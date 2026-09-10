#!/usr/bin/env python3

import argparse
import bisect
from collections import defaultdict
import os
import re
import sys

import numpy as np
import pandas as pd
import pysam

from genotype_utils import run_process_jobs, sample_name_from_bam


DNA_BASES = ("A", "C", "G", "T")
BASE_INDEX = {base: idx for idx, base in enumerate(DNA_BASES)}


def parse_args():
    ap = argparse.ArgumentParser(description="Discover segregating candidate SNPs from tagged BAMs inside assembled transcript loci.")
    ap.add_argument("--bams", nargs="+", required=True, help="Input BAMs")
    ap.add_argument("--reference-fa", required=True, help="Reference FASTA")
    ap.add_argument("--gtf", required=True, help="Assembler GTF used to define transcribed loci")
    ap.add_argument("--out-tsv", required=True, help="Output candidate SNP TSV")
    ap.add_argument("--min-alt-reads", type=int, default=4)
    ap.add_argument("--min-second-alt-reads", type=int, default=4,
                    help="Absolute read floor for the SECOND alt in the multiallelic drop. Kept SEPARATE "
                         "from --min-alt-reads (default matches its default): sharing one value made the "
                         "caller non-monotone -- raising --min-alt-reads TIGHTENED the first-alt floor "
                         "but simultaneously LOOSENED the multiallelic gate, so sites appeared only at a "
                         "stricter setting.")
    ap.add_argument("--min-total-cov", type=int, default=8)
    ap.add_argument("--min-alt-frac", type=float, default=0.10)
    ap.add_argument("--max-alt-frac", type=float, default=0.90)
    ap.add_argument("--multiallelic-frac", type=float, default=0.10,
                    help="A site is dropped as multiallelic only if its SECOND-most-common alt is both "
                         ">= min-alt-reads AND >= this FRACTION of coverage. A pure absolute count wrongly "
                         "discards deep clean biallelic hets, where a few %% third-base basecall error "
                         "accumulates past min-alt-reads at high depth. Set 0 for the old absolute-only rule.")
    ap.add_argument("--min-baseq", type=int, default=20)
    ap.add_argument("--min-mapq", type=int, default=10)
    ap.add_argument("--jobs", type=int, default=1, help="Number of scan shards to run in parallel")
    ap.add_argument("--window-bp", type=int, default=1_000_000,
                    help="Split each chromosome's exon intervals into shards spanning at most this "
                         "many bp, so parallelism isn't capped at one shard per chromosome (R2).")
    ap.add_argument("--threads", type=int, default=4,
                    help="(accepted for compatibility; the per-BAM prefilter copy was removed)")
    ap.add_argument("--primary-only", action="store_true")
    ap.add_argument("--verbose", action="store_true")
    return ap.parse_args()


def parse_attrs(attrs: str):
    out = {}
    for m in re.finditer(r'(\S+)\s+"([^"]*)"', attrs):
        out[m.group(1)] = m.group(2)
    return out


def load_gtf_exons(gtf_path: str):
    exon_records = defaultdict(list)
    merged_intervals = defaultdict(list)
    with open(gtf_path) as fh:
        for ln in fh:
            if ln.startswith("#") or not ln.strip():
                continue
            parts = ln.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            chrom, _, feature, start, end, _, strand, _, attrs = parts
            if feature != "exon":
                continue
            a = parse_attrs(attrs)
            start = int(start)
            end = int(end)
            exon_records[chrom].append({
                "start": start,
                "end": end,
                "strand": strand,
                "gene_id": a.get("gene_id", ""),
                "gene_name": a.get("ref_gene_name", a.get("gene_name", a.get("gene_id", ""))),
                "metagene_index": a.get("metagene_index", ""),
                "zt_label": a.get("zt_label", a.get("transcript_id", "")),
            })
            merged_intervals[chrom].append((start, end))

    for chrom, intervals in merged_intervals.items():
        intervals = sorted(intervals)
        merged = []
        for s, e in intervals:
            if not merged or s > merged[-1][1] + 1:
                merged.append([s, e])
            else:
                merged[-1][1] = max(merged[-1][1], e)
        merged_intervals[chrom] = [(s, e) for s, e in merged]
        exon_records[chrom].sort(key=lambda x: (x["start"], x["end"], x["gene_name"], x["zt_label"]))
    return exon_records, merged_intervals


_EXON_SCAN = {}


def _exon_window(recs, pos1):
    """Exon records that can contain pos1 (bisect on start + max exon length), in list order --
    replaces a linear scan from the chromosome start for every candidate."""
    if not recs:
        return recs
    idx = _EXON_SCAN.get(id(recs))
    if idx is None or idx[0] is not recs:
        starts = [r["start"] for r in recs]
        max_len = max((r["end"] - r["start"] for r in recs), default=0)
        idx = (recs, starts, max_len)
        _EXON_SCAN[id(recs)] = idx
    _, starts, max_len = idx
    lo = bisect.bisect_left(starts, pos1 - max_len)
    hi = bisect.bisect_right(starts, pos1)
    return recs[lo:hi]


def annotate_site(chrom: str, pos1: int, exon_records):
    genes = set()
    gene_ids = set()
    metagenes = set()
    zts = set()
    for rec in _exon_window(exon_records.get(chrom, []), pos1):
        if rec["start"] > pos1:
            break
        if rec["end"] < pos1:
            continue
        gene_ids.add(rec["gene_id"])
        genes.add(rec["gene_name"])
        if rec["metagene_index"]:
            metagenes.add(str(rec["metagene_index"]))
        if rec["zt_label"]:
            zts.add(rec["zt_label"])
    return {
        "gene_ids": ";".join(sorted(g for g in gene_ids if g)),
        "gene_names": ";".join(sorted(g for g in genes if g)),
        "metagene_indices": ";".join(sorted(m for m in metagenes if m)),
        "zt_labels": ";".join(sorted(z for z in zts if z)),
    }


def init_count_record():
    return {
        "counts": {b: 0 for b in DNA_BASES},
        "per_sample": defaultdict(lambda: {b: 0 for b in DNA_BASES}),
    }


def iter_window_shards(merged_intervals, window_bp):
    """R2: split each chrom's (sorted, merged, non-overlapping) intervals into shards spanning at
    most ~window_bp, keeping every interval whole. Because no interval is split, positions are
    partitioned across shards (each position lands in exactly one shard), so the cross-sample site
    reduce merges byte-identically to a single per-chromosome scan. Yields (chrom, [intervals])."""
    window_bp = max(1, int(window_bp))
    for chrom, ivs in merged_intervals.items():
        cur = []
        cur_start = None
        for s, e in ivs:
            if cur and (e - cur_start) > window_bp:
                yield chrom, cur
                cur = []
                cur_start = None
            if cur_start is None:
                cur_start = s
            cur.append((s, e))
        if cur:
            yield chrom, cur


def chroms_with_reads(bam_path):
    """R4: chromosomes with >=1 mapped read (idxstats col 3). Returns None on failure (never skip)."""
    try:
        out = set()
        for line in pysam.idxstats(bam_path).splitlines():
            f = line.split("\t")
            if len(f) >= 3 and f[0] != "*" and int(f[2]) > 0:
                out.add(f[0])
        return out
    except Exception:
        return None


# One AlignmentFile per BAM per worker (workers are long-lived; re-opening parses the index).
_BAMS = {}


def _bam(path):
    fh = _BAMS.get(path)
    if fh is None:
        fh = pysam.AlignmentFile(path, "rb")
        _BAMS[path] = fh
    return fh


def _make_read_filter(exclude_flag, min_mapq):
    """The old design copied every BAM through `samtools view -F <flag> -q <mapq>` (31 x ~35 GB of
    transient copies at cohort scale) so count_coverage could run callback-free. The zt_tagged BAMs
    already contain no secondary/supplementary records, so the filter is MAPQ + QCFAIL/DUP: a per-read
    callback reproduces it exactly, at ~1 us/read, with no copy."""
    def _ok(r):
        return r.mapping_quality >= min_mapq and not (r.flag & exclude_flag)
    return _ok


def scan_shard(chrom, intervals, bam_specs, reference_fa, min_baseq, filt, verbose=False):
    """Count bases for ONE genomic shard over ALL samples and reduce to candidate rows inside the
    worker. Returns (rows, n_positions_with_coverage, drops).

    Positions are partitioned across shards, so per-shard reduction is exactly the old global
    reduce -- but the parent never holds every covered position x every sample (34 GiB at 2
    samples, ~300 GiB projected at 31). The per-base Python loop is replaced by numpy over the
    (4 x L) count arrays; only candidate positions are touched in Python."""
    exclude_flag, min_mapq, primary_or_mapq, args = filt
    cb = _make_read_filter(exclude_flag, min_mapq) if primary_or_mapq else "all"
    fasta = pysam.FastaFile(reference_fa)
    rows = []
    n_cov_positions = 0
    drops = {"low_total_cov": 0, "low_alt_reads": 0, "low_alt_frac": 0, "high_alt_frac": 0, "multiallelic": 0}
    try:
        for start1, end1 in intervals:
            ref_seq = fasta.fetch(chrom, start1 - 1, end1).upper()
            L = len(ref_seq)
            if L <= 0:
                continue
            ref_idx = np.full(L, -1, dtype=np.int64)
            rb = np.frombuffer(ref_seq.encode("ascii"), dtype=np.uint8)
            for k, base in enumerate(DNA_BASES):
                ref_idx[rb == ord(base)] = k
            total = np.zeros((4, L), dtype=np.int64)
            per_sample = []   # (sample, (4,L) int array) -- one window at a time, bounded
            for sample, bam_path, chroms in bam_specs:
                if chroms is not None and chrom not in chroms:
                    continue
                fh = _bam(bam_path)
                if chrom not in fh.references:
                    continue
                # NB: count_coverage filters on raw base quality only -- it does NOT apply BAQ (see the
                # ADVANCED_USAGE note, finding N); counts match `samtools mpileup -B`.
                counts = fh.count_coverage(chrom, start1 - 1, end1, quality_threshold=int(min_baseq),
                                           read_callback=cb)
                arr = np.array(counts, dtype=np.int64)
                if not arr.any():
                    continue
                total += arr
                per_sample.append((sample, arr))
            if not per_sample:
                continue
            ar = np.arange(L)
            total_cov = total.sum(axis=0)
            valid = (ref_idx >= 0) & (total_cov > 0)
            n_cov_positions += int(valid.sum())
            if not valid.any():
                continue
            ridx = np.where(valid, ref_idx, 0)
            ref_count = total[ridx, ar]
            masked = total.copy(); masked[ridx, ar] = -1
            alt_idx = masked.argmax(axis=0)                      # ties -> smallest base index = alphabetical
            alt_count = masked[alt_idx, ar]
            masked2 = masked.copy(); masked2[alt_idx, ar] = -1
            second_alt_count = masked2.max(axis=0)
            tc = total_cov.astype(float)
            alt_frac = np.where(total_cov > 0, alt_count / np.where(total_cov > 0, tc, 1), 0.0)
            second_alt_frac = np.where(total_cov > 0, second_alt_count / np.where(total_cov > 0, tc, 1), 0.0)
            # sequential filters with drop accounting (same order + semantics as the old per-position loop)
            m = valid.copy()
            f1 = m & (total_cov < args.min_total_cov); drops["low_total_cov"] += int(f1.sum()); m &= ~f1
            f2 = m & (alt_count < args.min_alt_reads); drops["low_alt_reads"] += int(f2.sum()); m &= ~f2
            f3 = m & (alt_frac < args.min_alt_frac); drops["low_alt_frac"] += int(f3.sum()); m &= ~f3
            f4 = m & (alt_frac > args.max_alt_frac); drops["high_alt_frac"] += int(f4.sum()); m &= ~f4
            f5 = m & (second_alt_count >= args.min_second_alt_reads) & (second_alt_frac >= args.multiallelic_frac)
            drops["multiallelic"] += int(f5.sum()); m &= ~f5
            cand = np.flatnonzero(m)
            if cand.size == 0:
                continue
            for off in cand.tolist():
                pos1 = start1 + off
                ref = ref_seq[off]
                alt = DNA_BASES[int(alt_idx[off])]
                total_cov_i = int(total_cov[off]); ref_count_i = int(ref_count[off]); alt_count_i = int(alt_count[off])
                sample_summaries = []
                samples_with_alt = 0
                alt_i = int(alt_idx[off])
                for sample, arr in sorted(per_sample, key=lambda t: t[0]):
                    sc = arr[:, off]
                    if not sc.any():
                        continue   # this sample has no coverage here (the old code stored no entry)
                    if sc[alt_i] > 0:
                        samples_with_alt += 1
                    sample_summaries.append(f"{sample}:A={int(sc[0])},C={int(sc[1])},G={int(sc[2])},T={int(sc[3])}")
                ann = annotate_site(chrom, pos1, args._exon_records)
                rows.append({
                    "snp_id": f"{chrom}:{pos1}:{ref}>{alt}",
                    "chrom": chrom,
                    "pos1": pos1,
                    "start0": pos1 - 1,
                    "end0": pos1,
                    "ref": ref,
                    "alt": alt,
                    "total_cov": total_cov_i,
                    "ref_count": ref_count_i,
                    "ref_frac": round((ref_count_i / total_cov_i) if total_cov_i > 0 else 0.0, 6),
                    "alt_count": alt_count_i,
                    "second_alt_count": int(second_alt_count[off]),
                    "alt_frac": round((alt_count_i / total_cov_i) if total_cov_i > 0 else 0.0, 6),
                    "site_class": "segregating",
                    "samples_with_alt": samples_with_alt,
                    "gene_ids": ann["gene_ids"],
                    "gene_names": ann["gene_names"],
                    "metagene_indices": ann["metagene_indices"],
                    "zt_labels": ann["zt_labels"],
                    "sample_base_counts": "|".join(sample_summaries),
                })
    finally:
        fasta.close()
    if verbose:
        print(f"[info] SNP scan shard {chrom}:{intervals[0][0]}-{intervals[-1][1]}: "
              f"{n_cov_positions} covered positions, {len(rows)} candidates", file=sys.stderr, flush=True)
    return rows, n_cov_positions, drops


_SHARD_G = {}


def _init_shard_worker(exon_records):
    _SHARD_G["exon_records"] = exon_records


def _scan_shard_task(task):
    chrom, intervals, bam_specs, reference_fa, min_baseq, filt = task
    filt[3]._exon_records = _SHARD_G["exon_records"]
    return scan_shard(chrom, intervals, bam_specs, reference_fa, min_baseq, filt)


def main():
    args = parse_args()
    exon_records, merged_intervals = load_gtf_exons(args.gtf)

    # UNMAP | QCFAIL | DUP  + (SECONDARY|SUPPLEMENTARY if primary_only), applied per read (see
    # _make_read_filter); Mode B (no primary_only, no mapq) keeps pysam's default 'all' filter.
    use_filter = bool(args.primary_only) or int(args.min_mapq) > 0
    exclude_flag = 0x4 | 0x200 | 0x400 | (0x900 if args.primary_only else 0)

    bam_specs = []
    for b in args.bams:
        bam_specs.append((sample_name_from_bam(b), b, chroms_with_reads(b)))
    all_chroms = set()
    for _s, _b, ch in bam_specs:
        if ch is None:
            all_chroms = None
            break
        all_chroms |= ch

    # One task per genomic window (all samples), so the reduction to candidates happens inside the
    # worker and the parent only ever receives candidate rows (R2 sharding; R4 chrom skip).
    filt = (exclude_flag, int(args.min_mapq), use_filter, args)
    tasks = []
    for chrom, ivs in iter_window_shards(merged_intervals, args.window_bp):
        if all_chroms is not None and chrom not in all_chroms:
            continue
        tasks.append((chrom, ivs, bam_specs, args.reference_fa, args.min_baseq, filt))
    jobs = max(1, min(int(args.jobs), len(tasks))) if tasks else 1

    rows = []
    n_candidate_positions = 0
    drops = {"low_total_cov": 0, "low_alt_reads": 0, "low_alt_frac": 0, "high_alt_frac": 0, "multiallelic": 0}

    def _take(res):
        nonlocal n_candidate_positions
        r, n, d = res
        rows.extend(r)
        n_candidate_positions += n
        for k, v in d.items():
            drops[k] += v

    if jobs == 1:
        _init_shard_worker(exon_records)
        for t in tasks:
            _take(_scan_shard_task(t))
    else:
        import multiprocessing
        from concurrent.futures import ProcessPoolExecutor
        ctx = multiprocessing.get_context("fork")
        _init_shard_worker(exon_records)   # inherited by forked workers (no pickling of the GTF index)
        with ProcessPoolExecutor(max_workers=jobs, mp_context=ctx) as ex:
            for res in ex.map(_scan_shard_task, tasks, chunksize=1):
                _take(res)

    rows.sort(key=lambda r: (r["chrom"], r["pos1"], r["ref"]))
    df = pd.DataFrame(rows)
    # de-silence the scan: report how many candidate positions each filter removed (was invisible).
    if args.verbose:
        print(f"[candidate_snps] {n_candidate_positions:,} positions with coverage -> {len(rows):,} kept; "
              f"dropped: low_total_cov={drops['low_total_cov']:,} low_alt_reads={drops['low_alt_reads']:,} "
              f"low_alt_frac={drops['low_alt_frac']:,} high_alt_frac={drops['high_alt_frac']:,} "
              f"multiallelic={drops['multiallelic']:,}", file=sys.stderr, flush=True)
    out_dir = os.path.dirname(args.out_tsv) or "."
    os.makedirs(out_dir, exist_ok=True)
    if df.empty:
        df = pd.DataFrame(columns=[
            "snp_id", "chrom", "pos1", "start0", "end0", "ref", "alt", "total_cov", "ref_count",
            "ref_frac", "alt_count", "second_alt_count", "alt_frac", "site_class", "samples_with_alt", "gene_ids",
            "gene_names", "metagene_indices", "zt_labels", "sample_base_counts"
        ])
    _tmp = args.out_tsv + ".tmp"                # atomic write (see build_read_assignment_table)
    df.to_csv(_tmp, sep="\t", index=False)
    os.replace(_tmp, args.out_tsv)


if __name__ == "__main__":
    main()
