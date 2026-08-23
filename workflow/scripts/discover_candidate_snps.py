#!/usr/bin/env python3

import argparse
from collections import defaultdict
import os
import re
import shutil
import subprocess
import sys

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
                    help="samtools threads for the one-time per-BAM prefilter (R3).")
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


def annotate_site(chrom: str, pos1: int, exon_records):
    genes = set()
    gene_ids = set()
    metagenes = set()
    zts = set()
    for rec in exon_records.get(chrom, []):
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


def prefilter_bam(sample, in_bam, out_bam, exclude_flag, min_mapq, threads):
    """R3: one-time native prefilter so the per-shard count_coverage can run callback-free.
    `samtools view -F exclude_flag -q min_mapq` reproduces the old per-read callback EXACTLY:
    exclude_flag = UNMAP(0x4) [+ SECONDARY(0x100)|SUPPLEMENTARY(0x800) when primary_only], plus
    mapq>=min_mapq; QCFAIL/DUP are kept (the old callback never tested them). count_coverage then
    counts every remaining read (read_callback='nofilter'), so base counts stay byte-identical while
    the per-read Python callback (the hot path) is eliminated. Returns (sample, out_bam, chroms)."""
    threads = max(1, int(threads))
    subprocess.run(
        ["samtools", "view", "-b", "-F", str(int(exclude_flag)), "-q", str(max(0, int(min_mapq))),
         "-@", str(threads), "-o", out_bam, in_bam],
        check=True,
    )
    subprocess.run(["samtools", "index", "-@", str(threads), out_bam], check=True)
    return sample, out_bam, chroms_with_reads(out_bam)


def scan_bam_counts(
    sample: str,
    bam: str,
    reference_fa: str,
    merged_intervals,
    min_baseq: int,
    count_mode,
    verbose: bool = False,
):
    if verbose:
        print(f"[info] SNP scan start: {sample}", file=sys.stderr, flush=True)

    sample_counts = {}
    fasta = pysam.FastaFile(reference_fa)
    covered_positions = 0
    try:
        with pysam.AlignmentFile(bam, "rb") as fh:
            for chrom, intervals in merged_intervals.items():
                for start1, end1 in intervals:
                    ref_seq = fasta.fetch(chrom, start1 - 1, end1).upper()
                    counts = fh.count_coverage(
                        chrom,
                        start1 - 1,
                        end1,
                        quality_threshold=int(min_baseq),
                        read_callback=count_mode,
                    )
                    for offset, ref in enumerate(ref_seq):
                        if ref not in BASE_INDEX:
                            continue
                        base_counts = tuple(int(arr[offset]) for arr in counts)
                        if not any(base_counts):
                            continue
                        covered_positions += 1
                        sample_counts[(chrom, start1 + offset, ref)] = base_counts
    finally:
        fasta.close()

    if verbose:
        print(
            f"[info] SNP scan done: {sample} positions_with_coverage={covered_positions}",
            file=sys.stderr,
            flush=True,
        )
    return sample, sample_counts


def main():
    args = parse_args()
    exon_records, merged_intervals = load_gtf_exons(args.gtf)

    # R3: reproduce the old per-read callback with a one-time native prefilter, so each shard's
    # count_coverage runs callback-free ('nofilter'). Mode A (primary_only or min_mapq>0) prefilters;
    # the degenerate Mode B (no primary_only, no mapq) keeps pysam's default 'all' filter on the raw
    # BAM -- exactly what make_read_callback returned as "all" before this refactor.
    use_prefilter = bool(args.primary_only) or int(args.min_mapq) > 0
    # UNMAP | QCFAIL | DUP  + (SECONDARY|SUPPLEMENTARY if primary_only). QCFAIL(0x200) and DUP(0x400)
    # MUST be excluded: a marked duplicate is the same molecule (not independent allele support) and a
    # QC-failed read should not vote. Previously Mode A (the prefilter path, which the pipeline always
    # takes via --primary-only) counted both, while Mode B (pysam 'all') dropped them -- so the shipped
    # callset over-counted duplicate/QCFAIL reads AND the two modes disagreed (non-monotone in min_mapq).
    exclude_flag = 0x4 | 0x200 | 0x400 | (0x900 if args.primary_only else 0)
    # Per-PROCESS temp dir: a fixed name was shared by two concurrent CLI runs writing to the same
    # out dir, and the finally-cleanup then deleted the other run's prefiltered bams mid-scan.
    tmp_dir = os.path.join(os.path.dirname(args.out_tsv) or ".", f"._snp_prefilter_{os.getpid()}")

    bam_specs = []  # (sample, scan_bam_path, chroms_with_reads_or_None, count_mode)
    try:
        if use_prefilter:
            os.makedirs(tmp_dir, exist_ok=True)
            nb = max(1, len(args.bams))
            pf_threads = max(1, int(args.threads) // nb)
            pf_tasks = [
                (
                    sample_name_from_bam(b),
                    b,
                    # index-prefix the prefiltered path so two input BAMs sharing a basename do not
                    # collapse onto one file (their scans would otherwise clobber each other).
                    os.path.join(tmp_dir, f"{i:03d}_{sample_name_from_bam(b)}.prefiltered.bam"),
                    exclude_flag,
                    int(args.min_mapq),
                    pf_threads,
                )
                for i, b in enumerate(args.bams)
            ]
            pf_jobs = max(1, min(len(pf_tasks), int(args.jobs)))
            pf_results = run_process_jobs(
                prefilter_bam, pf_tasks, pf_jobs, verbose=args.verbose, label="prefilter_bam"
            )
            for sample, out_bam, chroms in pf_results:
                bam_specs.append((sample, out_bam, chroms, "nofilter"))
        else:
            for b in args.bams:
                sample = sample_name_from_bam(b)
                bam_specs.append((sample, b, chroms_with_reads(b), "all"))

        # R2 + R4: shard by genomic window, skipping chroms with no reads for a given BAM. The reduce
        # below keys by site, so per-shard partial counts merge identically to a per-(sample x chrom) scan.
        worker_args = []
        for sample, scan_bam, chroms, count_mode in bam_specs:
            for chrom, ivs in iter_window_shards(merged_intervals, args.window_bp):
                if chroms is not None and chrom not in chroms:
                    continue
                worker_args.append(
                    (sample, scan_bam, args.reference_fa, {chrom: ivs}, args.min_baseq, count_mode, args.verbose)
                )
        jobs = max(1, min(int(args.jobs), len(worker_args))) if worker_args else 1

        if jobs == 1:
            results = [scan_bam_counts(*item) for item in worker_args]
        else:
            results = run_process_jobs(
                scan_bam_counts,
                worker_args,
                jobs,
                verbose=args.verbose,
                label="discover_candidate_snps",
            )
    finally:
        if use_prefilter:
            shutil.rmtree(tmp_dir, ignore_errors=True)

    site_counts = {}
    for sample, sample_counts in results:
        for key, sample_base_counts in sample_counts.items():
            rec = site_counts.setdefault(key, {"counts": [0, 0, 0, 0], "per_sample": {}})
            rec["per_sample"][sample] = sample_base_counts
            for idx, value in enumerate(sample_base_counts):
                rec["counts"][idx] += int(value)

    rows = []
    # per-filter drop accounting: every candidate position removed by a filter is counted, so the scan
    # can report WHY sites were dropped instead of them vanishing silently.
    drops = {"low_total_cov": 0, "low_alt_reads": 0, "low_alt_frac": 0, "high_alt_frac": 0, "multiallelic": 0}
    n_candidate_positions = len(site_counts)
    for (chrom, pos1, ref), rec in sorted(site_counts.items()):
        total_counts = rec["counts"]
        counts = {base: int(total_counts[idx]) for idx, base in enumerate(DNA_BASES)}
        total_cov = sum(total_counts)
        ref_count = counts.get(ref, 0)
        alts = sorted(((b, c) for b, c in counts.items() if b != ref), key=lambda x: (-x[1], x[0]))
        alt, alt_count = alts[0] if alts else ("", 0)
        second_alt_count = alts[1][1] if len(alts) > 1 else 0
        alt_frac = (alt_count / total_cov) if total_cov > 0 else 0.0
        ref_frac = (ref_count / total_cov) if total_cov > 0 else 0.0
        if total_cov < args.min_total_cov:
            drops["low_total_cov"] += 1
            continue
        if alt_count < args.min_alt_reads:
            drops["low_alt_reads"] += 1
            continue
        if alt_frac < args.min_alt_frac:
            drops["low_alt_frac"] += 1
            continue
        if alt_frac > args.max_alt_frac:
            drops["high_alt_frac"] += 1
            continue
        # multiallelic: the SECOND alt must be a real fraction of coverage, not merely past the absolute
        # floor -- otherwise a deep clean biallelic het is discarded once accumulated basecall error on a
        # third base clears the floor (worse at higher depth, i.e. the highest-power sites). The second-alt
        # floor is its OWN parameter (--min-second-alt-reads), decoupled from --min-alt-reads.
        second_alt_frac = (second_alt_count / total_cov) if total_cov > 0 else 0.0
        if second_alt_count >= args.min_second_alt_reads and second_alt_frac >= args.multiallelic_frac:
            drops["multiallelic"] += 1
            continue

        ann = annotate_site(chrom, pos1, exon_records)
        sample_summaries = []
        samples_with_alt = 0
        alt_idx = BASE_INDEX.get(alt, -1)
        for sample in sorted(rec["per_sample"]):
            sc = rec["per_sample"][sample]
            if alt_idx >= 0 and sc[alt_idx] > 0:
                samples_with_alt += 1
            sample_summaries.append(
                f"{sample}:A={sc[0]},C={sc[1]},G={sc[2]},T={sc[3]}"
            )
        rows.append({
            "snp_id": f"{chrom}:{pos1}:{ref}>{alt}",
            "chrom": chrom,
            "pos1": pos1,
            "start0": pos1 - 1,
            "end0": pos1,
            "ref": ref,
            "alt": alt,
            "total_cov": total_cov,
            "ref_count": ref_count,
            "ref_frac": round(ref_frac, 6),
            "alt_count": alt_count,
            "second_alt_count": second_alt_count,
            "alt_frac": round(alt_frac, 6),
            "site_class": "segregating",
            "samples_with_alt": samples_with_alt,
            "gene_ids": ann["gene_ids"],
            "gene_names": ann["gene_names"],
            "metagene_indices": ann["metagene_indices"],
            "zt_labels": ann["zt_labels"],
            "sample_base_counts": "|".join(sample_summaries),
        })

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
    os.makedirs(os.path.dirname(args.out_tsv) or ".", exist_ok=True)
    _tmp = args.out_tsv + ".tmp"                # atomic write (see build_read_assignment_table)
    df.to_csv(_tmp, sep="\t", index=False)
    os.replace(_tmp, args.out_tsv)


if __name__ == "__main__":
    main()
