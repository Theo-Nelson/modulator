#!/usr/bin/env python3

import argparse
import os
import shutil
import sys
import tempfile

import pandas as pd
import pysam

from genotype_utils import run_process_jobs, sample_name_from_bam, safe_int, write_chrom_index

OUT_COLS = [
    "sample", "qname", "snp_id", "chrom", "pos1", "start0", "end0", "ref", "alt",
    "observed_base", "allele_class", "baseq", "mapq", "strand",
    "ZT", "ZG", "ZN", "ZM", "gene_names", "gene_ids", "metagene_indices"
]


def parse_args():
    ap = argparse.ArgumentParser(description="Build a per-read candidate SNP table from tagged BAMs.")
    ap.add_argument("--bams", nargs="+", required=True, help="Input BAMs")
    ap.add_argument("--candidate-snps", required=True, help="Candidate SNP TSV from discover_candidate_snps.py")
    ap.add_argument("--out-tsv", required=True, help="Output molecule SNP TSV")
    ap.add_argument("--min-baseq", type=int, default=20)
    ap.add_argument("--min-mapq", type=int, default=10)
    ap.add_argument("--jobs", type=int, default=1, help="Number of BAMs to scan in parallel")
    ap.add_argument("--primary-only", action="store_true")
    ap.add_argument("--verbose", action="store_true")
    return ap.parse_args()


def safe_get_tag(aln, tag, default=""):
    try:
        return aln.get_tag(tag)
    except Exception:
        return default


def build_windows(positions, max_gap=5000):
    if not positions:
        return []
    positions = sorted(set(int(p) for p in positions))
    windows = []
    start = positions[0]
    prev = positions[0]
    for pos in positions[1:]:
        if pos - prev > max_gap:
            windows.append((start, prev))
            start = pos
        prev = pos
    windows.append((start, prev))
    return windows


def extract_rows_from_bam(
    bam: str,
    cand_by_chrom,
    min_baseq: int,
    min_mapq: int,
    primary_only: bool,
    verbose: bool = False,
    shard_dir=None,
):
    """Scan one (BAM x chromosome) and write its rows to a pickle shard; returns (chrom, path, n).
    Rows are never returned to the parent (162 GiB of dicts on an 8-sample mosquito run); the parent
    assembles the output one chromosome at a time."""
    sample = sample_name_from_bam(bam)
    if verbose:
        print(f"[info] molecule SNP start: {sample}", file=sys.stderr, flush=True)

    rows = []
    chrom_done = ""
    with pysam.AlignmentFile(bam, "rb") as fh:
        bam_contigs = set(fh.references)
        for chrom, pos_map in cand_by_chrom.items():
            chrom_done = chrom
            # A candidate SNP's contig may be absent from THIS sample's BAM header (e.g. samples aligned
            # to different references); pileup() on an unknown contig raises and would abort the scan.
            # Nothing to extract for this (bam, chrom), so skip it.
            if chrom not in bam_contigs:
                continue
            windows = build_windows(pos_map.keys())
            for win_start1, win_end1 in windows:
                # Per-alignment cache of the decoded sequence / qualities / tags: pileup hands the same
                # alignment back at every column it covers, and each `query_sequence` /
                # `query_qualities` access re-decodes the whole read (2 kb) -- at a 10k-deep SNP that is
                # 20 MB of decoding per column. Keyed by (qname, flag, start) so distinct records of one
                # read (primary_only=False) never share an entry. Bounded to one window.
                cache = {}
                for col in fh.pileup(
                    chrom,
                    win_start1 - 1,
                    win_end1,
                    truncate=True,
                    stepper="samtools",
                    min_base_quality=min_baseq,
                    min_mapping_quality=min_mapq,
                    max_depth=10_000_000,   # pysam defaults to 8000 -> a high-depth candidate SNP found
                ):                          # by discover_candidate_snps was silently capped here

                    pos1 = int(col.reference_pos) + 1
                    cand_row = pos_map.get(pos1)
                    if cand_row is None:
                        continue
                    ref = str(cand_row["ref"])
                    alt = str(cand_row["alt"])
                    for pr in col.pileups:
                        aln = pr.alignment
                        if pr.is_del or pr.is_refskip or pr.query_position is None:
                            continue
                        if primary_only and (aln.is_secondary or aln.is_supplementary):
                            continue
                        if aln.mapping_quality < min_mapq:
                            continue
                        ck = (aln.query_name, aln.flag, aln.reference_start)
                        c = cache.get(ck)
                        if c is None:
                            c = cache[ck] = (
                                aln.query_sequence, aln.query_qualities, int(aln.mapping_quality),
                                "-" if aln.is_reverse else "+",
                                str(safe_get_tag(aln, "ZT", "")), safe_int(safe_get_tag(aln, "ZG", "")),
                                safe_int(safe_get_tag(aln, "ZN", "")), safe_int(safe_get_tag(aln, "ZM", "")),
                            )
                        seq, quals, mapq, strand, zt, zg, zn, zm = c
                        qpos = pr.query_position
                        bq = int(quals[qpos]) if quals is not None else 0
                        if bq < min_baseq:
                            continue
                        base = seq[qpos].upper()
                        if len(base) != 1:
                            continue
                        allele_class = "other"
                        if base == ref:
                            allele_class = "ref"
                        elif base == alt:
                            allele_class = "alt"
                        rows.append({
                            "sample": sample,
                            "qname": ck[0],
                            "snp_id": cand_row["snp_id"],
                            "chrom": chrom,
                            "pos1": pos1,
                            "start0": pos1 - 1,
                            "end0": pos1,
                            "ref": ref,
                            "alt": alt,
                            "observed_base": base,
                            "allele_class": allele_class,
                            "baseq": bq,
                            "mapq": mapq,
                            "strand": strand,
                            "ZT": zt,
                            "ZG": zg,
                            "ZN": zn,
                            "ZM": zm,
                            "gene_names": cand_row.get("gene_names", ""),
                            "gene_ids": cand_row.get("gene_ids", ""),
                            "metagene_indices": cand_row.get("metagene_indices", ""),
                        })

    if verbose:
        print(f"[info] molecule SNP done: {sample} rows={len(rows)}", file=sys.stderr, flush=True)
    if not rows:
        return chrom_done, None, 0
    path = os.path.join(shard_dir, f"{sample}.{chrom_done}.pkl")
    pd.DataFrame(rows)[OUT_COLS].to_pickle(path)
    return chrom_done, path, len(rows)


def main():
    args = parse_args()
    cand = pd.read_csv(args.candidate_snps, sep="\t", low_memory=False)
    os.makedirs(os.path.dirname(args.out_tsv) or ".", exist_ok=True)
    if cand.empty:
        out = pd.DataFrame(columns=OUT_COLS)
        _tmp = args.out_tsv + ".tmp"           # atomic write (see build_read_assignment_table)
        out.to_csv(_tmp, sep="\t", index=False)
        os.replace(_tmp, args.out_tsv)
        write_chrom_index(args.out_tsv, [])
        return

    cand_by_chrom = {}
    for chrom, sub in cand.groupby("chrom", sort=False):
        pos_map = {}
        for row in sub.to_dict("records"):
            pos_map[int(row["pos1"])] = row
        cand_by_chrom[chrom] = pos_map

    shard_dir = tempfile.mkdtemp(prefix=".molsnp_shards.", dir=os.path.dirname(args.out_tsv) or ".")
    try:
        # Shard per (BAM x chromosome) for genome-level parallelism; results are assembled per
        # chromosome (sort + dedup), which reproduces the global sort exactly because chrom is the
        # primary sort key.
        task_args = [
            (bam, {chrom: pos_map}, args.min_baseq, args.min_mapq, args.primary_only, args.verbose, shard_dir)
            for bam in args.bams
            for chrom, pos_map in cand_by_chrom.items()
        ]
        jobs = max(1, min(int(args.jobs), len(task_args)))
        if jobs == 1:
            results = [extract_rows_from_bam(*item) for item in task_args]
        else:
            results = run_process_jobs(
                extract_rows_from_bam, task_args, jobs, verbose=args.verbose, label="build_molecule_snp_table")

        by_chrom = {}
        for chrom, path, n in results:
            if path is not None:
                by_chrom.setdefault(chrom, []).append(path)

        _tmp = args.out_tsv + ".tmp"               # atomic write (see build_read_assignment_table)
        blocks = []
        with open(_tmp, "w") as out:
            wrote_header = False
            if not by_chrom:
                pd.DataFrame(columns=OUT_COLS).to_csv(out, sep="\t", index=False)
                wrote_header = True
            for chrom in sorted(by_chrom):
                df = pd.concat([pd.read_pickle(p) for p in by_chrom[chrom]], ignore_index=True)
                # Sort BEFORE dedup so keep="first" is deterministic. Parallel (BAM x chrom) shards
                # return rows in nondeterministic order, and a read with >1 alignment over a SNP yields
                # >1 row (possibly with different observed_base); sorting -- with the allele as a
                # tiebreak -- makes the retained call reproducible instead of order-dependent.
                df = (df.sort_values(["chrom", "pos1", "snp_id", "sample", "qname", "observed_base"])
                        .drop_duplicates(["sample", "qname", "snp_id"], keep="first")
                        .reset_index(drop=True))
                off = out.tell()
                df.to_csv(out, sep="\t", index=False, header=not wrote_header)
                wrote_header = True
                out.flush()
                blocks.append((chrom, off, out.tell() - off))
                del df
        os.replace(_tmp, args.out_tsv)
        write_chrom_index(args.out_tsv, blocks, header_in_first_block=(len(blocks) > 0))
    finally:
        shutil.rmtree(shard_dir, ignore_errors=True)


if __name__ == "__main__":
    main()
