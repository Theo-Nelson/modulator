#!/usr/bin/env python3

import argparse
import gzip
import os
import subprocess
import sys
import tempfile

import pandas as pd

from genotype_utils import load_read_assignments, normalize_string_series, run_process_jobs, sample_name_from_bam, safe_float, safe_int


def parse_args():
    ap = argparse.ArgumentParser(description="Build a per-read mod call table at candidate modulator sites.")
    ap.add_argument("--bams", nargs="+", required=True, help="Input BAMs with MM/ML tags")
    ap.add_argument("--candidate-sites-tsv", required=True, help="Candidate mod site TSV")
    ap.add_argument("--candidate-bed", required=True, help="Candidate mod BED for modkit include-bed")
    ap.add_argument("--read-assignments", required=True, help="Read assignment TSV")
    ap.add_argument("--reference-fa", required=True, help="Reference FASTA")
    ap.add_argument("--out-tsv", required=True, help="Output TSV")
    ap.add_argument("--modkit-bin", default="modkit", help="modkit executable")
    ap.add_argument("--threads", type=int, default=2)
    ap.add_argument("--jobs", type=int, default=1, help="Number of BAMs to process in parallel")
    ap.add_argument("--verbose", action="store_true")
    return ap.parse_args()


def parse_bool_text(x) -> bool:
    return str(x).strip().lower() in {"1", "true", "t", "yes", "y"}


def extract_rows_from_bam(
    bam: str,
    candidate_bed: str,
    reference_fa: str,
    modkit_bin: str,
    threads_per_job: int,
    lookup,
    verbose: bool = False,
    region=None,
):
    sample = sample_name_from_bam(bam)
    rows = []
    with tempfile.NamedTemporaryFile(prefix=f"{sample}.extract_calls.", suffix=".tsv.bgz", delete=False) as tmp:
        tmp_path = tmp.name

    try:
        cmd = [
            modkit_bin, "extract", "calls", bam, tmp_path,
            "--bgzf",
            "--force",
            *(["--region", str(region)] if region else []),
            # Cap modkit's per-chunk read buffering. The default 100kb interval over deep
            # direct-RNA piles up huge memory at highly-expressed loci (chr1/chr19 etc.) and
            # OOM'd even a 1TB node when many shards ran concurrently. Smaller chunks bound
            # peak RSS (more overhead, identical output).
            "--interval-size", "20000",
            # Don't estimate a pass-threshold by sampling reads: on sparse inputs
            # (e.g. region subsets, low-coverage samples) modkit aborts with
            # "Error! not enough datapoints" when there are too few mod calls over
            # the candidate-site BED. All calls are emitted; downstream genotype
            # logic applies its own coverage/quality filters.
            "--no-filtering",
            "--include-bed", candidate_bed,
            "--reference", reference_fa,
            "--mapped-only",
            "--threads", str(max(1, int(threads_per_job))),
            "--out-threads", "1",
            "--suppress-progress",
        ]
        if verbose:
            print(f"[info] mod extract start: {sample} threads={max(1, int(threads_per_job))}", file=sys.stderr, flush=True)
        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode != 0:
            raise SystemExit(f"modkit extract calls failed for {bam}:\n{proc.stderr}")

        header = None
        with gzip.open(tmp_path, "rt") as fh:
            for line in fh:
                line = line.rstrip("\n")
                if not line:
                    continue
                if header is None:
                    header = line.split("\t")
                    continue
                parts = line.split("\t")
                if len(parts) != len(header):
                    continue
                rec = dict(zip(header, parts))
                chrom = str(rec.get("chrom", ""))
                start0 = safe_int(rec.get("ref_position", -1), default=-1)
                qname = str(rec.get("read_id", ""))
                call_code = str(rec.get("call_code", ""))
                ref_strand = str(rec.get("ref_strand", ""))
                key = (chrom, start0)
                if key not in lookup:
                    continue
                for site in lookup[key]:
                    site_strand = str(site.get("strand", ""))
                    if site_strand and ref_strand and ref_strand not in {".", "?"} and site_strand != ref_strand:
                        continue
                    target_mod = str(site["mod_code"])
                    if call_code == target_mod:
                        state_detail = "modified"
                        target_modified = 1
                    elif call_code == "-":
                        state_detail = "canonical"
                        target_modified = 0
                    else:
                        state_detail = "other_mod"
                        target_modified = 0
                    rows.append({
                        "sample": sample,
                        "qname": qname,
                        "mod_site_id": site["mod_site_id"],
                        "chrom": chrom,
                        "start0": start0,
                        "end0": safe_int(site.get("end0", start0 + 1), default=start0 + 1),
                        "strand": site_strand or ref_strand,
                        "target_mod_code": target_mod,
                        "call_code": call_code,
                        "state_detail": state_detail,
                        "target_modified": target_modified,
                        "call_prob": safe_float(rec.get("call_prob", 0.0)),
                        "canonical_base": str(rec.get("canonical_base", "")),
                        "modified_primary_base": str(rec.get("modified_primary_base", "")),
                        "fail": parse_bool_text(rec.get("fail", False)),
                        "within_alignment": parse_bool_text(rec.get("within_alignment", True)),
                        "gene_id": str(site.get("gene_id", "")),
                        "gene_name": str(site.get("gene_name", "")),
                        "metagene_index": str(site.get("metagene_index", "")),
                    })
    finally:
        if os.path.exists(tmp_path):
            os.unlink(tmp_path)

    if verbose:
        print(f"[info] mod extract done: {sample} rows={len(rows)}", file=sys.stderr, flush=True)
    return rows


def main():
    args = parse_args()
    cand = pd.read_csv(args.candidate_sites_tsv, sep="\t", low_memory=False)
    if cand.empty:
        out = pd.DataFrame(columns=[
            "sample", "qname", "mod_site_id", "chrom", "start0", "end0", "strand",
            "target_mod_code", "call_code", "state_detail", "target_modified",
            "call_prob", "canonical_base", "modified_primary_base", "fail",
            "within_alignment", "gene_id", "gene_name", "metagene_index"
        ])
        os.makedirs(os.path.dirname(args.out_tsv) or ".", exist_ok=True)
        out.to_csv(args.out_tsv, sep="\t", index=False)
        return

    lookup = {}
    for row in cand.to_dict("records"):
        key = (str(row["chrom"]), int(row["start0"]))
        lookup.setdefault(key, []).append(row)

    # Shard per (BAM x chromosome) for genome-level parallelism: modkit extract runs
    # per chrom (--region) over its candidate sites; rows concatenate identically.
    chroms = sorted({k[0] for k in lookup.keys()})
    if not chroms:
        chroms = [None]
    lookup_by_chrom = {c: ({k: v for k, v in lookup.items() if k[0] == c} if c is not None else lookup) for c in chroms}
    n_tasks = len(args.bams) * len(chroms)
    jobs = max(1, min(int(args.jobs), n_tasks))
    threads_per_job = max(1, int(args.threads) // jobs)
    task_args = [
        (bam, args.candidate_bed, args.reference_fa, args.modkit_bin, threads_per_job, lookup_by_chrom[c], args.verbose, c)
        for bam in args.bams
        for c in chroms
    ]
    rows = []
    if jobs == 1:
        for item in task_args:
            rows.extend(extract_rows_from_bam(*item))
    else:
        for result in run_process_jobs(
            extract_rows_from_bam,
            task_args,
            jobs,
            verbose=args.verbose,
            label="build_molecule_mod_table",
        ):
            rows.extend(result)

    df = pd.DataFrame(rows)
    if df.empty:
        df = pd.DataFrame(columns=[
            "sample", "qname", "mod_site_id", "chrom", "start0", "end0", "strand",
            "target_mod_code", "call_code", "state_detail", "target_modified",
            "call_prob", "canonical_base", "modified_primary_base", "fail",
            "within_alignment", "gene_id", "gene_name", "metagene_index"
        ])
    else:
        # Load read assignments AFTER the parallel modkit extraction: read_assignments is
        # ~22GB on deep runs (~70GB in pandas), and loading it before the worker pool made
        # the fork copy-on-write that buffer across all jobs -> OOM. Workers never use it.
        assignments = load_read_assignments(args.read_assignments)
        keep_assign_cols = [c for c in [
            "sample", "qname", "ZT", "ZG", "ZN", "ZM", "assigned", "gene_id", "gene_name",
            "gene_index", "transcript_index", "metagene_index", "classification"
        ] if c in assignments.columns]
        assignments = assignments[keep_assign_cols].drop_duplicates(["sample", "qname"])
        assignments = assignments.rename(columns={
            col: f"assignment_{col}"
            for col in ["gene_id", "gene_name", "metagene_index"]
            if col in assignments.columns
        })
        df = df.merge(assignments, on=["sample", "qname"], how="left")
        # Keep site-derived context columns stable for downstream joins while
        # retaining assignment-derived metadata as explicit fallback columns.
        for col in ["gene_id", "gene_name", "metagene_index"]:
            assign_col = f"assignment_{col}"
            if col in df.columns and assign_col in df.columns:
                primary = normalize_string_series(df[col])
                fallback = normalize_string_series(df[assign_col])
                df[col] = primary.where(primary.ne(""), fallback)
        df["usable"] = (~df["fail"].fillna(True)) & df["within_alignment"].fillna(False)

    if not df.empty:
        # Deterministic on-disk order across (BAM x chrom) shards.
        df = df.sort_values(["chrom", "start0", "mod_site_id", "sample", "qname"]).reset_index(drop=True)

    os.makedirs(os.path.dirname(args.out_tsv) or ".", exist_ok=True)
    df.to_csv(args.out_tsv, sep="\t", index=False)


if __name__ == "__main__":
    main()
