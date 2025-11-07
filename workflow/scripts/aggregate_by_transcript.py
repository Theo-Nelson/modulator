#!/usr/bin/env python3
"""
Aggregate modkit bedMethyl outputs partitioned by ZT (transcript code).

v4 changes:
  • Compatible with assembler v7 ZT labels like "GENE.GENEID.G# .T#" (dots, no underscores)
    so we can still split <sample>_<ZT>.bed on the last underscore safely even when
    sample names contain underscores.
  • Adds optional parsing of .bed.gz and defensive filename handling.
  • Writes the same four outputs as v3: *_long.tsv, *_frac_pivot.tsv, *_cov_pivot.tsv, *_Nmod_pivot.tsv
  • Joins optional classification summary and recognizes 'zt_label', 'gene_index', 'transcript_index'.

Usage is unchanged except you should run modkit with --partition-tag ZT after using
assembler v7, which sets ZT to a deterministic transcript label.
"""

import os, sys, argparse, glob, gzip
from collections import defaultdict
from typing import Dict, Tuple, List

try:
    import pandas as pd
except ImportError:
    sys.exit("This script requires pandas. Install it (e.g. `micromamba install pandas`).")

BED_COLS = [
    "chrom", "start0", "end0", "mod_code", "score", "strand",
    "start0_compat", "end0_compat", "rgb",
    "Nvalid_cov", "frac_modified",
    "Nmod", "Ncanonical", "Nother_mod",
    "Ndelete", "Nfail", "Ndiff", "Nnocall"
]


def parse_args():
    ap = argparse.ArgumentParser(description="Aggregate modkit bedMethyl files partitioned by ZT")
    ap.add_argument("--modkit-dir", required=True,
                    help="Directory with <sample>_<ZT>.bed from `modkit pileup --partition-tag ZT`")
    ap.add_argument("--summary-tsv",
                    help="Optional assembler classification summary TSV to join on 'code'/'zt_label' (header may start with '#code')")
    ap.add_argument("--out-prefix", required=True, help="Prefix for output TSVs")
    ap.add_argument("--min-cov", type=int, default=0,
                    help="If >0, zero-out frac_modified where summed Nvalid_cov < MIN_COV (row kept)")
    ap.add_argument("--verbose", action="store_true")
    return ap.parse_args()


def is_header_line(line: str) -> bool:
    s = line.strip()
    return (not s) or s.startswith("#") or s.startswith("track") or s.startswith("browser")


def sample_and_code_from_path(path: str) -> Tuple[str, str]:
    """Infer sample and partition code from filename by splitting at the last underscore.
    Works if partition (ZT) contains no underscores. assembler v7 ensures this by using dots.
    Returns (sample, code) or (None, None) if the pattern doesn't match.
    """
    base = os.path.basename(path)
    if base.endswith(".bed.gz"):
        base = base[:-7]
    elif base.endswith(".bed"):
        base = base[:-4]
    # Ignore trailing _ungrouped if present
    if base.endswith("_ungrouped"):
        return None, None
    if "_" not in base:
        return None, None
    sample, code = base.rsplit("_", 1)
    return sample, code


def open_textmaybe_gzip(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def safe_int(x, default=0):
    try:
        return int(x)
    except Exception:
        try:
            return int(float(x))
        except Exception:
            return default


def read_bed_summing_counts(path: str) -> Dict[str, Dict[str, int]]:
    """
    Sum counts per mod_code across the entire bed file.
    Returns: mod_code -> counts dict
    """
    sums: Dict[str, Dict[str, int]] = defaultdict(lambda: {
        "Nvalid_cov": 0, "Nmod": 0, "Ncanonical": 0, "Nother_mod": 0,
        "Ndelete": 0, "Nfail": 0, "Ndiff": 0, "Nnocall": 0, "sites": 0
    })
    with open_textmaybe_gzip(path) as f:
        for line in f:
            if is_header_line(line):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 18:
                parts = line.strip().split()
                if len(parts) < 18:
                    continue
            parts = parts[:18]
            d = dict(zip(BED_COLS, parts))
            mod = d["mod_code"]
            sums[mod]["Nvalid_cov"] += safe_int(d["Nvalid_cov"])
            sums[mod]["Nmod"]       += safe_int(d["Nmod"])
            sums[mod]["Ncanonical"] += safe_int(d["Ncanonical"])
            sums[mod]["Nother_mod"] += safe_int(d["Nother_mod"])
            sums[mod]["Ndelete"]    += safe_int(d["Ndelete"])
            sums[mod]["Nfail"]      += safe_int(d["Nfail"])
            sums[mod]["Ndiff"]      += safe_int(d["Ndiff"])
            sums[mod]["Nnocall"]    += safe_int(d["Nnocall"])
            sums[mod]["sites"]      += 1
    return sums


def robust_load_summary(path: str, verbose=False) -> pd.DataFrame:
    if not path or not os.path.exists(path):
        return pd.DataFrame()
    lines = []
    with open(path) as f:
        for ln in f:
            if ln.strip() == "":
                continue
            lines.append(ln.rstrip("\n"))

    header_idx = None
    header = None
    for i, ln in enumerate(lines):
        h = ln.lstrip("#")
        if "\t" in h and any(c.strip() == "code" for c in h.split("\t")):
            header_idx = i
            header = [c.strip() for c in h.split("\t")]
            break
    if header is None:
        header_idx = 0
        header = [c.strip() for c in lines[0].lstrip("#").split("\t")]
        if verbose:
            print("[warn] Could not find explicit 'code' header; using first line as header", file=sys.stderr)

    rows = []
    expected = len(header)
    for ln in lines[header_idx+1:]:
        if ln.startswith("#"):
            continue
        parts = ln.split("\t")
        if len(parts) < expected:
            parts = parts + [""] * (expected - len(parts))
        elif len(parts) > expected:
            parts = parts[:expected]
        rows.append({header[j]: parts[j] for j in range(expected)})

    df = pd.DataFrame(rows)

    # Normalize 'code' column name if it's '#code' or similar
    if "code" not in df.columns:
        for c in list(df.columns):
            if c.lstrip("#").strip() == "code":
                df = df.rename(columns={c: "code"})
                break

    # Prefer 'zt_label' if present
    if "zt_label" in df.columns and "code" in df.columns:
        df["code"] = df["zt_label"].fillna(df["code"]).astype(str)
    elif "zt_label" in df.columns and "code" not in df.columns:
        df = df.rename(columns={"zt_label": "code"})

    if "code" in df.columns:
        df["code"] = df["code"].astype(str).str.strip()

    return df


def main():
    args = parse_args()

    # Look for .bed files either directly in modkit_dir or in sample subdirectories
    beds = sorted(glob.glob(os.path.join(args.modkit_dir, "*.bed"))) + \
           sorted(glob.glob(os.path.join(args.modkit_dir, "*.bed.gz"))) + \
           sorted(glob.glob(os.path.join(args.modkit_dir, "*", "*.bed"))) + \
           sorted(glob.glob(os.path.join(args.modkit_dir, "*", "*.bed.gz")))
    if not beds:
        sys.exit(f"No .bed or .bed.gz files found in {args.modkit_dir}")

    rows_long: List[dict] = []

    for bed in beds:
        sample, code = sample_and_code_from_path(bed)
        if sample is None or code is None:
            if args.verbose:
                print(f"[skip] {os.path.basename(bed)} (unparseable filename)", file=sys.stderr)
            continue

        sums_by_mod = read_bed_summing_counts(bed)
        for mod, cnts in sums_by_mod.items():
            Ncov = cnts["Nvalid_cov"]
            Nmod = cnts["Nmod"]
            frac = (Nmod / Ncov) if Ncov > 0 else 0.0
            if args.min_cov and Ncov < args.min_cov:
                frac = 0.0
            rows_long.append({
                "code": code,
                "sample": sample,
                "mod_code": mod,
                "Nvalid_cov": Ncov,
                "Nmod": Nmod,
                "Ncanonical": cnts["Ncanonical"],
                "Nother_mod": cnts["Nother_mod"],
                "Ndelete": cnts["Ndelete"],
                "Nfail": cnts["Nfail"],
                "Ndiff": cnts["Ndiff"],
                "Nnocall": cnts["Nnocall"],
                "n_sites": cnts["sites"],
                "frac_modified": round(frac, 6),
            })

    if not rows_long:
        sys.exit("No data parsed from bed files; check filenames and contents.")

    df_long = pd.DataFrame(rows_long)

    # Optional join to classification summary
    if args.summary_tsv:
        summ = robust_load_summary(args.summary_tsv, verbose=args.verbose)
        if not summ.empty and "code" in summ.columns:
            keep_cols = [c for c in ["gtf_gene_name","gtf_gene_id","gtf_transcript_id",
                                     "classification","match_source","iso_tes","iso_chain_tx",
                                     "gtf_tes","gtf_chain_tx","tes_delta_bp","exon_overlap_bp",
                                     "read_support","frac_global","polya_support_frac","sample_counts",
                                     "gene_index","transcript_index"]
                         if c in summ.columns]
            summ_sub = summ[["code"] + keep_cols].drop_duplicates("code")
            pre = len(df_long)
            df_long = df_long.merge(summ_sub, on="code", how="left")
            if args.verbose:
                print(f"[ok] merged summary on 'code' (rows={pre} -> {len(df_long)})", file=sys.stderr)
        else:
            if args.verbose:
                print("[warn] Summary join skipped (no 'code' column or empty).", file=sys.stderr)

    # Write long table
    out_long = f"{args.out_prefix}_long.tsv"
    df_long.to_csv(out_long, sep="\t", index=False)
    if args.verbose:
        print(f"[ok] wrote {out_long}", file=sys.stderr)

    # Build pivots
    def pivot_metric(metric: str, fname_suffix: str):
        piv = df_long.pivot_table(index=["code","mod_code"], columns="sample", values=metric, aggfunc="first")
        piv = piv.fillna(0).reset_index()
        outp = f"{args.out_prefix}_{fname_suffix}.tsv"
        piv.to_csv(outp, sep="\t", index=False)
        if args.verbose:
            print(f"[ok] wrote {outp}", file=sys.stderr)

    pivot_metric("frac_modified", "frac_pivot")
    pivot_metric("Nvalid_cov",   "cov_pivot")
    pivot_metric("Nmod",         "Nmod_pivot")

    print("[OK] Aggregation complete.")

if __name__ == "__main__":
    main()

