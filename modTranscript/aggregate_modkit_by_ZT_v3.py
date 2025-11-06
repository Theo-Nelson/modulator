#!/usr/bin/env python3
"""
Aggregate modkit bedMethyl outputs partitioned by ZT (transcript code).

Changes vs v2:
  • Robust summary loader: treats '#code' as header (not a comment),
    tolerates occasional missing tabs by padding/truncating rows to header length.
  • Will never crash if summary lacks 'code' — it falls back to no-join with a warning.

Outputs:
  <out_prefix>_long.tsv
  <out_prefix>_frac_pivot.tsv
  <out_prefix>_cov_pivot.tsv
  <out_prefix>_Nmod_pivot.tsv
"""

import os, sys, argparse, glob
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
                    help="Directory with <sample>_<CODE>.bed from `modkit pileup --partition-tag ZT`")
    ap.add_argument("--summary-tsv",
                    help="Optional assembler classification summary TSV to join on 'code' (header may start with '#code')")
    ap.add_argument("--out-prefix", required=True, help="Prefix for output TSVs")
    ap.add_argument("--min-cov", type=int, default=0,
                    help="If >0, zero-out frac_modified where summed Nvalid_cov < MIN_COV (row kept)")
    ap.add_argument("--verbose", action="store_true")
    return ap.parse_args()

def is_header_line(line: str) -> bool:
    s = line.strip()
    return (not s) or s.startswith("#") or s.startswith("track") or s.startswith("browser")

def sample_and_code_from_path(path: str) -> Tuple[str, str]:
    """Infer sample and code by splitting at the last underscore."""
    base = os.path.basename(path)
    if base.endswith(".bed"):
        base = base[:-4]
    if base.endswith("_ungrouped"):
        return None, None
    if "_" not in base:
        return None, None
    sample, code = base.rsplit("_", 1)
    return sample, code

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
    with open(path) as f:
        for line in f:
            if is_header_line(line): 
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 18:
                # Some modkit builds use mixed delim beyond col10; fallback to any whitespace.
                parts = line.strip().split()
                if len(parts) < 18:
                    continue
            # Only take first 18 fields to be safe
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
    """
    Read the assembler classification summary even if header begins with '#code'
    and rows occasionally miss a tab. Pads/truncates rows to header length.
    """
    if not os.path.exists(path):
        raise FileNotFoundError(path)
    lines = []
    with open(path) as f:
        for ln in f:
            if ln.strip() == "":
                continue
            lines.append(ln.rstrip("\n"))

    # Find header line: prefer one containing 'code' (ignoring leading '#')
    header_idx = None
    header = None
    for i, ln in enumerate(lines):
        h = ln.lstrip("#")
        if "\t" in h and "code" in [c.strip() for c in h.split("\t")]:
            header_idx = i
            header = [c.strip() for c in h.split("\t")]
            break
    if header is None:
        # Fallback: use the first non-empty as header
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
        # pad/truncate to header length (prevents crashes on a missing tab)
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

    if "code" not in df.columns:
        if verbose:
            print(f"[warn] Summary '{path}' has no 'code' column after parsing; join will be skipped.", file=sys.stderr)
        return pd.DataFrame()

    # Basic cleanup
    df["code"] = df["code"].astype(str).str.strip()
    return df

def main():
    args = parse_args()

    beds = sorted(glob.glob(os.path.join(args.modkit_dir, "*.bed")))
    if not beds:
        sys.exit(f"No .bed files found in {args.modkit_dir}")

    rows_long: List[dict] = []

    for bed in beds:
        sample, code = sample_and_code_from_path(bed)
        if sample is None or code is None:
            if args.verbose:
                print(f"[skip] {os.path.basename(bed)} (no code or ungrouped)", file=sys.stderr)
            continue

        sums_by_mod = read_bed_summing_counts(bed)
        for mod, cnts in sums_by_mod.items():
            Ncov = cnts["Nvalid_cov"]
            Nmod = cnts["Nmod"]
            frac = (Nmod / Ncov) if Ncov > 0 else 0.0
            if args.min_cov and Ncov < args.min_cov:
                frac = 0.0  # keep row; explicit zero for low coverage
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
            # Keep a useful subset if present
            keep_cols = [c for c in ["gtf_gene_name","gtf_gene_id","gtf_transcript_id",
                                     "classification","match_source","iso_tes","iso_chain_tx",
                                     "gtf_tes","gtf_chain_tx","tes_delta_bp","exon_overlap_bp",
                                     "read_support","frac_global","polya_support_frac","sample_counts"]
                         if c in summ.columns]
            summ_sub = summ[["code"] + keep_cols].drop_duplicates("code")
            pre = len(df_long)
            df_long = df_long.merge(summ_sub, on="code", how="left")
            if args.verbose:
                n_joined = df_long["gtf_gene_name"].notna().sum() if "gtf_gene_name" in df_long.columns else "NA"
                print(f"[ok] merged summary on 'code' (rows={pre} -> {len(df_long)}; annotated rows={n_joined})", file=sys.stderr)
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

