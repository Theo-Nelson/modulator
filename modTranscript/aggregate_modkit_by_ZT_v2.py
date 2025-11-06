#!/usr/bin/env python3
"""
Aggregate modkit bedMethyl outputs partitioned by ZT (transcript code).

Expected inputs:
  modkit_out/<SAMPLE>_<CODE>.bed
    where CODE looks like "GENE:abcdef12" and SAMPLE can contain underscores.
    We split on the LAST underscore to get (SAMPLE, CODE).

What it does:
  - For every bed file, sum the bedMethyl counts per modified-base code:
      Nvalid_cov, Nmod, Ncanonical, Nother_mod, Ndelete, Nfail, Ndiff, Nnocall
    and compute:
      frac_modified = Nmod / Nvalid_cov  (0 if Nvalid_cov==0)
  - Optionally join to a classification summary TSV from the assembler
    (columns like: code, gtf_gene_name, classification, iso_tes, iso_chain_tx, ...).
  - Write:
      <out_prefix>_long.tsv            (all metrics, long format)
      <out_prefix>_frac_pivot.tsv      (rows: code+mod, cols: samples, values: frac_modified)
      <out_prefix>_cov_pivot.tsv       (rows: code+mod, cols: samples, values: Nvalid_cov)
      <out_prefix>_Nmod_pivot.tsv      (rows: code+mod, cols: samples, values: Nmod)

Notes:
  - Skips files ending with '_ungrouped.bed'
  - Skips header/track lines starting with '#' or 'track' or 'browser'
  - Works with both letter (e.g., 'm','a') and numeric mod codes (e.g., '17596').

"""

import os, sys, argparse, csv, glob
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
    ap.add_argument("--modkit-dir", required=True, help="Directory containing <sample>_<code>.bed from modkit --partition-tag ZT")
    ap.add_argument("--summary-tsv", help="Optional assembler classification summary TSV to join on 'code'")
    ap.add_argument("--out-prefix", required=True, help="Prefix for output TSVs")
    ap.add_argument("--min-cov", type=int, default=0, help="If >0, zero-out frac_modified where summed Nvalid_cov < MIN_COV (kept in output)")
    ap.add_argument("--verbose", action="store_true")
    return ap.parse_args()

def is_header_line(line: str) -> bool:
    s = line.strip()
    return (not s) or s.startswith("#") or s.startswith("track") or s.startswith("browser")

def sample_and_code_from_path(path: str) -> Tuple[str, str]:
    base = os.path.basename(path)
    if base.endswith(".bed"):
        base = base[:-4]
    # skip *_ungrouped
    if base.endswith("_ungrouped"):
        return None, None
    # split on the last underscore
    if "_" not in base:
        return None, None
    sample, code = base.rsplit("_", 1)
    return sample, code

def safe_int(x, default=0):
    try:
        return int(x)
    except Exception:
        try:
            # modkit sometimes has floats in frac col only; others should be ints
            return int(float(x))
        except Exception:
            return default

def safe_float(x, default=0.0):
    try:
        return float(x)
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
                # some modkit builds may use mixed delim after col10; try to split on whitespace past 10th col
                parts = line.strip().split()
                if len(parts) < 18:
                    continue
            d = dict(zip(BED_COLS, parts[:18]))
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

def load_summary(summary_tsv: str) -> pd.DataFrame:
    # Robust reader; some summaries start with '#code' etc.
    df = pd.read_csv(summary_tsv, sep="\t", comment="#")
    if "code" not in df.columns:
        # some versions have a leading '#code' header
        for c in df.columns:
            if c.lstrip("#") == "code":
                df = df.rename(columns={c: "code"})
                break
    # Keep a useful subset; pass-thru other cols as well
    return df

def main():
    args = parse_args()
    in_dir = args.modkit_dir
    beds = sorted(glob.glob(os.path.join(in_dir, "*.bed")))
    if not beds:
        sys.exit(f"No .bed files found in {in_dir}")

    rows_long: List[dict] = []
    seen_samples = set()
    seen_codes = set()
    seen_mods = set()

    for bed in beds:
        sample, code = sample_and_code_from_path(bed)
        if sample is None or code is None:
            if args.verbose:
                print(f"[skip] {os.path.basename(bed)} (no code or ungrouped)", file=sys.stderr)
            continue
        seen_samples.add(sample)
        seen_codes.add(code)

        sums_by_mod = read_bed_summing_counts(bed)
        for mod, cnts in sums_by_mod.items():
            seen_mods.add(mod)
            Ncov = cnts["Nvalid_cov"]
            Nmod = cnts["Nmod"]
            frac = (Nmod / Ncov) if Ncov > 0 else 0.0
            if args.min_cov and Ncov < args.min_cov:
                # keep the row but blank out frac to 0 (explicit)
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
    if args.summary_tsv and os.path.exists(args.summary_tsv):
        summ = load_summary(args.summary_tsv)
        # Keep useful annotation columns if present
        keep_cols = [c for c in ["gtf_gene_name","gtf_gene_id","gtf_transcript_id",
                                 "classification","match_source","iso_tes","iso_chain_tx",
                                 "gtf_tes","gtf_chain_tx","tes_delta_bp","exon_overlap_bp",
                                 "read_support","frac_global","polya_support_frac","sample_counts"]
                     if c in summ.columns]
        summ_sub = summ[["code"] + keep_cols].drop_duplicates("code")
        df_long = df_long.merge(summ_sub, on="code", how="left")
    else:
        if args.verbose:
            print("[warn] No summary TSV provided; output will omit annotation columns.", file=sys.stderr)

    # Write long table
    out_long = f"{args.out_prefix}_long.tsv"
    df_long.to_csv(out_long, sep="\t", index=False)
    if args.verbose:
        print(f"[ok] wrote {out_long}")

    # Build pivots
    def pivot_metric(metric: str, fname_suffix: str):
        piv = df_long.pivot_table(index=["code","mod_code"], columns="sample", values=metric, aggfunc="first").reset_index()
        piv = piv.fillna(0)
        outp = f"{args.out_prefix}_{fname_suffix}.tsv"
        piv.to_csv(outp, sep="\t", index=False)
        if args.verbose:
            print(f"[ok] wrote {outp}")

    pivot_metric("frac_modified", "frac_pivot")
    pivot_metric("Nvalid_cov",   "cov_pivot")
    pivot_metric("Nmod",         "Nmod_pivot")

    print("[OK] Aggregation complete.")

if __name__ == "__main__":
    main()

