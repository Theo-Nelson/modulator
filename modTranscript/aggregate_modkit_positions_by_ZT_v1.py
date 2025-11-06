#!/usr/bin/env python3
"""
Make a per-position table from modkit bedMethyl outputs partitioned by ZT (transcript code).

Output rows: one per (sample, code, chrom, 1-based position, strand, mod_code)
Only positions with Nvalid_cov >= --min-cov are kept.

Columns written:
  code, sample, chrom, pos1, strand, mod_code,
  Nvalid_cov, Nmod, frac_modified,
  Ncanonical, Nother_mod, Ndelete, Nfail, Ndiff, Nnocall, n_sites,
  pos_key  (chrom:pos1:strand)
+ (optional) annotation from assembler summary if --summary-tsv is provided.

Notes:
- Robust to bedMethyl "mixed delimiters" (tabs first 10 fields, spaces after).
- Skips *_ungrouped.bed and files without a trailing "_<CODE>.bed".
- Tolerant summary loader: accepts header starting with '#code' and pads/truncates ragged rows.
"""

import os, sys, glob, argparse
from collections import defaultdict
from typing import Tuple, Dict, List

try:
    import pandas as pd
except ImportError:
    sys.exit("This script requires pandas. Install it (e.g., `micromamba install pandas`).")

BED_COLS = [
    "chrom","start0","end0","mod_code","score","strand",
    "start0_compat","end0_compat","rgb",
    "Nvalid_cov","frac_modified","Nmod","Ncanonical","Nother_mod",
    "Ndelete","Nfail","Ndiff","Nnocall"
]

def parse_args():
    ap = argparse.ArgumentParser(description="Aggregate modkit bedMethyl positions by ZT transcript code")
    ap.add_argument("--modkit-dir", required=True, help="Directory with *.bed from `modkit pileup --partition-tag ZT`")
    ap.add_argument("--summary-tsv", help="Optional assembler summary TSV (header may begin with #code)")
    ap.add_argument("--out-prefix", required=True, help="Prefix for output TSV")
    ap.add_argument("--min-cov", type=int, default=5, help="Min Nvalid_cov to keep a position")
    ap.add_argument("--verbose", action="store_true")
    return ap.parse_args()

def is_header_line(s: str) -> bool:
    s = s.strip()
    return (not s) or s.startswith("#") or s.startswith("track") or s.startswith("browser")

def sample_and_code_from_path(path: str) -> Tuple[str, str]:
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

def robust_load_summary(path: str, verbose=False) -> pd.DataFrame:
    if not path or not os.path.exists(path):
        return pd.DataFrame()
    lines = [ln.rstrip("\n") for ln in open(path) if ln.strip() != ""]
    header_idx, header = None, None
    for i, ln in enumerate(lines):
        cand = ln.lstrip("#")
        if "\t" in cand and "code" in [c.strip() for c in cand.split("\t")]:
            header_idx = i
            header = [c.strip() for c in cand.split("\t")]
            break
    if header is None:
        header_idx = 0
        header = [c.strip() for c in lines[0].lstrip("#").split("\t")]
        if verbose:
            print("[warn] summary: no explicit 'code' header; using first line", file=sys.stderr)
    rows = []
    expected = len(header)
    for ln in lines[header_idx+1:]:
        if ln.startswith("#"):
            continue
        parts = ln.split("\t")
        if len(parts) < expected:
            parts = parts + [""]*(expected-len(parts))
        elif len(parts) > expected:
            parts = parts[:expected]
        rows.append({header[j]: parts[j] for j in range(expected)})
    df = pd.DataFrame(rows)
    # normalize '#code' -> 'code'
    if "code" not in df.columns:
        for c in list(df.columns):
            if c.lstrip("#").strip() == "code":
                df = df.rename(columns={c: "code"})
                break
    if "code" in df.columns:
        df["code"] = df["code"].astype(str).str.strip()
    else:
        if verbose:
            print("[warn] summary has no 'code' column; annotation join will be skipped", file=sys.stderr)
        return pd.DataFrame()
    return df

def parse_bed_positions(path: str, min_cov: int) -> List[dict]:
    """Return list of per-position dicts for positions with Nvalid_cov >= min_cov."""
    out = []
    with open(path) as f:
        for line in f:
            if is_header_line(line):
                continue
            # mixed-delim safe: split on any whitespace
            parts = line.strip().split()
            if len(parts) < 18:
                # try tab split fallback just in case
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 18:
                    continue
            parts = parts[:18]
            d = dict(zip(BED_COLS, parts))
            Ncov = safe_int(d["Nvalid_cov"])
            if Ncov < min_cov:
                continue
            Nmod = safe_int(d["Nmod"])
            # 1-based pos center = start0+1 (bedMethyl is single-base intervals)
            pos1 = safe_int(d["start0"]) + 1
            frac = (Nmod/Ncov) if Ncov > 0 else 0.0
            out.append({
                "chrom": d["chrom"],
                "pos1": pos1,
                "strand": d["strand"],
                "mod_code": d["mod_code"],
                "Nvalid_cov": Ncov,
                "Nmod": Nmod,
                "frac_modified": round(frac, 6),
                "Ncanonical": safe_int(d["Ncanonical"]),
                "Nother_mod": safe_int(d["Nother_mod"]),
                "Ndelete": safe_int(d["Ndelete"]),
                "Nfail": safe_int(d["Nfail"]),
                "Ndiff": safe_int(d["Ndiff"]),
                "Nnocall": safe_int(d["Nnocall"]),
                "n_sites": 1  # each bed row is one site
            })
    return out

def main():
    args = parse_args()

    beds = sorted(glob.glob(os.path.join(args.modkit_dir, "*.bed")))
    if not beds:
        sys.exit(f"No .bed files found in {args.modkit_dir}")

    rows = []
    for bed in beds:
        sample, code = sample_and_code_from_path(bed)
        if sample is None or code is None:
            if args.verbose:
                print(f"[skip] {os.path.basename(bed)} (no code or ungrouped)", file=sys.stderr)
            continue
        perpos = parse_bed_positions(bed, args.min_cov)
        if not perpos:
            if args.verbose:
                print(f"[info] {os.path.basename(bed)} yielded 0 positions at cov>={args.min_cov}", file=sys.stderr)
            continue
        for r in perpos:
            r["sample"] = sample
            r["code"] = code
            r["pos_key"] = f"{r['chrom']}:{r['pos1']}:{r['strand']}"
            rows.append(r)

    if not rows:
        sys.exit("No positions passed the coverage filter. Try lowering --min-cov or check inputs.")

    df = pd.DataFrame(rows)

    # Optional join to assembler summary for convenience
    if args.summary_tsv:
        summ = robust_load_summary(args.summary_tsv, verbose=args.verbose)
        if not summ.empty and "code" in summ.columns:
            keep_cols = [c for c in [
                "gtf_gene_name","gtf_gene_id","gtf_transcript_id",
                "classification","match_source",
                "iso_tes","iso_chain_tx","gtf_tes","gtf_chain_tx",
                "tes_delta_bp","exon_overlap_bp",
                "read_support","frac_global","polya_support_frac","sample_counts"
            ] if c in summ.columns]
            df = df.merge(summ[["code"] + keep_cols].drop_duplicates("code"), on="code", how="left")
        elif args.verbose:
            print("[warn] Summary join skipped (no 'code' column or empty).", file=sys.stderr)

    out_tsv = f"{args.out_prefix}_positions_long.tsv"
    df.sort_values(["chrom","pos1","strand","mod_code","sample","code"]).to_csv(out_tsv, sep="\t", index=False)
    if args.verbose:
        print(f"[ok] wrote {out_tsv}  (rows={len(df)})", file=sys.stderr)

    print("[OK] Position-level aggregation complete.")

if __name__ == "__main__":
    main()

