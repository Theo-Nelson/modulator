#!/usr/bin/env python3
"""
Aggregate modkit bedMethyl outputs partitioned by ZT (transcript code) into per-site tables.

Semantics (important):
- RAW = all parsed rows (after duplicate collapsing), no site filter.
- FILTERED = keep a site (chrom,start0,end0,strand,mod_code) if **ANY** row at that site
             passes row-level checks; if kept, include **all** rows at that site.

Outputs (when enabled):
  <out_prefix>_RAW_long.tsv
  <out_prefix>_RAW_frac_pivot.tsv
  <out_prefix>_RAW_cov_pivot.tsv
  <out_prefix>_RAW_Nmod_pivot.tsv
  <out_prefix>_FILTERED_long.tsv
  <out_prefix>_FILTERED_frac_pivot.tsv
  <out_prefix>_FILTERED_cov_pivot.tsv
  <out_prefix>_FILTERED_Nmod_pivot.tsv
"""

import os, sys, argparse, glob, gzip
from typing import Tuple

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
    ap = argparse.ArgumentParser(description="Aggregate ZT-partitioned modkit outputs per-site with filtering")
    ap.add_argument("--modkit-dir", required=True)
    ap.add_argument("--summary-tsv", help="Assembler summary to join on 'code'/'zt_label'")
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--min-cov", type=int, default=0,
                    help="Zero frac_modified for display when Nvalid_cov < MIN_COV (does NOT affect filtering)")

    # filtering knobs (row-level checks)
    ap.add_argument("--filter-enable", action="store_true")
    ap.add_argument("--count-diff-factor", type=float, default=3.0,
                    help="Row FAIL if Ndiff > factor * Nvalid_cov")
    ap.add_argument("--mod-fail-margin", type=int, default=1,
                    help="Row FAIL if Nmod <= Nfail + margin")

    # outputs
    ap.add_argument("--emit-raw", dest="emit_raw", action="store_true")
    ap.add_argument("--no-emit-raw", dest="emit_raw", action="store_false"); ap.set_defaults(emit_raw=True)
    ap.add_argument("--emit-filtered", dest="emit_filt", action="store_true")
    ap.add_argument("--no-emit-filtered", dest="emit_filt", action="store_false"); ap.set_defaults(emit_filt=True)
    ap.add_argument("--write-long", dest="write_long", action="store_true")
    ap.add_argument("--no-write-long", dest="write_long", action="store_false"); ap.set_defaults(write_long=True)
    ap.add_argument("--write-pivots", dest="write_pivots", action="store_true")
    ap.add_argument("--no-write-pivots", dest="write_pivots", action="store_false"); ap.set_defaults(write_pivots=True)

    # lightweight debug (summary counts only)
    ap.add_argument("--debug-summary", action="store_true",
                    help="Print counts of rows/sites passing/failing to stderr")

    ap.add_argument("--verbose", action="store_true")
    return ap.parse_args()

def is_header_line(line: str) -> bool:
    s = line.strip()
    return (not s) or s.startswith("#") or s.startswith("track") or s.startswith("browser")

def sample_and_code_from_path(path: str) -> Tuple[str, str]:
    """
    Support both:
      - Nested: <modkit_dir>/<SAMPLE>/<CODE>.bed(.gz)
      - Flat:   <modkit_dir>/<SAMPLE>_<CODE>.bed(.gz)
    """
    base = os.path.basename(path)
    if base.endswith(".bed.gz"): base = base[:-7]
    elif base.endswith(".bed"):  base = base[:-4]
    if base == "ungrouped": return None, None

    parent = os.path.basename(os.path.dirname(path))
    # nested: parent is sample, filename is code
    if parent and parent != os.path.basename(os.path.abspath(path)):
        return parent, base

    # flat fallback
    if "_" in base:
        return base.rsplit("_", 1)
    return None, None

def open_textmaybe_gzip(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")

def safe_int(x, default=0):
    try:
        return int(x)
    except Exception:
        try:
            return int(float(x))
        except Exception:
            return default

def parse_bed_row(parts):
    parts = parts[:18]
    d = dict(zip(BED_COLS, parts))
    d["start0"] = safe_int(d["start0"])
    d["end0"]   = safe_int(d["end0"])
    for k in ["Nvalid_cov","Nmod","Ncanonical","Nother_mod","Ndelete","Nfail","Ndiff","Nnocall"]:
        d[k] = safe_int(d[k])
    try:
        d["frac_modified"] = float(d["frac_modified"])
    except Exception:
        d["frac_modified"] = 0.0
    return d

def robust_load_summary(path: str, verbose=False) -> pd.DataFrame:
    if not path or not os.path.exists(path):
        return pd.DataFrame()
    lines = []
    with open(path) as f:
        for ln in f:
            if ln.strip(): lines.append(ln.rstrip("\n"))
    header = None; header_idx = 0
    for i, ln in enumerate(lines):
        h = ln.lstrip("#")
        if "\t" in h and any(c.strip() == "code" for c in h.split("\t")):
            header_idx = i; header = [c.strip() for c in h.split("\t")]
            break
    if header is None:
        header = [c.strip() for c in lines[0].lstrip("#").split("\t")]
        if verbose: print("[warn] Using first line as header for summary.", file=sys.stderr)

    rows = []
    for ln in lines[header_idx+1:]:
        if ln.startswith("#"): continue
        parts = ln.split("\t")
        if len(parts) < len(header): parts += [""]*(len(header)-len(parts))
        elif len(parts) > len(header): parts = parts[:len(header)]
        rows.append({header[j]: parts[j] for j in range(len(header))})

    df = pd.DataFrame(rows)
    if "code" not in df.columns:
        for c in list(df.columns):
            if c.lstrip("#").strip() == "code":
                df = df.rename(columns={c: "code"})
                break
    if "zt_label" in df.columns and "code" in df.columns:
        df["code"] = df["zt_label"].fillna(df["code"]).astype(str)
    elif "zt_label" in df.columns and "code" not in df.columns:
        df = df.rename(columns={"zt_label": "code"})
    if "code" in df.columns:
        df["code"] = df["code"].astype(str).str.strip()
    return df

def main():
    args = parse_args()

    beds = (
        sorted(glob.glob(os.path.join(args.modkit_dir, "*.bed"))) +
        sorted(glob.glob(os.path.join(args.modkit_dir, "*.bed.gz"))) +
        sorted(glob.glob(os.path.join(args.modkit_dir, "*", "*.bed"))) +
        sorted(glob.glob(os.path.join(args.modkit_dir, "*", "*.bed.gz")))
    )
    if not beds:
        sys.exit(f"No .bed or .bed.gz files found in {args.modkit_dir}")

    rows = []
    for bed in beds:
        sample, code = sample_and_code_from_path(bed)
        if sample is None or code is None:
            continue
        with open_textmaybe_gzip(bed) as f:
            for line in f:
                if is_header_line(line): continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 18:
                    parts = line.strip().split()
                    if len(parts) < 18: continue
                d = parse_bed_row(parts)
                rows.append({
                    "code": code, "sample": sample,
                    "chrom": d["chrom"], "start0": d["start0"], "end0": d["end0"], "strand": d["strand"],
                    "mod_code": d["mod_code"],
                    "Nvalid_cov": d["Nvalid_cov"], "Nmod": d["Nmod"], "Ncanonical": d["Ncanonical"],
                    "Nother_mod": d["Nother_mod"], "Ndelete": d["Ndelete"], "Nfail": d["Nfail"],
                    "Ndiff": d["Ndiff"], "Nnocall": d["Nnocall"],
                })

    if not rows:
        sys.exit("No data parsed from bed files; check filenames and contents.")

    df = pd.DataFrame(rows)

    # collapse duplicates within same genomic+code+sample+mod
    key = ["code","sample","mod_code","chrom","start0","end0","strand"]
    sumcols = ["Nvalid_cov","Nmod","Ncanonical","Nother_mod","Ndelete","Nfail","Ndiff","Nnocall"]
    df = df.groupby(key, as_index=False)[sumcols].sum()

    # frac (display only; min-cov does NOT affect filtering)
    df["frac_modified"] = (df["Nmod"] / df["Nvalid_cov"].where(df["Nvalid_cov"]>0, 1)).fillna(0.0)
    if args.min_cov:
        df.loc[df["Nvalid_cov"] < args.min_cov, "frac_modified"] = 0.0
    df["frac_modified"] = df["frac_modified"].round(6)

    # attach summary metadata if available
    if args.summary_tsv:
        summ = robust_load_summary(args.summary_tsv, verbose=args.verbose)
        if not summ.empty and "code" in summ.columns:
            keep_cols = [c for c in ["gtf_gene_name","gtf_gene_id","gtf_transcript_id",
                                     "classification","match_source","iso_tes","iso_chain_tx",
                                     "gtf_tes","gtf_chain_tx","tes_delta_bp","exon_overlap_bp",
                                     "read_support","frac_global","polya_support_frac","sample_counts",
                                     "gene_index","transcript_index","code"]
                         if c in summ.columns]
            summ_sub = summ[keep_cols].drop_duplicates("code")
            df = df.merge(summ_sub, on="code", how="left")

    # ---- SITE FILTERING (ANY-pass keeps site) ----
    # define site key
    df["__site_key__"] = list(zip(df["chrom"], df["start0"].astype(int), df["end0"].astype(int), df["strand"], df["mod_code"]))

    # row-level pass/fail
    if args.filter_enable:
        pass_row = (~(df["Ndiff"] > (args.count_diff_factor * df["Nvalid_cov"]))) & \
                   (df["Nmod"] > (df["Nfail"] + args.mod_fail_margin))
        passing_sites = set(df.loc[pass_row, "__site_key__"])
        df_filt = df[df["__site_key__"].isin(passing_sites)].copy()
    else:
        pass_row = pd.Series([True]*len(df), index=df.index)
        df_filt = df.copy()

    # debug summary, if requested
    if args.debug_summary:
        n_rows = len(df)
        n_rows_pass = int(pass_row.sum())
        n_sites = df["__site_key__"].nunique()
        n_sites_pass = len(set(df.loc[pass_row, "__site_key__"]))
        print(f"[debug] rows total={n_rows}, row-pass={n_rows_pass}", file=sys.stderr)
        print(f"[debug] sites total={n_sites}, site-pass(any)={n_sites_pass}", file=sys.stderr)

    base = args.out_prefix

    def write_long_and_pivots(sub: pd.DataFrame, tag: str):
        if args.write_long:
            out_long = f"{base}_{tag}_long.tsv"
            sub.drop(columns=["__site_key__"], errors="ignore").to_csv(out_long, sep="\t", index=False)
        if args.write_pivots:
            def pivot_metric(metric: str, fname_suffix: str):
                piv = sub.pivot_table(
                    index=["chrom","start0","end0","strand","code","mod_code"],
                    columns="sample", values=metric, aggfunc="first"
                ).fillna(0).reset_index()
                piv.to_csv(f"{base}_{tag}_{fname_suffix}.tsv", sep="\t", index=False)
            pivot_metric("frac_modified", "frac_pivot")
            pivot_metric("Nvalid_cov",   "cov_pivot")
            pivot_metric("Nmod",         "Nmod_pivot")

    # write RAW (pre-filter)
    if args.emit_raw:
        write_long_and_pivots(df, "RAW")

    # write FILTERED (post-site keep)
    if args.emit_filt:
        write_long_and_pivots(df_filt, "FILTERED")

    print("[OK] ZT aggregation complete.", file=sys.stderr)

if __name__ == "__main__":
    main()

