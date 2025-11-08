#!/usr/bin/env python3
"""
Aggregate modkit ZN-partitioned bedMethyl outputs into per-site (genomic) and per-gene tables,
with duplicate-site collapsing and explicit transcript index (ZN).

Key features:
- Reads numbered ZN partition files (e.g., .../<sample>/<mod>_filtered_mod.bed/1.bed).
  We ignore 'ungrouped.bed' and any flat '*_filtered_mod.bed' files.
- Each numbered file contains bedMethyl rows with mod_code in column 4; we keep that column.
- Deduplication: sums counts for identical keys:
  (sample, mod_code, ZN_transcript_index, chrom, start0, end0, strand).
- Site-level filtering with hyperparameters:
    FAIL if (Ndiff > count_diff_factor * Nvalid_cov) OR (Nmod <= Nfail + mod_fail_margin).
  A SITE IS KEPT if **any** row at that site (across samples/transcripts) passes;
  when a site is kept, we keep **all rows** for that site.
- Outputs (controlled by flags):
    RAW  : before site filtering
    FILTERED: after site filtering rule above
  For each of RAW / FILTERED (as enabled):
    1) <out_prefix>_<TAG>_sites_long.tsv
    2) Optional per-gene×mod tables with pivots under
       <out_prefix>_<TAG>__per_gene_mod/
"""

import os, sys, re, argparse, gzip
from collections import defaultdict, namedtuple
from typing import List, Dict, Tuple

try:
    import pandas as pd
except ImportError:
    sys.exit("This script requires pandas. (e.g., `micromamba install pandas`)")

BED_COLS = [
    "chrom","start0","end0","mod_code","score","strand",
    "start0_compat","end0_compat","rgb",
    "Nvalid_cov","frac_modified",
    "Nmod","Ncanonical","Nother_mod","Ndelete","Nfail","Ndiff","Nnocall",
]

def parse_args():
    ap = argparse.ArgumentParser(description="Aggregate ZN-partitioned modkit outputs per gene/mod with site-level filtering")
    ap.add_argument("--modkit-dir", required=True, help="Parent dir with per-sample subdirs containing numbered ZN .bed files")
    ap.add_argument("--gtf", required=True, help="Assembler GTF (with gene coordinates). Exon or transcript features work.")
    ap.add_argument("--out-prefix", required=True, help="Prefix for outputs (no extension)")
    ap.add_argument("--min-cov", type=int, default=0, help="Zero frac_modified if Nvalid_cov < MIN_COV (row kept)")
    ap.add_argument("--filter-enable", action="store_true", help="Enable site-level filtering")
    ap.add_argument("--count-diff-factor", type=float, default=3.0, help="FAIL if Ndiff > factor * Nvalid_cov (default: 3)")
    ap.add_argument("--mod-fail-margin", type=int, default=1, help="FAIL if Nmod <= Nfail + margin (default: 1)")

    # output toggles
    ap.add_argument("--emit-raw", dest="emit_raw", action="store_true")
    ap.add_argument("--no-emit-raw", dest="emit_raw", action="store_false"); ap.set_defaults(emit_raw=True)
    ap.add_argument("--emit-filtered", dest="emit_filt", action="store_true")
    ap.add_argument("--no-emit-filtered", dest="emit_filt", action="store_false"); ap.set_defaults(emit_filt=True)

    ap.add_argument("--write-long", dest="write_long", action="store_true")
    ap.add_argument("--no-write-long", dest="write_long", action="store_false"); ap.set_defaults(write_long=True)
    ap.add_argument("--write-pivots", dest="write_pivots", action="store_true")
    ap.add_argument("--no-write-pivots", dest="write_pivots", action="store_false"); ap.set_defaults(write_pivots=True)

    # RAW vs FILTERED per-gene tables
    ap.add_argument("--write-raw-per-gene", dest="write_raw_per_gene", action="store_true")
    ap.add_argument("--no-write-raw-per-gene", dest="write_raw_per_gene", action="store_false"); ap.set_defaults(write_raw_per_gene=False)
    ap.add_argument("--write-filtered-per-gene", dest="write_filtered_per_gene", action="store_true")
    ap.add_argument("--no-write-filtered-per-gene", dest="write_filtered_per_gene", action="store_false"); ap.set_defaults(write_filtered_per_gene=True)

    ap.add_argument("--verbose", action="store_true")
    return ap.parse_args()

def open_text(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")

def is_header(line:str)->bool:
    s=line.strip()
    return (not s) or s.startswith("#") or s.startswith("track") or s.startswith("browser")

def safe_int(x, default=0):
    try:
        return int(x)
    except Exception:
        try:
            return int(float(x))
        except Exception:
            return default

# --- GTF interval indexing ---
Interval = namedtuple("Interval", ["start","end","gene_id","gene_name","strand"])

def load_gene_intervals_from_gtf(gtf_path:str, verbose=False)->Dict[Tuple[str,str], List[Interval]]:
    """Union per-gene exon/transcript spans to coarse intervals for site→gene mapping."""
    gene_bounds: Dict[Tuple[str,str,str], Tuple[int,int]] = {}
    gene_name_map: Dict[str,str] = {}
    with open_text(gtf_path) as f:
        for ln in f:
            if ln.startswith("#") or not ln.strip(): continue
            parts = ln.rstrip("\n").split("\t")
            if len(parts) < 9: continue
            chrom, source, feature, start, end, score, strand, frame, attrs = parts
            if feature not in ("exon","transcript","gene"): continue
            a = {}
            for kv in re.finditer(r'(\S+)\s+"([^"]*)"', attrs):
                a[kv.group(1)] = kv.group(2)
            gene_id = a.get("gene_id") or a.get("gtf_gene_id") or a.get("gene") or ""
            gene_name = a.get("ref_gene_name") or a.get("gene_name") or a.get("gtf_gene_name") or gene_id
            if not gene_id: continue
            s = int(start); e = int(end)
            key = (chrom, strand, gene_id)
            if key not in gene_bounds: gene_bounds[key] = (s, e)
            else:
                mn, mx = gene_bounds[key]
                gene_bounds[key] = (min(mn, s), max(mx, e))
            gene_name_map[gene_id] = gene_name

    by_cs: Dict[Tuple[str,str], List[Interval]] = defaultdict(list)
    for (chrom, strand, gid), (s,e) in gene_bounds.items():
        gname = gene_name_map.get(gid, gid)
        by_cs[(chrom, strand)].append(Interval(s, e, gid, gname, strand))
    for k in by_cs:
        by_cs[k].sort(key=lambda iv: (iv.start, iv.end))
    if verbose:
        print(f"[info] loaded {sum(len(v) for v in by_cs.values())} gene intervals from {gtf_path}", file=sys.stderr)
    return by_cs

def assign_gene(chrom:str, pos_start:int, pos_end:int, strand:str, gene_index:Dict[Tuple[str,str], List[Interval]]):
    """Return (gene_id, gene_name) by overlap; choose max-overlap; tie → first; try opposite strand if empty."""
    ivs = gene_index.get((chrom, strand), [])
    best = None; best_ov = -1
    for iv in ivs:
        if iv.start > pos_end: break
        if iv.end < pos_start: continue
        ov = min(iv.end, pos_end) - max(iv.start, pos_start) + 1
        if ov > best_ov: best_ov = ov; best = iv
    if best: return best.gene_id, best.gene_name
    ivs2 = gene_index.get((chrom, "+" if strand == "-" else "-"), [])
    for iv in ivs2:
        if iv.start > pos_end: break
        if iv.end < pos_start: continue
        return iv.gene_id, iv.gene_name
    return "", ""

# --- readers ---
def iter_numbered_beds(modkit_dir:str)->List[Tuple[str, str, str, int]]:
    """
    Yield (root, sample_name, bed_path, ZN_index) for files like '<sample>/<something>/<N>.bed'.
    Skip ungrouped and flat '*_filtered_mod.bed'.
    """
    out = []
    for root, dirs, files in os.walk(modkit_dir):
        rel = os.path.relpath(root, modkit_dir)
        if rel == ".":             # expect sample subdirs
            continue
        sample_name = rel.split(os.sep)[0]
        for fname in files:
            if fname.endswith("_filtered_mod.bed") or fname.endswith("_filtered_mod.bed.gz"):
                continue
            base = fname[:-3] if fname.endswith(".gz") else fname
            if base.lower() == "ungrouped.bed": continue
            m = re.fullmatch(r"(\d+)\.bed", base)
            if not m: continue
            zn = int(m.group(1))
            out.append((root, sample_name, os.path.join(root, fname), zn))
    return sorted(out)

def parse_bed_line(line:str):
    parts = line.rstrip("\n").split("\t")
    if len(parts) < 18:
        parts = line.strip().split()
        if len(parts) < 18:
            return None
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

def site_key_from_row(r):
    # group per genomic site + mod type (independent of transcript or sample)
    return (r["chrom"], int(r["start0"]), int(r["end0"]), r["strand"], r["mod_code"])

def main():
    args = parse_args()
    beds = iter_numbered_beds(args.modkit_dir)
    if not beds:
        sys.exit(f"No numbered ZN partition files found under {args.modkit_dir}")

    gene_index = load_gene_intervals_from_gtf(args.gtf, verbose=args.verbose)

    rows = []
    for _, sample_name, bed_path, zn in beds:
        with open_text(bed_path) as f:
            for ln in f:
                if is_header(ln): continue
                rec = parse_bed_line(ln)
                if not rec: continue
                gid, gname = assign_gene(rec["chrom"], rec["start0"], rec["end0"], rec["strand"], gene_index)
                rows.append({
                    "sample": sample_name,
                    "ZN_transcript_index": zn,
                    "chrom": rec["chrom"],
                    "start0": rec["start0"],
                    "end0": rec["end0"],
                    "strand": rec["strand"],
                    "mod_code": rec["mod_code"],
                    "Nvalid_cov": rec["Nvalid_cov"],
                    "Nmod": rec["Nmod"],
                    "Ncanonical": rec["Ncanonical"],
                    "Nother_mod": rec["Nother_mod"],
                    "Ndelete": rec["Ndelete"],
                    "Nfail": rec["Nfail"],
                    "Ndiff": rec["Ndiff"],
                    "Nnocall": rec["Nnocall"],
                    "gene_id": gid,
                    "gene_name": gname,
                })

    if not rows:
        sys.exit("Parsed zero rows from numbered ZN beds.")

    df = pd.DataFrame(rows)

    # collapse duplicates within same genomic+partition key
    key = ["sample","mod_code","ZN_transcript_index","chrom","start0","end0","strand","gene_id","gene_name"]
    sumcols = ["Nvalid_cov","Nmod","Ncanonical","Nother_mod","Ndelete","Nfail","Ndiff","Nnocall"]
    pre = len(df)
    df = df.groupby(key, as_index=False)[sumcols].sum()
    if args.verbose:
        print(f"[info] collapsed {pre} → {len(df)} rows", file=sys.stderr)

    # frac_modified (+ min-cov zeroing only affects the displayed fraction, not filtering)
    df["frac_modified"] = (df["Nmod"] / df["Nvalid_cov"].where(df["Nvalid_cov"]>0, 1)).fillna(0.0)
    if args.min_cov:
        df.loc[df["Nvalid_cov"] < args.min_cov, "frac_modified"] = 0.0
    df["frac_modified"] = df["frac_modified"].round(6)

    # --- site-level filtering ---
    df["__site_key__"] = df.apply(site_key_from_row, axis=1)

    if args.filter_enable:
        # row-pass logic
        pass_row = (~(df["Ndiff"] > (args.count_diff_factor * df["Nvalid_cov"]))) & \
                   (df["Nmod"] > (df["Nfail"] + args.mod_fail_margin))
        # sites that pass if ANY row at the site passes
        passing_sites = set(df.loc[pass_row, "__site_key__"])
        df_filt = df[df["__site_key__"].isin(passing_sites)].copy()
    else:
        df_filt = df.copy()

    df_raw = df.copy()  # before filtering
    base = args.out_prefix

    def write_all(sub: pd.DataFrame, tag: str, write_per_gene: bool):
        # long
        if args.write_long:
            out_long = f"{base}_{tag}_sites_long.tsv"
            sub.drop(columns=["__site_key__"], errors="ignore").to_csv(out_long, sep="\t", index=False)
            if args.verbose: print(f"[ok] wrote {out_long}", file=sys.stderr)

        # per-gene × mod tables + pivots
        if write_per_gene or args.write_pivots:
            out_dir = f"{base}_{tag}__per_gene_mod"
            os.makedirs(out_dir, exist_ok=True)
            prefix_base = os.path.basename(args.out_prefix)
            for (gname, mod), gsub in sub.groupby(["gene_name","mod_code"], dropna=False):
                safe_g = re.sub(r"[^A-Za-z0-9._+-]", "_", gname if gname else "NA")
                safe_mod = re.sub(r"[^A-Za-z0-9._+-]", "_", str(mod))
                fn_base = os.path.join(out_dir, f"{prefix_base}__{safe_g}__{safe_mod}")

                cols = ["gene_name","gene_id","mod_code","chrom","start0","end0","strand",
                        "ZN_transcript_index","sample","Nvalid_cov","Nmod","Ncanonical",
                        "Nother_mod","Ndelete","Nfail","Ndiff","Nnocall","frac_modified"]
                gsub2 = gsub[cols].sort_values(["chrom","start0","ZN_transcript_index","sample"])

                if write_per_gene:
                    gsub2.to_csv(f"{fn_base}.tsv", sep="\t", index=False)

                if args.write_pivots:
                    def piv(metric, suf):
                        p = gsub2.pivot_table(
                            index=["chrom","start0","end0","strand","ZN_transcript_index"],
                            columns="sample", values=metric, aggfunc="first"
                        ).fillna(0).reset_index()
                        p.to_csv(f"{fn_base}_{suf}.tsv", sep="\t", index=False)
                    piv("Nvalid_cov","cov_pivot")
                    piv("frac_modified","frac_pivot")
                    piv("Nmod","Nmod_pivot")

    if args.emit_raw:
        write_all(df_raw, "RAW", write_per_gene=args.write_raw_per_gene)
    if args.emit_filt:
        write_all(df_filt, "FILTERED", write_per_gene=args.write_filtered_per_gene)

    print("[OK] ZN aggregation complete.")

if __name__ == "__main__":
    main()

