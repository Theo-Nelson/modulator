#!/usr/bin/env python3
"""
Aggregate modkit ZN-partitioned bedMethyl outputs into per-gene, per-mod tables,
with duplicate-site collapsing and explicit transcript index (ZN).

What it handles:
- Directory layout like:
    modkit_out_ZN/<sample>/
      <sample>_<mod>_filtered_mod.bed            # flat per-mod file (ignored)
      <sample>_<mod>_filtered_mod_bed_log
      1.bed 2.bed 3.bed ... ungrouped.bed        # ZN partitions (we ONLY read N.bed)
- Each numbered file contains bedMethyl rows with the mod code in column 4.
- Deduplication: sums counts for identical keys:
    (sample, mod_code, ZN_transcript_index, chrom, start0, end0, strand)
- Gene mapping: intersects sites to GTF features to annotate gene_name/gene_id.
  (If multiple genes overlap, picks the one with largest overlap span; tie → first.)
- Outputs:
  1) <out_prefix>_sites_long.tsv                 (all sites, annotated, deduped)
  2) Per gene × mod:
       <out_prefix>__<geneName>__<mod>.tsv
       <out_prefix>__<geneName>__<mod>_cov_pivot.tsv
       <out_prefix>__<geneName>__<mod>_frac_pivot.tsv
       <out_prefix>__<geneName>__<mod>_Nmod_pivot.tsv

Minimal dependencies: pandas
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
    ap = argparse.ArgumentParser(description="Aggregate ZN-partitioned modkit outputs per gene/mod")
    ap.add_argument("--modkit-dir", required=True, help="Parent dir with per-sample subdirs containing numbered ZN .bed files")
    ap.add_argument("--gtf", required=True, help="Assembler GTF (with gene coordinates). Exon or transcript features work.")
    ap.add_argument("--out-prefix", required=True, help="Prefix for outputs")
    ap.add_argument("--min-cov", type=int, default=0, help="Zero frac_modified if Nvalid_cov < MIN_COV (row kept)")
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

# --- GTF interval indexing (simple, robust) ---
Interval = namedtuple("Interval", ["start","end","gene_id","gene_name","strand"])

def load_gene_intervals_from_gtf(gtf_path:str, verbose=False)->Dict[Tuple[str,str], List[Interval]]:
    """
    Build per-(chrom,strand) interval list from any 'exon' or 'transcript' (fallback 'gene') rows.
    We union per gene by simple min(start), max(end) to get a coarse gene span, good for site→gene mapping.
    """
    gene_bounds: Dict[Tuple[str,str,str], Tuple[int,int]] = {}  # (chrom,strand,gene_id)->(minS,maxE)
    gene_name_map: Dict[str,str] = {}

    with open_text(gtf_path) as f:
        for ln in f:
            if ln.startswith("#") or not ln.strip():
                continue
            parts = ln.rstrip("\n").split("\t")
            if len(parts) < 9: continue
            chrom, source, feature, start, end, score, strand, frame, attrs = parts
            if feature not in ("exon","transcript","gene"):
                continue
            # parse attrs
            a = {}
            for kv in re.finditer(r'(\S+)\s+"([^"]*)"', attrs):
                a[kv.group(1)] = kv.group(2)
            gene_id = a.get("gene_id") or a.get("gtf_gene_id") or a.get("gene") or ""
            gene_name = a.get("ref_gene_name") or a.get("gene_name") or a.get("gtf_gene_name") or gene_id
            if not gene_id:
                continue
            s = int(start)
            e = int(end)
            key = (chrom, strand, gene_id)
            if key not in gene_bounds:
                gene_bounds[key] = (s, e)
            else:
                mn, mx = gene_bounds[key]
                gene_bounds[key] = (min(mn, s), max(mx, e))
            # keep last-seen name
            gene_name_map[gene_id] = gene_name

    by_cs: Dict[Tuple[str,str], List[Interval]] = defaultdict(list)
    for (chrom, strand, gid), (s,e) in gene_bounds.items():
        gname = gene_name_map.get(gid, gid)
        by_cs[(chrom, strand)].append(Interval(s, e, gid, gname, strand))
    # sort intervals for fast scan
    for k in by_cs:
        by_cs[k].sort(key=lambda iv: (iv.start, iv.end))
    if verbose:
        print(f"[info] loaded {sum(len(v) for v in by_cs.values())} gene intervals from {gtf_path}", file=sys.stderr)
    return by_cs

def assign_gene(chrom:str, pos_start:int, pos_end:int, strand:str, gene_index:Dict[Tuple[str,str], List[Interval]]):
    """
    Return (gene_id, gene_name) by overlap; choose max overlap; tie → first.
    """
    ivs = gene_index.get((chrom, strand), [])
    best = None
    best_ov = -1
    for iv in ivs:
        # early break if iv.start already beyond site
        if iv.start > pos_end: break
        if iv.end < pos_start: continue
        ov = min(iv.end, pos_end) - max(iv.start, pos_start) + 1
        if ov > best_ov:
            best_ov = ov
            best = iv
    if best:
        return best.gene_id, best.gene_name
    # try opposite strand if stranded annotation is sparse
    ivs2 = gene_index.get((chrom, "+" if strand == "-" else "-"), [])
    for iv in ivs2:
        if iv.start > pos_end: break
        if iv.end < pos_start: continue
        return iv.gene_id, iv.gene_name
    return "", ""

# --- core readers ---
def iter_numbered_beds(modkit_dir:str)->List[Tuple[str, str, str, int]]:
    """
    Yields tuples: (sample_dir, sample_name, bed_path, ZN_index) for files like '.../<sample>/*/<N>.bed'
    We also support files directly under the sample dir ('<sample>/1.bed').
    Skips 'ungrouped.bed' and non-numeric basenames. Ignores flat '*_filtered_mod.bed'.
    """
    out = []
    for root, dirs, files in os.walk(modkit_dir):
        # sample_name is the immediate folder under modkit_dir
        rel = os.path.relpath(root, modkit_dir)
        if rel == ".":
            # top level; expect per-sample subdirs here—skip
            continue
        sample_name = rel.split(os.sep)[0]
        for fname in files:
            if fname.endswith("_filtered_mod.bed") or fname.endswith("_filtered_mod.bed.gz"):
                continue  # flat per-mod file → ignore for ZN aggregation
            base = fname
            if base.endswith(".gz"): base = base[:-3]
            if base.lower() == "ungrouped.bed":
                continue
            m = re.fullmatch(r"(\d+)\.bed", base)
            if not m:
                continue
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
    # coerce numerics
    d["start0"] = safe_int(d["start0"])
    d["end0"]   = safe_int(d["end0"])
    for k in ["Nvalid_cov","Nmod","Ncanonical","Nother_mod","Ndelete","Nfail","Ndiff","Nnocall"]:
        d[k] = safe_int(d[k])
    d["frac_modified"] = float(d["frac_modified"]) if d["frac_modified"] not in ("", ".") else 0.0
    return d

def main():
    args = parse_args()
    beds = iter_numbered_beds(args.modkit_dir)
    if not beds:
        sys.exit(f"No numbered ZN partition files found under {args.modkit_dir}")

    gene_index = load_gene_intervals_from_gtf(args.gtf, verbose=args.verbose)

    rows = []
    for sample_dir, sample_name, bed_path, zn in beds:
        with open_text(bed_path) as f:
            for ln in f:
                if is_header(ln): continue
                rec = parse_bed_line(ln)
                if not rec: continue
                gid, gname = assign_gene(rec["chrom"], rec["start0"], rec["end0"], rec["strand"], gene_index)
                # min-cov rule applied after collapsing; store raw now
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

    # Collapse duplicates within same genomic+partition key
    key = ["sample","mod_code","ZN_transcript_index","chrom","start0","end0","strand","gene_id","gene_name"]
    sumcols = ["Nvalid_cov","Nmod","Ncanonical","Nother_mod","Ndelete","Nfail","Ndiff","Nnocall"]
    pre = len(df)
    df = df.groupby(key, as_index=False)[sumcols].sum()
    if args.verbose:
        print(f"[info] collapsed {pre} → {len(df)} rows", file=sys.stderr)

    # Compute frac_modified with min-cov logic
    df["frac_modified"] = (df["Nmod"] / df["Nvalid_cov"].where(df["Nvalid_cov"]>0, 1)).fillna(0.0)
    if args.min_cov:
        df.loc[df["Nvalid_cov"] < args.min_cov, "frac_modified"] = 0.0
    df["frac_modified"] = df["frac_modified"].round(6)

    # Write master long table
    out_long = f"{args.out_prefix}_sites_long.tsv"
    df.to_csv(out_long, sep="\t", index=False)
    if args.verbose:
        print(f"[ok] wrote {out_long}", file=sys.stderr)

    # Per gene × mod outputs
    os.makedirs(f"{args.out_prefix}__per_gene_mod", exist_ok=True)
    for (gname, mod), sub in df.groupby(["gene_name","mod_code"], dropna=False):
        safe_g = re.sub(r"[^A-Za-z0-9._+-]", "_", gname if gname else "NA")
        fn_base = f"{args.out_prefix}__per_gene_mod/{args.out_prefix}__{safe_g}__{mod}"
        # tidy columns and sort
        cols = ["gene_name","gene_id","mod_code","chrom","start0","end0","strand",
                "ZN_transcript_index","sample","Nvalid_cov","Nmod","Ncanonical",
                "Nother_mod","Ndelete","Nfail","Ndiff","Nnocall","frac_modified"]
        sub = sub[cols].sort_values(["chrom","start0","ZN_transcript_index","sample"])
        sub.to_csv(f"{fn_base}.tsv", sep="\t", index=False)

        # pivots by sample
        def piv(metric, suf):
            p = sub.pivot_table(index=["chrom","start0","end0","strand","ZN_transcript_index"],
                                columns="sample", values=metric, aggfunc="first").fillna(0).reset_index()
            p.to_csv(f"{fn_base}_{suf}.tsv", sep="\t", index=False)
        piv("Nvalid_cov","cov_pivot")
        piv("frac_modified","frac_pivot")
        piv("Nmod","Nmod_pivot")

    print("[OK] ZN per-gene per-mod aggregation complete.")

if __name__ == "__main__":
    main()

