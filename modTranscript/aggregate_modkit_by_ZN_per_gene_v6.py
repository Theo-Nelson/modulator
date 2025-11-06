#!/usr/bin/env python3
"""
Aggregate modkit outputs written with --partition-tag ZN.

v6:
  • Walks subfolders AND handles the modkit layout where
    <sample>_<MOD>_filtered_mod.bed/ is a DIRECTORY containing ZN parts:
      1.bed, 2.bed, ..., ungrouped.bed
  • Extracts:
      sample  = parent folder name above the <sample>_<MOD>_filtered_mod.bed dir
      mod     = from '<MOD>' in '<sample>_<MOD>_filtered_mod.bed'
      ZN      = transcript_index from inner filename (e.g., '3.bed' -> 3)
  • Maps each site to a gene using the provided GTF (ReadBacked GTF is fine).
  • Writes:
      <out_prefix>_sites_long.tsv
      <out_prefix>.per_gene_per_mod/<geneIndex>_<geneName>_<mod>.tsv
"""

import os, sys, re, gzip, argparse
from collections import defaultdict
from typing import Dict, List, Tuple, Iterable

try:
    import pandas as pd
except ImportError:
    sys.exit("This script needs pandas. (micromamba install pandas)")

# columns in modkit bedMethyl pileup
BED_COLS = [
    "chrom","start0","end0","mod_code","score","strand",
    "start0_compat","end0_compat","rgb",
    "Nvalid_cov","frac_modified",
    "Nmod","Ncanonical","Nother_mod",
    "Ndelete","Nfail","Ndiff","Nnocall",
]

def parse_args():
    ap = argparse.ArgumentParser(description="Aggregate modkit (partitioned by ZN) into per-gene per-mod tables.")
    ap.add_argument("--modkit-dir", required=True, help="Root of modkit outputs")
    ap.add_argument("--gtf", required=True, help="ReadBacked GTF from assembler (has gene_index/transcript_index)")
    ap.add_argument("--out-prefix", required=True, help="Output prefix")
    ap.add_argument("--min-cov", type=int, default=0, help="Zero frac_modified where Nvalid_cov < MIN")
    ap.add_argument("--verbose", action="store_true")
    return ap.parse_args()

def open_text(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")

def is_header(line: str) -> bool:
    s = line.strip()
    return (not s) or s.startswith("#") or s.startswith("track") or s.startswith("browser")

def read_bed_rows(path: str) -> Iterable[Dict[str,str]]:
    with open_text(path) as f:
        for ln in f:
            if is_header(ln): 
                continue
            parts = ln.rstrip("\n").split("\t")
            if len(parts) < 18:
                continue
            d = dict(zip(BED_COLS, parts[:18]))
            yield d

# ---------------- GTF mapping (gene spans) ----------------

def parse_attrs(attr: str) -> Dict[str,str]:
    out = {}
    for kv in re.finditer(r'(\S+)\s+"([^"]*)"', attr):
        out[kv.group(1)] = kv.group(2)
    return out

def load_gene_spans(gtf_path: str, verbose=False):
    """
    Build per (chrom,strand) list of gene spans:
      [(start,end,gene_id,gene_name,gene_index), ...]
    We derive gene_name from ref_gene_name if present, else gene_id.
    We take transcript rows to collect gene spans.
    """
    per_key = defaultdict(list)  # (chrom,strand) -> list of (start,end,meta)
    by_gene = {}  # gene_id -> [minStart, maxEnd, name, index, chrom, strand]

    with open_text(gtf_path) as f:
        for ln in f:
            if ln.startswith("#"): 
                continue
            parts = ln.rstrip("\n").split("\t")
            if len(parts) < 9: 
                continue
            chrom, _, feature, start, end, _, strand, _, attrs = parts
            if feature != "transcript":
                continue
            a = parse_attrs(attrs)
            gene_id  = a.get("gene_id", "")
            gene_nm  = a.get("ref_gene_name", gene_id)
            gidx     = a.get("gene_index", "")
            s = int(start); e = int(end)
            if gene_id not in by_gene:
                by_gene[gene_id] = [s, e, gene_nm, gidx, chrom, strand]
            else:
                by_gene[gene_id][0] = min(by_gene[gene_id][0], s)
                by_gene[gene_id][1] = max(by_gene[gene_id][1], e)

    for gid, (s,e,gnm,gidx,chrom,strand) in by_gene.items():
        per_key[(chrom,strand)].append( (s,e, {"gene_id":gid, "gene_name":gnm, "gene_index":gidx}) )

    # sort for binary search
    for key in per_key:
        per_key[key].sort(key=lambda x: (x[0], x[1]))

    if verbose:
        n = sum(len(v) for v in per_key.values())
        print(f"[gtf] loaded {n} gene spans", file=sys.stderr)
    return per_key

def map_site_to_gene(chrom: str, pos1: int, strand: str, spans) -> Dict[str,str]:
    lst = spans.get((chrom, strand), [])
    # binary search by start
    lo, hi = 0, len(lst)-1
    hit = None
    while lo <= hi:
        mid = (lo+hi)//2
        s,e,meta = lst[mid]
        if pos1 < s:
            hi = mid - 1
        elif pos1 > e:
            lo = mid + 1
        else:
            hit = meta; break
    return hit  # dict or None

# ---------------- path parsing ----------------

MODDIR_RE = re.compile(r'(.+?)_([A-Za-z0-9]+)_filtered_mod\.bed$')

def enumerate_partitions(modkit_root: str, verbose=False):
    """
    Yields tuples for every ZN partition bed found:
      (bed_path, sample, mod_code, zn_index)
    Supports two layouts:
      A)  <root>/<sample>/<sample>_<MOD>_filtered_mod.bed/<ZN>.bed
      B)  <root>/**/<sample>_<MOD>_filtered_mod.bed (flat file; then zn_index=None)
    """
    out = []
    for dirpath, dirnames, filenames in os.walk(modkit_root):
        base = os.path.basename(dirpath)
        m = MODDIR_RE.match(base)
        if m:
            sample = m.group(1)
            mod = m.group(2)
            # This is a partition DIRECTORY; inner files are ZN parts
            inner = [f for f in os.listdir(dirpath) if f.endswith(".bed") or f.endswith(".bed.gz")]
            if not inner and verbose:
                print(f"[warn] empty partition dir: {dirpath}", file=sys.stderr)
            for f in inner:
                zn = os.path.splitext(os.path.basename(f))[0]
                zn_idx = int(zn) if zn.isdigit() else None  # 'ungrouped' -> None
                out.append( (os.path.join(dirpath, f), sample, mod, zn_idx) )
        else:
            # Also include flat .bed files (no partitions)
            for f in filenames:
                if not (f.endswith(".bed") or f.endswith(".bed.gz")):
                    continue
                # try to get sample/mod from filename itself
                stem = f[:-7] if f.endswith(".bed.gz") else f[:-4]
                mm = MODDIR_RE.match(stem)
                if mm:
                    out.append( (os.path.join(dirpath, f), mm.group(1), mm.group(2), None) )
    if verbose:
        print(f"[scan] found {len(out)} bed partitions", file=sys.stderr)
    return out

# ---------------- main aggregation ----------------

def safe_int(x):
    try:
        return int(x)
    except Exception:
        try:
            return int(float(x))
        except Exception:
            return 0

def main():
    args = parse_args()
    spans = load_gene_spans(args.gtf, verbose=args.verbose)

    parts = enumerate_partitions(args.modkit_dir, verbose=args.verbose)
    if not parts:
        sys.exit(f"No .bed/.bed.gz partitions found under {args.modkit_dir}")

    rows: List[dict] = []
    for bed_path, sample, mod, zn in parts:
        for rec in read_bed_rows(bed_path):
            chrom = rec["chrom"]
            strand = rec["strand"]
            start0 = safe_int(rec["start0"])
            end0   = safe_int(rec["end0"])
            pos1   = start0 + 1  # 1-based
            gene = map_site_to_gene(chrom, pos1, strand, spans)
            if not gene:
                continue  # skip sites outside any assembled gene on this strand
            Ncov = safe_int(rec["Nvalid_cov"])
            Nmod = safe_int(rec["Nmod"])
            frac = float(rec.get("frac_modified", 0))
            if args.min_cov and Ncov < args.min_cov:
                frac = 0.0
            rows.append({
                "sample": sample,
                "mod_code": mod,
                "ZN_transcript_index": zn,   # None for ungrouped/flat
                "chrom": chrom,
                "start0": start0,
                "end0": end0,
                "strand": strand,
                "Nvalid_cov": Ncov,
                "Nmod": Nmod,
                "frac_modified": round(frac, 6),
                # gene meta from GTF
                "gene_index": gene.get("gene_index", ""),
                "gene_name": gene.get("gene_name", gene.get("gene_id","")),
                "gene_id": gene.get("gene_id",""),
            })

    if not rows:
        sys.exit("No rows produced. Check GTF mapping and inputs.")

    df = pd.DataFrame(rows)
    out_long = f"{args.out_prefix}_sites_long.tsv"
    df.to_csv(out_long, sep="\t", index=False)

    # per-gene per-mod spreadsheets
    outdir = f"{args.out_prefix}.per_gene_per_mod"
    os.makedirs(outdir, exist_ok=True)
    for (gidx, gname, mod), sub in df.groupby(["gene_index","gene_name","mod_code"]):
        tag = f"{gidx}_{gname}_{mod}".replace("/", "_")
        sub.sort_values(["ZN_transcript_index","sample","chrom","start0","end0"], inplace=True)
        sub.to_csv(os.path.join(outdir, f"{tag}.tsv"), sep="\t", index=False)

    if args.verbose:
        print(f"[ok] wrote {out_long}", file=sys.stderr)
        print(f"[ok] wrote per-gene per-mod tables under {outdir}", file=sys.stderr)
    print("[OK] ZN aggregation complete.")

if __name__ == "__main__":
    main()

