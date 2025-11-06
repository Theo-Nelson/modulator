#!/usr/bin/env python3
"""
Aggregate modkit bedMethyl outputs partitioned by ZN or ZT into
per-gene, per-mod *site-level* spreadsheets, with the transcript index column.

v4 highlights (tailored for your current layout):
  • Supports BOTH partition filename styles:
      A) suffix code:  <sample>_<ZT>.bed[.gz]  (where ZT like ALCAM.ALCAM.G1.T2)
      B) numeric siblings:  <N>.bed[.gz] next to a base
         <sample>_<mod>_filtered_mod.bed  (ZN = N)
  • Assigns each site in numeric files to a gene by intersecting with
    the assembler GTF transcript ranges (per chromosome/strand).
  • Derives gene_name/gene_id and (gene_index, transcript_index) from:
        - GTF attributes (preferred: gene_index, transcript_index, zt_label), OR
        - transcript_id like NAME.G#.T# if attributes missing, OR
        - ZT code in filename for style A.
  • Writes one TSV per gene **per modification** with all sites and a
    `transcript_index` column (== ZN for numeric files, == T# for ZT files).
  • Also emits a single combined long table if you want to scan everything.

Outputs (under --out-prefix):
  <out_prefix>_sites_long.tsv                                     # optional global long table
  <out_prefix>.per_gene_sites/<gene_index>_<gene_name>_<gene_id>_<mod>.tsv

Usage example:
  python3 aggregate_modkit_by_ZN_per_gene_v4.py \
    --modkit-dir modkit_out_ZN \
    --gtf fivegenes_readbacked_annot.gtf \
    --summary-tsv fivegenes_readbacked_annot_classification_summary.tsv \
    --out-prefix modkit_by_transcript_ZN \
    --min-cov 5 \
    --verbose

Notes:
  • Files named "ungrouped.bed" and unpartitioned base beds (…_filtered_mod.bed) are ignored.
  • If summary lacks gene fields, names are recovered from GTF or ZT label.
  • If multiple transcripts overlap a site, we choose the gene with the
    GTF interval giving the largest overlap with the site window.
"""

import os, sys, argparse, glob, gzip, re
from collections import defaultdict
from typing import Dict, List, Tuple, Optional

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

# ------------------------- helpers -------------------------

def open_textmaybe_gzip(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")

def is_header_line(line: str) -> bool:
    s = line.strip()
    return (not s) or s.startswith("#") or s.startswith("track") or s.startswith("browser")

def safe_int(x, default=0):
    try:
        return int(x)
    except Exception:
        try:
            return int(float(x))
        except Exception:
            return default

def parse_zt(zt: str) -> Tuple[Optional[str], Optional[str], Optional[int], Optional[int]]:
    """Parse ZT like NAME.G<g>.T<t> with optional leading NAME.NAME.
       Returns (gene_name, gene_id, gene_index, transcript_index)."""
    if not zt:
        return None, None, None, None
    core = zt.split(".")
    # Accept NAME.GENEID.G#.T# (length >=4) or compact G#.T#
    gene_name, gene_id, g_idx, t_idx = None, None, None, None
    # Find last two tokens containing G# and T#
    for i in range(len(core)-1):
        if core[i].startswith("G") and core[i+1].startswith("T"):
            gene_name = core[0] if len(core) >= 1 else None
            gene_id = core[1] if len(core) >= 2 else gene_name
            try:
                g_idx = int(core[i][1:])
            except Exception:
                g_idx = None
            try:
                t_idx = int(core[i+1][1:])
            except Exception:
                t_idx = None
            return gene_name, gene_id, g_idx, t_idx
    # Fallback: try transcript_id style NAME.G#.T#
    m = re.search(r"^(?P<name>[^.]+)\.G(?P<g>\d+)\.T(?P<t>\d+)$", zt)
    if m:
        return m.group("name"), m.group("name"), int(m.group("g")), int(m.group("t"))
    return None, None, None, None

# ---------------------- GTF parsing & index ----------------------

def parse_attrs(attr_field: str) -> Dict[str,str]:
    d = {}
    for m in re.finditer(r'(\S+)\s+"([^"]+)"', attr_field or ""):
        d[m.group(1)] = m.group(2)
    # normalize
    if "gene_name" not in d and "ref_gene_name" in d:
        d["gene_name"] = d["ref_gene_name"]
    if "gene_name" not in d and "gene_id" in d:
        d["gene_name"] = d["gene_id"]
    return d

class TxRec:
    __slots__ = ("chrom","strand","start","end","gene_id","gene_name","gene_index","transcript_index","zt_label")
    def __init__(self, chrom, strand, start, end, attrs):
        self.chrom = chrom; self.strand = strand; self.start = start; self.end = end
        self.gene_id = attrs.get("gene_id") or attrs.get("ref_gene_name")
        self.gene_name = attrs.get("gene_name") or self.gene_id
        self.gene_index = None
        self.transcript_index = None
        self.zt_label = attrs.get("zt_label")
        # Prefer explicit indices
        if attrs.get("gene_index"):
            try: self.gene_index = int(attrs.get("gene_index"))
            except Exception: pass
        if attrs.get("transcript_index"):
            try: self.transcript_index = int(attrs.get("transcript_index"))
            except Exception: pass
        # Derive from transcript_id if needed
        if (self.gene_index is None or self.transcript_index is None) and attrs.get("transcript_id"):
            gname, gid, gi, ti = parse_zt(attrs.get("transcript_id"))
            if gi is not None and self.gene_index is None: self.gene_index = gi
            if ti is not None and self.transcript_index is None: self.transcript_index = ti
            if not self.gene_name and gname: self.gene_name = gname
            if not self.gene_id and gid: self.gene_id = gid
        # Derive from zt_label if still missing
        if (self.gene_index is None or self.transcript_index is None) and self.zt_label:
            gname, gid, gi, ti = parse_zt(self.zt_label)
            if gi is not None and self.gene_index is None: self.gene_index = gi
            if ti is not None and self.transcript_index is None: self.transcript_index = ti
            if not self.gene_name and gname: self.gene_name = gname
            if not self.gene_id and gid: self.gene_id = gid


def load_gtf_transcripts(gtf_path: str) -> List[TxRec]:
    txs: List[TxRec] = []
    with open(gtf_path) as f:
        for line in f:
            if not line or line.startswith("#"): continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9: continue
            chrom, source, feature, start, end, score, strand, frame, attrs = parts
            if feature != "transcript":
                continue
            a = parse_attrs(attrs)
            try:
                txs.append(TxRec(chrom, strand, int(start), int(end), a))
            except Exception:
                continue
    return [t for t in txs if t.gene_index is not None]

# Build simple interval bins per (chrom,strand) for fast overlap
class IntervalBin:
    def __init__(self):
        self.edges: List[int] = []
        self.items: List[Tuple[int,int,Tuple]] = []  # (start,end,key)

    def add(self, start:int, end:int, key:Tuple):
        self.edges.append(start); self.edges.append(end)
        self.items.append((start,end,key))

    def finalize(self):
        self.edges = sorted(set(self.edges))

    def query(self, pos_start:int, pos_end:int) -> List[Tuple[int,int,Tuple]]:
        # naive scan is fine for small five-gene test; could add bisect for scale
        out = []
        for s,e,key in self.items:
            if e < pos_start or s > pos_end: continue
            out.append((s,e,key))
        return out

# --------------------------- I/O ---------------------------

def sample_and_partition_from_path(path: str) -> Tuple[Optional[str], Optional[str], Optional[int]]:
    """Return (sample, ZT_code, ZN_index) based on filename.
       - If name has an underscore and ends with .bed(.gz), assume suffix code style A and extract ZT.
       - If basename is N.bed(.gz), return ZN=N.
    """
    base = os.path.basename(path)
    if base.endswith(".bed.gz"): base = base[:-7]
    elif base.endswith(".bed"): base = base[:-4]
    if base == "ungrouped":
        return None, None, None
    # Style B: numeric
    if re.fullmatch(r"\d+", base):
        return None, None, int(base)
    # Style A: suffix after last underscore
    if "_" in base:
        sample, code = base.rsplit("_", 1)
        return sample, code, None
    return None, None, None


def read_bed_sites(path: str) -> List[dict]:
    rows = []
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
            # coerce ints/floats
            d["start0"] = safe_int(d["start0"]) ; d["end0"] = safe_int(d["end0"]) ;
            for k in ["Nvalid_cov","Nmod","Ncanonical","Nother_mod","Ndelete","Nfail","Ndiff","Nnocall"]:
                d[k] = safe_int(d[k])
            try:
                d["frac_modified"] = float(d["frac_modified"]) if d["frac_modified"] != "." else 0.0
            except Exception:
                d["frac_modified"] = 0.0
            rows.append(d)
    return rows

# --------------------------- main ---------------------------

def parse_args():
    ap = argparse.ArgumentParser(description="Aggregate ZN/ZT modkit site beds into per-gene per-mod TSVs")
    ap.add_argument("--modkit-dir", required=True)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--summary-tsv", help="Optional assembler summary for zt_label recovery")
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--min-cov", type=int, default=0, help="Zero-out frac_modified if Nvalid_cov < MIN_COV")
    ap.add_argument("--verbose", action="store_true")
    return ap.parse_args()


def main():
    args = parse_args()

    # Load GTF transcript ranges
    txs = load_gtf_transcripts(args.gtf)
    by_cs: Dict[Tuple[str,str], IntervalBin] = defaultdict(IntervalBin)
    key_to_meta = {}  # key -> (gene_index, gene_name, gene_id)
    for t in txs:
        key = (t.gene_index, t.gene_name or t.gene_id or "NA", t.gene_id or t.gene_name or "NA")
        by_cs[(t.chrom, t.strand)].add(t.start, t.end, key)
        key_to_meta[key] = key
    for b in by_cs.values():
        b.finalize()
    if args.verbose:
        print(f"[info] built interval index for {sum(len(b.items) for b in by_cs.values())} transcript ranges", file=sys.stderr)

    out_dir = f"{args.out_prefix}.per_gene_sites"
    os.makedirs(out_dir, exist_ok=True)

    # Optional summary join for Style A ZT -> indices
    zt_to_indices = {}
    if args.summary_tsv and os.path.exists(args.summary_tsv):
        try:
            df = pd.read_csv(args.summary_tsv, sep="\t", comment=None, dtype=str)
            # normalize header
            cols = {c: c.lstrip('#') for c in df.columns}
            df.rename(columns=cols, inplace=True)
            if 'zt_label' in df.columns and 'code' not in df.columns:
                df['code'] = df['zt_label']
            if 'code' in df.columns:
                for _, r in df.iterrows():
                    zt = str(r.get('code') or '').strip()
                    gi = r.get('gene_index'); ti = r.get('transcript_index')
                    if zt:
                        try: gi_i = int(gi) if pd.notna(gi) else None
                        except Exception: gi_i = None
                        try: ti_i = int(ti) if pd.notna(ti) else None
                        except Exception: ti_i = None
                        zt_to_indices[zt] = (gi_i, ti_i)
        except Exception as e:
            if args.verbose:
                print(f"[warn] could not parse summary: {e}", file=sys.stderr)

    beds = sorted(glob.glob(os.path.join(args.modkit_dir, "*.bed"))) + \
           sorted(glob.glob(os.path.join(args.modkit_dir, "*.bed.gz")))
    if not beds:
        sys.exit(f"No .bed/.bed.gz files found in {args.modkit_dir}")

    long_rows: List[dict] = []
    per_gene_mod_rows: Dict[Tuple[int,str,str,str], List[dict]] = defaultdict(list)
    # key: (gene_index, gene_name, gene_id, mod_code)

    for bed in beds:
        base = os.path.basename(bed)
        if base.startswith("ungrouped"):  # ignore
            if args.verbose:
                print(f"[skip] {base} (ungrouped)", file=sys.stderr)
            continue

        sample, code, zn = sample_and_partition_from_path(bed)
        rows = read_bed_sites(bed)
        if not rows:
            if args.verbose:
                print(f"[skip] {base} (no data)", file=sys.stderr)
            continue

        if code:  # Style A: ZT label in filename
            gname, gid, gidx, tidx = parse_zt(code)
            if (gidx is None or tidx is None) and code in zt_to_indices:
                gidx2, tidx2 = zt_to_indices.get(code, (None,None))
                gidx = gidx if gidx is not None else gidx2
                tidx = tidx if tidx is not None else tidx2
            # fallbacks
            if gidx is None: gidx = -1
            if not gname: gname = gid or "NA"
            if not gid: gid = gname or "NA"
            for d in rows:
                if args.min_cov and d["Nvalid_cov"] < args.min_cov:
                    d["frac_modified"] = 0.0
                rec = {
                    **d,
                    "sample": sample or "",
                    "gene_index": gidx,
                    "gene_name": gname,
                    "gene_id": gid,
                    "transcript_index": tidx,
                    "partition": code,
                }
                long_rows.append(rec)
                key = (gidx, gname, gid, d["mod_code"])
                per_gene_mod_rows[key].append(rec)

        elif zn is not None:  # Style B: numeric ZN -> need to assign to a gene by interval overlap
            for d in rows:
                if args.min_cov and d["Nvalid_cov"] < args.min_cov:
                    d["frac_modified"] = 0.0
                # find gene by overlap with transcript spans on (chrom,strand)
                cs = (d["chrom"], d["strand"])
                start = d["start0"]; end = d["end0"]
                best = None; best_ov = -1
                if cs in by_cs:
                    for s,e,key in by_cs[cs].query(start, end):
                        ov = min(e, end) - max(s, start) + 1
                        if ov > best_ov:
                            best_ov = ov; best = key
                if best is None:
                    # could not assign — skip or send to NA bucket
                    if args.verbose:
                        print(f"[warn] unassigned site {d['chrom']}:{start}-{end}{d['strand']} in {base}", file=sys.stderr)
                    continue
                gidx, gname, gid = best
                rec = {
                    **d,
                    "sample": sample or "",
                    "gene_index": gidx,
                    "gene_name": gname,
                    "gene_id": gid,
                    "transcript_index": zn,
                    "partition": str(zn),
                }
                long_rows.append(rec)
                key = (gidx, gname, gid, d["mod_code"])
                per_gene_mod_rows[key].append(rec)
        else:
            if args.verbose:
                print(f"[skip] {base} (unparseable partition)", file=sys.stderr)
            continue

    if not long_rows:
        sys.exit("No site rows collected — check inputs.")

    # Write global long table
    out_long = f"{args.out_prefix}_sites_long.tsv"
    pd.DataFrame(long_rows).to_csv(out_long, sep='\t', index=False)

    # Write per-gene per-mod tables
    for (gidx, gname, gid, mod), rows in per_gene_mod_rows.items():
        df = pd.DataFrame(rows)
        # stable ordering
        df = df[[
            "chrom","start0","end0","strand","mod_code","sample",
            "gene_index","gene_name","gene_id","transcript_index",
            "Nvalid_cov","Nmod","Ncanonical","Nother_mod","Ndelete","Nfail","Ndiff","Nnocall","frac_modified"
        ]]
        fname = f"{gidx}_{gname}_{gid}_{mod}.tsv".replace("/","_")
        df.to_csv(os.path.join(out_dir, fname), sep='\t', index=False)

    if args.verbose:
        print(f"[ok] wrote {out_long}", file=sys.stderr)
        print(f"[ok] wrote per-gene site tables in {out_dir}", file=sys.stderr)

    print("[OK] ZN/ZT per-gene site aggregation complete.")

if __name__ == "__main__":
    main()

