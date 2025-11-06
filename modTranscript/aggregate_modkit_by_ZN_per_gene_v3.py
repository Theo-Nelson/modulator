#!/usr/bin/env python3
"""
Aggregate modkit bedMethyl outputs produced with `--partition-tag ZN` and emit:
  1) <out_prefix>_ZN_long.tsv  (all rows, with gene_index/transcript_index and gene name/id)
  2) <out_prefix>.per_gene/<Gidx>_<GeneName>_<GeneID>.tsv (tidy per-gene)
  3) optional per-gene pivots (frac/cov/Nmod) with --write-gene-pivots

v3 (this file) vs v2:
  * Understands BOTH filename styles that modkit can emit when partitioning:
      A) "<sample>_<mod>_filtered_mod_<ZN>.bed[.gz]"  (expected)
      B) Numeric siblings: "<ZN>.bed" (and "ungrouped.bed") written next to
         the base file "<sample>_<mod>_filtered_mod.bed". We infer:
            - sample := parent directory name (your workflow uses per-sample folders)
            - mod    := first data row's mod_code in the file
            - zn     := integer filename stem
    Unpartitioned base files are skipped by design.
  * More defensive parsing; clearer verbose output.

Requires: pandas
"""

import os, sys, argparse, glob, gzip, re
from collections import defaultdict
from typing import Dict, Tuple, List, Iterable, Optional

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

# ----------------- CLI -----------------

def parse_args():
    ap = argparse.ArgumentParser(description="Aggregate ZN-partitioned modkit beds into per-gene tables (v3)")
    ap.add_argument("--modkit-dir", required=True,
                    help="Directory tree containing per-ZN .bed/.bed.gz files (either *_filtered_mod_<ZN>.bed or <ZN>.bed)")
    ap.add_argument("--summary-tsv", required=False,
                    help="Assembler classification summary TSV (optional, improves gene labels)")
    ap.add_argument("--gtf", required=True,
                    help="Assembler v8 GTF (contains exon blocks and gene/transcript indices)")
    ap.add_argument("--out-prefix", required=True,
                    help="Prefix for outputs (a '.per_gene' folder will be created alongside)")
    ap.add_argument("--min-cov", type=int, default=0,
                    help="If >0, zero frac_modified where Nvalid_cov < MIN_COV (row kept)")
    ap.add_argument("--write-gene-pivots", action="store_true",
                    help="Also emit wide pivot tables per gene (frac, cov, Nmod)")
    ap.add_argument("--verbose", action="store_true")
    return ap.parse_args()

# ----------------- utils -----------------

def is_header_line(line: str) -> bool:
    s = line.strip()
    return (not s) or s.startswith("#") or s.startswith("track") or s.startswith("browser")

def open_text(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")

def safe_int(x, default=0):
    try:
        return int(x)
    except Exception:
        try:
            return int(float(x))
        except Exception:
            return default

# ----------------- GTF -----------------

def parse_gtf_attrs(attr_field: str) -> Dict[str, str]:
    d: Dict[str, str] = {}
    for m in re.finditer(r'(\S+)\s+"([^"]+)"', attr_field):
        d[m.group(1)] = m.group(2)
    return d

def parse_gtf(gtf_path: str) -> List[dict]:
    rows = []
    with open(gtf_path) as f:
        for ln in f:
            if not ln or ln.startswith("#"): continue
            parts = ln.rstrip("\n").split("\t")
            if len(parts) < 9: continue
            chrom, source, feature, start, end, score, strand, frame, attrs = parts
            if feature not in ("transcript", "exon"): continue
            start, end = int(start), int(end)
            a = parse_gtf_attrs(attrs)
            rows.append({
                "feature": feature,
                "chrom": chrom,
                "start": start,
                "end": end,
                "strand": strand,
                **a
            })
    return rows

def build_gene_tx_interval_index(gtf_rows: List[dict]):
    per_tid = {}
    for r in gtf_rows:
        if r["feature"] == "transcript":
            tid = r.get("transcript_id", f"{r.get('gene_id','NA')}.{r.get('gene_index','NA')}.{r.get('transcript_index','NA')}")
            per_tid[tid] = {
                "chrom": r["chrom"],
                "strand": r["strand"],
                "start": r["start"],
                "end": r["end"],
                "gene_index": safe_int(r.get("gene_index"), -1),
                "transcript_index": safe_int(r.get("transcript_index"), -1),
                "gene_id": r.get("gene_id", "NA"),
                "gene_name": r.get("ref_gene_name", r.get("gene_id", "NA")),
                "zt_label": r.get("zt_label", ""),
                "exons": []
            }
        elif r["feature"] == "exon":
            tid = r.get("transcript_id")
            if tid in per_tid:
                per_tid[tid]["exons"].append((r["start"], r["end"]))

    by_key = defaultdict(list)
    for tid, info in per_tid.items():
        if info["gene_index"] < 0 or info["transcript_index"] < 0:
            continue
        chrom = info["chrom"]; strand = info["strand"]
        gidx = info["gene_index"]; txi = info["transcript_index"]
        exons = sorted(info["exons"]) or [(info["start"], info["end"])]
        tes_1based = info["end"] if strand == "+" else info["start"]
        by_key[(chrom, strand, gidx, txi)].append({
            "tid": tid, "exons": exons, "tes": tes_1based,
            "gene_id": info["gene_id"], "gene_name": info["gene_name"],
            "zt_label": info["zt_label"],
        })
    return by_key

# ----------------- BED parsing -----------------

def read_bed_rows(path: str):
    with open_text(path) as f:
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
            d["start0"], d["end0"] = int(d["start0"]), int(d["end0"])
            for k in ("Nvalid_cov","Nmod","Ncanonical","Nother_mod","Ndelete","Nfail","Ndiff","Nnocall"):
                d[k] = safe_int(d[k], 0)
            try:
                d["frac_modified"] = float(d["frac_modified"])
            except Exception:
                d["frac_modified"] = 0.0
            yield d

# ----------------- ZN filename parsing -----------------

def parse_style_A(path: str) -> Tuple[Optional[str], Optional[str], Optional[int]]:
    """<sample>_<mod>_filtered_mod_<ZN>.bed[.gz]"""
    base = os.path.basename(path)
    if base.endswith(".bed.gz"): base = base[:-7]
    elif base.endswith(".bed"): base = base[:-4]
    if "_filtered_mod_" not in base:
        return None, None, None
    left, zn = base.split("_filtered_mod_", 1)
    zn_i = safe_int(zn, None)
    if zn_i is None:
        return None, None, None
    if "_" not in left:
        return None, None, None
    sample, mod = left.rsplit("_", 1)
    return sample, mod, zn_i


def infer_mod_from_file(path: str) -> Optional[str]:
    for row in read_bed_rows(path):
        mod = row.get("mod_code")
        if mod:
            return mod
    return None


def parse_style_B(path: str) -> Tuple[Optional[str], Optional[str], Optional[int]]:
    """Numeric sibling files: <ZN>.bed (and ungrouped.bed) next to base *_filtered_mod.bed.
       We infer sample from parent dir name; mod from first row's mod_code; ZN from filename stem.
    """
    base = os.path.basename(path)
    if base.endswith(".bed.gz"): base = base[:-7]
    elif base.endswith(".bed"): base = base[:-4]
    if base == "ungrouped":
        return None, None, None
    if not base.isdigit():
        return None, None, None
    zn_i = safe_int(base, None)
    if zn_i is None:
        return None, None, None
    sample = os.path.basename(os.path.dirname(path))  # parent directory name
    mod = infer_mod_from_file(path)
    return sample, mod, zn_i


def parse_sample_mod_zn(path: str):
    a = parse_style_A(path)
    if a != (None, None, None):
        return a
    b = parse_style_B(path)
    if b != (None, None, None):
        return b
    return (None, None, None)

# ----------------- Summary loader -----------------

def robust_load_summary(path: str, verbose=False) -> pd.DataFrame:
    if not path or not os.path.exists(path):
        return pd.DataFrame()
    lines = []
    with open(path) as f:
        for ln in f:
            if ln.strip() == "":
                continue
            lines.append(ln.rstrip("\n"))
    header = None; header_idx = None
    for i, ln in enumerate(lines):
        h = ln.lstrip("#")
        if "\t" in h and any(c.strip() == "code" for c in h.split("\t")):
            header_idx = i
            header = [c.strip() for c in h.split("\t")]
            break
    if header is None:
        header_idx = 0
        header = [c.strip() for c in lines[0].lstrip("#").split("\t")]
    rows = []
    expected = len(header)
    for ln in lines[header_idx+1:]:
        if ln.startswith("#"): continue
        parts = ln.split("\t")
        if len(parts) < expected:
            parts += [""] * (expected - len(parts))
        elif len(parts) > expected:
            parts = parts[:expected]
        rows.append({header[j]: parts[j] for j in range(expected)})
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
    return df

# ----------------- interval helpers -----------------

def interval_overlap_len(a_start: int, a_end: int, b_start: int, b_end: int) -> int:
    lo = max(a_start, b_start)
    hi = min(a_end, b_end)
    return max(0, hi - lo + 1)


def any_exon_overlap(exons: List[Tuple[int,int]], start0: int, end0: int) -> int:
    s, e = start0 + 1, end0
    best = 0
    for xs, xe in exons:
        best = max(best, interval_overlap_len(xs, xe, s, e))
        if best and best >= (e - s + 1):
            break
    return best

# ----------------- main -----------------

def main():
    args = parse_args()

    # Build GTF interval index
    gtf_rows = parse_gtf(args.gtf)
    idx = build_gene_tx_interval_index(gtf_rows)
    if args.verbose:
        print(f"[info] built interval index for {len(idx)} (chrom,strand,gene_index,transcript_index) keys", file=sys.stderr)

    # Summary (optional) → gene labels
    summ = robust_load_summary(args.summary_tsv, verbose=args.verbose) if args.summary_tsv else pd.DataFrame()

    # gene_index → names (summary preferred, else GTF)
    if (not summ.empty) and set(["gene_index","gtf_gene_name","gtf_gene_id"]).issubset(summ.columns):
        gene_meta = (summ.dropna(subset=["gene_index"]).sort_values(["gene_index"]) \
                        [["gene_index","gtf_gene_name","gtf_gene_id"]] \
                        .drop_duplicates("gene_index").set_index("gene_index").to_dict(orient="index"))
    else:
        gene_meta = {}
        for (c,s,gidx,txi), tlist in idx.items():
            if gidx not in gene_meta and tlist:
                t0 = tlist[0]
                gene_meta[gidx] = {"gtf_gene_name": t0.get("gene_name","NA"), "gtf_gene_id": t0.get("gene_id","NA")}

    # Find all bed files
    beds = sorted(glob.glob(os.path.join(args.modkit_dir, "**", "*.bed"), recursive=True)) + \
           sorted(glob.glob(os.path.join(args.modkit_dir, "**", "*.bed.gz"), recursive=True))
    if not beds:
        sys.exit(f"No .bed/.bed.gz files found under {args.modkit_dir}")

    rows_long = []
    per_gene_accum: Dict[Tuple[int,int,str,str], Dict[str,int]] = defaultdict(lambda: {
        "Nvalid_cov":0, "Nmod":0, "Ncanonical":0, "Nother_mod":0, "Ndelete":0, "Nfail":0, "Ndiff":0, "Nnocall":0, "n_sites":0
    })

    for path in beds:
        sample, mod, zn = parse_sample_mod_zn(path)
        if (sample, mod, zn) == (None, None, None):
            if args.verbose:
                print(f"[skip] cannot parse ZN from {os.path.basename(path)}", file=sys.stderr)
            continue
        if mod is None:
            # infer_mod_from_file failed (empty bed?)
            if args.verbose:
                print(f"[skip] no mod_code in {os.path.basename(path)} (empty?)", file=sys.stderr)
            continue

        # consume rows
        for row in read_bed_rows(path):
            chrom = row["chrom"]; strand = row["strand"]
            best = None; best_key = None
            for (c, s, gidx, txi), tlist in idx.items():
                if c != chrom or s != strand or txi != zn:
                    continue
                for tinfo in tlist:
                    ov = any_exon_overlap(tinfo["exons"], row["start0"], row["end0"])  # bp
                    if ov <= 0:
                        continue
                    site_mid_1b = (row["start0"] + row["end0"])//2 + 1
                    tes_dist = abs((tinfo["tes"]) - site_mid_1b)
                    score = (ov, -tes_dist)
                    if (best is None) or (score > best):
                        best = score
                        best_key = (gidx, tinfo["gene_id"], tinfo["gene_name"], zn)
            if best_key is None:
                if args.verbose:
                    print(f"[warn] site {chrom}:{row['start0']}-{row['end0']} ZN={zn} had no gene overlap; skipping", file=sys.stderr)
                continue
            gidx, gid, gname, txi = best_key

            acc_key = (gidx, txi, sample, mod)
            acc = per_gene_accum[acc_key]
            for k in ("Nvalid_cov","Nmod","Ncanonical","Nother_mod","Ndelete","Nfail","Ndiff","Nnocall"):
                acc[k] += row[k]
            acc["n_sites"] += 1

            frac = (row["Nmod"] / row["Nvalid_cov"]) if row["Nvalid_cov"] > 0 else 0.0
            if args.min_cov and row["Nvalid_cov"] < args.min_cov:
                frac = 0.0

            rows_long.append({
                "code": str(zn),
                "gene_index": gidx,
                "transcript_index": txi,
                "gtf_gene_id": gid,
                "gtf_gene_name": gname,
                "sample": sample,
                "mod_code": row["mod_code"],
                "Nvalid_cov": row["Nvalid_cov"],
                "Nmod": row["Nmod"],
                "Ncanonical": row["Ncanonical"],
                "Nother_mod": row["Nother_mod"],
                "Ndelete": row["Ndelete"],
                "Nfail": row["Nfail"],
                "Ndiff": row["Ndiff"],
                "Nnocall": row["Nnocall"],
                "n_sites": 1,
                "frac_modified": round(frac, 6),
            })

    # Write combined long table
    out_long = f"{args.out_prefix}_ZN_long.tsv"
    pd.DataFrame(rows_long).to_csv(out_long, sep="\t", index=False)
    if args.verbose:
        print(f"[ok] wrote {out_long}", file=sys.stderr)

    # Emit per-gene tables
    per_gene_dir = f"{args.out_prefix}.per_gene"
    os.makedirs(per_gene_dir, exist_ok=True)

    def safe_name(s):
        return str(s).replace("/","-").replace(" ","_")

    # Convert accumulators to DataFrame
    pg_rows = []
    for (gidx, txi, sample, mod), acc in per_gene_accum.items():
        meta = gene_meta.get(gidx, {"gtf_gene_name": "NA", "gtf_gene_id": "NA"})
        gname = meta.get("gtf_gene_name", "NA"); gid = meta.get("gtf_gene_id", "NA")
        frac = (acc["Nmod"] / acc["Nvalid_cov"]) if acc["Nvalid_cov"] > 0 else 0.0
        if args.min_cov and acc["Nvalid_cov"] < args.min_cov:
            frac = 0.0
        pg_rows.append({
            "gene_index": gidx,
            "gtf_gene_name": gname,
            "gtf_gene_id": gid,
            "transcript_index": txi,
            "sample": sample,
            "mod_code": mod,
            **acc,
            "frac_modified": round(frac,6),
        })

    if pg_rows:
        df_pg = pd.DataFrame(pg_rows)
        for gidx, df_g in df_pg.groupby("gene_index"):
            meta = gene_meta.get(gidx, {"gtf_gene_name": "NA", "gtf_gene_id": "NA"})
            gname = meta.get("gtf_gene_name", "NA"); gid = meta.get("gtf_gene_id", "NA")
            fname = os.path.join(per_gene_dir, f"{int(gidx)}_{safe_name(gname)}_{safe_name(gid)}.tsv")
            df_g.sort_values(["transcript_index","mod_code","sample"]).to_csv(fname, sep='\t', index=False)
            if args.verbose:
                print(f"[ok] wrote {fname}", file=sys.stderr)

        if args.write_gene_pivots:
            for gidx, df_g in df_pg.groupby("gene_index"):
                meta = gene_meta.get(gidx, {"gtf_gene_name": "NA", "gtf_gene_id": "NA"})
                gname = meta.get("gtf_gene_name", "NA"); gid = meta.get("gtf_gene_id", "NA")
                base = os.path.join(per_gene_dir, f"{int(gidx)}_{safe_name(gname)}_{safe_name(gid)}")
                def piv(metric, suffix):
                    piv = df_g.pivot_table(index=["transcript_index","mod_code"], columns="sample", values=metric, aggfunc="first").fillna(0).reset_index()
                    piv.to_csv(f"{base}_{suffix}.tsv", sep='\t', index=False)
                piv("frac_modified","frac_pivot")
                piv("Nvalid_cov","cov_pivot")
                piv("Nmod","Nmod_pivot")
                if args.verbose:
                    print(f"[ok] wrote {base}_*_pivot.tsv", file=sys.stderr)

    print("[OK] ZN aggregation complete.")

if __name__ == "__main__":
    main()

