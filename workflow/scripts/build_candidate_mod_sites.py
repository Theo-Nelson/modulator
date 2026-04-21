#!/usr/bin/env python3

import argparse
import os

import pandas as pd


def parse_args():
    ap = argparse.ArgumentParser(description="Build unique candidate mod sites from aggregated modulator outputs.")
    ap.add_argument("--zn-long", default="", help="ZN filtered long TSV")
    ap.add_argument("--zt-long", default="", help="ZT filtered long TSV fallback")
    ap.add_argument("--out-tsv", required=True, help="Output candidate mod site TSV")
    ap.add_argument("--out-bed", required=True, help="Output BED for modkit extract include-bed")
    ap.add_argument("--min-total-cov", type=int, default=1, help="Minimum aggregated coverage to keep a site")
    return ap.parse_args()


def load_input(path: str) -> pd.DataFrame:
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        return pd.DataFrame()
    return pd.read_csv(path, sep="\t", low_memory=False)


def main():
    args = parse_args()
    df = load_input(args.zn_long)
    if df.empty:
        df = load_input(args.zt_long)

    cols = [
        "chrom", "start0", "end0", "strand", "mod_code",
        "gene_id", "gene_name", "metagene_index"
    ]
    if df.empty:
        out = pd.DataFrame(columns=["mod_site_id"] + cols + ["total_cov", "total_nmod", "supporting_rows", "n_samples"])
    else:
        for col in cols:
            if col not in df.columns:
                df[col] = ""
        if "Nvalid_cov" not in df.columns:
            df["Nvalid_cov"] = 0
        if "Nmod" not in df.columns:
            df["Nmod"] = 0
        if "sample" not in df.columns:
            df["sample"] = ""
        grouped = (
            df.groupby(cols, dropna=False, as_index=False)
              .agg(
                  total_cov=("Nvalid_cov", "sum"),
                  total_nmod=("Nmod", "sum"),
                  supporting_rows=("chrom", "size"),
                  n_samples=("sample", lambda x: len({str(v) for v in x if str(v)})),
              )
        )
        grouped = grouped[grouped["total_cov"] >= int(args.min_total_cov)].copy()
        grouped["mod_site_id"] = grouped.apply(
            lambda r: f"{r['chrom']}:{int(r['start0'])}-{int(r['end0'])}:{r['strand']}:{r['mod_code']}",
            axis=1,
        )
        out = grouped[["mod_site_id"] + cols + ["total_cov", "total_nmod", "supporting_rows", "n_samples"]].sort_values(
            ["chrom", "start0", "end0", "strand", "mod_code"]
        )

    os.makedirs(os.path.dirname(args.out_tsv) or ".", exist_ok=True)
    out.to_csv(args.out_tsv, sep="\t", index=False)

    with open(args.out_bed, "w") as bed:
        if not out.empty:
            for row in out.itertuples(index=False):
                bed.write(
                    "\t".join([
                        str(row.chrom),
                        str(int(row.start0)),
                        str(int(row.end0)),
                        str(row.mod_site_id),
                        "0",
                        str(row.strand),
                    ]) + "\n"
                )


if __name__ == "__main__":
    main()
