#!/usr/bin/env python3

import argparse
import os

import pandas as pd

from genotype_utils import context_key_from_row, context_keys_from_snp_row, normalize_text_token


def parse_args():
    ap = argparse.ArgumentParser(description="Build unique candidate mod sites from aggregated modulator outputs.")
    ap.add_argument("--zn-long", default="", help="ZN filtered long TSV")
    ap.add_argument("--zt-long", default="", help="ZT filtered long TSV fallback")
    ap.add_argument("--candidate-snps", default="", help=(
        "Candidate SNP TSV. When given, keep only mod sites whose context_key "
        "(metagene/gene/chrom) matches a candidate SNP's context_key -- i.e. sites that can "
        "actually pair with a SNP on a shared read in snp_mod_assoc / "
        "haplotype_mod_assoc. Drops genome-wide mod sites with no linked SNP (lossless for "
        "those outputs) and keeps the per-read mod table tractable on deep data."))
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

    # Keep only SNP-linked mod sites: downstream pairing is by equal context_key on a shared
    # read, so a mod site whose context has no candidate SNP can never appear in
    # snp_mod_assoc / haplotype_mod_assoc. Dropping those is lossless
    # for those outputs and keeps the per-read mod-call table at SNP scale.
    if args.candidate_snps and os.path.exists(args.candidate_snps) and os.path.getsize(args.candidate_snps) and not out.empty:
        snps = load_input(args.candidate_snps)
        if not snps.empty:
            # A mod site loaded from the ZN long table has NO metagene_index column, so its
            # context_key is GENE:{gene_name}; the SNP side, at a single-metagene locus, yields
            # MG:{metagene}. Comparing those directly (the old code) made isin() False and dropped
            # EVERY linkable mod site -- violating this filter's "lossless" promise and silently
            # emptying snp_mod_assoc / haplotype_mod_assoc on real data.
            # Match at GENE granularity on both sides: also register each SNP's gene(s) as a
            # GENE: key. This is a lossless superset -- it never drops a mod site that could pair
            # with a SNP, and may keep a few extra same-gene/different-metagene sites (safe).
            snp_keys = set()
            for r in snps.to_dict("records"):
                # Fan a SNP out to ALL its metagene (MG:) contexts -- the singular
                # context_key_from_snp_row collapsed a multi-metagene SNP to a single CHR: key and never
                # emitted its MG: keys, dropping linkable mod sites at this UPSTREAM filter (invisible to
                # the already-fixed snp_mod_assoc / haplotype consumers). Matches 715916d.
                snp_keys.update(context_keys_from_snp_row(r))
                for g in str(r.get("gene_names", "")).split(";"):
                    g = normalize_text_token(g)
                    if g:
                        snp_keys.add(f"GENE:{g}")
            before = len(out)
            mod_keys = out.apply(context_key_from_row, axis=1)
            out = out[mod_keys.isin(snp_keys)].copy()
            print(
                f"[info] mod-site SNP-link filter: {before} -> {len(out)} sites "
                f"({len(snp_keys)} SNP context_keys)",
                flush=True,
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
