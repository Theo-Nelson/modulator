#!/usr/bin/env python3

import argparse
import os

import numpy as np
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
    ap.add_argument("--min-total-nmod", type=int, default=0,
                    help="Minimum POOLED modified-read count (sum of Nmod over samples x fragmentforms) for a "
                         "site to be a candidate. A position with a handful of modified calls across a whole "
                         "cohort can never yield an association but costs one per-read row per covering read "
                         "in the molecule table (a 31-library cohort: 1.4 M such sites, a >1 TB table). 0 = off.")
    ap.add_argument("--min-site-frac", type=float, default=0.0,
                    help="Keep a site only if at least ONE (sample x fragmentform) row with coverage >= "
                         "--min-frac-cov has modified fraction >= this. Per-row, so a modification private to "
                         "one fragmentform or one sample is kept. 0 = off.")
    ap.add_argument("--min-frac-cov", type=int, default=10, help="Row coverage floor for --min-site-frac.")
    ap.add_argument("--max-rows", type=int, default=0,
                    help="Row budget for the per-read modification table (rows ~ 2 x summed site coverage). When "
                         "the surviving sites would exceed it, sites are ranked by evidence of modification -- the "
                         "strongest (sample x fragmentform) modified fraction, then pooled modified reads -- and "
                         "kept from the top until the budget is filled. Runs that fit (a few samples) keep every "
                         "site; a large cohort is bounded automatically. 0 = no budget.")
    ap.add_argument("--max-snp-distance-bp", type=int, default=0,
                    help="Keep only sites within this distance of a candidate SNP on the same chromosome "
                         "(cis window for the SNP x modification tests). 0 = off (any SNP-linked site).")
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
        out = pd.DataFrame(columns=["mod_site_id"] + cols + ["total_cov", "total_nmod", "supporting_rows", "n_samples", "max_row_frac"])
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
        _cov = pd.to_numeric(df["Nvalid_cov"], errors="coerce").fillna(0)
        _nmod = pd.to_numeric(df["Nmod"], errors="coerce").fillna(0)
        df["_row_frac"] = np.where(_cov >= int(args.min_frac_cov), _nmod / _cov.clip(lower=1), 0.0)
        grouped = (
            df.groupby(cols, dropna=False, as_index=False)
              .agg(
                  total_cov=("Nvalid_cov", "sum"),
                  total_nmod=("Nmod", "sum"),
                  supporting_rows=("chrom", "size"),
                  n_samples=("sample", lambda x: len({str(v) for v in x if str(v)})),
                  max_row_frac=("_row_frac", "max"),
              )
        )
        n0 = len(grouped)
        grouped = grouped[grouped["total_cov"] >= int(args.min_total_cov)].copy()
        n1 = len(grouped)
        if int(args.min_total_nmod) > 0:
            grouped = grouped[grouped["total_nmod"] >= int(args.min_total_nmod)].copy()
        n2 = len(grouped)
        if float(args.min_site_frac) > 0:
            grouped = grouped[grouped["max_row_frac"] >= float(args.min_site_frac)].copy()
        print(f"[info] mod-site filters: {n0} aggregated sites -> {n1} (cov>={args.min_total_cov}) -> "
              f"{n2} (pooled Nmod>={args.min_total_nmod}) -> {len(grouped)} (a row with frac>={args.min_site_frac} "
              f"at cov>={args.min_frac_cov})", flush=True)
        grouped["max_row_frac"] = grouped["max_row_frac"].round(4)
        grouped["mod_site_id"] = grouped.apply(
            lambda r: f"{r['chrom']}:{int(r['start0'])}-{int(r['end0'])}:{r['strand']}:{r['mod_code']}",
            axis=1,
        )
        out = grouped[["mod_site_id"] + cols + ["total_cov", "total_nmod", "supporting_rows", "n_samples", "max_row_frac"]].sort_values(
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
            if int(args.max_snp_distance_bp) > 0 and not out.empty:
                # cis window: nearest candidate SNP on the same chromosome (bisect per chromosome)
                import bisect
                snp_pos = {c: np.sort(pd.to_numeric(g["pos1"], errors="coerce").dropna().astype(int).values)
                           for c, g in snps.groupby("chrom")}
                lim = int(args.max_snp_distance_bp)
                keep = []
                for c, s0 in zip(out["chrom"].astype(str), out["start0"].astype(int)):
                    a = snp_pos.get(c)
                    ok = False
                    if a is not None and a.size:
                        i = bisect.bisect_left(a, s0)
                        for j in (i - 1, i):
                            if 0 <= j < a.size and abs(int(a[j]) - s0) <= lim:
                                ok = True
                                break
                    keep.append(ok)
                before = len(out)
                out = out[np.array(keep, dtype=bool)].copy()
                print(f"[info] mod-site SNP-distance filter (<= {lim} bp): {before} -> {len(out)} sites", flush=True)
        else:
            # We were given a candidate-SNP file (require-snp-link is on) but it is header-only: there
            # are NO candidate SNPs, so NO mod site can be SNP-linked. Keep none -- the old `if not
            # snps.empty` skipped the filter entirely and silently kept EVERY mod site.
            print(f"[info] mod-site SNP-link filter: {len(out)} -> 0 sites (no candidate SNPs)", flush=True)
            out = out.iloc[0:0].copy()

    if int(args.max_rows) > 0 and not out.empty:
        proj = 2.0 * pd.to_numeric(out["total_cov"], errors="coerce").fillna(0)
        if proj.sum() > int(args.max_rows):
            ranked = out.assign(_proj=proj).sort_values(
                ["max_row_frac", "total_nmod", "total_cov"], ascending=[False, False, True], kind="stable")
            keep = ranked["_proj"].cumsum() <= int(args.max_rows)
            before = len(out)
            cut = ranked.loc[keep]
            out = cut.drop(columns=["_proj"]).sort_values(["chrom", "start0", "end0", "strand", "mod_code"]).copy()
            print(f"[info] mod-site row budget ({int(args.max_rows):,} rows): {before} -> {len(out)} sites "
                  f"(projected {int(proj.sum()):,} -> {int(cut['_proj'].sum()):,} rows; weakest kept site: "
                  f"max row fraction {float(cut['max_row_frac'].min()):.3f}, pooled modified reads "
                  f"{int(cut['total_nmod'].min())})", flush=True)
        else:
            print(f"[info] mod-site row budget ({int(args.max_rows):,} rows): all {len(out)} sites fit "
                  f"(projected {int(proj.sum()):,} rows)", flush=True)
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
