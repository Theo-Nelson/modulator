"""Differential USAGE between conditions: isoform, APA site, or splice junction.

All three are the same question -- "what fraction of this gene's reads use feature X, and does that
fraction change between conditions?" -- so they share one engine. The only difference is how a
fragmentform maps to a feature:

  --feature isoform   feature = the fragmentform itself      (which isoform does the gene use?)
  --feature apa       feature = the fragmentform's TES       (which 3' end / APA site does it use?)
  --feature junction  feature = each intron it contains      (which junctions does it splice?)
                      (a fragmentform carries many junctions, so it contributes to each)

Counts come from the per-sample tx_counts matrix; the denominator is the gene's total reads in that
sample. That makes each test a (successes, total) per sample comparison -> the same replicate-aware
beta-binomial LRT with dispersion shrinkage used for differential modification (see diffstats.py).
Never pools reads across replicates (that would be pseudoreplication; measured 62% FPR).
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

import diffstats
from genotype_utils import benjamini_hochberg

OUT_COLS = ["contrast", "feature_type", "feature", "gene_name",
            "n_reference", "n_test", "reads_reference", "reads_test",
            "mu_reference", "mu_test", "delta", "dispersion", "lrt_stat", "p_value", "p_adj_bh"]


def parse_args():
    ap = argparse.ArgumentParser(description="Replicate-aware differential isoform/APA/junction usage between conditions.")
    ap.add_argument("--tx-counts", required=True, help="*_tx_counts.tsv (fragmentform x sample counts)")
    ap.add_argument("--sample-metadata", required=True)
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--feature", choices=["isoform", "apa", "junction"], required=True)
    ap.add_argument("--classification-summary", default="", help="required for --feature apa (zt_label -> iso_tes)")
    ap.add_argument("--splice-junctions", default="", help="required for --feature junction (zt_label -> introns)")
    ap.add_argument("--contrast-name", default="")
    ap.add_argument("--column", default="condition")
    ap.add_argument("--test", required=True)
    ap.add_argument("--reference", required=True)
    ap.add_argument("--min-gene-reads", type=int, default=20, help="min gene reads per sample to test its features")
    ap.add_argument("--min-samples-per-group", type=int, default=2)
    ap.add_argument("--prior-weight", type=float, default=20.0)
    ap.add_argument("--ref-df", type=int, default=diffstats.REF_DF)
    ap.add_argument("--site-weight", default="auto", help="dispersion-shrinkage per-site weight; 'auto'=N_site-2 (scales with cohort)")
    ap.add_argument("--verbose", action="store_true")
    return ap.parse_args()


def _gene_of(zt: str) -> str:
    return str(zt).split(".")[0]


def _feature_map(args, zts) -> pd.DataFrame:
    """fragmentform -> (feature, gene). One row per (zt_label, feature); junctions fan out."""
    if args.feature == "isoform":
        return pd.DataFrame({"zt_label": zts, "feature": zts, "gene_name": [_gene_of(z) for z in zts]})
    if args.feature == "apa":
        if not args.classification_summary:
            sys.exit("--feature apa needs --classification-summary")
        s = pd.read_csv(args.classification_summary, sep="\t", low_memory=False)
        s.columns = [str(c).lstrip("#") for c in s.columns]
        gene_col = "gtf_gene_name" if "gtf_gene_name" in s.columns else "gene_name"
        s = s[["zt_label", "iso_tes", "chrom", "strand", gene_col]].dropna(subset=["zt_label", "iso_tes"])
        s["feature"] = (s["chrom"].astype(str) + ":" + s["iso_tes"].astype(int).astype(str)
                        + ":" + s["strand"].astype(str))
        s["gene_name"] = s[gene_col].fillna("").astype(str).replace("", np.nan)
        s["gene_name"] = s["gene_name"].fillna(s["zt_label"].map(_gene_of))
        return s[["zt_label", "feature", "gene_name"]]
    if not args.splice_junctions:
        sys.exit("--feature junction needs --splice-junctions")
    j = pd.read_csv(args.splice_junctions, sep="\t", low_memory=False)
    j.columns = [str(c).lstrip("#") for c in j.columns]
    j = j.dropna(subset=["zt_label", "intron_start1", "intron_end1"])
    j["feature"] = (j["chrom"].astype(str) + ":" + j["intron_start1"].astype(int).astype(str)
                    + "-" + j["intron_end1"].astype(int).astype(str) + ":" + j["strand"].astype(str))
    j["gene_name"] = j.get("gene_name", pd.Series(index=j.index, dtype=str)).fillna(
        j["zt_label"].map(_gene_of))
    return j[["zt_label", "feature", "gene_name"]].drop_duplicates()


def main():
    args = parse_args()
    name = args.contrast_name or f"{args.test}_vs_{args.reference}"

    meta = pd.read_csv(args.sample_metadata, sep="\t", low_memory=False, keep_default_na=False)
    grp = dict(zip(meta["sample"].astype(str), meta[args.column].astype(str)))
    ref_s = [s for s, g in grp.items() if g == args.reference]
    test_s = [s for s, g in grp.items() if g == args.test]

    tx = pd.read_csv(args.tx_counts, sep="\t", index_col=0, low_memory=False)
    tx.index = tx.index.astype(str)
    samples = [s for s in tx.columns if s in ref_s or s in test_s]
    ref_in = [s for s in samples if s in ref_s]; test_in = [s for s in samples if s in test_s]
    if len(ref_in) < args.min_samples_per_group or len(test_in) < args.min_samples_per_group:
        print(f"[condition_usage] {name}: need >={args.min_samples_per_group}/group "
              f"(got {len(ref_in)} vs {len(test_in)})", file=sys.stderr, flush=True)
        pd.DataFrame(columns=OUT_COLS).to_csv(args.out_tsv, sep="\t", index=False)
        return
    tx = tx[samples].fillna(0)

    fmap = _feature_map(args, list(tx.index))
    fmap = fmap[fmap["zt_label"].isin(tx.index)]
    if fmap.empty:
        pd.DataFrame(columns=OUT_COLS).to_csv(args.out_tsv, sep="\t", index=False)
        return

    # feature counts = sum of its fragmentforms; gene totals = sum of ALL the gene's fragmentforms.
    counts = tx.loc[fmap["zt_label"]].to_numpy(dtype=float)
    fc = pd.DataFrame(counts, columns=samples)
    fc["feature"] = fmap["feature"].to_numpy(); fc["gene_name"] = fmap["gene_name"].to_numpy()
    feat = fc.groupby(["gene_name", "feature"], sort=False)[samples].sum()

    gene_of_zt = pd.Series([_gene_of(z) for z in tx.index], index=tx.index)
    gene_tot = tx.groupby(gene_of_zt, sort=False)[samples].sum()

    genes = feat.index.get_level_values(0)
    tot = gene_tot.reindex(genes)[samples].to_numpy(dtype=float)
    K = feat[samples].to_numpy(dtype=float)
    gidx = np.array([0 if s in ref_s else 1 for s in samples], dtype=int)

    keep = (tot.min(axis=1) >= args.min_gene_reads)
    # a feature that is the gene's ONLY one is 100% by construction -> nothing to compare. Use an EXACT
    # integer compare (K == tot): np.isclose's default rtol=1e-5 would wrongly drop a feature that is
    # 99.999% (off by 1-2 reads) at very high coverage.
    varies = ~(K == tot).all(axis=1)
    keep &= varies
    K, tot = K[keep], tot[keep]
    idx = feat.index[keep]
    if args.verbose:
        print(f"[condition_usage:{args.feature}] {name}: {len(idx):,} features testable "
              f"({len(ref_in)} {args.reference} vs {len(test_in)} {args.test})", flush=True)
    if not len(idx):
        pd.DataFrame(columns=OUT_COLS).to_csv(args.out_tsv, sep="\t", index=False)
        return

    sites = [(i, K[i], tot[i], gidx) for i in range(K.shape[0])]
    res = diffstats.beta_binomial_diff(sites, prior_weight=args.prior_weight,
                                       min_group_samples=args.min_samples_per_group,
                                       ref_df=args.ref_df, calibrate=False,
                                       site_weight=diffstats.parse_site_weight(args.site_weight))
    rows = []
    for r in res:
        gene, feature = idx[r["key"]]
        rows.append({
            "contrast": name, "feature_type": args.feature, "feature": feature, "gene_name": gene,
            "n_reference": r["n_reference"], "n_test": r["n_test"],
            "reads_reference": r["reads_reference"], "reads_test": r["reads_test"],
            "mu_reference": round(r["mu_reference"], 5), "mu_test": round(r["mu_test"], 5),
            "delta": round(r["delta"], 5), "dispersion": round(r["dispersion"], 6),
            "lrt_stat": round(r["lrt_stat"], 4), "p_value": r["p_value"],
        })
    out = pd.DataFrame(rows)
    if out.empty:
        pd.DataFrame(columns=OUT_COLS).to_csv(args.out_tsv, sep="\t", index=False)
        return
    out["p_adj_bh"] = benjamini_hochberg(out["p_value"].values)
    out["_abs"] = out["delta"].abs()
    out = out.sort_values(["p_adj_bh", "_abs"], ascending=[True, False]).drop(columns="_abs")
    out = out[OUT_COLS].reset_index(drop=True)
    os.makedirs(os.path.dirname(args.out_tsv) or ".", exist_ok=True)
    out.to_csv(args.out_tsv, sep="\t", index=False)
    if args.verbose:
        print(f"[condition_usage:{args.feature}] {name}: {len(out):,} tested, "
              f"{int((out['p_adj_bh'] < 0.05).sum()):,} at FDR<0.05 -> {args.out_tsv}", flush=True)


if __name__ == "__main__":
    main()
