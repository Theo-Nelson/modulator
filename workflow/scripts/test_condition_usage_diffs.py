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
import json
import os
import re
import sys

import numpy as np
import pandas as pd

import diffstats
from genotype_utils import benjamini_hochberg

OUT_COLS = ["contrast", "feature_type", "feature", "gene_name",
            "n_reference", "n_test", "reads_reference", "reads_test",
            "mu_reference", "mu_test", "delta", "dispersion", "lrt_stat", "p_value", "p_adj_bh",
            "per_replicate_json"]


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


def _fallback_gene_of(zt: str) -> str:
    """Last-resort gene key when the authoritative gtf_gene_name is unavailable (e.g. a NOVEL_LOCUS
    with gtf_gene_name = NA). The zt_label is `{gene}.{gene_id}.G<n>.T<n>` and BOTH gene and gene_id
    can contain dots (GENCODE clone names like CTC-338M12.4, versioned Ensembl ids), so split(".")[0]
    (M4) truncated dotted names and MERGED distinct genes' denominators (CTC-338M12.4 and .3 -> one
    "CTC-338M12"). Strip only the `.G<n>.T<n>` suffix, keeping `{gene}.{gene_id}` -- unique per gene,
    so nothing is merged. Used only as a fallback; a real gtf_gene_name always wins."""
    return re.sub(r"\.G\d+\.T\d+$", "", str(zt))


def _load_gene_map(args) -> dict:
    """zt_label -> AUTHORITATIVE gtf_gene_name from the classification summary. Blank/NA are dropped so
    the string fallback handles those. Empty dict if no summary was passed."""
    if not getattr(args, "classification_summary", ""):
        return {}
    s = pd.read_csv(args.classification_summary, sep="\t", low_memory=False)
    s.columns = [str(c).lstrip("#") for c in s.columns]
    gcol = "gtf_gene_name" if "gtf_gene_name" in s.columns else ("gene_name" if "gene_name" in s.columns else None)
    if "zt_label" not in s.columns or gcol is None:
        return {}
    g = s[["zt_label", gcol]].dropna(subset=["zt_label"]).copy()
    g[gcol] = g[gcol].astype(str).str.strip()
    g = g[~g[gcol].str.lower().isin(("", "na", "nan", "none"))]
    return dict(zip(g["zt_label"].astype(str), g[gcol]))


def _feature_map(args, zts, gene_of) -> pd.DataFrame:
    """fragmentform -> (feature, gene). One row per (zt_label, feature); junctions fan out.
    `gene_of(zt)` returns the authoritative gtf_gene_name (or the no-merge fallback)."""
    if args.feature == "isoform":
        return pd.DataFrame({"zt_label": zts, "feature": zts, "gene_name": [gene_of(z) for z in zts]})
    if args.feature == "apa":
        if not args.classification_summary:
            sys.exit("--feature apa needs --classification-summary")
        s = pd.read_csv(args.classification_summary, sep="\t", low_memory=False)
        s.columns = [str(c).lstrip("#") for c in s.columns]
        s = s[["zt_label", "iso_tes", "chrom", "strand"]].dropna(subset=["zt_label", "iso_tes"])
        s["feature"] = (s["chrom"].astype(str) + ":" + s["iso_tes"].astype(int).astype(str)
                        + ":" + s["strand"].astype(str))
        s["gene_name"] = s["zt_label"].map(gene_of)
        return s[["zt_label", "feature", "gene_name"]]
    if not args.splice_junctions:
        sys.exit("--feature junction needs --splice-junctions")
    j = pd.read_csv(args.splice_junctions, sep="\t", low_memory=False)
    j.columns = [str(c).lstrip("#") for c in j.columns]
    j = j.dropna(subset=["zt_label", "intron_start1", "intron_end1"])
    j["feature"] = (j["chrom"].astype(str) + ":" + j["intron_start1"].astype(int).astype(str)
                    + "-" + j["intron_end1"].astype(int).astype(str) + ":" + j["strand"].astype(str))
    j["gene_name"] = j["zt_label"].map(gene_of)
    return j[["zt_label", "feature", "gene_name"]].drop_duplicates()


def main():
    args = parse_args()
    # Ensure the output directory exists BEFORE any early-return writes an empty table -- the
    # pipeline does not pre-create between_conditions/, so a no-result contrast otherwise crashes.
    os.makedirs(os.path.dirname(args.out_tsv) or ".", exist_ok=True)
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

    # Authoritative zt_label -> gene map (gtf_gene_name); fall back to a NO-MERGE key for anything
    # absent (novel loci). Used for BOTH the feature->gene assignment and the gene denominators, so a
    # dotted gene name can no longer be truncated or two distinct genes merged (M4).
    _gene_map = _load_gene_map(args)
    def gene_of(zt):
        return _gene_map.get(str(zt)) or _fallback_gene_of(zt)

    fmap = _feature_map(args, list(tx.index), gene_of)
    fmap = fmap[fmap["zt_label"].isin(tx.index)]
    if fmap.empty:
        pd.DataFrame(columns=OUT_COLS).to_csv(args.out_tsv, sep="\t", index=False)
        return

    # feature counts = sum of its fragmentforms; gene totals = sum of ALL the gene's fragmentforms.
    counts = tx.loc[fmap["zt_label"]].to_numpy(dtype=float)
    fc = pd.DataFrame(counts, columns=samples)
    fc["feature"] = fmap["feature"].to_numpy(); fc["gene_name"] = fmap["gene_name"].to_numpy()
    feat = fc.groupby(["gene_name", "feature"], sort=False)[samples].sum()

    gene_of_zt = pd.Series([gene_of(z) for z in tx.index], index=tx.index)
    gene_tot = tx.groupby(gene_of_zt, sort=False)[samples].sum()

    genes = feat.index.get_level_values(0)
    tot = gene_tot.reindex(genes)[samples].to_numpy(dtype=float)
    K = feat[samples].to_numpy(dtype=float)
    K = np.minimum(K, tot)   # defensive: a feature count can never exceed its gene total
    gidx = np.array([0 if s in ref_s else 1 for s in samples], dtype=int)

    # NOTE: requiring EVERY sample >= min_gene_reads (not >= min_samples_per_group covered per group)
    # silently drops features covered in most-but-not-all replicates; documented in the stress report.
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
        i = r["key"]
        gene, feature = idx[i]
        with np.errstate(divide="ignore", invalid="ignore"):
            frac_i = np.where(tot[i] > 0, K[i] / tot[i], np.nan)
        per_rep = {"reference": {s: round(float(frac_i[j]), 4)
                                 for j, s in enumerate(samples) if s in ref_s and np.isfinite(frac_i[j])},
                   "test": {s: round(float(frac_i[j]), 4)
                            for j, s in enumerate(samples) if s in test_s and np.isfinite(frac_i[j])}}
        rows.append({
            "contrast": name, "feature_type": args.feature, "feature": feature, "gene_name": gene,
            "n_reference": r["n_reference"], "n_test": r["n_test"],
            "reads_reference": r["reads_reference"], "reads_test": r["reads_test"],
            "mu_reference": round(r["mu_reference"], 5), "mu_test": round(r["mu_test"], 5),
            "delta": round(r["delta"], 5), "dispersion": round(r["dispersion"], 6),
            "lrt_stat": round(r["lrt_stat"], 4), "p_value": r["p_value"],
            "per_replicate_json": json.dumps(per_rep, separators=(",", ":")),
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
