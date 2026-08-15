#!/usr/bin/env python3
"""
Check the modulator pipeline outputs on the synthetic 3-exon dataset against the
known ground truth. Prints one PASS/FAIL line per feature.

Usage (from the repo root, after run_pipeline.sh):
    <modulator-env>/bin/python resources/synthetic_3exon/validate_outputs.py \
        --results results --prefix syn3exon
"""
from __future__ import annotations
import argparse
from pathlib import Path
import pandas as pd

SITES = {1080: "P", 1160: "Q", 1300: "R", 1380: "S", 2100: "X", 2260: "Y",
         3050: "D", 3200: "Msnp", 3350: "Cmod", 3450: "COND"}


def lbl(x):
    try:
        return SITES.get(int(x), str(int(x)))
    except Exception:
        return str(x)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default="results")
    ap.add_argument("--prefix", default="syn3exon")
    a = ap.parse_args()
    R = Path(a.results); P = a.prefix
    n_pass = n_fail = 0

    def check(name, cond, detail=""):
        nonlocal n_pass, n_fail
        tag = "PASS" if cond else "FAIL"
        if cond: n_pass += 1
        else: n_fail += 1
        print(f"  [{tag}] {name}" + (f"  --  {detail}" if detail else ""))

    def rd(rel):
        return pd.read_csv(R / rel, sep="\t")

    print("== ASSEMBLY / STRUCTURE ==")
    cs = rd(f"assemble/{P}_classification_summary.tsv")
    genes = cs["zt_label"].str.split(".").str[0]
    check("GENE_A/B/C/OV1/OV2/TR/TA all assembled",
          {"GENE_A", "GENE_B", "GENE_C", "GENE_OV1", "GENE_OV2", "GENE_TR", "GENE_TA"} <= set(genes))
    check("all fragmentforms classified EXACT vs reference",
          (cs["classification"] == "EXACT").all(),
          f"classes={sorted(cs['classification'].unique())}")
    check("GENE_A has exactly 2 fragmentforms (boundary-TES duplication bug fixed)",
          (genes == "GENE_A").sum() == 2, f"observed {(genes=='GENE_A').sum()}")

    print("== MULTIGENE OVERLAP FILTER ==")
    import glob as _glob
    kept = 0
    for f in _glob.glob(str(R / "assemble/zt_scrap/*.multigene_filter_summary.tsv")):
        m = pd.read_csv(f, sep="\t")
        row = m[m["metric"] == "multi_gene_kept_by_zt"]
        if len(row):
            kept += int(row.iloc[0]["value"])
    check("overlapping GENE_OV1/OV2 reads flagged multi-gene & kept by ZT (>0)",
          kept > 0, f"multi_gene_kept_by_zt total = {kept}")

    print("== TRUNCATION-AWARE STOICHIOMETRY (hierarchical) ==")
    hs = rd(f"test_diffs/{P}_hierarchical_stoich.tsv")
    tr = hs[hs.gene_name == "GENE_TR"]
    check("GENE_TR drops 5'-truncated reads as uninformative (>0)",
          len(tr) and (tr["reads_dropped_as_uninformative"] > 0).any(),
          f"dropped={int(tr['reads_dropped_as_uninformative'].max())}" if len(tr) else "no GENE_TR row")

    sj = rd(f"assemble/{P}_gene_splice_summary.tsv")
    check("all genes ALL_CANONICAL splice junctions (incl. minus-strand GENE_C)",
          (sj["intron_category"] == "ALL_CANONICAL").all())

    apa = rd(f"assemble/{P}_apa_motifs.tsv")
    prox = apa[apa.tes == 8001]
    dist = apa[apa.tes == 8481]
    check("GENE_B proximal TES -> PAS_NONE_INTERNAL_PRIMING",
          len(prox) and prox.iloc[0]["apa_motif_class"] == "PAS_NONE_INTERNAL_PRIMING")
    check("GENE_B distal TES -> PAS_CANONICAL",
          len(dist) and dist.iloc[0]["apa_motif_class"] == "PAS_CANONICAL")

    print("== BETWEEN-ISOFORM DIFFERENTIAL MODIFICATION (test_diffs) ==")
    td = rd(f"test_diffs/{P}__ZN_site_diff_results.tsv")
    sc = [c for c in td.columns if "start" in c.lower()][0]
    td["site"] = td[sc].map(lbl)
    drow = td[td.site == "D"].sort_values("p_adj_bh").head(1)
    check("site D flagged differential between isoforms (p_adj<0.05)",
          len(drow) and drow.iloc[0]["p_adj_bh"] < 0.05,
          f"p_adj={drow.iloc[0]['p_adj_bh']:.2e}" if len(drow) else "missing")
    for s in ("P", "Q", "R", "S"):
        row = td[td.site == s]
        check(f"site {s} NOT differential between isoforms (p_adj>=0.05)",
              len(row) and (row["p_adj_bh"] >= 0.05).all())

    print("== GENOTYPE: SNP + SNP-MOD + MECHANISM + HAPLOTYPE ==")
    snp = rd(f"genotype/{P}_candidate_snps.tsv")
    check("both designed SNPs discovered",
          {"chrSyn:1251:A>T", "chrSyn:3202:C>G"} <= set(snp["snp_id"]))
    sm = rd(f"genotype/{P}_snp_mod_assoc.tsv")
    ms = [c for c in sm.columns if "mod_start" in c.lower()][0]
    sm["site"] = sm[ms].map(lbl)
    snp1_msnp = sm[(sm.snp_id == "chrSyn:3202:C>G") & (sm.site == "Msnp")]
    check("SNP1 x Msnp significant (p_adj<0.05, large effect)",
          len(snp1_msnp) and snp1_msnp.iloc[0]["p_adj_bh"] < 0.05
          and snp1_msnp.iloc[0]["effect_abs_delta_mod_frac"] > 0.5,
          f"p_adj={snp1_msnp.iloc[0]['p_adj_bh']:.2e}" if len(snp1_msnp) else "missing")
    # No truly-independent site should show a STRONG SNP association. (We gate on effect
    # size, not FDR: with finite reads a null site can cross p<0.05 by chance -- exactly the
    # pipeline's own "rank by |delta|, not FDR alone" guidance -- but its effect stays small,
    # far below the designed SNP1 x Msnp effect. See GROUND_TRUTH.md.)
    indep = sm[sm.site.isin(["P", "Q", "R", "S", "D", "Cmod", "COND"])]
    max_indep = float(indep["effect_abs_delta_mod_frac"].max())
    check("no independent mod site has a strong SNP association (effect < 0.25)",
          max_indep < 0.25,
          f"max independent effect={max_indep:.2f}  vs designed SNP1xMsnp="
          f"{snp1_msnp.iloc[0]['effect_abs_delta_mod_frac']:.2f}")
    me = rd(f"genotype/{P}_snp_mod_mechanism.tsv")
    mrow = me[(me.snp_id == "chrSyn:3202:C>G") & (me.positional_class == "IN_MOTIF_CORE")]
    check("SNP1 mechanism = IN_MOTIF_CORE + MOTIF_DISRUPTED + CONCORDANT",
          len(mrow) and mrow.iloc[0]["motif_effect"] == "MOTIF_DISRUPTED"
          and mrow.iloc[0]["observed_direction"] == mrow.iloc[0]["predicted_direction"])
    hb = rd(f"genotype/{P}_haplotype_blocks.tsv")
    check("one haplotype block with 2 SNPs reconstructed",
          len(hb) == 1 and hb.iloc[0]["n_snps"] == 2)

    print("== CO-LOCALIZED MODIFICATIONS (mod x mod dependency) -- headline ==")
    mm = rd(f"genotype/{P}_mod_mod_assoc.tsv")
    mm["A"] = mm["start0_a"].map(lbl); mm["B"] = mm["start0_b"].map(lbl)

    def pair(x, y):
        r = mm[((mm.A == x) & (mm.B == y)) | ((mm.A == y) & (mm.B == x))]
        return r.iloc[0] if len(r) else None
    pq, rs, xy = pair("P", "Q"), pair("R", "S"), pair("X", "Y")
    check("P x Q -> CONCORDANT + significant (co-dependent)",
          pq is not None and pq["direction"] == "CONCORDANT" and pq["p_adj_bh"] < 0.05,
          f"OR={pq['odds_ratio']:.1f} p_adj={pq['p_adj_bh']:.1e}" if pq is not None else "missing")
    check("R x S -> INDEPENDENT + non-significant",
          rs is not None and rs["direction"] == "INDEPENDENT" and rs["p_adj_bh"] >= 0.05,
          f"OR={rs['odds_ratio']:.2f} p_adj={rs['p_adj_bh']:.2f}" if rs is not None else "missing")
    check("X x Y -> MUTUALLY_EXCLUSIVE + significant",
          xy is not None and xy["direction"] == "MUTUALLY_EXCLUSIVE" and xy["p_adj_bh"] < 0.05,
          f"OR={xy['odds_ratio']:.2f} p_adj={xy['p_adj_bh']:.1e}" if xy is not None else "missing")

    print("== POLY(A) ==")
    ptd = rd(f"polya/{P}_taillength_diffs.tsv")
    gb = ptd[ptd.gene_name == "GENE_B"]
    check("GENE_B differential tail (B1 vs B2) significant",
          len(gb) and gb.iloc[0]["p_adj_bh"] < 0.05)
    tm = rd(f"polya/{P}_taillength_mod.tsv")
    check("tail x modification detects the D coupling (>=1 site p_adj<0.05)",
          (tm["p_adj_bh"] < 0.05).any())

    print("== BETWEEN-CONDITIONS ==")
    md = rd(f"between_conditions/{P}_zikv_vs_mock_mod_diffs.tsv")
    sc2 = [c for c in md.columns if "start" in c.lower()][0]
    md["site"] = md[sc2].map(lbl)
    crow = md[md.site == "COND"]
    check("COND site differentially modified between conditions",
          len(crow) and crow.iloc[0]["p_adj_bh"] < 0.05,
          f"delta={crow.iloc[0]['delta']:.2f} p_adj={crow.iloc[0]['p_adj_bh']:.1e}" if len(crow) else "missing")
    iso = rd(f"between_conditions/{P}_zikv_vs_mock_isoform_usage_diffs.tsv")
    ga = iso[iso.gene_name == "GENE_A"]
    check("GENE_A isoform usage shifts between conditions (>=1 form p_adj<0.05)",
          (ga["p_adj_bh"] < 0.05).any())
    tl = rd(f"between_conditions/{P}_zikv_vs_mock_tail_diffs.tsv")
    check("GENE_B tail length differs between conditions",
          (tl[tl.gene_name == "GENE_B"]["p_adj_bh"] < 0.05).any())

    print("== STRUCTURAL CLASSIFICATION (classify_diffs) ==")
    cl = rd(f"test_diffs/{P}__ZN_site_classified.tsv")
    ta = cl[cl.gene_name == "GENE_TA"]
    check("GENE_TA differential site classified TANDEM_APA",
          len(ta) and (ta["structural_category"] == "TANDEM_APA").any())
    ga = cl[(cl.gene_name == "GENE_A")]
    check("GENE_A co-terminal isoform site classified SHARED_TERMINAL_EXON",
          len(ga) and (ga["structural_category"] == "SHARED_TERMINAL_EXON").any())

    print("== SEQUENCE ELEMENTS (unbiased element x modification) ==")
    se = rd(f"assemble/{P}_sequence_elements.tsv")
    check("sequence_elements produced with PAS detected (synthetic AATAAA)",
          len(se) and (se["element_type"] == "PAS").any(),
          f"element types: {sorted(se['element_type'].unique())}")
    ss = rd(f"assemble/{P}_sequence_elements_summary.tsv")
    check("sequence_elements summary has an unbiased mod-code column",
          "mod_codes_seen" in ss.columns)

    print("== REPORT ==")
    check("HTML report produced", (R / f"report/{P}_report.html").exists())
    check("interactive gene browser produced", (R / f"report/{P}_gene_browser.html").exists())

    print(f"\nSUMMARY: {n_pass} PASS / {n_fail} FAIL")
    import sys as _sys
    _sys.exit(1 if n_fail else 0)


if __name__ == "__main__":
    main()
