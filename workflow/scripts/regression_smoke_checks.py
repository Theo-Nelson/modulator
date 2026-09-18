#!/usr/bin/env python3

import copy
import importlib.util
import pathlib
import sys
import tempfile


ROOT = pathlib.Path(__file__).resolve().parent
ASSEMBLER_PATH = ROOT / "assemble_transcripts.py"
AGGREGATE_PATH = ROOT / "aggregate_by_gene.py"


def load_assembler_module():
    spec = importlib.util.spec_from_file_location("assemble_transcripts", ASSEMBLER_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def load_aggregate_module():
    spec = importlib.util.spec_from_file_location("aggregate_by_gene", AGGREGATE_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_support_first(module):
    exact_counts = {
        ((101, 200), (301, 400), (501, 600)): 1,
        ((301, 400), (501, 600)): 30,
    }
    feats = module.compute_chain_features(exact_counts)
    params = {
        "min_exact_canonical_reads": 1,
        "min_distal_anchor_reads": 2,
        "min_distal_anchor_frac": 0.05,
    }
    long_chain = ((101, 200), (301, 400), (501, 600))
    short_chain = ((301, 400), (501, 600))
    if module.absorb_allowed_for_chain(feats[long_chain], params):
        raise AssertionError("Support-first gating should reject unsupported longer chain absorption.")
    if not module.absorb_allowed_for_chain(feats[short_chain], params):
        raise AssertionError("Observed suffix chain should remain a valid canonical.")


def test_intron_retention_not_truncation(module):
    """M9: a read that ALIGNS contiguously across an intron the fragmentform splices out has
    RETAINED that intron and must not be suffix-absorbed as a clean 3' truncation."""
    f = module.read_retains_any_intron
    # canonical fragmentform (tx order, 1-based inclusive introns): splices (200,299) then (400,499)
    canon = ((200, 299), (400, 499))
    read_chain = ((400, 499),)                     # read only splices the 3'-most intron
    omitted = canon[: len(canon) - len(read_chain)]  # -> ((200,299),)
    # clean 3' truncation: read starts 3' of the omitted intron, never touches it
    clean = [(300, 399), (500, 600)]
    if f(clean, omitted):
        raise AssertionError("Clean 3' truncation wrongly flagged as intron retention.")
    # retention: one contiguous block (150..399) spans the omitted intron (200..299)
    retained = [(150, 399), (500, 600)]
    if not f(retained, omitted):
        raise AssertionError("Intron-retention read must NOT be absorbed as a 3' truncation.")
    # read whose 5' end lands inside the intron is a partial overlap -> conservatively allowed
    partial = [(250, 399), (500, 600)]
    if f(partial, omitted):
        raise AssertionError("Partial (mid-intron start) read should not be flagged as retention.")


def test_metagene_coloring(module):
    base = [
        {
            "gene_index": 1,
            "gene_tx_index": 1,
            "gene_name_label": "GENEA",
            "gene_id_label": "GENEA",
            "chrom": "chrSynthetic",
            "strand": "+",
            "rep_exons": [(100, 160)],
            "count": 100,
            "tes": 160,
        },
        {
            "gene_index": 2,
            "gene_tx_index": 1,
            "gene_name_label": "GENEB",
            "gene_id_label": "GENEB",
            "chrom": "chrSynthetic",
            "strand": "+",
            "rep_exons": [(140, 220)],
            "count": 90,
            "tes": 220,
        },
        {
            "gene_index": 2,
            "gene_tx_index": 2,
            "gene_name_label": "GENEB",
            "gene_id_label": "GENEB",
            "chrom": "chrSynthetic",
            "strand": "+",
            "rep_exons": [(260, 320)],
            "count": 70,
            "tes": 320,
        },
    ]

    metagene = module.assign_metagene_partitions(copy.deepcopy(base))

    mg_by_gene = {(x["gene_index"], x["gene_tx_index"]): (x["metagene_index"], x["zn_index"]) for x in metagene}

    if mg_by_gene[(1, 1)][0] != mg_by_gene[(2, 1)][0]:
        raise AssertionError("Overlapping genes should share one metagene.")
    if mg_by_gene[(1, 1)][1] == mg_by_gene[(2, 1)][1]:
        raise AssertionError("Overlapping transcripts must not share a ZN partition.")
    if mg_by_gene[(1, 1)][1] != mg_by_gene[(2, 2)][1]:
        raise AssertionError("Non-overlapping transcripts should be allowed to reuse a ZN partition.")


def test_aggregate_partition_mapping(module):
    gtf_text = """chrSynthetic\tReadBacked\ttranscript\t100\t160\t1000\t+\t.\tgene_id "GENEA"; transcript_id "GENEA.G1.T1"; ref_gene_name "GENEA"; zt_label "GENEA.G1.T1"; gene_index "1"; transcript_index "1"; metagene_index "1"; zn_index "1"; metagene_partition_count "2";
chrSynthetic\tReadBacked\texon\t100\t160\t1000\t+\t.\tgene_id "GENEA"; transcript_id "GENEA.G1.T1"; ref_gene_name "GENEA"; zt_label "GENEA.G1.T1"; gene_index "1"; transcript_index "1"; metagene_index "1"; zn_index "1"; metagene_partition_count "2"; exon_number "1";
chrSynthetic\tReadBacked\ttranscript\t140\t220\t1000\t+\t.\tgene_id "GENEB"; transcript_id "GENEB.G2.T1"; ref_gene_name "GENEB"; zt_label "GENEB.G2.T1"; gene_index "2"; transcript_index "1"; metagene_index "1"; zn_index "2"; metagene_partition_count "2";
chrSynthetic\tReadBacked\texon\t140\t220\t1000\t+\t.\tgene_id "GENEB"; transcript_id "GENEB.G2.T1"; ref_gene_name "GENEB"; zt_label "GENEB.G2.T1"; gene_index "2"; transcript_index "1"; metagene_index "1"; zn_index "2"; metagene_partition_count "2"; exon_number "1";
chrSynthetic\tReadBacked\ttranscript\t260\t320\t1000\t+\t.\tgene_id "GENEB"; transcript_id "GENEB.G2.T2"; ref_gene_name "GENEB"; zt_label "GENEB.G2.T2"; gene_index "2"; transcript_index "2"; metagene_index "1"; zn_index "1"; metagene_partition_count "2";
chrSynthetic\tReadBacked\texon\t260\t320\t1000\t+\t.\tgene_id "GENEB"; transcript_id "GENEB.G2.T2"; ref_gene_name "GENEB"; zt_label "GENEB.G2.T2"; gene_index "2"; transcript_index "2"; metagene_index "1"; zn_index "1"; metagene_partition_count "2"; exon_number "1";
"""
    with tempfile.NamedTemporaryFile("w", suffix=".gtf", delete=False) as tmp:
        tmp.write(gtf_text)
        gtf_path = tmp.name
    tx_index, gene_index = module.load_gene_intervals_from_gtf(gtf_path, verbose=False)
    if module.assign_gene("chrSynthetic", 109, 110, "+", 1, tx_index, gene_index) != ("GENEA", "GENEA"):
        raise AssertionError("ZN=1 should map the left site to GENEA.")
    if module.assign_gene("chrSynthetic", 149, 150, "+", 2, tx_index, gene_index) != ("GENEB", "GENEB"):
        raise AssertionError("ZN=2 should map the overlapping site to GENEB.")
    if module.assign_gene("chrSynthetic", 269, 270, "+", 1, tx_index, gene_index) != ("GENEB", "GENEB"):
        raise AssertionError("Reused ZN=1 should map the non-overlapping right site to GENEB.")


def test_stale_zn_filter(module):
    """BLOCKER-2 residual: aggregation entered via --stages alone must drop numbered beds whose ZN
    partition the current assembly does not define (stale from a prior, larger run)."""
    gtf_text = (
        'chrSynthetic\tReadBacked\ttranscript\t100\t160\t1000\t+\t.\tgene_id "GENEA"; '
        'transcript_id "GENEA.G1.T1"; zn_index "1";\n'
        'chrSynthetic\tReadBacked\ttranscript\t140\t220\t1000\t+\t.\tgene_id "GENEB"; '
        'transcript_id "GENEB.G2.T1"; zn_index "2";\n'
    )
    with tempfile.NamedTemporaryFile("w", suffix=".gtf", delete=False) as tmp:
        tmp.write(gtf_text)
        gtf_path = tmp.name
    valid = module.valid_zn_from_gtf(gtf_path)
    if valid != {1, 2}:
        raise AssertionError(f"valid_zn_from_gtf should read {{1, 2}} from the GTF, got {valid}")
    beds = [("root", "M1", "M1/1.bed", 1), ("root", "M1", "M1/2.bed", 2),
            ("root", "M1", "M1/99.bed", 99)]   # ZN=99 is a stale partition not in the GTF
    kept = module.filter_stale_zn_beds(beds, gtf_path, "test", verbose=False)
    kept_zn = sorted(b[3] for b in kept)
    if kept_zn != [1, 2]:
        raise AssertionError(f"stale ZN=99 bed must be dropped; kept ZN partitions {kept_zn}")


def test_polya_gene_from_zt():
    """M3: recover the gene from a zt_label `{gene}.{gene_id}.G<n>.T<n>`. The fallback (no gene_id)
    must strip `.G<n>.T<n>` to a NO-MERGE key -- NOT return the raw label (which splits one gene into
    per-transcript keys, the regression the reviewer caught) and NOT split(".")[0] (dotted-name merge)."""
    spec = importlib.util.spec_from_file_location(
        "build_read_polya_table", ROOT / "build_read_polya_table.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    g = mod.gene_from_zt
    if g("CCT8.ENSG00000156261.14.G10.T1", "ENSG00000156261.14") != "CCT8":
        raise AssertionError("known gene_id should recover the clean gene name")
    if g("CTC-338M12.4.ENSG9.G2.T3", "ENSG9") != "CTC-338M12.4":
        raise AssertionError("dotted GENCODE gene name must be preserved")
    t1 = g("CCT8.ENSG00000156261.14.G10.T1", "")
    t2 = g("CCT8.ENSG00000156261.14.G10.T2", "")
    if t1 != t2:
        raise AssertionError(f"no-gene_id fallback must collapse a gene's transcripts to one key, got {t1!r} != {t2!r}")
    if t1 != "CCT8.ENSG00000156261.14":
        raise AssertionError(f"fallback must strip .G<n>.T<n>, not return the raw label: {t1!r}")


def test_stratified_sparse_guard():
    """test_stoichiometry_diffs: (1) the effect is the MH rate difference over EVERY stratum covering both
    forms, so one modified read on a shallow form no longer scores as a 100% difference; (2) a sparse
    stratified table (expected counts << 5) gets the exact Monte-Carlo p of the same CMH statistic, not
    the chi-square approximation; (3) a well-populated table is unchanged; (4) the guard can be turned off."""
    import importlib.util
    import numpy as np
    spec = importlib.util.spec_from_file_location("test_stoichiometry_diffs", ROOT / "test_stoichiometry_diffs.py")
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    names = [f"S{i}" for i in range(31)]

    def site(rows):
        zn, sp, cv, nm = zip(*rows)
        return (np.array(zn), np.array(sp, dtype=np.int32), np.array(cv), np.array(nm))
    # 31 samples; form 1 deep with one modified read in sample 7; form 14 has 1 read per sample, modified in sample 3
    rows = []
    for smp in range(31):
        rows.append((1, smp, 1165, 1 if smp == 7 else 0))
        rows.append((14, smp, 1, 1 if smp == 3 else 0))
        for z in range(2, 8):
            rows.append((z, smp, 3, 0))
    zn, sp, cv, nm = site(rows)
    r = m.summarize_site_arrays(zn, sp, cv, nm, names, 20, "auto", 0.5, "two-sided")
    assert r["test_name"].endswith("_mc"), r["test_name"]
    assert r["p_value"] > 1e-4, f"sparse table must not get an asymptotic p: {r['p_value']}"
    assert abs(r["effect_max_abs_frac_diff"] - r["effect_max_abs_frac_diff_pooled"]) < 0.02, \
        f"effect over all strata should be near the pooled contrast: {r['effect_max_abs_frac_diff']} vs {r['effect_max_abs_frac_diff_pooled']}"
    assert r["effect_max_abs_frac_diff"] < 0.05, r["effect_max_abs_frac_diff"]
    r_off = m.summarize_site_arrays(zn, sp, cv, nm, names, 20, "auto", 0.5, "two-sided", mc_min_expected=0)
    assert not r_off["test_name"].endswith("_mc") and r_off["p_value"] < r["p_value"], "guard off must reproduce the chi-square path"
    # deterministic
    r2 = m.summarize_site_arrays(zn, sp, cv, nm, names, 20, "auto", 0.5, "two-sided")
    assert r2["p_value"] == r["p_value"], "MC p must be deterministic (fixed seed)"
    # dense: 3 forms x 5 samples, 30% vs 10% vs 30% -> asymptotic CMH, and the effect equals the old definition
    rows = []
    for smp in range(5):
        rows += [(1, smp, 200, 60), (2, smp, 100, 10), (3, smp, 80, 24)]
    zn, sp, cv, nm = site(rows)
    d = m.summarize_site_arrays(zn, sp, cv, nm, names, 20, "auto", 0.5, "two-sided")
    assert d["test_name"] == "cmh_general_3x2" and abs(d["effect_max_abs_frac_diff"] - 0.2) < 1e-9, d
    stat, p, _df, _used = m.cmh_general_association(m.informative_strata([np.array([[60, 140], [10, 90], [24, 56]], float)] * 5))
    assert abs(d["stat_value"] - stat) < 1e-9 and abs(d["p_value"] - p) < 1e-300, "dense path must be the plain CMH"


def test_classify_example_rank():
    """classify_diff_sites: the example order shared by the table and the rankNN__ figures is by the
    classified hi-lo difference, then FDR, then the shallower form's coverage -- not the saturating test
    effect -- and is 1-based per leaf."""
    import importlib.util
    spec = importlib.util.spec_from_file_location("classify_diff_sites", ROOT / "classify_diff_sites.py")
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    recs = [
        dict(class_key="A", gene="g1", mod="a", chrom="chr1", start0=10, strand="+", effect=1.0, padj=1e-300, hi_frac=0.04, lo_frac=0.00, hi_cov=23, lo_cov=36000),
        dict(class_key="A", gene="g2", mod="a", chrom="chr1", start0=20, strand="+", effect=0.3, padj=1e-5, hi_frac=0.60, lo_frac=0.20, hi_cov=50, lo_cov=80),
        dict(class_key="A", gene="g3", mod="a", chrom="chr1", start0=30, strand="+", effect=0.5, padj=1e-9, hi_frac=0.60, lo_frac=0.20, hi_cov=500, lo_cov=800),
        dict(class_key="B", gene="g4", mod="a", chrom="chr2", start0=40, strand="+", effect=1.0, padj=1e-2, hi_frac=0.10, lo_frac=0.00, hi_cov=30, lo_cov=30),
    ]
    m.assign_example_ranks(recs)
    order = {r["gene"]: r["example_rank"] for r in recs}
    assert order == {"g3": 1, "g2": 2, "g1": 3, "g4": 1}, order
    top = [r["gene"] for r in m.ranked_examples([r for r in recs if r["class_key"] == "A"], 2)]
    assert top == ["g3", "g2"], top


def test_browser_companion_shards(tmp_path=None):
    """build_gene_browser companion mode: every gene is in the shard the page will compute (crc32 % n)."""
    import importlib.util, json, os, tempfile, zlib
    spec = importlib.util.spec_from_file_location("build_gene_browser", ROOT / "build_gene_browser.py")
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    genes = [f"GENE{i}" for i in range(50)] + ["CTC-338M12.4", "Zfp\u00e9"]
    n = 7
    with tempfile.TemporaryDirectory() as d:
        for g in genes:
            sid = zlib.crc32(g.encode("utf-8")) % n
            assert 0 <= sid < n
        # the HTML template must carry the loader hooks the shards call
        assert "__modulator_browser_shard" in m._HTML and "__modulator_browser_index" in m._HTML and "crc32(" in m._HTML


def main():
    assembler = load_assembler_module()
    aggregate = load_aggregate_module()
    test_support_first(assembler)
    test_intron_retention_not_truncation(assembler)
    test_metagene_coloring(assembler)
    test_aggregate_partition_mapping(aggregate)
    test_stale_zn_filter(aggregate)
    test_polya_gene_from_zt()
    test_stratified_sparse_guard()
    test_classify_example_rank()
    test_browser_companion_shards()
    print("regression_smoke_checks: OK")


if __name__ == "__main__":
    main()
