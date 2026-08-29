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


def main():
    assembler = load_assembler_module()
    aggregate = load_aggregate_module()
    test_support_first(assembler)
    test_intron_retention_not_truncation(assembler)
    test_metagene_coloring(assembler)
    test_aggregate_partition_mapping(aggregate)
    print("regression_smoke_checks: OK")


if __name__ == "__main__":
    main()
