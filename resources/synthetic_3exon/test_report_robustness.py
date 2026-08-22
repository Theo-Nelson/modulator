#!/usr/bin/env python3
"""Regression tests for report-generator robustness edge cases found in stress-test campaign #2.

A malformed or missing field in ONE input row must never abort the whole HTML report, and modification
codes must stay labelled even when a TSV coerces them to float. Run:
  <modulator-env>/bin/python resources/synthetic_3exon/test_report_robustness.py
"""
import os
import sys
import tempfile

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "workflow", "scripts"))
import generate_html_report as R  # noqa: E402


def _no_crash(fn, *a, **k):
    try:
        fn(*a, **k)
        return True
    except Exception as e:  # noqa: BLE001
        print(f"    !! raised {type(e).__name__}: {e}")
        return False


def checks():
    out = []

    # A-F1: numeric mod codes coerced to float64 ('17802.0') still label + drive the glossary
    out.append(("label_mod_code(17802.0) -> pseudouridine",
                R.label_mod_code(17802.0) == "17802 (pseudouridine)"))
    out.append(("df_to_html relabels float mod_code",
                "pseudouridine" in R.df_to_html(pd.DataFrame({"mod_code": [17802.0, 17596.0, np.nan]}))))
    g = R.build_glossary_section({"17802.0", "17596.0"})
    out.append(("glossary(float codes) names pseU, no 'unknown code'",
                ("pseudouridine" in g) and ("unknown code" not in g)))

    # B-F1: malformed per_replicate_json (list / scalar / non-numeric value) must not raise
    bad = pd.DataFrame({
        "gene_name": ["G1", "G2", np.nan], "ZN_transcript_index": [0, 1, 2],
        "chrom": ["c", "c", "c"], "start0": [100, 200, 300], "mod_code": ["a", "a", "a"],
        "delta": [0.5, 0.4, 0.3],
        "per_replicate_json": ['{"reference":[1,2,3],"test":{"r":0.5}}',
                               '{"reference":{"a":0.1},"test":{"b":"abc"}}',
                               'not json at all'],
    })
    out.append(("between_cond_topn_png tolerates malformed per_replicate_json",
                _no_crash(R.between_cond_topn_png, bad, "delta", "t", "x", "mock", "zikv", "c")))

    # B-F5: NaN gene name must not render the literal 'nan'
    lab = R._bc_feature_label(pd.Series({"gene_name": np.nan, "chrom": "c", "start0": 100,
                                         "ZN_transcript_index": 1, "mod_code": "a"}))
    out.append(("_bc_feature_label(NaN gene) has no 'nan'", "nan" not in lab.lower()))

    # B-F7: malformed mods_json (list of non-dicts / scalar) must not raise
    out.append(("_elem_mods_string list-of-ints -> ''", R._elem_mods_string("[1,2,3]") == ""))
    out.append(("_elem_mods_string scalar -> ''", R._elem_mods_string("42") == ""))

    # B-F2: a contrast TSV with no p_adj_bh column must not abort the section
    d = tempfile.mkdtemp()
    pd.DataFrame({"contrast": ["z_vs_m"], "feature": ["F"], "gene_name": ["G"], "delta": [0.3]}).to_csv(
        os.path.join(d, "x_z_vs_m_isoform_usage_diffs.tsv"), sep="\t", index=False)
    out.append(("build_between_conditions_section w/o p_adj_bh",
                _no_crash(R.build_between_conditions_section, d, 10)))

    # ---- F7-F10 audit (schema-drift crashes + silent total-loss paths) ----
    # F2/F3: _scan_perff_by_allele must not silently return {} when target_modified is missing or
    # when 'usable' is float-encoded ('1.0' from a NaN-forced float64 column).
    def _scan(usable, tmod_col=True):
        dd = tempfile.mkdtemp(); mp = os.path.join(dd, "m.tsv"); sp = os.path.join(dd, "s.tsv")
        cols = {"sample": ["S"] * 3, "qname": ["r0", "r1", "r2"], "chrom": ["chr1"] * 3, "start0": [100] * 3,
                "target_mod_code": ["a"] * 3, "ZN": [5, 5, 5], "usable": usable}
        if tmod_col:
            cols["target_modified"] = [1, 0, 1]
        pd.DataFrame(cols).to_csv(mp, sep="\t", index=False)
        pd.DataFrame({"sample": ["S"] * 3, "qname": ["r0", "r1", "r2"], "snp_id": ["snp1"] * 3,
                      "allele_class": ["ref", "ref", "alt"]}).to_csv(sp, sep="\t", index=False)
        return R._scan_perff_by_allele(mp, sp, [("GENEA", "chr1", 100, "a", "snp1")])
    out.append(("_scan_perff_by_allele: missing target_modified still scans", bool(_scan(["True"] * 3, tmod_col=False))))
    out.append(("_scan_perff_by_allele: float64 usable ('1.0') still scans", bool(_scan([1, 1, np.nan]))))

    # F4: build_polya_section must not crash when frag_df lacks median_tail
    out.append(("build_polya_section w/o median_tail",
                _no_crash(R.build_polya_section, pd.DataFrame({"ZT": ["G.T1"], "n_reads": [10]}), None, None, 10)))
    # F5/F6: build_snp_mechanism_section must not crash w/o p_adj_bh / direction_concordance
    out.append(("build_snp_mechanism_section w/o p_adj_bh & direction_concordance",
                _no_crash(R.build_snp_mechanism_section, pd.DataFrame({"positional_class": ["AT_MOD_BASE"], "gene_names": ["G"]}), 10)))
    # F7: build_classification_section must not crash w/o an 'event' column
    _cls = pd.DataFrame({"bucket": ["SHARED_LOCAL"], "gene_name": ["G"], "chrom": ["c"], "start0": [1],
                         "strand": ["+"], "mod_code": ["a"]})
    out.append(("build_classification_section w/o event col", _no_crash(R.build_classification_section, _cls, pd.DataFrame(), "", "", 0, 5)))
    # F1(report): a PRIVATE row that lives only in class_df surfaces; one in both frames isn't double-shown
    _priv = pd.DataFrame({"bucket": ["PRIVATE"], "event": ["ALT_LAST_EXON"], "direction": ["PROXIMAL_HIGHER"],
                          "gene_name": ["INBOTH"], "chrom": ["c"], "start0": [1], "strand": ["+"], "mod_code": ["a"]})
    _cls2 = pd.DataFrame({"bucket": ["PRIVATE", "PRIVATE"], "event": ["SKIPPED_EXON", "SKIPPED_EXON"],
                          "direction": ["", ""], "gene_name": ["ONLYCLASS", "INBOTH"], "chrom": ["c", "c"],
                          "start0": [9, 1], "strand": ["+", "+"], "mod_code": ["a", "a"]})
    _h = R.build_classification_section(_cls2, _priv, "", "", 0, 20)
    out.append(("class-only PRIVATE row surfaces (not silently dropped)", "ONLYCLASS" in _h))
    out.append(("PRIVATE row in both frames shown once (no double-count)", _h.count("INBOTH") == 1))

    return out


def main():
    results = checks()
    n_pass = sum(1 for _, ok in results if ok)
    for name, ok in results:
        print(f"  [{'PASS' if ok else '**FAIL**'}] {name}")
    print(f"\nreport robustness: {n_pass}/{len(results)} checks passed")
    sys.exit(0 if n_pass == len(results) else 1)


if __name__ == "__main__":
    main()
