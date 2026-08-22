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
