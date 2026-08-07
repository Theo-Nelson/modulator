#!/usr/bin/env python3
"""
Unit-validate the FULL classify_diff_sites.py structural taxonomy.

The end-to-end synthetic run only reaches the classify categories whose
differential site is covered by >=2 isoforms (shared-coverage categories:
SHARED_TERMINAL_EXON, TANDEM_APA, EJC_SPLICING...). The "private-exon"
categories (CASSETTE_EXON, INTRONIC_POLYADENYLATION, INTERGENIC_TERMINAL_EXON,
LAST_EXON_DISTAL_ONLY) require the site to be intronic/absent in the low
isoform, so real reads never produce a 2-isoform test there. To validate those,
this test drives classify_diff_sites.py directly: it hand-writes a GTF + a
differential-site table (with per_transcript_json) engineered so each site lands
in a specific category, runs the REAL classifier, and checks the mechanism it
returns (the `structural_category` column = STRUCTURAL_OF[raw 14-label]).

Usage:  <modulator-env>/bin/python resources/synthetic_3exon/test_classify_categories.py
"""
from __future__ import annotations
import json
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd

CHROM = "chrUnit"
SCRIPT = Path(__file__).resolve().parents[2] / "workflow" / "scripts" / "classify_diff_sites.py"

# Each case: gene -> isoforms {zn: (exons[list of 1-based (s,e)], tes, read_support)},
# the differential-site 1-based base, which zn is hi (high m6A) vs lo, and the
# expected STRUCTURAL category (the mapped mechanism the classifier outputs).
CASES = [
    # name, isoforms, site(1based), hi, lo, expected_structural
    ("CASSETTE_EXON",              # exon-skip, site in skipped cassette
     {1: ([(1000, 1200), (2000, 2200), (3000, 3500)], 3500, 40),
      2: ([(1000, 1200), (3000, 3500)], 3500, 200)},           # skip = anchor (higher rs)
     2100, 1, 2, "CASSETTE_EXON"),

    ("INTRONIC_POLYADENYLATION_UNIQUE",  # hi is an IPA form, site in its intron-derived terminal exon
     {1: ([(1000, 1200), (2000, 2200), (3000, 3500)], 3500, 200),
      2: ([(1000, 1200), (2000, 2600)], 2600, 60)},
     2400, 2, 1, "INTRONIC_POLYADENYLATION"),

    ("INTRONIC_POLYADENYLATION_SHARED_EJC",  # IPA form terminalizes a base internal to the anchor
     {1: ([(1000, 1200), (2000, 2200), (3000, 3500)], 3500, 200),
      2: ([(1000, 1200), (2000, 2400)], 2400, 60)},
     2100, 2, 1, "INTRONIC_POLYADENYLATION"),

    ("EJC_SPLICING",               # read-through hi vs spliced lo, junction near site in lo
     {1: ([(1000, 1200), (2000, 3000)], 3000, 60),
      2: ([(1000, 1200), (2000, 2400), (2550, 3000)], 3000, 200)},
     2600, 1, 2, "EJC_SPLICING"),

    ("TANDEM_APA_distal_only",     # site only in the distal isoform's extended 3'UTR
     {1: ([(1000, 1200), (2000, 3000)], 3000, 200),
      2: ([(1000, 1200), (2000, 2500)], 2500, 60)},
     2800, 1, 2, "TANDEM_APA"),

    ("TANDEM_APA_proximal_favored",  # same acceptor, different TES, proximal isoform hi
     {1: ([(1000, 1200), (2000, 2500)], 2500, 100),
      2: ([(1000, 1200), (2000, 3000)], 3000, 100)},
     2300, 1, 2, "TANDEM_APA"),

    ("TANDEM_APA_distal_favored",    # same geometry, distal isoform hi
     {1: ([(1000, 1200), (2000, 2500)], 2500, 100),
      2: ([(1000, 1200), (2000, 3000)], 3000, 100)},
     2300, 2, 1, "TANDEM_APA"),

    ("ALTERNATIVE_LAST_EXON",      # different (overlapping) last-exon acceptors
     {1: ([(1000, 1200), (2000, 2600)], 2600, 100),
      2: ([(1000, 1200), (2100, 2600)], 2600, 100)},
     2400, 1, 2, "ALTERNATIVE_LAST_EXON"),

    ("INTERGENIC_TERMINAL_EXON",   # hi terminal exon far downstream & disjoint
     {1: ([(1000, 1200), (2000, 2500), (5000, 6000)], 6000, 200),
      2: ([(1000, 1200), (2000, 3500)], 3500, 60)},
     5500, 1, 2, "INTERGENIC_TERMINAL_EXON"),

    ("SHARED_TERMINAL_EXON",       # same terminal exon + same TES, differ only 5'
     {1: ([(1000, 1200), (2000, 2200), (3000, 3500)], 3500, 100),
      2: ([(1500, 1700), (2000, 2200), (3000, 3500)], 3500, 100)},
     3200, 1, 2, "SHARED_TERMINAL_EXON"),

    ("SHARED_INTERNAL_EXON",       # site in a shared internal exon, symmetric junctions
     {1: ([(1000, 1200), (2000, 2600), (3000, 3500)], 3500, 100),
      2: ([(1000, 1200), (2000, 2600), (3000, 3200)], 3200, 100)},
     2300, 1, 2, "SHARED_INTERNAL_EXON"),

    ("UNEXPLAINED",                # terminal in hi & anchor, internal in lo, no nearby junction
     {1: ([(1000, 1200), (2000, 3500)], 3500, 100),
      2: ([(1000, 1200), (2000, 3000)], 3000, 100),
      3: ([(1000, 1200), (2000, 2700), (2900, 3200)], 3200, 100)},
     2500, 2, 3, "UNEXPLAINED"),

    ("ARTIFACT",                   # highest-m6A isoform does not contain the base (intronic)
     {1: ([(1000, 1200), (3000, 3500)], 3500, 100),
      2: ([(1000, 1200), (2000, 2200), (3000, 3500)], 3500, 100)},
     2000, 1, 2, "ARTIFACT"),
]


def write_gtf(path: Path):
    with open(path, "w") as fh:
        for name, isos, *_ in CASES:
            for zn, (exons, tes, rs) in isos.items():
                gstart = min(s for s, _ in exons)
                gend = max(e for _, e in exons)
                attrs = (f'gene_id "{name}"; ref_gene_name "{name}"; '
                         f'transcript_id "{name}.T{zn}"; zn_index "{zn}"; '
                         f'tes "{tes}"; read_support "{rs}";')
                fh.write(f"{CHROM}\tunit\ttranscript\t{gstart}\t{gend}\t.\t+\t.\t{attrs}\n")
                for (s, e) in exons:
                    fh.write(f"{CHROM}\tunit\texon\t{s}\t{e}\t.\t+\t.\t{attrs}\n")


def write_diff(path: Path):
    rows = []
    for name, isos, site1, hi, lo, _exp in CASES:
        start0 = site1 - 1
        per_tx = []
        for zn in isos:
            frac = 0.90 if zn == hi else 0.10 if zn == lo else 0.50
            per_tx.append({"ZN": zn, "Ncov": 100, "Nmod": int(round(frac * 100)), "frac": frac})
        rows.append({
            "gene_name": name, "mod_code": "a", "chrom": CHROM,
            "start0": start0, "end0": site1, "strand": "+",
            "n_tx_tested": len(isos), "test_name": "unit", "stat_name": "unit",
            "stat_value": 0.0, "p_value": 1e-6,
            "effect_max_abs_frac_diff": 0.80, "per_transcript_json": json.dumps(per_tx),
            "p_adj_bh": 1e-4,
        })
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)


def main():
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        gtf, diff, out = td / "unit.gtf", td / "diff.tsv", td / "classified.tsv"
        write_gtf(gtf)
        write_diff(diff)
        subprocess.run([sys.executable, str(SCRIPT), "--diff-tsv", str(diff),
                        "--gtf", str(gtf), "--out-tsv", str(out)],
                       check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        got = pd.read_csv(out, sep="\t").set_index("gene_name")["structural_category"].to_dict()

    n_pass = n_fail = 0
    print(f"  {'case':<34} {'expected':<26} {'got':<26} ok")
    for name, _isos, _s, _hi, _lo, exp in CASES:
        g = got.get(name, "<missing>")
        ok = (g == exp)
        n_pass += ok
        n_fail += (not ok)
        print(f"  {name:<34} {exp:<26} {g:<26} {'PASS' if ok else '**FAIL**'}")
    print(f"\nclassify taxonomy: {n_pass}/{len(CASES)} categories correct"
          + ("" if not n_fail else f"  ({n_fail} FAILED)"))
    sys.exit(1 if n_fail else 0)


if __name__ == "__main__":
    main()
