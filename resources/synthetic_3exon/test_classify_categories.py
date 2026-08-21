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
# the differential-site 1-based base, which zn is hi (high m6A) vs lo, and the expected
# (bucket, event) the tree classifier (classify_tree) should return. All cases are + strand.
CASES = [
    # name, isoforms, site(1based), hi, lo, expected_bucket, expected_event
    ("PRIVATE_SKIPPED_EXON",        # base in a cassette exon present in hi, skipped (intronic) in lo
     {1: ([(1000, 1200), (2000, 2200), (3000, 3500)], 3500, 100),
      2: ([(1000, 1200), (3000, 3500)], 3500, 100)},
     2100, 1, 2, "PRIVATE", "SKIPPED_EXON"),

    ("PRIVATE_ALT_LAST_EXON",       # base in hi's distal terminal exon, absent from lo entirely
     {1: ([(1000, 1200), (2000, 2500), (5000, 6000)], 6000, 100),
      2: ([(1000, 1200), (2000, 3500)], 3500, 100)},
     5500, 1, 2, "PRIVATE", "ALT_LAST_EXON"),

    ("SHARED_LOCAL_ALT_ACCEPTOR",   # base's exon shares its donor but uses a different acceptor (5' start)
     {1: ([(1000, 1200), (2000, 2500), (3000, 3500)], 3500, 100),
      2: ([(1000, 1200), (2100, 2500), (3000, 3500)], 3500, 100)},
     2300, 1, 2, "SHARED_LOCAL", "ALT_ACCEPTOR"),

    ("SHARED_LOCAL_ALT_DONOR",      # base's exon shares its acceptor but uses a different donor (3' end)
     {1: ([(1000, 1200), (2000, 2500), (3000, 3500)], 3500, 100),
      2: ([(1000, 1200), (2000, 2400), (3000, 3500)], 3500, 100)},
     2200, 1, 2, "SHARED_LOCAL", "ALT_DONOR"),

    ("SHARED_LOCAL_ALT_POLYA_SITE", # both in the last exon, same acceptor, different poly(A) site
     {1: ([(1000, 1200), (2000, 2500)], 2500, 100),
      2: ([(1000, 1200), (2000, 3000)], 3000, 100)},
     2300, 1, 2, "SHARED_LOCAL", "ALT_POLYA_SITE"),

    ("SHARED_LOCAL_IPA_EXTENSION",  # hi reads into the intron and polyadenylates there (IPA); lo splices on
     {1: ([(1000, 1200), (2000, 2600)], 2600, 100),
      2: ([(1000, 1200), (2000, 2300), (3000, 3500)], 3500, 100)},
     2100, 1, 2, "SHARED_LOCAL", "IPA_EXTENSION"),

    ("SHARED_LOCAL_RETAINED_INTRON", # co-terminal; hi retains a 3'UTR intron that lo splices out
     {1: ([(1000, 1200), (2000, 3000)], 3000, 100),
      2: ([(1000, 1200), (2000, 2300), (2600, 3000)], 3000, 100)},
     2150, 1, 2, "SHARED_LOCAL", "RETAINED_INTRON"),

    ("SHARED_LOCAL_ALT_DONOR_INTERNAL", # base in an INTERNAL exon whose donor sits near (not at) the
     {1: ([(1000, 1200), (2000, 2480), (2490, 2500)], 2500, 100),   # tes; a real terminal exon follows,
      2: ([(1000, 1200), (2000, 2400), (2490, 2700)], 2700, 100)},  # so this is ALT_DONOR, not IPA
     2200, 1, 2, "SHARED_LOCAL", "ALT_DONOR"),

    ("SHARED_DISTAL_SPLICING",      # base's exon identical, SAME 3' end (tes 3500); forms differ only in the 5' exon
     {1: ([(1000, 1200), (2000, 2200), (3000, 3500)], 3500, 100),
      2: ([(1500, 1700), (2000, 2200), (3000, 3500)], 3500, 100)},
     3200, 1, 2, "SHARED_DISTAL", "DISTAL_SPLICING"),

    ("SHARED_DISTAL_APA",           # base's exon (middle) identical in both; forms differ only in 3' end / poly(A) site
     {1: ([(1000, 1200), (2000, 2200), (3000, 3500)], 3500, 100),
      2: ([(1000, 1200), (2000, 2200), (3000, 4200)], 4200, 100)},
     2100, 1, 2, "SHARED_DISTAL", "DISTAL_APA"),

    ("UNEXPLAINABLE_INTRON_READ",   # the higher-m6A form does not structurally contain the base (intronic)
     {1: ([(1000, 1200), (3000, 3500)], 3500, 100),
      2: ([(1000, 1200), (2000, 2200), (3000, 3500)], 3500, 100)},
     2000, 1, 2, "UNEXPLAINABLE", "INTRON_READ_ARTIFACT"),
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
    for name, isos, site1, hi, lo, _bkt, _evt in CASES:
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
        df = pd.read_csv(out, sep="\t").set_index("gene_name")
        got = {g: (df.loc[g, "bucket"], df.loc[g, "event"]) for g in df.index}

    n_pass = n_fail = 0
    print(f"  {'case':<30} {'expected (bucket/event)':<34} {'got':<34} ok")
    for name, _isos, _s, _hi, _lo, exp_b, exp_e in CASES:
        g = got.get(name, ("<missing>", ""))
        ok = (g == (exp_b, exp_e))
        n_pass += ok
        n_fail += (not ok)
        print(f"  {name:<30} {exp_b+'/'+exp_e:<34} {g[0]+'/'+str(g[1]):<34} {'PASS' if ok else '**FAIL**'}")
    print(f"\nclassify tree taxonomy: {n_pass}/{len(CASES)} cases correct"
          + ("" if not n_fail else f"  ({n_fail} FAILED)"))
    sys.exit(1 if n_fail else 0)


if __name__ == "__main__":
    main()
