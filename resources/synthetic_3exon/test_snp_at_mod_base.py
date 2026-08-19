#!/usr/bin/env python3
"""
Unit-validate the SNP-at-modified-base classifier (find_snp_at_mod_base.py) with
hand-crafted candidate-SNP + modification-site tables covering the diagnostic
cases (self-reporting A-to-I / pseudouridine, modified-base ablation), on both
strands. The 5-gene demo rarely contains clean cases, so this exercises the logic
directly.

Usage: <modulator-env>/bin/python resources/synthetic_3exon/test_snp_at_mod_base.py
"""
from __future__ import annotations
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd

SCRIPT = Path(__file__).resolve().parents[2] / "workflow" / "scripts" / "find_snp_at_mod_base.py"

# name, mod (start0, strand, mod_code), snp (ref, alt, at same pos0), expected class
CASES = [
    ("AtoI_selfreport_plus",  (100, "+", "17596"), ("A", "G"), "EDITING_SELF_REPORT"),
    ("AtoI_selfreport_minus", (200, "-", "17596"), ("T", "C"), "EDITING_SELF_REPORT"),  # rc(A>G)
    ("pseU_selfreport_plus",  (300, "+", "17802"), ("T", "C"), "PSEU_SELF_REPORT"),
    ("pseU_selfreport_minus", (350, "-", "17802"), ("A", "G"), "PSEU_SELF_REPORT"),      # rc(T>C)
    ("m6A_ablated_plus",      (400, "+", "a"),     ("A", "C"), "MOD_BASE_ABLATED"),
    ("m6A_ablated_minus",     (450, "-", "a"),     ("T", "G"), "MOD_BASE_ABLATED"),       # rc(A>C)
    ("m5C_ablated",           (500, "+", "m"),     ("C", "T"), "MOD_BASE_ABLATED"),
    ("inosine_ablated_notG",  (600, "+", "17596"), ("A", "T"), "MOD_BASE_ABLATED"),
]
CHROM = "chrU"


def main():
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        snps, mods, out = td / "snps.tsv", td / "mods.tsv", td / "out.tsv"
        # candidate SNPs
        srows = []
        for name, (s0, strand, code), (ref, alt), _exp in CASES:
            srows.append(dict(snp_id=f"{CHROM}:{s0+1}:{ref}>{alt}", chrom=CHROM, pos1=s0 + 1,
                              ref=ref, alt=alt, alt_frac=0.5))
        pd.DataFrame(srows).to_csv(snps, sep="\t", index=False)
        # modification sites (minimal FILTERED_sites_long columns)
        mrows = []
        for name, (s0, strand, code), _snp, _exp in CASES:
            mrows.append(dict(sample="S1", ZN_transcript_index=1, chrom=CHROM, start0=s0,
                              end0=s0 + 1, strand=strand, mod_code=code, Nvalid_cov=50, Nmod=30,
                              frac_modified=0.6, gene_id="G", gene_name="GENE", Ncanonical=20,
                              Nother_mod=0, Ndelete=0, Nfail=0, Ndiff=0, Nnocall=0))
        pd.DataFrame(mrows).to_csv(mods, sep="\t", index=False)

        subprocess.run([sys.executable, str(SCRIPT), "--candidate-snps", str(snps),
                        "--mod-sites", str(mods), "--out-tsv", str(out)],
                       check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        got = pd.read_csv(out, sep="\t").set_index("pos0")["snp_at_mod_base_class"].to_dict()

    n_pass = n_fail = 0
    print(f"  {'case':<26} {'expected':<22} {'got':<22} ok")
    for name, (s0, _st, _code), _snp, exp in CASES:
        g = got.get(s0, "<missing>")
        ok = (g == exp)
        n_pass += ok; n_fail += (not ok)
        print(f"  {name:<26} {exp:<22} {g:<22} {'PASS' if ok else '**FAIL**'}")
    print(f"\nsnp-at-mod-base classifier: {n_pass}/{len(CASES)} correct"
          + ("" if not n_fail else f"  ({n_fail} FAILED)"))
    sys.exit(1 if n_fail else 0)


if __name__ == "__main__":
    main()
