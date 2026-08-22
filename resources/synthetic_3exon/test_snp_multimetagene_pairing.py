#!/usr/bin/env python3
"""Regression test for F1 (stress campaign #2 follow-up): a SNP overlapping >1 metagene must still
reach snp_mod_assoc. It used to collapse to a single CHR: context key that could never match the mod
side (which keys on MG:<metagene_index>), so 0/123 multi-metagene SNPs paired on real HG002 data.

Drives the REAL pairing function (_pairs_for_one_chrom) with a hand-built multi-metagene SNP + a mod
site in one of its metagenes, sharing reads, and asserts an association row is produced.
Run:  <modulator-env>/bin/python resources/synthetic_3exon/test_snp_multimetagene_pairing.py
"""
import os
import sys
import tempfile
from types import SimpleNamespace

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "workflow", "scripts"))
import test_snp_mod_assoc as T  # noqa: E402


def _write(tmp):
    reads = [f"r{i}" for i in range(6)]           # r0..r5 shared between the SNP and the mod site
    # mod site in metagene 5: r0-2 modified, r3-5 canonical
    mod_rows = []
    for i, r in enumerate(reads):
        mod_rows.append({"sample": "S1", "qname": r, "mod_site_id": "chr1:200-201:+:a",
                         "chrom": "chr1", "start0": 200, "end0": 201, "target_mod_code": "a",
                         "state_detail": "modified" if i < 3 else "canonical",
                         "gene_name": "GENEA", "metagene_index": 5, "usable": True})
    mod_p = os.path.join(tmp, "mods.tsv")
    pd.DataFrame(mod_rows).to_csv(mod_p, sep="\t", index=False)

    # ONE SNP that overlaps TWO metagenes (5 and 7) -> multi-metagene -> used to get CHR: context.
    # ref allele on r0,r1,r3 ; alt on r2,r4,r5 (both alleles carry modified + canonical reads).
    alleles = {"r0": "ref", "r1": "ref", "r3": "ref", "r2": "alt", "r4": "alt", "r5": "alt"}
    snp_rows = [{"sample": "S1", "qname": r, "snp_id": "chr1:150:A>G", "chrom": "chr1", "pos1": 150,
                 "start0": 149, "end0": 150, "allele_class": alleles[r],
                 "gene_names": "GENEA;GENEB", "metagene_indices": "5;7"} for r in reads]
    snp_p = os.path.join(tmp, "snps.tsv")
    pd.DataFrame(snp_rows).to_csv(snp_p, sep="\t", index=False)
    return mod_p, snp_p


def main():
    with tempfile.TemporaryDirectory() as tmp:
        mod_p, snp_p = _write(tmp)
        args = SimpleNamespace(min_allele_reads=2, min_total_reads=4, pseudocount=0.5,
                               test="fisher", molecule_mods=mod_p, molecule_snps=snp_p,
                               out_tsv=os.path.join(tmp, "out.tsv"))
        rows = T._pairs_for_one_chrom(mod_p, snp_p, args)
        got = [r for r in rows if r.get("snp_id") == "chr1:150:A>G"]
        checks = [
            ("multi-metagene SNP produces >=1 association row", len(got) >= 1),
            ("paired against the metagene-5 mod site",
             any(r.get("mod_site_id") == "chr1:200-201:+:a" for r in got)),
        ]
    n_pass = sum(1 for _, ok in checks if ok)
    for name, ok in checks:
        print(f"  [{'PASS' if ok else '**FAIL**'}] {name}")
    print(f"\nSNP multi-metagene pairing: {n_pass}/{len(checks)} checks passed")
    sys.exit(0 if n_pass == len(checks) else 1)


if __name__ == "__main__":
    main()
