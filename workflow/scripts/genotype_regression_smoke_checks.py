#!/usr/bin/env python3

import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parent


def run(cmd):
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise SystemExit(f"Command failed: {' '.join(cmd)}\nSTDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}")
    return proc


def main():
    with tempfile.TemporaryDirectory(prefix="genotype_smoke_") as tmpdir:
        tmp = Path(tmpdir)
        off_context_mod_site = "chr1:300-301:+:a"

        snp_rows = []
        mod_rows = []
        reads = []
        for i in range(1, 9):
            q = f"ref_t1_{i}"
            reads.append(("S1", q, "TX1", "ref", "A", "A"))
        for i in range(1, 3):
            q = f"ref_t2_{i}"
            reads.append(("S1", q, "TX2", "ref", "A", "A"))
        for i in range(1, 3):
            q = f"alt_t1_{i}"
            reads.append(("S1", q, "TX1", "alt", "G", "G"))
        for i in range(1, 9):
            q = f"alt_t2_{i}"
            reads.append(("S1", q, "TX2", "alt", "G", "G"))

        hap_alleles = {}
        for idx, (sample, qname, zt, allele_class, obs1, obs2) in enumerate(reads, start=1):
            for snp_id, pos1, ref, alt, obs in [
                ("chr1:101:A>G", 101, "A", "G", obs1),
                ("chr1:151:C>T", 151, "C", "T", obs2),
            ]:
                snp_rows.append({
                    "sample": sample,
                    "qname": qname,
                    "snp_id": snp_id,
                    "chrom": "chr1",
                    "pos1": pos1,
                    "start0": pos1 - 1,
                    "end0": pos1,
                    "ref": ref,
                    "alt": alt,
                    "observed_base": obs,
                    "allele_class": "ref" if obs == ref else "alt",
                    "baseq": 40,
                    "mapq": 60,
                    "strand": "+",
                    "ZT": zt,
                    "ZG": 1,
                    "ZN": 1 if zt == "TX1" else 2,
                    "ZM": 1,
                    "gene_names": "GENE1",
                    "gene_ids": "GENE1",
                    "metagene_indices": "1",
                })

            target_modified = 0
            if allele_class == "alt" and zt == "TX2":
                target_modified = 1
            elif allele_class == "alt" and idx % 2 == 0:
                target_modified = 1
            mod_rows.append({
                "sample": sample,
                "qname": qname,
                "mod_site_id": "chr1:200-201:+:a",
                "chrom": "chr1",
                "start0": 200,
                "end0": 201,
                "strand": "+",
                "target_mod_code": "a",
                "call_code": "a" if target_modified else "-",
                "state_detail": "modified" if target_modified else "canonical",
                "target_modified": target_modified,
                "call_prob": 0.99,
                "canonical_base": "A",
                "modified_primary_base": "A",
                "fail": False,
                "within_alignment": True,
                "gene_id": "GENE1",
                "gene_name": "GENE1",
                "metagene_index": "1",
                "ZT": zt,
                "ZG": 1,
                "ZN": 1 if zt == "TX1" else 2,
                "ZM": 1,
                "assigned": True,
                "assignment_gene_id": "GENE1",
                "assignment_gene_name": "GENE1",
                "assignment_metagene_index": "1",
                "usable": True,
            })
            mod_rows.append({
                "sample": sample,
                "qname": qname,
                "mod_site_id": off_context_mod_site,
                "chrom": "chr1",
                "start0": 300,
                "end0": 301,
                "strand": "+",
                "target_mod_code": "a",
                "call_code": "a" if idx % 3 == 0 else "-",
                "state_detail": "modified" if idx % 3 == 0 else "canonical",
                "target_modified": 1 if idx % 3 == 0 else 0,
                "call_prob": 0.99,
                "canonical_base": "A",
                "modified_primary_base": "A",
                "fail": False,
                "within_alignment": True,
                "gene_id": "GENE_OFF",
                "gene_name": "GENE_OFF",
                "metagene_index": "99",
                "ZT": zt,
                "ZG": 1,
                "ZN": 1 if zt == "TX1" else 2,
                "ZM": 1,
                "assigned": True,
                "assignment_gene_id": "GENE_OFF",
                "assignment_gene_name": "GENE_OFF",
                "assignment_metagene_index": "99",
                "usable": True,
            })

        snp_path = tmp / "molecule_snps.tsv"
        mod_path = tmp / "molecule_mods.tsv"
        pd.DataFrame(snp_rows).to_csv(snp_path, sep="\t", index=False)
        pd.DataFrame(mod_rows).to_csv(mod_path, sep="\t", index=False)

        snp_tx_out = tmp / "snp_tx.tsv"
        snp_mod_out = tmp / "snp_mod.tsv"
        hap_blocks = tmp / "hap_blocks.tsv"
        hap_mols = tmp / "hap_molecules.tsv"
        hap_tx_out = tmp / "hap_tx.tsv"
        hap_mod_out = tmp / "hap_mod.tsv"

        run([sys.executable, str(ROOT / "test_snp_transcript_assoc.py"), "--molecule-snps", str(snp_path), "--out-tsv", str(snp_tx_out), "--min-allele-reads", "2", "--min-transcript-reads", "2"])
        run([sys.executable, str(ROOT / "test_snp_mod_assoc.py"), "--molecule-snps", str(snp_path), "--molecule-mods", str(mod_path), "--out-tsv", str(snp_mod_out), "--min-allele-reads", "2", "--min-total-reads", "4"])
        run([sys.executable, str(ROOT / "build_haplotype_blocks.py"), "--molecule-snps", str(snp_path), "--out-blocks-tsv", str(hap_blocks), "--out-molecules-tsv", str(hap_mols), "--min-alt-reads", "2", "--min-cocover-reads", "2", "--max-block-snps", "4", "--min-haplotype-reads", "2"])
        run([sys.executable, str(ROOT / "test_haplotype_associations.py"), "--molecule-haplotypes", str(hap_mols), "--molecule-mods", str(mod_path), "--out-haplotype-transcript", str(hap_tx_out), "--out-haplotype-mod", str(hap_mod_out), "--min-haplotype-reads", "2", "--min-transcript-reads", "2", "--min-total-reads", "4"])

        snp_tx = pd.read_csv(snp_tx_out, sep="\t")
        snp_mod = pd.read_csv(snp_mod_out, sep="\t")
        hap_blocks_df = pd.read_csv(hap_blocks, sep="\t")
        hap_tx = pd.read_csv(hap_tx_out, sep="\t")
        hap_mod = pd.read_csv(hap_mod_out, sep="\t")

        if snp_tx.empty:
            raise AssertionError("Expected non-empty SNP to transcript associations.")
        if snp_mod.empty:
            raise AssertionError("Expected non-empty SNP to mod associations.")
        if hap_blocks_df.empty:
            raise AssertionError("Expected at least one haplotype block.")
        if hap_tx.empty:
            raise AssertionError("Expected non-empty haplotype to transcript associations.")
        if hap_mod.empty:
            raise AssertionError("Expected non-empty haplotype to mod associations.")

        if float(snp_tx.iloc[0]["effect_max_abs_tx_frac_diff"]) <= 0.0:
            raise AssertionError("Expected positive SNP transcript effect size.")
        if float(snp_mod.iloc[0]["effect_abs_delta_mod_frac"]) <= 0.0:
            raise AssertionError("Expected positive SNP mod effect size.")
        if off_context_mod_site in set(snp_mod.get("mod_site_id", [])):
            raise AssertionError("Off-context mod site should not appear in SNP-mod associations.")
        if off_context_mod_site in set(hap_mod.get("mod_site_id", [])):
            raise AssertionError("Off-context mod site should not appear in haplotype-mod associations.")

    print("genotype_regression_smoke_checks: OK")


if __name__ == "__main__":
    main()
