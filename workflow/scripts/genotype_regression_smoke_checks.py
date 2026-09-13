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

    # ---- A9: a SNP whose exon annotation names metagene 1 only, observed on reads ASSIGNED to metagene
    # 2 (e.g. the SNP lies in a retained intron / 5' region absent from metagene-2 fragmentform models).
    # The mod calls on those same reads are keyed MG:2 (assignment metagene), so the pair must be found.
    with tempfile.TemporaryDirectory(prefix="genotype_smoke_a9_") as tmpdir:
        tmp = Path(tmpdir)
        snp_rows, mod_rows = [], []
        for i in range(1, 17):
            q = f"r{i}"
            allele = "alt" if i % 2 == 0 else "ref"
            modified = 1 if allele == "alt" else 0           # perfectly allele-linked modification
            snp_rows.append({"sample": "S1", "qname": q, "snp_id": "chr1:501:A>G", "chrom": "chr1", "pos1": 501,
                             "start0": 500, "end0": 501, "ref": "A", "alt": "G",
                             "observed_base": "G" if allele == "alt" else "A", "allele_class": allele,
                             "baseq": 40, "mapq": 60, "strand": "+", "ZT": "GENE2.GENE2.G2.T1", "ZG": 2, "ZN": 1,
                             "ZM": 2, "gene_names": "GENE1", "gene_ids": "GENE1", "metagene_indices": "1"})
            mod_rows.append({"sample": "S1", "qname": q, "mod_site_id": "chr1:900-901:+:a", "chrom": "chr1",
                             "start0": 900, "end0": 901, "strand": "+", "target_mod_code": "a",
                             "call_code": "a" if modified else "-", "state_detail": "modified" if modified else "canonical",
                             "target_modified": modified, "call_prob": 0.99, "canonical_base": "A",
                             "modified_primary_base": "A", "fail": False, "within_alignment": True,
                             "gene_id": "GENE2", "gene_name": "GENE2", "metagene_index": "2",
                             "ZT": "GENE2.GENE2.G2.T1", "ZG": 2, "ZN": 1, "ZM": 2, "assigned": True,
                             "assignment_gene_id": "GENE2", "assignment_gene_name": "GENE2",
                             "assignment_metagene_index": "2", "usable": True})
        snp_path = tmp / "molecule_snps.tsv"; mod_path = tmp / "molecule_mods.tsv"
        pd.DataFrame(snp_rows).to_csv(snp_path, sep="\t", index=False)
        pd.DataFrame(mod_rows).to_csv(mod_path, sep="\t", index=False)
        out = tmp / "snp_mod.tsv"
        run([sys.executable, str(ROOT / "test_snp_mod_assoc.py"), "--molecule-snps", str(snp_path),
             "--molecule-mods", str(mod_path), "--out-tsv", str(out), "--min-allele-reads", "2", "--min-total-reads", "4"])
        res = pd.read_csv(out, sep="\t")
        if res.empty or "chr1:900-901:+:a" not in set(res["mod_site_id"]):
            raise AssertionError("A9: SNP on reads assigned to metagene 2 must pair with the mod calls on those reads "
                                 "even though the SNP's exon annotation names metagene 1 only.")
        if float(res.iloc[0]["effect_abs_delta_mod_frac"]) < 0.9:
            raise AssertionError("A9: the perfectly allele-linked modification should show a ~1.0 effect.")

    # ---- depth cap + per-window checkpoint in discover_candidate_snps: a 3,000-read het locus and a
    # 20-read het locus. Cap off must be exact; a cap must bound the deep locus, leave the shallow one
    # untouched and keep the allele fraction ~0.5; a checkpointed window must be reloaded on a re-run.
    _check_depth_cap()

    print("genotype_regression_smoke_checks: OK")


def _check_depth_cap():
    import os
    import pickle
    import random
    import pysam
    with tempfile.TemporaryDirectory(prefix="genotype_smoke_depth_") as tmpdir:
        tmp = Path(tmpdir)
        random.seed(7)
        ref = "".join(random.choice("ACGT") for _ in range(400))
        ref = ref[:149] + "A" + ref[150:]
        ref = ref[:249] + "C" + ref[250:]
        (tmp / "ref.fa").write_text(">chrT\n" + ref + "\n")
        pysam.faidx(str(tmp / "ref.fa"))
        attrs = ('gene_id "G1"; transcript_id "T1"; gene_name "G1"; ref_gene_name "G1"; zt_label "G1.G1.G1.T1"; '
                 'gene_index "1"; transcript_index "1"; metagene_index "1"; zn_index "1";')
        (tmp / "genes.gtf").write_text(f"chrT\tt\ttranscript\t1\t400\t.\t+\t.\t{attrs}\nchrT\tt\texon\t1\t400\t.\t+\t.\t{attrs}\n")
        hdr = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chrT", "LN": 400}]}

        def mk(name, n_deep, n_shallow):
            reads = []
            for i in range(n_deep):
                seq = list(ref[100:200])
                if i % 2 == 0:
                    seq[49] = "G"
                reads.append((100, "".join(seq), f"deep{i}"))
            for i in range(n_shallow):
                seq = list(ref[220:280])
                if i % 2 == 0:
                    seq[29] = "T"
                reads.append((220, "".join(seq), f"shal{i}"))
            with pysam.AlignmentFile(str(name), "wb", header=hdr) as out:
                for start, seq, qn in sorted(reads):
                    a = pysam.AlignedSegment()
                    a.query_name = qn; a.query_sequence = seq; a.flag = 0; a.reference_id = 0
                    a.reference_start = start; a.mapping_quality = 60; a.cigar = [(0, len(seq))]
                    a.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
                    out.write(a)
            pysam.index(str(name))

        mk(tmp / "s1.bam", 3000, 20)
        mk(tmp / "s2.bam", 1500, 20)
        base = [sys.executable, str(ROOT / "discover_candidate_snps.py"), "--bams", str(tmp / "s1.bam"), str(tmp / "s2.bam"),
                "--reference-fa", str(tmp / "ref.fa"), "--gtf", str(tmp / "genes.gtf"), "--min-alt-reads", "4",
                "--min-total-cov", "8", "--jobs", "1", "--primary-only", "--verbose"]
        run(base + ["--out-tsv", str(tmp / "off.tsv"), "--max-depth", "0"])
        run(base + ["--out-tsv", str(tmp / "cap.tsv"), "--max-depth", "300", "--depth-chunk-bp", "200"])
        off = pd.read_csv(tmp / "off.tsv", sep="\t")
        cap = pd.read_csv(tmp / "cap.tsv", sep="\t")
        if list(off["snp_id"]) != ["chrT:150:A>G", "chrT:250:C>T"] or list(cap["snp_id"]) != list(off["snp_id"]):
            raise AssertionError(f"depth cap: candidate set changed: {list(off['snp_id'])} vs {list(cap['snp_id'])}")
        if int(off.loc[0, "total_cov"]) != 4500 or int(off.loc[1, "total_cov"]) != 40:
            raise AssertionError("depth cap off must give exact counts (4500 / 40)")
        if not (int(cap.loc[0, "total_cov"]) <= 2 * 300 * 1.1) or abs(float(cap.loc[0, "alt_frac"]) - 0.5) > 0.06:
            raise AssertionError(f"depth cap must bound the deep locus near 2x300 reads with alt_frac ~0.5 "
                                 f"(got {int(cap.loc[0, 'total_cov'])}, {cap.loc[0, 'alt_frac']})")
        if int(cap.loc[1, "total_cov"]) != 40 or float(cap.loc[1, "alt_frac"]) != 0.5:
            raise AssertionError("depth cap must leave a shallow locus untouched")
        # checkpoint: plant a completed-window checkpoint, replace s1 by an empty BAM, rerun -> capped result
        sys.path.insert(0, str(ROOT))
        import discover_candidate_snps as dcs
        _exon, merged = dcs.load_gtf_exons(str(tmp / "genes.gtf"))
        chrom, ivs = next(dcs.iter_window_shards(merged, 1_000_000))
        shards = tmp / "shards"; shards.mkdir()
        drops = {"low_total_cov": 0, "low_alt_reads": 0, "low_alt_frac": 0, "high_alt_frac": 0, "multiallelic": 0,
                 "chunks_subsampled": 0, "max_reads_per_chunk": 0}
        with open(dcs._shard_ckpt_path(str(shards), chrom, ivs), "wb") as fh:
            pickle.dump((cap.to_dict("records"), 0, drops), fh)
        os.remove(tmp / "s1.bam"); os.remove(tmp / "s1.bam.bai")
        mk(tmp / "s1.bam", 0, 0)
        run(base + ["--out-tsv", str(tmp / "ck.tsv"), "--max-depth", "300", "--depth-chunk-bp", "200", "--shard-dir", str(shards)])
        ck = pd.read_csv(tmp / "ck.tsv", sep="\t")
        if list(ck["snp_id"]) != list(cap["snp_id"]) or int(ck.loc[0, "total_cov"]) != int(cap.loc[0, "total_cov"]):
            raise AssertionError("per-window checkpoint was not reused on re-run")
        if shards.exists():
            raise AssertionError("shard dir must be removed after a successful run")


if __name__ == "__main__":
    main()
