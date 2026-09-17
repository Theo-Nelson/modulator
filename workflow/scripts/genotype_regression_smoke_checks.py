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

    # ---- parquet per-read tables: TSV -> parquet -> per-chromosome reads must equal the TSV reads
    # (values, NaN pattern, and read_csv-style dtypes), and the reverse conversion must round-trip.
    _check_parquet_roundtrip()

    # ---- between_conditions at scale: the streamed sample-filtered read must equal a whole-file read
    # (rows, order, dtypes incl. a mixed text/numeric column), and the parallel per-site fits must give
    # the same table as the sequential ones (threads=1 vs threads=3), for both site and per-transcript.
    _check_between_conditions()

    print("genotype_regression_smoke_checks: OK")


def _check_parquet_roundtrip():
    import numpy as np
    sys.path.insert(0, str(ROOT))
    from genotype_utils import open_chrom_table
    with tempfile.TemporaryDirectory(prefix="genotype_smoke_parquet_") as tmpdir:
        tmp = Path(tmpdir)
        rows = []
        for chrom in ("chr2", "chr10"):
            for i in range(6):
                rows.append({"sample": "S1", "qname": f"r{chrom}_{i}", "mod_site_id": f"{chrom}:{100+i}-{101+i}:+:a",
                             "chrom": chrom, "start0": 100 + i, "end0": 101 + i, "strand": "+", "target_mod_code": "a",
                             "call_code": "a" if i % 2 else "-", "state_detail": "modified" if i % 2 else "canonical",
                             "target_modified": i % 2, "call_prob": 0.5 + i / 20, "canonical_base": "A",
                             "modified_primary_base": "A", "fail": False, "within_alignment": True,
                             "gene_id": "G1", "gene_name": "G1", "metagene_index": 1 if chrom == "chr2" else np.nan,
                             "ZT": "G1.G1.G1.T1", "ZG": 1, "ZN": 1, "ZM": 1, "assigned": True,
                             "assignment_gene_id": "G1", "assignment_gene_name": "G1", "gene_index": 1,
                             "transcript_index": 1, "assignment_metagene_index": 1, "classification": "EXACT",
                             "usable": True})
        df = pd.DataFrame(rows)
        tsv = tmp / "t.tsv"; df.to_csv(tsv, sep="\t", index=False)
        run([sys.executable, str(ROOT / "convert_table_format.py"), "--in", str(tsv), "--out", str(tmp / "t.parquet")])
        run([sys.executable, str(ROOT / "convert_table_format.py"), "--in", str(tmp / "t.parquet"), "--out", str(tmp / "back.tsv")])
        a, b = open_chrom_table(str(tsv)), open_chrom_table(str(tmp / "t.parquet"))
        if a.chroms != b.chroms or a.header_cols != b.header_cols:
            raise AssertionError("parquet: chromosome index or columns differ from the TSV")
        for chrom in a.chroms:
            x, y = a.read(chrom), b.read(chrom)
            for c in x.columns:
                if str(x[c].dtype) != str(y[c].dtype) and c != "metagene_index":
                    raise AssertionError(f"parquet: dtype of {c} on {chrom}: {x[c].dtype} vs {y[c].dtype}")
                try:   # 1.0 (TSV float inference for an int column with NaN elsewhere) == 1 (typed int)
                    pd.testing.assert_series_equal(x[c].reset_index(drop=True), y[c].reset_index(drop=True),
                                                   check_dtype=False, check_names=False, check_exact=False)
                except AssertionError as e:
                    raise AssertionError(f"parquet: values of {c} on {chrom} differ: {e}") from None
        # a block larger than pyarrow's default 1,048,576-row group must not desynchronise the index
        from genotype_utils import ParquetBlockWriter
        import pyarrow.parquet as pq
        big = pd.DataFrame({"sample": "S", "qname": np.arange(1_300_000).astype(str), "chrom": "chrB",
                            "start0": np.arange(1_300_000), "target_modified": 0})
        w = ParquetBlockWriter(str(tmp / "big.parquet"), list(big.columns)); w.write("chrB", big); w.close()
        pf = pq.ParquetFile(str(tmp / "big.parquet"))
        got = open_chrom_table(str(tmp / "big.parquet")).read("chrB")
        if pf.metadata.num_row_groups != 1 or len(got) != 1_300_000 or int(got["start0"].iloc[-1]) != 1_299_999:
            raise AssertionError("parquet: a >1M-row block must be exactly one recorded row group")
        try:   # typed ints come back as "1" where the pandas-written TSV had "1.0": numerically equal
            _key = ["chrom", "start0", "qname"]   # the converter writes chromosomes in sorted order
            pd.testing.assert_frame_equal(
                pd.read_csv(tmp / "back.tsv", sep="\t").sort_values(_key).reset_index(drop=True),
                pd.read_csv(tsv, sep="\t").sort_values(_key).reset_index(drop=True),
                check_dtype=False, check_exact=False)
        except AssertionError as e:
            raise AssertionError(f"parquet -> tsv round trip differs: {e}") from None


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



def _check_between_conditions():
    import random
    import numpy as np
    sys.path.insert(0, str(ROOT))
    from genotype_utils import read_tsv_for_samples
    import diffstats
    with tempfile.TemporaryDirectory(prefix="between_conditions_smoke_") as tmpdir:
        tmp = Path(tmpdir)
        rng = random.Random(11)
        samples = [f"S{i}" for i in range(1, 7)]
        meta = pd.DataFrame({"sample": samples, "condition": ["ref"] * 3 + ["test"] * 3})
        meta.to_csv(tmp / "meta.tsv", sep="\t", index=False)
        rows = []
        for si in range(500):                      # 500 sites x up to 2 ZN partitions x 6 samples
            chrom = "chr1" if si < 400 else "chr2"
            start0 = 1000 + si * 7
            mod = rng.choice(["a", "m", "17802"])    # mixed text/numeric codes -> object column
            gene = rng.choice(["G1", "G2", ""])       # "" -> NaN after read_csv
            for zn in (1, 2)[: rng.choice((1, 2))]:
                p_ref = rng.uniform(0.05, 0.6)
                p_test = p_ref + rng.choice((0.0, 0.0, 0.25))
                for sample in samples:
                    if rng.random() < 0.1:
                        continue                      # uncovered (sample x site) -> NaN in the pivot
                    n = rng.randint(5, 80)
                    p = p_ref if sample in samples[:3] else p_test
                    k = sum(1 for _ in range(n) if rng.random() < p)
                    rows.append((sample, zn, chrom, start0, start0 + 1, "+", mod, n, k, gene))
        long = pd.DataFrame(rows, columns=["sample", "ZN_transcript_index", "chrom", "start0", "end0",
                                           "strand", "mod_code", "Nvalid_cov", "Nmod", "gene_name"])
        long.to_csv(tmp / "long.tsv", sep="\t", index=False)
        # (a) streamed sample-filtered read == whole-file read + filter (small chunks force many chunks,
        #     including chunks where mod_code parses as all-numeric or gene_name as all-empty)
        want = ["sample", "ZN_transcript_index", "chrom", "start0", "end0", "strand", "mod_code",
                "Nvalid_cov", "Nmod", "gene_name"]
        keep = {"S1", "S2", "S4", "S5"}
        whole = pd.read_csv(tmp / "long.tsv", sep="\t", low_memory=False, usecols=want)
        whole = whole[whole["sample"].astype(str).isin(keep)].reset_index(drop=True)
        streamed = read_tsv_for_samples(str(tmp / "long.tsv"), want, "sample", keep, chunksize=37)
        if list(streamed.columns) != list(whole.columns):
            raise AssertionError(f"streamed read: column order differs {list(streamed.columns)}")
        for c in want:
            if str(streamed[c].dtype) != str(whole[c].dtype):
                raise AssertionError(f"streamed read: dtype of {c} is {streamed[c].dtype}, whole-file {whole[c].dtype}")
        pd.testing.assert_frame_equal(streamed, whole)
        # (b) sequential vs parallel fits must be identical, row for row
        for by_tx in (False, True):
            outs = []
            for th in (1, 3):
                out = tmp / f"mod_{int(by_tx)}_{th}.tsv"
                run([sys.executable, str(ROOT / "test_condition_mod_diffs.py"), "--in-tsv", str(tmp / "long.tsv"),
                     "--sample-metadata", str(tmp / "meta.tsv"), "--out-tsv", str(out), "--test", "test",
                     "--reference", "ref", "--min-cov", "5", "--threads", str(th), "--chunk-rows", "101"]
                    + (["--by-transcript"] if by_tx else []))
                outs.append(out.read_bytes())
            if outs[0] != outs[1]:
                raise AssertionError(f"between_conditions: threads=1 and threads=3 tables differ (by_tx={by_tx})")
            n = len(pd.read_csv(tmp / f"mod_{int(by_tx)}_1.tsv", sep="\t"))
            if n < 100:
                raise AssertionError(f"between_conditions: expected >=100 tested sites, got {n} (by_tx={by_tx})")
        # (c) the tail test must not depend on string-hash order (it iterated set-valued groups, so
        #     Welch's t summed the replicate medians in PYTHONHASHSEED order -> last-digit p differences)
        tail_rows = []
        for si in range(60):
            zt = f"G{si}.G{si}.G1.T1"
            for sample in samples:
                base = 60 + si + (15 if sample in samples[3:] and si % 2 == 0 else 0)
                for r in range(12):
                    tail_rows.append((sample, f"r{si}_{sample}_{r}", base + rng.randint(-20, 20), zt, f"G{si}"))
        pd.DataFrame(tail_rows, columns=["sample", "qname", "tail_len", "ZT", "gene_name"]).to_csv(
            tmp / "tails.tsv", sep="\t", index=False)
        import os
        outs = []
        for seed in ("1", "2"):
            out = tmp / f"tail_{seed}.tsv"
            run_env = dict(os.environ, PYTHONHASHSEED=seed)
            proc = subprocess.run([sys.executable, str(ROOT / "test_condition_tail_diffs.py"), "--tail-tsv", str(tmp / "tails.tsv"),
                                   "--sample-metadata", str(tmp / "meta.tsv"), "--out-tsv", str(out), "--test", "test",
                                   "--reference", "ref", "--min-reads-per-sample", "5", "--chunk-rows", "50"],
                                  capture_output=True, text=True, env=run_env)
            if proc.returncode != 0:
                raise SystemExit(f"tail test failed:\n{proc.stderr}")
            outs.append(out.read_bytes())
        if outs[0] != outs[1]:
            raise AssertionError("test_condition_tail_diffs: output depends on PYTHONHASHSEED")
        if len(pd.read_csv(tmp / "tail_1.tsv", sep="\t")) < 50:
            raise AssertionError("test_condition_tail_diffs: expected >=50 tested fragmentforms")
        # (d) position-block streaming (report + gene browser): no (chrom, start0) block may be split
        #     across chunks, and per-block aggregation must equal the whole-table aggregation, in order
        from genotype_utils import iter_tsv_position_blocks
        whole_all = pd.read_csv(tmp / "long.tsv", sep="\t", low_memory=False,
                                dtype={c: str for c in ("sample", "chrom", "strand", "mod_code", "gene_name")})
        keys = ["gene_name", "chrom", "start0", "strand", "mod_code", "ZN_transcript_index"]
        g_whole = whole_all.groupby(keys, sort=False)[["Nvalid_cov", "Nmod"]].sum().reset_index()
        parts, n_rows, last_key = [], 0, None
        for blk in iter_tsv_position_blocks(str(tmp / "long.tsv"), list(whole_all.columns), chunksize=53,
                                            dtype={c: str for c in ("sample", "chrom", "strand", "mod_code", "gene_name")}):
            first_key = (blk["chrom"].iloc[0], int(blk["start0"].iloc[0]))
            if last_key is not None and first_key == last_key:
                raise AssertionError("iter_tsv_position_blocks split a (chrom, start0) block across chunks")
            last_key = (blk["chrom"].iloc[-1], int(blk["start0"].iloc[-1]))
            n_rows += len(blk)
            parts.append(blk.groupby(keys, sort=False)[["Nvalid_cov", "Nmod"]].sum().reset_index())
        if n_rows != len(whole_all):
            raise AssertionError(f"iter_tsv_position_blocks yielded {n_rows} rows, table has {len(whole_all)}")
        pd.testing.assert_frame_equal(pd.concat(parts, ignore_index=True), g_whole)
        # (e) the engine itself: n_workers must not change any statistic
        sites = []
        for i in range(300):
            n = np.array([rng.randint(5, 60) for _ in range(6)], dtype=float)
            k = np.array([sum(1 for _ in range(int(nn)) if rng.random() < (0.2 if j < 3 else 0.35))
                          for j, nn in enumerate(n)], dtype=float)
            sites.append((i, k, n, np.array([0, 0, 0, 1, 1, 1])))
        r1 = diffstats.beta_binomial_diff(sites, n_workers=1)
        r3 = diffstats.beta_binomial_diff(sites, n_workers=3)
        if [(r["key"], r["p_value"], r["lrt_stat"], r["theta_shrunk"]) for r in r1] != \
           [(r["key"], r["p_value"], r["lrt_stat"], r["theta_shrunk"]) for r in r3]:
            raise AssertionError("diffstats.beta_binomial_diff: parallel result differs from sequential")


if __name__ == "__main__":
    main()
