#!/usr/bin/env python3
"""M4: the between-condition usage test must not truncate dotted gene names or merge distinct genes.

zt_label is `{gene}.{gene_id}.G<n>.T<n>` and BOTH parts can contain dots (GENCODE clone names like
CTC-338M12.4, versioned Ensembl ids). The old `_gene_of = zt.split(".")[0]` truncated the gene and
collapsed CTC-338M12.4 and CTC-338M12.3 into one "CTC-338M12", dropping features and rescoring the
survivors against a merged denominator. The gene now comes from the authoritative gtf_gene_name in
the classification summary, with a no-merge fallback (`{gene}.{gene_id}`) when it is absent.
"""
import csv
import os
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
SCRIPT = HERE.parent.parent / "workflow" / "scripts" / "test_condition_usage_diffs.py"

# two DISTINCT dotted-name genes that split(".")[0] merges into "CTC-338M12", plus a normal gene
ZT_GENE = {
    "CTC-338M12.4.ENSG1.G1.T1": "CTC-338M12.4", "CTC-338M12.4.ENSG1.G1.T2": "CTC-338M12.4",
    "CTC-338M12.3.ENSG2.G2.T1": "CTC-338M12.3", "CTC-338M12.3.ENSG2.G2.T2": "CTC-338M12.3",
    "EEF2.ENSG3.G3.T1": "EEF2", "EEF2.ENSG3.G3.T2": "EEF2",
}
COUNTS = {
    "CTC-338M12.4.ENSG1.G1.T1": [80, 82, 20, 18], "CTC-338M12.4.ENSG1.G1.T2": [20, 18, 80, 82],
    "CTC-338M12.3.ENSG2.G2.T1": [50, 52, 50, 48], "CTC-338M12.3.ENSG2.G2.T2": [50, 48, 50, 52],
    "EEF2.ENSG3.G3.T1": [70, 72, 30, 28], "EEF2.ENSG3.G3.T2": [30, 28, 70, 72],
}
SAMPLES = ["m1", "m2", "z1", "z2"]


def _write(td):
    with open(td / "tx.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t"); w.writerow(["zt_label"] + SAMPLES)
        for z in ZT_GENE:
            w.writerow([z] + COUNTS[z])
    with open(td / "cls.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t"); w.writerow(["zt_label", "gtf_gene_name", "chrom", "strand", "iso_tes"])
        for z, g in ZT_GENE.items():
            w.writerow([z, g, "chr1", "+", 1000])
    with open(td / "meta.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t"); w.writerow(["sample", "condition"])
        for s in SAMPLES:
            w.writerow([s, "mock" if s.startswith("m") else "zikv"])


def _run(td, with_summary):
    out = td / ("out_%s.tsv" % ("cls" if with_summary else "nocls"))
    cmd = [sys.executable, str(SCRIPT), "--tx-counts", str(td / "tx.tsv"),
           "--sample-metadata", str(td / "meta.tsv"), "--out-tsv", str(out),
           "--feature", "isoform", "--test", "zikv", "--reference", "mock",
           "--min-gene-reads", "20", "--min-samples-per-group", "2"]
    if with_summary:
        cmd += ["--classification-summary", str(td / "cls.tsv")]
    subprocess.run(cmd, check=True)
    with open(out) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    return {r["gene_name"] for r in rows}


def main():
    rc = 0
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        _write(td)

        genes = _run(td, with_summary=True)
        if {"CTC-338M12.4", "CTC-338M12.3"} <= genes:
            print("  PASS  authoritative gtf_gene_name keeps CTC-338M12.4 and .3 as DISTINCT genes")
        elif "CTC-338M12" in genes:
            print(f"  FAIL  dotted gene names merged into 'CTC-338M12' (M4 not fixed): {sorted(genes)}"); rc = 1
        else:
            print(f"  FAIL  unexpected gene set: {sorted(genes)}"); rc = 1

        # fallback (no summary): must NOT merge distinct genes either (keeps {gene}.{gene_id})
        genes_fb = _run(td, with_summary=False)
        merged = any(g == "CTC-338M12" for g in genes_fb)
        n_ctc = len({g for g in genes_fb if g.startswith("CTC-338M12")})
        if not merged and n_ctc == 2:
            print("  PASS  fallback (no summary) keeps the two CTC-338M12.* genes separate (no merge)")
        else:
            print(f"  FAIL  fallback merged distinct genes: {sorted(genes_fb)}"); rc = 1

    print("usage dotted-gene: " + ("OK" if rc == 0 else "FAILURES"))
    sys.exit(rc)


if __name__ == "__main__":
    main()
