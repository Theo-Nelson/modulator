#!/usr/bin/env python3
"""BLOCKER-4 (continued): the gene-level between-condition tail test must not turn a pure
isoform-USAGE shift into a false tail difference.

The gene summary is the mean of per-fragmentform median tails. A fragmentform only enters that
mean where it clears --min-reads-per-sample. If a usage shift pushes a form below the threshold in
one condition, that form silently leaves the mean there, and the gene mean jumps -- a STEP at the
threshold, not a tail change. The fix restricts the mean to fragmentforms that clear the threshold
in EVERY contributing sample, so both conditions average the same form set.

Two genes, both with a long form (~130 nt) and a short form (~70/others) whose OWN tails are flat:
  GX  usage-shift-only  : short form drops below threshold in the test condition -> must read ~0 nt
  GY  genuine gene shift : BOTH forms move +30 nt, staying above threshold -> must survive
"""
import subprocess, sys, tempfile, csv
from pathlib import Path

HERE = Path(__file__).resolve().parent
SCRIPT = HERE.parent.parent / "workflow" / "scripts" / "test_condition_tail_diffs.py"

REF_SAMPLES = ["m1", "m2", "m3"]
TEST_SAMPLES = ["z1", "z2", "z3"]


def _reads(rows, gene, zt, sample, median_nt, n):
    # n identical reads -> the per-(gene,zt,sample) median is exactly median_nt.
    for _ in range(n):
        rows.append({"sample": sample, "tail_len": median_nt, "ZT": zt, "gene_name": gene})


def build_tail_tsv(path):
    rows = []
    # --- GX: short form T2 present in mock (20 reads) but sub-threshold in zikv (3 reads) ---
    for i, s in enumerate(REF_SAMPLES):
        _reads(rows, "GX", "GX.G1.T1", s, 130 + (i - 1), 20)   # long, flat, qualifies
        _reads(rows, "GX", "GX.G1.T2", s, 70 + (i - 1), 20)    # short, flat, qualifies in mock
    for i, s in enumerate(TEST_SAMPLES):
        _reads(rows, "GX", "GX.G1.T1", s, 130 + (i - 1), 20)   # long, flat, qualifies
        _reads(rows, "GX", "GX.G1.T2", s, 70 + (i - 1), 3)     # short, flat, SUB-THRESHOLD
    # --- GY: genuine gene-wide shift, both forms +30 nt, both above threshold everywhere ---
    for i, s in enumerate(REF_SAMPLES):
        _reads(rows, "GY", "GY.G1.T1", s, 100 + (i - 1), 20)
        _reads(rows, "GY", "GY.G1.T2", s, 60 + (i - 1), 20)
    for i, s in enumerate(TEST_SAMPLES):
        _reads(rows, "GY", "GY.G1.T1", s, 130 + (i - 1), 20)
        _reads(rows, "GY", "GY.G1.T2", s, 90 + (i - 1), 20)
    # --- GZ: coverage in every replicate, but NO form clears the threshold in all of them (T1 dips
    # in one zikv rep, T2 in one mock rep) -> no common form -> must appear as an UNTESTABLE row,
    # not silently vanish. reads: qualifying=20, sub-threshold=3.
    gz = {"m1": (20, 20), "m2": (20, 20), "m3": (20, 3),
          "z1": (20, 20), "z2": (3, 20), "z3": (20, 20)}
    for s, (n1, n2) in gz.items():
        _reads(rows, "GZ", "GZ.G1.T1", s, 110, n1)
        _reads(rows, "GZ", "GZ.G1.T2", s, 80, n2)
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["sample", "tail_len", "ZT", "gene_name"], delimiter="\t")
        w.writeheader()
        w.writerows(rows)


def build_meta(path):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["sample", "condition"])
        for s in REF_SAMPLES:
            w.writerow([s, "mock"])
        for s in TEST_SAMPLES:
            w.writerow([s, "zikv"])


def main():
    rc = 0
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        tail_tsv, meta_tsv, out_tsv = td / "tails.tsv", td / "meta.tsv", td / "out.tsv"
        build_tail_tsv(tail_tsv)
        build_meta(meta_tsv)
        subprocess.run([sys.executable, str(SCRIPT),
                        "--tail-tsv", str(tail_tsv), "--sample-metadata", str(meta_tsv),
                        "--out-tsv", str(out_tsv), "--level", "gene",
                        "--test", "zikv", "--reference", "mock",
                        "--min-reads-per-sample", "10", "--min-samples-per-group", "2",
                        "--min-tail", "1"], check=True)
        with open(out_tsv) as fh:
            out = {r["gene_name"]: r for r in csv.DictReader(fh, delimiter="\t")}

    gx = out.get("GX")
    if gx is None:
        print("  FAIL: GX produced no row"); rc = 1
    else:
        d = abs(float(gx["delta_nt"]))
        if d < 5.0:
            print(f"  PASS  GX usage-shift-only reads ~0 nt (no threshold-step false positive): delta={gx['delta_nt']} nt")
        else:
            print(f"  FAIL  GX usage-shift-only reports {gx['delta_nt']} nt (threshold-step confound NOT removed)"); rc = 1

    gy = out.get("GY")
    if gy is None:
        print("  FAIL: GY produced no row"); rc = 1
    else:
        d = float(gy["delta_nt"])
        if d > 20.0:
            print(f"  PASS  GY genuine gene-wide shift survives: delta={gy['delta_nt']} nt")
        else:
            print(f"  FAIL  GY real +30 nt effect lost after the fix: delta={gy['delta_nt']} nt"); rc = 1

    # GZ: no common fragmentform -> must be VISIBLE as an untestable row (NaN stats + status),
    # never silently dropped (the silent-data-loss MAJOR).
    gz = out.get("GZ")
    if gz is None:
        print("  FAIL  GZ (no common fragmentform) silently DROPPED -- should be an untestable row"); rc = 1
    elif gz["p_value"].strip() != "" or "untestable_no_common_fragmentform" not in gz["per_replicate_json"]:
        print(f"  FAIL  GZ present but not marked untestable: p_value={gz['p_value']!r} json={gz['per_replicate_json']!r}"); rc = 1
    else:
        print("  PASS  GZ (no common fragmentform) emitted as a visible untestable row, not dropped")

    print("condition-tail threshold-step: " + ("OK" if rc == 0 else "FAILURES"))
    sys.exit(rc)


if __name__ == "__main__":
    main()
