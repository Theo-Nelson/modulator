#!/usr/bin/env python3
"""mod x mod co-occurrence must EXCLUDE same-base pairs.

Two modification codes assayed at the SAME genomic base (e.g. A->inosine 17596 and A->m6A 'a') are
mutually exclusive BY CONSTRUCTION -- a single molecule's base carries at most one -- so pairing them
yields an OR~0 artifact that tops any mutually-exclusive ranking. The test must skip pairs at the same
(chrom, start0, strand) while still testing genuine different-base pairs.

Fixture: base 100 carries codes 'a' and '17596'; base 140 carries 'a'. The only same-base pair is
(100:a, 100:17596); the (100,140) and (100:17596,140) pairs are legitimate different-base pairs.
"""
import subprocess, sys, tempfile, csv
from pathlib import Path

HERE = Path(__file__).resolve().parent
SCRIPT = HERE.parent.parent / "workflow" / "scripts" / "test_mod_mod_assoc.py"

SITES = [(100, "a"), (100, "17596"), (140, "a")]   # (start0, mod_code)


def build_mods(path):
    fields = ["sample", "qname", "mod_site_id", "chrom", "start0", "strand", "target_mod_code",
              "state_detail", "gene_name", "metagene_index", "usable", "fail", "within_alignment"]
    rows = []
    for i in range(60):
        qn = f"r{i}"
        # deterministic, decorrelated-ish states so every 2x2 cell is populated for each pair
        st = {(100, "a"): (i % 2 == 0),
              (100, "17596"): (i % 5 == 0),
              (140, "a"): (i % 3 != 0)}
        for (pos, code), modified in st.items():
            rows.append({
                "sample": "s1", "qname": qn, "mod_site_id": f"chr1:{pos}:{code}",
                "chrom": "chr1", "start0": pos, "strand": "+", "target_mod_code": code,
                "state_detail": "modified" if modified else "canonical",
                "gene_name": "GENE", "metagene_index": 1,
                "usable": "True", "fail": "False", "within_alignment": "True",
            })
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(rows)


def main():
    rc = 0
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        mods, out = td / "mods.tsv", td / "out.tsv"
        build_mods(mods)
        subprocess.run([sys.executable, str(SCRIPT), "--molecule-mods", str(mods),
                        "--out-tsv", str(out), "--max-distance", "1000", "--min-state-reads", "4"],
                       check=True)
        with open(out) as fh:
            pairs = list(csv.DictReader(fh, delimiter="\t"))

    same_base = [p for p in pairs if p["start0_a"] == p["start0_b"]]
    diff_base = [p for p in pairs if p["start0_a"] != p["start0_b"]]
    if same_base:
        print(f"  FAIL  {len(same_base)} same-base pair(s) present (mutually exclusive by construction, "
              f"should be excluded): {[(p['mod_code_a'], p['mod_code_b'], p['start0_a']) for p in same_base]}")
        rc = 1
    else:
        print("  PASS  same-base pair (two mod codes on one base) excluded")
    if diff_base:
        print(f"  PASS  genuine different-base pairs still tested ({len(diff_base)} pair(s))")
    else:
        print("  FAIL  different-base pairs missing -- exclusion was too broad")
        rc = 1

    print("mod_mod same-base exclusion: " + ("OK" if rc == 0 else "FAILURES"))
    sys.exit(rc)


if __name__ == "__main__":
    main()
