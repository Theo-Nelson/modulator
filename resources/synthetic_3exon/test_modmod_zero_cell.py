#!/usr/bin/env python3
"""
Guard the single-molecule co-localized-modification (mod x mod) test on its most
important cases: a 2x2 co-occurrence table with a ZERO in a discordant cell.

The 2x2 is [[both_modified, a_only], [b_only, neither]]. The two discordant
cells -- a_only (A modified, B not) and b_only (B modified, A not) -- going to
zero means one modification essentially never occurs without the other: a strong,
asymmetric mechanistic dependency. Those pairs MUST (a) not crash the test and
(b) rank at the very top of statistical significance. This test hand-builds a
per-read modification-call table so a pair lands with a_only == 0, runs the REAL
test_mod_mod_assoc.py, and asserts the pair is emitted, called CONCORDANT, highly
significant, and ranked above an independent control pair.

Usage: <modulator-env>/bin/python resources/synthetic_3exon/test_modmod_zero_cell.py
"""
from __future__ import annotations
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd

SCRIPT = Path(__file__).resolve().parents[2] / "workflow" / "scripts" / "test_mod_mod_assoc.py"
CHROM = "chrUnit"


def _rows(qbase, site_id, start0, states):
    """One mod-call record per read for a site. states: list of 1(mod)/0(canonical)."""
    out = []
    for i, s in enumerate(states):
        out.append({
            "sample": "S1", "qname": f"{qbase}{i:04d}", "mod_site_id": site_id,
            "chrom": CHROM, "start0": start0, "strand": "+", "target_mod_code": "a",
            "state_detail": "modified" if s else "canonical",
            "gene_name": "GENEZ", "gene_names": "GENEZ", "metagene_index": 7,
            "usable": True, "fail": False, "within_alignment": True,
        })
    return out


def build_mods(path: Path):
    rows = []
    # --- Pair 1 (ZERO discordant cell): a_only == 0 -----------------------------
    # 100 shared reads: 40 both-modified, 0 A-only, 15 B-only, 45 neither.
    #   siteA modified on reads 0..39   (40 modified, 60 canonical)
    #   siteB modified on reads 0..54   (55 modified, 45 canonical)
    # => A modified never happens without B modified -> a_only = 0.
    a_states = [1] * 40 + [0] * 60
    b_states = [1] * 55 + [0] * 45
    rows += _rows("dep_", f"{CHROM}:100:+:a", 100, a_states)
    rows += _rows("dep_", f"{CHROM}:150:+:a", 150, b_states)
    # --- Pair 2 (INDEPENDENT control): balanced 2x2, no dependency ---------------
    c_states = ([1, 0] * 25) + ([1, 0] * 25)          # 50 modified interleaved
    d_states = ([1, 1, 0, 0] * 25)                    # 50 modified, independent of C
    rows += _rows("ind_", f"{CHROM}:400:+:a", 400, c_states)
    rows += _rows("ind_", f"{CHROM}:450:+:a", 450, d_states)
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)


def main():
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        mods, out = td / "mods.tsv", td / "modmod.tsv"
        build_mods(mods)
        subprocess.run([sys.executable, str(SCRIPT), "--molecule-mods", str(mods),
                        "--out-tsv", str(out), "--min-pair-reads", "8", "--min-state-reads", "4"],
                       check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        df = pd.read_csv(out, sep="\t")

    checks = []

    def check(name, ok):
        checks.append((name, ok))
        print(f"  {'PASS' if ok else '**FAIL**'}  {name}")

    # the zero-cell dependent pair is the one spanning start0 100 <-> 150
    dep = df[(df["start0_a"] == 100) & (df["start0_b"] == 150)]
    ind = df[(df["start0_a"] == 400) & (df["start0_b"] == 450)]
    check("dependent pair (a_only=0) is emitted (no crash)", len(dep) == 1)
    if len(dep) == 1:
        d = dep.iloc[0]
        check("a_only cell is exactly 0", int(d["n_a_only"]) == 0)
        check("dependent pair called CONCORDANT", d["direction"] == "CONCORDANT")
        check("dependent pair is significant (p_adj < 1e-6)", float(d["p_adj_bh"]) < 1e-6)
    if len(dep) == 1 and len(ind) == 1:
        check("dependent (zero-cell) pair ranks above the independent control",
              float(dep.iloc[0]["p_adj_bh"]) < float(ind.iloc[0]["p_adj_bh"]))
        check("dependent pair is the top-ranked row",
              df.sort_values("p_adj_bh").iloc[0][["start0_a", "start0_b"]].tolist() == [100, 150])

    n_fail = sum(1 for _, ok in checks if not ok)
    print(f"\nmod-mod zero-cell: {len(checks) - n_fail}/{len(checks)} checks passed"
          + ("" if not n_fail else f"  ({n_fail} FAILED)"))
    sys.exit(1 if n_fail else 0)


if __name__ == "__main__":
    main()
