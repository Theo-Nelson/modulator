#!/usr/bin/env python3
"""
Guard the ZN-track collision case in classify_diff_sites.load_isoforms.

A ``zn_index`` is a display/aggregation TRACK, and metagene partitioning is allowed to place two
NON-OVERLAPPING fragmentforms on the SAME track (metagene_partition_count < number of transcripts).
Keying isoform exon models by (gene, zn_index) would then MERGE two disjoint transcripts into one
structurally-meaningless model, corrupting every classification for that gene. load_isoforms must keep
the track members separate (primary + variants) and _iso_at must resolve the (gene, zn, pos) triple to
the fragmentform that actually contains the base.

Usage: <modulator-env>/bin/python resources/synthetic_3exon/test_zn_collision.py
"""
from __future__ import annotations
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "workflow" / "scripts"))
import classify_diff_sites as C  # noqa: E402

CHROM = "chrUnit"


def write_gtf(path):
    # gene G: two fragmentforms on the SAME zn track (zn_index=1), disjoint spans; plus a long form.
    #   T1 (zn 1): exons 1000-1200, 1400-1600      (the low-coord region)
    #   T4 (zn 1): exons 5000-5200, 5400-5600      (the high-coord region)  <- shares zn 1 with T1
    #   T3 (zn 2): exons 1000-1200, ... , 5400-5600 (spans everything; the anchor)
    tx = [
        ("1", "1", [(1000, 1200), (1400, 1600)], 1600, 100),
        ("4", "1", [(5000, 5200), (5400, 5600)], 5600, 40),
        ("3", "2", [(1000, 1200), (3000, 3100), (5400, 5600)], 5600, 80),
    ]
    with open(path, "w") as fh:
        for tid, zn, exons, tes, rs in tx:
            gmin, gmax = min(s for s, _ in exons), max(e for _, e in exons)
            attrs = (f'gene_id "G"; ref_gene_name "G"; transcript_id "G.T{tid}"; '
                     f'transcript_index "{tid}"; zn_index "{zn}"; tes "{tes}"; read_support "{rs}";')
            fh.write(f"{CHROM}\tunit\ttranscript\t{gmin}\t{gmax}\t.\t+\t.\t{attrs}\n")
            for s, e in exons:
                fh.write(f"{CHROM}\tunit\texon\t{s}\t{e}\t.\t+\t.\t{attrs}\n")


def main():
    checks = []

    def check(name, ok):
        checks.append(ok)
        print(f"  {'PASS' if ok else '**FAIL**'}  {name}")

    with tempfile.TemporaryDirectory() as td:
        gtf = Path(td) / "u.gtf"
        write_gtf(gtf)
        iso, genes = C.load_isoforms(str(gtf), tes_tol=25, inside_tol=50)

        d = iso[("G", "1")]
        # the merged model would have 4 exons (T1's 2 + T4's 2); the correct primary has just T1's 2.
        check("zn-1 track NOT merged (primary has 2 exons, not 4)", len(d["exons"]) == 2)
        check("zn-1 track keeps both members as variants", len(d.get("variants", [])) == 2)

        # a base at 5500 lives in T4 (5400-5600), NOT in the primary T1 -> _iso_at must pick T4.
        m = C._iso_at(iso, "G", "1", 5500)
        check("_iso_at(5500) resolves to the T4 member (span 5000-5600)",
              m["exons"][0][0] == 5000 and m["exons"][-1][1] == 5600)
        # a base at 1500 lives in T1 (1400-1600) -> _iso_at must pick T1.
        m2 = C._iso_at(iso, "G", "1", 1500)
        check("_iso_at(1500) resolves to the T1 member (span 1000-1600)",
              m2["exons"][0][0] == 1000 and m2["exons"][-1][1] == 1600)
        # status of 5500 in the RESOLVED zn-1 model is exonic (would be 'absent' in the merged model's
        # gap, i.e. misclassified) -> proves the collision no longer corrupts the call.
        check("base 5500 is exonic in the resolved zn-1 model",
              C.status_in(m["exons"], "+", 5500) in ("exonic_terminal", "exonic_internal"))

    n_fail = checks.count(False)
    print(f"\nzn-collision: {len(checks) - n_fail}/{len(checks)} checks passed"
          + ("" if not n_fail else f"  ({n_fail} FAILED)"))
    sys.exit(1 if n_fail else 0)


if __name__ == "__main__":
    main()
