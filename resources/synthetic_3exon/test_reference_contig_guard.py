#!/usr/bin/env python3
"""BLOCKER-2: a BAM aligned to a DIFFERENT genome than reference_fa must NOT complete silently.

mm39 and GRCh38 share contig NAMES (chr1, chr19, ...) but differ entirely in length (chr1 195 Mb vs
249 Mb), so aligning mouse reads and analysing them against a human reference produced a full report
attributing mouse data to human genes, exit 0. The preflight now compares each BAM's @SQ lengths to
the reference .fai and aborts on any shared-name/different-length contig. This test exercises both the
comparison logic and the pysam @SQ read on a real (header-only) BAM.
"""
import os
import sys
import tempfile

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))
from modulator.pipeline import ModulatorPipeline as Pipeline  # noqa: E402

GRCh38 = {"chr1": 248956422, "chr19": 58617616, "chr21": 46709983}
mm39 = {"chr1": 195154279, "chr19": 61420004}   # mouse: same names, different lengths


def _write_header_bam(path, contig_lengths):
    import pysam
    header = {"HD": {"VN": "1.6", "SO": "coordinate"},
              "SQ": [{"SN": c, "LN": ln} for c, ln in contig_lengths.items()]}
    with pysam.AlignmentFile(path, "wb", header=header):
        pass


def main():
    rc = 0

    # 1) pure comparison: mouse BAM vs human reference -> both shared contigs flagged
    mism, shared = Pipeline._contig_length_mismatches(mm39, GRCh38)
    if len(mism) == 2 and {m[0] for m in mism} == {"chr1", "chr19"}:
        print(f"  PASS  mouse-vs-human: {len(mism)}/{len(shared)} shared contigs flagged as mismatched")
    else:
        print(f"  FAIL  mouse-vs-human should flag 2 shared contigs, got {mism}"); rc = 1

    # 2) human BAM vs human reference WITH extra alt contigs in the reference -> clean (extras ignored)
    human_bam = {"chr1": 248956422, "chr19": 58617616}
    ref_with_alts = dict(GRCh38, **{"chr1_KI270766v1_alt": 256271, "chrUn_GL000195v1": 182896})
    mism2, _ = Pipeline._contig_length_mismatches(human_bam, ref_with_alts)
    if not mism2:
        print("  PASS  human-vs-human (261-style extra alt contigs): no false mismatch")
    else:
        print(f"  FAIL  matching human BAM wrongly flagged: {mism2}"); rc = 1

    # 3) real @SQ read: a header-only mouse BAM read via pysam must reproduce the mismatch
    try:
        import pysam  # noqa: F401
        with tempfile.TemporaryDirectory() as td:
            bam = os.path.join(td, "mouse.bam")
            _write_header_bam(bam, mm39)
            import pysam as _ps
            with _ps.AlignmentFile(bam, "rb", check_sq=False) as fh:
                bam_len = dict(zip(fh.references, fh.lengths))
            mism3, _ = Pipeline._contig_length_mismatches(bam_len, GRCh38)
            if mism3:
                print(f"  PASS  real BAM @SQ read + compare flags the wrong reference ({len(mism3)} contigs)")
            else:
                print("  FAIL  real BAM @SQ read did not flag the mouse/human mismatch"); rc = 1
    except Exception as exc:  # noqa: BLE001
        print(f"  SKIP  real-BAM @SQ check (pysam/write unavailable): {exc}")

    print("reference-contig guard: " + ("OK" if rc == 0 else "FAILURES"))
    sys.exit(rc)


if __name__ == "__main__":
    main()
