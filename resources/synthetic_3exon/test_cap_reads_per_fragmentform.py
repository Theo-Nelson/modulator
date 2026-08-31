#!/usr/bin/env python3
"""Round-trip test for cap_reads_per_fragmentform.py -- the per-fragmentform depth cap that is the
memory safety valve for deep loci. It ships OFF by default (max_reads_per_fragmentform=0), so absent a
test it would first execute on a user's hardest dataset. This builds a tiny BAM with known ZT tags and
asserts the cap's contract: (1) each fragmentform keeps <= N distinct qnames; (2) a fragmentform UNDER
the cap is left fully intact; (3) untagged reads are all kept; (4) all alignments of a kept qname are
kept (qname-level, not alignment-level -- so a read's secondaries survive with it); (5) the kept set is
deterministic for a fixed seed and responds to the seed; (6) the output is indexed (.bai written)."""
import os
import subprocess
import sys
import tempfile

import pysam

HERE = os.path.dirname(os.path.abspath(__file__))
CAP = os.path.join(HERE, "..", "..", "workflow", "scripts", "cap_reads_per_fragmentform.py")

T1 = "GENE.ENSG1.G1.T1"   # 10 qnames -> capped to 5
T2 = "GENE.ENSG1.G1.T2"   # 3 qnames -> under cap, untouched
N_CAP = 5


def _seg(header, qname, pos, tag, secondary=False):
    a = pysam.AlignedSegment(header)
    a.query_name = qname
    a.flag = 0x100 if secondary else 0
    a.reference_id = 0
    a.reference_start = pos
    a.mapping_quality = 60
    a.query_sequence = "ACGTACGTAC"
    a.cigarstring = "10M"
    a.query_qualities = pysam.qualitystring_to_array("IIIIIIIIII")
    if tag is not None:
        a.set_tag("ZT", tag, value_type="Z")
    return a


def build_bam(path):
    header = pysam.AlignmentHeader.from_dict(
        {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 100000}]})
    recs = []
    pos = 100
    for i in range(10):                       # T1: 10 distinct qnames
        recs.append((pos, f"t1_{i:02d}", T1, False)); pos += 5
    for i in range(3):                        # T2: 3 distinct qnames (under cap)
        recs.append((pos, f"t2_{i:02d}", T2, False)); pos += 5
    recs.append((pos, "t2_00", T2, True)); pos += 5   # a SECONDARY alignment of a T2 qname (always kept)
    for i in range(4):                        # untagged reads
        recs.append((pos, f"u_{i:02d}", None, False)); pos += 5
    recs.sort(key=lambda r: r[0])             # coordinate-sorted so the capped output is indexable
    with pysam.AlignmentFile(path, "wb", header=header) as out:
        for p, q, tag, sec in recs:
            out.write(_seg(header, q, p, tag, sec))


def run_cap(in_bam, out_bam, seed, n=N_CAP):
    subprocess.run([sys.executable, CAP, "--in-bam", in_bam, "--out-bam", out_bam,
                    "--max-per-tag", str(n), "--seed", str(seed)],
                   check=True, capture_output=True)


def qnames_by_tag(bam):
    d = {}
    n_untagged = set()
    align_counts = {}
    with pysam.AlignmentFile(bam, "rb") as f:
        for a in f.fetch(until_eof=True):
            align_counts[a.query_name] = align_counts.get(a.query_name, 0) + 1
            try:
                t = str(a.get_tag("ZT"))
            except KeyError:
                n_untagged.add(a.query_name); continue
            d.setdefault(t, set()).add(a.query_name)
    return d, n_untagged, align_counts


def main():
    checks = []
    with tempfile.TemporaryDirectory() as td:
        in_bam = os.path.join(td, "in.bam")
        build_bam(in_bam)
        o1 = os.path.join(td, "o1.bam")
        o1b = os.path.join(td, "o1b.bam")
        o2 = os.path.join(td, "o2.bam")
        run_cap(in_bam, o1, seed=1)
        run_cap(in_bam, o1b, seed=1)   # same seed -> identical
        run_cap(in_bam, o2, seed=7)    # different seed

        d1, unt1, ac1 = qnames_by_tag(o1)
        d1b, _, _ = qnames_by_tag(o1b)
        d2, _, _ = qnames_by_tag(o2)

        checks.append((f"T1 capped to exactly {N_CAP} distinct qnames (from 10)", len(d1.get(T1, set())) == N_CAP))
        checks.append(("T1 kept qnames are a subset of the originals", d1.get(T1, set()) <= {f"t1_{i:02d}" for i in range(10)}))
        checks.append(("T2 (under cap) left fully intact: 3 qnames", len(d1.get(T2, set())) == 3))
        checks.append(("all 4 untagged reads kept", unt1 == {f"u_{i:02d}" for i in range(4)}))
        checks.append(("secondary alignment of a kept qname survives (qname-level): t2_00 has 2 records", ac1.get("t2_00", 0) == 2))
        checks.append(("deterministic: same seed -> identical kept set", d1.get(T1) == d1b.get(T1)))
        checks.append(("seed matters: a different seed selects a different subset", d1.get(T1) != d2.get(T1)))
        checks.append(("still exactly N under the other seed", len(d2.get(T1, set())) == N_CAP))
        checks.append(("output BAM is indexed (.bai written)", os.path.exists(o1 + ".bai")))

    n_pass = 0
    for name, ok in checks:
        print(f"  [{'PASS' if ok else '**FAIL**'}] {name}")
        n_pass += bool(ok)
    print(f"\ncap_reads_per_fragmentform round-trip: {n_pass}/{len(checks)} checks passed")
    sys.exit(0 if n_pass == len(checks) else 1)


if __name__ == "__main__":
    main()
