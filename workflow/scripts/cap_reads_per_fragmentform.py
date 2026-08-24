#!/usr/bin/env python3
"""Depth-cap a BAM to at most N reads per fragmentform (ZT tag), by seeded random subsampling.

Modulator's per-molecule genotype analysis is per-fragmentform (each read's ZT tag = its assembled
transcript isoform). Highly-expressed loci pile up hundreds of thousands of reads on a few isoforms,
which makes the per-window `modkit extract` in build_molecule_mod_table saturate memory. Capping each
fragmentform to <=N reads collapses those deep isoforms to a bounded size while GUARANTEEING every
isoform keeps up to N reads (a flat per-region cap could let one dominant isoform starve a minor one).

Deterministic: with a fixed --seed the kept read set is identical across runs. Reads with no ZT tag
are kept as-is (they are sparse after multigene filtering and are not the deep-locus driver).
"""

import argparse
import random
import sys
from collections import defaultdict

import pysam


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--in-bam", required=True)
    ap.add_argument("--out-bam", required=True)
    ap.add_argument("--tag", default="ZT", help="Fragmentform tag to group on (default ZT).")
    ap.add_argument("--max-per-tag", type=int, default=400, help="Max reads (distinct qnames) per fragmentform.")
    ap.add_argument("--seed", type=int, default=12345)
    ap.add_argument("--threads", type=int, default=2)
    ap.add_argument("--verbose", action="store_true")
    return ap.parse_args()


def main():
    args = parse_args()
    tag = args.tag
    n = max(1, int(args.max_per_tag))
    io_threads = max(1, int(args.threads))

    # Pass 1: collect the distinct qnames per fragmentform tag.
    qname_tag = {}                      # qname -> its (first-seen) tag; membership = "this qname is tagged"
    tag_qnames = defaultdict(set)       # tag -> {qnames}
    with pysam.AlignmentFile(args.in_bam, "rb", threads=io_threads) as inp:
        for aln in inp.fetch(until_eof=True):
            try:
                t = str(aln.get_tag(tag))
            except KeyError:
                continue
            q = aln.query_name
            if q not in qname_tag:
                qname_tag[q] = t
            tag_qnames[t].add(q)

    # Seeded subsample: keep <=N qnames per tag. sorted() everywhere so the result depends only on
    # the seed (not on BAM read order or set iteration order).
    rng = random.Random(int(args.seed))
    kept = set()
    n_capped = 0
    for t in sorted(tag_qnames):
        qs = sorted(tag_qnames[t])
        if len(qs) > n:
            qs = rng.sample(qs, n)
            n_capped += 1
        kept.update(qs)

    # Pass 2: write kept-qname reads for tagged qnames; keep all untagged reads.
    n_in = 0
    n_out = 0
    with pysam.AlignmentFile(args.in_bam, "rb", threads=io_threads) as inp, \
         pysam.AlignmentFile(args.out_bam, "wb", header=inp.header, threads=io_threads) as outp:
        for aln in inp.fetch(until_eof=True):
            n_in += 1
            q = aln.query_name
            if q in qname_tag:
                if q in kept:
                    outp.write(aln)
                    n_out += 1
            else:
                outp.write(aln)
                n_out += 1
    # Do NOT swallow an index failure: every downstream reader fetch()es the capped BAM by region and
    # would fail (or silently read nothing) without a .bai. A cap run that could not index its output
    # is not a usable result -- fail loudly here instead of hours later.
    pysam.index(args.out_bam)

    if args.verbose:
        print(f"[cap] {args.in_bam}: fragmentforms={len(tag_qnames)} capped={n_capped} "
              f"(>{n} reads) seed={args.seed} reads {n_in}->{n_out}", file=sys.stderr, flush=True)


if __name__ == "__main__":
    main()
