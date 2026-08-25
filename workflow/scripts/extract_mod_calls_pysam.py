#!/usr/bin/env python3
"""Streaming pysam replacement for `modkit extract calls`.

WARNING (BLOCKER-4): this standalone utility is NOT used by the pipeline (build_molecule_mod_table has
its own, fixed, pysam backend). It still has the implicit-MM bug: read.modified_bases returns only the
positions LISTED in the MM tag, so for an implicit ('.') MM group every unlisted canonical base -- a
real UNMODIFIED observation that modkit emits -- is dropped, inflating every modified fraction on real
ONT data. Do not use it on implicit-MM BAMs without porting the parse_mm_groups / implicit-canonical
logic from build_molecule_mod_table.extract_rows_pysam.

modkit's `extract` is either memory-bounded-but-serial-slow (`--ignore-index`) or
fast-but-OOM (parallel index scan accumulates unboundedly, esp. on acrocentric/rDNA
chromosomes). This reproduces the per-(read, candidate-site) call table that
build_molecule_mod_table consumes, but streams one read at a time -> O(1) memory, one
pass, never OOMs, on any chromosome or the whole genome. No reference FASTA needed
(modkit only uses it for the ref_kmer column, which is not consumed downstream).

Output columns are exactly the ones build_molecule_mod_table.parse_extracted_calls reads.
Semantics validated byte-identical vs `modkit extract calls --no-filtering --mapped-only`
(688,932 rows chr10, 0 diffs):
  call_prob = (ML + 0.5) / 256                       (as float32, modkit's formatting)
  canonical_prob = 1 - sum(mod_probs)
  call_code = argmax over {'-': canonical} U {mod_code: mod_prob}; call_prob = that max
  fail=false, inferred=false, within_alignment=true  (constant with these flags)
  ref_strand = '-' if read.is_reverse else '+'; --include-bed matching is strand-aware
"""
import argparse
import gzip
from collections import defaultdict

import numpy as np
import pysam

COLS = ["read_id", "chrom", "ref_position", "ref_strand", "call_prob", "call_code",
        "canonical_base", "modified_primary_base", "fail", "inferred", "within_alignment"]


def load_candidates(bed_path):
    """(chrom, start0) -> set of strands, from a BED6 candidate-site file."""
    d = defaultdict(set)
    op = gzip.open if bed_path.endswith(".gz") else open
    with op(bed_path, "rt") as fh:
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            if len(f) < 6:
                continue
            d[(f[0], int(f[1]))].add(f[5])
    return d


def extract(bam_path, cand, out_fh, region=None):
    bam = pysam.AlignmentFile(bam_path, "rb")
    itr = bam.fetch(region) if region else bam.fetch(until_eof=True)
    out_fh.write("\t".join(COLS) + "\n")
    n = 0
    f32 = np.float32
    for read in itr:
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        mb = read.modified_bases
        if not mb:
            continue
        chrom = read.reference_name
        ref_strand = "-" if read.is_reverse else "+"
        # per-query-position reference coordinate (0-based); None where not aligned to ref
        refpos = read.get_reference_positions(full_length=True)
        nrp = len(refpos)
        # gather {read_pos -> {mod_code: ml}} and the canonical base at that read_pos
        pos_mods = defaultdict(dict)
        pos_base = {}
        for (base, _strand, mod_code), calls in mb.items():
            code = str(mod_code)
            for read_pos, ml in calls:
                pos_mods[read_pos][code] = ml
                pos_base[read_pos] = base
        for read_pos, mods in pos_mods.items():
            if read_pos >= nrp:
                continue
            rp = refpos[read_pos]
            if rp is None:
                continue
            strands = cand.get((chrom, rp))
            if not strands or ref_strand not in strands:
                continue
            # argmax over canonical + mods
            mod_sum = 0.0
            best_code = None
            best_prob = -1.0
            for code, ml in mods.items():
                p = (ml + 0.5) / 256.0
                mod_sum += p
                if p > best_prob:
                    best_prob = p
                    best_code = code
            canon = 1.0 - mod_sum
            if canon >= best_prob:
                call_code = "-"
                call_prob = str(f32(canon))    # str() -> "0.9941406"; repr() wraps as np.float32(...)
            else:
                call_code = best_code
                call_prob = str(f32(best_prob))
            base = str(pos_base[read_pos])
            out_fh.write(f"{read.query_name}\t{chrom}\t{rp}\t{ref_strand}\t{call_prob}\t"
                         f"{call_code}\t{base}\t{base}\tfalse\tfalse\ttrue\n")
            n += 1
    return n


def main():
    ap = argparse.ArgumentParser(description="Streaming pysam replacement for `modkit extract calls`.")
    ap.add_argument("--bam", required=True, help="modBAM with MM/ML tags")
    ap.add_argument("--candidate-bed", required=True, help="BED6 candidate mod sites (strand-aware)")
    ap.add_argument("--out", required=True, help="Output TSV (.gz/.bgz -> gzip-compressed)")
    ap.add_argument("--region", default=None, help="Restrict to a region/chrom (default: whole BAM)")
    args = ap.parse_args()
    cand = load_candidates(args.candidate_bed)
    op = gzip.open if args.out.endswith((".gz", ".bgz")) else open
    with op(args.out, "wt") as out_fh:
        n = extract(args.bam, cand, out_fh, region=args.region)
    import sys
    print(f"[pysam-extract] wrote {n} rows -> {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
