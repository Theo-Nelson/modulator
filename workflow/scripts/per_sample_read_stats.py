#!/usr/bin/env python3
"""
per_sample_read_stats.py

Compute per-sample read count + read-length distribution stats for:
  1) ALL reads in the original BAM (total)
  2) Reads that PASS the same "considered" filters used by assemble_transcripts.py
  3) Reads ASSIGNED to an assembled transcript (have ZT tag) in the zt_tagged BAM
  4) (Optional but useful) Reads that appear in zt_tagged BAM but are UNASSIGNED (no ZT)

Also reports a breakdown of how many reads FAIL the "considered" filters by first-failure reason.

Notes:
- "considered" is evaluated on the *original* BAM, matching the assembler's filter order:
    unmapped -> secondary/supp -> low mapq -> low introns -> low softclip3p
- "assigned" is evaluated from zt_tagged BAM by presence of ZT tag.
  This assumes you generated zt_tagged bams via assemble_transcripts.py.
"""

import argparse
import array
import os
import glob
import sys
from concurrent.futures import ProcessPoolExecutor
from typing import Dict, List, Tuple

import numpy as np
import pysam


def new_len_buf() -> "array.array":
    """Read lengths accumulate into a C int32 array rather than a Python list: ~4 bytes/read
    instead of ~36 (pointer + int object). At 33M reads/sample that is ~132 MB vs ~1.2 GB, which
    is what makes it safe to process several samples in parallel. Values are unchanged, so every
    downstream statistic is bit-identical to the previous list-based implementation."""
    return array.array("i")


# ----------------------------- stats helpers -----------------------------

def qstats(lengths) -> Dict[str, int]:
    """Return {n, min, p25, p50, p75, p90, max, mean} for a sequence of lengths."""
    if len(lengths) == 0:
        return dict(n=0, min=0, p25=0, p50=0, p75=0, p90=0, max=0, mean=0.0)
    arr = np.asarray(lengths, dtype=np.int64)
    # numpy quantile method kw differs across versions; use a safe fallback
    try:
        p25 = int(np.quantile(arr, 0.25, method="linear"))
        p50 = int(np.quantile(arr, 0.50, method="linear"))
        p75 = int(np.quantile(arr, 0.75, method="linear"))
        p90 = int(np.quantile(arr, 0.90, method="linear"))
    except TypeError:
        # older numpy
        p25 = int(np.quantile(arr, 0.25, interpolation="linear"))
        p50 = int(np.quantile(arr, 0.50, interpolation="linear"))
        p75 = int(np.quantile(arr, 0.75, interpolation="linear"))
        p90 = int(np.quantile(arr, 0.90, interpolation="linear"))

    return dict(
        n=int(arr.size),
        min=int(arr.min()),
        p25=p25,
        p50=p50,
        p75=p75,
        p90=p90,
        max=int(arr.max()),
        mean=round(float(arr.mean()), 2),
    )


def intron_count(aln: pysam.AlignedSegment) -> int:
    """Count N ops in CIGAR (introns)."""
    ct = aln.cigartuples or []
    return sum(1 for op, ln in ct if op == 3)


def softclip3p_len(aln: pysam.AlignedSegment) -> int:
    """Return 3' softclip length in transcript space, using genomic strand from FLAG."""
    ct = aln.cigartuples or []
    if not ct:
        return 0
    tx = "-" if aln.is_reverse else "+"
    if tx == "+":
        return ct[-1][1] if ct and ct[-1][0] == 4 else 0
    else:
        return ct[0][1] if ct and ct[0][0] == 4 else 0


def get_len(aln: pysam.AlignedSegment) -> int:
    """Best-effort query length."""
    L = aln.query_length
    if L is not None:
        return int(L)
    qs = aln.query_sequence
    return int(len(qs)) if qs else 0


# ------------------------ considered filter breakdown ------------------------

FAIL_UNMAPPED = "FAIL_unmapped"
FAIL_SECONDARY_SUPP = "FAIL_secondary_or_supplementary"
FAIL_LOW_MAPQ = "FAIL_low_mapq"
FAIL_LOW_INTRONS = "FAIL_low_introns"
FAIL_LOW_SOFTCLIP3P = "FAIL_low_softclip3p"
PASS = "PASS"


def considered_fail_reason(
    aln: pysam.AlignedSegment,
    primary_only: bool,
    min_mapq: int,
    min_introns_read: int,
    require_softclip3p: int,
) -> str:
    """
    Return PASS if the read would be "considered" by the assembler,
    else return the FIRST failure reason (matching assembler ordering).
    """
    if aln.is_unmapped:
        return FAIL_UNMAPPED
    if primary_only and (aln.is_secondary or aln.is_supplementary):
        return FAIL_SECONDARY_SUPP
    if aln.mapping_quality < min_mapq:
        return FAIL_LOW_MAPQ
    if intron_count(aln) < min_introns_read:
        return FAIL_LOW_INTRONS
    if require_softclip3p > 0 and softclip3p_len(aln) < require_softclip3p:
        return FAIL_LOW_SOFTCLIP3P
    return PASS


def has_ZT(aln: pysam.AlignedSegment) -> bool:
    try:
        aln.get_tag("ZT")
        return True
    except KeyError:
        return False


def safe_fetch_all(fh: pysam.AlignmentFile):
    """
    Robust iterator over all records, even if BAM isn't indexed.
    Uses until_eof=True which streams through.
    """
    return fh.fetch(until_eof=True)


FAIL_KEYS = [
    FAIL_UNMAPPED,
    FAIL_SECONDARY_SUPP,
    FAIL_LOW_MAPQ,
    FAIL_LOW_INTRONS,
    FAIL_LOW_SOFTCLIP3P,
]


def _stats_for_one_sample(task):
    """Scan one sample's original BAM + its zt_tagged BAM and return its stats row.

    Picklable ProcessPool entry point. Each sample is independent (reads only its own BAMs and
    returns a small dict of scalars), so results are identical to the serial path.
    """
    (bam, zt_tagged_dir, primary_only, min_mapq, min_introns_read, require_softclip3p) = task

    sample = os.path.basename(bam).replace(".bam", "")
    zt_bam = os.path.join(zt_tagged_dir, f"{sample}.zt_tagged.bam")

    total_lens = new_len_buf()
    considered_lens = new_len_buf()
    assigned_lens = new_len_buf()
    zt_mapped_unassigned_lens = new_len_buf()

    total_n = 0
    total_mapped = 0
    total_unmapped = 0
    considered_n = 0
    fail_counts = {k: 0 for k in FAIL_KEYS}

    # Scan original BAM
    with pysam.AlignmentFile(bam, "rb") as fh:
        for aln in safe_fetch_all(fh):
            total_n += 1
            L = get_len(aln)
            if L > 0:
                total_lens.append(L)

            if aln.is_unmapped:
                total_unmapped += 1
            else:
                total_mapped += 1

            reason = considered_fail_reason(
                aln,
                primary_only=primary_only,
                min_mapq=min_mapq,
                min_introns_read=min_introns_read,
                require_softclip3p=require_softclip3p,
            )

            if reason == PASS:
                considered_n += 1
                if L > 0:
                    considered_lens.append(L)
            else:
                fail_counts[reason] = fail_counts.get(reason, 0) + 1

    # Scan zt_tagged BAM for assigned vs unassigned (within zt_tagged universe)
    assigned_n = 0
    zt_total_records = 0
    zt_unmapped_records = 0
    zt_mapped_records = 0
    zt_mapped_unassigned_n = 0

    if os.path.exists(zt_bam):
        with pysam.AlignmentFile(zt_bam, "rb") as fh:
            for aln in safe_fetch_all(fh):
                zt_total_records += 1
                if aln.is_unmapped:
                    zt_unmapped_records += 1
                    continue
                zt_mapped_records += 1
                L = get_len(aln)
                if has_ZT(aln):
                    assigned_n += 1
                    if L > 0:
                        assigned_lens.append(L)
                else:
                    zt_mapped_unassigned_n += 1
                    if L > 0:
                        zt_mapped_unassigned_lens.append(L)

    s_total = qstats(total_lens)
    s_cons = qstats(considered_lens)
    s_asgn = qstats(assigned_lens)
    s_zt_un = qstats(zt_mapped_unassigned_lens)

    def _frac(num, den):
        return round(num / den, 6) if den else 0.0

    failed_total = sum(fail_counts.values())

    return dict(
        sample=sample,

        # Original BAM totals
        total_reads_bam=total_n,
        total_mapped=total_mapped,
        total_unmapped=total_unmapped,

        # Considered universe (original BAM + filters)
        considered_reads=considered_n,

        # Failure breakdown (original BAM)
        failed_unmapped=fail_counts[FAIL_UNMAPPED],
        failed_secondary_or_supp=fail_counts[FAIL_SECONDARY_SUPP],
        failed_low_mapq=fail_counts[FAIL_LOW_MAPQ],
        failed_low_introns=fail_counts[FAIL_LOW_INTRONS],
        failed_low_softclip3p=fail_counts[FAIL_LOW_SOFTCLIP3P],

        # zt_tagged summaries (if present)
        zt_tagged_exists=int(os.path.exists(zt_bam)),
        zt_total_records=zt_total_records,
        zt_unmapped_records=zt_unmapped_records,
        zt_mapped_records=zt_mapped_records,

        assigned_reads=assigned_n,
        zt_mapped_unassigned_reads=zt_mapped_unassigned_n,

        # Length summaries
        total_len_min=s_total["min"], total_len_p25=s_total["p25"], total_len_p50=s_total["p50"],
        total_len_p75=s_total["p75"], total_len_p90=s_total["p90"], total_len_max=s_total["max"],

        considered_len_min=s_cons["min"], considered_len_p25=s_cons["p25"], considered_len_p50=s_cons["p50"],
        considered_len_p75=s_cons["p75"], considered_len_p90=s_cons["p90"], considered_len_max=s_cons["max"],

        assigned_len_min=s_asgn["min"], assigned_len_p25=s_asgn["p25"], assigned_len_p50=s_asgn["p50"],
        assigned_len_p75=s_asgn["p75"], assigned_len_p90=s_asgn["p90"], assigned_len_max=s_asgn["max"],

        zt_unassigned_len_min=s_zt_un["min"], zt_unassigned_len_p25=s_zt_un["p25"], zt_unassigned_len_p50=s_zt_un["p50"],
        zt_unassigned_len_p75=s_zt_un["p75"], zt_unassigned_len_p90=s_zt_un["p90"], zt_unassigned_len_max=s_zt_un["max"],

        # --- Faithful retention funnel: the fractions the counts above imply, plus mean lengths.
        # total -> mapped -> considered (assembler filters) -> assigned (has ZT).
        frac_mapped_of_total=_frac(total_mapped, total_n),
        frac_considered_of_total=_frac(considered_n, total_n),
        frac_considered_of_mapped=_frac(considered_n, total_mapped),
        frac_assigned_of_considered=_frac(assigned_n, considered_n),
        frac_assigned_of_total=_frac(assigned_n, total_n),
        frac_failed_of_total=_frac(failed_total, total_n),
        failed_reads_total=failed_total,
        total_len_mean=s_total["mean"],
        considered_len_mean=s_cons["mean"],
        assigned_len_mean=s_asgn["mean"],
        zt_unassigned_len_mean=s_zt_un["mean"],
    )


# ----------------------------- main logic -----------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bams-dir", required=True)
    ap.add_argument("--bam-glob", default="*.bam")
    ap.add_argument("--zt-tagged-dir", required=True)
    ap.add_argument("--out", required=True)

    # must match assembler config
    ap.add_argument("--primary-only", action="store_true")
    ap.add_argument("--min-mapq", type=int, default=10)
    ap.add_argument("--min-introns-read", type=int, default=1)
    ap.add_argument("--require-softclip3p", type=int, default=0)
    ap.add_argument("--jobs", type=int, default=1,
                    help="Number of samples to scan in parallel (default 1 = serial, single-core safe).")

    args = ap.parse_args()

    bam_paths = sorted(glob.glob(os.path.join(args.bams_dir, args.bam_glob)))
    if not bam_paths:
        raise SystemExit(f"No BAMs found under {args.bams_dir} with glob {args.bam_glob}")

    tasks = [
        (bam, args.zt_tagged_dir, bool(args.primary_only), int(args.min_mapq),
         int(args.min_introns_read), int(args.require_softclip3p))
        for bam in bam_paths
    ]

    # Samples are independent -> scan them in parallel. jobs<=1 (or a single sample) stays
    # strictly serial, so the stage still runs on one core. ex.map preserves input order, so the
    # output row order always matches the sorted BAM order regardless of completion order.
    jobs = max(1, int(args.jobs))
    rows: List[Dict[str, object]] = []
    if jobs > 1 and len(tasks) > 1:
        try:
            with ProcessPoolExecutor(max_workers=min(jobs, len(tasks))) as ex:
                rows = list(ex.map(_stats_for_one_sample, tasks))
        except Exception as exc:
            print(f"[warn] falling back to serial read stats: {exc}", file=sys.stderr)
            rows = []
    if not rows:
        rows = [_stats_for_one_sample(t) for t in tasks]

    # write TSV
    if not rows:
        raise SystemExit("No rows produced (unexpected).")

    cols = list(rows[0].keys())
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out, "w") as out:
        out.write("\t".join(cols) + "\n")
        for r in rows:
            out.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")


if __name__ == "__main__":
    main()

