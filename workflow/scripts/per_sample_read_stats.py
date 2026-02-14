#!/usr/bin/env python3
import argparse, os, glob
import numpy as np
import pysam

def qstats(lengths):
    if not lengths:
        return dict(n=0, min=0, p25=0, p50=0, p75=0, p90=0, max=0)
    arr = np.asarray(lengths, dtype=np.int64)
    return dict(
        n=int(arr.size),
        min=int(arr.min()),
        p25=int(np.quantile(arr, 0.25, method="linear")),
        p50=int(np.quantile(arr, 0.50, method="linear")),
        p75=int(np.quantile(arr, 0.75, method="linear")),
        p90=int(np.quantile(arr, 0.90, method="linear")),
        max=int(arr.max()),
    )

def intron_count(aln):
    # count N ops in CIGAR
    ct = aln.cigartuples or []
    return sum(1 for op, ln in ct if op == 3)

def softclip3p_len(aln):
    ct = aln.cigartuples or []
    if not ct:
        return 0
    # determine tx strand from FLAG (same as assembler: genomic strand only)
    tx = "-" if aln.is_reverse else "+"
    if tx == "+":
        return ct[-1][1] if ct and ct[-1][0] == 4 else 0
    else:
        return ct[0][1] if ct and ct[0][0] == 4 else 0

def get_len(aln):
    L = aln.query_length
    if L is not None:
        return int(L)
    qs = aln.query_sequence
    return int(len(qs)) if qs else 0

def passes_considered_filters(aln, primary_only, min_mapq, min_introns_read, require_softclip3p):
    if aln.is_unmapped:
        return False
    if primary_only and (aln.is_secondary or aln.is_supplementary):
        return False
    if aln.mapping_quality < min_mapq:
        return False
    if intron_count(aln) < min_introns_read:
        return False
    if require_softclip3p > 0 and softclip3p_len(aln) < require_softclip3p:
        return False
    return True

def has_ZT(aln):
    try:
        aln.get_tag("ZT")
        return True
    except KeyError:
        return False

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

    args = ap.parse_args()

    bam_paths = sorted(glob.glob(os.path.join(args.bams_dir, args.bam_glob)))
    if not bam_paths:
        raise SystemExit(f"No BAMs found under {args.bams_dir} with glob {args.bam_glob}")

    rows = []
    for bam in bam_paths:
        sample = os.path.basename(bam).replace(".bam", "")
        zt_bam = os.path.join(args.zt_tagged_dir, f"{sample}.zt_tagged.bam")

        total_lens = []
        considered_lens = []

        total_n = 0
        total_mapped = 0
        total_unmapped = 0
        considered_n = 0

        with pysam.AlignmentFile(bam, "rb") as fh:
            for aln in fh.fetch(until_eof=True):
                total_n += 1
                L = get_len(aln)
                if L > 0:
                    total_lens.append(L)

                if aln.is_unmapped:
                    total_unmapped += 1
                else:
                    total_mapped += 1

                if passes_considered_filters(
                    aln,
                    primary_only=args.primary_only,
                    min_mapq=args.min_mapq,
                    min_introns_read=args.min_introns_read,
                    require_softclip3p=args.require_softclip3p,
                ):
                    considered_n += 1
                    if L > 0:
                        considered_lens.append(L)

        assigned_n = 0
        assigned_lens = []
        if os.path.exists(zt_bam):
            with pysam.AlignmentFile(zt_bam, "rb") as fh:
                for aln in fh.fetch(until_eof=True):
                    if has_ZT(aln):
                        assigned_n += 1
                        L = get_len(aln)
                        if L > 0:
                            assigned_lens.append(L)
        else:
            # If you didn’t write zt_tagged bams, you can’t do “assigned” this way.
            assigned_n = 0
            assigned_lens = []

        s_total = qstats(total_lens)
        s_cons  = qstats(considered_lens)
        s_asgn  = qstats(assigned_lens)

        rows.append(dict(
            sample=sample,

            total_reads_bam=total_n,
            total_mapped=total_mapped,
            total_unmapped=total_unmapped,

            considered_reads=considered_n,
            assigned_reads=assigned_n,

            total_len_min=s_total["min"], total_len_p25=s_total["p25"], total_len_p50=s_total["p50"],
            total_len_p75=s_total["p75"], total_len_p90=s_total["p90"], total_len_max=s_total["max"],

            considered_len_min=s_cons["min"], considered_len_p25=s_cons["p25"], considered_len_p50=s_cons["p50"],
            considered_len_p75=s_cons["p75"], considered_len_p90=s_cons["p90"], considered_len_max=s_cons["max"],

            assigned_len_min=s_asgn["min"], assigned_len_p25=s_asgn["p25"], assigned_len_p50=s_asgn["p50"],
            assigned_len_p75=s_asgn["p75"], assigned_len_p90=s_asgn["p90"], assigned_len_max=s_asgn["max"],
        ))

    # write TSV
    cols = list(rows[0].keys())
    with open(args.out, "w") as out:
        out.write("\t".join(cols) + "\n")
        for r in rows:
            out.write("\t".join(str(r[c]) for c in cols) + "\n")

if __name__ == "__main__":
    main()
