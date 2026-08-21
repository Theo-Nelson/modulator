#!/usr/bin/env python3
"""
Streaming ZN aggregator.

modkit's per-ZN bedMethyl files are already position-sorted and tabix-indexed, so
we can aggregate WITHOUT aggregate_by_gene.py's expensive front-end (normalize to a
giant unsorted TSV, then external-sort ~2e9 lines). Instead, per chromosome we
k-way merge the per-(sample,ZN) beds by genomic position in a single streaming
pass. This is dramatically less I/O, is parallel across chromosomes, and is
resumable (a per-chromosome checkpoint), which the sort-based path was not.

Per chromosome:
  - tabix-fetch each numbered bed (one (sample, ZN) each) for the chromosome
  - k-way merge by start; gather all rows at a start; sub-group into sites
    (chrom, start0, end0, strand, mod); for each site collect every (sample, ZN)
    measurement, assign_gene, and emit:
        * dedup.RAW rows + RAW long rows
        * FILTERED rows (whole site) iff ANY row passes row_pass_filter
  - write per-chrom partial files + a .done.<chrom> marker

Then concatenate the per-chrom partials into the final dedup/long files and reuse
aggregate_by_gene.py's (order-agnostic, internally-resorting) per-sample-stats and
per-gene/pivot writers, so those outputs are byte-for-byte the same logic.

Output is content-identical to aggregate_by_gene.py (validate with a sorted diff).
"""
from __future__ import annotations

import argparse
import os
import shutil
import sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed

import pysam

import aggregate_by_gene as agg  # same directory; reuse validated helpers


# ----- worker globals (populated once per worker via the pool initializer) -----
_TX = None
_GENE = None
_BEDS = None
_CFG = None


def _init_worker(gtf_path, beds, cfg):
    global _TX, _GENE, _BEDS, _CFG
    _TX, _GENE = agg.load_gene_intervals_from_gtf(gtf_path, verbose=False)
    _BEDS = beds
    _CFG = cfg


def _bed_rows_for_chrom(path, chrom):
    """Yield parsed bed rows for `chrom` from a tabix-indexed bed, in file order."""
    try:
        tbx = pysam.TabixFile(path)
    except Exception:
        return
    try:
        if chrom not in tbx.contigs:
            return
        for line in tbx.fetch(chrom):
            rec = agg.parse_bed_line(line if line.endswith("\n") else line + "\n")
            if rec:
                yield rec
    finally:
        tbx.close()


def process_chrom(chrom):
    """Stream-merge all beds for one chromosome -> per-chrom partial outputs."""
    tx, gene, beds, cfg = _TX, _GENE, _BEDS, _CFG
    workdir = cfg["workdir"]
    marker = os.path.join(workdir, f".done.{agg.sanitize_filename_token(chrom)}")
    raw_long_p = os.path.join(workdir, f"raw_long.{agg.sanitize_filename_token(chrom)}.tsv")
    filt_long_p = os.path.join(workdir, f"filt_long.{agg.sanitize_filename_token(chrom)}.tsv")
    dedup_raw_p = os.path.join(workdir, f"dedup_raw.{agg.sanitize_filename_token(chrom)}.tsv")
    dedup_filt_p = os.path.join(workdir, f"dedup_filt.{agg.sanitize_filename_token(chrom)}.tsv")
    if os.path.exists(marker):
        return (chrom, "skipped", 0)

    min_cov = cfg["min_cov"]
    cdf = cfg["count_diff_factor"]
    k_default = cfg.get("k_default", 1.0)
    k_per_mod = cfg.get("k_per_mod", {})
    filter_enable = cfg["filter_enable"]

    # one position-sorted iterator per bed that has this chrom
    streams = []
    for (_root, sample, path, zn) in beds:
        it = _bed_rows_for_chrom(path, chrom)
        head = next(it, None)
        if head is not None:
            streams.append([head, it, sample, str(int(zn))])

    n_sites = 0
    emit_raw = cfg["emit_raw"]
    # When site-filtering is OFF, the FILTERED outputs alias RAW (filtered == everything), matching the
    # sort engine (aggregate_by_gene.py). So we must still produce the RAW partials in that case, even
    # if emit_raw is off, because they are the source the FILTERED outputs are copied from.
    need_raw = emit_raw or (cfg.get("emit_filtered", True) and not filter_enable)
    # RAW partials are large and feed only the (unconsumed) RAW outputs; skip them
    # entirely when they are not needed so the on-disk footprint is just FILTERED.
    raw_long = open(raw_long_p, "w") if need_raw else None
    filt_long = open(filt_long_p, "w") if filter_enable else None
    dedup_raw = open(dedup_raw_p, "w") if need_raw else None
    dedup_filt = open(dedup_filt_p, "w") if filter_enable else None

    try:
        while streams:
            min_start = min(s[0]["start0"] for s in streams)
            rows_at_start = []  # (rec, sample, zn)
            still = []
            for s in streams:
                head, it, sample, zn = s
                if head is not None and head["start0"] == min_start:
                    while head is not None and head["start0"] == min_start:
                        rows_at_start.append((head, sample, zn))
                        head = next(it, None)
                    s[0] = head
                if s[0] is not None:
                    still.append(s)
            streams = still

            # sub-group rows at this start into sites = (end0, strand, mod_code)
            sites = defaultdict(list)
            for rec, sample, zn in rows_at_start:
                sites[(rec["end0"], rec["strand"], rec["mod_code"])].append((rec, sample, zn))

            for (end0, strand, mod), members in sites.items():
                # Sum any duplicate rows per (sample, ZN) to match the original
                # dedup_reduce_sorted (its key is sample,zn,chrom,start,end,strand,mod,gid,gname;
                # gid/gname are a deterministic function of chrom,pos,strand,zn).
                ag_counts = {}
                for rec, sample, zn in members:
                    v = ag_counts.setdefault((sample, zn), [0] * 8)
                    v[0] += int(rec["Nvalid_cov"]); v[1] += int(rec["Nmod"]); v[2] += int(rec["Ncanonical"])
                    v[3] += int(rec["Nother_mod"]); v[4] += int(rec["Ndelete"]); v[5] += int(rec["Nfail"])
                    v[6] += int(rec["Ndiff"]); v[7] += int(rec["Nnocall"])
                site_pass = False
                # k-ratio threshold is per mod_code (constant across this site's sample/ZN rows)
                site_k = agg.resolve_nfail_score_k(mod, k_default, k_per_mod)
                buf_dedup = []
                buf_long = []
                s0 = str(min_start); e0 = str(end0)
                for (sample, zn) in sorted(ag_counts, key=lambda t: (t[0], int(t[1]))):
                    cov, nmod, ncan, nother, ndel, nfail, ndiff, nnoc = ag_counts[(sample, zn)]
                    gid, gname = agg.assign_gene(chrom, min_start, end0, strand, int(zn), tx, gene)
                    frac = agg.frac_modified(nmod, cov, min_cov)
                    buf_dedup.append("\t".join([
                        sample, zn, chrom, s0, e0, strand, mod, gid, gname,
                        str(cov), str(nmod), str(ncan), str(nother), str(ndel),
                        str(nfail), str(ndiff), str(nnoc), f"{frac:.6f}",
                    ]) + "\n")
                    buf_long.append(
                        f"{sample}\t{zn}\t{chrom}\t{s0}\t{e0}\t{strand}\t{mod}\t"
                        f"{cov}\t{nmod}\t{frac:.6f}\t{gid}\t{gname}\t"
                        f"{ncan}\t{nother}\t{ndel}\t{nfail}\t{ndiff}\t{nnoc}\n"
                    )
                    if filter_enable and agg.row_pass_filter(cov, nmod, nfail, ndiff, cdf, site_k):
                        site_pass = True
                if need_raw:
                    dedup_raw.writelines(buf_dedup)
                    raw_long.writelines(buf_long)
                if filter_enable and site_pass:
                    dedup_filt.writelines(buf_dedup)
                    filt_long.writelines(buf_long)
                n_sites += 1
    finally:
        if raw_long: raw_long.close()
        if dedup_raw: dedup_raw.close()
        if filt_long: filt_long.close()
        if dedup_filt: dedup_filt.close()

    with open(marker, "w") as m:
        m.write(f"{chrom}\t{n_sites}\n")
    return (chrom, "done", n_sites)


def _chroms_from_beds(beds):
    chroms = set()
    for (_root, _sample, path, _zn) in beds:
        try:
            tbx = pysam.TabixFile(path)
            chroms.update(tbx.contigs)
            tbx.close()
        except Exception:
            continue
    return sorted(chroms)


def _concat(parts, out_path, header=None):
    with open(out_path, "w") as out:
        if header is not None:
            out.write("\t".join(header) + "\n")
        for p in parts:
            if os.path.exists(p) and os.path.getsize(p) > 0:
                with open(p, "r") as f:
                    shutil.copyfileobj(f, out, length=1024 * 1024)


def parse_args():
    ap = argparse.ArgumentParser(description="Streaming ZN aggregator (merge pre-sorted, tabix-indexed beds).")
    ap.add_argument("--modkit-dir", required=True)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--min-cov", type=int, default=0)
    ap.add_argument("--tmpdir", default=None)
    ap.add_argument("--chunk-lines", type=int, default=2_000_000)
    ap.add_argument("--count-diff-factor", type=float, default=3.0)
    ap.add_argument("--mod-fail-margin", type=int, default=1,
                    help="DEPRECATED / no-op: superseded by --nfail-score-k (k=1 reproduces margin=0).")
    ap.add_argument("--jobs", type=int, default=8, help="chromosomes processed in parallel")
    ap.add_argument("--verbose", action="store_true")
    ap.add_argument("--filter-enable", dest="filter_enable", action="store_true", default=False)
    ap.add_argument("--nfail-score-k", type=str, default="1.0",
                    help="NFail-SCORE k-ratio confident-call filter: FAIL if Nmod < k*(Nfail+1). Either a single "
                         "value ('1.0') or a per-mod-code map ('a=0.4,17802=1.0,default=1.0'). 0 disables.")
    for flag in ("emit-raw", "emit-filtered", "write-long",
                 "write-raw-per-gene", "write-filtered-per-gene"):
        dest = flag.replace("-", "_")
        ap.add_argument(f"--{flag}", dest=dest, action="store_true", default=True)
        ap.add_argument(f"--no-{flag}", dest=dest, action="store_false")
    # Per-gene pivots are optional inspection outputs (3 dense files per gene x mod group);
    # nothing downstream reads them. 'auto' skips them past --pivot-max-groups so a whole-
    # transcriptome run does not explode into hundreds of thousands of tiny files; 'on' forces
    # them even at scale; 'off' never writes them. Legacy --write-pivots/--no-write-pivots alias.
    ap.add_argument("--pivot-mode", dest="pivot_mode", choices=["auto", "on", "off"],
                    default="auto")
    ap.add_argument("--pivot-max-groups", dest="pivot_max_groups", type=int, default=2000)
    ap.add_argument("--write-pivots", dest="pivot_mode", action="store_const", const="on")
    ap.add_argument("--no-write-pivots", dest="pivot_mode", action="store_const", const="off")
    return ap.parse_args()


def main():
    args = parse_args()
    base = args.out_prefix
    agg.ensure_dir(os.path.dirname(base) or ".")
    # Stable per-prefix workdir so re-runs resume completed chromosomes.
    workdir = f"{base}__zn_stream_work"
    agg.ensure_dir(workdir)

    beds = agg.iter_numbered_beds(args.modkit_dir)
    if not beds:
        raise SystemExit(f"No numbered ZN beds under {args.modkit_dir}")
    chroms = _chroms_from_beds(beds)
    if args.verbose:
        print(f"[stream] {len(beds)} beds, {len(chroms)} chromosomes, jobs={args.jobs}", file=sys.stderr, flush=True)

    _k_default, _k_per_mod = agg.parse_nfail_score_k(args.nfail_score_k)
    cfg = dict(workdir=workdir, min_cov=args.min_cov, count_diff_factor=args.count_diff_factor,
               k_default=_k_default, k_per_mod=_k_per_mod,
               filter_enable=args.filter_enable, emit_raw=args.emit_raw,
               emit_filtered=args.emit_filtered)

    # ---- Phase 1: per-chromosome streaming merge (parallel, resumable) ----
    todo = [c for c in chroms if not os.path.exists(os.path.join(workdir, f".done.{agg.sanitize_filename_token(c)}"))]
    if args.verbose:
        print(f"[stream] phase1: {len(chroms)-len(todo)} chrom(s) already done, {len(todo)} to do", file=sys.stderr, flush=True)
    if todo:
        jobs = max(1, min(args.jobs, len(todo)))
        if jobs == 1:
            _init_worker(args.gtf, beds, cfg)
            for c in todo:
                r = process_chrom(c)
                if args.verbose:
                    print(f"[stream] {r[0]}: {r[1]} ({r[2]} sites)", file=sys.stderr, flush=True)
        else:
            with ProcessPoolExecutor(max_workers=jobs, initializer=_init_worker,
                                     initargs=(args.gtf, beds, cfg)) as ex:
                futs = {ex.submit(process_chrom, c): c for c in todo}
                for fut in as_completed(futs):
                    r = fut.result()
                    if args.verbose:
                        print(f"[stream] {r[0]}: {r[1]} ({r[2]} sites)", file=sys.stderr, flush=True)

    # ---- Phase 2: concatenate per-chrom partials ----
    def parts(kind):
        return [os.path.join(workdir, f"{kind}.{agg.sanitize_filename_token(c)}.tsv") for c in chroms]

    # When site-filtering is OFF the FILTERED outputs alias RAW (filtered == everything), exactly as the
    # sort engine (aggregate_by_gene.py) does — otherwise downstream stages that key on the FILTERED
    # tables (e.g. test_diffs) would be silently skipped under the default stream engine.
    filt_aliases_raw = args.emit_filtered and not args.filter_enable
    need_raw = args.emit_raw or filt_aliases_raw
    dedup_raw = os.path.join(workdir, "dedup.RAW.tsv")
    dedup_filt = os.path.join(workdir, "dedup.FILTERED.tsv")
    if need_raw:
        _concat(parts("dedup_raw"), dedup_raw)
    if args.filter_enable:
        _concat(parts("dedup_filt"), dedup_filt)
    elif filt_aliases_raw:
        dedup_filt = dedup_raw   # FILTERED per-gene/stats read the RAW dedup

    if args.emit_raw and args.write_long:
        _concat(parts("raw_long"), f"{base}_RAW_sites_long.tsv", header=agg.LONG_HEADER)
    if args.emit_filtered and args.write_long:
        # FILTERED long from the filtered partials when filtering is on; from the raw partials when it
        # is off (so the FILTERED long is present-and-complete rather than missing).
        _concat(parts("filt_long" if args.filter_enable else "raw_long"),
                f"{base}_FILTERED_sites_long.tsv", header=agg.LONG_HEADER)
    if args.verbose:
        print("[stream] phase2: concatenated long + dedup outputs", file=sys.stderr, flush=True)

    # ---- Phase 3: reuse the (order-agnostic) per-sample-stats + per-gene/pivot writers ----
    # The per-gene/pivot writers are parallelized across genes (jobs>1); jobs<=1 stays serial.
    if args.emit_raw:
        agg.compute_per_sample_mod_stats_from_dedup(dedup_raw, base, "RAW", workdir, args.chunk_lines, args.verbose)
        agg.generate_per_gene_outputs_from_dedup(dedup_raw, base, "RAW", args.write_raw_per_gene,
                                                 args.pivot_mode, workdir, args.chunk_lines, args.verbose,
                                                 jobs=args.jobs, pivot_max_groups=args.pivot_max_groups)
    if args.emit_filtered:
        agg.compute_per_sample_mod_stats_from_dedup(dedup_filt, base, "FILTERED", workdir, args.chunk_lines, args.verbose)
        agg.generate_per_gene_outputs_from_dedup(dedup_filt, base, "FILTERED", args.write_filtered_per_gene,
                                                 args.pivot_mode, workdir, args.chunk_lines, args.verbose,
                                                 jobs=args.jobs, pivot_max_groups=args.pivot_max_groups)

    # success -> drop the (large) workdir
    shutil.rmtree(workdir, ignore_errors=True)
    if args.verbose:
        print("[stream] done.", file=sys.stderr, flush=True)


if __name__ == "__main__":
    main()
