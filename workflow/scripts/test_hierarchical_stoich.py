#!/usr/bin/env python3
"""Truncation-aware differential stoichiometry between fragmentforms.

THE PROBLEM. Direct-RNA reads are sequenced 3'->5' and truncate at the 5' end. A read that stops
before a gene's 5' region is still ASSIGNED to some fragmentform -- but at any feature the read never
reached, that assignment is INFERRED, not observed. Its modification calls then contribute to a
fragmentform it was only guessed into, so a differential test between fragmentforms that diverge
5'-ward is confounded by reads that carry no evidence about the thing being tested. Measured on this
data: ~35% of assigned reads fall >500nt short of their fragmentform's 5' end, ~22% fall >1kb short.

THE FIX. For each pair of fragmentforms (A, B) find their DIVERGENCE POINT -- the position closest
to the 3' end at which their exon structures differ. Only reads spanning that point can distinguish
A from B by evidence. Restrict the comparison to those reads. As the divergence moves 5'-ward fewer
reads qualify, so power decays; every row reports how many informative reads remained, making an
underpowered 5' feature explicitly underpowered rather than silently "no effect".

WHY IT IS CHEAP. Because all reads share the 3' end and differ only in 5' reach, the informative
sets are NESTED: reach(t') subset-of reach(t) for t' more 5' than t. So "reads spanning position p"
is a PREFIX of the reads sorted by reach, and the whole hierarchy collapses to prefix bitmaps +
`intersection_cardinality` -- the same roaring machinery the association engines already use. Per
gene: one sort, then 4 intersections per (pair, site).

With --also-naive each row additionally carries the unrestricted (all-assigned-reads) result, so the
size of the confound is measurable per site rather than assumed.
"""
import argparse
import os
import re
import shutil
import sys
import tempfile
from bisect import bisect_right

import numpy as np
import pandas as pd
from pyroaring import BitMap

from genotype_utils import (benjamini_hochberg, cmh_stratified_test, mh_stratified_effect,
                            run_contingency_test, shard_tsv_by_chrom, tsv_header)

OUT_COLS = [
    "gene_name", "chrom", "strand", "fragmentform_a", "fragmentform_b",
    "divergence_pos", "divergence_from_3p_nt", "mod_site_id", "site_pos", "site_from_3p_nt",
    "n_informative", "n_informative_a", "n_informative_b", "n_strata_informative",
    "a_modified", "a_unmodified", "b_modified", "b_unmodified",
    "frac_modified_a", "frac_modified_b", "delta", "test_name", "stat_value", "p_value",
    "test_name_pooled", "stat_value_pooled", "p_value_pooled",
    "naive_n_a", "naive_n_b", "naive_frac_a", "naive_frac_b", "naive_delta", "naive_p_value",
    "reads_dropped_as_uninformative", "p_adj_bh",
]
_MOD_WANT = ["sample", "qname", "mod_site_id", "chrom", "start0", "strand", "target_mod_code",
             "state_detail", "gene_name", "ZT", "usable", "fail", "within_alignment"]
_RA_WANT = ["sample", "qname", "chrom", "start0", "end0", "strand", "ZT", "assigned", "gene_name"]


def parse_args():
    ap = argparse.ArgumentParser(description="Differential stoichiometry between fragmentforms using only reads that demonstrably span their divergence point.")
    ap.add_argument("--read-assignments", required=True, help="*_read_assignments*.tsv (per-read extent + ZT)")
    ap.add_argument("--molecule-mods", required=True, help="*_molecule_mod_calls.tsv")
    ap.add_argument("--gtf", required=True, help="Assembled fragmentform GTF (exon structures)")
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--sites", default="", help="Optional TSV of preselected sites (needs mod_site_id) to restrict testing -- the fast path")
    ap.add_argument("--min-informative-reads", type=int, default=10, help="Min informative reads per fragmentform for a pair to be tested")
    ap.add_argument("--min-state-reads", type=int, default=3, help="Min reads in each modification state per fragmentform")
    ap.add_argument("--max-fragmentforms-per-gene", type=int, default=12, help="Cap pairs: keep the N best-supported fragmentforms per gene")
    ap.add_argument("--min-divergence-from-3p", type=int, default=0,
                    help="Only test pairs whose divergence is at least this far (nt) from the 3' end. "
                         "This is the targeted fast path: measured genome-wide, pairs diverging <1kb from "
                         "the 3' end lose reads in 0.2%% of tests and change zero calls (the cheap "
                         "test_diffs already answers those correctly), while >20kb pairs lose reads in "
                         "~55%% of tests. Setting ~5000 skips ~42%% of pairs and keeps essentially all "
                         "of the correction.")
    ap.add_argument("--also-naive", action="store_true", help="Also compute the unrestricted result, to quantify the confound")
    ap.add_argument("--test", choices=["auto", "fisher", "chi2"], default="auto")
    ap.add_argument("--pseudocount", type=float, default=0.5)
    ap.add_argument("--verbose", action="store_true")
    return ap.parse_args()


def load_gtf_exons(path):
    """transcript_id -> (chrom, strand, [(start0,end0), ...]) from the assembled GTF."""
    exons = {}
    meta = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "exon":
                continue
            m = re.search(r'transcript_id "([^"]+)"', f[8])
            if not m:
                continue
            tid = m.group(1)
            exons.setdefault(tid, []).append((int(f[3]) - 1, int(f[4])))
            meta[tid] = (f[0], f[6])
    return {t: (meta[t][0], meta[t][1], sorted(v)) for t, v in exons.items()}


def _covered(intervals, lo, hi):
    """Is [lo,hi) inside any interval? (intervals sorted, non-overlapping)"""
    i = bisect_right([s for s, _ in intervals], lo) - 1
    return i >= 0 and intervals[i][1] >= hi


def divergence_point(ex_a, ex_b, strand):
    """3'-MOST INTERNAL position where the two exon structures differ, or None if none exists.

    A read distinguishes A from B only by spanning a position where they make different exon/intron
    claims, and ONT direct-RNA reads sequence 3'->5' and truncate on the 5' side. Two competing
    concerns set which differing position to gate on:

      * Returning the DEEPEST (5'-most) difference OVER-FILTERS: it drops every read that stopped short
        of a distal difference even though a more 3'-proximal difference had already classified it.
      * Returning the 3'-most difference NAIVELY is a no-op when that difference is a TERMINAL OVERHANG
        (a region beyond one form's 5'/3' transcript bound): there a non-spanning read is genuinely
        ambiguous -- it could be the short form or a truncated long form -- and at the 3' end every read
        trivially spans it (the old measured failure: 380/467 rows reported reads_dropped=0).

    The position that satisfies both is the 3'-MOST difference that is INTERNAL to BOTH forms'
    transcript extents (so spanning it genuinely discriminates), which keeps the most reads informative:
    the largest such coord on '+', the smallest on '-'. Pairs whose only differences are terminal
    overhangs (e.g. pure alt-TSS / alt-TES length) return None and are left to the cheap test_diffs.
    """
    lo_both = max(ex_a[0][0], ex_b[0][0])          # 5'/3'-clipped region internal to BOTH forms
    hi_both = min(ex_a[-1][1], ex_b[-1][1])         # (exons are sorted, so [0][0]=min start, [-1][1]=max end)
    if hi_both <= lo_both:
        return None
    bounds = sorted({p for iv in (ex_a, ex_b) for se in iv for p in se})
    if len(bounds) < 2:
        return None
    diffs = []
    for lo, hi in zip(bounds[:-1], bounds[1:]):
        if hi <= lo:
            continue
        if lo < lo_both or hi > hi_both:            # terminal overhang of one form -> not a valid gate
            continue
        if _covered(ex_a, lo, hi) != _covered(ex_b, lo, hi):
            diffs.append((lo, hi))
    if not diffs:
        return None
    # 3'-most internal difference -> the 5' edge a read must reach: largest lo on '+', smallest hi on '-'
    return max(l for l, _ in diffs) if strand == "+" else min(h for _, h in diffs)


def main():
    args = parse_args()
    gtf = load_gtf_exons(args.gtf)
    if args.verbose:
        print(f"[hier] GTF fragmentforms: {len(gtf):,}", flush=True)

    keep_sites = None
    if args.sites:
        s = pd.read_csv(args.sites, sep="\t", low_memory=False)
        col = "mod_site_id" if "mod_site_id" in s.columns else None
        if col:
            keep_sites = set(s[col].astype(str))
            if args.verbose:
                print(f"[hier] restricting to {len(keep_sites):,} preselected sites", flush=True)

    ra_hdr = tsv_header(args.read_assignments)
    ra_all = pd.read_csv(args.read_assignments, sep="\t", low_memory=False,
                         usecols=[c for c in _RA_WANT if c in ra_hdr])
    ra_all = ra_all[ra_all.get("assigned", True).astype(bool)] if "assigned" in ra_all.columns else ra_all
    ra_all = ra_all[ra_all["ZT"].astype(str).ne("")]
    # ZT is "<gene>.<gene_id>.G#.T#"; GTF transcript_id drops the leading gene name.
    ra_all["gtf_id"] = ra_all["ZT"].astype(str).str.split(".").str[1:].str.join(".")
    ra_all = ra_all[ra_all["gtf_id"].isin(gtf)]

    # The mod table is genotype-scale (tens of GB), so stream it one chromosome at a time -- the same
    # O(1)-RAM sharding the association engines use. Genes live inside a chromosome, so per-chrom
    # processing is exact, not an approximation.
    rows = []
    n_pairs = n_skipped = 0
    tmp = tempfile.mkdtemp(prefix=".hier_", dir=os.path.dirname(args.out_tsv) or ".")
    try:
        shards = shard_tsv_by_chrom(args.molecule_mods, os.path.join(tmp, "mod"))
        if args.verbose:
            print(f"[hier] assigned reads: {len(ra_all):,} | mod shards: {len(shards)}", flush=True)
        for chrom in sorted(shards):
            ra = ra_all[ra_all["chrom"].astype(str) == str(chrom)]
            if ra.empty:
                continue
            mod_hdr = tsv_header(shards[chrom])
            mods = pd.read_csv(shards[chrom], sep="\t", low_memory=False,
                               usecols=[c for c in _MOD_WANT if c in mod_hdr])
            if "usable" in mods.columns:
                mods = mods[mods["usable"].fillna(False)]
            elif "fail" in mods.columns:
                mods = mods[(~mods["fail"].fillna(True)) & mods["within_alignment"].fillna(False)]
            mods = mods[mods["state_detail"].isin(["modified", "canonical", "other_mod"])].copy()
            if mods.empty:
                continue
            mods["target_state"] = mods["state_detail"].eq("modified").astype(int)
            if keep_sites is not None:
                mods = mods[mods["mod_site_id"].astype(str).isin(keep_sites)]
            if mods.empty:
                continue
            r, p, s = _process_chrom(ra, mods, gtf, args)
            rows.extend(r); n_pairs += p; n_skipped += s
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    _finish(rows, args, n_pairs, n_skipped)


def _process_chrom(ra, mods, gtf, args):
    """All genes on one chromosome. Returns (rows, n_pairs, n_skipped)."""
    rows = []
    n_pairs = n_skipped = 0
    mods_by_gene = {g: d for g, d in mods.groupby("gene_name", sort=False)} if "gene_name" in mods.columns else {}

    for gene, rg in ra.groupby("gene_name", sort=False):
        md = mods_by_gene.get(gene)
        if md is None or md.empty:
            continue
        # index this gene's reads
        rg = rg.drop_duplicates(["sample", "qname"])
        key = rg["sample"].astype(str) + "\x00" + rg["qname"].astype(str)
        ridx = {k: i for i, k in enumerate(key)}
        rg = rg.assign(_ri=np.arange(len(rg)))
        strand = str(rg["strand"].iloc[0])
        chrom = str(rg["chrom"].iloc[0])
        # per-sample read bitmaps over the same _ri space, for sample-stratified CMH (sorted for a
        # stable stratum order). Each _ri is one (sample, qname) read, so samples partition _ri.
        sample_bm = {s: BitMap(sg["_ri"].to_numpy().astype(np.uint32))
                     for s, sg in rg.groupby("sample", sort=True)}

        zt_counts = rg["gtf_id"].value_counts()
        zts = list(zt_counts.index[:int(args.max_fragmentforms_per_gene)])
        if len(zts) < 2:
            continue
        zt_bm = {z: BitMap(rg.loc[rg["gtf_id"].eq(z), "_ri"].to_numpy().astype(np.uint32)) for z in zts}

        # reach = how far 5' each read got. Sort so "reads spanning p" is a prefix.
        if strand == "+":
            reach = rg["start0"].to_numpy()          # smaller = further 5'
            order = np.argsort(reach, kind="stable")  # ascending: most-5' first
            sorted_reach = reach[order]
        else:
            reach = rg["end0"].to_numpy()             # larger = further 5'
            order = np.argsort(-reach, kind="stable")
            sorted_reach = reach[order]
        ri_sorted = rg["_ri"].to_numpy()[order]

        def spanning(pos):
            """Bitmap of reads whose 5' reach passes `pos` -- a prefix of the reach-sorted order."""
            if strand == "+":
                k = bisect_right(sorted_reach.tolist(), pos)      # start0 <= pos
            else:
                k = bisect_right((-sorted_reach).tolist(), -pos)  # end0   >= pos
            return BitMap(ri_sorted[:k].astype(np.uint32)), k

        span_cache = {}
        # per-site bitmaps over this gene's reads
        mk = md["sample"].astype(str) + "\x00" + md["qname"].astype(str)
        md = md.assign(_ri=mk.map(ridx).to_numpy())
        md = md[md["_ri"].notna()]
        if md.empty:
            continue
        md = md.assign(_ri=md["_ri"].astype(int))
        site_bm = {}
        for sid, sg in md.groupby("mod_site_id", sort=False):
            ri = sg["_ri"].to_numpy(); ts = sg["target_state"].to_numpy()
            site_bm[sid] = (BitMap(ri[ts == 1].astype(np.uint32)),
                            BitMap(ri[ts == 0].astype(np.uint32)),
                            int(sg["start0"].iloc[0]))
        tes = rg["end0"].max() if strand == "+" else rg["start0"].min()

        for i in range(len(zts)):
            for j in range(i + 1, len(zts)):
                a, b = zts[i], zts[j]
                d = divergence_point(gtf[a][2], gtf[b][2], strand)
                if d is None:
                    continue
                # 3'-proximal pairs are exactly where every read is already informative, so the cheap
                # test answers them identically -- skip them when running the targeted 5' path.
                if abs(tes - d) < int(args.min_divergence_from_3p):
                    continue
                n_pairs += 1
                for sid, (mbm, ubm, spos) in site_bm.items():
                    # a read must span BOTH the divergence point and the site -> the more 5' of the two
                    thr = min(d, spos) if strand == "+" else max(d, spos)
                    if thr not in span_cache:
                        span_cache[thr] = spanning(thr)
                    inf_bm, _ = span_cache[thr]
                    ia = zt_bm[a] & inf_bm
                    ib = zt_bm[b] & inf_bm
                    am = ia.intersection_cardinality(mbm); au = ia.intersection_cardinality(ubm)
                    bm_ = ib.intersection_cardinality(mbm); bu = ib.intersection_cardinality(ubm)
                    # --min-informative-reads applies to the reads that actually ENTER the 2x2 (span the
                    # divergence AND carry a mod call at the site) -- NOT the spanning count len(ia),
                    # which also includes no-call reads that contribute no power. Gating on len(ia) let
                    # pairs with only 6-9 test reads pass a threshold of 10; n_informative_* is now the
                    # test-read count so it matches what was actually tested.
                    n_a = am + au; n_b = bm_ + bu
                    if n_a < args.min_informative_reads or n_b < args.min_informative_reads:
                        n_skipped += 1
                        continue
                    # require min_state_reads in EACH modification state of EACH fragmentform (per the
                    # help text + the sibling scripts) -- checking the per-fragmentform TOTAL let a zero
                    # cell through (e.g. 6 mod / 0 unmod) and inflated the BH family.
                    if min(am, au, bm_, bu) < args.min_state_reads:
                        continue
                    tab = [[float(am), float(au)], [float(bm_), float(bu)]]
                    # POOLED test (kept as *_pooled): pools informative reads across samples, so a
                    # per-sample difference in both modification rate AND a/b composition manufactures a
                    # fragmentform stoichiometry difference where the within-sample effect is ~0.
                    ptname, _psn, psval, ppval = run_contingency_test(tab, test=args.test, pseudocount=args.pseudocount)
                    # SAMPLE-STRATIFIED CMH (primary): one [[am,au],[bm,bu]] per sample among the
                    # informative reads, combined by Cochran-Mantel-Haenszel.
                    ia_m = ia & mbm; ia_u = ia & ubm; ib_m = ib & mbm; ib_u = ib & ubm
                    strata = [[[float(ia_m.intersection_cardinality(sbm)), float(ia_u.intersection_cardinality(sbm))],
                               [float(ib_m.intersection_cardinality(sbm)), float(ib_u.intersection_cardinality(sbm))]]
                              for sbm in sample_bm.values()]
                    if len(sample_bm) > 1:
                        tname, _sn, sval, pval, n_strata = cmh_stratified_test(strata)
                    else:
                        tname, sval, pval, n_strata = ptname, psval, ppval, 0
                    fa = am / (am + au); fb = bm_ / (bm_ + bu)
                    row = {
                        "gene_name": gene, "chrom": chrom, "strand": strand,
                        "fragmentform_a": a, "fragmentform_b": b,
                        "divergence_pos": int(d), "divergence_from_3p_nt": int(abs(tes - d)),
                        "mod_site_id": sid, "site_pos": int(spos), "site_from_3p_nt": int(abs(tes - spos)),
                        "n_informative": int(n_a + n_b),
                        "n_informative_a": int(n_a), "n_informative_b": int(n_b),
                        "n_strata_informative": int(n_strata),
                        "a_modified": am, "a_unmodified": au, "b_modified": bm_, "b_unmodified": bu,
                        "frac_modified_a": round(fa, 5), "frac_modified_b": round(fb, 5),
                        "delta": round(fb - fa, 5), "test_name": tname,
                        "stat_value": sval, "p_value": pval,
                        "test_name_pooled": ptname, "stat_value_pooled": psval, "p_value_pooled": ppval,
                    }
                    if args.also_naive:
                        na, nb = zt_bm[a], zt_bm[b]
                        nam = na.intersection_cardinality(mbm); nau = na.intersection_cardinality(ubm)
                        nbm = nb.intersection_cardinality(mbm); nbu = nb.intersection_cardinality(ubm)
                        nfa = nam / (nam + nau) if (nam + nau) else float("nan")
                        nfb = nbm / (nbm + nbu) if (nbm + nbu) else float("nan")
                        np_ = float("nan")
                        if (nam + nau) and (nbm + nbu):
                            np_ = run_contingency_test([[float(nam), float(nau)], [float(nbm), float(nbu)]],
                                                       test=args.test, pseudocount=args.pseudocount)[3]
                        row.update({"naive_n_a": nam + nau, "naive_n_b": nbm + nbu,
                                    "naive_frac_a": round(nfa, 5), "naive_frac_b": round(nfb, 5),
                                    "naive_delta": round(nfb - nfa, 5), "naive_p_value": np_,
                                    "reads_dropped_as_uninformative": int((nam + nau + nbm + nbu) - (am + au + bm_ + bu))})
                    rows.append(row)

    return rows, n_pairs, n_skipped


def _finish(rows, args, n_pairs, n_skipped):
    out = pd.DataFrame(rows)
    if out.empty:
        out = pd.DataFrame(columns=OUT_COLS)
    else:
        out["p_adj_bh"] = benjamini_hochberg(out["p_value"].values)
        out["_a"] = out["delta"].abs()
        out = out.sort_values(["p_adj_bh", "_a"], ascending=[True, False]).drop(columns="_a")
        for c in OUT_COLS:
            if c not in out.columns:
                out[c] = pd.NA
        out = out[OUT_COLS].reset_index(drop=True)
    os.makedirs(os.path.dirname(args.out_tsv) or ".", exist_ok=True)
    out.to_csv(args.out_tsv, sep="\t", index=False)
    if args.verbose:
        print(f"[hier] pairs considered={n_pairs:,} | (pair,site) skipped for too few informative reads={n_skipped:,}", flush=True)
        print(f"[hier] wrote {len(out):,} tested (pair, site) rows -> {args.out_tsv}", flush=True)


if __name__ == "__main__":
    main()
