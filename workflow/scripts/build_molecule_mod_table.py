#!/usr/bin/env python3

import argparse
import gzip
import os
import re
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict

import numpy as np
import pandas as pd
import pysam

from genotype_utils import (load_read_assignments, normalize_string_series, robust_load_summary,
                            run_process_jobs, sample_name_from_bam, safe_float, safe_int, write_chrom_index)

# Column order of the populated table (== what the read_assignments join used to produce).
FINAL_COLUMNS = [
    "sample", "qname", "mod_site_id", "chrom", "start0", "end0", "strand",
    "target_mod_code", "call_code", "state_detail", "target_modified",
    "call_prob", "canonical_base", "modified_primary_base", "fail",
    "within_alignment", "gene_id", "gene_name", "metagene_index",
    "ZT", "ZG", "ZN", "ZM", "assigned", "assignment_gene_id", "assignment_gene_name",
    "gene_index", "transcript_index", "assignment_metagene_index", "classification", "usable",
]


OUTPUT_COLUMNS = [
    "sample", "qname", "mod_site_id", "chrom", "start0", "end0", "strand",
    "target_mod_code", "call_code", "state_detail", "target_modified",
    "call_prob", "canonical_base", "modified_primary_base", "fail",
    "within_alignment", "gene_id", "gene_name", "metagene_index",
    # The populated path joins the read->fragmentform assignment (below) and derives `usable`; the
    # empty-output header must advertise the SAME columns so downstream readers see a consistent
    # schema whether or not the table has rows (they read by name, so order is immaterial here).
    "ZT", "ZG", "ZN", "ZM", "assigned", "gene_index", "transcript_index", "classification",
    "assignment_gene_id", "assignment_gene_name", "assignment_metagene_index", "usable",
]

# transcript-oriented base each modification sits on (canonical_base in modkit's output is already
# strand-adjusted -- verified: m6A rows are 'A' on both strands). A read whose canonical_base is NOT
# this base carries a variant at the site and CANNOT carry the modification, so recording it as a
# (usable, unmodified) observation manufactures false negative allele-specific-modification signal.
MOD_BASE = {"a": "A", "17596": "A", "69426": "A", "m": "C", "19228": "C",
            "17802": "T", "19227": "T", "19229": "G", "h": "C", "f": "C", "c": "C"}


def _base_mismatch(target_mod, canonical_base):
    """True iff the read's canonical base cannot carry `target_mod` (a variant at the modified base)."""
    exp = MOD_BASE.get(str(target_mod))
    cb = str(canonical_base or "").upper()
    return exp is not None and cb != "" and cb != exp


_COMP = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N", "U": "A"}
_MM_ENTRY_RE = re.compile(r'^([ACGTUNacgtun])([-+])([a-z]+|[0-9]+)([.?]?)')


def parse_mm_groups(read):
    """Parse the MM/Mm tag. Returns (groups, has_listed_deltas):
      groups            = {(canonical_base, mod_code): is_implicit}
      has_listed_deltas = True iff ANY group declares >=1 delta (a listed position)

    BLOCKER-4: a modBAM in IMPLICIT MM mode (flag '.' or absent) declares that every canonical base of
    a group NOT listed in the deltas is an implicitly-canonical (unmodified) observation. read.modified_bases
    returns ONLY the listed positions, so the pysam backend silently dropped those unmodified calls and
    inflated every modified fraction. Explicit mode ('?') means unlisted positions carry NO call.

    Keyed by (base, code) NOT code alone: a code can appear on two canonical bases (e.g. C+m. and A+m?),
    and the implicit pass must know the base to emit a canonical call ONLY at that base -- and to reject
    codes whose base the read does not carry, for ANY code (the old MOD_BASE allowlist could not).
    has_listed_deltas lets the caller detect an htslib PARSE FAILURE (deltas present but modified_bases
    empty) and skip the read instead of calling it entirely canonical."""
    mm = None
    for tag in ("MM", "Mm"):
        try:
            mm = read.get_tag(tag)
            break
        except KeyError:
            continue
    if not mm:
        return {}, False
    groups = {}
    has_listed = False
    for entry in str(mm).split(";"):
        entry = entry.strip()
        if not entry:
            continue
        m = _MM_ENTRY_RE.match(entry)
        if not m:
            continue
        base, strand, mods, flag = m.groups()
        # The MM strand marker is NO LONGER discarded: '+' means the calls are relative to the read as
        # basecalled (the dRNA norm), '-' means they refer to the complementary strand (rare, duplex).
        # Implicit-canonical synthesis walks the read's OWN bases, so it is only valid for '+' groups;
        # synthesizing for a '-' group would emit canonical calls at the wrong base. A '-' group's LISTED
        # calls still arrive via pysam's read.modified_bases (which carries strand) and are unaffected.
        implicit = (flag != "?") and (strand == "+")   # '.'/absent -> implicit; '?' -> explicit; '-' -> never implicit
        base = base.upper()
        if entry[m.end():].lstrip().startswith(","):   # a delta list follows the header -> listed positions
            has_listed = True
        for code in re.findall(r'[0-9]+|[a-z]', mods):   # 'mh' -> m,h ; '17802' -> 17802
            key = (base, code)
            # OR-merge so a later '-'/explicit group cannot clobber a '+' implicit True for the same
            # (base, code) (and vice-versa the two-base case C+m./A+m? keeps its own per-key flag).
            groups[key] = groups.get(key, False) or implicit
    return groups, has_listed


def parse_args():
    ap = argparse.ArgumentParser(description="Build a per-read mod call table at candidate modulator sites.")
    ap.add_argument("--bams", nargs="+", required=True, help="Input BAMs with MM/ML tags")
    ap.add_argument("--candidate-sites-tsv", required=True, help="Candidate mod site TSV")
    ap.add_argument("--candidate-bed", required=True, help="Candidate mod BED for modkit include-bed")
    ap.add_argument("--read-assignments", required=True, help="Read assignment TSV")
    ap.add_argument("--summary-tsv", default="",
                    help="Classification summary TSV (zt_label -> gene/metagene/classification). With the "
                         "pysam backend the per-read ZT/ZG/ZN/ZM tags are read from the BAM and joined to "
                         "this small table instead of loading the whole read-assignment table.")
    ap.add_argument("--reference-fa", required=True, help="Reference FASTA")
    ap.add_argument("--out-tsv", required=True, help="Output TSV")
    ap.add_argument("--modkit-bin", default="modkit", help="modkit executable")
    ap.add_argument("--threads", type=int, default=2)
    ap.add_argument("--jobs", type=int, default=1, help="Number of extract shards to run in parallel")
    ap.add_argument("--window-bp", type=int, default=1_000_000,
                    help="Split each chromosome's candidate sites into shards spanning at most this "
                         "many bp so a single heavy chromosome no longer serializes the whole run (B2).")
    ap.add_argument("--interval-size", type=int, default=20000,
                    help="modkit extract --interval-size (bounds modkit's per-extract RSS).")
    ap.add_argument("--chunk-rows", type=int, default=250000,
                    help="Flush a worker's rows to a numbered pickle part every this many rows, so a "
                         "single deep-locus window can't blow up the worker's RSS (keeps the pool alive).")
    ap.add_argument("--pre-extracted", nargs="*", default=None, metavar="SAMPLE=CALLS_TSV",
                    help="If given, skip running modkit and instead parse these already-extracted "
                         "per-sample `modkit extract calls` TSVs (one per subset BAM, produced by a "
                         "separate per-sample sbatch across nodes). Each entry is SAMPLE=path. The "
                         "per-site parsing/join/sort is identical to the in-line modkit path.")
    ap.add_argument("--pysam", action="store_true",
                    help="Extract per-molecule calls with the built-in pysam streaming reader instead "
                         "of modkit (no external modkit, no reference FASTA). Parallelised per "
                         "(BAM x chromosome) via --jobs; each task streams one read at a time so peak "
                         "RSS is ~100MB and it never OOMs (chr15 included) -- no windowing/interval-size "
                         "needed. Emits implicit-canonical calls (parse_mm_groups) so it matches modkit on "
                         "IMPLICIT-MM BAMs -- validated on real chrEBV: identical row count, canonical count "
                         "identical, 4/686765 rows differ at a float32 argmax tie-break (Jaccard 1.0000).")
    ap.add_argument("--verbose", action="store_true")
    return ap.parse_args()


def parse_bool_text(x) -> bool:
    return str(x).strip().lower() in {"1", "true", "t", "yes", "y"}


def iter_site_windows(lookup, window_bp):
    """B2: split each chrom's candidate positions into coordinate windows spanning at most window_bp.
    Every candidate position lands in exactly one window's lookup, so shards never double-count a call
    (the per-window lookup filter is the guard even if two windows' modkit --regions touched). Yields
    (chrom, region_1based, window_lookup)."""
    window_bp = max(1, int(window_bp))
    by_chrom = defaultdict(list)
    for (chrom, pos) in lookup.keys():
        by_chrom[chrom].append(pos)
    for chrom in sorted(by_chrom):
        positions = sorted(set(by_chrom[chrom]))
        cur = []
        cur_lo = None
        for p in positions:
            if cur and (p - cur_lo) >= window_bp:
                region = f"{chrom}:{cur_lo + 1}-{cur[-1] + 1}"
                yield chrom, region, {(chrom, x): lookup[(chrom, x)] for x in cur}
                cur = []
                cur_lo = None
            if cur_lo is None:
                cur_lo = p
            cur.append(p)
        if cur:
            region = f"{chrom}:{cur_lo + 1}-{cur[-1] + 1}"
            yield chrom, region, {(chrom, x): lookup[(chrom, x)] for x in cur}


def extract_rows_from_bam(
    bam: str,
    candidate_bed: str,
    reference_fa: str,
    modkit_bin: str,
    threads_per_job: int,
    lookup,
    region,
    shard_path: str,
    interval_size: int,
    chunk_rows: int,
    verbose: bool = False,
):
    """A: write this (BAM x window) shard's rows straight to disk (pickle) instead of returning a
    Python list -- eliminates the millions-of-dicts IPC back to the parent. Rows are flushed to
    numbered pickle parts every `chunk_rows`, so a single deep-locus window can NOT blow up the
    worker's RSS (the bug that wedged the 16-way pool at 256 GiB). Returns (chrom, [parts], nrows)."""
    sample = sample_name_from_bam(bam)
    chrom = region.split(":", 1)[0] if region else ""
    chunk_rows = max(1, int(chunk_rows))
    rows = []
    parts = []
    total = 0

    def _flush():
        if rows:
            p = f"{shard_path}.{len(parts)}.pkl"
            pd.DataFrame(rows).to_pickle(p)
            parts.append(p)
            rows.clear()
    with tempfile.NamedTemporaryFile(prefix=f"{sample}.extract_calls.", suffix=".tsv.bgz", delete=False) as tmp:
        tmp_path = tmp.name

    try:
        cmd = [
            modkit_bin, "extract", "calls", bam, tmp_path,
            "--bgzf",
            "--force",
            *(["--region", str(region)] if region else []),
            # Cap modkit's per-chunk read buffering. The default 100kb interval over deep
            # direct-RNA piles up huge memory at highly-expressed loci (chr1/chr19 etc.) and
            # OOM'd even a 1TB node when many shards ran concurrently. Smaller chunks bound
            # peak RSS (more overhead, identical output).
            "--interval-size", str(max(1, int(interval_size))),
            # Don't estimate a pass-threshold by sampling reads: on sparse inputs
            # (e.g. region subsets, low-coverage samples) modkit aborts with
            # "Error! not enough datapoints" when there are too few mod calls over
            # the candidate-site BED. All calls are emitted; downstream genotype
            # logic applies its own coverage/quality filters.
            "--no-filtering",
            "--include-bed", candidate_bed,
            "--reference", reference_fa,
            "--mapped-only",
            "--threads", str(max(1, int(threads_per_job))),
            "--out-threads", "1",
            "--suppress-progress",
        ]
        if verbose:
            print(f"[info] mod extract start: {sample} region={region} threads={max(1, int(threads_per_job))}", file=sys.stderr, flush=True)
        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode != 0:
            raise SystemExit(f"modkit extract calls failed for {bam}:\n{proc.stderr}")

        header = None
        with gzip.open(tmp_path, "rt") as fh:
            for line in fh:
                line = line.rstrip("\n")
                if not line:
                    continue
                if header is None:
                    header = line.split("\t")
                    continue
                fields = line.split("\t")
                if len(fields) != len(header):
                    continue
                rec = dict(zip(header, fields))
                rchrom = str(rec.get("chrom", ""))
                start0 = safe_int(rec.get("ref_position", -1), default=-1)
                qname = str(rec.get("read_id", ""))
                call_code = str(rec.get("call_code", ""))
                ref_strand = str(rec.get("ref_strand", ""))
                key = (rchrom, start0)
                if key not in lookup:
                    continue
                for site in lookup[key]:
                    site_strand = str(site.get("strand", ""))
                    if site_strand and ref_strand and ref_strand not in {".", "?"} and site_strand != ref_strand:
                        continue
                    target_mod = str(site["mod_code"])
                    if _base_mismatch(target_mod, rec.get("canonical_base", "")):
                        continue   # variant at the modified base -> read cannot carry the mod (not an observation)
                    if call_code == target_mod:
                        state_detail = "modified"
                        target_modified = 1
                    elif call_code == "-":
                        state_detail = "canonical"
                        target_modified = 0
                    else:
                        state_detail = "other_mod"
                        target_modified = 0
                    rows.append({
                        "sample": sample,
                        "qname": qname,
                        "mod_site_id": site["mod_site_id"],
                        "chrom": rchrom,
                        "start0": start0,
                        "end0": safe_int(site.get("end0", start0 + 1), default=start0 + 1),
                        "strand": site_strand or ref_strand,
                        "target_mod_code": target_mod,
                        "call_code": call_code,
                        "state_detail": state_detail,
                        "target_modified": target_modified,
                        "call_prob": safe_float(rec.get("call_prob", 0.0)),
                        "canonical_base": str(rec.get("canonical_base", "")),
                        "modified_primary_base": str(rec.get("modified_primary_base", "")),
                        "fail": parse_bool_text(rec.get("fail", False)),
                        "within_alignment": parse_bool_text(rec.get("within_alignment", True)),
                        "gene_id": str(site.get("gene_id", "")),
                        "gene_name": str(site.get("gene_name", "")),
                        "metagene_index": str(site.get("metagene_index", "")),
                    })
                    total += 1
                    if len(rows) >= chunk_rows:
                        _flush()
    finally:
        if os.path.exists(tmp_path):
            os.unlink(tmp_path)

    _flush()  # remaining rows below the last chunk boundary
    if verbose:
        print(f"[info] mod extract done: {sample} region={region} rows={total}", file=sys.stderr, flush=True)
    return chrom, parts, total


def parse_extracted_calls(sample, calls_tsv, lookup, shard_dir, chunk_rows, verbose=False):
    """Parse a pre-extracted per-sample `modkit extract calls` TSV (bgzip/plain) into per-chrom
    pickle shards, using the identical per-site expansion + strand filter as extract_rows_from_bam.
    modkit output is coordinate-sorted, so a shard is flushed whenever the chromosome changes or the
    buffer reaches chunk_rows; each shard holds rows from a single chromosome (assembly concatenates
    all parts per chrom regardless of order). Returns list of (chrom, [pickle parts], nrows)."""
    from collections import defaultdict as _dd
    by_chrom = _dd(lambda: {"parts": [], "n": 0})
    part_ct = _dd(int)
    rows = []
    cur = [None]

    def _flush():
        if not rows:
            return
        ch = cur[0]
        p = os.path.join(shard_dir, f"{sample}.{ch}.{part_ct[ch]}.pkl")
        part_ct[ch] += 1
        pd.DataFrame(rows).to_pickle(p)
        by_chrom[ch]["parts"].append(p)
        by_chrom[ch]["n"] += len(rows)
        rows.clear()

    total = 0
    opener = gzip.open if str(calls_tsv).endswith((".gz", ".bgz")) else open
    with opener(calls_tsv, "rt") as fh:
        header = None
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if header is None:
                header = line.split("\t")
                continue
            fields = line.split("\t")
            if fields == header:
                # per-chromosome extracts are concatenated into one per-sample file, so the modkit
                # header line reappears at each chunk boundary -- skip those repeats.
                continue
            if len(fields) != len(header):
                continue
            rec = dict(zip(header, fields))
            rchrom = str(rec.get("chrom", ""))
            start0 = safe_int(rec.get("ref_position", -1), default=-1)
            qname = str(rec.get("read_id", ""))
            call_code = str(rec.get("call_code", ""))
            ref_strand = str(rec.get("ref_strand", ""))
            key = (rchrom, start0)
            if key not in lookup:
                continue
            if cur[0] is not None and rchrom != cur[0]:
                _flush()
            cur[0] = rchrom
            for site in lookup[key]:
                site_strand = str(site.get("strand", ""))
                if site_strand and ref_strand and ref_strand not in {".", "?"} and site_strand != ref_strand:
                    continue
                target_mod = str(site["mod_code"])
                if _base_mismatch(target_mod, rec.get("canonical_base", "")):
                    continue   # variant at the modified base -> read cannot carry the mod (not an observation)
                if call_code == target_mod:
                    state_detail = "modified"
                    target_modified = 1
                elif call_code == "-":
                    state_detail = "canonical"
                    target_modified = 0
                else:
                    state_detail = "other_mod"
                    target_modified = 0
                rows.append({
                    "sample": sample,
                    "qname": qname,
                    "mod_site_id": site["mod_site_id"],
                    "chrom": rchrom,
                    "start0": start0,
                    "end0": safe_int(site.get("end0", start0 + 1), default=start0 + 1),
                    "strand": site_strand or ref_strand,
                    "target_mod_code": target_mod,
                    "call_code": call_code,
                    "state_detail": state_detail,
                    "target_modified": target_modified,
                    "call_prob": safe_float(rec.get("call_prob", 0.0)),
                    "canonical_base": str(rec.get("canonical_base", "")),
                    "modified_primary_base": str(rec.get("modified_primary_base", "")),
                    "fail": parse_bool_text(rec.get("fail", False)),
                    "within_alignment": parse_bool_text(rec.get("within_alignment", True)),
                    "gene_id": str(site.get("gene_id", "")),
                    "gene_name": str(site.get("gene_name", "")),
                    "metagene_index": str(site.get("metagene_index", "")),
                })
                total += 1
                if len(rows) >= chunk_rows:
                    _flush()
    _flush()
    if verbose:
        print(f"[info] parsed pre-extracted calls: {sample} rows={total} chroms={len(by_chrom)}",
              file=sys.stderr, flush=True)
    return [(ch, d["parts"], d["n"]) for ch, d in by_chrom.items()]


def _aligned_qr(read):
    """Query and reference positions of the aligned (M/=/X) bases as sorted int64 arrays -- the same
    (query index -> reference position) map as get_reference_positions(full_length=True) restricted to
    aligned bases, built per CIGAR op instead of per base."""
    q = 0
    r = read.reference_start
    qs = []
    rs = []
    for op, ln in (read.cigartuples or ()):
        if op == 0 or op == 7 or op == 8:
            qs.append(np.arange(q, q + ln)); rs.append(np.arange(r, r + ln)); q += ln; r += ln
        elif op == 1 or op == 4:      # I, S consume the query only
            q += ln
        elif op == 2 or op == 3:      # D, N consume the reference only
            r += ln
    if not qs:
        return None, None
    return np.concatenate(qs), np.concatenate(rs)


def _read_tags(read):
    def _g(tag, default=""):
        try:
            return read.get_tag(tag)
        except Exception:
            return default
    return (str(_g("ZT", "")), safe_int(_g("ZG", "")), safe_int(_g("ZN", "")), safe_int(_g("ZM", "")))


def extract_rows_pysam(bam, chrom_lookup, chrom, shard_path, chunk_rows, verbose=False, window_bp=1_000_000):
    """Stream one chromosome of a modBAM with pysam and emit the same per-(read, candidate site) rows
    the modkit path produces -- one read at a time, so peak RSS is ~100MB regardless of BAM/chrom size
    (never OOMs). Reproduces modkit `extract calls --no-filtering --mapped-only` semantics -- INCLUDING
    implicit-canonical calls (unlisted canonical bases of an implicit MM group), which the previous
    version dropped, inflating every modified fraction on real ONT data (BLOCKER-4). Namely:
    call_prob=(ML+0.5)/256 (float32), canonical=1-sum(mod_probs), call_code=argmax, strand-aware.

    Only the candidate sites inside each read's span are visited (searchsorted into the chromosome's
    sorted candidate positions + a per-CIGAR-op query/reference map), instead of every base of every
    read. The read's ZT/ZG/ZN/ZM tags are carried on each row so the read-assignment table is no
    longer joined. Rows are flushed to pickle parts bucketed by `window_bp` of start0, so the parent
    can sort one window (all samples) at a time. Returns (chrom, {window: [parts]}, nrows)."""
    sample = sample_name_from_bam(bam)
    window_bp = max(1, int(window_bp))
    rows = []
    parts_by_win = defaultdict(list)

    def _flush():
        if not rows:
            return
        df = pd.DataFrame(rows)
        wins = df["start0"].to_numpy() // window_bp
        for w in np.unique(wins):
            w = int(w)
            p = f"{shard_path}.w{w}.{len(parts_by_win[w])}.pkl"
            df[wins == w].to_pickle(p)
            parts_by_win[w].append(p)
        rows.clear()

    f32 = np.float32
    total = 0
    n_mm_parse_fail = 0
    cand_arr = np.array(sorted(p for (_c, p) in chrom_lookup.keys()), dtype=np.int64)
    bamf = pysam.AlignmentFile(bam, "rb")
    # A candidate site's contig may be absent from THIS sample's BAM header (multi-sample runs where a
    # contig -- e.g. a viral/alt chrom -- has reads in one sample but not another). fetch() on an
    # unknown contig raises ValueError and would abort the whole genotype stage; there is simply nothing
    # to extract for this (bam, chrom), so skip it.
    if chrom not in bamf.references or cand_arr.size == 0:
        bamf.close()
        return chrom, {}, 0
    call_prob1 = float(str(f32(1.0)))
    for read in bamf.fetch(chrom):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        r_end = read.reference_end
        if r_end is None:
            continue
        lo = int(np.searchsorted(cand_arr, read.reference_start))
        hi = int(np.searchsorted(cand_arr, r_end))
        if lo == hi:
            continue          # no candidate site inside this read's span -> nothing to emit
        # Parse MM FIRST (not read.modified_bases): a read that is entirely canonical for an implicit
        # group still has MM but an empty modified_bases, and its unmodified observations must be kept.
        mm_groups, mm_has_listed = parse_mm_groups(read)
        if not mm_groups:
            continue  # no MM/Mm tag -> no modification information at all
        mb = read.modified_bases or {}
        # BLOCKER: an htslib MM PARSE FAILURE (e.g. a low-complexity read missing a canonical base, whose
        # zero-delta group makes htslib reject the whole tag) returns modified_bases == {} even though the
        # tag declares LISTED positions. The implicit pass would then read "no listed calls" as "entirely
        # canonical" and write confident unmodified calls from unparsed data. Detect the mismatch and SKIP
        # the read (count it) rather than fabricate calls -- these were missing data before the B4 fix, and
        # must stay missing, not become wrong.
        if mm_has_listed and not mb:
            n_mm_parse_fail += 1
            continue
        seq = read.query_sequence
        if seq is None:
            continue
        qs, rsa = _aligned_qr(read)
        if qs is None:
            continue
        cand = cand_arr[lo:hi]
        pi = np.searchsorted(rsa, cand)
        ok = pi < rsa.size
        pi = pi[ok]; cand = cand[ok]
        hit = rsa[pi] == cand
        if not hit.any():
            continue
        cand_q = dict(zip(qs[pi[hit]].tolist(), cand[hit].tolist()))   # query pos -> ref pos (ascending)
        qname = read.query_name
        ref_strand = "-" if read.is_reverse else "+"
        zt, zg, zn, zm = _read_tags(read)
        nseq = len(seq)
        pos_mods = {}
        pos_base = {}
        for (base, _mstrand, mod_code), calls in mb.items():
            code = str(mod_code)
            for read_pos, ml in calls:
                if read_pos in cand_q:
                    d = pos_mods.get(read_pos)
                    if d is None:
                        d = pos_mods[read_pos] = {}
                    d[code] = ml
                    pos_base[read_pos] = base
        for read_pos, start0 in cand_q.items():
            sites = chrom_lookup.get((chrom, start0))
            if not sites:
                continue
            mods = pos_mods.get(read_pos)
            if mods is not None:
                # ---- LISTED call at a candidate site ----
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
                # Round-trip through float32's short repr so the emitted double matches what the
                # modkit path writes+parses (modkit prints float32; build parses it back to a double).
                if canon >= best_prob:
                    call_code = "-"
                    call_prob = float(str(f32(canon)))
                else:
                    call_code = best_code
                    call_prob = float(str(f32(best_prob)))
                base = str(pos_base[read_pos])
                for site in sites:
                    site_strand = str(site.get("strand", ""))
                    if site_strand and ref_strand and ref_strand not in {".", "?"} and site_strand != ref_strand:
                        continue
                    target_mod = str(site["mod_code"])
                    # BLOCKER-5: emit a row for target_mod ONLY if the read actually ASSESSED target_mod at
                    # this position -- either it LISTED target_mod here (target_mod in mods), or it declared an
                    # IMPLICIT group for target_mod on this base (so an unlisted position is a real canonical
                    # observation, mm_groups[(base,target_mod)] is True). A read that declared only A+a. must
                    # NOT get a fabricated 17596 "canonical" row at an A it never assessed for inosine.
                    if target_mod not in mods and not mm_groups.get((base.upper(), target_mod), False):
                        continue
                    if call_code == target_mod:
                        state_detail = "modified"
                        target_modified = 1
                    elif call_code == "-":
                        state_detail = "canonical"
                        target_modified = 0
                    else:
                        state_detail = "other_mod"
                        target_modified = 0
                    rows.append({
                        "sample": sample,
                        "qname": qname,
                        "mod_site_id": site["mod_site_id"],
                        "chrom": chrom,
                        "start0": start0,
                        "end0": safe_int(site.get("end0", start0 + 1), default=start0 + 1),
                        "strand": site_strand or ref_strand,
                        "target_mod_code": target_mod,
                        "call_code": call_code,
                        "state_detail": state_detail,
                        "target_modified": target_modified,
                        "call_prob": call_prob,
                        "canonical_base": base,
                        "modified_primary_base": base,
                        "fail": False,
                        "within_alignment": True,
                        "gene_id": str(site.get("gene_id", "")),
                        "gene_name": str(site.get("gene_name", "")),
                        "metagene_index": str(site.get("metagene_index", "")),
                        "ZT": zt, "ZG": zg, "ZN": zn, "ZM": zm,
                    })
                    total += 1
                    if len(rows) >= chunk_rows:
                        _flush()
            else:
                # ---- IMPLICIT-CANONICAL call (BLOCKER-4) ----
                # A candidate site the read covers at an UNLISTED position of an IMPLICIT MM group is a real
                # unmodified observation that modkit emits. read_pos indexes query_sequence and mb `base` is
                # transcript-oriented, so for a reverse read the stored base is complemented.
                qb = seq[read_pos].upper() if read_pos < nseq else ""
                if not qb:
                    continue
                tb = _COMP.get(qb, qb) if read.is_reverse else qb   # transcript-oriented read base
                for site in sites:
                    site_strand = str(site.get("strand", ""))
                    if site_strand and ref_strand and ref_strand not in {".", "?"} and site_strand != ref_strand:
                        continue
                    target_mod = str(site["mod_code"])
                    # MAJOR-1: look the group up by (read-base, code). Emit only if the read has an IMPLICIT
                    # group for this mod ON this read's actual base; a mismatched base / explicit group /
                    # unassessed mod all fall through with no call.
                    if not mm_groups.get((tb, target_mod), False):
                        continue
                    rows.append({
                        "sample": sample,
                        "qname": qname,
                        "mod_site_id": site["mod_site_id"],
                        "chrom": chrom,
                        "start0": start0,
                        "end0": safe_int(site.get("end0", start0 + 1), default=start0 + 1),
                        "strand": site_strand or ref_strand,
                        "target_mod_code": target_mod,
                        "call_code": "-",
                        "state_detail": "canonical",
                        "target_modified": 0,
                        "call_prob": call_prob1,
                        "canonical_base": tb,
                        "modified_primary_base": tb,
                        "fail": False,
                        "within_alignment": True,
                        "gene_id": str(site.get("gene_id", "")),
                        "gene_name": str(site.get("gene_name", "")),
                        "metagene_index": str(site.get("metagene_index", "")),
                        "ZT": zt, "ZG": zg, "ZN": zn, "ZM": zm,
                    })
                    total += 1
                    if len(rows) >= chunk_rows:
                        _flush()
    _flush()
    bamf.close()
    if verbose:
        print(f"[info] pysam extract done: {sample} {chrom} rows={total}", file=sys.stderr, flush=True)
    if n_mm_parse_fail:
        print(f"[warn] pysam extract: {sample} {chrom}: skipped {n_mm_parse_fail} read(s) whose MM tag "
              f"declared listed calls but htslib returned none (unparsable MM, e.g. low-complexity reads "
              f"missing a canonical base) -- NOT called canonical", file=sys.stderr, flush=True)
    return chrom, dict(parts_by_win), total


def _meta_by_zt(summary_tsv):
    """zt_label -> the assignment columns the read-assignment join used to supply (same names, and the
    integer columns as float64 so the CSV renders exactly as the old left-join did)."""
    summ = robust_load_summary(summary_tsv) if summary_tsv else pd.DataFrame()
    if summ.empty or "zt_label" not in summ.columns:
        return None
    keep = [c for c in ["zt_label", "gtf_gene_id", "gtf_gene_name", "gene_index", "transcript_index",
                        "metagene_index", "classification"] if c in summ.columns]
    meta = summ[keep].drop_duplicates("zt_label").rename(columns={
        "zt_label": "ZT", "gtf_gene_id": "assignment_gene_id", "gtf_gene_name": "assignment_gene_name",
        "metagene_index": "assignment_metagene_index"})
    for c in ("gene_index", "transcript_index", "assignment_metagene_index"):
        if c in meta.columns and pd.api.types.is_integer_dtype(meta[c]):
            meta[c] = meta[c].astype("float64")
    return meta.set_index("ZT")


def _empty_output(out_tsv):
    os.makedirs(os.path.dirname(out_tsv) or ".", exist_ok=True)
    pd.DataFrame(columns=OUTPUT_COLUMNS).to_csv(out_tsv, sep="\t", index=False)


def main():
    args = parse_args()
    cand = pd.read_csv(args.candidate_sites_tsv, sep="\t", low_memory=False)
    if cand.empty:
        _empty_output(args.out_tsv)
        return

    lookup = {}
    for row in cand.to_dict("records"):
        key = (str(row["chrom"]), int(row["start0"]))
        lookup.setdefault(key, []).append(row)

    shard_dir = tempfile.mkdtemp(prefix="molmod_shards.", dir=os.path.dirname(args.out_tsv) or ".")
    try:
        if args.pre_extracted:
            # Pre-extracted mode: `modkit extract calls` was run ONCE per subset BAM in a separate
            # per-sample sbatch (across nodes). Parse those TSVs here -- identical per-site logic,
            # just sourced from disk instead of re-running modkit windowed on one node.
            results = []
            for entry in args.pre_extracted:
                if "=" not in entry:
                    raise SystemExit(f"--pre-extracted entry must be SAMPLE=path, got: {entry}")
                sample, path = entry.split("=", 1)
                results.extend(parse_extracted_calls(
                    sample, path, lookup, shard_dir, args.chunk_rows, args.verbose))
        elif args.pysam:
            # pysam streaming backend: one task per (BAM x chromosome). Each streams reads one at a
            # time (peak RSS ~100MB, never OOMs -- chr15 included) and applies the identical per-site
            # expansion. No windowing / interval-size / reference / modkit needed.
            chrom_lookups = defaultdict(dict)
            for (c, p), sites in lookup.items():
                chrom_lookups[c][(c, p)] = sites
            chroms = sorted(chrom_lookups)
            if not chroms:
                _empty_output(args.out_tsv)
                return
            n_tasks = len(args.bams) * len(chroms)
            jobs = max(1, min(int(args.jobs), n_tasks))
            task_args = []
            for bam in args.bams:
                sample = sample_name_from_bam(bam)
                for chrom in chroms:
                    shard_path = os.path.join(shard_dir, f"{sample}.{chrom}.pkl")
                    task_args.append((bam, chrom_lookups[chrom], chrom, shard_path,
                                      args.chunk_rows, args.verbose, args.window_bp))
            if jobs == 1:
                results = [extract_rows_pysam(*item) for item in task_args]
            else:
                results = run_process_jobs(
                    extract_rows_pysam, task_args, jobs,
                    verbose=args.verbose, label="build_molecule_mod_table[pysam]",
                )
        else:
            # B: shard per (BAM x candidate-site window) so heavy chromosomes split into many balanced
            # tasks instead of one monolithic per-chrom extract that serializes the tail.
            windows = list(iter_site_windows(lookup, args.window_bp))
            if not windows:
                _empty_output(args.out_tsv)
                return
            n_tasks = len(args.bams) * len(windows)
            jobs = max(1, min(int(args.jobs), n_tasks))
            threads_per_job = max(1, int(args.threads) // jobs)
            task_args = []
            for bam in args.bams:
                sample = sample_name_from_bam(bam)
                for wi, (chrom, region, wlookup) in enumerate(windows):
                    shard_path = os.path.join(shard_dir, f"{sample}.{chrom}.{wi}.pkl")
                    task_args.append((bam, args.candidate_bed, args.reference_fa, args.modkit_bin,
                                      threads_per_job, wlookup, region, shard_path, args.interval_size,
                                      args.chunk_rows, args.verbose))

            if jobs == 1:
                results = [extract_rows_from_bam(*item) for item in task_args]
            else:
                results = run_process_jobs(
                    extract_rows_from_bam, task_args, jobs,
                    verbose=args.verbose, label="build_molecule_mod_table",
                )

        # Group shard files by chromosome and window. total_rows tells us whether anything survived.
        shards_by_chrom = defaultdict(lambda: defaultdict(list))   # chrom -> window -> [parts]
        total_rows = 0
        for chrom, parts, nrows in results:
            total_rows += int(nrows or 0)
            if isinstance(parts, dict):
                for w, pl in parts.items():
                    shards_by_chrom[chrom][int(w)].extend(pl)
            else:
                shards_by_chrom[chrom][0].extend(parts or [])

        if total_rows == 0:
            _empty_output(args.out_tsv)
            return

        # The pysam backend carries each read's ZT/ZG/ZN/ZM on the row, so the (genome-scale)
        # read-assignment table is never loaded; the per-fragmentform metadata comes from the small
        # classification summary instead. The modkit / pre-extracted paths keep the old join.
        use_tags = bool(args.pysam) and not args.pre_extracted and bool(args.summary_tsv)
        assignments = None
        meta = None
        if use_tags:
            meta = _meta_by_zt(args.summary_tsv)
        else:
            # Load read assignments ONCE (kept out of the fork so the table never gets COW-copied across
            # the worker pool). Index by (sample, qname) so the per-chrom join is an index lookup.
            assignments = load_read_assignments(args.read_assignments)
            keep_assign_cols = [c for c in [
                "sample", "qname", "ZT", "ZG", "ZN", "ZM", "assigned", "gene_id", "gene_name",
                "gene_index", "transcript_index", "metagene_index", "classification"
            ] if c in assignments.columns]
            assignments = assignments[keep_assign_cols].drop_duplicates(["sample", "qname"])
            assignments = assignments.rename(columns={
                col: f"assignment_{col}"
                for col in ["gene_id", "gene_name", "metagene_index"]
                if col in assignments.columns
            })
            assignments = assignments.set_index(["sample", "qname"])

        os.makedirs(os.path.dirname(args.out_tsv) or ".", exist_ok=True)
        tmp_out = args.out_tsv + ".tmp"
        wrote_header = False
        blocks = []
        with open(tmp_out, "w") as out_fh:
            # Stream chromosome-by-chromosome, window-by-window, in sorted order. Each (chrom, window)
            # block's rows are joined + sorted in isolation and appended; because (chrom, start0) lead the
            # sort key and windows partition start0, appending per-window-sorted blocks reproduces the
            # global sort byte-for-byte -- while holding one window of all samples in memory, never a
            # whole chromosome (chr1 of a 31-sample cohort is tens of GB of rows).
            for chrom in sorted(shards_by_chrom):
                off_chrom = out_fh.tell()
                for w in sorted(shards_by_chrom[chrom]):
                    parts = [pd.read_pickle(p) for p in shards_by_chrom[chrom][w]]
                    df = pd.concat(parts, ignore_index=True)
                    if use_tags:
                        df["assigned"] = normalize_string_series(df["ZT"]).ne("")
                        if meta is not None:
                            df = df.join(meta, on="ZT")
                        for c in ("assignment_gene_id", "assignment_gene_name", "gene_index",
                                  "transcript_index", "assignment_metagene_index", "classification"):
                            if c not in df.columns:
                                df[c] = np.nan
                    else:
                        df = df.join(assignments, on=["sample", "qname"], how="left")
                    for col in ["gene_id", "gene_name", "metagene_index"]:
                        assign_col = f"assignment_{col}"
                        if col in df.columns and assign_col in df.columns:
                            primary = normalize_string_series(df[col])
                            fallback = normalize_string_series(df[assign_col])
                            df[col] = primary.where(primary.ne(""), fallback)
                    df["usable"] = (~df["fail"].fillna(True)) & df["within_alignment"].fillna(False)
                    df = df.sort_values(["chrom", "start0", "mod_site_id", "sample", "qname"]).reset_index(drop=True)
                    if use_tags:
                        df = df[[c for c in FINAL_COLUMNS if c in df.columns]]
                    df.to_csv(out_fh, sep="\t", index=False, header=not wrote_header)
                    wrote_header = True
                    del df, parts
                out_fh.flush()
                blocks.append((chrom, off_chrom, out_fh.tell() - off_chrom))
        os.replace(tmp_out, args.out_tsv)
        write_chrom_index(args.out_tsv, blocks)
    finally:
        shutil.rmtree(shard_dir, ignore_errors=True)


if __name__ == "__main__":
    main()
