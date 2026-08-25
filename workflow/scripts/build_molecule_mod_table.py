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

from genotype_utils import load_read_assignments, normalize_string_series, run_process_jobs, sample_name_from_bam, safe_float, safe_int


OUTPUT_COLUMNS = [
    "sample", "qname", "mod_site_id", "chrom", "start0", "end0", "strand",
    "target_mod_code", "call_code", "state_detail", "target_modified",
    "call_prob", "canonical_base", "modified_primary_base", "fail",
    "within_alignment", "gene_id", "gene_name", "metagene_index",
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
    """Parse the MM/Mm tag into {mod_code: (canonical_base, is_implicit)}.

    BLOCKER-4: a modBAM in IMPLICIT MM mode (flag '.' or absent) declares that every canonical base of
    a group NOT listed in the deltas is an implicitly-canonical (unmodified) observation. read.modified_bases
    returns ONLY the listed positions, so the pysam backend silently dropped those unmodified calls and
    inflated every modified fraction. Explicit mode ('?') means unlisted positions carry NO call. This
    parser recovers the flag so the extractor can emit the implicit-canonical rows modkit already emits."""
    mm = None
    for tag in ("MM", "Mm"):
        try:
            mm = read.get_tag(tag)
            break
        except KeyError:
            continue
    if not mm:
        return {}
    out = {}
    for entry in str(mm).split(";"):
        m = _MM_ENTRY_RE.match(entry.strip())
        if not m:
            continue
        base, _strand, mods, flag = m.groups()
        implicit = (flag != "?")   # '.' or absent -> implicit; '?' -> explicit (unlisted = no call)
        base = base.upper()
        for code in re.findall(r'[0-9]+|[a-z]', mods):   # 'mh' -> m,h ; '17802' -> 17802
            out[code] = (base, implicit)
    return out


def parse_args():
    ap = argparse.ArgumentParser(description="Build a per-read mod call table at candidate modulator sites.")
    ap.add_argument("--bams", nargs="+", required=True, help="Input BAMs with MM/ML tags")
    ap.add_argument("--candidate-sites-tsv", required=True, help="Candidate mod site TSV")
    ap.add_argument("--candidate-bed", required=True, help="Candidate mod BED for modkit include-bed")
    ap.add_argument("--read-assignments", required=True, help="Read assignment TSV")
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


def extract_rows_pysam(bam, chrom_lookup, chrom, shard_path, chunk_rows, verbose=False):
    """Stream one chromosome of a modBAM with pysam and emit the same per-(read, candidate site) rows
    the modkit path produces -- one read at a time, so peak RSS is ~100MB regardless of BAM/chrom size
    (never OOMs). Reproduces modkit `extract calls --no-filtering --mapped-only` semantics -- INCLUDING
    implicit-canonical calls (unlisted canonical bases of an implicit MM group), which the previous
    version dropped, inflating every modified fraction on real ONT data (BLOCKER-4). Namely:
    call_prob=(ML+0.5)/256 (float32), canonical=1-sum(mod_probs), call_code=argmax, strand-aware.
    Flushes to numbered pickle parts every chunk_rows. Returns (chrom, [parts], nrows)."""
    sample = sample_name_from_bam(bam)
    rows = []
    parts = []

    def _flush():
        if not rows:
            return
        p = f"{shard_path}.{len(parts)}.pkl"
        pd.DataFrame(rows).to_pickle(p)
        parts.append(p)
        rows.clear()

    f32 = np.float32
    total = 0
    bamf = pysam.AlignmentFile(bam, "rb")
    # A candidate site's contig may be absent from THIS sample's BAM header (multi-sample runs where a
    # contig -- e.g. a viral/alt chrom -- has reads in one sample but not another). fetch() on an
    # unknown contig raises ValueError and would abort the whole genotype stage; there is simply nothing
    # to extract for this (bam, chrom), so skip it.
    if chrom not in bamf.references:
        bamf.close()
        return chrom, [], 0
    for read in bamf.fetch(chrom):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        # Parse MM FIRST (not read.modified_bases): a read that is entirely canonical for an implicit
        # group still has MM but an empty modified_bases, and its unmodified observations must be kept.
        mm_groups = parse_mm_groups(read)
        if not mm_groups:
            continue  # no MM/Mm tag -> no modification information at all
        mb = read.modified_bases or {}
        seq = read.query_sequence
        if seq is None:
            continue
        qname = read.query_name
        ref_strand = "-" if read.is_reverse else "+"
        refpos = read.get_reference_positions(full_length=True)  # per query_sequence position -> ref pos
        nrp = len(refpos)
        pos_mods = {}
        pos_base = {}
        for (base, _mstrand, mod_code), calls in mb.items():
            code = str(mod_code)
            for read_pos, ml in calls:
                d = pos_mods.get(read_pos)
                if d is None:
                    d = pos_mods[read_pos] = {}
                d[code] = ml
                pos_base[read_pos] = base
        emitted = set()  # query positions handled as a LISTED call (so the implicit pass skips them)
        for read_pos, mods in pos_mods.items():
            if read_pos >= nrp:
                continue
            start0 = refpos[read_pos]
            if start0 is None:
                continue
            sites = chrom_lookup.get((chrom, start0))
            if not sites:
                continue
            emitted.add(read_pos)
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
                # base-mismatch guard (parity with the modkit path): the read's base at this site cannot
                # carry target_mod -> it is a variant, not a usable unmodified observation. Skip it.
                if _base_mismatch(target_mod, base):
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
                })
                total += 1
                if len(rows) >= chunk_rows:
                    _flush()

        # ---- IMPLICIT-CANONICAL pass (BLOCKER-4) ----
        # A candidate site the read covers at an UNLISTED position of an IMPLICIT MM group is a real
        # unmodified observation that modkit emits and the old pysam backend dropped. Iterate the read's
        # aligned positions; for each candidate site not already emitted as a listed call, if the read's
        # transcript-oriented base is the mod's canonical base and that mod's MM group is implicit, emit a
        # canonical row. read_pos indexes query_sequence and mb `base` is transcript-oriented, so for a
        # reverse read the stored base is complemented to get the transcript base.
        for q in range(nrp):
            if q in emitted:
                continue
            start0 = refpos[q]
            if start0 is None:
                continue
            sites = chrom_lookup.get((chrom, start0))
            if not sites:
                continue
            qb = seq[q].upper() if q < len(seq) else ""
            if not qb:
                continue
            tb = _COMP.get(qb, qb) if read.is_reverse else qb   # transcript-oriented read base
            call_prob1 = float(str(f32(1.0)))
            for site in sites:
                site_strand = str(site.get("strand", ""))
                if site_strand and ref_strand and ref_strand not in {".", "?"} and site_strand != ref_strand:
                    continue
                target_mod = str(site["mod_code"])
                grp = mm_groups.get(target_mod)
                if grp is None or not grp[1]:
                    continue  # this mod was not assessed on the read, OR is explicit ('?') -> no call
                if _base_mismatch(target_mod, tb):
                    continue  # variant at the modified base -> not a canonical observation
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
                })
                total += 1
                if len(rows) >= chunk_rows:
                    _flush()
    _flush()
    if verbose:
        print(f"[info] pysam extract done: {sample} {chrom} rows={total}", file=sys.stderr, flush=True)
    return chrom, parts, total


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
                                      args.chunk_rows, args.verbose))
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

        # A: group shard files by chromosome. total_rows tells us whether anything survived.
        shards_by_chrom = defaultdict(list)
        total_rows = 0
        for chrom, parts, nrows in results:
            total_rows += int(nrows or 0)
            for p in (parts or []):
                shards_by_chrom[chrom].append(p)

        if total_rows == 0:
            _empty_output(args.out_tsv)
            return

        # Load read assignments ONCE (kept out of the fork so the ~22GB/70GB-in-pandas table never
        # gets COW-copied across the worker pool). Index by (sample, qname) so the per-chrom join is
        # an index lookup, not a re-hash of the whole table each chromosome.
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
        with open(tmp_out, "w") as out_fh:
            # A: stream chromosome-by-chromosome in sorted order. Each chrom's rows are joined + sorted
            # in isolation and appended; because `chrom` is the primary sort key and the full sort key
            # (chrom,start0,mod_site_id,sample,qname) is unique per (site,read,sample), appending
            # per-chrom-sorted blocks reproduces the global sort byte-for-byte -- without ever holding
            # every chromosome's rows (or a whole-table sort) in memory at once.
            for chrom in sorted(shards_by_chrom):
                parts = [pd.read_pickle(p) for p in shards_by_chrom[chrom]]
                df = pd.concat(parts, ignore_index=True)
                df = df.join(assignments, on=["sample", "qname"], how="left")
                for col in ["gene_id", "gene_name", "metagene_index"]:
                    assign_col = f"assignment_{col}"
                    if col in df.columns and assign_col in df.columns:
                        primary = normalize_string_series(df[col])
                        fallback = normalize_string_series(df[assign_col])
                        df[col] = primary.where(primary.ne(""), fallback)
                df["usable"] = (~df["fail"].fillna(True)) & df["within_alignment"].fillna(False)
                df = df.sort_values(["chrom", "start0", "mod_site_id", "sample", "qname"]).reset_index(drop=True)
                df.to_csv(out_fh, sep="\t", index=False, header=not wrote_header)
                wrote_header = True
                del df, parts
        os.replace(tmp_out, args.out_tsv)
    finally:
        shutil.rmtree(shard_dir, ignore_errors=True)


if __name__ == "__main__":
    main()
