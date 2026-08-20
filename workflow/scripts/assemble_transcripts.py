#!/usr/bin/env python3
"""
Read-backed TES/APA assembler (v9, regionized)

Key changes vs. prior v9 drafts
-------------------------------
- Fast, deterministic regionization by streaming each BAM once per chromosome
  to find all-sample zero-coverage gaps over coarse bins (--cov-bin-bp).
- Optional deterministic random breakpoints injected INSIDE all-zero gaps
  (--max-breakpoints-per-chrom, --rand-seed).
- Cores = intervals between gaps; each core is processed independently
  (parallel with --threads). Reads are fetched with ±pad, but an isoform is
  KEPT only if its TES falls INSIDE the core (prevents duplicates).
- Deterministic ordering preserved (stable sorts, seeded RNG only for gap injection).

Outputs
-------
- <out>.gtf
- <prefix>_metrics.tsv
- <prefix>_classification_summary.tsv
- <prefix>_tx_counts.tsv, <prefix>_tx_counts.pca.png
- <prefix>_per_sample_stats.tsv
- zt_tagged/*.zt_tagged.bam  (optional)
- zt_bams/*.bam + modkit_manifest.tsv (optional)

Typical use
-----------
python assemble_transcripts.py \
  --dir <bam_dir> --glob "*.bam" \
  --gtf hg38.ncbiRefSeq.gtf \
  --out-gtf results/assemble/readbacked_annot.gtf \
  --out-prefix results/assemble/readbacked_annot \
  --threads 64 \
  --primary-only --min-mapq 10 --min-introns-read 1 \
  --min-reads 40 --min-frac 0.00 --min-introns 1 \
  --min-polya-length 12 --min-polya-purity 0.5 --polya-support-frac 0.5 \
  --tes-match-tol 25 --exact-tes-tol 10 \
  --write-zt-tagged-sample-bams \
  --min-reads-per-sample-for-mod 5 \
  --min-total-reads-for-mod 20
"""

import argparse, os, sys, glob, re, gzip, math, random
from collections import defaultdict, Counter
from concurrent.futures import ProcessPoolExecutor, as_completed

import pysam


# --- HPC page-cache control -------------------------------------------------
# Sequentially streaming very large BAMs populates the OS page cache, which on
# cgroup-v2 HPC nodes is charged to the job's memory and can grow to the size of
# the data read (hundreds of GB) -> node-level OOM, even though the process's
# real RSS stays tiny. Periodically advising the kernel to drop already-read
# pages keeps that cache flat and does not change any results.
def _drop_page_cache(path, upto_bytes=0):
    """Best-effort POSIX_FADV_DONTNEED over [0, upto_bytes) of `path` (0 => whole file)."""
    if not path:
        return
    try:
        if isinstance(path, bytes):
            path = path.decode()
        fd = os.open(path, os.O_RDONLY)
        try:
            os.posix_fadvise(fd, 0, int(upto_bytes), os.POSIX_FADV_DONTNEED)
        finally:
            os.close(fd)
    except Exception:
        pass


def _bgzf_coffset(fh):
    """Compressed file offset (bytes) consumed so far by a pysam BGZF reader, else 0."""
    try:
        return max(0, int(fh.tell()) >> 16)
    except Exception:
        return 0

# --------------------------- utils ---------------------------

def median(vals):
    v = [x for x in vals if x is not None]
    if not v: return 0
    try:
        import statistics
        return statistics.median(v)
    except Exception:
        v.sort(); n=len(v)
        return v[n//2] if n%2 else 0.5*(v[n//2-1]+v[n//2])

def get_tx_strand(aln):
    # Genomic strand ONLY from FLAG (0x10).
    return "-" if aln.is_reverse else "+"

def tes_pos1(aln, tx):
    return aln.reference_end if tx=="+" else aln.reference_start+1

def intron_chain_1based(aln):
    chain = []
    ref = aln.reference_start
    for op, ln in (aln.cigartuples or []):
        if op == 3:  # N
            donor = ref + 1
            acceptor = ref + ln
            chain.append((donor, acceptor))
            ref += ln
        elif op in (0,2,7,8):  # M/D/=/X
            ref += ln
    return tuple(chain)

def chain_tx_order(chain, strand):
    return chain if strand=="+" else tuple(reversed(chain))

def exon_blocks_from_aln(aln):
    exons = []
    ref = aln.reference_start
    cur_start = ref
    for op, ln in (aln.cigartuples or []):
        if op == 3:
            exons.append((cur_start+1, ref))
            ref += ln
            cur_start = ref
        elif op in (0,2,7,8):
            ref += ln
    exons.append((cur_start+1, ref))
    return exons

def softclip3p_len_and_seq(aln, tx):
    ct = aln.cigartuples or []
    if tx == "+":
        if ct and ct[-1][0] == 4:
            L = ct[-1][1]
            seq = (aln.query_sequence or "")[-L:] if L>0 else ""
            return L, seq
        return 0, ""
    else:
        if ct and ct[0][0] == 4:
            L = ct[0][1]
            seq = (aln.query_sequence or "")[:L] if L>0 else ""
            return L, seq
        return 0, ""

def polya_purity(seq, tx):
    # + expects polyA; - expects polyT
    if not seq: return 0.0
    s = seq.upper()
    base = "A" if tx == "+" else "T"
    return s.count(base) / len(s)

def query_len(aln):
    L = aln.query_length
    if L is not None:
        return int(L)
    seq = aln.query_sequence
    return int(len(seq)) if seq else 0

def cluster_positions(sorted_positions, window):
    clusters = []
    cur = [sorted_positions[0]]
    for p in sorted_positions[1:]:
        if p - cur[-1] <= window:
            cur.append(p)
        else:
            clusters.append(cur); cur = [p]
    clusters.append(cur)
    out = []
    for cl in clusters:
        c = Counter(cl)
        rep = sorted(c.items(), key=lambda x: (x[1], x[0]))[-1][0]
        out.append({"positions": cl, "rep": rep, "count": sum(c.values())})
    return out

def is_suffix(longer, shorter):
    if len(shorter) > len(longer): return False
    if len(shorter) == 0: return True
    return longer[-len(shorter):] == shorter

def chain_to_str(chain):
    return "." if not chain else ";".join(f"{d}-{a}" for d,a in chain)

def merge_intervals(intervals):
    if not intervals:
        return []
    merged = [list(sorted(intervals)[0])]
    for s, e in sorted(intervals)[1:]:
        if s <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    return [(s, e) for s, e in merged]

def first_overlap_bp(exons_a, exons_b):
    return exon_overlap_len(exons_a, exons_b)

def smallest_unused_positive(used_vals):
    x = 1
    while x in used_vals:
        x += 1
    return x

def compute_chain_features(exact_counts):
    observed = list(exact_counts.keys())
    features = {}
    for chain in observed:
        shorter_suffixes = [
            other for other in observed
            if other != chain and len(other) < len(chain) and is_suffix(chain, other)
        ]
        max_suffix_len = max((len(other) for other in shorter_suffixes), default=0)
        unique_prefix = tuple(chain[: max(0, len(chain) - max_suffix_len)])
        reachable_count = sum(
            exact_counts[other] for other in observed
            if is_suffix(chain, other)
        )
        exact_count = exact_counts[chain]
        anchor_frac = (exact_count / reachable_count) if reachable_count else 0.0
        features[chain] = dict(
            exact_count=exact_count,
            reachable_count=reachable_count,
            has_shorter_suffix=bool(shorter_suffixes),
            unique_prefix=unique_prefix,
            unique_prefix_len=len(unique_prefix),
            distal_unique_5p_junction=(unique_prefix[0] if unique_prefix else None),
            anchor_reads=exact_count,
            anchor_frac=anchor_frac,
        )
    return features

def absorb_allowed_for_chain(feature, iso_params):
    if feature["exact_count"] < iso_params["min_exact_canonical_reads"]:
        return False
    if not feature["has_shorter_suffix"]:
        return True
    if feature["anchor_reads"] < iso_params["min_distal_anchor_reads"]:
        return False
    if feature["anchor_frac"] < iso_params["min_distal_anchor_frac"]:
        return False
    return True

def canonical_order(exact_counts, features):
    chains = list(exact_counts.keys())
    return sorted(
        chains,
        key=lambda ch: (
            features[ch]["exact_count"],
            features[ch]["anchor_frac"],
            features[ch]["reachable_count"],
            len(ch),
        ),
        reverse=True,
    )

# --------------------------- GTF parsing ---------------------------

_attr_re = re.compile(r'(\S+)\s+"([^"]+)"')

def parse_gtf_attrs(attr_field):
    d = {}
    for m in _attr_re.finditer(attr_field):
        d[m.group(1)] = m.group(2)
    if "gene_name" not in d and "gene_id" in d: d["gene_name"] = d["gene_id"]
    return d

class GTFTranscript:
    __slots__ = ("chrom","strand","start","end","tes_1based","exons","gene_id","gene_name","transcript_id","chain_tx")
    def __init__(self, chrom, strand, start, end, attrs):
        self.chrom=chrom; self.strand=strand; self.start=start; self.end=end
        self.tes_1based = end if strand=="+" else start
        self.exons = []
        self.gene_id = attrs.get("gene_id","NA")
        self.gene_name = attrs.get("gene_name", self.gene_id)
        self.transcript_id = attrs.get("transcript_id", "NA")
        self.chain_tx = tuple()

def intron_chain_from_exons(exons, strand):
    if not exons or len(exons)<2: return tuple()
    exons_sorted = sorted(exons)
    introns = []
    for i in range(len(exons_sorted)-1):
        e1 = exons_sorted[i][1]
        s2 = exons_sorted[i+1][0]
        introns.append((e1+1, s2-1))
    chain = tuple(introns)
    return chain if strand=="+" else tuple(reversed(chain))

def load_gtf(gtf_path):
    txs = {}
    opener = gzip.open if gtf_path.endswith('.gz') else open
    mode = 'rt' if gtf_path.endswith('.gz') else 'r'
    with opener(gtf_path, mode) as f:
        for line in f:
            if not line or line.startswith("#"): continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9: continue
            chrom, source, feature, start, end, score, strand, frame, attrs = parts
            start, end = int(start), int(end)
            a = parse_gtf_attrs(attrs)
            if feature == "transcript":
                tid = a.get("transcript_id")
                if not tid: continue
                txs[tid] = GTFTranscript(chrom, strand, start, end, a)
            elif feature == "exon":
                tid = a.get("transcript_id")
                if tid and tid in txs:
                    txs[tid].exons.append((start, end))
    for tx in txs.values():
        tx.chain_tx = intron_chain_from_exons(tx.exons, tx.strand)
    return list(txs.values())

def build_gtf_index(gtf_txs):
    by_cs = defaultdict(list)
    for tx in gtf_txs:
        by_cs[(tx.chrom, tx.strand)].append(tx)
    tes_sorted = {}
    for k, lst in by_cs.items():
        tes_sorted[k] = sorted([(t.tes_1based, t) for t in lst], key=lambda x: x[0])
    return by_cs, tes_sorted

# --------------------------- annotation ---------------------------

def exon_overlap_len(ex1, ex2):
    tot=0
    for s1,e1 in ex1:
        for s2,e2 in ex2:
            lo=max(s1,s2); hi=min(e1,e2)
            if hi>=lo: tot += (hi-lo+1)
    return tot

def annotate_isoform(iso, gtf_by_cs, tes_sorted_index, tes_match_tol, exact_tes_tol):
    chrom, strand, tes = iso["chrom"], iso["strand"], iso["tes"]
    key = (chrom, strand)
    cands = []
    if key in tes_sorted_index:
        arr = tes_sorted_index[key]
        lo, hi = 0, len(arr)
        left = 0
        # lower bound for tes - tol
        x = tes - tes_match_tol
        while lo < hi:
            mid = (lo+hi)//2
            if arr[mid][0] < x: lo = mid+1
            else: hi = mid
        left = lo
        i = left
        upper = tes + tes_match_tol
        while i < len(arr) and arr[i][0] <= upper:
            cands.append(arr[i][1]); i += 1
    match_source = "TES_TOL" if cands else "NONE"
    best = None
    if cands:
        def rank(tx):
            same_chain = 0 if tx.chain_tx == iso["chain_tx"] else 1
            return (abs(tx.tes_1based - tes), same_chain, -exon_overlap_len(iso["rep_exons"], tx.exons), -len(tx.chain_tx))
        best = min(cands, key=rank)
    else:
        overlaps = []
        for tx in gtf_by_cs.get(key, []):
            if exon_overlap_len(iso["rep_exons"], tx.exons) > 0:
                overlaps.append(tx)
        if overlaps:
            best = max(overlaps, key=lambda tx: exon_overlap_len(iso["rep_exons"], tx.exons))
            match_source = "OVERLAP"
    if best:
        tes_delta = abs(best.tes_1based - tes)
        same_chain = (best.chain_tx == iso["chain_tx"])
        if same_chain and tes_delta < exact_tes_tol:
            cls = "EXACT"
        elif same_chain:
            cls = "NOVEL_APA"
        else:
            cls = "NOVEL_CHAIN"
        gene_id, gene_name, tid = best.gene_id, best.gene_name, best.transcript_id
        ov_bp = exon_overlap_len(iso["rep_exons"], best.exons)
        return dict(classification=cls, gene_id=gene_id, gene_name=gene_name,
                    matched_tid=tid, gtf_tes=best.tes_1based, gtf_chain_tx=best.chain_tx,
                    tes_delta_bp=tes_delta, exon_overlap_bp=ov_bp, match_source=match_source)
    else:
        return dict(classification="NOVEL_LOCUS", gene_id="NA", gene_name="NA",
                    matched_tid="NA", gtf_tes="NA", gtf_chain_tx=tuple(), tes_delta_bp="NA", exon_overlap_bp=0, match_source=match_source)

# --------------------------- regionization ---------------------------

def _bins_for_chrom(chrom_len, bin_bp):
    n = max(1, math.ceil(chrom_len / bin_bp))
    starts = [i*bin_bp for i in range(n)]
    ends = [min(chrom_len, (i+1)*bin_bp) for i in range(n)]
    return list(zip(starts, ends))  # 0-based half-open

def _presence_bins_for_chrom(handles, chrom, chrom_len, bin_bp, status_every_reads=0):
    import numpy as np
    bins = _bins_for_chrom(chrom_len, bin_bp)
    nb = len(bins)
    union_presence = np.zeros(nb, dtype=np.bool_)
    for hi, fh in enumerate(handles):
        seen = 0
        bam_path = getattr(fh, "filename", None)
        try:
            for aln in fh.fetch(contig=chrom):
                seen += 1
                if aln.is_unmapped: continue
                s = aln.reference_start
                e = aln.reference_end
                if e <= 0: continue
                if s < 0: s = 0
                if e > chrom_len: e = chrom_len
                if e <= s: continue
                bs = s // bin_bp
                be = (e - 1) // bin_bp
                union_presence[bs:be+1] = True
                if seen % 1_000_000 == 0:
                    _drop_page_cache(bam_path, _bgzf_coffset(fh))
                if status_every_reads and (seen % status_every_reads == 0):
                    print(f"[REGIONIZE] {chrom} reads streamed (bam#{hi+1}): {seen:,}", file=sys.stderr)
        except KeyError:
            continue
        # Release the page cache consumed by this chromosome's scan of this BAM.
        _drop_page_cache(bam_path, _bgzf_coffset(fh))
    all_zero = ~union_presence
    return bins, all_zero

def _zero_runs_to_gaps(bins, all_zero, min_gap_bins):
    gaps = []
    i = 0
    n = len(all_zero)
    while i < n:
        if not all_zero[i]:
            i += 1; continue
        j = i
        while j < n and all_zero[j]:
            j += 1
        if (j - i) >= min_gap_bins:
            gaps.append((bins[i][0], bins[j-1][1]))  # [s,e)
        i = j
    return gaps

def _inject_random_breaks(chrom_len, gaps, bin_bp, max_breaks, seed=1337):
    if max_breaks <= 0: return []
    rng = random.Random(seed)
    # sample candidate bins uniformly; keep only those landing in a gap
    picked = []
    if not gaps: return picked
    # Flatten gaps into cumulative lengths for sampling proportional to gap size
    total_gap_bp = sum(e-s for s,e in gaps)
    if total_gap_bp <= 0: return picked
    # draw up to max_breaks distinct positions
    tried = 0
    seen = set()
    while len(picked) < max_breaks and tried < max_breaks * 20:
        tried += 1
        # choose a random point in the union of gaps
        r = rng.randrange(total_gap_bp)
        acc = 0
        for (s,e) in gaps:
            span = e - s
            if r < acc + span:
                pos = s + (r - acc)
                # snap to bin boundary for stability
                pos = (pos // bin_bp) * bin_bp
                if pos <= 0 or pos >= chrom_len: break
                if pos not in seen:
                    seen.add(pos)
                    picked.append((pos, pos))  # point break
                break
            acc += span
    # convert points to degenerate gaps so they act as cuts
    return [(p, p) for (p, _) in picked]

def _cores_from_gaps(chrom_len, gaps):
    # gaps are [s,e) disjoint (not guaranteed, but OK). Build complement intervals
    pts = [0]
    for s,e in gaps:
        s = max(0, min(chrom_len, s))
        e = max(0, min(chrom_len, e))
        if e <= s: continue
        pts.append(s); pts.append(e)
    pts.append(chrom_len)
    pts = sorted(set(pts))
    cores = []
    for i in range(len(pts)-1):
        s, e = pts[i], pts[i+1]
        if e > s:
            cores.append((s,e))
    # merge tiny slivers created by point breaks (zero-length removed above).
    return cores


def _regionize_chrom_core(handles, chrom, clen, cov_bin_bp, status_every, min_gap_bins,
                          max_breaks, rand_seed, c0, c1):
    """Regionize ONE chromosome given already-open BAM handles: stream its reads to build
    zero-coverage bins, derive gaps -> cores, clip to the requested subregion.

    Deterministic: _inject_random_breaks uses a LOCAL random.Random(seed) and
    _zero_runs_to_gaps/_cores_from_gaps are pure, so a chromosome's cores never depend on which
    worker runs it or in what order.
    """
    bins, all_zero = _presence_bins_for_chrom(
        handles, chrom, clen, cov_bin_bp, status_every_reads=status_every
    )
    gaps = _zero_runs_to_gaps(bins, all_zero, min_gap_bins)
    injected = _inject_random_breaks(clen, gaps, cov_bin_bp, max_breaks, seed=rand_seed)
    cores = _cores_from_gaps(clen, gaps + injected)
    kept = []
    for (s, e) in cores:
        if e <= c0 or s >= c1: continue
        ss = max(s, c0); ee = min(e, c1)
        if ee > ss:
            kept.append((chrom, ss, ee))
    return chrom, len(cores), kept


def _regionize_chrom_worker(task):
    """Picklable ProcessPool entry point: open this worker's own BAM handles, regionize one
    chromosome, close them. Result is identical to the serial shared-handle path."""
    (chrom, clen, bam_paths, cov_bin_bp, status_every, min_gap_bins,
     max_breaks, rand_seed, c0, c1) = task
    handles = [pysam.AlignmentFile(bp, "rb") for bp in bam_paths]
    try:
        return _regionize_chrom_core(handles, chrom, clen, cov_bin_bp, status_every,
                                     min_gap_bins, max_breaks, rand_seed, c0, c1)
    finally:
        for fh in handles:
            try: fh.close()
            except Exception: pass

def _write_zt_tagged_sample(task):
    """Write one sample's ZT/ZN-tagged BAM (picklable ProcessPool entry point). Each sample re-reads
    its own original BAM and writes its own tagged BAM, so samples are independent -> parallelizable.
    BGZF read/write use compression threads (`threads=`), which is the main per-sample speedup."""
    (in_path, out_path, sample_assign, primary_only, io_threads) = task
    io_threads = max(1, int(io_threads))
    _seen_w = 0
    with pysam.AlignmentFile(in_path, "rb", threads=io_threads) as inp, \
         pysam.AlignmentFile(out_path, "wb", header=inp.header, threads=io_threads) as outw:
        for aln in inp.fetch():
            _seen_w += 1
            if _seen_w % 1_000_000 == 0:
                _drop_page_cache(in_path, _bgzf_coffset(inp))
            if aln.is_unmapped:
                outw.write(aln); continue
            if primary_only and (aln.is_secondary or aln.is_supplementary):
                continue
            tup = sample_assign.get(aln.query_name) if sample_assign else None
            if tup:
                zt, zg, zn, zm = tup
                try:
                    aln.set_tag("ZT", zt, value_type="Z", replace=True)
                except TypeError:
                    aln.set_tag("ZT", zt)
                for tag, val in (("ZG", int(zg)), ("ZN", int(zn)), ("ZM", int(zm))):
                    try:
                        aln.set_tag(tag, val, value_type="i", replace=True)
                    except TypeError:
                        aln.set_tag(tag, val)
            outw.write(aln)
    try:
        pysam.index(out_path)
    except Exception:
        pass
    return out_path


# --------------------------- per-core processing ---------------------------

def _process_core(core_args):
    """
    Fetch reads for a (chrom, s, e) core with padding; cluster TES; assign truncations.
    Only keep isoforms whose TES is within [s,e] (core span) to avoid duplicates.
    Return list of isoform dicts (without annotation yet).
    """

    (chrom, s, e, pad_bp, bam_paths, filters, iso_params) = core_args
    apa_window = iso_params["apa_window"]

    # -----------------------------
    # Collect reads
    # -----------------------------
    mem = []
    for bam in bam_paths:
        sample = os.path.basename(bam).replace(".bam", "")
        if not os.path.exists(bam):
            continue

        with pysam.AlignmentFile(bam, "rb") as fh:
            s_fetch = max(0, s - pad_bp)
            e_fetch = e + pad_bp

            try:
                it = fh.fetch(contig=chrom, start=s_fetch, end=e_fetch)
            except ValueError:
                continue

            for aln in it:
                if aln.is_unmapped:
                    continue
                if filters["primary_only"] and (aln.is_secondary or aln.is_supplementary):
                    continue
                if aln.mapping_quality < filters["min_mapq"]:
                    continue

                tx = get_tx_strand(aln)

                sclen, tail = softclip3p_len_and_seq(aln, tx)
                if filters["require_softclip3p"] > 0 and sclen < filters["require_softclip3p"]:
                    continue

                purity = polya_purity(tail, tx) if sclen > 0 else 0.0

                chain = intron_chain_1based(aln)
                if len(chain) < filters["min_introns_read"]:
                    continue

                chain_tx = chain_tx_order(chain, tx)

                rd_chrom = fh.get_reference_name(aln.reference_id)
                if rd_chrom != chrom:
                    continue

                mem.append(dict(
                    chrom=chrom,
                    strand=tx,
                    tes=tes_pos1(aln, tx),
                    chain=chain,
                    chain_tx=chain_tx,
                    n_introns=len(chain),
                    exons=exon_blocks_from_aln(aln),
                    bam=os.path.basename(bam),
                    sample=sample,
                    qname=aln.query_name,
                    mapq=aln.mapping_quality,
                    read_length=query_len(aln),
                    sclen=sclen,
                    purity=purity
                ))

    if not mem:
        return []

    # -----------------------------
    # Group by chrom + strand
    # -----------------------------
    by_cs = defaultdict(list)
    for r in mem:
        by_cs[(r["chrom"], r["strand"])].append(r)

    isoforms = []

    # -----------------------------
    # Process each strand separately
    # -----------------------------
    for (c, strand), rlist in by_cs.items():

        positions = sorted(r["tes"] for r in rlist)
        if not positions:
            continue

        clusters = cluster_positions(positions, apa_window)

        for cl in clusters:

            rep_pos = cl["rep"]

            if not (s <= rep_pos <= e):
                continue

            # Assign every read whose TES falls in THIS cluster's span, not just those within
            # apa_window of the mode. cluster_positions() single-linkage-chains positions, so a
            # cluster can be wider than apa_window; selecting by |tes-rep_pos|<=apa_window would
            # silently drop the tail of such a cluster (those reads then match no other cluster
            # either, since inter-cluster gaps exceed the window) -- losing reads and, if the tail
            # is a real distinct 3' end, an entire APA isoform. Clusters are non-overlapping, so a
            # read's TES belongs to exactly one span.
            cl_lo, cl_hi = min(cl["positions"]), max(cl["positions"])
            members = [r for r in rlist if cl_lo <= r["tes"] <= cl_hi]
            if not members:
                continue

            # --------------------------------------------------------
            # TRUE 3′-ANCHORED SUFFIX COLLAPSE
            # --------------------------------------------------------

            chain_to_idxs = defaultdict(list)
            for i, m in enumerate(members):
                chain_to_idxs[tuple(m["chain_tx"])].append(i)

            exact_counts = Counter(tuple(m["chain_tx"]) for m in members)
            chain_features = compute_chain_features(exact_counts)
            for canon in chain_features:
                chain_features[canon]["absorb_allowed"] = absorb_allowed_for_chain(
                    chain_features[canon], iso_params
                )
            canons = canonical_order(
                exact_counts,
                chain_features,
            )

            assigned = set()

            for canon in canons:

                feat = chain_features[canon]
                allow_suffix_absorb = feat["absorb_allowed"]
                idxs = []

                for ch, idxlist in chain_to_idxs.items():

                    # strict suffix collapse (true 3′ anchoring)
                    compatible = len(ch) <= len(canon) and canon[-len(ch):] == ch
                    if not compatible:
                        continue
                    if ch != canon and not allow_suffix_absorb:
                        continue
                    if compatible:
                        for i in idxlist:
                            if i not in assigned:
                                idxs.append(i)

                if not idxs:
                    continue

                for i in idxs:
                    assigned.add(i)

                grp = [members[i] for i in idxs]

                full_len_members = [
                    m for m in grp if tuple(m["chain_tx"]) == canon
                ]

                rep = max(
                    full_len_members or grp,
                    key=lambda m: (m["exons"][-1][1] - m["exons"][0][0])
                )

                rep_exons = list(rep["exons"])
                tes = rep_pos

                # enforce TES boundary -- but never past the terminal exon's OTHER end, or we
                # would emit an inverted (start > end) exon. rep_pos (the cluster mode) can be up
                # to apa_window from the rep read's own TES, so a terminal exon shorter than the
                # window could otherwise be flipped. Clamp so the exon stays non-empty.
                if strand == "+":
                    if rep_exons[-1][1] != tes and tes > rep_exons[-1][0]:
                        rep_exons[-1] = (rep_exons[-1][0], tes)
                else:
                    if rep_exons[0][0] != tes and tes < rep_exons[0][1]:
                        rep_exons[0] = (tes, rep_exons[0][1])

                polya_ok = sum(
                    1 for m in grp
                    if m["sclen"] >= iso_params["min_polya_length"]
                    and m["purity"] >= iso_params["min_polya_purity"]
                )

                polya_frac = polya_ok / len(grp) if grp else 0.0

                isoforms.append(dict(
                    chrom=chrom,
                    strand=strand,
                    tes=tes,
                    chain_tx=canon,
                    n_introns=len(canon),
                    members=grp,
                    rep_exons=rep_exons,
                    polya_frac=polya_frac,
                    exact_chain_reads=exact_counts[canon],
                    trunc_assigned_reads=len(grp) - exact_counts[canon],
                    family_reachable_reads=feat["reachable_count"],
                    anchor_reads=feat["anchor_reads"],
                    anchor_frac=feat["anchor_frac"],
                    absorb_allowed=int(allow_suffix_absorb),
                    distal_unique_prefix=feat["unique_prefix"],
                    distal_unique_5p_junction=feat["distal_unique_5p_junction"],
                    exact_sample_ct=Counter(m["sample"] for m in full_len_members),
                    assignment_mode="support_first",
                ))

    return isoforms

def assign_metagene_partitions(final_kept):
    gene_records = {}
    for iso in final_kept:
        gidx = iso["gene_index"]
        rec = gene_records.setdefault(gidx, dict(
            gene_index=gidx,
            chrom=iso["chrom"],
            strand=iso["strand"],
            gene_name=iso["gene_name_label"],
            gene_id=iso["gene_id_label"],
            exons=[],
            isos=[],
        ))
        rec["exons"].extend(iso["rep_exons"])
        rec["isos"].append(iso)

    for rec in gene_records.values():
        rec["exon_union"] = merge_intervals(rec["exons"])
        rec["span_start"] = rec["exon_union"][0][0]
        rec["span_end"] = rec["exon_union"][-1][1]

    parent = {gidx: gidx for gidx in gene_records}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    by_cs = defaultdict(list)
    for rec in gene_records.values():
        by_cs[(rec["chrom"], rec["strand"])].append(rec)

    for key in by_cs:
        lst = sorted(by_cs[key], key=lambda g: (g["span_start"], g["span_end"], g["gene_index"]))
        for i, g1 in enumerate(lst):
            j = i + 1
            while j < len(lst) and lst[j]["span_start"] <= g1["span_end"]:
                g2 = lst[j]
                if first_overlap_bp(g1["exon_union"], g2["exon_union"]) > 0:
                    union(g1["gene_index"], g2["gene_index"])
                j += 1

    components = defaultdict(list)
    for rec in gene_records.values():
        components[find(rec["gene_index"])].append(rec)

    ordered_components = sorted(
        components.values(),
        key=lambda comp: (
            comp[0]["chrom"],
            comp[0]["strand"],
            min(g["span_start"] for g in comp),
            min(g["gene_index"] for g in comp),
        ),
    )

    for midx, comp in enumerate(ordered_components, 1):
        comp_isos = []
        component_gene_indexes = sorted(g["gene_index"] for g in comp)
        for grec in comp:
            comp_isos.extend(grec["isos"])
        comp_isos = sorted(
            comp_isos,
            key=lambda iso: (-iso["count"], iso["gene_index"], iso["gene_tx_index"], iso["tes"]),
        )

        colored = []
        max_partition = 0
        for iso in comp_isos:
            used = {
                other["zn_index"]
                for other in colored
                if first_overlap_bp(iso["rep_exons"], other["rep_exons"]) > 0
            }
            zn_index = smallest_unused_positive(used)
            max_partition = max(max_partition, zn_index)
            iso["metagene_index"] = midx
            iso["zn_index"] = zn_index
            iso["metagene_gene_indexes"] = "|".join(str(x) for x in component_gene_indexes)
            colored.append(iso)

        for iso in comp_isos:
            iso["metagene_partition_count"] = max_partition

    return final_kept
  
# --------------------------- main ---------------------------

def main():
    ap = argparse.ArgumentParser(description="TES/APA assembler v9 (regionized)")
    # Inputs
    ap.add_argument("--bams", nargs="+")
    ap.add_argument("--glob")
    ap.add_argument("--dir")
    ap.add_argument("--region", help="Limit to chrom[:start-end]")
    ap.add_argument("--gtf")
    ap.add_argument("--threads", type=int, default=0)

    # Read filters
    ap.add_argument("--primary-only", action="store_true")
    ap.add_argument("--min-mapq", type=int, default=1)
    ap.add_argument("--min-introns-read", type=int, default=0)
    ap.add_argument("--require-softclip3p", type=int, default=0)

    # Isoform support thresholds
    ap.add_argument("--apa-window", type=int, default=20)
    ap.add_argument("--tes-window", type=int, default=None)
    ap.add_argument("--min-reads", type=int, default=10)
    ap.add_argument("--min-frac", type=float, default=0.05)
    ap.add_argument("--min-introns", type=int, default=0)

    # Poly(A) heuristics
    ap.add_argument("--min-polya-length", type=int, default=12)
    ap.add_argument("--min-polya-purity", type=float, default=0.7)
    ap.add_argument("--polya-support-frac", type=float, default=0.6)

    # Annotation behavior
    ap.add_argument("--tes-match-tol", type=int, default=25)
    ap.add_argument("--exact-tes-tol", type=int, default=10)
    ap.add_argument("--min-distal-anchor-reads", type=int, default=2)
    ap.add_argument("--min-distal-anchor-frac", type=float, default=0.05)
    ap.add_argument("--min-exact-canonical-reads", type=int, default=1)

    # Outputs
    ap.add_argument("--write-zt-bams", action="store_true")
    ap.add_argument("--write-zt-tagged-sample-bams", action="store_true")
    ap.add_argument("--emit-modkit-manifest", action="store_true")
    ap.add_argument("--min-reads-per-sample-for-mod", type=int, default=5)
    ap.add_argument("--min-total-reads-for-mod", type=int, default=20)
    ap.add_argument("--out-gtf", default="readbacked_annot.gtf")
    ap.add_argument("--out-prefix", default=None)

    # Regionization (new)
    ap.add_argument("--cov-bin-bp", type=int, default=500, help="Coarse bin size for gap finding")
    ap.add_argument("--min-gap-bins", type=int, default=2, help="Min consecutive all-zero bins to call a gap")
    ap.add_argument("--pad-fetch-bp", type=int, default=2000, help="Padding when fetching per-core reads")
    ap.add_argument("--max-breakpoints-per-chrom", type=int, default=50,
                    help="Deterministic random breaks inside all-zero gaps")
    ap.add_argument("--rand-seed", type=int, default=1337)
    ap.add_argument("--status-every", type=int, default=0,
                    help="Print regionization read-count status every N reads per BAM (0=off)")

    args = ap.parse_args()

    # Collect BAMs
    bams = []
    if args.bams: bams += args.bams
    if args.dir and args.glob:
        bams += glob.glob(os.path.join(args.dir, args.glob))
    elif args.dir:
        bams += glob.glob(os.path.join(args.dir, "*.bam"))
    bams = sorted(set(bams))
    if not bams:
        sys.exit("No BAMs found")

    # Load GTF (+ index)
    gtf_txs = load_gtf(args.gtf) if args.gtf else []
    if args.gtf:
        print(f"[INFO] Loaded {len(gtf_txs)} transcripts from {args.gtf}", file=sys.stderr)
        gtf_by_cs, gtf_tes_sorted = build_gtf_index(gtf_txs)
    else:
        print("[INFO] No GTF supplied", file=sys.stderr)
        gtf_by_cs, gtf_tes_sorted = {}, {}

    # Parse optional region limit
    only_chrom = None; only_start0=None; only_end0=None
    if args.region:
        s = args.region
        if ":" in s:
            chrom, rest = s.split(":",1)
            only_chrom = chrom
            if "-" in rest:
                a,b = rest.split("-",1)
                try:
                    only_start0 = int(a.replace(",",""))
                    only_end0 = int(b.replace(",",""))
                except Exception:
                    only_start0 = None; only_end0 = None
        else:
            only_chrom = s

    # Chrom sizes from first BAM
    chrom_sizes = []
    with pysam.AlignmentFile(bams[0], "rb") as fh0:
        sq = fh0.header.get("SQ", [])
        for ent in sq:
            chrom = ent.get("SN")
            ln = int(ent.get("LN", 0))
            if chrom and ln>0:
                chrom_sizes.append((chrom, ln))

    # Regionize -- one task per chromosome. This pre-pass streams every read of every BAM, so
    # running it single-threaded (as it used to) left C-1 cores idle for an O(total_reads) scan.
    # Each chromosome is independent and deterministic, so fan it out across processes; with
    # threads<=1 (or a single chromosome) it runs serially exactly as before.
    print("[INFO] Regionizing via all-sample zero-coverage gaps ...", file=sys.stderr)
    n_threads = max(1, int(args.threads or 0))
    regionize_tasks = []
    for chrom, clen in chrom_sizes:
        if only_chrom and chrom != only_chrom:
            continue
        c0 = only_start0 if (only_chrom == chrom and only_start0 is not None) else 0
        c1 = only_end0 if (only_chrom == chrom and only_end0 is not None) else clen
        regionize_tasks.append((chrom, clen, bams, args.cov_bin_bp, args.status_every,
                                args.min_gap_bins, args.max_breakpoints_per_chrom,
                                args.rand_seed, c0, c1))

    region_results = {}
    if n_threads > 1 and len(regionize_tasks) > 1:
        # Each worker holds one open handle (+index) per BAM, so bound concurrency to keep the
        # total open-file count sane on many-sample runs (e.g. 64 threads x 6 BAMs = 384 fds).
        max_workers = min(n_threads, len(regionize_tasks), max(1, 256 // max(1, len(bams))))
        try:
            with ProcessPoolExecutor(max_workers=max_workers) as ex:
                for chrom, n_cores, kept in ex.map(_regionize_chrom_worker, regionize_tasks):
                    region_results[chrom] = (n_cores, kept)
        except Exception as exc:
            print(f"[WARN] Falling back to serial regionization: {exc}", file=sys.stderr)
            region_results = {}
    if not region_results:
        # Serial path: open the BAM handles ONCE and reuse them across chromosomes (avoids
        # re-opening per contig on references with hundreds of contigs).
        handles = [pysam.AlignmentFile(bp, "rb") for bp in bams]
        try:
            for task in regionize_tasks:
                (chrom, clen, _bam_paths, cov_bin_bp, status_every, min_gap_bins,
                 max_breaks, rand_seed, c0, c1) = task
                chrom, n_cores, kept = _regionize_chrom_core(
                    handles, chrom, clen, cov_bin_bp, status_every, min_gap_bins,
                    max_breaks, rand_seed, c0, c1)
                region_results[chrom] = (n_cores, kept)
        finally:
            for fh in handles:
                try: fh.close()
                except Exception: pass

    # Reassemble in the ORIGINAL chromosome order so cores_all is order-identical to the
    # serial path (workers may finish in any order).
    cores_all = []
    for chrom, _clen in chrom_sizes:
        if chrom not in region_results:
            continue
        n_cores, kept = region_results[chrom]
        cores_all.extend(kept)
        print(f"[INFO] {chrom}: cores_total={n_cores} kept_in_region={len(kept)}", file=sys.stderr)

    if not cores_all:
        sys.exit("No cores to process (after regionization)")

    # Prepare args for workers
    filters = dict(
        primary_only=bool(args.primary_only),
        min_mapq=int(args.min_mapq),
        require_softclip3p=int(args.require_softclip3p),
        min_introns_read=int(args.min_introns_read),
    )
    iso_params = dict(
        apa_window=int(args.tes_window if args.tes_window is not None else args.apa_window),
        min_polya_length=int(args.min_polya_length),
        min_polya_purity=float(args.min_polya_purity),
        min_distal_anchor_reads=int(args.min_distal_anchor_reads),
        min_distal_anchor_frac=float(args.min_distal_anchor_frac),
        min_exact_canonical_reads=int(args.min_exact_canonical_reads),
    )
    worker_args = [
        (chrom, s, e, int(args.pad_fetch_bp), bams, filters, iso_params)
        for (chrom, s, e) in cores_all
    ]

    # Process cores (parallel or serial)
    kept_isoforms = []
    n_threads = max(1, int(args.threads or 0))
    print(f"[INFO] Processing {len(worker_args)} cores with threads={n_threads}", file=sys.stderr)
    if n_threads > 1:
        try:
            with ProcessPoolExecutor(max_workers=n_threads) as ex:
                futs = [ex.submit(_process_core, wa) for wa in worker_args]
                for i, fut in enumerate(as_completed(futs), 1):
                    out = fut.result()
                    kept_isoforms.extend(out)
                    if i % max(1, len(worker_args)//20) == 0:
                        print(f"[INFO] cores done: {i}/{len(worker_args)}", file=sys.stderr)
        except PermissionError as exc:
            print(f"[WARN] Falling back to serial core processing: {exc}", file=sys.stderr)
            for i, wa in enumerate(worker_args, 1):
                out = _process_core(wa)
                kept_isoforms.extend(out)
                if i % max(1, len(worker_args)//20) == 0:
                    print(f"[INFO] cores done: {i}/{len(worker_args)}", file=sys.stderr)
    else:
        for i, wa in enumerate(worker_args, 1):
            out = _process_core(wa)
            kept_isoforms.extend(out)
            if i % max(1, len(worker_args)//20) == 0:
                print(f"[INFO] cores done: {i}/{len(worker_args)}", file=sys.stderr)

    if not kept_isoforms:
        sys.exit("No candidate isoforms found")

    # De-duplicate fragmentforms emitted by more than one regionization core.
    # Cores tile the contig with shared boundaries, and each core fetches its reads with
    # +/- pad_fetch_bp of padding, so a fragmentform whose TES lands exactly on a core
    # boundary (e.g. a shared 3' end at the edge of a coverage gap) is assembled by BOTH
    # adjacent cores. Those copies are the SAME fragmentform -- identical
    # (chrom, strand, tes, chain_tx) -- so collapse them, keeping the copy with the most
    # members. The core that actually contains the locus always sees the full read set,
    # while a padded neighbour may see only the reads within pad_fetch_bp of the boundary,
    # so "max members" never drops reads. (Without this, a shared-3'-end isoform pair is
    # double-emitted and its reads are double-counted in the metrics / tx_counts / usage
    # tables; the ZN read tags are unaffected because each read is tagged once.)
    _dedup: dict[tuple, dict] = {}
    for iso in kept_isoforms:
        key = (iso["chrom"], iso["strand"], iso["tes"], tuple(iso["chain_tx"]))
        prev = _dedup.get(key)
        if prev is None or len(iso["members"]) > len(prev["members"]):
            _dedup[key] = iso
    if len(_dedup) != len(kept_isoforms):
        print(f"[INFO] collapsed {len(kept_isoforms) - len(_dedup)} duplicate fragmentform(s) "
              f"straddling a regionization core boundary", file=sys.stderr)
    kept_isoforms = list(_dedup.values())

    # Global filtering + metrics (same as single-pass)
    total_reads_used = sum(len(iso["members"]) for iso in kept_isoforms)
    final_kept = []
    metrics_rows = []
    for iso in kept_isoforms:
        chrom = iso["chrom"]; strand = iso["strand"]; tes = iso["tes"]
        chain_tx = iso["chain_tx"]; n_introns = iso["n_introns"]
        members = iso["members"]; rep_exons = iso["rep_exons"]; polya_frac = iso["polya_frac"]
        count = len(members)
        frac_global = count/total_reads_used if total_reads_used else 0.0
        med_mapq = median([m["mapq"] for m in members])
        med_sclen = median([m["sclen"] for m in members])
        med_purity = median([m["purity"] for m in members])
        sample_ct = Counter(m["sample"] for m in members)
        chain_counts = Counter(chain_to_str(tuple(m["chain_tx"])) for m in members)
        n_full_len_reads = chain_counts.get(chain_to_str(chain_tx), 0)
        included_trunc = 1 if any(len(tuple(m["chain_tx"])) < len(chain_tx) for m in members) else 0

        keep = True
        if count < args.min_reads or frac_global < args.min_frac: keep=False
        if n_introns < args.min_introns: keep=False
        if polya_frac < args.polya_support_frac: keep=False

        metrics_rows.append([
            chrom, tes, strand,
            chain_to_str(chain_tx), n_introns,
            count, f"{frac_global:.4f}", f"{polya_frac:.4f}",
            f"{med_mapq:.1f}", f"{med_sclen:.1f}", f"{med_purity:.3f}",
            len(chain_counts), n_full_len_reads, included_trunc,
            iso["exact_chain_reads"], iso["trunc_assigned_reads"], iso["family_reachable_reads"],
            iso["anchor_reads"], f"{iso['anchor_frac']:.4f}", iso["absorb_allowed"],
            chain_to_str(iso["distal_unique_prefix"]),
            (f"{iso['distal_unique_5p_junction'][0]}-{iso['distal_unique_5p_junction'][1]}"
             if iso["distal_unique_5p_junction"] else "."),
            "|".join(f"{k}:{v}" for k,v in sorted(iso["exact_sample_ct"].items())),
            iso["assignment_mode"],
            "|".join(f"{k}:{v}" for k,v in sorted(sample_ct.items())),
            "|".join(f"{k}:{v}" for k,v in sorted(chain_counts.items())),
            int(keep)
        ])

        if keep:
            final_kept.append(dict(
                chrom=chrom, strand=strand, tes=tes, chain_tx=chain_tx,
                n_introns=n_introns, members=members, rep_exons=rep_exons,
                count=count, frac_global=frac_global, polya_frac=polya_frac,
                med_mapq=med_mapq, med_sclen=med_sclen, med_purity=med_purity,
                sample_ct=sample_ct, chain_counts=chain_counts,
                exact_chain_reads=iso["exact_chain_reads"],
                trunc_assigned_reads=iso["trunc_assigned_reads"],
                family_reachable_reads=iso["family_reachable_reads"],
                anchor_reads=iso["anchor_reads"],
                anchor_frac=iso["anchor_frac"],
                absorb_allowed=iso["absorb_allowed"],
                distal_unique_prefix=iso["distal_unique_prefix"],
                distal_unique_5p_junction=iso["distal_unique_5p_junction"],
                exact_sample_ct=iso["exact_sample_ct"],
                assignment_mode=iso["assignment_mode"],
            ))

    prefix = args.out_prefix if args.out_prefix else args.out_gtf.replace(".gtf","")
    os.makedirs(os.path.dirname(args.out_gtf) or ".", exist_ok=True)

    metrics_path = f"{prefix}_metrics.tsv"
    with open(metrics_path, "w") as m:
        m.write("#chrom\ttes_1based\tstrand\tintron_chain_tx_order\tn_introns\tread_support\tfrac_global\tpolya_support_frac\tmedian_mapq\tmedian_tail_len\tmedian_tail_purity\tn_unique_chains\tn_full_length_reads\tincluded_trunc\texact_chain_reads\ttrunc_assigned_reads\tfamily_reachable_reads\tanchor_reads\tanchor_frac\tabsorb_allowed\tdistal_unique_prefix_tx\tdistal_unique_5p_junction\texact_sample_counts\tassignment_mode\tsample_counts\tchain_counts\tkept\n")
        for row in metrics_rows:
            m.write("\t".join(map(str,row))+"\n")

    if not final_kept:
        sys.exit(f"No isoforms passed filters. See metrics: {metrics_path}")

    # Annotate & bucket
    def annotate_isoform_safe(iso):
        return annotate_isoform(iso, gtf_by_cs, gtf_tes_sorted, args.tes_match_tol, args.exact_tes_tol) if gtf_txs else dict(
            classification="NOVEL_LOCUS", gene_id="NA", gene_name="NA", matched_tid="NA",
            gtf_tes="NA", gtf_chain_tx=tuple(), tes_delta_bp="NA", exon_overlap_bp=0, match_source="NONE")

    for iso in final_kept:
        iso["annotation"] = annotate_isoform_safe(iso)

    # --- Novel loci: group by genomic OVERLAP, and name uniquely + deterministically. ---
    # The old key was (f"NOVEL_{chrom}_{strand}", f"{chrom}:{strand}:{tes}"), which was wrong in
    # both directions: the display half is identical for every novel locus on a (chrom, strand)
    # -- so distinct loci collide downstream (the report groups by gene_name; aggregate_by_gene
    # names per-gene files by gene_name and silently overwrites one locus with another) -- while
    # the id half keys on the exact TES, so ONE novel locus with APA is split into several
    # "genes". Instead, merge novel fragmentforms whose exon spans overlap on the same
    # (chrom, strand) into one locus, and name it after its span.
    novel_isos = [iso for iso in final_kept
                  if iso["annotation"]["gene_name"] == "NA" and iso["annotation"]["gene_id"] == "NA"]
    novel_by_cs = defaultdict(list)
    for iso in novel_isos:
        novel_by_cs[(iso["chrom"], iso["strand"])].append(iso)
    for (chrom, strand), isos in novel_by_cs.items():
        spans = sorted((i["rep_exons"][0][0], i["rep_exons"][-1][1], k) for k, i in enumerate(isos))
        comps = []  # [start, end, [iso indices]] -- connected components by span overlap
        for s, e, k in spans:
            if comps and s <= comps[-1][1]:
                comps[-1][1] = max(comps[-1][1], e)
                comps[-1][2].append(k)
            else:
                comps.append([s, e, [k]])
        tok = {"+": "plus", "-": "minus"}.get(strand, "na")
        for ordinal, (s, e, ks) in enumerate(comps, 1):
            name = f"NOVEL_{chrom}_{tok}_{ordinal:03d}_{s}_{e}"
            lid = f"{chrom}:{strand}:{s}-{e}"
            for k in ks:
                isos[k]["_novel_locus"] = (name, lid)

    buckets = defaultdict(list)
    for iso in final_kept:
        ann = iso["annotation"]
        if ann["gene_name"] != "NA" and ann["gene_id"] != "NA":
            gkey = (ann["gene_name"], ann["gene_id"])
        else:
            gkey = iso["_novel_locus"]
        buckets[gkey].append(iso)

    gene_keys_sorted = sorted(buckets.keys(), key=lambda x: (x[0], x[1]))
    gene_index = {gk: i+1 for i, gk in enumerate(gene_keys_sorted)}

    for gk in gene_keys_sorted:
        iso_list = buckets[gk]
        iso_list_sorted = sorted(iso_list, key=lambda x: (-x["count"], x["tes"]))
        for tidx, iso in enumerate(iso_list_sorted, 1):
            gn, gid = gk
            gidx = gene_index[gk]
            iso["gene_name_label"], iso["gene_id_label"] = gn, gid
            iso["gene_index"] = gidx
            iso["gene_tx_index"] = tidx
            iso["zt_label"] = f"{gn}.{gid}.G{gidx}.T{tidx}"

    assign_metagene_partitions(final_kept)

    # Write GTF
    with open(args.out_gtf, "w") as out:
        for iso in sorted(final_kept, key=lambda x: (x["gene_index"], x["gene_tx_index"])):
            chrom=iso["chrom"]; strand=iso["strand"]; tes=iso["tes"]
            rep_exons=iso["rep_exons"]; chain_tx=iso["chain_tx"]
            t_start, t_end = rep_exons[0][0], rep_exons[-1][1]
            ann = iso["annotation"]
            tid = f"{iso['gene_id_label']}.G{iso['gene_index']}.T{iso['gene_tx_index']}"
            attrs = (
                f'gene_id "{iso["gene_id_label"]}"; '
                f'transcript_id "{tid}"; '
                f'ref_gene_name "{iso["gene_name_label"]}"; '
                f'zt_label "{iso["zt_label"]}"; gene_index "{iso["gene_index"]}"; transcript_index "{iso["gene_tx_index"]}"; '
                f'metagene_index "{iso["metagene_index"]}"; zn_index "{iso["zn_index"]}"; metagene_partition_count "{iso["metagene_partition_count"]}"; '
                f'read_support "{iso["count"]}"; frac_support_global "{iso["frac_global"]:.4f}"; polya_support_frac "{iso["polya_frac"]:.4f}"; '
                f'intron_chain "{chain_to_str(chain_tx)}"; tes "{tes}"; classification "{ann["classification"]}"; '
                f'matched_tid "{ann["matched_tid"]}"; matched_gid "{ann["gene_id"]}"; match_source "{ann["match_source"]}";'
            )
            out.write(f"{chrom}\tReadBacked\ttranscript\t{t_start}\t{t_end}\t1000\t{strand}\t.\t{attrs}\n")
            for j,(s,e) in enumerate(rep_exons,1):
                out.write(f"{chrom}\tReadBacked\texon\t{s}\t{e}\t1000\t{strand}\t.\t{attrs} exon_number \"{j}\";\n")

    # Classification summary
    summary_path = f"{prefix}_classification_summary.tsv"
    with open(summary_path, "w") as s:
        s.write("#code\tzt_label\tgene_index\ttranscript_index\tmetagene_index\tzn_index\tmetagene_partition_count\tchrom\tstrand\tiso_tes\tiso_chain_tx\tgtf_gene_id\tgtf_gene_name\tgtf_transcript_id\tgtf_tes\tgtf_chain_tx\ttes_delta_bp\texon_overlap_bp\tmatch_source\tclassification\tread_support\texact_chain_reads\ttrunc_assigned_reads\tfamily_reachable_reads\tanchor_reads\tanchor_frac\tabsorb_allowed\tdistal_unique_prefix_tx\texact_sample_counts\tassignment_mode\tfrac_global\tpolya_support_frac\tsample_counts\n")
        for iso in sorted(final_kept, key=lambda x: (x["gene_index"], x["gene_tx_index"])):
            ann = iso["annotation"]
            s.write("\t".join(map(str,[
                iso["zt_label"], iso["zt_label"], iso["gene_index"], iso["gene_tx_index"], iso["metagene_index"], iso["zn_index"], iso["metagene_partition_count"], iso["chrom"], iso["strand"], iso["tes"],
                chain_to_str(iso["chain_tx"]),
                ann["gene_id"], ann["gene_name"], ann["matched_tid"],
                ann["gtf_tes"], chain_to_str(ann["gtf_chain_tx"]), ann["tes_delta_bp"], ann["exon_overlap_bp"], ann["match_source"],
                ann["classification"], iso["count"], iso["exact_chain_reads"], iso["trunc_assigned_reads"], iso["family_reachable_reads"],
                iso["anchor_reads"], f"{iso['anchor_frac']:.4f}", iso["absorb_allowed"], chain_to_str(iso["distal_unique_prefix"]),
                "|".join(f"{k}:{v}" for k,v in sorted(iso["exact_sample_ct"].items())), iso["assignment_mode"],
                f"{iso['frac_global']:.3f}", f"{iso['polya_frac']:.4f}",
                "|".join(f"{k}:{v}" for k,v in sorted(iso["sample_ct"].items()))
            ]))+"\n")

    # === Extra outputs: counts, PCA, per-sample stats ===
    import pandas as pd
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from plot_utils import save_figure

    all_samples = sorted({s for iso in final_kept for s in iso["sample_ct"].keys()})
    tx_ids = [iso["zt_label"] for iso in final_kept]
    data = []
    for iso in final_kept:
        data.append([iso["sample_ct"].get(s, 0) for s in all_samples])
    counts_df = pd.DataFrame(data, index=tx_ids, columns=all_samples)
    counts_path = f"{prefix}_tx_counts.tsv"
    counts_df.to_csv(counts_path, sep="\t")

    X = np.log1p(counts_df.values.T)
    pca_png = f"{prefix}_tx_counts.pca.png"
    if X.size and X.shape[0] >= 1 and X.shape[1] >= 1:
        Xc = X - X.mean(axis=0, keepdims=True)
        U, S, Vt = np.linalg.svd(Xc, full_matrices=False)
        pcs = min(2, S.size)
        PC = U[:, :pcs] * S[:pcs] if pcs else np.zeros((X.shape[0], 0))
        var_expl = (S**2) / (S**2).sum() if S.size else np.zeros_like(S)

        fig, ax = plt.subplots(figsize=(8.6, 6.8), layout="constrained")
        cmap = plt.get_cmap("tab10")
        colors = [cmap(i % cmap.N) for i in range(len(all_samples))]
        x = PC[:, 0] if pcs >= 1 else np.zeros(X.shape[0])
        y = PC[:, 1] if pcs >= 2 else np.zeros(X.shape[0])

        for idx, name in enumerate(all_samples):
            ax.scatter(
                x[idx],
                y[idx],
                s=90,
                color=colors[idx],
                edgecolors="white",
                linewidths=1.2,
                label=name,
                zorder=3,
            )

        ax.axhline(0.0, color="#c9d1d8", linewidth=1.0, zorder=1)
        ax.axvline(0.0, color="#c9d1d8", linewidth=1.0, zorder=1)
        ax.grid(True, linestyle="--", linewidth=0.6, alpha=0.35)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_color("#4d5761")
        ax.spines["bottom"].set_color("#4d5761")
        ax.tick_params(colors="#2f353a")

        ax.set_xlabel(f"PC1 ({var_expl[0]*100:.1f}% variance)" if pcs >= 1 else "PC1")
        ax.set_ylabel(f"PC2 ({var_expl[1]*100:.1f}% variance)" if pcs >= 2 else "PC2")
        ax.set_title("Sample PCA of log1p fragmentform counts")

        if all_samples:
            # 'outside lower center' -> constrained_layout reserves a band below the axes for the
            # legend, so it never collides with the x-axis label at the enlarged house font sizes.
            fig.legend(
                loc="outside lower center",
                ncol=min(2, max(1, len(all_samples))),
                frameon=False,
                fontsize=8,
                handletextpad=0.5,
                columnspacing=1.2,
            )

        save_figure(fig, pca_png, dpi=300, bbox_inches="tight")   # PNG + PDF + SVG
        plt.close(fig)
    else:
        fig, ax = plt.subplots(figsize=(8.6, 6.8))
        ax.set_title("Sample PCA (no data)")
        ax.axis("off")
        save_figure(fig, pca_png, dpi=300, bbox_inches="tight")   # PNG + SVG
        plt.close(fig)

    rows = []
    for sname in all_samples:
        col = counts_df[sname]
        detected = (col > 0)
        total_reads = int(col.sum())
        n_tx = int(detected.sum())
        med_reads = float(col[detected].median()) if n_tx > 0 else 0.0
        rows.append(dict(sample=sname, total_reads=total_reads, n_transcripts=n_tx, median_reads_per_tx=med_reads))
    stats_df = pd.DataFrame(rows).sort_values("sample")
    stats_path = f"{prefix}_per_sample_stats.tsv"
    stats_df.to_csv(stats_path, sep="\t", index=False)

    tx_read_length_rows = []
    partition_rows = []
    for iso in sorted(final_kept, key=lambda x: (x["gene_index"], x["gene_tx_index"])):
        read_lengths = [m["read_length"] for m in iso["members"] if m.get("read_length", 0) > 0]
        mean_len = (sum(read_lengths) / len(read_lengths)) if read_lengths else 0.0
        tx_read_length_rows.append(dict(
            code=iso["zt_label"],
            zt_label=iso["zt_label"],
            gene_index=iso["gene_index"],
            transcript_index=iso["gene_tx_index"],
            metagene_index=iso["metagene_index"],
            zn_index=iso["zn_index"],
            metagene_partition_count=iso["metagene_partition_count"],
            gene_id=iso["gene_id_label"],
            gene_name=iso["gene_name_label"],
            chrom=iso["chrom"],
            strand=iso["strand"],
            assigned_reads=len(read_lengths),
            mean_read_length=round(mean_len, 2),
            median_read_length=round(float(median(read_lengths)), 2) if read_lengths else 0.0,
            min_read_length=min(read_lengths) if read_lengths else 0,
            max_read_length=max(read_lengths) if read_lengths else 0,
        ))
        partition_rows.append(dict(
            code=iso["zt_label"],
            zt_label=iso["zt_label"],
            gene_index=iso["gene_index"],
            transcript_index=iso["gene_tx_index"],
            metagene_index=iso["metagene_index"],
            zn_index=iso["zn_index"],
            metagene_partition_count=iso["metagene_partition_count"],
            gene_id=iso["gene_id_label"],
            gene_name=iso["gene_name_label"],
            chrom=iso["chrom"],
            strand=iso["strand"],
            tes=iso["tes"],
            read_support=iso["count"],
            exact_chain_reads=iso["exact_chain_reads"],
            trunc_assigned_reads=iso["trunc_assigned_reads"],
            anchor_reads=iso["anchor_reads"],
            anchor_frac=round(float(iso["anchor_frac"]), 4),
        ))
    tx_read_length_df = pd.DataFrame(tx_read_length_rows)
    tx_read_lengths_path = f"{prefix}_tx_assigned_read_lengths.tsv"
    tx_read_length_df.to_csv(tx_read_lengths_path, sep="\t", index=False)
    partition_map_path = f"{prefix}_partition_map.tsv"
    pd.DataFrame(partition_rows).to_csv(partition_map_path, sep="\t", index=False)

    # Optional ZT outputs
    need_assign = bool(args.write_zt_bams or args.emit_modkit_manifest or args.write_zt_tagged_sample_bams)
    assign = None
    if need_assign:
        assign = defaultdict(dict)
        for iso in final_kept:
            for m in iso["members"]:
                assign[m["sample"]][m["qname"]] = (
                    iso["zt_label"],
                    iso["gene_index"],
                    iso["zn_index"],
                    iso["metagene_index"],
                )

    if args.write_zt_bams or args.emit_modkit_manifest:
        out_dir = os.path.join(os.path.dirname(args.out_gtf) or ".", "zt_bams")
        os.makedirs(out_dir, exist_ok=True)
        manifest_rows = []
        for iso in final_kept:
            if sum(iso["sample_ct"].values()) < args.min_total_reads_for_mod:
                continue
            for sample, n in iso["sample_ct"].items():
                if n < args.min_reads_per_sample_for_mod:
                    continue
                iso.setdefault("eligible_samples", set()).add(sample)
        for bam in bams:
            sample = os.path.basename(bam).replace(".bam","")
            elig_codes = {iso["zt_label"] for iso in final_kept if "eligible_samples" in iso and sample in iso["eligible_samples"]}
            if not elig_codes:
                continue
            with pysam.AlignmentFile(bam, "rb") as inp:
                writers = {}
                try:
                    for code in sorted(elig_codes):
                        out_bam = os.path.join(out_dir, f"{sample}.{code}.bam")
                        writers[code] = pysam.AlignmentFile(out_bam, "wb", header=inp.header)
                    for aln in (inp.fetch()):
                        if aln.is_unmapped: continue
                        if args.primary_only and (aln.is_secondary or aln.is_supplementary): continue
                        qn = aln.query_name
                        tup = assign.get(sample, {}).get(qn) if assign else None
                        if not tup: continue
                        zt, zg, zn, zm = tup
                        if zt not in writers: continue
                        try:
                            aln.set_tag("ZT", zt, value_type="Z", replace=True)
                        except TypeError:
                            aln.set_tag("ZT", zt)
                        for tag, val in (("ZG", int(zg)), ("ZN", int(zn)), ("ZM", int(zm))):
                            try:
                                aln.set_tag(tag, val, value_type="i", replace=True)
                            except TypeError:
                                aln.set_tag(tag, val)
                        writers[zt].write(aln)
                finally:
                    for w in writers.values():
                        w.close()
                for code in sorted(elig_codes):
                    out_bam = os.path.join(out_dir, f"{sample}.{code}.bam")
                    if os.path.exists(out_bam) and os.path.getsize(out_bam) > 0:
                        try: pysam.index(out_bam)
                        except Exception: pass
                        manifest_rows.append([sample, code, out_bam])
        if args.emit_modkit_manifest:
            mani = os.path.join(out_dir, "modkit_manifest.tsv")
            with open(mani, "w") as f:
                f.write("sample\tcode\tbam\n")
                for row in manifest_rows:
                    f.write("\t".join(row)+"\n")
            print(f"[OK] Wrote modkit manifest: {mani}", file=sys.stderr)

    if args.write_zt_tagged_sample_bams:
        out_dir2 = os.path.join(os.path.dirname(args.out_gtf) or ".", "zt_tagged")
        os.makedirs(out_dir2, exist_ok=True)
        # Parallelize the tagged-BAM writing across samples (was a serial for-loop; on genome-wide
        # data each sample re-reads a ~100 GB BAM and writes a ~40 GB tagged BAM single-threaded --
        # the dominant assemble tail). Samples are independent, so fan them out; each writer also
        # gets BGZF compression threads. Output is byte-identical to the serial path.
        n_threads = max(1, int(args.threads or 0))
        tasks = []
        for bam in bams:
            sample = os.path.basename(bam).replace(".bam", "")
            out_path = os.path.join(out_dir2, f"{sample}.zt_tagged.bam")
            tasks.append((bam, out_path, (assign.get(sample, {}) if assign else {}),
                          bool(args.primary_only), max(2, n_threads // max(1, 2 * len(bams)))))
        n_workers = min(len(tasks), n_threads)
        wrote = []
        if n_workers > 1 and len(tasks) > 1:
            try:
                with ProcessPoolExecutor(max_workers=n_workers) as ex:
                    wrote = list(ex.map(_write_zt_tagged_sample, tasks))
            except Exception as exc:
                print(f"[WARN] parallel zt_tagged write failed ({exc}); serial fallback", file=sys.stderr)
                wrote = []
        if not wrote:
            wrote = [_write_zt_tagged_sample(t) for t in tasks]
        for out_path in wrote:
            print(f"[OK] Wrote ZT/ZN-tagged sample BAM: {out_path}", file=sys.stderr)

    print(f"[OK] Wrote GTF: {args.out_gtf}", file=sys.stderr)
    print(f"[OK] Metrics: {metrics_path}", file=sys.stderr)
    print(f"[OK] Classification summary: {summary_path}", file=sys.stderr)
    print(f"[OK] Counts: {counts_path}", file=sys.stderr)
    print(f"[OK] PCA plot: {pca_png}", file=sys.stderr)
    print(f"[OK] Per-sample stats: {stats_path}", file=sys.stderr)
    print(f"[OK] Transcript assigned read lengths: {tx_read_lengths_path}", file=sys.stderr)
    print(f"[OK] Partition map: {partition_map_path}", file=sys.stderr)

if __name__ == "__main__":
    main()
