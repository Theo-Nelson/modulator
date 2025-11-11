#!/usr/bin/env python3
"""
Read-backed TES/APA assembler (v9)

What's new vs v8
----------------
1) New outputs (besides GTF, metrics, and classification summary):
   - <prefix>_tx_counts.tsv               : transcript (rows) x sample (cols) counts for KEPT transcripts
   - <prefix>_tx_counts.pca.png           : per-sample PCA scatter from log1p(counts)
   - <prefix>_per_sample_stats.tsv        : basic per-sample stats from the counts table
2) New CLI:
   - --out-prefix: base path for all auxiliary outputs (defaults to out-gtf without .gtf)
3) Keeps v8 features:
   - Robust writer guard for per-transcript BAMs
   - Deterministic 5' truncation assignment prioritizing strongest isoform
   - ZT/ZG/ZN tagging, per-sample tagged BAMs

Tags on assigned reads
----------------------
- ZN (i): transcript_index within its gene (1..k)
- ZG (i): gene_index (1..G) – deterministic within run
- ZT (Z): "{gene_name}.{gene_id}.G{ZG}.T{ZN}" (human-readable)

Standard Outputs
----------------
- <out>.gtf and <out>_classification_summary.tsv
- <prefix>_metrics.tsv
- zt_tagged/<sample>.zt_tagged.bam  (ALL reads retained; tags added where available)
- zt_bams/<sample>.<ZT>.bam         (optional per-transcript subset BAMs)
"""

import argparse, os, sys, glob, re, gzip
from collections import defaultdict, Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
import pysam

# ---------- utils ----------

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
    if aln.has_tag("ts"):
        v = aln.get_tag("ts")
        if v in ("+","-"): return v
    if aln.has_tag("XS"):
        v = aln.get_tag("XS")
        if v in ("+","-"): return v
    return "-" if aln.is_reverse else "+"

def tes_pos1(aln, tx):
    return aln.reference_end if tx=="+" else aln.reference_start+1

def intron_chain_1based(aln):
    chain = []
    ref = aln.reference_start
    for op, ln in (aln.cigartuples or []):
        if op == 3:  # N (splice)
            donor = ref + 1
            acceptor = ref + ln
            chain.append((donor, acceptor))
            ref += ln
        elif op in (0,2,7,8):  # M/D/=/X advance reference
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
    if not seq: return 0.0
    s = seq.upper()
    base = "A" if tx=="+" else "T"
    p1 = s.count(base)/len(s)
    comp = "T" if base=="A" else "A"
    p2 = s.count(comp)/len(s)
    return max(p1, p2)

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

# ---------- Streaming Cluster Manager ----------

class StreamingClusterManager:
    """Manages active TES clusters with automatic flushing as coordinates advance."""
    def __init__(self, apa_window):
        self.apa_window = apa_window
        self.active_clusters = defaultdict(list)  # key: (chrom, strand) -> list of clusters
        self.finalized_clusters = []  # clusters ready for isoform processing
        self.last_pos = defaultdict(lambda: -1)  # key: (chrom, strand) -> last TES seen
        
    def add_read(self, read_dict):
        """Add a read to appropriate cluster or create new cluster."""
        chrom = read_dict["chrom"]
        strand = read_dict["strand"]
        tes = read_dict["tes"]
        key = (chrom, strand)
        
        # Flush old clusters when we've moved far enough ahead
        if self.last_pos[key] >= 0:
            flush_threshold = tes - (self.apa_window * 3)  # conservative margin
            self._flush_clusters_before(key, flush_threshold)
        
        self.last_pos[key] = max(self.last_pos[key], tes)
        
        # Find or create cluster
        cluster_list = self.active_clusters[key]
        matched = None
        for cluster in cluster_list:
            if abs(cluster["rep_tes"] - tes) <= self.apa_window:
                matched = cluster
                break
        
        if matched:
            matched["members"].append(read_dict)
            # Update rep_tes to median (simple approach: just track)
            matched["positions"].append(tes)
        else:
            # New cluster
            new_cluster = {
                "chrom": chrom,
                "strand": strand,
                "rep_tes": tes,
                "positions": [tes],
                "members": [read_dict]
            }
            cluster_list.append(new_cluster)
    
    def _flush_clusters_before(self, key, threshold):
        """Move clusters with max TES before threshold to finalized."""
        cluster_list = self.active_clusters[key]
        remaining = []
        for cluster in cluster_list:
            max_tes = max(cluster["positions"])
            if max_tes < threshold:
                # Finalize this cluster
                self.finalized_clusters.append(cluster)
            else:
                remaining.append(cluster)
        self.active_clusters[key] = remaining
    
    def flush_all(self):
        """Finalize all remaining active clusters."""
        for key, cluster_list in self.active_clusters.items():
            self.finalized_clusters.extend(cluster_list)
        self.active_clusters.clear()
    
    def get_finalized_clusters(self):
        """Return list of finalized clusters grouped by (chrom, strand)."""
        groups = defaultdict(list)
        for cluster in self.finalized_clusters:
            # Recompute representative TES as most common position
            positions = cluster["positions"]
            c = Counter(positions)
            rep = sorted(c.items(), key=lambda x: (x[1], x[0]))[-1][0]
            groups[(cluster["chrom"], cluster["strand"])].append({
                "tes": rep,
                "members": cluster["members"]
            })
        return groups

# ---------- Chunked parallel helpers ----------

def _parse_region_to_parts(region_str):
    """Parse region like 'chr:start-end' or 'chr' into (chrom, start, end or None)."""
    if region_str is None:
        return None
    s = str(region_str)
    if ":" in s:
        chrom, rest = s.split(":", 1)
        if "-" in rest:
            a, b = rest.split("-", 1)
            try:
                return chrom, int(a.replace(",","")), int(b.replace(",",""))
            except Exception:
                return chrom, None, None
        else:
            return chrom, None, None
    return s, None, None

def _make_regions_from_header(bam_path, chunk_bases, overlap_bases, only_region=None):
    """Generate region strings from BAM header with fixed-size chunks and overlap."""
    import pysam
    regs = []
    with pysam.AlignmentFile(bam_path, "rb") as fh:
        sq = fh.header.get("SQ", [])
        if only_region:
            ochr, ost, oend = _parse_region_to_parts(only_region)
        else:
            ochr = ost = oend = None
        for ent in sq:
            chrom = ent.get("SN")
            ln = int(ent.get("LN", 0))
            if not chrom or ln <= 0:
                continue
            if ochr and chrom != ochr:
                continue
            start0 = 1 if (ost is None) else max(1, ost)
            end0 = ln if (oend is None) else min(ln, oend)
            if chunk_bases <= 0:
                regs.append(f"{chrom}:{start0}-{end0}")
            else:
                step = int(chunk_bases)
                ov = int(max(0, overlap_bases))
                s = start0
                while s <= end0:
                    e = min(end0, s + step - 1)
                    # expand with overlap inside available bounds
                    ss = max(start0, s - ov)
                    ee = min(end0, e + ov)
                    regs.append(f"{chrom}:{ss}-{ee}")
                    s = e + 1
    return regs

def _process_region_worker(args):
    """Worker: process a region across all BAMs, return list of clusters dicts.
    Returns a list of {chrom,strand,tes,members} clusters.
    """
    (region, bam_paths, filters, apa_window) = args
    import pysam
    cm = StreamingClusterManager(apa_window)
    for bam in bam_paths:
        if not os.path.exists(bam):
            continue
        sample = os.path.basename(bam).replace(".bam", "")
        with pysam.AlignmentFile(bam, "rb") as fh:
            it = fh.fetch(region=region)
            for aln in it:
                if aln.is_unmapped: continue
                if filters["primary_only"] and (aln.is_secondary or aln.is_supplementary): continue
                if aln.mapping_quality < filters["min_mapq"]: continue
                tx = get_tx_strand(aln)
                sclen, tail = softclip3p_len_and_seq(aln, tx)
                if filters["require_softclip3p"] > 0 and sclen < filters["require_softclip3p"]: continue
                purity = polya_purity(tail, tx) if sclen > 0 else 0.0
                chain = intron_chain_1based(aln)
                if len(chain) < filters["min_introns_read"]: continue
                chain_tx = chain_tx_order(chain, tx)
                chrom = fh.get_reference_name(aln.reference_id)
                read_dict = dict(
                    chrom=chrom, strand=tx, tes=tes_pos1(aln, tx),
                    chain=chain, chain_tx=chain_tx, n_introns=len(chain),
                    exons=exon_blocks_from_aln(aln),
                    bam=os.path.basename(bam), sample=sample, qname=aln.query_name,
                    mapq=aln.mapping_quality, sclen=sclen, purity=purity
                )
                cm.add_read(read_dict)
    cm.flush_all()
    # convert to flat list of clusters
    out = []
    for (chrom, strand), lst in cm.get_finalized_clusters().items():
        for cl in lst:
            out.append({"chrom": chrom, "strand": strand, "tes": cl["tes"], "members": cl["members"]})
    return out

def _merge_region_clusters(cluster_lists, apa_window):
    """Merge clusters from multiple regions by TES proximity within apa_window."""
    by_cs = defaultdict(list)
    for clist in cluster_lists:
        for c in clist:
            by_cs[(c["chrom"], c["strand"])].append(c)
    merged_groups = defaultdict(list)
    for (chrom, strand), lst in by_cs.items():
        if not lst:
            continue
        # positions for clustering
        positions = sorted(c["tes"] for c in lst)
        clusters = cluster_positions(positions, apa_window)
        for cl in clusters:
            rep_pos = cl["rep"]
            members = []
            for c in lst:
                if abs(c["tes"] - rep_pos) <= apa_window:
                    members.extend(c["members"])  # combine
            merged_groups[(chrom, strand)].append(dict(tes=rep_pos, members=members))
    return merged_groups

# ---------- GTF parsing ----------
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

# Build fast lookup structures for GTF transcripts to accelerate annotation
def build_gtf_index(gtf_txs):
    """Return two indices:
    - by_cs: dict[(chrom,strand)] -> list of transcripts for that key
    - tes_sorted: dict[(chrom,strand)] -> list of (tes_1based, transcript) sorted by tes
    """
    by_cs = defaultdict(list)
    for tx in gtf_txs:
        by_cs[(tx.chrom, tx.strand)].append(tx)
    tes_sorted = {}
    for k, lst in by_cs.items():
        tes_sorted[k] = sorted([(t.tes_1based, t) for t in lst], key=lambda x: x[0])
    return by_cs, tes_sorted

# ---------- annotation ----------

def exon_overlap_len(ex1, ex2):
    tot=0
    for s1,e1 in ex1:
        for s2,e2 in ex2:
            lo=max(s1,s2); hi=min(e1,e2)
            if hi>=lo: tot += (hi-lo+1)
    return tot

def annotate_isoform(iso, gtf_by_cs, tes_sorted_index, tes_match_tol, exact_tes_tol):
    chrom, strand, tes = iso["chrom"], iso["strand"], iso["tes"]
    # Candidates by TES proximity using pre-sorted TES list
    key = (chrom, strand)
    cands = []
    if key in tes_sorted_index:
        arr = tes_sorted_index[key]
        # binary search window around tes
        # find left bound
        lo, hi = 0, len(arr)
        while lo < hi:
            mid = (lo + hi) // 2
            if arr[mid][0] < tes - tes_match_tol:
                lo = mid + 1
            else:
                hi = mid
        left = lo
        # iterate until beyond upper bound
        i = left
        upper = tes + tes_match_tol
        while i < len(arr) and arr[i][0] <= upper:
            cands.append(arr[i][1])
            i += 1
    match_source = "TES_TOL" if cands else "NONE"
    best = None
    if cands:
        def rank(tx):
            same_chain = 0 if tx.chain_tx == iso["chain_tx"] else 1
            return (abs(tx.tes_1based - tes), same_chain, -exon_overlap_len(iso["rep_exons"], tx.exons), -len(tx.chain_tx))
        best = min(cands, key=rank)
    else:
        # fallback: overlap search within chrom/strand bin only
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

# ---------- main ----------

def main():
    ap = argparse.ArgumentParser(description="TES/APA assembler v9 with extra counts/PCA/stats outputs")
    # Inputs
    ap.add_argument("--bams", nargs="+")
    ap.add_argument("--glob")
    ap.add_argument("--dir")
    ap.add_argument("--region")
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

    # Outputs
    ap.add_argument("--write-zt-bams", action="store_true",
                    help="Write per-transcript subset BAMs (tags added)")
    ap.add_argument("--write-zt-tagged-sample-bams", action="store_true",
                    help="Write one tagged BAM per input BAM (ALL reads retained)")
    ap.add_argument("--emit-modkit-manifest", action="store_true")
    ap.add_argument("--min-reads-per-sample-for-mod", type=int, default=5)
    ap.add_argument("--min-total-reads-for-mod", type=int, default=20)
    ap.add_argument("--out-gtf", default="readbacked_annot.gtf")
    ap.add_argument("--out-prefix", default=None,
                    help="Prefix for extra outputs (counts/PCA/stats). If omitted, uses out-gtf without .gtf")
    
    # Performance options
    ap.add_argument("--streaming", action="store_true", default=True,
                    help="Use streaming cluster mode (lower memory, default: True)")
    ap.add_argument("--no-streaming", dest="streaming", action="store_false",
                    help="Disable streaming mode (load all reads into memory)")
    ap.add_argument("--progress-interval", type=int, default=100000,
                    help="Print progress every N reads (default: 100000, 0 to disable)")
    ap.add_argument("--chunk-bases", type=int, default=0,
                    help="If >0 and threads>1, process in parallel genomic chunks of this size (bp) with overlap.")
    ap.add_argument("--chunk-overlap-bases", type=int, default=None,
                    help="Optional explicit chunk overlap in bp (default: 3*apa_window).")

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

    # Load GTF + indices
    gtf_txs = load_gtf(args.gtf) if args.gtf else []
    if args.gtf:
        gtf_by_cs, gtf_tes_sorted = build_gtf_index(gtf_txs)
        print(f"[INFO] Loaded {len(gtf_txs)} transcripts from {args.gtf}", file=sys.stderr)
    else:
        gtf_by_cs, gtf_tes_sorted = {}, {}
        print("[INFO] No GTF supplied", file=sys.stderr)

    # Effective window for clustering
    apa_window = args.apa_window if args.apa_window is not None else (args.tes_window if args.tes_window else 20)

    # Decide execution mode: single-process streaming vs chunked parallel
    use_parallel_chunks = bool(args.streaming and args.chunk_bases and args.chunk_bases > 0 and args.threads and args.threads > 1)

    if use_parallel_chunks:
        # Build regions with overlap
        overlap_bases = args.chunk_overlap_bases if args.chunk_overlap_bases is not None else (3 * apa_window)
        # Regions derived from the first BAM header; applied to all BAMs
        regions = _make_regions_from_header(bams[0], args.chunk_bases, overlap_bases, only_region=args.region)
        print(f"[INFO] Parallel chunk mode: {len(regions)} regions, chunk={args.chunk_bases}bp, overlap={overlap_bases}bp, threads={args.threads}", file=sys.stderr)
        # Prepare filters bundle for worker
        filters = dict(primary_only=bool(args.primary_only), min_mapq=int(args.min_mapq),
                       require_softclip3p=int(args.require_softclip3p), min_introns_read=int(args.min_introns_read))
        cluster_lists = []
        with ProcessPoolExecutor(max_workers=args.threads) as ex:
            futs = []
            for reg in regions:
                futs.append(ex.submit(_process_region_worker, (reg, bams, filters, apa_window)))
            for i, fut in enumerate(as_completed(futs), 1):
                cl = fut.result()
                cluster_lists.append(cl)
                if i % max(1, len(regions)//10) == 0:
                    print(f"[INFO] Completed {i}/{len(regions)} regions", file=sys.stderr)
        # Merge clusters from regions
        groups = _merge_region_clusters(cluster_lists, apa_window)
        total_reads_processed = sum(len(c["members"]) for cl in cluster_lists for c in cl)
        if total_reads_processed == 0:
            sys.exit("No usable reads found after filters")
    else:
        # Single-process streaming over whole BAM set
        from pysam import AlignmentFile
        cluster_manager = StreamingClusterManager(apa_window)
        total_reads_processed = 0
        print(f"[INFO] Processing reads in streaming mode (window={apa_window})", file=sys.stderr)
        for bam in bams:
            if not os.path.exists(bam):
                print(f"[WARN] Missing BAM: {bam}", file=sys.stderr); continue
            sample = os.path.basename(bam).replace(".bam","")
            with AlignmentFile(bam, "rb") as fh:
                it = fh.fetch(region=args.region) if args.region else fh.fetch()
                for aln in it:
                    if aln.is_unmapped: continue
                    if args.primary_only and (aln.is_secondary or aln.is_supplementary): continue
                    if aln.mapping_quality < args.min_mapq: continue
                    tx = get_tx_strand(aln)
                    sclen, tail = softclip3p_len_and_seq(aln, tx)
                    if args.require_softclip3p > 0 and sclen < args.require_softclip3p: continue
                    purity = polya_purity(tail, tx) if sclen > 0 else 0.0
                    chain = intron_chain_1based(aln)
                    if len(chain) < args.min_introns_read: continue
                    chain_tx = chain_tx_order(chain, tx)
                    chrom = fh.get_reference_name(aln.reference_id)
                    read_dict = dict(
                        chrom=chrom, strand=tx, tes=tes_pos1(aln, tx),
                        chain=chain, chain_tx=chain_tx, n_introns=len(chain),
                        exons=exon_blocks_from_aln(aln),
                        bam=os.path.basename(bam), sample=sample, qname=aln.query_name,
                        mapq=aln.mapping_quality, sclen=sclen, purity=purity
                    )
                    cluster_manager.add_read(read_dict)
                    total_reads_processed += 1
                    if args.progress_interval > 0 and total_reads_processed % args.progress_interval == 0:
                        n_active = sum(len(clusters) for clusters in cluster_manager.active_clusters.values())
                        n_finalized = len(cluster_manager.finalized_clusters)
                        print(f"[INFO] Processed {total_reads_processed:,} reads (active clusters: {n_active}, finalized: {n_finalized})", file=sys.stderr)
        cluster_manager.flush_all()
        print(f"[INFO] Finalized {len(cluster_manager.finalized_clusters)} total clusters from {total_reads_processed:,} reads", file=sys.stderr)
        if total_reads_processed == 0:
            sys.exit("No usable reads found after filters")
        groups = cluster_manager.get_finalized_clusters()

    # Collapse & assign 5' truncations with priority to strongest isoform
    isoforms = []
    for (chrom, strand), cluster_list in groups.items():
        for cl in cluster_list:
            mem = cl["members"]
            if not mem: continue

            # Count exact-chain (full-length) support per canon
            exact_counts = Counter(tuple(m["chain_tx"]) for m in mem)
            # Also track total presence of each chain (including truncations) for tiebreaks
            total_chain_counts = exact_counts.copy()

            # Candidate canon chains are unique exact chains seen
            canons = list(exact_counts.keys())

            # Rank canons: (full-length support desc, total support desc, chain length desc)
            canon_rank = sorted(
                canons,
                key=lambda ch: (exact_counts[ch], total_chain_counts[ch], len(ch)),
                reverse=True,
            )

            # Map reads to best-ranked canon where read chain is suffix of the canon chain
            assigned = set()
            chain_to_idxs = defaultdict(list)
            for i, m in enumerate(mem):
                chain_to_idxs[tuple(m["chain_tx"])].append(i)

            for canon in canon_rank:
                idxs = []
                for ch, idxlist in chain_to_idxs.items():
                    if is_suffix(canon, ch):
                        for i in idxlist:
                            if i not in assigned:
                                idxs.append(i)
                if not idxs:
                    continue
                for i in idxs: assigned.add(i)

                members = [mem[i] for i in idxs]
                full_len_members = [m for m in members if tuple(m["chain_tx"]) == canon]
                rep = max(full_len_members or members,
                          key=lambda m: (m["exons"][-1][1]-m["exons"][0][0]))
                rep_exons = list(rep["exons"])
                tes = cl["tes"]
                if strand == "+":
                    if rep_exons[-1][1] != tes:
                        rep_exons[-1] = (rep_exons[-1][0], tes)
                else:
                    if rep_exons[0][0] != tes:
                        rep_exons[0] = (tes, rep_exons[0][1])

                polya_ok = sum(1 for m in members if m["sclen"] >= args.min_polya_length and m["purity"] >= args.min_polya_purity)
                polya_frac = polya_ok/len(members) if members else 0.0

                isoforms.append(dict(
                    chrom=chrom, strand=strand, tes=tes,
                    chain_tx=canon, n_introns=len(canon),
                    members=members, rep_exons=rep_exons,
                    polya_frac=polya_frac,
                ))

    # Global counts & filtering
    total_reads_used = sum(len(iso["members"]) for iso in isoforms)
    kept = []
    metrics_rows = []
    for iso in isoforms:
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
            "|".join(f"{k}:{v}" for k,v in sorted(sample_ct.items())),
            "|".join(f"{k}:{v}" for k,v in sorted(chain_counts.items())),
            int(keep)
        ])

        if keep:
            kept.append(dict(
                chrom=chrom, strand=strand, tes=tes, chain_tx=chain_tx,
                n_introns=n_introns, members=members, rep_exons=rep_exons,
                count=count, frac_global=frac_global, polya_frac=polya_frac,
                med_mapq=med_mapq, med_sclen=med_sclen, med_purity=med_purity,
                sample_ct=sample_ct, chain_counts=chain_counts
            ))

    # Resolve prefix & write metrics
    prefix = args.out_prefix if args.out_prefix else args.out_gtf.replace(".gtf","")
    os.makedirs(os.path.dirname(args.out_gtf) or ".", exist_ok=True)

    metrics_path = f"{prefix}_metrics.tsv"
    with open(metrics_path, "w") as m:
        m.write("#chrom\ttes_1based\tstrand\tintron_chain_tx_order\tn_introns\tread_support\tfrac_global\tpolya_support_frac\tmedian_mapq\tmedian_tail_len\tmedian_tail_purity\tn_unique_chains\tn_full_length_reads\tincluded_trunc\tsample_counts\tchain_counts\tkept\n")
        for row in metrics_rows:
            m.write("\t".join(map(str,row))+"\n")

    if not kept:
        sys.exit(f"No isoforms passed filters. See metrics: {metrics_path}")

    # Annotate and bucket by gene
    def annotate_isoform_safe(iso):
        return annotate_isoform(iso, gtf_by_cs, gtf_tes_sorted, args.tes_match_tol, args.exact_tes_tol) if gtf_txs else dict(
            classification="NOVEL_LOCUS", gene_id="NA", gene_name="NA", matched_tid="NA",
            gtf_tes="NA", gtf_chain_tx=tuple(), tes_delta_bp="NA", exon_overlap_bp=0, match_source="NONE")

    for iso in kept:
        iso["annotation"] = annotate_isoform_safe(iso)

    buckets = defaultdict(list)
    for iso in kept:
        ann = iso["annotation"]
        if ann["gene_name"] != "NA" and ann["gene_id"] != "NA":
            gkey = (ann["gene_name"], ann["gene_id"])
        else:
            gkey = (f"NOVEL_{iso['chrom']}_{iso['strand']}", f"{iso['chrom']}:{iso['strand']}:{iso['tes']}")
        buckets[gkey].append(iso)

    gene_keys_sorted = sorted(buckets.keys(), key=lambda x: (x[0], x[1]))
    gene_index = {gk: i+1 for i, gk in enumerate(gene_keys_sorted)}

    for gk in gene_keys_sorted:
        iso_list = buckets[gk]
        # order by support desc then TES position
        iso_list_sorted = sorted(iso_list, key=lambda x: (-x["count"], x["tes"]))
        for tidx, iso in enumerate(iso_list_sorted, 1):
            gn, gid = gk
            gidx = gene_index[gk]
            iso["gene_name_label"], iso["gene_id_label"] = gn, gid
            iso["gene_index"] = gidx        # ZG
            iso["tx_index"] = tidx          # ZN (per gene)
            iso["zt_label"] = f"{gn}.{gid}.G{gidx}.T{tidx}"

    # Write GTF
    with open(args.out_gtf, "w") as out:
        for iso in sorted(kept, key=lambda x: (x["gene_index"], x["tx_index"])):
            chrom=iso["chrom"]; strand=iso["strand"]; tes=iso["tes"]
            rep_exons=iso["rep_exons"]; chain_tx=iso["chain_tx"]
            t_start, t_end = rep_exons[0][0], rep_exons[-1][1]
            ann = iso["annotation"]
            tid = f"{iso['gene_id_label']}.G{iso['gene_index']}.T{iso['tx_index']}"
            attrs = (
                f'gene_id "{iso["gene_id_label"]}"; '
                f'transcript_id "{tid}"; '
                f'ref_gene_name "{iso["gene_name_label"]}"; '
                f'zt_label "{iso["zt_label"]}"; gene_index "{iso["gene_index"]}"; transcript_index "{iso["tx_index"]}"; '
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
        s.write("#code\tzt_label\tgene_index\ttranscript_index\tchrom\tstrand\tiso_tes\tiso_chain_tx\tgtf_gene_id\tgtf_gene_name\tgtf_transcript_id\tgtf_tes\tgtf_chain_tx\ttes_delta_bp\texon_overlap_bp\tmatch_source\tclassification\tread_support\tfrac_global\tpolya_support_frac\tsample_counts\n")
        for iso in sorted(kept, key=lambda x: (x["gene_index"], x["tx_index"])):
            ann = iso["annotation"]
            s.write("\t".join(map(str,[
                iso["zt_label"], iso["zt_label"], iso["gene_index"], iso["tx_index"], iso["chrom"], iso["strand"], iso["tes"],
                chain_to_str(iso["chain_tx"]),
                ann["gene_id"], ann["gene_name"], ann["matched_tid"],
                ann["gtf_tes"], chain_to_str(ann["gtf_chain_tx"]), ann["tes_delta_bp"], ann["exon_overlap_bp"], ann["match_source"],
                ann["classification"], iso["count"], f"{iso['frac_global']:.3f}", f"{iso['polya_frac']:.4f}",
                "|".join(f"{k}:{v}" for k,v in sorted(iso["sample_ct"].items()))
            ]))+"\n")

    # === Extra outputs: transcript counts, PCA (samples), per-sample stats ===
    import pandas as pd
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # collect samples and transcript IDs
    all_samples = sorted({s for iso in kept for s in iso["sample_ct"].keys()})
    tx_ids = [iso["zt_label"] for iso in kept]

    # counts matrix (KEPT transcripts)
    data = []
    for iso in kept:
        row = [iso["sample_ct"].get(s, 0) for s in all_samples]
        data.append(row)
    counts_df = pd.DataFrame(data, index=tx_ids, columns=all_samples)
    counts_path = f"{prefix}_tx_counts.tsv"
    counts_df.to_csv(counts_path, sep="\t")

    # PCA on samples (features = transcripts), using log1p counts
    X = np.log1p(counts_df.values.T)  # shape: (n_samples, n_transcripts)
    if X.size and X.shape[0] >= 1 and X.shape[1] >= 1:
        Xc = X - X.mean(axis=0, keepdims=True)
        # SVD
        U, S, Vt = np.linalg.svd(Xc, full_matrices=False)
        pcs = min(2, S.size)
        PC = U[:, :pcs] * S[:pcs] if pcs else np.zeros((X.shape[0], 0))
        var_expl = (S**2) / (S**2).sum() if S.size else np.zeros_like(S)

        pca_png = f"{prefix}_tx_counts.pca.png"
        plt.figure(figsize=(6,5))
        if pcs >= 2:
            plt.scatter(PC[:,0], PC[:,1])
            for i, name in enumerate(all_samples):
                plt.text(PC[i,0], PC[i,1], name, fontsize=8, ha="left", va="bottom")
            plt.xlabel(f"PC1 ({var_expl[0]*100:.1f}% var)")
            plt.ylabel(f"PC2 ({var_expl[1]*100:.1f}% var)")
        else:
            # fall back: 1D PC or all zeros
            x = PC[:,0] if pcs >= 1 else np.zeros(X.shape[0])
            plt.scatter(x, np.zeros_like(x))
            for i, name in enumerate(all_samples):
                plt.text(x[i], 0.0, name, fontsize=8, ha="left", va="bottom")
            xl = f"PC1 ({var_expl[0]*100:.1f}% var)" if pcs>=1 else "PC1"
            plt.xlabel(xl); plt.ylabel("PC2")
        plt.title("Sample PCA (log1p transcript counts)")
        plt.tight_layout()
        plt.savefig(pca_png, dpi=150)
        plt.close()
    else:
        pca_png = f"{prefix}_tx_counts.pca.png"
        # write an empty placeholder figure
        plt.figure(figsize=(6,5))
        plt.title("Sample PCA (no data)")
        plt.savefig(pca_png, dpi=150)
        plt.close()

    # per-sample stats
    rows = []
    for s in all_samples:
        col = counts_df[s]
        detected = (col > 0)
        total_reads = int(col.sum())
        n_tx = int(detected.sum())
        med_reads = float(col[detected].median()) if n_tx > 0 else 0.0
        rows.append(dict(sample=s, total_reads=total_reads, n_transcripts=n_tx, median_reads_per_tx=med_reads))
    stats_df = pd.DataFrame(rows).sort_values("sample")
    stats_path = f"{prefix}_per_sample_stats.tsv"
    stats_df.to_csv(stats_path, sep="\t", index=False)

    # Build read -> (ZT, ZG, ZN) map only if needed by downstream tagging/exports
    need_assign = bool(args.write_zt_bams or args.emit_modkit_manifest or args.write_zt_tagged_sample_bams)
    assign = None
    if need_assign:
        assign = defaultdict(dict)
        for iso in kept:
            for m in iso["members"]:
                assign[m["sample"]][m["qname"]] = (iso["zt_label"], iso["gene_index"], iso["tx_index"])  # ZT, ZG, ZN

    # Optional: per-transcript subset BAMs
    if args.write_zt_bams or args.emit_modkit_manifest:
        out_dir = os.path.join(os.path.dirname(args.out_gtf) or ".", "zt_bams")
        os.makedirs(out_dir, exist_ok=True)
        manifest_rows = []

        for iso in kept:
            if sum(iso["sample_ct"].values()) < args.min_total_reads_for_mod:
                continue
            for sample, n in iso["sample_ct"].items():
                if n < args.min_reads_per_sample_for_mod:
                    continue
                iso.setdefault("eligible_samples", set()).add(sample)

        for bam in bams:
            sample = os.path.basename(bam).replace(".bam","")
            elig_codes = {iso["zt_label"] for iso in kept if "eligible_samples" in iso and sample in iso["eligible_samples"]}
            if not elig_codes:
                continue
            with pysam.AlignmentFile(bam, "rb") as inp:
                writers = {}
                try:
                    for code in sorted(elig_codes):
                        out_bam = os.path.join(out_dir, f"{sample}.{code}.bam")
                        writers[code] = pysam.AlignmentFile(out_bam, "wb", header=inp.header)
                    for aln in (inp.fetch(region=args.region) if args.region else inp.fetch()):
                        if aln.is_unmapped: continue
                        if args.primary_only and (aln.is_secondary or aln.is_supplementary): continue
                        if not assign: continue
                        qn = aln.query_name
                        tup = assign.get(sample, {}).get(qn)
                        if not tup:
                            continue
                        zt, zg, zn = tup
                        if zt not in writers:
                            # Guard: this transcript not eligible for this sample
                            continue
                        try:
                            aln.set_tag("ZT", zt, value_type="Z", replace=True)
                        except TypeError:
                            aln.set_tag("ZT", zt)
                        for tag, val in (("ZG", int(zg)), ("ZN", int(zn))):
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
                        # Indexing can be expensive; keep default behavior but allow opting out via CLI
                        if getattr(args, "index_output_bams", True):
                            try:
                                pysam.index(out_bam)
                            except Exception:
                                pass
                        manifest_rows.append([sample, code, out_bam])
        if args.emit_modkit_manifest:
            mani = os.path.join(out_dir, "modkit_manifest.tsv")
            with open(mani, "w") as f:
                f.write("sample\tcode\tbam\n")
                for row in manifest_rows:
                    f.write("\t".join(row)+"\n")
            print(f"[OK] Wrote modkit manifest: {mani}", file=sys.stderr)

    # ZT/ZN-tagged copies of each input BAM (ALL reads) — same count as inputs
    if args.write_zt_tagged_sample_bams:
        out_dir2 = os.path.join(os.path.dirname(args.out_gtf) or ".", "zt_tagged")
        os.makedirs(out_dir2, exist_ok=True)
        for bam in bams:
            sample = os.path.basename(bam).replace(".bam","")
            in_path = bam
            out_path = os.path.join(out_dir2, f"{sample}.zt_tagged.bam")
            with pysam.AlignmentFile(in_path, "rb") as inp, \
                 pysam.AlignmentFile(out_path, "wb", header=inp.header) as outw:
                for aln in (inp.fetch(region=args.region) if args.region else inp.fetch()):
                    if aln.is_unmapped:
                        outw.write(aln); continue
                    if args.primary_only and (aln.is_secondary or aln.is_supplementary):
                        continue
                    qn = aln.query_name
                    tup = assign.get(sample, {}).get(qn) if assign else None
                    if tup:
                        zt, zg, zn = tup
                        try:
                            aln.set_tag("ZT", zt, value_type="Z", replace=True)
                        except TypeError:
                            aln.set_tag("ZT", zt)
                        for tag, val in (("ZG", int(zg)), ("ZN", int(zn))):
                            try:
                                aln.set_tag(tag, val, value_type="i", replace=True)
                            except TypeError:
                                aln.set_tag(tag, val)
                    outw.write(aln)
            if getattr(args, "index_output_bams", True):
                try:
                    pysam.index(out_path)
                except Exception:
                    pass
            print(f"[OK] Wrote ZT/ZN-tagged sample BAM: {out_path}", file=sys.stderr)

    print(f"[OK] Wrote GTF: {args.out_gtf}", file=sys.stderr)
    print(f"[OK] Metrics: {metrics_path}", file=sys.stderr)
    print(f"[OK] Classification summary: {summary_path}", file=sys.stderr)
    print(f"[OK] Counts: {counts_path}", file=sys.stderr)
    print(f"[OK] PCA plot: {pca_png}", file=sys.stderr)
    print(f"[OK] Per-sample stats: {stats_path}", file=sys.stderr)

if __name__ == "__main__":
    main()

