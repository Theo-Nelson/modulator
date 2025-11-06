#!/usr/bin/env python3
"""
Read-backed TES/APA assembler for ONT Direct RNA with 5' truncation handling,
APA clustering, safe GTF-based annotation, and ZT-tagged per-isoform BAMs.

Key features
------------
• TES/APA clustering: reads within --apa-window (default 20bp) are one APA site.
  The site's coordinate is the mode (most frequent TES) of its reads.
• 5' truncation handling: collapse reads whose intron chains are suffixes of a
  canonical chain chosen within each TES cluster.
• Poly(A) gating (per-read and per-isoform) is supported.
• GTF annotation (safe): candidates must be same chr/strand AND TES within
  --tes-match-tol; else overlap-based fallback; else NOVEL_LOCUS.
• Classification: EXACT (same intron chain & TES w/in tol), NOVEL_APA
  (same chain, TES shift), NOVEL_CHAIN (different chain in same locus), NOVEL_LOCUS.
• ZT tagging: optional per-sample, per-isoform BAMs + manifest for modkit.

Typical use (your 5-gene subsets):
  python assemble_tes_variants_from_reads_v3.py \
    --dir test_bams/ALCAM_NHSL1_SERAC1_MXD1_RIOK3_reads \
    --glob "*_5genes.bam" \
    --gtf hg38.ncbiRefSeq.gtf \
    --apa-window 20 --tes-match-tol 25 \
    --primary-only --min-mapq 10 \
    --min-introns-read 1 --min-introns 1 \
    --min-reads 10 --min-frac 0.05 \
    --min-polya-length 12 --min-polya-purity 0.5 --polya-support-frac 0.5 \
    --require-softclip3p 12 \
    --write-zt-bams --emit-modkit-manifest \
    --min-reads-per-sample-for-mod 5 --min-total-reads-for-mod 20 \
    --threads 8 \
    --out-gtf fivegenes_readbacked_annot.gtf
"""

import argparse, os, sys, glob, statistics, hashlib, re
from collections import defaultdict, Counter
import pysam

# ---------- small utils ----------

def median(vals):
    v = [x for x in vals if x is not None]
    if not v: return 0
    try:
        return statistics.median(v)
    except statistics.StatisticsError:
        v.sort(); n=len(v)
        return v[n//2] if n%2 else 0.5*(v[n//2-1]+v[n//2])

def get_tx_strand(aln):
    # Prefer minimap2 'ts', then 'XS', else flag
    if aln.has_tag("ts"):
        v = aln.get_tag("ts")
        if v in ("+","-"): return v
    if aln.has_tag("XS"):
        v = aln.get_tag("XS")
        if v in ("+","-"): return v
    return "-" if aln.is_reverse else "+"

def tes_pos1(aln, tx):
    # 1-based TES (3' end of transcript)
    return aln.reference_end if tx=="+" else aln.reference_start+1

def intron_chain_1based(aln):
    """Return tuple of (donor,acceptor) introns in genomic coords, 1-based, from CIGAR N ops."""
    chain = []
    ref = aln.reference_start  # 0-based
    for op, ln in (aln.cigartuples or []):
        if op == 3:  # N (intron)
            donor = ref + 1
            acceptor = ref + ln
            chain.append((donor, acceptor))
            ref += ln
        elif op in (0,2,7,8):  # M/D/=/X consume ref
            ref += ln
        else:
            pass
    return tuple(chain)

def chain_tx_order(chain, strand):
    return chain if strand=="+" else tuple(reversed(chain))

def exon_blocks_from_aln(aln):
    """List of (start,end) 1-based exons from CIGAR, split on N."""
    exons = []
    ref = aln.reference_start
    cur_start = ref
    for op, ln in (aln.cigartuples or []):
        if op == 3:  # intron
            exons.append((cur_start+1, ref))
            ref += ln
            cur_start = ref
        elif op in (0,2,7,8):
            ref += ln
        else:
            pass
    exons.append((cur_start+1, ref))
    return exons

def softclip3p_len_and_seq(aln, tx):
    """Length & sequence of transcript-3' soft-clip."""
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
    """Purity ∈[0,1] for A (tx '+') or T (tx '-') in soft-clip; robust to orientation quirks."""
    if not seq: return 0.0
    s = seq.upper()
    base = "A" if tx=="+" else "T"
    p1 = s.count(base)/len(s)
    comp = "T" if base=="A" else "A"
    p2 = s.count(comp)/len(s)
    return max(p1, p2)

def cluster_positions(sorted_positions, window):
    """Greedy cluster of sorted ints; return clusters -> list of dicts:
       {"positions": [...], "rep": mode_position, "count": total}
       Mode (most frequent pos) is chosen as representative APA/TES.
    """
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
        rep = sorted(c.items(), key=lambda x: (x[1], x[0]))[-1][0]  # mode; tie -> max pos
        out.append({"positions": cl, "rep": rep, "count": sum(c.values())})
    return out

def is_suffix(longer, shorter):
    """True if 'shorter' is a suffix of 'longer' (both tuples)."""
    if len(shorter) > len(longer): return False
    if len(shorter) == 0: return True
    return longer[-len(shorter):] == shorter

def chain_to_str(chain):
    return "." if not chain else ";".join(f"{d}-{a}" for d,a in chain)

def short_code(s):
    return hashlib.sha1(s.encode()).hexdigest()[:8]

# ---------- GTF parsing ----------

_attr_re = re.compile(r'(\S+)\s+"([^"]+)"')

def parse_gtf_attrs(attr_field):
    d = {}
    for m in _attr_re.finditer(attr_field):
        d[m.group(1)] = m.group(2)
    # tolerate RefSeq styles: sometimes gene_name missing, use gene_id fallback
    if "gene_name" not in d and "gene_id" in d: d["gene_name"] = d["gene_id"]
    return d

class GTFTranscript:
    __slots__ = ("chrom","strand","start","end","tes_1based","exons","gene_id","gene_name","transcript_id","chain_tx")
    def __init__(self, chrom, strand, start, end, attrs):
        self.chrom=chrom; self.strand=strand; self.start=start; self.end=end
        self.tes_1based = end if strand=="+" else start
        self.exons = []  # fill later
        self.gene_id = attrs.get("gene_id","NA")
        self.gene_name = attrs.get("gene_name", self.gene_id)
        self.transcript_id = attrs.get("transcript_id", "NA")
        self.chain_tx = tuple()  # fill later (from exons)

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
    with open(gtf_path) as f:
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
    # fill chains
    for tx in txs.values():
        tx.chain_tx = intron_chain_from_exons(tx.exons, tx.strand)
    return list(txs.values())

# ---------- main ----------

def exon_overlap_len(ex1, ex2):
    tot=0
    for s1,e1 in ex1:
        for s2,e2 in ex2:
            lo=max(s1,s2); hi=min(e1,e2)
            if hi>=lo: tot += (hi-lo+1)
    return tot

def annotate_isoform(iso, gtf_txs, tes_match_tol):
    chrom, strand, tes = iso["chrom"], iso["strand"], iso["tes"]
    # candidates: same chr & strand AND TES within tol
    cands = [tx for tx in gtf_txs
             if tx.chrom==chrom and tx.strand==strand and abs(int(tx.tes_1based)-int(tes)) <= tes_match_tol]

    def rank(tx):
        # prefer smallest TES distance, then chain equality, then larger exon overlap, then longer chain
        same_chain = 0 if tx.chain_tx == iso["chain_tx"] else 1
        return (abs(tx.tes_1based - tes), same_chain, -exon_overlap_len(iso["rep_exons"], tx.exons), -len(tx.chain_tx))

    best = min(cands, key=rank) if cands else None

    if best and best.chain_tx == iso["chain_tx"]:
        cls = "EXACT" if abs(best.tes_1based - tes) == 0 else "NOVEL_APA"
        gene_id, gene_name, tid = best.gene_id, best.gene_name, best.transcript_id
    elif best:
        cls = "NOVEL_CHAIN"
        gene_id, gene_name, tid = best.gene_id, best.gene_name, best.transcript_id
    else:
        # overlap fallback
        overlaps = [tx for tx in gtf_txs if tx.chrom==chrom and tx.strand==strand
                    and exon_overlap_len(iso["rep_exons"], tx.exons) > 0]
        if overlaps:
            best_ov = max(overlaps, key=lambda tx: exon_overlap_len(iso["rep_exons"], tx.exons))
            cls = "NOVEL_CHAIN"
            gene_id, gene_name, tid = best_ov.gene_id, best_ov.gene_name, best_ov.transcript_id
        else:
            cls = "NOVEL_LOCUS"
            gene_id = gene_name = tid = "NA"

    return dict(classification=cls, gene_id=gene_id, gene_name=gene_name, matched_tid=tid)

def main():
    ap = argparse.ArgumentParser(description="TES/APA assembler with GTF annotation & ZT tagging")
    # Inputs
    ap.add_argument("--bams", nargs="+", help="Input BAM(s)")
    ap.add_argument("--glob", help="Glob pattern for BAMs (used with --dir)")
    ap.add_argument("--dir", help="Directory of BAMs")
    ap.add_argument("--region", help="Optional region chr:start-end")
    ap.add_argument("--gtf", help="GTF for annotation")
    ap.add_argument("--threads", type=int, default=0)

    # Read filters
    ap.add_argument("--primary-only", action="store_true", help="Keep only primary alignments")
    ap.add_argument("--min-mapq", type=int, default=1)
    ap.add_argument("--min-introns-read", type=int, default=0, help="Drop reads with < this many introns")
    ap.add_argument("--require-softclip3p", type=int, default=0, help="Require ≥ this many 3' soft-clip bases per read")

    # TES clustering and isoform support
    ap.add_argument("--apa-window", type=int, default=20, help="Cluster TES within this window (bp) to define APA sites (mode used as site)")
    ap.add_argument("--tes-window", type=int, default=None, help="Deprecated; use --apa-window")
    ap.add_argument("--min-reads", type=int, default=10, help="Min reads per isoform")
    ap.add_argument("--min-frac", type=float, default=0.05, help="Min fraction of all reads")
    ap.add_argument("--min-introns", type=int, default=0, help="Require ≥ this many introns in canonical chain")

    # Poly(A) evidence
    ap.add_argument("--min-polya-length", type=int, default=12, help="Per-read: min tail length")
    ap.add_argument("--min-polya-purity", type=float, default=0.7, help="Per-read: min A/T purity")
    ap.add_argument("--polya-support-frac", type=float, default=0.6, help="Per-isoform: fraction of reads passing length+purity")

    # Annotation behavior
    ap.add_argument("--tes-match-tol", type=int, default=25, help="Max distance (bp) between isoform TES and GTF transcript TES to consider a match")

    # ZT-tag outputs / modkit prep
    ap.add_argument("--write-zt-bams", action="store_true", help="Write per-sample per-isoform BAMs with ZT tag")
    ap.add_argument("--emit-modkit-manifest", action="store_true", help="Emit a manifest TSV for modkit runs")
    ap.add_argument("--min-reads-per-sample-for-mod", type=int, default=5)
    ap.add_argument("--min-total-reads-for-mod", type=int, default=20)

    # Output labeling
    ap.add_argument("--out-gtf", default="readbacked_annot.gtf")

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

    # Load GTF (optional but recommended)
    gtf_txs = load_gtf(args.gtf) if args.gtf else []
    gtf_loaded = f"[INFO] Loaded {len(gtf_txs)} transcripts from {args.gtf}" if args.gtf else "[INFO] No GTF supplied; all isoforms will be novel-classified"
    print(gtf_loaded, file=sys.stderr)

    apa_window = args.apa_window if args.apa_window is not None else (args.tes_window if args.tes_window else 20)

    # Gather per-read features
    reads = []  # items: dict below
    for bam in bams:
        if not os.path.exists(bam):
            print(f"[WARN] Missing BAM: {bam}", file=sys.stderr); continue
        sample = os.path.basename(bam).replace(".bam","")
        with pysam.AlignmentFile(bam, "rb") as fh:
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
                reads.append(dict(
                    chrom=chrom,
                    strand=tx,
                    tes=tes_pos1(aln, tx),
                    chain=chain, chain_tx=chain_tx,
                    n_introns=len(chain),
                    exons=exon_blocks_from_aln(aln),
                    bam=os.path.basename(bam),
                    sample=sample,
                    qname=aln.query_name,
                    mapq=aln.mapping_quality,
                    sclen=sclen, purity=purity
                ))

    if not reads:
        sys.exit("No usable reads found after filters")

    # Group by (chrom, strand) -> APA clusters
    groups = defaultdict(list)  # key: (chrom,strand) -> list of clusters
    by_cs = defaultdict(list)
    for r in reads:
        by_cs[(r["chrom"], r["strand"])].append(r)

    for (chrom, strand), rlist in by_cs.items():
        positions = sorted(r["tes"] for r in rlist)
        if not positions: 
            continue
        clusters = cluster_positions(positions, apa_window)  # APA clustering (mode TES)
        # assign reads to nearest cluster rep within window
        for cl in clusters:
            rep_pos = cl["rep"]
            members = [r for r in rlist if abs(r["tes"] - rep_pos) <= apa_window]
            groups[(chrom, strand)].append(dict(tes=rep_pos, members=members))

    # Within each APA cluster, collapse 5' differences by suffix-of canonical chain
    isoforms = []  # canonical isoforms across all clusters
    for (chrom, strand), cluster_list in groups.items():
        for cl in cluster_list:
            mem = cl["members"]
            if not mem: continue

            # unique chains (in tx order) with counts
            chain_counts = Counter(tuple(m["chain_tx"]) for m in mem)
            sorted_chains = sorted(chain_counts.keys(),
                                   key=lambda ch: (len(ch), chain_counts[ch]),
                                   reverse=True)

            assigned = set()  # indices in mem
            chain_to_idxs = defaultdict(list)
            for i, m in enumerate(mem):
                chain_to_idxs[tuple(m["chain_tx"])].append(i)

            for canon in sorted_chains:
                # indices whose chain is a suffix of canon (collapse 5' truncation)
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
                # choose representative alignment (prefer full-length chain)
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

                # polyA support within this isoform
                polya_ok = sum(1 for m in members
                               if m["sclen"] >= args.min_polya_length and m["purity"] >= args.min_polya_purity)
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
    for idx, iso in enumerate(isoforms, 1):
        chrom = iso["chrom"]; strand = iso["strand"]
        tes = iso["tes"]; chain_tx = iso["chain_tx"]; n_introns = iso["n_introns"]
        members = iso["members"]; rep_exons = iso["rep_exons"]; polya_frac = iso["polya_frac"]
        count = len(members)
        frac_global = count/total_reads_used if total_reads_used else 0.0
        med_mapq = median([m["mapq"] for m in members])
        med_sclen = median([m["sclen"] for m in members])
        med_purity = median([m["purity"] for m in members])
        sample_ct = Counter(m["sample"] for m in members)

        chain_counts = Counter(chain_to_str(tuple(m["chain_tx"])) for m in members)
        n_unique_chains = len(chain_counts)
        n_full_len_reads = chain_counts.get(chain_to_str(chain_tx), 0)
        included_trunc = 1 if any(len(tuple(m["chain_tx"])) < len(chain_tx) for m in members) else 0

        keep = True
        if count < args.min_reads or frac_global < args.min_frac: keep=False
        if n_introns < args.min_introns: keep=False
        if polya_frac < args.polya_support_frac: keep=False

        metrics_rows.append([
            chrom, tes, strand,
            chain_to_str(chain_tx), n_introns,
            count, f"{frac_global:.4f}",
            f"{polya_frac:.4f}",
            f"{med_mapq:.1f}", f"{med_sclen:.1f}", f"{med_purity:.3f}",
            n_unique_chains, n_full_len_reads, included_trunc,
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

    metrics_path = args.out_gtf.replace(".gtf","") + "_metrics.tsv"
    os.makedirs(os.path.dirname(args.out_gtf) or ".", exist_ok=True)
    with open(metrics_path, "w") as m:
        m.write("#chrom\ttes_1based\tstrand\tintron_chain_tx_order\tn_introns\tread_support\tfrac_global\tpolya_support_frac\tmedian_mapq\tmedian_tail_len\tmedian_tail_purity\tn_unique_chains\tn_full_length_reads\tincluded_trunc\tsample_counts\tchain_counts\tkept\n")
        for row in metrics_rows:
            m.write("\t".join(map(str,row))+"\n")

    if not kept:
        sys.exit(f"No isoforms passed filters. See metrics: {metrics_path}")

    # Annotation
    for iso in kept:
        if gtf_txs:
            ann = annotate_isoform(iso, gtf_txs, args.tes_match_tol)
        else:
            ann = dict(classification="NOVEL_LOCUS", gene_id="NA", gene_name="NA", matched_tid="NA")
        iso["annotation"] = ann
        # stable code: gene_or_novel + hash(chrom,strand,chain_tx,tes)
        base = ann["gene_name"] if ann["gene_name"]!="NA" else f"NOVEL_{iso['chrom']}_{iso['strand']}"
        sig = f"{iso['chrom']}|{iso['strand']}|{chain_to_str(iso['chain_tx'])}|{iso['tes']}"
        iso["code"] = f"{base}:{short_code(sig)}"

    # Write GTF
    with open(args.out_gtf, "w") as out:
        for i, iso in enumerate(sorted(kept, key=lambda x: (-x["count"], x["tes"])), 1):
            chrom=iso["chrom"]; strand=iso["strand"]; tes=iso["tes"]
            rep_exons=iso["rep_exons"]; chain_tx=iso["chain_tx"]
            t_start, t_end = rep_exons[0][0], rep_exons[-1][1]
            ann = iso["annotation"]
            tid = f"{ann['gene_id'] if ann['gene_id']!='NA' else 'NOVEL'}.T{i}"
            attrs = (
                f'gene_id "{ann["gene_id"] if ann["gene_id"]!="NA" else base_gene_from_code(iso["code"])}"; '
                f'transcript_id "{tid}"; '
                f'ref_gene_name "{ann["gene_name"] if ann["gene_name"]!="NA" else base_gene_from_code(iso["code"])}"; '
                f'read_support "{iso["count"]}"; frac_support_global "{iso["frac_global"]:.4f}"; polya_support_frac "{iso["polya_frac"]:.4f}"; '
                f'intron_chain "{chain_to_str(chain_tx)}"; tes "{tes}"; classification "{ann["classification"]}"; '
                f'code "{iso["code"]}"; matched_tid "{ann["matched_tid"]}"; matched_gid "{ann["gene_id"]}";'
            )
            out.write(f"{chrom}\tReadBacked\ttranscript\t{t_start}\t{t_end}\t1000\t{strand}\t.\t{attrs}\n")
            for j,(s,e) in enumerate(rep_exons,1):
                out.write(f"{chrom}\tReadBacked\texon\t{s}\t{e}\t1000\t{strand}\t.\t{attrs} exon_number \"{j}\";\n")

    # ZT-tag BAMs and manifest
    if args.write_zt_bams or args.emit_modkit_manifest:
        # map (sample, qname) -> code
        assign = defaultdict(dict)  # sample -> qname -> code
        for iso in kept:
            code = iso["code"]
            for m in iso["members"]:
                assign[m["sample"]][m["qname"]] = code

        out_dir = os.path.join(os.path.dirname(args.out_gtf) or ".", "zt_bams")
        os.makedirs(out_dir, exist_ok=True)
        manifest_rows = []

        # eligibility by thresholds
        # only build per-sample per-code BAMs that pass thresholds
        for iso in kept:
            if sum(iso["sample_ct"].values()) < args.min_total_reads_for_mod: 
                continue
            for sample, n in iso["sample_ct"].items():
                if n < args.min_reads_per_sample_for_mod: 
                    continue
                iso.setdefault("eligible_samples", set()).add(sample)

        for bam in bams:
            sample = os.path.basename(bam).replace(".bam","")
            elig_codes = {iso["code"] for iso in kept if "eligible_samples" in iso and sample in iso["eligible_samples"]}
            if not elig_codes: 
                continue
            # open one writer per code for this sample
            with pysam.AlignmentFile(bam, "rb") as inp:
                writers = {}
                try:
                    for code in sorted(elig_codes):
                        out_bam = os.path.join(out_dir, f"{sample}.{code}.bam")
                        writers[code] = pysam.AlignmentFile(out_bam, "wb", header=inp.header)
                    # stream once
                    for aln in inp.fetch(region=args.region) if args.region else inp.fetch():
                        if aln.is_unmapped: continue
                        if args.primary_only and (aln.is_secondary or aln.is_supplementary): continue
                        qn = aln.query_name
                        code = assign.get(sample, {}).get(qn)
                        if not code or code not in writers: 
                            continue
                        # set ZT tag
                        try:
                            aln.set_tag("ZT", code, value_type="Z", replace=True)
                        except TypeError:
                            aln.set_tag("ZT", code)
                        writers[code].write(aln)
                finally:
                    for w in writers.values():
                        w.close()
                # index & add to manifest
                for code in sorted(elig_codes):
                    out_bam = os.path.join(out_dir, f"{sample}.{code}.bam")
                    if os.path.exists(out_bam) and os.path.getsize(out_bam) > 0:
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

    print(f"[OK] Wrote {args.out_gtf} and metrics {metrics_path}", file=sys.stderr)

def base_gene_from_code(code):
    # "GENE:abcd1234" -> "GENE"
    return code.split(":",1)[0] if ":" in code else code

if __name__ == "__main__":
    main()

