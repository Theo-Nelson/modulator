#!/usr/bin/env python3
"""
Read-backed TES assembler for ONT Direct RNA with:
  • TES clustering into APA sites (window = --apa-window; site = mode TES)
  • 5' truncation collapse (suffix-of canonical chain)
  • Poly(A) gating per read & per isoform
  • GTF-aware annotation & classification (fixed logic order)
  • Stable transcript "code" + ZT tagging of reads into per-sample BAMs
  • Modkit manifest for per-code analysis

Classification (final order):
  KNOWN_EXACT      : intron chain == annotated AND TES within --tes-match-tol
  KNOWN_COMPAT_5P  : intron chain is suffix of annotated AND TES within tol
  NOVEL_APA        : chain equal/suffix but TES outside tol for all matches
  NOVEL_CHAIN      : different splicing but stranded exonic overlap with a gene
  NOVEL_LOCUS      : no stranded exonic overlap (assigned --novel-prefix:N)

Notes
-----
• “Suffix-of” uses transcript order (5'→3' along transcript).
• TES comparison uses the APA site representative (mode within cluster).
• ZT-tagged BAMs contain only reads assigned to kept isoforms; they are sorted & indexed.
• Manifest lists (code, sample, bam, n_reads_sample, n_reads_total) filtered by
  --min-reads-per-sample-for-mod and --min-total-reads-for-mod.

Dependencies: pysam
"""

import argparse, os, sys, glob, statistics, re, hashlib, tempfile, shutil
from collections import defaultdict, Counter
import pysam

# -------------------- small utils --------------------

def median(vals):
    v = [x for x in vals if x is not None]
    if not v: return 0
    return statistics.median(v)

def get_tx_strand(aln):
    if aln.has_tag("ts"):
        v = aln.get_tag("ts")
        if v in ("+","-"): return v
    if aln.has_tag("XS"):
        v = aln.get_tag("XS")
        if v in ("+","-"): return v
    return "-" if aln.is_reverse else "+"

def tes_pos1(aln, tx):
    # 1-based TES at transcript 3' end
    return aln.reference_end if tx=="+" else aln.reference_start+1

def intron_chain_from_cigar_1based(aln):
    """Return tuple of (donor,acceptor) 1-based for each 'N' in CIGAR."""
    chain = []
    ref = aln.reference_start  # 0-based
    for op, ln in (aln.cigartuples or []):
        if op == 3:  # N -> intron
            donor = ref + 1
            acceptor = ref + ln
            chain.append((donor, acceptor))
            ref += ln
        elif op in (0,2,7,8):  # M/D/=/X consume ref
            ref += ln
        else:
            # I/S/H/P don't consume ref
            pass
    return tuple(chain)

def chain_tx_order(chain, strand):
    return chain if strand == "+" else tuple(reversed(chain))

def exon_blocks_from_aln(aln):
    """List of (start,end) 1-based exons from CIGAR, split on 'N'."""
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

def cluster_positions_by_window(sorted_positions, window):
    """Greedy cluster of sorted ints by <= window gaps; representative = mode."""
    if not sorted_positions: return []
    clusters, cur = [], [sorted_positions[0]]
    for p in sorted_positions[1:]:
        if p - cur[-1] <= window:
            cur.append(p)
        else:
            clusters.append(cur); cur = [p]
    clusters.append(cur)
    reps = []
    for cl in clusters:
        c = Counter(cl)
        # representative TES = mode (most frequent); tie → max position
        rep = sorted(c.items(), key=lambda x: (x[1], x[0]))[-1][0]
        reps.append((rep, sum(c.values()), cl))
    return reps  # list of tuples (rep_pos, total_count, member_positions)

def is_suffix(longer, shorter):
    if len(shorter) > len(longer): return False
    if len(shorter) == 0: return True
    return longer[-len(shorter):] == shorter

def chain_to_str(chain):
    return "." if not chain else ";".join(f"{d}-{a}" for d,a in chain)

def stable_code(prefix, chrom, strand, tes, chain_tx):
    h = hashlib.md5(f"{chrom}|{strand}|{tes}|{chain_to_str(chain_tx)}".encode()).hexdigest()[:8]
    return f"{prefix}:{h}"

def intervals_overlap_len(a, b):
    """Total overlap length between two exon lists [(s,e),...] 1-based inclusive."""
    tot = 0
    ia = ib = 0
    a = sorted(a); b = sorted(b)
    while ia < len(a) and ib < len(b):
        sa, ea = a[ia]; sb, eb = b[ib]
        s = max(sa, sb); e = min(ea, eb)
        if s <= e:
            tot += (e - s + 1)
        if ea < eb:
            ia += 1
        else:
            ib += 1
    return tot

# -------------------- GTF parsing --------------------

def parse_gtf(gtf_path):
    """
    Return:
      transcripts: tid -> dict(chrom,strand,gene_id,gene_name,exons:list[(s,e)], tes:int, chain_tx:tuple)
      genes: gid -> dict(gene_name, chrom, strand, span:(min,max), exons_union:list[(s,e)], tids:list)
      genes_by_chr_strand: (chrom,strand) -> list of gid
    """
    transcripts = {}
    gene_exons = defaultdict(list)  # gid -> list of exon intervals
    gene_span = {}
    gene_name = {}
    tid_to_gid = {}

    attr_re = re.compile(r'(\w+)\s+"([^"]+)"')
    with open(gtf_path, "r") as f:
        for line in f:
            if not line or line.startswith("#"): continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9: continue
            chrom, source, feature, start, end, score, strand, frame, attrs = parts
            if feature != "exon" and feature != "transcript":
                # We only need exon for chains/TES; transcript lines help when genes without exons (rare)
                pass
            m = dict(attr_re.findall(attrs))
            gid = m.get("gene_id", None)
            gname = m.get("gene_name", gid or "NA")
            tid = m.get("transcript_id", None)

            s = int(start); e = int(end)
            if feature == "exon" and tid and gid:
                # collect per-transcript
                if tid not in transcripts:
                    transcripts[tid] = dict(chrom=chrom, strand=strand, gene_id=gid,
                                            gene_name=gname, exons=[])
                transcripts[tid]["exons"].append((s,e))
                # collect per-gene
                gene_exons[gid].append((s,e))
                gene_name[gid] = gname
                if gid not in gene_span:
                    gene_span[gid] = [s,e]
                else:
                    gene_span[gid][0] = min(gene_span[gid][0], s)
                    gene_span[gid][1] = max(gene_span[gid][1], e)
                tid_to_gid[tid] = gid

    # finalize transcripts (TES + chain)
    for tid, t in list(transcripts.items()):
        ex = sorted(t["exons"], key=lambda x: x[0])
        if not ex:
            del transcripts[tid]
            continue
        strand = t["strand"]
        # TES
        tes = ex[-1][1] if strand == "+" else ex[0][0]
        # chain (from adjacent exons)
        chain = []
        for i in range(len(ex)-1):
            left = ex[i][1]
            right = ex[i+1][0]
            donor = left + 1
            acceptor = right - 1
            chain.append((donor, acceptor))
        chain = tuple(chain)
        chain_tx = chain if strand == "+" else tuple(reversed(chain))
        t.update(exons=ex, tes=tes, chain_tx=chain_tx)

    # finalize genes: merge exon intervals (simple merge)
    genes = {}
    for gid, exs in gene_exons.items():
        exs = sorted(exs)
        merged = []
        for s,e in exs:
            if not merged or s > merged[-1][1] + 1:
                merged.append([s,e])
            else:
                merged[-1][1] = max(merged[-1][1], e)
        merged = [tuple(x) for x in merged]
        # find chrom/strand from any transcript for this gene
        chrom = None; strand=None
        for tid, t in transcripts.items():
            if t["gene_id"] == gid:
                chrom = t["chrom"]; strand = t["strand"]; break
        if chrom is None:
            # fallback: skip genes without transcripts/exons (rare)
            continue
        span = (gene_span[gid][0], gene_span[gid][1])
        tids = [tid for tid, t in transcripts.items() if t["gene_id"] == gid]
        genes[gid] = dict(gene_name=gene_name.get(gid, gid), chrom=chrom, strand=strand,
                          span=span, exons_union=merged, tids=tids)

    genes_by_chr_strand = defaultdict(list)
    for gid, g in genes.items():
        genes_by_chr_strand[(g["chrom"], g["strand"])].append(gid)

    return transcripts, genes, genes_by_chr_strand

# -------------------- core --------------------

def classify_isoform(chrom, strand, tes, exons, chain_tx, transcripts, genes, genes_by_chr_strand, tes_tol, novel_prefix):
    """
    Returns dict with:
      class, matched_tid, matched_gid, gene_id, gene_name
    """
    # candidate genes on same chrom/strand whose span overlaps isoform span
    span_iso = (min(s for s,_ in exons), max(e for _,e in exons))
    cand_gids = []
    for gid in genes_by_chr_strand.get((chrom, strand), []):
        g = genes[gid]
        s1,e1 = span_iso; s2,e2 = g["span"]
        if not (e1 < s2 or e2 < s1):  # span overlap
            # require exonic overlap (not just span)
            if intervals_overlap_len(exons, g["exons_union"]) > 0:
                cand_gids.append(gid)

    # gather overlapping transcripts
    overlap_tids = []
    for gid in cand_gids:
        overlap_tids.extend(genes[gid]["tids"])

    # build equal/suffix sets among overlapping transcripts
    equal = []
    suffix = []
    for tid in overlap_tids:
        t = transcripts[tid]
        if t["chain_tx"] == chain_tx:
            equal.append(tid)
        elif is_suffix(t["chain_tx"], chain_tx):  # annotated is longer, isoform could be 5' truncated
            suffix.append(tid)

    # helper: closest TES within a set of tids
    def tes_match(tids):
        best = None
        for tid in tids:
            t = transcripts[tid]
            diff = abs(t["tes"] - tes)
            if diff <= tes_tol:
                if (best is None) or (diff < best[0]): best = (diff, tid)
        return best  # (diff, tid) or None

    # 1) equal chain
    if equal:
        m = tes_match(equal)
        if m is not None:
            tid = m[1]; gid = transcripts[tid]["gene_id"]; gname = transcripts[tid]["gene_name"]
            return dict(classification="KNOWN_EXACT", matched_tid=tid, matched_gid=gid,
                        gene_id=gid, gene_name=gname)
        else:
            # pick closest TES to report as matched_tid
            tid = min(equal, key=lambda x: abs(transcripts[x]["tes"] - tes))
            gid = transcripts[tid]["gene_id"]; gname = transcripts[tid]["gene_name"]
            return dict(classification="NOVEL_APA", matched_tid=tid, matched_gid=gid,
                        gene_id=gid, gene_name=gname)

    # 2) suffix-of chain
    if suffix:
        m = tes_match(suffix)
        if m is not None:
            tid = m[1]; gid = transcripts[tid]["gene_id"]; gname = transcripts[tid]["gene_name"]
            return dict(classification="KNOWN_COMPAT_5P", matched_tid=tid, matched_gid=gid,
                        gene_id=gid, gene_name=gname)
        else:
            tid = min(suffix, key=lambda x: abs(transcripts[x]["tes"] - tes))
            gid = transcripts[tid]["gene_id"]; gname = transcripts[tid]["gene_name"]
            return dict(classification="NOVEL_APA", matched_tid=tid, matched_gid=gid,
                        gene_id=gid, gene_name=gname)

    # 3) different chain but overlapping gene
    if cand_gids:
        # choose the gene with max exonic overlap
        best_gid = max(cand_gids, key=lambda gid: intervals_overlap_len(exons, genes[gid]["exons_union"]))
        g = genes[best_gid]
        return dict(classification="NOVEL_CHAIN", matched_tid=".", matched_gid=best_gid,
                    gene_id=best_gid, gene_name=g["gene_name"])

    # 4) No overlap -> novel locus
    return dict(classification="NOVEL_LOCUS", matched_tid=".", matched_gid=".",
                gene_id=None, gene_name=None)

def write_gtf(isoforms, out_gtf):
    with open(out_gtf, "w") as out:
        for i, iso in enumerate(sorted(isoforms, key=lambda x: (-x["count"], x["tes"])), 1):
            chrom=iso["chrom"]; strand=iso["strand"]; tes=iso["tes"]
            exons=iso["rep_exons"]
            t_start, t_end = exons[0][0], exons[-1][1]
            attrs = (
                f'gene_id "{iso["gene_id_out"]}"; transcript_id "{iso["tid_out"]}"; '
                f'ref_gene_name "{iso["gene_name_out"]}"; read_support "{iso["count"]}"; '
                f'frac_support_global "{iso["frac_global"]:.4f}"; polya_support_frac "{iso["polya_frac"]:.4f}"; '
                f'intron_chain "{chain_to_str(iso["chain_tx"])}"; tes "{tes}"; '
                f'classification "{iso["class"]}"; code "{iso["code"]}"; '
                f'matched_tid "{iso["matched_tid"]}"; matched_gid "{iso["matched_gid"]}";'
            )
            out.write(f"{chrom}\tReadBacked\ttranscript\t{t_start}\t{t_end}\t1000\t{strand}\t.\t{attrs}\n")
            for j,(s,e) in enumerate(exons,1):
                out.write(f"{chrom}\tReadBacked\texon\t{s}\t{e}\t1000\t{strand}\t.\t{attrs} exon_number \"{j}\";\n")

def sort_and_index_bam(in_bam, out_bam, threads=1):
    tmp_sorted = out_bam
    pysam.sort("-@", str(threads), "-o", tmp_sorted, in_bam)
    pysam.index(tmp_sorted)

# -------------------- main --------------------

def main():
    ap = argparse.ArgumentParser(description="Assemble TES isoforms, annotate with GTF, tag ZT, emit modkit manifest.")
    # Inputs
    ap.add_argument("--bams", nargs="+")
    ap.add_argument("--glob")
    ap.add_argument("--dir")
    ap.add_argument("--region")
    ap.add_argument("--gtf", required=True)

    # Read-level filters
    ap.add_argument("--primary-only", action="store_true")
    ap.add_argument("--min-mapq", type=int, default=1)
    ap.add_argument("--min-introns-read", type=int, default=0)
    ap.add_argument("--require-softclip3p", type=int, default=0)

    # TES clustering & isoform filters
    ap.add_argument("--apa-window", type=int, default=20, help="Cluster TES within this bp window; site = mode TES")
    ap.add_argument("--tes-match-tol", type=int, default=25, help="Tolerance to count an APA site as matching annotation TES")
    ap.add_argument("--min-reads", type=int, default=10, help="Min reads per isoform")
    ap.add_argument("--min-frac", type=float, default=0.05, help="Min fraction of all reads")
    ap.add_argument("--min-introns", type=int, default=0, help="Require >= this many introns in canonical chain")

    # Poly(A) gating
    ap.add_argument("--min-polya-length", type=int, default=12)
    ap.add_argument("--min-polya-purity", type=float, default=0.5)
    ap.add_argument("--polya-support-frac", type=float, default=0.5)

    # Output / tagging
    ap.add_argument("--novel-prefix", default="NOVEL")
    ap.add_argument("--write-zt-bams", action="store_true")
    ap.add_argument("--emit-modkit-manifest", action="store_true")
    ap.add_argument("--min-reads-per-sample-for-mod", type=int, default=5)
    ap.add_argument("--min-total-reads-for-mod", type=int, default=20)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--out-gtf", default="readbacked_annot.gtf")

    args = ap.parse_args()

    # Collect BAMs
    bams = []
    if args.bams: bams += args.bams
    if args.dir and args.glob:
        bams += glob.glob(os.path.join(args.dir, args.glob))
    elif args.dir and not args.glob:
        bams += glob.glob(os.path.join(args.dir, "*.bam"))
    bams = sorted(set(bams))
    if not bams:
        sys.exit("No BAMs found")

    # Parse GTF
    transcripts, genes, genes_by_chr_strand = parse_gtf(args.gtf)

    # Gather per-read features
    reads = []  # items: dict(... below ...)
    total_read_counter = 0
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
                chain = intron_chain_from_cigar_1based(aln)
                if len(chain) < args.min_introns_read: continue
                chain_tx = chain_tx_order(chain, tx)
                exons = exon_blocks_from_aln(aln)
                reads.append(dict(
                    chrom=fh.get_reference_name(aln.reference_id),
                    strand=tx,
                    tes=tes_pos1(aln, tx),
                    chain=chain, chain_tx=chain_tx,
                    n_introns=len(chain),
                    exons=exons,
                    bam=bam, sample=sample,
                    qname=aln.query_name,
                    mapq=aln.mapping_quality,
                    sclen=sclen, purity=purity
                ))
                total_read_counter += 1

    if not reads:
        sys.exit("No usable reads found after filters")

    # Group by (chrom, strand) then APA TES clusters
    groups = []  # list of dict(chrom,strand,tes_rep,members)
    by_cs = defaultdict(list)
    for r in reads:
        by_cs[(r["chrom"], r["strand"])].append(r)

    for (chrom, strand), rlist in by_cs.items():
        positions = sorted(r["tes"] for r in rlist)
        reps = cluster_positions_by_window(positions, args.apa_window)
        for rep_pos, _cnt, members_positions in reps:
            members = [r for r in rlist if abs(r["tes"] - rep_pos) <= args.apa_window]
            groups.append(dict(chrom=chrom, strand=strand, tes=rep_pos, members=members))

    # Within each APA cluster, collapse 5' differences by suffix-of canonical chain
    isoforms = []
    for cl in groups:
        mem = cl["members"]
        if not mem: continue
        # unique chains (in tx order) with counts
        chain_counts = Counter(tuple(m["chain_tx"]) for m in mem)
        sorted_chains = sorted(chain_counts.keys(),
                               key=lambda ch: (len(ch), chain_counts[ch]),
                               reverse=True)
        assigned = set()
        chain_to_idxs = defaultdict(list)
        for i, m in enumerate(mem):
            chain_to_idxs[tuple(m["chain_tx"])].append(i)

        for canon in sorted_chains:
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

            # representative with the longest genomic span
            rep = max(members, key=lambda m: (m["exons"][-1][1]-m["exons"][0][0]))
            rep_exons = list(rep["exons"])
            tes = cl["tes"]
            if cl["strand"] == "+":
                if rep_exons[-1][1] != tes:
                    rep_exons[-1] = (rep_exons[-1][0], tes)
            else:
                if rep_exons[0][0] != tes:
                    rep_exons[0] = (tes, rep_exons[0][1])

            polya_ok = sum(1 for m in members
                           if m["sclen"] >= args.min_polya_length and m["purity"] >= args.min_polya_purity)
            polya_frac = polya_ok/len(members) if members else 0.0

            isoforms.append(dict(
                chrom=cl["chrom"], strand=cl["strand"], tes=tes,
                chain_tx=canon, n_introns=len(canon),
                members=members, rep_exons=rep_exons,
                polya_frac=polya_frac
            ))

    # Global counts & filtering
    total_reads_used = sum(len(iso["members"]) for iso in isoforms)
    kept = []
    for idx, iso in enumerate(isoforms, 1):
        count = len(iso["members"])
        frac_global = count/total_reads_used if total_reads_used else 0.0
        if count < args.min_reads: continue
        if frac_global < args.min_frac: pass  # allow 0 if user set it
        if iso["n_introns"] < args.min_introns: continue
        if iso["polya_frac"] < args.polya_support_frac: continue
        iso["count"] = count
        iso["frac_global"] = frac_global
        kept.append(iso)

    if not kept:
        sys.exit("No isoforms passed filters; try relaxing thresholds.")

    # Classification against GTF and code assignment
    for i, iso in enumerate(kept, 1):
        chrom=iso["chrom"]; strand=iso["strand"]; tes=iso["tes"]
        exons=iso["rep_exons"]; chain_tx=iso["chain_tx"]
        cls = classify_isoform(chrom, strand, tes, exons, chain_tx,
                               transcripts, genes, genes_by_chr_strand,
                               args.tes_match_tol, args.novel_prefix)
        iso["class"] = cls["classification"]
        iso["matched_tid"] = cls["matched_tid"]
        iso["matched_gid"] = cls["matched_gid"]

        # Output gene labeling
        if iso["class"] == "NOVEL_LOCUS":
            # create a synthetic locus name
            iso["gene_name_out"] = f"{args.novel_prefix}"
            iso["gene_id_out"]   = f"{args.novel_prefix}_{i}"
            code_prefix = iso["gene_name_out"]
        else:
            gname = cls["gene_name"] or (transcripts[cls["matched_tid"]]["gene_name"]
                                         if cls["matched_tid"] in transcripts else args.novel_prefix)
            gid = cls["matched_gid"] if cls["matched_gid"] != "." else (transcripts[cls["matched_tid"]]["gene_id"]
                                                                         if cls["matched_tid"] in transcripts else f"{args.novel_prefix}_{i}")
            iso["gene_name_out"] = gname
            iso["gene_id_out"] = gid
            code_prefix = gname

        iso["tid_out"] = f'{iso["gene_name_out"]}.T{i}'
        iso["code"] = stable_code(code_prefix, chrom, strand, tes, chain_tx)

        # per-sample counts
        iso["sample_ct"] = Counter(m["sample"] for m in iso["members"])

    # Write GTF
    out_gtf = args.out_gtf
    write_gtf(kept, out_gtf)
    print(f"[OK] Wrote {out_gtf} (n={len(kept)} isoforms)")

    # Build read→code map per sample for ZT tagging (only reads from kept isoforms)
    sample_read_to_code = defaultdict(dict)  # sample -> {qname: code}
    code_totals = Counter()
    for iso in kept:
        code = iso["code"]
        for m in iso["members"]:
            sample_read_to_code[m["sample"]][m["qname"]] = code
        code_totals[code] += len(iso["members"])

    # Write ZT-tagged BAMs (one per input BAM), sorted & indexed
    zt_bams = {}  # sample -> path
    if args.write_zt_bams:
        for bam in bams:
            sample = os.path.basename(bam).replace(".bam","")
            if sample not in sample_read_to_code:
                print(f"[INFO] No assigned reads for sample {sample}; skipping ZT BAM.")
                continue
            out_unsorted = os.path.splitext(bam)[0] + ".ZT.unsorted.bam"
            out_sorted   = os.path.splitext(bam)[0] + ".ZT.bam"
            with pysam.AlignmentFile(bam, "rb") as fin, \
                 pysam.AlignmentFile(out_unsorted, "wb", template=fin) as fout:
                # add PG line? (optional)
                for aln in fin.fetch():
                    qn = aln.query_name
                    code = sample_read_to_code[sample].get(qn)
                    if code is None: 
                        continue  # write only reads that were assigned to an isoform
                    aln.set_tag("ZT", code, value_type="Z")
                    fout.write(aln)
            # sort & index
            sort_and_index_bam(out_unsorted, out_sorted, threads=args.threads)
            os.remove(out_unsorted)
            zt_bams[sample] = out_sorted
            print(f"[OK] Wrote {out_sorted} (+ .bai)")

    # Emit modkit manifest
    if args.emit_modkit_manifest:
        manifest_path = os.path.splitext(out_gtf)[0] + "_modkit_manifest.tsv"
        with open(manifest_path, "w") as mf:
            mf.write("code\tsample\tzt_bam\tn_reads_sample\tn_reads_total\n")
            # per code, per sample counts from kept isoforms
            per_code_per_sample = defaultdict(Counter)
            for iso in kept:
                code = iso["code"]
                for samp, n in iso["sample_ct"].items():
                    per_code_per_sample[code][samp] += n
            for code, samp_ct in per_code_per_sample.items():
                total = sum(samp_ct.values())
                if total < args.min_total_reads_for_mod:
                    continue
                for samp, n in samp_ct.items():
                    if n < args.min_reads_per_sample_for_mod:
                        continue
                    ztbam = zt_bams.get(samp, os.path.splitext([b for b in bams if os.path.basename(b).startswith(samp)][0])[0] + ".ZT.bam")
                    mf.write(f"{code}\t{samp}\t{ztbam}\t{n}\t{total}\n")
        print(f"[OK] Wrote modkit manifest: {manifest_path}")

if __name__ == "__main__":
    main()

