#!/usr/bin/env python3
"""
Read-backed TES assembler (ONT DRS) with:
 - 5' truncation collapse (unchanged from v1)
 - GTF-based classification (3' match within tolerance; ignore 5')
 - Gene naming (use GTF gene when possible; otherwise novel locus name)
 - Stable per-transcript "code" assignment
 - Output: annotated GTF, metrics, read→code lists, eligibility table
 - Optional: write ZT tag into BAMs; emit per-code BAM extract & modkit manifest

Usage (MXD1 slice example)
--------------------------
python assemble_tes_variants_from_reads_v2.py \
  --dir test_bams/MXD1_reads \
  --glob "*_MXD1_chr2_69915109_69942945.bam" \
  --region chr2:69915109-69942945 \
  --primary-only --min-mapq 10 \
  --min-introns-read 1 --min-introns 1 \
  --tes-window 12 --min-reads 10 --min-frac 0.05 \
  --min-polya-length 12 --min-polya-purity 0.7 --polya-support-frac 0.7 \
  --require-softclip3p 12 \
  --gene-name MXD1 --gene-id MXD1_custom \
  --gtf hg38.ncbiRefSeq.gtf \
  --tes-match-tol 25 \
  --novel-prefix NOVEL \
  --out-gtf MXD1_readbacked_annot.gtf \
  --write-zt-bams \
  --min-reads-per-sample-for-mod 5 \
  --min-total-reads-for-mod 20 \
  --emit-modkit-manifest
"""

import argparse, os, sys, glob, statistics, re, hashlib, itertools, json
from collections import defaultdict, Counter
import pysam

# ---------- small utils (kept from your v1, with a few additions) ----------

def median(vals):
    v = [x for x in vals if x is not None]
    if not v: return 0
    try:
        return statistics.median(v)
    except statistics.StatisticsError:
        v.sort()
        n = len(v)
        return v[n//2] if n % 2 == 1 else 0.5*(v[n//2-1] + v[n//2])

def get_tx_strand(aln):
    if aln.has_tag("ts"):
        v = aln.get_tag("ts")
        if v in ("+","-"): return v
    if aln.has_tag("XS"):
        v = aln.get_tag("XS")
        if v in ("+","-"): return v
    return "-" if aln.is_reverse else "+"

def tes_pos1(aln, tx):
    return aln.reference_end if tx=="+" else aln.reference_start+1  # 1-based

def intron_chain_1based(aln):
    chain = []
    ref = aln.reference_start  # 0-based
    for op, ln in (aln.cigartuples or []):
        if op == 3:  # N
            donor = ref + 1
            acceptor = ref + ln
            chain.append((donor, acceptor))
            ref += ln
        elif op in (0,2,7,8):
            ref += ln
        else:
            pass
    return tuple(chain)

def chain_tx_order(chain, strand):
    return chain if strand == "+" else tuple(reversed(chain))

def exon_blocks_from_aln(aln):
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

def cluster_positions(sorted_positions, window):
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
        rep = sorted(c.items(), key=lambda x: (x[1], x[0]))[-1][0]  # mode; tie -> max pos
        reps.append((rep, sum(c.values())))
    return reps

def is_suffix(longer, shorter):
    if len(shorter) > len(longer): return False
    if len(shorter) == 0: return True
    return longer[-len(shorter):] == shorter

def chain_to_str(chain):
    return "." if not chain else ";".join(f"{d}-{a}" for d,a in chain)

def chains_compatible_ignoring_5prime(a, b):
    """True if a is suffix of b or b is suffix of a."""
    return is_suffix(a, b) or is_suffix(b, a)

def exons_overlap(ex1, ex2):
    # ex1, ex2 are lists of (s,e), 1-based inclusive
    i, j = 0, 0
    A, B = sorted(ex1), sorted(ex2)
    while i < len(A) and j < len(B):
        s1,e1 = A[i]; s2,e2 = B[j]
        if e1 < s2: i += 1
        elif e2 < s1: j += 1
        else: return True
    return False

# ---------- GTF parsing ----------

def parse_gtf(gtf_path):
    """
    Return:
      transcripts_by_cs[(chrom,strand)] = list of dicts {
        tid, gid, gname, chrom, strand, exons[(s,e)], chain_tx[(d,a)], tes(int)
      }
    """
    attr_re = re.compile(r'(\S+)\s+"([^"]*)"')
    tx_exons = defaultdict(list)
    tx_meta  = {}
    with open(gtf_path) as fh:
        for line in fh:
            if not line or line.startswith("#"): continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9: continue
            chrom, source, feature, start, end, score, strand, frame, attrs = f
            if feature != "exon": continue
            start = int(start); end = int(end)
            ad = dict(attr_re.findall(attrs))
            tid = ad.get("transcript_id")
            gid = ad.get("gene_id", "")
            gname = ad.get("gene_name", gid if gid else "")
            if not tid:  # skip if no transcript_id
                continue
            tx_exons[tid].append((start,end))
            tx_meta[tid] = (chrom, strand, gid, gname)
    transcripts_by_cs = defaultdict(list)
    for tid, exons in tx_exons.items():
        chrom, strand, gid, gname = tx_meta[tid]
        ex = sorted(exons)
        # make intron chain from exons (1-based)
        introns = []
        for i in range(len(ex)-1):
            donor = ex[i][1]     # exon i end
            acceptor = ex[i+1][0]  # exon i+1 start
            introns.append((donor, acceptor))
        chain = tuple(introns)
        chain_tx = chain if strand == "+" else tuple(reversed(chain))
        tes = ex[-1][1] if strand == "+" else ex[0][0]
        transcripts_by_cs[(chrom,strand)].append(dict(
            tid=tid, gid=gid, gname=gname, chrom=chrom, strand=strand,
            exons=ex, chain_tx=chain_tx, tes=tes
        ))
    return transcripts_by_cs

# ---------- code assignment ----------

def code_for_isoform(iso, gene_label, style="short"):
    """
    Deterministic code per isoform across samples.
    base string includes gene, strand, TES, and intron chain (tx order).
    """
    base = f"{gene_label}|{iso['strand']}|TES{iso['tes']}|{chain_to_str(iso['chain_tx'])}"
    if style == "long":
        return base
    h = hashlib.sha1(base.encode()).hexdigest()[:8]
    return f"{gene_label}:{h}"

# ---------- main ----------

def main():
    ap = argparse.ArgumentParser(description="TES assembler + GTF classification + codes + optional ZT tagging")
    # Inputs
    ap.add_argument("--bams", nargs="+", help="Input BAM(s)")
    ap.add_argument("--glob", help="Glob pattern for BAMs (used with --dir)")
    ap.add_argument("--dir", help="Directory of BAMs")
    ap.add_argument("--region", help="Optional region chr:start-end")
    ap.add_argument("--primary-only", action="store_true", help="Keep only primary alignments")
    ap.add_argument("--min-mapq", type=int, default=1)

    # Per-read gating
    ap.add_argument("--min-introns-read", type=int, default=0)
    ap.add_argument("--require-softclip3p", type=int, default=0)

    # TES clustering and isoform support
    ap.add_argument("--tes-window", type=int, default=12)
    ap.add_argument("--min-reads", type=int, default=10)
    ap.add_argument("--min-frac", type=float, default=0.05)
    ap.add_argument("--min-introns", type=int, default=0)

    # Poly(A) evidence
    ap.add_argument("--min-polya-length", type=int, default=12)
    ap.add_argument("--min-polya-purity", type=float, default=0.7)
    ap.add_argument("--polya-support-frac", type=float, default=0.6)

    # GTF + classification
    ap.add_argument("--gtf", help="Reference GTF for classification")
    ap.add_argument("--tes-match-tol", type=int, default=25, help="Abs distance allowed for 3' end match")
    ap.add_argument("--novel-prefix", default="NOVEL", help="Prefix for novel locus gene names")
    ap.add_argument("--code-style", choices=["short","long"], default="short")

    # Output labeling
    ap.add_argument("--gene-name", default="GENE")
    ap.add_argument("--gene-id", default="GENE_custom")
    ap.add_argument("--out-gtf", default="readbacked_annot.gtf")

    # ZT + mod prep
    ap.add_argument("--write-zt-bams", action="store_true", help="Write new *_ZT.bam with per-read ZT tag")
    ap.add_argument("--min-reads-per-sample-for-mod", type=int, default=5)
    ap.add_argument("--min-total-reads-for-mod", type=int, default=20)
    ap.add_argument("--emit-modkit-manifest", action="store_true", help="Emit a shell script with per-code BAM extracts + modkit placeholders")
    ap.add_argument("--threads", type=int, default=4)

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
    samples = [os.path.basename(p).replace(".bam","") for p in bams]
    sample_by_bam = {os.path.basename(p): p for p in bams}

    # Load GTF (optional)
    tx_by_cs = defaultdict(list)
    if args.gtf:
        print(f"[INFO] Loading GTF: {args.gtf}", file=sys.stderr)
        tx_by_cs = parse_gtf(args.gtf)
        print(f"[INFO] Parsed {sum(len(v) for v in tx_by_cs.values())} transcripts from GTF.", file=sys.stderr)

    # Gather per-read features
    reads = []  # dict per read
    for bam in bams:
        if not os.path.exists(bam):
            print(f"[WARN] Missing BAM: {bam}", file=sys.stderr); continue
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
                reads.append(dict(
                    chrom=fh.get_reference_name(aln.reference_id),
                    strand=tx,
                    tes=tes_pos1(aln, tx),
                    chain=chain, chain_tx=chain_tx,
                    n_introns=len(chain),
                    exons=exon_blocks_from_aln(aln),
                    bam=os.path.basename(bam),
                    qname=aln.query_name,
                    mapq=aln.mapping_quality,
                    sclen=sclen, purity=purity
                ))

    if not reads:
        sys.exit("No usable reads found after filters")

    # Group by (chrom, strand) then TES clusters
    groups = defaultdict(list)
    by_cs = defaultdict(list)
    for r in reads:
        by_cs[(r["chrom"], r["strand"])].append(r)

    for (chrom, strand), rlist in by_cs.items():
        positions = sorted(r["tes"] for r in rlist)
        if not positions:
            continue
        reps = cluster_positions(positions, args.tes_window)
        for rep_pos, _cnt in reps:
            members = [r for r in rlist if abs(r["tes"] - rep_pos) <= args.tes_window]
            groups[(chrom, strand)].append(dict(tes=rep_pos, members=members))

    # Within each TES cluster, collapse 5' differences by suffix-of canonical chain
    isoforms = []  # canonical across clusters
    for (chrom, strand), cluster_list in groups.items():
        for cl in cluster_list:
            mem = cl["members"]
            if not mem: continue

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

                polya_ok = sum(1 for m in members
                               if m["sclen"] >= args.min_polya_length and m["purity"] >= args.min_polya_purity)
                polya_frac = polya_ok/len(members) if members else 0.0

                isoforms.append(dict(
                    chrom=chrom, strand=strand, tes=tes,
                    chain_tx=canon, n_introns=len(canon),
                    members=members, rep_exons=rep_exons,
                    polya_frac=polya_frac,
                ))

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
        sample_ct = Counter(m["bam"].replace(".bam","") for m in members)

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
    with open(metrics_path, "w") as m:
        m.write("#chrom\ttes_1based\tstrand\tintron_chain_tx_order\tn_introns\tread_support\tfrac_global\tpolya_support_frac\tmedian_mapq\tmedian_tail_len\tmedian_tail_purity\tn_unique_chains\tn_full_length_reads\tincluded_trunc\tsample_counts\tchain_counts\tkept\n")
        for row in metrics_rows:
            m.write("\t".join(map(str,row))+"\n")

    if not kept:
        sys.exit(f"No isoforms passed filters. See metrics: {metrics_path}")

    # ---------- classification vs GTF ----------
    novel_counter = itertools.count(1)
    annotated = []
    for iso in kept:
        chrom, strand, tes = iso["chrom"], iso["strand"], iso["tes"]
        rep_exons = iso["rep_exons"]
        candidates = tx_by_cs.get((chrom,strand), [])
        # rank matches
        exact_3p = []
        compat = []
        overlap_any = []
        for t in candidates:
            if exons_overlap(rep_exons, t["exons"]):
                overlap_any.append(t)
            # splice comparison ignoring 5'
            if chains_compatible_ignoring_5prime(iso["chain_tx"], t["chain_tx"]):
                compat.append(t)
                if abs(tes - t["tes"]) <= args.tes_match_tol:
                    exact_3p.append(t)
        cls = "NOVEL_LOCUS"
        gene_label = None
        match_tid = match_gid = match_gname = ""
        # pick best
        if exact_3p:
            # prefer exact intron chain equality if available
            exact_chain = [t for t in exact_3p if t["chain_tx"] == iso["chain_tx"]]
            tbest = exact_chain[0] if exact_chain else exact_3p[0]
            match_tid, match_gid, match_gname = tbest["tid"], tbest["gid"], tbest["gname"] or tbest["gid"]
            gene_label = match_gname or match_gid
            cls = "KNOWN_EXACT" if tbest["chain_tx"] == iso["chain_tx"] else "KNOWN_COMPAT_5P"
        elif compat:
            # Same splicing (ignoring 5'), different TES => APA
            tbest = compat[0]
            match_tid, match_gid, match_gname = tbest["tid"], tbest["gid"], tbest["gname"] or tbest["gid"]
            gene_label = match_gname or match_gid
            cls = "NOVEL_APA"
        elif overlap_any:
            # Different splicing within gene space
            tbest = overlap_any[0]
            match_tid, match_gid, match_gname = tbest["tid"], tbest["gid"], tbest["gname"] or tbest["gid"]
            gene_label = match_gname or match_gid
            cls = "NOVEL_SPLICE"
        else:
            # entirely new locus
            gene_label = f"{args.novel_prefix}_G{next(novel_counter):04d}"

        code = code_for_isoform(iso, gene_label, style=args.code-style if hasattr(args, "code-style") else args.code_style)
        iso.update(dict(classification=cls, gene_label=gene_label,
                        gtf_tid=match_tid, gtf_gid=match_gid, gtf_gname=match_gname,
                        code=code))
        annotated.append(iso)

    # ---------- write annotated GTF ----------
    with open(args.out_gtf, "w") as out:
        for i, iso in enumerate(sorted(annotated, key=lambda x: (-x["count"], x["tes"])), 1):
            chrom=iso["chrom"]; strand=iso["strand"]; tes=iso["tes"]
            rep_exons=iso["rep_exons"]; chain_tx=iso["chain_tx"]
            t_start, t_end = rep_exons[0][0], rep_exons[-1][1]
            tid = f'{iso["gene_label"]}.T{i}'
            attrs = (
                f'gene_id "{iso.get("gtf_gid", iso["gene_label"])}"; transcript_id "{tid}"; '
                f'ref_gene_name "{iso.get("gtf_gname", iso["gene_label"])}"; '
                f'read_support "{iso["count"]}"; frac_support_global "{iso["frac_global"]:.4f}"; polya_support_frac "{iso["polya_frac"]:.4f}"; '
                f'intron_chain "{chain_to_str(chain_tx)}"; tes "{tes}"; classification "{iso["classification"]}"; '
                f'code "{iso["code"]}"; matched_tid "{iso.get("gtf_tid","")}"; matched_gid "{iso.get("gtf_gid","")}";'
            )
            out.write(f"{chrom}\tReadBacked\ttranscript\t{t_start}\t{t_end}\t1000\t{strand}\t.\t{attrs}\n")
            for j,(s,e) in enumerate(rep_exons,1):
                out.write(f"{chrom}\tReadBacked\texon\t{s}\t{e}\t1000\t{strand}\t.\t{attrs} exon_number \"{j}\";\n")

    # ---------- read→code maps, eligibility for mod ----------
    # per-sample readname lists per code
    out_prefix = args.out_gtf.replace(".gtf","")
    readlists_dir = out_prefix + "_readlists"
    os.makedirs(readlists_dir, exist_ok=True)

    # counts and lists
    code_info = {}  # code -> dict(gene_label, class, chain, tes)
    counts_per_code_sample = defaultdict(lambda: Counter())
    readnames_per_code_sample = defaultdict(lambda: defaultdict(list))

    for iso in annotated:
        code = iso["code"]
        code_info[code] = dict(gene_label=iso["gene_label"], classification=iso["classification"],
                               chain=chain_to_str(iso["chain_tx"]), tes=iso["tes"])
        for m in iso["members"]:
            samp = m["bam"].replace(".bam","")
            counts_per_code_sample[code][samp] += 1
            readnames_per_code_sample[code][samp].append(m["qname"])

    # write per-code per-sample readname files
    for code, by_sample in readnames_per_code_sample.items():
        for samp, qnames in by_sample.items():
            fn = os.path.join(readlists_dir, f"{samp}__{code}.reads.txt")
            with open(fn, "w") as fh:
                fh.write("\n".join(sorted(set(qnames))) + "\n")

    # eligibility table
    elig_path = out_prefix + "_eligibility_for_mod.tsv"
    with open(elig_path, "w") as fh:
        header = ["code","gene_label","classification","chain_tx","tes"] + samples + ["total","eligible"]
        fh.write("\t".join(header) + "\n")
        for code in sorted(code_info.keys()):
            row = [code,
                   code_info[code]["gene_label"],
                   code_info[code]["classification"],
                   code_info[code]["chain"],
                   str(code_info[code]["tes"])]
            total = 0
            per_sample = []
            for s in samples:
                c = counts_per_code_sample[code][s]
                per_sample.append(str(c)); total += c
            eligible = (total >= args.min_total_reads_for_mod) and all(
                counts_per_code_sample[code][s] >= args.min_reads_per_sample_for_mod
                or counts_per_code_sample[code][s] == 0  # allow zeros in some samples
                for s in samples
            )
            fh.write("\t".join(row + per_sample + [str(total), "1" if eligible else "0"]) + "\n")

    # ---------- optional: add ZT tag to BAMs ----------
    if args.write_zt_bams:
        print("[INFO] Writing ZT-tagged BAMs...", file=sys.stderr)
        # build qname->code per sample
        q2code_by_sample = defaultdict(dict)
        for code, by_sample in readnames_per_code_sample.items():
            for samp, qnames in by_sample.items():
                # if a read appears in multiple codes (shouldn't), last wins
                for q in qnames:
                    q2code_by_sample[samp][q] = code

        for bam_path in bams:
            samp = os.path.basename(bam_path).replace(".bam","")
            mapping = q2code_by_sample.get(samp, {})
            if not mapping:
                print(f"[INFO] No reads to tag for {samp}", file=sys.stderr)
                continue
            outbam = os.path.splitext(bam_path)[0] + "_ZT.bam"
            with pysam.AlignmentFile(bam_path, "rb") as inp, \
                 pysam.AlignmentFile(outbam, "wb", header=inp.header) as outp:
                it = inp.fetch()  # full file; fast enough for test slices
                for aln in it:
                    q = aln.query_name
                    if q in mapping:
                        aln.set_tag("ZT", mapping[q], value_type="Z")
                    outp.write(aln)
            pysam.index(outbam)
            print(f"[INFO] Wrote {outbam}", file=sys.stderr)

    # ---------- optional: emit per-code BAM extraction + modkit manifest ----------
    if args.emit_modkit_manifest:
        sh = out_prefix + "_modkit_manifest.sh"
        with open(sh, "w") as fh:
            fh.write("#!/usr/bin/env bash\nset -euo pipefail\n\n")
            fh.write(f"THREADS={args.threads}\n\n")
            fh.write("# Extract per-code BAMs and (placeholder) run modkit\n")
            for bam_path in bams:
                samp = os.path.basename(bam_path).replace(".bam","")
                fh.write(f"# --- sample: {samp} ---\n")
                for code in sorted(code_info.keys()):
                    readlist = os.path.join(readlists_dir, f"{samp}__{code}.reads.txt")
                    outbam = f"{samp}__{code}.bam"
                    # only extract if readlist exists & non-empty
                    fh.write(f'if [[ -s "{readlist}" ]]; then\n')
                    fh.write(f'  echo "[INFO] {samp} {code}: extracting BAM"\n')
                    fh.write(f'  samtools view -@ ${{THREADS}} -N "{readlist}" -b "{bam_path}" -o "{outbam}"\n')
                    fh.write(f'  samtools index -@ ${{THREADS}} "{outbam}"\n')
                    # placeholder modkit command — adapt as needed
                    fh.write(f'  # modkit pileup --threads ${{THREADS}} --output {samp}__{code}_mods.bed "{outbam}"\n')
                    fh.write("fi\n")
                fh.write("\n")
        os.chmod(sh, 0o755)
        print(f"[OK] Wrote modkit manifest: {sh}", file=sys.stderr)

    print(f"[OK] Wrote {args.out_gtf}, metrics {metrics_path}, readlists dir {readlists_dir}, and eligibility {elig_path}")

if __name__ == "__main__":
    main()

