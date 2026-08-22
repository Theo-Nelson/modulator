#!/usr/bin/env python3
"""
Annotate SEQUENCE-based cis-elements on each assembled fragmentform's mature mRNA
and report EVERY modification that falls inside each element -- unbiased across mod
codes (m6A, 5mC, pseudouridine, inosine, and every other code in the aggregate
table). This generalizes the PAS check in check_apa_motifs.py to a family of
motifs anchored to the transcript's 3' end, start codon, stop codon, and 5' end.

Only elements that live in the MATURE (spliced) mRNA are scanned, because that is
what the direct-RNA reads -- and therefore the modification calls -- actually cover.
Intronic elements (branch point, polypyrimidine tract) are deliberately excluded:
there is no modification signal in spliced-out introns.

Elements (all sequence-defined, transcript-oriented, strand-aware):
  PAS          AATAAA + 11 variants, in the window upstream of the cleavage site
  ARE          AU-rich element AUUUA (+ the WWAUUUAWW nonamer), 3'UTR
  CPE          cytoplasmic polyadenylation element UUUUAU, 3'UTR
  GRE          GU-rich element (GU repeats), 3'UTR
  G4           RNA G-quadruplex sequence motif G>=3(N1-7 G>=3)x3, anywhere
  KOZAK        translation-initiation context around the start codon (-6..+4)
  UORF         upstream ORF (AUG..stop) in the 5'UTR
  TOP          5' terminal oligopyrimidine tract at the transcript 5' end
  STOP_CONTEXT stop codon + downstream base (readthrough-prone UGA-C etc.)
  M6AM         cap-adjacent first transcribed A (potential m6Am)

The modification join is purely positional: a mod site (any code) overlaps an
element iff its genomic base is one of the element's bases on the same strand.

Coverage note: direct-RNA reads truncate at the 5' end, so KOZAK / UORF / TOP /
M6AM sit where read (hence modification) coverage is sparsest. Each element row
therefore carries the number of overlapping mod sites AND their coverage, so a
5'-proximal element with no mods reads as "not covered", not "no modification".
"""
from __future__ import annotations
import argparse, json, re, sys
from collections import defaultdict
from pathlib import Path

PAS_CANON = "AATAAA"
PAS_VARIANTS = ["ATTAAA", "TATAAA", "AGTAAA", "AAGAAA", "AATATA", "AATACA",
                "CATAAA", "GATAAA", "AATGAA", "ACTAAA", "AATAGA"]
STOP_CODONS = {"TAA", "TAG", "TGA"}
G4_RE = re.compile(r"G{3,}[ACGTN]{1,7}G{3,}[ACGTN]{1,7}G{3,}[ACGTN]{1,7}G{3,}")
COMP = str.maketrans("ACGTRYKMSWBDHVNacgtrykmswbdhvn", "TGCAYRMKSWVHDBNtgcayrmkswvhdbn")  # full IUPAC


def log(*a):
    print("[seq_elem]", *a, file=sys.stderr, flush=True)


# --------------------------------------------------------------------------- #
#  parsing
# --------------------------------------------------------------------------- #
def _attr(s, key):
    m = re.search(rf'{key} "([^"]+)"', s)
    return m.group(1) if m else None


def load_assembled(gtf_path):
    """zt_label -> dict(chrom, strand, exons=[(s0,e0)], tes, matched_tid, gene_name)."""
    tx = {}
    with open(gtf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            typ, attr = f[2], f[8]
            zt = _attr(attr, "zt_label")
            if not zt:
                continue
            d = tx.setdefault(zt, dict(chrom=f[0], strand=f[6], exons=[],
                                       tes=_attr(attr, "tes"),
                                       matched_tid=_attr(attr, "matched_tid"),
                                       gene_name=_attr(attr, "ref_gene_name") or _attr(attr, "gene_id")))
            if typ == "exon":
                d["exons"].append((int(f[3]) - 1, int(f[4])))
    for d in tx.values():
        d["exons"].sort()
    return tx


def load_ref_codons(ref_gtf):
    """transcript_id -> dict(start=(s0,e0)|None, stop=(s0,e0)|None) genomic 0-based."""
    codons = defaultdict(lambda: dict(start=None, stop=None))
    with open(ref_gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] not in ("start_codon", "stop_codon"):
                continue
            tid = _attr(f[8], "transcript_id")
            if not tid:
                continue
            key = "start" if f[2] == "start_codon" else "stop"
            # GENCODE emits TWO (or three) features for a codon that straddles an exon junction;
            # assigning would let the LAST fragment win and anchor cds_start_t/cds_end_t to the 2nd/3rd
            # base. MERGE to the genomic extent (min start, max end): downstream uses only the codon's
            # 5'-most and 3'-most bases, both of which are real codon bases at these extremes.
            s0, e0 = int(f[3]) - 1, int(f[4])
            prev = codons[tid][key]
            codons[tid][key] = (s0, e0) if prev is None else (min(prev[0], s0), max(prev[1], e0))
    return codons


def load_mod_sites(path):
    """(chrom,strand) -> {pos0: [(mod_code, frac, cov)]}, pooled over sample/ZN."""
    agg = defaultdict(lambda: [0, 0])  # (chrom,pos0,strand,code) -> [Nmod, Ncov]
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(header)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            try:
                ch = f[idx["chrom"]]; pos = int(f[idx["start0"]]); st = f[idx["strand"]]
                code = f[idx["mod_code"]]
                nmod = float(f[idx["Nmod"]]); ncov = float(f[idx["Nvalid_cov"]])
            except (KeyError, ValueError, IndexError):
                continue
            a = agg[(ch, pos, st, code)]
            a[0] += nmod; a[1] += ncov
    out = defaultdict(lambda: defaultdict(list))
    for (ch, pos, st, code), (nmod, ncov) in agg.items():
        frac = (nmod / ncov) if ncov else 0.0
        out[(ch, st)][pos].append((code, round(frac, 4), int(ncov)))
    return out


# --------------------------------------------------------------------------- #
#  mature-mRNA construction + coordinate mapping
# --------------------------------------------------------------------------- #
def build_mrna(fa, chrom, strand, exons):
    """Return (mrna_seq_sense_5to3, tpos_to_g list). T stands for U."""
    seq = []
    tpos_to_g = []
    if strand == "+":
        for (s0, e0) in exons:
            block = fa.fetch(chrom, s0, e0).upper()
            seq.append(block)
            tpos_to_g.extend(range(s0, e0))
    else:
        for (s0, e0) in reversed(exons):
            block = fa.fetch(chrom, s0, e0).upper().translate(COMP)[::-1]
            seq.append(block)
            tpos_to_g.extend(range(e0 - 1, s0 - 1, -1))
    return "".join(seq), tpos_to_g


def genomic_index_of(tpos_to_g, gpos):
    """transcript index whose genomic base == gpos (or None)."""
    # linear scan is fine (transcripts are short); build a dict once per transcript instead
    return tpos_to_g.get(gpos) if isinstance(tpos_to_g, dict) else None


# --------------------------------------------------------------------------- #
#  element scanners  ->  list of (subclass, t_start, t_end, extra)
# --------------------------------------------------------------------------- #
def scan_pas(mrna, pas_window):
    hits = []
    L = len(mrna)
    region = mrna[max(0, L - pas_window):]
    off = max(0, L - pas_window)
    for motif in [PAS_CANON] + PAS_VARIANTS:
        for m in re.finditer(motif, region):
            t0 = off + m.start()
            dist = L - (t0 + 6)  # nt from hexamer 3' end to the transcript 3' end
            sub = "canonical" if motif == PAS_CANON else "variant"
            hits.append((f"{sub}:{motif}", t0, t0 + 6, dict(distance_to_3p=dist)))
    return hits


def _scan_regex(seq, pat, lo, hi, subclass):
    hits = []
    for m in re.finditer(pat, seq[lo:hi]):
        hits.append((subclass, lo + m.start(), lo + m.end(), {}))
    return hits


def scan_are(mrna, lo, hi):
    return _scan_regex(mrna, r"ATTTA", lo, hi, "AUUUA")


def scan_cpe(mrna, lo, hi):
    return _scan_regex(mrna, r"TTTTA[AT]", lo, hi, "UUUUAW")


def scan_gre(mrna, lo, hi):
    return _scan_regex(mrna, r"(?:TG){3,}", lo, hi, "GU_repeat")


def scan_g4(mrna, lo, hi):
    hits = []
    for m in G4_RE.finditer(mrna[lo:hi]):
        hits.append(("rG4", lo + m.start(), lo + m.end(), {}))
    return hits


def scan_kozak(mrna, cds_start_t):
    if cds_start_t is None or cds_start_t < 6 or cds_start_t + 4 > len(mrna):
        return []
    win = mrna[cds_start_t - 6:cds_start_t + 4]  # -6..+4 (A of AUG = +1 at index 6)
    m3 = mrna[cds_start_t - 3]      # -3 position
    p4 = mrna[cds_start_t + 3]      # +4 position
    strong = (m3 in "AG") and (p4 == "G")
    partial = (m3 in "AG") or (p4 == "G")
    sub = "strong" if strong else "adequate" if partial else "weak"
    return [(f"kozak_{sub}", cds_start_t - 6, cds_start_t + 4,
             dict(minus3=m3, plus4=p4, context=win))]


def scan_uorf(mrna, cds_start_t):
    if cds_start_t is None or cds_start_t < 3:
        return []
    hits = []
    utr5 = mrna[:cds_start_t]
    for m in re.finditer("ATG", utr5):
        s = m.start()
        # in-frame stop within the 5'UTR (upstream ORF fully in UTR) or extending toward CDS
        stop_t = None
        for j in range(s, len(mrna) - 2, 3):
            if mrna[j:j + 3] in STOP_CODONS:
                stop_t = j + 3
                break
        kind = "uORF" if (stop_t is not None and stop_t <= cds_start_t) else "uORF_overlapping"
        hits.append((kind, s, (stop_t if stop_t else s + 3), dict(atg_t=s)))
    return hits


def scan_top(mrna):
    # 5' TOP: C at +1 followed by a pyrimidine (C/T) tract
    if not mrna or mrna[0] != "C":
        return []
    i = 1
    while i < len(mrna) and mrna[i] in "CT":
        i += 1
    if i >= 5:  # C + >=4 pyrimidines
        return [("TOP", 0, i, dict(tract_len=i))]
    return []


def scan_m6am(mrna):
    return [("cap_adjacent_A", 0, 1, {})] if mrna[:1] == "A" else []


def scan_stop(mrna, cds_end_t):
    if cds_end_t is None or cds_end_t < 3 or cds_end_t > len(mrna):
        return []
    stop = mrna[cds_end_t - 3:cds_end_t]
    plus4 = mrna[cds_end_t] if cds_end_t < len(mrna) else "."
    readthrough = "readthrough_prone" if (stop == "TGA" and plus4 == "C") else "canonical"
    return [(f"stop_{stop}_{readthrough}", cds_end_t - 3, cds_end_t,
             dict(stop=stop, plus4=plus4))]


# --------------------------------------------------------------------------- #
#  main
# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--assembled-gtf", required=True)
    ap.add_argument("--reference-fa", required=True)
    ap.add_argument("--reference-gtf", required=True, help="for start/stop codons (Kozak/stop/uORF)")
    ap.add_argument("--mod-sites", required=True, help="aggregate_zn <prefix>_FILTERED_sites_long.tsv")
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--out-summary", default="")
    ap.add_argument("--pas-window", type=int, default=60)
    ap.add_argument("--utr3-window", type=int, default=400,
                    help="if no annotated stop codon, scan 3'UTR-type motifs in the last N nt")
    ap.add_argument("--verbose", action="store_true")
    a = ap.parse_args()

    import pysam
    fa = pysam.FastaFile(a.reference_fa)
    tx = load_assembled(a.assembled_gtf)
    codons = load_ref_codons(a.reference_gtf)
    mods = load_mod_sites(a.mod_sites)
    if a.verbose:
        log(f"{len(tx)} fragmentforms, {sum(len(v) for v in mods.values())} mod sites")

    cols = ["zt_label", "gene_name", "chrom", "strand", "element_type", "element_subclass",
            "elem_gstart0", "elem_gend0", "spans_junction", "matched_seq", "region", "detail",
            "n_mod_sites", "mod_codes", "mods_json"]
    rows = []
    per_type = defaultdict(lambda: [0, 0])       # element_type -> [instances, instances_with_mod]
    per_type_codes = defaultdict(lambda: defaultdict(int))  # element_type -> code -> n

    for zt, d in tx.items():
        chrom, strand, exons = d["chrom"], d["strand"], d["exons"]
        if not exons:
            continue
        if chrom not in fa.references:
            continue
        mrna, tg = build_mrna(fa, chrom, strand, exons)
        g2t = {g: i for i, g in enumerate(tg)}
        L = len(mrna)

        # CDS anchors from the matched reference transcript, mapped into this mRNA
        cds_start_t = cds_end_t = None
        c = codons.get(d["matched_tid"] or "", None)
        if c:
            if c["start"]:
                g5 = c["start"][0] if strand == "+" else c["start"][1] - 1  # 5'-most base of start codon
                cds_start_t = g2t.get(g5)
            if c["stop"]:
                g3 = c["stop"][1] - 1 if strand == "+" else c["stop"][0]     # 3'-most base of stop codon
                st_first = g2t.get(c["stop"][0] if strand == "+" else c["stop"][1] - 1)
                cds_end_t = (st_first + 3) if st_first is not None else None

        utr3_lo = cds_end_t if cds_end_t is not None else max(0, L - a.utr3_window)
        utr3_region = "3UTR" if cds_end_t is not None else "3p_window"

        elements = []
        elements += [("PAS", h, "3UTR") for h in scan_pas(mrna, a.pas_window)]
        elements += [("ARE", h, utr3_region) for h in scan_are(mrna, utr3_lo, L)]
        elements += [("CPE", h, utr3_region) for h in scan_cpe(mrna, utr3_lo, L)]
        elements += [("GRE", h, utr3_region) for h in scan_gre(mrna, utr3_lo, L)]
        elements += [("G4", h, "tx") for h in scan_g4(mrna, 0, L)]
        elements += [("KOZAK", h, "start") for h in scan_kozak(mrna, cds_start_t)]
        elements += [("UORF", h, "5UTR") for h in scan_uorf(mrna, cds_start_t)]
        elements += [("TOP", h, "5p") for h in scan_top(mrna)]
        elements += [("M6AM", h, "cap") for h in scan_m6am(mrna)]
        elements += [("STOP_CONTEXT", h, "stop") for h in scan_stop(mrna, cds_end_t)]

        site_map = mods.get((chrom, strand), {})
        for etype, (subclass, t0, t1, extra), region in elements:
            gpos = [tg[i] for i in range(t0, min(t1, L))]
            gs, ge = (min(gpos), max(gpos) + 1) if gpos else (-1, -1)
            # gs..ge is the min/max genomic span; when the element crosses an exon-exon junction the
            # spliced-out intron sits inside it, so [gs,ge) is NOT a contiguous element locus. Flag it
            # so downstream BED/overlap consumers don't treat the intron as part of the element.
            spans_junction = bool(gpos) and (max(gpos) - min(gpos) + 1) != len(gpos)
            hit_mods = []
            for g in gpos:
                for (code, frac, cov) in site_map.get(g, []):
                    hit_mods.append(dict(code=code, gpos=g, frac=frac, cov=cov))
            codes = sorted({m["code"] for m in hit_mods})
            per_type[etype][0] += 1
            if hit_mods:
                per_type[etype][1] += 1
                for cc in codes:
                    per_type_codes[etype][cc] += 1
            rows.append([
                zt, d["gene_name"], chrom, strand, etype, subclass, gs, ge, int(spans_junction),
                mrna[t0:min(t1, L)], region,
                json.dumps({k: v for k, v in extra.items()}) if extra else "",
                len(hit_mods), ",".join(codes), json.dumps(hit_mods) if hit_mods else "",
            ])

    Path(a.out_tsv).parent.mkdir(parents=True, exist_ok=True)
    with open(a.out_tsv, "w") as o:
        o.write("\t".join(cols) + "\n")
        for r in rows:
            o.write("\t".join(str(x) for x in r) + "\n")
    log(f"wrote {len(rows)} element instances -> {a.out_tsv}")

    if a.out_summary:
        with open(a.out_summary, "w") as o:
            o.write("element_type\tn_instances\tn_with_modification\tfrac_with_mod\tmod_codes_seen\n")
            for et in sorted(per_type):
                n, nm = per_type[et]
                codes = ",".join(f"{c}:{k}" for c, k in sorted(per_type_codes[et].items()))
                o.write(f"{et}\t{n}\t{nm}\t{nm/n if n else 0:.3f}\t{codes}\n")
        log(f"wrote summary -> {a.out_summary}")


if __name__ == "__main__":
    main()
