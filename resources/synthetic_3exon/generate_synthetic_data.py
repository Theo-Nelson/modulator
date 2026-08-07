#!/usr/bin/env python3
"""
Generate a synthetic dRNA-seq dataset with KNOWN GROUND TRUTH to validate every
feature of the `modulator` pipeline end to end.

Two contigs:
  chrSyn  (15 kb) -- GENE_A (the requested 3x500nt-exon / 500nt-intron model,
                     co-dependency + SNP + haplotype + condition structure),
                     GENE_B (tandem APA + PAS + internal priming + poly(A)),
                     GENE_C (minus strand).
  chrSyn2 (12 kb) -- GENE_OV1/GENE_OV2 (overlapping genes -> multigene filter),
                     GENE_TR (5' truncation -> hierarchical_stoich correction),
                     GENE_TA (tandem APA WITH a differential mod -> classify TANDEM_APA).

Deterministic (seeded). Also writes a ground-truth table describing what each
feature SHOULD report, so the outputs can be checked mechanically (validate_outputs.py).

Requires: pysam, samtools on PATH (both in the `modulator` conda env).
"""
from __future__ import annotations
import argparse, array, random, subprocess
from pathlib import Path

# ============================================================================= #
#  Contig 1 : chrSyn  (all coordinates 0-based, half-open [start, end))
# ============================================================================= #
CONTIG = "chrSyn"
CONTIG_LEN = 15000

# GENE_A : 3 exons x 500 nt, introns 500 nt   (the requested structure)
A_e1 = (1000, 1500); A_i1 = (1500, 2000); A_e2 = (2000, 2500)
A_i2 = (2500, 3000); A_e3 = (3000, 3500); A_TES = 3499

# GENE_B : 2 exons, long last exon -> tandem APA (distal vs proximal TES)
B_e1 = (6500, 7000); B_i1 = (7000, 7500); B_e2 = (7500, 8500)
B_TES_distal = 8480        # B1, has canonical PAS
B_TES_prox   = 8000        # B2, NO pas + genomic A-rich downstream -> internal priming

# GENE_C : 3 exons on the MINUS strand
C_e1 = (11000, 11400); C_i1 = (11400, 11800); C_e2 = (11800, 12200)
C_i2 = (12200, 12600); C_e3 = (12600, 13000); C_TES = 11000

# GENE_A modification sites (0-based genomic pos of the modified base)
P = 1080; Q = 1160          # exon1 co-dependent pair
R = 1300; S = 1380          # exon1 independent pair
X = 2100; Y = 2260          # exon2 mutually-exclusive pair (A1 only)
D    = 3050                 # exon3 isoform-differential (A1 hi / A2 lo)
Msnp = 3200                 # exon3 SNP1-controlled
Cmod = 3350                 # exon3 5mC (second mod code)
COND = 3450                 # exon3 condition-differential
SNP1 = 3201; SNP2 = 1250
A_M6A_SITES = [P, Q, R, S, X, Y, D, Msnp, COND]
A_5MC_SITES = [Cmod]

# ============================================================================= #
#  Contig 2 : chrSyn2
# ============================================================================= #
CONTIG2 = "chrSyn2"
CONTIG2_LEN = 12000

# --- overlapping genes (multigene filter): OV1 e3 overlaps OV2 e1 on the SAME strand
OV1_e1 = (500, 900); OV1_e2 = (1400, 1800); OV1_e3 = (2300, 2800); OV1_TES = 2799
OV2_e1 = (2500, 3000); OV2_e2 = (3500, 4000); OV2_TES = 3999      # e1 overlaps OV1_e3 in [2500,2800]

# --- 5' truncation gene (hierarchical_stoich): 4 exons, cassette-skip + truncated reads
TR_e1 = (5500, 5800); TR_e2 = (6300, 6600); TR_e3 = (7100, 7400); TR_e4 = (7900, 8400)
TR_TES = 8399
TR_site = 8100          # m6A in the shared terminal exon e4 (3' of the divergence)
TR_trunc_start = 7200   # 5'-truncated reads start mid-e3 (3' of the e2 divergence ~6600)

# --- tandem-APA gene WITH a differential mod (classify -> TANDEM_APA)
TA_e1 = (10000, 10300)
TA_last_acc = 10600
TA_TES_prox = 10999; TA_TES_dist = 11499
TA_site = 10800         # m6A in the shared part of the last exon (< prox TES)

RNG = random.Random(20260807)


def bern(p): return 1 if RNG.random() < p else 0
def revcomp(s):
    t = str.maketrans("ACGTN", "TGCAN"); return s.translate(t)[::-1]


# ----------------------------------------------------------------------------- #
#  Reference sequences
# ----------------------------------------------------------------------------- #
def build_references() -> dict[str, str]:
    seq1 = [RNG.choice("ACGT") for _ in range(CONTIG_LEN)]
    seq2 = [RNG.choice("ACGT") for _ in range(CONTIG2_LEN)]

    def put(seq, pos, base):
        for i, b in enumerate(base):
            seq[pos + i] = b

    def splice_plus(seq, introns):     # genomic GT..AG (transcript-oriented canonical)
        for (s, e) in introns:
            put(seq, s, "GT"); put(seq, e - 2, "AG")

    def drach(seq, pos):               # A A A C A around the modified A -> valid DRACH
        put(seq, pos - 2, "AAAC"); seq[pos + 2] = "A"

    # ---- chrSyn ----
    splice_plus(seq1, (A_i1, A_i2, B_i1))
    for (s, e) in (C_i1, C_i2):        # minus-strand gene: genomic CT..AC
        put(seq1, s, "CT"); put(seq1, e - 2, "AC")
    for pos in A_M6A_SITES:
        drach(seq1, pos)
    put(seq1, Cmod - 1, "TCT")
    seq1[SNP1] = "C"; seq1[SNP2] = "A"
    put(seq1, B_TES_distal - 21, "AATAAA")
    for i in range(B_TES_prox - 50, B_TES_prox + 1):
        seq1[i] = RNG.choice("CGT")                    # AT-poor upstream -> no PAS
    put(seq1, B_TES_prox + 1, "A" * 20)                # A-rich downstream -> internal priming
    put(seq1, A_TES - 20, "AATAAA")
    put(seq1, C_e3[1] + 5, revcomp("AATAAA"))          # minus-strand PAS (transcript-upstream = higher coord)

    # ---- chrSyn2 ----
    splice_plus(seq2, [
        (OV1_e1[1], OV1_e2[0]), (OV1_e2[1], OV1_e3[0]),           # OV1 introns
        (OV2_e1[1], OV2_e2[0]),                                    # OV2 intron
        (TR_e1[1], TR_e2[0]), (TR_e2[1], TR_e3[0]), (TR_e3[1], TR_e4[0]),  # TR_L introns
        (TR_e1[1], TR_e3[0]),                                      # TR_S cassette-skip intron
        (TA_e1[1], TA_last_acc),                                   # TA intron
    ])
    drach(seq2, TR_site)
    drach(seq2, TA_site)
    for tes in (OV1_TES, OV2_TES, TR_TES, TA_TES_prox, TA_TES_dist):
        put(seq2, tes - 20, "AATAAA")                             # canonical PAS upstream of each TES

    return {CONTIG: "".join(seq1), CONTIG2: "".join(seq2)}


# ----------------------------------------------------------------------------- #
#  Read construction
# ----------------------------------------------------------------------------- #
class ReadSpec:
    __slots__ = ("name", "contig", "strand", "blocks", "tail_softclip", "pt", "mods", "snps")
    def __init__(self, name, contig, strand, blocks, tail_softclip, pt, mods, snps):
        self.name = name; self.contig = contig; self.strand = strand
        self.blocks = blocks; self.tail_softclip = tail_softclip
        self.pt = pt; self.mods = mods; self.snps = snps


def emit_read(out, refs, contig_ids, spec: ReadSpec):
    ref = refs[spec.contig]
    blocks = spec.blocks
    g2q = {}; qbases = []; qi = 0
    for (gs, ge) in blocks:
        for g in range(gs, ge):
            g2q[g] = qi; qbases.append(ref[g]); qi += 1
    for gpos, alt in spec.snps.items():
        if gpos in g2q:
            qbases[g2q[gpos]] = alt

    a = pysam.AlignedSegment()
    a.query_name = spec.name
    a.reference_id = contig_ids[spec.contig]
    a.mapping_quality = 60

    if spec.strand == "+":
        a.flag = 0; a.reference_start = blocks[0][0]
        seqstr = "".join(qbases) + ("A" * spec.tail_softclip)
        cig = []
        for i, (gs, ge) in enumerate(blocks):
            if i > 0:
                cig.append((3, gs - blocks[i - 1][1]))
            cig.append((0, ge - gs))
        cig.append((4, spec.tail_softclip))
        off = 0
    else:
        a.flag = 16; a.reference_start = blocks[0][0]
        seqstr = ("T" * spec.tail_softclip) + "".join(qbases)
        cig = [(4, spec.tail_softclip)]
        for i, (gs, ge) in enumerate(blocks):
            if i > 0:
                cig.append((3, gs - blocks[i - 1][1]))
            cig.append((0, ge - gs))
        off = spec.tail_softclip

    a.query_sequence = seqstr
    a.cigartuples = cig
    a.query_qualities = pysam.qualitystring_to_array("I" * len(seqstr))

    if spec.mods:
        ml_at = {}
        for gpos, (code, ml) in spec.mods.items():
            if gpos in g2q:
                ml_at[off + g2q[gpos]] = (code, ml)
        mm_records = []; ml_bytes = []
        for base, code in (("A", "a"), ("C", "m")):
            idxs = [i for i, b in enumerate(seqstr) if b == base]
            if not idxs:
                continue
            vals = [ml_at[i][1] if (i in ml_at and ml_at[i][0] == code) else 0 for i in idxs]
            mm_records.append(f"{base}+{code}.," + ",".join("0" for _ in idxs) + ";")
            ml_bytes.extend(vals)
        if mm_records:
            a.set_tag("MM", "".join(mm_records), "Z")
            a.set_tag("ML", array.array("B", ml_bytes))

    a.set_tag("pt", int(spec.pt), "i")
    out.write(a)


# ----------------------------------------------------------------------------- #
#  Per-sample read generators (ground-truth-encoded)
# ----------------------------------------------------------------------------- #
SAMPLES = [("M1", "mock"), ("M2", "mock"), ("M3", "mock"), ("M4", "mock"),
           ("Z1", "zikv"), ("Z2", "zikv"), ("Z3", "zikv"), ("Z4", "zikv")]


def gene_a_reads(sample, cond, rid):
    specs = []
    n_a1, n_a2 = (45, 18) if cond == "mock" else (18, 45)
    for iso, n in (("A1", n_a1), ("A2", n_a2)):
        blocks = [A_e1, A_e2, A_e3] if iso == "A1" else [A_e1, A_e3]
        for _ in range(n):
            rid[0] += 1
            name = f"{sample}_A_{iso}_{rid[0]}"
            mods = {}; snps = {}
            H = bern(0.5)
            if (H == 1) ^ (bern(0.05) == 1): snps[SNP1] = "G"
            if (H == 1) ^ (bern(0.05) == 1): snps[SNP2] = "T"
            b_pq = bern(0.5)
            mods[P] = ("a", 255 if b_pq else 0); mods[Q] = ("a", 255 if b_pq else 0)
            mods[R] = ("a", 255 if bern(0.5) else 0); mods[S] = ("a", 255 if bern(0.5) else 0)
            d_mod = bern(0.9) if iso == "A1" else bern(0.1)
            mods[D] = ("a", 255 if d_mod else 0)
            mods[COND] = ("a", 255 if (bern(0.85) if cond == "mock" else bern(0.15)) else 0)
            mods[Cmod] = ("m", 255 if bern(0.5) else 0)
            mods[Msnp] = ("a", 255 if (bern(0.9) if SNP1 not in snps else bern(0.1)) else 0)
            if iso == "A1":
                b_xy = bern(0.5)
                mods[X] = ("a", 255 if b_xy else 0); mods[Y] = ("a", 255 if (1 - b_xy) else 0)
            pt = int(RNG.gauss(130 if d_mod else 85, 8))
            specs.append(ReadSpec(name, CONTIG, "+", blocks, 18, max(1, pt), mods, snps))
    return specs


def gene_b_reads(sample, cond, rid):
    specs = []
    for iso, tes, n in (("B1", B_TES_distal, 40), ("B2", B_TES_prox, 40)):
        blocks = [B_e1, (B_e2[0], tes + 1)]
        for _ in range(n):
            rid[0] += 1
            pt = int(RNG.gauss(180 if cond == "zikv" else 120, 10)) if iso == "B1" else int(RNG.gauss(50, 8))
            specs.append(ReadSpec(f"{sample}_B_{iso}_{rid[0]}", CONTIG, "+", blocks, 18, max(1, pt), {}, {}))
    return specs


def gene_c_reads(sample, cond, rid):
    specs = []
    for _ in range(40):
        rid[0] += 1
        specs.append(ReadSpec(f"{sample}_C_C1_{rid[0]}", CONTIG, "-", [C_e1, C_e2, C_e3],
                              18, max(1, int(RNG.gauss(100, 10))), {}, {}))
    return specs


def gene_ov_reads(sample, cond, rid):
    """Two overlapping genes on chrSyn2 -> multigene filter (multi_gene_kept_by_zt)."""
    specs = []
    for _ in range(45):
        rid[0] += 1
        specs.append(ReadSpec(f"{sample}_OV1_{rid[0]}", CONTIG2, "+", [OV1_e1, OV1_e2, OV1_e3],
                              18, max(1, int(RNG.gauss(90, 8))), {}, {}))
    for _ in range(45):
        rid[0] += 1
        specs.append(ReadSpec(f"{sample}_OV2_{rid[0]}", CONTIG2, "+", [OV2_e1, OV2_e2],
                              18, max(1, int(RNG.gauss(90, 8))), {}, {}))
    return specs


def gene_tr_reads(sample, cond, rid):
    """5'-truncation gene. TR_L (e1-e2-e3-e4), TR_S (e1-e3-e4, cassette-skip), and
    truncated reads (mid-e3 -> e4) that are ASSIGNED (>=1 intron, suffix of TR_L) yet
    start 3' of the TR_L/TR_S divergence -> dropped by the truncation-aware test."""
    specs = []
    for _ in range(35):                                # TR_L (canonical, most reads)
        rid[0] += 1
        mods = {TR_site: ("a", 255 if bern(0.85) else 0)}
        specs.append(ReadSpec(f"{sample}_TRL_{rid[0]}", CONTIG2, "+",
                              [TR_e1, TR_e2, TR_e3, TR_e4], 18, max(1, int(RNG.gauss(110, 8))), mods, {}))
    for _ in range(22):                                # TR_S (cassette skip)
        rid[0] += 1
        mods = {TR_site: ("a", 255 if bern(0.15) else 0)}
        specs.append(ReadSpec(f"{sample}_TRS_{rid[0]}", CONTIG2, "+",
                              [TR_e1, TR_e3, TR_e4], 18, max(1, int(RNG.gauss(110, 8))), mods, {}))
    for _ in range(18):                                # 5'-truncated (mid-e3 -> e4), one intron
        rid[0] += 1
        mods = {TR_site: ("a", 255 if bern(0.85) else 0)}
        specs.append(ReadSpec(f"{sample}_TRt_{rid[0]}", CONTIG2, "+",
                              [(TR_trunc_start, TR_e3[1]), TR_e4], 18, max(1, int(RNG.gauss(110, 8))), mods, {}))
    return specs


def gene_ta_reads(sample, cond, rid):
    """Tandem APA with a differential mod in the shared last-exon region -> classify TANDEM_APA."""
    specs = []
    for iso, tes, p_mod, n in (("prox", TA_TES_prox, 0.2, 40), ("dist", TA_TES_dist, 0.8, 40)):
        blocks = [TA_e1, (TA_last_acc, tes + 1)]
        for _ in range(n):
            rid[0] += 1
            mods = {TA_site: ("a", 255 if bern(p_mod) else 0)}
            specs.append(ReadSpec(f"{sample}_TA_{iso}_{rid[0]}", CONTIG2, "+",
                                  blocks, 18, max(1, int(RNG.gauss(100, 10))), mods, {}))
    return specs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=str(Path(__file__).resolve().parent))
    ap.add_argument("--samtools", default="samtools")
    args = ap.parse_args()

    global pysam
    import pysam  # noqa

    outdir = Path(args.outdir)
    refdir = outdir / "reference"; refdir.mkdir(parents=True, exist_ok=True)
    bamdir = outdir / "bams"; bamdir.mkdir(parents=True, exist_ok=True)
    gtdir = outdir / "ground_truth"; gtdir.mkdir(parents=True, exist_ok=True)

    refs = build_references()
    fa = refdir / "synthetic_ref.fa"
    with open(fa, "w") as fh:
        for name in (CONTIG, CONTIG2):
            fh.write(f">{name}\n")
            s = refs[name]
            for i in range(0, len(s), 70):
                fh.write(s[i:i + 70] + "\n")
    subprocess.run([args.samtools, "faidx", str(fa)], check=True)

    gtf = refdir / "synthetic_ref.gtf"
    with open(gtf, "w") as fh:
        def tx(contig, gene_id, gene_name, tid, strand, exons):
            gstart = min(e[0] for e in exons) + 1
            gend = max(e[1] for e in exons)
            attr = f'gene_id "{gene_id}"; gene_name "{gene_name}"; transcript_id "{tid}";'
            fh.write(f"{contig}\tsyn\ttranscript\t{gstart}\t{gend}\t.\t{strand}\t.\t{attr}\n")
            for (s, e) in exons:
                fh.write(f"{contig}\tsyn\texon\t{s+1}\t{e}\t.\t{strand}\t.\t{attr}\n")
        tx(CONTIG, "GENEA", "GENE_A", "GENEA.1", "+", [A_e1, A_e2, A_e3])
        tx(CONTIG, "GENEA", "GENE_A", "GENEA.2", "+", [A_e1, A_e3])
        tx(CONTIG, "GENEB", "GENE_B", "GENEB.1", "+", [B_e1, (B_e2[0], B_TES_distal + 1)])
        tx(CONTIG, "GENEB", "GENE_B", "GENEB.2", "+", [B_e1, (B_e2[0], B_TES_prox + 1)])
        tx(CONTIG, "GENEC", "GENE_C", "GENEC.1", "-", [C_e1, C_e2, C_e3])
        tx(CONTIG2, "GENEOV1", "GENE_OV1", "GENEOV1.1", "+", [OV1_e1, OV1_e2, OV1_e3])
        tx(CONTIG2, "GENEOV2", "GENE_OV2", "GENEOV2.1", "+", [OV2_e1, OV2_e2])
        tx(CONTIG2, "GENETR", "GENE_TR", "GENETR.1", "+", [TR_e1, TR_e2, TR_e3, TR_e4])
        tx(CONTIG2, "GENETR", "GENE_TR", "GENETR.2", "+", [TR_e1, TR_e3, TR_e4])
        tx(CONTIG2, "GENETA", "GENE_TA", "GENETA.1", "+", [TA_e1, (TA_last_acc, TA_TES_prox + 1)])
        tx(CONTIG2, "GENETA", "GENE_TA", "GENETA.2", "+", [TA_e1, (TA_last_acc, TA_TES_dist + 1)])

    hdr = pysam.AlignmentHeader.from_dict(
        {"HD": {"VN": "1.6", "SO": "coordinate"},
         "SQ": [{"SN": CONTIG, "LN": CONTIG_LEN}, {"SN": CONTIG2, "LN": CONTIG2_LEN}]})
    contig_ids = {CONTIG: 0, CONTIG2: 1}

    for sample, cond in SAMPLES:
        rid = [0]
        specs = (gene_a_reads(sample, cond, rid) + gene_b_reads(sample, cond, rid)
                 + gene_c_reads(sample, cond, rid) + gene_ov_reads(sample, cond, rid)
                 + gene_tr_reads(sample, cond, rid) + gene_ta_reads(sample, cond, rid))
        sam = bamdir / f"{sample}.sam"
        with pysam.AlignmentFile(str(sam), "w", header=hdr) as out:
            for sp in specs:
                emit_read(out, refs, contig_ids, sp)
        bam = bamdir / f"{sample}.bam"
        subprocess.run([args.samtools, "sort", "-o", str(bam), str(sam)], check=True)
        subprocess.run([args.samtools, "index", str(bam)], check=True)
        sam.unlink()
        print(f"[gen] {sample} ({cond}): {len(specs)} reads -> {bam}")

    ss = outdir / "samples.tsv"
    with open(ss, "w") as fh:
        fh.write("sample\tbam\tcondition\treplicate\n")
        for i, (sample, cond) in enumerate(SAMPLES):
            fh.write(f"{sample}\t{sample}.bam\t{cond}\t{i % 4 + 1}\n")

    write_ground_truth(gtdir)
    print(f"[gen] done. reference={fa}  gtf={gtf}  bams={bamdir}  samplesheet={ss}")


def write_ground_truth(gtdir: Path):
    rows = [
        ("assemble", "GENE_A", "2 fragmentforms (A1=e1e2e3, A2=e1e3), EXACT"),
        ("assemble", "GENE_B", "2 fragmentforms (tandem APA: distal 8480, proximal 8000)"),
        ("assemble", "GENE_C", "1 fragmentform on MINUS strand"),
        ("assemble", "GENE_OV1/OV2", "two OVERLAPPING genes assembled separately"),
        ("assemble", "GENE_TR", "TR_L (e1e2e3e4) + TR_S (e1e3e4, cassette skip)"),
        ("assemble", "GENE_TA", "2 fragmentforms (tandem APA prox 10999 / dist 11499)"),
        ("splice_junctions", "all genes", "ALL_CANONICAL (incl. minus-strand GENE_C)"),
        ("apa_motifs", f"GENE_B prox {B_TES_prox}", "PAS_NONE_INTERNAL_PRIMING"),
        ("apa_motifs", f"GENE_B distal {B_TES_distal}", "PAS_CANONICAL"),
        ("multigene_filter", "GENE_OV1/OV2 reads", "multi_gene_kept_by_zt > 0 (reads overlap both genes, resolved by ZT)"),
        ("test_diffs", f"GENE_A D@{D}", "SIGNIFICANT between-isoform (A1 hi ~90% vs A2 lo ~10%)"),
        ("classify_diffs", f"GENE_A D@{D}", "SHARED_TERMINAL_EXON (A1/A2 co-terminal)"),
        ("classify_diffs", f"GENE_TA @{TA_site}", "TANDEM_APA (prox vs dist, distal favored)"),
        ("classify_diffs", "taxonomy", "all 13 categories -> test_classify_categories.py"),
        ("genotype.snp", f"SNP1 {CONTIG}:{SNP1+1} C>G / SNP2 {CONTIG}:{SNP2+1} A>T", "both discovered (~50% alt)"),
        ("genotype.snp_mod", f"SNP1 x Msnp@{Msnp}", "SIGNIFICANT (ref->mod ~90%, alt->mod ~10%)"),
        ("genotype.mechanism", f"SNP1 vs Msnp@{Msnp}", "IN_MOTIF_CORE + MOTIF_DISRUPTED, CONCORDANT"),
        ("genotype.mod_mod", f"P@{P} x Q@{Q}", "SIGNIFICANT CONCORDANT (co-dependent)"),
        ("genotype.mod_mod", f"R@{R} x S@{S}", "NOT significant, INDEPENDENT"),
        ("genotype.mod_mod", f"X@{X} x Y@{Y}", "SIGNIFICANT MUTUALLY_EXCLUSIVE"),
        ("genotype.haplotype", "SNP1+SNP2", "one haplotype block (2 SNPs in LD)"),
        ("hierarchical_stoich", "GENE_TR", "reads_dropped_as_uninformative > 0 (5'-truncated reads dropped)"),
        ("polya", "GENE_B B1 vs B2", "SIGNIFICANT differential tail (B1~120-180 vs B2~50)"),
        ("polya.tail_x_mod", f"GENE_A D@{D}", "modified reads longer tail (~130 vs ~85)"),
        ("between_conditions.mod", f"COND@{COND}", "SIGNIFICANT mock(hi~85%) vs zikv(lo~15%)"),
        ("between_conditions.isoform", "GENE_A A1 vs A2", "SIGNIFICANT usage shift (mock A1, zikv A2)"),
        ("between_conditions.tail", "GENE_B B1", "SIGNIFICANT longer tail in zikv"),
        ("calibration", "ref_df", "within-mock 2v2 null -> calibrate_between_conditions.py"),
        ("report", "html + gene_browser", "both HTML files produced"),
    ]
    with open(gtdir / "ground_truth.tsv", "w") as fh:
        fh.write("feature\tentity\texpectation\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")


if __name__ == "__main__":
    main()
