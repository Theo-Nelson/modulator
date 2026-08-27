# Methods §2 — modulator (draft prose, one subsection per stage)

> **Revision note (working draft).** Prose updated to match the current codebase. Three subsections
> changed materially from the previous draft and are flagged inline: **§2.7** (site-filtering reformulated
> as the per-modification NFail-SCORE k-ratio), **§2.12** (the Cochran–Mantel–Haenszel test was retired;
> SNP discovery and SNP↔modification pairing were made more conservative), and **§2.15** (between-condition
> differential modification is now resolved per transcript). Illustrative *Demo* numbers are carried over
> from the previous run and marked **[re-derive]** wherever a code change since then could have moved them;
> regenerate all reported counts from a single fresh run before submission.

**Demo used throughout:** the **14-gene benchmark panel** — ATG7, ALCAM, C8A, SERAC1, MXD1, NHSL1,
RANBP3, IP6K2, RIOK3, PCID2, VGLL4, PROZ, ASGR1, ASGR2 — from a ZIKV-infection Nanopore direct-RNA
experiment, **3 mock + 3 ZIKV replicates** (~11–17k reads/sample). These are the same genes and reads used
for the APA truthset (80 curated APA sites), so the running example and the tool-comparison benchmark (Fig 7)
refer to one identical panel.

---

## 2.0 Features and implementation

modulator runs as a single resumable command of sixteen ordered stages that reconstruct **fragmentforms**
(read-backed partial transcripts) from long-read direct-RNA alignments, quantify base modifications per
fragmentform, identify and mechanistically explain modification-stoichiometry differences between the
fragmentforms of a gene, and optionally connect sequence variants, poly(A) tails, and cis-regulatory
elements to those modifications (Fig. 1A). It wraps modkit v0.5.0 and samtools v1.22.1 and runs in Python
3.13 (pysam 0.23.3, scipy 1.16.3, pandas 2.3.3, numpy 2.3.3, matplotlib 3.10.6). Input is coordinate-sorted
genomic alignments of Nanopore direct-RNA reads carrying basecaller `MM`/`ML` modification tags and, where
available, `pt:i` poly(A) estimates. Every stage checkpoints, so a run can be stopped and resumed.

**Reference preparation.** Alignment uses a splice-aware minimap2 index of the target genome; direct-RNA
reads are converted to FASTQ with `samtools fastq -T MM,ML,ts,pt` (preserving modification and poly(A)
tags) and mapped with `minimap2 -ax splice -uf -k14 -y`, then coordinate-sorted and indexed. For samples
whose transcriptome includes a non-host contribution — e.g. Epstein–Barr virus (EBV) in immortalized
lymphoblastoid lines — the viral genome is appended to the reference as an additional contig (EBV:
NC_007605.1, added as `chrEBV`) together with a viral gene-model track, so viral transcripts are quantified
in the same fragmentform framework rather than mismapped to the host or lost. All reported metadata
(basecaller/model/version, modification code set) is taken verbatim from the input BAM `@PG`/`MM` records.

## 2.1 Fragmentform reconstruction

modulator groups reads that share an intron chain and a clustered 3′ end into one **fragmentform**. 3′ ends
within `apa_window` (default 20 nt) are merged so basecaller imprecision at the exon–poly(A) transition does
not split a fragmentform, and a shorter intron chain matching the 3′ portion of a longer one is collapsed as
a compatible truncation. Poly(A) is evaluated from the 3′ soft-clip; a fragmentform is kept only if enough
supporting reads carry a high-purity tail (`polya_support_frac`, `min_polya_purity`) and it clears
`min_reads` (default 40). Its 5′ boundary is the longest intron chain with sufficient full-length support
(`min_distal_anchor_reads`, default 2), limiting run-through across loci. Each fragmentform is assigned to a
gene by exon overlap; genes whose exons overlap on the same strand are merged into a **metagene**, and within
a metagene overlapping fragmentforms receive distinct partition indices (the `ZN` aggregation track) so their
pileups are never pooled. Non-overlapping fragmentforms may share a `ZN` track; every downstream stage
resolves a track back to the specific fragmentform that contains the base under test.
*Demo:* 126 fragmentforms across the 14 genes — 22 exact matches to the reference and 104 novel forms (17
novel-APA, 87 novel-chain).

## 2.2 Read accounting

A per-sample funnel reports total → considered → assigned reads under the same filters used for assembly,
making the denominator of every downstream rate explicit. Every read falls into exactly one PASS/FAIL bucket
and the buckets sum to the mapped total. *Demo:* e.g. sample Z1, 17,171 reads → 14,801 considered → 9,840
assigned; across the six libraries 82–86% of reads pass QC and 55–66% are assigned.

## 2.3 Splice-junction classification

Every assembled intron's donor/acceptor dinucleotides are read from the genome and classified (GT-AG
canonical / GC-AG / AT-AC minor U12 / non-canonical), summarized per gene. *Demo:* of 14 genes, 8 are
all-canonical, 4 canonical-with-GC-AG, and 2 carry a non-canonical junction.

## 2.4 Polyadenylation-signal detection

For each fragmentform 3′ end, modulator scans `upstream` (60 nt) for a polyadenylation signal (canonical
`AATAAA` or one of 11 variants) within `pas_max_distance`, scores the downstream U/GU element, and flags
**internal priming** (no PAS + A-rich genomic downstream → likely oligo-dT artifact). *Demo:* of 126
cleavage sites, **115 (91%) carry a PAS** (50 canonical `AATAAA`, 65 variant); only 11 lack one.

## 2.5 Multi-gene read resolution

Because same-strand overlapping genes are merged into one metagene (§2.1), a read's exons can overlap more
than one gene and be ambiguous as to which gene it should be counted toward. This ambiguity does **not**
reach the per-fragmentform modification counts: the modification pileup partitions by fragmentform
(`ZN`; §2.6), a read carries a `ZN` tag only if it was assigned to a reported fragmentform, and reads
without a `ZN` tag are dropped by the partitioned pileup. The ambiguity *does* reach the genotype layer
(§2.12), whose SNP scan is a tag-agnostic pileup over every read and would otherwise count an ambiguous
read toward the allele frequencies of each gene it overlaps. modulator therefore resolves multi-gene
reads to a single gene immediately before genotyping: a read is assigned by its fragmentform (`ZT`) tag
and kept when that gene is among the genes it overlaps; reads whose `ZT` gene is not among the overlaps
are treated as misassignments and dropped; and reads that overlap several genes with no `ZT` tag are
unresolvable and, by default (`multi_gene_action = scrap_unresolved`), dropped rather than double-counted.
Resolution operates purely by read inclusion/exclusion — fragmentform tags are never altered — producing a
multigene-cleaned BAM for the genotype SNP scan (and the poly(A) tables); the per-sample resolution funnel
(kept-by-`ZT` / single-gene / dropped) is reported.

## 2.6 Per-fragmentform modification calling

`modkit pileup` is run **partitioned by fragmentform** (`--partition-tag ZN`) directly on the `ZT`-tagged
BAMs, so calls are never pooled across overlapping isoforms. The partitioning also excludes reads not
assigned to a reported fragmentform (they carry no `ZN` and are placed in modkit's `ungrouped` output,
which is not aggregated) — which is why the multi-gene resolution of §2.5 is scoped to the genotype layer
and is not needed here. All modkit calling parameters are exposed (`filter_thresholds`, `mod_thresholds`
at 0.99, `max_depth`, …); modkit is pinned to 0.5.0 (0.6.x removed `--partition-tag`).

## 2.7 Aggregation and site filtering  **[revised]**

Per-sample, per-fragmentform modBEDs are merged into one table (modified / canonical / failed /
different-base counts per position per fragmentform). Two filters decide whether a (fragmentform, sample)
row is trusted. (1) A **variant/misalignment guard** drops positions dominated by a different base rather
than a modification: fail if `Ndiff > count_diff_factor × coverage`. (2) A **confident-call guard**, the
NFail-SCORE k-ratio (after Nelson et al.): fail if `Nmod < k·(Nfail + 1)`, i.e. the confident modified
calls must outweigh the low-confidence (failed) calls by a factor `k`. `k = 1` requires simply that `Nmod`
exceed `Nfail`; because the optimal threshold depends on the modification, basecaller, model and version,
`k` is **calibrated per modification** (a lookup table of empirically-derived values is provided; e.g.
Dorado SUP m6A-DRACH v2.0.0 → 0.4, SUP pseU v2.0.0 → 1.0) and may be set globally or as a per-modification
map. A site is kept if any of its rows passes both guards, and all of that site's rows are then retained.
This replaces the previous fixed `Nmod > Nfail + margin` rule, of which `k = 1` is the special case. *Demo:*
13,685 fragmentform×site×sample rows over ~427 unique modified positions **[re-derive]**, dominated by m6A
(`a`) with pseudouridine (`17802`), inosine (`17596`), 5mC and the four 2′-O-methyls — the full RNA004 code
set (each modification labelled by name in all outputs).

## 2.8 Between-fragmentform differential modification

At each position, modulator tests whether modification stoichiometry differs between the fragmentforms of a
gene (Fisher exact for 2×2, chi-square for r×2, chosen by table size; BH-FDR). *Demo:* of 364 testable
sites, ~95 differ significantly between fragmentforms (p_adj<0.05) **[re-derive]** — roughly a quarter of
quantified positions are fragmentform-specific, the phenomenon per-site callers collapse away.

## 2.9 Structural classification

Each significant site (passing `min_effect`=0.10, `fdr`=0.05) is assigned one mutually-exclusive structural
cause (Table 1), anchored to the gene's longest-3′UTR fragmentform. modulator resolves the base's status
(terminal / internal / intronic / absent) independently within the high- and low-stoichiometry fragmentforms
and the anchor — using the specific fragmentform that contains the base for each `ZN` track — then labels the
difference by a 3′-end or splicing mechanism at three levels (bucket → event → direction). *Demo:* ~72
classified sites, dominated by alternative-3′-end usage (tandem-APA + intronic polyadenylation) **[re-derive]**.

## 2.10 Novel loci

Read-backed fragmentforms matching no reference gene are rolled up into named novel loci with their sites and
junctions; intronic vs intergenic origin is recorded. *Demo:* none on this fully-annotated panel (expected).

## 2.11 Sequence-element × modification annotation

On each fragmentform's mature mRNA, modulator scans sequence cis-elements — PAS, AU-rich element (ARE), CPE,
GU-rich, RNA G-quadruplex, Kozak, uORF, 5′TOP, stop-codon context, m6Am — and reports **every overlapping
modification, unbiased across mod codes**. *Demo:* modifications concentrate where biology predicts — e.g.
~24% of uORFs modified, m6A in stop-codon-context and PAS windows, ARE carrying both m6A and pseudouridine
**[re-derive]**; 5′-anchored elements (Kozak/TOP/m6Am) show few calls, reflecting direct-RNA 5′ truncation.

## 2.12 Genotype layer *(optional)*  **[revised]**

From the same alignments, modulator calls read-backed candidate SNPs and tests their relationship to
modification and fragmentform usage. **SNP discovery** requires alt depth ≥ `min_alt_reads` (4), total
coverage ≥ 8, pooled alt fraction 0.10–0.90, base quality ≥ 20 and mapping quality ≥ 10. A site is discarded
as **multiallelic** only when its second-most-common alternative allele is *both* ≥ `min_alt_reads` *and* ≥
a coverage fraction (`multiallelic_frac`, default 0.10); gating on a fraction (rather than an absolute count)
avoids discarding clean deep heterozygous sites where a few percent third-base basecall error accumulates
past the absolute floor at high depth. The scan reports a per-filter drop tally so no positions are removed
silently.

modulator then tests three relationships, each BH-controlled: (i) **allele-specific fragmentform usage**
(allele × fragmentform table); (ii) **allele-specific modification** (allele × modified-state 2×2 on reads
covering both the SNP and the modified base), reported with a per-fragmentform breakdown so an allelic
effect can be read against the fragmentform composition of each allele; and it also tests **modification
co-dependency** (modification × modification on co-covering reads). A SNP overlapping more than one metagene
(i.e. overlapping genes) is paired against the modifications of *each* metagene it spans rather than being
collapsed to a chromosome-level context, so cis effects at overlapping loci are not lost. Co-occurring SNPs
are phased into local haplotype blocks and re-tested (allele strings only; the pooled "OTHER" bucket of
sub-threshold haplotypes is not itself tested). A mechanism classifier places each SNP relative to the
modified base (at-base / 5-mer motif core / 9-mer motif extended / proximal / distal), flags motif
disruption/creation, and flags **self-reporting** A-to-I (inosine) and pseudouridine variants — where the
apparent A>G / U>C "SNP" *is* the modification's basecall signature, not an independent genetic variant,
classified from the transcript-oriented allele change so both strands are handled; a companion table flags
every modified base coinciding with a SNP so genotype-confounded differential sites can be excluded.

*Previous draft note:* a Cochran–Mantel–Haenszel test that stratified allele-specific modification by
fragmentform has been **removed**; the per-fragmentform breakdown above (and, where conditions exist, the
per-transcript between-condition test in §2.15) serves that role more transparently.

*Demo:* ~33 read-backed SNPs with a handful of allele-specific-modification sites and modification
co-dependencies, several phased haplotype blocks, and the SNP-at-modified-base flags **[re-derive]**.

## 2.13 Poly(A) tail length

The dorado `pt:i` estimate is read per assigned read → per-fragmentform tail distributions, differential
tail between a gene's fragmentforms (Mann-Whitney for two forms, Kruskal-Wallis for more), and tail ×
modification (tail of modified vs unmodified reads, resolved **per fragmentform**). The per-fragmentform
tail×modification comparison reports how many fragmentforms reproduce the pooled direction (concordance) and
warns when the effect rests on a single fragmentform. Reads whose fragmentform is not covered by the
classification summary are back-filled from their `ZM` (metagene) tag rather than dropped. *Demo:* of the
tested genes, most show significant between-fragmentform tail differences, and a minority of tail×mod sites
show modified reads carrying systematically different tails **[re-derive]**.

## 2.14 Truncation-aware stoichiometry *(optional, off by default)*

The 5′ complement to §2.8: for fragmentforms that diverge 5′-ward, the comparison is restricted to reads
that demonstrably span their divergence point, so a read assigned to a fragmentform it never reached does not
contribute — correcting the direct-RNA 3′→5′ truncation confound.

## 2.15 Between-condition testing *(needs replicates)*  **[revised]**

Given a samplesheet `condition` column with replicates, modulator performs replicate-aware between-condition
tests, **without pooling reads across replicates** (pooling millions of reads at n=3 per group is
pseudoreplication and inflates false positives). Count-based analyses — differential modification, and
fragmentform / APA-site / splice-junction **usage** — use a beta-binomial likelihood-ratio test with
dispersion shrinkage across features, referenced to an F(1, `ref_df`) null (default `ref_df`=10); the
shrinkage weight scales with cohort size (`site_weight=auto`) so large heterogeneous cohorts retain their
true per-site dispersion. Poly(A) tail length, being continuous, is compared with Welch's t-test across
per-replicate median tails. **Differential modification is resolved per transcript partition (`ZN`) by
default:** each site is tested for the *same* fragmentform across conditions, so a between-condition change is
attributed to the specific transcript that carries it (rather than summed over all fragmentforms at the
position). Every analysis exports the per-replicate values and renders a top-by-effect figure showing the
two condition means and the individual replicate spread, so effect size and replicate agreement are visible
together. Direction is fixed by the contrast name (`{test}_vs_{reference}`; every reported Δ is test −
reference). *Demo (mock vs ZIKV, 3v3):* infection reshapes splicing and 3′-end usage strongly (a majority of
isoform/APA/junction-usage tests significant) while moving modification stoichiometry more selectively and
poly(A) tail length least; per-transcript resolution yields fewer, more-specific modification hits than the
pooled test at equal FDR **[re-derive]**.

## 2.16 Outputs

Flat tables (per-site stoichiometry, differential + classification results, genotype tables, sequence-
element table) plus one self-contained HTML report — with a sidebar table of contents, a glossary of terms
and modification codes, modification codes labelled by name throughout, and an interactive gene browser. The
per-site stoichiometry table is the substrate for between-sample/condition comparison.

## 2.17 Validation and robustness

Correctness is guarded by a synthetic-ground-truth test suite (structural-classification taxonomy, the
aggregation-engine parity between the streaming and sort back-ends, the site filters, and report-generation
edge cases) and by simulation-based calibration of the between-condition null (measured false-positive rate
vs. replicate dispersion and cohort size). Two independent adversarial audits of the codebase drove the
current design choices: reads are never pooled across replicates for between-condition inference; the
confident-call filter is calibrated per modification rather than fixed; SNPs overlapping multiple metagenes
are paired per metagene; and each filtering stage reports its drops rather than removing data silently.

---

## Suggested figures

- **Fig 1** — schematic: reads → fragmentforms (metagene / `ZN` partitioning) → the 16-stage flow,
  block-grouped, ending in the two headline outputs.
- **Fig 2** — a worked demo gene (e.g. SERAC1 or ATG7): fragmentform architecture map + per-fragmentform,
  per-sample stoichiometry heatmap with the differential site highlighted.
- **Fig 3** — classification donut/bar of the demo sites, each category with an architecture cartoon
  (pairs with Table 1).
- **Fig 4** — metagene/landmark profile: stoichiometry vs distance to stop codon / TES per mod code.
- **Fig 5** — sequence-element × modification bar (modified vs total) with the uORF / stop-context / ARE hits.
- **Fig 6** — genotype panel: an allele-specific modification example (allele × mod 2×2 + the per-fragmentform
  breakdown), and a self-report (A-to-I / pseU) flag illustration.
- **Fig 7** — APA-truthset benchmark on these same 14 genes: recovery of the 80 curated APA sites, modulator
  vs FLAIR / StringTie3 / bambu / IsoQuant.
- **Fig 8** — validation/scaling: synthetic ground-truth scorecard + dispersion-shrinkage FPR-vs-cohort-size,
  and the per-modification NFail-SCORE calibration curve.
- **Fig 9 (optional)** — cross-species application: EBV transcripts quantified in a lymphoblastoid line on
  the EBV-augmented reference (per-fragmentform modification on a named viral locus, e.g. BHRF1 / RPMS1).

---

## Table 1 — structural classification categories
*(insert your existing Table 1 here; note the current taxonomy is three-level — bucket → event → direction —
with buckets PRIVATE / SHARED_LOCAL / SHARED_DISTAL / UNEXPLAINABLE and events including ALT_DONOR,
ALT_ACCEPTOR, ALT_EXON, ALT_POLYA_SITE, RETAINED_INTRON, IPA_EXTENSION, IPA_NO_EXTENSION under SHARED_LOCAL,
and DISTAL_APA / DISTAL_SPLICING under SHARED_DISTAL.)*
</content>
