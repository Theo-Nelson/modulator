# Ground truth & validation results

Coordinates are 0-based genomic. Mod codes follow dorado/modkit: `a` = m6A (on A),
`m` = 5mC (on C). A base is encoded modified with `ML`=255 (p≈0.996), canonical with
`ML`=0 — unambiguous for both the modkit pileup path (0.99 threshold) and the
genotype per-read argmax path.

## GENE_A modification-site map (`chrSyn`, `+` strand)

| site | pos | mod | exon | designed behaviour |
|------|-----|-----|------|--------------------|
| **P** | 1080 | m6A | e1 | **co-dependent with Q** |
| **Q** | 1160 | m6A | e1 | **co-dependent with P** |
| **R** | 1300 | m6A | e1 | **independent** |
| **S** | 1380 | m6A | e1 | **independent** |
| **X** | 2100 | m6A | e2 | **mutually exclusive with Y** (e2 = A1 only) |
| **Y** | 2260 | m6A | e2 | **mutually exclusive with X** |
| **D** | 3050 | m6A | e3 | **differential between isoforms** (A1 ≈ 90 % vs A2 ≈ 10 %) |
| **Msnp** | 3200 | m6A | e3 | **controlled by SNP1** (ref ≈ 90 % vs alt ≈ 10 %); DRACH |
| **Cmod** | 3350 | 5mC | e3 | independent; second mod code |
| **COND** | 3450 | m6A | e3 | **differential between conditions** (mock ≈ 85 % vs zikv ≈ 15 %) |

SNPs (≈50 % alt, in **linkage** → one haplotype block): **SNP1** `chrSyn:3202 C>G`
(+1 of Msnp, disrupts DRACH, controls Msnp); **SNP2** `chrSyn:1251 A>T`.

Isoform usage is condition-linked (mock favours A1, zikv favours A2). Poly(A):
GENE_B B1 ≈ 120 nt (mock) / 180 (zikv), B2 ≈ 50 nt; GENE_A tail coupled to site-D
modification (≈130 vs ≈85) for a tail×mod signal.

## `chrSyn2` scenario genes

- **GENE_OV1 / GENE_OV2** — same strand, OV1's last exon overlaps OV2's first exon, so
  every full read of each gene exonically overlaps *both* gene unions → multigene filter
  routes them (`multi_gene_kept_by_zt`).
- **GENE_TR** — `TR_L` (e1e2e3e4) vs `TR_S` (e1e3e4, cassette skip), diverging ~1800 nt
  from the 3′ end. **5′-truncated reads** start mid-e3 (3′ of the divergence), are
  assigned to `TR_L` by suffix collapse, and carry a mod call in the shared terminal exon
  → the truncation-aware test drops them as uninformative.
- **GENE_TA** — tandem APA (prox TES 10999 / dist TES 11499) with a differential mod in
  the shared last-exon region → `classify_diffs` → `TANDEM_APA`.

## Feature-by-feature expected output vs. observed

`validate_outputs.py` = **30 PASS / 0 FAIL**; `test_classify_categories.py` = **13/13**.

| feature | expectation | observed |
|---------|-------------|----------|
| **assemble** | 7 genes; GENE_A = **2** fragmentforms; all EXACT | matches |
| **read funnel** | all reads mapped/considered/assigned | matches |
| **splice junctions** | all `ALL_CANONICAL` incl. minus-strand GENE_C | matches |
| **apa_motifs** | GENE_B distal→`PAS_CANONICAL`, proximal→`PAS_NONE_INTERNAL_PRIMING` | matches |
| **multigene filter** | OV1/OV2 reads → `multi_gene_kept_by_zt` > 0 | 720 total |
| **test_diffs** | site D top hit; P/Q/R/S not significant | D p_adj 3.8e-80 |
| **classify_diffs** | GENE_TA→`TANDEM_APA`; GENE_A→`SHARED_TERMINAL_EXON`; **all 13 categories** (unit) | matches |
| **SNP discovery** | both SNPs, correct ref/alt/≈50 % | both found |
| **snp_mod** | SNP1×Msnp strong (effect 0.83); no independent site strong (effect<0.25) | 0.83 vs max-indep 0.15 |
| **snp_mod mechanism** | SNP1 → `IN_MOTIF_CORE`/`MOTIF_DISRUPTED`, CONCORDANT | matches |
| **haplotype** | 1 block, 2 SNPs, LD haplotypes A\|C & T\|G (+recombinants) | matches |
| **mod_mod (co-dependency)** | **P×Q CONCORDANT sig; R×S INDEPENDENT n.s.; X×Y MUT.EXCL. sig** | OR 2.5e5 / 1.03 / 0.0 |
| **hierarchical_stoich** | GENE_TR `reads_dropped_as_uninformative` > 0 | 144 dropped |
| **polya** | GENE_B B1 vs B2 tail sig; tail×mod at D sig | matches |
| **between-cond mod** | COND sig (mock hi vs zikv lo) | delta −0.74, p_adj 7e-8 |
| **between-cond isoform / tail** | GENE_A A1↔A2 usage sig; GENE_B tail sig | matches |
| **calibration** | within-mock 2v2 null → ref_df grid report | runs (see README) |
| **report + gene browser** | both HTML files | produced |

**Headline (mod_mod):** the co-dependent pair is `CONCORDANT` at p_adj ≈ 1e-149, the
independent pair is `INDEPENDENT` and non-significant (p_adj ≈ 1.0), the
mutually-exclusive pair is `MUTUALLY_EXCLUSIVE` at p_adj ≈ 1e-74.

### Emergent / expected correlations (correct, not bugs)

- **D × COND** co-localize significantly: the design makes `condition → isoform → D`
  and `condition → COND`, so they genuinely co-vary. Co-occurrence ≠ direct mechanism.
- **SNP2 × Msnp** is significant via **LD** (SNP2 is linked to SNP1, which controls Msnp).
  The `ground_truth.tsv` label "SNP2 controls no modification" is about *causation*, not the
  (real) LD association.
- **A control m6A/5mC site can cross p<0.05 by chance.** With ~250 reads per allele a
  null site can show a ~3-SD (≈0.13) fluctuation and a well-powered test flags it — exactly
  why the pipeline says *rank by |delta|, not FDR alone*. The validator therefore gates the
  controls on **effect size** (all < 0.25) versus the designed SNP1×Msnp effect (0.83), not
  on FDR.

## Fixed bug (was "Bug 1")

**Fragmentform duplication when a TES lands on a regionization core boundary** — FIXED
(`assemble_transcripts.py`, commit *"collapse duplicate fragmentforms straddling a core
boundary"*). Cores tile the contig with shared boundaries and fetch reads with ±pad; a
fragmentform whose TES sat exactly on a boundary was emitted by both adjacent cores,
duplicating it (GENE_A → 4 forms) and double-counting reads in the metrics / tx_counts /
usage tables. The fix collapses identical `(chrom, strand, tes, chain_tx)` fragmentforms
(keeping the copy with the most members). GENE_A now assembles as 2 forms; usage fractions
are correct. `validate_outputs.py` guards against regression.

## Observation (not a bug)

With dorado-style `A+a.` MM tags, modkit reports **every** covered A/C as a site (implicit
canonical), so the RAW aggregate scales with base count; the FILTERED table discards the
0 %-modified noise, so nothing downstream is affected.
