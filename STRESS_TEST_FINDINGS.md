# modulator stress-test findings (pre-v2.0.0)

Adversarial stress-testing: degenerate-input fuzzing (empty/malformed/edge inputs run
through each stage) + a five-way parallel code review of the whole `workflow/scripts/`
tree and `src/modulator/`. Synthetic regression suite stayed **33 PASS / 0 FAIL** after
every fix.

Two buckets: **FIXED** (safe, validated, committed) and **FLAGGED FOR REVIEW** (real but
subtle — a fix would change association/SNP-calling semantics I can't validate against
ground truth before release, so they're documented for the authors to decide).

---

## FIXED (committed on `test`)

### Crashes on realistic inputs
| # | file | bug | trigger |
|---|------|-----|---------|
| 1 | `test_stoichiometry_diffs.py` | empty/0-row input crashed at `df.apply(axis=1)` (multi-col→single-col) | no site passed the aggregation filter, or a gene/mod filter matched nothing → now writes empty result |
| 2 | `test_taillength_mod.py` | `_site_rows_for_chrom` returned bare `[]` but caller unpacks `r,c=…` → `ValueError` | any empty contig shard (chrM, unplaced) → now `[], []` |
| 3 | `test_taillength_diffs.py` | `kruskal(*groups)` raised on all-identical tails (≥3 fragmentforms) | degenerate low-diversity gene → now guarded, p=1.0 |
| 4 | `classify_diff_sites.py` | unguarded `json.loads(per_transcript_json)` aborted the whole stage | one malformed row → now skips the row |
| 5 | `aggregate_scrap_tx_counts.py` | bare `int("12.0")` crashed; duplicate `(code,sample)` overwrote instead of summing | non-int upstream value / duplicate rows → now tolerant + accumulates |

### Correctness
| # | file | bug | impact |
|---|------|-----|--------|
| 6 | `assemble_transcripts.py` | TES clusters are single-linkage-chained (can exceed `apa_window`), but members were selected by `|tes-rep_pos|≤apa_window` → the tail of a wide cluster was silently dropped and assigned to no isoform | **read loss / whole APA isoforms lost** on real data with TES spread; now assign by cluster span |
| 7 | `assemble_transcripts.py` | TES-boundary enforcement could write an inverted (`start>end`) exon when a terminal exon is shorter than `apa_window` | corrupt GTF exon; now clamped |
| 8 | `diffstats.py` (`continuous_diff`) | Welch t returns p=0.0 (finite!) when both groups have zero within-group variance but different means → passed the `isfinite` guard as the #1 hit | spurious top differential-tail hit on tied per-replicate medians; now p=1.0 for zero pooled variance |

### Robustness / performance
| # | file | bug | impact |
|---|------|-----|--------|
| 9 | `classify_splice_junctions.py` | GTF-vs-FASTA contig-name mismatch made every `fa.fetch` fail silently → **every** junction `NONCANONICAL`, exit 0 | now warns loudly with the offending contigs |
| 10 | `build_candidate_regions_bed.py` | a `#`-commented SNP header dropped **every** SNP interval → all SNP/haplotype outputs silently empty | now strips a leading `#` (matches the mod-BED reader) |
| 11 | `build_haplotype_blocks.py` | per-haplotype-member metadata fetch was a full boolean scan of the context frame → **O(members × reads)**, effectively a hang on a deep polymorphic locus | now a one-time `(sample,qname)` dict; **output identical** (protects the genome-wide genotype step) |

---

## FLAGGED FOR REVIEW (real, but a fix changes semantics — authors decide)

### HIGH
- **`genotype_utils.py` — asymmetric `context_key` resolution drops multi-metagene SNPs.**
  `context_key_from_row` (mod side) commits to `MG:{metagene}` whenever present;
  `context_key_from_snp_row` (SNP side) only uses `MG:` when *all* metagene tokens agree,
  else falls back to `GENE:`/`CHR:`. A SNP at an overlapping-gene locus
  (`metagene_indices="3;4"`) then never matches its mod calls (each `MG:`), and the pair is
  **silently dropped** from `snp_mod`, `snp_tx_mod_dependency`, and `haplotype_mod`. Since
  each read belongs to one metagene, this is a real false-negative class at overlapping loci
  (common in human), not conservatism. *Fix requires deciding the intended per-read grouping;
  needs an overlapping-locus ground-truth case to validate.*

### MEDIUM
- **`discover_candidate_snps.py` vs `build_molecule_snp_table.py` — inconsistent read filters.**
  Discovery counts SECONDARY/DUP/QCFAIL (`count_coverage`, `nofilter`); the molecule table
  excludes them (`pileup stepper="samtools"`). A SNP whose alt support is mostly
  duplicate/secondary reads is discovered but has ~0 molecule rows → vanishes from every
  association. Allele counts are non-reproducible between the two stages. *Fix = align the
  masks; changes which SNPs are called, so validate before shipping.*
- **`test_snp_tx_mod_dependency.py` — single-stratum pairs inflate the CMH FDR.**
  `benjamini_hochberg` is computed over all rows including `n_transcripts_tested<2` ones
  (which still got a 1-stratum CMH p-value), inflating `m` and demoting genuine multi-stratum
  hits. *Fix = BH only over ≥2-stratum rows.*
- **`discover_candidate_snps.py` — prefilter copies whole BAMs to a temp dir.**
  On N deep dRNA BAMs (~100 GB each) this adds ~N× transient disk before scanning; can
  exhaust a small working FS. (Fine on this 2.4 PB filesystem.) *Fix = stream instead of copy.*
- **`aggregate_by_transcript.py` — flat-layout sample/code detection broken (`:86`).**
  The nested-vs-flat guard compares parent-dir name to file name (always unequal) so flat
  `<sample>_<code>.bed` files collapse to one bogus sample. Only affects the legacy ZT flat
  layout (disabled by default). *Fix = compare directories, not dir-vs-file.*
- **`generate_html_report.py` / `build_gene_browser.py` — report crashes / unescaped values.**
  Several sections index optional columns guarded only on a *different* column
  (`pas_motif`, `end0`, `has_noncanonical` NaN→int, `p_adj_bh`) → `KeyError`/`ValueError`
  crash the report if an upstream table is missing a column. The gene browser interpolates
  `condition`/gene values into HTML/JS unescaped (a `</script>` in a condition name breaks the
  page). Latent on clean pipeline output; harden before exposing to arbitrary inputs.

### LOW (latent / edge)
- `discover_candidate_snps.py`: `--min-mapq 0` vs `1` flips which SAM flags are filtered (not just MAPQ).
- `build_molecule_mod_table.py` / `extract_mod_calls_pysam.py`: a `"."`-strand candidate mod site drops all its calls (strands are normally `+`/`-`).
- `test_mod_mod_assoc.py`: the `min_state_reads` gate checks only site A's marginal, so a monomorphic site B yields uninformative zero-column pairs that still enter BH.
- `genotype_utils.py` (`shard_tsv_by_chrom`): if `chrom` were the last column, the shard key keeps a trailing newline (latent; chrom isn't last today).
- `assemble_transcripts.py`: the *script* default `--min-frac 0.05` is a genome-global fraction that wipes all output at genome scale (the pipeline overrides to `0.00`; direct-CLI footgun).
- Non-determinism: `assemble_transcripts.py` collects cores via `as_completed`, so co-terminal isoforms with equal read counts can swap `T1`/`T2` labels run-to-run (contradicts the determinism claim; add a stable tiebreak on `chain_tx`).

---

## Verified NOT bugs (spot-checked, for the record)
- All six genotype 2×2/CMH contingency constructions are cell-consistent (no swapped cells).
- BH-FDR ranks only finite p-values, preserves NaN, handles ties/empty.
- The aggregation `Ndiff/Nmod` keep-filter direction/threshold matches the docs in all three implementations.
- MM/ML `(ml+0.5)/256` + canonical-wins-ties argmax matches modkit; reverse-strand `modified_bases` coordinate mapping is consistent.
- The beta-binomial LRT (null pooled-μ vs full per-group-μ, shrunk θ, `F(1,ref_df)` reference) is correctly formed.
