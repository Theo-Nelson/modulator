# Changelog

All notable changes to Modulator are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the project uses
[Semantic Versioning](https://semver.org/).

## [2.0.0] - 2026-10-14

Archived at Zenodo: https://doi.org/10.5281/zenodo.23381616

This is the release described in the Modulator manuscript (Nelson et al. 2026)
and used for every analysis in it (HG002, Huh7 Mock and ZIKV-infected). It is a
rewrite of the 1.x Snakemake workflow as a Python package with a single
`modulator` command. The 1.x history is preserved on the `legacy-v1` branch and
at tag `v1.0.0`; 2.0.0 shares no code with it.

### Added

- `modulator` command line interface with `run`, `validate-config`, and `demo`
  subcommands; `--stages` to run a subset of the sixteen stages and `--resume`
  to restart from existing stage outputs.
- Fragmentform assembly directly from spliced genomic alignments, without a
  reference transcriptome. Every read, 5'-truncated or not, is assigned to a
  fragmentform from its own 3' end and junction chain (`ZN` tag). Fragmentforms
  are classified against the annotation as exact match, new junction chain,
  new 3' end, or new locus.
- Per-fragmentform modification calls with `modkit pileup` (modkit 0.5.0),
  with the `Nfail + 1 <= Nmod` and `Ndiff <= 3 x Nvalid_cov` site filters and
  per-fragmentform read caps.
- Between-fragmentform stoichiometry test: library-stratified Cochran Mantel
  Haenszel test with a Monte Carlo exact fallback (999 / 9,999 / 99,990 draws)
  whenever an expected count is below 5, Benjamini Hochberg correction across
  all sites, and a coverage-weighted effect size.
- Classification of every between-fragmentform difference into twelve classes
  (four private, six shared local, two shared distal) or "unexplainable", plus
  a table of private sites and the distance from each site to the separating
  junction or poly(A) site along the spliced RNA.
- Genotype stage: SNP discovery from the reads (minimum alternate fraction
  0.30, maximum 0.90, base quality 20, mapping quality 10), haplotype blocks
  of 2 to 4 SNPs, SNP-to-modification association (library-stratified CMH,
  Fisher when one library is informative) with a distance-class annotation
  (replaced base, motif core, extended motif, proximal, distal), and
  SNP-to-fragmentform usage association. Per-read molecule tables are written
  as Parquet.
- Between-conditions stage: replicate-aware beta-binomial likelihood ratio
  test with dispersion shrinkage for modification stoichiometry, and the same
  engine for fragmentform, poly(A) site, and junction usage. Sites where the
  modification itself is read as a SNP are flagged.
- Poly(A) tail stage: per-fragmentform tail length distributions, differential
  tail length, and tail by modification.
- Supporting stages: read statistics, splice junction and poly(A) motif
  annotation, novel loci, sequence elements, multigene read filter, and a
  hierarchical stoichiometry summary.
- Self-contained HTML report and gene browser, streamed in chunks so that
  30-sample runs render.
- Bundled demo and test BAMs (`resources/test_bams/`, eight datasets over 14
  genes from the six Huh7 libraries) re-cut from the manuscript's uniform
  processing (dorado 2.1.0 sup v6.0.0, all-context modification models v1,
  minimap2 2.29, samtools 1.22.1), with their provenance in
  `resources/test_bams/README.md`; the published sample report in
  `docs/sample_report` is regenerated from them on this version. A synthetic
  three-exon fixture is used by the regression tests.
- `ADVANCED_USAGE.md` with the full parameter, output, and HPC reference.

### Changed

- The Snakemake workflow of 1.x is replaced by the `modulator` package; stage
  scripts live in `workflow/scripts/` and are launched by the runner.
- Reads are grouped by 3' end and junction chain (fragmentforms) instead of by
  annotated transcript, so a reference transcriptome is no longer required.
- Default SNP alternate-allele floor raised from 0.10 to 0.30 after a Genome in
  a Bottle HG002 check (precision 0.82 to 0.93 at recall 0.93).
- Intronic poly(A) classification uses the donor position rather than exon
  length, which fixes an alternative-donor versus intronic-poly(A) inversion.
- Gene-level condition summaries use the mean of fragmentform medians rather
  than pooled reads, so they do not change with fragmentform usage.

### Fixed

- Stale `modkit` partitions from a previous run are no longer aggregated.
- The modification table no longer fabricates unmodified observations for
  modifications a read never assessed.
- Contingency tests reduce empty rows and columns before testing, the exact
  single-stratum r x 2 test is used at all six call sites, and the Monte Carlo
  p value never falls back to the chi-square approximation.
- Usage tests no longer truncate dotted gene names or merge genes.
- Memory and speed work for 30-sample runs: parallel per-site fits, streamed
  reads, per-window checkpoints, and bounded per-read tables.

## [1.0.0] - 2026-03-21

Snakemake workflow release. See the `legacy-v1` branch.
