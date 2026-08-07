![modulator](workflow/images/modulator_banner.png)

# modulator

> Transcript-fragment-resolved RNA-modification analysis for long-read direct-RNA sequencing.

**modulator** takes aligned long-read (e.g. ONT direct-RNA) BAMs that carry
base-modification tags and, in a single command:

1. **assembles** fragment transcripts and partitions reads per transcript (`ZN`/`ZT` tags) — here *fragment* means the fragmented / partial transcript a long read represents, not an assumed full-length molecule;
2. **calls modifications per transcript** (via [`modkit`](https://github.com/nanoporetech/modkit) pileup, partitioned by transcript);
3. finds sites that are **differentially modified _between_ transcripts** of the same gene, and **classifies _why_** (alternative polyadenylation, intronic polyadenylation, EJC, splicing, …), anchored to each gene's longest-3′UTR transcript;
4. connects to **genotype** — read-backed SNPs, SNP-modification and SNP-fragment associations, fragment-conditioned dependency tests, and local haplotype blocks;
5. writes a single self-contained **HTML report**.

The supported interface is the `modulator` Python CLI.

**Full parameter, output, and HPC reference -> [ADVANCED_USAGE.md](ADVANCED_USAGE.md)**

---

## Contents

- [Installation](#installation)
- [Input requirements](#input-requirements)
- [Quick start](#quick-start)
- [What modulator does](#what-modulator-does)
- [Key outputs](#key-outputs)
- [Command-line interface](#command-line-interface)
- [Configuration](#configuration)
- [Running on a cluster](#running-on-a-cluster)
- [Citation & contributors](#citation--contributors)

## Installation

```bash
# 1. micromamba (https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html)
# 2. clone
git clone https://github.com/Theo-Nelson/modulator.git
cd modulator
# 3. create + activate the environment
micromamba env create -n modulator -f workflow/envs/modulator.yaml
micromamba activate modulator
# 4. install the CLI 
python -m pip install -e .
```

This provides the `modulator` command (Python ≥ 3.11; bundles `modkit` 0.5.0,
`samtools`, `pysam`, `pandas`, `scipy`, `matplotlib`). No-install fallback:
`PYTHONPATH=src python -m modulator …`.

Quick code-path check:

```bash
python workflow/scripts/regression_smoke_checks.py
python workflow/scripts/genotype_regression_smoke_checks.py
```

## Input requirements

- **BAMs**: coordinate-sorted and indexed (`.bai`), aligned to your genome
  (e.g. `minimap2 -ax splice -uf` for ONT dRNA), one per sample.
- **Modification tags**: `MM`/`ML` tags from a modification-aware basecaller
  (e.g. Dorado). modulator does not call modifications from the signal — it
  summarizes the modification tags already present in the reads.
- **Reference**: matching genome FASTA (+ `.fai`) and a GTF annotation.

## Quick start

```bash
# bundled demo (fast, single gene)
modulator demo --reference-fa ref.fa --reference-gtf ref.gtf

# a real run over a directory of BAMs
modulator run \
  --config config/config.yaml \
  --jobs 8 \
  --set \
    reference_fa=/path/to/ref.fa \
    reference_gtf=/path/to/ref.gtf \
    bams_dir=/path/to/bams \
    bam_glob='*.bam' \
    prefix=my_run
```

Outputs land in `results/` under the working directory.

## What modulator does

Fifteen stages, run in order and individually **resumable** (`--resume`):

| # | stage | what it does | headline output |
|---|-------|--------------|-----------------|
| 1 | `assemble` | build transcript models from intron chains; tag reads with metagene-aware `ZN` partitions | `<prefix>.gtf`, `<prefix>_partition_map.tsv` |
| 2 | `read_stats` | per-sample read-retention funnel + length summaries | `<prefix>_per_sample_read_stats.tsv` |
| 3 | `splice_junctions` | donor/acceptor dinucleotides of every assembled intron | `<prefix>_splice_junctions.tsv` |
| 4 | `apa_motifs` | polyadenylation-signal check per APA site + internal-priming flag | `<prefix>_apa_motifs.tsv` |
| 5 | `multigene_filter` | resolve/keep reads over overlapping genes -> cleaned tagged BAMs | `zt_filtered/*.bam` |
| 6 | `modkit_zn` | `modkit pileup` partitioned by `ZN` (per-transcript modification calls) | `modkit_zn/<sample>/*.bed` |
| 7 | `aggregate_zn` | merge into per-site × transcript × sample stoichiometry | `<prefix>_FILTERED_sites_long.tsv` |
| 8 | `novel_loci` | read-backed loci matching no reference gene | `<prefix>_novel_loci.tsv` |
| 9 | `test_diffs` | between-transcript differential-modification test per site | `<prefix>__ZN_site_diff_results.tsv` |
| 10 | `classify_diffs` | assign a structural category to each significant site (+ figures) | `<prefix>__ZN_site_classified.tsv` |
| 11 | `genotype` *(optional)* | SNP discovery, SNP-mod/fragment association, dependency, haplotypes, **and why each SNP changes a modification** | `genotype/<prefix>_*.tsv` |
| 12 | `polya` | dorado poly(A) tail length per read -> per-fragmentform distributions, differential tail, tail × modification | `polya/<prefix>_*.tsv` |
| 13 | `hierarchical_stoich` *(optional)* | truncation-aware differential stoichiometry: compares fragmentforms using only reads that demonstrably span their divergence point (the 5' complement to `test_diffs`) | `<prefix>_hierarchical_stoich.tsv` |
| 14 | `between_conditions` *(needs a samplesheet)* | replicate-aware differential modification / isoform / APA / junction usage / tail length between conditions | `between_conditions/<prefix>_<contrast>_*.tsv` |
| 15 | `report` | self-contained HTML report **+ interactive gene browser** | `report/<prefix>_report.html`, `report/<prefix>_gene_browser.html` |

### Comparing conditions (samplesheet)

By default samples are discovered by globbing `bams_dir`, and every "differential"
stage compares *transcripts within a pooled population* — there is no notion of a
condition. To compare **conditions** (e.g. infected vs mock), point `samplesheet:`
at a TSV:

```tsv
sample	bam	condition	replicate
M1	HornerLab_M1pA_*.bam	mock	1
Z1	HornerLab_Z1pA_*.bam	zikv	1
```

The samplesheet then becomes the **sample source** (it replaces `bam_glob`): each
BAM is symlinked to `<sample>.bam`, so your short ids — not 60-character BAM stems —
are the sample names in every output. `condition` unlocks the `between_conditions`
stage. See [ADVANCED_USAGE.md -> Samplesheet](ADVANCED_USAGE.md#samplesheet--between-condition-comparisons).

## Key outputs

Where to look first (all under `results/`):

- **`report/<prefix>_report.html`** — start here. Self-contained, with
  per-modification structural-category graphs, per-sample stoichiometry, and the
  genotype/SNP tables. Images are written to a sidecar `<prefix>_report_files/`
  folder so the report stays small and opens in any browser.
- **`aggregate_zn/<prefix>_FILTERED_sites_long.tsv`** — per-site, per-transcript,
  per-sample modification stoichiometry (the quantitative core).
- **`test_diffs/<prefix>__ZN_site_diff_results.tsv`** — sites differentially
  modified between transcripts (effect size + BH-FDR).
- **`test_diffs/<prefix>__ZN_site_classified.tsv`** (+ `__figs_by_category_arch/`)
  — the structural reason each site differs, with transcript architecture maps.
- **`genotype/<prefix>_*.tsv`** — SNPs, SNP-mod / SNP-fragment / dependency,
  and haplotype blocks (when `genotype.enable=true`). Includes
  `<prefix>_snp_mod_mechanism.tsv`: *why* each SNP changes a modification (at the
  modified base / in the DRACH 5-mer / 9-mer / proximal / distal), whether the alt
  allele breaks the m6A consensus, and whether the data agree with that prediction.
  It also flags **self-reporting** variants — A-to-I and pseudouridine change the
  basecall, so they get called as SNPs at their own site and the association is
  circular (~12% of significant hits in test data).
- **`polya/<prefix>_*.tsv`** — poly(A) tail length per fragmentform, differential
  tail between a gene's isoforms, and tail × modification.
- **`between_conditions/<prefix>_<contrast>_*.tsv`** — replicate-aware condition
  comparisons (needs a `samplesheet` with a `condition` column).
- Every figure is written as **both PNG and SVG** (SVG for publication).

See [ADVANCED_USAGE.md](ADVANCED_USAGE.md#outputs) for the complete file list.

## Command-line interface

```bash
modulator run             --config config/config.yaml [--jobs N] [--stages …] [--resume] [--set k=v …]
modulator demo            --reference-fa ref.fa --reference-gtf ref.gtf [--dataset …] [--mode full]
modulator validate-config --config config/config.yaml [--set k=v …]
```

- `--set` takes simple `key=value` or `nested.key=value` overrides.
- `--stages` runs a comma-separated subset (e.g. `--stages assemble,read_stats`).
- `--resume` skips stages whose outputs already exist.

Full flag reference: [ADVANCED_USAGE.md -> CLI](ADVANCED_USAGE.md#command-line-interface).

## Configuration

Nested settings live in `config/config.yaml` under sections `assembler`,
`multigene`, `modkit`, `aggregation`, `test_diffs`, `classify_diffs`,
`genotype`, and `report`; override any of them on the command line with
`--set nested.key=value`. Every knob is documented in
[ADVANCED_USAGE.md -> Stages & parameters](ADVANCED_USAGE.md#stages--parameters).

## Running on a cluster

modulator is designed to run one isolated project directory per sample-set and
to checkpoint between stages. A minimal Slurm pattern (full recipe, resume/clean
re-runs, and runtime collection in
[ADVANCED_USAGE.md -> HPC](ADVANCED_USAGE.md#running-on-an-hpc-cluster-slurm)):

```bash
modulator run --workdir "$RUNDIR" --config config/config.yaml --jobs "$JOBS" \
  --set reference_fa=ref.fa reference_gtf=ref.gtf bams_dir="$BAMS" bam_glob='*.bam' prefix="$NAME"
```

## Citation & contributors

- [Theodore Nelson](https://github.com/Theo-Nelson), Weill Cornell Medicine
- [Michael Goneos](https://github.com/mgoneos), Weill Cornell Medicine

Supported by NSF ACCESS Allocation BIO240371. T.M.N. was supported by an MSTP
grant (NIGMS/NIH T32GM152349) to the Weill Cornell/Rockefeller/Sloan Kettering
Tri-Institutional MD-PhD Program.
