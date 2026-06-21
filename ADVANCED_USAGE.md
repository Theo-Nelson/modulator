# modulator — Advanced Usage

Complete reference for the command-line interface, configuration, every stage's
parameters, outputs, the genotype module, HPC execution, and troubleshooting.
For installation and a quick start, see [README.md](README.md).

## Contents

- [Command-line interface](#command-line-interface)
- [Configuration file](#configuration-file)
- [Stages & parameters](#stages--parameters)
  - [1. Assembly](#1-assembly)
  - [2. Multigene filter](#2-multigene-filter)
  - [3. modkit pileup (modkit_zn)](#3-modkit-pileup-modkit_zn)
  - [4. Aggregation (aggregate_zn)](#4-aggregation-aggregate_zn)
  - [5. Differential test (test_diffs)](#5-differential-test-test_diffs)
  - [6. Classification (classify_diffs)](#6-classification-classify_diffs)
  - [7. Genotype (optional)](#7-genotype-optional)
  - [8. Report](#8-report)
- [Outputs](#outputs)
- [Running on an HPC cluster (Slurm)](#running-on-an-hpc-cluster-slurm)
- [Performance & troubleshooting](#performance--troubleshooting)
- [Recent changes](#recent-changes)

---

## Command-line interface

```bash
modulator run --config config/config.yaml [--workdir .] [--jobs N] [--stages all] [--resume] [--set k=v …]
```

| flag | default | meaning |
|------|---------|---------|
| `--config` | `config/config.yaml` | YAML config (relative to `--workdir` or absolute). |
| `--workdir` | `.` | Project dir containing `workflow/`, `config/`, `results/`, `resources/`. |
| `--jobs` | `1` | Independent sample-level jobs in parallel. (`genotype.jobs` is separate.) |
| `--stages` | `all` | Comma-separated subset of: `assemble, read_stats, multigene_filter, modkit_zn, aggregate_zn, test_diffs, classify_diffs, genotype, report`. |
| `--resume` | off | Skip stages whose outputs (or `results/.checkpoints/<stage>.done` markers) already exist. |
| `--set` | — | Simple overrides: `key=value` / `nested.key=value` (e.g. `genotype.enable=true`). |

```bash
modulator demo --reference-fa ref.fa --reference-gtf ref.gtf [--dataset NAME] [--mode full] [--prefix …] [--jobs 2] [--set …]
modulator validate-config --config config/config.yaml [--set …]
```

- `demo` runs a fast bundled dataset (default: MXD1-only, genotype off). Datasets
  include `RPL13_reads` (single gene) and `ALCAM_NHSL1_SERAC1_MXD1_RIOK3_reads`
  (five genes); `--mode full` enables the heavier path.
- `validate-config` loads the config, applies `--set`, and prints the resolved
  project root — use it to dry-check overrides before a long run.

The package runner executes the existing `workflow/scripts/*` directly, so the
algorithms are unchanged; only the launch UX differs.

## Configuration file

`config/config.yaml` holds top-level inputs and the nested per-stage sections.
Top-level keys (usually set with `--set`): `reference_fa`, `reference_gtf`,
`bams_dir`, `bam_glob`, `prefix`, `threads`. Nested sections: `assembler`,
`multigene`, `modkit` (`common` / `zn` / `zt`), `aggregation` (`zn` / `zt`),
`test_diffs`, `classify_diffs`, `genotype`, `report`. Override anything with
`--set nested.key=value`.

---

## Stages & parameters

### 1. Assembly

Builds isoform models from read intron chains and tags reads with a
**metagene-aware `ZN` partition index** (non-overlapping transcripts within one
overlapping metagene can share a partition without discarding overlapping-locus
reads).

| Parameter | Type | Default | Range / Options | Description |
|-----------|------|---------|-----------------|-------------|
| `primary_only` | bool | `true` | `true`/`false` | Use only primary alignments. |
| `min_mapq` | int | `10` | `0–60` | Minimum MAPQ to keep a read. |
| `min_introns_read` | int | `1` | `0–10` | Require ≥ this many introns per read (`0` keeps unspliced reads — e.g. mitochondrial). |
| `require_softclip3p` | int (nt) | `0` | `0–40` | Minimum 3′ soft-clip length (poly(A/T) proxy). |
| `apa_window` | int (nt) | `20` | `10–100` | Window to cluster 3′ ends into APA groups. |
| `tes_window` | int/null | `null` | — | Overrides `apa_window` for TES clustering if set. |
| `min_reads` | int | `40` | `1–200` | Minimum total reads to retain an isoform. |
| `min_frac` | float | `0.00` | `0.0–0.2` | Minimum fraction of all reads supporting an isoform. |
| `min_introns` | int | `1` | `0–10` | Minimum introns per assembled transcript (`0` keeps single-exon/intronless models). |
| `min_polya_length` | int (nt) | `12` | `8–30` | Minimum soft-clip length for poly(A/T) evidence. |
| `min_polya_purity` | float | `0.5` | `0.4–0.95` | A/T fraction in the soft-clipped tail to count as poly(A/T). |
| `polya_support_frac` | float | `0.5` | `0.3–0.9` | Minimum fraction of an isoform's reads with poly(A/T) support. |
| `tes_match_tol` | int (nt) | `25` | `5–50` | Tolerance for matching transcript end sites to reference TES. |
| `exact_tes_tol` | int (nt) | `10` | `2–20` | Distance to label a match `EXACT` vs `NOVEL_APA`. |
| `min_distal_anchor_reads` | int | `2` | `1–20` | Exact-chain reads needed before a longer chain absorbs shorter suffix-compatible chains. |
| `min_distal_anchor_frac` | float | `0.05` | `0.0–1.0` | Exact-chain fraction within a suffix family required for absorption. |
| `min_exact_canonical_reads` | int | `1` | `1–20` | Minimum exact-chain reads for a canonical to participate in suffix collapse. |
| `write_zt_bams` | bool | `false` | — | Write one BAM per transcript per sample (many files). |
| `write_zt_tagged_sample_bams` | bool | `true` | — | Write one ZT/ZN-tagged BAM per sample (enables downstream modification calling). |
| `emit_modkit_manifest` | bool | `false` | — | Also write `zt_bams/modkit_manifest.tsv`. |
| `min_reads_per_sample_for_mod` | int | `5` | `1–50` | Minimum per-sample read support to make a per-transcript BAM. |
| `min_total_reads_for_mod` | int | `20` | `10–200` | Minimum total read support for modkit BAM output. |
| `status_every` | int | `0` | `0+` | Print assembly progress every N reads (`0` = off). |

Notes: a longer suffix-compatible chain absorbs shorter chains only with direct
exact-chain support from its own distal unique 5′ structure; a read provides
poly(A/T) support if 3′ soft-clip ≥ `min_polya_length` **and** purity ≥
`min_polya_purity`; isoforms must pass **all** filters; `ZN` is a metagene-aware
partition index, not a within-gene transcript index.

> **Intronless transcriptomes (e.g. chrM/mitochondria):** set
> `assembler.min_introns=0` **and** `assembler.min_introns_read=0`, or assembly
> will discard every (intronless) read.

### 2. Multigene filter

Scans each tagged BAM, summarizes overlap against assembled gene exon-unions,
and writes a cleaned BAM for `modkit`. Overlapping reads are kept and categorized.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `enable` | bool | `true` | Run the overlap filter between assembly and modkit. |
| `zero_gene_action` | string | `"keep"` | Reads overlapping zero assembled genes: `keep` or `scrap`. |

### 3. modkit pileup (modkit_zn)

`modkit pileup` is run **per `ZN` partition** (`--partition-tag ZN`) on the
cleaned tagged BAMs. Keys under `modkit.common` map directly to modkit flags
(see ONT's [Advanced Usage](https://github.com/nanoporetech/modkit/blob/master/book/src/advanced_usage.md)):

| Key | Type | Default | Maps to |
|-----|------|--------:|---------|
| `log_file_template` | str | `results/{which}/{sample}.log` | `--log-filepath` |
| `region` | str/null | `null` | `--region` |
| `max_depth` | int | `1000` | `--max-depth` |
| `include_bed` | str/null | `null` | `--include-bed` |
| `include_unmapped` | bool | `false` | `--include-unmapped` |
| `edge_filter` / `invert_edge_filter` | str/bool | `null`/`false` | `--edge-filter` / `--invert-edge-filter` |
| `threads` | int | top-level `threads` | `-t/--threads` |
| `interval_size` | int | `100000` | `--interval-size` |
| `queue_size` | int | `1000` | `--queue-size` |
| `chunk_size` | int/null | `null` | `--chunk-size` |
| `num_reads` | int | `10042` | `--num-reads` |
| `sampling_frac` / `seed` | float/int | `null`/`null` | `--sampling-frac` / `--seed` |
| `sample_region` / `sampling_interval_size` | str/int | `null`/`1000000` | `--sample-region` / `--sampling-interval-size` |
| `no_filtering` | bool | `false` | `--no-filtering` |
| `filter_thresholds` | list[str] | `["A:0.8","C:0.8","G:0.8","T:0.8"]` | `--filter-threshold` (repeatable) |
| `mod_thresholds` | list[str] | eight codes at `0.99` | `--mod-threshold` (repeatable) |
| `ignore` | list[str] | `[]` | `--ignore` |
| `motif` / `cpg` | list/bool | `[]`/`false` | `--motif <m> <off>` / `--cpg` |
| `combine_mods` / `combine_strands` | bool | `false` | `--combine-mods` / `--combine-strands` |
| `only_tabs` / `mixed_delim` | bool | `false` | output delimiter control |
| `suppress_progress` | bool | `true` | `--suppress-progress` |

Partition tags are hard-coded to split reads by transcript assignment:

| Section | Key | Default | Meaning |
|---------|-----|--------:|---------|
| `zn` | `partition_tag` | `"ZN"` | metagene-aware partition index from the assembler |
| `zn` | `per_mod_bed` | `true` | one BED per modification |
| `zt` | `partition_tag` | `"ZT"` | gene+transcript human-readable code |

> `modkit` is pinned to **0.5.0**: 0.6.x removed `--partition-tag`, which the
> per-isoform pileup relies on.

### 4. Aggregation (aggregate_zn)

Merges the per-`ZN` modkit bedMethyl into long tables and per-gene/per-mod pivots.

> **Site-keeping rule (when `filter_enable`):** a row **fails** if
> `Ndiff > count_diff_factor·Nvalid_cov` **or** `Nmod ≤ Nfail + mod_fail_margin`.
> A **site is kept** in *FILTERED* outputs if **any** row at that site passes;
> when kept, **all** rows for that site (every ZN/sample) are retained. `min_cov`
> only zeroes the displayed `frac_modified` (`Nvalid_cov < min_cov`); it does not
> affect pass/fail.

| Parameter | Type | Default | Scope | Description |
|-----------|------|--------:|-------|-------------|
| `engine` | str | `"stream"` | ZN | `stream` (per-chrom k-way merge of pre-sorted, tabix-indexed beds; parallel + resumable) or `sort` (legacy). |
| `jobs` | int | `12` | ZN | Per-chromosome parallelism for the streaming engine. |
| `filter_enable` | bool | `true` | ZN/ZT | Enable the site-keeping rule above. |
| `count_diff_factor` | float | `3.0` | ZN/ZT | Factor for the `Ndiff` fail term. |
| `mod_fail_margin` | int | `1` | ZN/ZT | Extra margin on `Nfail` in the `Nmod` fail rule. |
| `emit_raw` | bool | `false`* | ZN/ZT | Write *RAW* (pre-filter) outputs. *Disabled by default to bound disk; the RAW long table is unconsumed downstream.* |
| `emit_filtered` | bool | `true` | ZN/ZT | Write *FILTERED* outputs. |
| `write_long` | bool | `true` | ZN/ZT | Emit the long TSV (one row per site × sample × transcript × mod). |
| `write_pivots` | bool | `true` | ZN/ZT | Emit per-gene × mod pivots (coverage, fraction, Nmod). |
| `write_raw_per_gene` | bool | `false` | ZN | Per-gene tables for *RAW*. |
| `write_filtered_per_gene` | bool | `true` | ZN | Per-gene tables for *FILTERED*. |
| `min_cov` | int | `0` | ZN/ZT | Zero `frac_modified` when `Nvalid_cov < min_cov` (row kept). |
| `gtf` | path | assembled GTF | ZN | GTF to map sites → genes (uses `zn_index` first, gene-exon fallback). |

### 5. Differential test (test_diffs)

Identifies sites that differ in modification stoichiometry **across transcripts
within the same gene locus**.

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `test_diffs.min_cov` | int | `20` | Minimum pooled coverage per ZN at a site to include that ZN; a site is tested only if ≥2 ZN pass. |
| `test_diffs.topk` | int | `10` | Top sites to plot. |
| `test_diffs.test` | str | `"auto"` | `auto` picks Fisher (2×2) or Chi-square (r×2); or force `fisher`/`chi2`. |
| `test_diffs.pseudocount` | float | `0.5` | Chi-square cell pseudocount (ignored for Fisher). |
| `test_diffs.alternative` | str | `"two-sided"` | Fisher alternative: `two-sided`/`greater`/`less`. |
| `test_diffs.gene_filter` | list/null | `null` | Optional gene_name whitelist. |
| `test_diffs.mod_filter` | list/null | `null` | Optional mod_code whitelist (e.g. `["a","m"]`). |

**How it works:** for each site `(gene_name, mod_code, chrom, start0, end0,
strand)`, counts are pooled across samples within each ZN (`Ncov = ΣNvalid_cov`,
`Nmod = ΣNmod`), a contingency table (rows = ZN, cols = `[Nmod, Nunmod]`) is
tested, and p-values are BH-adjusted (`p_adj_bh`). If no site has ≥2 covered
isoforms (e.g. single-isoform / mitochondrial genes), an empty result is written
and the pipeline continues.

### 6. Classification (classify_diffs)

Runs after `test_diffs`. For every significant between-isoform site (passing
`fdr` **and** `min_effect`, the “>10% absolute stoichiometry” rule) it assigns
one MECE structural category explaining *why* the isoforms differ, anchored to
the gene's longest-3′UTR isoform. Covers **all** detected `mod_code`s.

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `classify_diffs.enable` | bool | `true` | Run classification. |
| `classify_diffs.min_effect` | float | `0.10` | Minimum max\|Δ stoichiometry\|. |
| `classify_diffs.fdr` | float | `0.05` | Maximum BH-FDR. |
| `classify_diffs.mod_filter` | list/null | `null` | mod_code(s) to classify (`null` = all detected). |
| `classify_diffs.min_cov` | int | `0` | Extra per-isoform `Ncov` floor on JSON entries. |
| `classify_diffs.tes_tol` | int | `200` | TES match tolerance for architecture/APA calls (the run scripts set `25` to match the assembler). |
| `classify_diffs.inside_tol` | int | `50` | Last-exon acceptor match tolerance. |
| `classify_diffs.ejc_nt` | int | `150` | EJC suppression zone (nt). |
| `classify_diffs.intergenic_gap` | int | `1000` | Minimum gap (bp) to call `INTERGENIC_TERMINAL_EXON`. |
| `classify_diffs.figures` | bool | `true` | Render per-category isoform **architecture maps** (`__figs_by_category_arch/`) + 2-panel stoichiometry figures (`__figs_by_category/`). |
| `classify_diffs.figs_per_category` | int | `10` | Top sites (by effect) per category to plot. |

**Categories:** `IPA_UNIQUE`, `SPLICED_EXON_UNIQUE`, `LAST_EXON_DISTAL_ONLY`,
`IPA_SHARED_EJC`, `SPLICING_EJC`, `LAST_EXON_PROXIMAL_APA_FAVORED`,
`LAST_EXON_DISTAL_APA_FAVORED`, `ALTERNATIVE_LAST_EXON`,
`INTERGENIC_TERMINAL_EXON`, `SHARED_TERMINAL_EXON`, `SHARED_INTERNAL_EXON`,
`UNEXPLAINED_SHARED`, `HI_INTRONIC_ARTIFACT`, `UNCLASSIFIED`.

### 7. Genotype (optional)

Uses the cleaned tagged BAMs to add read-level SNP, modification, and haplotype
association layers. Enable with `--set genotype.enable=true`.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `enable` | bool | `false` | Turn on the genotype module. |
| `jobs` | int | `2` | BAMs processed in parallel (the BAM-heavy steps shard per BAM × chromosome — raise on compute nodes). |
| `min_alt_reads` | int | `4` | Minimum alt-base support for a candidate SNP. |
| `min_total_cov` | int | `8` | Minimum total depth at a candidate SNP. |
| `min_alt_frac` / `max_alt_frac` | float | `0.10` / `0.90` | Alt-allele fraction window for segregating SNPs. |
| `min_baseq` / `min_mapq` | int | `20` / `10` | Quality floors for discovery / molecule extraction. |
| `mod_sites_require_snp_link` | bool | `true` | Keep only mod sites whose context_key matches a candidate SNP (lossless for SNP↔mod / dependency / haplotype-mod; bounds memory on deep genome-wide data). |
| `mod_jobs` | int | `8` | Concurrent `modkit extract calls` in the mod-table step (each streams a chromosome — bounds memory). |
| `min_mod_site_cov` | int | `1` | Minimum aggregated mod-site coverage for SNP-mod testing. |
| `min_group_reads` | int | `4` | Minimum group support for association/dependency tests. |
| `min_haplotype_reads` | int | `4` | Minimum read support for a haplotype to be tested vs collapsed into `OTHER`. |
| `max_haplotype_snps` | int | `4` | Maximum SNPs per local haplotype block before chunking. |
| `test` | str | `"auto"` | `auto`/`fisher`/`chi2`. |
| `pseudocount` | float | `0.5` | Pseudocount for non-2×2 chi-square. |

The BAM-heavy genotype steps shard the work unit per **(BAM × chromosome)** for
parallelism beyond the sample count, and producer tables are written in a
deterministic order so the haplotype blocks and dependency outputs are
reproducible run-to-run.

### 8. Report

`results/report/<prefix>_report.html` is self-contained. Relevant knobs:

| Key | Default | Description |
|-----|--------:|-------------|
| `report.enable` | `true` | Generate the HTML report. |
| `report.top_transcripts` / `report.top_genes` | `20` | Rows shown in the various tables. |
| `report.max_diff_figs` | `6` | Differential-site figures embedded. |
| `report.max_class_figs_per_category` | `10` | Per-category classification figures embedded. |
| `report.max_snp_figs` | `12` | Per-example SNP/haplotype figures per genotype section. |

The Site-Classification section shows **one structural-category distribution
graph per detected modification** (m6A, 5mC, pseudoU, inosine, the 2′-O-methyls,
…) plus a combined overview, the per-category top sites, and isoform
architecture maps. **All figures are written to a sidecar
`<prefix>_report_files/` directory and referenced by relative path** (not inlined
as base64), so the HTML stays small and opens in any browser; a built-in
lightbox enlarges figures on click. Haplotype-block tables include gene names,
genomic `region`, `span_bp`, and per-SNP `chrom:pos ref>alt` coordinates.

---

## Outputs

All paths are under `results/`.

**Assembly** — read tags `ZN` (partition index), `ZG` (gene index), `ZM`
(metagene index), `ZT` (`gene_name.gene_id.G{ZG}.T{tx}`):
- `assemble/<prefix>.gtf` — assembled isoform models
- `assemble/<prefix>_metrics.tsv` — isoform-level assembly metrics
- `assemble/<prefix>_tx_counts.tsv`, `…_tx_counts.pca.png` — transcript × sample counts + PCA
- `assemble/<prefix>_per_sample_stats.tsv`, `…_per_sample_read_stats.tsv` — per-sample summaries + read-retention funnel
- `assemble/<prefix>_tx_assigned_read_lengths.tsv`, `…_partition_map.tsv`
- `assemble/zt_tagged/*.bam`, `zt_filtered/*.bam` (cleaned, used downstream), `zt_scrap/*.bam`, `zt_bams/*.bam` (optional)

**Aggregation** (`aggregate_zn/`): `<prefix>_FILTERED_sites_long.tsv` (and
`_RAW_sites_long.tsv` when `emit_raw=true`); `<prefix>_FILTERED__per_gene_mod/`
pivots. FILTERED long tables keep every ZN/sample row at kept sites.

**Differential & classification** (`test_diffs/`):
- `<prefix>__ZN_site_diff_results.tsv` — per-site between-ZN test (`effect`, `p_adj_bh`, `per_transcript_json`)
- `<prefix>__ZN_site_classified.tsv` — structural category per significant site (high/low isoform, 3′ architecture, junction distances)
- `<prefix>__figs/`, `__figs_by_category_arch/{CATEGORY}/…png`, `__figs_by_category/{CATEGORY}/…png`

**Genotype** (`genotype/`, when enabled):
- `<prefix>_candidate_snps.tsv` — segregating SNP candidates
- `<prefix>_molecule_snps.tsv` — one row per read per candidate SNP
- `<prefix>_candidate_mod_sites.tsv`, `<prefix>_molecule_mod_calls.tsv`
- `<prefix>_snp_transcript_assoc.tsv`, `<prefix>_snp_mod_assoc.tsv`
- `<prefix>_snp_tx_mod_dependency.tsv` — transcript-conditioned SNP-mod dependency
- `<prefix>_haplotype_blocks.tsv` (gene_names / region / span_bp / snp_coords), `<prefix>_haplotype_transcript_assoc.tsv`, `<prefix>_haplotype_mod_assoc.tsv`
- `<prefix>__snp_figs/` — per-example figures

**Report**: `report/<prefix>_report.html` + `report/<prefix>_report_files/`.

**Timing**: `stage_timings.tsv` (per-stage seconds + `TOTAL`) is written at the
end of every run.

---

## Running on an HPC cluster (Slurm)

modulator is built to run **one isolated project directory per sample-set** and
to **checkpoint between stages**, which makes Slurm execution and re-runs simple.

A self-contained example wrapper:

```bash
#!/usr/bin/env bash
#SBATCH --cpus-per-task=48 --mem=128G --time=12:00:00 --job-name=mod_run
set -uo pipefail
RUNDIR=/path/to/runs/my_run          # has config/ (a copy) + symlinks to workflow/ src/ resources/
cd "$RUNDIR"
export TMPDIR="$RUNDIR/results/tmp"; mkdir -p "$TMPDIR"   # keep temp on the shared FS, not node /tmp
[ "${CLEAN:-0}" = 1 ] && rm -rf "$RUNDIR/results"          # CLEAN=1 -> from-scratch (fully timed) run
t0=$SECONDS
micromamba run -n modulator modulator run ${RESUME:+--resume} \
  --workdir "$RUNDIR" --config config/config.yaml --jobs 48 \
  --set reference_fa=ref.fa reference_gtf=ref.gtf \
        bams_dir=/path/to/bams bam_glob='*.bam' prefix=my_run \
        genotype.enable=true genotype.jobs=48
echo "wall_seconds: $((SECONDS - t0))"
```

Recommendations:

- **One dir per run**: copy `config/` into each run dir (so per-run tweaks don't
  collide) and symlink `workflow/`, `src/`, `resources/` to a single canonical
  clone.
- **Resume**: submit with `RESUME=1` to skip completed stages after a failure or
  to re-run only the cheap tail (`--stages classify_diffs,report`).
- **From-scratch timing**: `CLEAN=1` wipes `results/` first so `stage_timings.tsv`
  reflects a full run. Collect totals from `sacct` (Elapsed) or the
  `wall_seconds` log line; per-stage from each `results/stage_timings.tsv`; read
  counts from `…_per_sample_read_stats.tsv` (`total_reads_bam` processed,
  `considered_reads` kept, `total − considered` filtered).
- **Sizing**: `--jobs` / `threads` / `modkit.common.threads` are independent;
  keep their product within `--cpus-per-task`. `genotype.jobs` parallelizes the
  BAM-heavy SNP/mod steps separately.

## Performance & troubleshooting

- **Memory**: the assembler streams BAMs and drops consumed page cache
  (`posix_fadvise`), so its working set is small; OS page cache can still inflate
  cgroup-accounted RSS on shared nodes — give generous `--mem` (it is advisory
  where `ConstrainRAMSpace=no`).
- **TMPDIR**: point it at the shared filesystem, not small node-local `/tmp`
  (the aggregation sort/temp can be large).
- **Disk / quota**: `aggregation.zn.emit_raw=false` (default) keeps the big RAW
  intermediates off disk; report images are externalized to a sidecar folder.
- **modkit**: pinned to 0.5.0 (`--partition-tag` was removed in 0.6.x).
- **Single-isoform / mitochondrial genes**: set
  `assembler.min_introns=0 assembler.min_introns_read=0`; `test_diffs` emits an
  empty result (nothing to contrast) and the run continues to genotype/report.
- **Determinism**: genotype producer tables are sorted before writing, so
  haplotype blocks and dependency outputs are byte-reproducible across runs.

## Recent changes

- **Report**: one structural-category graph **per modification** (not a single
  hard-labeled chart); images **externalized** to a sidecar folder + lightbox
  (small HTML that opens in Chrome); haplotype blocks annotated with gene names
  and genomic coordinates; corrected modification-code labels.
- **Genotype**: per-(BAM × chromosome) sharding; deterministic producer ordering;
  SNP-linked mod-site filter (`mod_sites_require_snp_link`) and `mod_jobs` to
  bound memory at genome scale.
- **Aggregation**: streaming per-chromosome ZN aggregator (`aggregation.engine:
  stream`); RAW outputs off by default.
- **Pipeline**: stage-level checkpointing (`--resume`); per-stage timing
  (`results/stage_timings.tsv`); intronless-transcriptome support; `test_diffs`
  no longer aborts when there is nothing to test.
