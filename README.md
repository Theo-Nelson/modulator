![modulator](workflow/images/modulator_banner.png)


A Snakemake pipeline for analyzing RNA modifications from aligned BAM files.

## Installation

1. Set up micromamba in your HPC environment: https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html
2. Clone the repository: `git clone https://github.com/Theo-Nelson/modulator.git`
3. Create the environment from the bundled YAML:

```bash
micromamba env create -n modulator -f workflow/envs/modulator.yaml
```

4. Activate the environment:

```bash
micromamba activate modulator
```

5. Install the package in editable mode so the `modulator` command is available:

```bash
python -m pip install -e .
```

## Input File Requirements

## Pipeline Parameters and Usage

Package-native quick demo on the bundled MXD1 slice:

```bash
modulator demo \
  --reference-fa /path/to/ref.fa \
  --reference-gtf /path/to/ref.gtf
```

Standard full run with explicit config overrides:

```bash
modulator run \
  --config config/config.yaml \
  --jobs 8 \
  --set \
    reference_fa=/path/to/ref.fa \
    reference_gtf=/path/to/ref.gtf \
    bams_dir=/path/to/your/bams \
    bam_glob='*.bam' \
    prefix=my_run
```

If you prefer not to install the console script yet, the exact same interface is available through:

```bash
PYTHONPATH=src python -m modulator run --config config/config.yaml --set reference_fa=/path/to/ref.fa reference_gtf=/path/to/ref.gtf
```

Useful CLI helpers:

```bash
modulator validate-config --config config/config.yaml --set reference_fa=/path/to/ref.fa reference_gtf=/path/to/ref.gtf
modulator run --config config/config.yaml --stages assemble,read_stats,multigene_filter
modulator demo --reference-fa /path/to/ref.fa --reference-gtf /path/to/ref.gtf --dataset ALCAM_NHSL1_SERAC1_MXD1_RIOK3_reads --mode full
```

A few notes on the new interface:

- `config/config.yaml` is now the normal place for nested settings like `assembler`, `modkit`, `aggregation`, `genotype`, and `report`.
- CLI overrides are intentionally simple `key=value` or `nested.key=value` items passed with `--set`.
- `modulator demo` defaults to a fast MXD1-only run and keeps genotype disabled unless you override it.
- The package runner executes the existing workflow scripts directly, so the algorithms stay the same while the launch UX is much more robust on systems like ACES.

### Assembly Parameters

The following table explains the different parameter functions available for transcript assembly. 

| Parameter | Type | Default | Typical Range / Options | Description |
|------------|------|----------|--------------------------|--------------|
| `primary_only` | boolean | `true` | `true` / `false` | Use only primary alignments (ignore secondary/supplementary). Keeps one best alignment per read for assembly. |
| `min_mapq` | integer | `10` | `0–60` | Minimum MAPQ (mapping quality) required to keep a read. Filters low-confidence alignments. |
| `min_introns_read` | integer | `1` | `0–10` | Require at least this many introns in a read. Filters unspliced reads unless `0`. |
| `require_softclip3p` | integer (nt) | `0` | `0–40` | Minimum 3′ soft-clip length for a read to be included (proxy for poly(A/T) tail evidence). |
| `apa_window` | integer (nt) | `20` | `10–100` | Window size around 3′ ends to cluster reads into APA (alternative polyadenylation) groups. |
| `tes_window` | integer (nt) or `null` | `null` | `null` or numeric | Overrides `apa_window` if provided; used for TES clustering or matching tolerance. |
| `min_reads` | integer | `40` | `1–200` | Minimum total reads supporting an isoform to retain it. |
| `min_frac` | float | `0.00` | `0.0–0.2` | Minimum fraction of all reads supporting an isoform. Acts as a global rarity filter. |
| `min_introns` | integer | `1` | `0–10` | Minimum introns per *assembled transcript*. Filters single-exon transcripts unless set to `0`. |
| `min_polya_length` | integer (nt) | `12` | `8–30` | Minimum soft-clip length required for poly(A/T) evidence. |
| `min_polya_purity` | float | `0.5` | `0.4–0.95` | Fraction of A/T bases in the soft-clipped tail to count as poly(A/T). |
| `polya_support_frac` | float | `0.5` | `0.3–0.9` | Minimum fraction of reads in an isoform with poly(A/T) support. |
| `tes_match_tol` | integer (nt) | `25` | `5–50` | Tolerance (bp) for matching transcript end sites to reference TES. |
| `exact_tes_tol` | integer (nt) | `10` | `2–20` | Distance threshold to label a match as `EXACT` (vs. `NOVEL_APA`). |
| `assignment_mode` | string | `"support_first"` | `"support_first"` / `"longest_first"` | `support_first` ranks TES-cluster canonicals by direct distal 5′ splice-junction support; `longest_first` reproduces the previous longest-chain-first suffix collapse. |
| `zn_mode` | string | `"metagene_colored"` | `"metagene_colored"` / `"gene_local"` | `metagene_colored` assigns ZN by greedy coloring within overlapping metagenes so non-overlapping transcripts can reuse partition indices. `gene_local` reproduces the previous per-gene transcript index semantics. |
| `min_distal_anchor_reads` | integer | `2` | `1–20` | Minimum number of exact-chain reads needed before a longer chain can absorb shorter suffix-compatible chains in `support_first` mode. |
| `min_distal_anchor_frac` | float | `0.05` | `0.0–1.0` | Minimum exact-chain fraction within a suffix family required for longer-chain absorption in `support_first` mode. |
| `min_exact_canonical_reads` | integer | `1` | `1–20` | Minimum exact-chain read count for any canonical to participate in suffix collapse. |
| `write_zt_bams` | boolean | `false` | `true` / `false` | Write one BAM per transcript per sample (ZT-tagged). Useful for per-transcript analyses; produces many files. |
| `write_zt_tagged_sample_bams` | boolean | `true` | `true` / `false` | Write one ZT/ZN-tagged BAM per sample (all reads). Enables downstream modification calling. |
| `emit_modkit_manifest` | boolean | `false` | `true` / `false` | Also write a manifest (`zt_bams/modkit_manifest.tsv`) for modkit processing. |
| `min_reads_per_sample_for_mod` | integer | `5` | `1–50` | Minimum per-sample read support for creating a per-transcript BAM. |
| `min_total_reads_for_mod` | integer | `20` | `10–200` | Minimum total read support (across samples) for a transcript to be eligible for modkit BAM output. |
| `status_every` | integer | `0` | `0+` | Print assembly progress every N reads (`0` disables status logging). |

A few notes on how these parameters cooperate:

- **TES logic:** `tes_window` (if not `null`) overrides `apa_window` for 3′ end clustering.  
- **Support-first assignment:** in `support_first` mode, a longer suffix-compatible chain can absorb shorter chains only when it has direct exact-chain support from its own distal unique 5′ splice structure.  
- **Poly(A) evidence:** A read provides poly(A/T) support if its 3′ soft-clip length ≥ `min_polya_length` **and** purity ≥ `min_polya_purity`.  
- **Filtering:** Isoforms must pass *all* filters (`min_reads`, `min_frac`, `min_introns`, `polya_support_frac`) to be retained.  
- **ZN semantics:** in `metagene_colored` mode, `ZN` is a metagene-aware partition index rather than a simple transcript index within one gene.

### Multigene Filter Parameters

This optional post-assembly module scans each `zt_tagged` BAM, summarizes overlap against assembled gene exon unions, and writes a cleaned BAM for downstream `modkit`. In the default `resolve` mode, overlapping reads are retained and categorized instead of being discarded.

| Parameter | Type | Default | Description |
|------------|------|----------|--------------|
| `enable` | boolean | `true` | Run the multigene-overlap filter between assembly and `modkit`. |
| `mode` | string | `"resolve"` | `"resolve"` keeps overlapping reads and records how they were resolved; `"legacy_scrap"` reproduces the previous behavior of removing multi-gene overlaps. |
| `zero_gene_action` | string | `"keep"` | Whether reads overlapping zero assembled genes should stay in the cleaned BAM or be sent to scrap (`keep` / `scrap`). |

### Modkit Parameters

The following table explains the different parameter functions available for modkit, which map to the parameters in the ONT's [Advanced Usage Guide](https://github.com/nanoporetech/modkit/blob/master/book/src/advanced_usage.md). 

| Key | Type | Default | Range/Options | Maps to |
|---|---|---:|---|---|
| `log_file_template` | str | `"results/{which}/{sample}.log"` | any path | `--log-filepath` |
| `region` | str/null | `null` | `chr`, `chr:start-end` | `--region` |
| `max_depth` | int | `1000` | `1–2,147,483,647` | `--max-depth` |
| `include_bed` | str/null | `null` | path | `--include-bed` |
| `include_unmapped` | bool | `false` | `true/false` | `--include-unmapped` |
| `edge_filter` | str/int/null | `null` | `N` or `"N,M"` | `--edge-filter` |
| `invert_edge_filter` | bool | `false` | `true/false` | `--invert-edge-filter` |
| `threads` | int | top-level `threads` config | `1+` | `-t/--threads` |
| `interval_size` | int | `100000` | `1+` | `--interval-size` |
| `queue_size` | int | `1000` | `1+` | `--queue-size` |
| `chunk_size` | int/null | `null` | `1+` | `--chunk-size` |
| `num_reads` | int | `10042` | `1+` | `--num-reads` |
| `sampling_frac` | float/null | `null` | `0–1` | `--sampling-frac` |
| `seed` | int/null | `null` | any int | `--seed` |
| `sample_region` | str/null | `null` | `chr`/`chr:start-end` | `--sample-region` |
| `sampling_interval_size` | int | `1000000` | `1+` | `--sampling-interval-size` |
| `no_filtering` | bool | `false` | `true/false` | `--no-filtering` |
| `filter_percentile` | float/null | (not an option in modulator) | `0–1` | `--filter-percentile` |
| `filter_thresholds` | list[str] | `["A:0.8","C:0.8","G:0.8","T:0.8"]` | per-base | `--filter-threshold` (repeatable) |
| `mod_thresholds` | list[str] | eight bases at `0.99` | per-mod | `--mod-threshold` (repeatable) |
| `ignore` | list[str] | `[]` | e.g., `["h"]` | `--ignore` (repeatable) |
| `force_allow_implicit` | bool | `false` | `true/false` | `--force-allow-implicit` |
| `motif` | list[str] | `[]` | e.g., `["CG:0","CGCG:2"]` | `--motif <motif> <offset>` (repeatable) |
| `cpg` | bool | `false` | `true/false` | `--cpg` |
| `ref_mask` | bool | `false` | `true/false` | `--mask` |
| `combine_mods` | bool | `false` | `true/false` | `--combine-mods` |
| `combine_strands` | bool | `false` | `true/false` | `--combine-strands` |
| `mixed_delim` | bool | `false` | `true/false` | `--mixed-delim` |
| `only_tabs` | bool | `false` | `true/false` | `--only-tabs` |
| `bedgraph` | bool | `false` | `true/false` | `--bedgraph` |
| `header` | bool | `false` | `true/false` | `--header` |
| `prefix` | str/null | `null` | any | `--prefix` (bedGraph mode) |
| `suppress_progress` | bool | `true` | `true/false` | `--suppress-progress` |

These partition tags are hard-coded into the pipeline to split reads according to their transcript assignments. 

| Section | Key | Default | Meaning |
|---|---|---:|---|
| `zn` | `partition_tag` | `"ZN"` | Partition by ZN (metagene-aware partition index in `metagene_colored` mode; legacy per-gene transcript index in `gene_local` mode) |
| `zn` | `per_mod_bed` | `true` | Emit one BED per mod (pipeline behavior) |
| `zt` | `partition_tag` | `"ZT"` | Partition by ZT (gene+tx human-readable code) |

### Aggregation Parameters

The parameters below control how *modulator* aggregates **ZN** (per-transcript index) and **ZT** (per-transcript code) modkit bedMethyl outputs into long tables and per‑gene/per‑mod pivots.

> **Site keeping rule (when `filter_enable` is true):**  
> A row **fails** if `(Ndiff > count_diff_factor * Nvalid_cov)` **or** `(Nmod <= Nfail + mod_fail_margin)`.  
> A **site is kept** in *FILTERED* outputs if **any** row at that site passes (considering all samples and transcripts). **When a site is kept, _all rows_ for that site are retained** (i.e., you keep every ZN/sample line for that genomic position+mod).  
> `min_cov` affects only the displayed `frac_modified` (set to 0 when `Nvalid_cov < min_cov`) and **does not** influence pass/fail logic.

| Parameter                   | Type    | Default | Typical Range / Options | Scope | Description |
|----------------------------|---------|---------|--------------------------|-------|-------------|
| `filter_enable`            | bool    | `true`  | `true`/`false`           | ZN, ZT | Enable site-level filtering using the rule above. |
| `count_diff_factor`        | float   | `3.0`   | `1–10`                   | ZN, ZT | Threshold factor for the `Ndiff` term in the fail rule. |

### Genotype Parameters

This optional module uses the same cleaned tagged BAMs as the rest of the pipeline and adds read-level SNP, mod, and haplotype association layers.

| Parameter | Type | Default | Description |
|------------|------|----------|--------------|
| `enable` | boolean | `false` | Turn on genotype-aware molecule tables and association tests. |
| `min_alt_reads` | integer | `4` | Minimum alternative-base read support for a candidate SNP. |
| `min_total_cov` | integer | `8` | Minimum total depth for a candidate SNP position. |
| `min_alt_frac` | float | `0.10` | Minimum alternative allele fraction for candidate SNP discovery. |
| `max_alt_frac` | float | `0.90` | Maximum alternative allele fraction for sites to remain in the segregating SNP set used for association tests. Sites above this are treated as near-fixed reference-discordant sites rather than informative segregating SNPs. |
| `min_baseq` | integer | `20` | Minimum base quality for SNP molecule extraction. |
| `min_mapq` | integer | `10` | Minimum alignment MAPQ for SNP discovery and molecule extraction. |
| `min_mod_site_cov` | integer | `1` | Minimum aggregated mod-site coverage needed before a site is used in SNP-mod testing. |
| `min_group_reads` | integer | `4` | Minimum group support used by SNP-transcript, SNP-mod, and conditional dependency tests. |
| `min_haplotype_reads` | integer | `4` | Minimum read support for a haplotype to be tested directly instead of collapsed into `OTHER`. |
| `max_haplotype_snps` | integer | `4` | Maximum number of SNPs retained in one local read-backed haplotype block before chunking. |
| `test` | string | `"auto"` | Association test mode: `"auto"`, `"fisher"`, or `"chi2"`. |
| `pseudocount` | float | `0.5` | Pseudocount used for non-2x2 chi-square contingency tests. |

The genotype module writes:

- `results/genotype/<prefix>_candidate_snps.tsv`: discovered segregating SNP candidates with read support and locus annotations
- `results/genotype/<prefix>_molecule_snps.tsv`: one row per read per candidate SNP
- `results/genotype/<prefix>_candidate_mod_sites.tsv`: modulator-derived candidate mod sites used for SNP-mod testing
- `results/genotype/<prefix>_molecule_mod_calls.tsv`: one row per read per candidate mod site
- `results/genotype/<prefix>_snp_transcript_assoc.tsv`: SNP to transcript association results
- `results/genotype/<prefix>_snp_mod_assoc.tsv`: SNP to epitranscriptome association results
- `results/genotype/<prefix>_snp_tx_mod_dependency.tsv`: transcript-conditioned SNP-mod dependency results
- `results/genotype/<prefix>_haplotype_blocks.tsv`: local read-backed haplotype blocks
- `results/genotype/<prefix>_haplotype_transcript_assoc.tsv`: haplotype to transcript association results
- `results/genotype/<prefix>_haplotype_mod_assoc.tsv`: haplotype to epitranscriptome association results
| `mod_fail_margin`          | int     | `1`     | `0–5`                    | ZN, ZT | Additional margin on `Nfail` for the `Nmod` fail rule. |
| `emit_raw`                 | bool    | `true`  | `true`/`false`           | ZN, ZT | Write *RAW* outputs (pre-filter). |
| `emit_filtered`            | bool    | `true`  | `true`/`false`           | ZN, ZT | Write *FILTERED* outputs (site-kept logic). |
| `write_long`               | bool    | `true`  | `true`/`false`           | ZN, ZT | Emit the long TSV (one row per site × sample × transcript × mod). |
| `write_pivots`             | bool    | `true`  | `true`/`false`           | ZN, ZT | Emit per‑gene × mod pivoted tables (coverage, fraction, Nmod). |
| `write_raw_per_gene`       | bool    | `false` | `true`/`false`           | ZN     | Also write per‑gene tables for *RAW*. |
| `write_filtered_per_gene`  | bool    | `true`  | `true`/`false`           | ZN     | Also write per‑gene tables for *FILTERED*. |
| `min_cov`                  | int     | `0`     | `0–20`                   | ZN, ZT | If `Nvalid_cov < min_cov`, set `frac_modified = 0` (row kept). Does **not** affect pass/fail. |
| `out_prefix`               | str     | —       | path                     | ZN, ZT | Prefix for all output files (both RAW and FILTERED variants). |
| `gtf`                      | path    | —       | path to GTF              | ZN     | GTF used to map sites → genes using transcript-aware `zn_index` assignments first, with gene-exon fallback. |

### Differential Site-Level Per-Transcript Within Locus Test Parameters

These parameters allow for the identification of sites with differences across transcripts within the same gene locus. 

| Key | Type | Default | Description |
|---|---|---|---|
| `test_diffs.min_cov` | int | `20` | Minimum **pooled** coverage per ZN at a site to include that ZN in testing. A site is tested only if ≥2 ZN pass. (This mirrors `min_cov_test`.) |
| `test_diffs.topk` | int | `10` | Number of top sites to plot as figures. |
| `test_diffs.test` | str | `"auto"` | `"auto"` chooses Fisher (2×2) or Chi-square (r×2) automatically; `"fisher"` forces Fisher (requires exactly 2 ZN); `"chi2"` forces Chi-square for any r≥2. |
| `test_diffs.pseudocount` | float | `0.5` | Pseudocount added to each cell for Chi-square stability (ignored for Fisher). |
| `test_diffs.alternative` | str | `"two-sided"` | Fisher alternative hypothesis: `"two-sided"`, `"greater"`, or `"less"`. |
| `test_diffs.gene_filter` | list[str] / null | `null` | Optional gene_name whitelist. |
| `test_diffs.mod_filter` | list[str] / null | `null` | Optional mod_code whitelist (e.g., `["a","m"]`). |

**How the test works:** For each site `(gene_name, mod_code, chrom, start0, end0, strand)`, counts are **pooled across samples** within each ZN: `Ncov = Σ Nvalid_cov`, `Nmod = Σ Nmod`, `Nunmod = Ncov − Nmod`. Build a contingency table (rows=ZN, cols=[Nmod, Nunmod]) and test differences in stoichiometry across ZN; adjust p-values with Benjamini–Hochberg (`p_adj_bh`).

## Outputs

### Assembly Outputs

- **Within-Bam File Tagging:**  
  - `ZN`: metagene-aware partition index (or legacy per-gene transcript index if `zn_mode: gene_local`)  
  - `ZG`: gene index (run deterministic)  
  - `ZM`: metagene index  
  - `ZT`: string label of the form `"gene_name.gene_id.G{ZG}.T{gene-local transcript index}"`  

- **Output Summary:**  
  - `<prefix>_metrics.tsv`: isoform-level assembly metrics including exact-chain reads, truncation-assigned reads, anchor support, and assignment mode  
  - `<prefix>_tx_counts.tsv`: transcript × sample read counts  
  - `<prefix>_tx_counts.pca.png`: PCA of samples (log1p counts)  
  - `<prefix>_per_sample_stats.tsv`: summary per sample (reads, transcripts, median per transcript)  
  - `<prefix>_tx_assigned_read_lengths.tsv`: mean/median/min/max assigned read length per transcript  
  - `<prefix>_partition_map.tsv`: transcript ↔ gene/metagene/ZN lookup table for downstream interpretation  
  - `<prefix>_multigene_scrap_tx_counts.tsv`: transcript × sample counts for reads removed by the multigene filter  
  - `zt_tagged/*.bam`: one tagged BAM per sample  
  - `zt_filtered/*.bam`: cleaned tagged BAMs used for downstream `modkit` runs  
  - `zt_scrap/*.bam`: legacy-scrapped reads plus per-sample overlap-resolution detail TSVs  
  - `zt_bams/*.bam`: per-transcript BAMs (optional)  
  - `results/report/<prefix>_report.html`: HTML summary report

### Modkit Outputs

### Aggregation Outputs
 
- **ZN**:  
  - `{out_prefix}_RAW_sites_long.tsv` and `{out_prefix}_FILTERED_sites_long.tsv`  
  - `{out_prefix}_RAW__per_gene_mod/` and `{out_prefix}_FILTERED__per_gene_mod/` (pivoted and/or per‑gene files)
- **ZT**:  
  - `{out_prefix}_RAW_long.tsv`, `{out_prefix}_RAW_*_pivot.tsv`  
  - `{out_prefix}_FILTERED_long.tsv`, `{out_prefix}_FILTERED_*_pivot.tsv`

> **Note**: *FILTERED* long tables still include every ZN/sample row at kept sites, preserving per‑transcript stoichiometries while removing entire sites that fail across all rows.

## Citation

## Contributors 

* [Theodore Nelson](https://github.com/Theo-Nelson), Weill Cornell Medicine
* [Michael Goneos](https://github.com/mgoneos), Weill Cornell Medicine

This project was supported with computational resources from the National Science Foundation ACCESS Allocation Request BIO240371. 

T.M.N. was supported by a Medical Scientist Training Program grant from the National Institute of General Medical Sciences of the National Institutes of Health under award number: T32GM152349 to the Weill Cornell/Rockefeller/Sloan Kettering Tri-Institutional MD-PhD Program. 

## Version Changes

- TES clustering can now use a `support_first` assignment mode that ranks canonicals by direct distal 5' splice-junction support instead of always letting the longest compatible chain absorb shorter suffixes.
- `ZN` can now be assigned in a `metagene_colored` mode that reuses the same partition index for non-overlapping transcripts within an overlapping metagene, reducing `modkit` partition count without throwing overlapping-locus reads away.
- The multigene stage now supports `mode: resolve` so overlapping reads are retained and summarized, while `mode: legacy_scrap` reproduces the older discard behavior for comparisons.
- A workflow-integrated HTML report is now generated by default.
- An optional genotype-aware module can now discover read-supported SNPs, test SNP to transcript associations, test SNP to epitranscriptome associations, evaluate SNP-transcript-mod dependency, and build local read-backed haplotype blocks.

For a quick code-path sanity check without running the full workflow, you can run:

```bash
python workflow/scripts/regression_smoke_checks.py
python workflow/scripts/genotype_regression_smoke_checks.py
```


