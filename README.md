![modulator](workflow/images/modulator_banner.png)



A Snakemake pipeline for analyzing RNA modifications from aligned BAM files.

## Installation

1. Set up micromamba in your HPC environment: https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html
2. Set up a custom micromamba environment for modulator: `micromamba create -y -n modulator -c conda-forge -c bioconda python=3.13.7 pandas=2.3.3 numpy=2.3.3 matplotlib=3.10.6 pysam=0.23.3 samtools=1.22.1 scipy snakemake ont-modkit=0.5.0`  
3. Activate the environment: `micromamba activate modulator`
4. Clone the repository: `git clone https://github.com/Theo-Nelson/modulator.git`
5. Configure `config/config.yaml` with your samples and references.
6. Run `snakemake` from the workflow directory.

## Input File Requirements

## Pipeline Parameters and Usage

Minimal command to run on demo data:

```bash
snakemake -j 8 \
  --configfile ../config/config.yaml \
  --config \
    reference_fa=/path/to/ref.fa \
    reference_gtf=/path/to/ref.gtf \
  --rerun-incomplete --printshellcmds
```

Expanded command to specify samples and references:

```bash
snakemake -j 8 \
  --configfile ../config/config.yaml \
  --config \
    reference_fa=/path/to/ref.fa \
    reference_gtf=/path/to/ref.gtf \
    bams_dir=/path/to/your/bams \
    bam_glob='*.bam' \
    prefix=fivegenes_readbacked_annot \
    threads=64 \
    mods='["17596","a","m","17802","69426","19228","19229","19227"]' \
    ref_bases='["A","A","C","T","A","C","G","T"]' \
    min_cov=0 \
    min_cov_test=20 \
    topk=10 \
    assembler='{primary_only: true,
                min_mapq: 10,
                min_introns_read: 1,
                require_softclip3p: 0,
                apa_window: 20,
                tes_window: null,
                min_reads: 40,
                min_frac: 0.00,
                min_introns: 1,
                min_polya_length: 12,
                min_polya_purity: 0.5,
                polya_support_frac: 0.5,
                tes_match_tol: 25,
                exact_tes_tol: 10,
                write_zt_bams: false,
                write_zt_tagged_sample_bams: true,
                emit_modkit_manifest: false,
                min_reads_per_sample_for_mod: 5,
                min_total_reads_for_mod: 20}' \
    modkit='{common: {
                log_file_template: "results/{which}/{sample}.log",
                region: null,
                max_depth: 1000,
                include_bed: null,
                include_unmapped: false,
                edge_filter: null,
                invert_edge_filter: false,
                threads: "{threads}",
                interval_size: 100000,
                queue_size: 1000,
                chunk_size: null,
                num_reads: 10042,
                sampling_frac: null,
                seed: null,
                sample_region: null,
                sampling_interval_size: 1000000,
                no_filtering: false,
                filter_thresholds: ["A:0.8","C:0.8","G:0.8","T:0.8"],
                mod_thresholds: ["17596:0.99","a:0.99","m:0.99","17802:0.99","69426:0.99","19228:0.99","19229:0.99","19227:0.99"],
                ignore: [],
                force_allow_implicit: false,
                motif: [],
                cpg: false,
                ref_mask: false,
                combine_mods: false,
                combine_strands: false,
                mixed_delim: false,
                only_tabs: false,
                bedgraph: false,
                header: false,
                prefix: null,
                suppress_progress: true
             },
             zn: { partition_tag: "ZN", per_mod_bed: true },
             zt: { partition_tag: "ZT" }}' \
    aggregation='{zn: { filter_enable: true,
                        count_diff_factor: 3,
                        mod_fail_margin: 1,
                        emit_raw: true,
                        emit_filtered: true,
                        write_long: true,
                        write_pivots: true,
                        write_raw_per_gene: true,
                        write_filtered_per_gene: true },
                  zt: { filter_enable: true,
                        count_diff_factor: 3,
                        mod_fail_margin: 1,
                        emit_raw: true,
                        emit_filtered: true,
                        write_long: true,
                        write_pivots: true }}' \
    test_diffs='{min_cov: 20,
                 topk: 10,
                 test: "auto",
                 pseudocount: 0.5,
                 alternative: "two-sided",
                 gene_filter: null,
                 mod_filter: null}' \
  --rerun-incomplete --printshellcmds
```

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
| `write_zt_bams` | boolean | `false` | `true` / `false` | Write one BAM per transcript per sample (ZT-tagged). Useful for per-transcript analyses; produces many files. |
| `write_zt_tagged_sample_bams` | boolean | `true` | `true` / `false` | Write one ZT/ZN-tagged BAM per sample (all reads). Enables downstream modification calling. |
| `emit_modkit_manifest` | boolean | `false` | `true` / `false` | Also write a manifest (`zt_bams/modkit_manifest.tsv`) for modkit processing. |
| `min_reads_per_sample_for_mod` | integer | `5` | `1–50` | Minimum per-sample read support for creating a per-transcript BAM. |
| `min_total_reads_for_mod` | integer | `20` | `10–200` | Minimum total read support (across samples) for a transcript to be eligible for modkit BAM output. |

A few notes on how these parameters cooperate:

- **TES logic:** `tes_window` (if not `null`) overrides `apa_window` for 3′ end clustering.  
- **Poly(A) evidence:** A read provides poly(A/T) support if its 3′ soft-clip length ≥ `min_polya_length` **and** purity ≥ `min_polya_purity`.  
- **Filtering:** Isoforms must pass *all* filters (`min_reads`, `min_frac`, `min_introns`, `polya_support_frac`) to be retained.  

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
| `threads` | int/str | `"{threads}"` | `1+` | `-t/--threads` |
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
| `zn` | `partition_tag` | `"ZN"` | Partition by ZN (transcript index per gene) |
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
| `mod_fail_margin`          | int     | `1`     | `0–5`                    | ZN, ZT | Additional margin on `Nfail` for the `Nmod` fail rule. |
| `emit_raw`                 | bool    | `true`  | `true`/`false`           | ZN, ZT | Write *RAW* outputs (pre-filter). |
| `emit_filtered`            | bool    | `true`  | `true`/`false`           | ZN, ZT | Write *FILTERED* outputs (site-kept logic). |
| `write_long`               | bool    | `true`  | `true`/`false`           | ZN, ZT | Emit the long TSV (one row per site × sample × transcript × mod). |
| `write_pivots`             | bool    | `true`  | `true`/`false`           | ZN, ZT | Emit per‑gene × mod pivoted tables (coverage, fraction, Nmod). |
| `write_raw_per_gene`       | bool    | `false` | `true`/`false`           | ZN     | Also write per‑gene tables for *RAW*. |
| `write_filtered_per_gene`  | bool    | `true`  | `true`/`false`           | ZN     | Also write per‑gene tables for *FILTERED*. |
| `min_cov`                  | int     | `0`     | `0–20`                   | ZN, ZT | If `Nvalid_cov < min_cov`, set `frac_modified = 0` (row kept). Does **not** affect pass/fail. |
| `out_prefix`               | str     | —       | path                     | ZN, ZT | Prefix for all output files (both RAW and FILTERED variants). |
| `gtf`                      | path    | —       | path to GTF              | ZN     | GTF used to map sites → genes (union-of-exons span per gene). |

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
  - `ZN`: transcript index within gene (1..k)  
  - `ZG`: gene index (run deterministic)  
  - `ZT`: string label of the form `"gene_name.gene_id.G{ZG}.T{ZN}"`  

- **Output Summary:**  
  - `<prefix>_tx_counts.tsv`: transcript × sample read counts  
  - `<prefix>_tx_counts.pca.png`: PCA of samples (log1p counts)  
  - `<prefix>_per_sample_stats.tsv`: summary per sample (reads, transcripts, median per transcript)  
  - `zt_tagged/*.bam`: one tagged BAM per sample  
  - `zt_bams/*.bam`: per-transcript BAMs (optional)

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

