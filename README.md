# modulator: Transcript-Specific Modification Calling for Nanopore Direct-RNA Sequencing Data

A Snakemake pipeline for analyzing RNA modifications from BAM files.

## Installation

1. Set up micromamba in your HPC environment: https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html
2. Set up a custom micromamba environment for modulator: `micromamba create -y -n modulator -c conda-forge -c bioconda python=3.13.7 pandas=2.3.3 numpy=2.3.3 matplotlib=3.10.6 pysam=0.23.3 samtools=1.22.1 scipy snakemake`  
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
    min_cov=5 \
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

## Outputs

- **Tagging:**  
  - `ZN`: transcript index within gene (1..k)  
  - `ZG`: gene index (run deterministic)  
  - `ZT`: string label of the form `"gene_name.gene_id.G{ZG}.T{ZN}"`  

- **Output Summary:**  
  - `<prefix>_tx_counts.tsv`: transcript × sample read counts  
  - `<prefix>_tx_counts.pca.png`: PCA of samples (log1p counts)  
  - `<prefix>_per_sample_stats.tsv`: summary per sample (reads, transcripts, median per transcript)  
  - `zt_tagged/*.bam`: one tagged BAM per sample  
  - `zt_bams/*.bam`: per-transcript BAMs (optional) 

## Citation

## Contributors 

* [Theodore Nelson](https://github.com/Theo-Nelson), Weill Cornell Medicine
* [Michael Goneos](https://github.com/mgoneos), Weill Cornell Medicine 

