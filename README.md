# modulator: Transcript-Specific Modification Calling for Nanopore Direct-RNA Sequencing Data

A Snakemake pipeline for analyzing RNA modifications from BAM files.

## Installation

1. Set up micromamba in your HPC environment: https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html
2. Set up a custom micromamba environment for modulator: ``micromamba create -y -n modulator -c conda-forge -c bioconda python=3.13.7 pandas=2.3.3 numpy=2.3.3 matplotlib=3.10.6 pysam=0.23.3 samtools=1.22.1 scipy snakemake` ` 
3. Activate the environment: `micromamba activate modulator`
4. Clone the repository: `git clone https://github.com/Theo-Nelson/modulator.git`
5. Install Snakemake and conda. 
6. Configure `config/config.yaml` with your samples and references.
7. Run `snakemake` from the workflow directory.

## Sample Preparation

## Pipeline Parameters and Usage

```bash
cd workflow
snakemake
```

## Citation

This pipeline assembles transcripts from BAM files, runs modkit pileup for modifications, aggregates results per gene and transcript, and tests for differences in modification stoichiometry between transcripts.

## Contributors 

