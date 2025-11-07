# Modulator

A Snakemake pipeline for analyzing RNA modifications from BAM files.

## Installation

1. Set up micromamba in your HPC environment: https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html
2. Clone the repository.
3. Install Snakemake and conda. 
4. Configure `config/config.yaml` with your samples and references.
5. Run `snakemake` from the workflow directory.

## Usage

```bash
cd workflow
snakemake
```

## Description

This pipeline assembles transcripts from BAM files, runs modkit pileup for modifications, aggregates results per gene and transcript, and tests for differences in modification stoichiometry between transcripts.
