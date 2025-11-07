# Modulator

A Snakemake pipeline for analyzing RNA modifications from BAM files.

## Installation

1. Clone the repository.
2. Install Snakemake and conda.
3. Configure `config/config.yaml` with your samples and references.
4. Run `snakemake` from the workflow directory.

## Usage

```bash
cd workflow
snakemake
```

## Description

This pipeline assembles transcripts from BAM files, runs modkit pileup for modifications, aggregates results per gene and transcript, and tests for differences in modification stoichiometry between transcripts.