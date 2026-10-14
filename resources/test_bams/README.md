# Bundled test and demo BAMs

Six Huh7 direct-RNA libraries (three mock-infected, `M1pA`, `M2pA`, `M3pA`; three
ZIKV-infected, `Z1pA`, `Z2pA`, `Z3pA`; 48 h post infection) subset to a few genes.
These are the same libraries, and the same processing, as the Huh7 analyses in the
Modulator manuscript.

Processing (identical for every library): dorado 2.1.0 `basecaller sup`
(`rna004_130bps_sup@v6.0.0`) with the four all-context modification models
(`inosine_m6A_2OmeA`, `2OmeG`, `m5C_2OmeC`, `pseU_2OmeU`, all v1) and
`--estimate-poly-a`, no `--modified-bases-threshold`; minimap2 2.29
`-ax splice -uf -k14 -y --secondary=no` against hg38 plus the Epstein-Barr virus
genome (NC_007605.1); samtools 1.22.1 sort. Each subset keeps every alignment
record overlapping the region (`samtools view -b -h -M <bam> <regions>`), so the
BAM headers list all 456 contigs of the reference.

| Dataset directory | Regions (hg38, 1-based) | Use |
|---|---|---|
| `MXD1_reads` | chr2:69915109-69942945 (MXD1) | `modulator demo` default |
| `ALCAM_NHSL1_SERAC1_MXD1_RIOK3_reads` | chr3:105366909-105576900 (ALCAM), chr6:138422043-138499494 (NHSL1, 3' part), chr6:158109519-158168262 (SERAC1), chr2:69915109-69942945 (MXD1), chr18:23453287-23483140 (RIOK3) | `modulator demo --dataset ALCAM_NHSL1_SERAC1_MXD1_RIOK3_reads`; default `bams_dir` in `config/config.yaml` |
| `ATG7` | chr3:11272397-11576353 (ATG7, VGLL4) | regression tests |
| `C8A` | chr1:56854797-56918223 | regression tests |
| `IP6K2` | chr3:48688003-48717221 | regression tests |
| `PCID2` | chr13:113165002-113208669 (PCID2, PROZ) | regression tests |
| `RANBP3` | chr19:5916139-5978140 | regression tests |
| `chr17_7083406_7198744_subbams` | chr17:7083406-7198744 (ASGR1, ASGR2, DLG4) | regression tests |

The published sample report (`docs/sample_report`, "demo14") is a run over the union
of all regions above (14 genes), cut from the same six libraries with the same
command, using a samplesheet with `condition` set so the between-conditions stage
runs. To reproduce it, extract the union of the regions into one BAM per library,
write a samplesheet (`sample`, `bam`, `condition`, `replicate`), and run:

```bash
modulator run --config config/config.yaml --jobs 6 \
  --set reference_fa=hg38.fa reference_gtf=hg38.ncbiRefSeq.gtf \
        samplesheet=samples.tsv bams_dir=<dir of the six BAMs> bam_glob='*.bam' \
        prefix=demo14 threads=8 genotype.enable=true genotype.jobs=6
```

The raw reads are deposited in the European Nucleotide Archive under PRJEB86065.
