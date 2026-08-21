#!/usr/bin/env python3

import argparse
import base64
import glob
import hashlib
import html
import itertools
import json
import os
import re
import sys
from pathlib import Path

import pandas as pd


# Display names for modkit/dorado modification codes (ChEBI ids -> friendly names).
# Dorado RNA004 v5.2.0 set: m6A, inosine, pseudoU, 5mC, and the 2'-O-methyls (Nm).
# Unknown codes are shown verbatim.
MOD_DISPLAY = {
    "a": "m6A", "m": "5mC", "h": "5hmC", "C": "4mC",
    "17596": "inosine", "17802": "pseudoU",
    "69426": "Am", "19227": "Um", "19228": "Cm", "19229": "Gm",
}


def mod_display(code):
    return MOD_DISPLAY.get(str(code), str(code))


CARD_DEFINITIONS = {
    "Fragmentforms": "Read-backed partial transcripts kept after assembly, read-support, and poly(A) filters. "
                     "One fragmentform = one distinct (intron chain + clustered 3' end) supported by at least "
                     "`min_reads` reads whose 3' ends fall within `apa_window` of each other "
                     "(see input parameters `min_reads`, `apa_window`).",
    "Genes": "Distinct genes represented across the retained fragmentforms (reference gene where the fragmentform "
             "matched the annotation, otherwise an assembled locus).",
    "Metagenes": "Connected groups of genes whose exons overlap on the same strand, merged so their overlapping "
                 "fragmentforms can be placed on separate ZN partitions. Every gene belongs to exactly one metagene, "
                 "so the number of metagenes is always ≤ the number of genes (a gene overlapping nothing is its "
                 "own metagene of one).",
    "Max ZN partitions / metagene": "The largest number of ZN partitions any single metagene needed in order to "
                                    "separate all of its overlapping genes and fragmentforms onto their own track, so "
                                    "their modification pileups are never pooled together.",
    "Assigned reads": "Total reads assigned to retained fragmentforms in the classification summary.",
    "Exact-chain reads": "Assigned reads that exactly match a retained fragmentform's canonical intron chain.",
    "Truncation-assigned reads": "Assigned reads absorbed into a retained fragmentform from suffix-compatible shorter "
                                 "chains (a read that is a 3' truncation of the fragmentform).",
    "Segregating SNPs": "Candidate SNP loci (read-backed variants passing the genotype discovery thresholds) retained "
                        "for genotype-aware association testing (see input parameters `genotype.min_alt_reads`, "
                        "`genotype.min_alt_frac`).",
    "Haplotype blocks": "Local read-backed SNP blocks retained for haplotype association testing.",
}


COLUMN_DEFINITIONS = {
    "zt_label": "Human-readable fragmentform code used in fragmentform-level outputs and BAM tags.",
    "code": "Fragmentform code or grouping identifier written by the relevant upstream stage.",
    "gtf_gene_name": "Reference gene name matched to the assembled fragmentform when available.",
    "gene_name": "Gene label assigned to the reported row or site.",
    "gene_names": "Gene labels associated with the reported SNP or molecule row.",
    "gene_ids": "Reference or assembled gene identifiers associated with the reported SNP or molecule row.",
    "gene_index": "Assembly-local gene index assigned during fragmentform summarization.",
    "transcript_index": "Within-gene fragmentform rank assigned after fragmentform sorting.",
    "metagene_index": "Index of the metagene (a connected group of same-strand overlapping genes) this row belongs to; genes that do not overlap each form their own metagene, so metagenes ≤ genes.",
    "metagene_indices": "Metagene labels associated with the reported SNP or molecule row.",
    "zn_index": "ZN partition index within the metagene — the track this fragmentform is placed on so its modification pileup is never pooled with an overlapping fragmentform.",
    "metagene_partition_count": "Number of ZN partitions this metagene needed to separate all of its overlapping genes/fragmentforms onto their own tracks.",
    "read_support": "Total reads assigned to the fragmentform or partition.",
    "exact_chain_reads": "Assigned reads whose intron chain exactly matches the retained fragmentform.",
    "trunc_assigned_reads": "Reads assigned to the fragmentform after suffix-compatible truncation absorption.",
    "anchor_reads": "Number of reads that exactly span this fragmentform's distal (5'-most) anchor point. "
                    "Support-first assignment: a fragmentform's 5' boundary is only extended to the longest intron "
                    "chain that has at least `min_distal_anchor_reads` reads demonstrably reaching it — the boundary "
                    "is driven by read support, not by the single longest read — which limits run-through across loci "
                    "(see input parameter `min_distal_anchor_reads`).",
    "anchor_frac": "Fraction of reachable suffix-family reads that exactly support the retained canonical fragmentform.",
    "assignment_mode": "Recorded suffix-family assignment policy used during fragmentform assembly.",
    "single_gene_kept": "Reads in an overlapping-gene region where only ONE gene's fragmentforms survived the "
                        "assembly/support filters, so the read is unambiguously assigned to that single gene — there "
                        "was no cross-gene conflict left to resolve.",
    "multi_gene_no_zt": "Reads spanning ≥2 overlapping genes' fragmentforms that carried NO fragmentform (ZT) tag: "
                        "the read could not be pinned to a specific fragmentform, so it could not be resolved to one "
                        "gene and was scrapped (`multi_gene_action`).",
    "zero_gene_kept": "Reads in an overlapping-gene region where NO gene's fragmentforms passed the filters, so the "
                      "read maps to zero retained fragmentforms and is scrapped.",
    # --- per-sample read funnel (per_sample_read_stats) ---
    "total_reads_bam": "Total alignment records in the input BAM for this sample.",
    "total_mapped": "Reads with a mapped primary alignment.",
    "total_unmapped": "Reads with no mapped alignment.",
    "considered_reads": "Reads that passed every primary QC filter (mapped, primary, MAPQ, intron count, 3' soft-clip) and are eligible for fragmentform assignment.",
    "failed_unmapped": "Number of reads filtered because they were unmapped.",
    "failed_secondary_or_supp": "Number of reads filtered because they were secondary or supplementary alignments.",
    "failed_low_mapq": "Number of reads filtered because their mapping quality was below the minimum (see input parameter min_mapq).",
    "failed_low_introns": "Number of reads filtered because they had too few introns (see input parameter assembler.min_introns_read).",
    "failed_low_softclip3p": "Number of reads filtered because their 3' soft-clip (candidate poly(A)) was too short (see input parameter assembler.require_softclip3p).",
    "failed_reads_total": "Total reads removed by the primary QC filters (sum of the failed_* categories).",
    "zt_tagged_exists": "Whether a ZT-tagged BAM (fragmentform assignments) was written for this sample.",
    "zt_total_records": "Total alignment records in the ZT-tagged BAM.",
    "zt_unmapped_records": "Unmapped records in the ZT-tagged BAM.",
    "zt_mapped_records": "Mapped records in the ZT-tagged BAM.",
    "zt_mapped_unassigned_reads": "Mapped reads in the ZT BAM that were not assigned to any retained fragmentform.",
    "frac_mapped_of_total": "Fraction of BAM reads that are mapped (total_mapped / total_reads_bam).",
    "frac_considered_of_total": "Fraction of BAM reads that passed QC (considered_reads / total_reads_bam).",
    "frac_considered_of_mapped": "Fraction of mapped reads that passed QC (considered_reads / total_mapped).",
    "frac_assigned_of_considered": "Fraction of QC-passing reads assigned to a fragmentform (assigned_reads / considered_reads).",
    "frac_assigned_of_total": "Fraction of BAM reads assigned to a fragmentform (assigned_reads / total_reads_bam).",
    "frac_failed_of_total": "Fraction of BAM reads removed by the QC filters.",
    # --- per-gene splice-junction summary (gene_splice_summary) ---
    "n_fragmentforms": "Number of assembled fragmentforms for this gene.",
    "n_distinct_junctions": "Number of distinct splice junctions (donor–acceptor pairs) across the gene's fragmentforms.",
    "n_canonical_GT_AG": "Junctions with canonical GT-AG dinucleotides (major U2 spliceosome).",
    "n_semi_canonical_GC_AG": "Junctions with GC-AG dinucleotides (semi-canonical).",
    "n_minor_AT_AC": "Junctions with AT-AC dinucleotides (minor U12 spliceosome).",
    "n_noncanonical": "Junctions whose dinucleotides match none of the above — worth inspecting.",
    "frac_canonical": "Fraction of the gene's junctions that are canonical GT-AG.",
    "has_noncanonical": "1 if the gene carries at least one non-canonical junction, else 0.",
    "intron_category": "Overall call for the gene: ALL_CANONICAL, CANONICAL_WITH_GC_AG, or HAS_NONCANONICAL.",
    # --- coverage-independent PRIVATE-scan table ---
    "carry_ZN": "The fragmentform (ZN) that carries the modification at this private site (highest modified fraction among the forms that contain the base).",
    "carry_arch": "3' architecture label of the carrying fragmentform (IPA / TANDEM_APA / FULL_LENGTH / REFERENCE / ...).",
    "carry_frac": "Pooled modified fraction of the site in the carrying fragmentform.",
    "carry_cov": "Pooled coverage of the site in the carrying fragmentform.",
    "n_forms_present": "Number of the gene's fragmentforms whose model contains the base (exonic).",
    "n_forms_absent": "Number of the gene's fragmentforms whose model LACKS the base (intronic/absent, not the 5' blind spot) — the private evidence.",
    "absent_in_ZN": "The fragmentform (ZN) indices that lack the base entirely (comma-separated).",
    # --- sequence-element tables ---
    "element_type": "The cis-element class scanned (PAS, ARE, CPE, GRE, G4, KOZAK, UORF, TOP, STOP_CONTEXT, M6AM).",
    "element_subclass": "Finer subclass of the element (e.g. the specific PAS hexamer or ARE variant).",
    "n_instances": "Number of instances of this element type found across all fragmentforms' mature mRNAs.",
    "n_with_modification": "How many of those instances overlap at least one modification (any mod code).",
    "frac_with_mod": "n_with_modification / n_instances.",
    "mod_codes_seen": "Modification codes found within this element type, with a count per code.",
    "n_mod_sites": "Number of distinct modification sites overlapping this element instance.",
    "matched_seq": "The mature-mRNA sequence matched for this element instance.",
    "mod_codes": "Modification codes overlapping this element instance.",
    "modifications": "Per-modification detail for the element: code@genomic-position (stoichiometry).",
    # --- poly(A) tail tables ---
    "median_tail": "Median poly(A) tail length (nt) across the fragmentform's reads.",
    "n_fragmentforms_tested": "Number of the gene's fragmentforms with enough tailed reads to compare.",
    "min_median_tail": "Smallest per-fragmentform median tail length (nt) among the gene's fragmentforms.",
    "max_median_tail": "Largest per-fragmentform median tail length (nt) among the gene's fragmentforms.",
    "effect_median_range_nt": "Spread (nt) of per-fragmentform median tail lengths within a gene — how differently its isoforms are tailed.",
    "effect_median_diff_nt": "Median tail of modified reads minus unmodified reads at a site (negative = modification associates with shorter tails).",
    "test_name": "Statistical test used: Mann-Whitney U (2 groups) or Kruskal-Wallis (>2 groups).",
    "median_tail_modified": "Median tail length (nt) of reads modified at the target site.",
    "median_tail_unmodified": "Median tail length (nt) of reads unmodified at the target site.",
    "n_unmodified": "Number of reads unmodified at the target site.",
    "sample": "Sample identifier derived from the BAM filename.",
    "chrom": "Reference chromosome or contig containing the reported feature.",
    "pos1": "1-based genomic coordinate of the reported SNP locus.",
    "mod_start0": "0-based inclusive start coordinate of the reported modification site.",
    "mod_end0": "0-based exclusive end coordinate of the reported modification site.",
    "total_reads": "Total assigned fragmentform reads for the sample.",
    "n_transcripts": "Number of fragmentforms detected in the sample.",
    "median_reads_per_tx": "Median assigned read count across fragmentforms detected in the sample.",
    "input_total_reads": "Total BAM alignments encountered before fragmentform assignment filtering.",
    "primary_reads": "Reads retained after removing secondary and supplementary alignments.",
    "mapq_reads": "Reads retained after applying the minimum MAPQ filter.",
    "intronic_reads": "Reads retained after the minimum intron-count filter.",
    "tagged_reads": "Reads written to the ZT-tagged BAM after fragmentform assignment.",
    "assigned_fraction": "Fraction of qualifying reads assigned to a retained fragmentform.",
    "assigned_reads": "Number of reads contributing to the reported fragmentform-length summary.",
    "mean_read_length": "Mean assigned read length in nucleotides.",
    "median_read_length": "Median assigned read length in nucleotides.",
    "min_read_length": "Shortest assigned read length in nucleotides.",
    "max_read_length": "Longest assigned read length in nucleotides.",
    "tes": "3' end (TES) reported for the retained fragmentform.",
    "mod_code": "Modification code reported by modkit or downstream aggregation.",
    "n_sites": "Number of unique genomic modification sites observed for the gene and modification code.",
    "p_value": "Nominal p-value from the reported hypothesis test.",
    "p_adj_bh": "Benjamini-Hochberg false-discovery-rate adjusted p-value.",
    "effect_max_abs_frac_diff": "Maximum absolute difference in pooled modified fraction across tested fragmentform partitions.",
    "effect_max_abs_tx_frac_diff": "Maximum absolute difference in fragmentform usage or stoichiometry between tested groups.",
    "effect_abs_delta_mod_frac": "Absolute difference in modified-site rate between the tested allele groups.",
    "weighted_within_tx_effect": "Coverage-weighted within-fragmentform SNP/mod effect size after fragmentform conditioning.",
    "classification": "How this fragmentform compares to the reference annotation (`reference_gtf`): EXACT (matches "
                      "an annotated transcript's chain + 3' end), NOVEL_APA (annotated chain, novel 3' end / APA site), "
                      "or NOVEL_CHAIN (intron chain not in the annotation). See input parameter `reference_gtf`.",
    "snp_id": "Canonical SNP identifier in `chrom:pos:ref>alt` format.",
    "mod_site_id": "Canonical modification-site identifier in `chrom:start-end:strand:mod` format.",
    "target_mod_code": "Modification code tested at the reported site.",
    "n_alt_reads": "Number of reads carrying the alternative allele in the tested contingency table.",
    "n_ref_reads": "Number of reads carrying the reference allele in the tested contingency table.",
    "n_reads": "Total reads contributing to the reported test.",
    "n_modified": "Reads classified as modified for the target mod code.",
    "n_not_target": "Reads classified as canonical or as another modification state for the target site.",
    "n_transcripts_tested": "Number of fragmentform partitions retained in the reported test.",
    "complete_reads": "Reads covering every SNP in the reported haplotype block.",
    "support_reads": "Reads overlapping at least one SNP in the reported haplotype block.",
    "n_snps": "Number of SNPs represented in the reported haplotype block.",
    "block_id": "Identifier for the reported haplotype block.",
    "block_region": "Genomic coordinates (chrom:start-end) spanned by the haplotype block.",
    "context_key": "Gene- or metagene-aware context key used to restrict genotype and modification joins to the same local feature family.",
    "haplotypes": "Observed allele strings for the retained read-backed haplotypes in the reported block.",
    "region": "Genomic span of the haplotype block as chrom:start-end (1-based; first to last SNP position).",
    "start1": "1-based coordinate of the first (leftmost) SNP in the haplotype block.",
    "end1": "1-based coordinate of the last (rightmost) SNP in the haplotype block.",
    "span_bp": "Distance in base pairs between the first and last SNP in the haplotype block.",
    "snp_coords": "Per-SNP coordinates and alleles in the block (chrom:pos ref>alt), in genomic order.",
    "alt_frac": "Alternative-allele fraction across all supporting reads at the candidate SNP locus.",
    "total_cov": "Total coverage accumulated across samples at the reported SNP or mod site.",
    "ref_count": "Reference-base support across all reads for the candidate SNP locus.",
    "alt_count": "Alternative-base support across all reads for the candidate SNP locus.",
    "samples_with_alt": "Samples contributing at least one alternative-allele read at the candidate SNP locus.",
    "effect_max_abs_mod_rate_diff": "Maximum absolute difference in target modified fraction across the haplotypes retained in the reported test.",
    "class_key": "Primary classification key = structural_category__stoich_direction (the structural mechanism plus which fragmentform carries more modification).",
    "structural_category": "The structural mechanism that makes the isoforms differ at the site (tandem APA, intronic polyadenylation, EJC/splicing, cassette/terminal/internal exon, …).",
    "stoich_direction": "Which fragmentform is more modified, by 3'UTR length: PROXIMAL_HIGHER (shorter form) / DISTAL_HIGHER (longer) / CO_TERMINAL (same 3' end).",
    "stoich_tier": "Magnitude of the between-isoform stoichiometry gap: T1_MARGINAL (.10–.25) / T2_MODERATE / T3_STRONG / T4_NEAR_BINARY (≥.75).",
    "hi_stoich_level": "Absolute modification level of the favored fragmentform: HI_HYPER (≥.66) / HI_INTERMED / HI_HYPO (<.33).",
    "start0": "0-based inclusive start coordinate of the reported modification site.",
    "end0": "0-based exclusive end coordinate of the reported modification site.",
    "strand": "Genomic strand of the reported feature.",
    "n_tx_tested": "Number of fragmentform (ZN) partitions retained in the differential test at the site.",
    "hi_ZN": "ZN partition index of the higher-stoichiometry isoform in the classified contrast.",
    "hi_arch": "3' architecture of the higher-stoichiometry isoform versus the gene's longest-3'UTR anchor (IPA / TANDEM_APA / FULL_LENGTH / DISTAL_EXT / REFERENCE / AMBIGUOUS).",
    "hi_frac": "Pooled modified fraction (stoichiometry) of the higher isoform at the site.",
    "lo_ZN": "ZN partition index of the lower-stoichiometry isoform in the classified contrast.",
    "lo_arch": "3' architecture of the lower-stoichiometry isoform versus the gene's longest-3'UTR anchor.",
    "lo_frac": "Pooled modified fraction (stoichiometry) of the lower isoform at the site.",
    "anchor_ZN": "ZN partition index of the anchor (longest-3'UTR) isoform used as the structural reference.",
    "status_hi": "Position status of the site within the high isoform (exonic_terminal / exonic_internal / intronic / absent).",
    "status_lo": "Position status of the site within the low isoform.",
    "status_anchor": "Position status of the site within the anchor (longest) isoform.",
    "jd_hi": "Distance (nt) from the site to the nearest spliced junction in the high isoform.",
    "jd_lo": "Distance (nt) from the site to the nearest spliced junction in the low isoform.",
}


CATEGORY_DEFINITIONS = {
    # structural_category (the primary axis): the 8 mechanisms + 3 residual labels the fine
    # 14-label `category` collapses into.
    "TANDEM_APA": "Tandem alternative polyadenylation: both isoforms share the same last exon (same acceptor) but cleave at different 3' ends; the site differs between the shorter- and longer-3'UTR forms.",
    "ALTERNATIVE_LAST_EXON": "The two isoforms are both terminal at the site but their terminal exons begin at DIFFERENT acceptor sites (distinct last exons).",
    "INTERGENIC_TERMINAL_EXON": "A non-IPA terminal exon that is genomically disjoint and far (≥ intergenic-gap) from the comparator's terminal exon.",
    "INTRONIC_POLYADENYLATION": "The high isoform terminates within an intron of the longer form (IPA); the modified base exists only in the mature IPA isoform (IPA-private, or IPA-vs-full-length with EJC relief).",
    "EJC_SPLICING": "A shared base whose modification tracks removal of a nearby splice junction / exon-junction-complex (EJC) footprint in one isoform.",
    "CASSETTE_EXON": "The base sits in an internal/cassette exon included in one isoform and spliced out (or absent) in the other.",
    "SHARED_TERMINAL_EXON": "Both isoforms share the same terminal exon AND the same 3' end; modification tracks isoform identity, not APA/EJC.",
    "SHARED_INTERNAL_EXON": "The base is in a constitutive internal exon with no junction asymmetry between the isoforms.",
    "UNEXPLAINED": "Residual: terminal in the high isoform, internal in the low, with no nearby differential junction.",
    "ARTIFACT": "The high isoform does not structurally contain the base (intronic/absent) — the 'high' stoichiometry is likely intron-read noise.",
    "IPA_UNIQUE": "High isoform is an IPA (intronic polyadenylation) isoform; the A is exonic-terminal in it but intronic/absent in the longer anchor — the modified A only exists in the mature IPA isoform. Cleavage-dependent, IPA-private.",
    "SPLICED_EXON_UNIQUE": "The A sits in an internal/cassette exon present in the high isoform but spliced out (intronic) or absent in the comparator/anchor.",
    "LAST_EXON_DISTAL_ONLY": "The A is in the anchor's distal 3'UTR but the low (proximal) isoform's cleavage site is upstream of it — a distal-3'UTR-private A.",
    "IPA_SHARED_EJC": "High isoform is IPA; the A is shared (exonic in both) but exonic-terminal in IPA versus exonic-internal in the long anchor — the A gains m6A in IPA because the downstream exon-junction complex is removed. Cleavage-dependent.",
    "SPLICING_EJC": "Shared A, non-IPA: terminalized, or a junction within the EJC window in the low/anchor is removed in the high isoform — EJC relief.",
    "LAST_EXON_PROXIMAL_APA_FAVORED": "Same terminal exon (same acceptor), different cleavage site; the PROXIMAL (shorter 3'UTR) isoform carries more m6A. Tandem 3'UTR APA, cleavage-independent.",
    "LAST_EXON_DISTAL_APA_FAVORED": "Same terminal-exon geometry; the DISTAL (longer 3'UTR) isoform carries more m6A. Tandem 3'UTR APA.",
    "ALTERNATIVE_LAST_EXON": "High and low isoforms are both terminal at the site, but their terminal exons begin at different (nearby/overlapping) acceptors — alternative last exons.",
    "INTERGENIC_TERMINAL_EXON": "The site is exonic-terminal in a non-IPA high isoform whose terminal exon is genomically disjoint from the comparator's and separated by a large gap (>= --intergenic-gap) — spatially separated / read-through / intergenic-scale alternative last exons.",
    "SHARED_TERMINAL_EXON": "High and low isoforms share the same terminal exon AND the same cleavage site; m6A tracks isoform identity, not APA or EJC.",
    "SHARED_INTERNAL_EXON": "The site is in a constitutive internal exon with no junction asymmetry — not attributable to 3' architecture.",
    "UNEXPLAINED_SHARED": "Rare residual: terminal in the high isoform, internal in the low, with no nearby differential junction.",
    "HI_INTRONIC_ARTIFACT": "The high-m6A isoform does not structurally contain the A (intronic/absent) — the 'high' stoichiometry is intron-read noise, not a real isoform-specific site.",
    "UNCLASSIFIED": "Fewer than two covered isoform models at the site, or no anchor isoform — cannot be assigned a structural category.",
}


def parse_args():
    ap = argparse.ArgumentParser(description="Generate a lightweight HTML report for modulator outputs.")
    ap.add_argument("--classification", required=True)
    ap.add_argument("--metrics", required=True)
    ap.add_argument("--tx-counts", required=True)
    ap.add_argument("--pca-png", required=True)
    ap.add_argument("--sample-stats", required=True)
    ap.add_argument("--read-stats", required=True)
    ap.add_argument("--tx-lengths", required=True)
    ap.add_argument("--partition-map", required=True)
    ap.add_argument("--out-html", required=True)
    ap.add_argument("--title", required=True)
    ap.add_argument("--run-manifest", default="",
                    help="Run manifest TXT (command line, resolved inputs, full config) to embed.")
    ap.add_argument("--sample-metadata", default="",
                    help="Sample metadata TSV (sample, condition, replicate) for the replicate-concordance note.")
    ap.add_argument("--zn-long", default="")
    ap.add_argument("--zt-long", default="")
    ap.add_argument("--diff-results", default="")
    ap.add_argument("--diff-figs-dir", default="")
    ap.add_argument("--private-sites", default="", help="{prefix}__ZN_site_private.tsv from the coverage-independent PRIVATE scan.")
    ap.add_argument("--classified-sites", default="",
                    help="{prefix}__ZN_site_classified.tsv from classify_diff_sites.py")
    ap.add_argument("--class-figs-dir", default="",
                    help="{prefix}__figs_by_category directory (one subdir per category) "
                         "holding the 2-panel per-sample stoichiometry / pooled-coverage figures")
    ap.add_argument("--arch-figs-dir", default="",
                    help="{prefix}__figs_by_category_arch directory (one subdir per category) "
                         "holding the isoform architecture-map (exon/intron locus-track) figures, "
                         "featured as the primary per-category figure")
    ap.add_argument("--max-class-figs-per-category", type=int, default=10,
                    help="max per-category figures to embed in the report. Default 10.")
    ap.add_argument("--multigene-summary-glob", default="")
    ap.add_argument("--splice-junctions", default="", help="*_splice_junctions.tsv (per-intron donor/acceptor classes)")
    ap.add_argument("--splice-genes", default="", help="*_gene_splice_summary.tsv (per-gene junction repertoire)")
    ap.add_argument("--novel-loci", default="", help="*_novel_loci.tsv (read-backed novel loci)")
    ap.add_argument("--novel-fragmentforms", default="", help="*_novel_fragmentforms.tsv")
    ap.add_argument("--mod-mod-assoc", default="", help="*_mod_mod_assoc.tsv (co-localized modification pairs)")
    ap.add_argument("--candidate-snps", default="")
    ap.add_argument("--snp-tx-assoc", default="")
    ap.add_argument("--snp-mod-assoc", default="")
    ap.add_argument("--assembled-gtf", default="",
                    help="assembled GTF (per-fragmentform exon models) used to break the cis-SNP→mod "
                         "stoichiometry hits down per fragmentform and flag which forms cover the site")
    ap.add_argument("--molecule-mod-calls", default="",
                    help="per-read modification-call table (with ZN fragmentform tag); read-backs the "
                         "per-fragmentform stoichiometry graphs in the cis-SNP→mod section")
    ap.add_argument("--molecule-snps", default="",
                    help="per-read SNP genotype table (with ZN + allele_class); splits the "
                         "per-fragmentform stoichiometry graphs by SNP allele")
    ap.add_argument("--hap-blocks", default="")
    ap.add_argument("--hap-tx-assoc", default="")
    ap.add_argument("--hap-mod-assoc", default="")
    ap.add_argument("--between-conditions-dir", default="", help="results/between_conditions dir (per-contrast differential TSVs)")
    ap.add_argument("--apa-motifs", default="", help="*_apa_motifs.tsv (PAS motif check per APA site)")
    ap.add_argument("--sequence-elements", default="", help="*_sequence_elements.tsv (cis-elements x modifications, per instance)")
    ap.add_argument("--sequence-elements-summary", default="", help="*_sequence_elements_summary.tsv (per element-type counts)")
    ap.add_argument("--snp-mod-mechanism", default="", help="*_snp_mod_mechanism.tsv (why a SNP changes a modification)")
    ap.add_argument("--polya-fragmentform", default="", help="*_polya_fragmentform.tsv (per-fragmentform tail-length distributions)")
    ap.add_argument("--taillength-diffs", default="", help="*_taillength_diffs.tsv (differential tail length between fragmentforms of a gene)")
    ap.add_argument("--taillength-mod", default="", help="*_taillength_mod.tsv (tail length vs modification state)")
    ap.add_argument("--taillength-diff-figs", default="", help="Directory of top-K per-gene tail-distribution PNGs")
    ap.add_argument("--taillength-mod-figs", default="", help="Directory of top-K per-site tail-vs-modification PNGs")
    ap.add_argument("--snp-figs-dir", default="", help="Directory to write per-example SNP/haplotype figures into (also embedded inline).")
    ap.add_argument("--max-snp-figs", type=int, default=12, help="Max per-example figures per SNP/haplotype section.")
    ap.add_argument("--max-diff-figs", type=int, default=6)
    ap.add_argument("--top-transcripts", type=int, default=20)
    ap.add_argument("--top-genes", type=int, default=20)
    return ap.parse_args()


def clean_columns(df):
    df.columns = [str(c).lstrip("#") for c in df.columns]
    return df


def read_tsv(path):
    """Tolerant TSV reader for the report. A malformed or unexpectedly-compressed input must never
    abort the whole report, so we sniff gzip magic bytes (a bgzipped file named .tsv) and skip
    (with a warning) any ragged rows rather than raising."""
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        return pd.DataFrame()
    comp = "infer"
    try:
        with open(path, "rb") as fh:
            if fh.read(2) == b"\x1f\x8b":
                comp = "gzip"
    except Exception:
        pass
    try:
        return clean_columns(pd.read_csv(path, sep="\t", compression=comp, on_bad_lines="warn"))
    except Exception as e:
        try:
            return clean_columns(pd.read_csv(path, sep="\t", compression=comp,
                                             on_bad_lines="skip", engine="python"))
        except Exception:
            sys.stderr.write(f"[report] warning: could not read {path}: {e}\n")
            return pd.DataFrame()


def read_summary_metrics(paths):
    rows = []
    for path in sorted(paths):
        sample = os.path.basename(path).split(".")[0]
        metrics = {}
        with open(path) as fh:
            for line in fh:
                line = line.rstrip("\n")
                if not line:
                    break
                if line == "metric\tvalue":
                    continue
                key, value = line.split("\t", 1)
                metrics[key] = value
        metrics["sample"] = sample
        rows.append(metrics)
    return pd.DataFrame(rows)


def embed_png(path):
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        return ""
    with open(path, "rb") as fh:
        encoded = base64.b64encode(fh.read()).decode("ascii")
    return f"data:image/png;base64,{encoded}"


def replicate_concordance(zn_long_df, meta_df, min_cov=20):
    """How reproducible are the replicates WITHIN each condition? For every site×fragmentform
    covered (Nvalid_cov >= min_cov) in >=2 replicates of a condition, take the range (max-min) of
    frac_modified across those replicates; the median of those ranges is the typical per-site
    replicate disagreement. Returns {condition: {n_reps, n_sites, median_pp}}."""
    if (zn_long_df is None or zn_long_df.empty or meta_df is None or meta_df.empty
            or "condition" not in meta_df.columns or "sample" not in meta_df.columns):
        return {}
    cond = dict(zip(meta_df["sample"].astype(str), meta_df["condition"].astype(str)))
    df = zn_long_df.copy()
    df["_cond"] = df["sample"].astype(str).map(cond)
    df = df[df["_cond"].notna()]
    df["_cov"] = pd.to_numeric(df.get("Nvalid_cov"), errors="coerce")
    df["_frac"] = pd.to_numeric(df.get("frac_modified"), errors="coerce")
    df = df[(df["_cov"] >= min_cov) & df["_frac"].notna()]
    site_keys = [k for k in ["chrom", "start0", "strand", "mod_code", "ZN_transcript_index"] if k in df.columns]
    out = {}
    for condition, cdf in df.groupby("_cond"):
        n_reps = cdf["sample"].nunique()
        entry = {"n_reps": int(n_reps), "n_sites": 0, "median_pp": None}
        if n_reps >= 2 and site_keys:
            g = cdf.groupby(site_keys)["_frac"]
            rng = (g.max() - g.min())[g.count() >= 2]
            if len(rng):
                entry["n_sites"] = int(len(rng))
                entry["median_pp"] = float(rng.median() * 100.0)
        out[condition] = entry
    return out


def significance_note_box(concordance):
    """A prominent, data-driven callout placed above the differential sections: reports THIS run's
    per-condition replicate reproducibility and warns that tight replicates make significance cheap,
    so hits should be ranked by effect size, not p-value."""
    parts = [f"within <b>{html.escape(str(c))}</b> ({d['n_reps']} replicates) the typical site "
             f"differs by only <b>~{d['median_pp']:.1f} pp</b> between replicates"
             for c, d in sorted(concordance.items()) if d.get("median_pp") is not None]
    if not parts:
        return ""
    return (
        "<div class='callout-warn'>"
        "<b>Read this first — please consider both effect size and p-value.</b> "
        "Your replicates are highly reproducible: " + "; ".join(parts) + ". "
        "When replicates agree this tightly the replicate-aware tests have very high power, so a shift "
        "of only a couple of percentage points can clear FDR even though it may be biologically less "
        "meaningful."
        "</div>"
    )


from plot_utils import save_figure, setup_matplotlib_style, bump_fonts

try:  # per-fragmentform structural coverage for the cis-SNP→mod breakdown graphs
    from classify_diff_sites import load_isoforms, status_in
except Exception:
    load_isoforms = status_in = None

_REPORT_FIGS_DIR = None  # set by main(); when set, inline summary charts also write PNG/PDF/SVG here


def _first_gene(g):
    """First gene name out of a possibly multi-gene field ('A,B' / 'A;B' / 'A|B' / 'A B')."""
    s = str(g or "")
    for sep in (",", ";", "|", "/", " "):
        s = s.replace(sep, ",")
    toks = [t for t in s.split(",") if t and t.lower() != "nan"]
    return toks[0] if toks else ""


def _parse_snp_alleles(snp_id):
    """'chr17:7115293:A>G' -> ('A', 'G'); ('', '') if unparseable."""
    try:
        ref, alt = str(snp_id).split(":")[-1].split(">")
        return ref.strip(), alt.strip()
    except Exception:
        return "", ""


def _zn_int(z):
    try:
        return int(float(z))
    except (TypeError, ValueError):
        return None


def _scan_perff_by_allele(mol_mods_path, mol_snps_path, sites):
    """For each cis-SNP→mod hit, split every fragmentform's reads at the modified site BY the SNP allele
    the read carries. Returns {(chrom,start0,mod_code): {ZN: {allele: [n_modified, n_reads]}}} with
    allele in {'ref','alt','na'} ('na' = read covers the modification but not the SNP, or no SNP table).

    Two chunked, site-filtered passes -- the per-read mod-call table (ZN + target_modified) and the
    per-read SNP table (allele_class) -- joined by (sample, read id). Bounded at genome scale because
    only the top-N hit sites are retained."""
    acc = {}
    mod_want = {(str(c), int(s), str(mc)) for (_g, c, s, mc, _snp) in sites}
    site_snp = {(str(c), int(s), str(mc)): str(snp) for (_g, c, s, mc, snp) in sites}
    snp_want = {str(snp) for (_g, _c, _s, _mc, snp) in sites if snp}
    if not mol_mods_path or not os.path.exists(mol_mods_path) or not mod_want:
        return acc
    # pass 1: reads at each mod site -> {(sample,qname): (ZN_int, modified)}
    reads_at = {k: {} for k in mod_want}
    keep = ["sample", "qname", "chrom", "start0", "target_mod_code", "target_modified", "ZN", "usable"]
    try:
        for chunk in pd.read_csv(mol_mods_path, sep="\t", usecols=lambda c: c in keep,
                                 chunksize=200000, low_memory=False):
            if not {"chrom", "start0", "target_mod_code", "ZN", "qname"}.issubset(chunk.columns):
                break
            chunk = chunk.copy()
            chunk["start0"] = pd.to_numeric(chunk["start0"], errors="coerce")
            chunk = chunk.dropna(subset=["start0"])
            chunk["start0"] = chunk["start0"].astype(int)
            # coerce target_modified up front so a NaN/blank can never raise inside the row loop
            # (a raise there would drop the whole result and silently erase every allele-split figure)
            chunk["_tmod"] = pd.to_numeric(chunk.get("target_modified", 0), errors="coerce").fillna(0).astype(int)
            if "usable" in chunk.columns:
                chunk = chunk[chunk["usable"].astype(str).str.lower().isin(("true", "1"))]
            for (c, s, mc), grp in chunk.groupby(["chrom", "start0", "target_mod_code"]):
                key = (str(c), int(s), str(mc))
                if key not in reads_at:
                    continue
                for _, r in grp.iterrows():
                    zi = _zn_int(r["ZN"])
                    if zi is None:
                        continue
                    sq = (str(r.get("sample", "")), str(r["qname"]))
                    reads_at[key][sq] = (zi, int(r["_tmod"]))
    except Exception:
        return acc
    # pass 2: allele_class per read at each wanted SNP -> {snp_id: {(sample,qname): 'ref'/'alt'}}
    allele_at = {snp: {} for snp in snp_want}
    if mol_snps_path and os.path.exists(mol_snps_path) and snp_want:
        skeep = ["sample", "qname", "snp_id", "allele_class"]
        try:
            for chunk in pd.read_csv(mol_snps_path, sep="\t", usecols=lambda c: c in skeep,
                                     chunksize=200000, low_memory=False):
                if not {"snp_id", "qname", "allele_class"}.issubset(chunk.columns):
                    break
                chunk = chunk[chunk["snp_id"].astype(str).isin(snp_want)]
                for _, r in chunk.iterrows():
                    # record the genotype for EVERY read covering the SNP (incl. third alleles), so the
                    # join can tell "not covered" (→ na) apart from "covered but off-allele" (→ dropped).
                    allele_at[str(r["snp_id"])][(str(r.get("sample", "")), str(r["qname"]))] = str(r["allele_class"]).lower()
        except Exception:
            allele_at = {snp: {} for snp in snp_want}
    # join
    for key, rmap in reads_at.items():
        amap = allele_at.get(site_snp.get(key, ""), {})
        per_zn = acc.setdefault(key, {})
        for sq, (zi, mod) in rmap.items():
            geno = amap.get(sq)                    # None = read did not cover the SNP
            if geno is None:
                allele = "na"
            elif geno in ("ref", "alt"):
                allele = geno
            else:
                continue                           # covered but a third/off allele -> not ref/alt, drop
            cell = per_zn.setdefault(zi, {}).setdefault(allele, [0.0, 0])
            cell[0] += float(mod); cell[1] += 1
    return acc


_ALLELE_STYLE = {"ref": ("#3b6ea5", "ref"), "alt": ("#c1121f", "alt"), "na": ("#17807f", "no SNP call")}


def snp_mod_fragmentform_png(gene, chrom, mod_start0, mod_code, iso, perff, ref_base="", alt_base=""):
    """Per-fragmentform modification stoichiometry at ONE cis-SNP→mod hit, each fragmentform's reads
    SPLIT by the SNP allele they carry (grouped horizontal bars: ref vs alt, plus 'no SNP call' for
    reads that miss the SNP). ``perff`` is {ZN: {allele: [n_modified, n_reads]}}. Fragmentforms with no
    reads at the site are collapsed into one summary row. base64 PNG or ''."""
    if not iso or status_in is None:
        return ""
    zns = [zn for (g, zn) in iso.keys() if g == gene]
    if not zns:
        return ""
    perff = perff or {}
    # order plotted fragmentforms by total read support; collapse the rest
    covered = []
    n_other = 0
    for zn in zns:
        d = perff.get(_zn_int(zn), {})
        tot = sum(v[1] for v in d.values())
        if tot:
            covered.append((zn, d, tot))
        else:
            n_other += 1
    covered.sort(key=lambda t: t[2], reverse=True)
    if not covered and not n_other:
        return ""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from io import BytesIO
    except Exception:
        return ""
    lab = {"ref": f"ref ({ref_base})" if ref_base else "ref",
           "alt": f"alt ({alt_base})" if alt_base else "alt", "na": "no SNP call"}
    # build a flat row list: for each fragmentform, one bar per allele present (ref, alt, na)
    bars = []  # (y, frac, color, text, ztick_label_or_None)
    yticks, ylabels = [], []
    y = 0.0
    order = ["ref", "alt", "na"]
    for zn, d, _tot in covered:
        present = [a for a in order if d.get(a, [0, 0])[1] > 0]
        y0 = y
        for a in present:
            nmod, nreads = d[a]
            frac = nmod / nreads if nreads else 0.0
            bars.append((y, frac, _ALLELE_STYLE[a][0], f"{frac*100:.0f}%  (n={nreads})", lab[a]))
            y += 1.0
        yticks.append((y0 + y - 1.0) / 2.0); ylabels.append(f"ZN{zn}")
        y += 0.6  # gap between fragmentforms
    if n_other:
        bars.append((y, 1.0, "#ededed", f"{n_other} form{'s' if n_other != 1 else ''}: no reads here", None))
        yticks.append(y); ylabels.append(f"{n_other} other")
        y += 1.0
    height = max(2.2, 0.34 * len(bars) + 0.5 * len(covered) + 1.2)
    fig, ax = plt.subplots(figsize=(8.8, height))
    for (yy, frac, color, text, _al) in bars:
        hatch = "///" if color == "#ededed" else None
        ax.barh(yy, frac if color != "#ededed" else 1.0, height=0.8, color=color,
                edgecolor="#cfcfcf" if color == "#ededed" else "white", linewidth=0.6, hatch=hatch)
        ax.text(min((frac + 0.02) if color != "#ededed" else 0.02, 0.72), yy, text,
                va="center", ha="left", fontsize=7.5, color="#4a3f35" if color != "#ededed" else "#9a9a9a")
    ax.set_yticks(yticks); ax.set_yticklabels(ylabels, fontsize=9)
    ax.invert_yaxis()  # most-supported fragmentform at the top
    ax.set_xlim(0, 1.0); ax.set_xlabel("Modified fraction (stoichiometry)")
    pos = int(mod_start0) + 1
    ax.set_title(f"{gene}  {chrom}:{pos} ({mod_code}) — per-fragmentform stoichiometry, split by SNP allele")
    # allele legend
    from matplotlib.patches import Patch
    handles = [Patch(facecolor=_ALLELE_STYLE[a][0], label=lab[a]) for a in ("ref", "alt", "na")]
    ax.legend(handles=handles, fontsize=7.5, frameon=False, loc="lower right", ncol=3)
    ax.grid(True, axis="x", linestyle="--", linewidth=0.6, alpha=0.4)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    fig.tight_layout()
    _save_report_chart(fig, f"snp_mod_fragmentform_{_first_gene(gene)}_{pos}_{mod_code}")
    buf = BytesIO()
    fig.savefig(buf, format="png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    return "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode("ascii")


def build_snp_mod_fragmentform_html(snp_mod_df, iso, mol_mods_path, mol_snps_path, top_n):
    """For the top cis-SNP→modification hits, a per-fragmentform stoichiometry breakdown at each hit's
    modified position, split by the SNP allele each read carries (read-backed; ZN = fragmentform)."""
    if snp_mod_df is None or snp_mod_df.empty or not iso or status_in is None:
        return ""
    if not {"chrom", "mod_start0", "target_mod_code"}.issubset(snp_mod_df.columns):
        return ""
    srt = snp_mod_df.sort_values("p_adj_bh", ascending=True) if "p_adj_bh" in snp_mod_df.columns else snp_mod_df
    iso_genes = {g for (g, _z) in iso.keys()}
    sites, snp_of, seen = [], {}, set()
    for _, r in srt.iterrows():
        ms0 = r.get("mod_start0")
        if pd.isna(ms0):
            continue
        chrom = str(r.get("chrom")); mc = str(r.get("target_mod_code"))
        gene = _first_gene(r.get("gene_names", ""))
        if gene not in iso_genes:
            continue  # no assembled models for this gene → cannot lay out its fragmentforms
        key = (gene, chrom, int(ms0), mc)
        if key in seen:
            continue
        seen.add(key)
        snp_id = str(r.get("snp_id", ""))
        sites.append((gene, chrom, int(ms0), mc, snp_id))
        snp_of[(gene, chrom, int(ms0), mc)] = snp_id
        if len(sites) >= top_n:
            break
    if not sites:
        return ""
    perff_all = _scan_perff_by_allele(mol_mods_path, mol_snps_path, sites)
    cards = []
    for (gene, chrom, ms0, mc, snp_id) in sites:
        perff = perff_all.get((str(chrom), int(ms0), str(mc)), {})
        ref_b, alt_b = _parse_snp_alleles(snp_id)
        img = snp_mod_fragmentform_png(gene, chrom, ms0, mc, iso, perff, ref_b, alt_b)
        if not img:
            continue
        cards.append(clickable_image_html(
            img, f"Per-fragmentform stoichiometry at {gene} {chrom}:{ms0+1} ({mc}) split by {snp_id}",
            caption=f"{gene} {chrom}:{ms0+1} ({mc}), SNP {snp_id}: within each fragmentform (ZN), the "
                    "target modification's read-backed stoichiometry is split by the SNP allele the read "
                    "carries — ref vs alt (bar = modified fraction, n = reads). Comparing ref vs alt WITHIN "
                    "a fragmentform isolates the direct allelic effect from the allele merely shifting "
                    "fragmentform usage. 'no SNP call' = reads covering the modification but not the SNP; "
                    "fragmentforms with no reads at the site are collapsed into one row."))
    if not cards:
        return ""
    return ("<h3>Per-fragmentform breakdown at top cis-SNP-linked modification sites (split by SNP allele)</h3>"
            "<p class='section-intro'>For each top hit, the target modification's stoichiometry is broken "
            "down across every fragmentform (ZN partition) of the gene AND split by SNP allele, computed "
            "directly from the per-read modification calls joined to the per-read SNP genotypes. Comparing "
            "ref vs alt within the same fragmentform is the direct allelic-effect read-out — it tells apart "
            "a variant that changes modification from one that merely shifts which fragmentform is "
            "used.</p>" + "".join(cards))


def _save_report_chart(fig, name):
    """Apply the house style (Arial-like + enlarged fonts) to an inline report summary chart, and
    persist it to disk as PNG + PDF + SVG. Called before the chart's own base64 embed on the SAME
    figure, so the inline copy inherits the identical font family + sizes."""
    setup_matplotlib_style()
    bump_fonts(fig)
    if _REPORT_FIGS_DIR:
        try:
            save_figure(fig, os.path.join(_REPORT_FIGS_DIR, name + ".png"), dpi=200, bbox_inches="tight")
        except Exception:
            pass


def category_distribution_png(counts, mod_label=None):
    """Horizontal bar chart of classified-site counts per category. Returns a
    base64 data URI (or "" if matplotlib/data unavailable). ``mod_label`` names the
    modification in the title/axis (e.g. "m6A"); None -> generic wording."""
    if not counts:
        return ""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from io import BytesIO
    except Exception:
        return ""
    items = sorted(counts.items(), key=lambda kv: kv[1])  # ascending -> largest on top
    labels = [k for k, _ in items]
    values = [v for _, v in items]
    total = sum(values) or 1
    height = max(2.4, 0.42 * len(labels) + 1.1)
    fig, ax = plt.subplots(figsize=(8.6, height))
    ypos = range(len(labels))
    bars = ax.barh(list(ypos), values, color="#c98a5e", edgecolor="#7d3c1f", linewidth=0.9)
    ax.set_yticks(list(ypos))
    ax.set_yticklabels(labels, fontsize=9)
    if mod_label:
        _ml = str(mod_label).strip()
        ax.set_xlabel(f"Significant {_ml} sites")
        ax.set_title(f"Differential {_ml} sites by structural category")
    else:
        ax.set_xlabel("Significant sites")
        ax.set_title("Differential modification sites by structural category")
    pad = max(values) * 0.01 if values else 0.1
    for rect, val in zip(bars, values):
        ax.text(rect.get_width() + pad, rect.get_y() + rect.get_height() / 2.0,
                f"{val:,} ({100.0 * val / total:.1f}%)", va="center", ha="left",
                fontsize=8, color="#4a3f35")
    ax.set_xlim(0, max(values) * 1.18 if values else 1)
    ax.grid(True, axis="x", linestyle="--", linewidth=0.6, alpha=0.4)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    fig.tight_layout()
    _save_report_chart(fig, "site_classification_distribution")
    buf = BytesIO()
    fig.savefig(buf, format="png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    encoded = base64.b64encode(buf.getvalue()).decode("ascii")
    return f"data:image/png;base64,{encoded}"


def sequence_elements_png(summ_df):
    """Stacked horizontal bar per element type: instances carrying a modification
    (dark) vs. not (light). base64 data URI, or '' if unavailable."""
    if summ_df is None or summ_df.empty or "element_type" not in summ_df.columns:
        return ""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from io import BytesIO
    except Exception:
        return ""
    df = summ_df.copy()
    df["n_instances"] = pd.to_numeric(df["n_instances"], errors="coerce").fillna(0).astype(int)
    df["n_with_modification"] = pd.to_numeric(df["n_with_modification"], errors="coerce").fillna(0).astype(int)
    df = df.sort_values("n_instances")  # largest on top
    labels = df["element_type"].tolist()
    total = df["n_instances"].tolist()
    withmod = df["n_with_modification"].tolist()
    rest = [t - w for t, w in zip(total, withmod)]
    height = max(2.4, 0.42 * len(labels) + 1.1)
    fig, ax = plt.subplots(figsize=(8.6, height))
    ypos = list(range(len(labels)))
    ax.barh(ypos, withmod, color="#c98a5e", edgecolor="#7d3c1f", linewidth=0.9, label="with a modification")
    ax.barh(ypos, rest, left=withmod, color="#ecdccf", edgecolor="#c9b7a6", linewidth=0.7, label="no modification")
    ax.set_yticks(ypos); ax.set_yticklabels(labels, fontsize=9)
    ax.set_xlabel("Element instances")
    ax.set_title("Sequence elements: instances carrying a modification (any code)")
    for y, (w, t) in enumerate(zip(withmod, total)):
        ax.text(t + max(total) * 0.01, y, f"{w:,}/{t:,}", va="center", ha="left", fontsize=8, color="#4a3f35")
    ax.set_xlim(0, max(total) * 1.20 if total else 1)
    ax.grid(True, axis="x", linestyle="--", linewidth=0.6, alpha=0.4)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.legend(loc="lower right", fontsize=8, frameon=False)
    fig.tight_layout()
    _save_report_chart(fig, "sequence_elements_distribution")
    buf = BytesIO()
    fig.savefig(buf, format="png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    return f"data:image/png;base64,{base64.b64encode(buf.getvalue()).decode('ascii')}"


def mod_mod_concordance_png(df, top_n=12):
    """Panel of 2x2 co-modification tables for the top pairs (by BH-FDR). Each panel is the
    observed 2x2 (rows = site A modified/unmodified, cols = site B modified/unmodified), coloured
    by log2(observed / expected-under-independence): red = enriched, blue = depleted. Concordant
    pairs light up the diagonal (both-modified top-left, both-unmodified bottom-right) -- i.e. they
    tend to be modified together and unmodified together. Cell text is 'obs (exp)'. base64 or ""."""
    if df is None or df.empty:
        return ""
    need = ["n_both_modified", "n_a_only", "n_b_only", "n_neither",
            "n_a_modified", "n_a_unmodified", "n_b_modified", "n_reads"]
    if any(c not in df.columns for c in need):
        return ""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
        from io import BytesIO
    except Exception:
        return ""
    sub = df.head(int(top_n))
    n = len(sub)
    if n == 0:
        return ""
    # 3 columns with generously-sized panels + constrained layout: at the enlarged house fonts the
    # per-panel titles need the extra width/height and the reflow to avoid overlapping neighbours.
    ncol = min(3, n)
    nrow = (n + ncol - 1) // ncol
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.6 * ncol, 4.2 * nrow), squeeze=False,
                             layout="constrained")
    for idx, (_, r) in enumerate(sub.iterrows()):
        ax = axes[idx // ncol][idx % ncol]
        nreads = float(r["n_reads"]) or 1.0
        a_mod = float(r["n_a_modified"]); a_un = float(r["n_a_unmodified"])
        b_mod = float(r["n_b_modified"]); b_un = nreads - b_mod
        obs = np.array([[float(r["n_both_modified"]), float(r["n_a_only"])],
                        [float(r["n_b_only"]), float(r["n_neither"])]])
        exp = np.array([[a_mod * b_mod, a_mod * b_un],
                        [a_un * b_mod, a_un * b_un]]) / nreads
        with np.errstate(divide="ignore", invalid="ignore"):
            l2 = np.log2(np.where((obs > 0) & (exp > 0), obs / np.where(exp > 0, exp, 1), 1.0))
        l2 = np.clip(np.nan_to_num(l2), -2, 2)
        ax.imshow(l2, cmap="RdBu_r", vmin=-2, vmax=2, aspect="equal")
        for i in range(2):
            for j in range(2):
                ax.text(j, i, f"{int(obs[i, j])}\n({exp[i, j]:.0f})", ha="center", va="center",
                        fontsize=8, color="#111" if abs(l2[i, j]) < 1.2 else "#fff")
        ax.set_xticks([0, 1]); ax.set_xticklabels(["B+", "B-"], fontsize=8)
        ax.set_yticks([0, 1]); ax.set_yticklabels(["A+", "A-"], fontsize=8)
        gene = str(r.get("gene_names", "")).split(";")[0]
        gene = "" if gene.strip().lower() in ("", "nan", "none") else gene[:12]
        ca, cb = str(r.get("mod_code_a", "")), str(r.get("mod_code_b", ""))
        orr = pd.to_numeric(r.get("odds_ratio", None), errors="coerce")
        orr_s = f"{orr:.1f}" if pd.notna(orr) else "n/a"
        dist = r.get("distance_bp", "")
        ttl = f"{gene} {ca}x{cb}\nd={dist}bp OR={orr_s}" if gene else f"{ca}x{cb}  d={dist}bp OR={orr_s}"
        ax.set_title(ttl, fontsize=8)
    for k in range(n, nrow * ncol):
        axes[k // ncol][k % ncol].axis("off")
    fig.suptitle("Co-localized modifications: observed vs expected 2x2 (red = enriched over independence)",
                 fontsize=10)
    _save_report_chart(fig, "mod_mod_concordance")
    buf = BytesIO()
    fig.savefig(buf, format="png", dpi=170, bbox_inches="tight")
    plt.close(fig)
    return "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode("ascii")


def clickable_image_html(src, alt, *, caption="", figure_class="image-card"):
    escaped_alt = html.escape(str(alt))
    caption_html = f"<figcaption>{html.escape(str(caption))}</figcaption>" if caption else ""
    return (
        f"<figure class='{figure_class}'>"
        f"<a class='image-link' href='{src}' target='_blank' rel='noopener noreferrer' title='Open full-size image in a new tab'>"
        f"<img src='{src}' alt='{escaped_alt}' loading='lazy' />"
        "<span class='expand-badge' aria-hidden='true'>↗</span>"
        "</a>"
        f"{caption_html}"
        "</figure>"
    )


def fmt_int(val):
    try:
        return f"{int(val):,}"
    except Exception:
        return "0"


def visible_definitions(labels):
    return [(label, CARD_DEFINITIONS[label]) for label in labels if label in CARD_DEFINITIONS]


_LEN_GROUP = {
    "total": "all alignment records",
    "considered": "reads passing the primary QC filters",
    "assigned": "reads assigned to a retained fragmentform",
    "zt_unassigned": "ZT-tagged reads not assigned to any fragmentform",
}
_LEN_STAT = {
    "min": "Minimum", "max": "Maximum", "mean": "Mean",
    "p25": "25th percentile", "p50": "Median (50th percentile)",
    "p75": "75th percentile", "p90": "90th percentile",
}


def _derived_column_def(col):
    """Definition for columns that follow a regular naming pattern (read-length percentiles, etc.)
    so they need not be enumerated one-by-one. Returns None if there is no derived definition."""
    for grp, gdesc in _LEN_GROUP.items():
        pref = f"{grp}_len_"
        if col.startswith(pref):
            stat = col[len(pref):]
            if stat in _LEN_STAT:
                return f"{_LEN_STAT[stat]} of the read-length distribution (nt) over {gdesc}."
    return None


def column_definitions(columns):
    # Only emit definitions we actually have (explicit or derived). Columns without one are left
    # out entirely rather than shown with a filler "carried through" line.
    defs = []
    seen = set()
    for col in columns:
        if col in seen:
            continue
        seen.add(col)
        desc = COLUMN_DEFINITIONS.get(col) or _derived_column_def(col)
        if desc:
            defs.append((col, desc))
    return defs


def definitions_html(items, *, summary="Definitions", open_by_default=False):
    if not items:
        return ""
    open_attr = " open" if open_by_default else ""
    defs = "".join(
        f"<dt>{html.escape(str(term))}</dt><dd>{html.escape(str(desc))}</dd>"
        for term, desc in items
    )
    return (
        f"<details class='definitions'{open_attr}>"
        f"<summary>{html.escape(summary)}</summary>"
        f"<dl>{defs}</dl>"
        "</details>"
    )


def df_to_html(df, max_rows=25):
    if df is None or df.empty:
        return "<p class='muted'>No data available.</p>"
    table = df.head(max_rows).to_html(index=False, escape=True, classes="datatable", border=0)
    return f"<div class='table-wrap'>{table}</div>"


def section(title, body, *, intro="", definitions=""):
    intro_html = f"<p class='section-intro'>{intro}</p>" if intro else ""
    # Each section is a collapsible open/close panel: the header (h2) lives in the <summary>, and
    # everything after it collapses. Defaults open.
    return (
        "<section><details class='report-section' open>"
        f"<summary><h2>{html.escape(title)}</h2></summary>"
        f"<div class='section-body'>{intro_html}{definitions}{body}</div>"
        "</details></section>"
    )


def subsection(title, body, *, definitions=""):
    return f"<div class='subsection'><h3>{html.escape(title)}</h3>{definitions}{body}</div>"


def category_figure_gallery(figs_dir, category, max_figs,
                            summary_label="site figure(s) — per-sample stoichiometry &amp; pooled coverage",
                            open_by_default=False):
    """Collapsible gallery of the per-category figures under figs_dir/<category>/*.png."""
    if not figs_dir or max_figs <= 0:
        return ""
    cat_dir = os.path.join(figs_dir, category)
    if not os.path.isdir(cat_dir):
        return ""
    fig_paths = sorted(glob.glob(os.path.join(cat_dir, "*.png")))[:max_figs]
    if not fig_paths:
        return ""
    pieces = []
    for path in fig_paths:
        img = embed_png(path)
        if img:
            pieces.append(clickable_image_html(img, os.path.basename(path), caption=os.path.basename(path)))
    if not pieces:
        return ""
    gallery = "<div class='gallery'>" + "".join(pieces) + "</div>"
    open_attr = " open" if open_by_default else ""
    return (
        f"<details class='definitions'{open_attr}>"
        f"<summary>Top {len(pieces)} {summary_label}</summary>"
        f"{gallery}</details>"
    )


def polya_distribution_png(frag_df):
    """Two-panel poly(A) overview: (left) histogram of per-fragmentform median tail length;
    (right) median tail length by fragmentform classification. base64 data URI or ""."""
    if frag_df is None or frag_df.empty or "median_tail" not in frag_df.columns:
        return ""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from io import BytesIO
    except Exception:
        return ""
    med = pd.to_numeric(frag_df["median_tail"], errors="coerce").dropna()
    if med.empty:
        return ""
    has_cls = "classification" in frag_df.columns and frag_df["classification"].notna().any()
    # constrained layout + a taller canvas so the enlarged house fonts / rotated tick labels don't
    # overlap the axis labels (build-time tight_layout can't account for the +12 bump).
    fig, axes = plt.subplots(1, 2 if has_cls else 1, figsize=(13 if has_cls else 7.0, 6.0),
                             layout="constrained")
    axes = axes if has_cls else [axes]
    axes[0].hist(med.values, bins=40, color="#3b6ea5", edgecolor="white", linewidth=0.4)
    axes[0].axvline(float(med.median()), color="#c1121f", ls="--", lw=1.2,
                    label=f"median {med.median():.0f} nt")
    axes[0].set_xlabel("fragmentform median poly(A) tail (nt)")
    axes[0].set_ylabel("fragmentforms")
    axes[0].legend(frameon=False, fontsize=8)
    if has_cls:
        order = [c for c in ["EXACT", "NOVEL_APA", "NOVEL_CHAIN", "NOVEL_LOCUS"]
                 if c in set(frag_df["classification"])]
        order += [c for c in frag_df["classification"].dropna().unique() if c not in order]
        data = [pd.to_numeric(frag_df.loc[frag_df["classification"] == c, "median_tail"],
                              errors="coerce").dropna().values for c in order]
        axes[1].boxplot(data, tick_labels=order, showfliers=False)  # 'labels=' removed in mpl 3.11
        axes[1].set_ylabel("median tail (nt)")  # short label so the +12 font bump doesn't clip it
        axes[1].set_xlabel("fragmentform class")
        for t in axes[1].get_xticklabels():
            t.set_rotation(25); t.set_ha("right")
    _save_report_chart(fig, "polya_tail_distribution")
    buf = BytesIO()
    fig.savefig(buf, format="png", dpi=130, bbox_inches="tight")
    plt.close(fig)
    import base64 as _b64
    return "data:image/png;base64," + _b64.b64encode(buf.getvalue()).decode("ascii")


def build_between_conditions_section(bc_dir, top_n, top_note=""):
    """Replicate-aware between-condition results, one subsection per contrast."""
    if not bc_dir or not os.path.isdir(bc_dir):
        return ""
    kinds = [("mod_diffs", "Differential modification", ["gene_name", "chrom", "start0", "strand", "mod_code",
                                                         "mu_reference", "mu_test", "delta", "p_adj_bh"], "delta"),
             ("isoform_usage_diffs", "Differential isoform usage", ["gene_name", "feature", "mu_reference",
                                                                    "mu_test", "delta", "p_adj_bh"], "delta"),
             ("apa_usage_diffs", "Differential APA-site usage", ["gene_name", "feature", "mu_reference",
                                                                 "mu_test", "delta", "p_adj_bh"], "delta"),
             ("junction_usage_diffs", "Differential junction usage", ["gene_name", "feature", "mu_reference",
                                                                      "mu_test", "delta", "p_adj_bh"], "delta"),
             ("tail_diffs", "Differential poly(A) tail length", ["gene_name", "feature", "median_tail_reference",
                                                                 "median_tail_test", "delta_nt", "p_adj_bh"], "delta_nt")]
    found = {}
    for path in sorted(glob.glob(os.path.join(bc_dir, "*.tsv"))):
        base = os.path.basename(path)
        for suffix, _, _, _ in kinds:
            if base.endswith(f"_{suffix}.tsv"):
                df = read_tsv(path)
                if df.empty or "contrast" not in df.columns:
                    continue
                found.setdefault(str(df["contrast"].iloc[0]), {})[suffix] = df
                break
    if not found:
        return ""
    parts = []
    for contrast, by_kind in found.items():
        summary, blocks = [], []
        for suffix, title, cols, eff in kinds:
            df = by_kind.get(suffix)
            if df is None or df.empty:
                continue
            padj = pd.to_numeric(df.get("p_adj_bh", pd.Series(dtype=float)), errors="coerce")
            e = pd.to_numeric(df.get(eff, pd.Series(dtype=float)), errors="coerce").abs()
            # Effect thresholds differ by units, so name the column consistently and report the
            # threshold in its own column (otherwise each analysis makes its own NaN-filled column).
            thr = 10.0 if eff == "delta_nt" else 0.10     # 10 nt of tail, or 10 percentage points
            summary.append({"analysis": title, "tested": len(df),
                            "FDR<0.05": int((padj < 0.05).sum()),
                            "effect_threshold": f"|{eff}| >= {thr:g}" + (" nt" if eff == "delta_nt" else ""),
                            "FDR<0.05 & above threshold": int(((padj < 0.05) & (e >= thr)).sum())})
            sig = df[(padj < 0.05)].copy()
            if not sig.empty:
                sig["_e"] = pd.to_numeric(sig[eff], errors="coerce").abs()
                sig = sig.sort_values("_e", ascending=False).drop(columns="_e")
                blocks.append(subsection(f"{title} — top by effect size",
                                         df_to_html(sig[[c for c in cols if c in sig.columns]], max_rows=top_n)))
        if summary:
            parts.append(subsection(f"Contrast: {contrast}",
                                    df_to_html(pd.DataFrame(summary), max_rows=10) + "".join(blocks)))
    return section(
        "Between-Condition Changes in Fragment-form Usage, Junction Usage, APA Usage, polyA Tail Length, and Modification Stoichiometries",
        top_note + "".join(parts),
        intro="Replicate-aware comparisons between conditions, from the samplesheet's `condition` column. "
              "Counts (modification and isoform/APA/junction usage) use a beta-binomial likelihood-ratio "
              "test with dispersion shrinkage across features; poly(A) tail length is continuous and is "
              "compared with Welch's t-test across per-replicate medians. Reads are never pooled across "
              "replicates -- the biological unit is the replicate, and pooling would make trivial "
              "differences look overwhelmingly significant.",
        definitions=definitions_html([
            ("mu_reference / mu_test", "Fitted modified (or usage) fraction in the reference and test condition."),
            ("delta", "mu_test - mu_reference: the effect size, in fraction units (delta_nt = nucleotides of tail)."),
            ("p_adj_bh", "Benjamini-Hochberg FDR over all features in that analysis."),
        ], summary="Column definitions"),
    )


def build_apa_motif_section(apa_df, top_n):
    """PAS motif support for every APA site, and the internal-priming artifact flag."""
    if apa_df is None or apa_df.empty or "apa_motif_class" not in apa_df.columns:
        return section("APA Motifs (Polyadenylation Signals) in Fragmentforms",
                       "<p class='muted'>No APA motif results available.</p>",
                       intro="Polyadenylation-signal support for each detected APA site.")
    counts = apa_df["apa_motif_class"].value_counts().to_dict()
    total = int(sum(counts.values())) or 1
    parts = []
    fig = category_distribution_png(counts)
    if fig:
        parts.append(clickable_image_html(fig, "APA PAS class distribution",
                                          caption="Polyadenylation-signal class per APA site."))
    n_pas = counts.get("PAS_CANONICAL", 0) + counts.get("PAS_VARIANT", 0)
    n_ip = counts.get("PAS_NONE_INTERNAL_PRIMING", 0)
    parts.append("<h3>Summary Across All Measured APA Sites</h3>")
    parts.append(
        "<ul>"
        f"<li><b>{total:,}</b> APA sites checked — <b>{n_pas:,} ({100.0 * n_pas / total:.1f}%)</b> carry a polyadenylation signal</li>"
        f"<li><b>{counts.get('PAS_CANONICAL', 0):,}</b> canonical <code>AATAAA</code>; "
        f"<b>{counts.get('PAS_VARIANT', 0):,}</b> a variant hexamer</li>"
        f"<li><b>{n_ip:,}</b> flagged <b>likely internal priming</b> (no PAS + A-rich downstream genome)</li>"
        "</ul>"
    )
    summ = pd.DataFrame({"apa_motif_class": list(counts.keys()), "n_sites": list(counts.values())})
    summ["pct"] = (100.0 * summ["n_sites"] / total).round(2)
    parts.append(subsection("PAS class summary", df_to_html(summ.sort_values("n_sites", ascending=False), max_rows=10)))
    hx = (apa_df[apa_df["pas_motif"].astype(str).ne("") & apa_df["pas_motif"].notna()]
          if "pas_motif" in apa_df.columns else apa_df.iloc[0:0])
    if not hx.empty and "pas_distance_nt" in hx.columns:
        hu = hx["pas_motif"].value_counts().reset_index()
        hu.columns = ["pas_motif", "n_sites"]
        med = hx.groupby("pas_motif")["pas_distance_nt"].median()
        hu["median_distance_upstream_nt"] = hu["pas_motif"].map(med)
        parts.append(subsection("PAS hexamer usage (and distance upstream of the cleavage site)",
                                df_to_html(hu, max_rows=12)))
    ip = apa_df[apa_df["apa_motif_class"].eq("PAS_NONE_INTERNAL_PRIMING")]
    if not ip.empty:
        cols = [c for c in ["gene_name", "zt_label", "chrom", "strand", "tes", "downstream_a_frac",
                            "fragmentform_class", "read_support"] if c in ip.columns]
        parts.append(subsection("Sites flagged as likely internal priming", df_to_html(ip[cols], max_rows=top_n)))
    return section(
        "APA Motifs (Polyadenylation Signals) in Fragmentforms",
        "".join(parts),
        intro="Each fragmentform's TES is a cleavage/polyadenylation site. For every one, the genomic "
              "sequence around it is read in sense (fragmentform) orientation and scanned for a polyadenylation "
              "signal (canonical AATAAA or a known variant hexamer, normally 10-30 nt upstream) plus a "
              "downstream U/GU-rich element. This both annotates the site and filters artifacts: a site "
              "with no PAS whose downstream genome is A-rich was most likely produced by an oligo-dT "
              "primer annealing to an internal A-stretch (internal priming), not by real cleavage.",
        definitions=definitions_html([
            ("PAS_CANONICAL", "Canonical AATAAA hexamer found upstream — the textbook signal."),
            ("PAS_VARIANT", "One of the 11 known variant hexamers (ATTAAA, TATAAA, AGTAAA, ...)."),
            ("PAS_NONE_INTERNAL_PRIMING", "No PAS AND A-rich downstream genome — likely an oligo-dT internal-priming artifact, not a real APA site."),
            ("PAS_NONE", "No PAS but not A-rich — a non-canonical or novel site worth inspecting."),
            ("downstream_a_frac", "A fraction of the genomic sequence just downstream of the cleavage site; high values are the internal-priming signature."),
        ], summary="Category definitions"),
    )


def _elem_mods_string(mods_json):
    """Compact 'code@pos frac' summary from a sequence_elements mods_json cell."""
    try:
        ms = json.loads(mods_json) if isinstance(mods_json, str) and mods_json.strip() else []
    except (ValueError, TypeError):
        return ""
    return "; ".join(f"{m.get('code')}@{m.get('gpos')} ({m.get('frac')})" for m in ms)


def build_sequence_elements_section(se_df, summ_df, top_n):
    """Sequence-based cis-elements on each fragmentform's mature mRNA, and every
    overlapping modification -- unbiased across mod codes (m6A, 5mC, pseudouridine, …)."""
    intro = ("Every fragmentform's mature (spliced) mRNA is scanned for sequence-defined cis-elements "
             "anchored to its 3' end (PAS, AU-rich element, CPE, GU-rich element, rG4), start codon "
             "(Kozak, uORF, 5'TOP, m6Am) and stop codon (readthrough context). Each element is then "
             "joined to EVERY overlapping modification with no mod-code bias. Only mature-mRNA (exonic) "
             "elements are scanned; 5'-anchored elements (Kozak/uORF/TOP/m6Am) sit where direct-RNA read "
             "coverage is sparsest, so their modification calls are coverage-limited.")
    if se_df is None or se_df.empty or "element_type" not in se_df.columns:
        return section("Modifications within RNA Sequence Elements",
                       "<p class='muted'>No sequence-element results available.</p>", intro=intro)
    parts = []
    n_inst = len(se_df)
    se_df = se_df.copy()
    se_df["n_mod_sites"] = pd.to_numeric(se_df.get("n_mod_sites", 0), errors="coerce").fillna(0).astype(int)
    with_mod = se_df[se_df["n_mod_sites"] > 0]
    codes = sorted({c for cell in with_mod.get("mod_codes", pd.Series(dtype=str)).dropna()
                    for c in str(cell).split(",") if c})
    fig = sequence_elements_png(summ_df)
    if fig:
        parts.append(clickable_image_html(fig, "Sequence elements: modification-carrying instances per type",
                                          caption="Per element type, how many instances overlap a modification (any code)."))
    parts.append("<h3>Summary Across All Sequence Elements</h3>")
    parts.append(
        "<ul>"
        f"<li><b>{n_inst:,}</b> element instances across <b>{se_df['element_type'].nunique()}</b> element types</li>"
        f"<li><b>{len(with_mod):,}</b> carry ≥1 modification — mod codes seen: "
        f"<code>{', '.join(codes) if codes else '—'}</code> (no code is filtered out)</li>"
        "</ul>"
    )
    # per element-type summary (prefer the precomputed summary table)
    if summ_df is not None and not summ_df.empty:
        parts.append(subsection("Elements by type (and how many carry a modification)",
                                df_to_html(summ_df, max_rows=20),
                                definitions=definitions_html(column_definitions(list(summ_df.columns)), summary="Column definitions")))
    # per-gene rollup: how many of each gene's fragmentforms have a MODIFIED element of each type
    if not with_mod.empty and "gene_name" in with_mod.columns:
        tx_col = "zt_label" if "zt_label" in with_mod.columns else None
        if tx_col:
            piv = (with_mod.groupby(["gene_name", "element_type"])[tx_col]
                   .nunique().unstack(fill_value=0))
        else:
            piv = with_mod.groupby(["gene_name", "element_type"]).size().unstack(fill_value=0)
        piv = piv.reset_index()
        parts.append(subsection(
            "Per-gene rollup — fragmentforms with a MODIFIED element of each type",
            df_to_html(piv, max_rows=max(top_n, 25)),
            definitions=definitions_html([
                ("gene_name", "The gene."),
                ("(each element-type column)", "Number of that gene's fragmentforms whose mature mRNA carries a MODIFIED element of that type. Columns are element types (defined under Element definitions on the section header)."),
            ], summary="Column definitions")))
    # the modification-carrying elements themselves (deduped across repeated fragmentforms)
    if not with_mod.empty:
        w = with_mod.copy()
        w["modifications"] = w["mods_json"].map(_elem_mods_string) if "mods_json" in w.columns else ""
        cols = [c for c in ["gene_name", "element_type", "element_subclass", "chrom", "strand",
                            "matched_seq", "mod_codes", "modifications"] if c in w.columns]
        w = w[cols].drop_duplicates()
        parts.append(subsection("Elements carrying a modification (any code)",
                                df_to_html(w, max_rows=max(top_n, 25)),
                                definitions=definitions_html(column_definitions(cols), summary="Column definitions")))
    return section("Modifications within RNA Sequence Elements", "".join(parts), intro=intro,
                   definitions=definitions_html([
                       ("PAS", "Polyadenylation signal hexamer (AATAAA + variants) upstream of the 3' end."),
                       ("ARE", "AU-rich element (AUUUA), a 3'UTR stability determinant."),
                       ("CPE", "Cytoplasmic polyadenylation element (UUUUAU)."),
                       ("GRE", "GU-rich element (GU repeats)."),
                       ("G4", "RNA G-quadruplex sequence motif (G≥3 runs)."),
                       ("KOZAK", "Translation-initiation context around the start codon (−3 purine / +4 G)."),
                       ("UORF", "Upstream ORF (AUG…stop) in the 5'UTR."),
                       ("TOP", "5' terminal oligopyrimidine tract."),
                       ("STOP_CONTEXT", "Stop codon + downstream base (readthrough-prone UGA-C, …)."),
                       ("M6AM", "Cap-adjacent first transcribed A (potential m6Am)."),
                   ], summary="Element definitions"))


def build_snp_mechanism_section(mech_df, top_n):
    """Why a SNP changes a modification: positional ladder, m6A motif effect, direction concordance."""
    if mech_df is None or mech_df.empty or "positional_class" not in mech_df.columns:
        return section("Distance of cis SNP to Affected Modifications",
                       "<p class='muted'>No SNP-to-modification mechanism results available.</p>",
                       intro="Positional and motif-level explanation of each SNP x modification association.")
    parts = []
    counts = mech_df["positional_class"].value_counts().to_dict()
    fig = category_distribution_png(counts)
    if fig:
        parts.append(clickable_image_html(fig, "SNP positional class distribution",
                                          caption="Where each SNP sits relative to the modified base."))
    padj = pd.to_numeric(mech_df.get("p_adj_bh"), errors="coerce")
    sig = mech_df[padj < 0.05]
    n_art_sig = int(sig["artifact_flag"].ne("NONE").sum()) if "artifact_flag" in sig.columns else 0
    clean = mech_df[mech_df["artifact_flag"].eq("NONE")] if "artifact_flag" in mech_df.columns else mech_df
    conc = clean[clean["direction_concordance"].isin(["CONCORDANT", "DISCORDANT"])]
    rate = (100.0 * conc["direction_concordance"].eq("CONCORDANT").mean()) if len(conc) else float("nan")
    parts.append(
        "<ul>"
        f"<li><b>{len(mech_df):,}</b> SNP x modification pairs classified</li>"
        f"<li><b>{n_art_sig:,}</b> of the <b>{len(sig):,}</b> significant (FDR&lt;0.05) pairs are "
        f"<b>self-reporting artifacts</b> — the variant and the modification are potentially the same physical event</li>"
        + (f"<li>Motif prediction matches the observed allelic direction in <b>{rate:.1f}%</b> of "
           f"testable motif-bearing pairs (n={len(conc)}) — the mechanism is internally consistent</li>" if len(conc) else "")
        + "</ul>"
    )
    if "artifact_flag" in mech_df.columns:
        af = mech_df["artifact_flag"].value_counts().reset_index()
        af.columns = ["artifact_flag", "n_pairs"]
        parts.append(subsection("Self-reporting / definitional artifacts", df_to_html(af, max_rows=8)))
    me = mech_df[mech_df["motif_effect"].ne("NOT_APPLICABLE")] if "motif_effect" in mech_df.columns else pd.DataFrame()
    if not me.empty:
        mt = me["motif_effect"].value_counts().reset_index()
        mt.columns = ["motif_effect", "n_pairs"]
        parts.append(subsection("Motif effect (SNPs inside the 5-mer)", df_to_html(mt, max_rows=8)))
    causal = clean[clean["motif_effect"].eq("MOTIF_DISRUPTED") & clean["direction_concordance"].eq("CONCORDANT")] \
        if "motif_effect" in clean.columns else pd.DataFrame()
    if not causal.empty:
        cols = [c for c in ["snp_id", "gene_names", "mod_site_id", "distance_bp", "ref_5mer", "alt_5mer",
                            "ref_mod_rate", "alt_mod_rate", "p_adj_bh"] if c in causal.columns]
        causal = causal.assign(_p=pd.to_numeric(causal["p_adj_bh"], errors="coerce")).sort_values("_p")
        parts.append(subsection("Top causal cis variants (motif disrupted, direction concordant)",
                                df_to_html(causal[cols], max_rows=top_n)))
    return section(
        "Distance of cis SNP to Affected Modifications",
        "".join(parts),
        intro="A SNP whose alleles carry different modification rates gets explained on three axes: WHERE "
              "it sits relative to the modified base (at it / inside the 5-mer / inside the 9-mer / "
              "proximal / distal cis), whether the alt allele DISRUPTS or CREATES the modification's "
              "sequence motif, and whether the observed allelic direction matches what the motif predicts. "
              "A sequence motif is used only where the modification has one (the DRACH consensus for m6A); "
              "modifications without a comparable motif (e.g. pseudouridine, inosine) get positional context "
              "only.",
        definitions=definitions_html([
            ("AT_MOD_BASE", "The SNP is the modified base itself."),
            ("IN_MOTIF_CORE / IN_MOTIF_EXTENDED", "Inside the 5-mer (<=2 nt) / the 9-mer (<=4 nt) around the modified base."),
            ("PROXIMAL_CIS / DISTAL_CIS", "Near (default <=50 nt) vs far from the modified base, same locus."),
            ("MOTIF_DISRUPTED / MOTIF_CREATED", "The alt allele breaks / creates the modification's sequence motif — predicting less / more modification on alt."),
            ("MOTIF_ABSENT_BOTH", "Neither allele matches the modification's motif, which makes the modification call itself suspect."),
            ("EDITING_SELF_REPORT / PSEU_SELF_REPORT", "The 'SNP' may itself BE the modification: A-to-I editing basecalls as G, and pseudouridine causes a U-to-C basecall error, so each can get called as a variant at its own site. The association is then circular rather than regulatory."),
            ("MOD_BASE_ABLATED", "SNP at a modified base whose alt allele cannot carry the modification (e.g. an m6A site whose alt allele is not an A) — there is no substrate to modify, so the association is definitional."),
            ("direction_concordance", "Does the observed allelic direction match the motif's prediction? CONCORDANT = coherent causal cis variant."),
        ], summary="Category definitions"),
    )


def _flat_figure_gallery(figs_dir, max_figs, summary_label):
    """Collapsible gallery of figs_dir/*.png (flat, rank-sorted). For the poly(A) top-K figures."""
    if not figs_dir or max_figs <= 0 or not os.path.isdir(figs_dir):
        return ""
    fig_paths = sorted(glob.glob(os.path.join(figs_dir, "*.png")))[:max_figs]
    pieces = []
    for path in fig_paths:
        img = embed_png(path)
        if img:
            pieces.append(clickable_image_html(img, os.path.basename(path), caption=os.path.basename(path)))
    if not pieces:
        return ""
    return (f"<details class='definitions' open><summary>Top {len(pieces)} {summary_label}</summary>"
            f"<div class='gallery'>{''.join(pieces)}</div></details>")


def build_polya_section(frag_df, diffs_df, mod_df, top_n, diff_figs_dir="", mod_figs_dir="", max_figs=10, top_note=""):
    """Poly(A) tail length as a first-class readout: distribution across fragmentforms,
    differential tail length between the fragmentforms of a gene, and tail x modification."""
    if (frag_df is None or frag_df.empty) and (diffs_df is None or diffs_df.empty) and (mod_df is None or mod_df.empty):
        return section("Differential polyA Tail Length by Fragmentform and Modification Status",
                       "<p class='muted'>No poly(A) tail-length data. Requires reads basecalled with "
                       "dorado <code>--estimate-poly-a</code> (the <code>pt:i</code> tag).</p>",
                       intro="Per-read dorado poly(A) tail-length estimates, grouped by fragmentform.")
    parts = []
    fig = polya_distribution_png(frag_df)
    if fig:
        parts.append(clickable_image_html(fig, "poly(A) tail-length distribution",
                                          caption="Left: per-fragmentform median tail length. "
                                                  "Right: median tail by fragmentform class."))
    if frag_df is not None and not frag_df.empty:
        med_all = pd.to_numeric(frag_df["median_tail"], errors="coerce").dropna()
        n_sig_diff = int((pd.to_numeric(diffs_df.get("p_adj_bh"), errors="coerce") < 0.05).sum()) if diffs_df is not None and not diffs_df.empty else 0
        n_sig_mod = int((pd.to_numeric(mod_df.get("p_adj_bh"), errors="coerce") < 0.05).sum()) if mod_df is not None and not mod_df.empty else 0
        parts.append("<h3>Summary Across All Poly(A) Tail Measurements</h3>")
        _med = med_all.median()
        _med_s = f"{_med:.0f} nt" if pd.notna(_med) else "n/a"
        parts.append(
            "<ul>"
            f"<li><b>{len(frag_df):,}</b> fragmentforms with a tail-length distribution; "
            f"overall median <b>{_med_s}</b></li>"
            f"<li><b>{n_sig_diff:,}</b> genes with differential tail length between their fragmentforms (FDR&lt;0.05)</li>"
            f"<li><b>{n_sig_mod:,}</b> modification sites where tail length differs by modification state (FDR&lt;0.05)</li>"
            "</ul>"
        )
    if diffs_df is not None and not diffs_df.empty:
        cols = [c for c in ["gene_name", "n_fragmentforms_tested", "n_reads", "effect_median_range_nt",
                            "min_median_tail", "max_median_tail", "test_name", "p_value", "p_adj_bh"] if c in diffs_df.columns]
        gallery = _flat_figure_gallery(diff_figs_dir, max_figs, "per-gene tail-distribution figure(s) — tail length by fragmentform")
        parts.append(subsection("Differential tail length between fragmentforms of a gene",
                                df_to_html(diffs_df[cols], max_rows=top_n) + gallery,
                                definitions=definitions_html(column_definitions(cols), summary="Column definitions")))
    if mod_df is not None and not mod_df.empty:
        cols = [c for c in ["gene_name", "mod_site_id", "target_mod_code", "n_modified", "n_unmodified",
                            "median_tail_modified", "median_tail_unmodified", "effect_median_diff_nt",
                            "p_value", "p_adj_bh"] if c in mod_df.columns]
        gallery = _flat_figure_gallery(mod_figs_dir, max_figs, "per-site figure(s) — modified vs unmodified tail length")
        parts.append(subsection("Poly(A) tail length vs modification state",
                                df_to_html(mod_df[cols], max_rows=top_n) + gallery,
                                definitions=definitions_html(column_definitions(cols), summary="Column definitions")))
    return section(
        "Differential polyA Tail Length by Fragmentform and Modification Status",
        top_note + "".join(parts),
        intro="Direct per-read poly(A) tail-length estimates from the dorado basecaller (pt:i tag), "
              "grouped by fragmentform (isoform). Tail length is mechanistically tied to m6A, RNA "
              "stability, and translation, so it is reported as a first-class readout: its distribution "
              "across fragmentforms, whether the fragmentforms of a gene differ in tail length, and "
              "whether modification at a site tracks with a read's tail length.",
        definitions=definitions_html([
            ("effect_median_range_nt", "Spread (nt) of per-fragmentform median tail lengths within a gene — how differently its isoforms are tailed."),
            ("effect_median_diff_nt", "Median tail length of modified reads minus unmodified reads at a site (negative = modification associates with shorter tails)."),
            ("test", "Mann-Whitney U (2 groups) or Kruskal-Wallis (>2) on tail-length distributions; p_adj_bh is BH-FDR."),
        ], summary="Column definitions"),
    )


CLASS_TAXONOMY = {
    "PRIVATE":       ["SKIPPED_EXON", "INTRONIC_POLYA", "ALT_LAST_EXON"],
    "SHARED_LOCAL":  ["ALT_DONOR", "ALT_ACCEPTOR", "ALT_POLYA_SITE", "RETAINED_INTRON",
                      "IPA_EXTENSION", "NEAR_ALT_JUNCTION"],
    "SHARED_DISTAL": ["DISTAL_APA", "DISTAL_SPLICING"],
    "UNEXPLAINABLE": ["FIVE_PRIME_UNCERTAIN", "INTRON_READ_ARTIFACT", "NO_MODEL", "UNRESOLVED"],
}

# Level 3 of the tree: the direction sub-categories each event can resolve to. Every direction listed
# here is always shown for its event (with a zero count when none landed there), so the tree is a
# complete, self-documenting enumeration rather than only whatever the run happened to produce.
EVENT_DIRECTIONS = {
    "SKIPPED_EXON":         ["WITH_EXON_HIGHER"],
    "INTRONIC_POLYA":       ["IPA_TRANSCRIPT_HIGHER"],
    "ALT_LAST_EXON":        ["LONGER_EXON_HIGHER", "SHORTER_EXON_HIGHER"],
    "ALT_DONOR":            ["LONGER_EXON_HIGHER", "SHORTER_EXON_HIGHER"],
    "ALT_ACCEPTOR":         ["LONGER_EXON_HIGHER", "SHORTER_EXON_HIGHER"],
    "ALT_POLYA_SITE":       ["PROXIMAL_HIGHER", "DISTAL_HIGHER"],
    "RETAINED_INTRON":      ["INTRON_RETAINED_HIGHER", "INTRON_RETAINED_LOWER"],
    "IPA_EXTENSION":        ["EXTENDS_TO_PA_HIGHER", "EXTENDS_TO_PA_LOWER"],
    "NEAR_ALT_JUNCTION":    ["JUNCTION_REMOVED_HIGHER", "JUNCTION_PRESENT_HIGHER"],
    "DISTAL_APA":           ["PROXIMAL_HIGHER", "DISTAL_HIGHER"],
    "DISTAL_SPLICING":      ["CO_TERMINAL_HIGHER"],
    "FIVE_PRIME_UNCERTAIN": [],
    "INTRON_READ_ARTIFACT": [],
    "NO_MODEL":             [],
    "UNRESOLVED":           [],
}
DIRECTION_DEFINITIONS = {
    "WITH_EXON_HIGHER":        "The more-modified fragmentform is the one that splices the cassette exon IN.",
    "IPA_TRANSCRIPT_HIGHER":   "The more-modified fragmentform is the intronic-polyadenylation isoform (base lives in its retained-intron 3'UTR).",
    "LONGER_EXON_HIGHER":      "The more-modified fragmentform carries the LONGER version of the base's exon (structural_delta_nt = the extra length).",
    "SHORTER_EXON_HIGHER":     "The more-modified fragmentform carries the SHORTER version of the base's exon.",
    "PROXIMAL_HIGHER":         "The more-modified fragmentform uses the PROXIMAL (upstream, shorter-3'UTR) poly(A) site.",
    "DISTAL_HIGHER":           "The more-modified fragmentform uses the DISTAL (downstream, longer-3'UTR) poly(A) site.",
    "INTRON_RETAINED_HIGHER":  "The more-modified fragmentform is the one that RETAINS the intron (its exon spans it); the other splices it out.",
    "INTRON_RETAINED_LOWER":   "The LESS-modified fragmentform is the one that retains the intron; the more-modified one splices it out.",
    "EXTENDS_TO_PA_HIGHER":    "The more-modified fragmentform is the one that reads into the intron and polyadenylates there (its exon is the extended, intronic-polyadenylation exon).",
    "EXTENDS_TO_PA_LOWER":     "The LESS-modified fragmentform is the one that reads into the intron and polyadenylates there.",
    "JUNCTION_REMOVED_HIGHER": "The more-modified fragmentform is the one in which the nearby splice junction has been REMOVED (EJC relief).",
    "JUNCTION_PRESENT_HIGHER": "The more-modified fragmentform is the one that still has the nearby splice junction.",
    "CO_TERMINAL_HIGHER":      "The two fragmentforms share the same 3' end; the more-modified one differs only in an internal / 5' splicing choice elsewhere.",
}
CLASS_BUCKET_ORDER = ["PRIVATE", "SHARED_LOCAL", "SHARED_DISTAL", "UNEXPLAINABLE"]
BUCKET_DEFINITIONS = {
    "PRIVATE": "The modified base exists (exonically) ONLY in the more-modified fragmentform — it is intronic or absent in the other form's assembled model, so the difference is structurally trivial: the base is physically missing from one isoform. (Called from transcript-model coordinates, guarded against the 5' blind spot.)",
    "SHARED_LOCAL": "The base is exonic in BOTH fragmentforms, and the structural difference between them lies ON / adjacent to the base's own exon — a shifted splice site (ALT_DONOR/ALT_ACCEPTOR), an alternative poly(A) site, a retained intron, an intronic-polyadenylation extension, or a nearby differential junction. This holds even when the base's exon is the terminal exon in one form.",
    "SHARED_DISTAL": "The base is exonic in both, in an IDENTICAL local context; the fragmentforms differ only ELSEWHERE in the transcript, so the modification difference tracks isoform identity rather than any feature at the base.",
    "UNEXPLAINABLE": "The difference cannot be attributed to structure: the 5' blind spot (direct-RNA truncation), an intron / soft-clip read artifact, or too few covered isoforms to compare.",
}
EVENT_DEFINITIONS = {
    "SKIPPED_EXON": "The base sits in a cassette exon spliced INTO the higher form and OUT of the lower form.",
    "INTRONIC_POLYA": "The base is in the retained-intron 3'UTR of an intronic-polyadenylation isoform — it does not exist in the spliced full-length form.",
    "ALT_LAST_EXON": "The base is in a mutually-exclusive alternative last exon used by the higher form.",
    "ALT_DONOR": "The base's exon uses a different 5' splice site (donor) in the two forms, changing that exon's length (structural_delta_nt).",
    "ALT_ACCEPTOR": "The base's exon uses a different 3' splice site (acceptor) in the two forms, changing that exon's length (structural_delta_nt).",
    "ALT_POLYA_SITE": "Same last exon, different poly(A) cleavage sites (tandem APA); the base is in the shared body of the last exon.",
    "RETAINED_INTRON": "The base's own exon spans an intron in one form that the other form splices out — intron retention right at the base (both forms share the same 3' end).",
    "IPA_EXTENSION": "The base's own exon is extended in one form because it reads into the downstream intron and polyadenylates there (intronic polyadenylation), while the other form splices on to a more distal 3' end.",
    "NEAR_ALT_JUNCTION": "The base sits within the exon-junction-complex footprint of a splice junction present in one form and removed in the other (EJC relief).",
    "DISTAL_APA": "The base's local exon is identical in both forms, but the forms end at DIFFERENT poly(A) sites (a distal alternative-polyadenylation choice); methylation tracks which 3' end the isoform uses.",
    "DISTAL_SPLICING": "The base's local exon is identical AND the forms share the same 3' end; they differ only in an internal / 5' splicing choice ELSEWHERE, so methylation tracks that distal splicing decision.",
    "FIVE_PRIME_UNCERTAIN": "The base falls in (or 5' of) the 5'-most exon, where direct-RNA 5' truncation makes the assembled model unreliable — no confident structural call.",
    "INTRON_READ_ARTIFACT": "The more-modified form does not structurally contain the base (it is intronic there) — likely an intron / soft-clip read, not a real exonic modification.",
    "NO_MODEL": "Fewer than two adequately-covered fragmentforms at the site, so there is nothing to compare.",
    "UNRESOLVED": "Passes the filters but fits no structural rule.",
}


def build_classification_section(class_df, private_df, class_figs_dir, arch_figs_dir, max_figs_per_category, top_n):
    """Tree view of the differential-site classification. The PRIVATE bucket is sourced from the
    coverage-INDEPENDENT private scan (private_df), which finds sites where the base is absent from
    a fragmentform even when that form has no coverage there; the SHARED / UNEXPLAINABLE buckets come
    from the differential test (class_df). Every bucket and event is always listed (with "No sites"
    for empty ones)."""
    title = ("Classification of Transcript Architecture Changes Associated with Differential "
             "Epitranscriptomic Modification")
    intro = (
        "Sites are placed in a <b>three-level tree</b>. <b>Level 1 (bucket)</b>: does the base exist in "
        "the mature RNA of both forms — <b>PRIVATE</b> (base only in the higher form) vs <b>SHARED</b>. "
        "<b>Level 2 (event)</b> for shared sites: is the distinguishing structural change ON the base's "
        "exon (<b>SHARED_LOCAL</b> → ALT_DONOR / ALT_ACCEPTOR / ALT_POLYA_SITE / NEAR_ALT_JUNCTION) or "
        "ELSEWHERE in the transcript (<b>SHARED_DISTAL</b>)? SHARED_DISTAL now splits into <b>DISTAL_APA</b> "
        "(the forms differ in their 3' end / poly(A) site) vs <b>DISTAL_SPLICING</b> (same 3' end, they "
        "differ in an internal / 5' splicing choice). <b>Level 3 (direction)</b>: which form is the "
        "more-modified one, worded per event (e.g. LONGER_EXON_HIGHER, PROXIMAL_HIGHER). A 4th bucket, "
        "<b>UNEXPLAINABLE</b>, holds sites with no confident structural cause. Every bucket, event and "
        "direction is ALWAYS listed — with a zero count where nothing landed there — so the taxonomy is "
        "reported completely. PRIVATE sites are detected by a COVERAGE-INDEPENDENT structural scan (every "
        "modified site checked against all fragmentform models of its gene) — so they are found even when "
        "the form lacking the base has no coverage there, which the differential test cannot see; SHARED / "
        "UNEXPLAINABLE come from the between-fragmentform differential test (BH-FDR + the &gt;10% effect "
        "rule, tunable via classify_diffs.min_effect)."
    )
    class_ok = class_df is not None and not class_df.empty and "bucket" in class_df.columns
    priv_ok = private_df is not None and not private_df.empty and "bucket" in private_df.columns
    if not class_ok and not priv_ok:
        return section(title, "<p class='muted'>No classified differential sites available.</p>", intro=intro)

    SHARED_COLS = [c for c in [
        "gene_name", "mod_code", "chrom", "start0", "strand", "direction", "structural_delta_nt",
        "hi_ZN", "hi_arch", "hi_frac", "lo_ZN", "lo_arch", "lo_frac", "stoich_tier", "hi_stoich_level",
        "effect_max_abs_frac_diff", "p_adj_bh",
    ] if class_ok and c in class_df.columns]
    PRIV_COLS = [c for c in [
        "gene_name", "mod_code", "chrom", "start0", "strand", "direction", "structural_delta_nt",
        "carry_ZN", "carry_arch", "carry_frac", "carry_cov", "n_forms_present", "n_forms_absent", "absent_in_ZN",
    ] if priv_ok and c in private_df.columns]

    def src_for(bucket):
        return private_df if bucket == "PRIVATE" else class_df

    def count(bucket, event=None):
        df = src_for(bucket)
        if df is None or df.empty or "bucket" not in df.columns:
            return 0
        m = df["bucket"] == bucket
        if event is not None:
            m &= df["event"] == event
        return int(m.sum())

    # Which rows the known bucket/event enumeration will actually render (PRIVATE from private_df,
    # the rest from class_df). Anything NOT covered here is a taxonomy-drift row (e.g. the classifier
    # emitted a newer label than this report knows) — surface it in a catch-all rather than let it
    # vanish and skew the denominator.
    def _rendered_mask(df, buckets):
        if df is None or df.empty or "bucket" not in df.columns:
            return None
        m = pd.Series(False, index=df.index)
        for b in buckets:
            for e in CLASS_TAXONOMY.get(b, []):
                sub = (df["bucket"] == b)
                if "event" in df.columns:
                    sub &= (df["event"] == e)
                m |= sub
        return m

    leftover_frames = []
    pr = _rendered_mask(private_df, ["PRIVATE"])
    if pr is not None and (~pr).any():
        leftover_frames.append(private_df[~pr])
    cr = _rendered_mask(class_df, [b for b in CLASS_BUCKET_ORDER if b != "PRIVATE"])
    if cr is not None and (~cr).any():
        leftover_frames.append(class_df[~cr])
    leftover_df = pd.concat(leftover_frames, ignore_index=True) if leftover_frames else pd.DataFrame()

    known_total = sum(count(b, e) for b in CLASS_BUCKET_ORDER for e in CLASS_TAXONOMY[b])
    total = (known_total + len(leftover_df)) or 1

    def count_dir(bucket, event, direction):
        df = src_for(bucket)
        if df is None or df.empty or "bucket" not in df.columns:
            return 0
        m = (df["bucket"] == bucket) & (df["event"] == event)
        if "direction" in df.columns:
            m &= df["direction"].fillna("") == direction
        elif direction:
            return 0
        return int(m.sum())

    # ---- overview: complete bucket x event x direction count grid (all leaves, zeros included) ----
    grid_rows = []
    for b in CLASS_BUCKET_ORDER:
        for e in CLASS_TAXONOMY[b]:
            dlist = EVENT_DIRECTIONS.get(e, [])
            if dlist:
                for d in dlist:
                    n = count_dir(b, e, d)
                    grid_rows.append({"bucket": b, "event": e, "direction": d,
                                      "n_sites": n, "pct": f"{100.0 * n / total:.1f}%"})
            else:
                n = count(b, e)
                grid_rows.append({"bucket": b, "event": e, "direction": "—",
                                  "n_sites": n, "pct": f"{100.0 * n / total:.1f}%"})
    grid = pd.DataFrame(grid_rows)
    bucket_counts = {b: count(b) for b in CLASS_BUCKET_ORDER}
    hero = category_distribution_png(bucket_counts, mod_label=None)
    hero_html = clickable_image_html(hero, "Classification by top-level bucket", figure_class="hero-figure",
                                     caption="Sites per top-level bucket.") if hero else ""
    overview = (
        "<div class='overview-layout'>"
        f"<div>{df_to_html(grid, max_rows=len(grid))}"
        f"{definitions_html(list(BUCKET_DEFINITIONS.items()), summary='Level 1-2: bucket definitions', open_by_default=False)}"
        f"{definitions_html(list(EVENT_DEFINITIONS.items()), summary='Level 2: event definitions (all events)', open_by_default=False)}"
        f"{definitions_html(list(DIRECTION_DEFINITIONS.items()), summary='Level 3: direction definitions (all directions)', open_by_default=False)}"
        "</div>"
        f"<div class='hero'>{hero_html or ''}</div>"
        "</div>"
    )

    def render_sites(sub, bucket, cols):
        """Table (+ architecture/mechanism figures for shared buckets) for a leaf set of sites."""
        sort_col = "carry_frac" if bucket == "PRIVATE" else "effect_max_abs_frac_diff"
        srt = sub.sort_values(sort_col, ascending=False) if sort_col in sub.columns else sub
        tbl = df_to_html(srt[cols] if cols else srt, max_rows=top_n)
        figs = ""
        if bucket != "PRIVATE" and "class_key" in sub.columns:
            for ck in sub["class_key"].dropna().unique():
                figs += category_figure_gallery(arch_figs_dir, ck, max_figs_per_category,
                                                summary_label="isoform architecture map(s) — exon/intron tracks, site marked",
                                                open_by_default=False)
                figs += category_figure_gallery(class_figs_dir, ck, max_figs_per_category)
        defs = definitions_html(column_definitions(cols), summary="Column definitions") if cols else ""
        return defs + tbl + figs

    # ---- the tree: bucket -> event -> direction -> (sites [+ figures]) ----
    tree = []
    for bucket in CLASS_BUCKET_ORDER:
        df = src_for(bucket)
        cols = PRIV_COLS if bucket == "PRIVATE" else SHARED_COLS
        bdf = df[df["bucket"] == bucket] if (df is not None and not df.empty and "bucket" in df.columns) else pd.DataFrame()
        n_b = len(bdf)
        events_html = []
        for event in CLASS_TAXONOMY[bucket]:
            edf = bdf[bdf["event"] == event] if not bdf.empty else bdf
            n_e = len(edf)
            defn = EVENT_DEFINITIONS.get(event, "")
            dir_list = EVENT_DIRECTIONS.get(event, [])
            if dir_list:
                # Level 3: one collapsible subheader per direction, ALL enumerated (incl. zeros)
                seen = set(dir_list)
                extra = ([d for d in edf["direction"].fillna("").unique() if d and d not in seen]
                         if (not edf.empty and "direction" in edf.columns) else [])
                dirs_html = []
                for direction in list(dir_list) + list(extra):
                    ddf = (edf[edf["direction"].fillna("") == direction]
                           if (not edf.empty and "direction" in edf.columns) else edf.iloc[0:0])
                    n_d = len(ddf)
                    ddef = DIRECTION_DEFINITIONS.get(direction, "")
                    body = (render_sites(ddf, bucket, cols) if n_d
                            else "<p class='muted'>No sites of this direction in this run.</p>")
                    dirs_html.append(
                        "<details class='report-subtree'>"
                        f"<summary><h4>{html.escape(direction)} — {n_d} site{'s' if n_d != 1 else ''}</h4></summary>"
                        f"<div class='section-body'><p class='section-intro'>{html.escape(ddef)}</p>{body}</div>"
                        "</details>"
                    )
                inner = "".join(dirs_html)
            else:
                inner = (render_sites(edf, bucket, cols) if n_e
                         else "<p class='muted'>No sites in this category in this run.</p>")
            events_html.append(
                "<details class='report-subtree'>"
                f"<summary><h3>{html.escape(event)} — {n_e} site{'s' if n_e != 1 else ''}</h3></summary>"
                f"<div class='section-body'><p class='section-intro'>{html.escape(defn)}</p>{inner}</div>"
                "</details>"
            )
        bdef = BUCKET_DEFINITIONS.get(bucket, "")
        tree.append(
            f"<details class='report-subtree'{' open' if n_b else ''}>"
            f"<summary><h2 class='bucket'>{html.escape(bucket)} — {n_b} site{'s' if n_b != 1 else ''}</h2></summary>"
            f"<div class='section-body'><p class='section-intro'>{html.escape(bdef)}</p>{''.join(events_html)}</div>"
            "</details>"
        )

    # catch-all: taxonomy-drift rows not in the known bucket/event constants — listed, never dropped
    if not leftover_df.empty:
        nlo = len(leftover_df)
        lo_cols = [c for c in ["gene_name", "mod_code", "chrom", "start0", "strand",
                               "bucket", "event", "direction", "structural_delta_nt"]
                   if c in leftover_df.columns]
        tree.append(
            "<details class='report-subtree' open>"
            f"<summary><h2 class='bucket'>UNRECOGNIZED — {nlo} site{'s' if nlo != 1 else ''}</h2></summary>"
            "<div class='section-body'><p class='section-intro'>Sites whose bucket/event is not in this "
            "report's taxonomy constants (e.g. the classifier emitted a newer label than the report "
            "knows). They are listed here rather than silently dropped, so the counts and percentages "
            "above stay complete.</p>"
            + df_to_html(leftover_df[lo_cols] if lo_cols else leftover_df, max_rows=top_n)
            + "</div></details>")

    return section(title, overview + "".join(tree), intro=intro)


def externalize_data_uris(html_doc, out_html):
    """Replace inline base64 image data: URIs with sidecar files referenced by relative
    path. Chrome blocks top-level navigation to data: URLs (so clicking a figure to
    open it full-size only worked in Firefox), and inlining many multi-MB PNGs bloats
    the HTML past 100 MB. Writing images next to the report fixes both; identical
    images are de-duplicated. Returns (new_html, n_images, rel_dir)."""
    out_dir = os.path.dirname(os.path.abspath(out_html)) or "."
    stem = os.path.splitext(os.path.basename(out_html))[0]
    rel_dir = f"{stem}_files"
    abs_dir = os.path.join(out_dir, rel_dir)
    seen = {}
    counter = itertools.count(1)
    made = {"n": 0}
    pattern = re.compile(r"data:image/(png|jpeg|jpg|gif|svg\+xml);base64,([A-Za-z0-9+/=]+)")

    def repl(match):
        fmt, b64 = match.group(1), match.group(2)
        digest = hashlib.md5(b64.encode("ascii")).hexdigest()
        rel = seen.get(digest)
        if rel is None:
            try:
                raw = base64.b64decode(b64)
            except Exception:
                return match.group(0)
            os.makedirs(abs_dir, exist_ok=True)
            ext = {"jpeg": "jpg", "svg+xml": "svg"}.get(fmt, fmt)
            name = f"fig_{next(counter):04d}.{ext}"
            with open(os.path.join(abs_dir, name), "wb") as fh:
                fh.write(raw)
            rel = f"{rel_dir}/{name}"
            seen[digest] = rel
            made["n"] += 1
        return rel

    return pattern.sub(repl, html_doc), made["n"], rel_dir


def main():
    args = parse_args()

    # Persist the report's inline summary charts as PNG + SVG files next to the report.
    global _REPORT_FIGS_DIR
    _REPORT_FIGS_DIR = os.path.join(os.path.dirname(os.path.abspath(args.out_html)) or ".", "figures")
    try:
        os.makedirs(_REPORT_FIGS_DIR, exist_ok=True)
    except OSError:
        _REPORT_FIGS_DIR = None

    class_df = read_tsv(args.classification)
    metrics_df = read_tsv(args.metrics)
    tx_counts_df = read_tsv(args.tx_counts)
    sample_stats_df = read_tsv(args.sample_stats)
    read_stats_df = read_tsv(args.read_stats)
    tx_lengths_df = read_tsv(args.tx_lengths)
    partition_map_df = read_tsv(args.partition_map)
    zn_long_df = read_tsv(args.zn_long)
    zt_long_df = read_tsv(args.zt_long)
    # Per-fragmentform (per-ZN) exon models, for the cis-SNP→mod stoichiometry per-fragmentform graphs.
    snp_mod_iso = {}
    if getattr(args, "assembled_gtf", "") and load_isoforms is not None and os.path.exists(args.assembled_gtf):
        try:
            snp_mod_iso, _snp_mod_genes = load_isoforms(args.assembled_gtf, tes_tol=25, inside_tol=50)
        except Exception:
            snp_mod_iso = {}
    diff_df = read_tsv(args.diff_results)
    # Data-driven "significance is cheap" callout, shown above every differential section.
    meta_df = read_tsv(args.sample_metadata) if getattr(args, "sample_metadata", "") else pd.DataFrame()
    sig_box = significance_note_box(replicate_concordance(zn_long_df, meta_df))
    classified_df = read_tsv(args.classified_sites)
    overlap_df = read_summary_metrics(glob.glob(args.multigene_summary_glob)) if args.multigene_summary_glob else pd.DataFrame()
    candidate_snps_df = read_tsv(args.candidate_snps)
    snp_tx_assoc_df = read_tsv(args.snp_tx_assoc)
    snp_mod_assoc_df = read_tsv(args.snp_mod_assoc)
    hap_blocks_df = read_tsv(args.hap_blocks)
    hap_tx_assoc_df = read_tsv(args.hap_tx_assoc)
    hap_mod_assoc_df = read_tsv(args.hap_mod_assoc)

    # Per-example (per-row) figures for the genotype/SNP sections.
    try:
        import snp_report_figures
        snp_galleries = snp_report_figures.build_snp_galleries(
            snp_tx=snp_tx_assoc_df, snp_mod=snp_mod_assoc_df,
            hap_tx=hap_tx_assoc_df, hap_mod=hap_mod_assoc_df,
            figs_dir=args.snp_figs_dir, max_figs=args.max_snp_figs,
        )
    except Exception:
        snp_galleries = {}

    overview_cards = []
    n_tx = len(class_df)
    n_genes = class_df["gene_index"].nunique() if "gene_index" in class_df.columns and not class_df.empty else 0
    n_metagenes = class_df["metagene_index"].nunique() if "metagene_index" in class_df.columns and not class_df.empty else 0
    def _num_col(col, reduce="sum"):
        # numeric reduction that tolerates non-numeric cells (coerce->NaN) instead of aborting the
        # whole report on a stray 'NA'; also avoids a lexical max on an object-typed count column.
        if col not in class_df.columns or class_df.empty:
            return 0
        s = pd.to_numeric(class_df[col], errors="coerce")
        v = s.max() if reduce == "max" else s.sum()
        return 0 if pd.isna(v) else v
    max_partitions = _num_col("metagene_partition_count", "max")
    total_reads = _num_col("read_support")
    trunc_reads = _num_col("trunc_assigned_reads")
    exact_reads = _num_col("exact_chain_reads")
    overview_cards.extend([
        ("Fragmentforms", fmt_int(n_tx)),
        ("Genes", fmt_int(n_genes)),
        ("Metagenes", fmt_int(n_metagenes)),
        ("Max ZN partitions / metagene", fmt_int(max_partitions)),
        ("Assigned reads", fmt_int(total_reads)),
        ("Exact-chain reads", fmt_int(exact_reads)),
        ("Truncation-assigned reads", fmt_int(trunc_reads)),
    ])
    if not candidate_snps_df.empty:
        overview_cards.extend([
            ("Segregating SNPs", fmt_int(len(candidate_snps_df))),
            ("Haplotype blocks", fmt_int(len(hap_blocks_df) if not hap_blocks_df.empty else 0)),
        ])

    pca_img = embed_png(args.pca_png)
    pca_html = clickable_image_html(pca_img, "PCA plot", figure_class="hero-figure") if pca_img else "<p class='muted'>PCA plot unavailable.</p>"

    top_tx_cols = [
        c for c in [
            "zt_label", "gtf_gene_name", "read_support", "exact_chain_reads",
            "trunc_assigned_reads", "anchor_reads", "anchor_frac",
            "metagene_index", "zn_index", "classification"
        ] if c in class_df.columns
    ]
    top_tx_df = class_df.sort_values(["read_support", "exact_chain_reads"], ascending=False) if not class_df.empty else class_df

    top_gene_sites_df = pd.DataFrame()
    if not zn_long_df.empty and "gene_name" in zn_long_df.columns:
        top_gene_sites_df = (
            zn_long_df.assign(site_key=zn_long_df[["chrom", "start0", "end0", "strand", "mod_code"]].astype(str).agg(":".join, axis=1))
            .groupby(["gene_name", "mod_code"], as_index=False)["site_key"].nunique()
            .rename(columns={"site_key": "n_sites"})
            .sort_values("n_sites", ascending=False)
        )

    diff_html = "<p class='muted'>No differential-site results available.</p>"
    diff_cols = []
    if not diff_df.empty:
        diff_cols = [
            c for c in [
                "gene_name", "mod_code", "chrom", "start0", "end0", "strand",
                "p_value", "p_adj_bh", "effect_max_abs_frac_diff"
            ] if c in diff_df.columns
        ]
        diff_html = df_to_html(
            diff_df[diff_cols].sort_values(["p_adj_bh", "effect_max_abs_frac_diff"], ascending=[True, False]),
            max_rows=20,
        )

    diff_fig_html = ""
    if args.diff_figs_dir and os.path.isdir(args.diff_figs_dir):
        fig_paths = sorted(glob.glob(os.path.join(args.diff_figs_dir, "*.png")))[: args.max_diff_figs]
        if fig_paths:
            pieces = []
            for path in fig_paths:
                img = embed_png(path)
                if img:
                    pieces.append(clickable_image_html(img, os.path.basename(path), caption=os.path.basename(path)))
            diff_fig_html = "<div class='gallery'>" + "".join(pieces) + "</div>"

    read_stats_counts_cols = [
        c for c in [
            "sample", "total_reads_bam", "total_mapped", "total_unmapped", "considered_reads",
            "failed_unmapped", "failed_secondary_or_supp", "failed_low_mapq", "failed_low_introns",
            "failed_low_softclip3p", "zt_tagged_exists", "zt_total_records", "zt_unmapped_records",
            "zt_mapped_records", "assigned_reads", "zt_mapped_unassigned_reads"
        ] if c in read_stats_df.columns
    ]
    read_stats_length_cols = ["sample"] + [c for c in read_stats_df.columns if "_len_" in c]
    read_stats_length_cols = [c for c in read_stats_length_cols if c in read_stats_df.columns]

    candidate_snps_df_view = pd.DataFrame()
    if not candidate_snps_df.empty:
        keep = [
            c for c in [
                "snp_id", "gene_names", "metagene_indices", "total_cov",
                "ref_count", "alt_count", "alt_frac", "samples_with_alt"
            ] if c in candidate_snps_df.columns
        ]
        candidate_snps_df_view = candidate_snps_df[keep].sort_values(["alt_count", "alt_frac"], ascending=False)

    snp_tx_df_view = pd.DataFrame()
    if not snp_tx_assoc_df.empty:
        keep = [
            c for c in [
                "snp_id", "gene_names", "metagene_indices", "n_alt_reads",
                "n_transcripts_tested", "p_value", "p_adj_bh", "effect_max_abs_tx_frac_diff"
            ] if c in snp_tx_assoc_df.columns
        ]
        snp_tx_df_view = snp_tx_assoc_df[keep].sort_values(["p_adj_bh", "effect_max_abs_tx_frac_diff"], ascending=[True, False])

    snp_mod_df_view = pd.DataFrame()
    if not snp_mod_assoc_df.empty:
        keep = [
            c for c in [
                "snp_id", "mod_site_id", "target_mod_code", "gene_names", "n_alt_reads",
                "n_modified", "p_value", "p_adj_bh", "effect_abs_delta_mod_frac"
            ] if c in snp_mod_assoc_df.columns
        ]
        snp_mod_df_view = snp_mod_assoc_df[keep].sort_values(["p_adj_bh", "effect_abs_delta_mod_frac"], ascending=[True, False])

    hap_blocks_df_view = pd.DataFrame()
    if not hap_blocks_df.empty:
        keep = [c for c in ["block_id", "gene_names", "context_key", "chrom", "region", "span_bp", "n_snps", "snp_coords", "support_reads", "complete_reads", "haplotypes"] if c in hap_blocks_df.columns]
        hap_blocks_df_view = hap_blocks_df[keep].sort_values(["complete_reads", "n_snps"], ascending=False)

    # block_id -> genomic coordinate string, so the association tables can carry the hapblock's location
    block_region = {}
    if not hap_blocks_df.empty and "block_id" in hap_blocks_df.columns:
        for _, b in hap_blocks_df.iterrows():
            if "region" in hap_blocks_df.columns and pd.notna(b.get("region")):
                block_region[b["block_id"]] = str(b["region"])
            elif {"chrom", "start1", "end1"}.issubset(hap_blocks_df.columns):
                block_region[b["block_id"]] = f"{b['chrom']}:{int(b['start1'])}-{int(b['end1'])}"

    def _with_region(df):
        if "block_id" not in df.columns:   # malformed input without the join key -> pass through
            return df.copy()
        out = df.copy()
        out.insert(1, "block_region", out["block_id"].map(block_region).fillna(""))
        return out

    hap_tx_df_view = pd.DataFrame()
    if not hap_tx_assoc_df.empty:
        keep = [c for c in ["block_id", "n_transcripts_tested", "p_value", "p_adj_bh", "effect_max_abs_tx_frac_diff"] if c in hap_tx_assoc_df.columns]
        hap_tx_df_view = _with_region(hap_tx_assoc_df[keep]).sort_values(["p_adj_bh", "effect_max_abs_tx_frac_diff"], ascending=[True, False])

    hap_mod_df_view = pd.DataFrame()
    if not hap_mod_assoc_df.empty:
        keep = [c for c in ["block_id", "mod_site_id", "target_mod_code", "p_value", "p_adj_bh", "effect_max_abs_mod_rate_diff"] if c in hap_mod_assoc_df.columns]
        hap_mod_df_view = _with_region(hap_mod_assoc_df[keep]).sort_values(["p_adj_bh", "effect_max_abs_mod_rate_diff"], ascending=[True, False])

    cards_html = "".join(
        f"<div class='card'><div class='label'>{html.escape(label)}</div><div class='value'>{html.escape(value)}</div></div>"
        for label, value in overview_cards
    )

    overview_defs = definitions_html(
        visible_definitions([label for label, _ in overview_cards]),
        summary="Definitions for the overview metrics",
        open_by_default=True,
    )

    # ---- build each section into a variable, then assemble in the intended reading order ----
    sec_overview = section(
        "Overview",
        (
            "<div class='overview-layout'>"
            "<div>"
            f"<div class='cards'>{cards_html}</div>{overview_defs}"
            "</div>"
            f"<div class='hero'>{pca_html}</div>"
            "</div>"
        ),
        intro="Total counts of the objects modulator assembled and tested in this run — fragmentforms, genes, "
              "metagenes, segregating SNPs — and the reads supporting them.",
    )

    # Read-length summaries, reordered to MATCH the "Most Expressed Fragmentforms" table (same
    # fragmentforms, same highest-expression order) so the two line up row-for-row.
    tx_lengths_ordered = tx_lengths_df
    if (not tx_lengths_df.empty and "zt_label" in tx_lengths_df.columns
            and not top_tx_df.empty and "zt_label" in top_tx_df.columns):
        _have = set(tx_lengths_df["zt_label"])
        _order = [z for z in top_tx_df["zt_label"].tolist() if z in _have]
        if _order:
            tx_lengths_ordered = tx_lengths_df.set_index("zt_label").reindex(_order).reset_index()
    sec_top_ff = section(
        "Most Expressed Fragmentforms",
        (df_to_html(top_tx_df[top_tx_cols], max_rows=args.top_transcripts) if top_tx_cols
         else "<p class='muted'>No fragmentform summary columns available.</p>")
        + (subsection(
            "Read Lengths for the Most Expressed Fragmentforms",
            df_to_html(tx_lengths_ordered, max_rows=args.top_transcripts),
            definitions=definitions_html(column_definitions(list(tx_lengths_ordered.columns)), summary="Column definitions"),
          ) if not tx_lengths_ordered.empty else ""),
        intro="The highest-support fragmentforms retained after assembly and reference classification, ordered by "
              "assigned reads.",
        definitions=definitions_html(column_definitions(top_tx_cols), summary="Column definitions"),
    )

    sec_sample_stats = section(
        "Per-Sample Fragmentform Statistics",
        df_to_html(sample_stats_df, max_rows=100),
        intro="Per-sample fragmentform assignment totals and detection breadth derived from the final retained fragmentforms.",
        definitions=definitions_html(column_definitions(list(sample_stats_df.columns)), summary="Column definitions"),
    )

    sec_read_funnel = section(
        "Per-Sample Read Usage & Length Statistics",
        (
            subsection(
                "Read Retention Counts",
                df_to_html(read_stats_df[read_stats_counts_cols], max_rows=100) if read_stats_counts_cols else "<p class='muted'>No read-count funnel columns available.</p>",
                definitions=definitions_html(column_definitions(read_stats_counts_cols), summary="Column definitions") if read_stats_counts_cols else "",
            ) +
            subsection(
                "Read Length Summaries",
                df_to_html(read_stats_df[read_stats_length_cols], max_rows=100) if len(read_stats_length_cols) > 1 else "<p class='muted'>No read-length summary columns available.</p>",
                definitions=definitions_html(column_definitions(read_stats_length_cols), summary="Column definitions") if len(read_stats_length_cols) > 1 else "",
            )
        ),
        intro="Sample-level retention across the primary read filters and fragmentform-assignment workflow.",
    )

    sec_gene_sites = section(
        "Genes with the Most Number of Modified Sites",
        df_to_html(top_gene_sites_df, max_rows=args.top_genes),
        intro="Counts of unique genomic modification sites observed per gene and modification code in the ZN aggregation.",
        definitions=definitions_html(column_definitions(list(top_gene_sites_df.columns)), summary="Column definitions"),
    )

    # Overlap resolution: put `sample` first and drop the pseudo-metric header rows.
    if not overlap_df.empty:
        overlap_view = overlap_df.drop(columns=[c for c in ["removed_reads_per_gene_id", "removed_reads_per_zt_label"] if c in overlap_df.columns])
        _ocols = (["sample"] if "sample" in overlap_view.columns else []) + [c for c in overlap_view.columns if c != "sample"]
        overlap_view = overlap_view[_ocols].copy()
        overlap_body = df_to_html(overlap_view.fillna("0"), max_rows=100)
        overlap_defs2 = definitions_html(column_definitions(list(overlap_view.columns)), summary="Column definitions")
    else:
        overlap_body = "<p class='muted'>No overlap-resolution summaries available.</p>"
        overlap_defs2 = ""
    sec_overlap = section(
        "Read Usage within Overlapping Gene Loci",
        overlap_body,
        intro="Where genes overlap on the same strand, modulator merges them into one metagene and resolves each "
              "read to a single gene by its fragmentform (ZT) assignment. This shows, per sample, how the reads in "
              "those overlapping-gene regions were used or scrapped (see input parameter multi_gene_action).",
        definitions=overlap_defs2,
    )

    sec_diff = section(
        "Sites with Differential Epitranscriptomic Modification Between Fragmentforms",
        diff_html + diff_fig_html,
        intro="Positions where the modification stoichiometry differs between the fragmentforms of a gene "
              "(across all detected mod codes; Fisher exact / chi-square with Benjamini-Hochberg FDR, keeping "
              "sites whose absolute stoichiometry difference clears the minimum effect). The effect threshold "
              "(default 10%) and the FDR are tunable — see input parameters classify_diffs.min_effect and "
              "classify_diffs.fdr.",
        definitions=definitions_html(column_definitions(diff_cols), summary="Result-column definitions") if diff_cols else "",
    )

    private_df = read_tsv(args.private_sites) if getattr(args, "private_sites", "") else pd.DataFrame()
    sec_class = build_classification_section(
        classified_df, private_df, args.class_figs_dir, args.arch_figs_dir,
        args.max_class_figs_per_category, args.max_class_figs_per_category,
    )

    # --- Splice-junction repertoire (canonical vs non-canonical), read-derived ---
    splice_genes_df = read_tsv(args.splice_genes) if args.splice_genes else pd.DataFrame()
    splice_junc_df = read_tsv(args.splice_junctions) if args.splice_junctions else pd.DataFrame()
    sec_splice = None
    if not splice_junc_df.empty and "junction_class" in splice_junc_df.columns:
        cls_counts = splice_junc_df["junction_class"].value_counts()
        total_j = int(cls_counts.sum())
        overview = (
            "<h3>Summary Across All Measured Junctions</h3><ul>" + "".join(
                f"<li><b>{html.escape(str(k))}</b>: {int(v):,} ({100.0 * int(v) / total_j:.2f}%)</li>"
                for k, v in cls_counts.items()
            ) + "</ul>"
            + "<p class='muted'>The exact coordinates and donor/acceptor dinucleotide of every junction — and "
              "the specific genes carrying non-canonical junctions — are in the per-gene summary and "
              "non-canonical tables below.</p>"
        )
        noncanon_genes = splice_genes_df[pd.to_numeric(
            splice_genes_df.get("has_noncanonical", 0), errors="coerce").fillna(0).astype(int) == 1] \
            if not splice_genes_df.empty and "has_noncanonical" in splice_genes_df.columns else pd.DataFrame()
        noncanon_html = (
            subsection("Genes carrying non-canonical junctions", df_to_html(noncanon_genes, max_rows=args.top_genes))
            if not noncanon_genes.empty
            else "<p class='muted'>Every gene's junctions are canonical (GT-AG) or semi-canonical (GC-AG).</p>"
        )
        sec_splice = section(
            "Splice Junctions in Fragmentforms",
            overview
            + (subsection("Per-gene summary", df_to_html(splice_genes_df, max_rows=args.top_genes)) if not splice_genes_df.empty else "")
            + noncanon_html,
            intro="Donor/acceptor dinucleotides of every assembled fragmentform intron, in fragmentform "
                  "orientation. GT-AG is the major (U2) spliceosome; GC-AG is semi-canonical; AT-AC is the "
                  "minor (U12) spliceosome; anything else is non-canonical and worth inspecting. Because the "
                  "intron chains are read-derived, these are the junctions the reads actually support.",
            definitions=definitions_html(column_definitions(list(splice_genes_df.columns)), summary="Column definitions") if not splice_genes_df.empty else "",
        )

    # --- Novel gene loci (fragmentforms matching no reference gene) ---
    novel_loci_df = read_tsv(args.novel_loci) if args.novel_loci else pd.DataFrame()
    novel_ff_df = read_tsv(args.novel_fragmentforms) if args.novel_fragmentforms else pd.DataFrame()
    sec_novel = None
    if args.novel_loci:
        if not novel_loci_df.empty:
            nl_body = (
                subsection("Novel loci", df_to_html(novel_loci_df, max_rows=args.top_genes))
                + (subsection("Fragmentforms of novel loci", df_to_html(novel_ff_df, max_rows=args.top_transcripts)) if not novel_ff_df.empty else "")
            )
        else:
            nl_body = "<p class='muted'>No novel loci: every assembled fragmentform matched a reference gene.</p>"
        sec_novel = section(
            "Novel Gene Loci",
            nl_body,
            intro="Read-backed loci whose fragmentforms overlap no annotated gene in the reference GTF. Each locus "
                  "is formed by merging overlapping novel fragmentforms on one strand and is given a unique, "
                  "deterministic, coordinate-anchored name (NOVEL_<chrom>_<strand>_<n>_<start>_<end>), so two "
                  "distinct novel loci on the same chromosome and strand can never be conflated.",
            definitions=definitions_html(column_definitions(list(novel_loci_df.columns)), summary="Column definitions") if not novel_loci_df.empty else "",
        )

    # --- APA motifs (conditional) ---
    sec_apa = build_apa_motif_section(read_tsv(args.apa_motifs), args.top_genes) if args.apa_motifs else None

    # --- Sequence elements (conditional) ---
    sec_seq = build_sequence_elements_section(
        read_tsv(args.sequence_elements),
        read_tsv(args.sequence_elements_summary) if args.sequence_elements_summary else pd.DataFrame(),
        args.top_genes) if args.sequence_elements else None

    # --- Poly(A) tail length (conditional) ---
    polya_frag_df = read_tsv(args.polya_fragmentform) if args.polya_fragmentform else pd.DataFrame()
    taillength_diffs_df = read_tsv(args.taillength_diffs) if args.taillength_diffs else pd.DataFrame()
    taillength_mod_df = read_tsv(args.taillength_mod) if args.taillength_mod else pd.DataFrame()
    sec_polya = None
    if args.polya_fragmentform or args.taillength_diffs or args.taillength_mod:
        sec_polya = build_polya_section(polya_frag_df, taillength_diffs_df, taillength_mod_df, args.top_genes,
                                        diff_figs_dir=args.taillength_diff_figs, mod_figs_dir=args.taillength_mod_figs,
                                        max_figs=int(getattr(args, "max_snp_figs", 12)))

    # --- Between-condition comparisons (conditional) ---
    sec_between = build_between_conditions_section(args.between_conditions_dir, args.top_genes, top_note=sig_box) if args.between_conditions_dir else None

    sec_snp_cand = section(
        "Segregating cis SNP Candidates",
        df_to_html(candidate_snps_df_view, max_rows=args.top_genes) if not candidate_snps_df_view.empty else "<p class='muted'>No segregating SNP candidates available.</p>",
        intro="Read-supported non-reference loci discovered in the cleaned tagged BAMs.",
        definitions=definitions_html(column_definitions(list(candidate_snps_df_view.columns)), summary="Column definitions") if not candidate_snps_df_view.empty else "",
    )
    sec_snp_ff = section(
        "cis SNP to Fragmentform Usage Associations",
        (df_to_html(snp_tx_df_view, max_rows=args.top_genes) if not snp_tx_df_view.empty else "<p class='muted'>No cis SNP to fragmentform-usage associations available.</p>") + snp_galleries.get("snp_tx", ""),
        intro="Associations between segregating SNP alleles and fragmentform-partition usage.",
        definitions=definitions_html(column_definitions(list(snp_tx_df_view.columns)), summary="Column definitions") if not snp_tx_df_view.empty else "",
    )
    snp_mod_ff_html = build_snp_mod_fragmentform_html(snp_mod_assoc_df, snp_mod_iso,
                                                      getattr(args, "molecule_mod_calls", ""),
                                                      getattr(args, "molecule_snps", ""), args.max_snp_figs)
    sec_snp_mod = section(
        "cis SNP to Modification Stoichiometry Association",
        (df_to_html(snp_mod_df_view, max_rows=args.top_genes) if not snp_mod_df_view.empty else "<p class='muted'>No cis SNP to modification-stoichiometry associations available.</p>") + snp_galleries.get("snp_mod", "") + snp_mod_ff_html,
        intro="Associations between segregating SNP alleles and target modification states on the same molecules.",
        definitions=definitions_html(column_definitions(list(snp_mod_df_view.columns)), summary="Column definitions") if not snp_mod_df_view.empty else "",
    )

    # --- Distance of cis SNP to affected modifications (positional + motif mechanism) ---
    sec_snp_mech = build_snp_mechanism_section(read_tsv(args.snp_mod_mechanism), args.top_genes) if args.snp_mod_mechanism else None

    # --- Single-molecule co-localized modifications (mod x mod dependency on shared molecules) ---
    mod_mod_df = read_tsv(args.mod_mod_assoc) if args.mod_mod_assoc else pd.DataFrame()
    sec_modmod = None
    if args.mod_mod_assoc:
        if not mod_mod_df.empty:
            mm_view = mod_mod_df
            if "p_adj_bh" in mm_view.columns:
                n_sig = int((pd.to_numeric(mm_view["p_adj_bh"], errors="coerce") < 0.05).sum())
                header = f"<p><b>{len(mm_view):,}</b> co-localized mod-site pairs tested; <b>{n_sig:,}</b> significant at BH-FDR &lt; 0.05.</p>"
            else:
                header = ""
            if "direction" in mm_view.columns:
                dc = mm_view["direction"].value_counts().to_dict()
                header += ("<p>direction: "
                           + " · ".join(f"<b>{int(dc.get(k,0)):,}</b> {k.lower().replace('_',' ')}"
                                        for k in ("CONCORDANT", "MUTUALLY_EXCLUSIVE", "INDEPENDENT"))
                           + "</p>")
            conc_img = mod_mod_concordance_png(mm_view, top_n=args.max_snp_figs)
            fig_html = clickable_image_html(
                conc_img, "Co-localized modification 2x2 concordance panels",
                caption="Top pairs by FDR. Each 2x2 is observed count (expected under independence); "
                        "red = enriched. A concordant pair reddens the diagonal (both-modified / both-unmodified).",
                figure_class="hero-figure") if conc_img else ""
            mm_body = header + fig_html + df_to_html(mm_view, max_rows=args.top_genes)
        else:
            mm_body = "<p class='muted'>No co-localized modification pairs passed the coverage thresholds.</p>"
        sec_modmod = section(
            "Single-Molecule cis Co-Localized Modifications",
            mm_body,
            intro="For every pair of nearby modification sites seen on the SAME reads, does the state of one predict "
                  "the state of the other? The 2x2 is tested conditional on reads covering both sites, and enrichment "
                  "is judged against the expected co-occurrence given each site's own marginal stoichiometry — so a "
                  "concordant call means the two are modified together MORE than their individual rates alone would "
                  "predict, not merely that both are often modified. effect_abs_delta_mod_frac is "
                  "|P(B modified | A modified) - P(B modified | A unmodified)|; jaccard_both is co-modified reads / "
                  "reads modified at either.",
            definitions=definitions_html(column_definitions(list(mod_mod_df.columns)), summary="Column definitions") if not mod_mod_df.empty else "",
        )

    sec_hap_blocks = section(
        "Association of multiple cis SNPs into Haplotypes",
        df_to_html(hap_blocks_df_view, max_rows=args.top_genes) if not hap_blocks_df_view.empty else "<p class='muted'>No haplotype blocks available.</p>",
        intro="Local read-backed SNP blocks (co-occurring cis SNPs phased onto the same reads) retained for haplotype association testing.",
        definitions=definitions_html(column_definitions(list(hap_blocks_df_view.columns)), summary="Column definitions") if not hap_blocks_df_view.empty else "",
    )

    hap_sections = []
    if not hap_tx_df_view.empty:
        hap_sections.append(
            subsection(
                "Haplotype to Fragmentform Usage",
                df_to_html(hap_tx_df_view, max_rows=args.top_genes) + snp_galleries.get("hap_tx", ""),
                definitions=definitions_html(column_definitions(list(hap_tx_df_view.columns)), summary="Column definitions"),
            )
        )
    if not hap_mod_df_view.empty:
        hap_sections.append(
            subsection(
                "Haplotype to Modification Stoichiometry",
                df_to_html(hap_mod_df_view, max_rows=args.top_genes) + snp_galleries.get("hap_mod", ""),
                definitions=definitions_html(column_definitions(list(hap_mod_df_view.columns)), summary="Column definitions"),
            )
        )
    sec_hap_assoc = section(
        "Association between Haplotypes and Fragmentform Usage + Modification Stoichiometry",
        "".join(hap_sections) if hap_sections else "<p class='muted'>No haplotype associations available.</p>",
        intro="Associations between local haplotype blocks and fragmentform usage or modification stoichiometry.",
    )

    body = [s for s in [
        sec_overview, sec_top_ff, sec_sample_stats, sec_read_funnel, sec_gene_sites, sec_overlap,
        sec_novel, sec_splice, sec_apa, sec_diff, sec_class, sec_seq, sec_polya, sec_between,
        sec_snp_cand, sec_snp_ff, sec_snp_mod, sec_snp_mech, sec_modmod,
        sec_hap_blocks, sec_hap_assoc,
    ] if s]

    # --- header: modulator logo banner + a run manifest (command line + resolved inputs + config)
    logo_html = ""
    try:
        _logo_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "images",
                                  "modulator_banner.png")
        with open(_logo_path, "rb") as _lf:
            _logo_uri = "data:image/png;base64," + base64.b64encode(_lf.read()).decode("ascii")
        logo_html = f"<img class='report-logo' src='{_logo_uri}' alt='modulator' />"
    except Exception:
        logo_html = "<h1>modulator</h1>"

    manifest_html = ""
    _mtext = ""
    if getattr(args, "run_manifest", "") and os.path.exists(args.run_manifest):
        try:
            with open(args.run_manifest) as _mf:
                _mtext = _mf.read()
        except Exception:
            _mtext = ""
    if _mtext.strip():
        manifest_html = (
            "<section><details class='run-manifest' open>"
            "<summary>Run inputs &amp; parameters</summary>"
            f"<pre class='manifest-pre'>{html.escape(_mtext)}</pre>"
            "</details></section>"
        )

    html_doc = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>{html.escape(args.title)}</title>
  <style>
    :root {{
      --bg:#f7f9fb; --panel:#ffffff; --raised:#fdfefe;
      --ink:#131a22; --ink-soft:#3d4b5a; --muted:#65778a;
      --line:#e2e8ef; --line-soft:#eef2f6;
      --accent:#1f6feb; --accent-soft:#e8f0fd; --accent-ink:#0d4ba0;
      --pos:#1a7f52; --neg:#c0392b; --warn:#b7791f;
      --stripe:#f6f9fb; --shadow:0 1px 2px rgba(16,24,40,.04),0 1px 3px rgba(16,24,40,.06);
      --mono:ui-monospace,SFMono-Regular,"SF Mono",Menlo,Consolas,monospace;
    }}
    @media (prefers-color-scheme: dark) {{
      :root {{
        --bg:#0d1117; --panel:#151b23; --raised:#1a222c;
        --ink:#e6edf3; --ink-soft:#c2cdd8; --muted:#8b98a5;
        --line:#232c37; --line-soft:#1c242e;
        --accent:#58a6ff; --accent-soft:#132132; --accent-ink:#79b8ff;
        --pos:#3fb950; --neg:#f85149; --warn:#d29922;
        --stripe:#131a22; --shadow:0 1px 2px rgba(0,0,0,.3);
      }}
    }}
    * {{ box-sizing:border-box; }}
    body {{
      margin:0; background:var(--bg); color:var(--ink);
      font:15px/1.6 -apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,"Helvetica Neue",Arial,sans-serif;
      -webkit-font-smoothing:antialiased; text-rendering:optimizeLegibility;
    }}
    main {{ max-width:1240px; margin:0 auto; padding:28px 26px 80px; }}
    header {{ margin:8px 0 26px; padding-bottom:20px; border-bottom:1px solid var(--line); }}
    h1 {{
      font-size:clamp(26px,3.4vw,36px); line-height:1.1; margin:0 0 6px;
      font-weight:680; letter-spacing:-.022em; text-wrap:balance;
    }}
    /* Logo banner + "Report" wordmark, styled to vibe with the modulator logo (teal #17807f /
       gold #d4a017, rounded-geometric face). */
    header {{ text-align:center; }}
    .report-logo {{ display:block; width:100%; max-width:620px; height:auto; margin:2px auto 6px; }}
    h1.report-title {{
      font-family:"Century Gothic","Futura","Avenir Next","Trebuchet MS","Segoe UI",system-ui,sans-serif;
      font-size:clamp(30px,4.4vw,48px); font-weight:700; letter-spacing:.10em;
      color:#17807f; line-height:1; margin:0; text-transform:uppercase; text-align:center;
    }}
    h1.report-title::after {{
      content:""; display:block; width:64px; height:5px; border-radius:3px;
      margin:10px auto 0; background:#d4a017;
    }}
    /* run-manifest sits under the centered header but reads left-aligned like the section bodies */
    details.run-manifest {{ text-align:left; }}
    /* collapsible section headers (open/close), enlarged +12 over the base h2 size */
    details.report-section {{ margin:0; }}
    details.report-section > summary {{
      cursor:pointer; user-select:none; list-style:none;
      display:flex; align-items:center; gap:10px;
      margin:0 0 14px; padding-bottom:12px; border-bottom:1px solid var(--line-soft);
    }}
    details.report-section > summary::-webkit-details-marker {{ display:none; }}
    details.report-section > summary::before {{
      content:"\\25B8"; color:#d4a017; font-size:.8em; line-height:1;
    }}
    details.report-section[open] > summary::before {{ content:"\\25BE"; }}
    details.report-section > summary > h2 {{
      display:inline; margin:0; padding:0; border:0;
      font-size:28px; color:#17807f; font-weight:660;  /* green, matches the run-manifest header */
    }}
    /* nested classification tree: buckets (h2.bucket) and events (h3), collapsible with the same
       green summary + gold arrow, progressively smaller and indented. */
    details.report-subtree {{ margin:6px 0; }}
    details.report-subtree > summary {{
      cursor:pointer; user-select:none; list-style:none; display:flex; align-items:center; gap:9px;
      padding:6px 0; margin:0;
    }}
    details.report-subtree > summary::-webkit-details-marker {{ display:none; }}
    details.report-subtree > summary::before {{ content:"\\25B8"; color:#d4a017; font-size:.8em; }}
    details.report-subtree[open] > summary::before {{ content:"\\25BE"; }}
    details.report-subtree > summary > h2.bucket {{
      display:inline; margin:0; padding:0; border:0; font-size:23px; color:#17807f; font-weight:680;
    }}
    details.report-subtree > summary > h3 {{
      display:inline; margin:0; padding:0; font-size:17px; color:#17807f; font-weight:640;
    }}
    details.report-subtree > summary > h4 {{
      display:inline; margin:0; padding:0; font-size:14px; color:#1f9a8f; font-weight:600; letter-spacing:.01em;
    }}
    details.report-subtree > .section-body {{ padding-left:20px; border-left:2px solid var(--line-soft); margin-left:4px; }}
    details.run-manifest {{ margin:16px 0 4px; }}
    details.run-manifest > summary {{
      cursor:pointer; font-weight:660; color:#17807f; letter-spacing:-.01em;
      font-size:28px; list-style:none; user-select:none;
    }}
    details.run-manifest > summary::-webkit-details-marker {{ display:none; }}
    details.run-manifest > summary::before {{ content:"\\25B8 "; color:#d4a017; font-size:.8em; }}
    details.run-manifest[open] > summary::before {{ content:"\\25BE "; }}
    pre.manifest-pre {{
      max-height:600px; overflow:auto; margin:10px 0 0; padding:14px 16px;
      background:var(--panel); border:1px solid var(--line); border-radius:10px;
      font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace;
      font-size:12px; line-height:1.5; white-space:pre; color:var(--ink-soft);
    }}
    h2 {{
      font-size:19px; font-weight:660; letter-spacing:-.012em; margin:0 0 4px;
      padding-bottom:10px; border-bottom:1px solid var(--line-soft);
    }}
    h3 {{ font-size:14.5px; font-weight:640; margin:22px 0 8px; color:var(--ink-soft); }}
    a {{ color:var(--accent); text-decoration:none; }}
    a:hover {{ text-decoration:underline; }}
    section {{
      background:var(--panel); border:1px solid var(--line); border-radius:12px;
      padding:20px 22px; margin:0 0 18px; box-shadow:var(--shadow);
    }}
    .section-intro {{ color:var(--muted); font-size:14px; max-width:74ch; margin:10px 0 16px; }}
    .callout-warn {{
      border:1px solid #d8b24a; border-left:5px solid #d4a017; background:#fdf6e3;
      color:#5c4a12; border-radius:9px; padding:12px 15px; margin:0 0 16px; font-size:13.5px;
      line-height:1.5;
    }}
    @media (prefers-color-scheme: dark) {{
      .callout-warn {{ background:#2a2410; color:#e7d9a8; border-color:#5c4d1c; border-left-color:#d4a017; }}
    }}
    :root[data-theme="dark"] .callout-warn {{ background:#2a2410; color:#e7d9a8; border-color:#5c4d1c; border-left-color:#d4a017; }}
    :root[data-theme="light"] .callout-warn {{ background:#fdf6e3; color:#5c4a12; border-color:#d8b24a; border-left-color:#d4a017; }}
    .muted {{ color:var(--muted); font-style:normal; }}
    .subsection {{ margin-top:20px; padding-top:4px; }}
    /* ---- cards / stat tiles ---- */
    .cards, .overview {{
      display:grid; grid-template-columns:repeat(auto-fit,minmax(168px,1fr));
      gap:12px; margin:16px 0;
    }}
    .card {{
      background:var(--raised); border:1px solid var(--line); border-radius:10px; padding:13px 15px;
    }}
    .card .value, .card b {{
      display:block; font-size:23px; font-weight:670; letter-spacing:-.02em;
      font-variant-numeric:tabular-nums; line-height:1.15;
    }}
    .card .label, .card small {{ display:block; color:var(--muted); font-size:12px; margin-top:5px; }}
    /* ---- tables ---- */
    .table-wrap {{
      overflow-x:auto; border:1px solid var(--line); border-radius:10px; margin:12px 0; background:var(--panel);
    }}
    table {{ border-collapse:separate; border-spacing:0; width:100%; font-size:13px; }}
    table.datatable {{ font-variant-numeric:tabular-nums; }}
    th {{
      position:sticky; top:0; z-index:1; background:var(--accent-soft); color:var(--accent-ink);
      text-align:left; font-weight:640; font-size:11.5px; letter-spacing:.03em; text-transform:uppercase;
      padding:9px 12px; white-space:nowrap; border-bottom:1px solid var(--line);
    }}
    td {{
      padding:8px 12px; border-bottom:1px solid var(--line-soft); white-space:nowrap;
      font-family:var(--mono); font-size:12px; color:var(--ink-soft);
    }}
    tbody tr:nth-child(even) td {{ background:var(--stripe); }}
    tbody tr:hover td {{ background:var(--accent-soft); }}
    tbody tr:last-child td {{ border-bottom:none; }}
    /* ---- figures ---- */
    figure {{ margin:16px 0; }}
    img {{ max-width:100%; height:auto; border:1px solid var(--line); border-radius:9px; background:var(--panel); }}
    figcaption {{ color:var(--muted); font-size:12px; margin-top:7px; }}
    .gallery {{
      display:grid; grid-template-columns:repeat(auto-fit,minmax(290px,1fr)); gap:14px; margin:12px 0;
    }}
    .gallery figure {{ margin:0; }}
    /* ---- collapsible definitions ---- */
    details.definitions {{
      margin:14px 0; border:1px solid var(--line); border-radius:10px; background:var(--raised); overflow:hidden;
    }}
    details.definitions > summary {{
      cursor:pointer; padding:10px 15px; font-weight:600; font-size:13px; color:var(--accent);
      list-style:none; user-select:none;
    }}
    details.definitions > summary::-webkit-details-marker {{ display:none; }}
    details.definitions > summary::before {{ content:"▸ "; color:var(--muted); }}
    details.definitions[open] > summary::before {{ content:"▾ "; }}
    details.definitions[open] > summary {{ border-bottom:1px solid var(--line); }}
    details.definitions dl {{ margin:0; padding:12px 16px 15px; }}
    details.definitions dt {{ font-weight:640; font-size:12.5px; margin-top:11px; font-family:var(--mono); color:var(--ink); }}
    details.definitions dt:first-child {{ margin-top:0; }}
    details.definitions dd {{ margin:3px 0 0; color:var(--muted); font-size:13px; line-height:1.5; }}
    ul {{ padding-left:20px; }} li {{ margin:5px 0; }}
    li b {{ font-variant-numeric:tabular-nums; }}
    code {{ font-family:var(--mono); font-size:12.5px; background:var(--accent-soft); padding:1px 6px; border-radius:4px; }}
    footer {{ margin-top:36px; padding-top:16px; border-top:1px solid var(--line); color:var(--muted); font-size:12.5px; }}
    @media (max-width:720px) {{ main {{ padding:18px 14px 60px; }} section {{ padding:16px; }} }}
  </style>
</head>
<body>
  <main>
    <header>
      {logo_html}
      <h1 class="report-title">Report</h1>
    </header>
    {manifest_html}
    {''.join(body)}
  </main>
  <!-- Figures are plain <a class="image-link" target="_blank"> anchors: a left-click opens the
       full-size image file (written next to this report by externalize_data_uris) in a NEW TAB.
       No in-page lightbox/JS interception -- data: URIs cannot be opened as a top-level tab in
       modern browsers, so the previous JS overlay is intentionally gone. -->
</body>
</html>
"""

    os.makedirs(os.path.dirname(args.out_html) or ".", exist_ok=True)
    html_doc, n_imgs, rel_dir = externalize_data_uris(html_doc, args.out_html)
    with open(args.out_html, "w") as out:
        out.write(html_doc)
    print(f"[report] externalized {n_imgs} image(s) into {rel_dir}/", file=sys.stderr)


if __name__ == "__main__":
    main()
