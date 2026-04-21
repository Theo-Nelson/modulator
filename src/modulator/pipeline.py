from __future__ import annotations

import glob
import os
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from modulator.runtime import (
    as_bool,
    ensure_mpl_config_dir,
    ensure_parent,
    find_project_root,
    format_command,
    is_set,
    require_tools,
    resolve_path,
    run_command,
    run_parallel,
)


STAGE_ORDER = [
    "assemble",
    "read_stats",
    "multigene_filter",
    "modkit_zn",
    "modkit_zt",
    "aggregate_zn",
    "aggregate_zt",
    "test_diffs",
    "genotype",
    "report",
]


@dataclass
class PipelinePaths:
    root: Path
    prefix: str

    @property
    def results(self) -> Path:
        return self.root / "results"

    @property
    def assemble(self) -> Path:
        return self.results / "assemble"

    @property
    def aggregate_zn(self) -> Path:
        return self.results / "aggregate_zn"

    @property
    def aggregate_zt(self) -> Path:
        return self.results / "aggregate_zt"

    @property
    def modkit_zn(self) -> Path:
        return self.results / "modkit_zn"

    @property
    def modkit_zt(self) -> Path:
        return self.results / "modkit_zt"

    @property
    def test_diffs(self) -> Path:
        return self.results / "test_diffs"

    @property
    def genotype(self) -> Path:
        return self.results / "genotype"

    @property
    def report(self) -> Path:
        return self.results / "report"

    @property
    def out_gtf(self) -> Path:
        return self.assemble / f"{self.prefix}.gtf"

    @property
    def classification_summary(self) -> Path:
        return self.assemble / f"{self.prefix}_classification_summary.tsv"

    @property
    def metrics(self) -> Path:
        return self.assemble / f"{self.prefix}_metrics.tsv"

    @property
    def tx_counts(self) -> Path:
        return self.assemble / f"{self.prefix}_tx_counts.tsv"

    @property
    def pca_png(self) -> Path:
        return self.assemble / f"{self.prefix}_tx_counts.pca.png"

    @property
    def sample_stats(self) -> Path:
        return self.assemble / f"{self.prefix}_per_sample_stats.tsv"

    @property
    def per_sample_read_stats(self) -> Path:
        return self.assemble / f"{self.prefix}_per_sample_read_stats.tsv"

    @property
    def tx_assigned_read_lengths(self) -> Path:
        return self.assemble / f"{self.prefix}_tx_assigned_read_lengths.tsv"

    @property
    def partition_map(self) -> Path:
        return self.assemble / f"{self.prefix}_partition_map.tsv"

    @property
    def multigene_scrap_tx_counts(self) -> Path:
        return self.assemble / f"{self.prefix}_multigene_scrap_tx_counts.tsv"

    @property
    def zt_tagged_dir(self) -> Path:
        return self.assemble / "zt_tagged"

    @property
    def zt_filtered_dir(self) -> Path:
        return self.assemble / "zt_filtered"

    @property
    def zt_scrap_dir(self) -> Path:
        return self.assemble / "zt_scrap"

    @property
    def zn_filtered_long(self) -> Path:
        return self.aggregate_zn / f"{self.prefix}_FILTERED_sites_long.tsv"

    @property
    def zt_filtered_long(self) -> Path:
        return self.aggregate_zt / f"{self.prefix}_FILTERED_long.tsv"

    @property
    def zn_diff_results(self) -> Path:
        return self.test_diffs / f"{self.prefix}__ZN_site_diff_results.tsv"

    @property
    def zn_diff_figs(self) -> Path:
        return self.test_diffs / f"{self.prefix}__figs"

    @property
    def report_html(self) -> Path:
        return self.report / f"{self.prefix}_report.html"

    @property
    def geno_read_assignments(self) -> Path:
        return self.genotype / f"{self.prefix}_read_assignments.tsv"

    @property
    def geno_candidate_snps(self) -> Path:
        return self.genotype / f"{self.prefix}_candidate_snps.tsv"

    @property
    def geno_molecule_snps(self) -> Path:
        return self.genotype / f"{self.prefix}_molecule_snps.tsv"

    @property
    def geno_candidate_mod_sites(self) -> Path:
        return self.genotype / f"{self.prefix}_candidate_mod_sites.tsv"

    @property
    def geno_candidate_mod_bed(self) -> Path:
        return self.genotype / f"{self.prefix}_candidate_mod_sites.bed"

    @property
    def geno_molecule_mod_calls(self) -> Path:
        return self.genotype / f"{self.prefix}_molecule_mod_calls.tsv"

    @property
    def geno_snp_tx(self) -> Path:
        return self.genotype / f"{self.prefix}_snp_transcript_assoc.tsv"

    @property
    def geno_snp_mod(self) -> Path:
        return self.genotype / f"{self.prefix}_snp_mod_assoc.tsv"

    @property
    def geno_joint(self) -> Path:
        return self.genotype / f"{self.prefix}_snp_tx_mod_dependency.tsv"

    @property
    def geno_hap_blocks(self) -> Path:
        return self.genotype / f"{self.prefix}_haplotype_blocks.tsv"

    @property
    def geno_molecule_haps(self) -> Path:
        return self.genotype / f"{self.prefix}_molecule_haplotypes.tsv"

    @property
    def geno_hap_tx(self) -> Path:
        return self.genotype / f"{self.prefix}_haplotype_transcript_assoc.tsv"

    @property
    def geno_hap_mod(self) -> Path:
        return self.genotype / f"{self.prefix}_haplotype_mod_assoc.tsv"

    def zt_tagged_bam(self, sample: str) -> Path:
        return self.zt_tagged_dir / f"{sample}.zt_tagged.bam"

    def clean_bam(self, sample: str) -> Path:
        return self.zt_filtered_dir / f"{sample}.zt_tagged.clean.bam"

    def scrap_bam(self, sample: str) -> Path:
        return self.zt_scrap_dir / f"{sample}.zt_tagged.multigene_scrap.bam"

    def multigene_summary(self, sample: str) -> Path:
        return self.zt_scrap_dir / f"{sample}.multigene_filter_summary.tsv"

    def multigene_removed(self, sample: str) -> Path:
        return self.zt_scrap_dir / f"{sample}.multigene_removed_reads.tsv"

    def scrap_tx_counts(self, sample: str) -> Path:
        return self.zt_scrap_dir / f"{sample}.multigene_scrap_tx_counts.tsv"

    def modkit_dir(self, which: str, sample: str) -> Path:
        return (self.modkit_zn if which == "zn" else self.modkit_zt) / sample


class ModulatorPipeline:
    def __init__(self, config: dict[str, Any], *, workdir: str | Path, jobs: int = 1, verbose: bool = True):
        self.root = find_project_root(workdir)
        ensure_mpl_config_dir(self.root)
        self.config = config
        self.jobs = max(1, int(jobs))
        self.verbose = verbose
        self.prefix = str(config.get("prefix", "modulator_run"))
        self.paths = PipelinePaths(self.root, self.prefix)
        self.samples = self._discover_samples()
        self.reference_fa = self._resolve_reference("reference_fa", ("reference", "fasta"))
        self.reference_gtf = self._resolve_reference("reference_gtf", ("reference", "gtf"))
        self.top_threads = int(config.get("threads", 1))
        self._validate_config()

    def _discover_samples(self) -> list[str]:
        bams_dir = self.bams_dir
        bam_glob = self.bam_glob
        found = sorted(glob.glob(str(bams_dir / bam_glob)))
        samples = [Path(path).stem for path in found]
        if not samples:
            raise FileNotFoundError(f"No BAMs matched {bam_glob!r} under {bams_dir}")
        return samples

    @property
    def bams_dir(self) -> Path:
        return resolve_path(self.root, self.config.get("bams_dir", "resources/test_bams/ALCAM_NHSL1_SERAC1_MXD1_RIOK3_reads"))

    @property
    def bam_glob(self) -> str:
        return str(self.config.get("bam_glob", "*.bam"))

    def _resolve_reference(self, flat_key: str, nested_key: tuple[str, str]) -> Path | None:
        if is_set(self.config.get(flat_key)):
            return resolve_path(self.root, self.config.get(flat_key))
        parent, child = nested_key
        nested = self.config.get(parent, {})
        if isinstance(nested, dict) and is_set(nested.get(child)):
            return resolve_path(self.root, nested.get(child))
        return None

    def _validate_config(self) -> None:
        toggles = self.config.get("toggles", {})
        downstream_needs_sample_bams = (
            as_bool(self.config.get("multigene_filter", {}).get("enable", True), True)
            or as_bool(toggles.get("enable_zn_pileup", True), True)
            or as_bool(toggles.get("enable_zt_pileup", True), True)
            or as_bool(self.config.get("genotype", {}).get("enable", False), False)
        )
        assembler_cfg = self.config.setdefault("assembler", {})
        if downstream_needs_sample_bams and not as_bool(assembler_cfg.get("write_zt_tagged_sample_bams", True), True):
            assembler_cfg["write_zt_tagged_sample_bams"] = True
            if self.verbose:
                print("[modulator] enabled assembler.write_zt_tagged_sample_bams because downstream stages require per-sample tagged BAMs.", flush=True)

    def _require_existing_file(self, path: Path | None, label: str) -> Path:
        if path is None:
            raise ValueError(f"{label} is required for this stage but was not configured.")
        if not path.exists():
            raise FileNotFoundError(f"{label} does not exist: {path}")
        if not path.is_file():
            raise FileNotFoundError(f"{label} is not a file: {path}")
        return path

    def _require_existing_dir(self, path: Path | None, label: str) -> Path:
        if path is None:
            raise ValueError(f"{label} is required for this stage but was not configured.")
        if not path.exists():
            raise FileNotFoundError(f"{label} does not exist: {path}")
        if not path.is_dir():
            raise FileNotFoundError(f"{label} is not a directory: {path}")
        return path

    def _require_reference_fa(self) -> Path:
        return self._require_existing_file(self.reference_fa, "reference FASTA")

    def _require_reference_gtf(self) -> Path:
        return self._require_existing_file(self.reference_gtf, "reference GTF")

    def _require_modkit_outputs(self, which: str) -> Path:
        base = self.paths.modkit_zn if which == "zn" else self.paths.modkit_zt
        self._require_existing_dir(base, f"modkit {which.upper()} output directory")
        has_partitions = any(base.rglob("*.bed")) or any(base.rglob("*.bed.gz"))
        if not has_partitions:
            raise FileNotFoundError(
                f"No partition BED outputs were found under {base}. Run the modkit {which.upper()} stage first."
            )
        return base

    def script_path(self, name: str) -> Path:
        path = self.root / "workflow" / "scripts" / name
        if not path.exists():
            raise FileNotFoundError(f"Required workflow script is missing: {path}")
        return path

    def run_python_script(self, script_name: str, args: list[str], *, label: str) -> None:
        cmd = [sys.executable, str(self.script_path(script_name)), *args]
        run_command(cmd, cwd=self.root, label=label, verbose=self.verbose)

    def run(self, stages: list[str] | None = None) -> None:
        selected = STAGE_ORDER if not stages else [stage for stage in STAGE_ORDER if stage in stages]
        for stage in selected:
            if self.verbose:
                print(f"[modulator] stage: {stage}", flush=True)
            getattr(self, f"stage_{stage}")()

    def stage_assemble(self) -> None:
        cfg = self.config.get("assembler", {})
        reference_gtf = self._require_reference_gtf()
        args = [
            "--dir", str(self.bams_dir),
            "--glob", self.bam_glob,
            "--gtf", str(reference_gtf),
            "--out-gtf", str(self.paths.out_gtf),
            "--out-prefix", str(self.paths.assemble / self.prefix),
            "--threads", str(self.top_threads),
            "--min-mapq", str(cfg.get("min_mapq", 10)),
            "--min-introns-read", str(cfg.get("min_introns_read", 1)),
            "--require-softclip3p", str(cfg.get("require_softclip3p", 0)),
            "--apa-window", str(cfg.get("apa_window", 20)),
            "--min-reads", str(cfg.get("min_reads", 40)),
            "--min-frac", str(cfg.get("min_frac", 0.00)),
            "--min-introns", str(cfg.get("min_introns", 1)),
            "--min-polya-length", str(cfg.get("min_polya_length", 12)),
            "--min-polya-purity", str(cfg.get("min_polya_purity", 0.5)),
            "--polya-support-frac", str(cfg.get("polya_support_frac", 0.5)),
            "--tes-match-tol", str(cfg.get("tes_match_tol", 25)),
            "--exact-tes-tol", str(cfg.get("exact_tes_tol", 10)),
            "--assignment-mode", str(cfg.get("assignment_mode", "support_first")),
            "--zn-mode", str(cfg.get("zn_mode", "metagene_colored")),
            "--min-distal-anchor-reads", str(cfg.get("min_distal_anchor_reads", 2)),
            "--min-distal-anchor-frac", str(cfg.get("min_distal_anchor_frac", 0.05)),
            "--min-exact-canonical-reads", str(cfg.get("min_exact_canonical_reads", 1)),
            "--min-reads-per-sample-for-mod", str(cfg.get("min_reads_per_sample_for_mod", 5)),
            "--min-total-reads-for-mod", str(cfg.get("min_total_reads_for_mod", 20)),
        ]
        if as_bool(cfg.get("primary_only", True), True):
            args.append("--primary-only")
        if cfg.get("tes_window") is not None:
            args.extend(["--tes-window", str(cfg.get("tes_window"))])
        if as_bool(cfg.get("write_zt_bams", False), False):
            args.append("--write-zt-bams")
        if as_bool(cfg.get("write_zt_tagged_sample_bams", True), True):
            args.append("--write-zt-tagged-sample-bams")
        if as_bool(cfg.get("emit_modkit_manifest", False), False):
            args.append("--emit-modkit-manifest")
        if cfg.get("status_every") is not None:
            args.extend(["--status-every", str(cfg.get("status_every"))])
        self.run_python_script("assemble_transcripts.py", args, label="assemble_transcripts")

    def stage_read_stats(self) -> None:
        cfg = self.config.get("assembler", {})
        self._require_existing_dir(self.paths.zt_tagged_dir, "ZT-tagged sample BAM directory")
        args = [
            "--bams-dir", str(self.bams_dir),
            "--bam-glob", self.bam_glob,
            "--zt-tagged-dir", str(self.paths.zt_tagged_dir),
            "--out", str(self.paths.per_sample_read_stats),
            "--min-mapq", str(cfg.get("min_mapq", 10)),
            "--min-introns-read", str(cfg.get("min_introns_read", 1)),
            "--require-softclip3p", str(cfg.get("require_softclip3p", 0)),
        ]
        if as_bool(cfg.get("primary_only", True), True):
            args.append("--primary-only")
        self.run_python_script("per_sample_read_stats.py", args, label="per_sample_read_stats")

    def stage_multigene_filter(self) -> None:
        cfg = self.config.get("multigene_filter", {})
        if not as_bool(cfg.get("enable", True), True):
            return
        self._require_existing_file(self.paths.out_gtf, "assembled GTF")

        tasks = []
        for sample in self.samples:
            input_bam = self.paths.zt_tagged_bam(sample)
            self._require_existing_file(input_bam, f"ZT-tagged BAM for sample {sample}")
            output_clean = self.paths.clean_bam(sample)
            output_scrap = self.paths.scrap_bam(sample)
            output_summary = self.paths.multigene_summary(sample)
            output_removed = self.paths.multigene_removed(sample)
            output_counts = self.paths.scrap_tx_counts(sample)
            for path in [output_clean, output_scrap, output_summary, output_removed, output_counts]:
                ensure_parent(path)

            args = [
                "--bam", str(input_bam),
                "--gtf", str(self.paths.out_gtf),
                "--sample", sample,
                "--out-clean-bam", str(output_clean),
                "--out-scrap-bam", str(output_scrap),
                "--out-summary-tsv", str(output_summary),
                "--out-removed-tsv", str(output_removed),
                "--out-scrap-tx-counts-tsv", str(output_counts),
                "--zero-gene-action", str(cfg.get("zero_gene_action", "keep")),
                "--mode", str(cfg.get("mode", "resolve")),
                "--same-strand-only",
            ]
            tasks.append((
                f"multigene_filter[{sample}]",
                lambda sample_name=sample, sample_args=args: self.run_python_script(
                    "filter_multigene_reads_from_zt_bam.py",
                    sample_args,
                    label=f"filter_multigene_reads_from_zt_bam[{sample_name}]",
                ),
            ))
        run_parallel(tasks, jobs=self.jobs)

        counts = [str(self.paths.scrap_tx_counts(sample)) for sample in self.samples]
        args = ["--counts", *counts, "--out", str(self.paths.multigene_scrap_tx_counts)]
        self.run_python_script("aggregate_scrap_tx_counts.py", args, label="aggregate_scrap_tx_counts")

    def _modkit_input_bam(self, sample: str) -> Path:
        if as_bool(self.config.get("multigene_filter", {}).get("enable", True), True):
            return self.paths.clean_bam(sample)
        return self.paths.zt_tagged_bam(sample)

    def _format_common_modkit_flags(self, common: dict[str, Any], *, sample: str, which: str, threads: int) -> list[str]:
        flags: list[str] = []
        log_tmpl = common.get("log_file_template")
        if is_set(log_tmpl):
            log_path = resolve_path(self.root, str(log_tmpl).format(sample=sample, which=which))
            ensure_parent(log_path)
            flags += ["--log-filepath", str(log_path)]
        if is_set(common.get("region")):
            flags += ["--region", str(common["region"])]
        if common.get("max_depth") is not None:
            flags += ["--max-depth", str(common["max_depth"])]
        if is_set(common.get("include_bed")):
            flags += ["--include-bed", str(resolve_path(self.root, common["include_bed"]))]
        if as_bool(common.get("include_unmapped"), False):
            flags += ["--include-unmapped"]
        if is_set(common.get("edge_filter")):
            flags += ["--edge-filter", str(common["edge_filter"])]
        if as_bool(common.get("invert_edge_filter"), False):
            flags += ["--invert-edge-filter"]
        flags += ["-t", str(threads)]
        if is_set(common.get("interval_size")):
            flags += ["--interval-size", str(common["interval_size"])]
        if is_set(common.get("queue_size")):
            flags += ["--queue-size", str(common["queue_size"])]
        if is_set(common.get("chunk_size")):
            flags += ["--chunk-size", str(common["chunk_size"])]
        if is_set(common.get("num_reads")):
            flags += ["--num-reads", str(common["num_reads"])]
        if is_set(common.get("sampling_frac")):
            flags += ["--sampling-frac", str(common["sampling_frac"])]
        seed_val = common.get("seed")
        if seed_val is not None and not (isinstance(seed_val, str) and seed_val.strip().lower() in {"none", "null", ""}):
            flags += ["--seed", str(seed_val)]
        if is_set(common.get("sample_region")):
            flags += ["--sample-region", str(common["sample_region"])]
        if is_set(common.get("sampling_interval_size")):
            flags += ["--sampling-interval-size", str(common["sampling_interval_size"])]
        if as_bool(common.get("no_filtering"), False):
            flags += ["--no-filtering"]
        if common.get("filter_percentile") is not None:
            flags += ["--filter-percentile", str(common["filter_percentile"])]
        for item in common.get("filter_thresholds") or []:
            if is_set(item):
                flags += ["--filter-threshold", str(item)]
        for item in common.get("mod_thresholds") or []:
            if is_set(item) and ":" in str(item):
                flags += ["--mod-threshold", str(item)]
        for item in common.get("ignore") or []:
            if is_set(item):
                flags += ["--ignore", str(item)]
        if as_bool(common.get("force_allow_implicit"), False):
            flags += ["--force-allow-implicit"]
        for item in common.get("motif") or []:
            if is_set(item) and ":" in str(item):
                motif, off = str(item).split(":", 1)
                flags += ["--motif", motif, str(off)]
        if as_bool(common.get("cpg"), False):
            flags += ["--cpg"]
        if as_bool(common.get("ref_mask"), False):
            flags += ["--mask"]
        if as_bool(common.get("combine_mods"), False):
            flags += ["--combine-mods"]
        if as_bool(common.get("combine_strands"), False):
            flags += ["--combine-strands"]
        if as_bool(common.get("only_tabs"), False):
            flags += ["--only-tabs"]
        if as_bool(common.get("mixed_delim"), False):
            flags += ["--mixed-delim"]
        if as_bool(common.get("bedgraph"), False):
            flags += ["--bedgraph"]
        if as_bool(common.get("header"), False):
            flags += ["--header"]
        if is_set(common.get("prefix")):
            flags += ["--prefix", str(common["prefix"])]
        if as_bool(common.get("suppress_progress", True), True):
            flags += ["--suppress-progress"]
        return flags

    def _run_modkit_pileup(self, *, sample: str, which: str, partition_tag: str) -> None:
        require_tools(["modkit", "bgzip", "tabix"])
        reference_fa = self._require_reference_fa()
        modkit_cfg = self.config.get("modkit", {})
        common = modkit_cfg.get("common", {})
        output_dir = self.paths.modkit_dir(which, sample)
        output_dir.mkdir(parents=True, exist_ok=True)
        input_bam = self._modkit_input_bam(sample)
        self._require_existing_file(input_bam, f"modkit input BAM for sample {sample}")
        threads = int(common.get("threads", self.top_threads if self.top_threads > 0 else 4))
        flags = self._format_common_modkit_flags(common, sample=sample, which=f"modkit_{which}", threads=threads)
        cmd = [
            "modkit",
            "pileup",
            str(input_bam),
            str(output_dir),
            "--ref",
            str(reference_fa),
            *flags,
            "--partition-tag",
            partition_tag,
        ]
        run_command(cmd, cwd=self.root, label=f"modkit_pileup_{which}[{sample}]", verbose=self.verbose)

        for bed_path in sorted(output_dir.rglob("*.bed")):
            if bed_path.name == "ungrouped.bed":
                continue
            run_command(["bgzip", "-f", "-@", str(threads), str(bed_path)], cwd=self.root, label=f"bgzip[{bed_path.name}]", verbose=self.verbose)
            run_command(["tabix", "-f", "-p", "bed", str(bed_path) + ".gz"], cwd=self.root, label=f"tabix[{bed_path.name}.gz]", verbose=self.verbose)

    def stage_modkit_zn(self) -> None:
        if not as_bool(self.config.get("toggles", {}).get("enable_zn_pileup", True), True):
            return
        zn_cfg = self.config.get("modkit", {}).get("zn", {})
        partition_tag = str(zn_cfg.get("partition_tag", "ZN"))
        tasks = [
            (
                f"modkit_zn[{sample}]",
                lambda sample=sample: self._run_modkit_pileup(sample=sample, which="zn", partition_tag=partition_tag),
            )
            for sample in self.samples
        ]
        run_parallel(tasks, jobs=self.jobs)

    def stage_modkit_zt(self) -> None:
        if not as_bool(self.config.get("toggles", {}).get("enable_zt_pileup", True), True):
            return
        zt_cfg = self.config.get("modkit", {}).get("zt", {})
        partition_tag = str(zt_cfg.get("partition_tag", "ZT"))
        tasks = [
            (
                f"modkit_zt[{sample}]",
                lambda sample=sample: self._run_modkit_pileup(sample=sample, which="zt", partition_tag=partition_tag),
            )
            for sample in self.samples
        ]
        run_parallel(tasks, jobs=self.jobs)

    def stage_aggregate_zn(self) -> None:
        if not as_bool(self.config.get("toggles", {}).get("enable_zn_aggregate", True), True):
            return
        agg_cfg = self.config.get("aggregation", {}).get("zn", {})
        filters_cfg = self.config.get("filters", {})
        out_prefix = self.paths.aggregate_zn / self.prefix
        out_prefix.parent.mkdir(parents=True, exist_ok=True)
        self._require_modkit_outputs("zn")
        self._require_existing_file(self.paths.out_gtf, "assembled GTF")
        tmpdir = (
            self.config.get("aggregation_tmpdir")
            or self.config.get("aggregation", {}).get("tmpdir")
            or os.environ.get("TMPDIR")
            or str((self.paths.results / "tmp").resolve())
        )
        tmpdir_path = resolve_path(self.root, tmpdir)
        if tmpdir_path is None:
            raise ValueError("Could not resolve aggregation tmpdir.")
        tmpdir_path.mkdir(parents=True, exist_ok=True)
        args = [
            "--modkit-dir", str(self.paths.modkit_zn),
            "--gtf", str(self.paths.out_gtf),
            "--out-prefix", str(out_prefix),
            "--min-cov", str(self.config.get("min_cov", 5)),
            "--tmpdir", str(tmpdir_path),
            "--chunk-lines", str(int(self.config.get("aggregation_chunk_lines") or self.config.get("aggregation", {}).get("chunk_lines", 2000000))),
            "--count-diff-factor", str(float(agg_cfg.get("count_diff_factor", filters_cfg.get("count_diff_factor", 3)))),
            "--mod-fail-margin", str(int(agg_cfg.get("mod_fail_margin", filters_cfg.get("mod_fail_margin", 1)))),
            "--verbose",
        ]
        if as_bool(agg_cfg.get("filter_enable", filters_cfg.get("enable_site_filter", True)), True):
            args.append("--filter-enable")
        args.append("--emit-raw" if as_bool(agg_cfg.get("emit_raw", self.config.get("aggregate_outputs", {}).get("emit_raw", True)), True) else "--no-emit-raw")
        args.append("--emit-filtered" if as_bool(agg_cfg.get("emit_filtered", self.config.get("aggregate_outputs", {}).get("emit_filtered", True)), True) else "--no-emit-filtered")
        args.append("--write-long" if as_bool(agg_cfg.get("write_long", self.config.get("aggregate_outputs", {}).get("write_long", True)), True) else "--no-write-long")
        args.append("--write-pivots" if as_bool(agg_cfg.get("write_pivots", self.config.get("aggregate_outputs", {}).get("write_pivots", True)), True) else "--no-write-pivots")
        args.append("--write-raw-per-gene" if as_bool(agg_cfg.get("write_raw_per_gene", self.config.get("aggregate_outputs", {}).get("write_raw_per_gene", False)), False) else "--no-write-raw-per-gene")
        args.append("--write-filtered-per-gene" if as_bool(agg_cfg.get("write_filtered_per_gene", self.config.get("aggregate_outputs", {}).get("write_filtered_per_gene", True)), True) else "--no-write-filtered-per-gene")
        self.run_python_script("aggregate_by_gene.py", args, label="aggregate_by_gene")

    def stage_aggregate_zt(self) -> None:
        if not as_bool(self.config.get("toggles", {}).get("enable_zt_aggregate", True), True):
            return
        agg_cfg = self.config.get("aggregation", {}).get("zt", {})
        filters_cfg = self.config.get("filters", {})
        out_prefix = self.paths.aggregate_zt / self.prefix
        out_prefix.parent.mkdir(parents=True, exist_ok=True)
        self._require_modkit_outputs("zt")
        self._require_existing_file(self.paths.classification_summary, "classification summary TSV")
        args = [
            "--modkit-dir", str(self.paths.modkit_zt),
            "--summary-tsv", str(self.paths.classification_summary),
            "--out-prefix", str(out_prefix),
            "--min-cov", str(self.config.get("min_cov", 5)),
            "--count-diff-factor", str(float(agg_cfg.get("count_diff_factor", filters_cfg.get("count_diff_factor", 3)))),
            "--mod-fail-margin", str(int(agg_cfg.get("mod_fail_margin", filters_cfg.get("mod_fail_margin", 1)))),
            "--debug-summary",
            "--verbose",
        ]
        if as_bool(agg_cfg.get("filter_enable", filters_cfg.get("enable_site_filter", True)), True):
            args.append("--filter-enable")
        args.append("--emit-raw" if as_bool(agg_cfg.get("emit_raw", self.config.get("aggregate_outputs", {}).get("emit_raw", True)), True) else "--no-emit-raw")
        args.append("--emit-filtered" if as_bool(agg_cfg.get("emit_filtered", self.config.get("aggregate_outputs", {}).get("emit_filtered", True)), True) else "--no-emit-filtered")
        args.append("--write-long" if as_bool(agg_cfg.get("write_long", self.config.get("aggregate_outputs", {}).get("write_long", True)), True) else "--no-write-long")
        args.append("--write-pivots" if as_bool(agg_cfg.get("write_pivots", self.config.get("aggregate_outputs", {}).get("write_pivots", True)), True) else "--no-write-pivots")
        self.run_python_script("aggregate_by_transcript.py", args, label="aggregate_by_transcript")

    def stage_test_diffs(self) -> None:
        if not as_bool(self.config.get("toggles", {}).get("enable_test_diffs", True), True):
            return
        if not self.paths.zn_filtered_long.exists() or self.paths.zn_filtered_long.stat().st_size == 0:
            if self.verbose:
                print(
                    f"[modulator] skipping test_diffs because the ZN filtered long table is missing or empty: {self.paths.zn_filtered_long}",
                    flush=True,
                )
            return
        args = [
            "--in-tsv", str(self.paths.zn_filtered_long),
            "--out-prefix", str(self.paths.test_diffs / self.prefix),
            "--min-cov", str(self.config.get("test_diffs", {}).get("min_cov", self.config.get("min_cov_test", 20))),
            "--topk", str(self.config.get("test_diffs", {}).get("topk", self.config.get("topk", 10))),
            "--verbose",
        ]
        td = self.config.get("test_diffs", {})
        if is_set(td.get("test")):
            args.extend(["--test", str(td["test"])])
        if td.get("pseudocount") is not None:
            args.extend(["--pseudocount", str(td["pseudocount"])])
        if is_set(td.get("alternative")):
            args.extend(["--alternative", str(td["alternative"])])
        for gene in td.get("gene_filter") or []:
            args.extend(["--gene-filter", str(gene)])
        for mod in td.get("mod_filter") or []:
            args.extend(["--mod-filter", str(mod)])
        self.run_python_script("test_stoichiometry_diffs.py", args, label="test_stoichiometry_diffs")

    def stage_genotype(self) -> None:
        geno = self.config.get("genotype", {})
        if not as_bool(geno.get("enable", False), False):
            return
        self._require_reference_fa()
        geno_jobs = max(1, min(len(self.samples), int(geno.get("jobs", 2))))
        sample_bams = [
            str(self._require_existing_file(self._modkit_input_bam(sample), f"genotype input BAM for sample {sample}"))
            for sample in self.samples
        ]
        self._require_existing_file(self.paths.classification_summary, "classification summary TSV")
        self.paths.genotype.mkdir(parents=True, exist_ok=True)

        self.run_python_script(
            "build_read_assignment_table.py",
            [
                "--bams", *sample_bams,
                "--summary-tsv", str(self.paths.classification_summary),
                "--out-tsv", str(self.paths.geno_read_assignments),
                "--jobs", str(geno_jobs),
                "--primary-only",
                "--verbose",
            ],
            label="build_read_assignment_table",
        )
        self.run_python_script(
            "discover_candidate_snps.py",
            [
                "--bams", *sample_bams,
                "--reference-fa", str(self.reference_fa),
                "--gtf", str(self.paths.out_gtf),
                "--out-tsv", str(self.paths.geno_candidate_snps),
                "--min-alt-reads", str(int(geno.get("min_alt_reads", 4))),
                "--min-total-cov", str(int(geno.get("min_total_cov", 8))),
                "--min-alt-frac", str(float(geno.get("min_alt_frac", 0.10))),
                "--max-alt-frac", str(float(geno.get("max_alt_frac", 0.90))),
                "--min-baseq", str(int(geno.get("min_baseq", 20))),
                "--min-mapq", str(int(geno.get("min_mapq", self.config.get("assembler", {}).get("min_mapq", 10)))),
                "--jobs", str(geno_jobs),
                "--primary-only",
                "--verbose",
            ],
            label="discover_candidate_snps",
        )
        self.run_python_script(
            "build_molecule_snp_table.py",
            [
                "--bams", *sample_bams,
                "--candidate-snps", str(self.paths.geno_candidate_snps),
                "--out-tsv", str(self.paths.geno_molecule_snps),
                "--min-baseq", str(int(geno.get("min_baseq", 20))),
                "--min-mapq", str(int(geno.get("min_mapq", self.config.get("assembler", {}).get("min_mapq", 10)))),
                "--jobs", str(geno_jobs),
                "--primary-only",
                "--verbose",
            ],
            label="build_molecule_snp_table",
        )
        zn_long = str(self.paths.zn_filtered_long) if self.paths.zn_filtered_long.exists() else ""
        zt_long = str(self.paths.zt_filtered_long) if self.paths.zt_filtered_long.exists() else ""
        self.run_python_script(
            "build_candidate_mod_sites.py",
            [
                "--zn-long", zn_long,
                "--zt-long", zt_long,
                "--out-tsv", str(self.paths.geno_candidate_mod_sites),
                "--out-bed", str(self.paths.geno_candidate_mod_bed),
                "--min-total-cov", str(int(geno.get("min_mod_site_cov", 1))),
            ],
            label="build_candidate_mod_sites",
        )
        self.run_python_script(
            "build_molecule_mod_table.py",
            [
                "--bams", *sample_bams,
                "--candidate-sites-tsv", str(self.paths.geno_candidate_mod_sites),
                "--candidate-bed", str(self.paths.geno_candidate_mod_bed),
                "--read-assignments", str(self.paths.geno_read_assignments),
                "--reference-fa", str(self.reference_fa),
                "--out-tsv", str(self.paths.geno_molecule_mod_calls),
                "--threads", str(max(1, min(8, self.top_threads))),
                "--jobs", str(geno_jobs),
                "--verbose",
            ],
            label="build_molecule_mod_table",
        )
        self.run_python_script(
            "test_snp_transcript_assoc.py",
            [
                "--molecule-snps", str(self.paths.geno_molecule_snps),
                "--out-tsv", str(self.paths.geno_snp_tx),
                "--min-allele-reads", str(int(geno.get("min_group_reads", 4))),
                "--min-transcript-reads", str(int(geno.get("min_group_reads", 4))),
                "--test", str(geno.get("test", "auto")),
                "--pseudocount", str(float(geno.get("pseudocount", 0.5))),
            ],
            label="test_snp_transcript_assoc",
        )
        self.run_python_script(
            "test_snp_mod_assoc.py",
            [
                "--molecule-snps", str(self.paths.geno_molecule_snps),
                "--molecule-mods", str(self.paths.geno_molecule_mod_calls),
                "--out-tsv", str(self.paths.geno_snp_mod),
                "--min-allele-reads", str(int(geno.get("min_group_reads", 4))),
                "--min-total-reads", str(int(geno.get("min_group_reads", 4))),
                "--test", str(geno.get("test", "auto")),
                "--pseudocount", str(float(geno.get("pseudocount", 0.5))),
            ],
            label="test_snp_mod_assoc",
        )
        self.run_python_script(
            "test_snp_tx_mod_dependency.py",
            [
                "--molecule-snps", str(self.paths.geno_molecule_snps),
                "--molecule-mods", str(self.paths.geno_molecule_mod_calls),
                "--snp-transcript-assoc", str(self.paths.geno_snp_tx),
                "--snp-mod-assoc", str(self.paths.geno_snp_mod),
                "--out-tsv", str(self.paths.geno_joint),
                "--min-stratum-reads", str(int(geno.get("min_group_reads", 4))),
            ],
            label="test_snp_tx_mod_dependency",
        )
        self.run_python_script(
            "build_haplotype_blocks.py",
            [
                "--molecule-snps", str(self.paths.geno_molecule_snps),
                "--out-blocks-tsv", str(self.paths.geno_hap_blocks),
                "--out-molecules-tsv", str(self.paths.geno_molecule_haps),
                "--min-alt-reads", str(int(geno.get("min_alt_reads", 4))),
                "--min-cocover-reads", str(int(geno.get("min_haplotype_reads", 4))),
                "--max-block-snps", str(int(geno.get("max_haplotype_snps", 4))),
                "--min-haplotype-reads", str(int(geno.get("min_haplotype_reads", 4))),
            ],
            label="build_haplotype_blocks",
        )
        self.run_python_script(
            "test_haplotype_associations.py",
            [
                "--molecule-haplotypes", str(self.paths.geno_molecule_haps),
                "--molecule-mods", str(self.paths.geno_molecule_mod_calls),
                "--out-haplotype-transcript", str(self.paths.geno_hap_tx),
                "--out-haplotype-mod", str(self.paths.geno_hap_mod),
                "--min-haplotype-reads", str(int(geno.get("min_haplotype_reads", 4))),
                "--min-transcript-reads", str(int(geno.get("min_group_reads", 4))),
                "--min-total-reads", str(int(geno.get("min_group_reads", 4))),
                "--test", str(geno.get("test", "auto")),
                "--pseudocount", str(float(geno.get("pseudocount", 0.5))),
            ],
            label="test_haplotype_associations",
        )

    def stage_report(self) -> None:
        report_cfg = self.config.get("report", {})
        if not as_bool(report_cfg.get("enable", True), True):
            return
        for required_path, label in [
            (self.paths.classification_summary, "classification summary TSV"),
            (self.paths.metrics, "metrics TSV"),
            (self.paths.tx_counts, "transcript counts TSV"),
            (self.paths.pca_png, "PCA PNG"),
            (self.paths.sample_stats, "per-sample stats TSV"),
            (self.paths.per_sample_read_stats, "per-sample read stats TSV"),
            (self.paths.tx_assigned_read_lengths, "assigned read lengths TSV"),
            (self.paths.partition_map, "partition map TSV"),
        ]:
            self._require_existing_file(required_path, label)
        args = [
            "--classification", str(self.paths.classification_summary),
            "--metrics", str(self.paths.metrics),
            "--tx-counts", str(self.paths.tx_counts),
            "--pca-png", str(self.paths.pca_png),
            "--sample-stats", str(self.paths.sample_stats),
            "--read-stats", str(self.paths.per_sample_read_stats),
            "--tx-lengths", str(self.paths.tx_assigned_read_lengths),
            "--partition-map", str(self.paths.partition_map),
            "--out-html", str(self.paths.report_html),
            "--title", str(report_cfg.get("title", f"modulator report: {self.prefix}")),
            "--max-diff-figs", str(int(report_cfg.get("max_diff_figs", 6))),
            "--top-transcripts", str(int(report_cfg.get("top_transcripts", 20))),
            "--top-genes", str(int(report_cfg.get("top_genes", 20))),
            "--zn-long", str(self.paths.zn_filtered_long) if self.paths.zn_filtered_long.exists() else "",
            "--zt-long", str(self.paths.zt_filtered_long) if self.paths.zt_filtered_long.exists() else "",
            "--diff-results", str(self.paths.zn_diff_results) if self.paths.zn_diff_results.exists() else "",
            "--diff-figs-dir", str(self.paths.zn_diff_figs) if self.paths.zn_diff_figs.exists() else "",
            "--multigene-summary-glob", str(self.paths.zt_scrap_dir / "*.multigene_filter_summary.tsv") if self.paths.zt_scrap_dir.exists() else "",
            "--candidate-snps", str(self.paths.geno_candidate_snps) if self.paths.geno_candidate_snps.exists() else "",
            "--snp-tx-assoc", str(self.paths.geno_snp_tx) if self.paths.geno_snp_tx.exists() else "",
            "--snp-mod-assoc", str(self.paths.geno_snp_mod) if self.paths.geno_snp_mod.exists() else "",
            "--snp-tx-mod-assoc", str(self.paths.geno_joint) if self.paths.geno_joint.exists() else "",
            "--hap-blocks", str(self.paths.geno_hap_blocks) if self.paths.geno_hap_blocks.exists() else "",
            "--hap-tx-assoc", str(self.paths.geno_hap_tx) if self.paths.geno_hap_tx.exists() else "",
            "--hap-mod-assoc", str(self.paths.geno_hap_mod) if self.paths.geno_hap_mod.exists() else "",
        ]
        self.run_python_script("generate_html_report.py", args, label="generate_html_report")
