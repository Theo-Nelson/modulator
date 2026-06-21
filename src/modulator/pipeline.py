from __future__ import annotations

import glob
import os
import sys
import time
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
    "aggregate_zn",
    "test_diffs",
    "classify_diffs",
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
    def zn_site_classified(self) -> Path:
        return self.test_diffs / f"{self.prefix}__ZN_site_classified.tsv"

    @property
    def zn_diff_figs(self) -> Path:
        return self.test_diffs / f"{self.prefix}__figs"

    @property
    def zn_class_figs(self) -> Path:
        return self.test_diffs / f"{self.prefix}__figs_by_category"

    @property
    def zn_class_figs_arch(self) -> Path:
        return self.test_diffs / f"{self.prefix}__figs_by_category_arch"

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
    def __init__(self, config: dict[str, Any], *, workdir: str | Path, jobs: int = 1, verbose: bool = True, resume: bool = False):
        self.root = find_project_root(workdir)
        ensure_mpl_config_dir(self.root)
        self.config = config
        self.jobs = max(1, int(jobs))
        self.verbose = verbose
        self.resume = bool(resume)
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

    @property
    def _checkpoint_dir(self) -> Path:
        return self.paths.results / ".checkpoints"

    def _stage_marker(self, stage: str) -> Path:
        return self._checkpoint_dir / f"{stage}.done"

    def _mark_stage_done(self, stage: str) -> None:
        try:
            self._checkpoint_dir.mkdir(parents=True, exist_ok=True)
            self._stage_marker(stage).write_text(f"done\t{stage}\n")
        except Exception:
            pass

    def _nonempty(self, path: Path) -> bool:
        try:
            return path.exists() and path.stat().st_size > 0
        except Exception:
            return False

    def _geno_reuse(self, out_path: Path, label: str) -> bool:
        """Within the genotype stage, reuse an existing non-empty output on --resume so a
        re-run after a later-step failure doesn't redo the expensive per-(sample x chrom)
        BAM scans. Applied only to the deterministic SNP-side tables, so reuse is exact;
        everything from candidate_mod_sites onward always re-runs."""
        if self.resume and self._nonempty(out_path):
            if self.verbose:
                print(f"[modulator]   genotype substep {label}: reusing existing output, skipping", flush=True)
            return True
        return False

    def _all_samples(self, fn) -> bool:
        return bool(self.samples) and all(self._nonempty(fn(sample)) for sample in self.samples)

    def _modkit_done(self, which: str) -> bool:
        base = self.paths.modkit_zn if which == "zn" else self.paths.modkit_zt
        if not base.exists() or not self.samples:
            return False
        for sample in self.samples:
            d = self.paths.modkit_dir(which, sample)
            if not d.exists() or not any(d.rglob("*.bed.gz")):
                return False
        return True

    def _outputs_present(self, stage: str) -> bool:
        """Best-effort check that a stage's outputs already exist, so --resume
        works on a results folder produced before checkpoint markers existed.
        Disabled/toggled-off stages count as already satisfied (they no-op)."""
        cfg = self.config
        p = self.paths
        toggles = cfg.get("toggles", {})
        if stage == "assemble":
            ok = all(self._nonempty(x) for x in (
                p.out_gtf, p.classification_summary, p.metrics, p.tx_counts,
                p.partition_map, p.sample_stats, p.tx_assigned_read_lengths))
            if ok and as_bool(cfg.get("assembler", {}).get("write_zt_tagged_sample_bams", True), True):
                ok = self._all_samples(p.zt_tagged_bam)
            return ok
        if stage == "read_stats":
            return self._nonempty(p.per_sample_read_stats)
        if stage == "multigene_filter":
            if not as_bool(cfg.get("multigene_filter", {}).get("enable", True), True):
                return True
            return self._all_samples(p.clean_bam) and self._nonempty(p.multigene_scrap_tx_counts)
        if stage == "modkit_zn":
            if not as_bool(toggles.get("enable_zn_pileup", True), True):
                return True
            return self._modkit_done("zn")
        if stage == "aggregate_zn":
            if not as_bool(toggles.get("enable_zn_aggregate", True), True):
                return True
            return self._nonempty(p.zn_filtered_long)
        if stage == "test_diffs":
            if not as_bool(toggles.get("enable_test_diffs", True), True):
                return True
            return self._nonempty(p.zn_diff_results)
        if stage == "classify_diffs":
            if not as_bool(cfg.get("classify_diffs", {}).get("enable", True), True):
                return True
            return self._nonempty(p.zn_site_classified)
        if stage == "genotype":
            if not as_bool(cfg.get("genotype", {}).get("enable", False), False):
                return True
            return self._nonempty(p.geno_hap_mod)
        if stage == "report":
            if not as_bool(cfg.get("report", {}).get("enable", True), True):
                return True
            return self._nonempty(p.report_html)
        return False

    def _stage_done(self, stage: str) -> bool:
        return self._stage_marker(stage).exists() or self._outputs_present(stage)

    def run(self, stages: list[str] | None = None) -> None:
        selected = STAGE_ORDER if not stages else [stage for stage in STAGE_ORDER if stage in stages]
        timings: list[tuple[str, float]] = []
        for stage in selected:
            if self.resume and self._stage_done(stage):
                if self.verbose:
                    print(f"[modulator] stage: {stage} — checkpoint found, skipping", flush=True)
                continue
            if self.verbose:
                print(f"[modulator] stage: {stage}", flush=True)
            t0 = time.perf_counter()
            getattr(self, f"stage_{stage}")()
            dt = time.perf_counter() - t0
            timings.append((stage, dt))
            print(f"[modulator] stage {stage} finished in {dt:.1f}s", flush=True)
            self._mark_stage_done(stage)
        if timings:
            try:
                tpath = self.paths.results / "stage_timings.tsv"
                tpath.parent.mkdir(parents=True, exist_ok=True)
                with open(tpath, "w") as fh:
                    fh.write("stage\tseconds\n")
                    for stg, secs in timings:
                        fh.write(f"{stg}\t{secs:.2f}\n")
                    fh.write(f"TOTAL\t{sum(s for _, s in timings):.2f}\n")
                print(f"[modulator] stage timings -> {tpath}", flush=True)
            except OSError:
                pass

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
        # Default to the streaming engine: it k-way-merges the already-sorted,
        # tabix-indexed per-ZN beds (no normalize.tsv, no genome-wide external sort),
        # is parallel across chromosomes, and is per-chromosome resumable. Output is
        # content-identical to the sort engine (validated). Set aggregation.engine=sort
        # to fall back to aggregate_by_gene.py.
        engine = str(self.config.get("aggregation", {}).get("engine", "stream")).strip().lower()
        if engine == "stream":
            agg_jobs = max(1, int(self.config.get("aggregation", {}).get("jobs", min(self.top_threads, 12) if self.top_threads else 8)))
            args.extend(["--jobs", str(agg_jobs)])
            self.run_python_script("aggregate_zn_stream.py", args, label="aggregate_zn_stream")
        else:
            self.run_python_script("aggregate_by_gene.py", args, label="aggregate_by_gene")

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

    def stage_classify_diffs(self) -> None:
        cfg = self.config.get("classify_diffs", {})
        if not as_bool(cfg.get("enable", True), True):
            return
        if not self.paths.zn_diff_results.exists() or self.paths.zn_diff_results.stat().st_size == 0:
            if self.verbose:
                print(
                    f"[modulator] skipping classify_diffs because the ZN diff results table is missing or empty: {self.paths.zn_diff_results}",
                    flush=True,
                )
            return
        if not self.paths.out_gtf.exists():
            if self.verbose:
                print(
                    f"[modulator] skipping classify_diffs because the assembled GTF is missing: {self.paths.out_gtf}",
                    flush=True,
                )
            return
        args = [
            "--diff-tsv", str(self.paths.zn_diff_results),
            "--gtf", str(self.paths.out_gtf),
            "--out-tsv", str(self.paths.zn_site_classified),
            "--min-effect", str(cfg.get("min_effect", 0.10)),
            "--fdr", str(cfg.get("fdr", 0.05)),
            "--min-cov", str(cfg.get("min_cov", 0)),
            "--tes-tol", str(cfg.get("tes_tol", 200)),
            "--inside-tol", str(cfg.get("inside_tol", 50)),
            "--ejc-nt", str(cfg.get("ejc_nt", 150)),
            "--intergenic-gap", str(cfg.get("intergenic_gap", 1000)),
            "--verbose",
        ]
        if as_bool(cfg.get("figures", True), True):
            # Isoform architecture-map figures are the PRIMARY per-category figure.
            # They are built from the GTF isoform models + the diff table's
            # per_transcript_json, so they need no --zn-long.
            args.extend([
                "--arch-figs-dir", str(self.paths.zn_class_figs_arch),
                "--figs-per-category", str(int(cfg.get("figs_per_category", 10))),
            ])
            # The 2-panel per-sample stoichiometry figures additionally need --zn-long.
            if self.paths.zn_filtered_long.exists():
                args.extend([
                    "--zn-long", str(self.paths.zn_filtered_long),
                    "--figs-dir", str(self.paths.zn_class_figs),
                ])
        mod_filter = cfg.get("mod_filter")
        if mod_filter is None:
            mod_filter = self.config.get("test_diffs", {}).get("mod_filter")
        # An empty/None mod_filter means classify ALL modifications -- the diff table
        # already carries every mod_code emitted upstream. Only restrict when set.
        if mod_filter:
            args.append("--mod-filter")
            args.extend(str(m) for m in mod_filter)
        self.run_python_script("classify_diff_sites.py", args, label="classify_diff_sites")

    def stage_genotype(self) -> None:
        geno = self.config.get("genotype", {})
        if not as_bool(geno.get("enable", False), False):
            return
        self._require_reference_fa()
        # The genotype scripts now shard BAM scans per (sample x chromosome), so
        # parallelism is no longer capped at the sample count -- size it to the CPU
        # budget instead.
        geno_jobs = max(1, min(self.top_threads if self.top_threads > 0 else 8, int(geno.get("jobs", self.top_threads or 8))))
        sample_bams = [
            str(self._require_existing_file(self._modkit_input_bam(sample), f"genotype input BAM for sample {sample}"))
            for sample in self.samples
        ]
        self._require_existing_file(self.paths.classification_summary, "classification summary TSV")
        self.paths.genotype.mkdir(parents=True, exist_ok=True)

        if not self._geno_reuse(self.paths.geno_read_assignments, "build_read_assignment_table"):
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
        if not self._geno_reuse(self.paths.geno_candidate_snps, "discover_candidate_snps"):
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
        if not self._geno_reuse(self.paths.geno_molecule_snps, "build_molecule_snp_table"):
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
        mod_site_args = [
            "--zn-long", zn_long,
            "--zt-long", zt_long,
            "--out-tsv", str(self.paths.geno_candidate_mod_sites),
            "--out-bed", str(self.paths.geno_candidate_mod_bed),
            "--min-total-cov", str(int(geno.get("min_mod_site_cov", 1))),
        ]
        # Restrict candidate mod sites to those that can pair with a candidate SNP (same
        # context_key on a shared read). Lossless for snp_mod_assoc/snp_tx_mod_dependency/
        # haplotype_mod_assoc and keeps the per-read mod table tractable on deep genome-wide
        # data (otherwise 100k+ sites -> molecule_mod_calls OOM). Toggle off to keep all sites.
        if as_bool(geno.get("mod_sites_require_snp_link", True), True):
            mod_site_args += ["--candidate-snps", str(self.paths.geno_candidate_snps)]
        self.run_python_script(
            "build_candidate_mod_sites.py",
            mod_site_args,
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
                "--threads", str(max(1, self.top_threads)),
                # Cap concurrent modkit `extract calls` separately from the (cheap, pysam)
                # SNP steps: each modkit process streams a chromosome's reads, so too many
                # at once is memory-heavy. Defaults to 24, never more than geno_jobs.
                "--jobs", str(max(1, min(geno_jobs, int(geno.get("mod_jobs", 8))))),
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
            "--classified-sites", str(self.paths.zn_site_classified) if self.paths.zn_site_classified.exists() else "",
            "--class-figs-dir", str(self.paths.zn_class_figs) if self.paths.zn_class_figs.exists() else "",
            "--arch-figs-dir", str(self.paths.zn_class_figs_arch) if self.paths.zn_class_figs_arch.exists() else "",
            "--max-class-figs-per-category", str(int(report_cfg.get("max_class_figs_per_category", 10))),
            "--multigene-summary-glob", str(self.paths.zt_scrap_dir / "*.multigene_filter_summary.tsv") if self.paths.zt_scrap_dir.exists() else "",
            "--candidate-snps", str(self.paths.geno_candidate_snps) if self.paths.geno_candidate_snps.exists() else "",
            "--snp-tx-assoc", str(self.paths.geno_snp_tx) if self.paths.geno_snp_tx.exists() else "",
            "--snp-mod-assoc", str(self.paths.geno_snp_mod) if self.paths.geno_snp_mod.exists() else "",
            "--snp-tx-mod-assoc", str(self.paths.geno_joint) if self.paths.geno_joint.exists() else "",
            "--hap-blocks", str(self.paths.geno_hap_blocks) if self.paths.geno_hap_blocks.exists() else "",
            "--hap-tx-assoc", str(self.paths.geno_hap_tx) if self.paths.geno_hap_tx.exists() else "",
            "--hap-mod-assoc", str(self.paths.geno_hap_mod) if self.paths.geno_hap_mod.exists() else "",
            "--snp-figs-dir", str(self.paths.genotype / f"{self.prefix}__snp_figs"),
            "--max-snp-figs", str(int(report_cfg.get("max_snp_figs", 12))),
        ]
        self.run_python_script("generate_html_report.py", args, label="generate_html_report")
