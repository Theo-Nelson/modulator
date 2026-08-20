from __future__ import annotations

import glob
import os
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from modulator.runtime import (
    PeakRSSSampler,
    as_bool,
    ensure_mpl_config_dir,
    ensure_parent,
    find_project_root,
    format_command,
    is_set,
    normalize_pivot_mode,
    require_tools,
    resolve_path,
    run_command,
    run_parallel,
)
from modulator.samplesheet import (
    read_samplesheet,
    resolve_contrasts,
    sample_metadata,
    stage_bams,
    write_metadata_tsv,
)


STAGE_ORDER = [
    "assemble",
    "read_stats",
    "splice_junctions",
    "apa_motifs",
    "multigene_filter",
    "modkit_zn",
    "aggregate_zn",
    "novel_loci",
    "sequence_elements",
    "test_diffs",
    "classify_diffs",
    "genotype",
    "polya",
    "hierarchical_stoich",
    "between_conditions",
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
    def splice_junctions(self) -> Path:
        return self.assemble / f"{self.prefix}_splice_junctions.tsv"

    @property
    def gene_splice_summary(self) -> Path:
        return self.assemble / f"{self.prefix}_gene_splice_summary.tsv"

    @property
    def apa_motifs(self) -> Path:
        return self.assemble / f"{self.prefix}_apa_motifs.tsv"

    @property
    def geno_snp_mod_mechanism(self) -> Path:
        return self.genotype / f"{self.prefix}_snp_mod_mechanism.tsv"

    @property
    def geno_snp_at_mod_base(self) -> Path:
        return self.genotype / f"{self.prefix}_snp_at_mod_base.tsv"

    @property
    def novel_loci_tsv(self) -> Path:
        return self.assemble / f"{self.prefix}_novel_loci.tsv"

    @property
    def novel_fragmentforms(self) -> Path:
        return self.assemble / f"{self.prefix}_novel_fragmentforms.tsv"

    @property
    def sequence_elements(self) -> Path:
        return self.assemble / f"{self.prefix}_sequence_elements.tsv"

    @property
    def sequence_elements_summary(self) -> Path:
        return self.assemble / f"{self.prefix}_sequence_elements_summary.tsv"

    @property
    def geno_mod_mod(self) -> Path:
        return self.genotype / f"{self.prefix}_mod_mod_assoc.tsv"

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
    def zn_site_private(self) -> Path:
        return self.test_diffs / f"{self.prefix}__ZN_site_private.tsv"

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
    def gene_browser_html(self) -> Path:
        return self.report / f"{self.prefix}_gene_browser.html"

    @property
    def geno_read_assignments(self) -> Path:
        return self.genotype / f"{self.prefix}_read_assignments.tsv"

    @property
    def geno_read_assignments_regions(self) -> Path:
        """Read-assignment table restricted to reads overlapping candidate SNP/mod sites.
        Deliberately named differently from the genome-wide table so the narrower scope is explicit."""
        return self.genotype / f"{self.prefix}_read_assignments_candidate_regions.tsv"

    @property
    def geno_candidate_regions_bed(self) -> Path:
        return self.genotype / f"{self.prefix}_candidate_regions.bed"

    @property
    def geno_subset_dir(self) -> Path:
        return self.genotype / "subset_bams"

    def geno_subset_bam(self, sample_bam: Path) -> Path:
        # Keep the SAME basename so genotype_utils.sample_name_from_bam() derives an identical
        # sample name from the subset BAM as from the original.
        return self.geno_subset_dir / Path(sample_bam).name

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

    @property
    def sample_metadata(self) -> Path:
        return self.results / f"{self.prefix}_sample_metadata.tsv"

    @property
    def staged_bams_dir(self) -> Path:
        return self.results / "staged_bams"

    @property
    def hierarchical_stoich(self) -> Path:
        return self.test_diffs / f"{self.prefix}_hierarchical_stoich.tsv"

    @property
    def between_conditions(self) -> Path:
        return self.results / "between_conditions"

    def cond_mod_diffs(self, contrast: str) -> Path:
        return self.between_conditions / f"{self.prefix}_{contrast}_mod_diffs.tsv"

    def cond_usage_diffs(self, contrast: str, feature: str) -> Path:
        return self.between_conditions / f"{self.prefix}_{contrast}_{feature}_usage_diffs.tsv"

    def cond_tail_diffs(self, contrast: str) -> Path:
        return self.between_conditions / f"{self.prefix}_{contrast}_tail_diffs.tsv"

    @property
    def polya(self) -> Path:
        return self.results / "polya"

    @property
    def polya_read_tails(self) -> Path:
        return self.polya / f"{self.prefix}_read_tail_lengths.tsv"

    @property
    def polya_fragmentform(self) -> Path:
        return self.polya / f"{self.prefix}_polya_fragmentform.tsv"

    @property
    def polya_taillength_diffs(self) -> Path:
        return self.polya / f"{self.prefix}_taillength_diffs.tsv"

    @property
    def polya_taillength_mod(self) -> Path:
        return self.polya / f"{self.prefix}_taillength_mod.tsv"

    @property
    def polya_diff_figs(self) -> Path:
        return self.polya / f"{self.prefix}__taillength_diff_figs"

    @property
    def polya_mod_figs(self) -> Path:
        return self.polya / f"{self.prefix}__taillength_mod_figs"

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
        # When set to a list (only during the genotype stage), run_python_script records
        # each substep's peak RSS into it, so per-substep memory can be written out.
        self._substep_mem: list[tuple[str, float]] | None = None
        self.prefix = str(config.get("prefix", "modulator_run"))
        self.paths = PipelinePaths(self.root, self.prefix)
        # Samplesheet (optional) is the sample SOURCE + the condition metadata; it stages the BAMs
        # and therefore must run before sample discovery.
        self._staged_bams_dir: Path | None = None
        self.sample_meta: dict[str, dict] = {}
        self.contrasts: list[dict] = []
        self._load_samplesheet()
        self.samples = self._discover_samples()
        self.reference_fa = self._resolve_reference("reference_fa", ("reference", "fasta"))
        self.reference_gtf = self._resolve_reference("reference_gtf", ("reference", "gtf"))
        self.top_threads = int(config.get("threads", 1))
        self._validate_config()

    def _load_samplesheet(self) -> None:
        """If a samplesheet is configured it becomes the sample source and the metadata.

        Each BAM is symlinked to ``<results>/staged_bams/<sample>.bam`` so the sample id IS the BAM
        stem for every downstream script (see samplesheet.py). Without a samplesheet the pipeline
        keeps its original bams_dir + bam_glob behaviour.
        """
        ss = self.config.get("samplesheet")
        if not is_set(ss):
            return
        rows = read_samplesheet(resolve_path(self.root, str(ss)))
        stage_bams(rows, self.paths.staged_bams_dir, self._config_bams_dir)
        self._staged_bams_dir = self.paths.staged_bams_dir
        self.sample_meta = sample_metadata(rows)
        self.contrasts = resolve_contrasts(self.config.get("contrasts"), rows)
        write_metadata_tsv(rows, self.paths.sample_metadata)
        if self.verbose:
            groups: dict[str, list[str]] = {}
            for r in rows:
                groups.setdefault(r.get("condition") or "(no condition)", []).append(r["sample"])
            print(f"[modulator] samplesheet: {len(rows)} sample(s) staged -> {self._staged_bams_dir}", flush=True)
            for cond, names in groups.items():
                print(f"[modulator]   {cond}: {', '.join(names)}", flush=True)
            if self.contrasts:
                print(f"[modulator]   contrasts: {', '.join(c['name'] for c in self.contrasts)}", flush=True)

    def _discover_samples(self) -> list[str]:
        bams_dir = self.bams_dir
        bam_glob = self.bam_glob
        found = sorted(glob.glob(str(bams_dir / bam_glob)))
        samples = [Path(path).stem for path in found]
        if not samples:
            raise FileNotFoundError(f"No BAMs matched {bam_glob!r} under {bams_dir}")
        return samples

    @property
    def _config_bams_dir(self) -> Path:
        """The bams_dir as configured -- samplesheet 'bam' values resolve relative to this."""
        return resolve_path(self.root, self.config.get("bams_dir", "resources/test_bams/ALCAM_NHSL1_SERAC1_MXD1_RIOK3_reads"))

    @property
    def bams_dir(self) -> Path:
        # With a samplesheet the staged symlink dir IS the bams dir, so every stage sees {sample}.bam.
        return self._staged_bams_dir or self._config_bams_dir

    @property
    def bam_glob(self) -> str:
        # Staged BAMs are named <sample>.bam, so the configured glob no longer applies.
        if self._staged_bams_dir is not None:
            return "*.bam"
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
        self._preflight_bams()

    @staticmethod
    def _bam_has_index(bam_path: str) -> bool:
        """True if a coordinate index (.bai/.csi) exists for this BAM, using every name modkit and
        pysam will try: <bam>.bai, <bam>.csi, and the samtools <stem>.bai/.csi variant."""
        p = str(bam_path)
        cands = [p + ".bai", p + ".csi"]
        if p.endswith(".bam"):
            cands += [p[:-4] + ".bai", p[:-4] + ".csi"]
        return any(os.path.exists(c) for c in cands)

    def _preflight_bams(self) -> None:
        """Fail fast, before any stage runs, on two input mistakes that otherwise surface as a
        cryptic mid-run modkit error or (worse) silently wrong statistics:

          * duplicate BAMs -- two samples resolving to the SAME physical file, which double-counts
            one library as two replicates (pseudo-replication). Errors. Distinct files of identical
            size are flagged as a warning (likely an accidental copy) but allowed.
          * a missing .bai/.csi index -- modkit pileup and pysam both require one. Errors with the
            exact `samtools index` command to fix, or builds them if preflight.auto_index is set.

        Controlled by the `preflight` config block; set preflight.enable: false to skip entirely.
        """
        pf = self.config.get("preflight", {}) or {}
        if not as_bool(pf.get("enable", True), True):
            return
        paths = sorted(glob.glob(str(self.bams_dir / self.bam_glob)))
        if not paths:
            return  # empty discovery is reported elsewhere (_discover_samples)

        problems: list[str] = []
        warnings: list[str] = []

        if as_bool(pf.get("check_duplicate_bams", True), True):
            by_inode: dict[tuple[int, int], str] = {}
            by_size: dict[int, list[str]] = {}
            for p in paths:
                try:
                    st = os.stat(os.path.realpath(p))
                except OSError as exc:
                    problems.append(f"cannot stat BAM '{Path(p).name}': {exc}")
                    continue
                key = (st.st_dev, st.st_ino)
                if key in by_inode:
                    problems.append(
                        f"duplicate BAM: '{Path(p).name}' and '{Path(by_inode[key]).name}' are the "
                        f"SAME physical file ({os.path.realpath(p)}); this double-counts one library "
                        f"as two samples (pseudo-replication)."
                    )
                else:
                    by_inode[key] = p
                by_size.setdefault(st.st_size, []).append(p)
            for sz, plist in by_size.items():
                realpaths = {os.path.realpath(x) for x in plist}
                if len(realpaths) > 1:
                    names = ", ".join(sorted(Path(x).name for x in plist))
                    warnings.append(
                        f"{len(realpaths)} distinct BAMs share an identical size ({sz} bytes): "
                        f"{names}; verify these are not accidental copies of one library."
                    )

        if as_bool(pf.get("check_bam_index", True), True):
            missing = [p for p in paths if not self._bam_has_index(p)]
            if missing and as_bool(pf.get("auto_index", False), False):
                require_tools(["samtools"])
                threads = str(max(1, self.top_threads))
                still: list[str] = []
                for p in missing:
                    try:
                        run_command(["samtools", "index", "-@", threads, p],
                                    cwd=self.root, label=f"preflight_index[{Path(p).name}]",
                                    verbose=self.verbose)
                    except Exception as exc:  # noqa: BLE001
                        warnings.append(f"auto-index failed for '{Path(p).name}': {exc}")
                    if not self._bam_has_index(p):
                        still.append(p)
                missing = still
            if missing:
                listed = "\n".join(f"      - {Path(p).name}" for p in missing)
                fix = "\n".join(f"      samtools index {os.path.realpath(p)}" for p in missing)
                problems.append(
                    "BAM(s) missing a .bai/.csi index (modkit pileup and pysam require one):\n"
                    + listed
                    + "\n    create them with:\n" + fix
                    + "\n    or set preflight.auto_index: true to build them automatically."
                )

        for w in warnings:
            print(f"[modulator] preflight WARNING: {w}", flush=True)
        if problems:
            raise ValueError(
                "modulator preflight failed for the input BAMs:\n  - "
                + "\n  - ".join(problems)
                + "\n(If this is intentional, set preflight.enable: false in the config.)"
            )
        if self.verbose:
            print(f"[modulator] preflight: {len(paths)} BAM(s) OK (unique files, indexes present).",
                  flush=True)

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
        sink = self._substep_mem
        on_peak = (lambda lbl, gib: sink.append((lbl, gib))) if sink is not None else None
        run_command(cmd, cwd=self.root, label=label, verbose=self.verbose, on_peak=on_peak)

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
            if not d.exists():
                return False
            beds_gz = list(d.rglob("*.bed.gz"))
            if not beds_gz:
                return False
            # Every compressed partition must be tabix-indexed, and no raw .bed may be
            # left behind (other than the intentionally-skipped ungrouped.bed) -- a leftover
            # raw .bed means the bgzip/tabix loop was interrupted, so aggregate_zn_stream
            # would silently drop those ZN partitions. Require a fully-finished directory.
            if any(not Path(str(b) + ".tbi").exists() for b in beds_gz):
                return False
            if any(b.name != "ungrouped.bed" for b in d.rglob("*.bed")):
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
        if stage == "splice_junctions":
            if not as_bool(cfg.get("splice_junctions", {}).get("enable", True), True):
                return True
            return self._nonempty(p.gene_splice_summary)
        if stage == "apa_motifs":
            if not as_bool(cfg.get("apa_motifs", {}).get("enable", True), True):
                return True
            return self._nonempty(p.apa_motifs)
        if stage == "novel_loci":
            if not as_bool(cfg.get("novel_loci", {}).get("enable", True), True):
                return True
            return self._nonempty(p.novel_loci_tsv)
        if stage == "sequence_elements":
            if not as_bool(cfg.get("sequence_elements", {}).get("enable", True), True):
                return True
            return self._nonempty(p.sequence_elements)
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
        if stage == "polya":
            if not as_bool(cfg.get("polya", {}).get("enable", True), True):
                return True
            return self._nonempty(p.polya_read_tails) and self._nonempty(p.polya_taillength_diffs)
        if stage == "hierarchical_stoich":
            if not as_bool(cfg.get("hierarchical_stoich", {}).get("enable", False), False):
                return True
            return self._nonempty(p.hierarchical_stoich)
        if stage == "between_conditions":
            # No samplesheet/contrasts -> the stage is a no-op, so count it as satisfied.
            if not as_bool(cfg.get("between_conditions", {}).get("enable", True), True) or not self.contrasts:
                return True
            return all(self._nonempty(p.cond_mod_diffs(c["name"])) for c in self.contrasts)
        if stage == "report":
            if not as_bool(cfg.get("report", {}).get("enable", True), True):
                return True
            return self._nonempty(p.report_html)
        return False

    def _stage_done(self, stage: str) -> bool:
        return self._stage_marker(stage).exists() or self._outputs_present(stage)

    def _stage_disabled(self, stage: str) -> bool:
        """True when the stage is toggled off (or a no-op, like between_conditions with no
        contrasts) rather than actually completed. Used only to report an accurate skip reason:
        a disabled stage must not be described as "checkpoint found / outputs reused"."""
        cfg = self.config
        toggles = cfg.get("toggles", {})
        checks = {
            "splice_junctions":    lambda: not as_bool(cfg.get("splice_junctions", {}).get("enable", True), True),
            "apa_motifs":          lambda: not as_bool(cfg.get("apa_motifs", {}).get("enable", True), True),
            "novel_loci":          lambda: not as_bool(cfg.get("novel_loci", {}).get("enable", True), True),
            "sequence_elements":   lambda: not as_bool(cfg.get("sequence_elements", {}).get("enable", True), True),
            "multigene_filter":    lambda: not as_bool(cfg.get("multigene_filter", {}).get("enable", True), True),
            "modkit_zn":           lambda: not as_bool(toggles.get("enable_zn_pileup", True), True),
            "aggregate_zn":        lambda: not as_bool(toggles.get("enable_zn_aggregate", True), True),
            "test_diffs":          lambda: not as_bool(toggles.get("enable_test_diffs", True), True),
            "classify_diffs":      lambda: not as_bool(cfg.get("classify_diffs", {}).get("enable", True), True),
            "genotype":            lambda: not as_bool(cfg.get("genotype", {}).get("enable", False), False),
            "polya":               lambda: not as_bool(cfg.get("polya", {}).get("enable", True), True),
            "hierarchical_stoich": lambda: not as_bool(cfg.get("hierarchical_stoich", {}).get("enable", False), False),
            "between_conditions":  lambda: (not as_bool(cfg.get("between_conditions", {}).get("enable", True), True)) or not self.contrasts,
            "report":              lambda: not as_bool(cfg.get("report", {}).get("enable", True), True),
        }
        fn = checks.get(stage)
        return bool(fn()) if fn else False

    def _merge_stage_table(self, path: Path, value_header: str,
                           new_rows: list[tuple[str, float]], fmt: str,
                           summary_label: str, summary_fn) -> None:
        """Write a per-stage TSV, MERGING with any existing file so a partial
        ``--stages`` run updates only the stages it actually ran instead of
        clobbering the rows of the stages it skipped. The summary row
        (TOTAL/MAX) is recomputed over the merged set."""
        merged: dict[str, float] = {}
        if path.exists():
            try:
                for i, line in enumerate(path.read_text().splitlines()):
                    if i == 0 or not line.strip():
                        continue
                    stg, _, val = line.partition("\t")
                    if stg == summary_label or not stg or not val:
                        continue
                    try:
                        merged[stg] = float(val)
                    except ValueError:
                        continue
            except OSError:
                pass
        for stg, val in new_rows:
            merged[stg] = val
        ordered = ([s for s in STAGE_ORDER if s in merged]
                   + [s for s in merged if s not in STAGE_ORDER])
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, "w") as fh:
            fh.write(f"stage\t{value_header}\n")
            for stg in ordered:
                fh.write(f"{stg}\t{format(merged[stg], fmt)}\n")
            if merged:
                fh.write(f"{summary_label}\t{format(summary_fn(list(merged.values())), fmt)}\n")

    def _write_run_manifest(self) -> "Path | None":
        """Write a human-readable record of exactly how this run was invoked -- timestamp, command
        line, resolved inputs (with the paths they came from), the sample sheet, and the fully
        merged configuration -- so the report can show precisely where every input and parameter
        originated. Returns the manifest path (or None if it could not be written)."""
        try:
            import datetime
            stamp = datetime.datetime.now().astimezone().strftime("%Y-%m-%d %H:%M:%S %Z").strip()
        except Exception:
            stamp = "(unavailable)"
        L = ["modulator run manifest", "=" * 64,
             f"Run date / time : {stamp}",
             f"Prefix          : {self.prefix}",
             f"Project root    : {self.root}",
             f"Command         : {sys.executable} {' '.join(sys.argv)}",
             "",
             "Resolved inputs", "-" * 64,
             f"reference_fa    : {self.reference_fa}",
             f"reference_gtf   : {self.reference_gtf}",
             f"bams_dir        : {self.bams_dir}",
             f"bam_glob        : {self.bam_glob}",
             f"samplesheet     : {self.config.get('samplesheet') or '(none — bams_dir/bam_glob discovery)'}",
             f"threads         : {self.config.get('threads')}",
             f"jobs            : {self.config.get('jobs', '(default)')}",
             ""]
        try:
            paths = sorted(glob.glob(str(self.bams_dir / self.bam_glob)))
            meta = getattr(self, "sample_meta", {}) or {}
            L += [f"Samples ({len(paths)})", "-" * 64]
            for p in paths:
                s = Path(p).stem
                m = meta.get(s, {})
                extra = ("  " + "  ".join(f"{k}={v}" for k, v in m.items())) if m else ""
                L.append(f"  {s}{extra}  ->  {os.path.realpath(p)}")
            L.append("")
        except Exception:
            pass
        L += ["Fully-resolved configuration (base config + every --set override)", "-" * 64]
        try:
            import yaml
            L.append(yaml.safe_dump(self.config, sort_keys=False, default_flow_style=False).rstrip())
        except Exception:
            import json
            L.append(json.dumps(self.config, indent=2, default=str))
        try:
            out = self.paths.results / f"{self.prefix}_run_manifest.txt"
            out.parent.mkdir(parents=True, exist_ok=True)
            out.write_text("\n".join(L) + "\n")
            return out
        except Exception:
            return None

    def run(self, stages: list[str] | None = None) -> None:
        self._write_run_manifest()
        selected = STAGE_ORDER if not stages else [stage for stage in STAGE_ORDER if stage in stages]
        timings: list[tuple[str, float]] = []
        mem_peaks: list[tuple[str, float]] = []
        geno_substeps: list[tuple[str, float]] = []
        for stage in selected:
            if self.resume and self._stage_done(stage):
                if self.verbose:
                    # Report the ACCURATE reason: a disabled/no-op stage was never run (so nothing
                    # was "reused"); a real checkpoint marker means outputs were reused; otherwise
                    # pre-existing outputs (from a pre-marker run) satisfied the resume check.
                    if self._stage_disabled(stage):
                        reason = "disabled (off) — skipping"
                    elif self._stage_marker(stage).exists():
                        reason = "checkpoint found — reusing existing outputs"
                    else:
                        reason = "existing outputs found (no checkpoint) — skipping"
                    print(f"[modulator] stage: {stage} — {reason}", flush=True)
                continue
            if self.verbose:
                print(f"[modulator] stage: {stage}", flush=True)
            # Per-substep memory is only meaningful for the (sequential) genotype scripts.
            self._substep_mem = [] if stage == "genotype" else None
            t0 = time.perf_counter()
            with PeakRSSSampler(os.getpid()) as sampler:  # peak over this process's whole subtree
                getattr(self, f"stage_{stage}")()
            dt = time.perf_counter() - t0
            timings.append((stage, dt))
            mem_peaks.append((stage, sampler.peak_gib))
            if self._substep_mem:
                geno_substeps.extend(self._substep_mem)
            self._substep_mem = None
            print(f"[modulator] stage {stage} finished in {dt:.1f}s  peak={sampler.peak_gib:.2f} GiB", flush=True)
            # Only checkpoint a stage that actually produced its outputs. A stage that
            # early-returned as a no-op on missing/empty input (e.g. test_diffs with an empty
            # ZN long table) must NOT be marked done, or --resume would skip it forever even
            # after the upstream input appears. Disabled/toggled-off stages report present.
            if self._outputs_present(stage):
                self._mark_stage_done(stage)
        if timings:
            try:
                tpath = self.paths.results / "stage_timings.tsv"
                self._merge_stage_table(tpath, "seconds", timings, ".2f", "TOTAL", sum)
                print(f"[modulator] stage timings -> {tpath}", flush=True)
            except OSError:
                pass
        if mem_peaks:
            try:
                mpath = self.paths.results / "stage_memory.tsv"
                self._merge_stage_table(mpath, "peak_rss_gib", mem_peaks, ".3f", "MAX", max)
                print(f"[modulator] stage memory -> {mpath}", flush=True)
            except OSError:
                pass
        if geno_substeps:
            try:
                gpath = self.paths.results / "genotype_memory.tsv"
                gpath.parent.mkdir(parents=True, exist_ok=True)
                with open(gpath, "w") as fh:
                    fh.write("substep\tpeak_rss_gib\n")
                    for lbl, gib in geno_substeps:
                        fh.write(f"{lbl}\t{gib:.3f}\n")
                print(f"[modulator] genotype substep memory -> {gpath}", flush=True)
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
            # Samples are scanned in parallel (jobs<=1 stays serial / single-core safe).
            "--jobs", str(self.jobs),
        ]
        if as_bool(cfg.get("primary_only", True), True):
            args.append("--primary-only")
        self.run_python_script("per_sample_read_stats.py", args, label="per_sample_read_stats")

    def stage_splice_junctions(self) -> None:
        """Classify every fragmentform intron as canonical (GT-AG) / semi-canonical (GC-AG) /
        minor U12 (AT-AC) / non-canonical, and summarize per gene."""
        cfg = self.config.get("splice_junctions", {})
        if not as_bool(cfg.get("enable", True), True):
            return
        reference_fa = self._require_reference_fa()
        self._require_existing_file(self.paths.out_gtf, "assembled GTF")
        self.run_python_script(
            "classify_splice_junctions.py",
            [
                "--gtf", str(self.paths.out_gtf),
                "--reference-fa", str(reference_fa),
                "--out-junctions", str(self.paths.splice_junctions),
                "--out-genes", str(self.paths.gene_splice_summary),
                "--verbose",
            ],
            label="classify_splice_junctions",
        )

    def stage_novel_loci(self) -> None:
        """Roll up read-backed NOVEL_LOCUS fragmentforms into uniquely-named loci, with their
        fragmentforms, modification sites, and splice-junction category."""
        cfg = self.config.get("novel_loci", {})
        if not as_bool(cfg.get("enable", True), True):
            return
        self._require_existing_file(self.paths.out_gtf, "assembled GTF")
        args = [
            "--gtf", str(self.paths.out_gtf),
            "--classification", str(self.paths.classification_summary),
            "--out-loci", str(self.paths.novel_loci_tsv),
            "--out-fragmentforms", str(self.paths.novel_fragmentforms),
            "--verbose",
        ]
        if self.paths.zn_filtered_long.exists():
            args.extend(["--zn-long", str(self.paths.zn_filtered_long)])
        if self.paths.gene_splice_summary.exists():
            args.extend(["--splice-genes", str(self.paths.gene_splice_summary)])
        self.run_python_script("summarize_novel_loci.py", args, label="summarize_novel_loci")

    def stage_sequence_elements(self) -> None:
        """Annotate sequence-based cis-elements (PAS, ARE, CPE, GRE, rG4, Kozak, uORF,
        5'TOP, stop context, m6Am) on each fragmentform's mature mRNA and report EVERY
        overlapping modification, unbiased across mod codes. Needs the assembled GTF, the
        reference FASTA+GTF (for start/stop codons), and the ZN modification table."""
        cfg = self.config.get("sequence_elements", {})
        if not as_bool(cfg.get("enable", True), True):
            return
        if not self._nonempty(self.paths.zn_filtered_long):
            if self.verbose:
                print("[modulator] sequence_elements: no ZN modification table (needs aggregate_zn); skipping",
                      flush=True)
            return
        self._require_existing_file(self.paths.out_gtf, "assembled GTF")
        self.run_python_script("annotate_sequence_elements.py", [
            "--assembled-gtf", str(self.paths.out_gtf),
            "--reference-fa", str(self._require_reference_fa()),
            "--reference-gtf", str(self._require_reference_gtf()),
            "--mod-sites", str(self.paths.zn_filtered_long),
            "--out-tsv", str(self.paths.sequence_elements),
            "--out-summary", str(self.paths.sequence_elements_summary),
            "--pas-window", str(int(cfg.get("pas_window", 60))),
            "--utr3-window", str(int(cfg.get("utr3_window", 400))),
            "--verbose",
        ], label="annotate_sequence_elements")

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
                "--multi-gene-action", str(cfg.get("multi_gene_action", "scrap_conflict")),
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
            or str((self.paths.results / "tmp" / self.prefix).resolve())
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
            "--nfail-score-k", str(float(agg_cfg.get("nfail_score_k", filters_cfg.get("nfail_score_k", 0.0)))),
            "--verbose",
        ]
        if as_bool(agg_cfg.get("filter_enable", filters_cfg.get("enable_site_filter", True)), True):
            args.append("--filter-enable")
        args.append("--emit-raw" if as_bool(agg_cfg.get("emit_raw", self.config.get("aggregate_outputs", {}).get("emit_raw", True)), True) else "--no-emit-raw")
        args.append("--emit-filtered" if as_bool(agg_cfg.get("emit_filtered", self.config.get("aggregate_outputs", {}).get("emit_filtered", True)), True) else "--no-emit-filtered")
        args.append("--write-long" if as_bool(agg_cfg.get("write_long", self.config.get("aggregate_outputs", {}).get("write_long", True)), True) else "--no-write-long")
        # Pivots are optional inspection outputs (3 dense files per gene x mod group; nothing
        # downstream reads them). Tri-state: 'auto' (default) writes them unless the run exceeds
        # pivot_max_groups, avoiding a hundreds-of-thousands-of-tiny-files explosion on whole-
        # transcriptome runs; 'on'/true forces them even at scale; 'off'/false never writes them.
        pivot_raw = agg_cfg.get("write_pivots", self.config.get("aggregate_outputs", {}).get("write_pivots", "auto"))
        pivot_mode = normalize_pivot_mode(pivot_raw)
        pivot_max_groups = int(agg_cfg.get("pivot_max_groups", self.config.get("aggregate_outputs", {}).get("pivot_max_groups", 2000)))
        args.extend(["--pivot-mode", pivot_mode, "--pivot-max-groups", str(pivot_max_groups)])
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
            "--tes-tol", str(cfg.get("tes_tol", 25)),
            "--inside-tol", str(cfg.get("inside_tol", 50)),
            "--ejc-nt", str(cfg.get("ejc_nt", 150)),
            "--intergenic-gap", str(cfg.get("intergenic_gap", 1000)),
            "--verbose",
        ]
        # Coverage-independent PRIVATE-site scan (needs the FILTERED long table; independent of the
        # differential test and of figures).
        if self.paths.zn_filtered_long.exists():
            args.extend([
                "--zn-long", str(self.paths.zn_filtered_long),
                "--private-out-tsv", str(self.paths.zn_site_private),
                "--private-min-frac", str(cfg.get("private_min_frac", 0.10)),
                "--private-min-cov", str(cfg.get("private_min_cov", 20)),
            ])
        if as_bool(cfg.get("figures", True), True):
            # Isoform architecture-map figures are the PRIMARY per-category figure.
            # They are built from the GTF isoform models + the diff table's
            # per_transcript_json, so they need no --zn-long.
            args.extend([
                "--arch-figs-dir", str(self.paths.zn_class_figs_arch),
                "--figs-per-category", str(int(cfg.get("figs_per_category", 10))),
            ])
            # The 2-panel per-sample stoichiometry figures additionally need --figs-dir (--zn-long
            # is already added above for the PRIVATE scan).
            if self.paths.zn_filtered_long.exists():
                args.extend(["--figs-dir", str(self.paths.zn_class_figs)])
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
        # discover_candidate_snps is embarrassingly parallel and (with the native prefilter) cheap
        # per shard, so give it its OWN, larger job budget than the memory-sensitive later substeps
        # (raising the shared geno_jobs would also inflate their peak RSS). snp_scan_jobs=0/unset ->
        # use the full thread budget.
        _snp_cfg = int(geno.get("snp_scan_jobs", 0) or 0)
        snp_jobs = max(1, min(self.top_threads if self.top_threads > 0 else 8,
                              _snp_cfg if _snp_cfg > 0 else (self.top_threads or 8)))
        sample_bams = [
            str(self._require_existing_file(self._modkit_input_bam(sample), f"genotype input BAM for sample {sample}"))
            for sample in self.samples
        ]
        self._require_existing_file(self.paths.classification_summary, "classification summary TSV")
        self.paths.genotype.mkdir(parents=True, exist_ok=True)

        # ---- Step 1: discover candidate SNPs (needs the full BAMs; already locus-restricted). ----
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
                    "--jobs", str(snp_jobs),
                    "--threads", str(self.top_threads or 8),
                    "--window-bp", str(int(geno.get("snp_scan_window_bp", 1_000_000))),
                    "--primary-only",
                    "--verbose",
                ],
                label="discover_candidate_snps",
            )

        # ---- Step 2: candidate mod sites (needs zn_long + candidate SNPs; no BAM scan). ----
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
        # context_key on a shared read). Lossless for snp_mod_assoc/
        # haplotype_mod_assoc and keeps the per-read mod table tractable on deep genome-wide
        # data (otherwise 100k+ sites -> molecule_mod_calls OOM). Toggle off to keep all sites.
        if as_bool(geno.get("mod_sites_require_snp_link", True), True):
            mod_site_args += ["--candidate-snps", str(self.paths.geno_candidate_snps)]
        self.run_python_script(
            "build_candidate_mod_sites.py",
            mod_site_args,
            label="build_candidate_mod_sites",
        )

        # ---- Step 3: pre-subset the BAMs to the candidate regions. ----
        # Every read that can contribute to a snp/mod/haplotype association overlaps at least one
        # candidate SNP or mod site, so the per-molecule scans below only ever need those reads.
        # Restricting them here is lossless for every genotype output and removes the genome-wide
        # off-target reads that otherwise force build_read_assignment_table to materialize an
        # all-reads table (~22 GB on disk / ~70 GB in pandas on deep runs).
        subset_cfg = geno.get("subset_bams", {})
        scan_bams = sample_bams
        read_assignments_path = self.paths.geno_read_assignments
        used_subset = False
        if as_bool(subset_cfg.get("enable", True), True):
            require_tools(["samtools"])
            fai = f"{self.reference_fa}.fai"
            self.run_python_script(
                "build_candidate_regions_bed.py",
                [
                    "--candidate-snps", str(self.paths.geno_candidate_snps),
                    "--candidate-mod-bed", str(self.paths.geno_candidate_mod_bed),
                    "--fai", fai if Path(fai).exists() else "",
                    "--pad", str(int(subset_cfg.get("pad", 0))),
                    "--out-bed", str(self.paths.geno_candidate_regions_bed),
                ],
                label="build_candidate_regions_bed",
            )
            bed = self.paths.geno_candidate_regions_bed
            if self._nonempty(bed):
                self.paths.geno_subset_dir.mkdir(parents=True, exist_ok=True)
                sam_threads = max(1, min(4, self.top_threads or 1))
                # Per-fragmentform depth cap (0 = off). Highly-expressed loci pile hundreds of
                # thousands of reads on a few isoforms, which makes the per-window modkit extract in
                # build_molecule_mod_table saturate memory; capping each ZT fragmentform to max_ff
                # reads (seeded) collapses those isoforms to a bounded size while keeping every isoform
                # represented. Deterministic; results are a fixed-seed subsample (not full-depth).
                max_ff = int(subset_cfg.get("max_reads_per_fragmentform", 0) or 0)
                ff_seed = int(subset_cfg.get("subsample_seed", 12345))
                cap_script = str(self.script_path("cap_reads_per_fragmentform.py"))

                def _subset(in_bam: str) -> None:
                    out_bam = self.paths.geno_subset_bam(Path(in_bam))
                    if max_ff > 0:
                        regions_bam = f"{out_bam}.regions.bam"
                        run_command(
                            ["samtools", "view", "-b", "-M", "-L", str(bed), "-@", str(sam_threads),
                             str(in_bam), "-o", regions_bam],
                            cwd=self.root, label=f"samtools_subset[{out_bam.name}]", verbose=self.verbose,
                        )
                        run_command(
                            [sys.executable, cap_script, "--in-bam", regions_bam, "--out-bam", str(out_bam),
                             "--tag", "ZT", "--max-per-tag", str(max_ff), "--seed", str(ff_seed),
                             "--threads", str(sam_threads), "--verbose"],
                            cwd=self.root, label=f"cap_fragmentform[{out_bam.name}]", verbose=self.verbose,
                        )
                        try:
                            os.remove(regions_bam)
                        except OSError:
                            pass
                    else:
                        run_command(
                            ["samtools", "view", "-b", "-M", "-L", str(bed), "-@", str(sam_threads),
                             str(in_bam), "-o", str(out_bam)],
                            cwd=self.root, label=f"samtools_subset[{out_bam.name}]", verbose=self.verbose,
                        )
                        run_command(["samtools", "index", "-@", str(sam_threads), str(out_bam)],
                                    cwd=self.root, label=f"samtools_index[{out_bam.name}]", verbose=self.verbose)

                run_parallel(
                    [(f"subset[{Path(b).name}]", (lambda b=b: _subset(b))) for b in sample_bams],
                    jobs=self.jobs,
                )
                scan_bams = [str(self.paths.geno_subset_bam(Path(b))) for b in sample_bams]
                read_assignments_path = self.paths.geno_read_assignments_regions
                used_subset = True
            elif self.verbose:
                print("[modulator]   genotype: no candidate regions -- skipping BAM subsetting, "
                      "scanning the full BAMs.", flush=True)
        if self.verbose:
            print(f"[modulator]   genotype: per-molecule scans use "
                  f"{'candidate-region subset BAMs' if used_subset else 'full BAMs'}", flush=True)

        # ---- Step 4: per-molecule tables, over the (subset) BAMs. ----
        if not self._geno_reuse(read_assignments_path, "build_read_assignment_table"):
            self.run_python_script(
                "build_read_assignment_table.py",
                [
                    "--bams", *scan_bams,
                    "--summary-tsv", str(self.paths.classification_summary),
                    "--out-tsv", str(read_assignments_path),
                    "--jobs", str(geno_jobs),
                    "--primary-only",
                    "--verbose",
                ],
                label="build_read_assignment_table",
            )
        if not self._geno_reuse(self.paths.geno_molecule_snps, "build_molecule_snp_table"):
            self.run_python_script(
                "build_molecule_snp_table.py",
                [
                    "--bams", *scan_bams,
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
        mod_args = [
            # Subset BAMs: modkit only ever emits calls at --candidate-bed positions, and the
            # read-assignments table it joins against now covers exactly the candidate-region
            # reads -- so this is lossless, and modkit streams far fewer reads.
            "--bams", *scan_bams,
            "--candidate-sites-tsv", str(self.paths.geno_candidate_mod_sites),
            "--candidate-bed", str(self.paths.geno_candidate_mod_bed),
            "--read-assignments", str(read_assignments_path),
            "--reference-fa", str(self.reference_fa),
            "--out-tsv", str(self.paths.geno_molecule_mod_calls),
            "--threads", str(max(1, self.top_threads)),
            # A+B: this substep now streams shards to disk (bounded parent RSS) and shards per
            # (BAM x site-window), so heavy chromosomes no longer serialize and it can use far
            # more concurrency than the memory-bound per-molecule SNP steps -- decouple from
            # geno_jobs and cap at the thread budget. window-bp splits big chroms; interval-size
            # bounds each modkit process's own RSS.
            "--jobs", str(max(1, min(self.top_threads, int(geno.get("mod_jobs", 8))))),
            "--window-bp", str(int(geno.get("mod_scan_window_bp", 1_000_000))),
            "--interval-size", str(int(geno.get("mod_interval_size", 20000))),
            "--verbose",
        ]
        # Fast path: if `modkit extract calls` was pre-computed once per subset BAM (a separate
        # per-sample sbatch array across nodes, ~30-60 min/sample vs ~20 h windowed-on-one-node),
        # parse those TSVs instead of re-running modkit. Identical output; resumable per sample.
        mod_calls_dir = self.paths.genotype / "mod_calls"

        def _sample_of(bam_path):
            base = os.path.basename(str(bam_path))
            for suf in (".zt_tagged.clean.bam", ".zt_tagged.bam", ".bam"):
                if base.endswith(suf):
                    return base[: -len(suf)]
            return os.path.splitext(base)[0]

        pre = []
        for b in scan_bams:
            s = _sample_of(b)
            cp = mod_calls_dir / f"{s}.calls.tsv.bgz"
            if self._nonempty(cp):
                pre.append(f"{s}={cp}")
        if pre and len(pre) == len(scan_bams):
            mod_args += ["--pre-extracted", *pre]
            if self.verbose:
                print(f"[modulator] build_molecule_mod_table: using {len(pre)} pre-extracted "
                      f"per-sample call set(s) from {mod_calls_dir}", flush=True)
        else:
            # Default backend: built-in pysam streaming reader (per BAM x chromosome), ~100MB RSS,
            # never OOMs (chr15 included), no modkit / reference needed. Validated equivalent to modkit
            # (identical row set + call_code/target_modified).
            mod_args += ["--pysam"]
            if self.verbose:
                print("[modulator] build_molecule_mod_table: pysam streaming backend", flush=True)
        self.run_python_script("build_molecule_mod_table.py", mod_args, label="build_molecule_mod_table")
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
        # Why does a SNP change a modification? Positional ladder (at the modified base / inside the
        # DRACH 5-mer / 9-mer / proximal / distal cis), m6A motif disruption, and whether the observed
        # allelic direction matches the motif's prediction. Also flags SELF-REPORTING variants: SNPs are
        # called from RNA reads, so A-to-I editing (reads as G) and pseudouridine (U-to-C basecall error)
        # get called as variants at their own modified base, making those associations circular.
        if as_bool(geno.get("snp_mod_mechanism", True), True):
            self.run_python_script(
                "classify_snp_mod_mechanism.py",
                [
                    "--snp-mod-assoc", str(self.paths.geno_snp_mod),
                    "--reference-fa", str(self._require_reference_fa()),
                    "--out-tsv", str(self.paths.geno_snp_mod_mechanism),
                    "--proximal-bp", str(int(geno.get("mechanism_proximal_bp", 50))),
                    "--verbose",
                ],
                label="classify_snp_mod_mechanism",
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
        # Co-localized modifications: the mod x mod analogue of snp_mod_assoc -- do two nearby
        # mod sites co-occur on the same molecule more/less than expected? Reuses the per-read
        # mod-call table, so it adds no BAM scanning.
        colo = geno.get("colocalized_mods", {})
        if as_bool(colo.get("enable", True), True):
            self.run_python_script(
                "test_mod_mod_assoc.py",
                [
                    "--molecule-mods", str(self.paths.geno_molecule_mod_calls),
                    "--out-tsv", str(self.paths.geno_mod_mod),
                    "--max-distance", str(int(colo.get("max_distance", 1000))),
                    "--min-pair-reads", str(int(colo.get("min_pair_reads", 8))),
                    "--min-state-reads", str(int(colo.get("min_state_reads", 4))),
                    "--max-sites-per-read", str(int(colo.get("max_sites_per_read", 200))),
                    "--test", str(geno.get("test", "auto")),
                    "--pseudocount", str(float(geno.get("pseudocount", 0.5))),
                ],
                label="test_mod_mod_assoc",
            )

        # SNPs at the modified base: enumerate every modified site that coincides with a
        # candidate SNP, classify it (self-reporting A-to-I / pseU, or modified-base ablation),
        # and FLAG those sites in the between-isoform differential table -- so a segregating
        # variant at the modified base can be scrutinised/excluded rather than silently
        # confounding the SNP-blind differential test. (No test recalibration.)
        if as_bool(geno.get("snp_at_mod_base", True), True) and self._nonempty(self.paths.zn_filtered_long):
            args = [
                "--candidate-snps", str(self.paths.geno_candidate_snps),
                "--mod-sites", str(self.paths.zn_filtered_long),
                "--out-tsv", str(self.paths.geno_snp_at_mod_base),
                "--verbose",
            ]
            if self._nonempty(self.paths.zn_diff_results):
                args += ["--annotate", str(self.paths.zn_diff_results)]
            self.run_python_script("find_snp_at_mod_base.py", args, label="find_snp_at_mod_base")

    def stage_apa_motifs(self) -> None:
        """Polyadenylation-signal check for every APA site (each fragmentform's TES): canonical
        AATAAA / variant hexamer and its distance, downstream U/GU-richness, and an internal-priming
        flag (no PAS + genomic A-rich downstream = a likely oligo-dT artifact rather than a real site)."""
        cfg = self.config.get("apa_motifs", {})
        if not as_bool(cfg.get("enable", True), True):
            return
        if not self._nonempty(self.paths.classification_summary):
            if self.verbose:
                print("[modulator] apa_motifs: no classification summary (needs assemble), skipping", flush=True)
            return
        self.run_python_script("check_apa_motifs.py", [
            "--classification-summary", str(self.paths.classification_summary),
            "--reference-fa", str(self._require_reference_fa()),
            "--out-tsv", str(self.paths.apa_motifs),
            "--upstream", str(int(cfg.get("upstream", 60))),
            "--downstream", str(int(cfg.get("downstream", 40))),
            "--pas-max-distance", str(int(cfg.get("pas_max_distance", 40))),
            "--internal-priming-a-frac", str(float(cfg.get("internal_priming_a_frac", 0.65))),
            "--internal-priming-window", str(int(cfg.get("internal_priming_window", 20))),
            "--verbose",
        ], label="check_apa_motifs")

    def stage_polya(self) -> None:
        """Poly(A) tail length as a first-class readout: per-read dorado pt:i tail length ->
        per-fragmentform distributions, differential tail length between fragmentforms of a gene,
        and tail length x modification. Reads the ZT-tagged BAMs (which retain pt:i); tail x mod
        also needs the genotype molecule_mod_calls table."""
        cfg = self.config.get("polya", {})
        if not as_bool(cfg.get("enable", True), True):
            return
        # Prefer the multigene-cleaned tagged BAMs (final assignment set); fall back to the tagged BAMs.
        bams = [self.paths.clean_bam(s) for s in self.samples]
        if not all(self._nonempty(b) for b in bams):
            bams = [self.paths.zt_tagged_bam(s) for s in self.samples]
        bams = [b for b in bams if self._nonempty(b)]
        if not bams:
            if self.verbose:
                print("[modulator] polya: no ZT-tagged BAMs found (needs the assemble stage), skipping", flush=True)
            return
        pjobs = int(cfg.get("jobs", 0)) or self.top_threads or 1  # 0 -> use `threads`
        # 1) per-read tail table (pt:i joined to fragmentform / gene / metagene)
        args = ["--bams", *[str(b) for b in bams],
                "--out-tsv", str(self.paths.polya_read_tails),
                "--jobs", str(pjobs)]
        if self._nonempty(self.paths.classification_summary):
            args += ["--summary-tsv", str(self.paths.classification_summary)]
        if as_bool(self.config.get("assembler", {}).get("primary_only", True), True):
            args.append("--primary-only")
        self.run_python_script("build_read_polya_table.py", args, label="build_read_polya_table")
        # 2) per-fragmentform distributions + between-fragmentform differential tail length
        top_figs = int(cfg.get("top_figures", 10))
        self.run_python_script("test_taillength_diffs.py", [
            "--in-tsv", str(self.paths.polya_read_tails),
            "--out-prefix", str(self.paths.polya / self.prefix),
            "--min-reads", str(int(cfg.get("min_fragmentform_reads", 10))),
            "--min-total-reads", str(int(cfg.get("min_total_reads", 20))),
            "--min-tail", str(int(cfg.get("min_tail", 1))),
            "--figs-dir", str(self.paths.polya_diff_figs),
            "--top-k", str(top_figs),
        ], label="test_taillength_diffs")
        # 3) tail length x modification (needs the genotype per-read mod-call table)
        if self._nonempty(self.paths.geno_molecule_mod_calls):
            self.run_python_script("test_taillength_mod.py", [
                "--tail-tsv", str(self.paths.polya_read_tails),
                "--mod-tsv", str(self.paths.geno_molecule_mod_calls),
                "--out-tsv", str(self.paths.polya_taillength_mod),
                "--min-state-reads", str(int(cfg.get("min_state_reads", 10))),
                "--min-tail", str(int(cfg.get("min_tail", 1))),
                "--figs-dir", str(self.paths.polya_mod_figs),
                "--top-k", str(top_figs),
            ], label="test_taillength_mod")
        elif self.verbose:
            print("[modulator] polya: no molecule_mod_calls (genotype disabled/not run); "
                  "skipping tail x modification", flush=True)

    def stage_hierarchical_stoich(self) -> None:
        """Truncation-aware differential stoichiometry between fragmentforms -- the 5' complement to
        test_diffs, NOT a replacement.

        Direct-RNA reads truncate at the 5' end, so a read assigned to a fragmentform it never
        reached carries no evidence about features there. This restricts each fragmentform pair to
        the reads that demonstrably span their divergence point. Measured genome-wide (3,086 genes,
        304k tests): pairs diverging <1kb from the 3' end lose reads in 0.2% of tests and change
        ZERO calls -- test_diffs already answers those correctly and far more cheaply -- while pairs
        diverging >20kb lose reads in ~55% of tests. Hence `min_divergence_from_3p`: run the
        expensive engine only where it can actually change the answer.

        Off by default: it needs the genotype stage's per-read tables and costs ~15 min genome-wide.
        """
        cfg = self.config.get("hierarchical_stoich", {})
        if not as_bool(cfg.get("enable", False), False):
            return
        ra = self.paths.geno_read_assignments_regions
        if not self._nonempty(ra):
            ra = self.paths.geno_read_assignments
        for path, what in ((ra, "read-assignment table"),
                           (self.paths.geno_molecule_mod_calls, "molecule mod-call table"),
                           (self.paths.out_gtf, "assembled GTF")):
            if not self._nonempty(path):
                if self.verbose:
                    print(f"[modulator] hierarchical_stoich: missing {what} (needs the genotype stage); skipping",
                          flush=True)
                return
        args = [
            "--read-assignments", str(ra),
            "--molecule-mods", str(self.paths.geno_molecule_mod_calls),
            "--gtf", str(self.paths.out_gtf),
            "--out-tsv", str(self.paths.hierarchical_stoich),
            "--min-informative-reads", str(int(cfg.get("min_informative_reads", 10))),
            "--min-state-reads", str(int(cfg.get("min_state_reads", 3))),
            "--max-fragmentforms-per-gene", str(int(cfg.get("max_fragmentforms_per_gene", 12))),
            "--min-divergence-from-3p", str(int(cfg.get("min_divergence_from_3p", 0))),
            "--test", str(self.config.get("genotype", {}).get("test", "auto")),
            "--pseudocount", str(float(self.config.get("genotype", {}).get("pseudocount", 0.5))),
            "--verbose",
        ]
        # Preselected sites keep this cheap -- default to the sites test_diffs already flagged.
        sites = cfg.get("sites", "auto")
        if is_set(sites) and str(sites) != "auto":
            args += ["--sites", str(resolve_path(self.root, str(sites)))]
        elif str(sites) == "auto" and self._nonempty(self.paths.zn_diff_results):
            args += ["--sites", str(self.paths.zn_diff_results)]
        if as_bool(cfg.get("also_naive", True), True):
            args.append("--also-naive")
        self.run_python_script("test_hierarchical_stoich.py", args, label="test_hierarchical_stoich")

    def stage_between_conditions(self) -> None:
        """Replicate-aware BETWEEN-CONDITION comparisons, for every configured contrast.

        Needs a samplesheet (it supplies sample -> condition) with >=2 levels. The count-based
        analyses -- modification, and isoform / APA / junction usage -- all share the beta-binomial
        LRT with dispersion shrinkage in diffstats.py; poly(A) tail length is continuous and is
        compared across replicate summaries with Welch. NOTHING here pools reads across replicates:
        with millions of reads and n=3 per group that is pseudoreplication (measured 62% false
        positives on simulated nulls). See diffstats.py for the model and its calibration.
        """
        cfg = self.config.get("between_conditions", {})
        if not as_bool(cfg.get("enable", True), True):
            return
        if not self.contrasts:
            if self.verbose:
                print("[modulator] between_conditions: no contrasts — needs a samplesheet with a "
                      "'condition' column and >=2 levels; skipping", flush=True)
            return
        if not self._nonempty(self.paths.sample_metadata):
            return
        min_grp = str(int(cfg.get("min_samples_per_group", 2)))
        common_stat = ["--prior-weight", str(float(cfg.get("prior_weight", 20.0))),
                       "--ref-df", str(int(cfg.get("ref_df", 10))),
                       "--site-weight", str(cfg.get("site_weight", "auto")),
                       "--min-samples-per-group", min_grp]
        mod_filter = [str(m) for m in (cfg.get("mod_filter") or [])]
        for c in self.contrasts:
            name = c["name"]
            common = ["--sample-metadata", str(self.paths.sample_metadata), "--column", c["column"],
                      "--test", c["test"], "--reference", c["reference"],
                      "--contrast-name", name, "--verbose"]
            # 1) differential modification (per-sample site counts from the ZN long table)
            if as_bool(cfg.get("mod_diffs", True), True) and self._nonempty(self.paths.zn_filtered_long):
                args = ["--in-tsv", str(self.paths.zn_filtered_long),
                        "--out-tsv", str(self.paths.cond_mod_diffs(name)),
                        "--min-cov", str(int(cfg.get("min_cov", 20))), *common, *common_stat]
                if mod_filter:
                    args += ["--mod-filter", *mod_filter]
                self.run_python_script("test_condition_mod_diffs.py", args,
                                       label=f"condition_mod_diffs:{name}")
                # Flag between-condition sites that sit on a segregating SNP at the modified base
                # (genotype confounder). Needs candidate SNPs from the genotype stage.
                if self._nonempty(self.paths.geno_candidate_snps) and self._nonempty(self.paths.cond_mod_diffs(name)):
                    self.run_python_script("find_snp_at_mod_base.py", [
                        "--candidate-snps", str(self.paths.geno_candidate_snps),
                        "--mod-sites", str(self.paths.zn_filtered_long),
                        "--out-tsv", str(self.paths.geno_snp_at_mod_base),
                        "--annotate", str(self.paths.cond_mod_diffs(name)),
                    ], label=f"flag_snp_at_mod_base:{name}")
            # 2) differential usage: isoform / APA site / splice junction (one engine, three maps)
            for feature in ("isoform", "apa", "junction"):
                if not as_bool(cfg.get(f"{feature}_usage", True), True):
                    continue
                if not self._nonempty(self.paths.tx_counts):
                    break
                args = ["--tx-counts", str(self.paths.tx_counts),
                        "--out-tsv", str(self.paths.cond_usage_diffs(name, feature)),
                        "--feature", feature,
                        "--min-gene-reads", str(int(cfg.get("min_gene_reads", 20))), *common, *common_stat]
                if feature == "apa":
                    if not self._nonempty(self.paths.classification_summary):
                        continue
                    args += ["--classification-summary", str(self.paths.classification_summary)]
                if feature == "junction":
                    if not self._nonempty(self.paths.splice_junctions):
                        continue
                    args += ["--splice-junctions", str(self.paths.splice_junctions)]
                self.run_python_script("test_condition_usage_diffs.py", args,
                                       label=f"condition_{feature}_usage:{name}")
            # 3) differential poly(A) tail length (continuous -> Welch across replicates)
            if as_bool(cfg.get("tail_diffs", True), True) and self._nonempty(self.paths.polya_read_tails):
                self.run_python_script("test_condition_tail_diffs.py", [
                    "--tail-tsv", str(self.paths.polya_read_tails),
                    "--out-tsv", str(self.paths.cond_tail_diffs(name)),
                    "--level", str(cfg.get("tail_level", "fragmentform")),
                    "--min-reads-per-sample", str(int(cfg.get("min_tail_reads_per_sample", 10))),
                    "--min-samples-per-group", min_grp,
                    *common,
                ], label=f"condition_tail_diffs:{name}")

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
            "--run-manifest", str(self.paths.results / f"{self.prefix}_run_manifest.txt"),
            "--sample-metadata", str(self.paths.sample_metadata) if self._nonempty(self.paths.sample_metadata) else "",
            "--max-diff-figs", str(int(report_cfg.get("max_diff_figs", 6))),
            "--top-transcripts", str(int(report_cfg.get("top_transcripts", 20))),
            "--top-genes", str(int(report_cfg.get("top_genes", 20))),
            "--zn-long", str(self.paths.zn_filtered_long) if self.paths.zn_filtered_long.exists() else "",
            "--zt-long", str(self.paths.zt_filtered_long) if self.paths.zt_filtered_long.exists() else "",
            "--diff-results", str(self.paths.zn_diff_results) if self.paths.zn_diff_results.exists() else "",
            "--diff-figs-dir", str(self.paths.zn_diff_figs) if self.paths.zn_diff_figs.exists() else "",
            "--classified-sites", str(self.paths.zn_site_classified) if self.paths.zn_site_classified.exists() else "",
            "--private-sites", str(self.paths.zn_site_private) if self.paths.zn_site_private.exists() else "",
            "--class-figs-dir", str(self.paths.zn_class_figs) if self.paths.zn_class_figs.exists() else "",
            "--arch-figs-dir", str(self.paths.zn_class_figs_arch) if self.paths.zn_class_figs_arch.exists() else "",
            "--max-class-figs-per-category", str(int(report_cfg.get("max_class_figs_per_category", 10))),
            "--multigene-summary-glob", str(self.paths.zt_scrap_dir / "*.multigene_filter_summary.tsv") if self.paths.zt_scrap_dir.exists() else "",
            "--splice-junctions", str(self.paths.splice_junctions) if self.paths.splice_junctions.exists() else "",
            "--splice-genes", str(self.paths.gene_splice_summary) if self.paths.gene_splice_summary.exists() else "",
            "--novel-loci", str(self.paths.novel_loci_tsv) if self.paths.novel_loci_tsv.exists() else "",
            "--novel-fragmentforms", str(self.paths.novel_fragmentforms) if self.paths.novel_fragmentforms.exists() else "",
            "--mod-mod-assoc", str(self.paths.geno_mod_mod) if self.paths.geno_mod_mod.exists() else "",
            "--candidate-snps", str(self.paths.geno_candidate_snps) if self.paths.geno_candidate_snps.exists() else "",
            "--snp-tx-assoc", str(self.paths.geno_snp_tx) if self.paths.geno_snp_tx.exists() else "",
            "--snp-mod-assoc", str(self.paths.geno_snp_mod) if self.paths.geno_snp_mod.exists() else "",
            "--assembled-gtf", str(self.paths.out_gtf) if self.paths.out_gtf.exists() else "",
            "--molecule-mod-calls", str(self.paths.geno_molecule_mod_calls) if self.paths.geno_molecule_mod_calls.exists() else "",
            "--hap-blocks", str(self.paths.geno_hap_blocks) if self.paths.geno_hap_blocks.exists() else "",
            "--hap-tx-assoc", str(self.paths.geno_hap_tx) if self.paths.geno_hap_tx.exists() else "",
            "--hap-mod-assoc", str(self.paths.geno_hap_mod) if self.paths.geno_hap_mod.exists() else "",
            "--between-conditions-dir", str(self.paths.between_conditions) if self.paths.between_conditions.is_dir() else "",
            "--apa-motifs", str(self.paths.apa_motifs) if self.paths.apa_motifs.exists() else "",
            "--sequence-elements", str(self.paths.sequence_elements) if self.paths.sequence_elements.exists() else "",
            "--sequence-elements-summary", str(self.paths.sequence_elements_summary) if self.paths.sequence_elements_summary.exists() else "",
            "--snp-mod-mechanism", str(self.paths.geno_snp_mod_mechanism) if self.paths.geno_snp_mod_mechanism.exists() else "",
            "--polya-fragmentform", str(self.paths.polya_fragmentform) if self.paths.polya_fragmentform.exists() else "",
            "--taillength-diffs", str(self.paths.polya_taillength_diffs) if self.paths.polya_taillength_diffs.exists() else "",
            "--taillength-mod", str(self.paths.polya_taillength_mod) if self.paths.polya_taillength_mod.exists() else "",
            "--taillength-diff-figs", str(self.paths.polya_diff_figs) if self.paths.polya_diff_figs.is_dir() else "",
            "--taillength-mod-figs", str(self.paths.polya_mod_figs) if self.paths.polya_mod_figs.is_dir() else "",
            "--snp-figs-dir", str(self.paths.genotype / f"{self.prefix}__snp_figs"),
            "--max-snp-figs", str(int(report_cfg.get("max_snp_figs", 12))),
        ]
        self.run_python_script("generate_html_report.py", args, label="generate_html_report")

        # Companion interactive browser: search a gene/fragmentform, click an exon to filter its
        # modification sites and differential results. Self-contained HTML, built from the same tables.
        if as_bool(report_cfg.get("gene_browser", True), True):
            gb = [
                "--gtf", str(self.paths.out_gtf),
                "--out-html", str(self.paths.gene_browser_html),
                "--title", f"{self.prefix} — gene browser",
                "--max-genes", str(int(report_cfg.get("browser_max_genes", 4000))),
                "--verbose",
            ]
            for flag, path in (
                ("--sites-long", self.paths.zn_filtered_long),
                ("--diff-results", self.paths.zn_diff_results),
                ("--classification-summary", self.paths.classification_summary),
                ("--apa-motifs", self.paths.apa_motifs),
                ("--polya-fragmentform", self.paths.polya_fragmentform),
                ("--hierarchical-stoich", self.paths.hierarchical_stoich),
            ):
                if self._nonempty(path):
                    gb += [flag, str(path)]
            cm = sorted(self.paths.between_conditions.glob(f"{self.prefix}_*_mod_diffs.tsv")) \
                if self.paths.between_conditions.is_dir() else []
            if cm:
                gb += ["--condition-mod-diffs", str(cm[0])]
            if self._nonempty(self.paths.out_gtf):
                self.run_python_script("build_gene_browser.py", gb, label="build_gene_browser")
