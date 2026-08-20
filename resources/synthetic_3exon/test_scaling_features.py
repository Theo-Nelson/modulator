#!/usr/bin/env python3
"""
Guard the two v2.0.0 scale-up features, directly (no full pipeline run needed):

  A. aggregate_zn pivot tri-state (aggregate_by_gene.generate_per_gene_outputs_from_dedup):
     per-gene pivots are optional inspection outputs (3 dense files per gene x mod group, read by
     nothing downstream). 'auto' must write them below pivot_max_groups and SKIP them above it
     (avoiding a tiny-file explosion on whole-transcriptome runs) while STILL writing the per-gene
     long tables; 'on' must force them even past the ceiling; 'off' must never write them.

  B. BAM preflight (pipeline.ModulatorPipeline._preflight_bams / _bam_has_index): fail fast on
     duplicate BAMs (two samples -> same physical file = pseudo-replication) and on a missing
     .bai/.csi index (modkit/pysam require one), and pass a clean, unique, indexed set.

Usage: <modulator-env>/bin/python resources/synthetic_3exon/test_scaling_features.py
"""
from __future__ import annotations

import glob
import os
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "workflow" / "scripts"))
sys.path.insert(0, str(REPO / "src"))

import aggregate_by_gene as agg  # noqa: E402


# One dedup row = the 18 tab-columns generate_per_gene_outputs_from_dedup expects, in order:
# sample, ZN, chrom, start0, end0, strand, mod, gid, gname, cov, Nmod, Ncan, Nother, Ndel, Nfail,
# Ndiff, Nnocall, frac
def _dedup_row(sample, zn, chrom, start0, mod, gid, gname):
    return "\t".join(str(x) for x in [
        sample, zn, chrom, start0, start0 + 1, "+", mod, gid, gname,
        50, 30, 20, 0, 0, 0, 0, 0, 0.6,
    ])


def _write_dedup(path, n_genes, mods=("a", "17802"), samples=("S1", "S2")):
    """n_genes x len(mods) groups; each group one site x len(samples) rows."""
    with open(path, "w") as f:
        for g in range(n_genes):
            for mi, mod in enumerate(mods):
                for s in samples:
                    f.write(_dedup_row(s, 1, "chr1", 1000 + g * 10 + mi, mod,
                                       f"GID{g}", f"GENE{g}") + "\n")
    return n_genes * len(mods)


def _count_pivots(out_prefix, tag):
    return len(glob.glob(f"{out_prefix}_{tag}__per_gene_mod/*_pivot.tsv"))


def _count_pergene(out_prefix, tag):
    d = f"{out_prefix}_{tag}__per_gene_mod"
    return len([p for p in glob.glob(f"{d}/*.tsv") if not p.endswith("_pivot.tsv")])


def _run(td, tag, pivot_mode, pivot_max_groups, n_genes=10, write_per_gene=True):
    dedup = td / f"dedup.{tag}.tsv"
    n_groups = _write_dedup(dedup, n_genes)
    base = str(td / f"out.{tag}")
    work = td / f"work.{tag}"
    work.mkdir(exist_ok=True)
    agg.generate_per_gene_outputs_from_dedup(
        str(dedup), base, tag, write_per_gene, pivot_mode, str(work),
        chunk_lines=100000, verbose=False, jobs=1, pivot_max_groups=pivot_max_groups,
    )
    return n_groups, _count_pivots(base, tag), _count_pergene(base, tag)


def test_pivots():
    checks = []
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        # 10 genes x 2 mods = 20 groups. Each group writes 3 pivots when enabled.
        ng, piv, pg = _run(td, "AUTO_ON", "auto", pivot_max_groups=1000)
        checks.append((f"auto BELOW ceiling -> pivots ON ({ng} groups <= 1000)", piv == 3 * ng))
        checks.append(("auto ON also writes per-gene tables", pg == ng))

        ng, piv, pg = _run(td, "AUTO_OFF", "auto", pivot_max_groups=5)
        checks.append((f"auto ABOVE ceiling -> pivots OFF ({ng} groups > 5)", piv == 0))
        checks.append(("auto OFF STILL writes per-gene long tables (analysis intact)", pg == ng))

        ng, piv, _ = _run(td, "FORCE_ON", "on", pivot_max_groups=5)
        checks.append((f"'on' forces pivots even past ceiling ({ng} groups > 5)", piv == 3 * ng))

        ng, piv, pg = _run(td, "FORCE_OFF", "off", pivot_max_groups=10_000)
        checks.append(("'off' never writes pivots (below ceiling)", piv == 0))
        checks.append(("'off' still writes per-gene tables when requested", pg == ng))

        ng, piv, _ = _run(td, "LEGACY_TRUE", True, pivot_max_groups=5)
        checks.append(("legacy bool True maps to 'on' (writes at scale)", piv == 3 * ng))
        ng, piv, _ = _run(td, "LEGACY_FALSE", False, pivot_max_groups=10_000)
        checks.append(("legacy bool False maps to 'off'", piv == 0))
    return checks


def test_preflight():
    """Drive ModulatorPipeline._preflight_bams via a tiny stub (no real BAM decoding needed --
    the check only stats files + looks for a sibling .bai/.csi)."""
    from modulator import pipeline as pl  # noqa: E402

    class Stub:
        # borrow the real methods unbound
        _preflight_bams = pl.ModulatorPipeline._preflight_bams
        _bam_has_index = staticmethod(pl.ModulatorPipeline._bam_has_index)

        def __init__(self, bams_dir, cfg):
            self._bams_dir = Path(bams_dir)
            self.config = cfg
            self.verbose = False
            self.top_threads = 1
            self.root = Path(bams_dir)

        @property
        def bams_dir(self):
            return self._bams_dir

        @property
        def bam_glob(self):
            return "*.bam"

    def make_bam(p, content=b"BAMSTUB", index=True):
        p.write_bytes(content)
        if index:
            Path(str(p) + ".bai").write_bytes(b"IDX")

    checks = []
    cfg = {"preflight": {"enable": True, "check_duplicate_bams": True, "check_bam_index": True}}

    # 1. clean, unique, indexed -> passes
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        make_bam(td / "M1.bam", b"AAA", index=True)
        make_bam(td / "M2.bam", b"BBB", index=True)
        try:
            Stub(td, cfg)._preflight_bams()
            checks.append(("clean unique+indexed set passes preflight", True))
        except Exception as e:  # noqa: BLE001
            checks.append((f"clean set passes ({e})", False))

    # 2. missing index -> raises, message names the fix
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        make_bam(td / "M1.bam", b"AAA", index=True)
        make_bam(td / "M2.bam", b"BBB", index=False)
        try:
            Stub(td, cfg)._preflight_bams()
            checks.append(("missing-index BAM is rejected", False))
        except ValueError as e:
            checks.append(("missing-index BAM is rejected", "index" in str(e).lower()
                           and "samtools index" in str(e)))

    # 3. duplicate physical file (hardlink) -> raises pseudo-replication error
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        make_bam(td / "A.bam", b"SAME", index=True)
        os.link(td / "A.bam", td / "B.bam")          # same inode -> duplicate
        Path(str(td / "B.bam") + ".bai").write_bytes(b"IDX")
        try:
            Stub(td, cfg)._preflight_bams()
            checks.append(("duplicate (same physical file) is rejected", False))
        except ValueError as e:
            checks.append(("duplicate (same physical file) is rejected",
                           "same physical file" in str(e).lower()))

    # 4. disabled -> skips even a broken set
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        make_bam(td / "M1.bam", b"AAA", index=False)  # unindexed
        try:
            Stub(td, {"preflight": {"enable": False}})._preflight_bams()
            checks.append(("preflight.enable:false skips the checks", True))
        except Exception:  # noqa: BLE001
            checks.append(("preflight.enable:false skips the checks", False))

    return checks


def test_nfail_score_k():
    """NFail-SCORE k-ratio site filter (aggregate_by_gene.row_pass_filter): keep a site iff
    k = Nmod/(Nfail+1) >= nfail_score_k, on top of the Ndiff and mod_fail_margin guards."""
    import aggregate_by_gene as agg
    P = agg.row_pass_filter
    checks = [
        ("k=0 disables the k-ratio filter (legacy behaviour)", P(100, 10, 4, 0, 3, 1, 0.0) is True),
        ("clean site passes at k=0.4 (k=10/5=2.0)", P(100, 10, 4, 0, 3, 1, 0.4) is True),
        ("error-prone site (Nmod=3,Nfail=40 -> k=0.07) fails at k=0.4", P(100, 3, 40, 0, 3, 1, 0.4) is False),
        ("stricter k=1.5 rejects a borderline site (Nmod=6,Nfail=4 -> k=1.2)", P(100, 6, 4, 0, 3, 1, 1.5) is False),
        ("Ndiff guard still applies regardless of k", P(100, 50, 0, 400, 3, 1, 0.4) is False),
    ]
    return checks


def test_stage_skip_reason():
    """A disabled/no-op stage must be reported as disabled, not "checkpoint found" (pipeline
    ._stage_disabled), so a --resume log never claims work was reused when a stage was simply off."""
    from modulator import pipeline as pl

    class Stub:
        _stage_disabled = pl.ModulatorPipeline._stage_disabled

        def __init__(self, config, contrasts):
            self.config = config
            self.contrasts = contrasts

    off = Stub({"genotype": {"enable": False}}, contrasts=[])
    on = Stub({"between_conditions": {"enable": True}}, contrasts=[{"name": "a_vs_b"}])
    return [
        ("genotype off -> disabled", off._stage_disabled("genotype") is True),
        ("hierarchical_stoich (default off) -> disabled", off._stage_disabled("hierarchical_stoich") is True),
        ("between_conditions with no contrasts -> disabled", off._stage_disabled("between_conditions") is True),
        ("between_conditions WITH contrasts -> not disabled", on._stage_disabled("between_conditions") is False),
        ("assemble (always on) -> not disabled", off._stage_disabled("assemble") is False),
    ]


def main():
    all_checks = [("PIVOT TRI-STATE", test_pivots()), ("BAM PREFLIGHT", test_preflight()),
                  ("NFAIL-SCORE k-RATIO", test_nfail_score_k()),
                  ("STAGE SKIP REASON", test_stage_skip_reason())]
    n_fail = 0
    for header, checks in all_checks:
        print(f"== {header} ==")
        for msg, ok in checks:
            print(f"  [{'PASS' if ok else '**FAIL**'}] {msg}")
            n_fail += (not ok)
    total = sum(len(c) for _, c in all_checks)
    print(f"\nscaling features: {total - n_fail}/{total} checks passed"
          + ("" if not n_fail else f"  ({n_fail} FAILED)"))
    sys.exit(1 if n_fail else 0)


if __name__ == "__main__":
    main()
