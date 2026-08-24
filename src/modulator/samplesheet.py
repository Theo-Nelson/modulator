"""Samplesheet: the sample source and the metadata for between-condition comparisons.

A TSV (or CSV) with one row per sample::

    sample   bam                             condition   replicate
    M1       HornerLab_M1pA_...sorted.bam    mock        1
    Z1       HornerLab_Z1pA_...sorted.bam    zikv        1

``sample`` and ``bam`` are required; ``condition`` is required to run any between-condition test.
Any further columns are carried through as covariates.

WHY STAGING: the pipeline names samples by the BAM stem (``sample_name_from_bam``), and that name
flows into every output table, so honouring a friendly ``sample`` id would otherwise mean touching
every script that derives a sample name. Instead each BAM is SYMLINKED to ``<staging>/<sample>.bam``
(plus its index) and the pipeline points ``bams_dir``/``bam_glob`` at the staging directory. The stem
then *is* the friendly id, everywhere, with no script changes and no copying.

stdlib only (the package's sole dependency is PyYAML).
"""
from __future__ import annotations

import csv
import glob as _glob
import os
from pathlib import Path

REQUIRED_COLUMNS = ("sample", "bam")
_BAM_INDEX_SUFFIXES = (".bam.bai", ".bai")


class SamplesheetError(ValueError):
    """Raised for a malformed samplesheet -- always names the offending row/value."""


def _sniff_delimiter(path: Path) -> str:
    with open(path, newline="") as fh:
        first = fh.readline()
    return "\t" if "\t" in first else ","


def read_samplesheet(path: str | Path) -> list[dict]:
    """Parse + validate a samplesheet. Returns a list of row dicts (order preserved)."""
    path = Path(path)
    if not path.exists():
        raise SamplesheetError(f"samplesheet not found: {path}")
    rows: list[dict] = []
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter=_sniff_delimiter(path))
        if reader.fieldnames is None:
            raise SamplesheetError(f"samplesheet {path} is empty")
        fields = [(f or "").strip().lstrip("#") for f in reader.fieldnames]
        missing = [c for c in REQUIRED_COLUMNS if c not in fields]
        if missing:
            raise SamplesheetError(
                f"samplesheet {path} is missing required column(s): {', '.join(missing)}. "
                f"Found: {', '.join(fields)}")
        for i, raw in enumerate(reader, start=2):  # header is line 1
            row = {(k or "").strip().lstrip("#"): (v.strip() if isinstance(v, str) else v)
                   for k, v in raw.items() if k is not None}
            if not any((v or "") for v in row.values()):
                continue  # blank line
            sample = row.get("sample") or ""
            if not sample:
                raise SamplesheetError(f"{path} line {i}: empty 'sample'")
            if any(ch in sample for ch in "/\\ \t"):
                raise SamplesheetError(
                    f"{path} line {i}: sample id {sample!r} must not contain spaces or path separators "
                    "(it becomes a filename)")
            if not (row.get("bam") or ""):
                raise SamplesheetError(f"{path} line {i}: empty 'bam' for sample {sample!r}")
            rows.append(row)
    if not rows:
        raise SamplesheetError(f"samplesheet {path} has no data rows")
    seen: dict[str, int] = {}
    for idx, r in enumerate(rows):
        s = r["sample"]
        if s in seen:
            raise SamplesheetError(f"{path}: duplicate sample id {s!r}")
        seen[s] = idx
    return rows


def resolve_bam(entry: str, bams_dir: Path) -> Path:
    """Resolve a samplesheet 'bam' value: absolute path, path relative to bams_dir, or a glob."""
    cand = Path(entry)
    if not cand.is_absolute():
        cand = Path(bams_dir) / entry
    if cand.exists():
        return cand.resolve()
    hits = sorted(_glob.glob(str(cand)))
    if len(hits) == 1:
        return Path(hits[0]).resolve()
    if not hits:
        raise SamplesheetError(f"no BAM matched {entry!r} (looked at {cand})")
    raise SamplesheetError(f"{entry!r} matched {len(hits)} BAMs; it must identify exactly one: "
                           f"{', '.join(os.path.basename(h) for h in hits[:4])}...")


def _find_index(bam: Path) -> Path | None:
    stem = str(bam.with_suffix(""))            # ".../x"  for  ".../x.bam"
    # accept both BAI and CSI (CSI is required for chromosomes > 512 Mbp), in either naming convention
    for cand in (Path(str(bam) + ".bai"), Path(stem + ".bai"),
                 Path(str(bam) + ".csi"), Path(stem + ".csi")):
        if cand.exists():
            return cand
    return None


def _link(src: Path, dst: Path) -> None:
    """Symlink src -> dst, replacing a stale link. Never copies (BAMs are huge)."""
    if dst.is_symlink() or dst.exists():
        try:
            if dst.is_symlink() and Path(os.readlink(dst)) == src:
                return
            dst.unlink()
        except OSError:
            return
    dst.symlink_to(src)


def stage_bams(rows: list[dict], staging_dir: str | Path, bams_dir: str | Path) -> list[str]:
    """Symlink each row's BAM (and index) to <staging_dir>/<sample>.bam. Returns the sample ids.

    This is what makes the samplesheet the sample SOURCE: downstream, the BAM stem is the sample id.
    """
    staging = Path(staging_dir)
    staging.mkdir(parents=True, exist_ok=True)
    samples: list[str] = []
    for r in rows:
        sample, bam = r["sample"], resolve_bam(r["bam"], Path(bams_dir))
        _link(bam, staging / f"{sample}.bam")
        idx = _find_index(bam)
        if idx is not None:
            # preserve the index type (.csi must not be linked as .bam.bai, or readers mis-handle it)
            ext = ".bam.csi" if idx.suffix == ".csi" else ".bam.bai"
            _link(idx, staging / f"{sample}{ext}")
        samples.append(sample)
    # Prune stale links from a PRIOR samplesheet: a re-run with a sample removed would otherwise leave
    # its <sample>.bam symlink behind, and anything that globs the staged dir would silently re-include
    # the dropped sample. Only remove symlinks (never real files a user may have dropped in) whose
    # sample stem is not in the current set.
    keep = set(samples)
    for link in staging.glob("*.bam*"):
        if not link.is_symlink():
            continue
        stem = link.name.split(".bam")[0]
        if stem not in keep:
            try:
                link.unlink()
            except OSError:
                pass
    return samples


def sample_metadata(rows: list[dict]) -> dict[str, dict]:
    """{sample: {condition: ..., <covariates>}} -- everything except the bam path."""
    return {r["sample"]: {k: v for k, v in r.items() if k not in ("sample", "bam")} for r in rows}


def write_metadata_tsv(rows: list[dict], out_path: str | Path) -> Path:
    """Write the resolved metadata for the differential scripts to consume."""
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    cols = ["sample"] + [c for c in rows[0].keys() if c not in ("sample", "bam")]
    with open(out_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    return out_path


def _safe_contrast_name(name: str) -> str:
    """Make a contrast name safe to embed in a filename.

    Contrast names become path components (``{prefix}_{name}_mod_diffs.tsv`` etc.), so a name with a
    ``/`` would silently redirect the output into a subdirectory (or fail), and stray whitespace makes
    brittle paths. Keep alnum / ``-`` / ``.``; collapse every other run to a single ``_``; trim
    leading/trailing ``._``. Collisions AFTER sanitizing are caught by the caller.
    """
    import re
    s = re.sub(r"[^A-Za-z0-9._-]+", "_", str(name).strip())
    s = s.strip("._")
    return s or "contrast"


def resolve_contrasts(config_contrasts, rows: list[dict], column: str = "condition") -> list[dict]:
    """Normalise the configured contrasts, or derive every pairwise one from ``column``.

    Each contrast -> {name, column, test, reference}. Validated against the samplesheet so a typo
    fails loudly here rather than silently producing an empty comparison. Names are sanitized to be
    filesystem-safe (they become output-filename components) and checked for post-sanitize collisions.
    """
    levels: list[str] = []
    for r in rows:
        v = (r.get(column) or "").strip()
        if v and v not in levels:
            levels.append(v)
    out: list[dict] = []
    if config_contrasts:
        for c in config_contrasts:
            if not isinstance(c, dict):
                raise SamplesheetError(f"contrast {c!r} must be a mapping with test/reference")
            col = c.get("column", column)
            test, ref = c.get("test"), c.get("reference")
            if not test or not ref:
                raise SamplesheetError(f"contrast {c!r} needs both 'test' and 'reference'")
            have = {(r.get(col) or "").strip() for r in rows}
            for side, val in (("test", test), ("reference", ref)):
                if val not in have:
                    raise SamplesheetError(
                        f"contrast {c.get('name', '')!r}: {side}={val!r} not found in samplesheet "
                        f"column {col!r} (levels: {', '.join(sorted(v for v in have if v))})")
            out.append({"name": _safe_contrast_name(c.get("name") or f"{test}_vs_{ref}"),
                        "column": col, "test": test, "reference": ref})
        _reject_name_collisions(out)
        return out
    for i, ref in enumerate(levels):          # no explicit contrasts -> all pairwise
        for test in levels[i + 1:]:
            out.append({"name": _safe_contrast_name(f"{test}_vs_{ref}"),
                        "column": column, "test": test, "reference": ref})
    _reject_name_collisions(out)
    return out


def _reject_name_collisions(contrasts: list[dict]) -> None:
    """Two contrasts whose sanitized names collide would overwrite each other's output tables."""
    seen: dict[str, dict] = {}
    for c in contrasts:
        n = c["name"]
        if n in seen:
            raise SamplesheetError(
                f"contrast name {n!r} is not unique after sanitizing (also used by "
                f"{seen[n]['test']}_vs_{seen[n]['reference']} vs {c['test']}_vs_{c['reference']}); "
                f"give each contrast a distinct 'name'")
        seen[n] = c
