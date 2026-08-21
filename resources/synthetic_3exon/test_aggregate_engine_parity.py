#!/usr/bin/env python3
"""
Guard the stream-vs-sort aggregation engines against the FILTERED-output divergence.

Both engines (aggregate_zn_stream.py = default, aggregate_by_gene.py = sort) must produce the same
FILTERED_sites_long table for a given config. The bug this guards: when site-filtering is turned OFF,
the sort engine writes FILTERED == RAW while the stream engine used to write NO filtered table at all,
which silently skipped every downstream stage (test_diffs, ...) that keys on the FILTERED table.

This is an INTEGRATION test that reuses the demo's real ZN bedMethyl input (bgzipped + tabixed), so it
is skipped (exit 0) when that input is not present.

Usage: <modulator-env>/bin/python resources/synthetic_3exon/test_aggregate_engine_parity.py
"""
from __future__ import annotations
import os
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = ROOT / "workflow" / "scripts"
# demo run sits next to the repo checkout
DEMO = ROOT.parent / "demo14genes_run"
MODKIT = DEMO / "results" / "modkit_zn"
GTF = DEMO / "results" / "assemble" / "demo14.gtf"


def _rows(long_tsv):
    """Return the FILTERED long table as a sorted set of data lines (order-agnostic comparison)."""
    if not long_tsv.exists():
        return None
    with open(long_tsv) as fh:
        lines = fh.read().splitlines()
    return set(lines[1:]) if lines else set()


def run_engine(script, out_prefix, filter_on):
    args = [sys.executable, str(SCRIPTS / script),
            "--modkit-dir", str(MODKIT), "--gtf", str(GTF), "--out-prefix", str(out_prefix),
            "--min-cov", "0", "--count-diff-factor", "3.0", "--mod-fail-margin", "1",
            "--nfail-score-k", "0.0", "--emit-raw", "--emit-filtered", "--write-long"]
    if "stream" in script:
        args += ["--jobs", "4"]   # the sort engine has no --jobs flag
    if filter_on:
        args.append("--filter-enable")
    subprocess.run(args, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return Path(f"{out_prefix}_FILTERED_sites_long.tsv"), Path(f"{out_prefix}_RAW_sites_long.tsv")


def main():
    if not (MODKIT.is_dir() and GTF.exists()):
        print(f"  SKIP: demo ZN bedMethyl input not found at {MODKIT} — run the demo first.")
        sys.exit(0)

    checks = []

    def check(name, ok):
        checks.append(ok)
        print(f"  {'PASS' if ok else '**FAIL**'}  {name}")

    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        # filter ON: the two engines must agree
        s_filt_on, _ = run_engine("aggregate_zn_stream.py", td / "stream_on", filter_on=True)
        g_filt_on, _ = run_engine("aggregate_by_gene.py", td / "sort_on", filter_on=True)
        so, go = _rows(s_filt_on), _rows(g_filt_on)
        check("filter ON: stream FILTERED long is written", so is not None)
        check("filter ON: sort FILTERED long is written", go is not None)
        check("filter ON: stream FILTERED == sort FILTERED (engine parity)", so == go)

        # filter OFF: both engines must still write FILTERED, and it must equal RAW (== all sites)
        s_filt_off, s_raw_off = run_engine("aggregate_zn_stream.py", td / "stream_off", filter_on=False)
        g_filt_off, g_raw_off = run_engine("aggregate_by_gene.py", td / "sort_off", filter_on=False)
        s_off, s_raw = _rows(s_filt_off), _rows(s_raw_off)
        g_off, g_raw = _rows(g_filt_off), _rows(g_raw_off)
        check("filter OFF: stream FILTERED long is PRESENT (was silently skipped)", s_off is not None)
        check("filter OFF: stream FILTERED == RAW (filtered == everything)", s_off == s_raw)
        check("filter OFF: sort FILTERED == RAW", g_off == g_raw)
        check("filter OFF: stream FILTERED == sort FILTERED (engine parity)", s_off == g_off)
        check("filter OFF: FILTERED is non-empty (downstream test_diffs not skipped)", bool(s_off))

    n_fail = checks.count(False)
    print(f"\naggregate engine parity: {len(checks) - n_fail}/{len(checks)} checks passed"
          + ("" if not n_fail else f"  ({n_fail} FAILED)"))
    sys.exit(1 if n_fail else 0)


if __name__ == "__main__":
    main()
