#!/usr/bin/env python3

from concurrent.futures import ProcessPoolExecutor, as_completed
import gzip
import math
import os
import sys
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np
import pandas as pd
from scipy.stats import chi2, chi2_contingency, fisher_exact


def safe_int(x, default=0) -> int:
    try:
        return int(x)
    except Exception:
        try:
            return int(float(x))
        except Exception:
            return default


def safe_float(x, default=0.0) -> float:
    try:
        return float(x)
    except Exception:
        return default


def sample_name_from_bam(path: str) -> str:
    base = os.path.basename(path)
    for suffix in (".zt_tagged.clean.bam", ".zt_tagged.bam", ".bam"):
        if base.endswith(suffix):
            return base[: -len(suffix)]
    return os.path.splitext(base)[0]


def benjamini_hochberg(pvals: Iterable[float]) -> np.ndarray:
    p = np.asarray(list(pvals), dtype=float)
    if p.size == 0:
        return np.asarray([], dtype=float)
    # Rank only the finite p-values (statsmodels.multipletests semantics). A single NaN
    # would otherwise sort last and, via minimum.accumulate on the reversed array, poison
    # every adjusted p-value with NaN -> total silent loss of significance. NaN in stays NaN.
    out = np.full(p.size, np.nan, dtype=float)
    idx = np.flatnonzero(np.isfinite(p))
    m = idx.size
    if m == 0:
        return out
    pf = p[idx]
    order = np.argsort(pf)
    ranks = np.empty(m, dtype=int)
    ranks[order] = np.arange(1, m + 1)
    adj = pf * m / ranks
    adj_sorted = np.minimum.accumulate(adj[order][::-1])[::-1]
    adj_final = np.empty(m, dtype=float)
    adj_final[order] = adj_sorted
    out[idx] = np.clip(adj_final, 0.0, 1.0)
    return out


def max_abs_distribution_shift(table: np.ndarray) -> float:
    table = np.asarray(table, dtype=float)
    if table.size == 0 or table.shape[0] < 2 or table.shape[1] < 1:
        return 0.0
    row_sums = table.sum(axis=1, keepdims=True)
    frac = np.divide(table, row_sums, out=np.zeros_like(table, dtype=float), where=row_sums > 0)
    max_diff = 0.0
    for i in range(frac.shape[0]):
        for j in range(i + 1, frac.shape[0]):
            max_diff = max(max_diff, float(np.max(np.abs(frac[i] - frac[j]))))
    return round(max_diff, 6)


def binary_rate_delta(table_2x2: np.ndarray) -> float:
    table = np.asarray(table_2x2, dtype=float)
    if table.shape != (2, 2):
        return 0.0
    rate0 = table[0, 0] / table[0].sum() if table[0].sum() > 0 else 0.0
    rate1 = table[1, 0] / table[1].sum() if table[1].sum() > 0 else 0.0
    return round(float(abs(rate0 - rate1)), 6)


def run_contingency_test(
    table: np.ndarray,
    test: str = "auto",
    pseudocount: float = 0.5,
) -> Tuple[str, str, float, float]:
    tab = np.asarray(table, dtype=float)
    if tab.size == 0 or tab.shape[0] < 2 or tab.shape[1] < 2:
        return "none", "none", 0.0, 1.0

    # A table with a zero marginal (a group or an outcome with no reads at all -- e.g. 100%/100%,
    # 0%/0%, or a monomorphic group) is structurally UNTESTABLE: there is no variation to test on one
    # axis. fisher returns an undefined (nan) odds ratio and chi2_contingency RAISES when pseudocount=0.
    # Return NaN p so it is excluded from the BH family (it can never be significant), rather than
    # crashing, writing a nan/inf statistic, or padding the multiple-testing burden with p==1 rows.
    if not ((tab.sum(axis=0) > 0).all() and (tab.sum(axis=1) > 0).all()):
        return "untestable", "none", float("nan"), float("nan")

    def do_fisher_2x2(tt):
        odds, p = fisher_exact(tt.astype(int))
        # a single zero CELL (not a zero margin) is still testable; its odds ratio is genuinely
        # infinite -> keep inf. nan cannot occur here (degenerate margins handled above).
        odds = float(odds) if math.isfinite(odds) else float("inf")
        return "fisher_exact_2x2", "fisher_odds", odds, float(p)

    def do_chi2(tt):
        tt_pc = tt + float(pseudocount)
        stat, p, _, _ = chi2_contingency(tt_pc, correction=False)
        return f"chi2_{tt.shape[0]}x{tt.shape[1]}_pc{pseudocount:g}", "chi2", float(stat), float(p)

    if test == "fisher":
        if tab.shape == (2, 2):
            return do_fisher_2x2(tab)
        return do_chi2(tab)
    if test == "chi2":
        return do_chi2(tab)
    if tab.shape == (2, 2):
        return do_fisher_2x2(tab)
    return do_chi2(tab)


def robust_load_summary(path: str) -> pd.DataFrame:
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        return pd.DataFrame()
    df = pd.read_csv(path, sep="\t", low_memory=False)
    df.columns = [str(c).lstrip("#") for c in df.columns]
    return df


def load_read_assignments(path: str) -> pd.DataFrame:
    df = robust_load_summary(path)
    if df.empty:
        return df
    need = ["sample", "qname"]
    missing = [c for c in need if c not in df.columns]
    if missing:
        raise SystemExit(f"Missing columns in read assignments table {path}: {missing}")
    return df


def normalize_text_token(value, *, numeric: bool = False) -> str:
    if value is None or pd.isna(value):
        return ""
    text = str(value).strip()
    if text.lower() in {"", "nan", "none", "null"}:
        return ""
    if numeric:
        try:
            num = float(text)
        except Exception:
            return text
        if math.isfinite(num) and num.is_integer():
            return str(int(num))
    return text


def first_present_token(row, keys: Iterable[str], *, numeric: bool = False) -> str:
    for key in keys:
        token = normalize_text_token(row.get(key, ""), numeric=numeric)
        if token:
            return token
    return ""


def build_context_key(chrom: str, *, metagene: str = "", gene: str = "") -> str:
    mg = normalize_text_token(metagene, numeric=True)
    if mg:
        return f"MG:{mg}"
    gene_name = normalize_text_token(gene)
    if gene_name:
        return f"GENE:{gene_name}"
    chrom_name = normalize_text_token(chrom)
    return f"CHR:{chrom_name}" if chrom_name else "CHR:"


def context_key_from_snp_row(row) -> str:
    metagenes = [
        normalize_text_token(token, numeric=True)
        for token in str(row.get("metagene_indices", "")).split(";")
    ]
    metagenes = [token for token in metagenes if token]
    if len(set(metagenes)) == 1 and metagenes:
        return f"MG:{metagenes[0]}"

    genes = [normalize_text_token(token) for token in str(row.get("gene_names", "")).split(";")]
    genes = [token for token in genes if token]
    if len(set(genes)) == 1 and genes:
        return f"GENE:{genes[0]}"

    return build_context_key(str(row.get("chrom", "")))


def context_keys_from_snp_row(row) -> list:
    """ALL context keys a SNP should pair against (a superset of context_key_from_snp_row).

    A SNP overlapping several metagenes (i.e. overlapping genes) is cis to modifications in EACH of
    them, so it must be registered under every MG: track it spans. Collapsing such a SNP to a single
    CHR: key -- as context_key_from_snp_row does -- leaves it unmatchable against the mod side, which
    always keys on one metagene_index -> MG:x, silently dropping it from snp_mod_assoc / haplotype
    associations. Falls back to per-gene GENE: keys, then a single CHR:, mirroring the single-key form
    (a single-metagene SNP returns exactly [MG:x], so behaviour is unchanged for the common case)."""
    metagenes = sorted({normalize_text_token(t, numeric=True)
                        for t in str(row.get("metagene_indices", "")).split(";")
                        if normalize_text_token(t, numeric=True)})
    if metagenes:
        return [f"MG:{m}" for m in metagenes]
    genes = sorted({normalize_text_token(t)
                   for t in str(row.get("gene_names", "")).split(";")
                   if normalize_text_token(t)})
    if genes:
        return [f"GENE:{g}" for g in genes]
    return [build_context_key(str(row.get("chrom", "")))]


def context_key_from_row(
    row,
    *,
    chrom_key: str = "chrom",
    metagene_keys: Iterable[str] = ("metagene_index",),
    gene_keys: Iterable[str] = ("gene_name",),
) -> str:
    return build_context_key(
        str(row.get(chrom_key, "")),
        metagene=first_present_token(row, metagene_keys, numeric=True),
        gene=first_present_token(row, gene_keys),
    )


def normalize_string_series(series: pd.Series, fill_value: str = "") -> pd.Series:
    return series.fillna(fill_value).astype(str).replace({"nan": fill_value, "None": fill_value, "null": fill_value})


# --------------------------------------------------------------------------------------
# Read-key prefiltering for the pairing tests (snp x mod, hap x mod).
#
# All of them inner-join a LARGE per-read table (molecule_snps: 7.5M rows / 1.7GB on Huh7 mock;
# molecule_haplotypes) against a SMALL one (molecule_mod_calls: ~100k rows over ~53k reads), then
# keep only rows sharing (sample, qname). Loading the large table whole costs GiB and a row-wise
# apply over every row -- yet an inner join can never keep a row whose read is absent from the
# small table. So: read the small table first, collect its read keys, then stream the large table
# in chunks and retain only matching rows. Exactly lossless, and peak memory drops to
# O(matching rows) instead of O(whole table).
# --------------------------------------------------------------------------------------

def tsv_header(path: str) -> List[str]:
    with open(path) as fh:
        return fh.readline().rstrip("\n").split("\t")


def shard_tsv_by_chrom(path: str, out_dir: str, chrom_col: str = "chrom") -> Dict[str, str]:
    """Route each raw data line of a TSV to a per-chromosome shard file, preserving exact bytes (so a
    loader parses a shard identically to a chrom-subset of the original). This lets the association
    tests process one chromosome at a time in bounded memory instead of loading the whole many-GB
    table -- peak RAM becomes one chromosome's data, so the pipeline scales to many samples.

    Lossless for these tests: a read maps to a single locus, so all of its mod calls, SNP
    observations and haplotype membership share one chromosome, and every context_key already
    embeds chrom -- no (snp, mod), mod-pair or haplotype group ever spans two chromosomes.

    O(#contigs) open handles + O(1) per-line RAM. Returns {chrom: shard_path} ordered by chrom.
    Handles .gz/.bgz input."""
    os.makedirs(out_dir, exist_ok=True)
    opener = gzip.open if str(path).endswith((".gz", ".bgz")) else open
    writers: Dict[str, object] = {}
    paths: Dict[str, str] = {}
    with opener(path, "rt") as fh:
        header = fh.readline()
        if not header:
            return {}
        cols = header.rstrip("\n").split("\t")
        try:
            ci = cols.index(chrom_col)
        except ValueError:
            raise ValueError(f"shard_tsv_by_chrom: no {chrom_col!r} column in {path}")
        for line in fh:
            parts = line.split("\t", ci + 1)
            if len(parts) <= ci:
                continue
            chrom = parts[ci]
            w = writers.get(chrom)
            if w is None:
                sp = os.path.join(out_dir, "shard_" + chrom.replace("/", "_") + ".tsv")
                w = open(sp, "wt")
                w.write(header)
                writers[chrom] = w
                paths[chrom] = sp
            w.write(line)
    for w in writers.values():
        w.close()
    return dict(sorted(paths.items()))


def read_keys_of(df: pd.DataFrame) -> set:
    """Set of 'sample\\x00qname' keys (a vectorized stand-in for tuple(sample, qname))."""
    if df.empty:
        return set()
    return set(df["sample"].astype(str) + "\x00" + df["qname"].astype(str))


def stream_filter_by_read_keys(
    path: str,
    usecols: List[str],
    read_keys: set,
    *,
    chunksize: int = 500_000,
    row_filter=None,
) -> pd.DataFrame:
    """Chunked read of a large per-read table, keeping only rows whose (sample, qname) appears in
    `read_keys` (and that pass `row_filter`, applied per chunk before the key test).

    `usecols` MUST still include every column that collides with the other table's columns, so the
    downstream merge's ("_snp", "_mod") suffixing is unchanged. Row order is preserved.
    """
    if not read_keys:
        return pd.DataFrame(columns=usecols)
    kept = []
    for chunk in pd.read_csv(path, sep="\t", usecols=usecols, low_memory=False, chunksize=chunksize):
        if row_filter is not None:
            chunk = chunk[row_filter(chunk)]
            if chunk.empty:
                continue
        keys = chunk["sample"].astype(str) + "\x00" + chunk["qname"].astype(str)
        chunk = chunk[keys.isin(read_keys)]
        if not chunk.empty:
            kept.append(chunk)
    if not kept:
        return pd.DataFrame(columns=usecols)
    return pd.concat(kept, ignore_index=True)


def drop_unassigned_reads(mod_df: pd.DataFrame) -> pd.DataFrame:
    """Keep only mod calls on reads assigned to a fragmentform (metagene_index populated).

    An unassigned read (assigned=False) has no fragmentform, so build_molecule_mod_table leaves its
    metagene_index empty and context_key_from_row falls back to GENE:{gene}. The SNP side always
    carries a GTF metagene (MG:{metagene}), so those reads can never pair in snp_mod / snp_tx / hap --
    they are already dropped there by the MG:/GENE: context mismatch, silently. mod_mod, whose pair
    key is context-agnostic, is the ONLY test that counts co-occurrences on these unassigned scrap
    reads. This filter makes the fragmentform scope explicit and CONSISTENT across all four tests
    (a no-op for the SNP-based ones, a correction for mod_mod)."""
    if mod_df.empty or "metagene_index" not in mod_df.columns:
        return mod_df
    mg = mod_df["metagene_index"].astype(str).str.strip()
    return mod_df[mg.ne("") & mg.ne("nan")].copy()


def load_molecule_mods_for_pairing(path: str, extra_cols: Optional[List[str]] = None) -> pd.DataFrame:
    """Load molecule_mod_calls with only the columns the pairing tests need, apply the usable /
    state_detail filters, and add target_state. Same filtering as before, just column-pruned."""
    header = tsv_header(path)
    want = ["sample", "qname", "mod_site_id", "chrom", "start0", "end0",
            "target_mod_code", "state_detail", "gene_name", "metagene_index"]
    if "usable" in header:
        want.append("usable")
    else:
        want.extend([c for c in ("fail", "within_alignment") if c in header])
    for c in (extra_cols or []):
        if c in header and c not in want:
            want.append(c)
    usecols = [c for c in want if c in header]

    mod_df = pd.read_csv(path, sep="\t", usecols=usecols, low_memory=False)
    if mod_df.empty:
        return mod_df
    if "usable" in mod_df.columns:
        mod_df = mod_df[mod_df["usable"].fillna(False)].copy()
    else:
        mod_df = mod_df[(~mod_df["fail"].fillna(True)) & mod_df["within_alignment"].fillna(False)].copy()
    mod_df = mod_df[mod_df["state_detail"].isin(["modified", "canonical", "other_mod"])].copy()
    mod_df = drop_unassigned_reads(mod_df)
    if not mod_df.empty:
        mod_df["target_state"] = mod_df["state_detail"].eq("modified").astype(int)
    return mod_df


def run_process_jobs(fn, task_args: List[tuple], jobs: int, *, verbose: bool = False, label: str = "parallel jobs"):
    if not task_args:
        return []
    jobs = max(1, min(int(jobs), len(task_args)))
    if jobs <= 1 or len(task_args) == 1:
        return [fn(*args) for args in task_args]

    try:
        results = []
        with ProcessPoolExecutor(max_workers=jobs) as executor:
            future_map = {executor.submit(fn, *args): args for args in task_args}
            for future in as_completed(future_map):
                results.append(future.result())
        return results
    except Exception as exc:
        if verbose:
            print(f"[warn] Falling back to serial {label}: {exc}", file=sys.stderr, flush=True)
        return [fn(*args) for args in task_args]
