#!/usr/bin/env python3

from concurrent.futures import ProcessPoolExecutor, as_completed
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
    order = np.argsort(p)
    ranks = np.empty(p.size, dtype=int)
    ranks[order] = np.arange(1, p.size + 1)
    adj = p * p.size / ranks
    adj_sorted = np.minimum.accumulate(adj[order][::-1])[::-1]
    out = np.empty_like(adj)
    out[order] = adj_sorted
    return np.minimum(out, 1.0)


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

    def do_fisher_2x2(tt):
        odds, p = fisher_exact(tt.astype(int))
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


def cmh_test_2x2xk(strata: List[np.ndarray]) -> Tuple[float, float, float]:
    valid = []
    for tab in strata:
        tt = np.asarray(tab, dtype=float)
        if tt.shape != (2, 2):
            continue
        if tt.sum() <= 1:
            continue
        valid.append(tt)
    if not valid:
        return 0.0, 1.0, 0.0

    num = 0.0
    den = 0.0
    or_num = 0.0
    or_den = 0.0
    for tt in valid:
        a, b = tt[0, 0], tt[0, 1]
        c, d = tt[1, 0], tt[1, 1]
        n = a + b + c + d
        row1 = a + b
        row2 = c + d
        col1 = a + c
        col2 = b + d
        expected_a = (row1 * col1) / n if n > 0 else 0.0
        if n > 1:
            var_a = (row1 * row2 * col1 * col2) / (n * n * (n - 1))
        else:
            var_a = 0.0
        num += (a - expected_a)
        den += var_a
        or_num += (a * d) / n if n > 0 else 0.0
        or_den += (b * c) / n if n > 0 else 0.0

    if den <= 0:
        return 0.0, 1.0, 0.0
    stat = (num * num) / den
    p_value = float(chi2.sf(stat, 1))
    common_or = float(or_num / or_den) if or_den > 0 else float("inf")
    return float(stat), p_value, common_or


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
