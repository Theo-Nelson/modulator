#!/usr/bin/env python3

from concurrent.futures import ProcessPoolExecutor, as_completed
import gzip
import math
import os
import sys
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np
import pandas as pd
from scipy.stats import chi2, chi2_contingency, fisher_exact, random_table, rankdata as _rankdata


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


def stratified_max_distribution_shift(strata) -> float:
    """Sample-stratified analogue of max_abs_distribution_shift for r x c tables: for each row pair
    (i, j) and each column, the coverage-weighted mean over strata of (row_i_frac - row_j_frac); the
    reported effect is the max |.| over columns and row pairs. Weights w_k = R_i R_j / N_k (matching
    mh_stratified_effect). This is the effect size consistent with cmh_stratified_test -- within-sample
    fractions, not pooled -- so a stratified CMH p is never paired with a confounded pooled effect."""
    strata = [np.asarray(T, dtype=float) for T in strata if np.asarray(T).size]
    if not strata:
        return 0.0
    r, c = strata[0].shape
    if r < 2 or c < 1:
        return 0.0
    best = 0.0
    for i in range(r):
        for j in range(i + 1, r):
            num = np.zeros(c)
            den = 0.0
            for T in strata:
                if T.shape != (r, c):
                    continue
                Ri, Rj, N = T[i].sum(), T[j].sum(), T.sum()
                if Ri <= 0 or Rj <= 0 or N <= 0:
                    continue
                w = Ri * Rj / N
                num += w * (T[i] / Ri - T[j] / Rj)
                den += w
            if den > 0:
                best = max(best, float(np.max(np.abs(num / den))))
    return round(best, 6)


def binary_rate_delta(table_2x2: np.ndarray) -> float:
    table = np.asarray(table_2x2, dtype=float)
    if table.shape != (2, 2):
        return 0.0
    rate0 = table[0, 0] / table[0].sum() if table[0].sum() > 0 else 0.0
    rate1 = table[1, 0] / table[1].sum() if table[1].sum() > 0 else 0.0
    return round(float(abs(rate0 - rate1)), 6)


def montecarlo_exact_test(table, seed: int = 12345, n_resamples: int = 9999):
    """Monte-Carlo EXACT test of independence for an r x c contingency table -- Patefield resampling of
    tables with BOTH margins fixed (R's chisq.test(simulate.p.value=TRUE)). Deterministic (fixed seed).

    This is the exact replacement for the asymptotic chi2 in the single-informative-stratum primary,
    which is ~2-5x anti-conservative on sparse / unequal-coverage / low-rate tables -- the m6A regime,
    and EVERY site on single-sample data (where every row takes that branch). Fixing BOTH margins makes
    it correct for 2 x c as well as r x 2 (a row-margin-only sampler returns an incoherent number on a
    2 x c table). The MC p has a resolution floor of 1/(n_resamples+1); when the observed statistic
    exceeds every resample -- a strong, well-conditioned table beyond MC resolution -- the asymptotic
    chi2 p is accurate there and is reported instead, so the floor never (a) makes a clean study with a
    few strong effects report NONE under BH nor (b) ties all strong hits at 1e-4 and destroys ranking.
    Robust to zero margins: an all-zero ROW (a transcript with no reads) or COLUMN (an outcome with
    none) carries no information and would put a 0 in the expected table -> 0/0 = NaN p, silently
    dropping a genuine hit from the BH family. run_contingency_test screens zero margins, but DIRECT
    callers do not (test_stoichiometry_diffs' single-informative-stratum path), and _stratum_informative
    only requires >=2 nonzero rows -- so at r>2 a zero row survives into this table. Drop all-zero rows
    and columns here; if <2x2 remains the table is genuinely untestable (NaN, leaves the BH family)."""
    tab = np.asarray(table, dtype=float)
    tab = tab[tab.sum(axis=1) > 0]        # drop all-zero rows
    tab = tab[:, tab.sum(axis=0) > 0]     # drop columns left all-zero (recomputed after the row drop)
    if tab.shape[0] < 2 or tab.shape[1] < 2:
        return "untestable", "none", float("nan"), float("nan")
    r, c = tab.shape
    row = tab.sum(axis=1); col = tab.sum(axis=0); N = tab.sum()
    exp = np.outer(row, col) / N
    obs_stat = float(np.sum((tab - exp) ** 2 / exp))
    dof = (r - 1) * (c - 1)
    asymptotic_p = float(chi2.sf(obs_stat, dof))
    rng = np.random.default_rng(seed)
    samples = random_table(row.astype(int), col.astype(int)).rvs(size=n_resamples, random_state=rng)
    stats = np.sum((samples - exp) ** 2 / exp, axis=(1, 2))
    hits = int(np.sum(stats >= obs_stat - 1e-9))
    p = asymptotic_p if hits == 0 else (hits + 1) / (n_resamples + 1)
    return f"montecarlo_exact_{r}x{c}", "chi2", obs_stat, float(p)


def run_contingency_test(
    table: np.ndarray,
    test: str = "auto",
    pseudocount: float = 0.5,
) -> Tuple[str, str, float, float]:
    tab = np.asarray(table, dtype=float)
    if tab.size == 0 or tab.shape[0] < 2 or tab.shape[1] < 2:
        return "none", "none", 0.0, 1.0

    # Drop all-zero ROWS/COLUMNS (a group or outcome with no reads carries no information) BEFORE
    # deciding testability. Strata are built over a FIXED row set -- a group absent from one sample is
    # a zero row by construction, routine at 3v3 -- and _stratum_informative keeps such tables (it
    # requires >=2 nonzero rows, not every row nonzero). The OLD zero-margin screen instead marked the
    # WHOLE table untestable (NaN), silently dropping a genuinely-significant reduced table from the BH
    # family -- and disagreeing with montecarlo_exact_test, which reduces internally for its direct
    # callers. Reduce here so the dispatch (2x2->Fisher, else->MC exact) sees the informative survivors.
    tab = tab[tab.sum(axis=1) > 0][:, tab.sum(axis=0) > 0]
    # After the reduction there is no all-zero row/column left; the only remaining untestable case is a
    # table that collapsed below 2x2 (all variation was on a single removed axis -- e.g. 100%/100%).
    if tab.shape[0] < 2 or tab.shape[1] < 2:
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

    # Non-2x2 tables use the Monte-Carlo EXACT test, NOT the asymptotic chi2 (M2: chi2 is ~2-5x
    # anti-conservative on the sparse/low-rate tables this branch actually sees, and on single-sample
    # data EVERY site takes it). `test="chi2"` still forces the asymptotic test as an explicit override.
    if test == "chi2":
        return do_chi2(tab)
    if tab.shape == (2, 2):
        return do_fisher_2x2(tab)
    return montecarlo_exact_test(tab)


def cmh_stratified_test(strata):
    """Generalized Cochran-Mantel-Haenszel GENERAL-ASSOCIATION test over a list of r x 2 tables, one
    per SAMPLE (rows in a FIXED order; cols = [positive, negative]). Reduces to the standard 2x2 CMH.

    This is the sample-stratified replacement for a sample-POOLED Fisher/chi2 on read-level molecule
    data: pooling reads across technical/biological replicates lets per-sample rate + composition
    imbalance manufacture a Simpson's-paradox association where the within-sample effect is zero. A
    stratum with no positive OR no negative outcome, or <2 covered rows, carries no information and is
    dropped (the intended power cost -- a sample loading only one row can't separate row from sample).

    Returns (test_name, stat_name, stat, p_value, n_informative_strata); NaN p when no stratum is
    informative (so it is excluded from the BH family, like run_contingency_test)."""
    strata = [np.asarray(T, dtype=float) for T in strata if np.asarray(T).size]
    if not strata:
        return "cmh_untestable", "cmh_chi2", float("nan"), float("nan"), 0
    r, c = strata[0].shape
    if r < 2 or c < 2:
        return "cmh_untestable", "cmh_chi2", float("nan"), float("nan"), 0
    # Landis-Koch general-association statistic: obs-expected over the first (r-1)x(c-1) cells, with
    # covariance kron(A_r, A_c)/(N-1). For c==2 this reduces exactly to the standard 2x2/rx2 CMH.
    d = (r - 1) * (c - 1)
    A = np.zeros(d)
    V = np.zeros((d, d))
    used = 0
    for T in strata:
        R = T.sum(axis=1); C = T.sum(axis=0); N = T.sum()
        if N < 2 or (R > 0).sum() < 2 or (C > 0).sum() < 2:
            continue
        Rr, Cc = R[:r - 1], C[:c - 1]
        obs = T[:r - 1, :c - 1].reshape(-1)
        m = np.outer(Rr, Cc).reshape(-1) / N
        Ar = np.diag(Rr) - np.outer(Rr, Rr) / N
        Ac = np.diag(Cc) - np.outer(Cc, Cc) / N
        A += obs - m
        V += np.kron(Ar, Ac) / (N - 1)
        used += 1
    if used == 0:
        return "cmh_untestable", "cmh_chi2", float("nan"), float("nan"), 0
    try:
        Q = float(A @ np.linalg.solve(V, A))
    except np.linalg.LinAlgError:
        Q = float(A @ np.linalg.pinv(V) @ A)
    if not np.isfinite(Q) or Q < 0:
        return "cmh_untestable", "cmh_chi2", float("nan"), float("nan"), used
    name = "cmh_2x2" if (r == 2 and c == 2) else f"cmh_general_{r}x{c}"
    return name, "cmh_chi2", Q, float(chi2.sf(Q, d)), used


def _stratum_informative(T) -> bool:
    """A stratum carries within-stratum information iff N>=2 AND >=2 rows and >=2 columns have reads --
    exactly the strata cmh_stratified_test keeps. NOTE: this requires >=2 NONZERO rows/columns, NOT
    that EVERY row/column is nonzero, so at r>2 an informative stratum can still contain a zero row (a
    transcript with no reads in that sample); montecarlo_exact_test drops such rows/columns itself."""
    T = np.asarray(T, dtype=float)
    if T.size == 0:
        return False
    R = T.sum(axis=1); C = T.sum(axis=0); N = T.sum()
    return N >= 2 and int((R > 0).sum()) >= 2 and int((C > 0).sum()) >= 2


def informative_strata(strata):
    """The subset of `strata` that carry within-stratum information (see _stratum_informative), with
    LEVELS (rows/columns) that are zero in EVERY informative stratum dropped. Filtering once here lets
    the caller reuse the SAME informative set for the primary test, the odds ratio, the effect, and the
    heterogeneity test -- so none of them silently reverts to a table that pools reads back in.

    Why the level-reduction (fixing the common point once, so CMH + heterogeneity + the exact test all
    get it): a level with no reads in ANY informative stratum is a phantom dimension. It costs
    cmh_stratified_test a spurious degree of freedom -> p inflated (12.6x per empty level, real hits
    lost -- CONSERVATIVE), and it makes stratum_heterogeneity rank-deficient relative to its Agresti-
    smoothed covariance so W = inv(Cov) blows up -> false heterogeneity flags (ANTI-conservative:
    null type-I 0.56-0.81 when the last column is empty). run_contingency_test already reduced for the
    single-stratum exact path; doing it HERE gives every consumer the reduced strata, so a site can no
    longer get the honest answer at 1 stratum and an inflated one at >=2. CMH needs a fixed row/col set,
    so a level is dropped only when it is zero across ALL strata (kept if nonzero in any)."""
    inf = [np.asarray(T, dtype=float) for T in strata if _stratum_informative(T)]
    if not inf or len({T.shape for T in inf}) != 1:   # ragged strata (shouldn't happen) -> don't reduce
        return inf
    stack = np.stack(inf)                              # (K, r, c)
    row_keep = stack.sum(axis=(0, 2)) > 0             # rows nonzero in >=1 stratum
    col_keep = stack.sum(axis=(0, 1)) > 0             # cols nonzero in >=1 stratum
    if row_keep.all() and col_keep.all():
        return inf
    return [T[row_keep][:, col_keep] for T in inf]


def stratified_primary(inf_strata, exact_test, min_strata=2):
    """Primary test computed from the INFORMATIVE strata only (pre-filter with informative_strata):

    - >= min_strata informative strata -> sample-stratified CMH.
    - exactly 1 informative stratum -> the EXACT test on THAT stratum alone. NOT a table that pools the
      dropped (non-informative) samples' reads back in: below 2 informative strata the asymptotic CMH is
      anti-conservative ([[6,1],[1,6]] -> 0.010 vs Fisher 0.029), but the fully-pooled table is the very
      Simpson statistic the stratification removed -- one exactly-independent informative sample plus two
      constant samples pools to Fisher p=5e-81 / OR=78, where the honest answer is p=1.0 / OR=1.0.
    - 0 informative strata -> NaN, so the row leaves the BH family (like any untestable row).

    `exact_test(table)` returns (name, stat_name, stat, p). The odds ratio / effect the CALLER reports
    should be computed from inf_strata too (Mantel-Haenszel over the informative strata), so MAJOR-B's
    'primary is sample-adjusted' guarantee holds on the fallback path as well. Returns
    (test_name, stat_name, stat, p_value, n_informative, mode) with mode in {cmh, exact_single, none}."""
    n = len(inf_strata)
    if n >= min_strata:
        name, sn, stat, p, _n = cmh_stratified_test(inf_strata)
        return (name, sn, stat, p, n, "cmh")
    if n == 1:
        name, sn, stat, p = exact_test(np.asarray(inf_strata[0], dtype=float))
        return (name, sn, stat, p, n, "exact_single")
    return ("untestable", "none", float("nan"), float("nan"), 0, "none")


def stratum_heterogeneity(inf_strata, min_row=10, min_col=5, smooth=1.0):
    """Cochran's Q test of effect HOMOGENEITY across the INFORMATIVE strata, on the row-conditional
    PROPORTION (RISK-DIFFERENCE) scale -- the SAME scale as the reported effect (mh_stratified_effect /
    stratified_max_distribution_shift). Testing homogeneity on the odds-ratio scale (the earlier
    generalized-CMH decomposition) FALSE-FLAGS the common real case where the risk difference is constant
    across samples but the baseline rate varies, so the OR moves: a visibly constant effect column would
    then sit beside "not homogeneous". Here per stratum the row-0-referenced proportion differences
    theta_k = {pi_ij - pi_0j : i=1..r-1, j=0..c-2}, pi_ij = T_ij / rowsum_i, are combined by
    inverse-variance into Q = sum (theta_k - theta_bar)' W_k (theta_k - theta_bar) ~ chi2_{(K-1)(r-1)(c-1)}.
    For 2x2 this is exactly Cochran's Q on the per-stratum rate differences, so a constant risk difference
    gives Q=0 regardless of baseline.

    LOW-COUNT ROBUSTNESS: the plug-in Wald variance p(1-p)/n collapses to ~0 when an observed proportion
    hits 0 or 1 (routine at low depth), so a bare ridge let Q -> huge and the flag fired at p=0.0 on every
    row (measured null 0.49-0.65 at 2-4 reads/row). Two guards: (1) the VARIANCE uses Agresti-smoothed
    proportions (x+smooth)/(n+2*smooth) so no cell has 0 variance (the point estimate theta stays RAW, so
    the risk-difference scale is unchanged); (2) a stratum enters only if EVERY row has >= min_row reads --
    below that the normal approximation is unreliable. With both, the null is ~nominal (0.037 on HG002-like
    K=2 depths, 0.06-0.08 at 16-32 reads/row), the p-values are finite (so add_heterogeneity_flag's BH
    adjustment gives real protection: post-BH false-flag ~0 on homogeneous tables), and >=98% of rows are
    still tested. Returns (hetero_stat, hetero_p, hetero_df, n_usable); NaN when <2 usable strata.

    NEAR-EMPTY COLUMN GUARD (min_col): the column-wise mirror of the per-stratum min_row guard.
    informative_strata already drops any response level empty in EVERY stratum, but a column with just
    a few reads POOLED across strata survives and still gives a near-degenerate row-conditional
    proportion whose Agresti-smoothed variance under-states the true noise -> Q inflation -> false
    heterogeneity (measured null ~0.26 at <5 pooled reads, well above nominal). So for c>=3 any column
    whose reads pooled across strata total < min_col is dropped before the test; the response is never
    reduced below 2 columns (untestable -> NaN). 2x2 tables are left untouched (nothing to drop)."""
    inf = [np.asarray(T, dtype=float) for T in inf_strata]
    if len(inf) < 2:
        return float("nan"), float("nan"), 0, len(inf)
    r, c = inf[0].shape
    inf = [T for T in inf if T.shape == (r, c)]                # common-shape strata only (mirrors in-loop guard)
    if len(inf) < 2:
        return float("nan"), float("nan"), 0, len(inf)
    if c > 2:                                                  # pooled near-empty column guard (c>=3 only)
        # Pool the column totals over the strata the test WILL ACTUALLY USE -- those passing the
        # min_row guard applied in the loop below -- NOT all informative strata. A column whose reads
        # live entirely in a min_row-SKIPPED stratum would otherwise clear a pooled-over-everything
        # guard yet be empty in every used stratum, and the near-empty degeneracy returns at full
        # magnitude (same shape as the globally- and near-empty column bugs, one level deeper:
        # informative-but-skipped). min_row routinely discards strata informative_strata kept.
        used = [T for T in inf if not (T.sum(axis=1) < min_row).any()]
        if len(used) >= 2:
            col_pooled = np.sum([T.sum(axis=0) for T in used], axis=0)
            keep = col_pooled >= min_col
            if int(keep.sum()) < 2:
                return float("nan"), float("nan"), 0, len(used)  # <2 testable response levels left
            if not keep.all():
                inf = [T[:, keep] for T in inf]
                r, c = inf[0].shape
    d = (r - 1) * (c - 1)
    if d < 1:
        return float("nan"), float("nan"), 0, len(inf)
    thetas = []
    Ws = []
    for T in inf:
        if T.shape != (r, c):
            continue
        R = T.sum(axis=1)
        if (R < min_row).any():                              # per-stratum minimum-count guard
            continue
        pi = T / R[:, None]                                   # RAW row-conditional proportions (effect scale)
        pv = (T + smooth) / (R[:, None] + 2.0 * smooth)       # Agresti-smoothed props for the VARIANCE only
        theta = (pi[1:, :c - 1] - pi[0:1, :c - 1]).reshape(-1)  # (r-1)(c-1) row-0-referenced diffs

        def _rowcov(i):                                       # multinomial cov of row i's first c-1 (smoothed)
            p = pv[i, :c - 1]
            return (np.diag(p) - np.outer(p, p)) / R[i]
        V0 = _rowcov(0)
        m = r - 1
        Cov = np.zeros((m * (c - 1), m * (c - 1)))
        for a in range(m):
            Va = _rowcov(a + 1)
            for b in range(m):
                blk = V0 + (Va if a == b else 0.0)            # shared row-0 term couples the blocks
                Cov[a * (c - 1):(a + 1) * (c - 1), b * (c - 1):(b + 1) * (c - 1)] = blk
        try:
            W = np.linalg.inv(Cov)
        except np.linalg.LinAlgError:
            W = np.linalg.pinv(Cov)
        thetas.append(theta)
        Ws.append(W)
    if len(thetas) < 2:
        return float("nan"), float("nan"), 0, len(thetas)
    SW = sum(Ws)
    try:
        SWinv = np.linalg.inv(SW)
    except np.linalg.LinAlgError:
        SWinv = np.linalg.pinv(SW)
    theta_bar = SWinv @ sum(W @ th for W, th in zip(Ws, thetas))
    Q = 0.0
    for W, th in zip(Ws, thetas):
        dif = th - theta_bar
        Q += float(dif @ W @ dif)
    Q = max(0.0, Q)
    df = (len(thetas) - 1) * d
    return Q, float(chi2.sf(Q, df)), df, len(thetas)


def add_heterogeneity_flag(out_df, alpha=0.05):
    """BH-adjust the per-row strata_heterogeneity_p across the whole table and set strata_heterogeneous =
    (adjusted < alpha). Every other per-row p in these outputs is BH-adjusted; a RAW p<0.05 heterogeneity
    flag would fire on ~5% of perfectly homogeneous rows (thousands of false flags on a 260k-row table).
    NaN heterogeneity p (rows whose strata never passed the heterogeneity test's own count guard, even if
    n_strata_informative reports otherwise) is excluded from the family and never flagged.

    strata_heterogeneous is a NULLABLE boolean: True = tested + heterogeneous, False = tested +
    homogeneous, <NA> = heterogeneity NOT assessed for this row. Collapsing the untested rows to plain
    False would read as "checked and clean" when nothing was checked -- a table where no row was testable
    would otherwise look uniformly homogeneous. Adds strata_heterogeneity_p_adj and (re)writes
    strata_heterogeneous. No-op if the column is absent."""
    if out_df is None or "strata_heterogeneity_p" not in getattr(out_df, "columns", []):
        return out_df
    padj = benjamini_hochberg(pd.to_numeric(out_df["strata_heterogeneity_p"], errors="coerce").values)
    out_df["strata_heterogeneity_p_adj"] = padj
    padj_num = pd.to_numeric(out_df["strata_heterogeneity_p_adj"], errors="coerce")
    het = (padj_num < alpha).astype("boolean")   # NaN comparison -> False; fix below
    out_df["strata_heterogeneous"] = het.where(padj_num.notna(), other=pd.NA)  # untested -> <NA>, not False
    return out_df


def mh_common_odds_ratio(strata):
    """Mantel-Haenszel common odds ratio across 2x2 strata [[a,b],[c,d]] (one per sample):
    OR_MH = sum(a*d/N_k) / sum(b*c/N_k). This is the stratum-ADJUSTED odds ratio consistent with
    cmh_stratified_test -- unlike a pooled OR it is not inflated by between-sample rate/composition
    imbalance, so it is the odds ratio (and hence the CONCORDANT/MUTUALLY_EXCLUSIVE direction) a reader
    should see next to a stratified CMH p. Returns NaN when no stratum is informative or the denominator
    is 0 (no discordant mass)."""
    num = den = 0.0
    for T in strata:
        T = np.asarray(T, dtype=float)
        if T.shape != (2, 2):
            continue
        N = T.sum()
        if N <= 0:
            continue
        num += T[0, 0] * T[1, 1] / N
        den += T[0, 1] * T[1, 0] / N
    if den <= 0:
        return float("nan")
    return num / den


def mh_stratified_effect(strata):
    """Mantel-Haenszel coverage-weighted rate difference (positive-fraction), max over row pairs -- the
    effect size consistent with cmh_stratified_test. Weights w_k = R_i R_j / N_k over strata covering
    both rows."""
    strata = [np.asarray(T, dtype=float) for T in strata if np.asarray(T).size]
    if not strata:
        return 0.0
    r = strata[0].shape[0]
    best = 0.0
    for i in range(r):
        for j in range(i + 1, r):
            num = den = 0.0
            for T in strata:
                Ri, Rj, N = T[i].sum(), T[j].sum(), T.sum()
                if N <= 0 or Ri <= 0 or Rj <= 0:
                    continue
                w = Ri * Rj / N
                num += w * (T[i, 0] / Ri - T[j, 0] / Rj); den += w
            if den > 0:
                best = max(best, abs(num / den))
    return best


def van_elteren_kw(strata):
    """Generalized van Elteren stratified rank test (stratified Kruskal-Wallis) for k>=2 groups.

    `strata` is a list of strata; each stratum is a list of 1D arrays, ONE PER GROUP in a FIXED order
    across every stratum (a group absent from a stratum is an empty array). Within each stratum the
    observations are midrank-scored and each group's rank-sum deviation from its null expectation is
    accumulated, weighted across strata by the design-free w_k = 1/(N_k + 1) -- the same weight the
    2-group van_elteren_stratified uses -- into a Landis-Koch quadratic form ~ chi2_{g-1}.

    This is the k-group, sample-stratified replacement for a sample-POOLED Kruskal/Mann-Whitney on
    read-level continuous data (e.g. tail length by fragmentform): pooling reads across replicates lets
    per-sample level + composition imbalance manufacture a Simpson's-paradox difference where the
    within-sample effect is ~0. A single stratum reduces to the tie-corrected Kruskal-Wallis H; two
    groups reduce to the 2-group van Elteren (chi2_1 = z^2). A stratum with <2 non-empty groups, <2
    observations, or no rank variation (all tied) carries no information and is dropped.

    Returns (stat, p_value, n_informative_strata, df); NaN when nothing is informative."""
    strata = [st for st in strata if st]
    if not strata:
        return float("nan"), float("nan"), 0, 0
    g = len(strata[0])
    if g < 2 or any(len(st) != g for st in strata):
        return float("nan"), float("nan"), 0, 0
    A = np.zeros(g)
    V = np.zeros((g, g))
    present = np.zeros(g, dtype=bool)
    used = 0
    for st in strata:
        arrs = [np.asarray(a, dtype=float).ravel() for a in st]
        ns = np.array([a.size for a in arrs], dtype=float)
        N = int(ns.sum())
        if N < 2 or int((ns > 0).sum()) < 2:
            continue
        allv = np.concatenate([a for a in arrs if a.size])
        ranks = _rankdata(allv)
        rbar = (N + 1) / 2.0
        sigma2 = float(((ranks - rbar) ** 2).sum()) / N   # population variance of midranks
        if sigma2 <= 0:                                    # every value tied -> no rank information
            continue
        w = 1.0 / (N + 1)
        Rk = np.zeros(g)
        off = 0
        for i, a in enumerate(arrs):
            if a.size:
                Rk[i] = float(ranks[off:off + a.size].sum())
                off += a.size
        Ek = ns * rbar
        cov = sigma2 / (N - 1) * (N * np.diag(ns) - np.outer(ns, ns))
        A += w * (Rk - Ek)
        V += w * w * cov
        present |= ns > 0
        used += 1
    if used == 0:
        return float("nan"), float("nan"), 0, 0
    idx = np.where(present)[0]
    if idx.size < 2:
        return float("nan"), float("nan"), used, 0
    keep = idx[:-1]                                        # drop one group (rows are linearly dependent)
    Ar = A[keep]
    Vr = V[np.ix_(keep, keep)]
    try:
        Q = float(Ar @ np.linalg.solve(Vr, Ar))
    except np.linalg.LinAlgError:
        Q = float(Ar @ np.linalg.pinv(Vr) @ Ar)
    df = int(keep.size)
    if not np.isfinite(Q) or Q < 0:
        return float("nan"), float("nan"), used, df
    return Q, float(chi2.sf(Q, df)), used, df


def van_elteren_heterogeneity(strata):
    """Rank-scale test of HOMOGENEITY of the stratified tail effect across strata -- the tail-length
    analogue of stratum_heterogeneity, so test_taillength_diffs gets a heterogeneity diagnostic too.

    It is a Cochran's Q on the NONPARAMETRIC RELATIVE EFFECT p_ik = mean-midrank_ik / (N_k + 1) of each
    fragmentform (dropping one), NOT on the raw rank-sum scores: the raw score scales with stratum size,
    so a common non-zero effect at UNEQUAL read depths (the normal case -- it is why van Elteren weights
    strata) makes the scores differ and a score-based Q false-flags (measured: 0.03 at equal size ->
    0.70+ at 5x). The relative effect is on [0,1] regardless of N -- the rank analog of the count test's
    risk-difference scale -- so per-stratum p_ik are directly comparable and inverse-variance Cochran's Q
    stays calibrated (verified flat ~0.03 to 5x depth ratio and k=3) while keeping power (reversed
    ordering -> 1.0, a +6-vs-+14 nt magnitude difference -> 0.9).

    Every p_ik must be on ONE reference frame: the ranking is done SEPARATELY per stratum but restricted
    to the fragmentforms present in EVERY informative stratum (the "common" forms). Ranking over all forms
    present in a stratum first and then comparing the common ones would let a common form's relative
    effect move when its NON-common competitors change -- forms A,B identical across two samples but a
    third form C present in only one sample (or with a different tail) shifted A,B's ranks and false-
    flagged. Re-ranking each stratum on the common forms alone removes that dependence. `strata` has
    van_elteren_kw's shape. Returns (hetero_stat, hetero_p, hetero_df, n_informative); NaN when <2
    informative strata or <2 common forms."""
    strata = [st for st in strata if st]
    if not strata:
        return float("nan"), float("nan"), 0, 0
    g = len(strata[0])
    if g < 2 or any(len(st) != g for st in strata):
        return float("nan"), float("nan"), 0, 0
    # informative strata (>=2 groups with reads, N>=2) and their per-group presence
    cand = []
    for st in strata:
        arrs = [np.asarray(a, dtype=float).ravel() for a in st]
        ns = np.array([a.size for a in arrs], dtype=float)
        if int(ns.sum()) < 2 or int((ns > 0).sum()) < 2:
            continue
        cand.append(arrs)
    if len(cand) < 2:
        return float("nan"), float("nan"), 0, len(cand)
    common = np.all([[a.size > 0 for a in arrs] for arrs in cand], axis=0)  # forms present in EVERY stratum
    cidx = np.where(common)[0]
    if cidx.size < 2:
        return float("nan"), float("nan"), 0, len(cand)
    keep = np.arange(cidx.size - 1)                          # local index over common forms; drop one
    es = []
    Ws = []
    used = 0
    for arrs in cand:
        cg = [arrs[i] for i in cidx]                         # rank ONLY the common forms -> one reference
        ns = np.array([a.size for a in cg], dtype=float)
        N = int(ns.sum())
        if N < 2 or int((ns > 0).sum()) < 2:
            continue
        ranks = _rankdata(np.concatenate([a for a in cg if a.size]))
        rbar = (N + 1) / 2.0
        sigma2 = float(((ranks - rbar) ** 2).sum()) / N
        if sigma2 <= 0:
            continue
        Rc = np.zeros(cidx.size)
        off = 0
        for j, a in enumerate(cg):
            if a.size:
                Rc[j] = float(ranks[off:off + a.size].sum())
                off += a.size
        cov = sigma2 / (N - 1) * (N * np.diag(ns) - np.outer(ns, ns))
        with np.errstate(divide="ignore", invalid="ignore"):
            e = ((Rc / ns) / (N + 1))[keep]                  # relative effect = (R_j / n_j) / (N+1)
            C = (cov / (np.outer(ns, ns) * (N + 1) ** 2))[np.ix_(keep, keep)]
        Ws.append(np.linalg.pinv(C))
        es.append(e)
        used += 1
    if used < 2:
        return float("nan"), float("nan"), 0, used
    SW = sum(Ws)
    SWinv = np.linalg.pinv(SW)
    e_bar = SWinv @ sum(W @ e for W, e in zip(Ws, es))
    Q = 0.0
    for W, e in zip(Ws, es):
        dif = e - e_bar
        Q += float(dif @ W @ dif)
    Q = max(0.0, Q)
    df = (used - 1) * int(keep.size)
    if df < 1:
        return float("nan"), float("nan"), 0, used
    return Q, float(chi2.sf(Q, df)), df, used


def weighted_within_stratum_median_range(strata):
    """Effect size consistent with van_elteren_kw: within each stratum the range (max-min) of the
    per-group medians, averaged over informative strata weighted by stratum size. `strata` has the same
    shape as van_elteren_kw's argument."""
    num = den = 0.0
    for st in strata:
        arrs = [np.asarray(a, dtype=float).ravel() for a in st]
        meds = [float(np.median(a)) for a in arrs if a.size]
        if len(meds) < 2:
            continue
        w = float(sum(a.size for a in arrs))
        num += w * (max(meds) - min(meds))
        den += w
    return num / den if den > 0 else 0.0


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
