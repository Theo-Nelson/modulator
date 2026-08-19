#!/usr/bin/env python3
"""Replicate-aware between-condition statistics.

Every between-sample comparison in this pipeline is one of two shapes:

  * COUNTS -- (successes, total) per sample: modified/covered reads at a site, reads on an isoform /
    APA site / junction out of its gene's reads. Handled by ``beta_binomial_diff`` below.
  * CONTINUOUS -- one value per sample (e.g. median poly(A) tail). Handled by ``continuous_diff``.

THE POINT OF THIS MODULE IS TO NOT COMMIT PSEUDOREPLICATION. With millions of reads but n=3 per
group, the biological unit is the REPLICATE, not the read. Pooling reads across replicates and
running Fisher gives p=1e-300 for trivial differences. Two noise sources must stay separated:

  * binomial sampling noise  -- k ~ Binomial(n, p_i); scales with coverage, so 2/5 is nearly
    uninformative while 400/1000 is tight. (A test on collapsed fractions is blind to this: it
    cannot tell 0.4 from 2/5 apart from 0.4 from 400/1000, and 2/5-vs-1/5 across three replicates
    has ZERO fraction-variance -> t -> infinity -> a spurious "hit" built on 5 reads.)
  * beta biological noise    -- p_i ~ Beta(mu_group, dispersion); replicate-to-replicate spread.

Beta-binomial log-likelihood, mean/precision parameterised (alpha = mu*theta, beta = (1-mu)*theta):

    ll(k; n, mu, theta) = betaln(k + mu*theta, n - k + (1-mu)*theta) - betaln(mu*theta, (1-mu)*theta)

Dispersion phi = 1/(theta+1). At n=3 per group a PER-SITE theta is badly estimated (~4 residual df),
so -- exactly as DSS/edgeR do -- theta is SHRUNK toward a trend fitted across all sites. That
shrinkage is where the power at n=3 comes from. Test = likelihood-ratio (2*(ll_full - ll_null)),
chi2 with 1 df, on the shrunk theta.

scipy only (no statsmodels).
"""
import numpy as np
from scipy.optimize import minimize, minimize_scalar
from scipy.special import betaln, expit, logit
from scipy.stats import f as fdist
from scipy.stats import ttest_ind

_EPS = 1e-9
_MIN_MU, _MAX_MU = 1e-6, 1 - 1e-6
_LOG_THETA_LO, _LOG_THETA_HI = np.log(1e-2), np.log(1e6)
# Reference distribution for the LRT: F(1, REF_DF).
#
# chi2(1) is NOT valid here. At n=3v3 the LRT is inflated purely because theta is ESTIMATED (with
# the TRUE theta the LRT is already chi2(1)-calibrated at n=3: mean 1.09, KS p=0.49). The reference
# therefore has to absorb the uncertainty in theta -- exactly edgeR's quasi-likelihood F idea, where
# df reflects how well the dispersion is known.
#
# REF_DF=10 was calibrated against the OBSERVED dispersion of this data (see DISPERSION note below),
# by simulation from the empirical (phi, coverage, mu) joint of the mock replicates:
#     NULL p<0.05 = 0.050 (nominal)   |   delta=0.10 -> 69% sensitivity at ~5% FDR
#     (chi2(1) -> 0.081 null, 24% FDR;  F(1,4) -> 0.021 null but only 14% sensitivity)
#
# DISPERSION note: the 3 mock replicates are near-BINOMIALLY reproducible -- median per-site
# overdispersion ratio ~1 (implied phi ~ 0), 90th pct phi ~ 0.009. theta is therefore well
# determined, which is what justifies a reference this light. Data with genuinely overdispersed
# replicates needs a LOWER ref_df (F(1,4) was right for phi 0.01-0.08); ref_df is exposed so it can
# be re-calibrated. Sanity-check any new dataset by running a within-condition (null) contrast and
# confirming p<0.05 lands near 0.05.
REF_DF = 10
_MIN_SITES_FOR_CALIBRATION = 50   # below this the empirical null median is unreliable -> no scaling
_MIN_CALIBRATION = 1.0            # never DEFLATE the statistic (stay on the conservative side)


def _bb_ll(k, n, mu, theta):
    """Beta-binomial log-likelihood (binomial coefficient dropped -- constant across models)."""
    mu = min(max(mu, _MIN_MU), _MAX_MU)
    a, b = mu * theta, (1.0 - mu) * theta
    return float(np.sum(betaln(k + a, n - k + b) - betaln(a, b)))


def _fit_mu(k, n, theta):
    """MLE of mu for one group with theta fixed (1-D, bounded)."""
    if len(k) == 0:
        return np.nan, 0.0
    tot = float(np.sum(n))
    p0 = float(np.sum(k)) / tot if tot > 0 else 0.5
    if p0 <= _MIN_MU or p0 >= _MAX_MU:          # degenerate: all-0 or all-1, MLE is at the bound
        mu = min(max(p0, _MIN_MU), _MAX_MU)
        return mu, _bb_ll(k, n, mu, theta)
    r = minimize_scalar(lambda x: -_bb_ll(k, n, expit(x), theta),
                        bracket=None, bounds=(logit(_MIN_MU), logit(_MAX_MU)), method="bounded",
                        options={"xatol": 1e-6})
    return float(expit(r.x)), float(-r.fun)


def _fit_theta(k, n, gidx, ngroups=2):
    """MLE of a site's theta with a mu per group -- ONE joint fit over (logit mu0, logit mu1, log theta).

    Profiling mu inside a theta search would nest an optimizer in an optimizer (~5k likelihood evals
    per site, ~80 min on 100k sites); the joint fit is ~1 optimization per site instead.
    """
    m0 = gidx == 0
    def _p0(mask):
        tot = float(np.sum(n[mask]))
        p = float(np.sum(k[mask])) / tot if tot > 0 else 0.5
        return min(max(p, 1e-3), 1 - 1e-3)

    def negll(params):
        mu0, mu1, theta = expit(params[0]), expit(params[1]), float(np.exp(params[2]))
        return -(_bb_ll(k[m0], n[m0], mu0, theta) + _bb_ll(k[~m0], n[~m0], mu1, theta))

    x0 = np.array([logit(_p0(m0)), logit(_p0(~m0)), np.log(100.0)])
    r = minimize(negll, x0, method="L-BFGS-B",
                 bounds=[(logit(_MIN_MU), logit(_MAX_MU)),
                         (logit(_MIN_MU), logit(_MAX_MU)),
                         (_LOG_THETA_LO, _LOG_THETA_HI)])
    return float(np.exp(r.x[2])), float(-r.fun)


def _site_lrt(k, n, gidx, theta):
    """LRT statistic for mu_0 != mu_1 with theta fixed. Returns (mu0, mu1, stat).

    The p-value is NOT computed here: the reference distribution is calibrated across sites (see
    beta_binomial_diff), because chi2(1) is invalid at this sample size.
    """
    m0, m1 = gidx == 0, gidx == 1
    mu0, ll0 = _fit_mu(k[m0], n[m0], theta)
    mu1, ll1 = _fit_mu(k[m1], n[m1], theta)
    _, ll_null = _fit_mu(k, n, theta)
    stat = 2.0 * ((ll0 + ll1) - ll_null)
    return mu0, mu1, max(stat, 0.0)


def parse_site_weight(x):
    """CLI '--site-weight' -> 'auto' or a float. Anything non-numeric (incl. 'auto'/'') -> 'auto'."""
    s = str(x).strip().lower() if x is not None else "auto"
    if s in ("auto", "", "none"):
        return "auto"
    try:
        return float(s)
    except ValueError:
        return "auto"


def beta_binomial_diff(sites, prior_weight=20.0, min_group_samples=2, ref_df=REF_DF,
                       calibrate=False, site_weight="auto"):
    """Replicate-aware differential test for count data, with dispersion shrinkage across sites.

    ``sites``: list of (key, k, n, gidx) where k/n/gidx are equal-length arrays over samples and
    gidx is 0/1 (reference/test). Returns a list of dicts (one per tested site).

    Two passes: (1) per-site theta MLE -> (2) shrink log-theta toward the across-site trend, then LRT
    with the shrunk theta. ``prior_weight`` is how many "pseudo-sites" of prior to mix in: higher =
    more shrinkage toward the global trend (more power, less per-site adaptivity).

    CALIBRATION: the LRT is referenced to F(1, ref_df), NOT chi2(1) -- see the REF_DF note above for
    why and for the validation numbers. ``calibrate=True`` additionally rescales the statistic by the
    observed LRT median (genomic-control style); it is OFF by default because real signal inflates
    that median too, so it over-corrects and destroys power (validated: at 10% true DE it called
    zero sites). Leave it off unless you have a reason.
    """
    prepared = []
    for key, k, n, gidx in sites:
        k = np.asarray(k, dtype=float); n = np.asarray(n, dtype=float); gidx = np.asarray(gidx, dtype=int)
        ok = n > 0
        k, n, gidx = k[ok], n[ok], gidx[ok]
        if (gidx == 0).sum() < min_group_samples or (gidx == 1).sum() < min_group_samples:
            continue
        theta_s, _ = _fit_theta(k, n, gidx, 2)
        prepared.append((key, k, n, gidx, theta_s))
    if not prepared:
        return []

    # Shrink log-theta toward the across-site median (the "trend"); this buys power at small n.
    # The per-site weight w in the shrinkage balance (w*logθ_s + prior_weight*log_prior)/(w+prior_weight)
    # must SCALE WITH how well the site's own dispersion is determined -- i.e. with the number of
    # replicates covering it. With w hardcoded to 1 against prior_weight=20, EVERY site (even in a
    # 400-sample cohort) got only ~5% of its own dispersion and was forced ~95% onto the global
    # median. That is fine when replicates are few and near-binomial (the site's θ ≈ the prior, so
    # shrinkage is inert), but on a large or heterogeneous cohort it crushes genuinely overdispersed
    # sites onto the near-binomial bulk and reads their replicate scatter as signal -> false positives.
    # site_weight="auto" sets w = max(1, N_site - 2) (the residual df for dispersion at that site), so
    # shrinkage automatically fades as the cohort grows; a numeric value forces a fixed w (w=1
    # reproduces the legacy behaviour). Near-binomial null sites are unaffected either way.
    log_thetas = np.array([np.log(p[4]) for p in prepared])
    log_prior = float(np.median(log_thetas))
    auto_w = (site_weight == "auto")
    fixed_w = None if auto_w else max(0.0, float(site_weight))
    out = []
    for key, k, n, gidx, theta_s in prepared:
        w = max(1.0, float(len(k) - 2)) if auto_w else fixed_w
        denom = w + prior_weight
        log_shrunk = (w * np.log(theta_s) + prior_weight * log_prior) / denom if denom > 0 else np.log(theta_s)
        theta = float(np.exp(log_shrunk))
        mu0, mu1, stat = _site_lrt(k, n, gidx, theta)
        out.append({
            "key": key, "mu_reference": mu0, "mu_test": mu1, "delta": mu1 - mu0,
            "theta_site": theta_s, "theta_shrunk": theta, "dispersion": 1.0 / (theta + 1.0),
            "lrt_stat": stat,
            "n_reference": int((gidx == 0).sum()), "n_test": int((gidx == 1).sum()),
            "reads_reference": float(n[gidx == 0].sum()), "reads_test": float(n[gidx == 1].sum()),
        })

    # Empirical-null calibration: scale the LRT so its null bulk matches F(1, ref_df).
    stats_arr = np.array([r["lrt_stat"] for r in out])
    c = 1.0
    if calibrate and len(stats_arr) >= _MIN_SITES_FOR_CALIBRATION:
        med = float(np.median(stats_arr))
        ref_med = float(fdist.ppf(0.5, 1, ref_df))
        if med > 0 and ref_med > 0:
            c = max(med / ref_med, _MIN_CALIBRATION)
    for r, s in zip(out, stats_arr):
        r["calibration_factor"] = c
        r["p_value"] = float(fdist.sf(s / c, 1, ref_df))
    return out


def continuous_diff(values_a, values_b):
    """Between-condition test for a continuous per-sample summary (e.g. median poly(A) tail).

    One value per REPLICATE (never per read -- that would be pseudoreplication). Welch's t-test.
    """
    a, b = np.asarray(values_a, dtype=float), np.asarray(values_b, dtype=float)
    a, b = a[np.isfinite(a)], b[np.isfinite(b)]
    if a.size < 2 or b.size < 2:
        return {"stat": np.nan, "p_value": np.nan, "mean_reference": np.nan,
                "mean_test": np.nan, "delta": np.nan, "n_reference": int(a.size), "n_test": int(b.size)}
    # Welch's t returns pvalue=0.0 (t=+/-inf) when BOTH groups have zero within-group variance but
    # different means (e.g. per-replicate medians [50,50,50] vs [52,52,52]). 0.0 is finite, so the
    # caller's isfinite guard would let it through as the #1 hit despite zero replicate variation.
    # Treat exactly-zero pooled variance as non-significant.
    if a.var(ddof=1) == 0.0 and b.var(ddof=1) == 0.0:
        return {"stat": np.nan, "p_value": 1.0, "mean_reference": float(a.mean()),
                "mean_test": float(b.mean()), "delta": float(b.mean() - a.mean()),
                "n_reference": int(a.size), "n_test": int(b.size)}
    t, p = ttest_ind(b, a, equal_var=False)
    return {"stat": float(t), "p_value": float(p), "mean_reference": float(a.mean()),
            "mean_test": float(b.mean()), "delta": float(b.mean() - a.mean()),
            "n_reference": int(a.size), "n_test": int(b.size)}
