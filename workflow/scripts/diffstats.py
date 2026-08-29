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
from scipy.special import betaln, expit, logit, polygamma
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
# DISPERSION note: per-site theta is estimated with a Cox-Reid ADJUSTED profile likelihood (see
# _fit_theta), which removes the downward dispersion bias of plain ML. WITHOUT that adjustment the
# test was ~2x anti-conservative on even mildly overdispersed cohorts (null type-I 0.11 at phi>=0.005,
# 0.26 at 2v2) while looking calibrated only at phi~0 -- the regime the 3 mock replicates happen to
# sit in. WITH Cox-Reid, at ref_df=10 the null type-I is ~nominal across dispersion (phi=0 -> 0.026;
# phi=0.005/0.02/0.05 -> 0.061/0.052/0.056 at 3v3) and only conservative at higher replicate counts --
# never anti-conservative except the 2v2 corner (0.093). So the fix was the dispersion estimator, NOT
# ref_df (ref_df is still exposed for re-calibration). Sanity-check any new dataset with a
# within-condition (null) contrast and confirm p<0.05 lands near/below 0.05.
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


def _mu_info(k, n, mu, theta):
    """Observed Fisher information of one group's mean mu at fixed theta: j = -d2/dmu2 sum_i ll_i.

    Analytic (trigamma) -- with a=mu*theta, b=(1-mu)*theta and betaln(x,y)=lgamma(x)+lgamma(y)-lgamma(x+y),
    d2 ll/dmu2 = theta^2 * sum[psi1(k+a)+psi1(n-k+b)-psi1(a)-psi1(b)] <= 0, so j = -that >= 0. This is the
    curvature of the mean nuisance parameter that Cox-Reid penalises when estimating theta."""
    mu = min(max(mu, _MIN_MU), _MAX_MU)
    a, b = mu * theta, (1.0 - mu) * theta
    j = theta * theta * float(np.sum(polygamma(1, a) + polygamma(1, b)
                                     - polygamma(1, k + a) - polygamma(1, n - k + b)))
    return max(j, _EPS)


def _fit_theta(k, n, gidx, ngroups=2):
    """Cox-Reid ADJUSTED profile estimate of a site's theta with a mu per group.

    Plain ML of theta jointly with the two group means is biased LOW in dispersion (theta too high):
    it ignores the degrees of freedom spent estimating mu0 and mu1 from few replicates, so at n=3 per
    group the between-condition LRT is ~2x anti-conservative once replicates are even mildly
    overdispersed. Exactly as edgeR/DSS do, we maximise the Cox-Reid adjusted profile log-likelihood
    over theta:  APL(theta) = l(theta, mu_hat(theta)) - 1/2 * log|j_mu(mu_hat; theta)|,  where j_mu is
    the observed information of the (nuisance) group means; with two independent groups |j_mu| factors
    as j0*j1. This is a bounded 1-D search over log theta (mu profiled analytically per group), so it
    is only a few group-mean fits per site -- fine at the site counts a between-condition test sees."""
    m0 = gidx == 0
    k0, n0, k1, n1 = k[m0], n[m0], k[~m0], n[~m0]

    def neg_apl(log_theta):
        theta = float(np.exp(log_theta))
        mu0, ll0 = _fit_mu(k0, n0, theta)
        mu1, ll1 = _fit_mu(k1, n1, theta)
        cr = 0.5 * (np.log(_mu_info(k0, n0, mu0, theta)) + np.log(_mu_info(k1, n1, mu1, theta)))
        return -((ll0 + ll1) - cr)

    r = minimize_scalar(neg_apl, bounds=(_LOG_THETA_LO, _LOG_THETA_HI), method="bounded",
                        options={"xatol": 1e-3})
    theta = float(np.exp(r.x))
    _, ll0 = _fit_mu(k0, n0, theta)
    _, ll1 = _fit_mu(k1, n1, theta)
    return theta, float(ll0 + ll1)


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
        # Flag sites that carry NO usable dispersion information: a degenerate pooled mean (all reads
        # modified or all unmodified -- e.g. the all-zero sites that flood the universe when the site
        # filter is disabled) leaves theta unidentified, so the optimiser pins it at a bound. Such sites
        # must NOT enter the across-site prior median, or they collapse the prior and destroy a genuine
        # effect (prior-collapse finding).
        #
        # MAJOR-4: exclude only the LOWER bound (theta=1e-2, extreme overdispersion -- the pathological
        # pin an unidentified/degenerate site runs to). A site at the UPPER bound (theta=1e6) is
        # genuinely NEAR-BINOMIAL (under-dispersed): that is real, legitimate low-dispersion evidence and
        # MUST stay in the prior median. Dropping it (the old both-bounds test) biased the prior toward
        # overdispersion, so every other site was over-shrunk to a too-dispersed prior and real effects
        # were washed out -- and on a near-binomial batch it emptied the informative set below the
        # shrink threshold, disabling shrinkage entirely.
        _mu_hat = float(k.sum() / n.sum()) if n.sum() > 0 else 0.0
        _lt = float(np.log(theta_s))
        _at_bound = abs(_lt - _LOG_THETA_LO) < 1e-2
        _degenerate = (_mu_hat <= 1e-9) or (_mu_hat >= 1.0 - 1e-9)
        informative = not (_at_bound or _degenerate)
        prepared.append((key, k, n, gidx, theta_s, informative))
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
    #
    # KNOWN LIMITATION (finding S1): because the prior is the across-site median dispersion of the batch
    # tested TOGETHER, a site's shrunk theta -- and therefore its p-value -- depends on its batch-mates.
    # At the pipeline's design point (n=3 vs 3) w = max(1, N-2) = 4 against prior_weight=20, so a site
    # keeps only ~17% of its own dispersion and "auto" is nearly inert -- fine for near-binomial
    # replicates (the calibrated, validated case) but it means a genuinely overdispersed site in an
    # otherwise near-binomial batch is under-shrunk-corrected and can inflate type-I error. This is the
    # standard cost of empirical-Bayes dispersion shrinkage (edgeR/DSS share it); it is documented rather
    # than re-engineered so the calibration in the REF_DF note stays valid. Lower prior_weight (or a
    # larger site_weight) to trust each site's own dispersion more when your replicates are heterogeneous.
    # Prior trend from INFORMATIVE sites only (degenerate/bound-pinned sites excluded above). If too few
    # remain to estimate a trend, do not shrink at all -- fall back to each site's own theta, which is
    # strictly safer than shrinking toward a prior built from junk.
    _info_lt = np.array([np.log(p[4]) for p in prepared if p[5]])
    _n_info = len(_info_lt)
    # MAJOR-4: the old `>= 10` was a HARD CLIFF -- 9 informative sites shrank NOTHING, 10 shrank FULLY,
    # so one borderline site flipped every p-value in the batch. Replace with a floor + smooth ramp: a
    # median needs a few sites to mean anything (floor = 3), and the prior earns full weight only once it
    # rests on enough sites. `prior_conf` ramps 0->1 linearly across [3, 10] informative sites and scales
    # the prior's pseudo-count, so shrinkage fades in gradually instead of switching on.
    _PRIOR_FLOOR, _PRIOR_FULL = 3, 10
    _shrink = _n_info >= _PRIOR_FLOOR
    log_prior = float(np.median(_info_lt)) if _shrink else 0.0
    prior_conf = min(1.0, (_n_info - _PRIOR_FLOOR) / float(_PRIOR_FULL - _PRIOR_FLOOR)) if _shrink else 0.0
    eff_prior_weight = prior_weight * prior_conf
    auto_w = (site_weight == "auto")
    fixed_w = None if auto_w else max(0.0, float(site_weight))
    out = []
    for key, k, n, gidx, theta_s, _informative in prepared:
        if not _shrink or eff_prior_weight <= 0:
            log_shrunk = np.log(theta_s)
        else:
            w = max(1.0, float(len(k) - 2)) if auto_w else fixed_w
            denom = w + eff_prior_weight
            log_shrunk = (w * np.log(theta_s) + eff_prior_weight * log_prior) / denom if denom > 0 else np.log(theta_s)
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


_PFLOOR_CACHE = {}


def design_pfloor(n_ref, n_test, ref_df=REF_DF, prior_weight=20.0, coverage=10000):
    """Minimum achievable p-value of beta_binomial_diff at a given replicate design -- the STRUCTURAL
    ceiling on significance set by the replicate count alone, independent of the data.

    Computed from the maximum conceivable effect (0% modified in every reference replicate, 100% in
    every test replicate) at saturating coverage. This is why a between-condition table can be all-
    non-significant purely because of small n per group: e.g. ~0.035 at 2v2, ~0.0066 at 4v4, ~0.0007
    at 8v8. Under Benjamini-Hochberg a p-floor F only yields a rejection when a fraction >= F/alpha of
    the whole family sits AT the floor, so at 2v2 (F~0.035, alpha 0.05) ~71% of tests must be
    maximum-effect before ANY call is made. Returns NaN if the design is below the minimum group size.
    """
    key = (int(n_ref), int(n_test), int(ref_df), float(prior_weight), int(coverage))
    if key in _PFLOOR_CACHE:
        return _PFLOOR_CACHE[key]
    nr, nt = int(n_ref), int(n_test)
    if nr < 2 or nt < 2:
        _PFLOOR_CACHE[key] = float("nan")
        return float("nan")
    k = np.array([0.0] * nr + [float(coverage)] * nt, dtype=float)
    n = np.array([float(coverage)] * (nr + nt), dtype=float)
    gidx = np.array([0] * nr + [1] * nt, dtype=int)
    res = beta_binomial_diff([("floor", k, n, gidx)], prior_weight=prior_weight,
                             min_group_samples=2, ref_df=ref_df, calibrate=False)
    val = float(res[0]["p_value"]) if res else float("nan")
    _PFLOOR_CACHE[key] = val
    return val


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
