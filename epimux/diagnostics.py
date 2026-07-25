"""Diagnostics that localise *which* sample or model assumption is the problem.

`audit` answers "is this result trustworthy?".  This module answers the follow-up:
"which replicate is dragging it, and is the model even appropriate?".

* :func:`outlier_replicates` — ranks samples by how far they sit from their own
  group. A single deviant replicate is the usual cause of a non-empty null
  contrast, and it is invisible in a group-mean comparison.
* :func:`pvalue_diagnostic` — the shape of the p-value histogram is the fastest
  test of whether a differential model fits. Uniform-plus-spike is healthy;
  U-shaped or hump-shaped means the variance model is wrong and the FDR is not
  what it claims.
* :func:`confounding_check` — is the contrast confounded with a covariate?
  If every KO was processed in batch 2, no statistical adjustment can separate
  them, and the design must be reported as confounded.
* :func:`power_analysis` — how many replicates would be needed for the observed
  effect size, and what effect size is detectable with the ones you have.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats as ss

from .audit import AuditResult
from .utils import get_logger, log2cpm

LOG = get_logger()

__all__ = ["outlier_replicates", "pvalue_diagnostic", "confounding_check",
           "power_analysis", "detectable_effect"]


# --------------------------------------------------------------------------
def outlier_replicates(values: pd.DataFrame, groups: dict,
                       max_drop: float = 0.10, z_threshold: float = 2.0,
                       method: str = "spearman") -> AuditResult:
    """Find replicates that sit far from the rest of their own group.

    For each sample the mean correlation to its group-mates is computed, then
    two criteria are applied:

    * **drop** -- how far the sample sits below the *median of the others*.
      This is the criterion that works at small n and is the primary test.
    * **z-score** -- only applied when a group has >= 5 replicates.

    A z-score alone is useless here: with ``n`` samples the most extreme
    z (population SD) is bounded by ``(n-1)/sqrt(n)``, which is 1.15 at n=3 and
    1.50 at n=4 -- a threshold of 2 can never be reached, so a deviant replicate
    in a typical 3-replicate design would always be missed.
    """
    rows = []
    for g, samples in groups.items():
        samples = [s for s in samples if s in values.columns]
        if len(samples) < 3:
            continue
        M = values[samples].to_numpy(dtype=float)
        ok = np.isfinite(M).all(axis=1)
        M = M[ok]
        if M.shape[0] < 50:
            continue
        # .to_numpy() can return a read-only view on newer pandas/numpy, and
        # np.corrcoef output is not guaranteed writable either -- copy before
        # touching the diagonal.
        C = np.array(pd.DataFrame(M, columns=samples).corr(method=method).to_numpy()
                     if method == "spearman" else np.corrcoef(M.T),
                     dtype=float, copy=True)
        np.fill_diagonal(C, np.nan)
        mean_r = np.nanmean(C, axis=1)
        n = len(samples)
        mu, sd = mean_r.mean(), mean_r.std(ddof=0)
        for k, (s, r) in enumerate(zip(samples, mean_r)):
            others = np.delete(mean_r, k)
            drop = float(np.median(others) - r)
            z = float((r - mu) / sd) if sd > 1e-9 else 0.0
            rows.append({"group": g, "sample": s, "mean_corr": float(r),
                         "drop_vs_others": drop, "z": z, "n_replicates": n})
    if not rows:
        return AuditResult("outlier_replicates", "warn",
                           "need >=3 replicates in at least one group")
    df = pd.DataFrame(rows).sort_values("drop_vs_others", ascending=False)
    flag = df["drop_vs_others"] > max_drop
    flag |= (df["n_replicates"] >= 5) & (df["z"] < -z_threshold)
    bad = df[flag]
    detail = {"table": df.to_dict("records"),
              "z_bound_for_n": {int(n): float((n - 1) / np.sqrt(n))
                                for n in sorted(df["n_replicates"].unique())}}
    if len(bad):
        worst = bad.iloc[0]
        return AuditResult(
            "outlier_replicates", "warn",
            f"{len(bad)} deviant replicate(s); worst is {worst['sample']} "
            f"(mean r = {worst['mean_corr']:.3f}, {worst['drop_vs_others']:.3f} below its "
            "group-mates). Re-run the contrast without it to check the result is not "
            "driven by one sample.", detail)
    return AuditResult("outlier_replicates", "pass",
                       f"no replicate sits more than {max_drop:.2f} below its group-mates "
                       f"(min mean r = {df['mean_corr'].min():.3f})", detail)


# --------------------------------------------------------------------------
def pvalue_diagnostic(result: pd.DataFrame, n_bins: int = 20) -> AuditResult:
    """Classify the p-value histogram: healthy, conservative, or mis-specified.

    Under the null p-values are uniform.  Real data should be uniform with a
    spike near zero.  A hump in the middle or a rise toward one means the model
    is wrong (usually variance under- or over-estimated), and the reported FDR
    cannot be trusted.
    """
    p = result["pvalue"].dropna().to_numpy()
    if p.size < 200:
        return AuditResult("pvalue_diagnostic", "warn",
                           f"only {p.size} p-values; histogram shape is uninformative")
    hist, _ = np.histogram(p, bins=n_bins, range=(0, 1))
    frac = hist / hist.sum()
    first, last = frac[0], frac[-1]
    flat = float(np.mean(frac[n_bins // 2:]))          # right half ~ null level
    mid_max = float(frac[2:-2].max()) if n_bins > 6 else flat
    # rough null proportion (Storey-style, conservative)
    pi0 = min(1.0, float(np.mean(p > 0.5) * 2))
    detail = {"first_bin": float(first), "last_bin": float(last),
              "flat_level": flat, "pi0": pi0}
    if last > 1.6 * flat:
        return AuditResult("pvalue_diagnostic", "fail",
                           f"p-values pile up near 1 (last bin {last:.3f} vs flat {flat:.3f}) — "
                           "the variance model is mis-specified; FDR is not interpretable", detail)
    if mid_max > 1.8 * flat and first < mid_max:
        return AuditResult("pvalue_diagnostic", "fail",
                           "hump-shaped p-value histogram — model mis-specification or "
                           "unmodelled batch structure; FDR is not interpretable", detail)
    if first < flat * 1.1:
        return AuditResult("pvalue_diagnostic", "pass",
                           f"p-values are essentially uniform (pi0 ~ {pi0:.2f}) — "
                           "consistent with little or no true signal", detail)
    return AuditResult("pvalue_diagnostic", "pass",
                       f"healthy shape: spike near zero ({first:.3f}) over a flat null "
                       f"({flat:.3f}); estimated null fraction {pi0:.2f}", detail)


# --------------------------------------------------------------------------
def confounding_check(contrast, covariates: pd.DataFrame) -> AuditResult:
    """Is the contrast confounded with any covariate?"""
    samples = contrast.ref_samples + contrast.test_samples
    cov = covariates.reindex(samples)
    cond = pd.Series(["ref"] * len(contrast.ref_samples) + ["test"] * len(contrast.test_samples),
                     index=samples)
    problems, detail = [], {}
    for c in cov.columns:
        tab = pd.crosstab(cond, cov[c].astype(str))
        detail[c] = tab.to_dict()
        # perfectly confounded: no level shared between the two conditions
        shared = (tab > 0).sum(axis=0)
        if (shared >= 2).sum() == 0:
            problems.append(f"'{c}' is perfectly confounded with the contrast")
        elif (tab > 0).all(axis=0).sum() < tab.shape[1] / 2:
            problems.append(f"'{c}' is largely imbalanced across the contrast")
    if any("perfectly" in p for p in problems):
        return AuditResult("confounding_check", "fail",
                           "; ".join(problems) + " — no adjustment can separate them; "
                           "report the design as confounded", detail)
    if problems:
        return AuditResult("confounding_check", "warn", "; ".join(problems), detail)
    return AuditResult("confounding_check", "pass",
                       f"covariates {list(cov.columns)} are balanced across the contrast", detail)


# --------------------------------------------------------------------------
def power_analysis(result: pd.DataFrame, n_per_group: int,
                   fc: float = 1.5, fdr: float = 0.1,
                   n_range=(2, 3, 4, 5, 6, 8, 10)) -> pd.DataFrame:
    """Replicates needed to detect the observed effect sizes.

    Uses the dispersion implied by the fitted results: for each candidate n, the
    two-sample t power at the target log2 fold-change, given the observed
    within-group variability.
    """
    d = result.dropna(subset=["log2FC", "stat", "pvalue"])
    if len(d) < 50:
        raise ValueError("too few results to estimate power")
    # back out an effective SD from |log2FC| / |t| on the tested features
    with np.errstate(invalid="ignore", divide="ignore"):
        se = np.abs(d["log2FC"] / d["stat"])
    se = se[np.isfinite(se) & (se > 0)]
    sd = float(np.median(se) * np.sqrt(n_per_group / 2))   # se = sd*sqrt(2/n)
    target = np.log2(fc)
    n_tests = int(d["pvalue"].notna().sum())
    alpha = fdr / max(n_tests, 1) * max(int((d["padj"] < fdr).sum()), 1)   # BH-ish
    alpha = float(np.clip(alpha, 1e-8, 0.05))
    rows = []
    for n in n_range:
        semi = sd * np.sqrt(2 / n)
        ncp = target / semi
        df = 2 * n - 2
        crit = ss.t.ppf(1 - alpha / 2, df)
        power = float(ss.nct.sf(crit, df, ncp) + ss.nct.cdf(-crit, df, ncp))
        rows.append({"n_per_group": n, "power": power,
                     "detectable_log2FC_at_80pct": float(
                         (ss.t.ppf(1 - alpha / 2, df) + ss.t.ppf(0.8, df)) * semi)})
    out = pd.DataFrame(rows)
    out.attrs["implied_sd"] = sd
    out.attrs["alpha_used"] = alpha
    LOG.info(f"power analysis: implied per-feature SD = {sd:.3f} (log2 units), "
             f"effective alpha = {alpha:.2e}")
    return out


def detectable_effect(result: pd.DataFrame, n_per_group: int,
                      power: float = 0.8, fdr: float = 0.1) -> float:
    """Smallest fold-change detectable at the given power with the current n."""
    tab = power_analysis(result, n_per_group, fdr=fdr, n_range=(n_per_group,))
    return float(2 ** tab["detectable_log2FC_at_80pct"].iloc[0]) if power == 0.8 else float("nan")
