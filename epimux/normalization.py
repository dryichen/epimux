"""Normalization, including the spike-in case.

Standard size-factor methods (median-of-ratios, TMM, quantile) all assume that
*most features do not change*.  When that assumption breaks -- a global gain or
loss of a mark -- they silently absorb the very effect you are trying to measure
and the result is a flat, confidently-wrong null.

Two consequences are implemented here:

* :func:`spike_in_factors` uses exogenous chromatin (or any external reference)
  to derive size factors that survive a global shift;
* :func:`assess_global_shift` estimates whether a global shift is *plausible*
  from the data alone, so a study without spike-ins at least knows whether its
  normalization assumption is at risk.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from .utils import get_logger

LOG = get_logger()

__all__ = ["median_of_ratios", "tmm", "quantile_normalize", "spike_in_factors",
           "apply_factors", "assess_global_shift", "reference_normalize"]


# --------------------------------------------------------------------------
def median_of_ratios(counts: pd.DataFrame, min_count: int = 1) -> pd.Series:
    """DESeq2-style size factors (geometric-mean reference)."""
    X = counts.to_numpy(dtype=float)
    keep = (X >= min_count).all(axis=1)
    if keep.sum() == 0:
        LOG.warning("no feature detected in every sample; falling back to library size")
        lib = X.sum(0)
        return pd.Series(lib / np.exp(np.mean(np.log(lib))), index=counts.columns)
    logX = np.log(X[keep])
    ref = logX.mean(axis=1, keepdims=True)
    sf = np.exp(np.median(logX - ref, axis=0))
    return pd.Series(sf / np.exp(np.mean(np.log(sf))), index=counts.columns)


def tmm(counts: pd.DataFrame, ref_col: str | None = None,
        log_ratio_trim: float = 0.3, sum_trim: float = 0.05) -> pd.Series:
    """edgeR-style trimmed mean of M-values."""
    X = counts.to_numpy(dtype=float)
    lib = X.sum(0)
    if ref_col is None:
        # sample whose upper-quartile is closest to the mean upper-quartile
        uq = np.array([np.percentile(X[:, j] / lib[j], 75) for j in range(X.shape[1])])
        ref_idx = int(np.argmin(np.abs(uq - uq.mean())))
    else:
        ref_idx = list(counts.columns).index(ref_col)
    r = X[:, ref_idx] / lib[ref_idx]
    factors = []
    for j in range(X.shape[1]):
        o = X[:, j] / lib[j]
        ok = (X[:, j] > 0) & (X[:, ref_idx] > 0)
        if ok.sum() < 10:
            factors.append(1.0)
            continue
        M = np.log2(o[ok] / r[ok])
        A = 0.5 * np.log2(o[ok] * r[ok])
        w = (lib[j] - X[ok, j]) / (lib[j] * X[ok, j]) + \
            (lib[ref_idx] - X[ok, ref_idx]) / (lib[ref_idx] * X[ok, ref_idx])
        mlo, mhi = np.quantile(M, [log_ratio_trim, 1 - log_ratio_trim])
        alo, ahi = np.quantile(A, [sum_trim, 1 - sum_trim])
        sel = (M >= mlo) & (M <= mhi) & (A >= alo) & (A <= ahi) & np.isfinite(w) & (w > 0)
        factors.append(2 ** (np.sum(M[sel] / w[sel]) / np.sum(1 / w[sel])) if sel.sum() else 1.0)
    f = np.array(factors) * lib
    return pd.Series(f / np.exp(np.mean(np.log(f))), index=counts.columns)


def quantile_normalize(values: pd.DataFrame) -> pd.DataFrame:
    """Force every sample to share one distribution (use with care: this is the
    most aggressive way to erase a genuine global shift)."""
    X = values.to_numpy(dtype=float)
    order = np.argsort(X, axis=0)
    ranks = np.empty_like(order)
    for j in range(X.shape[1]):
        ranks[order[:, j], j] = np.arange(X.shape[0])
    mean_sorted = np.sort(X, axis=0).mean(axis=1)
    return pd.DataFrame(mean_sorted[ranks], index=values.index, columns=values.columns)


# --------------------------------------------------------------------------
def spike_in_factors(spike_counts: pd.Series | dict,
                     target_lib: pd.Series | dict | None = None) -> pd.Series:
    """Size factors from exogenous spike-in reads.

    ``spike_counts`` : reads mapping to the spike-in genome, per sample.
    ``target_lib``   : optional total reads on the target genome; when given,
    factors are corrected for differences in sequencing depth so that they
    reflect *biological* occupancy rather than how deeply each library was run.

    Unlike median-of-ratios these factors are unaffected by a genome-wide gain
    or loss of the mark -- which is precisely why spike-ins are required to
    claim a global change.
    """
    s = pd.Series(spike_counts, dtype=float)
    if (s <= 0).any():
        raise ValueError("spike-in counts must be positive")
    f = s / np.exp(np.mean(np.log(s)))
    if target_lib is not None:
        t = pd.Series(target_lib, dtype=float).reindex(f.index)
        depth = t / np.exp(np.mean(np.log(t)))
        f = f / depth
        f = f / np.exp(np.mean(np.log(f)))
    LOG.info("spike-in size factors: " + ", ".join(f"{k}={v:.3f}" for k, v in f.items()))
    return f


def reference_normalize(counts: pd.DataFrame, reference_features) -> pd.Series:
    """Size factors from an internal set assumed invariant (e.g. constitutive
    CTCF sites, housekeeping promoters).  A poor substitute for spike-ins, but
    far better than assuming the whole genome is invariant."""
    idx = counts.index.intersection(pd.Index(reference_features))
    if len(idx) < 20:
        raise ValueError(f"need >=20 reference features, got {len(idx)}")
    sub = counts.loc[idx].to_numpy(dtype=float)
    lib = sub.sum(0)
    f = lib / np.exp(np.mean(np.log(lib)))
    LOG.info(f"reference-set normalization on {len(idx):,} features")
    return pd.Series(f, index=counts.columns)


def apply_factors(counts: pd.DataFrame, factors: pd.Series) -> pd.DataFrame:
    f = factors.reindex(counts.columns)
    if f.isna().any():
        raise KeyError(f"missing size factors for {list(f.index[f.isna()])}")
    return counts.divide(f, axis=1)


# --------------------------------------------------------------------------
def assess_global_shift(counts: pd.DataFrame, contrast, frip: dict | None = None,
                        quantiles=(0.25, 0.5, 0.75, 0.9)) -> dict:
    """Is a genome-wide shift plausible, and would normalization hide it?

    Compares raw per-sample signal distributions (depth-normalized only) between
    the two groups.  A consistent offset across quantiles is what a true global
    change looks like -- and is exactly what size-factor normalization removes.
    """
    ref, test = contrast.ref_samples, contrast.test_samples
    X = counts.to_numpy(dtype=float)
    lib = X.sum(0)
    cpm = X / lib * 1e6
    cols = list(counts.columns)
    ri = [cols.index(s) for s in ref]
    ti = [cols.index(s) for s in test]
    out = {}
    for q in quantiles:
        a = np.mean([np.quantile(cpm[:, j], q) for j in ri])
        b = np.mean([np.quantile(cpm[:, j], q) for j in ti])
        out[f"q{q}"] = float(np.log2((b + 1e-9) / (a + 1e-9)))
    shifts = np.array(list(out.values()))
    consistent = bool(np.all(shifts > 0.15) or np.all(shifts < -0.15))
    out["consistent_shift"] = consistent
    out["mean_log2_shift"] = float(shifts.mean())
    if frip:
        fr = np.mean([frip[s] for s in ref if s in frip])
        ft = np.mean([frip[s] for s in test if s in frip])
        out["frip_ref"], out["frip_test"] = float(fr), float(ft)
        out["frip_ratio"] = float(fr / ft) if ft else np.inf
    if consistent:
        LOG.warning(
            "a consistent shift is present across all quantiles "
            f"(mean log2 {shifts.mean():+.2f}). Size-factor normalization will remove it. "
            "Without spike-ins the global magnitude cannot be recovered; only the "
            "relative redistribution of signal is interpretable.")
    out["interpretation"] = (
        "global shift plausible -- spike-ins required to quantify it; "
        "report relative redistribution only" if consistent else
        "no consistent global shift detected; standard normalization is reasonable")
    return out
