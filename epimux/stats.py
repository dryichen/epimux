"""Differential engines, one per data type.

Every engine returns the same schema::

    baseMean  log2FC  stat  pvalue  padj

and ``log2FC`` is **always** ``log2(test / ref)``.  This is enforced by passing a
:class:`~epimux.utils.Contrast` rather than a factor, and is re-verified against
raw normalised values by :func:`epimux.audit.check_direction`.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats as ss

from .utils import Contrast, cpm, get_logger, log2cpm

LOG = get_logger()

__all__ = ["deseq2_de", "moderated_t_de", "methylation_de", "bh_fdr", "RESULT_COLS"]

RESULT_COLS = ["baseMean", "log2FC", "stat", "pvalue", "padj"]


def bh_fdr(p: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg, NaN-safe."""
    p = np.asarray(p, dtype=float)
    out = np.full_like(p, np.nan)
    ok = np.isfinite(p)
    if ok.sum() == 0:
        return out
    q = p[ok]
    n = q.size
    order = np.argsort(q)
    ranked = q[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    adj = np.empty(n)
    adj[order] = np.clip(ranked, 0, 1)
    out[ok] = adj
    return out


# --------------------------------------------------------------------------
# count data  (ATAC / ChIP / RNA)
# --------------------------------------------------------------------------
def deseq2_de(counts: pd.DataFrame, contrast: Contrast,
              min_count: int = 10, min_frac: float = 0.5,
              covariates: pd.DataFrame | None = None,
              n_cpus: int = 4) -> pd.DataFrame:
    """Negative-binomial DE via PyDESeq2, with an exact-sign contrast.

    ``counts`` : features x samples (raw integer counts).
    ``covariates`` : optional samples x variables frame (e.g. batch, sex, donor).
    Each column is added to the design as ``~ cov1 + cov2 + condition`` so the
    contrast is adjusted for it; the tested coefficient is still condition.
    Features failing the expression filter are returned with NaN statistics so
    the output always aligns with the input index.
    """
    samples = contrast.ref_samples + contrast.test_samples
    missing = [s for s in samples if s not in counts.columns]
    if missing:
        raise KeyError(f"samples not in count matrix: {missing}")
    sub = counts.loc[:, samples]
    keep = (sub >= min_count).sum(axis=1) >= max(2, int(np.floor(min_frac * len(samples))))
    res = pd.DataFrame(np.nan, index=counts.index, columns=RESULT_COLS)
    if keep.sum() == 0:
        LOG.warning("no features pass the expression filter")
        return res

    try:
        from pydeseq2.dds import DeseqDataSet
        from pydeseq2.ds import DeseqStats
    except ImportError as e:  # pragma: no cover
        raise ImportError("pydeseq2 is required for count-based DE") from e

    meta = contrast.design().loc[samples]
    design = "~condition"
    if covariates is not None:
        cov = covariates.reindex(samples)
        if cov.isna().any().any():
            raise ValueError("covariates must cover every sample in the contrast")
        for c in cov.columns:
            meta[c] = cov[c].astype(str).to_numpy()
        design = "~" + " + ".join(list(cov.columns) + ["condition"])
        LOG.info(f"adjusted design: {design}")
    # PyDESeq2 wants samples x genes
    dds = DeseqDataSet(counts=sub.loc[keep].T.astype(int), metadata=meta,
                       design=design, refit_cooks=True, quiet=True)
    dds.deseq2()
    st = DeseqStats(dds, contrast=["condition", contrast.test, contrast.ref], quiet=True)
    st.summary()
    r = st.results_df
    res.loc[r.index, "baseMean"] = r["baseMean"].to_numpy()
    res.loc[r.index, "log2FC"] = r["log2FoldChange"].to_numpy()
    res.loc[r.index, "stat"] = r["stat"].to_numpy()
    res.loc[r.index, "pvalue"] = r["pvalue"].to_numpy()
    res.loc[r.index, "padj"] = r["padj"].to_numpy()
    res.attrs["engine"] = "pydeseq2"
    res.attrs["design"] = design
    res.attrs["contrast"] = repr(contrast)
    return res


# --------------------------------------------------------------------------
# continuous signal  (bigWig-derived) -- moderated t, limma style
# --------------------------------------------------------------------------
def moderated_t_de(signal: pd.DataFrame, contrast: Contrast,
                   min_var: float = 1e-8) -> pd.DataFrame:
    """Empirical-Bayes moderated t-test on continuous signal.

    This is the right test for averaged bigWig signal, but note the package
    deliberately warns when you use it for a genotype contrast: averaged signal
    tracks are markedly less sensitive than read counts, and a null result from
    this engine should never be reported as "no change" without a positive
    control (see :func:`epimux.audit.positive_control`).
    """
    A = signal.loc[:, contrast.ref_samples].to_numpy(dtype=float)
    B = signal.loc[:, contrast.test_samples].to_numpy(dtype=float)
    na, nb = A.shape[1], B.shape[1]
    ok = np.isfinite(A).all(1) & np.isfinite(B).all(1)
    res = pd.DataFrame(np.nan, index=signal.index, columns=RESULT_COLS)
    if ok.sum() == 0:
        return res
    a, b = A[ok], B[ok]
    ma, mb = a.mean(1), b.mean(1)
    diff = mb - ma                      # test - ref  (log space -> log2FC)
    df = na + nb - 2
    s2 = ((a - ma[:, None]) ** 2).sum(1) + ((b - mb[:, None]) ** 2).sum(1)
    s2 = np.maximum(s2 / df, min_var)
    # empirical Bayes shrinkage of the variance toward the trend (limma-lite)
    log_s2 = np.log(s2)
    s2_prior = np.exp(np.median(log_s2))
    d0 = 4.0                            # moderate prior df
    s2_post = (d0 * s2_prior + df * s2) / (d0 + df)
    se = np.sqrt(s2_post * (1.0 / na + 1.0 / nb))
    t = diff / se
    p = 2 * ss.t.sf(np.abs(t), df + d0)
    res.loc[ok, "baseMean"] = (ma + mb) / 2
    res.loc[ok, "log2FC"] = diff
    res.loc[ok, "stat"] = t
    res.loc[ok, "pvalue"] = p
    res.loc[ok, "padj"] = bh_fdr(p)
    res.attrs["engine"] = "moderated_t"
    res.attrs["contrast"] = repr(contrast)
    res.attrs["warning"] = ("averaged signal is less sensitive than counts; "
                            "validate any null with a positive control")
    return res


# --------------------------------------------------------------------------
# methylation  (rates with coverage)
# --------------------------------------------------------------------------
def methylation_de(meth: pd.DataFrame, cov: pd.DataFrame, contrast: Contrast,
                   min_cov: int = 10, test: str = "beta") -> pd.DataFrame:
    """Differential methylation on rates, coverage-weighted.

    ``meth`` : methylated read counts, features x samples.
    ``cov``  : total read counts, same shape.
    ``log2FC`` is replaced by an **absolute rate difference** (``test - ref``);
    the column keeps the name ``log2FC`` for schema consistency but the
    ``value_kind`` attr records that it is a difference, not a ratio, and the
    plotting/threshold code honours that.
    """
    rs, ts = contrast.ref_samples, contrast.test_samples
    mA, cA = meth[rs].to_numpy(float), cov[rs].to_numpy(float)
    mB, cB = meth[ts].to_numpy(float), cov[ts].to_numpy(float)
    res = pd.DataFrame(np.nan, index=meth.index, columns=RESULT_COLS)

    with np.errstate(invalid="ignore", divide="ignore"):
        rA = np.where(cA >= min_cov, mA / cA, np.nan)
        rB = np.where(cB >= min_cov, mB / cB, np.nan)
    ok = (np.isfinite(rA).sum(1) >= 2) & (np.isfinite(rB).sum(1) >= 2)
    if ok.sum() == 0:
        return res
    a, b = rA[ok], rB[ok]
    mean_a = np.nanmean(a, 1)
    mean_b = np.nanmean(b, 1)
    diff = mean_b - mean_a

    if test == "beta":
        # arcsine (variance-stabilising) transform then Welch t
        ta, tb = np.arcsin(np.sqrt(np.clip(a, 0, 1))), np.arcsin(np.sqrt(np.clip(b, 0, 1)))
        stat, p = ss.ttest_ind(tb, ta, axis=1, nan_policy="omit", equal_var=False)
    elif test == "fisher":
        # pooled-count Fisher exact per feature (fast approximation via chi2)
        sm_a, sc_a = np.nansum(mA[ok], 1), np.nansum(cA[ok], 1)
        sm_b, sc_b = np.nansum(mB[ok], 1), np.nansum(cB[ok], 1)
        tbl = np.stack([sm_b, sc_b - sm_b, sm_a, sc_a - sm_a], axis=1)
        with np.errstate(invalid="ignore", divide="ignore"):
            n = tbl.sum(1)
            exp = (tbl[:, 0] + tbl[:, 1]) * (tbl[:, 0] + tbl[:, 2]) / n
            stat = (tbl[:, 0] - exp) ** 2 / np.maximum(exp, 1e-9)
            p = ss.chi2.sf(stat, 1)
    else:
        raise ValueError("test must be 'beta' or 'fisher'")

    stat = np.asarray(np.ma.filled(stat, np.nan), dtype=float)
    p = np.asarray(np.ma.filled(p, np.nan), dtype=float)
    res.loc[ok, "baseMean"] = (mean_a + mean_b) / 2
    res.loc[ok, "log2FC"] = diff
    res.loc[ok, "stat"] = stat
    res.loc[ok, "pvalue"] = p
    res.loc[ok, "padj"] = bh_fdr(p)
    res.attrs["engine"] = f"methylation:{test}"
    res.attrs["value_kind"] = "difference"   # NOT a log ratio
    res.attrs["contrast"] = repr(contrast)
    return res
