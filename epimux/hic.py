"""Hi-C analytics beyond compartments and insulation.

Contact maps are usually the deepest, least normalization-fragile assay in a
multi-omic study, and they need no spike-in.  When a ChIP is too shallow to
quantify a binding change, the functional consequence is often still measurable
here -- boundary insulation, loop strength, compartment segregation.

Everything takes a cooler URI or an :class:`~epimux.assays.HiCAssay`.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from .utils import as_intervals, get_logger

LOG = get_logger()

__all__ = ["contact_decay", "extrusion_shoulder", "saddle", "compartment_strength",
           "pileup", "apa", "boundary_strength", "differential_insulation"]


def _clr(source, resolution):
    import cooler
    if isinstance(source, str):
        return cooler.Cooler(f"{source}::/resolutions/{resolution}"
                             if "::" not in source else source)
    return source


def _view(clr):
    import bioframe
    return bioframe.make_viewframe(
        [(c, 0, clr.chromsizes[c]) for c in clr.chromnames
         if c in clr.chromsizes and str(c).startswith("chr") and c not in ("chrM", "chrY")])


# --------------------------------------------------------------------------
def contact_decay(source, resolution: int = 20_000) -> pd.DataFrame:
    """P(s): contact probability vs genomic separation."""
    import cooltools
    clr = _clr(source, resolution)
    cvd = cooltools.expected_cis(clr, view_df=_view(clr), smooth=True,
                                 aggregate_smoothed=True, nproc=1)
    cvd = cvd[cvd["dist"] > 0].copy()
    col = "balanced.avg.smoothed.agg" if "balanced.avg.smoothed.agg" in cvd else "balanced.avg"
    out = (cvd.groupby("dist")[col].mean().reset_index()
              .rename(columns={col: "contact"}))
    out["s"] = out["dist"] * resolution
    return out.dropna()


def extrusion_shoulder(decay: pd.DataFrame, s_min: int = 50_000,
                       s_max: int = 5_000_000) -> pd.DataFrame:
    """Log-derivative of P(s).

    Loop extrusion produces a characteristic shoulder; its depth and position
    track processivity, so comparing derivatives between conditions is far more
    informative than overlaying the raw curves (where large differences look
    tiny on a log-log plot).
    """
    d = decay[(decay["s"] >= s_min) & (decay["s"] <= s_max)].copy()
    d = d[d["contact"] > 0]
    ls, lc = np.log10(d["s"].to_numpy()), np.log10(d["contact"].to_numpy())
    d["slope"] = np.gradient(lc, ls)
    return d[["s", "contact", "slope"]]


# --------------------------------------------------------------------------
def saddle(source, eigenvector: pd.DataFrame, sample_col: str,
           resolution: int = 160_000, n_bins: int = 50,
           q_lo: float = 0.025, q_hi: float = 0.975):
    """Compartment saddle: observed/expected binned by eigenvector rank."""
    import cooltools
    clr = _clr(source, resolution)
    view = _view(clr)
    ev = eigenvector.rename(columns={sample_col: "E1"})[["chrom", "start", "end", "E1"]].copy()
    ev["chrom"] = ev["chrom"].astype(str)
    exp = cooltools.expected_cis(clr, view_df=view, nproc=1)
    q = np.linspace(q_lo, q_hi, n_bins)
    inter, sums, counts = cooltools.saddle(
        clr, exp, ev, "cis", n_bins=n_bins, qrange=(q_lo, q_hi), view_df=view)
    with np.errstate(invalid="ignore", divide="ignore"):
        S = sums / counts
    return S, (sums, counts)


def compartment_strength(S: np.ndarray, frac: float = 0.2) -> float:
    """(AA + BB) / 2AB from a saddle matrix -- higher = more segregated."""
    A = np.asarray(S, dtype=float)
    n = A.shape[0]
    k = max(1, int(n * frac))
    inner = A[1:-1, 1:-1] if n > 4 else A
    m = inner.shape[0]
    kk = max(1, int(m * frac))
    BB = np.nanmean(inner[:kk, :kk])
    AA = np.nanmean(inner[-kk:, -kk:])
    AB = np.nanmean([np.nanmean(inner[:kk, -kk:]), np.nanmean(inner[-kk:, :kk])])
    return float((AA + BB) / (2 * AB))


# --------------------------------------------------------------------------
def pileup(source, features: pd.DataFrame, resolution: int = 10_000,
           flank: int = 200_000, expected: bool = True) -> np.ndarray:
    """On-diagonal pileup (average map centered on each feature)."""
    clr = _clr(source, resolution)
    w = flank // resolution
    size = 2 * w + 1
    acc = np.zeros((size, size))
    n = 0
    feats = as_intervals(features)
    for chrom, grp in feats.groupby("chrom", sort=False):
        if chrom not in clr.chromnames:
            continue
        M = clr.matrix(balance=True).fetch(chrom)
        if expected:
            d = np.arange(M.shape[0])
            exp = np.array([np.nanmean(np.diagonal(M, k)) if M.shape[0] > abs(k) else np.nan
                            for k in range(M.shape[0])])
        mids = ((grp["start"] + grp["end"]) // 2 // resolution).to_numpy()
        for b in mids:
            lo, hi = b - w, b + w + 1
            if lo < 0 or hi > M.shape[0]:
                continue
            sub = M[lo:hi, lo:hi].copy()
            if expected:
                ii, jj = np.indices(sub.shape)
                e = exp[np.abs(ii - jj)]
                with np.errstate(invalid="ignore", divide="ignore"):
                    sub = sub / e
            if np.isfinite(sub).sum() < size:
                continue
            acc += np.nan_to_num(sub)
            n += 1
    if n == 0:
        raise ValueError("no usable features for pileup")
    LOG.info(f"pileup over {n:,} features")
    return acc / n


def apa(source, loops: pd.DataFrame, resolution: int = 10_000,
        flank: int = 100_000) -> tuple:
    """Aggregate peak analysis over off-diagonal loop anchors.

    ``loops`` needs chrom1/start1/end1/chrom2/start2/end2 (or chrom/start/end
    twice).  Returns (matrix, score) where score is center / corner mean.
    """
    clr = _clr(source, resolution)
    w = flank // resolution
    size = 2 * w + 1
    acc = np.zeros((size, size))
    n = 0
    lp = loops.copy()
    c1 = "chrom1" if "chrom1" in lp else "chrom"
    for chrom, grp in lp.groupby(c1, sort=False):
        if chrom not in clr.chromnames:
            continue
        M = clr.matrix(balance=True).fetch(chrom)
        a = ((grp["start1"] + grp["end1"]) // 2 // resolution).to_numpy() if "start1" in grp \
            else ((grp["start"] + grp["end"]) // 2 // resolution).to_numpy()
        b = ((grp["start2"] + grp["end2"]) // 2 // resolution).to_numpy()
        for i, j in zip(a, b):
            if i - w < 0 or j - w < 0 or i + w + 1 > M.shape[0] or j + w + 1 > M.shape[0]:
                continue
            sub = M[i - w:i + w + 1, j - w:j + w + 1]
            if np.isfinite(sub).sum() < size:
                continue
            acc += np.nan_to_num(sub)
            n += 1
    if n == 0:
        raise ValueError("no usable loops for APA")
    mat = acc / n
    c = size // 2
    k = max(1, size // 6)
    center = np.nanmean(mat[c - k:c + k + 1, c - k:c + k + 1])
    corner = np.nanmean(mat[-2 * k:, :2 * k])
    LOG.info(f"APA over {n:,} loops; score = {center / corner:.2f}")
    return mat, float(center / corner)


# --------------------------------------------------------------------------
def boundary_strength(insulation: pd.DataFrame, sample_col: str,
                      quantile: float = 0.05) -> pd.DataFrame:
    """Flag the strongest boundaries (lowest insulation score)."""
    d = insulation.dropna(subset=[sample_col]).copy()
    thr = d[sample_col].quantile(quantile)
    d["is_boundary"] = d[sample_col] < thr
    return d


def differential_insulation(ins_ref: pd.DataFrame, ins_test: pd.DataFrame,
                            ref_col: str, test_col: str,
                            anchors: pd.DataFrame | None = None,
                            resolution: int = 20_000) -> pd.DataFrame:
    """Insulation change per bin, optionally annotated by anchor overlap.

    Positive ``dIns`` = higher score = *weaker* boundary in the test condition.
    """
    a = ins_ref[["chrom", "start", "end", ref_col]].copy()
    b = ins_test[["chrom", "start", "end", test_col]].copy()
    for f in (a, b):
        f["chrom"] = f["chrom"].astype(str)
    m = a.merge(b, on=["chrom", "start", "end"]).dropna()
    m["dIns"] = m[test_col] - m[ref_col]
    if anchors is not None:
        an = as_intervals(anchors)
        key = an["chrom"].astype(str) + "_" + ((an["start"] + an["end"]) // 2 // resolution * resolution).astype(str)
        have = set(key)
        m["anchor"] = (m["chrom"] + "_" + (m["start"] // resolution * resolution).astype(str)).isin(have)
    return m
