"""Multi-omic module discovery over reference elements.

Elements are described by a stacked, per-assay z-scored matrix and partitioned
with k-means (fast, reproducible) or NMF (parts-based, non-negative).  Modules
are then characterised by their per-assay profile and tested for enrichment
among a set of elements of interest (e.g. those that changed in the KO).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy import stats as ss

from .utils import get_logger

LOG = get_logger()

__all__ = ["ModuleResult", "find_modules", "module_profile", "module_enrichment"]


@dataclass
class ModuleResult:
    labels: pd.Series
    profile: pd.DataFrame
    inertia: float | None = None
    method: str = "kmeans"

    @property
    def sizes(self) -> pd.Series:
        return self.labels.value_counts().sort_index()

    def __repr__(self):
        return (f"ModuleResult({self.method}, k={self.labels.nunique()})\n"
                + self.profile.round(2).to_string())


def _stack(layers: dict, scale: str = "zscore") -> tuple:
    """Stack per-assay element matrices into one feature matrix."""
    idx = None
    for v in layers.values():
        i = v.dropna().index
        idx = i if idx is None else idx.intersection(i)
    cols, names = [], []
    for name, v in layers.items():
        x = v.loc[idx].to_numpy(dtype=float)
        if scale == "zscore":
            sd = np.nanstd(x)
            x = (x - np.nanmean(x)) / (sd if sd > 0 else 1.0)
        elif scale == "rank":
            x = ss.rankdata(x) / len(x)
        cols.append(x.reshape(-1, 1))
        names.append(name)
    return np.hstack(cols), idx, names


def find_modules(layers: dict, k: int = 6, method: str = "kmeans",
                 scale: str = "zscore", seed: int = 0,
                 max_elements: int | None = None) -> ModuleResult:
    """Partition elements by their multi-omic profile.

    ``layers`` maps assay name -> element-level Series (e.g. WT signal).
    """
    X, idx, names = _stack(layers, scale=scale)
    keep = np.isfinite(X).all(axis=1)
    X, idx = X[keep], idx[keep]
    sub = np.arange(len(X))
    if max_elements and len(X) > max_elements:
        rng = np.random.default_rng(seed)
        sub = rng.choice(len(X), max_elements, replace=False)

    if method == "nmf":
        from sklearn.decomposition import NMF
        Xn = X - X.min(axis=0, keepdims=True)
        m = NMF(n_components=k, random_state=seed, init="nndsvda", max_iter=500)
        W = m.fit_transform(Xn)
        labels = W.argmax(axis=1)
        inertia = float(m.reconstruction_err_)
    else:
        from sklearn.cluster import KMeans
        km = KMeans(n_clusters=k, random_state=seed, n_init=10)
        km.fit(X[sub])
        labels = km.predict(X)
        inertia = float(km.inertia_)

    lab = pd.Series(labels, index=idx, name="module")
    prof = module_profile({n: pd.Series(X[:, j], index=idx) for j, n in enumerate(names)}, lab)
    LOG.info(f"modules: k={k} ({method}) over {len(idx):,} elements")
    return ModuleResult(labels=lab, profile=prof, inertia=inertia, method=method)


def module_profile(layers: dict, labels: pd.Series) -> pd.DataFrame:
    """Mean per-assay value in each module."""
    rows = {}
    for name, v in layers.items():
        idx = labels.index.intersection(v.dropna().index)
        rows[name] = v.loc[idx].groupby(labels.loc[idx]).mean()
    return pd.DataFrame(rows)


def module_enrichment(labels: pd.Series, selected: pd.Index,
                      background: pd.Index | None = None) -> pd.DataFrame:
    """Fisher enrichment of ``selected`` elements in each module."""
    bg = labels.index if background is None else labels.index.intersection(background)
    sel = labels.index.intersection(selected)
    rows = []
    for m in sorted(labels.unique()):
        in_m = labels.loc[bg] == m
        sel_m = labels.loc[sel] == m
        a = int(sel_m.sum())
        b = int(len(sel) - a)
        c = int(in_m.sum() - a)
        d = int(len(bg) - len(sel) - c)
        orr, p = ss.fisher_exact([[a, b], [c, d]])
        rows.append({"module": m, "n_selected": a, "n_module": int(in_m.sum()),
                     "odds_ratio": float(orr), "pvalue": float(p)})
    out = pd.DataFrame(rows)
    from .stats import bh_fdr
    out["padj"] = bh_fdr(out["pvalue"].to_numpy())
    return out.sort_values("pvalue")
