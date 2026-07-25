"""Element annotation: genomic context, TSS distance, super-enhancers.

Keeps the interpretation of an element list honest -- "distal enhancers" that
turn out to be promoters, or a super-enhancer set that is really a peak-density
artifact, change conclusions more often than people expect.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from .utils import as_intervals, get_logger, overlap, sort_intervals

LOG = get_logger()

__all__ = ["annotate_tss", "classify_context", "stitch_super_enhancers",
           "context_composition", "distance_enrichment"]


# --------------------------------------------------------------------------
def annotate_tss(elements: pd.DataFrame, tss: pd.DataFrame) -> pd.DataFrame:
    """Signed distance to the nearest TSS (negative = upstream of the gene)."""
    el = as_intervals(elements).reset_index(drop=True)
    ts = as_intervals(tss).reset_index(drop=True)
    if "gene" not in ts.columns:
        ts = ts.assign(gene=[f"g{i}" for i in range(len(ts))])
    strand = ts["strand"].to_numpy() if "strand" in ts.columns else np.array(["+"] * len(ts))
    out = pd.DataFrame({"nearest_gene": pd.NA, "tss_distance": np.nan}, index=el.index)
    by = {c: g.sort_values("start") for c, g in ts.groupby("chrom", sort=False)}
    for c, grp in el.groupby("chrom", sort=False):
        g = by.get(c)
        if g is None or len(g) == 0:
            continue
        pos = g["start"].to_numpy()
        names = g["gene"].to_numpy()
        st = strand[g.index.to_numpy()]
        mids = ((grp["start"] + grp["end"]) // 2).to_numpy()
        j = np.searchsorted(pos, mids)
        for k, (idx, mid) in enumerate(zip(grp.index.to_numpy(), mids)):
            cand = [q for q in (j[k] - 1, j[k]) if 0 <= q < len(pos)]
            if not cand:
                continue
            q = min(cand, key=lambda q: abs(pos[q] - mid))
            d = int(mid - pos[q])
            if st[q] == "-":
                d = -d
            out.loc[idx, "nearest_gene"] = names[q]
            out.loc[idx, "tss_distance"] = d
    return out


def classify_context(elements: pd.DataFrame, tss: pd.DataFrame,
                     promoter: int = 2_000, proximal: int = 10_000,
                     features: dict | None = None) -> pd.DataFrame:
    """Label elements promoter / proximal / distal, plus optional feature overlap.

    ``features`` maps a label -> interval frame (exons, CpG islands, repeats...).
    """
    ann = annotate_tss(elements, tss)
    d = ann["tss_distance"].abs()
    ann["context"] = np.select(
        [d <= promoter, d <= proximal],
        ["promoter", "proximal"], default="distal")
    ann.loc[d.isna(), "context"] = "unassigned"
    if features:
        el = as_intervals(elements).reset_index(drop=True)
        for name, feat in features.items():
            pairs = overlap(el, as_intervals(feat))
            hit = np.zeros(len(el), dtype=bool)
            if not pairs.empty:
                hit[pairs["idx_a"].astype(int).to_numpy()] = True
            ann[f"in_{name}"] = hit
    return ann


def context_composition(ann: pd.DataFrame, selected: pd.Index | None = None) -> pd.DataFrame:
    """Composition of a selected set vs background, with enrichment."""
    bg = ann["context"].value_counts(normalize=True)
    if selected is None:
        return bg.to_frame("background")
    sel = ann.loc[ann.index.intersection(selected), "context"].value_counts(normalize=True)
    out = pd.DataFrame({"background": bg, "selected": sel}).fillna(0.0)
    out["enrichment"] = out["selected"] / out["background"].replace(0, np.nan)
    return out


# --------------------------------------------------------------------------
def stitch_super_enhancers(peaks: pd.DataFrame, signal: pd.Series,
                           stitch_distance: int = 12_500,
                           tss: pd.DataFrame | None = None,
                           exclude_promoter: int = 2_000) -> pd.DataFrame:
    """ROSE-style super-enhancer calling.

    Peaks within ``stitch_distance`` are merged, signal is summed, stitched
    regions are ranked, and the cut-off is the point of slope 1 on the ranked
    signal curve (the standard ROSE geometric definition).
    """
    pk = as_intervals(peaks).reset_index(drop=True)
    pk["signal"] = np.asarray(signal, dtype=float)
    if tss is not None:
        ann = annotate_tss(pk, tss)
        keep = ~(ann["tss_distance"].abs() <= exclude_promoter)
        LOG.info(f"super-enhancers: dropped {int((~keep).sum()):,} promoter-proximal peaks")
        pk = pk[keep.to_numpy()].reset_index(drop=True)
    pk = sort_intervals(pk)

    rows, cur = [], None
    for r in pk.itertuples():
        if cur and r.chrom == cur["chrom"] and r.start - cur["end"] <= stitch_distance:
            cur["end"] = max(cur["end"], r.end)
            cur["signal"] += r.signal
            cur["n_peaks"] += 1
        else:
            if cur:
                rows.append(cur)
            cur = {"chrom": r.chrom, "start": r.start, "end": r.end,
                   "signal": float(r.signal), "n_peaks": 1}
    if cur:
        rows.append(cur)
    st = pd.DataFrame(rows).sort_values("signal", ascending=False).reset_index(drop=True)

    y = st["signal"].to_numpy()[::-1]                 # ascending
    x = np.arange(len(y))
    if len(y) < 3:
        st["super_enhancer"] = False
        return st
    xs = x / x.max()
    ys = (y - y.min()) / (y.max() - y.min() + 1e-12)
    slope = np.gradient(ys, xs)
    cut = int(np.argmax(slope >= 1.0)) if (slope >= 1.0).any() else int(0.95 * len(y))
    n_se = len(y) - cut
    st["super_enhancer"] = np.arange(len(st)) < n_se
    st["rank"] = np.arange(1, len(st) + 1)
    LOG.info(f"super-enhancers: {n_se:,} of {len(st):,} stitched regions")
    return st


# --------------------------------------------------------------------------
def distance_enrichment(ann: pd.DataFrame, selected: pd.Index,
                        bins=(0, 2_000, 10_000, 50_000, 200_000, 1_000_000, np.inf)) -> pd.DataFrame:
    """Are the selected elements distributed differently with TSS distance?"""
    from scipy import stats as ss
    d = ann["tss_distance"].abs()
    lab = pd.cut(d, bins=list(bins), right=False)
    bg = lab.value_counts().sort_index()
    sel = lab.loc[ann.index.intersection(selected)].value_counts().sort_index()
    rows = []
    for b in bg.index:
        a = int(sel.get(b, 0))
        rest_sel = int(sel.sum() - a)
        c = int(bg.get(b, 0)) - a
        rest_bg = int(bg.sum() - sel.sum() - c)
        orr, p = ss.fisher_exact([[a, rest_sel], [max(c, 0), max(rest_bg, 0)]])
        rows.append({"bin": str(b), "n_selected": a, "n_background": int(bg.get(b, 0)),
                     "odds_ratio": float(orr), "pvalue": float(p)})
    out = pd.DataFrame(rows)
    from .stats import bh_fdr
    out["padj"] = bh_fdr(out["pvalue"].to_numpy())
    return out
