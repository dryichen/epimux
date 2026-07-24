"""Element -> gene linking.

Nearest-gene assignment is the default in most pipelines and it is wrong often
enough to change conclusions: enhancers routinely skip their neighbour.  When a
contact map is available, epimux links by *contact*, using an Activity-by-Contact
style score (Fulco et al. 2019):

    score(e, g) = activity(e) * contact(e, g) / sum_over_e' in window

``activity`` defaults to the geometric mean of accessibility and activity marks,
which is what ABC uses; any element-level vector can be substituted.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from .utils import as_intervals, get_logger, map_to_reference, overlap

LOG = get_logger()

__all__ = ["nearest_gene", "abc_link", "aggregate_to_genes"]


def nearest_gene(elements: pd.DataFrame, tss: pd.DataFrame,
                 max_distance: int = 1_000_000) -> pd.DataFrame:
    """Baseline nearest-TSS assignment (kept for comparison, not recommended)."""
    el = as_intervals(elements).reset_index(drop=True)
    ts = as_intervals(tss).reset_index(drop=True)
    if "gene" not in ts.columns:
        raise ValueError("tss frame needs a 'gene' column")
    rows = []
    ts_by = {c: g.sort_values("start") for c, g in ts.groupby("chrom", sort=False)}
    for i, (c, s, e) in enumerate(zip(el["chrom"], el["start"], el["end"])):
        g = ts_by.get(c)
        if g is None:
            continue
        mid = (s + e) // 2
        pos = g["start"].to_numpy()
        j = np.searchsorted(pos, mid)
        cand = [k for k in (j - 1, j) if 0 <= k < len(pos)]
        if not cand:
            continue
        k = min(cand, key=lambda k: abs(pos[k] - mid))
        d = int(abs(pos[k] - mid))
        if d <= max_distance:
            rows.append({"element": i, "gene": g["gene"].to_numpy()[k], "distance": d})
    return pd.DataFrame(rows)


def abc_link(elements: pd.DataFrame, tss: pd.DataFrame, hic, sample: str,
             activity: pd.Series | None = None, resolution: int = 20_000,
             window: int = 1_000_000, min_score: float = 0.02,
             exclude_promoter: int = 2_000) -> pd.DataFrame:
    """Activity-by-Contact linking using a real contact map.

    Parameters
    ----------
    hic
        an :class:`~epimux.assays.HiCAssay` (or anything exposing ``.coolers``).
    activity
        element-level activity; if ``None`` a uniform value is used, which
        reduces the score to pure contact frequency.

    Returns a frame ``element, gene, contact, activity, score``.
    """
    import cooler

    el = as_intervals(elements).reset_index(drop=True)
    ts = as_intervals(tss).reset_index(drop=True)
    if activity is None:
        activity = pd.Series(1.0, index=np.arange(len(el)))
    activity = pd.Series(np.asarray(activity, dtype=float), index=np.arange(len(el)))

    clr = cooler.Cooler(f"{hic.coolers[sample]}::/resolutions/{resolution}")
    rows = []
    ts_by = {c: g for c, g in ts.groupby("chrom", sort=False)}

    for chrom, grp in el.groupby("chrom", sort=False):
        if chrom not in clr.chromnames:
            continue
        genes = ts_by.get(chrom)
        if genes is None or len(genes) == 0:
            continue
        M = clr.matrix(balance=True).fetch(chrom)
        n = M.shape[0]
        g_bin = ((genes["start"]) // resolution).to_numpy()
        g_name = genes["gene"].to_numpy()
        e_mid = ((grp["start"] + grp["end"]) // 2).to_numpy()
        e_bin = e_mid // resolution
        w = window // resolution

        for pos, emid, eb in zip(grp.index.to_numpy(), e_mid, e_bin):
            if eb < 0 or eb >= n:
                continue
            sel = np.abs(g_bin - eb) <= w
            if not sel.any():
                continue
            gb = np.clip(g_bin[sel], 0, n - 1)
            contact = M[eb, gb]
            contact = np.nan_to_num(contact, nan=0.0)
            # drop the element's own promoter-proximal gene bin if requested
            dist = np.abs(genes["start"].to_numpy()[sel] - emid)
            contact = np.where(dist < exclude_promoter, 0.0, contact)
            act = float(activity.get(pos, 0.0))
            num = act * contact
            denom = num.sum()
            if denom <= 0:
                continue
            score = num / denom
            keep = score >= min_score
            for gn, ct, sc in zip(g_name[sel][keep], contact[keep], score[keep]):
                rows.append({"element": int(pos), "gene": gn,
                             "contact": float(ct), "activity": act, "score": float(sc)})
    out = pd.DataFrame(rows)
    LOG.info(f"ABC: {len(out):,} element-gene links over {out['element'].nunique() if len(out) else 0:,} elements")
    return out


def aggregate_to_genes(links: pd.DataFrame, element_values: pd.Series,
                       weight: str = "score", how: str = "weighted_mean") -> pd.Series:
    """Push element-level effect sizes onto genes through the link table."""
    if links.empty:
        return pd.Series(dtype=float)
    df = links.copy()
    df["val"] = np.asarray(element_values)[df["element"].to_numpy()]
    df = df.dropna(subset=["val"])
    if how == "sum":
        return df.groupby("gene")["val"].sum()
    if how == "max_abs":
        df["a"] = df["val"].abs()
        return df.sort_values("a", ascending=False).drop_duplicates("gene").set_index("gene")["val"]
    g = df.groupby("gene")
    return g.apply(lambda d: np.average(d["val"], weights=np.maximum(d[weight], 1e-9)))
