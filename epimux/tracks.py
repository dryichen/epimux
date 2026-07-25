"""Genome-browser-style figure panels.

A locus figure is where multi-omic claims are usually checked by eye, so it is
worth drawing properly: shared x-axis, per-track group means with replicate
spread, and an explicit statement of what is being averaged.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from .plotting import FP, PALETTE, _clean, save  # noqa: F401
from .utils import as_intervals, get_logger

LOG = get_logger()

__all__ = ["extract_profile", "locus_plot", "metaprofile", "heatmap_profile"]


def extract_profile(bigwig, chrom: str, start: int, end: int, nbins: int = 500):
    """Binned mean signal over a locus."""
    import pyBigWig
    bw = bigwig if hasattr(bigwig, "stats") else pyBigWig.open(str(bigwig))
    try:
        v = bw.stats(chrom, int(start), int(end), type="mean", nBins=nbins)
    except RuntimeError:
        v = [None] * nbins
    finally:
        if not hasattr(bigwig, "stats"):
            bw.close()
    return np.array([np.nan if x is None else x for x in v], dtype=float)


def locus_plot(tracks: dict, chrom: str, start: int, end: int,
               groups: dict | None = None, features: dict | None = None,
               nbins: int = 600, colors: dict | None = None,
               title: str | None = None, height_per_track: float = 0.9,
               show_spread: bool = True):
    """Stacked browser view.

    ``tracks``  : {track name: bigWig path}  or  {track name: {sample: path}}
    ``groups``  : {group: [sample names]} -- when a track maps samples, the mean
                  of each group is drawn with the replicate range shaded.
    ``features``: {label: interval frame} drawn as blocks under the tracks.
    """
    import matplotlib.pyplot as plt
    n = len(tracks) + (1 if features else 0)
    fig, axes = plt.subplots(n, 1, figsize=(9.0, height_per_track * n + 1.0),
                             sharex=True, gridspec_kw={"hspace": 0.25})
    if n == 1:
        axes = [axes]
    x = np.linspace(start, end, nbins) / 1e6
    palette = colors or {"WT": PALETTE["navy"], "KO": PALETTE["red"]}

    for ax, (name, src) in zip(axes, tracks.items()):
        if isinstance(src, dict) and groups:
            for g, samples in groups.items():
                mats = [extract_profile(src[s], chrom, start, end, nbins)
                        for s in samples if s in src]
                if not mats:
                    continue
                M = np.vstack(mats)
                mu = np.nanmean(M, 0)
                c = palette.get(g, PALETTE["grey"])
                ax.plot(x, mu, color=c, lw=1.2, label=g)
                if show_spread and M.shape[0] > 1:
                    ax.fill_between(x, np.nanmin(M, 0), np.nanmax(M, 0), color=c, alpha=.20, lw=0)
        else:
            path = src if not isinstance(src, dict) else list(src.values())[0]
            v = extract_profile(path, chrom, start, end, nbins)
            ax.fill_between(x, 0, v, color=PALETTE["navy"], alpha=.85, lw=0)
        ax.set_ylabel(name, fontproperties=FP, fontsize=9, rotation=0,
                      ha="right", va="center", labelpad=8)
        ax.set_yticks([])
        for s in ("top", "right", "left"):
            ax.spines[s].set_visible(False)
        ax.margins(x=0)

    if features:
        ax = axes[-1]
        for k, (label, feat) in enumerate(features.items()):
            f = as_intervals(feat)
            f = f[(f["chrom"] == chrom) & (f["end"] > start) & (f["start"] < end)]
            for _, r in f.iterrows():
                ax.add_patch(plt.Rectangle((r["start"] / 1e6, -k - 0.35),
                                           (r["end"] - r["start"]) / 1e6, 0.7,
                                           color=PALETTE["teal"]))
            ax.text(start / 1e6, -k, f" {label}", va="center", ha="left",
                    fontproperties=FP, fontsize=8, color=PALETTE["dgrey"])
        ax.set_ylim(-len(features), 1)
        ax.set_yticks([])
        for s in ("top", "right", "left"):
            ax.spines[s].set_visible(False)

    axes[-1].set_xlabel(f"{chrom} (Mb)", fontproperties=FP, fontsize=10)
    for lab in axes[-1].get_xticklabels():
        lab.set_fontproperties(FP)
    if groups:
        axes[0].legend(prop=FP, fontsize=8, frameon=False, ncol=len(groups), loc="upper right")
    if title:
        fig.suptitle(title, fontproperties=FP, fontsize=12, y=1.0)
    return fig, axes


# --------------------------------------------------------------------------
def metaprofile(bigwigs: dict, regions: pd.DataFrame, flank: int = 2_000,
                nbins: int = 100, groups: dict | None = None,
                colors: dict | None = None, ax=None, title: str | None = None,
                center: str = "midpoint"):
    """Average signal centred on a set of regions."""
    import matplotlib.pyplot as plt
    reg = as_intervals(regions)
    if center == "midpoint":
        mids = ((reg["start"] + reg["end"]) // 2).to_numpy()
    else:
        mids = reg["start"].to_numpy()
    prof = {}
    for sample, path in bigwigs.items():
        import pyBigWig
        bw = pyBigWig.open(str(path))
        acc = np.zeros(nbins)
        n = 0
        for c, m in zip(reg["chrom"], mids):
            if c not in bw.chroms():
                continue
            a, b = int(m - flank), int(m + flank)
            if a < 0 or b > bw.chroms()[c]:
                continue
            try:
                v = bw.stats(c, a, b, type="mean", nBins=nbins)
            except RuntimeError:
                continue
            v = np.array([np.nan if q is None else q for q in v], dtype=float)
            if np.isfinite(v).sum() < nbins // 2:
                continue
            acc += np.nan_to_num(v)
            n += 1
        bw.close()
        prof[sample] = acc / max(n, 1)
        LOG.info(f"metaprofile {sample}: {n:,} regions")
    if ax is None:
        fig, ax = plt.subplots(figsize=(4.6, 3.8))
    else:
        fig = ax.figure
    x = np.linspace(-flank, flank, nbins) / 1000
    palette = colors or {"WT": PALETTE["navy"], "KO": PALETTE["red"]}
    if groups:
        for g, samples in groups.items():
            M = np.vstack([prof[s] for s in samples if s in prof])
            mu, sd = M.mean(0), M.std(0)
            c = palette.get(g, PALETTE["grey"])
            ax.plot(x, mu, color=c, lw=1.8, label=g)
            ax.fill_between(x, mu - sd, mu + sd, color=c, alpha=.20, lw=0)
        ax.legend(prop=FP, fontsize=9, frameon=False)
    else:
        for s, v in prof.items():
            ax.plot(x, v, lw=1.4, label=s)
        ax.legend(prop=FP, fontsize=8, frameon=False)
    ax.axvline(0, c="k", lw=.6, ls=":")
    _clean(ax, title or f"n = {len(reg):,} regions", "distance from centre (kb)", "mean signal")
    return fig, ax


def heatmap_profile(bigwig, regions: pd.DataFrame, flank: int = 2_000,
                    nbins: int = 100, sort_by: pd.Series | None = None,
                    vmax: float | None = None, ax=None, title: str | None = None):
    """Per-region signal heatmap (deepTools-style), rows optionally sorted."""
    import matplotlib.pyplot as plt
    import pyBigWig
    from .plotting import DIVERGING  # noqa: F401
    reg = as_intervals(regions)
    mids = ((reg["start"] + reg["end"]) // 2).to_numpy()
    bw = pyBigWig.open(str(bigwig))
    rows = []
    for c, m in zip(reg["chrom"], mids):
        if c not in bw.chroms():
            rows.append(np.full(nbins, np.nan)); continue
        a, b = int(m - flank), int(m + flank)
        if a < 0 or b > bw.chroms()[c]:
            rows.append(np.full(nbins, np.nan)); continue
        try:
            v = bw.stats(c, a, b, type="mean", nBins=nbins)
        except RuntimeError:
            v = [None] * nbins
        rows.append(np.array([np.nan if q is None else q for q in v], dtype=float))
    bw.close()
    M = np.vstack(rows)
    if sort_by is not None:
        M = M[np.argsort(-np.asarray(sort_by))]
    else:
        M = M[np.argsort(-np.nanmean(M, 1))]
    if ax is None:
        fig, ax = plt.subplots(figsize=(3.4, 5.0))
    else:
        fig = ax.figure
    im = ax.imshow(M, aspect="auto", cmap="Reds",
                   vmax=vmax or np.nanpercentile(M, 98), interpolation="nearest")
    ax.set_xticks([0, nbins // 2, nbins - 1])
    ax.set_xticklabels([f"-{flank//1000}kb", "0", f"+{flank//1000}kb"],
                       fontproperties=FP, fontsize=8)
    ax.set_yticks([])
    ax.set_ylabel(f"{len(M):,} regions", fontproperties=FP, fontsize=9)
    if title:
        ax.set_title(title, fontproperties=FP, fontsize=11)
    fig.colorbar(im, ax=ax, fraction=0.045, pad=0.03)
    return fig, ax
