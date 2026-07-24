"""Publication-grade plotting.

House style: NPG-like palette, Helvetica when available, no top/right spines,
embedded fonts (Type 42) so PDFs are editable in Illustrator, rasterised scatter
layers so vector files stay small.
"""
from __future__ import annotations

import os

import matplotlib
import numpy as np
import pandas as pd
from scipy import stats as ss

matplotlib.use("Agg")
import matplotlib as mpl                    # noqa: E402
import matplotlib.pyplot as plt             # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402
from matplotlib.font_manager import FontProperties     # noqa: E402

__all__ = ["set_style", "PALETTE", "ma_plot", "volcano", "coupling_plot",
           "signal_heatmap", "audit_plot", "state_barplot", "pca_plot", "save"]

PALETTE = {
    "navy": "#3C5488", "teal": "#00A087", "red": "#E64B35", "orange": "#F39B7F",
    "blue": "#4DBBD5", "grey": "#B8B8B8", "dgrey": "#5A5A5A", "dark": "#1A2433",
    "purple": "#8491B4", "gold": "#E6A817",
}
# Helvetica is the usual journal requirement but is rarely installed; set
# EPIMUX_FONT to a .ttf/.ttc to use it, otherwise a sane sans-serif is used.
_FONT_CANDIDATES = [
    os.environ.get("EPIMUX_FONT", ""),
    "/usr/share/fonts/truetype/helvetica/Helvetica.ttf",
    os.path.expanduser("~/.fonts/Helvetica.ttf"),
    "/System/Library/Fonts/Helvetica.ttc",
    "/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf",
]
FP = FontProperties(family="DejaVu Sans")
DIVERGING = LinearSegmentedColormap.from_list("epimux", [PALETTE["navy"], "#f4f4f4", PALETTE["red"]])


def set_style(font: str | None = None):
    """Apply the house style; returns the FontProperties actually used."""
    global FP
    cands = ([font] if font else []) + _FONT_CANDIDATES
    for c in cands:
        if c and os.path.exists(c):
            try:
                fp = FontProperties(fname=c)
                fig = plt.figure(figsize=(0.4, 0.4))
                fig.text(0.1, 0.1, "Ag", fontproperties=fp)
                fig.canvas.draw()
                plt.close(fig)
                FP = fp
                try:
                    mpl.font_manager.fontManager.addfont(c)
                    mpl.rcParams["font.family"] = fp.get_name()
                except Exception:
                    pass
                break
            except Exception:
                continue
    mpl.rcParams.update({
        "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
        "axes.linewidth": 0.8, "xtick.major.width": 0.8, "ytick.major.width": 0.8,
        "xtick.major.size": 3.5, "ytick.major.size": 3.5,
        "axes.edgecolor": "#333333", "figure.dpi": 150, "savefig.dpi": 300,
    })
    return FP


set_style()


def _clean(ax, title=None, xlabel=None, ylabel=None, color=None, size=11):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.tick_params(labelsize=9)
    for lab in ax.get_xticklabels() + ax.get_yticklabels():
        lab.set_fontproperties(FP)
    if xlabel:
        ax.set_xlabel(xlabel, fontproperties=FP, fontsize=size)
    if ylabel:
        ax.set_ylabel(ylabel, fontproperties=FP, fontsize=size)
    if title:
        ax.set_title(title, fontproperties=FP, fontsize=size + 1,
                     color=color or PALETTE["dark"])
    return ax


def save(fig, path, formats=("png", "pdf")):
    base = os.path.splitext(str(path))[0]
    out = []
    for f in formats:
        p = f"{base}.{f}"
        fig.savefig(p, bbox_inches="tight", facecolor="white")
        out.append(p)
    plt.close(fig)
    return out


# --------------------------------------------------------------------------
def ma_plot(res: pd.DataFrame, fc=1.5, fdr=0.1, title=None, ax=None,
            up_label=None, down_label=None, ylim=(-4, 4)):
    """Effect size vs abundance, significance coloured by direction."""
    d = res.dropna(subset=["log2FC", "baseMean"])
    lfc = np.log2(fc) if res.attrs.get("value_kind") != "difference" else 0.1
    x = np.log10(d["baseMean"] + 1)
    y = d["log2FC"]
    sig = (d["padj"] < fdr) & (y.abs() > lfc)
    up, dn = sig & (y > 0), sig & (y < 0)
    if ax is None:
        fig, ax = plt.subplots(figsize=(4.6, 4.2))
    else:
        fig = ax.figure
    ax.scatter(x[~sig], y[~sig], s=2.5, c=PALETTE["grey"], alpha=.25, lw=0, rasterized=True)
    ax.scatter(x[dn], y[dn], s=7, c=PALETTE["navy"], alpha=.78, lw=0, rasterized=True)
    ax.scatter(x[up], y[up], s=7, c=PALETTE["red"], alpha=.78, lw=0, rasterized=True)
    ax.axhline(0, c="k", lw=.6)
    if ylim:
        ax.set_ylim(*ylim)
    _clean(ax, title, "log10 mean abundance", "log2 fold-change (test / ref)")
    ax.text(.97, .93, up_label or f"{int(up.sum())} up", transform=ax.transAxes,
            ha="right", fontproperties=FP, fontsize=10, color=PALETTE["red"], fontweight="bold")
    ax.text(.97, .05, down_label or f"{int(dn.sum())} down", transform=ax.transAxes,
            ha="right", fontproperties=FP, fontsize=10, color=PALETTE["navy"], fontweight="bold")
    return fig, ax


def volcano(res: pd.DataFrame, fc=1.5, fdr=0.1, title=None, ax=None, label_top=0):
    d = res.dropna(subset=["log2FC", "padj"])
    lfc = np.log2(fc)
    y = -np.log10(d["padj"].clip(lower=1e-300))
    sig = (d["padj"] < fdr) & (d["log2FC"].abs() > lfc)
    if ax is None:
        fig, ax = plt.subplots(figsize=(4.4, 4.4))
    else:
        fig = ax.figure
    ax.scatter(d["log2FC"][~sig], y[~sig], s=3, c=PALETTE["grey"], alpha=.3, lw=0, rasterized=True)
    up, dn = sig & (d["log2FC"] > 0), sig & (d["log2FC"] < 0)
    ax.scatter(d["log2FC"][up], y[up], s=8, c=PALETTE["red"], alpha=.8, lw=0, rasterized=True)
    ax.scatter(d["log2FC"][dn], y[dn], s=8, c=PALETTE["navy"], alpha=.8, lw=0, rasterized=True)
    ax.axvline(0, c="k", lw=.5)
    for v in (-lfc, lfc):
        ax.axvline(v, c=PALETTE["dgrey"], lw=.6, ls=":")
    ax.axhline(-np.log10(fdr), c=PALETTE["dgrey"], lw=.6, ls=":")
    _clean(ax, title, "log2 fold-change (test / ref)", "-log10 FDR")
    return fig, ax


def coupling_plot(res_x, res_y, name_x="X", name_y="Y", fc=1.5, fdr=0.1,
                  lim=3, ax=None, title=None):
    """Cross-assay effect-size scatter with the dual-significant set highlighted."""
    idx = res_x.index.intersection(res_y.index)
    x = res_x.loc[idx, "log2FC"]
    y = res_y.loc[idx, "log2FC"]
    ok = x.notna() & y.notna()
    x, y = x[ok], y[ok]
    rho = ss.spearmanr(x, y)[0]
    if ax is None:
        fig, ax = plt.subplots(figsize=(4.8, 4.4))
    else:
        fig = ax.figure
    ax.hexbin(x, y, gridsize=55, bins="log", cmap="Greys",
              extent=(-lim, lim, -lim, lim), mincnt=1)
    lf = np.log2(fc)
    dual = ((res_x.loc[x.index, "padj"] < fdr) & (x.abs() > lf) &
            (res_y.loc[y.index, "padj"] < fdr) & (y.abs() > lf))
    ax.scatter(x[dual], y[dual], s=16, c=PALETTE["red"], alpha=.9, lw=0)
    b1, b0 = np.polyfit(x, y, 1)
    xs = np.array([-lim, lim])
    ax.plot(xs, b1 * xs + b0, c=PALETTE["teal"], lw=2)
    ax.axhline(0, c="k", lw=.5)
    ax.axvline(0, c="k", lw=.5)
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    _clean(ax, title or f"rho = {rho:+.2f}",
           f"{name_x} log2FC", f"{name_y} log2FC", color=PALETTE["teal"])
    ax.text(.04, .93, f"n = {len(x):,}\ndual-sig = {int(dual.sum())}",
            transform=ax.transAxes, va="top", fontproperties=FP,
            fontsize=9, color=PALETTE["dgrey"])
    return fig, ax


def signal_heatmap(matrix: pd.DataFrame, groups: dict | None = None,
                   sort_by: pd.Series | None = None, vmax=1.5,
                   title=None, ylabel=None, ax=None):
    """Row z-scored replicate heatmap of selected elements."""
    M = matrix.to_numpy(dtype=float)
    z = (M - M.mean(1, keepdims=True)) / (M.std(1, keepdims=True) + 1e-9)
    if sort_by is not None:
        z = z[np.argsort(-np.asarray(sort_by))]
    if ax is None:
        fig, ax = plt.subplots(figsize=(4.2, 5.0))
    else:
        fig = ax.figure
    im = ax.imshow(z, aspect="auto", cmap=DIVERGING, vmin=-vmax, vmax=vmax,
                   interpolation="nearest")
    ax.set_xticks(range(matrix.shape[1]))
    ax.set_xticklabels(matrix.columns, fontproperties=FP, fontsize=8, rotation=45, ha="right")
    ax.set_yticks([])
    if groups:
        n = 0
        for g, s in list(groups.items())[:-1]:
            n += len([x for x in s if x in matrix.columns])
            ax.axvline(n - 0.5, c="k", lw=1.4)
    if ylabel:
        ax.set_ylabel(ylabel, fontproperties=FP, fontsize=9.5)
    if title:
        ax.set_title(title, fontproperties=FP, fontsize=12)
    cb = fig.colorbar(im, ax=ax, fraction=0.045, pad=0.03)
    cb.set_label("row z-score", fontproperties=FP, fontsize=8)
    cb.ax.tick_params(labelsize=7)
    return fig, ax


def state_barplot(states: pd.Series, ax=None, title=None):
    """Counts of multi-layer states (coordinated / discordant / ...)."""
    order = ["coordinated_up", "coordinated_down", "discordant", "single_layer"]
    vals = [int(states.get(k, 0)) for k in order]
    cols = [PALETTE["red"], PALETTE["navy"], PALETTE["gold"], PALETTE["grey"]]
    if ax is None:
        fig, ax = plt.subplots(figsize=(4.6, 4.2))
    else:
        fig = ax.figure
    ax.bar(range(len(order)), vals, color=cols, width=.66, edgecolor="white")
    for i, v in enumerate(vals):
        ax.text(i, v + max(vals) * 0.02, f"{v:,}", ha="center",
                fontproperties=FP, fontsize=11, fontweight="bold", color=cols[i])
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels([o.replace("_", "\n") for o in order], fontproperties=FP, fontsize=9)
    ax.set_ylim(0, max(max(vals), 1) * 1.2)
    _clean(ax, title or "multi-layer response", None, "elements")
    return fig, ax


def pca_plot(matrix: pd.DataFrame, groups: dict, n_top=2000, ax=None, title=None,
             markers=None, colors=None):
    X = np.log2(matrix.to_numpy(dtype=float) /
                (matrix.to_numpy(dtype=float).sum(0, keepdims=True) / 1e6) + 1)
    X = X[np.argsort(X.var(1))[-n_top:]]
    X = X - X.mean(1, keepdims=True)
    U, S, _ = np.linalg.svd(X.T - X.T.mean(0), full_matrices=False)
    pcs = U[:, :2] * S[:2]
    ve = (S ** 2 / (S ** 2).sum())[:2] * 100
    g_of = {s: g for g, ss_ in groups.items() for s in ss_}
    palette = colors or {g: c for g, c in zip(groups, [PALETTE["navy"], PALETTE["teal"],
                                                       PALETTE["red"], PALETTE["gold"]])}
    if ax is None:
        fig, ax = plt.subplots(figsize=(4.4, 4.0))
    else:
        fig = ax.figure
    for s, (px, py) in zip(matrix.columns, pcs):
        g = g_of.get(s, "other")
        ax.scatter(px, py, c=palette.get(g, PALETTE["grey"]),
                   marker=(markers or {}).get(g, "o"), s=90,
                   edgecolor="white", lw=1, zorder=3)
    _clean(ax, title or "sample PCA", f"PC1 ({ve[0]:.0f}%)", f"PC2 ({ve[1]:.0f}%)")
    from matplotlib.lines import Line2D
    ax.legend(handles=[Line2D([0], [0], marker="o", color="w", markersize=9,
                              markerfacecolor=palette.get(g, PALETTE["grey"]), label=g)
                       for g in groups],
              prop=FP, fontsize=8, frameon=False, loc="best")
    return fig, ax


def audit_plot(report, ax=None):
    """Traffic-light summary of the audit battery."""
    df = report.to_frame()
    if df.empty:
        raise ValueError("empty audit report")
    cmap = {"pass": PALETTE["teal"], "warn": PALETTE["gold"], "fail": PALETTE["red"]}
    if ax is None:
        fig, ax = plt.subplots(figsize=(7.2, 0.42 * len(df) + 1.0))
    else:
        fig = ax.figure
    y = np.arange(len(df))[::-1]
    ax.barh(y, 1, color=[cmap[s] for s in df["status"]], height=.7)
    for yi, (_, r) in zip(y, df.iterrows()):
        ax.text(1.05, yi, r["summary"][:96], va="center",
                fontproperties=FP, fontsize=8.5, color=PALETTE["dark"])
    ax.set_yticks(y)
    ax.set_yticklabels(df["check"], fontproperties=FP, fontsize=9)
    ax.set_xlim(0, 6)
    ax.set_xticks([])
    for s in ("top", "right", "bottom"):
        ax.spines[s].set_visible(False)
    ax.set_title("epimux audit", fontproperties=FP, fontsize=12, loc="left")
    return fig, ax
