"""Enrichment: pathways for linked genes, motifs for elements, overlap for sets.

The recurring mistake this module guards against is a background set that does
not match the foreground.  Testing "genes near changed enhancers" against *all
genes* conflates enhancer-density with biology; the correct background is genes
near *tested* enhancers.  Every function here takes an explicit background and
refuses to invent one.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats as ss

from .stats import bh_fdr
from .utils import as_intervals, get_logger, overlap

LOG = get_logger()

__all__ = ["pathway_enrichment", "motif_enrichment", "overlap_enrichment",
           "gc_matched_background"]


# --------------------------------------------------------------------------
def pathway_enrichment(genes, background=None, gene_sets="MSigDB_Hallmark_2020",
                       organism: str = "mouse", top: int = 15) -> pd.DataFrame:
    """Pathway enrichment via gseapy/Enrichr.

    ``background`` should be the genes linked to *tested* elements, not the
    whole genome.  Enrichr ignores custom backgrounds for its hosted libraries,
    so when one is supplied the function additionally reports a Fisher test
    against it -- use that column, not the raw Enrichr P.
    """
    genes = [g for g in pd.unique(pd.Series(list(genes)).dropna().astype(str)) if g]
    if len(genes) < 5:
        LOG.warning(f"only {len(genes)} genes -- enrichment will be unstable")
    try:
        import gseapy as gp
    except ImportError as e:  # pragma: no cover
        raise ImportError("gseapy required for pathway enrichment") from e
    er = gp.enrichr(gene_list=genes, gene_sets=[gene_sets] if isinstance(gene_sets, str)
                    else list(gene_sets), organism=organism, outdir=None)
    res = er.results.copy()
    res = res.sort_values("Adjusted P-value").head(top)
    res["n_query"] = len(genes)
    if background is not None:
        bg = set(pd.Series(list(background)).dropna().astype(str))
        q = set(genes)
        rows = []
        for _, r in res.iterrows():
            members = set(str(r["Genes"]).split(";"))
            a = len(members & q)
            b = len(q) - a
            c = len(members & bg) - a
            d = len(bg) - len(q) - c
            orr, p = ss.fisher_exact([[a, b], [max(c, 0), max(d, 0)]])
            rows.append((orr, p))
        res["OR_vs_background"] = [r[0] for r in rows]
        res["P_vs_background"] = [r[1] for r in rows]
        res["padj_vs_background"] = bh_fdr(res["P_vs_background"].to_numpy())
        LOG.info("reported P_vs_background uses your background; prefer it over the Enrichr P")
    return res


# --------------------------------------------------------------------------
def gc_matched_background(elements: pd.DataFrame, selected: pd.Index,
                          gc: pd.Series, n_per: int = 5, tol: float = 0.02,
                          seed: int = 0) -> pd.Index:
    """Sample a background matched on GC content (and thus roughly on width).

    Motif enrichment without GC matching mostly rediscovers GC-rich motifs.
    """
    rng = np.random.default_rng(seed)
    gc = pd.Series(gc, index=elements.index).dropna()
    sel = gc.index.intersection(selected)
    pool = gc.index.difference(sel)
    out = []
    pool_gc = gc.loc[pool]
    for i in sel:
        target = gc.loc[i]
        cand = pool_gc.index[(pool_gc - target).abs() <= tol]
        if len(cand) == 0:
            continue
        out.extend(rng.choice(cand, size=min(n_per, len(cand)), replace=False))
    return pd.Index(pd.unique(pd.Series(out)))


def motif_enrichment(counts_fg: pd.DataFrame, counts_bg: pd.DataFrame) -> pd.DataFrame:
    """Fisher enrichment of motif hits, foreground vs an explicit background.

    Inputs are boolean/count matrices (elements x motifs); a motif is "present"
    when its value is > 0.  Supply a GC-matched background
    (:func:`gc_matched_background`) unless you know you do not need one.
    """
    fg = (counts_fg > 0).sum(axis=0)
    bg = (counts_bg > 0).sum(axis=0)
    nfg, nbg = len(counts_fg), len(counts_bg)
    rows = []
    for m in fg.index:
        a, c = int(fg[m]), int(bg.get(m, 0))
        orr, p = ss.fisher_exact([[a, nfg - a], [c, nbg - c]])
        rows.append({"motif": m, "n_fg": a, "n_bg": c,
                     "frac_fg": a / nfg, "frac_bg": c / max(nbg, 1),
                     "odds_ratio": float(orr), "pvalue": float(p)})
    out = pd.DataFrame(rows)
    out["padj"] = bh_fdr(out["pvalue"].to_numpy())
    return out.sort_values("pvalue")


# --------------------------------------------------------------------------
def overlap_enrichment(query: pd.DataFrame, target: pd.DataFrame,
                       universe: pd.DataFrame, n_shuffle: int = 1000,
                       seed: int = 0) -> dict:
    """Is ``query`` enriched for overlap with ``target``, relative to ``universe``?

    Significance by shuffling query labels within the universe -- which keeps
    the chromosome and width distribution of the query, unlike shuffling
    coordinates across the genome.
    """
    rng = np.random.default_rng(seed)
    uni = as_intervals(universe).reset_index(drop=True)
    q = as_intervals(query).reset_index(drop=True)
    t = as_intervals(target)
    hit_u = np.zeros(len(uni), dtype=bool)
    pr = overlap(uni, t)
    if not pr.empty:
        hit_u[pr["idx_a"].astype(int).to_numpy()] = True
    pr_q = overlap(q, t)
    obs = int(pd.unique(pr_q["idx_a"]).size) if not pr_q.empty else 0
    k = len(q)
    null = np.array([hit_u[rng.choice(len(uni), size=k, replace=False)].sum()
                     for _ in range(n_shuffle)])
    p = float((np.sum(null >= obs) + 1) / (n_shuffle + 1))
    exp = float(null.mean())
    return {"observed": obs, "expected": exp, "n_query": k,
            "fold_enrichment": obs / exp if exp else np.inf,
            "pvalue": p, "null_sd": float(null.std())}
