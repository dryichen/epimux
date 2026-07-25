"""Comparing effects across contrasts.

The question "is the KO effect the same in LSK and in GMP?" is not answered by
comparing two significance lists.  Two lists differ for reasons that have
nothing to do with biology -- different replicate counts, different depth,
different dispersion -- so an element significant in one and not the other is
weak evidence of a real difference.

This module compares **effect sizes with their uncertainty**:

* :func:`compare_contrasts` tests the *interaction* per element (is the effect
  in A different from the effect in B?), which is the question people mean.
* :func:`concordance_summary` reports how much of the apparent difference is
  explained by power rather than biology.
* :func:`meta_analyse` combines contrasts by inverse-variance weighting when the
  intent is to pool rather than contrast.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats as ss

from .stats import bh_fdr
from .utils import get_logger

LOG = get_logger()

__all__ = ["compare_contrasts", "concordance_summary", "meta_analyse",
           "replication_rate"]


def _se(res: pd.DataFrame) -> pd.Series:
    """Standard error of log2FC, from the Wald statistic when not stored."""
    if "lfcSE" in res.columns:
        return res["lfcSE"]
    with np.errstate(invalid="ignore", divide="ignore"):
        se = (res["log2FC"] / res["stat"]).abs()
    return se.replace([np.inf, -np.inf], np.nan)


def compare_contrasts(res_a: pd.DataFrame, res_b: pd.DataFrame,
                      name_a: str = "A", name_b: str = "B") -> pd.DataFrame:
    """Per-element interaction test: does the effect differ between contrasts?

    Uses a two-sample Wald test on the difference of log2 fold-changes,
    ``(b - a) / sqrt(se_a^2 + se_b^2)``.  Requires both results to share an
    index and orientation.
    """
    for r, n in ((res_a, name_a), (res_b, name_b)):
        if "log2FC" not in r:
            raise ValueError(f"{n}: not a differential result")
    ca, cb = res_a.attrs.get("contrast"), res_b.attrs.get("contrast")
    if ca and cb:
        oa = ca.split("log2(")[1].split(")")[0] if "log2(" in ca else None
        ob = cb.split("log2(")[1].split(")")[0] if "log2(" in cb else None
        if oa and ob and oa.split("/")[0] != ob.split("/")[0]:
            LOG.warning(f"contrast orientations differ ({oa} vs {ob}); "
                        "the interaction sign will follow these definitions")

    idx = res_a.index.intersection(res_b.index)
    a, b = res_a.loc[idx], res_b.loc[idx]
    sea, seb = _se(a), _se(b)
    diff = b["log2FC"] - a["log2FC"]
    se = np.sqrt(sea ** 2 + seb ** 2)
    with np.errstate(invalid="ignore", divide="ignore"):
        z = diff / se
    p = 2 * ss.norm.sf(np.abs(z))
    out = pd.DataFrame({
        f"log2FC_{name_a}": a["log2FC"], f"log2FC_{name_b}": b["log2FC"],
        f"padj_{name_a}": a["padj"], f"padj_{name_b}": b["padj"],
        "difference": diff, "se_difference": se, "z": z,
        "pvalue_interaction": p,
    }, index=idx)
    out["padj_interaction"] = bh_fdr(out["pvalue_interaction"].to_numpy())
    out.attrs["comparison"] = f"{name_b} - {name_a}"
    return out


def concordance_summary(cmp: pd.DataFrame, name_a: str = "A", name_b: str = "B",
                        fc: float = 1.5, fdr: float = 0.1) -> dict:
    """How much of an apparent difference between contrasts is real?

    Splits elements into shared / A-only / B-only by significance, then asks how
    many of the "only" calls are actually supported by an interaction test.  A
    large gap between the two numbers means the difference is mostly power.
    """
    lf = np.log2(fc)
    sa = (cmp[f"padj_{name_a}"] < fdr) & (cmp[f"log2FC_{name_a}"].abs() > lf)
    sb = (cmp[f"padj_{name_b}"] < fdr) & (cmp[f"log2FC_{name_b}"].abs() > lf)
    sa, sb = sa.fillna(False), sb.fillna(False)
    interaction = (cmp["padj_interaction"] < fdr).fillna(False)
    a_only, b_only, shared = sa & ~sb, sb & ~sa, sa & sb
    same_dir = shared & (np.sign(cmp[f"log2FC_{name_a}"]) == np.sign(cmp[f"log2FC_{name_b}"]))
    out = {
        f"significant_{name_a}": int(sa.sum()),
        f"significant_{name_b}": int(sb.sum()),
        "shared": int(shared.sum()),
        "shared_same_direction": int(same_dir.sum()),
        f"{name_a}_only": int(a_only.sum()),
        f"{name_b}_only": int(b_only.sum()),
        f"{name_a}_only_with_interaction": int((a_only & interaction).sum()),
        f"{name_b}_only_with_interaction": int((b_only & interaction).sum()),
        "elements_with_interaction": int(interaction.sum()),
        "spearman_effect_sizes": float(
            ss.spearmanr(cmp[f"log2FC_{name_a}"], cmp[f"log2FC_{name_b}"],
                         nan_policy="omit")[0]),
    }
    only = out[f"{name_a}_only"] + out[f"{name_b}_only"]
    supported = out[f"{name_a}_only_with_interaction"] + out[f"{name_b}_only_with_interaction"]
    out["fraction_of_differences_supported"] = float(supported / only) if only else np.nan
    if only and supported / only < 0.25:
        out["interpretation"] = (
            f"only {supported}/{only} ({100*supported/only:.0f}%) of the contrast-specific "
            "calls survive an interaction test -- most of the apparent difference between "
            "these contrasts is power, not biology")
    else:
        out["interpretation"] = (
            f"{supported}/{only} contrast-specific calls are supported by an interaction "
            "test -- the contrasts genuinely differ at these elements")
    LOG.info(out["interpretation"])
    return out


def meta_analyse(results: dict, method: str = "inverse_variance") -> pd.DataFrame:
    """Pool several contrasts into one effect size per element.

    Use when the contrasts are replicates of the same question (e.g. two cohorts),
    NOT when they are different conditions you want to contrast -- for that use
    :func:`compare_contrasts`.
    """
    names = list(results)
    idx = None
    for r in results.values():
        idx = r.index if idx is None else idx.intersection(r.index)
    L = pd.DataFrame({n: results[n].loc[idx, "log2FC"] for n in names})
    S = pd.DataFrame({n: _se(results[n]).loc[idx] for n in names})
    if method == "inverse_variance":
        W = 1 / (S ** 2)
        W = W.replace([np.inf, -np.inf], np.nan)
        eff = (L * W).sum(axis=1) / W.sum(axis=1)
        se = np.sqrt(1 / W.sum(axis=1))
    elif method == "mean":
        eff, se = L.mean(axis=1), S.mean(axis=1) / np.sqrt(len(names))
    else:
        raise ValueError("method must be 'inverse_variance' or 'mean'")
    z = eff / se
    p = 2 * ss.norm.sf(np.abs(z))
    # Cochran's Q for heterogeneity
    W = 1 / (S ** 2)
    Q = (W.mul((L.sub(eff, axis=0)) ** 2)).sum(axis=1)
    df = len(names) - 1
    p_het = ss.chi2.sf(Q, df) if df > 0 else np.full(len(idx), np.nan)
    out = pd.DataFrame({"log2FC": eff, "se": se, "stat": z, "pvalue": p,
                        "Q": Q, "pvalue_heterogeneity": p_het}, index=idx)
    out["padj"] = bh_fdr(out["pvalue"].to_numpy())
    out["baseMean"] = pd.concat([results[n].loc[idx, "baseMean"] for n in names],
                                axis=1).mean(axis=1)
    out.attrs["engine"] = f"meta:{method}"
    out.attrs["contrast"] = results[names[0]].attrs.get("contrast", "")
    return out


def replication_rate(discovery: pd.DataFrame, validation: pd.DataFrame,
                     fc: float = 1.5, fdr: float = 0.1,
                     loose_p: float = 0.05) -> dict:
    """Do discovery hits replicate, in direction and at a relaxed threshold?"""
    idx = discovery.index.intersection(validation.index)
    d, v = discovery.loc[idx], validation.loc[idx]
    lf = np.log2(fc)
    hits = ((d["padj"] < fdr) & (d["log2FC"].abs() > lf)).fillna(False)
    if hits.sum() == 0:
        return {"n_discovery_hits": 0}
    same_dir = np.sign(d.loc[hits, "log2FC"]) == np.sign(v.loc[hits, "log2FC"])
    loose = (v.loc[hits, "pvalue"] < loose_p) & same_dir
    strict = ((v.loc[hits, "padj"] < fdr) & (v.loc[hits, "log2FC"].abs() > lf) & same_dir)
    return {
        "n_discovery_hits": int(hits.sum()),
        "same_direction": int(same_dir.sum()),
        "same_direction_rate": float(same_dir.mean()),
        "replicated_loose": int(loose.sum()),
        "replicated_strict": int(strict.sum()),
        "replication_rate_loose": float(loose.mean()),
        "replication_rate_strict": float(strict.mean()),
    }
