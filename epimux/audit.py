"""The audit layer -- automated checks that catch the failure modes that
silently corrupt differential genomics.

Every check in here exists because it caught a real, published-grade error:

``check_direction``
    A contrast built from an R-style factor silently computed ``log2(WT/KO)``
    because factor levels sort alphabetically.  Every direction in the analysis
    was reversed and nothing downstream complained.  This check compares the
    reported ``log2FC`` against raw normalized means and fails loudly.

``positive_control``
    A bigWig/limma pipeline reported "no change" for a histone mark.  It was a
    false negative of a weak method: the same pipeline could not detect a
    cell-type difference either.  A differential pipeline that cannot find a
    difference you *know* exists cannot support a null result.

``null_contrast``
    Splitting replicates *within* one group must yield ~no hits.  If it does
    not, the "significant" features in the real contrast are noise.

``efficiency_balance``
    Systematically different ChIP efficiency (FRiP) between genotypes can
    manufacture a directional bias that survives depth normalization.  This
    reports the imbalance and, crucially, whether the observed effect runs
    *with* the technical bias (dangerous) or *against* it (conservative).

``replicate_reliability``
    Spearman-Brown reliability of a group mean; low reliability means a
    difference of two group means is mostly noise, and cross-assay correlations
    will be attenuated toward zero.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from scipy import stats as ss

from .utils import Contrast, cpm, get_logger

LOG = get_logger()

__all__ = ["AuditResult", "check_direction", "positive_control", "null_contrast",
           "efficiency_balance", "replicate_reliability", "AuditReport"]


@dataclass
class AuditResult:
    name: str
    status: str                       # "pass" | "warn" | "fail"
    summary: str
    detail: dict = field(default_factory=dict)

    @property
    def symbol(self):
        return {"pass": "PASS", "warn": "WARN", "fail": "FAIL"}[self.status]

    def __repr__(self):
        return f"[{self.symbol}] {self.name}: {self.summary}"


@dataclass
class AuditReport:
    results: list = field(default_factory=list)

    def add(self, r: AuditResult):
        self.results.append(r)
        LOG.info(repr(r))
        return r

    @property
    def failed(self):
        return [r for r in self.results if r.status == "fail"]

    def to_frame(self) -> pd.DataFrame:
        return pd.DataFrame([{"check": r.name, "status": r.status,
                              "summary": r.summary} for r in self.results])

    def __repr__(self):
        if not self.results:
            return "AuditReport(empty)"
        return "\n".join(repr(r) for r in self.results)


# --------------------------------------------------------------------------
def check_direction(result: pd.DataFrame, counts: pd.DataFrame,
                    contrast: Contrast, kind: str = "count",
                    top_n: int = 2000, min_corr: float = 0.5) -> AuditResult:
    """Verify the sign of ``log2FC`` against raw normalized values.

    This is the single most important check in the package.  It recomputes a
    naive ``log2(mean(test) / mean(ref))`` from the raw matrix and correlates it
    with the reported effect size.  A strong *negative* correlation means the
    contrast was built backwards.
    """
    sig = result.dropna(subset=["log2FC"])
    if "padj" in sig:
        ranked = sig.reindex(sig["padj"].sort_values().index)
    else:
        ranked = sig.reindex(sig["log2FC"].abs().sort_values(ascending=False).index)
    idx = ranked.index[:top_n]
    idx = [i for i in idx if i in counts.index]
    if len(idx) < 10:
        return AuditResult("check_direction", "warn",
                           "too few overlapping features to verify direction")

    sub = counts.loc[idx]
    if kind == "count":
        norm = pd.DataFrame(cpm(sub.to_numpy(dtype=float)),
                            index=sub.index, columns=sub.columns)
        ref = norm[contrast.ref_samples].mean(1)
        test = norm[contrast.test_samples].mean(1)
        naive = np.log2((test + 0.5) / (ref + 0.5))
    else:                                    # already on a linear/ratio scale
        ref = sub[contrast.ref_samples].mean(1)
        test = sub[contrast.test_samples].mean(1)
        naive = test - ref

    rep = result.loc[idx, "log2FC"].to_numpy(dtype=float)
    ok = np.isfinite(rep) & np.isfinite(naive.to_numpy())
    if ok.sum() < 10:
        return AuditResult("check_direction", "warn", "not enough finite values")
    r = float(ss.spearmanr(rep[ok], naive.to_numpy()[ok])[0])
    detail = {"spearman": r, "n": int(ok.sum()), "contrast": repr(contrast)}

    if r <= -min_corr:
        return AuditResult(
            "check_direction", "fail",
            f"log2FC sign is REVERSED (corr with raw = {r:+.3f}). "
            f"Reported values are log2({contrast.ref}/{contrast.test}), not "
            f"log2({contrast.test}/{contrast.ref}).", detail)
    if r < min_corr:
        return AuditResult("check_direction", "warn",
                           f"weak agreement with raw values (corr {r:+.3f}); "
                           "check normalization", detail)
    return AuditResult("check_direction", "pass",
                       f"sign verified against raw values (corr {r:+.3f})", detail)


# --------------------------------------------------------------------------
def positive_control(counts: pd.DataFrame, contrast: Contrast, de_fn,
                     fc: float = 1.5, fdr: float = 0.1,
                     min_frac: float = 0.02) -> AuditResult:
    """Run the pipeline on a comparison where a difference certainly exists.

    ``contrast`` should be something like cell type, not genotype.  If this
    finds essentially nothing, the pipeline lacks power and *no* null result
    from it is interpretable.
    """
    res = de_fn(counts, contrast)
    lfc = np.log2(fc)
    sig = ((res["padj"] < fdr) & (res["log2FC"].abs() > lfc)).sum()
    tested = int(res["log2FC"].notna().sum())
    frac = sig / max(tested, 1)
    detail = {"significant": int(sig), "tested": tested, "fraction": float(frac),
              "contrast": repr(contrast)}
    if frac < min_frac / 4:
        return AuditResult("positive_control", "fail",
                           f"only {sig}/{tested} ({frac:.2%}) differ for a comparison that "
                           "must differ -- pipeline lacks power; null results are "
                           "uninterpretable", detail)
    if frac < min_frac:
        return AuditResult("positive_control", "warn",
                           f"weak power: {sig}/{tested} ({frac:.2%})", detail)
    return AuditResult("positive_control", "pass",
                       f"detects {sig:,}/{tested:,} ({frac:.1%}) -- pipeline has power",
                       detail)


# --------------------------------------------------------------------------
def null_contrast(counts: pd.DataFrame, samples: list, de_fn,
                  fc: float = 1.5, fdr: float = 0.1,
                  max_hits: int = 25, seed: int = 0) -> AuditResult:
    """Split replicates *within* one biological group and run the same test.

    Any hit here is a false positive by construction.
    """
    if len(samples) < 3:
        return AuditResult("null_contrast", "warn",
                           f"need >=3 replicates in a group, got {len(samples)}")
    rng = np.random.default_rng(seed)
    s = list(samples)
    rng.shuffle(s)
    k = len(s) // 2
    a, b = s[:k], s[k:]
    ctr = Contrast(ref="nullA", test="nullB", group={"nullA": a, "nullB": b})
    res = de_fn(counts, ctr)
    lfc = np.log2(fc)
    sig = int(((res["padj"] < fdr) & (res["log2FC"].abs() > lfc)).sum())
    detail = {"significant": sig, "split": {"A": a, "B": b}}
    if sig > max_hits * 4:
        return AuditResult("null_contrast", "fail",
                           f"{sig} 'significant' features between replicates of the "
                           "same group -- real hits are not distinguishable from noise",
                           detail)
    if sig > max_hits:
        return AuditResult("null_contrast", "warn",
                           f"{sig} hits in a null contrast; interpret modest effects "
                           "with care", detail)
    return AuditResult("null_contrast", "pass",
                       f"only {sig} hits within-group -- false-positive rate is low",
                       detail)


# --------------------------------------------------------------------------
def efficiency_balance(frip: dict, contrast: Contrast,
                       observed_direction: str | None = None,
                       max_ratio: float = 1.2) -> AuditResult:
    """Compare signal-to-background (e.g. FRiP) between the two groups.

    ``frip`` maps sample -> fraction of reads in peaks (0-1 or percent).
    ``observed_direction`` is ``"up"``/``"down"``: the dominant direction of the
    real contrast.  The check reports whether the technical bias would *create*
    that direction (dangerous) or *oppose* it (the result is conservative).
    """
    ref = np.array([frip[s] for s in contrast.ref_samples if s in frip], dtype=float)
    test = np.array([frip[s] for s in contrast.test_samples if s in frip], dtype=float)
    if ref.size == 0 or test.size == 0:
        return AuditResult("efficiency_balance", "warn", "FRiP not available for both groups")
    mr, mt = float(ref.mean()), float(test.mean())
    ratio = mr / mt if mt else np.inf
    detail = {"ref_mean": mr, "test_mean": mt, "ref_over_test": ratio,
              "ref_values": ref.tolist(), "test_values": test.tolist()}

    balanced = (1 / max_ratio) <= ratio <= max_ratio
    if balanced:
        return AuditResult("efficiency_balance", "pass",
                           f"efficiency balanced ({contrast.ref} {mr:.3g} vs "
                           f"{contrast.test} {mt:.3g}; ratio {ratio:.2f}x)", detail)

    bias = "down" if ratio > 1 else "up"   # lower efficiency in test -> signal looks lower
    if observed_direction is None:
        return AuditResult("efficiency_balance", "warn",
                           f"efficiency imbalance {ratio:.2f}x ({contrast.ref} {mr:.3g} vs "
                           f"{contrast.test} {mt:.3g}); could bias the {bias} direction",
                           detail)
    detail["technical_bias_direction"] = bias
    if observed_direction == bias:
        return AuditResult("efficiency_balance", "fail",
                           f"efficiency imbalance {ratio:.2f}x biases toward '{bias}', which is "
                           f"the SAME direction as the observed effect -- the result may be "
                           "technical. Re-test on an efficiency-matched subset.", detail)
    return AuditResult("efficiency_balance", "pass",
                       f"efficiency imbalance {ratio:.2f}x biases toward '{bias}', OPPOSITE to "
                       f"the observed '{observed_direction}' effect -- the result is conservative",
                       detail)


# --------------------------------------------------------------------------
def replicate_reliability(values: pd.DataFrame, groups: dict,
                          min_reliability: float = 0.5) -> AuditResult:
    """Spearman-Brown reliability of each group mean.

    Low reliability attenuates every downstream correlation by
    ``sqrt(rel_x * rel_y)`` -- the reason single-replicate tracks produce
    correlations with the wrong magnitude and sometimes the wrong sign.
    """
    out = {}
    for g, samples in groups.items():
        samples = [s for s in samples if s in values.columns]
        if len(samples) < 2:
            out[g] = np.nan
            continue
        M = values[samples].to_numpy(dtype=float)
        ok = np.isfinite(M).all(1)
        M = M[ok]
        cors = []
        for i in range(len(samples)):
            for j in range(i + 1, len(samples)):
                cors.append(ss.spearmanr(M[:, i], M[:, j])[0])
        rbar = float(np.nanmean(cors))
        k = len(samples)
        out[g] = k * rbar / (1 + (k - 1) * rbar) if np.isfinite(rbar) else np.nan
    worst = np.nanmin(list(out.values())) if out else np.nan
    detail = {"reliability": out}
    if not np.isfinite(worst):
        return AuditResult("replicate_reliability", "warn", "could not estimate", detail)
    txt = ", ".join(f"{g} {v:.2f}" for g, v in out.items())
    if worst < min_reliability:
        return AuditResult("replicate_reliability", "warn",
                           f"low reliability ({txt}); cross-assay correlations will be "
                           f"attenuated by ~{np.sqrt(worst):.2f}x", detail)
    return AuditResult("replicate_reliability", "pass",
                       f"group means are reliable ({txt})", detail)
