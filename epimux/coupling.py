"""Cross-assay coupling, with sign conventions enforced rather than assumed.

The headline number in a multi-omic paper is usually a cross-layer correlation
("do accessibility and activity change together?").  That number is trivially
easy to get wrong:

* if the two assays were contrasted in opposite directions, the correlation
  flips sign and reads as "decoupling";
* if replicate reliability is low, the correlation is attenuated toward zero
  (regression dilution) and can be reported as "no relationship";
* a positive correlation can be induced by shared shrinkage if effect sizes are
  compared without stratifying on expression.

:func:`couple` refuses to run unless both results were produced from contrasts
with the same reference/test orientation, reports the attenuation-corrected
estimate, and can stratify by expression as a built-in control.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from scipy import stats as ss

from .utils import get_logger

LOG = get_logger()

__all__ = ["CouplingResult", "couple", "classify_elements", "concordance"]


@dataclass
class CouplingResult:
    assay_x: str
    assay_y: str
    spearman: float
    pearson: float
    pvalue: float
    n: int
    dual_significant: int
    same_direction: int
    opposite_direction: int
    strata: pd.DataFrame | None = None
    attenuation_corrected: float | None = None
    notes: list = field(default_factory=list)

    @property
    def interpretation(self) -> str:
        if self.n == 0:
            return "no overlapping features"
        if abs(self.spearman) < 0.05:
            return "no detectable relationship"
        word = "COUPLED (same direction)" if self.spearman > 0 else "DECOUPLED (opposite)"
        return f"{word}, rho = {self.spearman:+.3f}"

    def __repr__(self):
        s = [f"Coupling({self.assay_x} vs {self.assay_y})",
             f"  {self.interpretation}   P = {self.pvalue:.2e}, n = {self.n:,}",
             f"  dual-significant: {self.dual_significant}  "
             f"(same direction {self.same_direction}, opposite {self.opposite_direction})"]
        if self.attenuation_corrected is not None:
            s.append(f"  attenuation-corrected rho ~ {self.attenuation_corrected:+.3f}")
        for n in self.notes:
            s.append(f"  ! {n}")
        return "\n".join(s)


def _orientation(res: pd.DataFrame) -> str | None:
    c = res.attrs.get("contrast")
    if not c:
        return None
    # Contrast repr embeds "log2(test/ref)"
    if "log2(" in c:
        return c.split("log2(")[1].split(")")[0]
    return None


def couple(res_x: pd.DataFrame, res_y: pd.DataFrame,
           name_x: str = "X", name_y: str = "Y",
           fc: float = 1.5, fdr: float = 0.1,
           reliability: tuple | None = None,
           stratify: bool = True, n_strata: int = 4) -> CouplingResult:
    """Correlate two differential results element-wise.

    Both inputs must be indexed by the *same* element identifiers (that is what
    :class:`epimux.core.Dataset` guarantees).
    """
    notes = []
    ox, oy = _orientation(res_x), _orientation(res_y)
    if ox and oy and ox != oy:
        raise ValueError(
            f"contrast orientation mismatch: {name_x} is {ox} but {name_y} is {oy}. "
            "Correlating these would flip the sign and read as decoupling. "
            "Rebuild both with the same ref/test.")

    idx = res_x.index.intersection(res_y.index)
    x = res_x.loc[idx, "log2FC"]
    y = res_y.loc[idx, "log2FC"]
    ok = x.notna() & y.notna()
    x, y = x[ok], y[ok]
    if len(x) < 10:
        return CouplingResult(name_x, name_y, np.nan, np.nan, np.nan, len(x), 0, 0, 0,
                              notes=["too few overlapping elements"])

    rho, p = ss.spearmanr(x, y)
    pear = float(ss.pearsonr(x, y)[0])

    lfc = np.log2(fc)
    sx = (res_x.loc[x.index, "padj"] < fdr) & (x.abs() > lfc)
    sy = (res_y.loc[y.index, "padj"] < fdr) & (y.abs() > lfc)
    dual = sx & sy
    same = int((dual & (np.sign(x) == np.sign(y))).sum())
    opp = int((dual & (np.sign(x) != np.sign(y))).sum())

    strata = None
    if stratify and len(x) >= 200:
        base = (res_x.loc[x.index, "baseMean"].rank(pct=True))
        bins = np.clip((base * n_strata).astype(int), 0, n_strata - 1)
        rows = []
        for b in range(n_strata):
            m = bins == b
            if m.sum() < 30:
                continue
            r_b = ss.spearmanr(x[m], y[m])[0]
            rows.append({"stratum": b, "n": int(m.sum()), "spearman": r_b})
        strata = pd.DataFrame(rows)
        if len(strata) > 1:
            spread = strata["spearman"].max() - strata["spearman"].min()
            if np.sign(strata["spearman"]).nunique() > 1:
                notes.append("correlation changes sign across expression strata -- "
                             "the global estimate is not a stable summary")
            elif spread > 0.25:
                notes.append(f"correlation varies with expression (spread {spread:.2f}); "
                             "report the stratified values")

    corrected = None
    if reliability:
        rx, ry = reliability
        denom = np.sqrt(max(rx, 1e-6) * max(ry, 1e-6))
        corrected = float(rho / denom)
        if denom < 0.7:
            notes.append(f"low reliability attenuates rho by {denom:.2f}x; "
                         "the observed value is a lower bound")

    if opp == 0 and dual.sum() > 0:
        notes.append(f"no element changes in opposite directions -- "
                     f"there is no decoupled subpopulation at FC>{fc}/FDR<{fdr}")

    return CouplingResult(name_x, name_y, float(rho), pear, float(p), int(len(x)),
                          int(dual.sum()), same, opp, strata, corrected, notes)


def classify_elements(results: dict, fc: float = 1.5, fdr: float = 0.1) -> pd.DataFrame:
    """Label each element by its multi-layer response.

    Categories are derived from *significance*, never from the sign of a single
    noisy difference -- sign-only state calls were the original source of
    irreproducible "decoupled enhancer" lists.
    """
    lfc = np.log2(fc)
    idx = None
    for r in results.values():
        idx = r.index if idx is None else idx.intersection(r.index)
    out = pd.DataFrame(index=idx)
    for name, r in results.items():
        v = r.loc[idx, "log2FC"]
        s = (r.loc[idx, "padj"] < fdr) & (v.abs() > lfc)
        out[f"{name}_lfc"] = v
        out[f"{name}_sig"] = s.fillna(False)
        out[f"{name}_dir"] = np.where(~s.fillna(False), "ns",
                                      np.where(v > 0, "up", "down"))
    dirs = out[[c for c in out.columns if c.endswith("_dir")]]
    n_up = (dirs == "up").sum(1)
    n_dn = (dirs == "down").sum(1)
    out["n_sig"] = n_up + n_dn
    out["state"] = np.select(
        [(n_up >= 2) & (n_dn == 0), (n_dn >= 2) & (n_up == 0),
         (n_up >= 1) & (n_dn >= 1), out["n_sig"] == 1],
        ["coordinated_up", "coordinated_down", "discordant", "single_layer"],
        default="unchanged")
    return out


def concordance(cls: pd.DataFrame) -> pd.Series:
    """Summary counts of the multi-layer states."""
    return cls["state"].value_counts()
