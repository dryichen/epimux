"""The Dataset: one reference element set, many assays, one sign convention.

    ds = ep.Dataset(reference="enhancers.bed", genome="mm10")
    ds.add_counts("ATAC",    "atac_counts.txt")
    ds.add_counts("H3K27ac", "h3k_counts.txt")
    ds.add_methyl("WGBS",    "wgbs/*.cov")
    ds.add_hic("HiC", {"WT": "wt.mcool", "KO": "ko.mcool"})

    ds.set_design({"WT": [...], "KO": [...]})
    res = ds.differential(ref="WT", test="KO")     # log2FC = log2(KO/WT), always
    ds.audit()                                     # direction / power / null / efficiency
    ds.coupling("ATAC", "H3K27ac")

Design notes
------------
* Results are keyed by **reference element index**, so every layer is directly
  comparable without re-running interval arithmetic.
* The contrast is an object, not a factor.  R's ``factor()`` sorts levels
  alphabetically, which silently reverses ``c("WT","KO")`` contrasts; that class
  of bug is impossible here and is additionally checked at runtime.
* ``differential()`` runs :func:`epimux.audit.check_direction` by default and
  raises if the reported sign disagrees with the raw data.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from . import audit as _audit
from . import coupling as _coupling
from . import modules as _modules
from .assays import (Assay, CountAssay, HiCAssay, MethylAssay, SignalAssay,
                     read_featurecounts)
from .utils import Contrast, as_intervals, get_logger, read_bed

LOG = get_logger()

__all__ = ["Dataset"]


class Dataset:
    """Multi-assay container anchored on a shared element reference."""

    def __init__(self, reference, genome: str = "unknown", name: str = "dataset"):
        if isinstance(reference, (str, bytes)) or hasattr(reference, "__fspath__"):
            ref = read_bed(reference)
        else:
            ref = as_intervals(reference)
        self.reference = ref.reset_index(drop=True)
        self.genome = genome
        self.name = name
        self.assays: dict[str, Assay] = {}
        self.design: dict[str, list] = {}
        self.results: dict[str, pd.DataFrame] = {}
        self.audit_report = _audit.AuditReport()
        self._contrast: Contrast | None = None
        LOG.info(f"Dataset '{name}': {len(self.reference):,} reference elements ({genome})")

    # ------------------------------------------------------------------ IO
    def add_assay(self, assay: Assay):
        assay.map_to(self.reference)
        self.assays[assay.name] = assay
        return assay

    def add_counts(self, name, path=None, intervals=None, matrix=None, **kw):
        if path is not None:
            a = CountAssay.from_featurecounts(name, path, **kw)
        else:
            a = CountAssay(name=name, intervals=as_intervals(intervals), matrix=matrix, **kw)
        return self.add_assay(a)

    def add_methyl(self, name, pattern, **kw):
        a = MethylAssay.from_bismark_cov(name, pattern, self.reference, **kw)
        self.assays[name] = a           # already element-indexed
        a._ref_map = pd.Series(np.arange(len(self.reference)))
        return a

    def add_signal(self, name, files, **kw):
        a = SignalAssay.from_bigwigs(name, files, self.reference, **kw)
        self.assays[name] = a
        a._ref_map = pd.Series(np.arange(len(self.reference)))
        return a

    def add_hic(self, name, files):
        a = HiCAssay.from_mcools(name, files, self.reference)
        self.assays[name] = a
        return a

    # -------------------------------------------------------------- design
    def set_design(self, groups: dict):
        """``{"WT": [samples...], "KO": [samples...]}``"""
        self.design = {k: list(v) for k, v in groups.items()}
        return self

    def contrast(self, ref: str, test: str) -> Contrast:
        for g in (ref, test):
            if g not in self.design:
                raise KeyError(f"group '{g}' not in design {list(self.design)}")
        return Contrast(ref=ref, test=test, group=self.design)

    # --------------------------------------------------------- differential
    def differential(self, ref: str, test: str, assays: list | None = None,
                     verify_direction: bool = True, strict: bool = True,
                     **kw) -> dict:
        """Run per-assay differential analysis with one pinned direction."""
        ctr = self.contrast(ref, test)
        self._contrast = ctr
        LOG.info(f"differential: {ctr}")
        names = assays or [n for n, a in self.assays.items() if a.kind != "hic"]
        out = {}
        for n in names:
            a = self.assays[n]
            sub = Contrast(ref=ref, test=test,
                           group={g: [s for s in v if s in a.samples]
                                  for g, v in self.design.items()})
            if not sub.ref_samples or not sub.test_samples:
                LOG.warning(f"{n}: no samples for this contrast, skipped")
                continue
            r = a.differential(sub, **kw)
            r.attrs.setdefault("contrast", repr(sub))

            if verify_direction and a.kind in ("count", "methylation"):
                mat = a.matrix if a.kind == "count" else a.rates()
                chk = _audit.check_direction(r, mat, sub,
                                             kind="count" if a.kind == "count" else "rate")
                chk.name = f"check_direction[{n}]"
                self.audit_report.add(chk)
                if chk.status == "fail":
                    if strict:
                        raise RuntimeError(
                            f"{n}: {chk.summary} -- refusing to return reversed results. "
                            "Pass strict=False to override.")
                    LOG.warning(f"{n}: {chk.summary}")

            out[n] = self._to_elements(a, r)
        self.results.update(out)
        return out

    def _to_elements(self, assay: Assay, res: pd.DataFrame) -> pd.DataFrame:
        """Collapse a feature-level result onto reference elements."""
        if assay.kind in ("methylation", "signal"):
            r = res.copy()
            r.index = np.arange(len(self.reference))
            return r
        cols = {}
        rank = res["padj"].to_numpy()
        for c in ["baseMean", "log2FC", "stat", "pvalue", "padj"]:
            cols[c] = assay.to_elements(res[c], self.reference,
                                        agg="strongest", rank_by=rank)
        out = pd.DataFrame(cols)
        out.attrs.update(res.attrs)
        return out

    # ---------------------------------------------------------------- audit
    def audit(self, positive_control: tuple | None = None,
              null_group: str | None = None, frip: dict | None = None,
              assays: list | None = None) -> _audit.AuditReport:
        """Run the full battery of checks.

        ``positive_control`` -- ``(ref, test)`` for a comparison that must differ
        (e.g. two cell types).  ``null_group`` -- a group whose replicates are
        split against each other.
        """
        names = assays or [n for n, a in self.assays.items() if a.kind == "count"]
        for n in names:
            a = self.assays[n]
            if positive_control:
                pr, pt = positive_control
                sub = Contrast(ref=pr, test=pt,
                               group={g: [s for s in v if s in a.samples]
                                      for g, v in self.design.items()})
                if sub.ref_samples and sub.test_samples:
                    r = _audit.positive_control(
                        a.matrix, sub, lambda c, ct: a.differential(ct))
                    r.name = f"positive_control[{n}]"
                    self.audit_report.add(r)
            if null_group:
                samples = [s for s in self.design.get(null_group, []) if s in a.samples]
                r = _audit.null_contrast(a.matrix, samples,
                                         lambda c, ct: a.differential(ct))
                r.name = f"null_contrast[{n}:{null_group}]"
                self.audit_report.add(r)
            f = frip or a.meta.get("frip")
            if f and self._contrast:
                obs = None
                if n in self.results:
                    lf = self.results[n]["log2FC"]
                    sig = self.results[n]["padj"] < 0.1
                    if sig.sum() > 0:
                        obs = "up" if (lf[sig] > 0).mean() > 0.5 else "down"
                sub = Contrast(ref=self._contrast.ref, test=self._contrast.test,
                               group={g: [s for s in v if s in a.samples]
                                      for g, v in self.design.items()})
                r = _audit.efficiency_balance(f, sub, observed_direction=obs)
                r.name = f"efficiency_balance[{n}]"
                self.audit_report.add(r)
            r = _audit.replicate_reliability(
                a.norm() if hasattr(a, "norm") else a.matrix,
                {g: [s for s in v if s in a.samples] for g, v in self.design.items()})
            r.name = f"replicate_reliability[{n}]"
            self.audit_report.add(r)
        return self.audit_report

    # -------------------------------------------------------------- analysis
    def coupling(self, x: str, y: str, **kw) -> _coupling.CouplingResult:
        if x not in self.results or y not in self.results:
            raise KeyError("run differential() first")
        rel = None
        rr = [r for r in self.audit_report.results
              if r.name.startswith("replicate_reliability")]
        if len(rr) >= 1:
            vals = [v for r in rr for v in r.detail.get("reliability", {}).values()
                    if np.isfinite(v)]
            if vals:
                rel = (float(np.mean(vals)), float(np.mean(vals)))
        return _coupling.couple(self.results[x], self.results[y], x, y,
                                reliability=rel, **kw)

    def classify(self, assays: list | None = None, **kw) -> pd.DataFrame:
        names = assays or list(self.results)
        return _coupling.classify_elements({n: self.results[n] for n in names}, **kw)

    def modules(self, layers: dict, k: int = 6, **kw) -> _modules.ModuleResult:
        return _modules.find_modules(layers, k=k, **kw)

    def link_genes(self, tss, method: str = "abc", hic: str | None = None,
                   sample: str | None = None, activity: pd.Series | None = None, **kw):
        from . import linking
        if method == "nearest":
            return linking.nearest_gene(self.reference, tss, **kw)
        hic_assay = self.assays[hic] if hic else next(
            (a for a in self.assays.values() if a.kind == "hic"), None)
        if hic_assay is None:
            raise ValueError("ABC linking needs a Hi-C assay; use method='nearest'")
        sample = sample or next(iter(hic_assay.coolers))
        return linking.abc_link(self.reference, tss, hic_assay, sample,
                                activity=activity, **kw)

    # ----------------------------------------------------------------- misc
    def summary(self) -> pd.DataFrame:
        rows = []
        for n, a in self.assays.items():
            rows.append({"assay": n, "kind": a.kind,
                         "features": len(a.intervals),
                         "samples": len(a.samples) if a.kind != "hic" else len(a.coolers)})
        return pd.DataFrame(rows)

    def significant(self, assay: str, fc: float = 1.5, fdr: float = 0.1) -> pd.DataFrame:
        r = self.results[assay]
        thr = np.log2(fc) if r.attrs.get("value_kind") != "difference" else 0.1
        m = (r["padj"] < fdr) & (r["log2FC"].abs() > thr)
        return r[m.fillna(False)]

    def report(self, path="epimux_report.html", **kw):
        from .report import html_report
        return html_report(self, path, **kw)

    def __repr__(self):
        return (f"Dataset('{self.name}', {len(self.reference):,} elements, "
                f"assays={list(self.assays)})")
