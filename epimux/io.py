"""Import/export and interoperability.

Analyses outlive the session that produced them, so a Dataset must be able to
leave epimux without losing the two things that make it trustworthy: the pinned
contrast, and the audit record.  Every export carries both.
"""
from __future__ import annotations

import json
import os

import numpy as np
import pandas as pd

from .utils import as_intervals, get_logger

LOG = get_logger()

__all__ = ["to_anndata", "to_mudata", "export_bed", "export_results",
           "save_dataset", "load_results"]


# --------------------------------------------------------------------------
def to_anndata(dataset, assay: str):
    """One assay as AnnData (obs = samples, var = reference elements)."""
    try:
        import anndata as ad
    except ImportError as e:  # pragma: no cover
        raise ImportError("anndata required; pip install anndata") from e
    a = dataset.assays[assay]
    X = a.matrix.T.to_numpy(dtype=float)
    obs = pd.DataFrame(index=pd.Index(a.samples, name="sample"))
    for g, samples in dataset.design.items():
        obs.loc[[s for s in a.samples if s in samples], "group"] = g
    var = a.intervals.reset_index(drop=True).copy()
    var.index = var["chrom"].astype(str) + ":" + var["start"].astype(str) + "-" + var["end"].astype(str)
    adata = ad.AnnData(X=X, obs=obs, var=var)
    if assay in dataset.results:
        r = dataset.results[assay]
        adata.uns["differential"] = r.reset_index(drop=True).to_dict("list")
        adata.uns["contrast"] = r.attrs.get("contrast", "")
    adata.uns["audit"] = dataset.audit_report.to_frame().to_dict("list")
    adata.uns["epimux_version"] = _version()
    return adata


def to_mudata(dataset):
    """All assays as MuData, one modality per assay."""
    try:
        import mudata as md
    except ImportError as e:  # pragma: no cover
        raise ImportError("mudata required; pip install mudata") from e
    mods = {n: to_anndata(dataset, n) for n, a in dataset.assays.items()
            if a.kind != "hic"}
    mdata = md.MuData(mods)
    mdata.uns["audit"] = dataset.audit_report.to_frame().to_dict("list")
    mdata.uns["contrast"] = repr(dataset._contrast) if dataset._contrast else ""
    return mdata


# --------------------------------------------------------------------------
def export_bed(dataset, assay: str, path: str, fc: float = 1.5, fdr: float = 0.1,
               direction: str | None = None, name_prefix: str = "") -> str:
    """Significant elements as BED, score = 1000*min(1,|log2FC|/3), strand-free.

    ``direction`` : ``"up"``, ``"down"`` or ``None`` for both.
    """
    r = dataset.results[assay]
    thr = np.log2(fc) if r.attrs.get("value_kind") != "difference" else 0.1
    m = (r["padj"] < fdr) & (r["log2FC"].abs() > thr)
    if direction == "up":
        m &= r["log2FC"] > 0
    elif direction == "down":
        m &= r["log2FC"] < 0
    m = m.fillna(False)
    ref = dataset.reference.loc[m[m].index].copy()
    lfc = r.loc[m[m].index, "log2FC"]
    ref["name"] = [f"{name_prefix}{assay}_{i}" for i in range(len(ref))]
    ref["score"] = (1000 * np.minimum(1.0, lfc.abs() / 3)).astype(int).to_numpy()
    ref["strand"] = "."
    ref["log2FC"] = lfc.to_numpy()
    ref["padj"] = r.loc[m[m].index, "padj"].to_numpy()
    ref[["chrom", "start", "end", "name", "score", "strand", "log2FC", "padj"]] \
        .to_csv(path, sep="\t", header=False, index=False)
    LOG.info(f"wrote {len(ref):,} {direction or 'changed'} elements to {path}")
    return os.path.abspath(path)


def export_results(dataset, out_dir: str, fc: float = 1.5, fdr: float = 0.1) -> dict:
    """Write every result table, the audit record and a machine-readable manifest."""
    os.makedirs(out_dir, exist_ok=True)
    written = {}
    for name, r in dataset.results.items():
        p = os.path.join(out_dir, f"{name}_differential.tsv")
        out = dataset.reference.join(r)
        out.to_csv(p, sep="\t", index=False)
        written[name] = p
    ap = os.path.join(out_dir, "audit.tsv")
    dataset.audit_report.to_frame().to_csv(ap, sep="\t", index=False)
    written["audit"] = ap

    manifest = {
        "name": dataset.name,
        "genome": dataset.genome,
        "n_reference_elements": int(len(dataset.reference)),
        "epimux_version": _version(),
        "contrast": repr(dataset._contrast) if dataset._contrast else None,
        "sign_convention": "log2FC = log2(test / ref), verified against raw values",
        "thresholds": {"fold_change": fc, "fdr": fdr},
        "assays": {n: {"kind": a.kind,
                       "n_features": int(len(a.intervals)),
                       "samples": list(a.samples) if a.kind != "hic" else list(a.coolers)}
                   for n, a in dataset.assays.items()},
        "design": dataset.design,
        "audit": {r.name: {"status": r.status, "summary": r.summary}
                  for r in dataset.audit_report.results},
        "results": {n: _summarise(r, fc, fdr) for n, r in dataset.results.items()},
    }
    mp = os.path.join(out_dir, "manifest.json")
    with open(mp, "w") as fh:
        json.dump(manifest, fh, indent=2, default=str)
    written["manifest"] = mp
    LOG.info(f"exported {len(dataset.results)} result table(s) + manifest to {out_dir}")
    return written


def _summarise(r: pd.DataFrame, fc: float, fdr: float) -> dict:
    thr = np.log2(fc) if r.attrs.get("value_kind") != "difference" else 0.1
    sig = (r["padj"] < fdr) & (r["log2FC"].abs() > thr)
    up = int((sig & (r["log2FC"] > 0)).sum())
    dn = int((sig & (r["log2FC"] < 0)).sum())
    return {"tested": int(r["log2FC"].notna().sum()), "significant": up + dn,
            "up": up, "down": dn,
            "percent_up": round(100 * up / (up + dn), 1) if up + dn else None,
            "engine": r.attrs.get("engine"), "design": r.attrs.get("design")}


def save_dataset(dataset, path: str) -> str:
    """Serialise results + audit + design to a single JSON (data stays on disk)."""
    payload = {
        "name": dataset.name, "genome": dataset.genome,
        "design": dataset.design,
        "contrast": repr(dataset._contrast) if dataset._contrast else None,
        "audit": dataset.audit_report.to_frame().to_dict("records"),
        "results": {n: r.reset_index().to_dict("list") for n, r in dataset.results.items()},
        "epimux_version": _version(),
    }
    with open(path, "w") as fh:
        json.dump(payload, fh, default=str)
    return os.path.abspath(path)


def load_results(path: str) -> dict:
    """Read back the result tables written by :func:`save_dataset`."""
    with open(path) as fh:
        payload = json.load(fh)
    return {n: pd.DataFrame(v) for n, v in payload.get("results", {}).items()}


def _version():
    from . import __version__
    return __version__
