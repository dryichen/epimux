"""Peak set construction.

How a consensus peak set is built silently determines what a differential
analysis can find:

* union of per-sample peaks maximises sensitivity but inflates the test burden;
* "present in >= k replicates" is the common default, and is **biased whenever
  the groups have unequal replicate numbers** — the group with more replicates
  wins more peaks, which manufactures apparent gains;
* a fixed reference set (an atlas) avoids both but may miss condition-specific
  elements.

:func:`consensus_peaks` implements all three and *refuses* the biased case
silently: with unequal replicates it either balances by subsampling or raises.
"""
from __future__ import annotations

import glob
import os

import numpy as np
import pandas as pd

from .utils import as_intervals, get_logger, sort_intervals

LOG = get_logger()

__all__ = ["read_narrowpeak", "merge_intervals", "consensus_peaks",
           "peak_overlap_matrix", "saf_from_intervals", "jaccard"]

NARROWPEAK_COLS = ["chrom", "start", "end", "name", "score", "strand",
                   "signalValue", "pValue", "qValue", "peak"]


def read_narrowpeak(path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", header=None, comment="#")
    df.columns = NARROWPEAK_COLS[:df.shape[1]]
    return as_intervals(df)


def merge_intervals(df: pd.DataFrame, gap: int = 0) -> pd.DataFrame:
    """Merge overlapping/nearby intervals (bedtools merge equivalent)."""
    d = sort_intervals(as_intervals(df))
    out = []
    cur = None
    for r in d.itertuples():
        if cur and r.chrom == cur[0] and r.start - cur[2] <= gap:
            cur[2] = max(cur[2], r.end)
        else:
            if cur:
                out.append(cur)
            cur = [r.chrom, r.start, r.end]
    if cur:
        out.append(cur)
    return pd.DataFrame(out, columns=["chrom", "start", "end"])


def consensus_peaks(peak_files: dict, method: str = "union", min_replicates: int = 2,
                    gap: int = 0, balance: str = "raise") -> pd.DataFrame:
    """Build a consensus peak set.

    Parameters
    ----------
    peak_files
        ``{group: [paths]}`` or ``{sample: path}``.
    method
        ``"union"``      every peak from every sample;
        ``"replicated"`` present in >= ``min_replicates`` samples of a group,
                         then the union across groups;
        ``"intersect"``  present in every sample.
    balance
        Only relevant for ``"replicated"`` with unequal replicate counts:
        ``"raise"``      refuse (default) — the naive version biases toward the
                         group with more replicates;
        ``"subsample"``  use the smallest group size for every group;
        ``"ignore"``     proceed anyway (records a warning in the result attrs).
    """
    groups = {}
    for k, v in peak_files.items():
        paths = [v] if isinstance(v, (str, bytes)) else list(v)
        groups[k] = paths

    if method == "union":
        allp = pd.concat([read_narrowpeak(p)[["chrom", "start", "end"]]
                          for ps in groups.values() for p in ps], ignore_index=True)
        out = merge_intervals(allp, gap=gap)
        out.attrs["method"] = "union"
        LOG.info(f"consensus (union): {len(out):,} peaks from "
                 f"{sum(len(v) for v in groups.values())} samples")
        return out

    sizes = {g: len(ps) for g, ps in groups.items()}
    if method == "replicated" and len(set(sizes.values())) > 1:
        msg = (f"unequal replicate counts {sizes}: a 'present in >={min_replicates}' rule "
               "favours the larger group and can manufacture apparent gains")
        if balance == "raise":
            raise ValueError(msg + ". Pass balance='subsample' or balance='ignore'.")
        if balance == "subsample":
            k = min(sizes.values())
            groups = {g: ps[:k] for g, ps in groups.items()}
            LOG.warning(msg + f"; subsampled every group to {k} replicates")
        else:
            LOG.warning(msg + "; proceeding as requested")

    per_group = []
    for g, ps in groups.items():
        pool = pd.concat([read_narrowpeak(p)[["chrom", "start", "end"]] for p in ps],
                         ignore_index=True)
        cand = merge_intervals(pool, gap=gap)
        M = peak_overlap_matrix(cand, {os.path.basename(p): read_narrowpeak(p) for p in ps})
        need = len(ps) if method == "intersect" else min_replicates
        keep = M.sum(axis=1) >= need
        LOG.info(f"consensus [{g}]: {int(keep.sum()):,}/{len(cand):,} peaks in >={need} of {len(ps)}")
        per_group.append(cand[keep.to_numpy()])
    out = merge_intervals(pd.concat(per_group, ignore_index=True), gap=gap)
    out.attrs["method"] = method
    out.attrs["group_sizes"] = sizes
    LOG.info(f"consensus ({method}): {len(out):,} peaks")
    return out


def peak_overlap_matrix(reference: pd.DataFrame, peak_sets: dict) -> pd.DataFrame:
    """Boolean matrix: which reference peak is supported by which sample."""
    from .utils import overlap
    ref = as_intervals(reference).reset_index(drop=True)
    cols = {}
    for name, pk in peak_sets.items():
        hit = np.zeros(len(ref), dtype=bool)
        pr = overlap(ref, as_intervals(pk))
        if not pr.empty:
            hit[pr["idx_a"].astype(int).to_numpy()] = True
        cols[name] = hit
    return pd.DataFrame(cols, index=ref.index)


def saf_from_intervals(df: pd.DataFrame, path: str | None = None) -> pd.DataFrame:
    """SAF table for featureCounts (1-based, inclusive)."""
    d = as_intervals(df)
    saf = pd.DataFrame({
        "GeneID": d["chrom"].astype(str) + ":" + d["start"].astype(str) + "-" + d["end"].astype(str),
        "Chr": d["chrom"], "Start": d["start"] + 1, "End": d["end"], "Strand": "+",
    })
    if path:
        saf.to_csv(path, sep="\t", index=False)
        LOG.info(f"wrote SAF with {len(saf):,} features to {path}")
    return saf


def jaccard(a: pd.DataFrame, b: pd.DataFrame) -> float:
    """Jaccard index of two interval sets (by covered bases)."""
    from .utils import overlap
    A, B = merge_intervals(a), merge_intervals(b)
    inter = overlap(A, B)["ovl"].sum() if len(overlap(A, B)) else 0
    la = int((A["end"] - A["start"]).sum())
    lb = int((B["end"] - B["start"]).sum())
    union = la + lb - inter
    return float(inter / union) if union else 0.0
