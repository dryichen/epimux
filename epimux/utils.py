"""Interval primitives, logging and small helpers.

The whole package is anchored on genomic intervals: every assay, whatever its
native feature (ATAC peak, ChIP peak, CpG, Hi-C bin), is mapped onto one shared
reference element set.  These helpers implement that mapping with numpy only, so
the core works without bioframe/pyranges installed.
"""
from __future__ import annotations

import logging
import sys
from dataclasses import dataclass

import numpy as np
import pandas as pd

__all__ = [
    "get_logger", "BED_COLS", "read_bed", "as_intervals", "sort_intervals",
    "overlap", "map_to_reference", "midpoint_bin", "cpm", "log2cpm",
]

BED_COLS = ["chrom", "start", "end"]


def get_logger(name: str = "epimux") -> logging.Logger:
    log = logging.getLogger(name)
    if not log.handlers:
        h = logging.StreamHandler(sys.stderr)
        h.setFormatter(logging.Formatter("[epimux] %(message)s"))
        log.addHandler(h)
        log.setLevel(logging.INFO)
    return log


LOG = get_logger()


# --------------------------------------------------------------------------
# reading / normalising interval tables
# --------------------------------------------------------------------------
def read_bed(path, name_col: bool = True) -> pd.DataFrame:
    """Read a BED-like file into chrom/start/end[/name]."""
    df = pd.read_csv(path, sep="\t", header=None, comment="#")
    ncol = df.shape[1]
    cols = BED_COLS + ([f"c{i}" for i in range(3, ncol)])
    df.columns = cols
    if name_col and ncol > 3:
        df = df.rename(columns={"c3": "name"})
    return as_intervals(df)


def as_intervals(df: pd.DataFrame) -> pd.DataFrame:
    """Coerce a frame to canonical interval dtypes."""
    df = df.copy()
    lower = {c.lower(): c for c in df.columns}
    ren = {}
    for want, alts in (("chrom", ("chrom", "chr", "seqnames", "chromosome")),
                       ("start", ("start", "chromstart")),
                       ("end", ("end", "stop", "chromend"))):
        for a in alts:
            if a in lower and lower[a] != want:
                ren[lower[a]] = want
                break
    df = df.rename(columns=ren)
    missing = [c for c in BED_COLS if c not in df.columns]
    if missing:
        raise ValueError(f"missing interval columns: {missing}; got {list(df.columns)[:8]}")
    df["chrom"] = df["chrom"].astype(str)
    df["start"] = df["start"].astype(np.int64)
    df["end"] = df["end"].astype(np.int64)
    return df


def sort_intervals(df: pd.DataFrame) -> pd.DataFrame:
    return df.sort_values(["chrom", "start", "end"], kind="mergesort").reset_index(drop=True)


# --------------------------------------------------------------------------
# overlap engine (sweep per chromosome, numpy searchsorted)
# --------------------------------------------------------------------------
def overlap(a: pd.DataFrame, b: pd.DataFrame, min_overlap: int = 1) -> pd.DataFrame:
    """All overlapping pairs between two interval sets.

    Returns a frame with columns ``idx_a`` / ``idx_b`` (positional indices into
    the *input order* of ``a`` and ``b``) plus the overlap width ``ovl``.
    """
    a = as_intervals(a).reset_index(drop=True)
    b = as_intervals(b).reset_index(drop=True)
    out_a, out_b, out_o = [], [], []
    b_by_chrom = {c: g for c, g in b.groupby("chrom", sort=False)}
    for chrom, ga in a.groupby("chrom", sort=False):
        gb = b_by_chrom.get(chrom)
        if gb is None or len(gb) == 0:
            continue
        gb = gb.sort_values("start", kind="mergesort")
        bs = gb["start"].to_numpy()
        be = gb["end"].to_numpy()
        bidx = gb.index.to_numpy()
        # candidate window: b.start < a.end  and  b.end > a.start
        max_len = int((be - bs).max()) if len(be) else 0
        for ai, astart, aend in zip(ga.index.to_numpy(),
                                    ga["start"].to_numpy(),
                                    ga["end"].to_numpy()):
            lo = np.searchsorted(bs, astart - max_len, side="left")
            hi = np.searchsorted(bs, aend, side="left")
            if hi <= lo:
                continue
            cs, ce = bs[lo:hi], be[lo:hi]
            ov = np.minimum(ce, aend) - np.maximum(cs, astart)
            keep = ov >= min_overlap
            if not keep.any():
                continue
            k = np.nonzero(keep)[0]
            out_a.append(np.full(k.size, ai, dtype=np.int64))
            out_b.append(bidx[lo:hi][k])
            out_o.append(ov[k])
    if not out_a:
        return pd.DataFrame({"idx_a": [], "idx_b": [], "ovl": []}).astype(np.int64)
    return pd.DataFrame({
        "idx_a": np.concatenate(out_a),
        "idx_b": np.concatenate(out_b),
        "ovl": np.concatenate(out_o),
    })


def map_to_reference(features: pd.DataFrame, reference: pd.DataFrame,
                     how: str = "best", weight: str | None = None) -> pd.Series:
    """Map assay features onto reference elements.

    Parameters
    ----------
    how
        ``"best"``   keep, for each feature, the single reference element with
                     the largest overlap (a feature contributes once);
        ``"all"``    keep every overlapping pair.
    weight
        optional column in ``features`` used to break ties (largest wins) when
        ``how="best"`` and overlaps are equal.

    Returns
    -------
    Series indexed by feature position, values = reference positional index
    (``-1`` where the feature does not overlap the reference).
    """
    pairs = overlap(features, reference)
    out = pd.Series(np.full(len(features), -1, dtype=np.int64),
                    index=np.arange(len(features)), name="ref_idx")
    if pairs.empty:
        return out
    if how == "all":
        return pairs
    order = ["idx_a", "ovl"]
    asc = [True, False]
    if weight is not None and weight in features.columns:
        pairs = pairs.assign(_w=features[weight].to_numpy()[pairs["idx_a"].to_numpy()])
        order.append("_w")
        asc.append(False)
    best = pairs.sort_values(order, ascending=asc, kind="mergesort").drop_duplicates("idx_a")
    out.loc[best["idx_a"].to_numpy()] = best["idx_b"].to_numpy()
    return out


def midpoint_bin(df: pd.DataFrame, resolution: int) -> pd.Series:
    """chrom_binstart key at a fixed resolution (for Hi-C style binning)."""
    df = as_intervals(df)
    mid = ((df["start"] + df["end"]) // 2 // resolution) * resolution
    return df["chrom"].astype(str) + "_" + mid.astype(str)


# --------------------------------------------------------------------------
# normalisation helpers
# --------------------------------------------------------------------------
def cpm(counts: np.ndarray) -> np.ndarray:
    lib = counts.sum(axis=0, keepdims=True).astype(float)
    lib[lib == 0] = 1.0
    return counts / lib * 1e6


def log2cpm(counts: np.ndarray, prior: float = 1.0) -> np.ndarray:
    return np.log2(cpm(counts) + prior)


@dataclass
class Contrast:
    """A comparison, with the direction pinned down explicitly.

    ``log2FC`` produced anywhere in epimux is **always** ``log2(test / ref)``.
    Storing the contrast as an object (rather than relying on factor level
    ordering, which sorts alphabetically in R and bit us badly) is deliberate.
    """
    ref: str
    test: str
    group: dict  # {group_label: [sample names]}

    @property
    def ref_samples(self):
        return list(self.group[self.ref])

    @property
    def test_samples(self):
        return list(self.group[self.test])

    def design(self) -> pd.DataFrame:
        rows = [(s, self.ref) for s in self.ref_samples] + \
               [(s, self.test) for s in self.test_samples]
        return pd.DataFrame(rows, columns=["sample", "condition"]).set_index("sample")

    def __repr__(self):
        return (f"Contrast(log2FC = log2({self.test}/{self.ref}); "
                f"n_ref={len(self.ref_samples)}, n_test={len(self.test_samples)})")
