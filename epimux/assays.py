"""Assay adapters.

Each adapter owns its native features, knows how to load them, which statistical
engine is appropriate, and how to project itself onto the shared reference
elements.  Adding a new data type means adding one subclass.
"""
from __future__ import annotations

import glob
import os
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from . import stats as st
from .utils import (Contrast, as_intervals, get_logger, log2cpm,
                    map_to_reference, read_bed)

LOG = get_logger()

__all__ = ["Assay", "CountAssay", "MethylAssay", "SignalAssay", "HiCAssay",
           "read_featurecounts"]


# --------------------------------------------------------------------------
def read_featurecounts(path, strip=(".mLb.clN.sorted.bam", ".target.markdup.sorted.bam",
                                    ".sorted.bam", ".bam", "./")):
    """Read a featureCounts table into (intervals, counts).

    Sample names are cleaned of the usual nf-core BAM suffixes, which otherwise
    make design dictionaries unreadable.
    """
    df = pd.read_csv(path, sep="\t", comment="#")
    cols = list(df.columns)
    id_col = cols[0]
    ivl = df.rename(columns={"Chr": "chrom", "Start": "start", "End": "end"})
    ivl = as_intervals(ivl[["chrom", "start", "end"]])
    ivl.index = df[id_col].astype(str)
    sample_cols = [c for c in cols if c.lower().endswith(".bam")]
    if not sample_cols:
        sample_cols = cols[6:]
    counts = df[sample_cols].copy()
    clean = []
    for c in sample_cols:
        n = os.path.basename(str(c))
        for s in strip:
            n = n.replace(s, "")
        clean.append(n)
    counts.columns = clean
    counts.index = df[id_col].astype(str)
    return ivl, counts


# --------------------------------------------------------------------------
@dataclass
class Assay:
    """Base class: an assay is features + a matrix + a mapping to the reference."""
    name: str
    intervals: pd.DataFrame
    matrix: pd.DataFrame
    kind: str = "generic"
    meta: dict = field(default_factory=dict)
    _ref_map: pd.Series | None = field(default=None, repr=False)

    # -- projection onto the shared element set -----------------------------
    def map_to(self, reference: pd.DataFrame, how: str = "best",
               weight: str | None = None):
        feats = self.intervals.reset_index(drop=True)
        if weight is None and "baseMean" in self.matrix.columns:
            weight = None
        self._ref_map = map_to_reference(feats, reference, how="best", weight=weight)
        n = int((self._ref_map >= 0).sum())
        LOG.info(f"{self.name}: mapped {n:,}/{len(feats):,} features to reference elements")
        return self._ref_map

    def to_elements(self, values: pd.Series, reference: pd.DataFrame,
                    agg: str = "strongest", rank_by: pd.Series | None = None) -> pd.Series:
        """Collapse feature-level values onto reference elements.

        ``agg`` -- ``"strongest"`` keeps the feature with the smallest ``rank_by``
        (typically padj); ``"mean"`` averages; ``"max_abs"`` keeps the largest
        magnitude.
        """
        if self._ref_map is None:
            self.map_to(reference)
        m = self._ref_map
        ok = m >= 0
        df = pd.DataFrame({"ref": m[ok].to_numpy(),
                           "val": np.asarray(values)[ok.to_numpy()]})
        if agg == "mean":
            g = df.groupby("ref")["val"].mean()
        elif agg == "max_abs":
            df["a"] = df["val"].abs()
            g = df.sort_values("a", ascending=False).drop_duplicates("ref").set_index("ref")["val"]
        else:
            if rank_by is None:
                df["r"] = -df["val"].abs()
            else:
                df["r"] = np.asarray(rank_by)[ok.to_numpy()]
            g = df.sort_values("r").drop_duplicates("ref").set_index("ref")["val"]
        out = pd.Series(np.nan, index=np.arange(len(reference)))
        out.loc[g.index] = g.to_numpy()
        return out

    def differential(self, contrast: Contrast, **kw) -> pd.DataFrame:
        raise NotImplementedError

    @property
    def samples(self):
        return list(self.matrix.columns)


# --------------------------------------------------------------------------
@dataclass
class CountAssay(Assay):
    """Read-count assay (ATAC, ChIP-seq, CUT&RUN, RNA). Engine: PyDESeq2."""
    kind: str = "count"

    @classmethod
    def from_featurecounts(cls, name, path, **kw):
        ivl, cnt = read_featurecounts(path)
        summary = f"{path}.summary"
        meta = {}
        if os.path.exists(summary):
            meta["frip"] = frip_from_summary(summary)
        return cls(name=name, intervals=ivl, matrix=cnt, meta=meta, **kw)

    def differential(self, contrast: Contrast, **kw) -> pd.DataFrame:
        return st.deseq2_de(self.matrix, contrast, **kw)

    def norm(self) -> pd.DataFrame:
        return pd.DataFrame(log2cpm(self.matrix.to_numpy(dtype=float)),
                            index=self.matrix.index, columns=self.matrix.columns)


def frip_from_summary(path) -> dict:
    """Fraction of reads in peaks from a featureCounts .summary file."""
    s = pd.read_csv(path, sep="\t", index_col=0)
    s.columns = [os.path.basename(str(c)).replace(".mLb.clN.sorted.bam", "")
                 .replace(".target.markdup.sorted.bam", "").replace(".bam", "")
                 for c in s.columns]
    tot = s.sum(axis=0)
    return (s.loc["Assigned"] / tot.replace(0, np.nan)).to_dict()


# --------------------------------------------------------------------------
@dataclass
class MethylAssay(Assay):
    """Bisulfite methylation. Holds methylated and total counts per element."""
    kind: str = "methylation"
    cov: pd.DataFrame | None = None

    @classmethod
    def from_bismark_cov(cls, name, pattern, reference: pd.DataFrame,
                         sample_from=None):
        """Aggregate per-CpG Bismark .cov files onto reference elements.

        Bismark cov columns: chrom, start, end, meth%, n_meth, n_unmeth.
        """
        files = sorted(glob.glob(pattern)) if isinstance(pattern, str) else list(pattern)
        if not files:
            raise FileNotFoundError(f"no .cov files matched {pattern}")
        ref = as_intervals(reference).reset_index(drop=True)
        meth, cov = {}, {}
        for f in files:
            sample = (sample_from(f) if sample_from else
                      os.path.basename(f).split(".")[0])
            d = pd.read_csv(f, sep="\t", header=None,
                            names=["chrom", "start", "end", "pct", "m", "u"])
            d["chrom"] = d["chrom"].astype(str)
            d["start"] = d["start"] - 1
            d["end"] = d["start"] + 1
            d["tot"] = d["m"] + d["u"]
            mp = map_to_reference(d[["chrom", "start", "end"]], ref, how="best")
            ok = mp >= 0
            grp = pd.DataFrame({"ref": mp[ok].to_numpy(),
                                "m": d["m"].to_numpy()[ok.to_numpy()],
                                "t": d["tot"].to_numpy()[ok.to_numpy()]}).groupby("ref").sum()
            mm = pd.Series(0.0, index=np.arange(len(ref)))
            cc = pd.Series(0.0, index=np.arange(len(ref)))
            mm.loc[grp.index] = grp["m"].to_numpy()
            cc.loc[grp.index] = grp["t"].to_numpy()
            meth[sample], cov[sample] = mm, cc
            LOG.info(f"{name}: {sample} -> {len(d):,} CpGs aggregated")
        M, C = pd.DataFrame(meth), pd.DataFrame(cov)
        return cls(name=name, intervals=ref, matrix=M, cov=C)

    def differential(self, contrast: Contrast, **kw) -> pd.DataFrame:
        if self.cov is None:
            raise ValueError("MethylAssay needs coverage counts")
        return st.methylation_de(self.matrix, self.cov, contrast, **kw)

    def rates(self, min_cov: int = 10) -> pd.DataFrame:
        with np.errstate(invalid="ignore", divide="ignore"):
            r = self.matrix.to_numpy(float) / self.cov.to_numpy(float)
        r[self.cov.to_numpy(float) < min_cov] = np.nan
        return pd.DataFrame(r, index=self.matrix.index, columns=self.matrix.columns)


# --------------------------------------------------------------------------
@dataclass
class SignalAssay(Assay):
    """Continuous signal extracted from bigWig tracks.

    Provided for completeness and for legacy comparisons -- but the class warns
    on construction, because averaged track signal is markedly less sensitive
    than read counts and has historically produced both false negatives and
    sign errors in genotype contrasts.
    """
    kind: str = "signal"

    @classmethod
    def from_bigwigs(cls, name, files: dict, reference: pd.DataFrame,
                     stat: str = "mean", transform: str = "log1p"):
        try:
            import pyBigWig
        except ImportError as e:  # pragma: no cover
            raise ImportError("pyBigWig required for SignalAssay") from e
        LOG.warning(f"{name}: averaged bigWig signal is less sensitive than read counts; "
                    "prefer CountAssay for genotype contrasts and always run a positive control")
        ref = as_intervals(reference).reset_index(drop=True)
        out = {}
        for sample, path in files.items():
            bw = pyBigWig.open(str(path))
            chroms = bw.chroms()
            vals = np.full(len(ref), np.nan)
            for i, (c, s, e) in enumerate(zip(ref["chrom"], ref["start"], ref["end"])):
                if c not in chroms:
                    continue
                a, b = max(0, int(s)), min(int(e), chroms[c])
                if b <= a:
                    continue
                try:
                    v = bw.stats(c, a, b, type=stat)[0]
                except RuntimeError:
                    v = None
                if v is not None:
                    vals[i] = v
            bw.close()
            out[sample] = np.log1p(vals) if transform == "log1p" else vals
            LOG.info(f"{name}: extracted {sample}")
        return cls(name=name, intervals=ref, matrix=pd.DataFrame(out))

    def differential(self, contrast: Contrast, **kw) -> pd.DataFrame:
        return st.moderated_t_de(self.matrix, contrast, **kw)


# --------------------------------------------------------------------------
@dataclass
class HiCAssay(Assay):
    """Hi-C contact maps: compartments, insulation, local contact support.

    Hi-C needs no spike-in and is usually deeply sequenced, which makes it the
    right tool to ask about the *functional consequence* of a binding change
    when the corresponding ChIP is too shallow to quantify.
    """
    kind: str = "hic"
    coolers: dict = field(default_factory=dict)

    @classmethod
    def from_mcools(cls, name, files: dict, reference: pd.DataFrame):
        return cls(name=name, intervals=as_intervals(reference).reset_index(drop=True),
                   matrix=pd.DataFrame(index=range(len(reference))), coolers=dict(files))

    # -- compartments ------------------------------------------------------
    def eigenvector(self, sample: str, resolution: int = 160_000,
                    phasing_track: pd.DataFrame | None = None) -> pd.DataFrame:
        import bioframe
        import cooler
        import cooltools
        clr = cooler.Cooler(f"{self.coolers[sample]}::/resolutions/{resolution}")
        bins = clr.bins()[:][["chrom", "start", "end"]].copy()
        if phasing_track is not None:
            cov = bioframe.count_overlaps(bins, as_intervals(phasing_track)[["chrom", "start", "end"]])
            bins["phase"] = cov["count"].to_numpy(dtype=float)
        view = bioframe.make_viewframe(
            [(c, 0, clr.chromsizes[c]) for c in clr.chromnames
             if c in clr.chromsizes and str(c).startswith("chr") and c not in ("chrM", "chrY")])
        _, evec = cooltools.eigs_cis(clr, phasing_track=bins if phasing_track is not None else None,
                                     view_df=view, n_eigs=1)
        out = evec[["chrom", "start", "end", "E1"]].copy()
        out["chrom"] = out["chrom"].astype(str)
        return out.rename(columns={"E1": sample})

    def insulation(self, sample: str, resolution: int = 20_000,
                   window: int = 100_000) -> pd.DataFrame:
        import bioframe
        import cooler
        import cooltools
        clr = cooler.Cooler(f"{self.coolers[sample]}::/resolutions/{resolution}")
        view = bioframe.make_viewframe(
            [(c, 0, clr.chromsizes[c]) for c in clr.chromnames
             if c in clr.chromsizes and str(c).startswith("chr") and c not in ("chrM", "chrY")])
        t = cooltools.insulation(clr, [window], view_df=view, verbose=False)
        col = f"log2_insulation_score_{window}"
        out = t[["chrom", "start", "end", col]].copy()
        out["chrom"] = out["chrom"].astype(str)
        return out.rename(columns={col: sample})

    def local_contact(self, sample: str, targets: pd.DataFrame,
                      resolution: int = 20_000, flank: int = 200_000) -> pd.Series:
        """Summed balanced contacts within +/-flank of each target midpoint."""
        import cooler
        clr = cooler.Cooler(f"{self.coolers[sample]}::/resolutions/{resolution}")
        tg = as_intervals(targets).reset_index(drop=True)
        w = flank // resolution
        out = pd.Series(np.nan, index=np.arange(len(tg)))
        for chrom, grp in tg.groupby("chrom", sort=False):
            if chrom not in clr.chromnames:
                continue
            M = clr.matrix(balance=True).fetch(chrom)
            n = M.shape[0]
            mids = ((grp["start"] + grp["end"]) // 2 // resolution).to_numpy()
            for pos, b in zip(grp.index.to_numpy(), mids):
                if b < 0 or b >= n:
                    continue
                lo, hi = max(0, b - w), min(n, b + w + 1)
                out.loc[pos] = np.nansum(M[b, lo:hi])
        return out
