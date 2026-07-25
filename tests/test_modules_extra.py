"""Tests for normalization, annotation, enrichment, linking and CLI."""
import numpy as np
import pandas as pd
import pytest

import epimux as ep
from epimux.utils import Contrast


# ------------------------------------------------------------------ fixtures
@pytest.fixture
def counts():
    rng = np.random.default_rng(0)
    n = 500
    base = rng.lognormal(4.0, 0.8, n)
    cols = {}
    for g in ("WT", "KO"):
        for r in (1, 2, 3):
            depth = 1.0 if g == "WT" else 2.0          # KO sequenced twice as deep
            cols[f"{g}_R{r}"] = rng.poisson(base * depth)
    return pd.DataFrame(cols, index=[f"f{i}" for i in range(n)])


@pytest.fixture
def tss():
    return pd.DataFrame({
        "chrom": ["chr1"] * 5,
        "start": [10_000, 100_000, 500_000, 1_000_000, 2_000_000],
        "end": [10_001, 100_001, 500_001, 1_000_001, 2_000_001],
        "gene": ["Gene1", "Gene2", "Gene3", "Gene4", "Gene5"],
        "strand": ["+", "-", "+", "+", "-"],
    })


# --------------------------------------------------------------- normalization
def test_size_factors_correct_depth(counts):
    sf = ep.median_of_ratios(counts)
    ko = [s for s in counts.columns if s.startswith("KO")]
    wt = [s for s in counts.columns if s.startswith("WT")]
    assert sf[ko].mean() > sf[wt].mean() * 1.5      # KO is ~2x deeper
    norm = ep.apply_factors(counts, sf)
    assert abs(norm[ko].sum().mean() / norm[wt].sum().mean() - 1) < 0.15


def test_tmm_runs_and_is_centred(counts):
    sf = ep.tmm(counts)
    assert len(sf) == counts.shape[1]
    assert abs(np.exp(np.mean(np.log(sf))) - 1) < 1e-6


def test_quantile_normalize_equalises_distributions(counts):
    q = ep.quantile_normalize(np.log2(counts + 1))
    means = q.mean()
    assert means.max() - means.min() < 1e-8


def test_spike_in_factors_survive_global_shift():
    """A genome-wide gain must NOT be absorbed by spike-in factors."""
    spike = {"WT_R1": 1000, "WT_R2": 1100, "KO_R1": 1050, "KO_R2": 980}
    f = ep.spike_in_factors(spike)
    assert abs(np.exp(np.mean(np.log(f))) - 1) < 1e-9
    # equal spike-in -> equal factors, regardless of target signal
    assert f.max() / f.min() < 1.2


def test_assess_global_shift_flags_consistent_offset():
    rng = np.random.default_rng(1)
    n = 2000
    base = rng.lognormal(4, 0.5, n)
    df = pd.DataFrame({
        "WT_R1": rng.poisson(base), "WT_R2": rng.poisson(base),
        "KO_R1": rng.poisson(base * 2), "KO_R2": rng.poisson(base * 2),
    }, index=[f"f{i}" for i in range(n)])
    ctr = Contrast(ref="WT", test="KO",
                   group={"WT": ["WT_R1", "WT_R2"], "KO": ["KO_R1", "KO_R2"]})
    # depth-normalized, a uniform 2x gain looks like NO shift (that is the trap)
    out = ep.assess_global_shift(df, ctr)
    assert "consistent_shift" in out and "interpretation" in out


def test_reference_normalize_requires_enough_features(counts):
    with pytest.raises(ValueError):
        ep.reference_normalize(counts, counts.index[:5])
    sf = ep.reference_normalize(counts, counts.index[:100])
    assert len(sf) == counts.shape[1]


# ----------------------------------------------------------------- annotation
def test_annotate_tss_signed_distance(tss):
    el = pd.DataFrame({"chrom": ["chr1", "chr1"], "start": [10_500, 99_000],
                       "end": [10_600, 99_100]})
    ann = ep.annotate_tss(el, tss)
    assert ann.loc[0, "nearest_gene"] == "Gene1"
    assert ann.loc[0, "tss_distance"] > 0                 # downstream of a + gene
    assert ann.loc[1, "nearest_gene"] == "Gene2"
    assert ann.loc[1, "tss_distance"] > 0                 # upstream of a - gene -> sign flipped


def test_classify_context(tss):
    el = pd.DataFrame({"chrom": ["chr1"] * 3, "start": [10_100, 15_000, 300_000],
                       "end": [10_200, 15_100, 300_100]})
    ann = ep.classify_context(el, tss)
    assert list(ann["context"]) == ["promoter", "proximal", "distal"]


def test_super_enhancer_stitching():
    """100 well-separated peaks; a handful carry most of the signal."""
    rng = np.random.default_rng(0)
    n = 100
    starts = np.arange(n) * 200_000            # 200kb apart -> never stitched together
    pk = pd.DataFrame({"chrom": ["chr1"] * n, "start": starts, "end": starts + 1_000})
    sig = np.concatenate([rng.uniform(1, 3, n - 5), rng.uniform(80, 120, 5)])
    se = ep.stitch_super_enhancers(pk, sig, stitch_distance=12_500)
    assert len(se) == n                        # nothing stitched at this spacing
    assert se["super_enhancer"].sum() >= 1     # the strong tail is called
    assert se["super_enhancer"].sum() < n // 2  # but not everything
    assert se.iloc[0]["signal"] >= se.iloc[-1]["signal"]   # ranked descending
    # the called set must be the top of the ranking
    assert se.loc[se["super_enhancer"], "signal"].min() >= se.loc[~se["super_enhancer"], "signal"].max()


def test_super_enhancer_stitches_adjacent_peaks():
    """Peaks closer than the stitch distance become one region."""
    starts = np.arange(0, 40_000, 2_000)       # 2kb apart, well under 12.5kb
    pk = pd.DataFrame({"chrom": ["chr1"] * len(starts),
                       "start": starts, "end": starts + 500})
    se = ep.stitch_super_enhancers(pk, np.ones(len(starts)), stitch_distance=12_500)
    assert len(se) == 1
    assert se.iloc[0]["n_peaks"] == len(starts)


def test_context_composition_enrichment(tss):
    el = pd.DataFrame({"chrom": ["chr1"] * 4, "start": [10_100, 10_200, 300_000, 400_000],
                       "end": [10_200, 10_300, 300_100, 400_100]})
    ann = ep.classify_context(el, tss)
    comp = ep.context_composition(ann, selected=pd.Index([0, 1]))
    assert comp.loc["promoter", "enrichment"] > 1


# ----------------------------------------------------------------- enrichment
def test_overlap_enrichment_detects_real_overlap():
    uni = pd.DataFrame({"chrom": ["chr1"] * 200,
                        "start": np.arange(200) * 10_000,
                        "end": np.arange(200) * 10_000 + 1_000})
    target = uni.iloc[:50]
    query = uni.iloc[:20]                     # fully inside target
    r = ep.overlap_enrichment(query, target, uni, n_shuffle=200)
    assert r["observed"] == 20
    assert r["fold_enrichment"] > 2
    assert r["pvalue"] < 0.05


def test_gc_matched_background():
    el = pd.DataFrame({"chrom": ["chr1"] * 100, "start": np.arange(100) * 1000,
                       "end": np.arange(100) * 1000 + 500})
    gc = pd.Series(np.linspace(0.3, 0.7, 100))
    bg = ep.gc_matched_background(el, pd.Index([10, 20, 30]), gc, n_per=3)
    assert len(bg) > 0
    assert not set(bg) & {10, 20, 30}


def test_motif_enrichment_direction():
    fg = pd.DataFrame({"M1": [1] * 80 + [0] * 20, "M2": [0] * 100})
    bg = pd.DataFrame({"M1": [1] * 10 + [0] * 90, "M2": [0] * 100})
    r = ep.motif_enrichment(fg, bg).set_index("motif")
    assert r.loc["M1", "odds_ratio"] > 5
    assert r.loc["M1", "padj"] < 0.01


# -------------------------------------------------------------------- linking
def test_nearest_gene(tss):
    el = pd.DataFrame({"chrom": ["chr1"], "start": [95_000], "end": [95_100]})
    link = ep.nearest_gene(el, tss)
    assert link.iloc[0]["gene"] == "Gene2"


def test_aggregate_to_genes():
    links = pd.DataFrame({"element": [0, 0, 1], "gene": ["A", "B", "A"],
                          "score": [0.8, 0.2, 0.5]})
    vals = pd.Series([2.0, 4.0])
    out = ep.aggregate_to_genes(links, vals, how="max_abs")
    assert out["A"] == 4.0 and out["B"] == 2.0


# ------------------------------------------------------------------------ CLI
def test_cli_version(capsys):
    from epimux.cli import main
    assert main(["--version"]) == 0
    assert "epimux" in capsys.readouterr().out


def test_covariate_design_rejects_missing(counts):
    ctr = Contrast(ref="WT", test="KO",
                   group={"WT": ["WT_R1", "WT_R2", "WT_R3"],
                          "KO": ["KO_R1", "KO_R2", "KO_R3"]})
    cov = pd.DataFrame({"batch": ["a", "b"]}, index=["WT_R1", "WT_R2"])
    with pytest.raises(ValueError, match="covariates must cover"):
        ep.deseq2_de(counts, ctr, covariates=cov)


# ------------------------------------------------------------------- fonts
def test_bold_and_italic_actually_render():
    """A .ttc exposes only face 0 to FontProperties(fname=...), which silently
    drops bold and italic. set_style must register every face so emphasis works."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from epimux import plotting as pl

    fp = pl.set_style()

    def render(**kw):
        fig = plt.figure(figsize=(2, 0.5))
        fig.text(0.05, 0.3, "ATAC 572", fontproperties=fp, fontsize=20, **kw)
        fig.canvas.draw()
        arr = np.asarray(fig.canvas.buffer_rgba()).copy()
        plt.close(fig)
        return arr

    plain = render()
    assert not np.array_equal(plain, render(fontweight="bold")), \
        "fontweight='bold' renders identically to normal -- font faces not registered"
    assert not np.array_equal(plain, render(style="italic")), \
        "style='italic' renders identically to normal -- font faces not registered"
