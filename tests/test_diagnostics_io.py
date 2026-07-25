"""Tests for diagnostics, peak handling and IO."""
import json
import os

import numpy as np
import pandas as pd
import pytest

import epimux as ep
from epimux.utils import Contrast


# ----------------------------------------------------------------- fixtures
@pytest.fixture
def signal_with_outlier():
    """Six samples; WT_R3 is deliberately deviant."""
    rng = np.random.default_rng(0)
    n = 800
    base = rng.normal(5, 2, n)
    cols = {}
    for g in ("WT", "KO"):
        for r in (1, 2, 3):
            cols[f"{g}_R{r}"] = base + rng.normal(0, 0.3, n)
    cols["WT_R3"] = rng.normal(5, 2, n)          # unrelated -> outlier
    return pd.DataFrame(cols), {"WT": ["WT_R1", "WT_R2", "WT_R3"],
                                "KO": ["KO_R1", "KO_R2", "KO_R3"]}


@pytest.fixture
def narrowpeaks(tmp_path):
    """Two groups, unequal replicate numbers (2 vs 3) — the biased case."""
    def write(path, starts):
        pd.DataFrame({"chrom": ["chr1"] * len(starts), "start": starts,
                      "end": [s + 500 for s in starts],
                      "name": [f"p{i}" for i in range(len(starts))],
                      "score": 100, "strand": ".", "signalValue": 5.0,
                      "pValue": 10.0, "qValue": 5.0, "peak": 250}
                     ).to_csv(path, sep="\t", header=False, index=False)
        return str(path)
    common = list(range(0, 50_000, 5_000))
    files = {
        "WT": [write(tmp_path / "wt1.narrowPeak", common),
               write(tmp_path / "wt2.narrowPeak", common)],
        "KO": [write(tmp_path / "ko1.narrowPeak", common + [200_000]),
               write(tmp_path / "ko2.narrowPeak", common + [200_000]),
               write(tmp_path / "ko3.narrowPeak", common + [200_000])],
    }
    return files


# -------------------------------------------------------------- diagnostics
def test_outlier_replicates_finds_the_bad_sample(signal_with_outlier):
    vals, groups = signal_with_outlier
    r = ep.outlier_replicates(vals, groups)
    assert r.status == "warn"
    assert "WT_R3" in r.summary


def test_outlier_replicates_clean_data_passes():
    rng = np.random.default_rng(1)
    base = rng.normal(0, 1, 500)
    vals = pd.DataFrame({f"WT_R{i}": base + rng.normal(0, .1, 500) for i in (1, 2, 3)})
    r = ep.outlier_replicates(vals, {"WT": ["WT_R1", "WT_R2", "WT_R3"]})
    assert r.status == "pass"


def test_pvalue_diagnostic_uniform_and_signal():
    rng = np.random.default_rng(0)
    # pure null -> uniform
    null = pd.DataFrame({"pvalue": rng.uniform(0, 1, 5000)})
    assert ep.pvalue_diagnostic(null).status == "pass"
    # signal -> spike at zero
    p = np.concatenate([rng.uniform(0, 0.001, 500), rng.uniform(0, 1, 4500)])
    sig = pd.DataFrame({"pvalue": p})
    r = ep.pvalue_diagnostic(sig)
    assert r.status == "pass" and "spike" in r.summary


def test_pvalue_diagnostic_flags_misspecification():
    rng = np.random.default_rng(0)
    # pile-up near 1 -> variance over-estimated
    p = np.concatenate([rng.uniform(0.9, 1.0, 3000), rng.uniform(0, 1, 2000)])
    r = ep.pvalue_diagnostic(pd.DataFrame({"pvalue": p}))
    assert r.status == "fail"


def test_confounding_check_detects_perfect_confounding():
    ctr = Contrast(ref="WT", test="KO",
                   group={"WT": ["w1", "w2"], "KO": ["k1", "k2"]})
    bad = pd.DataFrame({"batch": ["b1", "b1", "b2", "b2"]},
                       index=["w1", "w2", "k1", "k2"])
    assert ep.confounding_check(ctr, bad).status == "fail"
    good = pd.DataFrame({"batch": ["b1", "b2", "b1", "b2"]},
                        index=["w1", "w2", "k1", "k2"])
    assert ep.confounding_check(ctr, good).status == "pass"


def test_power_analysis_monotone():
    rng = np.random.default_rng(0)
    n = 2000
    lfc = rng.normal(0, 0.5, n)
    res = pd.DataFrame({"log2FC": lfc, "stat": lfc / 0.2,
                        "pvalue": rng.uniform(0, 1, n), "padj": rng.uniform(0, 1, n)})
    tab = ep.power_analysis(res, n_per_group=3)
    assert (np.diff(tab["power"]) >= -1e-9).all()          # power rises with n
    assert (np.diff(tab["detectable_log2FC_at_80pct"]) <= 1e-9).all()


# -------------------------------------------------------------------- peaks
def test_consensus_union(narrowpeaks):
    out = ep.consensus_peaks(narrowpeaks, method="union")
    assert len(out) == 11                       # 10 common + 1 KO-specific
    assert out.attrs["method"] == "union"


def test_consensus_replicated_refuses_unequal_replicates(narrowpeaks):
    """The default must refuse the setup that manufactures apparent gains."""
    with pytest.raises(ValueError, match="unequal replicate counts"):
        ep.consensus_peaks(narrowpeaks, method="replicated", min_replicates=2)


def test_consensus_replicated_subsample(narrowpeaks):
    out = ep.consensus_peaks(narrowpeaks, method="replicated",
                             min_replicates=2, balance="subsample")
    assert len(out) >= 10


def test_merge_and_jaccard():
    a = pd.DataFrame({"chrom": ["chr1", "chr1"], "start": [0, 400], "end": [500, 900]})
    m = ep.merge_intervals(a)
    assert len(m) == 1 and m.iloc[0]["end"] == 900
    b = pd.DataFrame({"chrom": ["chr1"], "start": [0], "end": [900]})
    assert ep.jaccard(a, b) == pytest.approx(1.0, abs=1e-6)


def test_saf_is_one_based(tmp_path):
    d = pd.DataFrame({"chrom": ["chr1"], "start": [100], "end": [200]})
    saf = ep.saf_from_intervals(d, str(tmp_path / "x.saf"))
    assert saf.iloc[0]["Start"] == 101 and saf.iloc[0]["End"] == 200
    assert os.path.exists(tmp_path / "x.saf")


def test_peak_overlap_matrix():
    ref = pd.DataFrame({"chrom": ["chr1", "chr1"], "start": [0, 10_000],
                        "end": [500, 10_500]})
    sets = {"s1": pd.DataFrame({"chrom": ["chr1"], "start": [100], "end": [200]})}
    M = ep.peak_overlap_matrix(ref, sets)
    assert bool(M.iloc[0, 0]) and not bool(M.iloc[1, 0])


# ----------------------------------------------------------------------- io
def _small_dataset():
    rng = np.random.default_rng(0)
    n = 300
    base = rng.lognormal(4, .8, n)
    cols, design = {}, {"WT": [], "KO": []}
    for g in ("WT", "KO"):
        for r in (1, 2, 3):
            mu = base.copy()
            if g == "KO":
                mu[:50] *= 3
            cols[f"{g}_R{r}"] = rng.poisson(mu)
            design[g].append(f"{g}_R{r}")
    counts = pd.DataFrame(cols, index=[f"e{i}" for i in range(n)])
    el = pd.DataFrame({"chrom": "chr1", "start": np.arange(n) * 5000,
                       "end": np.arange(n) * 5000 + 1000})
    ds = ep.Dataset(el, genome="test", name="io-test")
    ds.add_counts("ATAC", intervals=el.set_index(counts.index), matrix=counts)
    ds.set_design(design)
    ds.differential(ref="WT", test="KO")
    return ds


def test_export_results_writes_manifest(tmp_path):
    ds = _small_dataset()
    written = ep.export_results(ds, str(tmp_path))
    assert os.path.exists(written["manifest"])
    man = json.load(open(written["manifest"]))
    assert man["sign_convention"].startswith("log2FC = log2(test / ref)")
    assert man["results"]["ATAC"]["up"] > 0
    assert man["results"]["ATAC"]["up"] > man["results"]["ATAC"]["down"]


def test_export_bed_direction(tmp_path):
    ds = _small_dataset()
    p = ep.export_bed(ds, "ATAC", str(tmp_path / "up.bed"), direction="up")
    bed = pd.read_csv(p, sep="\t", header=None)
    assert (bed[6] > 0).all()                    # log2FC column
    assert bed.shape[1] == 8


def test_save_and_load_roundtrip(tmp_path):
    ds = _small_dataset()
    p = ep.save_dataset(ds, str(tmp_path / "ds.json"))
    back = ep.load_results(p)
    assert "ATAC" in back
    assert len(back["ATAC"]) == len(ds.results["ATAC"])


def test_to_anndata_carries_contrast_and_audit():
    ds = _small_dataset()
    ad = ep.to_anndata(ds, "ATAC")
    assert ad.n_obs == 6
    assert "log2(KO/WT)" in ad.uns["contrast"]
    assert "epimux_version" in ad.uns


# ------------------------------------------------------------- cross-contrast
def _paired_results(seed=0, effect_b=1.0, n=1000, n_true=100, se_val=0.2):
    """Two contrasts over the same elements; `effect_b` scales the second one.

    The first `n_true` elements carry a real effect; the rest are null, so any
    correlation computed over ALL elements is diluted by design.
    """
    from scipy import stats as sstats
    rng = np.random.default_rng(seed)
    true = np.concatenate([rng.normal(1.5, .3, n_true), np.zeros(n - n_true)])

    def mk(scale):
        lfc = true * scale + rng.normal(0, .15, n)
        se = np.full(n, se_val)
        stat = lfc / se
        p = 2 * sstats.norm.sf(np.abs(stat))
        return pd.DataFrame({"baseMean": 100.0, "log2FC": lfc, "lfcSE": se,
                             "stat": stat, "pvalue": p, "padj": ep.bh_fdr(p)})
    return mk(1.0), mk(effect_b)


def test_compare_contrasts_finds_no_interaction_when_equal():
    a, b = _paired_results(effect_b=1.0)
    cmp = ep.compare_contrasts(a, b, "LSK", "GMP")
    assert (cmp["padj_interaction"] < 0.1).sum() < 30       # ~nothing differs


def test_compare_contrasts_detects_real_interaction():
    a, b = _paired_results(effect_b=0.0)                    # effect absent in B
    cmp = ep.compare_contrasts(a, b, "LSK", "GMP")
    assert (cmp["padj_interaction"] < 0.1).sum() > 50
    # and the interactions must sit on the elements that actually carry an effect
    top = cmp.nsmallest(50, "padj_interaction").index
    assert (top < 100).mean() > 0.8


def test_concordance_summary_separates_power_from_biology():
    """Identical effects: hits should be shared, with few real interactions."""
    a, b = _paired_results(effect_b=1.0)
    cmp = ep.compare_contrasts(a, b, "LSK", "GMP")
    s = ep.concordance_summary(cmp, "LSK", "GMP")
    # correlation over the elements that carry an effect (not the null 90%)
    has_effect = cmp.index < 100
    from scipy import stats as sstats
    r = sstats.spearmanr(cmp.loc[has_effect, "log2FC_LSK"], cmp.loc[has_effect, "log2FC_GMP"])[0]
    assert r > 0.8
    assert s["shared"] > 0
    assert s["shared_same_direction"] == s["shared"]
    assert s["elements_with_interaction"] < 30


def test_concordance_summary_flags_real_difference():
    a, b = _paired_results(effect_b=0.0)
    s = ep.concordance_summary(ep.compare_contrasts(a, b, "LSK", "GMP"), "LSK", "GMP")
    assert s["LSK_only"] > 0
    assert s["fraction_of_differences_supported"] > 0.5
    assert "genuinely differ" in s["interpretation"]


def test_meta_analyse_shrinks_standard_error():
    a, b = _paired_results(effect_b=1.0, se_val=0.2)
    m = ep.meta_analyze({"one": a, "two": b})
    # inverse-variance pooling of two SE=0.2 estimates -> 0.2/sqrt(2) ~ 0.141
    assert m["se"].median() == pytest.approx(0.2 / np.sqrt(2), rel=0.05)
    assert "pvalue_heterogeneity" in m


def test_meta_analyse_detects_heterogeneity():
    a, b = _paired_results(effect_b=0.0)
    m = ep.meta_analyze({"one": a, "two": b})
    het = (m["pvalue_heterogeneity"] < 0.05)
    assert het.iloc[:100].mean() > het.iloc[100:].mean()   # true-effect elements disagree


def test_replication_rate_direction():
    a, b = _paired_results(effect_b=1.0)
    r = ep.replication_rate(a, b)
    assert r["n_discovery_hits"] > 0
    assert r["same_direction_rate"] > 0.9
    assert r["replication_rate_loose"] > 0.8
