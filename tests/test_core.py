"""Unit tests: synthetic data with a known ground truth.

The critical test is `test_audit_catches_reversed_sign` -- it reproduces the
exact bug (an R-style alphabetical factor contrast) and asserts the audit layer
reports FAIL.
"""
import numpy as np
import pandas as pd
import pytest

import epimux as ep
from epimux.utils import Contrast


# --------------------------------------------------------------------------
def make_data(n=800, n_rep=3, n_change=120, effect=2.0, seed=0):
    """Counts where the first `n_change` features are truly UP in KO."""
    rng = np.random.default_rng(seed)
    base = rng.lognormal(4.2, 0.9, n)
    samples, groups = [], {"WT": [], "KO": []}
    cols = {}
    for g in ("WT", "KO"):
        for r in range(1, n_rep + 1):
            s = f"{g}_R{r}"
            mu = base.copy()
            if g == "KO":
                mu[:n_change] *= effect
            cols[s] = rng.poisson(mu)
            groups[g].append(s)
            samples.append(s)
    counts = pd.DataFrame(cols, index=[f"f{i}" for i in range(n)])
    ivl = pd.DataFrame({
        "chrom": "chr1",
        "start": np.arange(n) * 10_000,
        "end": np.arange(n) * 10_000 + 1_000,
    }, index=counts.index)
    return ivl, counts, groups


def make_dataset(**kw):
    ivl, counts, groups = make_data(**kw)
    ref = ivl.reset_index(drop=True)
    ds = ep.Dataset(ref, genome="test", name="synthetic")
    ds.add_counts("ATAC", intervals=ivl, matrix=counts)
    ds.set_design(groups)
    return ds, counts, groups


# --------------------------------------------------------------------------
def test_interval_overlap():
    a = pd.DataFrame({"chrom": ["chr1", "chr1"], "start": [0, 500], "end": [100, 600]})
    b = pd.DataFrame({"chrom": ["chr1", "chr1"], "start": [50, 5000], "end": [150, 5100]})
    ov = ep.overlap(a, b)
    assert len(ov) == 1
    assert ov.iloc[0]["idx_a"] == 0 and ov.iloc[0]["idx_b"] == 0
    assert ov.iloc[0]["ovl"] == 50


def test_contrast_direction_is_explicit():
    c = Contrast(ref="WT", test="KO", group={"WT": ["a"], "KO": ["b"]})
    assert c.ref_samples == ["a"] and c.test_samples == ["b"]
    # alphabetical ordering must NOT influence the contrast
    c2 = Contrast(ref="KO", test="WT", group={"WT": ["a"], "KO": ["b"]})
    assert c2.test_samples == ["a"]


def test_bh_fdr_monotone():
    p = np.array([0.001, 0.01, 0.02, 0.5, np.nan])
    q = ep.bh_fdr(p)
    assert np.isnan(q[-1])
    assert np.all(np.diff(q[:4]) >= -1e-12)
    assert np.all(q[:4] >= p[:4] - 1e-12)


def test_differential_recovers_truth():
    ds, counts, groups = make_dataset()
    res = ds.differential(ref="WT", test="KO")["ATAC"]
    sig = (res["padj"] < 0.1) & (res["log2FC"] > np.log2(1.5))
    # the planted features should dominate the up set
    assert sig.iloc[:120].sum() > 60
    assert sig.iloc[120:].sum() < 30
    # and the direction must be UP (KO higher)
    assert res.loc[sig, "log2FC"].median() > 0


def test_audit_catches_reversed_sign():
    """The regression test for the real bug: a backwards contrast must FAIL."""
    ds, counts, groups = make_dataset()
    good = Contrast(ref="WT", test="KO", group=groups)
    bad = Contrast(ref="KO", test="WT", group=groups)   # what alphabetical levels gave us
    res_good = ep.deseq2_de(counts, good)
    ok = ep.check_direction(res_good, counts, good)
    assert ok.status == "pass", ok.summary
    # feed the correct result but claim the opposite contrast -> must be caught
    caught = ep.check_direction(res_good, counts, bad)
    assert caught.status == "fail"
    assert "REVERSED" in caught.summary


def test_null_contrast_is_empty():
    ds, counts, groups = make_dataset(n_rep=4)
    r = ep.null_contrast(counts, groups["WT"],
                         lambda c, ct: ep.deseq2_de(c, ct))
    assert r.status == "pass", r.summary


def test_positive_control_detects_planted_effect():
    ds, counts, groups = make_dataset(n_change=400, effect=3.0)
    ctr = Contrast(ref="WT", test="KO", group=groups)
    r = ep.positive_control(counts, ctr, lambda c, ct: ep.deseq2_de(c, ct))
    assert r.status == "pass", r.summary


def test_coupling_sign_and_orientation_guard():
    ds, counts, groups = make_dataset()
    ds.differential(ref="WT", test="KO")
    a = ds.results["ATAC"]
    # a perfectly concordant second layer
    b = a.copy()
    b.attrs.update(a.attrs)
    c = ep.couple(a, b, "A", "B")
    assert c.spearman > 0.9
    assert c.opposite_direction == 0
    # opposite orientation must raise instead of silently flipping the sign
    flipped = a.copy()
    flipped["log2FC"] = -flipped["log2FC"]
    flipped.attrs["contrast"] = repr(Contrast(ref="KO", test="WT", group=groups))
    with pytest.raises(ValueError, match="orientation mismatch"):
        ep.couple(a, flipped, "A", "B")


def test_efficiency_balance_direction_logic():
    ctr = Contrast(ref="WT", test="KO", group={"WT": ["w1", "w2"], "KO": ["k1", "k2"]})
    frip = {"w1": 0.42, "w2": 0.43, "k1": 0.33, "k2": 0.34}   # KO less efficient
    # observed effect is DOWN -> same direction as the bias -> dangerous
    bad = ep.efficiency_balance(frip, ctr, observed_direction="down")
    assert bad.status == "fail"
    # observed effect is UP -> against the bias -> conservative
    good = ep.efficiency_balance(frip, ctr, observed_direction="up")
    assert good.status == "pass"


def test_classify_states():
    ds, counts, groups = make_dataset()
    ds.differential(ref="WT", test="KO")
    r = ds.results["ATAC"]
    cls = ep.classify_elements({"ATAC": r, "H3K": r})
    assert "coordinated_up" in cls["state"].unique()
    assert (cls["state"] == "discordant").sum() == 0


def test_modules_partition():
    rng = np.random.default_rng(1)
    n = 600
    layers = {"a": pd.Series(rng.normal(size=n)), "b": pd.Series(rng.normal(size=n))}
    layers["a"][:200] += 4
    m = ep.find_modules(layers, k=3, seed=0)
    assert m.labels.nunique() == 3
    assert m.profile.shape == (3, 2)
