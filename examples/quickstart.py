"""Quickstart -- runs anywhere, no external data needed.

Generates a small synthetic two-assay dataset in which a known set of elements
is coordinately activated, then runs the full epimux workflow: differential,
audit, coupling, states, figures, HTML report.

    python examples/quickstart.py
"""
import os

import numpy as np
import pandas as pd

import epimux as ep

OUT = os.path.dirname(os.path.abspath(__file__))
rng = np.random.default_rng(0)

# ---------------------------------------------------------------- synthetic
N, N_CHANGE, N_REP = 4000, 300, 3
elements = pd.DataFrame({
    "chrom": "chr1",
    "start": np.arange(N) * 5_000,
    "end": np.arange(N) * 5_000 + 1_000,
})


def make_counts(effect, seed):
    """Three groups: WT, KO (first N_CHANGE elements up by `effect`), and OTHER,
    a different cell type used as the positive control.

    The cell-type shift is drawn ONCE per group so its replicates agree with each
    other. A positive control whose replicates disagree tests nothing.
    """
    r = np.random.default_rng(seed)
    base = r.lognormal(4.4, 0.9, N)
    cell_type_shift = r.lognormal(0, 0.8, N)      # drawn once, not per replicate
    cols, design = {}, {"WT": [], "KO": [], "OTHER": []}
    for g in ("WT", "KO", "OTHER"):
        mu = base.copy()
        if g == "KO":
            mu[:N_CHANGE] *= effect               # <-- the ground truth
        elif g == "OTHER":
            mu = mu * cell_type_shift
        for i in range(1, N_REP + 1):
            s = f"{g}_R{i}"
            cols[s] = r.poisson(mu)
            design[g].append(s)
    return pd.DataFrame(cols, index=[f"e{i}" for i in range(N)]), design


atac, design = make_counts(effect=2.2, seed=1)
h3k, _ = make_counts(effect=2.0, seed=2)      # same elements change -> coupled

# ------------------------------------------------------------------ epimux
ds = ep.Dataset(elements, genome="synthetic", name="quickstart")
ds.add_counts("ATAC", intervals=elements.set_index(atac.index), matrix=atac)
ds.add_counts("H3K27ac", intervals=elements.set_index(h3k.index), matrix=h3k)
ds.set_design(design)

# log2FC = log2(KO / WT) -- pinned by the contrast, verified against raw counts
ds.differential(ref="WT", test="KO")

# the audit battery. The positive control is a DIFFERENT comparison, one we know
# differs -- using the contrast of interest as its own control is circular.
ds.audit(positive_control=("WT", "OTHER"), null_group="WT")
print("\n" + repr(ds.audit_report))

print("\n" + repr(ds.coupling("ATAC", "H3K27ac")))

cls = ds.classify()
print("\nmulti-layer states:")
print(ep.concordance(cls).to_string())

# ----------------------------------------------------------------- figures
from epimux import plotting as pl

for name, r in ds.results.items():
    fig, _ = pl.ma_plot(r, title=f"{name}  KO vs WT")
    pl.save(fig, f"{OUT}/quickstart_{name}_MA", formats=("png",))
fig, _ = pl.coupling_plot(ds.results["ATAC"], ds.results["H3K27ac"], "ATAC", "H3K27ac")
pl.save(fig, f"{OUT}/quickstart_coupling", formats=("png",))

print("\nreport:", ds.report(f"{OUT}/quickstart_report.html"))
print("\nGround truth: elements 0-299 were made UP in KO in both assays.")
