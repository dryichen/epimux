"""End-to-end example on the STAG2 LSK dataset.

Reproduces the corrected analysis in ~30 lines: two count assays over a shared
enhancer reference, one pinned contrast, the full audit battery, cross-layer
coupling, and an HTML report.
"""
import glob
import os

import pandas as pd

import epimux as ep

SP = os.environ.get("STAG2_DATA", "./data")   # set STAG2_DATA to your data directory
OUT = os.path.dirname(os.path.abspath(__file__))

# ---- 1. reference elements ------------------------------------------------
ds = ep.Dataset(f"{SP}/enh.bed", genome="mm10", name="STAG2 LSK")

# ---- 2. register assays ---------------------------------------------------
ds.add_counts("ATAC", f"{SP}/atac_rawcounts.txt")
ds.add_counts("H3K27ac", f"{SP}/h3k_counts.txt")

# ---- 3. design ------------------------------------------------------------
ds.set_design({
    "WT":  ["LSK_WT_REP1", "LSK_WT_REP2", "LSK_WT_REP3",
            "LSK_WT_R1", "LSK_WT_R2", "LSK_WT_R3"],
    "KO":  ["LSK_KO_REP1", "LSK_KO_REP2", "LSK_KO_REP3",
            "LSK_KO_R1", "LSK_KO_R2", "LSK_KO_R3"],
    "GMP": ["GMP_WT_REP1", "GMP_WT_REP2", "GMP_WT_R1", "GMP_WT_R2", "GMP_WT_R3"],
})

# ---- 4. differential: log2FC = log2(KO / WT), verified against raw counts --
res = ds.differential(ref="WT", test="KO")

# ---- 5. audit -------------------------------------------------------------
ds.audit(positive_control=("WT", "GMP"), null_group="WT")
print("\n" + repr(ds.audit_report))

# ---- 6. cross-layer coupling ---------------------------------------------
print("\n" + repr(ds.coupling("ATAC", "H3K27ac")))

# ---- 7. multi-layer states -----------------------------------------------
cls = ds.classify()
print("\nmulti-layer states:")
print(ep.concordance(cls).to_string())

# ---- 8. figures + report --------------------------------------------------
from epimux import plotting as pl

for name, r in ds.results.items():
    fig, _ = pl.ma_plot(r, title=f"{name}  KO vs WT (LSK)")
    pl.save(fig, f"{OUT}/fig_{name}_MA")
fig, _ = pl.coupling_plot(ds.results["ATAC"], ds.results["H3K27ac"], "ATAC", "H3K27ac")
pl.save(fig, f"{OUT}/fig_coupling")
fig, _ = pl.state_barplot(ep.concordance(cls))
pl.save(fig, f"{OUT}/fig_states")
fig, _ = pl.audit_plot(ds.audit_report)
pl.save(fig, f"{OUT}/fig_audit")

print("\nreport:", ds.report(f"{OUT}/stag2_lsk_report.html"))
