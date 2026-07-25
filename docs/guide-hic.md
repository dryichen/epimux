# Hi-C analyses

Contact maps are usually the deepest and least normalization-fragile assay in a
multi-omic study, and they need **no spike-in**. That makes Hi-C the right tool
when a ChIP is too shallow to quantify but the biological question still stands.

```python
ds.add_hic("HiC", {"WT": "wt.mcool", "KO": "ko.mcool"})
hic = ds.assays["HiC"]
```

---

## Compartments

```python
ev_wt = hic.eigenvector("WT", resolution=160_000, phasing_track=active_peaks)
ev_ko = hic.eigenvector("KO", resolution=160_000, phasing_track=active_peaks)

m = ev_wt.merge(ev_ko, on=["chrom", "start", "end"]).dropna()
m["dE1"] = m["KO"] - m["WT"]
print("identity agreement:", ((m["WT"] > 0) == (m["KO"] > 0)).mean())
```

`phasing_track` orients the sign so positive = A. Peak density of an active mark
works well; GC content is the standard alternative.

**Identity and strength are different questions.** Compartment *identity* can be
completely preserved (r ≈ 1.0) while compartment *strength* changes a lot —
reporting only the eigenvector correlation hides the effect.

```python
S, _ = ep.hic.saddle("wt.mcool", ev_wt, "WT")
strength = ep.hic.compartment_strength(S)      # (AA + BB) / 2AB
```

Compute per chromosome and compare paired — a change on 19/19 chromosomes is far
more convincing than a genome-wide average with no error bar.

---

## Contact decay and the extrusion shoulder

```python
decay = ep.hic.contact_decay("wt.mcool", resolution=20_000)
shoulder = ep.hic.extrusion_shoulder(decay, s_min=50_000, s_max=5_000_000)
```

Always compare the **log-derivative**, not the raw curves. On a log-log plot a
large difference in extrusion looks like two nearly identical lines; the
derivative makes the shoulder — and its change — visible.

---

## Insulation and boundaries

```python
ins_wt = hic.insulation("WT", resolution=20_000, window=100_000)
ins_ko = hic.insulation("KO", resolution=20_000, window=100_000)

d = ep.hic.differential_insulation(ins_wt, ins_ko, "WT", "KO", anchors=ctcf_peaks)
strong = ep.hic.boundary_strength(ins_wt, "WT", quantile=0.05)
```

**Sign convention:** the insulation score is a log ratio, so **lower = stronger
boundary**. A positive `dIns` therefore means a *weaker* boundary in the test
condition. Label it explicitly in figures — this is the most common
misinterpretation.

Restricting to the strongest boundaries (`quantile=0.05`) is usually where a
real effect is detectable; genome-wide medians average it away.

---

## Pileups and APA

```python
mat = ep.hic.pileup("wt.mcool", boundaries, flank=200_000)      # on-diagonal
mat, score = ep.hic.apa("wt.mcool", loops, flank=100_000)       # off-diagonal
```

`pileup` divides by the distance-dependent expected by default, so the result is
observed/expected rather than raw contact. `apa` returns a center/corner ratio —
compare that number between conditions on the **same loop set**, never on
per-condition loop calls, or you are measuring the calling threshold rather than
the biology.

---

## Linking elements to genes by contact

```python
links = ds.link_genes(tss, method="abc", hic="HiC", sample="WT",
                      activity=element_activity)
gene_effect = ep.aggregate_to_genes(links, ds.results["ATAC"]["log2FC"])
```

Nearest-gene assignment (`method="nearest"`) is available for comparison and is
clearly labeled as not recommended: enhancers routinely skip their neighbour,
often enough to change the conclusion of a pathway analysis.
