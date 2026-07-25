# FAQ and troubleshooting

### `check_direction` failed. What now?

Swap `ref` and `test` in your contrast — then re-check every figure and sentence
you have already written. The p-values were always right; only the direction was
wrong, which is why nothing else complained.

### My positive control finds almost nothing

Then no null result from that pipeline is interpretable. In order:

1. Are you using **read counts** rather than averaged bigWig signal?
2. Are the library sizes comparable to your other assays? A ChIP 30–100× shallower
   than your ATAC cannot be rescued by any statistical choice.
3. Do the samples cluster as expected in PCA? Mislabelling is common.

### `null_contrast` returns hits. Is my result dead?

Not necessarily. Compare scales: 100 null hits against 1,300 real hits means the
assay is noisier than you assumed and modest effects need care — it does not
invalidate a large, replicated effect. Run `outlier_replicates` next; a single
deviant library is the usual cause. Report the number either way.

### `couple()` raised "contrast orientation mismatch"

The two results were built from opposite ref/test assignments. Correlating them
would flip the sign and read as decoupling. Rebuild both with the same
orientation. This is deliberate — the error exists because that mistake produced
a published-grade "decoupling" claim that reversed on correction.

### Why is my methylation `log2FC` not a log ratio?

It is an absolute rate difference (`test - ref`). The column keeps the name for
schema consistency; `result.attrs["value_kind"] == "difference"` records it, and
thresholding and plotting honour it. A fold-change on a bounded rate is not
meaningful.

### `consensus_peaks` refuses to run

You asked for `method="replicated"` with unequal replicate numbers. That rule
favors the group with more replicates and manufactures apparent gains. Either
`balance="subsample"` (recommended) or `balance="ignore"` if you have a reason.

### Hi-C functions raise ImportError

Install the extras: `pip install -e ".[hic]"` (cooler, cooltools, bioframe).
bigWig needs `".[bigwig]"`, pathway enrichment `".[enrich]"`, everything
`".[all]"`.

### Sample names don't match my design

epimux strips common BAM suffixes from featureCounts columns. Check what it
produced:

```python
ds.assays["ATAC"].samples
```

and use exactly those strings in `set_design`.

### Figures use the wrong font

Set `EPIMUX_FONT` to a `.ttf`/`.ttc`, or call
`ep.plotting.set_style(font="/path/to/Helvetica.ttc")`. Without it a sane
sans-serif is used. PDFs embed Type 42 fonts so text stays editable in
Illustrator.

### How do I use spike-in normalization with `deseq2_de`?

Normalize first, then pass the adjusted matrix:

```python
sf = ep.spike_in_factors(spike_counts, target_lib)
adjusted = ep.apply_factors(counts, sf).round().astype(int)
ep.deseq2_de(adjusted, contrast)
```

### Is a "coordinated" or "discordant" call about a single element reliable?

Only if the element is significant in **both** layers. `classify_elements` calls
states from significance, never from the sign of a single noisy difference —
sign-only state calls are the classic source of irreproducible "discordant
element" lists that change membership when you shuffle replicates.

### Can I use this for single-cell data?

No. epimux is deliberately bulk. For single-cell multi-omics use `muon`,
`ArchR` or `Signac`.
