# Input formats

Everything in epimux is anchored on one **reference element set**. Each assay
keeps its own native features and is projected onto that reference, so you never
repeat interval arithmetic and every layer becomes directly comparable.

---

## Reference elements

A BED file (or any frame with `chrom`/`start`/`end`) defining the units you want
to reason about — enhancers, promoters, an atlas, or a consensus peak set.

```
chr1	4417942	4418246	enh0
chr1	4622164	4624427	enh1
```

```python
ds = ep.Dataset("enhancers.bed", genome="mm10")
ds = ep.Dataset(dataframe, genome="mm10")     # or pass a frame directly
```

Column names are normalised, so `chr`/`seqnames`/`chromosome` and
`stop`/`chromEnd` are all accepted. Coordinates are **0-based, half-open** (BED
convention). The only place this changes is `saf_from_intervals`, which converts
to featureCounts' 1-based inclusive convention for you.

If you do not have a reference yet, build one:

```python
peaks = ep.consensus_peaks({"WT": [...], "KO": [...]}, method="union")
ep.saf_from_intervals(peaks, "consensus.saf")     # then run featureCounts
```

---

## Count assays — ATAC, ChIP-seq, CUT&RUN, RNA

**featureCounts output**, directly:

```python
ds.add_counts("ATAC", "atac_counts.txt")
```

Expected layout (the standard featureCounts table):

```
# Program:featureCounts ...
Geneid	Chr	Start	End	Strand	Length	./WT_R1.mLb.clN.sorted.bam	...
peak_1	chr1	4417943	4418246	+	304	152	...
```

Sample names are cleaned of the usual nf-core suffixes
(`.mLb.clN.sorted.bam`, `.target.markdup.sorted.bam`, …) so your design
dictionary stays readable. Check what epimux ended up with:

```python
ds.assays["ATAC"].samples
```

If the file has a `.summary` companion (featureCounts writes one), epimux reads
it and derives **FRiP** automatically, which feeds `efficiency_balance`.

**Already have a matrix?**

```python
ds.add_counts("ATAC", intervals=peak_frame, matrix=count_dataframe)
```

`intervals` and `matrix` must share an index, one row per feature, raw integer
counts (not CPM, not logged — the negative binomial model needs counts).

---

## Methylation — Bismark

```python
ds.add_methyl("WGBS", "wgbs/*.cov")
```

Bismark `.cov` columns, no header:

```
chrom  start(1-based)  end  methylation%  count_methylated  count_unmethylated
chr1   3000827         3000827  100       2                 0
```

Per-CpG records are aggregated onto reference elements (methylated and total
counts summed), so coverage is preserved and the test can weight by it. Sample
names default to the part of the filename before the first dot; override with
`sample_from=lambda path: ...`.

Methylation results use an **absolute rate difference**, not a log ratio. The
column is still called `log2FC` for schema consistency, but
`result.attrs["value_kind"] == "difference"` records this and the thresholding
and plotting code honours it.

---

## Continuous signal — bigWig

```python
ds.add_signal("H3K27ac_bw", {"WT_R1": "wt1.bw", "KO_R1": "ko1.bw"})
```

Supported, but epimux warns when you use it, and you should listen: **averaged
track signal is markedly less sensitive than read counts**. A "no change" result
from bigWig signal is very often a false negative of the method rather than a
biological null. If BAMs exist, count them and use `add_counts`.

If you must use signal, always pair it with a positive control
(`ds.audit(positive_control=...)`) before reporting any null.

---

## Hi-C — cooler

```python
ds.add_hic("HiC", {"WT": "wt.mcool", "KO": "ko.mcool"})
```

Multi-resolution `.mcool` files. Resolution is chosen per call:

```python
hic = ds.assays["HiC"]
hic.insulation("WT", resolution=20_000, window=100_000)
hic.eigenvector("WT", resolution=160_000, phasing_track=peaks)
```

`phasing_track` orients the eigenvector so that positive = A compartment. Peak
density of an active mark works well; GC content is the usual alternative.

---

## Design and contrast

```python
ds.set_design({
    "WT":  ["LSK_WT_R1", "LSK_WT_R2", "LSK_WT_R3"],
    "KO":  ["LSK_KO_R1", "LSK_KO_R2", "LSK_KO_R3"],
    "GMP": ["GMP_WT_R1", "GMP_WT_R2", "GMP_WT_R3"],   # for the positive control
})
ds.differential(ref="WT", test="KO")
```

A group may list samples that only exist in some assays — epimux intersects per
assay and skips any assay with no samples on one side of the contrast.

Covariates, when you have batch structure:

```python
cov = pd.DataFrame({"batch": ["b1", "b2", "b1", "b2", "b1", "b2"]},
                   index=[...samples...])
ds.differential(ref="WT", test="KO", covariates=cov)   # design becomes ~batch + condition
ds.audit(covariates=cov)                               # also checks confounding
```

---

## What comes back

```python
res = ds.differential(ref="WT", test="KO")
res["ATAC"].head()
```

| column | meaning |
|---|---|
| `baseMean` | mean normalised abundance |
| `log2FC` | **always** `log2(test / ref)` — verified against raw values |
| `stat` | Wald statistic (or moderated *t*) |
| `pvalue`, `padj` | raw and Benjamini–Hochberg adjusted |

Results are indexed by **reference element position**, so `res["ATAC"]` and
`res["H3K27ac"]` line up row for row and can be correlated directly.

---

## Exporting

```python
ds.export("results/")                       # tables + audit + manifest.json
ep.export_bed(ds, "ATAC", "up.bed", direction="up")
ad = ds.to_anndata("ATAC")                  # obs = samples, var = elements
mdata = ep.to_mudata(ds)                    # one modality per assay
```

`manifest.json` records the contrast, the sign convention, the audit outcome and
per-assay summaries, so the analysis can be checked later without re-running it.
