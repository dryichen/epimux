# epimux

**Bulk epigenome integration anchored on genomic intervals — with an audit layer that catches the errors that silently invalidate multi-omic studies.**

```python
import epimux as ep

ds = ep.Dataset("enhancers.bed", genome="mm10")
ds.add_counts("ATAC",    "atac_counts.txt")      # featureCounts
ds.add_counts("H3K27ac", "h3k_counts.txt")
ds.add_methyl("WGBS",    "wgbs/*.cov")           # Bismark
ds.add_hic("HiC", {"WT": "wt.mcool", "KO": "ko.mcool"})

ds.set_design({"WT": [...], "KO": [...], "GMP": [...]})
ds.differential(ref="WT", test="KO")             # log2FC = log2(KO/WT), verified
ds.audit(positive_control=("WT", "GMP"), null_group="WT")
print(ds.coupling("ATAC", "H3K27ac"))
ds.report("report.html")
```

---

## Why this exists

The multi-omics tooling landscape is crowded — but it is crowded in *single-cell*
space (`muon`, `Seurat WNN`, `ArchR`, `Signac`) or in *feature-agnostic* space
(`MOFA+`, `mixOmics`). For **bulk epigenomics**, where the integration key is a
genomic coordinate and each assay needs its own statistics, Python has no
equivalent of Bioconductor's `MultiAssayExperiment` — and nothing that closes the
loop from intervals → per-assay statistics → cross-layer inference → 3D-aware
gene links.

`epimux` fills that gap, and adds something none of them have: **an audit layer.**

## The audit layer

Every check exists because it caught a real error that had already been written
into a figure and a slide deck:

| Check | The error it catches |
|---|---|
| `check_direction` | A contrast built from an R factor computed `log2(WT/KO)` because factor levels sort **alphabetically**. Every direction in the paper was reversed; nothing downstream complained. This check recomputes the effect from raw normalised values and **fails loudly** on a sign flip. |
| `positive_control` | A bigWig/limma pipeline reported "no change" for a histone mark. It was a **false negative of a weak method** — the same pipeline could not detect a cell-type difference either. A pipeline that cannot find a difference you *know* exists cannot support a null result. |
| `null_contrast` | Splitting replicates *within* one group must yield ~no hits. If it does not, the "significant" features in the real contrast are noise. |
| `efficiency_balance` | Unequal ChIP efficiency (FRiP) between genotypes can manufacture a directional bias that survives depth normalisation. Reports whether the observed effect runs **with** the technical bias (dangerous) or **against** it (conservative). |
| `replicate_reliability` | Low reliability attenuates every cross-assay correlation by `sqrt(r_x · r_y)` — the reason single-replicate tracks yield correlations with the wrong magnitude and sometimes the wrong sign. |

`differential()` runs `check_direction` **by default** and raises rather than
returning reversed results.

## Design principles

**The contrast is an object, not a factor.** `Contrast(ref="WT", test="KO")`
makes `log2FC = log2(test/ref)` structurally impossible to get backwards, and it
is re-verified against the raw matrix at runtime.

**One reference, many assays.** Every assay is projected onto a shared element
set, so ATAC peaks, ChIP peaks, CpGs and Hi-C bins become directly comparable
without repeated interval arithmetic.

**Assay-appropriate statistics.** Counts → PyDESeq2 (negative binomial);
methylation → coverage-weighted, variance-stabilised rates; averaged signal →
moderated *t* (and the class warns you that it is less sensitive than counts).

**Coupling refuses to guess.** `couple()` raises if the two results were built
from opposite orientations, reports the attenuation-corrected estimate, and
stratifies by abundance to expose shrinkage-induced correlation.

**States from significance, not signs.** `classify_elements()` labels elements by
significance in each layer. Sign-only state calls were the original source of
irreproducible "decoupled element" lists.

**Genes by contact, not proximity.** `link_genes(method="abc")` uses the contact
map (Fulco-style Activity-by-Contact). Nearest-gene is available for comparison
and clearly labelled as not recommended.

## Install

```bash
pip install -e ".[all]"      # core + hic + bigwig + enrichment
pip install -e ".[dev]" && pytest
```

Core needs only numpy/pandas/scipy/scikit-learn/matplotlib/pydeseq2.
Hi-C features need `cooler`, `cooltools`, `bioframe`; bigWig needs `pyBigWig`.

## Module map

| Module | Contents |
|---|---|
| `core` | `Dataset` — reference, assays, design, contrast, orchestration |
| `assays` | `CountAssay`, `MethylAssay`, `SignalAssay`, `HiCAssay` |
| `stats` | `deseq2_de`, `methylation_de`, `moderated_t_de`, `bh_fdr` |
| `audit` | the five checks + `AuditReport` |
| `coupling` | `couple`, `classify_elements`, `concordance` |
| `linking` | `abc_link`, `nearest_gene`, `aggregate_to_genes` |
| `modules` | `find_modules`, `module_profile`, `module_enrichment` |
| `plotting` | MA, volcano, coupling, heatmap, PCA, state bars, audit traffic-lights |
| `report` | one self-contained, theme-aware HTML file |

## Hi-C without spike-in

When a ChIP is too shallow to quantify a binding change, ask about its
*functional consequence* instead. `HiCAssay.insulation()` and `.eigenvector()`
need no spike-in and are usually deeply sequenced:

```python
hic = ds.assays["HiC"]
wt, ko = hic.insulation("WT"), hic.insulation("KO")
```

## Citation

If `epimux` is useful, please cite the accompanying study (in preparation) and
this repository.

## License

MIT
