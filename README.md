# epimux

[![CI](https://github.com/dryichen/epimux/actions/workflows/ci.yml/badge.svg)](https://github.com/dryichen/epimux/actions/workflows/ci.yml)
[![Python](https://img.shields.io/badge/python-3.10%20%7C%203.11%20%7C%203.12-blue)](https://www.python.org)
[![License: MIT](https://img.shields.io/badge/license-MIT-green)](LICENSE)

**Bulk epigenome integration anchored on genomic intervals — with an audit layer that catches the errors that silently invalidate multi-omic studies.**

```python
import epimux as ep

ds = ep.Dataset("enhancers.bed", genome="mm10")
ds.add_counts("ATAC",    "atac_counts.txt")          # featureCounts
ds.add_counts("H3K27ac", "h3k_counts.txt")
ds.add_methyl("WGBS",    "wgbs/*.cov")               # Bismark
ds.add_hic("HiC", {"WT": "wt.mcool", "KO": "ko.mcool"})

ds.set_design({"WT": [...], "KO": [...], "GMP": [...]})
ds.differential(ref="WT", test="KO")                 # log2FC = log2(KO/WT), verified
ds.audit(positive_control=("WT", "GMP"), null_group="WT")
print(ds.coupling("ATAC", "H3K27ac"))
ds.report("report.html")
```

---

## Documentation

| | |
|---|---|
| [**Tutorial**](docs/tutorial.ipynb) | full walkthrough on data with a known ground truth |
| [**Input formats**](docs/guide-inputs.md) | what each assay expects, design, covariates, exports |
| [**Reading the audit**](docs/guide-audit.md) | what each check tests and how to act on it |
| [**Choosing a normalisation**](docs/guide-normalisation.md) | spike-in, internal reference, and when a magnitude is unrecoverable |
| [**Hi-C analyses**](docs/guide-hic.md) | compartments, P(s), insulation, APA, ABC linking |
| [**Command line**](docs/guide-cli.md) | `epimux run/audit --config study.yaml` |
| [**FAQ**](docs/faq.md) | troubleshooting |
| [**API reference**](docs/api.md) | every public function |

## Install

```bash
pip install -e ".[all]"      # core + hic + bigwig + enrichment
pip install -e ".[dev]" && pytest
```

Core needs only numpy / pandas / scipy / scikit-learn / matplotlib / pydeseq2.
Hi-C needs `cooler`, `cooltools`, `bioframe`; bigWig needs `pyBigWig`; pathway
enrichment needs `gseapy`.

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
| `check_direction` | A contrast built from an R factor computed `log2(WT/KO)` because factor levels sort **alphabetically**. Every direction in the analysis was reversed; nothing downstream complained. This check recomputes the effect from raw normalised values and **fails loudly** on a sign flip. |
| `positive_control` | A bigWig/limma pipeline reported "no change" for a histone mark. It was a **false negative of a weak method** — the same pipeline could not detect a cell-type difference either. A pipeline that cannot find a difference you *know* exists cannot support a null result. |
| `null_contrast` | Splitting replicates *within* one group must yield ~no hits. If it does not, the "significant" features in the real contrast are noise. |
| `efficiency_balance` | Unequal ChIP efficiency (FRiP) between genotypes can manufacture a directional bias that survives depth normalisation. Reports whether the observed effect runs **with** the technical bias (dangerous) or **against** it (conservative). |
| `replicate_reliability` | Low reliability attenuates every cross-assay correlation by `sqrt(r_x · r_y)` — the reason single-replicate tracks yield correlations with the wrong magnitude and sometimes the wrong sign. |
| `outlier_replicates` | Which single library is driving the result. Uses the drop relative to group-mates, because the most extreme z-score with `n` samples is bounded by `(n-1)/sqrt(n)` — **1.15 at n=3** — so a z-threshold can never fire in a three-replicate design. |
| `pvalue_diagnostic` | A mis-specified variance model, read off the p-value histogram, before the FDR is believed. |
| `confounding_check` | A design where a covariate cannot be separated from the contrast — no adjustment can fix it, and it should be reported as confounded. |

`differential()` runs `check_direction` **by default** and raises rather than
returning reversed results. See [docs/guide-audit.md](docs/guide-audit.md).

## Design principles

**The contrast is an object, not a factor.** `Contrast(ref="WT", test="KO")`
makes `log2FC = log2(test/ref)` structurally impossible to get backwards, and it
is re-verified against the raw matrix at runtime.

**One reference, many assays.** Every assay is projected onto a shared element
set, so ATAC peaks, ChIP peaks, CpGs and Hi-C bins become directly comparable
without repeated interval arithmetic.

**Assay-appropriate statistics.** Counts → PyDESeq2 (negative binomial);
methylation → coverage-weighted, variance-stabilised rates; averaged signal →
moderated *t* (and the class warns you it is less sensitive than counts).

**Coupling refuses to guess.** `couple()` raises if the two results were built
from opposite orientations, reports the attenuation-corrected estimate, and
stratifies by abundance to expose shrinkage-induced correlation.

**States from significance, not signs.** `classify_elements()` labels elements by
significance in each layer. Sign-only state calls were the original source of
irreproducible "discordant element" lists.

**Genes by contact, not proximity.** `link_genes(method="abc")` uses the contact
map (Fulco-style Activity-by-Contact). Nearest-gene is available for comparison
and clearly labelled as not recommended.

**Defaults refuse the biased option.** `consensus_peaks(method="replicated")`
raises on unequal replicate numbers rather than quietly producing a peak set
that favours the larger group.

---

## What you can do with it

### Differential analysis with a pinned direction

```python
ds.differential(ref="WT", test="KO")                     # verified log2(KO/WT)
ds.differential(ref="WT", test="KO", covariates=batch)   # ~batch + condition
ds.significant("ATAC", fc=1.5, fdr=0.1)
```

### Audit and diagnose

```python
ds.audit(positive_control=("WT", "GMP"), null_group="WT", covariates=batch)
ep.outlier_replicates(log_cpm, groups)      # which library is deviant?
ep.pvalue_diagnostic(result)                # does the model fit?
ds.power("ATAC")                            # what effect size can this n detect?
```

### Normalisation when a global shift is possible

```python
ep.assess_global_shift(counts, contrast, frip=frip)   # is the assumption at risk?
sf = ep.spike_in_factors(spike_counts, target_lib)    # survives a global shift
counts_norm = ep.apply_factors(counts, sf)
```

Without spike-ins, `reference_normalize` against a set you believe is invariant
is a weaker but honest fallback — and `assess_global_shift` at least tells you
when to stop claiming a magnitude.
See [docs/guide-normalisation.md](docs/guide-normalisation.md).

### Cross-layer and cross-contrast

```python
ds.coupling("ATAC", "H3K27ac")               # coordinated or discordant?
states = ds.classify(); ep.concordance(states)

cmp = ep.compare_contrasts(res_lsk, res_gmp, "LSK", "GMP")   # interaction test
ep.concordance_summary(cmp, "LSK", "GMP")    # power difference or real difference?
ep.meta_analyse({"cohort1": r1, "cohort2": r2})
```

`compare_contrasts` exists because comparing two significance lists does **not**
answer "is the effect the same in A and B?" — lists differ because of replicate
count, depth and dispersion.

### Hi-C without spike-in

When a ChIP is too shallow to quantify a binding change, ask about its
*functional consequence* instead. Hi-C needs no spike-in and is usually deep:

```python
wt, ko = hic.insulation("WT"), hic.insulation("KO")
d = ep.hic.differential_insulation(wt, ko, "WT", "KO", anchors=ctcf_peaks)

decay = ep.hic.contact_decay("wt.mcool", resolution=20_000)
sh    = ep.hic.extrusion_shoulder(decay)          # log-derivative of P(s)
S, _  = ep.hic.saddle("wt.mcool", eigenvector, "WT")
print(ep.hic.compartment_strength(S))             # (AA+BB)/2AB
mat, score = ep.hic.apa("wt.mcool", loops)
```

See [docs/guide-hic.md](docs/guide-hic.md).

### Peak sets, annotation, enrichment

```python
peaks = ep.consensus_peaks(files, method="replicated", balance="subsample")
ep.saf_from_intervals(peaks, "consensus.saf")        # then featureCounts

ann = ep.classify_context(elements, tss)             # promoter / proximal / distal
se  = ep.stitch_super_enhancers(peaks, signal, tss=tss)

ep.pathway_enrichment(genes, background=tested_genes)   # background-aware
ep.motif_enrichment(fg, ep.gc_matched_background(...))  # GC-matched
ep.overlap_enrichment(query, target, universe)          # label-shuffling null
```

### Figures and export

```python
from epimux import plotting as pl
pl.ma_plot(res); pl.volcano(res); pl.coupling_plot(a, b)
pl.pca_plot(counts, groups); pl.audit_plot(ds.audit_report)

ep.tracks.locus_plot(bigwigs, "chr1", 3_000_000, 3_200_000, groups=groups)
ep.tracks.metaprofile(bigwigs, regions, flank=2000, groups=groups)

ds.report("report.html")          # self-contained, theme-aware
ds.export("results/")             # tables + audit + manifest.json
ep.export_bed(ds, "ATAC", "up.bed", direction="up")
ad = ds.to_anndata("ATAC")
```

Plots use a journal-style palette, embed Type 42 fonts so PDFs stay editable,
and rasterise scatter layers so vector files stay small.

---

## Module map

| Module | Contents |
|---|---|
| `core` | `Dataset` — reference, assays, design, contrast, orchestration |
| `assays` | `CountAssay`, `MethylAssay`, `SignalAssay`, `HiCAssay` |
| `stats` | `deseq2_de`, `methylation_de`, `moderated_t_de`, `bh_fdr` |
| `audit` | the five core checks + `AuditReport` |
| `diagnostics` | outlier replicates, p-value shape, confounding, power analysis |
| `normalization` | median-of-ratios, TMM, quantile, **spike-in**, internal reference, global-shift assessment |
| `coupling` | `couple`, `classify_elements`, `concordance` |
| `meta` | `compare_contrasts`, `concordance_summary`, `meta_analyse`, `replication_rate` |
| `modules` | `find_modules`, `module_profile`, `module_enrichment` |
| `linking` | `abc_link`, `nearest_gene`, `aggregate_to_genes` |
| `annotation` | TSS distance, genomic context, ROSE-style super-enhancers |
| `enrichment` | pathways (background-aware), motifs (GC-matched), interval overlap |
| `peaks` | consensus construction, SAF, Jaccard, overlap matrix |
| `hic` | P(s) + log-derivative, saddle, compartment strength, pileup, APA, insulation |
| `plotting` | MA, volcano, coupling, heatmaps, PCA, state bars, audit traffic-lights |
| `tracks` | browser-style locus panels, metaprofiles, per-region heatmaps |
| `io` | AnnData/MuData export, BED export, result manifest |
| `report` | one self-contained, theme-aware HTML file |
| `cli` | `epimux run/audit --config study.yaml` |

## Command line

```bash
epimux audit --config study.yaml          # checks only; non-zero exit on failure
epimux run   --config study.yaml --out results/
```

Both return exit code 1 when an audit check fails, so they can gate a pipeline.
Config reference in [docs/guide-cli.md](docs/guide-cli.md).

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md). The short version: new differential
engines must take a `Contrast` and return `log2(test/ref)`, every audit check
needs a test that makes it fail, enrichment functions take an explicit
background, and weak methods must warn at runtime.

## Citation

If `epimux` is useful, please cite the accompanying study (in preparation) and
this repository.

## License

MIT
