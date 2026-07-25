# Changelog

All notable changes to epimux are documented here.
Format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [0.1.0] — 2026-07-24

First release.

### Core
- `Dataset`: one genomic-interval reference, many assays, one pinned contrast.
- `Contrast` as an object rather than a factor — `log2FC = log2(test/ref)` is
  structurally unable to invert, and is re-verified against raw values at runtime.
- Assay adapters: `CountAssay` (featureCounts), `MethylAssay` (Bismark `.cov`),
  `SignalAssay` (bigWig, with a sensitivity warning), `HiCAssay` (mcool).
- Differential engines: PyDESeq2 negative binomial, coverage-weighted
  variance-stabilised methylation, empirical-Bayes moderated *t*.
- Covariate-adjusted designs (`~batch + condition`).

### Audit layer
- `check_direction`, `positive_control`, `null_contrast`, `efficiency_balance`,
  `replicate_reliability`, and an `AuditReport` with traffic-light rendering.
- `differential()` verifies direction by default and raises on a sign flip.

### Normalisation
- `median_of_ratios`, `tmm`, `quantile_normalize`.
- `spike_in_factors` and `reference_normalize` for the case where a genuine
  global shift would otherwise be absorbed by size-factor normalisation.
- `assess_global_shift` warns when that assumption is at risk.

### Analysis
- Cross-layer `couple()` (refuses mismatched contrast orientations, reports
  attenuation-corrected estimates, stratifies by abundance).
- `classify_elements` — multi-layer states from significance, not signs.
- ABC gene linking through a contact map; nearest-gene retained for comparison.
- Multi-omic modules (k-means / NMF) with Fisher enrichment.
- Annotation: TSS distance, genomic context, ROSE-style super-enhancers.
- Enrichment: pathways (background-aware), motifs (GC-matched background),
  interval overlap with a label-shuffling null.
- Hi-C: P(s) and its log-derivative, saddle plots and compartment strength,
  on-diagonal pileup, APA, boundary strength, differential insulation.

### Output
- Publication-grade plotting (MA, volcano, coupling, heatmaps, PCA, audit).
- Browser-style locus tracks, metaprofiles and per-region heatmaps.
- Self-contained, theme-aware HTML report.
- CLI: `epimux run|audit --config study.yaml`.

## [Unreleased]

### Added
- `diagnostics`: `outlier_replicates` (drop-based, because a z-score is bounded
  by `(n-1)/sqrt(n)` and can never fire at n=3), `pvalue_diagnostic` for model
  mis-specification, `confounding_check`, `power_analysis`, `detectable_effect`.
  Wired into `Dataset.audit()`.
- `peaks`: consensus peak construction that **refuses** the replicate-imbalanced
  "present in >= k" default, plus `merge_intervals`, `peak_overlap_matrix`,
  `saf_from_intervals`, `jaccard`.
- `io`: AnnData/MuData export, BED export by direction, and `export_results`,
  which writes a manifest recording the contrast, the sign convention and the
  audit outcome alongside the tables.
- `Dataset.export()`, `Dataset.to_anndata()`, `Dataset.power()`.
- `docs/tutorial.ipynb` — full walkthrough on data with a known ground truth.
- CI badge; 16 further tests (45 total).

### Fixed
- `outlier_replicates` could never flag a deviant replicate in a three-replicate
  design: the z-score criterion is mathematically bounded below the threshold at
  small n. Now uses the drop relative to group-mates, with the z-score applied
  only when a group has >= 5 replicates.
- `meta`: `compare_contrasts` (per-element interaction test — the question people
  actually mean by "is the effect the same in A and B?"), `concordance_summary`
  (how much of an apparent difference is power rather than biology),
  `meta_analyse` with Cochran's Q, `replication_rate`.
- Documentation: `docs/api.md` (auto-generated, 102 entries), plus guides for
  input formats, reading the audit, choosing a normalisation, Hi-C analyses and
  the CLI, and a troubleshooting FAQ. README expanded with worked examples for
  every module. All documented APIs are verified to exist.
