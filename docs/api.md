# API reference

Auto-generated from the source. Every effect size in epimux is
`log2FC = log2(test / ref)`; see [guide-audit.md](guide-audit.md) for why that is
enforced rather than assumed.

## `epimux.core` — Dataset — the container

*The Dataset: one reference element set, many assays, one sign convention.*

### `class Dataset`

Multi-assay container anchored on a shared element reference.

- **`add_assay(self, assay: 'Assay')`** — 
- **`add_counts(self, name, path=None, intervals=None, matrix=None, **kw)`** — 
- **`add_hic(self, name, files)`** — 
- **`add_methyl(self, name, pattern, **kw)`** — 
- **`add_signal(self, name, files, **kw)`** — 
- **`audit(self, positive_control: 'tuple | None' = None, null_group: 'str | None' = None, frip: 'dict | None' = None, assays: 'list | None' = None, covariates: 'pd.DataFrame | None' = None) -> '_audit.AuditReport'`** — Run the full battery of checks.
- **`classify(self, assays: 'list | None' = None, **kw) -> 'pd.DataFrame'`** — 
- **`contrast(self, ref: 'str', test: 'str') -> 'Contrast'`** — 
- **`coupling(self, x: 'str', y: 'str', **kw) -> '_coupling.CouplingResult'`** — 
- **`differential(self, ref: 'str', test: 'str', assays: 'list | None' = None, verify_direction: 'bool' = True, strict: 'bool' = True, **kw) -> 'dict'`** — Run per-assay differential analysis with one pinned direction.
- **`export(self, out_dir: 'str', **kw)`** — 
- **`link_genes(self, tss, method: 'str' = 'abc', hic: 'str | None' = None, sample: 'str | None' = None, activity: 'pd.Series | None' = None, **kw)`** — 
- **`modules(self, layers: 'dict', k: 'int' = 6, **kw) -> '_modules.ModuleResult'`** — 
- **`power(self, assay: 'str', **kw)`** — 
- **`report(self, path='epimux_report.html', **kw)`** — 
- **`set_design(self, groups: 'dict')`** — ``{"WT": [samples...], "KO": [samples...]}``
- **`significant(self, assay: 'str', fc: 'float' = 1.5, fdr: 'float' = 0.1) -> 'pd.DataFrame'`** — 
- **`summary(self) -> 'pd.DataFrame'`** — 
- **`to_anndata(self, assay: 'str')`** — 


## `epimux.assays` — Assays — one adapter per data type

*Assay adapters.*

### `class Assay`

Base class: an assay is features + a matrix + a mapping to the reference.

- **`differential(self, contrast: 'Contrast', **kw) -> 'pd.DataFrame'`** — 
- **`map_to(self, reference: 'pd.DataFrame', how: 'str' = 'best', weight: 'str | None' = None)`** — 
- **`to_elements(self, values: 'pd.Series', reference: 'pd.DataFrame', agg: 'str' = 'strongest', rank_by: 'pd.Series | None' = None) -> 'pd.Series'`** — Collapse feature-level values onto reference elements.

### `class CountAssay`

Read-count assay (ATAC, ChIP-seq, CUT&RUN, RNA). Engine: PyDESeq2.

- **`differential(self, contrast: 'Contrast', **kw) -> 'pd.DataFrame'`** — 
- **`map_to(self, reference: 'pd.DataFrame', how: 'str' = 'best', weight: 'str | None' = None)`** — 
- **`norm(self) -> 'pd.DataFrame'`** — 
- **`to_elements(self, values: 'pd.Series', reference: 'pd.DataFrame', agg: 'str' = 'strongest', rank_by: 'pd.Series | None' = None) -> 'pd.Series'`** — Collapse feature-level values onto reference elements.

### `class MethylAssay`

Bisulfite methylation. Holds methylated and total counts per element.

- **`differential(self, contrast: 'Contrast', **kw) -> 'pd.DataFrame'`** — 
- **`map_to(self, reference: 'pd.DataFrame', how: 'str' = 'best', weight: 'str | None' = None)`** — 
- **`rates(self, min_cov: 'int' = 10) -> 'pd.DataFrame'`** — 
- **`to_elements(self, values: 'pd.Series', reference: 'pd.DataFrame', agg: 'str' = 'strongest', rank_by: 'pd.Series | None' = None) -> 'pd.Series'`** — Collapse feature-level values onto reference elements.

### `class SignalAssay`

Continuous signal extracted from bigWig tracks.

- **`differential(self, contrast: 'Contrast', **kw) -> 'pd.DataFrame'`** — 
- **`map_to(self, reference: 'pd.DataFrame', how: 'str' = 'best', weight: 'str | None' = None)`** — 
- **`to_elements(self, values: 'pd.Series', reference: 'pd.DataFrame', agg: 'str' = 'strongest', rank_by: 'pd.Series | None' = None) -> 'pd.Series'`** — Collapse feature-level values onto reference elements.

### `class HiCAssay`

Hi-C contact maps: compartments, insulation, local contact support.

- **`differential(self, contrast: 'Contrast', **kw) -> 'pd.DataFrame'`** — 
- **`eigenvector(self, sample: 'str', resolution: 'int' = 160000, phasing_track: 'pd.DataFrame | None' = None) -> 'pd.DataFrame'`** — 
- **`insulation(self, sample: 'str', resolution: 'int' = 20000, window: 'int' = 100000) -> 'pd.DataFrame'`** — 
- **`local_contact(self, sample: 'str', targets: 'pd.DataFrame', resolution: 'int' = 20000, flank: 'int' = 200000) -> 'pd.Series'`** — Summed balanced contacts within +/-flank of each target midpoint.
- **`map_to(self, reference: 'pd.DataFrame', how: 'str' = 'best', weight: 'str | None' = None)`** — 
- **`to_elements(self, values: 'pd.Series', reference: 'pd.DataFrame', agg: 'str' = 'strongest', rank_by: 'pd.Series | None' = None) -> 'pd.Series'`** — Collapse feature-level values onto reference elements.

### `read_featurecounts(path, strip=('.mLb.clN.sorted.bam', '.target.markdup.sorted.bam', '.sorted.bam', '.bam', './'))`

Read a featureCounts table into (intervals, counts).


## `epimux.stats` — Differential engines

*Differential engines, one per data type.*

### `deseq2_de(counts: 'pd.DataFrame', contrast: 'Contrast', min_count: 'int' = 10, min_frac: 'float' = 0.5, covariates: 'pd.DataFrame | None' = None, n_cpus: 'int' = 4) -> 'pd.DataFrame'`

Negative-binomial DE via PyDESeq2, with an exact-sign contrast.

### `moderated_t_de(signal: 'pd.DataFrame', contrast: 'Contrast', min_var: 'float' = 1e-08) -> 'pd.DataFrame'`

Empirical-Bayes moderated t-test on continuous signal.

### `methylation_de(meth: 'pd.DataFrame', cov: 'pd.DataFrame', contrast: 'Contrast', min_cov: 'int' = 10, test: 'str' = 'beta') -> 'pd.DataFrame'`

Differential methylation on rates, coverage-weighted.

### `bh_fdr(p: 'np.ndarray') -> 'np.ndarray'`

Benjamini-Hochberg, NaN-safe.


## `epimux.audit` — Audit — is this result trustworthy?

*The audit layer -- automated checks that catch the failure modes that silently corrupt differential genomics.*

### `class AuditResult`

AuditResult(name: 'str', status: 'str', summary: 'str', detail: 'dict' = <factory>)


### `check_direction(result: 'pd.DataFrame', counts: 'pd.DataFrame', contrast: 'Contrast', kind: 'str' = 'count', top_n: 'int' = 2000, min_corr: 'float' = 0.5) -> 'AuditResult'`

Verify the sign of ``log2FC`` against raw normalized values.

### `positive_control(counts: 'pd.DataFrame', contrast: 'Contrast', de_fn, fc: 'float' = 1.5, fdr: 'float' = 0.1, min_frac: 'float' = 0.02) -> 'AuditResult'`

Run the pipeline on a comparison where a difference certainly exists.

### `null_contrast(counts: 'pd.DataFrame', samples: 'list', de_fn, fc: 'float' = 1.5, fdr: 'float' = 0.1, max_hits: 'int' = 25, seed: 'int' = 0) -> 'AuditResult'`

Split replicates *within* one biological group and run the same test.

### `efficiency_balance(frip: 'dict', contrast: 'Contrast', observed_direction: 'str | None' = None, max_ratio: 'float' = 1.2) -> 'AuditResult'`

Compare signal-to-background (e.g. FRiP) between the two groups.

### `replicate_reliability(values: 'pd.DataFrame', groups: 'dict', min_reliability: 'float' = 0.5) -> 'AuditResult'`

Spearman-Brown reliability of each group mean.

### `class AuditReport`

AuditReport(results: 'list' = <factory>)

- **`add(self, r: 'AuditResult')`** — 
- **`to_frame(self) -> 'pd.DataFrame'`** — 


## `epimux.diagnostics` — Diagnostics — what exactly is wrong?

*Diagnostics that localise *which* sample or model assumption is the problem.*

### `outlier_replicates(values: 'pd.DataFrame', groups: 'dict', max_drop: 'float' = 0.1, z_threshold: 'float' = 2.0, method: 'str' = 'spearman') -> 'AuditResult'`

Find replicates that sit far from the rest of their own group.

### `pvalue_diagnostic(result: 'pd.DataFrame', n_bins: 'int' = 20) -> 'AuditResult'`

Classify the p-value histogram: healthy, conservative, or mis-specified.

### `confounding_check(contrast, covariates: 'pd.DataFrame') -> 'AuditResult'`

Is the contrast confounded with any covariate?

### `power_analysis(result: 'pd.DataFrame', n_per_group: 'int', fc: 'float' = 1.5, fdr: 'float' = 0.1, n_range=(2, 3, 4, 5, 6, 8, 10)) -> 'pd.DataFrame'`

Replicates needed to detect the observed effect sizes.

### `detectable_effect(result: 'pd.DataFrame', n_per_group: 'int', power: 'float' = 0.8, fdr: 'float' = 0.1) -> 'float'`

Smallest fold-change detectable at the given power with the current n.


## `epimux.normalization` — Normalization

*Normalization, including the spike-in case.*

### `median_of_ratios(counts: 'pd.DataFrame', min_count: 'int' = 1) -> 'pd.Series'`

DESeq2-style size factors (geometric-mean reference).

### `tmm(counts: 'pd.DataFrame', ref_col: 'str | None' = None, log_ratio_trim: 'float' = 0.3, sum_trim: 'float' = 0.05) -> 'pd.Series'`

edgeR-style trimmed mean of M-values.

### `quantile_normalize(values: 'pd.DataFrame') -> 'pd.DataFrame'`

Force every sample to share one distribution (use with care: this is the most aggressive way to erase a genuine global shift).

### `spike_in_factors(spike_counts: 'pd.Series | dict', target_lib: 'pd.Series | dict | None' = None) -> 'pd.Series'`

Size factors from exogenous spike-in reads.

### `apply_factors(counts: 'pd.DataFrame', factors: 'pd.Series') -> 'pd.DataFrame'`



### `assess_global_shift(counts: 'pd.DataFrame', contrast, frip: 'dict | None' = None, quantiles=(0.25, 0.5, 0.75, 0.9)) -> 'dict'`

Is a genome-wide shift plausible, and would normalization hide it?

### `reference_normalize(counts: 'pd.DataFrame', reference_features) -> 'pd.Series'`

Size factors from an internal set assumed invariant (e.g. constitutive CTCF sites, housekeeping promoters).  A poor substitute for spike-ins, but far better than assuming the whole genome is invariant.


## `epimux.coupling` — Cross-layer coupling

*Cross-assay coupling, with sign conventions enforced rather than assumed.*

### `class CouplingResult`

CouplingResult(assay_x: 'str', assay_y: 'str', spearman: 'float', pearson: 'float', pvalue: 'float', n: 'int', dual_significant: 'int', same_direction: 'int', opposite_direction: 'int', strata: 'pd.DataFrame | None' = None, attenuation_corrected: 'float | None' = None, notes: 'list' = <factory>)


### `couple(res_x: 'pd.DataFrame', res_y: 'pd.DataFrame', name_x: 'str' = 'X', name_y: 'str' = 'Y', fc: 'float' = 1.5, fdr: 'float' = 0.1, reliability: 'tuple | None' = None, stratify: 'bool' = True, n_strata: 'int' = 4) -> 'CouplingResult'`

Correlate two differential results element-wise.

### `classify_elements(results: 'dict', fc: 'float' = 1.5, fdr: 'float' = 0.1) -> 'pd.DataFrame'`

Label each element by its multi-layer response.

### `concordance(cls: 'pd.DataFrame') -> 'pd.Series'`

Summary counts of the multi-layer states.


## `epimux.meta` — Cross-contrast comparison

*Comparing effects across contrasts.*

### `compare_contrasts(res_a: 'pd.DataFrame', res_b: 'pd.DataFrame', name_a: 'str' = 'A', name_b: 'str' = 'B') -> 'pd.DataFrame'`

Per-element interaction test: does the effect differ between contrasts?

### `concordance_summary(cmp: 'pd.DataFrame', name_a: 'str' = 'A', name_b: 'str' = 'B', fc: 'float' = 1.5, fdr: 'float' = 0.1) -> 'dict'`

How much of an apparent difference between contrasts is real?

### `meta_analyze(results: 'dict', method: 'str' = 'inverse_variance') -> 'pd.DataFrame'`

Pool several contrasts into one effect size per element.

### `replication_rate(discovery: 'pd.DataFrame', validation: 'pd.DataFrame', fc: 'float' = 1.5, fdr: 'float' = 0.1, loose_p: 'float' = 0.05) -> 'dict'`

Do discovery hits replicate, in direction and at a relaxed threshold?


## `epimux.modules` — Multi-omic modules

*Multi-omic module discovery over reference elements.*

### `class ModuleResult`

ModuleResult(labels: 'pd.Series', profile: 'pd.DataFrame', inertia: 'float | None' = None, method: 'str' = 'kmeans')


### `find_modules(layers: 'dict', k: 'int' = 6, method: 'str' = 'kmeans', scale: 'str' = 'zscore', seed: 'int' = 0, max_elements: 'int | None' = None) -> 'ModuleResult'`

Partition elements by their multi-omic profile.

### `module_profile(layers: 'dict', labels: 'pd.Series') -> 'pd.DataFrame'`

Mean per-assay value in each module.

### `module_enrichment(labels: 'pd.Series', selected: 'pd.Index', background: 'pd.Index | None' = None) -> 'pd.DataFrame'`

Fisher enrichment of ``selected`` elements in each module.


## `epimux.linking` — Element to gene linking

*Element -> gene linking.*

### `nearest_gene(elements: 'pd.DataFrame', tss: 'pd.DataFrame', max_distance: 'int' = 1000000) -> 'pd.DataFrame'`

Baseline nearest-TSS assignment (kept for comparison, not recommended).

### `abc_link(elements: 'pd.DataFrame', tss: 'pd.DataFrame', hic, sample: 'str', activity: 'pd.Series | None' = None, resolution: 'int' = 20000, window: 'int' = 1000000, min_score: 'float' = 0.02, exclude_promoter: 'int' = 2000) -> 'pd.DataFrame'`

Activity-by-Contact linking using a real contact map.

### `aggregate_to_genes(links: 'pd.DataFrame', element_values: 'pd.Series', weight: 'str' = 'score', how: 'str' = 'weighted_mean') -> 'pd.Series'`

Push element-level effect sizes onto genes through the link table.


## `epimux.annotation` — Annotation

*Element annotation: genomic context, TSS distance, super-enhancers.*

### `annotate_tss(elements: 'pd.DataFrame', tss: 'pd.DataFrame') -> 'pd.DataFrame'`

Signed distance to the nearest TSS (negative = upstream of the gene).

### `classify_context(elements: 'pd.DataFrame', tss: 'pd.DataFrame', promoter: 'int' = 2000, proximal: 'int' = 10000, features: 'dict | None' = None) -> 'pd.DataFrame'`

Label elements promoter / proximal / distal, plus optional feature overlap.

### `stitch_super_enhancers(peaks: 'pd.DataFrame', signal: 'pd.Series', stitch_distance: 'int' = 12500, tss: 'pd.DataFrame | None' = None, exclude_promoter: 'int' = 2000) -> 'pd.DataFrame'`

ROSE-style super-enhancer calling.

### `context_composition(ann: 'pd.DataFrame', selected: 'pd.Index | None' = None) -> 'pd.DataFrame'`

Composition of a selected set vs background, with enrichment.

### `distance_enrichment(ann: 'pd.DataFrame', selected: 'pd.Index', bins=(0, 2000, 10000, 50000, 200000, 1000000, inf)) -> 'pd.DataFrame'`

Are the selected elements distributed differently with TSS distance?


## `epimux.enrichment` — Enrichment

*Enrichment: pathways for linked genes, motifs for elements, overlap for sets.*

### `pathway_enrichment(genes, background=None, gene_sets='MSigDB_Hallmark_2020', organism: 'str' = 'mouse', top: 'int' = 15) -> 'pd.DataFrame'`

Pathway enrichment via gseapy/Enrichr.

### `motif_enrichment(counts_fg: 'pd.DataFrame', counts_bg: 'pd.DataFrame') -> 'pd.DataFrame'`

Fisher enrichment of motif hits, foreground vs an explicit background.

### `overlap_enrichment(query: 'pd.DataFrame', target: 'pd.DataFrame', universe: 'pd.DataFrame', n_shuffle: 'int' = 1000, seed: 'int' = 0) -> 'dict'`

Is ``query`` enriched for overlap with ``target``, relative to ``universe``?

### `gc_matched_background(elements: 'pd.DataFrame', selected: 'pd.Index', gc: 'pd.Series', n_per: 'int' = 5, tol: 'float' = 0.02, seed: 'int' = 0) -> 'pd.Index'`

Sample a background matched on GC content (and thus roughly on width).


## `epimux.peaks` — Peak sets

*Peak set construction.*

### `read_narrowpeak(path) -> 'pd.DataFrame'`



### `merge_intervals(df: 'pd.DataFrame', gap: 'int' = 0) -> 'pd.DataFrame'`

Merge overlapping/nearby intervals (bedtools merge equivalent).

### `consensus_peaks(peak_files: 'dict', method: 'str' = 'union', min_replicates: 'int' = 2, gap: 'int' = 0, balance: 'str' = 'raise') -> 'pd.DataFrame'`

Build a consensus peak set.

### `peak_overlap_matrix(reference: 'pd.DataFrame', peak_sets: 'dict') -> 'pd.DataFrame'`

Boolean matrix: which reference peak is supported by which sample.

### `saf_from_intervals(df: 'pd.DataFrame', path: 'str | None' = None) -> 'pd.DataFrame'`

SAF table for featureCounts (1-based, inclusive).

### `jaccard(a: 'pd.DataFrame', b: 'pd.DataFrame') -> 'float'`

Jaccard index of two interval sets (by covered bases).


## `epimux.hic` — Hi-C analytics

*Hi-C analytics beyond compartments and insulation.*

### `contact_decay(source, resolution: 'int' = 20000) -> 'pd.DataFrame'`

P(s): contact probability vs genomic separation.

### `extrusion_shoulder(decay: 'pd.DataFrame', s_min: 'int' = 50000, s_max: 'int' = 5000000) -> 'pd.DataFrame'`

Log-derivative of P(s).

### `saddle(source, eigenvector: 'pd.DataFrame', sample_col: 'str', resolution: 'int' = 160000, n_bins: 'int' = 50, q_lo: 'float' = 0.025, q_hi: 'float' = 0.975)`

Compartment saddle: observed/expected binned by eigenvector rank.

### `compartment_strength(S: 'np.ndarray', frac: 'float' = 0.2) -> 'float'`

(AA + BB) / 2AB from a saddle matrix -- higher = more segregated.

### `pileup(source, features: 'pd.DataFrame', resolution: 'int' = 10000, flank: 'int' = 200000, expected: 'bool' = True) -> 'np.ndarray'`

On-diagonal pileup (average map centered on each feature).

### `apa(source, loops: 'pd.DataFrame', resolution: 'int' = 10000, flank: 'int' = 100000) -> 'tuple'`

Aggregate peak analysis over off-diagonal loop anchors.

### `boundary_strength(insulation: 'pd.DataFrame', sample_col: 'str', quantile: 'float' = 0.05) -> 'pd.DataFrame'`

Flag the strongest boundaries (lowest insulation score).

### `differential_insulation(ins_ref: 'pd.DataFrame', ins_test: 'pd.DataFrame', ref_col: 'str', test_col: 'str', anchors: 'pd.DataFrame | None' = None, resolution: 'int' = 20000) -> 'pd.DataFrame'`

Insulation change per bin, optionally annotated by anchor overlap.


## `epimux.plotting` — Plotting

*Publication-grade plotting.*

### `set_style(font: 'str | None' = None)`

Apply the house style; returns the FontProperties actually used.

### `ma_plot(res: 'pd.DataFrame', fc=1.5, fdr=0.1, title=None, ax=None, up_label=None, down_label=None, ylim=(-4, 4))`

Effect size vs abundance, significance colored by direction.

### `volcano(res: 'pd.DataFrame', fc=1.5, fdr=0.1, title=None, ax=None, label_top=0)`



### `coupling_plot(res_x, res_y, name_x='X', name_y='Y', fc=1.5, fdr=0.1, lim=3, ax=None, title=None)`

Cross-assay effect-size scatter with the dual-significant set highlighted.

### `signal_heatmap(matrix: 'pd.DataFrame', groups: 'dict | None' = None, sort_by: 'pd.Series | None' = None, vmax=1.5, title=None, ylabel=None, ax=None)`

Row z-scored replicate heatmap of selected elements.

### `audit_plot(report, ax=None)`

Traffic-light summary of the audit battery.

### `state_barplot(states: 'pd.Series', ax=None, title=None)`

Counts of multi-layer states (coordinated / discordant / ...).

### `pca_plot(matrix: 'pd.DataFrame', groups: 'dict', n_top=2000, ax=None, title=None, markers=None, colors=None)`



### `save(fig, path, formats=('png', 'pdf'))`




## `epimux.tracks` — Browser tracks

*Genome-browser-style figure panels.*

### `extract_profile(bigwig, chrom: 'str', start: 'int', end: 'int', nbins: 'int' = 500)`

Binned mean signal over a locus.

### `locus_plot(tracks: 'dict', chrom: 'str', start: 'int', end: 'int', groups: 'dict | None' = None, features: 'dict | None' = None, nbins: 'int' = 600, colors: 'dict | None' = None, title: 'str | None' = None, height_per_track: 'float' = 0.9, show_spread: 'bool' = True)`

Stacked browser view.

### `metaprofile(bigwigs: 'dict', regions: 'pd.DataFrame', flank: 'int' = 2000, nbins: 'int' = 100, groups: 'dict | None' = None, colors: 'dict | None' = None, ax=None, title: 'str | None' = None, center: 'str' = 'midpoint')`

Average signal centered on a set of regions.

### `heatmap_profile(bigwig, regions: 'pd.DataFrame', flank: 'int' = 2000, nbins: 'int' = 100, sort_by: 'pd.Series | None' = None, vmax: 'float | None' = None, ax=None, title: 'str | None' = None)`

Per-region signal heatmap (deepTools-style), rows optionally sorted.


## `epimux.io` — Import / export

*Import/export and interoperability.*

### `to_anndata(dataset, assay: 'str')`

One assay as AnnData (obs = samples, var = reference elements).

### `to_mudata(dataset)`

All assays as MuData, one modality per assay.

### `export_bed(dataset, assay: 'str', path: 'str', fc: 'float' = 1.5, fdr: 'float' = 0.1, direction: 'str | None' = None, name_prefix: 'str' = '') -> 'str'`

Significant elements as BED, score = 1000*min(1,|log2FC|/3), strand-free.

### `export_results(dataset, out_dir: 'str', fc: 'float' = 1.5, fdr: 'float' = 0.1) -> 'dict'`

Write every result table, the audit record and a machine-readable manifest.

### `save_dataset(dataset, path: 'str') -> 'str'`

Serialise results + audit + design to a single JSON (data stays on disk).

### `load_results(path: 'str') -> 'dict'`

Read back the result tables written by :func:`save_dataset`.


## `epimux.report` — HTML report

*Self-contained HTML report: audit traffic-lights, per-assay results, coupling.*

### `html_report(ds, path='epimux_report.html', fc=1.5, fdr=0.1, figures=True)`

Render a Dataset into a single self-contained HTML file.


## `epimux.utils` — Interval utilities

*Interval primitives, logging and small helpers.*

### `get_logger(name: 'str' = 'epimux') -> 'logging.Logger'`



### `read_bed(path, name_col: 'bool' = True) -> 'pd.DataFrame'`

Read a BED-like file into chrom/start/end[/name].

### `as_intervals(df: 'pd.DataFrame') -> 'pd.DataFrame'`

Coerce a frame to canonical interval dtypes.

### `sort_intervals(df: 'pd.DataFrame') -> 'pd.DataFrame'`



### `overlap(a: 'pd.DataFrame', b: 'pd.DataFrame', min_overlap: 'int' = 1) -> 'pd.DataFrame'`

All overlapping pairs between two interval sets.

### `map_to_reference(features: 'pd.DataFrame', reference: 'pd.DataFrame', how: 'str' = 'best', weight: 'str | None' = None) -> 'pd.Series'`

Map assay features onto reference elements.

### `midpoint_bin(df: 'pd.DataFrame', resolution: 'int') -> 'pd.Series'`

chrom_binstart key at a fixed resolution (for Hi-C style binning).

### `cpm(counts: 'np.ndarray') -> 'np.ndarray'`



### `log2cpm(counts: 'np.ndarray', prior: 'float' = 1.0) -> 'np.ndarray'`




## `epimux.cli` — Command line

*Command-line interface.*

### `build_dataset(cfg)`



### `cmd_audit(args)`



### `cmd_run(args)`



### `get_logger(name: 'str' = 'epimux') -> 'logging.Logger'`



### `main(argv=None)`



