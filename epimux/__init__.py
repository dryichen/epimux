"""epimux -- bulk epigenome integration anchored on genomic intervals.

Multi-assay differential analysis (ATAC, ChIP/CUT&RUN, WGBS, Hi-C, RNA) over one
shared element reference, with the sign convention pinned and an audit layer that
catches the failure modes that silently invalidate this kind of study.

Quick start
-----------
>>> import epimux as ep
>>> ds = ep.Dataset("enhancers.bed", genome="mm10")
>>> ds.add_counts("ATAC", "atac_counts.txt")
>>> ds.add_counts("H3K27ac", "h3k_counts.txt")
>>> ds.set_design({"WT": [...], "KO": [...]})
>>> ds.differential(ref="WT", test="KO")     # log2FC = log2(KO/WT), verified
>>> ds.audit(positive_control=("LSK", "GMP"), null_group="WT")
>>> print(ds.coupling("ATAC", "H3K27ac"))
>>> ds.report("report.html")

Why the audit layer exists
--------------------------
Each check corresponds to an error that produced a wrong, confidently-reported
result in real analysis: a contrast reversed by alphabetical factor levels; a
"no change" conclusion from an underpowered signal-averaging method; a
"decoupling" that was two assays contrasted in opposite directions; a
directional bias manufactured by unequal ChIP efficiency.
"""
from .assays import (Assay, CountAssay, HiCAssay, MethylAssay, SignalAssay,
                     read_featurecounts)
from .audit import (AuditReport, AuditResult, check_direction, efficiency_balance,
                    null_contrast, positive_control, replicate_reliability)
from .core import Dataset
from .coupling import classify_elements, concordance, couple
from .diagnostics import (confounding_check, detectable_effect, outlier_replicates,
                          power_analysis, pvalue_diagnostic)
from .io import (export_bed, export_results, load_results, save_dataset,
                 to_anndata, to_mudata)
from .peaks import (consensus_peaks, jaccard, merge_intervals, peak_overlap_matrix,
                    read_narrowpeak, saf_from_intervals)
from .linking import abc_link, aggregate_to_genes, nearest_gene
from .modules import find_modules, module_enrichment, module_profile
from .annotation import (annotate_tss, classify_context, context_composition,
                         distance_enrichment, stitch_super_enhancers)
from .enrichment import (gc_matched_background, motif_enrichment,
                         overlap_enrichment, pathway_enrichment)
from .normalization import (apply_factors, assess_global_shift, median_of_ratios,
                            quantile_normalize, reference_normalize,
                            spike_in_factors, tmm)
from .stats import bh_fdr, deseq2_de, methylation_de, moderated_t_de
from .utils import Contrast, overlap, read_bed

from . import hic, plotting, tracks

__version__ = "0.1.0"

__all__ = [
    "Dataset", "Contrast",
    "CountAssay", "MethylAssay", "SignalAssay", "HiCAssay", "Assay",
    "read_featurecounts", "read_bed", "overlap",
    "deseq2_de", "moderated_t_de", "methylation_de", "bh_fdr",
    "check_direction", "positive_control", "null_contrast",
    "efficiency_balance", "replicate_reliability", "AuditReport", "AuditResult",
    "couple", "classify_elements", "concordance",
    "find_modules", "module_profile", "module_enrichment",
    "abc_link", "nearest_gene", "aggregate_to_genes",
    # normalisation (incl. spike-in)
    "median_of_ratios", "tmm", "quantile_normalize", "spike_in_factors",
    "reference_normalize", "apply_factors", "assess_global_shift",
    # annotation
    "annotate_tss", "classify_context", "context_composition",
    "stitch_super_enhancers", "distance_enrichment",
    # enrichment
    "pathway_enrichment", "motif_enrichment", "overlap_enrichment",
    "gc_matched_background",
    # diagnostics
    "outlier_replicates", "pvalue_diagnostic", "confounding_check",
    "power_analysis", "detectable_effect",
    # peaks
    "consensus_peaks", "read_narrowpeak", "merge_intervals",
    "peak_overlap_matrix", "saf_from_intervals", "jaccard",
    # io
    "to_anndata", "to_mudata", "export_bed", "export_results",
    "save_dataset", "load_results",
    "plotting", "tracks", "hic", "__version__",
]
