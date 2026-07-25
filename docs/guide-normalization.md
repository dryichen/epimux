# Choosing a normalization

Every size-factor method — median-of-ratios, TMM, quantile — rests on one
assumption:

> **most features do not change.**

When that holds they are excellent. When it breaks — a genome-wide gain or loss
of a mark — they absorb the very effect you are measuring and hand back a
confident null. This is the single most consequential methodological choice in a
ChIP-style differential analysis, and it is usually made by default.

---

## Decide in this order

**1. Could the change be global?**

```python
ep.assess_global_shift(counts, ds.contrast("WT", "KO"), frip=frip)
```

It compares depth-normalized signal quantiles between groups. A consistent
offset across *all* quantiles is what a true global change looks like — and is
exactly what size-factor normalization removes.

```python
{'q0.25': 0.31, 'q0.5': 0.29, 'q0.75': 0.28, 'q0.9': 0.27,
 'consistent_shift': True,
 'interpretation': 'global shift plausible -- spike-ins required to quantify it;
                    report relative redistribution only'}
```

**2. Do you have spike-ins?**

```python
sf = ep.spike_in_factors(spike_counts, target_lib)
counts_norm = ep.apply_factors(counts, sf)
```

Spike-in factors are derived from exogenous chromatin, so a genome-wide change
in the target does not move them. Passing `target_lib` additionally corrects for
sequencing depth, so the factors reflect biological occupancy rather than how
deeply each library was run.

**3. No spike-ins — what can still be said?**

Two honest options.

*Internal reference.* Normalize on a set you have reason to believe is
invariant (constitutive CTCF sites, housekeeping promoters):

```python
sf = ep.reference_normalize(counts, invariant_features)
```

Weaker than a spike-in — the assumption has just moved from "most of the genome"
to "this set" — but it is explicit and testable rather than hidden.

*Report the relative question.* Standard normalization still answers
**"which sites changed relative to the average site?"** perfectly well. That is
a real, publishable question. What you cannot claim without spike-ins is the
**global magnitude**. Say which one you are answering.

---

## When the ChIP is simply too shallow

Normalization cannot rescue depth. If a positive control fails and library sizes
are far below the other assays, the binding question is unanswerable with that
data at any normalization.

Ask about the **functional consequence** instead, using an assay that needs no
spike-in and is deeply sequenced:

```python
wt, ko = hic.insulation("WT"), hic.insulation("KO")
d = ep.hic.differential_insulation(wt, ko, "WT", "KO", anchors=ctcf_peaks)
```

If binding at a class of sites changed, boundary insulation there should change.
This converts an unanswerable question into a measurable one.

---

## Reference

| function | use |
|---|---|
| `median_of_ratios` | default; DESeq2-style geometric reference |
| `tmm` | edgeR-style; more robust when composition differs |
| `quantile_normalize` | forces identical distributions — the most aggressive way to erase a real global shift; use deliberately |
| `spike_in_factors` | exogenous reference; survives a global change |
| `reference_normalize` | internal invariant set; explicit fallback |
| `assess_global_shift` | is the assumption at risk? |
| `apply_factors` | divide counts by size factors |

`deseq2_de` applies median-of-ratios internally. To use spike-in factors,
normalize first and pass the adjusted matrix, or supply the factors to the
engine directly.
