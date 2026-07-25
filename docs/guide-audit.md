# Reading the audit

`ds.audit()` returns a traffic-light report. This guide explains what each check
actually tests, what a failure means, and what to do about it — because a check
you cannot act on is decoration.

```python
ds.differential(ref="WT", test="KO")
ds.audit(positive_control=("WT", "GMP"), null_group="WT")
print(ds.audit_report)
ds.audit_report.to_frame()
```

---

## `check_direction` — is the sign right?

**Tests:** recomputes a naive `log2(mean(test) / mean(ref))` from the raw matrix
and correlates it with the reported `log2FC`.

**FAIL means:** the contrast was built backwards. Every direction in your
analysis is inverted.

This runs automatically inside `differential()` and **raises** rather than
returning reversed results (`strict=False` to override).

**Why it exists.** In R, `factor(c("WT","KO"))` sorts its levels alphabetically
to `("KO","WT")`. A contrast written as "WT vs KO" then computes `log2(WT/KO)` —
the reverse of what the author intended. Nothing downstream complains: the
p-values are identical, the volcano plot looks normal, the peak counts are
right. Only the direction is wrong, and direction is the biology. In epimux the
contrast is an object (`Contrast(ref=..., test=...)`), so this class of bug
cannot be expressed — and it is still verified against the data at runtime.

**If it fails:** swap `ref` and `test`. Then check every figure and sentence you
have already written.

---

## `positive_control` — does the pipeline have power?

**Tests:** the same pipeline on a comparison where a difference certainly
exists — two cell types, a treated vs untreated pair, anything you already know
differs.

```python
ds.audit(positive_control=("WT", "GMP"))     # cell type, not genotype
```

**FAIL means:** the pipeline detects almost nothing even where it must. **No
null result from it is interpretable.**

**Why it exists.** A histone-mark analysis once reported "no change" between
genotypes. The number was real; the conclusion was not. The same pipeline could
not detect the difference between two cell types either — it was a false
negative of an underpowered method (averaged bigWig signal instead of read
counts), dressed up as a biological finding.

**Rule of thumb:** if the positive control finds 30–40% of features different
and your genotype contrast finds 1–2%, that contrast is a *specific* effect. If
the positive control finds ~0, stop and fix the data or the method.

**If it fails:** use read counts rather than averaged signal; check depth; check
that the samples are what you think they are (PCA).

---

## `null_contrast` — what is the false-positive rate?

**Tests:** splits replicates *within* one biological group and runs the same
test. Every hit is a false positive by construction.

**FAIL means:** the "significant" features in your real contrast are not
distinguishable from noise.

**Interpreting the number.** Zero to a handful is healthy. A hundred hits
against a real contrast of a thousand is a *warning*, not a disqualification —
it tells you the assay has more within-group variability than you assumed and
that modest effects deserve care. Compare across assays: if ATAC gives 0 and the
same design in H3K27ac gives 100, the difference is assay noise, not biology.

**If it fires:** run `outlier_replicates` next — a single deviant sample is the
usual cause, and it is invisible in a group-mean comparison.

---

## `efficiency_balance` — could a technical bias explain the direction?

**Tests:** compares signal-to-background (FRiP) between the two groups, and —
this is the useful part — reports whether the technical bias runs **with** or
**against** your observed effect.

**FAIL means:** the bias would produce the direction you are reporting. The
result may be technical.

**PASS with an imbalance means something stronger than "no problem":** the bias
runs opposite to your effect, so the true effect is *larger* than measured. This
is a sentence worth putting in a paper:

> KO libraries have lower ChIP efficiency (FRiP 33.7% vs 42.3%), which biases
> toward apparent loss; the observed gain therefore runs against the technical
> bias and the estimate is conservative.

**If it fails:** re-run on an efficiency-matched subset, drop the worst library,
or use spike-in normalisation (`ep.spike_in_factors`).

---

## `replicate_reliability` — are group means stable?

**Tests:** Spearman–Brown reliability of each group mean.

**Why it matters:** low reliability attenuates every cross-assay correlation by
`sqrt(r_x · r_y)`. This is the reason single-replicate tracks produce
correlations with the wrong *magnitude*, and — when two layers are contrasted
inconsistently — sometimes the wrong *sign*. `couple()` reports an
attenuation-corrected estimate when reliability is available.

---

## `outlier_replicates` — which sample is the problem?

**Tests:** how far each sample sits from its own group-mates, measured as the
**drop** relative to the median of the others.

**Note on method:** a z-score alone is useless here. With `n` samples the most
extreme z (population SD) is bounded by `(n-1)/sqrt(n)` — **1.15 at n = 3**. A
threshold of 2 can never fire in a typical three-replicate design, so a deviant
replicate would always be missed. epimux uses the drop criterion below n = 5 and
adds the z-score only above it.

**If it warns:** re-run the contrast without the flagged sample and compare. If
the result holds, say so; if it collapses, the finding was one library.

---

## `pvalue_diagnostic` — does the model fit?

**Tests:** the shape of the p-value histogram.

| shape | meaning |
|---|---|
| flat with a spike near 0 | healthy — signal over a proper null |
| flat, no spike | no signal; the null is well behaved |
| **pile-up near 1** | variance over-estimated — **FDR not interpretable** |
| **hump in the middle** | mis-specification or unmodelled batch structure |

**If it fails:** add the missing covariate (`covariates=`), check for batch
structure in PCA, or reconsider the test. Do not report FDR from a
mis-specified model — the number is not what it claims to be.

---

## `confounding_check` — can the design answer the question?

**Tests:** whether a covariate is separable from the contrast.

**FAIL means:** e.g. every KO was processed in batch 2. **No statistical
adjustment can separate them.** The honest response is to report the design as
confounded, not to "adjust for batch" and present the result as clean.

---

## Acting on the report

```python
rep = ds.audit(positive_control=("WT", "GMP"), null_group="WT")
if rep.failed:
    for r in rep.failed:
        print(r.name, r.summary, r.detail)
```

Every result carries a `detail` dict with the underlying numbers, so a check can
be turned into a supplementary table. `ds.report("report.html")` renders the
whole thing, and the CLI (`epimux audit --config study.yaml`) exits non-zero
when anything failed — which makes it usable in CI for an analysis pipeline.

---

## A worked example

From a real study, after the audit was run properly:

```
[PASS] check_direction[ATAC]       sign verified against raw values (corr +0.994)
[PASS] check_direction[H3K27ac]    sign verified against raw values (corr +1.000)
[PASS] positive_control[ATAC]      detects 44,834/106,839 (42.0%) -- pipeline has power
[PASS] null_contrast[ATAC:WT]      only 0 hits within-group
[WARN] null_contrast[H3K27ac:WT]   100 hits in a null contrast; interpret modest effects with care
[PASS] efficiency_balance[H3K27ac] imbalance 1.26x biases toward 'down', OPPOSITE to the
                                   observed 'up' effect -- the result is conservative
[PASS] replicate_reliability       group means are reliable (WT 0.94, KO 0.93)
[PASS] outlier_replicates          no replicate sits more than 0.10 below its group-mates
```

Read as a whole: the directions are verified, the pipeline demonstrably has
power, ATAC has essentially no false positives while H3K27ac is noisier, no
single library is responsible, and the technical bias makes the estimate
conservative rather than inflated. That is a result you can defend — and the
`WARN` belongs in the paper, not hidden.
