# Command line

A study is described by one config file, so an analysis is reproducible from the
shell and diffable in version control.

```bash
epimux audit --config study.yaml          # checks only; non-zero exit on failure
epimux run   --config study.yaml --out results/
epimux --version
```

`epimux run` writes per-assay tables, `audit.tsv`, `coupling.tsv`, a
`manifest.json` and a self-contained `report.html`.

**Exit codes.** Both commands return `1` when any audit check failed, which
makes them usable as a gate in a pipeline or in CI:

```bash
epimux audit --config study.yaml || { echo "audit failed"; exit 1; }
```

---

## Config reference

```yaml
name: STAG2 LSK                 # free text, appears in the report
genome: mm10                    # free text
reference: data/enhancers.bed   # BED-like; the shared element set

assays:
  ATAC:
    type: counts                # counts | methyl | signal | hic
    path: data/atac_counts.txt  # featureCounts output
  H3K27ac:
    type: counts
    path: data/h3k_counts.txt
  WGBS:
    type: methyl
    path: "data/wgbs/*.cov"     # glob of Bismark .cov files
  HiC:
    type: hic
    files:
      WT: data/wt.mcool
      KO: data/ko.mcool

design:                         # group -> sample names, as they appear in the matrices
  WT:  [LSK_WT_R1, LSK_WT_R2, LSK_WT_R3]
  KO:  [LSK_KO_R1, LSK_KO_R2, LSK_KO_R3]
  GMP: [GMP_WT_R1, GMP_WT_R2, GMP_WT_R3]

contrast:
  ref: WT                       # log2FC = log2(test / ref)
  test: KO

audit:
  positive_control: [WT, GMP]   # a comparison that MUST differ
  null_group: WT                # replicates split against each other

thresholds:
  fc: 1.5                       # fold-change (not log2)
  fdr: 0.1
```

JSON is accepted too — the extension decides the parser.

**Flags.** `--out DIR` sets the output directory (default `epimux_out`).
`--no-strict` lets `run` continue when a direction check fails instead of
aborting; use it only to inspect a broken analysis, never to publish one.
