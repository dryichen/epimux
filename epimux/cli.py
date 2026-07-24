"""Command-line interface.

    epimux audit  --config study.yaml
    epimux run    --config study.yaml --out results/

A study is described by one YAML/JSON file so an analysis is reproducible from
the shell and diffable in version control::

    name: STAG2 LSK
    genome: mm10
    reference: enhancers.bed
    assays:
      ATAC:    {type: counts, path: atac_counts.txt}
      H3K27ac: {type: counts, path: h3k_counts.txt}
      WGBS:    {type: methyl, path: "wgbs/*.cov"}
    design:
      WT:  [LSK_WT_R1, LSK_WT_R2, LSK_WT_R3]
      KO:  [LSK_KO_R1, LSK_KO_R2, LSK_KO_R3]
      GMP: [GMP_WT_R1, GMP_WT_R2, GMP_WT_R3]
    contrast: {ref: WT, test: KO}
    audit:
      positive_control: [WT, GMP]
      null_group: WT
    thresholds: {fc: 1.5, fdr: 0.1}
"""
from __future__ import annotations

import argparse
import json
import os
import sys

from .utils import get_logger

LOG = get_logger()


def _load_config(path):
    with open(path) as fh:
        text = fh.read()
    if path.endswith((".yaml", ".yml")):
        try:
            import yaml
        except ImportError:
            sys.exit("PyYAML needed for YAML configs; use JSON or `pip install pyyaml`")
        return yaml.safe_load(text)
    return json.loads(text)


def build_dataset(cfg):
    from .core import Dataset
    ds = Dataset(cfg["reference"], genome=cfg.get("genome", "unknown"),
                 name=cfg.get("name", "study"))
    for name, spec in cfg.get("assays", {}).items():
        t = spec.get("type", "counts")
        if t == "counts":
            ds.add_counts(name, spec["path"])
        elif t == "methyl":
            ds.add_methyl(name, spec["path"])
        elif t == "signal":
            ds.add_signal(name, spec["files"])
        elif t == "hic":
            ds.add_hic(name, spec["files"])
        else:
            LOG.warning(f"unknown assay type '{t}' for {name}, skipped")
    ds.set_design(cfg["design"])
    return ds


def cmd_run(args):
    cfg = _load_config(args.config)
    thr = cfg.get("thresholds", {})
    fc, fdr = thr.get("fc", 1.5), thr.get("fdr", 0.1)
    out = args.out or "epimux_out"
    os.makedirs(out, exist_ok=True)

    ds = build_dataset(cfg)
    ctr = cfg.get("contrast", {})
    ds.differential(ref=ctr["ref"], test=ctr["test"], strict=not args.no_strict)

    au = cfg.get("audit", {})
    ds.audit(positive_control=tuple(au["positive_control"]) if au.get("positive_control") else None,
             null_group=au.get("null_group"))

    for name, r in ds.results.items():
        r.to_csv(os.path.join(out, f"{name}_differential.tsv"), sep="\t")
    ds.audit_report.to_frame().to_csv(os.path.join(out, "audit.tsv"), sep="\t", index=False)

    names = list(ds.results)
    rows = []
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            try:
                c = ds.coupling(names[i], names[j], fc=fc, fdr=fdr)
                rows.append({"x": c.assay_x, "y": c.assay_y, "spearman": c.spearman,
                             "pvalue": c.pvalue, "n": c.n,
                             "dual_significant": c.dual_significant,
                             "same_direction": c.same_direction,
                             "opposite_direction": c.opposite_direction})
                print(repr(c))
            except Exception as e:
                LOG.warning(f"coupling {names[i]}/{names[j]}: {e}")
    if rows:
        import pandas as pd
        pd.DataFrame(rows).to_csv(os.path.join(out, "coupling.tsv"), sep="\t", index=False)

    p = ds.report(os.path.join(out, "report.html"), fc=fc, fdr=fdr)
    print(f"\nreport: {p}")
    return 1 if ds.audit_report.failed else 0


def cmd_audit(args):
    cfg = _load_config(args.config)
    ds = build_dataset(cfg)
    ctr = cfg.get("contrast", {})
    ds.differential(ref=ctr["ref"], test=ctr["test"], strict=False)
    au = cfg.get("audit", {})
    rep = ds.audit(
        positive_control=tuple(au["positive_control"]) if au.get("positive_control") else None,
        null_group=au.get("null_group"))
    print("\n" + repr(rep))
    return 1 if rep.failed else 0


def main(argv=None):
    p = argparse.ArgumentParser(prog="epimux", description=__doc__.split("\n")[0])
    p.add_argument("--version", action="store_true")
    sub = p.add_subparsers(dest="cmd")

    r = sub.add_parser("run", help="full analysis + report")
    r.add_argument("--config", required=True)
    r.add_argument("--out")
    r.add_argument("--no-strict", action="store_true",
                   help="do not abort when a direction check fails")
    r.set_defaults(fn=cmd_run)

    a = sub.add_parser("audit", help="run only the audit battery")
    a.add_argument("--config", required=True)
    a.set_defaults(fn=cmd_audit)

    args = p.parse_args(argv)
    if args.version:
        from . import __version__
        print(f"epimux {__version__}")
        return 0
    if not getattr(args, "fn", None):
        p.print_help()
        return 0
    return args.fn(args)


if __name__ == "__main__":
    sys.exit(main())
