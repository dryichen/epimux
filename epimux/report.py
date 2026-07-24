"""Self-contained HTML report: audit traffic-lights, per-assay results, coupling."""
from __future__ import annotations

import base64
import datetime
import io
import os

import numpy as np
import pandas as pd

__all__ = ["html_report"]

_CSS = """
:root{--bg:#fff;--fg:#1A2433;--mut:#5A5A5A;--line:#e6e6e6;
      --pass:#00A087;--warn:#E6A817;--fail:#E64B35;--navy:#3C5488}
@media (prefers-color-scheme:dark){:root{--bg:#14171c;--fg:#e8eaed;--mut:#9aa0a6;--line:#2a2f37}}
*{box-sizing:border-box}
body{margin:0;padding:2rem 1.25rem;background:var(--bg);color:var(--fg);
     font:15px/1.6 -apple-system,'Helvetica Neue',Helvetica,Arial,sans-serif}
.wrap{max-width:1040px;margin:0 auto}
h1{font-size:1.7rem;margin:0 0 .2rem}h2{font-size:1.15rem;margin:2rem 0 .6rem;
   padding-bottom:.3rem;border-bottom:1px solid var(--line)}
.sub{color:var(--mut);margin:0 0 1.5rem;font-size:.92rem}
table{border-collapse:collapse;width:100%;font-size:.88rem;margin:.5rem 0}
th,td{text-align:left;padding:.45rem .6rem;border-bottom:1px solid var(--line)}
th{color:var(--mut);font-weight:600}
td.num{text-align:right;font-variant-numeric:tabular-nums}
.pill{display:inline-block;padding:.1rem .55rem;border-radius:999px;
      font-size:.75rem;font-weight:700;color:#fff}
.pass{background:var(--pass)}.warn{background:var(--warn)}.fail{background:var(--fail)}
.card{border:1px solid var(--line);border-radius:10px;padding:1rem;margin:.6rem 0}
.note{color:var(--mut);font-size:.85rem}
.scroll{overflow-x:auto}
code{background:rgba(128,128,128,.14);padding:.1rem .35rem;border-radius:4px;font-size:.85em}
"""


def _fig_b64(fig) -> str:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=140, bbox_inches="tight", facecolor="white")
    buf.seek(0)
    return base64.b64encode(buf.read()).decode()


def html_report(ds, path="epimux_report.html", fc=1.5, fdr=0.1, figures=True):
    """Render a Dataset into a single self-contained HTML file."""
    from . import plotting as pl
    import matplotlib.pyplot as plt

    parts = []
    A = parts.append
    A(f"<h1>epimux report &mdash; {ds.name}</h1>")
    A(f"<p class='sub'>{len(ds.reference):,} reference elements &middot; genome {ds.genome} "
      f"&middot; generated {datetime.datetime.now():%Y-%m-%d %H:%M}</p>")

    # ---- contrast ----
    if ds._contrast:
        A("<div class='card'><b>Contrast</b><br>")
        A(f"<code>{ds._contrast}</code><br>")
        A("<span class='note'>All effect sizes are log2(test / ref). "
          "The direction is pinned by the contrast object and verified against raw values.</span></div>")

    # ---- assays ----
    A("<h2>Assays</h2><div class='scroll'>")
    A(ds.summary().to_html(index=False, border=0))
    A("</div>")

    # ---- audit ----
    rep = ds.audit_report
    if rep.results:
        A("<h2>Audit</h2><div class='scroll'><table><tr><th>check</th><th>status</th><th>summary</th></tr>")
        for r in rep.results:
            A(f"<tr><td>{r.name}</td><td><span class='pill {r.status}'>{r.symbol}</span></td>"
              f"<td>{r.summary}</td></tr>")
        A("</table></div>")
        if rep.failed:
            A(f"<p class='note' style='color:var(--fail)'><b>{len(rep.failed)} check(s) failed — "
              "results should not be reported until these are resolved.</b></p>")
        if figures:
            try:
                fig, _ = pl.audit_plot(rep)
                A(f"<img style='max-width:100%' src='data:image/png;base64,{_fig_b64(fig)}'>")
            except Exception:
                pass

    # ---- results ----
    if ds.results:
        A("<h2>Differential results</h2><div class='scroll'><table>")
        A("<tr><th>assay</th><th class='num'>tested</th><th class='num'>significant</th>"
          "<th class='num'>up</th><th class='num'>down</th><th class='num'>% up</th></tr>")
        for n, r in ds.results.items():
            thr = np.log2(fc) if r.attrs.get("value_kind") != "difference" else 0.1
            sig = (r["padj"] < fdr) & (r["log2FC"].abs() > thr)
            up = int((sig & (r["log2FC"] > 0)).sum())
            dn = int((sig & (r["log2FC"] < 0)).sum())
            tot = up + dn
            A(f"<tr><td>{n}</td><td class='num'>{int(r['log2FC'].notna().sum()):,}</td>"
              f"<td class='num'>{tot:,}</td><td class='num'>{up:,}</td>"
              f"<td class='num'>{dn:,}</td>"
              f"<td class='num'>{(100*up/tot if tot else float('nan')):.0f}%</td></tr>")
        A("</table></div>")
        if figures:
            for n, r in ds.results.items():
                try:
                    fig, _ = pl.ma_plot(r, fc=fc, fdr=fdr, title=n)
                    A(f"<img style='max-width:48%' src='data:image/png;base64,{_fig_b64(fig)}'>")
                except Exception:
                    pass

    # ---- coupling ----
    names = list(ds.results)
    if len(names) >= 2:
        A("<h2>Cross-layer coupling</h2>")
        for i in range(len(names)):
            for j in range(i + 1, len(names)):
                try:
                    c = ds.coupling(names[i], names[j], fc=fc, fdr=fdr)
                except Exception as e:
                    A(f"<p class='note'>{names[i]} vs {names[j]}: {e}</p>")
                    continue
                A("<div class='card'>")
                A(f"<b>{c.assay_x} vs {c.assay_y}</b> &mdash; {c.interpretation}<br>")
                A(f"<span class='note'>n={c.n:,} &middot; P={c.pvalue:.2e} &middot; "
                  f"dual-significant {c.dual_significant} "
                  f"(same {c.same_direction}, opposite {c.opposite_direction})</span>")
                for note in c.notes:
                    A(f"<br><span class='note'>! {note}</span>")
                A("</div>")
                if figures:
                    try:
                        fig, _ = pl.coupling_plot(ds.results[names[i]], ds.results[names[j]],
                                                  names[i], names[j], fc=fc, fdr=fdr)
                        A(f"<img style='max-width:48%' src='data:image/png;base64,{_fig_b64(fig)}'>")
                    except Exception:
                        pass

    html = (f"<!doctype html><html><head><meta charset='utf-8'>"
            f"<meta name='viewport' content='width=device-width,initial-scale=1'>"
            f"<title>epimux &mdash; {ds.name}</title><style>{_CSS}</style></head>"
            f"<body><div class='wrap'>{''.join(parts)}</div></body></html>")
    with open(path, "w") as fh:
        fh.write(html)
    return os.path.abspath(path)
