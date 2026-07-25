"""The tutorial must stay executable.

It was shipped broken once: the notebook was written with `source` lines that
had no trailing newlines, so nbformat joined every line into one and nothing
could run. Nothing caught it because no test ever executed the notebook.
"""
import json
import os
import subprocess
import sys

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
NB = os.path.join(ROOT, "docs", "tutorial.ipynb")


def _load():
    with open(NB) as fh:
        return json.load(fh)


def test_notebook_is_wellformed():
    nb = _load()
    assert nb["nbformat"] == 4
    assert nb["cells"], "no cells"
    for i, c in enumerate(nb["cells"]):
        assert c["cell_type"] in ("code", "markdown"), i
        src = c["source"]
        assert isinstance(src, list), f"cell {i}: source must be a list of lines"
        # every line except the last must end in a newline, or the notebook
        # collapses into a single unrunnable line
        for line in src[:-1]:
            assert line.endswith("\n"), f"cell {i}: line without newline: {line[:40]!r}"


def test_notebook_runs(tmp_path):
    """Execute every code cell in order and require a clean exit."""
    nb = _load()
    code = "\n".join("".join(c["source"]) for c in nb["cells"] if c["cell_type"] == "code")
    script = tmp_path / "tutorial.py"
    script.write_text(code)
    env = dict(os.environ, PYTHONPATH=ROOT, MPLBACKEND="Agg")
    r = subprocess.run([sys.executable, str(script)], cwd=tmp_path,
                       capture_output=True, text=True, env=env, timeout=900)
    assert r.returncode == 0, f"tutorial failed:\n{r.stdout[-3000:]}\n{r.stderr[-3000:]}"


def test_quickstart_runs(tmp_path):
    env = dict(os.environ, PYTHONPATH=ROOT, MPLBACKEND="Agg")
    r = subprocess.run([sys.executable, os.path.join(ROOT, "examples", "quickstart.py")],
                       cwd=tmp_path, capture_output=True, text=True, env=env, timeout=900)
    assert r.returncode == 0, f"quickstart failed:\n{r.stdout[-3000:]}\n{r.stderr[-3000:]}"
    assert "COUPLED" in r.stdout


def test_positive_control_is_not_circular():
    """A positive control must be a DIFFERENT comparison from the contrast."""
    nb = _load()
    code = "\n".join("".join(c["source"]) for c in nb["cells"] if c["cell_type"] == "code")
    qs = open(os.path.join(ROOT, "examples", "quickstart.py")).read()
    for name, src in (("tutorial", code), ("quickstart", qs)):
        assert 'positive_control=("WT", "KO")' not in src, (
            f"{name}: uses the contrast of interest as its own positive control")


def test_documented_api_exists():
    """Every epimux name referenced in the docs must actually exist."""
    import re
    import epimux as ep
    names = set()
    for fn in ["README.md"] + [f"docs/{f}" for f in os.listdir(os.path.join(ROOT, "docs"))
                               if f.endswith(".md") and f != "api.md"]:
        text = open(os.path.join(ROOT, fn)).read()
        names |= set(re.findall(r"\bep\.([a-z_][a-z0-9_]*)\s*\(", text))
    missing = sorted(n for n in names if not hasattr(ep, n))
    assert not missing, f"documented but missing from epimux: {missing}"
