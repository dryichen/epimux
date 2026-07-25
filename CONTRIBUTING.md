# Contributing

## Development setup

```bash
git clone https://github.com/dryichen/epimux
cd epimux
pip install -e ".[all,dev]"
pytest -q
```

## Ground rules

**Directions are pinned, never inferred.** Any new differential engine must take
a `Contrast` and return `log2FC = log2(test / ref)`. If you add one, also add a
`check_direction` test that feeds it a reversed contrast and asserts `FAIL`.

**Every audit check needs a failing case.** A check that cannot fail on
constructed bad input is not a check. See `tests/test_core.py` for the pattern.

**Backgrounds are explicit.** Enrichment functions take a background and must
not invent one from the whole genome.

**Warn on weak methods.** If a code path is less sensitive than the alternative
(averaged signal vs read counts, nearest-gene vs contact-based linking), say so
in the docstring and at runtime.

## Tests

New code needs tests with a known ground truth — plant an effect, then assert it
is recovered with the right sign and magnitude. Synthetic fixtures live at the
top of each test module.

## Style

Follow the surrounding code: NumPy-style docstrings that explain *why*, type
hints on public signatures, no abbreviations in public names.
