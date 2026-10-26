# Contributing to EBBEF2p

Bug reports, questions and pull requests are welcome on
[GitHub](https://github.com/bteodoru/ebbef2p-python/issues).

## Reporting a problem

Please include a minimal script that reproduces it, the result you expected
(with its source, if it comes from a publication or a hand calculation)
and the output of

```
python -c "import ebbef2p, numpy, scipy, sys; print(ebbef2p.__version__, numpy.__version__, scipy.__version__, sys.version)"
```

## Development setup

The project is managed with [uv](https://docs.astral.sh/uv/):

```
git clone https://github.com/bteodoru/ebbef2p-python.git
cd ebbef2p-python
uv sync
uv run pytest
uv run ruff check .
uv run sphinx-build -W -b html docs/source docs/build
```

## Guidelines

- **Tests.** Every change comes with tests. `tests/unit` holds tests of
  the API, input validation and internal consistency; `tests/verification`
  holds comparisons with independent references (closed-form solutions,
  published results, numerical integration, equilibrium). A new feature
  that computes something physical needs at least one verification test,
  with the reference cited in the test. Prefer exact values printed in a
  source over values read from a figure.
- **Conventions.** `w` is positive downwards, `theta = dw/dx`,
  `M = -EI w''`, `V = dM/dx` (beam shear, without the shear-layer force),
  and the second foundation parameter is `g` (`2t` for Vlasov, `G_p` for
  Pasternak). See `docs/source/user_guide/conventions.rst` before changing
  anything that touches signs.
- **Style.** `uv run ruff check .` must pass. Type hints on public
  functions, NumPy-style docstrings.
- **Documentation.** Update the docs and `CHANGELOG.md` (section
  *Unreleased*) for user-visible changes. The docs must build with `-W`.

## Releasing (maintainer)

1. Set the version in `src/ebbef2p/__init__.py` and in `CITATION.cff`
   (`version`, `date-released`), and move the *Unreleased* entries of
   `CHANGELOG.md` under the new version with its date.
2. Merge to `main` and wait for CI.
3. Publish a GitHub release with the tag `vX.Y.Z` on `main`. The `release`
   workflow checks that the tag matches the version, runs the tests and
   publishes to PyPI (trusted publishing); Zenodo archives the release and
   mints a DOI.

By contributing you agree that your contribution is licensed under the MIT
License of the project.
