# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

EBBEF2p: finite element analysis of Euler-Bernoulli beams on two-parameter elastic
foundations (Winkler, Pasternak, modified Vlasov). The package backs an educational web app
(FastAPI + Vue, separate repo) and is being prepared for publication in SoftwareX, followed by an
education study in *Computer Applications in Engineering Education*. The manuscript lives outside
this repo.

Communicate with the user in Romanian; code, comments, docstrings, docs and commit messages are in
English.

## Commands

The project is managed with uv (Python >= 3.12). `uv sync` creates `.venv` with the dev group.

```
uv run pytest                                   # all tests (unit + verification)
uv run pytest tests/verification                # only comparisons with references
uv run pytest tests/unit/test_model.py::test_model_validation   # single test
uv run pytest --cov --cov-report=term-missing
uv run ruff check .                             # lint (no formatter is used, see Style)
uv run sphinx-build -W -b html docs/source docs/build   # docs must build without warnings
uv build                                        # wheel + sdist
```

Minimum dependency versions (numpy 1.26, scipy 1.11.3 on Python 3.12) are checked with:

```
uv run --isolated --python 3.12 --resolution lowest-direct --no-default-groups --group test pytest -p no:cacheprovider
```

After changing dependencies in `pyproject.toml`, run `uv lock` (CI uses `--locked`).
pytest runs with `filterwarnings = error`, so any new warning fails the suite.

## Architecture

The pipeline is `Model` → `build_mesh` → `assemble` → sparse solve → `Solution`; every stage is a
separate module in `src/ebbef2p/` and data flows one way.

- `model.py` holds only the problem definition (frozen dataclasses for segments, loads, supports;
  chainable `add_*` builders returning `Self`). It does no numerics. Overlapping foundation
  segments add up; that is how jumps in `k` and `g` are expressed.
- `mesh.py` turns every key point (segment ends, loads, supports) into a node and stores per-element
  arrays: `EI (n,)`, `k, g, q (n, 2)` (values at both ends, linear in between). `Mesh.element_dofs`
  `(n, 4)` is the index map used everywhere.
- `elements.py` gives exact Galerkin matrices for linearly varying `k`, `g`, `q` with cubic Hermite
  shape functions. All functions broadcast: arrays of element properties give stacks
  `(n, 4, 4)`; shape functions broadcast `xi` against `l` and return the 4 functions on axis 0.
  Keep new code vectorised in the same way; no Python loops over elements.
- `solver.py` assembles a CSR matrix, adds end springs `sqrt(k*g)` when
  `model.soil_beyond_ends`, applies prescribed DOFs, checks rigid-body stability and solves.
- `results.py` (`Solution`) post-processes. `element_field(xi, order)` evaluates `w` or its
  derivatives at local points of every element and underlies pressure, deflection and the Vlasov
  update. Internal forces come from element equilibrium, `F = K_e d_e - f_e`, not from
  differentiating `w`. Quantities at element ends are returned as `ElementEndValues`
  (`2*n` values, left and right of every node) so that jumps survive.
- `vlasov.py` iterates: `gamma → (k, g)` from `VlasovSoil.parameters` → `solve` →
  `gamma_from_deflection` → repeat. It builds a copy of the model with a uniform foundation, so the
  input model must have no foundation of its own.
- `exceptions.py`: `ModelError`, `UnstableModelError`, `ConvergenceError` derive from
  `EBBEF2pError` and from `ValueError` / `LinAlgError` / `RuntimeError` respectively; keep both
  bases.

The public API is what `src/ebbef2p/__init__.py` re-exports.

## Conventions that must not change silently

Full text in `docs/source/user_guide/conventions.rst`.

- `w` positive downwards, `theta = dw/dx`; loads positive downwards; reactions positive upwards,
  so `applied = soil + supports` in `Solution.equilibrium()`.
- `M = -EI w''` (positive sagging). `V = dM/dx` is the shear in the beam alone and excludes the
  shear-layer force `g w'`. In element terms: `V_i = -F_1 - g_i theta_i`, `V_j = F_3 - g_j theta_j`.
  Razaqpur & Shah's "vertical shear" (`-EI w''' + k1 w'`) is a different quantity.
- `g` is the second parameter: `g = 2t` (Vlasov) `= G_p` (Pasternak). The pre-1.0 code called it
  `t` although it meant `2t`.
- Reactive pressure `p = k w - (g w')'` at element ends, plus concentrated forces: `-g w'` at the
  start, `+g w'` at the end, the jump of `g w'` where `g` jumps, and `sqrt(k g) w` at the ends when
  `soil_beyond_ends`.
- `soil_beyond_ends` is read only from `Model` and must enter both the stiffness matrix (end
  springs) and the gamma update (outside integrals). Treating them inconsistently was the main bug
  of the old code (gamma 0.498 instead of 0.395 for Vallabhan & Das, uniform load).
- `VlasovSoil.parameters` switches from the closed forms to their Taylor series below
  `gamma = 0.3` because the closed form of `g` cancels catastrophically as `gamma → 0`.

## Tests

- `tests/unit`: API behaviour, validation, internal consistency, examples, plotting.
- `tests/verification`: comparisons with independent references (closed-form beams, Hetényi,
  Razaqpur & Shah 1991, Vallabhan & Das 1991, quadrature, global equilibrium). This folder backs
  the verification table in `docs/source/theory/verification.rst`; keep the two in sync.
- Verification references must be exact values printed in a primary source or computed
  independently, never values read off a figure. Cite the source in the test.
- `tests/unit/test_examples.py` runs every script in `examples/`, so examples must print something
  and stay fast.

## Style

- ruff lint only. There is no auto-formatter: the code uses compact arithmetic (`EI/l**3`,
  `0.5*l`) close to the formulas, so do not run `ruff format`. Line length 100.
- `l` is the element length throughout (E741 is ignored on purpose).
- NumPy-style docstrings; type hints on public functions (`NDArray[np.float64]`, `Self`).
- User-visible changes go into `CHANGELOG.md`; the version lives in `src/ebbef2p/__init__.py`
  (read by hatchling).
