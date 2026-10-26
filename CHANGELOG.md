# Changelog

All notable changes to this project are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the project uses
[semantic versioning](https://semver.org/).

## [1.0.0] - unreleased

A rewrite with a new API. Code written for 0.0.x does not run unchanged;
see *Migrating from 0.0.x* below.

### Added

- `Model` with chainable builders (`add_beam`, `add_foundation`,
  `add_point_load`, `add_distributed_load`, `add_support`) and validation of
  the input when it is added and before solving.
- Both foundation parameters may vary linearly along every segment;
  overlapping segments add up, so jumps in `k` and `g` are allowed.
- `Model(soil_beyond_ends=True)`: the foundation continues beyond the beam
  ends, represented by equivalent springs `sqrt(k*g)`.
- `Solution` with nodal `w` and `theta`, `deflection(x)` from the Hermite
  interpolation, `bending_moment()` and `shear_force()` from element
  equilibrium (values on both sides of every node, so jumps are kept),
  `reactive_pressure()` including the concentrated forces of the shear
  layer, `support_reactions()` and an `equilibrium()` check.
- `VlasovSoil` (Young's modulus varying linearly with depth) and
  `solve_vlasov()` with the iteration for the decay parameter gamma.
- Package exceptions `ModelError`, `UnstableModelError`,
  `ConvergenceError`, all derived from `EBBEF2pError` and from the
  matching built-in exception.
- Type hints (`py.typed`).
- Test suite split into `tests/unit` and `tests/verification`; the
  verification tests compare with closed-form beam solutions, Hetényi's
  long beam, the exact results of Razaqpur & Shah (1991), the gamma values
  of Vallabhan & Das (1991), numerical integration of the element matrices
  and Vlasov parameters, and global equilibrium.

### Fixed

- Bending moments and shear forces were computed from the bending
  stiffness alone, ignoring the foundation and the distributed load on the
  element. The shear force did not converge with mesh refinement: it
  contained an error of the order of the shear-layer force `g*w'`.
- Two concentrated loads at the same node: only the last one was applied.
- The gamma update of the Vlasov model included the soil beyond the beam
  ends while the stiffness matrix had no end springs. For the uniform-load
  case of Vallabhan & Das this gave gamma = 0.498 instead of 0.395.
- The Vlasov parameters lost accuracy for small gamma because of
  cancellation (error in `g` of 1e-5 at gamma = 1e-3, and `g = 0` instead
  of `E*H/(6(1+nu))` for a rigidly settling beam). A Taylor series is now
  used for gamma < 0.3.
- Compatibility with NumPy 2 (`np.trapz` and `np.asmatrix` were removed).

### Changed

- The second foundation parameter is called `g` (`g = 2t` in Vlasov's
  notation, `G_p` in Pasternak's). The old code called it `t` although it
  meant `2t`.
- Sparse assembly and solution (SciPy); element computations are
  vectorised.
- Matplotlib is optional (`pip install "ebbef2p-python[plot]"`).
- Requires Python 3.12 or newer, NumPy 1.26 or newer and SciPy 1.11.3 or
  newer.
- Packaging with `pyproject.toml`; the distribution is `ebbef2p-python`
  (`setup.py` called it `ebbef2py`).

### Removed

- The `Structure` class and the modules `beam`, `element`,
  `elastic_foundation`, `distributed_load`, `nodal_load`, `nodal_support`,
  `vlasov_foundation_parameters` and `helpers`.

### Migrating from 0.0.x

| 0.0.x | 1.0 |
|---|---|
| `s = Structure("name")` | `m = Model()` |
| `s.add_beam([x0, x1], E, I)` | `m.add_beam(x0, x1, E, I)` |
| `s.add_elastic_foundation([k0, k1], [x0, x1], 'k')` | `m.add_foundation(x0, x1, k=(k0, k1))` |
| `s.add_elastic_foundation([t0, t1], [x0, x1], 't')` | `m.add_foundation(x0, x1, g=(t0, t1))` (the old `t` was already `2t`) |
| `s.add_nodal_load(P, x, 'fz')` / `(M, x, 'my')` | `m.add_point_load(x, force=P, moment=M)` |
| `s.add_distributed_load([q0, q1], [x0, x1])` | `m.add_distributed_load(x0, x1, q0, q1)` |
| `s.add_nodal_support({'uz': 0, 'ur': 'NaN'}, x)` | `m.add_support(x, w=0, theta=None)` |
| `add_nodes`, `add_elements`, `build_global_matrix`, `solve(...)` | `solution = solve(m, n_elements)` |
| `get_displacements()` (sign flipped) | `solution.w`, `solution.theta` (positive downwards) |
| `get_bending_moments()`, `get_shear_forces()` | `solution.bending_moment()`, `solution.shear_force()` |
| `compute_elastic_foundation_parameters({'E1':…, 'E2':…}, nu, H)` | `solve_vlasov(m, VlasovSoil(E_top, nu, H, E_bottom))` |

## [0.0.3.dev0] - 2020-05-25

- Foundation input reworked: `soil_conditions` became `elastic_foundation`,
  `elastic_foundation_support` was removed; `Structure` and the Vlasov
  parameters revised; the iteration for gamma
  (`Structure.compute_elastic_foundation_parameters`) and `example.py`
  added.
- Tests for the simply supported beam, the long beam on a Winkler
  foundation, a beam on a two-parameter foundation and the Vlasov
  parameters (replacing the first tests).

## [0.0.2.dev1] - 2020-04-22

- First tests (`Beam`, `Structure`), Travis CI with coverage,
  `requirements.txt`, documentation updates.

## [0.0.1.dev1] - 2020-04-21

First Python version: the `Structure` class, beams on Winkler and
two-parameter foundations, and the modified Vlasov formulas for `k`, `t`
and gamma.

The 0.0.x versions were tagged retroactively in 2026 on the commits whose
code matches the sdists published on PyPI; the version number in
`setup.py` was not always committed at the time.

[1.0.0]: https://github.com/bteodoru/ebbef2p-python/compare/v0.0.3.dev0...HEAD
[0.0.3.dev0]: https://github.com/bteodoru/ebbef2p-python/compare/v0.0.2.dev1...v0.0.3.dev0
[0.0.2.dev1]: https://github.com/bteodoru/ebbef2p-python/compare/v0.0.1.dev1...v0.0.2.dev1
[0.0.1.dev1]: https://github.com/bteodoru/ebbef2p-python/releases/tag/v0.0.1.dev1
