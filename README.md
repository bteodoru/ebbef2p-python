
# EBBEF2p for Python

[![PyPI version](https://badge.fury.io/py/ebbef2p-python.svg)](https://badge.fury.io/py/ebbef2p-python)
[![GitHub license](https://img.shields.io/github/license/bteodoru/ebbef2p-python)](https://github.com/bteodoru/ebbef2p-python/blob/master/LICENSE)
[![Documentation Status](https://readthedocs.org/projects/ebbef2p-python/badge/?version=latest)](https://ebbef2p-python.readthedocs.io/en/latest/?badge=latest)

EBBEF2p (Euler-Bernoulli Beams on Two-Parameter Elastic Foundations) is a
finite element package for static analysis of beams on Winkler, Pasternak
and modified Vlasov foundations.

* both foundation parameters may vary linearly along the beam, with jumps
  between segments;
* the soil may continue beyond the beam ends (equivalent end springs);
* modified Vlasov model (Vallabhan & Das, 1991) with iterated decay
  parameter and depth-dependent soil modulus;
* bending moments and shear forces from element equilibrium;
* soil pressure including the concentrated forces of the shear layer, with
  a global equilibrium check.

EBBEF2p was originally written as a MATLAB toolbox.

## Installation

```
pip install ebbef2p-python
```

## Quick start

```python
from ebbef2p import Model, solve

model = Model(soil_beyond_ends=True)
model.add_beam(0, 10, E=30e6, I=0.01)
model.add_foundation(0, 10, k=5000, g=2000)   # g = 2t (Vlasov) = Gp (Pasternak)
model.add_point_load(5, force=100)

solution = solve(model, n_elements=100)
print(solution.w.max(), solution.bending_moment().max())
print(solution.equilibrium())
```

A modified Vlasov foundation is defined by the soil layer instead:

```python
from ebbef2p import Model, VlasovSoil, solve_vlasov

model = Model().add_beam(0, 100, E=432000, I=2.25).add_distributed_load(0, 100, 3)
result = solve_vlasov(model, VlasovSoil(E_top=144, nu=0.2, depth=25))
print(result.gamma, result.k, result.g)
```

More in the [documentation](https://ebbef2p-python.readthedocs.io) and in
the `examples` folder.

## Development

```
pip install -e .[plot]
pytest
```
