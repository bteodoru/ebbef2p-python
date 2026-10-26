"""Element routines evaluated for many elements at once."""

import numpy as np
import pytest

from ebbef2p import elements as el

RNG = np.random.default_rng(0)
N = 7
L = RNG.uniform(0.1, 2.0, N)
A, B = RNG.uniform(0, 1e4, (2, N))


@pytest.mark.parametrize("matrix, args", [
    (el.bending_stiffness, (A, L)),
    (el.winkler_stiffness, (A, B, L)),
    (el.shear_stiffness, (A, B, L)),
    (el.distributed_load_vector, (A, B, L)),
])
def test_stacked_matrices_equal_element_by_element(matrix, args):
    stacked = matrix(*args)
    single = np.stack([matrix(*(a[e] for a in args)) for e in range(N)])
    assert stacked.shape == single.shape
    assert stacked == pytest.approx(single, rel=1e-14)


@pytest.mark.parametrize("order", [0, 1, 2, 3])
def test_shape_functions_broadcast_over_points_and_elements(order):
    xi = np.linspace(0, 1, 5)

    def f(x, l):
        if order == 0:
            return el.shape_functions(x, l)
        return el.shape_function_derivatives(x, l, order)

    grid = f(xi, L[:, None])
    assert grid.shape == (4, N, 5)
    for e in range(N):
        assert grid[:, e, :] == pytest.approx(f(xi, L[e]))


def test_derivative_order_is_checked():
    with pytest.raises(ValueError, match="order"):
        el.shape_function_derivatives(0.5, 1.0, 4)
