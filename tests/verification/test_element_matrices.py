"""Element matrices against direct numerical integration."""

import numpy as np
import pytest

from ebbef2p import elements as el

l = 1.7
GX, GW = np.polynomial.legendre.leggauss(6)


def integrate(f):
    xi = 0.5*(GX + 1)
    return 0.5*l*sum(w*f(x) for x, w in zip(xi, GW, strict=True))


def test_winkler_matrix():
    k0, k1 = 3.0, 11.0
    ref = integrate(lambda xi: (k0 + (k1 - k0)*xi)*np.outer(el.shape_functions(xi, l),
                                                            el.shape_functions(xi, l)))
    assert el.winkler_stiffness(k0, k1, l) == pytest.approx(ref)


def test_shear_matrix():
    g0, g1 = 5.0, 2.0
    ref = integrate(lambda xi: (g0 + (g1 - g0)*xi)*np.outer(
        el.shape_function_derivatives(xi, l, 1), el.shape_function_derivatives(xi, l, 1)))
    assert el.shear_stiffness(g0, g1, l) == pytest.approx(ref)


def test_bending_matrix():
    EI = 7.0
    ref = integrate(lambda xi: EI*np.outer(el.shape_function_derivatives(xi, l, 2),
                                           el.shape_function_derivatives(xi, l, 2)))
    assert el.bending_stiffness(EI, l) == pytest.approx(ref)


def test_load_vector():
    q0, q1 = 4.0, -1.0
    ref = integrate(lambda xi: (q0 + (q1 - q0)*xi)*el.shape_functions(xi, l))
    assert el.distributed_load_vector(q0, q1, l) == pytest.approx(ref)


def test_derivatives_are_consistent():
    h = 1e-6
    for xi in (0.2, 0.7):
        num = (el.shape_functions(xi + h, l) - el.shape_functions(xi - h, l))/(2*h*l)
        assert el.shape_function_derivatives(xi, l, 1) == pytest.approx(num, rel=1e-6)
