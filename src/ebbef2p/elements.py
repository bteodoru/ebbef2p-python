"""Element matrices for an Euler-Bernoulli beam on a two-parameter foundation.

Two-node element, degrees of freedom ``(w_i, theta_i, w_j, theta_j)`` with
``theta = dw/dx``, cubic Hermite interpolation. All foundation and load
quantities vary linearly along the element, so every matrix below is the
exact Galerkin integral for that variation.

The matrix functions accept scalars or arrays of element properties; for
arrays of shape ``(n,)`` they return a stack of shape ``(n, 4, 4)`` (or
``(n, 4)`` for load vectors).
"""

import numpy as np
from numpy.typing import ArrayLike, NDArray


def _stack(rows: list[list[ArrayLike]]) -> NDArray[np.float64]:
    """Matrix (or stack of matrices) from a nested list of broadcastable entries."""
    entries = np.broadcast_arrays(*(np.asarray(v, dtype=float) for row in rows for v in row))
    m = np.stack(entries, axis=-1)
    return m.reshape((*m.shape[:-1], len(rows), len(rows[0])))


def _rows(*entries: ArrayLike) -> NDArray[np.float64]:
    """Stack broadcastable entries along a new first axis."""
    return np.stack(np.broadcast_arrays(*(np.asarray(v, dtype=float) for v in entries)))


def shape_functions(xi: ArrayLike, l: ArrayLike) -> NDArray[np.float64]:
    """Hermite shape functions at the local coordinate ``xi`` in [0, 1].

    ``xi`` and ``l`` broadcast against each other; the result has the four
    functions along the first axis.
    """
    xi, l = np.asarray(xi, dtype=float), np.asarray(l, dtype=float)
    return _rows(
        1 - 3*xi**2 + 2*xi**3,
        l*xi*(1 - xi)**2,
        3*xi**2 - 2*xi**3,
        l*(xi**3 - xi**2),
    )


def shape_function_derivatives(xi: ArrayLike, l: ArrayLike, order: int) -> NDArray[np.float64]:
    """Derivatives of the shape functions with respect to ``x``."""
    xi, l = np.asarray(xi, dtype=float), np.asarray(l, dtype=float)
    if order == 1:
        return _rows(
            (-6*xi + 6*xi**2)/l,
            1 - 4*xi + 3*xi**2,
            (6*xi - 6*xi**2)/l,
            3*xi**2 - 2*xi,
        )
    if order == 2:
        return _rows(
            (-6 + 12*xi)/l**2,
            (-4 + 6*xi)/l,
            (6 - 12*xi)/l**2,
            (-2 + 6*xi)/l,
        )
    if order == 3:
        z = np.zeros_like(xi)
        return _rows(12/l**3 + z, 6/l**2 + z, -12/l**3 + z, 6/l**2 + z)
    raise ValueError("order must be 1, 2 or 3")


def bending_stiffness(EI: ArrayLike, l: ArrayLike) -> NDArray[np.float64]:
    """Flexural stiffness matrix, int EI N''^T N'' dx."""
    EI, l = np.asarray(EI, dtype=float), np.asarray(l, dtype=float)
    c = EI/l**3
    return _stack([
        [12*c, 6*l*c, -12*c, 6*l*c],
        [6*l*c, 4*l**2*c, -6*l*c, 2*l**2*c],
        [-12*c, -6*l*c, 12*c, -6*l*c],
        [6*l*c, 2*l**2*c, -6*l*c, 4*l**2*c],
    ])


def winkler_stiffness(k_i: ArrayLike, k_j: ArrayLike, l: ArrayLike) -> NDArray[np.float64]:
    """First-parameter foundation matrix, int k(x) N^T N dx, k linear."""
    a, b, l = (np.asarray(v, dtype=float) for v in (k_i, k_j, l))
    return _stack([
        [l*(10*a + 3*b)/35, l**2*(15*a + 7*b)/420,
         9*l*(a + b)/140, -l**2*(7*a + 6*b)/420],
        [l**2*(15*a + 7*b)/420, l**3*(5*a + 3*b)/840,
         l**2*(6*a + 7*b)/420, -l**3*(a + b)/280],
        [9*l*(a + b)/140, l**2*(6*a + 7*b)/420,
         l*(3*a + 10*b)/35, -l**2*(7*a + 15*b)/420],
        [-l**2*(7*a + 6*b)/420, -l**3*(a + b)/280,
         -l**2*(7*a + 15*b)/420, l**3*(3*a + 5*b)/840],
    ])


def shear_stiffness(g_i: ArrayLike, g_j: ArrayLike, l: ArrayLike) -> NDArray[np.float64]:
    """Second-parameter foundation matrix, int g(x) N'^T N' dx, g linear.

    ``g`` is the shear parameter of the foundation: 2t in Vlasov's notation,
    G_p in Pasternak's.
    """
    a, b, l = (np.asarray(v, dtype=float) for v in (g_i, g_j, l))
    return _stack([
        [3*(a + b)/(5*l), b/10, -3*(a + b)/(5*l), a/10],
        [b/10, l*(3*a + b)/30, -b/10, -l*(a + b)/60],
        [-3*(a + b)/(5*l), -b/10, 3*(a + b)/(5*l), -a/10],
        [a/10, -l*(a + b)/60, -a/10, l*(a + 3*b)/30],
    ])


def distributed_load_vector(q_i: ArrayLike, q_j: ArrayLike, l: ArrayLike) -> NDArray[np.float64]:
    """Consistent nodal loads, int q(x) N dx, q linear."""
    q_i, q_j, l = (np.asarray(v, dtype=float) for v in (q_i, q_j, l))
    return _stack([[
        l*(7*q_i + 3*q_j)/20,
        l**2*(3*q_i + 2*q_j)/60,
        l*(3*q_i + 7*q_j)/20,
        -l**2*(2*q_i + 3*q_j)/60,
    ]])[..., 0, :]
