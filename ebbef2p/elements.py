"""Element matrices for an Euler-Bernoulli beam on a two-parameter foundation.

Two-node element, degrees of freedom ``(w_i, theta_i, w_j, theta_j)`` with
``theta = dw/dx``, cubic Hermite interpolation. All foundation and load
quantities vary linearly along the element, so every matrix below is the
exact Galerkin integral for that variation.
"""

import numpy as np


def shape_functions(xi, l):
    """Hermite shape functions at the local coordinate ``xi`` in [0, 1]."""
    xi = np.asarray(xi, dtype=float)
    return np.array([
        1 - 3*xi**2 + 2*xi**3,
        l*xi*(1 - xi)**2,
        3*xi**2 - 2*xi**3,
        l*(xi**3 - xi**2),
    ])


def shape_function_derivatives(xi, l, order):
    """Derivatives of the shape functions with respect to ``x``."""
    xi = np.asarray(xi, dtype=float)
    if order == 1:
        return np.array([
            (-6*xi + 6*xi**2)/l,
            1 - 4*xi + 3*xi**2,
            (6*xi - 6*xi**2)/l,
            3*xi**2 - 2*xi,
        ])
    if order == 2:
        return np.array([
            (-6 + 12*xi)/l**2,
            (-4 + 6*xi)/l,
            (6 - 12*xi)/l**2,
            (-2 + 6*xi)/l,
        ])
    if order == 3:
        one = np.ones_like(xi)
        return np.array([12*one/l**3, 6*one/l**2, -12*one/l**3, 6*one/l**2])
    raise ValueError("order must be 1, 2 or 3")


def bending_stiffness(EI, l):
    """Flexural stiffness matrix, int EI N''^T N'' dx."""
    return EI/l**3*np.array([
        [12, 6*l, -12, 6*l],
        [6*l, 4*l**2, -6*l, 2*l**2],
        [-12, -6*l, 12, -6*l],
        [6*l, 2*l**2, -6*l, 4*l**2],
    ])


def winkler_stiffness(k_i, k_j, l):
    """First-parameter foundation matrix, int k(x) N^T N dx, k linear."""
    a, b = k_i, k_j
    return np.array([
        [l*(10*a + 3*b)/35, l**2*(15*a + 7*b)/420,
         9*l*(a + b)/140, -l**2*(7*a + 6*b)/420],
        [l**2*(15*a + 7*b)/420, l**3*(5*a + 3*b)/840,
         l**2*(6*a + 7*b)/420, -l**3*(a + b)/280],
        [9*l*(a + b)/140, l**2*(6*a + 7*b)/420,
         l*(3*a + 10*b)/35, -l**2*(7*a + 15*b)/420],
        [-l**2*(7*a + 6*b)/420, -l**3*(a + b)/280,
         -l**2*(7*a + 15*b)/420, l**3*(3*a + 5*b)/840],
    ])


def shear_stiffness(g_i, g_j, l):
    """Second-parameter foundation matrix, int g(x) N'^T N' dx, g linear.

    ``g`` is the shear parameter of the foundation: 2t in Vlasov's notation,
    G_p in Pasternak's.
    """
    a, b = g_i, g_j
    return np.array([
        [3*(a + b)/(5*l), b/10, -3*(a + b)/(5*l), a/10],
        [b/10, l*(3*a + b)/30, -b/10, -l*(a + b)/60],
        [-3*(a + b)/(5*l), -b/10, 3*(a + b)/(5*l), -a/10],
        [a/10, -l*(a + b)/60, -a/10, l*(a + 3*b)/30],
    ])


def distributed_load_vector(q_i, q_j, l):
    """Consistent nodal loads, int q(x) N dx, q linear."""
    return np.array([
        l*(7*q_i + 3*q_j)/20,
        l**2*(3*q_i + 2*q_j)/60,
        l*(3*q_i + 7*q_j)/20,
        -l**2*(2*q_i + 3*q_j)/60,
    ])
