"""Assembly and solution of the finite element system."""

import warnings

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from . import elements as el
from .mesh import build_mesh, end_spring_stiffness
from .results import Solution


def element_matrices(mesh):
    """Total element stiffness matrices and consistent load vectors.

    Returns arrays of shape (n, 4, 4) and (n, 4).
    """
    n = mesh.n_elements
    ke = np.empty((n, 4, 4))
    fe = np.empty((n, 4))
    for e, l in enumerate(mesh.lengths):
        ke[e] = (el.bending_stiffness(mesh.EI[e], l)
                 + el.winkler_stiffness(*mesh.k[e], l)
                 + el.shear_stiffness(*mesh.g[e], l))
        fe[e] = el.distributed_load_vector(*mesh.q[e], l)
    return ke, fe


def assemble(model, mesh):
    """Global stiffness matrix (CSR), load vector and element contributions."""
    n = mesh.n_elements
    ndof = 2*(n + 1)
    ke, fe = element_matrices(mesh)

    dofs = 2*np.arange(n)[:, None] + np.arange(4)[None, :]
    rows = np.repeat(dofs, 4, axis=1).ravel()
    cols = np.tile(dofs, (1, 4)).ravel()
    K = sp.coo_matrix((ke.ravel(), (rows, cols)), shape=(ndof, ndof)).tolil()

    f = np.zeros(ndof)
    np.add.at(f, dofs.ravel(), fe.ravel())

    if model.soil_beyond_ends:
        k_left, k_right = end_spring_stiffness(mesh)
        K[0, 0] += k_left
        K[ndof - 2, ndof - 2] += k_right

    for p in model.point_loads:
        i = mesh.node_index(p.x)
        f[2*i] += p.force
        f[2*i + 1] += p.moment

    return K.tocsr(), f, ke, fe


def constraints(model, mesh):
    """Prescribed degrees of freedom and their values."""
    prescribed = {}
    for s in model.supports:
        i = mesh.node_index(s.x)
        for offset, value in ((0, s.w), (1, s.theta)):
            if value is None:
                continue
            dof = 2*i + offset
            if dof in prescribed and prescribed[dof] != value:
                raise ValueError(f"conflicting supports at x = {s.x}")
            prescribed[dof] = float(value)
    dofs = np.array(sorted(prescribed), dtype=int)
    values = np.array([prescribed[d] for d in dofs], dtype=float)
    return dofs, values


def check_stability(model, mesh, fixed):
    """Raise if the beam can move as a rigid body.

    A foundation with k > 0 over a finite length restrains both
    translation and rotation. Without it, the supports must prescribe w at
    two distinct points, or w and theta at the same point.
    """
    if np.any(mesh.k > 0):
        return
    w_nodes = {d // 2 for d in fixed if d % 2 == 0}
    theta_nodes = {d // 2 for d in fixed if d % 2 == 1}
    if len(w_nodes) >= 2 or (w_nodes and theta_nodes):
        return
    raise np.linalg.LinAlgError(
        "the beam is not sufficiently supported: without a Winkler "
        "parameter k > 0 it needs w prescribed at two points, or w and "
        "theta at one point")


def solve(model, n_elements=100):
    """Solve the static problem.

    Parameters
    ----------
    model : Model
    n_elements : int
        Approximate number of elements along the beam.

    Returns
    -------
    Solution
    """
    mesh = build_mesh(model, n_elements)
    K, f, ke, fe = assemble(model, mesh)
    ndof = K.shape[0]

    fixed, values = constraints(model, mesh)
    check_stability(model, mesh, fixed)
    free = np.setdiff1d(np.arange(ndof), fixed)

    u = np.zeros(ndof)
    u[fixed] = values
    rhs = f[free] - K[free][:, fixed] @ values
    with warnings.catch_warnings():
        warnings.simplefilter("error", spla.MatrixRankWarning)
        try:
            u[free] = spla.spsolve(K[free][:, free].tocsc(), rhs)
        except (spla.MatrixRankWarning, RuntimeError):
            u[free] = np.nan
    if not np.all(np.isfinite(u)):
        raise np.linalg.LinAlgError(
            "the system is singular: the beam is not sufficiently supported "
            "(add a foundation or supports)")

    return Solution(model=model, mesh=mesh, u=u, K=K, f=f,
                    element_stiffness=ke, element_loads=fe,
                    fixed_dofs=fixed)
