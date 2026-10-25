"""Discretisation of a :class:`~ebbef2p.model.Model` into beam elements."""

from dataclasses import dataclass
from itertools import pairwise

import numpy as np
from numpy.typing import NDArray

from .exceptions import ModelError
from .model import Model


@dataclass(frozen=True)
class Mesh:
    """Nodes and per-element properties.

    Attributes
    ----------
    nodes : (n+1,) array
    EI : (n,) array
    k, g : (n, 2) arrays
        Foundation parameters at the start and end of every element.
    q : (n, 2) array
        Distributed load at the start and end of every element.
    """
    nodes: NDArray[np.float64]
    EI: NDArray[np.float64]
    k: NDArray[np.float64]
    g: NDArray[np.float64]
    q: NDArray[np.float64]

    @property
    def n_elements(self) -> int:
        return len(self.nodes) - 1

    @property
    def lengths(self) -> NDArray[np.float64]:
        return np.diff(self.nodes)

    @property
    def element_dofs(self) -> NDArray[np.intp]:
        """Global degrees of freedom of every element, shape (n, 4)."""
        return 2*np.arange(self.n_elements)[:, None] + np.arange(4)[None, :]

    def node_index(self, x: float, tol: float = 1e-9) -> int:
        """Index of the node at ``x`` (raises if there is none)."""
        span = self.nodes[-1] - self.nodes[0]
        i = int(np.argmin(np.abs(self.nodes - x)))
        if abs(self.nodes[i] - x) > tol*max(1.0, span):
            raise ModelError(f"no node at x = {x}")
        return i

    def element_at(self, x: float | NDArray[np.float64]) -> NDArray[np.intp]:
        """Index of the element containing ``x`` (right-continuous)."""
        i = np.searchsorted(self.nodes, x, side='right') - 1
        return np.clip(i, 0, self.n_elements - 1)


def end_spring_stiffness(mesh: Mesh) -> tuple[float, float]:
    """Springs representing the foundation beyond each beam end.

    For a two-parameter foundation extending to infinity, the surface
    deflection outside the beam decays as w_end*exp(-sqrt(k/g)*s), and that
    soil reacts on the beam end with a force sqrt(k*g)*w_end. Returns the
    spring stiffnesses at the start and end of the beam.
    """
    return (float(np.sqrt(mesh.k[0, 0]*mesh.g[0, 0])),
            float(np.sqrt(mesh.k[-1, 1]*mesh.g[-1, 1])))


def _interp_segment(x, x0, x1, v0, v1):
    return v0 + (v1 - v0)*(x - x0)/(x1 - x0)


def build_mesh(model: Model, n_elements: int = 100) -> Mesh:
    """Mesh the model with roughly ``n_elements`` elements.

    Every key point of the model (segment ends, loads, supports) becomes a
    node; each interval between key points is divided uniformly so that the
    element size is close to ``beam length / n_elements``.
    """
    model.validate()
    x_start, x_end = model.extent
    h = (x_end - x_start)/max(int(n_elements), 1)

    pts = model.key_points()
    nodes = [pts[0]]
    for a, b in pairwise(pts):
        n = max(1, int(np.ceil((b - a)/h - 1e-9)))
        nodes.extend(np.linspace(a, b, n + 1)[1:])
    nodes = np.array(nodes)

    n = len(nodes) - 1
    xi, xj = nodes[:-1], nodes[1:]
    xm = 0.5*(xi + xj)

    EI = np.empty(n)
    for b in model.beams:
        inside = (xm > b.x0) & (xm < b.x1)
        EI[inside] = b.EI

    k = np.zeros((n, 2))
    g = np.zeros((n, 2))
    for f in model.foundation:
        inside = (xm > f.x0) & (xm < f.x1)
        for col, x in enumerate((xi, xj)):
            k[inside, col] += _interp_segment(x[inside], f.x0, f.x1, *f.k)
            g[inside, col] += _interp_segment(x[inside], f.x0, f.x1, *f.g)

    q = np.zeros((n, 2))
    for d in model.distributed_loads:
        inside = (xm > d.x0) & (xm < d.x1)
        for col, x in enumerate((xi, xj)):
            q[inside, col] += _interp_segment(x[inside], d.x0, d.x1, d.q0, d.q1)

    return Mesh(nodes=nodes, EI=EI, k=k, g=g, q=q)
