"""Post-processing of a solved model.

Sign conventions follow :mod:`ebbef2p.model`: ``w`` positive downwards,
``theta = dw/dx``, bending moment ``M = -EI w''`` (positive when sagging),
shear force ``V = dM/dx``. Reactions (soil pressure and support reactions)
are positive when they act upwards on the beam, i.e. against positive loads.
"""

from dataclasses import dataclass, field
from typing import Literal

import numpy as np
import scipy.sparse as sp
from numpy.typing import ArrayLike, NDArray

from . import elements as el
from .mesh import Mesh, end_spring_stiffness
from .model import Model

_GAUSS_X, _GAUSS_W = np.polynomial.legendre.leggauss(4)
_GAUSS_XI = 0.5*(_GAUSS_X + 1)


@dataclass(frozen=True)
class ElementEndValues:
    """A quantity evaluated at both ends of every element.

    ``x`` and ``values`` have length ``2 * n_elements``; consecutive pairs
    belong to one element, so jumps at nodes (under point loads, at
    changes of stiffness) are preserved.
    """
    x: NDArray[np.float64]
    values: NDArray[np.float64]

    def max(self) -> float:
        return float(np.max(self.values))

    def min(self) -> float:
        return float(np.min(self.values))

    def at(self, x: float, side: Literal['left', 'right'] = 'left') -> float:
        """Value at a node, taken from the element on its left or right."""
        idx = np.flatnonzero(np.isclose(self.x, x))
        if idx.size == 0:
            raise ValueError(f"no element end at x = {x}")
        if side == 'left':
            return float(self.values[idx[0]])
        if side == 'right':
            return float(self.values[idx[-1]])
        raise ValueError("side must be 'left' or 'right'")

    def interpolate(self, x: ArrayLike,
                    side: Literal['left', 'right'] = 'left') -> float | NDArray[np.float64]:
        """Value at any ``x`` on the beam, linear between the ends of each element.

        Between nodes this is a linear interpolation of the element-end
        values, so its accuracy depends on the mesh. At a node, or within
        ``np.isclose`` of one as in :meth:`at` (relative tolerance 1e-5,
        absolute 1e-8, so the snapping distance grows with ``|x|``), ``side``
        picks the element on the left or on the right, so jumps are kept; at
        the beam ends it is ignored.

        Returns a float for a scalar ``x`` and an array of the same shape
        otherwise. Raises ``ValueError`` for points outside the beam and for
        NaN or infinite ``x``.
        """
        if side not in ('left', 'right'):
            raise ValueError("side must be 'left' or 'right'")
        x = np.asarray(x, dtype=float)
        if not np.all(np.isfinite(x)):
            raise ValueError("x must be finite")
        starts, ends = self.x[0::2], self.x[1::2]
        nodes = np.append(starts, ends[-1])
        # snap to the nearest node when within isclose of it
        i = np.clip(np.searchsorted(nodes, x), 1, nodes.size - 1)
        nearest = np.where(np.abs(x - nodes[i - 1]) <= np.abs(nodes[i] - x),
                           nodes[i - 1], nodes[i])
        x = np.where(np.isclose(x, nearest), nearest, x)
        outside = (x < starts[0]) | (x > ends[-1])
        if np.any(outside):
            raise ValueError(f"x = {x[outside].ravel()[0]} lies outside the beam "
                             f"[{starts[0]}, {ends[-1]}]")
        if side == 'left':
            e = np.searchsorted(ends, x, side='left')
        else:
            e = np.searchsorted(starts, x, side='right') - 1
        e = np.clip(e, 0, starts.size - 1)
        t = (x - starts[e])/(ends[e] - starts[e])
        value = (1 - t)*self.values[2*e] + t*self.values[2*e + 1]
        return float(value) if value.ndim == 0 else value


@dataclass(frozen=True)
class ReactivePressure(ElementEndValues):
    """Soil reaction under the beam.

    Attributes
    ----------
    x, values
        Distributed pressure p = k w - (g w')' at element ends.
    concentrated : list of (x, force)
        Forces where the shear layer is cut: at the beam ends
        (-g w' at the start, +g w' at the end, plus sqrt(k g) w if the soil
        continues beyond the ends) and at interior nodes where g jumps.
    resultant : float
        Integral of p plus all concentrated forces (total soil reaction).

    :meth:`~ElementEndValues.at` and :meth:`~ElementEndValues.interpolate`
    give the distributed pressure only; the concentrated forces are in
    ``concentrated``.
    """
    concentrated: list[tuple[float, float]] = field(default_factory=list)
    resultant: float = 0.0


@dataclass
class Solution:
    model: Model
    mesh: Mesh
    u: NDArray[np.float64]
    K: sp.csr_matrix
    f: NDArray[np.float64]
    element_stiffness: NDArray[np.float64]
    element_loads: NDArray[np.float64]
    fixed_dofs: NDArray[np.intp]

    # -- nodal values ---------------------------------------------------
    @property
    def x(self) -> NDArray[np.float64]:
        """Node coordinates."""
        return self.mesh.nodes

    @property
    def w(self) -> NDArray[np.float64]:
        """Nodal displacements (positive downwards)."""
        return self.u[0::2]

    @property
    def theta(self) -> NDArray[np.float64]:
        """Nodal rotations dw/dx."""
        return self.u[1::2]

    @property
    def element_displacements(self) -> NDArray[np.float64]:
        """Degrees of freedom of every element, shape (n, 4)."""
        return self.u[self.mesh.element_dofs]

    def element_field(self, xi: ArrayLike, order: int = 0) -> NDArray[np.float64]:
        """``w`` (order 0) or its ``order``-th derivative at local points.

        ``xi`` holds local coordinates in [0, 1], the same in every element;
        the result has shape (n_elements, len(xi)).
        """
        xi = np.atleast_1d(np.asarray(xi, dtype=float))
        l = self.mesh.lengths[:, None]
        N = (el.shape_functions(xi, l) if order == 0
             else el.shape_function_derivatives(xi, l, order))
        return np.einsum('iep,ei->ep', N, self.element_displacements)

    def deflection(self, x: ArrayLike) -> NDArray[np.float64]:
        """Displacement at arbitrary points, from the Hermite interpolation."""
        x = np.atleast_1d(np.asarray(x, dtype=float))
        e = self.mesh.element_at(x)
        l = self.mesh.lengths[e]
        xi = (x - self.mesh.nodes[e])/l
        return np.einsum('im,mi->m', el.shape_functions(xi, l), self.element_displacements[e])

    # -- internal forces --------------------------------------------------
    def _element_end_forces(self) -> NDArray[np.float64]:
        """Forces exerted on each element by its nodes: K_e d_e - f_e."""
        return (np.einsum('eij,ej->ei', self.element_stiffness, self.element_displacements)
                - self.element_loads)

    def _end_coords(self) -> NDArray[np.float64]:
        return np.column_stack([self.mesh.nodes[:-1], self.mesh.nodes[1:]]).ravel()

    def bending_moment(self) -> ElementEndValues:
        """Bending moment at element ends, from element equilibrium."""
        F = self._element_end_forces()
        M = np.column_stack([F[:, 1], -F[:, 3]]).ravel()
        return ElementEndValues(self._end_coords(), M)

    def shear_force(self) -> ElementEndValues:
        """Shear force in the beam at element ends, from element equilibrium.

        The force carried by the foundation's shear layer (g w') is
        excluded, so this is the shear force in the beam itself.
        """
        F = self._element_end_forces()
        th = self.theta
        g = self.mesh.g
        V = np.column_stack([-F[:, 0] - g[:, 0]*th[:-1],
                             F[:, 2] - g[:, 1]*th[1:]]).ravel()
        return ElementEndValues(self._end_coords(), V)

    # -- reactions --------------------------------------------------------
    def reactive_pressure(self) -> ReactivePressure:
        """Soil reaction: distributed pressure and concentrated forces."""
        mesh = self.mesh
        l = mesh.lengths
        k, g = mesh.k, mesh.g
        dg = (g[:, 1] - g[:, 0])/l
        w_end = self.element_field([0.0, 1.0])
        theta_end = self.element_field([0.0, 1.0], order=1)
        curvature_end = self.element_field([0.0, 1.0], order=2)
        values = (k*w_end - g*curvature_end - dg[:, None]*theta_end).ravel()

        flux = g*theta_end                     # g w' at both ends of every element
        k_gauss = k[:, :1] + (k[:, 1:] - k[:, :1])*_GAUSS_XI
        w_gauss = self.element_field(_GAUSS_XI)
        resultant = float((0.5*l) @ ((k_gauss*w_gauss) @ _GAUSS_W)
                          - np.sum(flux[:, 1] - flux[:, 0]))

        x0, x1 = mesh.nodes[0], mesh.nodes[-1]
        f0, f1 = -flux[0, 0], flux[-1, 1]
        if self.model.soil_beyond_ends:
            k_left, k_right = end_spring_stiffness(mesh)
            f0 += k_left*self.w[0]
            f1 += k_right*self.w[-1]

        concentrated = [(float(x0), float(f0))]
        jumps = flux[:-1, 1] - flux[1:, 0]
        scale = max(1e-300, np.max(np.abs(flux[:, 1])))
        for e in np.flatnonzero(np.abs(jumps) > 1e-10*scale):
            concentrated.append((float(mesh.nodes[e + 1]), float(jumps[e])))
        concentrated.append((float(x1), float(f1)))
        resultant += sum(force for _, force in concentrated)

        return ReactivePressure(self._end_coords(), values,
                                concentrated=concentrated,
                                resultant=resultant)

    def support_reactions(self) -> list[tuple[float, float, float]]:
        """Reactions at the supports as a list of (x, force, moment).

        Positive when acting against positive loads (force upwards).
        """
        r = self.f - self.K @ self.u
        out: dict[float, tuple[float, float]] = {}
        for dof in self.fixed_dofs:
            node = dof // 2
            x = float(self.mesh.nodes[node])
            force, moment = out.get(x, (0.0, 0.0))
            if dof % 2 == 0:
                force = float(r[dof])
            else:
                moment = float(r[dof])
            out[x] = (force, moment)
        return [(x, F, M) for x, (F, M) in sorted(out.items())]

    def equilibrium(self) -> dict[str, float]:
        """Vertical force balance: applied loads vs soil and support reactions."""
        applied = sum(p.force for p in self.model.point_loads)
        applied += sum(0.5*(d.q0 + d.q1)*(d.x1 - d.x0)
                       for d in self.model.distributed_loads)
        soil = self.reactive_pressure().resultant
        supports = sum(F for _, F, _ in self.support_reactions())
        return dict(applied=applied, soil=soil, supports=supports,
                    residual=applied - soil - supports)
