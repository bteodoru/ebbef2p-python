"""Post-processing of a solved model.

Sign conventions follow :mod:`ebbef2p.model`: ``w`` positive downwards,
``theta = dw/dx``, bending moment ``M = -EI w''`` (positive when sagging),
shear force ``V = dM/dx``. Reactions (soil pressure and support reactions)
are positive when they act upwards on the beam, i.e. against positive loads.
"""

from dataclasses import dataclass, field
from typing import List, Tuple

import numpy as np

from . import elements as el
from .mesh import end_spring_stiffness


@dataclass(frozen=True)
class ElementEndValues:
    """A quantity evaluated at both ends of every element.

    ``x`` and ``values`` have length ``2 * n_elements``; consecutive pairs
    belong to one element, so jumps at nodes (under point loads, at
    changes of stiffness) are preserved.
    """
    x: np.ndarray
    values: np.ndarray

    def max(self):
        return float(np.max(self.values))

    def min(self):
        return float(np.min(self.values))

    def at(self, x, side='left'):
        """Value at a node, taken from the element on its left or right."""
        idx = np.flatnonzero(np.isclose(self.x, x))
        if idx.size == 0:
            raise ValueError(f"no element end at x = {x}")
        if side == 'left':
            return float(self.values[idx[0]])
        if side == 'right':
            return float(self.values[idx[-1]])
        raise ValueError("side must be 'left' or 'right'")


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
    """
    concentrated: List[Tuple[float, float]] = field(default_factory=list)
    resultant: float = 0.0


@dataclass
class Solution:
    model: object
    mesh: object
    u: np.ndarray
    K: object
    f: np.ndarray
    element_stiffness: np.ndarray
    element_loads: np.ndarray
    fixed_dofs: np.ndarray

    # -- nodal values ---------------------------------------------------
    @property
    def x(self):
        """Node coordinates."""
        return self.mesh.nodes

    @property
    def w(self):
        """Nodal displacements (positive downwards)."""
        return self.u[0::2]

    @property
    def theta(self):
        """Nodal rotations dw/dx."""
        return self.u[1::2]

    def _element_dofs(self, e):
        return self.u[2*e:2*e + 4]

    def deflection(self, x):
        """Displacement at arbitrary points, from the Hermite interpolation."""
        x = np.atleast_1d(np.asarray(x, dtype=float))
        out = np.empty_like(x)
        for n, xv in enumerate(x):
            e = self.mesh.element_at(xv)
            l = self.mesh.lengths[e]
            xi = (xv - self.mesh.nodes[e])/l
            out[n] = el.shape_functions(xi, l) @ self._element_dofs(e)
        return out

    # -- internal forces --------------------------------------------------
    def _element_end_forces(self):
        """Forces exerted on each element by its nodes: K_e d_e - f_e."""
        d = np.stack([self._element_dofs(e) for e in range(self.mesh.n_elements)])
        return np.einsum('eij,ej->ei', self.element_stiffness, d) - self.element_loads

    def _end_coords(self):
        return np.column_stack([self.mesh.nodes[:-1], self.mesh.nodes[1:]]).ravel()

    def bending_moment(self):
        """Bending moment at element ends, from element equilibrium."""
        F = self._element_end_forces()
        M = np.column_stack([F[:, 1], -F[:, 3]]).ravel()
        return ElementEndValues(self._end_coords(), M)

    def shear_force(self):
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
    def reactive_pressure(self):
        """Soil reaction: distributed pressure and concentrated forces."""
        mesh = self.mesh
        gx, gw = np.polynomial.legendre.leggauss(4)
        values = np.empty(2*mesh.n_elements)
        flux_start = np.empty(mesh.n_elements)
        flux_end = np.empty(mesh.n_elements)
        resultant = 0.0

        for e, l in enumerate(mesh.lengths):
            d = self._element_dofs(e)
            ki, kj = mesh.k[e]
            gi, gj = mesh.g[e]
            dg = (gj - gi)/l
            wi, thi, wj, thj = d
            values[2*e] = (ki*wi - gi*(el.shape_function_derivatives(0.0, l, 2) @ d)
                           - dg*thi)
            values[2*e + 1] = (kj*wj - gj*(el.shape_function_derivatives(1.0, l, 2) @ d)
                               - dg*thj)
            flux_start[e] = gi*thi
            flux_end[e] = gj*thj

            xi = 0.5*(gx + 1)
            kq = ki + (kj - ki)*xi
            wq = el.shape_functions(xi, l).T @ d
            resultant += 0.5*l*np.sum(gw*kq*wq)
            resultant -= flux_end[e] - flux_start[e]

        x0, x1 = mesh.nodes[0], mesh.nodes[-1]
        f0, f1 = -flux_start[0], flux_end[-1]
        if self.model.soil_beyond_ends:
            k_left, k_right = end_spring_stiffness(mesh)
            f0 += k_left*self.w[0]
            f1 += k_right*self.w[-1]

        concentrated = [(float(x0), float(f0))]
        jumps = flux_end[:-1] - flux_start[1:]
        scale = max(1e-300, np.max(np.abs(flux_end)))
        for e in np.flatnonzero(np.abs(jumps) > 1e-10*scale):
            concentrated.append((float(mesh.nodes[e + 1]), float(jumps[e])))
        concentrated.append((float(x1), float(f1)))
        resultant += sum(force for _, force in concentrated)

        return ReactivePressure(self._end_coords(), values,
                                concentrated=concentrated,
                                resultant=float(resultant))

    def support_reactions(self):
        """Reactions at the supports as a list of (x, force, moment).

        Positive when acting against positive loads (force upwards).
        """
        r = self.f - self.K @ self.u
        out = {}
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

    def equilibrium(self):
        """Vertical force balance: applied loads vs soil and support reactions."""
        applied = sum(p.force for p in self.model.point_loads)
        applied += sum(0.5*(d.q0 + d.q1)*(d.x1 - d.x0)
                       for d in self.model.distributed_loads)
        soil = self.reactive_pressure().resultant
        supports = sum(F for _, F, _ in self.support_reactions())
        return dict(applied=applied, soil=soil, supports=supports,
                    residual=applied - soil - supports)
