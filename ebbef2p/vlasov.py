"""Modified Vlasov foundation model (Vallabhan & Das, 1991).

The soil is an elastic layer of depth ``H`` on a rigid base, with Poisson's
ratio ``nu`` and Young's modulus varying linearly from ``E_top`` at the
surface to ``E_bottom`` at the base. The vertical displacement in the layer
is ``w(x) * phi(z)`` with

    phi(z) = sinh(gamma*(1 - z/H)) / sinh(gamma),

and the foundation parameters follow from the layer properties and the
non-dimensional decay parameter ``gamma``:

    k = (1-nu)/((1+nu)(1-2nu)) * int_0^H E(z) phi'(z)^2 dz
    g = 1/(2(1+nu))            * int_0^H E(z) phi(z)^2  dz     (g = 2t)

``gamma`` itself depends on the beam deflection, so it is found iteratively.
"""

from dataclasses import dataclass, field, replace
from typing import List, Optional

import numpy as np

from . import elements as el
from .solver import solve

# Taylor coefficients in powers of gamma**2 (up to gamma**14) of the four
# gamma-dependent factors in VlasovSoil.parameters: k = C_k (E1 kA + dE kB),
# g = C_g (E1 gA + dE gB). Derived with sympy.
_SERIES = (
    (1, 0, 1/45, -4/945, 1/1575, -8/93555, 1382/127702575, -8/6081075),
    (1/2, -1/12, 1/60, -1/378, 1/2700, -1/20790, 691/116093250, -1/1403325),
    (1/6, -1/45, 1/315, -2/4725, 1/18711, -1382/212837625, 2/2606175,
     -14468/162820783125),
    (1/24, -1/120, 1/756, -1/5400, 1/41580, -691/232186500, 1/2806650,
     -3617/86837751000),
)
_SERIES_GAMMA = 0.3


@dataclass(frozen=True)
class VlasovSoil:
    """Elastic layer for the modified Vlasov model.

    Parameters
    ----------
    E_top : float
        Young's modulus at the surface.
    nu : float
        Poisson's ratio, 0 <= nu < 0.5.
    depth : float
        Layer thickness H.
    E_bottom : float, optional
        Young's modulus at the base; defaults to ``E_top`` (homogeneous layer).
    """
    E_top: float
    nu: float
    depth: float
    E_bottom: Optional[float] = None

    def __post_init__(self):
        if not 0 <= self.nu < 0.5:
            raise ValueError("nu must satisfy 0 <= nu < 0.5")
        if self.depth <= 0 or self.E_top <= 0:
            raise ValueError("depth and E_top must be positive")
        if self.E_bottom is not None and self.E_bottom <= 0:
            raise ValueError("E_bottom must be positive")

    def parameters(self, gamma):
        """Foundation parameters ``(k, g)`` for a given ``gamma`` (g = 2t).

        Below ``_SERIES_GAMMA`` the closed forms lose accuracy to
        cancellation (the error in g grows like 1e-16/gamma**4), so their
        Taylor series in gamma**2 is used instead; both agree to 1e-15 at
        the switch.
        """
        E1 = self.E_top
        E2 = self.E_top if self.E_bottom is None else self.E_bottom
        nu, H = self.nu, self.depth
        if gamma < _SERIES_GAMMA:
            kA, kB, gA, gB = (np.polynomial.polynomial.polyval(gamma**2, c)
                              for c in _SERIES)
        else:
            s2 = np.sinh(gamma)**2
            kA = (2*gamma*np.sinh(2*gamma) + 4*gamma**2)/(8*s2)
            kB = (np.cosh(2*gamma) - 1 + 2*gamma**2)/(8*s2)
            gA = (2*gamma*np.sinh(2*gamma) - 4*gamma**2)/(16*gamma**2*s2)
            gB = (np.cosh(2*gamma) - 1 - 2*gamma**2)/(16*gamma**2*s2)
        k = (1 - nu)/((1 + nu)*(1 - 2*nu)*H)*(E1*kA + (E2 - E1)*kB)
        g = H/(1 + nu)*(E1*gA + (E2 - E1)*gB)
        return float(k), float(g)

    def gamma_from_deflection(self, solution, k, g, soil_beyond_ends):
        """Update gamma from a deflected shape (Vallabhan & Das, 1991).

            (gamma/H)^2 = (1-2nu)/(2(1-nu)) * int w'^2 dx / int w^2 dx

        with the integrals extended beyond the beam ends when the soil
        continues there (w decays as w_end*exp(-sqrt(k/g)*s)).
        """
        mesh = solution.mesh
        gx, gw = np.polynomial.legendre.leggauss(4)
        xi = 0.5*(gx + 1)
        int_w2 = int_dw2 = 0.0
        for e, l in enumerate(mesh.lengths):
            d = solution.u[2*e:2*e + 4]
            w = el.shape_functions(xi, l).T @ d
            dw = el.shape_function_derivatives(xi, l, 1).T @ d
            int_w2 += 0.5*l*np.sum(gw*w**2)
            int_dw2 += 0.5*l*np.sum(gw*dw**2)
        if soil_beyond_ends:
            ends = solution.w[0]**2 + solution.w[-1]**2
            alpha = np.sqrt(k/g)
            int_w2 += ends/(2*alpha)
            int_dw2 += alpha*ends/2
        nu = self.nu
        return float(self.depth*np.sqrt((1 - 2*nu)/(2*(1 - nu))*int_dw2/int_w2))


@dataclass
class VlasovResult:
    """Outcome of :func:`solve_vlasov`."""
    solution: object
    soil: VlasovSoil
    gamma: float
    k: float
    g: float
    history: List[float] = field(default_factory=list)

    @property
    def t(self):
        """Vlasov's t (the second parameter is g = 2t)."""
        return 0.5*self.g


def solve_vlasov(model, soil, n_elements=100, gamma0=1.0, tol=1e-6,
                 max_iter=100):
    """Solve a beam on a modified Vlasov foundation.

    The foundation parameters are uniform along the beam and are found
    iteratively together with ``gamma``. Any foundation already defined in
    ``model`` is not allowed, since the Vlasov layer replaces it.

    Whether the layer continues beyond the beam ends is taken from
    ``model.soil_beyond_ends``; the same assumption is used in the
    stiffness matrix and in the gamma update.

    Parameters
    ----------
    model : Model
        Beam, loads and supports.
    soil : VlasovSoil

    Returns
    -------
    VlasovResult
    """
    if model.foundation:
        raise ValueError("solve_vlasov defines the foundation itself; "
                         "remove the foundation segments from the model")
    x0, x1 = model.extent
    beyond = model.soil_beyond_ends

    def solve_with(gamma):
        k, g = soil.parameters(gamma)
        trial = replace(model, foundation=[])
        trial.add_foundation(x0, x1, k=k, g=g)
        return solve(trial, n_elements), k, g

    gamma = float(gamma0)
    history = [gamma]
    for _ in range(max_iter):
        solution, k, g = solve_with(gamma)
        new_gamma = soil.gamma_from_deflection(solution, k, g, beyond)
        history.append(new_gamma)
        if abs(new_gamma - gamma) < tol*max(1.0, gamma):
            solution, k, g = solve_with(new_gamma)
            return VlasovResult(solution, soil, new_gamma, k, g, history)
        gamma = new_gamma
    raise RuntimeError(f"gamma did not converge in {max_iter} iterations "
                       f"(last values: {history[-3:]})")
