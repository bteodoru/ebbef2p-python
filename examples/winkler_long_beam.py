"""Long beam on a Winkler foundation under a central point load.

Compares the finite element result with Hetenyi's closed-form solution.
"""

from ebbef2p import Model, solve

L, E, I = 2.0, 1.0, 1.0
k, P = 10000.0, 100.0

model = Model()
model.add_beam(0, L, E=E, I=I)
model.add_foundation(0, L, k=k)
model.add_point_load(L/2, force=P)

solution = solve(model, n_elements=100)

lam = (k/(4*E*I))**0.25
print(f"w max     FE: {solution.w.max():.6f}   Hetenyi: {P*lam/(2*k):.6f}")
print(f"M max     FE: {solution.bending_moment().max():.4f}   Hetenyi: {P/(4*lam):.4f}")
print(f"residual of vertical equilibrium: {solution.equilibrium()['residual']:.2e}")
