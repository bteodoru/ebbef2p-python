"""Free beam on a modified Vlasov foundation (Vallabhan & Das, 1991).

The foundation parameters are not given; they follow from the elastic
layer (E, nu, depth) and the iterated decay parameter gamma.
"""

from ebbef2p import Model, VlasovSoil, solve_vlasov

L, E, I = 100.0, 432000.0, 2.25

model = Model()
model.add_beam(0, L, E=E, I=I)
model.add_distributed_load(0, L, 3.0)

soil = VlasovSoil(E_top=E/3000, nu=0.2, depth=L/4)
result = solve_vlasov(model, soil, n_elements=100)

print(f"gamma = {result.gamma:.4f}  (Vallabhan & Das: 0.395)")
print(f"k = {result.k:.4f}   2t = {result.g:.4f}")
print("gamma iterations:", ", ".join(f"{g:.4f}" for g in result.history))

s = result.solution
print(f"settlement at the ends  : {s.w[0]:.4f}")
print(f"settlement at mid-span  : {s.w[len(s.w)//2]:.4f}")
