"""Strip footing crossing a weak zone.

The first foundation parameter drops over a soft lens between x = 6 and
x = 10, and the shear parameter changes abruptly at x = 8. The soil is
assumed to continue beyond both ends of the footing.
"""

from ebbef2p import Model, solve

L = 20.0
model = Model(soil_beyond_ends=True)
model.add_beam(0, L, E=30e6, I=0.02)

model.add_foundation(0, 6, k=20000)
model.add_foundation(6, 8, k=(20000, 6000))
model.add_foundation(8, 10, k=(6000, 20000))
model.add_foundation(10, L, k=20000)
model.add_foundation(0, 8, g=4000)
model.add_foundation(8, L, g=8000)

model.add_distributed_load(0, L, 150)
model.add_point_load(5, force=400)
model.add_point_load(15, force=400)

solution = solve(model, n_elements=200)

pressure = solution.reactive_pressure()
print(f"max settlement   : {solution.w.max():.5f} at x = {solution.x[solution.w.argmax()]:.2f}")
print(f"max moment       : {solution.bending_moment().max():.1f}")
print(f"min moment       : {solution.bending_moment().min():.1f}")
print("concentrated soil forces:")
for x, force in pressure.concentrated:
    print(f"  x = {x:5.2f}   F = {force:8.2f}")
eq = solution.equilibrium()
print(f"applied {eq['applied']:.2f} = soil {eq['soil']:.2f} + supports {eq['supports']:.2f}")
