"""Quick-look plots of a :class:`~ebbef2p.results.Solution` (needs matplotlib)."""


def plot_solution(solution, axes=None):
    """Deflection, bending moment, shear force and soil pressure, stacked.

    Returns the array of matplotlib axes.
    """
    import matplotlib.pyplot as plt

    if axes is None:
        _, axes = plt.subplots(4, 1, sharex=True, figsize=(7, 9))

    import numpy as np
    xs = np.linspace(solution.x[0], solution.x[-1], 400)
    panels = [
        (xs, solution.deflection(xs), "w", True),
        (*_pairs(solution.bending_moment()), "M", True),
        (*_pairs(solution.shear_force()), "V", False),
        (*_pairs(solution.reactive_pressure()), "p", True),
    ]
    for ax, (x, y, label, invert) in zip(axes, panels):
        ax.plot(x, y)
        ax.axhline(0, color="0.6", lw=0.8)
        ax.set_ylabel(label)
        if invert:
            ax.invert_yaxis()
    axes[-1].set_xlabel("x")
    return axes


def _pairs(values):
    return values.x, values.values
