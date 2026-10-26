"""Quick-look plots."""

import pytest

from ebbef2p import Model, solve

pytest.importorskip("matplotlib").use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


def test_plot_solution_draws_four_panels():
    from ebbef2p.plotting import plot_solution

    m = Model().add_beam(0, 10, 30e6, 0.01).add_foundation(0, 10, k=5000, g=2000)
    m.add_point_load(5, force=100)
    try:
        axes = plot_solution(solve(m, 50))
        assert [ax.get_ylabel() for ax in axes] == ["w", "M", "V", "p"]
        assert all(ax.lines for ax in axes)
    finally:
        plt.close("all")
