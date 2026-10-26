"""Saving and loading models and Vlasov soils as JSON-ready dictionaries."""

import contextlib
import functools
import io
import json
import pathlib
import re
import runpy

import numpy as np
import pytest

from ebbef2p import Model, ModelError, VlasovSoil, solve, solve_vlasov
from ebbef2p.io import (
    model_from_dict,
    model_to_dict,
    vlasov_soil_from_dict,
    vlasov_soil_to_dict,
)

EXAMPLES = pathlib.Path(__file__).parents[2]/"examples"


def razaqpur_example_2():
    """Razaqpur & Shah (1991), example 2."""
    m = Model().add_beam(0, 15, E=2000, I=1).add_foundation(0, 15, k=64, g=800)
    m.add_point_load(7.5, force=80).add_point_load(2, moment=100)
    m.add_distributed_load(10, 15, 0, 18)
    m.add_support(0, w=0, theta=0)
    for x in (5, 10, 15):
        m.add_support(x)
    return m


def stepped_beam():
    m = Model().add_beam(0, 4, 30e6, 0.01).add_beam(4, 10, 60e6, 0.02)
    m.add_foundation(0, 10, k=(4000, 9000))
    return m.add_point_load(4, force=300, moment=-20)


def jump_in_g():
    m = Model().add_beam(0, 20, 30e6, 0.01)
    m.add_foundation(0, 20, k=5000, g=1000).add_foundation(8, 20, g=(3000, 500))
    return m.add_distributed_load(0, 20, 10)


def soil_beyond_ends():
    m = Model(soil_beyond_ends=True).add_beam(0, 10, 30e6, 0.02)
    return m.add_foundation(0, 10, k=8000, g=4000).add_distributed_load(0, 10, 15)


def prescribed_settlement():
    m = Model().add_beam(0, 6, 2e8, 1e-4).add_distributed_load(0, 6, 5, 0)
    return m.add_support(0, w=0.01, theta=None).add_support(6, w=None, theta=0).add_support(6)


@functools.cache
def example(name):
    with contextlib.redirect_stdout(io.StringIO()):
        return runpy.run_path(str(EXAMPLES/f"{name}.py"), run_name="__main__")


MODELS = {
    "razaqpur_shah": razaqpur_example_2,
    "stepped_beam": stepped_beam,
    "jump_in_g": jump_in_g,
    "soil_beyond_ends": soil_beyond_ends,
    "prescribed_settlement": prescribed_settlement,
    **{f"example_{p.stem}": functools.partial(lambda s: example(s)["model"], p.stem)
       for p in sorted(EXAMPLES.glob("*.py"))},
}


@pytest.fixture(params=list(MODELS))
def model(request):
    return MODELS[request.param]()


def test_round_trip_is_identical(model):
    assert model_from_dict(model_to_dict(model)) == model


def test_json_text_round_trip(model):
    text = json.dumps(model_to_dict(model))
    assert model_from_dict(json.loads(text)) == model


def test_round_trip_gives_same_solution(model):
    loaded = model_from_dict(model_to_dict(model))
    if model.foundation:
        assert np.array_equal(solve(model).u, solve(loaded).u)
    else:
        soil = example("modified_vlasov")["soil"]
        assert solve_vlasov(loaded, soil).gamma == solve_vlasov(model, soil).gamma


def test_output_format():
    assert model_to_dict(razaqpur_example_2()) == {
        "format": "ebbef2p-model", "version": 1,
        "soil_beyond_ends": False,
        "beams": [{"x0": 0.0, "x1": 15.0, "E": 2000.0, "I": 1.0}],
        "foundation": [{"x0": 0.0, "x1": 15.0, "k": [64.0, 64.0], "g": [800.0, 800.0]}],
        "point_loads": [{"x": 7.5, "force": 80.0, "moment": 0.0},
                        {"x": 2.0, "force": 0.0, "moment": 100.0}],
        "distributed_loads": [{"x0": 10.0, "x1": 15.0, "q0": 0.0, "q1": 18.0}],
        "supports": [{"x": 0.0, "w": 0.0, "theta": 0.0}, {"x": 5.0, "w": 0.0, "theta": None},
                     {"x": 10.0, "w": 0.0, "theta": None}, {"x": 15.0, "w": 0.0, "theta": None}],
    }


def test_scalar_k_and_missing_header_accepted():
    m = model_from_dict({"beams": [{"x0": 0, "x1": 15, "E": 2000, "I": 1}],
                         "foundation": [{"x0": 0, "x1": 15, "k": 64}]})
    assert m.foundation[0].k == (64.0, 64.0)
    assert m.foundation[0].g == (0.0, 0.0)
    assert m.soil_beyond_ends is False


def base():
    return model_to_dict(Model().add_beam(0, 5, 1, 1).add_beam(5, 10, 1, 1))


def with_beam(index, **changes):
    data = base()
    data["beams"][index].update(changes)
    return data


@pytest.mark.parametrize("data, message", [
    (with_beam(0, EE=3), "beams[0].EE: unknown key"),
    (with_beam(1, E=-1), "beams[1]: E and I must be positive"),
    (with_beam(0, x1=4), "model: beam segments leave a gap"),
    ({**base(), "version": 2}, "version"),
    ({**base(), "format": "ebbef2p-soil"}, "format"),
    (with_beam(0, E="2000"), "beams[0].E: expected a number"),
    ({**base(), "loads": []}, "loads: unknown key"),
    ({"beams": [{"x0": 0, "x1": 1, "E": 1}]}, "beams[0].I: missing"),
    ({**base(), "soil_beyond_ends": 1}, "soil_beyond_ends: expected true or false"),
])
def test_errors_name_the_field(data, message):
    with pytest.raises(ModelError, match=re.escape(message)):
        model_from_dict(data)


@pytest.mark.parametrize("value", [True, float("nan"), float("inf")])
def test_bool_and_non_finite_are_not_numbers(value):
    with pytest.raises(ModelError, match=re.escape("beams[0].E: expected a number")):
        model_from_dict(with_beam(0, E=value))


def test_partial_model_needs_validate_false():
    for partial in (Model(), Model().add_point_load(3, force=1)):
        data = model_to_dict(partial)
        assert model_from_dict(data, validate=False) == partial
        with pytest.raises(ModelError, match="model: "):
            model_from_dict(data)


@pytest.mark.parametrize("data, message", [
    ([], "expected an object"),
    ({"beams": {}}, "beams: expected a list"),
    ({"beams": [1]}, "beams[0]: expected an object"),
    ({"beams": [{"x0": 0, "x1": 1, "E": 1, "I": 1}],
      "foundation": [{"x0": 0, "x1": 1, "k": [1, 2, 3]}]},
     "foundation[0].k: expected a number or a pair of numbers"),
    ({"supports": [{"x": 0, "w": [0]}]}, "supports[0].w: expected a number or null"),
])
def test_wrong_containers(data, message):
    with pytest.raises(ModelError, match=re.escape(message)):
        model_from_dict(data)


SOILS = [VlasovSoil(144.0, 0.2, 25.0), VlasovSoil(144.0, 0.25, 10.0, E_bottom=300.0)]


@pytest.mark.parametrize("soil", SOILS)
def test_vlasov_soil_round_trip(soil):
    data = vlasov_soil_to_dict(soil)
    assert data["format"] == "ebbef2p-vlasov-soil" and data["version"] == 1
    assert vlasov_soil_from_dict(json.loads(json.dumps(data))) == soil


def test_example_vlasov_soil_gives_same_gamma():
    ns = example("modified_vlasov")
    soil = vlasov_soil_from_dict(vlasov_soil_to_dict(ns["soil"]))
    assert soil == ns["soil"]
    assert solve_vlasov(ns["model"], soil).gamma == solve_vlasov(ns["model"], ns["soil"]).gamma


@pytest.mark.parametrize("data, message", [
    ({"E_top": 144, "nu": 0.5, "depth": 25}, "soil: nu must satisfy"),
    ({"E_top": 144, "nu": 0.2, "depth": 25, "H": 3}, "soil.H: unknown key"),
    ({"E_top": 144, "nu": 0.2}, "soil.depth: missing"),
    ({"E_top": 144, "nu": 0.2, "depth": 25, "format": "ebbef2p-model"}, "format"),
])
def test_vlasov_soil_errors(data, message):
    with pytest.raises(ModelError, match=re.escape(message)):
        vlasov_soil_from_dict(data)
