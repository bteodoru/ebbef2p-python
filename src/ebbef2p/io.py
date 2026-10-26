"""Saving and loading models as JSON-ready dictionaries.

The format of a model is::

    {
      "format": "ebbef2p-model", "version": 1,
      "soil_beyond_ends": false,
      "beams": [{"x0": 0, "x1": 15, "E": 2000, "I": 1}],
      "foundation": [{"x0": 0, "x1": 15, "k": [64, 64], "g": [800, 800]}],
      "point_loads": [{"x": 7.5, "force": 80, "moment": 0}],
      "distributed_loads": [{"x0": 10, "x1": 15, "q0": 0, "q1": 18}],
      "supports": [{"x": 0, "w": 0, "theta": 0}]
    }

The keys of every item are the fields of the matching dataclass of
:mod:`ebbef2p.model`; missing optional keys take their defaults, unknown
keys are errors. ``k`` and ``g`` may be a number or a pair on input and are
always written as pairs. Models are rebuilt through the ``add_*`` methods,
so they are validated exactly as when built in code; errors are
:class:`~ebbef2p.exceptions.ModelError` with the path of the offending item
or field in the message, e.g. ``beams[1].E: expected a number, got str``.

Use :mod:`json` to turn the dictionaries into text. Inside the package import
this module as ``from . import io as model_io`` so that it is never confused
with the standard library module.
"""

import dataclasses
import math
import typing
from collections.abc import Callable, Mapping
from typing import Any, Literal

from .exceptions import ModelError
from .model import (
    BeamSegment,
    DistributedLoad,
    FoundationSegment,
    Model,
    Pair,
    PointLoad,
    Support,
)

MODEL_FORMAT = "ebbef2p-model"
FORMAT_VERSION = 1

_SECTIONS: dict[str, tuple[type, Callable[..., Model]]] = {
    "beams": (BeamSegment, Model.add_beam),
    "foundation": (FoundationSegment, Model.add_foundation),
    "point_loads": (PointLoad, Model.add_point_load),
    "distributed_loads": (DistributedLoad, Model.add_distributed_load),
    "supports": (Support, Model.add_support),
}

Kind = Literal["number", "number_or_null", "pair", "boolean"]


def _kind(hint: Any) -> Kind:
    """Classify the type hint of a serialized field."""
    if hint is float:
        return "number"
    if hint == float | None:
        return "number_or_null"
    if hint == Pair:
        return "pair"
    if hint is bool:
        return "boolean"
    raise TypeError(f"no JSON form for {hint!r}")


def _serialized_fields(cls: type) -> list[tuple[dataclasses.Field, Kind]]:
    hints = typing.get_type_hints(cls)
    return [(f, _kind(hints[f.name])) for f in dataclasses.fields(cls)]


def _required(f: dataclasses.Field) -> bool:
    return f.default is dataclasses.MISSING and f.default_factory is dataclasses.MISSING


# -- writing --------------------------------------------------------------
def _item_to_dict(item: Any) -> dict[str, Any]:
    return {f.name: list(v) if isinstance(v := getattr(item, f.name), tuple) else v
            for f, _ in _serialized_fields(type(item))}


def model_to_dict(model: Model) -> dict[str, Any]:
    """JSON-ready dictionary of ``model`` (see the module docstring for the format)."""
    data: dict[str, Any] = {"format": MODEL_FORMAT, "version": FORMAT_VERSION,
                            "soil_beyond_ends": model.soil_beyond_ends}
    for key in _SECTIONS:
        data[key] = [_item_to_dict(item) for item in getattr(model, key)]
    return data


# -- reading --------------------------------------------------------------
def _object(data: Any, path: str) -> Mapping[str, Any]:
    if not isinstance(data, Mapping):
        raise ModelError(f"{path}: expected an object, got {type(data).__name__}")
    return data


def _is_number(value: Any) -> bool:
    return (isinstance(value, (int, float)) and not isinstance(value, bool)
            and math.isfinite(value))


def _value(kind: Kind, value: Any, path: str) -> Any:
    if kind == "number" and _is_number(value):
        return value
    if kind == "number_or_null" and (value is None or _is_number(value)):
        return value
    if kind == "pair":
        if _is_number(value):
            return value
        if isinstance(value, (list, tuple)) and len(value) == 2 and all(map(_is_number, value)):
            return tuple(value)
    if kind == "boolean" and isinstance(value, bool):
        return value
    expected = {"number": "a number", "number_or_null": "a number or null",
                "pair": "a number or a pair of numbers", "boolean": "true or false"}[kind]
    raise ModelError(f"{path}: expected {expected}, got {value!r}")


def _check_keys(data: Mapping[str, Any], allowed: typing.Iterable[str], prefix: str) -> None:
    allowed = set(allowed)
    for key in data:
        if key not in allowed:
            raise ModelError(f"{prefix}{key}: unknown key")


def _fields_from_dict(cls: type, data: Any, path: str) -> dict[str, Any]:
    """Checked keyword arguments for ``cls`` (or its builder) from ``data``."""
    data = _object(data, path)
    spec = _serialized_fields(cls)
    _check_keys(data, (f.name for f, _ in spec), f"{path}.")
    kwargs = {}
    for f, kind in spec:
        if f.name in data:
            kwargs[f.name] = _value(kind, data[f.name], f"{path}.{f.name}")
        elif _required(f):
            raise ModelError(f"{path}.{f.name}: missing")
    return kwargs


def _check_header(data: Mapping[str, Any], expected_format: str) -> None:
    fmt = data.get("format", expected_format)
    if fmt != expected_format:
        raise ModelError(f"format: expected {expected_format!r}, got {fmt!r}")
    version = data.get("version", FORMAT_VERSION)
    if type(version) is not int or version != FORMAT_VERSION:
        raise ModelError(f"version: unsupported version {version!r} "
                         f"(this ebbef2p reads version {FORMAT_VERSION})")


def model_from_dict(data: Mapping[str, Any], validate: bool = True) -> Model:
    """Build a :class:`~ebbef2p.model.Model` from a dictionary.

    Parameters
    ----------
    data : mapping
        As written by :func:`model_to_dict`, or by hand. Without ``format``
        and ``version`` it is read as version 1.
    validate : bool
        If True, also run :meth:`Model.validate` (segments tile the beam,
        everything lies on it). False allows an incomplete model to be
        reloaded, e.g. one that has no beam yet.

    Raises
    ------
    ModelError
        With the path of the item or field in the message.
    """
    data = _object(data, "model")
    _check_keys(data, ("format", "version", "soil_beyond_ends", *_SECTIONS), "")
    _check_header(data, MODEL_FORMAT)
    model = Model(soil_beyond_ends=_value("boolean", data.get("soil_beyond_ends", False),
                                          "soil_beyond_ends"))
    for key, (cls, add) in _SECTIONS.items():
        items = data.get(key, [])
        if not isinstance(items, list):
            raise ModelError(f"{key}: expected a list, got {type(items).__name__}")
        for i, item in enumerate(items):
            path = f"{key}[{i}]"
            kwargs = _fields_from_dict(cls, item, path)
            try:
                add(model, **kwargs)
            except ModelError as err:
                raise ModelError(f"{path}: {err}") from err
    if validate:
        try:
            model.validate()
        except ModelError as err:
            raise ModelError(f"model: {err}") from err
    return model
