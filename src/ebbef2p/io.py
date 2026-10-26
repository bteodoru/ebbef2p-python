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
keys are errors. A key is required only if the matching ``add_*`` argument
has no default, so ``q1`` may be omitted for a uniform load. ``k`` and ``g``
may be a number or a pair on input and are always written as pairs. Models
are rebuilt through the ``add_*`` methods, so they are validated exactly as
when built in code; errors are :class:`~ebbef2p.exceptions.ModelError`
whose message starts with the path of the offending item or field, which is
also in their ``path`` attribute, e.g. ``beams[1].E: E must be positive``.

A :class:`~ebbef2p.vlasov.VlasovSoil` is written in the same way::

    {"format": "ebbef2p-vlasov-soil", "version": 1,
     "E_top": 144, "nu": 0.2, "depth": 25, "E_bottom": null}

:func:`model_json_schema` and :func:`vlasov_soil_json_schema` give JSON Schemas
(draft 2020-12) of both formats, generated from the dataclasses: every field
is described by its ``doc`` and its dimension (see :mod:`ebbef2p.dimensions`).
They check the structure only; values are checked when the model is built.

Use :mod:`json` to turn the dictionaries into text. Inside the package import
this module as ``from . import io as model_io`` so that it is never confused
with the standard library module.
"""

import copy
import dataclasses
import inspect
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
from .vlasov import VlasovSoil

MODEL_FORMAT = "ebbef2p-model"
VLASOV_SOIL_FORMAT = "ebbef2p-vlasov-soil"
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


def _required_args(builder: Callable[..., Any]) -> set[str]:
    """Arguments of an ``add_*`` method (or class) that have no default."""
    return {name for name, p in inspect.signature(builder).parameters.items()
            if name != "self" and p.default is inspect.Parameter.empty}


def _error(path: str, message: str) -> ModelError:
    return ModelError(f"{path}: {message}", path=path)


def _located(err: ModelError, path: str) -> ModelError:
    """``err`` from a builder, located at ``path`` (and its field, if known)."""
    where = f"{path}.{err.field}" if err.field else path
    return _error(where, str(err))


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
        raise _error(path, f"expected an object, got {type(data).__name__}")
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
    raise _error(path, f"expected {expected}, got {value!r}")


def _check_keys(data: Mapping[str, Any], allowed: typing.Iterable[str], prefix: str) -> None:
    allowed = set(allowed)
    for key in data:
        if key not in allowed:
            raise _error(f"{prefix}{key}", "unknown key")


def _fields_from_dict(cls: type, builder: Callable[..., Any], data: Any,
                      path: str) -> dict[str, Any]:
    """Checked keyword arguments for ``builder`` (which builds ``cls``) from ``data``."""
    data = _object(data, path)
    spec = _serialized_fields(cls)
    _check_keys(data, (f.name for f, _ in spec), f"{path}.")
    required = _required_args(builder)
    kwargs = {}
    for f, kind in spec:
        if f.name in data:
            kwargs[f.name] = _value(kind, data[f.name], f"{path}.{f.name}")
        elif f.name in required:
            raise _error(f"{path}.{f.name}", "missing")
    return kwargs


def _check_header(data: Mapping[str, Any], expected_format: str) -> None:
    fmt = data.get("format", expected_format)
    if fmt != expected_format:
        raise _error("format", f"expected {expected_format!r}, got {fmt!r}")
    version = data.get("version", FORMAT_VERSION)
    # any JSON number equal to 1 (also 1.0), as in the schema; not a bool or a string
    if isinstance(version, bool) or not isinstance(version, (int, float)) \
            or version != FORMAT_VERSION:
        raise _error("version", f"unsupported version {version!r} "
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
        With the path of the item or field at the start of the message and
        in its ``path`` attribute.
    """
    data = _object(data, "model")
    _check_keys(data, ("format", "version", "soil_beyond_ends", *_SECTIONS), "")
    _check_header(data, MODEL_FORMAT)
    model = Model(soil_beyond_ends=_value("boolean", data.get("soil_beyond_ends", False),
                                          "soil_beyond_ends"))
    for key, (cls, add) in _SECTIONS.items():
        items = data.get(key, [])
        if not isinstance(items, list):
            raise _error(key, f"expected a list, got {type(items).__name__}")
        for i, item in enumerate(items):
            path = f"{key}[{i}]"
            kwargs = _fields_from_dict(cls, add, item, path)
            try:
                add(model, **kwargs)
            except ModelError as err:
                raise _located(err, path) from err
    if validate:
        try:
            model.validate()
        except ModelError as err:
            raise _error("model", str(err)) from err
    return model


def vlasov_soil_to_dict(soil: VlasovSoil) -> dict[str, Any]:
    """JSON-ready dictionary of a :class:`~ebbef2p.vlasov.VlasovSoil`."""
    return {"format": VLASOV_SOIL_FORMAT, "version": FORMAT_VERSION, **_item_to_dict(soil)}


def vlasov_soil_from_dict(data: Mapping[str, Any]) -> VlasovSoil:
    """Build a :class:`~ebbef2p.vlasov.VlasovSoil` from a dictionary.

    Raises
    ------
    ModelError
        With the path of the field in the message (``soil.depth: missing``).
    """
    data = _object(data, "soil")
    header = {key: data[key] for key in ("format", "version") if key in data}
    _check_header(header, VLASOV_SOIL_FORMAT)
    fields = {key: value for key, value in data.items() if key not in header}
    kwargs = _fields_from_dict(VlasovSoil, VlasovSoil, fields, "soil")
    try:
        return VlasovSoil(**kwargs)
    except ModelError as err:
        raise _located(err, "soil") from err


# -- JSON Schema ------------------------------------------------------------
_SCHEMA_DIALECT = "https://json-schema.org/draft/2020-12/schema"
_NUMBER = {"type": "number"}
_JSON_TYPES: dict[Kind, dict[str, Any]] = {
    "number": _NUMBER,
    "number_or_null": {"type": ["number", "null"]},
    "pair": {"oneOf": [_NUMBER, {"type": "array", "items": _NUMBER,
                                 "minItems": 2, "maxItems": 2}]},
    "boolean": {"type": "boolean"},
}


def _description(f: dataclasses.Field) -> str:
    dimension = f.metadata.get("dimension")
    return f"{f.metadata['doc']} [{dimension}]" if dimension else f.metadata["doc"]


def _object_schema(cls: type, builder: Callable[..., Any]) -> dict[str, Any]:
    properties = {}
    for f, kind in _serialized_fields(cls):
        prop = {**copy.deepcopy(_JSON_TYPES[kind]), "description": _description(f)}
        if f.default is not dataclasses.MISSING:
            prop["default"] = list(f.default) if isinstance(f.default, tuple) else f.default
        properties[f.name] = prop
    return {"type": "object", "description": cls.__doc__.splitlines()[0],
            "properties": properties,
            "required": [f.name for f, _ in _serialized_fields(cls)
                         if f.name in _required_args(builder)],
            "additionalProperties": False}


def _header_schema(fmt: str) -> dict[str, Any]:
    return {"format": {"const": fmt, "description": "format name (optional on input)"},
            "version": {"const": FORMAT_VERSION,
                        "description": "format version (optional on input)"}}


def model_json_schema() -> dict[str, Any]:
    """JSON Schema (draft 2020-12) of the dictionaries read by :func:`model_from_dict`."""
    beyond = next(f for f in dataclasses.fields(Model) if f.name == "soil_beyond_ends")
    properties = {**_header_schema(MODEL_FORMAT),
                  "soil_beyond_ends": {"type": "boolean", "default": beyond.default,
                                       "description": _description(beyond)}}
    for key, (cls, add) in _SECTIONS.items():
        properties[key] = {"type": "array", "items": _object_schema(cls, add), "default": []}
    return {"$schema": _SCHEMA_DIALECT, "title": "EBBEF2p model", "type": "object",
            "description": "Beam on a two-parameter elastic foundation. No units are built "
                           "in: use any consistent set (e.g. kN and m). Dimensions are given "
                           "in brackets as powers of force F and length L. Numbers must be "
                           "finite.",
            "properties": properties, "additionalProperties": False}


def vlasov_soil_json_schema() -> dict[str, Any]:
    """JSON Schema (draft 2020-12) of the dictionaries read by :func:`vlasov_soil_from_dict`."""
    schema = _object_schema(VlasovSoil, VlasovSoil)
    schema["properties"] = {**_header_schema(VLASOV_SOIL_FORMAT), **schema["properties"]}
    schema["description"] += " Numbers must be finite."
    return {"$schema": _SCHEMA_DIALECT, "title": "EBBEF2p modified Vlasov soil", **schema}
