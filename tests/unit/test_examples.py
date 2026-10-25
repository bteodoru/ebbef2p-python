"""The example scripts used in the documentation must run."""

import pathlib
import runpy

import pytest

EXAMPLES = sorted((pathlib.Path(__file__).parents[2]/"examples").glob("*.py"))


@pytest.mark.parametrize("script", EXAMPLES, ids=lambda p: p.name)
def test_example_runs(script, capsys):
    runpy.run_path(str(script), run_name="__main__")
    assert capsys.readouterr().out
