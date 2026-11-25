r"""G04i: every space constructor states its variance, and the bare ``FunctionSpace(name, shape)`` call is gone (N04, spec 7.5).

``TensorSpace.declared(name, shape)`` (no variance), ``of_measure(measure, label)`` and ``of_basis(basis)`` raise ``TypeError``; the same with a variance builds.
``ast`` over ``orpheus/`` and ``tests/``: 0 calls of the bare ``FunctionSpace(`` (a call whose callee is the name ``FunctionSpace`` or an attribute of that name).

First red at the base tip: ``TensorSpace`` does not exist (import-only for the behavioural legs); the ``ast`` leg reads the 4 production and 291 test
bare calls of the opening (re-measured here by the same parse). Mutation: give ``variance`` a default in ``declared`` reddens the first leg. Exact on any host.
"""

from __future__ import annotations

import ast
from pathlib import Path

import numpy as np
import pytest

pytestmark = pytest.mark.foundation

ROOT = Path(__file__).resolve().parents[4]


def _require(cond: object, msg: str) -> None:
    if not cond:
        pytest.fail(msg)


def test_a_constructor_without_a_variance_is_refused() -> None:
    from orpheus.numerics import factor as F
    from orpheus.numerics.space import TensorSpace

    with pytest.raises(TypeError):
        TensorSpace.declared("x", (3,))  # type: ignore[call-arg]
    with pytest.raises(TypeError):
        TensorSpace.of_measure(None, "x")  # type: ignore[call-arg]
    with pytest.raises(TypeError):
        TensorSpace.of_basis(None)  # type: ignore[call-arg]
    built = TensorSpace.declared("x", (3,), F.Variance(F.Kind.FUNCTION, F.VarianceChart.DENSITY))
    _require(built.shape == (3,), f"declared(name, shape, variance) reads shape {built.shape}")
    _require(np.array_equal(built.apply_metric(np.array([1.0, 2.0, 3.0])), np.array([1.0, 2.0, 3.0])), "a declared space with no form is not counting")


def _bare_calls(source: str) -> int:
    n = 0
    for node in ast.walk(ast.parse(source)):
        if isinstance(node, ast.Call):
            f = node.func
            if (isinstance(f, ast.Name) and f.id == "FunctionSpace") or (isinstance(f, ast.Attribute) and f.attr == "FunctionSpace"):
                n += 1
    return n


def test_no_bare_function_space_call_remains() -> None:
    _require(_bare_calls("x = FunctionSpace('a', (3,))\ny = m.FunctionSpace('b', (2,))") == 2, "the census does not find a planted bare call of either spelling")
    per_file: dict[str, int] = {}
    population = 0
    for tree in ("orpheus", "tests"):
        for path in sorted((ROOT / tree).rglob("*.py")):
            population += 1
            n = _bare_calls(path.read_text())
            if n:
                per_file[str(path.relative_to(ROOT))] = n
    _require(population > 1000, f"the census read {population} files")
    _require(not per_file, f"{sum(per_file.values())} bare FunctionSpace( calls in {len(per_file)} of {population} files: {sorted(per_file.items())[:8]}")
