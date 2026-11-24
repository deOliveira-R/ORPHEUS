r"""G01b: ``numerics/metric.py`` holds no numerical cutoff (N01, spec 7.2, 3.3).

``ast`` over the module: no call named ``pinv`` or ``allclose`` and no keyword named ``rcond``, ``rtol`` or ``atol``, with no
exemption. The symmetry admission ``_DENSE_METRIC_SYMMETRY_RTOL`` is a module constant, neither a call nor a keyword, so the
predicate does not see it. The predicate is also run on the module source with a forbidden construct appended (the positive
control of the instrument): each of three forms must read red. First red today: ``pinv`` at two sites and ``rcond`` at two.
"""

from __future__ import annotations

import ast
import importlib.util
from pathlib import Path

import pytest

pytestmark = pytest.mark.foundation

FORBIDDEN_CALLS = ("pinv", "allclose")
FORBIDDEN_KEYWORDS = ("rcond", "rtol", "atol")


def _require(cond: object, msg: str) -> None:
    if not cond:
        pytest.fail(msg)


def metric_source_path() -> Path:
    spec = importlib.util.find_spec("orpheus.numerics.metric")
    _require(spec is not None and spec.origin is not None, "orpheus.numerics.metric is not found")
    assert spec is not None and spec.origin is not None
    return Path(spec.origin)


def cutoff_hits(source: str) -> list[tuple[str, int]]:
    """``(name, line)`` of each call named a forbidden call and each keyword named a forbidden keyword, in a parse of ``source``."""
    hits: list[tuple[str, int]] = []
    for node in ast.walk(ast.parse(source)):
        if isinstance(node, ast.Call):
            f = node.func
            name = f.attr if isinstance(f, ast.Attribute) else f.id if isinstance(f, ast.Name) else None
            if name in FORBIDDEN_CALLS:
                hits.append((str(name), node.lineno))
            for kw in node.keywords:
                if kw.arg in FORBIDDEN_KEYWORDS:
                    hits.append((str(kw.arg), kw.value.lineno))
    return sorted(hits)


def test_metric_module_holds_no_cutoff() -> None:
    hits = cutoff_hits(metric_source_path().read_text())
    _require(not hits, f"numerics/metric.py holds {len(hits)} cutoff constructs (name, line): {hits}")


@pytest.mark.parametrize(
    "appended,expected",
    [
        ("np.linalg.pinv(None, rcond=1e-12)", {"pinv", "rcond"}),
        ("np.allclose(1.0, 1.0)", {"allclose"}),
        ("f(1.0, rtol=1e-3)", {"rtol"}),
        ("f(1.0, atol=1e-3)", {"atol"}),
    ],
)
def test_the_predicate_reads_red_on_a_forbidden_construct(appended: str, expected: set[str]) -> None:
    """Positive control: the same predicate on the real source with one forbidden construct appended reads exactly its names."""
    source = metric_source_path().read_text()
    before = cutoff_hits(source)
    after = cutoff_hits(source + "\n" + appended + "\n")
    new = {name for name, _ in after} - {name for name, _ in before} if before else {name for name, _ in after}
    _require(len(after) == len(before) + len(expected), f"appending {appended!r} added {len(after) - len(before)} hits, expected {len(expected)}")
    _require(expected <= {name for name, _ in after}, f"appending {appended!r} did not read {expected}: {sorted(new)}")
