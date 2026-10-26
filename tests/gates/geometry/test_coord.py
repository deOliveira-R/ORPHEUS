r"""The one measure of a 1-D coordinate system (P1 step 3a, round 2; re-filed into geometry by the layer DAG, ruling 23).

The measure of the cells of an interval is owned by the coordinate system,
:math:`m_j = c\,(T(r_{j+1}) - T(r_j))` with :math:`T(r) = r^d` evaluated by
numpy (``CoordSystem.measure``, :mod:`orpheus.geometry.coord`). This file pins
that definition (``TestTheOneMeasure``): the measure is ``c * diff(T(edges))``
bitwise, ``compute_volumes_1d`` is a delegate of it, the measure coordinate
``T`` is the system's, and a ``StructuredGeometry`` asks it. It also holds the
witness that the tests' oracle of the measure is written once
(``TestTheMeasureOracleHasOneHome``).

The class lived in ``tests/gates/mesh/test_partition.py`` until the layer-DAG
re-filing (``.claude/plans/layer_dag.md``, 2026-10-26, candidate 1 of the
Phase 3 sub-cone report): its claims are the coordinate system's and the
geometry's, so geometry is the lowest layer that states them. Its ROUTE row,
that every MESH reader (the interval rules' equal shares and realised shells)
asks the one definition, names mesh objects and stays there
(``test_partition.py::TestTheOneMeasure::test_every_reader_asks_the_one_definition``).

Gate ids S3.1 to S3.7 are those of the P1 verification specification
(``.claude/plans/reference_p1_spec.md`` §1.3a). Every test is ``foundation``
(no theory-page label, so no ``verifies``); claim kind THEOREM unless marked
RECORD.
"""
from __future__ import annotations

import ast
import itertools
import re
from pathlib import Path

import numpy as np
import pytest

import orpheus.geometry.coord as coord_module
from orpheus.geometry import CoordSystem, compute_volumes_1d
from orpheus.geometry.coord import MeasureCoordinate
from tests.gates.geometry._measure import CONSTANT as _CONSTANT
from tests.gates.geometry._measure import EXPONENT as _EXPONENT
from tests.gates.geometry._measure import LAW_INTERVALS as _LAW_INTERVALS
from tests.gates.geometry._measure import measure as _measure
from tests.gates.geometry._measure import r2_edges as _r2_edges
from tests.gates.geometry._measure import report as _report
from tests.gates.geometry._measure import uniform_geometry as _geometry

pytestmark = pytest.mark.foundation

_ROOT = Path(__file__).resolve().parents[3]

_SLAB = CoordSystem.CARTESIAN
_CYLINDER = CoordSystem.CYLINDRICAL
_SPHERE = CoordSystem.SPHERICAL
_COORDS = (_SLAB, _CYLINDER, _SPHERE)


class TestTheOneMeasure:
    r"""The measure has one definition, owned by the coordinate system (the
    ruling of 2026-09-29): ``coord.measure(edges) = c * diff(T(edges))``,
    ``T`` numpy's power; the geometry asks it; ``compute_volumes_1d`` is it.

    Claim kind: THEOREM. The route row, which swaps the definition for a
    decoy and requires every MESH reader to move, is
    ``tests/gates/mesh/test_partition.py::TestTheOneMeasure::test_every_reader_asks_the_one_definition``.
    """

    @pytest.mark.parametrize("coord", _COORDS, ids=lambda c: c.name.lower())
    def test_the_measure_is_c_times_the_difference_of_T(self, coord):
        failures = []
        for (a, b), n in itertools.product(_LAW_INTERVALS, (1, 5, 64)):
            edges = _r2_edges(1, a, b, n)
            if not np.array_equal(coord.measure(edges), _measure(coord, edges)):
                failures.append(f"[{a}, {b}] n={n}")
        _report(failures, 3 * len(_LAW_INTERVALS), f"{coord.name} measure")

    def test_compute_volumes_1d_is_a_delegate(self):
        """RECORD of the retirement's residue: the legacy name survives only
        as a one-line delegate (3b's migration retires it)."""
        import inspect

        body = inspect.getsource(coord_module.compute_volumes_1d)
        assert "coord.measure(edges)" in body
        assert "match" not in body

    @pytest.mark.parametrize("coord", _COORDS, ids=lambda c: c.name.lower())
    def test_the_measure_coordinate_is_the_systems(self, coord):
        T = coord.measure_coordinate
        assert T == MeasureCoordinate(_EXPONENT[coord])
        r = np.array([0.0, 0.3, 1.7, 2.0])
        np.testing.assert_array_equal(T(r), r ** _EXPONENT[coord])
        np.testing.assert_allclose(T.inverse(T(r)), r, rtol=4 * np.finfo(float).eps, atol=0.0)
        assert coord.measure_constant == _CONSTANT[coord]

    def test_the_identity_inverse_does_not_alias(self):
        t = np.array([0.0, 1.0])
        assert MeasureCoordinate(1).inverse(t) is not t

    @pytest.mark.parametrize("exponent", [0, 4, -1, 1.5], ids=["0", "4", "-1", "1.5"])
    def test_a_measure_coordinate_is_r_r2_or_r3(self, exponent):
        with pytest.raises(ValueError, match=re.escape("a measure coordinate is r, r**2 or r**3")):
            MeasureCoordinate(exponent)  # type: ignore[arg-type]  # a refusal input

    def test_the_geometry_s_intervals_and_measure(self):
        g = _geometry(_CYLINDER, (0.0, 0.5, 1.5, 2.0))
        assert g.intervals == ((0.0, 0.5), (0.5, 1.5), (1.5, 2.0))
        edges = np.array([0.0, 0.5, 1.5, 2.0])
        np.testing.assert_array_equal(g.measure(edges), _CYLINDER.measure(edges))


# ─────────────────────────────────────────────────────────────────────
# The tests' oracle of the measure is written once (layer DAG ruling 31)
# ─────────────────────────────────────────────────────────────────────

#: The files that read the measure's oracle: the definition pin (this file) and the three mesh consumers.
_CONSUMERS = (
    "tests/gates/geometry/test_coord.py",
    "tests/gates/mesh/test_partition.py",
    "tests/gates/mesh/test_mesher.py",
    "tests/gates/mesh/test_mesh1d.py",
)
_ORACLE = "tests/gates/geometry/_measure.py"
_ORACLE_MODULE = "tests.gates.geometry._measure"
_COORD_MEMBERS = {"CARTESIAN": _SLAB, "CYLINDRICAL": _CYLINDER, "SPHERICAL": _SPHERE}
_LOCAL_ALIASES = {"_SLAB": _SLAB, "_CYLINDER": _CYLINDER, "_SPHERE": _SPHERE, "_CART": _SLAB, "_CYL": _CYLINDER, "_SPH": _SPHERE}


def _key_member(key: ast.expr | None) -> CoordSystem | None:
    if isinstance(key, ast.Attribute) and key.attr in _COORD_MEMBERS:
        return _COORD_MEMBERS[key.attr]
    if isinstance(key, ast.Name) and key.id in _LOCAL_ALIASES:
        return _LOCAL_ALIASES[key.id]
    return None


_OPERATORS = {ast.Add: float.__add__, ast.Sub: float.__sub__, ast.Mult: float.__mul__, ast.Div: float.__truediv__}


def _value(node: ast.expr) -> float | None:
    """The number a value expression spells (literals, ``np.pi`` or ``math.pi``, ``+ - * /``); ``None`` for anything else
    (a lambda, a call, a string): such a value is not a written copy of the measure's numbers."""
    if isinstance(node, ast.Constant) and isinstance(node.value, int | float) and not isinstance(node.value, bool):
        return float(node.value)
    if isinstance(node, ast.Attribute) and node.attr == "pi" and isinstance(node.value, ast.Name) and node.value.id in {"np", "math"}:
        return float(np.pi)
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub):
        inner = _value(node.operand)
        return None if inner is None else -inner
    if isinstance(node, ast.BinOp) and type(node.op) in _OPERATORS:
        left, right = _value(node.left), _value(node.right)
        return None if left is None or right is None else _OPERATORS[type(node.op)](left, right)
    return None


def _measure_copies(source: str) -> list[str]:
    """Each dict literal keyed by the three coordinate systems whose values are the measure's exponents or constants."""
    found = []
    for node in ast.walk(ast.parse(source)):
        if not isinstance(node, ast.Dict) or len(node.keys) != 3:
            continue
        members = [_key_member(k) for k in node.keys]
        if None in members or len(set(members)) != 3:
            continue
        spelled = {m: _value(v) for m, v in zip(members, node.values, strict=True)}
        values = {m: v for m, v in spelled.items() if v is not None}
        if len(values) != 3:
            continue
        for written in (_EXPONENT, _CONSTANT):
            if all(np.isclose(values[m], float(written[m]), rtol=1e-12, atol=0.0) for m in written):  # any spelling of the number
                found.append(f"line {node.lineno}: {ast.unparse(node)}")
    return found


#: The production attributes and names the oracle must never read.
_PRODUCTION_MEASURE = {"measure", "measure_constant", "measure_coordinate", "compute_volumes_1d", "compute_volumes"}


def _production_reads(source: str) -> list[str]:
    reads = []
    for node in ast.walk(ast.parse(source)):
        if isinstance(node, ast.Attribute) and node.attr in _PRODUCTION_MEASURE:
            reads.append(f"line {node.lineno}: .{node.attr}")
        elif isinstance(node, ast.Name) and node.id in _PRODUCTION_MEASURE:
            reads.append(f"line {node.lineno}: {node.id}")
        elif isinstance(node, ast.alias) and node.name.split(".")[-1] in _PRODUCTION_MEASURE:
            reads.append(f"line {getattr(node, 'lineno', '?')}: import {node.name}")
    return reads


class TestTheMeasureOracleHasOneHome:
    r"""The X4 witness of ruling 31. The mechanism that keeps the tests'
    measure single is one module, ``tests/gates/geometry/_measure.py``,
    imported by every consumer; these rows are what reddens when a consumer
    writes its own copy again (the state of the tree until 2026-10-26: two
    full copies and a third exponent table) or when the oracle starts reading
    the production definition it pins (which would demote every gate that
    compares with it, ``retirement-audit`` D.14).

    Population: the four consumer files of ``_CONSUMERS``, not the tree.
    Claim kind: THEOREM over that population, each scan with its positive
    control.
    """

    def test_the_consumers_define_no_copy(self) -> None:
        planted = (
            "_P = {CoordSystem.CARTESIAN: 1, CoordSystem.CYLINDRICAL: 2, CoordSystem.SPHERICAL: 3}\n"
            "_Q = {_SLAB: 1.0, _CYLINDER: np.pi, _SPHERE: 4 * math.pi / 3}\n"
            "_R = {_SLAB: 'slab', _CYLINDER: 'cylinder', _SPHERE: 'sphere'}\n"
        )
        planted_found = _measure_copies(planted)
        assert len(planted_found) == 2, planted_found  # the positive control: both shapes found, the labels not
        assert len(_measure_copies((_ROOT / _ORACLE).read_text())) == 2  # the one allowed definition is seen
        copies, missing_import = [], []
        for rel in _CONSUMERS:
            source = (_ROOT / rel).read_text()
            copies += [f"{rel} {hit}" for hit in _measure_copies(source)]
            imported = {
                node.module for node in ast.walk(ast.parse(source))
                if isinstance(node, ast.ImportFrom) and node.module
            }
            if _ORACLE_MODULE not in imported:
                missing_import.append(rel)
        assert len(_CONSUMERS) == 4
        assert not copies, "a consumer writes the measure again:\n  " + "\n  ".join(copies)
        assert not missing_import, f"consumers not reading {_ORACLE_MODULE}: {missing_import}"

    def test_the_oracle_reads_no_production_measure(self) -> None:
        assert _production_reads("def f(coord, e):\n    return coord.measure(e)\n") == ["line 2: .measure"]  # the positive control
        assert _production_reads("from orpheus.geometry import compute_volumes_1d\n")  # an import is a read
        reads = _production_reads((_ROOT / _ORACLE).read_text())
        assert not reads, f"{_ORACLE} reads the production measure it pins: {reads}"
