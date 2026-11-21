r"""The retirement of wave 2 (WF2.4, unit R): the typed-law vocabulary and its three homes are gone, and nothing names them.

The carve (slices S1 to S5) replaced the law classes of ``orpheus/geometry/boundary/`` (16 modules), the realization seam
``orpheus/transport/boundary_realization.py`` and the SN shim ``orpheus/sn/boundary/bound_compat.py``. A retirement that only
deletes loses coverage and leaves orphans (``retirement-audit`` C and G), so this file was the done-when of the deletion commit:
every check it held was RED at the base tip (``8c5578c5`` plus the S5 checkpoint) and green exactly when the retirement is complete.
Plan: ``.claude/plans/wave2_boundary_carve.md`` section 7 (7.1 census, 7.2 dependency table, 7.3 the edits).

The name and module census of this retirement is the shared engine's (``tests/gates/test_retirements.py``, rows in
``tests/gates/_retired/rows_wave2.py``); this file keeps the capture-script, cold-import and method-space checks.

The claims kept here, each with its kind (``vv-principles``) and the input that reddens it (X1):

* the layer gate's cold-import list names no retired module (a retired module imported in a fresh interpreter is an
  ``ImportError``, so the gate would turn red by the deletion alone).
* ERR-041's index sets retire: ``SNMethodSpace`` carries no ``inflow_indices`` / ``outflow_indices`` (the typed operators read
  their half-spaces off the trace ladder, and ERR-041's invariant is carried there by
  ``test_angular_trace_space.py::TestTheInflowHalfSpaceIsTheIncomingOrdinates``).

**No exemption, and the two capture scripts.** ``tests/gates/transport/_capture_carry_head.py`` and
``tests/gates/derivations/_capture_reference_carry_head.py`` spelled the retired API by design (frozen records of how
``carry_head.npz`` and ``reference_carry_head.json`` were captured from ``5d04c296``). The ruling of 2026-10-24 deletes them: the
fixtures stay, each reading gate's docstring names the capture commit and the script's git path (recoverable with ``git show``),
and until the deletion the scripts were live users the census read like any other file.
"""

from __future__ import annotations

import ast
import dataclasses
from pathlib import Path

import pytest

from tests.gates._retired import RetiredModule, rows_wave2

_REPO = Path(__file__).resolve().parents[2]


_CAPTURE_SCRIPTS: tuple[str, ...] = (
    "tests/gates/transport/_capture_carry_head.py",
    "tests/gates/derivations/_capture_reference_carry_head.py",
)
_CARRY_FIXTURES: tuple[str, ...] = ("tests/gates/transport/carry_head.npz", "tests/gates/derivations/reference_carry_head.json")


@pytest.mark.foundation
def test_the_capture_scripts_are_deleted_and_their_fixtures_stay() -> None:
    r"""The two frozen capture scripts no longer exist (ruling of 2026-10-24: no exemption remains), and the fixtures they were run
    to produce stay (the carry gates read them). Red input at the base tip: both scripts are present. Green when both are deleted
    and both fixtures are files."""
    present = [rel for rel in _CAPTURE_SCRIPTS if (_REPO / rel).exists()]
    assert not present, f"the frozen capture scripts still exist: {present}"
    for rel in _CARRY_FIXTURES:
        assert (_REPO / rel).is_file(), rel


@pytest.mark.foundation
def test_the_cold_import_list_of_the_layer_gate_names_no_retired_module() -> None:
    r"""``test_entry_point_imports_in_a_fresh_interpreter`` imports each listed entry in a cold interpreter; a retired module in
    the list is an ``ImportError`` the day it is deleted. The list is read as the string constants of that test's decorators (AST).

    Positive control: the same list holds the live entry ``orpheus.mesh`` and at least 30 entries. Red input at the base tip: the
    three entries ``orpheus.geometry.boundary``, ``orpheus.transport.boundary_realization`` and ``orpheus.sn.boundary.bound_compat``.
    """
    tree = ast.parse((_REPO / "tests" / "gates" / "test_layer_imports.py").read_text(encoding="utf-8"))
    fn = next(n for n in ast.walk(tree) if isinstance(n, ast.FunctionDef) and n.name == "test_entry_point_imports_in_a_fresh_interpreter")
    entries = [c.value for d in fn.decorator_list for c in ast.walk(d) if isinstance(c, ast.Constant) and isinstance(c.value, str)]
    assert len(entries) >= 30 and "orpheus.mesh" in entries, entries
    retired_modules = tuple(r.module for r in rows_wave2.ROWS if isinstance(r, RetiredModule))
    retired = [e for e in entries if e in retired_modules]
    assert not retired, f"the cold-import list names retired modules: {retired}"


@pytest.mark.foundation
def test_the_sn_method_space_carries_no_index_sets() -> None:
    r"""ERR-041's index sets are retired: ``SNMethodSpace`` has no ``inflow_indices`` and no ``outflow_indices`` field. The typed
    operators read their half-spaces off the trace ladder (``trace.inflow_space(face)``), so the pair had 0 production readers
    outside the zero map's own cross-check, and the invariant that guard checked is carried by the ladder's own witness
    (``tests/gates/numerics/test_angular_trace_space.py::TestTheInflowHalfSpaceIsTheIncomingOrdinates``).

    The method space also names no ``face`` and holds no ``mesh``: the face is the ``FaceSite`` beside the space, so a realization
    has one face value (the space's own, once cross-checked against the site's, was two encodings of one orientation).

    Read from the dataclass fields (a name grep would read a homonym local). Positive control: the same class still carries
    ``quadrature`` and ``trace``. Red input at the base tip: the index fields, ``face`` and ``mesh`` exist.
    """
    from orpheus.sn.mesh.method_space import SNMethodSpace

    names = {f.name for f in dataclasses.fields(SNMethodSpace)}
    assert {"quadrature", "trace"} <= names, sorted(names)
    assert not {"inflow_indices", "outflow_indices", "face", "mesh"} & names, sorted(names)
