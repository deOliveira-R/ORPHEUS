r"""The carry of wave 2 for CP: the six eigenvalues the carve must not move (law L-K; plan sections 6.0, 6.4 gate 2, 8).

``tests/gates/transport/carry_head.npz`` was captured FROM THE HEAD TREE (``main`` ``5d04c296``, a detached worktree, before any edit) by
``tests/gates/transport/_capture_carry_head.py``: ``solve_cp`` at ``CPParams()`` on the three one-dimensional charts, a slab ``(0, 1.5, 3)``
and a cylinder and a sphere ``(0, 1, 2)``, each of two regions (materials ``1`` then ``0``) with six cells apiece and the two-group library,
under two declarations of the outer wall: vacuum and white (6 rows, key ``keff_cp/<chart>/<declaration>``). The slab's left face is the
mirror CP's kernel images (section 4.2): the HEAD capture declared ``BC.reflective`` there. The same problems, spelled in the NEW
vocabulary (``tests/gates/transport/_wave2_s4_helpers.cp_mesh`` and ``cp_solve``: a ``MirrorDeck`` at the slab's left face, the outer face
a tag bound to ``WallCondition.vacuum()`` or ``white(1.0)``, rule R4), must give the same ``k`` BIT-EQUAL (``==`` on the float; a non-bit
difference is a finding to investigate, never a tolerance to widen).

What each row can and cannot see
================================

* The solve reads every layer between the declaration and ``k``: the condition, the binding, the table, the realized transform
  (``P_cell -> P_inf``: the identity for vacuum, the white closure for the diffuse chain), CP's kernels and the power iteration. A move of
  either transform, of the closure's weight, or of the left-face handling moves a ``k`` here.
* Blind: a defect that is the same for every wall (the kernel), or that needs a heterogeneous body to show: the body is two-region
  two-group, so redistribution is in view.

The fixture's own control (``vv-principles`` anti-pattern 19: a positive reading is no evidence of loading): the row
``test_the_fixture_holds_the_cp_rows`` asserts the ORDER the captured numbers have, ``k_vacuum > k_white`` on 3 of 3 charts. That is the
inverted order of physics, #546 as captured (CP's ``compute_keff`` has no leakage term), and a flat or truncated column fails it. So the
carry freezes the defect together with the numbers, and the carve cannot fix or move it unseen: the day #546 is repaired this row and
the carry rows redden TOGETHER and are re-baselined together (the ``_VACUUM_ABOVE_WHITE`` rows XPASS).

First red at the base tip: ``ImportError`` from ``orpheus.cp.solver`` (the retired ``BC``) for each carry row; the instrument row reads
the file only and is green today. Mutations for the reviewers: M-vacuum-as-white (the empty network realized as the white closure:
reddens the 3 vacuum rows), M-white-as-identity (the diffuse chain realized as the identity: reddens the 3 white rows),
M-left-face-white (the slab's left face realized as the white closure: reddens the slab rows).
"""

from __future__ import annotations

from functools import cache
from typing import Any

import numpy as np
import pytest

from tests.gates.transport import _wave2_s4_helpers as H
from tests.gates.transport import _carry_cases as carry

_THIS = "tests/gates/cp/test_cp_boundary_carry.py"
_FIXTURE_TEST = "tests/gates/sn/architecture/test_wave2_carry.py::test_the_fixture_is_the_head_capture"
_CASES = [(chart, d) for chart in carry.CHARTS for d in carry.CP_DECLARATIONS]
_IDS = [f"{chart}-{carry.tag(d)}" for chart, d in _CASES]


@cache
def _fixture() -> Any:
    assert carry.CARRY_FIXTURE.exists(), f"{carry.CARRY_FIXTURE} is missing: the carry fixture is the instrument of this file"
    return np.load(carry.CARRY_FIXTURE)


def _condition(declaration: tuple) -> Any:
    W = H.bc().WallCondition
    match declaration:
        case ("vacuum",):
            return W.vacuum()
        case ("white", a):
            return W.white(a)
    raise ValueError(declaration)


@pytest.mark.foundation
@pytest.mark.rests_on(_FIXTURE_TEST)
def test_the_fixture_holds_the_cp_rows() -> None:
    r"""The fixture holds the 6 ``keff_cp/*`` rows (3 charts x 2 declarations) captured at ``5d04c296``, and they carry the defect #546
    as captured: on each chart ``k_vacuum > k_white`` (the inverted order; physics says ``<=``). A fixture of a flat column, a truncated
    one, or one captured from another tree fails here.

    Red input: none today (the capture is complete); it reddens on a truncated fixture or one captured from another tree.
    """
    fixture = _fixture()
    assert str(fixture["meta/commit"]) == carry.HEAD_COMMIT
    keys = [carry.other_key("cp", chart, d) for chart, d in _CASES]
    assert len(keys) == 6 and all(k in fixture.files for k in keys), [k for k in keys if k not in fixture.files]
    for chart in carry.CHARTS:
        vacuum, white = (float(fixture[carry.other_key("cp", chart, d)]) for d in carry.CP_DECLARATIONS)
        assert vacuum > white, f"{chart}: k_vacuum {vacuum!r} is not above k_white {white!r}: the capture lost #546"


@pytest.mark.l1
@pytest.mark.rests_on(f"{_THIS}::test_the_fixture_holds_the_cp_rows")
@pytest.mark.parametrize("chart,declaration", _CASES, ids=_IDS)
def test_the_carry_of_the_cp_solves(chart: str, declaration: tuple) -> None:
    r"""``solve_cp`` in the new spelling equals the HEAD ``k`` bit for bit, on the chart and wall of the row.

    Red input: ``ImportError`` at the base tip (every row). Mutations: M-vacuum-as-white, M-white-as-identity, M-left-face-white.
    """
    k = H.cp_solve(H.cp_mesh(chart), {0: _condition(declaration)}).keff
    expected = H.fixture_k("cp", chart, declaration)
    assert k == expected, f"{chart} / {carry.tag(declaration)}: k = {k!r}, the HEAD k = {expected!r}"


@pytest.mark.l1
def test_the_default_pin_cell_outer_face_is_bound_to_white() -> None:
    r"""``solve_cp(materials)`` with no mesh solves the default Wigner-Seitz pin cell (10 + 3 + 7 equal-volume cells) whose outer tag
    ``solve_cp`` binds to ``WallCondition.white()`` itself (at the base tip the white law sat on the geometry). The row pins that
    binding: the default solve equals, bit for bit, the solve on the same mesh with ``{0: WallCondition.white()}`` declared by the
    caller, on a 2-group library of three mixtures, one per material id of the cell.

    Red input: none today (the carry rows above run on bodies the caller builds, and nothing else calls ``solve_cp`` without a mesh).
    Mutation: ``cp-preset-vacuum`` (the preset ``white()`` replaced by ``vacuum()``): the vacuum k is above the white k (#546), so
    the equality fails by the full gap of the inverted order, not by rounding.
    """
    from orpheus.cp.solver import CPParams, solve_cp
    from orpheus.derivations.common.xs_library import get_mixture
    from orpheus.geometry import StructuredGeometry
    from orpheus.mesh import CellsByCount, Mesher
    from orpheus.transport.boundary_conditions import WallCondition

    mesh = Mesher(StructuredGeometry.wigner_seitz_pin_cell()).partition(
        (CellsByCount.uniform_volume(10), CellsByCount.uniform_volume(3), CellsByCount.uniform_volume(7))
    ).mesh
    ids = sorted({int(i) for i in mesh.mat_ids})
    assert len(ids) == 3, ids
    materials = {i: get_mixture(name, "2g") for i, name in zip(ids, "ABC")}
    with H.quiet():
        preset = solve_cp(materials).keff
        declared = solve_cp(materials, mesh, CPParams(), boundary_conditions={0: WallCondition.white()}).keff
        vacuum = solve_cp(materials, mesh, CPParams(), boundary_conditions={0: WallCondition.vacuum()}).keff
    assert preset == declared, f"default solve k = {preset!r}, the white declaration's k = {declared!r}"
    assert vacuum != declared, "the control is blind: vacuum and white gave one k"
