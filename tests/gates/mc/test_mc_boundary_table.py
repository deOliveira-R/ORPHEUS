r"""MC's boundary table: ``MC_REALIZERS``, the torus refusal and the #547 strict xfails (wave 2, WF2.3 slice S4; plan sections 2.8, 4.2, 6.4).

Spec: ``.claude/plans/wave2_boundary_carve.md`` section 4.2 (the MC columns), 4.3 delta 9 (decision 31, the user's ruling of the thirteenth
exchange), 6.4. MC is NOT improved on its own terms (``CLAUDE.md``, "The direction of development"): each entry is today's behaviour
re-keyed by type, a harmonisation onto the shared machinery. The MC rows of the shared admission pin carry the grid; THIS file carries the
table's identity, the torus refusal at construction, the defect records and the carry.

What MC computes
================

``MCMesh`` reads a ``Mesh1D`` as a unit cell of side ``pitch`` whose faces are PERIODIC: a particle's position is wrapped by modulo
arithmetic (the registry's one kind, ``periodic``, applied to every face). On a slab that is the translation of the slab's two faces
(``TranslationDeck()`` on each: ``Realized(True)``). On a SOLID CYLINDER it is a periodic equal-area SQUARE lattice cell, a torus
wrap, that no one-dimensional declaration describes: a specular wall declared on the cylinder (the same declaration MoC realizes on its
square cell) is ``Refused(TOPOLOGY, "MC realizes a periodic equal-area square lattice cell (a torus wrap), not the declared cylinder; no
1-D declaration describes it", #547)`` for EVERY weight (the user's ruling, thirteenth exchange; decision 31). The axis is ``TOPOLOGY``
(section 4.3 delta 9, ratified 2026-10-22): the carrier MC samples is a torus whose opposite faces ARE coupled by an image lattice, and
the declared body is not.

The ladder
==========

::

    test_the_mc_table_is_sealed_and_the_old_registry_is_gone    foundation: the table exists; BC_REGISTRY, bc_kind and the factory do not
        ^
    test_mc_refuses_a_specular_wall_on_a_cylinder_for_every_weight   gate 6 (MC half): TOPOLOGY, #547, "torus", at construction
    test_mc_realizes_the_translation_deck_on_a_slab                  the control of the refusal: the same builder admits what MC computes
    test_mc_does_not_compare_the_pitch_with_the_extent              RECORD (#655): green today, reddens the day the defect is fixed
        ^
    test_the_mc_cylinder_sites_are_strict_xfails_naming_547          foundation: the companion of the two xfail marks (#547)
    test_the_carry_of_the_mc_solve                                   bit-equal to the HEAD capture (slab, periodic, fixed seed)

The #547 strict xfails
======================

The 2 MC cylinder sites (``tests/gates/mc/test_gaps.py:692`` and ``tests/gates/mc/test_properties.py:157``; Appendix B rows
``:691`` and ``:156``; the plan's third row, ``:204``, is a SPHERE site and gates the coordinate-system refusal instead) declare the cylinder's outer face as a tag bound to ``WallCondition.specular(1.0)``. They are marked
``@pytest.mark.xfail(strict=True, raises=BoundaryRefused, reason="#547: ...")``: ``raises=BoundaryRefused`` makes any OTHER failure (a
different exception, a wrong number) a hard failure and not an expected one, and strictness turns them red (XPASS) the day MC takes the
two-dimensional geometry value that describes what it computes. The companion row below asserts the marks, by introspection of the two
test functions (a mark that lost ``strict`` or ``raises`` is a silent pass).

First reds at the base tip ``d4da4328``: ``ImportError`` from ``orpheus.mc.solver`` (the retired ``BC``) for every row but the fixture
row; the companion row fails on the two sites (their marks do not exist before the re-spelling lands). Mutations for the reviewers:
M-mc-cylinder-approximated (the specular cell on a cylinder admitted as ``Approximated``: the refusal rows redden and the two xfail
sites XPASS, strictly red), M-mc-weight-gate (the refusal only at weight 1.0: the 0.3 leg reddens), M-bc-kind-restored.
"""

from __future__ import annotations

import ast
from functools import cache
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from orpheus.mc import solver
from tests.gates.mc import test_gaps, test_properties
from tests.gates.transport import _carry_cases as carry
from tests.gates.transport import _wave2_s4_helpers as H

_THIS = "tests/gates/mc/test_mc_boundary_table.py"
_TABLE = f"{_THIS}::test_the_mc_table_is_sealed_and_the_old_registry_is_gone"
_FIXTURE_TEST = "tests/gates/sn/architecture/test_wave2_carry.py::test_the_fixture_is_the_head_capture"
_HERE = Path(__file__).resolve().parent


@pytest.mark.foundation
def test_the_mc_table_is_sealed_and_the_old_registry_is_gone() -> None:
    r"""``MC_REALIZERS`` is a sealed ``RealizerTable`` of method ``"mc"`` on the non-negative-flux cone whose unregistered-unit door cites
    #654; ``MCMesh.BC_REGISTRY``, the instance attribute ``bc_kind`` (0 readers) and the module function ``_mc_bc_periodic`` are gone in the
    same commit as the table.

    Red input at the base tip: ``orpheus.mc.solver`` cannot be imported (the retired ``BC``); once it can, the table does not exist while
    ``BC_REGISTRY`` does. Mutation: M-bc-kind-restored.
    """
    B = H.bc()
    table = H.mc_table()
    assert isinstance(table, B.RealizerTable)
    assert table.method == "mc" and table.state_cone is B.StateCone.NONNEGATIVE_FLUX and table.unregistered_issue == 654
    table.catalogue()
    assert not hasattr(solver.MCMesh, "BC_REGISTRY") and not hasattr(solver, "_mc_bc_periodic")
    from orpheus.geometry.deck import TranslationDeck

    built = H.mc_head(H.mc_slab_mesh1d(left=TranslationDeck(), right=TranslationDeck()), 3.0)
    assert not hasattr(built, "bc_kind"), "bc_kind has 0 readers and retires (section 2.8)"


@pytest.mark.l1
@pytest.mark.rests_on(_TABLE)
@pytest.mark.parametrize("weight", [0.3, 1.0], ids=["w0.3", "w1.0"])
def test_mc_refuses_a_specular_wall_on_a_cylinder_for_every_weight(weight: float) -> None:
    r"""Gate 6 (the MC half): ``MCMesh`` on a solid cylinder whose outer tag is bound to ``WallCondition.specular(w)`` raises
    ``BoundaryRefused(TOPOLOGY, ..., #547)`` at construction, for ``w = 0.3`` AND ``w = 1.0`` (the user's ruling: every weight), the reason
    naming the torus. The faces' name is ``xmax``.

    Red input at the base tip: ``ImportError`` (and, on the HEAD tree, a constructed ``MCMesh`` that ran as periodic whatever law the face
    held). Mutations: M-mc-cylinder-approximated, M-mc-weight-gate.
    """
    B = H.bc()
    mesh = H.chart_mesh("cylinder", outer=0, mat_ids=(2, 0))
    with pytest.raises(B.BoundaryRefused) as info:
        H.mc_head(mesh, 4.0, {0: B.WallCondition.specular(weight)})
    admission = info.value.admission
    assert admission.axis is B.RefusalAxis.TOPOLOGY and admission.issue == 547, f"w = {weight}: {admission!r}"
    assert H.MC_TORUS_FRAGMENT in admission.reason, f"the refusal must name what MC computes: {admission.reason!r}"
    assert info.value.face == "xmax"


@pytest.mark.l1
@pytest.mark.rests_on(f"{_THIS}::test_mc_refuses_a_specular_wall_on_a_cylinder_for_every_weight")
def test_mc_realizes_the_translation_deck_on_a_slab() -> None:
    r"""The control of the refusal (X1: a refusal row whose builder cannot succeed refuses everything): the same slab builder with a
    ``TranslationDeck`` on each face (what MC computes) is ``Realized(True)`` and constructs, and the solve returns a finite eigenvalue.
    A mirror deck on the slab is ``Refused(PLUMBING, #654)`` (the registry held only ``periodic``).

    Red input: ``ImportError`` at the base tip. Mutation: a table that refuses the translation deck.
    """
    import math

    from orpheus.geometry.deck import MirrorDeck, TranslationDeck

    B = H.bc()
    table = H.mc_table()
    assert H.summary(table.admit(TranslationDeck(), H.site("CART_XMAX"))) == ("R", True)
    assert H.summary(table.admit(MirrorDeck(), H.site("CART_XMAX"))) == ("F", "PLUMBING", 654)
    geometry = H.mc_head(H.mc_slab_mesh1d(left=TranslationDeck(), right=TranslationDeck()), 3.0)
    assert math.isfinite(H.mc_solve(geometry).keff)
    with pytest.raises(B.BoundaryRefused) as info:
        H.mc_head(H.mc_slab_mesh1d(left=MirrorDeck(), right=MirrorDeck()), 3.0)
    assert info.value.admission.issue == 654


@pytest.mark.l1
@pytest.mark.rests_on(f"{_THIS}::test_mc_realizes_the_translation_deck_on_a_slab")
def test_mc_does_not_compare_the_pitch_with_the_extent() -> None:
    r"""RECORD of #655 (a defect filed and NOT fixed in this slice): MC's lattice pitch is a separate parameter that nothing compares with
    the declared geometry's extent. A slab of extent 3.0 with a pitch of 7.0 constructs today. The row is green today and reddens the day
    #655 is fixed (the constructor then refuses a pitch that disagrees with the extent), which is the signal to re-pose it as the refusal
    gate; it is a RECORD row (``vv-principles`` test kinds), designed to move at that carve, and carries no ``verifies``.

    Red input: ``ImportError`` at the base tip (the row reads production). Mutation: a pitch/extent comparison added to ``MCMesh``.
    """
    from orpheus.geometry.deck import TranslationDeck

    H.mc_head(H.mc_slab_mesh1d(left=TranslationDeck(), right=TranslationDeck()), 7.0)


#: The two cylinder test functions of decision 31, by name (Appendix B rows ``test_gaps.py:691`` and ``test_properties.py:156``; the lines of
#: ``5d04c296`` plus one are 692 and 157). The plan's third row, ``test_properties.py:204``, declares a SPHERE: its subject is the
#: coordinate-system rejection, which precedes the binding, so it is an ordinary ``pytest.raises`` row and not an xfail.
_XFAIL_SITES = (
    (test_gaps, "test_mcmesh_vs_concentric_keff"),
    (test_properties, "test_mcmesh_cylindrical_lookup"),
)


@pytest.mark.foundation
def test_the_mc_cylinder_sites_are_strict_xfails_naming_547() -> None:
    r"""The companion of decision 31. The two test functions of :data:`_XFAIL_SITES` (the MC cylinder sites that declare a curved face's outer
    wall as the specular tag the cell of MoC uses) carry ONE ``xfail`` each with ``strict=True``, ``raises=BoundaryRefused`` and a reason
    beginning ``#547:``, read through the imported modules (the marks as pytest sees them). Population: exactly 2 functions, each found by
    name, each building an ``MCMesh`` (read by AST: a renamed function or a site that lost its ``MCMesh`` fails the population leg).

    Red input: at the base tip ``d4da4328`` the two test modules cannot be imported (they import the retired ``BC``) and the functions
    carry no mark (the gates commit precedes the re-spelling commit, which adds the marks and turns this row green). Mutation: a mark without ``raises``; without ``strict``; a
    reason without ``#547:``.
    """
    B = H.bc()
    assert len(_XFAIL_SITES) == 2 and len(set(_XFAIL_SITES)) == 2
    for module, func in _XFAIL_SITES:
        module_name = module.__name__.rpartition(".")[2]
        path = _HERE / f"{module_name}.py"
        defs = {n.name: ast.unparse(n) for n in ast.walk(ast.parse(path.read_text())) if isinstance(n, ast.FunctionDef)}
        assert func in defs and "MCMesh(" in defs[func], f"{module_name}::{func}: not found, or no MCMesh construction in its body"
        target = getattr(module, func)
        marks = [m for m in getattr(target, "pytestmark", []) if m.name == "xfail"]
        assert len(marks) == 1, f"{module_name}::{func}: {len(marks)} xfail marks"
        kwargs = marks[0].kwargs
        assert kwargs.get("strict") is True and kwargs.get("raises") is B.BoundaryRefused, f"{module_name}::{func}: {kwargs}"
        assert str(kwargs.get("reason", "")).startswith("#547:"), f"{module_name}::{func}: {kwargs.get('reason')!r}"


@cache
def _fixture() -> Any:
    assert carry.CARRY_FIXTURE.exists(), f"{carry.CARRY_FIXTURE} is missing: the carry fixture is the instrument of this file"
    return np.load(carry.CARRY_FIXTURE)


@pytest.mark.foundation
@pytest.mark.rests_on(_FIXTURE_TEST)
def test_the_fixture_holds_the_mc_row() -> None:
    r"""The fixture holds ``keff_mc/slab/periodic`` captured at ``5d04c296`` and it is physical: the periodic fuel-moderator slab is
    within 0.1 of CP's white-closed slab (``1.2208`` against ``1.2113`` at capture; an MC estimate with 200 neutrons and 20 active cycles
    carries a few percent of statistical noise, so the band is a sanity of the body and not a value claim): a flat or foreign column fails.

    Red input: none today; it reddens on a truncated fixture or one captured from another tree.
    """
    fixture = _fixture()
    assert str(fixture["meta/commit"]) == carry.HEAD_COMMIT
    mc, cp_white = H.fixture_k("mc", "slab", ("periodic",)), H.fixture_k("cp", "slab", ("white", 1.0))
    assert abs(mc - cp_white) < 0.1, f"MC {mc!r} against CP white {cp_white!r}"


@pytest.mark.l1
@pytest.mark.rests_on(f"{_THIS}::test_the_fixture_holds_the_mc_row")
def test_the_carry_of_the_mc_solve() -> None:
    r"""``solve_monte_carlo`` in the new spelling (a ``TranslationDeck`` on each slab face, rule R3) with the fixed seed equals the HEAD
    ``k`` bit for bit: slab ``(0, 1.5, 3)``, materials ``(2, 0)``, 6 + 6 cells, ``pitch = 3.0``, ``seed = 12345``.

    Red input: ``ImportError`` at the base tip. Mutations: a wrap that reads another pitch; a face handled as reflective.
    """
    from orpheus.geometry.deck import TranslationDeck

    geometry = H.mc_head(H.mc_slab_mesh1d(left=TranslationDeck(), right=TranslationDeck()), carry.MC_PARAMS["pitch"])
    k = H.mc_solve(geometry).keff
    expected = H.fixture_k("mc", "slab", ("periodic",))
    assert k == expected, f"k = {k!r}, the HEAD k = {expected!r}"
