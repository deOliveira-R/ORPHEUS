r"""CP's boundary table: ``CP_REALIZERS`` and what it refuses (wave 2, WF2.3 slice S4; plan sections 2.8, 4.2, 6.4 gates 4 and 5).

Spec: ``.claude/plans/wave2_boundary_carve.md`` section 4.2 (the CP columns), 4.3 deltas 3 and 5, section 6.4. CP is NOT improved on its
own terms (``CLAUDE.md``, "The direction of development"): each entry is today's behaviour re-keyed by type and by face, a
harmonisation onto the shared machinery. The CP rows of the shared admission pin (``tests/gates/transport/test_admission_matrix.py``)
carry the grid cell by cell; THIS file carries the rest of CP's contract.

The ladder
==========

::

    test_the_cp_table_is_sealed_and_the_old_registry_is_gone     foundation: the table exists; BC_REGISTRY and its factories do not
        ^
    test_cp_decks_and_the_hollow_inner_face_are_pinned           reference: the deck and the hollow-inner cells the shared grid cannot hold
        ^
    test_cp_left_face_is_a_forced_mirror                         gate 5: the xmin face of a CART site admits MirrorDeck() only
    test_cp_albedo_below_one_is_refused                          gate 4: white(0.5) is Refused(PLUMBING, #654), a silent drop today
    test_cp_white_at_unit_weight_and_vacuum_are_admitted         the control of gate 4: the same builder passes for what CP realizes

The behaviours the rows pin, each read at the line of today's code that makes it so
===================================================================================

* **The left face (gate 5).** CP's slab kernel images a mirror at the left face whatever is declared there (``[M]`` 2026-09-25 to
  2026-09-30: a white, vacuum or undeclared left law gave ``k = 1.8749980808246423`` to all 16 digits; CP declared white|white
  1.212883 against SN reflective|white 1.212884 and SN white|white 1.212537), so the slab CP computes is half of a symmetric slab. Today
  ``_refuse_a_law_cp_drops`` admits only a reflective left law and raises ``NotImplementedError`` otherwise (#513). After the carve
  the refusal is admission: the ``xmin`` face of a ``CART`` site admits ``MirrorDeck()`` only, a tag bound to ANY condition, and a
  ``TranslationDeck``, are ``Refused(PLUMBING, #513)`` (the migrated ``tests/gates/mesh/test_dropped_laws_are_refused.py::TestCP``,
  whose ``NotImplementedError`` expectations become ``BoundaryRefused``; the ``wrong-axis-mirror`` row has no successor, because a
  ``MirrorDeck`` has no axis, rule R5). The right face of the slab and the outer face of a radial body: ``white(1.0)`` and the empty
  network are ``Realized(True)``, everything else is ``Refused(PLUMBING, #654)``.
* **White below one (gate 4).** ``WhiteBoundary(albedo=0.5)`` is realized TODAY with the albedo dropped (the registry key was
  ``white`` whatever alpha): ``k = 0.6 = k_inf`` on 3 of 3 geometries on the homogeneous one-group probe (``_probe_admission_num.out``),
  against SN's 0.536. After the carve ``WallCondition.white(0.5)`` is ``Refused(PLUMBING, #654)`` (section 4.3 delta 3, announced at the
  S4 checkpoint); the diffuse chain at weight 1.0 is ``Realized(True)``.
* **The hollow inner face.** Section 4.2 states no CP ``RAD_INNER`` column; hand-off (b) of S3 says CP does not call
  ``hollow_inner_refusal`` (it is the SN and diffusion refusal, #511). TODAY a hollow CP mesh with a declared inner law is refused by
  ``_refuse_a_law_cp_drops`` with ``NotImplementedError`` citing #513 ("the declared inner law would be dropped"), for every law
  (``tests/gates/mesh/test_dropped_laws_are_refused.py::TestCP::test_a_hollow_body_with_an_inner_law_is_refused``). Today's behaviour
  re-keyed by face is ``Refused(PLUMBING, #513)`` for every unit at CP's ``RAD_INNER`` site, and the hollow-inner row below pins exactly
  that. It is the builder's to confirm or amend (the orchestrator's ruling: NEEDS 1 of the S4 gates report); the row is separate from the
  admission pin so that a ruling edits one test and not the grid.

First reds at the base tip ``d4da4328``: ``ImportError`` for every row (``orpheus.cp.solver`` imports the retired ``BC`` from
``orpheus.geometry``), each at the accessor of its own first production read. Mutations for the reviewers: M-xmin-admits-white (the
left-face rule dropped: the five ``test_cp_left_face_is_a_forced_mirror`` condition rows redden, the deck rows stay); M-albedo-dropped
(``white(w)`` admitted at every ``w``, the old registry's coverage: the three gate-4 rows redden and the pin rows ``white 0.5``);
M-inner-honoured (CP's inner face admitting a condition: the hollow-inner rows redden).
"""

from __future__ import annotations

import importlib

import pytest

from tests.gates.transport import _wave2_s4_helpers as H

_THIS = "tests/gates/cp/test_cp_realizer_table.py"
_TABLE = f"{_THIS}::test_the_cp_table_is_sealed_and_the_old_registry_is_gone"
_CHARTS = ("slab", "cylinder", "sphere")

#: The names the carve retires from ``orpheus.cp.solver`` in the SAME commit as their replacements (section 2.8).
_RETIRED = ("_refuse_a_law_cp_drops", "_cp_bc_white", "_cp_bc_vacuum")


@pytest.mark.foundation
def test_the_cp_table_is_sealed_and_the_old_registry_is_gone() -> None:
    r"""``CP_REALIZERS`` is a sealed ``RealizerTable`` of method ``"cp"`` on the non-negative-flux cone whose unregistered-unit door cites
    #654 (section 4.2: CP, MoC and MC keep the old string registries' coverage); ``CPMesh.BC_REGISTRY``, ``CPMesh._resolve_bc`` and the
    three module functions of ``_RETIRED`` are gone in the same commit.

    Red input at the base tip: ``orpheus.cp.solver`` cannot be imported (the retired ``BC``), and once it can, the table does not exist
    while ``BC_REGISTRY`` does. Mutation: the table built with ``unregistered_issue=659`` (SN's); ``BC_REGISTRY`` restored.
    """
    solver = importlib.import_module("orpheus.cp.solver")
    B = H.bc()
    table = H.cp_table()
    assert isinstance(table, B.RealizerTable)
    assert table.method == "cp" and table.state_cone is B.StateCone.NONNEGATIVE_FLUX and table.unregistered_issue == 654
    table.catalogue()  # raises TableNotSealedError iff the table is not sealed
    assert not hasattr(solver.CPMesh, "BC_REGISTRY") and not hasattr(solver.CPMesh, "_resolve_bc"), "the string registry must be retired"
    stale = [name for name in _RETIRED if hasattr(solver, name)]
    assert not stale, f"retired in the same commit as the table, still defined: {stale}"


@pytest.mark.foundation
@pytest.mark.rests_on(_TABLE)
def test_cp_decks_and_the_hollow_inner_face_are_pinned() -> None:
    r"""The cells of CP's column the shared value grid cannot hold, typed by hand from section 4.2: the two decks at the slab's two faces
    (``MirrorDeck``: ``Realized(True)`` at ``xmin`` only, ``Refused(PLUMBING, #654)`` elsewhere; ``TranslationDeck``: ``Refused(PLUMBING,
    #513)`` at ``xmin``, #654 elsewhere), and every unit at the hollow inner face (``Refused(PLUMBING, #513)``, today's refusal: see the
    module docstring; the section states no cell).

    Population: 2 decks x 3 Cartesian/radial-outer sites + 6 value classes x 1 inner site = 12 cells.
    Red input: ``ImportError`` at the base tip. Mutation: M-inner-honoured, M-xmin-admits-white.
    """
    B, table = H.bc(), H.cp_table()
    from orpheus.geometry.deck import MirrorDeck, TranslationDeck

    deck_cells = {
        ("mirror deck", "CART_XMIN"): ("R", True), ("mirror deck", "CART_XMAX"): ("F", "PLUMBING", 654),
        ("translation deck", "CART_XMIN"): ("F", "PLUMBING", 513), ("translation deck", "CART_XMAX"): ("F", "PLUMBING", 654),
    }
    decks = {"mirror deck": MirrorDeck(), "translation deck": TranslationDeck()}
    for (name, where), expected in deck_cells.items():
        assert H.summary(table.admit(decks[name], H.site(where))) == expected, f"{name} @ {where}"
    W = B.WallCondition
    inner = {
        "vacuum": W.vacuum(), "white 1.0": W.white(1.0), "specular 1.0": W.specular(1.0), "specular 0.5": W.specular(0.5),
        "both 0.3 0.5": W.from_amplitudes(0.3, 0.5), "inflow 1.0": W.uniform(1.0),
    }
    for label, declared in inner.items():
        assert H.summary(table.admit(declared, H.site("RAD_INNER"))) == ("F", "PLUMBING", 513), f"{label} @ RAD_INNER"


def _slab_with_left(left: object, conditions: dict):
    return H.cp_head(H.cp_mesh("slab", left=left, outer=0), conditions)


@pytest.mark.l1
@pytest.mark.rests_on(_TABLE)
@pytest.mark.parametrize(
    "label", ["vacuum", "white 1.0", "white 0.5", "specular 0.7", "specular 1.0"],
)
def test_cp_left_face_is_a_forced_mirror(label: str) -> None:
    r"""Gate 5 of 6.4, the migrated #513 guard. The slab's left face carries a TAG bound to ``label``'s condition (the right face a
    second tag bound to ``white(1.0)``, the one wall CP realizes there): ``CPMesh`` refuses at the face ``xmin``,
    ``BoundaryRefused`` whose admission is ``Refused(PLUMBING, #513)``. The five conditions are the migrated rows of
    ``TestCP::test_a_slab_whose_left_law_is_not_the_mirror_is_refused`` (vacuum, white, partial specular; the wrong-axis mirror dissolves,
    rule R5), plus the unit specular wall (the permutation CP would replace by a full mirror) and ``white(0.5)``.

    Red input at the base tip: ``ImportError``. After the carve without the left-face rule (a table keyed on the site class only): the
    ``specular 1.0`` and ``white 1.0`` rows build, so they redden (M-xmin-admits-white).
    """
    W = H.bc().WallCondition
    declared = {"vacuum": W.vacuum(), "white 1.0": W.white(1.0), "white 0.5": W.white(0.5), "specular 0.7": W.specular(0.7), "specular 1.0": W.specular(1.0)}
    with pytest.raises(H.bc().BoundaryRefused) as info:
        _slab_with_left(1, {1: declared[label], 0: W.white(1.0)})
    assert info.value.face == "xmin", f"{label}: refused at {info.value.face!r}, the left face is xmin"
    admission = info.value.admission
    assert admission.axis is H.bc().RefusalAxis.PLUMBING and admission.issue == 513, f"{label}: {admission!r}"


@pytest.mark.l1
@pytest.mark.rests_on(f"{_THIS}::test_cp_left_face_is_a_forced_mirror")
def test_cp_left_face_admits_the_mirror_deck_and_refuses_the_translation() -> None:
    r"""The positive leg of gate 5 (``vv-principles`` anti-pattern 11): a ``MirrorDeck`` on the left face builds, with the right face
    ``white(1.0)`` or vacuum (the migrated ``test_a_slab_with_the_mirror_on_its_left_builds``); a ``TranslationDeck`` pair is refused at
    ``xmin`` with #513 (a slab with a translation deck needs one on each face, so the right face is the partner deck).

    Red input: ``ImportError`` at the base tip. Mutation: M-xmin-admits-white does not touch this row; a table that refuses the mirror
    at ``xmin`` (the deck realizer returning the plumbing refusal) reddens the first leg.
    """
    B = H.bc()
    from orpheus.geometry.deck import TranslationDeck

    for right in (B.WallCondition.white(1.0), B.WallCondition.vacuum()):
        _slab_with_left(None, {0: right})  # builds: the default left of cp_mesh is MirrorDeck()
    with pytest.raises(B.BoundaryRefused) as info:
        H.cp_head(H.cp_mesh("slab", left=TranslationDeck(), outer=TranslationDeck()))
    assert info.value.face == "xmin" and info.value.admission.issue == 513 and info.value.admission.axis is B.RefusalAxis.PLUMBING


@pytest.mark.l1
@pytest.mark.rests_on(_TABLE)
@pytest.mark.parametrize("chart", _CHARTS)
def test_cp_albedo_below_one_is_refused(chart: str) -> None:
    r"""Gate 4 of 6.4. ``WallCondition.white(0.5)`` on CP raises ``BoundaryRefused(PLUMBING, #654)`` at the outer face, on 3 of 3
    charts, at the head AND at the solve entry. Today ``WhiteBoundary(albedo=0.5)`` is realized with ``alpha`` dropped
    (the registry key was ``white`` whatever alpha): ``k = 0.6 = k_inf`` on 3 of 3 geometries (``_probe_admission_num.out``, a silent
    wrong answer), the red input of the carve's behaviour change (section 4.3 delta 3). ``specular(0.5)`` and ``specular(1.0)`` are
    refused the same way (the principled axis is angular resolution, ``[R]``; the issue is #654).

    Red input at the base tip: ``ImportError`` (and, on the HEAD tree, a realized ``k``, not a refusal). Mutation: M-albedo-dropped.
    """
    B, W = H.bc(), H.bc().WallCondition
    for label, declared in (("white 0.5", W.white(0.5)), ("specular 0.5", W.specular(0.5)), ("specular 1.0", W.specular(1.0))):
        with pytest.raises(B.BoundaryRefused) as at_head:
            H.cp_head(H.cp_mesh(chart), {0: declared})
        assert at_head.value.admission.axis is B.RefusalAxis.PLUMBING and at_head.value.admission.issue == 654, f"{chart} / {label}: {at_head.value.admission!r}"
        assert at_head.value.face == "xmax"
        with pytest.raises(B.BoundaryRefused):
            H.cp_solve(H.cp_mesh(chart), {0: declared})


@pytest.mark.l1
@pytest.mark.rests_on(f"{_THIS}::test_cp_albedo_below_one_is_refused")
@pytest.mark.parametrize("chart", _CHARTS)
def test_cp_white_at_unit_weight_and_vacuum_are_admitted(chart: str) -> None:
    r"""The control of gate 4 (X1: a refusal gate whose builder cannot succeed refuses everything): on the same charts and the same
    builder, ``white(1.0)`` (the weight at which CP's closure is the diffuse chain) and the vacuum build and ``solve_cp`` returns a
    finite eigenvalue. The spelling of the same physics as a shaped albedo, refused today as ``kind == "albedo"``, is admitted at weight
    1.0 by the fold of decision D1 (``white(1.0)`` IS the diffuse chain at unit weight).

    Red input: ``ImportError`` at the base tip. Mutation: a table that refuses ``white(1.0)`` reddens it.
    """
    import math

    W = H.bc().WallCondition
    for declared in (W.white(1.0), W.vacuum()):
        k = H.cp_solve(H.cp_mesh(chart), {0: declared}).keff
        assert math.isfinite(k) and k > 0.0, f"{chart}: k = {k!r}"
