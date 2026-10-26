r"""Gates of the SN consumers of the conditions package, posed through ``SNProblem`` and the ``solve_sn*`` entries (wave 2, S2b).

Spec: ``.claude/plans/wave2_boundary_carve.md`` sections 2.10 and 6.2 (S2b gates 2, 5, 6, 7, 8, 10 and 12) and section 4.3
(deltas 4 and 6). Every problem is spelled in the NEW vocabulary (``tests/gates/sn/architecture/_sn_consumers_helpers.py``):
a face carries a role, an ``int`` tag or a deck, and the conditions are a ``{tag: WallCondition}`` mapping handed to the
entry as ``boundary_conditions=``.

First red at the base tip of slice S2b, for every row: ``ImportError`` from ``import orpheus.sn`` (``orpheus.transport.method``
imports the retired ``orpheus.geometry.BC``) before any row can read its own line; the per-row reds are described in each
docstring as the INPUT that reddens it once the carve lands, and as the mutation the reviewers apply.

The ladder::

    test_bc_holds_realized_boundaries_and_the_binding   foundation: the problem's surface (bc, binding, reflective_axes)
    test_the_empty_default_is_safe                      foundation: no condition is filled in for the caller (decision 6)
        ^
    test_the_binding_through_sn                         edge: law L-C through the real binding
    test_a_certification_failure_names_its_axis         edge: ANGULAR_RESOLUTION at construction
    test_curvilinear_refusals_are_at_admit              edge: the radial corner action is refused AT admit (#653)
    test_the_hollow_inner_witness                       edge: the cavity (carry) and its refusals (#511)
        ^
    test_the_vacuum_question_is_the_realized_map        interior: DSA reads the type
    test_the_adjoint_ignores_q                          interior: law L-A (iii)
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pytest

from tests.gates.sn.architecture import _sn_consumers_helpers as S
from tests.gates.transport import _carry_cases as carry

_THIS = "tests/gates/sn/architecture/test_sn_boundary_binding.py"
_CARRY = "tests/gates/sn/architecture/test_wave2_carry.py"

pytestmark = pytest.mark.l1


def _one_material() -> dict[int, Any]:
    return {0: carry.materials()[1]}


def _slab(left: Any, right: Any, *, mat_ids: tuple[int, ...] = (0, 0), breakpoints: tuple[float, ...] = (0.0, 1.0, 2.0), cells: tuple[int, ...] = (3, 3)) -> Any:
    return S.mesh_of("slab", left=left, right=right, mat_ids=mat_ids, breakpoints=breakpoints, cells=cells)


# ═════════════════════════════════════════════════════════════════════════════
# The problem's surface
# ═════════════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
def test_bc_holds_realized_boundaries_and_the_binding() -> None:
    r"""``SNProblem.bc`` holds one ``RealizedBoundary`` per face with a trace (a vacuum tag is a ``ZeroMap``, a mirror
    deck a ``ScaledOrdinatePermutation`` with ``alpha = 1`` that permutes ordinates, a diffuse tag a ``RankOneReemission``);
    ``SNProblem.binding`` is the ``Binding`` of the carrier; ``reflective_axes`` reads the DECKS of the binding (a slab
    mirrored on both faces is reflective on its axis; vacuum tags, and a deck on one face only, are not); a hollow body's
    inner face is a tag with no trace and so no ``bc`` entry; the class carries no ``BOUNDARY_OPERATOR_REGISTRY``.

    Red input at the base tip: ``problem.bc`` holds the ``_BoundBoundaryOperator`` shim and ``reflective_axes`` reads
    ``isinstance(bound.law, ReflectiveBoundary)``. Mutation: ``reflective_axes`` reading the tags (a ``specular(1.0)``
    pair would then be reflective); ``bc`` wrapping the operators in the retired shim.
    """
    from orpheus.geometry.deck import MirrorDeck

    B = S.bc()
    from tests.gates.sn.architecture import _realized_boundary_helpers as H

    R = H.ro()
    quad = S.quadrature_of("slab", 8)
    mats = _one_material()
    vac = S.problem(_slab(0, 0), quad, {0: B.WallCondition.vacuum()}, mats)
    assert set(vac.bc) == {"xmin", "xmax"} and all(isinstance(op, R.ZeroMap) for op in vac.bc.values())
    assert isinstance(vac.binding, B.Binding) and set(vac.binding.tags) == {"xmin", "xmax"} and not vac.binding.decks
    assert vac.reflective_axes == ()

    mirrored = S.problem(_slab(MirrorDeck(), MirrorDeck()), quad, None, mats)
    assert all(isinstance(op, R.ScaledOrdinatePermutation) and op.return_fraction == 1.0 and op.permutes_ordinates for op in mirrored.bc.values())
    assert set(mirrored.binding.decks) == {"xmin", "xmax"} and mirrored.reflective_axes == (0,)

    half = S.problem(_slab(MirrorDeck(), 0), quad, {0: B.WallCondition.vacuum()}, mats)
    assert half.reflective_axes == (), "a deck on one face only does not make the axis reflective"

    diffuse = S.problem(_slab(0, 0), quad, {0: B.WallCondition.white(0.5)}, mats)
    assert all(isinstance(op, R.RankOneReemission) for op in diffuse.bc.values())

    hollow = S.problem(S.mesh_of("hollow_cylinder", outer=0, inner=1), S.quadrature_of("hollow_cylinder"), {0: B.WallCondition.vacuum(), 1: B.WallCondition.specular(1.0)})
    assert set(hollow.bc) == {"xmax"}, "the hollow inner face has no trace and no bc entry"
    assert set(hollow.binding.tags) == {"xmin", "xmax"} and hollow.binding.tags["xmin"] == 1
    for problem in (vac, mirrored, hollow):
        assert not hasattr(type(problem), "BOUNDARY_OPERATOR_REGISTRY") and not hasattr(problem, "realize_boundary_law")


@pytest.mark.foundation
def test_the_empty_default_is_safe() -> None:
    r"""Decision 6's witness: no condition is filled in for the caller. A TAGGED SN problem built with no
    ``boundary_conditions`` raises ``UnboundTagError`` naming the tag and the face; each tag is named in its turn (the
    other declared); a declaration for a tag no face carries raises ``UnusedConditionError``; a DECK-ONLY problem (the PWR
    half cell, mirror decks on both faces) constructs with none.

    Red input at the base tip: ``boundary_conditions`` is not a keyword of ``SNProblem`` and the slab's laws ride the mesh.
    Mutation: an empty default that binds every tag to ``vacuum()`` (the silent fill-in of a wall).
    """
    B = S.bc()
    quad = S.quadrature_of("slab", 8)
    mats = _one_material()
    mesh = _slab(0, 1)
    with pytest.raises(B.UnboundTagError) as info:
        S.problem(mesh, quad, None, mats)
    assert info.value.tag in (0, 1)
    with pytest.raises(B.UnboundTagError) as only_left:
        S.problem(mesh, quad, {1: B.WallCondition.vacuum()}, mats)
    assert (only_left.value.tag, only_left.value.face) == (0, "xmin")
    with pytest.raises(B.UnboundTagError) as only_right:
        S.problem(mesh, quad, {0: B.WallCondition.vacuum()}, mats)
    assert (only_right.value.tag, only_right.value.face) == (1, "xmax")
    S.problem(mesh, quad, {0: B.WallCondition.vacuum(), 1: B.WallCondition.vacuum()}, mats)
    with pytest.raises(B.UnusedConditionError) as unused:
        S.problem(mesh, quad, {0: B.WallCondition.vacuum(), 1: B.WallCondition.vacuum(), 5: B.WallCondition.vacuum()}, mats)
    assert unused.value.tag == 5

    from orpheus.data.macro_xs.synthetic import get_mixture
    from orpheus.geometry import StructuredGeometry
    from orpheus.mesh import CellsByCount, Mesher

    half_cell = Mesher(StructuredGeometry.pwr_slab_half_cell()).partition(tuple(CellsByCount.uniform_volume(3) for _ in range(3))).mesh
    deck_only = S.problem(half_cell, S.quadrature_of("slab", 8), None, {2: get_mixture("A", "2g"), 1: get_mixture("B", "2g"), 0: get_mixture("C", "2g")})
    assert not deck_only.binding.tags and set(deck_only.binding.decks) == {"xmin", "xmax"}


# ═════════════════════════════════════════════════════════════════════════════
# Gate 2: law L-C through SN's real binding
# ═════════════════════════════════════════════════════════════════════════════


@pytest.mark.rests_on(f"{_THIS}::test_the_empty_default_is_safe")
def test_the_binding_through_sn() -> None:
    r"""Law L-C: a lab-frame source term on a symmetric two-face tag raises ``NonInvariantConditionError`` through
    ``SNProblem``; the controls bind. The authored first red: the symmetric slab ``(0, 1, 2)`` of material ``(0, 0)``, tag 0
    on both faces, ``WallCondition.framed(LabFrameInflow(f))``. Controls that bind: the same with a uniform source or a
    face-frame source (``face_local`` implies invariant on a 1-D slab), an ASYMMETRIC slab (two materials), tags 0 and 1
    (one face each), a one-face tag (a deck on the other face).

    Red input at the base tip: the equivalent typed ``PrescribedInflow`` on both faces is realized silently. Mutation:
    ``invariant_under`` replaced by ``True`` in a copy of ``binding.py``.
    """
    from orpheus.geometry.deck import MirrorDeck

    from tests.gates.transport._boundary_conditions_helpers import ProfileSampler

    B = S.bc()
    quad = S.quadrature_of("slab", 8)
    lab = B.WallCondition.framed(B.LabFrameInflow(ProfileSampler(2.0)))
    face = B.WallCondition.framed(B.FaceFrameInflow(ProfileSampler(2.0)))
    one = _one_material()
    both = carry.materials()

    with pytest.raises(B.NonInvariantConditionError) as info:
        S.problem(_slab(0, 0), quad, {0: lab}, one)
    assert info.value.tag == 0 and set(info.value.faces) == {"xmin", "xmax"}

    S.problem(_slab(0, 0), quad, {0: B.WallCondition.uniform(1.0)}, one)
    S.problem(_slab(0, 0), quad, {0: face}, one)
    S.problem(_slab(0, 0, mat_ids=(0, 1)), quad, {0: lab}, both)
    S.problem(_slab(0, 1), quad, {0: lab, 1: lab}, one)
    S.problem(_slab(0, MirrorDeck()), quad, {0: lab}, one)


# ═════════════════════════════════════════════════════════════════════════════
# Gate 5: a certification failure names its axis, at construction
# ═════════════════════════════════════════════════════════════════════════════


@pytest.mark.rests_on(f"{_THIS}::test_bc_holds_realized_boundaries_and_the_binding")
def test_a_certification_failure_names_its_axis() -> None:
    r"""An odd-``n_phi`` product rule has no x-mirror closure: a mirror deck on a 2-D Cartesian mesh raises
    ``BoundaryRefused`` whose admission is ``Refused(ANGULAR_RESOLUTION)`` at ``SNProblem`` construction, the face named;
    the even rule of the same shape certifies (the control), so a refusal that fired on every product rule cannot pass.

    Red input at the base tip: the construction fails with the ``BoundaryError`` of the realizer's own certification
    (``law=...``), not the typed admission. Mutation: a table whose ``admit`` skips the pairing (the error moves to ``realize``).
    """
    from orpheus.geometry.deck import MirrorDeck
    from orpheus.mesh import Mesh2D

    B = S.bc()
    mats = _one_material()

    def mesh2d() -> Any:
        return Mesh2D(np.array([0.0, 1.0, 2.0]), np.array([0.0, 1.0, 2.0]), np.zeros((2, 2), dtype=int), face_roles={f: MirrorDeck() for f in ("xmin", "xmax", "ymin", "ymax")})

    from tests.gates.sn.architecture import _realized_boundary_helpers as H

    Q = H.quadrature_class()
    S.problem(mesh2d(), Q.product(n_mu=2, n_phi=4), None, mats)
    with pytest.raises(B.BoundaryRefused) as info:
        S.problem(mesh2d(), Q.product(n_mu=2, n_phi=3), None, mats)
    assert info.value.admission.axis is B.RefusalAxis.ANGULAR_RESOLUTION and info.value.face in {"xmin", "xmax"}


# ═════════════════════════════════════════════════════════════════════════════
# Gate 6: the curvilinear refusals are at admit
# ═════════════════════════════════════════════════════════════════════════════

_CURVED = ("cylinder", "sphere")
#: the declarations a radial outer face refuses at admit (the corner action is ruled only for a permuting wall or a zero
#: return, #653): the diffuse chain, a mix with it, a source alone, a response over a source.
_REFUSED_ON_A_RADIAL_OUTER = {
    "white 0.5": lambda w: w.white(0.5),
    "both 0.3 0.5": lambda w: w.from_amplitudes(0.3, 0.5),
    "inflow 1.0": lambda w: w.uniform(1.0),
    "specular 0.5 over inflow 1.0": lambda w: w.specular(0.5) + w.uniform(1.0),
}
_ADMITTED_ON_A_RADIAL_OUTER = {"vacuum": lambda w: w.vacuum(), "specular 0.7": lambda w: w.specular(0.7), "specular 1.0": lambda w: w.specular(1.0)}


@pytest.mark.rests_on(f"{_THIS}::test_bc_holds_realized_boundaries_and_the_binding")
@pytest.mark.parametrize("label", list(_REFUSED_ON_A_RADIAL_OUTER))
@pytest.mark.parametrize("body", _CURVED)
def test_curvilinear_refusals_are_at_admit(body: str, label: str) -> None:
    r"""A cylinder or a sphere with a diffuse, mixed or sourced outer wall raises ``BoundaryRefused(PLUMBING, #653)`` at
    CONSTRUCTION. The vacuum, specular 0.7 and mirror walls of the same body construct and solve (the controls, so a
    construction that refused every curved body cannot pass). Section 4.3 delta 4.

    Red input at the base tip (``[M]`` the plan's measurement, 12 of 12 declarations over the two bodies): the
    construction succeeds and ``solve_sn`` raises ``NotImplementedError`` at the first sweep. Mutation: ``SNDiffuseChain.admit``
    returns ``Realized`` on a radial outer site (reddens this and the SN rows of ``test_admission_matrix.py``).
    """
    B = S.bc()
    mesh = S.mesh_of(body, outer=0)
    quad = S.quadrature_of(body)
    for control in _ADMITTED_ON_A_RADIAL_OUTER.values():
        S.problem(mesh, quad, {0: control(B.WallCondition)})
    with pytest.raises(B.BoundaryRefused) as info:
        S.problem(mesh, quad, {0: _REFUSED_ON_A_RADIAL_OUTER[label](B.WallCondition)})
    assert info.value.admission.axis is B.RefusalAxis.PLUMBING and info.value.admission.issue == 653 and info.value.face == "xmax"


# ═════════════════════════════════════════════════════════════════════════════
# Gate 10: the hollow inner witness
# ═════════════════════════════════════════════════════════════════════════════


@pytest.mark.rests_on(f"{_THIS}::test_bc_holds_realized_boundaries_and_the_binding", f"{_CARRY}::test_the_carry_of_the_solves")
def test_the_hollow_inner_witness() -> None:
    r"""A hollow cylinder with the inner tag bound to ``WallCondition.specular(1.0)`` builds and solves to the SAME ``k``
    as today's reflective cavity (the carry fixture's ``hollow_cylinder`` entries, bit-equal), and an inner tag bound to a
    vacuum, white, partial specular or sourced condition raises ``BoundaryRefused(PLUMBING, #511)`` at construction with the
    inner face named (no trace and so no ``realize``: ``table.require`` is the only door there).

    Red input at the base tip: the refusal is ``NotImplementedError`` from the adapter (15 of 15 constructible
    non-reflective typed laws). Mutation: the inner site admitting every weight (the cavity is the closure's own).
    """
    B = S.bc()
    fixture = np.load(_carry_path())
    mesh = S.mesh_of("hollow_cylinder", outer=0, inner=1)
    quad = S.quadrature_of("hollow_cylinder")
    cavity = {0: B.WallCondition.vacuum(), 1: B.WallCondition.specular(1.0)}
    S.problem(mesh, quad, cavity)

    from orpheus.sn.solver import solve_sn

    solved = solve_sn(carry.materials(), mesh, quad, boundary_conditions=cavity, **carry.EIGEN)
    assert solved.outcome.keff == float(fixture["keff/hollow_cylinder/vacuum"])

    for label, inner in (("vacuum", B.WallCondition.vacuum()), ("white 0.5", B.WallCondition.white(0.5)), ("specular 0.5", B.WallCondition.specular(0.5)), ("inflow", B.WallCondition.uniform(1.0))):
        with pytest.raises(B.BoundaryRefused) as info:
            S.problem(mesh, quad, {0: B.WallCondition.vacuum(), 1: inner})
        assert info.value.admission.axis is B.RefusalAxis.PLUMBING and info.value.admission.issue == 511 and info.value.face == "xmin", label


def _carry_path() -> str:
    return str(carry.CARRY_FIXTURE)


# ═════════════════════════════════════════════════════════════════════════════
# Gate 7: the vacuum question is the realized map
# ═════════════════════════════════════════════════════════════════════════════


@pytest.mark.rests_on(f"{_THIS}::test_bc_holds_realized_boundaries_and_the_binding")
def test_the_vacuum_question_is_the_realized_map() -> None:
    r"""Attack A1: ``WallCondition.specular(0.0)``, ``WallCondition.white(0.0)`` and ``WallCondition.vacuum()`` are ONE
    value (equal, one digest), realize a ``ZeroMap``, ``is_vacuum`` is True, and the DSA admits the face and builds the
    Marshak row: the three low-order systems are ``array_equal``. A unit mirror gets the (39) row (a different system),
    a partial specular wall (ERR-094) and a diffuse wall are refused (``NotImplementedError``): the DSA reads
    ``face.is_vacuum or (face.permutes_ordinates and face.return_fraction == 1.0)`` off the TYPE.

    Red input at the base tip: ``AlbedoBoundary(0, S)`` and ``VacuumInflow()`` are two values and the DSA asks
    ``law.response_kernel.is_zero``. Mutation: the DSA admitting on ``permutes_ordinates`` alone (a 0.5 mirror would build
    the mirror's row, ERR-094).
    """
    from orpheus.sn.acceleration.dsa import DSALowOrderSystem
    from orpheus.values.content import content_digest

    B = S.bc()
    W = B.WallCondition
    spellings = [W.specular(0.0), W.white(0.0), W.vacuum()]
    assert spellings[0] == spellings[1] == spellings[2] and len({content_digest(x) for x in spellings}) == 1

    mesh = S.mesh_of("slab", left=0, right=0)
    quad = S.quadrature_of("slab", 8)
    systems = []
    for wall in spellings:
        problem = S.problem(mesh, quad, {0: wall})
        assert all(op.is_vacuum for op in problem.bc.values())
        systems.append(DSALowOrderSystem.from_problem(problem))
    np.testing.assert_array_equal(systems[0].a_low, systems[1].a_low)
    np.testing.assert_array_equal(systems[1].a_low, systems[2].a_low)

    mirror = DSALowOrderSystem.from_problem(S.problem(mesh, quad, {0: W.specular(1.0)}))
    assert not np.array_equal(mirror.a_low, systems[0].a_low), "the unit mirror's row (39) is not Marshak's (38)"
    for refused in (W.specular(0.5), W.white(1.0), W.white(0.5)):
        with pytest.raises(NotImplementedError):
            DSALowOrderSystem.from_problem(S.problem(mesh, quad, {0: refused}))


# ═════════════════════════════════════════════════════════════════════════════
# Gate 8: the adjoint ignores q
# ═════════════════════════════════════════════════════════════════════════════


@pytest.mark.rests_on(f"{_THIS}::test_the_binding_through_sn")
def test_the_adjoint_ignores_q() -> None:
    r"""Law L-A (iii): ``q`` is not adjointed. ``solve_sn_adjoint_fixed_source`` returns the SAME adjoint (importance)
    solution whatever source is bound (a vacuum face, ``uniform(1.0)``, ``uniform(50.0)``), while ``solve_sn`` and
    ``solve_sn_adjoint`` refuse a bound source (the eigenvalue question has no source term, ``_refuse_a_boundary_source``,
    ERR-095); the source-free declaration solves (the control). HEAD already agrees on the first leg; the gate pins the
    silent asymmetry.

    Red input at the base tip: none for the first leg; the refusal legs read the retired ``law.source``. Mutation: an
    adjoint entry that adjoints the boundary source into the right-hand side.
    """
    from orpheus.sn.solver import solve_sn, solve_sn_adjoint, solve_sn_adjoint_fixed_source

    B = S.bc()
    W = B.WallCondition
    mats = carry.materials()
    quad = S.quadrature_of("slab", 8)
    mesh = S.mesh_of("slab", left=0, right=1)
    ng, nx = 2, sum(carry.BODIES["slab"].cells)
    detector = np.zeros((ng, nx))
    detector[1, 2:5] = 0.7

    flux = []
    for q in (W.vacuum(), W.uniform(1.0), W.uniform(50.0)):
        adjoint = solve_sn_adjoint_fixed_source(mats, mesh, quad, detector, boundary_conditions={0: q, 1: W.vacuum()}, inner_tol=1e-12, max_inner=4000)
        flux.append(np.asarray(adjoint.scalar_flux.values))
    np.testing.assert_array_equal(flux[0], flux[1])
    np.testing.assert_array_equal(flux[1], flux[2])

    for entry in (solve_sn, solve_sn_adjoint):
        with pytest.raises(ValueError, match="eigenvalue question"):
            entry(mats, mesh, quad, boundary_conditions={0: W.uniform(1.0), 1: W.vacuum()})
        assert np.isfinite(entry(mats, mesh, quad, boundary_conditions={0: W.vacuum(), 1: W.vacuum()}).outcome.keff)
