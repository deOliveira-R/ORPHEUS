r"""W2 gates — the reified splitting matrix ``M = (L+C−B_lower)`` (#226 step 2).

Spec: ``issue_226_inverse_operator_verification.md`` §13.1–§13.3.  The
dissolved ``_GaussSeidelResolvent`` paired ``apply = (L+C)ψ`` with
``solve = (L+C−B_lower)⁻¹`` — inverses of DIFFERENT operators; its round-trip
defect was O(1) (2.667, §17 falsifier-3).  These gates pin the reification:
the forward and the inverse are now the two faces of ONE
:class:`~orpheus.sn.operators.scheduled_invertible.ScheduledInvertibleOperator`.

Domain note (implementation finding, recorded in the spec §13):
the sweep substrate re-derives the outflow-definition rows (``shed`` writes
``z.out = streamed``), so the walk realizes :math:`M^{-1}` exactly on the
SOURCE SUBSPACE ``{y : y.outflow-rows = 0}`` — which contains every
production rhs (``q + S·ψ + B_upper·ψ`` all write bulk/inflow only) — whose
:math:`M`-preimage is the TRACE-CONSISTENT states (``x.out =
streamed(x.interior)``, i.e. actual transport states; solve outputs by
construction).  The round-trip gate therefore round-trips a consistent
state and asserts machine precision on BOTH the bulk and the trace — a
STRONGER claim than the bulk-only falsifier.  The discriminator keeps its
teeth there: the confused pairing still fails at O(1) on a consistent state
(its apply lacks the ``B_lower`` subtraction entirely).

Mutation redefinition (spec §13 note): the masked ``B_lower`` single-sources
the row split for BOTH ``M.apply`` and the in-sweep reflect, so the spec'd
M-SPLIT ("mask ≠ walk fold") is UNREPRESENTABLE by construction.  Its
replacement pair, one per gate:

* **M-SPLIT-DIR** — flip the split direction (upper-as-lower) in
  ``SweepSchedule.lower_inflow_rows`` → the flipped rows are read BEFORE
  their face's reflect, so the late fold never reaches the reader →
  W2-round-trip REDs (exercised in
  :func:`test_mutation_split_direction_reddens_round_trip`).
* **M-SPLIT-PART** — doctor one half's rows after the split (partition no
  longer complementary) → W2-split REDs
  (:func:`test_mutation_partition_break_reddens_split`).
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.common.xs_library import get_mixture
from orpheus.geometry import BC
from orpheus.geometry.mesh import Mesh2D
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.loss_representation.sweep_schedule import (
    SweepSchedule,
    reflective_faces,
)
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.sn.operators.boundary import SNBoundaryOperator
from orpheus.sn.operators.scheduled_invertible import (
    ScheduledInvertibleOperator,
)
from orpheus.sn.operators.streaming import StreamingOperator
from orpheus.sn.operators.sweep_operator import SweepOperator
from orpheus.sn.solver import _select_si_splitting, solve_sn_fixed_source
from orpheus.transport.operators.multiplication_operator import (
    MultiplicationOperator,
)
from orpheus.transport.timed_full_field import TimedFullField
from tests.sn.operators.test_removal_form_matvec_sweep import (
    _cart2d,
    _random_state,
    _removal_sigmas,
)

pytestmark = pytest.mark.foundation


def _as_timed(field, like):
    return TimedFullField(
        interior=field.interior, boundary=field.boundary,
        _history=(), history_depth=like.history_depth,
    )


def _reified(bc: str = "reflective", *, seed: int = 7):
    """(LC, B, schedule, parts, M) on the het 2-D all-``bc`` mesh.

    All-reflective maximizes the ``B_lower`` support (every face reflected
    in-sweep); the removal σ set is heterogeneous 2G (config-blindness —
    a uniform σ nulls the streaming asymmetry the fold couples).
    """
    sn = _cart2d(nx=4, ny=5, ng=2, bc=bc)
    _, sig_r = _removal_sigmas(sn, seed=seed)
    LC = StreamingOperator.pose(sn) + MultiplicationOperator.from_mesh(sig_r, sn)
    B = SNBoundaryOperator(sn)
    schedule = SweepSchedule.gauss_seidel(
        sn.ndim, sn.quad.octants, reflective_faces(sn))
    parts = B.split(schedule)
    return sn, LC, B, schedule, parts, LC - parts.lower


def _consistent_state(sn, LC, *, seed: int):
    """A trace-CONSISTENT random state (``x.out = streamed(x.interior)``): a
    solve output on an INFLOW-ONLY rhs — the walk's honest domain (module
    docstring).

    ERR-071: the honest inverse emits ``ψ_out = streamed − rhs_out``, so a
    solve output is trace-consistent IFF the rhs carries ZERO outflow rows
    (exactly the production idiom).  Pre-fix the clobber made ANY solve
    output consistent — passing the raw random trace through relied on the
    dropped-row bug, so the outflow rows are zeroed here.
    """
    raw = _random_state(sn, seed=seed)
    trace = sn.angular_trace
    for face in raw.boundary.layout.faces:
        out_rows = trace.outflow_indices_for_face(face)
        raw.boundary.face_view(face)[out_rows] = 0.0
    return LC.solve(raw)


# ─────────────────────────────────────────────────────────────────────────
# §13.1 — W2-round-trip: machine precision on M's OWN forward.
# ─────────────────────────────────────────────────────────────────────────


def test_w2_round_trip_machine_precision():
    r"""``M.inverse().apply(M.apply(x)) == x`` at machine precision, bulk AND
    trace — what the dissolved resolvent could not do (2.667, O(1) RED).

    ≥2 reflective faces (here all four — maximal non-empty ``B_lower``),
    heterogeneous 2G removal σ, trace-consistent x (module docstring).
    De-risk: measured 5.2e-16 bulk / 4.4e-16 trace at authoring.
    """
    sn, LC, _B, _sched, parts, M = _reified()
    n_lower = sum(int(r.size) for r in parts.lower.rows.values())
    if n_lower == 0:
        pytest.fail(
            "B_lower support is EMPTY on an all-reflective 2-D mesh — the "
            "round-trip would be vacuously the Jacobi degenerate."
        )
    if not isinstance(M, ScheduledInvertibleOperator):
        pytest.fail(f"(L+C) - B_lower fused to {type(M).__name__}")

    x = _consistent_state(sn, LC, seed=13)
    inv = M.inverse()
    if not isinstance(inv, SweepOperator):
        pytest.fail(f"M.inverse() returned {type(inv).__name__}")

    rt = inv.apply(_as_timed(M.apply(x), x))
    np.testing.assert_allclose(
        np.asarray(rt.interior.values), np.asarray(x.interior.values),
        rtol=1e-10, atol=1e-12 * float(np.abs(x.interior.values).max()),
        err_msg="M⁻¹(Mx) ≠ x (bulk) — the reified G-S forward substitution "
                "does not invert its OWN forward",
    )
    for face in x.boundary.layout.faces:
        np.testing.assert_allclose(
            np.asarray(rt.boundary.face_view(face)),
            np.asarray(x.boundary.face_view(face)),
            rtol=1e-10, atol=1e-12,
            err_msg=f"M⁻¹(Mx) ≠ x on the {face} trace",
        )
    # The involution is an OBJECT-IDENTITY fact (taxonomy §13 I2).
    if inv.inverse() is not M:
        pytest.fail("M.inverse().inverse() is not M")


@pytest.mark.catches("ERR-071")
def test_off_domain_outflow_rhs_is_not_completed_yet():
    r"""Honest-scope witness: the SCHEDULED walk realizes ``M⁻¹`` on the
    source subspace ``{y : y.outflow-rows = 0}`` ONLY.

    A structurally nonzero outflow row leaves the lower-coupled inflow
    rows off by ``B(y_out)``: the mid-march reflect consumes UN-restored
    streamed values (the ERR-071 outflow-defect restore fires only after
    the whole march — correct for the bare unscheduled sweep, which IS
    exact on the whole space, but too late for a mid-march reader).  The
    full-space completion needs the restore interleaved per-group with
    the schedule — a walk-internal carve deferred until a full-space
    consumer exists (e.g. a G-S-preconditioned Krylov posture).

    This pin measures the documented gap so the scope claim is a fact,
    not prose.  THE TRIPWIRE: if the per-group completion lands, this
    test REDS — flip it to the exactness assert (round-trip ≡ y at
    machine precision) and lift the docstring restriction in
    :meth:`ScheduledInvertibleOperator._solve_timed_full_field`.
    """
    sn, LC, _B, _sched, _parts, M = _reified()
    x = _consistent_state(sn, LC, seed=13)
    y = M.apply(x)
    # Structurally populate the outflow rows (an off-domain rhs — the
    # role-conflated iterate-trace cast ERR-071 documents).
    trace = sn.angular_trace
    for face in y.boundary.layout.faces:
        out_rows = trace.outflow_indices_for_face(face)
        y.boundary.face_view(face)[out_rows] = 1.0
    z = M.inverse().apply(_as_timed(y, x))
    defect = M.apply(z)
    worst = 0.0
    for face in y.boundary.layout.faces:
        in_rows = trace.inflow_indices_for_face(face)
        worst = max(worst, float(np.abs(
            np.asarray(defect.boundary.face_view(face))[in_rows]
            - np.asarray(y.boundary.face_view(face))[in_rows]
        ).max()))
    if not worst > 0.1:
        pytest.fail(
            f"the scheduled walk's off-domain round-trip defect on the "
            f"lower-coupled inflow rows measured {worst:.2e} — the "
            f"per-group restore completion appears to have LANDED. "
            f"Flip this pin to the full-space exactness assert and lift "
            f"the documented source-subspace restriction."
        )


def test_factory_returns_reified_pair():
    r"""``_select_si_splitting(gauss_seidel)`` returns the splitting pair
    ``(M, (S, B_upper))`` — congruent with the Jacobi arm ``(LC, (S, B))``."""
    sn, LC, B, _sched, _parts, _M = _reified()
    sentinel_S = object()
    resolvent, gains = _select_si_splitting(
        LC, sentinel_S, B, sn, "gauss_seidel",
    )
    if not isinstance(resolvent, ScheduledInvertibleOperator):
        pytest.fail(f"G-S arm returned {type(resolvent).__name__}")
    np.testing.assert_equal(len(gains), 2)
    if gains[0] is not sentinel_S:
        pytest.fail("gains[0] must be the scattering gain, passed through")
    from orpheus.sn.operators.boundary import SNMaskedBoundaryOperator

    if not isinstance(gains[1], SNMaskedBoundaryOperator):
        pytest.fail(f"gains[1] must be B_upper; got {type(gains[1]).__name__}")
    # Jacobi arm unchanged: (LC, (S, B)).
    resolvent_j, gains_j = _select_si_splitting(
        LC, sentinel_S, B, sn, "jacobi",
    )
    if resolvent_j is not LC or gains_j != (sentinel_S, B):
        pytest.fail("the Jacobi arm must stay (LC, (S, B))")


# ─────────────────────────────────────────────────────────────────────────
# §13.2 — W2-split: B == B_lower + B_upper, bit-exact per face.
# ─────────────────────────────────────────────────────────────────────────


def test_w2_split_exactness():
    r"""The schedule split is an exact partition: ``B·ψ == B_lower·ψ +
    B_upper·ψ`` bit-identically per face (disjoint row writes; no
    arithmetic change).  The specular map has no octant-diagonal."""
    sn, _LC, B, _sched, parts, _M = _reified()
    psi = _random_state(sn, seed=13)  # raw random — a pure partition claim
    whole = B.apply(psi)
    summed = parts.lower.apply(psi) + parts.upper.apply(psi)
    for face in psi.boundary.layout.faces:
        np.testing.assert_array_equal(
            np.asarray(whole.boundary.face_view(face)),
            np.asarray(summed.boundary.face_view(face)),
            err_msg=f"B ≠ B_lower + B_upper on face {face} — the octant "
                    f"split is not an exact partition",
        )
    # Row supports are complementary within each face's inflow.
    trace = sn.angular_trace
    for face in psi.boundary.layout.faces:
        lo = np.asarray(parts.lower.rows.get(face, np.empty(0, dtype=np.intp)))
        up = np.asarray(parts.upper.rows[face])
        inflow = trace.inflow_indices_for_face(face)
        np.testing.assert_array_equal(
            np.sort(np.concatenate([lo, up])), np.sort(inflow),
            err_msg=f"lower ∪ upper ≠ inflow rows on {face}",
        )
        np.testing.assert_equal(np.intersect1d(lo, up).size, 0)


def test_jacobi_schedule_split_is_degenerate():
    r"""The Jacobi schedule's lower support is EMPTY (``B_lower = 0``,
    ``B_upper = B``) — the degenerate that recovers the plain lagged-B
    iteration."""
    sn, _LC, B, _sched, _parts, _M = _reified()
    jac = SweepSchedule.jacobi(sn.ndim, sn.quad.octants)
    np.testing.assert_equal(jac.lower_inflow_rows(sn), {})
    parts = B.split(jac)
    psi = _random_state(sn, seed=5)
    lower_out = parts.lower.apply(psi)
    whole = B.apply(psi)
    upper_out = parts.upper.apply(psi)
    for face in psi.boundary.layout.faces:
        np.testing.assert_array_equal(
            np.asarray(lower_out.boundary.face_view(face)), 0.0,
            err_msg=f"Jacobi B_lower must be zero on {face}",
        )
        np.testing.assert_array_equal(
            np.asarray(upper_out.boundary.face_view(face)),
            np.asarray(whole.boundary.face_view(face)),
            err_msg=f"Jacobi B_upper must equal B on {face}",
        )


# ─────────────────────────────────────────────────────────────────────────
# §13.3 — W2-FP: fixed-point equivalence vs Jacobi, Mode-9-SAFE config.
# ─────────────────────────────────────────────────────────────────────────


def test_w2_fixed_point_equivalence_diagonal_cubature():
    r"""The converged fixed point is splitting-invariant (``vv`` Mode 9) —
    G-S ≡ Jacobi to solver tolerance on a config that BREAKS the degenerate
    coincidences: a DIAGONAL cubature (level-symmetric — shared faces,
    ERR-056's regime; an axis-aligned product quad makes octant-G-S
    accidentally exact) on a HETEROGENEOUS 2-D box with vacuum-x /
    reflective-y (anisotropic flux — the fully-reflective isotropic box is
    the Mode-9 degenerate).

    NECESSARY but not sufficient (§13.3): the fixed point cannot even
    distinguish G-S from Jacobi — §13.1 is the load-bearing correctness
    gate; this pins the SPLITTING claim (same ψ*, only the rate differs).
    """
    fuel, mod = get_mixture("A", "2g"), get_mixture("B", "2g")
    nx, ny = 8, 4
    mat = np.zeros((nx, ny), dtype=int)
    mat[:4, :] = 2
    mesh = Mesh2D(
        edges_x=np.linspace(0.0, 2.0, nx + 1),
        edges_y=np.linspace(0.0, 1.0, ny + 1),
        mat_map=mat,
        bc_xmin=BC("vacuum"), bc_xmax=BC("vacuum"),
        bc_ymin=BC("reflective"), bc_ymax=BC("reflective"),
    )
    quad = Quadrature.level_symmetric(sn_order=4)
    q_ext = np.full((quad.N, 2, nx, ny), 1.0 / float(quad.weights.sum()))
    inner_tol = 1e-12
    kwargs = dict(
        scattering_order=1, inner_solver="source_iteration",
        max_inner=3000, inner_tol=inner_tol,
    )
    materials = {2: fuel, 0: mod}
    phi_gs = solve_sn_fixed_source(
        materials, mesh, quad, q_ext, inner_schedule="gauss_seidel", **kwargs,
    ).scalar_flux
    phi_jac = solve_sn_fixed_source(
        materials, mesh, quad, q_ext, inner_schedule="jacobi", **kwargs,
    ).scalar_flux
    # SAFETY × conv_tol (L7): the two iterations stop within inner_tol of
    # the SAME ψ*; 50× covers the stopping-criterion slack on both sides.
    np.testing.assert_allclose(
        np.asarray(phi_gs.values), np.asarray(phi_jac.values),
        rtol=50 * inner_tol,
        err_msg="G-S and Jacobi converged to DIFFERENT fixed points — the "
                "splitting moved ψ* (vv Mode 9; cf. ERR-056)",
    )


# ─────────────────────────────────────────────────────────────────────────
# Mutations — each gate's teeth (monkeypatch only; -O-proof asserts).
# ─────────────────────────────────────────────────────────────────────────


def test_mutation_split_direction_reddens_round_trip(monkeypatch):
    r"""M-SPLIT-DIR: flip the split direction (upper-as-lower).  The flipped
    rows are read BEFORE their face's reflect fires, so the in-sweep fold
    never reaches the reader while ``apply`` subtracts it — the round-trip
    defect is O(1).  Proves §13.1 bites the split-direction convention
    (Mode-1-family: a ``>`` vs ``<`` flip)."""
    real = SweepSchedule.lower_inflow_rows

    def flipped(self, sn_mesh):
        rows = real(self, sn_mesh)
        trace = sn_mesh.angular_trace
        return {
            face: np.setdiff1d(trace.inflow_indices_for_face(face), lo)
            for face, lo in rows.items()
        }

    monkeypatch.setattr(SweepSchedule, "lower_inflow_rows", flipped)
    sn, LC, _B, _sched, parts, M = _reified()
    if sum(int(r.size) for r in parts.lower.rows.values()) == 0:
        pytest.fail("flipped split produced an empty lower — vacuous mutation")
    x = _consistent_state(sn, LC, seed=13)
    rt = M.inverse().apply(_as_timed(M.apply(x), x))
    defect = float(
        np.abs(np.asarray(rt.interior.values) - np.asarray(x.interior.values)).max()
        / np.abs(np.asarray(x.interior.values)).max()
    )
    if defect < 1e-3:
        pytest.fail(
            f"split-direction flip NOT caught: round-trip defect {defect:.2e} "
            f"— §13.1 has no teeth on the split convention"
        )


def test_mutation_partition_break_reddens_split():
    r"""M-SPLIT-PART: doctor one half's rows post-split (drop one face's
    rows from ``B_lower`` without re-complementing) — the partition claim
    ``B == B_lower + B_upper`` REDs.  Proves §13.2 bites."""
    sn, _LC, B, _sched, parts, _M = _reified()
    doctored_rows = dict(parts.lower.rows)
    face_with_rows = next(
        f for f, r in doctored_rows.items() if np.asarray(r).size
    )
    doctored_rows[face_with_rows] = np.empty(0, dtype=np.intp)
    from orpheus.sn.operators.boundary import SNMaskedBoundaryOperator

    doctored = SNMaskedBoundaryOperator(B, doctored_rows, parts.lower.schedule)
    psi = _random_state(sn, seed=13)
    whole = np.asarray(B.apply(psi).boundary.face_view(face_with_rows))
    summed = np.asarray(
        (doctored.apply(psi) + parts.upper.apply(psi))
        .boundary.face_view(face_with_rows)
    )
    if np.array_equal(whole, summed):
        pytest.fail(
            "partition break NOT caught: dropping a face's lower rows left "
            "B == B_lower + B_upper — §13.2 has no teeth"
        )
