r"""The composite sweep-inverse identity — ``(L+C) ∘ (L+C)⁻¹ ≡ I``
on the FULL composite space, outflow-trace rows included (ERR-071).

The forward's boundary block is a sibling of the bulk: the outflow-trace
row is the self-consistency DEFECT ``streamed − ψ_out`` (the
authoritative forward's row, loss_representation) and the inflow
row the identity on the given inflow. The exact inverse must therefore
emit ``ψ_out = streamed − rhs_out`` (the sign is pinned by this gate's
round-trip) — and until 2026-07-26 the sweep
DROPPED the rhs's outflow-row content (seeded it into the mutable
boundary buffer, then let the march clobber it). Every PHYSICAL rhs
carries zero there (builders populate inflow slots only; outflow rows
are 0 = 0 identities at the fixed point), so all SI/eigenvalue paths
were blind — but a GMRES preconditioner ``M = (I + 𝒞)∘(L+C)⁻¹``
exercises the full composite space, where the dropped term made M
SINGULAR on the outflow-trace subspace (measured ‖M q‖/‖q‖ = 1e-15 on
a pure outflow-row Krylov residual; GMRES stalled at an O(1) true
residual and the end-of-solve claim check refused the claim — the
catch that exposed the class). The P1-DSA (d₁) Krylov posture (#2)
excited it deterministically.

This gate pins the identity the exact-inverse pair owes: round-trip a
RANDOM composite (interior + boundary, all rows live), and the pure
outflow-row leg that was the singular subspace.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.numerics.coupled_system import CoupledField, CoupledOperator
from orpheus.transport.radial_characteristic_field import (
    RadialCharacteristicField,
)
from orpheus.sn.splitting import Splitting, resolve_schedule
from orpheus.sn.coupled_system import build_within_group_system
from orpheus.transport.fields.angular_boundary_flux import AngularBoundaryFlux
from orpheus.transport.fields.angular_flux import AngularFlux
from orpheus.transport.full_field import FullField
from tests.gates.sn.operators._full_space_states import (
    GEOMS as _GEOMS,
    MESHES as _MESHES,
    random_state as _random_state,
    system_a as _system_a,
    zero_source_composite as _zero_source_composite,
)

pytestmark = [
    pytest.mark.foundation,
    pytest.mark.catches("ERR-071"),
    # ERR-078: the ψ½ march's solve dropped the outflow-row rhs — the
    # System-B twin of ERR-071, caught by this file's coupled rows.
    pytest.mark.catches("ERR-078"),
]


def _lc_pair(geom: str):
    """The production within-group forward + its exact inverse on a
    small heterogeneous mesh (every trace row is live: inflow rows are
    identities, outflow rows are defects).  Slabs carry the bare
    ``(L+C)``; the cylinder — carrying since the 6.3 flip — carries the
    upper-triangular coupled ``[[LC, Seeding], [None, march]]``, whose
    ``inverse()`` is the block back-substitution."""
    problem = _MESHES[geom]()
    system = build_within_group_system(
        problem, problem.mat_xs,
    )
    lc = Splitting.from_schedule(system, resolve_schedule(problem, "jacobi")).implicit
    if geom in ("cyl_folded", "sphere_gl"):
        if not isinstance(lc, CoupledOperator):
            pytest.fail(
                f"{geom}: a carrying mesh's implicit operator must be "
                f"the coupled (bulk ⊕ ψ½) composite"
            )
    elif isinstance(lc, CoupledOperator):
        pytest.fail(
            f"{geom}: a non-carrying mesh must carry the bare (L+C) arm"
        )
    return problem, lc, lc.inverse()


class TestSweepInverseIdentity:
    @pytest.mark.verifies("sn-dsa-sweep-inverse-identity")
    @pytest.mark.parametrize("geom", _GEOMS)
    def test_forward_of_inverse_is_identity_on_a_random_composite(
        self, geom,
    ):
        """``(L+C)((L+C)⁻¹ rhs) ≡ rhs`` with EVERY block populated —
        bulk, inflow-trace, and the outflow-trace rows the old sweep
        dropped.

        Honest scope on the trace: the identity is claimed on the
        inflow ∪ outflow rows.  A DEGENERATE pure-azimuthal ordinate
        (``μ_r = 0``, the folded n_φ ≡ 2 (mod 4) rule — excluded from
        BOTH selectors) has NO streaming coupling to the face: its
        trace slot is a free DOF of the composite (#284), where the
        forward is a structural ZERO row and the inverse completes
        with the identity (seed passthrough).  Both halves of that
        pair are asserted explicitly; on a slab the degenerate set is
        empty and the claim is the full-trace identity.  On the
        carrying cylinder the round-trip runs the COUPLED composite —
        the identity is additionally claimed on the ψ½ System-B
        block."""
        problem, lc, sweep = _lc_pair(geom)
        trace = problem.angular_trace
        rhs = _random_state(problem, lc, seed=17)
        psi = sweep.apply(rhs)
        back = lc.apply(psi)
        rhs_a, psi_a, back_a = _system_a(rhs), _system_a(psi), _system_a(back)
        np.testing.assert_allclose(
            np.asarray(back_a.interior.values),
            np.asarray(rhs_a.interior.values),
            rtol=1e-12, atol=1e-12,
            err_msg=f"{geom}: (L+C)∘(L+C)⁻¹ must be the identity on "
                    f"the bulk",
        )
        if isinstance(rhs, CoupledField):
            np.testing.assert_allclose(
                np.asarray(back.systems[1].to_flat()),
                np.asarray(rhs.systems[1].to_flat()),
                rtol=1e-12, atol=1e-12,
                err_msg=f"{geom}: the coupled round-trip must be the "
                        f"identity on the ψ½ System-B block",
            )
        n_live = 0
        n_degenerate = 0
        for face in rhs_a.boundary.layout.faces:
            live = np.union1d(
                trace.inflow_indices_for_face(face),
                trace.outflow_indices_for_face(face),
            )
            degenerate = np.setdiff1d(
                np.arange(problem.quad.N), live,
            )
            n_live += live.size
            n_degenerate += degenerate.size
            np.testing.assert_allclose(
                np.asarray(back_a.boundary.face_view(face))[live],
                np.asarray(rhs_a.boundary.face_view(face))[live],
                rtol=1e-12, atol=1e-12,
                err_msg=f"{geom}/{face}: (L+C)∘(L+C)⁻¹ must be the "
                        f"identity on the live trace — inflow "
                        f"identities AND outflow defect rows",
            )
            if degenerate.size:
                np.testing.assert_allclose(
                    np.asarray(back_a.boundary.face_view(face))[degenerate],
                    0.0, atol=1e-12,
                    err_msg=f"{geom}/{face}: the forward must be a "
                            f"structural zero row on the degenerate "
                            f"(μ_r = 0) free-DOF trace slots",
                )
                np.testing.assert_allclose(
                    np.asarray(psi_a.boundary.face_view(face))[degenerate],
                    np.asarray(rhs_a.boundary.face_view(face))[degenerate],
                    rtol=1e-12, atol=1e-12,
                    err_msg=f"{geom}/{face}: the inverse must complete "
                            f"the free-DOF slots with the identity "
                            f"(seed passthrough)",
                )
        if not n_live > 0:
            pytest.fail(f"{geom}: no live trace rows — vacuous gate")
        if geom == "cyl_folded" and n_degenerate == 0:
            # This is the ONLY assertion in the tree that the forward is a
            # structural zero on the tangential (μ_r = 0) slots — measured
            # 2026-08-03 by mutating those rows with a LINEAR bug
            # (``out[tan] = ±ψ[tan]``): exactly 1 of 5076 tests reddened,
            # this one. The trace metric ``G = |Ω·n|·w_n`` is EXACTLY zero
            # there, so every G-weighted or solver-level gate is Mode-12
            # designed-green and structurally cannot see it.
            #
            # The ``if degenerate.size:`` branch above is therefore
            # load-bearing but self-silencing: swap this fixture to an
            # n_φ ≡ 0 (mod 4) folded rule (no tangential ordinate) and
            # the branch simply stops executing — green, with the
            # property unasserted anywhere. Fail loudly instead.
            pytest.fail(
                "cyl_folded carries no tangential ordinates — the tree's "
                "only catcher for the structural-zero trace row has gone "
                "vacuous. Restore an n_phi ≡ 2 (mod 4) folded rule here "
                "(see the #280 MANDATORY config on _full_space_states.mesh_cyl) rather than "
                "deleting this guard."
            )

    @pytest.mark.parametrize("geom", _GEOMS)
    def test_pure_outflow_rhs_round_trips(self, geom):
        """The previously-singular subspace: a rhs living ONLY on the
        outflow-trace rows must round-trip exactly (the old sweep
        mapped it to ZERO — the singular preconditioner's kernel)."""
        problem, lc, sweep = _lc_pair(geom)
        trace = problem.angular_trace
        boundary = AngularBoundaryFlux.zeros(problem.angular_trace)
        rng = np.random.default_rng(23)
        for face in boundary.layout.faces:
            out_rows = trace.outflow_indices_for_face(face)
            view = boundary.face_view(face)
            view[out_rows] = rng.normal(size=view[out_rows].shape)
        if isinstance(lc, CoupledOperator):
            rhs_a = _zero_source_composite(problem)
            for face in rhs_a.boundary.layout.faces:
                out_rows = trace.outflow_indices_for_face(face)
                rhs_a.boundary.face_view(face)[out_rows] = (
                    np.asarray(boundary.face_view(face))[out_rows]
                )
            rhs = CoupledField(systems=(
                rhs_a, RadialCharacteristicField.source_zeros(problem.radial_characteristic_field_space),
            ))
        else:
            rhs_a = FullField(
                interior=AngularFlux.zeros(problem.angular_bulk_space), boundary=boundary,
            )
            rhs = rhs_a
        psi = sweep.apply(rhs)
        psi_a = _system_a(psi)
        # Zero bulk source + zero inflow ⟹ the marched interior is zero;
        # the outflow trace carries the rhs's defect content (as −rhs_out
        # under the forward's streamed − ψ_out row convention).  On the
        # coupled arm the ψ½ rhs is zero too, so the marched System-B
        # state is zero and the seed feeds nothing into the bulk.
        np.testing.assert_allclose(
            np.asarray(psi_a.interior.values), 0.0, atol=1e-14,
            err_msg="a pure outflow-row rhs drives no interior flux",
        )
        if isinstance(psi, CoupledField):
            np.testing.assert_allclose(
                np.asarray(psi.systems[1].to_flat()), 0.0, atol=1e-14,
                err_msg="a pure outflow-row rhs drives no ψ½ state",
            )
        norm_out = float(np.abs(np.asarray(psi_a.boundary.values)).max())
        if not norm_out > 0.1:
            pytest.fail(
                f"the sweep must carry the outflow-row rhs through "
                f"(max |trace| = {norm_out:.2e}) — the ERR-071 dropped "
                f"term has regressed and the Krylov preconditioner is "
                f"singular again"
            )
        back = lc.apply(psi)
        back_a = _system_a(back)
        np.testing.assert_allclose(
            np.asarray(back_a.boundary.values),
            np.asarray(rhs_a.boundary.values),
            rtol=1e-12, atol=1e-13,
            err_msg="the outflow defect rows must round-trip",
        )

    def test_sentinel_scanner_tooth(self, monkeypatch):
        """The identity gate's own tooth: re-introduce the drop (an
        empty outflow selector reproduces the pre-fix clobber) and the
        pure-outflow round-trip must break."""
        from orpheus.numerics.spaces.angular_trace_space import (
            AngularTraceSpace,
        )

        monkeypatch.setattr(
            AngularTraceSpace,
            "outflow_indices_for_face",
            lambda self, face: np.array([], dtype=int),
        )
        problem, _lc, sweep = _lc_pair("slab_vacuum")
        boundary = AngularBoundaryFlux.zeros(problem.angular_trace)
        # populate what WOULD be the outflow rows (computed from the
        # quadrature directly — the monkeypatched selector is the
        # production read under mutation)
        mu = np.asarray(problem.quad.mu_x)
        rows = {"xmin": np.flatnonzero(mu < 0), "xmax": np.flatnonzero(mu > 0)}
        for face, out_rows in rows.items():
            boundary.face_view(face)[out_rows] = 1.0
        rhs = FullField(
            interior=AngularFlux.zeros(problem.angular_bulk_space), boundary=boundary,
        )
        psi = sweep.apply(rhs)
        norm_out = float(np.abs(np.asarray(psi.boundary.values)).max())
        if not norm_out < 1e-12:
            pytest.fail(
                f"the mutation (dropped outflow restore) must reproduce "
                f"the pre-fix annihilation (got max |trace| = "
                f"{norm_out:.2e}) — the gate is not sensing the fix"
            )
