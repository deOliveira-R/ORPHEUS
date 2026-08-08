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
residual and the end-of-solve certificate refused the claim — the
catch that exposed the class). The P1-DSA (d₁) Krylov posture (#2)
excited it deterministically.

This gate pins the identity the exact-inverse pair owes: round-trip a
RANDOM composite (interior + boundary, all rows live), and the pure
outflow-row leg that was the singular subspace.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.common.xs_library import make_mixture
from orpheus.geometry import BC, CoordSystem
from orpheus.geometry.mesh import Mesh1D
from orpheus.numerics.coupled_system import CoupledField, CoupledOperator
from orpheus.numerics.quadrature import Quadrature
from orpheus.transport.radial_characteristic_field import (
    RadialCharacteristicField,
)
from orpheus.sn.coupled_system import build_within_group_system
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.transport.fields.angular_boundary_flux import AngularBoundaryFlux
from orpheus.transport.fields.angular_flux import AngularFlux
from orpheus.transport.full_field import FullField

pytestmark = [
    pytest.mark.foundation,
    pytest.mark.catches("ERR-071"),
    # ERR-078: the ψ½ march's solve dropped the outflow-row rhs — the
    # System-B twin of ERR-071, caught by this file's coupled rows.
    pytest.mark.catches("ERR-078"),
]


def _mixtures():
    mix_a = make_mixture(
        sig_t=np.array([1.0]), sig_c=np.array([0.1]),
        sig_f=np.array([0.0]), nu=np.array([0.0]),
        chi=np.array([0.0]), sig_s=np.array([[0.9]]),
    )
    mix_b = make_mixture(
        sig_t=np.array([2.0]), sig_c=np.array([1.0]),
        sig_f=np.array([0.0]), nu=np.array([0.0]),
        chi=np.array([0.0]), sig_s=np.array([[1.0]]),
    )
    return {0: mix_a, 1: mix_b}


def _mesh_slab(bc_left: str) -> SNMesh:
    return SNMesh(
        Mesh1D(
            edges=np.array([0.0, 0.5, 1.5, 3.0, 5.0, 6.0, 8.0]),
            mat_ids=np.array([0, 1, 1, 0, 1, 0]),
            bc_left=BC(bc_left),
            bc_right=BC("vacuum"),
        ),
        Quadrature.gauss_legendre(n_ordinates=4),
        _mixtures(),
    )


def _mesh_cyl() -> SNMesh:
    # The #280 MANDATORY cylinder config, re-posed at the 6.3 flip onto
    # the admitted family: ``folded_product(4, 6)`` — the staggered
    # parent at n_φ ≡ 2 (mod 4) places φ = π/2 exactly, and the
    # roots-of-unity circle (E3) makes those ordinates' μ_r = 0.0
    # BIT-EXACT — so the rule carries degenerate pure-azimuthal
    # ordinates AND (like every admitted cylinder rule) a live ψ½
    # System B.  xmax-only trace layout — exercises the restore's
    # per-face membership loop on the curvilinear face set.
    return SNMesh(
        Mesh1D(
            edges=np.array([0.0, 0.3, 0.8, 1.0]),
            mat_ids=np.array([0, 1, 0]),
            bc_left=BC("reflective"),
            bc_right=BC("vacuum"),
            coord=CoordSystem.CYLINDRICAL,
        ),
        Quadrature.folded_product(n_mu=4, n_phi=6),
        _mixtures(),
    )


def _mesh_sphere() -> SNMesh:
    # The carrying SPHERE row — GL-4, one seed level; the same coupled
    # round-trip as the cylinder row.  The sphere's seed-carrying
    # inverse reciprocity was explicitly deferred pre-6.3 (the "#29
    # domain" note in test_loss_transpose_solve.G3); ERR-078's fix
    # covers both curvilinear arms, so both are gated here.
    return SNMesh(
        Mesh1D(
            edges=np.array([0.0, 0.3, 0.8, 1.0]),
            mat_ids=np.array([0, 1, 0]),
            bc_left=BC("reflective"),
            bc_right=BC("vacuum"),
            coord=CoordSystem.SPHERICAL,
        ),
        Quadrature.gauss_legendre(n_ordinates=4),
        _mixtures(),
    )


_MESHES = {
    # vacuum walls — every trace row live: inflow identities, outflow defects
    "slab_vacuum": lambda: _mesh_slab("vacuum"),
    # the identity is bc-INDEPENDENT: B is a coupling GAIN, never inside
    # the bare (L+C) — a reflective wall must not change the round-trip
    "slab_reflective": lambda: _mesh_slab("reflective"),
    "cyl_folded": _mesh_cyl,
    "sphere_gl": _mesh_sphere,
}

_GEOMS = list(_MESHES)


def _lc_pair(geom: str):
    """The production within-group forward + its exact inverse on a
    small heterogeneous mesh (every trace row is live: inflow rows are
    identities, outflow rows are defects).  Slabs carry the bare
    ``(L+C)``; the cylinder — carrying since the 6.3 flip — carries the
    upper-triangular coupled ``[[LC, Seeding], [None, march]]``, whose
    ``inverse()`` is the block back-substitution."""
    sn_mesh = _MESHES[geom]()
    system = build_within_group_system(
        sn_mesh, sn_mesh.material_xs_field(),
    )
    lc = system.implicit_operator
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
    return sn_mesh, lc, lc.inverse()


def _zero_source_composite(sn_mesh: SNMesh) -> FullField:
    """A zero SOURCE-role System-A carrier for the coupled arm's rhs.

    Role-honest member algebra: a solve's rhs is a SOURCE, and the
    substitution computes ``q_A − Seeding·ψ_B`` — cross-role arithmetic
    is forbidden by the typed fields, so a flux-role rhs raises at the
    block boundary (#289-F2)."""
    from orpheus.transport.source_sinks import (
        AngularBoundarySourceSink,
        AngularSourceSink,
    )

    return FullField(
        interior=AngularSourceSink.from_mesh(
            np.zeros((sn_mesh.quad.N, sn_mesh.ng, *sn_mesh.spatial_shape)),
            sn_mesh,
        ),
        boundary=AngularBoundarySourceSink.zeros_on(sn_mesh),
    )


def _random_state(sn_mesh: SNMesh, lc, seed: int):
    """A random rhs in ``lc``'s domain — the bare composite, or the
    coupled (bulk ⊕ ψ½) SOURCE-role state with EVERY member block
    populated (randomized through the coupled ``from_flat``, so the
    ψ½ member's interior ⊕ boundary blocks are live too)."""
    if not isinstance(lc, CoupledOperator):
        return _random_composite(sn_mesh, seed)
    template = CoupledField(systems=(
        _zero_source_composite(sn_mesh),
        RadialCharacteristicField.source_zeros_on(sn_mesh),
    ))
    flat = np.asarray(template.to_flat())
    rng = np.random.default_rng(seed + 1)
    return CoupledField.from_flat(rng.normal(size=flat.shape), template)


def _system_a(x):
    """Project the bulk (System-A) member of a possibly-coupled field."""
    return x.systems[0] if isinstance(x, CoupledField) else x


def _random_composite(sn_mesh: SNMesh, seed: int) -> FullField:
    """Every block populated — bulk, inflow-trace, AND the outflow-trace
    rows the old sweep dropped — with shapes read off the mesh (so the
    same builder serves slab and the xmax-only curvilinear layout)."""
    rng = np.random.default_rng(seed)
    interior = AngularFlux.from_mesh(
        rng.normal(size=(sn_mesh.quad.N, sn_mesh.ng, *sn_mesh.spatial_shape)),
        sn_mesh,
    )
    boundary = AngularBoundaryFlux.zeros_on(sn_mesh)
    for face in boundary.layout.faces:
        view = boundary.face_view(face)
        view[...] = rng.normal(size=view.shape)
    return FullField(interior=interior, boundary=boundary)


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
        sn_mesh, lc, sweep = _lc_pair(geom)
        trace = sn_mesh.angular_trace
        rhs = _random_state(sn_mesh, lc, seed=17)
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
                np.arange(sn_mesh.quad.N), live,
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
                "(see the #280 MANDATORY config on _mesh_cyl) rather than "
                "deleting this guard."
            )

    @pytest.mark.parametrize("geom", _GEOMS)
    def test_pure_outflow_rhs_round_trips(self, geom):
        """The previously-singular subspace: a rhs living ONLY on the
        outflow-trace rows must round-trip exactly (the old sweep
        mapped it to ZERO — the singular preconditioner's kernel)."""
        sn_mesh, lc, sweep = _lc_pair(geom)
        trace = sn_mesh.angular_trace
        boundary = AngularBoundaryFlux.zeros_on(sn_mesh)
        rng = np.random.default_rng(23)
        for face in boundary.layout.faces:
            out_rows = trace.outflow_indices_for_face(face)
            view = boundary.face_view(face)
            view[out_rows] = rng.normal(size=view[out_rows].shape)
        if isinstance(lc, CoupledOperator):
            rhs_a = _zero_source_composite(sn_mesh)
            for face in rhs_a.boundary.layout.faces:
                out_rows = trace.outflow_indices_for_face(face)
                rhs_a.boundary.face_view(face)[out_rows] = (
                    np.asarray(boundary.face_view(face))[out_rows]
                )
            rhs = CoupledField(systems=(
                rhs_a, RadialCharacteristicField.source_zeros_on(sn_mesh),
            ))
        else:
            rhs_a = FullField(
                interior=AngularFlux.zeros_on(sn_mesh), boundary=boundary,
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
        sn_mesh, _lc, sweep = _lc_pair("slab_vacuum")
        boundary = AngularBoundaryFlux.zeros_on(sn_mesh)
        # populate what WOULD be the outflow rows (computed from the
        # quadrature directly — the monkeypatched selector is the
        # production read under mutation)
        mu = np.asarray(sn_mesh.quad.mu_x)
        rows = {"xmin": np.flatnonzero(mu < 0), "xmax": np.flatnonzero(mu > 0)}
        for face, out_rows in rows.items():
            boundary.face_view(face)[out_rows] = 1.0
        rhs = FullField(
            interior=AngularFlux.zeros_on(sn_mesh), boundary=boundary,
        )
        psi = sweep.apply(rhs)
        norm_out = float(np.abs(np.asarray(psi.boundary.values)).max())
        if not norm_out < 1e-12:
            pytest.fail(
                f"the mutation (dropped outflow restore) must reproduce "
                f"the pre-fix annihilation (got max |trace| = "
                f"{norm_out:.2e}) — the gate is not sensing the fix"
            )
