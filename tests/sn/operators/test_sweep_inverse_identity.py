r"""The composite sweep-inverse identity — ``(L+C) ∘ (L+C)⁻¹ ≡ I``
on the FULL composite space, outflow-trace rows included (ERR-071).

The forward's boundary block is a sibling of the bulk: the outflow-trace
row is the self-consistency DEFECT ``streamed − ψ_out`` (the sign is
the authoritative forward's, loss_representation, pinned by this gate's
round-trip) and the inflow
row the identity on the given inflow. The exact inverse must therefore
emit ``ψ_out = streamed − rhs_out`` (the forward's row is ``streamed −
ψ_out``; the sign is pinned by this gate's round-trip) — and until
2026-07-26 the sweep
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
from orpheus.numerics.coupled_system import CoupledOperator
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.coupled_system import build_within_group_system
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.transport.fields.angular_boundary_flux import AngularBoundaryFlux
from orpheus.transport.fields.angular_flux import AngularFlux
from orpheus.transport.full_field import FullField

pytestmark = [pytest.mark.foundation, pytest.mark.catches("ERR-071")]


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
    # The #280 MANDATORY cylinder config: a PRODUCT quadrature carries
    # degenerate pure-azimuthal ordinates AND a live seed-fold — the two
    # cylinder-specific terms level_symmetric nulls.  Non-carrying (no
    # ψ½ System B), xmax-only trace layout — exercises the restore's
    # per-face membership loop on the curvilinear face set.
    return SNMesh(
        Mesh1D(
            edges=np.array([0.0, 0.3, 0.8, 1.0]),
            mat_ids=np.array([0, 1, 0]),
            bc_left=BC("reflective"),
            bc_right=BC("vacuum"),
            coord=CoordSystem.CYLINDRICAL,
        ),
        Quadrature.product(n_mu=4, n_phi=8),
        _mixtures(),
    )


_MESHES = {
    # vacuum walls — every trace row live: inflow identities, outflow defects
    "slab_vacuum": lambda: _mesh_slab("vacuum"),
    # the identity is bc-INDEPENDENT: B is a coupling GAIN, never inside
    # the bare (L+C) — a reflective wall must not change the round-trip
    "slab_reflective": lambda: _mesh_slab("reflective"),
    "cyl_product": _mesh_cyl,
}

_GEOMS = list(_MESHES)


def _lc_pair(geom: str):
    """The production ``(L+C)`` forward + its sweep inverse on a small
    heterogeneous non-carrying mesh (every trace row is live: inflow
    rows are identities, outflow rows are defects)."""
    sn_mesh = _MESHES[geom]()
    system = build_within_group_system(
        sn_mesh, sn_mesh.material_xs_field(),
    )
    lc = system.resolvent
    if isinstance(lc, CoupledOperator):
        pytest.fail(
            f"{geom}: a non-carrying mesh must carry the bare (L+C) arm"
        )
    return sn_mesh, lc, lc.inverse()


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
    @pytest.mark.parametrize("geom", _GEOMS)
    def test_forward_of_inverse_is_identity_on_a_random_composite(
        self, geom,
    ):
        """``(L+C)((L+C)⁻¹ rhs) ≡ rhs`` with EVERY block populated —
        bulk, inflow-trace, and the outflow-trace rows the old sweep
        dropped.

        Honest scope on the trace: the identity is claimed on the
        inflow ∪ outflow rows.  A DEGENERATE pure-azimuthal ordinate
        (``μ_r = 0``, cylinder product quadrature — excluded from BOTH
        selectors) has NO streaming coupling to the face: its trace
        slot is a free DOF of the composite (#284), where the forward
        is a structural ZERO row and the inverse completes with the
        identity (seed passthrough).  Both halves of that pair are
        asserted explicitly; on a slab the degenerate set is empty and
        the claim is the full-trace identity."""
        sn_mesh, lc, sweep = _lc_pair(geom)
        trace = sn_mesh.angular_trace
        rhs = _random_composite(sn_mesh, seed=17)
        psi = sweep.apply(rhs)
        back = lc.apply(psi)
        np.testing.assert_allclose(
            np.asarray(back.interior.values),
            np.asarray(rhs.interior.values),
            rtol=1e-12, atol=1e-12,
            err_msg=f"{geom}: (L+C)∘(L+C)⁻¹ must be the identity on "
                    f"the bulk",
        )
        n_live = 0
        for face in rhs.boundary.layout.faces:
            live = np.union1d(
                trace.inflow_indices_for_face(face),
                trace.outflow_indices_for_face(face),
            )
            degenerate = np.setdiff1d(
                np.arange(sn_mesh.quad.N), live,
            )
            n_live += live.size
            np.testing.assert_allclose(
                np.asarray(back.boundary.face_view(face))[live],
                np.asarray(rhs.boundary.face_view(face))[live],
                rtol=1e-12, atol=1e-12,
                err_msg=f"{geom}/{face}: (L+C)∘(L+C)⁻¹ must be the "
                        f"identity on the live trace — inflow "
                        f"identities AND outflow defect rows",
            )
            if degenerate.size:
                np.testing.assert_allclose(
                    np.asarray(back.boundary.face_view(face))[degenerate],
                    0.0, atol=1e-12,
                    err_msg=f"{geom}/{face}: the forward must be a "
                            f"structural zero row on the degenerate "
                            f"(μ_r = 0) free-DOF trace slots",
                )
                np.testing.assert_allclose(
                    np.asarray(psi.boundary.face_view(face))[degenerate],
                    np.asarray(rhs.boundary.face_view(face))[degenerate],
                    rtol=1e-12, atol=1e-12,
                    err_msg=f"{geom}/{face}: the inverse must complete "
                            f"the free-DOF slots with the identity "
                            f"(seed passthrough)",
                )
        if not n_live > 0:
            pytest.fail(f"{geom}: no live trace rows — vacuous gate")

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
        rhs = FullField(
            interior=AngularFlux.zeros_on(sn_mesh), boundary=boundary,
        )
        psi = sweep.apply(rhs)
        # Zero bulk source + zero inflow ⟹ the marched interior is zero;
        # the outflow trace carries the rhs's defect content (as −rhs_out
        # under the forward's streamed − ψ_out row convention).
        np.testing.assert_allclose(
            np.asarray(psi.interior.values), 0.0, atol=1e-14,
            err_msg="a pure outflow-row rhs drives no interior flux",
        )
        norm_out = float(np.abs(np.asarray(psi.boundary.values)).max())
        if not norm_out > 0.1:
            pytest.fail(
                f"the sweep must carry the outflow-row rhs through "
                f"(max |trace| = {norm_out:.2e}) — the ERR-071 dropped "
                f"term has regressed and the Krylov preconditioner is "
                f"singular again"
            )
        back = lc.apply(psi)
        np.testing.assert_allclose(
            np.asarray(back.boundary.values),
            np.asarray(rhs.boundary.values),
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
