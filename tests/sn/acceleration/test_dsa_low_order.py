r"""The 3b object gates — the production low-order build (D7/D8 + ties).

Three foundation clusters:

1. **The production tie**: :class:`DSALowOrderSystem` (the SN-side
   production build) equals the derivation-side reference builder
   ``orpheus.derivations.discrete.sn.dsa.build_consistent_dd_system``
   entry-for-entry, per group, on the heterogeneous non-uniform slab —
   both BC variants. The derivation is the algebra of record (its rows
   are theorems); a drift in either spelling reds here instead of
   forking silently.
2. **Admission teeth**: the loud seams (non-DD scheme, unsupported
   walls, the Σw = 2 convention boundary, the D-positivity guard)
   actually refuse.
3. **D7/D8 — the R/P object laws**: R conservation (hand-posed
   ``⟨1, R r⟩ = ⟨1, r⟩`` with explicit :math:`w_n, V_i`), the R
   identifications (``integrate_angular`` ≡ the ℓ=0 frame row ≡ the
   displacement tangent map, 0-ULP), and the exact round-trip
   ``R ∘ P = I`` with the frame's :math:`Y_0^0 = 1` table backing the
   normalized-injection spelling of P.

Levels: foundation (object identities; no physics reference).
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.common.xs_library import get_mixture
from orpheus.derivations.discrete.sn import dsa as dsa_reference
from orpheus.geometry import BC
from orpheus.geometry.mesh import Mesh1D
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.acceleration import DSACorrection, DSALowOrderSystem
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.transport.fields.angular_flux import AngularFlux

pytestmark = pytest.mark.foundation


def _slab(bc_left: str = "vacuum", bc_right: str = "vacuum") -> SNMesh:
    """Heterogeneous, non-uniform 4-cell slab, S4, 2 groups (the 3a tie
    fixture — mixtures carry real P1 data, so the (23c) D is exercised
    beyond the bare-P0 coincidence)."""
    mesh1d = Mesh1D(
        edges=np.array([0.0, 0.5, 1.5, 3.0, 5.0]),
        mat_ids=np.array([0, 1, 1, 0]),
        bc_left=BC(bc_left),
        bc_right=BC(bc_right),
    )
    quad = Quadrature.gauss_legendre(n_ordinates=4)
    return SNMesh(
        mesh1d, quad, {0: get_mixture("A", "2g"), 1: get_mixture("B", "2g")}
    )


def _reference_inputs(sn_mesh: SNMesh):
    """The reference builder's per-group inputs, gathered the same way
    the production build gathers them (data path shared; the FORMULAS
    are what the tie discriminates)."""
    h = np.diff(np.asarray(getattr(sn_mesh.mesh, "edges"), float))
    xs = sn_mesh.material_xs_field()
    sigma_t = np.asarray(xs.total_cross_section_field.values, float)
    ng = sigma_t.shape[0]
    mat_ids = np.asarray(sn_mesh.mat_map, int).ravel()
    fold = xs.foldable_sigma()
    sigma_s0 = np.stack([fold[int(m)] for m in mat_ids], axis=1)
    residual = xs.residual_sig_s()
    s1 = {
        mid: (
            np.asarray(np.diag(mats[1]), float)
            if len(mats) > 1
            else np.zeros(ng)
        )
        for mid, mats in residual.items()
    }
    sigma_s1 = np.stack([s1[int(m)] for m in mat_ids], axis=1)
    mu = np.asarray(sn_mesh.quad.mu_x, float)
    w = np.asarray(sn_mesh.quad.weights, float)
    return h, sigma_t, sigma_s0, sigma_s1, mu, w


class TestProductionTie:
    """Weld: the production build ≡ the derivation reference builder."""

    @pytest.mark.parametrize(
        "bc", [("vacuum", "vacuum"), ("reflective", "vacuum")]
    )
    @pytest.mark.verifies("sn-dsa-consistent-low-order")
    def test_low_order_matches_reference_builder(self, bc):
        """At ``scattering_order=1`` so the mixtures' real P1 rows
        exercise the (23c) transport-corrected D, not the bare-P0
        coincidence."""
        sn_mesh = _slab(*bc)
        system = DSALowOrderSystem.from_sn_mesh(sn_mesh, scattering_order=1)
        h, sigma_t, sigma_s0, sigma_s1, mu, w = _reference_inputs(sn_mesh)
        for g in range(sigma_t.shape[0]):
            a_ref, g_ref = dsa_reference.build_consistent_dd_system(
                h, sigma_t[g], sigma_s0[g], sigma_s1[g], mu, w, bc=bc
            )
            np.testing.assert_allclose(
                system.a_low[g], a_ref, rtol=0, atol=1e-15,
                err_msg=f"A_low group {g} must equal the proven reference",
            )
            np.testing.assert_allclose(
                system.g_map[g], g_ref, rtol=0, atol=1e-15,
                err_msg=f"G group {g} must equal the proven reference",
            )

    def test_solve_correction_realizes_the_reference_solve(self):
        """f0 = A⁻¹(G d) and the (28a) cell average, against a direct
        dense solve of the reference system."""
        sn_mesh = _slab()
        system = DSALowOrderSystem.from_sn_mesh(sn_mesh, scattering_order=1)
        h, sigma_t, sigma_s0, sigma_s1, mu, w = _reference_inputs(sn_mesh)
        rng = np.random.default_rng(7)
        d0 = rng.standard_normal((sigma_t.shape[0], h.shape[0]))
        f0 = system.solve_correction(d0)
        for g in range(sigma_t.shape[0]):
            a_ref, g_ref = dsa_reference.build_consistent_dd_system(
                h, sigma_t[g], sigma_s0[g], sigma_s1[g], mu, w
            )
            d = np.concatenate([d0[g], np.zeros_like(d0[g])])
            f_ref = np.linalg.solve(a_ref, g_ref @ d)
            np.testing.assert_allclose(f0[g], f_ref, rtol=1e-12, atol=1e-14)
            np.testing.assert_allclose(
                system.cell_update(f0)[g],
                0.5 * (f_ref[:-1] + f_ref[1:]),
                rtol=1e-12, atol=1e-14,
            )


class TestAdmissionTeeth:
    """The loud seams refuse — each guard actually bites."""

    def test_unsupported_boundary_refused(self):
        """The albedo/white seam guard. SNMesh's own registry pre-refuses
        such kinds today, so the guard is defense-in-depth for a future
        registry entry — exercised through a structural stub carrying
        the admission surface (geometry, scheme, bc kinds)."""
        from types import SimpleNamespace

        real = _slab()
        stub = SimpleNamespace(
            is_cartesian=True,
            ndim=1,
            mesh=real.mesh,
            scheme=SimpleNamespace(key="diamond_difference"),
            bc={
                "xmin": SimpleNamespace(kind="white"),
                "xmax": SimpleNamespace(kind="vacuum"),
            },
        )
        with pytest.raises(NotImplementedError, match="Marshak-albedo"):
            DSALowOrderSystem.from_sn_mesh(stub)  # type: ignore[arg-type]

    def test_non_dd_scheme_refused(self):
        from orpheus.transport.spatial.linear_discontinuous import (
            LinearDiscontinuous,
        )

        mesh1d = Mesh1D(
            edges=np.array([0.0, 0.5, 1.5, 3.0, 5.0]),
            mat_ids=np.array([0, 1, 1, 0]),
            bc_left=BC("vacuum"),
            bc_right=BC("vacuum"),
        )
        sn_mesh = SNMesh(
            mesh1d,
            Quadrature.gauss_legendre(n_ordinates=4),
            {0: get_mixture("A", "2g"), 1: get_mixture("B", "2g")},
            scheme=LinearDiscontinuous(),
        )
        with pytest.raises(NotImplementedError, match="diamond"):
            DSALowOrderSystem.from_sn_mesh(sn_mesh)

    def test_quadrature_convention_guard_fires(self):
        sn_mesh = _slab()
        h, sigma_t, sigma_s0, sigma_s1, mu, w = _reference_inputs(sn_mesh)
        with pytest.raises(ValueError, match="Σw = 2"):
            DSALowOrderSystem._build(
                h, sigma_t, sigma_s0, sigma_s1, mu, w / 2.0,
                ("vacuum", "vacuum"),
            )

    def test_d_positivity_guard_fires(self):
        sn_mesh = _slab()
        h, sigma_t, sigma_s0, _sigma_s1, mu, w = _reference_inputs(sn_mesh)
        bad_s1 = np.full_like(sigma_t, 2.0) * sigma_t  # σ_s1 > σ_t
        with pytest.raises(ValueError, match="positive"):
            DSALowOrderSystem._build(
                h, sigma_t, sigma_s0, bad_s1, mu, w, ("vacuum", "vacuum"),
            )


class TestRestrictionProlongation:
    """D7/D8 — the R/P object laws on the DSA carrier."""

    @pytest.fixture()
    def psi(self):
        sn_mesh = _slab()
        rng = np.random.default_rng(11)
        values = rng.standard_normal((4, 2, 4))
        return AngularFlux.from_mesh(values, sn_mesh)

    @pytest.mark.verifies("sn-dsa-restriction")
    def test_d7_restriction_conserves_particles(self, psi):
        r"""⟨1, R r⟩ = ⟨1, r⟩ — hand-posed with explicit w_n and V_i
        (structurally independent of the einsum body)."""
        sn_mesh = psi.mesh
        w = np.asarray(sn_mesh.quad.weights, float)
        volumes = np.diff(np.asarray(sn_mesh.mesh.edges, float))
        reduced = psi.integrate_angular().values  # (ng, nx)
        lhs = 0.0
        rhs = 0.0
        for g in range(reduced.shape[0]):
            for i in range(reduced.shape[1]):
                lhs += volumes[i] * reduced[g, i]
                for n in range(w.shape[0]):
                    rhs += volumes[i] * w[n] * psi.values[n, g, i]
        np.testing.assert_allclose(lhs, rhs, rtol=1e-14, atol=0)

    @pytest.mark.verifies("sn-dsa-restriction")
    def test_d8_restriction_is_the_frame_moment_row(self, psi):
        r"""``integrate_angular`` ≡ the ℓ=0 analysis row of
        ``Quadrature.angular_frame(0)`` (Y⁰₀ = 1 under the no-prefactor
        SH convention ⟹ the row IS the weight vector) — 0-ULP."""
        sn_mesh = psi.mesh
        frame = sn_mesh.quad.angular_frame(0)
        table = np.asarray(frame.table)  # (N, 1, 1); Y00 ≡ 1
        np.testing.assert_array_equal(table.ravel(), np.ones(4))
        w = np.asarray(sn_mesh.quad.weights, float)
        frame_row = np.einsum("n,ng...->g...", w * table.ravel(), psi.values)
        np.testing.assert_array_equal(
            psi.integrate_angular().values, frame_row,
        )

    def test_d8_displacement_reduction_is_the_tangent_map(self, psi):
        r"""``AngularDisplacement.integrate_angular`` carries the SAME
        body (a linear map is its own tangent map) — 0-ULP, and typed
        as a ScalarDisplacement."""
        from orpheus.transport.displacements.scalar_displacement import (
            ScalarDisplacement,
        )

        zero = AngularFlux.from_mesh(np.zeros_like(psi.values), psi.mesh)
        delta = psi - zero  # AngularDisplacement with psi's values
        reduced = delta.integrate_angular()
        if not isinstance(reduced, ScalarDisplacement):
            pytest.fail(
                f"the displacement reduction must stay in the "
                f"displacement role; got {type(reduced).__name__}"
            )
        np.testing.assert_array_equal(
            reduced.values, psi.integrate_angular().values,
        )

    def test_d8_prolongation_round_trip_is_identity(self, psi):
        r"""R ∘ P = I exactly: the normalized isotropic injection's
        moment-0 reproduces the input scalar values (Σ w · (x/Σw) = x
        up to one product-sum; pinned at 1-ULP-scale rtol)."""
        sn_mesh = psi.mesh
        corr = DSACorrection.from_sn_mesh(sn_mesh)
        phi = psi.integrate_angular().values
        sum_w = float(np.asarray(sn_mesh.quad.weights).sum())
        injected = AngularFlux.from_mesh(
            np.broadcast_to(phi[None] / sum_w, psi.values.shape).copy(),
            sn_mesh,
        )
        np.testing.assert_allclose(
            injected.integrate_angular().values, phi, rtol=1e-15, atol=0,
        )
        # and the operator's own injection normalization agrees
        np.testing.assert_array_equal(corr._sum_w, sum_w)
