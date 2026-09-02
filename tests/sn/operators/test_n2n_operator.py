r"""Gates for :class:`~orpheus.transport.operators.n2n.N2NOperator` — the
first-class :math:`(n,2n)` source operator (CS4c step 3, design record
§14.1).

The extraction's claims, each with its own witness:

* the composite action IS the isotropic lift —
  ``N2N.apply(ψ) ≡ (1/W)·frame.conjugate(N2NMomentOperator).apply(ψ)``
  (the ℓ=0 conjugation the fast path realizes without a moment tensor);
* the Euclidean transpose closes reciprocity
  ``⟨N2N ψ, χ⟩ = ⟨ψ, N2Nᵀ χ⟩`` on the full per-ordinate contraction;
* the carrier arms mirror S's (composite bulk-only; moment arm ≡ angular
  arm on ``φ = Mψ``; ScalarFlux refused toward the ENERGY binding);
* the EXTRACTION preserves the value: ``S.apply + N2N.apply`` reproduces
  the pre-extraction fused source (pinned by the pre-T3 snapshots in
  ``test_scattering_operator.py::TestAnisoMomentSourcePath``, which now
  compare the composed pair against the frozen artifact).

Fixture: asymmetric ``Σ_s0`` + non-zero, cross-group ``Σ_2n`` (L27 — a
group flip or a dropped channel is detectable).
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry import BC, CoordSystem, Mesh1D
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.sn.solver import SNSolver
from orpheus.transport.fields.angular_flux import AngularFlux
from orpheus.transport.fields.harmonic_moment_flux import HarmonicMomentFlux
from orpheus.transport.fields.scalar_flux import ScalarFlux
from orpheus.transport.full_field import FullField
from orpheus.transport.operators.n2n import N2NOperator
from orpheus.transport.operators.scattering import N2NMomentOperator
from orpheus.transport.source_sinks import AngularSourceSink

from tests.sn._test_helpers import material_xs_from_raw
from orpheus.numerics.basis.spherical_harmonic_basis import SphericalHarmonicBasis

pytestmark = pytest.mark.foundation

_SIGS0 = np.array([[0.38, 0.10], [0.05, 0.90]])
_SIG2 = np.array([[0.00, 0.03], [0.01, 0.00]])
_NX, _NG = 4, 2


def _solver():
    mat_xs = material_xs_from_raw(
        sig_s={0: [_SIGS0]}, sig2={0: _SIG2},
        cells_by_mat={0: (np.arange(_NX), np.zeros(_NX, dtype=int))},
        ng=_NG, nx=_NX, ny=1,
    )
    mesh = Mesh1D(
        edges=np.linspace(0.0, 1.0, _NX + 1),
        mat_ids=np.zeros(_NX, dtype=int),
        coord=CoordSystem.CARTESIAN,
        bc_left=BC("vacuum"),
        bc_right=BC("vacuum"),
    )
    sn = SNMesh(mesh, Quadrature.gauss_legendre(n_ordinates=4), mat_xs.materials)
    # SNMesh re-derives its own mat_xs field; use the SOLVER path so the
    # operator pair is the production mint (injection-consistent).
    return SNSolver(sn), sn


def _psi(sn, seed=3):
    rng = np.random.default_rng(seed)
    values = rng.uniform(0.05, 1.0, size=sn.angular_bulk_space.shape)
    return AngularFlux(values=values, space=sn.angular_bulk_space)


class TestLiftIsTheConjugation:
    def test_apply_equals_l0_conjugation(self):
        r"""``N2N.apply(ψ)`` ≡ ``(1/W)·conjugate(N2NMoment).apply(ψ)`` —
        the fast path IS the ℓ=0 frame conjugation (algebra eager,
        performance lazy; the c₀ = 1 theorem the P0 crosscheck pins)."""
        solver, sn = _solver()
        n2n = solver.n2n_op
        psi = _psi(sn)
        got = n2n.apply(psi).values

        S = solver.scattering_op
        # ⛔ RE-KEYED 2026-09-02 (#429). The moment operator was minted from a
        # hand-named ``SphericalHarmonicBasis``; on this 1-D fixture the
        # frame now binds the FLAT Legendre basis, so the conjugation's ends
        # no longer composed (`[M]` ``A.domain=SphericalHarmonicSpace(1,1)``
        # vs ``B.codomain=LegendreSpace(1,)``). ⭐ Reading the basis OFF THE
        # FRAME is not merely the repair — it is the property this gate is
        # about: the conjugation ``R ∘ X ∘ M`` is well-posed exactly when
        # ``X``'s ends are the frame's own coefficient space, and the
        # composition guard says so.
        moment = N2NMomentOperator.from_material_xs(
            mat_xs=solver.mat_xs, basis=S.frame.basis,
        )
        assert S.frame.basis.L == S.scattering_order
        conjugated = S.frame.conjugate(moment).apply(psi.values)
        np.testing.assert_allclose(
            got, np.asarray(conjugated) / n2n.total_weight,
            rtol=1e-13, atol=1e-16,
            err_msg="the (n,2n) lift drifted from its own ℓ=0 conjugation",
        )

    def test_euclidean_reciprocity(self):
        r"""``⟨N2N ψ, χ⟩ = ⟨ψ, N2Nᵀ χ⟩`` (full per-ordinate contraction)."""
        solver, sn = _solver()
        n2n = solver.n2n_op
        psi, chi = _psi(sn, 5), _psi(sn, 6)
        lhs = float(np.sum(n2n.apply(psi).values * chi.values))
        rhs = float(np.sum(psi.values * n2n.apply_transpose(chi.values)))
        np.testing.assert_allclose(lhs, rhs, rtol=1e-12)

    def test_transpose_reds_on_group_flip(self):
        r"""Negative leg (vv #11): a hand-flipped transpose breaks
        reciprocity on the asymmetric fixture — the identity has teeth."""
        solver, sn = _solver()
        n2n = solver.n2n_op
        psi, chi = _psi(sn, 5), _psi(sn, 6)
        lhs = float(np.sum(n2n.apply(psi).values * chi.values))
        # The WRONG transpose: forward applied to χ (un-transposed K).
        wrong = float(np.sum(psi.values * n2n.apply(chi).values))
        if abs(lhs - wrong) <= 1e-10 * abs(lhs):
            pytest.fail(
                "group-flip control did not move reciprocity — the "
                "fixture cannot see a transposed kernel",
            )


class TestCarrierArms:
    def test_composite_arm_is_bulk_only(self):
        solver, sn = _solver()
        state = FullField(
            interior=_psi(sn, 9),
            boundary=__import__(
                "orpheus.transport.fields.angular_boundary_flux",
                fromlist=["AngularBoundaryFlux"],
            ).AngularBoundaryFlux.zeros(sn.angular_trace),
        )
        out = solver.n2n_op.apply(state)
        assert isinstance(out, FullField)
        np.testing.assert_array_equal(out.boundary.values, 0.0)
        np.testing.assert_array_equal(
            out.interior.values, solver.n2n_op.apply(state.interior).values,
        )

    def test_moment_arm_equals_angular_arm_on_projected_flux(self):
        r"""``N2N.apply(Mψ) ≡ N2N.apply(ψ)`` — the ℓ=0 moment IS the
        scalar flux, so windowing loses nothing (mirrors S's guard)."""
        solver, sn = _solver()
        S = solver.scattering_op
        psi = _psi(sn, 11)
        moments = S.frame.analysis.apply(psi.values)
        phi_field = HarmonicMomentFlux.from_mesh_and_L(
            moments, sn, S.scattering_order,
        )
        np.testing.assert_allclose(
            solver.n2n_op.apply(phi_field).values,
            solver.n2n_op.apply(psi).values,
            rtol=1e-13, atol=1e-16,
        )

    def test_scalar_flux_refused_toward_the_energy_binding(self):
        solver, sn = _solver()
        phi = ScalarFlux(
            values=np.ones(sn.bulk_space.shape),
            space=sn.bulk_space,
        )
        with pytest.raises(TypeError, match="ENERGY binding"):
            solver.n2n_op.apply(phi)

    def test_unknown_carrier_refused_with_named_operator(self):
        solver, _ = _solver()
        with pytest.raises(TypeError, match="N2NOperator.apply"):
            solver.n2n_op.apply(object())


class TestAdmission:
    def test_mismatched_frame_refused(self):
        """CS4c step-4 harmonization: the retained state is the L=0 FRAME
        (was a bare ``weights`` array — this row re-keyed with it); a
        frame minted from a DIFFERENT quadrature than the posed axis is
        refused at construction (the F guard, mirrored)."""
        solver, sn = _solver()
        from dataclasses import replace

        from orpheus.numerics.quadrature import Quadrature
        from orpheus.transport.frames.harmonic_frame import HarmonicFrame

        other = Quadrature.gauss_legendre(n_ordinates=2)
        wrong_frame = HarmonicFrame.from_galerkin(other.angular_frame(0))
        with pytest.raises(TypeError, match="mint the frame"):
            replace(solver.n2n_op, frame=wrong_frame)

    def test_from_solver_data_refuses_a_bare_space(self):
        solver, _ = _solver()
        from orpheus.numerics.space import FunctionSpace

        with pytest.raises(TypeError, match="composite FullFieldSpace"):
            N2NOperator.from_solver_data(
                mat_xs=solver.mat_xs,
                space=FunctionSpace("bare", (2, 4, 1)),
            )

    def test_total_weight_is_the_measure_mass(self):
        solver, sn = _solver()
        np.testing.assert_allclose(
            solver.n2n_op.total_weight, float(sn.quad.weights.sum()),
        )
