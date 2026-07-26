r"""``IntegratedReactionRate`` — :math:`\int_V \langle \Sigma_x,\phi\rangle\,dV`.

The volume-integrated reaction rate (the k-eigenvalue numerator/denominator,
and the :math:`\phi^\dagger\!=\!1` degenerate of the homogenization PG bilinear).
Three legs:

* **A — hand volume-weighted loop** on a **non-uniform-volume** mesh: the
  per-cell ``V`` weighting is only constrained when cells have distinct volumes
  (a uniform mesh nulls it — Mode-10). The reference is an explicit Python
  triple-loop ``Σ_cells V_cell Σ_g Σx_g φ_g`` (no numpy reduction shared with
  the SUT).
* **B — closed-form k∞ ratio**: on the homogeneous Mixture-A medium,
  ``IRR(νΣf)/IRR(Σa) == k∞`` (the volume measure CANCELS in the ratio — so this
  pins the fold is consistent between numerator and denominator). 4G mandatory
  (A 2G's spectrum is coincidentally flat).
* **C — (n,2n) activation in the SN production rate**: every ``xs_library``
  mixture has ``Sig2 = 0``, so a library fixture would null the (n,2n) term
  vacuously — this uses ``_make_mixture_with_n2n`` to genuinely activate it and
  assert it ADDS to the IRR fission production.

``foundation`` — L0 value correctness + software invariants.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.common.eigenvalue import kinf_and_spectrum_homogeneous
from orpheus.derivations.common.xs_library import get_mixture
from orpheus.geometry import BC, CoordSystem, Mesh1D, Mesh2D
from orpheus.numerics.functional import Functional
from orpheus.numerics.operator import LinearOperator
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.sn.solver import SNSolver
from orpheus.transport.reaction_rate_functional import IntegratedReactionRate

from tests.sn._test_helpers import placeholder_materials
from tests.transport._functional_helpers import cross_section_field, require, slab_mesh

pytestmark = pytest.mark.foundation


def _non_uniform_slab(ng: int) -> SNMesh:
    """A 1-D slab with NON-UNIFORM cell widths → distinct cell volumes.

    A uniform mesh makes the per-cell ``V`` a constant that the ratio/sum can't
    distinguish from 1 (Mode-10); distinct widths genuinely constrain the
    ``∫·dV`` weighting.
    """
    edges = np.array([0.0, 0.1, 0.3, 0.6, 1.0])  # widths 0.1 / 0.2 / 0.3 / 0.4
    mesh = Mesh1D(
        edges=edges, mat_ids=np.zeros(4, dtype=int), coord=CoordSystem.CARTESIAN,
        bc_left=BC("vacuum"), bc_right=BC("vacuum"),
    )
    return SNMesh(mesh, Quadrature.gauss_legendre(n_ordinates=4), placeholder_materials(ng=ng))


def _hand_volume_integral(sigx: np.ndarray, phi: np.ndarray, volumes: np.ndarray) -> float:
    """``Σ_cells V_cell Σ_g Σx_g φ_g`` by explicit Python loops (no numpy reduction)."""
    ng, ncells = sigx.shape[0], sigx.shape[1]
    vols = np.asarray(volumes).ravel()
    total = 0.0
    for c in range(ncells):
        cell = 0.0
        for g in range(ng):
            cell += float(sigx[g, c]) * float(phi[g, c])
        total += cell * float(vols[c])
    return total


class TestVolumeIntegralAgainstHandLoop:
    def test_evaluate_matches_hand_volume_weighted_loop_nonuniform(self):
        sn = _non_uniform_slab(ng=3)
        nx = sn.spatial_shape[0]
        rng = np.random.default_rng(2026)
        sigx = rng.uniform(0.1, 1.0, size=(3, nx))
        phi = rng.uniform(0.05, 1.0, size=(3, nx))
        got = IntegratedReactionRate(cross_section_field(sigx, sn)).evaluate(phi)
        ref = _hand_volume_integral(sigx, phi, sn.volumes)
        np.testing.assert_allclose(got, ref, rtol=1e-13)

    def test_mesh_has_distinct_cell_volumes(self):
        """The non-uniform fixture genuinely constrains ``V`` (Mode-10 guard)."""
        sn = _non_uniform_slab(ng=2)
        require(
            len(set(np.round(np.asarray(sn.volumes).ravel(), 12))) > 1,
            "the test mesh must have DISTINCT cell volumes — else ∫·dV is nulled.",
        )


class TestClosedFormKinfRatio:
    @pytest.mark.parametrize("ng_key, kinf", [("2g", 1.8750000000), ("4g", 1.4877619048)])
    def test_kinf_is_ratio_of_integrated_rates(self, ng_key, kinf):
        mix = get_mixture("A", ng_key)
        sigp = np.asarray(mix.SigP)
        siga = np.asarray(mix.absorption_xs)
        sigt = np.asarray(mix.SigT)
        sigs0 = np.asarray(mix.SigS[0].todense())
        chi = np.asarray(mix.chi)
        k, phi = kinf_and_spectrum_homogeneous(sigt, sigs0, sigp, chi)
        ng = len(sigp)
        sn = slab_mesh(nx=4, ng=ng)  # uniform OK — the volume CANCELS in the ratio
        nx = sn.spatial_shape[0]
        phi_field = np.repeat(phi[:, None], nx, axis=1)  # flat in space

        prod = IntegratedReactionRate(
            cross_section_field(np.repeat(sigp[:, None], nx, axis=1), sn)
        ).evaluate(phi_field)
        absn = IntegratedReactionRate(
            cross_section_field(np.repeat(siga[:, None], nx, axis=1), sn)
        ).evaluate(phi_field)

        require(prod > 0 and absn > 0, "integrated rates must be positive.")
        np.testing.assert_allclose(prod / absn, k, rtol=1e-10)
        np.testing.assert_allclose(k, kinf, rtol=1e-9)

    def test_category_membership(self):
        sn = slab_mesh(nx=2, ng=2)
        irr = IntegratedReactionRate(cross_section_field(np.array([[0.025, 0.025], [0.2, 0.2]]), sn))
        require(isinstance(irr, Functional), "IntegratedReactionRate must satisfy the Functional Protocol.")
        require(
            not isinstance(irr, LinearOperator),
            "a Functional must NOT satisfy LinearOperator (no apply/capabilities).",
        )


class TestAdjointWeightedBilinear:
    r"""The P6 (#281) bilinear arm: ``evaluate(phi, adjoint=φ*)`` = ∫⟨φ*⊙Σx, φ⟩dV.

    The vector-channel worth numerator of the eigenvalue-consistent collapse
    (theorem T1, :mod:`orpheus.derivations.common.homogenization`). Legs:
    the hand triple-loop reference on the non-uniform mesh (the ``∫·dV``
    weighting genuinely constrained), the EXACT flat-adjoint degenerate law,
    linearity in the adjoint slot (exact under a power-of-two scale), the
    ⟨φ*,Σφ⟩ = ⟨φ,Σφ*⟩ symmetry (diagonal Σ), and the mesh-identity guard.
    """

    @staticmethod
    def _hand_bilinear_integral(
        sigx: np.ndarray, phis: np.ndarray, phi: np.ndarray, volumes: np.ndarray,
    ) -> float:
        """``Σ_cells V_cell Σ_g φ*_g Σx_g φ_g`` by explicit Python loops."""
        ng, ncells = sigx.shape[0], sigx.shape[1]
        vols = np.asarray(volumes).ravel()
        total = 0.0
        for c in range(ncells):
            cell = 0.0
            for g in range(ng):
                cell += float(phis[g, c]) * float(sigx[g, c]) * float(phi[g, c])
            total += cell * float(vols[c])
        return total

    def test_bilinear_matches_hand_loop_nonuniform(self):
        from orpheus.transport.fields.scalar_flux import ScalarFlux

        sn = _non_uniform_slab(ng=3)
        nx = sn.spatial_shape[0]
        rng = np.random.default_rng(2027)
        sigx = rng.uniform(0.1, 1.0, size=(3, nx))
        phi = rng.uniform(0.05, 1.0, size=(3, nx))
        phis = rng.uniform(0.2, 2.0, size=(3, nx))
        got = IntegratedReactionRate(cross_section_field(sigx, sn)).evaluate(
            phi, adjoint=ScalarFlux.from_mesh(phis, sn),
        )
        ref = self._hand_bilinear_integral(sigx, phis, phi, sn.volumes)
        np.testing.assert_allclose(got, ref, rtol=1e-13)

    def test_flat_adjoint_degenerates_to_unweighted_exactly(self):
        """φ* = 1 reproduces the unweighted call BIT-FOR-BIT (1.0·x == x in
        IEEE floats; the contraction body is shared) — the φ†=1 degenerate
        law the class docstring claims."""
        from orpheus.transport.fields.scalar_flux import ScalarFlux

        sn = _non_uniform_slab(ng=2)
        nx = sn.spatial_shape[0]
        rng = np.random.default_rng(7)
        sigx = rng.uniform(0.1, 1.0, size=(2, nx))
        phi = rng.uniform(0.05, 1.0, size=(2, nx))
        irr = IntegratedReactionRate(cross_section_field(sigx, sn))
        ones = ScalarFlux.from_mesh(np.ones((2, nx)), sn)
        require(
            irr.evaluate(phi, adjoint=ones) == irr.evaluate(phi),
            "adjoint=1 must reproduce the unweighted integral EXACTLY "
            "(the φ†=1 degenerate law).",
        )

    def test_adjoint_slot_scales_linearly_exact(self):
        """Doubling φ* doubles the bilinear EXACTLY (scaling by 2.0 is
        round-off-free per IEEE, so linearity in the adjoint slot admits an
        exact assertion, not a tolerance)."""
        from orpheus.transport.fields.scalar_flux import ScalarFlux

        sn = _non_uniform_slab(ng=2)
        nx = sn.spatial_shape[0]
        rng = np.random.default_rng(11)
        sigx = rng.uniform(0.1, 1.0, size=(2, nx))
        phi = rng.uniform(0.05, 1.0, size=(2, nx))
        phis = rng.uniform(0.2, 2.0, size=(2, nx))
        irr = IntegratedReactionRate(cross_section_field(sigx, sn))
        base = irr.evaluate(phi, adjoint=ScalarFlux.from_mesh(phis, sn))
        doubled = irr.evaluate(phi, adjoint=ScalarFlux.from_mesh(2.0 * phis, sn))
        require(
            doubled == 2.0 * base,
            f"bilinear must be exactly linear in the adjoint slot under a "
            f"power-of-two scale: {doubled!r} != 2·{base!r}.",
        )

    def test_bilinear_symmetry_under_slot_swap(self):
        """⟨φ*, Σφ⟩ == ⟨φ, Σφ*⟩ — the diagonal-Σ symmetry of the bilinear
        (the two slots pair through the same multiplicative kernel)."""
        from orpheus.transport.fields.scalar_flux import ScalarFlux

        sn = _non_uniform_slab(ng=3)
        nx = sn.spatial_shape[0]
        rng = np.random.default_rng(13)
        sigx = rng.uniform(0.1, 1.0, size=(3, nx))
        a = rng.uniform(0.05, 1.0, size=(3, nx))
        b = rng.uniform(0.2, 2.0, size=(3, nx))
        irr = IntegratedReactionRate(cross_section_field(sigx, sn))
        ab = irr.evaluate(a, adjoint=ScalarFlux.from_mesh(b, sn))
        ba = irr.evaluate(b, adjoint=ScalarFlux.from_mesh(a, sn))
        np.testing.assert_allclose(ab, ba, rtol=1e-13)

    def test_foreign_mesh_adjoint_raises(self):
        """The σ↔geometry pairing tier: the adjoint must live on the SAME
        mesh OBJECT — an equal-shaped foreign mesh raises."""
        from orpheus.transport.fields.scalar_flux import ScalarFlux

        sn = _non_uniform_slab(ng=2)
        other = _non_uniform_slab(ng=2)          # equal-shaped, different object
        nx = sn.spatial_shape[0]
        rng = np.random.default_rng(17)
        irr = IntegratedReactionRate(
            cross_section_field(rng.uniform(0.1, 1.0, size=(2, nx)), sn)
        )
        foreign = ScalarFlux.from_mesh(rng.uniform(0.2, 2.0, size=(2, nx)), other)
        with pytest.raises(ValueError, match="different mesh"):
            irr.evaluate(rng.uniform(0.05, 1.0, size=(2, nx)), adjoint=foreign)


class TestN2NActivationInProductionRate:
    """SN ``compute_production_rate`` = IRR(νΣf) + (n,2n) — the (n,2n) genuinely adds."""

    def test_production_includes_n2n_when_sig2_nonzero(self):
        from tests.cp.test_verification import _make_mixture_with_n2n

        mat = _make_mixture_with_n2n("2g")  # Sig2 ≠ 0 (xs_library mixtures are all Sig2=0)
        mesh = Mesh2D(
            edges_x=np.linspace(0.0, 0.8, 5), edges_y=np.linspace(0.0, 0.6, 4),
            mat_map=np.zeros((4, 3), dtype=int),
        )
        solver = SNSolver(SNMesh(mesh, Quadrature.lebedev(order=17), {0: mat}))
        ng = solver.ng
        nx, ny = solver.sn_mesh.spatial_shape
        flux = np.random.default_rng(3).uniform(0.1, 1.0, size=(ng, nx, ny))

        fission_only = IntegratedReactionRate(
            solver.mat_xs.fission_production_field
        ).evaluate(flux)
        total = solver.compute_production_rate(flux)
        require(
            total > fission_only * (1.0 + 1e-9),
            f"(n,2n) must ADD to production when Sig2≠0: total={total!r}, "
            f"fission-only IRR={fission_only!r}. Guards the Mode-10 vacuity — every "
            f"xs_library mixture has Sig2=0, so a library fixture would null the term.",
        )
