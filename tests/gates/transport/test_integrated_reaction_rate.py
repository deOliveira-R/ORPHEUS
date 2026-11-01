r"""``IntegratedReactionRate`` — :math:`\int_V \langle \Sigma_x,\phi\rangle\,dV`.

The volume-integrated reaction rate (the k-eigenvalue numerator/denominator,
and the :math:`\phi^\dagger\!=\!1` degenerate of the homogenization PG bilinear).
Two legs here, a third in SN:

* **A — hand volume-weighted loop** on a **non-uniform-volume** mesh: the
  per-cell ``V`` weighting is only constrained when cells have distinct volumes
  (a uniform mesh nulls it — Mode-10). The reference is an explicit Python
  triple-loop ``Σ_cells V_cell Σ_g Σx_g φ_g`` (no numpy reduction shared with
  the SUT).
* **B — closed-form k∞ ratio**: on the homogeneous Mixture-A medium,
  ``IRR(νΣf)/IRR(Σa) == k∞`` (the volume measure CANCELS in the ratio — so this
  pins the fold is consistent between numerator and denominator). 4G mandatory
  (A 2G's spectrum is coincidentally flat).
* **C — (n,2n) activation in the SN production rate** moved to
  ``tests/gates/sn/operators/test_production_rate_counts_n2n.py`` (layer DAG Phase 4: an SN solver claim).

``foundation`` — L0 value correctness + software invariants.

Layer DAG Phase 4 (2026-11-01): the file builds material meshes and spells each space it reads from the raw
constructors of orpheus.transport.phase_spaces (the SN problem's and the diffusion mesh's cached mints are those
constructors on the problem's own ingredients, gated in
tests/gates/sn/mesh/test_the_mints_are_the_raw_constructors.py and
tests/gates/diffusion/test_scalar_trace_is_the_raw_constructor.py).
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.common.eigenvalue import kinf_and_spectrum_homogeneous
from orpheus.data.macro_xs.synthetic import get_mixture
from orpheus.geometry import StructuredGeometry
from orpheus.mesh import CellEdges, Mesher
from orpheus.numerics.functional import Functional
from orpheus.numerics.operator import LinearOperator
from orpheus.transport.mesh.material_mesh import MaterialMesh
from orpheus.transport.reaction_rate_functional import IntegratedReactionRate

from tests.gates.transport._materials import placeholder_materials
from tests.gates.transport._functional_helpers import cross_section_field, require, slab_mesh


pytestmark = pytest.mark.foundation


def _non_uniform_slab(ng: int) -> MaterialMesh:
    """A 1-D slab with NON-UNIFORM cell widths → distinct cell volumes.

    A uniform mesh makes the per-cell ``V`` a constant that the ratio/sum can't
    distinguish from 1 (Mode-10); distinct widths genuinely constrain the
    ``∫·dV`` weighting.
    """
    edges = np.array([0.0, 0.1, 0.3, 0.6, 1.0])  # widths 0.1 / 0.2 / 0.3 / 0.4
    geometry = StructuredGeometry.from_homogeneous(1.0, 0)
    mesh = Mesher(geometry).partition(CellEdges(edges)).mesh
    return MaterialMesh(mesh, placeholder_materials(ng=ng))

def _stretched_nonuniform_slab(ng: int) -> MaterialMesh:
    """The non-uniform slab at doubled width — same shape, different cell
    volumes (the F2 content discriminator)."""
    edges = 2.0 * np.array([0.0, 0.1, 0.3, 0.6, 1.0])
    geometry = StructuredGeometry.from_homogeneous(2.0, 0)
    mesh = Mesher(geometry).partition(CellEdges(edges)).mesh
    return MaterialMesh(mesh, placeholder_materials(ng=ng))


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
            phi, adjoint=ScalarFlux(values=phis, space=sn.bulk_space),
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
        ones = ScalarFlux(values=np.ones((2, nx)), space=sn.bulk_space)
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
        base = irr.evaluate(phi, adjoint=ScalarFlux(values=phis, space=sn.bulk_space))
        doubled = irr.evaluate(phi, adjoint=ScalarFlux(values=2.0 * phis, space=sn.bulk_space))
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
        ab = irr.evaluate(a, adjoint=ScalarFlux(values=b, space=sn.bulk_space))
        ba = irr.evaluate(b, adjoint=ScalarFlux(values=a, space=sn.bulk_space))
        np.testing.assert_allclose(ab, ba, rtol=1e-13)

    def test_foreign_mesh_adjoint_raises(self):
        """The σ↔geometry pairing tier: the adjoint must live on the SAME
        mesh OBJECT — an equal-shaped foreign mesh raises."""
        from orpheus.transport.fields.scalar_flux import ScalarFlux

        sn = _non_uniform_slab(ng=2)
        nx = sn.spatial_shape[0]
        rng = np.random.default_rng(17)
        irr = IntegratedReactionRate(
            cross_section_field(rng.uniform(0.1, 1.0, size=(2, nx)), sn)
        )
        # CS4b S3 (F2): a twin carrier's adjoint pairs legally (content).
        twin = ScalarFlux(values=rng.uniform(0.2, 2.0, size=(2, nx)), space=_non_uniform_slab(ng=2).bulk_space)
        irr.evaluate(rng.uniform(0.05, 1.0, size=(2, nx)), adjoint=twin)
        # A carrier whose volumes differ refuses — the σ↔geometry pairing.
        stretched = _stretched_nonuniform_slab(ng=2)
        foreign = ScalarFlux(values=rng.uniform(0.2, 2.0, size=(2, nx)), space=stretched.bulk_space)
        with pytest.raises(ValueError, match="not a member of"):
            irr.evaluate(rng.uniform(0.05, 1.0, size=(2, nx)), adjoint=foreign)


class TestTheReactionRateIsTheFlatOfItsCrossSection:
    """The space-membership unit U2c: the reaction rate is the canonical pairing of the co-vector ``V.flat(σ)`` with the
    flux, the metric read once at the mint. Against the old ``Σ_i V_i Σ_g σφ`` this re-associates one sum, so the claim
    is principled equivalence, not bit identity: leg A above is the independent pin (a hand triple loop), and the
    re-baselines it moved are declared in ``tests/gates/transport/carry_rebaselines.json`` (event U2c)."""

    def test_the_reaction_rate_pairs_the_flat_of_its_cross_section(self):
        sn = _non_uniform_slab(ng=3)
        sigx = np.random.default_rng(7).uniform(0.1, 1.0, size=(3, sn.spatial_shape[0]))
        field = cross_section_field(sigx, sn)
        covector = IntegratedReactionRate(field).covector
        assert covector.space == field.space.dual()
        np.testing.assert_array_equal(covector.values, field.space.flat(field).values)
