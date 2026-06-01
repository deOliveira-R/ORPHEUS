"""Foundation tests for the Issue #197 PR-TYPED-4 HarmonicMomentField.

Pins the structural contract of
:class:`~orpheus.sn.harmonic_moment_field.HarmonicMomentField` — shape
validation, named slicing / decomposition primitives, dunder algebra,
truncation, ``scalar_flux()`` agreement with the bare-ordinate
reduction, and the SN-side typed ``R \\cdot \\Lambda \\cdot M``
round-trip.

These are foundation tests — they verify software invariants
(constructor shape check, frozen-ness, dunder algebra) and the
algebraic identities the type carries; they do NOT make physics
claims, so they carry ``@pytest.mark.foundation`` per the V&V harness
convention.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry import BC, CoordSystem, Mesh1D, Mesh2D
from orpheus.sn.geometry import SNMesh
from orpheus.transport.fields.angular_flux import AngularFlux
from orpheus.transport.fields.harmonic_moment_field import HarmonicMomentField
from orpheus.numerics.quadrature import Quadrature
from orpheus.transport.fields.scalar_flux import ScalarFlux

from tests.sn._test_helpers import placeholder_materials

pytestmark = pytest.mark.foundation


# ── Fixtures ─────────────────────────────────────────────────────────


def _slab_mesh(nx: int = 4, ng: int = 2) -> SNMesh:
    """Build a small slab :class:`SNMesh` for unit testing."""
    mesh = Mesh1D(
        edges=np.linspace(0.0, 1.0, nx + 1),
        mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.CARTESIAN,
        bc_left=BC("vacuum"),
        bc_right=BC("vacuum"),
    )
    quad = Quadrature.gauss_legendre(n_ordinates=4)
    return SNMesh(mesh, quad, placeholder_materials(ng=ng))


def _2d_mesh(nx: int = 3, ny: int = 3, ng: int = 1) -> SNMesh:
    """Build a small 2-D Cartesian :class:`SNMesh`."""
    mesh = Mesh2D(
        edges_x=np.linspace(0, 1, nx + 1),
        edges_y=np.linspace(0, 1, ny + 1),
        mat_map=np.zeros((nx, ny), dtype=int),
    )
    quad = Quadrature.level_symmetric(sn_order=4)
    return SNMesh(mesh, quad, placeholder_materials(ng=ng))


# ════════════════════════════════════════════════════════════════════
# Construction + shape validation
# ════════════════════════════════════════════════════════════════════


class TestHarmonicMomentFieldConstruction:
    def test_construct_from_factory(self) -> None:
        m = _slab_mesh()
        L = 2
        phi = m.zeros_harmonic_moments(L)
        # (L+1, 2L+1, ng, nx, ny) = (3, 5, 2, 4, 1)
        assert phi.values.shape == (L + 1, 2 * L + 1, m.ng, m.nx, m.ny)
        assert np.all(phi.values == 0.0)
        assert isinstance(phi, HarmonicMomentField)
        assert phi.mesh is m
        assert phi.L == L

    def test_construct_explicit(self) -> None:
        m = _2d_mesh(ng=1)
        L = 1
        # shape (2, 3, 1, 3, 3)
        vals = np.arange(2 * 3 * 1 * 3 * 3, dtype=float).reshape(
            (2, 3, 1, 3, 3),
        )
        phi = HarmonicMomentField.from_mesh_and_L(vals, m, L)
        assert phi.L == 1
        np.testing.assert_array_equal(phi.values, vals)

    def test_shape_validation_wrong_L_block_size(self) -> None:
        m = _slab_mesh()
        L = 2
        # Wrong shape: (L+1, 2L, ...) instead of (L+1, 2L+1, ...).
        with pytest.raises(ValueError, match="HarmonicMomentField.*does not match"):
            HarmonicMomentField.from_mesh_and_L(
                np.zeros((L + 1, 2 * L, m.ng, m.nx, m.ny)), m, L,
            )

    def test_shape_validation_wrong_mesh_dims(self) -> None:
        m = _slab_mesh()
        L = 1
        with pytest.raises(ValueError, match="HarmonicMomentField.*does not match"):
            HarmonicMomentField.from_mesh_and_L(
                np.zeros((L + 1, 2 * L + 1, m.ng + 1, m.nx, m.ny)), m, L,
            )

    def test_metadata_read_throughs(self) -> None:
        m = _slab_mesh()
        phi = m.zeros_harmonic_moments(L=0)
        assert phi.ng == m.ng
        assert phi.nx == m.nx
        assert phi.ny == m.ny

    def test_space_is_tensor_product_of_sh_and_cell_group(self) -> None:
        r"""D-E invariant: the typed field's ``space`` is a
        :class:`TensorProductSpace` whose factors are
        ``SphericalHarmonicSpace(L)`` and a cell-group
        :class:`FunctionSpace`.

        This is the FIRST production consumer of D-B's
        :class:`TensorProductSpace` primitive in a typed transport
        Field. The moment-axis structure is type-visible via the
        composition tree — consumers asking for ``L`` via the
        composition-tree query (Issue #207) get the right answer:
        ``phi.space.factors[0].L == phi.L``.
        """
        from orpheus.numerics.space import (
            FunctionSpace, TensorProductSpace,
        )
        from orpheus.numerics.spaces.spherical_harmonic_space import (
            SphericalHarmonicSpace,
        )
        m = _slab_mesh()
        L = 2
        phi = m.zeros_harmonic_moments(L=L)
        # Space typing invariants.
        assert isinstance(phi.space, TensorProductSpace)
        assert len(phi.space.factors) == 2
        # Factor 0: SphericalHarmonicSpace carrying the L parameter.
        sh_factor = phi.space.factors[0]
        assert isinstance(sh_factor, SphericalHarmonicSpace)
        assert sh_factor.L == L
        assert sh_factor.shape == (L + 1, 2 * L + 1)
        # Factor 1: cell-group FunctionSpace.
        cg_factor = phi.space.factors[1]
        assert isinstance(cg_factor, FunctionSpace)
        assert cg_factor.shape == (m.ng, m.nx, m.ny)
        # Combined shape matches the field's data shape.
        assert phi.space.shape == phi.values.shape


# ════════════════════════════════════════════════════════════════════
# Slicing / decomposition primitives (Pattern 3 — named intermediates)
# ════════════════════════════════════════════════════════════════════


class TestHarmonicMomentFieldSlicing:
    def test_l_block_returns_view_with_right_shape(self) -> None:
        m = _slab_mesh()
        L = 2
        phi = m.zeros_harmonic_moments(L)
        for l in range(L + 1):
            block = phi.l_block(l)
            assert block.shape == (2 * l + 1, m.ng, m.nx, m.ny)

    def test_l_block_values_match_underlying_slice(self) -> None:
        m = _slab_mesh()
        L = 2
        rng = np.random.default_rng(seed=0)
        vals = rng.standard_normal(
            (L + 1, 2 * L + 1, m.ng, m.nx, m.ny),
        )
        phi = HarmonicMomentField.from_mesh_and_L(vals, m, L)
        for l in range(L + 1):
            np.testing.assert_array_equal(
                phi.l_block(l), vals[l, : 2 * l + 1],
            )

    def test_l_block_out_of_range_raises(self) -> None:
        m = _slab_mesh()
        phi = m.zeros_harmonic_moments(L=1)
        with pytest.raises(ValueError):
            phi.l_block(2)
        with pytest.raises(ValueError):
            phi.l_block(-1)

    def test_isotropic_part_zeros_l_ge_1(self) -> None:
        m = _slab_mesh()
        L = 2
        rng = np.random.default_rng(seed=1)
        vals = rng.standard_normal(
            (L + 1, 2 * L + 1, m.ng, m.nx, m.ny),
        )
        phi = HarmonicMomentField.from_mesh_and_L(vals, m, L)
        iso = phi.isotropic_part()
        assert isinstance(iso, HarmonicMomentField)
        assert iso.L == L
        # ℓ=0 m=0 slot preserved.
        np.testing.assert_array_equal(iso.values[0, 0], vals[0, 0])
        # Every other slot zeroed.
        np.testing.assert_array_equal(iso.values[1:], 0.0)
        # The ℓ=0 row's m≥1 slots (zero-padding) stay zero.
        np.testing.assert_array_equal(iso.values[0, 1:], 0.0)

    def test_anisotropic_part_zeros_l_eq_0(self) -> None:
        m = _slab_mesh()
        L = 2
        rng = np.random.default_rng(seed=2)
        vals = rng.standard_normal(
            (L + 1, 2 * L + 1, m.ng, m.nx, m.ny),
        )
        phi = HarmonicMomentField.from_mesh_and_L(vals, m, L)
        aniso = phi.anisotropic_part()
        assert isinstance(aniso, HarmonicMomentField)
        # ℓ=0 m=0 slot zeroed.
        np.testing.assert_array_equal(aniso.values[0, 0], 0.0)
        # ℓ ≥ 1 preserved.
        np.testing.assert_array_equal(aniso.values[1:], vals[1:])

    def test_iso_plus_aniso_recovers_self(self) -> None:
        m = _slab_mesh()
        L = 2
        rng = np.random.default_rng(seed=3)
        vals = rng.standard_normal(
            (L + 1, 2 * L + 1, m.ng, m.nx, m.ny),
        )
        # Zero the ℓ=0 m≥1 padding so iso+aniso reproduces vals
        # exactly (the padding is not part of the legitimate
        # moment-space content).
        vals[0, 1:] = 0.0
        phi = HarmonicMomentField.from_mesh_and_L(vals, m, L)
        recombined = phi.isotropic_part() + phi.anisotropic_part()
        np.testing.assert_array_equal(recombined.values, vals)


# ════════════════════════════════════════════════════════════════════
# scalar_flux() — extracts ℓ=0, m=0 moment as a ScalarFlux
# ════════════════════════════════════════════════════════════════════


class TestHarmonicMomentFieldScalarFlux:
    def test_scalar_flux_returns_ScalarFlux_with_isotropic_slot(self) -> None:
        m = _slab_mesh()
        L = 1
        rng = np.random.default_rng(seed=4)
        vals = rng.standard_normal(
            (L + 1, 2 * L + 1, m.ng, m.nx, m.ny),
        )
        phi = HarmonicMomentField.from_mesh_and_L(vals, m, L)
        sf = phi.scalar_flux()
        assert isinstance(sf, ScalarFlux)
        assert sf.values.shape == (m.ng, m.nx, m.ny)
        np.testing.assert_array_equal(sf.values, vals[0, 0])

    def test_scalar_flux_agrees_with_integrate_angular(self) -> None:
        """``M(\\psi).scalar_flux() == \\psi.integrate_angular()``.

        Under the no-prefactor SH convention (``Y_0^0 = 1``), the
        isotropic moment :math:`\\sum_n w_n Y_0^0 \\psi_n` reduces to
        :math:`\\sum_n w_n \\psi_n` — the bare angular reduction.  This
        identity makes the moment-space and ordinate-space scalar
        fluxes algebraically equivalent.
        """
        from orpheus.numerics.projection import MomentProjection
        m = _slab_mesh()
        L = 1
        rng = np.random.default_rng(seed=5)
        N = m.quad.N
        psi_values = rng.standard_normal((N, m.ng, m.nx, m.ny))
        psi = AngularFlux.from_mesh(psi_values, m)

        # Direct angular reduction.
        sf_direct = psi.integrate_angular()

        # Via moment projection + scalar_flux extraction.
        Y = m.quad.spherical_harmonics(L)
        M = MomentProjection(
            weights=m.quad.weights, Y=Y, L=L,
        )
        moments_values = M.apply(psi.values)
        moments = HarmonicMomentField.from_mesh_and_L(moments_values, m, L)
        sf_via_moments = moments.scalar_flux()

        np.testing.assert_allclose(
            sf_via_moments.values, sf_direct.values, rtol=1e-13,
        )


# ════════════════════════════════════════════════════════════════════
# truncate
# ════════════════════════════════════════════════════════════════════


class TestHarmonicMomentFieldTruncate:
    def test_truncate_preserves_lower_blocks(self) -> None:
        m = _slab_mesh()
        L = 3
        rng = np.random.default_rng(seed=6)
        vals = rng.standard_normal(
            (L + 1, 2 * L + 1, m.ng, m.nx, m.ny),
        )
        phi = HarmonicMomentField.from_mesh_and_L(vals, m, L)
        for L_new in range(L + 1):
            trunc = phi.truncate(L_new)
            assert trunc.L == L_new
            assert trunc.values.shape == (
                L_new + 1, 2 * L_new + 1, m.ng, m.nx, m.ny,
            )
            # For each ℓ ≤ L_new, the (ℓ, m) entries inside |m|≤ℓ
            # match the source.
            for l in range(L_new + 1):
                np.testing.assert_array_equal(
                    trunc.l_block(l), phi.l_block(l),
                )

    def test_truncate_rejects_L_new_greater_than_L(self) -> None:
        m = _slab_mesh()
        phi = m.zeros_harmonic_moments(L=1)
        with pytest.raises(ValueError, match="truncate"):
            phi.truncate(2)

    def test_truncate_rejects_negative(self) -> None:
        m = _slab_mesh()
        phi = m.zeros_harmonic_moments(L=1)
        with pytest.raises(ValueError, match="truncate"):
            phi.truncate(-1)


# ════════════════════════════════════════════════════════════════════
# Dunder algebra (Pattern 1 — match the math)
# ════════════════════════════════════════════════════════════════════


class TestHarmonicMomentFieldAlgebra:
    def test_add_within_type(self) -> None:
        m = _slab_mesh()
        L = 1
        shape = (L + 1, 2 * L + 1, m.ng, m.nx, m.ny)
        a = HarmonicMomentField.from_mesh_and_L(np.ones(shape), m, L)
        b = HarmonicMomentField.from_mesh_and_L(2.0 * np.ones(shape), m, L)
        c = a + b
        assert isinstance(c, HarmonicMomentField)
        np.testing.assert_array_equal(c.values, 3.0 * np.ones(shape))
        # Frozen — originals unchanged.
        np.testing.assert_array_equal(a.values, np.ones(shape))

    def test_sub_within_type(self) -> None:
        m = _slab_mesh()
        L = 1
        shape = (L + 1, 2 * L + 1, m.ng, m.nx, m.ny)
        a = HarmonicMomentField.from_mesh_and_L(3.0 * np.ones(shape), m, L)
        b = HarmonicMomentField.from_mesh_and_L(np.ones(shape), m, L)
        c = a - b
        assert isinstance(c, HarmonicMomentField)
        np.testing.assert_array_equal(c.values, 2.0 * np.ones(shape))

    def test_scalar_mul_left_and_right(self) -> None:
        m = _slab_mesh()
        L = 1
        shape = (L + 1, 2 * L + 1, m.ng, m.nx, m.ny)
        a = HarmonicMomentField.from_mesh_and_L(np.ones(shape), m, L)
        np.testing.assert_array_equal(
            (3.0 * a).values, (a * 3.0).values,
        )
        np.testing.assert_array_equal(
            (3.0 * a).values, 3.0 * np.ones(shape),
        )

    def test_div(self) -> None:
        m = _slab_mesh()
        L = 0
        shape = (L + 1, 2 * L + 1, m.ng, m.nx, m.ny)
        a = HarmonicMomentField.from_mesh_and_L(2.0 * np.ones(shape), m, L)
        np.testing.assert_array_equal(
            (a / 2.0).values, np.ones(shape),
        )

    def test_neg(self) -> None:
        m = _slab_mesh()
        L = 1
        shape = (L + 1, 2 * L + 1, m.ng, m.nx, m.ny)
        a = HarmonicMomentField.from_mesh_and_L(np.ones(shape), m, L)
        np.testing.assert_array_equal((-a).values, -np.ones(shape))

    def test_partner_must_be_same_type(self) -> None:
        m = _slab_mesh()
        L = 1
        shape = (L + 1, 2 * L + 1, m.ng, m.nx, m.ny)
        a = HarmonicMomentField.from_mesh_and_L(np.ones(shape), m, L)
        with pytest.raises(TypeError):
            a + 5  # not a HarmonicMomentField

    def test_partner_must_share_mesh(self) -> None:
        m1 = _slab_mesh()
        m2 = _slab_mesh()  # distinct instance — same shape
        L = 1
        shape = (L + 1, 2 * L + 1, m1.ng, m1.nx, m1.ny)
        a = HarmonicMomentField.from_mesh_and_L(np.ones(shape), m1, L)
        b = HarmonicMomentField.from_mesh_and_L(np.ones(shape), m2, L)
        with pytest.raises(ValueError, match="mesh-bound"):
            a + b

    def test_partner_must_share_L(self) -> None:
        m = _slab_mesh()
        a = m.zeros_harmonic_moments(L=1)
        b = m.zeros_harmonic_moments(L=2)
        # Post-D-E: L mismatch surfaces as a space-equality error,
        # because different L values produce different
        # SphericalHarmonicSpace shapes, which propagates to different
        # TensorProductSpace identities. Field._check_partner's space
        # check fires before the explicit L check in
        # HarmonicMomentField._check_partner. Both gate the same
        # invariant; the space-level message is more general.
        with pytest.raises(ValueError, match="equal space"):
            a + b


# ════════════════════════════════════════════════════════════════════
# R · Λ · M · ψ round-trip through the typed pipeline
# ════════════════════════════════════════════════════════════════════


class TestRLambdaMRoundTrip:
    def test_aniso_part_zero_under_isotropic_psi(self) -> None:
        r"""For isotropic :math:`\psi_n = c \forall n`, the
        Pℓ≥1 reconstruction :math:`R \cdot \Lambda \cdot M \cdot \psi`
        is identically zero (the anisotropic moments vanish by
        construction).  This test pins the typed-pipeline output as
        a :class:`AngularSource` matching that algebraic claim.
        """
        from orpheus.derivations.common.xs_library import get_mixture
        from orpheus.sn.scattering import ScatteringOperator
        from orpheus.sn.solver import SNSolver
        from orpheus.transport.sources import AngularSource
        mix = get_mixture("A", "2g")
        if len(mix.SigS) < 2:
            pytest.skip("No P1 data in test mixture")

        nx, ny = 2, 2
        mesh = Mesh2D(
            edges_x=np.linspace(0, 1, nx + 1),
            edges_y=np.linspace(0, 1, ny + 1),
            mat_map=np.zeros((nx, ny), dtype=int),
        )
        quad = Quadrature.level_symmetric(sn_order=4)
        sn_mesh = SNMesh(mesh, quad, {0: mix})
        solver = SNSolver(sn_mesh, scattering_order=1)
        op = solver.scattering_op

        # Build a typed AngularFlux with isotropic content.
        N = quad.N
        psi_values = np.ones((N, mix.ng, nx, ny))
        psi = AngularFlux.from_mesh(psi_values, sn_mesh)

        # Typed pipeline output.
        out = op.build_aniso_source(psi)
        # Must be AngularSource (not bare ndarray) under typed-in.
        assert isinstance(out, AngularSource)
        assert out.values.shape == (N, mix.ng, nx, ny)
        np.testing.assert_allclose(out.values, 0.0, atol=1e-12)

    def test_lambda_apply_typed_in_typed_out(self) -> None:
        """``LegendreMomentScattering.apply(HarmonicMomentField)``
        returns ``HarmonicMomentField`` with matching mesh + L."""
        from orpheus.derivations.common.xs_library import get_mixture
        from orpheus.sn.scattering import (
            LegendreMomentScattering, ScatteringOperator,
        )
        from orpheus.sn.solver import SNSolver
        mix = get_mixture("A", "2g")
        if len(mix.SigS) < 2:
            pytest.skip("No P1 data in test mixture")

        nx, ny = 2, 2
        mesh = Mesh2D(
            edges_x=np.linspace(0, 1, nx + 1),
            edges_y=np.linspace(0, 1, ny + 1),
            mat_map=np.zeros((nx, ny), dtype=int),
        )
        quad = Quadrature.level_symmetric(sn_order=4)
        sn_mesh = SNMesh(mesh, quad, {0: mix})
        solver = SNSolver(sn_mesh, scattering_order=1)
        op = solver.scattering_op

        L = 1
        rng = np.random.default_rng(seed=7)
        moments_values = rng.standard_normal(
            (L + 1, 2 * L + 1, mix.ng, nx, ny),
        )
        moments = HarmonicMomentField.from_mesh_and_L(moments_values, sn_mesh, L)

        Lam = LegendreMomentScattering(
            mat_xs=op.mat_xs, L=L, skip_l0=True,
        )
        out = Lam.apply(moments)
        assert isinstance(out, HarmonicMomentField)
        assert out.mesh is sn_mesh
        assert out.L == L
        assert out.values.shape == moments.values.shape

    def test_lambda_bare_in_bare_out_legacy_path(self) -> None:
        """Bare-ndarray path is preserved for legacy probe tests."""
        from orpheus.derivations.common.xs_library import get_mixture
        from orpheus.sn.scattering import LegendreMomentScattering
        from orpheus.sn.solver import SNSolver
        mix = get_mixture("A", "2g")
        if len(mix.SigS) < 2:
            pytest.skip("No P1 data in test mixture")

        nx, ny = 2, 2
        mesh = Mesh2D(
            edges_x=np.linspace(0, 1, nx + 1),
            edges_y=np.linspace(0, 1, ny + 1),
            mat_map=np.zeros((nx, ny), dtype=int),
        )
        quad = Quadrature.level_symmetric(sn_order=4)
        sn_mesh = SNMesh(mesh, quad, {0: mix})
        solver = SNSolver(sn_mesh, scattering_order=1)
        op = solver.scattering_op

        L = 1
        moments_values = np.zeros(
            (L + 1, 2 * L + 1, mix.ng, nx, ny),
        )
        Lam = LegendreMomentScattering(
            mat_xs=op.mat_xs, L=L, skip_l0=True,
        )
        out = Lam.apply(moments_values)
        assert isinstance(out, np.ndarray)
        assert out.shape == moments_values.shape


# ════════════════════════════════════════════════════════════════════
# Factory + zeros_harmonic_moments
# ════════════════════════════════════════════════════════════════════


class TestSNMeshFactory:
    def test_factory_for_L_zero(self) -> None:
        m = _slab_mesh()
        phi = m.zeros_harmonic_moments(L=0)
        assert phi.L == 0
        assert phi.values.shape == (1, 1, m.ng, m.nx, m.ny)

    def test_factory_returns_owned_ndarray(self) -> None:
        m = _slab_mesh()
        phi1 = m.zeros_harmonic_moments(L=1)
        phi2 = m.zeros_harmonic_moments(L=1)
        # Independent allocations.
        assert phi1.values is not phi2.values
        phi1.values.flags.writeable
        # Mesh shared — by reference.
        assert phi1.mesh is phi2.mesh

    def test_copy_creates_independent(self) -> None:
        m = _slab_mesh()
        phi = m.zeros_harmonic_moments(L=1)
        phi_copy = phi.copy()
        assert phi.values is not phi_copy.values
        np.testing.assert_array_equal(phi.values, phi_copy.values)
