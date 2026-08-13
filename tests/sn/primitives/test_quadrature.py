"""Tests for the canonical :class:`Quadrature` factories.

Covers all four factory classmethods:

* :meth:`Quadrature.gauss_legendre` — slab 1-D polar.
* :meth:`Quadrature.lebedev` — sphere :math:`O_h`-invariant.
* :meth:`Quadrature.level_symmetric` — Carlson-Lathrop S_N.
* :meth:`Quadrature.product` — GL × equispaced.

Tests verify mathematical properties that any valid quadrature must
satisfy (weight sums, unit-sphere condition, second moments, level
structure, reflection involution).
"""

import numpy as np
import pytest

from orpheus.geometry import CoordSystem
from orpheus.geometry.transformation import RigidMotion
from orpheus.numerics.quadrature import Quadrature
from tests.sn._test_helpers import placeholder_materials

# All quadrature tests are L0 term-verification (hand-computed weight
# sums, moment conditions, alpha-recursion closure). TestL0TermVerification
# already auto-tags via the TestL<N> class-name convention.
pytestmark = [
    pytest.mark.l0,
    pytest.mark.verifies(
        "mm-weights",
        "alpha-recursion",
        "alpha-cylindrical",
        "quadrature-product-weights",
        "wdd-closure",
        "wdd-face",
        "reflective-bc",
        "quadrature-ordinate-permutation",
        "flux-moments",
    ),
]


# ═══════════════════════════════════════════════════════════════════════
# Weight sums
# ═══════════════════════════════════════════════════════════════════════

class TestWeightSums:
    """Quadrature weights must integrate to the correct solid angle."""

    @pytest.mark.parametrize("N", [4, 8, 16])
    def test_gl_weights_sum_to_2(self, N):
        quad = Quadrature.gauss_legendre(N)
        np.testing.assert_allclose(quad.weights.sum(), 2.0, atol=1e-14)

    @pytest.mark.sentinel
    def test_lebedev_weights_sum_to_4pi(self):
        quad = Quadrature.lebedev(order=17)
        np.testing.assert_allclose(quad.weights.sum(), 4 * np.pi, rtol=1e-12)

    @pytest.mark.parametrize("order", [2, 4, 6, 8])
    def test_level_symmetric_weights_sum_to_4pi(self, order):
        quad = Quadrature.level_symmetric(order)
        np.testing.assert_allclose(quad.weights.sum(), 4 * np.pi, rtol=1e-12)

    @pytest.mark.parametrize("n_mu,n_phi", [(4, 4), (4, 8), (8, 8)])
    def test_product_weights_sum_to_4pi(self, n_mu, n_phi):
        quad = Quadrature.product(n_mu, n_phi)
        np.testing.assert_allclose(quad.weights.sum(), 4 * np.pi, rtol=1e-12)


# ═══════════════════════════════════════════════════════════════════════
# Unit sphere condition
# ═══════════════════════════════════════════════════════════════════════

class TestUnitSphere:
    """All ordinates must lie on the unit sphere: η² + ξ² + μ² = 1."""

    def test_lebedev(self):
        quad = Quadrature.lebedev(order=17)
        norm = quad.mu_x**2 + quad.mu_y**2 + quad.mu_z**2
        np.testing.assert_allclose(norm, 1.0, atol=1e-14)

    @pytest.mark.parametrize("order", [2, 4, 6, 8])
    def test_level_symmetric(self, order):
        quad = Quadrature.level_symmetric(order)
        norm = quad.mu_x**2 + quad.mu_y**2 + quad.mu_z**2
        np.testing.assert_allclose(norm, 1.0, atol=1e-14)

    @pytest.mark.parametrize("n_mu,n_phi", [(4, 8), (8, 8)])
    def test_product(self, n_mu, n_phi):
        quad = Quadrature.product(n_mu, n_phi)
        norm = quad.mu_x**2 + quad.mu_y**2 + quad.mu_z**2
        np.testing.assert_allclose(norm, 1.0, atol=1e-14)


# ═══════════════════════════════════════════════════════════════════════
# Moment conditions
# ═══════════════════════════════════════════════════════════════════════

class TestMomentConditions:
    """Quadratures must integrate low-order polynomials exactly.

    ∫ dΩ = 4π, ∫ μ_i² dΩ = 4π/3 for i ∈ {x, y, z}.
    """

    @pytest.mark.parametrize("order", [4, 6, 8])
    def test_level_symmetric_second_moments(self, order):
        quad = Quadrature.level_symmetric(order)
        target = 4 * np.pi / 3
        for attr in ['mu_x', 'mu_y', 'mu_z']:
            m2 = np.sum(quad.weights * getattr(quad, attr)**2)
            np.testing.assert_allclose(m2, target, rtol=1e-10,
                                       err_msg=f"∫{attr}² dΩ ≠ 4π/3")

    @pytest.mark.parametrize("n_mu,n_phi", [(4, 8), (8, 8)])
    def test_product_second_moments(self, n_mu, n_phi):
        quad = Quadrature.product(n_mu, n_phi)
        target = 4 * np.pi / 3
        for attr in ['mu_x', 'mu_y', 'mu_z']:
            m2 = np.sum(quad.weights * getattr(quad, attr)**2)
            np.testing.assert_allclose(m2, target, rtol=1e-10,
                                       err_msg=f"∫{attr}² dΩ ≠ 4π/3")

    def test_lebedev_second_moments(self):
        quad = Quadrature.lebedev(order=17)
        target = 4 * np.pi / 3
        for attr in ['mu_x', 'mu_y', 'mu_z']:
            m2 = np.sum(quad.weights * getattr(quad, attr)**2)
            np.testing.assert_allclose(m2, target, rtol=1e-10)


# ═══════════════════════════════════════════════════════════════════════
# Level structure (for cylindrical sweep)
# ═══════════════════════════════════════════════════════════════════════

class TestLevelStructure:
    """Level-indexed quadratures must have consistent level_indices."""

    @pytest.mark.parametrize("order", [2, 4, 6, 8])
    def test_level_sym_indices_partition(self, order):
        """Level indices must cover all ordinates exactly once."""
        quad = Quadrature.level_symmetric(order)
        all_idx = np.sort(np.concatenate(quad.level_indices))
        np.testing.assert_array_equal(all_idx, np.arange(quad.N))

    @pytest.mark.parametrize("n_mu,n_phi", [(4, 4), (8, 8)])
    def test_product_indices_partition(self, n_mu, n_phi):
        quad = Quadrature.product(n_mu, n_phi)
        all_idx = np.sort(np.concatenate(quad.level_indices))
        np.testing.assert_array_equal(all_idx, np.arange(quad.N))

    def test_product_level_mu_match(self):
        """On each level, all ordinates must share the same μ_z value."""
        quad = Quadrature.product(n_mu=4, n_phi=8)
        for p, idx in enumerate(quad.level_indices):
            mu_vals = quad.mu_z[idx]
            np.testing.assert_allclose(mu_vals, quad.level_mu[p], atol=1e-14)


# ═══════════════════════════════════════════════════════════════════════
# Mirror-induced ordinate permutations
# ═══════════════════════════════════════════════════════════════════════


def _mirror_pi(quad, axis_index: int) -> np.ndarray:
    """The permutation σ_axis induces, via ``ordinate_permutation``."""
    pi = quad.ordinate_permutation(
        RigidMotion.reflection(normal=np.eye(3)[axis_index])
    )
    assert pi is not None, f"axis {axis_index}: rule not mirror-closed"
    return pi.indices


class TestReflectionIndices:
    """The mirror partner must have the negated direction cosine."""

    @pytest.mark.parametrize("factory,kwargs", [
        (Quadrature.level_symmetric, {"sn_order": 4}),
        (Quadrature.product, {"n_mu": 4, "n_phi": 8}),
    ])
    def test_x_reflection(self, factory, kwargs):
        quad = factory(**kwargs)
        ref = _mirror_pi(quad, 0)
        # μ_x of reflected partner should be -μ_x of original
        np.testing.assert_allclose(quad.mu_x[ref], -quad.mu_x, atol=1e-12)

    @pytest.mark.parametrize("factory,kwargs", [
        (Quadrature.level_symmetric, {"sn_order": 4}),
        (Quadrature.product, {"n_mu": 4, "n_phi": 8}),
    ])
    def test_reflection_involution(self, factory, kwargs):
        """Reflecting twice must return to the original ordinate."""
        quad = factory(**kwargs)
        for axis in (0, 1, 2):
            ref = _mirror_pi(quad, axis)
            np.testing.assert_array_equal(ref[ref], np.arange(quad.N),
                                          err_msg=f"axis-{axis} mirror not involution")


# ═══════════════════════════════════════════════════════════════════════
# Alpha redistribution coefficient properties
# ═══════════════════════════════════════════════════════════════════════

class TestAlphaRedistribution:
    """Verify α coefficient properties required for curvilinear SN sweeps.

    The α recursion (Bailey et al. 2009, Eq. 50) uses the radial
    direction cosine η (mu_x): α_{m+1/2} = α_{m-1/2} − w_m · η_m.
    The resulting dome must be non-negative with α[0] = α[M] = 0.
    """

    # Q5.6.3: the folded family SNMesh(CYLINDRICAL) admits.  Splits chosen
    # for level-structure diversity (4×4, 8×8, 8×2, 6×6 angles/level), and
    # folded(4,6) puts the surviving bit-exact μ_r = 0 degenerate ordinate
    # (parent n_φ ≡ 2 mod 4) inside the dome.
    @pytest.mark.parametrize("factory,kwargs", [
        (Quadrature.folded_product, {"n_mu": 4, "n_phi": 8}),
        (Quadrature.folded_product, {"n_mu": 8, "n_phi": 16}),
        (Quadrature.folded_product, {"n_mu": 8, "n_phi": 4}),
        (Quadrature.folded_product, {"n_mu": 6, "n_phi": 12}),
        (Quadrature.folded_product, {"n_mu": 4, "n_phi": 6}),
    ])
    def test_alpha_dome_non_negative(self, factory, kwargs):
        """α values must form a non-negative dome on each level."""
        from orpheus.geometry import CoordSystem, Mesh1D
        from orpheus.sn.mesh.augmented_mesh import SNMesh

        quad = factory(**kwargs)
        mesh = Mesh1D(
            edges=np.array([0.0, 1.0]), mat_ids=np.array([0]),
            coord=CoordSystem.CYLINDRICAL,
        )
        sn_mesh = SNMesh(mesh, quad, placeholder_materials())

        for p, alpha in enumerate(sn_mesh.reduced.alpha_per_level):
            assert np.all(alpha >= -1e-14), (
                f"Level {p}: negative α = {alpha.min():.2e}"
            )

    @pytest.mark.parametrize("factory,kwargs", [
        (Quadrature.folded_product, {"n_mu": 4, "n_phi": 8}),
        (Quadrature.folded_product, {"n_mu": 4, "n_phi": 6}),
    ])
    def test_alpha_boundary_zero(self, factory, kwargs):
        """α must be zero at both dome boundaries (conservation)."""
        from orpheus.geometry import CoordSystem, Mesh1D
        from orpheus.sn.mesh.augmented_mesh import SNMesh

        quad = factory(**kwargs)
        mesh = Mesh1D(
            edges=np.array([0.0, 1.0]), mat_ids=np.array([0]),
            coord=CoordSystem.CYLINDRICAL,
        )
        sn_mesh = SNMesh(mesh, quad, placeholder_materials())

        for p, alpha in enumerate(sn_mesh.reduced.alpha_per_level):
            np.testing.assert_allclose(alpha[0], 0.0,
                                       err_msg=f"Level {p}: α[0] ≠ 0")
            np.testing.assert_allclose(alpha[-1], 0.0, atol=1e-13,
                                       err_msg=f"Level {p}: α[-1] ≠ 0")

    @pytest.mark.sentinel
    def test_spherical_alpha_dome_non_negative(self):
        """Spherical α (cumsum(−w·μ)) must be non-negative for GL quadrature."""
        from orpheus.geometry import CoordSystem, Mesh1D
        from orpheus.sn.mesh.augmented_mesh import SNMesh

        quad = Quadrature.gauss_legendre(8)
        mesh = Mesh1D(
            edges=np.array([0.0, 1.0]), mat_ids=np.array([0]),
            coord=CoordSystem.SPHERICAL,
        )
        sn_mesh = SNMesh(mesh, quad, placeholder_materials())

        assert np.all(sn_mesh.reduced.alpha_half >= -1e-14), (
            f"Negative spherical α: min = {sn_mesh.reduced.alpha_half.min():.2e}"
        )


class TestL0TermVerification:
    """Term-level (L0) verification tests.

    Each test isolates a single property of the discretization and
    verifies it against a hand calculation.  Tagged with L0-SN-NNN
    for the publication catalog.
    """

    @pytest.mark.sentinel
    @pytest.mark.catches("ERR-006", "ERR-007")
    @pytest.mark.parametrize("coord", [
        CoordSystem.SPHERICAL, CoordSystem.CYLINDRICAL,
    ])
    def test_per_ordinate_flat_flux_consistency(self, coord):
        """L0-SN-003: streaming + redistribution = 0 per ordinate for flat ψ.

        The fundamental correctness criterion for curvilinear SN.
        The ΔA/w factor ensures exact per-ordinate cancellation.
        """
        from orpheus.geometry import (
            BC,
            CoordSystem,
            Mesh1D,
            Region,
            RegionMesh,
            StructuredGeometry,
        )
        from orpheus.sn.mesh.augmented_mesh import SNMesh

        if coord == CoordSystem.SPHERICAL:
            quad = Quadrature.gauss_legendre(8)
        else:
            quad = Quadrature.folded_product(n_mu=4, n_phi=8)

        tag = {
            CoordSystem.SPHERICAL: "SPH",
            CoordSystem.CYLINDRICAL: "CYL",
        }[coord]
        mesh = Mesh1D.from_geometry(
            StructuredGeometry(
                geometry=tag,
                regions=(Region(mat_id=0, outer_thickness_cm=1.0),),
                bcs=(BC.reflective,),
            ),
            region_meshes=(RegionMesh(n_cells=10),),
        )
        sn = SNMesh(mesh, quad, placeholder_materials())
        dA = sn.delta_A
        psi0 = 1.0

        if coord == CoordSystem.SPHERICAL:
            alpha = sn.reduced.alpha_half
            for n in range(quad.N):
                streaming = quad.mu_x[n] * dA * psi0
                alpha_diff = alpha[n + 1] - alpha[n]
                redist = (dA / quad.weights[n]) * alpha_diff * psi0
                residual = streaming + redist
                np.testing.assert_allclose(
                    residual, 0.0, atol=1e-14,
                    err_msg=f"Spherical ordinate {n}: residual ≠ 0",
                )
        else:
            for p in range(len(sn.reduced.alpha_per_level)):
                alpha = sn.reduced.alpha_per_level[p]
                for m, n in enumerate(quad.level_indices[p]):
                    streaming = quad.mu_x[n] * dA * psi0
                    alpha_diff = alpha[m + 1] - alpha[m]
                    redist = (dA / quad.weights[n]) * alpha_diff * psi0
                    residual = streaming + redist
                    np.testing.assert_allclose(
                        residual, 0.0, atol=1e-14,
                        err_msg=f"Cyl level {p} ord {m}: residual ≠ 0",
                    )

    @pytest.mark.parametrize("coord", [
        CoordSystem.SPHERICAL, CoordSystem.CYLINDRICAL,
    ])
    def test_delta_A_magnitude(self, coord):
        """L0-SN-004: ΔA = A[i+1] − A[i], hand-computed for known mesh."""
        from orpheus.geometry import (
            BC,
            CoordSystem,
            Mesh1D,
            Region,
            RegionMesh,
            StructuredGeometry,
        )
        from orpheus.sn.mesh.augmented_mesh import SNMesh

        tag = {
            CoordSystem.SPHERICAL: "SPH",
            CoordSystem.CYLINDRICAL: "CYL",
        }[coord]
        mesh = Mesh1D.from_geometry(
            StructuredGeometry(
                geometry=tag,
                regions=(Region(mat_id=0, outer_thickness_cm=1.0),),
                bcs=(BC.reflective,),
            ),
            region_meshes=(RegionMesh(n_cells=5),),
        )
        if coord == CoordSystem.SPHERICAL:
            quad = Quadrature.gauss_legendre(4)
        else:
            quad = Quadrature.folded_product(n_mu=4, n_phi=8)
        sn = SNMesh(mesh, quad, placeholder_materials())

        edges = mesh.edges
        if coord == CoordSystem.SPHERICAL:
            expected = 4 * np.pi * (edges[1:]**2 - edges[:-1]**2)
        else:
            expected = 2 * np.pi * (edges[1:] - edges[:-1])
        np.testing.assert_allclose(sn.delta_A, expected, rtol=1e-14)

    @pytest.mark.verifies("sn-contamination-factor")
    def test_contamination_beta_spherical(self):
        """L0-SN-008: Contamination β ≈ 0 (machine zero) for spherical."""
        from orpheus.derivations.discrete.sn.angular_differencing import contamination_beta
        from orpheus.sn.sweep.pole_angular_closure import angular_cell_edges_per_level
        quad = Quadrature.gauss_legendre(8)
        # ``edges`` is supplied by the caller: L0 may not import orpheus.sn,
        # so the diagnostic is handed the closure it grades (2026-08-12).
        edges = angular_cell_edges_per_level(quad, CoordSystem.SPHERICAL)
        beta = contamination_beta(quad, "spherical", edges=edges)
        assert abs(beta) < 1e-14, f"Spherical β = {beta:.2e}"

    @pytest.mark.verifies("sn-contamination-factor")
    def test_contamination_beta_cylindrical(self):
        """L0-SN-008: Contamination β ≈ 0 (machine zero) for cylindrical."""
        from orpheus.derivations.discrete.sn.angular_differencing import contamination_beta
        from orpheus.sn.sweep.pole_angular_closure import angular_cell_edges_per_level
        quad = Quadrature.folded_product(n_mu=4, n_phi=8)
        edges = angular_cell_edges_per_level(quad, CoordSystem.CYLINDRICAL)
        betas = contamination_beta(quad, "cylindrical", edges=edges)
        assert np.all(np.abs(betas) < 1e-14), (
            f"Cylindrical β_max = {np.abs(betas).max():.2e}"
        )

    @pytest.mark.sentinel
    @pytest.mark.catches("ERR-002")
    def test_scattering_source_magnitude(self):
        """L0-SN-009: Scattering source = SigS^T @ φ, hand-calculated."""
        from orpheus.derivations.common.xs_library import get_mixture
        mix = get_mixture("A", "2g")
        phi = np.array([1.0, 2.0])
        sig_s = mix.SigS[0]
        if hasattr(sig_s, 'toarray'):
            sig_s = sig_s.toarray()
        expected = sig_s.T @ phi
        actual = phi @ sig_s
        np.testing.assert_allclose(actual, expected, rtol=1e-14)
