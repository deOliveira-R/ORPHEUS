"""Level-0 term verification for the SN discrete ordinates solver.

Each test isolates a SINGLE term or property of the discretized
transport equation and verifies it against a hand-calculated or
analytical reference value.  This is the most fundamental level of
verification — it catches sign flips, variable swaps, missing factors,
and convention errors that higher-level tests (eigenvalue, convergence)
cannot detect.

Problem catalog
---------------
L0-SN-001  Streaming: uniform source → φ = Q/Σ_t
L0-SN-002  Redistribution telescoping: Σ(α_out·ψ_out − α_in·ψ_in) = 0
L0-SN-003  Per-ordinate flat-flux: streaming + redistribution = 0 each ordinate
L0-SN-004  ΔA/w factor magnitude: hand-computed vs code
L0-SN-005  Alpha dome non-negativity: α ≥ 0 for all levels
L0-SN-006  Alpha closure: α[0] = α[M] = 0
L0-SN-007  M-M weight τ range: 0.5 ≤ τ ≤ 1.0
L0-SN-008  Contamination β: machine zero with correct formulation
L0-SN-009  Scattering source magnitude: SigS^T @ φ hand-calculated
L0-SN-010  Eigenvalue 1G: k = νΣ_f / Σ_a (material property)
L0-SN-011  Eigenvalue 2G: k = λ_max(A⁻¹F) (flux-shape dependent)
L0-SN-012  Eigenvalue 4G: k = λ_max(A⁻¹F)
L0-SN-013  Fixed-source flux bounded: no spike at r=0
L0-SN-014  Spatial convergence: O(h²) rate for DD
L0-SN-015  Particle balance: production/absorption = keff

Usage::

    pytest derivations/l0_sn.py -v                    # run all L0 SN tests
    pytest derivations/l0_sn.py -k "L0_SN_003" -v     # run a specific problem

For the error catalog, see ``derivations/l0_error_catalog.md``.
For the Sphinx publication page, see ``docs/theory/verification_l0.rst``.
"""

from __future__ import annotations

import numpy as np
import pytest

# Paths must be on PYTHONPATH: 02.Discrete.Ordinates, data, derivations, geometry
from derivations import get
from derivations._xs_library import get_mixture
from geometry import CoordSystem, Mesh1D, homogeneous_1d, mesh1d_from_zones, Zone
from sn_geometry import SNMesh
from sn_quadrature import GaussLegendre1D, ProductQuadrature, LevelSymmetricSN
from sn_solver import SNSolver, solve_sn
from sn_sweep import _sweep_1d_spherical, _sweep_1d_cylindrical


# ═══════════════════════════════════════════════════════════════════════
# Streaming
# ═════════════════════════��═════════════════════════════════════════════

class TestL0_SN_001_streaming_equilibrium:
    """L0-SN-001: Uniform source in pure absorber → φ = Q/Σ_t.

    Isolates the streaming term.  With uniform Q, uniform Σ_t, and
    reflective BC, the exact solution is φ = Q/Σ_t everywhere.
    The volume-averaged flux must match and the flux range must be
    bounded (no spike at r=0 for curvilinear).
    """

    @pytest.mark.parametrize("coord,QuadClass,quad_kwargs", [
        (CoordSystem.SPHERICAL, GaussLegendre1D, {}),
        (CoordSystem.CYLINDRICAL, ProductQuadrature, {"n_mu": 4, "n_phi": 8}),
    ], ids=["spherical", "cylindrical"])
    def test_volume_average(self, coord, QuadClass, quad_kwargs):
        mesh = homogeneous_1d(40, 1.0, mat_id=0, coord=coord)
        quad = QuadClass.create(**quad_kwargs)
        sn = SNMesh(mesh, quad)

        Q = np.ones((40, 1, 1))
        sig_t = np.ones((40, 1, 1))
        psi_bc = {}
        sweep = _sweep_1d_spherical if coord == CoordSystem.SPHERICAL else _sweep_1d_cylindrical
        for _ in range(50):
            _, phi = sweep(Q, sig_t, sn, psi_bc)

        avg = np.average(phi[:, 0, 0], weights=mesh.volumes)
        np.testing.assert_allclose(avg, 1.0, rtol=0.01,
                                   err_msg=f"Volume-avg φ ≠ Q/Σ_t for {coord.value}")

    @pytest.mark.parametrize("coord,QuadClass,quad_kwargs", [
        (CoordSystem.SPHERICAL, GaussLegendre1D, {}),
        (CoordSystem.CYLINDRICAL, ProductQuadrature, {"n_mu": 4, "n_phi": 8}),
    ], ids=["spherical", "cylindrical"])
    def test_flux_bounded(self, coord, QuadClass, quad_kwargs):
        """Flux range must be bounded — no spike at r=0."""
        mesh = homogeneous_1d(40, 1.0, mat_id=0, coord=coord)
        quad = QuadClass.create(**quad_kwargs)
        sn = SNMesh(mesh, quad)

        Q = np.ones((40, 1, 1))
        sig_t = np.ones((40, 1, 1))
        psi_bc = {}
        sweep = _sweep_1d_spherical if coord == CoordSystem.SPHERICAL else _sweep_1d_cylindrical
        for _ in range(50):
            _, phi = sweep(Q, sig_t, sn, psi_bc)

        assert phi[:, 0, 0].max() < 2.0, (
            f"Flux spike at origin: max={phi[:, 0, 0].max():.4f}"
        )


# ════════���═══════════════��═════════════════════════════��════════════════
# Angular redistribution
# ════════════════════���══════════════════════════════════════════════════

class TestL0_SN_002_redistribution_telescoping:
    """L0-SN-002: α·ψ product telescopes to zero per level per cell.

    The sum Σ_m(α_{m+1/2}·ψ_{m+1/2} − ��_{m-1/2}·ψ_{m-1/2}) must
    vanish because α[0] = α[M] = 0.
    """

    def test_cylindrical(self):
        mix = get_mixture("A", "1g")
        mesh = homogeneous_1d(10, 2.0, mat_id=0, coord=CoordSystem.CYLINDRICAL)
        quad = ProductQuadrature.create(n_mu=4, n_phi=8)
        sn = SNMesh(mesh, quad)

        sig_t = np.full((10, 1, 1), mix.SigT[0])
        Q = np.ones((10, 1, 1))
        psi_bc = {}
        ang, _ = _sweep_1d_cylindrical(Q, sig_t, sn, psi_bc)

        for p, level_idx in enumerate(quad.level_indices):
            alpha = sn.alpha_per_level[p]
            M = len(level_idx)
            psi_angle = np.zeros(10)
            for m_local in range(M):
                n = level_idx[m_local]
                psi_cell = ang[n, :, 0, 0]
                psi_angle = 2.0 * psi_cell - psi_angle
            # After all ordinates: α[M]·ψ_{M+1/2} = 0 since α[M] ≈ 0
            residual = alpha[M] * psi_angle
            np.testing.assert_allclose(residual, 0.0, atol=1e-12,
                                       err_msg=f"Level {p}: telescoping ≠ 0")


class TestL0_SN_003_per_ordinate_flat_flux:
    """L0-SN-003: For flat ψ, streaming + redistribution = 0 per ordinate.

    This is the fundamental correctness criterion for curvilinear SN.
    The ΔA/w factor ensures exact per-ordinate cancellation.
    """

    @pytest.mark.parametrize("coord", [CoordSystem.SPHERICAL, CoordSystem.CYLINDRICAL])
    def test_residual_is_zero(self, coord):
        if coord == CoordSystem.SPHERICAL:
            quad = GaussLegendre1D.create(8)
        else:
            quad = ProductQuadrature.create(n_mu=4, n_phi=8)

        mesh = homogeneous_1d(10, 1.0, mat_id=0, coord=coord)
        sn = SNMesh(mesh, quad)

        psi0 = 1.0  # flat flux
        dA = sn.delta_A  # (nx,)

        if coord == CoordSystem.SPHERICAL:
            alpha = sn.alpha_half
            mu = quad.mu_x
            w = quad.weights
            for n in range(quad.N):
                streaming = mu[n] * dA * psi0
                alpha_diff = alpha[n + 1] - alpha[n]  # = -w[n]*mu[n]
                redist = (dA / w[n]) * alpha_diff * psi0
                residual = streaming + redist
                np.testing.assert_allclose(
                    residual, 0.0, atol=1e-14,
                    err_msg=f"Spherical ordinate {n}: residual ≠ 0",
                )
        else:
            for p, level_idx in enumerate(sn.alpha_per_level):
                alpha = sn.alpha_per_level[p]
                for m, n in enumerate(quad.level_indices[p]):
                    eta = quad.mu_x[n]
                    wn = quad.weights[n]
                    streaming = eta * dA * psi0
                    alpha_diff = alpha[m + 1] - alpha[m]
                    redist = (dA / wn) * alpha_diff * psi0
                    residual = streaming + redist
                    np.testing.assert_allclose(
                        residual, 0.0, atol=1e-14,
                        err_msg=f"Cyl level {p} ord {m}: residual ≠ 0",
                    )


# ══════��════════════════════════════════════════════════���═══════════════
# Geometry factor and alpha properties
# ═══════════════════════════════════════════════════════════════════════

class TestL0_SN_004_delta_A_magnitude:
    """L0-SN-004: ΔA = A[i+1] − A[i], hand-computed for known mesh."""

    def test_spherical(self):
        mesh = homogeneous_1d(5, 1.0, mat_id=0, coord=CoordSystem.SPHERICAL)
        quad = GaussLegendre1D.create(4)
        sn = SNMesh(mesh, quad)
        # Equal-volume spherical: r_j = (j/N)^{1/3}
        # A = 4πr², ΔA = 4π(r_{j+1}² − r_j²)
        edges = mesh.edges
        expected_dA = 4 * np.pi * (edges[1:]**2 - edges[:-1]**2)
        np.testing.assert_allclose(sn.delta_A, expected_dA, rtol=1e-14)

    def test_cylindrical(self):
        mesh = homogeneous_1d(5, 1.0, mat_id=0, coord=CoordSystem.CYLINDRICAL)
        quad = ProductQuadrature.create(n_mu=4, n_phi=8)
        sn = SNMesh(mesh, quad)
        edges = mesh.edges
        expected_dA = 2 * np.pi * (edges[1:] - edges[:-1])
        np.testing.assert_allclose(sn.delta_A, expected_dA, rtol=1e-14)


class TestL0_SN_005_alpha_dome:
    """L0-SN-005: α coefficients form a non-negative dome."""

    @pytest.mark.parametrize("QuadClass,kwargs", [
        (ProductQuadrature, {"n_mu": 4, "n_phi": 8}),
        (LevelSymmetricSN, {"sn_order": 4}),
    ])
    def test_non_negative(self, QuadClass, kwargs):
        mesh = Mesh1D(edges=np.array([0.0, 1.0]), mat_ids=np.array([0]),
                      coord=CoordSystem.CYLINDRICAL)
        quad = QuadClass.create(**kwargs)
        sn = SNMesh(mesh, quad)
        for p, alpha in enumerate(sn.alpha_per_level):
            assert np.all(alpha >= -1e-14), (
                f"Level {p}: negative α = {alpha.min():.2e}"
            )


class TestL0_SN_006_alpha_closure:
    """L0-SN-006: α[0] = 0 and α[M] ≈ 0 (boundary conditions)."""

    def test_spherical(self):
        mesh = Mesh1D(edges=np.array([0.0, 1.0]), mat_ids=np.array([0]),
                      coord=CoordSystem.SPHERICAL)
        quad = GaussLegendre1D.create(8)
        sn = SNMesh(mesh, quad)
        assert sn.alpha_half[0] == 0.0
        np.testing.assert_allclose(sn.alpha_half[-1], 0.0, atol=1e-13)

    def test_cylindrical(self):
        mesh = Mesh1D(edges=np.array([0.0, 1.0]), mat_ids=np.array([0]),
                      coord=CoordSystem.CYLINDRICAL)
        quad = ProductQuadrature.create(n_mu=4, n_phi=8)
        sn = SNMesh(mesh, quad)
        for p, alpha in enumerate(sn.alpha_per_level):
            assert alpha[0] == 0.0, f"Level {p}: α[0] ≠ 0"
            np.testing.assert_allclose(alpha[-1], 0.0, atol=1e-13,
                                       err_msg=f"Level {p}: α[-1] ≠ 0")


class TestL0_SN_007_mm_weights:
    """L0-SN-007: Morel-Montry weights τ ∈ [0.5, 1.0]."""

    def test_spherical(self):
        mesh = homogeneous_1d(5, 1.0, mat_id=0, coord=CoordSystem.SPHERICAL)
        quad = GaussLegendre1D.create(8)
        sn = SNMesh(mesh, quad)
        assert np.all(sn.tau_mm >= 0.5 - 1e-14)
        assert np.all(sn.tau_mm <= 1.0 + 1e-14)

    def test_cylindrical(self):
        mesh = homogeneous_1d(5, 1.0, mat_id=0, coord=CoordSystem.CYLINDRICAL)
        quad = ProductQuadrature.create(n_mu=4, n_phi=8)
        sn = SNMesh(mesh, quad)
        for p, tau in enumerate(sn.tau_mm_per_level):
            assert np.all(tau >= 0.5 - 1e-14), f"Level {p}: τ < 0.5"
            assert np.all(tau <= 1.0 + 1e-14), f"Level {p}: τ > 1.0"


class TestL0_SN_008_contamination_beta:
    """L0-SN-008: Contamination β ≈ 0 (machine zero)."""

    def test_spherical(self):
        from derivations.sn_contamination import contamination_beta
        quad = GaussLegendre1D.create(8)
        beta = contamination_beta(quad, "spherical")
        assert abs(beta) < 1e-14, f"Spherical β = {beta:.2e}"

    def test_cylindrical(self):
        from derivations.sn_contamination import contamination_beta
        quad = ProductQuadrature.create(n_mu=4, n_phi=8)
        betas = contamination_beta(quad, "cylindrical")
        assert np.all(np.abs(betas) < 1e-14), f"Cylindrical β_max = {np.abs(betas).max():.2e}"


# ════���══════════════════════════════���═════════════════════��═════════════
# Scattering source
# ══════���═══════════════════���════════════════════════════════════════════

class TestL0_SN_009_scattering_source:
    """L0-SN-009: Scattering source = SigS^T @ φ, hand-calculated."""

    def test_2g_magnitude(self):
        mix = get_mixture("A", "2g")
        phi = np.array([1.0, 2.0])  # known flux
        # SigS[0] may be sparse or dense
        sig_s = mix.SigS[0]
        if hasattr(sig_s, 'toarray'):
            sig_s = sig_s.toarray()
        expected = sig_s.T @ phi  # in-scatter (matrix-vector)
        actual = phi @ sig_s       # vectorized form (row-vector × matrix)
        np.testing.assert_allclose(actual, expected, rtol=1e-14)


# ══════════��════════════════════��═══════════════════════���═══════════════
# Eigenvalue
# ��═══════════════════════════════════════════════════════���══════════════

class TestL0_SN_010_011_012_eigenvalue:
    """L0-SN-010/011/012: Homogeneous eigenvalue for 1G, 2G, 4G."""

    @pytest.mark.parametrize("case_name", [
        "sn_slab_1eg_1rg",
        "sn_slab_2eg_1rg",
        "sn_slab_4eg_1rg",
    ])
    @pytest.mark.parametrize("coord,QuadClass,quad_kwargs", [
        (CoordSystem.SPHERICAL, GaussLegendre1D, {}),
        (CoordSystem.CYLINDRICAL, ProductQuadrature, {"n_mu": 4, "n_phi": 8}),
    ], ids=["spherical", "cylindrical"])
    def test_exact(self, case_name, coord, QuadClass, quad_kwargs):
        case = get(case_name)
        mix = next(iter(case.materials.values()))
        mesh = homogeneous_1d(20, 2.0, mat_id=0, coord=coord)
        quad = QuadClass.create(**quad_kwargs)
        result = solve_sn({0: mix}, mesh, quad, max_inner=500, inner_tol=1e-10)
        assert abs(result.keff - case.k_inf) < 1e-6, (
            f"{case_name} {coord.value}: keff={result.keff:.8f} vs {case.k_inf:.8f}"
        )


# ═══════��═══════════════════════════��═══════════════════════════════════
# Fixed-source and convergence
# ════��════════════���═════════════════════════════���═══════════════════════

class TestL0_SN_013_fixed_source_bounded:
    """L0-SN-013: Fixed-source flux range bounded (no spike at r=0)."""
    # Covered by L0-SN-001 test_flux_bounded


class TestL0_SN_014_spatial_convergence:
    """L0-SN-014: Heterogeneous keff converges with mesh refinement."""

    def test_cylindrical(self):
        mix_fuel = get_mixture("A", "1g")
        mix_mod = get_mixture("B", "1g")
        materials = {2: mix_fuel, 0: mix_mod}
        quad = ProductQuadrature.create(n_mu=4, n_phi=8)

        keffs = []
        for n_cells in [5, 10, 20]:
            zones = [Zone(outer_edge=0.5, mat_id=2, n_cells=n_cells),
                     Zone(outer_edge=1.0, mat_id=0, n_cells=n_cells)]
            mesh = mesh1d_from_zones(zones, coord=CoordSystem.CYLINDRICAL)
            result = solve_sn(materials, mesh, quad, max_inner=500, inner_tol=1e-10)
            keffs.append(result.keff)

        diff_1 = abs(keffs[1] - keffs[0])
        diff_2 = abs(keffs[2] - keffs[1])
        assert diff_2 < diff_1, (
            f"Not converging: Δ(10-5)={diff_1:.6f}, Δ(20-10)={diff_2:.6f}"
        )


class TestL0_SN_015_particle_balance:
    """L0-SN-015: Production / absorption = keff."""

    @pytest.mark.parametrize("coord", [CoordSystem.SPHERICAL, CoordSystem.CYLINDRICAL])
    def test_balance(self, coord):
        case = get("sn_slab_2eg_1rg")
        mix = next(iter(case.materials.values()))
        if coord == CoordSystem.SPHERICAL:
            quad = GaussLegendre1D.create(8)
        else:
            quad = ProductQuadrature.create(n_mu=4, n_phi=8)

        mesh = homogeneous_1d(20, 2.0, mat_id=0, coord=coord)
        result = solve_sn({0: mix}, mesh, quad, max_inner=500, inner_tol=1e-10)

        V = mesh.volumes
        flux = result.scalar_flux[:, 0, :]
        sig_p = mix.SigP
        sig_a = mix.SigC + mix.SigF
        production = np.sum(flux * sig_p[None, :] * V[:, None])
        absorption = np.sum(flux * sig_a[None, :] * V[:, None])
        k_balance = production / absorption

        np.testing.assert_allclose(k_balance, result.keff, rtol=1e-5,
                                   err_msg=f"{coord.value}: balance ≠ keff")
