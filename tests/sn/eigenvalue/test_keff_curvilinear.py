"""Curvilinear (cylinder + sphere) SN k-eigenvalue verification.

Split from the legacy ``tests/sn/test_cylindrical.py`` /
``test_spherical.py`` (SN taxonomy reorg). This is the standalone home
the spec mandates for curvilinear k-eff — it previously existed ONLY
inside the mixed cyl/sph files and would have been lost in a naive
split. It carries the only ``l2`` markers in the 1-D suite.

The cylindrical eigenvalue tests are below; the spherical ones are
appended by the test_spherical split (phase 4). The per-sweep
regression guards moved to ``tests/sn/sweep/curvilinear/``; the CP
cross-checks to ``tests/sn/verification/analytical/``.
"""

import numpy as np
import pytest

from orpheus.derivations import get
from orpheus.derivations.common.xs_library import get_mixture
from orpheus.geometry import CoordSystem
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.geometry import SNMesh
from orpheus.sn.solver import SNSolver, solve_sn
from tests.sn._test_helpers import (
    curvilinear_homogeneous_mesh as _homogeneous_mesh,
    curvilinear_two_region_mesh as _two_region_mesh,
)

# Cylinder and sphere carry DIFFERENT equation-coverage lists (the
# legacy modules did too), so the verifies(...) sets are applied
# per-section, NOT as a single module-level pytestmark — a module mark
# would falsely stamp the cylinder labels onto the sphere tests (and
# vice versa). Each constant is the verbatim list from its legacy
# module so no verifies(...) edge is lost in the split.
_CYL_VERIFIES = pytest.mark.verifies(
    "transport-cylindrical",
    "alpha-cylindrical",
    "alpha-recursion",
    "wdd-closure",
    "wdd-face",
    "mm-weights",
    "multigroup",
    "reflective-bc",
    "one-group-kinf",
    "matrix-eigenvalue",
    "mg-balance",
    "balance-general",
)
_SPH_VERIFIES = pytest.mark.verifies(
    "transport-spherical",
    "alpha-recursion",
    "wdd-closure",
    "wdd-face",
    "multigroup",
    "reflective-bc",
    "one-group-kinf",
    "matrix-eigenvalue",
    "mg-balance",
    "balance-general",
)


# ═══════════════════════════════════════════════════════════════════════
# Cylindrical eigenvalue
# ═══════════════════════════════════════════════════════════════════════

@_CYL_VERIFIES
@pytest.mark.l1
@pytest.mark.parametrize("case_name", [
    "sn_slab_1eg_1rg",
    "sn_slab_2eg_1rg",
    "sn_slab_4eg_1rg",
])
@pytest.mark.parametrize("quad_factory", [
    lambda: Quadrature.product(n_mu=4, n_phi=8),
    lambda: Quadrature.level_symmetric(4),
], ids=["product", "level_sym"])
def test_homogeneous_exact(case_name, quad_factory):
    """Cylindrical SN on a homogeneous cylinder with reflective BC must
    match the analytical infinite-medium eigenvalue."""
    case = get(case_name)
    mix = next(iter(case.materials.values()))
    mesh = _homogeneous_mesh(20, 2.0, mat_id=0, coord=CoordSystem.CYLINDRICAL)
    quad = quad_factory()
    result = solve_sn({0: mix}, mesh, quad,
                      max_inner=500, inner_tol=1e-10)

    assert abs(result.keff - case.k_inf) < 1e-6, (
        f"keff={result.keff:.8f} vs analytical={case.k_inf:.8f}"
    )


@_CYL_VERIFIES
@pytest.mark.l1
@pytest.mark.parametrize("quad_factory", [
    lambda: Quadrature.product(n_mu=4, n_phi=8),
    lambda: Quadrature.level_symmetric(4),
], ids=["product", "level_sym"])
def test_particle_balance(quad_factory):
    """For reflective BCs (no leakage), production / absorption = keff."""
    case = get("sn_slab_2eg_1rg")
    mix = next(iter(case.materials.values()))
    mesh = _homogeneous_mesh(20, 2.0, mat_id=0, coord=CoordSystem.CYLINDRICAL)
    quad = quad_factory()
    result = solve_sn({0: mix}, mesh, quad,
                      max_inner=500, inner_tol=1e-10)

    V = mesh.volumes
    flux = result.scalar_flux.values[:, :, 0].T  # PR-INDEX-5
    sig_p = mix.SigP
    sig_a = mix.SigC + mix.SigF

    production = np.sum(flux * sig_p[None, :] * V[:, None])
    absorption = np.sum(flux * sig_a[None, :] * V[:, None])

    k_balance = production / absorption
    np.testing.assert_allclose(
        k_balance, result.keff, rtol=1e-5,
        err_msg=f"Particle balance: prod/abs={k_balance:.8f} ≠ keff={result.keff:.8f}",
    )


@_CYL_VERIFIES
@pytest.mark.l1
def test_flux_non_negative():
    mix = get_mixture("A", "1g")
    mesh = _homogeneous_mesh(10, 2.0, mat_id=0, coord=CoordSystem.CYLINDRICAL)
    quad = Quadrature.product(n_mu=4, n_phi=8)
    result = solve_sn({0: mix}, mesh, quad, max_inner=500, inner_tol=1e-10)

    assert np.all(result.scalar_flux.values >= 0), (
        f"Negative flux: min={result.scalar_flux.values.min():.4e}"
    )


@_CYL_VERIFIES
@pytest.mark.l2
class TestCylinderMultiGroupMultiRegion:
    """Multi-group heterogeneous cylindrical eigenvalue integration.

    The minimum problems that catch normalization, scattering, and
    eigenvector-distortion bugs simultaneously — invisible in 1G keff.
    """

    def test_2g_heterogeneous_fuel_moderator(self):
        """2G fuel+moderator cylinder — minimum problem catching
        normalization, scattering, and redistribution bugs at once."""
        fuel = get_mixture("A", "2g")
        mod = get_mixture("B", "2g")
        materials = {2: fuel, 0: mod}

        mesh = _two_region_mesh(
            outers=(0.5, 1.0), mat_ids=(2, 0), n_cells=(10, 10),
            coord=CoordSystem.CYLINDRICAL,
        )
        quad = Quadrature.product(n_mu=4, n_phi=8)
        result = solve_sn(materials, mesh, quad,
                          max_inner=500, inner_tol=1e-10)

        assert np.isfinite(result.keff), "keff is NaN/Inf"
        assert result.keff > 0, f"keff is non-positive: {result.keff}"
        assert np.all(np.isfinite(result.scalar_flux.values)), "Non-finite flux"
        assert 0.5 < result.keff < 3.0, f"keff={result.keff:.4f} out of physical range"

    def test_2g_heterogeneous_product_different_resolutions(self):
        """Product quadrature at two resolutions must give close keff."""
        fuel = get_mixture("A", "2g")
        mod = get_mixture("B", "2g")
        materials = {2: fuel, 0: mod}

        keffs = {}
        for label, quad in [
            ("4×8", Quadrature.product(n_mu=4, n_phi=8)),
            ("8×8", Quadrature.product(n_mu=8, n_phi=8)),
        ]:
            mesh = _two_region_mesh(
                outers=(0.5, 1.0), mat_ids=(2, 0), n_cells=(10, 10),
                coord=CoordSystem.CYLINDRICAL,
            )
            result = solve_sn(materials, mesh, quad,
                              max_inner=500, inner_tol=1e-10)
            keffs[label] = result.keff

        assert abs(keffs["4×8"] - keffs["8×8"]) < 0.05, (
            f"Product resolutions disagree: "
            f"4×8={keffs['4×8']:.6f}, 8×8={keffs['8×8']:.6f}"
        )

    def test_4g_homogeneous_scattering_convergence(self):
        """4G homogeneous with strong scattering must converge.

        4-group has the richest scattering matrix (10 nonzero entries)
        and is the most sensitive to iteration divergence.
        """
        mix = get_mixture("A", "4g")
        mesh = _homogeneous_mesh(20, 2.0, mat_id=0, coord=CoordSystem.CYLINDRICAL)
        quad = Quadrature.product(n_mu=4, n_phi=8)
        sn_mesh = SNMesh(mesh, quad, {0: mix})
        solver = SNSolver(sn_mesh, max_inner=500, inner_tol=1e-10)

        phi = solver.initial_flux_distribution()
        keff = 1.0
        for _ in range(5):
            fs = solver.compute_fission_source(phi, keff)
            phi = solver.solve_fixed_source(fs, phi)
            keff = solver.compute_keff(phi)

        assert np.all(np.isfinite(phi)), "4G scattering iteration diverged"
        assert phi.max() < 1e10, f"4G flux blew up to {phi.max():.2e}"

    def test_multigroup_eigenvector_not_flat(self):
        """For multi-group heterogeneous, the flux spectrum must vary
        between fuel and moderator — a flat spectrum indicates the
        multi-group coupling is broken.
        """
        fuel = get_mixture("A", "2g")
        mod = get_mixture("B", "2g")
        materials = {2: fuel, 0: mod}

        mesh = _two_region_mesh(
            outers=(0.5, 1.0), mat_ids=(2, 0), n_cells=(10, 10),
            coord=CoordSystem.CYLINDRICAL,
        )
        quad = Quadrature.product(n_mu=4, n_phi=8)
        result = solve_sn(materials, mesh, quad,
                          max_inner=500, inner_tol=1e-10)

        flux = result.scalar_flux.values[:, :, 0].T  # PR-INDEX-5  # (nx, ng)
        V = mesh.volumes
        mat_ids = mesh.mat_ids

        fuel_flux = np.average(flux[mat_ids == 2], axis=0, weights=V[mat_ids == 2])
        mod_flux = np.average(flux[mat_ids == 0], axis=0, weights=V[mat_ids == 0])

        fuel_ratio = fuel_flux[0] / fuel_flux[1]
        mod_ratio = mod_flux[0] / mod_flux[1]

        assert abs(fuel_ratio - mod_ratio) > 0.01, (
            f"Flux spectrum identical in fuel and moderator — "
            f"multi-group coupling may be broken: "
            f"fuel ratio={fuel_ratio:.4f}, mod ratio={mod_ratio:.4f}"
        )

    def test_particle_balance_heterogeneous(self):
        """Particle balance must hold for heterogeneous multi-region."""
        fuel = get_mixture("A", "2g")
        mod = get_mixture("B", "2g")
        materials = {2: fuel, 0: mod}

        mesh = _two_region_mesh(
            outers=(0.5, 1.0), mat_ids=(2, 0), n_cells=(10, 10),
            coord=CoordSystem.CYLINDRICAL,
        )
        quad = Quadrature.product(n_mu=4, n_phi=8)
        sn_mesh = SNMesh(mesh, quad, materials)
        solver = SNSolver(sn_mesh, max_inner=500, inner_tol=1e-10)

        phi = solver.initial_flux_distribution()
        keff = 1.0
        for _ in range(100):
            fs = solver.compute_fission_source(phi, keff)
            phi = solver.solve_fixed_source(fs, phi)
            keff = solver.compute_keff(phi)

        vol = solver.volume[None, :, :]
        production = np.sum(solver.mat_xs.fission_production * phi * vol)
        absorption = np.sum(solver.mat_xs.absorption_cross_section * phi * vol)
        k_balance = production / absorption

        np.testing.assert_allclose(
            k_balance, keff, rtol=1e-4,
            err_msg=f"Heterogeneous particle balance: {k_balance:.6f} ≠ {keff:.6f}",
        )

    def test_heterogeneous_1g_spatial_convergence(self):
        """keff must converge monotonically with mesh refinement."""
        mix_fuel = get_mixture("A", "1g")
        mix_mod = get_mixture("B", "1g")
        materials = {2: mix_fuel, 0: mix_mod}
        quad = Quadrature.product(n_mu=4, n_phi=8)

        keffs = []
        for n_cells in [5, 10, 20]:
            mesh = _two_region_mesh(
                outers=(0.5, 1.0), mat_ids=(2, 0), n_cells=(n_cells, n_cells),
                coord=CoordSystem.CYLINDRICAL,
            )
            result = solve_sn(materials, mesh, quad,
                              max_inner=500, inner_tol=1e-10)
            keffs.append(result.keff)

        diff_1 = abs(keffs[1] - keffs[0])
        diff_2 = abs(keffs[2] - keffs[1])
        assert diff_2 < diff_1, (
            f"keff not converging: Δ(10−5)={diff_1:.6f}, Δ(20−10)={diff_2:.6f}, "
            f"keffs={[f'{k:.6f}' for k in keffs]}"
        )


# ═══════════════════════════════════════════════════════════════════════
# Spherical eigenvalue (split from the legacy test_spherical.py)
# ═══════════════════════════════════════════════════════════════════════
#
# Grouped in classes (collision-free with the cylinder module functions)
# and decorated with _SPH_VERIFIES — the sphere's own equation-coverage
# list, distinct from the cylinder's.

@_SPH_VERIFIES
@pytest.mark.l1
class TestSphereEigenvalue:
    """Spherical homogeneous-exact / balance / convergence eigenvalue gates."""

    @pytest.mark.parametrize("case_name", [
        "sn_slab_1eg_1rg",
        "sn_slab_2eg_1rg",
        "sn_slab_4eg_1rg",
    ])
    def test_homogeneous_exact(self, case_name):
        """Spherical SN on a homogeneous sphere with reflective BC must
        match the analytical infinite-medium eigenvalue."""
        case = get(case_name)
        mix = next(iter(case.materials.values()))
        materials = {0: mix}
        mesh = _homogeneous_mesh(20, 2.0, mat_id=0, coord=CoordSystem.SPHERICAL)
        quad = Quadrature.gauss_legendre(8)
        result = solve_sn(materials, mesh, quad,
                          max_inner=500, inner_tol=1e-10)

        # Spherical DD has larger discretization error than Cartesian
        # due to angular redistribution coupling. 1G is exact (keff
        # independent of flux shape); multi-group has ~1% error on S8/20-cell.
        tol = 1e-6 if case.n_groups == 1 else 0.02
        assert abs(result.keff - case.k_inf) < tol, (
            f"keff={result.keff:.8f} vs analytical={case.k_inf:.8f} "
            f"err={abs(result.keff - case.k_inf):.2e}"
        )

    def test_particle_balance(self):
        """For reflective BCs (no leakage), production / absorption = keff."""
        case = get("sn_slab_2eg_1rg")
        mix = next(iter(case.materials.values()))
        materials = {0: mix}
        mesh = _homogeneous_mesh(20, 2.0, mat_id=0, coord=CoordSystem.SPHERICAL)
        quad = Quadrature.gauss_legendre(8)
        result = solve_sn(materials, mesh, quad,
                          max_inner=500, inner_tol=1e-10)

        V = mesh.volumes
        flux = result.scalar_flux.values[:, :, 0].T  # PR-INDEX-5  # (nx, ng)
        sig_p = mix.SigP
        sig_a = mix.SigC + mix.SigF

        production = np.sum(flux * sig_p[None, :] * V[:, None])
        absorption = np.sum(flux * sig_a[None, :] * V[:, None])

        k_balance = production / absorption
        np.testing.assert_allclose(
            k_balance, result.keff, rtol=1e-5,
            err_msg=f"Particle balance: prod/abs={k_balance:.8f} ≠ keff={result.keff:.8f}",
        )

    @pytest.mark.slow
    def test_spatial_convergence(self):
        """Diamond-difference on spherical mesh must show O(h²) convergence."""
        fuel = get_mixture("A", "1g")
        mod = get_mixture("B", "1g")
        materials = {2: fuel, 0: mod}

        keffs = []
        drs = []
        for n_per in [5, 10, 20, 40]:
            mesh = _two_region_mesh(
                outers=(0.5, 1.0), mat_ids=(2, 0), n_cells=(n_per, n_per),
                coord=CoordSystem.SPHERICAL,
            )
            quad = Quadrature.gauss_legendre(16)
            result = solve_sn(
                materials, mesh, quad,
                max_outer=300, max_inner=500, inner_tol=1e-10,
            )
            keffs.append(result.keff)
            drs.append(0.5 / n_per)

        k_ref = keffs[-1] + (keffs[-1] - keffs[-2]) / 3.0

        orders = []
        for i in range(1, len(keffs)):
            err_prev = abs(keffs[i - 1] - k_ref)
            err_curr = abs(keffs[i] - k_ref)
            if err_prev > 0 and err_curr > 0:
                orders.append(
                    np.log(err_prev / err_curr)
                    / np.log(drs[i - 1] / drs[i])
                )

        assert orders[-1] > 1.5, (
            f"Expected O(h²) convergence, got order {orders[-1]:.2f}"
        )

    def test_flux_non_negative(self):
        """Converged scalar flux must be non-negative everywhere."""
        mix = get_mixture("A", "1g")
        mesh = _homogeneous_mesh(10, 2.0, mat_id=0, coord=CoordSystem.SPHERICAL)
        quad = Quadrature.gauss_legendre(8)
        result = solve_sn({0: mix}, mesh, quad, max_inner=500, inner_tol=1e-10)

        assert np.all(result.scalar_flux.values >= 0), (
            f"Negative flux: min={result.scalar_flux.values.min():.4e}"
        )


@_SPH_VERIFIES
@pytest.mark.l2
class TestMultiGroupMultiRegionSpherical:
    """Multi-group heterogeneous spherical eigenvalue integration.

    Spherical-specific: angular redistribution + multi-group scattering
    is the combination most likely to expose normalization / coupling bugs.
    """

    def test_2g_heterogeneous_converges(self):
        """2G fuel+moderator sphere must converge to finite keff."""
        fuel = get_mixture("A", "2g")
        mod = get_mixture("B", "2g")
        materials = {2: fuel, 0: mod}

        mesh = _two_region_mesh(
            outers=(0.5, 1.0), mat_ids=(2, 0), n_cells=(10, 10),
            coord=CoordSystem.SPHERICAL,
        )
        quad = Quadrature.gauss_legendre(8)
        result = solve_sn(materials, mesh, quad,
                          max_inner=500, inner_tol=1e-10)

        assert np.isfinite(result.keff), "keff is NaN/Inf"
        assert 0.1 < result.keff < 3.0, f"keff={result.keff:.4f} out of range"
        assert np.all(np.isfinite(result.scalar_flux.values)), "Non-finite flux"

    def test_4g_scattering_convergence(self):
        """4G homogeneous must converge (richest scattering matrix)."""
        mix = get_mixture("A", "4g")
        mesh = _homogeneous_mesh(20, 2.0, mat_id=0, coord=CoordSystem.SPHERICAL)
        quad = Quadrature.gauss_legendre(8)
        sn_mesh = SNMesh(mesh, quad, {0: mix})
        solver = SNSolver(sn_mesh, max_inner=500, inner_tol=1e-10)

        phi = solver.initial_flux_distribution()
        keff = 1.0
        for _ in range(5):
            fs = solver.compute_fission_source(phi, keff)
            phi = solver.solve_fixed_source(fs, phi)
            keff = solver.compute_keff(phi)

        assert np.all(np.isfinite(phi)), "4G scattering iteration diverged"
        assert phi.max() < 1e10, f"4G flux blew up to {phi.max():.2e}"

    def test_multigroup_eigenvector_not_flat(self):
        """Flux spectrum must differ between fuel and moderator."""
        fuel = get_mixture("A", "2g")
        mod = get_mixture("B", "2g")
        materials = {2: fuel, 0: mod}

        mesh = _two_region_mesh(
            outers=(0.5, 1.0), mat_ids=(2, 0), n_cells=(10, 10),
            coord=CoordSystem.SPHERICAL,
        )
        quad = Quadrature.gauss_legendre(8)
        result = solve_sn(materials, mesh, quad,
                          max_inner=500, inner_tol=1e-10)

        flux = result.scalar_flux.values[:, :, 0].T  # PR-INDEX-5
        V = mesh.volumes
        mat_ids = mesh.mat_ids

        fuel_flux = np.average(flux[mat_ids == 2], axis=0, weights=V[mat_ids == 2])
        mod_flux = np.average(flux[mat_ids == 0], axis=0, weights=V[mat_ids == 0])

        fuel_ratio = fuel_flux[0] / fuel_flux[1]
        mod_ratio = mod_flux[0] / mod_flux[1]

        assert abs(fuel_ratio - mod_ratio) > 0.01, (
            f"Spectrum identical in fuel/mod — coupling broken: "
            f"fuel={fuel_ratio:.4f}, mod={mod_ratio:.4f}"
        )

    def test_particle_balance_heterogeneous(self):
        """Particle balance on 2G heterogeneous sphere."""
        fuel = get_mixture("A", "2g")
        mod = get_mixture("B", "2g")
        materials = {2: fuel, 0: mod}

        mesh = _two_region_mesh(
            outers=(0.5, 1.0), mat_ids=(2, 0), n_cells=(10, 10),
            coord=CoordSystem.SPHERICAL,
        )
        quad = Quadrature.gauss_legendre(8)
        sn_mesh = SNMesh(mesh, quad, materials)
        solver = SNSolver(sn_mesh, max_inner=500, inner_tol=1e-10)

        phi = solver.initial_flux_distribution()
        keff = 1.0
        for _ in range(100):
            fs = solver.compute_fission_source(phi, keff)
            phi = solver.solve_fixed_source(fs, phi)
            keff = solver.compute_keff(phi)

        vol = solver.volume[None, :, :]  # PR-INDEX-5
        production = np.sum(solver.mat_xs.fission_production * phi * vol)
        absorption = np.sum(solver.mat_xs.absorption_cross_section * phi * vol)
        k_balance = production / absorption

        np.testing.assert_allclose(
            k_balance, keff, rtol=1e-4,
            err_msg=f"Heterogeneous balance: {k_balance:.6f} ≠ {keff:.6f}",
        )

    def test_fixed_source_flux_bounded(self):
        """Fixed-source flux range must be bounded near r=0.

        Without the ΔA/w geometry factor, the flux spikes to ~5x at
        the origin.  With the fix, the range should be bounded.
        """
        from orpheus.sn.sweep import transport_sweep
        from orpheus.transport.sources import PerOrdinateSource
        from tests.sn._test_helpers import placeholder_materials

        mesh = _homogeneous_mesh(40, 1.0, mat_id=0, coord=CoordSystem.SPHERICAL)
        quad = Quadrature.gauss_legendre(8)
        sn_mesh = SNMesh(mesh, quad, placeholder_materials())

        sig_t = np.ones((1, 40, 1))             # (ng, nx, ny)
        Q_iso = np.ones((1, 40, 1))             # (ng, nx, ny)
        source = PerOrdinateSource.from_isotropic(Q_iso, sn_mesh)
        boundary_flux = sn_mesh.zeros_boundary_flux()
        phi = None
        for _ in range(50):
            _, phi = transport_sweep(source, sig_t, sn_mesh, boundary_flux)

        phi_avg = np.average(phi[0, :, 0], weights=mesh.volumes)
        np.testing.assert_allclose(phi_avg, 1.0, rtol=0.01,
                                   err_msg="Volume-avg φ ≠ Q/Σ_t")
        assert phi[0, :, 0].max() < 2.0, (
            f"Flux spike at origin: max={phi[0, :, 0].max():.4f}"
        )

    def test_heterogeneous_1g_spatial_convergence(self):
        """keff must converge with mesh refinement for fuel+moderator."""
        mix_fuel = get_mixture("A", "1g")
        mix_mod = get_mixture("B", "1g")
        materials = {2: mix_fuel, 0: mix_mod}
        quad = Quadrature.gauss_legendre(8)

        keffs = []
        for n_cells in [5, 10, 20]:
            mesh = _two_region_mesh(
                outers=(0.5, 1.0), mat_ids=(2, 0), n_cells=(n_cells, n_cells),
                coord=CoordSystem.SPHERICAL,
            )
            result = solve_sn(materials, mesh, quad,
                              max_inner=500, inner_tol=1e-10)
            keffs.append(result.keff)

        diff_1 = abs(keffs[1] - keffs[0])
        diff_2 = abs(keffs[2] - keffs[1])
        assert diff_2 < diff_1, (
            f"keff not converging: Δ(10−5)={diff_1:.6f}, Δ(20−10)={diff_2:.6f}, "
            f"keffs={[f'{k:.6f}' for k in keffs]}"
        )
