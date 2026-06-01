"""2-D Cartesian SN k-eigenvalue verification.

Split from the legacy ``tests/sn/test_solver_components.py`` (SN
taxonomy reorg). This file is the 2-D eigenvalue home — it fills the
gap the reorg spec flagged ("NO 2-D k-eigenvalue test exists" as a
standalone). It carries the homogeneous-exact k_inf checks, the
multigroup-eigenvector recovery, the quadrature-normalisation
eigenvalue invariants (BiCGSTAB GL-vs-Lebedev, ERR-004), the Pn
eigenvalue effects, and the SI-vs-Krylov agreement.

R-1 Step E note: every test that drives the 2-D Cartesian
*source-iteration* inner solver currently raises NotImplementedError
(the 2-D SI carve is deferred; the B1'' face block is 1-D-only). Those
tests are marked ``xfail(raises=NotImplementedError, strict=False)`` —
they preserve the coverage and will xpass when the carve lands. Tests
that drive only the *Krylov* 2-D path (carved separately) pass today.
"""

from pathlib import Path

import numpy as np
import pytest

from orpheus.derivations import get
from orpheus.derivations.common.xs_library import get_mixture
from orpheus.geometry import Mesh1D, Mesh2D
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.geometry import SNMesh
from orpheus.sn.solver import SNSolver, solve_sn

pytestmark = pytest.mark.l0  # SN 2-D eigenvalue component checks

# R-1 Step E: the 2-D Cartesian source-iteration inner solver is
# deferred (sn_mesh.reduced is None ⇒ NotImplementedError). Shared
# xfail mark for the tests that exercise that path. strict=False so the
# tests xpass — alerting us — once the 2-D SI carve lands.
_DEFERRED_2D_SI = pytest.mark.xfail(
    reason="R-1 Step E: 2-D Cartesian source-iteration carve deferred "
    "(NotImplementedError). The B1'' face block is 1-D-only; 2-D needs "
    "the 4-face (xmin,xmax,ymin,ymax) layout. Coverage preserved; xpass "
    "when the carve lands.",
    raises=NotImplementedError,
    strict=False,
)


def _uniform_2d(nx, ny, delta, mat_map):
    """Build a uniform Mesh2D."""
    return Mesh2D(
        edges_x=np.linspace(0, nx * delta, nx + 1),
        edges_y=np.linspace(0, ny * delta, ny + 1),
        mat_map=np.asarray(mat_map, dtype=int),
    )


class TestHomogeneousExact:
    """2D SN on homogeneous infinite medium must match analytical k_inf.

    Uses small 2×2 mesh with reflective BCs. k_inf is the max eigenvalue
    of A⁻¹F where A = Σ_t - Σ_s^T, F = χ⊗(νΣ_f).
    """

    @_DEFERRED_2D_SI
    @pytest.mark.parametrize("ng_key,label", [("1g", "1G"), ("2g", "2G"), ("4g", "4G")])
    def test_homogeneous_exact(self, ng_key, label):
        case = get(f"sn_slab_{ng_key[0]}eg_1rg")
        mix = next(iter(case.materials.values()))
        materials = {0: mix}

        mesh = _uniform_2d(2, 2, 0.5, np.zeros((2, 2), dtype=int))
        quad = Quadrature.lebedev(order=17)
        solver = SNSolver(SNMesh(mesh, quad, materials), max_inner=500, inner_tol=1e-10)

        phi = solver.initial_flux_distribution()
        keff = 1.0
        for _ in range(100):
            fs = solver.compute_fission_source(phi, keff)
            phi = solver.solve_fixed_source(fs, phi)
            keff = solver.compute_keff(phi)
            phi = phi / np.linalg.norm(phi)

        assert abs(keff - case.k_inf) < 1e-8, (
            f"{label}: keff={keff:.10f} vs analytical={case.k_inf:.10f}"
        )


class TestMultiGroupEigenvector:
    """The converged flux group ratio must match the analytical eigenvector.

    1-group tests hide bugs because k = νΣf/Σa is flux-shape independent.
    Multi-group problems have a specific eigenvector (group ratio).
    """

    @_DEFERRED_2D_SI
    def test_2g_eigenvector(self):
        case = get("sn_slab_2eg_1rg")
        mix = next(iter(case.materials.values()))

        # Analytical eigenvector from (Σ_t - Σ_s^T)^{-1} · χ⊗(νΣf)
        sig_s = mix.SigS[0].toarray()
        A = np.diag(mix.SigT) - sig_s.T
        F = np.outer(mix.chi, mix.SigP)
        _, vecs = np.linalg.eig(np.linalg.solve(A, F))
        idx = np.argmax(np.real(np.linalg.eigvals(np.linalg.solve(A, F))))
        phi_expected = np.real(vecs[:, idx])
        phi_expected /= phi_expected.sum()

        materials = {0: mix}
        mesh = _uniform_2d(2, 2, 0.5, np.zeros((2, 2), dtype=int))
        quad = Quadrature.lebedev(order=17)
        solver = SNSolver(SNMesh(mesh, quad, materials), max_inner=500, inner_tol=1e-10)

        phi = solver.initial_flux_distribution()
        keff = 1.0
        for _ in range(100):
            fs = solver.compute_fission_source(phi, keff)
            phi = solver.solve_fixed_source(fs, phi)
            keff = solver.compute_keff(phi)
            phi = phi / np.linalg.norm(phi)

        phi_cell = phi[:, 0, 0]
        phi_ratio = phi_cell / phi_cell.sum()

        np.testing.assert_allclose(phi_ratio, phi_expected, rtol=1e-6,
                                   err_msg="Converged group ratio ≠ analytical eigenvector")


class TestBicgstabNormalization:
    """BiCGSTAB must give the same keff regardless of quadrature type.

    build_rhs hardcoded 4π for the angular normalization, but GL weights
    sum to 2, not 4π. The normalization must use sum(weights),
    quadrature-dependent. (ERR-004.)
    """

    @pytest.mark.catches("ERR-004")
    def test_1d_gl_homogeneous_exact(self):
        """BiCGSTAB with GL quadrature on 1D slab must match analytical k_inf."""
        case = get("sn_slab_2eg_1rg")
        mix = next(iter(case.materials.values()))

        mesh = Mesh1D(edges=np.linspace(0, 2, 5), mat_ids=np.zeros(4, dtype=int))
        gl = Quadrature.gauss_legendre(8)
        solver = SNSolver(SNMesh(mesh, gl, {0: mix}), inner_solver="krylov", max_inner=2000, inner_tol=1e-6)

        phi = solver.initial_flux_distribution()
        keff = 1.0
        for _ in range(50):
            fs = solver.compute_fission_source(phi, keff)
            phi = solver.solve_fixed_source(fs, phi)
            keff = solver.compute_keff(phi)
            phi /= np.linalg.norm(phi)

        assert abs(keff - case.k_inf) < 1e-4, (
            f"1D GL BiCGSTAB keff={keff:.8f} vs analytical={case.k_inf:.8f}"
        )

    def test_2d_lebedev_homogeneous_exact(self):
        """BiCGSTAB with Lebedev quadrature on 2D mesh must match analytical k_inf."""
        case = get("sn_slab_2eg_1rg")
        mix = next(iter(case.materials.values()))

        mesh = _uniform_2d(2, 2, 0.5, np.zeros((2, 2), dtype=int))
        quad = Quadrature.lebedev(order=17)
        solver = SNSolver(SNMesh(mesh, quad, {0: mix}), inner_solver="krylov", max_inner=2000, inner_tol=1e-6)

        phi = solver.initial_flux_distribution()
        keff = 1.0
        for _ in range(50):
            fs = solver.compute_fission_source(phi, keff)
            phi = solver.solve_fixed_source(fs, phi)
            keff = solver.compute_keff(phi)
            phi /= np.linalg.norm(phi)

        assert abs(keff - case.k_inf) < 1e-4, (
            f"2D Lebedev BiCGSTAB keff={keff:.8f} vs analytical={case.k_inf:.8f}"
        )

    def test_gl_and_lebedev_agree(self):
        """BiCGSTAB keff must not depend on which quadrature is used."""
        case = get("sn_slab_2eg_1rg")
        mix = next(iter(case.materials.values()))

        results = {}
        for label, mesh, quad in [
            ("GL", Mesh1D(edges=np.linspace(0, 2, 5), mat_ids=np.zeros(4, dtype=int)),
             Quadrature.gauss_legendre(8)),
            ("Lebedev", _uniform_2d(2, 2, 0.5, np.zeros((2, 2), dtype=int)),
             Quadrature.lebedev(order=17)),
        ]:
            solver = SNSolver(SNMesh(mesh, quad, {0: mix}), inner_solver="krylov", max_inner=2000, inner_tol=1e-6)
            phi = solver.initial_flux_distribution()
            keff = 1.0
            for _ in range(50):
                fs = solver.compute_fission_source(phi, keff)
                phi = solver.solve_fixed_source(fs, phi)
                keff = solver.compute_keff(phi)
                phi /= np.linalg.norm(phi)
            results[label] = keff

        assert abs(results["GL"] - results["Lebedev"]) < 1e-3, (
            f"GL keff={results['GL']:.6f} vs Lebedev keff={results['Lebedev']:.6f}"
        )


@pytest.mark.verifies("pn-scatter")
class TestAnisotropicScatteringKeff:
    """Pn-scattering EIGENVALUE effects on a 2-D mesh.

    Split from the legacy ``TestAnisotropicScattering`` — these are the
    members that drive a full 2-D eigenvalue solve (not the harmonic
    orthogonality / clamp checks, which stayed with the operator-level
    component tests in operators/test_solver_components.py). Both
    exercise the Legendre-expanded scattering source from
    Eq. :label:`pn-scatter`.
    """

    @_DEFERRED_2D_SI
    def test_p0_gives_identical_keff(self):
        """scattering_order=0 must give the exact same keff as the default."""
        case = get("sn_slab_2eg_1rg")
        mix = next(iter(case.materials.values()))
        mesh = _uniform_2d(2, 2, 0.5, np.zeros((2, 2), dtype=int))
        quad = Quadrature.lebedev(order=17)

        solver_default = SNSolver(SNMesh(mesh, quad, {0: mix}), max_inner=500, inner_tol=1e-10)
        phi = solver_default.initial_flux_distribution()
        keff = 1.0
        for _ in range(50):
            fs = solver_default.compute_fission_source(phi, keff)
            phi = solver_default.solve_fixed_source(fs, phi)
            keff = solver_default.compute_keff(phi)
            phi /= np.linalg.norm(phi)
        keff_p0 = keff

        solver_explicit = SNSolver(SNMesh(mesh, quad, {0: mix}), scattering_order=0, max_inner=500, inner_tol=1e-10)
        phi = solver_explicit.initial_flux_distribution()
        keff = 1.0
        for _ in range(50):
            fs = solver_explicit.compute_fission_source(phi, keff)
            phi = solver_explicit.solve_fixed_source(fs, phi)
            keff = solver_explicit.compute_keff(phi)
            phi /= np.linalg.norm(phi)
        keff_explicit = keff

        assert abs(keff_p0 - keff_explicit) < 1e-14, (
            f"P0 default {keff_p0:.10f} != explicit P0 {keff_explicit:.10f}"
        )

    @_DEFERRED_2D_SI
    def test_p1_changes_heterogeneous_keff(self):
        """P1 scattering must produce a different keff than P0 on a
        heterogeneous problem where anisotropy matters at interfaces."""
        fuel = get_mixture("A", "2g")
        mod = get_mixture("B", "2g")  # B has mu_bar=0.6, strongly anisotropic
        materials = {2: fuel, 0: mod}

        mat = np.zeros((6, 2), dtype=int)
        mat[:3, :] = 2
        mesh = _uniform_2d(6, 2, 0.2, mat)
        quad = Quadrature.lebedev(order=17)

        keffs = {}
        for L in [0, 1]:
            solver = SNSolver(SNMesh(mesh, quad, materials), scattering_order=L, max_inner=500, inner_tol=1e-10)
            phi = solver.initial_flux_distribution()
            keff = 1.0
            for _ in range(50):
                fs = solver.compute_fission_source(phi, keff)
                phi = solver.solve_fixed_source(fs, phi)
                keff = solver.compute_keff(phi)
                phi /= np.linalg.norm(phi)
            keffs[L] = keff

        assert abs(keffs[0] - keffs[1]) > 1e-4, (
            f"P0 keff={keffs[0]:.6f} and P1 keff={keffs[1]:.6f} should differ"
        )


@pytest.mark.verifies("pn-scatter")
class TestBicgstabPnScattering:
    """BiCGSTAB path must handle Pn scattering consistently with source iteration.

    Verifies the Legendre-expanded scattering source (Eq. :label:`pn-scatter`)
    is assembled identically on the BiCGSTAB and source-iteration paths.
    """

    @_DEFERRED_2D_SI
    def test_bicgstab_p0_matches_si_p0(self):
        """BiCGSTAB and source iteration must agree at P0."""
        case = get("sn_slab_2eg_1rg")
        mix = next(iter(case.materials.values()))
        mesh = _uniform_2d(2, 2, 0.5, np.zeros((2, 2), dtype=int))
        quad = Quadrature.lebedev(order=17)

        keffs = {}
        for label, solver_type in [("SI", "source_iteration"), ("BC", "krylov")]:
            solver = SNSolver(SNMesh(mesh, quad, {0: mix}), inner_solver=solver_type, scattering_order=0, max_inner=500 if solver_type == "source_iteration" else 2000, inner_tol=1e-10 if solver_type == "source_iteration" else 1e-6)
            phi = solver.initial_flux_distribution()
            keff = 1.0
            for _ in range(50):
                fs = solver.compute_fission_source(phi, keff)
                phi = solver.solve_fixed_source(fs, phi)
                keff = solver.compute_keff(phi)
                phi /= np.linalg.norm(phi)
            keffs[label] = keff

        assert abs(keffs["SI"] - keffs["BC"]) < 1e-4, (
            f"P0 SI keff={keffs['SI']:.8f} vs BC keff={keffs['BC']:.8f}"
        )

    def test_bicgstab_p1_homogeneous_same_as_p0(self):
        """BiCGSTAB with P1 on homogeneous must match P0 (isotropic flux)."""
        mix = get_mixture("A", "2g")
        mesh = _uniform_2d(2, 2, 0.5, np.zeros((2, 2), dtype=int))
        quad = Quadrature.lebedev(order=17)

        keffs = {}
        for L in [0, 1]:
            solver = SNSolver(SNMesh(mesh, quad, {0: mix}), inner_solver="krylov", scattering_order=L, max_inner=2000, inner_tol=1e-6)
            phi = solver.initial_flux_distribution()
            keff = 1.0
            for _ in range(50):
                fs = solver.compute_fission_source(phi, keff)
                phi = solver.solve_fixed_source(fs, phi)
                keff = solver.compute_keff(phi)
                phi /= np.linalg.norm(phi)
            keffs[L] = keff

        assert abs(keffs[0] - keffs[1]) < 1e-4, (
            f"BiCGSTAB P0 keff={keffs[0]:.6f} vs P1 keff={keffs[1]:.6f} "
            f"should be equal on homogeneous"
        )

    @_DEFERRED_2D_SI
    def test_bicgstab_p1_matches_si_p1_homogeneous(self):
        """BiCGSTAB and source iteration must agree at P1 on homogeneous."""
        mix = get_mixture("A", "2g")
        mesh = _uniform_2d(2, 2, 0.5, np.zeros((2, 2), dtype=int))
        quad = Quadrature.lebedev(order=17)

        keffs = {}
        for label, solver_type in [("SI", "source_iteration"), ("BC", "krylov")]:
            solver = SNSolver(SNMesh(mesh, quad, {0: mix}), inner_solver=solver_type, scattering_order=1, max_inner=500 if solver_type == "source_iteration" else 2000, inner_tol=1e-10 if solver_type == "source_iteration" else 1e-6)
            phi = solver.initial_flux_distribution()
            keff = 1.0
            for _ in range(50):
                fs = solver.compute_fission_source(phi, keff)
                phi = solver.solve_fixed_source(fs, phi)
                keff = solver.compute_keff(phi)
                phi /= np.linalg.norm(phi)
            keffs[label] = keff

        assert abs(keffs["SI"] - keffs["BC"]) < 1e-3, (
            f"P1 SI keff={keffs['SI']:.8f} vs BC keff={keffs['BC']:.8f}"
        )


class TestSolveFixedSource:
    """Integration test: one outer iteration must reduce residual."""

    @_DEFERRED_2D_SI
    def test_source_iteration_converges(self):
        fuel = get_mixture("A", "2g")
        mod = get_mixture("B", "2g")
        materials = {2: fuel, 0: mod}
        nx, ny = 6, 4
        mat = np.zeros((nx, ny), dtype=int)
        mat[:3, :] = 2
        mesh = _uniform_2d(nx, ny, 0.2, mat)
        quad = Quadrature.lebedev(order=17)
        solver = SNSolver(SNMesh(mesh, quad, materials))

        phi = solver.initial_flux_distribution()
        keff = 1.0

        fission_src = solver.compute_fission_source(phi, keff)
        phi_new = solver.solve_fixed_source(fission_src, phi)

        assert not np.allclose(phi, phi_new), "No update from solve_fixed_source"
        assert np.all(np.isfinite(phi_new)), "NaN/Inf in solve output"

    @_DEFERRED_2D_SI
    def test_bicgstab_matches_source_iteration(self):
        """BiCGSTAB and source iteration must converge to the same keff."""
        case = get("sn_slab_2eg_1rg")
        mix = next(iter(case.materials.values()))
        mesh = _uniform_2d(2, 2, 0.5, np.zeros((2, 2), dtype=int))
        quad = Quadrature.lebedev(order=17)

        solver_si = SNSolver(SNMesh(mesh, quad, {0: mix}), inner_solver="source_iteration", max_inner=500, inner_tol=1e-10)
        phi = solver_si.initial_flux_distribution()
        keff = 1.0
        for _ in range(50):
            fs = solver_si.compute_fission_source(phi, keff)
            phi = solver_si.solve_fixed_source(fs, phi)
            keff = solver_si.compute_keff(phi)
            phi /= np.linalg.norm(phi)
        keff_si = keff

        solver_bc = SNSolver(SNMesh(mesh, quad, {0: mix}), inner_solver="krylov", max_inner=2000, inner_tol=1e-6)
        phi = solver_bc.initial_flux_distribution()
        keff = 1.0
        for _ in range(50):
            fs = solver_bc.compute_fission_source(phi, keff)
            phi = solver_bc.solve_fixed_source(fs, phi)
            keff = solver_bc.compute_keff(phi)
            phi /= np.linalg.norm(phi)
        keff_bc = keff

        assert abs(keff_si - keff_bc) < 1e-5, (
            f"BiCGSTAB keff={keff_bc:.8f} vs SI keff={keff_si:.8f}"
        )
