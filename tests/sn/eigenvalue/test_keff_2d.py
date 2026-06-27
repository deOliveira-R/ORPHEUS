"""2-D Cartesian SN k-eigenvalue verification.

Split from the legacy ``tests/sn/test_solver_components.py`` (SN
taxonomy reorg). This file is the 2-D eigenvalue home — it fills the
gap the reorg spec flagged ("NO 2-D k-eigenvalue test exists" as a
standalone). It carries the homogeneous-exact k_inf checks, the
multigroup-eigenvector recovery, the quadrature-normalisation
eigenvalue invariants (BiCGSTAB GL-vs-Lebedev, ERR-004), the Pn
eigenvalue effects, and the SI-vs-Krylov agreement.

Wave O "2-D SI Phase A" (2026-06-04): the 2-D Cartesian
source-iteration inner solver is now LIVE. The historical
``NotImplementedError`` guard (``SNSolver._solve_source_iteration``)
was stale — its "B1'' face block is 1-D-only" reason never existed as
code; B1'' was a legacy 1-D boundary closure superseded by the L2
``BoundaryFlux`` + ``SNBoundaryOperator`` bare-boundary architecture.
The SI inner is the structural twin of the Krylov inner (same operator
triple + RHS, only the driver differs). The previously-deferred 2-D SI
tests (formerly ``xfail(NotImplementedError)``) now run; the dedicated
:class:`TestSIKrylov2DEquivalence` gate pins SI ≡ Krylov ≡ closed-form
``k_inf`` on the repaired DEFAULT ``solve_sn`` entry.
"""

from pathlib import Path

import numpy as np
import pytest

from orpheus.derivations import get
from orpheus.derivations.common.xs_library import get_mixture
from orpheus.geometry import Mesh1D, Mesh2D
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.sn.solver import SNSolver, solve_sn

pytestmark = pytest.mark.l0  # SN 2-D eigenvalue component checks


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

    @pytest.mark.parametrize("ng_key,label", [("1g", "1G"), ("2g", "2G"), ("4g", "4G")])
    def test_homogeneous_exact(self, ng_key, label):
        case = get(f"sn_slab_{ng_key[0]}eg_1rg")
        mix = next(iter(case.materials.values()))
        materials = {0: mix}

        mesh = _uniform_2d(2, 2, 0.5, np.zeros((2, 2), dtype=int))
        # Homogeneous infinite-medium k_inf = νΣ_f/Σ_a is flux-SHAPE-INDEPENDENT
        # → every quadrature gives the SAME eigenvalue. Genuine 2-D Cartesian box
        # ⟹ O_h symmetry, so the SN-canonical level-symmetric set is the right
        # tool; Lebedev (O_h N=110 doe=17) is the SO(3) moment cubature, overkill
        # here. level_symmetric(4) (SAME O_h group, N=24 doe=3). Verified:
        # err ≤ 2.76e-12 (1G exact, 2G/4G to round-off) at the EXISTING 1e-8 tol.
        quad = Quadrature.level_symmetric(sn_order=4)
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
        # The group-spectrum eigenvector of a HOMOGENEOUS medium is the
        # quadrature-INDEPENDENT (Σ_t−Σ_s^T)⁻¹·χ⊗(νΣ_f) eigenvector — the spatial
        # flux is flat, so the angular cubature does not enter the group ratio.
        # Genuine 2-D Cartesian box ⟹ O_h; level_symmetric(4) (O_h N=24 doe=3)
        # replaces Lebedev (O_h N=110 doe=17 moment cubature). Verified: group
        # ratio matches the analytical eigenvector to 1.33e-11 (rtol 1e-6).
        quad = Quadrature.level_symmetric(sn_order=4)
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

    def test_p0_gives_identical_keff(self):
        """scattering_order=0 must give the exact same keff as the default."""
        case = get("sn_slab_2eg_1rg")
        mix = next(iter(case.materials.values()))
        mesh = _uniform_2d(2, 2, 0.5, np.zeros((2, 2), dtype=int))
        # Default-P0 vs explicit-scattering_order=0 EQUIVALENCE on the SAME quad
        # (the two solvers share this one `quad`) → the agreement holds for ANY
        # quadrature. Genuine 2-D Cartesian box ⟹ O_h; level_symmetric(4) (O_h
        # N=24 doe=3) replaces Lebedev (O_h N=110 doe=17). Verified: |Δ|=0.0
        # (bit-identical) at the EXISTING 1e-14 tol.
        quad = Quadrature.level_symmetric(sn_order=4)

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

    def test_p1_changes_heterogeneous_keff(self):
        """P1 scattering must produce a different keff than P0 on a
        heterogeneous problem where anisotropy matters at interfaces."""
        fuel = get_mixture("A", "2g")
        mod = get_mixture("B", "2g")  # B has mu_bar=0.6, strongly anisotropic
        materials = {2: fuel, 0: mod}

        mat = np.zeros((6, 2), dtype=int)
        mat[:3, :] = 2
        mesh = _uniform_2d(6, 2, 0.2, mat)
        # Genuine 2-D Cartesian (nx>1, ny>1, x-heterogeneous fuel|mod) ⟹ O_h
        # symmetry, so the SN-canonical level-symmetric set is the right tool;
        # Lebedev is the SO(3) spherical-harmonic-MOMENT cubature (over-quad
        # here: O_h N=110 doe=17). level_symmetric(4) is the SAME O_h group at
        # N=24 doe=3, and doe=3 EXACTLY integrates the degree-2 Y₁·Y₁ moment
        # products the P1 scattering source needs, so the anisotropy effect is
        # genuine — not vacuous. The assertion is EFFECT-PRESENCE (P1 keff
        # measurably ≠ P0 keff, |Δ|>1e-4), not a value vs a quadrature-dependent
        # reference, and both P0/P1 legs use the SAME swapped quad. Verified:
        # |Δ|=4.10e-3 (33× the 1e-4 bar; Lebedev gives 3.34e-3 — same effect).
        # 45.2s → 8.7s.
        quad = Quadrature.level_symmetric(sn_order=4)

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

    def test_bicgstab_p0_matches_si_p0(self):
        """BiCGSTAB and source iteration must agree at P0."""
        case = get("sn_slab_2eg_1rg")
        mix = next(iter(case.materials.values()))
        mesh = _uniform_2d(2, 2, 0.5, np.zeros((2, 2), dtype=int))
        # SI-vs-Krylov EQUIVALENCE at P0 on the SAME quad (both inners solve the
        # identical operator) → holds for ANY quadrature; swap both legs. Genuine
        # 2-D Cartesian box ⟹ O_h; level_symmetric(4) (O_h N=24 doe=3) replaces
        # Lebedev (O_h N=110 doe=17). Verified: |SI−BC|=1.36e-11 (<1e-4).
        quad = Quadrature.level_symmetric(sn_order=4)

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
        # P0-vs-P1 EQUIVALENCE on a HOMOGENEOUS medium: the flux is isotropic, so
        # the P1 (ℓ=1) moments vanish and P1 ≡ P0 for ANY quadrature; swap both
        # legs (same `quad`). Genuine 2-D Cartesian box ⟹ O_h; level_symmetric(4)
        # (O_h N=24 doe=3) replaces Lebedev (O_h N=110 doe=17). Verified:
        # |k0−k1|=3.88e-13 (<1e-4).
        quad = Quadrature.level_symmetric(sn_order=4)

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

    def test_bicgstab_p1_matches_si_p1_homogeneous(self):
        """BiCGSTAB and source iteration must agree at P1 on homogeneous."""
        mix = get_mixture("A", "2g")
        mesh = _uniform_2d(2, 2, 0.5, np.zeros((2, 2), dtype=int))
        # SI-vs-Krylov EQUIVALENCE at P1 on a HOMOGENEOUS medium, SAME quad (both
        # inners solve the identical P1 operator) → holds for ANY quadrature; swap
        # both legs. Genuine 2-D Cartesian box ⟹ O_h; level_symmetric(4) (O_h
        # N=24 doe=3) replaces Lebedev (O_h N=110 doe=17). Verified:
        # |SI−BC|=1.40e-11 (<1e-3).
        quad = Quadrature.level_symmetric(sn_order=4)

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

    def test_source_iteration_converges(self):
        fuel = get_mixture("A", "2g")
        mod = get_mixture("B", "2g")
        materials = {2: fuel, 0: mod}
        nx, ny = 6, 4
        mat = np.zeros((nx, ny), dtype=int)
        mat[:3, :] = 2
        mesh = _uniform_2d(nx, ny, 0.2, mat)
        # The assertion is quadrature-AGNOSTIC: solve_fixed_source must produce a
        # finite, non-identical update (one SI step moves the flux) — true for any
        # consistent SN set. Genuine 2-D Cartesian (nx>1, ny>1) ⟹ O_h;
        # level_symmetric(4) (O_h N=24 doe=3) replaces Lebedev (O_h N=110 doe=17).
        quad = Quadrature.level_symmetric(sn_order=4)
        solver = SNSolver(SNMesh(mesh, quad, materials))

        phi = solver.initial_flux_distribution()
        keff = 1.0

        fission_src = solver.compute_fission_source(phi, keff)
        phi_new = solver.solve_fixed_source(fission_src, phi)

        assert not np.allclose(phi, phi_new), "No update from solve_fixed_source"
        assert np.all(np.isfinite(phi_new)), "NaN/Inf in solve output"

    def test_bicgstab_matches_source_iteration(self):
        """BiCGSTAB and source iteration must converge to the same keff."""
        case = get("sn_slab_2eg_1rg")
        mix = next(iter(case.materials.values()))
        mesh = _uniform_2d(2, 2, 0.5, np.zeros((2, 2), dtype=int))
        # SI-vs-Krylov EQUIVALENCE on the SAME quad (both inners solve the
        # identical operator) → holds for ANY quadrature; swap both legs. Genuine
        # 2-D Cartesian box ⟹ O_h; level_symmetric(4) (O_h N=24 doe=3) replaces
        # Lebedev (O_h N=110 doe=17). Verified: |SI−BC|=1.36e-11 (<1e-5).
        quad = Quadrature.level_symmetric(sn_order=4)

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


class TestSIKrylov2DEquivalence:
    """2-D Cartesian eigenvalue source-iteration ≡ Krylov ≡ closed-form k_inf.

    The gate for the Wave O "2-D SI Phase A" carve (the stale
    ``NotImplementedError`` guard at ``SNSolver._solve_source_iteration``
    deleted).  ``solve_sn`` defaults ``inner_solver="source_iteration"`` for
    every geometry, so the stale guard broke the DEFAULT 2-D Cartesian
    eigenvalue entry point — these tests drive the production ``solve_sn``
    entry directly (not a hand-rolled power-iteration loop), so they pin the
    thing users actually call.

    Structural-independence discipline (vv-principles L11 / L14): SI ≡ Krylov
    alone is twin-path agreement — NECESSARY, NOT SUFFICIENT (both could
    share a defect).  It becomes correctness evidence ONLY because the same
    production path is independently pinned to the closed-form
    ``k_inf = λ_max(A⁻¹F)`` on the homogeneous leg
    (:meth:`test_default_entry_hits_kinf`).  The heterogeneous leg has no 2-D
    closed-form eigenvalue reference, so its bar is SI≡Krylov flux-SHAPE
    agreement + convergence-to-the-same-limit, anchored by the k_inf leg.
    The heterogeneous flux is genuinely NON-FLAT (≥2G fuel|moderator), so the
    angular / wavefront redistribution terms are active — not a
    homogeneous/1G degenerate (vv-principles anti-patterns #3 / #4).
    """

    @pytest.mark.parametrize("ng_key", ["1g", "2g"])
    @pytest.mark.verifies("transport-cartesian-2d", "matrix-eigenvalue")
    def test_default_entry_hits_kinf(self, ng_key):
        """The repaired DEFAULT entry: ``solve_sn`` (source-iteration inner by
        default) on a 2-D mesh must not raise and must hit closed-form k_inf.

        Pre-carve this raised ``NotImplementedError`` (the stale 2-D SI guard).
        """
        case = get(f"sn_slab_{ng_key[0]}eg_1rg")
        mix = next(iter(case.materials.values()))
        mesh = _uniform_2d(2, 2, 0.5, np.zeros((2, 2), dtype=int))
        # Homogeneous infinite-medium k_inf is flux-SHAPE-INDEPENDENT → the
        # default-entry SI eigenvalue is the SAME for any quadrature. Genuine 2-D
        # Cartesian box ⟹ O_h; level_symmetric(4) (O_h N=24 doe=3) replaces
        # Lebedev (O_h N=110 doe=17). Verified: err ≤ 1.68e-11 at the 1e-8 tol.
        sol = solve_sn(
            {0: mix}, mesh, Quadrature.level_symmetric(sn_order=4),
            keff_tol=1e-12, flux_tol=1e-10, max_inner=500, inner_tol=1e-10,
        )
        assert np.isfinite(sol.keff)
        assert abs(sol.keff - case.k_inf) < 1e-8, (
            f"{ng_key}: default-entry SI keff={sol.keff:.10f} "
            f"vs closed-form k_inf={case.k_inf:.10f}"
        )

    @pytest.mark.slow
    @pytest.mark.verifies(
        "transport-cartesian-2d", "matrix-eigenvalue", "multigroup",
    )
    def test_si_krylov_heterogeneous_2g_nonflat_flux(self):
        """SI ≡ Krylov on a HETEROGENEOUS, non-flat, 2G 2-D problem.

        Flux SHAPE (not just the shape-independent eigenvalue) must agree
        between the two inner solvers.  The k_inf anchor (above) supplies the
        structural-independence leg; here the bar is the twin agreement on a
        problem where the redistribution terms are genuinely exercised.
        """
        materials = {2: get_mixture("A", "2g"), 0: get_mixture("B", "2g")}
        nx, ny = 8, 4
        mat = np.zeros((nx, ny), dtype=int)
        mat[:4, :] = 2  # fuel | moderator split across x → non-flat flux
        mesh = _uniform_2d(nx, ny, 0.25, mat)
        # SI-vs-Krylov EQUIVALENCE on a genuine 2-D Cartesian (8×4, fuel|mod)
        # heterogeneous case ⟹ O_h symmetry. Both inners solve the identical
        # (L+C−S−B)ψ=q operator on the SAME quadrature, so the equivalence holds
        # for ANY quadrature — swap both sides together. level_symmetric(4) (O_h,
        # N=24 doe=3) replaces Lebedev (O_h, N=110 doe=17 — the SO(3) moment
        # cubature, over-quadrature for a plain 2-D sweep). Verified: non-flat
        # guard max/min=1.455 (>1.2 fires), SI≡Krylov keff Δ=3.76e-10 (<1e-7) and
        # flux SHAPE max-diff 2.64e-9 (within rtol 1e-6/atol 1e-8). 138s → 14.6s.
        quad = Quadrature.level_symmetric(sn_order=4)

        sol_si = solve_sn(
            materials, mesh, quad, inner_solver="source_iteration",
            keff_tol=1e-12, flux_tol=1e-10, max_inner=500, inner_tol=1e-10,
        )
        sol_kry = solve_sn(
            materials, mesh, quad, inner_solver="krylov",
            keff_tol=1e-12, flux_tol=1e-10, max_inner=4000, inner_tol=1e-8,
        )

        phi_si = np.asarray(sol_si.scalar_flux.values, dtype=np.float64)
        phi_kry = np.asarray(sol_kry.scalar_flux.values, dtype=np.float64)

        # Degenerate-gate guard: the flux MUST be genuinely non-flat, else the
        # redistribution terms are not exercised and SI≡Krylov is vacuous.
        prof = phi_si[0].mean(axis=1)  # group-0 profile across x
        assert prof.max() / prof.min() > 1.2, (
            f"flux too flat (max/min={prof.max() / prof.min():.3f}); "
            "redistribution not exercised — gate degenerate"
        )

        # Eigenvalue twin agreement (shape-independent — necessary not enough).
        assert abs(sol_si.keff - sol_kry.keff) < 1e-7, (
            f"SI keff={sol_si.keff:.10f} vs Krylov keff={sol_kry.keff:.10f}"
        )
        # Flux SHAPE twin agreement (eigenvectors are scale-free → mean-norm).
        phi_si_n = phi_si / phi_si.mean()
        phi_kry_n = phi_kry / phi_kry.mean()
        np.testing.assert_allclose(
            phi_si_n, phi_kry_n, rtol=1e-6, atol=1e-8,
            err_msg="SI vs Krylov heterogeneous 2-D flux SHAPE diverged",
        )

    @pytest.mark.verifies(
        "transport-cartesian-2d", "matrix-eigenvalue", "multigroup",
    )
    def test_eigenvalue_jacobi_gauss_seidel_equivalence(self):
        """#218: the eigenvalue SI inner gives the SAME eigenvalue + eigenmode
        under the Jacobi (default) and boundary-Gauss-Seidel schedules.

        vv-principles Mode 9: the two boundary splittings reach the SAME
        within-group fixed point each outer step (only the inner SI spectral
        rate differs), so the outer power iteration converges to the same
        eigenvalue.  The eigenvalue inner being unable to reach the G-S
        accelerator (which the fixed-source path got in Phase 3) was the #218
        gap; the shared ``_within_group_si`` builder + ``inner_schedule``
        plumbing close it.  G-S stays OPT-IN (default Jacobi) — a schedule
        change shifts the converged k_eff by ~inner_tol (the inner SI stops at
        a slightly different point, NOT keff_tol-tight), which is why the
        keff_tol-tight regression snapshots stay on Jacobi.

        Reflective het 2-D config: ``B`` is non-trivial (G-S genuinely folds
        it) and the flux is non-flat (the equivalence is not vacuous).
        """
        materials = {2: get_mixture("A", "2g"), 0: get_mixture("B", "2g")}
        nx, ny = 8, 4
        mat = np.zeros((nx, ny), dtype=int)
        mat[:4, :] = 2  # fuel | moderator split across x → non-flat flux
        mesh = _uniform_2d(nx, ny, 0.25, mat)
        # Jacobi-vs-boundary-G-S schedule EQUIVALENCE on a genuine 2-D Cartesian
        # (8×4, fuel|mod) heterogeneous case ⟹ O_h symmetry. Both schedules reach
        # the SAME within-group fixed point on the SAME quadrature (only the inner
        # SI spectral rate differs), so the equivalence holds for ANY quadrature —
        # swap both legs together. level_symmetric(4) (O_h, N=24 doe=3) replaces
        # Lebedev (O_h, N=110 doe=17). Verified: non-flat guard max/min=1.455
        # (>1.2 fires, B coupling exercised), Jacobi≡G-S keff Δ=2.59e-11 (<1e-8)
        # and flux SHAPE max-diff 3.92e-10 (within rtol 1e-6/atol 1e-8). 37s → 7.7s.
        quad = Quadrature.level_symmetric(sn_order=4)
        kw = dict(keff_tol=1e-12, flux_tol=1e-10, max_inner=500, inner_tol=1e-10)

        sol_jac = solve_sn(materials, mesh, quad, inner_schedule="jacobi", **kw)
        sol_gs = solve_sn(
            materials, mesh, quad, inner_schedule="gauss_seidel", **kw,
        )

        # Non-degeneracy: the flux MUST be genuinely non-flat, else B +
        # redistribution are not exercised and the Jacobi≡G-S claim is vacuous.
        phi_jac = np.asarray(sol_jac.scalar_flux.values, dtype=np.float64)
        prof = phi_jac[0].mean(axis=1)
        assert prof.max() / prof.min() > 1.2, (
            f"flux too flat (max/min={prof.max() / prof.min():.3f}); "
            "B coupling not exercised — gate degenerate"
        )

        # Mode 9: same eigenvalue to within ~inner_tol (the schedule change
        # shifts the inner SI stopping point, NOT the fixed point — so the
        # k_eff agreement scales with inner_tol=1e-10, not keff_tol=1e-12).
        assert abs(sol_jac.keff - sol_gs.keff) < 1e-8, (
            f"#218: eigenvalue Jacobi keff={sol_jac.keff:.12f} vs "
            f"boundary-G-S keff={sol_gs.keff:.12f} (Δ exceeds the inner_tol-"
            "scale schedule shift — investigate as a real FP discrepancy)"
        )
        # Same eigenmode shape (eigenvectors are scale-free → mean-norm).
        phi_gs = np.asarray(sol_gs.scalar_flux.values, dtype=np.float64)
        np.testing.assert_allclose(
            phi_jac / phi_jac.mean(), phi_gs / phi_gs.mean(),
            rtol=1e-6, atol=1e-8,
            err_msg="#218: eigenvalue Jacobi vs boundary-G-S flux SHAPE diverged",
        )

    @pytest.mark.slow
    @pytest.mark.verifies(
        "transport-cartesian-2d", "matrix-eigenvalue", "multigroup",
    )
    def test_si_2d_keff_converges_under_refinement(self):
        """Leg (d): the 2-D SI keff is Cauchy under mesh refinement.

        Convergence is NECESSARY, NOT SUFFICIENT (vv-principles anti-pattern
        #5) — a monotone approach to the WRONG limit would still be Cauchy.
        The value is anchored by the k_inf + SI≡Krylov legs above; this leg
        catches a consistency regression (the discrete operator drifting off a
        single fixed point under refinement).
        """
        materials = {2: get_mixture("A", "2g"), 0: get_mixture("B", "2g")}
        # Genuine 2-D Cartesian (n×n, x-heterogeneous fuel|mod) ⟹ O_h symmetry,
        # so the SN-canonical level-symmetric set is the right tool; Lebedev is
        # the SO(3) moment cubature (over-quad: O_h N=110 doe=17). The assertion
        # is a CONVERGENCE TREND (the keff Cauchy diffs shrink, d2<d1, under mesh
        # refinement), which is quadrature-AGNOSTIC: any consistent SN set
        # converges, so the trend holds for any quadrature. level_symmetric(4)
        # (O_h, N=24 doe=3). Verified: d1(4→8)=4.34e-4, d2(8→16)=1.07e-4 (d2<d1).
        # 77s → 48s.
        quad = Quadrature.level_symmetric(sn_order=4)
        keffs = []
        for n in (4, 8, 16):
            mat = np.zeros((n, n), dtype=int)
            mat[: n // 2, :] = 2  # fuel in x<0.5 of the fixed [0,1]² domain
            mesh = _uniform_2d(n, n, 1.0 / n, mat)
            sol = solve_sn(
                materials, mesh, quad, inner_solver="source_iteration",
                keff_tol=1e-12, flux_tol=1e-10, max_inner=500, inner_tol=1e-10,
            )
            keffs.append(sol.keff)
        d1 = abs(keffs[1] - keffs[0])
        d2 = abs(keffs[2] - keffs[1])
        assert d2 < d1, (
            f"2-D SI keff not converging under refinement: "
            f"|Δ|(4→8)={d1:.3e}, |Δ|(8→16)={d2:.3e}; keffs={keffs}"
        )
