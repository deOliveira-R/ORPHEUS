"""Regression tests for the SN sweep + SNMesh stencil.

These tests cover bugs and edge cases found during the geometry
migration (2026-04-04).  Each test targets a specific failure mode.

Historical context — Gotcha #5 (resolved by Issue #196 Phase G Step 2.5)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
An algebraically-equivalent rewrite of the DD recurrence
(``0.5*(psi_in + psi_out)`` → ``2*(psi_in*(1-a)/(1-a+eps) + s) - psi_in``)
once produced divergent multi-group scattering source iteration
(flux grew by ~1e34 per outer iteration); 1-group was unaffected.
The recurrence-formulation helpers (``_solve_recurrence``,
``_outgoing``) were retired in Phase G Step 2.5 along with the slab
cumprod sweep path, so the original unit-level regressions for those
helpers are no longer applicable.  The end-to-end multi-group
scattering convergence regression (catches ERR-005) is preserved at
:class:`TestScatteringConvergence` below — it exercises the same
failure mode via :class:`~orpheus.sn.solver.SNSolver`, which is the
right level for a structural-bug guard now that the helpers are gone.
"""

import numpy as np
import pytest

from orpheus.geometry import (
    BC,
    Mesh1D,
    Mesh2D,
    Region,
    RegionMesh,
    StructuredGeometry,
)
from orpheus.sn.mesh.augmented_mesh import SNMesh


def _homogeneous_slab_mesh(n_cells: int, total_width: float, mat_id: int = 0) -> Mesh1D:
    """Single-region Cartesian mesh helper (replaces legacy ``homogeneous_1d``)."""
    geom = StructuredGeometry(
        geometry="SLB",
        regions=(Region(mat_id=mat_id, outer_thickness_cm=total_width),),
        bcs=(BC.reflective, BC.reflective),
    )
    return Mesh1D.from_geometry(geom, region_meshes=(RegionMesh(n_cells=n_cells),))
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.loss_representation import transport_sweep
from tests.sn._test_helpers import placeholder_materials

pytestmark = pytest.mark.l0  # SN sweep + SNMesh stencil regressions


# ═══════════════════════════════════════════════════════════════════════
# Multi-group scattering convergence regression (ERR-005)
# ═══════════════════════════════════════════════════════════════════════

class TestScatteringConvergence:
    """End-to-end regression for multi-group scattering convergence.

    Phase G Step 2.5 retired the slab cumprod recurrence (``_solve_recurrence``,
    ``_outgoing``); the unit-level tests that imported those helpers are gone.
    This class preserves the structural regression for ERR-005 at the
    :class:`~orpheus.sn.solver.SNSolver` level — the right level once the
    helper functions no longer exist.
    """

    @pytest.mark.catches("ERR-005")
    def test_regression_multigroup_scattering_convergence(self):
        """Scattering source iteration must converge for multi-group.

        Gotcha #5 (historical, see module docstring): an algebraically-
        equivalent rewrite of the slab DD recurrence caused multi-group
        scattering iteration to diverge while 1-group was unaffected.
        Catches that failure mode at the SNSolver end-to-end level.
        """
        from orpheus.derivations.common.xs_library import get_mixture
        from orpheus.sn.solver import SNSolver

        mix = get_mixture("A", "2g")
        mesh = _homogeneous_slab_mesh(20, 2.0, mat_id=0)
        quad = Quadrature.gauss_legendre(8)
        sn_mesh = SNMesh(mesh, quad, {0: mix})
        solver = SNSolver(sn_mesh, max_inner=500, inner_tol=1e-10)

        # One outer iteration: flux must remain bounded
        phi = solver.initial_flux_distribution()
        fission = solver.compute_fission_source(phi, 1.0)
        phi_new = solver.solve_fixed_source(fission, phi)

        assert np.all(np.isfinite(phi_new)), "Non-finite flux after one outer iteration"
        assert phi_new.max() < 100, (
            f"Flux blew up to {phi_new.max():.2e} — "
            f"scattering iteration may have regressed"
        )


# ═══════════════════════════════════════════════════════════════════════
# SNMesh stencil and shape tests
# ═══════════════════════════════════════════════════════════════════════

class TestSNMesh:
    """Tests for the SNMesh augmented geometry."""

    def test_stencil_values_cartesian(self):
        """streaming(0)[n,i] is the RAW down-face streaming |μ_x[n]| / dx[i]
        (#240) — the diamond 2 = 1/w_DD is the scheme's, applied in the cell
        kernel, NOT baked into this geometric accessor."""
        mesh = Mesh1D(edges=np.array([0.0, 0.1, 0.3, 0.6]),
                      mat_ids=np.array([0, 1, 2]))
        quad = Quadrature.gauss_legendre(4)
        sn_mesh = SNMesh(mesh, quad, placeholder_materials(mat_ids=(0, 1, 2)))

        for n in range(quad.N):
            for i in range(sn_mesh.nx):
                expected = abs(quad.mu_x[n]) / mesh.widths[i]
                np.testing.assert_allclose(
                    sn_mesh.streaming(0)[n, i], expected, rtol=1e-14,
                )

    def test_stencil_exactly_equals_axis_cosines_build(self):
        """G-b3 (C3.6): EXACT equality (not allclose) of the d-generic
        stencil against the hand-spelled ``2|axis_cosines(a)|/Δa`` —
        the bit-identity guard.  ``mu_x``/``mu_y`` are property VIEWS of
        ``axis_cosines(0)/(1)``, so the ``range(ndim)`` build must be
        bit-identical to the retired hand-listed x/y pair; if a future
        refactor turns the view into a rounded copy, the rtol gates above
        stay green while this one trips."""
        mesh = Mesh2D(
            edges_x=np.linspace(0, 1, 4),
            edges_y=np.linspace(0, 0.5, 3),
            mat_map=np.zeros((3, 2), dtype=int),
        )
        quad = Quadrature.lebedev(order=17)
        sn_mesh = SNMesh(mesh, quad, placeholder_materials())
        for a, widths in enumerate((mesh.dx, mesh.dy)):
            np.testing.assert_array_equal(
                sn_mesh.streaming(a),
                np.abs(quad.axis_cosines(a))[:, None] / widths[None, :],
            )

    def test_streaming_axis_out_of_range_raises(self):
        """G-b5 (C3.6): ``streaming(2)`` on a 2-D mesh is an IndexError —
        no phantom third axis."""
        mesh = Mesh2D(
            edges_x=np.linspace(0, 1, 3),
            edges_y=np.linspace(0, 1, 3),
            mat_map=np.zeros((2, 2), dtype=int),
        )
        sn_mesh = SNMesh(
            mesh, Quadrature.lebedev(order=5), placeholder_materials(),
        )
        with pytest.raises(IndexError, match="out of range for ndim=2"):
            sn_mesh.streaming(2)

    def test_stencil_dd_denom_equivalence(self):
        """Precomputed stencil must reproduce the original DD denominator.

        Original: denom = Σ_t + 2|μ_x|/dx + 2|μ_y|/dy
        Stencil:  denom = Σ_t + 2·streaming(0)[n,i] + 2·streaming(1)[n,j]
                  (#240: ``streaming`` is the RAW g = |μ|/Δ; the DD scheme
                  applies the diamond 2 = 1/w_DD, so the denom doubles each
                  streaming term).
        """
        mesh = Mesh2D(
            edges_x=np.linspace(0, 1, 4),  # 3 cells, dx=1/3
            edges_y=np.linspace(0, 0.5, 3),  # 2 cells, dy=0.25
            mat_map=np.zeros((3, 2), dtype=int),
        )
        quad = Quadrature.lebedev(order=17)
        sn_mesh = SNMesh(mesh, quad, placeholder_materials())

        sig_t = 0.5  # scalar for simplicity
        for n in range(quad.N):
            for i in range(sn_mesh.nx):
                for j in range(sn_mesh.spatial_shape[1]):
                    old = sig_t + 2*abs(quad.mu_x[n])/mesh.dx[i] + 2*abs(quad.mu_y[n])/mesh.dy[j]
                    new = sig_t + 2*sn_mesh.streaming(0)[n, i] + 2*sn_mesh.streaming(1)[n, j]
                    np.testing.assert_allclose(new, old, rtol=1e-14)

    def test_mesh1d_shapes(self):
        """SNMesh from Mesh1D must have rank-1 (N,) shaped mat_map and volumes."""
        mesh = Mesh1D(edges=np.linspace(0, 1, 6), mat_ids=np.array([0,1,2,1,0]))
        quad = Quadrature.gauss_legendre(4)
        sn_mesh = SNMesh(mesh, quad, placeholder_materials(mat_ids=(0, 1, 2)))

        assert sn_mesh.nx == 5
        assert sn_mesh.spatial_shape == (5,)
        assert sn_mesh.mat_map.shape == sn_mesh.spatial_shape
        assert sn_mesh.volumes.shape == sn_mesh.spatial_shape
        assert sn_mesh.is_1d is True

    def test_mesh2d_shapes(self):
        """SNMesh from Mesh2D preserves shapes."""
        mesh = Mesh2D(
            edges_x=np.linspace(0, 1, 4),
            edges_y=np.linspace(0, 1, 3),
            mat_map=np.zeros((3, 2), dtype=int),
        )
        quad = Quadrature.lebedev(order=17)
        sn_mesh = SNMesh(mesh, quad, placeholder_materials())

        assert sn_mesh.nx == 3
        assert sn_mesh.spatial_shape[1] == 2
        assert sn_mesh.mat_map.shape == (3, 2)
        assert sn_mesh.volumes.shape == (3, 2)
        assert sn_mesh.is_1d is False

    def test_cylindrical_requires_level_quadrature(self):
        """Cylindrical coords require a quadrature with level structure."""
        from orpheus.geometry import CoordSystem

        mesh = Mesh1D(edges=np.array([0.0, 1.0]), mat_ids=np.array([0]),
                      coord=CoordSystem.CYLINDRICAL)
        quad = Quadrature.gauss_legendre(4)
        with pytest.raises(ValueError, match="level structure"):
            SNMesh(mesh, quad, placeholder_materials())

    def test_spherical_setup(self):
        """Spherical SNMesh must precompute face areas and α coefficients."""
        from orpheus.geometry import CoordSystem

        mesh = Mesh1D(edges=np.array([0.0, 0.5, 1.0]), mat_ids=np.array([0, 1]),
                      coord=CoordSystem.SPHERICAL)
        quad = Quadrature.gauss_legendre(4)
        sn_mesh = SNMesh(mesh, quad, placeholder_materials(mat_ids=(0, 1)))

        assert sn_mesh.curvature == "spherical"
        assert sn_mesh.face_areas is not None
        assert sn_mesh.reduced.alpha_half is not None
        assert len(sn_mesh.reduced.alpha_half) == quad.N + 1
        # α_{1/2} = 0 and α_{N+1/2} ≈ 0
        np.testing.assert_allclose(sn_mesh.reduced.alpha_half[0], 0.0)
        np.testing.assert_allclose(sn_mesh.reduced.alpha_half[-1], 0.0, atol=1e-14)

    def test_sweep_1d_2d_consistency(self):
        r"""1D and 2D sweeps on equivalent meshes must produce the same keff.

        A ``Mesh1D`` slab (Gauss-Legendre) and a degenerate ``ny=1``
        ``Mesh2D`` must both recover the closed-form ``k_inf`` — a
        regression for issue #214 (the 2-D wavefront on a ``ny=1`` mesh
        no longer crashes once the ``(N, ng, *spatial)`` rank-d carve
        made the extent-1 axis handling clean).

        **Quadrature choice (by symmetry-group affinity, #154 G-tags).**
        The 2-D leg uses a **product** quadrature
        (``Quadrature.product(n_mu, n_phi)`` = Gauss-Legendre on
        :math:`\mu` × equispaced on :math:`\phi`), NOT Lebedev. A
        degenerate ``ny=1`` mesh is physically a slab in :math:`x`, whose
        symmetry is **SO(2)** (axial) — exactly the ``invariance_group``
        of the 1-D ``gauss_legendre`` reference. The product rule is
        SO(2)-invariant and *is* the GL polar rule lifted to 2-D, so the
        1-D-vs-2-D comparison is apples-to-apples in the polar angle.
        Lebedev is an :math:`O_h`-invariant ``S^2`` cubature built for
        SO(3)/spherical-harmonic-moment integration (PN, anisotropic
        scattering); on this flux-flat ``k_inf`` sweep it imposes a
        spurious cubic symmetry and ``degree_of_exactness=17`` over-kill
        at ``N=110`` ordinates — ~500× slower than ``product(8,4)``'s
        ``N=32`` for identical (exact) ``k_inf``. ``n_phi=4`` (not 2) is
        deliberate: it keeps genuine :math:`\Omega_y\neq 0` ordinates so
        the 2-D y-wavefront on the ``ny=1`` mesh is actually exercised —
        the exact code path #214 was about.
        """
        from orpheus.derivations.common.xs_library import get_mixture
        from orpheus.sn.solver import SNSolver

        mix = get_mixture("A", "1g")

        # 1D with GL
        mesh_1d = _homogeneous_slab_mesh(10, 1.0, mat_id=0)
        quad_gl = Quadrature.gauss_legendre(8)
        solver_1d = SNSolver(SNMesh(mesh_1d, quad_gl, {0: mix}), max_inner=500, inner_tol=1e-10)
        phi = solver_1d.initial_flux_distribution()
        keff_1d = 1.0
        for _ in range(50):
            fs = solver_1d.compute_fission_source(phi, keff_1d)
            phi = solver_1d.solve_fixed_source(fs, phi)
            keff_1d = solver_1d.compute_keff(phi)

        # 2D degenerate ny=1 — SO(2)-affine product quadrature (see docstring):
        # GL(8) polar (matches the 1-D reference) × 4 azimuthal angles, so
        # Ω_y≠0 ordinates genuinely exercise the y-wavefront on ny=1.
        mesh_2d = Mesh2D(
            edges_x=np.linspace(0, 1, 11),
            edges_y=np.array([0.0, 1.0]),
            mat_map=np.zeros((10, 1), dtype=int),
        )
        quad_2d = Quadrature.product(n_mu=8, n_phi=4)
        solver_2d = SNSolver(SNMesh(mesh_2d, quad_2d, {0: mix}), max_inner=500, inner_tol=1e-10)
        phi = solver_2d.initial_flux_distribution()
        keff_2d = 1.0
        for _ in range(50):
            fs = solver_2d.compute_fission_source(phi, keff_2d)
            phi = solver_2d.solve_fixed_source(fs, phi)
            keff_2d = solver_2d.compute_keff(phi)

        # Both must match the analytical k_inf (homogeneous, 1G)
        from orpheus.derivations import get
        k_ref = get("sn_slab_1eg_1rg").k_inf
        assert abs(keff_1d - k_ref) < 1e-8
        assert abs(keff_2d - k_ref) < 1e-8
