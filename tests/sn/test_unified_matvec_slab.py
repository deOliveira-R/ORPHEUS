r"""Comparison tests — unified slab matvec (PR-TYPED-6c Step 4).

Issue #197 PR-TYPED-6c Step 4 verifies the 1-D Cartesian (slab) branch
of :func:`~orpheus.sn.operator.transport_operator_matvec_unified`.

Slab is the M-M-neutral 1-level case of the per-level body that
sphere and cylinder also flow through:
  - level_indices = [arange(N)]
  - alpha = 0, redist_dAw = 0, tau = 1 (M-M closure is the identity)
  - A = 1 (slab face area)
  - psi_angular_upstream = None (no M-M angular recurrence)

The per-cell algebra :func:`cell_balance_for_streaming` collapses to
the slab form ``2|μ|·1 + Σ_t·V`` with neutral closure constants — no
algebraic shortcut, just neutral data.

Phase 1 (pole/inner) seed difference from curvilinear: slab uses
``bc_left.apply(ψ_view[:, :, 0, 0])`` (BC-applied cell-centre proxy
at x=0) where curvilinear uses ``ψ_view[:, :, 0, 0]`` directly
(cell-centre IS the pole-face proxy for r=0). The Phase 2 inward IC
uses ``bc_right.apply(ψ_view[:, :, -1, 0])`` symmetrically.

Important: the unified Cartesian matvec is WDD-based; the legacy
:func:`~orpheus.sn.operator.transport_operator_matvec` is FD-based.
These are DIFFERENT discretisations — they converge to the same
continuous solution under refinement but DIFFER by O(h) at any fixed
mesh.  The L1 verification here uses analytical k_∞ (shape-invariant,
so insensitive to discretisation order) and the FD-vs-WDD
characterisation is informational (NOT a pass/fail gate).

The full Case singular-eigenfunction Krylov verification is deferred
to PR-TYPED-6c Step 5 (when ``SNStreamingOperator.apply`` is rewired
to call ``transport_operator_matvec_unified`` and the production
``test_heterogeneous_transport.py`` runs through the unified path).
The standalone Krylov + Case combination is unusably slow at the
n_per=320 production mesh; that production test uses
``inner_solver="source_iteration"`` (the sweep, not the matvec).
"""
from __future__ import annotations

import contextlib

import numpy as np
import pytest

from orpheus.derivations.common.eigenvalue import kinf_homogeneous
from orpheus.derivations.common.xs_library import make_mixture
from orpheus.geometry import BC, CoordSystem, Mesh1D
from orpheus.sn import operator as sn_op
from orpheus.sn import solve_sn
from orpheus.sn.geometry import SNMesh
from orpheus.sn.operator import (
    build_equation_map,
    solution_to_angular_flux,
    transport_operator_matvec,
    transport_operator_matvec_unified,
)
from orpheus.sn.quadrature import GaussLegendre1D
from tests.sn._test_helpers import legacy_proxy_matvec, placeholder_materials


# ═══════════════════════════════════════════════════════════════════════
# L0 — sanity (zero, constant, 2D-raises)
# ═══════════════════════════════════════════════════════════════════════


def _build_slab(n_cells: int = 5, n_ord: int = 4) -> SNMesh:
    mesh = Mesh1D(
        edges=np.linspace(0.0, 2.0, n_cells + 1),
        mat_ids=np.zeros(n_cells, dtype=int),
        coord=CoordSystem.CARTESIAN,
        bc_left=BC("reflective"),
        bc_right=BC("reflective"),
    )
    quad = GaussLegendre1D.create(n_ordinates=n_ord)
    return SNMesh(mesh, quad, placeholder_materials())


@pytest.mark.l0
def test_unified_slab_zero_psi_gives_zero() -> None:
    """Linear operator: zero input → zero output."""
    sn_mesh = _build_slab(n_cells=5, n_ord=4)
    ng = 1
    sigma_t = np.full((ng, sn_mesh.nx, 1), 2.0)
    psi_view = np.zeros((sn_mesh.quad.N, ng, sn_mesh.nx, 1))

    m_unified = legacy_proxy_matvec(psi_view, sn_mesh, sigma_t)
    np.testing.assert_array_equal(m_unified, np.zeros_like(m_unified))


@pytest.mark.l0
def test_unified_slab_constant_psi_gives_sigma_t() -> None:
    """At ψ = constant on homogeneous reflective slab, unified matvec
    returns σ_t · ψ everywhere. Flat flux activates only the collision
    term (streaming cancels via WDD; M-M redistribution is identity
    on the slab-neutral data)."""
    sn_mesh = _build_slab(n_cells=5, n_ord=4)
    ng = 1
    sigma_t_val = 2.0
    sigma_t = np.full((ng, sn_mesh.nx, 1), sigma_t_val)
    psi_view = np.ones((sn_mesh.quad.N, ng, sn_mesh.nx, 1))

    m_unified = legacy_proxy_matvec(psi_view, sn_mesh, sigma_t)
    np.testing.assert_allclose(
        m_unified, sigma_t_val, rtol=1e-13, atol=1e-14,
    )


# ═══════════════════════════════════════════════════════════════════════
# L1 — homogeneous 2G k_∞ via the unified matvec (Krylov inner)
# ═══════════════════════════════════════════════════════════════════════
#
# k_∞ = (νΣ_f)·φ / Σ_a·φ is FLUX-SHAPE INDEPENDENT on homogeneous
# reflective.  Any consistent matvec converges to the same k_∞.  This
# pins the unified matvec's basic correctness (does it converge to a
# reasonable eigenvalue at all?) without the slow heterogeneous Case
# Krylov path.  Heterogeneous verification arrives at Step 5 when
# SNStreamingOperator.apply is rewired and the production
# test_heterogeneous_transport.py exercises the unified path.
# ═══════════════════════════════════════════════════════════════════════


def _make_2g_mixture(sigma_t, sig_s_matrix, nu_sigma_f, chi):
    sigma_t = np.asarray(sigma_t, dtype=float)
    sig_s = np.asarray(sig_s_matrix, dtype=float)
    nu_sig_f = np.asarray(nu_sigma_f, dtype=float)
    chi = np.asarray(chi, dtype=float)
    sig_a = sigma_t - sig_s.sum(axis=1)
    nu = np.ones_like(nu_sig_f)
    sig_f = nu_sig_f.copy()
    sig_c = sig_a - sig_f
    return make_mixture(
        sig_t=sigma_t, sig_c=sig_c, sig_f=sig_f, nu=nu, chi=chi, sig_s=sig_s,
    )


@contextlib.contextmanager
def _patch_apply_to_unified():
    """Route :meth:`SNStreamingOperator.apply` Cartesian path through the
    unified matvec (sphere/cylinder fall through to legacy).  The shim
    decodes packed ψ via :func:`solution_to_angular_flux` (which BC-fills
    the boundary cell-centre slots), calls the unified matvec on the
    canonical 4-D view, then re-packs at the unknown slots.
    """
    orig_apply = sn_op.SNStreamingOperator.apply

    def _unified_apply(self, psi_packed: np.ndarray) -> np.ndarray:
        sn_mesh = self.sn_mesh
        curv = getattr(sn_mesh, "curvature", None)
        if curv is not None:
            return orig_apply(self, psi_packed)

        eq_map = self._ensure_eq_map()
        quad = sn_mesh.quad
        ng = self.sig_t.shape[0]
        psi_view = solution_to_angular_flux(
            psi_packed, eq_map, quad, sn_mesh.nx, sn_mesh.ny, ng,
            bc_xmin=sn_mesh.bc_xmin, bc_xmax=sn_mesh.bc_xmax,
            bc_ymin=sn_mesh.bc_ymin, bc_ymax=sn_mesh.bc_ymax,
        )
        m_view = legacy_proxy_matvec(psi_view, sn_mesh, self.sig_t)
        flux = np.empty((ng, eq_map.n_eq))
        for k in range(eq_map.n_eq):
            flux[:, k] = m_view[
                eq_map.ordinate[k], :, eq_map.ix[k], eq_map.iy[k],
            ]
        return flux.ravel(order="F")

    sn_op.SNStreamingOperator.apply = _unified_apply
    try:
        yield
    finally:
        sn_op.SNStreamingOperator.apply = orig_apply


@pytest.mark.l1
@pytest.mark.parametrize("nx", [10, 20])
def test_unified_slab_l1_homogeneous_kinf_2g(nx: int) -> None:
    r"""L1 — 2G homogeneous reflective slab via the unified matvec.

    Drives :func:`solve_sn` with ``inner_solver="krylov"`` and a
    monkey-patched :meth:`SNStreamingOperator.apply` that routes the
    matvec through :func:`transport_operator_matvec_unified` for
    Cartesian meshes.  The converged ``k_eff`` is compared to the
    analytical ``k_∞ = (νΣ_f^T φ) / (Σ_a^T φ)`` reference at rtol < 5e-4
    — the same tolerance the production
    ``test_sn_spherical_homogeneous_kinf_recovery_2g`` uses.

    Per ``.claude/lessons.md`` L2: k_∞ is flux-shape independent on
    homogeneous reflective, so any consistent matvec hits the same
    eigenvalue.  Heterogeneous (shape-sensitive) verification lives at
    Step 5 when the production sweep is rewired.
    """
    sigma_t = [0.5, 1.0]
    sig_s = [[0.3, 0.05], [0.0, 0.7]]
    nu_sig_f = [0.4, 0.6]
    chi = [1.0, 0.0]
    mat = _make_2g_mixture(sigma_t, sig_s, nu_sig_f, chi)
    k_analytical = kinf_homogeneous(
        np.asarray(sigma_t), np.asarray(sig_s),
        np.asarray(nu_sig_f), np.asarray(chi),
    )

    mesh = Mesh1D(
        edges=np.linspace(0.0, 2.0, nx + 1),
        mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.CARTESIAN,
        bc_left=BC("reflective"),
        bc_right=BC("reflective"),
    )
    quad = GaussLegendre1D.create(n_ordinates=8)

    with _patch_apply_to_unified():
        sol = solve_sn(
            materials={0: mat},
            mesh=mesh,
            quadrature=quad,
            inner_solver="krylov",
            max_outer=200, keff_tol=1e-9, flux_tol=1e-8,
            max_inner=200, inner_tol=1e-10,
        )

    rel = abs(sol.keff - k_analytical) / k_analytical
    assert rel < 5e-4, (
        f"unified slab k_∞ recovery violated: "
        f"k_analytical={k_analytical:.10f}, k_unified={sol.keff:.10f}, "
        f"rel={rel:.2e}, nx={nx}"
    )


# ═══════════════════════════════════════════════════════════════════════
# Characterization — FD-vs-WDD order-of-accuracy delta (informational)
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.l0
def test_unified_slab_differs_from_legacy_fd_O_h() -> None:
    r"""Characterize the FD-vs-WDD delta on a homogeneous reflective slab.

    The legacy :func:`transport_operator_matvec` uses upwind finite
    differences via :func:`_compute_gradients` (1st-order at the
    interface, but on uniform homogeneous reflective the FD is
    consistent with WDD on flat ψ).  The unified
    :func:`transport_operator_matvec_unified` uses WDD via
    :func:`cell_balance_for_streaming`.

    On flat ψ both should produce ``σ_t · ψ`` exactly.  On a
    structured ψ they diverge by an O(magnitude) order-of-accuracy
    difference.  This test pins the divergence as a CHARACTERISATION,
    not a correctness gate — it just confirms the two methods are
    actually computing different discretisations as designed.
    """
    n_cells, n_ord, ng = 5, 4, 1
    sn_mesh = _build_slab(n_cells=n_cells, n_ord=n_ord)
    sigma_t = np.full((ng, n_cells, 1), 2.0)
    eq_map = build_equation_map(n_cells, 1, sn_mesh.quad, ng)

    rng = np.random.default_rng(0)
    psi_packed = rng.standard_normal(eq_map.n_unknowns)

    # Legacy FD
    m_legacy = transport_operator_matvec(
        psi_packed, eq_map, sn_mesh.quad, sigma_t, n_cells, 1, ng,
        sn_mesh.dx, sn_mesh.dy,
        bc_xmin=sn_mesh.bc_xmin, bc_xmax=sn_mesh.bc_xmax,
        bc_ymin=sn_mesh.bc_ymin, bc_ymax=sn_mesh.bc_ymax,
    )

    # Unified WDD
    psi_view = solution_to_angular_flux(
        psi_packed, eq_map, sn_mesh.quad, n_cells, 1, ng,
        bc_xmin=sn_mesh.bc_xmin, bc_xmax=sn_mesh.bc_xmax,
        bc_ymin=sn_mesh.bc_ymin, bc_ymax=sn_mesh.bc_ymax,
    )
    m_view = legacy_proxy_matvec(psi_view, sn_mesh, sigma_t)
    flux = np.empty((ng, eq_map.n_eq))
    for k in range(eq_map.n_eq):
        flux[:, k] = m_view[
            eq_map.ordinate[k], :, eq_map.ix[k], eq_map.iy[k],
        ]
    m_unified = flux.ravel(order="F")

    # Both finite, both produce the same result for FLAT ψ (collision term
    # only); both produce DIFFERENT results for non-flat ψ (different
    # discretisations).
    assert np.all(np.isfinite(m_legacy))
    assert np.all(np.isfinite(m_unified))
    diff = np.linalg.norm(m_legacy - m_unified)
    assert diff > 1e-3, (
        f"Expected legacy FD vs unified WDD to differ on random ψ "
        f"(different discretisation orders); got ||Δ||={diff:.3e}"
    )

    # Sanity: M·1 = σ_t = 2.0 for BOTH (collision-only on flat ψ).
    ones_packed = np.ones(eq_map.n_unknowns)
    m_legacy_flat = transport_operator_matvec(
        ones_packed, eq_map, sn_mesh.quad, sigma_t, n_cells, 1, ng,
        sn_mesh.dx, sn_mesh.dy,
        bc_xmin=sn_mesh.bc_xmin, bc_xmax=sn_mesh.bc_xmax,
        bc_ymin=sn_mesh.bc_ymin, bc_ymax=sn_mesh.bc_ymax,
    )
    psi_view_flat = solution_to_angular_flux(
        ones_packed, eq_map, sn_mesh.quad, n_cells, 1, ng,
        bc_xmin=sn_mesh.bc_xmin, bc_xmax=sn_mesh.bc_xmax,
        bc_ymin=sn_mesh.bc_ymin, bc_ymax=sn_mesh.bc_ymax,
    )
    m_view_flat = legacy_proxy_matvec(psi_view_flat, sn_mesh, sigma_t)
    flux_flat = np.empty((ng, eq_map.n_eq))
    for k in range(eq_map.n_eq):
        flux_flat[:, k] = m_view_flat[
            eq_map.ordinate[k], :, eq_map.ix[k], eq_map.iy[k],
        ]
    m_unified_flat = flux_flat.ravel(order="F")

    np.testing.assert_allclose(m_legacy_flat, 2.0, rtol=1e-13, atol=1e-14)
    np.testing.assert_allclose(m_unified_flat, 2.0, rtol=1e-13, atol=1e-14)
