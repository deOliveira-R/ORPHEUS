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
    transport_operator_matvec_unified,
)
from orpheus.numerics.quadrature import Quadrature
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
    quad = Quadrature.gauss_legendre(n_ordinates=n_ord)
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


# D-K.5 (2026-05-29): ``_patch_apply_to_unified`` retired together with
# :class:`SNStreamingOperator`.  Pre-D-K.1, this monkey-patched
# ``SNStreamingOperator.apply`` Cartesian path through the unified
# matvec.  Post-D-K.1 (commit ``400ca33``), ``SNSolver.L`` is
# ``StreamingOperator + CollisionOperator`` (= :class:`InvertibleOperator`),
# which routes through the unified matvec natively for 1-D and through
# :meth:`StreamingOperator._apply_2d_cartesian` for 2-D.  The
# monkey-patch had no effect on the call path post-D-K.1; deleting it
# removes the no-op + the dependency on the retiring
# :class:`SNStreamingOperator` class.


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
    quad = Quadrature.gauss_legendre(n_ordinates=8)

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


# D-H.2-C4e.6 (2026-05-29): retired
# ``test_unified_slab_differs_from_legacy_fd_O_h`` together with the
# legacy ``transport_operator_matvec`` it characterised against.
# The test pinned the FD-vs-WDD discretisation-order delta on a
# homogeneous reflective slab as informational; with the legacy FD
# retired (only the L2-native :meth:`StreamingOperator._apply_2d_
# cartesian_l2` survives), the cross-comparison has no architectural
# target.  Slab correctness is pinned by ``test_unified_slab_l1_
# homogeneous_kinf_2g`` (L1 k_∞ recovery) at L1, the C4a unit-source
# tests (``test_2d_l2_matvec_correctness.py``) at L0 for the 2-D
# matvec, and ``test_l1_standoff_slab_cylinder.py`` for cross-method
# slab anchoring at L1.
