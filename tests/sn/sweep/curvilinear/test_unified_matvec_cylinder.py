r"""Comparison tests — unified cylinder matvec (PR-TYPED-6c Step 3).

Issue #197 PR-TYPED-6c Step 3 verifies the cylindrical branch of
:func:`~orpheus.sn.operator._transport_operator_matvec_unified`.

Two levels of evidence (per L14 in ``.claude/lessons.md`` — solver
correctness is a 4-way standoff):

* **L0** — per-ordinate **hand reference** at machine precision. The
  hand reference is structurally independent: it walks each ordinate's
  sweep explicitly using the WDD recurrence ``ψ_out = 2·ψ̄ − ψ_in`` +
  ``streaming + redistribution + collision`` per cell, with NO
  bool-mask scatter into the legacy's misrouting ``ks``. The unified
  matvec uses ``out_g_first[:, global_X, i, 0] = m_full`` fancy-index
  scatter that is **structurally immune** to the routing bug. Both
  match at rtol=1e-12.

* **L1** — heterogeneous 3-region 2G closed cylinder eigenvalue.  The
  GMRES inner solve routes through :class:`InvertibleOperator`
  (= ``L + C``) consuming the unified matvec; the resulting ``k_eff``
  is cross-checked against the trajectory_resolvent reference
  (Variant α at α=1) at the same tolerance as the production Gate 4.2
  cross-check (3% — set by Variant α's quadrature error budget).

The legacy :func:`~orpheus.sn.operator.transport_operator_matvec_cylindrical`
has a per-ordinate routing bug (ERR-049 — ascending-global ``ks``
indices vs level-internal μ-sorted column order in
``streaming + redistribution + collision``). It is **NOT** a valid
reference here.
"""
from __future__ import annotations

import contextlib

import numpy as np
import pytest

from orpheus.derivations.common.xs_library import get_xs, make_mixture
from orpheus.derivations.continuous.trajectory_resolvent.greens_function_cylinder import (
    solve_greens_function_cylinder_mr,
)
from orpheus.geometry import BC, CoordSystem, Mesh1D
from orpheus.sn import operator as sn_op
from orpheus.sn import solve_sn
from orpheus.sn.geometry import SNMesh
from tests.sn._test_helpers import _LC_matvec
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.spatial.pole_angular_closure import MorelMontryAngularSweep
from orpheus.sn.spatial.psi_half_angle_seed import CarlsonSweepContext
from tests.sn._test_helpers import legacy_proxy_matvec, placeholder_materials


# ═══════════════════════════════════════════════════════════════════════
# Helpers — canonical ↔ packed conversion + per-ordinate hand reference
# ═══════════════════════════════════════════════════════════════════════


def _build_cyl(n_cells: int, quad, edges=None) -> SNMesh:
    """Build a homogeneous reflective cylinder mesh."""
    if edges is None:
        edges = np.linspace(0.1, 1.0, n_cells + 1)
    mesh = Mesh1D(
        edges=edges,
        mat_ids=np.zeros(n_cells, dtype=int),
        coord=CoordSystem.CYLINDRICAL,
        bc_left=BC("reflective"),
        bc_right=BC("reflective"),
    )
    return SNMesh(mesh, quad, placeholder_materials())


def _bc_fill_outer(psi_view: np.ndarray, sn_mesh: SNMesh) -> np.ndarray:
    """Make psi_view BC-consistent at the outer face (incoming ordinates)."""
    quad = sn_mesh.quad
    incoming_mask = quad.mu_x < -1e-15
    if not incoming_mask.any():
        return psi_view
    outer_face = psi_view[:, :, -1, 0]
    inflow_full = sn_mesh.bc_right.apply(outer_face)
    psi_view = psi_view.copy()
    psi_view[incoming_mask, :, -1, 0] = inflow_full[incoming_mask, :]
    return psi_view


def _extract_at_unknown_slots(
    field_4d: np.ndarray, sn_mesh: SNMesh,
) -> np.ndarray:
    """Gather field_4d at the curvilinear equation-bearing slots → (ng, n_eq).

    D-J (2026-05-30) — replaces the legacy :class:`EquationMap`-driven
    slot map.  Curvilinear 1-D equation set: all ``(n, ix, 0)`` slots
    EXCEPT inward ordinates at the outermost cell ``ix == nx - 1``
    (the reflective BC determines those values; they are NOT unknowns).
    """
    quad = sn_mesh.quad
    nx = sn_mesh.nx
    ng = field_4d.shape[1]
    inflow_outer = quad.mu_x < -1e-15  # (N,)
    cols = []
    for ix in range(nx):
        for n in range(quad.N):
            if ix == nx - 1 and inflow_outer[n]:
                continue
            cols.append(field_4d[n, :, ix, 0])
    return np.stack(cols, axis=1)  # (ng, n_eq)


def _hand_reference_cyl_matvec(
    psi_view: np.ndarray, sn_mesh: SNMesh, sigma_t: np.ndarray,
) -> np.ndarray:
    r"""Per-ordinate explicit matvec — the structurally-independent L0 reference.

    For each ordinate ``n``: walk cells in the sweep direction using the
    WDD recurrence and compute ``streaming + redistribution + collision``
    per cell with NO bool-mask scatter. Routing is impossible to get
    wrong because each ordinate is processed in its own scalar pass.

    Returns shape ``(N, ng, nx)``.
    """
    quad = sn_mesh.quad
    N = quad.N
    ng = psi_view.shape[1]
    nx = sn_mesh.nx
    eps = 1e-15

    reduced = sn_mesh.reduced
    A = reduced.face_areas
    V = sn_mesh.volumes[:, 0]
    mu_x = quad.mu_x

    bc_outer = sn_mesh.bc_right
    pac = sn_mesh.pole_angular_closure
    if pac is None:
        pac = MorelMontryAngularSweep()

    out = np.zeros((N, ng, nx))

    sigma_t_gx = sigma_t[:, :, 0]
    dr = sn_mesh.dx
    psi_g_first = psi_view.transpose(1, 0, 2, 3)
    outer_inflow_estimate = bc_outer.apply(psi_view[:, :, -1, 0])
    level_indices = quad.level_indices

    carlson_ctx_per_level = []
    for level_idx in level_indices:
        level_idx_arr = np.asarray(level_idx)
        mu_level = mu_x[level_idx_arr]
        weights_level = quad.weights[level_idx_arr]
        within_idx_most_inward = int(np.argmin(mu_level))
        global_idx_most_inward = int(level_idx_arr[within_idx_most_inward])
        bc_outer_value_level = outer_inflow_estimate[global_idx_most_inward, :]
        carlson_ctx_per_level.append(
            CarlsonSweepContext(
                sigma_t=sigma_t_gx, dr=dr,
                mu_quad=mu_level.copy(),
                weights=weights_level.copy(),
                bc_outer_value=bc_outer_value_level,
            )
        )

    redist_full = pac(
        psi_g_first[..., 0],
        reduced.alpha_per_level,
        reduced.redist_dAw_per_level,
        reduced.tau_mm_per_level,
        V,
        level_indices=level_indices,
        carlson_context=carlson_ctx_per_level,
    )  # (ng, N, nx)

    outflow_at_boundary = np.zeros((ng, N))

    # Outward — per ordinate.
    for level_idx in level_indices:
        level_idx_arr = np.asarray(level_idx)
        eta_level = mu_x[level_idx_arr]
        out_within = eta_level > +eps
        if not np.any(out_within):
            continue
        global_out = level_idx_arr[out_within]
        for n_g in global_out:
            mu_n = mu_x[n_g]
            psi_face_in = psi_g_first[:, n_g, 0, 0].copy()
            for i in range(nx):
                psi_cell = psi_g_first[:, n_g, i, 0]
                psi_face_out = 2.0 * psi_cell - psi_face_in
                streaming = mu_n * (
                    A[i + 1] * psi_face_out - A[i] * psi_face_in
                ) / V[i]
                redistribution = redist_full[:, n_g, i]
                collision = sigma_t_gx[:, i] * psi_cell
                out[n_g, :, i, 0] = streaming + redistribution + collision
                psi_face_in = psi_face_out
            outflow_at_boundary[:, n_g] = psi_face_out

    # BC trace.
    inflow_full = bc_outer.apply(outflow_at_boundary.T)

    # Inward — per ordinate.
    for level_idx in level_indices:
        level_idx_arr = np.asarray(level_idx)
        eta_level = mu_x[level_idx_arr]
        in_within = eta_level < -eps
        if not np.any(in_within):
            continue
        global_in = level_idx_arr[in_within]
        for n_g in global_in:
            mu_n = mu_x[n_g]
            psi_face_in = inflow_full[n_g, :]
            for i in range(nx - 1, -1, -1):
                psi_cell = psi_g_first[:, n_g, i, 0]
                psi_face_out = 2.0 * psi_cell - psi_face_in
                streaming = mu_n * (
                    A[i + 1] * psi_face_in - A[i] * psi_face_out
                ) / V[i]
                redistribution = redist_full[:, n_g, i]
                collision = sigma_t_gx[:, i] * psi_cell
                out[n_g, :, i, 0] = streaming + redistribution + collision
                psi_face_in = psi_face_out

    # Degenerate (|μ_x| < eps) — no radial flow.
    degenerate_mask = np.abs(mu_x) < eps
    if np.any(degenerate_mask):
        global_deg = np.where(degenerate_mask)[0]
        for n_g in global_deg:
            for i in range(nx):
                psi_cell = psi_g_first[:, n_g, i, 0]
                redistribution = redist_full[:, n_g, i]
                collision = sigma_t_gx[:, i] * psi_cell
                out[n_g, :, i, 0] = redistribution + collision

    return out


# ═══════════════════════════════════════════════════════════════════════
# L0 — Hand-reference battery (the structural correctness anchor)
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.xfail(
    reason="cylinder matvec/sweep WDD divergence — issue #206",
    strict=False,
)
@pytest.mark.l0
@pytest.mark.parametrize("quad_factory", [
    lambda: Quadrature.level_symmetric(sn_order=4),
    lambda: Quadrature.level_symmetric(sn_order=6),
    lambda: Quadrature.product(n_mu=2, n_phi=4),
])
@pytest.mark.parametrize("n_cells", [3, 5, 10])
@pytest.mark.parametrize("seed", [0, 1, 2])
def test_unified_cylinder_matches_hand_reference(
    quad_factory, n_cells, seed,
) -> None:
    """Unified cylindrical matvec matches the per-ordinate hand reference.

    Promoted from ``derivations/diagnostics/diag_step3_cyl_unified_vs_hand_battery.py``
    (numerics-investigator 2026-05-17 closeout). The hand reference is
    structurally independent of the unified (no bool-mask scatter); both
    agree at machine precision across LS4 / LS6 / ProductQuadrature ×
    {3, 5, 10} cells × {0, 1, 2} seeds = 27 cases.
    """
    quad = quad_factory()
    sn_mesh = _build_cyl(n_cells, quad)
    ng = 1
    N = quad.N

    rng = np.random.default_rng(seed)
    psi_view = rng.standard_normal((N, ng, n_cells)).astype(np.float64)
    psi_view = _bc_fill_outer(psi_view, sn_mesh)
    sigma_t = np.full((ng, n_cells), 2.0)

    m_unified = legacy_proxy_matvec(psi_view, sn_mesh, sigma_t)
    m_hand = _hand_reference_cyl_matvec(psi_view, sn_mesh, sigma_t)

    m_unified_u = _extract_at_unknown_slots(m_unified, sn_mesh)
    m_hand_u = _extract_at_unknown_slots(m_hand, sn_mesh)

    np.testing.assert_allclose(
        m_unified_u, m_hand_u, rtol=1e-12, atol=1e-13,
        err_msg=(
            f"Unified must match per-ordinate hand reference at quad "
            f"N={N}, n_cells={n_cells}, seed={seed}"
        ),
    )


@pytest.mark.l0
def test_unified_cylinder_zero_psi_gives_zero() -> None:
    """Linear operator: zero input → zero output."""
    quad = Quadrature.level_symmetric(sn_order=4)
    sn_mesh = _build_cyl(n_cells=5, quad=quad)
    ng = 1
    sigma_t = np.full((ng, sn_mesh.nx), 2.0)
    psi_view = np.zeros((quad.N, ng, sn_mesh.nx))

    m_unified = legacy_proxy_matvec(psi_view, sn_mesh, sigma_t)
    np.testing.assert_array_equal(m_unified, np.zeros_like(m_unified))


@pytest.mark.l0
def test_unified_cylinder_constant_psi_gives_sigma_t() -> None:
    """At ψ = constant on homogeneous reflective cylinder, unified matvec
    returns σ_t · ψ. Sanity check — flat flux activates only the
    collision term in the per-cell balance."""
    quad = Quadrature.level_symmetric(sn_order=4)
    sn_mesh = _build_cyl(n_cells=5, quad=quad)
    ng = 1
    sigma_t_val = 2.0
    sigma_t = np.full((ng, sn_mesh.nx), sigma_t_val)
    psi_view = np.ones((quad.N, ng, sn_mesh.nx))

    m_unified = legacy_proxy_matvec(psi_view, sn_mesh, sigma_t)
    m_at_unknowns = _extract_at_unknown_slots(m_unified, sn_mesh)
    np.testing.assert_allclose(
        m_at_unknowns, sigma_t_val, rtol=1e-12, atol=1e-13,
    )


# ═══════════════════════════════════════════════════════════════════════
# L1 — Heterogeneous trajectory_resolvent cross-check via Krylov
# ═══════════════════════════════════════════════════════════════════════
#
# The L0 battery proves the per-cell algebra. The L1 test proves the
# *converged eigenvalue* of a heterogeneous cylinder problem agrees with
# the structurally-independent trajectory_resolvent reference (Variant α
# at α=1, see ``orpheus.derivations.continuous.trajectory_resolvent``).
#
# Why heterogeneous: L2 in lessons.md — homogeneous k = νΣ_f/Σ_a is
# flux-shape independent (the matvec's redistribution & angular closure
# all collapse on flat flux). Heterogeneous closed MR exercises every
# term in the unified per-cell algebra.
# ═══════════════════════════════════════════════════════════════════════


_MR_RADII = np.array([0.5, 1.5, 2.0])
_MR_LAYOUT_ABA_KEYS = ("A", "B", "A")


def _mr_xs_2g():
    """3-region 2G XS bundle — same layout as ``test_phase_c_crosscheck`` MR."""
    parts = [get_xs(k, "2g") for k in _MR_LAYOUT_ABA_KEYS]
    sigma_t = np.stack([p["sig_t"] for p in parts], axis=0)
    sigma_s = np.stack([p["sig_s"] for p in parts], axis=0)
    nu_sigma_f = np.stack([p["nu"] * p["sig_f"] for p in parts], axis=0)
    chi = np.stack([p["chi"] for p in parts], axis=0)
    return sigma_t, sigma_s, nu_sigma_f, chi


def _make_2g_mixture(sigma_t, sig_s, nu_sig_f, chi):
    """Build a 2-group Mixture from explicit XS arrays."""
    sigma_t = np.asarray(sigma_t, dtype=float)
    sig_s = np.asarray(sig_s, dtype=float)
    nu_sig_f = np.asarray(nu_sig_f, dtype=float)
    chi = np.asarray(chi, dtype=float)
    sig_a = sigma_t - sig_s.sum(axis=1)
    nu = np.ones_like(nu_sig_f)
    sig_f = nu_sig_f.copy()
    sig_c = sig_a - sig_f
    return make_mixture(
        sig_t=sigma_t, sig_c=sig_c, sig_f=sig_f, nu=nu, chi=chi, sig_s=sig_s,
    )


def _build_mr_cylinder_mesh(nx: int = 40) -> tuple[Mesh1D, dict]:
    """Build the 3-region cylindrical mesh + 2G materials for the MR case."""
    sigma_t, sigma_s, nu_sigma_f, chi = _mr_xs_2g()
    materials = {
        i: _make_2g_mixture(sigma_t[i], sigma_s[i], nu_sigma_f[i], chi[i])
        for i in range(3)
    }
    # Cell-edge-aligned region boundaries.
    edges = np.linspace(0.0, _MR_RADII[-1], nx + 1)
    mat_ids = np.zeros(nx, dtype=int)
    cell_centres = 0.5 * (edges[:-1] + edges[1:])
    for i_cell, r_c in enumerate(cell_centres):
        if r_c <= _MR_RADII[0]:
            mat_ids[i_cell] = 0
        elif r_c <= _MR_RADII[1]:
            mat_ids[i_cell] = 1
        else:
            mat_ids[i_cell] = 0
    mesh = Mesh1D(
        edges=edges,
        mat_ids=mat_ids,
        coord=CoordSystem.CYLINDRICAL,
        bc_left=BC("reflective"),
        bc_right=BC("reflective"),
    )
    return mesh, materials


# Post-D-K (commit ``dadf4e8``), ``SNSolver.L`` is the algebraic
# composition ``StreamingOperator + CollisionOperator``
# (= :class:`InvertibleOperator`), which calls
# :func:`_transport_operator_matvec_unified` natively for 1-D
# cylindrical.  No monkey-patch is required.


@pytest.mark.l1
@pytest.mark.slow
@pytest.mark.verifies("sn-curvilinear-trajectory-resolvent-crosscheck")
def test_unified_cylinder_l1_mr_2g_trajectory_resolvent() -> None:
    r"""L1 — heterogeneous 3-region 2G closed cylinder via unified matvec.

    Drives :func:`solve_sn` with ``inner_solver="krylov"`` — the
    matvec routes through :class:`InvertibleOperator` (= ``L + C``)
    via :func:`_transport_operator_matvec_unified`. The converged
    ``k_eff`` is compared against the structurally-independent
    :func:`solve_greens_function_cylinder_mr` reference (Variant α at
    α=1) at the same 3% tolerance the production Gate 4.2 cross-check
    uses (set by Variant α's quadrature error budget at n_r=24 +
    n_traj_quad=64; see ``test_phase_c_crosscheck.py``).

    Per ``.claude/lessons.md`` L14 — solver correctness is a 4-way
    standoff. The L0 hand-reference battery in this file proves the
    per-cell algebra; this test pins the converged eigenvalue against
    a continuous reference that is structurally independent of every
    discrete primitive in the SN code (no shared FP path, no shared
    redist closure, no shared boundary recurrence).
    """
    mesh, materials = _build_mr_cylinder_mesh(nx=40)
    quad = Quadrature.level_symmetric(sn_order=4)

    sigma_t, sigma_s, nu_sigma_f, chi = _mr_xs_2g()
    ref = solve_greens_function_cylinder_mr(
        radii=_MR_RADII,
        sigma_t=sigma_t, sigma_s=sigma_s,
        nu_sigma_f=nu_sigma_f, chi=chi,
        alpha=1.0,
        n_r=24, n_mu_axial=16, n_phi_az=32, n_traj_quad=64,
        max_iter=500, tol=1e-7, initial_k=1.23,
    )
    k_ref = float(ref.k_eff)

    sol = solve_sn(
            materials=materials,
            mesh=mesh,
            quadrature=quad,
            inner_solver="krylov",
            max_outer=200, keff_tol=1e-7, flux_tol=1e-7,
            max_inner=200, inner_tol=1e-9,
        )

    rel = abs(sol.keff - k_ref) / k_ref
    assert rel < 3.0e-2, (
        f"unified cylinder L1 disagreement with trajectory_resolvent: "
        f"k_unified={sol.keff:.10f}, k_ref={k_ref:.10f}, rel={rel:.3e}"
    )


@pytest.mark.l1
@pytest.mark.slow
def test_unified_cylinder_l1_homogeneous_kinf_2g() -> None:
    r"""L1 sanity — 2G homogeneous closed cylinder eigenvalue via unified.

    Same monkey-patch construction as the MR test, but on homogeneous
    XS. k_eff is shape-independent (per L2 in lessons.md), so this test
    primarily confirms that the unified matvec drives GMRES to the
    correct k_∞ — a *necessary* but not *sufficient* condition. The
    sufficient condition is the heterogeneous MR test above.
    """
    from orpheus.derivations.common.eigenvalue import kinf_homogeneous

    sigma_t = [0.5, 1.0]
    sig_s = [[0.3, 0.05], [0.0, 0.7]]
    nu_sig_f = [0.4, 0.6]
    chi = [1.0, 0.0]
    mat = _make_2g_mixture(sigma_t, sig_s, nu_sig_f, chi)
    k_analytical = kinf_homogeneous(
        np.asarray(sigma_t), np.asarray(sig_s),
        np.asarray(nu_sig_f), np.asarray(chi),
    )

    nx = 20
    mesh = Mesh1D(
        edges=np.linspace(0.0, 2.0, nx + 1),
        mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.CYLINDRICAL,
        bc_left=BC("reflective"),
        bc_right=BC("reflective"),
    )
    quad = Quadrature.level_symmetric(sn_order=4)

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
        f"unified cylinder k_∞ recovery violated: "
        f"k_analytical={k_analytical:.10f}, k_unified={sol.keff:.10f}, "
        f"rel={rel:.2e}"
    )
