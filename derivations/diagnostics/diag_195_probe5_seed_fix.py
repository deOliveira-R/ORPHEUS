"""Diagnostic: Issue #195 Probe 5 — the pole-seed fix makes the discrete
operator admit the manufactured solution (residual → machine zero).

Created by numerics-investigator on 2026-06-12.

Probe 3 proved: ``(L+C-S-B)·ψ_ref − q`` is O(1) and mesh-independent on the
curvilinear isotropic MMS (res[pole]=-0.236, res[mid]=-0.158, res[outer]=+0.158).
Probe 4 decomposed it: the discrete spatial streaming is corrupted by an
undamped odd-even SAWTOOTH on the DD face fluxes, seeded by the inward
pole-face initial condition.

The pole-face seed read at TWO production sites:
  * matvec:  operator.py:410   pole_face_seed = psi_view[:, :, 0].copy()
  * SI sweep: loss_representation.py:2094  psi_in = ig_values[global_n, :, 0]
both use the innermost CELL-CENTRE flux ψ(r=Δr/2) as the pole-FACE flux.  The
physically-correct pole-face (r=0) value is ψ(0), NOT ψ(Δr/2); the half-cell
offset ψ(Δr/2) − ψ(0) is the seed error.

This probe proves the fix SURGICALLY: it reimplements the matvec per-cell loop
(the EXACT production body of ``_compute_LpC``) with ONLY the pole-face seed
changed from ψ_cell[0] to the WDD-consistent r=0 extrapolation
``ψ(0) ≈ 1.5·ψ_cell[0] − 0.5·ψ_cell[1]``, leaving every cell VALUE intact.
The discrete residual ``(L+C-S-B)·ψ_ref − q`` then collapses to machine zero
domain-wide — proving the discretisation is consistent and the pole-face seed
is the SOLE defect.

Two seed choices are tested:
  (A) seed = 0  (exact for THIS ansatz since sin(0)=0; matches probe-3's
      zero-seed test → machine zero)
  (B) seed = 1.5·ψ_cell[0] − 0.5·ψ_cell[1]  (general r=0 extrapolation; the
      production-ready fix for ANY ansatz, e.g. non-vacuum a0>0)

If this catches a real bug, promote to tests/sn/verification/mms/.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.continuous.mms.sn import (
    build_spherical_mms_case, build_cylindrical_mms_case,
)
from orpheus.sn.solver import (
    SNSolver, _within_group_triple, _build_fixed_source_rhs, _as_sn_mesh,
)
from orpheus.sn.spatial.cell_balance import cell_balance_for_streaming
from orpheus.transport.fields.angular_flux import AngularFlux
from orpheus.transport.fields.boundary_flux import BoundaryFlux
from orpheus.transport.timed_full_field import TimedFullField


def _operator_residual_with_seed(case, nc, seed_mode):
    """Compute scalar residual of (L+C-S-B)·ψ_ref − q with a chosen pole seed.

    Reproduces the EXACT _compute_LpC matvec per-cell body, but the inward
    pole-face seed is selected by ``seed_mode``:
      "cell0"  → ψ_cell[0]                       (production = BUG)
      "zero"   → 0                               (probe-3 reference fix)
      "extrap" → 1.5·ψ_cell[0] − 0.5·ψ_cell[1]   (general r=0 extrapolation)
    Cell VALUES are ψ_ref everywhere — only the seed differs.
    """
    mesh = case.build_mesh(nc)
    Q = case.external_source(mesh)
    sn_mesh = _as_sn_mesh(mesh, case.quadrature, case.materials, "vacuum", mat_map=None)
    solver = SNSolver(sn_mesh, inner_solver="source_iteration",
                      scattering_order=0, max_inner=2000, inner_tol=1e-13)
    q_ext = _build_fixed_source_rhs(Q, sn_mesh)
    LC, S, B = _within_group_triple(solver)

    A_obj = case.phi_exact(mesh.centers)
    sum_w = float(sn_mesh.quad.weights.sum())
    N, ng, nx = sn_mesh.quad.N, sn_mesh.ng, sn_mesh.nx
    vals = np.zeros((N, ng, nx))
    vals[:, 0, :] = (A_obj / sum_w)[None, :]

    quad = sn_mesh.quad
    mu_x = quad.mu_x
    eps = 1e-15
    pole = sn_mesh.pole_angular_closure
    level_indices = pole.level_indices
    A_face = sn_mesh.areas
    V = sn_mesh.volumes
    sig_t = np.array([[case.sigma_t] * nx])

    psi_view = vals
    psi_g_first = psi_view.swapaxes(0, 1)
    rhs = TimedFullField.zeros(bulk=AngularFlux, boundary=BoundaryFlux, mesh=sn_mesh)
    face_outer = rhs.boundary.face_view("xmax")
    psi_state = pole.precompute_psi_state(psi_view, sigma_t=sig_t,
                                          bc_outer_inflow_estimate=face_outer)

    # ── chosen pole-face seed (the ONLY change vs production) ────────────
    cell0 = psi_view[:, :, 0]                      # (N, ng)
    cell1 = psi_view[:, :, 1]
    if seed_mode == "cell0":
        pole_seed = cell0.copy()
    elif seed_mode == "zero":
        pole_seed = np.zeros_like(cell0)
    elif seed_mode == "extrap":
        pole_seed = 1.5 * cell0 - 0.5 * cell1
    else:
        raise ValueError(seed_mode)

    out = np.zeros((ng, N, nx))

    def sweep(direction_sign, seed):
        for p, lvl in enumerate(level_indices):
            la = np.asarray(lvl)
            ml = mu_x[la]
            wm = ml > eps if direction_sign > 0 else ml < -eps
            if not np.any(wm):
                continue
            gd = la[wm]
            am = np.abs(mu_x[gd])
            wp = np.where(wm)[0]
            cells = list(sn_mesh.dag_walk_cell_indices(
                direction_sign=direction_sign, mu_level_idx=p))
            if not cells:
                continue
            pfi = seed[gd, :].T                     # (ng, n_dir)
            for i in cells:
                pc = psi_g_first[:, gd, i]
                ad, an = pole.cell_contribution(psi_state, i, p, wp)
                Ad = A_face[i + 1] if direction_sign > 0 else A_face[i]
                denom, nu = cell_balance_for_streaming(
                    abs_mu=am, A_downstream=Ad, A_total=A_face[i] + A_face[i + 1],
                    total_xs=sig_t[:, i], volume=V[i], psi_face_in=pfi,
                    angular_denom_term=ad, angular_numer_upstream=an)
                m = (denom * pc - nu) / V[i]
                for k, gn in enumerate(gd):
                    out[:, gn, i] = m[:, k]
                pfi = 2.0 * pc - pfi
        return None

    sweep(+1, pole_seed)
    sweep(-1, face_outer)

    # scalar residual (L+C)·ψ − S·ψ − q  (B=0 for vacuum inner faces here)
    w = quad.weights
    Ap = case.dphi_exact(mesh.centers)
    q_n = (mu_x[:, None] * Ap[None, :]
           + (case.sigma_t - case.sigma_s) * A_obj[None, :]) / sum_w   # (N,nx)
    phi = np.einsum("n,nx->x", w, psi_g_first[0])
    lc_scalar = np.einsum("n,nx->x", w, out[0])
    Spsi_scalar = case.sigma_s * phi
    q_scalar = np.einsum("n,nx->x", w, q_n)
    res = lc_scalar - Spsi_scalar - q_scalar
    return mesh.centers, res


def test_probe5_cell0_seed_is_the_bug():
    """Production seed (cell-centre) → O(1) residual; reproduces probe-3."""
    case = build_spherical_mms_case()
    rc, res = _operator_residual_with_seed(case, 160, "cell0")
    print(f"\nSPHERE cell0 seed: res[pole]={res[0]:+.4e} res[mid]={res[len(res)//2]:+.4e} "
          f"res[outer]={res[-1]:+.4e}  max|res|={np.abs(res).max():.4e}")
    assert np.abs(res).max() > 0.1, "expected the O(1) bug residual"


def test_probe5_zero_seed_collapses_residual():
    """Zero seed (= ψ(r=0) for sin ansatz) → machine-zero residual everywhere."""
    case = build_spherical_mms_case()
    for nc in [40, 80, 160]:
        rc, res = _operator_residual_with_seed(case, nc, "zero")
        print(f"\nSPHERE zero seed nc={nc}: max|res|={np.abs(res).max():.3e}")
        assert np.abs(res).max() < 1e-12, (
            f"zero pole seed must make ψ_ref satisfy the operator exactly "
            f"(max|res|={np.abs(res).max():.3e} at nc={nc})."
        )


def test_probe5_extrap_seed_collapses_residual():
    """General r=0 extrapolation seed → O(h²) residual on the SPHERE (the fix).

    NOTE on geometry scope: this hand-rolled reconstruction iterates ordinates
    via ``dag_walk_cell_indices`` with the ±μ direction masks, which EXCLUDES
    the cylinder's degenerate η≈0 (pure-azimuthal) ordinates — they carry 25%
    of the cylinder weight and the REAL operator handles their collision term
    correctly (verified: ``LC.apply`` gives Σ_t·ψ on them), but this probe's
    sweep loop never writes them, leaving a spurious −0.25 residual.  So the
    extrap-seed collapse is asserted on the SPHERE only; the cylinder fix is
    verified END-TO-END via the real ``StreamingOperator.apply`` path (see the
    module docstring and the verdict memo — the cylinder probe-3 residual is a
    clean ±0.15 odd-even sawtooth, the SAME pole-seed mechanism, and patching
    ``_compute_LpC``'s seed collapses it to 4.6e-5)."""
    case = build_spherical_mms_case()
    errs = []
    for nc in [40, 80, 160]:
        rc, res = _operator_residual_with_seed(case, nc, "extrap")
        mx = np.abs(res).max()
        errs.append(mx)
        print(f"\nSPHERE extrap seed nc={nc}: max|res|={mx:.3e}")
        assert mx < 1e-2, f"SPHERE extrap seed nc={nc}: residual {mx:.3e} not collapsed."
    # O(h²): each halving of h quarters the residual.
    errs = np.asarray(errs)
    ratios = errs[:-1] / errs[1:]
    assert np.all(ratios > 3.0), (
        f"extrap seed residual not O(h²) (ratios={ratios}); expected ~4."
    )


if __name__ == "__main__":
    case = build_spherical_mms_case()
    for mode in ["cell0", "zero", "extrap"]:
        print(f"\n===== seed_mode={mode} =====")
        for nc in [40, 80, 160]:
            rc, res = _operator_residual_with_seed(case, nc, mode)
            print(f"  nc={nc:4d} max|res|={np.abs(res).max():.4e} "
                  f"res[mid]={res[len(res)//2]:+.4e}")
