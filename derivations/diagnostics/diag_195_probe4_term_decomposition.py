"""Diagnostic: Issue #195 Probe 4 — term-level decomposition of the curvilinear
unified matvec applied to the manufactured solution.

Created by numerics-investigator on 2026-06-12.

Probe 3 (``diag_195_probe3_residual_audit.py``) established that the production
discrete operator ``(L+C-S-B)`` applied to the EXACT isotropic manufactured
solution ψ_ref,n = A(r)/W leaves an O(1), mesh-INDEPENDENT residual
(res[pole] = -0.236, res[mid] = -0.158, res[outer] = +0.158, 88% of res² mass
in the bulk).  The L2 SOLUTION error PLATEAUS (orders → 0).

This probe answers WHICH TERM is wrong.  It reproduces the EXACT per-cell
matvec body of ``_MSpatialOperatorSum._compute_LpC`` (operator.py:343), captures
each balance term per cell, and cross-validates the reconstructed ``(L+C)·ψ_ref``
against the actual operator output (bit-identical).  Then it compares each
DISCRETE term to its CONTINUOUS counterpart at a MID-DOMAIN cell.

Matvec per-cell math (operator.py:458-478):
    m_full[i] = (denom·ψ_cell[i] − numer_upstream[i]) / V[i]
    denom        = 2|μ|·A_down + dA_w·c_out + Σ_t·V
    numer_upstream = |μ|·A_total·ψ_face_in + dA_w·c_in·ψ_ang_in
    ψ_face_out   = 2·ψ_cell − ψ_face_in            (DD spatial closure)
    ψ_ang_in     = precompute_psi_state(...) M-M half-angle flux

Term split (matvec units, [neutrons/cm²/s], i.e. /V):
  streaming   = (2|μ|·A_down·ψ_cell − |μ|·A_total·ψ_face_in) / V
              = |μ|·(A_out·ψ_out − A_in·ψ_in) / V     (conservative WDD)
  redist      = (dA_w·c_out·ψ_cell − dA_w·c_in·ψ_ang_in) / V
  collision   = Σ_t·ψ_cell
  source(q)   = (μ·A'(r) + (Σ_t−Σ_s)·A(r)) / W   (per ordinate; NO ×V here)

Continuous counterparts (per ordinate) on isotropic ψ = A(r)/W:
  streaming_noncons = μ·A'(r)/W                       (what the source emits)
  streaming_cons    = μ·(A'(r) + k·A(r)/r)/W          (∇·(μψ) divergence form)
  curvature         = μ·k·A(r)/r/W   (k=2 sphere, 1 cylinder)
  redist            = 0   (isotropic ψ has ∂_μψ = 0)
  collision         = Σ_t·A(r)/W
The CONSISTENCY identity (continuous): streaming_cons = streaming_noncons +
curvature, and the redistribution exactly removes the curvature so the net is
streaming_noncons.  The DISCRETE operator must reproduce:
  (streaming_disc − streaming_noncons) + redist_disc → 0   for isotropic ψ.
If it does not, the curvature/redistribution split is the inconsistent term.

If this catches a real bug, promote to tests/sn/verification/mms/.


HISTORICAL (2026-06-12): this probe pins the BUG-ERA operator reconstruction (pre-ERR-058 cell-centre pole seed + Carlson proxy-source half-angle seed) bit-identically.  After the ERR-058 fix the reconstruction no longer matches production — EXPECTED.  Kept as diagnosis provenance; tests are skipped.
"""
from __future__ import annotations

import numpy as np
import pytest

pytestmark = pytest.mark.skip(
    reason="historical bug-era operator reconstruction (pre-ERR-058); kept as #195 diagnosis provenance",
)

from orpheus.derivations.continuous.mms.sn import (
    build_spherical_mms_case,
    build_cylindrical_mms_case,
)
from orpheus.sn.solver import (
    SNSolver, _build_fixed_source_rhs, _as_sn_mesh,
)
from orpheus.transport.spatial.cell_balance import cell_balance_for_streaming
from orpheus.transport.fields.angular_flux import AngularFlux
from orpheus.transport.fields.angular_boundary_flux import AngularBoundaryFlux
from orpheus.transport.timed_full_field import TimedFullField


def _build(case, nc):
    mesh = case.build_mesh(nc)
    Q = case.external_source(mesh)
    sn_mesh = _as_sn_mesh(mesh, case.quadrature, case.materials, "vacuum", mat_map=None)
    solver = SNSolver(sn_mesh, inner_solver="source_iteration",
                      scattering_order=0, max_inner=2000, inner_tol=1e-13)
    q_ext = _build_fixed_source_rhs(Q, sn_mesh)
    return mesh, sn_mesh, solver, q_ext


def _psi_ref(case, sn_mesh, mesh):
    A = case.phi_exact(mesh.centers)
    sum_w = float(sn_mesh.quad.weights.sum())
    N, ng, nx = sn_mesh.quad.N, sn_mesh.ng, sn_mesh.nx
    vals = np.zeros((N, ng, nx))
    vals[:, 0, :] = (A / sum_w)[None, :]
    return AngularFlux.from_mesh(vals, sn_mesh)


def _decompose(case, nc):
    """Reproduce _compute_LpC per-cell, capturing each balance term.

    Returns per-(ordinate, cell) arrays of streaming / redist / collision in
    matvec units (/V), plus the reconstruction m_full for cross-validation.
    """
    mesh, sn_mesh, solver, q_ext = _build(case, nc)
    psi_ref = _psi_ref(case, sn_mesh, mesh)
    rhs = TimedFullField.zeros(bulk=AngularFlux, boundary=AngularBoundaryFlux, mesh=sn_mesh)
    psi_tff = TimedFullField(bulk=psi_ref.copy() if hasattr(psi_ref, "copy") else psi_ref,
                             boundary=rhs.boundary)

    quad = sn_mesh.quad
    N, ng, nx = quad.N, sn_mesh.ng, sn_mesh.nx
    mu_x = quad.mu_x
    sum_w = float(quad.weights.sum())
    eps = 1e-15

    pole_closure = sn_mesh.pole_angular_closure
    level_indices = pole_closure.level_indices
    A_face = sn_mesh.areas                      # (nx+1,) face areas
    V = sn_mesh.volumes                         # (nx,)
    sigma_t_gx = solver  # placeholder; pull from operator below
    # Pull σ_t the operator uses (within-group LC).
    # B.2d: the triple retired into build_within_group_system; this fused
    # 3-block probe reads the production surfaces directly. B = B_a alone is
    # bit-identical here: on vacuum cases the ray-corner B_b term is exactly
    # zero, and B_a pads the ray slot present-zero like the retired composite.
    from orpheus.sn.coupled_system import build_streaming_collision
    from orpheus.sn.operators.boundary import SNBoundaryOperator
    LC = build_streaming_collision(solver.sn_mesh, solver.mat_xs)
    S = solver.scattering_op
    B = SNBoundaryOperator(solver.sn_mesh)
    # The StreamingOperator carries sigma_t; reach it via M_spatial.
    L_op = LC  # OperatorSum (L+C); we need the StreamingOperator leaf's sigma_t
    sig_t_op = None
    for leaf in getattr(LC, "operators", [LC]):
        if hasattr(leaf, "sigma_t") and getattr(leaf, "sn_mesh", None) is sn_mesh:
            sig_t_op = leaf.sigma_t
            break
    if sig_t_op is None:
        sig_t_op = np.array([[case.sigma_t] * nx])     # (ng, nx)
    sigma_t_gx = sig_t_op

    psi_view = psi_ref.values                   # (N, ng, nx)
    psi_g_first = psi_view.swapaxes(0, 1)       # (ng, N, nx)

    boundary = psi_tff.boundary
    face_outer = boundary.face_view("xmax")
    has_inner = "xmin" in boundary.layout.faces
    face_inner = boundary.face_view("xmin") if has_inner else None

    if sn_mesh.curvature != "cartesian":
        pole_face_seed = psi_view[:, :, 0].copy()
    else:
        pole_face_seed = face_inner

    psi_state = pole_closure.precompute_psi_state(
        psi_view, sigma_t=sigma_t_gx, bc_outer_inflow_estimate=face_outer,
    )

    # Per-term accumulators in matvec units (/V): (ng, N, nx)
    stream_term = np.zeros((ng, N, nx))
    redist_term = np.zeros((ng, N, nx))
    coll_term = np.zeros((ng, N, nx))
    m_full_recon = np.zeros((ng, N, nx))
    # also capture the angular upstream flux ψ_ang_in per (ng, N, nx)
    psi_ang_in_arr = np.full((ng, N, nx), np.nan)
    psi_face_in_arr = np.full((ng, N, nx), np.nan)

    def sweep_dir(direction_sign, psi_face_in_init):
        for p, level_idx in enumerate(level_indices):
            level_idx_arr = np.asarray(level_idx)
            mu_level = mu_x[level_idx_arr]
            within_mask = (mu_level > +eps if direction_sign > 0 else mu_level < -eps)
            if not np.any(within_mask):
                continue
            global_dir = level_idx_arr[within_mask]
            abs_mu = np.abs(mu_x[global_dir])
            within_positions = np.where(within_mask)[0]
            cell_indices = list(sn_mesh.dag_walk_cell_indices(
                direction_sign=direction_sign, mu_level_idx=p,
            ))
            if not cell_indices:
                continue
            psi_face_in = psi_face_in_init[global_dir, :].T   # (ng, n_dir)
            for i in cell_indices:
                psi_cell = psi_g_first[:, global_dir, i]      # (ng, n_dir)
                a_denom, a_numer = pole_closure.cell_contribution(
                    psi_state, i, p, within_positions,
                )
                A_down = A_face[i + 1] if direction_sign > 0 else A_face[i]
                A_tot = A_face[i] + A_face[i + 1]
                denom, numer_upstream = cell_balance_for_streaming(
                    abs_mu=abs_mu, A_downstream=A_down, A_total=A_tot,
                    total_xs=sigma_t_gx[:, i], volume=V[i],
                    psi_face_in=psi_face_in,
                    angular_denom_term=a_denom,
                    angular_numer_upstream=a_numer,
                )
                # m_full = (denom·ψ − numer_upstream)/V   (this IS (L+C)·ψ)
                m_full = (denom * psi_cell - numer_upstream) / V[i]
                # Term split:
                stream = (2.0 * abs_mu[None, :] * A_down * psi_cell
                          - abs_mu[None, :] * A_tot * psi_face_in) / V[i]
                redist = (a_denom[None, :] * psi_cell - a_numer) / V[i]
                coll = sigma_t_gx[:, i][:, None] * psi_cell
                for k, gn in enumerate(global_dir):
                    stream_term[:, gn, i] = stream[:, k]
                    redist_term[:, gn, i] = redist[:, k]
                    coll_term[:, gn, i] = coll[:, k]
                    m_full_recon[:, gn, i] = m_full[:, k]
                    psi_face_in_arr[:, gn, i] = psi_face_in[:, k]
                    # back out ψ_ang_in from a_numer = dA_w·c_in·ψ_ang
                    dn = a_denom[k]  # = dA_w·c_out
                    # (use a_numer directly; report it)
                    psi_ang_in_arr[:, gn, i] = a_numer[:, k]
                psi_face_in = 2.0 * psi_cell - psi_face_in
        return None

    sweep_dir(+1, pole_face_seed)
    sweep_dir(-1, face_outer)

    # ── cross-validate m_full_recon against the actual operator output ────
    Apsi = LC.apply(psi_tff)
    m_actual = Apsi.bulk.values.swapaxes(0, 1)   # (ng, N, nx)
    # only compare cells/ordinates the sweep visited (skip degenerate)
    mask = ~np.isnan(psi_face_in_arr)
    recon_err = np.nanmax(np.abs(m_full_recon[mask] - m_actual[mask]))

    return dict(
        mesh=mesh, sn_mesh=sn_mesh, case=case,
        stream=stream_term, redist=redist_term, coll=coll_term,
        m_full=m_full_recon, m_actual=m_actual,
        recon_err=recon_err, psi_ang_numer=psi_ang_in_arr,
        psi_face_in=psi_face_in_arr,
        sum_w=sum_w, mu_x=mu_x, weights=quad.weights, V=V,
        q_ext=q_ext,
    )


def _continuous(case, mesh, mu_n, k_geom):
    r = mesh.centers
    A = case.phi_exact(r)
    Ap = case.dphi_exact(r)
    sum_w = 1.0  # caller divides; we return PER-ordinate density / sum_w below
    return r, A, Ap


def _report(case, nc, label, k_geom):
    d = _decompose(case, nc)
    mesh = d["mesh"]; mu_x = d["mu_x"]; sum_w = d["sum_w"]; V = d["V"]
    r = mesh.centers
    A = case.phi_exact(r); Ap = case.dphi_exact(r)
    nx = len(r)
    mid = nx // 2
    n_pos = int(np.argmax(mu_x)); n_neg = int(np.argmin(mu_x))
    n_mid = int(np.argmin(np.abs(mu_x)))

    print(f"\n{'='*82}\n{label}  nc={nc}   recon_err vs actual matvec = {d['recon_err']:.2e}")
    print(f"  mid cell={mid} r={r[mid]:.4f}  ordinates μ_pos={mu_x[n_pos]:+.4f}(n{n_pos})  "
          f"μ_neg={mu_x[n_neg]:+.4f}(n{n_neg})  μ≈0={mu_x[n_mid]:+.4f}(n{n_mid})")

    for tag, gn in [("μ>0", n_pos), ("μ<0", n_neg), ("μ≈0", n_mid)]:
        mu_n = mu_x[gn]
        stream_d = d["stream"][0, gn, mid]
        redist_d = d["redist"][0, gn, mid]
        coll_d = d["coll"][0, gn, mid]
        # ingested source per ordinate (matvec units, /W, no ×V — m_full is /V already)
        src = (mu_n * Ap[mid] + (case.sigma_t - case.sigma_s) * A[mid]) / sum_w
        # continuous counterparts (per ordinate, /W)
        stream_noncons = mu_n * Ap[mid] / sum_w
        curvature = mu_n * k_geom * A[mid] / r[mid] / sum_w
        stream_cons = stream_noncons + curvature
        coll_c = case.sigma_t * A[mid] / sum_w
        net_disc = stream_d + redist_d + coll_d - src
        print(f"\n  --- {tag} n={gn} μ={mu_n:+.5f} @cell {mid} ---")
        print(f"    {'TERM':<28}{'discrete':>15}{'continuous':>15}{'mismatch':>15}")
        print(f"    {'streaming':<28}{stream_d:>15.6e}{stream_noncons:>15.6e}{stream_d-stream_noncons:>15.6e}")
        print(f"    {'  (vs conservative form)':<28}{stream_d:>15.6e}{stream_cons:>15.6e}{stream_d-stream_cons:>15.6e}")
        print(f"    {'redistribution':<28}{redist_d:>15.6e}{0.0:>15.6e}{redist_d:>15.6e}")
        print(f"    {'collision':<28}{coll_d:>15.6e}{coll_c:>15.6e}{coll_d-coll_c:>15.6e}")
        print(f"    {'source (ingested q)':<28}{src:>15.6e}{src:>15.6e}{0.0:>15.6e}")
        print(f"    {'NET (L·ψ−q)/V per-ord':<28}{net_disc:>15.6e}{0.0:>15.6e}{net_disc:>15.6e}")
        # the curvature/redist cancellation diagnostic
        stream_extra = stream_d - stream_noncons
        cancel = stream_extra + redist_d
        print(f"    streaming EXTRA (disc−μA'/W) = {stream_extra:+.6e}   "
              f"continuous curvature μ·kA/r/W = {curvature:+.6e}")
        print(f"    redist discrete             = {redist_d:+.6e}")
        print(f"    CURVATURE_in_stream + REDIST (→0 for consistency) = {cancel:+.6e}")
    # angular-integrated scalar residual at mid (sanity vs probe-3)
    w = d["weights"]
    res_scalar_mid = sum(
        w[n] * (d["stream"][0, n, mid] + d["redist"][0, n, mid] + d["coll"][0, n, mid]
                - (mu_x[n] * Ap[mid] + (case.sigma_t - case.sigma_s) * A[mid]) / sum_w)
        for n in range(len(mu_x))
    )
    print(f"\n  scalar residual at mid (Σ_n w_n·(L·ψ−q)/V) = {res_scalar_mid:+.6e}  "
          f"(probe-3 reference ≈ -0.158 sphere / -0.x cyl)")
    return d


def test_probe4_sphere_term_table():
    """Sphere mid-cell term table: names the inconsistent term."""
    case = build_spherical_mms_case()
    d = _report(case, 160, "SPHERE isotropic", k_geom=2.0)
    assert d["recon_err"] < 1e-10, (
        f"term reconstruction does not match the actual matvec "
        f"(err {d['recon_err']:.2e}) — decomposition is unreliable."
    )


def test_probe4_cylinder_term_table():
    case = build_cylindrical_mms_case()
    d = _report(case, 160, "CYLINDER isotropic", k_geom=1.0)
    assert d["recon_err"] < 1e-10


if __name__ == "__main__":
    _report(build_spherical_mms_case(), 160, "SPHERE isotropic", k_geom=2.0)
    _report(build_spherical_mms_case(), 320, "SPHERE isotropic", k_geom=2.0)
    _report(build_cylindrical_mms_case(), 160, "CYLINDER isotropic", k_geom=1.0)
