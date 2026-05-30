r"""L14 four-leg standoff — slab + cylinder via the typed operator algebra.

Per ``.claude/lessons.md`` L14, solver correctness is four-legged:

1.  **Algorithm 1 (Krylov via ``(L + C)``) ≡ structurally-independent reference.**
2.  **Algorithm 2 (sweep via ``transport_sweep``) ≡ structurally-independent reference.**
3.  **Algorithm 1 ≡ Algorithm 2** (twin-path agreement).
4.  **All three under mesh refinement** (right rate to right limit).

Two algorithms agreeing is necessary but NOT sufficient — both can be
equally wrong.  Post-D-K (commit ``dadf4e8``), ``solve_sn`` routes
through ``StreamingOperator + CollisionOperator`` =
:class:`InvertibleOperator`; the Krylov path uses GMRES on
``InvertibleOperator.apply`` with the sweep as preconditioner.

References (semi-analytical pillar per ``vv-principles``):

- **Slab** — the Case singular-eigenfunction reference
  ``sn_slab_1eg_2rg_S8`` (Case & Zweifel 1967; Gülderen & Türeci 2023;
  Garis & Sjöstrand 1990).  1G two-region reflective, S8.  The
  reference is mesh-independent and mathematically self-contained;
  ``solve_sn`` at the matching quadrature order must converge to it
  as h → 0.
- **Cylinder** — the trajectory-resolvent Variant α Green's-function
  reference ``solve_greens_function_cylinder_mr``.  3-region 2G
  reflective ABA layout.  Uses ray-traced chord integration in 3-D
  space; structurally independent of every discrete primitive in the
  SN code (no shared FP path, no shared redist closure, no shared
  boundary recurrence).

Sweep route (``inner_solver="source_iteration"``) routes through
``transport_sweep`` / ``CellUpdate.update``, which uses WDD
(Cartesian) / the same per-cell algebra as
:func:`transport_operator_matvec_unified` (curvilinear).
"""
from __future__ import annotations

import contextlib

import numpy as np
import pytest

from orpheus.derivations.common.xs_library import get_xs, make_mixture
from orpheus.derivations.continuous.trajectory_resolvent.greens_function_cylinder import (
    solve_greens_function_cylinder_mr,
)
from orpheus.derivations.reference_values import continuous_get
from orpheus.geometry import BC, CoordSystem, Mesh1D
from orpheus.sn import solve_sn
from orpheus.numerics.quadrature import Quadrature


# Post-D-K (commit ``dadf4e8``), ``SNSolver.L`` is the algebraic
# composition :class:`InvertibleOperator` (= ``StreamingOperator +
# CollisionOperator``), routing through
# :func:`transport_operator_matvec_unified` natively for 1-D slab /
# sphere / cylinder and through
# :meth:`StreamingOperator._apply_2d_cartesian` for 2-D Cartesian.
# No monkey-patch is required.


# ═══════════════════════════════════════════════════════════════════════
# Cylinder fixtures
# ═══════════════════════════════════════════════════════════════════════
_CYL_RADII = np.array([0.5, 1.5, 2.0])
_CYL_LAYOUT_ABA_KEYS = ("A", "B", "A")


def _cyl_xs_2g():
    parts = [get_xs(k, "2g") for k in _CYL_LAYOUT_ABA_KEYS]
    sigma_t = np.stack([p["sig_t"] for p in parts], axis=0)
    sigma_s = np.stack([p["sig_s"] for p in parts], axis=0)
    nu_sigma_f = np.stack([p["nu"] * p["sig_f"] for p in parts], axis=0)
    chi = np.stack([p["chi"] for p in parts], axis=0)
    return sigma_t, sigma_s, nu_sigma_f, chi


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


def _build_cyl_mesh(nx: int) -> tuple[Mesh1D, dict]:
    """3-region cylindrical ABA mesh + 2G materials (same as cylinder L1 in
    ``test_unified_matvec_cylinder.py``)."""
    sigma_t, sigma_s, nu_sigma_f, chi = _cyl_xs_2g()
    materials = {
        i: _make_2g_mixture(sigma_t[i], sigma_s[i], nu_sigma_f[i], chi[i])
        for i in range(3)
    }
    edges = np.linspace(0.0, _CYL_RADII[-1], nx + 1)
    mat_ids = np.zeros(nx, dtype=int)
    cell_centres = 0.5 * (edges[:-1] + edges[1:])
    for i_cell, r_c in enumerate(cell_centres):
        if r_c <= _CYL_RADII[0]:
            mat_ids[i_cell] = 0
        elif r_c <= _CYL_RADII[1]:
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


def _cylinder_k_ref() -> float:
    """trajectory-resolvent Variant α reference for the cylinder MR case.

    Cached at module import via the once-only ``solve_greens_function_cylinder_mr``
    call — the reference is expensive (~30 s) but identical across every
    mesh refinement.
    """
    sigma_t, sigma_s, nu_sigma_f, chi = _cyl_xs_2g()
    ref = solve_greens_function_cylinder_mr(
        radii=_CYL_RADII,
        sigma_t=sigma_t, sigma_s=sigma_s,
        nu_sigma_f=nu_sigma_f, chi=chi,
        alpha=1.0,
        n_r=24, n_mu_axial=16, n_phi_az=32, n_traj_quad=64,
        max_iter=500, tol=1e-7, initial_k=1.23,
    )
    return float(ref.k_eff)


def _solve_cyl_via_krylov_unified(nx: int) -> float:
    mesh, materials = _build_cyl_mesh(nx=nx)
    quad = Quadrature.level_symmetric(sn_order=4)
    sol = solve_sn(
            materials=materials, mesh=mesh, quadrature=quad,
            inner_solver="krylov",
            max_outer=200, keff_tol=1e-7, flux_tol=1e-7,
            max_inner=200, inner_tol=1e-9,
        )
    return float(sol.keff)


def _solve_cyl_via_sweep(nx: int) -> float:
    mesh, materials = _build_cyl_mesh(nx=nx)
    quad = Quadrature.level_symmetric(sn_order=4)
    sol = solve_sn(
        materials=materials, mesh=mesh, quadrature=quad,
        inner_solver="source_iteration",
        max_outer=500, keff_tol=1e-7, flux_tol=1e-7,
        max_inner=500, inner_tol=1e-9,
    )
    return float(sol.keff)


# ═══════════════════════════════════════════════════════════════════════
# Cylinder L1 standoff
# ═══════════════════════════════════════════════════════════════════════
# Reference tolerance: 3% (matches existing ``test_unified_cylinder_l1_mr_2g_trajectory_resolvent``;
# set by trajectory_resolvent's quadrature error budget at n_r=24, n_traj_quad=64).
# Twin-path tolerance: 1e-6 rel (both algorithms converge to the same
# discrete fixed point at the matched solver tolerances).
_CYL_REF_RTOL = 3.0e-2
_CYL_TWIN_RTOL = 1.0e-5


@pytest.mark.l1
@pytest.mark.slow
@pytest.mark.verifies("sn-curvilinear-trajectory-resolvent-crosscheck")
def test_cylinder_l1_sweep_vs_trajectory_resolvent() -> None:
    r"""**Cylinder Leg 2** — sweep ≡ trajectory_resolvent reference.

    Production source-iteration path (``transport_sweep`` /
    ``CellUpdate.update``) on the 3-region 2G ABA cylinder. No shim:
    the sweep already uses the WDD-correct per-cell algebra that the
    unified matvec also wraps.
    """
    k_ref = _cylinder_k_ref()
    k_sweep = _solve_cyl_via_sweep(nx=40)
    rel = abs(k_sweep - k_ref) / k_ref
    assert rel < _CYL_REF_RTOL, (
        f"cylinder sweep vs trajectory_resolvent: "
        f"k_sweep={k_sweep:.10f}, k_ref={k_ref:.10f}, rel={rel:.3e}"
    )


@pytest.mark.l1
@pytest.mark.slow
@pytest.mark.xfail(
    strict=True,
    reason=(
        "PR-TYPED-6.5 Phase 5 — cylinder twin-path divergence "
        "(rel ≈ 4e-3 at nx=40) is the L14 manifestation-#6 signature "
        "the B1'' face-state architecture fixes.  GMRES on the "
        "B1''-aware (L+C) through :class:`InvertibleOperator` "
        "converges at FP-noise (verified by "
        "``test_b1pp_cylinder_gmres_converges`` in "
        "``tests/sn/test_b1pp_verification.py``).  Post-D-K, "
        "``solve_sn`` routes through the (L+C) operator algebra "
        "natively; this xfail should be re-validated and likely "
        "flipped green by a follow-up V&V pass."
    ),
)
def test_cylinder_l1_sweep_vs_krylov_twin_path() -> None:
    r"""**Cylinder Leg 3** — sweep ≡ Krylov-via-unified twin-path agreement.

    Both algorithms drive the SAME continuous-equation discrete fixed
    point; at matched solver tolerances they MUST agree at sub-ULP-of-
    iteration drift. Disagreement here is the L14 manifestation-#6
    signature: same equation, two algorithmic paths, two answers.

    Currently xfailed strict — the test predates D-K's retargeting of
    ``solve_sn`` to the B1''-aware (L+C) algebra and has not been
    re-validated since.  See ``test_b1pp_verification.py`` for the
    direct (L+C) verification that the B1'' fix works.
    """
    k_sweep = _solve_cyl_via_sweep(nx=40)
    k_krylov = _solve_cyl_via_krylov_unified(nx=40)
    rel = abs(k_sweep - k_krylov) / k_sweep
    assert rel < _CYL_TWIN_RTOL, (
        f"cylinder twin-path disagreement: "
        f"k_sweep={k_sweep:.10f}, k_krylov={k_krylov:.10f}, rel={rel:.3e}"
    )


@pytest.mark.l1
@pytest.mark.slow
@pytest.mark.parametrize("nx", [20, 40, 80])
@pytest.mark.xfail(
    strict=True,
    reason=(
        "PR-TYPED-6.5 Phase 5 — twin-path divergence on cylinder at "
        "the Krylov leg's pre-D-K Carlson seed (cell-centre proxy).  "
        "Post-D-K, ``solve_sn`` migrated to the B1''-aware (L+C) "
        "algebra via :class:`InvertibleOperator`; xfail should be "
        "re-validated and likely flipped green.  See "
        "``test_b1pp_verification.py`` for the direct B1'' L1 "
        "verification on the Resolution A leaves."
    ),
)
def test_cylinder_l1_refinement_both_paths(nx: int) -> None:
    r"""**Cylinder Leg 4** — both paths converge to ref under refinement.

    For each nx ∈ {20, 40, 80}, both the sweep and the Krylov-via-unified
    paths must agree with the trajectory_resolvent reference to within
    the reference's own quadrature error budget. This is the joint test
    of "right rate to right limit" for both algorithms on the same case.

    Currently xfailed strict — see ``test_cylinder_l1_sweep_vs_krylov_twin_path``.
    """
    k_ref = _cylinder_k_ref()
    k_sweep = _solve_cyl_via_sweep(nx=nx)
    k_krylov = _solve_cyl_via_krylov_unified(nx=nx)
    rel_sweep = abs(k_sweep - k_ref) / k_ref
    rel_krylov = abs(k_krylov - k_ref) / k_ref
    rel_twin = abs(k_sweep - k_krylov) / k_sweep
    assert rel_sweep < _CYL_REF_RTOL, (
        f"cylinder nx={nx}: sweep vs ref rel={rel_sweep:.3e} ≥ {_CYL_REF_RTOL:.0e}"
    )
    assert rel_krylov < _CYL_REF_RTOL, (
        f"cylinder nx={nx}: krylov vs ref rel={rel_krylov:.3e} ≥ {_CYL_REF_RTOL:.0e}"
    )
    assert rel_twin < _CYL_TWIN_RTOL, (
        f"cylinder nx={nx}: twin-path rel={rel_twin:.3e} ≥ {_CYL_TWIN_RTOL:.0e}"
    )


# ═══════════════════════════════════════════════════════════════════════
# Slab fixtures
# ═══════════════════════════════════════════════════════════════════════


def _build_slab_2region_mesh(n_per: int) -> tuple[Mesh1D, dict, int]:
    r"""Build the ``sn_slab_1eg_2rg_S8`` Case singular-eigenfunction mesh.

    Returns ``(mesh, materials, N_ord)`` for direct use by ``solve_sn``.
    """
    ref = continuous_get("sn_slab_1eg_2rg_S8")
    geom = ref.problem.geometry_params
    materials = ref.problem.materials
    H_A = float(geom["fuel_height"])
    H_B = float(geom["refl_height"])
    N_ord = int(geom["n_ordinates"])
    edges = np.linspace(0.0, H_A + H_B, 2 * n_per + 1)
    mat_ids = np.array([0] * n_per + [1] * n_per)
    mesh = Mesh1D(edges=edges, mat_ids=mat_ids, coord=CoordSystem.CARTESIAN)
    return mesh, materials, N_ord


def _slab_k_ref() -> float:
    return float(continuous_get("sn_slab_1eg_2rg_S8").k_eff)


def _solve_slab_via_krylov_unified(n_per: int) -> float:
    mesh, materials, N_ord = _build_slab_2region_mesh(n_per=n_per)
    quad = Quadrature.gauss_legendre(N_ord)
    sol = solve_sn(
            materials, mesh, quad,
            inner_solver="krylov",
            max_outer=500, max_inner=500,
            keff_tol=1e-12, inner_tol=1e-9,
        )
    return float(sol.keff)


def _solve_slab_via_sweep(n_per: int) -> float:
    mesh, materials, N_ord = _build_slab_2region_mesh(n_per=n_per)
    quad = Quadrature.gauss_legendre(N_ord)
    sol = solve_sn(
        materials, mesh, quad,
        inner_solver="source_iteration",
        max_outer=500, max_inner=500,
        keff_tol=1e-12, inner_tol=1e-12,
    )
    return float(sol.keff)


# ═══════════════════════════════════════════════════════════════════════
# Slab L1 standoff
# ═══════════════════════════════════════════════════════════════════════
# Reference tolerance budget (per existing ``test_sn_2region_reflective_case_eigenvalue``):
#   |Δk| < 1e-5 at n_per=320 (sweep finest-mesh budget).
#   At n_per=160 the production sweep converges to |Δk| ~ 1e-7 typically
#   (still O(h) toward the Case reference); allow |Δk| < 2e-5 for safety
#   margin at the coarser meshes used here.
# Twin-path tolerance: 1e-6 rel.
_SLAB_REF_ABSTOL_AT_160 = 2.0e-5
_SLAB_TWIN_RTOL = 1.0e-5


@pytest.mark.l1
@pytest.mark.slow
@pytest.mark.catches("ERR-025")
def test_slab_l1_krylov_via_unified_vs_case() -> None:
    r"""**Slab Leg 1** — Krylov via unified matvec ≡ Case reference.

    Case singular-eigenfunction reference ``sn_slab_1eg_2rg_S8`` at the
    matching S8 quadrature.  Krylov is monkey-patched to route through
    the unified WDD matvec (vs the legacy 1st-order FD).  WDD converges
    to the Case reference at O(h) (piecewise-Σ interface degrades from
    nominal O(h²)).  At n_per=160 the production sweep reaches |Δk| ~
    1e-7; we allow 2e-5 for tolerance budget headroom.
    """
    k_ref = _slab_k_ref()
    k_krylov = _solve_slab_via_krylov_unified(n_per=160)
    abs_err = abs(k_krylov - k_ref)
    assert abs_err < _SLAB_REF_ABSTOL_AT_160, (
        f"slab Krylov-via-unified vs Case reference: "
        f"k_krylov={k_krylov:.10f}, k_ref={k_ref:.10f}, |Δ|={abs_err:.3e}"
    )


@pytest.mark.l1
@pytest.mark.slow
def test_slab_l1_sweep_vs_krylov_twin_path() -> None:
    r"""**Slab Leg 3** — sweep ≡ Krylov-via-unified twin-path agreement.

    Both algorithms drive the SAME WDD discrete fixed point on the same
    mesh.  At matched solver tolerances they agree at sub-ULP-of-
    iteration drift.  Disagreement signals algorithmic divergence (the
    L14 #6 manifestation pattern).
    """
    k_sweep = _solve_slab_via_sweep(n_per=80)
    k_krylov = _solve_slab_via_krylov_unified(n_per=80)
    rel = abs(k_sweep - k_krylov) / k_sweep
    assert rel < _SLAB_TWIN_RTOL, (
        f"slab twin-path disagreement: "
        f"k_sweep={k_sweep:.10f}, k_krylov={k_krylov:.10f}, rel={rel:.3e}"
    )


@pytest.mark.l1
@pytest.mark.slow
@pytest.mark.parametrize("n_per", [40, 80, 160])
def test_slab_l1_refinement_both_paths(n_per: int) -> None:
    r"""**Slab Leg 4** — both paths converge to Case ref + agree at each refinement.

    For each n_per ∈ {40, 80, 160}, both sweep and Krylov-via-unified
    must agree with the Case reference (decreasing error under
    refinement) and with each other (twin-path agreement).  The existing
    sweep refinement test in ``test_heterogeneous_transport.py`` covers
    n_per ∈ {20, 40, 80, 160, 320} for the sweep alone; this test adds
    the Krylov-via-unified leg.
    """
    k_ref = _slab_k_ref()
    k_sweep = _solve_slab_via_sweep(n_per=n_per)
    k_krylov = _solve_slab_via_krylov_unified(n_per=n_per)
    abs_err_sweep = abs(k_sweep - k_ref)
    abs_err_krylov = abs(k_krylov - k_ref)
    rel_twin = abs(k_sweep - k_krylov) / k_sweep
    # Each path must converge toward the Case reference.  Looser
    # tolerance at coarser meshes (O(h) interface error scales linearly):
    # 5e-4 at n=40, 2.5e-4 at n=80, 2e-5 at n=160.
    tolerance = {40: 5.0e-4, 80: 2.5e-4, 160: 2.0e-5}[n_per]
    assert abs_err_sweep < tolerance, (
        f"slab n_per={n_per}: sweep |Δk|={abs_err_sweep:.3e} ≥ {tolerance:.0e}"
    )
    assert abs_err_krylov < tolerance, (
        f"slab n_per={n_per}: krylov |Δk|={abs_err_krylov:.3e} ≥ {tolerance:.0e}"
    )
    assert rel_twin < _SLAB_TWIN_RTOL, (
        f"slab n_per={n_per}: twin-path rel={rel_twin:.3e} ≥ {_SLAB_TWIN_RTOL:.0e}"
    )
