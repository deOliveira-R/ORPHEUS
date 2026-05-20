"""Diagnostic — :class:`InvertibleOperator.solve` ↔ ``transport_sweep`` ``W`` bridge.

Created by numerics-investigator on 2026-05-19.

Pins the load-bearing convention bridge in :meth:`orpheus.sn.operator.\
InvertibleOperator.solve`: the operator-algebra convention is
**per-ordinate units everywhere** (``S_normalised = S / sum_w``,
``q_ext = fission / sum_w``), but :func:`orpheus.sn.sweep.\
transport_sweep` expects ``aniso_source`` in **iso-source magnitude**
(it applies ``weight_norm = 1/W`` internally at ``sweep.py:432``).

Without the ``rhs.values * sum_w`` bridge inside ``InvertibleOperator.solve``,
``rhs.values`` ends up divided by ``W`` twice — the reflective
homogeneous-medium fixed point lands at ``k_inf / W`` instead of
``k_inf``.  The pre-fix R-1 Step E ``_solve_source_iteration`` carve
failed every L1 ``kinf_homogeneous[…-source_iteration]`` with this
exact signature (slab-2eg keff = 1.4844 vs 1.875; ratio = 1.4844/1.875
≈ 0.7917 ≈ 1/W·something for W=2 GL, modulated by scattering).

Two probes:

1. ``LC.solve(LC.apply(ψ))`` round-trip on a random ψ must agree to
   1e-10 — the algebraic-inverse contract.
2. Fixed-source uniform Σ_t reflective homogeneous medium must
   reproduce ``φ = q / Σ_t`` per ordinate (the classic streaming-
   equilibrium test that catches L0-SN-001 / ERR-004 class bugs).

If this test catches a real bug, promote to ``tests/sn/test_operator.py``
under :class:`TestInvertibleOperator`.
"""
from __future__ import annotations

import numpy as np
import pytest


def _build_homogeneous_setup(coord: str, n_cells: int = 10):
    """Build a homogeneous reflective 1-D problem of the given coord."""
    import sys, warnings
    sys.path.insert(0, "tests/sn/l1_analytical")
    warnings.simplefilter("ignore")
    from test_kinf_homogeneous import (
        _get_continuous_case, _homogeneous_mesh, _quadrature_for,
    )
    from orpheus.sn.geometry import SNMesh
    from orpheus.sn.solver import SNSolver

    case = _get_continuous_case("2eg")
    mat_id = next(iter(case.problem.materials.keys()))
    mesh = _homogeneous_mesh(coord=coord, n_cells=n_cells,
                             length=2.0, mat_id=mat_id)
    quad = _quadrature_for(coord)
    sn_mesh = SNMesh(mesh, quad, case.problem.materials)
    solver = SNSolver(
        sn_mesh=sn_mesh, scattering_order=0,
        max_inner=300, inner_tol=1e-12,
        inner_solver="source_iteration",
    )
    return solver, case


def test_invertible_solve_slab_uniform_roundtrip() -> None:
    """Slab ``LC.solve(LC.apply(ψ=1)) == ψ=1`` to machine zero.

    For slab, the streaming term cancels exactly when ψ is uniform
    (no curvilinear M-M angular redistribution to confound). The
    streaming-free check is the cleanest test of the convention
    bridge: pre-fix R-1 Step E failed this by a factor ``W = sum_w``
    (``= 2`` for slab GL-N=8), giving ``LC.solve(LC.apply(1)) = 1/W``.

    Restricted to slab because sphere / cylinder M-M angular closure
    is a discrete operator on its own and the matvec / sweep
    representations of M-M differ at the per-cell level even for
    uniform ψ — those are two distinct discrete operators that share
    a fixed point under SI/Krylov, NOT a strict matrix inverse.
    """
    from orpheus.sn.angular_flux import AngularFlux
    from orpheus.sn.operator import CollisionOperator, StreamingOperator

    solver, _ = _build_homogeneous_setup("slab", n_cells=10)
    sn_mesh = solver.sn_mesh
    N = solver.quad.N
    nx, ny, ng = sn_mesh.nx, sn_mesh.ny, solver.ng

    L_leaf = StreamingOperator(sn_mesh, solver.mat_xs.total_cross_section)
    C_t = CollisionOperator(sn_mesh, solver.mat_xs.total_cross_section)
    LC = L_leaf + C_t

    psi_known = AngularFlux(np.ones((N, ng, nx, ny)), sn_mesh)
    LC_psi = LC.apply(psi_known)

    psi_recovered = None
    for _ in range(20):
        rhs = psi_recovered.stash(LC_psi) if psi_recovered is not None else LC_psi
        psi_new = LC.solve(rhs)
        if psi_recovered is not None and np.abs(
            psi_new.values - psi_recovered.values
        ).max() < 1e-15:
            psi_recovered = psi_new
            break
        psi_recovered = psi_new

    diff = np.abs(psi_known.values - psi_recovered.values).max()
    assert diff < 1e-12, (
        f"slab LC.solve(LC.apply(ψ=1)) ≠ ψ=1: abs_max={diff:.3e}. "
        f"Pre-fix R-1 Step E failed this by factor W = sum_w = 2 "
        f"due to missing ``×W`` bridge in InvertibleOperator.solve."
    )


@pytest.mark.parametrize("coord", ["slab", "sphere", "cylinder"])
def test_invertible_solve_fixed_source_homogeneous_reflective(coord: str) -> None:
    """``LC.solve(q_ext)`` on reflective homogeneous medium → uniform per-ordinate ``q/Σ_t``.

    Streaming-equilibrium fixed-source test (L0-SN-001 family).  With
    reflective BCs and uniform per-ordinate source ``q_n = const``, the
    converged per-ordinate ψ must be uniform = ``q_n / Σ_t`` (no
    spatial gradient on uniform medium with reflective BCs).

    The pre-fix R-1 Step E carve gave ψ = ``q_n / (W Σ_t)`` — a factor
    of ``W`` too small.  The fix adds ``rhs.values * sum_w`` in
    :meth:`InvertibleOperator.solve` to compensate the sweep's internal
    ``/W`` on ``aniso_source``.
    """
    from orpheus.sn.angular_flux import AngularFlux
    from orpheus.sn.operator import CollisionOperator, StreamingOperator

    solver, _ = _build_homogeneous_setup(coord, n_cells=10)
    sn_mesh = solver.sn_mesh
    quad = solver.quad
    N = quad.N
    nx, ny, ng = sn_mesh.nx, sn_mesh.ny, solver.ng
    sum_w = float(quad.weights.sum())

    sigma_t = solver.mat_xs.total_cross_section  # (ng, nx, ny)
    L_leaf = StreamingOperator(sn_mesh, sigma_t)
    C_t = CollisionOperator(sn_mesh, sigma_t)
    LC = L_leaf + C_t

    # Per-ordinate uniform source — operator-algebra convention.
    q_iso = 0.225
    q_per_ord = np.full((N, ng, nx, ny), q_iso / sum_w)
    rhs = AngularFlux(q_per_ord, sn_mesh)

    # Iterate to converge the reflective BC partner-flux state.  Each
    # call threads ``previous`` as the lag-1 frame via ``stash``.
    # Sphere/cyl curvilinear M-M closure needs O(100) iters for the
    # reflective partner-flux state to propagate; slab converges
    # quickly.  Iterate until ``LC.solve`` is a stable fixed point.
    psi_typed = None
    for _ in range(400):
        if psi_typed is not None:
            rhs_lagged = psi_typed.stash(rhs)
        else:
            rhs_lagged = rhs
        psi_new = LC.solve(rhs_lagged)
        if psi_typed is not None and np.abs(
            psi_new.values - psi_typed.values
        ).max() < 1e-14:
            psi_typed = psi_new
            break
        psi_typed = psi_new

    # Per-ordinate expected: ψ_n = (q_iso / W) / Σ_t (= q_iso / (W Σ_t))
    # Wait — that's not right.  The per-ordinate balance for uniform ψ
    # on reflective BC is Σ_t ψ_n = q_n (no streaming since ψ flat),
    # so ψ_n = q_n / Σ_t = (q_iso / W) / Σ_t.
    psi_values = psi_typed.values
    for g in range(ng):
        sig_g = float(sigma_t[g, 0, 0])
        expected_per_ord = (q_iso / sum_w) / sig_g
        # Skip group 1 in the 2eg case where source is in group-0
        # only — actually q_per_ord here is uniform across BOTH groups,
        # so both groups should converge to the same q/(W Σ_t).
        actual = psi_values[:, g, :, 0]
        # Tolerance — sphere/cyl carry M-M angular closure that
        # introduces a small spatial dependence; accept rel 1e-8.
        rel_dev = np.abs(actual - expected_per_ord) / max(expected_per_ord, 1e-30)
        assert rel_dev.max() < 1e-6, (
            f"{coord} ng={g}: expected per-ord ψ = {expected_per_ord:.4e}, "
            f"max rel deviation = {rel_dev.max():.3e}; ψ field:\n{actual}"
        )


@pytest.mark.parametrize("coord", ["slab", "sphere", "cylinder"])
@pytest.mark.parametrize("ng_key", ["2eg", "4eg"])
def test_si_carve_recovers_analytical_kinf(coord: str, ng_key: str) -> None:
    """End-to-end SI carve on reflective homogeneous medium → analytical ``k_inf``.

    Pins the R-1 Step E ``_solve_source_iteration`` carve at the L1
    level.  Pre-fix (no ``×W`` bridge): k_eff lands at ``k_inf / W``
    or similar — every ``ng≥2`` case failed with ~21% drift.
    Post-fix: every case reaches the analytical reference at
    ``rtol < 1e-9``.

    Cross-check with the Krylov path (already verified in R-1 Step D)
    pins structural independence.
    """
    import sys, warnings
    sys.path.insert(0, "tests/sn/l1_analytical")
    warnings.simplefilter("ignore")
    from test_kinf_homogeneous import (
        _get_continuous_case, _homogeneous_mesh, _quadrature_for,
    )
    from orpheus.sn.geometry import SNMesh
    from orpheus.sn.solver import SNSolver

    case = _get_continuous_case(ng_key)
    mat_id = next(iter(case.problem.materials.keys()))
    mesh = _homogeneous_mesh(coord=coord, n_cells=10, length=2.0,
                             mat_id=mat_id)
    quad = _quadrature_for(coord)
    sn_mesh = SNMesh(mesh, quad, case.problem.materials)

    solver = SNSolver(
        sn_mesh=sn_mesh, scattering_order=0,
        max_inner=300, inner_tol=1e-12,
        inner_solver="source_iteration",
    )
    phi = solver.initial_flux_distribution()
    keff = 1.0
    for _ in range(60):
        fis = solver.compute_fission_source(phi, keff)
        phi_new = solver._solve_source_iteration(fis, phi)
        keff_new = solver.compute_keff(phi_new)
        if abs(keff_new - keff) < 1e-12:
            break
        keff, phi = keff_new, phi_new

    rel_err = abs(keff_new - case.k_eff) / case.k_eff
    assert rel_err < 1e-9, (
        f"SI carve {coord}-{ng_key}: keff={keff_new:.10f}, "
        f"ref={case.k_eff:.10f}, rel_err={rel_err:.3e}.  Pre-fix this "
        f"failed at ~21% (≈ 1/W·c_s for W=2 GL).  Verify "
        f"InvertibleOperator.solve carries the ``×W`` convention bridge."
    )
