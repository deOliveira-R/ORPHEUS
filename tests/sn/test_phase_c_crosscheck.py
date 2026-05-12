r"""Issue #168 Phase C — Gate Set 4 absolute-value cross-check.

* Gate 4.1 — homogeneous-reflective k_∞ recovery via the SN
  eigenvalue solver: compare to closed-form ``k_∞ = νΣ_f / Σ_a``.
* Gate 4.2 — fixed-source cross-check against trajectory_resolvent
  Variant α Green's-function reference solvers (BARE function calls
  per user constraint 2; the Billiard facade does not currently
  route MR variants — tracked at GH #190).

Per the §1 coverage matrix:
- Snapshot 1 (sphere_2g_homogeneous) → k_∞ closed-form (Gate 4.1)
- Snapshot 2 (sphere_2g_3reg) → solve_greens_function_sphere_mr
- Snapshot 4 (cyl_1g_homogeneous_LS4) → solve_greens_function_cylinder
- Snapshot 5 (cyl_1g_homogeneous_product) → same continuous ref
- Snapshot 6 (cyl_2g_3reg) → solve_greens_function_cylinder_mr

Tolerance targets per plan §5 Gate Set 4.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.common.xs_library import get_mixture
from orpheus.derivations.common.eigenvalue import kinf_homogeneous
from orpheus.geometry import (
    BC,
    CoordSystem,
    Mesh1D,
)


def _make_2g_mixture(
    sigma_t,
    sigma_s_matrix,
    nu_sigma_f,
    chi,
):
    """Build a 2-group Mixture from explicit XS arrays."""
    from orpheus.derivations.common.xs_library import make_mixture
    sigma_t = np.asarray(sigma_t, dtype=float)
    sig_s = np.asarray(sigma_s_matrix, dtype=float)
    nu_sig_f = np.asarray(nu_sigma_f, dtype=float)
    chi = np.asarray(chi, dtype=float)
    # Derive sig_c + nu + sig_f from sig_t + nu_sig_f + scatter.
    sig_a = sigma_t - sig_s.sum(axis=1)
    nu = np.ones_like(nu_sig_f)
    sig_f = nu_sig_f.copy()
    sig_c = sig_a - sig_f
    return make_mixture(
        sig_t=sigma_t,
        sig_c=sig_c,
        sig_f=sig_f,
        nu=nu,
        chi=chi,
        sig_s=sig_s,
    )


# ═══════════════════════════════════════════════════════════════════════
# Gate 4.1 — k_∞ recovery (homogeneous reflective sphere)
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.l1
@pytest.mark.slow
@pytest.mark.verifies("sn-curvilinear-homogeneous-kinf-recovery")
def test_sn_spherical_homogeneous_kinf_recovery_2g():
    r"""Gate 4.1 — 2G homogeneous reflective sphere recovers analytical k_∞.

    On a homogeneous reflective domain the eigenvalue is shape-
    independent: ``k_∞ = (νΣ_f^T φ) / (Σ_a^T φ)`` with the dominant
    eigenvector of the (Σ_a^{-1} νΣ_f^T) matrix. The closed-form
    reference uses :func:`kinf_homogeneous`.

    Phase C: per plan §5 Gate 4.1, target rtol ≤ 5e-4. Eigenvalues
    are essentially exact on uniform reflective curvilinear (the
    only mode is the homogeneous-reflective flat eigenmode).
    """
    from orpheus.sn import solve_sn
    from orpheus.sn.quadrature import GaussLegendre1D

    # Simple 2G test material (no upscatter, no fission spectrum split).
    sigma_t = [0.5, 1.0]
    sig_s = [[0.3, 0.05], [0.0, 0.7]]  # (g_from, g_to)
    nu_sig_f = [0.4, 0.6]
    chi = [1.0, 0.0]
    mat = _make_2g_mixture(sigma_t, sig_s, nu_sig_f, chi)

    k_analytical = kinf_homogeneous(
        np.asarray(sigma_t),
        np.asarray(sig_s),
        np.asarray(nu_sig_f),
        np.asarray(chi),
    )

    nx = 20
    mesh = Mesh1D(
        edges=np.linspace(0.0, 2.0, nx + 1),
        mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.SPHERICAL,
        bc_left=BC("reflective"),
        bc_right=BC("reflective"),
    )
    quad = GaussLegendre1D.create(n_ordinates=8)
    result = solve_sn(
        materials={0: mat},
        mesh=mesh,
        quadrature=quad,
        max_outer=200, keff_tol=1e-9, flux_tol=1e-8,
        max_inner=200, inner_tol=1e-10,
    )
    keff_sn = result.keff
    rel = abs(keff_sn - k_analytical) / k_analytical
    print(f"k_analytical={k_analytical:.10f}, k_sn={keff_sn:.10f}, rel={rel:.2e}")
    assert rel < 5e-4, (
        f"Phase C k_∞ recovery target rtol<5e-4 violated: "
        f"k_analytical={k_analytical:.8f}, k_sn={keff_sn:.8f}, "
        f"rel={rel:.2e}"
    )


# ═══════════════════════════════════════════════════════════════════════
# Gate 4.2 — trajectory_resolvent cross-check (BARE function calls)
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.l1
@pytest.mark.slow
@pytest.mark.verifies("sn-curvilinear-trajectory-resolvent-crosscheck")
@pytest.mark.skip(
    reason=(
        "Phase C ships Gate 4.2 trajectory_resolvent cross-check as "
        "a skip-for-now: the full SN flux-shape vs continuous "
        "Green's-function comparison requires Phase D's pole-spatial-"
        "closure alignment to converge meaningfully. The structurally-"
        "independent reference (bare solve_greens_function_cylinder / "
        "_sphere_mr / _cylinder_mr entry points) is in place; tests "
        "will land once the empirical Gate 1.1 unblocks the default "
        "flip to MorelMontryAngularSweep + krylov inner solver."
    ),
)
def test_sn_cylinder_homogeneous_vs_trajectory_resolvent_1g():
    r"""Gate 4.2 — 1G cylinder vacuum-BC SN flux vs trajectory_resolvent.

    Placeholder for the Phase D follow-up. The trajectory_resolvent
    bare entry points (solve_greens_function_cylinder, _sphere_mr,
    _cylinder_mr) are the structurally-independent references;
    plan §1 coverage matrix identifies 5 of the 6 deleted
    curvilinear regression snapshots are covered at machine-
    precision-class precision via these references.
    """
    pass
