r"""L1 MMS tests for SN — 1D spherical and cylindrical geometries.

Phase 3.3 (spherical) and Phase 3.4 (cylindrical) of the verification
campaign.  Both use isotropic-in-angle ansatz :math:`\psi_n(r) = A(r)/W`
with :math:`A(r) = \sin(\pi r/R)`, so the angular redistribution terms
vanish and the spatial DD convergence rate is isolated.

Both tests are tagged ``@pytest.mark.xfail(strict=True)`` through Wave
E Round 3 — see :file:`tests/sn/l1_analytical/test_mms_curvilinear_aniso_dd_convergence.py`
for the full closure narrative; the legacy isotropic tests share the
same root cause (curvilinear FD operator boundary face-flux is
first-order on non-constant solutions, plus ERR-026 sweep WDD wrong
fixed point).

See :doc:`/theory/discrete_ordinates` (curvilinear MMS sections)
for the full derivation.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.continuous.mms.sn import (
    build_spherical_mms_case,
    build_cylindrical_mms_case,
)
from orpheus.sn import solve_sn_fixed_source


def _l2_1d(phi_num: np.ndarray, phi_ref: np.ndarray, volumes: np.ndarray) -> float:
    """Volume-weighted L2 norm for 1D curvilinear meshes."""
    diff = phi_num - phi_ref
    return float(np.sqrt(np.sum(volumes * diff * diff)))


# ═══════════════════════════════════════════════════════════════════════
# Phase 3.3 — 1D Spherical MMS
# ═══════════════════════════════════════════════════════════════════════

@pytest.mark.l1
@pytest.mark.xfail(
    strict=True,
    reason=(
        "ERR-026 PARTIAL — Phase D (commits 9512459..c44fe9b, 2026-05-12) "
        "shipped the canonical Hébert §3.9.4 Carlson coupled-pole "
        "inward μ=−1 sweep as the M-M angular recurrence's half-angle "
        "seed.  Per-ordinate flat-flux residual on Gate 1.1 sphere MMS "
        "collapsed to machine precision (≤ 1e-15), AND the empirical "
        "spatial convergence rate on this MMS ansatz under "
        "Krylov+Carlson is O(h²) at nx ∈ {20, 40, 80} (orders [3.33, "
        "2.46]).  However, the absolute-magnitude check at line ~96 "
        "(`1e-8 < errors[-1] < 1e-3`) fails at nx=160 because the L2 "
        "error is pre-asymptotic — at nx=80 the error is ~0.11 (probe), "
        "extrapolating with order 2.46 gives ~0.02 at nx=160, falling "
        "below 1e-3 only at nx ≥ 320-640.  This is either a benign "
        "pre-asymptotic transient (fix: refine mesh or relax magnitude "
        "bound) OR a coefficient bug (fix: investigate the MMS source "
        "discretisation vs the operator).  Issue #195 tracks the "
        "investigation + marker removal."
    ),
)
@pytest.mark.slow
@pytest.mark.verifies(
    "transport-spherical",
    "sn-mms-spherical-psi", "sn-mms-spherical-qext",
)
def test_sn_spherical_mms_converges_second_order():
    r"""Spherical SN with isotropic ansatz shows :math:`\mathcal{O}(h^2)`.

    The ansatz :math:`A(r) = \sin(\pi r/R)` vanishes at r=0 (symmetry)
    and r=R (vacuum).  Angular redistribution vanishes for isotropic
    flux, so only the radial DD closure drives the convergence rate.
    """
    case = build_spherical_mms_case()

    n_cells = [20, 40, 80, 160]
    errors = []
    for nc in n_cells:
        mesh = case.build_mesh(nc)
        Q = case.external_source(mesh)
        result = solve_sn_fixed_source(
            case.materials, mesh, case.quadrature, Q,
            max_inner=500, inner_tol=1e-13,
        )
        phi_num = result.scalar_flux.values[0, :]  # PR-INDEX-5: g=0 radial slice
        phi_ref = case.phi_exact(mesh.centers)
        errors.append(_l2_1d(phi_num, phi_ref, mesh.volumes))

    errors = np.asarray(errors)
    orders = np.log2(errors[:-1] / errors[1:])

    assert np.all(orders > 1.9), (
        f"Expected O(h^2), got orders={orders}, errors={errors}"
    )
    assert 1e-8 < errors[-1] < 1e-3


# ═══════════════════════════════════════════════════════════════════════
# Phase 3.4 — 1D Cylindrical MMS
# ═══════════════════════════════════════════════════════════════════════

@pytest.mark.l1
@pytest.mark.xfail(
    strict=True,
    reason=(
        "ERR-026 PARTIAL — coupled to the spherical isotropic case "
        "(this file, sphere variant).  Phase D's Carlson coupled-pole "
        "seed restored the canonical convergence rate (O(h²) under "
        "Krylov+Carlson on the empirical probe), but the absolute-"
        "magnitude check at nx=160 still fails because the L2 error "
        "is pre-asymptotic.  Issue #195 tracks the investigation + "
        "marker removal.  Cylindrical Gate 1.1 PASSES under canonical "
        "M-M angular closure (the per-level α-dome telescoping "
        "absorbed the wrong-seed discrepancy pre-Phase-D; the canonical "
        "seed now ships for cylindrical too).  This test stays xfail "
        "coupled to the sphere's #195 outcome — the convergence-magnitude "
        "investigation applies to both geometries."
    ),
)
@pytest.mark.slow
@pytest.mark.verifies(
    "transport-cylindrical",
    "sn-mms-cylindrical-psi", "sn-mms-cylindrical-qext",
)
def test_sn_cylindrical_mms_converges_second_order():
    r"""Cylindrical SN with isotropic ansatz shows :math:`\mathcal{O}(h^2)`.

    Same structure as the spherical test but on a cylindrical mesh
    with Product quadrature (polar × azimuthal).  Azimuthal
    redistribution vanishes for isotropic flux.
    """
    case = build_cylindrical_mms_case()

    n_cells = [20, 40, 80, 160]
    errors = []
    for nc in n_cells:
        mesh = case.build_mesh(nc)
        Q = case.external_source(mesh)
        result = solve_sn_fixed_source(
            case.materials, mesh, case.quadrature, Q,
            max_inner=500, inner_tol=1e-13,
        )
        phi_num = result.scalar_flux.values[0, :]  # PR-INDEX-5: g=0 radial slice
        phi_ref = case.phi_exact(mesh.centers)
        errors.append(_l2_1d(phi_num, phi_ref, mesh.volumes))

    errors = np.asarray(errors)
    orders = np.log2(errors[:-1] / errors[1:])

    assert np.all(orders > 1.9), (
        f"Expected O(h^2), got orders={orders}, errors={errors}"
    )
    assert 1e-8 < errors[-1] < 1e-3
