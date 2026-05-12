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
        "ERR-026 — same root cause as the anisotropic curvilinear MMS "
        "(see tests/sn/l1_analytical/test_mms_curvilinear_aniso_dd_convergence.py "
        "for the full closure narrative).  Wave H Phase A/B shipped the "
        "BC-aware FD operator infrastructure; Wave H Phase C "
        "(commits eae6f05..814a895) shipped the sweep-frame matvec "
        "rewrite + retired Phase A's BoundaryFaceFlux Protocol + "
        "honored the §16A.3 BC trace contract — but the spherical "
        "pole-face WDD initial condition interacts with the canonical "
        "Hébert §3.9.4 angular recurrence in a way that still breaks "
        "Gate 1.1 per-ordinate flat-flux invariant on sphere MMS "
        "(cylindrical Gate 1.1 PASSES; the cylindrical per-level "
        "α-dome telescoping happens to absorb the discrepancy).  "
        "Full ERR-026 closure on MMS depends on Issue #192 Phase D — "
        "the Carlson coupled-pole sweep where the outward-ordinate "
        "pole-face initial condition is determined by the inward-"
        "ordinate pole-face propagation (the continuous symmetry "
        "condition at r=0)."
    ),
)
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
        phi_num = result.scalar_flux[:, 0, 0]
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
        "ERR-026 — same root cause as the spherical isotropic MMS. "
        "Cylindrical Gate 1.1 PASSES under canonical M-M angular "
        "closure (the per-level α-dome telescoping absorbs the "
        "pole-recurrence discrepancy that breaks sphere); this test "
        "stays xfail because the underlying spherical-pole flux-"
        "shape ERR-026 bug is what governs the default flip decision "
        "for both geometries. Issue #192 Phase D ships the pole-face "
        "spatial closure refinement + the default flip + this test's "
        "marker removal as a single coordinated change. "
        "See tests/sn/l1_analytical/test_mms_curvilinear_aniso_dd_convergence.py "
        "for the full closure narrative."
    ),
)
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
        phi_num = result.scalar_flux[:, 0, 0]
        phi_ref = case.phi_exact(mesh.centers)
        errors.append(_l2_1d(phi_num, phi_ref, mesh.volumes))

    errors = np.asarray(errors)
    orders = np.log2(errors[:-1] / errors[1:])

    assert np.all(orders > 1.9), (
        f"Expected O(h^2), got orders={orders}, errors={errors}"
    )
    assert 1e-8 < errors[-1] < 1e-3
