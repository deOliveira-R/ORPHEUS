r"""L1 MMS tests for SN — 1D spherical and cylindrical geometries.

Phase 3.3 (spherical) and Phase 3.4 (cylindrical) of the verification
campaign.  Both use isotropic-in-angle ansatz :math:`\psi_n(r) = A(r)/W`
with :math:`A(r) = \sin(\pi r/R)`, so the angular redistribution terms
vanish and the spatial DD convergence rate is isolated.

Closure history (Issue #195, ERR-026 → ERR-058): these tests carried
``xfail(strict=True)`` from Wave E Round 3 (2026-05) through the
ERR-058 fix (2026-06-12).  The terminal root causes were two
self-referential closure seeds in the curvilinear within-group
operator, each exact on flat ψ (hence invisible to every flat-flux
gate, vv-principles Mode 7) and O(1)-wrong on non-flat fields:

* **ERR-058 a (spatial)** — the +μ pole-face seed read the innermost
  CELL-CENTRE flux as the r=0 face value (O(h) wrong); fixed by the
  Carlson coupled-pole seed ψ(0, +μ) = ψ(0, −μ) (the inward chain's
  pole outflow at the mirror ordinate).
* **ERR-058 b (angular)** — the M-M half-angle thread's Carlson seed
  solved the starting-direction ODE with the proxy source
  Σ_t·φ₀/Σw (exact only at flat-flux equilibrium); fixed by the
  ``AngularEdgeExtrapolation`` seed (the input field extrapolated in
  μ to the level's angular edge).  The proxy's per-ordinate O(1)
  residual was invisible to every scalar check because the α-dome
  telescopes under the angular weight sum (anti-pattern #8).

Post-fix measured ladders (2026-06-12, SI ≡ Krylov bit-identical):
sphere ``[1.49e-2, 3.73e-3, 9.28e-4, 2.31e-4, 5.74e-5]`` at
nx ∈ {20..320} (orders 2.00–2.01); cylinder ``[2.16e-3, 5.39e-4,
1.35e-4, 3.37e-5]`` at nx ∈ {20..160} (orders 2.00).  Both
assertions (rate AND magnitude band) now hold as plain tests.

See :doc:`/theory/discrete_ordinates` (curvilinear MMS sections)
for the full derivation, and ``error_catalog.md`` ERR-058.
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
@pytest.mark.slow
@pytest.mark.catches("ERR-058")
@pytest.mark.verifies(
    "transport-spherical",
    "sn-mms-spherical-psi", "sn-mms-spherical-qext",
)
def test_sn_spherical_mms_converges_second_order():
    r"""Spherical SN with isotropic ansatz shows :math:`\mathcal{O}(h^2)`.

    The ansatz :math:`A(r) = \sin(\pi r/R)` vanishes at r=0 (symmetry)
    and r=R (vacuum).  Angular redistribution vanishes for isotropic
    flux, so only the radial DD closure drives the convergence rate.

    xfail removed 2026-06-12 (Issue #195): the ERR-058 closure-seed
    fix (coupled-pole spatial seed + angular-edge-extrapolation
    half-angle seed) restored O(h²) collapse — measured errors
    ``[1.49e-2, 3.73e-3, 9.28e-4, 2.31e-4]`` at this ladder, orders
    2.00–2.01, magnitude band satisfied at nx=160.  This gate is the
    primary regression catcher for the ERR-058 closure-seed class
    (any seed that is exact-on-flat-ψ but O(1)-wrong on non-flat
    fields re-floors this ladder).
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
@pytest.mark.slow
@pytest.mark.catches("ERR-058")
@pytest.mark.verifies(
    "transport-cylindrical",
    "sn-mms-cylindrical-psi", "sn-mms-cylindrical-qext",
)
def test_sn_cylindrical_mms_converges_second_order():
    r"""Cylindrical SN with isotropic ansatz shows :math:`\mathcal{O}(h^2)`.

    Same structure as the spherical test but on a cylindrical mesh
    with Product quadrature (polar × azimuthal).  Azimuthal
    redistribution vanishes for isotropic flux.

    xfail removed 2026-06-12 (Issue #195 / ERR-058): measured errors
    ``[2.16e-3, 5.39e-4, 1.35e-4, 3.37e-5]`` at this ladder, orders
    2.00; magnitude band satisfied from nx=40.  See the sphere
    variant's docstring + the module docstring for the closure
    narrative.
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
