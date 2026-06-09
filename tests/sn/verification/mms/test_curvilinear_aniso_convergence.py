r"""L1 MMS convergence gates for curvilinear SN with the **anisotropic** ansatz.

Consolidates two legacy files that exercised the SAME
``build_{cylindrical,spherical}_anisotropic_mms_case`` references via
the angularly-non-trivial :math:`(A(r) + B(r)\,\zeta)/W` ansatz (the
Mode-7 / vv-principles redesign that activates the angular-redistribution
term the isotropic ansatz nulls):

* ``tests/sn/test_phase_c_mms.py`` — Issue #168 Phase C Gate Set 3
  (spatial convergence on sphere + cylinder, ``catches("ERR-026")``,
  plus the Gate 3.3 angular-convergence-at-fixed-mesh test).
* ``tests/sn/l1_analytical/test_mms_curvilinear_aniso_dd_convergence.py``
  — the strict-xfail spatial-convergence variant with the
  ``sn-mms-{sph,cyl}-aniso-{psi,qext}`` equation labels + the absolute-
  magnitude band assertion.

Both variants are kept (they verify DIFFERENT equation labels and the
phase-C pair carries the ERR-026 catcher) — the consolidation removes
the duplicated FILE, not the distinct coverage. All stay ``xfail``:
ERR-026's curvilinear-sweep WDD fixed point is closed at the
per-ordinate flat-flux level (Phase D Carlson coupled-pole seed) but
the L2 convergence magnitude is pre-asymptotic at nx=80 (Issue #195).

Pairing rationale (the failure-narrowing instrument): the anisotropic
ansatz differs from the isotropic one ONLY in the :math:`B(r)\zeta`
term that activates :math:`(1-\mu^2)B/r` (sphere) / :math:`\xi^2 B/r`
(cylinder). If the isotropic companion
(``test_mms_curvilinear.py``) passes and these fail, the bug is in the
angular-redistribution path; if both fail alike, it is upstream in the
DD spatial closure. The pair IS the diagnostic.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.continuous.mms.sn import (
    build_cylindrical_anisotropic_mms_case,
    build_spherical_anisotropic_mms_case,
)
from orpheus.sn import solve_sn_fixed_source

pytestmark = pytest.mark.l1


def _l2_1d(phi_num: np.ndarray, phi_ref: np.ndarray, volumes: np.ndarray) -> float:
    """Volume-weighted L2 norm for 1D curvilinear meshes."""
    diff = phi_num - phi_ref
    return float(np.sqrt(np.sum(volumes * diff * diff)))


# ═══════════════════════════════════════════════════════════════════════
# Phase C Gate Set 3 — spatial + angular convergence (catches ERR-026)
# ═══════════════════════════════════════════════════════════════════════

@pytest.mark.slow
@pytest.mark.verifies("sn-mms-spherical-aniso-spatial-convergence")
@pytest.mark.catches("ERR-026")
@pytest.mark.xfail(
    strict=False,
    reason=(
        "ERR-026 PARTIAL CLOSURE (Phase C 2026-05-12, updated "
        "PR-TYPED-6.5 2026-05-18): the sweep-frame apply matvec "
        "rewrite + B1'' face-state architecture closed the spatial "
        "side of ERR-026, but the Gate 1.1 empirical probe with "
        "``MorelMontryAngularSweep`` on sphere still does not preserve "
        "the per-ordinate flat-flux invariant by design (α-dome "
        "cancellation across pure-azimuthal degenerate ordinates "
        "happens to telescope cleanly on cylinder; sphere does not).  "
        "The MMS convergence rate stays at the ~O(h^1.3) profile until "
        "a follow-up aligns the pole-face spatial closure with the "
        "canonical Hébert §3.9.4 angular recurrence."
    ),
)
def test_sn_spherical_aniso_mms_spatial_convergence_phase_c():
    r"""Gate 3.1 — spherical spatial MMS, anisotropic ansatz.

    Refines nx ∈ {10, 20, 40, 80}. Asserts ``min(orders[-2:]) > 1.9``.
    """
    case = build_spherical_anisotropic_mms_case()
    n_cells = [10, 20, 40, 80]
    errors = []
    for nc in n_cells:
        mesh = case.build_mesh(nc)
        Q = case.external_source(mesh)
        result = solve_sn_fixed_source(
            case.materials, mesh, case.quadrature, Q,
            max_inner=500, inner_tol=1e-13,
        )
        phi_num = result.scalar_flux.values[0, :]  # (ng=1, nx, ny=1)
        phi_ref = case.phi_exact(mesh.centers)
        errors.append(_l2_1d(phi_num, phi_ref, mesh.volumes))

    errors = np.asarray(errors)
    orders = np.log2(errors[:-1] / errors[1:])
    print(f"spherical_aniso errors: {errors}")
    print(f"spherical_aniso orders: {orders}")
    assert np.all(orders[-2:] > 1.9), (
        f"Expected O(h^2) in the last 2 orders; got orders={orders}, "
        f"errors={errors}"
    )


@pytest.mark.slow
@pytest.mark.verifies("sn-mms-cylindrical-aniso-spatial-convergence")
@pytest.mark.catches("ERR-026")
@pytest.mark.xfail(
    strict=False,
    reason=(
        "ERR-026 PARTIAL CLOSURE — same rationale as the spherical "
        "Gate 3.1 test; see its docstring. Cylindrical follows the "
        "same deferral path until the pole-face spatial-closure fix."
    ),
)
def test_sn_cylindrical_aniso_mms_spatial_convergence_phase_c():
    r"""Gate 3.2 — cylindrical spatial MMS, anisotropic ansatz."""
    case = build_cylindrical_anisotropic_mms_case()
    n_cells = [10, 20, 40, 80]
    errors = []
    for nc in n_cells:
        mesh = case.build_mesh(nc)
        Q = case.external_source(mesh)
        result = solve_sn_fixed_source(
            case.materials, mesh, case.quadrature, Q,
            max_inner=500, inner_tol=1e-13,
        )
        phi_num = result.scalar_flux.values[0, :]  # (ng=1, nx, ny=1)
        phi_ref = case.phi_exact(mesh.centers)
        errors.append(_l2_1d(phi_num, phi_ref, mesh.volumes))

    errors = np.asarray(errors)
    orders = np.log2(errors[:-1] / errors[1:])
    print(f"cyl_aniso errors: {errors}")
    print(f"cyl_aniso orders: {orders}")
    assert np.all(orders[-2:] > 1.9), (
        f"Expected O(h^2) in the last 2 orders; got orders={orders}, "
        f"errors={errors}"
    )


@pytest.mark.slow
def test_sn_spherical_angular_convergence_at_fixed_mesh():
    r"""Gate 3.3 — angular convergence on spherical at fixed nx.

    Asserts monotone decrease of the L2 error with increasing
    n_ordinates, saturating to the spatial floor.
    """
    nx_fixed = 40
    n_ordinates_list = [4, 8, 16]
    errors = []
    for n_ord in n_ordinates_list:
        case = build_spherical_anisotropic_mms_case(n_ordinates=n_ord)
        mesh = case.build_mesh(nx_fixed)
        Q = case.external_source(mesh)
        result = solve_sn_fixed_source(
            case.materials, mesh, case.quadrature, Q,
            max_inner=500, inner_tol=1e-13,
        )
        phi_num = result.scalar_flux.values[0, :]  # (ng=1, nx, ny=1)
        phi_ref = case.phi_exact(mesh.centers)
        errors.append(_l2_1d(phi_num, phi_ref, mesh.volumes))

    print(f"angular convergence errors: {errors}")
    errors = np.asarray(errors)
    assert errors[1] <= errors[0] * 1.1, (
        f"angular convergence regression: e8={errors[1]} > 1.1·e4={errors[0]}"
    )


# ═══════════════════════════════════════════════════════════════════════
# Strict-xfail spatial variant (psi / qext equation labels + magnitude band)
# ═══════════════════════════════════════════════════════════════════════

@pytest.mark.xfail(
    strict=True,
    reason=(
        "ERR-026 PARTIAL — Phase D (commits 9512459..c44fe9b, 2026-05-12) "
        "shipped the canonical Hébert §3.9.4 Carlson coupled-pole "
        "inward μ=−1 sweep as the M-M angular recurrence's half-angle "
        "seed AND flipped the curvilinear inner solver default to "
        "Krylov.  Per-ordinate flat-flux residual on Gate 1.1 sphere "
        "MMS collapsed to machine precision (≤ 1e-15).  The absolute-"
        "magnitude check at nx=80 (`1e-8 < errors[-1] < 1e-3`) fails "
        "because the L2 error is pre-asymptotic.  Issue #195 tracks the "
        "investigation + marker removal."
    ),
)
@pytest.mark.slow
@pytest.mark.verifies(
    "transport-spherical",
    "sn-mms-spherical-aniso-psi",
    "sn-mms-spherical-aniso-qext",
)
def test_sn_spherical_aniso_mms_converges_second_order():
    r"""Spherical SN with anisotropic ansatz must show :math:`\mathcal{O}(h^2)`
    AND land in the absolute-magnitude band ``1e-8 < err[-1] < 1e-3``.

    Activates the :math:`(1-\mu^2) B(r)/r` angular-redistribution term
    via :math:`\psi_n(r) = (A(r) + B(r) \mu_n)/W`.
    """
    case = build_spherical_anisotropic_mms_case()

    n_cells = [10, 20, 40, 80]
    errors = []
    for nc in n_cells:
        mesh = case.build_mesh(nc)
        Q = case.external_source(mesh)
        result = solve_sn_fixed_source(
            case.materials, mesh, case.quadrature, Q,
            max_inner=500, inner_tol=1e-13,
        )
        phi_num = result.scalar_flux.values[:, 0]
        phi_ref = case.phi_exact(mesh.centers)
        errors.append(_l2_1d(phi_num, phi_ref, mesh.volumes))

    errors = np.asarray(errors)
    orders = np.log2(errors[:-1] / errors[1:])

    assert np.all(orders > 1.9), (
        f"Expected O(h^2), got orders={orders}, errors={errors}"
    )
    assert 1e-8 < errors[-1] < 1e-3


@pytest.mark.xfail(
    strict=True,
    reason=(
        "ERR-026 PARTIAL — coupled to the spherical anisotropic case. "
        "Phase D shipped the Carlson coupled-pole seed for both "
        "geometries; per-ordinate flat-flux identity is closed.  The "
        "absolute-magnitude check at nx=80 still fails because the L2 "
        "error is pre-asymptotic, same root cause as the sphere "
        "variant.  Issue #195 tracks the convergence-magnitude "
        "investigation + marker removal."
    ),
)
@pytest.mark.slow
@pytest.mark.verifies(
    "transport-cylindrical",
    "sn-mms-cylindrical-aniso-psi",
    "sn-mms-cylindrical-aniso-qext",
)
def test_sn_cylindrical_aniso_mms_converges_second_order():
    r"""Cylindrical SN with anisotropic ansatz must show :math:`\mathcal{O}(h^2)`
    AND land in the absolute-magnitude band.

    Activates the :math:`\xi_n^2 B(r)/r` cylindrical analog of the
    spherical redistribution term. :math:`\xi^2 \neq 1 - \eta^2` in
    general, so the cylindrical case exercises a structurally distinct
    quadrature evaluation from the sphere.
    """
    case = build_cylindrical_anisotropic_mms_case()

    n_cells = [10, 20, 40, 80]
    errors = []
    for nc in n_cells:
        mesh = case.build_mesh(nc)
        Q = case.external_source(mesh)
        result = solve_sn_fixed_source(
            case.materials, mesh, case.quadrature, Q,
            max_inner=500, inner_tol=1e-13,
        )
        phi_num = result.scalar_flux.values[:, 0]
        phi_ref = case.phi_exact(mesh.centers)
        errors.append(_l2_1d(phi_num, phi_ref, mesh.volumes))

    errors = np.asarray(errors)
    orders = np.log2(errors[:-1] / errors[1:])

    assert np.all(orders > 1.9), (
        f"Expected O(h^2), got orders={orders}, errors={errors}"
    )
    assert 1e-8 < errors[-1] < 1e-3
