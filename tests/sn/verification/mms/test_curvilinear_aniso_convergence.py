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
the duplicated FILE, not the distinct coverage. All stay ``xfail``,
but the reason changed on 2026-06-12: the ERR-058 closure-seed fix
(#195) closed the curvilinear wrong-fixed-point family (the isotropic
companions now pass plain and their markers are GONE), and the aniso
cases dropped ~50× in error.  What remains here is the
fixed-quadrature ANGULAR floor of the per-ordinate-imposed aniso
ansatz — a test-design limitation tracked at Issue #229 (the M-M
half-angle thread values are interpolated, not imposed; the floor
scales down with quadrature order).

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


def _aniso_sphere_l2_ladder(n_cells, *, n_ordinates: int) -> np.ndarray:
    """Run the spherical anisotropic MMS at each ``n_cells``; return the
    volume-weighted L2 error ladder."""
    case = build_spherical_anisotropic_mms_case(n_ordinates=n_ordinates)
    errors = []
    for nc in n_cells:
        mesh = case.build_mesh(nc)
        Q = case.external_source(mesh)
        result = solve_sn_fixed_source(
            case.materials, mesh, case.quadrature, Q,
            max_inner=500, inner_tol=1e-13,
        )
        phi_num = result.scalar_flux.values[0, :]
        phi_ref = case.phi_exact(mesh.centers)
        errors.append(_l2_1d(phi_num, phi_ref, mesh.volumes))
    return np.asarray(errors)


# ═══════════════════════════════════════════════════════════════════════
# Phase C Gate Set 3 — spatial + angular convergence (catches ERR-026)
# ═══════════════════════════════════════════════════════════════════════

@pytest.mark.slow
@pytest.mark.verifies("sn-mms-spherical-aniso-spatial-convergence")
@pytest.mark.catches("ERR-026")
@pytest.mark.xfail(
    strict=False,
    reason=(
        "ERR-058 (#195, 2026-06-12) CLOSED the curvilinear "
        "wrong-fixed-point family (closure-seed fix); this case's error "
        "dropped ~50x and the coarse-segment orders are now ~1.98.  The "
        "remaining failure is the fixed-quadrature ANGULAR floor of the "
        "per-ordinate-imposed aniso ansatz (the M-M half-angle thread "
        "values are interpolated, not imposed): the fine-segment rate "
        "degrades as the spatial error meets the floor (sphere S16 "
        "~7e-4; cylinder n_mu=4 ~1.9e-2; floor scales with quadrature). "
        "Test-design retune tracked at Issue #229."
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
        "ERR-058 (#195, 2026-06-12) CLOSED the curvilinear "
        "wrong-fixed-point family (closure-seed fix); this case's error "
        "dropped ~50x and the coarse-segment orders are now ~1.98.  The "
        "remaining failure is the fixed-quadrature ANGULAR floor of the "
        "per-ordinate-imposed aniso ansatz (the M-M half-angle thread "
        "values are interpolated, not imposed): the fine-segment rate "
        "degrades as the spatial error meets the floor (sphere S16 "
        "~7e-4; cylinder n_mu=4 ~1.9e-2; floor scales with quadrature). "
        "Test-design retune tracked at Issue #229."
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
        "ERR-058 (#195, 2026-06-12) CLOSED the curvilinear "
        "wrong-fixed-point family (closure-seed fix: coupled-pole "
        "spatial seed + angular-edge-extrapolation half-angle seed); "
        "this case's error dropped ~50x. The REMAINING failure is a "
        "TEST-DESIGN limitation, not a solver bug: the "
        "per-ordinate-imposed anisotropic ansatz leaves the M-M "
        "half-angle thread values INTERPOLATED (not imposed), so at "
        "fixed quadrature the solution converges spatially to an "
        "angular-discretization floor (sphere S16 ~7e-4, S32 ~2.9e-4; "
        "cylinder n_mu=4 ~1.9e-2, n_mu=8 ~7.4e-3 - floor scales with "
        "quadrature). The pure-spatial rate+band assertions cannot both "
        "hold in the tested window at the shipped quadrature. Issue "
        "#229 tracks the quadrature-aware retune + marker removal."
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
        phi_num = result.scalar_flux.values[0, :]  # (ng=1, nx)
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
        "ERR-058 (#195, 2026-06-12) CLOSED the curvilinear "
        "wrong-fixed-point family (closure-seed fix: coupled-pole "
        "spatial seed + angular-edge-extrapolation half-angle seed); "
        "this case's error dropped ~50x. The REMAINING failure is a "
        "TEST-DESIGN limitation, not a solver bug: the "
        "per-ordinate-imposed anisotropic ansatz leaves the M-M "
        "half-angle thread values INTERPOLATED (not imposed), so at "
        "fixed quadrature the solution converges spatially to an "
        "angular-discretization floor (sphere S16 ~7e-4, S32 ~2.9e-4; "
        "cylinder n_mu=4 ~1.9e-2, n_mu=8 ~7.4e-3 - floor scales with "
        "quadrature). The pure-spatial rate+band assertions cannot both "
        "hold in the tested window at the shipped quadrature. Issue "
        "#229 tracks the quadrature-aware retune + marker removal."
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
        phi_num = result.scalar_flux.values[0, :]  # (ng=1, nx)
        phi_ref = case.phi_exact(mesh.centers)
        errors.append(_l2_1d(phi_num, phi_ref, mesh.volumes))

    errors = np.asarray(errors)
    orders = np.log2(errors[:-1] / errors[1:])

    assert np.all(orders > 1.9), (
        f"Expected O(h^2), got orders={orders}, errors={errors}"
    )
    assert 1e-8 < errors[-1] < 1e-3


# ═══════════════════════════════════════════════════════════════════════
# W1 — spherical τ-clamp removal: the anisotropic improvement gate
# ═══════════════════════════════════════════════════════════════════════
#
# W1 unclamps the spherical M-M weight (``tau_mm = tau_raw``).  This
# CHANGES the anisotropic solution.  CRITICAL MEASURED FINDING
# (2026-06-13, test-architect): W1 does NOT give "clean O(h²), no floor"
# on this aniso sphere MMS.  The floor here is the #229 INTERPOLATED-
# half-angle-thread floor (the M-M thread values are interpolated, not
# imposed by the ansatz), which scales with QUADRATURE ORDER and is
# INDEPENDENT of the clamp.  Side-by-side at matched quadrature:
#
#   S16 nx[10,20,40,80,160]:
#     CLAMPED  err [5.94e-2,1.51e-2,3.82e-3,1.16e-3,7.29e-4] floor 7.3e-4
#     UNCLAMP  err [5.92e-2,1.48e-2,3.71e-3,1.40e-3,1.22e-3] floor 1.2e-3
#       → unclamp CLEANS the coarse orders (1.995/1.999 vs 1.979/1.978)
#         but RAISES the fine-mesh floor (1.2e-3 vs 7.3e-4).
#   S32 nx[10,20,40,80]:
#     UNCLAMP  orders [1.985, 2.015, 2.000]  (clean O(h²) full ladder)
#
# So W1's anisotropic "improvement" is a SLIGHTLY cleaner coarse-mesh
# rate, NOT floor removal.  The honest gate (vv anti-pattern #5/#17 — do
# not assert a claim that cannot hold; pin what IS true):
#
#   * S32 has a clean O(h²) window nx∈{10,20,40,80} post-W1 (the floor
#     drops below 1e-3 at this quadrature) — assert min(orders) > 1.9.
#   * S16 holds O(h²) only on the COARSE segment nx∈{10,20,40} post-W1
#     (the floor at nx≥80 degrades it) — assert min(orders[:2]) > 1.9
#     AND the coarse orders clear 1.99 (clamped reaches only ~1.978 —
#     a true W1 discriminator).
#   * The floor is a VERIFIED quantity (it SCALES with quadrature): the
#     S32 floor (nx=160) is below the S16 floor — assert it scales.
#
# These W1 gates do NOT flip the legacy strict-xfail "clean O(h²) on the
# WHOLE S16 ladder" claims above (which W1 does NOT satisfy — the floor is
# #229-interpolation, not the clamp).  The legacy xfails stay until #229
# (quadrature-aware retune) lands.  See the test-architect memo
# curvilinear_aniso_229_9_verification.md for the #229 framing.


@pytest.mark.slow
@pytest.mark.verifies("sn-mms-spherical-aniso-spatial-convergence")
@pytest.mark.catches("ERR-026")
def test_w1_aniso_sphere_S32_clean_o_h2_full_ladder():
    r"""W1 [aniso improvement] Post-unclamp, the spherical anisotropic MMS
    shows clean :math:`O(h^2)` in volume-weighted L2 across the FULL
    ladder ``nx∈{10,20,40,80}`` at S32.

    S32 is required: at S32 the #229 interpolation floor drops below 1e-3,
    so the O(h²) window covers the full nx∈{10,20,40,80}.  (At the default
    S16 the floor at nx=80 degrades the fine order — see the coarse-only
    gate below.)  Volume-weighted L2, NOT L∞ — the pole-cell O(h) defect
    (W2's domain) limits L∞ to ~O(h) but is suppressed in L2 by V∝r²
    (verified 2026-06-13: S32 L∞ orders ~0.82/0.91/0.98 vs L2
    ~1.985/2.015/2.000).

    Measured 2026-06-13 (UNCLAMPED, S32):
    err [5.94e-2, 1.50e-2, 3.71e-3, 9.29e-4], orders [1.985, 2.015, 2.000].

    PRE-W1 (clamped, S32) this ladder reads orders [1.981, 2.001, 1.987] —
    also clean at S32, so this gate is the standing O(h²) claim that W1
    must NOT break (the W1 DISCRIMINATOR is the S16 coarse-rate cleaning,
    next gate).  Stays GREEN on both clamped and W1 builds.
    """
    n_cells = [10, 20, 40, 80]
    errors = _aniso_sphere_l2_ladder(n_cells, n_ordinates=32)
    orders = np.log2(errors[:-1] / errors[1:])
    print(f"W1 aniso S32 errors={errors} orders={orders}")
    assert np.all(orders > 1.9), (
        f"Expected clean O(h^2) full ladder at S32; got orders={orders}, "
        f"errors={errors}"
    )


@pytest.mark.slow
@pytest.mark.catches("ERR-026")
def test_w1_aniso_sphere_S16_coarse_rate_cleaner_unclamped():
    r"""W1 [aniso discriminator] Post-unclamp, the COARSE-segment order of
    the S16 anisotropic sphere MMS is cleaner (closer to 2.0) than the
    clamped baseline.

    This is the gate that actually DISCRIMINATES W1 from the clamped
    baseline on the anisotropic MMS.  The unclamped weight is exact-on-
    linear-in-μ, so it threads the ``B(r)μ`` angular component more
    accurately, cleaning the coarse-mesh order:

      clamped   S16 orders[:2] = [1.979, 1.978]
      unclamped S16 orders[:2] = [1.995, 1.999]

    Asserts the coarse window nx∈{10,20,40} holds ``min > 1.9`` AND clears
    ``1.99`` (the clamped baseline reaches only ~1.978, so 1.99 FAILS on
    a clamped build — a true W1 discriminator).  The full ladder is NOT
    asserted: at S16 the #229 interpolation floor (~1.2e-3 unclamped)
    degrades the fine order to ~1.4, the standing #229 limitation, NOT a
    W1 regression.

    Measured 2026-06-13 (UNCLAMPED, S16, nx[10,20,40]):
    orders [1.99520, 1.99922].
    """
    n_cells = [10, 20, 40]
    errors = _aniso_sphere_l2_ladder(n_cells, n_ordinates=16)
    orders = np.log2(errors[:-1] / errors[1:])
    print(f"W1 aniso S16 coarse errors={errors} orders={orders}")
    assert np.all(orders > 1.9), (
        f"S16 coarse window lost O(h^2); orders={orders}, errors={errors}"
    )
    assert np.min(orders) > 1.99, (
        f"S16 coarse orders={orders} did not clear 1.99 — the unclamped "
        f"exact-on-linear-μ weight should clean the coarse rate past the "
        f"clamped baseline (~1.978)."
    )


@pytest.mark.slow
def test_w1_aniso_sphere_floor_scales_with_quadrature():
    r"""W1 [verified floor] The anisotropic sphere MMS floor SCALES with
    quadrature order (it is the #229 interpolation floor, NOT a fixed
    clamp artefact).

    At fixed-fine ``nx=160`` the volume-weighted L2 error must DROP when
    the quadrature order doubles S16→S32.  This pins the floor as a
    verified, quadrature-dependent quantity (vv-principles: a floor shown
    to scale is a CLAIM, not an unexplained limitation), and distinguishes
    the #229 interpolation floor (scales) from a hypothetical fixed
    closure-bug floor (would not scale).

    Measured 2026-06-13 (UNCLAMPED, nx=160): S16 floor 1.22e-3, S32 floor
    3.57e-4 → ratio 3.4× (gate uses a 2.0× margin).  Holds clamped too
    (clamped S16 7.29e-4, S32 2.89e-4 → 2.5×) — a floor-character gate,
    GREEN on both builds.
    """
    err_s16 = _aniso_sphere_l2_ladder([160], n_ordinates=16)[0]
    err_s32 = _aniso_sphere_l2_ladder([160], n_ordinates=32)[0]
    print(f"W1 aniso floor nx=160: S16={err_s16:.3e} S32={err_s32:.3e} "
          f"ratio={err_s16 / err_s32:.2f}")
    assert err_s32 < err_s16 / 2.0, (
        f"aniso sphere floor did NOT scale with quadrature: S16={err_s16:.3e}, "
        f"S32={err_s32:.3e} (ratio {err_s16 / err_s32:.2f} < 2.0) — the floor "
        f"is not the #229 interpolation floor; investigate a fixed closure "
        f"artefact"
    )
