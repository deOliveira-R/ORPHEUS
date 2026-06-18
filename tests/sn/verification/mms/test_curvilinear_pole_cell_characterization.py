r"""Characterization gate — curvilinear SN central-cell (r->0) spatial closure.

Issue #233 (WONTFIX, literature-accepted).  This gate pins a KNOWN
limitation of the plain weighted-diamond curvilinear cell-update
:eq:`dd-curvilinear-scalar` in TEST form.  It is NOT an xfail-pending-fix
and NOT a bug-catcher: it asserts what is TRUE today (the global-L2
guarantee) AND protects the regression floor (the pole order does not
DROP below first order), WITHOUT calcifying the limitation (no upper
bound that a future higher-order central-cell scheme would violate).

The headline — two norms, two orders, ONE defect
=================================================

The curvilinear SN scalar flux is **second-order O(h^2) in the
volume-weighted L2 norm** (``sqrt(sum V_i * diff_i^2)`` — the norm the
production gate ``test_sn_*_mms_converges_second_order`` uses) and
**first-order O(h) in the pointwise / L-infinity norm**, simultaneously.
There is no contradiction: the L-infinity error lives at the SINGLE
central cell (r->0), whose volume is V ~ h^3 (sphere) / h^2 (cylinder).
Its sqrt(V) ~ h^1.5 (sphere) weighting dilutes the pole's O(h)
pointwise error to an ~h^2.5 contribution in the L2 sum — subdominant
by construction.  This is precisely why issue #233 is INVISIBLE to the
production L2/keff gates and required an L-infinity / per-cell-rate
probe to surface.

Root cause (see #233 + ``docs/theory/discrete_ordinates.rst``,
:eq:`dd-curvilinear-scalar`): at the pole the inner face area
A_{i-1/2} = A(0) = 0, so the diamond closure psi_bar = 1/2(psi_in +
psi_out) over-predicts the pole OUTER face by exactly +50%
(mesh-independent), and the conservative balance itself is inconsistent
there (A_in=0 degenerates the streaming surface integral while V~h^3).
Hebert 2009 sec.3.9.4 and Stacey sec.9.9 BOTH use this plain-diamond +
Carlson-starting-direction scheme at the central cell with no special
O(h^2) closure, and neither flags reduced order there.  First order at
the single central cell is the accepted, unflagged behaviour of the
standard DD scheme.  The genuine O(h^2)-at-origin resolution is a
higher-order central-cell scheme (LD #6 / cell-update family #158 /
nodal); when one lands, the characterization tests below still PASS
(order 2.0 > 0.8 lower bound) and document the improvement.

The structurally-independent cross-check (why dual reference)
============================================================

The spherical DD discrete unknown IS the cell-VOLUME-average
(4pi/V) integral r^2 phi dr (Hebert 2009 Eq. 3.430 — the unknown is
DEFINED as the shell average, not a point value).  Comparing the solver
against that shell average is therefore the PRINCIPLED reference (the
discrete unknown compared to its own definition).  The production gate
instead compares against the midpoint point-value phi_exact(centers);
that midpoint comparison could, in principle, manufacture or mask an
apparent order.  This gate asserts O(h^2) under BOTH references: the
midpoint AND the Hebert-3.430 shell-volume-average.  Agreement on order
across two structurally-different references proves the L2 order is REAL
and not a comparison artifact (vv-principles: structural independence
applies to ALL test design, not just numerical cross-checks).  The
shell-average reference is built here from ``scipy.integrate.quad`` — a
trusted-library primitive, structurally independent of the SN solver.

The cylinder shares the IDENTICAL defect, MASKED
================================================

The cylinder pole-cell error vs the MIDPOINT reference is accidentally
O(h^2) (the r dr linear volume weight places the volume-centroid at
~2/3 h, while the diamond face ~1/2 A(h) ~ midpoint A(h/2), so the
midpoint comparison is fortuitously second-order).  But vs the
VOLUME-AVERAGE reference it is O(h) — the same diamond/balance
inconsistency as the sphere.  The earlier "cylinder pole is clean
O(h^2)" reading was an artifact of the midpoint comparison.  The
cylinder global L2 is clean O(h^2) under midpoint.

Measured numbers (this session, ``.venv/bin/python -O``, 2026-06-13;
reproduced exactly from ``diag_14`` / ``diag_30`` in job tmp
``negclosure/``).  Sphere, n_ordinates=16, ladder [40, 80, 160, 320]:

    L2 (midpoint)   = [3.734e-3, 9.279e-4, 2.306e-4, 5.740e-5]  orders 2.01 2.01 2.01
    L2 (shell-avg)  = [8.775e-3, 2.196e-3, 5.492e-4, 1.373e-4]  orders 2.00 2.00 2.00
    L-inf (pole)    = [2.123e-2, 1.132e-2, 5.874e-3, 2.999e-3]  orders 0.91 0.95 0.97
    interior (r/R~0.5)                                          orders 2.02 2.00 2.00
    pole fraction of total L-inf error                         1.00 at every mesh

Cylinder, n_mu=4 n_phi=8, ladder [40, 80, 160, 320]:

    pole vs midpoint    orders 1.94 1.97 1.98   (accidentally O(h^2))
    pole vs volume-avg  orders 0.99 0.99 1.00   (the REAL O(h) defect)
    global L2 (midpt)   = [5.389e-4, 1.347e-4, 3.366e-5, 8.414e-6]  orders 2.00 2.00 2.00

The pole order MEASURES ~1.0 (asymptotically); the assertions
LOWER-BOUND it at 0.8 deliberately, so a future LD/nodal fix that lifts
the pole to O(h^2) keeps this gate GREEN (2.0 > 0.8).

The pole-cell defect is catalogued as **ERR-059** (W5, 2026-06-13;
``.claude/skills/vv-principles/error_catalog.md``) — a DOCUMENTED
INHERENT LIMITATION (WONTFIX-for-DD).  All four tests below carry
``@pytest.mark.catches("ERR-059")``.

References
----------
* Issue #233 — the WONTFIX rationale + full cross-refs.
* ERR-059 — the catalogued pole-cell inherent-limitation entry.
* :eq:`dd-curvilinear-scalar` — the curvilinear DD cell-update whose
  pole behaviour this gate characterizes.
* Hebert 2009 Eq. 3.430 — the shell-volume-average definition of the
  curvilinear DD unknown (the principled reference).
* ``.claude/agent-memory/numerics-investigator/curvilinear_tau_clamp_vs_pole_floor.md``
  — the "REFINED VERDICT 2026-06-13b" pole-O(h) decomposition.
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy import integrate

from orpheus.derivations.continuous.mms.sn import (
    build_cylindrical_mms_case,
    build_spherical_mms_case,
)
from orpheus.sn import solve_sn_fixed_source

from tests.sn._test_helpers import volume_weighted_l2

# ── ladders / quadrature ────────────────────────────────────────────────
# [40, 80, 160, 320] is decisive for the rates (the production gate's
# [20..160] does not resolve the pole asymptote cleanly) and runs in
# ~0.3 s (sphere) / ~2.8 s (cylinder) — fast enough to NOT need @slow.
_LADDER = [40, 80, 160, 320]
_N_ORDINATES_SPHERE = 16


# ── shell-volume-average reference (Hebert 2009 Eq. 3.430) ───────────────
# Built from scipy.integrate.quad: a trusted-library primitive,
# structurally independent of the SN solver.  This is the PRINCIPLED
# reference — the curvilinear DD unknown IS the cell shell average, so
# comparing the solver against it compares the unknown to its definition.

def _shell_avg(profile, lo: float, hi: float, weight_power: int) -> float:
    r"""Volume-average of ``profile`` over the shell [lo, hi].

    ``weight_power=2`` -> sphere (r^2 dr); ``weight_power=1`` -> cylinder
    (r dr).  Returns ``(int profile * r^p dr) / (int r^p dr)``.
    """
    num, _ = integrate.quad(lambda r: profile(r) * r**weight_power, lo, hi)
    den, _ = integrate.quad(lambda r: r**weight_power, lo, hi)
    return num / den


def _shell_avg_reference(case, mesh, weight_power: int) -> np.ndarray:
    """Per-cell shell-volume-average of the manufactured exact profile."""
    radius = case.radius
    profile = lambda r: np.sin(np.pi * r / radius)  # noqa: E731 (the MMS ansatz A(r))
    return np.array([
        _shell_avg(profile, mesh.edges[i], mesh.edges[i + 1], weight_power)
        for i in range(mesh.N)
    ])


def _solve_sphere(nx: int) -> np.ndarray:
    case = build_spherical_mms_case(n_ordinates=_N_ORDINATES_SPHERE)
    mesh = case.build_mesh(nx)
    src = case.external_source(mesh)
    sol = solve_sn_fixed_source(
        case.materials, mesh, case.quadrature, src,
        boundary_condition="vacuum", scattering_order=0,
        max_inner=8000, inner_tol=1e-13,
    )
    return np.asarray(sol.scalar_flux.values)[0]


def _solve_cylinder(nx: int) -> np.ndarray:
    case = build_cylindrical_mms_case()
    mesh = case.build_mesh(nx)
    src = case.external_source(mesh)
    sol = solve_sn_fixed_source(
        case.materials, mesh, case.quadrature, src,
        boundary_condition="vacuum", scattering_order=0,
        max_inner=8000, inner_tol=1e-13,
    )
    return np.asarray(sol.scalar_flux.values)[0]


def _orders(errors: np.ndarray) -> np.ndarray:
    """Successive log2 ratios (uniform mesh-halving ladder)."""
    return np.log2(errors[:-1] / errors[1:])


# ═══════════════════════════════════════════════════════════════════════
# THE GUARANTEE — global volume-weighted L2 is O(h^2) (must hold)
# ═══════════════════════════════════════════════════════════════════════

@pytest.mark.l1
@pytest.mark.catches("ERR-059")
@pytest.mark.verifies(
    "dd-curvilinear-scalar",
    "transport-spherical",
    "sn-mms-spherical-psi", "sn-mms-spherical-qext",
)
def test_sphere_global_L2_second_order_dual_reference():
    r"""GUARANTEE: sphere global L2 is O(h^2) under BOTH references (#233).

    Asserts the production volume-weighted L2 (``sqrt(sum V * diff^2)``)
    converges O(h^2) for the sphere against BOTH the midpoint reference
    AND the Hebert-3.430 shell-volume-average reference.  The
    shell-average reference is the principled one (the curvilinear DD
    unknown IS the shell average); agreement on order across two
    structurally-different references proves the L2 order is REAL, not a
    midpoint artifact.  This is the FALSIFIABLE guarantee of #233: the
    pole O(h) defect does NOT leak into the global L2.

    Measured (n_ordinates=16, [40,80,160,320]):
      L2 midpoint  = [3.734e-3, 9.279e-4, 2.306e-4, 5.740e-5]  orders 2.01 2.01 2.01
      L2 shell-avg = [8.775e-3, 2.196e-3, 5.492e-4, 1.373e-4]  orders 2.00 2.00 2.00
    """
    err_mid = []
    err_va = []
    for nx in _LADDER:
        case = build_spherical_mms_case(n_ordinates=_N_ORDINATES_SPHERE)
        mesh = case.build_mesh(nx)
        src = case.external_source(mesh)
        sol = solve_sn_fixed_source(
            case.materials, mesh, case.quadrature, src,
            boundary_condition="vacuum", scattering_order=0,
            max_inner=8000, inner_tol=1e-13,
        )
        phi = np.asarray(sol.scalar_flux.values)[0]
        ref_mid = case.phi_exact(mesh.centers)
        ref_va = _shell_avg_reference(case, mesh, weight_power=2)
        err_mid.append(volume_weighted_l2(phi, ref_mid, mesh.volumes))
        err_va.append(volume_weighted_l2(phi, ref_va, mesh.volumes))

    err_mid = np.asarray(err_mid)
    err_va = np.asarray(err_va)
    ord_mid = _orders(err_mid)
    ord_va = _orders(err_va)

    assert np.all(ord_mid > 1.9), (
        f"Sphere L2 (midpoint) not O(h^2): orders={ord_mid}, errors={err_mid}"
    )
    assert np.all(ord_va > 1.9), (
        "Sphere L2 (shell-volume-average, Hebert 3.430) not O(h^2): "
        f"orders={ord_va}, errors={err_va} — the L2 order is a MIDPOINT "
        "ARTIFACT if this fails while midpoint passes (#233 falsified)."
    )


@pytest.mark.l1
@pytest.mark.catches("ERR-059")
@pytest.mark.verifies(
    "dd-curvilinear-scalar",
    "transport-cylindrical",
    "sn-mms-cylindrical-psi", "sn-mms-cylindrical-qext",
)
def test_cylinder_global_L2_second_order():
    r"""GUARANTEE: cylinder global L2 is O(h^2) (#233).

    The cylinder shares the identical pole defect as the sphere but its
    global volume-weighted L2 (midpoint reference) stays clean O(h^2) —
    the pole's single-cell O(h) error is diluted by the r dr volume
    weight just as on the sphere.

    Measured (n_mu=4, n_phi=8, [40,80,160,320]):
      L2 midpoint = [5.389e-4, 1.347e-4, 3.366e-5, 8.414e-6]  orders 2.00 2.00 2.00
    """
    errors = []
    for nx in _LADDER:
        case = build_cylindrical_mms_case()
        mesh = case.build_mesh(nx)
        src = case.external_source(mesh)
        sol = solve_sn_fixed_source(
            case.materials, mesh, case.quadrature, src,
            boundary_condition="vacuum", scattering_order=0,
            max_inner=8000, inner_tol=1e-13,
        )
        phi = np.asarray(sol.scalar_flux.values)[0]
        ref_mid = case.phi_exact(mesh.centers)
        errors.append(volume_weighted_l2(phi, ref_mid, mesh.volumes))

    errors = np.asarray(errors)
    orders = _orders(errors)
    assert np.all(orders > 1.9), (
        f"Cylinder global L2 not O(h^2): orders={orders}, errors={errors}"
    )


# ═══════════════════════════════════════════════════════════════════════
# THE CHARACTERIZATION — pole is first-order + L-inf-dominant (#233)
#   bounded BELOW (regression floor), NO upper bound (a future fix MUST
#   still PASS): if LD/nodal lifts the pole to O(h^2), order 2.0 > 0.8.
#   No verifies(...): these pin a LIMITATION, not a correctness claim.
# ═══════════════════════════════════════════════════════════════════════

@pytest.mark.l1
@pytest.mark.catches("ERR-059")
def test_sphere_pole_cell_first_order_and_Linf_dominant():
    r"""CHARACTERIZATION: sphere pole cell is first-order + L-inf-dominant (#233).

    Three pinned facts, all on the isotropic (angular-redistribution-free)
    sphere MMS so the radial DD closure is isolated:

    1. The central (r->0) cell is the dominant L-infinity contributor —
       it IS the maximum pointwise error at every mesh (fraction = 1.0).
    2. The pole error converges at FIRST order — measured ~0.91..0.97
       (approaching 1.0).  Asserted as a LOWER bound (> 0.8): "at least
       first order, does not regress".  NO upper bound is asserted, so a
       future higher-order central-cell scheme (LD #6 / nodal #158/#233)
       that lifts the pole to O(h^2) keeps this test GREEN (2.0 > 0.8).
       This is the regression floor without calcifying the limitation
       (vv-principles anti-pattern #5/#17 — do not assert a claim that
       blocks a legitimate improvement).
    3. The interior is clean O(h^2) — the pole is the SINGLE bad cell;
       everything off the pole is second-order (asserted > 1.8 to leave
       coarse-mesh headroom; measured ~1.84..1.96 on the max-over-interior
       norm, 2.00..2.02 at a fixed fractional position).

    Measured (n_ordinates=16, [40,80,160,320]):
      pole L-inf   = [2.123e-2, 1.132e-2, 5.874e-3, 2.999e-3]  orders 0.91 0.95 0.97
      pole fraction of total L-inf error: 1.00 at every mesh
      interior (max over r/R > 0.1)                           orders 1.84 1.92 1.96

    Catalogued as ERR-059 (DOCUMENTED INHERENT LIMITATION); this test
    carries ``@pytest.mark.catches("ERR-059")``.
    """
    pole_err = []
    interior_err = []
    pole_fraction = []
    for nx in _LADDER:
        case = build_spherical_mms_case(n_ordinates=_N_ORDINATES_SPHERE)
        mesh = case.build_mesh(nx)
        src = case.external_source(mesh)
        sol = solve_sn_fixed_source(
            case.materials, mesh, case.quadrature, src,
            boundary_condition="vacuum", scattering_order=0,
            max_inner=8000, inner_tol=1e-13,
        )
        phi = np.asarray(sol.scalar_flux.values)[0]
        abs_err = np.abs(phi - case.phi_exact(mesh.centers))
        frac = mesh.centers / case.radius
        pole_err.append(abs_err[0])
        interior_err.append(abs_err[frac > 0.1].max())
        pole_fraction.append(abs_err[0] / abs_err.max())

    pole_err = np.asarray(pole_err)
    interior_err = np.asarray(interior_err)

    # Fact 1: the pole cell IS the L-infinity-dominant cell.
    assert np.all(np.asarray(pole_fraction) > 0.99), (
        "Sphere pole cell is NOT the dominant L-inf contributor: "
        f"fractions={pole_fraction} — #233 mis-attributes the L-inf floor "
        "if the dominant cell moved off the pole."
    )

    # Fact 2: pole converges at LEAST first order (lower bound only).
    pole_orders = _orders(pole_err)
    assert np.all(pole_orders > 0.8), (
        f"Sphere pole REGRESSED below first order: orders={pole_orders}, "
        f"errors={pole_err}.  Measured floor is ~1.0; this asserts only "
        "the > 0.8 regression floor (a future O(h^2) fix must still PASS)."
    )

    # Fact 3: the interior is clean O(h^2) — the pole is the single bad cell.
    interior_orders = _orders(interior_err)
    assert np.all(interior_orders > 1.8), (
        f"Sphere interior is NOT clean O(h^2): orders={interior_orders}, "
        f"errors={interior_err} — the pole should be the ONLY first-order cell."
    )


@pytest.mark.l1
@pytest.mark.catches("ERR-059")
def test_cylinder_pole_first_order_vs_volume_average_masked_by_midpoint():
    r"""CHARACTERIZATION: cylinder shares the identical pole defect, MASKED (#233).

    The cylinder pole-cell error vs the MIDPOINT reference is accidentally
    O(h^2) (the r dr linear weight places the volume-centroid such that the
    diamond face ~ midpoint, fortuitously second-order).  But vs the
    VOLUME-AVERAGE reference (the principled Hebert-3.430 unknown) it is
    O(h) — the SAME diamond/balance pole inconsistency as the sphere.
    This test pins both halves so the "cylinder pole is clean O(h^2)"
    misreading cannot creep back in.

    The volume-average order is asserted as a LOWER bound (> 0.8), same
    regression-floor-not-calcification discipline as the sphere pole: a
    future higher-order central-cell scheme that lifts the cylinder pole
    to O(h^2) keeps this GREEN.

    Measured (n_mu=4, n_phi=8, [40,80,160,320]):
      pole vs midpoint     orders 1.94 1.97 1.98  (accidentally O(h^2))
      pole vs volume-avg   orders 0.99 0.99 1.00  (the REAL O(h) defect)
    """
    pole_mid = []
    pole_va = []
    for nx in _LADDER:
        case = build_cylindrical_mms_case()
        mesh = case.build_mesh(nx)
        src = case.external_source(mesh)
        sol = solve_sn_fixed_source(
            case.materials, mesh, case.quadrature, src,
            boundary_condition="vacuum", scattering_order=0,
            max_inner=8000, inner_tol=1e-13,
        )
        phi = np.asarray(sol.scalar_flux.values)[0]
        ref_mid = case.phi_exact(mesh.centers)
        ref_va = _shell_avg_reference(case, mesh, weight_power=1)
        pole_mid.append(abs(phi[0] - ref_mid[0]))
        pole_va.append(abs(phi[0] - ref_va[0]))

    pole_mid = np.asarray(pole_mid)
    pole_va = np.asarray(pole_va)
    ord_mid = _orders(pole_mid)
    ord_va = _orders(pole_va)

    # The midpoint comparison is accidentally second-order — pin it so the
    # masking mechanism is documented and cannot silently change.
    assert np.all(ord_mid > 1.8), (
        f"Cylinder pole vs MIDPOINT lost its accidental O(h^2): "
        f"orders={ord_mid}, errors={pole_mid}."
    )

    # The principled (volume-average) comparison exposes the REAL O(h)
    # defect — lower-bounded so a future O(h^2) fix still PASSES.
    assert np.all(ord_va > 0.8), (
        f"Cylinder pole vs VOLUME-AVERAGE regressed below first order: "
        f"orders={ord_va}, errors={pole_va} (the #233 defect floor)."
    )
    # The masking ITSELF (midpoint accidentally O(h²) WHILE volume-average is
    # genuinely O(h)) is the headline.  It is DOCUMENTED, not hard-asserted:
    # measured 2026-06-13 ord_mid≈[1.94,1.97,1.98] vs ord_va≈[0.99,0.99,1.00].
    # We deliberately do NOT gate on an upper bound for ``ord_va`` (e.g.
    # ``ord_va < 1.5``): that would CALCIFY the limitation (vv anti-pattern
    # #5/#17) — a future higher-order central-cell scheme (#6/#158/#233) that
    # lifts the real pole order to O(h²) must keep this gate GREEN.  The two
    # lower-bound assertions above (midpoint accidental O(h²); volume-average
    # ≥ O(h)) are both fix-survivable AND catch a regression; the current
    # first-order value lives in the printed line + the docstring.
    print(f"cyl pole masking: ord_mid={ord_mid} ord_va={ord_va} "
          f"(midpoint accidentally O(h²); volume-average genuinely O(h))")
