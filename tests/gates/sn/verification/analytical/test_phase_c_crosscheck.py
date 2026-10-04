r"""Issue #168 Phase C/D/E — Gate Set 4: the curvilinear discrete-ordinates solve against a semi-analytical reference.

* Gate 4.1 — homogeneous-reflective :math:`k_\infty` recovery by the SN
  eigenvalue solver, against the closed form :math:`k_\infty =
  \nu\Sigma_f/\Sigma_a`.
* Gate 4.2 — the SN eigenvalue (Phase D) and flux shape (Phase E) against the
  trajectory-resolvent Variant α Green's-function solvers, called bare (the
  Billiard facade routes the multi-region variants since P1 step 2b, #190,
  but these rows predate it and keep their direct calls).

The rows, with the claim each makes:

- ``sphere_2g_homogeneous_dd_n20`` and the two ``cyl_1g_homogeneous`` rows:
  on a uniform reflective medium both methods return :math:`k_\infty` exactly
  (the V_α1 / V_α1_cyl identities), so :math:`k` agrees to ``1e-9``. These are
  the edge rows the heterogeneous rows rest on. They are blind to the angular
  closure: a homogeneous medium's angular flux is near-flat, which nulls the
  redistribution the closure feeds (``[M]`` the cylinder rows' flux moves
  1.1e-10 under a deliberate ``tau := 0.7`` mutation, against 8.8e-2 for the
  2-group 3-region cylinder; ``vv-principles`` anti-pattern #3).
- ``sphere_2g_3reg``: a live SN solve at Gauss-Legendre 32, 40 cells,
  against the reference at :math:`(n_r, n_\mu) = (36, 96)`, eigenvalue and
  flux shape, at tolerances derived from both methods' measured ladders.
  Since #405 P2 step 7b.2.3 these are ``compare_uncertified`` comparisons,
  not verification: the reference's family derives no bound (#566).
- ``cyl_2g_3reg``: a LIVE SN solve at folded 16x32, 40 cells, against the
  reference at :math:`(n_r, n_{\mu,\rm axial}, n_\varphi) = (24, 16, 32)`. The
  eigenvalue and shape rows are strict ``xfail`` on the verification verbs'
  refusal of an uncertified reference (#566, #516). A RECORD row pins today's
  k readings.

What changed on 2026-09-26 (ERR-090, and two defects of this file):

1. The reference fitted one cubic spline to the emission density across the
   material interfaces; it now fits one per region. Its eigenvalue moved by
   up to 2 % (cylinder 1.20693 to 1.23104; sphere 1.35808 to 1.38374 at the
   old resolutions), and it converges in :math:`n_r`.
2. The Phase D sphere row compared the reference with a hand-typed
   :math:`k = 1.3578153` that had been stale since e30d8d14 (2026-05-12, the
   Phase F snapshot regeneration: the snapshot holds 1.3816447). It agreed at
   2e-4 only because the defective reference happened to read 1.35808. The SN
   eigenvalue is now read from the snapshot file.
3. The Phase E rows evaluated the reference at uniformly spaced cell centres,
   but the snapshot meshes are equal-volume (sphere) and equal-area (cylinder)
   (the default rule of the retired region-mesh spelling since b5e85c2d,
   ``CellsByCount.uniform_volume`` now). The reference is now
   compared as volume-weighted cell averages over the SN's own cells, read
   from the mesh that produced the SN flux.

The flux-shape metric changed with item 3. Each profile was normalised to its
own maximum per group, which discards the ratio between the groups and
amplifies an error in the maximum's cell. Both profiles are now scaled to unit
total fission production, one gauge for both groups, so the group ratio is
part of the claim; the metric is the largest cell-average difference relative
to the largest reference cell average.

trajectory_resolvent is the **semi-analytical pillar** (``vv-principles``,
the three pillars): chords integrated with scipy/numpy quadrature along
characteristics, sharing no project primitive with the SN sweep above the
trusted-library line. ORPHEUS SN is the production discretisation under test.
"""
from __future__ import annotations

import functools

import numpy as np
import pytest

from orpheus.derivations.common.xs_library import get_xs
from orpheus.derivations.common.eigenvalue import kinf_homogeneous
from orpheus.derivations.continuous.trajectory_resolvent.greens_function import (
    solve_greens_function_sphere_mg,
    solve_greens_function_sphere_mr,
)
from orpheus.derivations.continuous.trajectory_resolvent.greens_function_cylinder import (
    solve_greens_function_cylinder,
)
from orpheus.geometry import BC, CoordSystem, StructuredGeometry
from orpheus.mesh import CellsByCount, Mesher
from orpheus.numerics.observable import Eigenvalue
from orpheus.reference.verification import compare_uncertified, verify_agreement
from tests.gates.sn._test_helpers import mixture_from_transport_data
from tests.gates.derivations._trajectory_resolvent_ladders import (
    CYLINDER_3REG_SN_STEPS,
    SPHERE_3REG_SN_STEPS,
    sn_residual,
    sphere_3reg_reference_ladder_estimate,
    tolerance_for,
)
from tests.gates.sn.verification.analytical._aba_reference import (
    NO_ESTIMATOR,
    aba_reference,
    assert_cylinder_record,
    awaits_cylinder_bound,
    scaled_tolerance,
    shape_observables,
    verify_cylinder_k,
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
    from orpheus.numerics.quadrature import Quadrature

    # Simple 2G test material (no upscatter, no fission spectrum split).
    sigma_t = [0.5, 1.0]
    sig_s = [[0.3, 0.05], [0.0, 0.7]]  # (g_from, g_to)
    nu_sig_f = [0.4, 0.6]
    chi = [1.0, 0.0]
    mat = mixture_from_transport_data(sigma_t, sig_s, nu_sig_f, chi)

    k_analytical = kinf_homogeneous(
        np.asarray(sigma_t),
        np.asarray(sig_s),
        np.asarray(nu_sig_f),
        np.asarray(chi),
    )

    nx = 20
    mesh = Mesher(StructuredGeometry.sphere(
        (0.0, 2.0), (0,), outer=BC("reflective"),
    )).partition(CellsByCount.uniform_width(nx)).mesh
    quad = Quadrature.gauss_legendre(n_ordinates=8)
    result = solve_sn(
        materials={0: mat},
        mesh=mesh,
        quadrature=quad,
        max_outer=200, keff_tol=1e-9, flux_tol=1e-8,
        max_inner=200, inner_tol=1e-10,
    )
    keff_sn = result.outcome.keff
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

_THIS = "tests/gates/sn/verification/analytical/test_phase_c_crosscheck.py"
_REGIONWISE = (
    "tests/gates/derivations/test_trajectory_resolvent_regionwise_source.py"
    "::test_mr_oracle_first_leg_matches_the_line_integral"
)


def _snapshot(snapshot_id: str):
    r"""The frozen SN regression snapshot ``snapshot_id`` (``scalar_flux`` stored ``(ng, nx)``, ``keff``).

    The SN side of a snapshot row is READ from the file that
    ``test_dd_regression`` pins, never typed here: a hand-typed copy of the
    sphere eigenvalue went stale for four months (see the module docstring).
    """
    from tests.gates.sn._test_helpers import SN_TESTS_ROOT
    path = SN_TESTS_ROOT / "regression" / "snapshots" / f"{snapshot_id}.npz"
    if not path.exists():
        pytest.skip(f"snapshot {snapshot_id!r} not present at {path}")
    return np.load(path)


def _run_sphere_2g_homogeneous_closed() -> float:
    """Bare ``solve_greens_function_sphere_mg``, uniform A, R = 2 cm."""
    A2 = get_xs("A", "2g")
    res = solve_greens_function_sphere_mg(
        R=2.0,
        sigma_t=A2["sig_t"],
        sigma_s=A2["sig_s"],
        nu_sigma_f=A2["nu"] * A2["sig_f"],
        chi=A2["chi"],
        alpha=1.0,                                    # closed sphere
        n_r=16, n_mu=16, n_traj_quad=32,              # V_α1 exact at α=1
        max_iter=20, tol=1e-10,
    )
    return float(res.k_eff)


def _run_cyl_1g_homogeneous_closed() -> float:
    """Bare ``solve_greens_function_cylinder``, uniform A, R = 2 cm."""
    A1 = get_xs("A", "1g")
    res = solve_greens_function_cylinder(
        R=2.0,
        sigma_t=float(A1["sig_t"][0]),
        sigma_s=float(A1["sig_s"][0, 0]),
        nu_sigma_f=float(A1["nu"][0] * A1["sig_f"][0]),
        alpha=1.0,                                    # closed cylinder
        n_r=16, n_mu_axial=12, n_phi_az=24, n_traj_quad=32,
        max_iter=20, tol=1e-10,
    )
    return float(res.k_eff)


# (snapshot_id, runner, rtol, rationale): the homogeneous edge rows. On a
# uniform reflective medium k = k_inf = νΣ_f/Σ_a exactly for both methods
# (V_α1 / V_α1_cyl, derive_T00_equals_P_ss_sphere in the SymPy origins), so
# rtol 1e-9 is the two iterations' floor with headroom.
_GATE_4_2_CASES: tuple[tuple[str, object, float, str], ...] = (
    (
        "sphere_2g_homogeneous_dd_n20",
        _run_sphere_2g_homogeneous_closed,
        1e-9,
        "V_α1 algebraic identity — k=k_∞ exact",
    ),
    (
        "cyl_1g_homogeneous_folded_4x8_dd_n20",
        _run_cyl_1g_homogeneous_closed,
        1e-9,
        "V_α1_cyl algebraic identity — k=k_∞ exact",
    ),
    (
        "cyl_1g_homogeneous_folded_2x4_dd_n20",
        _run_cyl_1g_homogeneous_closed,
        1e-9,
        "V_α1_cyl algebraic identity — k=k_∞ exact (folded 2x4 split)",
    ),
)


@pytest.mark.l1
@pytest.mark.slow
@pytest.mark.verifies("sn-curvilinear-trajectory-resolvent-crosscheck")
@pytest.mark.parametrize(
    "snapshot_id, runner, rtol, rationale",
    _GATE_4_2_CASES,
    ids=[case[0] for case in _GATE_4_2_CASES],
)
def test_phase_d_trajectory_resolvent_crosscheck(
    snapshot_id, runner, rtol, rationale,
) -> None:
    r"""Gate 4.2, the edge rows: SN snapshot k against the trajectory resolvent on a homogeneous medium.

    Both sides must return :math:`k_\infty` of the uniform medium; the
    reference does so by the V_α1 / V_α1_cyl identities
    (:mod:`orpheus.derivations.continuous.trajectory_resolvent.origins.specular.greens_function`).
    These rows are 1-group or homogeneous and so blind to every spatial and
    angular operator (``vv-principles`` anti-pattern #3); they are the
    foundation the heterogeneous rows below rest on, not evidence about them.
    """
    expected_keff = float(_snapshot(snapshot_id)["keff"])
    k_ref = runner()
    rel = abs(k_ref - expected_keff) / expected_keff
    print(
        f"{snapshot_id}: k_sn={expected_keff:.10f}  "
        f"k_trajres={k_ref:.10f}  rel={rel:.2e}  target={rtol:.0e}  "
        f"({rationale})"
    )
    assert rel < rtol, (
        f"Gate 4.2 cross-check for {snapshot_id!r} exceeded tolerance: "
        f"k_sn_snapshot={expected_keff:.8f}, "
        f"k_trajectory_resolvent={k_ref:.8f}, rel={rel:.2e}, "
        f"target rtol={rtol:.0e}. Rationale: {rationale}"
    )


# ════════════════════════════════════════════════════════════════════════
# The heterogeneous rows — eigenvalue and flux shape, against an UNCERTIFIED reference (#405 P2 step 7b.2.3)
# ════════════════════════════════════════════════════════════════════════
#
# Since step 7b.2.3 both sides are read through the reading verb: production's
# ``Solution.read`` (cell averages paired with mesh-free weights) and the
# trajectory-resolvent ``ReferenceSolution`` (its emission density's transport
# integral, the user's ruling 1 of 2026-10-03). The reference's family derives
# no bound (#566; the cylinder also #516), so it reads ``Uncertified``: the
# sphere rows are explicit ``compare_uncertified`` comparisons at their
# 2026-09-26 tolerances, a weaker claim than verification and spelled as such,
# and carry no ``verifies`` marker; the cylinder rows are strict xfails on the
# verbs' own refusal (``ReferenceNotValid``), XPASSing only when P4 gives the
# family a certificate.
#
# The shape observable is, per SN cell and group, the cell average of φ gauged
# to unit total fission production: :func:`._aba_reference.shape_observables`,
# a ``Ratio`` of two flux integrals (the 2026-09-26 metric, the largest cell
# difference relative to the largest reference cell average, is exactly the
# largest per-cell difference against the absolute tolerance τ × M).
#
# Relative tolerances become absolute ones through a scale read off PRODUCTION
# and truncated to three figures (``_aba_reference.scaled_tolerance``), which
# only tightens: the verbs compare absolutely, and a cylinder xfail must refuse
# before the reference is read, so the scale cannot come from the reference.


@functools.cache
def _cylinder_3reg_sn_16x32():
    r"""A live SN solve of the snapshot's cylinder problem at folded 16x32: the ``Solution`` (and its mesh).

    The problem is ``cyl_2g_3reg_folded_4x8_dd_n40``'s (same materials, same
    40-cell equal-area mesh) with the angular grid refined from 4x8 to 16x32
    and the inner budget raised so the inner iteration converges (the
    snapshot's own 300 exits best-effort). The 4x8 snapshot stays pinned
    bit-exactly by ``test_dd_regression``, whose τ sensitivity (8.8e-2 in the
    flux under ``tau := 0.7``) makes it the cylinder angular-closure
    catcher; at 4x8 the SN solve's own angular error in this metric is about
    0.13, larger than any tolerance this comparison can hold, so the
    cross-check is posed where SN is angularly resolved.
    """
    from orpheus.numerics.quadrature import Quadrature
    from tests.gates.sn.regression import _generate_snapshots as generator
    config = {
        **generator._cylinder_3region("2g", 40, "folded_4x8"),
        "quadrature": Quadrature.folded_product(n_mu=16, n_phi=32),
        "max_inner": 2000,
    }
    return generator.run_case(config), config["mesh"]


@functools.cache
def _sphere_3reg_sn_gl32():
    r"""A live SN solve of the snapshot's sphere problem at Gauss-Legendre 32: the ``Solution`` (and its mesh).

    ``sphere_2g_3reg_dd_n40``'s problem (same materials, same 40-cell
    equal-volume mesh) with 32 ordinates instead of 8 and an inner budget
    the inner iteration converges within. The 8-ordinate snapshot stays
    pinned bit-exactly by ``test_dd_regression``; its own angular error in
    the shape metric (1.3e-2) would consume most of this comparison's tolerance.
    """
    from orpheus.numerics.quadrature import Quadrature
    from tests.gates.sn.regression import _generate_snapshots as generator
    config = {
        **generator._sphere_3region("2g", 40),
        "quadrature": Quadrature.gauss_legendre(n_ordinates=32),
        "max_inner": 2000,
    }
    return generator.run_case(config), config["mesh"]


def _largest_cell_average(solution, observables) -> float:
    """M, the largest gauged cell average production reads: the scale of the shape rows' tolerance."""
    return max(solution.read(ratio).value for _, _, ratio in observables)


# ── the sphere: an uncertified comparison ────────────────────────────────
#
# Tolerances are COMPUTED from the ladders in
# tests/gates/derivations/_trajectory_resolvent_ladders.py: tolerance_for(SN
# residual, the reference's ladder ESTIMATE). Since the step-5 ruling (no ladder
# certifies) the estimate is the tolerance's documented provenance, not a bound.
# [M] 2026-09-26: estimate 3.2e-4 (k) and 1.4e-3 (shape); SN residual 1.5e-5 (k)
# and 7.3e-3 (shape); tolerances 4e-3 and 2e-2 (relative).
_SPHERE_LADDER_ESTIMATE = sphere_3reg_reference_ladder_estimate()
_SPHERE_TOLERANCE = {
    observable: tolerance_for(sn_residual(SPHERE_3REG_SN_STEPS[observable]), _SPHERE_LADDER_ESTIMATE[observable])
    for observable in ("k", "shape")
}

_SPHERE_SUPPORTS = (
    f"{_REGIONWISE}[sphere]",
    "tests/gates/derivations/test_peierls_greens_function_mr.py::test_mr_sphere_k_converges_in_n_r",
    f"{_THIS}::test_phase_d_trajectory_resolvent_crosscheck[sphere_2g_homogeneous_dd_n20]",
    f"{_THIS}::test_sn_spherical_homogeneous_kinf_recovery_2g",
    "tests/gates/derivations/test_trajectory_resolvent_reference.py::test_r7b2_2_1_the_sphere_reading_against_an_unsplit_fine_angular_rule",
    "tests/gates/sn/test_solution_read.py::test_r7b2_9_2_the_cells_sum_to_the_whole_domain",
)


@pytest.mark.l1
@pytest.mark.slow
@pytest.mark.rests_on(*_SPHERE_SUPPORTS)
def test_sphere_3reg_k_against_trajectory_resolvent() -> None:
    r"""The heterogeneous closed sphere's eigenvalue: SN at GL32 against the UNCERTIFIED reference at (36, 96).

    Fuel A | moderator B | fuel A at 0.5, 1.5, 2.0 cm, 2 groups, reflective
    at r = R. Held to 4e-3 relative (derived above), as the absolute
    4e-3 × truncated(k_SN). ``[M]`` reading 6.5e-5 relative; the one-spline
    reference (ERR-090) read 7.9e-3. Not verification: the reference has no
    bound (#566), so this is ``compare_uncertified``.

    Blind to the SN angular-closure defect class: ``tau := 0.7`` reads 8.7e-4
    here, inside the tolerance; that defect's catcher is
    ``tests/gates/sn/regression/test_dd_regression.py``. The SN-defect
    witness of this row is a wrong boundary law (the reflective face realised
    as vacuum reads 2.6e1).

    Until 2026-09-26 this row compared the reference at (24, 24) with a
    hand-typed k four months stale (1.3578153, against the snapshot's
    1.3816447) and read 2e-4 under a 2e-2 bound: the stale number and the
    one-spline reference had drifted to the same value.
    """
    solution, _ = _sphere_3reg_sn_gl32()
    tolerance = scaled_tolerance(_SPHERE_TOLERANCE["k"], solution.read(Eigenvalue()).value)
    comparison = compare_uncertified(solution, Eigenvalue(), aba_reference(CoordSystem.SPHERICAL), tolerance)
    print(f"sphere 3-region k: sn={comparison.reading.value!r} ref={comparison.reference_reading.value!r} tol={tolerance:.3e}")
    comparison.require()


@pytest.mark.l1
@pytest.mark.slow
@pytest.mark.rests_on(
    *_SPHERE_SUPPORTS,
    "tests/gates/derivations/test_trajectory_resolvent_reference.py"
    "::test_r7b2_2_2_the_emission_density_is_the_solves_fixed_point[spherical]",
)
def test_sphere_3reg_flux_shape_against_trajectory_resolvent() -> None:
    r"""The heterogeneous closed sphere's flux shape: 80 fission-gauged cell averages over the SN cells.

    Activates the spatial and angular redistribution in both groups and the
    group ratio (one gauge for both groups); the material interfaces are
    inside the comparison. Each of the 80 ratios against the uncertified
    reference at the absolute 2e-2 × M, M the largest gauged cell average.

    ``[M]`` 2026-10-03, the RE-BASELINE of the reading (not of the tolerance):
    until step 7b.2.3 the reference was read as its NODAL φ through a
    per-region cubic spline (4.361e-3 against M); since then it is read as
    the transport integral of its emission density, its natural extension
    (the user's ruling 1 of 2026-10-03: one reference, one answer), and the
    row reads 4.325e-3. The two readings differ by 3.4e-4 in this metric.
    Blind to ``tau := 0.7`` (1.0e-2, measured on the nodal reading) like the
    eigenvalue row. ``[M]`` 2026-10-03, re-measured on this reading (the
    step-7b.2.3 battery): red under the vacuum-for-reflective law and under a
    1 % scaling of the solver's moderator emission density; green under the
    one-spline reference (worst 2.38e-3 against 4.98e-3, which its
    eigenvalue row catches) and under a 1 % scaling of the moderator density
    applied after the solve, at reading time only (worst 1.92e-3: a
    reading-only defect of that size is below this row's resolution). Its
    catcher is the reference's fixed-point identity, R7b2.2.2 (the density
    read reproduces the solve's angular flux to the rounding level), which
    this row rests on.
    """
    solution, mesh = _sphere_3reg_sn_gl32()
    reference = aba_reference(CoordSystem.SPHERICAL)
    observables = shape_observables(mesh)
    tolerance = scaled_tolerance(_SPHERE_TOLERANCE["shape"], _largest_cell_average(solution, observables))
    worst = 0.0
    for i, g, ratio in observables:
        comparison = compare_uncertified(solution, ratio, reference, tolerance)
        worst = max(worst, abs(comparison.reading.value - comparison.reference_reading.value))
        comparison.require()
    print(f"sphere 3-region shape: worst |m - v| = {worst:.3e} against {tolerance:.3e}")


# ── the cylinder: not yet certifiable (#516, #566) ───────────────────────
#
# The tolerances the cylinder rows are held to once the reference is
# certified: tolerance_for(SN residual, None), the reference assumed at the
# floor. [M] SN residual at folded 16x32, 40 cells: 3.0e-5 (k), 5.0e-3
# (shape); tolerances 8e-5 and 2e-2. The verbs refuse the uncertified
# reference before reading it: the expected failure.
_CYLINDER_TOLERANCE = {
    observable: tolerance_for(sn_residual(CYLINDER_3REG_SN_STEPS[observable]), None)
    for observable in ("k", "shape")
}

_CYLINDER_SUPPORTS = (
    f"{_REGIONWISE}[cylinder]",
    f"{_THIS}::test_phase_d_trajectory_resolvent_crosscheck[cyl_1g_homogeneous_folded_4x8_dd_n20]",
    "tests/gates/derivations/test_peierls_greens_function_cylinder_mr.py::test_mr_K3_uniform_reduces_to_mg_2g",
    "tests/gates/derivations/test_peierls_greens_function_cylinder_mr_xverif.py::test_mr_single_region_vacuum_matches_wm72",
)

@pytest.mark.l1
@pytest.mark.slow
@pytest.mark.rests_on(*_CYLINDER_SUPPORTS)
@awaits_cylinder_bound
def test_cylinder_3reg_k_against_trajectory_resolvent() -> None:
    r"""The heterogeneous closed cylinder's eigenvalue: live SN at folded 16x32 against the reference at (24, 16, 32).

    Held to 8e-5 (derived above). The reference has no certificate (#566;
    its azimuthal ladder is not monotone, #516), so ``verify_agreement``
    refuses it before reading: the expected failure. The row XPASSes when
    the family's factory returns a ``Valid`` certificate (P4); the RECORD
    row below keeps the reading live meanwhile.
    """
    verify_cylinder_k(_cylinder_3reg_sn_16x32()[0], _CYLINDER_TOLERANCE["k"])


@pytest.mark.l1
@pytest.mark.slow
@pytest.mark.rests_on(*_CYLINDER_SUPPORTS)
@awaits_cylinder_bound
def test_cylinder_3reg_flux_shape_against_trajectory_resolvent() -> None:
    r"""The heterogeneous closed cylinder's flux shape: the 80 gauged cell averages, each verified once certified.

    Re-posed on 2026-09-26 from the 4x8 snapshot onto a live 16x32 solve: at
    4x8 SN's own angular error in this metric is 6.1e-2, which no comparison
    can separate from a defect; the snapshot keeps its job in
    ``test_dd_regression``. Held to 2e-2 × M (derived above); the expected
    failure is the verbs' refusal, raised at the first ratio.
    """
    solution, mesh = _cylinder_3reg_sn_16x32()
    reference = aba_reference(CoordSystem.CYLINDRICAL)
    observables = shape_observables(mesh)
    tolerance = scaled_tolerance(_CYLINDER_TOLERANCE["shape"], _largest_cell_average(solution, observables))
    for _, _, ratio in observables:
        verify_agreement(solution, ratio, reference, tolerance, NO_ESTIMATOR).require()


def _cylinder_readings() -> dict[str, float]:
    """Today's cylinder k readings, keyed as ``_aba_reference.CYLINDER_3REG_RECORD`` keys them (no shape keys: the user's
    cost ruling of 2026-10-03, the 80-ratio set through the reference's extension costs about 38 min)."""
    solution, _ = _cylinder_3reg_sn_16x32()
    k_sn = solution.read(Eigenvalue()).value
    k_ref = aba_reference(CoordSystem.CYLINDRICAL).read(Eigenvalue()).value
    return {"k_ref": k_ref, "phase_c_k_sn": k_sn, "phase_c_k_gap": abs(k_ref - k_sn) / k_sn}


@pytest.mark.l1
@pytest.mark.slow
@pytest.mark.rests_on(*_CYLINDER_SUPPORTS)
def test_cylinder_3reg_crosscheck_record() -> None:
    r"""RECORD: today's cylinder k readings, green until either side moves.

    Not verification. Red under ``tau := 0.7`` (the SN k moves 7.8e-4), the
    vacuum-for-reflective law, the one-spline reference and a 1 % reference
    perturbation (measured 2026-09-26 with the shape keys; the k keys alone
    are re-measured by the step-7b.2.3 battery). The shape keys were dropped
    at step 7b.2.3 (the user's cost ruling).
    """
    readings = _cylinder_readings()
    print(f"cylinder 3-region readings: {readings}")
    assert_cylinder_record(readings)
