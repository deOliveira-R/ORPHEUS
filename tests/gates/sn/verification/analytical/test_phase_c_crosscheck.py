r"""Issue #168 Phase C/D/E — Gate Set 4: the curvilinear discrete-ordinates solve against a semi-analytical reference.

* Gate 4.1 — homogeneous-reflective :math:`k_\infty` recovery by the SN
  eigenvalue solver, against the closed form :math:`k_\infty =
  \nu\Sigma_f/\Sigma_a`.
* Gate 4.2 — the SN eigenvalue (Phase D) and flux shape (Phase E) against a
  semi-analytical reference: the characteristic reference
  (:func:`~orpheus.derivations.continuous.characteristic.characteristic_reference`)
  since P1 step (d) of ``.claude/plans/characteristic_reference_architecture.md``
  (2026-10-09), the trajectory-resolvent Variant α Green's-function solvers
  before it. The row names keep the old reference's name: the
  ``#516`` census (``test_crosscheck_harness.py``) keys on them.

The rows, with the claim each makes:

- ``sphere_2g_homogeneous_dd_n20`` and the two ``cyl_1g_homogeneous`` rows:
  on a uniform reflective medium both methods return :math:`k_\infty` exactly
  (the reference because a flat emission lies in its panel space and a closed
  body returns everything, ``test_characteristic_reference.py``'s closed
  homogeneous row), so :math:`k` agrees to each SN snapshot's own distance
  from :math:`k_\infty`, by ``tolerance_for``. These are
  the edge rows the heterogeneous rows rest on. They are blind to the angular
  closure: a homogeneous medium's angular flux is near-flat, which nulls the
  redistribution the closure feeds (``[M]`` the cylinder rows' flux moves
  1.1e-10 under a deliberate ``tau := 0.7`` mutation, against 8.8e-2 for the
  2-group 3-region cylinder; ``vv-principles`` anti-pattern #3).
- ``sphere_2g_3reg``: a live SN solve at Gauss-Legendre 32, 40 cells,
  against the reference at its working point (p = 5,
  ``_characteristic_ladders``), eigenvalue and flux shape, at tolerances
  derived from both methods' measured ladders. These are
  ``compare_uncertified`` comparisons, not verification: the reference derives
  no bound (#566).
- ``cyl_2g_3reg``: a LIVE SN solve at folded 16x32, 40 cells, against the
  reference at its working point (the door default, p = 3). The eigenvalue
  and shape rows are strict ``xfail`` on the verification verbs' refusal of
  an uncertified reference (#566; the old family's azimuthal ladder was also
  non-monotone, #516). A RECORD row pins today's k readings.

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

The characteristic reference is the **semi-analytical pillar**
(``vv-principles``, the three pillars): the transport equation integrated
along the body's lines, Galerkin over them, solved by a dense pencil, sharing
no project primitive with the SN sweep above the trusted-library line (spec
§9 of ``scratch/characteristic_architecture/p1_verification_spec.md`` walks
the axes). What both sides read from one object is the posed problem: the
cross sections, the geometry and the boundary law, whose response factor
(``SpecularReturn.kernel``) both read through ``law.response_kernel``. A
defect there moves both sides together (``[M]`` 2026-10-10, qa: alpha :=
alpha / 2 left the ERR-094 slab row green), so it is pinned outside this
file, by ``tests/gates/derivations/test_characteristic_walls.py`` (red 3
times under that mutation). ORPHEUS SN is the production discretisation under test.
"""
from __future__ import annotations

import functools

import numpy as np
import pytest

from orpheus.derivations.common.eigenvalue import kinf_homogeneous
from orpheus.geometry import BC, CoordSystem, StructuredGeometry
from orpheus.mesh import CellsByCount, Mesher
from orpheus.numerics.observable import Eigenvalue
from orpheus.reference.verification import compare_uncertified, verify_agreement
from tests.gates.sn._test_helpers import mixture_from_transport_data
from tests.gates.derivations._characteristic_ladders import (
    CYLINDER_3REG_SN_STEPS,
    EDGE_SN_RESIDUAL,
    SPHERE_3REG_SN_STEPS,
    WORKING_DEGREE,
    aba_cylinder_shape_error,
    aba_sphere_shape_error,
    edge_specification,
    reference_error,
    rung,
)
from tests.gates.derivations._ladder_rules import sn_residual, tolerance_for
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
# Gate 4.2 — the reference cross-check (the characteristic reference since P1 step (d))
# ═══════════════════════════════════════════════════════════════════════

_THIS = "tests/gates/sn/verification/analytical/test_phase_c_crosscheck.py"
#: The reference's line integral of a per-region source that jumps at both interfaces, against mpmath (until step
#: (d) the retired resolvent's ``test_trajectory_resolvent_regionwise_source.py::
#: test_mr_oracle_first_leg_matches_the_line_integral``).
_REGIONWISE = (
    "tests/gates/derivations/test_characteristic_transport.py"
    "::test_the_outflow_of_a_per_region_polynomial_is_its_line_integral_attenuated_to_the_exit"
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


def _edge_reference_k(body: str) -> float:
    """The characteristic reference's k on the edge body (``_characteristic_ladders.edge_specification``) at its
    working point."""
    from orpheus.derivations.continuous.characteristic import characteristic_reference

    resolution = rung(WORKING_DEGREE[f"edge_{body}"])
    return float(characteristic_reference(edge_specification(body), resolution).read(Eigenvalue()).value)


def _run_sphere_2g_homogeneous_closed() -> float:
    """The reference on the uniform A sphere, R = 2 cm, two groups, reflective."""
    return _edge_reference_k("sphere_2g")


def _run_cyl_1g_homogeneous_closed() -> float:
    """The reference on the uniform A cylinder, R = 2 cm, one group, reflective."""
    return _edge_reference_k("cylinder_1g")


def _edge_tolerance(snapshot_id: str, body: str) -> float:
    """``tolerance_for`` of the snapshot's own distance from k_inf and the reference's error on the body."""
    return tolerance_for(EDGE_SN_RESIDUAL[snapshot_id], reference_error(f"edge_{body}"))


# (snapshot_id, runner, rtol, rationale): the homogeneous edge rows. On a
# uniform reflective medium k = k_inf = νΣ_f/Σ_a exactly for both methods, so
# each side's error is its distance from the closed form: the snapshot's, frozen
# with it ([M] ``_characteristic_ladders.EDGE_SN_RESIDUAL``), and the reference's
# at its working point. Until step (d) the rtol was a typed 1e-9 ("the two
# iterations' floor with headroom").
_GATE_4_2_CASES: tuple[tuple[str, object, float, str], ...] = (
    (
        "sphere_2g_homogeneous_dd_n20",
        _run_sphere_2g_homogeneous_closed,
        _edge_tolerance("sphere_2g_homogeneous_dd_n20", "sphere_2g"),
        "closed homogeneous body — k=k_∞ exact",
    ),
    (
        "cyl_1g_homogeneous_folded_4x8_dd_n20",
        _run_cyl_1g_homogeneous_closed,
        _edge_tolerance("cyl_1g_homogeneous_folded_4x8_dd_n20", "cylinder_1g"),
        "closed homogeneous body — k=k_∞ exact",
    ),
    (
        "cyl_1g_homogeneous_folded_2x4_dd_n20",
        _run_cyl_1g_homogeneous_closed,
        _edge_tolerance("cyl_1g_homogeneous_folded_2x4_dd_n20", "cylinder_1g"),
        "closed homogeneous body — k=k_∞ exact (folded 2x4 split)",
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
    r"""Gate 4.2, the edge rows: SN snapshot k against the characteristic reference on a homogeneous medium.

    Both sides must return :math:`k_\infty` of the uniform medium; the
    reference does so because a flat emission lies in its panel space and the
    closed body returns everything it emits
    (``test_characteristic_reference.py::test_a_closed_homogeneous_body_reads_the_exact_infinite_mediums_k_and_flux``).
    The tolerance and the reading per row are in the ``_GATE_4_2_CASES``
    comment's source, ``_characteristic_ladders``. These rows are 1-group or homogeneous and so blind to every spatial and
    angular operator (``vv-principles`` anti-pattern #3); they are the
    foundation the heterogeneous rows below rest on, not evidence about them.
    """
    expected_keff = float(_snapshot(snapshot_id)["keff"])
    k_ref = runner()
    rel = abs(k_ref - expected_keff) / expected_keff
    print(
        f"{snapshot_id}: k_sn={expected_keff:.10f}  "
        f"k_ref={k_ref:.16f}  rel={rel:.2e}  target={rtol:.0e}  "
        f"({rationale})"
    )
    assert rel < rtol, (
        f"Gate 4.2 cross-check for {snapshot_id!r} exceeded tolerance: "
        f"k_sn_snapshot={expected_keff:.8f}, "
        f"k_reference={k_ref:.16f}, rel={rel:.2e}, "
        f"target rtol={rtol:.0e}. Rationale: {rationale}"
    )


# ════════════════════════════════════════════════════════════════════════
# The heterogeneous rows — eigenvalue and flux shape, against an UNCERTIFIED reference (#405 P2 step 7b.2.3)
# ════════════════════════════════════════════════════════════════════════
#
# Since step 7b.2.3 both sides are read through the reading verb: production's
# ``Solution.read`` (cell averages paired with mesh-free weights) and the
# reference's ``ReferenceSolution`` (since P1 step (d) the characteristic
# reference, its Galerkin flux paired with the weights). The reference derives
# no bound (#566), so it reads ``Uncertified``: the
# sphere rows are explicit ``compare_uncertified`` comparisons at their
# ladder-derived tolerances, a weaker claim than verification and spelled as such,
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
# tests/gates/derivations/_characteristic_ladders.py: tolerance_for(SN
# residual, the reference's ladder ESTIMATE). Since the step-5 ruling (no ladder
# certifies) the estimate is the tolerance's documented provenance, not a bound.
# [M] 2026-10-09: SN residual 1.5e-5 (k) and 7.3e-3 (shape); tolerances 4e-5
# and 2e-2 (relative), the reference's estimates at p = 5 three orders and more
# below the SN residual. Until step (d) the estimate was the trajectory
# resolvent's, 3.2e-4 (k) and 1.4e-3 (shape), and the k tolerance its 10 b
# floor, 4e-3.
_SPHERE_LADDER_ESTIMATE = {"k": reference_error("aba_sphere"), "shape": aba_sphere_shape_error()}
_SPHERE_TOLERANCE = {
    observable: tolerance_for(sn_residual(SPHERE_3REG_SN_STEPS[observable]), _SPHERE_LADDER_ESTIMATE[observable])
    for observable in ("k", "shape")
}

#: Until step (d): the retired resolvent's ``test_peierls_greens_function_mr.py::test_mr_sphere_k_converges_in_n_r``
#: (now the Rayleigh-Ritz nesting row) and ``test_trajectory_resolvent_reference.py::
#: test_r7b2_2_1_the_sphere_reading_against_an_unsplit_fine_angular_rule`` (now the flux integral of a step weight,
#: the shape observables' shape).
_SPHERE_SUPPORTS = (
    f"{_REGIONWISE}[sphere_solid_b0.3-1]",
    "tests/gates/derivations/test_characteristic_system.py::test_one_group_k_increases_on_nested_spaces[sphere]",
    f"{_THIS}::test_phase_d_trajectory_resolvent_crosscheck[sphere_2g_homogeneous_dd_n20]",
    f"{_THIS}::test_sn_spherical_homogeneous_kinf_recovery_2g",
    "tests/gates/derivations/test_characteristic_reference.py::test_a_flux_integral_is_the_volume_integral_of_the_galerkin_flux[step-symbolic]",
    "tests/gates/sn/test_solution_read.py::test_r7b2_9_2_the_cells_sum_to_the_whole_domain",
)


@pytest.mark.l1
@pytest.mark.slow
@pytest.mark.rests_on(*_SPHERE_SUPPORTS)
def test_sphere_3reg_k_against_trajectory_resolvent() -> None:
    r"""The heterogeneous closed sphere's eigenvalue: SN at GL32 against the UNCERTIFIED characteristic reference at p = 5.

    Fuel A | moderator B | fuel A at 0.5, 1.5, 2.0 cm, 2 groups, reflective
    at r = R. Held to 4e-5 relative (derived above: 2 e governs, the SN
    residual 1.5e-5; the reference's estimate 3.9e-11), as the absolute
    4e-5 × truncated(k_SN). ``[M]`` 2026-10-09 reading 1.8e-5 relative, a
    margin of 2.2. Not verification: the reference has no bound (#566), so
    this is ``compare_uncertified``.

    Until P1 step (d) the reference was the trajectory resolvent at
    (36, 96), its ladder estimate 3.2e-4 set the tolerance at 4e-3 (its 10 b
    floor), and the row read 6.5e-5; the one-spline resolvent (ERR-090) read
    7.9e-3. At that tolerance the row was blind to the SN angular-closure
    defect class (``tau := 0.7`` read 8.7e-4, 2026-09-26). At 4e-5 the same
    reading would be 22 times the tolerance, so the row is PROMOTED to a
    catcher of it ``[R]``, not re-measured at step (d); the class's catcher
    of record stays ``tests/gates/sn/regression/test_dd_regression.py``. The
    SN-defect witness of this row is a wrong boundary law (the reflective face
    realised as vacuum reads 2.6e1) and, at step (d), SN's k scaled by
    1 + 2 × 4e-5 (red, ``p1_step_d/ta_d/battery``).

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
    "tests/gates/derivations/test_characteristic_system.py"
    "::test_the_fundamental_mode_satisfies_its_pencil[closed-sphere]",
)
def test_sphere_3reg_flux_shape_against_trajectory_resolvent() -> None:
    r"""The heterogeneous closed sphere's flux shape: 80 fission-gauged cell averages over the SN cells.

    Activates the spatial and angular redistribution in both groups and the
    group ratio (one gauge for both groups); the material interfaces are
    inside the comparison. Each of the 80 ratios against the uncertified
    reference at the absolute 2e-2 × M, M the largest gauged cell average:
    the SN residual (7.3e-3) governs, the characteristic reference's estimate
    at p = 5 being 9.4e-6 (1.4e-3 for the trajectory resolvent it replaced at
    P1 step (d), which gave the same 2e-2).

    The readings below were taken on the trajectory resolvent; at step (d)
    the row reads the characteristic reference (``[M]`` 2026-10-09:
    ``p1_step_d/ta_d/``, the step-(d) runs), and its SN-side witness is SN's
    every ratio scaled by 1 + 2 × 2e-2 (red, ``p1_step_d/ta_d/battery``).

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


# ── the cylinder: not yet certifiable (#566) ─────────────────────────────
#
# The tolerances the cylinder rows are held to once the reference is
# certified: tolerance_for(SN residual, the characteristic reference's estimate
# at the door default). [M] SN residual at folded 16x32, 40 cells: 3.0e-5 (k),
# 5.0e-3 (shape); the reference's estimates 3.5e-8 (k) and 2.5e-4 (shape,
# _characteristic_ladders, 2026-10-10); tolerances 6e-5 and 2e-2. Until step
# (d) the trajectory resolvent had no estimate on the cylinder (#516), the
# reference was assumed at the floor, and the k tolerance was 8e-5. The verbs
# refuse the uncertified reference before reading it: the expected failure.
_CYLINDER_TOLERANCE = {
    "k": tolerance_for(sn_residual(CYLINDER_3REG_SN_STEPS["k"]), reference_error("aba_cylinder")),
    "shape": tolerance_for(sn_residual(CYLINDER_3REG_SN_STEPS["shape"]), aba_cylinder_shape_error()),
}

#: Until step (d): the retired resolvent's ``test_peierls_greens_function_cylinder_mr.py::
#: test_mr_K3_uniform_reduces_to_mg_2g`` (now: an interface between equal materials is invisible) and
#: ``test_peierls_greens_function_cylinder_mr_xverif.py::test_mr_single_region_vacuum_matches_wm72`` (now: the
#: cylinder's escape and transmission probabilities against their closed forms).
_CYLINDER_SUPPORTS = (
    f"{_REGIONWISE}[cylinder_solid_b0.7_wz0.8-1]",
    f"{_THIS}::test_phase_d_trajectory_resolvent_crosscheck[cyl_1g_homogeneous_folded_4x8_dd_n20]",
    "tests/gates/derivations/test_characteristic_assembly.py::test_an_interface_between_equal_materials_is_invisible[3]",
    "tests/gates/derivations/test_characteristic_assembly.py::test_the_escape_and_transmission_probabilities_are_the_closed_forms[cylinder-tau2.0]",
)

@pytest.mark.l1
@pytest.mark.slow
@pytest.mark.rests_on(*_CYLINDER_SUPPORTS)
@awaits_cylinder_bound
def test_cylinder_3reg_k_against_trajectory_resolvent() -> None:
    r"""The heterogeneous closed cylinder's eigenvalue: live SN at folded 16x32 against the reference at the door default.

    Held to 6e-5 (derived above). The reference has no certificate (#566),
    so ``verify_agreement`` refuses it before reading: the expected failure.
    The row XPASSes when the reference returns a ``Valid`` certificate (P4);
    the RECORD row below keeps the reading live meanwhile (``[M]``
    2026-10-10: the gap it records is 3.5e-5, inside 6e-5).
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
    ``test_dd_regression``. Held to 2e-2 × M (derived above, unchanged by
    step (d): the SN residual governs); the expected failure is the verbs'
    refusal, raised at the first ratio.
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

    Not verification. Since P1 step (d) ``k_ref`` is the characteristic
    reference's at the door default (re-baselined 2026-10-10, the three SN
    keys unmoved). Its SN-side witness at step (d): SN's k scaled by
    1 + 2 × 2e-5 (red, ``p1_step_d/ta_d/battery``). Earlier: red under ``tau := 0.7`` (the SN k moves 7.8e-4), the
    vacuum-for-reflective law, the one-spline reference and a 1 % reference
    perturbation (measured 2026-09-26 with the shape keys; the k keys alone
    are re-measured by the step-7b.2.3 battery). The shape keys were dropped
    at step 7b.2.3 (the user's cost ruling).
    """
    readings = _cylinder_readings()
    print(f"cylinder 3-region readings: {readings}")
    assert_cylinder_record(readings)
