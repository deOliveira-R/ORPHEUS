r"""Cross-method eigenvalue / critical-dimension regression gates.

Tests in this file exercise the
:class:`~tests.cross_method.protocol.SolverAdapter` protocol over
the populated case sets in :mod:`~tests.cross_method.cases`.

Relationship to direct math-heart-class construction (Phase D)
--------------------------------------------------------------

The math-heart classes
(:class:`~orpheus.derivations.continuous.trajectory_resolvent.Billiard`,
:class:`~orpheus.derivations.continuous.fn_method.moment_space.MomentSpace`)
are constructed directly with a :class:`StructuredGeometry` plus
``materials: dict[int, Mixture]``. **The tests in this file
deliberately keep their names + per-method bodies** — pytest
collection IDs are preserved (CI / pytest-xdist contract) and the
per-method adapter classes (``FNSlabAdapter``,
``TrajectoryResolventSphereAdapter``, ...) stay as the
unit-conversion layer. The agreement between the adapter route and
the direct-construction route is exercised by
:mod:`tests.cross_method.test_polymorphism` (foundation-tier
regression net, 5 tests).

Three classes of test:

1. **Truth gates** — each adapter reproduces its case's truth value
   to within the case's per-adapter tolerance. One parametrised
   test per (case × adapter) pair. These are the L1 backing —
   each method's agreement with closed-form truth.

2. **Cross-method agreement gates** — pairs of adapters that both
   support a case agree to within the larger of their two truth
   tolerances (per
   :func:`~tests.cross_method.protocol.agreement_tolerance` —
   tighter is reference contamination). These are the L4 cross-
   implementation gates, each backed by its case's L1 truth.

3. **Schema gates** (foundation) — every case has at least one
   declared tolerance; every adapter named in tolerances exists
   in :data:`ADAPTERS_BY_NAME`; the case sets enumerate disjoint
   case_id values. Catch silent drift in the case populations.

V&V tagging
-----------

* Foundation gates: ``@pytest.mark.foundation``. Software invariants
  on the protocol metadata.
* Truth gates: ``@pytest.mark.l1``. Each method matches a closed-
  form / semi-analytical truth value from primary literature.
* Cross-method agreement gates: also ``@pytest.mark.l1``. The
  conceptual level per :doc:`/skills/vv-principles` §"V&V level
  taxonomy" is L4 (code-to-code agreement) but the codebase's
  existing cross-method gates
  (``test_fn_la13511_slab_xverif.py``,
  ``test_fn_la13511_sphere_xverif.py``) tag these as L1 because
  the agreement is **L1-strength evidence** for either method
  when both methods' L1 truth-backing is established and
  structural independence is genuine. The L1 backing here is the
  per-adapter truth gates in the same file.

Slow tests
----------

The trajectory_resolvent slab vacuum-BC adapter at default
quadrature (``n_x=48, n_mu=128, n_traj_quad=96``) takes ~30 s per
solve. These tests carry ``@pytest.mark.slow``. Use ``pytest -m
"l1 and not slow"`` for the fast subset.
"""
from __future__ import annotations

import pytest

# Suppress the F_N bracket-scan divide-by-zero warnings (intermediate
# `a` values give near-singular matrices the bracket scan correctly
# brackets through; not a numerical pathology).
pytestmark = [
    pytest.mark.filterwarnings(
        "ignore:divide by zero encountered in det:RuntimeWarning"
    ),
    pytest.mark.filterwarnings(
        "ignore:invalid value encountered in det:RuntimeWarning"
    ),
]

from dataclasses import replace

from orpheus.geometry.structured_geometry import (
    Region,
    StructuredGeometry,
)

from .adapters import (
    ADAPTERS_BY_NAME,
    FNReflectedSlabAdapter,
    FNSlabAdapter,
    FNSphereAdapter,
    TrajectoryResolventSlabAdapter,
    TrajectoryResolventSphereAdapter,
    TrajectoryResolventSphereClosedAdapter,
    _extract_1g_xs,
)
from .cases import (
    ALL_CASES,
    BARE_CRITICAL_SLAB_CASES,
    BARE_CRITICAL_SPHERE_CASES,
    CLOSED_SPHERE_KINF_CASES,
    GRANDJEAN_SIEWERT_SLAB_PARAMETRIC,
    REFLECTED_SLAB_CASES,
)
from .protocol import (
    CrossMethodCase,
    agreement_tolerance,
)


def _shadow_with_thickness_mfp(
    case: CrossMethodCase,
    *,
    a_critical_mfp: float | None = None,
    R_critical_mfp: float | None = None,
) -> CrossMethodCase:
    """Return a copy of ``case`` whose ``structured_geometry`` encodes
    a different critical dimension (in mfp).

    Used by cross-method agreement tests to feed one method's
    predicted critical dimension into another method's adapter
    without altering the underlying truth or XS. Exactly one of
    ``a_critical_mfp`` (slab half-thickness) or ``R_critical_mfp``
    (sphere radius) must be provided. The cm value is derived via
    ``critical_dimension_mfp / sigma_t`` (mfp ↔ cm conversion using
    the case's own σ_t).

    The shadow case sets only ``structured_geometry`` (not
    ``materials``), so the registry-backed ``materials`` path is
    preserved — this is the protocol's "Override" path
    (registry_case + inline structured_geometry, materials=None).
    """
    if (a_critical_mfp is None) == (R_critical_mfp is None):
        raise ValueError(
            "Provide exactly one of a_critical_mfp or R_critical_mfp."
        )
    sigma_t, _, _ = _extract_1g_xs(case)
    if case.registry_case is not None and hasattr(
        case.registry_case, "to_geometry"
    ):
        base_geom = case.registry_case.to_geometry()
    elif case.structured_geometry is not None:
        base_geom = case.structured_geometry
    else:
        raise ValueError(
            f"Case {case.case_id!r} has no geometry to shadow."
        )

    cd_mfp = a_critical_mfp if a_critical_mfp is not None else R_critical_mfp
    cd_cm = float(cd_mfp) / sigma_t

    if base_geom.geometry == "SLB":
        # Slab: published critical dimension is the half-thickness;
        # the structured-geometry extent is the FULL slab width.
        full_width_cm = 2.0 * cd_cm
        new_geom = StructuredGeometry(
            geometry="SLB",
            regions=(Region(mat_id=0, outer_thickness_cm=full_width_cm),),
            bcs=base_geom.bcs,
        )
    else:
        # SPH / CYL: published critical dimension IS the radius.
        new_geom = StructuredGeometry(
            geometry=base_geom.geometry,
            regions=(Region(mat_id=0, outer_thickness_cm=cd_cm),),
            bcs=base_geom.bcs,
        )
    return replace(case, structured_geometry=new_geom)


# ═══════════════════════════════════════════════════════════════════
# Foundation — schema gates
# ═══════════════════════════════════════════════════════════════════


@pytest.mark.foundation
def test_case_inventory_disjoint_ids():
    """Every CrossMethodCase ``case_id`` is unique."""
    ids = [c.case_id for c in ALL_CASES]
    assert len(ids) == len(set(ids)), (
        f"Duplicate case_ids in ALL_CASES: "
        f"{[i for i in ids if ids.count(i) > 1]}"
    )


@pytest.mark.foundation
@pytest.mark.parametrize("case", ALL_CASES, ids=lambda c: c.case_id)
def test_case_has_at_least_one_tolerance(case: CrossMethodCase):
    """Every case opts in at least one adapter via tolerances."""
    assert len(case.tolerances) >= 1, (
        f"Case {case.case_id!r} has no adapter tolerances declared. "
        f"Either add an adapter or remove the case."
    )


@pytest.mark.foundation
@pytest.mark.parametrize("case", ALL_CASES, ids=lambda c: c.case_id)
def test_case_tolerance_adapters_exist(case: CrossMethodCase):
    """Every adapter named in case.tolerances is registered."""
    for adapter_name in case.tolerances:
        assert adapter_name in ADAPTERS_BY_NAME, (
            f"Case {case.case_id!r}: tolerance for {adapter_name!r} "
            f"but adapter not in ADAPTERS_BY_NAME = "
            f"{sorted(ADAPTERS_BY_NAME)}"
        )


@pytest.mark.foundation
@pytest.mark.parametrize("case", ALL_CASES, ids=lambda c: c.case_id)
def test_case_pillar_is_not_ancillary(case: CrossMethodCase):
    """Truth pillars must be closed-form / MMS / semi-analytical.

    "Ancillary" is reserved for cross-implementation references;
    those should not back a truth value. Per
    :doc:`/skills/vv-principles`.
    """
    assert case.pillar != "ancillary", (
        f"Case {case.case_id!r} pillar is 'ancillary' — that is the "
        f"L4 reference status, not a verification pillar. Backed-by "
        f"truth values must be closed-form / MMS / semi-analytical."
    )


# ═══════════════════════════════════════════════════════════════════
# Truth gates — bare-critical slab (fn_method)
# ═══════════════════════════════════════════════════════════════════


@pytest.mark.l1
@pytest.mark.parametrize(
    "case", BARE_CRITICAL_SLAB_CASES, ids=lambda c: c.case_id,
)
def test_fn_slab_matches_truth(case: CrossMethodCase):
    """F_N slab reproduces the case's truth ``a_critical_mfp``.

    Backed by Sood LA-13511 / Grandjean-Siewert / KLL via Sood
    transcription.
    """
    adapter = FNSlabAdapter()
    if adapter.name not in case.tolerances:
        pytest.skip(
            f"fn_slab not opted in for {case.case_id!r}"
        )
    res = adapter.solve(case)
    tol = case.tolerance_for(adapter)
    assert abs(res.value - case.truth_value) < tol, (
        f"{case.case_id}: fn_slab {res.value:.10f} vs truth "
        f"{case.truth_value} (source: {case.truth_source}) "
        f"diff={abs(res.value - case.truth_value):.3e} > tol={tol:.1e}"
    )


@pytest.mark.l1
@pytest.mark.parametrize(
    "case", GRANDJEAN_SIEWERT_SLAB_PARAMETRIC, ids=lambda c: c.case_id,
)
def test_fn_slab_grandjean_siewert_table_xi(case: CrossMethodCase):
    """F_N slab reproduces Grandjean-Siewert Table XI for the
    parametric c-sweep with unit XS.

    These cases extend the slab c-sweep beyond the Sood family
    (c=1.10, 1.70, 1.90 not in Sood). They are fn_method-only — no
    trajectory_resolvent counterpart since the unit-XS path doesn't
    flow through the registry.
    """
    # GS Table XI cases use the c parameter directly, not registry XS.
    # The FN slab solver takes c as input; we route via a special path.
    from orpheus.derivations.continuous.fn_method.slab import (
        solve_fn_slab_bare_critical,
    )
    # Extract c from case_id (encoded as "GS-Table-XI-slab-c1.10").
    c = float(case.case_id.split("-c")[-1])
    res = solve_fn_slab_bare_critical(c=c, n_modes=10)
    tol = case.tolerance_for("fn_slab")
    diff = abs(res.a_critical_mfp - case.truth_value)
    assert diff < tol, (
        f"{case.case_id}: fn_slab a_c={res.a_critical_mfp:.10f} vs "
        f"GS Table XI {case.truth_value:.10f}, diff={diff:.3e} > "
        f"tol={tol:.1e}"
    )


# ═══════════════════════════════════════════════════════════════════
# Truth gates — bare-critical sphere (fn_method)
# ═══════════════════════════════════════════════════════════════════


@pytest.mark.l1
@pytest.mark.parametrize(
    "case", BARE_CRITICAL_SPHERE_CASES, ids=lambda c: c.case_id,
)
def test_fn_sphere_matches_truth(case: CrossMethodCase):
    """F_N sphere reproduces the case's truth ``R_critical_mfp``.

    Backed by Sood LA-13511 / KLL via Sood transcription. F_N
    sphere at N=10 reaches ~5e-8 absolute on R_c.
    """
    adapter = FNSphereAdapter()
    if adapter.name not in case.tolerances:
        pytest.skip(
            f"fn_sphere not opted in for {case.case_id!r}"
        )
    res = adapter.solve(case)
    tol = case.tolerance_for(adapter)
    assert abs(res.value - case.truth_value) < tol, (
        f"{case.case_id}: fn_sphere {res.value:.10f} vs truth "
        f"{case.truth_value} (source: {case.truth_source}) "
        f"diff={abs(res.value - case.truth_value):.3e} > tol={tol:.1e}"
    )


# ═══════════════════════════════════════════════════════════════════
# Truth gates — bare-critical slab (trajectory_resolvent)
# ═══════════════════════════════════════════════════════════════════


@pytest.mark.l1
@pytest.mark.slow
@pytest.mark.parametrize(
    "case", BARE_CRITICAL_SLAB_CASES, ids=lambda c: c.case_id,
)
def test_trajectory_resolvent_slab_matches_truth_keff_one(
    case: CrossMethodCase,
):
    """trajectory_resolvent slab vacuum-BC at the case's truth half-
    thickness reproduces ``k_eff = 1.0``.

    The truth comes from F_N (or KLL via Sood). Agreement at
    ``k_eff = 1.0`` is the structural-independence pillar — the
    bouncing-trajectory operator agrees with the Case singular-
    eigenfunction representation that produced the truth.
    """
    adapter = TrajectoryResolventSlabAdapter()
    if adapter.name not in case.tolerances:
        pytest.skip(
            f"trajectory_resolvent_slab not opted in for {case.case_id!r}"
        )
    res = adapter.solve(case)
    tol = case.tolerance_for(adapter)
    assert res.metadata["converged"], (
        f"{case.case_id}: trajectory_resolvent_slab did not converge "
        f"in {res.metadata['iterations']} iter (n_x={res.metadata['n_x']}, "
        f"n_mu={res.metadata['n_mu']}, n_traj_quad="
        f"{res.metadata['n_traj_quad']})"
    )
    assert abs(res.value - 1.0) < tol, (
        f"{case.case_id}: trajectory_resolvent_slab k_eff="
        f"{res.value:.8f} at truth half-thickness "
        f"{case.truth_value} mfp, |k-1|={abs(res.value-1.0):.3e} > "
        f"tol={tol:.1e}. Truth source: {case.truth_source}"
    )


# ═══════════════════════════════════════════════════════════════════
# Truth gates — bare-critical sphere (trajectory_resolvent)
# ═══════════════════════════════════════════════════════════════════


@pytest.mark.l1
@pytest.mark.parametrize(
    "case", BARE_CRITICAL_SPHERE_CASES, ids=lambda c: c.case_id,
)
def test_trajectory_resolvent_sphere_matches_truth_keff_one(
    case: CrossMethodCase,
):
    """trajectory_resolvent sphere vacuum-BC at the case's truth radius
    reproduces ``k_eff = 1.0``.

    Sphere is NOT slow (no μ=0 cusp); fast quadrature suffices.
    """
    adapter = TrajectoryResolventSphereAdapter()
    if adapter.name not in case.tolerances:
        pytest.skip(
            f"trajectory_resolvent_sphere not opted in for "
            f"{case.case_id!r}"
        )
    res = adapter.solve(case)
    tol = case.tolerance_for(adapter)
    assert res.metadata["converged"], (
        f"{case.case_id}: trajectory_resolvent_sphere did not converge"
    )
    assert abs(res.value - 1.0) < tol, (
        f"{case.case_id}: trajectory_resolvent_sphere k_eff="
        f"{res.value:.8f} at truth R_c "
        f"{case.truth_value} mfp, |k-1|={abs(res.value-1.0):.3e} > "
        f"tol={tol:.1e}. Truth source: {case.truth_source}"
    )


# ═══════════════════════════════════════════════════════════════════
# Truth gates — reflected slab (fn_method only; one-sided coverage)
# ═══════════════════════════════════════════════════════════════════


@pytest.mark.l1
@pytest.mark.parametrize(
    "case", REFLECTED_SLAB_CASES, ids=lambda c: c.case_id,
)
def test_fn_reflected_slab_matches_truth(case: CrossMethodCase):
    """F_N reflected slab reproduces the case's truth ``tau_critical_mfp``.

    Backed by Sood LA-13511 Table 10 + NM 1980 Table 2 + Burkart
    1976 'Exact'. **No trajectory_resolvent counterpart** — this is
    one-sided coverage. The trajectory_resolvent slab has an
    asymmetric variant
    (``solve_greens_function_slab_asymmetric``) that COULD be
    extended to host a reflector via partial-α boundary conditions,
    but that work is not in this task's scope.
    """
    adapter = FNReflectedSlabAdapter()
    if adapter.name not in case.tolerances:
        pytest.skip(
            f"fn_reflected_slab not opted in for {case.case_id!r}"
        )
    res = adapter.solve(case)
    tol = case.tolerance_for(adapter)
    assert res.metadata["converged"], (
        f"{case.case_id}: fn_reflected_slab outer iter did not converge"
    )
    assert abs(res.value - case.truth_value) < tol, (
        f"{case.case_id}: fn_reflected_slab tau="
        f"{res.value:.5f} mfp vs truth "
        f"{case.truth_value} ({case.truth_source}), "
        f"diff={abs(res.value-case.truth_value):.3e} > tol={tol:.1e}"
    )


# ═══════════════════════════════════════════════════════════════════
# Truth gate — closed-sphere k_inf (trajectory_resolvent V_α1)
# ═══════════════════════════════════════════════════════════════════


@pytest.mark.l1
@pytest.mark.parametrize(
    "case", CLOSED_SPHERE_KINF_CASES, ids=lambda c: c.case_id,
)
def test_trajectory_resolvent_sphere_closed_matches_kinf(
    case: CrossMethodCase,
):
    r"""Closed-sphere (α=1) trajectory_resolvent gives ``k_eff = k_inf``
    to machine precision.

    V_α1 algebraic identity: at α=1 (perfect specular BC) the
    closed sphere has rank-1 isotropic eigenmode and ``k_eff =
    νΣ_f / Σ_a`` independent of R. This is the **multi-group
    cross-method gate** scaffolding — extending to 2G+ with this
    adapter is the natural next step (closed-sphere α=1 has
    ``k_eff = k_inf`` for any group structure).
    """
    adapter = TrajectoryResolventSphereClosedAdapter()
    if adapter.name not in case.tolerances:
        pytest.skip(
            f"trajectory_resolvent_sphere_closed not opted in for "
            f"{case.case_id!r}"
        )
    res = adapter.solve(case)
    tol = case.tolerance_for(adapter)
    assert res.metadata["converged"], (
        f"{case.case_id}: closed-sphere trajectory_resolvent did not converge"
    )
    assert abs(res.value - case.truth_value) < tol, (
        f"{case.case_id}: trajectory_resolvent_sphere_closed k_inf="
        f"{res.value:.16e} vs analytic k_inf={case.truth_value:.16e}, "
        f"diff={abs(res.value-case.truth_value):.3e} > tol={tol:.1e}"
    )


# ═══════════════════════════════════════════════════════════════════
# Cross-method agreement gates — L4 with L1 backing
# ═══════════════════════════════════════════════════════════════════


@pytest.mark.l1
@pytest.mark.slow
@pytest.mark.parametrize(
    "case", BARE_CRITICAL_SLAB_CASES, ids=lambda c: c.case_id,
)
def test_fn_slab_vs_trajectory_resolvent_slab(case: CrossMethodCase):
    r"""**Cross-method (L4) agreement gate**: F_N slab and
    trajectory_resolvent slab agree on the bare-critical configuration.

    Both methods are independently L1-verified against the same
    Sood/KLL truth (see the truth-gate tests above). This test
    pins their **mutual agreement** to within the larger of their
    two truth tolerances. The agreement is the
    structural-independence cross-check between two methods that
    share only ``numpy``/``scipy`` above the trusted-library line:

    * **F_N**: Case singular-eigenfunction representation, Wiener-
      Hopf factorisation, collocation system on µ-moments.
    * **trajectory_resolvent**: bouncing-trajectory angle-resolved
      Green's function, fixed-point iteration on surface inflow.

    These representations agree on the same physics by genuinely
    disjoint mathematical paths — their agreement is L1-strength
    evidence for either method, NOT L4 cross-implementation. The
    L4 tag here is the bookkeeping convention from
    :doc:`/skills/vv-principles` §"V&V level taxonomy" — every
    code-to-code comparison wears L4 even when the underlying
    methods are L1-grade independent.
    """
    fn = FNSlabAdapter()
    tr = TrajectoryResolventSlabAdapter()
    if not (
        fn.name in case.tolerances and tr.name in case.tolerances
    ):
        pytest.skip(
            f"Both fn_slab and trajectory_resolvent_slab must be "
            f"opted in for cross-check; case={case.case_id!r} "
            f"opts in: {list(case.tolerances)}"
        )

    # F_N predicts the critical half-thickness; trajectory_resolvent
    # at that half-thickness predicts k_eff = 1.0. The agreement is
    # the trajectory_resolvent k_eff vs 1.0, with tolerance set by
    # the larger of the two truth tolerances (= the
    # trajectory_resolvent floor).
    res_fn = fn.solve(case)
    # Build a trajectory_resolvent solve at F_N's predicted thickness
    # (NOT at the case's published truth, to test agreement of the
    # methods themselves rather than a triple agreement with truth).
    # The trajectory_resolvent slab adapter reads its width off
    # structured_geometry.domain_extent_cm; we shadow the registry
    # case with an inline structured_geometry whose extent reflects
    # F_N's prediction (full slab = 2 × half-thickness).
    case_at_fn_thickness = _shadow_with_thickness_mfp(
        case, a_critical_mfp=float(res_fn.value)
    )
    res_tr = tr.solve(case_at_fn_thickness)

    tol = agreement_tolerance(case, fn, tr)
    diff = abs(res_tr.value - 1.0)
    assert diff < tol, (
        f"{case.case_id}: F_N predicts a_c={res_fn.value:.10f} mfp; "
        f"trajectory_resolvent at that half-thickness gives k_eff="
        f"{res_tr.value:.8f}, |k-1|={diff:.3e} > pairwise tol "
        f"max({case.tolerance_for(fn):.1e}, "
        f"{case.tolerance_for(tr):.1e}) = {tol:.1e}. "
        f"Both methods backed by {case.truth_source}."
    )


@pytest.mark.l1
@pytest.mark.parametrize(
    "case", BARE_CRITICAL_SPHERE_CASES, ids=lambda c: c.case_id,
)
def test_fn_sphere_vs_trajectory_resolvent_sphere(case: CrossMethodCase):
    r"""**Cross-method (L4) agreement gate**: F_N sphere (Siewert-Thomas
    1986 Wiener-Hopf via Case eigenfunctions) vs trajectory_resolvent
    sphere (Sanchez 1986 bouncing-trajectory operator).

    Both methods are independently L1-verified against Sood/KLL
    truth. This test pins their pairwise agreement.
    """
    fn = FNSphereAdapter()
    tr = TrajectoryResolventSphereAdapter()
    if not (
        fn.name in case.tolerances and tr.name in case.tolerances
    ):
        pytest.skip(
            f"Both adapters must be opted in for cross-check; "
            f"case={case.case_id!r} opts in: {list(case.tolerances)}"
        )

    res_fn = fn.solve(case)
    # The trajectory_resolvent sphere adapter reads its radius off
    # structured_geometry.domain_extent_cm; shadow with an inline
    # geometry at F_N's predicted radius.
    case_at_fn_radius = _shadow_with_thickness_mfp(
        case, R_critical_mfp=float(res_fn.value)
    )
    res_tr = tr.solve(case_at_fn_radius)

    tol = agreement_tolerance(case, fn, tr)
    diff = abs(res_tr.value - 1.0)
    assert diff < tol, (
        f"{case.case_id}: F_N predicts R_c={res_fn.value:.10f} mfp; "
        f"trajectory_resolvent at that radius gives k_eff="
        f"{res_tr.value:.8f}, |k-1|={diff:.3e} > pairwise tol "
        f"max({case.tolerance_for(fn):.1e}, "
        f"{case.tolerance_for(tr):.1e}) = {tol:.1e}. "
        f"Both methods backed by {case.truth_source}."
    )


# ═══════════════════════════════════════════════════════════════════
# Coverage diagnostics (foundation) — print the agreement matrix
# ═══════════════════════════════════════════════════════════════════


@pytest.mark.foundation
def test_coverage_matrix_diagnostic(capsys):
    """Print the (case × adapter) agreement matrix for visibility.

    This test always passes — it exists to surface the cross-method
    coverage in the test output for the agreement-matrix renderer.
    The matrix shows which (case, adapter) combinations are opted
    in via tolerances; cells marked ``--`` are not opted in.
    """
    adapter_names = sorted(ADAPTERS_BY_NAME)
    header = f"{'case_id':<42} | " + " | ".join(
        f"{n[:14]:>14}" for n in adapter_names
    )
    print()
    print(header)
    print("-" * len(header))
    for case in ALL_CASES:
        row = f"{case.case_id:<42} | " + " | ".join(
            f"{case.tolerances.get(n, '--'):>14}"
            if isinstance(case.tolerances.get(n), float)
            else f"{'--':>14}"
            for n in adapter_names
        )
        print(row)

    # Sanity: every case has at least one float tolerance.
    for case in ALL_CASES:
        assert any(
            isinstance(t, float) for t in case.tolerances.values()
        ), case.case_id
