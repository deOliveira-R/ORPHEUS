r"""Cross-method eigenvalue / critical-dimension regression gates.

Tests in this file exercise the
:class:`~tests.gates.cross_method.protocol.SolverAdapter` protocol over
the populated case sets in :mod:`~tests.gates.cross_method.cases`.

Relationship to direct math-heart-class construction (Phase D)
--------------------------------------------------------------

The math-heart classes
(:class:`~orpheus.derivations.continuous.fn_method.moment_space.MomentSpace`)
are constructed directly with a :class:`StructuredGeometry` plus
``materials: dict[int, Mixture]``. **The tests in this file
deliberately keep their names + per-method bodies** — pytest
collection IDs are preserved (CI / pytest-xdist contract) and the
per-method adapter classes (``FNSlabAdapter``,
``CharacteristicAdapter``, ...) stay as the
unit-conversion layer. The agreement between the adapter route and
the direct-construction route is exercised by
:mod:`tests.gates.cross_method.test_polymorphism` (foundation-tier
regression net, 5 tests).

Three classes of test:

1. **Truth gates** — each adapter reproduces its case's truth value
   to within the case's per-adapter tolerance. One parametrised
   test per (case × adapter) pair. These are the L1 backing —
   each method's agreement with closed-form truth.

2. **Cross-method agreement gates** — pairs of adapters that both
   support a case agree to within the larger of their two truth
   tolerances (per
   :func:`~tests.gates.cross_method.protocol.agreement_tolerance` —
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
  (``test_fn_sood2003_slab_xverif.py``,
  ``test_fn_sood2003_sphere_xverif.py``) tag these as L1 because
  the agreement is **L1-strength evidence** for either method
  when both methods' L1 truth-backing is established and
  structural independence is genuine. The L1 backing here is the
  per-adapter truth gates in the same file.

Slow tests
----------

None today: the characteristic adapters' 15 rows take about 14 s in
all. Use ``pytest -m "l1 and not slow"`` for the fast subset.
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

from orpheus.geometry import CoordSystem

from .adapters import (
    ADAPTERS_BY_NAME,
    CHARACTERISTIC_SLAB,
    CHARACTERISTIC_SPHERE,
    CHARACTERISTIC_SPHERE_CLOSED,
    CharacteristicAdapter,
    FNReflectedSlabAdapter,
    FNSlabAdapter,
    FNSphereAdapter,
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

    # Slab: the published critical dimension is the half-thickness and
    # the geometry's extent the FULL slab width. Sphere / cylinder: the
    # published critical dimension IS the radius.
    extent_cm = 2.0 * cd_cm if base_geom.coord is CoordSystem.CARTESIAN else cd_cm
    new_geom = replace(base_geom, breakpoints=(0.0, extent_cm), mat_ids=(0,))
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

    Backed by Sood 2003 / Grandjean-Siewert / KLL via Sood
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
    characteristic-reference counterpart is registered for them.
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

    Backed by Sood 2003 / KLL via Sood transcription. F_N
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
# Truth gates — reflected slab (fn_method only; one-sided coverage)
# ═══════════════════════════════════════════════════════════════════


@pytest.mark.l1
@pytest.mark.parametrize(
    "case", REFLECTED_SLAB_CASES, ids=lambda c: c.case_id,
)
def test_fn_reflected_slab_matches_truth(case: CrossMethodCase):
    """F_N reflected slab reproduces the case's truth ``tau_critical_mfp``.

    Backed by Sood 2003 Table 7 (problem 4) + NM 1980 Table 2 + Burkart
    1976 'Exact'. **No characteristic-reference counterpart is
    registered** — this is one-sided coverage. The characteristic
    reference poses a reflected slab (two regions, any walls), so a
    counterpart is a case-set and adapter addition, not in this task's
    scope.
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
# The characteristic reference: the successors of the (deleted) trajectory-resolvent rows
# ═══════════════════════════════════════════════════════════════════
#
# P1 step (e1b) of ``.claude/plans/characteristic_reference_architecture.md`` (the user's ruling 3 of 2026-10-10):
# each trajectory-resolvent truth and agreement row had its successor here, on the characteristic adapters,
# with the tolerance ``cases.characteristic_tolerance`` computes from the reference's own ladder and the truth's
# resolution. Step (e2) deleted the predecessors. None is slow: a slab case at the working rung
# takes about 1.5 s.

_CHARACTERISTIC_SLAB = CHARACTERISTIC_SLAB
_CHARACTERISTIC_SPHERE = CHARACTERISTIC_SPHERE
_CHARACTERISTIC_SPHERE_CLOSED = CHARACTERISTIC_SPHERE_CLOSED


def _characteristic_k_against_one(adapter: CharacteristicAdapter, case: CrossMethodCase, at: str) -> None:
    res = adapter.solve(case)
    tol = case.tolerance_for(adapter.name)
    assert abs(res.value - 1.0) < tol, (
        f"{case.case_id}: {adapter.name} k_eff={res.value!r} {at}, |k-1|={abs(res.value - 1.0):.3e} > tol={tol:.1e}. "
        f"Truth source: {case.truth_source}"
    )


@pytest.mark.l1
@pytest.mark.parametrize("case", BARE_CRITICAL_SLAB_CASES, ids=lambda c: c.case_id)
def test_characteristic_slab_matches_truth_keff_one(case: CrossMethodCase):
    """The characteristic reference's vacuum slab at the case's published critical half-thickness reads k = 1.

    Successor of ``test_trajectory_resolvent_slab_matches_truth_keff_one``. Truth: Sood / Kaper-Lindeman-Leaf
    (F_N), a structurally independent method. First red: the slab's cosine rule integrating one hemisphere only
    (``consumers/cross_method_battery.py``).
    """
    _characteristic_k_against_one(_CHARACTERISTIC_SLAB, case, f"at the truth half-thickness {case.truth_value} mfp")


@pytest.mark.l1
@pytest.mark.parametrize("case", BARE_CRITICAL_SPHERE_CASES, ids=lambda c: c.case_id)
def test_characteristic_sphere_matches_truth_keff_one(case: CrossMethodCase):
    """The characteristic reference's vacuum sphere at the case's published critical radius reads k = 1.

    Successor of ``test_trajectory_resolvent_sphere_matches_truth_keff_one``.
    """
    _characteristic_k_against_one(_CHARACTERISTIC_SPHERE, case, f"at the truth radius {case.truth_value} mfp")


@pytest.mark.l1
@pytest.mark.parametrize("case", CLOSED_SPHERE_KINF_CASES, ids=lambda c: c.case_id)
def test_characteristic_sphere_closed_matches_kinf(case: CrossMethodCase):
    r"""The characteristic reference's closed sphere (a mirror at r = R) reads :math:`k_\infty = \nu\Sigma_f/\Sigma_a`.

    Successor of ``test_trajectory_resolvent_sphere_closed_matches_kinf`` and the cross-method protocol's only
    TR-keyed case, ``closed-sphere-1G-fuelA-tauR2.5``, which would otherwise have no tolerance after step (e).
    """
    res = _CHARACTERISTIC_SPHERE_CLOSED.solve(case)
    tol = case.tolerance_for(_CHARACTERISTIC_SPHERE_CLOSED.name)
    assert res.tag == "k_inf"
    assert abs(res.value - case.truth_value) < tol, (
        f"{case.case_id}: characteristic_sphere_closed k_inf={res.value!r} vs {case.truth_value!r}, "
        f"diff={abs(res.value - case.truth_value):.3e} > tol={tol:.1e}"
    )


@pytest.mark.l1
@pytest.mark.parametrize("case", BARE_CRITICAL_SLAB_CASES, ids=lambda c: c.case_id)
def test_fn_slab_vs_characteristic_slab(case: CrossMethodCase):
    """F_N's predicted critical half-thickness, read by the characteristic slab, gives k = 1 within the pairwise
    tolerance. Successor of ``test_fn_slab_vs_trajectory_resolvent_slab``: two methods sharing nothing above the
    trusted-library line (Case eigenfunctions and collocation against lines, panels and a dense pencil)."""
    fn = FNSlabAdapter()
    res_fn = fn.solve(case)
    res = _CHARACTERISTIC_SLAB.solve(_shadow_with_thickness_mfp(case, a_critical_mfp=float(res_fn.value)))
    tol = agreement_tolerance(case, fn.name, _CHARACTERISTIC_SLAB.name)
    assert abs(res.value - 1.0) < tol, (
        f"{case.case_id}: F_N a_c={res_fn.value:.10f} mfp; the characteristic slab there reads k={res.value!r}, "
        f"|k-1|={abs(res.value - 1.0):.3e} > {tol:.1e}"
    )


@pytest.mark.l1
@pytest.mark.parametrize("case", BARE_CRITICAL_SPHERE_CASES, ids=lambda c: c.case_id)
def test_fn_sphere_vs_characteristic_sphere(case: CrossMethodCase):
    """F_N's predicted critical radius, read by the characteristic sphere, gives k = 1 within the pairwise tolerance.
    Successor of ``test_fn_sphere_vs_trajectory_resolvent_sphere``."""
    fn = FNSphereAdapter()
    res_fn = fn.solve(case)
    res = _CHARACTERISTIC_SPHERE.solve(_shadow_with_thickness_mfp(case, R_critical_mfp=float(res_fn.value)))
    tol = agreement_tolerance(case, fn.name, _CHARACTERISTIC_SPHERE.name)
    assert abs(res.value - 1.0) < tol, (
        f"{case.case_id}: F_N R_c={res_fn.value:.10f} mfp; the characteristic sphere there reads k={res.value!r}, "
        f"|k-1|={abs(res.value - 1.0):.3e} > {tol:.1e}"
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
