r"""The reference cross-check's harness: its tolerance rules, its ladder estimates, its RECORD, its xfail mark.

:mod:`tests.gates.sn.verification.analytical._aba_reference` owns the A|B|A
problem, the reference, the RECORD table and the strict-xfail mark;
:mod:`tests.gates.derivations._ladder_rules` holds the rules and
:mod:`tests.gates.derivations._characteristic_ladders` the measured ladders
every tolerance is derived from (the characteristic reference's since P1 step
(d); the trajectory resolvent's before, whose tables remain in
``_trajectory_resolvent_ladders`` for the corroboration rows). All are
consulted by slow rows whose failures would read as physics; these rows pin
the harness itself, fast, one per branch (``foundation``).

Until #405 P2 step 7b.2.3 this file was ``test_certified_agreement.py`` and
also pinned the four branches of the test-side ``AgreementCertificate``; that
type retired with ``certify_agreement`` when ``verify_agreement`` and
``compare_uncertified`` (``orpheus.reference.verification``, gated in
``tests/gates/reference/``) took its place, and its row went with it.
"""

from __future__ import annotations

import importlib
from fractions import Fraction

import pytest

from orpheus.reference.verification import ReferenceNotValid
from tests.gates.derivations import _characteristic_ladders as ladders
from tests.gates.derivations._ladder_rules import (
    alternating_error,
    ceil_one_significant_figure,
    geometric_error,
    richardson_error,
    summed_tolerance,
    tolerance_for,
)
from tests.gates.sn.verification.analytical._aba_reference import (
    assert_record,
    awaits_cylinder_bound,
    truncated,
)


@pytest.mark.foundation
def test_tolerance_rule() -> None:
    """``tolerance_for`` is the smallest one-significant-figure T with T >= 10 b and T >= 2 (e + b), or 2.5 e with no
    estimate. The floor property T >= 10 b is asserted exactly (it was read through the retired certificate)."""
    assert ceil_one_significant_figure(1.5e-3) == pytest.approx(2e-3)
    assert ceil_one_significant_figure(3.2e-3) == pytest.approx(4e-3)
    assert ceil_one_significant_figure(4e-3) == pytest.approx(4e-3)
    assert tolerance_for(1.5e-5, 3.2e-4) == pytest.approx(4e-3)       # the floor governs
    assert tolerance_for(7.3e-3, 1.4e-3) == pytest.approx(2e-2)       # both errors govern
    assert tolerance_for(5.9e-4, None) == pytest.approx(2e-3)         # reference assumed at the floor
    for residual, estimate in ((1.5e-5, 3.2e-4), (7.3e-3, 1.4e-3), (1.43e-3, 4.6e-5)):
        assert Fraction(tolerance_for(residual, estimate)) >= 10 * Fraction(estimate)


@pytest.mark.foundation
def test_ladder_error_estimates() -> None:
    """Second order at refinement 2 inflates a step by 4/3; an alternating sequence's error is within its last step; a
    geometric one's is its step over one minus the contraction, refused when the steps do not contract; the summed
    band is the two errors' sum rounded up, with no floor."""
    assert richardson_error(-3e-4, 2.0, 2) == pytest.approx(4e-4)
    assert richardson_error(3e-4, 2.0, 1) == pytest.approx(6e-4)
    assert alternating_error(-1.3e-4) == pytest.approx(1.3e-4)
    assert geometric_error(-1e-6, 1e-5) == pytest.approx(1e-6 / 0.9)
    with pytest.raises(ValueError, match="the steps do not contract"):
        geometric_error(1e-6, 1e-6)
    assert summed_tolerance(6.0e-4, 6.05e-8) == pytest.approx(7e-4)
    assert summed_tolerance(3.0e-4, 1.1e-9) == pytest.approx(4e-4)


@pytest.mark.foundation
def test_each_reference_estimate_covers_the_finest_measured_rung() -> None:
    """Each estimate of the characteristic reference's error at a p = 5 working point is at least its distance to the
    finest rung measured (p = 8), and within a factor 2 of it: the estimator's model (geometric steps; for the A|B|A
    sphere's shape, pairs of degrees) checked against the ladder's own tail.

    The data are ``_characteristic_ladders``' ``[M]`` tables, so this row reads no reference. Declared: p = 8 is not
    the limit, so "covers" holds to the p = 8 rung's own error: ``[M]`` that error is at least 1.3 orders below each
p = 5 estimate (the A|B|A sphere's shape, its least, p = 7 -> 8 being 4.7e-7 against 9.4e-6). The row reds on
    a stride-1 shape estimate (3.6e-4 against 9.0e-6: the factor-2 leg) and on a stride-2 k estimate of the A|B|A
    sphere (2.77e-11 against 2.80e-11: the cover leg).
    """
    finest = 8
    pairs = {
        "aba_sphere": (ladders.reference_error("aba_sphere"), ladders.ABA_SPHERE_K),
        "partial_reflector_slab": (ladders.reference_error("partial_reflector_slab"), ladders.PARTIAL_REFLECTOR_K["slab"]),
        "partial_reflector_sphere": (ladders.reference_error("partial_reflector_sphere"), ladders.PARTIAL_REFLECTOR_K["sphere"]),
    }
    assert len(pairs) == 3
    for fixture, (estimate, k) in pairs.items():
        p = ladders.WORKING_DEGREE[fixture]
        distance = abs(k[p] - k[finest]) / abs(k[finest])
        assert distance <= estimate <= 2.0 * distance, (fixture, estimate, distance)
    shape = ladders.aba_sphere_shape_error()
    distance = ladders.ABA_SPHERE_SHAPE_STEP[(ladders.WORKING_DEGREE["aba_sphere"], finest)]
    assert distance <= shape <= 2.0 * distance, (shape, distance)


@pytest.mark.foundation
def test_truncation_only_tightens() -> None:
    """``truncated`` cuts toward zero to three figures and never lands above the value, so τ_rel × truncated(x) never
    exceeds τ_rel × x. 1.4 is the witness of qa F4 (the floor-times-power spelling returned 1.4000000000000001), and
    200 000 seeded draws over six decades of each sign are the population (it failed 45 035 of 2e6 at HEAD)."""
    import numpy as np

    assert truncated(1.4) == 1.4 and truncated(1.381079639349644) == 1.38 and truncated(0.24937654699132977) == 0.249
    rng = np.random.default_rng(20261003)
    draws = 10.0 ** rng.uniform(-3.0, 3.0, 100_000)
    for value in (*draws.tolist(), *(-draws).tolist(), 1.4, -1.4, 999.9999, 1e-300):
        cut = truncated(value)
        assert abs(cut) <= abs(value) and abs(value) - abs(cut) < 1e-2 * abs(value) and (cut > 0) == (value > 0), value


@pytest.mark.foundation
def test_record_moves_red() -> None:
    """``assert_record`` passes inside its band and raises outside it (relative and absolute legs)."""
    recorded = {"k": 1.2, "gap": 1e-3}
    assert_record({"k": 1.2 * (1 + 1e-6), "gap": 1e-3 + 1e-6}, recorded, 1e-5, relative=frozenset({"k"}))
    with pytest.raises(AssertionError, match="record 'k' moved"):
        assert_record({"k": 1.2 * (1 + 1e-4), "gap": 1e-3}, recorded, 1e-5, relative=frozenset({"k"}))
    with pytest.raises(AssertionError, match="record 'gap' moved"):
        assert_record({"k": 1.2, "gap": 1.1e-3}, recorded, 1e-5, relative=frozenset({"k"}))


#: Every row whose comparison needs the cylinder reference's certificate, by module.
_AWAITING_ROWS = {
    "tests.gates.sn.verification.analytical.test_phase_c_crosscheck": (
        "test_cylinder_3reg_k_against_trajectory_resolvent",
        "test_cylinder_3reg_flux_shape_against_trajectory_resolvent",
    ),
    "tests.gates.sn.sweep.curvilinear.test_unified_matvec_cylinder": (
        "test_unified_cylinder_l1_mr_2g_trajectory_resolvent",
    ),
    "tests.gates.sn.verification.analytical.test_l1_standoff_slab_cylinder": (
        "test_cylinder_l1_sweep_vs_trajectory_resolvent",
        "test_cylinder_l1_refinement_against_reference",
    ),
}


@pytest.mark.foundation
def test_cylinder_bound_rows_carry_the_shared_strict_xfail() -> None:
    """The five rows (seven collected cases) carry the ONE shared mark: strict, expecting only ``ReferenceNotValid``.

    A non-strict xfail would absorb an XPASS, and one expecting any exception
    would absorb a crash before the verbs' refusal (``vv-principles`` mode
    8(4)). The census runs the other way too: no other function in those
    modules carries a #516 xfail.
    """
    assert awaits_cylinder_bound.kwargs["strict"] is True
    assert awaits_cylinder_bound.kwargs["raises"] is ReferenceNotValid
    for module_name, expected in _AWAITING_ROWS.items():
        module = importlib.import_module(module_name)
        carriers = sorted(
            name for name, obj in vars(module).items()
            if name.startswith("test_") and any(
                m.name == "xfail" and "#516" in str(m.kwargs.get("reason", ""))
                for m in getattr(obj, "pytestmark", [])
            )
        )
        assert carriers == sorted(expected), f"{module_name}: #516 xfails on {carriers}, expected {sorted(expected)}"
        for name in expected:
            marks = [m for m in getattr(module, name).pytestmark if m.name == "xfail"]
            assert marks == [awaits_cylinder_bound.mark], f"{module_name}.{name}: not the shared mark"
