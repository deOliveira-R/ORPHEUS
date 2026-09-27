r"""The verification floor, its tolerance rule and its xfail mark: software invariants of the cross-check harness.

:mod:`tests.gates.sn.verification.analytical._certified_agreement` decides
whether an SN-against-reference comparison certifies anything, and
:mod:`tests.gates.derivations._trajectory_resolvent_ladders` derives every
bound and tolerance those comparisons use. Both are consulted by slow rows
whose failures would read as physics; these rows pin the harness itself,
fast, one per branch.
"""

from __future__ import annotations

import importlib

import pytest

from tests.gates.derivations._trajectory_resolvent_ladders import (
    alternating_error,
    ceil_one_significant_figure,
    richardson_error,
    tolerance_for,
)
from tests.gates.sn.verification.analytical._certified_agreement import (
    assert_record,
    awaits_cylinder_bound,
    certify_agreement,
)


@pytest.mark.foundation
def test_certificate_reaches_every_branch() -> None:
    """No bound; a bound above a tenth of the tolerance; a reading at the tolerance; agreement."""
    absent = certify_agreement("x", 1e-5, 1e-3, None)
    assert not absent.floor_holds and not absent.agrees
    with pytest.raises(AssertionError, match="carries no certified error bound"):
        absent.require()

    too_loose = certify_agreement("x", 1e-5, 1e-3, 1.1e-4)
    assert not too_loose.floor_holds
    with pytest.raises(AssertionError, match="exceeds a tenth of the tolerance"):
        too_loose.require()

    disagrees = certify_agreement("x", 1e-3, 1e-3, 1e-4)
    assert disagrees.floor_holds and not disagrees.agrees
    with pytest.raises(AssertionError, match="above the tolerance"):
        disagrees.require()

    agrees = certify_agreement("x", 9.9e-4, 1e-3, 1e-4)
    assert agrees.floor_holds and agrees.agrees
    agrees.require()


@pytest.mark.foundation
def test_tolerance_rule() -> None:
    """``tolerance_for`` is the smallest one-significant-figure T with T >= 10 b and T >= 2 (e + b), or 2.5 e with no bound."""
    assert ceil_one_significant_figure(1.5e-3) == pytest.approx(2e-3)
    assert ceil_one_significant_figure(3.2e-3) == pytest.approx(4e-3)
    assert ceil_one_significant_figure(4e-3) == pytest.approx(4e-3)
    assert tolerance_for(1.5e-5, 3.2e-4) == pytest.approx(4e-3)       # the floor governs
    assert tolerance_for(7.3e-3, 1.4e-3) == pytest.approx(2e-2)       # both errors govern
    assert tolerance_for(5.9e-4, None) == pytest.approx(2e-3)         # reference assumed at the floor
    for residual, bound in ((1.5e-5, 3.2e-4), (7.3e-3, 1.4e-3), (1.43e-3, 4.6e-5)):
        tolerance = tolerance_for(residual, bound)
        assert certify_agreement("x", 0.0, tolerance, bound).floor_holds


@pytest.mark.foundation
def test_ladder_error_estimates() -> None:
    """Second order at refinement 2 inflates a step by 4/3; an alternating sequence's error is within its last step."""
    assert richardson_error(-3e-4, 2.0, 2) == pytest.approx(4e-4)
    assert richardson_error(3e-4, 2.0, 1) == pytest.approx(6e-4)
    assert alternating_error(-1.3e-4) == pytest.approx(1.3e-4)


@pytest.mark.foundation
def test_record_moves_red() -> None:
    """``assert_record`` passes inside its band and raises outside it (relative and absolute legs)."""
    recorded = {"k": 1.2, "gap": 1e-3}
    assert_record({"k": 1.2 * (1 + 1e-6), "gap": 1e-3 + 1e-6}, recorded, 1e-5, relative=frozenset({"k"}))
    with pytest.raises(AssertionError, match="record 'k' moved"):
        assert_record({"k": 1.2 * (1 + 1e-4), "gap": 1e-3}, recorded, 1e-5, relative=frozenset({"k"}))
    with pytest.raises(AssertionError, match="record 'gap' moved"):
        assert_record({"k": 1.2, "gap": 1.1e-3}, recorded, 1e-5, relative=frozenset({"k"}))


#: Every row that needs the cylinder reference's bound, by module.
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
    """The five rows (seven collected cases) carry the ONE shared mark: strict, expecting only an ``AssertionError``.

    A non-strict xfail would absorb an XPASS, and one expecting any exception
    would absorb a crash before the floor assertion (``vv-principles`` mode
    8(4)). The census runs the other way too: no other function in those
    modules carries a #516 xfail.
    """
    assert awaits_cylinder_bound.kwargs["strict"] is True
    assert awaits_cylinder_bound.kwargs["raises"] is AssertionError
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
