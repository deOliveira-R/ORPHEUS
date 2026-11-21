"""``verify_order`` reads the algebraic ``Evidence`` of each answer exhaustively.

Third reader of the closed union ``Evidence = Measured | Asserted | NotApplicable | NotYet`` (the other two are
``measured_value`` and its callers, gated in ``tests/gates/numerics/test_convergence_evidence_reader_is_exhaustive.py``).
``parse_member`` admits ``get_args(Evidence)``, so a fifth arm added to the alias passes the parse; the ``match`` in
``verify_order`` must then end in ``case _: assert_never(item)`` (it survives ``python -O``), or the arm falls through the
loop as if established, an optimistic default.

What reddens it today (X1): a fifth arm, simulated by widening the parse's admitted set with ``NewKind`` (in process,
on the module's own ``get_args`` binding), is accepted and the verdict is returned instead of refused.

Mutation: replacing ``case _: assert_never(item)`` by ``case _: continue`` reddens the foreign-arm leg; the four
legitimate arms are unchanged.
"""

from __future__ import annotations

from typing import Any, get_args

import pytest

from orpheus.numerics.evidence import Asserted, Evidence, Measured
from orpheus.reference import verification as V
from tests.gates.numerics._evidence_arms import ARM_FIXTURES, NewKind
from tests.gates.reference import _step5 as s5
from tests.gates.reference import _step7 as s7

pytestmark = pytest.mark.foundation


def _order(errors: list[Any]) -> Any:
    answers = [(h, s7.FixedAnswer({s5.eigenvalue(): 1.0 + 0.5 * h**2})) for h in (0.4, 0.2, 0.1)]
    return s7.verify_order(answers, s5.eigenvalue(), s7.reference_with_bound(1.0, 0.0), 1e-3, errors, 2.0, 0.1)


def test_a_fifth_arm_is_refused_not_read_as_established(monkeypatch: pytest.MonkeyPatch) -> None:
    arms = get_args(Evidence)
    monkeypatch.setattr(V, "get_args", lambda _alias: (*arms, NewKind))
    with pytest.raises(AssertionError):
        _order([Measured(0.0), NewKind(), Measured(0.0)])


def test_control_the_widened_parse_admits_the_established_arms(monkeypatch: pytest.MonkeyPatch) -> None:
    """The widening alone does not refuse: with only the four arms the verdict is returned."""
    arms = get_args(Evidence)
    monkeypatch.setattr(V, "get_args", lambda _alias: (*arms, NewKind))
    assert _order([Measured(0.0), Asserted(1e-9, "t"), Measured(0.0)]) is not None


@pytest.mark.parametrize("arm", [a for a in ARM_FIXTURES if a is not Measured and a is not Asserted], ids=lambda a: a.__name__)
def test_the_unestablished_arms_are_still_refused_as_unestablished(arm: type) -> None:
    with pytest.raises(V.Unestablished):
        _order([Measured(0.0), ARM_FIXTURES[arm], Measured(0.0)])
