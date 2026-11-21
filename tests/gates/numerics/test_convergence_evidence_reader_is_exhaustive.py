"""The ``Evidence`` readers are exhaustive: the projection ``measured_value`` and the ``balance`` clause of ``warn_if_unconverged``.

``Evidence = Measured | Asserted | NotApplicable | NotYet`` is a closed union, so its reader is an
exhaustive ``match`` ending in ``assert_never`` (it survives ``python -O``). Three arms are LEGITIMATE
silent arms (``Asserted``, ``NotApplicable``, ``NotYet``); the red inputs are values OUTSIDE the union
(``NewKind()``, ``None``), which both readers read silently before the change.

Rungs: ``measured_value`` is the foundation; the converged leg and the stalled leg read it.

What reddens each leg today (X1): ``measured_value`` does not exist (``AttributeError``);
``NewKind()`` and ``None`` return silently from ``warn_if_unconverged``, on a stalled record and,
because the early ``fully_converged`` return precedes any read, on a converged one.

Mutations (in process; never patch ``orpheus.numerics.evidence``, ``convergence`` holds its own binding):
``monkeypatch.setattr(convergence, "measured_value", lambda e: e.value if isinstance(e, Measured) else None)``
reddens every foreign-kind leg; moving the call after the ``fully_converged`` guard reddens the
converged leg; replacing ``assert_never(e)`` by ``raise AssertionError()`` in ``measured_value`` keeps
every runtime leg green and reddens only the static predicate ``_reads_exhaustively``.
"""

from __future__ import annotations

import ast
import inspect
import textwrap
import typing
import warnings

import pytest

from orpheus.numerics import convergence, evidence
from orpheus.numerics.convergence import IterationBudget, IterationRecord, StoppingCriterion
from orpheus.numerics.evidence import Evidence, Measured, NotApplicable
from tests.gates.numerics._evidence_arms import ARM_FIXTURES, NewKind

pytestmark = pytest.mark.foundation

_ARM_NAMES = {c.__name__ for c in typing.get_args(Evidence)}
_LEGITIMATE_SILENT = [a for a in ARM_FIXTURES if a is not Measured]


def _stalled() -> IterationRecord:
    return IterationRecord(
        label="inner(probe)",
        criteria=(StoppingCriterion(name="residual", trajectory=tuple(1e-3 for _ in range(40)), tolerance=1e-10),),
        budget=IterationBudget(40, "max_inner"),
        iterations_run=40,
    )


def _converged() -> IterationRecord:
    record = IterationRecord(
        label="inner(probe)",
        criteria=(StoppingCriterion(name="residual", trajectory=(1e-3, 1e-12), tolerance=1e-10),),
        budget=IterationBudget(40, "max_inner"),
        iterations_run=2,
    )
    if not record.fully_converged:
        raise AssertionError("the converged fixture must be fully converged")
    return record


def _messages(record: IterationRecord, ev: object) -> list[str]:
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        convergence.warn_if_unconverged(record, where="probe", balance_defect=ev)  # type: ignore[arg-type]
    return [str(w.message) for w in caught]


def _reads_exhaustively(source: str) -> bool:
    """True when the ``match`` of ``source`` ends in ``case _: assert_never(<the parameter>)`` and its other
    patterns name exactly the arms of ``Evidence``."""
    func = ast.parse(textwrap.dedent(source)).body[0]
    if not isinstance(func, ast.FunctionDef) or not func.args.args:
        return False
    param = func.args.args[0].arg
    matches = [n for n in ast.walk(func) if isinstance(n, ast.Match)]
    if len(matches) != 1:
        return False
    cases = matches[0].cases
    last = cases[-1]
    if not (isinstance(last.pattern, ast.MatchAs) and last.pattern.pattern is None and last.pattern.name is None):
        return False
    if len(last.body) != 1 or not isinstance(last.body[0], ast.Expr):
        return False
    call = last.body[0].value
    if not (
        isinstance(call, ast.Call)
        and isinstance(call.func, ast.Name)
        and call.func.id == "assert_never"
        and len(call.args) == 1
        and isinstance(call.args[0], ast.Name)
        and call.args[0].id == param
    ):
        return False
    names: set[str] = set()
    for case in cases[:-1]:
        for node in ast.walk(case.pattern):
            if isinstance(node, ast.MatchClass) and isinstance(node.cls, ast.Name):
                names.add(node.cls.id)
    return names == _ARM_NAMES


def test_the_arm_fixtures_cover_the_union() -> None:
    """A fifth arm added to the alias reddens this until ``ARM_FIXTURES`` gains its fixture."""
    assert set(ARM_FIXTURES) == set(typing.get_args(Evidence))


class TestMeasuredValue:
    def test_a_measured_arm_projects_to_its_number(self) -> None:
        assert evidence.measured_value(Measured(0.5)) == 0.5  # type: ignore[attr-defined]

    @pytest.mark.parametrize("arm", _LEGITIMATE_SILENT, ids=lambda a: a.__name__)
    def test_the_arms_that_carry_no_number_project_to_none(self, arm: type) -> None:
        assert evidence.measured_value(ARM_FIXTURES[arm]) is None  # type: ignore[attr-defined]

    @pytest.mark.parametrize("foreign", [NewKind(), None], ids=["NewKind", "None"])
    def test_a_value_outside_the_union_is_refused(self, foreign: object) -> None:
        with pytest.raises(AssertionError):
            evidence.measured_value(foreign)  # type: ignore[attr-defined]


class TestTheStaticLaw:
    def test_measured_value_ends_in_assert_never(self) -> None:
        assert _reads_exhaustively(inspect.getsource(evidence.measured_value))  # type: ignore[attr-defined]

    def test_control_the_predicate_rejects_a_raise_in_place_of_assert_never(self) -> None:
        """The mutation that keeps every runtime leg green and silences pyright on a fifth arm."""
        source = inspect.getsource(evidence.measured_value)  # type: ignore[attr-defined]
        assert "assert_never(e)" in source
        assert not _reads_exhaustively(source.replace("assert_never(e)", "raise AssertionError()"))


class TestWarnIfUnconverged:
    def test_a_measured_balance_defect_renders(self) -> None:
        (message,) = _messages(_stalled(), Measured(1e-3))
        assert "balance defect" in message

    @pytest.mark.parametrize("arm", _LEGITIMATE_SILENT, ids=lambda a: a.__name__)
    def test_the_arms_without_a_number_leave_the_clause_absent(self, arm: type) -> None:
        (message,) = _messages(_stalled(), ARM_FIXTURES[arm])
        (reference,) = _messages(_stalled(), ARM_FIXTURES[NotApplicable])
        assert message == reference
        assert "balance defect" not in message

    @pytest.mark.parametrize("foreign", [NewKind(), None], ids=["NewKind", "None"])
    def test_a_value_outside_the_union_is_refused_on_a_stalled_record(self, foreign: object) -> None:
        with pytest.raises(AssertionError):
            _messages(_stalled(), foreign)

    @pytest.mark.parametrize("foreign", [NewKind(), None], ids=["NewKind", "None"])
    def test_a_value_outside_the_union_is_refused_on_a_converged_record_too(self, foreign: object) -> None:
        """The read precedes the ``fully_converged`` return: a converged exit must not hide a foreign kind."""
        with pytest.raises(AssertionError):
            _messages(_converged(), foreign)

    @pytest.mark.parametrize("arm", list(ARM_FIXTURES), ids=lambda a: a.__name__)
    def test_a_converged_record_is_silent_for_every_arm(self, arm: type) -> None:
        assert _messages(_converged(), ARM_FIXTURES[arm]) == []
