"""``warn_if_gauge_freedom`` reads its ``correction: Evidence`` exhaustively (via ``measured_value``).

Companion of ``tests/gates/numerics/test_convergence_evidence_reader_is_exhaustive.py`` (a numerics-directory
test may not import ``orpheus.sn``, so this reader's gate lives here).

First reds today (X1): ``NewKind()`` and ``None`` return silently, and with ``undetermined=True`` the
branch warns and returns BEFORE ``correction`` is read, so a foreign kind is silent there too.

Mutations (in process, on the reader's own binding, never on ``orpheus.numerics.evidence``):
``monkeypatch.setattr(K, "measured_value", lambda e: e.value if isinstance(e, Measured) else None)``
reddens the foreign-kind legs; moving the ``measured_value`` call after the ``undetermined`` branch
reddens the ``undetermined=True`` legs.
"""

from __future__ import annotations

import types
import warnings

import pytest

from orpheus.numerics.evidence import Measured
from orpheus.sn.operators import loss_kernel_gauge as K
from tests.gates.numerics._evidence_arms import ARM_FIXTURES, SILENT_ARMS, NewKind

pytestmark = pytest.mark.foundation



def _warnings(monkeypatch: pytest.MonkeyPatch, ev: object, *, undetermined: bool) -> list[warnings.WarningMessage]:
    verdict = types.SimpleNamespace(undetermined=undetermined, present=True, because="stub")
    monkeypatch.setattr(K, "gauge_freedom", lambda p: verdict)
    problem = types.SimpleNamespace(ndim=2)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        K.warn_if_gauge_freedom(problem, ev, where="probe")  # type: ignore[arg-type]
    return list(caught)


def test_a_measured_correction_above_the_floor_warns_once(monkeypatch: pytest.MonkeyPatch) -> None:
    assert len(_warnings(monkeypatch, Measured(0.5), undetermined=False)) == 1


@pytest.mark.parametrize("arm", SILENT_ARMS, ids=lambda a: a.__name__)
def test_the_arms_without_a_number_are_silent(monkeypatch: pytest.MonkeyPatch, arm: type) -> None:
    assert _warnings(monkeypatch, ARM_FIXTURES[arm], undetermined=False) == []


@pytest.mark.parametrize("undetermined", [False, True], ids=["classified", "undetermined"])
@pytest.mark.parametrize("foreign", [NewKind(), None], ids=["NewKind", "None"])
def test_a_value_outside_the_union_is_refused_whichever_the_verdict(
    monkeypatch: pytest.MonkeyPatch, foreign: object, undetermined: bool
) -> None:
    with pytest.raises(AssertionError):
        _warnings(monkeypatch, foreign, undetermined=undetermined)


def test_an_undetermined_closure_still_warns_for_a_legitimate_arm(monkeypatch: pytest.MonkeyPatch) -> None:
    """Control: the first-statement read does not silence the undetermined branch for the four arms."""
    assert len(_warnings(monkeypatch, ARM_FIXTURES[SILENT_ARMS[0]], undetermined=True)) == 1
