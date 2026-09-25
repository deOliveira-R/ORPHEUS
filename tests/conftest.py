"""Shared fixtures + V&V auto-tagging for ORPHEUS tests.

Every pytest collection pass populates :data:`tests._harness.registry.TEST_REGISTRY`
with one :class:`~tests._harness.registry.TestMetadata` entry per item.
Each entry carries the resolved V&V level, the source of that resolution
(``"explicit"`` / ``"class-name"`` / ``"func-name"`` /
``"case"`` / ``"unmarked"``), the Sphinx equation labels the test
verifies, the failure-mode / error-catalog tags it catches, and the
``VerificationCase`` names it parametrizes over.

The audit CLI (``python -m tests._harness.audit``) and a future
Sphinx verification-matrix generator read the registry directly.

See ``docs/theory/verification/harness.rst`` for the design rationale
and the full contributor guide.
"""

from __future__ import annotations

import re
from typing import Any

import pytest

from orpheus.derivations.common.withdrawal import (
    RUN_WITHDRAWN_VARIABLE,
    AllWithdrawals,
    GeneratorWithdrawn,
    Withdrawal,
    lifted_withdrawals,
)
from orpheus.derivations.reference_values import get as get_reference
from tests._harness import registry
from tests._harness.registry import TestMetadata


def pytest_addoption(parser) -> None:
    """``--capture-baseline`` — write (not assert) pre-carve snapshots.

    Used by the Wave O (#208) O.4a.2 BC-extraction matvec gate
    (``tests/gates/sn/operators/test_bc_extraction_matvec.py``): when present the
    snapshot tests WRITE the pre-extraction matvec output and skip the
    assert; absent (the default, incl. the post-carve gate) they READ the
    committed snapshots and assert byte-identity. ``pytest_addoption`` only
    fires from a ROOT ``conftest.py``, which is why it lives here rather
    than in the test module.
    """
    parser.addoption(
        "--capture-baseline",
        action="store_true",
        default=False,
        help="Wave O O.4a.2: write the pre-carve matvec snapshots instead "
             "of asserting against them (run BEFORE the BC-extraction carve).",
    )


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def ref():
    """Access analytical reference values by name.

    Usage in tests::

        def test_something(ref):
            case = ref("homo_1eg")
            assert abs(result.k_inf - case.k_inf) < 1e-10
    """
    return get_reference


# ---------------------------------------------------------------------------
# V&V marker auto-tagging + registry population
# ---------------------------------------------------------------------------

_LEVEL_CLASS_RE = re.compile(r"TestL([0-3])")
_LEVEL_FUNC_RE = re.compile(r"(?:^|_)l([0-3])(?:_|$)")

# Markers the harness recognises as V&V-level tags. ``l0``..``l3`` are
# the physics ladder (Cardinal Rule 4); ``foundation`` is the orthogonal
# software-invariant bucket described in
# ``docs/theory/verification/harness.rst``. All share the same resolution
# precedence and the same audit reporting path.
_LEVEL_MARKERS = ("l0", "l1", "l2", "l3", "foundation")


def _marker_to_level(marker_name: str) -> str:
    """Convert a pytest marker name to its registry ``VVLevel`` value.

    ``l0`` -> ``L0``, ``foundation`` -> ``foundation``. The asymmetry
    (uppercase for the ladder, lowercase for foundation) is load-bearing:
    the L0..L3 values are sorted and compared as strings elsewhere, and
    ``foundation`` deliberately sorts below them so the "numerically
    highest wins" tiebreak in :func:`_existing_level` never promotes a
    foundation tag over a conflicting physics-level tag. If both markers
    are stacked, the physics level wins and the foundation marker is
    surfaced as the conflict.
    """
    return marker_name.upper() if marker_name.startswith("l") else marker_name


def _existing_level(item: pytest.Item) -> tuple[str | None, bool]:
    """Return (level_str, was_already_marked).

    Inspects already-applied markers, including file-level ``pytestmark``
    and class-level markers. Recognises
    both the L0..L3 physics ladder and the orthogonal ``foundation``
    marker. If two different markers are present, the numerically highest
    wins (matching pytest's precedence for stacked markers); a
    foundation marker stacked with any L<N> marker yields the L<N>,
    never foundation (the physics ladder dominates because it is the
    stronger claim). A warning is emitted so duplicate tagging surfaces
    early.
    """
    present = [m.name for m in item.iter_markers() if m.name in _LEVEL_MARKERS]
    if not present:
        return None, False
    if len(set(present)) > 1:
        # Surface the conflict; choose deterministically. Sort ensures
        # ``foundation`` < ``l0`` < ``l1`` < ``l2`` < ``l3`` alphabetically,
        # and we pick the last (highest) entry — so an L<N> always wins
        # over foundation, which is the desired dominance.
        chosen = sorted(set(present))[-1]
        item.warn(
            pytest.PytestUnknownMarkWarning(
                f"{item.nodeid} has conflicting V&V level markers "
                f"{sorted(set(present))}; using {chosen!r}"
            )
        )
        return _marker_to_level(chosen), True
    return _marker_to_level(present[0]), True


def _level_from_class_name(item: pytest.Item) -> str | None:
    cls = getattr(item, "cls", None)
    if cls is None:
        return None
    m = _LEVEL_CLASS_RE.search(cls.__name__)
    return f"L{m.group(1)}" if m else None


def _level_from_func_name(item: pytest.Item) -> str | None:
    name = getattr(item, "originalname", None) or item.name
    m = _LEVEL_FUNC_RE.search(name)
    return f"L{m.group(1)}" if m else None


def _resolve_case(item: pytest.Item) -> object | None:
    """Return the VerificationCase a parametrized item consumes, or None.

    Two supported shapes:

    1. ``@pytest.mark.parametrize("case", [VerificationCase(...), ...])`` —
       the parameter value is the case object itself (any object with a
       ``vv_level`` attribute qualifies).
    2. ``@pytest.mark.parametrize("case_name", ["homo_1eg", ...])`` —
       the parameter value is a string that keys into the reference
       registry via :func:`orpheus.derivations.reference_values.get`.
       Legacy shape used by most pre-harness tests.

    Returns the resolved ``VerificationCase`` instance in both cases, or
    None if the item doesn't parametrize over a case at all.
    """
    callspec = getattr(item, "callspec", None)
    if callspec is None:
        return None
    params = callspec.params
    # Shape 1: case object directly.
    case_obj = params.get("case")
    if case_obj is not None and hasattr(case_obj, "vv_level"):
        return case_obj
    # Shape 2: case_name string → registry lookup.
    case_name = params.get("case_name")
    if isinstance(case_name, str):
        try:
            return get_reference(case_name)
        except KeyError:
            return None
    return None


def _case_names_from_parametrize(item: pytest.Item) -> tuple[str, ...]:
    """Return a single-element tuple with the case name the item uses, or ()."""
    case = _resolve_case(item)
    if case is None:
        return ()
    name = getattr(case, "name", None)
    return (name,) if isinstance(name, str) else ()


def _level_from_case(item: pytest.Item) -> tuple[str | None, tuple[str, ...]]:
    """Return (level, equation_labels) inherited from the item's VerificationCase.

    Walks both the ``case`` (object) and ``case_name`` (string → registry)
    parametrize shapes. Returns ``(None, ())`` for tests not parametrized
    over a case or for cases without an assigned ``vv_level``.
    """
    case = _resolve_case(item)
    if case is None:
        return None, ()
    level = getattr(case, "vv_level", None)
    labels = tuple(getattr(case, "equation_labels", ()) or ())
    return level, labels


def _collect_str_marker_args(item: pytest.Item, marker_name: str) -> tuple[str, ...]:
    """Flatten string args from every ``@pytest.mark.<name>(*args)`` on the item."""
    out: list[str] = []
    for m in item.iter_markers(name=marker_name):
        for arg in m.args:
            if isinstance(arg, str):
                out.append(arg)
            elif isinstance(arg, (list, tuple)):
                out.extend(str(x) for x in arg)
    return tuple(dict.fromkeys(out))  # dedupe preserving order


def _apply_level_marker(item: pytest.Item, level: str) -> None:
    # Physics levels are stored uppercase (``L0``..``L3``) but the
    # pytest marker name is lowercase (``l0``..``l3``). Foundation is
    # already lowercase in both forms.
    marker = "foundation" if level == "foundation" else level.lower()
    item.add_marker(getattr(pytest.mark, marker))


#: The ``ORPHEUS_RUN_WITHDRAWN`` lift, parsed ONCE per session in
#: :func:`pytest_configure` (a mistyped value is a usage error there, never
#: an internal error inside collection).
_LIFTED = pytest.StashKey[frozenset[int] | AllWithdrawals]()


def pytest_configure(config: pytest.Config) -> None:
    """Parse ``ORPHEUS_RUN_WITHDRAWN`` once; a mistyped lift is a ``UsageError``."""
    try:
        config.stash[_LIFTED] = lifted_withdrawals()
    except ValueError as exc:
        raise pytest.UsageError(
            f"{exc} (unset {RUN_WITHDRAWN_VARIABLE} to run with every withdrawal in force)"
        ) from exc


def _withdrawal_of(item: pytest.Item) -> Withdrawal | None:
    """Parse the item's ``@pytest.mark.withdrawn(reason, issue=N)``, or ``None``.

    The marker is resolved by pytest (function, then class, then module
    ``pytestmark``). A malformed marker, or one placed on a single
    ``pytest.param``, is a ``UsageError`` naming the test: a withdrawal
    is a property of a test FUNCTION or class (0 of the parametrised
    functions in the #506 set are mixed; ``.claude/plans/reference_p0_spec.md``
    §1.2), so a per-case withdrawal is a misplacement, never a need.
    """
    mark = item.get_closest_marker("withdrawn")
    if mark is None:
        return None
    callspec = getattr(item, "callspec", None)
    if callspec is not None and any(m.name == "withdrawn" for m in callspec.marks):
        raise pytest.UsageError(
            f"{item.nodeid}: @pytest.mark.withdrawn is placed on a pytest.param; "
            "withdraw the test function or class instead"
        )
    try:
        return Withdrawal.from_mark(mark, where=item.nodeid)
    except ValueError as exc:
        raise pytest.UsageError(str(exc)) from exc


def _apply_withdrawal(item: pytest.Item, withdrawal: Withdrawal | None) -> None:
    """Skip a withdrawn test, with its reason and issue, unless its issue is lifted.

    ``-rs`` prints the reason; the default summary counts the skip, so a
    withdrawal is never a silent deselection. The lift is the session's,
    parsed in :func:`pytest_configure`; the generator's lock
    (:func:`orpheus.derivations.common.withdrawal.withdrawn_generator`)
    reads the same variable again at call time.
    """
    if withdrawal is not None and withdrawal.issue not in item.config.stash[_LIFTED]:
        item.add_marker(
            pytest.mark.skip(
                reason=f"withdrawn (#{withdrawal.issue}): {withdrawal.reason}"
            )
        )


@pytest.hookimpl(wrapper=True)
def pytest_runtest_makereport(item: pytest.Item, call: pytest.CallInfo[None]):
    """No xfail absorbs a withdrawn generator's refusal, in any phase.

    The lock's refusal (:class:`~orpheus.derivations.common.withdrawal.GeneratorWithdrawn`)
    is a policy, never the expected failure an ``xfail`` documents, so a
    report that would read ``xfailed`` because the test body, or a fixture
    in its setup or teardown, raised it is turned into a failure carrying the
    lock's own message (which names the issue and the ``withdrawn`` marker to
    add). Without this, an unmarked xfail test that reaches a locked
    generator reads ``x`` and the placement is unenforced there (``[M]``
    2026-09-25: 12 of the 294 #506 cases were such xfails).

    **ELEGANCE-DEBT[guard] #506** — a report rewrite stands where a withdrawn
    generator can still be CALLED; it retires at phase P4 of
    ``.claude/plans/reference_cache.md``, when the ``ReferenceCertificate``
    carries the ``Withdrawn`` state and no withdrawn generator runs.
    """
    report = yield
    if (
        call.excinfo is not None
        and call.excinfo.errisinstance(GeneratorWithdrawn)
        and hasattr(report, "wasxfail")
    ):
        del report.wasxfail
        report.outcome = "failed"
        report.longrepr = (
            f"{item.nodeid}: an xfail absorbed a withdrawn generator's refusal; "
            f"the refusal is not an expected failure. {call.excinfo.value}"
        )
    return report


def pytest_collection_modifyitems(
    config: pytest.Config, items: list[pytest.Item]
) -> None:
    """Resolve V&V level per item and populate ``TEST_REGISTRY``.

    Precedence (most specific wins):

    1. Explicit marker already on the test (``@pytest.mark.lN``,
       file-level ``pytestmark``, or ``@verify.lN(...)`` which stamps
       an explicit marker).
    2. Class name matching ``TestL<N>Foo``.
    3. Function name matching ``test_l<N>_*``.
    4. ``VerificationCase.vv_level`` inherited through a parametrized
       ``case`` argument (requires PR-2 to populate case metadata; in
       PR-1 this branch is inert until cases are tagged).
    5. Unmarked — recorded in the registry with ``level=None`` so the
       audit tool can surface it.

    A test carrying ``@pytest.mark.withdrawn(reason, issue=N)`` is skipped
    unless ``ORPHEUS_RUN_WITHDRAWN`` names ``N`` (:func:`_apply_withdrawal`),
    and its :class:`~orpheus.derivations.common.withdrawal.Withdrawal` is
    recorded on its registry entry either way.
    """
    registry.clear()
    for item in items:
        level, had_explicit = _existing_level(item)
        source: str = "explicit" if had_explicit else "unmarked"

        if level is None:
            cls_level = _level_from_class_name(item)
            if cls_level is not None:
                level, source = cls_level, "class-name"

        if level is None:
            fn_level = _level_from_func_name(item)
            if fn_level is not None:
                level, source = fn_level, "func-name"

        inherited_equations: tuple[str, ...] = ()
        if level is None:
            case_level, case_labels = _level_from_case(item)
            if case_level is not None:
                level, source = case_level, "case"
                inherited_equations = case_labels

        if level is not None and not had_explicit:
            _apply_level_marker(item, level)

        explicit_equations = _collect_str_marker_args(item, "verifies")
        catches = _collect_str_marker_args(item, "catches")
        case_names = _case_names_from_parametrize(item)
        slow = any(m.name == "slow" for m in item.iter_markers())

        equations = tuple(
            dict.fromkeys(explicit_equations + inherited_equations)
        )

        withdrawal = _withdrawal_of(item)
        _apply_withdrawal(item, withdrawal)

        nodeid = item.nodeid
        file_path = nodeid.split("::", 1)[0]

        registry.record(
            TestMetadata(
                nodeid=nodeid,
                file=file_path,
                level=level,  # type: ignore[arg-type]
                level_source=source,  # type: ignore[arg-type]
                equations=equations,
                catches=catches,
                case_names=case_names,
                slow=slow,
                withdrawn=withdrawal,
            )
        )
