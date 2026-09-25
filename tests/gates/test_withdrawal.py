r"""The withdrawal mechanism: the value, the marker, the skip, the placement, the lock.

A withdrawn reference generator (``.claude/plans/reference_p0_spec.md`` §2;
#506, the Peierls Nyström solver half) is declared twice and the two
declarations must agree:

* statically, by ``@pytest.mark.withdrawn(reason, issue=N)`` on every test
  that consumes it, which ``tests/conftest.py`` parses with
  :meth:`~orpheus.derivations.common.withdrawal.Withdrawal.from_mark`, turns
  into a skip unless ``ORPHEUS_RUN_WITHDRAWN`` names ``N``, and records on
  the registry (so the V&V matrix and the error catalogue count the test as
  neither verifying nor catching);
* at run time, by :func:`~orpheus.derivations.common.withdrawal.withdrawn_generator`
  on every withdrawn symbol, which refuses the call with
  :class:`~orpheus.derivations.common.withdrawal.GeneratorWithdrawn` unless
  the issue is lifted.

The seven gates (the specification's M1 to M6, and M7 from the review):

* **M1** the :class:`Withdrawal` constructor law, one leg per clause;
* **M2** :meth:`Withdrawal.from_mark` on real ``pytest`` marks, and the
  hook's collection refusals (a malformed marker, a ``pytest.param``
  placement) in a child ``pytest``;
* **M3** the skip and the opt-in end to end, in a child ``pytest`` over
  ``TestMGInputValidation`` (4 cases): unset, ``506`` and ``999``; a
  mistyped opt-in is a usage error;
* **M4** the placement census: the tests withdrawn under #506 in the whole
  collected tree are exactly the committed list
  ``tests/gates/withdrawal_506_placement.txt``; each carries the one
  test-side mark :data:`tests._harness.withdrawals.PEIERLS_NYSTROM_WITHDRAWN`,
  minted from the package's :data:`PEIERLS_NYSTROM_WITHDRAWAL`, and no test
  file spells the marker by hand;
* **M5** the placement theorem: a full run of the 26 files that hold the
  withdrawn tests, plus ``test_peierls_assembly_drivers.py`` (the control
  file of the same package, which holds none), has no failure (a kept test
  that reached a locked generator would be red) and skips exactly the
  listed ids;
* **M6** the lock's own laws, in process and in a child interpreter; the
  refusal pickles; a mistyped lift refuses the call; and the production
  census, tree-wide over ``orpheus/`` by AST and at run time over the eight
  modules of ``peierls_nystrom``: the locked symbols are exactly the 34
  of #506;
* **M7** no ``xfail`` absorbs the lock, in the call phase or in a
  function- or module-scoped fixture's setup (``tests/conftest.py``,
  ``pytest_runtest_makereport``).

All seven are ``foundation`` (software invariants; no ``verifies``). M5
runs the 27 files again (about 35 s), so it is also ``slow``; every kept
test in those files enforces the same theorem in its own run, and M5 adds
the end-to-end count.
"""

from __future__ import annotations

import inspect
import os
import subprocess
import sys
import textwrap
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any

import pytest

from orpheus.derivations.common.withdrawal import (
    ALL_WITHDRAWALS,
    RUN_WITHDRAWN_VARIABLE,
    GeneratorWithdrawn,
    Withdrawal,
    lifted_withdrawals,
    withdrawal_of,
    withdrawn_generator,
)
from orpheus.derivations.continuous.peierls_nystrom import PEIERLS_NYSTROM_WITHDRAWAL
from tests._harness.withdrawals import PEIERLS_NYSTROM_WITHDRAWN

pytestmark = pytest.mark.foundation

REPO_ROOT = Path(__file__).resolve().parents[2]
PLACEMENT = Path(__file__).with_name("withdrawal_506_placement.txt")

#: The symbols #506 withdraws, by module of ``peierls_nystrom`` (spec §1.1).
#: ``BoundaryClosureOperator`` is a frozen dataclass whose ``__init__`` is
#: generated, so its construction is locked at ``__post_init__``.
WITHDRAWN_506 = {
    "geometry": {
        "K_vol_element_adaptive", "build_volume_kernel_adaptive",
        "build_volume_kernel", "build_white_bc_correction",
        "build_white_bc_correction_rank_n", "build_closure_operator",
        "_build_closure_operator_rank2_white", "_build_closure_operator_rank_n_white",
        "_build_slab_per_face_specular_PG", "_build_sphere_specular_mode_PG",
        "_build_cylinder_specular_mode_PG", "_build_white_rank1_mark_op",
        "_build_white_f4_op", "_build_white_hebert_op", "_build_specular_op",
        "_build_specular_multibounce_op", "_build_full_K_per_group",
        "solve_peierls_mg", "solve_peierls_1g",
        "BoundaryClosureOperator.__post_init__",
    },
    "slab": {
        "_build_kernel_matrix", "_build_system_matrices",
        "solve_peierls_eigenvalue", "_build_peierls_slab_case",
    },
    "cylinder": {"_build_peierls_cylinder_case", "_build_peierls_cylinder_hollow_f4_case"},
    "sphere": {"_build_peierls_sphere_case", "_build_peierls_sphere_hollow_f4_case"},
    "cases": {
        "build_two_surface_case", "_build_peierls_slab_case_via_unified",
        "build_one_surface_compact_case", "_build", "_class_a_cases",
        "continuous_cases",
    },
}


def _placement() -> set[str]:
    return {
        line.strip() for line in PLACEMENT.read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.startswith("#")
    }


def _env(lift: str | None = None) -> dict[str, str]:
    """The current environment with ``ORPHEUS_RUN_WITHDRAWN`` set to ``lift``,
    or removed when ``lift`` is ``None``."""
    env = {k: v for k, v in os.environ.items() if k != RUN_WITHDRAWN_VARIABLE}
    if lift is not None:
        env[RUN_WITHDRAWN_VARIABLE] = lift
    return env


# ---------------------------------------------------------------------------
# M1 — the constructor law
# ---------------------------------------------------------------------------


def test_m1_a_withdrawal_is_a_reason_and_a_positive_issue():
    """M1, the accepting leg: the #506 spelling constructs."""
    w = Withdrawal("r", issue=506)
    assert (w.reason, w.issue) == ("r", 506)


@pytest.mark.parametrize(
    "kwargs, match",
    [
        ({"reason": "", "issue": 506}, "reason must be a non-empty string"),
        ({"reason": "   ", "issue": 506}, "reason must be a non-empty string"),
        ({"reason": "r", "issue": "506"}, "issue must be an int"),
        ({"reason": "r", "issue": True}, "issue must be an int"),
        ({"reason": "r", "issue": 0}, "issue must be a positive issue number"),
        ({"reason": "r", "issue": -3}, "issue must be a positive issue number"),
    ],
    ids=["empty-reason", "blank-reason", "string-issue", "bool-issue", "zero-issue", "negative-issue"],
)
def test_m1_the_constructor_refuses_each_malformed_field(kwargs, match):
    """M1, one refusal leg per clause of ``Withdrawal.__post_init__``."""
    with pytest.raises(ValueError, match=match):
        Withdrawal(**kwargs)


def test_m1_an_issue_is_required():
    """M1: there is no default issue; a withdrawal names the ruling's record."""
    reason_only: dict[str, Any] = {"reason": "r"}
    with pytest.raises(TypeError, match="issue"):
        Withdrawal(**reason_only)


# ---------------------------------------------------------------------------
# M2 — the marker parser and the hook's collection refusals
# ---------------------------------------------------------------------------


def test_m2_from_mark_parses_a_real_pytest_mark():
    """M2, accepting leg: a real ``Mark`` parses to the value."""
    mark = pytest.mark.withdrawn("r", issue=506).mark
    assert Withdrawal.from_mark(mark, where="t::x") == Withdrawal("r", issue=506)


@pytest.mark.parametrize(
    "decorator, match",
    [
        (pytest.mark.withdrawn("r"), r"t::x: a withdrawn marker must name its issue"),
        (pytest.mark.withdrawn(issue=506), r"t::x: .*exactly one positional"),
        (pytest.mark.withdrawn("r", "s", issue=506), r"t::x: .*exactly one positional"),
        (pytest.mark.withdrawn("r", issue=506, why="x"), r"t::x: .*only issue="),
        (pytest.mark.withdrawn("r", issue="506"), r"t::x: .*issue must be an int"),
    ],
    ids=["no-issue", "no-reason", "two-reasons", "extra-kwarg", "string-issue"],
)
def test_m2_from_mark_refuses_a_malformed_marker_naming_the_test(decorator, match):
    """M2, refusal legs: each names the carrying test's node id."""
    with pytest.raises(ValueError, match=match):
        Withdrawal.from_mark(decorator.mark, where="t::x")


_CHILD_TESTS = {
    "test_param.py": """
        import pytest
        @pytest.mark.parametrize(
            "x", [1, pytest.param(2, marks=pytest.mark.withdrawn("r", issue=506))])
        def test_p(x):
            pass
    """,
    "test_noissue.py": """
        import pytest
        @pytest.mark.withdrawn("r")
        def test_n():
            pass
    """,
}


def _child_pytest(path: Path, env: dict[str, str], *extra: str) -> subprocess.CompletedProcess[str]:
    """A child ``python -O -m pytest`` on ``path`` with the repository's conftest
    loaded as a plugin (``path`` lives outside ``tests/``)."""
    return subprocess.run(
        [sys.executable, "-O", "-m", "pytest", "-p", "tests.conftest",
         "-c", str(REPO_ROOT / "pyproject.toml"), "-p", "no:cacheprovider",
         "--color=no", "-q", "-rs", *extra, str(path)],
        cwd=REPO_ROOT, env=env, capture_output=True, text=True, timeout=300,
    )


@pytest.mark.parametrize(
    "name, message",
    [
        ("test_param.py", "is placed on a pytest.param"),
        ("test_noissue.py", "a withdrawn marker must name its issue"),
    ],
)
def test_m2_the_hook_refuses_a_malformed_or_per_case_marker(tmp_path, name, message):
    """M2, the hook: a malformed marker, or one on a single ``pytest.param``,
    is a collection-time usage error (exit 4) naming the test."""
    path = tmp_path / name
    path.write_text(textwrap.dedent(_CHILD_TESTS[name]), encoding="utf-8")
    proc = _child_pytest(path, _env())
    assert proc.returncode == 4, (proc.returncode, proc.stdout[-2000:], proc.stderr[-2000:])
    assert message in proc.stdout + proc.stderr


# ---------------------------------------------------------------------------
# M3 — the skip and the opt-in, end to end
# ---------------------------------------------------------------------------

_M3_TARGET = "tests/gates/derivations/test_peierls_multigroup.py::TestMGInputValidation"


def _summary(proc: subprocess.CompletedProcess[str]) -> str:
    lines = [ln for ln in proc.stdout.splitlines() if ln.strip()]
    return lines[-1] if lines else proc.stderr[-500:]


def _outcomes(proc: subprocess.CompletedProcess[str]) -> dict[str, int]:
    """The summary line's counts, ``{"failed": 1, "errors": 2, ...}``
    (warnings excluded)."""
    import re

    return {
        word: int(n) for n, word in re.findall(r"(\d+) (\w+)", _summary(proc))
        if word not in ("warning", "warnings")
    }


@pytest.mark.parametrize(
    "lift, outcome, reason_printed",
    [
        (None, "4 skipped", True),
        ("506", "4 passed", False),
        ("999", "4 skipped", True),
    ],
    ids=["unset", "lift-506", "lift-999"],
)
def test_m3_a_withdrawn_test_is_skipped_unless_its_issue_is_lifted(lift, outcome, reason_printed):
    """M3: ``TestMGInputValidation`` (4 cases, cheap: each raises on bad input
    before any solve) is skipped with ``withdrawn (#506)`` in ``-rs`` unless
    ``ORPHEUS_RUN_WITHDRAWN`` names 506; lifting another issue lifts nothing."""
    proc = subprocess.run(
        [sys.executable, "-O", "-m", "pytest", "-p", "no:cacheprovider",
         "--color=no", "-q", "-rs", _M3_TARGET],
        cwd=REPO_ROOT, env=_env(lift), capture_output=True, text=True, timeout=300,
    )
    summary = _summary(proc)
    assert summary.startswith(outcome) and "failed" not in summary, summary
    assert ("withdrawn (#506)" in proc.stdout) is reason_printed, proc.stdout[-1500:]


def test_m3_a_mistyped_opt_in_is_a_usage_error():
    """M3: ``ORPHEUS_RUN_WITHDRAWN=5o6`` stops the session at configure time
    with a usage error naming the value (exit 4), never an internal error."""
    proc = subprocess.run(
        [sys.executable, "-O", "-m", "pytest", "-p", "no:cacheprovider",
         "--color=no", "-q", _M3_TARGET],
        cwd=REPO_ROOT, env=_env("5o6"), capture_output=True, text=True, timeout=300,
    )
    assert proc.returncode == 4, (proc.returncode, proc.stdout[-1500:], proc.stderr[-1500:])
    assert "'5o6' is not a positive issue number" in proc.stdout + proc.stderr
    assert "INTERNALERROR" not in proc.stdout + proc.stderr


# ---------------------------------------------------------------------------
# M4 — the placement census
# ---------------------------------------------------------------------------


def test_m4_the_tests_withdrawn_under_506_are_exactly_the_placement():
    """M4: the registry's #506 withdrawals, over the WHOLE collected tree,
    equal the committed placement list (294 node ids), and each parses to the
    package's one :data:`PEIERLS_NYSTROM_WITHDRAWAL`."""
    from tests._harness.audit import audit_payload

    withdrawn = audit_payload()["withdrawn_tests"]
    under_506 = {nid for nid, w in withdrawn.items() if w["issue"] == 506}
    expected = _placement()
    assert len(expected) == 294
    missing, extra = sorted(expected - under_506), sorted(under_506 - expected)
    assert not missing and not extra, (
        f"placement drift — listed but not marked: {missing}; marked but not "
        f"listed: {extra}. Update {PLACEMENT.name} with the marker change."
    )
    values = {(w["reason"], w["issue"]) for nid, w in withdrawn.items() if nid in under_506}
    assert values == {(PEIERLS_NYSTROM_WITHDRAWAL.reason, PEIERLS_NYSTROM_WITHDRAWAL.issue)}, values


def test_m4_the_test_side_mark_is_the_package_value():
    """M4: the one test-side mark parses to the package constant, so the
    reason sentence has one spelling."""
    parsed = Withdrawal.from_mark(PEIERLS_NYSTROM_WITHDRAWN.mark, where="PEIERLS_NYSTROM_WITHDRAWN")
    assert parsed == PEIERLS_NYSTROM_WITHDRAWAL


def _hand_spelled_withdrawn_marks(root: Path) -> list[str]:
    """Every ``pytest.mark.withdrawn(...)`` CALL under ``root`` (AST), as
    ``path:line``."""
    import ast

    sites: list[str] = []
    for path in sorted(root.rglob("*.py")):
        for node in ast.walk(ast.parse(path.read_text(encoding="utf-8"))):
            if (
                isinstance(node, ast.Call)
                and isinstance(node.func, ast.Attribute)
                and node.func.attr == "withdrawn"
                and isinstance(node.func.value, ast.Attribute)
                and node.func.value.attr == "mark"
            ):
                sites.append(f"{path.relative_to(REPO_ROOT)}:{node.lineno}")
    return sites


def test_m4_no_test_spells_the_marker_by_hand():
    """M4: outside the harness module that mints the mark (and this file's
    fixtures of malformed markers), no test spells ``pytest.mark.withdrawn``;
    each withdrawn test uses the minted mark. Positive control: the census
    finds the mint itself."""
    sites = _hand_spelled_withdrawn_marks(REPO_ROOT / "tests")
    assert any(s.startswith("tests/_harness/withdrawals.py:") for s in sites), sites
    stray = [
        s for s in sites
        if not s.startswith(("tests/_harness/withdrawals.py:", "tests/gates/test_withdrawal.py:"))
    ]
    assert not stray, f"hand-spelled withdrawn markers (use PEIERLS_NYSTROM_WITHDRAWN): {stray}"


# ---------------------------------------------------------------------------
# M5 — the placement theorem, end to end
# ---------------------------------------------------------------------------


def _junit_key(nodeid: str) -> tuple[str, str]:
    """A node id as JUnit spells it: ``(dotted module[.Class], name)``."""
    path, *rest = nodeid.split("::")
    return (".".join([path.removesuffix(".py").replace("/", "."), *rest[:-1]]), rest[-1])


@pytest.mark.slow
def test_m5_no_kept_test_reaches_a_locked_generator(tmp_path):
    """M5: the 26 files that hold a #506 withdrawal, plus the package's
    control file ``test_peierls_assembly_drivers.py`` (0 withdrawn), run as
    the canonical invocation, have 0 failures and 0 errors (the lock reds any kept test
    that reaches a withdrawn generator), and the cases skipped with the
    ``withdrawn (#506)`` reason are exactly the listed ids."""
    placement = _placement()
    files = sorted({nid.split("::", 1)[0] for nid in placement}
                   | {"tests/gates/derivations/test_peierls_assembly_drivers.py"})
    assert len(files) == 27
    report = tmp_path / "m5.xml"
    subprocess.run(
        [sys.executable, "-O", "-m", "pytest", "-p", "no:cacheprovider", "--color=no",
         "-q", f"--junitxml={report}", *files],
        cwd=REPO_ROOT, env=_env(), capture_output=True, text=True, timeout=1200,
    )
    red: list[str] = []
    skipped_withdrawn: set[tuple[str, str]] = set()
    for case in ET.parse(report).iter("testcase"):
        key = (case.get("classname", ""), case.get("name", ""))
        for child in case:
            text = (child.get("message") or "") + (child.text or "")
            if child.tag in ("failure", "error"):
                red.append(f"{'::'.join(key)}: {text[:300]}")
            elif child.tag == "skipped" and "withdrawn (#506)" in text:
                skipped_withdrawn.add(key)
    assert not red, "a kept test reached a locked generator, or failed:\n" + "\n".join(red)
    expected = {_junit_key(nid) for nid in placement}
    assert skipped_withdrawn == expected, (
        f"skipped but not listed: {sorted(skipped_withdrawn - expected)}; "
        f"listed but not skipped: {sorted(expected - skipped_withdrawn)}"
    )


_XFAIL_CHILD = """
import pytest
from orpheus.derivations.common.withdrawal import Withdrawal, withdrawn_generator

@withdrawn_generator(Withdrawal("stub", issue=506))
def locked():
    return 1

@pytest.fixture
def locked_function_fixture():
    return locked()

@pytest.fixture(scope="module")
def locked_module_fixture():
    return locked()

@pytest.mark.xfail(reason="an unrelated expected failure")
def test_unmarked_xfail_reaching_a_lock():
    locked()

@pytest.mark.xfail(reason="an unrelated expected failure")
def test_unmarked_xfail_whose_function_fixture_reaches_a_lock(locked_function_fixture):
    pass

@pytest.mark.xfail(reason="an unrelated expected failure")
def test_unmarked_xfail_whose_module_fixture_reaches_a_lock(locked_module_fixture):
    pass

@pytest.mark.xfail(reason="an ordinary expected failure")
def test_ordinary_xfail():
    raise ValueError("expected")
"""


def test_m7_no_xfail_absorbs_the_lock(tmp_path):
    """M7: an UNMARKED ``xfail`` test that reaches a locked generator reads
    red with the lock's message naming #506 and the marker to add
    (``tests/conftest.py``, ``pytest_runtest_makereport``): FAILED when its
    body reaches it, ERROR when a function- or module-scoped fixture does
    (a setup failure). An ordinary xfail beside them still reads xfailed
    (the control leg)."""
    path = tmp_path / "test_xfail_lock.py"
    path.write_text(_XFAIL_CHILD, encoding="utf-8")
    proc = _child_pytest(path, _env(), "-rfEx")
    assert _outcomes(proc) == {"failed": 1, "xfailed": 1, "errors": 2}, proc.stdout[-3000:]
    assert proc.stdout.count("an xfail absorbed a withdrawn generator's refusal") >= 3
    assert "@pytest.mark.withdrawn(" in proc.stdout and "issue=506" in proc.stdout


# ---------------------------------------------------------------------------
# M6 — the lock's own laws
# ---------------------------------------------------------------------------

_STUB_WITHDRAWAL = Withdrawal("stub", issue=506)


@withdrawn_generator(_STUB_WITHDRAWAL)
def _stub_generator(x: float, *, scale: float = 2.0) -> float:
    """A stub generator."""
    return scale * x


def test_m6_the_lock_refuses_unless_lifted_and_reads_the_variable_at_call_time(monkeypatch):
    """M6: refuses when unset, runs when 506 (set AFTER the decoration, at
    import), refuses under another issue, runs under ``all``."""
    monkeypatch.delenv(RUN_WITHDRAWN_VARIABLE, raising=False)
    with pytest.raises(GeneratorWithdrawn, match=r"_stub_generator is withdrawn \(#506\)"):
        _stub_generator(1.0)
    monkeypatch.setenv(RUN_WITHDRAWN_VARIABLE, "506")
    assert _stub_generator(1.5) == 3.0
    monkeypatch.setenv(RUN_WITHDRAWN_VARIABLE, "999")
    with pytest.raises(GeneratorWithdrawn):
        _stub_generator(1.0)
    monkeypatch.setenv(RUN_WITHDRAWN_VARIABLE, "all")
    assert _stub_generator(1.0, scale=3.0) == 3.0


@pytest.mark.parametrize(
    "value, lifted",
    [("", frozenset()), ("506", frozenset({506})), (" 506, 512 ", frozenset({506, 512})),
     ("all", ALL_WITHDRAWALS), ("ALL", ALL_WITHDRAWALS)],
)
def test_m6_the_opt_in_parses(monkeypatch, value, lifted):
    """M6: the opt-in's grammar, comma-separated issues or ``all``."""
    monkeypatch.setenv(RUN_WITHDRAWN_VARIABLE, value)
    assert lifted_withdrawals() == lifted


@pytest.mark.parametrize("value", ["abc", "506,", "0", "#506", "506;512"])
def test_m6_a_mistyped_opt_in_is_refused_not_read_as_nothing(monkeypatch, value):
    """M6: a mistyped lift raises, naming the variable (a silent "lift
    nothing" would leave the caller believing the generator ran)."""
    monkeypatch.setenv(RUN_WITHDRAWN_VARIABLE, value)
    with pytest.raises(ValueError, match=RUN_WITHDRAWN_VARIABLE):
        lifted_withdrawals()


def test_m6_the_lock_is_not_catchable_as_an_exception(monkeypatch):
    """M6: ``GeneratorWithdrawn`` is a ``BaseException``, so an ``except
    Exception`` fallback around a generator call cannot absorb it."""
    monkeypatch.delenv(RUN_WITHDRAWN_VARIABLE, raising=False)
    assert not issubclass(GeneratorWithdrawn, Exception)
    with pytest.raises(GeneratorWithdrawn):
        try:
            _stub_generator(1.0)
        except Exception:  # noqa: BLE001 — the fallback the base class defeats
            pytest.fail("an `except Exception` absorbed the lock")


def test_m6_the_lock_preserves_the_generators_identity():
    """M6: name, docstring, signature and ``__wrapped__`` survive, so the lazy
    registry, ``inspect`` and attribute ``monkeypatch`` see the generator."""
    assert _stub_generator.__name__ == "_stub_generator"
    assert _stub_generator.__doc__ == "A stub generator."
    assert str(inspect.signature(_stub_generator)) == "(x: 'float', *, scale: 'float' = 2.0) -> 'float'"
    assert inspect.unwrap(_stub_generator)(1.0) == 2.0
    assert withdrawal_of(_stub_generator) is _STUB_WITHDRAWAL
    assert withdrawal_of(inspect.unwrap(_stub_generator)) is None


def test_m6_a_mistyped_lift_refuses_the_call_naming_the_value(monkeypatch):
    """M6: inside the lock a lift that does not parse is a refusal (never a
    bare ``ValueError``, never a run), and it names the bad value."""
    monkeypatch.setenv(RUN_WITHDRAWN_VARIABLE, "5o6")
    with pytest.raises(GeneratorWithdrawn, match="'5o6' is not a positive issue number") as info:
        _stub_generator(1.0)
    assert info.value.invalid_lift is not None


def test_m6_the_refusal_pickles(monkeypatch):
    """M6: a refusal raised in a worker process reaches the parent as itself."""
    import pickle

    monkeypatch.delenv(RUN_WITHDRAWN_VARIABLE, raising=False)
    with pytest.raises(GeneratorWithdrawn) as info:
        _stub_generator(1.0)
    back = pickle.loads(pickle.dumps(info.value))
    assert type(back) is GeneratorWithdrawn
    assert (back.withdrawal, back.generator, back.invalid_lift, str(back)) == (
        info.value.withdrawal, info.value.generator, None, str(info.value)
    )


_CHILD_LOCK = """
import sys
from orpheus.derivations.common.withdrawal import GeneratorWithdrawn
from orpheus.derivations.continuous.peierls_nystrom.cases import continuous_cases
try:
    continuous_cases()
except GeneratorWithdrawn as exc:
    print("REFUSED", exc.withdrawal.issue)
    sys.exit(0)
print("RAN")
"""


@pytest.mark.parametrize("lift, expected", [(None, "REFUSED 506"), ("999", "REFUSED 506")])
def test_m6_a_child_interpreter_is_locked_through_the_environment(lift, expected):
    """M6, the subprocess leg: a child ``python`` that calls a real withdrawn
    generator is refused (no pytest, no marker: the lock alone)."""
    proc = subprocess.run(
        [sys.executable, "-O", "-c", _CHILD_LOCK],
        cwd=REPO_ROOT, env=_env(lift), capture_output=True, text=True, timeout=300,
    )
    assert proc.stdout.strip() == expected, (proc.stdout, proc.stderr[-1500:])


_CHILD_STUB = """
from orpheus.derivations.common.withdrawal import Withdrawal, withdrawn_generator
@withdrawn_generator(Withdrawal("stub", issue=506))
def g():
    return "RAN"
try:
    print(g())
except BaseException as exc:
    print(type(exc).__name__)
"""


@pytest.mark.parametrize(
    "lift, expected",
    [(None, "GeneratorWithdrawn"), ("506", "RAN"), ("all", "RAN")],
)
def test_m6_a_child_interpreter_inherits_the_lift(lift, expected):
    """M6, the subprocess leg: the lift travels in the environment."""
    proc = subprocess.run(
        [sys.executable, "-O", "-c", _CHILD_STUB],
        cwd=REPO_ROOT, env=_env(lift), capture_output=True, text=True, timeout=120,
    )
    assert proc.stdout.strip() == expected, (proc.stdout, proc.stderr[-1500:])


def _locked_symbols_by_ast(root: Path) -> set[str]:
    """``module.Qualname`` of every def decorated ``@withdrawn_generator(...)``
    under ``root`` (AST: the whole tree, nothing imported)."""
    import ast

    found: set[str] = set()
    for path in sorted(root.rglob("*.py")):
        module = ".".join(path.relative_to(REPO_ROOT).with_suffix("").parts)

        def visit(body: list[ast.stmt], prefix: str) -> None:
            for node in body:
                if isinstance(node, ast.ClassDef):
                    visit(node.body, f"{prefix}{node.name}.")
                elif isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                    for deco in node.decorator_list:
                        call = deco.func if isinstance(deco, ast.Call) else deco
                        if isinstance(call, ast.Name) and call.id == "withdrawn_generator":
                            found.add(f"{module}.{prefix}{node.name}")

        visit(ast.parse(path.read_text(encoding="utf-8")).body, "")
    return found


_PACKAGE = "orpheus.derivations.continuous.peierls_nystrom"


def test_m6_the_locked_symbols_of_the_tree_are_exactly_the_506_set():
    """M6, the production census over the whole of ``orpheus/`` by AST: the
    ``@withdrawn_generator`` sites are exactly the 34 of #506 (the escape
    primitives and their angular-assembly drivers stay in service).
    Positive control: the census finds ``solve_peierls_1g``."""
    expected = {f"{_PACKAGE}.{m}.{q}" for m, names in WITHDRAWN_506.items() for q in names}
    found = _locked_symbols_by_ast(REPO_ROOT / "orpheus")
    assert f"{_PACKAGE}.geometry.solve_peierls_1g" in found
    assert found == expected, (sorted(found - expected), sorted(expected - found))
    assert len(found) == 34


def test_m6_the_locked_symbols_of_the_package_are_exactly_the_506_set():
    """M6, the run-time census over the eight modules of ``peierls_nystrom``
    (``geometry``, ``slab``, ``cylinder``, ``sphere``, ``cases``, ``naming``,
    ``reference``, ``ps1982_reference``): every function and method locked
    at import is in the #506 set, carries the package's one constant, and
    every member of the set is locked."""
    import importlib

    found: dict[str, set[str]] = {}
    for module_name in ("geometry", "slab", "cylinder", "sphere", "cases", "naming",
                        "reference", "ps1982_reference"):
        module = importlib.import_module(f"{_PACKAGE}.{module_name}")
        locked: set[str] = set()
        for value in vars(module).values():
            if getattr(value, "__module__", None) != module.__name__:
                continue
            withdrawal = withdrawal_of(value)
            if withdrawal is not None:
                assert withdrawal is PEIERLS_NYSTROM_WITHDRAWAL, value
                # By the function's own name: ``cases = continuous_cases`` is
                # an alias binding, one symbol.
                locked.add(value.__qualname__)
            if inspect.isclass(value):
                for attr, member in vars(value).items():
                    member_withdrawal = withdrawal_of(member)
                    if member_withdrawal is not None:
                        assert member_withdrawal is PEIERLS_NYSTROM_WITHDRAWAL
                        locked.add(f"{value.__name__}.{attr}")
        if locked:
            found[module_name] = locked
    assert found == WITHDRAWN_506
    assert sum(len(v) for v in found.values()) == 34
