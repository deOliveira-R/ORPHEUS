"""Mutation battery for the #340/#342 convergence contract.

vv anti-#17: a gate suite that only ever sees the FIXED code proves nothing.
Each mutation below re-introduces one of the defects the contract forbids;
the suite must redden.  MA0 is the positive control and runs FIRST — an
all-blind verdict is a BROKEN HARNESS until the control says otherwise.

In-process monkeypatch only — never a git-level revert (the tree carries
uncommitted-by-policy state).

⛔ REPAIRED 2026-08-10 (#340 N6).  This file was TRACKED and could not
IMPORT at HEAD: line 24 bound ``sv._claims_convergence``, a predicate the
campaign retired on 2026-08-09.  A campaign that retires symbols breaks its
own instruments by MODULE-SCOPE BINDING, silently, because nobody runs the
harness between carves.  Two further staleness modes, both anti-#18
(over-powered: their reds are CRASHES, not property reds) —

* ``ma1`` built ``PowerIterationOutcome(converged=True)``, a keyword removed
  by N2b-ii and whose absence has its own gate.  RETIRED: the defect is now
  unspellable, so a mutation cannot express it.
* ``ma4`` / ``ma5`` declared ``(history, *, where, budget_name, budget,
  tol)``; N6 retired the last three, so the replacements raised ``TypeError``
  at every call site.  RE-SIGNED.

⭐ M1/M2/M3/M8/MN/MR are built by SOURCE TRANSFORMATION of the live
function rather than by hand-copying it: a hand copy is a twin path that
drifts from production, and a `str.replace` whose target is absent RAISES
(the instrument asserts its own installation, rather than printing a banner
nobody reads).
"""
from __future__ import annotations

import contextlib
import inspect
import io
import re
import textwrap
from unittest import mock

import pytest as _pytest

import orpheus.numerics.convergence as cv
import orpheus.numerics.eigenvalue as ev
import orpheus.numerics.iteration as it
import orpheus.sn.solver as sv

SUITE = [
    "tests/sn/solve/test_convergence_contract.py",
    "tests/numerics/test_iteration_record.py",
    "tests/numerics/test_power_iteration_record.py",
]
REAL_WARN = sv._warn_if_unconverged
REAL_SI_INIT = it.SourceIteration.__init__
REAL_KRY_INIT = it.KrylovAcceleration.__init__
REAL_PI = ev.power_iteration
REAL_FIRST_FAILURE = cv.IterationRecord.first_failure
REAL_POST_INIT = cv.IterationRecord.__post_init__


def run(paths=SUITE):
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(buf):
        _pytest.main(
            ["-q", "-rf", "--color=no", "-p", "no:cacheprovider", *paths]
        )
    out = buf.getvalue()
    summary = ""
    for line in out.splitlines():
        if re.search(r"\d+ (passed|failed|error)", line):
            summary = line.strip()
    failed = int(m.group(1)) if (m := re.search(r"(\d+) failed", out)) else 0
    return failed, sorted({i.split("::")[-1] for i in re.findall(r"FAILED (\S+)", out)}), summary


# ── source-transformed mutants of the live warning helper ───────────────
def _warn_mutant(*replacements: tuple[str, str]):
    """A copy of the LIVE ``_warn_if_unconverged`` with N lines swapped.

    Raises if any target is absent — the mutation cannot silently no-op,
    which is the failure mode that manufactures a false "gate is blind".
    """
    source = textwrap.dedent(inspect.getsource(REAL_WARN))
    for old, new in replacements:
        if old not in source:
            raise RuntimeError(f"mutation target absent from production: {old!r}")
        source = source.replace(old, new, 1)
    namespace = dict(sv.__dict__)
    exec(compile(source, "<n6-mutant>", "exec"), namespace)  # noqa: S102
    return namespace["_warn_if_unconverged"]


_NAVIGATE = "failing = history.record.first_failure or history.record"
_CRITERION = "criterion = failing.binding_criterion"
_KNOB = "budget_name = failing.budget_name"
_BUDGET = "budget = failing.budget"
_GUARD = "if history.converged:"


def report(name, expect, got, detail):
    print(f"  {name:<26} expect={expect:<28} got={got}")
    print(f"     {detail}")


def case(name, expect, patches, note=""):
    """Install ``patches`` (a list of context managers), run, report."""
    with contextlib.ExitStack() as stack:
        for patch in patches:
            stack.enter_context(patch)
        failed, names, summary = run()
    report(name, expect, f"{failed} failed", f"{summary}")
    if names:
        print(f"     reds: {', '.join(names[:12])}"
              + (f" (+{len(names) - 12} more)" if len(names) > 12 else ""))
    if note:
        print(f"     {note}")
    return failed, names


def main() -> None:
    print("=" * 78)
    print("#340 N6 mutation battery — control first (vv anti-#17)")
    print("=" * 78)

    failed, names, summary = run()
    report("M0 CONTROL unmutated", "0 failed", f"{failed} failed", summary)
    if failed:
        print("  !!!! CONTROL FAILED — the harness or the suite is broken; "
              "no negative below is trustworthy.")
        return

    print("\n-- the positive control: is the emission path reachable at all? --")
    case("MC warning suppressed", ">=5 (audibility legs)",
         [mock.patch.object(sv, "_warn_if_unconverged", lambda *a, **k: None)])

    print("\n-- the N6 claim: every fact comes from the FAILING level --")
    case("M1 no navigation", "keystone + nested only",
         [mock.patch.object(sv, "_warn_if_unconverged",
                            _warn_mutant((_NAVIGATE, "failing = history.record")))])
    case("M2 criterion left on top", "the tol/rate/advice legs",
         [mock.patch.object(sv, "_warn_if_unconverged", _warn_mutant(
             (_CRITERION, "criterion = history.record.binding_criterion")))])
    case("M3 knob+budget left on top", "the level/knob/budget legs",
         [mock.patch.object(sv, "_warn_if_unconverged", _warn_mutant(
             (_KNOB, "budget_name = history.record.budget_name"),
             (_BUDGET, "budget = history.record.budget")))])
    case("M7 first_failure SELF-first", "keystone + nested + record gates",
         [mock.patch.object(
             cv.IterationRecord, "first_failure",
             property(lambda self: None if self.converged else self))])

    print("\n-- the producer stamps --")
    def _si_forget(self, *a, budget_name="max_iter", **kw):
        return REAL_SI_INIT(self, *a, budget_name="max_iter", **kw)

    def _si_swap(self, *a, budget_name="max_iter", **kw):
        return REAL_SI_INIT(self, *a, budget_name="max_outer", **kw)

    def _kry_forget(self, *a, budget_name="max_iter", **kw):
        return REAL_KRY_INIT(self, *a, budget_name="max_iter", **kw)

    def _pi_forget(solver, *a, budget_name="max_iter", **kw):
        return REAL_PI(solver, *a, budget_name="max_iter", **kw)

    def _pi_swap(solver, *a, budget_name="max_iter", **kw):
        return REAL_PI(solver, *a, budget_name="max_inner", **kw)

    case("M4 SI+Krylov forget knob", "leaf msg legs + knob sweep",
         [mock.patch.object(it.SourceIteration, "__init__", _si_forget),
          mock.patch.object(it.KrylovAcceleration, "__init__", _kry_forget)])
    case("M5 SWAP the two knobs", "knob sweep exact-legs",
         [mock.patch.object(it.SourceIteration, "__init__", _si_swap),
          mock.patch.object(sv, "power_iteration", _pi_swap),
          mock.patch.object(ev, "power_iteration", _pi_swap)],
         note="the membership leg alone is BLIND here — solve_sn has both knobs")
    case("M6 power_iteration forgets", "eigenvalue rows ONLY (not leaf)",
         [mock.patch.object(sv, "power_iteration", _pi_forget),
          mock.patch.object(ev, "power_iteration", _pi_forget)])

    print("\n-- the guard, the field's invariant, and the retirement --")
    case("M8 guard -> fully_converged", "the commit-1 scope row + xpass",
         [mock.patch.object(sv, "_warn_if_unconverged",
                            _warn_mutant((_GUARD, "if history.fully_converged:")))])
    case("M9 drop the empty-knob guard", "the field's negative leg",
         [mock.patch.object(
             cv.IterationRecord, "__post_init__",
             lambda self: REAL_POST_INIT(
                 object.__setattr__(self, "budget_name",
                                    self.budget_name or "x") or self))],
         note="substitutes a placeholder so the refusal never fires")
    case("MR re-add a defaulted tol=", "the signature row ONLY",
         [mock.patch.object(sv, "_warn_if_unconverged", _warn_mutant(
             ("    where: str,\n", "    where: str,\n    tol: float = 1e-8,\n")))])

    print("\n-- labelled NON-CATCHERS (not coverage; recorded as such) --")
    case("MN drop the None fallback", "0 — provably unreachable",
         [mock.patch.object(sv, "_warn_if_unconverged",
                            _warn_mutant((" or history.record", "")))],
         note="first_failure is None IFF fully_converged, and the guard "
              "returns on converged; a 0 here is the CORRECT result")
    case("MX first_failure := None", "many, ALL via AttributeError",
         [mock.patch.object(sv, "_warn_if_unconverged",
                            _warn_mutant((_NAVIGATE, "failing = None")))],
         note="anti-#18 NEGATIVE EXAMPLE — breaks a type contract, not the "
              "sourcing property. Its reds prove nothing; do NOT count them")


if __name__ == "__main__":
    main()
