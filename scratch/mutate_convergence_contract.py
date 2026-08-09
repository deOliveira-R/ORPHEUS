"""Mutation battery for the #340/#342 convergence contract.

vv anti-#17: a gate suite that only ever sees the FIXED code proves nothing.
Each mutation below re-introduces one of the defects the contract forbids;
the suite must redden.  MA0 is the positive control and runs FIRST.

In-process monkeypatch only — never a git-level revert (the tree carries
uncommitted-by-policy state).
"""
from __future__ import annotations

import contextlib
import io
import re
from unittest import mock

import pytest as _pytest

import orpheus.numerics.eigenvalue as ev
import orpheus.sn.solver as sv

SUITE = ["tests/sn/solve/test_convergence_contract.py"]
REAL_PI = ev.power_iteration
REAL_CLAIMS = sv._claims_convergence
REAL_WARN = sv._warn_if_unconverged


def run(paths):
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(buf):
        _pytest.main(["-q", "-rf", "--color=no", "-p", "no:cacheprovider", *paths])
    out = buf.getvalue()
    m = re.search(r"(\d+) failed", out)
    return (int(m.group(1)) if m else 0), re.findall(r"FAILED (\S+)", out)


def report(name, expect, ok, detail):
    print(f"[{'PASS' if ok else '!!!! FAIL'}] {name}: {expect} — {detail}")
    return ok


# ── MA0 POSITIVE CONTROL: nothing mutated; the suite must be GREEN ──────
def ma0():
    fails, ids = run(SUITE)
    ok = fails == 0
    return report("MA0 control (unmutated)", "0 failures", ok, f"fails={fails} {ids}")


# ── MA1: #342 verbatim — the eigenvalue flag hardcoded True ────────────
def ma1():
    def mutated(solver, max_iter=500):
        out = REAL_PI(solver, max_iter=max_iter)
        return ev.PowerIterationOutcome(
            keff=out.keff, keff_history=out.keff_history,
            flux_distribution=out.flux_distribution,
            converged=True,                       # the lie
        )
    with mock.patch.object(sv, "power_iteration", mutated):
        fails, ids = run(SUITE)
    names = {i.split("::")[-1] for i in ids}
    ok = any("starved_budget" in n for n in names)
    return report("MA1 hardcoded converged=True", "the STARVED eigenvalue leg reds",
                  ok, f"fails={fails}: {sorted(names)}")


# ── MA2: the fixed-source predicate always claims convergence ──────────
def ma2():
    with mock.patch.object(sv, "_claims_convergence", lambda r, t: True):
        fails, ids = run(SUITE)
    names = {i.split("::")[-1] for i in ids}
    ok = any("starved" in n for n in names) and any("warns" in n or "audible" in n.lower()
                                                    or "escalat" in n or "message" in n
                                                    for n in names)
    return report("MA2 predicate always True",
                  "the starved fixed-source leg AND the warning legs red",
                  ok, f"fails={fails}: {sorted(names)}")


# ── MA3: the warning is never emitted (silent truncation returns) ──────
def ma3():
    with mock.patch.object(sv, "_warn_if_unconverged", lambda *a, **k: None):
        fails, ids = run(SUITE)
    names = {i.split("::")[-1] for i in ids}
    ok = any("warns" in n for n in names) and any("escalat" in n for n in names)
    return report("MA3 warning suppressed", "the audibility legs red",
                  ok, f"fails={fails}: {sorted(names)}")


# ── MA4: the warning fires ALWAYS (noise — would get filtered away) ────
def ma4():
    def always(history, *, where, budget_name, budget, tol):
        import warnings as w
        from orpheus.numerics.convergence import ConvergenceWarning
        w.warn(f"{where}: hit {budget_name}={budget} tol={tol:.3e} "
               f"(last residual 0.0) BEST-EFFORT", ConvergenceWarning, stacklevel=3)
    with mock.patch.object(sv, "_warn_if_unconverged", always):
        fails, ids = run(SUITE)
    names = {i.split("::")[-1] for i in ids}
    ok = any("SILENT" in n for n in names)
    return report("MA4 warning always fires", "the converged-is-SILENT leg reds",
                  ok, f"fails={fails}: {sorted(names)}")


# ── MA5: the message loses its actionable content ──────────────────────
def ma5():
    def terse(history, *, where, budget_name, budget, tol):
        if history.converged:
            return
        import warnings as w
        from orpheus.numerics.convergence import ConvergenceWarning
        w.warn("did not converge", ConvergenceWarning, stacklevel=3)
    with mock.patch.object(sv, "_warn_if_unconverged", terse):
        fails, ids = run(SUITE)
    names = {i.split("::")[-1] for i in ids}
    ok = any("message" in n for n in names)
    return report("MA5 uninformative message", "the message-content leg reds",
                  ok, f"fails={fails}: {sorted(names)}")


results = [("MA0", ma0())]
if not results[0][1]:
    print("!!!! CONTROL FAILED — the harness or the suite is broken; "
          "no negative below is trustworthy (vv anti-#17).")
for fn, tag in ((ma1, "MA1"), (ma2, "MA2"), (ma3, "MA3"), (ma4, "MA4"), (ma5, "MA5")):
    results.append((tag, fn()))

print()
n = sum(1 for _, ok in results if ok)
print(f"=== convergence-contract mutation battery: {n}/{len(results)} as designed ===")
