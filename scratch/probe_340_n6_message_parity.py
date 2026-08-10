"""#340 N6 — does the RECORD reproduce what the five callers pass by hand?

The N6 carve re-derives `_warn_if_unconverged`'s three hand-passed facts
(`budget_name`, `budget`, `tol`) from `record.first_failure`, so the message can
speak for the level that actually failed instead of for the caller's top level.

Commit 1 leaves the GUARD on `history.converged`, so only solves that already
warn are affected.  For those, the honest question is not "is it better" but
**"what changes"** — and the answer decides which of the four message gates in
`tests/sn/solve/test_convergence_contract.py` need re-baselining, and with what
rationale.

`first_failure` checks SELF before children (the MA3 lesson), so while the guard
stays on `converged` the failing level IS the top level.  Therefore:

    budget  : caller's `budget`   vs  record.budget
    tol     : caller's `tol`      vs  record.binding_criterion.tolerance

⚠ These are NOT the same question.  A level stops on a CONJUNCTION, so it has
one budget but N tolerances — and the caller can only pass one.  The eigenvalue
entries pass `keff_tol`, so whenever `dphi` is the binding criterion the current
message quotes the tolerance of a criterion that is NOT the one that failed.
That is the F2 disease inside the diagnostic itself.
"""
from __future__ import annotations

import numpy as np

from orpheus.derivations.common.xs_library import make_mixture
from orpheus.geometry import BC
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.solver import solve_sn, solve_sn_fixed_source
from orpheus.transport.mesh.axis import AxisMesh

R = BC("reflective")
EXTENTS, CELLS = (1.0, 2.0, 3.0), (3, 4, 5)


def axes():
    return tuple(
        AxisMesh(edges=np.linspace(0.0, e, n + 1), bc_low=R, bc_high=R)
        for e, n in zip(EXTENTS, CELLS)
    )


def absorber(sig_t=(0.8, 1.6)):
    sig_t = np.asarray(sig_t, dtype=float)
    return make_mixture(
        sig_t=sig_t, sig_c=sig_t, sig_f=np.zeros(2), nu=np.zeros(2),
        chi=np.zeros(2), sig_s=np.zeros((2, 2)),
    )


def fissile(c=0.5, sig_t=(0.8, 1.6)):
    sig_t = np.asarray(sig_t, dtype=float)
    s = np.zeros((2, 2))
    s[0, 0] = c * sig_t[0] * 0.7
    s[0, 1] = c * sig_t[0] * 0.3          # [from, to]: 0 -> 1 downscatter
    s[1, 1] = c * sig_t[1]
    sig_f = np.array([0.05, 0.30])
    return make_mixture(
        sig_t=sig_t, sig_c=sig_t - s.sum(axis=1) - sig_f, sig_f=sig_f,
        nu=np.array([2.4, 2.4]), chi=np.array([1.0, 0.0]), sig_s=s,
    )


quad = Quadrature.level_symmetric(sn_order=4)
W = float(np.sum(quad.weights))


def uniform_source():
    shape = (quad.N, 2, *CELLS)
    return np.broadcast_to(
        (np.array([1.0, 0.5]) / W)[None, :, *([None] * len(CELLS))], shape
    ).copy()


CASES = [
    # (label, thunk, caller's budget_name, caller's budget, caller's tol)
    (
        "solve_sn_fixed_source (leaf, starved)",
        lambda: solve_sn_fixed_source(
            {0: absorber()}, axes(), quad, external_source=uniform_source(),
            boundary_condition="reflective", inner_tol=1e-13, max_inner=50,
        ),
        "max_inner", 50, 1e-13,
    ),
    (
        "solve_sn (eigenvalue, starved OUTER)",
        lambda: solve_sn(
            {0: fissile()}, axes(), quad, inner_solver="source_iteration",
            inner_schedule="gauss_seidel", max_outer=3, keff_tol=1e-12,
            flux_tol=1e-12, max_inner=200, inner_tol=1e-8,
        ),
        "max_outer", 3, 1e-12,
    ),
]

print(f"{'case':<40} {'field':<10} {'caller':>14} {'record':>14}  match")
print("-" * 92)
for label, thunk, knob, budget, tol in CASES:
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sol = thunk()
    rec = sol.history.record
    fail = rec.first_failure
    if fail is None:
        print(f"{label:<40} <converged — this case does not warn; fixture drift>")
        continue
    crit = fail.binding_criterion
    rec_tol = None if crit is None else crit.tolerance
    rows = [
        ("level", "<top>", fail.label),
        ("budget", budget, fail.budget),
        ("tol", tol, rec_tol),
        ("binding", "n/a", None if crit is None else crit.name),
    ]
    for i, (field, a, b) in enumerate(rows):
        same = "" if field in ("level", "binding") else (
            "  OK" if a == b else "  ** DIFFERS **")
        head = label if i == 0 else ""
        fa = f"{a:.3e}" if isinstance(a, float) else str(a)
        fb = f"{b:.3e}" if isinstance(b, float) else str(b)
        print(f"{head:<40} {field:<10} {fa:>14} {fb:>14}{same}")
    print(f"{'':<40} all criteria: "
          + ", ".join(f"{c.name}(tol={c.tolerance:.1e}, "
                      f"cleared={c.cleared})" for c in fail.criteria))
    print()
