"""#340 — EXISTENCE-CHECK the plan's own "Done when", against the tree.

`plan-authoring` §1: a resume pointer owes an existence-check per DELIVERABLE,
not only a grep per symbol.  The #340 plan's done-when was written 2026-08-09,
BEFORE N1/N2a/N2b-i/N2b-ii/N3 landed.  Five steps later, the honest question is
not "how do we build it" but "does the tree already say it?".

The predicate, verbatim from `.claude/plans/nested_iteration_diagnostics.md`:

    for the measured d=3 all-reflective eigenvalue solve at max_inner=200, the
    returned object reports "inner(within-group) hit budget=200 at rho=0.985,
    needs ~849 for tol=1e-8" -- and a caller asking the one question "is this
    answer trustworthy?" gets False.

So there are FIVE separable claims, and each is checked on its own line below.
A partial pass is the interesting outcome: it names exactly what N4/N5 still owe.

    Q1  WHICH level failed          -> record.first_failure.label
    Q2  WHICH criterion bound it    -> .binding_criterion.name
    Q3  HOW FAR it got              -> .n_iterations vs .budget, .exhausted_budget
    Q4  WHAT BUDGET would suffice   -> .rate (rho) and .projected_iterations()
    Q5  IS IT TRUSTWORTHY           -> solution.history.fully_converged is False

CONFIGURATION (plan-authoring §4 -- a number without its fixture is not
reusable).  d=3, extents (1.0, 2.0, 3.0), cells (3, 4, 5), all six faces
reflective, LS-S4, 2 groups, `max_inner=200`, `inner_tol=1e-8`, source
iteration on a gauss_seidel schedule.  Zero leakage is the point: it is the
regime the plan measured rho ~ 0.985 in, and the regime both shipped budgets
are SHORT for (F4).
"""
from __future__ import annotations

import numpy as np

from orpheus.derivations.common.xs_library import make_mixture
from orpheus.geometry import BC
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.solver import solve_sn
from orpheus.transport.mesh.axis import AxisMesh

R = BC("reflective")
EXTENTS, CELLS = (1.0, 2.0, 3.0), (3, 4, 5)
MAX_INNER, INNER_TOL = 200, 1e-8


def axes():
    return tuple(
        AxisMesh(edges=np.linspace(0.0, e, n + 1), bc_low=R, bc_high=R)
        for e, n in zip(EXTENTS, CELLS)
    )


def fissile(c=0.5, sig_t=(0.8, 1.6)):
    r"""A 2-group fissile mixture with real down-scatter coupling.

    ⛔ **CORRECTED 2026-08-10 — this function shipped with the scattering
    matrix TRANSPOSED, and the probe's own headline number came from a
    hand-corrected variant that was never written back here.**  The bug and
    all three of its symptoms, measured:

    ``make_mixture`` reads ``sig_s[g_from, g_to]``; the original wrote
    ``[g_to, g_from]`` (its own comment said ``s[1,0] = 0 -> 1 downscatter``,
    which is the transposed reading) and removed out-scatter with
    ``s.sum(axis=0)`` — the correct row sum for THAT convention, so the two
    errors were self-consistent and neither showed up as an exception:

    * ``sigma_t = [0.8, 1.6]`` against ``sigma_c + sigma_f + outscatter =
      [0.68, 1.32]`` — **inconsistent by [+0.12, -0.52]**.  An inconsistent
      total makes the two zero-leakage expressions for :math:`k` *different
      numbers*, which is why the probe read ``0.2999999999999999`` while the
      0-D reference read ``0.2307692307692307``.  Neither solver was wrong.
    * ``SigS[0, 1] = 0.0`` — group 0 never feeds group 1, so ``phi_1 = 0``
      identically and the "2-group" fixture was **1-group in disguise**
      (`vv-principles` anti-pattern #3: a 1-group eigenvalue proves nothing
      about transport, because :math:`k = \nu\Sigma_f/\Sigma_a` is
      flux-shape independent).
    * ``fissile(c=0.9)`` gave ``sigma_c = -0.14`` — a negative capture
      cross-section, admitted silently.

    The repaired fixture DOES show the pole this probe exists to
    demonstrate: `[M]` ``keff = 0.43846153845055`` against
    ``solve_homogeneous_infinite(mix).k_inf = 0.43846153846154``, i.e.
    :math:`|\Delta k| = 1.1\times10^{-11}`, with 4 of 4 inners TRUNCATED
    200/200 at :math:`\rho \approx 0.985` — a truncation that is genuinely
    **benign for keff**.

    ⭐ The lesson is not "transposes are easy to get wrong".  It is that a
    mixture whose channels do not sum to its total is an ILLEGAL state that
    every consumer reads differently, and nothing refused it.  A three-line
    consistency check belongs in the test suite.
    """
    sig_t = np.asarray(sig_t, dtype=float)
    s = np.zeros((2, 2))                  # [g_from, g_to] — the library's order
    s[0, 0] = c * sig_t[0] * 0.7
    s[0, 1] = c * sig_t[0] * 0.3          # 0 -> 1 DOWNSCATTER: row = source group
    s[1, 1] = c * sig_t[1]
    sig_f = np.array([0.05, 0.30])
    sig_c = sig_t - s.sum(axis=1) - sig_f  # row sum = total out-scatter from g
    if not (sig_c > 0.0).all():
        raise ValueError(f"unphysical mixture: sigma_c = {sig_c}")
    return make_mixture(
        sig_t=sig_t, sig_c=sig_c,
        sig_f=sig_f, nu=np.array([2.4, 2.4]),
        chi=np.array([1.0, 0.0]), sig_s=s,
    )


sol = solve_sn(
    {0: fissile()}, axes(), Quadrature.level_symmetric(sn_order=4),
    inner_solver="source_iteration", inner_schedule="gauss_seidel",
    max_inner=MAX_INNER, inner_tol=INNER_TOL,
)

hist = sol.history
assert hist is not None
rec = hist.record

print("=" * 72)
print("THE RECORD, AS THE TREE RENDERS IT TODAY")
print("=" * 72)
print(rec.report())
print()
print(f"keff = {sol.keff!r}")
print()

print("=" * 72)
print("THE DONE-WHEN, CLAIM BY CLAIM")
print("=" * 72)

fail = rec.first_failure
print(f"Q1  which level failed      : {fail.label if fail else '<none — all converged>'}")

crit = fail.binding_criterion if fail else None
print(f"Q2  which criterion bound it: {crit.name if crit else '<none>'}")

if fail is not None:
    print(f"Q3  how far it got          : {fail.n_iterations} of budget "
          f"{fail.budget}   exhausted={fail.exhausted_budget}  "
          f"truncated={fail.truncated}  status={fail.status!r}")
    print(f"Q4  rate (rho)              : {fail.rate!r}")
    print(f"    projected_iterations()  : {fail.projected_iterations()!r}")
else:
    print("Q3  how far it got          : n/a — nothing failed")
    print("Q4  rate / projection       : n/a")

print(f"Q5  is it trustworthy?      : fully_converged = {hist.fully_converged!r}"
      f"   (flat `converged` = {hist.converged!r})")

print()
print("=" * 72)
print("WHAT A USER WOULD ACTUALLY SEE — is the warning emitted, and does it")
print("say what to set?  (re-run under `warnings.catch_warnings`)")
print("=" * 72)
import warnings

with warnings.catch_warnings(record=True) as caught:
    warnings.simplefilter("always")
    solve_sn(
        {0: fissile()}, axes(), Quadrature.level_symmetric(sn_order=4),
        inner_solver="source_iteration", inner_schedule="gauss_seidel",
        max_inner=MAX_INNER, inner_tol=INNER_TOL,
    )
if not caught:
    print("  *** NOTHING WARNED ***")
for w in caught:
    print(f"  [{w.category.__name__}] {w.message}")
