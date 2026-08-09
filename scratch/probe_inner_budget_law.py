"""#340 — is ``max_inner`` DERIVABLE, or is a magic integer the honest answer?

The tree carries five hardcoded budgets and the split is not random: it
co-varies exactly with ``inner_tol``.

    solver.py:1001 SNSolver.__init__            200   @ inner_tol 1e-8
    solver.py:2111 solve_sn                     200   @ inner_tol 1e-8
    solver.py:2541 solve_sn_adjoint             200   @ inner_tol 1e-8
    solver.py:2651 solve_sn_adjoint_fixed_src  1000   @ inner_tol 1e-12
    solver.py:3089 solve_sn_fixed_source       1000   @ inner_tol 1e-12

So somebody DID reason about it — four decades of extra tolerance bought
5x the budget.  The question this probe settles is whether 5x is the
right factor, and whether the whole number can come from a law instead
of a constant.

THE HYPOTHESIS.  Source iteration converges geometrically: the residual
decays as ``r_k = r_0 rho^k``, so reaching ``tol`` takes

    n(tol) = [log(r_0) - log(tol)] / [-log(rho)]                     (1)

which is AFFINE in ``|log tol|`` with slope ``1/|log rho|``.  If (1)
holds, then a budget is not a magic number: it is "the iterations needed
to reach the tolerance the caller asked for, at the worst rate we
promise to serve", and the only free parameter is that worst ``rho``.

Under (1) the tolerance factor between the two families is
``log(1e-12)/log(1e-8) = 1.5``, NOT 5 -- so the eigenvalue entries are
either relatively starved or the fixed-source ones are over-provisioned,
by ~3.3x.  Which one is a policy call; that the two disagree is not.

THE INSTRUMENT.  One run per configuration at a very tight tolerance,
then read ``history.flux_residuals`` and ask "at which k did the
residual first drop below X?".  That is EXACTLY the solver's own stop
test (``_claims_convergence`` is ``residual_history[-1] < tol``), so the
per-tolerance counts are exact re-reads of one descent, not predictions
from a fit -- and one run answers every tolerance at once.

vv anti-#17: row 1 is a CONTROL that must reproduce the independently
measured 1631 sweeps of the d3 budget study.  Read no other row until
it does.
"""
from __future__ import annotations

import numpy as np

from orpheus.derivations.common.xs_library import make_mixture
from orpheus.geometry import BC
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.solver import solve_sn_fixed_source
from orpheus.transport.mesh.axis import AxisMesh

V, R = BC("vacuum"), BC("reflective")
TOLS = (1e-6, 1e-8, 1e-10, 1e-12)
FLOOR_TOL = 1e-14          # run to here; read the looser tolerances off it
BUDGET = 40000


def axes(extents, cells, bcs):
    return tuple(
        AxisMesh(edges=np.linspace(0.0, e, n + 1), bc_low=lo, bc_high=hi)
        for e, n, (lo, hi) in zip(extents, cells, bcs)
    )


def absorber(sig_t=(0.8, 1.6)):
    sig_t = np.asarray(sig_t, dtype=float)
    return make_mixture(
        sig_t=sig_t, sig_c=sig_t, sig_f=np.zeros(2), nu=np.zeros(2),
        chi=np.zeros(2), sig_s=np.zeros((2, 2)),
    )


def scatterer(c=0.5, sig_t=(0.8, 1.6)):
    sig_t = np.asarray(sig_t, dtype=float)
    s = np.zeros((2, 2))
    s[0, 0] = c * sig_t[0] * 0.7
    s[1, 0] = c * sig_t[0] * 0.3          # 0 -> 1 downscatter: real coupling
    s[1, 1] = c * sig_t[1]
    return make_mixture(
        sig_t=sig_t, sig_c=sig_t - s.sum(axis=0), sig_f=np.zeros(2),
        nu=np.zeros(2), chi=np.zeros(2), sig_s=s,
    )


CONFIGS = [
    # (label, extents, cells, bcs, mixture)
    ("d3 all-refl absorber [CTRL n=1631]", (1.0, 2.0, 3.0), (3, 4, 5),
     [(R, R)] * 3, absorber()),
    ("d1 all-refl absorber", (1.0,), (3,), [(R, R)], absorber()),
    ("d2 all-refl absorber", (1.0, 2.0), (3, 4), [(R, R)] * 2, absorber()),
    ("d3 all-refl c=0.5", (1.0, 2.0, 3.0), (3, 4, 5),
     [(R, R)] * 3, scatterer()),
    ("d3 all-refl absorber St/4", (1.0, 2.0, 3.0), (3, 4, 5),
     [(R, R)] * 3, absorber(sig_t=(0.2, 0.4))),
    # The realistic regime the DEFAULT actually has to serve: some leakage.
    ("d3 x-vacuum c=0.5", (1.0, 2.0, 3.0), (3, 4, 5),
     [(V, V), (R, R), (R, R)], scatterer()),
    ("d3 all-vacuum c=0.5", (1.0, 2.0, 3.0), (3, 4, 5),
     [(V, V)] * 3, scatterer()),
    ("d2 all-vacuum c=0.9 (scatter-hard)", (1.0, 2.0), (3, 4),
     [(V, V)] * 2, scatterer(c=0.9)),
]

quad = Quadrature.level_symmetric(sn_order=4)
W = float(np.sum(quad.weights))


def first_below(residuals, tol):
    """The iteration index the solver's own stop test would have exited at."""
    for k, r in enumerate(residuals, start=1):
        if r < tol:
            return k
    return None


hdr = f"{'config':<36} " + " ".join(f"{f'n(1e-{p})':>9}"
                                    for p in (6, 8, 10, 12)) + \
      f" {'rho':>8} {'slope':>7} {'R2':>6}"
print(hdr)
print("-" * len(hdr))

rows = []
for label, extents, cells, bcs, mix in CONFIGS:
    shape = (quad.N, 2, *cells)
    q = np.broadcast_to(
        (np.array([1.0, 0.5]) / W)[None, :, *([None] * len(cells))], shape
    ).copy()
    sol = solve_sn_fixed_source(
        {0: mix}, axes(extents, cells, bcs), quad,
        external_source=q, inner_solver="source_iteration",
        inner_schedule="gauss_seidel", inner_tol=FLOOR_TOL, max_inner=BUDGET,
    )
    res = list(sol.history.flux_residuals)
    if not sol.history.converged:
        print(f"{label:<36} did NOT reach {FLOOR_TOL:.0e} in {BUDGET} "
              f"(last {res[-1]:.3e}) -- row unusable")
        continue

    counts = [first_below(res, t) for t in TOLS]

    # rho from the ASYMPTOTIC tail (last 60%), where the slowest mode
    # dominates; the early transient is a mixture of modes and would bias
    # the rate optimistic.
    tail = np.array(res[int(0.4 * len(res)):])
    k = np.arange(tail.size, dtype=float)
    lg = np.log(tail)
    A = np.vstack([k, np.ones_like(k)]).T
    (m, b), *_ = np.linalg.lstsq(A, lg, rcond=None)
    pred = A @ np.array([m, b])
    ss_res = float(np.sum((lg - pred) ** 2))
    ss_tot = float(np.sum((lg - lg.mean()) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    rho = float(np.exp(m))

    # Is n AFFINE in |log tol|, as (1) predicts?  Fit n vs |ln tol|; the
    # slope should equal 1/|ln rho|.
    ok = [(t, c) for t, c in zip(TOLS, counts) if c is not None]
    if len(ok) >= 2:
        x = np.array([-np.log(t) for t, _ in ok])
        y = np.array([float(c) for _, c in ok])
        B = np.vstack([x, np.ones_like(x)]).T
        (sl, ic), *_ = np.linalg.lstsq(B, y, rcond=None)
        py = B @ np.array([sl, ic])
        sst = float(np.sum((y - y.mean()) ** 2))
        nr2 = 1.0 - float(np.sum((y - py) ** 2)) / sst if sst > 0 else float("nan")
    else:
        sl, nr2 = float("nan"), float("nan")

    cells_txt = " ".join(f"{(c if c is not None else '--'):>9}" for c in counts)
    print(f"{label:<36} {cells_txt} {rho:>8.5f} {sl:>7.1f} {nr2:>6.4f}")
    rows.append((label, counts, rho, sl, 1.0 / abs(np.log(rho))))

print()
print("SLOPE CHECK -- law (1) says the n-vs-|ln tol| slope IS 1/|ln rho|.")
print(f"{'config':<36} {'fitted':>9} {'1/|ln rho|':>11} {'ratio':>7}")
print("-" * 66)
for label, _c, rho, sl, inv in rows:
    print(f"{label:<36} {sl:>9.1f} {inv:>11.1f} {sl / inv:>7.3f}")

print()
print("WHAT THE TREE'S TWO DEFAULTS BUY, per configuration:")
print(f"{'config':<36} {'need@1e-8':>10} {'have 200':>9} "
      f"{'need@1e-12':>11} {'have 1000':>10}")
print("-" * 80)
for label, counts, *_ in rows:
    n8, n12 = counts[1], counts[3]
    v8 = "OK" if n8 is not None and n8 <= 200 else "SHORT"
    v12 = "OK" if n12 is not None and n12 <= 1000 else "SHORT"
    print(f"{label:<36} {str(n8):>10} {v8:>9} {str(n12):>11} {v12:>10}")
