"""#340 — is the EIGENVALUE path's inner truncation visible at all?

The `d9b027d7` carve wired `ConvergenceWarning` at the four public
entries, but not symmetrically:

    solve_sn                    warns on max_outer  (tol = keff_tol)
    solve_sn_adjoint            warns on max_outer  (tol = keff_tol)
    solve_sn_adjoint_fixed_src  warns on max_inner  (tol = inner_tol)
    solve_sn_fixed_source       warns on max_inner  (tol = inner_tol)

The fixed-source entries have only ONE loop, so one warning covers
them.  The eigenvalue entries have TWO, and only the outer is watched.
`Solution.history.converged` on that path is `outcome.converged` -- the
POWER iteration's flag -- which is an honest claim about the outer and
says nothing about the inners underneath it.  Meanwhile
`_certify_within_group_exit` is a deliberate no-op precisely when an
inner made no claim.  So the question:

    can a k-eigenvalue solve truncate every inner solve, return a k
    that is WRONG by more than keff_tol, report converged=True, and
    emit NOTHING?

If yes that is the #340 defect one level in, and the first pass did not
close it.

FIXTURE CHOICE MATTERS (vv anti-#3).  An all-reflective homogeneous box
has k = k_inf = nu*Sigma_f / Sigma_a, which is FLUX-SHAPE INDEPENDENT --
inner truncation cannot move it, and the probe would report a
comfortable false negative.  So this uses a LEAKING box (x-vacuum,
y-reflective): k then depends on the shape the inner solve is
responsible for producing.  The all-reflective row is kept as a
deliberate contrast, to show the degeneracy rather than hide in it.
"""
from __future__ import annotations

import warnings

import numpy as np

from orpheus.derivations.common.xs_library import get_mixture
from orpheus.geometry import BC
from orpheus.numerics.convergence import ConvergenceWarning
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.solver import solve_sn
from orpheus.transport.mesh.axis import AxisMesh

V, R = BC("vacuum"), BC("reflective")
quad = Quadrature.level_symmetric(sn_order=4)
KEFF_TOL = 1e-7


def axes(bcs):
    return tuple(
        AxisMesh(edges=np.linspace(0.0, e, n + 1), bc_low=lo, bc_high=hi)
        for e, n, (lo, hi) in zip((2.0, 1.5), (4, 3), bcs)
    )


def run(bcs, max_inner):
    """Return (keff, converged_flag, n_ConvergenceWarnings)."""
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        sol = solve_sn(
            {0: get_mixture("A", "2g")}, axes(bcs), quad,
            max_outer=500, keff_tol=KEFF_TOL, flux_tol=1e-6,
            max_inner=max_inner, inner_tol=1e-10,
        )
    n_cw = sum(1 for w in caught if issubclass(w.category, ConvergenceWarning))
    return float(sol.keff), bool(sol.history.converged), n_cw


for label, bcs in (
    ("LEAKING  (x-vacuum, y-reflective)", [(V, V), (R, R)]),
    ("all-reflective (k = k_inf, shape-BLIND -- contrast row)", [(R, R)] * 2),
):
    print(f"\n=== {label} ===")
    k_ref, conv_ref, cw_ref = run(bcs, max_inner=2000)
    print(f"  reference   max_inner=2000 : k={k_ref:.12f} "
          f"converged={conv_ref} warnings={cw_ref}")
    for mi in (1, 2, 3, 5, 10):
        k, conv, cw = run(bcs, max_inner=mi)
        dk = abs(k - k_ref)
        flag = "SILENT+WRONG" if (conv and cw == 0 and dk > KEFF_TOL) else (
            "audible" if cw else ("honest" if not conv else "ok"))
        print(f"  starved     max_inner={mi:<5}: k={k:.12f} "
              f"converged={conv} warnings={cw} |dk|={dk:.3e} "
              f"({dk / KEFF_TOL:>8.1f} x keff_tol)  -> {flag}")
