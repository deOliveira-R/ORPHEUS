"""Is the d3 pure-absorber red a CONVERGENCE artifact or a real BIAS?

The gate asserts psi_n = Q/(W*Sigma_t) at rtol=1e-10 and measures 3.29e-10.
Discriminator: sweep inner_tol.  If the error tracks inner_tol, the gate is
measuring the solver's exit tolerance (its docstring's own "solver-tolerance-
tight" bound) and the quadrature change moved the iteration's rho.  If the
error sits on a FIXED floor independent of inner_tol, it is a real bias in
the discretization and a Cardinal-Rule-1 bug.
"""
from __future__ import annotations

import importlib.util
import sys

import numpy as np

from orpheus.derivations.common.xs_library import make_mixture
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.solver import solve_sn_fixed_source

spec = importlib.util.spec_from_file_location(
    "_d3mod", "tests/sn/solve/test_d3_admission.py",
)
mod = importlib.util.module_from_spec(spec)
sys.modules["_d3mod"] = mod
spec.loader.exec_module(mod)

quad = Quadrature.level_symmetric(sn_order=4)
mix = make_mixture(
    sig_t=np.array([0.8, 1.6]),
    sig_c=np.array([0.8, 1.6]),
    sig_f=np.array([0.0, 0.0]),
    nu=np.array([0.0, 0.0]),
    chi=np.zeros(2),
    sig_s=np.zeros((2, 2)),
)
Q_g = np.array([1.0, 0.5])
W = float(np.sum(quad.weights))
q = np.broadcast_to(
    (Q_g / W)[None, :, None, None, None], (quad.N, 2, 3, 4, 5),
).copy()

print(f"{'inner_tol':>12} {'iters':>7} {'max rel err g0':>16} {'g1':>14}")
for tol in (1e-9, 1e-11, 1e-13, 1e-15):
    sol = solve_sn_fixed_source(
        {0: mix}, mod._d3_axes(), quad,
        external_source=q, boundary_condition="reflective",
        inner_tol=tol,
    )
    psi = np.asarray(sol.angular_flux.interior.values)
    sig_t = np.asarray(mix.SigT)
    errs = []
    for g in range(2):
        want = Q_g[g] / (W * sig_t[g])
        errs.append(float(np.max(np.abs(psi[:, g] - want) / want)))
    hist = getattr(sol, "history", None)
    n_it = len(getattr(hist, "residuals", []) or []) if hist is not None else -1
    print(f"{tol:>12.0e} {n_it:>7d} {errs[0]:>16.3e} {errs[1]:>14.3e}")
