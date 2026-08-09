"""WHERE does the d3 pure-absorber bias live, and is the analytic psi a
solution of the DISCRETE system?

Three questions:
 1. Which ordinates carry the 3.287e-10 error, and what are their cosines?
 2. Is the returned psi the discrete solution (=> the discrete operator is
    inconsistent with the analytic answer) or is the analytic psi the
    discrete solution (=> the solve did not converge to it)?
 3. Does the same bias appear under the PRE-#337 level-symmetric S4 nodes?
    (monkeypatched node values, in-process only -- never a git revert)
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
assert spec is not None and spec.loader is not None
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

sol = solve_sn_fixed_source(
    {0: mix}, mod._d3_axes(), quad,
    external_source=q, boundary_condition="reflective", inner_tol=1e-13,
)
psi = np.asarray(sol.angular_flux.interior.values)      # (N, ng, 3, 4, 5)
sig_t = np.asarray(mix.SigT)

print("=== Q1: per-ordinate error structure, group 0 ===")
want0 = Q_g[0] / (W * sig_t[0])
per_ord = np.max(np.abs(psi[:, 0] - want0), axis=(1, 2, 3)) / want0
mu = np.asarray(quad.mu_x), np.asarray(quad.mu_y), np.asarray(quad.mu_z)
print(f"{'n':>3} {'mu_x':>11} {'mu_y':>11} {'mu_z':>11} {'max rel':>11}")
for n in range(quad.N):
    flag = "  <== OFF" if per_ord[n] > 1e-12 else ""
    print(f"{n:>3} {mu[0][n]:>11.6f} {mu[1][n]:>11.6f} {mu[2][n]:>11.6f} "
          f"{per_ord[n]:>11.3e}{flag}")

n_off = int(np.sum(per_ord > 1e-12))
print(f"\nordinates carrying the bias: {n_off} / {quad.N}")

print("\n=== Q2: sum of weights vs 4*pi (the W the source is normalised by) ===")
print(f"  sum(w)      = {W!r}")
print(f"  4*pi        = {4 * np.pi!r}")
print(f"  rel defect  = {abs(W - 4 * np.pi) / (4 * np.pi):.3e}")

print("\n=== Q3: is the flat analytic psi a FIXED POINT of the same solve? ===")
# Re-solve starting from a converged-looking state is not exposed; instead
# check the per-ordinate balance the analytic answer must satisfy:
#   grad psi = 0  =>  Sigma_t * psi == q  per ordinate, per cell.
for g in (0, 1):
    want = Q_g[g] / (W * sig_t[g])
    resid = np.abs(sig_t[g] * psi[:, g] - q[:, g])
    print(f"  g{g}: max |Sigma_t*psi - q| = {resid.max():.3e}"
          f"   (q = {q[0, g]:.6e}, psi_want = {want:.6e})")

print("\n=== Q4: cell-level spread -- is the bias uniform or spatial? ===")
worst = int(np.argmax(per_ord))
blk = psi[worst, 0]
print(f"  worst ordinate n={worst}, cosines "
      f"({mu[0][worst]:.6f}, {mu[1][worst]:.6f}, {mu[2][worst]:.6f})")
print(f"  psi block min={blk.min():.15e} max={blk.max():.15e}")
print(f"  want          ={want0:.15e}")
print(f"  spread across cells = {(blk.max() - blk.min()) / want0:.3e}")
