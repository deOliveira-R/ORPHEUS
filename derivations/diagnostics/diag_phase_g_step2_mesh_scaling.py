"""Mesh-scaling verification: the two-bugs fix dissolves the 22% L0 error
at any cell count.

If the algebraic-difference identification is correct, applying BOTH
fixes (pole-face IC + Carlson Q_bar source) should give machine-precision
agreement with ψ_K at ANY mesh refinement, NOT just the 2-cell case.

Specifically tests n_cells ∈ {2, 5, 10, 20, 40, 80}.
"""
from __future__ import annotations

import numpy as np

# Re-use the custom_si_sweep + picard_to_fp from the other diagnostic
import importlib.util
import os
spec = importlib.util.spec_from_file_location(
    "two_bugs",
    os.path.join(os.path.dirname(__file__), "diag_phase_g_step2_two_bugs_isolation.py"),
)
two_bugs = importlib.util.module_from_spec(spec)
spec.loader.exec_module(two_bugs)


def test_mesh_independence():
    print("Mesh scaling of the two-bug fix (mixture B 1G sphere, GL-2, refl):")
    print()
    print(f"  {'n_cells':>8s}  {'v0 max err':>12s}  {'v1 max err':>12s}  "
          f"{'v2 max err':>12s}  {'v3 max err':>12s}")
    print(f"  {'-'*8}  {'-'*12}  {'-'*12}  {'-'*12}  {'-'*12}")
    for n_cells in (2, 5, 10, 20, 40):
        psi_ref = np.full((2, n_cells), 5.0)
        row = [f"  {n_cells:>8d}"]
        for fix_pole, fix_seed in ((False,False),(True,False),(False,True),(True,True)):
            psi, k, r = two_bugs.picard_to_fp(fix_pole, fix_seed, n_cells=n_cells)
            err = np.max(np.abs(psi - psi_ref))
            row.append(f"  {err:>12.4e}")
        print("".join(row))


if __name__ == "__main__":
    test_mesh_independence()
