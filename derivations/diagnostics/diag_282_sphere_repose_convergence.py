"""Diagnostic: the #282 route-(a) sphere ψ½ seed re-pose is PRINCIPLED — the
heterogeneous keff converges in BOTH h (space) and N (angular order), and the seed
choice changes only the low-angular-order truncation, not the converged value.

Created by numerics-investigator on 2026-07-04.  (Route (a) landed as commit
``a29ab2d`` — "2.5d d3 — #282 route (a) direct starting-direction seed + carrier" —
DURING this investigation; the pre-carve OLD scheme was ``5170f20`` with the seed
arms dormant = edge-extrapolation.)

Context.  Route (a) replaced the OLD edge-extrapolation ψ½ seed with a direct
inward Carlson march of the closed μ=−1 starting-direction equation.  The two seed
treatments give DIFFERENT keff at a fixed low angular order (GL8): the heterogeneous
1g fuel|mod reflective sphere h→0 limit is 0.73825 (NEW) vs 0.73654 (OLD) — a 1.7e-3
gap that does NOT close under h-refinement (so the naive "do they share keff(h→0)?"
test FALSELY reads regression).  But sweeping GL order at a fixed fine mesh, OLD and
NEW CONVERGE TO THE SAME transport answer (measured n=80: NEW/OLD agree to 1.5e-6 at
N=32, 2.7e-6 at N=64, both → ~0.73681; the dd_regression 2g/3reg config agrees to 8e-8
at N=64).  The seed is an angular closure; changing it changes O(N) truncation, not
the fully-refined value — so the re-pose is principled, not a regression, and the
dd_regression ``sphere_2g_3reg_dd_n40`` snapshot move (N=8) is a legitimate re-baseline.

MMS is BLIND to this (every curvilinear MMS ansatz is ≤ linear-in-μ = the seed's
EXACT regime, vv-principles Mode 7 / numerics-investigator lessons L14/L15), so MMS
O(h²) does NOT certify the seed — the angular-order sweep below does.

test_a is the ROBUSTIFIED replacement for the fragile
``test_heterogeneous_1g_spatial_convergence`` ladder (its diff_2<diff_1 check on
n∈[5,10,20] trips on the REAL n=5→10 near-coincidence Δ≈8e-7, which persists at
keff_tol=1e-12 so it is a discretization feature, not iteration noise).  test_b is
the general angular-consistency property.  Promote to
``tests/sn/eigenvalue/test_keff_curvilinear.py``.
"""
import numpy as np
import pytest

from orpheus.geometry import CoordSystem
from orpheus.numerics.quadrature import Quadrature
from orpheus.derivations.common.xs_library import get_mixture
from orpheus.sn.solver import solve_sn
from tests.sn._test_helpers import curvilinear_two_region_mesh as _two_region_mesh


def _sphere_keff(n_cells, n_gl):
    materials = {2: get_mixture("A", "1g"), 0: get_mixture("B", "1g")}
    mesh = _two_region_mesh(
        outers=(0.5, 1.0), mat_ids=(2, 0), n_cells=(n_cells, n_cells),
        coord=CoordSystem.SPHERICAL,
    )
    return solve_sn(materials, mesh, Quadrature.gauss_legendre(n_gl),
                    max_outer=800, max_inner=500, inner_tol=1e-10,
                    keff_tol=1e-12, flux_tol=1e-11).keff


@pytest.mark.slow
def test_a_spatial_ladder_monotone_from_n10():
    """Robustified spatial-convergence ladder: from n=10 the keff increments shrink
    monotonically (the n=5→10 near-coincidence that trips the legacy n∈[5,10,20]
    ladder is a real coarse-mesh feature, not non-convergence)."""
    keffs = [_sphere_keff(n, 8) for n in (10, 20, 40)]
    d1, d2 = abs(keffs[1] - keffs[0]), abs(keffs[2] - keffs[1])
    assert d2 < d1, (
        f"keff not converging from n=10: Δ(20−10)={d1:.3e}, Δ(40−20)={d2:.3e}, "
        f"keffs={[f'{k:.8f}' for k in keffs]}"
    )


@pytest.mark.slow
def test_b_angular_order_consistency():
    """The seed re-pose is a CONSISTENT angular closure: at a fixed spatial mesh the
    keff increments shrink as the angular order N grows (converging to a limit).  A
    seed that did NOT converge in N would be an inconsistent (regression) closure."""
    n_cells = 20
    keffs = [_sphere_keff(n_cells, n_gl) for n_gl in (8, 16, 32)]
    d1, d2 = abs(keffs[1] - keffs[0]), abs(keffs[2] - keffs[1])
    assert d2 < d1 and d2 < 5e-4, (
        f"keff not converging in angular order N: Δ(16−8)={d1:.3e}, "
        f"Δ(32−16)={d2:.3e}, keffs={[f'{k:.8f}' for k in keffs]}"
    )
