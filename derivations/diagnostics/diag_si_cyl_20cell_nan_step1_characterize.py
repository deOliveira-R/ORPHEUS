"""Diagnostic: characterize the SI-CYL-20cell-LS-S8 NaN bug.

Created by numerics-investigator on 2026-05-28.

Reproduces the user-reported failure: solve_sn on a homogeneous
reflective 1-D cylinder with thickness=2.0, n_cells=20, level-symmetric
S8 quadrature, and inner_solver='source_iteration' returns
keff = NaN.  Krylov on the SAME configuration returns the analytical
k_inf = 1.875.

If this test catches a real bug (it does — see
``derivations/diagnostics/diag_si_cyl_20cell_nan_step5_root_cause.py``
for the closed-form root-cause proof), promote to
``tests/sn/test_si_cyl_20cell_nan_regression.py`` as the permanent
regression catcher.

ROOT CAUSE (proven in step5 and step6):
the Blelloch closed-form ``ordinate_scan`` in
``orpheus/sn/spatial/scan.py:138`` computes
``cumprod_a · (psi_0 + cumsum(b / cumprod_a))`` which divides by
``cumprod_a``.  At the pole cell of the cylindrical mesh, ``A_down = 0``
exactly; this drives the per-cell attenuation
``a = 2|μ|·A_total / denom − 1`` toward zero, and at the specific
algebraic resonance ``2|μ|·A_total = dA_w·c_out + Σ_t·V`` (which
holds for μ_x = 1/√20, dr = 0.1, Σ_t = 1.0 of mixture A group 1),
``a`` becomes EXACTLY 0.  The cumprod collapses; the divide produces
NaN; the NaN propagates through ``np.cumsum`` and out of the sweep.

The math is fine — the explicit per-cell recurrence
``psi[i+1] = a·psi[i] + b[i]`` is well-defined even at a=0
(yielding ``psi = b``).  The bug is in the chosen *closed-form*
algorithm, not the math.  Krylov uses
``transport_operator_matvec`` and does NOT go through ``ordinate_scan``,
which is why it works on the same problem.

This is a numerical-stability failure of the Blelloch reduction, not a
discrete-ordinates failure.
"""
import warnings
warnings.filterwarnings('ignore')

import numpy as np
import pytest

from orpheus.geometry import BC, Mesh1D, Region, RegionMesh, StructuredGeometry
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.solver import solve_sn
from orpheus.derivations.common.xs_library import get_mixture


def _build_failing_case(n_cells: int, thickness: float = 2.0):
    fuel = get_mixture('A', '2g')
    geom = StructuredGeometry(
        geometry='CYL',
        regions=(Region(mat_id=0, outer_thickness_cm=thickness),),
        bcs=(BC.reflective,),
    )
    mesh = Mesh1D.from_geometry(geom, region_meshes=(RegionMesh(n_cells=n_cells),))
    quad = Quadrature.level_symmetric(sn_order=8)
    return fuel, mesh, quad


def test_si_cyl_20cell_ls_s8_produces_nan_keff():
    """The reported bug: SI on (CYL, thick=2.0, n=20, LS-S8, 2G refl) → NaN.

    This is the canonical reproducer.  When the underlying bug is fixed
    this test should FAIL (the keff will become 1.875), which is the
    signal to promote this test with an inverted assertion.
    """
    fuel, mesh, quad = _build_failing_case(n_cells=20)
    res = solve_sn(
        materials={0: fuel}, mesh=mesh, quadrature=quad,
        inner_solver='source_iteration',
        keff_tol=1e-12, flux_tol=1e-10,
        max_outer=5, max_inner=5,    # NaN appears in first inner iteration
    )
    # Bug signature: keff is NaN.  Documenting current buggy behaviour.
    assert np.isnan(res.keff), (
        f"keff was {res.keff}, expected NaN (bug signature).  If this "
        f"assertion fires, the bug is FIXED — invert this assertion and "
        f"promote to tests/sn/test_si_cyl_20cell_nan_regression.py."
    )


def test_si_cyl_20cell_krylov_is_correct():
    """Cross-check: SAME configuration with Krylov inner solver works."""
    fuel, mesh, quad = _build_failing_case(n_cells=20)
    res = solve_sn(
        materials={0: fuel}, mesh=mesh, quadrature=quad,
        inner_solver='krylov',
        keff_tol=1e-12, flux_tol=1e-10,
        max_outer=50, max_inner=50,
    )
    assert np.isfinite(res.keff), f"Krylov should not return NaN; got {res.keff}"
    assert abs(res.keff - 1.875) < 1e-6, (
        f"Krylov k_eff = {res.keff}, expected k_inf = 1.875"
    )


@pytest.mark.parametrize("n_cells,expected_finite", [
    (10, True),   # works
    (15, True),   # works
    (19, True),   # works (just below the resonance)
    (20, False),  # BUG — NaN
    (21, True),   # works (just above)
    (40, True),   # works
])
def test_si_cyl_nan_is_sharp_resonance_at_n_cells_20(n_cells, expected_finite):
    """The NaN appears only at the exact (thickness=2.0, n_cells=20)
    combination.  Larger or smaller meshes converge cleanly.  This rules
    out spectral-radius divergence and pins the failure as an algebraic
    cancellation at one cell."""
    fuel, mesh, quad = _build_failing_case(n_cells=n_cells)
    res = solve_sn(
        materials={0: fuel}, mesh=mesh, quadrature=quad,
        inner_solver='source_iteration',
        keff_tol=1e-12, flux_tol=1e-10,
        max_outer=3, max_inner=3,
    )
    is_finite = np.isfinite(res.keff)
    assert is_finite == expected_finite, (
        f"n_cells={n_cells}: keff finite = {is_finite}, expected {expected_finite}"
    )
