"""Diagnostic Step 6: tighten outer/inner tolerances across mesh refinement.

Created by numerics-investigator on 2026-05-28.

Step 5 found Krylov error GROWS with refinement at defaults:
  n=10: 1.9153 (4%)
  n=20: 1.4019 (25%)

But ``test_kinf_homogeneous[sphere-2eg-krylov]`` passes at n_cells=10
using ``inner_tol=1e-12, max_inner=1000, keff_tol=1e-14, flux_tol=1e-12``
(the L1 anchor's ``_TIGHT_KW`` dict).

This step tests whether the L1 anchor's tight tolerances rescue the
failing test config (n_cells=20).  If YES → the issue is a slow but
correct convergence + an outer convergence-criterion bug that calls
"converged" too early.  If NO → there's a structural defect that
tightening doesn't fix.

Step 1 already established that ``inner_tol=1e-8 vs 1e-12`` gives the
SAME wrong keff (1.4019).  So tightening just inner_tol won't help.
But maybe tightening outer also matters?
"""
from __future__ import annotations

import warnings


def _solve(*, n_cells, inner_solver, **kw):
    from orpheus.derivations.common.xs_library import get_mixture
    from orpheus.geometry import (
        BC, Mesh1D, Region, RegionMesh, StructuredGeometry,
    )
    from orpheus.numerics.quadrature import Quadrature
    from orpheus.sn.solver import solve_sn

    fuel = get_mixture("A", "2g")
    geom = StructuredGeometry(
        geometry="SPH",
        regions=(Region(mat_id=0, outer_thickness_cm=2.0),),
        bcs=(BC.reflective,),
    )
    mesh = Mesh1D.from_geometry(
        geom, region_meshes=(RegionMesh(n_cells=n_cells),),
    )
    quad = Quadrature.gauss_legendre(n_ordinates=8)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        res = solve_sn(
            materials={0: fuel}, mesh=mesh, quadrature=quad,
            inner_solver=inner_solver, **kw,
        )
    return res.keff


def test_step6_tight_tol_sweep():
    """Test if L1 anchor's tight tolerances fix Krylov at n_cells=20."""
    kinf = 1.875
    # L1 anchor's TIGHT_KW
    tight_kw = dict(
        max_outer=1000, keff_tol=1e-14, flux_tol=1e-12,
        max_inner=1000, inner_tol=1e-12,
    )
    # Failing test's looser params
    loose_kw = dict(
        max_outer=500, keff_tol=1e-12, flux_tol=1e-10,
        max_inner=200, inner_tol=1e-8,
    )

    print()
    print("=== Krylov at n_cells=10 ===")
    for label, kw in [("tight", tight_kw), ("loose", loose_kw)]:
        keff = _solve(n_cells=10, inner_solver="krylov", **kw)
        print(f"  {label}: keff={keff:.10f}  err={abs(keff - kinf):.3e}")

    print("\n=== Krylov at n_cells=20 ===")
    for label, kw in [("tight", tight_kw), ("loose", loose_kw)]:
        keff = _solve(n_cells=20, inner_solver="krylov", **kw)
        print(f"  {label}: keff={keff:.10f}  err={abs(keff - kinf):.3e}")


if __name__ == "__main__":
    test_step6_tight_tol_sweep()
