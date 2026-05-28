"""Diagnostic Step 5: mesh-refinement scaling at fixed material/length.

Created by numerics-investigator on 2026-05-28.

Step 4 finding: L1 anchor flow with mixture A 2g at length=2.0 sphere:
  n_cells=10 → keff = 1.8750 (PASS, matches k_inf)
  n_cells=20 → keff = 1.5582 (FAIL, 17% drift)

The bug is mesh-refinement induced.  k_inf is geometry-/mesh-independent
for homogeneous reflective.  This is the smoking gun for a discretisation
inconsistency.

Step 5 establishes the convergence trend across n_cells {5, 8, 10, 12,
14, 16, 20, 30}.

The brief states SI converges correctly at n_cells=20.  So this is NOT a
"problem becomes hard at fine mesh"; it's "Krylov-on-sphere has a
mesh-dependent defect that SI doesn't share".
"""
from __future__ import annotations

import warnings

import numpy as np


def _kinf_via(*, n_cells: int, inner_solver: str, inner_tol: float = 1e-8):
    """solve_sn at given n_cells, return keff."""
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
            inner_solver=inner_solver,
            keff_tol=1e-12, flux_tol=1e-10,
            inner_tol=inner_tol,
        )
    return res.keff


def test_step5_mesh_refinement_si_vs_krylov():
    """Sweep n_cells; tabulate SI and Krylov keff."""
    kinf = 1.875
    n_list = [5, 8, 10, 12, 14, 16, 20]

    print()
    print(f"  {'n_cells':>8s} {'SI_keff':>15s} {'SI_err':>10s} {'KR_keff':>15s} {'KR_err':>10s}")
    for n in n_list:
        keff_si = _kinf_via(n_cells=n, inner_solver="source_iteration")
        keff_kr = _kinf_via(n_cells=n, inner_solver="krylov")
        err_si = abs(keff_si - kinf)
        err_kr = abs(keff_kr - kinf)
        print(f"  {n:>8d} {keff_si:>15.10f} {err_si:>10.3e} {keff_kr:>15.10f} {err_kr:>10.3e}")


if __name__ == "__main__":
    test_step5_mesh_refinement_si_vs_krylov()
