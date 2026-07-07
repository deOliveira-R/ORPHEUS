"""Diagnostic: the within-group Krylov restart MUST cover the full augmented
(bulk+trace+seed) composite ravel — not the bulk DOF count — else restarted GMRES
truncates and the sphere-eigenvalue within-group solve stalls (ERR-053, #282/#280).

Created by numerics-investigator on 2026-07-04.  REGRESSION GATE for a bug that was
introduced AND fixed within #282 route (a) (commit ``a29ab2d``):

* THE BUG (pre-fix working tree): ``_solve_krylov`` sized the gmres ``restart`` from
  the BULK formula ``N*ng*prod(spatial_shape)`` (solver.py ``n_dof=...``).  Route (a)
  grew the Krylov state to a 3-block ``TimedFullField(bulk ⊕ trace ⊕
  radial_characteristic)`` whose ``to_flat`` is LARGER (n=10 sphere GL8 1g: bulk 160 <
  composite 210 = +42 seed +8 trace).  Restarted GMRES(160) STAGNATED on the 210-dim
  augmented system (info=300, residual plateau, 868 s; keff best-effort — WRONG=0.865
  under a bounded outer cap).  Bit the eigenvalue path (moderator c=0.95 reflective =
  poorly conditioned, needs >restart iters); the fixed-source c=0.5 path fit within
  one restart cycle so it did not stall.  Distinct from #200 (the identity precond).

* THE FIX (committed ``a29ab2d``): ``n_dof=int(initial_guess.to_flat().size)`` — the
  restart is sized from the FULL composite ravel.  Verified end-to-end: restart=210,
  info=0, k_SI≡k_Krylov to 4.7e-11, 3.4 s.

This gate SPIES the gmres ``restart`` during a real sphere-eigenvalue Krylov solve and
asserts it covers the composite AND no call stalls — it PASSES on the fixed tree and
REDDENS if the sizing regresses to the bulk formula.  ``test_structural_*`` documents
WHY (bulk < composite) with no solve.  Promote to ``tests/sn/solve/`` (Krylov path).
"""
import numpy as np
import pytest
import scipy.sparse.linalg as _spla

from orpheus.geometry import CoordSystem
from orpheus.numerics.quadrature import Quadrature
from orpheus.derivations.common.xs_library import get_mixture
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.sn.solver import solve_sn
from orpheus.transport.timed_full_field import TimedFullField
from orpheus.transport.fields.angular_flux import AngularFlux
from orpheus.transport.fields.angular_boundary_flux import AngularBoundaryFlux
from orpheus.transport.fields.radial_characteristic_flux import RadialCharacteristicFlux
from tests.sn._test_helpers import curvilinear_two_region_mesh as _two_region_mesh


def _composite_dim_and_bulk(sn):
    bulk = sn.quad.N * sn.ng * int(np.prod(sn.spatial_shape))
    composite = TimedFullField.zeros(
        bulk=AngularFlux, boundary=AngularBoundaryFlux, mesh=sn,
        radial_characteristic=RadialCharacteristicFlux,
    )
    return bulk, len(composite.to_flat())


@pytest.mark.parametrize("n_cells", [5, 10, 20])
def test_structural_bulk_formula_undercounts_the_augmented_composite(n_cells):
    """Documents the ERR-053 hazard: the BULK DOF formula is strictly smaller than the
    3-block composite ravel on a carrying (sphere) mesh, so sizing restart from bulk
    would truncate.  No solve — pure structural fact."""
    mesh = _two_region_mesh(outers=(0.5, 1.0), mat_ids=(2, 0),
                            n_cells=(n_cells, n_cells), coord=CoordSystem.SPHERICAL)
    sn = SNMesh(mesh, Quadrature.gauss_legendre(8),
                {2: get_mixture("A", "1g"), 0: get_mixture("B", "1g")})
    bulk, composite = _composite_dim_and_bulk(sn)
    assert bulk < composite, (
        f"expected bulk {bulk} < composite {composite} (seed+trace) — the whole "
        f"point of the ERR-053 hazard; if equal, the seed carrier is absent"
    )


@pytest.mark.slow
def test_production_krylov_restart_covers_composite_and_does_not_stall():
    """REGRESSION GATE: the production sphere-eigenvalue Krylov solve sizes gmres
    restart from the composite ravel (>= its dimension) and NO gmres call stalls
    (info==0).  Reverting ``n_dof`` to the bulk formula drops restart below the
    composite and reddens this (info>0 stall + restart < composite)."""
    mesh = _two_region_mesh(outers=(0.5, 1.0), mat_ids=(2, 0),
                            n_cells=(5, 5), coord=CoordSystem.SPHERICAL)
    sn = SNMesh(mesh, Quadrature.gauss_legendre(8),
                {2: get_mixture("A", "1g"), 0: get_mixture("B", "1g")})
    _bulk, composite = _composite_dim_and_bulk(sn)

    seen = []
    orig = _spla.gmres

    def _spy(A, b, *a, **k):
        seen.append((k.get("restart"), None))
        x, info = orig(A, b, *a, **k)
        seen[-1] = (k.get("restart"), info)
        return x, info

    _spla.gmres = _spy
    try:
        result = solve_sn(
            {2: get_mixture("A", "1g"), 0: get_mixture("B", "1g")}, mesh,
            Quadrature.gauss_legendre(8), inner_solver="krylov",
            max_outer=200, max_inner=4000, inner_tol=1e-10,
            keff_tol=1e-9, flux_tol=1e-8,
        )
    finally:
        _spla.gmres = orig

    assert seen, "no gmres call recorded — the Krylov path did not run"
    bad_restart = [r for r, _ in seen if r is not None and r < composite]
    stalled = [info for _, info in seen if info not in (0, None)]
    assert not bad_restart, (
        f"gmres restart {min(bad_restart)} < composite {composite} — the bulk-sized "
        f"restart regressed (ERR-053 #282); GMRES will truncate the trace+seed DOFs"
    )
    assert not stalled, (
        f"{len(stalled)} within-group GMRES call(s) returned info!=0 (stall) — the "
        f"augmented composite was not fully spanned; keff={result.keff}"
    )
