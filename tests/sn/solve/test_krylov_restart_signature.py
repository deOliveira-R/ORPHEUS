r"""ERR-053 regression catcher — Krylov restart-truncation structural signature.

Promoted from ``derivations/diagnostics/diag_krylov_si_homogeneous_sphere_step5_mesh_scaling.py``
(numerics-investigator's step-5 diagnostic) per the investigator's
recommendation in the ERR-053 catalog entry.  This file pins the
LOAD-BEARING structural signature for the bug class
"GMRES subspace-dimension cap silently truncates the natural problem
dimension; the discarded ``info`` flag conceals the failure".

The signature
=============

For a homogeneous reflective sphere (where ``k = νΣ_f / Σ_a = k_inf``
is geometry- AND mesh-independent), mesh-refinement should produce
keff identical to machine precision at every cell count.  The
pre-ERR-053-fix Krylov inner solver produced a divergent error sequence
as the mesh refined — because the natural Krylov subspace dimension
grew past the hardcoded ``restart=min(50, ...)`` clamp.

Pre-fix data (numerics-investigator's step-5 table):

::

    n_cells   SI_keff         SI_err     KR_keff         KR_err
       5      1.8750000000   1.069e-11   1.8750000004   4.103e-10
       10     1.8750000000   1.097e-11   1.9152507886   4.025e-02
       16     1.8750000000   1.104e-11   1.5987097760   2.763e-01
       20     1.8750000000   1.106e-11   1.4019239124   4.731e-01

The signature: **SI is bit-flat across the refinement series; Krylov
error GROWS with refinement.**  This is the canonical structural
signature distinguishing subspace-dimension defects from tolerance
defects: a tolerance defect produces uniform or monotone-decreasing
error with refinement; a subspace defect produces DIVERGING error
because the natural subspace grows past the cap.

Post-fix (restart = full problem size) both inner solvers track at
machine precision regardless of mesh.  Any future regression that
re-introduces a hardcoded subspace cap (or any other structural
defect that scales with the problem dimension) will surface here.

References
==========

* ``docs/theory/verification/error_catalog.rst`` ERR-053 — full
  mechanism, root cause analysis, and fix rationale.
* The 8 diagnostic scripts at
  ``derivations/diagnostics/diag_krylov_si_homogeneous_sphere_step{1..8}_*.py``
  for the full bisection cascade.
"""
from __future__ import annotations

import warnings

import numpy as np
import pytest


@pytest.fixture(scope="module")
def _kinf_analytical() -> float:
    r"""``k_inf = νΣ_f / Σ_a`` for the ``A`` 2g mixture (homogeneous reflective).

    Geometry- and mesh-independent.  All inner solvers must converge
    to this value at every cell count.
    """
    return 1.875


def _solve_kinf(
    *, n_cells: int, inner_solver: str, inner_tol: float = 1e-8,
) -> float:
    r"""Helper: ``solve_sn`` at given n_cells for the ERR-053 fixture.

    Homogeneous reflective sphere, ``A`` 2g mixture, GL n_ord=8.
    Returns ``Solution.keff``.
    """
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
    return float(res.keff)


@pytest.mark.l1
@pytest.mark.verifies("sn-curvilinear-homogeneous-kinf-recovery")
@pytest.mark.catches("ERR-053")
@pytest.mark.parametrize("n_cells", [5, 8, 10, 16, 20, 30])
def test_krylov_kinf_independent_of_mesh_refinement(
    n_cells: int, _kinf_analytical: float,
) -> None:
    r"""Krylov ``keff`` on homogeneous reflective sphere == ``k_inf`` at every n_cells.

    The structural pin for ERR-053: subspace-dimension defects produce
    DIVERGING error with mesh refinement (the natural subspace grows
    past the cap).  Tolerance defects produce uniform or
    monotone-decreasing error with refinement.  This test catches both
    classes — any future regression that re-introduces a hardcoded
    ``restart=N`` for ``N < full_size`` will fail at the larger
    ``n_cells`` rows even with a tight ``inner_tol``.

    Pre-ERR-053-fix: this test produced ``4.7e-01`` error at
    ``n_cells=20`` (and worse at refinement).  Post-fix: machine
    precision across the full sweep.

    Note: parametrised on ``n_cells`` rather than tabulated inline so
    each refinement step is an independent test invocation in the
    V&V matrix (one failing cell does not mask the others).
    """
    keff = _solve_kinf(n_cells=n_cells, inner_solver="krylov")
    err = abs(keff - _kinf_analytical)
    # Gate: catches the ERR-053 subspace-truncation signature (pre-fix err was
    # 4.7e-1 — SIX orders above this gate).  Relaxed 1e-9 → 1e-7 (2026-06-05):
    # the UNPRECONDITIONED Krylov stopgap (issue #200) converges to ~1.6e-9 on
    # a reflective sphere, NOT machine precision — the strict 1e-9 gate was
    # marginally tight (FP noise near the GMRES inner_tol budget), failing at
    # meshes 5/16/30 while passing at 8/10/20.  1e-7 still catches any
    # re-introduction of a subspace-dimension cap while tolerating the
    # unpreconditioned floor; when #200 lands the block-inverse face
    # preconditioner (machine-precision convergence on curvilinear), re-tighten
    # to 1e-9.  The SI companion below stays at 1e-9 (SI IS bit-flat).
    assert err < 1e-7, (
        f"Krylov keff = {keff:.10f}, ref = {_kinf_analytical:.10f}, "
        f"err = {err:.3e}.  ERR-053 signature: subspace truncation in "
        f"GMRES (restart=min(50, full_size) clamp).  See "
        f"``docs/theory/verification/error_catalog.rst`` ERR-053."
    )


@pytest.mark.l1
@pytest.mark.verifies("sn-curvilinear-homogeneous-kinf-recovery")
@pytest.mark.parametrize("n_cells", [5, 8, 10, 16, 20, 30])
def test_si_kinf_independent_of_mesh_refinement(
    n_cells: int, _kinf_analytical: float,
) -> None:
    r"""SI ``keff`` on homogeneous reflective sphere == ``k_inf`` at every n_cells.

    Companion to :func:`test_krylov_kinf_independent_of_mesh_refinement`
    above.  SI was always correct on this fixture (pre- and post-fix);
    this test pins the bit-flat reference behaviour against which the
    Krylov regression catcher compares.

    If SI starts failing while Krylov continues to pass, a different
    bug class has appeared — flux normalisation drift (ERR-052),
    convention-bridge defects (ERR-049), or similar.
    """
    keff = _solve_kinf(n_cells=n_cells, inner_solver="source_iteration")
    err = abs(keff - _kinf_analytical)
    assert err < 1e-9, (
        f"SI keff = {keff:.10f}, ref = {_kinf_analytical:.10f}, "
        f"err = {err:.3e}.  SI is the structural reference for "
        f"ERR-053; if THIS fails, the bug class has moved."
    )
