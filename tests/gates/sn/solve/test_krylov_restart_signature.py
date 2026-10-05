r"""ERR-053 regression catcher — Krylov restart-truncation structural signature.

Promoted from ``derivations/diagnostics/diag_krylov_si_homogeneous_sphere_step5_mesh_scaling.py``
(numerics-investigator's step-5 diagnostic, retired at ``f36572c8``) per the investigator's
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
machine precision regardless of mesh.

Since #200 (2026-10-04) the within-group GMRES is preconditioned by the
sweep and needs 15 to 25 steps on this fixture at every mesh, so a
re-introduced ``restart=min(50, n_dof)`` clamp no longer bites here: the
mesh-refinement rows stay green under it (``[M]`` qa 2026-10-04,
``scratch/reference_architecture/p3/qa200/m1_restart50_full.log``; again
2026-10-05, ``scratch/reference_architecture/p3/gates200b/arm53_*.log``).
They keep the k_inf recovery claim and lost the ERR-053 one. ERR-053's
catchers are now the value row at the end of this file, one restart cycle
on a diffusive infinite medium whose honest solve needs more than 50
steps, and the restart site gates of
``test_krylov_curvilinear_precond_safety.py`` (``test_g_d3_3_*``).

References
==========

* ``docs/theory/verification/error_catalog.rst`` ERR-053 — full
  mechanism, root cause analysis, and fix rationale.
* The 8 diagnostic scripts of the full bisection cascade,
  ``diag_krylov_si_homogeneous_sphere_step{1..8}_*.py``, are no longer in
  the tree: steps 1, 2, 5 and 6 were retired at ``f36572c8`` and steps
  3, 4, 7 and 8 at ``d8843ba9``; recover one with
  ``git show <hash>^:derivations/diagnostics/<file>``.
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
    from orpheus.geometry import BC, StructuredGeometry
    from orpheus.mesh import CellsByCount, Mesher
    from orpheus.numerics.quadrature import Quadrature
    from orpheus.sn.solver import solve_sn

    fuel = get_mixture("A", "2g")
    geom = StructuredGeometry.sphere((0.0, 2.0), (0,), outer=BC.reflective)
    mesh = Mesher(geom).partition(CellsByCount.uniform_volume(n_cells)).mesh
    quad = Quadrature.gauss_legendre(n_ordinates=8)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        res = solve_sn(
            materials={0: fuel}, mesh=mesh, quadrature=quad,
            inner_solver=inner_solver,
            keff_tol=1e-12, flux_tol=1e-10,
            inner_tol=inner_tol,
        )
    return float(res.outcome.keff)


@pytest.mark.l1
@pytest.mark.verifies("sn-curvilinear-homogeneous-kinf-recovery")
@pytest.mark.parametrize("n_cells", [5, 8, 10, 16, 20, 30])
def test_krylov_kinf_independent_of_mesh_refinement(
    n_cells: int, _kinf_analytical: float,
) -> None:
    r"""Krylov ``keff`` on homogeneous reflective sphere == ``k_inf`` at every n_cells.

    Subspace-dimension defects produce DIVERGING error with mesh refinement
    (the natural subspace grows past a cap); tolerance defects produce
    uniform or monotone-decreasing error.  Pre-ERR-053-fix this test
    produced ``4.7e-01`` error at ``n_cells=20``, with the UNPRECONDITIONED
    GMRES of the time, whose step count grew with the mesh.

    No longer an ERR-053 catcher (marker removed 2026-10-05): the
    sweep-preconditioned solve needs 15 to 25 steps at every row, below
    any clamp at 50, so a re-introduced ``restart=min(50, n_dof)`` leaves
    all six rows green (``[M]``, the module docstring). What the rows
    still test is the k_inf recovery of the production Krylov path at
    every mesh, and a defect whose step count outgrows the restart at the
    finer rows. ERR-053's value catcher is
    :func:`test_one_full_restart_cycle_solves_a_diffusive_infinite_medium`.

    Note: parametrised on ``n_cells`` rather than tabulated inline so
    each refinement step is an independent test invocation in the
    V&V matrix (one failing cell does not mask the others).
    """
    keff = _solve_kinf(n_cells=n_cells, inner_solver="krylov")
    err = abs(keff - _kinf_analytical)
    # Gate: catches the ERR-053 subspace-truncation signature (pre-fix err was
    # 4.7e-1 — EIGHT orders above this gate).  Relaxed 1e-9 → 1e-7 on
    # 2026-06-05 for the then-unpreconditioned Krylov floor (~1.6e-9), and
    # re-tightened to 1e-9 with #200 (2026-10-04), as the relaxation asked.
    # ``[M]`` 2026-10-04 (scratch/reference_architecture/p3/krylov200/probes/
    # p3_restart_signature.log), |k − k_inf| at n_cells 5, 8, 10, 16, 20, 30:
    # sweep-preconditioned 8.4e-11, 2.2e-11, 2.4e-10, 9.9e-13, 6.8e-12, 5.1e-11;
    # the identity, on today's tree, 1.6e-14 to 8.3e-14. scipy's GMRES steers
    # its inner loop on the preconditioned residual, ‖M r‖ ≤ rtol·‖M b‖, and
    # accepts the solve on the TRUE residual, ‖b − A x‖ ≤ rtol·‖b‖ (``gmres``
    # in scipy/sparse/linalg/_isolve/iterative.py); the sweep-preconditioned
    # solve passes both at inner_tol = 1e-8 within 15 to 25 steps and stops
    # there, so it is not the machine-precision arm on this flat problem; the
    # 1e-9 band has 4x headroom over its worst row.
    assert err < 1e-9, (
        f"Krylov keff = {keff:.10f}, ref = {_kinf_analytical:.10f}, "
        f"err = {err:.3e}.  A Krylov k_inf error growing with refinement is "
        f"the subspace-truncation signature (ERR-053's class).  See "
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


#: The value row's diffusive slab: an infinite medium posed as a 100 cm slab (100 mean free paths) reflective on
#: both faces, one group, scattering ratio c = 0.9999, a flat isotropic source, 100 cells, Gauss-Legendre S8.
_DIFFUSIVE_C = 0.9999
_DIFFUSIVE_INNER_TOL = 1e-8
#: An iterative result is held to ten times the tolerance that stopped it.
_DIFFUSIVE_RTOL = 10.0 * _DIFFUSIVE_INNER_TOL
#: The clamp the defect re-introduces; the row is live only while the honest solve needs more Arnoldi steps.
_ERR053_CLAMP = 50


@pytest.mark.l1
@pytest.mark.catches("ERR-053")
@pytest.mark.rests_on(
    "tests/gates/sn/solve/test_krylov_sweep_preconditioner.py::"
    "test_p200_1_the_preconditioner_is_the_full_space_sweep_inverse[slab_reflective]",
    "tests/gates/sn/operators/test_sweep_inverse_identity.py::TestSweepInverseIdentity::"
    "test_forward_of_inverse_is_identity_on_a_random_composite[slab_reflective]",
)
def test_one_full_restart_cycle_solves_a_diffusive_infinite_medium() -> None:
    r"""ERR-053's VALUE gate under the sweep preconditioner: one GMRES restart cycle (``max_inner = 1``) solves the
    within-group system of a diffusive infinite medium to its closed form, :math:`\phi = W q / \Sigma_a` in every
    cell (``W`` the quadrature's weight sum), because the full-size restart spans the whole Krylov space.

    Since #200 the sweep-preconditioned solves of the rows above need 15 to 25 Arnoldi steps, so a
    ``restart = min(50, n_dof)`` clamp never bites there (``[M]`` qa 2026-10-04, and this session: the six k_inf rows
    and the consistency row stay green under it). Here the scattering ratio is near 1 on a slab 100 mean free paths
    thick, and the honest solve needs 56 steps (``[M]`` 2026-10-05,
    ``scratch/reference_architecture/p3/gates200b/probe_value53b.log``): relative error 4.9e-13 in one cycle. With
    the clamp re-dropped into ``_within_group_krylov`` the cycle stops at 50 steps, scipy returns ``info = 1`` and
    the flux is off by 1.3e-3, four orders above :data:`_DIFFUSIVE_RTOL`; with the discarded ``info`` added (the
    defect's second half) the GMRES warning is gone and the value leg still reds. Under the default cycle budget
    the clamp converges after 280 steps to 2.2e-9 (inside the tolerance), so the value moves only when the budget
    is counted in restart cycles that the clamp shortens, which is ERR-053's mechanism.

    One group and a flat solution null every spatial and group-coupling term; the row's threat is the Krylov
    truncation, which they do not touch. The infinite-medium answer is the structurally independent reference."""
    from orpheus.derivations.common.xs_library import make_mixture
    from orpheus.geometry import BC, StructuredGeometry
    from orpheus.mesh import CellsByCount, Mesher
    from orpheus.numerics.quadrature import Quadrature
    from orpheus.sn.solver import solve_sn_fixed_source

    sigma_t = 1.0
    sigma_a = sigma_t * (1.0 - _DIFFUSIVE_C)
    medium = make_mixture(
        sig_t=np.array([sigma_t]), sig_c=np.array([sigma_a]), sig_f=np.array([0.0]),
        nu=np.array([0.0]), chi=np.array([0.0]), sig_s=np.array([[sigma_t * _DIFFUSIVE_C]]),
    )
    n_cells = 100
    geometry = StructuredGeometry.slab((0.0, 100.0), (0,), left=BC.reflective, right=BC.reflective)
    mesh = Mesher(geometry).partition(CellsByCount.uniform_width(n_cells)).mesh
    quadrature = Quadrature.gauss_legendre(n_ordinates=8)
    source = 1.0
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        solution = solve_sn_fixed_source(
            materials={0: medium}, mesh=mesh, quadrature=quadrature,
            external_source=np.full((quadrature.N, 1, n_cells), source),
            inner_solver="krylov", inner_tol=_DIFFUSIVE_INNER_TOL, max_inner=1,
        )
    exact = float(np.sum(quadrature.weights)) * source / sigma_a
    phi = np.asarray(solution.scalar_flux.values, dtype=float)
    error = float(np.max(np.abs(phi - exact))) / exact
    assert error <= _DIFFUSIVE_RTOL, (
        f"one restart cycle left the infinite-medium flux off by {error:.2e} (relative), above "
        f"{_DIFFUSIVE_RTOL:.0e}: the restart no longer spans the Krylov space (ERR-053, a restart clamp)"
    )
    inner = [child for child in solution.record.children if child.label == "inner(gmres)"] or [solution.record]
    steps = max(record.n_iterations for record in inner)
    assert steps > _ERR053_CLAMP, (
        f"the honest solve needed only {steps} Arnoldi steps, not more than {_ERR053_CLAMP}: a clamp at "
        f"{_ERR053_CLAMP} would no longer bite and the row is blind"
    )
    convergence = [str(w.message)[:120] for w in caught if issubclass(w.category, RuntimeWarning)]
    assert not convergence, f"the one-cycle solve warned: {convergence}"
    assert solution.record.fully_converged, "the one-cycle solve does not read fully converged"
