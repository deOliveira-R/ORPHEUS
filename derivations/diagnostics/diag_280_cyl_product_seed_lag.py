"""Diagnostic: the product-cylinder (L+C).solve lags the Morel-Montry seed (#280 Phase 2.5b).

Created by numerics-investigator on 2026-07-05.

Feasibility investigation for #280 Phase 2.5b: the SN cylindrical within-group
loss operator ``(L+C)`` has an EXACTLY block-lower-triangular matvec (pinned by
``tests/sn/sweep/test_assembly_mode.py::test_282_augmented_walk_order_is_triangular``)
yet its cold ``.solve`` != ``(L+C)^-1`` for a PRODUCT quadrature.  This file
characterizes the root cause: the ONLY lagged element is the Morel-Montry
per-level starting seed, which for a product quadrature is the first-swept
ordinate's own flux ``psi_{m0}`` (edge_extrapolated_seed with t=0, #229), read
off the ITERATE.

Root cause CONFIRMED (this file):
  * cold solve != M^-1 (bulk err ~0.57), while seeding with the converged
    flux gives a SINGLE-PASS exact inverse -> the only lag is the seed;
  * the seed stencil is t=0 exactly (mu_start == mu_{m0}) and reads the
    SAME ordinate it feeds (first-swept) -> a per-ordinate self-reference.

If this test catches a real bug, promote to ``tests/sn/sweep/`` (the matching
per-module folder for the augmented sweep operator).  The
``test_cold_single_pass_equals_matvec_inverse`` xfail is the TARGET gate for
the Phase 2.5b fold: it flips GREEN when the direct single-pass reformulation
lands.
"""
import numpy as np
import pytest

from tests.sn.sweep.test_assembly_mode import _loss, _probe_augmented_matrix_one_group
from orpheus.geometry import BC, CoordSystem
from orpheus.geometry.mesh import Mesh1D
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.derivations.common.xs_library import get_mixture
from orpheus.transport.full_field import FullField
from orpheus.transport.source_sinks import AngularSourceSink, AngularBoundarySourceSink


def _cyl_product_mesh():
    return SNMesh(
        Mesh1D(edges=np.array([0.0, 0.3, 0.8, 1.0]), mat_ids=np.array([0, 1, 0]),
               bc_left=BC("reflective"), bc_right=BC("vacuum"),
               coord=CoordSystem.CYLINDRICAL),
        Quadrature.product(n_mu=4, n_phi=8),
        {0: get_mixture("A", "2g"), 1: get_mixture("B", "2g")},
    )


def _random_source(sn, N, ng):
    q = FullField(bulk=AngularSourceSink.zeros_on(sn),
                  boundary=AngularBoundarySourceSink.zeros_on(sn))
    q.bulk.values[:] = np.random.default_rng(5).random((N, ng, *sn.spatial_shape))
    return q


def test_cold_solve_lags_but_seeding_with_converged_is_single_pass():
    """Root cause: the ONLY lagged element is the seed.

    Cold ``.solve`` misses ``M^-1`` by O(0.5); seeding the SAME cold solve with
    the converged flux gives a single-pass exact inverse.  That isolates the
    lag to the seed (everything else in the sweep is feed-forward).
    """
    sn = _cyl_product_mesh()
    A = _loss(sn)
    g = 0
    N = sn.quad.n_ordinates
    M = _probe_augmented_matrix_one_group(sn, g)
    q = _random_source(sn, N, sn.ng)
    ref = np.linalg.solve(M, q.bulk.values[:, g].ravel())

    cold = np.asarray(A.solve(q).bulk.values)[:, g].ravel()
    # Converge by threading the iterate.
    ig = None
    for _ in range(60):
        psi = A.solve(q, initial_guess=ig)
        ig = psi
    one_pass_from_converged = np.asarray(
        A.solve(q, initial_guess=psi).bulk.values)[:, g].ravel()

    cold_err = np.abs(cold - ref).max()
    seeded_err = np.abs(one_pass_from_converged - ref).max()
    assert cold_err > 1e-2, (
        f"Diagnostic: cold solve error {cold_err:.3e} — expected the lag "
        "(>1e-2); if this is ~0 the product-cyl seed is no longer lagged "
        "(the Phase 2.5b fold has landed — retire this characterization)."
    )
    assert seeded_err < 1e-12, (
        f"Diagnostic: seeded-with-converged error {seeded_err:.3e} (expected "
        "<1e-12).  If this fails, the lag is NOT solely the seed — there is a "
        "second lagged element."
    )


def test_product_seed_is_first_swept_ordinate_flux_t_equals_zero():
    """The product-cyl seed stencil is t=0 exactly and reads the ordinate it feeds.

    ``edge_extrapolated_seed`` = (1-t)*psi[m0] + t*psi[m1].  For a product
    quadrature the level's starting direction coincides with the first-swept
    ordinate (mu_start == mu_{m0}, #229), so t=0 and the seed IS psi_{m0} —
    the very ordinate whose balance consumes the seed.  A clean self-reference.
    """
    sn = _cyl_product_mesh()
    closure = sn.pole_angular_closure
    for p, li in enumerate(sn.quad.level_indices):
        li = np.asarray(li)
        m0, m1, t = closure._edge_seed_stencil(p)
        # first-swept ordinate = level_indices[p][0] (the walk stacks in this order)
        assert t == pytest.approx(0.0, abs=1e-15), (
            f"Diagnostic: level {p} seed stencil t={t} (expected 0 for product)."
        )
        assert int(li[m0]) == int(li[0]), (
            f"Diagnostic: level {p} seed ordinate m0={int(li[m0])} != first-swept "
            f"{int(li[0])} — the seed must read the ordinate it feeds for the "
            "per-ordinate self-reference to be foldable."
        )


@pytest.mark.xfail(reason="Phase 2.5b fold not yet implemented: the cold product-cyl "
                          "solve lags the m0 seed. Flips GREEN when the direct "
                          "single-pass reformulation lands (#280 Phase 2.5b).",
                   strict=True)
def test_cold_single_pass_equals_matvec_inverse():
    """TARGET gate: the cold cyl solve reproduces the (triangular) matvec inverse
    in a single pass.  Fails today (seed lag); the Phase 2.5b fold makes it pass.
    """
    sn = _cyl_product_mesh()
    A = _loss(sn)
    N = sn.quad.n_ordinates
    q = _random_source(sn, N, sn.ng)
    for g in range(sn.ng):
        M = _probe_augmented_matrix_one_group(sn, g)
        ref = np.linalg.solve(M, q.bulk.values[:, g].ravel())
        cold = np.asarray(A.solve(q).bulk.values)[:, g].ravel()
        assert np.abs(cold - ref).max() < 1e-12
