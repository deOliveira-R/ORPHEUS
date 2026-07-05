"""Diagnostic: the product-cyl m0-seed self-coupling folds into the cell diagonal (#280 Phase 2.5b).

Created by numerics-investigator on 2026-07-05.

Feasibility proof for #280 Phase 2.5b — is the lagged product-cylinder seed
retirable by a DIRECT single-pass reformulation?  YES.  The evidence:

  * the one-group augmented ``(L+C)`` matrix is EXACTLY walk-order triangular
    (single-pass forward substitution IS its inverse);
  * the seed ordinate m0's self-block is a triangular nx x nx block, and m0's
    output depends ONLY on itself (no coupling from any other ordinate) — so
    ``psi_{m0}`` is a LOCAL triangular solve of its own self-block against its
    own source;
  * the seed's contribution to that self-block is PURE DIAGONAL (off-diagonal
    change is exactly 0) — i.e. it folds into the m0 cell diagonal exactly as
    self-scatter does, with coefficient ``dA_w[m0]*c_in[m0]``;
  * injecting the locally-solved seed (NOT the iterate) makes the full solve a
    single pass equal to ``M^-1`` to machine precision.

These are GENERAL structural properties of the curvilinear augmented operator
(per L6/L12 they are the correct permanent gates), so promote to
``tests/sn/sweep/`` even after the fold lands: they pin the fold's correctness
and would red if a future change re-introduced a back edge or an off-diagonal
seed coupling.
"""
import numpy as np
import pytest

from tests.sn.sweep.test_assembly_mode import (
    _loss, _probe_augmented_matrix_one_group, _augmented_sweep_order)
from orpheus.geometry import BC, CoordSystem
from orpheus.geometry.mesh import Mesh1D
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.derivations.common.xs_library import get_mixture
from orpheus.sn.spatial.sweep_cache import GeometryCoefficients
from orpheus.sn.spatial.pole_angular_closure import MorelMontryAngularSweep
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


def _m0_globals(sn):
    closure = sn.pole_angular_closure
    out = []
    for p, li in enumerate(sn.quad.level_indices):
        m0, _, _ = closure._edge_seed_stencil(p)
        out.append(int(np.asarray(li)[m0]))
    return out


def test_m0_seed_ordinate_is_a_local_triangular_solve():
    """m0 depends ONLY on itself and its self-block is triangular.

    The seed ordinate's output row reads no other ordinate as input, and its
    nx x nx self-block is (walk-order) triangular — so ``psi_{m0}`` is a local
    single-pass forward substitution of its own source, independent of the
    global iteration.
    """
    sn = _cyl_product_mesh()
    N, nx = sn.quad.N, sn.nx
    for g in range(sn.ng):
        Mb = _probe_augmented_matrix_one_group(sn, g).reshape(N, nx, N, nx)
        for m0 in _m0_globals(sn):
            other = max(np.abs(Mb[m0, :, o, :]).max() for o in range(N) if o != m0)
            assert other < 1e-12, (
                f"Diagnostic g{g} m0={m0}: output couples to another ordinate "
                f"(max {other:.3e}) — m0 is not a pure starting ordinate."
            )
            selfblk = Mb[m0, :, m0, :]
            # inward ordinate (mu<0) marches outer->inner: natural-order upper-tri.
            lower = np.abs(np.tril(selfblk, k=-1)).max()
            upper = np.abs(np.triu(selfblk, k=1)).max()
            assert min(lower, upper) < 1e-12, (
                f"Diagnostic g{g} m0={m0}: self-block is not triangular "
                f"(tril={lower:.3e}, triu={upper:.3e})."
            )


def test_seed_contribution_is_pure_diagonal_fold():
    """The seed's contribution to the m0 self-block is PURE diagonal.

    Zeroing the seed in the matvec changes ONLY the m0 self-block diagonal
    (off-diagonal change is exactly 0) — the fold is a per-cell diagonal
    modification with coefficient ``dA_w[m0]*c_in[m0]`` (the m0 redistribution
    coefficient), exactly as self-scatter folds into the cell diagonal.
    """
    sn = _cyl_product_mesh()
    N, nx = sn.quad.N, sn.nx
    m0 = _m0_globals(sn)[0]
    Mb = _probe_augmented_matrix_one_group(sn, 0).reshape(N, nx, N, nx)
    orig = MorelMontryAngularSweep.edge_extrapolated_seed
    MorelMontryAngularSweep.edge_extrapolated_seed = (
        lambda self, psi_level, p: np.zeros(psi_level[:, 0, :].shape))
    try:
        Mb0 = _probe_augmented_matrix_one_group(sn, 0).reshape(N, nx, N, nx)
    finally:
        MorelMontryAngularSweep.edge_extrapolated_seed = orig

    diff = Mb[m0, :, m0, :] - Mb0[m0, :, m0, :]
    offdiag = np.abs(diff - np.diag(np.diag(diff))).max()
    diag = np.abs(np.diag(diff)).max()
    assert diag > 1e-2, (
        f"Diagnostic: seed contributes {diag:.3e} to the diagonal — expected "
        "O(1) (the seed is LIVE for product, NOT 'telescoped away')."
    )
    assert offdiag < 1e-12, (
        f"Diagnostic: seed contributes {offdiag:.3e} OFF-diagonal (expected 0) "
        "— the fold is not pure-diagonal; a self-block solve is needed instead."
    )


def test_local_seed_solve_gives_single_pass_matvec_inverse():
    """POC: seed = local triangular solve of m0's self-block -> single pass = M^-1.

    Solve each level's ``psi_{m0}`` directly from its own self-block (the fold),
    inject as the seed, and confirm the cold solve reproduces the matvec inverse
    in ONE pass.  This is the direct single-pass reformulation, proven feasible.
    """
    sn = _cyl_product_mesh()
    A = _loss(sn)
    N, nx, ng = sn.quad.N, sn.nx, sn.ng
    q = FullField(bulk=AngularSourceSink.zeros_on(sn),
                  boundary=AngularBoundarySourceSink.zeros_on(sn))
    q.bulk.values[:] = np.random.default_rng(5).random((N, ng, *sn.spatial_shape))
    m0_globals = _m0_globals(sn)

    # Build the per-level seed by local solve of each m0 self-block.
    seeds = {p: np.zeros((ng, nx)) for p in range(len(m0_globals))}
    for g in range(ng):
        Mb = _probe_augmented_matrix_one_group(sn, g).reshape(N, nx, N, nx)
        qg = q.bulk.values[:, g].reshape(N, nx)
        for p, m0 in enumerate(m0_globals):
            seeds[p][g] = np.linalg.solve(Mb[m0, :, m0, :], qg[m0])

    patched = lambda self, psi_level, p: seeds[p].copy()
    orig = MorelMontryAngularSweep.edge_extrapolated_seed
    MorelMontryAngularSweep.edge_extrapolated_seed = patched
    try:
        one_pass = A.solve(q)  # cold start; seed is the local fold, not the iterate
    finally:
        MorelMontryAngularSweep.edge_extrapolated_seed = orig

    for g in range(ng):
        M = _probe_augmented_matrix_one_group(sn, g)
        ref = np.linalg.solve(M, q.bulk.values[:, g].ravel())
        got = np.asarray(one_pass.bulk.values)[:, g].ravel()
        err = np.abs(got - ref).max()
        assert err < 1e-12, (
            f"Diagnostic g{g}: folded single-pass error {err:.3e} (expected "
            "<1e-12) — the local-seed fold does NOT reproduce M^-1."
        )


def test_level_symmetric_seed_is_dead_no_fold_needed():
    """The level-symmetric cyl seed is dead (c_in[m0]=0) — the fold is inert.

    LS quadratures give the first-swept ordinate ``c_in[m0]=0`` (raw tau=1 ->
    (1-tau)/tau=0), so there is no self-coupling to fold; the cold solve is
    already a single-pass exact inverse.  The Phase 2.5b fold must be a no-op here.
    """
    snL = SNMesh(
        Mesh1D(edges=np.array([0.0, 0.3, 0.8, 1.0]), mat_ids=np.array([0, 1, 0]),
               bc_left=BC("reflective"), bc_right=BC("vacuum"),
               coord=CoordSystem.CYLINDRICAL),
        Quadrature.level_symmetric(8),
        {0: get_mixture("A", "2g"), 1: get_mixture("B", "2g")},
    )
    geom = GeometryCoefficients.from_mesh_and_quad(snL)
    closure = snL.pole_angular_closure
    for p, li in enumerate(snL.quad.level_indices):
        m0, _, _ = closure._edge_seed_stencil(p)
        g_m0 = int(np.asarray(li)[m0])
        assert abs(geom.c_in[g_m0]) < 1e-14, (
            f"Diagnostic: LS level {p} m0={g_m0} c_in={geom.c_in[g_m0]:.3e} != 0 "
            "— the LS seed is not dead; the fold would NOT be a no-op."
        )
    # And the cold solve is already single-pass exact.
    A = _loss(snL)
    N = snL.quad.n_ordinates
    q = FullField(bulk=AngularSourceSink.zeros_on(snL),
                  boundary=AngularBoundarySourceSink.zeros_on(snL))
    q.bulk.values[:, 0] = np.random.default_rng(5).random((N, *snL.spatial_shape))
    M = _probe_augmented_matrix_one_group(snL, 0)
    ref = np.linalg.solve(M, q.bulk.values[:, 0].ravel())
    cold = np.asarray(A.solve(q).bulk.values)[:, 0].ravel()
    assert np.abs(cold - ref).max() < 1e-12
