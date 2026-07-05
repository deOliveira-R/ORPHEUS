r"""The #280 Phase 2.5b direct-seed fold — the product-cylinder ``(L+C).solve``
is a single-pass direct inverse.

Before the fold the curvilinear within-group loss operator ``(L+C)`` had an
*exactly* block-lower-triangular matvec (pinned by
``test_assembly_mode.test_282_augmented_walk_order_is_triangular``) yet its cold
``.solve`` did **not** equal ``(L+C)^{-1}`` for a PRODUCT quadrature: the sole
lagged element was the Morel–Montry per-level starting seed ``ψ½``, which for a
product rule is the first-swept ordinate's OWN average ``ψ̄_{m0}`` (``t = 0``,
#229) read off the ITERATE.

The fold (``loss_representation/__init__.py`` ``_OneDimScanWalk._run``) folds that
self-coupling ``κ = (ΔA/w)·c_in`` into m0's cell diagonal — ``c_out → c_out −
c_in`` — exactly as self-scatter folds, and drops the seed source.  The cold
solve is then a single pass ``= M^{-1}`` for every geometry (the cylinder
analogue of the sphere's route-(a) starting-direction solve, #282 / 2.5d).

These gates pin, from the precondition down to the production result:

* the seed stencil IS ``t = 0`` and reads the ordinate it feeds (the
  self-reference that makes the fold pure-diagonal);
* m0's augmented self-block is a LOCAL triangular solve (couples to no other
  ordinate);
* the seed's matvec contribution is PURELY the m0 diagonal, coefficient
  ``(ΔA/w)·c_in`` (the fold coefficient — pins its sign/magnitude via a
  zero-the-seed mutation of the matvec);
* **the production cold solve reproduces the dense matvec inverse in one pass**
  (the headline gate — structurally-independent oracle
  ``_probe_augmented_matrix_one_group`` builds ``M`` from the FORWARD apply);
* the level-symmetric cylinder seed is dead (``c_in = 0``) so the fold is an
  exact no-op there.

vv Mode-8: ``np.testing`` / bare ``assert`` mixed — the module runs under ``-O``
(canonical), so the value gates use ``np.testing`` (fire under ``-O``); the
structural asserts guard shapes only.  The full §11/§14 fold mutation matrix
(``c_out+c_in`` wrong-sign, a'/b' folded coefficient, m0-not-first) lands with
the test-architect gate pass (#280 2.5b task #29).
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.common.xs_library import get_mixture
from orpheus.geometry import BC, CoordSystem
from orpheus.geometry.mesh import Mesh1D
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.sn.spatial.sweep_cache import GeometryCoefficients
from orpheus.transport.full_field import FullField
from orpheus.transport.source_sinks import (
    AngularBoundarySourceSink,
    AngularSourceSink,
)

# The augmented FORWARD-apply probe — the structurally-independent oracle.
from tests.sn.sweep.test_assembly_mode import (
    _loss,
    _probe_augmented_matrix_one_group,
)


def _cyl_product_mesh() -> SNMesh:
    return SNMesh(
        Mesh1D(edges=np.array([0.0, 0.3, 0.8, 1.0]), mat_ids=np.array([0, 1, 0]),
               bc_left=BC("reflective"), bc_right=BC("vacuum"),
               coord=CoordSystem.CYLINDRICAL),
        Quadrature.product(n_mu=4, n_phi=8),
        {0: get_mixture("A", "2g"), 1: get_mixture("B", "2g")},
    )


def _cyl_level_symmetric_mesh() -> SNMesh:
    return SNMesh(
        Mesh1D(edges=np.array([0.0, 0.3, 0.8, 1.0]), mat_ids=np.array([0, 1, 0]),
               bc_left=BC("reflective"), bc_right=BC("vacuum"),
               coord=CoordSystem.CYLINDRICAL),
        Quadrature.level_symmetric(8),
        {0: get_mixture("A", "2g"), 1: get_mixture("B", "2g")},
    )


def _random_source(sn: SNMesh) -> FullField:
    q = FullField(bulk=AngularSourceSink.zeros_on(sn),
                  boundary=AngularBoundarySourceSink.zeros_on(sn))
    q.bulk.values[:] = np.random.default_rng(5).random(
        (sn.quad.n_ordinates, sn.ng, *sn.spatial_shape)
    )
    return q


def _m0_globals(sn: SNMesh) -> list[int]:
    closure = sn.pole_angular_closure
    out = []
    for p, li in enumerate(sn.quad.level_indices):
        m0, _m1, _t = closure._edge_seed_stencil(p)
        out.append(int(np.asarray(li)[m0]))
    return out


# ── precondition: the product seed is a t=0 self-reference ─────────────


@pytest.mark.foundation
def test_product_seed_stencil_is_t_zero_and_first_swept():
    """``edge_extrapolated_seed`` = (1−t)ψ[m0] + t·ψ[m1]; product ⇒ t=0, m0 first.

    The level's starting direction coincides with the first-swept ordinate
    (μ_start == μ_{m0}, #229), so t=0 and the seed IS ψ_{m0} — the ordinate whose
    balance consumes the seed.  This self-reference is what makes the fold
    pure-diagonal; if a quadrature change broke t=0 the fold's ``_run`` guard
    would refuse the level (NotImplementedError) and this gate localizes why.
    """
    sn = _cyl_product_mesh()
    closure = sn.pole_angular_closure
    for p, li in enumerate(sn.quad.level_indices):
        li = np.asarray(li)
        m0, _m1, t = closure._edge_seed_stencil(p)
        # np.testing / pytest.fail — fire under -O (bare assert is stripped).
        np.testing.assert_allclose(
            t, 0.0, atol=1e-15,
            err_msg=f"level {p} seed stencil t={t} (expected 0 for product).",
        )
        if int(li[m0]) != int(li[0]):
            pytest.fail(
                f"level {p} seed ordinate m0={int(li[m0])} != first-swept "
                f"{int(li[0])} — the seed must read the ordinate it feeds."
            )


# ── structural: m0's self-block is a local triangular solve ────────────


@pytest.mark.foundation
def test_m0_self_block_is_local_triangular_solve():
    """m0 depends ONLY on itself and its nx×nx self-block is triangular.

    The seed ordinate's augmented-matrix row reads no other ordinate as input
    (other-ordinate coupling = 0), and its self-block is walk-order triangular —
    so ψ_{m0} is a local single-pass forward substitution of its own source,
    independent of any global iteration.  This is the structural licence for the
    in-place diagonal fold.
    """
    sn = _cyl_product_mesh()
    N, nx = sn.quad.N, sn.nx
    for g in range(sn.ng):
        Mb = _probe_augmented_matrix_one_group(sn, g).reshape(N, nx, N, nx)
        for m0 in _m0_globals(sn):
            other = max(np.abs(Mb[m0, :, o, :]).max() for o in range(N) if o != m0)
            np.testing.assert_array_less(
                other, 1e-12,
                err_msg=f"g{g} m0={m0}: output couples to another ordinate "
                        f"(max {other:.3e}) — m0 is not a pure starting ordinate.",
            )
            selfblk = Mb[m0, :, m0, :]
            lower = np.abs(np.tril(selfblk, k=-1)).max()
            upper = np.abs(np.triu(selfblk, k=1)).max()
            np.testing.assert_array_less(
                min(lower, upper), 1e-12,
                err_msg=f"g{g} m0={m0}: self-block is not triangular "
                        f"(tril={lower:.3e}, triu={upper:.3e}).",
            )


@pytest.mark.foundation
def test_seed_contribution_is_pure_diagonal_fold():
    """The seed's matvec contribution is PURE diagonal, coefficient (ΔA/w)·c_in.

    Zeroing the seed in the FORWARD matvec changes ONLY the m0 self-block
    diagonal (off-diagonal change is exactly 0) — so the fold is a per-cell
    diagonal modification ``κ = dA_w[m0]·c_in[m0]``, exactly as self-scatter
    folds.  This pins the fold coefficient the production ``_run`` subtracts
    (``c_out → c_out − c_in``): a wrong coefficient would leave a nonzero
    off-diagonal or a mismatched diagonal here.
    """
    from orpheus.sn.spatial.pole_angular_closure import MorelMontryAngularSweep

    sn = _cyl_product_mesh()
    N, nx = sn.quad.N, sn.nx
    m0 = _m0_globals(sn)[0]
    Mb = _probe_augmented_matrix_one_group(sn, 0).reshape(N, nx, N, nx)
    orig = MorelMontryAngularSweep.edge_extrapolated_seed
    MorelMontryAngularSweep.edge_extrapolated_seed = (  # type: ignore[method-assign]
        lambda self, psi_level, p: np.zeros(psi_level[:, 0, :].shape))
    try:
        Mb0 = _probe_augmented_matrix_one_group(sn, 0).reshape(N, nx, N, nx)
    finally:
        MorelMontryAngularSweep.edge_extrapolated_seed = orig  # type: ignore[method-assign]

    diff = Mb[m0, :, m0, :] - Mb0[m0, :, m0, :]
    offdiag = np.abs(diff - np.diag(np.diag(diff))).max()
    diag = np.abs(np.diag(diff)).max()
    np.testing.assert_array_less(
        1e-2, diag,
        err_msg=f"seed contributes {diag:.3e} to the diagonal — expected O(1) "
                "(the seed is LIVE for product, NOT 'telescoped away').",
    )
    np.testing.assert_array_less(
        offdiag, 1e-12,
        err_msg=f"seed contributes {offdiag:.3e} OFF-diagonal (expected 0) — "
                "the fold is not pure-diagonal; a self-block solve is needed.",
    )


# ── the headline: cold solve ≡ dense matvec inverse in one pass ────────


@pytest.mark.foundation
def test_cold_solve_equals_matvec_inverse_single_pass():
    """THE fold gate: the COLD product-cyl solve reproduces ``M^{-1}`` in one pass.

    ``M`` is built by column-probing the FORWARD ``apply``
    (``_probe_augmented_matrix_one_group``) — structurally independent of the
    solve's reverse-substitution code.  Before the fold this missed by O(0.5)
    (the seed lag); the direct-seed fold makes ``A.solve(q)`` (cold, no
    ``initial_guess``) equal ``np.linalg.solve(M, q)`` to machine precision.
    """
    sn = _cyl_product_mesh()
    A = _loss(sn)
    q = _random_source(sn)
    for g in range(sn.ng):
        M = _probe_augmented_matrix_one_group(sn, g)
        ref = np.linalg.solve(M, q.bulk.values[:, g].ravel())
        cold = np.asarray(A.solve(q).bulk.values)[:, g].ravel()
        np.testing.assert_allclose(
            cold, ref, rtol=1e-10, atol=1e-12,
            err_msg=f"g{g}: cold cyl solve != M^-1 single-pass (the fold lag).",
        )


@pytest.mark.foundation
def test_cold_solve_equals_si_converged_fixed_point():
    """Baseline-neutral: the fold's cold solve equals the SI-converged solve.

    The fold changes only HOW ψ_{m0} is obtained (iterate-lag → in-place direct
    solve), NOT the fixed point.  Threading the iterate to convergence reaches
    the SAME flux the cold single-pass now produces — so keff / MMS / converged
    snapshots do not move (the re-baseline is machine-identical, unlike the
    sphere's route-(a) closure change which moved keff(N) at fixed N).
    """
    sn = _cyl_product_mesh()
    A = _loss(sn)
    q = _random_source(sn)
    psi = A.solve(q)                       # cold single-pass
    cold = np.asarray(psi.bulk.values)
    for _ in range(40):                    # thread the iterate to convergence
        psi = A.solve(q, initial_guess=psi)
    converged = np.asarray(psi.bulk.values)
    np.testing.assert_allclose(
        cold, converged, rtol=1e-10, atol=1e-12,
        err_msg="cold single-pass solve != SI-converged fixed point.",
    )


# ── level-symmetric: the fold is an exact no-op (c_in = 0) ─────────────


@pytest.mark.foundation
def test_level_symmetric_cyl_seed_is_dead_fold_is_noop():
    """The level-symmetric cyl seed is dead (c_in[m0]=0) — the fold is inert.

    LS quadratures give the first-swept ordinate ``c_in[m0]=0`` (raw τ=1 →
    (1−τ)/τ=0), so ``c_out − c_in = c_out`` and there is no self-coupling to
    fold; the cold solve is already a single-pass exact inverse.
    """
    sn = _cyl_level_symmetric_mesh()
    geom = GeometryCoefficients.from_mesh_and_quad(sn)
    for m0 in _m0_globals(sn):
        np.testing.assert_array_less(
            abs(geom.c_in[m0]), 1e-14,
            err_msg=f"LS m0={m0} c_in={geom.c_in[m0]:.3e} != 0 — the LS seed is "
                    "not dead; the fold would not be a no-op.",
        )
    A = _loss(sn)
    q = _random_source(sn)
    for g in range(sn.ng):
        M = _probe_augmented_matrix_one_group(sn, g)
        ref = np.linalg.solve(M, q.bulk.values[:, g].ravel())
        cold = np.asarray(A.solve(q).bulk.values)[:, g].ravel()
        np.testing.assert_allclose(
            cold, ref, rtol=1e-10, atol=1e-12,
            err_msg=f"LS g{g}: cold cyl solve != M^-1 (should already be exact).",
        )
