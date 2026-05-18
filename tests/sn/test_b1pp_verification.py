r"""PR-TYPED-6.5 Phase 5 — B1'' face-state verification gates.

This module verifies the correctness of the B1'' face-state architecture
landed by PR-TYPED-6.5 Phase 3b on :class:`StreamingOperator` +
:class:`CollisionOperator` (the Resolution A leaves).

Why this lives alongside ``test_l1_standoff_slab_cylinder.py``
=============================================================

The existing L1 standoff test (``test_l1_standoff_slab_cylinder.py``)
exercises the cylinder + slab matvec-vs-sweep twin path through
``solve_sn → SNStreamingOperator``.  ``SNStreamingOperator`` is the
LEGACY bundle that retires at Step 7 of Issue #197 PR-TYPED-6c; until
then ``solve_sn`` keeps consuming it, and the Krylov leg there routes
through a patched ``apply`` that uses the cell-centre proxy at the
Carlson seed (the bug B1'' fixes).  That standoff therefore stays
``xfail strict`` on the twin-path and refinement tests until Step 7.

THIS module verifies B1'' directly on the new leaves
(``StreamingOperator + CollisionOperator``) without going through
``solve_sn`` or ``SNStreamingOperator`` at all.  GMRES on
``(L + C) ψ = q`` through the operator algebra is the direct
verification path the bug fix gates against.

After Step 7, ``solve_sn`` will consume the (L+C) algebra natively;
the existing L1 standoff's strict xfails flip green at that point,
and this file's tests subsume into the standoff naturally (or stay as
the focused B1''-only verification — TBD at retirement time).
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy.sparse.linalg import LinearOperator as SciLinearOperator
from scipy.sparse.linalg import gmres

from orpheus.geometry import BC, CoordSystem, Mesh1D
from orpheus.sn.geometry import SNMesh
from orpheus.sn.operator import (
    CollisionOperator,
    StreamingOperator,
    solution_to_angular_flux_with_traces,
)
from orpheus.sn.quadrature import GaussLegendre1D, LevelSymmetricSN
from tests.sn._test_helpers import placeholder_materials

# Per-test V&V level markers (see individual @pytest.mark.lN decorators).


# ═══════════════════════════════════════════════════════════════════════
# Fixtures
# ═══════════════════════════════════════════════════════════════════════


def _build_cylinder(nx: int, sn_order: int = 4) -> SNMesh:
    edges = np.linspace(0.01, 2.0, nx + 1)
    mesh = Mesh1D(
        edges=edges,
        mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.CYLINDRICAL,
        bc_left=BC("reflective"),
        bc_right=BC("reflective"),
    )
    quad = LevelSymmetricSN.create(sn_order=sn_order)
    return SNMesh(mesh, quad, placeholder_materials())


def _build_sphere(nx: int, n_ord: int = 8) -> SNMesh:
    edges = np.linspace(0.0, 2.0, nx + 1)
    mesh = Mesh1D(
        edges=edges,
        mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.SPHERICAL,
        bc_left=BC("reflective"),
        bc_right=BC("reflective"),
    )
    quad = GaussLegendre1D.create(n_ordinates=n_ord)
    return SNMesh(mesh, quad, placeholder_materials())


def _build_slab(nx: int, n_ord: int = 8) -> SNMesh:
    edges = np.linspace(0.0, 2.0, nx + 1)
    mesh = Mesh1D(
        edges=edges,
        mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.CARTESIAN,
        bc_left=BC("reflective"),
        bc_right=BC("reflective"),
    )
    quad = GaussLegendre1D.create(n_ordinates=n_ord)
    return SNMesh(mesh, quad, placeholder_materials())


_GEOMETRIES = [
    ("cylinder", _build_cylinder),
    ("sphere", _build_sphere),
    ("slab", _build_slab),
]


# ═══════════════════════════════════════════════════════════════════════
# L0 — sanity: (L + C) operator is full-rank and well-conditioned under B1''.
# ═══════════════════════════════════════════════════════════════════════
#
# The B1'' face-state packed vector adds genuine unknowns to the system
# (``ψ_face_outer`` for sphere/cylinder + slab, ``ψ_face_inner`` for
# slab and future hollow / annulus).  This test verifies the resulting
# (L + C) operator stays full-rank: the face residual equations
# ``WDD-propagated face − ψ_face = 0`` couple properly with the
# cell-balance equations.


@pytest.mark.l0
@pytest.mark.parametrize("name,builder", _GEOMETRIES)
def test_b1pp_lplusc_is_full_rank(name, builder):
    r"""(L + C) under B1'' is full-rank + well-conditioned on a small case.

    Builds (L + C) via :class:`StreamingOperator` + :class:`CollisionOperator`
    on a 5-cell homogeneous 1G case.  Materialises the dense matrix by
    probing with unit basis vectors; checks rank and condition number.
    """
    sn_mesh = builder(nx=5)
    ng = 1
    sigma_t = np.full((ng, sn_mesh.nx, 1), 0.4)
    L = StreamingOperator(sn_mesh, sigma_t)
    C = CollisionOperator(sn_mesh, sigma_t)
    n_unknowns = L.n_unknowns

    def matvec(psi: np.ndarray) -> np.ndarray:
        return L.apply(psi) + C.apply(psi)

    M = np.empty((n_unknowns, n_unknowns))
    for k in range(n_unknowns):
        e = np.zeros(n_unknowns)
        e[k] = 1.0
        M[:, k] = matvec(e)

    rank = np.linalg.matrix_rank(M)
    cond = np.linalg.cond(M)
    assert rank == n_unknowns, (
        f"{name}: (L + C) matrix is rank-deficient — "
        f"rank={rank}, n_unknowns={n_unknowns}"
    )
    assert cond < 1e8, (
        f"{name}: (L + C) condition number too large — cond={cond:.3e}"
    )


# ═══════════════════════════════════════════════════════════════════════
# L0 — constant flux is annihilated to σ_t · const under B1''.
# ═══════════════════════════════════════════════════════════════════════
#
# For a homogeneous reflective system, ψ = const at every cell-centre
# AND ψ_face_outer = const AT every outward face slot AND ψ_face_inner
# = const at every inner face slot (slab) is a meaningful state: every
# face equals every neighbour cell (no streaming gradient), and reflective
# BC maps const-on-outflow to const-on-inflow.  At this state:
#   - Streaming term vanishes (no gradient)
#   - Angular redistribution vanishes (homogeneous angular distribution)
#   - WDD face residual = (2·const − const) − const = 0 at every face
#   - Cell balance reduces to σ_t · const (collision only)


@pytest.mark.l0
@pytest.mark.parametrize("name,builder", _GEOMETRIES)
def test_b1pp_constant_flux_collapses_to_collision(name, builder):
    r"""For ψ = const at every B1'' slot, ``(L + C).apply(ψ) = σ_t·const``
    at cell slots and zero at face slots.

    Streaming cancels (WDD with flat ψ produces flat face flux);
    collision is the only surviving term in the cell-balance.  At face
    slots the WDD face residual is ``2·const − const_face_in − const_face = 0``.
    """
    sn_mesh = builder(nx=5)
    ng = 1
    sigma_t_val = 0.4
    sigma_t = np.full((ng, sn_mesh.nx, 1), sigma_t_val)
    L = StreamingOperator(sn_mesh, sigma_t)
    C = CollisionOperator(sn_mesh, sigma_t)
    n_unknowns = L.n_unknowns

    psi = np.ones(n_unknowns)
    out = L.apply(psi) + C.apply(psi)

    # Decode to inspect cell vs face slots.
    eq_map = L._ensure_eq_map(ng=ng)
    n_cell_scalars = eq_map.n_eq * ng

    # Cell slots: (L + C)·1 = σ_t · 1
    np.testing.assert_allclose(
        out[:n_cell_scalars], sigma_t_val, rtol=1e-12, atol=1e-13,
        err_msg=(
            f"{name}: (L + C) on flat ψ should reduce to σ_t at cell "
            f"slots; got out_cell with max diff "
            f"{np.max(np.abs(out[:n_cell_scalars] - sigma_t_val)):.3e}"
        ),
    )
    # Face slots: WDD residual = 0 at flat ψ.
    np.testing.assert_allclose(
        out[n_cell_scalars:], 0.0, rtol=1e-12, atol=1e-13,
        err_msg=(
            f"{name}: face residuals should be zero at flat ψ; "
            f"got max |face residual| = "
            f"{np.max(np.abs(out[n_cell_scalars:])):.3e}"
        ),
    )


# ═══════════════════════════════════════════════════════════════════════
# L1 — GMRES on (L + C) converges to FP-noise on the cylinder.
# ═══════════════════════════════════════════════════════════════════════
#
# Promotes the diagnostic at ``scratch/diag_b1pp_cylinder_gmres.py``
# into a permanent gate.  Pre-Phase-3b, the cylinder matvec carried an
# O(h) discretisation error (cell-centre-as-face proxy at the Carlson
# seed) that manifested as the cylinder twin-path divergence rel ≈ 4e-3.
# Post-Phase-3b, GMRES on (L + C) through B1'' converges at FP-noise
# (≈ 1e-13 rel residual) — the bug is fixed at its architectural source.


@pytest.mark.l1
@pytest.mark.parametrize("name,builder", _GEOMETRIES)
def test_b1pp_lplusc_gmres_converges_fp_noise(name, builder):
    r"""GMRES on ``(L + C) ψ = q`` through B1'' converges to FP-noise.

    Drives GMRES at ``rtol=1e-12``; verifies the relative residual
    ``||(L + C)ψ − q|| / ||q|| < 1e-10``.  Any residual above that
    bound signals the B1'' algebra is internally inconsistent.
    """
    sn_mesh = builder(nx=10)
    ng = 1
    sigma_t = np.full((ng, sn_mesh.nx, 1), 0.4)
    L = StreamingOperator(sn_mesh, sigma_t)
    C = CollisionOperator(sn_mesh, sigma_t)
    n_unknowns = L.n_unknowns

    def matvec(psi: np.ndarray) -> np.ndarray:
        return L.apply(psi) + C.apply(psi)

    LO = SciLinearOperator(
        shape=(n_unknowns, n_unknowns),
        matvec=matvec, dtype=np.float64,
    )

    rng = np.random.default_rng(seed=42)
    q = rng.standard_normal(n_unknowns)
    psi, info = gmres(
        LO, q,
        rtol=1e-12,
        maxiter=2000,
        restart=min(200, n_unknowns),
    )
    assert info == 0, (
        f"{name}: GMRES did not converge — info={info}, "
        f"n_unknowns={n_unknowns}"
    )
    r = matvec(psi) - q
    rel_residual = float(np.linalg.norm(r) / np.linalg.norm(q))
    assert rel_residual < 1e-10, (
        f"{name}: GMRES residual too large under B1'' — "
        f"rel_residual={rel_residual:.3e} (n_unknowns={n_unknowns})"
    )


# ═══════════════════════════════════════════════════════════════════════
# L0 — decoder round-trip is the identity (B1'' packed ↔ (cell, faces)).
# ═══════════════════════════════════════════════════════════════════════
#
# The B1'' packed vector layout is
#   ψ_packed = ψ_cell.ravel() ⊕ ψ_face_outer.ravel() ⊕ ψ_face_inner.ravel()
# Pure-reshape semantics demand decode ∘ encode == identity.  This is
# the Pattern-4 round-trip pin that catches index drift in the
# decoder / encoder pair.


@pytest.mark.l0
@pytest.mark.parametrize("name,builder", _GEOMETRIES)
def test_b1pp_decode_encode_roundtrip(name, builder):
    r"""Decoding a packed vector and re-packing the parts is the identity."""
    from orpheus.sn.operator import pack_with_traces

    sn_mesh = builder(nx=7)
    ng = 2
    sigma_t = np.full((ng, sn_mesh.nx, 1), 0.4)
    L = StreamingOperator(sn_mesh, sigma_t)
    eq_map = L._ensure_eq_map(ng=ng)

    rng = np.random.default_rng(seed=17)
    psi = rng.standard_normal(eq_map.n_unknowns)
    psi_cell, psi_face_outer, psi_face_inner = (
        solution_to_angular_flux_with_traces(
            psi, eq_map, sn_mesh.nx, ng, N=sn_mesh.quad.N,
        )
    )
    psi_back = pack_with_traces(
        psi_cell, psi_face_outer, psi_face_inner, eq_map,
    )
    np.testing.assert_array_equal(
        psi_back, psi,
        err_msg=f"{name}: pack ∘ decode is not the identity",
    )
