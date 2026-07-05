r"""PR-TYPED-6.5 Phase 5 — B1'' face-state verification gates.

This module verifies the correctness of the B1'' face-state architecture
landed by PR-TYPED-6.5 Phase 3b on :class:`StreamingOperator` + the
collision multiplier ``C = M[σ_t]``
(:class:`~orpheus.transport.operators.multiplication_operator.MultiplicationOperator`,
the Resolution A leaves).

Why this lives alongside ``test_l1_standoff_slab_cylinder.py``
=============================================================

The existing L1 standoff test (``test_l1_standoff_slab_cylinder.py``)
exercises the cylinder + slab matvec-vs-sweep twin path through
``solve_sn`` (now routing the Krylov leg through
:class:`InvertibleOperator` post-D-K).  Pre-D-K, that test's Krylov
leg used a cell-centre proxy at the Carlson seed (the bug B1'' fixes)
and stayed ``xfail strict`` on the twin-path and refinement gates.

THIS module verifies B1'' directly on the operator leaves
(``StreamingOperator + MultiplicationOperator`` — the collision
multiplier ``C = M[σ_t]``).  GMRES on
``(L + C) ψ = q`` through the operator algebra is the direct
verification path the bug fix gates against.
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy.sparse.linalg import LinearOperator as SciLinearOperator
from scipy.sparse.linalg import gmres

from orpheus.geometry import BC, CoordSystem, Mesh1D
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.sn.operators.streaming import (
    StreamingOperator,
)
from orpheus.transport.operators.multiplication_operator import MultiplicationOperator
from orpheus.numerics.quadrature import Quadrature
from orpheus.transport.timed_full_field import TimedFullField
from tests.sn._test_helpers import (
    placeholder_materials,
    starting_direction_edge_seed,
)
from orpheus.transport.fields.angular_flux import AngularFlux
from orpheus.transport.fields.angular_boundary_flux import AngularBoundaryFlux
from orpheus.transport.fields.starting_direction_flux import StartingDirectionFlux

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
    quad = Quadrature.level_symmetric(sn_order=sn_order)
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
    quad = Quadrature.gauss_legendre(n_ordinates=n_ord)
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
    quad = Quadrature.gauss_legendre(n_ordinates=n_ord)
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

    Builds (L + C) via :class:`StreamingOperator` + the collision
    multiplier ``C = M[σ_t]``
    (:class:`~orpheus.transport.operators.multiplication_operator.MultiplicationOperator`)
    on a 5-cell homogeneous 1G case.  Materialises the dense matrix by
    probing with unit basis vectors; checks rank and condition number.

    D-I.3c (2026-05-29) — migrated from the bare-ndarray packed-vector
    contract to :class:`TimedFullField`.  The unit-basis probing now
    bridges through :meth:`TimedFullField.to_flat` /
    :meth:`TimedFullField.from_flat` (the direct-sum flat
    representation ``concat(bulk.ravel(), boundary.flat)``).  The
    full-rank / condition-number claim is layout-agnostic; the typed
    layout's flat dimension may differ from the legacy B1'' packed
    layout but the math (matrix is invertible) is the same.
    """
    sn_mesh = builder(nx=5)
    ng = 1
    sigma_t = np.full((ng, sn_mesh.nx), 0.4)
    L = StreamingOperator(sn_mesh)
    C = MultiplicationOperator.from_mesh(sigma_t, sn_mesh)

    template = TimedFullField.zeros(bulk=AngularFlux, boundary=AngularBoundaryFlux, mesh=sn_mesh, starting_direction=StartingDirectionFlux)
    n_flat = template.to_flat().size

    def matvec_flat(flat: np.ndarray) -> np.ndarray:
        psi_typed = TimedFullField.from_flat(flat, template)
        out_typed = L.apply(psi_typed) + C.apply(psi_typed)
        return out_typed.to_flat()

    M = np.empty((n_flat, n_flat))
    for k in range(n_flat):
        e = np.zeros(n_flat)
        e[k] = 1.0
        M[:, k] = matvec_flat(e)

    # The typed direct-sum flat layout includes AngularBoundaryFlux slots for
    # inflow ordinates (face_view slots where no equation is defined —
    # implicit-zero in the operator output).  Those rows AND columns
    # of M are structurally zero (no equation at that slot, and that
    # slot does not feed any equation through L or C).  Restrict the
    # rank / conditioning check to the equation-bearing subspace by
    # excluding all-zero rows AND columns.  The math claim
    # ("(L+C) is full-rank on the equation subspace") is layout-agnostic.
    nonzero_rows = np.any(M != 0, axis=1)
    nonzero_cols = np.any(M != 0, axis=0)
    nontrivial_idx = nonzero_rows & nonzero_cols
    M_sub = M[np.ix_(nontrivial_idx, nontrivial_idx)]
    n_eq = int(nontrivial_idx.sum())

    rank = np.linalg.matrix_rank(M_sub)
    cond = np.linalg.cond(M_sub)
    assert rank == n_eq, (
        f"{name}: (L + C) restricted matrix is rank-deficient — "
        f"rank={rank}, n_eq={n_eq} (full flat size {n_flat})"
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
    r"""For ψ = const at every slot, ``(L + C).apply(ψ)`` reduces to
    ``σ_t·const`` at cell slots; on the boundary block it emits the bare-matvec
    rows — the OUTFLOW self-consistency defect ``streamed − ψ.outflow`` (= 0 for
    flat ψ) on the outflow ordinate slots, and the INFLOW identity ``ψ.inflow``
    (= const) on the inflow ordinate slots.

    Streaming cancels (WDD with flat ψ produces flat face flux); collision is the
    only surviving term in the cell-balance.  **Wave O O.4b bare-matvec
    semantics:** the boundary block is NOT a single "WDD residual = 0" — the
    BC-extraction made ``(L + C)`` bare (it reads the given inflow, no in-sweep
    ``bc.apply``), so its boundary block carries the OUTFLOW defect (zero for
    flat ψ, since the WDD-propagated outflow equals ψ.outflow) AND the INFLOW
    identity ``ψ.inflow`` (so the composed ``(L + C − B)`` inflow residual is
    ``ψ.inflow − B·ψ.outflow``).

    D-I.3c (2026-05-29) — migrated from the bare-ndarray packed-vector contract
    to :class:`TimedFullField` (typed ``out.bulk`` / ``out.boundary`` split,
    Pattern 4: illegal-states-unrepresentable).  2026-06-05 — the boundary-block
    assertion migrated from the pre-extraction "WDD residual = 0 at every face"
    to the O.4b inflow-identity / outflow-defect split (the pre-extraction
    expectation predated O.4a.2/O.4b).
    """
    from dataclasses import replace

    sn_mesh = builder(nx=5)
    ng = 1
    sigma_t_val = 0.4
    sigma_t = np.full((ng, sn_mesh.nx), sigma_t_val)
    L = StreamingOperator(sn_mesh)
    C = MultiplicationOperator.from_mesh(sigma_t, sn_mesh)

    # Build flat-ψ TimedFullField: bulk = 1 everywhere AND boundary
    # face_view = 1 at every face slot (the "ψ = const at every B1''
    # slot" condition the docstring describes).  The face-flat buffer
    # is filled by assigning to every face_view in turn.
    state = TimedFullField.zeros(bulk=AngularFlux, boundary=AngularBoundaryFlux, mesh=sn_mesh, starting_direction=StartingDirectionFlux)
    bulk_values = np.ones_like(state.bulk.values)
    new_bulk = replace(state.bulk, values=bulk_values)
    new_boundary = state.boundary
    for face in new_boundary.layout.faces:
        new_boundary.face_view(face)[:] = 1.0
    # #282 route (a): the CONSISTENT ψ½ seed of a per-ordinate FLAT field is the
    # same constant (edge-extrap is constant-preserving), so the pole march
    # reproduces the flat field and (L+C)·const collapses to σ_t·const; a
    # present-ZERO seed would break the collapse.  None on non-carrying meshes.
    state = replace(
        state, bulk=new_bulk, boundary=new_boundary,
        starting_direction=starting_direction_edge_seed(bulk_values, sn_mesh),
    )

    out = L.apply(state) + C.apply(state)

    # Cell slots: (L + C)·1.bulk = σ_t · 1.bulk
    np.testing.assert_allclose(
        out.bulk.values, sigma_t_val, rtol=1e-12, atol=1e-13,
        err_msg=(
            f"{name}: (L + C) on flat ψ should reduce to σ_t at cell "
            f"slots; got out.bulk with max diff "
            f"{np.max(np.abs(out.bulk.values - sigma_t_val)):.3e}"
        ),
    )
    # Boundary block (Wave O O.4b bare-matvec semantics): NOT a single
    # "WDD residual = 0".  On each face the OUTFLOW ordinate slots carry the
    # self-consistency defect ``streamed − ψ.outflow`` (= 0 for flat ψ) and the
    # INFLOW ordinate slots carry the identity ``ψ.inflow`` (= 1 here).
    trace = sn_mesh.angular_trace
    for face in out.boundary.layout.faces:
        fv = out.boundary.face_view(face)
        outflow = trace.outflow_indices_for_face(face)
        inflow = trace.inflow_indices_for_face(face)
        np.testing.assert_allclose(
            fv[outflow], 0.0, rtol=1e-12, atol=1e-13,
            err_msg=(
                f"{name}: OUTFLOW defect (streamed − ψ.outflow) should be zero "
                f"at flat ψ on face {face!r}; got max "
                f"{np.max(np.abs(fv[outflow])):.3e}"
            ),
        )
        np.testing.assert_allclose(
            fv[inflow], 1.0, rtol=1e-12, atol=1e-13,
            err_msg=(
                f"{name}: INFLOW slots should carry the identity ψ.inflow=1 "
                f"(the bare-matvec inflow row) at flat ψ on face {face!r}; "
                f"got {fv[inflow]}"
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

    D-I.3c (2026-05-29) — migrated from the bare-ndarray packed-vector
    contract to :class:`TimedFullField` via the direct-sum flat-bridge
    (Pattern 7 at the producer): the matvec closure lifts the GMRES
    iterate ``flat → TimedFullField`` via :meth:`from_flat`, applies
    ``(L + C)`` in the typed algebra, and packs back ``TimedFullField
    → flat`` via :meth:`to_flat` at the scipy interface boundary.
    GMRES sees a pure flat-vector contract; the producer (the
    operator algebra) sees the typed contract; the bridge collapses
    at the producer side, never at the consumer side.
    """
    sn_mesh = builder(nx=10)
    ng = 1
    sigma_t = np.full((ng, sn_mesh.nx), 0.4)
    L = StreamingOperator(sn_mesh)
    C = MultiplicationOperator.from_mesh(sigma_t, sn_mesh)

    template = TimedFullField.zeros(bulk=AngularFlux, boundary=AngularBoundaryFlux, mesh=sn_mesh, starting_direction=StartingDirectionFlux)
    n_flat = template.to_flat().size

    def matvec(flat: np.ndarray) -> np.ndarray:
        psi_typed = TimedFullField.from_flat(flat, template)
        out_typed = L.apply(psi_typed) + C.apply(psi_typed)
        return out_typed.to_flat()

    # The typed flat layout has zero-row/zero-column slots for inflow
    # ordinates on each AngularBoundaryFlux face (no-equation slots).  Identify
    # the equation-bearing subspace by probing every column and finding
    # which rows AND columns are non-trivial (Pattern 7 at the producer:
    # the equation subspace is determined ONCE up-front, then GMRES
    # operates on the well-posed system).  Without this restriction the
    # operator is singular on the inflow-slot null space and GMRES does
    # not converge.
    M_probe = np.empty((n_flat, n_flat))
    for k in range(n_flat):
        e = np.zeros(n_flat)
        e[k] = 1.0
        M_probe[:, k] = matvec(e)
    nonzero_rows = np.any(M_probe != 0, axis=1)
    nonzero_cols = np.any(M_probe != 0, axis=0)
    eq_idx = np.where(nonzero_rows & nonzero_cols)[0]
    n_eq = int(eq_idx.size)

    def matvec_eq(x_eq: np.ndarray) -> np.ndarray:
        flat = np.zeros(n_flat)
        flat[eq_idx] = x_eq
        return matvec(flat)[eq_idx]

    LO = SciLinearOperator(
        shape=(n_eq, n_eq),
        matvec=matvec_eq, dtype=np.float64,
    )

    rng = np.random.default_rng(seed=42)
    q = rng.standard_normal(n_eq)
    psi, info = gmres(
        LO, q,
        rtol=1e-12,
        maxiter=2000,
        restart=min(200, n_eq),
    )
    assert info == 0, (
        f"{name}: GMRES did not converge — info={info}, "
        f"n_eq={n_eq} (n_flat={n_flat})"
    )
    r = matvec_eq(psi) - q
    rel_residual = float(np.linalg.norm(r) / np.linalg.norm(q))
    assert rel_residual < 1e-10, (
        f"{name}: GMRES residual too large under B1'' — "
        f"rel_residual={rel_residual:.3e} (n_eq={n_eq})"
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


# D-J (2026-05-30): ``test_b1pp_decode_encode_roundtrip`` retired with
# the :func:`solution_to_angular_flux_with_traces` /
# :func:`pack_with_traces` codec family.  The legacy bare-ndarray
# packed-vector layout no longer exists; the codec roundtrip has
# nothing to invert.  The math claim it pinned (decoder ∘ encoder ==
# identity on the packed layout) was a property of the retired layout.
