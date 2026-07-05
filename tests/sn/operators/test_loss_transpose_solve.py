r"""The #280 Phase 2.5b reverse-scan ``(L+C)⁻ᵀ`` — G1 round-trip + G2 dense-Mᵀ.

``InvertibleOperator.solve_transpose`` is the transpose-solve
:math:`(L+C)^{-\mathsf T}` (the reverse-mode adjoint of the forward WDD
sweep-scan, sharing its ``ordinate_scan`` substrate via
:func:`~orpheus.sn.spatial.scan.ordinate_scan_transpose`).  These are the
object-level correctness gates:

* **G1 — round-trip.** ``solve_transpose`` inverts ``apply_transpose``
  (= the matvec transpose :math:`(L+C)^{\mathsf T}`, independently pinned
  by ``test_g_adjoint_reciprocity``): ``apply_transpose ∘ solve_transpose =
  I`` and its dual, on the **bulk** subspace where ``apply``/``solve`` are
  genuine inverses (the #284 source subspace — ``solve`` *computes* the
  boundary/seed OUTFLOW slots while ``apply`` treats them as free DOFs, so
  the two are inverse only off those slots).
* **G2 — dense-`Mᵀ` oracle (the structural keystone).** The augmented
  one-group matrix ``M`` is built by column-probing the FORWARD ``apply``
  (``_probe_augmented_matrix_one_group``); ``solve_transpose(b)`` reproduces
  ``np.linalg.solve(M.T, b)`` on the source-carried slots (bulk ⊕ seed
  cells ⊕ inflow corner).  Structurally independent of the reverse-walk
  code (the reference never touches it — L11 / the #284-discharge shape).
* **assembled-`Mᵀ`** (DD slab): the LAPACK upper-triangular back-substitution
  ``solve_triangular(permuted.T, …, lower=False)`` — a second
  structurally-independent realization, catches a wrong transposed-scan
  coefficient.

vv Mode-8: ``np.testing`` only (fires under ``-O``).  The full §11/§14
mutation suite (forward-DAG order, ±μ mirror, a'/b', the ×V/÷V two-denom
seam) + the cylinder leg land with the test-architect gate pass; this file
is the main-agent keystone the reverse-scan commit rides on.
"""
from __future__ import annotations

import numpy as np
import pytest
from scipy.linalg import solve_triangular

from orpheus.derivations.common.xs_library import get_mixture
from orpheus.geometry import BC, CoordSystem
from orpheus.geometry.mesh import Mesh1D
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.sn.operators.streaming import StreamingOperator
from orpheus.transport.fields.angular_boundary_flux import AngularBoundaryFlux
from orpheus.transport.fields.angular_flux import AngularFlux
from orpheus.transport.fields.starting_direction_flux import StartingDirectionFlux
from orpheus.transport.full_field import FullField
from orpheus.transport.operators.multiplication_operator import (
    MultiplicationOperator,
)

# The augmented forward-apply probe — the structurally-independent G2 oracle.
from tests.sn.sweep.test_assembly_mode import (
    _augmented_sweep_order,
    _probe_augmented_matrix_one_group,
)

_RTOL = 1e-10


def _slab() -> SNMesh:
    return SNMesh(
        Mesh1D(edges=np.array([0.0, 0.5, 1.5, 3.0, 5.0]),
               mat_ids=np.array([0, 1, 1, 0]),
               bc_left=BC("vacuum"), bc_right=BC("vacuum")),
        Quadrature.gauss_legendre(n_ordinates=4),
        {0: get_mixture("A", "2g"), 1: get_mixture("B", "2g")},
    )


def _sphere() -> SNMesh:
    return SNMesh(
        Mesh1D(edges=np.array([0.0, 0.3, 0.8, 1.0]),
               mat_ids=np.array([0, 1, 0]),
               bc_left=BC("reflective"), bc_right=BC("vacuum"),
               coord=CoordSystem.SPHERICAL),
        Quadrature.gauss_legendre(n_ordinates=4),
        {0: get_mixture("A", "2g"), 1: get_mixture("B", "2g")},
    )


_MESHES = {"slab": _slab, "sphere": _sphere}


def _loss(sn_mesh: SNMesh):
    mat_xs = sn_mesh.material_xs_field()
    return StreamingOperator(sn_mesh) + MultiplicationOperator(
        coefficient=mat_xs.total_cross_section_field,
        space=sn_mesh.full_field_space,
    )


def _fresh(sn_mesh: SNMesh) -> FullField:
    space = sn_mesh.starting_direction_space
    return FullField.zeros(
        bulk=AngularFlux, boundary=AngularBoundaryFlux, mesh=sn_mesh,
        starting_direction=None if space is None else StartingDirectionFlux,
    )


def _read_augmented(out, sn_mesh, g) -> np.ndarray:
    """The full augmented probe layout (incl. the outflow corner) — for the
    dense oracle's index bookkeeping."""
    bulk = np.asarray(out.bulk.values)[:, g].ravel()
    space = sn_mesh.starting_direction_space
    if space is None:
        return bulk
    seed = np.concatenate([
        np.concatenate([
            [space.corner_view(out.starting_direction.values, p, -1)[g]],
            space.cells_view(out.starting_direction.values, p, -1)[g][::-1],
            space.cells_view(out.starting_direction.values, p, +1)[g],
            [space.corner_view(out.starting_direction.values, p, +1)[g]],
        ]) for p in space.levels
    ])
    return np.concatenate([seed, bulk])


def _source_carried_mask(sn_mesh) -> np.ndarray:
    """Boolean mask (augmented layout) selecting source-carried slots — every
    slot EXCEPT each seed leg's trailing outflow corner."""
    space = sn_mesh.starting_direction_space
    N = sn_mesh.quad.n_ordinates
    nx = int(np.prod(sn_mesh.spatial_shape))
    if space is None:
        return np.ones(N * nx, dtype=bool)
    per = 2 * nx + 2  # seed leg: corner_in, cells⁻, cells⁺, corner_out
    mask = []
    for _ in space.levels:
        leg = np.ones(per, dtype=bool)
        leg[-1] = False                       # the outflow corner (free DOF)
        mask.append(leg)
    return np.concatenate([np.concatenate(mask), np.ones(N * nx, dtype=bool)])


# ── G1 — round-trip (bulk subspace) ────────────────────────────────────


@pytest.mark.foundation
@pytest.mark.parametrize("geom", list(_MESHES))
def test_g1_round_trip_bulk(geom):
    """``solve_transpose`` inverts ``apply_transpose`` on the bulk subspace
    (both directions)."""
    sn = _MESHES[geom]()
    A = _loss(sn)
    rng = np.random.default_rng(20260705)
    x = _fresh(sn)
    x.bulk.values[:] = rng.random(x.bulk.values.shape)
    if sn.starting_direction_space is not None:
        x.starting_direction.values[:] = rng.random(x.starting_direction.values.shape)
    back = A.solve_transpose(A.apply_transpose(x))
    np.testing.assert_allclose(
        np.asarray(back.bulk.values), np.asarray(x.bulk.values),
        rtol=_RTOL, atol=1e-11,
        err_msg=f"{geom}: solve_transpose∘apply_transpose ≠ I on the bulk",
    )
    b = _fresh(sn)                            # bulk-only source-subspace b
    b.bulk.values[:] = rng.random(b.bulk.values.shape)
    fwd = A.apply_transpose(A.solve_transpose(b))
    np.testing.assert_allclose(
        np.asarray(fwd.bulk.values), np.asarray(b.bulk.values),
        rtol=_RTOL, atol=1e-11,
        err_msg=f"{geom}: apply_transpose∘solve_transpose ≠ I on the bulk",
    )


# ── G2 — dense-Mᵀ oracle (source-carried slots) ────────────────────────


@pytest.mark.foundation
@pytest.mark.parametrize("geom", list(_MESHES))
def test_g2_dense_transpose_oracle(geom):
    """``solve_transpose(b)`` ≡ ``np.linalg.solve(Mᵀ, b)`` on the
    source-carried slots (M from the FORWARD apply — structurally
    independent of the reverse walk)."""
    sn = _MESHES[geom]()
    A = _loss(sn)
    rng = np.random.default_rng(20260706)
    mask = _source_carried_mask(sn)
    for g in range(sn.ng):
        M = _probe_augmented_matrix_one_group(sn, g)
        b = _fresh(sn)
        b.bulk.values[:, g] = rng.random((sn.quad.n_ordinates, *sn.spatial_shape))
        if sn.starting_direction_space is not None:
            b.starting_direction.values[:] = rng.random(
                b.starting_direction.values.shape
            )
        b_vec = _read_augmented(b, sn, g)
        got = _read_augmented(A.solve_transpose(b), sn, g)
        ref = np.linalg.solve(M.T, b_vec)
        np.testing.assert_allclose(
            got[mask], ref[mask], rtol=_RTOL, atol=1e-12,
            err_msg=f"{geom} g{g}: solve_transpose ≠ M⁻ᵀ b (source-carried)",
        )


# ── assembled-Mᵀ — the LAPACK back-substitution cross-check (slab) ──────


@pytest.mark.foundation
def test_assembled_transpose_lapack_slab():
    """DD slab: ``solve_transpose(b)`` ≡ upper-tri back-substitution of the
    walk-order-permuted ``Mᵀ`` (LAPACK ``dtrtrs`` — a second
    structurally-independent realization; catches a wrong a'/b')."""
    sn = _slab()
    A = _loss(sn)
    order = _augmented_sweep_order(sn)
    rng = np.random.default_rng(20260707)
    for g in range(sn.ng):
        M = _probe_augmented_matrix_one_group(sn, g)
        b = _fresh(sn)
        b.bulk.values[:, g] = rng.random((sn.quad.n_ordinates, *sn.spatial_shape))
        b_vec = np.asarray(b.bulk.values)[:, g].ravel()
        got = np.asarray(A.solve_transpose(b).bulk.values)[:, g].ravel()
        ref = np.empty_like(b_vec)
        ref[order] = solve_triangular(
            M[np.ix_(order, order)].T, b_vec[order], lower=False,
        )
        np.testing.assert_allclose(
            got, ref, rtol=_RTOL, atol=1e-12,
            err_msg=f"slab g{g}: solve_transpose ≠ LAPACK back-substitution",
        )
