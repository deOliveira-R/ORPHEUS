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
from orpheus.transport.radial_characteristic_field import RadialCharacteristicField
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


def _cyl_product() -> SNMesh:
    # The MANDATORY cylinder config (#280 2.5b cyl reverse-scan): a PRODUCT
    # quadrature has 8 pure-azimuthal DEGENERATE ordinates AND a LIVE seed-fold
    # (t=0, c_in≠0, mm_a_in[m0]=1) — the two cylinder-specific transpose terms
    # that ``level_symmetric`` simultaneously NULLS.
    return SNMesh(
        Mesh1D(edges=np.array([0.0, 0.3, 0.8, 1.0]),
               mat_ids=np.array([0, 1, 0]),
               bc_left=BC("reflective"), bc_right=BC("vacuum"),
               coord=CoordSystem.CYLINDRICAL),
        Quadrature.product(n_mu=4, n_phi=8),
        {0: get_mixture("A", "2g"), 1: get_mixture("B", "2g")},
    )


def _cyl_ls() -> SNMesh:
    # The everything-nulled CONTROL: level_symmetric cyl has 0 degenerate ords
    # and a DEAD seed-fold (c_in[m0]=0, mm_a_in[m0]=0) — it exercises only the
    # multi-level bulk M-M thread transpose.
    return SNMesh(
        Mesh1D(edges=np.array([0.0, 0.3, 0.8, 1.0]),
               mat_ids=np.array([0, 1, 0]),
               bc_left=BC("reflective"), bc_right=BC("vacuum"),
               coord=CoordSystem.CYLINDRICAL),
        Quadrature.level_symmetric(8),
        {0: get_mixture("A", "2g"), 1: get_mixture("B", "2g")},
    )


_MESHES = {
    "slab": _slab, "sphere": _sphere,
    "cyl_product": _cyl_product, "cyl_ls": _cyl_ls,
}


def _loss(sn_mesh: SNMesh):
    mat_xs = sn_mesh.material_xs_field()
    return StreamingOperator(sn_mesh) + MultiplicationOperator(
        coefficient=mat_xs.total_cross_section_field,
        space=sn_mesh.full_field_space,
    )


def _fresh(sn_mesh: SNMesh) -> FullField:
    return FullField.zeros(
        interior=AngularFlux, boundary=AngularBoundaryFlux, mesh=sn_mesh,
    )


def _seed_cot(sn_mesh: SNMesh, values=None):
    """The ray cotangent member for the transposed JOINT surfaces on a
    carrying mesh (step 6: the joint spelling is the grid — the cotangent
    is a CoupledField member, never a leg kwarg) — ``None`` seedless.

    SOURCE-role (the transposed substitution's member algebra is
    role-honest — the flux-role zeros the fused legs tolerated were a
    role leak; see ``test_inverse_adjoint_coherence._coupled_bulk_b``).
    ``values`` is a flat array of ``radial_characteristic_field_space``
    size (the split composite's ``to_flat`` layout: interior ⊕ boundary)."""
    if sn_mesh.radial_characteristic_field_space is None:
        return None
    cot = RadialCharacteristicField.source_zeros_on(sn_mesh)
    if values is not None:
        cot = RadialCharacteristicField.from_flat(
            np.asarray(values, dtype=float), cot,
        )
    return cot


def _fresh_source(sn_mesh: SNMesh) -> FullField:
    """A zero SOURCE-role composite — the codomain-side cotangent carrier
    the grid's transposed surfaces consume (role-honest member algebra)."""
    from orpheus.transport.source_sinks import (
        AngularBoundarySourceSink,
        AngularSourceSink,
    )

    return FullField(
        interior=AngularSourceSink.from_mesh(
            np.zeros((sn_mesh.quad.N, sn_mesh.ng, *sn_mesh.spatial_shape)),
            sn_mesh,
        ),
        boundary=AngularBoundarySourceSink.zeros_on(sn_mesh),
    )


def _read_augmented(out, sn_mesh, g, ray=None) -> np.ndarray:
    """The full augmented probe layout (incl. the outflow corner) — for the
    dense oracle's index bookkeeping. ``ray`` is the EXPLICIT ray-leg composite
    (B.2d / 4e — System B's split ``interior ⊕ boundary`` member)."""
    bulk = np.asarray(out.interior.values)[:, g].ravel()
    if sn_mesh.radial_characteristic_field_space is None or ray is None:
        return bulk
    seed = np.concatenate([
        np.concatenate([
            [ray.boundary.corner(p, -1)[g]],
            ray.interior.cells(p, -1)[g][::-1],
            ray.interior.cells(p, +1)[g],
            [ray.boundary.corner(p, +1)[g]],
        ]) for p in sn_mesh.radial_characteristic_levels
    ])
    return np.concatenate([seed, bulk])


def _source_carried_mask(sn_mesh) -> np.ndarray:
    """Boolean mask (augmented layout) selecting source-carried slots — every
    slot EXCEPT each seed leg's trailing outflow corner."""
    N = sn_mesh.quad.n_ordinates
    nx = int(np.prod(sn_mesh.spatial_shape))
    if sn_mesh.radial_characteristic_field_space is None:
        return np.ones(N * nx, dtype=bool)
    per = 2 * nx + 2  # seed leg: corner_in, cells⁻, cells⁺, corner_out
    mask = []
    for _ in sn_mesh.radial_characteristic_levels:
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
    carrying = sn.radial_characteristic_field_space is not None
    if carrying:
        # step 6: the joint Mᵀ round-trip is THE GRID's (presence
        # structural — the walk carries only the decoupled block);
        # cotangent members are SOURCE-role (role-honest algebra).
        from orpheus.numerics.coupled_system import CoupledField

        from tests.sn._test_helpers import joint_m_grid

        grid, _space = joint_m_grid(sn, A)
        x = _fresh_source(sn)
        x.interior.values[:] = rng.random(x.interior.values.shape)
        cot_in = _seed_cot(
            sn, rng.random(sn.radial_characteristic_field_space.shape[0]),
        )
        back = grid.solve_transpose(
            grid.apply_transpose(CoupledField(systems=(x, cot_in))),
        ).systems[0]
    else:
        x = _fresh(sn)
        x.interior.values[:] = rng.random(x.interior.values.shape)
        back = A.solve_transpose(A.apply_transpose(x))
    np.testing.assert_allclose(
        np.asarray(back.interior.values), np.asarray(x.interior.values),
        rtol=_RTOL, atol=1e-11,
        err_msg=f"{geom}: solve_transpose∘apply_transpose ≠ I on the bulk",
    )
    if carrying:
        from orpheus.numerics.coupled_system import CoupledField

        from tests.sn._test_helpers import joint_m_grid

        grid, _space = joint_m_grid(sn, A)
        b = _fresh_source(sn)                 # bulk-only source-subspace b
        b.interior.values[:] = rng.random(b.interior.values.shape)
        z_cot = RadialCharacteristicField.source_zeros_on(sn)
        fwd = grid.apply_transpose(
            grid.solve_transpose(CoupledField(systems=(b, z_cot))),
        ).systems[0]
    else:
        b = _fresh(sn)                        # bulk-only source-subspace b
        b.interior.values[:] = rng.random(b.interior.values.shape)
        fwd = A.apply_transpose(A.solve_transpose(b))
    np.testing.assert_allclose(
        np.asarray(fwd.interior.values), np.asarray(b.interior.values),
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
    carrying = sn.radial_characteristic_field_space is not None
    for g in range(sn.ng):
        M = _probe_augmented_matrix_one_group(sn, g)
        b = _fresh_source(sn) if carrying else _fresh(sn)
        b.interior.values[:, g] = rng.random((sn.quad.n_ordinates, *sn.spatial_shape))
        cot_in = _seed_cot(
            sn,
            rng.random(sn.radial_characteristic_field_space.shape[0])
            if carrying else None,
        )
        b_vec = _read_augmented(b, sn, g, cot_in)
        if carrying:
            # step 6: the joint M⁻ᵀ is the grid's transposed substitution.
            from orpheus.numerics.coupled_system import CoupledField

            from tests.sn._test_helpers import joint_m_grid

            grid, _space = joint_m_grid(sn, A)
            out = grid.solve_transpose(CoupledField(systems=(b, cot_in)))
            got = _read_augmented(out.systems[0], sn, g, out.systems[1])
        else:
            got = _read_augmented(A.solve_transpose(b), sn, g)
        ref = np.linalg.solve(M.T, b_vec)
        np.testing.assert_allclose(
            got[mask], ref[mask], rtol=_RTOL, atol=1e-12,
            err_msg=f"{geom} g{g}: solve_transpose ≠ M⁻ᵀ b (source-carried)",
        )


# ── G3 — full-field Euclidean solve-reciprocity (bulk ⊕ boundary) ──────


def _fdot(x, y, sn) -> float:
    """Plain Euclidean pairing over bulk ⊕ boundary.  ``solve_transpose`` is the
    EUCLIDEAN transpose-solve, so its reciprocity is metric-free (contrast the
    G-adjoint ``.H`` reciprocity in ``test_g_adjoint_reciprocity``)."""
    d = float(np.sum(np.asarray(x.interior.values) * np.asarray(y.interior.values)))
    for face in sn.angular_trace.layout.faces:
        d += float(np.sum(x.boundary.face_view(face) * y.boundary.face_view(face)))
    return d


@pytest.mark.foundation
@pytest.mark.parametrize("geom", ["slab", "cyl_product", "cyl_ls"])
def test_g3_full_field_solve_reciprocity(geom):
    r"""``⟨A.solve(q), p⟩ = ⟨q, A.solve_transpose(p)⟩`` over bulk ⊕ boundary.

    The ONLY gate that pins the boundary-cotangent passthrough: the bulk-only
    pairing fails at ~14–19 % (the boundary carries load-bearing inner-product
    mass), so a dropped/​wrong boundary cotangent reds here at O(1) while the
    bulk-only G1/G2 stay green (the #284 boundary subspace).  Caught the
    degenerate-ord boundary passthrough drop (cyl_product 2.8e-2 → 0 exact).

    Parametrized over the non-carrying configs + slab, where ``_fdot``
    (bulk ⊕ boundary) is the complete Euclidean pairing.  The sphere's
    seed-carrying solve-reciprocity (the #284 seed free-DOF, Mode-12
    seed-blind) is the test-architect's #29 domain.
    """
    sn = _MESHES[geom]()
    A = _loss(sn)
    rng = np.random.default_rng(20260705)
    q, p = _fresh(sn), _fresh(sn)
    q.interior.values[:] = rng.random(q.interior.values.shape)
    p.interior.values[:] = rng.random(p.interior.values.shape)
    for face in sn.angular_trace.layout.faces:
        q.boundary.face_view(face)[:] = rng.random(q.boundary.face_view(face).shape)
        p.boundary.face_view(face)[:] = rng.random(p.boundary.face_view(face).shape)
    lhs = _fdot(A.solve(q), p, sn)
    rhs = _fdot(q, A.solve_transpose(p), sn)
    rel = abs(lhs - rhs) / (abs(lhs) + abs(rhs) + 1e-300)
    if not rel < _RTOL:
        pytest.fail(
            f"{geom}: ⟨A.solve q,p⟩={lhs:.6e} ≠ ⟨q,A.solve_transpose p⟩="
            f"{rhs:.6e} (rel={rel:.2e}) — boundary-cotangent reciprocity broken"
        )


# ── G4 — the non-carrying-cylinder seed contract (m_seed = None) ───────


@pytest.mark.foundation
@pytest.mark.parametrize("geom", ["cyl_product", "cyl_ls"])
def test_g4_cyl_returns_no_seed_cotangent(geom):
    """The non-carrying cylinder's transpose-solve returns
    ``radial_characteristic = None`` (``m_seed = None``) — the R12a contract, the
    mirror of the forward's non-carrying refusal."""
    sn = _MESHES[geom]()
    if sn.radial_characteristic_field_space is not None:
        pytest.fail(f"{geom} unexpectedly carries a starting-direction space")
    A = _loss(sn)
    p = _fresh(sn)
    p.interior.values[:] = np.random.default_rng(1).random(p.interior.values.shape)
    out = A.solve_transpose(p)
    if type(out) is not FullField:
        pytest.fail(f"{geom}: solve_transpose emitted {type(out).__name__}")
    # B.2d: the R12a mirror is STRUCTURAL — the transposed ψ½ legs cannot
    # even be BUILT on a non-carrying mesh (the leaf factories raise; pinned
    # in test_radial_characteristic_carrier), so a silent None pass-through
    # is unspellable. The live claim here is the clean 2-block run above.


# ── G5 — the mandatory-config activation sentinel (Mode-7 pin) ─────────


def _m0_mm_a_in(sn) -> float:
    """The seed-ordinate M-M coefficient on level 0 (live iff a product seed-fold)."""
    from orpheus.sn.loss_representation import _OneDimScanWalk
    w = _OneDimScanWalk(sn)
    geom = w._ensure_geom_cache()
    m0_local = sn.pole_angular_closure._edge_seed_stencil(0)[0]
    m0 = int(np.asarray(sn.quad.level_indices[0])[m0_local])
    return float(geom.mm_a_in_coeff[m0])


@pytest.mark.foundation
def test_g5_mandatory_config_activates_cyl_terms():
    """``cyl_product`` activates BOTH cylinder-specific transpose terms
    (degenerate pure-azimuthal ords + a LIVE seed-fold); ``cyl_ls`` NULLS both.

    The Mode-7 activation pin: a future quadrature change cannot silently decay
    the mandatory config into the blind control (the ERR-066 blindness that hid
    the degenerate-drop in every ``level_symmetric`` row) without reddening
    here."""
    from orpheus.sn.loss_representation import _OneDimScanWalk
    sn_p = _cyl_product()
    n_deg_p = _OneDimScanWalk(sn_p)._degenerate_positions()[0].size
    if not (n_deg_p > 0 and _m0_mm_a_in(sn_p) != 0.0):
        pytest.fail(
            f"cyl_product must activate BOTH terms: degenerate ords={n_deg_p}, "
            f"seed mm_a_in[m0]={_m0_mm_a_in(sn_p)} (need >0 and ≠0)"
        )
    sn_l = _cyl_ls()
    n_deg_l = _OneDimScanWalk(sn_l)._degenerate_positions()[0].size
    if not (n_deg_l == 0 and _m0_mm_a_in(sn_l) == 0.0):
        pytest.fail(
            f"cyl_ls control must NULL BOTH terms: degenerate ords={n_deg_l}, "
            f"seed mm_a_in[m0]={_m0_mm_a_in(sn_l)} (need 0 and 0)"
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
        b.interior.values[:, g] = rng.random((sn.quad.n_ordinates, *sn.spatial_shape))
        b_vec = np.asarray(b.interior.values)[:, g].ravel()
        got = np.asarray(A.solve_transpose(b).interior.values)[:, g].ravel()
        ref = np.empty_like(b_vec)
        ref[order] = solve_triangular(
            M[np.ix_(order, order)].T, b_vec[order], lower=False,
        )
        np.testing.assert_allclose(
            got, ref, rtol=_RTOL, atol=1e-12,
            err_msg=f"slab g{g}: solve_transpose ≠ LAPACK back-substitution",
        )
