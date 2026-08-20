r"""Intrinsic-property gate for the generic ``Composite[Interior, Boundary]`` base.

:class:`~orpheus.transport.full_field.Composite` is the structural two-block
carrier the transport operator algebra acts on; :class:`FullField` (SN + the ψ½
block) is its specialization. Because ``FullField`` OVERRIDES every algebra hook
(``_map_binary`` / ``_map_unary`` / ``_recombine`` / ``_flat_parts`` /
``_from_flat``) to thread its third block, the base's OWN 2-block methods are
exercised by no SN test — so this file pins them directly.

Two things at once (the carrier-generalization gate, campaign Phase A2):

* **Intrinsic vector-space laws** — a math-bearing type ships a test of its
  defining laws ([[feedback-test-intrinsic-properties]]): commutativity, the
  additive inverse, scalar distributivity, the ``to_flat`` / ``from_flat``
  round-trip, ``copy`` independence, same-class partner rejection.
* **The multi-instantiation gate (N2)** — the composite is built with SCALAR
  leaves (``ScalarFlux`` / ``ScalarBoundaryFlux``), NOT the SN ``AngularFlux``.
  A hidden ``AngularFlux`` assumption anywhere in the generic body would red
  here while the SN wall stays green (L20 sharpened to the leaf-type parameter).
  The pure ``Composite`` carries NO ``radial_characteristic`` slot — the ψ½
  block is the SN subclass's, not the generic's.

Reuses the diffusion 2-material / 4-cell / 2-group slab (its scalar leaves are
the only non-SN ``BulkField`` / ``BoundaryField`` family in the tree today).
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.common.xs_library import make_mixture
from orpheus.diffusion import DiffusionMesh
from orpheus.geometry import BC
from orpheus.geometry.mesh import Mesh1D
from orpheus.transport.fields.scalar_boundary_flux import ScalarBoundaryFlux
from orpheus.transport.fields.scalar_flux import ScalarFlux
from orpheus.transport.full_field import Composite, FullField

pytestmark = pytest.mark.foundation

# A heterogeneous 2-material, non-uniform 4-cell, 2-group slab (the diffusion
# fixture geometry) — its scalar composite is a pure ``Composite`` with NO ray.
_SIG_T_A = np.array([0.2181, 0.7850])
_SIG_T_B = np.array([0.3416, 0.9431])
_SIG_S_A = np.array([[0.1900, 0.0160], [0.0, 0.4200]])
_SIG_S_B = np.array([[0.1000, 0.0020], [0.0, 0.0500]])
_SIG_F_A = np.array([0.0024, 0.0489])
_EDGES = np.array([0.0, 0.5, 1.5, 3.0, 5.0])
_MAT_IDS = np.array([0, 1, 1, 0])


def _materials() -> dict[int, object]:
    mix_a = make_mixture(
        sig_t=_SIG_T_A, sig_c=_SIG_T_A - _SIG_F_A - _SIG_S_A.sum(axis=1),
        sig_f=_SIG_F_A, nu=np.array([2.54, 2.47]),
        chi=np.array([1.0, 0.0]), sig_s=_SIG_S_A,
    )
    mix_b = make_mixture(
        sig_t=_SIG_T_B, sig_c=_SIG_T_B - _SIG_S_B.sum(axis=1),
        sig_f=np.array([0.0, 0.0]), nu=np.array([0.0, 0.0]),
        chi=np.zeros(2), sig_s=_SIG_S_B,
    )
    return {0: mix_a, 1: mix_b}


@pytest.fixture
def mesh() -> DiffusionMesh:
    mesh1d = Mesh1D(
        edges=_EDGES, mat_ids=_MAT_IDS,
        bc_left=BC("reflective"), bc_right=BC("vacuum"),
    )
    return DiffusionMesh(mesh1d, _materials())


def _scalar_composite(mesh: DiffusionMesh, seed: int) -> "Composite":
    """A pure scalar ``Composite`` (interior ⊕ boundary, NO ψ½) on ``mesh``."""
    rng = np.random.default_rng(seed)
    return Composite(
        interior=ScalarFlux.from_mesh(rng.random((2, 4)) + 0.5, mesh),
        boundary=ScalarBoundaryFlux.from_mesh(
            rng.random(mesh.scalar_trace.shape[0]) + 0.1, mesh,
        ),
    )


# ── The type is a pure Composite (multi-instantiation / structural) ──


def test_pure_composite_is_not_a_fullfield(mesh) -> None:
    """The scalar carrier is a ``Composite``, NOT the SN ``FullField`` — the
    generic base is directly instantiable with a non-angular leaf pair."""
    c = _scalar_composite(mesh, 0)
    assert isinstance(c, Composite)
    assert not isinstance(c, FullField)


def test_pure_composite_has_no_ray_slot(mesh) -> None:
    """The ψ½ block belongs to the SN subclass, not the generic base — the pure
    ``Composite`` carries no ``radial_characteristic`` attribute at all."""
    c = _scalar_composite(mesh, 0)
    assert not hasattr(c, "radial_characteristic")


def test_mesh_reads_off_either_leaf(mesh) -> None:
    c = _scalar_composite(mesh, 0)
    assert c.mesh is mesh


# ── Vector-space intrinsic laws (base hooks exercised via SCALAR leaves) ──


def test_addition_propagates_per_block(mesh) -> None:
    """Two flux composites ADD member-wise (flux lives in V — campaign 1
    CS3; the retired affine gate used to refuse this through the same
    delegation), propagated through BOTH blocks of the generic base."""
    a = _scalar_composite(mesh, 1)
    b = _scalar_composite(mesh, 2)
    s = a + b
    np.testing.assert_array_equal(
        s.interior.values, a.interior.values + b.interior.values,
    )
    np.testing.assert_array_equal(
        s.boundary.values, a.boundary.values + b.boundary.values,
    )


def test_subtraction_is_same_typed_per_block(mesh) -> None:
    """flux − flux returns the same leaf types, propagated to interior AND
    boundary (until CS3 this minted the displacement siblings)."""
    a = _scalar_composite(mesh, 3)
    b = _scalar_composite(mesh, 4)
    d = a - b
    np.testing.assert_array_equal(
        d.interior.values, a.interior.values - b.interior.values,
    )
    np.testing.assert_array_equal(
        d.boundary.values, a.boundary.values - b.boundary.values,
    )


def test_update_step_recovers_the_point(mesh) -> None:
    """a + (b − a) == b, propagated to both blocks — the plain V recovery
    (spelled through the displacement mint until the CS3 carve)."""
    a = _scalar_composite(mesh, 5)
    b = _scalar_composite(mesh, 6)
    recovered = a + (b - a)  # flux + displacement → flux
    np.testing.assert_allclose(recovered.interior.values, b.interior.values)
    np.testing.assert_allclose(recovered.boundary.values, b.boundary.values)


def test_scalar_scales_both_blocks(mesh) -> None:
    a = _scalar_composite(mesh, 7)
    out = 2.0 * a
    np.testing.assert_allclose(out.interior.values, 2.0 * a.interior.values)
    np.testing.assert_allclose(out.boundary.values, 2.0 * a.boundary.values)


def test_mul_and_rmul_agree(mesh) -> None:
    a = _scalar_composite(mesh, 8)
    left, right = 3.0 * a, a * 3.0
    np.testing.assert_array_equal(left.interior.values, right.interior.values)
    np.testing.assert_array_equal(left.boundary.values, right.boundary.values)


def test_truediv_is_reciprocal_mul(mesh) -> None:
    a = _scalar_composite(mesh, 9)
    half, scaled = a / 2.0, a * 0.5
    np.testing.assert_allclose(half.interior.values, scaled.interior.values)
    np.testing.assert_allclose(half.boundary.values, scaled.boundary.values)


# ── Flat protocol (the base's _flat_parts / _from_flat hooks) ──


def test_principal_bulk_leaf_is_the_interior(mesh) -> None:
    """The composite's convergence-diagnostic carrier is its INTERIOR leaf
    (CS3-R: the carrier owns the convention; the iteration layer reads only
    ``principal_bulk_leaf``. The norm-convention pin — interior space norm,
    not composite flat — is tests/numerics/test_si_diagnostic_trajectory.py)."""
    c = _scalar_composite(mesh, 3)
    assert c.principal_bulk_leaf is c.interior


def test_to_flat_layout_is_interior_then_boundary(mesh) -> None:
    a = _scalar_composite(mesh, 10)
    flat = a.to_flat()
    expected = np.concatenate([a.interior.values.ravel(), a.boundary.values])
    np.testing.assert_array_equal(flat, expected)


def test_from_flat_roundtrip_is_exact(mesh) -> None:
    a = _scalar_composite(mesh, 11)
    back = Composite.from_flat(a.to_flat(), a)
    assert isinstance(back, Composite) and not isinstance(back, FullField)
    np.testing.assert_array_equal(back.interior.values, a.interior.values)
    np.testing.assert_array_equal(back.boundary.values, a.boundary.values)


def test_from_flat_rejects_wrong_size(mesh) -> None:
    a = _scalar_composite(mesh, 12)
    with pytest.raises(ValueError, match="flat.size"):
        Composite.from_flat(a.to_flat()[:-1], a)


# ── copy / partner discipline ──


def test_copy_is_value_independent(mesh) -> None:
    a = _scalar_composite(mesh, 13)
    c = a.copy()
    np.testing.assert_array_equal(c.interior.values, a.interior.values)
    c.interior.values[0, 0] += 1.0
    # the copy owns its arrays — mutating it must not touch the original.
    assert a.interior.values[0, 0] != c.interior.values[0, 0]


def test_rejects_non_composite_partner(mesh) -> None:
    a = _scalar_composite(mesh, 14)
    with pytest.raises(TypeError, match="same-class partner"):
        _ = a + 42  # type: ignore[operator]
