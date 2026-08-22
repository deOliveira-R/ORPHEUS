r"""L0/L2 — :class:`TimedFullField` cross-method generic container.

Pins the post-D-H.1 contract for
:class:`orpheus.transport.timed_full_field.TimedFullField`:

* Bulk + boundary pair construction with locus-type validation (bulk
  must be a :class:`~orpheus.transport.fields._bases.BulkField`,
  boundary a :class:`~orpheus.transport.fields._bases.AngularBoundaryField`)
  and mesh-identity validation.
* Algebra propagates to bulk + boundary; history dropped on algebra
  results.
* Cross-method composition rejected (bulk-type-identity gate).
* History shift-register (``advance``, ``at_lag``, ``history_length``).
* Frozen-instance contract (assignment raises).

These are foundation tests — they pin the data-shape and algebra
contract of the cross-method generic container. Method-specific
factories (SN's :meth:`~orpheus.transport.timed_full_field.TimedFullField.zeros`) get their own
tests.

References
----------

* Depth B plan §3.8 (Container architecture).
* `.claude/agent-memory/test-architect/dg_boundary_flux_pure_field_verification.md`
  (D-G predecessor tests).
"""
from __future__ import annotations

from dataclasses import FrozenInstanceError

import numpy as np
import pytest

from orpheus.geometry import BC, CoordSystem, Mesh1D, Mesh2D
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.transport.fields.angular_flux import AngularFlux
from orpheus.transport.fields.angular_boundary_flux import AngularBoundaryFlux
from orpheus.transport.fields.scalar_flux import ScalarFlux
from orpheus.transport.timed_full_field import TimedFullField

from tests.sn._test_helpers import placeholder_materials


pytestmark = [pytest.mark.foundation]


# ───────────────────────────────────────────────────────────────────────
# Mesh fixtures
# ───────────────────────────────────────────────────────────────────────


def _slab_mesh(nx: int = 4, ng: int = 2) -> SNMesh:
    mesh = Mesh1D(
        edges=np.linspace(0.0, 1.0, nx + 1),
        mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.CARTESIAN,
        bc_left=BC("vacuum"),
        bc_right=BC("vacuum"),
    )
    quad = Quadrature.gauss_legendre(n_ordinates=4)
    return SNMesh(mesh, quad, placeholder_materials(ng=ng))


def _stretched_mesh(nx: int = 4, ng: int = 2) -> SNMesh:
    """Same shape as ``_slab_mesh``, doubled width — the cell VOLUMES differ,
    so the carrier mints an UNEQUAL space (the F2 content discriminator)."""
    mesh = Mesh1D(
        edges=np.linspace(0.0, 2.0, nx + 1),
        mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.CARTESIAN,
        bc_left=BC("vacuum"),
        bc_right=BC("vacuum"),
    )
    quad = Quadrature.gauss_legendre(n_ordinates=4)
    return SNMesh(mesh, quad, placeholder_materials(ng=ng))


def _cartesian_2d_mesh(nx: int = 3, ny: int = 2, ng: int = 2) -> SNMesh:
    mesh = Mesh2D(
        edges_x=np.linspace(0, 1, nx + 1),
        edges_y=np.linspace(0, 1, ny + 1),
        mat_map=np.zeros((nx, ny), dtype=int),
    )
    quad = Quadrature.level_symmetric(sn_order=4)
    return SNMesh(mesh, quad, placeholder_materials(ng=ng))


# ───────────────────────────────────────────────────────────────────────
# A. Construction
# ───────────────────────────────────────────────────────────────────────


class TestConstruction:
    def test_factory_slab(self) -> None:
        m = _slab_mesh()
        state = TimedFullField.zeros(interior=AngularFlux, boundary=AngularBoundaryFlux, mesh=m)
        assert isinstance(state, TimedFullField)
        assert isinstance(state.interior, AngularFlux)
        assert isinstance(state.boundary, AngularBoundaryFlux)
        assert state.interior.mesh is m
        assert state.boundary.mesh is m
        np.testing.assert_array_equal(state.interior.values, 0.0)
        np.testing.assert_array_equal(state.boundary.values, 0.0)
        assert state._history == ()
        assert state.history_depth == 2

    def test_factory_2d(self) -> None:
        m = _cartesian_2d_mesh()
        state = TimedFullField.zeros(interior=AngularFlux, boundary=AngularBoundaryFlux, mesh=m)
        assert isinstance(state, TimedFullField)
        assert state.interior.values.shape == (m.quad.N, m.ng, *m.spatial_shape)
        assert state.boundary.layout.total_size == state.boundary.values.size

    def test_factory_custom_history_depth(self) -> None:
        m = _slab_mesh()
        state = TimedFullField.zeros(interior=AngularFlux, boundary=AngularBoundaryFlux, mesh=m, history_depth=5)
        assert state.history_depth == 5
        assert state._history == ()

    def test_rejects_non_field_bulk(self) -> None:
        m = _slab_mesh()
        bf = AngularBoundaryFlux.zeros_on(m)
        with pytest.raises(TypeError, match="bulk must be a BulkField"):
            TimedFullField(interior="not a field", boundary=bf)  # type: ignore[arg-type]

    def test_rejects_non_field_boundary(self) -> None:
        m = _slab_mesh()
        psi = AngularFlux.zeros_on(m)
        with pytest.raises(TypeError, match="boundary must be a BoundaryField"):
            TimedFullField(interior=psi, boundary="not a field")  # type: ignore[arg-type]

    def test_rejects_boundary_field_as_bulk(self) -> None:
        """B.6 locus gate: a boundary-locus field cannot be slotted as bulk.

        ``AngularBoundaryFlux`` is a :class:`AngularBoundaryField`, NOT a
        :class:`BulkField`. Even though it is a perfectly valid
        ``Field``, slotting it into the bulk position is a
        locus-confusion bug the type now makes unrepresentable
        (illegal-states-unrepresentable, ``coding-elegance`` Pattern 4).
        """
        m = _slab_mesh()
        bf = AngularBoundaryFlux.zeros_on(m)
        with pytest.raises(TypeError, match="bulk must be a BulkField"):
            TimedFullField(interior=bf, boundary=bf)  # type: ignore[arg-type]

    def test_rejects_bulk_field_as_boundary(self) -> None:
        """B.6 locus gate: a bulk-locus field cannot be slotted as boundary.

        ``AngularFlux`` is a :class:`BulkField`, NOT a
        :class:`AngularBoundaryField`. The bulk slot accepts it; the boundary
        slot must reject it (the dual of
        :meth:`test_rejects_boundary_field_as_bulk`).
        """
        m = _slab_mesh()
        psi = AngularFlux.zeros_on(m)
        with pytest.raises(TypeError, match="boundary must be a BoundaryField"):
            TimedFullField(interior=psi, boundary=psi)  # type: ignore[arg-type]

    def test_twin_mesh_blocks_mix_and_derive_the_carrier_space(self) -> None:
        # CS4b S4 (F4 — the F2 doctrine's composite leg): twin carriers
        # mint content-equal spaces, so blocks assemble across them; the
        # pre-S4 mesh-IDENTITY refusal retired (cross-slot coherence is
        # an admission question at the seams that hold a carrier
        # reference, not a construction question).
        m1 = _slab_mesh()
        m2 = _slab_mesh()  # different instance, same structure
        psi = AngularFlux.zeros_on(m1)
        bf = AngularBoundaryFlux.zeros_on(m2)
        state = TimedFullField(interior=psi, boundary=bf)
        if state.space != m1.full_field_space:
            pytest.fail("twin-carrier composite must content-equal the mint")

    def test_rejects_negative_history_depth(self) -> None:
        m = _slab_mesh()
        psi = AngularFlux.zeros_on(m)
        bf = AngularBoundaryFlux.zeros_on(m)
        with pytest.raises(ValueError, match="history_depth"):
            TimedFullField(interior=psi, boundary=bf, history_depth=-1)


# ───────────────────────────────────────────────────────────────────────
# B. Algebra propagates to bulk + boundary
# ───────────────────────────────────────────────────────────────────────


def _filled_state(m: SNMesh, bulk_val: float, bound_val: float) -> TimedFullField:
    state = TimedFullField.zeros(interior=AngularFlux, boundary=AngularBoundaryFlux, mesh=m)
    state.interior.values[:] = bulk_val
    state.boundary.values[:] = bound_val
    return state


class TestAlgebraPropagation:
    def test_add_flux_flux_propagates_blockwise_slab(self) -> None:
        """Composite delegation of the V algebra (campaign 1 CS3): ψ + ψ'
        propagates member-wise to bulk + boundary, and the update
        ψ + (ψ' − ψ) recovers ψ' — the same delegation discipline that
        used to carry the affine gate now carries the plain vector ops."""
        m = _slab_mesh()
        a = _filled_state(m, bulk_val=1.0, bound_val=2.0)
        b = _filled_state(m, bulk_val=3.0, bound_val=4.0)
        s = a + b
        if not isinstance(s, TimedFullField):
            raise AssertionError("composite + composite left the composite type")
        np.testing.assert_array_equal(s.interior.values, 4.0)
        np.testing.assert_array_equal(s.boundary.values, 6.0)
        d = b - a            # composite difference (bulk + boundary)
        out = a + d          # the update step → composite flux
        assert isinstance(out, TimedFullField)
        np.testing.assert_array_almost_equal_nulp(out.interior.values, b.interior.values, nulp=4)
        np.testing.assert_array_almost_equal_nulp(
            out.boundary.values, b.boundary.values, nulp=4,
        )
        assert out._history == ()  # algebra drops history

    def test_sub_propagates_2d(self) -> None:
        m = _cartesian_2d_mesh()
        a = _filled_state(m, bulk_val=5.0, bound_val=10.0)
        b = _filled_state(m, bulk_val=2.0, bound_val=3.0)
        out = a - b
        np.testing.assert_array_equal(out.interior.values, 3.0)
        np.testing.assert_array_equal(out.boundary.values, 7.0)

    def test_neg_propagates(self) -> None:
        m = _slab_mesh()
        a = _filled_state(m, bulk_val=1.5, bound_val=2.5)
        out = -a
        np.testing.assert_array_equal(out.interior.values, -1.5)
        np.testing.assert_array_equal(out.boundary.values, -2.5)

    def test_scalar_mul_propagates(self) -> None:
        m = _slab_mesh()
        a = _filled_state(m, bulk_val=2.0, bound_val=4.0)
        out = a * 3.0
        np.testing.assert_array_equal(out.interior.values, 6.0)
        np.testing.assert_array_equal(out.boundary.values, 12.0)

    def test_rmul_propagates(self) -> None:
        m = _slab_mesh()
        a = _filled_state(m, bulk_val=2.0, bound_val=4.0)
        out = 3.0 * a
        np.testing.assert_array_equal(out.interior.values, 6.0)
        np.testing.assert_array_equal(out.boundary.values, 12.0)

    def test_truediv_propagates(self) -> None:
        m = _slab_mesh()
        a = _filled_state(m, bulk_val=8.0, bound_val=12.0)
        out = a / 4.0
        np.testing.assert_array_equal(out.interior.values, 2.0)
        np.testing.assert_array_equal(out.boundary.values, 3.0)

    def test_displacement_distributive_property(self) -> None:
        """Distributivity (d1 + d2)·c == d1·c + d2·c in the composite V —
        since the CS3 cone carve the differences are ordinary same-typed
        composites, so this is a law of the ONE composite vector space."""
        m = _slab_mesh()
        base = _filled_state(m, bulk_val=0.0, bound_val=0.0)
        a = _filled_state(m, bulk_val=1.0, bound_val=2.0)
        b = _filled_state(m, bulk_val=3.0, bound_val=4.0)
        d1 = a - base  # composite displacement (bulk + boundary)
        d2 = b - base
        left = (d1 + d2) * 2.0
        right = d1 * 2.0 + d2 * 2.0
        np.testing.assert_allclose(left.interior.values, right.interior.values, rtol=1e-15)
        np.testing.assert_allclose(left.boundary.values, right.boundary.values, rtol=1e-15)


# ───────────────────────────────────────────────────────────────────────
# C. Cross-class rejection
# ───────────────────────────────────────────────────────────────────────


class TestCrossClassRejection:
    def test_cross_carrier_discrimination_is_space_content(self) -> None:
        """CS4b S3 (F2): identity is space CONTENT, not carrier provenance.
        Twin carriers mint EQUAL spaces — composites mix member-wise (the
        EQUAL leg); a carrier whose cell volumes differ mints UNEQUAL bulk
        spaces and both binary ops refuse (the UNEQUAL leg)."""
        m1, m2 = _slab_mesh(), _slab_mesh()
        a = TimedFullField.zeros(interior=AngularFlux, boundary=AngularBoundaryFlux, mesh=m1)
        b = TimedFullField.zeros(interior=AngularFlux, boundary=AngularBoundaryFlux, mesh=m2)
        _ = a - b  # twin content — legal since the F2 re-key
        _ = a + b
        c = TimedFullField.zeros(
            interior=AngularFlux, boundary=AngularBoundaryFlux,
            mesh=_stretched_mesh(),
        )
        with pytest.raises(ValueError, match="equal space"):
            a - c
        with pytest.raises(ValueError, match="equal space"):
            a + c

    def test_wrong_type_rejected(self) -> None:
        m = _slab_mesh()
        state = TimedFullField.zeros(interior=AngularFlux, boundary=AngularBoundaryFlux, mesh=m)
        with pytest.raises(TypeError, match="same-class partner"):
            state + 42  # type: ignore[operator]


# ───────────────────────────────────────────────────────────────────────
# D. History shift-register
# ───────────────────────────────────────────────────────────────────────


class TestHistoryShiftRegister:
    def test_advance_pushes_current_to_lag_1(self) -> None:
        m = _slab_mesh()
        state = _filled_state(m, bulk_val=1.0, bound_val=2.0)
        new_bulk = AngularFlux.from_mesh(
            np.full((m.quad.N, m.ng, *m.spatial_shape), 10.0), m,
        )
        new_boundary = AngularBoundaryFlux.from_face_arrays(
            m, {
                "xmin": np.full((m.quad.N, m.ng), 20.0),
                "xmax": np.full((m.quad.N, m.ng), 20.0),
            },
        )
        advanced = state.advance(new_bulk, new_boundary)

        # Current frame = new values.
        np.testing.assert_array_equal(advanced.interior.values, 10.0)
        np.testing.assert_array_equal(advanced.boundary.values, 20.0)

        # Lag-1 = the pre-advance state.
        lag_1 = advanced.at_lag(1)
        np.testing.assert_array_equal(lag_1.interior.values, 1.0)
        np.testing.assert_array_equal(lag_1.boundary.values, 2.0)

    def test_advance_chained(self) -> None:
        """Multiple advances rotate the history correctly."""
        m = _slab_mesh()
        state = _filled_state(m, bulk_val=1.0, bound_val=10.0)
        for i in range(2, 5):
            new_bulk = AngularFlux.from_mesh(
                np.full((m.quad.N, m.ng, *m.spatial_shape), float(i)), m,
            )
            new_boundary = AngularBoundaryFlux.from_face_arrays(
                m, {
                    "xmin": np.full((m.quad.N, m.ng), float(i) * 10),
                    "xmax": np.full((m.quad.N, m.ng), float(i) * 10),
                },
            )
            state = state.advance(new_bulk, new_boundary)
        # Current = 4.0; lag-1 = 3.0; (history_depth=2 trims older).
        np.testing.assert_array_equal(state.interior.values, 4.0)
        np.testing.assert_array_equal(state.at_lag(1).interior.values, 3.0)
        np.testing.assert_array_equal(state.at_lag(2).interior.values, 2.0)
        assert state.history_length == 2  # trimmed to history_depth

    def test_advance_trims_to_history_depth(self) -> None:
        m = _slab_mesh()
        state = TimedFullField.zeros(interior=AngularFlux, boundary=AngularBoundaryFlux, mesh=m, history_depth=2)
        for i in range(5):
            nb = AngularFlux.zeros_on(m)
            nf = AngularBoundaryFlux.zeros_on(m)
            state = state.advance(nb, nf)
        assert state.history_length == 2

    def test_advance_type_mismatch_raises(self) -> None:
        m = _slab_mesh()
        state = TimedFullField.zeros(interior=AngularFlux, boundary=AngularBoundaryFlux, mesh=m)
        # Pass a ScalarFlux as bulk — type mismatch.
        sf = ScalarFlux.zeros_on(m)
        bf = AngularBoundaryFlux.zeros_on(m)
        with pytest.raises(TypeError, match="new_bulk type"):
            state.advance(sf, bf)  # type: ignore[arg-type]

    def test_at_lag_zero_returns_self(self) -> None:
        m = _slab_mesh()
        state = TimedFullField.zeros(interior=AngularFlux, boundary=AngularBoundaryFlux, mesh=m)
        assert state.at_lag(0) is state

    def test_at_lag_out_of_range_raises(self) -> None:
        m = _slab_mesh()
        state = TimedFullField.zeros(interior=AngularFlux, boundary=AngularBoundaryFlux, mesh=m)
        with pytest.raises(IndexError):
            state.at_lag(5)

    def test_at_lag_negative_raises(self) -> None:
        m = _slab_mesh()
        state = TimedFullField.zeros(interior=AngularFlux, boundary=AngularBoundaryFlux, mesh=m)
        with pytest.raises(ValueError, match="non-negative"):
            state.at_lag(-1)


# ───────────────────────────────────────────────────────────────────────
# E. Time-derivative use case (the load-bearing purpose of history)
# ───────────────────────────────────────────────────────────────────────


class TestTimeDerivativeStencil:
    """Pin the canonical use case for history: time-derivative stencils.

    The first-order backward difference is ``(state[n] - state[n-1]) / dt``.
    Higher-order BDF stencils read more lags.
    """

    def test_first_order_backward_difference(self) -> None:
        m = _slab_mesh()
        state_old = _filled_state(m, bulk_val=1.0, bound_val=2.0)
        new_bulk = AngularFlux.from_mesh(
            np.full((m.quad.N, m.ng, *m.spatial_shape), 1.5), m,
        )
        new_boundary = AngularBoundaryFlux.from_face_arrays(
            m, {
                "xmin": np.full((m.quad.N, m.ng), 2.5),
                "xmax": np.full((m.quad.N, m.ng), 2.5),
            },
        )
        state = state_old.advance(new_bulk, new_boundary)

        dt = 0.1
        # dpsi/dt ≈ (state[n] - state[n-1]) / dt
        dpsi_dt = (state.at_lag(0) - state.at_lag(1)) / dt
        # Bulk: (1.5 - 1.0) / 0.1 = 5.0
        np.testing.assert_allclose(dpsi_dt.interior.values, 5.0, rtol=1e-15)
        # Boundary: (2.5 - 2.0) / 0.1 = 5.0
        np.testing.assert_allclose(dpsi_dt.boundary.values, 5.0, rtol=1e-15)


# ───────────────────────────────────────────────────────────────────────
# F. Frozen-instance contract
# ───────────────────────────────────────────────────────────────────────


class TestFrozenInstance:
    def test_assign_bulk_raises(self) -> None:
        m = _slab_mesh()
        state = TimedFullField.zeros(interior=AngularFlux, boundary=AngularBoundaryFlux, mesh=m)
        with pytest.raises(FrozenInstanceError):
            state.interior = state.interior  # type: ignore[misc]

    def test_assign_boundary_raises(self) -> None:
        m = _slab_mesh()
        state = TimedFullField.zeros(interior=AngularFlux, boundary=AngularBoundaryFlux, mesh=m)
        with pytest.raises(FrozenInstanceError):
            state.boundary = state.boundary  # type: ignore[misc]


# ───────────────────────────────────────────────────────────────────────
# G. Copy
# ───────────────────────────────────────────────────────────────────────


class TestCopy:
    def test_copy_returns_owned_arrays(self) -> None:
        m = _slab_mesh()
        state = _filled_state(m, bulk_val=1.0, bound_val=2.0)
        copy = state.copy()
        # Same values
        np.testing.assert_array_equal(copy.interior.values, state.interior.values)
        # But independent arrays — modify original; copy unchanged.
        state.interior.values[:] = 99.0
        np.testing.assert_array_equal(copy.interior.values, 1.0)

    def test_copy_drops_history(self) -> None:
        m = _slab_mesh()
        state = TimedFullField.zeros(interior=AngularFlux, boundary=AngularBoundaryFlux, mesh=m)
        nb = AngularFlux.zeros_on(m)
        nf = AngularBoundaryFlux.zeros_on(m)
        advanced = state.advance(nb, nf)
        assert advanced.history_length == 1
        copy = advanced.copy()
        assert copy.history_length == 0


# ───────────────────────────────────────────────────────────────────────
# H. Flat-vector protocol (Krylov / scipy.gmres adapter)
# ───────────────────────────────────────────────────────────────────────


class TestFlatVectorProtocol:
    """Pin the to_flat / from_flat round-trip used by the Krylov adapter."""

    def test_to_flat_returns_concatenated_1d(self) -> None:
        m = _slab_mesh()
        state = _filled_state(m, bulk_val=1.0, bound_val=2.0)
        flat = state.to_flat()
        assert flat.ndim == 1
        expected_size = state.interior.values.size + state.boundary.values.size
        assert flat.size == expected_size

    def test_to_flat_preserves_values(self) -> None:
        m = _slab_mesh()
        state = _filled_state(m, bulk_val=3.0, bound_val=7.0)
        flat = state.to_flat()
        n_interior = state.interior.values.size
        np.testing.assert_array_equal(flat[:n_interior], 3.0)
        np.testing.assert_array_equal(flat[n_interior:], 7.0)

    def test_from_flat_round_trip_slab(self) -> None:
        m = _slab_mesh()
        rng = np.random.default_rng(42)
        state = _filled_state(m, bulk_val=0.0, bound_val=0.0)
        state.interior.values[:] = rng.standard_normal(state.interior.values.shape)
        state.boundary.values[:] = rng.standard_normal(state.boundary.values.shape)
        flat = state.to_flat()
        reconstructed = TimedFullField.from_flat(flat, state)
        np.testing.assert_array_equal(
            reconstructed.interior.values, state.interior.values,
        )
        np.testing.assert_array_equal(
            reconstructed.boundary.values, state.boundary.values,
        )

    def test_from_flat_round_trip_2d(self) -> None:
        m = _cartesian_2d_mesh()
        rng = np.random.default_rng(7)
        state = TimedFullField.zeros(interior=AngularFlux, boundary=AngularBoundaryFlux, mesh=m)
        state.interior.values[:] = rng.standard_normal(state.interior.values.shape)
        state.boundary.values[:] = rng.standard_normal(state.boundary.values.shape)
        flat = state.to_flat()
        reconstructed = TimedFullField.from_flat(flat, state)
        np.testing.assert_array_equal(
            reconstructed.interior.values, state.interior.values,
        )
        np.testing.assert_array_equal(
            reconstructed.boundary.values, state.boundary.values,
        )

    def test_from_flat_drops_history(self) -> None:
        m = _slab_mesh()
        state = TimedFullField.zeros(interior=AngularFlux, boundary=AngularBoundaryFlux, mesh=m)
        nb = AngularFlux.zeros_on(m)
        nf = AngularBoundaryFlux.zeros_on(m)
        advanced = state.advance(nb, nf)
        assert advanced.history_length == 1
        reconstructed = TimedFullField.from_flat(advanced.to_flat(), advanced)
        assert reconstructed.history_length == 0  # Krylov iterates carry no history

    def test_from_flat_preserves_history_depth(self) -> None:
        m = _slab_mesh()
        state = TimedFullField.zeros(interior=AngularFlux, boundary=AngularBoundaryFlux, mesh=m, history_depth=5)
        reconstructed = TimedFullField.from_flat(state.to_flat(), state)
        assert reconstructed.history_depth == 5

    def test_from_flat_wrong_size_raises(self) -> None:
        m = _slab_mesh()
        state = TimedFullField.zeros(interior=AngularFlux, boundary=AngularBoundaryFlux, mesh=m)
        with pytest.raises(ValueError, match="flat.size"):
            TimedFullField.from_flat(np.zeros(3), state)

    def test_iteration_protocol_detection(self) -> None:
        """The Krylov adapter's _is_ravellable detects TimedFullField."""
        from orpheus.numerics.iteration import (
            _is_ravellable, _ravel, _unravel_like, _zeros_like, _l2_norm,
        )
        m = _slab_mesh()
        state = _filled_state(m, bulk_val=1.0, bound_val=2.0)
        assert _is_ravellable(state)
        flat = _ravel(state)
        assert flat.ndim == 1
        rec = _unravel_like(state, flat)
        assert isinstance(rec, TimedFullField)
        zeros = _zeros_like(state)
        np.testing.assert_array_equal(zeros.interior.values, 0.0)
        np.testing.assert_array_equal(zeros.boundary.values, 0.0)
        # L2 norm reads the flat representation directly.
        l2 = _l2_norm(state)
        expected = float(np.linalg.norm(state.to_flat()))
        np.testing.assert_allclose(l2, expected, rtol=1e-15)
