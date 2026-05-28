r"""L0/L2 — :class:`TimedFullField` cross-method generic container.

Pins the post-D-H.1 contract for
:class:`orpheus.transport.timed_full_field.TimedFullField`:

* Bulk + boundary pair construction with mesh-identity validation.
* Algebra propagates to bulk + boundary; history dropped on algebra
  results.
* Cross-method composition rejected (bulk-type-identity gate).
* History shift-register (``advance``, ``at_lag``, ``history_length``).
* Frozen-instance contract (assignment raises).

These are foundation tests — they pin the data-shape and algebra
contract of the cross-method generic container. Method-specific
factories (SN's :meth:`SNMesh.zeros_timed_full_field`) get their own
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
from orpheus.sn.geometry import SNMesh
from orpheus.transport.fields.angular_flux import AngularFlux
from orpheus.transport.fields.boundary_flux import BoundaryFlux
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
        state = m.zeros_timed_full_field()
        assert isinstance(state, TimedFullField)
        assert isinstance(state.bulk, AngularFlux)
        assert isinstance(state.boundary, BoundaryFlux)
        assert state.bulk.mesh is m
        assert state.boundary.mesh is m
        np.testing.assert_array_equal(state.bulk.values, 0.0)
        np.testing.assert_array_equal(state.boundary.values, 0.0)
        assert state._history == ()
        assert state.history_depth == 2

    def test_factory_2d(self) -> None:
        m = _cartesian_2d_mesh()
        state = m.zeros_timed_full_field()
        assert isinstance(state, TimedFullField)
        assert state.bulk.values.shape == (m.quad.N, m.ng, m.nx, m.ny)
        assert state.boundary.layout.total_size == state.boundary.values.size

    def test_factory_custom_history_depth(self) -> None:
        m = _slab_mesh()
        state = m.zeros_timed_full_field(history_depth=5)
        assert state.history_depth == 5
        assert state._history == ()

    def test_rejects_non_field_bulk(self) -> None:
        m = _slab_mesh()
        bf = BoundaryFlux.zeros_for_sn_mesh(m)
        with pytest.raises(TypeError, match="bulk must be a Field"):
            TimedFullField(bulk="not a field", boundary=bf)  # type: ignore[arg-type]

    def test_rejects_non_field_boundary(self) -> None:
        m = _slab_mesh()
        psi = AngularFlux.zeros_for_sn_mesh(m)
        with pytest.raises(TypeError, match="boundary must be a Field"):
            TimedFullField(bulk=psi, boundary="not a field")  # type: ignore[arg-type]

    def test_rejects_mismatched_mesh(self) -> None:
        m1 = _slab_mesh()
        m2 = _slab_mesh()  # different instance, same structure
        psi = AngularFlux.zeros_for_sn_mesh(m1)
        bf = BoundaryFlux.zeros_for_sn_mesh(m2)
        with pytest.raises(ValueError, match="mesh identity"):
            TimedFullField(bulk=psi, boundary=bf)

    def test_rejects_negative_history_depth(self) -> None:
        m = _slab_mesh()
        psi = AngularFlux.zeros_for_sn_mesh(m)
        bf = BoundaryFlux.zeros_for_sn_mesh(m)
        with pytest.raises(ValueError, match="history_depth"):
            TimedFullField(bulk=psi, boundary=bf, history_depth=-1)


# ───────────────────────────────────────────────────────────────────────
# B. Algebra propagates to bulk + boundary
# ───────────────────────────────────────────────────────────────────────


def _filled_state(m: SNMesh, bulk_val: float, bound_val: float) -> TimedFullField:
    state = m.zeros_timed_full_field()
    state.bulk.values[:] = bulk_val
    state.boundary.values[:] = bound_val
    return state


class TestAlgebraPropagation:
    def test_add_propagates_slab(self) -> None:
        m = _slab_mesh()
        a = _filled_state(m, bulk_val=1.0, bound_val=2.0)
        b = _filled_state(m, bulk_val=3.0, bound_val=4.0)
        out = a + b
        assert isinstance(out, TimedFullField)
        np.testing.assert_array_equal(out.bulk.values, 4.0)
        np.testing.assert_array_equal(out.boundary.values, 6.0)
        assert out._history == ()  # algebra drops history

    def test_sub_propagates_2d(self) -> None:
        m = _cartesian_2d_mesh()
        a = _filled_state(m, bulk_val=5.0, bound_val=10.0)
        b = _filled_state(m, bulk_val=2.0, bound_val=3.0)
        out = a - b
        np.testing.assert_array_equal(out.bulk.values, 3.0)
        np.testing.assert_array_equal(out.boundary.values, 7.0)

    def test_neg_propagates(self) -> None:
        m = _slab_mesh()
        a = _filled_state(m, bulk_val=1.5, bound_val=2.5)
        out = -a
        np.testing.assert_array_equal(out.bulk.values, -1.5)
        np.testing.assert_array_equal(out.boundary.values, -2.5)

    def test_scalar_mul_propagates(self) -> None:
        m = _slab_mesh()
        a = _filled_state(m, bulk_val=2.0, bound_val=4.0)
        out = a * 3.0
        np.testing.assert_array_equal(out.bulk.values, 6.0)
        np.testing.assert_array_equal(out.boundary.values, 12.0)

    def test_rmul_propagates(self) -> None:
        m = _slab_mesh()
        a = _filled_state(m, bulk_val=2.0, bound_val=4.0)
        out = 3.0 * a
        np.testing.assert_array_equal(out.bulk.values, 6.0)
        np.testing.assert_array_equal(out.boundary.values, 12.0)

    def test_truediv_propagates(self) -> None:
        m = _slab_mesh()
        a = _filled_state(m, bulk_val=8.0, bound_val=12.0)
        out = a / 4.0
        np.testing.assert_array_equal(out.bulk.values, 2.0)
        np.testing.assert_array_equal(out.boundary.values, 3.0)

    def test_distributive_property(self) -> None:
        m = _slab_mesh()
        a = _filled_state(m, bulk_val=1.0, bound_val=2.0)
        b = _filled_state(m, bulk_val=3.0, bound_val=4.0)
        left = (a + b) * 2.0
        right = a * 2.0 + b * 2.0
        np.testing.assert_allclose(left.bulk.values, right.bulk.values, rtol=1e-15)
        np.testing.assert_allclose(left.boundary.values, right.boundary.values, rtol=1e-15)


# ───────────────────────────────────────────────────────────────────────
# C. Cross-class rejection
# ───────────────────────────────────────────────────────────────────────


class TestCrossClassRejection:
    def test_cross_mesh_rejected(self) -> None:
        m1 = _slab_mesh()
        m2 = _slab_mesh()
        a = m1.zeros_timed_full_field()
        b = m2.zeros_timed_full_field()
        # Same TimedFullField class, same bulk/boundary types, but
        # different mesh instances. AngularFlux._check_partner catches it.
        with pytest.raises(ValueError, match="mesh-bound"):
            a + b

    def test_wrong_type_rejected(self) -> None:
        m = _slab_mesh()
        state = m.zeros_timed_full_field()
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
            np.full((m.quad.N, m.ng, m.nx, m.ny), 10.0), m,
        )
        new_boundary = BoundaryFlux.from_face_arrays(
            m, {
                "xmin": np.full((m.quad.N, m.ng), 20.0),
                "xmax": np.full((m.quad.N, m.ng), 20.0),
            },
        )
        advanced = state.advance(new_bulk, new_boundary)

        # Current frame = new values.
        np.testing.assert_array_equal(advanced.bulk.values, 10.0)
        np.testing.assert_array_equal(advanced.boundary.values, 20.0)

        # Lag-1 = the pre-advance state.
        lag_1 = advanced.at_lag(1)
        np.testing.assert_array_equal(lag_1.bulk.values, 1.0)
        np.testing.assert_array_equal(lag_1.boundary.values, 2.0)

    def test_advance_chained(self) -> None:
        """Multiple advances rotate the history correctly."""
        m = _slab_mesh()
        state = _filled_state(m, bulk_val=1.0, bound_val=10.0)
        for i in range(2, 5):
            new_bulk = AngularFlux.from_mesh(
                np.full((m.quad.N, m.ng, m.nx, m.ny), float(i)), m,
            )
            new_boundary = BoundaryFlux.from_face_arrays(
                m, {
                    "xmin": np.full((m.quad.N, m.ng), float(i) * 10),
                    "xmax": np.full((m.quad.N, m.ng), float(i) * 10),
                },
            )
            state = state.advance(new_bulk, new_boundary)
        # Current = 4.0; lag-1 = 3.0; (history_depth=2 trims older).
        np.testing.assert_array_equal(state.bulk.values, 4.0)
        np.testing.assert_array_equal(state.at_lag(1).bulk.values, 3.0)
        np.testing.assert_array_equal(state.at_lag(2).bulk.values, 2.0)
        assert state.history_length == 2  # trimmed to history_depth

    def test_advance_trims_to_history_depth(self) -> None:
        m = _slab_mesh()
        state = m.zeros_timed_full_field(history_depth=2)
        for i in range(5):
            nb = AngularFlux.zeros_for_sn_mesh(m)
            nf = BoundaryFlux.zeros_for_sn_mesh(m)
            state = state.advance(nb, nf)
        assert state.history_length == 2

    def test_advance_type_mismatch_raises(self) -> None:
        m = _slab_mesh()
        state = m.zeros_timed_full_field()
        # Pass a ScalarFlux as bulk — type mismatch.
        sf = m.zeros_scalar_flux()
        bf = BoundaryFlux.zeros_for_sn_mesh(m)
        with pytest.raises(TypeError, match="new_bulk type"):
            state.advance(sf, bf)  # type: ignore[arg-type]

    def test_at_lag_zero_returns_self(self) -> None:
        m = _slab_mesh()
        state = m.zeros_timed_full_field()
        assert state.at_lag(0) is state

    def test_at_lag_out_of_range_raises(self) -> None:
        m = _slab_mesh()
        state = m.zeros_timed_full_field()
        with pytest.raises(IndexError):
            state.at_lag(5)

    def test_at_lag_negative_raises(self) -> None:
        m = _slab_mesh()
        state = m.zeros_timed_full_field()
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
            np.full((m.quad.N, m.ng, m.nx, m.ny), 1.5), m,
        )
        new_boundary = BoundaryFlux.from_face_arrays(
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
        np.testing.assert_allclose(dpsi_dt.bulk.values, 5.0, rtol=1e-15)
        # Boundary: (2.5 - 2.0) / 0.1 = 5.0
        np.testing.assert_allclose(dpsi_dt.boundary.values, 5.0, rtol=1e-15)


# ───────────────────────────────────────────────────────────────────────
# F. Frozen-instance contract
# ───────────────────────────────────────────────────────────────────────


class TestFrozenInstance:
    def test_assign_bulk_raises(self) -> None:
        m = _slab_mesh()
        state = m.zeros_timed_full_field()
        with pytest.raises(FrozenInstanceError):
            state.bulk = state.bulk  # type: ignore[misc]

    def test_assign_boundary_raises(self) -> None:
        m = _slab_mesh()
        state = m.zeros_timed_full_field()
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
        np.testing.assert_array_equal(copy.bulk.values, state.bulk.values)
        # But independent arrays — modify original; copy unchanged.
        state.bulk.values[:] = 99.0
        np.testing.assert_array_equal(copy.bulk.values, 1.0)

    def test_copy_drops_history(self) -> None:
        m = _slab_mesh()
        state = m.zeros_timed_full_field()
        nb = AngularFlux.zeros_for_sn_mesh(m)
        nf = BoundaryFlux.zeros_for_sn_mesh(m)
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
        expected_size = state.bulk.values.size + state.boundary.values.size
        assert flat.size == expected_size

    def test_to_flat_preserves_values(self) -> None:
        m = _slab_mesh()
        state = _filled_state(m, bulk_val=3.0, bound_val=7.0)
        flat = state.to_flat()
        n_bulk = state.bulk.values.size
        np.testing.assert_array_equal(flat[:n_bulk], 3.0)
        np.testing.assert_array_equal(flat[n_bulk:], 7.0)

    def test_from_flat_round_trip_slab(self) -> None:
        m = _slab_mesh()
        rng = np.random.default_rng(42)
        state = _filled_state(m, bulk_val=0.0, bound_val=0.0)
        state.bulk.values[:] = rng.standard_normal(state.bulk.values.shape)
        state.boundary.values[:] = rng.standard_normal(state.boundary.values.shape)
        flat = state.to_flat()
        reconstructed = TimedFullField.from_flat(flat, state)
        np.testing.assert_array_equal(
            reconstructed.bulk.values, state.bulk.values,
        )
        np.testing.assert_array_equal(
            reconstructed.boundary.values, state.boundary.values,
        )

    def test_from_flat_round_trip_2d(self) -> None:
        m = _cartesian_2d_mesh()
        rng = np.random.default_rng(7)
        state = m.zeros_timed_full_field()
        state.bulk.values[:] = rng.standard_normal(state.bulk.values.shape)
        state.boundary.values[:] = rng.standard_normal(state.boundary.values.shape)
        flat = state.to_flat()
        reconstructed = TimedFullField.from_flat(flat, state)
        np.testing.assert_array_equal(
            reconstructed.bulk.values, state.bulk.values,
        )
        np.testing.assert_array_equal(
            reconstructed.boundary.values, state.boundary.values,
        )

    def test_from_flat_drops_history(self) -> None:
        m = _slab_mesh()
        state = m.zeros_timed_full_field()
        nb = AngularFlux.zeros_for_sn_mesh(m)
        nf = BoundaryFlux.zeros_for_sn_mesh(m)
        advanced = state.advance(nb, nf)
        assert advanced.history_length == 1
        reconstructed = TimedFullField.from_flat(advanced.to_flat(), advanced)
        assert reconstructed.history_length == 0  # Krylov iterates carry no history

    def test_from_flat_preserves_history_depth(self) -> None:
        m = _slab_mesh()
        state = m.zeros_timed_full_field(history_depth=5)
        reconstructed = TimedFullField.from_flat(state.to_flat(), state)
        assert reconstructed.history_depth == 5

    def test_from_flat_wrong_size_raises(self) -> None:
        m = _slab_mesh()
        state = m.zeros_timed_full_field()
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
        np.testing.assert_array_equal(zeros.bulk.values, 0.0)
        np.testing.assert_array_equal(zeros.boundary.values, 0.0)
        # L2 norm reads the flat representation directly.
        l2 = _l2_norm(state)
        expected = float(np.linalg.norm(state.to_flat()))
        np.testing.assert_allclose(l2, expected, rtol=1e-15)
