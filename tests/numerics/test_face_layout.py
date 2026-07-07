r"""L0 — :class:`~orpheus.numerics.face_layout.FaceLayout` foundation tests.

The L1 primitive that powers post-D-G
:class:`~orpheus.transport.fields.angular_boundary_flux.AngularBoundaryFlux`. Pins:

* :class:`FaceSlot` invariants (``flat_size == prod(shape)``,
  ``offset >= 0``).
* :class:`FaceLayout` invariants (``total_size == sum(slot.flat_size)``,
  contiguous packing with no gap or overlap).
* :class:`FaceSlot.slice_view` returns a memory-shared view (no copy).
* :meth:`FaceLayout.from_named_shapes` produces contiguous offsets in
  iteration order.

References
----------

* Depth B plan §3.4 (Option Ω flat-buffer storage).
* Plan §6 step D-G prep.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.numerics.face_layout import FaceLayout, FaceSlot


pytestmark = [pytest.mark.foundation]


# ───────────────────────────────────────────────────────────────────────
# A. FaceSlot invariants
# ───────────────────────────────────────────────────────────────────────


class TestFaceSlotInvariants:
    def test_flat_size_matches_prod_shape(self) -> None:
        slot = FaceSlot(key="xmin", offset=0, shape=(4, 3), flat_size=12)
        assert slot.flat_size == 12

    def test_flat_size_mismatch_raises(self) -> None:
        with pytest.raises(ValueError, match="flat_size"):
            FaceSlot(key="xmin", offset=0, shape=(4, 3), flat_size=11)

    def test_negative_offset_raises(self) -> None:
        with pytest.raises(ValueError, match="offset"):
            FaceSlot(key="xmin", offset=-1, shape=(4,), flat_size=4)

    def test_frozen(self) -> None:
        from dataclasses import FrozenInstanceError
        slot = FaceSlot(key="xmin", offset=0, shape=(4,), flat_size=4)
        with pytest.raises(FrozenInstanceError):
            slot.offset = 1  # type: ignore[misc]


class TestFaceSlotSliceView:
    def test_slice_view_round_trip(self) -> None:
        flat = np.arange(20, dtype=float)
        slot = FaceSlot(key="middle", offset=4, shape=(2, 3), flat_size=6)
        view = slot.slice_view(flat)
        assert view.shape == (2, 3)
        np.testing.assert_array_equal(view, flat[4:10].reshape(2, 3))

    def test_slice_view_is_memory_shared(self) -> None:
        flat = np.zeros(10, dtype=float)
        slot = FaceSlot(key="seg", offset=2, shape=(4,), flat_size=4)
        view = slot.slice_view(flat)
        view[:] = 5.0
        np.testing.assert_array_equal(flat[2:6], np.full(4, 5.0))
        assert flat[0] == 0.0 and flat[6] == 0.0

    def test_slice_view_reshape_2d(self) -> None:
        flat = np.arange(24, dtype=float)
        slot = FaceSlot(key="square", offset=8, shape=(2, 2, 2), flat_size=8)
        view = slot.slice_view(flat)
        assert view.shape == (2, 2, 2)
        np.testing.assert_array_equal(view.ravel(), flat[8:16])


# ───────────────────────────────────────────────────────────────────────
# B. FaceLayout invariants
# ───────────────────────────────────────────────────────────────────────


class TestFaceLayoutInvariants:
    def test_total_size_matches_sum(self) -> None:
        faces = {
            "a": FaceSlot(key="a", offset=0, shape=(3,), flat_size=3),
            "b": FaceSlot(key="b", offset=3, shape=(4,), flat_size=4),
        }
        layout = FaceLayout(faces=faces, total_size=7)
        assert layout.total_size == 7

    def test_total_size_mismatch_raises(self) -> None:
        faces = {
            "a": FaceSlot(key="a", offset=0, shape=(3,), flat_size=3),
        }
        with pytest.raises(ValueError, match="total_size"):
            FaceLayout(faces=faces, total_size=5)

    def test_gap_raises(self) -> None:
        """Gap detection: total_size matches sum, but slots leave a hole."""
        faces = {
            "a": FaceSlot(key="a", offset=0, shape=(2,), flat_size=2),
            "b": FaceSlot(key="b", offset=5, shape=(3,), flat_size=3),  # GAP at 2..5
        }
        # Sum of flat_size = 5; pick total_size=5 so the sum check passes
        # and the gap/overlap walker fires instead.
        with pytest.raises(ValueError, match="gap or overlap"):
            FaceLayout(faces=faces, total_size=5)

    def test_overlap_raises(self) -> None:
        """Overlap detection: total_size matches sum, but slots overlap."""
        faces = {
            "a": FaceSlot(key="a", offset=0, shape=(4,), flat_size=4),
            "b": FaceSlot(key="b", offset=2, shape=(2,), flat_size=2),  # OVERLAP [2..4)
        }
        # Sum of flat_size = 6; pick total_size=6 so the sum check passes
        # and the walker detects slot b starting at 2 while cursor is at 4.
        with pytest.raises(ValueError, match="gap or overlap"):
            FaceLayout(faces=faces, total_size=6)


# ───────────────────────────────────────────────────────────────────────
# C. FaceLayout.from_named_shapes — the canonical constructor
# ───────────────────────────────────────────────────────────────────────


class TestFromNamedShapes:
    def test_two_face_slab_1g(self) -> None:
        layout = FaceLayout.from_named_shapes([
            ("xmin", (2, 1)),
            ("xmax", (2, 1)),
        ])
        assert layout.total_size == 4
        assert layout.faces["xmin"].offset == 0
        assert layout.faces["xmin"].flat_size == 2
        assert layout.faces["xmax"].offset == 2
        assert layout.faces["xmax"].flat_size == 2

    def test_one_face_sphere(self) -> None:
        layout = FaceLayout.from_named_shapes([("xmax", (4, 3))])
        assert layout.total_size == 12
        assert "xmin" not in layout.faces
        assert layout.faces["xmax"].offset == 0

    def test_four_face_2d_cartesian(self) -> None:
        # N=4, ng=2, nx=3, ny=2 — typical small 2-D mesh
        N, ng, nx, ny = 4, 2, 3, 2
        layout = FaceLayout.from_named_shapes([
            ("xmin", (N, ng, ny)),  # 16
            ("xmax", (N, ng, ny)),  # 16
            ("ymin", (N, ng, nx)),  # 24
            ("ymax", (N, ng, nx)),  # 24
        ])
        expected_total = N * ng * (2 * ny + 2 * nx)
        assert layout.total_size == expected_total == 80

        # Verify contiguous packing.
        assert layout.faces["xmin"].offset == 0
        assert layout.faces["xmax"].offset == 16
        assert layout.faces["ymin"].offset == 32
        assert layout.faces["ymax"].offset == 56

    def test_per_face_shapes_reachable_via_slice_view(self) -> None:
        N, ng, nx, ny = 4, 2, 3, 2
        layout = FaceLayout.from_named_shapes([
            ("xmin", (N, ng, ny)),
            ("xmax", (N, ng, ny)),
            ("ymin", (N, ng, nx)),
            ("ymax", (N, ng, nx)),
        ])
        flat = np.zeros(layout.total_size)

        xmin_view = layout.faces["xmin"].slice_view(flat)
        ymax_view = layout.faces["ymax"].slice_view(flat)
        assert xmin_view.shape == (N, ng, ny)
        assert ymax_view.shape == (N, ng, nx)

        xmin_view[:] = 1.0
        ymax_view[:] = 9.0
        # Pure writes through views must land in the right flat region.
        np.testing.assert_array_equal(flat[:16], np.full(16, 1.0))
        np.testing.assert_array_equal(flat[56:], np.full(24, 9.0))
        # Other regions untouched.
        np.testing.assert_array_equal(flat[16:56], np.zeros(40))
