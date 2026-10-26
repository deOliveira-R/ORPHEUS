"""The indicator basis keeps its partition under descending edges and reverses the cell order (the orientation trap).

The orientation trap (P5.0 memo FLAG 1)
---------------------------------------
``IndicatorBasis.evaluate`` buckets points with ``searchsorted(side=
"right")``, which assumes ASCENDING edges. The canonical ``eg`` is
DESCENDING. Feeding descending coarse edges to ``IndicatorBasis`` does
NOT break the partition (row-sums stay 1, because searchsorted+clip
always returns an in-range bin) — it silently REVERSES the coarse-group
order (fast fine groups land in the thermal coarse column). So the
``overlap_to`` law that MUST be pinned is the coarse-group ORDER
(fast-first), not merely partition completeness. ``TestOrientationTrap``
*demonstrates* the reversal directly against the live ``IndicatorBasis``
(it needs no SUT), so the implementer cannot ship a condensation whose
coarse groups are silently thermal-first.

These two rows exercise only :class:`~orpheus.numerics.basis.indicator_basis.IndicatorBasis` (numerics): they moved from
``tests/gates/data/test_energy_grid.py`` by the layer-DAG re-filing (Phase 1b, 2026-10-26), bodies unchanged, and the data-side law
they motivate (``EnergyGrid.overlap_to`` keeps the coarse groups fast-first) stays there as
``TestOverlapToIntrinsic::test_coarse_group_order_is_fast_first``. The descending grid ``_EG_DESC`` is defined here, the lower layer, and that file imports it.

vv Mode-8: every assertion is ``np.testing.*`` (the suite runs ``-O``; bare ``assert`` is stripped).
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry.manifold import HALF_LINE
from orpheus.numerics.basis.indicator_basis import IndicatorBasis

pytestmark = pytest.mark.foundation


# A small DESCENDING energy grid: 4 groups, fast-first (group 0 fastest).
_EG_DESC = np.array([1.0e7, 1.0e6, 1.0e3, 1.0e0, 1.0e-5])  # 5 edges → 4 groups


# ══════════════════════════════════════════════════════════════════════
# The orientation trap — demonstrated against the LIVE IndicatorBasis
# (needs NO SUT; pins the FLAG-1 silent coarse-group-order reversal)
# ══════════════════════════════════════════════════════════════════════


class TestOrientationTrap:
    """The descending-edge / ``searchsorted`` coarse-group-order reversal.

    These run unconditionally — they exercise the existing
    ``IndicatorBasis`` to prove WHY ``GroupCondensation`` must pin the
    coarse-group ORDER, not just partition completeness.
    """

    @staticmethod
    def _fine_points(eg_desc: np.ndarray) -> np.ndarray:
        """Geometric-midpoint representative energy of each fine group."""
        return np.sqrt(eg_desc[:-1] * eg_desc[1:])

    def test_descending_edges_keep_partition_but_reverse_order(self) -> None:
        """Descending coarse edges → valid partition BUT reversed coarse order.

        Fast fine groups (low fine-index, high energy) land in the HIGH
        coarse-column when edges are descending — the silent reversal the
        condensation MUST NOT inherit.
        """
        fine_pts = self._fine_points(_EG_DESC)            # 4 pts, descending energy
        # Coarse: fast {g0,g1}, thermal {g2,g3}. DESCENDING coarse edges.
        coarse_desc = np.array([1.0e7, 1.0e3, 1.0e-5])
        table = IndicatorBasis(
            edges_per_axis=(coarse_desc,), partition_of=HALF_LINE,
        ).evaluate(fine_pts)
        # Partition is "valid" (every fine group → exactly one coarse cell):
        np.testing.assert_array_equal(
            table.sum(axis=1), np.ones(4),
            err_msg="descending edges should still give a one-hot partition",
        )
        # ...but the column assignment is REVERSED: fast fine groups (0,1)
        # land in coarse column 1, thermal (2,3) in column 0.
        assignment = table.argmax(axis=1)
        np.testing.assert_array_equal(
            assignment, np.array([1, 1, 0, 0]),
            err_msg="descending edges silently reverse the coarse-group order",
        )

    def test_ascending_edges_required_for_fast_first(self) -> None:
        """Reversing edges to ASCENDING is what the searchsorted bucketing needs.

        With ascending coarse edges the fast fine groups map to the HIGH
        ascending-cell index; the implementer must reverse that cell axis
        back to recover fast-first coarse groups (or build the partition map
        explicitly, orientation-free).
        """
        fine_pts = self._fine_points(_EG_DESC)
        coarse_asc = np.array([1.0e-5, 1.0e3, 1.0e7])     # ascending
        table = IndicatorBasis(
            edges_per_axis=(coarse_asc,), partition_of=HALF_LINE,
        ).evaluate(fine_pts)
        np.testing.assert_array_equal(table.sum(axis=1), np.ones(4))
        # Fast fine groups (high energy) → high ascending-cell index (1).
        np.testing.assert_array_equal(
            table.argmax(axis=1), np.array([1, 1, 0, 0]),
            err_msg="ascending-edge cell index puts fast groups in the high cell",
        )
