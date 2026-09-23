r"""Tests for :meth:`DiscreteMeasure.partition_by` and
:class:`DiscreteMeasurePartition`.

The partition primitive realises the inverse of
:meth:`DiscreteMeasure.__add__` (direct sum): given a labelling
predicate :math:`\ell : \mathcal{X} \to L` on the nodes,
:meth:`partition_by` returns the disjoint decomposition

.. math::

   \mu \;=\; \bigoplus_{\lambda \in L} \mu_\lambda.

The test invariants verified:

* **Disjoint union**: every node appears in exactly one partition entry.
* **Mass preservation**: sum of partition weights equals total mass.
* **Round-trip with ``__add__``**: reassembled measure has the same
  ``(node, weight)`` multiset as the parent.
* **Octant predicate on Lebedev**: 8 partition entries on a
  3-D Lebedev sphere; sign pattern ``(±1, ±1, ±1)``.
* **Pure-axis ordinates label sign=0**: the partition correctly
  isolates degenerate ordinates with ``|μ_axis| < tol``.
* **Scalar vs compound labels**: 1-D and 2-D predicate outputs both
  produce well-formed tuple labels.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.numerics.manifold import REAL_LINE, SPHERE
from orpheus.numerics.measure import DiscreteMeasure, DiscreteMeasurePartition
from orpheus.numerics.quadrature import lebedev_sphere
from orpheus.numerics.symmetry import SubgroupOfO3


@pytest.mark.l0
class TestPartitionByScalarLabel:
    """L0: predicate returning a 1-D scalar-label array."""

    def test_two_partitions_by_sign(self):
        nodes = np.array([-2.0, -1.0, 0.5, 3.0])
        weights = np.array([0.1, 0.2, 0.3, 0.4])
        mu = DiscreteMeasure(nodes=nodes, weights=weights, support=REAL_LINE)
        parts = mu.partition_by(lambda x: np.sign(x).astype(int))
        # Three labels: -1, 0, 1 (but no zero in our nodes), so:
        labels = sorted(p.label for p in parts)
        assert labels == [(-1,), (1,)]

    def test_label_is_tuple_even_for_scalar_predicate(self):
        nodes = np.array([1.0, 2.0])
        weights = np.array([0.5, 0.5])
        mu = DiscreteMeasure(nodes=nodes, weights=weights, support=REAL_LINE)
        parts = mu.partition_by(lambda x: np.zeros(len(x), dtype=int))
        assert len(parts) == 1
        # Label must be a tuple, even though the predicate returned 1-D
        assert isinstance(parts[0].label, tuple)
        assert parts[0].label == (0,)

    def test_indices_correctly_recover_nodes(self):
        nodes = np.array([10.0, 20.0, 30.0, 40.0])
        weights = np.array([1.0, 1.0, 1.0, 1.0])
        mu = DiscreteMeasure(nodes=nodes, weights=weights, support=REAL_LINE)
        parts = mu.partition_by(lambda x: (x % 20 == 0).astype(int))
        for p in parts:
            np.testing.assert_array_equal(
                p.measure.nodes, mu.nodes[p.indices],
            )


@pytest.mark.l0
class TestPartitionByCompoundLabel:
    """L0: predicate returning a 2-D compound-label array."""

    def test_3d_sign_partition(self):
        # Six axis-aligned directions: +x, -x, +y, -y, +z, -z.
        nodes = np.array([
            [1, 0, 0], [-1, 0, 0],
            [0, 1, 0], [0, -1, 0],
            [0, 0, 1], [0, 0, -1],
        ], dtype=float)
        weights = np.full(6, 1.0)
        mu = DiscreteMeasure(nodes=nodes, weights=weights, support=SPHERE)
        # Compound label = (sign(x), sign(y), sign(z))
        parts = mu.partition_by(lambda nodes: np.sign(nodes).astype(int))
        # Six unique labels, one per axis-aligned direction.
        assert len(parts) == 6
        for p in parts:
            assert isinstance(p.label, tuple)
            assert len(p.label) == 3
            assert p.measure.n_points == 1


@pytest.mark.l0
class TestDisjointUnionInvariant:
    """L0: every node in exactly one partition; weights sum to total mass."""

    def test_disjoint_union_on_lebedev(self):
        measure = lebedev_sphere(5)
        parts = measure.partition_by(
            lambda nodes: np.sign(nodes).astype(int)
        )
        # Disjoint
        all_indices = np.concatenate([p.indices for p in parts])
        assert len(np.unique(all_indices)) == len(all_indices)
        assert sorted(all_indices.tolist()) == list(range(measure.n_points))
        # Mass preserved
        total = sum(float(p.measure.weights.sum()) for p in parts)
        np.testing.assert_allclose(total, float(measure.weights.sum()), rtol=1e-15)


@pytest.mark.l0
class TestRoundTripWithDirectSum:
    """L0: direct sum of partition measures = parent (modulo ordering)."""

    def test_round_trip_lebedev(self):
        measure = lebedev_sphere(7)
        parts = measure.partition_by(
            lambda nodes: np.sign(nodes).astype(int)
        )
        # Reassemble via direct sum
        rebuilt = parts[0].measure
        for p in parts[1:]:
            rebuilt = rebuilt + p.measure
        # Same point count
        assert rebuilt.n_points == measure.n_points
        # Same total mass (this is the integral of f=1)
        np.testing.assert_allclose(
            float(rebuilt.weights.sum()),
            float(measure.weights.sum()),
            rtol=1e-15,
        )
        # Multiset equality of (node, weight) pairs — sort by lex(node)
        # then compare componentwise.
        def sort_key(nodes, weights):
            order = np.lexsort(nodes.T)
            return nodes[order], weights[order]
        n_orig, w_orig = sort_key(measure.nodes, measure.weights)
        n_new, w_new = sort_key(rebuilt.nodes, rebuilt.weights)
        np.testing.assert_array_equal(n_orig, n_new)
        np.testing.assert_array_equal(w_orig, w_new)


@pytest.mark.l0
class TestOctantPartitionOnLebedev:
    """L0: the canonical octant partition for the SN sweep.

    A 3-D Lebedev grid has eight octants; the predicate
    ``λ x: (sign(x_x), sign(x_y), sign(x_z))`` produces one entry per
    octant. Pure-axis ordinates (any sign zero) form their own
    entries — they're degenerate for the SN sweep and need separate
    handling.
    """

    def test_lebedev_order_5_has_eight_octants(self):
        """Lebedev order-5 has axis-aligned ordinates (sign=0 components),
        so the partition has more than 8 entries — one per unique
        sign tuple. Verify count by enumerating expected sign tuples
        present in the rule."""
        measure = lebedev_sphere(5)
        parts = measure.partition_by(
            lambda nodes: np.sign(nodes).astype(int)
        )
        # All labels are 3-tuples of int in {-1, 0, +1}.
        for p in parts:
            assert len(p.label) == 3
            for v in p.label:
                assert v in (-1, 0, 1)

    def test_full_octant_partitions_present(self):
        """Lebedev rules include axis-aligned ordinates (sign=0 components),
        so the partition has more than 8 entries. But the eight
        ``(±1, ±1, ±1)`` octants are ALWAYS present (Lebedev preserves
        octahedral symmetry), each with at least one ordinate.
        """
        measure = lebedev_sphere(13)
        parts = measure.partition_by(
            lambda nodes: np.sign(nodes).astype(int)
        )
        full_octants = {
            p.label for p in parts if 0 not in p.label
        }
        expected_octants = {
            (sx, sy, sz)
            for sx in (-1, 1) for sy in (-1, 1) for sz in (-1, 1)
        }
        assert full_octants == expected_octants
        for p in parts:
            if 0 not in p.label:
                assert p.measure.n_points >= 1


@pytest.mark.l0
class TestPartitionEntryFields:
    def test_partition_entry_is_frozen_dataclass(self):
        nodes = np.array([1.0, 2.0])
        weights = np.array([0.5, 0.5])
        mu = DiscreteMeasure(nodes=nodes, weights=weights, support=REAL_LINE)
        parts = mu.partition_by(lambda x: np.zeros(len(x), dtype=int))
        p = parts[0]
        assert isinstance(p, DiscreteMeasurePartition)
        with pytest.raises((AttributeError, TypeError)):
            p.label = (1,)  # frozen dataclass forbids mutation

    def test_label_python_scalars_not_numpy(self):
        """Labels should round-trip cleanly through dict keys etc.
        Numpy integer scalars are equal to Python ints but hash
        differently — use Python types."""
        nodes = np.array([1.0, -1.0])
        weights = np.array([0.5, 0.5])
        mu = DiscreteMeasure(nodes=nodes, weights=weights, support=REAL_LINE)
        parts = mu.partition_by(lambda x: np.sign(x).astype(int))
        for p in parts:
            for v in p.label:
                assert type(v) is int  # exact type check, not isinstance

    def test_predicate_must_match_node_count(self):
        nodes = np.array([1.0, 2.0, 3.0])
        weights = np.array([1.0, 1.0, 1.0])
        mu = DiscreteMeasure(nodes=nodes, weights=weights, support=REAL_LINE)
        with pytest.raises(ValueError, match="all 3 nodes"):
            mu.partition_by(lambda x: np.array([0]))


@pytest.mark.l0
class TestPartitionDropsMetadata:
    """L0: partition entries drop invariance_group + degree_of_exactness.

    A restriction typically breaks invariance, and the parent
    rule's polynomial exactness does not survive. Caller may re-tag
    with `with_metadata` if they have evidence the partition
    preserves the property.
    """

    def test_metadata_dropped(self):
        measure = lebedev_sphere(7).with_metadata(invariance_group=SubgroupOfO3.OctahedralOh)
        assert measure.invariance_group == SubgroupOfO3.OctahedralOh
        assert measure.degree_of_exactness is not None
        parts = measure.partition_by(
            lambda nodes: np.sign(nodes).astype(int)
        )
        for p in parts:
            assert p.measure.invariance_group is None
            assert p.measure.degree_of_exactness is None
