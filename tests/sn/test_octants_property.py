r"""Tests for :attr:`AngularQuadrature.octants` cached property.

The :attr:`~orpheus.sn.quadrature._OctantsMixin.octants` property is the
SN-side adapter over
:meth:`orpheus.numerics.measure.DiscreteMeasure.partition_by` with the
sign-of-direction predicate. It realises the direct-sum decomposition

.. math::

    \mu_{S^2} \;=\; \bigoplus_\sigma \mu_\sigma,
    \qquad \sigma = (\mathrm{sign}\,\mu_x, \mathrm{sign}\,\mu_y,
                     \mathrm{sign}\,\mu_z)

on each :class:`AngularQuadrature` adapter. This test pins the
contract that the Wave-2 ``SweepDependencyGraph`` will consume:

* Disjoint coverage — every ordinate appears in exactly one entry.
* Sign correctness — for every entry with label ``σ``, every ordinate
  in ``entry.indices`` carries direction cosines whose signs match
  ``σ`` (component-wise).
* Mass conservation — the sum of partition weights equals the parent
  measure's total mass.
* Pure-axis ordinates form their own partition entry with one or more
  zero label components (verified explicitly on Lebedev 5).
* Caching identity — the property is computed once and re-used.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.numerics.measure import DiscreteMeasurePartition
from orpheus.sn.quadrature import (
    GaussLegendre1D,
    LebedevSphere,
    LevelSymmetricSN,
    ProductQuadrature,
    _OCTANT_SIGN_EPS,
    _octant_sign_predicate,
)


# Threshold the property uses to label "pure-axis" ordinates.
EPS = _OCTANT_SIGN_EPS


# ─────────────────────────────────────────────────────────────────────
# Predicate-only tests — pure helper, no quadrature in the loop
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.l0
class TestSignPredicate:
    """L0: ``_octant_sign_predicate`` returns the right shape and signs."""

    def test_1d_input(self):
        nodes = np.array([0.7, -0.3, 0.0, 1e-16, -1e-16, -1.0])
        labels = _octant_sign_predicate(nodes)
        np.testing.assert_array_equal(labels, [+1, -1, 0, 0, 0, -1])
        assert labels.shape == nodes.shape

    def test_3d_input(self):
        nodes = np.array([
            [+0.6, -0.4, +0.5],
            [-0.6, +0.4, -0.5],
            [0.0, 0.0, +1.0],
            [0.0, 0.0, -1.0],
            [+0.7, +0.7, 0.0],
        ])
        labels = _octant_sign_predicate(nodes)
        np.testing.assert_array_equal(labels, [
            [+1, -1, +1],
            [-1, +1, -1],
            [0, 0, +1],
            [0, 0, -1],
            [+1, +1, 0],
        ])
        assert labels.shape == nodes.shape

    def test_eps_threshold_is_strict(self):
        """Components within `EPS` of zero get label 0; outside it, ±1."""
        nodes = np.array([+EPS, -EPS, +EPS * 1.5, -EPS * 1.5])
        labels = _octant_sign_predicate(nodes)
        # |x| <= eps -> 0; |x| > eps -> sign(x)
        np.testing.assert_array_equal(labels, [0, 0, +1, -1])


# ─────────────────────────────────────────────────────────────────────
# Cross-adapter parametrised invariants
# ─────────────────────────────────────────────────────────────────────


_QUAD_FACTORIES = [
    pytest.param(lambda: GaussLegendre1D.create(8), id="gauss_legendre_1d_n8"),
    pytest.param(lambda: GaussLegendre1D.create(16), id="gauss_legendre_1d_n16"),
    pytest.param(lambda: LebedevSphere.create(5), id="lebedev_5"),
    pytest.param(lambda: LebedevSphere.create(7), id="lebedev_7"),
    pytest.param(lambda: LebedevSphere.create(11), id="lebedev_11"),
    pytest.param(lambda: LevelSymmetricSN.create(4), id="ls4"),
    pytest.param(lambda: LevelSymmetricSN.create(6), id="ls6"),
    pytest.param(lambda: ProductQuadrature.create(8, 8), id="product_8x8"),
]


@pytest.mark.l0
class TestOctantsInvariants:
    """L0: the four invariants every adapter's ``octants`` must satisfy."""

    @pytest.mark.parametrize("factory", _QUAD_FACTORIES)
    def test_returns_tuple_of_partitions(self, factory):
        quad = factory()
        octs = quad.octants
        assert isinstance(octs, tuple)
        assert len(octs) >= 2  # at least + and - hemisphere
        for entry in octs:
            assert isinstance(entry, DiscreteMeasurePartition)
            assert isinstance(entry.label, tuple)

    @pytest.mark.parametrize("factory", _QUAD_FACTORIES)
    def test_disjoint_and_complete_coverage(self, factory):
        """Every ordinate index appears in exactly one entry's `indices`."""
        quad = factory()
        all_indices = np.concatenate([e.indices for e in quad.octants])
        # Disjoint: no duplicates.
        assert len(np.unique(all_indices)) == len(all_indices)
        # Complete: union = {0, ..., N-1}.
        np.testing.assert_array_equal(
            np.sort(all_indices), np.arange(quad.N),
        )

    @pytest.mark.parametrize("factory", _QUAD_FACTORIES)
    def test_label_matches_direction_signs(self, factory):
        """For each entry with label σ, every ordinate's signs equal σ."""
        quad = factory()
        for entry in quad.octants:
            indices = entry.indices
            label = entry.label
            # Reconstruct signs from the cached mu_x/mu_y/mu_z arrays.
            if hasattr(quad, "mu_z"):
                # 3-D adapter: label is (s_x, s_y, s_z).
                actual = np.column_stack([
                    _octant_sign_predicate(quad.mu_x[indices]),
                    _octant_sign_predicate(quad.mu_y[indices]),
                    _octant_sign_predicate(quad.mu_z[indices]),
                ])
                expected = np.broadcast_to(
                    np.array(label, dtype=int), actual.shape,
                )
            else:
                # 1-D adapter: label is (s,) — only μ_x carries the sign.
                actual = _octant_sign_predicate(quad.mu_x[indices])
                expected = np.full(actual.shape, label[0], dtype=int)
            np.testing.assert_array_equal(actual, expected)

    @pytest.mark.parametrize("factory", _QUAD_FACTORIES)
    def test_weight_conservation(self, factory):
        """Sum of partition weights equals total parent mass.

        Tolerance: ``rtol=1e-14``. Each partition's ``weights.sum()`` and
        the parent's ``weights.sum()`` reduce in different orders, so
        IEEE-754 non-associativity introduces FP drift bounded by
        ``N × ULP × max_weight``. For Lebedev 11 (N=50, max_weight≈0.25)
        this is ≈ 2.8e-15 absolute / 2.2e-16 relative — well within
        ``rtol=1e-14``. Per ``vv-principles`` Bit-identity vs
        principled-equivalence: same value in real arithmetic, drift
        bounded by reduction depth, no algorithmic content lost.
        """
        quad = factory()
        partition_total = sum(
            float(e.measure.weights.sum()) for e in quad.octants
        )
        np.testing.assert_allclose(
            partition_total, float(quad.weights.sum()),
            rtol=1e-14, atol=0,
        )

    @pytest.mark.parametrize("factory", _QUAD_FACTORIES)
    def test_caching_identity(self, factory):
        """Two reads of the property return the SAME tuple object."""
        quad = factory()
        first = quad.octants
        second = quad.octants
        assert first is second


# ─────────────────────────────────────────────────────────────────────
# Per-adapter shape expectations
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.l0
class TestGaussLegendre1DOctants:
    """1-D: only ``(+1,)`` and ``(-1,)`` (GL has no μ=0 node)."""

    def test_two_partitions_only(self):
        quad = GaussLegendre1D.create(8)
        octs = quad.octants
        labels = {e.label for e in octs}
        assert labels == {(-1,), (+1,)}

    def test_label_is_one_tuple(self):
        quad = GaussLegendre1D.create(16)
        for e in quad.octants:
            assert len(e.label) == 1


@pytest.mark.l0
class TestLebedev5Octants:
    """Lebedev 5 has axis-aligned ``(0, 0, ±1)`` ordinates: 14 → 8 + 6 split."""

    def test_includes_pure_z_partitions(self):
        quad = LebedevSphere.create(5)
        labels = {e.label for e in quad.octants}
        # The two pure-z ordinates form their own partitions.
        assert (0, 0, +1) in labels
        assert (0, 0, -1) in labels

    def test_pure_z_partitions_have_one_ordinate_each(self):
        quad = LebedevSphere.create(5)
        per_label = {e.label: e for e in quad.octants}
        assert per_label[(0, 0, +1)].indices.size == 1
        assert per_label[(0, 0, -1)].indices.size == 1

    def test_full_octant_partitions_present(self):
        """The 8 full octants ``(±1, ±1, ±1)`` must all be present."""
        quad = LebedevSphere.create(5)
        labels = {e.label for e in quad.octants}
        for sx in (-1, +1):
            for sy in (-1, +1):
                for sz in (-1, +1):
                    assert (sx, sy, sz) in labels


@pytest.mark.l0
class TestLevelSymmetricSNOctants:
    """LS rules: only the 8 full octants (no axis-aligned ordinates)."""

    @pytest.mark.parametrize("order", [4, 6, 8])
    def test_eight_full_octants_only(self, order):
        quad = LevelSymmetricSN.create(order)
        octs = quad.octants
        # All 8 ``(±1, ±1, ±1)`` must be present.
        assert len(octs) == 8
        labels = {e.label for e in octs}
        for sx in (-1, +1):
            for sy in (-1, +1):
                for sz in (-1, +1):
                    assert (sx, sy, sz) in labels

    def test_partitions_have_equal_size(self):
        """LS quadrature is octant-symmetric: every partition has N/8 ordinates."""
        quad = LevelSymmetricSN.create(4)
        sizes = sorted(e.indices.size for e in quad.octants)
        assert sizes == [quad.N // 8] * 8


@pytest.mark.l0
class TestProtocolConformance:
    """The Protocol declares ``octants``; every adapter satisfies it."""

    @pytest.mark.parametrize("factory", _QUAD_FACTORIES)
    def test_protocol_has_octants(self, factory):
        quad = factory()
        # Duck-typed: accessing ``quad.octants`` succeeds; type-hint check
        # is implicit via the runtime_checkable Protocol.
        assert quad.octants is not None
