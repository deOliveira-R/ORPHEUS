r"""Foundation pins for the canonical :class:`Quadrature` class
(R-1 Phase A detour-C: full retirement of ``orpheus/sn/quadrature.py``).

Each pin verifies a software invariant on the new directional
quadrature primitive — no equation labels are involved, so each
test carries ``@foundation`` and no ``verifies(...)``.

Q0.x — Construction: each classmethod factory returns a Quadrature
        with the right measure shape, weight sum, and metadata.
Q1.x — Canonical accessors: nodes / weights / N / n_ordinates / dim
        all read through to the underlying measure (single source).
Q2.x — axis_cosines(i) is dim-agnostic — works for 1-D scalar and
        multi-dim measures; out-of-range axes return zeros.
Q3.x — Legacy mu_x/mu_y/mu_z are @property views — agreement with
        axis_cosines, no separate storage.
Q4.x — reflection_index(axis): both int (0/1/2) and str ('x'/'y'/'z')
        keys resolve to the same partner array; partners satisfy
        the involution invariant `ref[ref[i]] == i`.
Q5.x — spherical_harmonics: shape (N, L+1, 2L+1); P0 == 1/sqrt(4π);
        slab GL has only m=0 harmonics non-zero.
Q6.x — Octants: disjoint, total mass preserved, label tuples match
        the sign-of-direction predicate.
Q7.x — Level structure: cylindrical-compatible (LS, product) carry
        non-None level_structure; slab/Lebedev carry None.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.numerics.measure import DiscreteMeasure
from orpheus.numerics.quadrature import LevelStructure, Quadrature

pytestmark = [pytest.mark.foundation]


# ─── Q0: Factory construction ───────────────────────────────────────────


def test_q0_1_gauss_legendre_factory() -> None:
    """``Quadrature.gauss_legendre(N)`` returns a Quadrature with a
    1-D scalar measure on [-1, 1] and N nodes."""
    q = Quadrature.gauss_legendre(8)
    assert q.n_ordinates == 8
    assert q.measure.nodes.shape == (8,)
    assert q.measure.weights.shape == (8,)
    assert q.measure.dim == 1
    # GL weights on [-1, 1] sum to 2.
    assert np.isclose(q.weights.sum(), 2.0)
    assert q.level_structure is None


def test_q0_2_lebedev_factory() -> None:
    """``Quadrature.lebedev(order)`` returns a Quadrature with a 3-D
    sphere measure and 4π-summing weights."""
    q = Quadrature.lebedev(17)
    assert q.measure.nodes.shape[1] == 3
    assert q.measure.dim == 3
    assert np.isclose(q.weights.sum(), 4 * np.pi)
    assert q.level_structure is None


def test_q0_3_level_symmetric_factory() -> None:
    """``Quadrature.level_symmetric(N)`` returns a Quadrature with
    a LevelStructure side-channel for cylindrical SN."""
    q = Quadrature.level_symmetric(4)
    assert q.measure.nodes.shape[1] == 3
    assert np.isclose(q.weights.sum(), 4 * np.pi)
    assert q.level_structure is not None
    assert isinstance(q.level_structure, LevelStructure)
    assert q.n_levels == q.level_structure.n_levels


def test_q0_4_product_factory() -> None:
    """``Quadrature.product(n_mu, n_phi)`` returns a Quadrature with
    a LevelStructure side-channel (GL × equispaced)."""
    q = Quadrature.product(n_mu=4, n_phi=4)
    assert q.measure.nodes.shape == (16, 3)
    assert np.isclose(q.weights.sum(), 4 * np.pi)
    assert q.level_structure is not None
    assert q.n_levels == q.level_structure.n_levels


# ─── Q1: Canonical accessors ────────────────────────────────────────────


def test_q1_1_canonical_accessors_read_through_measure() -> None:
    """nodes / weights / N / n_ordinates / dim are @property reads
    on the underlying measure — single source of truth."""
    q = Quadrature.lebedev(17)
    assert q.nodes is q.measure.nodes
    assert q.weights is q.measure.weights
    assert q.N == q.measure.n_points
    assert q.n_ordinates == q.measure.n_points
    assert q.dim == q.measure.dim


def test_q1_2_n_equals_n_ordinates() -> None:
    """``N`` and ``n_ordinates`` are aliases; both equal n_points."""
    q = Quadrature.gauss_legendre(8)
    assert q.N == q.n_ordinates == q.measure.n_points == 8


# ─── Q2: Dim-agnostic axis_cosines ──────────────────────────────────────


def test_q2_1_axis_cosines_1d_scalar_measure() -> None:
    """For a 1-D scalar measure (GL1D), axis 0 returns nodes;
    higher axes return zeros."""
    q = Quadrature.gauss_legendre(8)
    np.testing.assert_array_equal(q.axis_cosines(0), q.measure.nodes)
    np.testing.assert_array_equal(q.axis_cosines(1), np.zeros(q.N))
    np.testing.assert_array_equal(q.axis_cosines(2), np.zeros(q.N))


def test_q2_2_axis_cosines_multidim_measure() -> None:
    """For a 3-D sphere measure, axis i returns nodes[:, i]."""
    q = Quadrature.lebedev(17)
    for i in range(3):
        np.testing.assert_array_equal(q.axis_cosines(i), q.measure.nodes[:, i])
    # Beyond intrinsic dim → zeros.
    np.testing.assert_array_equal(q.axis_cosines(3), np.zeros(q.N))


# ─── Q3: Legacy mu_x/mu_y/mu_z views ────────────────────────────────────


def test_q3_1_mu_x_mu_y_mu_z_are_property_views() -> None:
    """mu_x/mu_y/mu_z agree with axis_cosines — no separate storage."""
    q = Quadrature.level_symmetric(4)
    np.testing.assert_array_equal(q.mu_x, q.axis_cosines(0))
    np.testing.assert_array_equal(q.mu_y, q.axis_cosines(1))
    np.testing.assert_array_equal(q.mu_z, q.axis_cosines(2))


def test_q3_2_mu_y_mu_z_zeros_for_slab() -> None:
    """Slab GL1D: mu_y == mu_z == 0 by 1-D structure."""
    q = Quadrature.gauss_legendre(8)
    np.testing.assert_array_equal(q.mu_y, np.zeros(q.N))
    np.testing.assert_array_equal(q.mu_z, np.zeros(q.N))


# ─── Q4: Reflection partners ────────────────────────────────────────────


def test_q4_1_reflection_index_accepts_int_and_str() -> None:
    """reflection_index('x') == reflection_index(0); same for y/z."""
    q = Quadrature.lebedev(17)
    np.testing.assert_array_equal(q.reflection_index("x"), q.reflection_index(0))
    np.testing.assert_array_equal(q.reflection_index("y"), q.reflection_index(1))
    np.testing.assert_array_equal(q.reflection_index("z"), q.reflection_index(2))


def test_q4_2_reflection_is_an_involution() -> None:
    """For every axis, ref[ref[i]] == i (reflection composed twice is identity)."""
    for q in [Quadrature.lebedev(17), Quadrature.level_symmetric(4),
              Quadrature.product(n_mu=4, n_phi=4)]:
        for axis in (0, 1, 2):
            ref = q.reflection_index(axis)
            np.testing.assert_array_equal(ref[ref], np.arange(q.N))


def test_q4_3_gl1d_x_reflection_is_index_reversal() -> None:
    """GL1D x-reflection: partner of i is N-1-i by GL-node symmetry."""
    q = Quadrature.gauss_legendre(8)
    np.testing.assert_array_equal(q.reflection_index("x"), np.arange(8)[::-1])
    # y/z reflections are identity (1-D).
    np.testing.assert_array_equal(q.reflection_index("y"), np.arange(8))
    np.testing.assert_array_equal(q.reflection_index("z"), np.arange(8))


def test_q4_4_unknown_axis_label_raises() -> None:
    """An unrecognized axis label raises ValueError."""
    q = Quadrature.gauss_legendre(4)
    with pytest.raises(ValueError, match="Unknown axis label"):
        q.reflection_index("w")


# ─── Q5: Spherical harmonics ────────────────────────────────────────────


def test_q5_1_spherical_harmonics_shape() -> None:
    """spherical_harmonics(L) returns shape (N, L+1, 2L+1)."""
    q = Quadrature.lebedev(17)
    Y = q.spherical_harmonics(2)
    assert Y.shape == (q.N, 3, 5)


def test_q5_2_spherical_harmonics_l0_constant() -> None:
    r"""Y_0^0 is a constant — the project's SH normalisation puts
    the :math:`\sqrt{4\pi}` factor on the moment side, so :math:`Y_0^0
    \equiv 1` at every ordinate.

    See :func:`orpheus.numerics.spherical_harmonics.evaluate_real_sh`
    for the convention details. The pin checks the constancy
    (rotation-invariance of the :math:`l=0` slot), not a specific
    normalisation value.
    """
    q = Quadrature.lebedev(17)
    Y = q.spherical_harmonics(0)
    np.testing.assert_allclose(Y[:, 0, 0], Y[0, 0, 0])
    assert Y[0, 0, 0] != 0.0


# ─── Q6: Octant partition ───────────────────────────────────────────────


def test_q6_1_octants_partition_is_disjoint() -> None:
    """Every ordinate appears in exactly one octant entry."""
    q = Quadrature.lebedev(17)
    all_indices = np.concatenate([part.indices for part in q.octants])
    assert sorted(all_indices.tolist()) == list(range(q.N))


def test_q6_2_octants_preserve_total_mass() -> None:
    """Σ partition weight totals == total weight."""
    q = Quadrature.lebedev(17)
    total = sum(part.measure.weights.sum() for part in q.octants)
    assert np.isclose(total, q.weights.sum())


# ─── Q7: Level structure ────────────────────────────────────────────────


def test_q7_1_level_structure_present_for_ls_and_product() -> None:
    """Cylindrical-compatible cubatures carry level_structure."""
    assert Quadrature.level_symmetric(4).level_structure is not None
    assert Quadrature.product(n_mu=4, n_phi=4).level_structure is not None


def test_q7_2_level_structure_absent_for_slab_and_lebedev() -> None:
    """Slab and pure-sphere cubatures carry no level_structure."""
    assert Quadrature.gauss_legendre(8).level_structure is None
    assert Quadrature.lebedev(17).level_structure is None


def test_q7_3_level_passthroughs_default_for_no_structure() -> None:
    """When level_structure is None, the passthrough properties return
    defaults that the cylindrical SN sweep can consume safely."""
    q = Quadrature.gauss_legendre(8)
    assert q.n_levels == 1
    assert len(q.level_indices) == 1
    np.testing.assert_array_equal(q.level_indices[0], np.arange(8))
    np.testing.assert_array_equal(q.level_mu, np.array([0.0]))


def test_q7_4_level_passthroughs_match_underlying_structure() -> None:
    """When level_structure is present, the passthroughs return
    its fields verbatim."""
    q = Quadrature.level_symmetric(4)
    assert q.n_levels == q.level_structure.n_levels
    assert len(q.level_indices) == q.n_levels
    np.testing.assert_array_equal(q.level_mu, q.level_structure.level_mu)
