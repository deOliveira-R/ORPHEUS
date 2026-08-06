"""Foundation + L1 tests for ``orpheus.numerics.quadrature.rules_sphere``.

Lebedev rules (Lebedev 1976) and Carlson-Lathrop level-symmetric
:math:`S_N` rules (Carlson & Lathrop 1968) are both
:math:`O_h`-invariant by construction. Foundation tests pin shape /
metadata + the pass-through contract of
:meth:`~orpheus.numerics.quadrature.Quadrature.lebedev` and
:meth:`~orpheus.numerics.quadrature.Quadrature.level_symmetric`; L1
tests check polynomial-exactness for small Lebedev orders against
closed-form spherical integrals. The per-family SN adapters those
pass-through gates once compared against (``LebedevSphere`` /
``LevelSymmetricSN`` at ``orpheus.sn.quadrature``) were retired into
those classmethod factories in R-1 Phase A detour-C.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.numerics.measure import SPACE_SPHERE, DiscreteMeasure
from orpheus.numerics.quadrature import lebedev_sphere, level_symmetric_sn
from orpheus.numerics.quadrature import product_mu_phi
from orpheus.numerics.quadrature.rules_sphere import LevelStructure
from orpheus.numerics.quadrature import Quadrature
from orpheus.numerics.symmetry import SubgroupOfO3


# ---------------------------------------------------------------------------
# Lebedev — foundation
# ---------------------------------------------------------------------------


@pytest.mark.foundation
@pytest.mark.parametrize("order", [3, 5, 7, 11, 17, 23])
def test_lebedev_returns_discrete_measure(order: int) -> None:
    """Shape, space, ``invariance_group`` / ``degree_of_exactness``."""
    m = lebedev_sphere(order)
    assert isinstance(m, DiscreteMeasure)
    assert m.nodes.ndim == 2
    assert m.nodes.shape[1] == 3  # (N, 3)
    assert m.weights.shape == (m.n_points,)
    assert m.support == SPACE_SPHERE
    assert m.invariance_group == SubgroupOfO3.OctahedralOh
    assert m.degree_of_exactness == order


@pytest.mark.foundation
@pytest.mark.parametrize("order", [3, 5, 7, 11, 17])
def test_lebedev_weight_sum_is_4pi(order: int) -> None:
    """Weight sum must equal :math:`4\\pi`."""
    m = lebedev_sphere(order)
    assert m.weights.sum() == pytest.approx(4.0 * np.pi, abs=1e-13)


@pytest.mark.foundation
@pytest.mark.parametrize("order", [3, 5, 7, 11, 17])
def test_lebedev_bit_identical_to_legacy_adapter(order: int) -> None:
    """:meth:`~orpheus.numerics.quadrature.Quadrature.lebedev` passes
    this rule through unmodified — the ``mu_x`` / ``mu_y`` / ``mu_z``
    views read the right node columns of the wrapped measure.

    NARROWER than the test name (and the retired ``LebedevSphere``
    comparison it was written against) suggests: the factory *calls*
    ``lebedev_sphere``, so both sides move together by construction
    and this can NEVER detect node drift."""
    m = lebedev_sphere(order)
    legacy = Quadrature.lebedev(order)
    nodes = m.nodes
    assert np.array_equal(legacy.mu_x, nodes[:, 0])
    assert np.array_equal(legacy.mu_y, nodes[:, 1])
    assert np.array_equal(legacy.mu_z, nodes[:, 2])
    assert np.array_equal(legacy.weights, m.weights)


# ---------------------------------------------------------------------------
# Lebedev — L1 polynomial-exactness
# ---------------------------------------------------------------------------


@pytest.mark.l1
@pytest.mark.parametrize("order", [5, 7, 11])
def test_lebedev_integrates_constants(order: int) -> None:
    """:math:`\\int_{S^2} 1 \\, d\\Omega = 4\\pi` — the normalisation
    invariant. This is the degree-0 polynomial exactness."""
    m = lebedev_sphere(order)
    one = m.integrate(lambda x: np.ones(x.shape[0]))
    assert one == pytest.approx(4.0 * np.pi, abs=1e-13)


@pytest.mark.l1
def test_lebedev_integrates_quadratic_invariant() -> None:
    r"""For an :math:`O_h`-invariant rule of order :math:`\ge 2`,
    :math:`\int_{S^2} \mu_x^2 \, d\Omega = \int \mu_y^2 = \int \mu_z^2 =
    \frac{1}{3} \int |\Omega|^2 \, d\Omega = 4\pi/3`. The three
    components must match by symmetry."""
    m = lebedev_sphere(5)  # degree-5 exact
    val_x = m.integrate(lambda x: x[:, 0] ** 2)
    val_y = m.integrate(lambda x: x[:, 1] ** 2)
    val_z = m.integrate(lambda x: x[:, 2] ** 2)
    expected = 4.0 * np.pi / 3.0
    assert val_x == pytest.approx(expected, abs=1e-12)
    assert val_y == pytest.approx(expected, abs=1e-12)
    assert val_z == pytest.approx(expected, abs=1e-12)


# ---------------------------------------------------------------------------
# Level-symmetric S_N — foundation
# ---------------------------------------------------------------------------


@pytest.mark.foundation
@pytest.mark.parametrize("sn_order", [2, 4, 6, 8])
def test_level_symmetric_returns_measure_and_structure(sn_order: int) -> None:
    """Shape, space, ``invariance_group``, accompanying
    :class:`LevelStructure`."""
    m, s = level_symmetric_sn(sn_order)
    assert isinstance(m, DiscreteMeasure)
    assert m.nodes.ndim == 2
    assert m.nodes.shape[1] == 3
    assert m.support == SPACE_SPHERE
    assert m.invariance_group == SubgroupOfO3.OctahedralOh
    # ``max(3, N-1)``, not ``N-1``, since #327 gave the rule per-O_h-orbit
    # moment-matched weights: `[M]` S2 reaches degree 3, because its 8 nodes
    # are a SINGLE orbit whose one free weight is already fixed by
    # ``Σw = 4π`` — the rule over-delivers against the formula there and the
    # old ``N-1`` UNDER-claimed it. That under-claim is what identified the
    # tag as a formula for a different construction rather than a
    # mis-measurement of this one; the value is measured against the
    # closed-form monomial integral in
    # ``tests/numerics/test_advertised_degree_is_measured.py``.
    assert m.degree_of_exactness == max(3, sn_order - 1)

    assert isinstance(s, LevelStructure)
    assert s.n_levels == sn_order // 2
    assert len(s.level_indices) == sn_order // 2
    # Total ordinate count = sum of per-level counts (in the upper
    # hemisphere) × 2 (lower hemisphere mirror).
    levels_total = sum(len(li) for li in s.level_indices)
    assert levels_total == m.n_points


@pytest.mark.foundation
@pytest.mark.parametrize("sn_order", [2, 4, 6, 8])
def test_level_symmetric_weight_sum_is_4pi(sn_order: int) -> None:
    """Weight sum must equal :math:`4\\pi`."""
    m, _ = level_symmetric_sn(sn_order)
    assert m.weights.sum() == pytest.approx(4.0 * np.pi, abs=1e-13)


@pytest.mark.foundation
@pytest.mark.parametrize("sn_order", [2, 4, 6, 8])
def test_level_symmetric_bit_identical_to_legacy_adapter(sn_order: int) -> None:
    """:meth:`~orpheus.numerics.quadrature.Quadrature.level_symmetric`
    passes this rule through unmodified — the accessor views read the
    right node columns and the per-level grouping survives the wrap.

    NARROWER than the test name (and the retired ``LevelSymmetricSN``
    comparison it was written against) suggests: the factory *calls*
    ``level_symmetric_sn``, so both sides move together by
    construction and this can NEVER detect node drift. Node order and
    per-level grouping matter — if either drifts, the cylindrical
    sweep's α-recursion delivers different per-cell flux.

    ⚠ Nothing witnesses that the carve PRESERVED them. An earlier
    version of this docstring named the cylindrical regression
    snapshots (``cyl_*_LS4_*.npz``) as what "pins them against the
    pre-carve behaviour"; that is backwards, because both LS4
    snapshots were themselves re-captured by the consolidation
    (``81689a58``). They are a POST-carve regression floor going
    forward, not a witness to the carve."""
    m, s = level_symmetric_sn(sn_order)
    legacy = Quadrature.level_symmetric(sn_order)
    nodes = m.nodes
    assert np.array_equal(legacy.mu_x, nodes[:, 0])
    assert np.array_equal(legacy.mu_y, nodes[:, 1])
    assert np.array_equal(legacy.mu_z, nodes[:, 2])
    assert np.array_equal(legacy.weights, m.weights)
    assert legacy.n_levels == s.n_levels
    assert np.array_equal(legacy.level_mu, s.level_mu)
    for li_legacy, li_new in zip(legacy.level_indices, s.level_indices):
        assert np.array_equal(li_legacy, li_new)


@pytest.mark.foundation
def test_level_symmetric_rejects_odd_or_zero() -> None:
    """``sn_order`` must be a positive even integer."""
    with pytest.raises(ValueError):
        level_symmetric_sn(0)
    with pytest.raises(ValueError):
        level_symmetric_sn(3)


# ---------------------------------------------------------------------------
# Level-symmetric S_N — L1
# ---------------------------------------------------------------------------


@pytest.mark.l1
@pytest.mark.parametrize("sn_order", [4, 8])
def test_level_symmetric_integrates_constants(sn_order: int) -> None:
    """:math:`\\int_{S^2} 1 \\, d\\Omega = 4\\pi` — moment-0
    invariant satisfied by all level-symmetric rules."""
    m, _ = level_symmetric_sn(sn_order)
    one = m.integrate(lambda x: np.ones(x.shape[0]))
    assert one == pytest.approx(4.0 * np.pi, abs=1e-13)


@pytest.mark.l1
def test_level_symmetric_quadratic_isotropy() -> None:
    """For an :math:`O_h`-invariant rule, the second-moment integrals
    of the three Cartesian cosines must agree (each is :math:`4\\pi/3`).
    """
    m, _ = level_symmetric_sn(8)
    val_x = m.integrate(lambda x: x[:, 0] ** 2)
    val_y = m.integrate(lambda x: x[:, 1] ** 2)
    val_z = m.integrate(lambda x: x[:, 2] ** 2)
    expected = 4.0 * np.pi / 3.0
    # S_8 reaches degree 5 in moment-condition rules; 4π/3 is degree 2.
    assert val_x == pytest.approx(expected, rel=1e-10)
    assert val_y == pytest.approx(expected, rel=1e-10)
    assert val_z == pytest.approx(expected, rel=1e-10)


# ---------------------------------------------------------------------------
# Octahedral invariance is EXACT, not approximate
# ---------------------------------------------------------------------------
#
# A rule that advertises O_h invariance should realize it as an integer
# permutation of ordinate indices. That is what a reflecting boundary
# condition needs: the face reflection must map ordinate m onto ordinate
# m', not merely near it.


def _oh_exactness(nodes: "np.ndarray") -> "tuple[float, int, int]":
    """(worst landing distance, #operators exact, #operators) over O_h."""
    from orpheus.numerics.symmetry import SubgroupOfO3, _realized_ops

    ops = _realized_ops(SubgroupOfO3.OctahedralOh._tag)
    assert ops is not None and len(ops) == 48
    worst, exact = 0.0, 0
    for g in ops:
        # The element applies ITSELF to the points. The `nodes @ M.T` this
        # replaced put the row-vs-column convention at the call site, which
        # is one of the conventions the carve exists to seat in one place.
        moved = g.on_points(nodes)
        landing = np.linalg.norm(
            moved[:, None, :] - nodes[None, :, :], axis=2
        ).min(axis=1).max()
        worst = max(worst, float(landing))
        exact += int(landing == 0.0)
    return worst, exact, len(ops)


@pytest.mark.foundation
# S12 caps every level-symmetric sweep in this module since #327: the
# moment-matched solve has no POSITIVE solution above it (`[M]` min weight
# -0.027 at S14, -0.142 at S20), so the builder refuses and the ORDERS
# SIMPLY DO NOT EXIST. These rows gate NODE structure, which is unchanged
# by #327 -- the coverage is not weakened, its subject is gone. If a future
# node choice pushes the positivity frontier up, restore the orders here.
@pytest.mark.parametrize("sn_order", [4, 6, 8, 10, 12])
def test_level_symmetric_is_EXACTLY_octahedral(sn_order: int) -> None:
    r"""All 48 :math:`O_h` operators map the node set onto itself bit-exactly.

    Regression for a defect that made the rule advertise :math:`O_h` and
    realize only :math:`D_{2h}`. The construction used to recover the
    third direction cosine numerically as
    :math:`\sqrt{1 - \mu_z^2 - \eta^2}`, landing within ~1e-16 of the
    level value but not on it. `[M]` At :math:`N = 16` the y-axis then
    carried **22** distinct magnitudes where the level array has 8, and
    only the 8 pure sign flips were exact — negation is exact in IEEE,
    a coordinate swap only if the values match.

    The third level index is fixed by the triangular relation
    :math:`p + k + j = N/2 - 1`, so it is index arithmetic, not
    arithmetic.
    """
    measure, _ = level_symmetric_sn(sn_order)
    worst, exact, total = _oh_exactness(measure.nodes)
    assert exact == total, (
        f"S_{sn_order}: only {exact}/{total} O_h operators are exact "
        f"(worst landing {worst:.3e}) — the rule realizes a subgroup of "
        f"what it advertises"
    )


@pytest.mark.foundation
@pytest.mark.parametrize("sn_order", [4, 8, 10, 12])
def test_level_symmetric_axes_share_one_magnitude_set(sn_order: int) -> None:
    """The defining property of a *level*-symmetric rule.

    Every ordinate's three cosines are drawn from the SAME level array,
    so the three axes must carry identical magnitude sets — as sets of
    exact floats, not merely to a tolerance. This is what makes the
    axis-permuting half of :math:`O_h` a relabeling.
    """
    measure, structure = level_symmetric_sn(sn_order)
    per_axis = [set(np.abs(measure.nodes[:, k]).tolist()) for k in range(3)]
    assert per_axis[0] == per_axis[1] == per_axis[2], (
        f"S_{sn_order}: axes carry "
        f"{[len(s) for s in per_axis]} distinct magnitudes"
    )
    # And that shared set IS the level array.
    assert per_axis[0] == set(np.abs(structure.level_mu).tolist())


@pytest.mark.foundation
@pytest.mark.parametrize("sn_order", [4, 6, 8, 10, 12])
def test_lebedev_is_EXACTLY_octahedral(sn_order: int) -> None:
    """Lebedev already had this property — it is the reference case.

    Its tabulated orbits are sign-flip and coordinate-permutation images
    of a few representatives, so closure is bit-exact by construction.
    Pinned so a future table edit cannot quietly lose it.
    """
    order = {4: 3, 6: 5, 8: 7, 10: 9, 12: 11}[sn_order]
    measure = lebedev_sphere(order=order)
    worst, exact, total = _oh_exactness(measure.nodes)
    assert exact == total, (
        f"lebedev({order}): {exact}/{total} exact, worst {worst:.3e}"
    )


# ---------------------------------------------------------------------------
# A level is a FIBER of an invariant
# ---------------------------------------------------------------------------


@pytest.mark.foundation
def test_the_two_producers_fiber_over_different_invariants() -> None:
    """One type, two invariants — which is why the type must say which.

    ``product_mu_phi`` gives every Gauss-Legendre polar node its own
    level, so levels run over all of ``[-1, 1]`` and a level's fiber is
    ONE circle. ``level_symmetric_sn`` indexes levels by ``|mu_z|``, so
    each carries BOTH hemispheres and the fiber is two circles.

    A consumer reading ``level_indices`` without knowing which is
    reading an ambiguous object.
    """
    from orpheus.numerics.quadrature.rules_sphere import PolarInvariant

    _, prod = product_mu_phi(n_mu=4, n_phi=8)
    _, lsn = level_symmetric_sn(8)
    assert prod.polar_invariant is PolarInvariant.SIGNED_MU_Z
    assert lsn.polar_invariant is PolarInvariant.ABS_MU_Z

    # And the consequence is observable, not merely declared: signed
    # levels span both signs across levels; abs levels do so WITHIN one.
    m_prod, _ = product_mu_phi(n_mu=4, n_phi=8)
    m_lsn, _ = level_symmetric_sn(8)
    assert len(np.unique(np.sign(m_prod.nodes[prod.level_indices[0], 2]))) == 1
    assert len(np.unique(np.sign(m_lsn.nodes[lsn.level_indices[0], 2]))) == 2


@pytest.mark.foundation
@pytest.mark.parametrize(
    ("label", "build"),
    [
        ("product(4,8)", lambda: product_mu_phi(n_mu=4, n_phi=8)),
        ("product(3,7)", lambda: product_mu_phi(n_mu=3, n_phi=7)),
        ("level_symmetric(8)", lambda: level_symmetric_sn(8)),
        ("level_symmetric(12)", lambda: level_symmetric_sn(12)),
    ],
)
def test_fiber_coordinate_is_injective_on_every_level(
    label: str, build
) -> None:
    """``(hemisphere, azimuth)`` separates the ordinates of a level.

    This is the property the stored ``level_indices`` order does NOT
    have: its key ``eta = sin(theta) cos(phi)`` is even in ``phi``,
    hence 2-to-1 on the fiber. An ordering by a non-injective key is
    not an ordering of the fiber at all — it is an ordering of the
    fiber modulo the mirror, which is the mechanism behind #326.
    """
    _, s = build()
    for level in range(s.n_levels):
        members = s.fiber(level)
        pairs = np.stack([s.hemisphere[members], s.azimuth[members]], axis=1)
        assert len(np.unique(pairs, axis=0)) == len(members), (
            f"{label} level {level}: fiber coordinate repeats"
        )


@pytest.mark.foundation
def test_projection_order_is_2_to_1_where_the_fiber_order_is_not() -> None:
    """The measurement that makes the previous test's point concrete."""
    _, s = product_mu_phi(n_mu=2, n_phi=8)
    m, _ = product_mu_phi(n_mu=2, n_phi=8)
    members = s.level_indices[0]
    eta = m.nodes[members, 0]

    # The stored key collapses 8 ordinates onto 5 values...
    assert len(np.unique(np.round(eta, 12))) < len(members)
    # ...and cannot recover the sign of xi, which the sweep needs.
    xi_signs = np.sign(m.nodes[members, 1])
    assert len(np.unique(xi_signs)) > 1

    # The fiber coordinate keeps all 8 apart.
    fib = s.fiber(0)
    assert len(np.unique(s.azimuth[fib])) == len(fib)


@pytest.mark.foundation
@pytest.mark.parametrize("sn_order", [4, 8, 12])
def test_level_membership_is_exact_equality(sn_order: int) -> None:
    """Level assignment is an equality, not a neighbourhood test.

    The 8-fold sign replication copies ``mu_z`` straight out of the
    level array, so ``|mu_z|`` IS the level value bit-for-bit. The
    construction carried a ``tol = 1e-12`` window until 2026-08-02 —
    a float comparison answering a question the generating loop had
    already answered exactly.
    """
    measure, structure = level_symmetric_sn(sn_order)
    mu_z = measure.nodes[:, 2]
    covered = np.zeros(len(mu_z), dtype=bool)
    for level in range(structure.n_levels):
        members = structure.level_indices[level]
        assert np.all(np.abs(mu_z[members]) == structure.level_mu[level])
        covered[members] = True
    assert covered.all(), "every ordinate belongs to exactly one level"
