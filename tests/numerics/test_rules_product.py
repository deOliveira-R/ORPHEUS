"""Foundation + L1 tests for ``orpheus.numerics.quadrature.rules_product``.

The product rule :math:`\\mu_{\\text{GL}} \\times \\phi_{\\text{equispaced}}`
is the cylindrical-SN workhorse. Tests pin shape / metadata, weight
sum, the pass-through contract of
:meth:`~orpheus.numerics.quadrature.Quadrature.product`, and
polynomial-exactness in the polar factor. The per-family adapter
class the pass-through gate once compared against
(``orpheus.sn.quadrature.ProductQuadrature``) was retired into that
classmethod factory in R-1 Phase A detour-C.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.numerics.measure import SPACE_SPHERE, DiscreteMeasure
from orpheus.numerics.quadrature import product_mu_phi
from orpheus.numerics.quadrature.rules_sphere import LevelStructure
from orpheus.numerics.quadrature import Quadrature
from orpheus.numerics.symmetry import SubgroupOfO3


# ---------------------------------------------------------------------------
# Foundation
# ---------------------------------------------------------------------------


@pytest.mark.foundation
@pytest.mark.parametrize(
    ("n_mu", "n_phi"), [(4, 4), (8, 8), (8, 16), (16, 8), (12, 24)]
)
def test_product_returns_measure_and_structure(n_mu: int, n_phi: int) -> None:
    """Shape, metadata, accompanying :class:`LevelStructure`."""
    m, s = product_mu_phi(n_mu, n_phi)
    assert isinstance(m, DiscreteMeasure)
    assert m.nodes.ndim == 2
    assert m.nodes.shape == (n_mu * n_phi, 3)
    assert m.weights.shape == (n_mu * n_phi,)
    assert m.support == SPACE_SPHERE
    # D_{n_phi h}, not SO(2): the phi grid is a FINITE set, so no
    # continuous rotation preserves it. The rule carried SO2 until
    # 2026-08-02 — a claim no finite point set on S^2 can satisfy.
    assert m.invariance_group == SubgroupOfO3.Dnh(n_phi)
    assert m.degree_of_exactness == min(2 * n_mu - 1, n_phi - 1)

    assert isinstance(s, LevelStructure)
    assert s.n_levels == n_mu
    assert len(s.level_indices) == n_mu


@pytest.mark.foundation
@pytest.mark.parametrize(("n_mu", "n_phi"), [(4, 4), (8, 8), (16, 12)])
def test_product_weight_sum_is_4pi(n_mu: int, n_phi: int) -> None:
    """Weight sum :math:`= 4\\pi`."""
    m, _ = product_mu_phi(n_mu, n_phi)
    assert m.weights.sum() == pytest.approx(4.0 * np.pi, abs=1e-12)


@pytest.mark.foundation
@pytest.mark.parametrize(
    ("n_mu", "n_phi"), [(4, 4), (8, 8), (8, 16), (16, 8), (12, 24)]
)
def test_product_bit_identical_to_legacy_adapter(n_mu: int, n_phi: int) -> None:
    """:meth:`~orpheus.numerics.quadrature.Quadrature.product` passes
    this rule through unmodified — the accessor views read the right
    node columns and the level structure survives the wrap.

    A real contract, but NARROWER than the test name (and the retired
    ``ProductQuadrature.create`` comparison it was written against)
    suggests: the factory *calls* ``product_mu_phi``, so both sides
    move together by construction and this can NEVER detect node
    drift. Replacing ``product_mu_phi``'s body with random nodes
    leaves every assertion below GREEN.

    ⚠ And no pre-carve anchor exists elsewhere either. An earlier
    version of this docstring named the cylindrical regression
    snapshots (``cyl_*_product_*.npz``) as what "pins the ordinate
    order against the pre-carve behaviour". That is backwards: those
    snapshots were THEMSELVES re-captured by the consolidation
    (``81689a58``, "re-capture the SN baselines the quadrature
    consolidation moved"), so they pin POST-carve behaviour as a
    regression floor going forward. Nothing in the tree witnesses
    that the carve preserved the ordinate order."""
    m, s = product_mu_phi(n_mu, n_phi)
    legacy = Quadrature.product(n_mu=n_mu, n_phi=n_phi)
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
def test_product_rejects_zero() -> None:
    """``n_mu < 1`` or ``n_phi < 1`` must raise."""
    with pytest.raises(ValueError):
        product_mu_phi(0, 8)
    with pytest.raises(ValueError):
        product_mu_phi(8, 0)


# ---------------------------------------------------------------------------
# L1 — polynomial-exactness in the polar factor
# ---------------------------------------------------------------------------


@pytest.mark.l1
@pytest.mark.parametrize("n_mu", [4, 8, 16])
def test_product_polar_polynomial_exactness(n_mu: int) -> None:
    """The polar (Gauss-Legendre) factor of the product rule
    integrates :math:`\\mu^k` exactly for :math:`k \\le 2 n_\\mu - 1`,
    after multiplying out the trivial azimuthal integral
    :math:`\\int_0^{2\\pi} d\\phi = 2\\pi`.

    Closed-form: :math:`\\int_{S^2} \\mu_z^k \\, d\\Omega = 2\\pi \\cdot
    \\int_{-1}^{1} \\mu^k \\, d\\mu` which is :math:`0` for odd ``k``
    and :math:`4\\pi/(k+1)` for even ``k``.
    """
    n_phi = 8  # azimuthal exactness comfortable
    m, _ = product_mu_phi(n_mu, n_phi)
    for k in range(2 * n_mu):
        if k % 2 == 1:
            exact = 0.0
        else:
            exact = 4.0 * np.pi / (k + 1)
        approx = m.integrate(lambda x, k=k: x[:, 2] ** k)
        assert approx == pytest.approx(exact, abs=1e-12), (
            f"product({n_mu},{n_phi}) fails at polar degree {k}: "
            f"{approx} vs {exact}"
        )


# ---------------------------------------------------------------------------
# E4 — the spherical product theorem
# ---------------------------------------------------------------------------
#
# The degree used to be a hand-written ``min(2*n_mu - 1, n_phi - 1)`` over two
# bare integers. These pin it as a THEOREM over the two factors' typed claims:
# the value is unchanged, but it can no longer disagree with the factors, and
# the systems are checked rather than assumed.


@pytest.mark.foundation
@pytest.mark.parametrize(
    ("n_mu", "n_phi"), [(1, 1), (4, 4), (4, 8), (8, 4), (6, 12), (12, 24)]
)
def test_e4_the_product_claim_is_the_min_of_the_factors(
    n_mu: int, n_phi: int
) -> None:
    """The shipped claim equals the theorem applied to the factors' own
    claims — not to two integers that happen to match them today."""
    from orpheus.numerics.quadrature.rules_circle import (
        NODE_ALIGNED,
        periodic_trapezoid,
    )
    from orpheus.numerics.quadrature.rules_1d import gauss_legendre_on_mu
    from orpheus.numerics.quadrature.rules_product import (
        spherical_product_claim,
    )

    polar = gauss_legendre_on_mu(n_mu).exactness
    azimuthal = periodic_trapezoid(n_phi, shift=NODE_ALIGNED).exactness
    assert polar is not None and azimuthal is not None

    shipped = product_mu_phi(n_mu, n_phi)[0].exactness
    assert shipped == spherical_product_claim(polar, azimuthal)
    assert shipped is not None
    assert shipped.degree == min(2 * n_mu - 1, n_phi - 1)  # unchanged value


@pytest.mark.foundation
def test_e4_the_result_names_the_sphere_and_its_harmonics() -> None:
    """All three parts of the claim, since the whole point is that the
    integer alone was ambiguous between three different statements."""
    from orpheus.numerics.exactness import (
        UNIFORM_ON_SPHERE,
        OrthogonalSystem,
    )

    claim = product_mu_phi(4, 8)[0].exactness
    assert claim is not None
    assert claim.system is OrthogonalSystem.SPHERICAL_HARMONIC
    assert claim.reference == UNIFORM_ON_SPHERE
    assert claim.degree == 7


@pytest.mark.foundation
def test_e4_the_factors_do_NOT_compose_by_tensor_with() -> None:
    r"""Why the theorem is a separate function rather than
    ``ExactnessClaim.tensor_with``.

    A tensor product lands on the SQUARE :math:`[-1,1] \times S^1` and
    spans a mixed tensor system, so ``tensor_with`` returns ``None`` —
    correctly. This rule composes its factors through the EMBEDDING
    :math:`(\mu, \varphi) \mapsto \Omega` instead, so its claim is about
    :math:`S^2`. Different space, different theorem.
    """
    from orpheus.numerics.quadrature.rules_circle import (
        NODE_ALIGNED,
        periodic_trapezoid,
    )
    from orpheus.numerics.quadrature.rules_1d import gauss_legendre_on_mu

    polar = gauss_legendre_on_mu(4).exactness
    azimuthal = periodic_trapezoid(8, shift=NODE_ALIGNED).exactness
    assert polar is not None and azimuthal is not None
    assert polar.tensor_with(azimuthal) is None
    assert azimuthal.tensor_with(polar) is None


@pytest.mark.foundation
def test_e4_the_theorem_refuses_the_wrong_systems() -> None:
    """Two ALGEBRAIC claims are the square's tensor product, not this
    embedding — and a bare-integer ``min`` could not tell them apart.

    The azimuthal leg is the exact shape of the bug this carve exists to
    remove: ``equispaced`` on :math:`[0, 2\\pi]` carries the same nodes
    with an ALGEBRAIC degree of 1, and feeding it here would report a
    sphere degree of 1.
    """
    from orpheus.numerics.exactness import ExactnessClaim, UniformMeasure
    from orpheus.numerics.exactness import OrthogonalSystem as OS
    from orpheus.numerics.generating_measure import LEGENDRE
    from orpheus.numerics.measure import SPACE_CIRCLE, equispaced
    from orpheus.numerics.quadrature.rules_product import (
        spherical_product_claim,
    )

    algebraic = ExactnessClaim(LEGENDRE, 7)
    trigonometric = ExactnessClaim(
        UniformMeasure(SPACE_CIRCLE, OS.TRIGONOMETRIC), 7
    )

    with pytest.raises(ValueError, match="TRIGONOMETRIC"):
        spherical_product_claim(algebraic, algebraic)
    with pytest.raises(ValueError, match="ALGEBRAIC"):
        spherical_product_claim(trigonometric, trigonometric)

    # ...and the interval rule really does carry the degree-1 claim that
    # would have been silently accepted before the systems were typed.
    interval = equispaced(0.0, 2 * np.pi, 8).exactness
    assert interval is not None and interval.degree == 1
    with pytest.raises(ValueError, match="TRIGONOMETRIC"):
        spherical_product_claim(algebraic, interval)


@pytest.mark.foundation
@pytest.mark.parametrize("n_mu,n_phi", [(4, 4), (4, 8), (6, 12), (4, 5)])
def test_e4_sigma_has_one_size_now(n_mu: int, n_phi: int) -> None:
    r"""The mis-count closed.

    Before the azimuthal substitution the node set met the
    :math:`\xi = 0` axis twice per level but only ONE of the two was on
    it in exact arithmetic — ``np.sin(np.pi)`` is ``1.22e-16`` — so
    :math:`|\Sigma|` was 4 by equality and 8 by any tolerance. A fold
    whose well-posedness condition is :math:`\Sigma = \emptyset` cannot
    be decided on a set whose size depends on how you measure it.
    """
    xi = product_mu_phi(n_mu, n_phi)[0].nodes[:, 1]
    by_equality = int((xi == 0.0).sum())
    by_tolerance = int((np.abs(xi) < 1e-12).sum())
    assert by_equality == by_tolerance
    # Two axis crossings per level at even n_phi, one at odd.
    assert by_equality == n_mu * (2 if n_phi % 2 == 0 else 1)


@pytest.mark.foundation
@pytest.mark.parametrize("n_phi", [4, 8, 12, 24, 32])
def test_e4_the_level_order_tie_break_is_NAMED_not_round_off(
    n_phi: int,
) -> None:
    r"""η is 2-to-1 on the fiber, and with exact azimuths the ties are
    real.

    :math:`\eta = \sin\theta\cos\varphi` and
    :math:`\cos\varphi_m = \cos\varphi_{n-m}` holds bit-exactly for roots
    of unity, so a level carries only :math:`\lfloor n_\varphi/2 \rfloor
    + 1` distinct η values. Under the ``linspace``+``cos`` azimuths this
    rule used until 2026-08-02, round-off manufactured :math:`n_\varphi`
    fake distinctions and the intra-pair order was decided by noise.

    The order is increasing :math:`\varphi` within an η-tie, which this
    asserts by index.

    ⛔ **This gate was DEMOTED on 2026-08-13 and the demotion is
    correct** — read before crediting it with more than it has. It was
    written against ``kind="stable"``, and it had real teeth: `[M]`
    numpy's default (unstable) ``argsort`` AGREES with stable at
    :math:`n_\varphi \in \{4, 8, 12\}` (small arrays fall to insertion
    sort) and first DIVERGES at 24, which is why ``n_phi >= 24`` is in
    the parameterisation and why the first draft of this test, without
    it, asserted nothing. There is now no ``kind=`` to remove: the
    producer routes through
    :meth:`~orpheus.numerics.quadrature.rules_sphere.LevelStructure.from_level_membership`,
    whose key :math:`(\eta, \varphi, \operatorname{sign}\mu_z)` is
    injective, so no tie survives for any algorithm to break.

    The cost of that fix is this gate's teeth, and it is unavoidable
    rather than an oversight: ``np.lexsort`` is itself stable, so
    dropping a key component falls back to CONSTRUCTION order — and on
    this producer construction order within an η-tie already IS
    increasing :math:`\varphi` (the ``for m in range(n_phi)`` loop). So
    the order here is now over-determined three ways and `[M]` NO
    in-class mutation of the key can redden this: dropping
    :math:`\varphi`, dropping :math:`\operatorname{sign}\mu_z`, or
    dropping both leaves the stored order bit-identical at
    ``product(4,8)``, ``(4,24)``, ``(6,32)``, ``folded(4,16)``,
    ``(4,32)``.

    What it still tests, and what nothing else does: that **the ties are
    REAL** (the ``unique`` count — a claim about the rule's exact
    azimuths, which would fail again under a ``linspace``+``cos``
    regression) and that the level is η-sorted. The φ CONVENTION is
    gated where it is observable, on the other producer:
    ``tests/numerics/test_rules_sphere.py::
    test_the_azimuth_component_of_the_key_is_LOAD_BEARING_not_decorative``.
    """
    _, structure = product_mu_phi(2, n_phi)
    level = structure.level_indices[0]
    eta = product_mu_phi(2, n_phi)[0].nodes[level, 0]

    assert len(np.unique(eta)) == n_phi // 2 + 1  # the ties are REAL
    assert np.all(np.diff(eta) >= 0.0)            # ...and still sorted

    # Within each tie, the lower construction index (smaller φ) is first.
    ties = [i for i in range(len(level) - 1) if eta[i] == eta[i + 1]]
    assert ties, "no η-tie in this level — the gate would be vacuous"
    for i in ties:
        assert level[i] < level[i + 1]


# ---------------------------------------------------------------------------
# spherical_product — the combinator seam (campaign Q5.2)
# ---------------------------------------------------------------------------
#
# `product_mu_phi` is now a named family OF the combinator: the pair of
# factor rules arrives as ARGUMENTS, so the fold selects the staggered
# circle rule by substitution — no `half_range=True` flag exists. The
# result's group is DERIVED from three RigidMotion.preserves generator
# checks on the factors (this module shipped three false symmetry
# declarations: ERR-072/073/074), and the Hasse walk over the assembled
# nodes is the independent realization that confirms the derivation.


@pytest.mark.foundation
def test_the_named_family_is_a_spherical_product_of_its_two_rules() -> None:
    """Pass-through smoke: the wrapper IS the combinator at NODE_ALIGNED.

    Both sides move together (the wrapper calls the combinator), so this
    pins the API relationship, NOT the carve — the carve's independence
    pin was the pre-carve capture, bit-identical on four configs
    ([M] 2026-08-07, commit message of the Q5.2 landing).
    """
    from orpheus.numerics.quadrature import (
        NODE_ALIGNED,
        gauss_legendre_on_mu,
        periodic_trapezoid,
        spherical_product,
    )

    via_wrapper, s_wrapper = product_mu_phi(4, 8)
    via_combinator, s_combinator = spherical_product(
        gauss_legendre_on_mu(4), periodic_trapezoid(8, shift=NODE_ALIGNED)
    )
    np.testing.assert_array_equal(via_wrapper.nodes, via_combinator.nodes)
    np.testing.assert_array_equal(via_wrapper.weights, via_combinator.weights)
    assert via_wrapper.invariance_group == via_combinator.invariance_group
    assert via_wrapper.exactness == via_combinator.exactness
    for a, b in zip(s_wrapper.level_indices, s_combinator.level_indices):
        np.testing.assert_array_equal(a, b)


@pytest.mark.foundation
def test_the_seam_selects_the_staggered_fiber_and_the_fold_composes() -> None:
    r"""The Q5.2 -> Q5.1 handshake, with no flag anywhere.

    Substituting the STAGGERED circle rule empties the singular set
    (:math:`\Sigma = \emptyset`, T24's admissibility condition — by
    integer equality, since roots-of-unity azimuths make on-axis sines
    exactly 0.0), and the fold then composes straight through
    ``DiscreteMeasure.quotient``: 32 -> 16 at uniform doubled weights.
    NODE_ALIGNED, by contrast, puts two nodes per level ON the mirror.
    """
    from orpheus.numerics.quadrature import (
        NODE_ALIGNED,
        STAGGERED,
        gauss_legendre_on_mu,
        periodic_trapezoid,
        spherical_product,
    )
    from orpheus.numerics.symmetry import singular_set

    mirror = SubgroupOfO3.Mirror("y")
    polar = gauss_legendre_on_mu(4)

    aligned, _ = spherical_product(polar, periodic_trapezoid(8, shift=NODE_ALIGNED))
    assert singular_set(aligned, mirror).size == 2 * 4  # phi=0 and pi per level

    staggered, _ = spherical_product(polar, periodic_trapezoid(8, shift=STAGGERED))
    assert singular_set(staggered, mirror).size == 0  # free — the fold's condition

    folded = staggered.quotient(mirror)
    assert folded.n_points == staggered.n_points // 2
    # Weight multiset: every orbit is free, so each folded atom carries
    # exactly twice one member's weight. `[::2]` is a valid one-member-
    # per-orbit selection ONLY because weights are constant within a
    # level and levels are contiguous blocks in construction order.
    np.testing.assert_array_equal(
        np.sort(folded.weights), np.sort(2.0 * staggered.weights[::2])
    )
    np.testing.assert_array_almost_equal_nulp(
        np.array(folded.weights.sum()), np.array(staggered.weights.sum()), nulp=1
    )


@pytest.mark.foundation
def test_the_derived_group_and_the_hasse_walk_verify_each_other() -> None:
    """Two independent realizations of Sym(product): the construction
    derives D_nh from three generator checks on the FACTORS; the walk
    computes it from the ASSEMBLED nodes. Neither can be wrong alone
    without this reddening."""
    from orpheus.numerics.quadrature import (
        STAGGERED,
        gauss_legendre_on_mu,
        periodic_trapezoid,
        spherical_product,
    )
    from orpheus.numerics.symmetry import maximal_invariance_groups

    for n_mu, n_phi, shift in [(2, 4, None), (4, 8, STAGGERED)]:
        if shift is None:
            measure, _ = product_mu_phi(n_mu, n_phi)
        else:
            measure, _ = spherical_product(
                gauss_legendre_on_mu(n_mu),
                periodic_trapezoid(n_phi, shift=shift),
            )
        assert measure.invariance_group == SubgroupOfO3.Dnh(n_phi)
        walk = maximal_invariance_groups(measure)
        assert SubgroupOfO3.Dnh(n_phi) in walk, (n_mu, n_phi, walk)


@pytest.mark.foundation
def test_a_factor_pair_with_an_unrealized_group_family_is_refused() -> None:
    """An asymmetric polar factor would make the group C_nv, which
    ``SubgroupOfO3`` does not realize — refused naming the family,
    never tagged wrong-or-weaker (the field means MAXIMAL). Positive
    control: the symmetric pair constructs."""
    from orpheus.numerics.quadrature import (
        NODE_ALIGNED,
        gauss_legendre_on_mu,
        periodic_trapezoid,
        spherical_product,
    )

    azimuthal = periodic_trapezoid(8, shift=NODE_ALIGNED)
    lopsided_polar = DiscreteMeasure(
        nodes=np.array([-0.5, 0.2, 0.7]),
        weights=np.array([0.6, 0.7, 0.7]),
        support="[-1,1]",
        exactness=gauss_legendre_on_mu(3).exactness,  # borrowed, valid-typed
    )
    with pytest.raises(ValueError, match=r"C_8v.*does not realize"):
        spherical_product(lopsided_polar, azimuthal)

    measure, _ = spherical_product(gauss_legendre_on_mu(3), azimuthal)
    assert measure.invariance_group == SubgroupOfO3.Dnh(8)


@pytest.mark.foundation
def test_an_interval_rule_as_the_fiber_is_refused_by_its_space() -> None:
    """E3's same-nodes/wrong-space confusion, unrepresentable end to
    end: an interval rule carrying the azimuthal ANGLES is a different
    object from the circle rule carrying the POINTS, and the support
    tag is what discriminates them here."""
    from orpheus.numerics.measure import equispaced
    from orpheus.numerics.quadrature import gauss_legendre_on_mu, spherical_product

    angles_as_interval_rule = equispaced(0.0, 2.0 * np.pi, 8)
    with pytest.raises(ValueError, match="must be a circle rule"):
        spherical_product(gauss_legendre_on_mu(4), angles_as_interval_rule)


@pytest.mark.foundation
def test_off_circle_fiber_nodes_are_refused() -> None:
    """Non-unit circle nodes would embed to direction cosines with
    ||Omega|| != 1 — refused at the door."""
    from orpheus.numerics.measure import SPACE_CIRCLE
    from orpheus.numerics.quadrature import (
        NODE_ALIGNED,
        gauss_legendre_on_mu,
        periodic_trapezoid,
        spherical_product,
    )

    honest = periodic_trapezoid(8, shift=NODE_ALIGNED)
    inflated = DiscreteMeasure(
        nodes=1.1 * honest.nodes,
        weights=honest.weights,
        support=SPACE_CIRCLE,
        exactness=honest.exactness,
    )
    with pytest.raises(ValueError, match="unit"):
        spherical_product(gauss_legendre_on_mu(4), inflated)
