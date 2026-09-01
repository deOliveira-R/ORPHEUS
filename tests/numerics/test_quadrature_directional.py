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
Q4.x — ordinate_permutation(σ_a): the derived mirror permutations are
        certified (involution, weights, the ACTUAL reflection) and exist
        for exactly the axes each rule is closed under.
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

from orpheus.geometry.transformation import RigidMotion
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


def test_q3_3_cylindrical_eta_xi_aliases_agree_with_mu_x_mu_y() -> None:
    r"""``eta`` / ``xi`` aliases on :class:`Quadrature` are the
    cylindrical-frame names for ``mu_x`` / ``mu_y`` — same data,
    honest naming.

    The slab convention names ``mu_x``/``mu_y`` mislead in
    cylindrical SN: column 0 of ``measure.nodes`` is the radial
    cosine :math:`\eta`, NOT a Cartesian-X projection; column 1
    is the azimuthal cosine :math:`\xi`. The aliases let
    cylindrical-context call sites read in the right frame.
    """
    q = Quadrature.level_symmetric(4)
    np.testing.assert_array_equal(q.eta, q.mu_x)
    np.testing.assert_array_equal(q.xi, q.mu_y)
    # And they round-trip through axis_cosines(i).
    np.testing.assert_array_equal(q.eta, q.axis_cosines(0))
    np.testing.assert_array_equal(q.xi, q.axis_cosines(1))


# ─── Q4: Mirror-induced ordinate permutations ───────────────────────────
#
# Until G6.3 §7d.3 these rows gated the precomputed ``reflection_index``
# table; they now gate the DERIVATION that replaced it
# (``ordinate_permutation``), keeping each row's own claim. The axis-letter
# sugar (Q4.1) and the axis-vocabulary refusal (Q4.4) died with the
# accessor — the letter→mirror spelling lives on the deck tier
# (``_mirror_motion``, gated in ``tests/geometry/test_paired_deck.py``).


def _mirror(axis_index: int) -> RigidMotion:
    """σ_axis — the mirror whose plane is normal to the named axis."""
    return RigidMotion.reflection(normal=np.eye(3)[axis_index])


def test_q4_2_mirror_permutation_is_an_involution() -> None:
    """For every axis, pi[pi] == arange (a mirror composed twice is identity).

    ERR-044's quadrature-tier home: ``preserves`` certifies bijection and
    weight equality, NOT the involution — σ² = id is asserted ON the
    derived π.
    """
    for q in [Quadrature.lebedev(17), Quadrature.level_symmetric(4),
              Quadrature.product(n_mu=4, n_phi=4)]:
        for axis in (0, 1, 2):
            pi = q.ordinate_permutation(_mirror(axis))
            assert pi is not None
            ref = pi.indices
            np.testing.assert_array_equal(ref[ref], np.arange(q.N))


def test_q4_3_gl1d_x_mirror_is_index_reversal() -> None:
    """GL1D σ_x pairs i with N-1-i by GL-node symmetry; σ_y/σ_z fix every
    ordinate (the 1-D embedding ``(μ, 0, 0)`` has zero y/z components).

    Until §7d.3 these three maps were a stored closed form on the factory;
    they are now DERIVED, and this row pins that the derivation reproduces
    the closed form.
    """
    q = Quadrature.gauss_legendre(8)
    x = q.ordinate_permutation(_mirror(0))
    assert x is not None
    np.testing.assert_array_equal(x.indices, np.arange(8)[::-1])
    for axis in (1, 2):
        pi = q.ordinate_permutation(_mirror(axis))
        assert pi is not None
        np.testing.assert_array_equal(pi.indices, np.arange(8))


# ─── Q4.5+: the derived permutation is CERTIFIED, not merely involutive ─
#
# Campaign Q5.0.1 (re-posed at §7d.3). ``test_q4_2`` above checks only
# ``ref[ref[i]] == i`` — which ``geometry/boundary/_specular.py`` documents
# as insufficient for exactly this object: "the involution does not [catch
# it] (it can still be its own inverse)". That is ERR-042, catalogued one
# layer up for specular BC pairings; the pre-Q5.0.1 partner table carried
# none of that module's three checks and reproduced the failure.


_SHIPPED_RULES = [
    ("lebedev(17)", lambda: Quadrature.lebedev(17), (0, 1, 2)),
    ("level_symmetric(4)", lambda: Quadrature.level_symmetric(4), (0, 1, 2)),
    ("level_symmetric(8)", lambda: Quadrature.level_symmetric(8), (0, 1, 2)),
    ("product(4,4)", lambda: Quadrature.product(n_mu=4, n_phi=4), (0, 1, 2)),
    ("product(4,8)", lambda: Quadrature.product(n_mu=4, n_phi=8), (0, 1, 2)),
    ("product(4,12)", lambda: Quadrature.product(n_mu=4, n_phi=12), (0, 1, 2)),
    # ODD n_phi belongs in this list even though σ_x is refused on it:
    # `[M]` with only even rows, restoring the pre-Q5.0.1 bare-argmin body
    # left this gate GREEN — every map it inspected happened to be correct.
    # The odd rows are what give it teeth against the defect it documents,
    # because there the legacy body produced an axis-0 map whose reflection
    # residual was 0.33-0.58.
    ("product(4,5)", lambda: Quadrature.product(n_mu=4, n_phi=5), (1, 2)),
    ("product(4,7)", lambda: Quadrature.product(n_mu=4, n_phi=7), (1, 2)),
    ("gauss_legendre(8)", lambda: Quadrature.gauss_legendre(8), (0, 1, 2)),
]


@pytest.mark.parametrize(
    "name,build,expected_axes", _SHIPPED_RULES, ids=[r[0] for r in _SHIPPED_RULES],
)
def test_q4_5_every_derived_mirror_permutation_really_is_its_reflection(
    name, build, expected_axes,
) -> None:
    r"""Where σ_a induces a permutation it satisfies all three pairing
    checks — and it EXISTS for exactly the axes the rule is closed under.

    The three are independent — each passes a map the other two reject
    (``_specular.py``'s ERR-042 / ERR-044 / ERR-045 table):

    * **involution** — ``ref[ref[i]] == i``;
    * **the actual reflection** — ``R x_i == x_{ref[i]}``, which is the
      one a garbage map fails and the involution cannot see (checked on
      the full ``(N, 3)`` embedding, so degenerate axes are exercised
      rather than skipped);
    * **measure preservation** — ``w_i == w_{ref[i]}``.

    The per-row ``expected_axes`` makes the availability claim explicit —
    a derivation that started answering ``None`` everywhere would red
    here instead of silently vacating the loop.
    """
    q = build()
    nodes3 = np.column_stack([q.axis_cosines(a) for a in range(3)])
    available = []
    for axis in (0, 1, 2):
        pi = q.ordinate_permutation(_mirror(axis))
        if pi is None:
            continue
        available.append(axis)
        ref = pi.indices
        np.testing.assert_array_equal(
            ref[ref], np.arange(q.N), err_msg=f"{name} axis {axis}: not an involution",
        )
        np.testing.assert_array_equal(
            q.weights[ref], q.weights,
            err_msg=f"{name} axis {axis}: partners carry different weights",
        )
        reflected = nodes3.copy()
        reflected[:, axis] *= -1.0
        residual = float(np.max(np.abs(reflected - nodes3[ref])))
        assert residual < 1e-11, (
            f"{name} axis {axis}: the derived map is NOT the reflection — "
            f"max |R x_i - x_ref[i]| = {residual:.4e}. An involutive but "
            f"wrong map is ERR-042's signature."
        )
    assert tuple(available) == expected_axes, (
        f"{name}: σ-mirror availability {tuple(available)} != expected "
        f"{expected_axes}"
    )


@pytest.mark.parametrize("n_phi", [4, 6, 8, 12, 16])
def test_q4_6_even_n_phi_products_keep_all_three_mirrors(n_phi: int) -> None:
    """The certification is a NO-OP on everything shipped with even ``n_phi``.

    Pins that Q5.0.1 tightened the contract without narrowing any rule the
    tree actually uses — and that the §7d.3 derivation preserves that.
    """
    q = Quadrature.product(n_mu=4, n_phi=n_phi)
    available = [
        a for a in (0, 1, 2) if q.ordinate_permutation(_mirror(a)) is not None
    ]
    assert available == [0, 1, 2]


@pytest.mark.parametrize("n_phi", [5, 7, 9])
def test_q4_7_odd_n_phi_product_has_NO_x_mirror(n_phi: int) -> None:
    r"""An odd-:math:`n_\varphi` product is not :math:`\sigma_x`-closed, so
    the derivation answers ``None`` for σ_x.

    Mechanism: a product rule's mirror planes sit at
    :math:`k\pi/n_\varphi`. :math:`\sigma_x` is :math:`\varphi \to \pi -
    \varphi`, i.e. a plane at :math:`\pi/2`, needing :math:`k =
    n_\varphi/2` — an integer only for even :math:`n_\varphi`.
    :math:`\sigma_y` is the :math:`k = 0` plane and therefore survives at
    every :math:`n_\varphi`, which is why the cylindrical fold is
    unconditional while the centreline map is not.

    `[M]` Before Q5.0.1 the stored axis-0 map was wrong by ``0.58 / 0.42 /
    0.33`` in the direction cosines at ``n_phi = 5 / 7 / 9`` — **and was
    still an involution**, so ``test_q4_2`` passed on it. It feeds the
    :math:`r=0` pole continuation.

    The RAISE this row used to pin ("no precomputed reflection partner" —
    the table's lookup miss) retired with the table at §7d.3. ``None`` is
    the derivation's honest answer here; the LOUD refusals are pinned at
    the consumer tiers where they now live — "no specular pairing" at
    ``realize()`` (``test_sn_boundary_realizer``) and "cannot seed the
    r = 0 pole" at the curvilinear sweep
    (``test_coupled_pole_mu_level_invariant``).
    """
    q = Quadrature.product(n_mu=4, n_phi=n_phi)
    assert q.ordinate_permutation(_mirror(0)) is None
    # sigma_y survives — the fold does not depend on the parity.
    assert q.ordinate_permutation(_mirror(1)) is not None


def test_q4_8_derivation_refuses_when_closure_actually_breaks() -> None:
    """Mutation control: perturb ONE node off the mirror and σ_y answers
    ``None``.

    Without this the parametrized rows above could all be passing because
    the certification never refuses anything.
    """
    good = Quadrature.lebedev(17)
    assert good.ordinate_permutation(_mirror(1)) is not None

    broken_nodes = good.measure.nodes.copy()
    broken_nodes[0, 1] += 0.05  # move one node off the y-mirror orbit
    broken = Quadrature(
        measure=DiscreteMeasure(
            nodes=broken_nodes,
            weights=good.measure.weights,
            support=good.measure.support,
        )
    )
    assert broken.ordinate_permutation(_mirror(1)) is None, (
        "perturbing a node off the y-mirror still yielded a σ_y permutation "
        "— the closure check is not actually checking closure"
    )


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

    See :mod:`orpheus.numerics.basis.spherical_harmonic_basis` — the
    canonical home of the :math:`Y_\ell^m` evaluator and of the discrete
    Gram, which absorbed the free ``evaluate_real_sh`` function this
    docstring used to cite — for the convention details. The pin checks the constancy
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
    structure = q.level_structure
    assert structure is not None
    assert q.n_levels == structure.n_levels
    assert len(q.level_indices) == q.n_levels
    np.testing.assert_array_equal(q.level_mu, structure.level_mu)


# ─── Q8: the angular frame's measure is ROUTED, not rebuilt ─────────────
#
# Phase 0.1a/0.1c of `.claude/plans/angular_spaces_derived_from_symmetry.md`.
#
# ⚠ The keystone here is Q8.1, a ROUTE gate asserting OBJECT IDENTITY —
# deliberately NOT a value gate.  `[M]` 2026-08-31 (test-architect, 106 tool
# calls): inverting the routing predicate outright is bit-identical on
# end-to-end keff across slab/sphere/cylinder and reds 0 of 120 and 0 of 1913.
# A value gate is structurally blind to the routing decision this change IS
# (`plan-authoring` §10, third shape — the gate cannot detect its own
# campaign's success OR failure).  Identity, not equality: equality is what
# the REBUILT measure already satisfied, which is exactly why it hid the leak
# for as long as it did.

_ROUTED_RULES = [
    ("level_symmetric(4)", lambda: Quadrature.level_symmetric(4)),
    ("level_symmetric(8)", lambda: Quadrature.level_symmetric(8)),
    ("lebedev(11)", lambda: Quadrature.lebedev(11)),
    ("lebedev(17)", lambda: Quadrature.lebedev(17)),
    ("product(4,4)", lambda: Quadrature.product(n_mu=4, n_phi=4)),
    ("product(4,6)", lambda: Quadrature.product(n_mu=4, n_phi=6)),
    ("folded_product(2,4)", lambda: Quadrature.folded_product(n_mu=2, n_phi=4)),
    ("folded_product(4,8)", lambda: Quadrature.folded_product(n_mu=4, n_phi=8)),
]

_LIFTED_RULES = [
    ("gauss_legendre(2)", lambda: Quadrature.gauss_legendre(2)),
    ("gauss_legendre(8)", lambda: Quadrature.gauss_legendre(8)),
]


@pytest.mark.parametrize(
    "make", [r[1] for r in _ROUTED_RULES], ids=[r[0] for r in _ROUTED_RULES]
)
def test_q8_1_sphere_rules_hand_the_frame_their_OWN_measure(make) -> None:
    """⭐ KEYSTONE. A rule whose nodes already are three-component directions
    hands ``angular_frame`` its own measure — by IDENTITY, not by rebuilding
    an equal one.

    `[M]` 2026-09-01, before this landed: ``frame.measure is q.measure`` was
    ``False`` on 12 of 12 shipped rules; after, ``True`` on the 10 that route.
    """
    q = make()
    for L in (0, 1, 2):
        assert q.angular_frame(L).measure is q.measure


@pytest.mark.parametrize(
    "make", [r[1] for r in _ROUTED_RULES], ids=[r[0] for r in _ROUTED_RULES]
)
def test_q8_2_routing_carries_all_three_truths_the_rebuild_destroyed(make) -> None:
    """The rebuilt measure carried only nodes/weights/a literal support, so it
    silently dropped three things the rule knows. Identity restores all three
    at once — this pins each SEPARATELY, so a future partial copy cannot pass.

    `[M]` 2026-09-01 before the carve: support was falsified on 4 of 12 rules
    (2 slab + 2 fold), ``invariance_group`` on 10 of 12, ``exactness`` on 10 of 12.

    ⚠ **Not every clause bites on every row, and the count is stated rather
    than silently claimed** (``vv-principles`` #20). The two folded rules carry
    ``invariance_group = None`` and ``exactness = None`` on their OWN measure —
    deliberately, since a σ_y-quotient of an ``O_h``-invariant measure is not
    ``O_h``-invariant (pinned at ``test_measure.py:353`` and ``:949``) — so for
    those two rows those clauses read ``None is None`` and cannot fail. `[M]`
    the mutation arithmetic shows it exactly: stripping the group reds **6** of
    these 8 rows, forging the support reds the complementary **2**, and the
    full pre-carve rebuild reds all **8**.
    """
    q = make()
    m = q.angular_frame(2).measure
    assert m.support == q.measure.support
    assert m.invariance_group is q.measure.invariance_group
    assert m.exactness is q.measure.exactness


def test_q8_3_the_fold_keeps_its_quotient_tag_all_the_way_to_the_frame() -> None:
    """0.1c. ``folded_product`` declares the quotient :math:`S^2/\\sigma_y` —
    the measure layer has spoken quotients all along — and the rebuild
    overwrote it with :math:`S^2`, so the frame asserted a domain twice the
    size of the one its nodes cover.

    This is the one row of the carve that is NOT bit-identical in its tag:
    ``frame.measure_space`` moves ``L2[S^2]`` → ``L2[S^2/sigma_y]``.
    """
    q = Quadrature.folded_product(n_mu=4, n_phi=8)
    assert q.measure.support == "S^2/sigma_y"
    assert q.angular_frame(2).measure.support == "S^2/sigma_y"
    assert q.angular_frame(2).measure_space.name == "L2[S^2/sigma_y]"
    # …and an UNfolded sibling of the same family still says S^2, so the
    # assertion above is discriminating rather than a tautology on the tag.
    assert Quadrature.product(n_mu=4, n_phi=8).angular_frame(2).measure_space.name == "L2[S^2]"


@pytest.mark.parametrize(
    "make", [r[1] for r in _LIFTED_RULES], ids=[r[0] for r in _LIFTED_RULES]
)
def test_q8_4_the_1d_lift_is_still_a_FICTION_and_says_so(make) -> None:
    """⚠ SELF-RETIRING. A 1-D rule does NOT route, because the map it would
    need does not exist: a point of :math:`[-1,1]` is an ORBIT of the
    :math:`SO(2)` action, not a point of :math:`S^2`. The arrow that exists
    runs the other way (the quotient :math:`S^2 \\to [-1,1]`).

    So ``angular_frame`` still pads :math:`\\mu` to :math:`(\\mu, 0, 0)` and
    calls the result :math:`S^2` — ERR-080's construction, kept deliberately
    and named at
    :meth:`Quadrature._harmonic_frame_measure`.

    ⏏ **This test is the retirement trigger.** When Phase 3.4 gives the 1-D
    chart its trivial isotypic sub-basis, the branch disappears and this gate
    goes RED — which is what forces it to be rewritten rather than silently
    outlived.
    """
    q = make()
    frame_measure = q.angular_frame(2).measure
    assert frame_measure is not q.measure
    assert q.measure.support == "[-1,1]"
    assert frame_measure.support == "S^2"          # the fiction, stated
    assert q.measure.nodes.ndim == 1
    assert frame_measure.nodes.shape == (q.N, 3)
    # the fabricated azimuth: every node is padded onto the phi = 0 meridian
    np.testing.assert_array_equal(frame_measure.nodes[:, 1], np.zeros(q.N))
    np.testing.assert_array_equal(frame_measure.nodes[:, 2], np.zeros(q.N))


@pytest.mark.parametrize(
    "make", [r[1] for r in _ROUTED_RULES], ids=[r[0] for r in _ROUTED_RULES]
)
def test_q8_5_routing_moved_no_numbers(make) -> None:
    """A REGRESSION FLOOR, explicitly NOT the keystone (see the section header):
    it cannot see the routing decision, only that the decision moved no values.

    The nodes the rebuild produced were ``np.array_equal`` to the rule's own on
    every routed rule (`[M]` 10 of 12, exact shapes), so identity costs nothing
    numerically — that is what makes 0.1a landable ahead of 3.4.
    """
    q = make()
    m = q.angular_frame(2).measure
    rebuilt = np.column_stack([q.axis_cosines(0), q.axis_cosines(1), q.axis_cosines(2)])
    np.testing.assert_array_equal(m.nodes, rebuilt)
    np.testing.assert_array_equal(m.weights, q.weights)


@pytest.mark.parametrize(
    "make",
    [r[1] for r in _ROUTED_RULES + _LIFTED_RULES],
    ids=[r[0] for r in _ROUTED_RULES + _LIFTED_RULES],
)
def test_q8_6_frame_table_still_equals_spherical_harmonics_bit_for_bit(make) -> None:
    """``angular_frame(L).table`` == ``spherical_harmonics(L)``, bit for bit.

    ⚠ This gate exists BECAUSE of 0.1a, and it is the interesting kind of
    debt: the claim is old, the guarantee behind it is new and weaker.
    ``angular_frame``'s docstring has always asserted this equality, and until
    2026-09-01 it was true **by construction** — both spellings shared one
    literal ``column_stack(axis_cosines(0..2))`` expression, so no input could
    separate them. Routing the measure makes the two sides *independently
    assembled*: the frame reads ``measure.nodes``, ``spherical_harmonics``
    still column-stacks the cosines. They agree because those arrays are equal,
    which is a fact about the rules rather than a fact about the code.

    A claim that silently drops from by-construction to by-coincidence is
    exactly the demotion ``coding-standards`` warns about, and `[M]` nothing in
    the tree pinned this one — 0 tests compared the two. So it gets a gate.

    ⭐ NOT a tautology: both sides call the same ``SphericalHarmonicBasis``
    object, and that is fine — the independence lives in the **input
    assembly**, which is the criterion in
    ``feedback_verify_shared_primitive_pure_math``. Mutating either assembly
    path separates them; `[M]` 36 of 36 (12 rules × L ∈ {0,1,2}) agree today,
    against a positive control comparing mismatched degrees that reads False.
    """
    q = make()
    for L in (0, 1, 2):
        np.testing.assert_array_equal(q.angular_frame(L).table, q.spherical_harmonics(L))
