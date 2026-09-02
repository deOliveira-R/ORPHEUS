r"""Foundation tests for :mod:`orpheus.numerics.manifold`.

Every test here is tagged ``foundation``: the assertions verify the
intrinsic laws of the manifold type — dimension additivity, membership,
the quotient's dimension drop, and the recorded orbit-space derivation —
not a physics-equation claim with an L0..L3 ladder.  The mathematical
content is fixed by the standard invariant-theory construction (Procesi &
Schwarz, *Invent. Math.* **81** (1985) 539–554) and by the elementary
geometry of :math:`S^2`.

Three groups:

1. **The type's own laws** — the base is uninstantiable, variants are
   frozen and compare by value, and ``replace()`` re-runs the
   construction invariant rather than bypassing it.

2. **Membership, with both legs** — a contract predicate tested only
   against a broken instance validates the *raising*, not the *claim*
   (``vv-principles`` #11).  Every ``contains`` case therefore carries a
   positive leg (a genuine point set, must be admitted) beside its
   negative one.  The load-bearing negative is the **ERR-080 forgery**:
   a 1-D quadrature's ordinates :math:`(\mu, 0, 0)` declared to live on
   :math:`S^2`.

3. **The catalogued derivation, recomputed symbolically.**  ⭐ These are
   the tests ``D0.1`` requires of every catalogue entry, and they are
   deliberately written *before* the derivation engine exists: an engine
   ships when it reproduces them by computation instead of by lookup.  A
   specification written first cannot be shaped to flatter the
   implementation (``vv-principles`` #17).


4. **The arrows** (#429 tracker 2.3, 2026-09-02).  A :class:`ManifoldMap`
   is a frozen value with two endpoints; composition is functorial on the
   pushforward; and the two named maps — the Archimedes chart
   :math:`[-1,1] \times S^1 \to S^2` and the orbit barycentre
   :math:`S^2/SO(2)_a \to D^3` — each land where their codomain says,
   with the ERR-080 discriminator (the barycentre is OFF the sphere except
   at the poles) stated as a property of the map.

5. **The entry's own map and its reference** (#429 tracker 3.1,
   2026-09-02).  D0.1's last two derivation outputs are fields of the
   :class:`Quotient`: the quotient map :math:`\pi : M \to M/H` (an arrow
   INTO the entry, :math:`H`-invariant by the defining property, agreeing
   with the recorded symbolic invariants, and sitting between the two
   arrows of group 4 — :math:`\pi \circ \varphi_a = \mathrm{pr}_1`,
   :math:`\beta_a \circ \pi = ` the axial projection) and the pushforward
   reference (the hat-box measure on the axial entries; ``None`` where no
   shipped realization spells it), which the registry now READS.
"""

from __future__ import annotations

import dataclasses

import numpy as np
import pytest

from orpheus.numerics.manifold import (
    ambient_dim,
    CIRCLE,
    COSINE_INTERVAL,
    ENERGY,
    HALF_LINE,
    REAL_LINE,
    Ball,
    FundamentalDomain,
    SPHERE,
    UNIT_INTERVAL,
    Circle,
    EnergyGroups,
    IndexSet,
    Interval,
    Manifold,
    Product,
    Quotient,
    RealSpace,
    Sphere,
)
from orpheus.numerics.symmetry import SubgroupOfO3

pytestmark = pytest.mark.foundation


#: Every shipped variant, once.  A new member added to the module and not
#: to this tuple fails :func:`test_every_variant_is_reachable_from_this_list`.
_ALL: tuple[Manifold, ...] = (
    SPHERE,
    CIRCLE,
    COSINE_INTERVAL,
    UNIT_INTERVAL,
    HALF_LINE,
    REAL_LINE,
    RealSpace(1),
    RealSpace(3),
    ENERGY,
    IndexSet("angular", 8),
    SPHERE * COSINE_INTERVAL,
    Ball(2),
    FundamentalDomain(SPHERE, ((0.0, 1.0, 0.0),), "y>=0"),
)


# ---------------------------------------------------------------------------
# 1. The type's own laws
# ---------------------------------------------------------------------------


class TestTypeLaws:
    """The closed sum's construction invariants."""

    def test_the_base_is_not_a_manifold(self):
        """``Manifold()`` is uninstantiable — the base is a sum, not a member.

        Making the three TOTAL operations abstract is what forces every
        variant to answer them; a raising stub would let a half-built
        member ship and fail only when someone asked.
        """
        with pytest.raises(TypeError, match="abstract"):
            Manifold()  # type: ignore[abstract]

    @pytest.mark.parametrize("m", _ALL, ids=lambda m: m.name)
    def test_every_variant_answers_the_total_operations(self, m: Manifold):
        """dim, name and contains are TOTAL — no member may abstain."""
        assert isinstance(m.dim, int) and m.dim >= 0
        assert isinstance(m.name, str) and m.name

    def test_variants_are_frozen_and_compare_by_value(self):
        assert Sphere() == SPHERE
        assert Interval(-1.0, 1.0) == COSINE_INTERVAL
        assert Sphere() != Circle()
        assert len({Sphere(), SPHERE, CIRCLE}) == 2  # hashable, deduping
        with pytest.raises(dataclasses.FrozenInstanceError):
            COSINE_INTERVAL.a = 0.0  # type: ignore[misc]

    def test_replace_re_runs_the_construction_invariant(self):
        """``replace()`` must not be a hole in ``__post_init__``.

        The elegance rule (Pattern 4 ∩ Pattern 2): a same-type-producing
        path routes through construction so the law fires for free,
        rather than being restated — and drifting — at each site.
        """
        assert dataclasses.replace(COSINE_INTERVAL, b=2.0) == Interval(-1.0, 2.0)
        with pytest.raises(ValueError, match="requires a < b"):
            dataclasses.replace(COSINE_INTERVAL, a=2.0)

    def test_a_degenerate_interval_is_refused(self):
        """A point is not a 1-manifold."""
        with pytest.raises(ValueError, match="requires a < b"):
            Interval(1.0, 0.0)
        with pytest.raises(ValueError, match="requires a < b"):
            Interval(1.0, 1.0)

    def test_real_space_requires_a_positive_dimension(self):
        with pytest.raises(ValueError, match="requires d >= 1"):
            RealSpace(0)

    def test_the_names_reproduce_the_retired_string_tags(self):
        """The migration off ``Space = str`` must not re-word any message.

        Every one of these was a ``SPACE_*`` constant or an interpolated
        tag; a pinned error message or an ``L2[...]`` space name that
        quoted one keeps reading the same.
        """
        assert SPHERE.name == "S^2"
        assert CIRCLE.name == "S^1"
        assert COSINE_INTERVAL.name == "[-1,1]"
        assert UNIT_INTERVAL.name == "[0,1]"
        assert HALF_LINE.name == "[0,inf)"
        assert REAL_LINE.name == "R"
        assert ENERGY.name == "energy"
        assert RealSpace(1).name == "spatial_R1"
        assert IndexSet("angular").name == "index(angular)"


# ---------------------------------------------------------------------------
# 2. The product algebra
# ---------------------------------------------------------------------------


class TestProduct:
    r"""``M × N`` — the operation ``f"{a} × {b}"`` was performing."""

    def test_dimension_is_additive(self):
        r""":math:`\dim(M \times N) = \dim M + \dim N`."""
        for left in _ALL:
            for right in _ALL:
                assert (left * right).dim == left.dim + right.dim

    def test_the_name_reproduces_the_retired_interpolation(self):
        assert (COSINE_INTERVAL * CIRCLE).name == "[-1,1] × S^1"

    def test_the_product_is_a_manifold_so_it_composes(self):
        """Associativity of the concept, not of the constructor."""
        triple = (COSINE_INTERVAL * CIRCLE) * ENERGY
        other = COSINE_INTERVAL * (CIRCLE * ENERGY)
        assert isinstance(triple, Product) and isinstance(other, Product)
        assert triple.dim == other.dim == 2

    def test_multiplying_a_non_manifold_is_refused(self):
        with pytest.raises(TypeError):
            SPHERE * "S^1"  # type: ignore[operator]

    def test_membership_splits_the_coordinates(self):
        """A product admits a point iff each factor admits its block."""
        prod = COSINE_INTERVAL * CIRCLE  # ambient 1 + 2
        good = np.array([[0.5, 1.0, 0.0], [-0.5, 0.0, 1.0]])
        assert prod.contains(good)
        bad_left = np.array([[7.0, 1.0, 0.0]])  # 7 not in [-1,1]
        assert not prod.contains(bad_left)
        bad_right = np.array([[0.5, 2.0, 0.0]])  # not on S^1
        assert not prod.contains(bad_right)


# ---------------------------------------------------------------------------
# 3. Membership — positive AND negative legs (vv-principles #11)
# ---------------------------------------------------------------------------


class TestMembership:
    """``contains`` is what makes a support claim falsifiable."""

    def test_sphere_admits_unit_directions_and_refuses_the_err_080_forgery(self):
        r"""⭐ The load-bearing case: a 1-D rule's ordinates are not on :math:`S^2`.

        ``Quadrature.angular_frame`` builds an :math:`S^2` measure whose
        nodes are ``column_stack([mu_x, mu_y, mu_z])``.  On a 1-D rule
        ``mu_y = mu_z = 0`` means *"there is no azimuthal information"*,
        so the rows are :math:`(\mu, 0, 0)` with :math:`\|\Omega\| =
        |\mu| \ne 1` — points of :math:`[-1,1]`, not of the sphere.
        That forgery is the first link of ERR-080, and this predicate
        refuses it at construction, three hops before the wrong answer.

        Measured on ``gauss_legendre(8)``: the declared unit directions
        have norms in ``[0.1834, 0.9603]``, none within ``1e-12`` of 1.
        """
        mu = np.polynomial.legendre.leggauss(8)[0]
        forged = np.column_stack([mu, np.zeros_like(mu), np.zeros_like(mu)])

        # NEGATIVE leg — the forgery must be refused.
        assert not SPHERE.contains(forged)

        # POSITIVE leg — the same nodes, honestly normalised, are admitted.
        # Without this, the test would pin the raising and not the claim.
        honest = forged / np.linalg.norm(forged, axis=1, keepdims=True)
        assert SPHERE.contains(honest)

        # ...and the polar cosines DO live on the interval, which is the
        # manifold a 1-D rule's nodes actually belong to.
        assert COSINE_INTERVAL.contains(mu)

    def test_sphere_refuses_a_single_off_shell_node(self):
        """One bad row is enough — ``contains`` is a universal, not a mean."""
        good = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
        assert SPHERE.contains(good)
        assert not SPHERE.contains(np.array([[0.0, 0.0, 0.9]]))
        # and a node just OUTSIDE the tolerance, to pin where the edge is
        assert not SPHERE.contains(np.array([[1.0 + 1e-9, 0.0, 0.0]]))
        assert SPHERE.contains(np.array([[1.0 + 1e-14, 0.0, 0.0]]))

    def test_circle_admits_roots_of_unity(self):
        n = 8
        theta = 2.0 * np.pi * np.arange(n) / n
        assert CIRCLE.contains(np.column_stack([np.cos(theta), np.sin(theta)]))
        assert not CIRCLE.contains(np.array([[0.5, 0.5]]))

    def test_interval_includes_its_endpoints(self):
        """``[a, b]`` is closed — the endpoints are points of it."""
        assert COSINE_INTERVAL.contains(np.array([-1.0, 0.0, 1.0]))
        assert not COSINE_INTERVAL.contains(np.array([1.5]))
        assert not COSINE_INTERVAL.contains(np.array([-1.5]))

    def test_half_line_and_real_line(self):
        assert HALF_LINE.contains(np.array([0.0, 1.0, 1e12]))
        assert not HALF_LINE.contains(np.array([-1.0]))
        assert REAL_LINE.contains(np.array([-1e12, 0.0, 1e12]))
        assert not REAL_LINE.contains(np.array([np.inf]))

    def test_index_sets_admit_integers_only(self):
        assert IndexSet("angular", 8).contains(np.arange(8))
        assert not IndexSet("angular", 8).contains(np.array([8.0]))
        assert not IndexSet("angular", 8).contains(np.array([0.5]))
        assert not IndexSet("angular").contains(np.array([-1.0]))

    def test_energy_groups_admit_group_indices(self):
        assert EnergyGroups(2).contains(np.array([0.0, 1.0]))
        assert not EnergyGroups(2).contains(np.array([2.0]))

    def test_real_space_admits_any_finite_point(self):
        assert RealSpace(3).contains(np.zeros((4, 3)))
        assert not RealSpace(3).contains(np.full((1, 3), np.nan))

    def test_a_wrong_ambient_dimension_is_a_typed_refusal(self):
        """Not ``False`` — a shape mismatch is a caller error, not a verdict."""
        with pytest.raises(ValueError, match="ambient dimension 3"):
            SPHERE.contains(np.zeros((4, 2)))


# ---------------------------------------------------------------------------
# 4. The quotient, and its recorded derivation
# ---------------------------------------------------------------------------


class TestQuotient:
    r"""``S^2/SO(2) = [-1,1]``, and the fields that derive it."""

    @pytest.fixture(scope="class")
    def s2_mod_so2(self) -> Quotient:
        return SPHERE.quotient(SubgroupOfO3.SO2("x"))

    def test_the_realization_is_the_cosine_interval(self, s2_mod_so2):
        r""":math:`S^2/SO(2)_x \cong [-1,1]`, with :math:`\mu = \hat\Omega \cdot \hat x` — the slab\'s axis."""
        assert s2_mod_so2.realization == COSINE_INTERVAL
        assert s2_mod_so2.name == "S^2/SO2_x"

    def test_the_dimension_drops_by_the_group_dimension(self, s2_mod_so2):
        r""":math:`\dim S^2 = 2`, :math:`\dim SO(2) = 1`, quotient :math:`= 1`.

        The orbits are 1-dimensional (circles of constant latitude)
        everywhere except at the two poles — which is exactly the
        singular stratum below.
        """
        assert SPHERE.dim == 2
        assert s2_mod_so2.dim == 1

    def test_the_procesi_schwarz_matrix_matches_the_hand_derivation(
        self, s2_mod_so2
    ):
        r"""⭐ The symbolic regression test D0.1 requires of every entry.

        Invariants :math:`p_1 = z`, :math:`p_2 = x^2 + y^2` are
        algebraically independent (empty syzygy ideal), and

        .. math:: P_{ij} = \langle \nabla p_i, \nabla p_j \rangle
                  = \begin{pmatrix} 1 & 0 \\ 0 & 4p_2 \end{pmatrix},
                  \qquad \det P = 4 p_2 .

        This is the field an engine must compute rather than read, so
        the assertion is on the SYMBOLIC value, not on a string.
        """
        sp = pytest.importorskip("sympy")
        p2 = sp.Symbol("p2", real=True)

        assert s2_mod_so2.syzygy == ()
        assert sp.simplify(s2_mod_so2.det_gram - 4 * p2) == 0
        assert sp.simplify(s2_mod_so2.gram - sp.Matrix([[1, 0], [0, 4 * p2]])
                           ) == sp.zeros(2, 2)

    def test_det_P_vanishes_exactly_on_the_singular_stratum(self, s2_mod_so2):
        r"""⭐ The stratum is DERIVED, not declared — and this proves it.

        On :math:`S^2` the sphere's own ideal gives :math:`p_2 = 1 -
        \mu^2`, so :math:`\det P = 4(1-\mu^2)`, whose zeros are
        :math:`\mu = \pm 1`: the two points with full stabilizer, where
        the orbit degenerates from a circle to a point.  That is why the
        orbit space is an orbifold rather than a quotient manifold, and
        it is the same polynomial as the curvilinear
        angular-redistribution coefficient in
        :math:`(1/r)\,\partial_\mu[(1-\mu^2)\psi]`.
        """
        sp = pytest.importorskip("sympy")
        mu, p2 = sp.symbols("mu p2", real=True)

        det_on_sphere = sp.simplify(s2_mod_so2.det_gram.subs(p2, 1 - mu**2))
        assert sp.simplify(det_on_sphere - 4 * (1 - mu**2)) == 0

        roots = sorted(sp.solve(sp.Eq(det_on_sphere, 0), mu))
        assert [float(r) for r in roots] == [-1.0, 1.0]
        assert not s2_mod_so2.is_free

        # The stratum is a LOCUS in the realization's coordinate, and it
        # is the SAME locus det P cuts out — 1 - u0^2, vanishing at the
        # poles.  ⛔ It was a `tuple[float, ...]` until 2026-08-31, which
        # could not hold the sigma_y entry's stratum (a circle); the
        # first entry's SHAPE had become the field's type.
        u0 = sp.Symbol("u0", real=True)
        assert sp.simplify(s2_mod_so2.singular_stratum - (1 - u0**2)) == 0
        assert sorted(
            float(r) for r in sp.solve(sp.Eq(s2_mod_so2.singular_stratum, 0), u0)
        ) == [-1.0, 1.0]

    def test_the_stratum_lies_on_the_realization(self, s2_mod_so2):
        """A singular point is still a point of the orbit space.

        Solved from the locus rather than read from a literal, so the
        assertion survived the field's retyping from a float tuple
        (2026-08-31) without weakening: it still says the stratum's
        points are points of the realization.
        """
        sp = pytest.importorskip("sympy")
        u0 = sp.Symbol("u0", real=True)
        pts = [float(r) for r in sp.solve(sp.Eq(s2_mod_so2.singular_stratum, 0), u0)]
        assert s2_mod_so2.contains(np.asarray(sorted(pts)))

    def test_the_entry_records_who_derived_it(self, s2_mod_so2):
        """A MIXED hand/engine state must be expressible, or an
        incremental engine rollout would have to be all-or-nothing."""
        assert s2_mod_so2.derived_by == "hand"

    @pytest.mark.parametrize("base", [SPHERE, CIRCLE, COSINE_INTERVAL],
                             ids=lambda m: m.name)
    def test_the_trivial_quotient_is_derived_not_tabulated(self, base):
        r"""``M/{e} = M`` for EVERY manifold — a theorem, not a table row.

        Run the same procedure on the trivial group: its invariant ring
        is the whole polynomial ring, minimally generated by the
        coordinates, so :math:`P = I` and :math:`\det P = 1` — which
        vanishes nowhere, giving an empty singular stratum, i.e. a free
        action.  That is right vacuously (the only element is the
        identity), so this case doubles as a **positive control on the
        machinery**: the procedure reproduces a known answer.
        """
        sp = pytest.importorskip("sympy")
        q = base.quotient(SubgroupOfO3.Trivial)

        assert q.realization == base
        assert sp.simplify(q.det_gram - 1) == 0
        assert q.singular_stratum is None
        assert q.is_free
        assert q.dim == base.dim

    def test_the_derived_orbit_space_agrees_with_the_hand_written_table(self):
        r"""⭐ Two independent producers of one orbit space must agree.

        ``AngularSymmetry.support`` (``quadrature/registry.py``) was a
        **hand-written table** — ``SO2 -> COSINE_INTERVAL``,
        ``Trivial -> SPHERE``, raise otherwise.  :meth:`Manifold.quotient`
        **derives** the same orbit space from invariant theory (the
        Procesi--Schwarz matrix of the invariant ring).  Neither read the
        other, so this was a genuine two-implementation agreement and the
        only check that either was right about :math:`S^2/G^0`.

        ⛔ **RE-SCOPED 2026-09-01 by tracker 2.4, which single-sourced the
        SO(2) row.** The registry now DERIVES its domain through
        ``SPHERE.quotient(spent)`` — so for the axial row ``theirs`` is the
        very object ``mine`` was built from, and asking whether they agree
        is asking one call whether it equals itself (``coding-standards``:
        single-sourcing a duplicate demotes every gate that compared its
        copies, and the fix stays; the gate's DESCRIPTION is what moves).
        What survives, and is asserted below: (i) the hand-written
        ``expected`` column is still an INDEPENDENT pin on the derivation's
        realization — the one input that could make it fail is a wrong
        Procesi–Schwarz derivation; (ii) the ``Trivial`` row is still two
        producers, because the registry answers ``SPHERE`` explicitly while
        the catalogue re-derives ``M/{e}``; (iii) the axial row's registry
        answer must be the ORBIT SPACE, carrying the spent group, not its
        chart — the reading a spatial rule could share.

        ⭐ **Strengthened 2026-09-01 by the very retype it was written to
        pin** (tracker 2.0c).  It used to compare ``mine.name == theirs``
        because ``theirs`` was a *string* — the strongest claim the string
        vocabulary admitted.  Both sides are now ``Manifold``, so the gate
        asserts **object equality**: not merely that two producers spell the
        orbit space the same way, but that they produce the same point set.
        A name gate is satisfied by any self-consistent lie; this one is not.
        (``coding-standards``' mirror clause: a retirement can silently
        PROMOTE a gate's claim class, and the docstring must move with it or
        it goes on advertising the weaker claim.)

        ⚠ If this ever reddens, do NOT re-baseline: one of the two is wrong
        about a quotient, and which one is a mathematical question, not a
        test-maintenance one.
        """
        from orpheus.numerics.quadrature.registry import AngularSymmetry

        rows = {
            SubgroupOfO3.SO2("x"): COSINE_INTERVAL,
            SubgroupOfO3.SO2("z"): COSINE_INTERVAL,
            SubgroupOfO3.Trivial: SPHERE,
        }
        for group, expected in rows.items():
            derived = SPHERE.quotient(group)
            mine = derived.realization
            theirs = AngularSymmetry(
                continuous_isotropy=group,
                discrete_residual=SubgroupOfO3.Trivial,
            ).support
            # (i) the independent pin on the derivation
            assert mine == expected, f"{group.name}: mine gave {mine.name}"
            if group == SubgroupOfO3.Trivial:
                # (ii) still two producers: an explicit SPHERE vs M/{e}
                assert theirs == SPHERE and mine == theirs
            else:
                # (iii) the registry hands out the orbit space itself
                assert isinstance(theirs, Quotient) and theirs == derived, (
                    f"AngularSymmetry.support for S^2/{group.name} is "
                    f"{theirs!r}, not the derived orbit space"
                )
                assert theirs.by == group and theirs.realization == expected

    def test_an_uncatalogued_quotient_names_the_missing_WORK(self):
        """The engine's absence is a pick-up-able work item, not a wall.

        The message must name the manifold, the group, the derivation
        procedure, and where to register the result — so a fresh session
        can act on it without reading this plan.
        """
        with pytest.raises(NotImplementedError) as exc:
            SPHERE.quotient(SubgroupOfO3.OctahedralOh)
        msg = str(exc.value)
        assert "S^2/Oh" in msg
        assert "Procesi-Schwarz" in msg
        assert "_ORBIT_CATALOGUE" in msg
        # ...and says what IS catalogued, by class name (labelled
        # as such, since the request is named by MANIFOLD name).
        assert "Sphere/SO2_x" in msg and "Sphere/SO2_z" in msg
        assert "manifold CLASS" in msg



# ---------------------------------------------------------------------------
# 4b. The sigma_y fold — the entry that forced the two-slot ruling
# ---------------------------------------------------------------------------


class TestTheSigmaYFoldIsExpressibleAndDiscriminating:
    r"""``S^2/sigma_y``: the SHIPPED cylindrical fold, as a typed quotient.

    ⭐ This entry, not the first, is what forced ``Quotient`` to carry two
    coordinate systems.  ``sigma_y`` is FINITE, so ``dim`` does not drop
    (:math:`2 - 0 = 2`): the invariant chart buys no reduction and the
    section is canonical, which is the exact opposite of ``SO(2)``'s
    situation and is why one worked entry could not expose the fork.
    """

    @pytest.fixture(scope="class")
    def fold(self):
        return SPHERE.quotient(SubgroupOfO3.Mirror("y"))

    def test_the_derivation_reproduces_procesi_schwarz(self, fold):
        """``P = diag(1, 1, 4 p_3)``, ``det P = 4 p_3``, empty syzygy.

        The syzygy is empty by THEOREM, not by inspection:
        Chevalley--Shephard--Todd makes a reflection group's invariant
        ring polynomial, so the three generators are algebraically
        independent.
        """
        sp = pytest.importorskip("sympy")
        u2 = sp.Symbol("u2", real=True)

        assert fold.syzygy == ()
        assert sp.simplify(fold.det_gram - 4 * u2) == 0
        assert sp.simplify(
            fold.gram - sp.diag(1, 1, 4 * u2)
        ) == sp.zeros(3, 3)
        assert fold.realization == Ball(2)
        assert fold.dim == 2  # H is FINITE: no dimension drop
        assert not fold.is_free

    def test_the_stratum_is_a_CIRCLE_not_a_point_set(self, fold):
        """⛔ The field that could not hold this is why it was retyped.

        ``det P = 4 x_a^2`` vanishes on the mirror's fixed great circle,
        which in the realization's coordinates is the disk's BOUNDARY.
        A ``tuple[float, ...]`` cannot say that; the first catalogued
        entry's shape had become the field's type.
        """
        sp = pytest.importorskip("sympy")
        u0, u1 = sp.symbols("u0 u1", real=True)
        assert sp.simplify(fold.singular_stratum - (1 - u0**2 - u1**2)) == 0
        # It really is a curve: a one-parameter family of solutions.
        sols = sp.solve(sp.Eq(fold.singular_stratum, 0), u1)
        assert len(sols) == 2 and all(u0 in s.free_symbols for s in sols)

    def test_the_section_admits_the_SHIPPED_nodes_and_refuses_both_forgeries(
        self, fold
    ):
        r"""⭐ The load-bearing gate: three legs, on production data.

        ``DiscreteMeasure.quotient`` keeps ``nodes[representative]`` — a
        SELECTION applying no chart — so what the tree emits carries the
        base's three ambient columns.  Only the fundamental domain can
        validate that, and it must refuse BOTH wrong inputs:

        * the ORBIT TWINS (``mu_y -> -mu_y``): the wrong representative,
          which ``realization=SPHERE`` would have admitted, making
          ``Quotient.contains`` bit-for-bit ``SPHERE.contains`` — a
          predicate that cannot distinguish ``M/H`` from ``M``;
        * the ERR-080 forgery, which is not unit-norm.
        """
        from orpheus.numerics.quadrature.directional import Quadrature

        shipped = Quadrature.folded_product(n_mu=4, n_phi=8).measure.nodes
        assert shipped.shape == (16, 3)

        assert fold.contains(shipped)                                # POSITIVE
        twins = shipped * np.array([1.0, -1.0, 1.0])
        assert not fold.contains(twins)                              # NEGATIVE
        forgery = np.column_stack(
            [np.polynomial.legendre.leggauss(8)[0], np.zeros(8), np.zeros(8)]
        )
        assert not fold.contains(forgery)                            # NEGATIVE

    def test_the_chart_is_MODE_12_BLIND_to_the_forgery_and_the_section_is_not(
        self, fold
    ):
        r"""⛔ Why the chart alone could not have been the answer.

        The chart ``(x, y, z) -> (x, z)`` drops exactly the coordinate
        ERR-080 corrupts, so the forged row ``(mu, 0)`` is a PERFECTLY
        LEGAL point of the disk — it is the orbit of the real direction
        pair ``(mu, +-sqrt(1-mu^2), 0)``.  The mechanism is exact, not
        statistical: ``vv-principles`` Mode 12, the measured functional's
        stabiliser containing the error class.
        """
        forgery = np.column_stack(
            [np.polynomial.legendre.leggauss(8)[0], np.zeros(8), np.zeros(8)]
        )
        charted = forgery[:, [0, 2]]
        assert fold.realization.contains(charted)          # blind, by design
        assert not fold.fundamental_domain.contains(forgery)   # not blind

    def test_the_half_space_is_CLOSED_because_production_marches_from_it(
        self, fold
    ):
        r"""⛔ Strict ``> 0`` would refuse a direction the sweep uses.

        ``[M]`` the cylindrical march seeds each level at :math:`\xi = 0`
        exactly — ON the stratum — so the fundamental domain's
        inequality must be non-strict.  This is
        ``coding-elegance`` anti-pattern #18's half (ii): every LEGAL
        value must be admitted, which is a claim about the producers,
        and it is the half that gets skipped.

        ⚠ The shipped folded rule cannot witness this: its nodes are
        strictly positive (the even-``n_phi`` staggering makes the fold
        free on them), so the seed is the only witness available.
        """
        # BUILT, not transcribed.  [M] hand-typed to 12 dp these sit
        # 4.79e-13 off S^2 — 48 % of _MEMBERSHIP_ATOL's budget for no
        # reason, where production's are exact (0.0).  The seed is the
        # level's start point at xi = 0, so (-sqrt(1-mu^2), 0, mu) on the
        # rule's own level cosines reproduces it to machine precision
        # without coupling a numerics test to `orpheus.sn`.
        mu = np.polynomial.legendre.leggauss(4)[0]
        seeds = np.column_stack([-np.sqrt(1.0 - mu**2), np.zeros(4), mu])
        np.testing.assert_allclose(np.linalg.norm(seeds, axis=1), 1.0, atol=0.0)
        assert fold.contains(seeds)          # closed: ON the boundary, admitted

    def test_all_three_mirrors_share_ONE_derivation(self):
        """Three catalogue keys, one procedure — it reads the axis off
        the group.  A per-axis copy would be three places for the same
        math to drift."""
        from orpheus.numerics.manifold import _ORBIT_CATALOGUE

        builders = {
            _ORBIT_CATALOGUE[(Sphere, f"sigma_{a}")] for a in "xyz"
        }
        assert len(builders) == 1
        for a, ax in zip("xyz", range(3)):
            q = SPHERE.quotient(SubgroupOfO3.Mirror(a))
            fd = q.fundamental_domain
            assert isinstance(fd, FundamentalDomain)   # narrows, and IS a claim
            normal = fd.normals[0]
            assert normal[ax] == 1.0 and sum(normal) == 1.0

    def test_a_disk_is_not_its_bounding_square(self):
        """``Product(COSINE_INTERVAL, COSINE_INTERVAL)`` was the shipped
        2-D member nearest the answer, and it is the SQUARE — it admits
        ``(eta, mu)`` pairs corresponding to no direction at all."""
        square = COSINE_INTERVAL * COSINE_INTERVAL
        witness = np.array([[0.9, 0.9]])          # eta^2+mu^2 = 1.62 > 1
        assert square.contains(witness)
        assert not Ball(2).contains(witness)


class TestTheTwoCoordinateSystemsMustAgree:
    """The construction invariant, and the case it catches."""

    def test_a_mismatched_fundamental_domain_is_REFUSED_at_construction(self):
        r"""§6c: the gate lands with its own witness.

        A hemisphere (``dim`` 2, no antipodal pair) offered against a
        1-D realization is exactly the mistake the rule exists to catch —
        an equality normal forgotten, so the section describes a bigger
        object than the chart.
        """
        from orpheus.numerics.manifold import Quotient

        hemisphere = FundamentalDomain(SPHERE, ((0.0, 1.0, 0.0),), "y>=0")
        assert hemisphere.dim == 2 and COSINE_INTERVAL.dim == 1
        with pytest.raises(ValueError, match="must describe the same"):
            Quotient(
                base=SPHERE,
                by=SubgroupOfO3.Mirror("y"),
                realization=COSINE_INTERVAL,
                orbit_coordinates=lambda points: points[:, 1],
                fundamental_domain=hemisphere,
            )

    def test_the_shipped_entries_SATISFY_it(self):
        """The positive leg — without it the raise above is only evidence
        that the method raises when told to (``vv-principles`` #11)."""
        for q in (
            SPHERE.quotient(SubgroupOfO3.Mirror("y")),
            SPHERE.quotient(SubgroupOfO3.Trivial),
            SPHERE.quotient(SubgroupOfO3.SO2("x")),
        ):
            if q.fundamental_domain is not None:
                assert q.fundamental_domain.dim == q.realization.dim

    def test_an_antipodal_pair_spells_an_equality_and_drops_a_dimension(self):
        """One tuple expresses both a half-space and a hyperplane, which
        is what lets the sigma_y hemisphere and an SO(2) half-meridian
        share a type."""
        e_y, m_y, e_x = (0.0, 1.0, 0.0), (0.0, -1.0, 0.0), (1.0, 0.0, 0.0)
        assert FundamentalDomain(SPHERE, (e_y,), "half").dim == 2
        assert FundamentalDomain(SPHERE, (e_y, m_y, e_x), "meridian").dim == 1


# ---------------------------------------------------------------------------
# 5. Exhaustiveness — a new member must not be able to hide
# ---------------------------------------------------------------------------


def test_every_variant_is_reachable_from_this_modules_list():
    """A member added to the module but not exercised here fails HERE.

    The closed sum's whole benefit is that an operation can be checked
    against every member; that benefit is only real if the member list
    is itself pinned.  ``Manifold.__subclasses__`` is the tree's own
    answer, so this cannot drift from the module.
    """
    from orpheus.numerics import manifold as mod

    declared = {
        c.__name__
        for c in Manifold.__subclasses__()
        if c.__module__ == mod.__name__
    }
    exercised = {type(m).__name__ for m in _ALL} | {"Quotient"}
    assert declared == exercised, (
        f"variants declared but not exercised: {sorted(declared - exercised)}; "
        f"exercised but not declared: {sorted(exercised - declared)}"
    )


def test_ambient_dimension_is_defined_for_every_variant():
    """``ambient_dim`` decides how many columns ``contains`` consumes.

    Its ``match`` is deliberately exhaustive with a raising fallthrough,
    so a new member that forgets it fails loudly rather than silently
    mis-splitting a product's coordinates.

    Public since #429 tracker 2.1 (it was ``_ambient``): an
    :class:`~orpheus.numerics.basis.indicator_basis.IndicatorBasis` asks it
    from outside this module, to check that its per-axis partition has one
    axis per coordinate of the manifold it partitions.
    """
    for m in _ALL:
        assert ambient_dim(m) >= 1
    assert ambient_dim(SPHERE.quotient(SubgroupOfO3.SO2("x"))) == 1
    assert ambient_dim(SPHERE * CIRCLE) == 5


# ---------------------------------------------------------------------------
# 4. Maps between manifolds — the arrows (#429 tracker 2.3, 2026-09-02)
# ---------------------------------------------------------------------------


class TestManifoldMap:
    r"""A :class:`ManifoldMap` is a typed arrow, and the tree's three maps
    around the quotient are instances of it.

    The intrinsic laws (``feedback_test_intrinsic_properties``): a map is a
    frozen value with two endpoints; composition is refused across
    mismatched endpoints and is functorial on the pushforward,
    :math:`(\psi \circ \varphi)_*\mu = \psi_*(\varphi_*\mu)`.  Then the two
    named maps, each with a positive leg AND the refusal leg
    (``vv-principles`` #11): :func:`archimedes` lands on the sphere for
    every axis and collapses the fibre onto the stratum; :func:`barycentre`
    lands INSIDE the ball and on the sphere only at the poles — which is the
    ERR-080 discriminator stated as a property of a map rather than of a
    quadrature — and refuses any orbit space that is not an axial one.
    """

    def test_a_map_is_a_frozen_value_with_two_endpoints(self) -> None:
        from orpheus.numerics.manifold import ManifoldMap

        double = ManifoldMap(COSINE_INTERVAL, Interval(-2.0, 2.0), lambda x: 2.0 * x)
        assert double.domain == COSINE_INTERVAL
        assert double.codomain == Interval(-2.0, 2.0)
        assert np.array_equal(double(np.array([0.5, -1.0])), [1.0, -2.0])
        with pytest.raises(dataclasses.FrozenInstanceError):
            double.codomain = SPHERE  # type: ignore[misc]

    def test_composition_requires_matching_endpoints(self) -> None:
        from orpheus.numerics.manifold import ManifoldMap

        phi = ManifoldMap(COSINE_INTERVAL, UNIT_INTERVAL, lambda x: 0.5 * (x + 1.0))
        psi = ManifoldMap(UNIT_INTERVAL, UNIT_INTERVAL, lambda x: x * x)
        both = psi @ phi
        assert both.domain is COSINE_INTERVAL and both.codomain is UNIT_INTERVAL
        assert np.array_equal(both(np.array([1.0, -1.0, 0.0])), [1.0, 0.0, 0.25])

        # phi lands on [0,1]; a map out of [-1,1] cannot follow it.
        with pytest.raises(ValueError, match="cannot compose"):
            _ = phi @ phi
        with pytest.raises(TypeError):
            _ = phi @ "not a map"  # type: ignore[operator]

    def test_functoriality_the_pushforward_of_a_composite_is_the_composite_of_pushforwards(
        self,
    ) -> None:
        from orpheus.numerics.manifold import ManifoldMap
        from orpheus.numerics.measure import gauss_legendre

        phi = ManifoldMap(COSINE_INTERVAL, UNIT_INTERVAL, lambda x: 0.5 * (x + 1.0))
        psi = ManifoldMap(UNIT_INTERVAL, UNIT_INTERVAL, lambda x: x * x)
        mu = gauss_legendre(5)

        at_once = mu.pushforward(psi @ phi)
        stepwise = mu.pushforward(phi).pushforward(psi)
        assert np.array_equal(at_once.nodes, stepwise.nodes)
        assert np.array_equal(at_once.weights, stepwise.weights)
        assert at_once.support == stepwise.support == UNIT_INTERVAL

    @pytest.mark.parametrize("axis", ["x", "y", "z"])
    def test_archimedes_lands_on_the_sphere_and_collapses_the_fibre_at_the_poles(
        self, axis: str
    ) -> None:
        from orpheus.numerics.manifold import AXIS_INDEX, archimedes

        chart = archimedes(axis)
        assert chart.domain == COSINE_INTERVAL * CIRCLE
        assert chart.codomain is SPHERE

        mu = np.linspace(-1.0, 1.0, 7)
        phi = 2.0 * np.pi * np.arange(8) / 8
        # Cartesian order of the product measure: outer mu, inner phi.
        points = np.column_stack(
            [np.repeat(mu, 8), np.tile(np.cos(phi), 7), np.tile(np.sin(phi), 7)]
        )
        image = chart(points)
        assert SPHERE.contains(image)
        a = AXIS_INDEX[axis]
        assert np.array_equal(image[:, a], np.repeat(mu, 8))
        # Both poles are the image of a WHOLE fibre — the singular stratum.
        pole = np.zeros(3)
        pole[a] = 1.0
        assert np.array_equal(image[-8:], np.tile(pole, (8, 1)))
        assert np.array_equal(image[:8], np.tile(-pole, (8, 1)))

    def test_archimedes_about_z_is_the_labelled_equation_verbatim(self) -> None:
        r""":eq:`product-mu-phi-cosines`: :math:`\mu_x = \sin\theta\cos\varphi`,
        :math:`\mu_y = \sin\theta\sin\varphi`, :math:`\mu_z = \mu` — spelled
        by hand, and bit-identical to the chart."""
        from orpheus.numerics.manifold import archimedes

        mu = np.array([-0.9, -0.3, 0.0, 0.4, 0.75])
        phi = np.array([0.3, 1.1, 2.9, 4.0, 5.5])
        points = np.column_stack([mu, np.cos(phi), np.sin(phi)])
        sin_theta = np.sqrt(1.0 - mu**2)
        by_hand = np.column_stack([sin_theta * np.cos(phi), sin_theta * np.sin(phi), mu])
        assert np.array_equal(archimedes("z")(points), by_hand)
        assert archimedes("z") is archimedes("z")  # one object per axis
        with pytest.raises(ValueError, match="axis must be one of"):
            archimedes("w")

    @pytest.mark.parametrize("axis", ["x", "y", "z"])
    def test_barycentre_lands_in_the_ball_and_on_the_sphere_only_at_the_poles(
        self, axis: str
    ) -> None:
        from orpheus.numerics.manifold import AXIS_INDEX, barycentre

        orbit_space = SPHERE.quotient(SubgroupOfO3.SO2(axis))
        centre = barycentre(orbit_space)
        assert centre.domain is orbit_space  # the memoised entry itself
        assert centre.codomain == Ball(3)

        mu = np.linspace(-1.0, 1.0, 9)
        image = centre(mu)
        assert np.array_equal(image, centre(mu[:, None]))  # (N,) and (N,1) agree
        assert Ball(3).contains(image)
        # Off the sphere everywhere except the poles: this is ERR-080's
        # forgery, stated as what the map's image IS.
        assert not SPHERE.contains(image[1:-1])
        assert SPHERE.contains(image[[0, -1]])
        # 1 - |b|^2 is the squared orbit radius, det P / 4 = 1 - mu^2.
        assert np.allclose(1.0 - np.sum(image * image, axis=1), 1.0 - mu**2)
        a = AXIS_INDEX[axis]
        assert np.array_equal(image[:, a], mu)
        assert not np.any(image[:, [i for i in range(3) if i != a]])

    def test_barycentre_refuses_anything_but_an_axial_orbit_space(self) -> None:
        from orpheus.numerics.manifold import barycentre

        barycentre(SPHERE.quotient(SubgroupOfO3.SO2("x")))  # positive control
        for not_axial in (
            SPHERE.quotient(SubgroupOfO3.Mirror("y")),
            SPHERE.quotient(SubgroupOfO3.Trivial),
            COSINE_INTERVAL,
        ):
            with pytest.raises(ValueError, match="axial rotation group"):
                barycentre(not_axial)  # type: ignore[arg-type]

    @pytest.mark.parametrize("axis", ["x", "y", "z"])
    def test_barycentre_is_the_embedding_the_invariance_check_reads(
        self, axis: str
    ) -> None:
        """Pattern 2: the honest spelling of the map READS the map.

        ``symmetry._embedded_nodes`` is what ``is_invariant`` and the
        motion-preservation check embed a polar marginal through; since
        2.3 it is this map, so the two cannot drift.
        """
        from orpheus.numerics.manifold import barycentre
        from orpheus.numerics.quadrature.rules_1d import gauss_legendre_on_polar_orbit
        from orpheus.numerics.symmetry import _embedded_nodes

        rule = gauss_legendre_on_polar_orbit(8, axis)
        assert isinstance(rule.support, Quotient)
        assert np.array_equal(
            _embedded_nodes(rule), barycentre(rule.support)(rule.nodes)
        )


# ---------------------------------------------------------------------------
# 5. The entry's own map and its reference (#429 tracker 3.1)
# ---------------------------------------------------------------------------


def _random_directions(n: int, seed: int = 0) -> np.ndarray:
    rng = np.random.default_rng(seed)
    omega = rng.normal(size=(n, 3))
    return omega / np.linalg.norm(omega, axis=1, keepdims=True)


def _rotation_about(axis: int, theta: float) -> np.ndarray:
    """A proper rotation about coordinate axis ``axis`` — an element of
    :math:`SO(2)_a`, written by hand so the test does not read the group
    it is checking."""
    b, c = (axis + 1) % 3, (axis + 2) % 3
    R = np.eye(3)
    R[b, b] = R[c, c] = np.cos(theta)
    R[b, c], R[c, b] = -np.sin(theta), np.sin(theta)
    return R


#: Every shipped orbit space of the sphere: the six catalogue keys plus the
#: derived identity quotient.  A universal below is a universal over THIS.
_ROSTER = [
    SubgroupOfO3.SO2("x"), SubgroupOfO3.SO2("y"), SubgroupOfO3.SO2("z"),
    SubgroupOfO3.Mirror("x"), SubgroupOfO3.Mirror("y"), SubgroupOfO3.Mirror("z"),
    SubgroupOfO3.Trivial,
]


class TestTheEntryCarriesItsQuotientMapAndItsReference:
    r"""D0.1's last two derivation outputs are FIELDS of the entry.

    Until 2026-09-02 the quotient map :math:`\pi : S^2 \to S^2/SO(2)_a`
    *"was not a value anywhere"* and the pushforward reference lived on
    the registry as a tabulated twin.  Each test below pins one law the
    map must satisfy AS a quotient map — :math:`H`-invariance (its
    defining property), agreement with the symbolic invariants the entry
    records, the two witnesses against the arrows on either side of it,
    and the pushforward identity :eq:`discrete-measure-pushforward` on a
    real sphere rule — then the reference, and that the registry reads it
    off ONE object.
    """

    _AXES = ["x", "y", "z"]

    @pytest.mark.parametrize("axis", _AXES)
    def test_the_quotient_map_is_an_arrow_from_the_base_INTO_the_entry(
        self, axis: str
    ) -> None:
        """Codomain = the entry itself, never its chart codomain — the
        distinction tracker 2.4 made refusable, now carried by the map."""
        from orpheus.numerics.manifold import ManifoldMap

        q = SPHERE.quotient(SubgroupOfO3.SO2(axis))
        pi = q.quotient_map
        assert isinstance(pi, ManifoldMap)
        assert pi.domain is SPHERE
        assert pi.codomain is q
        assert pi.codomain != COSINE_INTERVAL
        mu = pi(_random_directions(64))
        assert mu.shape == (64,)
        assert q.contains(mu)

    @pytest.mark.parametrize("group", _ROSTER, ids=lambda g: g.name)
    def test_the_quotient_map_is_INVARIANT_under_the_group_it_quotients_by(
        self, group
    ) -> None:
        r"""The defining property of a quotient map: :math:`\pi(h\Omega) =
        \pi(\Omega)` for every :math:`h \in H` — bit-exactly here, because
        every shipped :math:`h` has exact 0/±1 entries in the rows
        :math:`\pi` reads.  With the negative leg (``vv-principles`` #11):
        an element OUTSIDE :math:`H` moves the image.

        ⚠ The outsider must be chosen: `[M]` 2026-09-02 (archivist)
        :math:`\pi_a` is ALSO invariant under :math:`\sigma_b`, :math:`b
        \neq a` — :math:`O(2)_a` and :math:`SO(2)_a` induce the same orbit
        partition, so a quotient map determines the PARTITION and the
        partition does not determine the group.  ``Quotient.by`` is a
        declaration, not something the map could recover.  The outsider
        below is a rotation about ANOTHER axis, which does move it."""
        q = SPHERE.quotient(group)
        pi = q.quotient_map
        omega = _random_directions(128, seed=1)
        before = pi(omega)
        if group.rotation_axis is not None:
            elements = [_rotation_about(group.rotation_axis, t) for t in (0.3, 1.7, 4.1)]
            outsider = _rotation_about((group.rotation_axis + 1) % 3, 0.3)
        elif group.mirror_axis is not None:
            elements = [np.diag([-1.0 if i == group.mirror_axis else 1.0 for i in range(3)])]
            outsider = np.diag([-1.0 if i != group.mirror_axis else 1.0 for i in range(3)])
        else:
            elements = [np.eye(3)]
            outsider = _rotation_about(0, 0.3)
        for h in elements:
            assert np.array_equal(pi(omega @ h.T), before), group.name
        assert not np.array_equal(pi(omega @ outsider.T), before), group.name

    @pytest.mark.parametrize("group", _ROSTER, ids=lambda g: g.name)
    def test_the_numeric_map_IS_the_recorded_symbolic_invariants(self, group) -> None:
        r"""⭐ The D0.1 tie: the engine would ``lambdify`` the surviving
        invariants; the hand entry spells them as a column selection.  The
        two must agree on random directions — bit-exactly, since every
        surviving invariant of the shipped entries is a coordinate
        function.  The convention pinned here: the realization's
        coordinates are the FIRST ``ambient_dim(realization)`` generators,
        each written as ``p_i - (invariant in x, y, z)`` (or the bare
        coordinate for the trivial entry)."""
        sp = pytest.importorskip("sympy")
        q = SPHERE.quotient(group)
        k = ambient_dim(q.realization)
        # Each generator is written `slot - invariant(ambient)`, the slot
        # being the ONE symbol it depends on with derivative exactly 1 — or
        # a bare ambient coordinate (the trivial entry).  Found structurally,
        # not by name: the sphere entries spell the ambient `x y z`, the
        # generic trivial derivation spells it `x0:n`.
        surviving = []
        for g in q.generators[:k]:
            slots = [sym for sym in g.free_symbols if sp.diff(g, sym) == 1]
            if isinstance(g, sp.Symbol):
                surviving.append(g)
            else:
                assert len(slots) == 1, (group.name, g)
                surviving.append(sp.solve(g, slots[0])[0])
        ambient = sorted(set().union(*(e.free_symbols for e in surviving)), key=str)
        assert len(ambient) <= 3
        # the derivations' coordinate order IS the name order (x < y < z,
        # x0 < x1 < x2), so the sorted symbols bind to the columns in order
        columns = [int(str(a)[1:]) if str(a)[1:].isdigit() else "xyz".index(str(a)) for a in ambient]
        symbolic = sp.lambdify(ambient, surviving, "numpy")
        omega = _random_directions(32, seed=2)
        expected = np.column_stack(np.broadcast_arrays(*symbolic(*omega[:, columns].T)))
        got = q.orbit_coordinates(omega)
        assert np.array_equal(got.reshape(32, -1), expected), group.name

    @pytest.mark.parametrize("axis", _AXES)
    def test_pi_after_archimedes_is_the_first_projection_BIT_EXACTLY(
        self, axis: str
    ) -> None:
        r""":math:`\pi_a \circ \varphi_a = \mathrm{pr}_1` — the arrow on the
        other side of the Archimedes parametrisation, composed through the
        typed ``@`` so the endpoints are checked too."""
        from orpheus.numerics.manifold import archimedes
        from orpheus.numerics.quadrature.rules_1d import gauss_legendre_on_mu

        q = SPHERE.quotient(SubgroupOfO3.SO2(axis))
        composite = q.quotient_map @ archimedes(axis)
        assert composite.domain == COSINE_INTERVAL * CIRCLE
        assert composite.codomain is q
        rng = np.random.default_rng(3)
        for n in (2, 4, 8, 16):
            mu = np.asarray(gauss_legendre_on_mu(n).nodes).reshape(-1)
            phi = rng.uniform(0.0, 2.0 * np.pi, size=mu.size)
            points = np.column_stack([mu, np.cos(phi), np.sin(phi)])
            assert np.array_equal(composite(points), mu), (axis, n)

    @pytest.mark.parametrize("axis", _AXES)
    def test_barycentre_after_pi_is_the_axial_projection(self, axis: str) -> None:
        r""":math:`\beta_a \circ \pi_a : \Omega \mapsto (\Omega \cdot \hat e_a)
        \hat e_a` — the orbit's centre, and the ERR-080 forgery's honest
        spelling: a point of the ball, on the sphere only at the poles."""
        from orpheus.numerics.manifold import AXIS_INDEX, barycentre

        q = SPHERE.quotient(SubgroupOfO3.SO2(axis))
        composite = barycentre(q) @ q.quotient_map
        assert composite.domain is SPHERE and composite.codomain == Ball(3)
        omega = _random_directions(256, seed=4)
        projection = np.zeros_like(omega)
        projection[:, AXIS_INDEX[axis]] = omega[:, AXIS_INDEX[axis]]
        assert np.array_equal(composite(omega), projection)
        assert Ball(3).contains(composite(omega))
        assert not SPHERE.contains(composite(omega))

    def test_a_sphere_rule_pushed_along_pi_lands_on_the_entry_and_obeys_the_pullback_identity(
        self,
    ) -> None:
        r"""The two induced arrows, on a REAL rule: :math:`\pi_*\mu` lives on
        :math:`S^2/SO(2)_x` because :math:`\pi` says so, and
        :math:`\int f \, d(\pi_*\mu) = \int (f \circ \pi) \, d\mu`
        (:eq:`discrete-measure-pushforward`) — checked on
        :math:`f(\mu) = \mu^2`, whose sphere integral is :math:`4\pi/3`.
        With the refusal leg: a rule on the chart cannot be pushed along a
        map out of the sphere."""
        from orpheus.numerics.quadrature.rules_1d import gauss_legendre_on_mu
        from orpheus.numerics.quadrature.rules_sphere import level_symmetric_sn

        q = SPHERE.quotient(SubgroupOfO3.SO2("x"))
        rule, _ = level_symmetric_sn(4)
        pushed = rule.pushforward(q.quotient_map)
        assert pushed.support is q
        assert np.array_equal(pushed.nodes, rule.nodes[:, 0])
        assert np.array_equal(pushed.weights, rule.weights)
        lhs = pushed.integrate(lambda mu: mu**2)
        rhs = rule.integrate(lambda omega: omega[:, 0] ** 2)
        assert lhs == rhs
        assert np.isclose(lhs, 4.0 * np.pi / 3.0)
        with pytest.raises(ValueError, match="map's domain must be the measure's support"):
            gauss_legendre_on_mu(4).pushforward(q.quotient_map)

    def test_the_axial_entries_carry_the_hat_box_reference_and_the_others_carry_None(
        self,
    ) -> None:
        r"""Archimedes: :math:`\mu_* d\Omega = 2\pi\,d\mu`, so the reference
        is Lebesgue on :math:`[-1,1]` — ``LEGENDRE``, identically, on all
        three axes.  ``None`` on the mirror (a weighted disk measure no
        shipped realization spells) and on the trivial quotient (Lebesgue
        on the BASE, whose orthogonal system the type does not carry)."""
        from orpheus.numerics.generating_measure import LEGENDRE

        for axis in self._AXES:
            assert SPHERE.quotient(SubgroupOfO3.SO2(axis)).reference is LEGENDRE
        assert LEGENDRE.support == COSINE_INTERVAL
        for axis in self._AXES:
            assert SPHERE.quotient(SubgroupOfO3.Mirror(axis)).reference is None
        assert SPHERE.quotient(SubgroupOfO3.Trivial).reference is None
        assert CIRCLE.quotient(SubgroupOfO3.Trivial).reference is None
        # ⚠ A reference lives on the REALIZATION — the chart codomain, where
        # the quotient map's image lands and where 2 pi d(mu) is written —
        # never on the entry, which would demand a measure carry an axis it
        # does not know.  `[M]` 2026-09-02 (archivist): nothing in orpheus/
        # reads `reference.support`, so this is the only pin on it.
        for group in _ROSTER:
            q = SPHERE.quotient(group)
            if q.reference is not None:
                assert q.reference.support == q.realization
                assert q.reference.support != q

    def test_the_map_is_DERIVED_so_it_takes_no_part_in_equality_and_survives_a_pickle(
        self,
    ) -> None:
        """Two honest constructions of one entry compare equal (a function
        has no value equality, and the map is a function of ``(base, by)``),
        and an entry — a support, hence a thing measures carry — pickles."""
        import pickle

        q = SPHERE.quotient(SubgroupOfO3.SO2("x"))
        twin = dataclasses.replace(q, orbit_coordinates=lambda points: points[:, 0])
        assert twin == q and hash(twin) == hash(q)
        assert twin.quotient_map.codomain is twin
        omega = _random_directions(8, seed=5)
        back = pickle.loads(pickle.dumps(q))
        assert back == q
        assert np.array_equal(back.quotient_map(omega), q.quotient_map(omega))
        assert back.reference == q.reference

    def test_no_quotient_is_minted_at_MODULE_scope_anywhere_in_the_package(self) -> None:
        r"""The function-scope import's safety condition, gated.

        ``_sphere_mod_so2`` imports ``LEGENDRE`` at FUNCTION scope; `[M]` the
        same line at module scope kills 7 of 7 fresh import orders (the cycle
        closes through ``measure``, which ``generating_measure`` imports, and
        ``numerics/__init__`` imports ``measure`` eagerly).  That is safe
        exactly while no quotient is built during package initialisation —
        a property of the CALL SITES, which a later module-scope
        ``SPHERE.quotient(...)`` anywhere in ``orpheus/`` would silently
        break.  So the condition is asserted here, by AST over the package,
        with a positive control proving the detector sees a violation.
        """
        import ast
        import pathlib

        import orpheus

        def module_scope_minting(tree: ast.AST) -> list[int]:
            hits: list[int] = []
            for node in getattr(tree, "body", []):
                if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
                    continue  # a body below this runs at CALL time, not import time
                for sub in ast.walk(node):
                    if (
                        isinstance(sub, ast.Call)
                        and isinstance(sub.func, ast.Attribute)
                        and sub.func.attr in ("quotient", "on_orbit_space")
                    ):
                        hits.append(sub.lineno)
            return hits

        # positive control: the detector flags a module-scope mint
        control = ast.parse("from m import SPHERE, G\nQ = SPHERE.quotient(G)\ndef f():\n    return SPHERE.quotient(G)\n")
        assert module_scope_minting(control) == [2]

        root = pathlib.Path(orpheus.__file__).parent
        offenders = {
            str(path.relative_to(root.parent)): lines
            for path in sorted(root.rglob("*.py"))
            if (lines := module_scope_minting(ast.parse(path.read_text())))
        }
        assert offenders == {}, (
            f"a quotient is minted at module scope in {offenders}; that runs "
            f"during package initialisation and re-opens the manifold -> "
            f"generating_measure -> measure import cycle"
        )
