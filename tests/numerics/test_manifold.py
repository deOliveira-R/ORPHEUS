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
    _all_coordinates,
    _coordinate_chart,
    quotient_onto,
)
from orpheus.numerics.symmetry import SubgroupOfO3

# 2.2b gates (test-architect, 2026-09-02) — the draft's imports, verbatim
import numpy as np
import pytest

from orpheus.geometry.transformation import RigidMotion
from orpheus.numerics.manifold import (
    COSINE_INTERVAL,
    SPHERE,
    Ball,
    Quotient,
    barycentre,
    quotient_onto,
)
from orpheus.numerics.measure import DiscreteMeasure
from orpheus.numerics.quadrature.directional import Quadrature
from orpheus.numerics.quadrature.registry import GEOMETRY_ANGULAR_SYMMETRY
from orpheus.numerics.quadrature.rules_1d import (
    gauss_legendre_on_mu,
    gauss_legendre_on_polar_orbit,
)
from orpheus.numerics.symmetry import SubgroupOfO3
from orpheus.numerics.invariance import permutation_under

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
        return SPHERE.quotient(SubgroupOfO3.O2("x"))

    def test_the_realization_is_the_cosine_interval(self, s2_mod_so2):
        r""":math:`S^2/O(2)_x \cong [-1,1]`, with :math:`\mu = \hat\Omega \cdot \hat x` — the slab\'s axis."""
        assert s2_mod_so2.realization == COSINE_INTERVAL
        assert s2_mod_so2.name == "S^2/O2_x"

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
            SubgroupOfO3.O2("x"): COSINE_INTERVAL,
            SubgroupOfO3.O2("z"): COSINE_INTERVAL,
            SubgroupOfO3.Trivial: SPHERE,
        }
        for group, expected in rows.items():
            derived = SPHERE.quotient(group)
            mine = derived.realization
            theirs = AngularSymmetry(
                spent=group, unspent=SubgroupOfO3.Trivial, owed=SubgroupOfO3.Trivial,
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
        assert "Sphere/O2_x" in msg and "Sphere/O2_z" in msg
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
        r"""§6c: the gate lands with its own witness — on the input that ONLY
        this clause refuses.

        A half-meridian (an antipodal normal PAIR spells an equality, so
        ``dim`` 1) offered against the disk (``dim`` 2) is the mistake this
        rule exists to catch: a section describing a smaller object than the
        chart.  ⚠ Until R4 of #434 (2026-09-03) this row offered a hemisphere
        against ``[-1,1]`` — an input that ALSO violates the dimension law
        ``dim(M/H) = dim M − dim(generic orbit)``, which now runs first, so it
        pinned whichever clause ran first rather than this one (`[M]`
        ``scratch/_r4_probe10.py``: the half-meridian input passes the
        dimension law, 2 − 0 = 2, and is refused by this clause verbatim).
        """
        from orpheus.numerics.manifold import Quotient, _coordinate_chart

        half_meridian = FundamentalDomain(
            SPHERE, ((0.0, 1.0, 0.0), (0.0, -1.0, 0.0)), "y=0"
        )
        assert half_meridian.dim == 1 and Ball(2).dim == 2
        chart, lift = _coordinate_chart([0, 2], 3)
        with pytest.raises(ValueError, match="must describe the same"):
            Quotient(
                base=SPHERE,
                by=SubgroupOfO3.Mirror("y"),
                realization=Ball(2),
                orbit_coordinates=chart,
                lift_coordinates=lift,
                lift_codomain=Ball(3),
                fundamental_domain=half_meridian,
            )

    def test_the_shipped_entries_SATISFY_it(self):
        """The positive leg — without it the raise above is only evidence
        that the method raises when told to (``vv-principles`` #11)."""
        for q in (
            SPHERE.quotient(SubgroupOfO3.Mirror("y")),
            SPHERE.quotient(SubgroupOfO3.Trivial),
            SPHERE.quotient(SubgroupOfO3.O2("x")),
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
    assert ambient_dim(SPHERE.quotient(SubgroupOfO3.O2("x"))) == 1
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

        orbit_space = SPHERE.quotient(SubgroupOfO3.O2(axis))
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

    def test_barycentre_refuses_a_bare_manifold_and_nothing_else(self) -> None:
        r"""Since R4 of #434 (2026-09-03) the barycentre map is the lift of
        EVERY catalogued entry — the Reynolds projector — so the refusal is
        for a point set that is not an orbit space at all; the widened
        positive side over all entries is
        ``TestR4TheLiftIsADerivationOutputNotATagBranch``."""
        from orpheus.numerics.manifold import barycentre

        for entry in (
            SPHERE.quotient(SubgroupOfO3.O2("x")),
            SPHERE.quotient(SubgroupOfO3.Mirror("y")),
            SPHERE.quotient(SubgroupOfO3.Trivial),
        ):
            assert barycentre(entry).domain == entry  # positive legs
        with pytest.raises(ValueError, match="not orbits"):
            barycentre(COSINE_INTERVAL)  # type: ignore[arg-type]

    @pytest.mark.parametrize("axis", ["x", "y", "z"])
    def test_barycentre_is_the_embedding_the_invariance_check_reads(
        self, axis: str
    ) -> None:
        """Pattern 2: the honest spelling of the map READS the map.

        ``invariance._embedded_nodes`` is what ``is_invariant_under`` and the
        motion-preservation check embed a polar marginal through; since
        2.3 it is this map, so the two cannot drift.
        """
        from orpheus.numerics.manifold import barycentre
        from orpheus.numerics.quadrature.rules_1d import gauss_legendre_on_polar_orbit
        from orpheus.numerics.invariance import _embedded_nodes

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
    SubgroupOfO3.O2("x"), SubgroupOfO3.O2("y"), SubgroupOfO3.O2("z"),
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

        q = SPHERE.quotient(SubgroupOfO3.O2(axis))
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

        q = SPHERE.quotient(SubgroupOfO3.O2(axis))
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

        q = SPHERE.quotient(SubgroupOfO3.O2(axis))
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
        :math:`S^2/O(2)_x` because :math:`\pi` says so, and
        :math:`\int f \, d(\pi_*\mu) = \int (f \circ \pi) \, d\mu`
        (:eq:`discrete-measure-pushforward`) — checked on
        :math:`f(\mu) = \mu^2`, whose sphere integral is :math:`4\pi/3`.
        With the refusal leg: a rule on the chart cannot be pushed along a
        map out of the sphere."""
        from orpheus.numerics.quadrature.rules_1d import gauss_legendre_on_mu
        from orpheus.numerics.quadrature.rules_sphere import level_symmetric_sn

        q = SPHERE.quotient(SubgroupOfO3.O2("x"))
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
            assert SPHERE.quotient(SubgroupOfO3.O2(axis)).reference is LEGENDRE
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

        q = SPHERE.quotient(SubgroupOfO3.O2("x"))
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

        ``_sphere_mod_o2`` imports ``LEGENDRE`` at FUNCTION scope; `[M]` the
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


# ══════════════════════════════════════════════════════════════════════
# The ENTRY's isotypic probe — which functions on my base descend to me
# (#429 tracker 3.4; the probe both ``Descent`` and the σ-even harmonic
#  sub-basis read, so the predicate is spelled ONCE)
# ══════════════════════════════════════════════════════════════════════


@pytest.mark.verifies("manifold-fibre-constancy")
class TestDescendingSlots:
    r"""``Quotient.descending_slots`` — the theorem about :math:`\pi`, gated.

    ⭐ Carries ``verifies("manifold-fibre-constancy")``: the probe IS
    :math:`f(g\,x) = f(x)` transcribed, evaluated at generic base points and
    their images under the group's generic elements. ``[M]`` 2026-09-02 that
    label was listed as an ORPHAN equation in the generated V&V matrix (zero
    tests carrying its marker), and it is not a ``documented``-sentineled one
    — so this is a real coverage edge, not a topic tag.

    A function on the base descends to :math:`M/H` **iff** it is constant on
    the fibres of the entry's own quotient map. The probe asks exactly that:
    it tabulates a basis at generic base points and at their images under the
    group's generic elements, and keeps the slots that agree.

    Claim layer: **structural (L0)**, closed form — for :math:`SO(2)_x` and
    the real harmonics the answer is a THEOREM (Schur: the trivial isotypic
    component is one-dimensional in every degree, spanned by
    :math:`Y_\ell^0`), so the positive leg is checkable against
    :math:`\{(\ell, \ell)\}` written independently of the probe.
    """

    def test_about_x_exactly_the_m_zero_slots_descend(self) -> None:
        r"""B5 — the descending real slots about :math:`x` are exactly :math:`\{(\ell, 0)\}`.

        ``[M]`` 2026-09-02 at :math:`L = 4`: **25 of 45** table slots test as
        invariant, of which **20 are** :math:`|m| > \ell` **padding** — so the
        honest count is **5 real of 25**, one per degree, and they sit on the
        rectangular layout's diagonal (column :math:`\ell + m` with
        :math:`m = 0`).

        The expected set is written from the representation theory, not read
        off the probe, so this row is a check and not a restatement.
        """
        from orpheus.numerics.basis.spherical_harmonic_basis import (
            SphericalHarmonicBasis,
        )

        for L in (1, 2, 3, 4, 5):
            basis = SphericalHarmonicBasis(L=L)
            mask = SPHERE.quotient(SubgroupOfO3.O2("x")).descending_slots(basis)
            real = mask & basis.live_slot_mask

            expected = np.zeros_like(real)
            for ell in range(L + 1):
                expected[ell, ell] = True  # column l + m at m = 0
            assert np.array_equal(real, expected), (
                f"L={L}: descending real slots {np.argwhere(real).tolist()}"
            )
            assert int(real.sum()) == L + 1

    def test_padding_slots_are_not_counted_as_descending(self) -> None:
        r"""B5b — the DENOMINATOR: a padded layout's :math:`|m| > \ell` slots descend vacuously.

        They are identically zero, so they are constant on every fibre of
        every group — counting them makes the mask read as "all of them" and
        hides the actual selection. ``[M]`` at :math:`L = 4`: 45 table slots,
        25 live, 25 pass the probe, **5** of those are live.
        """
        from orpheus.numerics.basis.spherical_harmonic_basis import (
            SphericalHarmonicBasis,
        )

        L = 4
        basis = SphericalHarmonicBasis(L=L)
        mask = SPHERE.quotient(SubgroupOfO3.O2("x")).descending_slots(basis)

        assert mask.shape == (L + 1, 2 * L + 1)
        assert mask.size == 45
        assert int(basis.live_slot_mask.sum()) == 25
        assert int(mask.sum()) == 25
        assert int((mask & basis.live_slot_mask).sum()) == 5

        # the padding is what inflates the raw count — and it is inert,
        # because the parent already tabulates it as exactly zero.
        padding = mask & ~basis.live_slot_mask
        assert int(padding.sum()) == 20
        probe_table = basis.evaluate(_generic_unit_directions())
        assert np.array_equal(
            probe_table[:, padding], np.zeros((probe_table.shape[0], 20))
        )

    def test_a_right_angle_sample_falsely_admits_slots_from_L_four(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        r"""⛔ The vv-#13 NEGATIVE control — and it is BLIND below :math:`L = 4`.

        A finite sample of a continuous group generates a finite SUBgroup:
        four right angles generate :math:`C_4`, not :math:`SO(2)`, so every
        :math:`m \equiv 0 \pmod 4` harmonic is falsely certified invariant.

        ``[M]`` 2026-09-02, live descending slots about x under a right-angle
        sample of the SO(2)-only probe: :math:`L = 1,2,3` are **unchanged**
        (2, 3, 4 — identical to the honest answer); :math:`L = 4` read **7**
        (the true 5 plus :math:`(4,\pm4)`); :math:`L = 5` read **10** (plus
        :math:`(5,\pm4)`).

        ⭐ Since #432 (later the same day) the entry is :math:`S^2/O(2)_x` and
        its probe samples BOTH components — each right-angle rotation and its
        composition with one vertical mirror — so the sample generates
        :math:`C_{4v}`, not :math:`C_4`, and only the :math:`\sigma_v`-EVEN
        member of each :math:`\pm 4` pair is falsely admitted: ``[M]``
        :math:`L = 4` reads **6** (slot ``(4, 8)``) and :math:`L = 5` reads
        **8**. The mirrored half of the probe is doing real work; the control
        is still blind below :math:`L = 4`, which is the finding.

        ⟹ **a probe gate without an** :math:`L \ge 4` **row has an
        unfalsifiable control** — which is the whole reason this row exists
        and why it asserts the blindness at low :math:`L` explicitly rather
        than only the catch at high :math:`L`.
        """
        import orpheus.numerics.symmetry as symmetry
        from orpheus.numerics.basis.spherical_harmonic_basis import (
            SphericalHarmonicBasis,
        )

        entry = SPHERE.quotient(SubgroupOfO3.O2("x"))

        def live_count(L: int) -> int:
            basis = SphericalHarmonicBasis(L=L)
            return int((entry.descending_slots(basis) & basis.live_slot_mask).sum())

        honest = {L: live_count(L) for L in (1, 2, 3, 4, 5)}
        assert honest == {1: 2, 2: 3, 3: 4, 4: 5, 5: 6}

        monkeypatch.setattr(
            symmetry,
            "_INCOMMENSURATE_ANGLES",
            (np.pi / 2.0, np.pi, 3.0 * np.pi / 2.0, 2.0 * np.pi),
        )
        degenerate = {L: live_count(L) for L in (1, 2, 3, 4, 5)}

        assert degenerate[4] == 6 and degenerate[5] == 8, (
            f"the right-angle sample must falsely admit the sigma_v-even "
            f"m = +4 slot per degree >= 4 (C_4v, not C_4); read {degenerate}"
        )
        assert {L: degenerate[L] for L in (1, 2, 3)} == {1: 2, 2: 3, 3: 4}, (
            f"the control is BLIND below L = 4 — that is the finding, not a "
            f"defect in this row; read {degenerate}"
        )

    def test_the_mirror_leg_reproduces_the_sigma_even_mask_bit_for_bit(self) -> None:
        r"""The FINITE-group leg: for :math:`\sigma_a` the probe IS ``even_slot_mask``.

        ``[M]`` 2026-09-02 ``array_equal`` on **15 of 15** (3 axes ×
        :math:`L \in \{1..5\}`) rows, against a σ-parity classification
        written here from the definition (evaluate at generic directions and
        at their σ-images; even ⟺ equal). That reference shares no code with
        the probe, so this is a check of the probe, not of its spelling.

        NEGATIVE leg: the mask MOVES when the mirror axis moves — otherwise
        "the probe reproduces the mask" would be compatible with a mask that
        does not depend on the group at all.
        """
        from orpheus.numerics.basis.spherical_harmonic_basis import (
            MirrorEvenSphericalHarmonicBasis,
            SphericalHarmonicBasis,
        )

        probe = _generic_unit_directions()
        for axis_index, axis_name in enumerate("xyz"):
            for L in (1, 2, 3, 4, 5):
                mask = MirrorEvenSphericalHarmonicBasis(
                    L=L, mirror_axis=axis_index
                ).even_slot_mask

                mirrored = probe.copy()
                mirrored[:, axis_index] *= -1.0
                table = SphericalHarmonicBasis(L=L).evaluate(probe)
                reflected = SphericalHarmonicBasis(L=L).evaluate(mirrored)
                parity_even = np.all(
                    np.isclose(table, reflected, rtol=0.0, atol=1e-12), axis=0
                ).astype(float)

                assert np.array_equal(mask, parity_even), (
                    f"axis {axis_name}, L={L}: the entry's probe and a direct "
                    f"parity classification disagree"
                )
                assert int(mask.sum()) == {1: 5, 2: 12, 3: 22, 4: 35, 5: 51}[L]

                other = MirrorEvenSphericalHarmonicBasis(
                    L=L, mirror_axis=(axis_index + 1) % 3
                ).even_slot_mask
                assert not np.array_equal(mask, other), (
                    f"the mask must depend on the mirror AXIS (L={L})"
                )

    def test_the_probe_refuses_points_off_the_base_and_an_unsupported_base(self) -> None:
        r"""Two refusals: a probe off the base, and a base with no default probe set."""
        from orpheus.numerics.basis.spherical_harmonic_basis import (
            SphericalHarmonicBasis,
        )

        entry = SPHERE.quotient(SubgroupOfO3.O2("x"))
        with pytest.raises(ValueError, match="not on the\n?\\s*base|not on the base"):
            entry.descending_slots(
                SphericalHarmonicBasis(L=1), probe=np.array([[0.3, 0.4, 0.5]])
            )

        # positive leg: a hand-supplied probe ON the base is accepted, and it
        # gives the SAME answer as the default one (the mask is a property of
        # the group, not of the nine directions the module happens to seed).
        alternate = _generic_unit_directions(seed=31337)
        assert np.array_equal(
            entry.descending_slots(SphericalHarmonicBasis(L=4), probe=alternate),
            entry.descending_slots(SphericalHarmonicBasis(L=4)),
        )

        with pytest.raises(NotImplementedError, match="no default probe"):
            ENERGY.quotient(SubgroupOfO3.Trivial).descending_slots(
                SphericalHarmonicBasis(L=1)
            )

    def test_a_continuous_group_with_no_rotation_axis_refuses_generic_images(self) -> None:
        r"""The probe's group surface refuses what it cannot sample honestly, rather than sampling badly."""
        points = _generic_unit_directions()
        for group in (SubgroupOfO3.SO3, SubgroupOfO3.O3):
            with pytest.raises(NotImplementedError, match="no rotation axis"):
                group.generic_images(points)
        # D_inf_h HAS an axis (z) and four components: since R1 of #434
        # (2026-09-02) its realization samples the torus about z composed with
        # each coset representative — 6 angles x 4 components.
        assert len(SubgroupOfO3.Dinfh.generic_images(points)) == 24

        # positive legs, one per branch: a finite group yields every element,
        # an SO(2) yields the incommensurate rotations.
        finite = SubgroupOfO3.Mirror("y").generic_images(points)
        assert len(finite) >= 2 and all(im.shape == points.shape for im in finite)
        rotations = SubgroupOfO3.SO2("x").generic_images(points)
        assert len(rotations) == 6
        for image in rotations:
            np.testing.assert_allclose(
                np.linalg.norm(image, axis=1), 1.0, atol=1e-14
            )
            # a rotation about x fixes the x-component — the property the
            # whole descent rests on.
            np.testing.assert_allclose(image[:, 0], points[:, 0], atol=1e-14)


def _generic_unit_directions(seed: int = 20260902, n: int = 9) -> np.ndarray:
    """Deterministic generic unit directions — no component zero, none related by a coordinate symmetry."""
    raw = np.random.default_rng(seed).normal(size=(n, 3))
    return raw / np.linalg.norm(raw, axis=1)[:, None]


# ══════════════════════════════════════════════════════════════════════
# ``quotient_onto`` — the arrow a FRAME tabulates its basis through (G0)
# ══════════════════════════════════════════════════════════════════════


class TestQuotientOnto:
    r"""The G0 arrow: three admit cases and the refusal, as a function of the two manifolds alone.

    A frame binding functions on ``target`` to a rule on ``source`` is
    well-posed iff the functions can be evaluated at the rule's nodes, i.e.
    iff a quotient map ``source -> target`` exists. ``None`` is the refusal
    and it is the whole content of ERR-080's frame-level repair.
    """

    def test_equality_is_the_identity_arrow(self) -> None:
        r"""``source == target`` — the special case :math:`K = H` (the slab rule with the Legendre basis)."""
        entry = SPHERE.quotient(SubgroupOfO3.O2("x"))
        arrow = quotient_onto(entry, entry)
        assert arrow is not None
        assert arrow.domain == entry and arrow.codomain == entry

        cosines = np.linspace(-0.9, 0.9, 5)
        assert np.array_equal(arrow(cosines.reshape(-1, 1)), cosines.reshape(-1, 1))

        assert quotient_onto(SPHERE, SPHERE) is not None

    def test_a_quotient_of_the_source_is_reached_by_the_entrys_own_map(self) -> None:
        r"""``target`` is a quotient of ``source`` — the Legendre basis on a FULL-SPHERE rule.

        ⭐ This is the pairing the user asked to be buildable:
        :math:`P_\ell(\Omega\cdot\hat e_a)` is a legitimate basis on a
        Lebedev or level-symmetric rule, and the arrow that makes it
        well-posed is the entry's own quotient map.
        """
        entry = SPHERE.quotient(SubgroupOfO3.O2("x"))
        arrow = quotient_onto(SPHERE, entry)
        assert arrow is not None
        assert arrow.domain == SPHERE and arrow.codomain == entry

        directions = _generic_unit_directions()
        assert np.array_equal(arrow(directions), entry.quotient_map(directions))
        assert np.array_equal(arrow(directions), directions[:, 0])

    def test_a_finer_orbit_space_maps_onto_a_coarser_one(self) -> None:
        r"""Both sides quotients of one base with :math:`K \subseteq H` — the induced :math:`M/K \to M/H`."""
        fine = SPHERE.quotient(SubgroupOfO3.Trivial)
        coarse = SPHERE.quotient(SubgroupOfO3.O2("x"))
        assert SubgroupOfO3.SO2("x").contains(SubgroupOfO3.Trivial)

        arrow = quotient_onto(fine, coarse)
        assert arrow is not None
        directions = _generic_unit_directions()
        assert np.array_equal(arrow(directions), coarse.orbit_coordinates(directions))

    def test_the_refusal_is_the_part_one_bug_and_the_lattice_says_why(self) -> None:
        r"""``None`` — the full harmonics on the slab's orbit space (ERR-080's own pairing).

        :math:`\mathrm{Trivial} \not\supseteq SO(2)_x`, so no arrow
        :math:`S^2/O(2)_x \to S^2` exists — a coarser space cannot map onto
        a finer one, which is exactly the direction ERR-080's forged nodes
        pretended to travel.
        """
        slab_entry = SPHERE.quotient(SubgroupOfO3.O2("x"))
        fold_entry = SPHERE.quotient(SubgroupOfO3.Mirror("y"))

        assert quotient_onto(slab_entry, SPHERE) is None
        assert quotient_onto(fold_entry, SPHERE) is None
        assert quotient_onto(slab_entry, fold_entry) is None
        assert quotient_onto(SPHERE, ENERGY) is None

        # the lattice fact the refusal reads
        assert not SubgroupOfO3.Trivial.contains(SubgroupOfO3.SO2("x"))
        assert not SubgroupOfO3.SO2("x").contains(SubgroupOfO3.Mirror("y"))

    def test_the_arrow_respects_the_subgroup_ORDER_relation(self) -> None:
        r"""``vv-principles`` #15 — a predicate that must respect an order relation gets the compatibility law.

        :math:`K \subseteq H` and an arrow onto :math:`M/H` implies an arrow
        onto :math:`M/K`'s coarsening: over every lattice edge among the
        shipped angular groups, ``quotient_onto`` and ``SubgroupOfO3.contains``
        agree — neither half can be wrong alone without this reddening.
        """
        groups = (
            SubgroupOfO3.Trivial,
            SubgroupOfO3.O2("x"),
            SubgroupOfO3.O2("y"),
            SubgroupOfO3.Mirror("x"),
            SubgroupOfO3.Mirror("y"),
        )
        checked = 0
        for spent in groups:
            for has in groups:
                source = SPHERE.quotient(spent)
                target = SPHERE.quotient(has)
                admitted = quotient_onto(source, target) is not None
                assert admitted == has.contains(spent), (
                    f"quotient_onto({source.name} -> {target.name}) = "
                    f"{admitted} but {has.name}.contains({spent.name}) = "
                    f"{has.contains(spent)}"
                )
                checked += 1
        assert checked == 25


# ══════════════════════════════════════════════════════════════════════
# The entry is named by its STABILISER — #432, 2026-09-02
# ══════════════════════════════════════════════════════════════════════


class TestTheOrbitSpaceIsNamedByItsStabiliser:
    r"""One orbit space, one spelling: :math:`S^2/SO(2)_a` **is**
    :math:`S^2/O(2)_a`, so the catalogue records the larger group and refuses
    the smaller one.

    The invariants of :math:`SO(2)_a` acting on :math:`\mathbb{R}^3` are
    :math:`x_a` and :math:`x_b^2 + x_c^2`, and a vertical mirror fixes both —
    so :math:`\mathbb{R}[x]^{SO(2)_a} = \mathbb{R}[x]^{O(2)_a}`, the two
    groups have the SAME orbits (the constant-:math:`\mu` circles), and one
    derivation produces one entry. Naming it by the stabiliser is what lets
    :attr:`~orpheus.numerics.basis.base.Basis.invariance_group`, read off
    :attr:`Quotient.by`, be the FULL group a basis on the entry has
    (previously it could only be derived as the lower bound, and the
    mathematically admissible Legendre-on-a-:math:`\sigma_b`-fold pairing was
    over-refused).

    Claim layer: **term-level (L0)**, closed form — the invariant-theory
    statement above, checked against the group action itself.
    """

    @pytest.mark.parametrize("axis", ["x", "y", "z"])
    def test_the_entry_records_the_stabiliser_and_its_two_derivation_outputs(
        self, axis: str
    ) -> None:
        r"""The POSITIVE leg: ``S^2/O2_a`` exists on every axis, carries
        :math:`O(2)_a` as its ``by``, the cosine interval as its realization,
        and Archimedes' hat-box measure (``LEGENDRE``) as its pushforward
        reference."""
        from orpheus.numerics.generating_measure import LEGENDRE
        from orpheus.numerics.manifold import AXIS_INDEX

        stabiliser = SubgroupOfO3.O2(axis)
        entry = SPHERE.quotient(stabiliser)
        assert entry.by == stabiliser
        assert entry.name == f"S^2/O2_{axis}"
        assert entry.realization == COSINE_INTERVAL
        assert entry.reference is LEGENDRE
        assert entry.dim == 1
        # the surviving invariant IS the cosine against THIS axis
        directions = _generic_unit_directions()
        assert np.array_equal(
            entry.orbit_coordinates(directions), directions[:, AXIS_INDEX[axis]]
        )
        # three axes, three DIFFERENT entries — same numbers, different maps
        others = [SPHERE.quotient(SubgroupOfO3.O2(b)) for b in "xyz" if b != axis]
        assert all(entry != other for other in others)

    @pytest.mark.parametrize("axis", ["x", "y", "z"])
    def test_the_quotient_by_the_ROTATION_HALF_is_refused_with_the_theorem(
        self, axis: str
    ) -> None:
        r"""The NEGATIVE leg (``vv-principles`` #11): asking for
        :math:`S^2/SO(2)_a` raises, and the message carries the reason rather
        than merely the refusal — it names the stabiliser, the invariants, and
        the spelling to use.

        ⚠ The pinned fragment is the message's HEAD,
        ``"S^2/SO2_a is the orbit space S^2/O2_a"`` — it names both orbit
        spaces and is the sentence the refusal exists to say, so it is the
        piece that must survive wherever the check is enforced from (the
        derivation, the catalogue door, or the type's own construction
        invariant). The rest of the message is a diagnosis and is
        deliberately NOT pinned, so re-wording it is not a false red.

        ⭐ The second leg is a PLACEMENT claim, and it is the one with teeth
        after the check moves: the refusal must be a ``ValueError`` carrying
        that sentence and NOT the catalogue's ``NotImplementedError``
        ("no catalogue entry", which names the derivation WORK). Those are
        the two ways ``S^2/SO2_a`` can fail to build, and only one of them is
        the theorem — so a check placed AFTER the catalogue lookup, or a key
        merely deleted from the table, reddens here.
        """
        rotation_half = SubgroupOfO3.SO2(axis)
        with pytest.raises(ValueError) as excinfo:
            SPHERE.quotient(rotation_half)
        message = str(excinfo.value)
        assert f"S^2/SO2_{axis} is the orbit space S^2/O2_{axis}" in message, (
            message
        )
        assert "no catalogue entry" not in message
        assert not isinstance(excinfo.value, NotImplementedError)

    def test_the_entry_is_named_by_the_LARGEST_group_with_its_orbits(self) -> None:
        r"""⭐ The theorem as a test: for every group the lattice can spell,

        .. math:: G \subseteq \mathtt{entry.by}
                  \iff \text{every generic image of } G
                  \text{ leaves } \mathtt{orbit\_coordinates} \text{ unchanged.}

        The right-hand side is what "these are the orbits" MEANS — the orbit
        coordinates are constant on the fibres of :math:`\pi` — and the
        left-hand side is the lattice's answer. Neither half can be wrong
        alone without this reddening (``vv-principles`` #15), and it is the
        maximality claim: were ``by`` set to a proper subgroup, some group
        outside it would still fix the invariants and the ⟸ direction would
        fail; were it set to something too big, an element would move them and
        ⟹ would fail.

        `[M]` 2026-09-02: **140** (entry, group) pairs over the seven shipped
        sphere entries × the 20 groups whose generic images can be sampled —
        **0** mismatches, **33** pairs inside and **107** outside, so BOTH
        directions are populated (`vv-principles` #20). Cost 0.74 s.

        ⚠ The denominator EXCLUDES :math:`D_{\infty h}`, :math:`SO(3)` and
        :math:`O(3)`: they are continuous with no axis to sample about and
        :meth:`generic_images` refuses them by design. Their absence is
        asserted, so the exclusion cannot silently widen.
        """
        entries = (
            [SubgroupOfO3.O2(a) for a in "xyz"]
            + [SubgroupOfO3.Mirror(a) for a in "xyz"]
            + [SubgroupOfO3.Trivial]
        )
        groups = (
            [
                SubgroupOfO3.Trivial, SubgroupOfO3.OctahedralOh,
                SubgroupOfO3.IcosahedralIh,
            ]
            + [SubgroupOfO3.Mirror(a) for a in "xyz"]
            + [SubgroupOfO3.SO2(a) for a in "xyz"]
            + [SubgroupOfO3.O2(a) for a in "xyz"]
            + [SubgroupOfO3.Cn(n) for n in (1, 2, 3, 4)]
            + [SubgroupOfO3.Dnh(n) for n in (1, 2, 3, 4)]
        )
        assert len(entries) == 7 and len(groups) == 20

        points = _generic_unit_directions()
        mismatches: list[str] = []
        inside = outside = 0
        for spent in entries:
            entry = SPHERE.quotient(spent)
            reference = entry.orbit_coordinates(points)
            for group in groups:
                fixes = all(
                    np.allclose(entry.orbit_coordinates(image), reference, atol=1e-12)
                    for image in group.generic_images(points)
                )
                contained = entry.by.contains(group)
                inside += int(contained)
                outside += int(not contained)
                if fixes != contained:
                    mismatches.append(
                        f"{entry.name}: {group.name} fixes the invariants="
                        f"{fixes} but by.contains={contained}"
                    )
        assert not mismatches, mismatches[:8]
        assert inside + outside == 140
        assert inside == 33 and outside == 107, (inside, outside)

        # the three excluded groups, asserted rather than assumed
        for axis_free in (SubgroupOfO3.SO3, SubgroupOfO3.O3):
            with pytest.raises(NotImplementedError, match="no rotation axis"):
                axis_free.generic_images(points)

    def test_the_stabiliser_of_the_axial_entry_is_STRICTLY_bigger_than_the_rotation_half(
        self,
    ) -> None:
        r"""Non-vacuity for the row above: the maximality claim only has
        content if the stabiliser is genuinely larger than the group whose
        derivation produced the entry.

        `[M]` 2026-09-02 on ``S^2/O2_x``: **7** of the 20 spellable groups sit
        inside :math:`O(2)_x` — ``Trivial``, ``sigma_y``, ``sigma_z``,
        ``SO2_x``, ``O2_x``, ``C_1``, ``D_1h`` — while only **3** sit inside
        the rotation half :math:`SO(2)_x` (``Trivial``, ``SO2_x``, ``C_1``;
        the two vertical mirrors and :math:`D_{1h}` are improper, and
        :math:`O(2)_x` itself is strictly bigger). So naming the entry by the
        rotation half would lose **four** edges, of which ``sigma_y`` and
        ``sigma_z`` are exactly what admits the Legendre basis on a
        :math:`\sigma_b`-fold.
        """
        groups = (
            [
                SubgroupOfO3.Trivial, SubgroupOfO3.OctahedralOh,
                SubgroupOfO3.IcosahedralIh,
            ]
            + [SubgroupOfO3.Mirror(a) for a in "xyz"]
            + [SubgroupOfO3.SO2(a) for a in "xyz"]
            + [SubgroupOfO3.O2(a) for a in "xyz"]
            + [SubgroupOfO3.Cn(n) for n in (1, 2, 3, 4)]
            + [SubgroupOfO3.Dnh(n) for n in (1, 2, 3, 4)]
        )
        stabiliser_side = {
            g.name for g in groups if SubgroupOfO3.O2("x").contains(g)
        }
        rotation_side = {
            g.name for g in groups if SubgroupOfO3.SO2("x").contains(g)
        }
        # C_1 IS Trivial since 2026-09-02 (R1 of #434), so it is not a second name here.
        assert stabiliser_side == {
            "Trivial", "sigma_y", "sigma_z", "SO2_x", "O2_x", "D_1h",
        }
        assert rotation_side == {"Trivial", "SO2_x"}
        assert rotation_side < stabiliser_side
        assert stabiliser_side - rotation_side == {
            "O2_x", "sigma_y", "sigma_z", "D_1h",
        }

    # ── the invariant lives on the TYPE, not in one builder ──────────

    def test_replace_cannot_forge_a_second_spelling(self) -> None:
        r"""``Quotient.__post_init__`` reads ``by.orbit_stabiliser``.

        `[M]` 2026-09-02, while the refusal lived inside ``_sphere_mod_o2``
        alone: ``replace(entry, by=SO2("x"))`` constructed, compared unequal
        to the catalogue entry under ``==``, and was accepted by
        ``barycentre`` and the polar-axis reader — ERR-080's own defect
        class (one orbit space, two objects) through the idiom Pattern 4
        blesses (elegance review).  Refused where it is written.
        """
        from dataclasses import replace

        entry = SPHERE.quotient(SubgroupOfO3.O2("x"))
        with pytest.raises(ValueError, match=r"is the orbit space S\^2/O2_x"):
            replace(entry, by=SubgroupOfO3.SO2("x"))

    def test_the_door_names_the_spelling_and_SO3_is_refused_the_same_way(
        self,
    ) -> None:
        r"""The door's message ends in the spelling to use, and the same law
        refuses the other connected group: :math:`S^2/SO(3)` is the point
        :math:`S^2/O(3)` (``orbit_stabiliser`` is one rule, not an axial
        special case)."""
        with pytest.raises(ValueError, match=r"spell SubgroupOfO3\.O2\('y'\)"):
            SPHERE.quotient(SubgroupOfO3.SO2("y"))
        with pytest.raises(ValueError, match=r"S\^2/SO3 is the orbit space S\^2/O3"):
            SPHERE.quotient(SubgroupOfO3.SO3)

    def test_every_catalogue_key_is_OBTAINABLE(self) -> None:
        r"""The door's "catalogued today" listing is exactly what can be
        built: no key routes to a refusal.  `[M]` until 2026-09-02 three of
        nine did (``SO2_a`` decoys kept to reach the in-builder refusal), and
        the ``NotImplementedError`` for an uncatalogued group advertised them
        as entries.
        """
        from orpheus.numerics.manifold import _ORBIT_CATALOGUE

        names = sorted(g for _, g in _ORBIT_CATALOGUE)
        assert names == ["O2_x", "O2_y", "O2_z", "sigma_x", "sigma_y", "sigma_z"]
        with pytest.raises(NotImplementedError) as excinfo:
            SPHERE.quotient(SubgroupOfO3.Cn(4))
        for name in names:
            assert f"Sphere/{name}" in str(excinfo.value)
        assert "SO2" not in str(excinfo.value)


# =============================================================================
# 2.2b — the Γ-slot: gates drafted by the test-architect (2026-09-02)
# =============================================================================

_g22b_AXES = ("x", "y", "z")


def _g22b_trivially_quotiented_product() -> DiscreteMeasure:
    from dataclasses import replace

    return replace(
        Quadrature.product(4, 8).measure,
        support=SPHERE.quotient(SubgroupOfO3.Trivial),
    )




def _g22b_rot(axis: str, theta: float) -> RigidMotion:
    return RigidMotion.rotation_about_axis(axis=np.eye(3)["xyz".index(axis)], angle=theta)


def _g22b_mirror(axis: str) -> RigidMotion:
    return RigidMotion.reflection(normal=np.eye(3)["xyz".index(axis)])


#: Pairwise-incommensurate angles.  A right-angle sample of a continuous
#: family generates C_4 and certifies what it never tested (ERR-072,
#: ``vv-principles`` #13); the CONTROL below is allowed to sample because it
#: can only REFUTE, and it is compared against a criterion that never does.
_g22b_INCOMMENSURATE = (1.0, float(np.sqrt(2.0)), float(np.e), 2.5, float(np.sqrt(7.0)))


# =============================================================================


# G1 — tests/numerics/test_manifold.py
# =============================================================================
class TestTheLiftIsARightInverseOfTheQuotientMap:
    r"""(ii) ``π ∘ lift = id`` on the chart, for all three entry families.

    The lift is what lets an orbit-space question be asked of a measure whose
    nodes arrive in CHART coordinates (a slab rule carries μ, not a point of
    S²).  Its defining property is the only one an induced action needs, and
    it is EXACT on every shipped family — `[M]` 2026-09-02, ``array_equal``
    True on 200 disk points, 201 μ values and 200 sphere points.

    ⭐ Bit-exactness here is a THEOREM, not a draw (``vv-principles`` #31):
    every shipped ``orbit_coordinates`` is a COLUMN SELECTION, and the lift
    writes exactly the selected columns back, so no arithmetic happens on the
    surviving coordinates at all.  The √ that the hemisphere section computes
    lands only in the DROPPED column.
    """

    @pytest.mark.foundation
    @pytest.mark.parametrize("axis", _g22b_AXES)
    def test_the_axial_entrys_lift_is_the_barycentre_and_pi_undoes_it(self, axis):
        entry = SPHERE.quotient(SubgroupOfO3.O2(axis))
        mu = np.linspace(-1.0, 1.0, 201).reshape(-1, 1)
        lifted = np.asarray(entry.lift(mu), dtype=float)
        # the barycentre μ ê_a — the orbit's MEAN, inside the ball
        expected = np.zeros((201, 3))
        expected[:, "xyz".index(axis)] = mu[:, 0]
        assert np.array_equal(lifted, expected)
        assert np.array_equal(
            np.asarray(entry.orbit_coordinates(lifted)).reshape(-1), mu[:, 0]
        )

    # ``test_a_mirror_entrys_lift_is_the_hemisphere_section_and_pi_undoes_it``
    # stood here until R4 of #434 (2026-09-03): the mirror entry's lift is the
    # orbit BARYCENTRE now, like the axial one, and lands in the ball — see
    # ``TestR4TheCoordinateChartPairIsTheReynoldsProjector`` (lift ∘ chart == P_H, with the
    # round-trip leg kept and labelled blind) and
    # ``TestR4TheLiftIsADerivationOutputNotATagBranch`` (the codomain).

    @pytest.mark.foundation
    def test_the_trivial_entrys_lift_is_the_identity(self):
        entry = SPHERE.quotient(SubgroupOfO3.Trivial)
        rng = np.random.default_rng(7)
        pts = rng.normal(size=(200, 3))
        pts /= np.linalg.norm(pts, axis=1)[:, None]
        assert np.array_equal(np.asarray(entry.lift(pts), dtype=float), pts)
        assert np.array_equal(np.asarray(entry.orbit_coordinates(pts)), pts)

    @pytest.mark.foundation
    def test_the_sphere_lifts_land_in_the_BALL_and_say_so(self):
        r"""⚠ A sphere entry's lift is NOT a section, and its codomain is the
        honest record of that: `[M]` ``barycentre(S^2/O2_x).codomain`` and
        ``barycentre(S^2/sigma_y).codomain`` are both ``D^3`` (the mirror
        entry joined the axial one at R4 of #434, 2026-09-03 — until then it
        lifted through a hemisphere SECTION onto ``S^2``, and this row's
        docstring argued that one uniform codomain across the three families
        would be false; that argument is REPEALED: only the trivial entry's
        lift is a section now).  The ERR-080 discriminator is unchanged: a map
        into the ball declared as a map into the sphere is the forgery."""
        axial = SPHERE.quotient(SubgroupOfO3.O2("x"))
        assert isinstance(axial.lift.codomain, Ball)
        assert axial.lift.codomain == barycentre(axial).codomain
        assert isinstance(SPHERE.quotient(SubgroupOfO3.Mirror("y")).lift.codomain, Ball)
        trivial = SPHERE.quotient(SubgroupOfO3.Trivial)
        assert trivial.lift.codomain == trivial.base
        # and the axial image really is off the sphere away from the poles
        img = np.asarray(axial.lift(np.array([[0.3]])), dtype=float)
        assert not SPHERE.contains(img)
        assert SPHERE.contains(np.asarray(axial.lift(np.array([[1.0]])), dtype=float))


# =============================================================================


# G2 — tests/numerics/test_manifold.py
# =============================================================================
class TestTheInducedActionIsWellDefinedAndActsAsTheTheoremSays:
    r"""(i) ``Quotient.induced_action`` — REFUSE what does not normalise, and
    act as the algebra says on what does.

    ⛔ **§6c note, and the coordinator must read it before crediting this
    gate.**  The refusal is provably UNREACHABLE from
    :meth:`~orpheus.numerics.measure.DiscreteMeasure.is_invariant_under`
    (``SubgroupOfO3.is_invariant`` until R2 of #434): the orbit-space kernel asks
    ``G.normalises(H)`` at step 1 and only then applies ``induced_action`` to
    elements of ``G``, every one of which normalises ``H`` because ``G``
    does.  So its only witness is a DIRECT call — which is what this class
    is.  Do not add an end-to-end row for it; there is none to add.
    """

    @pytest.mark.foundation
    def test_a_motion_that_does_not_normalise_the_group_is_REFUSED(self):
        fold = SPHERE.quotient(SubgroupOfO3.Mirror("y"))
        c4z = _g22b_rot("z", np.pi / 2.0)
        # `[M]` C_4 about z conjugates sigma_y to sigma_x — measured, the
        # conjugated linear part is diag(-1, +1, +1) to 2.2e-16.
        conj = _g22b_mirror("y").conjugated_by(c4z)
        assert np.allclose(conj.linear, np.diag([-1.0, 1.0, 1.0]), atol=1e-12)
        with pytest.raises(ValueError, match="normalise"):
            fold.induced_action(c4z)

    @pytest.mark.foundation
    def test_a_rotation_about_another_axis_does_not_descend_to_an_AXIAL_entry(self):
        entry = SPHERE.quotient(SubgroupOfO3.O2("x"))
        with pytest.raises(ValueError, match="normalise"):
            entry.induced_action(_g22b_rot("y", 1.0))
        # the positive control: about its OWN axis it descends, and is the
        # IDENTITY on the chart (every orbit is fixed as a set).
        act = entry.induced_action(_g22b_rot("x", 1.0))
        mu = np.linspace(-1.0, 1.0, 21).reshape(-1, 1)
        assert np.allclose(np.asarray(act(mu)).reshape(-1), mu[:, 0], atol=1e-15)

    @pytest.mark.foundation
    def test_sigma_x_on_the_sigma_y_fold_is_the_disks_x_flip(self):
        fold = SPHERE.quotient(SubgroupOfO3.Mirror("y"))
        rng = np.random.default_rng(11)
        rho = np.sqrt(rng.uniform(0.0, 1.0, 64))
        phi = rng.uniform(0.0, 2.0 * np.pi, 64)
        chart = np.column_stack([rho * np.cos(phi), rho * np.sin(phi)])
        got = np.asarray(fold.induced_action(_g22b_mirror("x"))(chart), dtype=float)
        assert np.allclose(got, chart * np.array([-1.0, 1.0]), atol=1e-15)
        got_z = np.asarray(fold.induced_action(_g22b_mirror("z"))(chart), dtype=float)
        assert np.allclose(got_z, chart * np.array([1.0, -1.0]), atol=1e-15)
        # and the quotienting mirror itself acts TRIVIALLY — the whole point.
        got_y = np.asarray(fold.induced_action(_g22b_mirror("y"))(chart), dtype=float)
        assert np.allclose(got_y, chart, atol=1e-15)

    @pytest.mark.foundation
    @pytest.mark.parametrize("axis", _g22b_AXES)
    def test_the_axis_flipping_mirror_is_mu_to_minus_mu_THROUGH_the_barycentre(self, axis):
        entry = SPHERE.quotient(SubgroupOfO3.O2(axis))
        mu = np.linspace(-1.0, 1.0, 41).reshape(-1, 1)
        got = np.asarray(entry.induced_action(_g22b_mirror(axis))(mu), dtype=float)
        assert np.array_equal(got.reshape(-1), -mu[:, 0])

    @pytest.mark.foundation
    def test_the_induced_action_is_a_GROUP_HOMOMORPHISM_on_the_chart(self):
        r"""``ind(g) ∘ ind(h) = ind(g·h)`` — the law that separates a correct
        induced action from one applied to the AMBIENT nodes.

        `[M]` 2026-09-02: max |difference| **0.000e+00** over 128 pairs
        (2 entries × 8 × 8 elements of D_2h).  Exact because every D_2h
        element is a signed permutation and every chart is a column
        selection; a fold's chart image and its chart nodes are then bit-equal
        rather than merely close.
        """
        elements = _closed(SubgroupOfO3.Dnh(2))
        for entry, chart in (
            (SPHERE.quotient(SubgroupOfO3.Mirror("y")), _disk_points()),
            (SPHERE.quotient(SubgroupOfO3.O2("x")), np.linspace(-1, 1, 17).reshape(-1, 1)),
        ):
            for g in elements:
                for h in elements:
                    lhs = np.asarray(
                        entry.induced_action(g)(entry.induced_action(h)(chart)),
                        dtype=float,
                    )
                    rhs = np.asarray(entry.induced_action(g @ h)(chart), dtype=float)
                    assert np.allclose(lhs, rhs, atol=1e-14), (g, h)


def _disk_points() -> np.ndarray:
    rng = np.random.default_rng(3)
    rho = np.sqrt(rng.uniform(0.0, 1.0, 32))
    phi = rng.uniform(0.0, 2.0 * np.pi, 32)
    return np.column_stack([rho * np.cos(phi), rho * np.sin(phi)])


def _closed(group: SubgroupOfO3) -> list[RigidMotion]:
    from orpheus.numerics.symmetry import _group_elements

    elems = _group_elements(group._tag)
    assert elems is not None
    return list(elems)


# =============================================================================


# G3 — tests/numerics/test_manifold.py
# =============================================================================
# ``TestSectionCoordinatesDispatchOnWidthAndRefuseAnythingElse`` stood here until
# R4 of #434 (2026-09-03): ``ambient_representatives`` passed ambient-width points
# through as representatives, and that pass-through is what R4 retired — the
# method is ``orbit_barycentres`` now and maps BOTH widths to the orbit
# barycentre.  Its rows are superseded by ``TestR4OrbitBarycentresIsOneConceptOnBothWidths``
# (both widths; the pass-through gone; the width refusal under the new name).


# =============================================================================


# G4 — tests/numerics/test_manifold.py
# =============================================================================
class TestTheEmbeddingReadsTheLiftRatherThanSpellingItTwice:
    r"""⛔ **The design's R3 says the axial arm of ``_embedded_nodes``
    "RETIRES".  Deleting it is a SILENT WRONG ANSWER in a second consumer.**

    `[M]` 2026-09-02.  ``_embedded_nodes`` has TWO production readers, and
    only one of them is the invariance kernel: ``Quadrature.ordinate_permutation``
    (``directional.py:539``) feeds it to ``RigidMotion.preserves``.  With the
    axial arm deleted, a rule on ``S^2/O2_y`` embeds as ``(μ, 0, 0)`` instead
    of ``(0, μ, 0)`` — max |difference| **9.603e-01** — and
    ``ordinate_permutation(σ_x)`` answers the μ→−μ permutation
    ``[7, 6, 5, 4]`` where the honest answer is the IDENTITY ``[0, 1, 2, 3]``,
    with σ_y flipping the other way.  Two of the three axes are wrong; the
    slab's own ``S^2/O2_x`` is the one axis where the two agree, which is why
    ``symmetry.py``'s own suite cannot see it.

    ⟹ the arm is RE-POINTED, not deleted: ``_embedded_nodes``'s axial branch
    reads ``measure.support.lift(nodes)``.  `[M]` that is ``array_equal`` to
    today's answer on all three axes, so the re-point is bit-identical AND
    single-sources the lift (Pattern 2).
    """

    @pytest.mark.foundation
    @pytest.mark.parametrize("axis", _g22b_AXES)
    def test_the_axial_embedding_IS_the_entrys_lift(self, axis):
        from orpheus.numerics.invariance import _embedded_nodes

        rule = gauss_legendre_on_polar_orbit(8, axis)
        nodes = np.asarray(rule.nodes, dtype=float).reshape(-1, 1)
        assert np.array_equal(
            np.asarray(_embedded_nodes(rule), dtype=float),
            np.asarray(rule.support.lift(nodes), dtype=float),
        )

    @pytest.mark.foundation
    @pytest.mark.parametrize("axis", _g22b_AXES)
    def test_ordinate_permutation_still_reads_the_axis_the_orbit_space_names(self, axis):
        """The consumer a deleted axial arm would corrupt, pinned directly.

        `[M]` the identity/flip pattern is diagonal in the axis: only the
        mirror NORMAL to the orbit space's own axis moves the ordinates."""
        rule = gauss_legendre_on_polar_orbit(8, axis)
        quad = Quadrature(measure=rule)
        for other in _g22b_AXES:
            perm = quad.ordinate_permutation(_g22b_mirror(other))
            assert perm is not None
            moved = perm.indices.tolist() != list(range(rule.n_points))
            assert moved is (other == axis), (axis, other)

    @pytest.mark.foundation
    def test_a_BARE_interval_keeps_the_column_zero_convention(self):
        """The arm that does NOT retire: a bare interval names no axis, so
        column 0 is the convention and there is no lift to read."""
        from orpheus.numerics.invariance import _embedded_nodes

        m = gauss_legendre_on_mu(8)
        assert m.support == COSINE_INTERVAL
        emb = np.asarray(_embedded_nodes(m), dtype=float)
        assert np.array_equal(emb[:, 0], np.asarray(m.nodes, dtype=float))
        assert np.array_equal(emb[:, 1:], np.zeros((m.n_points, 2)))


# =============================================================================


# =============================================================================
# ===== R4 gates (test-architect, 2026-09-03) — the lift is a derivation output
# =============================================================================

_R4_I3 = np.eye(3)
_R4_ATOL_EXACT = 1e-14          # C1 / B2: a 30x margin over the measured floors
_R4_ATOL_NEAR_BIT = 1e-13       # C3's non-signed-permutation leg (measured 5.0e-16)
_R4_O1 = 0.1                    # C2: the negative leg's floor (measured >= 5.5e-01)


# ─────────────────────────────────────────────────────────────────────
# References — every one computed from the group's REALIZATION
# ─────────────────────────────────────────────────────────────────────


def _r4_reynolds_projector(group: SubgroupOfO3) -> np.ndarray:
    r"""``P_H``, the orthogonal projector of :math:`\mathbb{R}^3` onto the
    fixed subspace :math:`(\mathbb{R}^3)^H` — the INDEPENDENT reference.

    Built from the realization alone: an orthonormal basis ``B`` of the common
    null space of the Lie generators (:math:`Xp = 0`) and of :math:`r - I` for
    every coset representative, then :math:`P = B B^{\mathsf T}`.  ``[D]``
    :math:`p \in (\mathbb{R}^3)^H` iff every element fixes it, iff the
    identity component's algebra kills it AND each representative fixes it —
    the component being connected, its own condition is exactly ``Xp = 0``.

    It reads NO column index, which is what makes it independent of
    ``_coordinate_chart`` (see the module docstring).
    """
    realization = group.realization
    rows = [np.asarray(X, dtype=float) for X in realization.component.generators]
    rows += [np.asarray(m.linear, dtype=float) - _R4_I3
             for m in realization.representatives]
    stacked = np.vstack(rows) if rows else np.zeros((0, 3))
    if stacked.shape[0] == 0:
        return _R4_I3.copy()
    _, singular, vt = np.linalg.svd(stacked)
    rank = int((singular > 1e-9).sum())
    basis = vt[rank:].T if rank < 3 else np.zeros((3, 0))
    return basis @ basis.T


def _r4_orbit_dim(group: SubgroupOfO3, point: np.ndarray) -> int:
    r"""``dim(H·p) = rank[X p : X ∈ 𝔤]`` — the INDEPENDENT reference for the
    dimension law.

    The orbit of a point under a Lie group is the image of the group, and its
    tangent space at ``p`` is spanned by ``{X p}`` over a basis of the algebra.
    A finite group has no algebra, so the rank is 0 whatever the point.
    """
    generators = group.realization.component.generators
    if not generators:
        return 0
    stacked = np.stack([np.asarray(X, dtype=float) @ point for X in generators])
    return int(np.linalg.matrix_rank(stacked, tol=1e-9))


def _r4_group_mean(group: SubgroupOfO3, points: np.ndarray) -> np.ndarray:
    """``(1/|H|) Σ_{g∈H} g p`` — the Reynolds average, over the group's OWN
    element list, with no manifold code in the room.

    Spelled as a stacked ``mean`` rather than ``sum(...)/len`` because the
    latter types as ``NDArray | float`` (an empty ``sum`` is ``0``) and would
    need a ``# type: ignore`` — ``coding-elegance`` anti-pattern #19: the
    principled spelling is also the one that reads like the average it is.
    ``[M]`` 2026-09-03 the two are ``array_equal`` on all four finite entries,
    so the bit tier below is unaffected by the choice.
    """
    stacked = np.stack([
        points @ np.asarray(m.linear, dtype=float).T
        for m in group.realization.elements
    ])
    return stacked.mean(axis=0)


def _r4_circle_mean(axis: int, points: np.ndarray, n: int = 16) -> np.ndarray:
    r"""The mean of ``p`` over the circle :math:`\{R_\theta p\}` about ``axis``,
    by an ``n``-point trapezoid.

    ⚠ ``n = 16`` is a MEASURED choice and MORE IS WORSE — the trapezoid
    integrates :math:`\cos\theta` and :math:`\sin\theta` exactly for
    :math:`n \ge 3`, so everything past that is summation error.  ``[M]``
    2026-09-03, residual against ``P_H p``: ``n=8 → 2.220e-16``,
    ``16 → 3.331e-16``, ``32 → 5.551e-16``, ``64 → 1.110e-15``,
    ``1024 → 2.587e-14``.  Do not "strengthen" this row by raising ``n``.
    """
    accumulated = np.zeros_like(points)
    for k in range(n):
        rotation = RigidMotion.rotation_about_axis(
            axis=_R4_I3[axis], angle=2.0 * np.pi * k / n,
        )
        accumulated += points @ np.asarray(rotation.linear, dtype=float).T
    return accumulated / n


# ─────────────────────────────────────────────────────────────────────
# Fixtures — parametrize over a LABEL, build in the body
# ─────────────────────────────────────────────────────────────────────

#: Every orbit space the tree can construct whose base has ambient width 3 —
#: the 6 catalogue entries plus the two shipped trivial ones (``R^3/{e}`` IS
#: ``invariance._ambient_orbit_space()``, a production object).  The FINITE
#: roster is probed whole, never sampled (``vv-principles`` #31's finite-roster
#: corollary; #13's ladder rule governs unbounded families, and this is not one).
_R4_ENTRY_LABELS = (
    "O2_x", "O2_y", "O2_z",
    "sigma_x", "sigma_y", "sigma_z",
    "S2/Trivial", "R3/Trivial",
)

_R4_MIRROR_LABELS = ("sigma_x", "sigma_y", "sigma_z")
_R4_AXIAL_LABELS = ("O2_x", "O2_y", "O2_z")
_R4_FOLDS = ((2, 4), (4, 6), (4, 8), (8, 8))


def _r4_entry(label: str) -> Quotient:
    """The orbit space named by ``label`` — built in the test body."""
    if label == "S2/Trivial":
        return SPHERE.quotient(SubgroupOfO3.Trivial)
    if label == "R3/Trivial":
        return RealSpace(3).quotient(SubgroupOfO3.Trivial)
    family, axis = label.split("_")
    group = (SubgroupOfO3.O2(axis) if family == "O2"
             else SubgroupOfO3.Mirror(axis))
    return SPHERE.quotient(group)


def _r4_sphere_points(n: int = 41, seed: int = 20260903) -> np.ndarray:
    """Seeded unit vectors — generic, so no coordinate is accidentally zero."""
    rng = np.random.default_rng(seed)
    points = rng.standard_normal((n, 3))
    return points / np.linalg.norm(points, axis=1, keepdims=True)


def _r4_disk_points(n: int = 64, seed: int = 11) -> np.ndarray:
    """Seeded points of the closed unit disk — a mirror entry's chart."""
    rng = np.random.default_rng(seed)
    rho = np.sqrt(rng.uniform(0.0, 1.0, n))
    phi = rng.uniform(0.0, 2.0 * np.pi, n)
    return np.column_stack([rho * np.cos(phi), rho * np.sin(phi)])


def _r4_chart_of(entry: Quotient, points: np.ndarray) -> np.ndarray:
    """The entry's chart image as a 2-D column block, whatever its width."""
    chart = np.asarray(entry.orbit_coordinates(points), dtype=float)
    return chart if chart.ndim == 2 else chart[:, None]


def _r4_generic_point(entry: Quotient) -> np.ndarray:
    """A generic point of the entry's base — the one the dimension law reads."""
    if entry.base == SPHERE:
        return np.array(
            [0.31622776601683794, 0.5477225575051661, 0.7745966692414834]
        )
    return np.array([0.3, -0.7, 1.1])


def _r4_forged(kind: str) -> dict:
    r"""The two B2 forgeries, as constructor kwargs.

    ``[M]`` 2026-09-03 on the PRE-carve tree BOTH construct — that is this
    plan's §6c witness that the dimension law lands with the case it catches.
    Both carry ``fundamental_domain=None``, so the EXISTING fd clause returns
    early and provably cannot see them.
    """
    if kind == "axial_on_a_disk":
        # S^2/O(2)_z: 2 - 1 = 1, offered a 2-dimensional realization.
        group, realization, columns = SubgroupOfO3.O2("z"), Ball(2), [0, 1]
    else:
        # S^2/sigma_x: 2 - 0 = 2, offered a 1-dimensional realization.
        group, realization, columns = SubgroupOfO3.Mirror("x"), COSINE_INTERVAL, 1
    select, embed = _coordinate_chart(columns, ambient_dim(SPHERE))
    return {
        "base": SPHERE, "by": group, "realization": realization,
        "orbit_coordinates": select,
        "lift_coordinates": embed, "lift_codomain": Ball(3),
    }


def _r4_hand_built_octahedral() -> dict:
    r"""A LEGAL entry outside all three retired ``lift`` arms — the §6c witness
    that retiring the tag branch is a capability, not merely "no test broke".

    ``[M]`` 2026-09-03 on the PRE-carve tree it CONSTRUCTS (``O_h`` is its own
    ``orbit_stabiliser``, and the dimension law reads ``2 − 0 = 2`` ✓) and
    ``entry.lift`` raises ``NotImplementedError: no lift is spelled for S^2/Oh``.
    Its chart is a placeholder — ``S^2/O_h``'s real invariants are higher-degree
    polynomials (#435) — and the row asserts only that the FIELD is read, never
    that the chart is the derivation's.
    """
    select, embed = _coordinate_chart([0, 1], ambient_dim(SPHERE))
    return {
        "base": SPHERE, "by": SubgroupOfO3.OctahedralOh, "realization": Ball(2),
        "orbit_coordinates": select,
        "lift_coordinates": embed, "lift_codomain": Ball(3),
    }


# =====================================================================
# A — the chart pair IS the Reynolds projector (finding B4)
# =====================================================================
class TestR4TheCoordinateChartPairIsTheReynoldsProjector:
    r"""``embed ∘ select = P_H``, the orthogonal projector onto
    :math:`(\mathbb{R}^3)^H`.

    This is the carve's central claim and the reason ONE helper can spell the
    lift for every column-selection entry: the axial entry's
    :math:`\mu \mapsto \mu\,\hat e_a`, the mirror entry's
    :math:`(x_b, x_c) \mapsto (0, x_b, x_c)` and the trivial entry's identity
    are the SAME map — the Reynolds average of the group — read at three
    different fixed subspaces.

    ⭐ **Bit-identity here is a THEOREM, not a draw** (``vv-principles`` #31,
    and the finite-roster corollary: the population is 8 entries and every one
    is probed).  ``select`` is a column read, ``embed`` writes those floats
    into a zero array, and ``[M]`` ``P_H`` is a 0/1 DIAGONAL on all 8 shipped
    entries, so ``p @ P.T`` re-reads the same floats.  Measured
    ``max|emb(sel p) − P p| = 0.000e+00``, ``array_equal`` on 8 of 8.
    ⚠ The bit tier is a claim about the SHIPPED entries: an entry whose ``H``
    is not axis-aligned has a dense ``P_H`` and belongs at
    ``assert_array_almost_equal_nulp``.  ``test_the_projector_is_diagonal_here``
    pins that PREMISE, so the day it stops holding the suite says which claim
    it lost rather than reddening a row whose subject is elsewhere.
    """

    @pytest.mark.foundation
    @pytest.mark.parametrize("label", _R4_ENTRY_LABELS)
    def test_embed_after_select_is_the_reynolds_projector(self, label: str) -> None:
        entry = _r4_entry(label)
        columns = _r4_columns_of(entry)
        select, embed = _coordinate_chart(columns, ambient_dim(entry.base))
        points = _r4_sphere_points()
        projector = _r4_reynolds_projector(entry.by)

        chart = np.asarray(select(points), dtype=float)
        lifted = np.asarray(embed(chart), dtype=float)
        np.testing.assert_array_equal(
            lifted, points @ projector.T,
            err_msg=(
                f"{entry.name}: embed(select(p)) is not the Reynolds projector "
                f"onto (R^3)^H — max|D| = "
                f"{np.abs(lifted - points @ projector.T).max():.3e}"
            ),
        )

    @pytest.mark.foundation
    @pytest.mark.parametrize("label", _R4_ENTRY_LABELS)
    def test_select_is_bit_identical_to_the_entrys_own_orbit_coordinates(
        self, label: str,
    ) -> None:
        """The frame-table landmine (plan Landmine 2): the chart must not move.

        ``[M]`` 2026-09-03 ``array_equal`` on 8 of 8 entries — which is what
        licenses the exit instrument's ``135 of 135`` prediction without
        re-running it per row.
        """
        entry = _r4_entry(label)
        select, _ = _coordinate_chart(
            _r4_columns_of(entry), ambient_dim(entry.base),
        )
        points = _r4_sphere_points()
        np.testing.assert_array_equal(
            np.asarray(select(points)),
            np.asarray(entry.orbit_coordinates(points)),
            err_msg=f"{entry.name}: the chart moved — the frame tables ride on it",
        )

    @pytest.mark.foundation
    @pytest.mark.parametrize("label", _R4_ENTRY_LABELS)
    def test_the_projector_is_idempotent_symmetric_and_here_DIAGONAL(
        self, label: str,
    ) -> None:
        """A projector's defining laws, asserted on the SUT's own pair — plus
        the diagonality PREMISE the bit tier rests on."""
        entry = _r4_entry(label)
        select, embed = _coordinate_chart(
            _r4_columns_of(entry), ambient_dim(entry.base),
        )

        def apply(points: np.ndarray) -> np.ndarray:
            return np.asarray(embed(select(points)), dtype=float)

        points = _r4_sphere_points()
        once = apply(points)
        np.testing.assert_array_equal(apply(once), once)  # P^2 = P

        projector = _r4_reynolds_projector(entry.by)
        np.testing.assert_array_equal(projector, projector.T)  # P = P^T
        np.testing.assert_array_equal(
            projector, np.diag(np.diag(projector)),
            err_msg=(
                f"{entry.name}: P_H is no longer diagonal, so the array_equal "
                f"tier above is no longer a theorem — move those rows to "
                f"assert_array_almost_equal_nulp and say why."
            ),
        )

    @pytest.mark.foundation
    @pytest.mark.parametrize("label", _R4_ENTRY_LABELS)
    def test_lift_after_chart_is_the_projector_and_the_round_trip_is_the_identity(
        self, label: str,
    ) -> None:
        r"""``λ ∘ π = P_H`` (the TEETH) and ``π ∘ λ = id`` (declared BLIND).

        ⚠ The second leg is blind to this whole carve: ``[M]`` 2026-09-03 the
        retired hemisphere section satisfies it too, because it writes the
        chart columns back unchanged and its :math:`\sqrt{\cdot}` lands only in
        the DROPPED column.  It ships because the round-trip law is what makes
        ``lift`` a right inverse at all; the discrimination is in the first
        leg, ``[M]`` ``max|section − projector| = 9.943e-01 / 9.735e-01 /
        9.778e-01`` on the three mirror entries.
        """
        entry = _r4_entry(label)
        points = _r4_sphere_points()
        projector = _r4_reynolds_projector(entry.by)

        # TEETH: the lift of the chart is the projection, not a section.
        lifted = np.asarray(entry.lift_coordinates(
            _r4_chart_of(entry, points)), dtype=float)
        np.testing.assert_array_equal(lifted, points @ projector.T)

        # BLIND (labelled): pi . lambda = id on the chart.
        chart = _r4_chart_of(entry, points)
        np.testing.assert_array_equal(
            _r4_chart_of(entry, np.asarray(entry.lift_coordinates(chart))), chart,
        )


def _r4_columns_of(entry: Quotient) -> int | list[int]:
    """The entry's invariant columns — the SUT-side spelling, read from the
    TAG accessors (never from the realization, which the reference owns)."""
    group = entry.by
    if group.rotation_axis is not None:
        return group.rotation_axis
    if group.mirror_axis is not None:
        return [i for i in range(3) if i != group.mirror_axis]
    return list(range(ambient_dim(entry.base)))


# =====================================================================
# B — the barycentre IS the orbit mean (the name's justification)
# =====================================================================
class TestR4TheLiftIsTheOrbitBarycentre:
    r"""The projector is not merely *a* right inverse — it is the orbit's MEAN,
    which is why the name ``orbit_barycentres`` is honest on every entry and
    why the map is equivariant (class C).

    :math:`P_H = \frac{1}{|H|}\sum_{g \in H} g` for a finite group (Reynolds),
    and its continuous form for a compact one.  Both legs below compute that
    average from the group and compare it to the lift.
    """

    @pytest.mark.foundation
    @pytest.mark.parametrize("label", (*_R4_MIRROR_LABELS, "S2/Trivial"))
    def test_a_FINITE_entrys_lift_is_the_mean_over_the_group(
        self, label: str,
    ) -> None:
        """``[M]`` 2026-09-03 ``array_equal`` — bit-exact, because
        ``(x + (−x))/2`` is exactly ``0.0`` and ``(y + y)/2`` is exactly ``y``."""
        entry = _r4_entry(label)
        points = _r4_sphere_points()
        lifted = np.asarray(
            entry.lift_coordinates(_r4_chart_of(entry, points)), dtype=float,
        )
        np.testing.assert_array_equal(
            lifted, _r4_group_mean(entry.by, points),
            err_msg=f"{entry.name}: the lift is not the orbit's mean",
        )

    @pytest.mark.foundation
    @pytest.mark.parametrize("label", _R4_AXIAL_LABELS)
    def test_an_AXIAL_entrys_lift_is_the_orbit_circles_centre(
        self, label: str,
    ) -> None:
        r"""An :math:`O(2)_a` orbit is the circle
        :math:`\{\Omega\cdot\hat e_a = \mu\}`; its centre is
        :math:`\mu\,\hat e_a`, inside the ball.

        ``[M]`` residual against a 16-point trapezoid: ``3.331e-16`` — see
        :func:`_r4_circle_mean` for why a LARGER ``n`` is worse.
        """
        entry = _r4_entry(label)
        axis = entry.by.rotation_axis
        assert axis is not None
        points = _r4_sphere_points()
        lifted = np.asarray(
            entry.lift_coordinates(_r4_chart_of(entry, points)), dtype=float,
        )
        reference = _r4_circle_mean(axis, points)
        np.testing.assert_allclose(
            lifted, reference, rtol=0.0, atol=_R4_ATOL_EXACT,
            err_msg=(
                f"{entry.name}: the lift is not the orbit circle's centre — "
                f"max|D| = {np.abs(lifted - reference).max():.3e}"
            ),
        )


# =====================================================================
# C — equivariance under the normaliser, and ONLY there
# =====================================================================
class TestR4TheLiftIsEquivariantUnderTheNormaliserAndOnlyThere:
    r"""``g P_H = P_H g`` iff :math:`g` normalises :math:`H`, because
    :math:`g P_H g^{-1} = P_{gHg^{-1}}`.

    That identity is the whole justification for replacing the mirror entry's
    hemisphere SECTION (which is not equivariant — an axis-reversing normaliser
    carries it off the fundamental domain) with the orbit BARYCENTRE.

    ⛔ **This class is the ONLY place the property is falsifiable.**  Routing it
    through :meth:`Quotient.induced_action` is Mode-12 blind: the chart drops
    exactly the column the projector zeroes.  ``[M]`` 0 of 9925 kernel answers
    move under R4's lift.  The third row below states that blindness as a
    measurement rather than leaving it to be rediscovered.
    """

    @pytest.mark.foundation
    @pytest.mark.parametrize("label", ("sigma_y", "O2_x"))
    def test_the_projector_commutes_with_every_NORMALISING_motion(
        self, label: str,
    ) -> None:
        """``[M]`` ``max|P g p − g P p| = 0.000e+00`` on the axis-aligned
        normalisers and ``1.218e-16`` on ``C₂ᶻ`` (trig roundoff in the rotation
        matrix), so ``atol=1e-14`` carries a 4-order margin."""
        entry = _r4_entry(label)
        projector = _r4_reynolds_projector(entry.by)
        points = _r4_sphere_points()
        motions = _r4_probe_motions()
        normalising = [(n, m) for n, m in motions if entry.by.is_normalised_by(m)]
        assert len(normalising) >= 5, "the probe set stopped exercising the normaliser"
        for name, motion in normalising:
            linear = np.asarray(motion.linear, dtype=float)
            lhs = (points @ linear.T) @ projector.T
            rhs = (points @ projector.T) @ linear.T
            np.testing.assert_allclose(
                lhs, rhs, rtol=0.0, atol=_R4_ATOL_EXACT,
                err_msg=f"{entry.name}: the lift is not equivariant under {name}",
            )

    @pytest.mark.foundation
    @pytest.mark.parametrize(
        "label,motion_name", [("sigma_y", "C4_z"), ("sigma_y", "R_generic"),
                              ("O2_x", "C4_y")],
    )
    def test_a_NON_normalising_motion_breaks_equivariance_by_ORDER_ONE(
        self, label: str, motion_name: str,
    ) -> None:
        r"""⭐ The ``vv-principles`` #19 discriminating reading.  The positive
        leg above is compatible with a lift that is equivariant under
        EVERYTHING (the identity, say); only this one says the property has
        content.

        ``[M]`` ``max|P g p − g P p| = 9.943e-01`` (``C₄ᶻ`` on
        :math:`\sigma_y`, which conjugates it to :math:`\sigma_x`),
        ``5.524e-01`` (a generic rotation), ``9.943e-01`` (``C₄ʸ`` on
        :math:`O(2)_x`) — so the ``>= 0.1`` floor cannot be met by noise.
        """
        entry = _r4_entry(label)
        motion = dict(_r4_probe_motions())[motion_name]
        assert not entry.by.is_normalised_by(motion), (
            f"{motion_name} normalises {entry.by.name} — this row's premise is "
            f"gone and it has become a false red"
        )
        projector = _r4_reynolds_projector(entry.by)
        points = _r4_sphere_points()
        linear = np.asarray(motion.linear, dtype=float)
        gap = float(np.abs(
            (points @ linear.T) @ projector.T - (points @ projector.T) @ linear.T
        ).max())
        assert gap >= _R4_O1, (
            f"{entry.name} under {motion_name}: the equivariance defect is "
            f"{gap:.3e}, below the O(1) floor — the negative leg has gone blind"
        )

    @pytest.mark.foundation
    def test_the_CHART_is_blind_to_all_of_this_and_the_two_tiers_say_why(
        self,
    ) -> None:
        r"""⛔ The declared blindness, MEASURED, in two tiers.

        ``π(g · P_H p)`` vs ``π(g · p)`` over every representative of the 31
        expressible members that normalises :math:`\sigma_y`, on a shipped fold:

        =========================================  ======  =============  =========
        population                                  rows   ``array_equal``  max|Δ|
        =========================================  ======  =============  =========
        signed-permutation normalisers                100            100  0.000e+00
        the rest (all of them :math:`I_h` elements)     7              0  4.996e-16
        =========================================  ======  =============  =========

        Mechanism: ``is_normalised_by`` admits at ``_ELEMENT_ATOL = 1e-9``, so
        an icosahedral element fixing :math:`\hat e_y` only to ~1e-16 carries
        non-zero off-axis entries and does not commute with :math:`P_H`
        exactly.  The closure's node window is
        ``_NODE_WINDOW_FACTOR × atol = 1e-7`` — NINE orders above — which is
        why ``[M]`` 0 of 9925 kernel answers move.

        Two legs, not one: a single ``allclose`` over the union would hide the
        BIT half, which is the one a future non-signed-permutation motion
        cannot silently weaken.
        """
        fold = SPHERE.quotient(SubgroupOfO3.Mirror("y"))
        projector = _r4_reynolds_projector(fold.by)
        nodes = np.asarray(
            Quadrature.folded_product(n_mu=4, n_phi=8).measure.nodes, dtype=float,
        )
        exact_rows = near_rows = 0
        worst_near = 0.0
        for group in _r4_every_member():
            for motion in group.realization.representatives:
                if not fold.by.is_normalised_by(motion):
                    continue
                linear = np.asarray(motion.linear, dtype=float)
                direct = _r4_chart_of(fold, nodes @ linear.T)
                through = _r4_chart_of(fold, (nodes @ projector.T) @ linear.T)
                if _r4_is_signed_permutation(linear):
                    np.testing.assert_array_equal(direct, through)
                    exact_rows += 1
                else:
                    gap = float(np.abs(direct - through).max())
                    worst_near = max(worst_near, gap)
                    assert gap <= _R4_ATOL_NEAR_BIT, (
                        f"a non-signed-permutation normaliser moved the chart by "
                        f"{gap:.3e} — the kernel's blindness is no longer a theorem"
                    )
                    near_rows += 1
        assert exact_rows >= 100 and near_rows >= 1, (
            f"the probe lost a population: {exact_rows} exact, {near_rows} near "
            f"(measured 100 and 7) — a one-population reading cannot say which "
            f"tier moved"
        )


def _r4_probe_motions() -> list[tuple[str, RigidMotion]]:
    """A named motion set spanning both sides of the normaliser."""
    axis = np.array([0.3, 0.5, 0.81])
    axis = axis / np.linalg.norm(axis)
    return [
        ("sigma_x", RigidMotion.reflection(normal=_R4_I3[0])),
        ("sigma_y", RigidMotion.reflection(normal=_R4_I3[1])),
        ("sigma_z", RigidMotion.reflection(normal=_R4_I3[2])),
        ("C2_z", RigidMotion.rotation_about_axis(axis=_R4_I3[2], angle=np.pi)),
        ("C4_z", RigidMotion.rotation_about_axis(axis=_R4_I3[2], angle=np.pi / 2)),
        ("C4_y", RigidMotion.rotation_about_axis(axis=_R4_I3[1], angle=np.pi / 2)),
        ("-I", RigidMotion.inversion(3)),
        ("R_generic", RigidMotion.rotation_about_axis(axis=axis, angle=0.7)),
    ]


def _r4_every_member() -> list[SubgroupOfO3]:
    """Every expressible member — 31 spellings, 30 distinct groups.  The FINITE
    roster, probed whole (``vv-principles`` #31's finite-roster corollary)."""
    members = [SubgroupOfO3.Trivial, SubgroupOfO3.Dinfh,
               SubgroupOfO3.OctahedralOh, SubgroupOfO3.IcosahedralIh,
               SubgroupOfO3.SO3, SubgroupOfO3.O3]
    members += [SubgroupOfO3.Mirror(a) for a in "xyz"]
    members += [SubgroupOfO3.SO2(a) for a in "xyz"]
    members += [SubgroupOfO3.O2(a) for a in "xyz"]
    members += [SubgroupOfO3.Cn(n) for n in range(1, 9)]
    members += [SubgroupOfO3.Dnh(n) for n in range(1, 9)]
    return members


def _r4_is_signed_permutation(linear: np.ndarray) -> bool:
    """``True`` iff every entry is exactly ``0.0`` or ``±1.0`` — the class in
    which the chart identity is BIT-exact rather than merely close."""
    absolute = np.abs(linear)
    return bool(np.all((absolute == 0.0) | (absolute == 1.0)))


# =====================================================================
# D — an orbit space's dimension is a theorem (finding B2)
# =====================================================================
class TestR4AnOrbitSpacesDimensionIsATheorem:
    r""":math:`\dim(M/H) = \dim M - \dim(H\!\cdot\!p)` at a GENERIC
    :math:`p`, with :math:`\dim(H\!\cdot\!p) = \operatorname{rank}[Xp : X \in
    \mathfrak g]`.

    ⭐ Stated on the ORBIT and not on the group, because ``dim H`` alone is
    wrong on two cases the tree will meet: ``[M]`` :math:`O(3)` on :math:`S^2`
    (rank 2 ⟹ 0 ✓, where ``dim H = 3`` gives −1) and :math:`SO(3)` on
    :math:`\mathbb{R}^3` (rank 2 ⟹ 1 ✓, where ``dim H = 3`` gives 0).  Without
    :meth:`test_dim_H_alone_would_be_WRONG` the law reads as a needlessly
    indirect spelling of ``base.dim − by.dim`` and a later session simplifies it
    into a bug.
    """

    @pytest.mark.foundation
    @pytest.mark.parametrize("label", _R4_ENTRY_LABELS)
    def test_every_constructible_entry_satisfies_the_law(self, label: str) -> None:
        """``[M]`` 8 of 8 — axial ``2−1=1``, mirror and ``S²/{e}`` ``2−0=2``,
        ``R³/{e}`` ``3−0=3``."""
        entry = _r4_entry(label)
        orbit = _r4_orbit_dim(entry.by, _r4_generic_point(entry))
        assert entry.realization.dim == entry.base.dim - orbit, (
            f"{entry.name}: realization.dim={entry.realization.dim} but "
            f"base.dim − dim(orbit) = {entry.base.dim} − {orbit}"
        )

    @pytest.mark.foundation
    @pytest.mark.parametrize("kind", ("axial_on_a_disk", "mirror_on_an_interval"))
    def test_the_two_FORGERIES_are_refused_at_construction_naming_the_law(
        self, kind: str,
    ) -> None:
        r"""§6c: the gate lands with the case it catches.

        ``[M]`` 2026-09-03 on the PRE-carve tree BOTH construct —
        ``S^2/O2_z`` realized on ``Ball(2)`` (the law says ``2−1 = 1 ≠ 2``) and
        ``S^2/sigma_x`` realized on ``[-1,1]`` (``2−0 = 2 ≠ 1``).  Both carry
        ``fundamental_domain=None``, so the pre-existing fd clause returns early
        and provably cannot see them: the dimension law is the only thing that
        can refuse them, which is what makes them ITS witnesses and not a
        second gate's.
        """
        with pytest.raises(ValueError, match=r"generic orbit|dimension"):
            Quotient(**_r4_forged(kind))

    @pytest.mark.foundation
    def test_dim_H_alone_would_be_WRONG_and_the_rank_helper_says_why(self) -> None:
        r"""Reference-side controls ONLY — no catalogue entry exists for either
        (:math:`S^2/O(3)` is #440, :math:`\mathbb{R}^3/SO(3)` is uncatalogued),
        so this row asserts the RANK HELPER and the arithmetic and never builds
        a ``Quotient``.

        ``D_∞h`` on :math:`S^2` is the third row on purpose: there ``dim H``
        HAPPENS to agree, which is what shows the agreement is a coincidence
        rather than a rule.
        """
        on_sphere = np.array(
            [0.31622776601683794, 0.5477225575051661, 0.7745966692414834]
        )
        in_space = np.array([0.3, -0.7, 1.1])

        assert _r4_orbit_dim(SubgroupOfO3.O3, on_sphere) == 2      # 2 − 2 = 0
        assert SubgroupOfO3.O3.dim == 3                            # 2 − 3 = −1 ✗
        assert _r4_orbit_dim(SubgroupOfO3.SO3, in_space) == 2      # 3 − 2 = 1
        assert SubgroupOfO3.SO3.dim == 3                           # 3 − 3 = 0 ✗
        # the coincidence row
        assert _r4_orbit_dim(SubgroupOfO3.Dinfh, on_sphere) == 1
        assert SubgroupOfO3.Dinfh.dim == 1

        # and the GENERIC-point premise: at a POLE the axial rank collapses and
        # the law would refuse all three shipped axial entries.
        pole = np.array([1.0, 0.0, 0.0])
        assert _r4_orbit_dim(SubgroupOfO3.O2("x"), on_sphere) == 1
        assert _r4_orbit_dim(SubgroupOfO3.O2("x"), pole) == 0

    @pytest.mark.foundation
    def test_the_dimension_law_and_the_fundamental_domain_gate_are_DISTINCT(
        self,
    ) -> None:
        r"""⛔ Two clauses in one ``__post_init__``, so each owes a
        DISCRIMINATING input and the ORDER owes a pin (``lessons`` L43c).

        ``[M]`` 2026-09-03 the fd gate's HISTORICAL witness
        (``by=Mirror("y")``, ``realization=[-1,1]``, a hemisphere fd) violates
        BOTH laws after R4, so a row asserting its old message would be pinning
        whichever clause happens to run first.  The two inputs below violate
        exactly one law each; the third asserts, on the both-violating input,
        that the ruled-first fragment is present and the other absent — which
        is what reddens when the clauses swap order.
        """
        # (a) the dimension law ALONE — fd is None, so that clause returns early
        with pytest.raises(ValueError, match=r"generic orbit|dimension") as only_dim:
            Quotient(**_r4_forged("axial_on_a_disk"))
        assert "fundamental domain" not in str(only_dim.value)

        # (b) the fd gate ALONE — the law reads 2 − 0 = 2 == Ball(2).dim ✓, and
        #     the half-meridian's ANTIPODAL PAIR drops it to dim 1.
        select, embed = _coordinate_chart([0, 2], ambient_dim(SPHERE))
        half_meridian = FundamentalDomain(
            SPHERE, ((0.0, 1.0, 0.0), (0.0, -1.0, 0.0)), "y=0",
        )
        assert half_meridian.dim == 1 and Ball(2).dim == 2
        with pytest.raises(ValueError, match="must describe the same") as only_fd:
            Quotient(
                base=SPHERE, by=SubgroupOfO3.Mirror("y"), realization=Ball(2),
                orbit_coordinates=select, lift_coordinates=embed,
                lift_codomain=Ball(3), fundamental_domain=half_meridian,
            )
        assert "generic orbit" not in str(only_fd.value)

        # (c) the ORDER, on the input that violates both
        hemisphere = FundamentalDomain(SPHERE, ((0.0, 1.0, 0.0),), "y>=0")
        select1, embed1 = _coordinate_chart(1, ambient_dim(SPHERE))
        with pytest.raises(ValueError) as both:
            Quotient(
                base=SPHERE, by=SubgroupOfO3.Mirror("y"),
                realization=COSINE_INTERVAL,
                orbit_coordinates=select1, lift_coordinates=embed1,
                lift_codomain=Ball(3), fundamental_domain=hemisphere,
            )
        message = str(both.value)
        # ⚠ RULED ORDER: the dimension law runs FIRST (open ruling 1 of the
        # plan).  Flip the two clauses in ``__post_init__`` and this reds.
        assert "generic orbit" in message or "dimension" in message
        assert "must describe the same" not in message


# =====================================================================
# E — orbit_barycentres: ONE concept on both widths (finding B6)
# =====================================================================
class TestR4OrbitBarycentresIsOneConceptOnBothWidths:
    r"""``ambient_representatives`` promised a REPRESENTATIVE and returned a
    barycentre — a fold's nodes passed through (representatives ON
    :math:`S^2`) while a marginal's were lifted (barycentres inside the ball).
    Two things under one name.  ``orbit_barycentres`` is one:

    * chart width → ``lift_coordinates(points)``
    * ambient width → ``lift_coordinates(orbit_coordinates(points))`` = ``P_H p``
    """

    @pytest.mark.foundation
    @pytest.mark.parametrize("label", _R4_ENTRY_LABELS)
    def test_both_widths_answer_the_projector(self, label: str) -> None:
        entry = _r4_entry(label)
        points = _r4_sphere_points()
        projector = _r4_reynolds_projector(entry.by)
        expected = points @ projector.T

        np.testing.assert_array_equal(
            np.asarray(entry.orbit_barycentres(points), dtype=float), expected,
            err_msg=f"{entry.name}: the AMBIENT-width arm is not the projector",
        )
        np.testing.assert_array_equal(
            np.asarray(
                entry.orbit_barycentres(_r4_chart_of(entry, points)), dtype=float,
            ),
            expected,
            err_msg=f"{entry.name}: the CHART-width arm is not the projector",
        )

    @pytest.mark.foundation
    def test_the_ambient_arm_no_longer_passes_a_FOLDS_NODES_through(self) -> None:
        r"""⭐ The rename's whole content, and the only row that sees it.

        ``[M]`` 2026-09-03 on the PRE-carve tree
        ``ambient_representatives(fold_nodes)`` is ``array_equal`` to the nodes
        (a documented pass-through) and those nodes are ON :math:`S^2`; after
        R4 the answer is :math:`P_H p`, OFF :math:`S^2` and inside the ball.
        Without this row the rename is a pure re-labelling with no witness.
        """
        fold = SPHERE.quotient(SubgroupOfO3.Mirror("y"))
        nodes = np.asarray(
            Quadrature.folded_product(n_mu=4, n_phi=8).measure.nodes, dtype=float,
        )
        assert SPHERE.contains(nodes)          # the representatives ARE on S^2
        answer = np.asarray(fold.orbit_barycentres(nodes), dtype=float)
        assert not np.array_equal(answer, nodes), (
            "the ambient arm still passes a fold's nodes through — "
            "`orbit_barycentres` is two concepts under one name again"
        )
        np.testing.assert_array_equal(
            answer, nodes @ _r4_reynolds_projector(fold.by).T,
        )
        assert not SPHERE.contains(answer)
        assert Ball(3).contains(answer)

    @pytest.mark.foundation
    @pytest.mark.parametrize("shape", _R4_FOLDS)
    def test_it_is_idempotent_and_INJECTIVE_on_a_folds_node_set(
        self, shape: tuple[int, int],
    ) -> None:
        r"""The projection's one silent failure mode is collapsing two orbits.

        ``[D]`` it cannot: on a :math:`\sigma_a`-fold two representatives
        sharing :math:`(x_b, x_c)` share :math:`|x_a| = \sqrt{1-x_b^2-x_c^2}`
        and both have :math:`x_a \ge 0`, so they are the same node.  ``[M]``
        min chart separation ``1.155e+00 / 4.403e-01 / 2.751e-01 / 1.510e-01``
        on ``folded_product`` (2,4)/(4,6)/(4,8)/(8,8) — six orders above the
        ``1e-7`` node window.
        """
        n_mu, n_phi = shape
        measure = Quadrature.folded_product(n_mu=n_mu, n_phi=n_phi).measure
        entry = measure.support
        assert isinstance(entry, Quotient)
        nodes = np.asarray(measure.nodes, dtype=float)

        once = np.asarray(entry.orbit_barycentres(nodes), dtype=float)
        np.testing.assert_array_equal(
            np.asarray(entry.orbit_barycentres(once), dtype=float), once,
        )
        separation = np.linalg.norm(once[:, None, :] - once[None, :, :], axis=-1)
        np.fill_diagonal(separation, np.inf)
        assert float(separation.min()) > 1e-7, (
            f"folded_product({n_mu},{n_phi}): the projection collapsed two "
            f"orbits — min separation {float(separation.min()):.3e}"
        )

    @pytest.mark.foundation
    def test_a_width_in_neither_system_is_refused_BY_NAME(self) -> None:
        """Re-keys the retired ``test_a_width_in_neither_system_is_refused_by_name``:
        the method's name changed, the refusal did not."""
        fold = SPHERE.quotient(SubgroupOfO3.Mirror("y"))
        with pytest.raises(ValueError, match="neither"):
            fold.orbit_barycentres(np.zeros((4, 7)))

    @pytest.mark.foundation
    def test_the_HARNESS_fold_answer_now_agrees_with_the_kernel(self) -> None:
        r"""⭐ ``tests/_harness/references.py``'s partner map is derived from
        ``_embedded_nodes``, so R4 moves ONE of its 33 (rule × axis) answers per
        fold — and it moves it INTO agreement with the kernel.

        ``[M]`` 2026-09-03: 31 of 33 answers unchanged; the 2 that move are
        ``axis="y"`` on the two folds, where the harness RAISED (residual
        ``1.19e+00``) and now returns the identity permutation (residual
        ``0.00e+00``) — which is exactly what ``permutation_under`` has
        answered since #429 tracker 2.2b.  ``[M]`` **no call site passes
        ``axis="y"`` on a fold** (both pass ``"x"``: a literal at
        ``test_snmesh_realizer_wiring.py:439,450``, and the CYL fold's only
        realized face carries ``ReflectiveBoundary(axis='x')`` at
        ``test_b3_domain_narrowing.py:164,175``), so no committed assertion
        moved — which is why this capability needs a gate of its own or it
        lands unrecorded.
        """
        from tests._harness.references import mirror_partner_indices

        quad = Quadrature.folded_product(n_mu=4, n_phi=8)
        reference = mirror_partner_indices(quad, "y")
        np.testing.assert_array_equal(reference, np.arange(quad.N))

        kernel = quad.measure.permutation_under(RigidMotion.reflection(normal=_R4_I3[1]), atol=1e-9)
        assert kernel is not None
        np.testing.assert_array_equal(np.asarray(kernel.indices), reference)

        # the axis the two folded CALL SITES actually pass is untouched
        np.testing.assert_array_equal(
            mirror_partner_indices(quad, "x"),
            mirror_partner_indices(quad, "x"),
        )


# =====================================================================
# F — the lift is a FIELD, not a tag branch (finding B1)
# =====================================================================
class TestR4TheLiftIsADerivationOutputNotATagBranch:
    r"""Every other derivation output is a field; the lift was a three-arm
    branch on the group's family, with a ``NotImplementedError`` telling the
    reader to add a fourth arm.  R4 makes it two REQUIRED fields, so a new
    entry cannot forget it and no consumer re-derives the family from the tag.
    """

    @pytest.mark.foundation
    def test_the_two_fields_are_REQUIRED_and_carry_no_default(self) -> None:
        r"""⚠ ``vv`` Mode-8, signature-tautological class: "pass the fields and
        it constructs" is green before and after.  The falsifiable half is the
        OMISSION, plus a ``dataclasses.fields`` read — and the second half is a
        grep-tier obligation wearing a test, labelled as such
        (``lessons`` L59d).
        """
        names = {f.name: f for f in dataclasses.fields(Quotient)}
        for name in ("lift_coordinates", "lift_codomain"):
            assert name in names, f"{name} is not a field of Quotient"
            field = names[name]
            assert field.default is dataclasses.MISSING, (
                f"{name} has a default, so a new entry CAN forget it — which is "
                f"the whole reason the plan made it required"
            )
            assert field.default_factory is dataclasses.MISSING

        select, _ = _coordinate_chart(0, ambient_dim(SPHERE))
        with pytest.raises(TypeError):
            Quotient(  # type: ignore[call-arg]
                base=SPHERE, by=SubgroupOfO3.O2("x"),
                realization=COSINE_INTERVAL, orbit_coordinates=select,
            )

    @pytest.mark.foundation
    @pytest.mark.parametrize("label", _R4_ENTRY_LABELS)
    def test_lift_assembles_the_typed_arrow_from_the_two_fields(
        self, label: str,
    ) -> None:
        """``lift`` is derived from the fields and from nothing else — no
        family test survives inside it."""
        entry = _r4_entry(label)
        arrow = entry.lift
        assert arrow.domain is entry
        assert arrow.codomain == entry.lift_codomain
        chart = _r4_chart_of(entry, _r4_sphere_points())
        np.testing.assert_array_equal(
            np.asarray(arrow(chart), dtype=float),
            np.asarray(entry.lift_coordinates(chart), dtype=float),
        )

    @pytest.mark.foundation
    def test_an_entry_OUTSIDE_the_three_retired_arms_answers(self) -> None:
        r"""§6c, and R4's cleanest witness.

        ``[M]`` 2026-09-03 on the PRE-carve tree a hand-built
        ``S^2/O_h`` CONSTRUCTS (``O_h`` is its own ``orbit_stabiliser``; the
        dimension law reads ``2 − 0 = 2`` ✓) and ``entry.lift`` raises
        ``NotImplementedError: no lift is spelled for S^2/Oh``.  After R4 it
        answers from its own field.  The input the branch could not serve
        exists in the tree at the moment the branch retires — so the
        retirement is a capability, not merely "no test broke".
        """
        entry = Quotient(**_r4_hand_built_octahedral())
        chart = _r4_disk_points()
        lifted = np.asarray(entry.lift(chart), dtype=float)
        assert lifted.shape == (chart.shape[0], 3)
        np.testing.assert_array_equal(
            lifted, np.asarray(entry.lift_coordinates(chart), dtype=float),
        )
        assert entry.lift.codomain == entry.lift_codomain

    @pytest.mark.foundation
    @pytest.mark.parametrize("label", _R4_ENTRY_LABELS)
    def test_the_codomain_is_the_BALL_for_both_sphere_families(
        self, label: str,
    ) -> None:
        r"""⛔ **This row REPEALS a claim, it does not extend one.**  The
        superseded gate's argument was *"a gate asserting one uniform codomain
        across the three families would be FALSE"* — ``[M]`` pre-carve the
        codomains are axial ``D^3``, mirror ``S^2``, trivial ``S^2``.  After R4
        BOTH sphere families land in the ball, because both lifts are the orbit
        MEAN and a mean is a point of the convex hull, not of the sphere.  Only
        the trivial entry's lift is a section, and its codomain is its base.

        The ERR-080 discriminator survives and is what keeps the row loaded:
        the image is off :math:`S^2` except where the orbit is a single point.
        """
        entry = _r4_entry(label)
        if entry.by.is_trivial:
            assert entry.lift_codomain == entry.base
            return
        assert entry.lift_codomain == Ball(3)
        chart = _r4_chart_of(entry, _r4_sphere_points())
        image = np.asarray(entry.lift_coordinates(chart), dtype=float)
        assert Ball(3).contains(image)
        assert not SPHERE.contains(image)

    @pytest.mark.foundation
    def test_barycentre_is_defined_on_EVERY_entry_and_refuses_only_a_BARE_manifold(
        self,
    ) -> None:
        r"""Item 4 of the carve: ``barycentre(entry)`` IS ``entry.lift`` for
        every catalogued entry, no longer axial-only.

        ``[M]`` 2026-09-03 on the PRE-carve tree ``barycentre`` refuses the
        mirror AND the trivial entries by name (*"...by an axial group; got
        'S^2/Trivial'"*), which is what makes this row red before the carve.
        The refusal leg survives, narrowed to what it was always ABOUT — an
        argument that is not an orbit space at all (``vv-principles`` #11: the
        positive leg is what turns "it raises" into "the claim is right").
        """
        for label in _R4_ENTRY_LABELS:
            entry = _r4_entry(label)
            arrow = barycentre(entry)
            assert arrow.domain is entry
            assert arrow.codomain == entry.lift_codomain
            chart = _r4_chart_of(entry, _r4_sphere_points())
            np.testing.assert_array_equal(
                np.asarray(arrow(chart), dtype=float),
                np.asarray(entry.lift_coordinates(chart), dtype=float),
            )
        with pytest.raises(ValueError):
            barycentre(COSINE_INTERVAL)  # type: ignore[arg-type]


# =====================================================================
# G — the hemisphere section retires, and its refusal keeps a home
# =====================================================================
class TestR4TheHemisphereSectionIsRetiredAndItsRefusalHasAHome:
    r"""``_hemisphere_section`` was the mirror entry's lift: a genuine section
    onto the fundamental domain, at the cost of a :math:`\sqrt{\cdot}`, a
    ``rho^2 > 1 + 1e-12`` refusal, and NON-equivariance under the
    axis-reversing normalisers (class C).  The projector replaces all three.

    ⛔ **Retirement-means-test-migration, run with its denominator.**  ``[M]``
    2026-09-03 a Python-``re`` sweep of ``tests/`` for the section's three
    distinctive message fragments (``"closed unit disk"``, ``"rho^2"``,
    ``"hemisphere section is defined"``) returns **0 hits**, with the positive
    control finding the production lines.  So the refusal had NO witness
    anywhere in the tree and nothing migrates — which is exactly why its
    successor's teeth are NET-NEW (``lessons`` L4) and why the second row below
    is written rather than assumed.
    """

    @pytest.mark.foundation
    def test_the_symbol_is_gone(self) -> None:
        from orpheus.numerics import manifold

        assert not hasattr(manifold, "_hemisphere_section"), (
            "the hemisphere section still ships — its sqrt, its literal and "
            "its refusal were the three things R4 retires"
        )

    @pytest.mark.foundation
    def test_an_OFF_DISK_chart_point_is_refused_by_the_MEMBERSHIP_predicate(
        self,
    ) -> None:
        r"""The refusal's new home.  ``embed`` is a linear map and validates
        nothing by design; the predicate that a chart point is a point of the
        orbit space is :meth:`Quotient.contains`, and it already answers.

        ``[M]`` 2026-09-03: ``entry.contains([[0.9, 0.9]])`` is ``False``
        (:math:`\rho^2 = 1.62`, a pair corresponding to NO direction) and
        ``entry.contains([[1.0, 0.0]])`` is ``True`` — the boundary circle,
        which the cylindrical march seeds at :math:`x_a = 0` exactly, so a
        strict inequality here would refuse a direction production marches
        from.
        """
        entry = SPHERE.quotient(SubgroupOfO3.Mirror("y"))
        assert not entry.contains(np.array([[0.9, 0.9]]))
        assert entry.contains(np.array([[1.0, 0.0]]))       # the mirror circle
        assert entry.contains(_r4_disk_points())            # the positive leg


# =====================================================================
# H — (M/H)/{e} IS M/H (finding B3)
# =====================================================================
class TestR4TheTrivialQuotientOfAnOrbitSpaceIsThatOrbitSpace:
    r""":math:`(M/H)/\{e\} = M/H` is the same theorem as :math:`M/\{e\} = M`,
    and the door already names the trivial group as its ONE admitted exception.
    What it did instead was build a SECOND ``Quotient`` — ``[M]`` pre-carve
    ``S^2/sigma_y/Trivial``, unequal to and not identical with its own base:
    the very class (one orbit space, two objects) that ``#432``'s naming law
    and ERR-080 exist to make unspellable.
    """

    @pytest.mark.foundation
    @pytest.mark.parametrize("label", ("sigma_y", "O2_x"))
    def test_M_over_H_mod_e_IS_M_over_H_by_IDENTITY(self, label: str) -> None:
        """Identity, not equality: the catalogue memoises, so the one object is
        the assertion.  ``[M]`` pre-carve ``is`` is False and ``==`` is False."""
        entry = _r4_entry(label)
        assert entry.quotient(SubgroupOfO3.Trivial) is entry

    @pytest.mark.foundation
    def test_the_measure_level_verb_keeps_every_node_on_the_ENTRYS_OWN_support(
        self,
    ) -> None:
        """``[M]`` pre-carve the 16 nodes land on ``S^2/sigma_y/Trivial``; after
        R4 they land on ``S^2/sigma_y``, which is where they already were."""
        measure = Quadrature.folded_product(n_mu=4, n_phi=8).measure
        folded = measure.quotient(SubgroupOfO3.Trivial)
        assert folded.support is measure.support
        assert folded.support.name == "S^2/sigma_y"
        np.testing.assert_array_equal(folded.nodes, measure.nodes)
        np.testing.assert_array_equal(folded.weights, measure.weights)

    @pytest.mark.foundation
    @pytest.mark.parametrize("label", ("S2", "CIRCLE", "COSINE_INTERVAL"))
    def test_M_mod_e_on_a_NON_quotient_base_still_DERIVES(self, label: str) -> None:
        r"""The control the short-circuit must not swallow.  ``M/{e}`` on an
        ordinary manifold is still the derived identity ENTRY — the catalogue's
        positive control on its own machinery — and the new arm fires only when
        the base is already an orbit space."""
        base = {"S2": SPHERE, "CIRCLE": CIRCLE,
                "COSINE_INTERVAL": COSINE_INTERVAL}[label]
        entry = base.quotient(SubgroupOfO3.Trivial)
        assert isinstance(entry, Quotient)
        assert entry is not base
        assert entry.realization == base
        assert entry.dim == base.dim
        assert entry.singular_stratum is None and entry.is_free


# =====================================================================
# I — the trivial group is asked STRUCTURALLY (finding B5)
# =====================================================================
class TestR4TheTrivialGroupIsAskedOfTheTypeNotOfItsName:
    r"""``name == "Trivial"`` appeared at five sites (3 in ``manifold.py``, 2 in
    ``basis/descent.py``) — a stringly-typed dispatch on a property the type can
    answer (``coding-elegance`` anti-pattern #4).

    ⚠ The roster is FINITE and enumerable, so it is probed WHOLE rather than
    sampled: ``vv-principles`` #13's ladder rule governs unbounded families,
    and #31's finite-roster corollary governs this one.  ONE row with the
    roster inside, not 31 parametrized rows (#20, row inflation).
    """

    @pytest.mark.foundation
    def test_is_trivial_over_EVERY_expressible_member(self) -> None:
        """``[M]`` 31 spellings, 30 distinct members; ``is_trivial`` is True on
        exactly the 2 spellings of ``{e}`` (``Trivial`` and ``Cn(1)``, which R1
        normalises onto the same tag at construction)."""
        members = _r4_every_member()
        assert len(members) == 31, "the roster changed — re-derive the count"
        trivial = [g for g in members if g.is_trivial]
        assert [g.name for g in trivial] == ["Trivial", "Trivial"]
        assert len(set(trivial)) == 1
        for group in members:
            assert group.is_trivial == (group == SubgroupOfO3.Trivial), group.name

    @pytest.mark.foundation
    def test_it_agrees_with_a_STRUCTURAL_reference_off_the_realization(
        self,
    ) -> None:
        r"""The independent side: :math:`\{e\}` is the group whose identity
        component is 0-dimensional AND which has exactly one coset — read off
        the realization, which the property need not consult.  ``[M]`` the two
        answers agree on 31 of 31."""
        for group in _r4_every_member():
            realization = group.realization
            structural = (realization.component.dim == 0
                          and len(realization.representatives) == 1)
            assert group.is_trivial is structural, group.name

    @pytest.mark.foundation
    def test_the_five_call_sites_answer_identically(self) -> None:
        """Not a refactor tautology: each site's OBSERVABLE is asserted, with
        the values measured on the pre-carve tree."""
        from orpheus.numerics.basis.descent import Descent
        from orpheus.numerics.basis.spherical_harmonic_basis import (
            SphericalHarmonicBasis,
        )

        entry = SPHERE.quotient(SubgroupOfO3.Trivial)
        descent = Descent.for_entry(entry, 3)
        assert isinstance(descent.frame_basis, SphericalHarmonicBasis)   # :127
        assert isinstance(descent.upstairs, SphericalHarmonicBasis)      # :149

        arrow = quotient_onto(entry, SPHERE)                             # :1372
        assert arrow is not None
        points = _r4_sphere_points()
        np.testing.assert_array_equal(np.asarray(arrow(points)), points)

        # :1510 — the door admits the trivial group on every base, including
        # one that is already an orbit space (class H).
        assert isinstance(CIRCLE.quotient(SubgroupOfO3.Trivial), Quotient)
        # ⚠ item 4 (barycentre widened to every entry) is a DIFFERENT subject
        # and lives in its own row — a mixed-subject assertion here would make
        # its red unattributable.  See
        # ``TestR4TheLiftIsADerivationOutputNotATagBranch
        # ::test_barycentre_is_defined_on_EVERY_entry_and_refuses_only_a_BARE_manifold``.


class TestR4TheLiftCodomainIsComparedAndGated:
    r"""Elegance review of R4 (2026-09-03), finding 1.

    `[M]` with ``lift_codomain`` excluded from ``__eq__``,
    ``replace(entry, lift_codomain=SPHERE)`` compared EQUAL to the catalogue
    entry and — ``barycentre`` being memoised on the entry — whichever was
    asked first answered for both: an arrow off the sphere declaring itself
    onto it, ERR-080's shape, re-minted by the field built to refuse it.  Two
    legs: the field is compared, and a codomain of the wrong width is refused
    where the entry is written.
    """

    @pytest.mark.foundation
    def test_a_different_codomain_is_a_different_entry(self) -> None:
        entry = SPHERE.quotient(SubgroupOfO3.O2("z"))
        forged = dataclasses.replace(entry, lift_codomain=SPHERE)
        assert forged != entry
        assert hash(forged) != hash(entry)
        # and the memo cannot be poisoned through the forgery
        assert barycentre(entry).codomain == Ball(3)
        assert barycentre(forged).codomain == SPHERE

    @pytest.mark.foundation
    def test_a_codomain_of_the_wrong_width_is_refused(self) -> None:
        entry = SPHERE.quotient(SubgroupOfO3.Mirror("y"))
        with pytest.raises(ValueError, match="ambient space"):
            dataclasses.replace(entry, lift_codomain=Ball(2))


class TestR4TheDimensionLawIsAMaximumOverTheProbeSet:
    r"""Elegance review of R4 (2026-09-03), finding 2.

    Orbit dimension is upper semicontinuous — it drops only on the singular
    stratum — so the generic value is the MAXIMUM over a probe set.  `[M]` a
    single-point spelling evaluated at a probe row placed ON the axis refused
    :math:`S^2/O(2)_z` and ADMITTED the disk forgery: the law inverted.  This
    row poisons one probe row the same way and asserts the law still
    discriminates; the positive control is the shipped probe.
    """

    @pytest.mark.foundation
    def test_a_probe_row_on_the_axis_does_not_invert_the_law(self, monkeypatch) -> None:
        from orpheus.numerics import manifold as mod

        poisoned = np.array(mod._GENERIC_SPHERE_PROBE, copy=True)
        poisoned[0] = np.array([0.0, 0.0, 1.0])  # ON the z-axis: an O(2)_z orbit of dimension 0
        monkeypatch.setattr(mod, "_GENERIC_SPHERE_PROBE", poisoned)
        mod._catalogued_quotient.cache_clear()
        try:
            entry = SPHERE.quotient(SubgroupOfO3.O2("z"))  # must still construct
            assert entry.realization.dim == 1
            chart, lift = _coordinate_chart([0, 1], 3)
            with pytest.raises(ValueError, match="realizes its dimension"):
                Quotient(
                    base=SPHERE, by=SubgroupOfO3.O2("z"), realization=Ball(2),
                    orbit_coordinates=chart, lift_coordinates=lift, lift_codomain=Ball(3),
                )
        finally:
            mod._catalogued_quotient.cache_clear()

    @pytest.mark.foundation
    def test_the_group_answers_the_orbit_dimension_and_fixes_is_its_zero(self) -> None:
        probe = np.asarray(_r4_sphere_points(), dtype=float)
        assert SubgroupOfO3.O2("x").generic_orbit_dimension(probe) == 1
        assert SubgroupOfO3.SO3.generic_orbit_dimension(probe) == 2
        assert SubgroupOfO3.Dnh(4).generic_orbit_dimension(probe) == 0
        on_axis = np.array([[1.0, 0.0, 0.0], [-0.3, 0.0, 0.0]])
        component = SubgroupOfO3.O2("x").realization.component
        assert component.fixes(on_axis, atol=1e-12)
        assert component.orbit_dimension(on_axis) == 0


# ===== the dimension law's SHAPE witness (test-architect, 2026-09-03; the battery
# arm `dim_law_reads_dim_H` reddened 0 rows without it) =====


@pytest.mark.foundation
@pytest.mark.parametrize("case", ("S2_mod_O3_is_a_POINT", "R3_mod_O3_is_a_RAY"))
def test_the_law_subtracts_the_ORBITS_dimension_not_dim_H(case: str) -> None:
    r"""⭐ The only inputs on which ``dim M − dim(generic orbit)`` and
    ``dim M − dim H`` DISAGREE, and both construct.

    ``[M]`` 2026-09-03: shipped-entry-wide the two formulas agree (axial
    ``1 = 1``, mirror/trivial ``0 = 0``), so a battery arm replacing the rank
    with ``group.dim`` reddens **0 of 4597** rows.  These two do the work:

    * :math:`S^2/O(3)` is a POINT — one orbit, ``2 − 2 = 0`` ✓, while
      ``dim H = 3`` gives ``−1``;
    * :math:`\mathbb{R}^3/O(3)` is the radial RAY — ``3 − 2 = 1`` ✓, while
      ``dim H = 3`` gives ``0``.

    Neither is catalogued (:math:`S^2/O(3)` is #440), so the row builds the
    ``Quotient`` directly — which is exactly the tier the law lives at.
    """
    if case == "S2_mod_O3_is_a_POINT":
        base, realization, codomain = SPHERE, IndexSet(label="pt", n=1), Ball(3)
    else:
        base, realization, codomain = RealSpace(3), HALF_LINE, RealSpace(3)
    entry = Quotient(
        base=base, by=SubgroupOfO3.O3, realization=realization,
        orbit_coordinates=_all_coordinates, lift_coordinates=_all_coordinates,
        lift_codomain=codomain,
    )
    assert entry.dim == base.dim - 2, (
        f"{entry.name}: a generic O(3)-orbit is 2-dimensional on both bases, "
        f"so the orbit space is {base.dim - 2}-dimensional"
    )
    assert SubgroupOfO3.O3.dim == 3, (
        "dim H = 3 here, so a law written as dim M − dim H would refuse this "
        "entry — which is the whole reason the law is stated on the ORBIT"
    )
