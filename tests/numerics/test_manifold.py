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
"""

from __future__ import annotations

import dataclasses

import numpy as np
import pytest

from orpheus.numerics.manifold import (
    CIRCLE,
    COSINE_INTERVAL,
    ENERGY,
    HALF_LINE,
    REAL_LINE,
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
        return SPHERE.quotient(SubgroupOfO3.SO2)

    def test_the_realization_is_the_cosine_interval(self, s2_mod_so2):
        r""":math:`S^2/SO(2) \cong [-1,1]`, with :math:`\mu = \hat\Omega \cdot \hat z`."""
        assert s2_mod_so2.realization == COSINE_INTERVAL
        assert s2_mod_so2.name == "S^2/SO2"

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
        assert s2_mod_so2.singular_stratum == (-1.0, 1.0)
        assert not s2_mod_so2.is_free

    def test_the_stratum_lies_on_the_realization(self, s2_mod_so2):
        """A singular point is still a point of the orbit space."""
        assert s2_mod_so2.contains(np.asarray(s2_mod_so2.singular_stratum))

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
        assert q.singular_stratum == ()
        assert q.is_free
        assert q.dim == base.dim

    def test_the_trivial_answer_agrees_with_the_shipped_string_twin(self):
        r"""⭐ The mint is a RE-TYPING, and this pins the migration.

        ``AngularSymmetry.support`` (``quadrature/registry.py:869``) is
        already an orbit-space catalogue in the string vocabulary: it
        maps ``SO2 -> "[-1,1]"`` and ``Trivial -> "S^2"`` and raises for
        anything else, with the same "extend it when a geometry first
        spends it" shape as :meth:`Manifold.quotient`'s refusal.  It is
        the Pattern-2 twin this type exists to absorb, so the two must
        agree on every row both can answer — otherwise the eventual
        collapse (tracker 2.2) silently changes an answer.

        `[M]` both rows agree today.  ⚠ If this ever reddens, do NOT
        re-baseline: one of the two is wrong about a quotient, and which
        one is a mathematical question, not a test-maintenance one.
        """
        from orpheus.numerics.quadrature.registry import AngularSymmetry

        rows = {
            SubgroupOfO3.SO2: COSINE_INTERVAL,
            SubgroupOfO3.Trivial: SPHERE,
        }
        for group, expected in rows.items():
            mine = SPHERE.quotient(group).realization
            theirs = AngularSymmetry(
                continuous_isotropy=group,
                discrete_residual=SubgroupOfO3.Trivial,
            ).support
            assert mine == expected, f"{group.name}: mine gave {mine.name}"
            assert mine.name == theirs, (
                f"the mint disagrees with the shipped string twin for "
                f"S^2/{group.name}: Manifold says {mine.name!r}, "
                f"AngularSymmetry.support says {theirs!r}. One of them is "
                f"wrong about the quotient -- do not re-baseline this."
            )

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
        assert "Sphere/SO2" in msg
        assert "manifold CLASS" in msg


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
    """``_ambient`` decides how many columns ``contains`` consumes.

    Its ``match`` is deliberately exhaustive with a raising fallthrough,
    so a new member that forgets it fails loudly rather than silently
    mis-splitting a product's coordinates.
    """
    from orpheus.numerics.manifold import _ambient

    for m in _ALL:
        assert _ambient(m) >= 1
    assert _ambient(SPHERE.quotient(SubgroupOfO3.SO2)) == 1
    assert _ambient(SPHERE * CIRCLE) == 5
