r"""S5 of the numerics untangling: invariance is a composition — the gates.

Design: ``scratch/numerics_untangling/packages/S5_design.md``; rulings 42–44 of
``.claude/plans/numerics_untangling.md``; the spec, with every row's mutation and its
first red, is ``scratch/numerics_untangling/packages/S5_gate_spec.md``. Written by the
test-architect at ``c27830f5``, before any production change.

The claim under test: every invariance verb on a :class:`DiscreteMeasure` is the
composition

.. math::

   \pi_g \;=\; \mu.\mathrm{relabelling\_to}\big(\mu.\mathrm{pushforward}(
   \mu.\mathrm{support}.\mathrm{action}(g))\big),

the support carries the action (``GSpace``), a support with no action is a
``TypeError``, and ``orpheus.numerics.invariance`` is gone.

**Rows that need the S5 API.** No row imports a post-S5 name at module scope, so
collection never dies: a new name is fetched inside the row by :func:`_api`, which
FAILS the row with ``S5 API absent: <name>`` when the name does not exist. On today's
tree those rows are red for that reason; the spec's §6 lists each row's first red. Rows that read only
names stable across the carve (the measure verbs, ``Quadrature``, the manifolds, the
baseline JSON) run green today and must stay green.

**The baseline** (``data/s5_invariance_baseline.json``, generator beside it) is a
RECORD of the old kernel at ``c27830f5``: 104 rules, 3242 verdicts, 2306 finite
(rule, group) rows of which 914 carry a certificate, 42 product-factor rows. The
REFERENCE beside it (``TestTheCertificateAgainstAnIndependentMatcher``) builds every
permutation from the elements' raw matrices and ``scipy.optimize.linear_sum_assignment``,
so the record is anchored to something the old kernel did not compute.

**Markers.** Every row is a software invariant with no theory ``:label:``:
``foundation``, no ``verifies``. ``rests_on`` names the existing gates whose
green each family assumes.

**Discipline.** No production call in a ``parametrize`` list (rows are keyed by LABEL
and built in the body, so a ``-p`` battery reaches them); every ``pytest.raises``
carries a ``match=`` fragment keyed to the triggering argument; every spy asserts its
own installation; group answers are compared as SETS.
"""

from __future__ import annotations

import functools
import importlib
import json
import pathlib
import re

import numpy as np
import pytest
from scipy.optimize import linear_sum_assignment

from orpheus.geometry.deck import mirror_motion, translation_motion
from orpheus.geometry.transformation import RigidMotion
from orpheus.geometry.manifold import (
    CIRCLE,
    COSINE_INTERVAL,
    HALF_LINE,
    REAL_LINE,
    SPHERE,
    UNIT_INTERVAL,
    Ball,
    Interval,
    Quotient,
    RealSpace,
)
from orpheus.numerics.index_set import ENERGY, EnergyGroups, IndexSet
from orpheus.numerics.measure import DiscreteMeasure, equispaced
from orpheus.numerics.quadrature import Quadrature, quadrature_registry
from orpheus.numerics.quadrature.rules_1d import (
    gauss_legendre_on_mu,
    gauss_legendre_on_polar_orbit,
)
from orpheus.numerics.quadrature.rules_circle import NODE_ALIGNED, STAGGERED, periodic_trapezoid
from orpheus.geometry.symmetry import SubgroupOfO3

pytestmark = pytest.mark.foundation

DATA = pathlib.Path(__file__).parent / "data" / "s5_invariance_baseline.json"
BASELINE = json.loads(DATA.read_text())
RULES: dict = BASELINE["rules"]
FACTORS: dict = BASELINE["product_factors"]

#: The refusal fragments S5 owes (spec §3). Disjoint by construction, asserted once.
FRAG_NOT_GSPACE = "is not a G-space"
FRAG_NO_O3 = "no O(3) action"
FRAG_DOES_NOT_MAP = "does not map"
#: A motion of the wrong dimension. NOT "dimension": numpy's own shape error says
#: "core dimension", so a missing refusal would still match it (battery arm
#: ``no_dimension_check``, `[M]` the 3-D-motion-on-an-interval row stayed green).
FRAG_WRONG_DIMENSION = "does not act"
#: Pinned since #429 tracker 2.2b by test_manifold.py; the rename keeps it.
FRAG_NOT_NORMALISE = "normalise"
FRAG_SUPPORT = "support"
#: Ruling 47 / F6: ``spherical_product``'s own refusal of a polar factor off the chart.
FRAG_POLAR_SUPPORT = "polar factor must be a rule on"
#: The pushforward's refusal it pre-empts (pinned by test_rules_product.py today).
FRAG_PUSHFORWARD_DOMAIN = "map's domain must be the measure's support"

_MANIFOLD_ROW = (
    "tests/gates/geometry/test_manifold.py::"
    "TestTheInducedActionIsWellDefinedAndActsAsTheTheoremSays::"
)
_MATCHER_ROWS = (
    "tests/gates/geometry/test_transformation.py::"
    "test_o4_a_duplicated_point_is_refused_and_guard_two_is_the_one_that_fires",
    "tests/gates/geometry/test_transformation.py::"
    "test_o6_the_weighted_question_is_permutes_plus_exactly_one_weight_guard",
)
_HOMOMORPHISM = _MANIFOLD_ROW + "test_the_induced_action_is_a_GROUP_HOMOMORPHISM_on_the_chart"
_DOWNWARD = "tests/gates/geometry/test_symmetry.py::test_invariance_is_downward_closed"


# ---------------------------------------------------------------------------
# Helpers — nothing here calls production at collection time.
# ---------------------------------------------------------------------------


def _api(dotted: str):
    """Resolve a post-S5 name at RUN time, or fail the row naming it."""
    module, _, attr = dotted.rpartition(".")
    try:
        obj = importlib.import_module(module)
    except ModuleNotFoundError:
        pytest.fail(f"S5 API absent: module {module}")
    if not hasattr(obj, attr):
        pytest.fail(f"S5 API absent: {dotted}")
    return getattr(obj, attr)


def _method(obj, name: str):
    if not hasattr(obj, name):
        pytest.fail(f"S5 API absent: {type(obj).__name__}.{name}")
    return getattr(obj, name)


def _params(text: str) -> dict:
    return {k: (int(v) if v.lstrip("-").isdigit() else v)
            for k, v in (kv.split("=") for kv in text.split(","))}


@functools.cache
def _built(label: str):
    """(measure, quadrature or None) for a baseline label — rebuilt, never stored."""
    kind, *rest = label.split("|")
    if kind == "registry":
        name, params = rest
        spec = {s.name: s for s in quadrature_registry}[name]
        q = spec.factory(**_params(params))
        return q.measure, q
    if kind == "folded_product":
        q = Quadrature.folded_product(**_params(rest[0]))
        return q.measure, q
    if kind == "polar_orbit":
        p = _params(rest[0])
        return gauss_legendre_on_polar_orbit(p["n"], p["axis"]), None
    if kind == "mirror_fold":
        builder, params, by = rest
        q = getattr(Quadrature, builder)(**_params(params))
        return q.measure.quotient(SubgroupOfO3.Mirror(by.removeprefix("by=sigma_"))), None
    if kind == "chart_mu":
        return gauss_legendre_on_mu(_params(rest[0])["n"]), None
    raise KeyError(label)


class _GroupsByName(dict):
    """Name -> SubgroupOfO3, resolving the parametrised families (``C_n``, ``D_nh``) on demand."""

    def __missing__(self, name: str) -> SubgroupOfO3:
        m = re.fullmatch(r"C_(\d+)|D_(\d+)h", name)
        if m is None:
            raise KeyError(name)
        g = SubgroupOfO3.Cn(int(m[1])) if m[1] else SubgroupOfO3.Dnh(int(m[2]))
        assert g.name == name
        self[name] = g
        return g


@functools.cache
def _groups() -> _GroupsByName:
    universe = (
        [SubgroupOfO3.Trivial, SubgroupOfO3.OctahedralOh, SubgroupOfO3.IcosahedralIh,
         SubgroupOfO3.Dinfh, SubgroupOfO3.SO3, SubgroupOfO3.O3]
        + [f(a) for f in (SubgroupOfO3.Mirror, SubgroupOfO3.SO2, SubgroupOfO3.O2) for a in "xyz"]
        + [SubgroupOfO3.Cn(n) for n in range(2, 9)]
        + [SubgroupOfO3.Dnh(n) for n in range(1, 9)]
    )
    return _GroupsByName({g.name: g for g in universe})


def _motions() -> dict[str, RigidMotion]:
    out = {f"mirror_{a}": mirror_motion(a, 3) for a in "xyz"}
    out["pole_deck_x"] = mirror_motion("x", 3)
    out["quarter_turn_z"] = RigidMotion.rotation(
        plane=(np.array([1.0, 0, 0]), np.array([0, 1.0, 0])), angle=np.pi / 2,
    )
    return out


def _fingerprint(measure) -> str:
    import hashlib

    h = hashlib.sha256()
    h.update(np.ascontiguousarray(measure.nodes, dtype=float).tobytes())
    h.update(np.ascontiguousarray(measure.weights, dtype=float).tobytes())
    return h.hexdigest()[:16]


def _perm(pi) -> list[int] | None:
    return None if pi is None else [int(i) for i in pi.indices]


_ANSWERED = sorted(k for k, r in RULES.items() if not k.startswith("chart_mu"))
_REFUSED = sorted(k for k in RULES if k.startswith("chart_mu"))


def test_the_baseline_population_is_what_the_spec_states():
    """X2: the fixture's denominators, so a truncated file cannot read as a clean pin."""
    assert BASELINE["commit"] == "c27830f5"
    assert len(RULES) == BASELINE["n_rules"] == 104
    assert len(FACTORS) == BASELINE["n_product_factor_rows"] == 42
    assert sum(len(r["is_invariant_under"]) for r in RULES.values()) == 3242
    finite = sum(len(r["certificate_permutations"]) for r in RULES.values())
    certified = sum(
        v is not None for r in RULES.values() for v in r["certificate_permutations"].values()
    )
    assert (finite, certified) == (2306, 914)
    assert len(_ANSWERED) == 101 and len(_REFUSED) == 3


def test_the_refusal_fragments_are_disjoint():
    frags = [FRAG_NOT_GSPACE, FRAG_NO_O3, FRAG_DOES_NOT_MAP, FRAG_NOT_NORMALISE,
             FRAG_WRONG_DIMENSION]
    for a in frags:
        for b in frags:
            assert a == b or a not in b


# ---------------------------------------------------------------------------
# A. The discriminator and the refusals (design gates 1 and 2)
# ---------------------------------------------------------------------------


def _verbs():
    """The five measure verbs, each as an O(3) question, keyed by name."""
    return {
        "is_invariant_under": lambda m: m.is_invariant_under(SubgroupOfO3.SO2("x")),
        "certificate_under": lambda m: m.certificate_under(SubgroupOfO3.Mirror("x")),
        "permutation_under": lambda m: m.permutation_under(mirror_motion("x", 3)),
        "singular_set_under": lambda m: m.singular_set_under(SubgroupOfO3.Mirror("x")),
        "symmetry_groups": lambda m: m.symmetry_groups(),
    }


_VERB_NAMES = ("is_invariant_under", "certificate_under", "permutation_under",
               "singular_set_under", "symmetry_groups")


class TestTheChartRuleIsRefused:
    """Gate 1. ``gauss_legendre_on_mu(4)`` lives on ``Interval[-1,1]``, a G-space of the
    1-D motions only. Today the kernel zero-pads it into R^3 and answers (``SO2_x``
    True, ``SO2_y`` False: the x is the padding's column 0). A design that kept the
    padding passes every sphere-family row and fails these."""

    @pytest.mark.parametrize("verb", _VERB_NAMES)
    def test_every_verb_refuses_an_O3_question_on_the_chart(self, verb):
        m = gauss_legendre_on_mu(4)
        with pytest.raises(TypeError, match=r"\[-1,1\]") as info:
            _verbs()[verb](m)
        assert FRAG_NO_O3 in str(info.value)
        assert "Interval" in str(info.value)

    @pytest.mark.parametrize("label", ["chart_mu|n=2", "chart_mu|n=4", "chart_mu|n=8"])
    def test_the_padded_verdicts_of_the_baseline_are_refused_not_re_answered(self, label):
        """Every True the baseline recorded on a chart rule (9 of 25 per rule) is a padding
        artefact; after S5 none of the 25 questions has an answer."""
        m, _ = _built(label)
        assert sum(RULES[label]["is_invariant_under"].values()) == 9  # the record
        for name in RULES[label]["is_invariant_under"]:
            with pytest.raises(TypeError, match=re.escape(FRAG_NO_O3)):
                m.is_invariant_under(_groups()[name])

    @pytest.mark.parametrize("case", ["unit_interval", "real_line_1", "real_space_2", "circle"])
    def test_a_G_space_of_another_dimension_is_refused_naming_it(self, case):
        m = {
            "unit_interval": lambda: equispaced(0.0, 1.0, 4),
            "real_line_1": lambda: DiscreteMeasure(
                nodes=np.array([-0.5, 0.25, 0.5]), weights=np.ones(3), support=RealSpace(1)),
            "real_space_2": lambda: DiscreteMeasure(
                nodes=np.array([[1.0, 0], [0, 1.0], [-1.0, 0], [0, -1.0]]),
                weights=np.ones(4), support=RealSpace(2)),
            "circle": lambda: periodic_trapezoid(8, shift=NODE_ALIGNED),
        }[case]()
        with pytest.raises(TypeError, match=re.escape(FRAG_NO_O3)) as info:
            m.is_invariant_under(SubgroupOfO3.Mirror("x"))
        assert m.support.name in str(info.value)

    def test_the_control_the_slab_rule_on_its_orbit_space_still_answers(self):
        """The positive leg: the SAME nodes on S^2/O2_x are asked and answer, so the
        refusal above is keyed on the support, not on the nodes."""
        m = Quadrature.gauss_legendre(4).measure
        assert m.is_invariant_under(SubgroupOfO3.SO2("x"))
        assert not m.is_invariant_under(SubgroupOfO3.SO2("y"))


class TestTheNonGSpaceSupportsAreRefused:
    """Gate 2, ruling 43: EnergyGroups, IndexSet, Product and FundamentalDomain carry
    no action; asking raises ``TypeError`` naming the support."""

    @staticmethod
    def _measure(case):
        if case == "energy":
            return DiscreteMeasure(nodes=np.arange(4.0), weights=np.ones(4), support=EnergyGroups(4))
        if case == "index":
            return DiscreteMeasure(nodes=np.arange(3.0), weights=np.ones(3), support=IndexSet("cells", 3))
        if case == "product":
            return gauss_legendre_on_mu(2) * periodic_trapezoid(4, shift=NODE_ALIGNED)
        if case == "fundamental_domain":
            fd = SPHERE.quotient(SubgroupOfO3.Mirror("y")).fundamental_domain
            assert fd is not None
            return DiscreteMeasure(
                nodes=np.array([[1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]), weights=np.ones(2), support=fd)
        raise KeyError(case)

    @pytest.mark.parametrize("case", ["energy", "index", "product", "fundamental_domain"])
    @pytest.mark.parametrize("verb", _VERB_NAMES)
    def test_asking_raises_TypeError_naming_the_support(self, case, verb):
        m = self._measure(case)
        with pytest.raises(TypeError, match=FRAG_NOT_GSPACE) as info:
            _verbs()[verb](m)
        assert m.support.name in str(info.value)
        assert FRAG_NO_O3 not in str(info.value)

    @pytest.mark.parametrize("case", ["energy", "index", "product", "fundamental_domain"])
    def test_the_support_has_no_quotient_verb(self, case):
        """S5 ruling (a): a Quotient's base is a G-space, so ``quotient`` lives on
        ``GSpace``; asking an EnergyGroups for an orbit space is not spellable."""
        assert not hasattr(self._measure(case).support, "quotient")
        if case == "energy":
            assert not hasattr(ENERGY, "quotient")

    @pytest.mark.parametrize("case", ["energy", "index", "product", "fundamental_domain"])
    def test_the_support_has_no_action_verb(self, case):
        gspace = _api("orpheus.geometry.manifold.GSpace")
        support = self._measure(case).support
        assert not isinstance(support, gspace)
        assert not hasattr(support, "action")


# ---------------------------------------------------------------------------
# B. The G-space family (design §2.1, gate 5)
# ---------------------------------------------------------------------------


def _rotation(dim: int, angle: float) -> RigidMotion:
    e = np.eye(dim)
    return RigidMotion.rotation(plane=(e[0], e[1]), angle=angle)


def _x_mirror(dim: int) -> RigidMotion:
    return RigidMotion.reflection(normal=np.eye(dim)[0])


class TestTheGSpaceFamily:
    """Each member's action and refusals (spec §4). Supports: Sphere, Ball, RealSpace,
    Circle, Interval and Quotient are G-spaces; the action is a ManifoldMap from the
    support to itself."""

    @pytest.mark.parametrize("member", ["sphere", "ball", "real_1", "real_2", "real_3",
                                        "circle", "interval", "quotient"])
    def test_the_six_families_are_G_spaces(self, member):
        gspace = _api("orpheus.geometry.manifold.GSpace")
        support = {
            "sphere": SPHERE, "ball": Ball(3), "real_1": RealSpace(1), "real_2": RealSpace(2),
            "real_3": RealSpace(3), "circle": CIRCLE, "interval": COSINE_INTERVAL,
            "quotient": SPHERE.quotient(SubgroupOfO3.Mirror("y")),
        }[member]
        assert isinstance(support, gspace)

    @pytest.mark.parametrize("support", ["sphere", "ball"])
    def test_sphere_and_ball_apply_the_3D_motion(self, support):
        s = {"sphere": SPHERE, "ball": Ball(3)}[support]
        g = _rotation(3, 0.7)
        pts = Quadrature.lebedev(7).measure.nodes * (1.0 if support == "sphere" else 0.5)
        got = np.asarray(_method(s, "action")(g)(pts))
        np.testing.assert_array_equal(got, g.on_points(pts))

    @pytest.mark.parametrize("support", ["sphere", "ball"])
    @pytest.mark.parametrize("dim", [1, 2])
    def test_sphere_and_ball_refuse_a_motion_of_another_dimension(self, support, dim):
        s = {"sphere": SPHERE, "ball": Ball(3)}[support]
        with pytest.raises(ValueError, match=FRAG_WRONG_DIMENSION):
            _method(s, "action")(_x_mirror(dim))

    def test_the_sphere_refuses_a_translation(self):
        """A wrap is an isometry of R^3 that does not map S^2 onto itself. Not in the
        design's table (it names only the dimension): spec §4, boundary B-T."""
        wrap = translation_motion("xmax", 3)
        assert np.any(wrap.translation != 0)
        with pytest.raises(ValueError, match=FRAG_DOES_NOT_MAP):
            _method(SPHERE, "action")(wrap)

    @pytest.mark.parametrize("d", [1, 2, 3])
    def test_real_space_acts_by_O_d_and_refuses_another_dimension(self, d):
        s = RealSpace(d)
        pts = np.linspace(-1.0, 1.0, 5 * d).reshape(5, d)
        got = np.asarray(_method(s, "action")(_x_mirror(d))(pts))
        np.testing.assert_array_equal(got.reshape(5, d), pts * np.r_[-1.0, np.ones(d - 1)])
        other = 3 if d != 3 else 2
        with pytest.raises(ValueError, match=FRAG_WRONG_DIMENSION):
            _method(s, "action")(_x_mirror(other))

    def test_the_circle_acts_by_O2_and_refuses_another_dimension(self):
        nodes = periodic_trapezoid(8, shift=NODE_ALIGNED).nodes
        g = _rotation(2, 2 * np.pi / 8)
        got = np.asarray(_method(CIRCLE, "action")(g)(nodes))
        np.testing.assert_array_equal(got, g.on_points(nodes))
        for dim in (1, 3):
            with pytest.raises(ValueError, match=FRAG_WRONG_DIMENSION):
                _method(CIRCLE, "action")(_x_mirror(dim))

    def test_the_symmetric_interval_carries_mu_to_minus_mu(self):
        mu = gauss_legendre_on_mu(8).nodes
        got = np.asarray(_method(COSINE_INTERVAL, "action")(_x_mirror(1))(mu)).reshape(-1)
        np.testing.assert_array_equal(got, -mu)

    @pytest.mark.parametrize("interval", ["unit", "half_line", "lopsided"])
    def test_an_asymmetric_interval_refuses_x_to_minus_x_naming_itself(self, interval):
        s = {"unit": UNIT_INTERVAL, "half_line": HALF_LINE, "lopsided": Interval(-1.0, 2.0)}[interval]
        with pytest.raises(ValueError, match=FRAG_DOES_NOT_MAP) as info:
            _method(s, "action")(_x_mirror(1))
        assert s.name in str(info.value)

    def test_the_real_line_is_symmetric(self):
        x = np.array([-3.0, 0.5, 7.0])
        got = np.asarray(_method(REAL_LINE, "action")(_x_mirror(1))(x)).reshape(-1)
        np.testing.assert_array_equal(got, -x)

    @pytest.mark.parametrize("interval", ["cosine", "unit", "half_line", "real_line"])
    def test_the_identity_acts_on_every_interval(self, interval):
        s = {"cosine": COSINE_INTERVAL, "unit": UNIT_INTERVAL, "half_line": HALF_LINE,
             "real_line": REAL_LINE}[interval]
        x = np.array([0.0, 0.25, 0.5])
        got = np.asarray(_method(s, "action")(RigidMotion.identity(1))(x)).reshape(-1)
        np.testing.assert_array_equal(got, x)

    def test_the_unit_interval_carries_its_own_isometry(self):
        """[0,1]'s isometry group is {id, x -> 1 - x}, which IS a RigidMotion (linear
        -1, translation 1). The design's table names only x -> -x; the mathematically
        right rule admits this one (spec §4, boundary B-I; a ruling, NEEDS 3)."""
        flip = RigidMotion(linear=np.array([[-1.0]]), translation=np.array([1.0]))
        x = np.array([0.0, 0.25, 1.0])
        got = np.asarray(_method(UNIT_INTERVAL, "action")(flip)(x)).reshape(-1)
        np.testing.assert_array_equal(got, 1.0 - x)

    def test_an_interval_refuses_a_3D_motion(self):
        with pytest.raises(ValueError, match=FRAG_WRONG_DIMENSION):
            _method(COSINE_INTERVAL, "action")(mirror_motion("x", 3))

    @pytest.mark.rests_on(_MANIFOLD_ROW + "test_sigma_x_on_the_sigma_y_fold_is_the_disks_x_flip")
    def test_the_quotient_action_is_the_renamed_induced_action(self):
        """sigma_x on S^2/sigma_y is the chart's x-flip — a closed form, not the old method."""
        fold = SPHERE.quotient(SubgroupOfO3.Mirror("y"))
        rng = np.random.default_rng(11)
        rho, phi = np.sqrt(rng.uniform(0, 1, 32)), rng.uniform(0, 2 * np.pi, 32)
        chart = np.column_stack([rho * np.cos(phi), rho * np.sin(phi)])
        got = np.asarray(_method(fold, "action")(mirror_motion("x", 3))(chart))
        np.testing.assert_array_equal(got, chart * np.array([-1.0, 1.0]))

    def test_the_quotient_refuses_a_translation(self):
        """Boundary B-T on an orbit space: the normaliser test reads the LINEAR part, so
        a wrap 'normalises' O2_x and today's induced action translates the barycentres
        off the ball. The action must refuse it, as the sphere's does."""
        entry = SPHERE.quotient(SubgroupOfO3.O2("x"))
        with pytest.raises(ValueError, match=FRAG_DOES_NOT_MAP):
            _method(entry, "action")(translation_motion("xmax", 3))

    def test_the_quotient_refuses_a_non_normaliser_with_the_pinned_fragment(self):
        fold = SPHERE.quotient(SubgroupOfO3.Mirror("y"))
        with pytest.raises(ValueError, match=FRAG_NOT_NORMALISE):
            _method(fold, "action")(_motions()["quarter_turn_z"])

    def test_induced_action_is_retired(self):
        """One name for one verb across the family (ruling 43)."""
        _method(SPHERE.quotient(SubgroupOfO3.O2("x")), "action")
        assert not hasattr(Quotient, "induced_action")


class TestAmbientPoints:
    """Ruling 46: ``GSpace.ambient_points`` is TOTAL. The identity on every G-space that is
    not a quotient (an Interval's points in R^1, a Circle's in R^2); the orbit barycentre in
    the base's ambient space on a Quotient. Expected values are closed forms, not
    ``orbit_barycentres``."""

    @pytest.mark.parametrize("member", ["sphere", "ball", "real_1", "real_2", "real_3",
                                        "interval", "circle"])
    def test_the_identity_off_a_quotient(self, member):
        support, pts = {
            "sphere": (SPHERE, Quadrature.lebedev(7).measure.nodes),
            "ball": (Ball(3), 0.5 * Quadrature.lebedev(7).measure.nodes),
            "real_1": (RealSpace(1), np.array([[-2.0], [0.5], [3.0]])),
            "real_2": (RealSpace(2), np.array([[1.0, -2.0], [0.0, 3.0]])),
            "real_3": (RealSpace(3), np.array([[1.0, -2.0, 4.0]])),
            "interval": (COSINE_INTERVAL, gauss_legendre_on_mu(4).nodes),
            "circle": (CIRCLE, periodic_trapezoid(8, shift=NODE_ALIGNED).nodes),
        }[member]
        got = np.asarray(_method(support, "ambient_points")(pts))
        np.testing.assert_array_equal(got.reshape(np.shape(pts)), pts)

    def test_the_barycentre_on_S2_mod_sigma_y(self):
        """A representative p of an orbit {p, sigma_y p} has barycentre (x, 0, z)."""
        fold = SPHERE.quotient(SubgroupOfO3.Mirror("y"))
        reps = Quadrature.folded_product(4, 8).measure.nodes
        assert np.all(reps[:, 1] != 0)
        want = reps * np.array([1.0, 0.0, 1.0])
        np.testing.assert_array_equal(np.asarray(_method(fold, "ambient_points")(reps)), want)

    def test_the_barycentre_on_S2_mod_O2x(self):
        """A chart point mu of S^2/O2_x has barycentre (mu, 0, 0), the centre of its circle."""
        entry = SPHERE.quotient(SubgroupOfO3.O2("x"))
        mu = gauss_legendre_on_polar_orbit(4, "x").nodes
        want = np.column_stack([mu, np.zeros(4), np.zeros(4)])
        np.testing.assert_array_equal(np.asarray(_method(entry, "ambient_points")(mu)), want)


def test_the_polar_mirror_spelled_twice_agrees():
    """Ruling 47's commuting square (X4): on the chart of S^2/O2_x, the 1-D mirror acting
    on Interval(-1,1) and sigma_x acting on the quotient are one map, bit for bit."""
    entry = SPHERE.quotient(SubgroupOfO3.O2("x"))
    mu = gauss_legendre_on_polar_orbit(8, "x").nodes
    on_interval = np.asarray(_method(COSINE_INTERVAL, "action")(_x_mirror(1))(mu)).reshape(-1)
    on_quotient = np.asarray(_method(entry, "action")(mirror_motion("x", 3))(mu)).reshape(-1)
    assert np.any(on_interval != mu)          # the map moves the points: not vacuous
    np.testing.assert_array_equal(on_interval, on_quotient)


def test_spherical_product_refuses_a_polar_factor_off_the_chart_FIRST():
    """F6 / ruling 47, a discrimination row: ``spherical_product`` refuses a polar factor
    on S^2/O2_x with its OWN message, before the group derivation can raise the measure's
    TypeError, and before the pushforward's refusal. Positive leg: the chart factor builds."""
    from orpheus.numerics.quadrature.rules_product import spherical_product

    azimuthal = periodic_trapezoid(8, shift=NODE_ALIGNED)
    accepted, _ = spherical_product(gauss_legendre_on_mu(4), azimuthal)
    assert accepted.support == SPHERE
    with pytest.raises(ValueError, match=FRAG_POLAR_SUPPORT) as info:
        spherical_product(gauss_legendre_on_polar_orbit(4, "x"), azimuthal)
    assert COSINE_INTERVAL.name in str(info.value)
    assert FRAG_PUSHFORWARD_DOMAIN not in str(info.value)


class TestPointIdentityOnAQuotient:
    """``GSpace.coordinates``: two points are equal iff their coordinates are."""

    def test_on_a_bare_sphere_the_coordinates_are_the_points(self):
        pts = Quadrature.lebedev(7).measure.nodes
        np.testing.assert_array_equal(np.asarray(_method(SPHERE, "coordinates")(pts)), pts)

    def test_two_representatives_of_one_sigma_y_orbit_are_one_point(self):
        fold = SPHERE.quotient(SubgroupOfO3.Mirror("y"))
        reps = Quadrature.folded_product(4, 8).measure.nodes
        mirrored = reps * np.array([1.0, -1.0, 1.0])
        coords = _method(fold, "coordinates")
        np.testing.assert_array_equal(np.asarray(coords(reps)), np.asarray(coords(mirrored)))

    def test_distinct_orbits_stay_distinct(self):
        """Injectivity on the shipped fold: 16 representatives, 16 distinct coordinates."""
        fold = SPHERE.quotient(SubgroupOfO3.Mirror("y"))
        reps = Quadrature.folded_product(4, 8).measure.nodes
        c = np.asarray(_method(fold, "coordinates")(reps))
        assert len({tuple(row) for row in c.reshape(len(reps), -1)}) == len(reps) == 16

    def test_a_point_and_its_rotation_about_the_axis_are_one_point_of_S2_mod_O2x(self):
        entry = SPHERE.quotient(SubgroupOfO3.O2("x"))
        p = Quadrature.lebedev(7).measure.nodes
        rotated = p @ _rotation_about("x", 1.1).T
        coords = _method(entry, "coordinates")
        np.testing.assert_array_equal(np.asarray(coords(p)), np.asarray(coords(rotated)))


def _rotation_about(axis: str, angle: float) -> np.ndarray:
    i = "xyz".index(axis)
    j, k = [a for a in range(3) if a != i]
    m = np.eye(3)
    c, s = np.cos(angle), np.sin(angle)
    m[j, j], m[j, k], m[k, j], m[k, k] = c, -s, s, c
    return m


# ---------------------------------------------------------------------------
# C. relabelling_to — the measure's one group-free verb (design §2.2)
# ---------------------------------------------------------------------------


def _toy(nodes, weights, support=SPHERE) -> DiscreteMeasure:
    return DiscreteMeasure(nodes=np.asarray(nodes, float), weights=np.asarray(weights, float),
                           support=support)


_P = np.array([[1.0, 0, 0], [0, 1.0, 0], [0, 0, 1.0], [-1.0, 0, 0]])


@pytest.mark.rests_on(*_MATCHER_ROWS)
class TestRelabellingTo:
    """Convention: ``mu.relabelling_to(nu)`` is pi with ``nu``'s atom i equal to ``mu``'s
    atom pi(i), points compared in ``support.coordinates``, weights compared, None if no
    bijection. A 3-cycle is used throughout: it is not its own inverse, so a pi / pi^-1
    slip reddens."""

    def test_a_reordering_is_returned_in_the_stated_direction(self):
        w = np.array([0.1, 0.2, 0.3, 0.4])
        mu = _toy(_P, w)
        sigma = np.array([1, 2, 0, 3])          # nu_i = mu_sigma(i)
        nu = _toy(_P[sigma], w[sigma])
        pi = _method(mu, "relabelling_to")(nu)
        assert pi is not None
        np.testing.assert_array_equal(pi.indices, sigma)

    def test_the_weight_leg_compares_the_OTHER_measures_weights(self):
        """The design names ``permutation_preserving`` as the matcher; it compares
        ``w[pi]`` with ``w`` for ONE weight array (a pushforward's assumption). A general
        ``nu`` carries its own weights: with distinct weights and a 3-cycle, both naive
        spellings (``mu.weights`` or ``nu.weights`` as the one array) answer None here."""
        w = np.array([0.1, 0.2, 0.3, 0.4])
        sigma = np.array([1, 2, 0, 3])
        nu = _toy(_P[sigma], w[sigma])
        assert _method(_toy(_P, w), "relabelling_to")(nu) is not None

    def test_nodes_reordered_without_their_weights_is_not_a_relabelling(self):
        w = np.array([0.1, 0.2, 0.3, 0.4])
        sigma = np.array([1, 2, 0, 3])
        nu = _toy(_P[sigma], w)                  # weights left in place
        assert _method(_toy(_P, w), "relabelling_to")(nu) is None

    def test_the_bijection_guard_ERR_073(self):
        """mu = {a, a, b}, nu = {a, b, b}: every atom of nu has a partner in mu, no
        bijection exists."""
        a, b = [1.0, 0, 0], [0, 1.0, 0]
        mu, nu = _toy([a, a, b], [1, 1, 1]), _toy([a, b, b], [1, 1, 1])
        assert _method(mu, "relabelling_to")(nu) is None

    @pytest.mark.parametrize("delta,matches", [(1e-15, True), (1e-9, False)])
    def test_the_weight_window(self, delta, matches):
        w = np.array([0.1, 0.2, 0.3, 0.4])
        nu = _toy(_P, w + np.r_[delta, 0, 0, 0])
        assert (_method(_toy(_P, w), "relabelling_to")(nu) is not None) is matches

    @pytest.mark.parametrize("scale,matches", [(0.5, True), (2.0, False)])
    def test_the_node_window_is_the_weight_window_times_the_factor(self, scale, matches):
        window = _api("orpheus.numerics.measure.WEIGHT_ATOL") * 100.0
        shifted = _P + np.array([[scale * window, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0]])
        nu = DiscreteMeasure(nodes=shifted, weights=np.ones(4), support=RealSpace(3))
        mu = DiscreteMeasure(nodes=_P, weights=np.ones(4), support=RealSpace(3))
        assert (_method(mu, "relabelling_to")(nu) is not None) is matches

    def test_a_different_length_is_None(self):
        mu = _toy(_P, np.ones(4))
        assert _method(mu, "relabelling_to")(_toy(_P[:3], np.ones(3))) is None

    def test_a_different_support_is_refused(self):
        mu = _toy(_P, np.ones(4))
        nu = DiscreteMeasure(nodes=_P, weights=np.ones(4), support=RealSpace(3))
        with pytest.raises(ValueError, match=FRAG_SUPPORT):
            _method(mu, "relabelling_to")(nu)

    def test_on_a_fold_two_representatives_of_one_orbit_are_one_atom(self):
        """The point-identity leg: replacing a representative by its sigma_y mate gives
        the IDENTITY relabelling, where a raw-node matcher answers None."""
        mu = Quadrature.folded_product(4, 8).measure
        nodes = mu.nodes.copy()
        nodes[3] = nodes[3] * np.array([1.0, -1.0, 1.0])
        assert not np.array_equal(nodes, mu.nodes)
        nu = DiscreteMeasure(nodes=nodes, weights=mu.weights, support=mu.support)
        pi = _method(mu, "relabelling_to")(nu)
        assert pi is not None
        np.testing.assert_array_equal(pi.indices, np.arange(mu.n_points))


# ---------------------------------------------------------------------------
# D. Agreement with the record (design gate 3) — bit identity on permutations
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("label", _ANSWERED)
class TestAgreementWithTheRecord:
    """Every (rule, finite group) certificate, verdict, singular set, maximal group set
    and motion permutation equals the old kernel's, bit for bit (``==`` on index lists).
    A red in ``test_the_rule_did_not_move`` invalidates the label's other rows."""

    def test_the_rule_did_not_move(self, label):
        m, _ = _built(label)
        assert (m.support.name, m.n_points, _fingerprint(m)) == (
            RULES[label]["support"], RULES[label]["n_points"], RULES[label]["fingerprint"])

    def test_certificate_permutations(self, label):
        m, _ = _built(label)
        bad = []
        for name, want in RULES[label]["certificate_permutations"].items():
            cert = m.certificate_under(_groups()[name])
            got = None if cert is None else [[int(i) for i in p.indices] for p in cert.permutations]
            if got != want:
                bad.append(name)
        assert not bad, f"{label}: certificates differ under {bad}"

    def test_verdicts(self, label):
        m, _ = _built(label)
        got = {name: bool(m.is_invariant_under(_groups()[name]))
               for name in RULES[label]["is_invariant_under"]}
        diff = {k for k in got if got[k] != RULES[label]["is_invariant_under"][k]}
        assert not diff, f"{label}: verdicts differ under {sorted(diff)}"

    def test_singular_sets(self, label):
        m, _ = _built(label)
        for name, want in RULES[label]["singular_set_under"].items():
            assert [int(i) for i in m.singular_set_under(_groups()[name])] == want, name

    @pytest.mark.parametrize("method", ["walk", "bruteforce"])
    def test_symmetry_groups(self, label, method):
        m, _ = _built(label)
        got = sorted(g.name for g in m.symmetry_groups(method=method))
        assert got == RULES[label][f"symmetry_groups_{method}"]

    def test_motion_permutations(self, label):
        m, q = _built(label)
        for key, motion in _motions().items():
            assert _perm(m.permutation_under(motion.linear_part)) == RULES[label]["permutation_under"][key], key
            if q is not None:
                assert _perm(q.ordinate_permutation(motion)) == RULES[label]["ordinate_permutation"][key], key


@pytest.mark.parametrize("label", sorted(FACTORS))
def test_the_product_factor_checks_through_the_measure(label):
    """Gate 3, second half: ``_derived_product_group``'s three questions, asked through
    the measure verb on the factors' own G-spaces (S^1 and Interval[-1,1]), give the
    permutations ``RigidMotion.preserves`` gave. The window moves (node 1e-13 to 1e-11)."""
    kind, params = label.split("|")
    p = _params(params)
    shift = STAGGERED if kind.startswith("folded") else NODE_ALIGNED
    polar, azimuthal = gauss_legendre_on_mu(p["n_mu"]), periodic_trapezoid(p["n_phi"], shift=shift)
    step = RigidMotion.rotation(plane=(np.eye(2)[0], np.eye(2)[1]), angle=2.0 * np.pi / p["n_phi"])
    want = FACTORS[label]
    assert _perm(azimuthal.permutation_under(step)) == want["rotation"]
    assert _perm(azimuthal.permutation_under(RigidMotion.reflection(normal=np.eye(2)[1]))) == want["vertical_mirror"]
    assert _perm(polar.permutation_under(_x_mirror(1))) == want["polar_mirror"]


# ---------------------------------------------------------------------------
# E. The REFERENCE: an independent matcher (anchors the record)
# ---------------------------------------------------------------------------


def _hungarian(points: np.ndarray, images: np.ndarray, weights: np.ndarray) -> list[int] | None:
    """pi with images[i] == points[pi[i]], weights equal; scipy's assignment, no ORPHEUS code."""
    cost = np.linalg.norm(images[:, None, :] - points[None, :, :], axis=-1)
    cost = cost + 1e3 * (np.abs(weights[:, None] - weights[None, :]) > 1e-13)
    rows, cols = linear_sum_assignment(cost)
    if cost[rows, cols].max() > 1e-10:
        return None
    return [int(c) for c in cols]


def _invariant_coordinates(label: str, pts: np.ndarray) -> np.ndarray:
    """The orbit invariants, written from the group, not read off the Quotient."""
    support = RULES[label]["support"]
    if support == "S^2":
        return pts
    if support.startswith("S^2/sigma_"):
        a = "xyz".index(support[-1])
        out = pts.copy()
        out[:, a] = np.abs(out[:, a])
        return out
    raise KeyError(support)


_REFERENCE_LABELS = sorted(
    k for k in _ANSWERED if RULES[k]["support"] == "S^2" or RULES[k]["support"].startswith("S^2/sigma_"))
_POLAR_LABELS = sorted(k for k in _ANSWERED if RULES[k]["support"].startswith("S^2/O2_"))


@pytest.mark.parametrize("label", _REFERENCE_LABELS)
def test_the_certificate_against_an_independent_matcher(label):
    """REFERENCE: for every certified finite group, every element's permutation equals
    the Hungarian assignment of the element's raw matrix applied to the nodes, compared
    in the orbit invariants. Shared input (declared): the group's element matrices."""
    m, _ = _built(label)
    checked = 0
    for name, want in RULES[label]["certificate_permutations"].items():
        if want is None:
            continue
        cert = m.certificate_under(_groups()[name])
        assert cert is not None, name
        for g, perm in zip(cert.operators, cert.permutations):
            images = _invariant_coordinates(label, m.nodes @ g.linear.T)
            ref = _hungarian(_invariant_coordinates(label, m.nodes), images, m.weights)
            assert ref == [int(i) for i in perm.indices], (name, g)
            checked += 1
    assert checked > 0


@pytest.mark.parametrize("label", _POLAR_LABELS)
def test_the_polar_orbit_certificate_against_mu_to_sign_mu(label):
    """On S^2/O2_a a normaliser element maps e_a to +-e_a, so mu -> Q[a,a] mu."""
    m, _ = _built(label)
    a = "xyz".index(RULES[label]["support"][-1])
    mu = m.nodes.reshape(-1, 1)
    for name, want in RULES[label]["certificate_permutations"].items():
        if want is None:
            continue
        cert = m.certificate_under(_groups()[name])
        assert cert is not None, name
        for g, perm in zip(cert.operators, cert.permutations):
            ref = _hungarian(mu, g.linear[a, a] * mu, m.weights)
            assert ref == [int(i) for i in perm.indices], (name, g)


# ---------------------------------------------------------------------------
# F. The composition is the kernel (design gate 4) — counting spy and route
# ---------------------------------------------------------------------------


@pytest.fixture
def relabelling_spy(monkeypatch):
    original = _method(DiscreteMeasure, "relabelling_to")
    calls: list[int] = []

    def counting(self, other, *a, **k):
        calls.append(self.n_points)
        return original(self, other, *a, **k)

    monkeypatch.setattr(DiscreteMeasure, "relabelling_to", counting)
    assert getattr(DiscreteMeasure, "relabelling_to") is counting  # installed
    return calls


def test_ordinate_permutation_enters_relabelling(relabelling_spy):
    q = Quadrature.level_symmetric(4)
    assert q.ordinate_permutation(mirror_motion("x", 3)) is not None
    assert len(relabelling_spy) == 1


def test_the_fold_enters_relabelling_once_per_element(relabelling_spy):
    Quadrature.product(4, 8).quotient(SubgroupOfO3.Mirror("y"))
    assert len(relabelling_spy) >= 2


def test_derived_product_group_enters_relabelling_and_the_second_door_is_gone(relabelling_spy):
    """Route gate: the three factor questions go through the kernel (>= 3 entries), and
    the second door they used until S5 no longer exists. ``RigidMotion.preserves`` and
    ``permutation_preserving`` were retired with it (elegance review, 2026-10-19): with
    ``_derived_product_group`` re-pointed they had 0 production callers, and their weight
    leg compared ONE array (finding F1)."""
    from orpheus.geometry import transformation
    from orpheus.numerics.quadrature.rules_product import spherical_product

    spherical_product(gauss_legendre_on_mu(4), periodic_trapezoid(8, shift=NODE_ALIGNED))
    assert len(relabelling_spy) >= 3
    assert not hasattr(RigidMotion, "preserves")
    assert not hasattr(transformation, "permutation_preserving")


# ---------------------------------------------------------------------------
# G. The positive control (design gate 6)
# ---------------------------------------------------------------------------


def _ls8(perturb: float) -> DiscreteMeasure:
    m = Quadrature.level_symmetric(8).measure
    w = m.weights.copy()
    w[0] += perturb
    return DiscreteMeasure(nodes=m.nodes, weights=w, support=m.support)


@pytest.mark.parametrize("perturb,invariant", [(0.0, True), (1e-9, False)])
def test_one_LS8_weight_moved_by_1e_9_breaks_Oh(perturb, invariant):
    assert _ls8(perturb).is_invariant_under(SubgroupOfO3.OctahedralOh) is invariant


def test_the_control_reaches_the_composition(relabelling_spy):
    assert not _ls8(1e-9).is_invariant_under(SubgroupOfO3.OctahedralOh)
    assert len(relabelling_spy) >= 1


# ---------------------------------------------------------------------------
# H. Contracts the composition must keep (spec §5): None, not a raise
# ---------------------------------------------------------------------------


def test_a_non_normalising_motion_on_a_fold_is_None_not_a_raise():
    """``Quotient.action`` raises for a non-normaliser; ``permutation_under`` must still
    answer None (the ordinate_permutation callers own the refusal message)."""
    m = Quadrature.folded_product(4, 8).measure
    assert m.permutation_under(_motions()["quarter_turn_z"].linear_part) is None
    assert Quadrature.folded_product(4, 8).ordinate_permutation(_motions()["quarter_turn_z"]) is None


def test_a_non_normalising_group_is_False_not_a_raise():
    m = Quadrature.folded_product(4, 8).measure
    assert m.is_invariant_under(SubgroupOfO3.Cn(4)) is False
    assert m.certificate_under(SubgroupOfO3.Cn(4)) is None


@pytest.mark.parametrize("rule", ["gauss_legendre_8", "lebedev_7"])
def test_an_affine_motion_is_None_not_a_raise(rule):
    """The wrap translates every node off S^2: not a symmetry, answered None today, on
    the orbit space (GL8, S^2/O2_x) and on the bare sphere (Lebedev 7), where the
    sphere's action REFUSES the translation (``test_the_sphere_refuses_a_translation``)
    and ``permutation_under`` must turn that refusal into None."""
    m = {"gauss_legendre_8": lambda: Quadrature.gauss_legendre(8),
         "lebedev_7": lambda: Quadrature.lebedev(7)}[rule]().measure
    wrap = translation_motion("xmax", 3)
    assert m.permutation_under(wrap) is None


def test_relabelling_needs_no_group_two_permuted_energy_measures_relabel():
    """The user's ruling of 2026-10-19 (elegance review, blocker 1): the point identity
    is every manifold's, so ``relabelling_to`` is defined on a support no motion acts on.
    Until the fix it asked the support to be a G-space and raised here."""
    w = np.array([1.0, 2.0, 3.0, 4.0])
    sigma = np.array([2, 0, 3, 1])
    mu = DiscreteMeasure(nodes=np.arange(4.0), weights=w, support=EnergyGroups(4))
    nu = DiscreteMeasure(nodes=np.arange(4.0)[sigma], weights=w[sigma], support=EnergyGroups(4))
    pi = mu.relabelling_to(nu)
    assert pi is not None
    np.testing.assert_array_equal(pi.indices, sigma)


def test_an_internal_error_in_a_familys_law_is_not_read_as_no_symmetry(monkeypatch):
    """Elegance review, blocker 2: the refusal is a typed VALUE, so ``permutation_under``
    turns only a ``DoesNotMap`` / ``DoesNotNormalise`` into None. Until the fix an
    ``except ValueError`` around the action read ANY ValueError raised inside the law as
    "not a symmetry". The control: an honest refusal still answers None."""
    mirror = RigidMotion.reflection(normal=np.array([1.0]))
    assert equispaced(0.0, 1.0, 4).permutation_under(mirror) is None  # [0,1]: DoesNotMap

    def broken(self, motion):
        raise ValueError("an internal defect in the Interval law")

    monkeypatch.setattr(Interval, "_refusal", broken)
    with pytest.raises(ValueError, match="internal defect"):
        gauss_legendre_on_mu(4).permutation_under(mirror)


# ---------------------------------------------------------------------------
# I. The lattice in ``symmetry`` (design §2.3), measure-free
# ---------------------------------------------------------------------------


_LATTICE = (["Oh", "D_4h", "D_2h", "C_4", "C_2", "sigma_x", "Trivial"])


def _independent_maximal(true_names: set[str], cands: list[SubgroupOfO3]) -> set[str]:
    trues = [g for g in cands if g.name in true_names]
    return {g.name for g in trues if not any(h is not g and h.contains(g) and not g.contains(h) for h in trues)}


@pytest.mark.parametrize("true_names", [
    ("D_4h", "D_2h", "C_4", "C_2", "sigma_x", "Trivial"),       # a downward-closed ideal
    ("C_4", "C_2", "Trivial"),
    ("Trivial",),
])
def test_maximal_satisfying_returns_the_maximal_elements(true_names):
    maximal_satisfying = _api("orpheus.geometry.symmetry.maximal_satisfying")
    cands = [_groups()[n] for n in _LATTICE]
    truth = set(true_names)
    for method in ("walk", "bruteforce"):
        got = {g.name for g in maximal_satisfying(lambda g: g.name in truth, cands, method=method)}
        assert got == _independent_maximal(truth, cands), method


def test_walk_and_bruteforce_agree_on_a_predicate_that_is_NOT_downward_closed():
    """DECLARED BLINDNESS (spec §8, finding F2): the walk visits every candidate below a
    failing node, so walk == bruteforce for EVERY predicate (`[M]` 2000 of 2000 random
    non-downward-closed predicates on LS8's 32 candidates, 0 disagreements). Their
    agreement is therefore no evidence of downward closure; that law has its own gate
    (``test_symmetry.py::test_invariance_is_downward_closed``). This row pins the fact so
    nobody re-credits the agreement."""
    maximal_satisfying = _api("orpheus.geometry.symmetry.maximal_satisfying")
    cands = [_groups()[n] for n in _LATTICE]
    pred = lambda g: g.name in ("D_4h", "C_2")  # noqa: E731 — C_2 true, D_2h false: not an ideal
    walk = {g.name for g in maximal_satisfying(pred, cands, method="walk")}
    brute = {g.name for g in maximal_satisfying(pred, cands, method="bruteforce")}
    assert walk == brute == {"D_4h"}


@pytest.mark.rests_on(_DOWNWARD)
@pytest.mark.parametrize("label", _ANSWERED)
def test_the_record_is_downward_closed(label):
    """Anti-pattern #15's compatibility law on the record: H <= G and G-invariant imply
    H-invariant, over every (pair x rule) of the baseline. A RECORD-coherence row."""
    v = RULES[label]["is_invariant_under"]
    gs = _groups()
    for big, ok in v.items():
        if not ok:
            continue
        for small in v:
            if gs[big].contains(gs[small]):
                assert v[small], (label, big, small)


def test_candidate_groups_reads_ambient_points():
    """``symmetry.candidate_groups(points)`` on the ambient barycentres equals the record
    on every answered rule."""
    candidate_groups = _api("orpheus.geometry.symmetry.candidate_groups")
    # the node window the measure asks at: symmetry sits below the measure that
    # owns it, so the window is the caller's to state (orchestrator, 2026-10-19)
    window = _api("orpheus.numerics.measure.WEIGHT_ATOL") * 100.0
    for label in _ANSWERED:
        m, _ = _built(label)
        pts = _method(m.support, "ambient_points")(m.nodes)   # ruling 46
        assert sorted(g.name for g in candidate_groups(pts, atol=window)) == RULES[label]["candidates"], label


# ---------------------------------------------------------------------------
# J. The module is gone and the pieces are home (design §2.3)
# ---------------------------------------------------------------------------


def test_the_invariance_module_is_deleted():
    with pytest.raises(ModuleNotFoundError):
        importlib.import_module("orpheus.numerics.invariance")


@pytest.mark.parametrize("dotted", [
    "orpheus.geometry.symmetry.OrbitCertificate",
    "orpheus.geometry.symmetry.maximal_satisfying",
    "orpheus.geometry.symmetry.maximal_subgroups",
    "orpheus.geometry.symmetry.candidate_groups",
    "orpheus.numerics.measure.WEIGHT_ATOL",
    "orpheus.geometry.manifold.GSpace",
])
def test_each_piece_lives_where_the_design_puts_it(dotted):
    _api(dotted)


def test_one_weight_window_and_the_fold_defaults_to_it():
    import inspect

    weight_atol = _api("orpheus.numerics.measure.WEIGHT_ATOL")
    assert weight_atol == 1e-13
    assert inspect.signature(Quadrature.quotient).parameters["atol"].default == weight_atol
    assert inspect.signature(DiscreteMeasure.is_invariant_under).parameters["atol"].default == weight_atol
