r"""The manifold a measure lives on and a basis is defined over.

A :class:`Manifold` is the *point set* — level 1 of the three-level stack
that :doc:`the theory page </theory/foundations/spaces>` describes:

.. code-block:: text

    level 1   the manifold M            S^2, [-1,1], R^d, M/H, M x N
    level 2   fields on it              L^2(M)          <- FunctionSpace
    level 3   coefficients              R^K             <- a Basis's shape

Until 2026-08 level 1 was ``Space = str`` — an opaque tag whose own
docstring called the entries *"recommendations, not constraints"*.  Three
consequences, all measured, and all of them the reason this type exists:

* **Nothing could be refused.**  A 1-D quadrature could declare its nodes
  live on :math:`S^2` while supplying :math:`\mu_y = \mu_z = 0`, which is
  not a point of the sphere.  That forgery is the first link of
  **ERR-080** — a live wrong-answer defect in :math:`P_L` scattering.
  (Plain text, not a ``:ref:``: `[M]` the ``error-entry`` directive emits a
  ``container`` + ``rubric`` and **no** ``nodes.target``, so there is no
  label to resolve — ``vv-error-ERR-080`` occurs nowhere else in ``docs/``
  or ``orpheus/``.  Every one of the ~150 ERR citations in ``orpheus/`` is
  plain text for the same reason.)  :meth:`Manifold.contains` refuses it at
  construction, three hops before the symptom.
* **The algebra ran as string concatenation.**  ``f"{a} x {b}"`` was the
  product, ``f"{a}/{H.name}"`` the quotient, ``f"phi_*({a})"`` the
  pushforward codomain.  Those are :meth:`__mul__`, :meth:`quotient` and a
  :class:`ManifoldMap`'s codomain; the interpolation *was* the operation,
  unnamed.
* **The vocabulary drifted.**  ``'S^2/<sigma_y>'`` and ``'S^2/sigma_y'``
  both shipped, naming one quotient and unequal under ``==``.

The type is a **closed sum**: a small frozen variant set, with the
operations every manifold answers (:attr:`dim`, :attr:`ambient_dim`,
:meth:`contains`, :meth:`__mul__`, :attr:`name`) on the base, and the operations only a
*quotient* can answer on :class:`Quotient` alone.  A sphere has no syzygy
ideal, and under this shape it cannot be asked for one.

⚠ **The import direction is** ``manifold -> symmetry`` **, at module scope,
since R2 of #434 (2026-09-03).**  A :class:`Quotient` is a manifold and a
group, so this module imports :class:`~orpheus.geometry.symmetry.SubgroupOfO3`
like any other dependency (and, since layer DAG WP-B of 2026-10-22, the axis
table from :mod:`orpheus.geometry.faces`, its one home); ``symmetry`` imports
only :mod:`orpheus.geometry.transformation`, the axis table and the roots of
unity.  Until R2 the direction ran the other
way — ``symmetry`` held the invariance kernel and therefore imported
``measure``, which imports this module, so a module-scope ``manifold ->
symmetry`` edge closed a cycle and the group could be referenced here only
under ``TYPE_CHECKING``, ten of its members duck-typed (`[M]` by AST on a
``group`` / ``.by`` receiver, 24 sites — the count has been published as three,
seven and nine under narrower predicates).  From R2 the kernel lived in an
``invariance`` module that ``measure`` imported and that imported ``measure``
back; since S5 of the numerics untangling (2026-10-19) invariance is the
measure's own composition of its pushforward with its support's
:meth:`GSpace.action`, and no such module exists.  This module imports
nothing from the measure layer, at any scope: it says which point set a
measure lives on, and the measure a degree of exactness on an orbit space is
AGAINST (the pushforward of the base's Lebesgue measure) is the measure
layer's answer,
:func:`~orpheus.numerics.generating_measure.pushforward_reference`.

Why a closed sum rather than a polymorphic hierarchy: the members are
stable (about a dozen, and every new orbit space is another *instance* of
:class:`Quotient`, not another class) while the operations keep growing.
The deciding evidence is that all of the derivation fields a catalogued
quotient carries belong to :class:`Quotient` alone — on a hierarchy they
would sit on the base returning ``None`` for every non-quotient, which is
the tax :attr:`SubgroupOfO3.mirror_axis` already pays and which
``directional.py`` already branches on.

One intermediate layer is real and is typed: :class:`GSpace`, the members a
rigid motion acts on (:class:`Sphere`, :class:`Ball`, :class:`Circle`,
:class:`Interval`, :class:`RealSpace`, :class:`Quotient`), carrying
:meth:`GSpace.action` by one law, :math:`g` acts iff :math:`g(M) = M`
(rulings 43 and 47 of the numerics untangling, 2026-10-19).  It is the same
rule the base states for quotients — an operation lives on the families that
answer it — applied to the action: an energy-group index set has no action,
and asking it for one is a type error rather than a zero-padded answer.

**The members live in two modules** (layer DAG WP-B, 2026-10-22, ruling 8 of
``.claude/plans/layer_dag.md``).  The base and every geometric member are this
module, in geometry: a point set and the rigid motions acting on it are
statements about the domain, which numerics (measures, quadratures, bases)
imports.  The two finite index sets, :class:`~orpheus.numerics.index_set.IndexSet`
and :class:`~orpheus.numerics.index_set.EnergyGroups`, are the support of a
counting measure and stay in numerics (:mod:`orpheus.numerics.index_set`), which
subclasses this base across the boundary, the allowed direction.  This module
names neither: each member answers :attr:`Manifold.ambient_dim` itself, where
until WP-B a module-level ``match`` named every member, the index sets included.
Until 2026-10-22 this module was ``orpheus/numerics/manifold.py``.
"""

from __future__ import annotations

import functools
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from collections.abc import Callable
from typing import TYPE_CHECKING, Any

import numpy as np
from numpy.typing import ArrayLike, NDArray

from orpheus.geometry.faces import AXIS_INDEX, AXIS_NAMES
from orpheus.geometry.symmetry import SubgroupOfO3

if TYPE_CHECKING:  # RigidMotion is annotation-only here (transformation
    # imports only orpheus.values and symmetry already loads it, so a runtime
    # import would close no cycle: the guard is a choice)
    from orpheus.geometry.transformation import RigidMotion

__all__ = [
    "Manifold",
    "GSpace",
    "Refusal",
    "OutOfDimension",
    "DoesNotMap",
    "DoesNotNormalise",
    "Ball",
    "FundamentalDomain",
    "Sphere",
    "Circle",
    "Interval",
    "RealSpace",
    "Product",
    "Point",
    "Quotient",
    "SPHERE",
    "CIRCLE",
    "UNIT_INTERVAL",
    "COSINE_INTERVAL",
    "HALF_LINE",
    "REAL_LINE",
    "POINT",
    "ManifoldMap",
    "archimedes",
    "barycentre",
    "quotient_onto",
]

#: Tolerance for :meth:`Manifold.contains` on a curved manifold.  A node
#: set is admitted when it satisfies the defining equation to this
#: absolute tolerance.  Chosen to match the construction quality of the
#: shipped quadrature rules, whose nodes are exact to a few ULP; it is
#: NOT a physics tolerance and no caller should widen it to make a
#: measure fit.
_MEMBERSHIP_ATOL = 1e-12


@dataclass(frozen=True)
class Manifold(ABC):
    r"""The point set a measure is supported on — the closed sum's base.

    Only the operations **every** manifold answers live here.  Anything a
    particular family alone can answer belongs on that family, so that
    asking a sphere for a quotient's syzygy ideal is a type error rather
    than a ``None``.

    Subclassing outside this module and :mod:`orpheus.numerics.index_set`
    (the finite index sets) is not supported: the set is closed
    so that a new operation can be checked against every member.
    """

    @property
    @abstractmethod
    def dim(self) -> int:
        """Topological dimension of the manifold."""
        ...

    @property
    @abstractmethod
    def ambient_dim(self) -> int:
        r"""Ambient coordinate count: how many columns a point of this manifold occupies.

        The *embedding* width, distinct from :attr:`dim` (the *topological*
        dimension): :class:`Sphere` is a 2-manifold carried in 3 columns, and
        :class:`~orpheus.numerics.index_set.EnergyGroups` is 0-dimensional
        carried in 1. Both numbers are real and neither derives from the other,
        which is why this is a separate question rather than an arithmetic on
        :attr:`dim`.

        Abstract, so every member answers it: a member that forgets it cannot be
        instantiated. The laws: additive on :class:`Product`, the base's on
        :class:`FundamentalDomain`, the realization's on :class:`Quotient`.
        Until layer DAG WP-B (2026-10-22) this was a module-level ``match``
        naming every member, the index sets included; that match was the one
        place the geometric members named the numerics ones, and the property
        dissolved it (ruling 8 of ``.claude/plans/layer_dag.md``).
        """
        ...

    @property
    @abstractmethod
    def name(self) -> str:
        """Display name, used in messages and in space names.

        The names reproduce the retired ``SPACE_*`` string tags exactly,
        so the migration off ``Space = str`` does not silently re-word
        any message a test pins.
        """
        ...

    @abstractmethod
    def contains(self, points: ArrayLike) -> bool:
        r"""``True`` iff every row of ``points`` is a point OF this manifold.

        This is the predicate that makes a support claim falsifiable.  A
        measure whose nodes fail it is not a measure on this manifold,
        whatever its tag says.

        Parameters
        ----------
        points
            ``(n, d)`` array of candidate points, or ``(n,)`` for a
            1-dimensional ambient space.

        Notes
        -----
        The check is on the AMBIENT coordinates, so it catches the
        ERR-080 class directly: a 1-D quadrature's ordinates
        :math:`(\mu, 0, 0)` satisfy :math:`\|\Omega\| = 1` only at
        :math:`\mu = \pm 1`, so :class:`Sphere` refuses them.
        """
        ...

    def __mul__(self, other: object) -> "Product":
        r"""The product manifold :math:`M \times N`.

        Replaces ``new_space = f"{self.support} x {other.support}"``.
        """
        if not isinstance(other, Manifold):
            return NotImplemented
        return Product(self, other)

    def coordinates(self, points: ArrayLike) -> NDArray:
        """The point identity: two points of this manifold are one point iff
        these coordinates agree.

        The points themselves, except on a :class:`Quotient`, where two
        representatives of one orbit are one point.  On every manifold, not
        only a :class:`GSpace`: relabelling two measures needs no group (the
        user's ruling of 2026-10-19, amending ruling 43)."""
        return np.asarray(points, dtype=float)

    def _as_points(self, points: ArrayLike, ambient: int) -> NDArray:
        """Coerce ``points`` to ``(n, ambient)``, or raise saying why not."""
        arr = np.atleast_1d(np.asarray(points, dtype=float))
        if ambient == 1 and arr.ndim == 1:
            arr = arr[:, None]
        if arr.ndim != 2 or arr.shape[1] != ambient:
            raise ValueError(
                f"{self.name}.contains expects points of ambient dimension "
                f"{ambient} — an (n, {ambient}) array"
                f"{' or an (n,) array' if ambient == 1 else ''}; got shape "
                f"{arr.shape}."
            )
        return arr


@dataclass(frozen=True)
class GSpace(Manifold):
    r"""A manifold the rigid motions of its ambient space act on — a
    :math:`G`-space, :math:`G \le E(d)`.

    An invariance question ("does :math:`g` permute these atoms?") is asked of
    a measure ON its support, so the support must say how :math:`g` moves its
    points.  The families that can say it are the members of this class:
    :class:`Sphere`, :class:`Ball`, :class:`Circle`, :class:`Interval`,
    :class:`RealSpace` and :class:`Quotient`.
    :class:`~orpheus.numerics.index_set.EnergyGroups`,
    :class:`~orpheus.numerics.index_set.IndexSet`, :class:`Product` and :class:`FundamentalDomain` are not
    G-spaces and carry no :meth:`action` at all — asking a group index set how a
    reflection moves it is a type error, not a ``None`` (the base class's rule;
    ruling 43 of the numerics untangling, 2026-10-19).  Until then
    ``invariance._orbit_space_of`` asked EVERY support on
    :math:`\mathbb{R}^3/\{e\}` after zero-padding its nodes, so an
    energy-group measure answered ``SO2_x``-invariant: True, and the chart
    :math:`\mu`-rule answered ``SO2_x`` True and ``SO2_y`` False, the "x" being
    the padding's column 0 and not a property of the rule.

    **One law** (ruling 47, 2026-10-19): a motion :math:`g` acts on :math:`M`
    iff :math:`g(M) = M`, and the action is then :math:`x \mapsto g(x)`,
    :math:`M \to M`.  :meth:`action` states it once; each family answers only
    whether :math:`g` preserves the data that DEFINES it (:meth:`_refusal`):

    =====================  ==================================  =========================
    family                 defining data                       :math:`g(M) = M` iff
    =====================  ==================================  =========================
    Sphere, Ball, Circle   centred at the origin               :math:`g` fixes the origin
    RealSpace(d)           all of :math:`\mathbb{R}^d`          always (:math:`E(d)`)
    Interval [a, b]        the endpoint set :math:`\{a, b\}`   :math:`g` permutes it
    Quotient M/H           the base and :math:`H`              the base's law, and
                                                               :math:`g H g^{-1} = H`
    =====================  ==================================  =========================

    There is no special case: the polar mirror :math:`\mu \to -\mu` of
    ``COSINE_INTERVAL`` is the Interval law's answer at :math:`a + b = 0`
    (:math:`x \mapsto a + b - x`), and on :math:`[0, 1]` the same law admits
    :math:`x \mapsto 1 - x`.  What is special about :math:`[-1, 1]` is that it
    is the chart of :math:`S^2/O(2)_x`, where :math:`\mu \to -\mu` is the
    action :math:`\sigma_x` induces on the quotient — a fact a gate checks
    (the polar mirror spelled twice must agree), never a branch.

    Three readings of a point accompany the action:

    * :attr:`acting_dimension` — the :math:`d` of the motions in :math:`E(d)`
      that are asked to act: the base's ambient width on a quotient (its
      points are acted on through the base), the manifold's own elsewhere;
    * :meth:`~Manifold.coordinates` — the point identity, inherited from
      :class:`Manifold` (every manifold has one; a quotient overrides it);
    * :meth:`ambient_points` — where a point sits in the space the motions act
      on: the orbit barycentre on a quotient, the point itself elsewhere
      (ruling 46).
    """

    def quotient(self, group: SubgroupOfO3) -> "Quotient":
        r"""The orbit space :math:`M/H`, by catalogue lookup.

        On a G-space only (S5, 2026-10-19): an orbit space is the set of
        orbits of a group ACTING on the manifold, and a manifold no rigid
        motion acts on has none — until then this lived on :class:`Manifold`
        and ``ENERGY.quotient(Trivial)`` built an entry.

        Computing an orbit space from scratch — invariant ring, syzygy
        ideal by elimination, then the Procesi–Schwarz PSD condition — is
        mechanical but is a symbolic-algebra engine, and the groups that
        occur in transport number about a dozen.  So entries are derived
        once by hand and recorded; the engine is **deferred, not
        refused**, and when it is built it populates exactly these
        fields rather than introducing a parallel representation.

        Raises
        ------
        NotImplementedError
            When no entry exists, naming the missing WORK rather than
            the gap, so a fresh session can pick it up.
        ValueError
            When ``group`` is not the orbit space's STABILISER
            (:attr:`~orpheus.geometry.symmetry.SubgroupOfO3.orbit_stabiliser`)
            — the largest subgroup with its orbits.  A smaller orbit-equivalent
            group names the same point set, so it is not a second entry
            but a second spelling, and it is refused with the theorem
            (`[M]` ``SPHERE.quotient(SubgroupOfO3.SO2("x"))`` is refused
            naming ``O2("x")``; #432, 2026-09-02).  And when ``self`` is
            already an orbit space and ``group`` lies in the group it spent:
            such a group acts TRIVIALLY here, the quotient is the identity,
            and ``S^2/sigma_y/sigma_y`` would be a second spelling of
            ``S^2/sigma_y`` (#429 tracker 2.2b) — refused naming the orbit
            space the caller already has.

        Notes
        -----
        Memoised on ``(self, group)``: an orbit space is derived ONCE
        and recorded, which is the catalogue's own philosophy, and since
        2026-09-01 every slab quadrature carries one (tracker 2.4), so a
        symbolic derivation per rule construction would put ~6 ms of
        SymPy on the path of every slab solve. Every :class:`Manifold`
        is a frozen value type and hashable, which is what the memo
        requires. The returned object is shared, which is safe for the
        same reason.
        """
        return _catalogued_quotient(self, group)

    @property
    def acting_dimension(self) -> int:
        """The :math:`d` of the rigid motions in :math:`E(d)` asked to act here."""
        return self.ambient_dim

    @abstractmethod
    def _refusal(self, motion: "RigidMotion") -> "DoesNotMap | DoesNotNormalise | None":
        """The family's own law: why ``motion`` does not map this manifold onto
        itself, or ``None``.

        Called by :meth:`refusal` with a motion of :attr:`acting_dimension`
        only.  Every family answers :class:`DoesNotMap`; a quotient also
        answers :class:`DoesNotNormalise`.
        """
        ...

    def refusal(self, motion: "RigidMotion") -> "Refusal | None":
        r"""Why ``motion`` does not act on this manifold, or ``None`` if it does
        — the one door to the law, as a VALUE.

        A typed value rather than an exception, so a caller that turns "this
        motion is no symmetry here" into ``None`` branches on the KIND and
        cannot swallow an unrelated error (the S5 elegance review, 2026-10-19:
        an ``except ValueError`` around :meth:`action` read an internal bug as
        "not a symmetry").  The kinds are a closed sum:
        :class:`OutOfDimension`, :class:`DoesNotMap`, :class:`DoesNotNormalise`.
        """
        if motion.dimension != self.acting_dimension:
            return OutOfDimension(
                f"{motion!r} does not act on {self.name}: it is a motion of "
                f"R^{motion.dimension}, and {self.name} is acted on by the "
                f"motions of R^{self.acting_dimension}."
            )
        return self._refusal(motion)

    def _image_map(self, motion: "RigidMotion") -> Callable[[NDArray], NDArray]:
        """The action's ``apply``: the motion on the points, shape kept."""
        return functools.partial(_apply_motion, motion)

    def action(self, motion: "RigidMotion") -> ManifoldMap:
        r"""The action of ``motion`` on this manifold, :math:`x \mapsto g(x)`,
        as the typed arrow :math:`M \to M`.

        Raises
        ------
        ValueError
            With the :meth:`refusal`'s message: ``does not act`` (the motion is
            not one of :math:`E(d)` for this manifold's
            :attr:`acting_dimension`), ``does not map`` (it does not carry the
            manifold onto itself), ``normalise`` (on a quotient, it does not
            normalise the group spent there).
        """
        refusal = self.refusal(motion)
        if refusal is not None:
            raise ValueError(refusal.message)
        return ManifoldMap(domain=self, codomain=self, apply=self._image_map(motion))

    def ambient_points(self, points: ArrayLike) -> NDArray:
        """Each point's position in the space the motions act on.

        The points themselves, except on a :class:`Quotient`, where it is the
        orbit barycentre."""
        return np.asarray(points, dtype=float)


def _apply_motion(motion: "RigidMotion", points: NDArray) -> NDArray:
    """``motion`` on ``points``, keeping an ``(n,)`` interval array ``(n,)``.

    Module-level, partially applied, so the arrow is picklable (as
    :func:`_act_through` is)."""
    arr = np.asarray(points, dtype=float)
    if arr.ndim == 1:
        return motion.on_points(arr[:, None])[:, 0]
    return motion.on_points(arr)


@dataclass(frozen=True)
class Refusal:
    """Why a rigid motion does not act on a :class:`GSpace` — a closed sum of
    three kinds, each carrying the message :meth:`GSpace.action` raises."""

    message: str


@dataclass(frozen=True)
class OutOfDimension(Refusal):
    """The motion is not one of :math:`E(d)` for the space's acting dimension."""


@dataclass(frozen=True)
class DoesNotMap(Refusal):
    r"""The motion does not carry the manifold onto itself: :math:`g(M) \ne M`."""


@dataclass(frozen=True)
class DoesNotNormalise(Refusal):
    r"""On an orbit space :math:`M/H`: the motion does not normalise :math:`H`,
    so :math:`[p] \mapsto [g p]` is not a function of the orbit."""


def _origin_refusal(manifold: Manifold, motion: "RigidMotion") -> DoesNotMap | None:
    """The centred families' law: :math:`g(M) = M` iff :math:`g` fixes the origin."""
    if np.all(np.abs(motion.translation) <= _MEMBERSHIP_ATOL):
        return None
    return DoesNotMap(
        f"{motion!r} does not map {manifold.name} onto itself: it moves the "
        f"centre to {np.asarray(motion.translation).tolist()}. Only the motions "
        f"fixing the origin, O({motion.dimension}), carry {manifold.name} onto itself."
    )


# ---------------------------------------------------------------------------
# Curved members — the defining equation IS the membership predicate
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Sphere(GSpace):
    r""":math:`S^2 = \{\Omega \in \mathbb{R}^3 : \|\Omega\| = 1\}`."""

    @property
    def dim(self) -> int:
        return 2

    @property
    def ambient_dim(self) -> int:
        return 3

    @property
    def name(self) -> str:
        return "S^2"

    def contains(self, points: ArrayLike) -> bool:
        arr = self._as_points(points, 3)
        return bool(
            np.all(
                np.abs(np.linalg.norm(arr, axis=1) - 1.0) <= _MEMBERSHIP_ATOL
            )
        )


    def _refusal(self, motion: "RigidMotion") -> DoesNotMap | None:
        return _origin_refusal(self, motion)


@dataclass(frozen=True)
class Circle(GSpace):
    r""":math:`S^1 = \{u \in \mathbb{R}^2 : \|u\| = 1\}`.

    The circle names the MANIFOLD, not a chart of it — the tag was
    ``"[0,2pi)"`` until 2026-08-02, which asserted a coordinate.  Its
    nodes are the roots of unity as POINTS, because only in that
    representation is an on-axis node exactly ``0.0``.
    """

    @property
    def dim(self) -> int:
        return 1

    @property
    def ambient_dim(self) -> int:
        return 2

    @property
    def name(self) -> str:
        return "S^1"

    def contains(self, points: ArrayLike) -> bool:
        arr = self._as_points(points, 2)
        return bool(
            np.all(
                np.abs(np.linalg.norm(arr, axis=1) - 1.0) <= _MEMBERSHIP_ATOL
            )
        )


    def _refusal(self, motion: "RigidMotion") -> DoesNotMap | None:
        return _origin_refusal(self, motion)


# ---------------------------------------------------------------------------
# Flat members
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Interval(GSpace):
    r"""The real interval :math:`[a, b]`, endpoints included.

    One family, not three tags: the retired ``SPACE_INTERVAL_M11``,
    ``SPACE_INTERVAL_01``, ``SPACE_HALF_LINE`` and ``SPACE_R`` are four
    MEMBERS of it — which is what ``generating_measure.py``'s
    ``support=f"[{a},{b}]"`` was already saying, one interpolation at a
    time.
    """

    a: float
    b: float

    def __post_init__(self) -> None:
        if not self.a < self.b:
            raise ValueError(
                f"Interval requires a < b; got a={self.a!r}, b={self.b!r}. "
                f"A degenerate interval is a point, not a 1-manifold."
            )

    @property
    def dim(self) -> int:
        return 1

    @property
    def ambient_dim(self) -> int:
        return 1

    @property
    def name(self) -> str:
        # Reproduce the retired tags verbatim so no pinned message moves.
        if self.a == -np.inf and self.b == np.inf:
            return "R"
        if self.a == 0.0 and self.b == np.inf:
            return "[0,inf)"
        return f"[{_fmt(self.a)},{_fmt(self.b)}]"

    def contains(self, points: ArrayLike) -> bool:
        arr = self._as_points(points, 1)[:, 0]
        # Finiteness is a SEPARATE condition from the bounds: on an
        # unbounded interval `inf <= inf` holds, and infinity is not a
        # point of R.  A NaN fails every comparison and so is refused
        # by the bounds anyway; asserting finiteness says why.
        return bool(
            np.all(np.isfinite(arr))
            and np.all(arr >= self.a - _MEMBERSHIP_ATOL)
            and np.all(arr <= self.b + _MEMBERSHIP_ATOL)
        )


    def _refusal(self, motion: "RigidMotion") -> DoesNotMap | None:
        r""":math:`g([a, b]) = [a, b]` iff :math:`g` permutes the endpoint set
        :math:`\{a, b\}`: the identity always, :math:`x \mapsto a + b - x`
        on a bounded interval, every motion on :math:`\mathbb{R}`, only the
        identity on a half-line.  The endpoints may be infinite; ``g(±inf)``
        is read off the affine map, as the image interval's own ends."""
        q = float(motion.linear[0, 0])
        t = float(motion.translation[0])
        lo, hi = sorted((q * self.a + t, q * self.b + t))
        if np.isclose(lo, self.a, rtol=0.0, atol=_MEMBERSHIP_ATOL) and np.isclose(
            hi, self.b, rtol=0.0, atol=_MEMBERSHIP_ATOL
        ):
            return None
        return DoesNotMap(
            f"{motion!r} does not map {self.name} onto itself: it carries it "
            f"onto [{_fmt(lo)},{_fmt(hi)}]. A motion acts on an interval iff it "
            f"permutes the endpoints {{{_fmt(self.a)}, {_fmt(self.b)}}}."
        )


@dataclass(frozen=True)
class RealSpace(GSpace):
    r""":math:`\mathbb{R}^d` — the spatial manifold, indexed by ``d``.

    Replaces the interpolated ``support=f"spatial_R{ndim}"``.
    """

    d: int

    def __post_init__(self) -> None:
        if self.d < 1:
            raise ValueError(
                f"RealSpace requires d >= 1; got d={self.d!r}."
            )

    @property
    def dim(self) -> int:
        return self.d

    @property
    def ambient_dim(self) -> int:
        return self.d

    @property
    def name(self) -> str:
        return f"spatial_R{self.d}"

    def contains(self, points: ArrayLike) -> bool:
        arr = self._as_points(points, self.d)
        return bool(np.all(np.isfinite(arr)))


    def _refusal(self, motion: "RigidMotion") -> DoesNotMap | None:
        r"""Every motion of :math:`E(d)` maps :math:`\mathbb{R}^d` onto itself."""
        return None


# ---------------------------------------------------------------------------
# Recursive members
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Product(Manifold):
    r""":math:`M \times N`, with coordinates concatenated in that order."""

    left: Manifold
    right: Manifold

    @property
    def dim(self) -> int:
        return self.left.dim + self.right.dim

    @property
    def ambient_dim(self) -> int:
        """Additive: the factors' columns, concatenated."""
        return self.left.ambient_dim + self.right.ambient_dim

    @property
    def name(self) -> str:
        return f"{self.left.name} × {self.right.name}"

    def contains(self, points: ArrayLike) -> bool:
        arr = np.atleast_2d(np.asarray(points, dtype=float))
        split = self.left.ambient_dim
        if arr.shape[1] != split + self.right.ambient_dim:
            raise ValueError(
                f"{self.name}.contains expects ambient dimension "
                f"{split + self.right.ambient_dim}; got shape {arr.shape}."
            )
        return self.left.contains(arr[:, :split]) and self.right.contains(
            arr[:, split:]
        )


@dataclass(frozen=True)
class Ball(GSpace):
    r"""The closed unit ball :math:`D^d = \{p \in \mathbb{R}^d : \|p\| \le 1\}`.

    Minted 2026-08-31 because the derivation demanded it: ``S^2/sigma_y``
    IS the closed 2-disk in invariant coordinates, and no shipped member
    could say so.  ``Product(COSINE_INTERVAL, COSINE_INTERVAL)`` is the
    bounding SQUARE, and the discriminator is measured — ``(0.9, 0.9)``
    is in the square and not in the disk, and it corresponds to NO
    direction, since ``eta^2 + mu^2 = 1.62 > 1`` forces ``xi^2 < 0``.
    """

    d: int

    @property
    def dim(self) -> int:
        return self.d

    @property
    def ambient_dim(self) -> int:
        return self.d

    @property
    def name(self) -> str:
        return f"D^{self.d}"

    def contains(self, points: ArrayLike) -> bool:
        arr = self._as_points(points, self.d)
        return bool(
            np.all(np.isfinite(arr))
            and np.all(
                np.sum(arr * arr, axis=1) <= 1.0 + _MEMBERSHIP_ATOL
            )
        )


    def _refusal(self, motion: "RigidMotion") -> DoesNotMap | None:
        return _origin_refusal(self, motion)


@dataclass(frozen=True)
class Point(Manifold):
    r"""The one-point manifold: dimension 0, carried in 0 coordinates.

    The orbit space of a group acting TRANSITIVELY is a point: :math:`O(3)`
    carries any direction onto any other, so :math:`S^2/O(3)` has one orbit
    (:func:`_sphere_mod_o3`). It is the angular domain of a method that
    integrates direction away analytically: collision probability's
    isotropic angular domain is this point (#440).

    Zero coordinates, not one: a chart of an orbit space reads the
    invariants that SURVIVE the base's own ideal, and on :math:`S^2` the one
    :math:`O(3)` invariant, :math:`x^2 + y^2 + z^2`, is the constant 1, so
    none survives. A point is therefore an ``(n, 0)`` array, and every such
    row is the point.

    Not a finite index set (:class:`~orpheus.numerics.index_set.IndexSet`):
    that is a labelled set of degree-of-freedom indices, carried in one
    column, with no geometry; this is a point of space.
    """

    @property
    def dim(self) -> int:
        return 0

    @property
    def ambient_dim(self) -> int:
        return 0

    @property
    def name(self) -> str:
        return "point"

    def contains(self, points: ArrayLike) -> bool:
        self._as_points(points, 0)
        return True


@dataclass(frozen=True)
class FundamentalDomain(Manifold):
    r"""A strict fundamental domain of :math:`H` acting on ``base`` — the
    IMAGE of a section of :math:`M \to M/H`, in the BASE's coordinates.

    Cut from ``base`` by closed half-spaces
    :math:`\langle p, n_i \rangle \ge 0`.  One normal cuts a half-space;
    an antipodal PAIR :math:`\{n, -n\}` cuts the hyperplane
    :math:`\langle p, n \rangle = 0`, so equalities need no second slot —
    which is why one tuple expresses both the :math:`\sigma_y` hemisphere
    and the :math:`SO(2)` half-meridian.

    ⭐ **This is the OTHER half of a quotient, and the tree has always
    needed it.**  :meth:`~orpheus.numerics.measure.DiscreteMeasure.quotient`
    computes orbit representatives and then keeps ``nodes[representative]``
    — a SELECTION, applying no chart — so every measure it emits carries
    the base's ambient columns.  A chart codomain
    (:attr:`Quotient.realization`) cannot validate those points; only a
    fundamental domain can.

    ⛔ **The inequalities must be CLOSED.**  ``[M]`` the cylindrical march
    seeds a level at :math:`\xi = 0` exactly — on the stratum — so a strict
    :math:`\langle p, n \rangle > 0` would refuse a direction production
    marches from (`coding-elegance` anti-pattern #18's half (ii): every
    legal value must be admitted, which is a claim about the PRODUCERS).
    """

    base: Manifold
    #: Outward normals of the closed half-spaces, in the base's ambient
    #: coordinates.  An antipodal pair spells an equality.
    normals: tuple[tuple[float, ...], ...]
    #: Short label for :attr:`name` (e.g. ``"xi>=0"``), since a normal
    #: tuple makes an unreadable one.
    label: str

    @property
    def dim(self) -> int:
        """``base.dim`` less one per antipodal pair (each is an equality).

        An inequality carves a region WITH BOUNDARY and does not drop
        dimension; an equality does.  This is what lets the σ_y
        hemisphere read 2 and the SO(2) half-meridian read 1 from one
        rule — and :class:`Quotient` gates it against the chart.
        """
        seen = {tuple(n) for n in self.normals}
        pairs = sum(
            1
            for n in seen
            if tuple(-c for c in n) in seen and n < tuple(-c for c in n)
        )
        return self.base.dim - pairs

    @property
    def ambient_dim(self) -> int:
        """Invariant: a fundamental domain is cut from its base's columns."""
        return self.base.ambient_dim

    @property
    def name(self) -> str:
        return f"{self.base.name}|{self.label}"

    def contains(self, points: ArrayLike) -> bool:
        arr = self._as_points(points, self.base.ambient_dim)
        if not self.base.contains(arr):
            return False
        # Closed half-spaces — see the class docstring's ⛔ note on why
        # this may not be strict.
        return bool(
            all(
                np.all(arr @ np.asarray(n, dtype=float) >= -_MEMBERSHIP_ATOL)
                for n in self.normals
            )
        )


@dataclass(frozen=True)
class Quotient(GSpace):
    r"""The orbit space :math:`M/H`, and the derivation that produced it.

    ⭐ **The fields below ARE the derivation procedure's output**, not a
    human summary of its answer.  That is what makes this class the SEED
    of the deferred engine rather than a rival to it: an engine ships
    when it computes these same fields, introducing no new type and no
    translation layer.  The ~12 catalogued entries' symbolic regression
    tests are therefore the engine's specification, written before it.

    ⚠ An orbit space is generally an **orbifold**, not a quotient
    manifold: where the action is not free, the image of :math:`H`'s
    fixed-point set is a singular stratum.  Any design that assumes a
    quotient is a smooth submersion is wrong exactly there — and for
    :math:`S^2/O(2)_a` that locus is :math:`\mu = \pm 1`, the poles.
    """

    base: GSpace
    #: The group the orbit space is NAMED BY — its **stabiliser**: the
    #: largest subgroup of :math:`O(3)` whose orbits are these orbits,
    #: equivalently the largest one fixing every invariant in
    #: :attr:`generators`.  A smaller orbit-equivalent group gives the same
    #: point set and the same derivation output, so it is not a second
    #: entry — it is refused (`[M]` :math:`S^2/SO(2)_a` is
    #: :math:`S^2/O(2)_a`: a vertical mirror carries each constant-\ :math:`\mu`
    #: circle onto itself), and refused WHERE IT IS WRITTEN: ``__post_init__``
    #: reads ``by.orbit_stabiliser``, so a ``replace(entry, by=...)`` cannot
    #: forge a second spelling either, and the catalogue door reads it before
    #: any lookup so the answer is the theorem and never "no catalogue
    #: entry" (#432, 2026-09-02).  One orbit space, one spelling — which is what lets
    #: :attr:`~orpheus.numerics.basis.base.Basis.invariance_group`, read
    #: off this field, be the FULL group a basis on the entry has.
    by: SubgroupOfO3
    #: The orbit space realized as an honest manifold — what a chart of
    #: ``M/H`` maps ONTO.  For ``S^2/O(2)_a`` this is ``Interval(-1, 1)``.
    #: Its coordinates are the INVARIANTS', the same language as
    #: :attr:`generators`, :attr:`gram` and :attr:`det_gram`.
    realization: Manifold
    #: The quotient map's action on the base's AMBIENT coordinates,
    #: returning the orbit's coordinates in the :attr:`realization`'s
    #: language — D0.1's "chart" output, as the derivation returns it:
    #: the invariants that SURVIVE eliminating the base's own ideal
    #: (``p_1 = x_a`` for ``S^2/O(2)_a``; ``(p_1, p_2) = (x_b, x_c)`` for
    #: ``S^2/sigma_a``; every coordinate for ``M/{e}``).  An engine would
    #: ``lambdify`` them; the hand entries spell the column selection
    #: directly, and the entry's regression test pins the two against each
    #: other.  The typed arrow is :attr:`quotient_map`; this is only its
    #: ``apply``, kept as a plain callable because the arrow's codomain is
    #: the entry ITSELF, which does not exist until construction ends.
    #: Excluded from equality and hashing: a function has no value
    #: equality, and the map is DERIVED from ``(base, by)``, which the
    #: comparison already covers.
    orbit_coordinates: Callable[[NDArray], NDArray] = field(
        compare=False, repr=False
    )
    #: The lift's action on CHART coordinates — a right inverse of the quotient
    #: map, as the derivation returns it: the orbit BARYCENTRE, which is the
    #: Reynolds projector onto :math:`H`'s fixed subspace read from the chart's
    #: side (the axial entry: :math:`\mu \mapsto \mu\,\hat e_a`; a mirror
    #: entry: :math:`(x_b, x_c) \mapsto (0, x_b, x_c)`; the trivial entry: the
    #: identity).  A field for the reason :attr:`orbit_coordinates` is one — it
    #: is derivation OUTPUT — and until 2026-09-03 (#434 R4) it was a three-arm
    #: branch on the group's tag whose own message read "add it to
    #: ``Quotient.lift``": a second dispatch over the key the catalogue had
    #: already used to choose the builder, and one a seventh entry could forget.
    #: Excluded from equality and hashing as :attr:`orbit_coordinates` is.
    lift_coordinates: Callable[[NDArray], NDArray] = field(
        compare=False, repr=False
    )
    #: What the lift lands ON — the ball :math:`D^3` for a sphere entry, whose
    #: barycentres leave the sphere (on it only where the orbit is a point:
    #: the axial poles, the mirror's equator), the base itself when the lift is
    #: the identity.  The codomain a consumer reads to learn which it was
    #: handed.  COMPARED, unlike :attr:`orbit_coordinates` and
    #: :attr:`lift_coordinates`: those are excluded because a FUNCTION has no
    #: value equality, and a :class:`Manifold` has one.  A second entry that
    #: agrees on ``(base, by)`` and disagrees on where its lift lands is two
    #: orbit spaces, not one — `[M]` 2026-09-03 (elegance review) with the
    #: field excluded, ``replace(entry, lift_codomain=SPHERE)`` compared EQUAL
    #: to the catalogue entry and, :func:`barycentre` being memoised on the
    #: entry, whichever was asked FIRST answered for both: ERR-080's shape
    #: (an arrow off the sphere declaring itself onto it), re-minted by the
    #: field built to refuse it.  Gated in ``__post_init__`` against the base's
    #: ambient width, where the lift actually lands.
    lift_codomain: Manifold
    #: The section's IMAGE, in the BASE's coordinates — the other of the
    #: quotient's two honest coordinate systems (user ruling,
    #: 2026-08-31).  ``None`` when no canonical section exists, which for
    #: a positive-dimensional group is the normal case: any half-meridian
    #: sections ``S^2 -> S^2/O(2)_a`` and none is distinguished.
    #:
    #: ⭐ Why a SECOND slot rather than a wider ``realization``: the two
    #: answer different questions and the tree needs both.  ``[M]`` the
    #: chart is Mode-12 BLIND to the ERR-080 forgery — ``(x,y,z) -> (x,z)``
    #: drops exactly the coordinate the forgery corrupts, so the forged
    #: row ``(mu, 0)`` is a legal point of the disk — while the section
    #: refuses it.  And ``[M]`` ERR-080 IS a botched section of
    #: ``S^2/O(2)_a``: a consumer needed one, the realization is a chart,
    #: and the tree fabricated one by zero-padding to ``(mu, 0, 0)``,
    #: which is off ``S^2``.  With this slot that padding has nowhere to
    #: live.
    fundamental_domain: Manifold | None = None
    #: Minimal generators of the invariant ring, as SymPy expressions in
    #: the ambient coordinates.  ``()`` is a legal value only for the
    #: trivial group.
    #:
    #: The four SymPy fields below and :attr:`singular_stratum` are excluded
    #: from equality, hashing and the content digest: they are the
    #: derivation's OUTPUT from ``(base, by)``, which the comparison already
    #: covers, and a SymPy expression has no content encoding, so while they
    #: were compared an axis on an orbit space could not be keyed at all
    #: (S7 of the numerics untangling, ruling 71: the support is axis
    #: content). The theorems in ``__post_init__`` still read them.
    generators: tuple[Any, ...] = field(default=(), compare=False)
    #: Generators of the syzygy ideal ``I = ker(R[y] -> R[x])``.  Empty
    #: when the invariants are algebraically independent.
    syzygy: tuple[Any, ...] = field(default=(), compare=False)
    #: ``P_ij = <grad p_i, grad p_j>`` re-expressed in the invariants —
    #: the Procesi–Schwarz matrix whose ``P >= 0`` locus cuts out the
    #: orbit space's inequalities.
    gram: Any | None = field(default=None, compare=False)
    #: ``det P``.  Its vanishing locus IS the orbit-space boundary, and
    #: it is the squared orbit radius.
    det_gram: Any | None = field(default=None, compare=False)
    #: ``"hand"`` for an entry derived once by a human, ``"engine"`` once
    #: the derivation engine populates it.  A MIXED state is expressible
    #: on purpose: an incremental engine rollout is exactly that, and
    #: without this field the migration would have to be all-or-nothing.
    derived_by: str = "hand"

    def __post_init__(self) -> None:
        r"""Three theorems, refused where the entry is WRITTEN, not where it is read.

        1. **The group is the orbit space's STABILISER**
           (:func:`_assert_named_by_stabiliser`).  `[M]` 2026-09-02, before
           this check: ``replace(SPHERE.quotient(O2("x")), by=SO2("x"))``
           constructed, compared unequal to the catalogue entry, and was
           accepted by :func:`barycentre` and the polar-axis reader —
           ERR-080's own defect class (one orbit space, two objects) reopened
           by the very idiom Pattern 4 blesses, because the refusal then
           lived in ONE builder rather than on the type.
        2. **The chart REALIZES the orbit space's dimension**, :math:`\dim
           M/H = \dim M - \dim(H \cdot p)` at a generic :math:`p` — the
           generic ORBIT's dimension, asked of the group over the base's
           generic points (:meth:`SubgroupOfO3.generic_orbit_dimension
           <orpheus.geometry.symmetry.SubgroupOfO3.generic_orbit_dimension>`),
           and NOT :math:`\dim H`: :math:`O(3)` on :math:`S^2` has
           :math:`\dim 3` and 2-dimensional orbits (the sphere is one orbit),
           :math:`SO(3)` on :math:`\mathbb{R}^3` likewise, so a law spelled
           ``base.dim - by.dim`` would read −1 and 0.  `[M]` 2026-09-03
           (elegance review), before this check: a forged :math:`S^2/O(2)_z`
           realized on the DISK and a forged :math:`S^2/\sigma_x` realized on
           :math:`[-1,1]` both constructed and compared unequal to the
           catalogue entry — the same defect class, one field over from the
           one #432 closed.  A non-maximal ``by`` passes this clause
           (:math:`SO(2)_a` has :math:`O(2)_a`'s orbit dimension) and is
           refused by clause 1: each clause has its own discriminating input.
        3. **The lift lands in the base's AMBIENT space** —
           :attr:`lift_codomain`'s width is the base's, because
           :func:`_act_through` hands the lift's output straight to
           :meth:`RigidMotion.on_points <orpheus.geometry.transformation.RigidMotion.on_points>`.
        4. **Where a fundamental domain ships, the two coordinate systems
           describe the SAME orbit space** — their dimensions agree, and that
           is a real check, not a tautology: the fundamental domain derives
           its ``dim`` from the base less one per antipodal normal pair, while
           the realization states its own.  `[M]` the σ_y hemisphere reads 2
           against the disk's 2; a half-meridian (an antipodal PAIR) reads 1
           against the disk's 2 and is refused HERE — the discriminating
           input, because clause 2 runs FIRST and refuses its own two
           forgeries before this clause sees them (they ship
           ``fundamental_domain=None`` besides; a σ_x forgery carrying
           the x-hemisphere, ``dim`` 2 against the interval's 1, would
           violate this clause too and still be refused by clause 2).
        """
        if not isinstance(self.base, GSpace):
            raise TypeError(
                f"{self.base.name}/{self.by.name}: the base {self.base.name} is "
                f"not a G-space, so no group acts on it and it has no orbit space."
            )
        _assert_named_by_stabiliser(self.base, self.by)
        # The generic orbit's dimension is a property of the GROUP acting on
        # the base's generic points; a finite group has point orbits, so it
        # is asked of no point at all (some bases have none to offer).
        orbit_dimension = (
            0 if self.by.dim == 0
            else self.by.generic_orbit_dimension(_generic_points(self.base))
        )
        expected = self.base.dim - orbit_dimension
        if self.realization.dim != expected:
            raise ValueError(
                f"{self.name}: dim(M/H) = dim {self.base.name} - dim(generic "
                f"{self.by.name}-orbit) = {expected}, but the realization "
                f"{self.realization.name!r} has dim {self.realization.dim}. A "
                f"chart of an orbit space realizes its dimension; check the "
                f"elimination against {self.base.name}'s own ideal."
            )
        if self.lift_codomain.ambient_dim != self.base.ambient_dim:
            raise ValueError(
                f"{self.name}: the lift lands in {self.base.name}'s ambient "
                f"space ({self.base.ambient_dim} columns), but lift_codomain "
                f"{self.lift_codomain.name!r} has "
                f"{self.lift_codomain.ambient_dim}. A codomain a consumer "
                f"reads to learn what it was handed must be the one the map "
                f"actually lands in."
            )
        if self.fundamental_domain is None:
            return
        if self.fundamental_domain.dim != self.realization.dim:
            raise ValueError(
                f"{self.name}: the fundamental domain "
                f"{self.fundamental_domain.name!r} has dim "
                f"{self.fundamental_domain.dim} but the realization "
                f"{self.realization.name!r} has dim "
                f"{self.realization.dim} — the two must describe the same "
                f"orbit space. Check the normals: an antipodal PAIR spells "
                f"an equality and drops a dimension; a lone normal does not."
            )

    @property
    def dim(self) -> int:
        return self.realization.dim

    @property
    def ambient_dim(self) -> int:
        """The realization's columns: the orbit space is carried in its invariant chart."""
        return self.realization.ambient_dim

    @property
    def name(self) -> str:
        return f"{self.base.name}/{self.by.name}"

    @property
    def quotient_map(self) -> ManifoldMap:
        r"""The quotient map :math:`\pi : M \to M/H`, as a typed arrow.

        Its codomain is THIS entry — never the :attr:`realization`, which
        is where the same numbers land when read as a chart, and which is
        exactly the reading tracker 2.4 made refusable (a rule on
        :math:`[-1,1]` is not a rule on :math:`S^2/O(2)_x`).  A measure
        pushed forward along it lives on the orbit space; a function on
        the orbit space pulled back along it is :math:`H`-invariant by
        construction, which is what :class:`~orpheus.numerics.basis.descent.Descent`
        reads, and what a frame tabulates through (:func:`quotient_onto`).

        ⚠ Its image is in the REALIZATION's coordinates.  The orbit
        retraction that :meth:`DiscreteMeasure.quotient
        <orpheus.numerics.measure.DiscreteMeasure.quotient>` builds lands
        on the same codomain in the SECTION's coordinates (a representative
        per realized orbit, still in the base's ambient width); the two
        are the two honest coordinate systems :meth:`contains` accepts.

        Derived, not stored, for the reason
        :attr:`~orpheus.numerics.basis.base.Basis.invariance_group` is:
        the arrow is a function of fields the entry already carries
        (:attr:`base`, itself, :attr:`orbit_coordinates`), and a stored
        copy would be one more thing able to disagree with them.  `[M]`
        2026-09-02: :math:`\pi_a \circ \varphi_a = \mathrm{pr}_1`
        bit-exactly on 12 of 12 (three axes × four Gauss orders) and
        :math:`\beta_a \circ \pi_a` the axial projection on 3 of 3 —
        :func:`archimedes` and :func:`barycentre` are the arrows on either
        side of this one.
        """
        return ManifoldMap(
            domain=self.base, codomain=self, apply=self.orbit_coordinates
        )

    @property
    def lift(self) -> ManifoldMap:
        r"""The lift of the chart into the base's ambient space — the orbit
        BARYCENTRE, :math:`\lambda : M/H \to D^3`, a right inverse of the
        quotient map (:math:`\pi \circ \lambda = \mathrm{id}` on the
        :attr:`realization`).

        On every catalogued entry the lift is the Reynolds projector
        :math:`P_H = \int_H \rho(g)\,dg` onto :math:`H`'s fixed subspace, read
        from the chart's side: :math:`\mu \mapsto \mu\,\hat e_a` for the
        axial entry (the centre of the constant-:math:`\mu` circle),
        :math:`(x_b, x_c) \mapsto (0, x_b, x_c)` for a mirror entry (the
        midpoint of :math:`\{p, \sigma_a p\}`), the identity for the trivial
        one — ONE formula, :func:`_coordinate_chart`'s ``embed`` half, where
        until 2026-09-03 (#434 R4) the mirror entry lifted through a hemisphere
        section (a square root and a :math:`\rho > 1` refusal) and the axial
        one through a barycentre, chosen by a branch on the group's tag.  The
        barycentre is not a section — it lands ON the base only where the
        orbit is a point — but it is what an ORBIT question needs: canonical,
        hence EQUIVARIANT under every isometry that descends to the orbit
        space (:math:`g P_H g^{-1} = P_{gHg^{-1}}`), which is all
        :meth:`action` requires.  A consumer that needs a point of
        :math:`S^2` must ask for a section, and no shipped consumer does.

        Typed honestly: the domain is this entry in its chart coordinates, the
        codomain :attr:`lift_codomain` — the ball for a sphere entry, the base
        for the identity — so a reader learns which it was handed
        (:func:`barycentre` is this arrow under its public name).
        """
        return ManifoldMap(
            domain=self, codomain=self.lift_codomain, apply=self.lift_coordinates,
        )

    def orbit_barycentres(self, points: ArrayLike) -> NDArray:
        r"""The orbit barycentre of each point, in the base's AMBIENT space —
        :math:`P_H p` for a point given in the base's coordinates (a fold's
        nodes are representatives), :math:`\lambda(c)` for one given in the
        chart's (a polar marginal's :math:`\mu`).

        The two coordinate systems an entry's points arrive in — the base's
        and the chart's — told apart by width and mapped to ONE concept: the
        point of the ambient space that determines the orbit, canonical under
        the normaliser.  ⚠ A WIDER admission than :meth:`contains` makes: an
        axial entry has no fundamental domain, so ``contains`` knows only the
        chart and REFUSES an ambient-width point, while every entry can
        project one (:math:`P_H` is defined on all of the ambient space).
        Named for what it returns on every entry — until 2026-09-03 this was
        ``ambient_representatives``, which passed ambient-width points through
        as representatives and lifted chart-width ones to barycentres: two
        things under a name that promised a representative, with a ⚠ in its
        docstring retracting the promise for the axial entry.  ⚠ The output is
        NOT claimed to lie on the base; it lies on :attr:`lift_codomain`.
        """
        arr = np.asarray(points, dtype=float)
        if arr.ndim == 1:
            arr = arr[:, None]
        width = arr.shape[1]
        if width == self.base.ambient_dim:
            chart = self.orbit_coordinates(arr)
        elif width == self.realization.ambient_dim:
            chart = arr
        else:
            raise ValueError(
                f"{self.name}: points of width {width} are in neither the base's "
                f"coordinates ({self.base.ambient_dim}) nor the chart's "
                f"({self.realization.ambient_dim})."
            )
        lifted = np.asarray(self.lift_coordinates(chart), dtype=float)
        return lifted if lifted.ndim == 2 else lifted[:, None]

    @property
    def acting_dimension(self) -> int:
        """The BASE's acting dimension: an orbit space is acted on through its
        base, so the motions are the base's, whatever the chart's width
        (:attr:`~Manifold.ambient_dim` of a quotient is the chart's).  Read through the
        base's own property, so a quotient of a quotient asks the sphere."""
        return self.base.acting_dimension

    def _refusal(self, motion: "RigidMotion") -> DoesNotMap | DoesNotNormalise | None:
        r""":math:`g` acts on :math:`M/H` iff it maps the base onto itself AND
        NORMALISES :math:`H` (:math:`gHg^{-1} = H`): then :math:`g` carries
        every :math:`H`-orbit onto an :math:`H`-orbit and the image orbit is a
        function of the orbit alone.  Otherwise :math:`[g\,p]` depends on the
        representative and the map is not a map — refused with the reason,
        rather than evaluated on whatever representative arrived (`[M]`
        :math:`C_4` about :math:`z` does not normalise :math:`\sigma_y`; it
        conjugates it to :math:`\sigma_x`).  The base's law runs first: the
        normaliser test reads the LINEAR part only, so a translated motion
        (a deck wrap) "normalises" :math:`O(2)_x`, and until S5 the induced
        action applied its translation to the barycentres, off the ball."""
        base_refusal = self.base._refusal(motion)
        if base_refusal is not None:
            return base_refusal
        if not self.by.is_normalised_by(motion):
            return DoesNotNormalise(
                f"{motion!r} does not descend to {self.name}: it does not "
                f"normalise {self.by.name} (g H g^-1 != H), so [p] -> [g p] "
                f"is not a function of the orbit. Only the normaliser of the "
                f"quotienting group acts on an orbit space."
            )
        return None

    def _image_map(self, motion: "RigidMotion") -> Callable[[NDArray], NDArray]:
        r""":math:`[p] \mapsto [g\,p]`: barycentres in (either width, through
        :meth:`orbit_barycentres`), chart coordinates out — the convention of
        :attr:`quotient_map`.

        This is the arrow the residual symmetry :math:`\Gamma` of a geometry
        acts through on a folded rule (#429 tracker 2.2b, 2026-09-02): a
        fold's nodes are representatives, so asking :math:`\sigma_y` of the
        representatives of :math:`S^2/\sigma_y` in the ambient space reads
        *"not closed"* while :math:`\sigma_y` acts TRIVIALLY on the orbit
        space — the question was being asked in the wrong place.  Named
        ``induced_action`` until S5 (2026-10-19): one name for one verb across
        the G-space family (ruling 43), the word "induced" kept here."""
        return functools.partial(_act_through, self, motion)

    def coordinates(self, points: ArrayLike) -> NDArray:
        r"""The orbit's chart coordinates, from a point in EITHER of the
        quotient's coordinate systems: two representatives of one orbit
        (a fold's node and its :math:`\sigma_y` mate) are one point."""
        return np.asarray(self.orbit_coordinates(self.orbit_barycentres(points)), dtype=float)

    def ambient_points(self, points: ArrayLike) -> NDArray:
        """The orbit barycentre in the BASE's ambient space
        (:meth:`orbit_barycentres`, ruling 46)."""
        return self.orbit_barycentres(points)

    def descending_slots(
        self, basis: Any, *, probe: ArrayLike | None = None, atol: float = 1e-12,
    ) -> NDArray:
        r"""Which slots of a basis on my BASE descend to me — the isotypic probe.

        A function on the base descends to the orbit space iff it is
        constant on the fibres of :attr:`quotient_map`, i.e. iff it takes the
        same value at a point and at every image of that point under the
        group. So the entry asks exactly that: tabulate ``basis`` at generic
        base points and at their images under the group's GENERIC elements
        (:meth:`SubgroupOfO3.generic_images
        <orpheus.geometry.symmetry.SubgroupOfO3.generic_images>` — every
        element of a finite group; for the axial groups :math:`SO(2)_a` /
        :math:`O(2)_a`, rotations by
        INCOMMENSURATE angles, because a sample of right angles generates
        :math:`C_4` and falsely admits the :math:`m = \pm 4` slots at
        :math:`L \ge 4`; ``vv-principles`` #13) and keep the slots that agree
        to ``atol`` at every image. Returns a boolean mask over the table's
        mode axes (``basis.evaluate(points).shape[1:]``).

        The probe lives HERE (user-ruled 2026-09-02) because the predicate
        is a theorem about :math:`\pi`, which this entry owns, and its two
        readers — :class:`~orpheus.numerics.basis.descent.Descent` and the
        σ-even harmonic sub-basis
        (:attr:`~orpheus.numerics.basis.spherical_harmonic_basis.MirrorEvenSphericalHarmonicBasis.even_slot_mask`)
        — must not spell it twice. `[M]` 2026-09-02, about x at
        :math:`L = 4`: exactly :math:`\{(\ell, 0)\}`, 5 real slots of 25
        (the :math:`|m| > \ell` padding is identically zero and descends
        vacuously — count real slots, not table slots); about y and z the
        invariant subspace is one-dimensional per degree but NOT slot-aligned
        from :math:`\ell \ge 2` (the harmonics' polar axis is x), which is
        the descent's business, not this probe's — it answers honestly
        about whatever slots the basis has.

        ``basis`` is any object with ``evaluate(points) -> (N, *modes)`` on
        the base (duck-typed: this module imports nothing from
        ``numerics.basis``). Only sphere entries carry a default probe.
        """
        if probe is None:
            if self.base != SPHERE:
                raise NotImplementedError(
                    f"Quotient.descending_slots: no default probe set for a "
                    f"quotient of {self.base.name}; pass generic base points "
                    f"as probe=."
                )
            probe = _GENERIC_SPHERE_PROBE
        points = np.asarray(probe, dtype=float)
        if not self.base.contains(points):
            raise ValueError(
                f"Quotient.descending_slots: the probe points are not on the "
                f"base {self.base.name}."
            )
        reference = np.asarray(basis.evaluate(points))
        mask = np.ones(reference.shape[1:], dtype=bool)
        for image in self.by.generic_images(points):
            mask &= np.all(
                np.isclose(basis.evaluate(image), reference, rtol=0.0, atol=atol),
                axis=0,
            )
        return mask

    def contains(self, points: ArrayLike) -> bool:
        r"""Membership, in EITHER of the quotient's two coordinate systems.

        A point of :math:`M/H` may be given as chart coordinates (the
        :attr:`realization`'s language) or as a representative in the base
        (the :attr:`fundamental_domain`'s).  Both are honest and the type
        knows both, so it accepts both and dispatches on the ambient
        width — the one place the distinction is a genuine local split
        rather than a repeated tag test.

        ⚠ :attr:`~Manifold.ambient_dim` still reports the REALIZATION's width: that is
        the canonical coordinate for composition (a :class:`Product`
        factor must have one width).  This method is deliberately the
        wider of the two.
        """
        arr = np.atleast_1d(np.asarray(points, dtype=float))
        chart = self.realization.ambient_dim
        if self.fundamental_domain is not None:
            section = self.fundamental_domain.ambient_dim
            width = arr.shape[1] if arr.ndim == 2 else 1
            if width == section and section != chart:
                return self.fundamental_domain.contains(arr)
        return self.realization.contains(arr)

    @property
    def is_free(self) -> bool:
        """``True`` iff the action has no fixed points, i.e. no stratum."""
        return self.singular_stratum is None

    #: The singular stratum as a LOCUS: a SymPy expression in the
    #: realization's coordinates whose vanishing set is the stratum, or
    #: ``None`` for a free action.
    #:
    #: ⛔ This was ``tuple[float, ...]`` until 2026-08-31 and could not
    #: hold the second catalogued entry: ``S^2/sigma_y``'s stratum is the
    #: disk's boundary CIRCLE, not a finite point set.  The first entry's
    #: shape had become the field's type — a stratum is a locus, and two
    #: poles are a locus that happens to be finite.
    #:
    #: It is derivation OUTPUT, not a stored copy of :attr:`det_gram`:
    #: recovering it needs the BASE's own ideal (``det P = 4 p_2`` becomes
    #: ``4(1 - mu^2)`` only after substituting ``p_1^2 + p_2 = 1``), and a
    #: :class:`Quotient` does not carry that ideal.  So the type cannot
    #: recompute it, which is exactly when storing is right.  Excluded from
    #: equality like :attr:`generators`: it is fixed by ``(base, by)``.
    singular_stratum: Any | None = field(default=None, compare=False)


# ---------------------------------------------------------------------------
# The shipped members. Until 2026-09-01 (#429 tracker 2.0c) each of these had a
# string twin in ``measure.py`` — ``SPACE_SPHERE``, ``SPACE_INTERVAL_M11``, … —
# and ``DiscreteMeasure.support`` carried the twin rather than the member. Both
# vocabularies shipped at once, which is why the comments below name what each
# member replaced: the tags are gone, and every one of these reproduces its
# tag's spelling bit-identically through :attr:`Manifold.name`.
# ---------------------------------------------------------------------------

SPHERE = Sphere()
CIRCLE = Circle()
COSINE_INTERVAL = Interval(-1.0, 1.0)  # replaced SPACE_INTERVAL_M11
UNIT_INTERVAL = Interval(0.0, 1.0)  # replaced SPACE_INTERVAL_01
HALF_LINE = Interval(0.0, float(np.inf))  # replaced SPACE_HALF_LINE
REAL_LINE = Interval(float(-np.inf), float(np.inf))  # replaced SPACE_R
POINT = Point()


# ---------------------------------------------------------------------------
# Maps between manifolds — the ARROWS of the category the members above are
# the objects of (2026-09-02, #429 tracker 2.3)
# ---------------------------------------------------------------------------

# The axis table (``AXIS_NAMES`` / ``AXIS_INDEX``) lives in
# :mod:`orpheus.geometry.faces` and is imported above. It lived HERE from
# 2026-09-02 to 2026-09-03 because this module could not import ``symmetry``
# at module scope, and moved back to ``symmetry`` (as ``AXIS_LETTER``) with R2
# of #434 when the import direction reversed (`[M]` keeping it here under the
# new direction closed a 2-cycle that killed 6 of 9 entry points); since layer
# DAG WP-B (2026-10-22) it is the face naming's table, merged with the
# point groups' transcription of it.

#: Nine GENERIC unit directions — no component zero, no two related by a
#: coordinate symmetry — the base points at which an orbit-space entry asks
#: which functions of a basis are constant on its fibres
#: (:meth:`Quotient.descending_slots`).  Deterministic (a fixed seed) so the
#: answer is reproducible, normalised so they ARE points of :math:`S^2`
#: (`[M]` 2026-09-02: the retired mirror probe's five directions had norms
#: 0.83–0.998 and were refused the day the harmonics refused off-sphere
#: points; normalising them left every mask bit-identical, 15 of 15).
_GENERIC_SPHERE_PROBE: NDArray = (
    lambda g: g / np.linalg.norm(g, axis=1)[:, None]
)(np.random.default_rng(20260902).normal(size=(9, 3)))


@dataclass(frozen=True)
class ManifoldMap:
    r"""A map :math:`\varphi : M \to N` between manifolds — an arrow, TYPED.

    The point-level analogue of a
    :class:`~orpheus.numerics.operator.LinearOperator`: where an operator
    carries its ``domain`` and ``codomain`` as function spaces, a map
    carries its two point sets, and every construction that used to be
    *"apply a callable to the nodes and name the result's manifold by
    hand"* becomes an arrow whose codomain is READ.  The two induced arrows
    on the levels above are the reason the type exists:

    * the **pushforward** of a measure,
      :meth:`DiscreteMeasure.pushforward(phi) <orpheus.numerics.measure.DiscreteMeasure.pushforward>`
      — :math:`\varphi_*\mu` lives on :math:`N` **because** :math:`\varphi`
      says so, and a measure on :math:`X \neq M` is refused;
    * the **pullback** of a function, :math:`f \mapsto f \circ \varphi`,
      which is what :func:`~orpheus.numerics.measure.DiscreteMeasure.integrate`
      evaluates on the pushed measure
      (:eq:`discrete-measure-pushforward`) and what a frame applies to its
      basis along :func:`quotient_onto` (the G0 arrow, #429 tracker 2.2) —
      :class:`~orpheus.numerics.basis.descent.Descent` reads the same map.

    ⭐ **Why this is a repair and not decoration.**  `[M]` at 2.3's opener
    (2026-09-02) the tree spelled THREE such maps around the quotient, none
    typed: the orbit retraction inside
    :meth:`~orpheus.numerics.measure.DiscreteMeasure.quotient` (a lambda
    plus a hand-named ``new_space``), the Archimedes chart inside
    ``spherical_product`` (a loop plus the literal ``support=SPHERE``), and
    the orbit-barycentre map :math:`S^2/O(2)_a \to D^3` spelled **twice**
    — honestly by ``_embedded_nodes`` (then in ``symmetry``, in ``invariance``
    from R2 until S5, which retired it for :meth:`GSpace.ambient_points`), and
    dishonestly by the
    ERR-080 forgery arm, which applies the SAME :math:`\mu \mapsto \mu\,\hat
    e_a` and declares the result on :math:`S^2`.  A codomain that is a
    field of the map cannot be forged at the call site; that is the whole
    of what this type buys.  (`[M]` 2026-09-02: :class:`Ball` had no consumer
    outside this module — ``Ball(2)`` is the :math:`\sigma_y` entry's
    :attr:`~Quotient.realization` — and ``Ball(3)`` had never been
    constructed; :func:`barycentre` is the first arrow whose CODOMAIN is a
    ball, which is what makes the forged codomain refusable in kind.)

    **Composition** is ``psi @ phi`` for :math:`\psi \circ \varphi`, refused
    unless ``phi.codomain == psi.domain`` — the same guard, at the same
    tier, as :class:`~orpheus.numerics.operator.OperatorProduct`.  The law
    that makes it a functor, and the one intrinsic property worth testing
    (there is nothing else a map of finite point sets can get wrong that
    :meth:`Manifold.contains` does not already catch):

    .. math::

       (\psi \circ \varphi)_* \mu \;=\; \psi_* (\varphi_* \mu).

    Parameters
    ----------
    domain, codomain : Manifold
        The two point sets.  Compared by VALUE (every manifold is a frozen
        value type), so a map out of ``COSINE_INTERVAL * CIRCLE`` accepts a
        measure whose support was built by the same product — and refuses
        one on ``S^2/O2_x * S^1``, which is the same numbers meaning a
        different function of direction.
    apply : callable
        The map on AMBIENT coordinates, vectorised over an
        ``(n, domain.ambient_dim)`` array (``(n,)`` is accepted for a 1-D
        domain) and returning ``(n, codomain.ambient_dim)``.  It is called
        ONCE with the whole array, so an index-based map — the orbit
        retraction ``nodes -> nodes[representative]``, a well-defined map
        on the finite node set that happens to be spelled by index — is a
        map like any other.

    Notes
    -----
    The **Jacobian is not this type's business**.  ``pushforward`` is the
    :math:`\varphi`-image (weights preserved), and a change of variables
    against a target reference measure is the caller's — the same
    asymmetry the pushforward has always documented.  The pushforward
    REFERENCE measure of a catalogued orbit space (the :math:`2\pi\,d\mu`
    of Archimedes' hat-box) is a fact about measures, not about the point
    set or the map: it is
    :func:`~orpheus.numerics.generating_measure.pushforward_reference`.

    No membership check runs on the image: :meth:`Manifold.contains` at
    measure CONSTRUCTION (tracker 0.1b) is the ruled home of that refusal,
    and a map whose ``apply`` lands outside its declared codomain is caught
    there, on the measure it produced, rather than twice.
    """

    domain: Manifold
    codomain: Manifold
    apply: Callable[[NDArray], NDArray]

    def __call__(self, points: ArrayLike) -> NDArray:
        """``apply`` on a float array — the map, as a map."""
        return np.asarray(self.apply(np.asarray(points, dtype=float)))

    def __matmul__(self, other: object) -> "ManifoldMap":
        r"""Composition :math:`\psi \circ \varphi` as ``psi @ phi``.

        Refused unless ``phi.codomain == psi.domain`` — a map out of the
        wrong point set is a type error at the level of spaces, exactly as
        an operator product over unequal spaces is.
        """
        if not isinstance(other, ManifoldMap):
            return NotImplemented
        if other.codomain != self.domain:
            raise ValueError(
                f"cannot compose: the inner map lands on "
                f"{other.codomain.name!r} but the outer map is defined on "
                f"{self.domain.name!r}. A composition psi @ phi needs "
                f"phi.codomain == psi.domain."
            )
        inner, outer = other.apply, self.apply
        return ManifoldMap(
            domain=other.domain,
            codomain=self.codomain,
            apply=lambda points: outer(inner(points)),
        )


@functools.cache
def archimedes(axis: str = "z") -> ManifoldMap:
    r"""The Archimedes chart :math:`[-1,1] \times S^1 \to S^2` about ``axis``.

    Write :math:`a` for the axis and :math:`b, c` for its two cyclic
    successors (:math:`z \to x, y`; :math:`x \to y, z`; :math:`y \to z, x`
    — a right-handed frame in every case).  Then

    .. math::

       (\mu,\ (\cos\varphi, \sin\varphi)) \;\mapsto\;
       \mu\,\hat e_a + \sqrt{1-\mu^2}\,(\cos\varphi\,\hat e_b
       + \sin\varphi\,\hat e_c),

    which for ``axis="z"`` is :eq:`product-mu-phi-cosines` verbatim:
    :math:`\mu_x = \sin\theta\cos\varphi`, :math:`\mu_y = \sin\theta
    \sin\varphi`, :math:`\mu_z = \mu`.  The domain is the PRODUCT
    manifold ``COSINE_INTERVAL * CIRCLE`` — a polar factor times a circle
    factor, in that coordinate order — so
    :func:`~orpheus.numerics.quadrature.rules_product.spherical_product`
    reads as ``(polar * azimuthal).pushforward(archimedes("z"))`` and its
    support is the chart's codomain, not a literal.

    **Named for Archimedes' hat-box theorem**, which is the statement about
    this map that transport needs: the pushforward of the uniform measure
    :math:`d\Omega` along the polar coordinate :math:`\mu = \Omega \cdot
    \hat e_a` is :math:`2\pi\,d\mu` — uniform on :math:`[-1,1]` — so a
    Gauss-Legendre rule in :math:`\mu` times any circle rule is exact on
    the sphere against Lebesgue measure.  That fact is what
    ``AngularSymmetry.reference`` answers (``LEGENDRE`` for every axis), via
    :func:`~orpheus.numerics.generating_measure.pushforward_reference`; the
    map is the geometric half of it.

    ⚠ It is a *parametrisation*, not a chart in the strict sense: the
    circle collapses at :math:`\mu = \pm 1` (both poles are the image of
    a whole fibre), which is exactly the singular stratum of
    :math:`S^2/O(2)_a` — the map is one-to-one off the stratum and the
    stratum is where it is not.  Its inverse on :math:`S^2 \setminus
    \{\pm \hat e_a\}` is the :math:`(\mu, \varphi)` chart.

    Memoised so that two requests for the same axis are ONE object, which
    is what lets a derivation agreement be stated by identity.

    Raises
    ------
    ValueError
        For an ``axis`` other than ``"x"``, ``"y"``, ``"z"``.
    """
    try:
        a = AXIS_INDEX[axis]
    except KeyError:
        raise ValueError(
            f"axis must be one of {sorted(AXIS_INDEX)}, got {axis!r}"
        ) from None
    b, c = (a + 1) % 3, (a + 2) % 3

    def apply(points: NDArray) -> NDArray:
        mu, cos_phi, sin_phi = points[:, 0], points[:, 1], points[:, 2]
        sin_theta = np.sqrt(1.0 - mu**2)
        out = np.empty((points.shape[0], 3))
        out[:, a] = mu
        out[:, b] = sin_theta * cos_phi
        out[:, c] = sin_theta * sin_phi
        return out

    return ManifoldMap(
        domain=COSINE_INTERVAL * CIRCLE, codomain=SPHERE, apply=apply
    )


@functools.cache
def barycentre(orbit_space: Quotient) -> ManifoldMap:
    r"""The orbit-barycentre map of a catalogued orbit space — its
    :attr:`~Quotient.lift`, under the name that says what it is.

    For :math:`S^2/O(2)_a` this is :math:`\mu \mapsto \mu\,\hat e_a`: an
    orbit of the axial group is the circle :math:`\{\Omega : \Omega \cdot
    \hat e_a = \mu\}` of radius :math:`\sqrt{1-\mu^2}` about the point
    :math:`\mu\,\hat e_a` on the axis — the orbit's mean under the fibre's
    uniform measure, **inside the ball, off the sphere**: :math:`1 -
    \|\mu\hat e_a\|^2 = 1 - \mu^2 = \tfrac14 \det P`, the squared orbit
    radius the catalogue records as :attr:`Quotient.det_gram`.  For
    :math:`S^2/\sigma_a` it is :math:`(x_b, x_c) \mapsto (0, x_b, x_c)`, the
    midpoint of :math:`\{p, \sigma_a p\}`, on the sphere only on the equator;
    for :math:`M/\{e\}` the identity.  Every one is the Reynolds projector
    :math:`P_H` onto :math:`H`'s fixed subspace, which is why one formula
    serves (:func:`_coordinate_chart`) and why the codomain is
    :class:`Ball` for a sphere entry and can be nothing else.  Until
    2026-09-03 (#434 R4) this arrow existed for the axial entries alone and
    the mirror entries lifted through a hemisphere section instead.

    **This is the map ERR-080 forges.**  `[M]` 2026-09-02 the tree spelled it
    twice: ``_embedded_nodes`` (then in ``symmetry``, in ``invariance`` from
    R2; :meth:`Quotient.ambient_points` since S5 — honestly: an invariance
    check wants
    the barycentre, because a rotation about :math:`a` fixes it), and
    ``Quadrature._harmonic_frame_measure``'s 1-D arm, which computed the
    same :math:`(\mu, 0, 0)` and declared it a measure on ``SPHERE``.  The
    first reads this arrow; the second was RETIRED with #429's fused commit —
    a 1-D rule's frame reads the rule's own measure and binds the Legendre
    basis on the orbit space — and it could not have been re-spelled through
    this map without lying about the codomain, which was the point.

    Not a section: a section lands ON the base by picking a representative,
    and for a positive-dimensional group none is canonical
    (:attr:`Quotient.fundamental_domain` is ``None`` for every axial entry on
    purpose).  The barycentre is canonical precisely because it is not a
    representative.  Memoised, so ``barycentre(e) is barycentre(e)``; note
    that ``e.lift`` assembles a FRESH (and ``==``-equal) arrow on every read
    from the same two fields, so an agreement between the two spellings of
    this one arrow is stated by ``==``, never by ``is``.

    Raises
    ------
    ValueError
        If ``orbit_space`` is not a :class:`Quotient` — a point of a bare
        interval is not an orbit at all.
    """
    if not isinstance(orbit_space, Quotient):
        raise ValueError(
            f"the barycentre map is defined on an orbit space M/H; got "
            f"{getattr(orbit_space, 'name', orbit_space)!r}, whose points are "
            f"not orbits."
        )
    return orbit_space.lift


# ---------------------------------------------------------------------------
# The orbit-space catalogue
# ---------------------------------------------------------------------------


def quotient_onto(source: Manifold, target: Manifold) -> ManifoldMap | None:
    r"""The quotient map ``source -> target`` when ``target`` is a quotient of ``source`` — else ``None``.

    The arrow a FRAME tabulates its basis through (#429 tracker 2.2, G0,
    user-ruled 2026-09-02: ONE predicate, the lattice one): a frame binding
    functions on ``target`` to a rule on ``source`` is admissible iff this
    arrow exists, and its table is the basis pulled back along it. Three
    honest cases:

    * ``source == target`` — the identity (equality is the special case
      :math:`K = H`; the slab rule on :math:`S^2/O(2)_x` and the Legendre
      basis on the same entry);
    * ``target`` is a quotient of ``source`` itself — the entry's own
      :attr:`Quotient.quotient_map` (a Legendre basis on a full-sphere rule:
      :math:`P_\ell(\Omega\cdot\hat e_a)`, legitimate on a Lebedev or
      level-symmetric rule);
    * both are quotients of one base and the group ``source`` SPENT is
      contained in the group ``target`` was quotiented BY (:math:`K \subseteq H`,
      :meth:`SubgroupOfO3.contains <orpheus.geometry.symmetry.SubgroupOfO3.contains>`)
      — the induced map :math:`M/K \to M/H`, applied to the source's
      representatives (section coordinates: a folded rule's nodes are points
      of the base). ⚠ A source given in its REALIZATION's coordinates against
      a different target has no chart-to-chart map here yet; no shipped rule
      needs one.

    ``None`` is the refusal — the Part I bug (full harmonics on the slab's
    :math:`S^2/O(2)_x`: ``Trivial ⊉ O2('x')``), or the full harmonics on a
    mirror fold. The Legendre basis on a :math:`\sigma_b`-fold
    (:math:`b \ne a`) is ADMITTED by the third case (#432, 2026-09-02):
    the entry is named by its stabiliser :math:`O(2)_a` and
    :math:`\sigma_b \in O(2)_a`, so the induced :math:`S^2/\sigma_b \to
    S^2/O(2)_a` exists and the fold's section coordinates are read
    straight into :math:`\mu`. Until then the basis could only declare the
    derived lower bound :math:`SO(2)_a`, and that pairing — mathematically
    admissible — was over-refused.
    """
    if source == target:
        return ManifoldMap(domain=source, codomain=target, apply=_all_coordinates)
    if (
        isinstance(source, Quotient)
        and source.base == target
        and source.by.is_trivial
    ):
        # M/{e} -> M: the trivial quotient IS the base (`_mod_trivial`'s
        # theorem, `P = I`), so the arrow is the identity on coordinates.
        # Unreachable by any shipped unfolded rule (they declare the bare
        # sphere), recorded because the lattice says it exists.
        return ManifoldMap(domain=source, codomain=target, apply=_all_coordinates)
    if isinstance(target, Quotient):
        if source == target.base:
            return target.quotient_map
        if (
            isinstance(source, Quotient)
            and source.base == target.base
            and target.by.contains(source.by)
        ):
            return ManifoldMap(
                domain=source, codomain=target, apply=target.orbit_coordinates,
            )
    return None


def _act_through(orbit_space: Quotient, motion: "RigidMotion", points: NDArray) -> NDArray:
    r""":math:`[p] \mapsto [g\,p]` — :meth:`Quotient.action`'s apply,
    a module-level function partially applied so the arrow is picklable."""
    return orbit_space.orbit_coordinates(
        motion.on_points(orbit_space.orbit_barycentres(points))
    )


def _ambient_columns(columns: int | list[int], points: NDArray) -> NDArray:
    r"""``points[:, columns]`` — every shipped entry's orbit coordinates.

    They are a SELECTION of the base's ambient coordinates because each
    surviving invariant is a coordinate function (``x_a``; ``x_b, x_c``);
    an entry whose invariants are higher-degree polynomials (``S^2/C_n``)
    would carry a genuine polynomial map here instead.  A module-level
    function partially applied, rather than a lambda, so an entry stays
    picklable and the column is recoverable from the partial's ``args``
    (the field itself is ``repr=False``).  An ``int``
    selects one column and returns ``(n,)`` — the shape a 1-D measure
    carries its nodes in; a list returns ``(n, k)``.
    """
    return np.asarray(points)[:, columns]


def _embed_columns(columns: int | list[int], ambient: int, chart: NDArray) -> NDArray:
    r"""The chart columns scattered back into a zero vector of the ambient
    width — :func:`_ambient_columns`'s right inverse, and the ``embed`` half of
    :func:`_coordinate_chart`.  An ``int`` takes ``(n,)`` chart coordinates
    (a polar marginal's :math:`\mu`), a list takes ``(n, k)``."""
    arr = np.asarray(chart, dtype=float)
    # A 1-D input means n points of a ONE-column chart, or ONE point of a
    # k-column one — the same reading `_ambient_columns` produces on the way
    # in, inherited from numpy's indexing.
    arr = arr.reshape(-1) if isinstance(columns, int) else np.atleast_2d(arr)
    out = np.zeros((arr.shape[0], ambient))
    out[:, columns] = arr
    return out


def _coordinate_chart(
    columns: int | list[int], ambient: int,
) -> "tuple[Callable[[NDArray], NDArray], Callable[[NDArray], NDArray]]":
    r"""The chart and its lift, for an entry whose surviving invariants are
    COORDINATE FUNCTIONS — ``(select, embed)``.

    ``select`` reads the invariant columns off the base's ambient coordinates
    (:func:`_ambient_columns`: the chart :math:`\pi`); ``embed`` scatters them
    back into a zero vector (:func:`_embed_columns`: the lift
    :math:`\lambda`).  Their composite ``embed ∘ select`` is the orthogonal
    projector onto the span of those coordinate axes — which for a group
    whose fixed subspace IS that span is the Reynolds projector
    :math:`P_H = \int_H \rho(g)\,dg`, the orbit barycentre.  So the pair is
    spelled ONCE, here, and the two sphere builders call it with their
    invariant columns: the axial entry with ``axis`` (fixed subspace
    :math:`\mathbb{R}\hat e_a`), the mirror entry with the two kept columns
    (fixed subspace :math:`\{x_a = 0\}`).  An entry whose invariants are not
    linear — :math:`S^2/C_n`, #439's :math:`\mathbb{RP}^2` — supplies its own
    pair; the fields are REQUIRED, so it cannot forget the lift.  Both halves
    are module-level functions partially applied, so an entry stays picklable
    and the columns are recoverable from each partial's ``args`` (the two
    fields are ``repr=False``; the ``repr`` shows neither).
    """
    return (
        functools.partial(_ambient_columns, columns),
        functools.partial(_embed_columns, columns, ambient),
    )


def _all_coordinates(points: NDArray) -> NDArray:
    """The identity — ``M/{e}``'s orbit coordinates ARE the point's, and its
    lift is the same map."""
    return np.asarray(points)


def _generic_points(base: Manifold) -> NDArray:
    r"""A SET of generic points of ``base`` — points off every stratum, over
    which the generic orbit dimension is the MAXIMUM (orbit dimension is
    upper semicontinuous: it drops only on the singular stratum, so the
    largest orbit in a probe set is the generic one).  `[M]` 2026-09-03
    (elegance review): with ONE probe row placed on the axis, a
    single-point spelling of the law refused :math:`S^2/O(2)_z` and ADMITTED
    the disk forgery — the maximum over the set is what makes "generic" a
    computation rather than a claim about a seeded draw.  The sphere's set is
    :data:`_GENERIC_SPHERE_PROBE`; a flat ambient base gets four seeded
    vectors inside the unit ball (forward machinery for #440 — `[M]` no
    shipped entry reaches it); a base no continuous group acts on has none,
    and none is asked (a finite group's orbits are points).  A ``match`` with
    a raising fallthrough: a base with a continuous action and no spelled probe
    is refused by name rather than handed a plausible set.
    """
    match base:
        case Sphere():
            return _GENERIC_SPHERE_PROBE
        case RealSpace(d=d) | Ball(d=d):
            v = np.random.default_rng(20260903).normal(size=(4, d))
            return v / (2.0 * np.linalg.norm(v, axis=1)[:, None])
    raise NotImplementedError(
        f"no generic points are spelled for {base.name}; only the sphere and "
        f"the flat ambient bases carry a continuous action today."
    )


@functools.cache
def _catalogued_quotient(base: GSpace, group: SubgroupOfO3) -> Quotient:
    """The memoised body of :meth:`GSpace.quotient` — see its Notes.

    **SCOPE-BOUNDARY[guard]** — machinery: the Gröbner-basis orbit-space
    engine, which would derive every entry. ruling: D0.1, the user,
    2026-08-31, "deferred, not refused" (:ref:`manifold-engine-seed`).
    revisit: an orbit space requested outside the catalogue, or more
    than the dozen groups transport needs. The ``NotImplementedError``
    below is the boundary's refusal; each entry is the engine's output in
    its own data model, so building the engine develops the catalogue.
    """
    # M/{e} = M is a THEOREM for every G-space, not a table row —
    # so it is derived here rather than needing one entry per member.
    # Run the same procedure on the trivial group: its invariant ring
    # is the whole polynomial ring, generated by the coordinates
    # themselves, so P_ij = <grad x_i, grad x_j> = delta_ij, and
    # det P = 1 vanishes NOWHERE.  No vanishing locus means no
    # singular stratum, which means the action is free — and the
    # trivial action is free, vacuously, since its only element is
    # the identity.  The procedure reproduces the known answer, so
    # this case is also a positive control on the machinery.
    # The trivial group is admitted on EVERY G-space base, first: M/{e} is
    # the identity ENTRY (`_mod_trivial` — a derivation output and the
    # catalogue's positive control).  And (M/H)/{e} IS M/H — the same theorem the door below
    # states for every G ⊆ H, applied to G = {e}: until 2026-09-03 (#434 R4)
    # this built a SECOND Quotient named ``S^2/sigma_y/Trivial``, unequal to
    # the entry it is, which is the two-spellings defect the door exists to
    # refuse.  Named here so that no session "tightens" the door for
    # consistency and takes the ambient container with it.
    if group.is_trivial:
        return base if isinstance(base, Quotient) else _mod_trivial(base, group)
    # A group that acts TRIVIALLY on an orbit space — one inside the group
    # already spent — has nothing to quotient: (M/H)/G for G ⊆ H is M/H
    # itself, a second spelling of one orbit space.  Refused with the theorem
    # (#429 tracker 2.2b, 2026-09-02), by the same law that names an entry
    # by its stabiliser: never "derive S^2/sigma_y/sigma_y", which would
    # send a caller off to derive an orbit space that exists.
    if isinstance(base, Quotient) and base.by.contains(group):
        raise ValueError(
            f"{base.name}/{group.name} is {base.name}: {group.name} lies in "
            f"the spent group {base.by.name}, so it acts trivially on the "
            f"orbit space and the quotient is the identity — there is "
            f"nothing to quotient. The orbit space is {base.name} already."
        )
    # The door's half of the invariant `Quotient.__post_init__` enforces:
    # refuse a non-stabiliser group BEFORE the lookup, so the answer is the
    # theorem (this orbit space exists, under its stabiliser's name) and
    # never "no catalogue entry" for an entry that exists.
    _assert_named_by_stabiliser(base, group)

    entry = _ORBIT_CATALOGUE.get((type(base), group.name))
    if entry is None:
        raise NotImplementedError(
            f"no catalogue entry for {base.name}/{group.name}: derive "
            f"it once (minimal invariants p_1..p_k of R[x]^H; syzygy "
            f"ideal I by elimination; Procesi-Schwarz P_ij = "
            f"<grad p_i, grad p_j> with P >= 0; intersect with the "
            f"ideal of {base.name}) and register it in "
            f"orpheus/geometry/manifold.py's _ORBIT_CATALOGUE, or "
            f"implement the derivation engine. "
            f"Catalogued today (manifold CLASS / group): "
            f"{sorted(f'{m.__name__}/{g}' for m, g in _ORBIT_CATALOGUE)}."
        )
    return entry(base, group)


def _sphere_mod_o2(base: GSpace, group: SubgroupOfO3) -> Quotient:
    r"""``S^2 / O(2)_a = [-1, 1]``, derived per the standard procedure.

    Write :math:`a` for the axis and :math:`b, c` for the other two.
    Invariants of :math:`O(2)_a` — the rotations about :math:`a` AND the
    mirrors through it — acting on :math:`\mathbb{R}^3` are :math:`p_1 =
    x_a` and :math:`p_2 = x_b^2 + x_c^2`; they are algebraically
    independent, so the syzygy ideal is empty.  ⭐ They are exactly the
    invariants of the rotation half :math:`SO(2)_a` alone: a vertical
    mirror fixes :math:`x_a` and :math:`x_b^2 + x_c^2`, so
    :math:`\mathbb{R}[x]^{SO(2)_a} = \mathbb{R}[x]^{O(2)_a}`, the two
    groups have the same orbits (the constant-:math:`\mu` circles), and
    this ONE derivation is the orbit space of both.  **The entry is named
    by the stabiliser** — the largest group with these orbits, which is
    :math:`O(2)_a = \{g : g\hat e_a = \hat e_a\}` — and asking for the
    quotient by :math:`SO(2)_a` is refused with that theorem, at the
    catalogue door and at :class:`Quotient`'s construction, both reading
    :attr:`~orpheus.geometry.symmetry.SubgroupOfO3.orbit_stabiliser`
    (:func:`_assert_named_by_stabiliser`; this builder is a pure
    derivation and trusts ``group``), so one
    orbit space has one spelling (#432, 2026-09-02; until then the entry
    was keyed and named by the rotation half, and the Legendre basis on it
    could declare only that lower bound).  The Procesi–Schwarz matrix is

    .. math:: P = \begin{pmatrix} 1 & 0 \\ 0 & 4 p_2 \end{pmatrix},
              \qquad \det P = 4 p_2,

    so :math:`\mathbb{R}^3/O(2)_a = \{p_2 \ge 0\}`.  Adjoining the
    sphere's own ideal :math:`p_1^2 + p_2 = 1` and writing
    :math:`\mu = p_1`:

    .. math:: S^2/O(2)_a = \{\mu \in \mathbb{R} : 1 - \mu^2 \ge 0\}
              = [-1, 1].

    ⭐ The three axes are ONE derivation reading the axis off the group —
    three catalogue keys, not three procedures — exactly as the mirrors
    are (2026-09-01, tracker 2.4).  The three orbit spaces are isometric
    and their realizations identical, and they are still three different
    quotients: :math:`\mu` is the cosine against a DIFFERENT axis in
    each, so a rule on :math:`S^2/O(2)_x` (the slab's — `[M]` the real
    spherical-harmonic pole is :math:`x`) and a rule on
    :math:`S^2/O(2)_z` (the polar factor of a product rule) carry the
    same numbers and mean different functions of direction.  That is why
    :attr:`Quotient.name` spells the axis (``S^2/O2_x``) and why the two
    do not compare equal.

    ⭐ :math:`\det P = 4(1-\mu^2)` is the squared orbit radius, and it is
    the SAME polynomial as the curvilinear angular-redistribution
    coefficient in :math:`(1/r)\,\partial_\mu[(1-\mu^2)\psi]` — the
    quotient's connection term.  Its zeros, :math:`\mu = \pm 1`, are the
    poles: the two points with full stabilizer, hence the singular
    stratum, hence why the orbit space is an orbifold.

    ⭐ The entry's remaining derivation output (#429 tracker 3.1,
    2026-09-02).  The **quotient map** is the surviving invariant read as
    a function of direction, :math:`\pi_a : \Omega \mapsto \Omega \cdot
    \hat e_a = p_1` — `[M]` :math:`\pi_a \circ \varphi_a = \mathrm{pr}_1`
    bit-exactly against the Archimedes parametrisation, on every axis.
    The measure a degree of exactness here is against (Archimedes' hat-box
    theorem, ``LEGENDRE``) is the measure layer's:
    :func:`~orpheus.numerics.generating_measure.pushforward_reference`.

    SymPy is imported lazily: it costs ~250 ms and no module of
    :mod:`orpheus.geometry` imports it at module scope (`[M]` 2026-10-22: its
    3 imports there are function-local, in this module), so a session that never asks for a
    quotient never pays for it — and :meth:`GSpace.quotient` memoises,
    so a session that asks pays the derivation once.
    """
    import sympy as sp

    axis = group.rotation_axis
    if axis is None:
        # An admission contract, so a real raise: `-O` strips `assert`.
        # Reachable only by registering this builder under a non-axial key.
        raise ValueError(
            f"the axial orbit-space derivation needs an axis, and "
            f"{group.name!r} has none. Register _sphere_mod_o2 only against "
            f"O2 entries."
        )
    coords = sp.symbols("x y z", real=True)
    x_a = coords[axis]
    x_b, x_c = (c for i, c in enumerate(coords) if i != axis)
    p1, p2 = sp.symbols("p1 p2", real=True)

    invariants = (x_a, x_b**2 + x_c**2)
    grad = [[sp.diff(p, v) for v in coords] for p in invariants]
    gram_xyz = sp.Matrix(
        [[sum(gi * gj for gi, gj in zip(a, b)) for b in grad] for a in grad]
    )
    # Re-express in the invariants: <grad p2, grad p2> = 4(x_b^2+x_c^2) = 4 p2.
    gram = gram_xyz.subs({x_b**2 + x_c**2: p2}).applyfunc(sp.simplify)
    gram = gram.subs({x_b: sp.sqrt(p2), x_c: 0})  # on the orbit, wlog
    gram = gram.applyfunc(sp.simplify)
    det_gram = sp.simplify(gram.det())

    # pi_a: the surviving invariant p_1 = x_a, read off the base's ambient
    # coordinates ((n,) out, the shape a polar marginal's nodes have); and its
    # lift mu -> mu e_a, the orbit barycentre.  One pair, spelled once.
    ambient = base.ambient_dim
    chart, lift = _coordinate_chart(axis, ambient)
    return Quotient(
        base=base,
        by=group,
        realization=COSINE_INTERVAL,
        orbit_coordinates=chart,
        lift_coordinates=lift,
        lift_codomain=Ball(ambient),
        generators=(p1 - x_a, p2 - (x_b**2 + x_c**2)),
        syzygy=(),
        # Frozen, so the Quotient — a frozen value type — is hashable, which
        # the catalogue's memo and any set/dict keyed on a support need.
        gram=sp.ImmutableMatrix(gram),
        det_gram=det_gram,
        derived_by="hand",
        # No CANONICAL section: any half-meridian sections S^2 -> S^2/O(2)_a
        # and none is distinguished, which is the normal case for a
        # positive-dimensional group.  ⛔ ERR-080 is the tree having needed
        # one anyway and fabricated it by zero-padding to (mu, 0, 0), which
        # is off S^2.  Declaring one (the phi=0 half-meridian, mu on axis a,
        # sqrt(1-mu^2) on axis b, 0 on axis c) is a CHOICE, not a derivation
        # output: tracker 2.3 declined it ([M] 0 production readers of a
        # section) and 3.1 ships the QUOTIENT MAP instead, which is one.
        fundamental_domain=None,
        # det P = 4 p_2 = 4(1 - mu^2) after the sphere ideal p_1^2 + p_2 = 1;
        # it vanishes at the poles mu = +-1.  A locus, in the realization's
        # coordinate — the two poles are a stratum that happens to be finite.
        singular_stratum=sp.simplify(1 - sp.Symbol("u0", real=True) ** 2),
    )


def _sphere_mod_mirror(base: GSpace, group: SubgroupOfO3) -> Quotient:
    r"""``S^2 / <sigma_a> = D^2``, derived per the standard procedure.

    The shipped CYLINDRICAL FOLD.  ``[M]``
    :meth:`~orpheus.numerics.quadrature.directional.Quadrature.folded_product`
    ends in ``full.quotient(SubgroupOfO3.Mirror("y"))``, so this entry is
    what makes the fold's support expressible as a typed quotient at all.

    Write :math:`a` for the mirrored axis and :math:`b, c` for the other
    two.  The invariant ring :math:`\mathbb{R}[x,y,z]^{\langle \sigma_a
    \rangle}` is minimally generated by :math:`p_1 = x_b`, :math:`p_2 =
    x_c`, :math:`p_3 = x_a^2`.  ⭐ The syzygy ideal is **empty**, and
    predictably so rather than by luck: :math:`\sigma_a` is a REFLECTION,
    so Chevalley–Shephard–Todd forces the invariant ring to be a
    polynomial ring — the three generators are algebraically independent
    by a theorem, not by inspection.  Hence

    .. math:: P = \mathrm{diag}(1,\, 1,\, 4 p_3),
              \qquad \det P = 4 p_3,

    so :math:`\mathbb{R}^3/\langle\sigma_a\rangle = \{p_3 \ge 0\}`.
    Adjoining the sphere's ideal :math:`p_1^2 + p_2^2 + p_3 = 1` and
    eliminating :math:`p_3 = 1 - p_1^2 - p_2^2` leaves the CLOSED UNIT
    DISK :math:`\{p_1^2 + p_2^2 \le 1\} = D^2`.

    ⚠ The dimension does NOT drop — :math:`\dim = 2 - 0 = 2`, because
    :math:`H` is finite.  That single fact is what makes this entry
    structurally unlike ``S^2/O(2)_a`` (where :math:`2 - 1 = 1`) and is
    why it, and not the first entry, forced the chart-vs-section ruling:
    with no reduction, the chart buys nothing and the section is
    canonical.

    :math:`\det P = 4 x_a^2` vanishes exactly on the mirror's own
    fixed-point set, the great circle :math:`x_a = 0` — which in the
    realization's coordinates is the disk's BOUNDARY, and is a circle,
    not a finite point set.

    ⭐ The entry's remaining derivation output (#429 tracker 3.1,
    2026-09-02).  The **quotient map** is the pair of surviving invariants
    read as functions of direction, :math:`\pi_a : \Omega \mapsto (x_b,
    x_c)` — the disk coordinates — and it is :math:`\sigma_a`-invariant
    bit-exactly, since the reflection touches only the column it drops.
    The pushforward of :math:`d\Omega` along it is a weighted disk measure
    that no shipped reference measure spells; that fact, and why, is
    :func:`~orpheus.numerics.generating_measure.pushforward_reference`.

    Full derivation with SymPy output, the Molien completeness check and
    the ``dim(m/m^2)`` minimality check:
    ``scratch/sigma_y_orbit_derivation.md`` (untracked; the load-bearing
    content is reproduced in ``docs/theory/foundations/manifolds.rst``).
    """
    import sympy as sp

    axis = group.mirror_axis
    if axis is None:
        # An admission contract, so a real raise: `-O` strips `assert`,
        # and this is the one precondition the derivation cannot check
        # from its own arguments.  Reachable only by registering this
        # builder under a non-Mirror key.
        raise ValueError(
            f"the reflection orbit-space derivation needs a mirror axis, "
            f"and {group.name!r} has none. Register "
            f"_sphere_mod_mirror only against Mirror entries."
        )
    coords = sp.symbols("x y z", real=True)
    x_a = coords[axis]
    kept = [c for i, c in enumerate(coords) if i != axis]
    u = sp.symbols("u0:3", real=True)

    invariants = (kept[0], kept[1], x_a**2)
    grad = [[sp.diff(q, v) for v in coords] for q in invariants]
    gram_xyz = sp.Matrix(
        [[sum(gi * gj for gi, gj in zip(a, b)) for b in grad] for a in grad]
    )
    # Re-express in the invariants: <grad p_3, grad p_3> = 4 x_a^2 = 4 p_3.
    gram = gram_xyz.subs({x_a**2: u[2]}).applyfunc(sp.simplify)
    det_gram = sp.simplify(gram.det())

    # pi_a: the two surviving invariants (x_b, x_c) — the kept columns, in
    # the order the disk's coordinates (p_1, p_2) are named; and the lift
    # (x_b, x_c) -> (0, x_b, x_c), the barycentre of {p, sigma_a p}.
    ambient = base.ambient_dim
    chart, lift = _coordinate_chart([i for i in range(ambient) if i != axis], ambient)
    return Quotient(
        base=base,
        by=group,
        # The CHART codomain, in invariant coordinates (p_1, p_2).
        realization=Ball(2),
        orbit_coordinates=chart,
        lift_coordinates=lift,
        lift_codomain=Ball(ambient),
        # The SECTION's image, in the base's ambient coordinates — what
        # DiscreteMeasure.quotient actually emits.  Closed, not strict:
        # [M] the cylindrical march seeds a level at x_a = 0 exactly.
        fundamental_domain=FundamentalDomain(
            base=base,
            normals=(tuple(1.0 if i == axis else 0.0 for i in range(3)),),
            label=f"{AXIS_NAMES[axis]}>=0",
        ),
        generators=(
            u[0] - invariants[0],
            u[1] - invariants[1],
            u[2] - invariants[2],
        ),
        syzygy=(),
        gram=sp.ImmutableMatrix(gram),
        det_gram=det_gram,
        derived_by="hand",
        # det P = 4 p_3 = 4(1 - u0^2 - u1^2) after the sphere ideal: the
        # disk's boundary circle.
        singular_stratum=sp.simplify(1 - u[0] ** 2 - u[1] ** 2),
    )


def _mod_trivial(base: GSpace, group: SubgroupOfO3) -> Quotient:
    r"""``M/{e} = M`` — the identity quotient, derived not tabulated.

    The invariant ring of the trivial group is all of
    :math:`\mathbb{R}[x_1..x_n]`, minimally generated by the coordinates,
    so :math:`P_{ij} = \langle \nabla x_i, \nabla x_j \rangle =
    \delta_{ij}` and :math:`\det P = 1`.  A determinant that vanishes
    nowhere means an empty singular stratum, i.e. a free action — which
    is right, vacuously, because the only group element is the identity.

    ⭐ Who reads it.  Not the registry: `[M]` ``AngularSymmetry.support``
    hands out the bare ``SPHERE`` for a geometry that spends nothing and
    never this entry (0 production ``SPHERE.quotient(Trivial)`` sites).  It
    is the catalogue's positive control: the procedure reproducing a known
    answer.  Until S5 (2026-10-19) it was also the orbit space a measure that
    named no group was asked ON (``invariance._ambient_orbit_space``,
    ``RealSpace(3)/{e}``, after zero-padding the nodes into
    :math:`\mathbb{R}^3`); since S5 a measure is asked on its own support,
    which must be a :class:`GSpace`.
    (Until 2026-09-03 this paragraph claimed the registry read it — the
    twin it names had been absorbed the other way, by the registry handing
    out the base.)

    The quotient map is the identity (every coordinate survives).  The
    pushforward of Lebesgue along it is Lebesgue on the BASE, which the
    measure layer does not spell generically
    (:func:`~orpheus.numerics.generating_measure.pushforward_reference`),
    and the registry's trivial arm reads ``UNIFORM_ON_SPHERE`` off the bare
    sphere it hands out as the domain — that arm never sees this entry
    (#429 tracker 3.1, user-ruled 2026-09-02).
    """
    import sympy as sp

    n = base.ambient_dim
    coords = sp.symbols(f"x0:{n}", real=True)
    return Quotient(
        base=base,
        by=group,
        realization=base,
        orbit_coordinates=_all_coordinates,
        # M/{e} = M: the lift is the identity too, and it lands on the base.
        lift_coordinates=_all_coordinates,
        lift_codomain=base,
        generators=tuple(coords),
        syzygy=(),
        gram=sp.ImmutableMatrix(sp.eye(n)),
        det_gram=sp.Integer(1),
        derived_by="hand",
        # M/{e} = M, so the identity map is a section and its image is all
        # of M — a fundamental domain cut by NO half-spaces.
        fundamental_domain=FundamentalDomain(base, (), "all"),
        # det P = 1 vanishes nowhere: no stratum, a free action (vacuously,
        # the only element being the identity).
        singular_stratum=None,
    )


def _sphere_mod_o3(base: GSpace, group: SubgroupOfO3) -> Quotient:
    r"""``S^2 / O(3) = {pt}``, derived per the standard procedure (#440).

    The invariant ring :math:`\mathbb{R}[x,y,z]^{O(3)}` is generated by
    :math:`p_1 = x^2 + y^2 + z^2` alone (the first fundamental theorem for
    the orthogonal group), so the syzygy ideal is empty and

    .. math:: P = \langle \nabla p_1, \nabla p_1 \rangle = 4 p_1,
              \qquad \det P = 4 p_1,

    so :math:`\mathbb{R}^3/O(3) = \{p_1 \ge 0\}`, the radial ray.
    Adjoining the sphere's own ideal :math:`p_1 = 1` fixes the one
    invariant: no coordinate survives, and the orbit space is a POINT. The
    dimension law agrees: :math:`O(3)` is transitive on :math:`S^2`, so its
    generic orbit is the whole sphere, and
    :math:`\dim S^2/O(3) = 2 - 2 = 0` (not :math:`2 - \dim O(3) = -1`).

    The chart selects no column (``(n, 0)`` out); the lift is the orbit
    barycentre, the Reynolds projector onto :math:`O(3)`'s fixed subspace of
    :math:`\mathbb{R}^3`, which is :math:`\{0\}`: every orbit's barycentre
    is the origin, which lies in the ball, not on the sphere. The action is
    transitive with every stabiliser conjugate to :math:`O(2)`, so there is
    one orbit type and no singular stratum; :math:`\det P = 4` on the sphere
    vanishes nowhere.

    Who reads it: collision probability's isotropic angular domain is this
    point, today an unstated convention (#440 holds the consumer wiring and
    the one-node measure of weight :math:`4\pi`, which is the measure
    layer's).
    """
    import sympy as sp

    x, y, z = sp.symbols("x y z", real=True)
    p1 = sp.Symbol("p1", real=True)
    chart, lift = _coordinate_chart([], base.ambient_dim)
    return Quotient(
        base=base,
        by=group,
        realization=POINT,
        orbit_coordinates=chart,
        lift_coordinates=lift,
        lift_codomain=Ball(base.ambient_dim),
        generators=(p1 - (x**2 + y**2 + z**2),),
        syzygy=(),
        gram=sp.ImmutableMatrix([[4 * p1]]),
        det_gram=4 * p1,
        derived_by="hand",
        # The sphere itself is ONE orbit, so any single direction is a
        # section; none is distinguished, and none is shipped.
        fundamental_domain=None,
        singular_stratum=None,
    )


#: ``(manifold type, group name) -> builder``.  Keyed on the PAIR because
#: a quotient is binary dispatch: it is a property of neither operand
#: alone.  About a dozen entries are expected; the engine that would
#: compute them instead of reading them is deferred, not refused.
_ORBIT_CATALOGUE: dict[
    tuple[type, str], Callable[[GSpace, SubgroupOfO3], Quotient]
] = {
    # The three axial orbit spaces S^2/O(2)_a share ONE derivation — it
    # reads the axis off the group — so they are three keys, not three
    # procedures.  (The key was the bare "SO2" until 2026-09-01, see
    # symmetry.SO2, and the rotation half "SO2_a" until 2026-09-02.)  The
    # entry is NAMED BY ITS STABILISER O(2)_a (#432).  Every key here is a
    # DERIVATION: a smaller orbit-equivalent group is refused at the door
    # by `orbit_stabiliser`, never routed through a key, so the "catalogued
    # today" listing in the door's message is exactly what can be built.
    (Sphere, "O2_x"): _sphere_mod_o2,
    (Sphere, "O2_y"): _sphere_mod_o2,
    (Sphere, "O2_z"): _sphere_mod_o2,
    # And the three mirrors, the same way.
    (Sphere, "sigma_x"): _sphere_mod_mirror,
    (Sphere, "sigma_y"): _sphere_mod_mirror,
    (Sphere, "sigma_z"): _sphere_mod_mirror,
    # S^2/O(3) is a point (#440): O(3) is transitive on the sphere.
    (Sphere, "O3"): _sphere_mod_o3,
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _assert_named_by_stabiliser(base: Manifold, group: SubgroupOfO3) -> None:
    r"""Refuse a group that is not its own
    :attr:`~orpheus.geometry.symmetry.SubgroupOfO3.orbit_stabiliser`.

    An orbit space is named by its stabiliser — the largest subgroup of
    :math:`O(3)` with its orbits.  A smaller orbit-equivalent group
    (:math:`SO(2)_a` for :math:`O(2)_a`; :math:`SO(3)` for :math:`O(3)`)
    has the same invariant ring, hence the same orbits, hence the same
    orbit space: a second SPELLING of one entry, which is the defect class
    ERR-080 belongs to (one orbit space, two unequal objects).  Refused
    with the theorem and the spelling to use, at BOTH doors — construction
    (:meth:`Quotient.__post_init__`) and the catalogue lookup
    (:func:`_catalogued_quotient`) — so the message is never "no catalogue
    entry" for an entry that exists (#432, 2026-09-02).  The group arrives
    as an argument, so the read is duck-typed and this module imports
    nothing for it (module docstring).
    """
    stabiliser = group.orbit_stabiliser
    if group == stabiliser:
        return
    raise ValueError(
        f"{base.name}/{group.name} is the orbit space "
        f"{base.name}/{stabiliser.name}: {group.name} and its stabiliser "
        f"{stabiliser.name} — the largest subgroup of O(3) with its orbits — "
        f"have the same invariant ring, so the same orbits and ONE orbit "
        f"space. An orbit space is named by its stabiliser; spell "
        f"{stabiliser!r} (#432)."
    )


def _fmt(v: float) -> str:
    """Format an endpoint the way the retired string tags did."""
    return str(int(v)) if float(v).is_integer() else repr(v)
