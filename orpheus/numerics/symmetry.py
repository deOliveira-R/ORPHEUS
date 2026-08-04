r"""Subgroups of :math:`O(3)` and discrete-measure invariance.

Quadrature *selection* is fundamentally a symmetry-matching problem.
A geometry has a natural symmetry group :math:`G_{\text{geom}}`
(slab :math:`\to` :math:`O(2)` about the slab normal, sphere :math:`\to`
:math:`O(3)`, hexagonal lattice :math:`\to` :math:`D_{6h}`, …) and an
angular quadrature carries an *invariance group*
:math:`G_{\text{quad}}` — the maximal subgroup of :math:`O(3)` that
permutes the quadrature nodes among themselves while preserving
weights. The selection rule is

.. math::
   :label: subgroup-of-o3-containment

   G_{\text{geom}} \subseteq G_{\text{quad}}
   \;\Longleftrightarrow\;
   \forall\, g \in G_{\text{geom}},\; g \in G_{\text{quad}}.

The "if and only if" is a tautology of subgroup containment, but its
*application* is non-trivial: given a geometry-quadrature pair, decide
whether the quadrature respects every symmetry the geometry exhibits.
A failure (e.g. integrating an :math:`O_h`-symmetric flux against an
:math:`SO(2)`-only product quadrature) introduces spurious
azimuthal-symmetry breaking that pollutes higher-order moments — the
classical observation behind Lebedev grids in computational chemistry
(Lebedev 1976, Stiefel & Fässler 1979).

ORPHEUS's relevant subgroup lattice is **finite and small**. Every
geometry-quadrature pair the project supports today (slab, sphere,
1-D / 2-D Cartesian; Gauss-Legendre, Lebedev, level-symmetric,
product) is one of eight named groups plus the parameterised families
:math:`C_n` and :math:`D_{nh}`. Generator-based machinery for arbitrary
discrete groups (`character tables, double cosets, …`) is intentionally
**not** implemented: a static ``(G_a, G_b) -> bool`` lookup table for
containment, plus per-quadrature *fingerprint checks* for invariance,
suffices for every selection decision the project will face before
hex / triangular lattices land. When :math:`C_{6v}` and :math:`D_{6h}`
become consumed, this module will be extended with the parameterised
relations they need — *not* before. The project's experience is that
"unify after two instances" beats "unify upfront with a generic
framework" (Hamermesh 1962 §2.5 covers the abstract theory if it is
ever needed).

The :meth:`SubgroupOfO3.is_invariant` check uses the **fingerprint**
strategy rather than full orbit enumeration: a measure :math:`\mu` is
:math:`G`-invariant iff the multiset of node-weight pairs is closed
under every defining permutation of :math:`G`. For the discrete groups
that arise here (:math:`O_h`, :math:`I_h`, :math:`SO(2)`), this collapses
to (i) generating a representative orbit of each node, (ii) snapping
the orbit back onto the input nodes via nearest-neighbour matching,
and (iii) confirming the matched weights agree. This is the standard
verification approach for Lebedev grids — Lebedev's construction *guarantees*
:math:`O_h` symmetry by construction, so the test confirms what is
mathematically known rather than discovering it.

References
----------

* Hamermesh, M. (1962). *Group Theory and Its Application to Physical
  Problems*. Addison-Wesley. §2.5 (finite point groups), §9.4
  (octahedral and icosahedral groups).
* Stiefel, E. and Fässler, A. (1979). *Group Theoretical Methods and
  Their Applications*. Birkhäuser. Chapters 4-5 (crystallographic
  point groups, invariant theory of finite groups).
* Lebedev, V.I. (1976). "Quadratures on a sphere." *USSR Computational
  Mathematics and Mathematical Physics* **16**(2), 10-24. The
  octahedral-invariant construction this module's :math:`O_h` check
  validates.
* Carlson, B.G. and Lathrop, K.D. (1968). "Transport theory: the
  method of discrete ordinates." In *Computing Methods in Reactor
  Physics*, Greenspan, Kelber, Okrent, eds., Gordon & Breach.
  Level-symmetric :math:`S_N` quadratures, also :math:`O_h`-invariant
  by construction.

See Also
--------

:ref:`discrete-measures` (theory page) — section "Symmetry groups for
quadrature invariance" — narrates how invariance ties into the
quadrature-selection logic that Issue 5 of the SN reshape campaign
will install.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import ClassVar, Iterable

import numpy as np

from orpheus.geometry.transformation import (
    Permutation,
    RigidMotion,
    close_group as _close_rigid_motions,
)

from .measure import DiscreteMeasure
from .roots_of_unity import roots_of_unity

#: The standard setting's coordinate frame. Every realization below is built
#: in it — principal axis along z, vertex on the x-axis — so containment is
#: literal subgroup containment *in that setting*.
_X_AXIS, _Y_AXIS, _Z_AXIS = np.eye(3)


# ---------------------------------------------------------------------------
# Group enumeration
# ---------------------------------------------------------------------------


class _NamedSubgroup(Enum):
    """Named, parameter-free subgroups of :math:`O(3)`.

    Parameterised families (:class:`Cn`, :class:`Dnh`) are separate
    types so the order parameter ``n`` is type-safe without dragging
    parameters into the enum machinery.
    """

    Trivial = "Trivial"        # {e}
    SO2 = "SO2"                # axial rotations C_∞ (continuous 1-parameter)
    Dinfh = "Dinfh"            # D_∞h — the FULL cylindrical group (see below)
    OctahedralOh = "Oh"        # full octahedral group, 48 elements
    IcosahedralIh = "Ih"       # full icosahedral group, 120 elements
    SO3 = "SO3"                # proper rotations of the sphere
    O3 = "O3"                  # orthogonal group of R^3 (top of lattice)

    # ``Dinfh`` was named ``O2`` until 2026-08-02, and the name was
    # wrong in a load-bearing way. The entry is REALIZED as axial
    # rotations + σ_h, which is C_∞h; true O(2) embedded in 3-D is
    # C_∞v (rotations + VERTICAL mirrors). Neither contains D_{nh},
    # because D_{nh} carries C₂ axes lying IN the plane — so the
    # lattice's asserted ``D_nh ⊆ O2`` was false under both readings,
    # and a committed test pinned it.
    #
    # D_∞h — C_∞ rotations, the in-plane C₂ axes, σ_h, the vertical
    # mirrors and the improper products — DOES contain every D_{nh},
    # so the relation the test asserted is correct once the name is.
    # It is also the group a cylinder actually carries, which is what
    # the geometry table needs.


@dataclass(frozen=True)
class Cn:
    """Cyclic group :math:`C_n` of :math:`n`-fold proper rotations
    about a fixed axis.

    Parameters
    ----------
    n : int
        Order of the cyclic group; must be :math:`\\ge 1`. ``Cn(1)``
        is the trivial group; ``Cn(2)`` is the rotation-only sibling
        of :class:`Mirror` — ``Cn(2)`` is specifically the
        :math:`180^{\\circ}`-ROTATION group (:math:`\\det = +1`), where a
        mirror is a reflection (:math:`\\det = -1`). Same order, different
        group, and the two spellings exist to keep them apart.
    """

    n: int

    def __post_init__(self) -> None:  # pragma: no cover - guard
        if self.n < 1:
            raise ValueError(f"Cn requires n >= 1, got n={self.n}")


@dataclass(frozen=True)
class Dnh:
    """Dihedral group :math:`D_{nh}` of order :math:`4n`.

    :math:`D_{nh}` is the full symmetry group of an :math:`n`-prism:
    :math:`n` proper rotations about the principal axis, :math:`n`
    rotations through perpendicular :math:`C_2` axes, plus a horizontal
    mirror plane and the improper rotations that combine these.

    Parameters
    ----------
    n : int
        Order of the principal axis; must be :math:`\\ge 1`. Common
        cases for ORPHEUS: ``Dnh(2)`` for rectangular Cartesian 2-D
        geometries, ``Dnh(4)`` for square lattices, ``Dnh(6)`` for
        hexagonal lattices (forward-looking — not yet consumed).
    """

    n: int

    def __post_init__(self) -> None:  # pragma: no cover - guard
        if self.n < 1:
            raise ValueError(f"Dnh requires n >= 1, got n={self.n}")


@dataclass(frozen=True)
class Mirror:
    r"""The order-2 reflection group :math:`\{e, \sigma_a\}`.

    **Parameterised by the plane, because the plane is not a
    convention.** This family replaced a parameter-free ``Z2`` entry on
    2026-08-02, and the entry was wrong in a way no test could see.

    ``Z2`` was realized as :math:`\sigma_z` while its docstring said
    *"pick* :math:`z \to -z` *as a canonical representative; any single
    reflection works, the choice is convention"*. `[M]` That is false on
    a rule this tree ships: ``product(4, 3)`` is closed under
    :math:`\sigma_z` and **not** under :math:`\sigma_x`, and ``Z2``
    answered ``True``. The claim survived because every fixture that
    reached it was symmetric in all three planes, where the answer is
    ``True`` either way — a degenerate fixture masking an overload, one
    more instance of the pattern this campaign keeps finding.

    Worse, the tag meant two DIFFERENT things by node shape: the 3-D arm
    tested :math:`\sigma_z`, the 1-D arm tested :math:`x \to -x` on the
    single coordinate (plane-free). The tree's canonical embedding of a
    polar marginal is :math:`(\mu, 0, 0)` (see
    :meth:`~orpheus.numerics.quadrature.directional.Quadrature.axis_cosines`),
    so the slab mirror :math:`\mu \to -\mu` **is** :math:`\sigma_x` — a
    different matrix from the one the 3-D arm used. `[M]` On an
    ASYMMETRIC :math:`\mu`-set the two arms disagree: the 3-D arm
    certifies ``True`` for a set that violates :math:`\mu \to -\mu`,
    which is a false certification in the dangerous direction (the
    ERR-072 family).

    Naming the plane is what makes the two arms say the same thing. On a
    1-D measure embedded as :math:`(\mu, 0, 0)`, :math:`\sigma_y` and
    :math:`\sigma_z` fix every node **pointwise** — so they hold
    trivially, and :math:`\sigma_x` is the one real test. That is now
    derived from the embedding rather than asserted by a branch.

    Parameters
    ----------
    axis : str
        The mirror's NORMAL, one of ``"x"`` / ``"y"`` / ``"z"``; the
        plane is the one perpendicular to it through the origin. So
        ``Mirror("z")`` is :math:`\sigma_h` for an axial geometry and
        ``Mirror("x")`` is the slab's polar mirror.

    See Also
    --------
    Cn
        The PROPER order-2 sibling is ``Cn(2)`` — a rotation, not a
        reflection. Keeping the two spellings distinct is the point:
        ``Mirror`` has :math:`\det = -1` and therefore lies in no
        proper-rotation group.
    """

    axis: str

    def __post_init__(self) -> None:  # pragma: no cover - guard
        if self.axis not in ("x", "y", "z"):
            raise ValueError(
                f"Mirror requires axis in x/y/z, got {self.axis!r}. The "
                f"axis names the mirror's NORMAL; there is no unnamed "
                f"reflection, which is exactly what the retired Z2 tag "
                f"pretended there was."
            )


# Public type for "any subgroup tag". The selection logic below
# branches on type — Cn, Dnh and Mirror are kept as separate dataclasses
# (not enum entries) so their parameter is explicit.
SubgroupTag = "_NamedSubgroup | Cn | Dnh | Mirror"


# ---------------------------------------------------------------------------
# Static containment lattice
# ---------------------------------------------------------------------------

# Edges of the partial order. Reading direction: ``(A, B)`` means
# :math:`A \subseteq B`. Reflexivity (G ⊂ G) is handled by the lookup
# function, not stored here. Transitivity is *NOT* expanded
# automatically — every needed relation is listed explicitly so the
# table reads as a one-shot lookup.
#
# References for the relations:
#   - Hamermesh 1962 §9.4 (Oh, Ih)
#   - Stiefel & Fässler 1979 Ch. 4 (crystallographic subgroups)
#
# The improper-rotation status of Oh and Ih (both contain inversion
# and reflections) means they are NOT subgroups of SO(3); they sit
# directly under O(3). This matches the chemistry/Lebedev convention
# used for the project's quadratures.
_NAMED_LATTICE: set[tuple[_NamedSubgroup, _NamedSubgroup]] = {
    # Trivial ⊂ everything (handled implicitly + explicit edges below).
    (_NamedSubgroup.Trivial, _NamedSubgroup.SO2),
    (_NamedSubgroup.Trivial, _NamedSubgroup.Dinfh),
    (_NamedSubgroup.Trivial, _NamedSubgroup.OctahedralOh),
    (_NamedSubgroup.Trivial, _NamedSubgroup.IcosahedralIh),
    (_NamedSubgroup.Trivial, _NamedSubgroup.SO3),
    (_NamedSubgroup.Trivial, _NamedSubgroup.O3),

    # The reflection group's edges are NOT here, and the reason is worth
    # stating: it is `Mirror(axis)`, a parameterised dataclass, and this
    # table is typed enum-to-enum. Its relations live in `_contains`'s
    # isinstance arms, exactly as `Cn`/`Dnh`'s do.
    #
    # `[M]` That costs less than it sounds: of the five edges the old
    # parameter-free ``Z2`` entry carried here, THREE WERE ALREADY DEAD
    # CODE — `_contains` decides finite × finite by computed matrix
    # containment and only consults this table when one side is
    # continuous, so `Trivial ⊆ Z2`, `Z2 ⊆ Oh` and `Z2 ⊆ Ih` were never
    # read. A reflection family owes exactly TWO facts, `Mirror ⊆ D_∞h`
    # and `Mirror ⊆ O(3)`, and gets every finite relation — including
    # the n-dependent `D_nh ⊇ Mirror(x)` answers — for free.
    #
    # The improper-element argument that shaped those edges still holds
    # and now lives on the arms: det(σ) = −1, so a reflection lies in no
    # PROPER-rotation group. (A ``(Z2, SO3)`` edge asserted otherwise
    # until 2026-08-02 and broke monotonicity on any measure that is
    # SO3-invariant but not reflection-symmetric.) The proper order-2
    # sibling is ``Cn(2)``, which IS inside SO2 and SO3 — the distinction
    # the two spellings exist to carry.

    # SO2 ⊂ Dinfh, SO3, O3.
    (_NamedSubgroup.SO2, _NamedSubgroup.Dinfh),
    (_NamedSubgroup.SO2, _NamedSubgroup.SO3),
    (_NamedSubgroup.SO2, _NamedSubgroup.O3),

    # Dinfh ⊂ O3 (the projection R^2 → R^3 along a fixed axis).
    (_NamedSubgroup.Dinfh, _NamedSubgroup.O3),

    # Oh, Ih ⊂ O3 (NOT ⊂ SO3 — both contain improper rotations).
    (_NamedSubgroup.OctahedralOh, _NamedSubgroup.O3),
    (_NamedSubgroup.IcosahedralIh, _NamedSubgroup.O3),

    # SO3 ⊂ O3.
    (_NamedSubgroup.SO3, _NamedSubgroup.O3),
}


def _named_contains(outer: _NamedSubgroup, inner: _NamedSubgroup) -> bool:
    """Static lookup for ``outer ⊃ inner`` (equivalently ``inner ⊂ outer``).

    Reflexive (G ⊃ G is always true) and uses the explicit edges above
    for everything else. The lattice is small enough that exhaustive
    listing of strict edges is cleaner than a transitive-closure
    computation.
    """
    if outer is inner:
        return True
    return (inner, outer) in _NAMED_LATTICE


# ---------------------------------------------------------------------------
# Public façade
# ---------------------------------------------------------------------------


class SubgroupOfO3:
    """Subgroup-of-:math:`O(3)` lattice with containment + invariance checks.

    Class attributes (named entries) are pre-instantiated singletons:

    .. code-block:: python

        SubgroupOfO3.Trivial
        SubgroupOfO3.Mirror('z')
        SubgroupOfO3.SO2
        SubgroupOfO3.Dinfh
        SubgroupOfO3.OctahedralOh
        SubgroupOfO3.IcosahedralIh
        SubgroupOfO3.SO3
        SubgroupOfO3.O3

    Parameterised families are constructors:

    .. code-block:: python

        SubgroupOfO3.Cn(6)        # cyclic order 6
        SubgroupOfO3.Dnh(6)       # dihedral D_6h (hex lattice)

    Containment via :meth:`contains` (or its synonym
    :meth:`is_subgroup_of` for the reverse-direction reading)
    implements :eq:`subgroup-of-o3-containment` over the static
    lattice; invariance via :meth:`is_invariant` checks a
    :class:`~orpheus.numerics.measure.DiscreteMeasure` against the
    group's defining permutations.
    """

    # The actual storage: a tag (named enum or parameterised dataclass).
    __slots__ = ("_tag",)

    # Pre-instantiated named singletons — assigned at module scope below,
    # after the class and ``_NamedSubgroup`` are fully defined.  Declared
    # here as ClassVars so the public ``SubgroupOfO3.SO2`` / ``.OctahedralOh``
    # / ... access surface is statically known (these are class attributes,
    # not instance slots, so they coexist with ``__slots__``).
    Trivial: ClassVar[SubgroupOfO3]
    SO2: ClassVar[SubgroupOfO3]
    Dinfh: ClassVar[SubgroupOfO3]
    OctahedralOh: ClassVar[SubgroupOfO3]
    IcosahedralIh: ClassVar[SubgroupOfO3]
    SO3: ClassVar[SubgroupOfO3]
    O3: ClassVar[SubgroupOfO3]

    def __init__(self, tag: "_NamedSubgroup | Cn | Dnh | Mirror") -> None:
        self._tag = tag

    # --- Constructors ----------------------------------------------------

    @classmethod
    def _from_named(cls, name: _NamedSubgroup) -> "SubgroupOfO3":
        return cls(name)

    @classmethod
    def Cn(cls, n: int) -> "SubgroupOfO3":
        """Cyclic group :math:`C_n` of :math:`n`-fold proper rotations."""
        return cls(Cn(n))

    @classmethod
    def Dnh(cls, n: int) -> "SubgroupOfO3":
        """Dihedral group :math:`D_{nh}` of order :math:`4n`."""
        return cls(Dnh(n))

    @classmethod
    def Mirror(cls, axis: str) -> "SubgroupOfO3":
        r"""Reflection group :math:`\{e, \sigma_a\}` about the plane
        normal to ``axis``. Replaced the parameter-free ``Z2`` on
        2026-08-02 — see :class:`Mirror` for why the plane is not a
        convention."""
        return cls(Mirror(axis))

    # --- Equality / hashing / repr --------------------------------------

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, SubgroupOfO3):
            return NotImplemented
        return self._tag == other._tag

    def __hash__(self) -> int:
        return hash(self._tag)

    def __repr__(self) -> str:
        tag = self._tag
        if isinstance(tag, _NamedSubgroup):
            return f"SubgroupOfO3.{tag.name}"
        if isinstance(tag, Cn):
            return f"SubgroupOfO3.Cn({tag.n})"
        if isinstance(tag, Dnh):
            return f"SubgroupOfO3.Dnh({tag.n})"
        if isinstance(tag, Mirror):
            # MUST carry the axis: `maximal_invariance_groups` keys its
            # `visited` set and `_GROUP_CACHE` on repr(), so a repr that
            # dropped the plane would silently collapse the three mirrors
            # into one cache entry.
            return f"SubgroupOfO3.Mirror({tag.axis!r})"
        return f"SubgroupOfO3({tag!r})"

    @property
    def name(self) -> str:
        """Human-readable name for the group (used in error messages)."""
        tag = self._tag
        if isinstance(tag, _NamedSubgroup):
            return tag.value
        if isinstance(tag, Cn):
            return f"C_{tag.n}"
        if isinstance(tag, Dnh):
            return f"D_{tag.n}h"
        if isinstance(tag, Mirror):
            return f"sigma_{tag.axis}"
        return repr(tag)

    # --- Containment lattice --------------------------------------------

    def contains(self, other: "SubgroupOfO3") -> bool:
        r"""``True`` iff ``other ⊆ self`` in the :math:`O(3)` lattice.

        Implements :eq:`subgroup-of-o3-containment` via a static lookup
        over named entries plus structural rules for the parameterised
        families:

        - :math:`C_n \subseteq C_m` iff :math:`n \mid m`. (Cyclic
          subgroups by divisibility.)
        - :math:`C_n \subseteq D_{mh}` iff :math:`n \mid m`.
          (:math:`D_{mh}` contains :math:`C_m` as the principal-axis
          subgroup.)
        - :math:`D_{nh} \subseteq D_{mh}` iff :math:`n \mid m` *and*
          :math:`n` and :math:`m` share the same principal axis (which
          we assume implicitly — multi-axis containment is out of scope
          for the static-lattice approach).
        - :math:`C_n` and :math:`D_{nh}` both sit below :math:`O(3)`;
          neither sits inside the named compact subgroups :math:`O_h`,
          :math:`I_h` (those have specific finite cyclic subgroups —
          :math:`C_4` ⊂ :math:`O_h`, :math:`C_5` ⊂ :math:`I_h`, etc. —
          but exhaustive enumeration of which :math:`n` lie inside is
          out of the static-lattice scope and not consumed today).
        - :math:`SO(2)` is the union :math:`\bigcup_n C_n`, so
          :math:`C_n \subseteq SO(2)` for every :math:`n`. Likewise
          :math:`D_{nh} \subseteq D_{\infty h}` for every :math:`n` —
          note :math:`D_{nh} \not\subseteq O(2)` under either
          embedding, since :math:`D_{nh}` carries :math:`C_2` axes
          lying IN the plane.

        Reflexivity (``self.contains(self) == True``) is honoured for
        every group.

        Parameters
        ----------
        other : SubgroupOfO3
            The candidate subgroup.

        Returns
        -------
        bool
            ``True`` iff every element of ``other`` is also an element
            of ``self``.
        """
        return _contains(self._tag, other._tag)

    def is_subgroup_of(self, other: "SubgroupOfO3") -> bool:
        """Reverse-direction synonym: ``self ⊆ other``.

        Equivalent to ``other.contains(self)`` and provided for
        readability at consumer sites where the natural phrasing is
        "the geometry's group is a subgroup of the quadrature's group."
        """
        return other.contains(self)

    # --- Invariance check -----------------------------------------------

    def is_invariant(
        self,
        measure: DiscreteMeasure,
        *,
        atol: float = 1e-13,
    ) -> bool:
        r"""``True`` iff ``measure`` is invariant under every element of
        ``self``.

        For a finite group :math:`G` acting on the ambient space of the
        measure's nodes, invariance means: for every :math:`g \in G`,
        the action :math:`g` permutes the nodes among themselves and
        the matched node carries the same weight (within ``atol``).

        For continuous groups (:math:`SO(2)`, :math:`SO(3)`, :math:`O(3)`)
        the check uses a representative orbit (a non-trivial finite
        rotation) and asserts the measure is closed; this is *necessary
        but not sufficient* for general measures, but **sufficient by
        construction** for the rules ORPHEUS ships, all of which are
        either built to be group-invariant (Lebedev, level-symmetric)
        or live in a space where the group action is trivial
        (1-D measure on :math:`[-1, 1]` is trivially :math:`SO(2)`-invariant
        because there is no azimuthal coordinate to rotate).

        Strategy by group:

        - **Trivial**: always ``True`` (every measure is invariant
          under the identity).
        - **Mirror / Cn / SO2 / Dinfh / SO3 / O3 / Dnh**: handled via the
          1-D vs sphere dispatch — see :func:`_check_invariance` below.
        - **OctahedralOh**: requires (i) closure under the 6
          coordinate-axis sign flips (:math:`x \to -x`, etc.) and
          (ii) closure under coordinate permutations (24 even
          permutations + 24 odd = the symmetric group :math:`S_3`
          action on the axes). The fingerprint is the multiset of
          :math:`(|x|, |y|, |z|)` triples sorted, weighted.
        - **IcosahedralIh**: 120 elements; the multiset of node
          radii (all equal to 1 for an :math:`S^2` rule) plus the
          icosahedral symmetry of the polyhedron is checked via a
          12-element representative orbit on a regular icosahedron
          vertex set. Lebedev grids do **not** carry :math:`I_h`
          symmetry (they are :math:`O_h`-invariant by construction),
          so this check correctly rejects them.

        Parameters
        ----------
        measure : DiscreteMeasure
            The measure to test. Must live on either a 1-D space
            (``"[-1,1]"``, ``"[0,1]"``, …) or :math:`S^2`
            (``"S^2"``); other spaces fall back to a conservative
            "trivial group invariance only" semantics.
        atol : float, optional
            Tolerance for weight equality between matched orbit
            partners. Default ``1e-13`` matches the floating-point
            noise on a 1-D Gauss-Legendre weight computation.

        Returns
        -------
        bool
            ``True`` if ``measure`` is :math:`G`-invariant within
            ``atol``, ``False`` otherwise.
        """
        return _check_invariance(self._tag, measure, atol)


# ---------------------------------------------------------------------------
# Containment lookup
# ---------------------------------------------------------------------------


def _realized_ops(tag) -> "list[RigidMotion] | None":
    r"""The FINITE matrix realization of ``tag``, or ``None`` if continuous.

    Single-sources the operator sets from the very functions
    :func:`_invariance_on_points` applies, so ``contains`` and
    ``is_invariant`` answer questions about the *same* group. A
    hand-maintained containment table is a claim with no construction
    behind it, and this module already shipped two such claims that
    were false.

    Every realization is in the **standard setting**: principal axis
    along z, vertex on the x-axis. Containment is therefore literal
    subgroup containment *in that setting* — which is the question the
    quadrature-selection gate asks, because a rule and its geometry are
    used in one frame. The setting-independent relation (subconjugacy,
    :math:`\exists g: gHg^{-1} \le K`) is a different question and gets
    its own verb; conflating them is what let the two senses disagree
    without anything noticing.
    """
    if isinstance(tag, _NamedSubgroup):
        if tag is _NamedSubgroup.Trivial:
            return [RigidMotion.identity(3)]
        if tag is _NamedSubgroup.OctahedralOh:
            return _octahedral_ops()
        if tag is _NamedSubgroup.IcosahedralIh:
            return _icosahedral_ops()
        return None  # SO2 / Dinfh / SO3 / O3 are continuous — no finite set
    if isinstance(tag, Mirror):
        # One generator, so the closed group is {e, sigma_a}: order 2.
        return _reflections(tag.axis)

    if isinstance(tag, Cn):
        return _cyclic_ops(tag.n)
    if isinstance(tag, Dnh):
        return _cyclic_ops(tag.n) + _reflections("z") + _vertical_mirrors(tag.n)
    return None


def _close_group(ops: "list[RigidMotion]") -> "list[RigidMotion]":
    """The group generated by ``ops`` — closure under composition.

    The realizations above are GENERATING sets, not full groups
    (:math:`D_{nh}` ships :math:`2n+1` matrices for a group of order
    :math:`4n`). Orbit closure only needs generators, but containment
    needs the elements, so the closure is taken here.

    A thin wrapper on :func:`~orpheus.geometry.transformation.close_group`,
    kept because callers here always mean *the point group in the standard
    3-D setting*. The rounding-key trick, the quadratic-not-cubic membership
    test and the finite-group cap all live in the core now, where the cap
    raises a NAMED refusal rather than a bare ``ValueError``.
    """
    return list(_close_rigid_motions(ops, dimension=3))


_GROUP_CACHE: "dict[str, list[RigidMotion]]" = {}


def _group_elements(tag) -> "list[RigidMotion] | None":
    """Memoised :func:`_close_group` of ``tag``'s realization.

    The lattice is queried O(n^3) times by the transitivity law, and every
    query would otherwise re-derive the same closures. Tags are immutable
    value objects, so their ``repr`` is a sound cache key.
    """
    key = repr(tag)
    cached = _GROUP_CACHE.get(key)
    if cached is not None:
        return cached
    ops = _realized_ops(tag)
    if ops is None:
        return None
    elems = _close_group(ops)
    _GROUP_CACHE[key] = elems
    return elems


def _finite_contains(outer, inner) -> bool:
    """``True`` iff every element of ``<inner>`` lies in ``<outer>``."""
    outer_elems, inner_elems = _group_elements(outer), _group_elements(inner)
    assert outer_elems is not None and inner_elems is not None
    return all(
        any(M.approximately_equals(E, atol=1e-9) for E in outer_elems)
        for M in inner_elems
    )


def _contains(outer, inner) -> bool:
    """Dispatch table for :meth:`SubgroupOfO3.contains`.

    Finite-vs-finite is COMPUTED from the realized operator sets;
    anything involving a continuous group falls back to the static
    named-lattice, since a continuous group has no finite realization
    to compare against.
    """
    # Reflexivity for parameterised families.
    if outer == inner:
        return True

    # Both finite — decide by literal matrix-set containment.
    if _realized_ops(outer) is not None and _realized_ops(inner) is not None:
        return _finite_contains(outer, inner)

    # Both named — use the static lattice.
    if isinstance(outer, _NamedSubgroup) and isinstance(inner, _NamedSubgroup):
        return _named_contains(outer, inner)

    # ----- Inner is named, outer is parameterised --------------------
    if isinstance(outer, Cn) and isinstance(inner, _NamedSubgroup):
        # The only named subgroup of an arbitrary Cn is the trivial group.
        # (Cn(2) and Mirror share an order but not a determinant: Cn(2) is
        # a proper rotation, a reflection is not. They are different
        # groups and the two spellings exist to keep them apart.)
        return inner is _NamedSubgroup.Trivial

    if isinstance(outer, Dnh) and isinstance(inner, _NamedSubgroup):
        # D_nh contains the trivial group and (for n >= 1) the
        # principal-axis cyclic subgroup C_n. Higher named subgroups
        # (SO2, Dinfh, ...) are not contained in a finite dihedral group.
        #
        # `Mirror` is deliberately absent from this arm: both sides are
        # finite, so the branch above already decided it by computed
        # matrix containment — and the computed answer is n-DEPENDENT,
        # which is the whole reason the plane had to be named. `[M]`
        # Dnh(1) ⊇ Mirror(z) but NOT Mirror(x); Dnh(2) ⊇ both; Dnh(3)
        # again only z. The retired `Z2` arm returned an unconditional
        # True here, and a committed test asserted it as "a single
        # reflection sits inside every D_nh" — false for two of the four
        # orders it was parametrized over.
        return inner is _NamedSubgroup.Trivial

    # ----- Outer is named, inner is parameterised --------------------
    if isinstance(outer, _NamedSubgroup) and isinstance(inner, Cn):
        # Cn ⊂ SO2 for every n (SO(2) = union over n of C_n).
        # Cn ⊂ Dinfh, SO3, O3 transitively. Cn ⊂ Oh / Ih only for specific
        # n that match the polyhedral rotation axes — out of scope for
        # the static lattice (Cn ⊂ Oh would need n ∈ {1,2,3,4,6}; we do
        # not encode this until a consumer needs it).
        if outer in (
            _NamedSubgroup.SO2,
            _NamedSubgroup.Dinfh,
            _NamedSubgroup.SO3,
            _NamedSubgroup.O3,
        ):
            return True
        if outer is _NamedSubgroup.Trivial:
            return inner.n == 1
        return False

    if isinstance(outer, _NamedSubgroup) and isinstance(inner, Dnh):
        # Dnh ⊂ Dinfh, O3 always. Not ⊂ SO2 / SO3 (Dnh contains improper
        # rotations).
        return outer in (_NamedSubgroup.Dinfh, _NamedSubgroup.O3)

    if isinstance(outer, _NamedSubgroup) and isinstance(inner, Mirror):
        # The two facts a reflection family owes the CONTINUOUS groups —
        # every finite relation is computed above. det(sigma) = -1, so a
        # reflection lies in no proper-rotation group; D_inf_h carries
        # sigma_h AND every sigma_v, so it contains all three.
        #
        # Reaching the fallthrough here instead would be silent and
        # wrong, not an error: `_contains` ends in a bare `return False`,
        # and `O3 not-contains Mirror` would break the soundness
        # precondition `maximal_invariance_groups` states for the walk.
        return outer in (_NamedSubgroup.Dinfh, _NamedSubgroup.O3)

    if isinstance(outer, Mirror) and isinstance(inner, _NamedSubgroup):
        # {e, sigma} has exactly two subgroups, itself and {e}.
        return inner is _NamedSubgroup.Trivial

    # ----- Both parameterised ---------------------------------------
    if isinstance(outer, Cn) and isinstance(inner, Cn):
        return inner.n != 0 and (outer.n % inner.n == 0)

    if isinstance(outer, Dnh) and isinstance(inner, Cn):
        # C_n ⊂ D_mh iff n divides m (C_n is the principal-axis
        # rotation subgroup of D_mh).
        return inner.n != 0 and (outer.n % inner.n == 0)

    if isinstance(outer, Cn) and isinstance(inner, Dnh):
        # D_mh always contains a reflection — never inside a pure cyclic
        # rotation group, except for the degenerate case D_1h = {e, σ_h}
        # which has order 2 and is not C_n for any n.
        return False

    if isinstance(outer, Dnh) and isinstance(inner, Dnh):
        return inner.n != 0 and (outer.n % inner.n == 0)

    return False


# ---------------------------------------------------------------------------
# Invariance check
# ---------------------------------------------------------------------------


def _embedded_nodes(measure: DiscreteMeasure) -> np.ndarray:
    r"""The measure's nodes as points of :math:`\mathbb{R}^3`.

    The tree's canonical embedding, written down in
    :meth:`Quadrature.axis_cosines` and used by ``spherical_harmonics``
    internally: a polar marginal :math:`\mu` becomes :math:`(\mu, 0, 0)`, a
    planar rule :math:`(x, y)` becomes :math:`(x, y, 0)`. It is the *data*
    that is lifted, not the group — :class:`SubgroupOfO3`'s named entries
    (:math:`O_h`, :math:`I_h`) genuinely are three-dimensional, and there is
    nothing to restrict them to.
    """
    nodes = measure.nodes
    if nodes.ndim == 1:
        nodes = nodes[:, None]
    n, d = nodes.shape
    if d == 3:
        return nodes
    embedded = np.zeros((n, 3))
    embedded[:, :d] = nodes
    return embedded


def _check_invariance(tag, measure: DiscreteMeasure, atol: float) -> bool:
    r"""Dispatch for :meth:`SubgroupOfO3.is_invariant` — **one arm**.

    Every measure is embedded in :math:`\mathbb{R}^3` and asked the same
    question: does every element of the group carry this point set onto
    itself, weights and all?

    **This replaced a separate 1-D arm, and that arm was over-promising.**
    `[M]` It had exactly ONE discriminating branch — the
    :math:`\mu \to -\mu` reflection test — and waved every other group
    through: :math:`SO(2)` and :math:`C_n` returned ``True``
    unconditionally, :math:`\sigma_y` and :math:`\sigma_z` returned ``True``
    unconditionally, and :math:`D_{\infty h} / SO(3) / O(3) / O_h / I_h /
    D_{nh}` all ran *the same* reflection test. So its output was a **one-bit
    function of the input** wearing a nineteen-group lattice walk: a symmetric
    rule passed the one test, every group read ``True``, and the maximal
    element was necessarily the TOP of the lattice. `[M]`
    ``Sym(gauss_legendre_on_mu(8))`` was reported as :math:`O(3)`; it is
    :math:`\{e, \sigma_x\}` together with the mirrors that fix the embedded
    axis pointwise.

    Worse in the dangerous direction: `[M]` an ASYMMETRIC 1-D measure read
    :math:`\sigma_x`-non-invariant (right) but :math:`SO(2)`- and
    :math:`C_4`-invariant (wrong), certifying a CONTINUOUS group that was
    never tested — ERR-072's shape, and in direct contradiction with the
    exact criterion :func:`_is_axis_supported` that the same module applies
    to three-dimensional nodes.

    The cause was a **convention conflict inside one function**: the mirror
    branch was derived from the :math:`(\mu, 0, 0)` embedding (so
    :math:`\sigma_x` is the real test), while the rotational branch asserted
    that ":math:`C_n` rotates about z, which does not move the polar cosine"
    — true only if :math:`\mu` sits on **z**. Under :math:`(\mu,0,0)` a
    rotation about z moves the node, and the two branches were answering
    questions about different embeddings.

    Note what monotonicity could NOT see. The compatibility law
    :math:`A \subseteq B \wedge P(B) \Rightarrow P(A)` measured **zero**
    violations here — because when *everything* reads ``True`` the
    implication is vacuously satisfied. **A consistency law is blind to
    uniform over-certification**; only comparing against a computed answer
    catches it.
    """
    # Trivial is always invariant.
    if isinstance(tag, _NamedSubgroup) and tag is _NamedSubgroup.Trivial:
        return True
    if isinstance(tag, Cn) and tag.n == 1:
        return True
    return _invariance_on_points(
        tag, _embedded_nodes(measure), measure.weights, atol,
    )


def _invariance_on_points(
    tag,
    nodes: np.ndarray,
    weights: np.ndarray,
    atol: float,
) -> bool:
    """Invariance check for 3-D-node measures (typically on :math:`S^2`)."""
    if isinstance(tag, _NamedSubgroup):
        if tag is _NamedSubgroup.SO2:
            # DECIDED EXACTLY, never sampled (ERR-072).  A continuous
            # group cannot be tested by a finite sample: the sample
            # generates a finite SUBgroup, and closure under that
            # subgroup is strictly weaker than closure under G.
            return _is_axis_supported(nodes, atol)

        if tag is _NamedSubgroup.Dinfh:
            # D_∞h = C_∞ + in-plane C₂ axes + σ_h + σ_v + improper
            # products.  Its continuous factor C_∞ forces axis support
            # exactly as SO(2) does; on an axis-supported set the
            # in-plane C₂ and σ_v act identically to σ_h (both send
            # (0,0,z) -> (0,0,-z)), so σ_h closure is the whole of the
            # remaining condition.  Both conjuncts exact.
            return _is_axis_supported(nodes, atol) and (
                _orbit_closure(nodes, weights, _reflections("z"), atol) is not None
            )

        if tag is _NamedSubgroup.OctahedralOh:
            return _orbit_closure(nodes, weights, _octahedral_ops(), atol) is not None

        if tag is _NamedSubgroup.IcosahedralIh:
            return _orbit_closure(nodes, weights, _icosahedral_ops(), atol) is not None

        if tag is _NamedSubgroup.SO3:
            # DECIDED EXACTLY (ERR-072).  Every SO(3) orbit of a
            # non-origin point is a whole 2-sphere, so a FINITE set is
            # SO(3)-closed iff it is supported at the origin.
            return _is_origin_supported(nodes, atol)

        if tag is _NamedSubgroup.O3:
            # O(3) contains SO(3), so the same exact criterion binds.
            return _is_origin_supported(nodes, atol)

    if isinstance(tag, Mirror):
        # ONE named plane, not "a reflection". The retired Z2 arm hard-
        # coded _reflections("z") under a docstring calling the choice a
        # convention; `[M]` product(4, 3) is sigma_z-closed and NOT
        # sigma_x-closed, so it was not.
        return _orbit_closure(
            nodes, weights, _reflections(tag.axis), atol,
        ) is not None

    if isinstance(tag, Cn):
        # Cyclic group about z, n proper rotations.
        return _orbit_closure(
            nodes, weights, _cyclic_ops(tag.n), atol,
        ) is not None

    if isinstance(tag, Dnh):
        # Dihedral D_nh: C_n + horizontal mirror + n vertical mirrors.
        ops = _cyclic_ops(tag.n) + _reflections("z") + _vertical_mirrors(tag.n)
        return _orbit_closure(nodes, weights, ops, atol) is not None

    return False


# ---------------------------------------------------------------------------
# Group operation generators
# ---------------------------------------------------------------------------
#
# Each operator is a `RigidMotion` — an element of E(3), not a bare matrix —
# so it carries its own composition, inverse, determinant and action, and the
# orbit-closure check below asks IT whether it permutes a point set. The
# constructions all live in `orpheus.geometry.transformation`; what remains
# here is the translation from a TAG (`Mirror("x")`, `Dnh(6)`) into the
# generating set that realizes it in the standard 3-D setting.


#: The coordinate axes, by name. The mirror family names its plane by the
#: NORMAL, so ``_AXIS_INDEX["x"]`` selects the normal of :math:`\sigma_x`.
_AXIS_INDEX = {"x": 0, "y": 1, "z": 2}


def _axis_vector(axis: str) -> np.ndarray:
    """The unit vector along a named coordinate axis."""
    try:
        return np.eye(3)[_AXIS_INDEX[axis]]
    except KeyError:
        raise ValueError(f"axis must be x/y/z, got {axis!r}") from None


def _reflections(axis: str) -> list[RigidMotion]:
    r"""The single coordinate mirror whose NORMAL is ``axis``.

    :math:`\sigma_x` is the reflection in the plane :math:`x = 0`, i.e. the
    one with normal :math:`\hat e_x` — not "a reflection through the x axis",
    which in :math:`\mathbb{R}^3` would be :math:`\mathrm{diag}(1,-1,-1)`, a
    half-turn with :math:`\det = +1`. The core makes the wrong reading
    unspellable by taking ``normal`` keyword-only; this wrapper exists to
    translate the tag's axis LETTER into that normal.
    """
    return [RigidMotion.reflection(normal=_axis_vector(axis))]


def _inversion_op() -> RigidMotion:
    r"""Inversion through the origin, :math:`x \to -x`.

    :math:`\det = (-1)^3 = -1`, and it fixes only the origin — so in
    :math:`\mathbb{R}^3` it is neither a reflection (which would fix a plane)
    nor a rotation (which would fix a line). That is why the core builds it
    from its definition rather than from either family.
    """
    return RigidMotion.inversion(3)


def _is_axis_supported(nodes: np.ndarray, atol: float) -> bool:
    r"""``True`` iff every node lies on the :math:`z`-axis.

    The exact criterion for **SO(2)-invariance of a FINITE point set**.
    The :math:`SO(2)` orbit of a point at cylindrical radius
    :math:`\rho > 0` is a whole circle — an infinite set — so a finite
    set containing it cannot be closed. Points with :math:`\rho = 0`
    are fixed by every rotation about :math:`z`. Hence closure
    :math:`\iff` every node has :math:`\rho \le` ``atol``.

    Replaces the retired ``_so2_representatives()``, which sampled
    :math:`\{0, 90, 180, 270\}^\circ` and therefore tested closure under
    the finite subgroup :math:`C_4`, **not** under :math:`SO(2)`
    (ERR-072). That is unsound in the dangerous direction: it CERTIFIES
    non-invariant rules. Measured on the shipped product family, the
    old answer was a function of ``n_phi mod 4`` — ``True`` for
    ``n_phi`` in 4/8/12/16 and ``False`` for 2/3/5/6/7 — while the true
    azimuthal group is the finite :math:`D_{n_\varphi h}` in every case.

    Consequence, and it is intended: **no real angular cubature is
    SO(2)-invariant.** A consumer needing "this rule respects a
    continuous azimuthal symmetry" is asking a question about the
    rule's exactness space, not about node-set closure, and must not
    be answered here.
    """
    return bool(np.all(np.hypot(nodes[:, 0], nodes[:, 1]) <= atol))


def _is_origin_supported(nodes: np.ndarray, atol: float) -> bool:
    r"""``True`` iff every node sits at the origin.

    The exact criterion for **SO(3)- (and O(3)-) invariance of a FINITE
    point set**: every orbit of a non-origin point under :math:`SO(3)`
    is a whole 2-sphere, so only the origin can be carried by a finite
    invariant set.

    Replaces the retired icosahedral-sample check, which tested
    :math:`I_h` closure and called it :math:`SO(3)` (ERR-072). That
    conflation had a second edge: :math:`-I \in I_h`, so the
    ``SO3`` and ``O3`` branches ran the *same* operator set and were
    identically equal for every input — and 60 of those 120 matrices
    are improper, hence not in :math:`SO(3)` at all.
    """
    return bool(np.all(np.linalg.norm(nodes, axis=1) <= atol))


def _cyclic_ops(n: int) -> list[RigidMotion]:
    r""":math:`C_n`: the ``n`` proper rotations about the z-axis.

    The angles are the :math:`n`-th roots of unity, taken from
    :func:`~orpheus.numerics.roots_of_unity.roots_of_unity` as exact points
    on the circle rather than recomputed from :math:`\theta_k = 2\pi k/n`
    through :func:`numpy.cos`. Two reasons, and the second is the
    structural one:

    1. The symmetric angles have exactly-representable components
       (:math:`\cos(\pi/2) = 0`) that a rounded :math:`\pi/2` does not
       reproduce. `[M]` at :math:`n = 4` this takes the operator matrices
       from 30/36 exactly-representable entries to **36/36**, and the
       landing residual on ``product(4, 4)`` from ``2.1e-16`` to exactly
       ``0.0``.
    2. It is the SAME generator the azimuthal rules are built from
       (:func:`~orpheus.numerics.quadrature.periodic_trapezoid`). A checker
       that derives its own cosines is a second spelling of the very nodes
       it is checking — so the two sides could drift, and the tolerance
       absorbing that drift would be reported as the rule's own accuracy.
    """
    cos, sin = roots_of_unity(np.arange(n), n)
    return [
        RigidMotion.rotation_from_circle_point(
            plane=(_X_AXIS, _Y_AXIS), point=(cos[k], sin[k])
        )
        for k in range(n)
    ]


def _vertical_mirrors(n: int) -> list[RigidMotion]:
    r"""The n vertical mirror planes of :math:`D_{nh}` (containing the z-axis).

    Built as the **coset** :math:`C_n\,\sigma_0` — that is,
    :math:`\sigma_k = R(2\pi k/n)\,\sigma_0`, where :math:`\sigma_0` is the
    :math:`\varphi = 0` mirror (normal :math:`\hat e_y`). This is the group
    fact :math:`D_n = C_n \sqcup C_n\sigma` spelled directly, and it has
    three consequences worth naming:

    - The mirror set IS the rotation set carried by one reflection, not a
      second trigonometric construction of the same angles. There is one
      source of truth for the :math:`n`-fold angles in this module.
    - :math:`\sigma_0` is a coordinate mirror, hence a bit-exact signed
      diagonal, so composing with it is a column sign flip and introduces
      **no** round-off. `[M]` the mirrors' landing residual is therefore
      exactly EQUAL to the rotations' on every rule measured — the mirror
      half can no longer be the worse half, which it previously was by
      1.7--4.8x because it went through a half-angle and a normalisation.
    - The half-angle disappears. The normal at :math:`k\pi/n + \pi/2` is
      a root of unity of order :math:`4n`, i.e. a *different* generator
      from the rule's own order-:math:`n` one; the coset form needs only
      the order-:math:`n` roots the rules themselves use.

    The plane placement is the convention that matters, and it is fixed
    by the standard :math:`D_{nh}` setting: the principal axis along z
    with a vertex on the **x-axis**, so :math:`\sigma_v` through
    :math:`\varphi = 0` is a symmetry. This is the setting every
    azimuthal rule in the tree is built in —
    :func:`~orpheus.numerics.quadrature.periodic_trapezoid` puts a node
    at :math:`\varphi = 0` for the unshifted (``NODE_ALIGNED``) family.
    Composing on the RIGHT by :math:`\sigma_0` is what places the k-th
    plane at :math:`k\pi/n`: :math:`\sigma_k` sends
    :math:`\varphi \mapsto 2\pi k/n - \varphi`, whose fixed azimuths are
    :math:`k\pi/n` and :math:`k\pi/n + \pi`.

    Previously the k-th **normal** was placed at :math:`k\pi/n`, i.e.
    every plane was rotated by :math:`\pi/2`. For even ``n`` that maps
    the normal set onto itself (:math:`\pi/2` is then a multiple of
    :math:`\pi/n`) and is invisible; for **odd** ``n`` it is a
    genuinely different set of planes, and
    ``Dnh(n).is_invariant(product(4, n))`` read ``False`` for
    ``n = 3, 5, 7`` while the rule is demonstrably closed under the
    :math:`\varphi = 0` mirror at every ``n``. Orthogonality,
    determinant, closure and group order are all preserved by a rotated
    mirror set, so none of those checks can see it. Two gates can, by
    different mechanisms: the plane-placement gate in
    ``tests/numerics/test_symmetry.py`` compares the placement against
    the setting, and the mirror-accuracy gate in
    ``tests/numerics/test_symmetry_exactness.py`` sees it as a RESIDUAL
    — rotating the planes takes the mirrors off the odd-``n`` rule
    entirely, so the landing blows up from ``2.5e-16`` to ``1.97e-01``.
    `[M]` that second gate reds on odd ``n`` and stays green on even
    ``n``, which is this paragraph's own claim, measured.
    """
    sigma_0 = RigidMotion.reflection(normal=_Y_AXIS)
    return [rotation @ sigma_0 for rotation in _cyclic_ops(n)]


def _octahedral_ops() -> list[RigidMotion]:
    """Generator set for the full octahedral group :math:`O_h` (48 elements).

    The full group is generated by:

    - All sign flips on the three coordinate axes (8 elements,
      :math:`(\\pm 1)^3`).
    - All coordinate permutations of :math:`(x, y, z)` (6 elements,
      :math:`S_3`).

    The Cartesian product gives :math:`8 \\times 6 = 48` matrices —
    exactly :math:`|O_h|`. The orbit-closure check applies all 48; a
    measure invariant under all of them is :math:`O_h`-invariant.
    """
    from itertools import product, permutations
    return [
        RigidMotion.signed_permutation(permutation=perm, signs=signs)
        for signs in product([-1.0, 1.0], repeat=3)
        for perm in permutations([0, 1, 2])
    ]


def _icosahedral_ops() -> list[RigidMotion]:
    """Generator set for the icosahedral group :math:`I_h` (120 elements).

    The icosahedral group is generated by a 5-fold rotation about a
    vertex axis and a 3-fold rotation about a face axis, plus
    inversion. We construct a representative orbit by:

    1. Listing the 12 vertices of a regular icosahedron.
    2. Generating proper rotations that map the first vertex to each
       other vertex (12 cosets) — gives the 60 elements of :math:`I`.
    3. Multiplying by inversion to extend to :math:`I_h` (120 elements).

    This is the standard construction; see Hamermesh 1962 §9.4.
    """
    # Icosahedron vertices (golden-ratio convention, normalised).
    phi = (1.0 + np.sqrt(5.0)) / 2.0
    raw = np.array([
        [0,  1,  phi], [0,  1, -phi], [0, -1,  phi], [0, -1, -phi],
        [1,  phi, 0], [1, -phi, 0], [-1,  phi, 0], [-1, -phi, 0],
        [phi, 0,  1], [phi, 0, -1], [-phi, 0,  1], [-phi, 0, -1],
    ], dtype=float)
    raw /= np.linalg.norm(raw[0])

    # 5-fold rotation about the axis through raw[0] = (0, 1, phi)/|.|.
    axis5 = raw[0]
    R5 = _rotation_about_axis(axis5, 2.0 * np.pi / 5.0)

    # 3-fold rotation about a face center (choose center of face
    # containing raw[0], raw[2], raw[8]: average and normalise).
    face = (raw[0] + raw[2] + raw[8]) / 3.0
    axis3 = face / np.linalg.norm(face)
    R3 = _rotation_about_axis(axis3, 2.0 * np.pi / 3.0)

    # The group generated by the 5-fold, the 3-fold and inversion IS I_h —
    # 120 elements. The hand-rolled BFS that used to live here (with its own
    # dedup-by-allclose and its own iteration cap) was a second closure
    # implementation; the core's is the one that carries the cap as a NAMED
    # refusal.
    return list(_close_rigid_motions([R5, R3, _inversion_op()]))


def _rotation_about_axis(axis: np.ndarray, theta: float) -> RigidMotion:
    r"""Rotation about ``axis`` by ``theta``.

    Delegates to the core, which builds the rotation **plane** orthogonal to
    ``axis`` and rotates in it. "About an axis" is a three-dimensional
    convenience — the :math:`(d-2)`-dimensional fixed set happens to be a line
    — not a separate construction.
    """
    return RigidMotion.rotation_about_axis(axis=axis, angle=theta)


# ---------------------------------------------------------------------------
# Orbit-closure check, and the certificate it produces
# ---------------------------------------------------------------------------

#: The node-match window is this multiple of ``atol``, which is the WEIGHT
#: window. The asymmetry is deliberate and worth naming rather than leaving
#: as a bare ``* 100``: a node coordinate is the accumulated result of a
#: matrix product against a constructed direction cosine, so it carries more
#: round-off than a weight, which is usually read straight from a table. Both
#: windows are ABSOLUTE — for rules whose weights are O(1e-3) the weight test
#: is correspondingly stricter in relative terms, which is a known
#: characteristic of this check and not an accident of it.
_NODE_WINDOW_FACTOR = 100.0


@dataclass(frozen=True)
class OrbitCertificate:
    r"""The realized action of a group on a point set.

    Returned by the closure check instead of a bare ``True``. The check
    ALREADY computes, for each operator :math:`M` and node :math:`i`, the
    index :math:`\pi_M(i)` with :math:`M x_i = x_{\pi_M(i)}` — it simply
    discarded it. A ``-> bool`` predicate that internally builds the
    permutation IS the missing primitive; widening the return type is
    cheaper than minting a class to recompute it.

    Everything below is a *reading* of the permutations, not new work:

    * :attr:`singular_set` — the orbifold **singular set**
      :math:`\Sigma = \{x : \mathrm{Stab}(x) \neq \{e\}\}`. Because
      :math:`\pi_M(i) = i` means exactly :math:`M x_i = x_i`, i.e.
      :math:`x_i \in \mathrm{Fix}(M)`, membership is an **integer
      identity** — exact, with no tolerance. The only place a tolerance
      enters is matching nodes while BUILDING :math:`\pi`, which is the
      one place the question is honestly numerical.
    * :attr:`stabilizer_order` — :math:`|\mathrm{Stab}(x_i)|`, the orbit
      type; by orbit-stabilizer the orbit length is :math:`|G|` divided
      by it.
    * :meth:`orbits` — the :math:`G`-orbits themselves.

    A certificate exists only for a :math:`G`-INVARIANT set, so
    :math:`\Sigma` is unrepresentable without the closure proof — which is
    the precondition ("the quotient is defined only on a G-invariant
    measure") enforced by construction rather than by a comment.
    """

    operators: "tuple[RigidMotion, ...]"
    permutations: "tuple[Permutation, ...]"

    @property
    def n_points(self) -> int:
        return self.permutations[0].n if self.permutations else 0

    @property
    def _non_identity(self) -> "list[Permutation]":
        """The permutations of operators other than the identity.

        The identity fixes every point, so including it would report the
        whole set as singular. The test asks the ELEMENT whether it is the
        identity of its own dimension, rather than comparing against a
        hard-coded ``np.eye(3)`` — the certificate is a statement about a
        group acting on a point set, and nothing in it is three-dimensional.
        """
        return [
            perm
            for M, perm in zip(self.operators, self.permutations)
            if not M.approximately_equals(RigidMotion.identity(M.dimension))
        ]

    @property
    def singular_set(self) -> np.ndarray:
        r"""Indices of the **singular points** — :math:`\Sigma`.

        A point is singular iff some non-identity group element fixes it.
        Under a reflection the fixed locus is a **mirror**; a point fixed
        by two mirrors is a corner reflector; a point on a rotation axis
        with no mirror is a cone point.
        """
        fixed = np.zeros(self.n_points, dtype=bool)
        for perm in self._non_identity:
            fixed[perm.fixed_points] = True
        return np.flatnonzero(fixed)

    @property
    def stabilizer_order(self) -> np.ndarray:
        r""":math:`|\mathrm{Stab}(x_i)|` for every node ``i``.

        ``1`` exactly on the regular (free) points; the singular set is
        ``> 1``. Meaningful only when the certificate was built from the
        FULL group rather than a generating set.
        """
        order = np.zeros(self.n_points, dtype=np.int64)
        for perm in self.permutations:
            order[perm.fixed_points] += 1
        return order

    def orbits(self) -> "tuple[np.ndarray, ...]":
        """The :math:`G`-orbits, as arrays of node indices."""
        parent = list(range(self.n_points))

        def find(a: int) -> int:
            while parent[a] != a:
                parent[a] = parent[parent[a]]
                a = parent[a]
            return a

        for perm in self.permutations:
            for i, j in enumerate(perm.indices):
                ra, rb = find(i), find(int(j))
                if ra != rb:
                    parent[ra] = rb

        buckets: "dict[int, list[int]]" = {}
        for i in range(self.n_points):
            buckets.setdefault(find(i), []).append(i)
        return tuple(np.array(v, dtype=np.int64) for v in buckets.values())


def _orbit_closure(
    nodes: np.ndarray,
    weights: np.ndarray,
    ops: Iterable[RigidMotion],
    atol: float,
) -> "OrbitCertificate | None":
    """Check that applying every operator ``M ∈ ops`` to ``nodes``
    yields the same multiset of (node, weight) pairs.

    The per-element work is
    :meth:`~orpheus.geometry.transformation.RigidMotion.preserves`; what
    remains here is the *measure*-level question — that EVERY element
    preserves it, and the certificate assembled from the results. The
    matching, the bijectivity requirement (ERR-073) and the weight guard all
    live in the core, so this module no longer carries a second copy of them.

    The two windows: ``atol`` is the WEIGHT tolerance, and the node-match
    window is :data:`_NODE_WINDOW_FACTOR` times larger. That relationship is
    a claim about *measures* — a coordinate accumulates matrix-product
    round-off while a weight is usually read from a table — so it is stated
    here rather than inside the geometry primitive, which takes the two
    windows separately and refuses to guess.

    For each :math:`M`, compute :math:`M(\\text{nodes})` and find a
    permutation :math:`\\pi` such that
    :math:`M(\\text{nodes})_i = \\text{nodes}_{\\pi(i)}` (within
    ``atol``). The match is verified element-wise, then the weight
    equality :math:`w_i = w_{\\pi(i)}` is checked.

    The match is required to be a **bijection**. Nearest-neighbour
    matching alone only proves every image has *a* same-weight partner,
    which is strictly weaker than "the action permutes the nodes": two
    distinct sources may land on one target, leaving some node
    unmatched entirely. Such a set is NOT :math:`G`-invariant — its
    mass is distributed differently from its image — yet it passes an
    injectivity-free check. Appending a bit-identical duplicate of any
    node to an :math:`O_h`-invariant rule is the minimal witness: the
    duplicated position then carries twice the mass of its mirror
    image, and every one of the 48 match maps is non-injective.
    Since :math:`\\pi` maps an :math:`n`-set to an :math:`n`-set,
    injectivity is equivalent to bijectivity, so counting distinct
    targets suffices. (ERR-073.)

    Returns ``None`` at the first failure, else the
    :class:`OrbitCertificate` carrying every :math:`\\pi_M`.
    """
    kept_ops: "list[RigidMotion]" = []
    perms: "list[Permutation]" = []

    for M in ops:
        pi = M.preserves(
            nodes,
            weights,
            atol=atol * _NODE_WINDOW_FACTOR,
            weight_atol=atol,
        )
        if pi is None:
            return None
        kept_ops.append(M)
        perms.append(pi)
    return OrbitCertificate(operators=tuple(kept_ops), permutations=tuple(perms))


# ---------------------------------------------------------------------------
# Public attribute install
# ---------------------------------------------------------------------------

# Pre-instantiated singletons attached to the public class. These are
# the canonical way users obtain a SubgroupOfO3 for a named entry.
SubgroupOfO3.Trivial = SubgroupOfO3._from_named(_NamedSubgroup.Trivial)
SubgroupOfO3.SO2 = SubgroupOfO3._from_named(_NamedSubgroup.SO2)
SubgroupOfO3.Dinfh = SubgroupOfO3._from_named(_NamedSubgroup.Dinfh)
SubgroupOfO3.OctahedralOh = SubgroupOfO3._from_named(_NamedSubgroup.OctahedralOh)
SubgroupOfO3.IcosahedralIh = SubgroupOfO3._from_named(_NamedSubgroup.IcosahedralIh)
SubgroupOfO3.SO3 = SubgroupOfO3._from_named(_NamedSubgroup.SO3)
SubgroupOfO3.O3 = SubgroupOfO3._from_named(_NamedSubgroup.O3)


# ---------------------------------------------------------------------------
# The subgroup graph, and walking it to find a measure's symmetry
# ---------------------------------------------------------------------------
#
# Crystallography does not ASK a structure what its symmetry group is — it
# walks the subgroup lattice downward from high symmetry until it reaches
# the symmetry the structure actually has. The graph is the Hasse diagram of
# the lattice: nodes are groups, edges are MAXIMAL-subgroup relations
# (International Tables Vol. A1 renders it as a Bärnighausen tree).
#
# That is the machinery this module needs. A DECLARED invariance group is a
# claim with no construction behind it, and such claims have already shipped
# false here twice. A COMPUTED one cannot lie about the object it was
# computed from.


def _distinct_azimuths(nodes: np.ndarray, atol: float) -> int:
    r"""Number of distinct azimuthal angles among the off-axis nodes.

    Bounds the cyclic families: a :math:`C_n` rotation (``n > 1``) fixes no
    azimuth, so the azimuths split into FREE orbits of size ``n`` and
    therefore ``n`` divides this count. That turns two infinite families
    into a handful of divisors.
    """
    rho = np.hypot(nodes[:, 0], nodes[:, 1])
    off_axis = nodes[rho > atol]
    if off_axis.shape[0] == 0:
        return 0
    phi = np.sort(np.mod(np.arctan2(off_axis[:, 1], off_axis[:, 0]), 2.0 * np.pi))
    distinct = [phi[0]]
    for p in phi[1:]:
        if p - distinct[-1] > 1e-9:
            distinct.append(float(p))
    # The first and last can be the same angle seen either side of the branch.
    if len(distinct) > 1 and (2.0 * np.pi - distinct[-1] + distinct[0]) < 1e-9:
        distinct.pop()
    return len(distinct)


def candidate_groups(
    measure: DiscreteMeasure, *, atol: float = 1e-13,
) -> "tuple[SubgroupOfO3, ...]":
    """The expressible groups worth testing against ``measure``.

    The named entries always, plus :math:`C_n` / :math:`D_{nh}` for each
    divisor of the measure's distinct-azimuth count (see
    :func:`_distinct_azimuths` for why divisors suffice).
    """
    named = [
        SubgroupOfO3.Trivial, SubgroupOfO3.SO2,
        SubgroupOfO3.Dinfh, SubgroupOfO3.OctahedralOh,
        SubgroupOfO3.IcosahedralIh, SubgroupOfO3.SO3, SubgroupOfO3.O3,
    ]
    # All three mirrors, always. The parameter set of a reflection family
    # is {x, y, z} — FINITE BY CONSTRUCTION, so unlike Cn/Dnh (which need
    # the distinct-azimuth divisor bound to stay finite) it needs no bound
    # at all, and unlike the retired parameter-free Z2 it offers the walk
    # all three planes instead of silently only sigma_z.
    named += [SubgroupOfO3.Mirror(a) for a in ("x", "y", "z")]
    nodes = measure.nodes
    if nodes.ndim == 1 or nodes.shape[1] < 3:
        return tuple(named)
    n_az = _distinct_azimuths(nodes, atol)
    families: "list[SubgroupOfO3]" = []
    for d in (d for d in range(1, n_az + 1) if n_az % d == 0) if n_az else (1,):
        families.append(SubgroupOfO3.Cn(d))
        families.append(SubgroupOfO3.Dnh(d))
    return tuple(named + families)


def _maximal(groups: "Iterable[SubgroupOfO3]") -> "tuple[SubgroupOfO3, ...]":
    """Those members not strictly contained in another member."""
    items = list(groups)
    return tuple(
        g for g in items
        if not any(h != g and h.contains(g) for h in items)
    )


def maximal_subgroups(
    group: "SubgroupOfO3", candidates: "Iterable[SubgroupOfO3]",
) -> "tuple[SubgroupOfO3, ...]":
    """The Hasse edges below ``group`` — its maximal proper subgroups.

    Derived from :meth:`SubgroupOfO3.contains`, never tabulated: a
    hand-drawn Bärnighausen tree would be exactly the unverifiable claim
    the computed lattice exists to eliminate.
    """
    below = [h for h in candidates if h != group and group.contains(h)]
    return _maximal(below)


def maximal_invariance_groups(
    measure: DiscreteMeasure,
    *,
    atol: float = 1e-13,
    candidates: "Iterable[SubgroupOfO3] | None" = None,
    method: str = "walk",
) -> "tuple[SubgroupOfO3, ...]":
    r"""The maximal groups leaving ``measure`` invariant — its symmetry.

    Returns the maximal elements of
    :math:`\{G \in C : \mu \text{ is } G\text{-invariant}\}`. Several
    incomparable maxima can survive (a rule may carry both a rotation group
    and a reflection with neither inside the other); the true
    :math:`\mathrm{Sym}(\mu)` is the group they GENERATE, which need not be
    expressible in the candidate set — so the maximal elements, not a single
    tag, are the honest answer.

    **Soundness.** Invariance is DOWNWARD-CLOSED: if :math:`\mu` is
    :math:`G`-invariant and :math:`H \le G`, then :math:`\mu` is
    :math:`H`-invariant. So the invariant set is an order ideal, giving the
    walk both prunings — a failing node kills every supergroup, a passing
    node implies every subgroup. The walk therefore *requires* a correct
    lattice and is simultaneously a test of it.

    Parameters
    ----------
    method : {"walk", "bruteforce"}
        Two realizations of the same definition. ``"walk"`` descends the
        Hasse diagram with both prunings; ``"bruteforce"`` tests every
        candidate. They must agree on every input — that agreement is the
        verification instrument, not a mere fast-path check.

    Notes
    -----
    The answer is about the group's realization in the **standard setting**
    (principal axis along z), not up to conjugation. A rule whose symmetry
    axis is not z reports a smaller group, which is correct for a gate
    comparing against a geometry in the same frame.
    """
    cands = tuple(candidates) if candidates is not None else candidate_groups(
        measure, atol=atol
    )

    if method == "bruteforce":
        return _maximal(
            g for g in cands if g.is_invariant(measure, atol=atol)
        )
    if method != "walk":
        raise ValueError(
            f"method must be 'walk' or 'bruteforce', got {method!r}"
        )

    accepted: "list[SubgroupOfO3]" = []
    visited: set[str] = set()
    stack = list(_maximal(cands))
    while stack:
        group = stack.pop()
        key = repr(group)
        if key in visited:
            continue
        visited.add(key)
        if any(a.contains(group) for a in accepted):
            continue  # already implied by an accepted supergroup
        if group.is_invariant(measure, atol=atol):
            accepted.append(group)  # everything below is implied — do not descend
        else:
            stack.extend(maximal_subgroups(group, cands))
    return _maximal(accepted)


def orbit_certificate(
    measure: DiscreteMeasure,
    group: "SubgroupOfO3",
    *,
    atol: float = 1e-13,
) -> "OrbitCertificate | None":
    r"""The realized action of ``group`` on ``measure``, or ``None``.

    ``None`` when the measure is not ``group``-invariant, or when the group
    is CONTINUOUS and so has no finite realization to permute nodes by. The
    certificate is built from the group's ELEMENTS, not a generating set:
    orbit closure only needs generators, but a point's stabilizer may be
    generated by a non-generator, so :math:`\Sigma` needs them all.
    """
    ops = _realized_ops(group._tag)
    if ops is None:
        return None
    nodes = measure.nodes
    if nodes.ndim == 1 or nodes.shape[1] < 3:
        return None
    return _orbit_closure(nodes, measure.weights, _close_group(ops), atol)


def singular_set(
    measure: DiscreteMeasure,
    group: "SubgroupOfO3",
    *,
    atol: float = 1e-13,
) -> np.ndarray:
    r"""The **singular set** :math:`\Sigma` of ``measure`` under ``group``.

    Indices of the **singular points** — those whose stabilizer is
    non-trivial, i.e. fixed by some non-identity group element. Under a
    reflection the fixed locus is a **mirror**; a point on two mirrors is a
    corner reflector; a point on a rotation axis with no mirror through it
    is a cone point.

    Membership is decided by :math:`\pi_M(i) = i` — an **integer identity**,
    exact, with no tolerance. The three ad-hoc float comparisons the tree
    grew for this question (``_OCTANT_SIGN_EPS``, ``_MU_DIRECTION_EPS``,
    ``TANGENTIAL_EPS``) were all asking it numerically; measured across 29
    production rules, the separation between "zero" and "nonzero" cosines is
    a factor of :math:`2.7\times10^{13}`, so the tolerance was never doing
    real work. The only honestly-numerical step is matching nodes while
    BUILDING :math:`\pi`, which is ``atol``.

    Raises
    ------
    ValueError
        If ``measure`` is not ``group``-invariant. :math:`\Sigma` is defined
        only on an invariant set — a quotient needs something to quotient —
        so the illegal state is unrepresentable rather than silently
        returning a wrong answer.
    """
    cert = orbit_certificate(measure, group, atol=atol)
    if cert is None:
        raise ValueError(
            f"the singular set is defined only for a {group.name}-invariant "
            f"measure with a finite realization; this measure is not "
            f"{group.name}-invariant (or {group.name} is continuous, and a "
            f"continuous group has no finite node permutation)"
        )
    return cert.singular_set


__all__ = [
    "OrbitCertificate",
    "SubgroupOfO3",
    "candidate_groups",
    "maximal_invariance_groups",
    "maximal_subgroups",
    "orbit_certificate",
    "singular_set",
]
