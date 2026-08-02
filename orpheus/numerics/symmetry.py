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

from .measure import DiscreteMeasure, SPACE_INTERVAL_M11, SPACE_SPHERE


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
    Z2 = "Z2"                  # {e, σ} — single reflection / 180° rotation
    SO2 = "SO2"                # axial rotations (continuous 1-parameter)
    O2 = "O2"                  # SO(2) ⋊ Z_2 — axial rotations + reflection
    OctahedralOh = "Oh"        # full octahedral group, 48 elements
    IcosahedralIh = "Ih"       # full icosahedral group, 120 elements
    SO3 = "SO3"                # proper rotations of the sphere
    O3 = "O3"                  # orthogonal group of R^3 (top of lattice)


@dataclass(frozen=True)
class Cn:
    """Cyclic group :math:`C_n` of :math:`n`-fold proper rotations
    about a fixed axis.

    Parameters
    ----------
    n : int
        Order of the cyclic group; must be :math:`\\ge 1`. ``Cn(1)``
        is the trivial group; ``Cn(2)`` is the rotation-only sibling
        of :class:`_NamedSubgroup.Z2` (``Z2`` here is *any* order-2
        subgroup including pure reflections; ``Cn(2)`` is specifically
        the :math:`180^{\\circ}`-rotation group).
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


# Public type for "any subgroup tag". The selection logic below
# branches on type — Cn and Dnh are kept as separate dataclasses (not
# enum entries) so the order parameter is explicit.
SubgroupTag = "_NamedSubgroup | Cn | Dnh"


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
    (_NamedSubgroup.Trivial, _NamedSubgroup.Z2),
    (_NamedSubgroup.Trivial, _NamedSubgroup.SO2),
    (_NamedSubgroup.Trivial, _NamedSubgroup.O2),
    (_NamedSubgroup.Trivial, _NamedSubgroup.OctahedralOh),
    (_NamedSubgroup.Trivial, _NamedSubgroup.IcosahedralIh),
    (_NamedSubgroup.Trivial, _NamedSubgroup.SO3),
    (_NamedSubgroup.Trivial, _NamedSubgroup.O3),

    # Z2 ⊂ everything except SO2 (a continuous proper-rotation group
    # has no order-2 reflection element).
    (_NamedSubgroup.Z2, _NamedSubgroup.O2),
    (_NamedSubgroup.Z2, _NamedSubgroup.OctahedralOh),
    (_NamedSubgroup.Z2, _NamedSubgroup.IcosahedralIh),
    (_NamedSubgroup.Z2, _NamedSubgroup.SO3),
    (_NamedSubgroup.Z2, _NamedSubgroup.O3),

    # SO2 ⊂ O2, SO3, O3.
    (_NamedSubgroup.SO2, _NamedSubgroup.O2),
    (_NamedSubgroup.SO2, _NamedSubgroup.SO3),
    (_NamedSubgroup.SO2, _NamedSubgroup.O3),

    # O2 ⊂ O3 (the projection R^2 → R^3 along a fixed axis).
    (_NamedSubgroup.O2, _NamedSubgroup.O3),

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
        SubgroupOfO3.Z2
        SubgroupOfO3.SO2
        SubgroupOfO3.O2
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
    Z2: ClassVar[SubgroupOfO3]
    SO2: ClassVar[SubgroupOfO3]
    O2: ClassVar[SubgroupOfO3]
    OctahedralOh: ClassVar[SubgroupOfO3]
    IcosahedralIh: ClassVar[SubgroupOfO3]
    SO3: ClassVar[SubgroupOfO3]
    O3: ClassVar[SubgroupOfO3]

    def __init__(self, tag: "_NamedSubgroup | Cn | Dnh") -> None:
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
          :math:`D_{nh} \subseteq O(2)` for every :math:`n`.

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
        - **Z2 / Cn / SO2 / O2 / SO3 / O3 / Dnh**: handled via the
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


def _contains(outer, inner) -> bool:
    """Dispatch table for :meth:`SubgroupOfO3.contains`.

    Branches on the (outer, inner) tag types and falls back to the
    static named-lattice when both are named entries.
    """
    # Reflexivity for parameterised families.
    if outer == inner:
        return True

    # Both named — use the static lattice.
    if isinstance(outer, _NamedSubgroup) and isinstance(inner, _NamedSubgroup):
        return _named_contains(outer, inner)

    # ----- Inner is named, outer is parameterised --------------------
    if isinstance(outer, Cn) and isinstance(inner, _NamedSubgroup):
        # The only named subgroup of an arbitrary Cn is the trivial group.
        # (Cn(2) coincides with Z2 in the proper-rotation sense, but our
        # Z2 entry is the abstract order-2 group — we keep them distinct
        # so the named lattice stays unambiguous.)
        return inner is _NamedSubgroup.Trivial

    if isinstance(outer, Dnh) and isinstance(inner, _NamedSubgroup):
        # D_nh contains the trivial group, the Z_2 mirror reflection,
        # and (for n >= 1) the principal-axis cyclic subgroup C_n.
        # Higher named subgroups (SO2, O2, etc.) are not contained in a
        # finite dihedral group.
        if inner is _NamedSubgroup.Trivial or inner is _NamedSubgroup.Z2:
            return True
        return False

    # ----- Outer is named, inner is parameterised --------------------
    if isinstance(outer, _NamedSubgroup) and isinstance(inner, Cn):
        # Cn ⊂ SO2 for every n (SO(2) = union over n of C_n).
        # Cn ⊂ O2, SO3, O3 transitively. Cn ⊂ Oh / Ih only for specific
        # n that match the polyhedral rotation axes — out of scope for
        # the static lattice (Cn ⊂ Oh would need n ∈ {1,2,3,4,6}; we do
        # not encode this until a consumer needs it).
        if outer in (
            _NamedSubgroup.SO2,
            _NamedSubgroup.O2,
            _NamedSubgroup.SO3,
            _NamedSubgroup.O3,
        ):
            return True
        if outer is _NamedSubgroup.Trivial:
            return inner.n == 1
        return False

    if isinstance(outer, _NamedSubgroup) and isinstance(inner, Dnh):
        # Dnh ⊂ O2, O3 always. Not ⊂ SO2 / SO3 (Dnh contains improper
        # rotations).
        return outer in (_NamedSubgroup.O2, _NamedSubgroup.O3)

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


def _check_invariance(tag, measure: DiscreteMeasure, atol: float) -> bool:
    """Dispatch table for :meth:`SubgroupOfO3.is_invariant`."""
    # Trivial is always invariant.
    if isinstance(tag, _NamedSubgroup) and tag is _NamedSubgroup.Trivial:
        return True
    if isinstance(tag, Cn) and tag.n == 1:
        return True

    # Determine which "ambient space" the measure lives in. For 1-D
    # measures the rotation/reflection axis story is degenerate and
    # most groups are trivially invariant; for S^2 (or any 3-D node
    # array) the check is non-trivial.
    is_1d = measure.dim == 1 or measure.support in (
        SPACE_INTERVAL_M11,
        "[0,1]",
        "R",
    )

    if is_1d:
        return _check_invariance_1d(tag, measure, atol)

    # 3-D node arrays (S^2, R^3, …). For SN-relevant rules the nodes
    # lie on the unit sphere.
    nodes = measure.nodes
    if nodes.ndim == 1:
        return False  # 1-D nodes — handled above; safety net
    if nodes.shape[1] == 2:
        # 2-D ambient — promote to 3-D by appending z=0 so the same
        # transformation kernels work. (Useful for product quadrature
        # restricted to a plane, not consumed today but harmless.)
        nodes_full = np.column_stack([nodes, np.zeros(nodes.shape[0])])
    else:
        nodes_full = nodes
    return _check_invariance_3d(tag, nodes_full, measure.weights, atol)


def _check_invariance_1d(tag, measure: DiscreteMeasure, atol: float) -> bool:
    """Invariance check for 1-D measures.

    Most rotational groups act trivially on a 1-D measure (there is no
    azimuthal coordinate). The non-trivial check is reflection
    :math:`x \\to -x`, which is what Z2 / O2 / O3 require on
    :math:`[-1, 1]`.
    """
    # SO(2), SO(3), Cn, Dnh (with n >= 1): trivially invariant on a 1-D
    # measure — the rotation axis is the (missing) angular coordinate.
    rotational_only = (
        (isinstance(tag, _NamedSubgroup) and tag in (
            _NamedSubgroup.SO2, _NamedSubgroup.SO3,
        ))
        or isinstance(tag, Cn)
    )
    if rotational_only:
        return True

    # Z2, O2, O3, Oh, Ih, Dnh: require closure under x -> -x for a
    # measure on a symmetric interval. (For an asymmetric interval
    # the check is automatically False, which is the correct answer.)
    if isinstance(tag, _NamedSubgroup) and tag in (
        _NamedSubgroup.Z2,
        _NamedSubgroup.O2,
        _NamedSubgroup.O3,
        _NamedSubgroup.OctahedralOh,
        _NamedSubgroup.IcosahedralIh,
    ):
        return _is_reflection_invariant_1d(measure, atol)

    if isinstance(tag, Dnh):
        return _is_reflection_invariant_1d(measure, atol)

    return False


def _is_reflection_invariant_1d(measure: DiscreteMeasure, atol: float) -> bool:
    """Closure of a 1-D measure under :math:`x \\to -x`."""
    nodes = measure.nodes if measure.nodes.ndim == 1 else measure.nodes[:, 0]
    weights = measure.weights
    # For each node x, find the reflected partner -x and check weights match.
    sorted_idx = np.argsort(nodes)
    sorted_nodes = nodes[sorted_idx]
    sorted_weights = weights[sorted_idx]

    n = len(sorted_nodes)
    for i in range(n):
        # Find -sorted_nodes[i] in sorted_nodes via binary search.
        target = -sorted_nodes[i]
        # np.searchsorted gives an insertion index; use it to find the
        # nearest node and check distance.
        j = np.searchsorted(sorted_nodes, target)
        # Try both j and j-1 (insertion point may be on either side).
        candidates = []
        if 0 <= j < n:
            candidates.append(j)
        if 0 <= j - 1 < n:
            candidates.append(j - 1)
        if not candidates:
            return False
        best = min(candidates, key=lambda k: abs(sorted_nodes[k] - target))
        if abs(sorted_nodes[best] - target) > atol * 10:
            return False
        if abs(sorted_weights[best] - sorted_weights[i]) > atol:
            return False
    return True


def _check_invariance_3d(
    tag,
    nodes: np.ndarray,
    weights: np.ndarray,
    atol: float,
) -> bool:
    """Invariance check for 3-D-node measures (typically on :math:`S^2`)."""
    if isinstance(tag, _NamedSubgroup):
        if tag is _NamedSubgroup.Z2:
            # Single reflection — pick z -> -z as a canonical
            # representative. (Any single reflection works; the choice
            # is convention.)
            return _orbit_closure(nodes, weights, _reflections("z"), atol)

        if tag is _NamedSubgroup.SO2:
            # DECIDED EXACTLY, never sampled (ERR-072).  A continuous
            # group cannot be tested by a finite sample: the sample
            # generates a finite SUBgroup, and closure under that
            # subgroup is strictly weaker than closure under G.
            return _is_axis_supported(nodes, atol)

        if tag is _NamedSubgroup.O2:
            # O(2) = SO(2) x <sigma_h>.  Both conjuncts, the first
            # exact and the second an honest finite-orbit check.
            return _is_axis_supported(nodes, atol) and _orbit_closure(
                nodes, weights, _reflections("z"), atol
            )

        if tag is _NamedSubgroup.OctahedralOh:
            return _orbit_closure(nodes, weights, _octahedral_ops(), atol)

        if tag is _NamedSubgroup.IcosahedralIh:
            return _orbit_closure(nodes, weights, _icosahedral_ops(), atol)

        if tag is _NamedSubgroup.SO3:
            # DECIDED EXACTLY (ERR-072).  Every SO(3) orbit of a
            # non-origin point is a whole 2-sphere, so a FINITE set is
            # SO(3)-closed iff it is supported at the origin.
            return _is_origin_supported(nodes, atol)

        if tag is _NamedSubgroup.O3:
            # O(3) contains SO(3), so the same exact criterion binds.
            return _is_origin_supported(nodes, atol)

    if isinstance(tag, Cn):
        # Cyclic group about z, n proper rotations.
        return _orbit_closure(
            nodes, weights, _cyclic_ops(tag.n), atol,
        )

    if isinstance(tag, Dnh):
        # Dihedral D_nh: C_n + horizontal mirror + n vertical mirrors.
        ops = _cyclic_ops(tag.n) + _reflections("z") + _vertical_mirrors(tag.n)
        return _orbit_closure(nodes, weights, ops, atol)

    return False


# ---------------------------------------------------------------------------
# Group operation generators
# ---------------------------------------------------------------------------
#
# Each operator is a (3, 3) numpy matrix M; the action on a node array
# of shape (N, 3) is ``nodes @ M.T``. The orbit-closure check applies
# every M and verifies the resulting node set matches the input.


def _rotation_z(theta: float) -> np.ndarray:
    """Rotation matrix about the z-axis by angle :math:`\\theta`."""
    c, s = np.cos(theta), np.sin(theta)
    return np.array([
        [c, -s, 0.0],
        [s,  c, 0.0],
        [0.0, 0.0, 1.0],
    ])


def _rotation_x(theta: float) -> np.ndarray:
    c, s = np.cos(theta), np.sin(theta)
    return np.array([
        [1.0, 0.0, 0.0],
        [0.0, c, -s],
        [0.0, s, c],
    ])


def _rotation_y(theta: float) -> np.ndarray:
    c, s = np.cos(theta), np.sin(theta)
    return np.array([
        [c, 0.0, s],
        [0.0, 1.0, 0.0],
        [-s, 0.0, c],
    ])


def _reflections(axis: str) -> list[np.ndarray]:
    """Single-axis reflection matrices."""
    M = np.eye(3)
    if axis == "x":
        M[0, 0] = -1.0
    elif axis == "y":
        M[1, 1] = -1.0
    elif axis == "z":
        M[2, 2] = -1.0
    else:
        raise ValueError(f"axis must be x/y/z, got {axis!r}")
    return [M]


def _inversion_op() -> np.ndarray:
    """Inversion through the origin: :math:`x \\to -x`."""
    return -np.eye(3)


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


def _cyclic_ops(n: int) -> list[np.ndarray]:
    """C_n: n proper rotations about the z-axis."""
    return [_rotation_z(2.0 * np.pi * k / n) for k in range(n)]


def _vertical_mirrors(n: int) -> list[np.ndarray]:
    """n vertical mirror planes (passing through the z-axis).

    The k-th mirror normal makes angle :math:`k\\pi/n` with the x-axis.
    Reflection through the plane is :math:`x \\to x - 2 (x \\cdot n) n`.
    """
    out = []
    for k in range(n):
        theta = k * np.pi / n
        n_vec = np.array([np.cos(theta), np.sin(theta), 0.0])
        M = np.eye(3) - 2.0 * np.outer(n_vec, n_vec)
        out.append(M)
    return out


def _octahedral_ops() -> list[np.ndarray]:
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
    ops = []
    for signs in product([-1.0, 1.0], repeat=3):
        for perm in permutations([0, 1, 2]):
            M = np.zeros((3, 3))
            for i, p in enumerate(perm):
                M[i, p] = signs[i]
            ops.append(M)
    return ops


def _icosahedral_ops() -> list[np.ndarray]:
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

    # Generate the icosahedral rotation group I (60 elements) by
    # closure under R5 and R3.
    group: list[np.ndarray] = [np.eye(3)]
    seen: list[np.ndarray] = [np.eye(3)]
    frontier = [np.eye(3)]
    generators = [R5, R3]
    # BFS over the group, with deduplication via close-comparison.
    iteration_cap = 200  # safety: I has 60 elements
    for _ in range(iteration_cap):
        new_frontier = []
        for M in frontier:
            for g in generators:
                cand = g @ M
                if not any(np.allclose(cand, s, atol=1e-10) for s in seen):
                    seen.append(cand)
                    group.append(cand)
                    new_frontier.append(cand)
        if not new_frontier:
            break
        frontier = new_frontier

    # Extend to I_h by multiplying by inversion.
    inv = _inversion_op()
    full = list(group) + [inv @ M for M in group]
    return full


def _rotation_about_axis(axis: np.ndarray, theta: float) -> np.ndarray:
    """Rodrigues rotation matrix about a unit ``axis`` by angle ``theta``."""
    a = axis / np.linalg.norm(axis)
    c, s = np.cos(theta), np.sin(theta)
    K = np.array([
        [0.0, -a[2], a[1]],
        [a[2], 0.0, -a[0]],
        [-a[1], a[0], 0.0],
    ])
    return np.eye(3) * c + s * K + (1.0 - c) * np.outer(a, a)


# ---------------------------------------------------------------------------
# Orbit-closure check
# ---------------------------------------------------------------------------


def _orbit_closure(
    nodes: np.ndarray,
    weights: np.ndarray,
    ops: Iterable[np.ndarray],
    atol: float,
) -> bool:
    """Check that applying every operator ``M ∈ ops`` to ``nodes``
    yields the same multiset of (node, weight) pairs.

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
    injectivity is equivalent to bijectivity, so tracking claimed
    targets suffices. (ERR-073.)

    Returns ``False`` at the first failure.
    """
    n = nodes.shape[0]
    # Pre-build a KD-tree-like index via lexicographic sort on rounded
    # coordinates. For n in the 100s (Lebedev-17 has 110 nodes,
    # LS_4 has 24, …) an O(n^2) match is fine; we use a binary-search
    # on a sorted view for clarity.
    rounded = np.round(nodes / atol).astype(np.int64) * atol  # quantise
    # Build a dict keyed on the bytes of the quantised triple.
    idx_of: dict[bytes, list[int]] = {}
    for i in range(n):
        key = rounded[i].tobytes()
        idx_of.setdefault(key, []).append(i)

    for M in ops:
        moved = nodes @ M.T
        # ``claimed_by[j] = i`` once source i has matched target j. A
        # second claim on the same j means the match is not injective,
        # hence not a permutation — see the docstring (ERR-073).
        claimed_by = np.full(n, -1, dtype=np.int64)
        for i in range(n):
            target = moved[i]
            target_q = (np.round(target / atol).astype(np.int64) * atol).tobytes()
            cands = idx_of.get(target_q, [])
            if not cands:
                # Fall back to a brute-force scan within tolerance —
                # quantisation can put a near-grid-line node into the
                # wrong bucket.
                dists = np.linalg.norm(nodes - target, axis=1)
                j = int(np.argmin(dists))
                if dists[j] > atol * 100:
                    return False
            else:
                # Pick the closest candidate (handles bucket collisions).
                j = min(cands, key=lambda k: float(np.linalg.norm(nodes[k] - target)))
                if np.linalg.norm(nodes[j] - target) > atol * 100:
                    return False
            if claimed_by[j] != -1:
                return False
            claimed_by[j] = i
            if abs(weights[j] - weights[i]) > atol:
                return False
    return True


# ---------------------------------------------------------------------------
# Public attribute install
# ---------------------------------------------------------------------------

# Pre-instantiated singletons attached to the public class. These are
# the canonical way users obtain a SubgroupOfO3 for a named entry.
SubgroupOfO3.Trivial = SubgroupOfO3._from_named(_NamedSubgroup.Trivial)
SubgroupOfO3.Z2 = SubgroupOfO3._from_named(_NamedSubgroup.Z2)
SubgroupOfO3.SO2 = SubgroupOfO3._from_named(_NamedSubgroup.SO2)
SubgroupOfO3.O2 = SubgroupOfO3._from_named(_NamedSubgroup.O2)
SubgroupOfO3.OctahedralOh = SubgroupOfO3._from_named(_NamedSubgroup.OctahedralOh)
SubgroupOfO3.IcosahedralIh = SubgroupOfO3._from_named(_NamedSubgroup.IcosahedralIh)
SubgroupOfO3.SO3 = SubgroupOfO3._from_named(_NamedSubgroup.SO3)
SubgroupOfO3.O3 = SubgroupOfO3._from_named(_NamedSubgroup.O3)


__all__ = [
    "SubgroupOfO3",
]
