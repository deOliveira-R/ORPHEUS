r"""How well a method realizes a declaration: the admission lattice (``admission.py`` of the conditions package).

**The mathematics.** For a method :math:`M` with trace discretization :math:`\Pi_M`, realization is a PARTIAL map
:math:`\mathrm{real}_M : (\text{declared}, \text{site}, \text{method space}) \to B_f : \Gamma_{+,h}(f') \to \Gamma_{-,h}(f)`;
*admission* answers "is :math:`\mathrm{real}_M` defined at this value on this face, and how good is it". The vocabulary is
the theory page's (``docs/theory/foundations/boundary_conditions.rst``, ``bc-method-realizability``, ``bc-refusal-axes``):

* **tier 1**, exact and faithful: :math:`R = \alpha I` (a scalar commutes with every projection):
  ``Realized(faithful=True)``;
* **tier 2**, exact but not faithful: the projection identifies laws that differ upstairs (at P1 specular and Lambertian
  both give :math:`J^- = \alpha J^+`): ``Realized(faithful=False)``;
* **tier 3**, not exact: the law depends on structure the method does not represent: ``Approximated(names, issues)``
  (also a named DOMAIN approximation, #547);
* **refused**, on exactly one of three axes (:class:`RefusalAxis`: angular resolution, topology, state cone) or on the
  fourth thing that is not an axis, ``PLUMBING`` (a missing registry entry or an unwired channel, with no theory behind
  it). "Not built yet" is NOT an axis: a refusal with no theory carries ``PLUMBING`` and the issue that tracks the
  missing machinery, so :class:`Refused` REFUSES to be built as ``PLUMBING`` without one.

**The lattice, defined once.** ``Realized(True) > Realized(False) > Approximated > Refused``. ``Admission.meet(*members)``
is the only spelling of the meet: the lesser member; the meet of two :class:`Approximated` is the union of their
``(name, issue)`` pairs (a SET, because ``meet`` on a total order would drop one, and a pair, because which issue tracks
which approximation is the information a union of two separate sets would lose); the meet of two :class:`Refused` is the first argument's (documented, deterministic); the meet of NO members is
``Realized(True)`` (an empty network and no source need no realizer). A wall's admission is the meet over its nonzero
chains and its source.

**Why a class and not a ``TypeAlias``.** The specification writes ``Admission: TypeAlias = Realized | Approximated |
Refused`` AND ``Admission.meet``; a ``typing.Union`` carries no attributes, so the two statements cannot both hold. The
sum is therefore a sealed base class (:class:`Admission`) with exactly three subclasses: ``isinstance``, ``match``
and annotations read as for the union, and ``meet`` lives on the type.

**Where a declaration is placed.** :class:`FaceSite` is derived from the carrier, computable before a method space
exists: the ``face``, the ``coord`` system and the ``dimension`` (whether the face is the hollow inner face of a radial
mesh is derived from the first two). Its :attr:`~FaceSite.site_class` is one of the three GEOMETRIC classes of the admission table
(:class:`SiteClass`: Cartesian, radial outer, radial inner), derived once; the pin and the realizers read this name.
CP's slab-left-face rule is CP's own predicate on ``site.face`` (decision 56).

**Gotchas.**

* There are two levels: the TYPE level (the front end's drop-down; the upper bound) and the INSTANCE level (per face
  and per value, plus the quadrature's equivariance certified by ``ordinate_pairing``). Both read ONE table per method
  (:mod:`.realization`), so the query and the realization cannot disagree.
* ``Approximated.of(name, issue)`` is the one-member spelling the tables write; ``names`` and ``issues`` are derived
  projections of the stored pairs, never stored.
* The sum is sealed at the SUBCLASS level and by identity: ``Realized``, ``Approximated`` and ``Refused`` are the only
  members, each a DIRECT subclass of :class:`Admission`; once the module has defined them every further subclass is
  refused, a subclass of a member and a foreign class reusing a member's name alike. Each member parses its fields (a
  ``Realized`` holds a ``bool``, an ``Approximated`` a non-empty set of ``(str, int)`` pairs, a ``Refused`` an axis, a
  reason and an optional issue).
* :attr:`FaceSite.hollow_inner` is DERIVED, never stored: the face name and the coordinate system already say whether
  a face is the inner surface of a hollow radial body, because a solid body has no inner face (``face_inventory``
  omits it). A stored flag would be a second spelling of that fact, constructible in contradiction with it (the
  orchestrator's ruling, 2026-10-23, after the S2a review).

Theory: :ref:`bc-architecture-of-record`.

Layer: ``transport`` (L2). Imports ``geometry`` (``CoordSystem`` and the face vocabulary) only.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum

from orpheus.geometry.coord import CoordSystem
from orpheus.geometry.faces import face_normal

__all__ = ["Admission", "Approximated", "FaceSite", "Realized", "RefusalAxis", "Refused", "SiteClass", "hollow_inner_refusal"]


class RefusalAxis(Enum):
    """On which axis a method refuses a declaration; ``PLUMBING`` is the fourth thing, not an axis of theory."""

    ANGULAR_RESOLUTION = "angular resolution"  # can the method represent R's angular structure?
    TOPOLOGY = "topology"  # can its operator couple two faces (spatial / topological)?
    STATE_CONE = "state cone"  # is the value in the method's state cone (psi >= 0), or on its trace?
    PLUMBING = "plumbing"  # no theory refuses it: a missing registry entry or an unwired channel


class Admission:
    """The sum ``Realized | Approximated | Refused``, sealed, with the lattice's ``meet`` defined once."""

    __slots__ = ()

    def __init_subclass__(cls, **kwargs: object) -> None:
        super().__init_subclass__(**kwargs)
        # The sum is closed: the three members below are its only subclasses, each direct. The seal is by identity, not by
        # name: once the module has defined the three, ``_SUM_CLOSED`` refuses every further subclass, whatever it is called.
        if _SUM_CLOSED or cls.__bases__ != (Admission,):
            raise TypeError(f"{cls.__name__}: Admission is a closed sum of Realized, Approximated and Refused, each a direct subclass")

    @staticmethod
    def meet(*members: Admission) -> Admission:
        """The lesser of ``members`` in ``Realized(True) > Realized(False) > Approximated > Refused`` (see the module docstring)."""
        result: Admission = Realized(True)
        for member in members:
            if isinstance(result, Refused):
                break
            if isinstance(result, Approximated) and isinstance(member, Approximated):
                result = Approximated(result.approximations | member.approximations)
            elif _rank(member) < _rank(result):
                result = member
        return result


_SUM_CLOSED = False


@dataclass(frozen=True)
class Realized(Admission):
    """Exact: ``faithful`` is True at tier 1 and False at tier 2."""

    faithful: bool

    def __post_init__(self) -> None:
        if not isinstance(self.faithful, bool):
            raise TypeError(f"Realized: faithful is a bool, got {self.faithful!r}")


@dataclass(frozen=True)
class Approximated(Admission):
    """Tier 3, or a named domain approximation: the non-empty SET of ``(name, issue)`` pairs, each approximation with its tracking issue."""

    approximations: frozenset[tuple[str, int]]

    def __post_init__(self) -> None:
        pairs = frozenset(self.approximations)
        if not pairs:
            raise ValueError("Approximated: an approximation names at least one (name, issue) pair")
        for pair in pairs:
            if not (isinstance(pair, tuple) and len(pair) == 2 and isinstance(pair[0], str) and isinstance(pair[1], int) and not isinstance(pair[1], bool)):
                raise TypeError(f"Approximated: a member is a (name: str, issue: int) pair, got {pair!r}")
        object.__setattr__(self, "approximations", pairs)

    @classmethod
    def of(cls, name: str, issue: int) -> Approximated:
        """The one-member spelling the tables write."""
        return cls(frozenset({(name, issue)}))

    @property
    def names(self) -> frozenset[str]:
        """The approximations' names (derived)."""
        return frozenset(name for name, _ in self.approximations)

    @property
    def issues(self) -> frozenset[int]:
        """The tracking issues (derived)."""
        return frozenset(issue for _, issue in self.approximations)


@dataclass(frozen=True)
class Refused(Admission):
    """Not realized; ``PLUMBING`` owes the ``issue`` that tracks the missing machinery."""

    axis: RefusalAxis
    reason: str
    issue: int | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.axis, RefusalAxis):
            raise TypeError(f"Refused: axis is a RefusalAxis, got {self.axis!r}")
        if not isinstance(self.reason, str):
            raise TypeError(f"Refused: reason is a str, got {self.reason!r}")
        if self.issue is not None and (not isinstance(self.issue, int) or isinstance(self.issue, bool)):
            raise TypeError(f"Refused: issue is an int or None, got {self.issue!r}")
        if self.axis is RefusalAxis.PLUMBING and self.issue is None:
            raise ValueError(f"Refused(PLUMBING, {self.reason!r}): a refusal with no theory behind it names the issue that tracks the missing machinery")


_SUM_CLOSED = True


def _rank(member: Admission) -> int:
    """The position in the lattice: three for ``Realized(True)`` down to zero for ``Refused``."""
    match member:
        case Realized(faithful=True):
            return 3
        case Realized():
            return 2
        case Approximated():
            return 1
        case Refused():
            return 0
        case _:
            raise TypeError(f"not an Admission: {member!r}")


class SiteClass(Enum):
    """The three GEOMETRIC classes of a face of the admission table."""

    CART = "cartesian"
    RAD_OUTER = "radial outer"
    RAD_INNER = "radial inner"


@dataclass(frozen=True)
class FaceSite:
    """Where a declaration is placed; derived from the carrier, computable before a method space exists."""

    face: str
    coord: CoordSystem
    dimension: int

    @property
    def is_curved(self) -> bool:
        """The coordinate system is not Cartesian and the face is not a Cartesian-axis face of a radial mesh."""
        return self.coord is not CoordSystem.CARTESIAN and face_normal(self.face)[0] == 0

    @property
    def site_class(self) -> SiteClass:
        """The geometric class of the face, one of Cartesian, radial outer or radial inner.

        The inner face has a negative outward sign.
        """
        if not self.is_curved:
            return SiteClass.CART
        return SiteClass.RAD_INNER if face_normal(self.face)[1] < 0 else SiteClass.RAD_OUTER

    @property
    def hollow_inner(self) -> bool:
        """Whether the face is the inner surface of a hollow radial body (derived, since a solid body has no inner face at all)."""
        return self.site_class is SiteClass.RAD_INNER


def hollow_inner_refusal(reason: str) -> Refused:
    r"""The refusal of a declaration the inner face of a hollow radial body cannot honour, citing #511 (ONE definition).

    The radial axis (:class:`~orpheus.mesh.axis.RadialAxisMesh`) carries one trace slot, the outer surface, so a method assembled
    on the axes' trace inventory has nowhere to put an inner condition: it computes the reflective cavity (a unit specular return)
    whatever is declared. A table whose method reads that inventory (SN, diffusion) refuses every other inner declaration through
    this constructor, ``reason`` naming what the method computes there. It is a per-TABLE admission and not a site rule of the
    shared dispatch: the reason for refusing an inner face is the method's (CP refuses a hollow body's inner face on its own guard,
    #513, and MoC and MC refuse a hollow body outright, #514 and #513), so a method that gains an inner-surface trace lifts it alone.

    SCOPE-BOUNDARY[guard] machinery: an inner-surface trace on the radial axis (a trace slot, and for SN a corner action, on the inner face).
    ruling: wave2_boundary_carve.md section 4.2 (the SN and diffusion columns), revised at the S3 integration 2026-10-23.
    revisit: #511, when the radial axis gains an inner-surface trace.
    """
    return Refused(RefusalAxis.PLUMBING, reason, 511)
