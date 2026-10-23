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
is the only spelling of the meet: the lesser member; the meet of two :class:`Approximated` is
``Approximated(x.names | y.names, x.issues | y.issues)`` (a SET of names, because ``meet`` on a total order would drop
one); the meet of two :class:`Refused` is the first argument's (documented, deterministic); the meet of NO members is
``Realized(True)`` (an empty network and no source need no realizer). A wall's admission is the meet over its nonzero
chains and its source.

**Why a class and not a ``TypeAlias``.** The specification writes ``Admission: TypeAlias = Realized | Approximated |
Refused`` AND ``Admission.meet``; a ``typing.Union`` carries no attributes, so the two statements cannot both hold. The
sum is therefore a sealed base class (:class:`Admission`) with exactly three subclasses: ``isinstance``, ``match``
and annotations read as for the union, and ``meet`` lives on the type (reported in NEEDS).

**Where a declaration is placed.** :class:`FaceSite` is derived from the carrier, computable before a method space
exists: the ``face``, the ``coord`` system, the ``dimension`` and whether the face is the hollow inner face of a radial
mesh. Its :attr:`~FaceSite.site_class` is one of the three GEOMETRIC classes of the admission table
(:class:`SiteClass`: Cartesian, radial outer, radial inner), derived once; the pin and the realizers read this name.
CP's slab-left-face rule is CP's own predicate on ``site.face`` (decision 56).

**Gotchas.**

* There are two levels: the TYPE level (the front end's drop-down; the upper bound) and the INSTANCE level (per face
  and per value, plus the quadrature's equivariance certified by ``ordinate_pairing``). Both read ONE table per method
  (:mod:`.realization`), so the query and the realization cannot disagree.
* ``Approximated.of(name, issue)`` is the one-member spelling the tables write.

Layer: ``transport`` (L2). Imports ``geometry`` (``CoordSystem`` and the face vocabulary) only.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum

from orpheus.geometry.coord import CoordSystem
from orpheus.geometry.faces import face_normal

__all__ = ["Admission", "Approximated", "FaceSite", "Realized", "RefusalAxis", "Refused", "SiteClass"]


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
        if len(Admission.__subclasses__()) > 3:  # the sum is closed: Realized, Approximated, Refused
            raise TypeError(f"{cls.__name__}: Admission is a closed sum of Realized, Approximated and Refused")

    @staticmethod
    def meet(*members: Admission) -> Admission:
        """The lesser of ``members`` in ``Realized(True) > Realized(False) > Approximated > Refused`` (see the module docstring)."""
        result: Admission = Realized(True)
        for member in members:
            if isinstance(result, Refused):
                break
            if isinstance(result, Approximated) and isinstance(member, Approximated):
                result = Approximated(result.names | member.names, result.issues | member.issues)
            elif _rank(member) < _rank(result):
                result = member
        return result


@dataclass(frozen=True)
class Realized(Admission):
    """Exact: ``faithful`` is True at tier 1 and False at tier 2."""

    faithful: bool


@dataclass(frozen=True)
class Approximated(Admission):
    """Tier 3, or a named domain approximation: the SET of approximations' ``names`` and tracking ``issues``."""

    names: frozenset[str]
    issues: frozenset[int]

    @classmethod
    def of(cls, name: str, issue: int) -> Approximated:
        """The one-member spelling the tables write."""
        return cls(frozenset({name}), frozenset({issue}))


@dataclass(frozen=True)
class Refused(Admission):
    """Not realized; ``PLUMBING`` owes the ``issue`` that tracks the missing machinery."""

    axis: RefusalAxis
    reason: str
    issue: int | None = None

    def __post_init__(self) -> None:
        if self.axis is RefusalAxis.PLUMBING and self.issue is None:
            raise ValueError(f"Refused(PLUMBING, {self.reason!r}): a refusal with no theory behind it names the issue that tracks the missing machinery")


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
    hollow_inner: bool

    @property
    def is_curved(self) -> bool:
        """The coordinate system is not Cartesian and the face is not a Cartesian-axis face of a radial mesh."""
        return self.coord is not CoordSystem.CARTESIAN and face_normal(self.face)[0] == 0

    @property
    def site_class(self) -> SiteClass:
        """The geometric class: Cartesian, radial outer or radial inner (the inner face has a negative outward sign)."""
        if not self.is_curved:
            return SiteClass.CART
        return SiteClass.RAD_INNER if face_normal(self.face)[1] < 0 else SiteClass.RAD_OUTER
