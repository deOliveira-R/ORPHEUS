r"""The reference specification (#405 P1 step 8).

A specification is the QUESTION a reference answers, never an answer. It is
the key the reference cache stores under, so it is a content value
(:class:`~orpheus.numerics.content.ContentIdentity`), and the layer of the
posing filtration it is posed at is its TYPE:

* :class:`InfiniteMediumSpecification` is posed on the materials layer alone:
  one material and a question. The infinite medium is the point in position
  AND in direction, so its problem lives on energy alone; it has no geometry
  field, because admitting a geometry admits spatial dimension, and with it
  more than one material and a direction chart (the user's ruling of
  2026-10-02; the theory page ``docs/theory/foundations/infinite_medium.rst``
  derives why this is exact for every question it admits);
* :class:`GeometrySpecification` is posed on a finite geometry: the materials,
  a :class:`~orpheus.geometry.structured_geometry.StructuredGeometry` with its
  boundary laws, and a question.

:data:`Specification` is the closed set of the two. Each type answers for
itself what the shared admission asks (the problem's materials, the regions
a regionwise datum has, the coordinates a ``Symbolic`` cannot depend on, how
an extent key resolves), so no code path asks which layer it is on.

Both are admitted in a CANONICAL form at construction:

* **the materials** of a geometry specification are restricted to those the
  geometry assigns (:meth:`~orpheus.data.materials.Materials.restrict`): a
  spectator changes no answer, so it is not part of the key (the leak
  principle); the group count is read over the materials kept;
* **the question's keys** are resolved: an ``Eigen`` parameter and every point
  key must be a coordinate of this problem, a
  :class:`~orpheus.data.cells.CellCoefficient` over its materials or a
  :class:`~orpheus.geometry.extent.GeometryExtent` of its geometry, and each
  is replaced by its resolved form (explicit non-zero cells), so two
  spellings of one direction are one specification. ``spec.question`` is
  therefore the canonical question, which may be unequal to the value the
  caller passed. It follows that RE-POSING a specification is done from the
  caller's question, never from ``spec.question``: ``replace(spec,
  materials=other)`` keeps the cells resolved against the old materials, so a
  question written ``every(fission emission)`` would then name only the
  materials that carried fission before;
* **the datum** of a ``FixedSource`` or a ``Response`` (a mesh-free function)
  must fit: its group count is the materials', a ``RegionwiseConstant`` has
  one value per region, and a ``Symbolic`` depends only on coordinates the
  problem has (none on the infinite medium; not the azimuth on a sphere,
  which has no azimuth reference).

A value the encoder cannot key (a boundary law holding a function, the MMS
inflow) is refused at construction: a specification exists to be keyed.

What is not here: the coordinates' charts (zero, physical value, admissible
range) and the mode law that finds a pole, #529; the nuclide-density
coordinate, deferred until the materials keep number densities; the
direction-dependent spatial marginal of a body problem, which is the
retraction of a geometry problem's solution along space, not a specification.
"""

from __future__ import annotations

import dataclasses
from dataclasses import dataclass
from functools import cached_property
from typing import TYPE_CHECKING, Any, ClassVar, TypeAlias, assert_never, get_args

from orpheus.data.cells import CellCoefficient
from orpheus.data.macro_xs.mixture import Mixture
from orpheus.data.materials import Materials
from orpheus.geometry.extent import GeometryExtent
from orpheus.geometry.structured_geometry import StructuredGeometry
from orpheus.numerics.content import ContentIdentity, FrozenMapping, content_digest
from orpheus.numerics.mesh_free_function import MeshFreeFunction, RegionwiseConstant, Symbolic
from orpheus.numerics.observable import Eigenvalue, FluxIntegral, Observable, PointValue, Ratio
from orpheus.numerics.question import Eigen, FixedSource, Question, Response
from orpheus.numerics.scalars import parse_integer

if TYPE_CHECKING:
    import sympy

__all__ = ["Coordinate", "GeometrySpecification", "InfiniteMediumSpecification", "Specification", "admit_observable"]

Coordinate: TypeAlias = CellCoefficient | GeometryExtent
"""The coordinates a specification resolves a question's keys to."""


def _names(alias: Any) -> str:
    return " or ".join(kind.__name__ for kind in get_args(alias))


@dataclass(frozen=True, eq=False)
class InfiniteMediumSpecification(ContentIdentity):
    """One material and a question, posed on energy alone (no geometry: the infinite medium)."""

    material_id: int
    mixture: Mixture
    question: Question

    #: The infinite medium is one region.
    n_regions: ClassVar[int] = 1

    def __post_init__(self) -> None:
        object.__setattr__(self, "material_id", parse_integer(self.material_id, "InfiniteMediumSpecification", "the material id"))
        if not isinstance(self.mixture, Mixture):
            raise TypeError(f"InfiniteMediumSpecification: the mixture is a Mixture, got a {type(self.mixture).__name__}")
        _admit(self)

    @property
    def materials(self) -> Materials:
        """The problem's materials: the one material, under its id."""
        return Materials({self.material_id: self.mixture})

    @property
    def n_groups(self) -> int:
        """The one material's group count."""
        return self.mixture.ng

    @property
    def unreadable(self) -> tuple["sympy.Symbol", ...]:
        """Every coordinate: the infinite medium has no position and no direction chart."""
        return (Symbolic.r, Symbolic.mu, Symbolic.phi)

    def resolve_extent(self, extent: GeometryExtent) -> GeometryExtent:
        raise ValueError(f"{type(self).__name__}: {extent!r} names an extent, and the infinite medium has no geometry")


@dataclass(frozen=True, eq=False)
class GeometrySpecification(ContentIdentity):
    """The materials, a finite geometry and a question; the materials kept are those the geometry assigns."""

    materials: Materials
    geometry: StructuredGeometry
    question: Question

    def __post_init__(self) -> None:
        if not isinstance(self.materials, Materials):
            raise TypeError(f"GeometrySpecification: the materials are a Materials, got a {type(self.materials).__name__}")
        if not isinstance(self.geometry, StructuredGeometry):
            raise TypeError(f"GeometrySpecification: the geometry is a StructuredGeometry, got a {type(self.geometry).__name__}")
        object.__setattr__(self, "materials", self.materials.restrict(set(self.geometry.mat_ids)))
        self.n_groups  # the group-count rule, read once over the materials kept, before the keys resolve
        _admit(self)

    @cached_property
    def n_groups(self) -> int:
        """The group count of the materials the geometry assigns (refused if they disagree)."""
        return self.materials.uniform_group_count()

    @property
    def n_regions(self) -> int:
        """One region per interval of the geometry."""
        return len(self.geometry.intervals)

    @property
    def unreadable(self) -> tuple["sympy.Symbol", ...]:
        """The azimuth on a chart with no azimuth reference (the sphere); otherwise none."""
        return (Symbolic.phi,) if self.geometry.coord.angular_chart.azimuth_reference is None else ()

    def resolve_extent(self, extent: GeometryExtent) -> GeometryExtent:
        return extent.resolve(self.geometry)


Specification: TypeAlias = InfiniteMediumSpecification | GeometrySpecification
"""A reference specification: the layer it is posed at is its type."""


def _resolve(key: Any, where: str, spec: Specification) -> Coordinate:
    """The canonical form of one key: a coordinate of this problem, resolved."""
    match key:
        case CellCoefficient():
            return key.resolve(spec.materials)
        case GeometryExtent():
            return spec.resolve_extent(key)
        case _:
            raise TypeError(f"{where}: {key!r} (a {type(key).__name__}) is not a coordinate; a key is a {_names(Coordinate)}")


def _admit_datum(datum: MeshFreeFunction, role: str, spec: Specification) -> None:
    """The mesh-free datum fits the materials' groups, the problem's regions and its coordinates."""
    if datum.n_groups != spec.n_groups:
        raise ValueError(f"{type(spec).__name__}: the {role} has {datum.n_groups} groups; the materials have {spec.n_groups}")
    match datum:
        case RegionwiseConstant():
            if datum.n_regions != spec.n_regions:
                raise ValueError(f"{type(spec).__name__}: the {role} has {datum.n_regions} regions; the problem has {spec.n_regions}")
        case Symbolic():
            dependent = tuple(c for c in spec.unreadable if datum.depends_on(c))
            if dependent:
                names = ", ".join(c.name for c in dependent)
                raise ValueError(f"{type(spec).__name__}: the {role} depends on {names}, which this problem has no coordinate to read")
        case _:
            assert_never(datum)


def admit_observable(observable: Observable, specification: Specification) -> None:
    """Refuse an observable this specification cannot pose; the one admission every reader reuses.

    A flux integral's weight must fit the problem as a question's datum must
    (its groups, its regions, the coordinates a ``Symbolic`` weight may depend
    on); a ratio's two operands must each fit; a point value's group must be
    one of the problem's and its position must lie on the geometry, so the
    infinite medium, which has no position, refuses it; an eigenvalue exists
    only for an eigen question. A reader admits an observable once, at the
    point it accepts it (a
    :class:`~orpheus.reference.published.PublishedSolution` at construction,
    for every observable it prints; a reference solution's ``read`` at #405
    P2 step 6), so the refusal is decided once (#405 P2, the elegance review
    of step 3).
    """
    where = type(specification).__name__
    match observable:
        case FluxIntegral():
            _admit_datum(observable.weight, "weight", specification)
        case Ratio():
            admit_observable(observable.numerator, specification)
            admit_observable(observable.denominator, specification)
        case PointValue():
            if observable.group >= specification.n_groups:
                raise ValueError(
                    f"{where}: the point value's group {observable.group} is not among the problem's "
                    f"{specification.n_groups} groups"
                )
            match specification:
                case InfiniteMediumSpecification():
                    raise ValueError(f"{where}: a point value names a position, and the infinite medium has none")
                case GeometrySpecification():
                    low, high = specification.geometry.breakpoints[0], specification.geometry.breakpoints[-1]
                    if not low <= observable.position <= high:
                        raise ValueError(
                            f"{where}: the point value's position {observable.position!r} is off the geometry "
                            f"[{low!r}, {high!r}]"
                        )
                case _:
                    assert_never(specification)
        case Eigenvalue():
            if not isinstance(specification.question, Eigen):
                raise ValueError(
                    f"{where}: an eigenvalue is read off an eigen question's answer, and this question is a "
                    f"{type(specification.question).__name__}"
                )
        case _:
            assert_never(observable)


def _canonical_point(point: Any, spec: Specification) -> FrozenMapping[Coordinate, float]:
    """The point with every key resolved; two keys naming one coordinate are refused."""
    resolved: dict[Coordinate, tuple[Any, float]] = {}
    for key, offset in point.items():
        coordinate = _resolve(key, f"{type(spec).__name__}: the point key", spec)
        if coordinate in resolved:
            raise ValueError(f"{type(spec).__name__}: the point keys {resolved[coordinate][0]!r} and {key!r} name the same coordinate")
        resolved[coordinate] = (key, offset)
    return FrozenMapping((coordinate, offset) for coordinate, (_, offset) in resolved.items())


def _canonical_question(question: Any, spec: Specification) -> Question:
    """The question with its datum admitted and every key resolved."""
    match question:
        case Eigen():
            parameter = _resolve(question.parameter, f"{type(spec).__name__}: the parameter", spec)
            return dataclasses.replace(question, parameter=parameter, point=_canonical_point(question.point, spec))
        case FixedSource():
            _admit_datum(question.source, "source", spec)
            return dataclasses.replace(question, point=_canonical_point(question.point, spec))
        case Response():
            _admit_datum(question.detector, "detector", spec)
            return dataclasses.replace(question, point=_canonical_point(question.point, spec))
        case _:
            raise TypeError(f"{type(spec).__name__}: the question is an {_names(Question)}, got a {type(question).__name__}")


def _admit(spec: Specification) -> None:
    """The canonical question, then the digest: a specification exists to be keyed."""
    object.__setattr__(spec, "question", _canonical_question(spec.question, spec))
    content_digest(spec)
