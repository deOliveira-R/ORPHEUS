r"""The reference specification (#405 P1 step 8).

A :class:`Specification` is the QUESTION a reference answers, never an
answer: the materials, the geometry with its boundary laws (``None`` for the
infinite medium), and one question value of :mod:`orpheus.numerics.question`.
It is the key the reference cache stores under, so it is a content value
(:class:`~orpheus.numerics.content.ContentIdentity`) and is admitted in a
CANONICAL form at construction:

* **the materials** are restricted to those the geometry assigns
  (:meth:`~orpheus.data.materials.Materials.restrict`): a spectator material
  changes no answer, so it is not part of the key (the leak principle);
  with no geometry there is exactly one material, the infinite medium;
* **the question's keys** are resolved: an ``Eigen`` parameter and every point
  key must be a coordinate of this specification, a
  :class:`~orpheus.data.cells.CellCoefficient` over its materials or a
  :class:`~orpheus.geometry.extent.GeometryExtent` of its geometry, and each
  is replaced by its resolved form (explicit non-zero cells), so two
  spellings of one direction are one specification. ``spec.question`` is
  therefore the canonical question, which may be unequal to the value the
  caller passed;
* **the datum** of a ``FixedSource`` or a ``Response`` (a mesh-free function)
  must fit: its group count is the materials', a ``RegionwiseConstant`` has
  one value per interval (one for the infinite medium), and a ``Symbolic`` is
  read in a chart this geometry has (no dependence on position or direction
  with no geometry; none on the azimuth on a sphere, which has no azimuth
  reference).

A value the encoder cannot key (a boundary law holding a function, the MMS
inflow) is refused at construction: a specification exists to be keyed.

What is not here: the coordinates' charts (zero, physical value, admissible
range) and the mode law that finds a pole, #529; the nuclide-density
coordinate, deferred until the materials keep number densities.
"""

from __future__ import annotations

import dataclasses
from dataclasses import dataclass
from typing import Any

from orpheus.data.cells import CellCoefficient
from orpheus.data.materials import Materials
from orpheus.geometry.coord import CoordSystem
from orpheus.geometry.extent import GeometryExtent
from orpheus.geometry.structured_geometry import StructuredGeometry
from orpheus.numerics.content import ContentIdentity, FrozenMapping, content_digest
from orpheus.numerics.mesh_free_function import MeshFreeFunction, RegionwiseConstant, Symbolic
from orpheus.numerics.question import Eigen, FixedSource, Question, Response

__all__ = ["Specification"]


def _resolve(key: Any, where: str, materials: Materials, geometry: StructuredGeometry | None) -> Any:
    """The canonical form of one key: a coordinate of this specification, resolved."""
    match key:
        case CellCoefficient():
            return key.resolve(materials)
        case GeometryExtent():
            if geometry is None:
                raise ValueError(
                    f"{where}: a GeometryExtent names an extent, and the infinite medium (no geometry) has none"
                )
            return key.resolve(geometry)
        case _:
            raise TypeError(
                f"{where}: {key!r} (a {type(key).__name__}) is not a coordinate; "
                f"a key is a CellCoefficient or a GeometryExtent"
            )


def _admit_datum(datum: MeshFreeFunction, role: str, n_groups: int, geometry: StructuredGeometry | None) -> None:
    """The mesh-free datum fits the materials' groups and the geometry's chart."""
    if datum.n_groups != n_groups:
        raise ValueError(f"Specification: the {role} has {datum.n_groups} groups; the materials have {n_groups}")
    if isinstance(datum, RegionwiseConstant):
        n_regions = 1 if geometry is None else len(geometry.intervals)
        if datum.n_regions != n_regions:
            where = "the infinite medium has 1 region" if geometry is None else f"the geometry has {n_regions} intervals"
            raise ValueError(f"Specification: the {role} has {datum.n_regions} regions; {where}")
        return
    if geometry is None:
        if datum.depends_on(Symbolic.r, Symbolic.mu, Symbolic.phi):
            raise ValueError(
                f"Specification: the {role} depends on position or direction, and with no geometry "
                f"there is no coordinate system to read them in"
            )
    elif geometry.coord is CoordSystem.SPHERICAL and datum.depends_on(Symbolic.phi):
        raise ValueError(f"Specification: the {role} depends on the azimuth phi, which has no reference on a sphere")


def _canonical_question(question: Question, materials: Materials, geometry: StructuredGeometry | None) -> Question:
    """The question with every key resolved; two point keys naming one coordinate are refused."""
    changes: dict[str, Any] = {}
    if isinstance(question, Eigen):
        changes["parameter"] = _resolve(question.parameter, "Specification: the parameter", materials, geometry)
    resolved: dict[Any, tuple[Any, float]] = {}
    for key, offset in question.point.items():
        coordinate = _resolve(key, "Specification: the point key", materials, geometry)
        if coordinate in resolved:
            raise ValueError(
                f"Specification: the point keys {resolved[coordinate][0]!r} and {key!r} name the same coordinate"
            )
        resolved[coordinate] = (key, offset)
    changes["point"] = FrozenMapping((coordinate, offset) for coordinate, (_, offset) in resolved.items())
    return dataclasses.replace(question, **changes)


@dataclass(frozen=True, eq=False)
class Specification(ContentIdentity):
    """Materials, a geometry (``None``: the infinite medium) and a question, in canonical form."""

    materials: Materials
    geometry: StructuredGeometry | None
    question: Question

    def __post_init__(self) -> None:
        if not isinstance(self.materials, Materials):
            raise TypeError(f"Specification: the materials are a Materials, got a {type(self.materials).__name__}")
        geometry = self.geometry
        if geometry is None:
            if len(self.materials) != 1:
                raise ValueError(
                    f"Specification: with no geometry, the infinite medium has one material; "
                    f"{len(self.materials)} are declared ({sorted(self.materials.ids)})"
                )
            materials = self.materials
        elif isinstance(geometry, StructuredGeometry):
            materials = self.materials.restrict(sorted(set(geometry.mat_ids)))
        else:
            raise TypeError(f"Specification: the geometry is a StructuredGeometry or None, got a {type(geometry).__name__}")
        n_groups = materials.uniform_group_count()
        question = self.question
        if not isinstance(question, (Eigen, FixedSource, Response)):
            raise TypeError(
                f"Specification: the question is an Eigen, a FixedSource or a Response, got a {type(question).__name__}"
            )
        if isinstance(question, FixedSource):
            _admit_datum(question.source, "source", n_groups, geometry)
        elif isinstance(question, Response):
            _admit_datum(question.detector, "detector", n_groups, geometry)
        object.__setattr__(self, "materials", materials)
        object.__setattr__(self, "question", _canonical_question(question, materials, geometry))
        content_digest(self)  # admission: a specification exists to be keyed
