r"""Structured 1-D geometry — pure geometry layer, no mesh concerns.

A :class:`StructuredGeometry` is an interval of positions in one
coordinate system, cut into material intervals, with a boundary law at
every point of its boundary:

* **the coordinate system** (:class:`CoordSystem`: slab, cylinder,
  sphere), which fixes the measure of an interval and the topology of
  its boundary;
* **the breakpoints** :math:`r_0 < r_1 < \dots < r_R`, the positions
  where one material interval ends and the next begins; they are stored
  as given and never re-derived;
* **one material id per interval** :math:`[r_k, r_{k+1}]`;
* **one boundary law per boundary point**, in the order (inner, outer).

It is "structured" in the memory-layout sense: the intervals are
traversed in order with no lookup. This is the counterpart to a future
``ConstructiveSolidGeometry`` (CSG) that would describe unstructured,
Boolean-composed shapes, which need different machinery (BVH lookups,
ray–surface intersection).

The class lives at the **geometry layer**, not the mesh layer, and
knows nothing about cell counts, discretisation rules or numerics. The
same geometry is discretised in several ways for different studies
(mesh refinement, equal-width against equal-volume cells) at the mesh
layer, :mod:`orpheus.mesh`.

The boundary is derived, not declared
-------------------------------------

The boundary of the region :math:`[r_0, r_R]` is its topological
boundary in its coordinate system, so the number of laws a geometry
takes follows from the coordinate system and :math:`r_0`:

======================  =====================  =====================
coordinate system       :math:`r_0 = 0`        :math:`r_0 > 0`
======================  =====================  =====================
Cartesian (slab)        2 (left, right)        2 (left, right)
cylindrical, spherical  1 (outer)              2 (inner, outer)
======================  =====================  =====================

A slab admits any :math:`r_0` (a position on a line). On a cylinder or
a sphere :math:`r_0 \ge 0` is a radius, and the centre :math:`r = 0` of
a solid body is an interior point of the region, not a boundary point:
the regularity of the solution there is a property of the coordinate
chart, and no law is declared at it. A hollow body
(:math:`r_0 > 0`) has an inner surface, which carries a law like any
other boundary point. :attr:`StructuredGeometry.boundary_points` lists
the positions, paired one to one with
:attr:`StructuredGeometry.boundaries`.

Architectural role
------------------

* **Reference solution generators** (``Billiard``, ``MomentSpace``,
  ``Spectrum``, ``BasisSpace``) consume a :class:`StructuredGeometry`
  directly: no mesh, no cell counts.
* **Discrete production solvers** (``solve_cp``, ``solve_sn``,
  ``solve_moc``, ``solve_mc``) consume a :class:`~orpheus.mesh.Mesh1D`
  built from the geometry; the discretisation is supplied at that build
  step, never stored on the geometry.

Examples
--------

A bare-critical sphere, one region, one outer law::

    geom = StructuredGeometry(
        coord=CoordSystem.SPHERICAL,
        breakpoints=(0.0, 2.872),
        mat_ids=(0,),
        boundaries=(BC.vacuum,),
    )

A reflected slab (reflector | core | reflector), two laws::

    geom = StructuredGeometry.from_thicknesses(
        coord=CoordSystem.CARTESIAN,
        thicknesses=(0.5, 2.0, 0.5),
        mat_ids=(1, 0, 1),
        boundaries=(BC.vacuum, BC.vacuum),
    )

A hollow sphere, an inner and an outer law::

    geom = StructuredGeometry(
        coord=CoordSystem.SPHERICAL,
        breakpoints=(0.5, 1.0, 2.0),
        mat_ids=(1, 0),
        boundaries=(BC.reflective, BC.vacuum),
    )

A PWR pin cell through the Wigner–Seitz factory::

    geom = StructuredGeometry.wigner_seitz_pin_cell(
        r_fuel=0.9, r_clad=1.1, pitch=3.6,
    )
"""
from __future__ import annotations

import itertools
import math
from collections.abc import Iterable
from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

from .coord import CoordSystem
from .scalars import parse_entries, parse_integer, parse_positive_real, parse_real
from .boundary import BC

if TYPE_CHECKING:
    from .boundary import BoundaryTraceLaw


@dataclass(frozen=True, kw_only=True)
class StructuredGeometry:
    r"""A 1-D geometry: a coordinate system, breakpoints, one material id
    per interval and one boundary law per boundary point.

    Parameters
    ----------
    coord : CoordSystem
        The coordinate system. The retired string kind tags (``"SLB"``,
        ``"CYL"``, ``"SPH"``) are refused.
    breakpoints : sequence of real numbers
        :math:`r_0 < r_1 < \dots < r_R`, at least two, finite, strictly
        increasing; :math:`r_0 \ge 0` on a cylinder or a sphere. Stored
        as a tuple of ``float``, bit for bit as given: a breakpoint is
        never re-derived from widths (the thickness constructor is
        :meth:`from_thicknesses`).
    mat_ids : sequence of int
        One material id per interval :math:`[r_k, r_{k+1}]`, inside-out
        on a cylinder or a sphere and left to right on a slab. Adjacent
        intervals may share a material. The ids key the
        ``materials: dict[int, Mixture]`` that consumers receive beside
        the geometry.
    boundaries : tuple of BC or BoundaryTraceLaw
        One law per boundary point, paired one to one with
        :attr:`boundary_points` (inner first). A declaration is either a
        :class:`BC` tag or an already-typed ``BoundaryTraceLaw``; the
        law arm exists because a tag cannot carry a function (a
        prescribed inflow whose source is a manufactured solution), and
        declaring such a law on the geometry is what makes it survive
        the method-mesh rebuild every public solver entry performs.
        ``None`` is not a law.

    Notes
    -----
    **What is not here.** No cell counts and no discretisation rule (the
    mesh layer's); no critical dimension (a registry's truth record); no
    group count (the materials'); no infinite medium (a problem with no
    geometry, or a finite domain with reflective laws).

    **The geometry never interprets a law.** Each method's mesh resolves
    what a tag means through its own admission table.
    """

    coord: CoordSystem
    breakpoints: tuple[float, ...]
    mat_ids: tuple[int, ...]
    boundaries: "tuple[BC | BoundaryTraceLaw, ...]"

    def __post_init__(self) -> None:
        _parse_coord(self.coord)
        object.__setattr__(
            self, "breakpoints", _parse_breakpoints(self.coord, self.breakpoints),
        )
        object.__setattr__(
            self, "mat_ids", _parse_mat_ids(self.mat_ids, len(self.breakpoints) - 1),
        )
        object.__setattr__(
            self, "boundaries",
            parse_entries(self.boundaries, "StructuredGeometry.boundaries", "boundary laws"),
        )
        _check_boundaries(self)

    @property
    def is_hollow(self) -> bool:
        r"""A cylinder or sphere whose first breakpoint :math:`r_0 > 0`.

        The one place the centre's role is decided: on a solid cylinder
        or sphere (:math:`r_0 = 0`) the centre is an interior point of the
        region and carries no law; a hollow one has an inner surface. A
        slab has no centre and is never hollow.
        """
        return _is_hollow(self.coord, self.breakpoints)

    @property
    def boundary_points(self) -> tuple[float, ...]:
        r"""The positions of the boundary, inner first; one law each.

        :math:`(r_0, r_R)` on a slab and on a hollow cylinder or sphere;
        :math:`(r_R,)` on a solid cylinder or sphere, whose centre is an
        interior point.
        """
        return _boundary_points(self.coord, self.breakpoints)

    @property
    def intervals(self) -> tuple[tuple[float, float], ...]:
        r"""The material intervals :math:`[r_k, r_{k+1}]`, in order."""
        return tuple(itertools.pairwise(self.breakpoints))

    def measure(self, edges: np.ndarray) -> np.ndarray:
        r"""The measures of the cells between ``edges``, in this geometry's coordinate system.

        Lengths on a slab, areas per unit height on a cylinder, volumes on a
        sphere: :meth:`CoordSystem.measure`, the one definition. An
        interval's measure is the one-cell case, ``measure([r_k, r_{k+1}])``.
        """
        return self.coord.measure(edges)

    @property
    def domain_extent_cm(self) -> float:
        r"""The width :math:`r_R - r_0` of the interval of positions, in cm.

        The full slab width on a slab (the production convention: a user
        with "a slab of thickness L" means L); the outer radius on a
        solid cylinder or sphere; the shell thickness on a hollow one.
        """
        return self.breakpoints[-1] - self.breakpoints[0]

    # ─────────────────────────────────────────────────────────────────
    # Constructors that name each boundary point's law
    # ─────────────────────────────────────────────────────────────────

    @classmethod
    def slab(
        cls,
        breakpoints: Iterable[float],
        mat_ids: Iterable[int],
        *,
        left: "BC | BoundaryTraceLaw",
        right: "BC | BoundaryTraceLaw",
    ) -> "StructuredGeometry":
        r"""A slab, with the law at its left and at its right face."""
        return cls(
            coord=CoordSystem.CARTESIAN, breakpoints=tuple(breakpoints),
            mat_ids=tuple(mat_ids), boundaries=(left, right),
        )

    @classmethod
    def cylinder(
        cls,
        breakpoints: Iterable[float],
        mat_ids: Iterable[int],
        *,
        outer: "BC | BoundaryTraceLaw",
        inner: "BC | BoundaryTraceLaw | None" = None,
    ) -> "StructuredGeometry":
        r"""A cylinder, with the law at its outer surface, and at its inner one when hollow.

        ``inner`` is given exactly when the first breakpoint is positive: the
        centre of a solid cylinder is an interior point and carries no law.
        """
        return cls._radial(CoordSystem.CYLINDRICAL, breakpoints, mat_ids, inner, outer)

    @classmethod
    def sphere(
        cls,
        breakpoints: Iterable[float],
        mat_ids: Iterable[int],
        *,
        outer: "BC | BoundaryTraceLaw",
        inner: "BC | BoundaryTraceLaw | None" = None,
    ) -> "StructuredGeometry":
        r"""A sphere, with the law at its outer surface, and at its inner one when hollow.

        ``inner`` is given exactly when the first breakpoint is positive: the
        centre of a solid sphere is an interior point and carries no law.
        """
        return cls._radial(CoordSystem.SPHERICAL, breakpoints, mat_ids, inner, outer)

    @classmethod
    def _radial(
        cls,
        coord: CoordSystem,
        breakpoints: Iterable[float],
        mat_ids: Iterable[int],
        inner: "BC | BoundaryTraceLaw | None",
        outer: "BC | BoundaryTraceLaw",
    ) -> "StructuredGeometry":
        # An absent ``inner`` keyword is no law at all, never a stored None;
        # whether the body needs one is decided by the one boundary check.
        laws = (outer,) if inner is None else (inner, outer)
        return cls(
            coord=coord, breakpoints=tuple(breakpoints),
            mat_ids=tuple(mat_ids), boundaries=laws,
        )

    @classmethod
    def uniform_boundary(
        cls,
        coord: CoordSystem,
        breakpoints: Iterable[float],
        mat_ids: Iterable[int],
        law: "BC | BoundaryTraceLaw",
    ) -> "StructuredGeometry":
        r"""The geometry whose every boundary point carries ``law``.

        Coordinate-generic: two laws on a slab or a hollow body, one on a
        solid cylinder or sphere. For a test parametrised over coordinate
        systems whose body has one law on its whole boundary.
        """
        parsed = _parse_breakpoints(_parse_coord(coord), breakpoints)
        return cls(
            coord=coord, breakpoints=parsed, mat_ids=tuple(mat_ids),
            boundaries=(law,) * len(_boundary_points(coord, parsed)),
        )

    @classmethod
    def from_homogeneous(
        cls, width: float, boundary: "BC | BoundaryTraceLaw",
    ) -> "StructuredGeometry":
        r"""A slab ``[0, width]`` of material 0, with ``boundary`` on both faces.

        The infinite medium as a test realises it: one finite slab whose
        faces carry one law (reflective, for the infinite medium itself).
        Only a slab: a finite curvilinear body is not an infinite medium in
        continuous transport.
        """
        extent = parse_positive_real(
            width, "StructuredGeometry.from_homogeneous", "the width",
        )
        return cls.slab((0.0, extent), (0,), left=boundary, right=boundary)

    # ─────────────────────────────────────────────────────────────────
    # Constructors for the callers that speak another vocabulary
    # ─────────────────────────────────────────────────────────────────

    @classmethod
    def from_thicknesses(
        cls,
        *,
        coord: CoordSystem,
        thicknesses: Iterable[float],
        mat_ids: Iterable[int],
        boundaries: "tuple[BC | BoundaryTraceLaw, ...]",
        r_0: float = 0.0,
    ) -> "StructuredGeometry":
        r"""The geometry whose intervals have the given widths, from ``r_0``.

        For the registries that publish a layered configuration as
        thicknesses. The breakpoints are the left fold
        :math:`r_{k+1} = r_k + t_k` evaluated in order, the same
        sequential sum every mesh built from thicknesses has used, so a
        registry's geometry keeps its bits. A thickness :math:`\le 0` is
        refused as a non-increasing breakpoint pair.
        """
        breakpoints = tuple(itertools.accumulate(
            (
                parse_real(t, f"StructuredGeometry.from_thicknesses: thicknesses[{k}]")
                for k, t in enumerate(thicknesses)
            ),
            initial=parse_real(r_0, "StructuredGeometry.from_thicknesses: r_0"),
        ))
        return cls(
            coord=coord,
            breakpoints=breakpoints,
            mat_ids=tuple(mat_ids),
            boundaries=boundaries,
        )

    @classmethod
    def wigner_seitz_pin_cell(
        cls,
        *,
        r_fuel: float = 0.9,
        r_clad: float = 1.1,
        pitch: float = 3.6,
    ) -> "StructuredGeometry":
        r"""Wigner–Seitz equivalent pin-cell geometry.

        Replaces a square unit cell of side ``pitch`` with a cylinder of
        equal cross-sectional area,

        .. math::

            r_{\rm cell} = \frac{\rm pitch}{\sqrt{\pi}},

        cut at the fuel and cladding radii: breakpoints
        :math:`(0, r_{\rm fuel}, r_{\rm clad}, r_{\rm cell})`, material
        ids fuel ``2``, cladding ``1``, coolant ``0``. The outer law is
        ``BC("white")``, and it is part of the model, not a default: the
        Wigner–Seitz cell IS the cylindricalised lattice cell with isotropic
        re-entry at its outer surface, the assumption that maps a periodic
        lattice to one cell. The same radii under another law are a
        different body, built with :meth:`cylinder`.

        Parameters
        ----------
        r_fuel, r_clad : float
            The fuel-pellet and cladding outer radii (cm).
        pitch : float
            The square unit cell's side (cm).
        """
        r_cell = float(pitch / np.sqrt(np.pi))
        return cls(
            coord=CoordSystem.CYLINDRICAL,
            breakpoints=(0.0, float(r_fuel), float(r_clad), r_cell),
            mat_ids=(2, 1, 0),
            boundaries=(BC("white"),),
        )

    @classmethod
    def pwr_slab_half_cell(
        cls,
        *,
        fuel_half: float = 0.9,
        clad_thick: float = 0.2,
        cool_thick: float = 0.7,
    ) -> "StructuredGeometry":
        r"""Cartesian 1-D PWR half-cell geometry: fuel | clad | coolant.

        The symmetry of a square PWR unit cell about the fuel centreline:
        the slab starts at the symmetry plane :math:`x = 0` and crosses
        half the fuel, the cladding and the coolant to the unit-cell
        boundary. Material ids fuel ``2``, cladding ``1``, coolant ``0``.
        Both faces are reflective, and that is part of the model, not a
        default: the fuel centreline and the unit-cell boundary are both
        symmetry planes of the infinite lattice the half cell stands for.
        The same stack under other laws (an isolated cell, say, white on
        the coolant face) is a different body, built with :meth:`slab`.

        Parameters
        ----------
        fuel_half, clad_thick, cool_thick : float
            Half the fuel thickness, the cladding and the coolant
            thicknesses (cm).
        """
        return cls.from_thicknesses(
            coord=CoordSystem.CARTESIAN,
            thicknesses=(fuel_half, clad_thick, cool_thick),
            mat_ids=(2, 1, 0),
            boundaries=(BC("reflective"), BC("reflective")),
        )


def _parse_coord(coord: object) -> CoordSystem:
    """The coordinate system, or the keyed refusal of anything else."""
    if not isinstance(coord, CoordSystem):
        raise TypeError(
            f"StructuredGeometry.coord must be a CoordSystem member, got "
            f"{type(coord).__name__} {coord!r}; the string "
            f"kind tags ('SLB', 'CYL', 'SPH') are retired."
        )
    return coord


def _is_hollow(coord: CoordSystem, breakpoints: tuple[float, ...]) -> bool:
    """A cylinder or sphere whose first (parsed) breakpoint is positive."""
    return coord is not CoordSystem.CARTESIAN and breakpoints[0] > 0.0


def _boundary_points(
    coord: CoordSystem, breakpoints: tuple[float, ...],
) -> tuple[float, ...]:
    """The boundary's positions, inner first, from parsed breakpoints."""
    return coord.boundary_points(breakpoints[0], breakpoints[-1])


def _parse_breakpoints(
    coord: CoordSystem, breakpoints: object,
) -> tuple[float, ...]:
    """The breakpoints as a tuple of ``float``, or a keyed refusal."""
    entries = parse_entries(breakpoints, "StructuredGeometry.breakpoints", "real numbers")
    parsed = tuple(
        parse_real(value, f"StructuredGeometry.breakpoints[{k}]")
        for k, value in enumerate(entries)
    )
    if len(parsed) < 2:
        raise ValueError(
            f"StructuredGeometry needs at least 2 breakpoints (one interval); "
            f"got {len(parsed)}"
        )
    if not all(math.isfinite(value) for value in parsed):
        raise ValueError(
            f"StructuredGeometry.breakpoints must be finite; got {parsed}"
        )
    if any(b <= a for a, b in itertools.pairwise(parsed)):
        raise ValueError(
            f"StructuredGeometry.breakpoints must be strictly increasing "
            f"(every interval has positive width); got {parsed}"
        )
    if coord is not CoordSystem.CARTESIAN and parsed[0] < 0.0:
        raise ValueError(
            f"a radial coordinate starts at r_0 >= 0: the {coord.name.lower()} "
            f"geometry's first breakpoint is {parsed[0]!r}"
        )
    return parsed


def _parse_mat_ids(mat_ids: object, n_intervals: int) -> tuple[int, ...]:
    """The material ids as a tuple of ``int``, one per interval."""
    entries = tuple(
        parse_integer(value, f"StructuredGeometry.mat_ids[{k}]", "a material id")
        for k, value in enumerate(parse_entries(mat_ids, "StructuredGeometry.mat_ids", "int"))
    )
    if len(entries) != n_intervals:
        raise ValueError(
            f"StructuredGeometry takes one material id per interval: "
            f"{n_intervals} interval(s), {len(entries)} material id(s)"
        )
    return entries


def _check_boundaries(geometry: StructuredGeometry) -> None:
    """One law per boundary point, each a ``BC`` tag or a typed law."""
    parse_boundary_laws(
        geometry.boundaries, geometry.coord, geometry.breakpoints[0],
        geometry.breakpoints[-1], "StructuredGeometry.boundaries",
    )


def parse_boundary_law(law: object, where: str) -> "BC | BoundaryTraceLaw":
    """One boundary law: a ``BC`` tag or a typed ``BoundaryTraceLaw``, never ``None``.

    The one check of a single boundary declaration, shared by the geometry,
    both meshes and both axis classes.

    **Why the law arm exists.** A ``BC`` tag is
    ``(kind: str, params: dict[str, float])``, structurally unable to carry a
    law whose content is a FUNCTION: a
    :class:`~orpheus.geometry.boundary.PrescribedInflow` whose source is a
    manufactured solution restricted to a face has no tag spelling. Declaring
    such a law on the mesh is what lets it survive the public solver entries,
    which rebuild the method mesh from the declaration: the shared
    :func:`~orpheus.transport.method.resolve_boundary_conditions` reads the
    declared law, so a law reaches the realizer through the path a tag does.
    The tag remains the spelling for everything expressible as one: it is
    serialisable, comparable, and the input-deck surface.
    """
    # Imported lazily: ``orpheus.geometry.boundary`` transitively loads
    # THIS module, so a top-level import cycles.
    from orpheus.geometry.boundary import BoundaryTraceLaw

    if law is None:
        raise TypeError(
            f"{where} is None, and None is not a boundary law: "
            f"declare the law the boundary carries."
        )
    if not isinstance(law, (BC, BoundaryTraceLaw)):
        raise TypeError(
            f"{where} must be a BC tag or a BoundaryTraceLaw "
            f"instance, got {type(law).__name__}"
        )
    return law


def parse_boundary_laws(
    laws: "tuple[object, ...]",
    coord: CoordSystem,
    r_0: float,
    r_R: float,
    where: str,
) -> "tuple[BC | BoundaryTraceLaw, ...]":
    """One law per boundary point of ``[r_0, r_R]`` in ``coord``, or a keyed refusal.

    The check of the geometry's boundary declaration: each law passes
    :func:`parse_boundary_law`, and there is one per
    :meth:`CoordSystem.boundary_points`. A mesh's face laws are checked by
    :meth:`~orpheus.mesh.face_laws.FaceLaws.over` against the same
    boundary points, named.
    """
    parsed = tuple(
        parse_boundary_law(law, f"{where}[{k}]") for k, law in enumerate(laws)
    )
    points = coord.boundary_points(r_0, r_R)
    if len(parsed) == len(points):
        return parsed
    if coord is CoordSystem.CARTESIAN:
        reason = "a slab has two boundary points (left, right)"
    elif len(points) == 1:
        reason = (
            f"the centre r = 0 of a solid {coord.name.lower()} body is an "
            f"interior point and carries no law; the only boundary point is "
            f"the outer surface"
        )
    else:
        reason = (
            f"a hollow {coord.name.lower()} body (r_0 = {points[0]!r} > 0) "
            f"has an inner surface, which needs its own law"
        )
    raise ValueError(
        f"{where}: {reason}; expected {len(points)} law(s) at r = {points}, "
        f"got {len(laws)}."
    )


__all__ = [
    "StructuredGeometry",
    "parse_boundary_law",
    "parse_boundary_laws",
]
