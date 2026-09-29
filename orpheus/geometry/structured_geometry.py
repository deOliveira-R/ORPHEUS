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
from numbers import Real
from typing import TYPE_CHECKING

import numpy as np

from .coord import CoordSystem
from .boundary import BC

if TYPE_CHECKING:
    from .boundary import BoundaryTraceLaw


def _parse_real(value: object, where: str) -> float:
    """A real scalar as a ``float``, or a keyed refusal.

    ``bool`` is refused (``True`` is an ``int``), and ``-0.0`` is
    canonicalised to ``+0.0``: the two compare equal, so they are one
    breakpoint, and a digest over the bits must see one value.
    """
    if not isinstance(value, Real) or isinstance(value, bool):
        raise TypeError(
            f"{where} must be a real number, got {type(value).__name__}"
        )
    return float(value) + 0.0


def _parse_integer(value: object, where: str) -> int:
    """An integer scalar as an ``int``, or a keyed refusal."""
    if not isinstance(value, (int, np.integer)) or isinstance(value, bool):
        raise TypeError(
            f"{where}: a material id is an int, got {type(value).__name__}"
        )
    return int(value)


def _entries(value: object, where: str, expected: str) -> tuple[object, ...]:
    """The entries of a sequence field as a tuple, or a keyed refusal."""
    if isinstance(value, (str, bytes)) or not isinstance(value, Iterable):
        raise TypeError(f"{where} must be a sequence of {expected}, got {type(value).__name__}")
    return tuple(value)


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
        if not isinstance(self.coord, CoordSystem):
            raise TypeError(
                f"StructuredGeometry.coord must be a CoordSystem member, got "
                f"{type(self.coord).__name__} {self.coord!r}; the string "
                f"kind tags ('SLB', 'CYL', 'SPH') are retired."
            )
        object.__setattr__(
            self, "breakpoints", _parse_breakpoints(self.coord, self.breakpoints),
        )
        object.__setattr__(
            self, "mat_ids", _parse_mat_ids(self.mat_ids, len(self.breakpoints) - 1),
        )
        object.__setattr__(
            self, "boundaries",
            _entries(self.boundaries, "StructuredGeometry.boundaries", "boundary laws"),
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
        return self.coord is not CoordSystem.CARTESIAN and self.breakpoints[0] > 0.0

    @property
    def boundary_points(self) -> tuple[float, ...]:
        r"""The positions of the boundary, inner first; one law each.

        :math:`(r_0, r_R)` on a slab and on a hollow cylinder or sphere;
        :math:`(r_R,)` on a solid cylinder or sphere, whose centre is an
        interior point.
        """
        r_0, r_R = self.breakpoints[0], self.breakpoints[-1]
        if self.coord is CoordSystem.CARTESIAN or self.is_hollow:
            return (r_0, r_R)
        return (r_R,)

    @property
    def domain_extent_cm(self) -> float:
        r"""The width :math:`r_R - r_0` of the interval of positions, in cm.

        The full slab width on a slab (the production convention: a user
        with "a slab of thickness L" means L); the outer radius on a
        solid cylinder or sphere; the shell thickness on a hollow one.
        """
        return self.breakpoints[-1] - self.breakpoints[0]

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
                _parse_real(t, f"StructuredGeometry.from_thicknesses: thicknesses[{k}]")
                for k, t in enumerate(thicknesses)
            ),
            initial=_parse_real(r_0, "StructuredGeometry.from_thicknesses: r_0"),
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
        boundaries: "tuple[BC | BoundaryTraceLaw, ...]" = (BC("white"),),
    ) -> "StructuredGeometry":
        r"""Wigner–Seitz equivalent pin-cell geometry.

        Replaces a square unit cell of side ``pitch`` with a cylinder of
        equal cross-sectional area,

        .. math::

            r_{\rm cell} = \frac{\rm pitch}{\sqrt{\pi}},

        cut at the fuel and cladding radii: breakpoints
        :math:`(0, r_{\rm fuel}, r_{\rm clad}, r_{\rm cell})`, material
        ids fuel ``2``, cladding ``1``, coolant ``0``. The default outer
        law is ``BC("white")``, the unit-cell symmetry assumption that
        maps a periodic lattice to one cell with isotropic re-entry.

        Parameters
        ----------
        r_fuel, r_clad : float
            The fuel-pellet and cladding outer radii (cm).
        pitch : float
            The square unit cell's side (cm).
        boundaries : tuple
            The law at the outer surface, a 1-tuple (the cylinder is
            solid).
        """
        r_cell = float(pitch / np.sqrt(np.pi))
        return cls(
            coord=CoordSystem.CYLINDRICAL,
            breakpoints=(0.0, float(r_fuel), float(r_clad), r_cell),
            mat_ids=(2, 1, 0),
            boundaries=boundaries,
        )

    @classmethod
    def pwr_slab_half_cell(
        cls,
        *,
        fuel_half: float = 0.9,
        clad_thick: float = 0.2,
        cool_thick: float = 0.7,
        boundaries: "tuple[BC | BoundaryTraceLaw, ...]" = (
            BC("reflective"), BC("reflective"),
        ),
    ) -> "StructuredGeometry":
        r"""Cartesian 1-D PWR half-cell geometry: fuel | clad | coolant.

        The symmetry of a square PWR unit cell about the fuel centreline:
        the slab starts at the symmetry plane :math:`x = 0` and crosses
        half the fuel, the cladding and the coolant to the unit-cell
        boundary. Material ids fuel ``2``, cladding ``1``, coolant ``0``.
        The default laws are reflective on both faces, the infinite-
        lattice convention; ``(reflective, white)`` models an isolated
        cell with isotropic re-entry on the coolant face.

        Parameters
        ----------
        fuel_half, clad_thick, cool_thick : float
            Half the fuel thickness, the cladding and the coolant
            thicknesses (cm).
        boundaries : tuple
            ``(left, right)`` laws.
        """
        return cls.from_thicknesses(
            coord=CoordSystem.CARTESIAN,
            thicknesses=(fuel_half, clad_thick, cool_thick),
            mat_ids=(2, 1, 0),
            boundaries=boundaries,
        )


def _parse_breakpoints(
    coord: CoordSystem, breakpoints: object,
) -> tuple[float, ...]:
    """The breakpoints as a tuple of ``float``, or a keyed refusal."""
    entries = _entries(breakpoints, "StructuredGeometry.breakpoints", "real numbers")
    parsed = tuple(
        _parse_real(value, f"StructuredGeometry.breakpoints[{k}]")
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
        _parse_integer(value, f"StructuredGeometry.mat_ids[{k}]")
        for k, value in enumerate(_entries(mat_ids, "StructuredGeometry.mat_ids", "int"))
    )
    if len(entries) != n_intervals:
        raise ValueError(
            f"StructuredGeometry takes one material id per interval: "
            f"{n_intervals} interval(s), {len(entries)} material id(s)"
        )
    return entries


def _check_boundaries(geometry: StructuredGeometry) -> None:
    """One law per boundary point, each a ``BC`` tag or a typed law."""
    # Imported lazily: ``orpheus.geometry.boundary`` transitively loads
    # THIS module, so a top-level import cycles.
    from orpheus.geometry.boundary import BoundaryTraceLaw

    boundaries = geometry.boundaries
    for k, law in enumerate(boundaries):
        if law is None:
            raise TypeError(
                f"StructuredGeometry.boundaries[{k}] is None, and None is not "
                f"a boundary law: declare the law the boundary point carries."
            )
        if not isinstance(law, (BC, BoundaryTraceLaw)):
            raise TypeError(
                f"StructuredGeometry.boundaries[{k}] must be a BC tag or a "
                f"BoundaryTraceLaw instance, got {type(law).__name__}"
            )
    points = geometry.boundary_points
    if len(boundaries) == len(points):
        return
    if geometry.coord is CoordSystem.CARTESIAN:
        reason = "a slab has two boundary points (left, right)"
    elif not geometry.is_hollow:
        reason = (
            f"the centre r = 0 of a solid {geometry.coord.name.lower()} "
            f"geometry is an interior point and carries no law; the only "
            f"boundary point is the outer surface"
        )
    else:
        reason = (
            f"a hollow {geometry.coord.name.lower()} geometry (r_0 = "
            f"{points[0]!r} > 0) has an inner surface, which needs its own law"
        )
    raise ValueError(
        f"StructuredGeometry: {reason}; expected {len(points)} law(s) at "
        f"r = {points}, got {len(boundaries)}."
    )


__all__ = [
    "StructuredGeometry",
]
