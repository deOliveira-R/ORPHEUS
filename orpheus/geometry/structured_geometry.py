r"""Structured 1-D geometry — pure geometry layer, no mesh concerns.

A :class:`StructuredGeometry` is an interval of positions in one
coordinate system, cut into material intervals, with a role at every
point of its boundary:

* **the coordinate system** (:class:`CoordSystem`: slab, cylinder,
  sphere), which fixes the measure of an interval and the topology of
  its boundary;
* **the breakpoints** :math:`r_0 < r_1 < \dots < r_R`, the positions
  where one material interval ends and the next begins; they are stored
  as given and never re-derived;
* **one material id per interval** :math:`[r_k, r_{k+1}]`;
* **one face role per boundary point**, in the order (inner, outer): a
  deck (the domain continues across the face: a mirror, a translation) or
  an ``int`` tag (the domain ends there; the problem binds the tag to an
  exterior condition where the method is known). The geometry names no
  response.

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
boundary in its coordinate system, so the number of roles a geometry
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
chart, and no role is declared at it. A hollow body
(:math:`r_0 > 0`) has an inner surface, which carries a role like any
other boundary point, a tag and never a deck (the quotient
:math:`W = N(G_c)/G_c` of a 1-D radial chart is trivial, so a reflecting
curved wall is a response a problem binds to a tag).
:attr:`StructuredGeometry.boundary_points` lists the positions, paired
one to one with :attr:`StructuredGeometry.face_roles`.

Architectural role
------------------

* **Reference solution generators** (``MomentSpace``,
  ``Spectrum``, ``BasisSpace``) consume a :class:`StructuredGeometry`
  directly: no mesh, no cell counts.
* **Discrete production solvers** (``solve_cp``, ``solve_sn``,
  ``solve_moc``, ``solve_mc``) consume a :class:`~orpheus.mesh.Mesh1D`
  built from the geometry; the discretisation is supplied at that build
  step, never stored on the geometry.

Examples
--------

A bare-critical sphere, one region, one outer tag::

    geom = StructuredGeometry(
        coord=CoordSystem.SPHERICAL,
        breakpoints=(0.0, 2.872),
        mat_ids=(0,),
        face_roles=(0,),
    )

A reflected slab (reflector | core | reflector), two tags::

    geom = StructuredGeometry.from_thicknesses(
        coord=CoordSystem.CARTESIAN,
        thicknesses=(0.5, 2.0, 0.5),
        mat_ids=(1, 0, 1),
        face_roles=(0, 0),
    )

A hollow sphere, an inner and an outer tag::

    geom = StructuredGeometry(
        coord=CoordSystem.SPHERICAL,
        breakpoints=(0.5, 1.0, 2.0),
        mat_ids=(1, 0),
        face_roles=(1, 0),
    )

A PWR pin cell through the Wigner–Seitz factory::

    geom = StructuredGeometry.wigner_seitz_pin_cell(
        r_fuel=0.9, r_clad=1.1, pitch=3.6,
    )
"""
from __future__ import annotations

import itertools
from collections.abc import Iterable
from dataclasses import dataclass

import numpy as np

from orpheus.values.content import ContentIdentity

from .coord import CoordSystem
from orpheus.values.scalars import parse_entries, parse_integer, parse_positions, parse_positive_real, parse_real
from .colouring import ColouringSymmetry, colouring_symmetry
from .deck import DeckTransformation, FaceRole, MirrorDeck, parse_face_roles
from .faces import face_name


@dataclass(frozen=True, kw_only=True, eq=False)
class StructuredGeometry(ContentIdentity):
    r"""A 1-D geometry: a coordinate system, breakpoints, one material id
    per interval and one face role per boundary point.

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
    face_roles : tuple of FaceRole
        One role per boundary point, paired one to one with
        :attr:`boundary_points` (inner first). A role is a deck
        (:class:`~orpheus.geometry.deck.MirrorDeck`,
        :class:`~orpheus.geometry.deck.TranslationDeck`: the domain
        continues across the face) or an ``int`` tag (the domain ends there;
        a ``bool`` is not a tag), the exact twin of a material id. A deck on
        a face of a 1-D radial chart is a
        :class:`~orpheus.geometry.deck.DeckError`, and a translation needs a
        translation on the opposite face. ``None`` is not a role.

        Notes
    -----
    **What is not here.** No cell counts and no discretisation rule (the
    mesh layer's); no critical dimension (a registry's truth record); no
    group count (the materials'); no infinite medium (a problem with no
    geometry, or a finite domain with reflective laws).

    **The geometry names no response.** A tag says which faces share an
    exterior model and never what the model is; a binding, supplied where
    the method is known, maps the tags to conditions.

    **Identity is content** (:class:`~orpheus.values.content.ContentIdentity`):
    two geometries with equal coordinate systems, breakpoints, material ids
    and face roles are equal and hash alike in every process, so a reference
    cache can key on one. A deck and a tag are different declarations
    (``MirrorDeck()`` is not tag ``1``), and the faces' order is content
    (inner first).
    """

    coord: CoordSystem
    breakpoints: tuple[float, ...]
    mat_ids: tuple[int, ...]
    face_roles: tuple[FaceRole, ...]

    def __post_init__(self) -> None:
        _parse_coord(self.coord)
        object.__setattr__(
            self, "breakpoints", _parse_breakpoints(self.coord, self.breakpoints),
        )
        object.__setattr__(
            self, "mat_ids", _parse_mat_ids(self.mat_ids, len(self.breakpoints) - 1),
        )
        object.__setattr__(self, "face_roles", _parse_face_roles(self))

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
        r"""The positions of the boundary, inner first; one role each.

        :math:`(r_0, r_R)` on a slab and on a hollow cylinder or sphere;
        :math:`(r_R,)` on a solid cylinder or sphere, whose centre is an
        interior point.
        """
        return _boundary_points(self.coord, self.breakpoints)

    @property
    def tags(self) -> frozenset[int]:
        r"""The tags its faces carry (the decks excluded)."""
        return frozenset(role for role in self.face_roles if not isinstance(role, DeckTransformation))

    @property
    def symmetry(self) -> ColouringSymmetry:
        r"""The colouring symmetry: the part of the mirror group that the material runs and the end roles respect.

        Computed from this geometry's breakpoints, material ids and end roles (never from a mesh's cells), by
        :func:`~orpheus.geometry.colouring.colouring_symmetry`, so one ``(geometry, binding)`` pair is admitted or
        refused for every mesh of it. A radial chart answers the identity only.
        """
        return colouring_symmetry(
            self.coord, (self.breakpoints[0], self.breakpoints[-1]), self.breakpoints, self.mat_ids,
            self.face_roles[0], self.face_roles[-1],
        )

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
    # Constructors that name each boundary point's role
    # ─────────────────────────────────────────────────────────────────

    @classmethod
    def slab(
        cls,
        breakpoints: Iterable[float],
        mat_ids: Iterable[int],
        *,
        left: FaceRole,
        right: FaceRole,
    ) -> "StructuredGeometry":
        r"""A slab, with the role at its left and at its right face."""
        return cls(
            coord=CoordSystem.CARTESIAN, breakpoints=tuple(breakpoints),
            mat_ids=tuple(mat_ids), face_roles=(left, right),
        )

    @classmethod
    def cylinder(
        cls,
        breakpoints: Iterable[float],
        mat_ids: Iterable[int],
        *,
        outer: FaceRole,
        inner: FaceRole | None = None,
    ) -> "StructuredGeometry":
        r"""A cylinder, with the role at its outer surface, and at its inner one when hollow.

        ``inner`` is given exactly when the first breakpoint is positive: the
        centre of a solid cylinder is an interior point and carries no role.
        A curved face carries a tag, never a deck.
        """
        return cls._radial(CoordSystem.CYLINDRICAL, breakpoints, mat_ids, inner, outer)

    @classmethod
    def sphere(
        cls,
        breakpoints: Iterable[float],
        mat_ids: Iterable[int],
        *,
        outer: FaceRole,
        inner: FaceRole | None = None,
    ) -> "StructuredGeometry":
        r"""A sphere, with the role at its outer surface, and at its inner one when hollow.

        ``inner`` is given exactly when the first breakpoint is positive: the
        centre of a solid sphere is an interior point and carries no role.
        A curved face carries a tag, never a deck.
        """
        return cls._radial(CoordSystem.SPHERICAL, breakpoints, mat_ids, inner, outer)

    @classmethod
    def _radial(
        cls,
        coord: CoordSystem,
        breakpoints: Iterable[float],
        mat_ids: Iterable[int],
        inner: FaceRole | None,
        outer: FaceRole,
    ) -> "StructuredGeometry":
        # An absent ``inner`` keyword is no role at all, never a stored None;
        # whether the body needs one is decided by the one boundary check.
        roles = (outer,) if inner is None else (inner, outer)
        return cls(
            coord=coord, breakpoints=tuple(breakpoints),
            mat_ids=tuple(mat_ids), face_roles=roles,
        )

    @classmethod
    def uniform_boundary(
        cls,
        coord: CoordSystem,
        breakpoints: Iterable[float],
        mat_ids: Iterable[int],
        role: FaceRole,
    ) -> "StructuredGeometry":
        r"""The geometry whose every boundary point carries ``role``.

        Coordinate-generic: two roles on a slab or a hollow body, one on a
        solid cylinder or sphere. For a test parametrised over coordinate
        systems whose body has one role on its whole boundary.
        """
        parsed = _parse_breakpoints(_parse_coord(coord), breakpoints)
        return cls(
            coord=coord, breakpoints=parsed, mat_ids=tuple(mat_ids),
            face_roles=(role,) * len(_boundary_points(coord, parsed)),
        )

    @classmethod
    def from_homogeneous(
        cls, width: float, boundary: FaceRole,
    ) -> "StructuredGeometry":
        r"""A slab ``[0, width]`` of material 0, with the role ``boundary`` on both faces.

        The infinite medium as a test realises it: one finite slab whose
        faces carry one role (a mirror deck, for the infinite medium itself).
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
        face_roles: tuple[FaceRole, ...],
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
            face_roles=face_roles,
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
        ids fuel ``2``, cladding ``1``, coolant ``0``. The outer face is tag ``0``.
        The Wigner–Seitz model binds tag ``0`` to ``WallCondition.white()``:
        the cell IS the cylindricalised lattice cell with isotropic re-entry
        at its outer surface, the assumption that maps a periodic lattice to
        one cell. The geometry cannot hold that (it names no response, and
        there is no problem-level preset type, #529), so the model is stated
        here and bound where the method is known. The same radii under
        another condition are the same body with another binding; a
        different body is built with :meth:`cylinder`.

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
            face_roles=(0,),
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
        Both faces are mirror decks, and that is part of the model, not a
        default: the fuel centreline and the unit-cell boundary are both
        symmetry planes of the infinite lattice the half cell stands for.
        The same stack under other faces (an isolated cell, say, a tag on
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
            face_roles=(MirrorDeck(), MirrorDeck()),
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
    parsed = tuple(float(b) for b in parse_positions(breakpoints, "StructuredGeometry.breakpoints"))
    if len(parsed) < 2:
        raise ValueError(
            f"StructuredGeometry needs at least 2 breakpoints (one interval); "
            f"got {len(parsed)}"
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


def _parse_face_roles(geometry: StructuredGeometry) -> tuple[FaceRole, ...]:
    """One role per boundary point, or a keyed refusal; a deck only where the chart admits one.

    The count rule is the boundary points' (a slab two, a hollow body two, a solid body one: the centre of a solid
    cylinder or sphere is an interior point and carries no role); each role is parsed, a deck on a radial face is a
    :class:`~orpheus.geometry.deck.DeckError`, and the pairings close, all through
    :func:`~orpheus.geometry.deck.parse_face_roles`.
    """
    where = "StructuredGeometry.face_roles"
    roles = parse_entries(geometry.face_roles, where, "face roles")
    points = geometry.boundary_points
    if len(roles) != len(points):
        if geometry.coord is CoordSystem.CARTESIAN:
            reason = "a slab has two boundary points (left, right)"
        elif len(points) == 1:
            reason = (
                f"the centre r = 0 of a solid {geometry.coord.name.lower()} body is an "
                f"interior point and carries no role; the only boundary point is "
                f"the outer surface"
            )
        else:
            reason = (
                f"a hollow {geometry.coord.name.lower()} body (r_0 = {points[0]!r} > 0) "
                f"has an inner surface, which needs its own role"
            )
        raise ValueError(f"{where}: {reason}; expected {len(points)} role(s) at r = {points}, got {len(roles)}.")
    faces = tuple(face_name(0, sign) for sign in (-1, +1))[-len(points):]
    return tuple(parse_face_roles(dict(zip(faces, roles, strict=True)), geometry.coord, where).values())


__all__ = [
    "StructuredGeometry",
]
