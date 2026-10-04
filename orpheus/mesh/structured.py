"""Mesh data structures for 1-D and 2-D geometries.

A mesh is an immutable description of a spatial domain: cell edges,
material assignments, and derived quantities (volumes, areas,
cell centres).  Solvers receive a mesh and build mutable,
solver-specific state on top of it.

Both :class:`Mesh1D` and :class:`Mesh2D` are frozen dataclasses --
once created, their fields cannot be reassigned.

Mesh construction from geometry
-------------------------------

A :class:`Mesh1D` is built by a :class:`~orpheus.mesh.mesher.Mesher`: load a
:class:`~orpheus.geometry.StructuredGeometry`, partition it by interval rules
(:mod:`orpheus.mesh.partition`), read the mesh::

    mesh = Mesher(geometry).partition(CellsByCount.uniform_width(8)).mesh

The mesh itself holds no geometry: its coordinate system, cells, the material
of each cell and the law on each boundary face are all it is. The mesher lifts
the geometry's regions and boundary laws onto the cells and faces.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import KW_ONLY, dataclass, field, replace
from typing import TYPE_CHECKING, assert_never, cast

import numpy as np

from orpheus.geometry.boundary import BC, BoundaryTraceLaw
from orpheus.geometry.coord import (
    CoordSystem,
    compute_areas_1d,
    compute_volumes_2d,
)
from orpheus.numerics.scalars import parse_entries, parse_index, parse_integer, parse_positions
from orpheus.mesh.face_laws import FaceLaws, face_inventory
from orpheus.numerics.content import ContentIdentity
from orpheus.numerics.mesh_free_function import MeshFreeFunction, RegionwiseConstant, Symbolic

if TYPE_CHECKING:
    import sympy

# ═══════════════════════════════════════════════════════════════════════
# Mesh1D
# ═══════════════════════════════════════════════════════════════════════

def _volume_ulps(coord: CoordSystem) -> int:
    r"""How far a stored cell volume may sit from its cell's measure, in ulp of :math:`c\,T(r_{j+1})`.

    The stored volume of an equal-measure cell is the share :math:`m/n`; the
    measure it is checked against is the shell :math:`c\,(T(r_{j+1}) - T(r_j))`
    between the realised edges, :math:`T(r) = r^p`. Each realised edge carries
    at most 1 ulp from forming :math:`t_j = T(a) + f_j\,\Delta T`, the root's
    rounding (at most 1 ulp, allowing a :math:`\sqrt[3]{\cdot}` that is not
    correctly rounded) amplified by :math:`p` in :math:`T(r)`, and ½ ulp from
    re-evaluating :math:`T`: :math:`p + 1.5` ulp per edge. Two edges, plus the
    subtraction, the constant and the share (about 2 ulp together), give
    :math:`2p + 5`: 7 on a slab, 9 on a cylinder, 11 on a sphere. `[M]` The
    largest gap over about 46 000 random legal equal-volume partitions is 3.05,
    6.53 and 7.62 ulp (the elegance review of P1 step 3b, macOS libm); a volume
    from the wrong coordinate system is about 6e15 ulp off. A volume given with
    its edges (``CellEdges``, a shell) is the measure itself: 0 ulp.
    """
    return 2 * coord.measure_coordinate.exponent + 5



def _require_polynomial_steps(expression: "sympy.Expr", r: "sympy.Symbol", group: int) -> None:
    """Refuse a non-smooth construct whose argument is not polynomial in ``r`` (see :meth:`Mesh1D.cell_integrals`)."""
    import sympy

    arguments: list[sympy.Basic] = []
    for piecewise in expression.atoms(sympy.Piecewise):
        for branch in piecewise.args:
            condition = branch.args[1]
            for relation in condition.atoms(sympy.core.relational.Relational):
                arguments += [relation.lhs, relation.rhs]  # both sides polynomial: the step's boundary is a root set
    for kind in (sympy.Heaviside, sympy.Abs, sympy.Max, sympy.Min, sympy.sign):
        for atom in expression.atoms(kind):
            arguments += list(atom.args)
    for argument in arguments:
        if not cast("sympy.Expr", argument).is_polynomial(r):
            raise ValueError(
                f"Mesh1D.cell_integrals: group {group}'s weight steps where {argument} changes sign; a step is "
                f"integrated only where its argument is polynomial in r (a scope boundary of the exact integration)"
            )


@dataclass(frozen=True, eq=False)
class Mesh1D(ContentIdentity):
    r"""A 1-D mesh: cells in one coordinate system, a region label on each, a law on each boundary face.

    The general constructor takes exactly what a 1-D discretisation is. A
    mesh is normally built by a :class:`~orpheus.mesh.mesher.Mesher` from a
    :class:`~orpheus.geometry.StructuredGeometry` and interval rules; the
    constructor is what the mesher, a relabelling or a refinement calls.

    Parameters
    ----------
    coord : CoordSystem
        The coordinate system; it gives the measure and the face areas.
    edges : array_like, shape (N+1,)
        The cell edges, strictly increasing and finite; :math:`r_0 \ge 0` on a
        cylinder or a sphere.
    volumes : array_like, shape (N,)
        The cell volumes (lengths on a slab, areas per unit height on a
        cylinder). Stored, not recomputed: an equal share :math:`m/n` of an
        interval is exact, and the shell between realised edges is not
        (ERR-020). Each must be the coordinate system's measure of its cell
        within :func:`_volume_ulps` ulp.
    region_ids : array_like of int, shape (N,)
        The region label of each cell, ``0 .. R-1``: the partition the mesh
        refines, which an external mesh carries as its cell labels (Gmsh's
        physical volume groups) and the :class:`~orpheus.mesh.mesher.Mesher`
        writes as each cell's interval. Required; an unlabelled mesh is one
        region. Labels are positions, the convention of a
        :class:`~orpheus.numerics.mesh_free_function.RegionwiseConstant`'s
        rows, so a region table reads out to the cells as :meth:`per_cell`.
    region_materials : tuple of int, length R
        The material id of each region, indexed by label; every entry labels
        some cell. The cells' material ids are derived from it
        (:attr:`mat_ids`, through :meth:`per_cell`), never given beside it. It moves to the system
        with the posing work (#522); the mesh does not hold its geometry,
        because an external mesh has none (the user's ruling of 2026-10-03).
    face_laws : mapping of face name to BC or BoundaryTraceLaw
        The law on each boundary face, over exactly the mesh's
        :func:`~orpheus.mesh.face_laws.face_inventory`: ``xmin`` and ``xmax``
        on a slab or a hollow body, ``xmax`` alone on a solid cylinder or
        sphere, whose centre carries no law. Stored as a
        :class:`~orpheus.mesh.face_laws.FaceLaws`, the value
        :class:`Mesh2D` stores too.

    Notes
    -----
    **Identity is content** (:class:`~orpheus.numerics.content.ContentIdentity`,
    #405 P1 step 5, 2026-10-02): the coordinate system, edges, volumes,
    region labels, the region → material map and the face laws (the
    material ASSIGNMENT, never the materials' data), so two meshes with the same cells, labels, map and laws
    are equal and hash alike in every process; this is the discretisation
    digest the reference cache keys on, and it never reads the materials.
    Until then the mesh compared bitwise by a hand-written ``__eq__`` and
    could not be hashed.
    """

    coord: CoordSystem
    edges: np.ndarray
    volumes: np.ndarray
    region_ids: np.ndarray
    region_materials: tuple[int, ...]
    face_laws: "Mapping[str, BC | BoundaryTraceLaw]"

    # Derived from the content, so not content (``compare=False``).
    mat_ids: np.ndarray = field(init=False, repr=False, compare=False)
    widths: np.ndarray = field(init=False, repr=False, compare=False)
    centers: np.ndarray = field(init=False, repr=False, compare=False)
    areas: np.ndarray = field(init=False, repr=False, compare=False)

    def __post_init__(self) -> None:
        if not isinstance(self.coord, CoordSystem):
            raise TypeError(
                f"Mesh1D.coord is a CoordSystem member, got {type(self.coord).__name__}"
            )
        edges = parse_positions(self.edges, "Mesh1D.edges")
        if len(edges) < 2 or not np.all(np.diff(edges) > 0):
            raise ValueError(
                f"Mesh1D.edges are at least two strictly increasing positions; got {edges}"
            )
        if self.coord is not CoordSystem.CARTESIAN and edges[0] < 0.0:
            raise ValueError(
                f"Mesh1D.edges: a radial coordinate starts at r_0 >= 0, got {edges[0]!r}"
            )
        n = len(edges) - 1
        volumes = parse_positions(self.volumes, "Mesh1D.volumes")
        if volumes.shape != (n,):
            raise ValueError(
                f"Mesh1D.volumes: {n} cell(s) need {n} volume(s), got shape {volumes.shape}"
            )
        if not np.all(volumes > 0.0):
            j = int(np.argmin(volumes))
            raise ValueError(
                f"Mesh1D.volumes[{j}] = {volumes[j]!r}: a cell volume is positive"
            )
        shells = self.coord.measure(edges)
        ulps = np.spacing(self.coord.measure_constant * self.coord.measure_coordinate(edges[1:]))
        off = np.abs(volumes - shells) / ulps
        band = _volume_ulps(self.coord)
        if not np.all(off <= band):
            j = int(np.argmax(off))
            raise ValueError(
                f"Mesh1D.volumes[{j}] = {volumes[j]!r} is not the {self.coord.name.lower()} "
                f"measure of the cell [{edges[j]!r}, {edges[j + 1]!r}], {shells[j]!r} "
                f"({off[j]:.3g} ulp off; at most {band})"
            )
        region_ids = np.array(
            [parse_index(r, f"Mesh1D.region_ids[{k}]", "a region label")
             for k, r in enumerate(parse_entries(self.region_ids, "Mesh1D.region_ids", "int"))],
            dtype=int,
        )
        if region_ids.shape != (n,):
            raise ValueError(
                f"Mesh1D.region_ids: {n} cell(s) need {n} region label(s), got {len(region_ids)}"
            )
        region_materials = tuple(
            parse_integer(m, f"Mesh1D.region_materials[{k}]", "a material id")
            for k, m in enumerate(parse_entries(
                self.region_materials, "Mesh1D.region_materials",
                "material ids indexed by region label (a tuple of material ids)",
            ))
        )
        cells_per_region = np.bincount(region_ids, minlength=len(region_materials))
        if len(cells_per_region) > len(region_materials):
            raise ValueError(
                f"Mesh1D.region_ids: the label {len(cells_per_region) - 1} has no material "
                f"({len(region_materials)} region(s) in region_materials)"
            )
        if not np.all(cells_per_region):
            raise ValueError(
                f"Mesh1D.region_materials: the region {int(np.argmin(cells_per_region))} labels no cell; "
                f"every region is some cells' label"
            )
        face_laws = FaceLaws.over(
            face_inventory(self.coord, edges, 1), self.face_laws,
            "Mesh1D.face_laws", self.coord,
        )
        region_ids.flags.writeable = False
        widths = np.diff(edges)
        centers = 0.5 * (edges[:-1] + edges[1:])
        areas = compute_areas_1d(self.coord, edges)
        for name, value in (
            ("edges", edges), ("volumes", volumes), ("region_ids", region_ids),
            ("region_materials", region_materials), ("face_laws", face_laws), ("widths", widths), ("centers", centers),
            ("areas", areas),
        ):
            object.__setattr__(self, name, value)
        mat_ids = self.per_cell(np.asarray(region_materials, dtype=int))
        mat_ids.flags.writeable = False
        object.__setattr__(self, "mat_ids", mat_ids)


    # ── Derived properties ────────────────────────────────────────────

    @property
    def N(self) -> int:
        """Number of cells."""
        return len(self.edges) - 1

    def per_cell(self, region_values: "np.typing.ArrayLike") -> np.ndarray:
        """A per-region table read out to the cells: row ``r`` of ``region_values`` for every cell labelled ``r``."""
        return np.asarray(region_values)[self.region_ids]

    def cell_integrals(self, weight: MeshFreeFunction) -> np.ndarray:
        r"""The weight integrated over each cell: ``[g, i]`` is :math:`\int_{V_i} w_g\,dV`, shape ``(G, N)``.

        The co-vector a flux integral :math:`\langle w, \phi\rangle` pairs with
        a cell-average flux, :math:`\sum_{g,i} \phi_{g,i} \int_{V_i} w_g\,dV`,
        which is exact for a piecewise-constant answer (#405 P2 step 7b.2.1).
        In frame vocabulary it is the analysis of the weight onto the mesh's
        indicator basis against the continuous measure
        :math:`dV = c\,dT(r)`, :math:`T(r) = r^p`. It is not the analysis
        against :attr:`volume_measure`, whose atoms sit at the cell centres:
        that would be the midpoint rule.

        * A :class:`~orpheus.numerics.mesh_free_function.RegionwiseConstant`
          is read through the region labels (:meth:`per_cell`; its rows are
          the regions, so its region count must be the mesh's) times the
          stored cell volumes. The stored volume is the cell's measure (it
          can sit a few ulp from :math:`c(T(r_{i+1}) - T(r_i))`, and on a thin
          cell far from the axis the two differ: `[M]` 1.9e-8 relative on a
          hollow cylinder at r = 100, qa of 7b.2.1).
        * A :class:`~orpheus.numerics.mesh_free_function.Symbolic` weight is
          integrated in :math:`r` by SymPy: one antiderivative per group
          (continuous for a ``Piecewise``, so a step inside a cell is
          integrated exactly), evaluated at the cells' edges as the binary
          rationals they are, then converted to a double and scaled by
          :math:`c` (two roundings, at most 2.7e-16 relative, `[M]`
          against mpmath). A weight that depends on the direction is refused
          (:meth:`~orpheus.numerics.mesh_free_function.Symbolic.without`).

        SCOPE-BOUNDARY[guard] machinery: a quadrature split at the roots of any step condition.
        ruling: the orchestrator, #405 P2 step 7b.2.1, spec §7b.2.1 (qa: `[M]` SymPy 1.14 integrates a sin step wrongly).
        revisit: when a weight with a transcendental step is needed.
        Qa's witness: ``Piecewise((1, sin(3r) > 0), (0, True))`` over [2, 3]
        integrates to 0. A non-smooth construct (a ``Piecewise`` condition, ``Heaviside``,
        ``Abs``, ``Max``, ``Min``, ``sign``) is admitted only with arguments
        polynomial in :math:`r`, where SymPy's antiderivative is reliable;
        any other failure of the integration is refused naming the weight.
        """
        match weight:
            case RegionwiseConstant():
                if weight.n_regions != len(self.region_materials):
                    raise ValueError(
                        f"Mesh1D.cell_integrals: the weight has {weight.n_regions} region(s); "
                        f"the mesh has {len(self.region_materials)}"
                    )
                return self.per_cell(weight.values).T * self.volumes
            case Symbolic():
                import sympy

                r = Symbolic.r
                measure = sympy.diff(r ** self.coord.measure_coordinate.exponent, r)
                edges = [sympy.Rational(float(e)) for e in self.edges]
                integrals = np.empty((weight.n_groups, self.N))
                for g, w in enumerate(weight.without(Symbolic.mu, Symbolic.phi).expressions):
                    _require_polynomial_steps(w, r, g)
                    try:
                        antiderivative = sympy.integrate(w * measure, r)
                        at_edges = [antiderivative.subs(r, e) for e in edges]
                        integrals[g] = [float(b - a) for a, b in zip(at_edges, at_edges[1:])]
                    except (TypeError, ValueError, NotImplementedError) as failure:
                        raise ValueError(
                            f"Mesh1D.cell_integrals: group {g}'s weight {w} cannot be integrated over the "
                            f"cells ({type(failure).__name__}: {failure})"
                        ) from None
                return self.coord.measure_constant * integrals
            case _:
                assert_never(weight)

    @property
    def total_width(self) -> float:
        """Total extent of the mesh (outer edge minus inner edge)."""
        return float(self.edges[-1] - self.edges[0])

    @property
    def boundary_points(self) -> tuple[float, ...]:
        """The positions of the boundary faces, inner first, in the order of :attr:`face_laws`."""
        return self.coord.boundary_points(float(self.edges[0]), float(self.edges[-1]))

    @property
    def outer_law(self) -> "BC | BoundaryTraceLaw":
        r"""The law on the outer face, :math:`r = r_R` (a slab's right face)."""
        return self.face_laws["xmax"]

    @property
    def volume_measure(self):
        r"""Cell-volume :class:`~orpheus.numerics.measure.DiscreteMeasure`.

        The natural integration measure :math:`\mu_V = \sum_i V_i \,
        \delta_{c_i}` whose atoms are the cell centres :math:`c_i =
        \tfrac12(x_{i-1/2} + x_{i+1/2})` carrying weight equal to the
        cell volume :math:`V_i`. Used at production integration
        sites where the historic ``np.sum(values * mesh.volumes)``
        idiom expresses an integral against the cell-volume measure
        — calling ``mesh.volume_measure(values)`` makes the integral
        structurally explicit and removes the manual ``np.sum`` /
        broadcasting boilerplate.

        See Also
        --------
        :meth:`orpheus.numerics.measure.DiscreteMeasure.integrate` —
        the array overload accepts the pre-evaluated value array
        directly (Issue 9.6 B4 extension).
        """
        # Local import to avoid a circular dependency at module
        # import time — :mod:`orpheus.numerics.measure` does not
        # import :mod:`orpheus.mesh.structured`, but the inverse
        # direction would force every consumer of mesh.py to bring
        # in the measure module.
        from orpheus.numerics.manifold import RealSpace
        from orpheus.numerics.measure import DiscreteMeasure
        return DiscreteMeasure(
            nodes=self.centers,
            weights=self.volumes,
            # The same ``RealSpace(1)`` ``indicator_basis`` partitions below.
            support=RealSpace(1),
        )

    def indicator_basis(self):
        r"""The mesh's cells AS a piecewise-constant :class:`~orpheus.numerics.basis.IndicatorBasis`.

        The **trial (synthesis) side** of a homogenisation / condensation
        :class:`~orpheus.numerics.frame.FrameBase` (the flux-weighted case is a
        :class:`~orpheus.numerics.frame.PetrovGalerkinFrame`) — the span of the
        cell indicators :math:`\{\mathbf{1}_R\}` is exactly the space of
        functions piecewise-constant on this mesh's cells, the coarse
        target of the flux-weighted projection (see
        :meth:`orpheus.sn.solution.Solution.homogenize`).

        Symmetric with :attr:`volume_measure`: the mesh **yields** its
        basis view, exactly as it yields its measure view. The mesh does
        **not** inherit :class:`~orpheus.numerics.basis.base.Basis` — a
        basis is the measure-*free* half of a frame, while the mesh
        carries the volume measure, so inheriting would conflate the two
        roles. The yielded :class:`IndicatorBasis` is geometry-free (it
        holds only this mesh's edge array), keeping
        :mod:`orpheus.numerics` free of any geometry dependency.
        """
        # Local import — same circular-dependency avoidance as
        # :attr:`volume_measure` (numerics does not import geometry).
        from orpheus.numerics.basis.indicator_basis import IndicatorBasis
        from orpheus.numerics.manifold import RealSpace
        return IndicatorBasis(
            edges_per_axis=(np.asarray(self.edges, dtype=float),),
            partition_of=RealSpace(1),
        )

    def with_distinct_cell_ids(self) -> "Mesh1D":
        r"""This mesh with every cell its **own** region and material id, ``0 .. N-1`` in cell order.

        The geometry the **homogenisation** result needs: a coarse
        :class:`~orpheus.transport.mesh.material_mesh.MaterialMesh` carries one
        fresh effective :class:`~orpheus.data.macro_xs.mixture.Mixture` per cell,
        so the cell-indexed homogenised materials key 1:1 into the cells.
        Polymorphic with :meth:`Mesh2D.with_distinct_cell_ids`, so
        ``Solution.homogenize`` stays dimension-agnostic. Edges, volumes and face
        laws carry through unchanged.
        """
        return replace(self, region_ids=np.arange(self.N, dtype=int), region_materials=tuple(range(self.N)))


# ═══════════════════════════════════════════════════════════════════════
# Mesh2D
# ═══════════════════════════════════════════════════════════════════════

@dataclass(frozen=True, eq=False)
class Mesh2D(ContentIdentity):
    """Two-dimensional mesh: Cartesian (x, y) or cylindrical (r, z).

    The boundary faces follow the one topology rule of
    :func:`~orpheus.mesh.face_laws.face_inventory`, as :class:`Mesh1D`'s do:
    along the first axis the :meth:`CoordSystem.boundary_points` of its edges
    (both ends on (x, y) and on a hollow (r, z); only the outer surface on a
    solid (r, z), whose axis :math:`r = 0` is interior and carries no law),
    and along the second axis (y or z) both ends, named ``xmin``, ``xmax``,
    ``ymin``, ``ymax``. A law per face, rather than per side, is the seed for
    a side whose faces carry different laws.

    Parameters
    ----------
    edges_x : ndarray, shape (Nx+1,)
        Edge positions in the first direction (x or r).
    edges_y : ndarray, shape (Ny+1,)
        Edge positions in the second direction (y or z).
    mat_map : ndarray, shape (Nx, Ny)
        Integer material ID for each cell.
    face_laws : mapping of face name to BC or BoundaryTraceLaw
        Keyword-only. One law per boundary face of the mesh's
        :func:`~orpheus.mesh.face_laws.face_inventory`, exactly: a missing
        face and an extra face are both refused, and ``None`` is not a law
        (:func:`~orpheus.geometry.structured_geometry.parse_boundary_law`).
        Stored as a :class:`~orpheus.mesh.face_laws.FaceLaws`, in inventory
        order, the value :class:`Mesh1D` stores too.
    coord : CoordSystem
        ``CARTESIAN`` for (x, y) or ``CYLINDRICAL`` for (r, z).

    Notes
    -----
    **A value with content identity**, as :class:`Mesh1D`
    (:class:`~orpheus.numerics.content.ContentIdentity`, #405 P1 step 5,
    2026-10-02): the arrays are stored as read-only copies, ``-0.0``
    canonicalised, and two meshes with equal edges, material map, laws and
    coordinate system are equal and hash alike. Until then the constructor
    aliased the caller's writeable arrays, and ``==`` raised on any pair of
    distinct meshes (the dataclass default compared ndarrays).
    """

    edges_x: np.ndarray
    edges_y: np.ndarray
    mat_map: np.ndarray
    _: KW_ONLY
    face_laws: "Mapping[str, BC | BoundaryTraceLaw]"
    coord: CoordSystem = CoordSystem.CARTESIAN

    def __post_init__(self) -> None:
        # The shared parsers, as Mesh1D: read-only copies (the mesh is a
        # value, and the caller's arrays stay theirs), ``-0.0`` canonicalised,
        # NaN, infinities, booleans and non-integral material ids refused.
        if np.ndim(self.edges_x) != 1 or np.ndim(self.edges_y) != 1:
            raise ValueError("edges_x and edges_y must be 1-D arrays")
        edges_x = parse_positions(self.edges_x, "Mesh2D.edges_x")
        edges_y = parse_positions(self.edges_y, "Mesh2D.edges_y")
        rows = parse_entries(self.mat_map, "Mesh2D.mat_map", "rows of int")
        mat_map = np.array(
            [[parse_integer(m, f"Mesh2D.mat_map[{i}][{j}]", "a material id")
              for j, m in enumerate(parse_entries(row, f"Mesh2D.mat_map[{i}]", "int"))]
             for i, row in enumerate(rows)],
            dtype=int,
        )
        if len(edges_x) < 2 or len(edges_y) < 2:
            raise ValueError("edge arrays must have at least 2 elements")
        if not np.all(np.diff(edges_x) > 0):
            raise ValueError("edges_x must be strictly monotonically increasing")
        if not np.all(np.diff(edges_y) > 0):
            raise ValueError("edges_y must be strictly monotonically increasing")

        nx = len(edges_x) - 1
        ny = len(edges_y) - 1
        if mat_map.shape != (nx, ny):
            raise ValueError(
                f"mat_map shape {mat_map.shape} must be ({nx}, {ny})"
            )
        if self.coord not in (CoordSystem.CARTESIAN, CoordSystem.CYLINDRICAL):
            raise ValueError(
                f"Mesh2D supports CARTESIAN or CYLINDRICAL, got {self.coord}"
            )

        face_laws = FaceLaws.over(
            face_inventory(self.coord, edges_x, 2), self.face_laws,
            "Mesh2D.face_laws", self.coord,
        )

        mat_map.flags.writeable = False
        object.__setattr__(self, "face_laws", face_laws)
        object.__setattr__(self, "edges_x", edges_x)
        object.__setattr__(self, "edges_y", edges_y)
        object.__setattr__(self, "mat_map", mat_map)

    # ── Derived properties ────────────────────────────────────────────

    @property
    def nx(self) -> int:
        """Number of cells in x (or r) direction."""
        return len(self.edges_x) - 1

    @property
    def ny(self) -> int:
        """Number of cells in y (or z) direction."""
        return len(self.edges_y) - 1

    @property
    def dx(self) -> np.ndarray:
        """Cell widths in x (or r) direction, shape (Nx,)."""
        return np.diff(self.edges_x)

    @property
    def dy(self) -> np.ndarray:
        """Cell widths in y (or z) direction, shape (Ny,)."""
        return np.diff(self.edges_y)

    @property
    def centers_x(self) -> np.ndarray:
        """Cell centres in x (or r) direction, shape (Nx,)."""
        return 0.5 * (self.edges_x[:-1] + self.edges_x[1:])

    @property
    def centers_y(self) -> np.ndarray:
        """Cell centres in y (or z) direction, shape (Ny,)."""
        return 0.5 * (self.edges_y[:-1] + self.edges_y[1:])

    @property
    def volumes(self) -> np.ndarray:
        """Cell volumes, shape (Nx, Ny).  Formula depends on *coord*."""
        return compute_volumes_2d(self.coord, self.edges_x, self.edges_y)

    @property
    def volume_measure(self):
        r"""2-D cell-volume :class:`~orpheus.numerics.measure.DiscreteMeasure`.

        Atoms are the ``(x_i, y_j)`` cell-centre pairs ordered with
        ``np.meshgrid(..., indexing='ij')`` — the same layout the
        ``volumes.ravel()`` exposes — and weights are the
        flattened cell volumes. Shape ``(Nx*Ny, 2)`` for the node
        array, ``(Nx*Ny,)`` for the weights.

        See Also
        --------
        :meth:`Mesh1D.volume_measure` — the 1-D analogue.
        """
        from orpheus.numerics.manifold import RealSpace
        from orpheus.numerics.measure import DiscreteMeasure
        cx = self.centers_x
        cy = self.centers_y
        X, Y = np.meshgrid(cx, cy, indexing="ij")
        nodes = np.stack([X.ravel(), Y.ravel()], axis=-1)  # (Nx*Ny, 2)
        weights = self.volumes.ravel()                     # (Nx*Ny,)
        return DiscreteMeasure(
            nodes=nodes,
            weights=weights,
            support=RealSpace(2),
        )

    @property
    def mat_ids(self) -> np.ndarray:
        """Flat material-ID array, shape (Nx*Ny,).

        Compatible with :func:`data.macro_xs.cell_xs.assemble_cell_xs`.
        """
        return self.mat_map.ravel()

    def indicator_basis(self):
        r"""The mesh's cells AS a 2-D :class:`~orpheus.numerics.basis.IndicatorBasis`.

        The 2-D analogue of :meth:`Mesh1D.indicator_basis` — the **trial (synthesis)
        side** of a homogenisation :class:`~orpheus.numerics.frame.FrameBase`. The
        indicator basis holds both axes' edges (``edges_per_axis = (edges_x,
        edges_y)``), so its membership table is built per axis and flattened in the
        ``"ij"`` / C order that matches this mesh's :attr:`volume_measure` nodes (and
        ``mat_map.ravel()``) — the same flat-cell ordering in any dimension.
        """
        from orpheus.numerics.basis.indicator_basis import IndicatorBasis
        from orpheus.numerics.manifold import RealSpace
        return IndicatorBasis(
            edges_per_axis=(
                np.asarray(self.edges_x, dtype=float),
                np.asarray(self.edges_y, dtype=float),
            ),
            partition_of=RealSpace(2),
        )

    def with_distinct_cell_ids(self) -> "Mesh2D":
        r"""A geometry-identical copy whose every cell is its **own** material id.

        Returns this mesh with ``mat_map = arange(Nx*Ny).reshape(Nx, Ny)`` — one
        distinct id per cell in ``"ij"`` / C order, matching the
        :meth:`indicator_basis` cell ordering — so the cell-indexed homogenised
        materials key 1:1 into the cells. Polymorphic with
        :meth:`Mesh1D.with_distinct_cell_ids`, so ``Solution.homogenize`` stays
        dimension-agnostic. The geometry (edges, coord, BCs) carries through.
        """
        mat_map = np.arange(self.nx * self.ny, dtype=int).reshape(self.nx, self.ny)
        return replace(self, mat_map=mat_map)

    def cell_integrals(self, weight: MeshFreeFunction) -> np.ndarray:
        """The 2-D twin of :meth:`Mesh1D.cell_integrals`; not built.

        ELEGANCE-DEBT[guard] #569: retires when ``Mesh2D`` carries region
        labels and ``Symbolic`` owns the second axis's coordinate (z, or x
        and y), so a weight can be read through the labels or integrated over
        a 2-D cell. The one refusal stands here, where the missing code goes,
        so no reader branches on the mesh's type.
        """
        raise NotImplementedError(
            "Mesh2D.cell_integrals: a 2-D cell integral needs the mesh's region labels and the second "
            "coordinate of a symbolic weight (#569)"
        )
