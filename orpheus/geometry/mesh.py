"""Mesh data structures for 1-D and 2-D geometries.

A mesh is an immutable description of a spatial domain: cell edges,
material assignments, and derived quantities (volumes, surfaces,
cell centres).  Solvers receive a mesh and build mutable,
solver-specific state on top of it.

Both :class:`Mesh1D` and :class:`Mesh2D` are frozen dataclasses --
once created, their fields cannot be reassigned.

Mesh construction from geometry
-------------------------------

The canonical path from a :class:`StructuredGeometry` to a
:class:`Mesh1D` is :meth:`Mesh1D.from_geometry`, which takes the
per-region discretization description as a tuple of
:class:`RegionMesh` instances. Discretization (cell counts, scheme)
is a mesh-layer concern — the geometry itself doesn't know or care
about cell counts.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Literal

import numpy as np

from .coord import (
    CoordSystem,
    compute_surfaces_1d,
    compute_volumes_1d,
    compute_volumes_2d,
)

if TYPE_CHECKING:
    from .structured_geometry import StructuredGeometry


# ═══════════════════════════════════════════════════════════════════════
# Boundary condition declaration
# ═══════════════════════════════════════════════════════════════════════

@dataclass(frozen=True)
class BC:
    """Solver-agnostic boundary condition declaration.

    A lightweight tag attached to geometry surfaces.  The geometry module
    makes no assumptions about what a given ``kind`` means — semantics
    are resolved by each solver's augmented mesh at construction time
    via its ``BC_REGISTRY``.

    Parameters
    ----------
    kind : str
        Boundary condition identifier (e.g. ``"vacuum"``,
        ``"reflective"``, ``"white"``).  Each solver defines which
        kinds it supports.
    params : dict[str, float]
        Optional numeric parameters (e.g. ``{"albedo": 0.7}``).
    """

    kind: str
    params: dict[str, float] = field(default_factory=dict)

    def __repr__(self) -> str:
        if self.params:
            return f"BC({self.kind!r}, {self.params!r})"
        return f"BC({self.kind!r})"

    def to_alpha(self) -> float:
        r"""Map this BC tag to a continuous specular albedo
        :math:`\alpha \in [0, 1]`.

        The trajectory_resolvent / Birkhoff–Sinai billiard family
        parametrises specular boundary conditions on a continuous
        albedo: :math:`\alpha = 0` is vacuum (no return), :math:`\alpha
        = 1` is perfect specular reflection, :math:`\alpha \in (0, 1)`
        is partial reflection. This method translates the production
        :class:`BC` tag-system to that scalar:

        * :data:`BC.vacuum` → ``0.0``
        * :data:`BC.reflective` → ``1.0``
        * ``BC("partial", {"albedo": x})`` → ``x``

        Other tags (``"white"``, Marshak diffuse, …) raise
        :class:`NotImplementedError` — they require a different
        closure structure that the specular-albedo parametrisation
        does not represent.

        Returns
        -------
        float
            The continuous-albedo equivalent of this BC tag.

        Raises
        ------
        ValueError
            If the BC kind is ``"partial"`` and the ``"albedo"`` key
            is missing from :attr:`params`.
        NotImplementedError
            If the BC kind has no specular-albedo equivalent
            (e.g. ``"white"``).
        """
        if self.kind == "vacuum":
            return 0.0
        if self.kind == "reflective":
            return 1.0
        if self.kind == "partial":
            try:
                return float(self.params["albedo"])
            except KeyError as exc:
                raise ValueError(
                    f"BC.to_alpha: BC('partial', ...) is missing the "
                    f"'albedo' parameter; got params={self.params!r}"
                ) from exc
        raise NotImplementedError(
            f"BC.to_alpha: BC kind {self.kind!r} has no specular-albedo "
            f"equivalent. Tags supported today: vacuum, reflective, partial."
        )


# Convenience instances — tab-completable, zero-import overhead.
BC.vacuum = BC("vacuum")  # type: ignore[attr-defined]
BC.reflective = BC("reflective")  # type: ignore[attr-defined]
BC.white = BC("white")  # type: ignore[attr-defined]


# ═══════════════════════════════════════════════════════════════════════
# RegionMesh — per-region discretization description (mesh-layer concern)
# ═══════════════════════════════════════════════════════════════════════

@dataclass(frozen=True)
class RegionMesh:
    r"""How to discretize one region of a :class:`StructuredGeometry`.

    Lives at the **mesh layer** because discretization is a mesh
    concern, not a geometry concern. The same
    :class:`StructuredGeometry` can be meshed with different
    :class:`RegionMesh` tuples for different studies (mesh refinement,
    uniform-vs-equal-volume comparison, future temperature-aware
    schemes).

    Parameters
    ----------
    n_cells : int
        Number of sub-cells inside this region. Must be ≥ 1.
    method : str
        Discretization scheme. Today's options:

        * ``"equal-volume"`` (default) — cells of equal volume.
          Yields uniform width for Cartesian; radially-graded width
          for cylindrical (``r ∝ √k``) and spherical (``r ∝ k^(1/3)``)
          per the equal-volume invariants in
          :func:`~orpheus.geometry.factories._subdivide_zone`.
        * ``"uniform"`` — cells of equal radial extent regardless of
          coordinate system. Useful when a method's accuracy depends
          on radial spacing rather than volume.

        Future: ``"temperature-graded"`` for inhomogeneous-temperature
        cases (Doppler-broadening grids), ``"adaptive"`` for output
        of refinement studies.

    Examples
    --------

    Default scheme (equal-volume), single region with 64 cells::

        RegionMesh(n_cells=64)

    Uniform discretization, fine outer mesh::

        (
            RegionMesh(n_cells=10),                          # equal-volume
            RegionMesh(n_cells=20, method="uniform"),        # uniform
        )

    See Also
    --------
    Mesh1D.from_geometry : Construct a Mesh1D from a
        :class:`StructuredGeometry` + tuple of :class:`RegionMesh`.
    Region : The geometry-layer per-region descriptor that this
        :class:`RegionMesh` pairs with at mesh-build time.
    """

    n_cells: int
    method: Literal["equal-volume", "uniform"] = "equal-volume"

    def __post_init__(self) -> None:
        if not isinstance(self.n_cells, (int, np.integer)):
            raise TypeError(
                f"RegionMesh.n_cells must be int, "
                f"got {type(self.n_cells).__name__}"
            )
        if self.n_cells < 1:
            raise ValueError(
                f"RegionMesh.n_cells must be ≥ 1; got {self.n_cells}"
            )
        if self.method not in ("equal-volume", "uniform"):
            raise ValueError(
                f"RegionMesh.method must be 'equal-volume' or 'uniform'; "
                f"got {self.method!r}"
            )


# ═══════════════════════════════════════════════════════════════════════
# Mesh1D
# ═══════════════════════════════════════════════════════════════════════

@dataclass(frozen=True)
class Mesh1D:
    """One-dimensional mesh in Cartesian, cylindrical, or spherical coordinates.

    Parameters
    ----------
    edges : ndarray, shape (N+1,)
        Monotonically increasing cell boundary positions.
        For cylindrical / spherical meshes these are radii.
    mat_ids : ndarray, shape (N,)
        Integer material ID for each cell.
    coord : CoordSystem
        Coordinate system (default: Cartesian).
    precomputed_volumes : ndarray or None, shape (N,)
        Optional override for cell volumes. When provided, the
        :attr:`volumes` property returns this array verbatim instead
        of recomputing it from ``edges``. This is the escape hatch for
        equal-volume subdivisions (cylindrical / spherical) where
        recomputation from edges loses ~1 ULP per cell through the
        ``sqrt→**2`` or ``cbrt→**3`` round trip and breaks invariants
        like "every cell in a zone has identical volume." Set by
        :func:`~orpheus.geometry.factories.mesh1d_from_zones` and
        friends; None for manually-constructed meshes with arbitrary
        edges, which still derive volumes from edges as before.
    """

    edges: np.ndarray
    mat_ids: np.ndarray
    coord: CoordSystem = CoordSystem.CARTESIAN
    precomputed_volumes: np.ndarray | None = None
    bc_left: BC | None = None
    bc_right: BC | None = None

    def __post_init__(self) -> None:
        edges = np.asarray(self.edges, dtype=float)
        mat_ids = np.asarray(self.mat_ids, dtype=int)

        if edges.ndim != 1 or mat_ids.ndim != 1:
            raise ValueError("edges and mat_ids must be 1-D arrays")
        if len(edges) < 2:
            raise ValueError("edges must have at least 2 elements (1 cell)")
        if len(mat_ids) != len(edges) - 1:
            raise ValueError(
                f"len(mat_ids)={len(mat_ids)} must equal "
                f"len(edges)-1={len(edges) - 1}"
            )
        if not np.all(np.diff(edges) > 0):
            raise ValueError("edges must be strictly monotonically increasing")

        precomputed = self.precomputed_volumes
        if precomputed is not None:
            precomputed = np.asarray(precomputed, dtype=float)
            if precomputed.shape != (len(edges) - 1,):
                raise ValueError(
                    f"precomputed_volumes shape {precomputed.shape} must be "
                    f"({len(edges) - 1},)"
                )
            if not np.all(precomputed > 0):
                raise ValueError("precomputed_volumes must be strictly positive")

        # Validate BC fields
        for attr in ("bc_left", "bc_right"):
            bc = getattr(self, attr)
            if bc is not None and not isinstance(bc, BC):
                raise TypeError(
                    f"{attr} must be a BC instance or None, got {type(bc).__name__}"
                )

        # Store validated arrays (frozen bypass via object.__setattr__)
        object.__setattr__(self, "edges", edges)
        object.__setattr__(self, "mat_ids", mat_ids)
        object.__setattr__(self, "precomputed_volumes", precomputed)

    # ── Derived properties ────────────────────────────────────────────

    @property
    def N(self) -> int:
        """Number of cells."""
        return len(self.edges) - 1

    @property
    def widths(self) -> np.ndarray:
        """Cell widths (edge-to-edge distance), shape (N,)."""
        return np.diff(self.edges)

    @property
    def centers(self) -> np.ndarray:
        """Cell centre positions, shape (N,)."""
        return 0.5 * (self.edges[:-1] + self.edges[1:])

    @property
    def volumes(self) -> np.ndarray:
        """Cell volumes, shape (N,).

        When ``precomputed_volumes`` was supplied at construction (the
        normal path via :func:`~orpheus.geometry.factories.mesh1d_from_zones`),
        those exact values are returned. Otherwise volumes are derived
        from edges via :func:`~orpheus.geometry.coord.compute_volumes_1d`
        — which is correct to ~1 ULP per cell but can break
        "all cells in an equal-volume zone are bit-identical"
        assertions for the cylindrical/spherical cases where the
        ``sqrt→**2`` / ``cbrt→**3`` edge round trip loses precision.
        """
        if self.precomputed_volumes is not None:
            return self.precomputed_volumes
        return compute_volumes_1d(self.coord, self.edges)

    @property
    def surfaces(self) -> np.ndarray:
        """Surface areas at each edge, shape (N+1,).  Formula depends on *coord*."""
        return compute_surfaces_1d(self.coord, self.edges)

    @property
    def total_width(self) -> float:
        """Total extent of the mesh (outer edge minus inner edge)."""
        return float(self.edges[-1] - self.edges[0])

    # ─────────────────────────────────────────────────────────────────
    # Construction from StructuredGeometry — the canonical entry point
    # ─────────────────────────────────────────────────────────────────

    @classmethod
    def from_geometry(
        cls,
        geometry: "StructuredGeometry",
        *,
        region_meshes: tuple[RegionMesh, ...],
        origin: float = 0.0,
    ) -> "Mesh1D":
        r"""Build a :class:`Mesh1D` from a :class:`StructuredGeometry`
        plus a per-region discretization description.

        This is the canonical geometry → mesh transition. Production
        solvers (CP, SN, MOC, MC) consume the resulting :class:`Mesh1D`
        directly; reference solvers do not need a mesh and consume the
        :class:`StructuredGeometry` instead.

        Each region in :attr:`geometry.regions <StructuredGeometry.regions>`
        is paired with the matching :class:`RegionMesh` at the same
        index. The lengths must match. For each pair:

        * The region's ``mat_id`` is broadcast across the resulting
          cells.
        * The region's ``outer_thickness_cm`` and the region-mesh's
          ``n_cells`` and ``method`` determine the cell edges and
          volumes.

        Edges are accumulated across regions starting from
        :paramref:`origin` (default 0). The first region starts at
        ``origin``; the last region ends at ``origin +
        geometry.domain_extent_cm``.

        BCs from :attr:`geometry.bcs <StructuredGeometry.bcs>` are
        propagated onto the resulting mesh's :attr:`bc_left` /
        :attr:`bc_right` fields:

        * ``"SLB"``: ``bcs[0] → bc_left``, ``bcs[1] → bc_right``.
        * ``"CYL"`` / ``"SPH"``: ``bcs[0] → bc_right`` (outer surface);
          ``bc_left`` is left ``None`` (the centreline is implicit
          reflective at the coordinate origin and is interpreted by
          each solver's augmented mesh).

        Parameters
        ----------
        geometry : StructuredGeometry
            The geometry to discretize.
        region_meshes : tuple[RegionMesh, ...]
            Per-region discretization descriptors. Length must equal
            ``len(geometry.regions)``.
        origin : float, optional
            Position of the inner-most edge in cm. Default 0.0.

        Returns
        -------
        Mesh1D
            A frozen mesh with cells, mat_ids, exact precomputed
            volumes, and BCs propagated from the geometry.

        Raises
        ------
        ValueError
            If ``len(region_meshes) != len(geometry.regions)``.

        Examples
        --------

        Three-region pin cell::

            geom = StructuredGeometry.wigner_seitz_pin_cell(
                r_fuel=0.9, r_clad=1.1, pitch=3.6,
            )
            mesh = Mesh1D.from_geometry(geom, region_meshes=(
                RegionMesh(n_cells=10),
                RegionMesh(n_cells=3),
                RegionMesh(n_cells=7),
            ))
        """
        # Local imports to avoid circular: structured_geometry imports
        # from this module (BC), and we import StructuredGeometry only
        # for type checking above.
        from .factories import _subdivide_zone

        if len(region_meshes) != len(geometry.regions):
            raise ValueError(
                f"Mesh1D.from_geometry: len(region_meshes)="
                f"{len(region_meshes)} must equal "
                f"len(geometry.regions)={len(geometry.regions)}."
            )

        coord = geometry.coord
        edges_list: list[np.ndarray] = []
        mat_ids_list: list[np.ndarray] = []
        volumes_list: list[np.ndarray] = []
        inner = float(origin)

        for region, region_mesh in zip(
            geometry.regions, region_meshes, strict=True,
        ):
            outer = inner + float(region.outer_thickness_cm)

            if region_mesh.method == "equal-volume":
                sub_edges, sub_volumes = _subdivide_zone(
                    inner, outer, region_mesh.n_cells, coord,
                )
            else:  # "uniform"
                sub_edges = np.linspace(
                    inner, outer, region_mesh.n_cells + 1,
                )
                sub_volumes = compute_volumes_1d(coord, sub_edges)

            # Skip the first sub-edge (== previous region's outer edge)
            # to avoid duplication at the inter-region boundary.
            edges_list.append(sub_edges[1:])
            mat_ids_list.append(
                np.full(region_mesh.n_cells, region.mat_id, dtype=int)
            )
            volumes_list.append(sub_volumes)
            inner = outer

        edges = np.concatenate([[float(origin)], *edges_list])
        mat_ids = np.concatenate(mat_ids_list)
        volumes = np.concatenate(volumes_list)

        # Map geometry BCs onto the mesh's bc_left / bc_right fields.
        if geometry.geometry == "SLB":
            bc_left, bc_right = geometry.bcs
        else:  # CYL / SPH — single outer BC; centreline implicit
            bc_left = None
            bc_right = geometry.bcs[0]

        return cls(
            edges=edges,
            mat_ids=mat_ids,
            coord=coord,
            precomputed_volumes=volumes,
            bc_left=bc_left,
            bc_right=bc_right,
        )


# ═══════════════════════════════════════════════════════════════════════
# Mesh2D
# ═══════════════════════════════════════════════════════════════════════

@dataclass(frozen=True)
class Mesh2D:
    """Two-dimensional mesh: Cartesian (x, y) or cylindrical (r, z).

    Parameters
    ----------
    edges_x : ndarray, shape (Nx+1,)
        Edge positions in the first direction (x or r).
    edges_y : ndarray, shape (Ny+1,)
        Edge positions in the second direction (y or z).
    mat_map : ndarray, shape (Nx, Ny)
        Integer material ID for each cell.
    coord : CoordSystem
        ``CARTESIAN`` for (x, y) or ``CYLINDRICAL`` for (r, z).
    """

    edges_x: np.ndarray
    edges_y: np.ndarray
    mat_map: np.ndarray
    coord: CoordSystem = CoordSystem.CARTESIAN
    bc_xmin: BC | None = None
    bc_xmax: BC | None = None
    bc_ymin: BC | None = None
    bc_ymax: BC | None = None

    def __post_init__(self) -> None:
        edges_x = np.asarray(self.edges_x, dtype=float)
        edges_y = np.asarray(self.edges_y, dtype=float)
        mat_map = np.asarray(self.mat_map, dtype=int)

        if edges_x.ndim != 1 or edges_y.ndim != 1:
            raise ValueError("edges_x and edges_y must be 1-D arrays")
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

        # Validate BC fields
        for attr in ("bc_xmin", "bc_xmax", "bc_ymin", "bc_ymax"):
            bc = getattr(self, attr)
            if bc is not None and not isinstance(bc, BC):
                raise TypeError(
                    f"{attr} must be a BC instance or None, got {type(bc).__name__}"
                )

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
    def mat_ids(self) -> np.ndarray:
        """Flat material-ID array, shape (Nx*Ny,).

        Compatible with :func:`data.macro_xs.cell_xs.assemble_cell_xs`.
        """
        return self.mat_map.ravel()
