r"""Dimension-agnostic axis primitive for SN phase-space construction.

The S\ :sub:`N` discretisation factors as a tensor product of per-axis
1-D meshes (grand report §15.1 ``L = Σ_axis (D_axis ⊗ Ω_axis ⊗ I_g)``).
This module declares the per-axis primitive used to build that product:

* :class:`Axis1D` — the protocol any 1-D mesh axis satisfies.
* :class:`AxisMesh` — Cartesian axis. Two endpoints (``min`` / ``max``),
  both BC-bearing.
* :class:`RadialAxisMesh` — solid radial axis (sphere / cylinder). One
  endpoint (``outer``); the pole at :math:`r=0` is a coordinate
  singularity treated by :class:`MorelMontryAngularSweep` (Hébert
  §3.9.4 Carlson coupled-pole inward sweep) — NOT a BC trace law.
  Conflating the pole with a face endpoint would force every
  consumer to handle a "BC that is not a BC", violating Pattern 4
  (illegal states unrepresentable).

The dim-agnostic shape primitives are exposed as **pure functions on
axis tuples** — :func:`spatial_shape`, :func:`face_labels`,
:func:`face_shape`, :func:`face_outflow_ordinates`,
:func:`n_unknowns_flat`. The :class:`~orpheus.sn.geometry.SNMesh` class
delegates to these so 3-D admission gates can exercise the shape
algebra on a synthetic axis tuple without instantiating a full SN
phase space (no ``Mesh3D`` exists today — D9 of
``r1_phase_a_dim_agnostic_ultraplan.md``).

See ``docs/theory/sn_dim_agnostic.rst`` (lands in C8) for the
architectural narrative; this module is its implementation.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum
from typing import Protocol, runtime_checkable

import numpy as np

from orpheus.geometry import BC


# ═══════════════════════════════════════════════════════════════════════
# Coordinate-system enum for an axis
# ═══════════════════════════════════════════════════════════════════════

class AxisCoord(StrEnum):
    """Coordinate system carried by a single 1-D axis.

    Distinct from :class:`orpheus.geometry.coord.CoordSystem` (which
    describes a *whole-mesh* coordinate system) because dimension-
    agnostic construction needs to mix coordinate systems on different
    axes — e.g. an (r, z) 2-D mesh has a ``RADIAL_CYLINDRICAL`` axis
    and a ``CARTESIAN`` axis.
    """

    CARTESIAN = "cartesian"
    RADIAL_SPHERICAL = "radial_spherical"
    RADIAL_CYLINDRICAL = "radial_cylindrical"


# ═══════════════════════════════════════════════════════════════════════
# FaceLabel — canonical key for face buffers, BCs, ordinate masks
# ═══════════════════════════════════════════════════════════════════════

#: Spatial axis names, positional-by-axis — the same axis order as
#: :attr:`SNMesh.axes`, ``OctantLabel.signs``, and the per-axis kernel
#: tuples. Re-exported from its C5.3 home
#: :data:`orpheus.numerics.face_layout.AXIS_NAMES` (moved down so the
#: geometry-blind trace space shares the crosswalk without an sn-ward
#: import); SN consumers keep importing it from here, next to the axis
#: primitives it names.
from orpheus.numerics.face_layout import AXIS_NAMES  # noqa: E402

#: Canonical ``FaceLabel.endpoint`` → face-name suffix crosswalk. A
#: solid radial axis's single ``"outer"`` endpoint renders as the
#: ``max`` face of its axis (the historical curvilinear convention:
#: the outer radius IS ``xmax``). Any other endpoint label — e.g. an
#: :class:`AxisMesh` with overridden ``label_low`` / ``label_high`` —
#: has NO face name and fails loud in :attr:`FaceLabel.face_name`,
#: rather than silently desynchronizing from the ``"{axis}{min|max}"``
#: world that :class:`~orpheus.numerics.face_layout.FaceLayout`,
#: the trace space, and the sweep schedule all key on.
_ENDPOINT_SUFFIX = {"min": "min", "max": "max", "outer": "max"}


@dataclass(frozen=True, slots=True)
class FaceLabel:
    r"""Canonical key identifying one boundary face of a D-dim SN mesh.

    A face is identified by *which axis it lies on* (``axis_index``)
    and *which endpoint of that axis* (``endpoint`` — ``"min"`` /
    ``"max"`` for Cartesian, ``"outer"`` for solid radial).

    This dataclass is the load-bearing key for every dim-agnostic
    boundary-keyed lookup in the SN module: ``BoundaryFlux.face_buffers``
    (C4), ``SNMesh.bc`` (C4), the outflow-ordinate mask cache (C1),
    and the sweep DAG's face-trace state (C5). It exists so that all
    those consumers share ONE definition of "which face" rather than
    each rolling its own ad-hoc identifier (Pattern 2 — single source
    of truth).

    Parameters
    ----------
    axis_index : int
        Position of the axis in :attr:`SNMesh.axes`.
    endpoint : str
        Endpoint label on that axis (``"min"``, ``"max"``, ``"outer"``).
    """

    axis_index: int
    endpoint: str

    def __str__(self) -> str:
        return f"face_{self.axis_index}_{self.endpoint}"

    @property
    def face_name(self) -> str:
        r"""Cross-module face-name rendering, ``"{axis}{min|max}"``.

        The single-sourced crosswalk from the structural identity
        ``(axis_index, endpoint)`` to the string key that
        :class:`~orpheus.numerics.face_layout.FaceLayout`, the trace
        space, :attr:`SNMesh.bc`, and the sweep schedule all share:
        axis name from :data:`AXIS_NAMES`, endpoint suffix from
        :data:`_ENDPOINT_SUFFIX` (``"outer"`` renders as ``max`` — a
        solid radial axis's outer surface IS its ``max`` face).

        Raises
        ------
        ValueError
            If ``endpoint`` is not one of the canonical labels
            (``"min"`` / ``"max"`` / ``"outer"``) — an overridden
            axis endpoint label has no face name and must fail loud
            rather than silently desynchronize the face-name world.
        IndexError
            If ``axis_index`` exceeds the named-axis inventory.
        """
        suffix = _ENDPOINT_SUFFIX.get(self.endpoint)
        if suffix is None:
            canonical = ", ".join(f"'{e}'" for e in _ENDPOINT_SUFFIX)
            raise ValueError(
                f"FaceLabel endpoint '{self.endpoint}' has no face name; "
                f"canonical endpoints are {canonical}."
            )
        return f"{AXIS_NAMES[self.axis_index]}{suffix}"


# ═══════════════════════════════════════════════════════════════════════
# Axis1D Protocol
# ═══════════════════════════════════════════════════════════════════════

@runtime_checkable
class Axis1D(Protocol):
    r"""Capability contract for a 1-D mesh axis.

    Implementations declare:

    * ``edges`` — strictly monotonically increasing cell-face positions,
      shape ``(n+1,)``.
    * ``coord`` — the :class:`AxisCoord` for this axis.
    * ``n`` — cell count along the axis (= ``len(edges) - 1``).
    * ``endpoints`` — tuple of endpoint labels that carry BCs. Cartesian
      axes have two (``"min"``, ``"max"``); solid radial axes have one
      (``"outer"``); the pole at :math:`r=0` is intentionally NOT an
      endpoint here — it carries an angular-closure regularity
      condition, not a BC trace law.
    * ``bc`` — mapping ``endpoint_label → BC | None``.

    Pole rationale (D3 of the ultraplan): a BC trace law is a linear
    operator ``bc.apply(outflow) → inflow`` that closes the transport
    equation at a Dirichlet / Neumann / reflective boundary. The pole
    at :math:`r=0` has a fundamentally different mathematical structure
    — the angular redistribution coefficient :math:`1-\mu^2` vanishes
    at :math:`\mu=\pm 1`, so the streaming-collision balance decouples
    from the :math:`\alpha`-cascade and reduces to a plain DD inward
    recurrence in radius (Hébert Eq. 3.434). There is no "BC trace
    law" to apply; the seed comes from a separate sweep against a
    moment-folded source at :math:`\mu=-1`. Conflating the two into a
    single ``face_buffers[label]`` interface would force every
    consumer to handle a "BC that is not a BC", violating Pattern 4.
    """

    # Read-only properties (not mutable attributes): implementers are frozen
    # dataclasses / expose ``coord`` as a ``@property``. A mutable Protocol
    # attribute is invariant and rejects a frozen field / read-only property
    # ("edges is not read-only in protocol"; "coord property not assignable");
    # a read-only Protocol member accepts both. Matches ``n``/``endpoints``/
    # ``bc`` below (already read-only) — this is a pure read contract.
    @property
    def edges(self) -> np.ndarray: ...

    @property
    def coord(self) -> AxisCoord: ...

    @property
    def n(self) -> int: ...

    @property
    def endpoints(self) -> tuple[str, ...]: ...

    @property
    def bc(self) -> dict[str, BC | None]: ...


# ═══════════════════════════════════════════════════════════════════════
# Concrete axis implementations
# ═══════════════════════════════════════════════════════════════════════

@dataclass(frozen=True, slots=True)
class AxisMesh:
    r"""Cartesian 1-D axis. Two endpoints (``min``, ``max``).

    ``edges`` is the source of truth — ``n`` is derived as
    ``len(edges) - 1``, so a mismatched cell count is unrepresentable
    (Pattern 4). The label strings (``min`` / ``max``) are kept as
    fields so a user can override them (e.g. ``"left"`` / ``"right"``
    for slab-conventional naming) without rebuilding the dataclass.

    Parameters
    ----------
    edges : ndarray, shape (n+1,)
        Strictly monotonically increasing cell-face positions.
    bc_low : BC or None
        Boundary condition at the low endpoint. ``None`` means "use
        the mesh-level default" (typically reflective).
    bc_high : BC or None
        Boundary condition at the high endpoint.
    label_low, label_high : str
        Endpoint labels; default ``"min"`` / ``"max"``.
    """

    edges: np.ndarray
    bc_low: BC | None = None
    bc_high: BC | None = None
    label_low: str = "min"
    label_high: str = "max"

    def __post_init__(self) -> None:
        edges = np.ascontiguousarray(self.edges, dtype=float)
        if edges.ndim != 1 or edges.size < 2:
            raise ValueError(
                f"AxisMesh.edges must be a 1-D array with ≥2 entries; "
                f"got shape {edges.shape}"
            )
        if not np.all(np.diff(edges) > 0):
            raise ValueError(
                "AxisMesh.edges must be strictly monotonically increasing"
            )
        object.__setattr__(self, "edges", edges)

    @property
    def coord(self) -> AxisCoord:
        return AxisCoord.CARTESIAN

    @property
    def n(self) -> int:
        return self.edges.size - 1

    @property
    def endpoints(self) -> tuple[str, ...]:
        return (self.label_low, self.label_high)

    @property
    def bc(self) -> dict[str, BC | None]:
        return {self.label_low: self.bc_low, self.label_high: self.bc_high}


@dataclass(frozen=True, slots=True)
class RadialAxisMesh:
    r"""Solid radial 1-D axis (sphere or cylinder). One endpoint (``outer``).

    The pole at :math:`r=0` is **not** an endpoint in this abstraction
    — see the :class:`Axis1D` docstring for the pole-treatment
    rationale.

    Parameters
    ----------
    edges : ndarray, shape (n+1,)
        Strictly monotonically increasing radial positions. The first
        entry MAY be 0 (solid-radial mesh that touches the pole); when
        it is, the pole carries an angular-closure regularity
        condition rather than a BC.
    coord : AxisCoord
        Must be :attr:`AxisCoord.RADIAL_SPHERICAL` or
        :attr:`AxisCoord.RADIAL_CYLINDRICAL`.
    bc_outer : BC or None
        Boundary condition at the outer radial surface.
    label_outer : str
        Endpoint label; default ``"outer"``.
    """

    edges: np.ndarray
    coord: AxisCoord
    bc_outer: BC | None = None
    label_outer: str = "outer"

    def __post_init__(self) -> None:
        edges = np.ascontiguousarray(self.edges, dtype=float)
        if edges.ndim != 1 or edges.size < 2:
            raise ValueError(
                f"RadialAxisMesh.edges must be a 1-D array with ≥2 entries; "
                f"got shape {edges.shape}"
            )
        if not np.all(np.diff(edges) > 0):
            raise ValueError(
                "RadialAxisMesh.edges must be strictly monotonically increasing"
            )
        if edges[0] < 0:
            raise ValueError(
                f"RadialAxisMesh.edges[0] must be ≥0; got {edges[0]}"
            )
        if self.coord not in (
            AxisCoord.RADIAL_SPHERICAL, AxisCoord.RADIAL_CYLINDRICAL,
        ):
            raise ValueError(
                f"RadialAxisMesh.coord must be RADIAL_SPHERICAL or "
                f"RADIAL_CYLINDRICAL; got {self.coord!r}"
            )
        object.__setattr__(self, "edges", edges)

    @property
    def n(self) -> int:
        return self.edges.size - 1

    @property
    def endpoints(self) -> tuple[str, ...]:
        return (self.label_outer,)

    @property
    def bc(self) -> dict[str, BC | None]:
        return {self.label_outer: self.bc_outer}


# ═══════════════════════════════════════════════════════════════════════
# Dim-agnostic shape primitives — pure functions on axis tuples
# ═══════════════════════════════════════════════════════════════════════

def spatial_shape(axes: tuple[Axis1D, ...]) -> tuple[int, ...]:
    """Per-axis cell counts ``(n_0, n_1, ...)``."""
    return tuple(a.n for a in axes)


def face_labels(axes: tuple[Axis1D, ...]) -> tuple[FaceLabel, ...]:
    """Canonical face inventory derived from each axis's endpoints.

    The iteration order is ``(axis_index ascending, endpoint in
    axis.endpoints order)``. This order is the canonical concatenation
    order for :meth:`AngularFlux.to_flat` (C3) and the canonical
    iteration order for :attr:`BoundaryFlux.face_buffers` (C4).
    """
    return tuple(
        FaceLabel(axis_index=i, endpoint=ep)
        for i, axis in enumerate(axes)
        for ep in axis.endpoints
    )


def face_shape(
    axes: tuple[Axis1D, ...], label: FaceLabel,
) -> tuple[int, ...]:
    """Spatial shape of one boundary face.

    The face on ``axes[label.axis_index]`` lies in the codimension-1
    hyperplane spanned by the other axes; its shape is the per-axis
    cell count of those axes in axis-index order.
    """
    return tuple(
        axes[j].n for j in range(len(axes)) if j != label.axis_index
    )


_OUTWARD_ENDPOINTS = frozenset({"max", "outer"})
_OUTFLOW_COSINE_TOL = 1e-15


def face_outflow_ordinates(
    axes: tuple[Axis1D, ...], label: FaceLabel, quad,
) -> np.ndarray:
    r"""Ordinate indices whose direction-cosine is OUTWARD at this face.

    At the ``max`` / ``outer`` endpoint of an axis, ordinates with
    :math:`\mu_{axis} > 0` are outflowing; at ``min``, outflow is
    :math:`\mu_{axis} < 0`. The strictly-positive cutoff (``> 1e-15``)
    excludes ordinates exactly tangent to the face — these contribute
    neither inflow nor outflow at the boundary and are skipped by
    every consumer (the sweep, the matvec, and the pack convention).
    """
    mu_axis = quad.axis_cosines(label.axis_index)
    sign = +1.0 if label.endpoint in _OUTWARD_ENDPOINTS else -1.0
    return np.where(sign * mu_axis > _OUTFLOW_COSINE_TOL)[0]


def n_unknowns_flat(
    axes: tuple[Axis1D, ...], quad, ng: int,
) -> int:
    r"""Total flat-vector size for typed :class:`AngularFlux` on this mesh.

    The pack convention (C3) is the direct-sum decomposition
    :math:`V = V_\text{cells} \oplus \bigoplus_\ell V_{\text{face}, \ell}`
    where the cells block has size
    :math:`N \cdot n_g \cdot \prod_i n_i` and each face block has size
    :math:`|\text{outflow}_\ell| \cdot n_g \cdot \prod_{i \ne \text{axis}(\ell)} n_i`.

    Parameters
    ----------
    axes : tuple of Axis1D
        The axis tuple (i.e. :attr:`SNMesh.axes`).
    quad : Quadrature
        Angular quadrature; needed for the per-face outflow mask.
    ng : int
        Number of energy groups.
    """
    shape = spatial_shape(axes)
    n_cells = int(quad.N) * int(ng) * int(np.prod(shape, dtype=np.int64))
    n_face = 0
    for label in face_labels(axes):
        outflow = face_outflow_ordinates(axes, label, quad)
        block = (
            int(outflow.size)
            * int(ng)
            * int(np.prod(face_shape(axes, label), dtype=np.int64))
        )
        n_face += block
    return n_cells + n_face


def coord_system(axes: tuple[Axis1D, ...]):
    r"""Whole-mesh :class:`~orpheus.geometry.CoordSystem` of an axis tuple.

    The single-axis coordinate map is

    * :attr:`AxisCoord.CARTESIAN` → ``CoordSystem.CARTESIAN``
    * :attr:`AxisCoord.RADIAL_SPHERICAL` → ``CoordSystem.SPHERICAL``
    * :attr:`AxisCoord.RADIAL_CYLINDRICAL` → ``CoordSystem.CYLINDRICAL``

    Multi-axis tuples must be all-Cartesian — a curvilinear axis only
    has meaning as the sole axis of a 1-D solid mesh (the reduced
    streaming operators and the pole angular closure are 1-D
    constructions); mixing it into a product mesh is unrepresentable.

    Raises
    ------
    NotImplementedError
        If a multi-axis tuple contains a non-Cartesian axis.
    """
    from orpheus.geometry import CoordSystem

    if len(axes) == 1:
        return {
            AxisCoord.CARTESIAN: CoordSystem.CARTESIAN,
            AxisCoord.RADIAL_SPHERICAL: CoordSystem.SPHERICAL,
            AxisCoord.RADIAL_CYLINDRICAL: CoordSystem.CYLINDRICAL,
        }[axes[0].coord]
    non_cartesian = [ax.coord for ax in axes if ax.coord is not AxisCoord.CARTESIAN]
    if non_cartesian:
        raise NotImplementedError(
            f"coord_system: multi-axis meshes must be all-Cartesian; "
            f"got non-Cartesian axis coords {non_cartesian!r}"
        )
    return CoordSystem.CARTESIAN


# ═══════════════════════════════════════════════════════════════════════
# Legacy-mesh → axis-tuple adapter
# ═══════════════════════════════════════════════════════════════════════

def axes_from_legacy_mesh(mesh) -> tuple[Axis1D, ...]:
    r"""Extract the canonical axis tuple from a legacy ``Mesh1D`` /
    ``Mesh2D``.

    Wraps the legacy mesh's per-axis edge arrays + BC declarations in
    :class:`AxisMesh` / :class:`RadialAxisMesh` so the
    :class:`~orpheus.sn.geometry.SNMesh` constructor can expose
    :attr:`SNMesh.axes` regardless of which constructor surface the
    caller used.

    Mesh1D mapping (per coordinate system):

    * ``CARTESIAN`` → ``(AxisMesh(edges=mesh.edges, bc_low=mesh.bc_left,
      bc_high=mesh.bc_right),)``.
    * ``SPHERICAL`` → ``(RadialAxisMesh(edges=mesh.edges,
      coord=RADIAL_SPHERICAL, bc_outer=mesh.bc_right),)``. The
      centreline ``mesh.bc_left`` is implicit-reflective per
      :meth:`Mesh1D.from_geometry`; it does NOT become a BC on the
      axis because the pole is a coordinate singularity, not an
      endpoint.
    * ``CYLINDRICAL`` → same as SPHERICAL with
      ``RADIAL_CYLINDRICAL``.

    Mesh2D mapping (Cartesian only — :class:`Mesh2D` cylindrical is
    not used by any SN sweep today):

    * ``CARTESIAN`` → ``(AxisMesh(edges=mesh.edges_x,
      bc_low=mesh.bc_xmin, bc_high=mesh.bc_xmax),
      AxisMesh(edges=mesh.edges_y, bc_low=mesh.bc_ymin,
      bc_high=mesh.bc_ymax))``.
    """
    # Local import to avoid a circular dependency at module import time
    # (mesh.py is imported by sn/geometry.py which imports this module).
    from orpheus.geometry import CoordSystem, Mesh1D, Mesh2D

    if isinstance(mesh, Mesh1D):
        if mesh.coord == CoordSystem.CARTESIAN:
            return (
                AxisMesh(
                    edges=mesh.edges,
                    bc_low=mesh.bc_left,
                    bc_high=mesh.bc_right,
                ),
            )
        if mesh.coord == CoordSystem.SPHERICAL:
            return (
                RadialAxisMesh(
                    edges=mesh.edges,
                    coord=AxisCoord.RADIAL_SPHERICAL,
                    bc_outer=mesh.bc_right,
                ),
            )
        if mesh.coord == CoordSystem.CYLINDRICAL:
            return (
                RadialAxisMesh(
                    edges=mesh.edges,
                    coord=AxisCoord.RADIAL_CYLINDRICAL,
                    bc_outer=mesh.bc_right,
                ),
            )
        raise ValueError(
            f"axes_from_legacy_mesh: unsupported Mesh1D.coord {mesh.coord!r}"
        )

    if isinstance(mesh, Mesh2D):
        if mesh.coord != CoordSystem.CARTESIAN:
            raise NotImplementedError(
                f"axes_from_legacy_mesh: 2-D non-Cartesian "
                f"(coord={mesh.coord!r}) has no SN sweep today; "
                f"axis-tuple extraction is undefined."
            )
        return (
            AxisMesh(
                edges=mesh.edges_x,
                bc_low=mesh.bc_xmin,
                bc_high=mesh.bc_xmax,
            ),
            AxisMesh(
                edges=mesh.edges_y,
                bc_low=mesh.bc_ymin,
                bc_high=mesh.bc_ymax,
            ),
        )

    raise TypeError(
        f"axes_from_legacy_mesh: expected Mesh1D or Mesh2D; "
        f"got {type(mesh).__name__}"
    )


# ═══════════════════════════════════════════════════════════════════════
# Axis-tuple → legacy-mesh ADAPTER builder (for SNMesh.from_axes)
# ═══════════════════════════════════════════════════════════════════════

def legacy_mesh_from_axes(
    axes: tuple[Axis1D, ...],
    mat_map: np.ndarray | None = None,
):
    r"""Build the legacy :class:`Mesh1D` / :class:`Mesh2D` ADAPTER for an axis tuple.

    C5.1 (#225): this is NO LONGER a round-trip source — the axes an
    :class:`~orpheus.sn.geometry.SNMesh` carries are the caller's
    objects verbatim. :meth:`SNMesh.from_axes` calls this builder only
    to synthesize the d≤2 ``SNMesh.mesh`` adapter for the consumers
    still reading through it (1-D reduced streaming construction,
    trace build, realizer metadata) — each dissolves across C5.2–C5.5,
    and this builder narrows/retires with them.

    Parameters
    ----------
    axes : tuple of Axis1D
        Axes describing the mesh. Length 1 → :class:`Mesh1D`;
        length 2 → :class:`Mesh2D` Cartesian. Length ≥3 is not
        supported in C1 (no ``Mesh3D`` dataclass exists today —
        D9 of the ultraplan).
    mat_map : ndarray or None
        Material-id assignment. Shape ``(n_0,)`` for 1-D or
        ``(n_0, n_1)`` for 2-D. Defaults to all-zeros (single
        material with id 0).
    """
    from orpheus.geometry import CoordSystem, Mesh1D, Mesh2D

    shape = spatial_shape(axes)
    if mat_map is None:
        mat_map = np.zeros(shape, dtype=int)
    else:
        mat_map = np.asarray(mat_map, dtype=int)
        # NOTE: SNMesh._init_core re-validates this shape — it is the
        # SURVIVING guard when this adapter builder retires with the
        # legacy mesh consumers (C5.2-C5.5); do not delete the
        # _init_core copy as "redundant" with this one.
        if mat_map.shape != shape:
            raise ValueError(
                f"legacy_mesh_from_axes: mat_map shape {mat_map.shape} "
                f"must match spatial_shape={shape}"
            )

    if len(axes) == 1:
        (ax,) = axes
        if ax.coord == AxisCoord.CARTESIAN:
            assert isinstance(ax, AxisMesh)
            return Mesh1D(
                edges=ax.edges,
                mat_ids=mat_map.ravel(),
                coord=CoordSystem.CARTESIAN,
                bc_left=ax.bc_low,
                bc_right=ax.bc_high,
            )
        if ax.coord == AxisCoord.RADIAL_SPHERICAL:
            assert isinstance(ax, RadialAxisMesh)
            return Mesh1D(
                edges=ax.edges,
                mat_ids=mat_map.ravel(),
                coord=CoordSystem.SPHERICAL,
                bc_left=None,
                bc_right=ax.bc_outer,
            )
        if ax.coord == AxisCoord.RADIAL_CYLINDRICAL:
            assert isinstance(ax, RadialAxisMesh)
            return Mesh1D(
                edges=ax.edges,
                mat_ids=mat_map.ravel(),
                coord=CoordSystem.CYLINDRICAL,
                bc_left=None,
                bc_right=ax.bc_outer,
            )
        raise ValueError(
            f"legacy_mesh_from_axes: unsupported AxisCoord {ax.coord!r}"
        )

    if len(axes) == 2:
        ax0, ax1 = axes
        if ax0.coord != AxisCoord.CARTESIAN or ax1.coord != AxisCoord.CARTESIAN:
            raise NotImplementedError(
                f"legacy_mesh_from_axes: 2-D non-Cartesian (axes coords "
                f"{ax0.coord!r}, {ax1.coord!r}) has no Mesh2D today."
            )
        assert isinstance(ax0, AxisMesh) and isinstance(ax1, AxisMesh)
        return Mesh2D(
            edges_x=ax0.edges,
            edges_y=ax1.edges,
            mat_map=mat_map,
            coord=CoordSystem.CARTESIAN,
            bc_xmin=ax0.bc_low,
            bc_xmax=ax0.bc_high,
            bc_ymin=ax1.bc_low,
            bc_ymax=ax1.bc_high,
        )

    raise NotImplementedError(
        f"legacy_mesh_from_axes: the legacy mesh ADAPTER is genuinely "
        f"d≤2 (Mesh1D / Mesh2D are the d≤2 user-facing dataclasses); "
        f"{len(axes)}-axis meshes are mesh-adapter-free by design — "
        f"construct via SNMesh.from_axes (C5.5, #225), which passes "
        f"mesh=None at d≥3 and never calls this builder."
    )
