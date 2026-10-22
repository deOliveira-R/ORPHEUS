r"""The legacy-mesh adapter: translation between an axis tuple and a structured mesh.

:func:`axes_from_legacy_mesh` reads the per-axis edges and endpoint laws of a
:class:`~orpheus.mesh.structured.Mesh1D` or
:class:`~orpheus.mesh.structured.Mesh2D` into a tuple of
:class:`~orpheus.mesh.axis.AxisMesh` / :class:`~orpheus.mesh.axis.RadialAxisMesh`;
:func:`legacy_mesh_from_axes` builds the structured mesh an axis tuple
describes, the adapter mesh :meth:`~orpheus.sn.problem.SNProblem.from_axes`
synthesizes for the consumers still reading ``SNProblem.mesh``.

The adapter knows both the axis primitive (:mod:`orpheus.mesh.axis`) and the
structured mesh (:mod:`orpheus.mesh.structured`), so it carries more knowledge
than either and lives above both: it imports the two, and neither imports it.
It sat in :mod:`orpheus.mesh.axis` until the layer-DAG campaign
(``.claude/plans/layer_dag.md``), where its function-local import of the
structured mesh closed the cycle ``axis -> structured -> face_laws -> axis``.
The adapter narrows and retires with the legacy-mesh consumers (C5.2-C5.5,
#225); retiring it is deleting this module.
"""

from __future__ import annotations

import numpy as np

from orpheus.geometry import CoordSystem
from orpheus.geometry.boundary import BC
from orpheus.mesh.axis import (
    Axis1D,
    AxisCoord,
    AxisMesh,
    RadialAxisMesh,
    coord_system,
    face_labels,
    spatial_shape,
)
from orpheus.mesh.structured import Mesh1D, Mesh2D


# ═══════════════════════════════════════════════════════════════════════
# Legacy-mesh → axis-tuple adapter
# ═══════════════════════════════════════════════════════════════════════

def _refuse_an_inner_law_the_radial_axis_drops(mesh) -> None:
    """Refuse a hollow curvilinear mesh whose inner law is not reflective.

    **SCOPE-BOUNDARY[guard]** — machinery: an inner-surface trace on :class:`RadialAxisMesh` (#511).
    ruling: the user, 2026-09-29, P1 step 2 of ``.claude/plans/reference_cache.md``.
    revisit: the curvilinear SN build, deferred by the user, which gives the radial axis its inner slot.

    The radial axis has one law slot, ``bc_outer``, so the declared inner
    law of a hollow cylinder or sphere is dropped and every method built
    on the axes (SN, diffusion) computes a reflective cavity instead. The
    guard lands with the step that makes a hollow geometry declarable.
    Reflective is admitted because it is what the methods compute.
    """
    if "xmin" not in mesh.face_laws:
        return
    inner_law = mesh.face_laws["xmin"]
    if inner_law.kind != "reflective":
        raise NotImplementedError(
            f"a hollow {mesh.coord.name.lower()} mesh (r_0 = {float(mesh.edges[0])!r}) "
            f"declares the inner law {inner_law!r}, which the radial axis "
            f"cannot carry: it has no inner-surface trace, so the law would "
            f"be dropped and a reflective cavity computed instead (#511). "
            f"Only a reflective inner law is admitted."
        )


def axes_from_legacy_mesh(mesh) -> tuple[Axis1D, ...]:
    r"""Extract the canonical axis tuple from a legacy ``Mesh1D`` /
    ``Mesh2D``.

    Wraps the legacy mesh's per-axis edge arrays + BC declarations in
    :class:`AxisMesh` / :class:`RadialAxisMesh` so the
    :class:`~orpheus.sn.problem.SNProblem` constructor can expose
    :attr:`SNProblem.axes` regardless of which constructor surface the
    caller used.

    Mesh1D mapping (per coordinate system):

    * ``CARTESIAN`` → ``(AxisMesh(edges=mesh.edges,
      bc_low=mesh.face_laws["xmin"], bc_high=mesh.face_laws["xmax"]),)``.
    * ``SPHERICAL`` → ``(RadialAxisMesh(edges=mesh.edges,
      coord=RADIAL_SPHERICAL, bc_outer=mesh.outer_law),)``. A solid body
      has one face (the centre is an interior point and carries no law).
      On a hollow body the inner face is an endpoint the radial axis has
      no slot for, so only a reflective inner law is admitted, the one
      the methods compute; any other is refused (#511).
    * ``CYLINDRICAL`` → same as SPHERICAL with
      ``RADIAL_CYLINDRICAL``.

    Mesh2D mapping (Cartesian only — :class:`Mesh2D` cylindrical is
    not used by any SN sweep today):

    * ``CARTESIAN`` → ``(AxisMesh(edges=mesh.edges_x,
      bc_low=face_laws["xmin"], bc_high=face_laws["xmax"]),
      AxisMesh(edges=mesh.edges_y, bc_low=face_laws["ymin"],
      bc_high=face_laws["ymax"]))``.
    """
    if isinstance(mesh, Mesh1D):
        if mesh.coord == CoordSystem.CARTESIAN:
            return (
                AxisMesh(
                    edges=mesh.edges,
                    bc_low=mesh.face_laws["xmin"],
                    bc_high=mesh.face_laws["xmax"],
                ),
            )
        if mesh.coord == CoordSystem.SPHERICAL:
            _refuse_an_inner_law_the_radial_axis_drops(mesh)
            return (
                RadialAxisMesh(
                    edges=mesh.edges,
                    coord=AxisCoord.RADIAL_SPHERICAL,
                    bc_outer=mesh.outer_law,
                ),
            )
        if mesh.coord == CoordSystem.CYLINDRICAL:
            _refuse_an_inner_law_the_radial_axis_drops(mesh)
            return (
                RadialAxisMesh(
                    edges=mesh.edges,
                    coord=AxisCoord.RADIAL_CYLINDRICAL,
                    bc_outer=mesh.outer_law,
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
        laws = mesh.face_laws
        return (
            AxisMesh(edges=mesh.edges_x, bc_low=laws["xmin"], bc_high=laws["xmax"]),
            AxisMesh(edges=mesh.edges_y, bc_low=laws["ymin"], bc_high=laws["ymax"]),
        )

    raise TypeError(
        f"axes_from_legacy_mesh: expected Mesh1D or Mesh2D; "
        f"got {type(mesh).__name__}"
    )


# ═══════════════════════════════════════════════════════════════════════
# Axis-tuple → legacy-mesh ADAPTER builder (for SNProblem.from_axes)
# ═══════════════════════════════════════════════════════════════════════

def _face_laws_of_axis(ax, coord, edges):
    """The adapter mesh's face laws: the laws SN resolves on this axis.

    ``SNProblem.from_axes`` resolves its laws from the axes themselves; the
    adapter mesh only needs a declaration, so it carries the axis's declared
    laws, and on the inner face of a hollow radial axis, which has no law
    slot, the reflective cavity SN computes.
    SCOPE-BOUNDARY[guard] — machinery: an inner-surface trace on :class:`RadialAxisMesh` (#511).
    ruling: the user, 2026-09-29, P1 step 2 of ``.claude/plans/reference_cache.md``.
    revisit: the curvilinear SN build, deferred by the user.
    """
    if isinstance(ax, AxisMesh):
        return {"xmin": ax.bc_low, "xmax": ax.bc_high}
    points = coord.boundary_points(float(edges[0]), float(edges[-1]))
    if len(points) == 2:
        return {"xmin": BC("reflective"), "xmax": ax.bc_outer}
    return {"xmax": ax.bc_outer}


def legacy_mesh_from_axes(
    axes: tuple[Axis1D, ...],
    mat_map: np.ndarray | None = None,
):
    r"""Build the legacy :class:`Mesh1D` / :class:`Mesh2D` ADAPTER for an axis tuple.

    C5.1 (#225): this is NO LONGER a round-trip source — the axes an
    :class:`~orpheus.sn.problem.SNProblem` carries are the caller's
    objects verbatim. :meth:`SNProblem.from_axes` calls this builder only
    to synthesize the d≤2 ``SNProblem.mesh`` adapter for the consumers
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
    shape = spatial_shape(axes)
    if mat_map is None:
        mat_map = np.zeros(shape, dtype=int)
    else:
        mat_map = np.asarray(mat_map, dtype=int)
        # NOTE: SNProblem._init_core re-validates this shape — it is the
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
        coord = coord_system(axes)
        edges = np.asarray(ax.edges, dtype=float)
        # No geometry: each maximal run of one material is a region (labelling by
        # material id would fuse disjoint pieces of one material into one region).
        materials = mat_map.ravel()
        run_starts = np.r_[True, materials[1:] != materials[:-1]]
        return Mesh1D(
            coord=coord,
            edges=edges,
            volumes=coord.measure(edges),
            region_ids=np.cumsum(run_starts) - 1,
            region_materials=tuple(materials[run_starts].tolist()),
            face_laws=_face_laws_of_axis(ax, coord, edges),
        )

    if len(axes) == 2:
        ax0, ax1 = axes
        if ax0.coord != AxisCoord.CARTESIAN or ax1.coord != AxisCoord.CARTESIAN:
            raise NotImplementedError(
                f"legacy_mesh_from_axes: 2-D non-Cartesian (axes coords "
                f"{ax0.coord!r}, {ax1.coord!r}) has no Mesh2D today."
            )
        return Mesh2D(
            ax0.edges, ax1.edges, mat_map,
            face_laws={
                label.face_name: axes[label.axis_index].bc[label.endpoint]
                for label in face_labels(axes)
            },
            coord=CoordSystem.CARTESIAN,
        )

    raise NotImplementedError(
        f"legacy_mesh_from_axes: the legacy mesh ADAPTER is genuinely "
        f"d≤2 (Mesh1D / Mesh2D are the d≤2 user-facing dataclasses); "
        f"{len(axes)}-axis meshes are mesh-adapter-free by design — "
        f"construct via SNProblem.from_axes (C5.5, #225), which passes "
        f"mesh=None at d≥3 and never calls this builder."
    )
