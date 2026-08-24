r"""Foundation: L2 :class:`AngularBoundaryFlux.face_view` writability contract.

Test-architect Test 2.6 (D-H.2-C4a verification spec).  Pins the
software-level invariant that
:meth:`orpheus.transport.fields.angular_boundary_flux.AngularBoundaryFlux.face_view`
returns a writable view into the flat backing buffer — writes
through the view MUST propagate to the stored ``values`` array.

If ``face_view`` ever silently returns a copy (e.g. via an upstream
NumPy refactor that breaks the slice/reshape view chain in
:class:`FaceSlot.slice_view`), the SN sweep's persistent partner-flux
plumbing fails silently — every BC apply on a "writable" view would
no-op the persistence, the reflective BC fixed point would drift to
vacuum, and the surface eventually catches it as a divergent k_eff
or a NaN, but not at the source.  This foundation test catches the
break immediately at the API boundary.

Why this catches the bug class
==============================

* Pattern 4 (illegal-states-unrepresentable) — ``face_view`` claims
  to return a writable view; the type system cannot enforce that.
  A positive integration test is the only mechanism.
* No legacy ``AngularBoundaryFlux`` analog — the legacy class exposed
  per-face arrays as mutable attributes (``bf.xmin_face = ...``),
  so the writability question was trivially answered by Python's
  attribute model.  The L2 carve replaces attribute mutation with
  ``view = bf.face_view(name); view[...] = ...``; this is the test
  that pins the new contract.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry import BC, CoordSystem, Mesh1D, Mesh2D
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.mesh.augmented_mesh import SNMesh
from tests.sn._test_helpers import placeholder_materials
from orpheus.transport.fields.angular_boundary_flux import AngularBoundaryFlux

pytestmark = pytest.mark.foundation


# ── Mesh fixtures spanning every supported layout ───────────────────────


def _slab_mesh() -> SNMesh:
    geom = Mesh1D(
        edges=np.linspace(0.0, 1.0, 5),
        mat_ids=np.zeros(4, dtype=int),
        coord=CoordSystem.CARTESIAN,
        bc_left=BC("vacuum"),
        bc_right=BC("vacuum"),
    )
    quad = Quadrature.gauss_legendre(n_ordinates=4)
    return SNMesh(geom, quad, placeholder_materials(ng=2))


def _spherical_mesh() -> SNMesh:
    geom = Mesh1D(
        edges=np.linspace(0.0, 1.0, 5),
        mat_ids=np.zeros(4, dtype=int),
        coord=CoordSystem.SPHERICAL,
        bc_left=BC("reflective"),
        bc_right=BC("vacuum"),
    )
    quad = Quadrature.gauss_legendre(n_ordinates=4)
    return SNMesh(geom, quad, placeholder_materials(ng=2))


def _cartesian_2d_mesh() -> SNMesh:
    geom = Mesh2D(
        edges_x=np.linspace(0.0, 1.0, 4),
        edges_y=np.linspace(0.0, 1.0, 4),
        mat_map=np.zeros((3, 3), dtype=int),
        bc_xmin=BC("vacuum"),
        bc_xmax=BC("vacuum"),
        bc_ymin=BC("vacuum"),
        bc_ymax=BC("vacuum"),
    )
    quad = Quadrature.level_symmetric(sn_order=4)
    return SNMesh(geom, quad, placeholder_materials(ng=2))


GEOMETRY_FACE_TABLE = [
    ("slab_1d", _slab_mesh, ("xmin", "xmax")),
    ("sphere_1d", _spherical_mesh, ("xmax",)),
    ("cartesian_2d", _cartesian_2d_mesh, ("xmin", "xmax", "ymin", "ymax")),
]


# ── The writability contract ────────────────────────────────────────────


@pytest.mark.parametrize("name,mesh_builder,faces", GEOMETRY_FACE_TABLE)
def test_face_view_inplace_mutation_propagates_to_backing(
    name: str, mesh_builder, faces: tuple[str, ...],
) -> None:
    r"""Writes through ``face_view`` must propagate to the flat backing.

    Per-face round-trip test: write a recognisable sentinel via
    ``face_view(face)[...] = sentinel``, then re-fetch via
    ``face_view(face)`` and verify the sentinel reads back.  Catches
    a regression where ``face_view`` returns a copy (writes lost) or
    where ``FaceSlot.slice_view`` doesn't preserve the view chain
    (writes silently propagate to scratch).
    """
    mesh = mesh_builder()
    boundary = AngularBoundaryFlux.zeros(mesh.angular_trace)

    sentinel = 7.5
    for face in faces:
        view = boundary.face_view(face)
        view[...] = sentinel
        view_read = boundary.face_view(face)
        assert np.allclose(view_read, sentinel), (
            f"face_view({face!r}) writability broken on {name}: "
            f"wrote {sentinel}, read back {view_read.ravel()[:3]}..."
        )


@pytest.mark.parametrize("name,mesh_builder,faces", GEOMETRY_FACE_TABLE)
def test_face_view_writes_appear_in_flat_values(
    name: str, mesh_builder, faces: tuple[str, ...],
) -> None:
    r"""Writes through ``face_view`` must propagate to ``boundary.values``.

    The flat backing is the single source of truth; ``face_view``
    returns a structured view into it.  Writes through any face's
    view must be observable in the flat ``values`` array — otherwise
    downstream consumers reading ``boundary.values`` directly (e.g.
    the ``_is_ravellable`` / ``to_flat`` adapter chain) see stale
    data.
    """
    mesh = mesh_builder()
    boundary = AngularBoundaryFlux.zeros(mesh.angular_trace)

    # Distinct sentinel per face to disambiguate the slot mappings.
    for i, face in enumerate(faces):
        sentinel = 1.0 + i * 0.25  # 1.0, 1.25, 1.5, 1.75 — bit-distinct
        boundary.face_view(face)[...] = sentinel

    # Re-read everything via the flat buffer; sentinel values must be
    # present at the layout's per-face offsets.
    layout = boundary.layout
    for i, face in enumerate(faces):
        sentinel = 1.0 + i * 0.25
        slot = layout.faces[face]
        flat_slice = boundary.values[slot.offset:slot.offset + slot.flat_size]
        assert np.allclose(flat_slice, sentinel), (
            f"Flat backing slice for face {face!r} on {name}: expected "
            f"{sentinel}, got {flat_slice[:3]}..."
        )


@pytest.mark.parametrize("name,mesh_builder,faces", GEOMETRY_FACE_TABLE)
def test_face_view_returns_distinct_views_per_face(
    name: str, mesh_builder, faces: tuple[str, ...],
) -> None:
    r"""Different face names must yield disjoint views into the flat backing.

    Catches a layout regression where two face slots accidentally
    alias (e.g. an off-by-one in ``FaceLayout.from_named_shapes``
    offset computation).  Write a distinct sentinel to each face;
    verify each face reads back ITS sentinel, not its neighbour's.
    """
    if len(faces) < 2:
        pytest.skip(f"{name} has only one face; alias check trivially holds")

    mesh = mesh_builder()
    boundary = AngularBoundaryFlux.zeros(mesh.angular_trace)

    sentinels = {face: 1.0 + i * 0.25 for i, face in enumerate(faces)}
    for face, sentinel in sentinels.items():
        boundary.face_view(face)[...] = sentinel

    for face, expected in sentinels.items():
        view = boundary.face_view(face)
        assert np.allclose(view, expected), (
            f"Face {face!r} on {name}: expected sentinel {expected}, "
            f"read {view.ravel()[:3]}... — slot aliasing or offset bug."
        )


def test_face_view_raises_on_unknown_face_name() -> None:
    r"""``face_view`` MUST reject unknown face names at the API boundary."""
    mesh = _slab_mesh()
    boundary = AngularBoundaryFlux.zeros(mesh.angular_trace)
    with pytest.raises(KeyError, match="no face keyed"):
        boundary.face_view("not_a_real_face")


def test_face_view_writable_independent_of_geometry_face_count() -> None:
    r"""1-D curvilinear (single face) and 1-D slab (two faces) both writable.

    Pins that the single-face layout (curvilinear) handles the
    writability contract identically to the multi-face layouts.
    A regression here would break SI on sphere/cylinder.
    """
    sphere_mesh = _spherical_mesh()
    boundary = AngularBoundaryFlux.zeros(sphere_mesh.angular_trace)
    view = boundary.face_view("xmax")
    view[...] = 3.14
    assert np.allclose(boundary.face_view("xmax"), 3.14)
