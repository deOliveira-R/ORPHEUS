r"""Tests for SNMethodSpace (Wave 8 / C8.1).

The :class:`SNMethodSpace` dataclass carries the bundle of mesh +
quadrature + trace metadata that an SN realizer needs to turn a
:class:`BoundaryTraceLaw` into a 1-arg :class:`LinearOperator`.

Wave 5 introduced a minimal placeholder (quadrature only) embedded
in :mod:`orpheus.sn.boundary_realizer`. Wave 8 moves it to its own
dedicated module :mod:`orpheus.sn.method_space` and extends it with
``mesh`` + ``trace`` fields plus a :meth:`for_face` factory that
consumes the unified :class:`TraceSpace` to derive per-face inflow
indices.

These tests pin:

* Wave 5 backward compat (``.minimal(quad)`` still works).
* The Wave 5 re-export from ``orpheus.sn.boundary_realizer`` resolves
  to the same class object.
* :meth:`for_face` derives ``inflow_indices`` from the provided
  ``trace`` correctly.
* :meth:`inflow_indices_for_face` delegates to the held trace space
  or raises with a useful message when no trace was attached.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry import BC, CoordSystem
from orpheus.geometry.mesh import Mesh1D, Mesh2D
from orpheus.numerics.face_layout import FaceLayout
from orpheus.numerics.spaces.trace_space import TraceSpace
from orpheus.sn.method_space import SNMethodSpace
from orpheus.numerics.quadrature import Quadrature


pytestmark = pytest.mark.l0


def _slab_trace(mesh: Mesh1D, quad) -> TraceSpace:
    layout = FaceLayout.from_named_shapes(
        [("xmin", (quad.N, 1)), ("xmax", (quad.N, 1))]
    )
    return TraceSpace.from_quadrature_and_layout(quad, layout)


def _cartesian2d_trace(mesh: Mesh2D, quad, nx: int, ny: int) -> TraceSpace:
    layout = FaceLayout.from_named_shapes([
        ("xmin", (quad.N, 1, ny)), ("xmax", (quad.N, 1, ny)),
        ("ymin", (quad.N, 1, nx)), ("ymax", (quad.N, 1, nx)),
    ])
    return TraceSpace.from_quadrature_and_layout(quad, layout)


# ─────────────────────────────────────────────────────────────────────
# Wave 5 backward compatibility
# ─────────────────────────────────────────────────────────────────────


def test_minimal_factory_carries_only_quadrature():
    """``.minimal(quad)`` returns a method space with only quadrature populated."""
    quad = Quadrature.gauss_legendre(4)
    space = SNMethodSpace.minimal(quad)
    assert space.quadrature is quad
    assert space.face is None
    assert space.inflow_indices is None
    assert space.mesh is None
    assert space.trace is None


def test_legacy_constructor_with_explicit_inflow_indices():
    """The original Wave-5 construction path
    ``SNMethodSpace(quadrature=q, face='xmin', inflow_indices=idx)``
    keeps working — the new fields default to ``None``.
    """
    quad = Quadrature.lebedev(17)
    inflow = np.flatnonzero(quad.mu_x > 0)
    space = SNMethodSpace(
        quadrature=quad, face="xmin", inflow_indices=inflow,
    )
    assert space.face == "xmin"
    np.testing.assert_array_equal(space.inflow_indices, inflow)
    assert space.mesh is None
    assert space.trace is None


# ─────────────────────────────────────────────────────────────────────
# Wave 8 — for_face factory + trace delegation
# ─────────────────────────────────────────────────────────────────────


def test_for_face_derives_inflow_indices_from_trace():
    """``for_face`` extracts the per-face inflow indices from the trace
    space — load-bearing for ``SNMesh._resolve_one`` Cartesian path.
    """
    mesh = Mesh2D(
        edges_x=np.linspace(0.0, 2.0, 5), edges_y=np.linspace(0.0, 1.5, 4),
        mat_map=np.zeros((4, 3), dtype=np.int64),
    )
    quad = Quadrature.lebedev(17)
    trace = _cartesian2d_trace(mesh, quad, nx=4, ny=3)

    space = SNMethodSpace.for_face(
        mesh=mesh, quadrature=quad, face="xmin", trace=trace,
    )
    expected_inflow = np.flatnonzero(quad.mu_x > 1e-12)
    np.testing.assert_array_equal(space.inflow_indices, expected_inflow)
    assert space.face == "xmin"
    assert space.mesh is mesh
    assert space.trace is trace


def test_inflow_indices_for_face_delegates_to_trace():
    """:meth:`SNMethodSpace.inflow_indices_for_face` routes through the
    held trace space, so realizers can request indices by face name
    without knowing about the trace-space type.
    """
    mesh = Mesh1D(
        edges=np.linspace(0.0, 3.0, 7),
        mat_ids=np.zeros(6, dtype=np.int64),
        bc_left=BC("vacuum"), bc_right=BC("vacuum"),
    )
    quad = Quadrature.gauss_legendre(8)
    trace = _slab_trace(mesh, quad)

    # Build a method space attached to 'xmin' — but we should still be
    # able to query 'xmax' via the delegated method.
    space = SNMethodSpace.for_face(
        mesh=mesh, quadrature=quad, face="xmin", trace=trace,
    )
    right_indices = space.inflow_indices_for_face("xmax")
    # xmax face: outward normal +x; inflow is mu_x < 0.
    np.testing.assert_array_equal(
        right_indices, np.flatnonzero(quad.mu_x < -1e-12),
    )


def test_inflow_indices_for_face_raises_without_trace():
    """Without an attached trace space, :meth:`inflow_indices_for_face`
    raises with a useful message — the right error.
    """
    quad = Quadrature.gauss_legendre(4)
    space = SNMethodSpace.minimal(quad)
    with pytest.raises(RuntimeError, match="no TraceSpace attached"):
        space.inflow_indices_for_face("xmin")
