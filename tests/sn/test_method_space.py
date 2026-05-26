r"""Tests for SNMethodSpace (Wave 8 / C8.1).

The :class:`SNMethodSpace` dataclass carries the bundle of mesh +
quadrature + trace metadata that an SN realizer needs to turn a
:class:`BoundaryTraceLaw` into a 1-arg :class:`LinearOperator`.

Wave 5 introduced a minimal placeholder (quadrature only) embedded
in :mod:`orpheus.sn.boundary_realizer`. Wave 8 moves it to its own
dedicated module :mod:`orpheus.sn.method_space` and extends it with
``mesh`` + ``inflow_trace`` + ``outflow_trace`` fields plus a
:meth:`for_face` factory that consumes a precomputed
:class:`InflowTraceSpace` to derive per-face inflow indices.

These tests pin:

* Wave 5 backward compat (``.minimal(quad)`` still works).
* The Wave 5 re-export from ``orpheus.sn.boundary_realizer`` resolves
  to the same class object.
* :meth:`for_face` derives ``inflow_indices`` from the provided
  ``inflow_trace`` correctly.
* :meth:`inflow_indices_for_face` delegates to the held trace space
  or raises with a useful message when no trace was attached.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry import BC, CoordSystem
from orpheus.geometry.mesh import Mesh1D, Mesh2D
from orpheus.numerics.trace_space import InflowTraceSpace, OutflowTraceSpace
from orpheus.sn.method_space import SNMethodSpace
from orpheus.numerics.quadrature import Quadrature


pytestmark = pytest.mark.l0


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
    assert space.inflow_trace is None
    assert space.outflow_trace is None


def test_realizer_module_reexports_same_class():
    """``from orpheus.sn.boundary_realizer import SNMethodSpace`` resolves
    to the same class object as the canonical module — Wave 5 import
    sites stay working unchanged.
    """
    from orpheus.sn.boundary_realizer import SNMethodSpace as ReexportedClass

    assert ReexportedClass is SNMethodSpace


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
    assert space.inflow_trace is None


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
    inflow_trace = InflowTraceSpace.from_mesh_and_quadrature(mesh, quad)
    outflow_trace = OutflowTraceSpace.from_mesh_and_quadrature(mesh, quad)

    space = SNMethodSpace.for_face(
        mesh=mesh, quadrature=quad, face="xmin",
        inflow_trace=inflow_trace, outflow_trace=outflow_trace,
    )
    expected_inflow = np.flatnonzero(quad.mu_x > 1e-12)
    np.testing.assert_array_equal(space.inflow_indices, expected_inflow)
    assert space.face == "xmin"
    assert space.mesh is mesh
    assert space.inflow_trace is inflow_trace
    assert space.outflow_trace is outflow_trace


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
    inflow_trace = InflowTraceSpace.from_mesh_and_quadrature(mesh, quad)

    # Build a method space attached to 'left' — but we should still be
    # able to query 'right' via the delegated method.
    space = SNMethodSpace.for_face(
        mesh=mesh, quadrature=quad, face="left", inflow_trace=inflow_trace,
    )
    right_indices = space.inflow_indices_for_face("right")
    # Right face: outward normal +x; inflow is mu_x < 0.
    np.testing.assert_array_equal(
        right_indices, np.flatnonzero(quad.mu_x < -1e-12),
    )


def test_inflow_indices_for_face_raises_without_trace():
    """Without an attached trace space, :meth:`inflow_indices_for_face`
    raises with a useful message — the right error.
    """
    quad = Quadrature.gauss_legendre(4)
    space = SNMethodSpace.minimal(quad)
    with pytest.raises(RuntimeError, match="no InflowTraceSpace attached"):
        space.inflow_indices_for_face("left")
