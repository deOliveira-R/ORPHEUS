r"""Tests for :mod:`orpheus.numerics.trace_space`.

Wave 2 of the ``transient-giggling-cake`` plan — trace function
spaces with per-face inflow / outflow masks.

L0 tests cover structural / type / equality / dtype / shape /
subset selection / error raising.

L1 tests cover the mathematical correctness of the directional
predicate :math:`\mathrm{sign}(\Omega \cdot \hat n_f)` mask
construction against hand-computed expectations.
"""

from __future__ import annotations

import dataclasses

import numpy as np
import pytest

from orpheus.geometry.coord import CoordSystem
from orpheus.geometry.mesh import Mesh1D, Mesh2D
from orpheus.numerics.space import FunctionSpace
from orpheus.numerics.trace_space import (
    InflowTraceSpace,
    OutflowTraceSpace,
    TraceSpace,
)
from orpheus.sn.quadrature import GaussLegendre1D, LebedevSphere


# ─────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────


def _mesh1d_cartesian(n: int = 4) -> Mesh1D:
    """Tiny 1-D slab mesh — Cartesian, n cells, uniform edges."""
    return Mesh1D(
        edges=np.linspace(0.0, 1.0, n + 1),
        mat_ids=np.zeros(n, dtype=int),
        coord=CoordSystem.CARTESIAN,
    )


def _mesh1d_spherical(n: int = 4) -> Mesh1D:
    """Tiny spherical mesh for the non-Cartesian-error test."""
    return Mesh1D(
        edges=np.linspace(0.0, 1.0, n + 1),
        mat_ids=np.zeros(n, dtype=int),
        coord=CoordSystem.SPHERICAL,
    )


def _mesh2d_cartesian(nx: int = 3, ny: int = 3) -> Mesh2D:
    """Tiny 2-D Cartesian mesh, (nx, ny) cells."""
    return Mesh2D(
        edges_x=np.linspace(0.0, 1.0, nx + 1),
        edges_y=np.linspace(0.0, 1.0, ny + 1),
        mat_map=np.zeros((nx, ny), dtype=int),
        coord=CoordSystem.CARTESIAN,
    )


# ─────────────────────────────────────────────────────────────────────
# Test 1 — L0: type construction + frozen invariance
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.l0
def test_inflow_trace_space_constructible_and_frozen():
    """L0: InflowTraceSpace constructs via factory; mutation raises."""
    mesh = _mesh2d_cartesian()
    quad = LebedevSphere.create(11)
    space = InflowTraceSpace.from_mesh_and_quadrature(mesh, quad, ng=1)

    # Subclass check.
    assert isinstance(space, InflowTraceSpace)
    assert isinstance(space, TraceSpace)
    assert isinstance(space, FunctionSpace)
    # Frozen: attempting to mutate any dataclass field raises.
    with pytest.raises(dataclasses.FrozenInstanceError):
        space.name = "evil"  # type: ignore[misc]
    with pytest.raises(dataclasses.FrozenInstanceError):
        space.shape = (1,)  # type: ignore[misc]


# ─────────────────────────────────────────────────────────────────────
# Test 2 — L0: equality + distinguishability
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.l0
def test_inflow_outflow_distinguishable_by_name():
    """L0: InflowTraceSpace and OutflowTraceSpace with same N/ng UNEQUAL."""
    mesh = _mesh2d_cartesian()
    quad = LebedevSphere.create(11)
    inflow = InflowTraceSpace.from_mesh_and_quadrature(mesh, quad, ng=1)
    outflow = OutflowTraceSpace.from_mesh_and_quadrature(mesh, quad, ng=1)
    assert inflow != outflow
    assert hash(inflow) != hash(outflow)


@pytest.mark.l0
def test_inflow_trace_space_equality_independent_of_mask():
    """L0: two InflowTraceSpaces with same N/ng compare equal — mask not in __eq__."""
    mesh = _mesh2d_cartesian()
    quad = LebedevSphere.create(11)
    a = InflowTraceSpace.from_mesh_and_quadrature(mesh, quad, ng=1)
    b = InflowTraceSpace.from_mesh_and_quadrature(mesh, quad, ng=1)
    assert a == b
    assert hash(a) == hash(b)
    # And the masks are present + equal-content.
    np.testing.assert_array_equal(a.inflow_mask, b.inflow_mask)


# ─────────────────────────────────────────────────────────────────────
# Test 3 — L1: from_mesh_and_quadrature on 2-D Cartesian + Lebedev
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.l1
def test_mesh2d_cartesian_lebedev_inflow_mask_per_face():
    """L1: per-face inflow predicate matches mu_x / mu_y signs."""
    mesh = _mesh2d_cartesian()
    quad = LebedevSphere.create(11)
    eps = 1e-12
    space = InflowTraceSpace.from_mesh_and_quadrature(mesh, quad, ng=1)
    mask = space.inflow_mask
    assert mask.shape == (4, quad.N)
    # Faces in canonical order: xmin, xmax, ymin, ymax.
    # bc_xmin: outward -x → inflow iff mu_x > eps
    np.testing.assert_array_equal(mask[0], quad.mu_x > eps)
    # bc_xmax: outward +x → inflow iff mu_x < -eps
    np.testing.assert_array_equal(mask[1], quad.mu_x < -eps)
    # bc_ymin: outward -y → inflow iff mu_y > eps
    np.testing.assert_array_equal(mask[2], quad.mu_y > eps)
    # bc_ymax: outward +y → inflow iff mu_y < -eps
    np.testing.assert_array_equal(mask[3], quad.mu_y < -eps)


# ─────────────────────────────────────────────────────────────────────
# Test 4 — L1: mask consistency with tangential ordinates excluded
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.l1
def test_lebedev_axis_aligned_ordinates_excluded_from_both_masks():
    """L1: pure-axis ordinates (e.g. (0,0,±1)) appear in NEITHER mask
    on a face perpendicular to them."""
    # Lebedev order=3 has 6 ordinates: (±1,0,0), (0,±1,0), (0,0,±1) —
    # all axis-aligned. Order=11 still includes the 6 axis ordinates.
    mesh = _mesh2d_cartesian()
    quad = LebedevSphere.create(11)
    inflow = InflowTraceSpace.from_mesh_and_quadrature(mesh, quad, ng=1)
    outflow = OutflowTraceSpace.from_mesh_and_quadrature(mesh, quad, ng=1)
    # For face "xmin" (perpendicular to z-axis ordinates), the
    # (0, 0, ±1) ordinates have Ω · n = 0 and must be excluded.
    z_axis_indices = np.where(
        (np.abs(quad.mu_x) < 1e-12) & (np.abs(quad.mu_y) < 1e-12)
        & (np.abs(quad.mu_z) > 0.5)
    )[0]
    # Sanity: Lebedev 11 has axis-aligned ordinates.
    assert z_axis_indices.size >= 2, "expected (0,0,±1) ordinates"
    # For face "xmin" (idx 0), these ordinates are tangential.
    assert not inflow.inflow_mask[0, z_axis_indices].any()
    assert not outflow.outflow_mask[0, z_axis_indices].any()


# ─────────────────────────────────────────────────────────────────────
# Test 5 — L0: 1-D slab edge case
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.l0
def test_mesh1d_cartesian_gausslegendre():
    """L0: Mesh1D Cartesian + GL produces (2, N) inflow mask;
    bc_left has the 4 positive-mu_x ordinates as inflow."""
    mesh = _mesh1d_cartesian()
    quad = GaussLegendre1D.create(n_ordinates=8)
    space = InflowTraceSpace.from_mesh_and_quadrature(mesh, quad, ng=1)
    assert space.inflow_mask.shape == (2, 8)
    # bc_left: inflow iff mu_x > 0. GL on (-1, 1) with N=8 has 4
    # positive nodes.
    inflow_left = space.inflow_mask[0]
    assert inflow_left.sum() == 4
    np.testing.assert_array_equal(inflow_left, quad.mu_x > 1e-12)
    # bc_right: inflow iff mu_x < 0.
    inflow_right = space.inflow_mask[1]
    assert inflow_right.sum() == 4
    np.testing.assert_array_equal(inflow_right, quad.mu_x < -1e-12)


# ─────────────────────────────────────────────────────────────────────
# Test 6 — L0: inflow_indices_for_face returns the right indices
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.l0
def test_inflow_indices_for_face_cross_check_against_mask():
    """L0: inflow_indices_for_face("xmin") matches np.flatnonzero of mask row."""
    mesh = _mesh2d_cartesian()
    quad = LebedevSphere.create(11)
    space = InflowTraceSpace.from_mesh_and_quadrature(mesh, quad, ng=1)
    indices = space.inflow_indices_for_face("xmin")
    expected = np.flatnonzero(space.inflow_mask[0])
    np.testing.assert_array_equal(indices, expected)
    # And for "ymax" (mask row 3).
    indices_ymax = space.inflow_indices_for_face("ymax")
    expected_ymax = np.flatnonzero(space.inflow_mask[3])
    np.testing.assert_array_equal(indices_ymax, expected_ymax)


# ─────────────────────────────────────────────────────────────────────
# Test 7 — L0: goodness of mask — sums bounded by N
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.l0
def test_inflow_outflow_sums_bounded_by_n_ordinates():
    """L0: inflow_mask[f].sum() + outflow_mask[f].sum() <= N (tangentials excluded)."""
    mesh = _mesh2d_cartesian()
    quad = LebedevSphere.create(11)
    inflow = InflowTraceSpace.from_mesh_and_quadrature(mesh, quad, ng=1)
    outflow = OutflowTraceSpace.from_mesh_and_quadrature(mesh, quad, ng=1)
    for f_idx in range(4):
        total = (
            int(inflow.inflow_mask[f_idx].sum())
            + int(outflow.outflow_mask[f_idx].sum())
        )
        assert total <= quad.N


# ─────────────────────────────────────────────────────────────────────
# Test 8 — L0: non-Cartesian Mesh1D raises NotImplementedError
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.l0
def test_non_cartesian_mesh1d_raises_notimplemented():
    """L0: spherical Mesh1D raises NotImplementedError with a
    grep-able deferral message."""
    mesh = _mesh1d_spherical()
    quad = GaussLegendre1D.create(n_ordinates=4)
    with pytest.raises(NotImplementedError, match="curvilinear"):
        InflowTraceSpace.from_mesh_and_quadrature(mesh, quad, ng=1)
    with pytest.raises(NotImplementedError, match="curvilinear"):
        OutflowTraceSpace.from_mesh_and_quadrature(mesh, quad, ng=1)


# ─────────────────────────────────────────────────────────────────────
# Test 9 — L0: subset of faces
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.l0
def test_subset_of_faces_produces_smaller_mask():
    """L0: faces=["xmin"] → inflow_mask shape (1, N)."""
    mesh = _mesh2d_cartesian()
    quad = LebedevSphere.create(11)
    space = InflowTraceSpace.from_mesh_and_quadrature(
        mesh, quad, faces=["xmin"], ng=1,
    )
    assert space.inflow_mask.shape == (1, quad.N)
    assert space.face_names == ("xmin",)


# ─────────────────────────────────────────────────────────────────────
# Test 10 — L0: ng > 1 sets shape correctly
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.l0
def test_ng_greater_than_one_in_shape():
    """L0: ng=4 yields shape == (n_ordinates, 4)."""
    mesh = _mesh2d_cartesian()
    quad = LebedevSphere.create(11)
    space = InflowTraceSpace.from_mesh_and_quadrature(mesh, quad, ng=4)
    assert space.shape == (quad.N, 4)


# ─────────────────────────────────────────────────────────────────────
# Test 11 — L0: inflow XOR outflow on non-tangential quadrature
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.l0
def test_inflow_xor_outflow_complementary_for_gl_1d():
    """L0: for GL 1-D (no axis-aligned-to-x ordinates at standard N),
    inflow XOR outflow is True on every ordinate per face."""
    # GL ordinates are strictly in (-1, 1), never 0 → no tangentials
    # against the x-aligned faces.
    mesh = _mesh1d_cartesian()
    quad = GaussLegendre1D.create(n_ordinates=8)
    inflow = InflowTraceSpace.from_mesh_and_quadrature(mesh, quad, ng=1)
    outflow = OutflowTraceSpace.from_mesh_and_quadrature(mesh, quad, ng=1)
    for f_idx in range(2):
        complement = inflow.inflow_mask[f_idx] ^ outflow.outflow_mask[f_idx]
        assert complement.all(), f"face {f_idx} fails XOR-complementarity"


# ─────────────────────────────────────────────────────────────────────
# Test 12 — L0: mask dtype is bool
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.l0
def test_mask_dtype_is_bool():
    """L0: inflow_mask / outflow_mask are bool arrays, not int."""
    mesh = _mesh2d_cartesian()
    quad = LebedevSphere.create(11)
    inflow = InflowTraceSpace.from_mesh_and_quadrature(mesh, quad, ng=1)
    outflow = OutflowTraceSpace.from_mesh_and_quadrature(mesh, quad, ng=1)
    assert inflow.inflow_mask.dtype == np.bool_
    assert outflow.outflow_mask.dtype == np.bool_


# ─────────────────────────────────────────────────────────────────────
# Bonus — L0: outflow factory + index method
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.l0
def test_outflow_indices_for_face_cross_check_against_mask():
    """L0: outflow_indices_for_face matches np.flatnonzero of mask."""
    mesh = _mesh2d_cartesian()
    quad = LebedevSphere.create(11)
    space = OutflowTraceSpace.from_mesh_and_quadrature(mesh, quad, ng=1)
    for f_idx, face in enumerate(("xmin", "xmax", "ymin", "ymax")):
        indices = space.outflow_indices_for_face(face)
        expected = np.flatnonzero(space.outflow_mask[f_idx])
        np.testing.assert_array_equal(indices, expected)


@pytest.mark.l0
def test_unknown_face_raises_value_error():
    """L0: inflow_indices_for_face('bogus') raises ValueError."""
    mesh = _mesh2d_cartesian()
    quad = LebedevSphere.create(11)
    inflow = InflowTraceSpace.from_mesh_and_quadrature(mesh, quad, ng=1)
    with pytest.raises(ValueError, match="Unknown face"):
        inflow.inflow_indices_for_face("bogus")
