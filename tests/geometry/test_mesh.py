"""Foundation tests for :mod:`orpheus.geometry.mesh` Issue 9.6 affordances.

Tests the ``Mesh1D.volume_measure`` and ``Mesh2D.volume_measure``
properties shipped in Phase B / Issue 9.6 — the natural integration
affordance that replaces hand-rolled ``np.sum(values * mesh.volumes)``
patterns at production integration sites.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry.coord import CoordSystem
from orpheus.geometry.mesh import Mesh1D, Mesh2D
from orpheus.numerics.measure import DiscreteMeasure


# ─────────────────────────────────────────────────────────────────────
# Mesh1D.volume_measure
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.foundation
def test_mesh1d_volume_measure_is_discrete_measure():
    """``mesh.volume_measure`` returns a DiscreteMeasure."""
    edges = np.linspace(0.0, 5.0, 6)
    mat_ids = np.zeros(5, dtype=int)
    mesh = Mesh1D(coord=CoordSystem.CARTESIAN, edges=edges, mat_ids=mat_ids)
    mu = mesh.volume_measure
    assert isinstance(mu, DiscreteMeasure)


@pytest.mark.foundation
def test_mesh1d_volume_measure_nodes_are_centers():
    edges = np.linspace(0.0, 5.0, 6)
    mat_ids = np.zeros(5, dtype=int)
    mesh = Mesh1D(coord=CoordSystem.CARTESIAN, edges=edges, mat_ids=mat_ids)
    np.testing.assert_array_equal(mesh.volume_measure.nodes, mesh.centers)


@pytest.mark.foundation
def test_mesh1d_volume_measure_weights_are_volumes():
    edges = np.linspace(0.0, 5.0, 6)
    mat_ids = np.zeros(5, dtype=int)
    mesh = Mesh1D(coord=CoordSystem.CARTESIAN, edges=edges, mat_ids=mat_ids)
    np.testing.assert_array_equal(mesh.volume_measure.weights, mesh.volumes)


@pytest.mark.foundation
def test_mesh1d_volume_measure_integrates_constant_to_total_volume():
    """``mu(1) == sum(volumes)``."""
    edges = np.linspace(0.0, 5.0, 6)
    mat_ids = np.zeros(5, dtype=int)
    mesh = Mesh1D(coord=CoordSystem.CARTESIAN, edges=edges, mat_ids=mat_ids)
    mu = mesh.volume_measure
    one = np.ones(mesh.N)
    assert mu(one) == pytest.approx(mesh.volumes.sum())


@pytest.mark.foundation
def test_mesh1d_volume_measure_array_overload_matches_explicit_sum():
    """``mu(values) == sum(volumes * values)`` for arbitrary values."""
    edges = np.linspace(0.0, 5.0, 6)
    mat_ids = np.zeros(5, dtype=int)
    mesh = Mesh1D(coord=CoordSystem.CARTESIAN, edges=edges, mat_ids=mat_ids)
    mu = mesh.volume_measure
    values = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    assert mu(values) == pytest.approx(np.sum(mesh.volumes * values))


@pytest.mark.foundation
def test_mesh1d_volume_measure_spherical():
    """Volume measure on a spherical mesh integrates ``1`` to total volume."""
    edges = np.linspace(0.0, 2.0, 5)
    mat_ids = np.zeros(4, dtype=int)
    mesh = Mesh1D(coord=CoordSystem.SPHERICAL, edges=edges, mat_ids=mat_ids)
    mu = mesh.volume_measure
    one = np.ones(mesh.N)
    assert mu(one) == pytest.approx(mesh.volumes.sum())


# ─────────────────────────────────────────────────────────────────────
# Mesh2D.volume_measure
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.foundation
def test_mesh2d_volume_measure_is_discrete_measure():
    edges_x = np.linspace(0.0, 3.0, 4)
    edges_y = np.linspace(0.0, 2.0, 3)
    mat_map = np.zeros((3, 2), dtype=int)
    mesh = Mesh2D(
        coord=CoordSystem.CYLINDRICAL,
        edges_x=edges_x,
        edges_y=edges_y,
        mat_map=mat_map,
    )
    mu = mesh.volume_measure
    assert isinstance(mu, DiscreteMeasure)


@pytest.mark.foundation
def test_mesh2d_volume_measure_node_layout_matches_meshgrid():
    """Nodes follow ``np.meshgrid(..., indexing='ij')`` order."""
    edges_x = np.linspace(0.0, 3.0, 4)
    edges_y = np.linspace(0.0, 2.0, 3)
    mat_map = np.zeros((3, 2), dtype=int)
    mesh = Mesh2D(
        coord=CoordSystem.CYLINDRICAL,
        edges_x=edges_x,
        edges_y=edges_y,
        mat_map=mat_map,
    )
    mu = mesh.volume_measure
    cx = mesh.centers_x
    cy = mesh.centers_y
    X, Y = np.meshgrid(cx, cy, indexing="ij")
    expected_nodes = np.stack([X.ravel(), Y.ravel()], axis=-1)
    np.testing.assert_array_equal(mu.nodes, expected_nodes)


@pytest.mark.foundation
def test_mesh2d_volume_measure_weights_match_flattened_volumes():
    edges_x = np.linspace(0.0, 3.0, 4)
    edges_y = np.linspace(0.0, 2.0, 3)
    mat_map = np.zeros((3, 2), dtype=int)
    mesh = Mesh2D(
        coord=CoordSystem.CYLINDRICAL,
        edges_x=edges_x,
        edges_y=edges_y,
        mat_map=mat_map,
    )
    mu = mesh.volume_measure
    np.testing.assert_array_equal(mu.weights, mesh.volumes.ravel())


@pytest.mark.foundation
def test_mesh2d_volume_measure_integrates_constant_to_total_volume():
    edges_x = np.linspace(0.0, 3.0, 4)
    edges_y = np.linspace(0.0, 2.0, 3)
    mat_map = np.zeros((3, 2), dtype=int)
    mesh = Mesh2D(
        coord=CoordSystem.CYLINDRICAL,
        edges_x=edges_x,
        edges_y=edges_y,
        mat_map=mat_map,
    )
    mu = mesh.volume_measure
    one = np.ones(mesh.volumes.size)
    assert mu(one) == pytest.approx(mesh.volumes.sum())
