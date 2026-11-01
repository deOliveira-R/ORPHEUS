r"""``SNProblem`` is-a ``MaterialMesh``: the data/behaviour split, seen from the SN problem.

Moved from ``tests/gates/transport/test_material_mesh.py`` by layer DAG Phase 4 (2026-11-01): the three rows read an
``SNProblem`` (Liskov substitutability, the inherited data block bit-identical to a standalone ``MaterialMesh``, and
the ``from_material_mesh`` promotion), so they are SN's claims; the transport file keeps the carrier's own data
contract. The two-material fixture is that file's, imported.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.problem import SNProblem
from orpheus.transport.mesh import MaterialMesh, MaterialXSField
from tests.gates.transport.test_material_mesh import _two_material_mesh

pytestmark = pytest.mark.foundation


def test_snmesh_is_a_material_mesh():
    assert issubclass(SNProblem, MaterialMesh)


def test_snmesh_data_block_bit_identical_to_standalone_carrier():
    """An SNProblem's inherited data block matches a standalone MaterialMesh
    built from the same mesh + materials (the split is bit-identical)."""
    mesh, mats = _two_material_mesh()
    quad = Quadrature.gauss_legendre(n_ordinates=8)
    snm = SNProblem(mesh, quad, mats)
    mm = MaterialMesh(mesh, mats)
    assert snm.ng == mm.ng
    assert snm.ndim == mm.ndim
    assert snm.spatial_shape == mm.spatial_shape
    np.testing.assert_array_equal(snm.mat_map, mm.mat_map)
    np.testing.assert_array_equal(snm.volumes, mm.volumes)
    np.testing.assert_array_equal(
        snm.volume_measure.weights, mm.volume_measure.weights,
    )
    # Every carrier accessor is callable on the SNProblem (substitutability).
    assert isinstance(snm.mat_xs, MaterialXSField)


def test_from_material_mesh_round_trips_to_solvable_snmesh():
    mesh, mats = _two_material_mesh()
    quad = Quadrature.gauss_legendre(n_ordinates=8)
    mm = MaterialMesh(mesh, mats)
    snm = SNProblem.from_material_mesh(mm, quad)
    assert isinstance(snm, SNProblem)
    # carries the carrier's data verbatim …
    np.testing.assert_array_equal(snm.mat_map, mm.mat_map)
    np.testing.assert_array_equal(snm.volumes, mm.volumes)
    assert snm.materials is mm.materials
    # … plus the SN method layer (quadrature + streaming stencil).
    assert snm.quad is quad
    assert snm.reduced is not None
