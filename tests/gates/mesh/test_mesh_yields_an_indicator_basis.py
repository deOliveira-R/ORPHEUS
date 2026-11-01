"""MOVED from tests/gates/numerics/test_indicator_basis.py (layer DAG Phase 2 refresh, specs_2r §5). Bodies unchanged."""

from __future__ import annotations
import numpy as np
import pytest
from orpheus.geometry import MirrorDeck, StructuredGeometry
from orpheus.mesh import CellsByCount, Mesher
from orpheus.numerics.basis import Basis


pytestmark = [pytest.mark.foundation]








_COARSE_EDGES = np.array([0.0, 2.0, 4.0])               # 2 coarse cells










def test_mesh_is_not_a_basis_but_yields_one():
    """The mesh YIELDS a basis view; it does NOT inherit Basis (role separation)."""
    geometry = StructuredGeometry.slab(
        _COARSE_EDGES, (0, 1), left=MirrorDeck(), right=MirrorDeck(),
    )
    mesh = Mesher(geometry).partition(CellsByCount.uniform_width(1)).mesh
    assert not isinstance(mesh, Basis)
    assert isinstance(mesh.indicator_basis(), Basis)
    # the yielded basis is built from the mesh's own edges
    np.testing.assert_array_equal(
        mesh.indicator_basis().edges_per_axis[0], _COARSE_EDGES,
    )




def test_the_mesh_frames_two_halves_name_one_manifold() -> None:
    """The mesh's basis view and its volume measure name ONE manifold (split from numerics'
    ``test_basis_domain.py::test_d6_a_frames_two_halves_name_ONE_manifold`` by the layer DAG Phase 2 refresh:
    the producer is the mesh's). Before 2.1 the measure named the manifold while the basis hard-coded a spatial one."""
    geometry = StructuredGeometry.slab((0.0, 0.5, 1.0), (0, 1), left=0, right=0)
    mesh = Mesher(geometry).partition(CellsByCount.uniform_width(1)).mesh
    basis, measure = mesh.indicator_basis(), mesh.volume_measure
    assert basis.domain == measure.support, (
        f"Mesh1D: the frame's basis says it lives on {basis.domain!r} and its measure says {measure.support!r}"
    )
