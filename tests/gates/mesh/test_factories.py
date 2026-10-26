"""The mesh factories of :mod:`orpheus.mesh.factories`.

:func:`~orpheus.mesh.pwr_pin_2d` builds the 2-D Cartesian PWR pin cell, a
:class:`~orpheus.mesh.Mesh2D`: its shape, its material map (fuel, cladding and
coolant all present), its flat material ids, and its refusal of a material-id
list whose length does not match the radii.

``TestPWRPin2D`` and ``TestPWRPin2DEdgeCases`` lived in
``tests/gates/geometry/test_geometry.py`` until the layer-DAG re-filing
(`.claude/plans/layer_dag.md`, ruling 21): a claim about an ``orpheus.mesh``
object is stated in ``mesh`` vocabulary. Every test is ``foundation``.
"""

from __future__ import annotations

import pytest

from orpheus.geometry import MirrorDeck
from orpheus.mesh import pwr_pin_2d

pytestmark = pytest.mark.foundation


# ═══════════════════════════════════════════════════════════════════════
# pwr_pin_2d — 2-D Cartesian PWR pin factory (post-Phase-F survivor)
# ═══════════════════════════════════════════════════════════════════════

class TestPWRPin2D:
    """:func:`pwr_pin_2d` 2-D factory invariants."""

    def test_pin_2d_shape(self):
        mesh = pwr_pin_2d(n_cells=10, role=MirrorDeck())
        assert mesh.nx == 10
        assert mesh.ny == 10
        assert mesh.mat_map.shape == (10, 10)

    def test_pin_2d_has_all_materials(self):
        mesh = pwr_pin_2d(n_cells=20, role=MirrorDeck())
        mats = set(mesh.mat_map.ravel())
        assert mats == {0, 1, 2}

    def test_pin_2d_mat_ids_flat(self):
        """mat_ids returns flat array for assemble_cell_xs."""
        mesh = pwr_pin_2d(n_cells=5, role=MirrorDeck())
        assert mesh.mat_ids.shape == (25,)


# ═══════════════════════════════════════════════════════════════════════
# pwr_pin_2d edge cases
# ═══════════════════════════════════════════════════════════════════════

class TestPWRPin2DEdgeCases:
    """Edge cases for :func:`pwr_pin_2d`."""

    def test_pin_2d_wrong_mat_ids_length_raises(self):
        with pytest.raises(ValueError, match="len\\(mat_ids\\)"):
            pwr_pin_2d(radii=[1.0], mat_ids=[0, 1, 2], role=MirrorDeck())
