"""Membership of the index sets: ``IndexSet.contains`` and ``EnergyGroups.contains``.

Re-filed 2026-10-25 from ``tests/gates/geometry/test_manifold.py`` (layer-DAG ruling 23):
:class:`~orpheus.numerics.index_set.IndexSet` and its energy specialisation are
numerics members of the ``Manifold`` sum, so their membership rows live beside
their module; the geometry members' membership rows stay in geometry's
``TestMembership``.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.numerics.index_set import EnergyGroups, IndexSet


@pytest.mark.foundation
class TestMembership:
    """The numerics half of geometry's ``TestMembership``.

    Re-filed 2026-10-25 from ``tests/gates/geometry/test_manifold.py`` (layer-DAG
    ruling 23): the rows below assert a numerics object's verb, so they live here; the class keeps its geometry name so the node id
    changes only in its path.  The manifold rows of the class stay in geometry.
    """

    def test_index_sets_admit_integers_only(self):
        assert IndexSet("angular", 8).contains(np.arange(8))
        assert not IndexSet("angular", 8).contains(np.array([8.0]))
        assert not IndexSet("angular", 8).contains(np.array([0.5]))
        assert not IndexSet("angular").contains(np.array([-1.0]))

    def test_energy_groups_admit_group_indices(self):
        assert EnergyGroups(2).contains(np.array([0.0, 1.0]))
        assert not EnergyGroups(2).contains(np.array([2.0]))
