r"""The finite index sets: the 0-dimensional members of the manifold closed sum.

An :class:`IndexSet` is ``n`` points with no metric structure; an
:class:`EnergyGroups` is the multigroup energy axis, the counting set of ``ng``
groups. Both are :class:`~orpheus.geometry.manifold.Manifold` members (a measure
is supported on them as on a sphere), and neither is geometry: no rigid motion
acts on a group index, so neither is a :class:`~orpheus.geometry.manifold.GSpace`.

Why this module is in numerics and is named ``index_set`` (layer DAG WP-B,
2026-10-22, ruling 8 of ``.claude/plans/layer_dag.md``): the ``Manifold`` ABC and
its geometric members moved down into :mod:`orpheus.geometry.manifold`, because
a point set and the group acting on it are statements about the domain. A finite
index set is not a statement about the domain; it is the support of a counting
measure, which is numerics' vocabulary. The module is named for what it holds,
the finite index sets, of which the energy groups are the instance with a
physical identity. Geometry never imports it (the layer gate); it subclasses
geometry's ``Manifold`` across the boundary, the allowed direction.

Each member answers :attr:`~orpheus.geometry.manifold.Manifold.ambient_dim`
itself (one column), which is what let the geometric half stop naming these
classes: until WP-B the width was a module-level ``match`` with an arm
``Interval() | IndexSet() | EnergyGroups()``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike

from orpheus.geometry.manifold import Manifold

__all__ = ["ENERGY", "EnergyGroups", "IndexSet"]


@dataclass(frozen=True)
class IndexSet(Manifold):
    r"""A finite index set — ``n`` points with no metric structure.

    Two shipped sites minted this concept under incompatible spellings:
    ``frame.py``'s ``f"index({axis_label})"`` and
    ``loss_kernel_gauge.py``'s ``f"sn_trace_orbit{orbit}_g{group}"``,
    whose "points" are trace DOF indices cast to float.  They are one
    family; ``label`` is what distinguishes the instances.

    ``n`` is optional: a label alone identifies the set, and the
    cardinality is often the measure's own ``n_points``.
    """

    label: str
    n: int | None = None

    @property
    def dim(self) -> int:
        return 0

    @property
    def ambient_dim(self) -> int:
        return 1

    @property
    def name(self) -> str:
        return f"index({self.label})"

    def contains(self, points: ArrayLike) -> bool:
        arr = self._as_points(points, 1)[:, 0]
        if not np.all(arr == np.floor(arr)):
            return False
        if self.n is None:
            return bool(np.all(arr >= 0))
        return bool(np.all((arr >= 0) & (arr < self.n)))


@dataclass(frozen=True)
class EnergyGroups(Manifold):
    r"""The multigroup energy axis — a counting set of ``ng`` groups.

    Distinct from a bare :class:`IndexSet` because the measure's PHASE
    depends on it: an energy group index and any other integer-noded
    counting rule are indistinguishable from their nodes alone, which is
    exactly why the tag had to supply the physical identity.
    """

    ng: int | None = None

    @property
    def dim(self) -> int:
        return 0

    @property
    def ambient_dim(self) -> int:
        return 1

    @property
    def name(self) -> str:
        return "energy"

    def contains(self, points: ArrayLike) -> bool:
        arr = self._as_points(points, 1)[:, 0]
        if not np.all(arr == np.floor(arr)):
            return False
        if self.ng is None:
            return bool(np.all(arr >= 0))
        return bool(np.all((arr >= 0) & (arr < self.ng)))


#: The energy axis of unspecified group count, the support a multigroup
#: measure declares before it knows ``ng``.
ENERGY = EnergyGroups()
