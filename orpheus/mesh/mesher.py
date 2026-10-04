r"""The mesher: a meshing session on one geometry.

Load a :class:`~orpheus.geometry.StructuredGeometry`, partition it by
interval rules, refine, and read the mesh::

    mesh = Mesher(geometry).partition(CellsByCount.uniform_width(8)).mesh

The mesher holds its current mesh, the way a meshing program holds the mesh
it previews: :meth:`Mesher.partition` builds it, :meth:`Mesher.refine`
replaces it, and :attr:`Mesher.mesh` returns it. Quality measures, a preview,
adaptive refinement on a field and a protocol for external meshers are #539.

The mesher lifts the geometry onto the cells: each cell is labelled with the
interval it lies in (its region; the geometry's region → material map comes
with the labels), and each boundary face takes the geometry's law at that
boundary point. Every breakpoint is a cell edge, so every cell lies in exactly
one interval: the mesh refines the geometry's regions.
"""
from __future__ import annotations

import numpy as np

from orpheus.geometry.structured_geometry import StructuredGeometry
from orpheus.mesh.face_laws import face_inventory
from orpheus.mesh.partition import IntervalRule
from orpheus.mesh.structured import Mesh1D


class Mesher:
    r"""A meshing session on one geometry: partition by rules, refine, read the mesh.

    Parameters
    ----------
    geometry : StructuredGeometry
        The geometry to mesh.
    """

    def __init__(self, geometry: StructuredGeometry) -> None:
        if not isinstance(geometry, StructuredGeometry):
            raise TypeError(
                f"Mesher loads a StructuredGeometry, got {type(geometry).__name__}"
            )
        self._geometry = geometry
        self._rules: tuple[IntervalRule, ...] | None = None
        self._mesh: Mesh1D | None = None

    @property
    def geometry(self) -> StructuredGeometry:
        """The geometry this mesher meshes."""
        return self._geometry

    def partition(self, rule: "IntervalRule | tuple[IntervalRule, ...]") -> "Mesher":
        r"""Partition every interval by ``rule``, or interval ``k`` by ``rule[k]``.

        Builds the current mesh and returns the mesher, so a session chains.
        """
        intervals = self._geometry.intervals
        rules = rule if isinstance(rule, tuple) else (rule,) * len(intervals)
        if len(rules) != len(intervals):
            raise ValueError(
                f"Mesher.partition: {len(rules)} rule(s) for {len(intervals)} "
                f"interval(s); give one rule, or one per interval"
            )
        for k, r in enumerate(rules):
            if not isinstance(r, IntervalRule):
                raise TypeError(
                    f"Mesher.partition: rule {k} is an interval rule (CellsByCount, "
                    f"CellsByMaxWidth, Refined, CellEdges), got {type(r).__name__}"
                )
        cells = [r.cells(self._geometry, iv) for r, iv in zip(rules, intervals, strict=True)]
        # Every breakpoint is a cell edge, so every cell lies in one interval:
        # asserted here, since an interval rule is an open protocol.
        for k, ((edges_k, _), (a, b)) in enumerate(zip(cells, intervals, strict=True)):
            if edges_k[0] != a or edges_k[-1] != b:
                raise ValueError(
                    f"Mesher.partition: rule {k} ({type(rules[k]).__name__}) placed cells "
                    f"spanning [{edges_k[0]!r}, {edges_k[-1]!r}], not the interval "
                    f"[{a!r}, {b!r}]: every breakpoint must be a cell edge, bit for bit"
                )
        edges = np.concatenate([cells[0][0], *(e[1:] for e, _ in cells[1:])])
        volumes = np.concatenate([v for _, v in cells])
        # Each cell's region is the interval it partitions: the label the mesh keeps.
        region_ids = np.repeat(np.arange(len(cells)), [len(v) for _, v in cells])
        self._mesh = Mesh1D(
            coord=self._geometry.coord,
            edges=edges,
            volumes=volumes,
            region_ids=region_ids,
            region_materials=tuple(self._geometry.mat_ids),
            face_laws=dict(zip(
                face_inventory(self._geometry.coord, edges, 1),
                self._geometry.boundaries, strict=True,
            )),
        )
        self._rules = rules
        return self

    def refine(self, factor: int) -> "Mesher":
        r"""Refine the current partition: ``factor`` (a power of two) times the cells of each rule."""
        if self._rules is None:
            raise ValueError("Mesher.refine: partition the geometry first")
        return self.partition(tuple(factor * r for r in self._rules))

    @property
    def mesh(self) -> Mesh1D:
        """The current mesh."""
        if self._mesh is None:
            raise ValueError("Mesher.mesh: partition the geometry first")
        return self._mesh


__all__ = ["Mesher"]
