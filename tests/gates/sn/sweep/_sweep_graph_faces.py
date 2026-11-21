r"""Test-side 2-D face-offset accessors of :class:`SweepDependencyGraph` (moved off production).

``face_in_x``, ``face_out_x``, ``face_in_y`` and ``face_out_y`` were four 2-D
convenience properties of the graph with no production reader: the generic walks read
``face_in[a]`` and ``face_out[a]`` per axis. Each is now a free function of the graph.
"""
from __future__ import annotations

from orpheus.sn.loss_representation.sweep_graph import SweepDependencyGraph


def face_in_x(sweep_dependency_graph: SweepDependencyGraph) -> int:
    """Incoming-face offset along axis 0 (0 for ``sx >= 0``, else 1)."""
    return sweep_dependency_graph.face_in[0]


def face_out_x(sweep_dependency_graph: SweepDependencyGraph) -> int:
    """Outgoing-face offset along axis 0 (1 for ``sx >= 0``, else 0)."""
    return sweep_dependency_graph.face_out[0]


def face_in_y(sweep_dependency_graph: SweepDependencyGraph) -> int:
    """Incoming-face offset along axis 1 (0 for ``sy >= 0``, else 1)."""
    return sweep_dependency_graph.face_in[1]


def face_out_y(sweep_dependency_graph: SweepDependencyGraph) -> int:
    """Outgoing-face offset along axis 1 (1 for ``sy >= 0``, else 0)."""
    return sweep_dependency_graph.face_out[1]
