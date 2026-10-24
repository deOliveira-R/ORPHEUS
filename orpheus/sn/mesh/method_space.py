r"""SN method-space carrying mesh + quadrature + trace metadata.

The :class:`SNMethodSpace` is the SN-specific argument to
:meth:`~orpheus.sn.boundary.realizer.SN_REALIZERS.realize`. It
carries everything an SN realizer needs to turn a declared unit (a chain
of the response network) into a
:class:`~orpheus.sn.boundary.realized_operators.RealizedBoundary`:

* :attr:`mesh` -- the spatial mesh (Mesh1D / Mesh2D). Carries
  face-normal information.
* :attr:`quadrature` -- the angular quadrature. Carries ordinate
  direction cosines + weights.
* :attr:`face` -- the face label (``"xmin"``, ``"xmax"``,
  ``"ymin"``, ``"ymax"``). Identifies WHICH face the realized BC
  operates on. (1-D meshes use ``"xmin"`` / ``"xmax"``; the radial
  axis IS the x-axis. A solid sphere / cylinder has only the outer
  ``"xmax"`` face — the pole at r=0 is the angular closure's
  regularity condition, not a BC face.)
* :attr:`trace` -- the unified
  :class:`~orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace` precomputed
  for the mesh+quad pair. The typed realized operators read every
  half-space of a face off it (``trace.inflow_space(face)``,
  ``trace.outflow_space(face)``); the space carries no index set of its own.

The class also exposes :meth:`inflow_indices_for_face` (delegating to
the held :attr:`trace`) so a consumer can request indices without
knowing about the trace-space type.

Backward-compatible construction
================================

Wave 5 introduced a minimal :meth:`SNMethodSpace.minimal` factory
that returns a method space carrying only the quadrature (with
``face=None``, ``mesh=None``, ``trace=None``). Wave 8 keeps that
factory functional: realizers / tests that don't need mesh+face
metadata can still call it. A realization that reads the trace ladder
fails loudly on a minimal method space -- that's the right error.

The Wave-5 ``SNMethodSpace`` lived at
:mod:`orpheus.sn.boundary.realizer`; Wave 8 moved it to this
dedicated module, and the realizer imports it from here.

References
----------

* ``.claude/plans/transient-giggling-cake.md`` Wave 8 -- C8.1 (this
  module's brief), Wave 5 -- C5.3 (the minimal placeholder that
  this module supersedes).
* Grand Report v3 §16A.4 lines 2864-2876 (``realize(law,
  method_space)`` vocabulary).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Optional

import numpy as np

if TYPE_CHECKING:
    from orpheus.mesh import Mesh1D, Mesh2D
    from orpheus.numerics.spaces.angular_trace_space import AngularTraceSpace
    from orpheus.numerics.quadrature import Quadrature


__all__ = ["SNMethodSpace"]


@dataclass(frozen=True)
class SNMethodSpace:
    r"""Full SN method space carrying mesh + quadrature + trace metadata.

    See module docstring for the per-field semantics.

    ``quadrature`` is the only required field; ``face``, ``mesh`` and ``trace`` default to ``None`` so a quadrature-only method
    space (:meth:`minimal`) stays constructible. The canonical construction site is :meth:`for_face` (used by
    ``SNProblem.realize_boundary_role``, the SN arm of the ``TransportMethod`` hook of the shared ``realize_boundary``). The
    typed realized operators read their half-spaces off ``trace``, so the space holds no index set (ERR-041's index sets and the
    cross-check that read them retired in WF2.4; the invariant is witnessed on the ladder itself,
    ``tests/gates/numerics/test_angular_trace_space.py::TestTheInflowHalfSpaceIsTheIncomingOrdinates``).
    """

    quadrature: "Quadrature"
    face: Optional[str] = None
    mesh: "Optional[Mesh1D | Mesh2D]" = None
    trace: "Optional[AngularTraceSpace]" = None

    @classmethod
    def minimal(cls, quadrature: "Quadrature") -> "SNMethodSpace":
        """Quadrature-only method space (Wave 5 backward-compat).

        Used by realizers / tests that don't need mesh+face
        metadata. A realization that reads the trace ladder on such
        a space raises :class:`BoundaryError` -- that's the
        intended behaviour.
        """
        return cls(quadrature=quadrature)

    @classmethod
    def for_face(
        cls,
        *,
        mesh: "Optional[Mesh1D | Mesh2D]" = None,
        quadrature: "Quadrature",
        face: str,
        trace: "Optional[AngularTraceSpace]" = None,
    ) -> "SNMethodSpace":
        r"""Build a method space for a specific face.

        This is the standard construction site inside ``SNProblem.realize_boundary_role`` (per face, driven by the shared
        :func:`~orpheus.transport.method.realize_boundary` body): the trace space is built once per quadrature+layout pair at
        ``SNProblem`` construction and passed through, and every half-space is read off it on demand.

        Parameters
        ----------
        mesh
            Spatial mesh, OPTIONAL metadata (C5.3, #225): nothing in the realizer chain reads it; an axis-native ``SNProblem``
            with no legacy mesh adapter passes ``None``.
        quadrature
            Angular quadrature.
        face
            Face name (``"{axis}{min|max}"``, e.g. ``"xmin"`` to ``"zmax"``).
        trace
            Optional precomputed unified
            :class:`~orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace`.
        """
        return cls(quadrature=quadrature, face=face, mesh=mesh, trace=trace)

    def inflow_indices_for_face(self, face: str) -> np.ndarray:
        """Delegate to the held :class:`AngularTraceSpace`.

        Raises :class:`RuntimeError` if no ``trace`` was attached --
        in that case the method space cannot serve per-face inflow
        lookup.
        """
        if self.trace is None:
            raise RuntimeError(
                f"SNMethodSpace.inflow_indices_for_face({face!r}): no "
                f"AngularTraceSpace attached to this method space. "
                f"Construct via .for_face(...) with trace= populated."
            )
        return self.trace.inflow_indices_for_face(face)
