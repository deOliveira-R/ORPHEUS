r"""SN method-space carrying mesh + quadrature + trace metadata.

The :class:`SNMethodSpace` is the SN-specific argument to
:meth:`~orpheus.sn.boundary_realizer.SNBoundaryRealizer.realize`. It
carries everything an SN realizer needs to turn a
:class:`~orpheus.geometry.boundary.BoundaryTraceLaw` into a 1-arg
:class:`~orpheus.numerics.operator.LinearOperator`:

* :attr:`mesh` -- the spatial mesh (Mesh1D / Mesh2D). Carries
  face-normal information.
* :attr:`quadrature` -- the angular quadrature. Carries ordinate
  direction cosines + weights.
* :attr:`face` -- the face label (``"left"``, ``"right"``,
  ``"xmin"``, ``"xmax"``, ``"ymin"``, ``"ymax"``). Identifies WHICH
  face the realized BC operates on.
* :attr:`inflow_trace` -- :class:`~orpheus.numerics.spaces.trace_space.InflowTraceSpace`
  precomputed for the mesh+quad pair. Used by the realizer's vacuum
  branch to extract per-face inflow indices.
* :attr:`outflow_trace` -- :class:`~orpheus.numerics.spaces.trace_space.OutflowTraceSpace`
  precomputed for the mesh+quad pair. Symmetric for any future law
  that needs outflow indices.

The class also exposes :meth:`inflow_indices_for_face` (delegating to
the held :attr:`inflow_trace`) so the realizer can request indices
without knowing about the trace-space type.

Backward-compatible construction
================================

Wave 5 introduced a minimal :meth:`SNMethodSpace.minimal` factory
that returns a method space carrying only the quadrature (with
``face=None``, ``inflow_indices=None``, ``mesh=None``,
``inflow_trace=None``, ``outflow_trace=None``). Wave 8 keeps that
factory functional: realizers / tests that don't need mesh+face
metadata can still call it. Vacuum realization (which requires
inflow_indices) will fail loudly on a minimal method space -- that's
the right error.

The Wave-5 ``SNMethodSpace`` lived at
:mod:`orpheus.sn.boundary_realizer`; Wave 8 moves it to this
dedicated module and re-exports the same name from
:mod:`~orpheus.sn.boundary_realizer` for backward compat.

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
    from orpheus.geometry.mesh import Mesh1D, Mesh2D
    from orpheus.numerics.spaces.trace_space import InflowTraceSpace, OutflowTraceSpace
    from orpheus.numerics.quadrature import Quadrature


__all__ = ["SNMethodSpace"]


@dataclass(frozen=True)
class SNMethodSpace:
    r"""Full SN method space carrying mesh + quadrature + trace metadata.

    See module docstring for the per-field semantics.

    Wave 5 (the minimal placeholder this supersedes) only required
    ``quadrature`` plus optional ``face`` / ``inflow_indices``. Those
    three fields stay first in the field ordering and keep their
    defaults so existing callers
    ``SNMethodSpace(quadrature=q, face=..., inflow_indices=...)``
    keep working unchanged.

    Wave 8 adds ``mesh``, ``inflow_trace``, ``outflow_trace`` -- all
    with ``None`` defaults so legacy construction paths are not
    forced to populate them. The canonical construction site is
    :meth:`for_face` (used by :meth:`SNMesh._resolve_bcs`) which
    derives ``inflow_indices`` from ``inflow_trace`` and keeps the
    references for any future consumer.
    """

    quadrature: "AngularQuadrature"
    face: Optional[str] = None
    inflow_indices: Optional[np.ndarray] = None
    mesh: "Optional[Mesh1D | Mesh2D]" = None
    inflow_trace: "Optional[InflowTraceSpace]" = None
    outflow_trace: "Optional[OutflowTraceSpace]" = None

    @classmethod
    def minimal(cls, quadrature: "AngularQuadrature") -> "SNMethodSpace":
        """Quadrature-only method space (Wave 5 backward-compat).

        Used by realizers / tests that don't need mesh+face
        metadata -- typically the Wave 5 unit tests and the Wave 6
        snapshot harness's realized-1-arg path. Vacuum realization
        on such a space raises :class:`BoundaryError` (since the
        per-face inflow indices are required) -- that's the
        intended behaviour.
        """
        return cls(quadrature=quadrature)

    @classmethod
    def for_face(
        cls,
        *,
        mesh: "Mesh1D | Mesh2D",
        quadrature: "AngularQuadrature",
        face: str,
        inflow_trace: "Optional[InflowTraceSpace]" = None,
        outflow_trace: "Optional[OutflowTraceSpace]" = None,
    ) -> "SNMethodSpace":
        r"""Build a method space for a specific face.

        If ``inflow_trace`` is provided, extract the inflow indices
        for ``face`` and store them. This is the standard
        construction site at :meth:`SNMesh._resolve_bcs` time --
        the trace space is built once per ``(mesh, quad)`` pair and
        passed through so per-face indices are derived on demand.

        Parameters
        ----------
        mesh
            Spatial mesh.
        quadrature
            Angular quadrature.
        face
            Face label (one of ``"left"``, ``"right"``, ``"xmin"``,
            ``"xmax"``, ``"ymin"``, ``"ymax"``).
        inflow_trace, outflow_trace
            Optional precomputed trace spaces. If ``inflow_trace``
            is non-``None`` the ``inflow_indices`` field is
            populated for this face.
        """
        inflow_indices: Optional[np.ndarray] = None
        if inflow_trace is not None:
            inflow_indices = inflow_trace.inflow_indices_for_face(face)
        return cls(
            quadrature=quadrature,
            face=face,
            inflow_indices=inflow_indices,
            mesh=mesh,
            inflow_trace=inflow_trace,
            outflow_trace=outflow_trace,
        )

    def inflow_indices_for_face(self, face: str) -> np.ndarray:
        """Delegate to the held :class:`InflowTraceSpace`.

        Raises :class:`RuntimeError` if no ``inflow_trace`` was
        attached -- in that case the method space cannot serve
        per-face inflow lookup.
        """
        if self.inflow_trace is None:
            raise RuntimeError(
                f"SNMethodSpace.inflow_indices_for_face({face!r}): no "
                f"InflowTraceSpace attached to this method space. "
                f"Construct via .for_face(...) with inflow_trace= "
                f"populated, or use the realizer's vacuum dispatch "
                f"branch which reads .inflow_indices directly."
            )
        return self.inflow_trace.inflow_indices_for_face(face)
