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
* :attr:`face` -- the face label (``"xmin"``, ``"xmax"``,
  ``"ymin"``, ``"ymax"``). Identifies WHICH face the realized BC
  operates on. (1-D meshes use ``"xmin"`` / ``"xmax"``; the radial
  axis IS the x-axis. A solid sphere / cylinder has only the outer
  ``"xmax"`` face — the pole at r=0 is the angular closure's
  regularity condition, not a BC face.)
* :attr:`trace` -- the unified
  :class:`~orpheus.numerics.spaces.trace_space.TraceSpace` precomputed
  for the mesh+quad pair. The realizer's vacuum branch extracts
  per-face inflow indices from it (and any future law needing outflow
  indices reads the same signed-:math:`\Omega\cdot\hat n` data).

The class also exposes :meth:`inflow_indices_for_face` (delegating to
the held :attr:`trace`) so the realizer can request indices without
knowing about the trace-space type.

Backward-compatible construction
================================

Wave 5 introduced a minimal :meth:`SNMethodSpace.minimal` factory
that returns a method space carrying only the quadrature (with
``face=None``, ``inflow_indices=None``, ``mesh=None``,
``trace=None``). Wave 8 keeps that factory functional: realizers /
tests that don't need mesh+face metadata can still call it. Vacuum
realization (which requires inflow_indices) will fail loudly on a
minimal method space -- that's the right error.

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
    from orpheus.numerics.spaces.trace_space import TraceSpace
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

    Wave 8 adds ``mesh`` + ``trace`` -- both with ``None`` defaults so
    legacy construction paths are not forced to populate them. The
    canonical construction site is :meth:`for_face` (used by
    :meth:`SNMesh._resolve_bcs`) which derives ``inflow_indices`` from
    ``trace`` and keeps the reference for any future consumer.
    """

    quadrature: "AngularQuadrature"
    face: Optional[str] = None
    inflow_indices: Optional[np.ndarray] = None
    mesh: "Optional[Mesh1D | Mesh2D]" = None
    trace: "Optional[TraceSpace]" = None

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
        mesh: "Optional[Mesh1D | Mesh2D]" = None,
        quadrature: "AngularQuadrature",
        face: str,
        trace: "Optional[TraceSpace]" = None,
    ) -> "SNMethodSpace":
        r"""Build a method space for a specific face.

        If ``trace`` is provided, extract the inflow indices for
        ``face`` and store them. This is the standard construction
        site at :meth:`SNMesh._resolve_bcs` time -- the trace space is
        built once per quadrature+layout pair and passed through so
        per-face indices are derived on demand.

        Parameters
        ----------
        mesh
            Spatial mesh — OPTIONAL metadata (C5.3, #225): nothing in
            the realizer chain reads it (inflow indices come from the
            trace); an axis-native ``SNMesh`` with no legacy mesh
            adapter passes ``None``.
        quadrature
            Angular quadrature.
        face
            Face name (``"{axis}{min|max}"`` — e.g. ``"xmin"`` …
            ``"zmax"``).
        trace
            Optional precomputed unified
            :class:`~orpheus.numerics.spaces.trace_space.TraceSpace`.
            If non-``None`` the ``inflow_indices`` field is populated
            for this face.
        """
        inflow_indices: Optional[np.ndarray] = None
        if trace is not None:
            inflow_indices = trace.inflow_indices_for_face(face)
        return cls(
            quadrature=quadrature,
            face=face,
            inflow_indices=inflow_indices,
            mesh=mesh,
            trace=trace,
        )

    def inflow_indices_for_face(self, face: str) -> np.ndarray:
        """Delegate to the held :class:`TraceSpace`.

        Raises :class:`RuntimeError` if no ``trace`` was attached --
        in that case the method space cannot serve per-face inflow
        lookup.
        """
        if self.trace is None:
            raise RuntimeError(
                f"SNMethodSpace.inflow_indices_for_face({face!r}): no "
                f"TraceSpace attached to this method space. "
                f"Construct via .for_face(...) with trace= populated, "
                f"or use the realizer's vacuum dispatch branch which "
                f"reads .inflow_indices directly."
            )
        return self.trace.inflow_indices_for_face(face)
