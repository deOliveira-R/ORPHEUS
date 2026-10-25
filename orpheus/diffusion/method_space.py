r"""Diffusion method-space carrying scalar-trace metadata.

The :class:`DiffusionMethodSpace` is the diffusion-specific argument to
:data:`~orpheus.diffusion.boundary_realizer.DIFFUSION_REALIZERS`
— the sibling of :class:`~orpheus.sn.mesh.method_space.SNMethodSpace`
(same ``trace`` field name — the per-method protocol surface kept at
#290 P2.5). SN's space is ``(quadrature, trace)`` with no face (the
face travels in the ``FaceSite`` beside it, WF2.4); this one still
carries a face label that no realizer reads (the same two-encodings
exposure; #661 holds the collapse trigger). What it
carries differs because the method does:

* :attr:`~orpheus.diffusion.method_space.DiffusionMethodSpace.face` -- the face label (``"xmin"`` / ``"xmax"`` / …).
  Identifies WHICH face the realized BC operates on.
* :attr:`trace` -- the
  :class:`~orpheus.numerics.spaces.scalar_trace_space.ScalarTraceSpace`
  (``DiffusionMesh.scalar_trace``, passed to :meth:`for_face` by the
  mesh, #290 P7a). Carries the face inventory, the
  ``(2, ng, *face_spatial)`` slot shapes, and the
  :attr:`~orpheus.numerics.spaces.scalar_trace_space.ScalarTraceSpace.OUTFLOW_ROW`
  / ``INFLOW_ROW`` component convention the P4 boundary operator
  assembles against.

There is deliberately NO quadrature and NO ``inflow_indices`` field:
the SN trace's inflow/outflow split is a quadrature-dependent sign
partition (read off its ladder), while the scalar trace's
inflow/outflow split is STRUCTURAL — the fixed ``(J⁺, J⁻)`` component
rows owned by :class:`ScalarTraceSpace`. The diffusion realizer's
albedo-family scalars read nothing from the method space at all; the fields exist for the P4
``DiffusionBoundaryOperator`` assembly and the uniform cross-method
construction surface (the ``TransportMethod`` Protocol, #290 P7).

There is no ``mesh`` field either. It was metadata nothing read (0
production reads), and it was the one import of
:class:`~orpheus.diffusion.augmented_mesh.DiffusionMesh` here, so a
method space below the mesh that builds it named the mesh back: the
diffusion package's only module cycle (layer DAG, 2026-10). The method
space is a statement about one face of a trace, and the trace is all it
needs.

References
----------

* ``.claude/plans/archive/diffusion_integration_290.md`` §P3 (this module's
  brief) + ``.claude/plans/diffusion_crosswalk.md`` (the ``(J⁺, J⁻)``
  convention contract).
* :mod:`orpheus.sn.mesh.method_space` — the mirrored template.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Optional

if TYPE_CHECKING:
    from orpheus.numerics.spaces.scalar_trace_space import ScalarTraceSpace


__all__ = ["DiffusionMethodSpace"]


@dataclass(frozen=True)
class DiffusionMethodSpace:
    r"""Diffusion method space carrying scalar-trace metadata.

    See the module docstring for the per-field semantics and for why
    there is no quadrature / inflow-index field (the scalar trace's
    inflow/outflow split is structural, not quadrature-dependent).

    All fields default to ``None`` so the realizer's rank-1
    albedo-family dispatch — which reads nothing from the method
    space — works on a :meth:`minimal` instance. The canonical
    construction site is :meth:`for_face` (the
    ``DiffusionMesh.realize_boundary_role`` wiring — the diffusion arm
    of the #290 P7b ``TransportMethod`` hook), which validates the face
    against the trace's inventory.
    """

    face: Optional[str] = None
    trace: "Optional[ScalarTraceSpace]" = None

    @classmethod
    def minimal(cls) -> "DiffusionMethodSpace":
        """Metadata-free method space.

        Sufficient for realizing any single albedo-family law — the
        realized operator is a face-slot-shape-agnostic scalar multiple
        of the identity, so no trace metadata is consulted.
        """
        return cls()

    @classmethod
    def for_face(
        cls,
        *,
        trace: "ScalarTraceSpace",
        face: str,
    ) -> "DiffusionMethodSpace":
        r"""Build the method space for one boundary face of ``trace``.

        The standard construction site — one method space per boundary
        face (each method-mesh's
        ``realize_boundary_role`` arm builds its method space here,
        driven by the shared
        :func:`~orpheus.transport.method.realize_boundary`
        body, #290 P7b). The mesh passes its own ``scalar_trace``, the
        single source. The scalar trace has nothing to derive per face
        (the ``(J⁺, J⁻)`` rows are fixed); the face is VALIDATED against
        the trace's face inventory so a method space for a nonexistent
        face is unrepresentable.

        Parameters
        ----------
        trace
            The mesh's
            :class:`~orpheus.numerics.spaces.scalar_trace_space.ScalarTraceSpace`,
            which supplies the face inventory.
        face
            Face name (``"{axis}{min|max}"`` — e.g. ``"xmin"``; a
            radial outer face renders ``"xmax"``; the pole is not a
            face).

        Raises
        ------
        ValueError
            If ``face`` is not in the trace's face inventory.
        """
        if face not in trace.face_names:
            raise ValueError(
                f"DiffusionMethodSpace.for_face: face {face!r} is not a "
                f"boundary face of the mesh's ScalarTraceSpace; "
                f"available: {trace.face_names}."
            )
        return cls(face=face, trace=trace)
