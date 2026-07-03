r"""Diffusion method-space carrying mesh + scalar-trace metadata.

The :class:`DiffusionMethodSpace` is the diffusion-specific argument to
:meth:`~orpheus.diffusion.boundary_realizer.DiffusionBoundaryRealizer.realize`
— the mirror of :class:`~orpheus.sn.mesh.method_space.SNMethodSpace`
(same construction surface: :meth:`minimal` / :meth:`for_face`; same
``trace`` field name — the per-method protocol surface kept at #290
P2.5). What it carries differs because the method does:

* :attr:`mesh` -- the :class:`~orpheus.diffusion.augmented_mesh.DiffusionMesh`
  the solve lives on. Metadata, mirroring SN's C5.3 ruling: nothing in
  the realizer chain reads it.
* :attr:`face` -- the face label (``"xmin"`` / ``"xmax"`` / …).
  Identifies WHICH face the realized BC operates on.
* :attr:`trace` -- the
  :class:`~orpheus.numerics.spaces.scalar_trace_space.ScalarTraceSpace`
  (``DiffusionMesh.scalar_trace`` — :meth:`for_face` reads it off the
  mesh, #290 P7a). Carries the face inventory, the
  ``(2, ng, *face_spatial)`` slot shapes, and the
  :attr:`~orpheus.numerics.spaces.scalar_trace_space.ScalarTraceSpace.OUTFLOW_ROW`
  / ``INFLOW_ROW`` component convention the P4 boundary operator
  assembles against.

There is deliberately NO quadrature and NO ``inflow_indices`` field:
the SN trace needs a per-face ordinate *selector* because inflow/outflow
is a quadrature-dependent sign partition, while the scalar trace's
inflow/outflow split is STRUCTURAL — the fixed ``(J⁺, J⁻)`` component
rows owned by :class:`ScalarTraceSpace`. The diffusion realizer's
albedo-family scalars read nothing from the method space at all
(mirroring SN's albedo branch); the fields exist for the P4
``DiffusionBoundaryOperator`` assembly and the uniform cross-method
construction surface (the ``TransportMethod`` Protocol, #290 P7).

References
----------

* ``.claude/plans/diffusion_integration_290.md`` §P3 (this module's
  brief) + ``.claude/plans/diffusion_crosswalk.md`` (the ``(J⁺, J⁻)``
  convention contract).
* :mod:`orpheus.sn.mesh.method_space` — the mirrored template.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Optional

if TYPE_CHECKING:
    from orpheus.diffusion.augmented_mesh import DiffusionMesh
    from orpheus.numerics.spaces.scalar_trace_space import ScalarTraceSpace


__all__ = ["DiffusionMethodSpace"]


@dataclass(frozen=True)
class DiffusionMethodSpace:
    r"""Diffusion method space carrying mesh + scalar-trace metadata.

    See the module docstring for the per-field semantics and for why
    there is no quadrature / inflow-index field (the scalar trace's
    inflow/outflow split is structural, not quadrature-dependent).

    All fields default to ``None`` so the realizer's rank-1
    albedo-family dispatch — which reads nothing from the method
    space — works on a :meth:`minimal` instance, exactly like SN's
    albedo branch on ``SNMethodSpace.minimal(quad)``. The canonical
    construction site is :meth:`for_face` (the
    ``DiffusionMesh.realize_boundary_law`` wiring — the diffusion arm
    of the #290 P7b ``TransportMethod`` hook), which reads the trace
    off the mesh and validates the face against its inventory.
    """

    mesh: "Optional[DiffusionMesh]" = None
    face: Optional[str] = None
    trace: "Optional[ScalarTraceSpace]" = None

    @classmethod
    def minimal(cls) -> "DiffusionMethodSpace":
        """Metadata-free method space (mirror of ``SNMethodSpace.minimal``).

        Sufficient for realizing any single albedo-family law — the
        realized operator is a face-slot-shape-agnostic scalar multiple
        of the identity, so no mesh / trace metadata is consulted.
        """
        return cls()

    @classmethod
    def for_face(
        cls,
        *,
        mesh: "DiffusionMesh",
        face: str,
    ) -> "DiffusionMethodSpace":
        r"""Build the method space for one of ``mesh``'s boundary faces.

        The standard construction site — one method space per boundary
        face, mirroring ``SNMethodSpace.for_face`` (each method-mesh's
        ``realize_boundary_law`` arm builds its method space here,
        driven by the shared
        :func:`~orpheus.transport.method.resolve_boundary_conditions`
        body, #290 P7b). The trace is READ
        OFF the mesh (``mesh.scalar_trace`` — the single source; no
        separate ``trace=`` parameter to drift). Where SN's
        ``for_face`` *derives* the per-face inflow indices from its
        trace, the scalar trace has nothing to derive (the ``(J⁺,
        J⁻)`` rows are fixed); instead the face is VALIDATED against
        the trace's face inventory so a method space for a nonexistent
        face is unrepresentable.

        Parameters
        ----------
        mesh
            The :class:`~orpheus.diffusion.augmented_mesh.DiffusionMesh`
            whose trace supplies the face inventory.
        face
            Face name (``"{axis}{min|max}"`` — e.g. ``"xmin"``; a
            radial outer face renders ``"xmax"``; the pole is not a
            face).

        Raises
        ------
        ValueError
            If ``face`` is not in the mesh trace's face inventory.
        """
        if face not in mesh.scalar_trace.face_names:
            raise ValueError(
                f"DiffusionMethodSpace.for_face: face {face!r} is not a "
                f"boundary face of the mesh's ScalarTraceSpace; "
                f"available: {mesh.scalar_trace.face_names}."
            )
        return cls(mesh=mesh, face=face, trace=mesh.scalar_trace)
