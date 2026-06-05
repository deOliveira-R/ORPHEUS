r"""The realized boundary law ``B`` as a whole-trace BOUNDARY-block operator.

Wave O (Issue #208) extracts the boundary conditions from the streaming
operator ``L`` so that ``(L_full + C − S − F − B)\psi = q`` is the canonical
transport algebra. ``B`` is the **boundary-block** operator on the direct-sum
transport state ``V = V_bulk ⊕ V_boundary``: a 2×2 block matrix with only the
``A_ss`` (boundary → boundary) block non-zero. It maps the **outflow** trace to
the **inflow** trace via the per-face realized boundary laws (reflective /
vacuum / white / albedo / periodic), with **no bulk action**.

Block structure
===============

On ``V = V_bulk ⊕ V_boundary`` the four operator families are::

    C, S, F  →  [ A_bb  0 ]   (BULK   — bulk → bulk only)
                [ 0     0 ]

    L_full   →  [ A_bb  A_bs ] (FULL   — streaming couples bulk ↔ boundary)
                [ A_sb  0    ]

    B        →  [ 0     0   ] (BOUNDARY — boundary → boundary only, ``A_ss``)
                [ 0     A_ss]

:class:`SNBoundaryOperator` is the ``A_ss`` leaf. As a sibling ``−B`` of ``L``
in the :class:`~orpheus.numerics.operator.OperatorSum` algebra it supplies the
reflective coupling that ``L`` previously absorbed inside its sweep (the
``inflow = bc.apply(outflow)`` re-apply); the outer Krylov / SI loop then drives
the boundary **consistency residual** ``ψ.inflow − B·ψ.outflow − q.inflow → 0``.

Construction
============

The per-face boundary laws already live on the
:class:`~orpheus.sn.geometry.SNMesh` as ``bc_xmin`` / ``bc_xmax`` / ``bc_ymin``
/ ``bc_ymax`` (each a :class:`~orpheus.geometry.boundary._bound_compat._BoundBoundaryOperator`
wrapping a realized law that carries :attr:`BlockRole.BOUNDARY`). The whole-trace
``B`` is the block-diagonal composition over the mesh's true boundary faces: for
each face present in the trace it applies that face's law to that face's slot.
``B`` is therefore block-diagonal over faces — it never mixes faces.

See :ref:`operator-algebra` and the Wave O plan
``.claude/plans/wave_o_operator_typing.md`` (step O.4a.2).
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Optional

import numpy as np

from orpheus.numerics.operator import (
    CAP_APPLY,
    CAP_APPLY_TRANSPOSE,
    BlockRole,
    LinearOperatorMixin,
)

if TYPE_CHECKING:
    from orpheus.numerics.space import FunctionSpace
    from orpheus.sn.geometry import SNMesh
    from orpheus.transport.fields.boundary_flux import BoundaryFlux
    from orpheus.transport.source_sinks import BoundarySourceSink
    from orpheus.transport.timed_full_field import TimedFullField


__all__ = ["SNBoundaryOperator"]


class SNBoundaryOperator(LinearOperatorMixin):
    r"""Whole-trace boundary law ``B`` — the ``A_ss`` block of the SN algebra.

    Block-diagonal over the mesh's true boundary faces: ``B.apply(ψ)`` returns a
    :class:`~orpheus.transport.timed_full_field.TimedFullField` with **zero bulk**
    and, on each face, ``bc_<face>.apply(ψ.boundary.face_view(<face>))`` — the
    per-face realized boundary law applied to that face's trace slot. It composes
    as ``−B`` in ``(L_full + C − S − F − B)`` (acting on the same
    :class:`TimedFullField` carrier as ``L``/``C``/``S``/``F``).

    The role is :attr:`BlockRole.BOUNDARY`; the domain and codomain are the
    mesh's unified :class:`~orpheus.numerics.spaces.trace_space.TraceSpace`
    (``sn_mesh.trace``) — the cosine-weighted ``|Ω·n|·w`` boundary metric on that
    space (Wave O step O.2) makes the Hilbert adjoint ``B.H`` the physically
    correct partial-current adjoint, which is the one channel by which the
    white-BC adjoint becomes available.

    Capabilities
    ------------

    ``apply`` always. ``apply_transpose`` is advertised iff EVERY per-face law
    advertises it — reflective (involution), vacuum, periodic and albedo do;
    **white does NOT** (it is self-adjoint only under the ``|Ω·n|·w`` metric, so
    its adjoint routes through ``B.H`` on the weighted trace space at O.2). The
    intersection rule keeps ``apply_transpose`` honest: it is reachable only when
    every block can honour it.

    Parameters
    ----------
    sn_mesh : SNMesh
        The augmented geometry — carries the per-face boundary laws
        (``bc_xmin``/...) and the unified trace space (same instance the
        composite carrier is bound to; the mesh-identity invariant of
        :class:`~orpheus.sn.operator.StreamingOperator` applies here too).
    """

    block_role = BlockRole.BOUNDARY

    def __init__(self, sn_mesh: "SNMesh") -> None:
        self.sn_mesh = sn_mesh

    @property
    def _face_laws(self) -> dict[str, object]:
        """Map each true boundary face → its per-face realized law.

        Read from ``sn_mesh.bc_<face>`` for the faces the trace carries
        (slab ``xmin``/``xmax``; curvilinear ``xmax`` only; 2-D Cartesian
        all four). Single source of truth — the laws are the same objects
        the sweep consumes, so ``B`` cannot drift from the realized BCs.
        """
        return {
            face: getattr(self.sn_mesh, f"bc_{face}")
            for face in self.sn_mesh.trace.layout.faces
        }

    @property
    def capabilities(self) -> frozenset[str]:
        caps = {CAP_APPLY}
        laws = self._face_laws.values()
        if laws and all(
            CAP_APPLY_TRANSPOSE in law.capabilities for law in laws
        ):
            caps.add(CAP_APPLY_TRANSPOSE)
        return frozenset(caps)

    @property
    def domain(self) -> Optional["FunctionSpace"]:
        return self.sn_mesh.trace

    @property
    def codomain(self) -> Optional["FunctionSpace"]:
        return self.sn_mesh.trace

    def _reflect_trace(
        self, boundary: "BoundaryFlux", method: str,
    ) -> "BoundarySourceSink":
        r"""Core ``A_ss`` action on the trace ALONE — apply each face's law
        (``method`` ∈ {apply, apply_transpose}) to that face's slot, project onto
        the codomain row, and return a boundary-only
        :class:`~orpheus.transport.source_sinks.BoundarySourceSink`.

        ``B`` is the ``A_ss`` block ``V_outflow → V_inflow``: it maps the
        **outflow** trace to the **inflow** trace, so the forward action must be
        non-zero **only on the inflow ordinate slots** of each face. The realized
        per-face law (a :class:`~orpheus.numerics.operator.PermutationOperator`
        for reflective, :class:`AngularAverageOperator` for white, …) is a
        *full-face* operator: e.g. the specular permutation also maps the input
        inflow slots onto the output **outflow** slots (``R·ψ.inflow``). The
        legacy sweep only ever read the inflow slots of ``bc.apply(face)`` so
        that spurious outflow output was harmless — but ``−B`` as a sibling of
        ``L`` in ``(L_full + C − S − F − B)`` reads the WHOLE boundary block, and
        a non-zero outflow emission would corrupt the outflow-definition residual
        ``ψ.outflow − streamed`` (it is supposed to carry no ``B`` term — see the
        block matrix in :mod:`orpheus.sn.boundary_operator`). So the forward
        action is projected onto ``inflow_indices_for_face`` (the consistency
        row); the Euclidean transpose ``Bᵀ: V_inflow → V_outflow`` is projected
        onto ``outflow_indices_for_face`` accordingly. (The metric-correct Hilbert
        adjoint ``B.H`` under ``|Ω·n|·w`` is Wave O step O.2; this Euclidean
        ``apply_transpose`` is the un-weighted shadow, not yet a live consumer.)

        This is the **single source of truth** for the boundary reflection: both
        the full-field :meth:`apply` (lifted onto a zero-bulk carrier) and the
        trace-only :meth:`reflect_into_inflow` (the direct-loop inflow seed) route
        through it, so the two cannot drift (Cardinal Rule 2).
        """
        from orpheus.transport.source_sinks import BoundarySourceSink

        # Single mesh source (mesh-identity invariant — see class docstring):
        # the output buffers, the trace selectors, and ``_face_laws`` ALL read
        # ``self.sn_mesh``, so a mismatched input trace cannot desync the
        # projection from the buffer geometry.
        mesh = self.sn_mesh
        trace = mesh.trace
        out_boundary = BoundarySourceSink.zeros_on(mesh)
        for face, law in self._face_laws.items():
            face_in = boundary.face_view(face)
            full = getattr(law, method)(face_in)
            target = (
                trace.inflow_indices_for_face(face)
                if method == "apply"
                else trace.outflow_indices_for_face(face)
            )
            out_boundary.face_view(face)[target] = full[target]
        return out_boundary

    def _apply_faces(
        self, psi: "TimedFullField", method: str,
    ) -> "TimedFullField":
        r"""Lift the trace-only :meth:`_reflect_trace` onto the full
        :class:`TimedFullField` carrier with **zero bulk** — the ``A_ss`` block
        as an operator on ``V = V_bulk ⊕ V_boundary``.
        """
        from orpheus.transport.source_sinks import AngularSourceSink
        from orpheus.transport.timed_full_field import TimedFullField

        mesh = self.sn_mesh
        if psi.bulk.mesh is not mesh:
            raise ValueError(
                "SNBoundaryOperator.apply: input field and operator must "
                "share the same SNMesh instance (mesh-identity invariant); "
                f"got field mesh {psi.bulk.mesh!r} vs operator mesh {mesh!r}."
            )
        return TimedFullField(
            bulk=AngularSourceSink.from_mesh(
                np.zeros_like(psi.bulk.values), mesh,
            ),
            boundary=self._reflect_trace(psi.boundary, method),
            _history=(),
            history_depth=psi.history_depth,
        )

    def apply(self, psi: "TimedFullField") -> "TimedFullField":
        r"""Forward action ``B·ψ`` — per-face boundary law on the trace, zero bulk."""
        return self._apply_faces(psi, "apply")

    def reflect_into_inflow(
        self, boundary: "BoundaryFlux",
    ) -> "BoundarySourceSink":
        r"""Trace-only forward reflection ``B·ψ.outflow`` projected onto the
        inflow row — the ``A_ss`` action expressed on the boundary trace ALONE.

        Returns a boundary-only
        :class:`~orpheus.transport.source_sinks.BoundarySourceSink` whose
        **inflow** ordinate slots carry the per-face reflected outflow (``R·G``
        for reflective, the angular average for white, zero for vacuum) and whose
        outflow slots are zero. It is :meth:`apply` without the zero-bulk carrier
        — the entry the direct fixed-source SI loop and the final eigenvalue
        reconstruction sweep use to seed ``ψ.inflow = B·ψ.outflow`` on a bare
        boundary buffer, without fabricating a throwaway zero-bulk field just to
        reach the ``A_ss`` block.
        """
        return self._reflect_trace(boundary, "apply")

    def apply_transpose(self, psi: "TimedFullField") -> "TimedFullField":
        r"""Euclidean transpose ``Bᵀ·ψ`` — per-face ``apply_transpose``, zero bulk.

        Reachable only when every per-face law advertises
        :data:`CAP_APPLY_TRANSPOSE` (see :attr:`capabilities`). The white BC has
        no Euclidean transpose; its physically-correct adjoint is ``B.H`` under
        the ``|Ω·n|·w`` trace metric (Wave O step O.2).
        """
        return self._apply_faces(psi, "apply_transpose")
