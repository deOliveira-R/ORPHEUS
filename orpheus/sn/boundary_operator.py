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

    def _apply_faces(
        self, psi: "TimedFullField", method: str,
    ) -> "TimedFullField":
        r"""Apply each face's law (``method`` ∈ {apply, apply_transpose}) to its
        slot; emit a zero-bulk, boundary-only :class:`TimedFullField`."""
        from orpheus.transport.fields.boundary_flux import BoundaryFlux
        from orpheus.transport.source_sinks import AngularSourceSink
        from orpheus.transport.timed_full_field import TimedFullField

        mesh = psi.bulk.mesh
        out_boundary = BoundaryFlux.zeros_on(mesh)
        for face, law in self._face_laws.items():
            face_in = psi.boundary.face_view(face)
            out_boundary.face_view(face)[:] = getattr(law, method)(face_in)
        return TimedFullField(
            bulk=AngularSourceSink.from_mesh(
                np.zeros_like(psi.bulk.values), mesh,
            ),
            boundary=out_boundary,
            _history=(),
            history_depth=psi.history_depth,
        )

    def apply(self, psi: "TimedFullField") -> "TimedFullField":
        r"""Forward action ``B·ψ`` — per-face boundary law on the trace, zero bulk."""
        return self._apply_faces(psi, "apply")

    def apply_transpose(self, psi: "TimedFullField") -> "TimedFullField":
        r"""Euclidean transpose ``Bᵀ·ψ`` — per-face ``apply_transpose``, zero bulk.

        Reachable only when every per-face law advertises
        :data:`CAP_APPLY_TRANSPOSE` (see :attr:`capabilities`). The white BC has
        no Euclidean transpose; its physically-correct adjoint is ``B.H`` under
        the ``|Ω·n|·w`` trace metric (Wave O step O.2).
        """
        return self._apply_faces(psi, "apply_transpose")
