r"""Physics-free spatial-moment layout policy (tensor-Legendre Kronecker).

L1 primitive (mathematics, knows no neutrons) — the moment-axis sibling of
:mod:`orpheus.numerics.face_layout`. The physics-free moment-axis primitives
every spatial-moment consumer keys on live HERE, in exactly one place: the
slot-0 cell/face **average** index (:data:`AVERAGE_MOMENT`), the "append a
trailing moment axis iff there is more than one moment" **tail** policy
(:func:`face_moment_tail`), and the rank-based "is this buffer moment-valued?"
discriminator (:func:`is_moment_valued_by_rank`).

Why ``numerics`` and not ``sn.spatial`` (#245)
==============================================

These conventions describe the Kronecker ordering of a tensor-Legendre
DG basis; they carry no transport physics (no :math:`\Sigma`, no
:math:`\mu`). Both halves of the spatial-moment iterate need them:

* the typed :class:`~orpheus.numerics.spaces.spatial_moment_space.SpatialMomentSpace`
  (the CAPABILITY half — surfaces ``average_moment_index`` and the cell
  tail ``spatial_moment_tail``), and
* the SN-side UBLD cell assembler :mod:`orpheus.sn.spatial._ubld` (the
  REALIZATION half — buffers, sweeps, the face cochain).

``numerics`` sits *below* ``sn``, so homing the policy here lets both
import it in the correct (downward) direction. Previously the constants
lived in ``sn.spatial._ubld`` and ``SpatialMomentSpace`` reached UP into
``sn.spatial`` for them via a *deferred* (call-time) import — a band-aid
over a layering inversion (a ``numerics`` space depending on the SN
package). Relocating here removes the deferral and the inversion; the
SN module re-exports these names downward (the
:data:`orpheus.numerics.face_layout.AXIS_NAMES` precedent) so SN
consumers keep importing them next to the UBLD primitives they name.

This module is leaf — no ``orpheus`` imports (numpy enters only under
``TYPE_CHECKING``, for the rank-predicate type hints) — so importing it can
never re-introduce a cycle.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import numpy as np


__all__ = ["AVERAGE_MOMENT", "face_moment_tail", "is_moment_valued_by_rank"]


#: Index of the cell/face AVERAGE moment in the tensor-Legendre Kronecker
#: layout (``[bar, …]``): the all-``P₀`` moment is first (d=2 cell order
#: ``[ψ̄, ψ̂_y, ψ̂_x, ψ̂_xy]``; per-axis face order ``[bar, slope]``).  Single
#: source for the slot-0 convention every moment consumer reduces on (#240 D5b)
#: — change the layout here, not at the scattered ``[..., 0]`` call sites.
AVERAGE_MOMENT = 0


def face_moment_tail(n_face_moments: int) -> tuple[int, ...]:
    r"""Trailing moment-axis shape suffix for a face-cochain buffer.

    A multi-moment closure (LD's bilinear face, ``n_face_moments > 1``) carries a
    trailing ``2^{d-1}``-moment axis; a cell-average closure (DD/Step,
    ``== 1``) leaves the face rank untouched (``()`` — NO length-1 axis appended)
    so its buffers stay byte-identical (#240 D5b — the backward-compat invariant).
    Single source for both storage policies (``_MovingFrontier`` window +
    ``FullFieldWavefront`` full cochain), which must agree on the tail shape.
    """
    return () if n_face_moments == 1 else (n_face_moments,)


def is_moment_valued_by_rank(array: "np.ndarray", reference: "np.ndarray") -> bool:
    r"""Does ``array`` carry the trailing spatial-moment axis, judged by RANK?

    ``True`` iff ``array`` has more than one axis beyond ``reference`` — the
    S4-safe discriminator for "is this a moment-valued buffer". A moment buffer
    is ``(N…, ng, *spatial, 2^d)`` while its scalar reference (``Σ_t`` /
    ``reaction_xs``) is the per-ordinate-stripped ``(ng, *spatial)``, so a genuine
    moment buffer carries one MORE leading (ordinate) axis PLUS the trailing
    ``2^d`` moment axis — net ``> reference.ndim + 1`` — whereas a flat
    ``(N…, ng, *spatial)`` buffer (a matvec-zero / flat external source) sits at
    exactly ``reference.ndim + 1``.

    RANK, not trailing-size: a coincidental ``n_diag == 2^d`` (a d=2 anti-diagonal
    of exactly 4 cells) mis-fires a ``shape[-1] == 2^d`` probe, but the rank of a
    flat buffer never collides with a moment buffer's (#246). The single source
    for the matvec moment-broadcast
    (:func:`orpheus.sn.loss_representation._moment_broadcast_sigma`) and the
    cell-solve source-reframe gate (``orpheus.sn.sweep_graph._CellSolve``).
    """
    return array.ndim > reference.ndim + 1
