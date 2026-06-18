r"""Physics-free spatial-moment layout policy (tensor-Legendre Kronecker).

L1 primitive (mathematics, knows no neutrons) — the moment-axis sibling of
:mod:`orpheus.numerics.face_layout`. The two layout conventions every
spatial-moment consumer keys on — the slot-0 cell/face **average** index
and the "append a trailing moment axis iff there is more than one moment"
**tail** policy — live HERE, in exactly one place.

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

This module is leaf (no imports beyond the standard library) — pure
``int`` / ``tuple`` arithmetic — so importing it can never re-introduce
a cycle.
"""

from __future__ import annotations


__all__ = ["AVERAGE_MOMENT", "face_moment_tail"]


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
