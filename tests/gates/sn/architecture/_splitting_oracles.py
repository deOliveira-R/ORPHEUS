r"""The test-side law instrument of :class:`orpheus.sn.splitting.Splitting` (moved off production).

``Splitting.law_residual`` had no production caller: it is the splitting law's certificate,
read by the gates that state ``A = M - N``. It is now a free function of the splitting.
"""
from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:  # pragma: no cover - typing only
    from orpheus.numerics.coupled_system import CoupledField
    from orpheus.sn.splitting import Splitting


def law_residual(splitting: "Splitting", state: "CoupledField") -> np.ndarray:
    r"""The residual of the splitting law, ``(M - sum N_i - A) x``, flat and aligned with ``A x``.

    The value's certificate against the Problem it was minted from (exactly ``0.0`` on every
    shipped labelling: the pieces are the loss's own terms, and the Gauss-Seidel split writes
    disjoint rows). The carrier bridge lives here: on a carrying mesh ``M`` and ``N`` consume the
    coupled state; on a seedless mesh they are System-A operators and the state's single member is
    applied.
    """
    loss_image = splitting.system.loss.apply(state).to_flat()
    if splitting.system.is_coupled:
        image = splitting.implicit.apply(state).to_flat()
        for gain in splitting.explicit:
            image = image - gain.apply(state).to_flat()
        return image - loss_image
    member = state.systems[0]
    image = splitting.implicit.apply(member).to_flat()
    for gain in splitting.explicit:
        image = image - gain.apply(member).to_flat()
    return image - loss_image
