r"""Test-side kernel of the Morel-Montry half-angle grid (moved off production).

``compute_psi_half_per_level`` of :class:`orpheus.sn.angular.closure._MMHalfGrid`'s
module had no production caller (the matvec reads ``upstream_per_ordinate`` and the
kernel ``_psi_half_grid_single_level`` directly); it is the hand-built-coefficient
verification surface of the M-M recurrence and lives here, in the gate directory of
the package it left. The two slices that moved with it on 2026-11-20
(``upstream`` and ``downstream_per_ordinate``) had no consumer once their own
indexing tests went, and were retired on 2026-11-21.
"""
from __future__ import annotations

import numpy as np

from orpheus.sn.angular.closure import _MMHalfGrid, _psi_half_grid_single_level




def compute_psi_half_per_level(
    psi_level: np.ndarray,
    tau_level: np.ndarray,
    *,
    psi_half_seed: np.ndarray | None = None,
) -> _MMHalfGrid:
    r"""Return the half-angle grid :math:`\phi_{m\pm 1/2, i, g}` for one level
    under the M-M recurrence, wrapped in the private :class:`_MMHalfGrid` accessor.

    ``tau_level`` arrives via argument, so algebraic-identity tests verify the
    recurrence under hand-built :math:`\tau` (no closure instance, no mesh).
    Production runs the SAME :func:`_psi_half_grid_single_level` kernel through
    ``MorelMontryAngularSweep.precompute_psi_state``.

    Parameters
    ----------
    psi_level :
        Shape ``(ng, M_p, nx)`` — the cell-centre angular flux at the ``M_p``
        ordinates of one level.
    tau_level :
        Shape ``(M_p,)``: the Morel-Montry angular closure weights for the level.
    psi_half_seed :
        The half-angle face flux seed VALUES, shape ``(ng, nx)``; ``None`` seeds
        the recurrence at zero.

    Returns
    -------
    _MMHalfGrid
        Accessor wrapping the half-angle grid ``faces`` of shape ``(ng, M_p+1, nx)``.
    """
    faces = _psi_half_grid_single_level(
        psi_level, tau_level, psi_half_seed=psi_half_seed,
    )
    return _MMHalfGrid(faces=faces)
