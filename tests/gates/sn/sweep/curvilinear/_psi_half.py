r"""Test-side accessors of the Morel-Montry half-angle grid (moved off production).

``compute_psi_half_per_level`` and the two slices ``downstream_per_ordinate`` and
``upstream`` of :class:`orpheus.sn.angular.closure._MMHalfGrid` had no production
caller (the matvec reads ``upstream_per_ordinate`` and the kernel
``_psi_half_grid_single_level`` directly); they are the hand-built-coefficient
verification surface of the M-M recurrence and live here, in the gate directory of
the package they left. Each was a method or a property of the grid and is now a free
function of it: ``half_grid.upstream(m)`` is ``upstream(half_grid, m)`` and
``half_grid.downstream_per_ordinate`` is ``downstream_per_ordinate(half_grid)``.
"""
from __future__ import annotations

import numpy as np

from orpheus.sn.angular.closure import _MMHalfGrid, _psi_half_grid_single_level


def downstream_per_ordinate(half_grid: _MMHalfGrid) -> np.ndarray:
    r"""Shape ``(ng, M, nx)`` — downstream half-face flux for each ordinate.

    ``[g, m, i]`` is :math:`\phi_{m+1/2, i, g}` — the downstream face of ordinate
    ``m`` (equivalently, the upstream face of ordinate ``m+1``). Excludes the
    leading face ``faces[:, 0, :]``, the Carlson seed.
    """
    return half_grid.faces[:, 1:, :]


def upstream(half_grid: _MMHalfGrid, m: int) -> np.ndarray:
    r"""Shape ``(ng, nx)`` — upstream half-face of ordinate ``m``.

    Returns :math:`\phi_{m-1/2, i, g}` per ``(g, i)``.
    """
    return half_grid.faces[:, m, :]


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
