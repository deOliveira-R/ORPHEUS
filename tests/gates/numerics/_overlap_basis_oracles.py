"""Test-side oracle of the re-binning provenance of an :class:`~orpheus.numerics.basis.OverlapBasis` (moved from production, dead-B row 43).

The table ``overlap_basis.overlap_table`` carries the provenance; production reads only the table, so the diagnostic lives here.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray


def fractional_columns(overlap_basis) -> NDArray:
    r"""Coarse columns that received a STRICTLY-FRACTIONAL contribution.

    ``(k,)`` indices of coarse cells/groups where some fine row straddled the boundary and contributed a weight strictly
    between 0 and 1, i.e. where the re-binning leaned on the within-group model rather than a clean (one-hot) nesting. Empty
    for a nested table; non-empty where the coarse partition is locally finer than the fine one. For the energy axis these
    are the coarse groups whose data leaned on the within-group flux model :math:`w(E)`.
    """
    table = overlap_basis.overlap_table
    frac = (table > 1e-12) & (table < 1.0 - 1e-12)
    return np.nonzero(frac.any(axis=0))[0]
