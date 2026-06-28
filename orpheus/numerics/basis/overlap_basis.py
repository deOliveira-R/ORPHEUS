r"""Fractional-overlap trial basis — the non-nested generalization of IndicatorBasis.

:class:`OverlapBasis` carries a precomputed **fractional** membership table
:math:`T[g,G]\in[0,1]` (a partition of unity, rows summing to 1) instead of the
point-sampled one-hot table of
:class:`~orpheus.numerics.basis.indicator_basis.IndicatorBasis`. A fine cell that
straddles a coarse boundary belongs *fractionally* to each coarse cell it overlaps —
the conservative re-binning that a NON-nested condensation needs (see
:class:`orpheus.data.energy_grid.GroupCondensation`). The one-hot
:class:`IndicatorBasis` is the **nested degenerate** (no straddles → every
:math:`T[g,G]\in\{0,1\}`).

Why a subclass (and why this is a no-op extension, not a new arm)
================================================================

Every table contraction on :class:`IndicatorBasis` — :meth:`analyze`,
:meth:`reconstruct`, :meth:`synthesize`, the transposes, and the diagonal Gram a
:class:`~orpheus.numerics.frame.FrameBase` builds via
``analysis(reconstruction(ones))`` — is a **pure function of the membership table**
and never assumes it was one-hot. So the fractional generalization reuses ALL of
them unchanged; only the table *production* differs. :class:`OverlapBasis` overrides
exactly one method, :meth:`evaluate`, to return the precomputed overlap table rather
than bucketing points with ``searchsorted``. The rate-preserving projection then
falls out of the existing ``frame.project = G⁻¹ M`` because a partition-of-unity
table gives ``reconstruction(ones) = 1`` → the diagonal Gram :math:`\Phi_G =
\sum_g T[g,G]\,\varphi_g`.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from orpheus.numerics.basis.indicator_basis import IndicatorBasis

__all__ = ["OverlapBasis"]


@dataclass(frozen=True, eq=False)
class OverlapBasis(IndicatorBasis):
    r"""Cell-indicator basis carrying a precomputed FRACTIONAL membership table.

    Parameters
    ----------
    edges_per_axis : tuple[NDArray, ...]
        The coarse-index partition (inherited from :class:`IndicatorBasis`), used
        only for the coarse cell count :attr:`n_cells` and the coefficient
        :attr:`space` — the membership table, not these edges, is the source of
        truth (:meth:`evaluate` is overridden). For the energy axis this is a
        length-1 tuple of the coarse-group index edges.
    overlap_table : NDArray, shape ``(n_fine, n_coarse)``
        The precomputed partition-of-unity table :math:`T[g,G]\in[0,1]`,
        ``rows.sum(axis=1) == 1``. One-hot recovers :class:`IndicatorBasis`.

    Notes
    -----
    Frozen, ``eq=False`` (identity equality; the fields are NumPy arrays).

    The inherited :meth:`~orpheus.numerics.basis.indicator_basis.IndicatorBasis.mass_matrix`
    docstring claims a **diagonal** Gram ("the indicators have disjoint support") —
    that claim does **not** hold for a fractional table (a straddling row makes two
    columns share support, so the cross Gram is non-diagonal). It is latent: no
    consumer calls ``mass_matrix`` (the frame's :meth:`~orpheus.numerics.frame.FrameBase.project`
    uses the partition-of-unity row-sum probe, not the full Gram — see
    :attr:`~orpheus.numerics.frame.FrameBase.gram`). A future least-squares consumer
    that needs the dense Gram must compute it for the fractional case, not trust the
    inherited diagonal claim.
    """

    overlap_table: NDArray

    def evaluate(self, points: NDArray, /) -> NDArray:
        r"""Return the precomputed ``(n_fine, n_coarse)`` partition-of-unity table.

        Interval overlap is computed once when the table is built (it depends on the
        fine AND coarse partitions plus the within-group weight, not on sample
        points), so this returns it directly. ``points`` (the fine-group nodes) are
        accepted for the :class:`~orpheus.numerics.basis.base.Basis` contract and
        validated against the table's fine-row count.
        """
        pts = np.asarray(points)
        n_rows = int(self.overlap_table.shape[0])
        if pts.shape[0] != n_rows:
            raise ValueError(
                f"OverlapBasis table has {n_rows} fine rows but got "
                f"{pts.shape[0]} sample points."
            )
        return self.overlap_table
