r"""The :math:`G`-orthonormal basis of the span of a set of generators, tabulated over an index set.

A basis computed from data, not from a formula: the caller hands a SPANNING set
of columns over a finite support (possibly over-complete, possibly not
:math:`G`-orthogonal), and :meth:`OrthonormalSpanBasis.of_generators` returns
the table :math:`\Phi` with

.. math::

    \operatorname{span}\Phi = \operatorname{span}(\text{generators}),
    \qquad
    \Phi^{\mathsf T} G \Phi = I ,

by ONE :math:`\sqrt G`-weighted SVD: :math:`\sqrt G\,B = U S V^{\mathsf T}`,
the rank cut on :math:`S`, then :math:`\Phi = U_{:,:r} / \sqrt G`, so
:math:`\Phi^{\mathsf T} G \Phi = U_{:,:r}^{\mathsf T} U_{:,:r} = I` exactly. The
one factorisation does both jobs a derivation would otherwise spell twice (a
rank reduction and an orthonormalisation).

Why orthonormal, and not the generators themselves: a frame's
:class:`~orpheus.numerics.frame_faces.CrossGramInverse` reads a ``DIAGONAL``
Gram by the row-sum probe ``analysis(reconstruction(ones))``, which equals the
true diagonal only if the Gram IS diagonal. Generators that share support are
not :math:`G`-orthogonal (``[M]`` the SN loss-kernel pair generators: 43 %
off-diagonal at ``(3,4,5)``), so declaring ``DIAGONAL`` on them would normalise
every coefficient by the wrong number, silently. The orthonormalisation is what
earns the declaration.

Its first consumer is the kernel gauge
(:class:`~orpheus.numerics.kernel_projection_gauge.KernelProjectionGauge`), one
summand per basis; the SN loss kernel's blocks were its origin (the space-
membership campaign's U3, 2026-10-31, moved the table here from
``orpheus.sn.operators.loss_kernel_gauge.LossKernelBasis``).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from orpheus.geometry.manifold import Manifold
from orpheus.numerics.basis.base import Basis, GramStructure
from orpheus.numerics.coefficient_layout import ModeLayout
from orpheus.numerics.index_set import IndexSet

__all__ = ["OrthonormalSpanBasis", "SPAN_RANK_RTOL"]

#: Relative singular-value cut for the rank of a span. The populations it
#: separates are far apart on every measured input, so the cut is not a tuning
#: knob: ``[M]`` (2026-10-30, ``scratch/membership_covector/u3/rank_gap.py``) the
#: over-complete numerics fixtures keep ``s/s0`` of ``0.598`` and ``0.727`` and
#: discard ``2.9e-17`` and ``2.0e-17``; the SN loss kernel's 50 blocks over six
#: reflective boxes keep at least ``3.68e-1`` and discard at most ``5.5e-16``.
#: Any cut in ``[1e-14, 1e-2]`` returns the same rank on all of them; the SN
#: rank is independently checked against its counting law
#: (``predicted_kernel_dimension``).
SPAN_RANK_RTOL = 1e-10


@dataclass(frozen=True, eq=False)
class OrthonormalSpanBasis(Basis):
    r"""A precomputed :math:`G`-orthonormal table over a finite support.

    ``eq=False`` (identity equality) follows
    :class:`~orpheus.numerics.basis.indicator_basis.IndicatorBasis`: the fields
    are NumPy arrays, which have no value equality a dataclass can use.

    Parameters
    ----------
    table
        ``(n_support, n_modes)``, :math:`G`-orthonormal by construction (mint it
        with :meth:`of_generators`).
    label
        The support's name: the basis's :attr:`domain` is the index set of that
        label, and the frame's measure takes its support from the basis, so the
        two halves of one frame name ONE manifold. The caller names it (an SN
        kernel block is ``sn_trace_orbit{orbit}_g{group}``), so two summands'
        coefficient spaces never compare equal.
    unit
        The coordinates of the constant function's :math:`G`-orthogonal
        projection onto the span (ruling 37), one per mode. Computed by
        :meth:`of_generators` as :math:`\Phi^{\mathsf T} G \mathbf 1` unless the
        caller declares the span orthogonal to the constants, in which case it
        is exactly zero (a computed projection would then be round-off, and the
        coefficient space's identity would hang on arithmetic order: the SN
        kernel's declared zero, the user's ruling of 2026-10-21).
    """

    table: NDArray
    label: str
    unit: NDArray

    @classmethod
    def of_generators(
        cls,
        generators: NDArray,
        weights: NDArray,
        *,
        label: str,
        orthogonal_to_constants: bool = False,
    ) -> "OrthonormalSpanBasis":
        r"""The :math:`G`-orthonormal basis of ``span(generators)``, :math:`G = \operatorname{diag}(w)`.

        ``generators`` is ``(n_support, n_generators)``; ``weights`` is the
        metric's diagonal on the support, ``(n_support,)``, every entry
        positive (a zero-metric row has no minimum-norm representative to
        choose, and the :math:`\sqrt G` division would fail). A rank-0 span
        returns a ``(n_support, 0)`` table.
        """
        columns = np.asarray(generators, dtype=float)
        metric = np.asarray(weights, dtype=float)
        if columns.ndim != 2 or metric.shape != (columns.shape[0],):
            raise ValueError(
                f"OrthonormalSpanBasis {label!r}: generators must be (n_support, n_generators) and weights "
                f"(n_support,), got {columns.shape} and {metric.shape}."
            )
        if not np.all(metric > 0.0):
            raise ValueError(
                f"OrthonormalSpanBasis {label!r}: {int(np.sum(metric <= 0.0))} support row(s) carry a vanishing "
                f"metric. Every value on such a row has the same (zero) G-norm, so there is no minimum-norm "
                f"member to choose, and the sqrt(G) orthonormalisation would divide by zero: exclude those rows "
                f"from the support."
            )
        root = np.sqrt(metric)
        left, singular, _ = np.linalg.svd(root[:, None] * columns, full_matrices=False)
        rank = int(np.sum(singular > SPAN_RANK_RTOL * singular[0])) if singular.size else 0
        table = left[:, :rank] / root[:, None]
        unit = np.zeros(rank) if orthogonal_to_constants else table.T @ metric
        return cls(table=table, label=label, unit=unit)

    # ── the projection-validity declaration ───────────────────────────
    @property
    def gram_structure(self) -> GramStructure:
        r"""``DIAGONAL``, a theorem here: :math:`\Phi^{\mathsf T} G \Phi = U^{\mathsf T} U = I`.

        The frame's row-sum probe is therefore not merely valid (all that
        ``DIAGONAL`` claims) but returns all ones, and :math:`G^{-1}` is the
        identity. It would be false for raw generators that share support.
        """
        return GramStructure.DIAGONAL

    # ── tabulation ────────────────────────────────────────────────────
    def evaluate(self, points: NDArray, /) -> NDArray:
        """Return the precomputed table, validating the caller's node count.

        The modes were fixed when the basis was minted; there is nothing to
        evaluate AT a point, so ``points`` is accepted for the :class:`Basis`
        contract and checked against the table's row count (the
        :class:`~orpheus.numerics.basis.overlap_basis.OverlapBasis` idiom). That
        check is what stops a support/basis mismatch from becoming a silently
        mis-shaped contraction.
        """
        rows = int(np.asarray(points).shape[0])
        if rows != self.table.shape[0]:
            raise ValueError(
                f"OrthonormalSpanBasis {self.label!r} was built for {self.table.shape[0]} support rows but was "
                f"evaluated at {rows} nodes: the measure and the basis describe different supports."
            )
        return self.table

    # ── table contractions ────────────────────────────────────────────
    # The canonical dual factor is the IDENTITY: for an orthonormal basis the
    # dual frame IS the frame, so `reconstruct` and `synthesize` coincide.
    def _synthesize(self, coefficients: NDArray, table: NDArray, /) -> NDArray:
        return np.einsum("nk,k...->n...", table, coefficients)

    def _analyze(self, values: NDArray, table: NDArray, weights: NDArray, /) -> NDArray:
        return np.einsum("n,nk,n...->k...", weights, table, values)

    def _analyze_transpose(self, coefficients: NDArray, table: NDArray, weights: NDArray, /) -> NDArray:
        return np.einsum("n,nk,k...->n...", weights, table, coefficients)

    def _reconstruct(self, coefficients: NDArray, table: NDArray, /) -> NDArray:
        return np.einsum("nk,k...->n...", table, coefficients)

    def _reconstruct_transpose(self, values: NDArray, table: NDArray, /) -> NDArray:
        return np.einsum("nk,n...->k...", table, values)

    @property
    def domain(self) -> Manifold:
        r"""The support, as an index set named by :attr:`label`.

        A finite set with no metric structure, which is exactly
        :class:`~orpheus.numerics.index_set.IndexSet`; the frame's measure takes
        its ``support`` from here, so the basis and the measure name one
        manifold.
        """
        return IndexSet(label=self.label)

    @property
    def layout(self) -> ModeLayout:
        """One amplitude per mode, with the unit the basis carries."""
        return ModeLayout(self.unit, "span_modes")

    @property
    def coefficient_weights(self) -> None:
        """``None``: the modes are :math:`G`-orthonormal, so the Gram is :math:`I`."""
        return None
