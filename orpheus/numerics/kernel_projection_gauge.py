r"""The kernel gauge: the :math:`G`-orthogonal projector onto an explicit kernel, and its section.

A source posing :math:`A\psi = q` with a non-trivial :math:`\ker A` and a
consistent :math:`q` has a COSET of solutions :math:`\psi_0 + \ker A`. Returning
one member means choosing a section of the quotient; the canonical one is the
minimum-:math:`G`-norm member

.. math::

    \sigma(\psi) \;=\; \psi - \Pi\psi ,
    \qquad
    \Pi \;=\; \sum_b \gamma_b^{\mathsf T}\, R_b\, G_b^{-1} M_b\, \gamma_b ,

with :math:`\Pi` the :math:`G`-orthogonal projector onto :math:`\ker A`, here a
sum of SUMMANDS. Each summand is the span of a set of generators over a finite
support: :math:`\gamma_b` gathers the support (a
:class:`~orpheus.numerics.elementary_operators.TraceRestrictionOperator`), the
span's :math:`G`-orthonormal basis
(:class:`~orpheus.numerics.basis.orthonormal_span_basis.OrthonormalSpanBasis`)
and the metric's restriction to the support make a
:class:`~orpheus.numerics.frame.GalerkinFrame`, and the summand's projector is
that frame's ``conjugate(frame, CrossGramInverse.of(frame))``
:math:`= R_b G_b^{-1} M_b`, the one spelling of "project onto a span" in the
package.

**When the sum is a projector.** A sum of :math:`G`-orthogonal projectors is a
projector iff their ranges are mutually :math:`G`-orthogonal. Disjoint supports
give that for free; overlapping supports are admitted only when the summands'
bases are measured :math:`G`-orthogonal (the union of a discretization kernel
and a physical one overlaps on the boundary trace, the space-membership
campaign's R3), and refused otherwise.

**What the gauge is, and what it is not.** This module owns the projector and
its section. Whether a posing is GAUGED at all (rather than solved, or refused
because :math:`q` is not orthogonal to :math:`\ker A^{\mathsf H}`) is the
posing's decision (the Fredholm alternative, the campaign's ruling R3), and
which basis spans the kernel is the method's: SN's diamond-difference loss
kernel supplies its pair generators
(:func:`orpheus.sn.operators.loss_kernel_gauge.predicted_discretization_kernel`),
inside the kernel the posing's decision reads
(:class:`orpheus.numerics.fredholm.PredictedKernel`).

Its laws are gated on a numerics fixture with a non-uniform :math:`G`
(``tests/gates/numerics/test_kernel_projection_gauge.py``): idempotence of
:math:`\Pi` and of :math:`\sigma`, :math:`G`-self-adjointness, range = span,
:math:`\sigma \perp_G \ker`, minimum norm, residual neutrality
(:math:`A\sigma(\psi) = A\psi`), and the overlap rule.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from functools import reduce
from operator import add
from typing import Callable, Optional

import numpy as np
from numpy.typing import NDArray

from orpheus.numerics import frame_faces
from orpheus.numerics.basis.orthonormal_span_basis import OrthonormalSpanBasis
from orpheus.numerics.elementary_operators import TraceRestrictionOperator
from orpheus.numerics.frame import GalerkinFrame
from orpheus.numerics.measure import DiscreteMeasure
from orpheus.numerics.operator import LinearOperator, ZeroMorphism
from orpheus.numerics.space import FunctionSpace, act_on_flat_coordinates

__all__ = ["KernelProjectionGauge", "KernelSummand", "SUMMAND_ORTHOGONALITY_TOL", "diagonal_metric"]

#: The largest ``|G-cosine|`` between two overlapping summands' orthonormal modes
#: that still counts as orthogonal. Orthogonal-by-theory summands read round-off
#: and non-orthogonal ones read O(1), so the cut is not a tuning knob: ``[M]``
#: (2026-10-30, ``scratch/membership_covector/u3/union_orthogonality.py``) the
#: R3 union kernel (the SN diamond-difference modes against the physical constant
#: mode, 28 modes over 4 boxes) reads at most ``1e-17``; the gate's verdict is
#: the same for any cut between about ``1e-14`` and ``1e-8``.
SUMMAND_ORTHOGONALITY_TOL = 1e-12


@dataclass(frozen=True, eq=False)
class KernelSummand:
    r"""One summand of a kernel: a spanning set of generators over a support.

    Parameters
    ----------
    support
        Sorted, unique flat indices into the gauge's space.
    generators
        ``(support.size, n_generators)`` columns spanning the summand's part of
        the kernel; over-complete and non-orthogonal sets are legal (the basis
        mint does the rank cut and the :math:`G`-orthonormalisation).
    label
        The support's name, which the summand's basis and measure share
        (:attr:`OrthonormalSpanBasis.label`).
    orthogonal_to_constants
        The caller's declaration that the span is :math:`G`-orthogonal to the
        constants, so the basis's unit is exactly zero rather than round-off
        (:class:`OrthonormalSpanBasis`). A declaration owes a witness in the
        caller's tests.
    """

    support: NDArray
    generators: NDArray
    label: str
    orthogonal_to_constants: bool = False

    def __post_init__(self) -> None:
        support = np.asarray(self.support)
        generators = np.asarray(self.generators, dtype=float)
        if support.ndim != 1 or not np.issubdtype(support.dtype, np.integer):
            raise ValueError(
                f"KernelSummand {self.label!r}: the support is a 1-D array of integer flat indices, got "
                f"{support.dtype} of shape {support.shape}."
            )
        if support.size and (np.any(np.diff(support) <= 0) or support[0] < 0):
            raise ValueError(
                f"KernelSummand {self.label!r}: the support must be sorted, unique and non-negative."
            )
        if generators.ndim != 2 or generators.shape[0] != support.size:
            raise ValueError(
                f"KernelSummand {self.label!r}: generators must be (support.size, n_generators) = "
                f"({support.size}, k), got {generators.shape}."
            )
        object.__setattr__(self, "support", support.astype(np.intp, copy=False))
        object.__setattr__(self, "generators", generators)


@dataclass(frozen=True, slots=True)
class _ProjectedSummand:
    r"""One term :math:`\gamma_b^{\mathsf T} \Pi_b \gamma_b` of the projector, built."""

    gather: TraceRestrictionOperator
    projector: LinearOperator
    basis: OrthonormalSpanBasis


def diagonal_metric(space: FunctionSpace) -> NDArray:
    r"""The diagonal of the space's metric on its flat chart (C order), :meth:`FunctionSpace.metric_diagonal`.

    SCOPE-BOUNDARY[guard]: a space whose metric is not diagonal is refused.
    machinery: a summand frame over a dense-metric measure; ruling: the space-membership spec §U3.1 (2026-10-30);
    revisit: a kernel posed on a space whose metric is a dense form on a modal block.

    A summand frame's measure is a diagonal weight vector, so a dense Gram has
    no node weights to restrict to a support. (A multi-axis space is admitted
    since U0, 2026-10-31: a summand's support is flat indices into the space's
    C-order flat chart, :func:`_flat_chart`.)
    """
    diagonal = space.metric_diagonal()
    if diagonal is None:
        raise ValueError(
            f"KernelProjectionGauge on {space!r}: the metric is not diagonal. A summand's frame measure is a "
            f"diagonal weight vector, so this is a declared scope boundary (SCOPE-BOUNDARY: the projector is "
            f"defined for diagonal metrics)."
        )
    return np.asarray(diagonal, dtype=float).reshape(-1)


def _flat_chart(space: FunctionSpace) -> FunctionSpace:
    r"""The space a summand's flat indices address: ``space`` itself when one-dimensional, else its C-order flat chart.

    The coordinates of a multi-axis leaf space are an array of its shape; their
    C-order ravel is a vector of :math:`\mathbb R^n` with the ravelled
    diagonal metric. The projector is built on that chart and acts on an
    element by ravelling it in and reshaping it out (:meth:`KernelProjectionGauge.apply`).
    """
    if len(space.shape) == 1:
        return space
    from orpheus.numerics.metric import DiagonalMetric

    return FunctionSpace(f"{space.name}#flat", (space.n_coordinates,),
                         metric=DiagonalMetric(diagonal_metric(space)))


class KernelProjectionGauge(LinearOperator):
    r"""The :math:`G`-orthogonal projector :math:`\Pi` onto a kernel given by spanning summands.

    An endomorphism of ``space``. :meth:`gauge` is the section
    :math:`\psi - \Pi\psi`, the minimum-:math:`G`-norm member of the coset;
    :attr:`dimension` is the projector's rank, ``0`` for the empty direct sum
    (whose section is the identity, so no caller needs a ``None`` branch).

    Not invertible (a projector onto a proper subspace never is): spelled by
    :attr:`is_invertible` and the absence of an ``inverse`` method, the house
    convention (``TraceRestrictionOperator``). Mint it with
    :meth:`direct_sum`.
    """

    def __init__(self, summands: tuple[_ProjectedSummand, ...], space: FunctionSpace) -> None:
        self._summands = tuple(summands)
        self._space = space
        self._chart = _flat_chart(space)
        #: :math:`\Pi = \sum_b \gamma_b^{\mathsf H} \Pi_b \gamma_b`, spelled in the operator algebra,
        #: so its transpose and adjoint follow by law; the zero map when there is no summand.
        self._projector: LinearOperator = (
            reduce(add, (s.gather.H @ s.projector @ s.gather for s in self._summands))
            if self._summands else ZeroMorphism(domain=self._chart, codomain=self._chart)
        )

    @classmethod
    def direct_sum(
        cls, space: FunctionSpace, summands: "tuple[KernelSummand, ...]",
    ) -> "KernelProjectionGauge":
        r"""Build :math:`\Pi` from spanning summands on ``space``.

        A summand of rank 0 contributes nothing and is dropped. Overlapping
        supports are admitted iff the summands' modes are :math:`G`-orthogonal
        to :data:`SUMMAND_ORTHOGONALITY_TOL`; otherwise the sum would not be a
        projector, and it is refused.
        """
        weights = diagonal_metric(space)
        n_total = weights.size
        labels = [summand.label for summand in summands]
        if len(set(labels)) != len(labels):
            raise ValueError(
                f"KernelProjectionGauge: summand labels must be distinct (each names its support's manifold, "
                f"so two equal labels would make two coefficient spaces compare equal); got {labels}."
            )
        built: list[_ProjectedSummand] = []
        for summand in summands:
            support = summand.support
            if support.size == 0:
                continue
            if support[-1] >= n_total:
                raise ValueError(
                    f"KernelSummand {summand.label!r}: support index {int(support[-1])} is outside the space's "
                    f"{n_total} coordinates."
                )
            basis = OrthonormalSpanBasis.of_generators(summand.generators, weights[support], label=summand.label)
            if basis.table.shape[1] == 0:
                continue
            if summand.orthogonal_to_constants:
                basis = _declared_orthogonal_to_constants(basis, weights[support])
            frame = GalerkinFrame(basis, DiscreteMeasure(
                nodes=support.astype(float), weights=weights[support], support=basis.domain,
            ))
            built.append(_ProjectedSummand(
                gather=TraceRestrictionOperator(
                    support, n_total=n_total, domain=_flat_chart(space), codomain=frame_faces.measure_space(frame)),
                projector=frame_faces.conjugate(frame, frame_faces.CrossGramInverse.of(frame)),
                basis=basis,
            ))
        _require_mutually_orthogonal(tuple(built), weights)
        return cls(tuple(built), space)

    # ── the arrow ─────────────────────────────────────────────────────
    @property
    def domain(self) -> Optional[FunctionSpace]:
        return self._space

    @property
    def codomain(self) -> Optional[FunctionSpace]:
        return self._space

    @property
    def summands(self) -> tuple[_ProjectedSummand, ...]:
        """The built summands, in the order given (rank-0 ones dropped)."""
        return self._summands

    @property
    def dimension(self) -> int:
        r"""The rank of :math:`\Pi`, :math:`\dim` of the kernel it spans."""
        return sum(int(s.basis.table.shape[1]) for s in self._summands)

    # ── the action ────────────────────────────────────────────────────
    def apply(self, x: NDArray) -> NDArray:
        r"""The component of ``x`` inside the kernel, :math:`\Pi x`."""
        return self._on_chart(x, self._projector.apply)

    def apply_transpose(self, x: NDArray) -> NDArray:
        r""":math:`\Pi^{\mathsf T} = \sum_b \gamma_b^{\mathsf T} \Pi_b^{\mathsf T} \gamma_b`, the Euclidean transpose.

        Not :math:`\Pi` itself: a :math:`G`-orthogonal projector
        :math:`\Phi\Phi^{\mathsf T}G` is :math:`G`-self-adjoint, and its
        Euclidean transpose :math:`G\Phi\Phi^{\mathsf T}` equals it only when
        :math:`G` is uniform (the SN predecessor returned :math:`\Pi`; nothing
        read it). It is the operator algebra's transpose of the sum, by law.
        """
        return self._on_chart(x, self._projector.dual().apply)

    def _on_chart(self, x: NDArray, action: "Callable[[NDArray], NDArray]") -> NDArray:
        r"""Parse ``x`` at the space's action verb, act on its flat chart, and return it in the space's shape."""
        return act_on_flat_coordinates(self._space, self._space, x, action)

    @property
    def is_adjointable(self) -> bool:
        return True

    @property
    def is_invertible(self) -> bool:
        """``False``: a projector onto a proper subspace kills its complement."""
        return False

    def gauge(self, state: NDArray) -> NDArray:
        r"""The section :math:`(I - \Pi)\,\psi`: the canonical, minimum-:math:`G`-norm member.

        Residual-neutral: :math:`A(\psi - \Pi\psi) = A\psi` because
        :math:`\Pi\psi \in \ker A`, so no convergence check can move.
        """
        values = np.asarray(self._space.leading_coordinates_of(state), dtype=float)
        return values - self.apply(values)

    def displacement(self, state: NDArray) -> float:
        r"""How far the section moves ``state``: :math:`\lVert\Pi\psi\rVert_G / \lVert\psi\rVert_G`, in the space's own norm.

        ``0.0`` for the zero state and for the empty direct sum. The exit
        report's gauge reading on every method that gauges. A scalar verb, so
        ``state`` is parsed exactly (:meth:`~orpheus.numerics.space.FunctionSpace.coordinates_of`):
        a batch has no single displacement.
        """
        values = np.asarray(self._space.coordinates_of(state), dtype=float).reshape(-1)
        weights = diagonal_metric(self._space)
        total = float(np.sqrt(np.sum(weights * values**2)))
        if total == 0.0 or not self._summands:
            return 0.0
        component = np.asarray(self.apply(values.reshape(self._space.shape)), dtype=float).reshape(-1)
        return float(np.sqrt(np.sum(weights * component**2))) / total

    def __repr__(self) -> str:
        return (
            f"KernelProjectionGauge(summands={len(self._summands)}, "
            f"dimension={self.dimension}, space={self._space!r})"
        )


def _require_mutually_orthogonal(summands: tuple[_ProjectedSummand, ...], weights: NDArray) -> None:
    r"""Refuse two summands that share support rows unless their modes are :math:`G`-orthogonal."""
    if len(summands) < 2:
        return
    owners: dict[int, list[int]] = {}
    counts = np.zeros(weights.size, dtype=np.intp)
    for s in summands:
        counts[s.gather.indices] += 1
    if counts.max() <= 1:
        return
    shared = set(np.flatnonzero(counts > 1).tolist())
    for k, s in enumerate(summands):
        for row in shared.intersection(s.gather.indices.tolist()):
            owners.setdefault(row, []).append(k)
    pairs = {(a, b) for ks in owners.values() for i, a in enumerate(ks) for b in ks[i + 1:]}
    for a, b in sorted(pairs):
        sa, sb = summands[a], summands[b]
        common, ia, ib = np.intersect1d(sa.gather.indices, sb.gather.indices, return_indices=True)
        cosine = sa.basis.table[ia].T @ (weights[common][:, None] * sb.basis.table[ib])
        worst = float(np.max(np.abs(cosine)))
        if worst > SUMMAND_ORTHOGONALITY_TOL:
            raise ValueError(
                f"KernelProjectionGauge: summands {sa.basis.label!r} and {sb.basis.label!r} overlap on "
                f"{common.size} row(s) and are not mutually G-orthogonal (largest |G-cosine| {worst:.3e} > "
                f"{SUMMAND_ORTHOGONALITY_TOL:.0e}); their projectors do not sum to a projector. Merge them "
                f"into one summand."
            )


def _declared_orthogonal_to_constants(basis: OrthonormalSpanBasis, weights: NDArray) -> OrthonormalSpanBasis:
    r"""Check the caller's declaration that the span is :math:`G`-orthogonal to the constants, then store the exact zero.

    The computed unit :math:`\Phi^{\mathsf T} G \mathbf 1` holds the
    modes' :math:`G`-inner products with the constant; divided by
    :math:`\lVert \mathbf 1\rVert_G` it is a :math:`G`-cosine, held to
    :data:`SUMMAND_ORTHOGONALITY_TOL` like an overlap. A declaration that
    the measurement contradicts is refused; one it confirms stores the
    exact zero, so the coefficient space's identity does not hang on
    round-off (the SN loss kernel's declared zero, the user's ruling of
    2026-10-21).
    """
    cosine = float(np.max(np.abs(basis.unit))) / float(np.sqrt(np.sum(weights)))
    if cosine > SUMMAND_ORTHOGONALITY_TOL:
        raise ValueError(
            f"KernelSummand {basis.label!r} is declared orthogonal to the constants, but its largest "
            f"|G-cosine| with the constant is {cosine:.3e} > {SUMMAND_ORTHOGONALITY_TOL:.0e}."
        )
    return replace(basis, unit=np.zeros(basis.table.shape[1]))
