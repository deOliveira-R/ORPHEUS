r"""The Fredholm alternative on a posed linear question :math:`A\psi = q`, decided from a PREDICTED kernel.

A finite-dimensional linear operator :math:`A : V \to V` answers a source
:math:`q` in exactly one of three ways (the Fredholm alternative):

* :math:`\ker A = \{0\}`: one solution. The verdict is :class:`Solve`.
* :math:`\ker A \neq \{0\}` and :math:`q` is annihilated by every co-vector of
  :math:`\ker A^{\mathsf T}` (equivalently :math:`q \in \operatorname{range} A`):
  a COSET of solutions :math:`\psi_0 + \ker A`. The verdict is :class:`Gauge`,
  carrying the kernel projector whose section
  (:meth:`~orpheus.numerics.kernel_projection_gauge.KernelProjectionGauge.gauge`)
  picks the minimum-:math:`G`-norm member.
* otherwise: no solution. The decision raises :class:`NoSolution`, so no
  outcome ever holds a refused posing.

The kernel is never measured in production. It is PREDICTED from structure by
the method (:class:`PredictedKernel`: the right generators spanning
:math:`\ker A` and the left generators spanning :math:`\ker A^{\mathsf T}`, the
annihilators of the range), from declared data alone: the absorbing-chain
predicate of the transport layer, SN's diamond-difference face transmission.
A predictor's claim is checked against a dense SVD in the gates, never at run
time. (The user's ruling R3, 2026-10-30: "Fredholm as the decision rule,
kernel measurement as witness".)

**The annihilators are co-vectors.** :math:`\ker A^{\mathsf T}` lives in the
dual :math:`V^*`, and the test :math:`\langle w, q\rangle = 0` is the canonical
pairing, coordinates contracted with no metric. Its scale-free reading is the
coordinate cosine ON THE ANNIHILATOR'S SUPPORT,

.. math::

    c(w, q) \;=\; \frac{|\langle w, q\rangle|}{\lVert w\rVert_2\,\lVert q|_{\operatorname{supp} w}\rVert_2},

the largest over the left generators' orthonormalised span: round-off inside
the support enters squared, and a source living off the support (a bulk
source against SN's trace-only DD annihilators) cannot dilute it.

**The adjoint question daggers the kernel.** :math:`\ker A^{\mathsf H} =
G^{-1}\ker A^{\mathsf T}` and the annihilators of :math:`\operatorname{range}
A^{\mathsf H}` are :math:`G\ker A`, so :attr:`PredictedKernel.H` swaps the two
generator sets through the Riesz legs; an adjoint posing gauges with the LEFT
kernel, never the forward projector.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import cached_property

import numpy as np
from numpy.typing import NDArray

from orpheus.numerics.evidence import Evidence, Measured
from orpheus.numerics.kernel_projection_gauge import KernelProjectionGauge, KernelSummand, _metric_diagonal
from orpheus.numerics.operator import LinearOperator
from orpheus.numerics.space import FunctionSpace

__all__ = [
    "CONSISTENCY_RTOL",
    "Gauge",
    "NoSolution",
    "PredictedKernel",
    "SEMISIMPLE_RTOL",
    "Solvable",
    "Solve",
    "fredholm_alternative",
]

#: The largest coordinate cosine :math:`c(w, q)` between a source and an
#: annihilator that still counts as orthogonal. Consistent sources read
#: round-off and inconsistent ones read O(1), so the cut is not a tuning knob:
#: ``[M]`` (2026-10-31, the test-architect's ``scratch/membership_covector/u0/probes/``,
#: ``p_cosine.py``, the cosine on the declared support) the 10 gauged cases
#: (diffusion, SN 1-D and 2-D, the SN adjoint) read at most ``2.4e-15`` and the
#: 6 refused ones at least ``0.092`` (a source on a closed white wall); the cut
#: sits at least 4.6 decades inside each side.
CONSISTENCY_RTOL = 1e-10

#: The smallest principal cosine between :math:`\ker A^{\mathsf T}` and
#: :math:`\ker A` for which the zero eigenvalue counts as semisimple, so that
#: :math:`A + \Pi` is invertible. ``[M]`` (2026-10-31,
#: ``scratch/membership_covector/u0/deflated_probe.py``) the absorbing-chain
#: kernels read 0.68 to 0.96 and SN's diamond-difference kernel reads exactly 0
#: (its left and right kernels are orthogonal: a Jordan chain).
SEMISIMPLE_RTOL = 1e-10


class NoSolution(ValueError):
    r"""The source is not annihilated by :math:`\ker A^{\mathsf T}`: :math:`A\psi = q` has no steady solution."""


@dataclass(frozen=True, eq=False)
class PredictedKernel:
    r""":math:`\ker A` and :math:`\ker A^{\mathsf T}` on ``space``, predicted from structure.

    Parameters
    ----------
    space
        The operator's domain (and codomain: the question is an endomorphism's).
    right
        Spanning summands of :math:`\ker A`, in the coordinates of ``space``.
    left
        Spanning summands of :math:`\ker A^{\mathsf T}`, the annihilators of
        :math:`\operatorname{range} A`, in the coordinates of :math:`V^*`.
    because
        The predictor's witness, in words a refusal or a record can cite
        (the absorption-free closed classes, the reflective axes).
    """

    space: FunctionSpace
    right: tuple[KernelSummand, ...]
    left: tuple[KernelSummand, ...]
    because: str

    def __post_init__(self) -> None:
        right, left = self._right_span.shape[1], self._annihilators.shape[1]
        if right != left:
            raise ValueError(
                f"PredictedKernel ({self.because}): {right} right kernel vectors against {left} left co-vectors. "
                f"An endomorphism of a finite-dimensional space has Fredholm index 0 (dim ker A = dim ker A^T), so "
                f"the predictor has lost or invented a mode."
            )

    @classmethod
    def trivial(cls, space: FunctionSpace, because: str) -> "PredictedKernel":
        r"""The prediction :math:`\ker A = \{0\}`, with its reason."""
        return cls(space=space, right=(), left=(), because=because)

    def __add__(self, other: "PredictedKernel") -> "PredictedKernel":
        r"""The sum of two predicted kernels on one space (SN: the discretization kernel and a physical one)."""
        if other.space != self.space:
            raise ValueError(f"PredictedKernel: cannot add kernels on {self.space!r} and {other.space!r}.")
        reasons = "; ".join(r for r in (self.because, other.because) if r)
        return PredictedKernel(self.space, self.right + other.right, self.left + other.left, reasons)

    # ── the right kernel: its projector ───────────────────────────────
    @cached_property
    def projector(self) -> KernelProjectionGauge:
        r"""The :math:`G`-orthogonal projector :math:`\Pi` onto :math:`\ker A`; its section is the gauge."""
        return KernelProjectionGauge.direct_sum(self.space, self.right)

    @property
    def dimension(self) -> int:
        r""":math:`\dim\ker A`, the projector's rank."""
        return self.projector.dimension

    # ── the left kernel: the annihilators ─────────────────────────────
    @cached_property
    def _annihilators(self) -> NDArray:
        r"""An orthonormal (Euclidean) basis of :math:`\ker A^{\mathsf T}`, embedded, ``(n, k)``."""
        return _embedded_span(self.left, self._n)

    @cached_property
    def _right_span(self) -> NDArray:
        return _embedded_span(self.right, self._n)

    @property
    def _n(self) -> int:
        return int(np.prod(self.space.shape))

    def consistency(self, source: NDArray) -> "tuple[float, str]":
        r"""The largest coordinate cosine :math:`c(w, q)` of ``source`` against the annihilators, and the summand that reads it.

        Per left summand, the cosine between the source restricted to the
        summand's declared support and the summand's span; ``(0.0, "")`` when
        there is no annihilator or the source vanishes on every support.
        """
        q = np.ravel(np.asarray(self.space.coordinates_of(source), dtype=float))
        worst, label = 0.0, ""
        for summand in self.left:
            basis = _orthonormal_columns(summand.generators)
            local = q[summand.support]
            scale = float(np.linalg.norm(local))
            if basis.shape[1] == 0 or scale == 0.0:
                continue
            cosine = float(np.linalg.norm(basis.T @ local)) / scale
            if cosine > worst:
                worst, label = cosine, summand.label
        return worst, label

    @property
    def is_semisimple(self) -> bool:
        r"""Whether :math:`\ker A^{\mathsf T}` pairs nondegenerately with :math:`\ker A`, so :math:`A + \Pi` is invertible.

        The smallest principal cosine between the two spans, against
        :data:`SEMISIMPLE_RTOL`. A trivial kernel is semisimple.
        """
        if self.dimension == 0:
            return True
        left, right = self._annihilators, self._right_span
        if left.shape[1] != right.shape[1]:
            return False
        return float(np.linalg.svd(left.T @ right, compute_uv=False).min()) > SEMISIMPLE_RTOL

    def deflate(self, loss: LinearOperator) -> LinearOperator:
        r""":math:`A + \Pi`, invertible when the kernel is semisimple: its inverse maps a consistent source to the minimum-:math:`G`-norm member.

        For :math:`q \in \operatorname{range} A`, :math:`x = (A + \Pi)^{-1} q`
        satisfies :math:`\Pi x = 0` and :math:`A x = q`: pairing the equation
        with a left co-vector :math:`y` gives :math:`y^{\mathsf T}\Pi x = 0`,
        and :math:`y^{\mathsf T} R` nonsingular forces :math:`\Pi x = 0`. So a
        factorising method factorises a NONSINGULAR operator and never meets a
        zero pivot. ``[M]`` (``scratch/membership_covector/u0/deflated_probe.py``)
        on three singular diffusion losses, cond :math:`(A + \Pi)` 35 to 42
        against :math:`10^{17}` for :math:`A`. A trivial kernel returns ``loss``
        unchanged; a non-semisimple one (a Jordan chain, SN's
        diamond-difference kernel) is refused, since :math:`A + \Pi` stays
        singular there.
        """
        if self.dimension == 0:
            return loss
        if not self.is_semisimple:
            raise ValueError(
                f"PredictedKernel ({self.because}): the kernel is not semisimple (its left co-vectors annihilate its "
                f"right vectors: a Jordan chain at 0), so A + Pi is singular and cannot be factorised; solve "
                f"iteratively and gauge the result instead."
            )
        return loss + self.projector

    # ── the adjoint question ──────────────────────────────────────────
    @property
    def H(self) -> "PredictedKernel":
        r"""The predicted kernel of :math:`A^{\mathsf H} = G^{-1}A^{\mathsf T}G`.

        :math:`\ker A^{\mathsf H} = \sharp\ker A^{\mathsf T}` and its
        annihilators are :math:`\flat\ker A`: the two generator sets swap
        through the Riesz legs, read from the metric's diagonal on each support.
        """
        weights = _metric_diagonal(self.space)
        sharp = tuple(_reweighted(s, 1.0 / weights[s.support], "sharp") for s in self.left)
        flat = tuple(_reweighted(s, weights[s.support], "flat") for s in self.right)
        return PredictedKernel(self.space, sharp, flat, f"the adjoint of: {self.because}")


def _reweighted(summand: KernelSummand, factor: NDArray, leg: str) -> KernelSummand:
    return KernelSummand(summand.support, factor[:, None] * summand.generators, f"{summand.label}_{leg}")


def _orthonormal_columns(generators: NDArray) -> NDArray:
    r"""An orthonormal basis of the generators' column span (Euclidean; the rank cut of the U3 basis mint)."""
    from orpheus.numerics.basis.orthonormal_span_basis import SPAN_RANK_RTOL

    if generators.size == 0:
        return np.zeros((generators.shape[0], 0))
    left, singular, _ = np.linalg.svd(generators, full_matrices=False)
    return left[:, : int(np.sum(singular > SPAN_RANK_RTOL * singular[0]))]


def _embedded_span(summands: tuple[KernelSummand, ...], n: int) -> NDArray:
    columns = []
    for summand in summands:
        basis = _orthonormal_columns(summand.generators)
        block = np.zeros((n, basis.shape[1]))
        block[summand.support] = basis
        columns.append(block)
    if not columns:
        return np.zeros((n, 0))
    span = np.hstack(columns)
    return _orthonormal_columns(span)


# ── the verdict ───────────────────────────────────────────────────────


@dataclass(frozen=True)
class Solve:
    r""":math:`\ker A = \{0\}`: the posing has one solution.

    ``because`` is the predictor's reason; ``conditioning`` what is known of
    :math:`A`'s conditioning (the factorisation's reciprocal-condition
    estimate on a factorising method, ``NotApplicable`` where the inverse is
    iterative).
    """

    because: str
    conditioning: Evidence


@dataclass(frozen=True)
class Gauge:
    r"""A kernel and a consistent source: a coset of solutions, and the section that picks one.

    ``projector`` is the kernel's :math:`G`-orthogonal projector (its
    :meth:`~orpheus.numerics.kernel_projection_gauge.KernelProjectionGauge.gauge`
    is the minimum-norm section); ``pairing`` the source's measured
    :func:`consistency <PredictedKernel.consistency>`, at or below
    :data:`CONSISTENCY_RTOL`.
    """

    projector: KernelProjectionGauge
    pairing: Measured
    because: str
    conditioning: Evidence


Solvable = Solve | Gauge


def fredholm_alternative(kernel: PredictedKernel, source: NDArray, *, conditioning: Evidence) -> Solvable:
    r"""Decide the posing :math:`A\psi = q` from its predicted kernel: :class:`Solve`, :class:`Gauge`, or raise :class:`NoSolution`."""
    if kernel.dimension == 0:
        return Solve(because=kernel.because, conditioning=conditioning)
    pairing, label = kernel.consistency(source)
    if pairing > CONSISTENCY_RTOL:
        raise NoSolution(
            f"the posing has no steady solution: the loss is singular ({kernel.because}), and the source is not "
            f"annihilated by its left kernel ({label!r} reads a coordinate cosine {pairing:.3g} > "
            f"{CONSISTENCY_RTOL:g}). On a closed "
            f"domain with an absorption-free class, every neutron a source puts into the class stays there; only a "
            f"balanced source (zero net emission into the class) has a solution."
        )
    return Gauge(projector=kernel.projector, pairing=Measured(pairing), because=kernel.because,
                 conditioning=conditioning)
