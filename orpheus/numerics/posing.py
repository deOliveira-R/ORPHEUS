r"""The POSINGS — layer 2 of the terminal object: the QUESTION asked of a pencil or an operator.

Shape (B) (``consumers_step2_design.md`` §3.5, RULED 2026-09-13): two kind-typed
types and ZERO Optionals —

* :class:`EigenPosing` ``(pencil, spectral_map)`` — the homogeneous question
  :math:`A\psi = \mu M\psi` with the physical eigenvalue read through a
  spectral map (k: :math:`\mu \mapsto 1/\mu`; α: :math:`\mu \mapsto -1/\mu`).
  The unknown is a RAY in the cone plus a scalar; the adjoint is NULLARY
  (:math:`k^\dagger = k`, the pencil's ``.H`` shares the spectrum).
* :class:`SourcePosing` ``(operator, source)`` — the affine question
  :math:`A\psi = q`.  The unknown is a COSET :math:`\psi_0 + \ker A`; the adjoint
  is UNARY — a source problem's dual needs a DETECTOR, extra data the generating
  data does not carry.

The ``(M and q)`` cell is a COMPOSITION, not a third type:
``SourcePosing(pencil.at(σ), q(σ))`` — the subcritical multiplying source at the
physical σ = 1, the Laplace/noise family as a path over one pencil.  Which cells
the generating data occupies is the PROBLEM's; which functional of the resolvent
is computed within a cell is the STRATEGY's.

ONE balance functional :math:`\beta(\psi, \lambda) = \langle 1, \mathcal{A}(\lambda)\psi - q\rangle`
underlies both kinds — the eigen kind SOLVES :math:`\beta = 0` for λ (a Rayleigh
quotient with the constant weight), the source kind EVALUATES it.  The pairing is
the explicit weight's: these numerics primitives take ``w``; a method's
volume-weighted member (SN's ``compute_keff``) is that quotient with the space's
measure as the weight — a REFERENCE-class agreement (``[M]`` the gap is the
convergence residual, 7e-10 → 7e-14 with the outer count), not one body.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Callable, Generic, TypeVar

import numpy as np

from orpheus.numerics.evidence import Evidence, Measured
from orpheus.numerics.fredholm import (
    CONSISTENCY_RTOL, Gauge, NoSolution, PredictedKernel, Solvable, Solve, fredholm_alternative,
)
from orpheus.numerics.operator import LinearOperator, Vector
from orpheus.numerics.pencil import OperatorPencil

V = TypeVar("V", bound=Vector)  # the operator algebra's own carrier bound


def _pair(w: Any, y: Any) -> float:
    """The explicit-weight pairing ⟨w, y⟩ over flattened arrays (``w = 1`` → the plain sum;
    a scalar or lower-rank ``w`` broadcasts to ``y``'s shape)."""
    y_arr = _flat(y)
    w_arr = np.broadcast_to(np.asarray(w, dtype=float).ravel() if np.ndim(w) else np.asarray(w, dtype=float), y_arr.shape)
    # ``np.sum`` (pairwise) rather than a BLAS dot: with ``w = 1`` this IS the
    # coordinate sum the method-tier estimators take (bit-identical to them).
    return float(np.sum(w_arr * y_arr))


def _flat(y: Any) -> np.ndarray:
    """A carrier-honest flatten: typed composites (``FullField``/``CoupledField``)
    expose ``to_flat``; bare arrays ravel."""
    flat = getattr(y, "to_flat", None)
    return np.asarray(flat() if callable(flat) else y, dtype=float).ravel()


def _carrier(space: Any, flat: np.ndarray) -> Any:
    r"""The element of ``space`` whose flat coordinates are ``flat``: its typed carrier when the space mints one, else the coordinate array."""
    zeros = getattr(space, "zeros", None)
    if callable(zeros):
        template: Any = zeros()
        return type(template).from_flat(flat, template)
    return flat.reshape(space.shape)


def _predicted(pencil: OperatorPencil[Any], question: str) -> PredictedKernel:
    if pencil.kernel is None:
        raise ValueError(
            f"{question}: the pencil carries no predicted kernel, so the question has no Fredholm verdict; pose it on a "
            f"hub's pencil (whose kernel the method predicts) to ask for one."
        )
    return pencil.kernel


def _pencil_annihilation(pencil: OperatorPencil[Any], kernel: PredictedKernel) -> tuple[float, str]:
    r"""How far the rhs :math:`M` fails to respect the kernel of :math:`A`: the largest of :math:`\lVert M^{\mathsf T}w\rVert`
    over each annihilator and :math:`\lVert M r\rVert` over each kernel mode, each relative to the same action on the
    entrywise absolute value, and the summand that reads it.

    Zero means :math:`\ker A \subseteq \ker M` and :math:`\operatorname{range} M \perp \ker A^{\mathsf T}`: then
    :math:`\ker A` is the kernel of every member :math:`A - \mu M` with the same annihilators, the eigen question is a
    solvable family on it, and gauging is residual-neutral at every :math:`\mu`. The relative scale assumes :math:`M`
    is entrywise non-negative, so the absolute pass bounds the signed one with no cancellation. Every production
    pencil's rhs is a fission emission :math:`\chi \otimes \nu\Sigma_f`, non-negative because the mixture parses its
    spectrum as a probability simplex at construction (:func:`~orpheus.data.emission_spectrum.enforce_emission_spectrum`
    refuses a negative entry); ``[M]`` qa, 2026-10-31: SN's assembled rhs, forward and daggered, has smallest entry
    0.0 on 4 of 4 fixtures with (n,2n) and fission (the (n,2n) term is the loss's, not :math:`M`'s).
    :math:`M^{\mathsf T}` is the coordinate transpose (the canonical pairing, no metric).
    """
    from orpheus.numerics.fredholm import _embedded_columns

    space = pencil.lhs.domain
    worst, label = 0.0, ""
    for summands, act in ((kernel.left, pencil.rhs.dual().apply), (kernel.right, pencil.rhs.apply)):
        for summand, column in _embedded_columns(summands, kernel.n_coordinates):
            signed = _flat(act(_carrier(space, column)))
            scale = float(np.linalg.norm(_flat(act(_carrier(space, np.abs(column))))))
            if scale == 0.0:
                continue
            ratio = float(np.linalg.norm(signed)) / scale
            if ratio > worst:
                worst, label = ratio, summand.label
    return worst, label


@dataclass(frozen=True)
class SpectralMap:
    r"""The bijection between the pencil's eigenvalue μ and the physical eigenvalue λ.

    ``of_quotient(a, m)`` states λ as ONE division of the two pairings
    :math:`a = \langle w, A\psi\rangle`, :math:`m = \langle w, M\psi\rangle`
    — ``forward(a / m)`` is the same number mathematically and NOT numerically
    (``[M]`` ``1/(41/6)`` vs ``6/41`` differ in the last digit): the Rayleigh
    estimator must be the method-tier ratio bit-for-bit.
    """

    forward: Callable[[float], float]  # μ ↦ λ
    inverse: Callable[[float], float]  # λ ↦ μ
    of_quotient: Callable[[float, float], float]  # (⟨w,Aψ⟩, ⟨w,Mψ⟩) ↦ λ, one division
    name: str

    def __call__(self, mu: float) -> float:
        return float(self.forward(mu))


K_MAP = SpectralMap(forward=lambda mu: 1.0 / mu, inverse=lambda k: 1.0 / k, of_quotient=lambda a, m: m / a, name="k")
r"""k-eigenvalue: :math:`A\psi = F\psi/k` ⟺ :math:`\mu = 1/k`; :math:`k = \langle w,F\psi\rangle/\langle w,A\psi\rangle`."""

ALPHA_MAP = SpectralMap(forward=lambda mu: -mu, inverse=lambda alpha: -alpha, of_quotient=lambda a, m: -a / m, name="alpha")
r"""α-eigenvalue: :math:`(L+C-S-F)\psi = -\alpha T\psi` ⟺ :math:`\mu = -\alpha` (a later posing; the map is stated now).

μ is the PENCIL's eigenvalue, :math:`A\psi = \mu M\psi` with :math:`M = T = 1/v`,
the convention ``K_MAP`` also reads (:math:`k = 1/\mu`). Where the corpus writes
μ for the eigenvalue of :math:`A^{-1}M` instead (the power iteration's), the
same α is :math:`-1/\mu`; this map once carried that formula (ERR-089)."""


@dataclass(frozen=True)
class EigenPosing(Generic[V]):
    """The eigenvalue question over a pencil, with its spectral map."""

    pencil: OperatorPencil[V]
    spectral_map: SpectralMap

    def residual(self, psi: Any, lam: float) -> Any:
        r"""The eigen-residual :math:`\mathcal{A}(\mu)\psi = A\psi - \mu M\psi` at the PHYSICAL eigenvalue λ."""
        return self.pencil.at(self.spectral_map.inverse(float(lam))).apply(psi)

    def rayleigh(self, psi: Any, w: Any = 1.0) -> float:
        r"""The PHYSICAL eigenvalue estimate :math:`\lambda(\mu)`, :math:`\mu = \langle w, A\psi\rangle / \langle w, M\psi\rangle`.

        ``w = 1`` is the balance-functional member (the eigen kind SOLVING
        β = 0); ``w = ψ†`` the stationary (adjoint-weighted) one.
        """
        a = _pair(w, self.pencil.lhs.apply(psi))
        m = _pair(w, self.pencil.rhs.apply(psi))
        return float(self.spectral_map.of_quotient(a, m))

    def balance(self, psi: Any, lam: float, w: Any = 1.0) -> float:
        r"""β(ψ, λ) = ⟨w, 𝒜(μ(λ))ψ⟩ — zero at the solution."""
        return _pair(w, self.residual(psi, lam))

    def H(self) -> "EigenPosing[V]":
        r"""NULLARY — the adjoint eigen-question needs no datum (:math:`k^\dagger = k`)."""
        return EigenPosing(self.pencil.H, self.spectral_map)

    def solvability(self, *, conditioning: Evidence) -> Solvable:
        r"""The Fredholm verdict on :math:`A\psi = \mu M\psi`, derived from the pencil's predicted kernel.

        A trivial kernel is ``Solve``. Otherwise the question is a solvable
        family iff the kernel is the WHOLE pencil's (:func:`_pencil_annihilation`
        reads zero: :math:`M` annihilates the kernel modes and its range is
        annihilated by the left kernel), and then it is ``Gauge``; else it
        raises :class:`~orpheus.numerics.fredholm.NoSolution` (fission emits
        into, or a kernel mode carries fission from, a class the loss never
        removes from). The condition is symmetric under :attr:`H`, so the
        forward question and its adjoint refuse together.
        """
        kernel = _predicted(self.pencil, "EigenPosing.solvability")
        if kernel.dimension == 0:
            return Solve(because=kernel.because, conditioning=conditioning)
        reading, label = _pencil_annihilation(self.pencil, kernel)
        if reading > CONSISTENCY_RTOL:
            raise NoSolution(
                f"the eigen posing has no steady solution: the loss is singular ({kernel.because}), and the pencil's "
                f"rhs does not respect its kernel ({label!r} reads {reading:.3g} > {CONSISTENCY_RTOL:g}): fission "
                f"emits into a class the loss never removes neutrons from."
            )
        return Gauge(projector=kernel.projector, pairing=Measured(reading), because=kernel.because,
                     conditioning=conditioning)


@dataclass(frozen=True)
class SourcePosing(Generic[V]):
    r"""The affine question :math:`A\psi = q` over one operator and one source."""

    operator: LinearOperator[V, V]
    source: Any
    #: The operator's predicted kernel (:class:`~orpheus.numerics.fredholm.PredictedKernel`), from which
    #: :meth:`solvability` derives the verdict; ``None`` when the operator's kernel was not predicted.
    #: Mint it from a pencil with :meth:`of_pencil`.
    kernel: "PredictedKernel | None" = field(default=None, kw_only=True, compare=False)

    @classmethod
    def of_pencil(cls, pencil: OperatorPencil[V], sigma: float, source: Any) -> "SourcePosing[V]":
        r"""The affine question :math:`(A - \sigma M)\psi = q` over ``pencil``, carrying its predicted kernel.

        At :math:`\sigma = 0` the operator is :math:`A` and its kernel is the
        pencil's. At :math:`\sigma \neq 0` the pencil's kernel is the
        member's only when :math:`M` respects it (the eigen question's
        condition, :meth:`EigenPosing.solvability`), which is checked here and
        refused otherwise; a kernel the member gains at a spectral coincidence
        of the data is a value question, the multiplying source's
        subcriticality.
        """
        sigma = float(sigma)
        if sigma != 0.0:
            kernel = _predicted(pencil, "SourcePosing.of_pencil")
            if kernel.dimension:
                reading, label = _pencil_annihilation(pencil, kernel)
                if reading > CONSISTENCY_RTOL:
                    raise NoSolution(
                        f"the posing (A - {sigma:g} M) psi = q has no predicted kernel: the pencil's rhs does not "
                        f"respect the loss's kernel ({label!r} reads {reading:.3g}), so the loss's singular class "
                        f"receives fission and the question has no steady solution."
                    )
        return cls(pencil.at(sigma), source, kernel=pencil.kernel)

    def __post_init__(self) -> None:
        codomain = self.operator.codomain
        if codomain is None:
            return
        space = getattr(self.source, "space", None)
        if space is None:
            # A bare source is parsed once, here, at the codomain's exact shape.
            # A typed source keeps the comparison below: a composite source is
            # not a flat member, and its identity leg is #666's.
            codomain.coordinates_of(self.source)
        elif space != codomain:
            raise ValueError(
                f"SourcePosing: the source lives on {space!r} but the operator's "
                f"codomain is {codomain!r} — pose q on the operator's ends."
            )

    def residual(self, psi: Any) -> Any:
        r""":math:`A\psi - q` — the LOSS-sign convention every SN residual uses.

        ⛔ Until step 3 U1 (2026-09-17) this read :math:`q - A\psi`, the sign
        twin of :func:`~orpheus.sn.solver.evaluate_residual` (``[M]`` bit-exactly
        its negative on both arms) while :meth:`balance` already reported
        :math:`\langle w, A\psi - q\rangle`; one convention now (Pattern 7),
        so the exit report's balance member, which reads this residual, and the typed full-system
        residual agree without a sign to remember.
        """
        return self.operator.apply(psi) - self.source

    def balance(self, psi: Any, w: Any = 1.0) -> float:
        r"""β(ψ) = ⟨w, Aψ − q⟩ — the exit balance defect (a signed number, zero at the solution)."""
        return _pair(w, self.residual(psi))

    def H(self, detector: Any) -> "SourcePosing[V]":
        r"""UNARY — the adjoint source problem :math:`A^\dagger\psi^\dagger = R` needs a DETECTOR; the kernel daggers with it."""
        return SourcePosing(self.operator.H, detector, kernel=None if self.kernel is None else self.kernel.H)

    def solvability(self, *, conditioning: Evidence) -> Solvable:
        r"""The Fredholm verdict on :math:`A\psi = q` from the predicted kernel: ``Solve``, ``Gauge``, or raise ``NoSolution``."""
        if self.kernel is None:
            raise ValueError(
                "SourcePosing.solvability: the posing carries no predicted kernel, so it has no Fredholm verdict; "
                "mint it with SourcePosing.of_pencil from a hub's pencil."
            )
        return fredholm_alternative(self.kernel, _flat(self.source), conditioning=conditioning)
