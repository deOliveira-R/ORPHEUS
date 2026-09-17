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

from dataclasses import dataclass
from typing import Any, Callable, Generic, TypeVar

import numpy as np

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

ALPHA_MAP = SpectralMap(forward=lambda mu: -1.0 / mu, inverse=lambda a: -1.0 / a, of_quotient=lambda a, m: -m / a, name="alpha")
r"""α-eigenvalue: :math:`(L+C-S-F)\psi = -\alpha T\psi` ⟺ :math:`\mu = -1/\alpha` (a later posing; the map is stated now)."""


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


@dataclass(frozen=True)
class SourcePosing(Generic[V]):
    r"""The affine question :math:`A\psi = q` over one operator and one source."""

    operator: LinearOperator[V, V]
    source: Any

    def __post_init__(self) -> None:
        space = getattr(self.source, "space", None)
        if space is not None and self.operator.codomain is not None and space != self.operator.codomain:
            raise ValueError(
                f"SourcePosing: the source lives on {space!r} but the operator's "
                f"codomain is {self.operator.codomain!r} — pose q on the operator's ends."
            )

    def residual(self, psi: Any) -> Any:
        r""":math:`A\psi - q` — the LOSS-sign convention every SN residual uses.

        ⛔ Until step 3 U1 (2026-09-17) this read :math:`q - A\psi`, the sign
        twin of :func:`~orpheus.sn.solver.evaluate_residual` (``[M]`` bit-exactly
        its negative on both arms) while :meth:`balance` already reported
        :math:`\langle w, A\psi - q\rangle`; one convention now (Pattern 7),
        so the certificate that reads this residual and the typed full-system
        residual agree without a sign to remember.
        """
        return self.operator.apply(psi) - self.source

    def balance(self, psi: Any, w: Any = 1.0) -> float:
        r"""β(ψ) = ⟨w, Aψ − q⟩ — the exit balance defect (a signed number, zero at the solution)."""
        return _pair(w, self.residual(psi))

    def H(self, detector: Any) -> "SourcePosing[V]":
        r"""UNARY — the adjoint source problem :math:`A^\dagger\psi^\dagger = R` needs a DETECTOR."""
        return SourcePosing(self.operator.H, detector)
