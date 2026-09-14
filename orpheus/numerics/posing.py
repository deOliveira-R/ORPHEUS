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
    y_arr = np.asarray(y, dtype=float)
    w_arr = np.broadcast_to(np.asarray(w, dtype=float), y_arr.shape)
    return float(np.vdot(w_arr.ravel(), y_arr.ravel()))


@dataclass(frozen=True)
class SpectralMap:
    r"""The bijection between the pencil's eigenvalue μ and the physical eigenvalue λ."""

    forward: Callable[[float], float]  # μ ↦ λ
    inverse: Callable[[float], float]  # λ ↦ μ
    name: str

    def __call__(self, mu: float) -> float:
        return float(self.forward(mu))


K_MAP = SpectralMap(forward=lambda mu: 1.0 / mu, inverse=lambda k: 1.0 / k, name="k")
r"""k-eigenvalue: :math:`A\psi = F\psi/k` ⟺ :math:`\mu = 1/k`."""

ALPHA_MAP = SpectralMap(forward=lambda mu: -1.0 / mu, inverse=lambda a: -1.0 / a, name="alpha")
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
        num = _pair(w, self.pencil.lhs.apply(psi))
        den = _pair(w, self.pencil.rhs.apply(psi))
        return self.spectral_map(num / den)

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
        r""":math:`q - A\psi`."""
        return self.source - self.operator.apply(psi)

    def balance(self, psi: Any, w: Any = 1.0) -> float:
        r"""β(ψ) = ⟨w, Aψ − q⟩ — the exit balance defect (a signed number, zero at the solution)."""
        return -_pair(w, self.residual(psi))

    def H(self, detector: Any) -> "SourcePosing[V]":
        r"""UNARY — the adjoint source problem :math:`A^\dagger\psi^\dagger = R` needs a DETECTOR."""
        return SourcePosing(self.operator.H, detector)
