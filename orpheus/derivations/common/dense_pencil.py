r"""The reference kernel's dense eigen and source solves, in weak form.

A continuous reference that has reduced its transport problem to a few
hundred unknowns asks two questions of the resulting matrices: the
eigenpairs of the pencil :math:`F c = k\,L c` (the fundamental mode, the
higher modes, the adjoint modes), and the solution of a source problem
:math:`M x = G x + s`, which is the pencil :math:`(M, G)` read as a source
problem. This module answers both with numpy and scipy only.

**Why the references carry their own.** Production's
:mod:`orpheus.numerics.eigenvalue` and :mod:`orpheus.numerics.pencil` are
built for versatility (any operator, any space, the operator algebra's
``.H``) and change as production grows. A reference is closed, purpose
built and limited in scope, so its machinery is small and changes only by
ruling, and production's churn never reaches it (the user's ruling of
2026-10-06, recorded in ``docs/architecture/layering.rst``). The
Perron--Frobenius refusal therefore exists twice, once on each side of the
branch line, by design: each copy is verified on its own, so agreement
between a reference and a production solver means something.

**Weak form.** Both matrices of a :class:`DensePencil` are bilinear forms
on one basis, :math:`L_{ij} = a_L(u_i, u_j)` and
:math:`F_{ij} = a_F(u_i, u_j)`, as a Galerkin assembly produces them. The
adjoint operator's bilinear forms are then the transposes, so
:meth:`DensePencil.adjoint` needs no metric. A nodal basis makes the
coefficients function values, which is what lets
:meth:`DensePencil.fundamental` read the single sign of the mode off its
coefficients.

**The refusals.** For the transport operator of a well-posed criticality
problem, the Krein--Rutman theorem (the infinite-dimensional
Perron--Frobenius theorem for a positive compact operator) makes the
fundamental eigenvalue real, positive, simple and strictly dominant in
modulus, with a single-signed eigenfunction. A discretisation that loses
any of these has not produced a fundamental mode, and
:meth:`DensePencil.fundamental` refuses instead of returning one
(:class:`NoFundamentalMode`). The higher modes carry no such contract and
:meth:`DensePencil.spectrum` returns them as they are, complex or
sign-changing.

The source solution is the sum of the Neumann series
:math:`x = \sum_n (M^{-1}G)^n M^{-1} s`, which for a positive gain is the
least non-negative solution; it exists for every source iff the spectral
radius :math:`\rho(M^{-1}G) < 1` on the unknowns the source reaches, and
:meth:`DensePencil.least_solution` refuses
otherwise (:class:`NoLeastSolution`), because a direct solve of a
supercritical system returns a finite, negative flux without complaint.
A part of the system the source never reaches (in the Frobenius normal form,
a class downstream of no source) carries zero whatever its spectral radius.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import scipy.linalg

#: The safety factor on the fundamental eigenvector's first-order error bound,
#: the band within which a negative coefficient counts as rounding of an exact
#: zero (a group with no source, a node on a vacuum wall) rather than a sign
#: change. The bound itself is the pencil's, :meth:`DensePencil.sign_band`: a
#: fixed band refused 9 of 20 pencils with cond(L) near 1e8, whose exact zeros
#: return at up to 5.6e-9 (``[M]`` qa, 2026-10-06).
SIGN_SAFETY: float = 10.0

#: Smallest relative gap in modulus between the fundamental eigenvalue and the
#: next one for the fundamental to count as strictly dominant.
DOMINANCE_GAP: float = 1e-10

#: Smallest distance below 1 of the gain's spectral radius for a source problem
#: to count as subcritical. The radius is computed, so an exactly critical
#: system reads 1 to within the eigensolver's backward error (``[M]``
#: 2026-10-06, 4-by-4 gains of radius exactly 1, 20 seeds of 200 draws: the
#: computed radius spans 1 - 7.2e-13 to 1 + 1.8e-13, and a bare comparison
#: with 1 admitted 473 of 4000, solving a singular system); and a
#: system nearer to critical than this has a least solution amplified by more
#: than ``1 / SUBCRITICAL_MARGIN``, which no reference value survives.
SUBCRITICAL_MARGIN: float = 1e-10


class NoFundamentalMode(ValueError):
    """The pencil has no fundamental mode in the Krein--Rutman sense."""


class NoLeastSolution(ValueError):
    """The source problem's Neumann series diverges: the system is not subcritical."""


def _frozen(array: np.ndarray) -> np.ndarray:
    a = np.array(array)
    a.flags.writeable = False
    return a


@dataclass(frozen=True)
class FundamentalMode:
    r"""The fundamental mode: :math:`k` real and positive, ``vector`` real, non-negative, unit norm.

    The invariant is the Krein--Rutman contract, checked at construction, so
    a value of this type is a fundamental mode whoever built it.
    """

    k: float
    vector: np.ndarray

    def __post_init__(self) -> None:
        vector = _frozen(np.asarray(self.vector, dtype=float))
        if not (np.isfinite(self.k) and self.k > 0.0):
            raise NoFundamentalMode(f"a fundamental eigenvalue is finite and positive; got {self.k!r}")
        if vector.ndim != 1 or not np.all(np.isfinite(vector)) or vector.min() < 0.0:
            raise NoFundamentalMode("a fundamental mode's vector is finite and non-negative")
        if abs(np.linalg.norm(vector) - 1.0) > 8 * np.finfo(float).eps * np.sqrt(vector.size):
            raise NoFundamentalMode("a fundamental mode's vector has unit Euclidean norm")
        object.__setattr__(self, "vector", vector)


@dataclass(frozen=True)
class SpectralMode:
    r"""One eigenpair of the pencil, complex in general; ``vector`` has unit Euclidean norm."""

    eigenvalue: complex
    vector: np.ndarray

    def __post_init__(self) -> None:
        object.__setattr__(self, "vector", _frozen(self.vector))


def _real(name: str, values: np.ndarray) -> np.ndarray:
    if np.iscomplexobj(values):
        raise ValueError(f"{name}: a real array is required; got a complex one")
    return np.array(values, dtype=float)


def _square_finite(name: str, matrix: np.ndarray) -> np.ndarray:
    a = _real(name, matrix)
    if a.ndim != 2 or a.shape[0] != a.shape[1]:
        raise ValueError(f"{name}: a square matrix is required; got shape {a.shape}")
    if not np.all(np.isfinite(a)):
        raise ValueError(f"{name}: every entry must be finite")
    a.flags.writeable = False
    return a


@dataclass(frozen=True)
class DensePencil:
    r"""The pencil :math:`F c = k\,L c`, both forms in weak form on one basis.

    Attributes
    ----------
    loss : (n, n) ndarray
        The loss form :math:`L`.
    production : (n, n) ndarray
        The production form :math:`F`.
    """

    loss: np.ndarray
    production: np.ndarray

    def __post_init__(self) -> None:
        loss = _square_finite("DensePencil.loss", self.loss)
        production = _square_finite("DensePencil.production", self.production)
        if loss.shape != production.shape:
            raise ValueError(
                f"DensePencil: loss {loss.shape} and production {production.shape} differ in shape"
            )
        object.__setattr__(self, "loss", loss)
        object.__setattr__(self, "production", production)

    def spectrum(self) -> tuple[SpectralMode, ...]:
        r"""Every eigenpair, by decreasing :math:`|k|`, complex and sign-changing modes included.

        Refuses a singular loss form, which shows as an infinite or
        undefined eigenvalue of the QZ reduction.
        """
        values, vectors = scipy.linalg.eig(self.production, self.loss)
        if not np.all(np.isfinite(values)):
            raise ValueError(
                "DensePencil.spectrum: the loss form is singular "
                f"(the QZ reduction returned {values[~np.isfinite(values)][:3]})"
            )
        order = np.argsort(-np.abs(values), kind="stable")
        return tuple(
            SpectralMode(complex(values[i]), vectors[:, i] / np.linalg.norm(vectors[:, i])) for i in order
        )

    def fundamental(self) -> FundamentalMode:
        r"""The fundamental mode, or :class:`NoFundamentalMode` naming the broken condition.

        The eigenvalue of largest modulus must be real (the QZ reduction
        returns a real eigenvalue with an exactly zero imaginary part),
        positive and strictly dominant, and its unit eigenvector
        single-signed within :meth:`sign_band`.
        """
        modes = self.spectrum()
        lead = modes[0]
        if lead.eigenvalue.imag != 0.0:
            raise NoFundamentalMode(
                f"the dominant eigenvalue {lead.eigenvalue:.6g} is complex"
            )
        if lead.eigenvalue.real <= 0.0:
            raise NoFundamentalMode(
                f"the dominant eigenvalue {lead.eigenvalue.real:.6g} is not positive"
            )
        if len(modes) > 1:
            gap = 1.0 - abs(modes[1].eigenvalue) / abs(lead.eigenvalue)
            if gap <= DOMINANCE_GAP:
                raise NoFundamentalMode(
                    f"the dominant eigenvalue {lead.eigenvalue.real:.6g} is not strictly dominant "
                    f"(the next has modulus {abs(modes[1].eigenvalue):.6g})"
                )
        vector = lead.vector.real
        vector = vector if vector.sum() >= 0.0 else -vector
        band = self.sign_band(lead.eigenvalue.real, modes[1].eigenvalue if len(modes) > 1 else 0.0)
        if vector.min() < -band:
            raise NoFundamentalMode(
                f"the dominant eigenvector is not single-signed "
                f"(most negative coefficient {vector.min():.3g} of a unit vector, beyond the "
                f"rounding band {band:.3g})"
            )
        # Within the band a negative coefficient is rounding of an exact zero.
        vector = np.maximum(vector, 0.0)
        return FundamentalMode(float(lead.eigenvalue.real), vector / np.linalg.norm(vector))

    def sign_band(self, k: float, next_eigenvalue: complex) -> float:
        r"""The rounding band on the unit fundamental eigenvector's coefficients.

        :data:`SIGN_SAFETY` times the first-order error bound of the
        eigenvector of :math:`A = L^{-1}F` under the QZ reduction's backward
        error, a relative :math:`\varepsilon` in each form:
        :math:`\|\delta A\| \le \varepsilon(\|L^{-1}\|\,\|F\| + \kappa(L)\,\|A\|)`,
        divided by the eigenvalue's separation from the rest of the spectrum,
        :math:`|k| - |k_2|` (Golub and Van Loan, *Matrix Computations*, 4th
        ed., the eigenvector perturbation bound of a simple eigenvalue). It
        grows with the conditioning of the loss form and shrinks with the
        dominance gap, so an exact zero is admitted and a real sign change is
        refused at every conditioning.
        """
        singular = np.linalg.svd(self.loss, compute_uv=False)
        loss_inverse_norm = 1.0 / singular[-1]
        resolvent_norm = np.linalg.norm(np.linalg.solve(self.loss, self.production), 2)
        backward = np.finfo(float).eps * (
            loss_inverse_norm * np.linalg.norm(self.production, 2)
            + singular[0] * loss_inverse_norm * resolvent_norm
        )
        return SIGN_SAFETY * backward / (abs(k) - abs(next_eigenvalue))

    def adjoint(self) -> DensePencil:
        r"""The adjoint pencil: in weak form, the transposed bilinear forms."""
        return DensePencil(self.loss.T, self.production.T)

    def reach(self, source: np.ndarray) -> np.ndarray:
        r"""The unknowns a source reaches through the pencil's couplings, as a boolean mask.

        Unknown :math:`j` enters equation :math:`i` when :math:`L_{ij}` or
        :math:`F_{ij}` is non-zero. The mask is the closure of the source's
        support under that relation: the union of the source's classes in the
        Frobenius normal form of the coupling and every class downstream of
        them. Nothing outside it can carry a non-zero least solution.
        """
        couples = (self.loss != 0.0) | (self.production != 0.0)
        reached = np.asarray(source) != 0.0
        frontier = reached.copy()
        while frontier.any():
            frontier = couples[:, frontier].any(axis=1) & ~reached
            reached |= frontier
        return reached

    def least_solution(self, source: np.ndarray) -> np.ndarray:
        r"""The solution of :math:`L x = F x + s` as the sum of its Neumann series.

        Read as a source problem, the pencil's loss is the form the unknown
        is measured in and its production the secondary production (the
        scattering and the fission). For a positive production this is the
        least non-negative solution.

        The series is zero outside :meth:`reach`, the unknowns the source
        reaches: those rows and columns decouple, block triangular in the
        Frobenius normal form, with a zero source. On the reached block it
        converges iff the block's spectral radius
        :math:`\rho(L_{RR}^{-1} F_{RR})` is below 1; a radius not below
        :math:`1 -` :data:`SUBCRITICAL_MARGIN` raises :class:`NoLeastSolution`
        naming it. A class with radius at least 1 that the source never
        reaches (a lossless trapped line with no source) is therefore not a
        refusal, and a zero source reaches nothing and gives zero.
        """
        s = _real("DensePencil.least_solution.source", source)
        if s.shape != self.loss.shape[:1] or not np.all(np.isfinite(s)):
            raise ValueError(
                f"DensePencil.least_solution: the source must be finite with shape "
                f"{self.loss.shape[:1]}; got {s.shape}"
            )
        reached = self.reach(s)
        x = np.zeros_like(s)
        if not reached.any():
            return x
        block = np.ix_(reached, reached)
        reached_pencil = DensePencil(self.loss[block], self.production[block])
        radius = abs(reached_pencil.spectrum()[0].eigenvalue)
        if not radius < 1.0 - SUBCRITICAL_MARGIN:
            raise NoLeastSolution(
                f"the Neumann series diverges or is not resolvable: the spectral radius of the gain "
                f"on the {int(reached.sum())} unknowns the source reaches is {radius:.12g}, "
                f"not below 1 - {SUBCRITICAL_MARGIN:g}"
            )
        x[reached] = np.linalg.solve(reached_pencil.loss - reached_pencil.production, s[reached])
        return x


__all__ = [
    "DOMINANCE_GAP",
    "SIGN_SAFETY",
    "SUBCRITICAL_MARGIN",
    "DensePencil",
    "FundamentalMode",
    "SpectralMode",
    "NoFundamentalMode",
    "NoLeastSolution",
]
