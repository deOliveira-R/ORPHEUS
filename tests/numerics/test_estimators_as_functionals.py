r"""#257 S5 → #259 P1 R8 — the KEigenvalue estimators, hardwired.

History: ``iteration.py`` used to carry two injectable estimator
callables threaded into ``KEigenvalue`` (``keff_estimator`` /
``production_estimator`` kwargs with ``_default_*`` module-function
defaults). #257 S5 examined their category (bare 4-arg callables, NOT
``Functional[V]``s) and kept the lighter-touch layout. **#259 P1 / R8
(2026-07-03) retired the injection seam entirely**: production solvers
implement the ``EigenvalueSolver`` protocol directly by design (the seam
had zero production consumers), and at a converged eigenpair every
estimator CONSISTENT with the posed problem agrees — so injection could
only introduce an inconsistent functional. The default spellings were
folded into the methods, arithmetic unchanged:

* ``KEigenvalue.compute_production_rate(ψ) = float(np.sum(F.apply(ψ)))``
* ``KEigenvalue.compute_keff(ψ)`` = the leakage-inclusive Rayleigh
  quotient ``Σ(Fψ) / (Σ(Aψ) − Σ(Sψ))`` — the operator-form spelling of
  the unified k discipline (when ``A`` carries streaming, the
  denominator IS absorption + leakage − scattering-family gains).

This module is the **bit-identity pin of the surviving spellings**: on a
synthetic L0 operator triple with a numpy ground truth, the hardwired
methods return the EXACT float the documented formulas give. The S5
category leg retired with the callables it examined (a method is plainly
not a mislabeled ``Functional``; there is nothing left to classify).

vv claim layer: regression/value claims (bit-identity against the
literal documented formula — a structurally trivial reference because
the formula is the spec). Zero eigenvalue claims here (KEigenvalue
eigenvalue correctness is pinned by the ``test_iteration.py`` gates).

vv Mode-8: structural assertions via ``_require``; value assertions via
``np.testing.*`` — both fire under ``python -O``.

``foundation`` — software invariants on the estimator surface.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.numerics.iteration import KEigenvalue
from orpheus.numerics.operator import (
    LinearOperator,
    ZeroOperator,
)

pytestmark = pytest.mark.foundation


def _require(condition: bool, message: str) -> None:
    """A -O-firing assertion (NOT a bare ``assert``)."""
    if not condition:
        pytest.fail(message)


# ───────────────────────────────────────────────────────────────────────
# Synthetic L0 operator triple (numpy ground truth) — same shape as the
# fixtures in test_iteration.py, kept self-contained here.
# ───────────────────────────────────────────────────────────────────────


class _MatrixOperator(LinearOperator):
    """Dense-matrix test operator with an honest invertibility pair.

    The estimators consume ``.apply``; ``KEigenvalue`` CONSTRUCTION
    additionally builds ``A.inverse()`` for its inner driver (taxonomy
    step 3 — the posing layer builds the inverse), so the A slot needs
    the ``is_invertible`` + ``inverse()`` pair. Diagonal test matrices
    make the inverse exact and trivial.
    """
    # S4-amendment: the base DEMANDS an answer from every subclass; this
    # double is a deliberately-unbound probe, so it DECLARES the unbound
    # state instead of inheriting a silent default (which no longer exists).
    domain = None
    codomain = None

    def __init__(self, matrix: np.ndarray) -> None:
        self.matrix = np.asarray(matrix, dtype=float)

    @property
    def is_invertible(self) -> bool:
        return True

    def inverse(self) -> "_MatrixOperator":
        return _MatrixOperator(np.linalg.inv(self.matrix))

    def apply(self, x: np.ndarray) -> np.ndarray:
        return self.matrix @ x


def _synthetic_triple():
    """An (A, S, F, ψ) quadruple with a hand-computable estimator value."""
    A = _MatrixOperator(np.diag([3.0, 5.0, 7.0]))
    S = _MatrixOperator(np.diag([0.5, 0.25, 0.1]))
    F = _MatrixOperator(np.diag([2.0, 1.0, 0.5]))
    psi = np.array([1.0, 2.0, 4.0])
    return A, S, F, psi


# ═══════════════════════════════════════════════════════════════════════
# The bit-identity pins — the hardwired methods ARE the documented
# formulas (the R8 fold must not perturb a single bit vs the retired
# module-function defaults, whose arithmetic these references restate).
# ═══════════════════════════════════════════════════════════════════════


class TestHardwiredEstimatorArithmetic:
    """R8 folded the defaults into methods — every bit is preserved."""

    def test_production_rate_bit_identical(self):
        """``compute_production_rate(ψ)`` == ``float(np.sum(F.apply(ψ)))``.

        The reference is the literal documented formula recomputed
        inline. 0 ULP expected — same ``np.sum`` over the same array.
        """
        A, S, F, psi = _synthetic_triple()
        ke = KEigenvalue(A, S, F)
        got = ke.compute_production_rate(psi)
        ref = float(np.sum(F.apply(psi)))  # F = diag(2,1,0.5); ψ=(1,2,4)
        # = 2·1 + 1·2 + 0.5·4 = 6.0
        _require(got == ref, f"production rate {got} != formula {ref}.")
        _require(got == 6.0, f"production rate {got} != hand value 6.0.")

    def test_keff_bit_identical(self):
        """``compute_keff(ψ)`` == ``Σ(Fψ)/(Σ(Aψ)−Σ(Sψ))``.

        Hand value: Σ(Fψ)=6.0; Σ(Aψ)=3+10+28=41; Σ(Sψ)=0.5+0.5+0.4=1.4;
        k = 6.0/(41−1.4) = 6.0/39.6.
        """
        A, S, F, psi = _synthetic_triple()
        ke = KEigenvalue(A, S, F)
        got = ke.compute_keff(psi)
        num = float(np.sum(F.apply(psi)))
        den = float(np.sum(A.apply(psi))) - float(np.sum(S.apply(psi)))
        ref = num / den
        _require(got == ref, f"keff estimator {got} != formula {ref}.")
        np.testing.assert_array_almost_equal_nulp(
            np.array(got), np.array(6.0 / 39.6), nulp=4
        )

    def test_keff_with_zero_S_unchanged(self):
        """The ``S = ZeroOperator`` path (the KEigenvalue default posture)."""
        A, _, F, psi = _synthetic_triple()
        ke = KEigenvalue(A, ZeroOperator(), F)
        got = ke.compute_keff(psi)
        num = float(np.sum(F.apply(psi)))
        den = float(np.sum(A.apply(psi)))  # Σ(Sψ)=0
        _require(got == num / den, f"keff with ZeroOperator S drifted: {got}.")

    def test_injection_kwargs_are_gone(self):
        """R8 teeth: the retired kwargs are NOT silently accepted.

        A ``keff_estimator=`` caller must get a loud ``TypeError`` — the
        seam is retired, not deprecated-and-ignored (a silently swallowed
        kwarg would let an inconsistent-estimator call site believe it
        took effect).
        """
        A, S, F, _ = _synthetic_triple()
        with pytest.raises(TypeError):
            KEigenvalue(A, S, F, keff_estimator=lambda a, s, f, p: 1.0)
        with pytest.raises(TypeError):
            KEigenvalue(A, S, F, production_estimator=lambda a, s, f, p: 1.0)
