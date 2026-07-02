r"""#257 S5 — Spec C: the keff / production estimators as ``Functional``s.

``iteration.py`` carries two estimators threaded into ``KEigenvalue``:

* ``_default_production_estimator(L,S,F,ψ) = float(np.sum(F.apply(ψ)))``
* ``_default_keff_estimator(L,S,F,ψ)`` = the Rayleigh quotient
  ``Σ(Fψ) / (Σ(Lψ) − Σ(Sψ))``.

S5 narrows their TYPE contract toward the ``Functional`` category; the
ARITHMETIC is UNCHANGED (bit-identical).

Two legs:

1. **Unchanged-arithmetic invariant (the load-bearing leg).** On a
   synthetic L0 operator triple with a numpy ground truth, the
   estimators return the EXACT same float they returned before S5.
   This is a regression pin: S5 must not perturb a single bit of the
   estimator output.

2. **Category-honesty leg (tolerant of the chosen layout).** The brief
   asks whether the bare ``Callable`` aliases are kept or wrapped in a
   ``Functional`` object, and recommends the lighter-touch option that
   still makes the category honest. RECOMMENDATION (see the module-level
   note): keep the bare ``(L,S,F,ψ)->float`` callables — they are NOT
   ``Functional[V]`` (whose surface is ``evaluate(x: V) -> float | V``,
   a ONE-argument field→scalar map; the estimators take the operator
   triple too). The honest ``Functional`` is the field→scalar CORE
   ``ψ -> Σ(Fψ)`` once ``F`` is bound. IF the method-implementer ships
   such a bound wrapper, this leg asserts it satisfies the Protocol AND
   reproduces the bare callable's float. IF no wrapper is shipped (the
   lighter touch), the leg SKIPS with a clear reason — the category is
   still honest because the estimators were never claimed to be
   ``Functional``s; they remain plain renormalisation callables.

vv claim layer (1.5 gate): leg 1 is a regression/value claim
(bit-identity inherited from the pre-S5 arithmetic — a structurally
trivial reference because the formula is unchanged); leg 2 is a
CATEGORY claim (Protocol membership). Zero eigenvalue claims here (the
KEigenvalue eigenvalue correctness is pinned by the EXISTING
``test_iteration.py`` gates, listed in the gate list).

vv Mode-8: structural assertions via ``require``; value assertions via
``np.testing.*`` — both fire under ``python -O``.

``foundation`` — software invariants on the estimator surface.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.numerics.iteration import (
    _default_keff_estimator,
    _default_production_estimator,
)
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
    """Dense-matrix apply-only test operator.

    Inherits the base ``False`` predicates — the estimators consume
    ONLY ``.apply``, so no invertibility/adjointability advertisement
    is needed (the pre-carve ``can_solve`` capability flag was pure
    advertisement over a class with no ``solve`` body; retired with
    the frozenset at carve P4).
    """

    def __init__(self, matrix: np.ndarray) -> None:
        self.matrix = np.asarray(matrix, dtype=float)

    def apply(self, x: np.ndarray) -> np.ndarray:
        return self.matrix @ x


def _synthetic_triple():
    """An (L, S, F, ψ) quadruple with a hand-computable estimator value."""
    L = _MatrixOperator(np.diag([3.0, 5.0, 7.0]))
    S = _MatrixOperator(np.diag([0.5, 0.25, 0.1]))
    F = _MatrixOperator(np.diag([2.0, 1.0, 0.5]))
    psi = np.array([1.0, 2.0, 4.0])
    return L, S, F, psi


# ═══════════════════════════════════════════════════════════════════════
# C.1 — Unchanged-arithmetic invariant (bit-identity, the regression pin).
# ═══════════════════════════════════════════════════════════════════════


class TestEstimatorArithmeticUnchanged:
    """S5 narrows the type, not the arithmetic — every bit is preserved."""

    def test_production_estimator_value_bit_identical(self):
        """``_default_production_estimator`` == ``float(np.sum(F.apply(ψ)))``.

        The reference is the literal documented formula recomputed inline
        (the arithmetic S5 must not touch). 0 ULP expected — same
        ``np.sum`` over the same array.
        """
        L, S, F, psi = _synthetic_triple()
        got = _default_production_estimator(L, S, F, psi)
        ref = float(np.sum(F.apply(psi)))  # F = diag(2,1,0.5); ψ=(1,2,4)
        # = 2·1 + 1·2 + 0.5·4 = 6.0
        _require(got == ref, f"production estimator {got} != formula {ref}.")
        _require(got == 6.0, f"production estimator {got} != hand value 6.0.")

    def test_keff_estimator_value_bit_identical(self):
        """``_default_keff_estimator`` == ``Σ(Fψ)/(Σ(Lψ)−Σ(Sψ))``.

        Hand value: Σ(Fψ)=6.0; Σ(Lψ)=3+10+28=41; Σ(Sψ)=0.5+0.5+0.4=1.4;
        k = 6.0/(41−1.4) = 6.0/39.6.
        """
        L, S, F, psi = _synthetic_triple()
        got = _default_keff_estimator(L, S, F, psi)
        num = float(np.sum(F.apply(psi)))
        den = float(np.sum(L.apply(psi))) - float(np.sum(S.apply(psi)))
        ref = num / den
        _require(got == ref, f"keff estimator {got} != formula {ref}.")
        np.testing.assert_array_almost_equal_nulp(
            np.array(got), np.array(6.0 / 39.6), nulp=4
        )

    def test_production_estimator_with_zero_S_unchanged(self):
        """The ``S = ZeroOperator`` path (the KEigenvalue default) is unchanged."""
        L, _, F, psi = _synthetic_triple()
        S = ZeroOperator()
        got = _default_keff_estimator(L, S, F, psi)
        num = float(np.sum(F.apply(psi)))
        den = float(np.sum(L.apply(psi)))  # Σ(Sψ)=0
        _require(got == num / den, f"keff with ZeroOperator S drifted: {got}.")


# ═══════════════════════════════════════════════════════════════════════
# C.2 — Category-honesty leg (tolerant of the chosen layout).
# ═══════════════════════════════════════════════════════════════════════


class TestEstimatorCategoryHonesty:
    """Make the category honest WITHOUT over-coupling to a layout choice.

    The recommended lighter-touch design keeps the bare ``(L,S,F,ψ)``
    callables (they are not field→scalar ``Functional[V]``s). This leg
    only fires if the method-implementer ALSO shipped a bound
    field→scalar wrapper; otherwise it skips, documenting that the
    category remains honest because the estimators never claimed to be
    Functionals.
    """

    def test_bare_callables_are_not_functionals(self):
        """A bare ``(L,S,F,ψ)`` callable does NOT satisfy ``Functional[V]``.

        ``Functional.evaluate`` is a ONE-argument field→scalar map; the
        estimators take the operator triple. So a runtime ``isinstance``
        against the Protocol must be False — confirming we did NOT
        mislabel a 4-arg renormalisation callable as a Functional.
        Skips if ``Functional`` is not yet on the tree.
        """
        functional_mod = pytest.importorskip(
            "orpheus.numerics.functional",
            reason="#257 S5 PRE-IMPL: Functional Protocol not yet written.",
        )
        Functional = getattr(functional_mod, "Functional", None)
        if Functional is None:
            pytest.skip("#257 S5 PRE-IMPL: `Functional` symbol not defined.")
        # The bare module-level callables are plain functions — they
        # carry no ``evaluate`` method, so they are not Functionals.
        _require(
            not isinstance(_default_production_estimator, Functional),
            "A bare (L,S,F,ψ) estimator callable must NOT be classified as "
            "a Functional[V] (whose surface is the 1-arg evaluate(x)->scalar).",
        )

    def test_bound_production_functional_if_shipped(self):
        """IF a bound ``ψ->Σ(Fψ)`` Functional wrapper exists, it is honest.

        The recommended optional wrapper binds ``F`` so the surface
        becomes the genuine field→scalar map ``evaluate(ψ) = Σ(Fψ)``.
        If shipped (probed below), assert it (a) satisfies the Functional
        Protocol AND (b) reproduces the bare callable's float bit-for-bit.
        If NOT shipped, skip — the lighter-touch design is the accepted
        recommendation and leaves nothing to assert here.
        """
        functional_mod = pytest.importorskip(
            "orpheus.numerics.functional",
            reason="#257 S5 PRE-IMPL: Functional Protocol not yet written.",
        )
        Functional = getattr(functional_mod, "Functional", None)
        wrapper_cls = getattr(functional_mod, "ProductionFunctional", None)
        if Functional is None or wrapper_cls is None:
            pytest.skip(
                "Optional bound ProductionFunctional wrapper not shipped — "
                "the recommended lighter-touch design keeps bare callables; "
                "nothing to assert."
            )
        L, S, F, psi = _synthetic_triple()
        wrapper = wrapper_cls(F)
        _require(
            isinstance(wrapper, Functional),
            "A bound ProductionFunctional must satisfy the Functional Protocol.",
        )
        bound_value = wrapper.evaluate(psi)
        bare_value = _default_production_estimator(L, S, F, psi)
        _require(
            float(bound_value) == bare_value,
            f"Bound Functional value {bound_value} must bit-match the bare "
            f"callable {bare_value}.",
        )
