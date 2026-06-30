"""Shared faithfulness assertion for the capability→predicate carve (#226).

During the inverse-as-operator carve's Phase 2–3 coexistence, the runtime
:attr:`~orpheus.numerics.operator.LinearOperator.is_invertible` /
:attr:`~orpheus.numerics.operator.LinearOperator.is_adjointable` predicates and
the legacy ``capabilities: frozenset[str]`` advertisement BOTH exist. The
carve's keystone invariant is that they agree EXACTLY for every operator —
``is_X ≡ CAP_X ∈ capabilities`` — which is what licenses deleting the frozenset
in Phase 4 (``operator_inverse_algebra_carve.md`` §4; verification spec §2.3).

This single assertion IS that invariant. Both the numerics-layer enumeration
(``tests/numerics/test_operator_capability_predicates.py``) and the
SN/transport-layer enumeration (``tests/sn/operators/test_capability_survival.py``)
call it, so the keystone contract lives in ONE place (Cardinal Rule 2).

**Why an explicit ``raise``, not ``assert``.** This is a plain helper module,
NOT a pytest-collected test module, so pytest does NOT AST-rewrite its
assertions. Under the canonical ``python -O`` invocation a bare ``assert`` here
would be stripped to a no-op — an inert gate (``vv-principles`` Mode 8). An
explicit ``raise AssertionError`` is a normal statement that fires under ``-O``.
"""
from __future__ import annotations

from orpheus.numerics.operator import (
    CAP_APPLY_TRANSPOSE,
    CAP_SOLVE,
    LinearOperator,
)


def assert_capability_faithful(op: LinearOperator) -> None:
    """Pin ``is_invertible ≡ CAP_SOLVE∈caps`` AND ``is_adjointable ≡ CAP_APPLY_TRANSPOSE∈caps``.

    The coexistence-era equivalence between the derived predicates and the
    legacy frozenset, asserted per operator instance. Raises
    :class:`AssertionError` (``-O``-safe — see module docstring) on any drift.
    """
    caps = op.capabilities
    if op.is_invertible != (CAP_SOLVE in caps):
        raise AssertionError(
            f"{type(op).__name__}: is_invertible={op.is_invertible} but "
            f"(CAP_SOLVE in caps)={CAP_SOLVE in caps}"
        )
    if op.is_adjointable != (CAP_APPLY_TRANSPOSE in caps):
        raise AssertionError(
            f"{type(op).__name__}: is_adjointable={op.is_adjointable} but "
            f"(CAP_APPLY_TRANSPOSE in caps)={CAP_APPLY_TRANSPOSE in caps}"
        )
