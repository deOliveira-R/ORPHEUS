"""The permanent two-axis inverse/adjoint contract (carve P4, spec §36).

The COEXISTENCE-era scaffold that used to live here
(``assert_capability_faithful``, pinning ``is_X ≡ CAP_X ∈ capabilities``)
DELETED with the frozenset at W2 — it could only run while
``op.capabilities`` existed, and its job (licensing the deletion) is
done. What remains is its permanent successor: the direct two-axis
contract, referencing NO frozenset, shared by the numerics enumeration
(``tests/numerics/test_operator_capability_predicates.py``) and the
SN/transport enumeration (``tests/sn/operators/test_capability_survival.py``)
so the keystone contract lives in ONE place (Cardinal Rule 2).

**Why explicit ``raise``, not ``assert``.** This is a plain helper module,
NOT a pytest-collected test module, so pytest does NOT AST-rewrite its
assertions. Under the canonical ``python -O`` invocation a bare ``assert``
here would be stripped to a no-op — an inert gate (``vv-principles``
Mode 8). An explicit ``raise AssertionError`` is a normal statement that
fires under ``-O``.
"""
from __future__ import annotations

# ───────────────────────────────────────────────────────────────────────
# Keystone v2 (carve P4, spec §36) — the standing faithfulness net.
#
# The inverse axis has a THREE-VALUED contract (the structural/value split
# IS part of the contract, Design C):

#: No inverse exists for the TYPE, mathematically — ``inverse()`` is not
#: declared; misuse is a STATIC error (``AttributeError`` only if forced).
STRUCTURAL_ABSENT = "structural-absent"
#: The type supports inversion; THIS instance refuses — ``inverse()`` is
#: declared and raises ``NotInvertible`` eagerly.
VALUE_RAISE = "value-raise"
#: ``inverse()`` returns an inverse operator (the I1 round-trip gates its
#: VALUE elsewhere — this contract pins returns-vs-raises-vs-absent only).
INVERTIBLE = "invertible"


def assert_inverse_adjoint_contract(
    op, *, invertible, adjointable, inverse_contract
):
    """Pin the two-axis inverse/adjoint contract on one operator (spec §36).

    invertible : bool
        Expected ``op.is_invertible``.
    adjointable : bool
        Expected ``op.is_adjointable`` — and the eager-``.H`` behaviour
        that must match it (returns the wrapper vs raises
        ``MissingAdjoint`` at CONSTRUCTION, never lazily at apply).
    inverse_contract : str
        One of :data:`INVERTIBLE` / :data:`VALUE_RAISE` /
        :data:`STRUCTURAL_ABSENT`.

    Explicit ``raise`` throughout (``-O``-safe — un-collected helper
    module, vv-principles Mode 8; see the module docstring).
    """
    from orpheus.numerics.operator import (
        MissingAdjoint,
        NotInvertible,
        adjointable as _adj_bridge,
        invertible as _inv_bridge,
    )

    # (a) inverse axis — predicate + the three-valued method contract.
    if op.is_invertible != invertible:
        raise AssertionError(
            f"{op!r}: is_invertible={op.is_invertible}, expected {invertible}"
        )
    if invertible:
        if inverse_contract != INVERTIBLE:
            raise AssertionError(
                f"{op!r}: is_invertible=True demands the INVERTIBLE contract"
            )
        op.inverse()  # MUST return; the I1 round-trip gates the value
    elif inverse_contract == VALUE_RAISE:
        if not hasattr(op, "inverse"):
            raise AssertionError(
                f"{op!r}: value-dependent type must DECLARE inverse()"
            )
        try:
            op.inverse()
        except NotInvertible:
            pass
        else:
            raise AssertionError(
                f"{op!r}: singular inverse() did not raise NotInvertible"
            )
    elif inverse_contract == STRUCTURAL_ABSENT:
        if hasattr(op, "inverse"):
            raise AssertionError(
                f"{op!r}: structurally non-invertible type must NOT "
                f"declare inverse()"
            )
    else:
        raise AssertionError(f"unknown inverse_contract {inverse_contract!r}")
    # (b) adjoint axis — uniformly eager return-or-raise (.H is base-level).
    if op.is_adjointable != adjointable:
        raise AssertionError(
            f"{op!r}: is_adjointable={op.is_adjointable}, expected {adjointable}"
        )
    if adjointable:
        op.H  # eager: returns the wrapper, no raise
    else:
        try:
            op.H
        except MissingAdjoint:
            pass
        else:
            raise AssertionError(
                f"{op!r}: non-adjointable .H did not raise MissingAdjoint "
                f"EAGERLY at construction"
            )
    # (c) bridge consistency — pins the one-line TypeGuard bodies (spec
    # §36 leg c; a bridge reading anything but the predicate drifts here).
    if _inv_bridge(op) != op.is_invertible:
        raise AssertionError(
            f"{op!r}: invertible() bridge drifted from is_invertible"
        )
    if _adj_bridge(op) != op.is_adjointable:
        raise AssertionError(
            f"{op!r}: adjointable() bridge drifted from is_adjointable"
        )
