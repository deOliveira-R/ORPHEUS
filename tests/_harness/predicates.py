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

import dataclasses
import inspect
import typing
from collections.abc import Callable

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


# ───────────────────────────────────────────────────────────────────────
# #340 N4 — the caller-facing knob contract, shared by every entry point.
#
# ``IterationRecord.budget.name`` is the ONE fact in a convergence warning
# that the record cannot derive for itself, so it is the one a producer can
# silently forget — and the default (``"max_iter"``) is a plausible-looking
# string that is a parameter of NO public ORPHEUS entry.  A forgotten stamp
# ships advice telling the reader to set something that does not exist, and
# no value gate anywhere can see it.
#
# The reference is ``inspect.signature`` of the entry itself — the API,
# structurally independent of the record under test.  Lives here so the SN
# enumeration and the CP/MoC/diffusion enumeration share ONE definition of
# "a knob the caller can reach" (Cardinal Rule 2).

def reachable_knobs(entry: "Callable[..., object]") -> frozenset[str]:
    """Every knob name a caller of ``entry`` can actually set.

    A parameter of ``entry`` is reachable under its own name.  A parameter
    whose annotation is a **dataclass** is a knob BUNDLE, so each of its
    fields is reachable as ``"<param>.<field>"``.

    That second arm is not hypothetical tidiness — it is how CP ships:
    ``solve_cp(materials, mesh, params: CPParams | None)`` exposes no
    ``max_outer`` of its own, because the knob is a field of
    :class:`~orpheus.cp.solver.CPParams`.  A reference built from
    ``signature().parameters`` alone would refuse CP's *correct* stamp, so
    the naive version of this check is not merely incomplete: it is wrong in
    the direction that punishes the honest producer.

    ``X | None`` is unwrapped, since an optional bundle is still reachable.
    """
    names: set[str] = set()
    try:
        hints = typing.get_type_hints(entry)
    except Exception:  # unresolvable forward ref — parameters still count
        hints = {}
    for param in inspect.signature(entry).parameters:
        names.add(param)
        for candidate in _annotation_members(hints.get(param)):
            if dataclasses.is_dataclass(candidate):
                names.update(
                    f"{param}.{f.name}" for f in dataclasses.fields(candidate)
                )
    return frozenset(names)


def _annotation_members(annotation: object) -> tuple[object, ...]:
    """``X`` -> ``(X,)``; ``X | None`` / ``Optional[X]`` -> ``(X, NoneType)``."""
    args = typing.get_args(annotation)
    return args if args else (annotation,)


def assert_every_budgeted_level_names_a_reachable_knob(
    record, entry: "Callable[..., object]",
) -> None:
    """Every level that can EXHAUST advises a knob ``entry`` really exposes.

    Levels with ``budget == 0`` are skipped, and the reason is a property of
    the primitive rather than a convenience: ``exhausted_budget`` is
    ``budget > 0 and ...``, so a budget-free level can never be ``TRUNCATED``
    and its ``budget_name`` can never reach a reader.  Diffusion's
    ``DIRECT`` inner (one LU back-substitution) is exactly that case — an
    exact solve has no knob to advise, and demanding it name one would be
    asserting a string with no consequence.
    """
    knobs = reachable_knobs(entry)
    for level in record.walk():
        if not level.budget.is_budgeted:
            continue
        if level.budget.name not in knobs:
            raise AssertionError(
                f"{getattr(entry, '__name__', entry)} returned level "
                f"{level.label!r} advising `set {level.budget.name}=...`, "
                f"which is not a knob its caller can reach. Reachable: "
                f"{sorted(knobs)}"
            )
