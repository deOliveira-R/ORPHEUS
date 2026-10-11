r"""The rules a cross-check's tolerance is computed by, from a reference's error and the system under test's.

The ladders hold measured rungs; these functions turn them into errors and
tolerances. They read no reference, so any reference's ladder uses them:
:mod:`tests.gates.derivations._characteristic_ladders` (the SN rows' reference
since P1 step (d) of ``.claude/plans/characteristic_reference_architecture.md``)
and the cross-method tolerances (``tests/gates/cross_method/cases.py``). They
lived apart from the ladders while the old trajectory-resolvent family's
ladders used them too (deleted with the family in step (e2)); the user's
ruling of 2026-10-10 keeps them here while the characteristic ladders use them.

Pinned by ``tests/gates/sn/verification/analytical/test_crosscheck_harness.py``.
"""

from __future__ import annotations

import math


def richardson_error(step: float, refinement_ratio: float, order: float) -> float:
    r"""The coarse rung's error from one step at a known order.

    With :math:`e(n) = C n^{-p}` and the step :math:`s = e(n) - e(\rho n)`,
    :math:`e(n) = s / (1 - \rho^{-p})`.
    """
    return abs(step) / (1.0 - refinement_ratio ** (-order))


def alternating_error(step: float) -> float:
    """The coarse rung's error when the sequence alternates with shrinking steps: the limit lies within one step."""
    return abs(step)


def geometric_error(step: float, previous_step: float) -> float:
    r"""The coarse rung's error when the steps shrink geometrically, from its step up and the step below it.

    With steps :math:`s_{n+j} = s_n r^j` (exponential convergence, an
    hp-ladder's), the rung's error is the sum of every step above it,
    :math:`|s_n| / (1 - r)`, with :math:`r = |s_n| / |s_{n-1}|` read off
    the two measured steps. The model needs :math:`r < 1`; a ladder that does
    not contract has no error estimate here, and is refused.
    """
    ratio = abs(step) / abs(previous_step)
    if not ratio < 1.0:
        raise ValueError(f"the steps do not contract (ratio {ratio:.3g}): no geometric error estimate")
    return abs(step) / (1.0 - ratio)


def ceil_one_significant_figure(x: float) -> float:
    """The smallest one-significant-figure value at least ``x`` (1.5e-3 -> 2e-3, 3.2e-3 -> 4e-3)."""
    magnitude = 10.0 ** math.floor(math.log10(x))
    return math.ceil(round(x / magnitude, 9)) * magnitude


def tolerance_for(sut_residual: float, reference_error: float | None) -> float:
    r"""The tolerance a comparison is held to.

    The smallest one-significant-figure value :math:`T` with
    :math:`T \ge 10\,b` (the verification floor) and
    :math:`T \ge 2\,(e + b)` (room for both errors), for the SUT's residual
    :math:`e` and the reference's error :math:`b`: a derived bound once the
    family is certified, a ladder ESTIMATE until then (every reference today,
    whose rows are uncertified comparisons, not verification). With no figure
    for the reference it is assumed at the floor, :math:`b = T/10`, which gives
    :math:`T \ge 2.5\,e`: the tolerance the row will hold once a reference
    certifies a tenth of it.
    """
    if reference_error is None:
        return ceil_one_significant_figure(2.5 * sut_residual)
    return ceil_one_significant_figure(
        max(10.0 * reference_error, 2.0 * (sut_residual + reference_error))
    )


def summed_tolerance(sut_error: float, reference_error: float) -> float:
    r"""The band of a comparison whose two errors are each measured against a common limit: :math:`e + b`, rounded up.

    The ERR-094 partial-reflector rows' rule (the user's ruling 2 of
    2026-10-09: they keep it). It has no factor 2 and no :math:`10\,b`
    floor, so it is tighter than :func:`tolerance_for` wherever both apply,
    and it holds only while each error is a MEASURED distance to the limit.
    """
    return ceil_one_significant_figure(sut_error + reference_error)


def sn_residual(steps: dict[str, float]) -> float:
    """The SUT's residual at its fixture: the sum of its per-axis steps."""
    return float(sum(steps.values()))
