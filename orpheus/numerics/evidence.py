r"""The EVIDENCE sum — what a solve can say about a property of its answer, typed.

A property of a returned state (its balance defect, its gauge displacement, its
solvability's conditioning) is either a number that was measured, a bound that
was asserted, a question that does not arise, or open work: never a ``None``
with several documented meanings. The four members are the exit report's
vocabulary (:class:`~orpheus.numerics.outcome.ExitReport`) and the Fredholm
verdict's (:mod:`orpheus.numerics.fredholm`), so they rank below both: they
import nothing.
"""

from __future__ import annotations

from dataclasses import dataclass

__all__ = ["Asserted", "Evidence", "Measured", "NotApplicable", "NotYet"]


@dataclass(frozen=True)
class Measured:
    """A number was measured on the returned state."""

    value: float

    def __truediv__(self, other: Measured) -> Measured:
        """The quotient of two measurements is a measurement (a ratio observable's production reading)."""
        if not isinstance(other, Measured):
            return NotImplemented
        return Measured(self.value / other.value)


@dataclass(frozen=True)
class Asserted:
    """No number was needed: the convergence-claim check ASSERTED the bound (raising on a defect beyond it)."""

    bound: float
    by: str


@dataclass(frozen=True)
class NotApplicable:
    """The question does not arise for this posing / configuration — say why."""

    reason: str


@dataclass(frozen=True)
class NotYet:
    """Open work — the issue that will make it measurable, and what is missing."""

    issue: int
    reason: str


Evidence = Measured | Asserted | NotApplicable | NotYet
