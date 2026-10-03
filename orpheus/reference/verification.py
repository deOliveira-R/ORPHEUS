r"""Verification: a production answer's reading put on trial against a valid reference (#405 P2 step 7a).

Verification compares production's SELF-REPORT of one observable (a
:data:`~orpheus.numerics.outcome.ProductionReading`) with a reference's
GUARANTEE about the exact answer of the same observable (an
:class:`~orpheus.numerics.enclosure.Enclosure`), and only a reference whose
certificate is ``Valid`` may stand as the reference (the user's rulings "G3"
and "G4", 2026-10-02; a published solution is an anchor of a reference
certificate, never the reference itself; Monte Carlo never is one, #505).

**The verb reads.** :func:`verify_agreement` and :func:`verify_order` call
``answer.read(observable)`` and ``reference.read(observable)`` themselves, so
the reading verb is the only producer of a production reading, and a
diagnostic that happens to be a ``Measured`` cannot be passed in its place
(the elegance review of step 1, carried note S4).

**Two verdicts, two floors** (spec §0 item 3, ruled 2026-10-03):

* AGREEMENT: :math:`|m - v_{\rm ref}| + b_{\rm ref} \le \tau`, the
  comparison bounding the total error whatever its split, with the floor
  :math:`b_{\rm ref} \le \tau/10` (a reference too loose for the tolerance
  verifies nothing). Production's algebraic-error evidence is recorded, not
  required.
* ORDER: a refinement of production answers whose observed orders are
  compared with the declared order. It attributes error to the
  discretisation, so the algebraic error of every answer must be established
  (``Measured`` or ``Asserted`` at most :math:`\tau/10`; ``NotYet`` and
  ``NotApplicable`` are refused as unestablished, #564), the reference's
  bound must sit below :math:`\tau/10`, and every error must be resolved
  above both floors (at least :math:`\tau`). The order and its band are the
  caller's declared inputs: this verdict is a claim about PRODUCTION's order
  of accuracy against a reference, not a reference certifying itself, so
  the user's ruling that no ladder certifies a reference does not reach it.

Both results are returned and never stored: they are a test's evidence, not
data. Declared limit: nothing yet checks that the answer answers the
reference's specification (production results hold no specification until
#405 P4's projection of a specification onto a mesh).
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass
from fractions import Fraction
from typing import Protocol, final, get_args, runtime_checkable

from orpheus.numerics.enclosure import Enclosure
from orpheus.numerics.observable import Observable
from orpheus.numerics.outcome import Asserted, Evidence, Measured, NotApplicable, NotYet, ProductionReading
from orpheus.numerics.scalars import parse_finite_real, parse_member, parse_positive_real
from orpheus.reference.certificate import Invalid, Valid
from orpheus.reference.solution import ReferenceSolution
from orpheus.reference.withdrawal import Withdrawal

__all__ = [
    "OrderVerification",
    "ProductionAnswer",
    "ReferenceNotValid",
    "Unestablished",
    "VerificationCertificate",
    "verify_agreement",
    "verify_order",
]


class ReferenceNotValid(ValueError):
    """The reference has no certificate, or its certificate is not ``Valid``."""


class Unestablished(ValueError):
    """An order verdict needs a floor (an algebraic error, a reference bound, a resolved error) that is not established."""


@runtime_checkable
class ProductionAnswer(Protocol):
    """A production answer that reads an observable as its self-report."""

    def read(self, observable: Observable) -> ProductionReading:
        """The answer's reading of ``observable``."""
        ...


def _require_valid(reference: ReferenceSolution) -> None:
    """Refuse a reference that may not stand as one, naming why, before anything is read."""
    parse_member(reference, (ReferenceSolution,), "verification", "the reference", "a reference solution")
    if reference.certificate is None:
        raise ReferenceNotValid("verification: the reference has no certificate, so it cannot stand as a reference")
    match reference.certificate.state:
        case Valid():
            return
        case Invalid(reasons=reasons):
            raise ReferenceNotValid(f"verification: the reference's certificate is invalid ({reasons[0]})")
        case Withdrawal(reason=reason, issue=issue):
            raise ReferenceNotValid(f"verification: the reference is withdrawn (#{issue}: {reason})")


def _production_reading(answer: ProductionAnswer, observable: Observable) -> Measured:
    if not isinstance(answer, ProductionAnswer):
        raise TypeError(f"verification: the answer reads an observable (a ProductionAnswer), got a {type(answer).__name__}")
    reading = answer.read(observable)
    parse_member(reading, (Measured,), "verification", "the answer's reading", "a production reading")
    return reading


def _algebraic(evidence: Evidence) -> Evidence:
    parse_member(evidence, get_args(Evidence), "verification", "the algebraic error", "an Evidence")
    return evidence


def _error(reading: Measured, reference: Enclosure) -> Fraction:
    """The production reading's distance to the reference's centre, exactly."""
    return abs(Fraction(reading.value) - Fraction(reference.value))


@final
@dataclass(frozen=True)
class VerificationCertificate:
    """One comparison of a production reading with a valid reference's enclosure, and its verdict."""

    observable: Observable
    reading: Measured
    reference_reading: Enclosure
    tolerance: float
    algebraic_error: Evidence

    @property
    def floor_holds(self) -> bool:
        """The reference's bound is at most a tenth of the tolerance."""
        return Fraction(self.reference_reading.bound) <= Fraction(self.tolerance) / 10

    @property
    def agrees(self) -> bool:
        """The floor holds and the total error, reference bound included, is within the tolerance."""
        error = _error(self.reading, self.reference_reading) + Fraction(self.reference_reading.bound)
        return self.floor_holds and error <= Fraction(self.tolerance)

    def require(self) -> None:
        """Raise unless the comparison verifies: the floor first, then the agreement."""
        if not self.floor_holds:
            raise AssertionError(
                f"{self.observable!r}: the reference's bound {self.reference_reading.bound!r} is above the floor "
                f"tolerance/10 = {self.tolerance / 10!r}; the reference is too loose to verify this tolerance"
            )
        if not self.agrees:
            raise AssertionError(
                f"{self.observable!r}: production reads {self.reading.value!r}, which disagrees with the reference "
                f"{self.reference_reading!r} at tolerance {self.tolerance!r}"
            )


def verify_agreement(
    answer: ProductionAnswer,
    observable: Observable,
    reference: ReferenceSolution,
    tolerance: float,
    algebraic_error: Evidence,
) -> VerificationCertificate:
    """Compare production's reading of ``observable`` with a valid reference's enclosure."""
    _require_valid(reference)
    algebraic = _algebraic(algebraic_error)
    bound = parse_positive_real(tolerance, "verification", "the tolerance")
    reference_reading = reference.read(observable)
    reading = _production_reading(answer, observable)
    return VerificationCertificate(observable, reading, reference_reading, bound, algebraic)


@final
@dataclass(frozen=True)
class OrderVerification:
    """Production's observed order of accuracy over a refinement, against a valid reference."""

    observable: Observable
    resolutions: tuple[float, ...]
    errors: tuple[Fraction, ...]
    order: float
    band: float

    @property
    def observed_orders(self) -> tuple[float, ...]:
        """One observed order per consecutive pair of resolutions."""
        return tuple(
            math.log(float(e0 / e1)) / math.log(h0 / h1)
            for (h0, e0), (h1, e1) in zip(zip(self.resolutions, self.errors), zip(self.resolutions[1:], self.errors[1:]))
        )

    @property
    def holds(self) -> bool:
        """Every observed order lies within the declared band about the declared order."""
        return all(abs(observed - self.order) <= self.band for observed in self.observed_orders)


def verify_order(
    answers: Sequence[tuple[float, ProductionAnswer]],
    observable: Observable,
    reference: ReferenceSolution,
    tolerance: float,
    algebraic_errors: Sequence[Evidence],
    order: float,
    band: float,
) -> OrderVerification:
    """Production's observed order over a refinement, every floor established."""
    _require_valid(reference)
    bound = parse_positive_real(tolerance, "verification", "the tolerance")
    declared = parse_finite_real(order, "verification: the declared order")
    width = parse_positive_real(band, "verification", "the band")
    pairs = tuple(answers)
    evidence = tuple(_algebraic(e) for e in algebraic_errors)
    if len(pairs) < 2 or len(evidence) != len(pairs):
        raise ValueError("verification: an order needs at least two answers, each with its algebraic error")
    floor = Fraction(bound) / 10
    for index, item in enumerate(evidence):
        match item:
            case Measured(value=value) if Fraction(abs(value)) <= floor:
                continue
            case Asserted(bound=asserted) if Fraction(asserted) <= floor:
                continue
            case Measured() | Asserted():
                raise Unestablished(f"verification: answer {index}'s algebraic error {item!r} is above tolerance/10")
            case NotYet() | NotApplicable():
                raise Unestablished(f"verification: answer {index}'s algebraic error is unestablished ({item!r})")
    reference_reading = reference.read(observable)
    if Fraction(reference_reading.bound) > floor:
        raise Unestablished(f"verification: the reference's bound {reference_reading.bound!r} is above tolerance/10")
    resolutions = tuple(parse_finite_real(h, "verification: a resolution") for h, _ in pairs)
    errors = tuple(_error(_production_reading(answer, observable), reference_reading) for _, answer in pairs)
    for index, error in enumerate(errors):
        if error < Fraction(bound):
            raise Unestablished(
                f"verification: answer {index}'s error {float(error)!r} is unresolved: below the tolerance, it is "
                f"not separated from the floors"
            )
    return OrderVerification(observable, resolutions, errors, declared, width)
