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

A certificate holds the answer and the reference and reads both itself at
construction, so a certificate against a reference that is not ``Valid``, or
over a reading the answer did not produce, cannot be built (the elegance
review of step 7a). Both results are returned and never stored: they are a test's evidence, not
data. Declared limit: nothing yet checks that the answer answers the
reference's specification (production results hold no specification until
#405 P4's projection of a specification onto a mesh).
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass, field
from fractions import Fraction
from typing import Protocol, final, get_args, runtime_checkable

from orpheus.numerics.enclosure import Enclosure
from orpheus.numerics.observable import Observable
from orpheus.numerics.outcome import PRODUCTION_READINGS, Asserted, Evidence, Measured, NotApplicable, NotYet, ProductionReading
from orpheus.numerics.scalars import parse_member, parse_positive_real
from orpheus.reference.certificate import Invalid, Valid, parse_resolutions
from orpheus.reference.solution import ReferenceSolution
from orpheus.reference.withdrawal import Withdrawal

__all__ = [
    "Disagreement",
    "OrderNotObserved",
    "OrderVerification",
    "ProductionAnswer",
    "ReferenceNotValid",
    "ReferenceTooLoose",
    "Unestablished",
    "VerificationCertificate",
    "observed_order",
    "verify_agreement",
    "verify_order",
]


class ReferenceNotValid(ValueError):
    """The reference has no certificate, or its certificate is not ``Valid``."""


class Unestablished(ValueError):
    """An order verdict needs a floor (an algebraic error, a reference bound, a resolved error) that is not established."""


class ReferenceTooLoose(AssertionError):
    """The reference's bound is above a tenth of the tolerance: it is too loose to verify that tolerance."""


class Disagreement(AssertionError):
    """Production's reading and the reference's enclosure disagree at the tolerance."""


class OrderNotObserved(AssertionError):
    """An observed order lies outside the declared band about the declared order."""


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


def _production_reading(answer: ProductionAnswer, observable: Observable) -> ProductionReading:
    if not isinstance(answer, ProductionAnswer):
        raise TypeError(f"verification: the answer reads an observable (a ProductionAnswer), got a {type(answer).__name__}")
    reading = answer.read(observable)
    parse_member(reading, PRODUCTION_READINGS, "verification", "the answer's reading", "a production reading")
    return reading


def _algebraic(evidence: Evidence) -> Evidence:
    parse_member(evidence, get_args(Evidence), "verification", "the algebraic error", "an Evidence")
    return evidence


def _error(reading: ProductionReading, reference: Enclosure) -> Fraction:
    """The production reading's distance to the reference's centre, exactly."""
    return abs(Fraction(reading.value) - Fraction(reference.value))


def _established_bound(evidence: Evidence) -> Fraction:
    """The bound an established algebraic error contributes (its measured size, or its asserted bound)."""
    match evidence:
        case Measured(value=value):
            return Fraction(abs(value))
        case Asserted(bound=bound):
            return Fraction(bound)
        case _:
            raise Unestablished(f"verification: the algebraic error {evidence!r} is unestablished")


def observed_order(coarse: tuple[float, Fraction], fine: tuple[float, Fraction]) -> float:
    """The order two (resolution, error) pairs suggest: log(e_coarse / e_fine) / log(h_coarse / h_fine)."""
    (h0, e0), (h1, e1) = coarse, fine
    return math.log(float(e0 / e1)) / math.log(h0 / h1)


@final
@dataclass(frozen=True)
class VerificationCertificate:
    """One comparison of a production answer's reading with a valid reference's enclosure, and its verdict.

    Built from the answer and the reference, it reads both at construction:
    the reference must be ``Valid`` before anything is read.
    """

    answer: ProductionAnswer
    observable: Observable
    reference: ReferenceSolution
    tolerance: float
    algebraic_error: Evidence
    reading: ProductionReading = field(init=False)
    reference_reading: Enclosure = field(init=False)

    def __post_init__(self) -> None:
        _require_valid(self.reference)
        _algebraic(self.algebraic_error)
        object.__setattr__(self, "tolerance", parse_positive_real(self.tolerance, "verification: the tolerance", "the tolerance"))
        object.__setattr__(self, "reference_reading", self.reference.read(self.observable))
        object.__setattr__(self, "reading", _production_reading(self.answer, self.observable))

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
        """Raise unless the comparison verifies: the floor first (:class:`ReferenceTooLoose`), then :class:`Disagreement`."""
        if not self.floor_holds:
            raise ReferenceTooLoose(
                f"{self.observable!r}: the reference's bound {self.reference_reading.bound!r} is above the floor "
                f"tolerance/10 = {self.tolerance / 10!r}; the reference is too loose to verify this tolerance"
            )
        if not self.agrees:
            raise Disagreement(
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
    return VerificationCertificate(answer, observable, reference, tolerance, algebraic_error)


@final
@dataclass(frozen=True)
class OrderVerification:
    r"""Production's observed order of accuracy over a refinement, against a valid reference.

    Each answer's error :math:`e_i = m_i - v_{\rm ref}` is known only within
    :math:`u_i`, the reference's bound plus that answer's established
    algebraic error, so each observed order is an INTERVAL: the range of
    :math:`\log(|e_i| / |e_{i+1}|) / \log(h_i / h_{i+1})` over the error
    intervals. The verdict holds only when every such interval lies inside the
    declared band about the declared order (qa of step 7a: point estimates let
    a true order 1.7 read 1.963, inside 2 ± 0.1, with both floors at their
    limit), and the errors keep one sign (a sign change, 1.04 then 0.99, is
    not monotone convergence, whatever order it suggests).
    """

    observable: Observable
    resolutions: tuple[float, ...]
    errors: tuple[Fraction, ...]
    uncertainties: tuple[Fraction, ...]
    order: float
    band: float

    def __post_init__(self) -> None:
        object.__setattr__(self, "resolutions", parse_resolutions(self.resolutions, "OrderVerification"))
        if not len(self.errors) == len(self.uncertainties) == len(self.resolutions):
            raise ValueError(
                f"OrderVerification: {len(self.errors)} errors and {len(self.uncertainties)} uncertainties for "
                f"{len(self.resolutions)} resolutions; one each"
            )
        if any(abs(error) <= uncertainty for error, uncertainty in zip(self.errors, self.uncertainties)):
            raise ValueError("OrderVerification: an error is not resolved above its uncertainty")
        object.__setattr__(self, "order", parse_positive_real(self.order, "OrderVerification: the order", "the order"))
        object.__setattr__(self, "band", parse_positive_real(self.band, "OrderVerification: the band", "the band"))

    @property
    def monotone(self) -> bool:
        """The errors keep one sign across the refinement."""
        return all(error > 0 for error in self.errors) or all(error < 0 for error in self.errors)

    @property
    def observed_orders(self) -> tuple[tuple[float, float], ...]:
        """Per consecutive pair, the interval of orders the error intervals admit (low, high)."""
        intervals = []
        for i in range(len(self.resolutions) - 1):
            h0, h1 = self.resolutions[i], self.resolutions[i + 1]
            e0, e1 = abs(self.errors[i]), abs(self.errors[i + 1])
            u0, u1 = self.uncertainties[i], self.uncertainties[i + 1]
            low = observed_order((h0, e0 - u0), (h1, e1 + u1))
            high = observed_order((h0, e0 + u0), (h1, e1 - u1))
            intervals.append((low, high))
        return tuple(intervals)

    @property
    def holds(self) -> bool:
        """The errors keep one sign and every observed-order interval lies inside the declared band."""
        return self.monotone and all(
            self.order - self.band <= low and high <= self.order + self.band for low, high in self.observed_orders
        )

    def require(self) -> None:
        """Raise :class:`OrderNotObserved` unless the verdict holds."""
        if not self.holds:
            raise OrderNotObserved(
                f"{self.observable!r}: the observed-order intervals {self.observed_orders} (errors of one sign: "
                f"{self.monotone}) are not within {self.band!r} of the declared order {self.order!r}"
            )


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
    bound = parse_positive_real(tolerance, "verification: the tolerance", "the tolerance")
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
    resolutions = parse_resolutions((h for h, _ in pairs), "verification")
    errors = tuple(
        Fraction(_production_reading(answer, observable).value) - Fraction(reference_reading.value) for _, answer in pairs
    )
    uncertainties = tuple(Fraction(reference_reading.bound) + _established_bound(item) for item in evidence)
    for index, error in enumerate(errors):
        if abs(error) < Fraction(bound):
            raise Unestablished(
                f"verification: answer {index}'s error {float(error)!r} is unresolved: below the tolerance, it is "
                f"not separated from the floors"
            )
    return OrderVerification(observable, resolutions, errors, uncertainties, order, band)
