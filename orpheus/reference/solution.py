r"""The reference solution: a reference method's answer, read on demand (#405 P2 step 6).

A :class:`ReferenceSolution` is the output of a reference method, never of a
production solver (the user's ruling "G5"): the specification it answers,
the :class:`Derivation` that evaluates its answer, and the
:class:`~orpheus.reference.certificate.ReferenceCertificate` that declares
targets and corroborates them, when the family has one.

**Read on demand.** ``read(observable)`` admits the observable once
(:func:`~orpheus.specification.specification.admit_observable`), then asks the
derivation to EVALUATE it: the reference's natural extension (for an
integral-equation reference, the transport integral of its emission density),
as an :data:`Evaluation`. Where the family derives a bound the evaluation is
an :data:`~orpheus.reference.certificate.Establishment` and reads as its
enclosure; where it cannot, it is an
:class:`~orpheus.reference.reading.Uncertified` value and reads as itself
(the user's ruling of 2026-10-03, step 7b: no float without a bound passes
for a guarantee, and none is withheld either). Nothing about the answer
is stored here: an observable a test poses later, the cell averages of its own
mesh, is read the same way as one the certificate declared (the user's rulings
of 2026-10-03: the reading is the natural extension; no stored field in P2;
P3's traced cache will make repeated readings cheap). A ratio reads as the
quotient of its operands' readings, its bound carried outward (uncertified
when either operand is), so a ratio is never evaluated, and never claimed, on
its own: one definition.

**The certificate's role.** A derived enclosure is a guarantee on its own
(G3). The certificate adds, for the observables it claims, a target and the
evidence that can refute the claim. A claim and the derivation's
establishment of the same observable are one quantity in two places (X4), so
they are checked once, at construction: a reference whose derivation
disagrees with its own claim cannot be built ("disagrees"), nor one that
claims an observable its derivation leaves uncertified, and ``read`` is a
pure evaluation.

**The same problem.** Construction refuses a claim the specification cannot
pose, and an anchor whose specification is not this one (qa of step 5, F2).

Not :class:`~orpheus.numerics.content.ContentIdentity` yet: P3 keys a
reference by its derivation's execution trace, so until then two reference
solutions are equal only if they are the same object.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol, TypeAlias, assert_never, final, get_args, runtime_checkable

from orpheus.numerics.enclosure import Enclosure, common_part
from orpheus.numerics.observable import Eigenvalue, Linear, Observable, Ratio
from orpheus.numerics.scalars import parse_member
from orpheus.reference.certificate import DerivedBound, Establishment, Exact, ReferenceCertificate
from orpheus.reference.reading import Uncertified
from orpheus.specification.specification import Specification, admit_observable

__all__ = ["Derivation", "Evaluation", "ReferenceSolution"]

Evaluation: TypeAlias = Establishment | Uncertified
"""A derivation's value of one observable: established (exactly or by a derived bound), or uncertified."""


@runtime_checkable
class Derivation(Protocol):
    """A reference method's natural extension, evaluated with its derived bound where the family has one."""

    def evaluate(self, observable: Eigenvalue | Linear) -> Evaluation:
        """The observable's value: established with its bound, or :class:`~orpheus.reference.reading.Uncertified`."""
        ...


@final
@dataclass(frozen=True, eq=False)
class ReferenceSolution:
    """A reference method's answer to a specification, read on demand through its derivation."""

    specification: Specification
    derivation: Derivation
    certificate: ReferenceCertificate | None

    def __post_init__(self) -> None:
        parse_member(self.specification, get_args(Specification), "ReferenceSolution", "the specification", "a specification")
        if not isinstance(self.derivation, Derivation):
            raise TypeError(
                f"ReferenceSolution: the derivation evaluates an observable (a Derivation), "
                f"got a {type(self.derivation).__name__}"
            )
        if self.certificate is None:
            return
        parse_member(self.certificate, (ReferenceCertificate,), "ReferenceSolution", "the certificate", "a reference certificate")
        for observable in self.certificate.claims:
            if isinstance(observable, Ratio):
                raise ValueError(
                    f"ReferenceSolution: {observable!r} is claimed, and a ratio reads as the quotient of its "
                    f"operands' readings, never on its own"
                )
            admit_observable(observable, self.specification)
            established = self._read(observable)
            if isinstance(established, Uncertified):
                raise ValueError(
                    f"ReferenceSolution: {observable!r} is claimed, and the derivation leaves it uncertified "
                    f"(no derived bound), so it cannot be claimed"
                )
            claimed = self.certificate.claims[observable].enclosure()
            if common_part([claimed, established]) is None:
                raise ValueError(
                    f"ReferenceSolution: the derivation establishes {established!r}, which disagrees with the "
                    f"claimed {claimed!r} for {observable!r}"
                )
        for corroboration in self.certificate.corroborations:
            if corroboration.anchor.specification != self.specification:
                raise ValueError(
                    "ReferenceSolution: an anchor answers another specification, so it cannot corroborate this one"
                )

    def read(self, observable: Observable) -> Enclosure | Uncertified:
        """The reference's reading of ``observable``, evaluated on demand (admitted once): enclosed, or uncertified."""
        admit_observable(observable, self.specification)
        return self._read(observable)

    def _read(self, observable: Observable) -> Enclosure | Uncertified:
        if isinstance(observable, Ratio):
            return observable.quotient(self._read)
        match self._evaluate(observable):
            case Exact() | DerivedBound() as established:
                return established.enclosure()
            case Uncertified() as uncertified:
                return uncertified
            case unreachable:
                assert_never(unreachable)

    def _evaluate(self, observable: Eigenvalue | Linear) -> Evaluation:
        """The derivation's evaluation, admitted as an :data:`Evaluation` (a Derivation is an open protocol)."""
        evaluation = self.derivation.evaluate(observable)
        parse_member(evaluation, get_args(Evaluation), "ReferenceSolution", "the derivation's evaluation", "an evaluation")
        return evaluation
