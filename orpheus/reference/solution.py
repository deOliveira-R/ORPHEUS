r"""The reference solution: a reference method's answer, read on demand (#405 P2 step 6).

A :class:`ReferenceSolution` is the output of a reference method, never of a
production solver (the user's ruling "G5"): the specification it answers,
the :class:`Derivation` that evaluates its answer, and the
:class:`~orpheus.reference.certificate.ReferenceCertificate` that declares
targets and corroborates them, when the family has one.

**Read on demand.** ``read(observable)`` admits the observable once
(:func:`~orpheus.specification.specification.admit_observable`), then asks the
derivation to ESTABLISH it: to evaluate the reference's natural extension (for
an integral-equation reference, the transport integral of its emission
density) together with a derived bound, as an
:data:`~orpheus.reference.certificate.Establishment`. Nothing about the answer
is stored here: an observable a test poses later, the cell averages of its own
mesh, is read the same way as one the certificate declared (the user's rulings
of 2026-10-03: the reading is the natural extension; no stored field in P2;
P3's traced cache will make repeated readings cheap). A ratio reads as the
quotient of its operands' readings, its bound carried outward, so a ratio is
never established, and never claimed, on its own: one definition.

**The certificate's role.** A derived enclosure is a guarantee on its own
(G3). The certificate adds, for the observables it claims, a target and the
evidence that can refute the claim. A claim and the derivation's
establishment of the same observable are one quantity in two places (X4), so
they are checked once, at construction: a reference whose derivation
disagrees with its own claim cannot be built ("disagrees"), and ``read`` is a
pure evaluation. A family whose derivation cannot
derive a bound raises :class:`NotCertified`; what such a reference returns
instead waits for the user's ruling on the uncertified reading.

**The same problem.** Construction refuses a claim the specification cannot
pose, and an anchor whose specification is not this one (qa of step 5, F2).

Not :class:`~orpheus.numerics.content.ContentIdentity` yet: P3 keys a
reference by its derivation's execution trace, so until then two reference
solutions are equal only if they are the same object.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol, final, get_args, runtime_checkable

from orpheus.numerics.enclosure import Enclosure, common_part
from orpheus.numerics.observable import Eigenvalue, Linear, Observable, Ratio
from orpheus.numerics.scalars import parse_member
from orpheus.reference.certificate import Establishment, ReferenceCertificate
from orpheus.specification.specification import Specification, admit_observable

__all__ = ["Derivation", "NotCertified", "ReferenceSolution"]


class NotCertified(LookupError):
    """The reference cannot establish a derived bound for the observable asked for."""


@runtime_checkable
class Derivation(Protocol):
    """A reference method's natural extension, evaluated with a derived bound."""

    def establish(self, observable: Eigenvalue | Linear) -> Establishment:
        """The observable's value and its derived bound; raises :class:`NotCertified` if none can be derived."""
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
                f"ReferenceSolution: the derivation establishes an observable (a Derivation), "
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
            established = self.derivation.establish(observable).enclosure()
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

    def read(self, observable: Observable) -> Enclosure:
        """The reference's enclosure of ``observable``, established on demand (admitted once)."""
        admit_observable(observable, self.specification)
        return self._establish(observable)

    def _establish(self, observable: Observable) -> Enclosure:
        if isinstance(observable, Ratio):
            return self._establish(observable.numerator) / self._establish(observable.denominator)
        return self.derivation.establish(observable).enclosure()
