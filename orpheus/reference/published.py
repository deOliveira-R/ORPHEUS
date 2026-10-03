r"""A published solution: values printed in a cited work, never recomputed (#405 P2 step 5).

A :class:`PublishedSolution` holds a specification and, per observable, the
value a work prints for it, each with its own citation
(:class:`~orpheus.reference.reading.Printed`): 8 of the 47 Sood cases print
their values in different places, so citations are per value, never pooled.
It reads what it prints and refuses anything else; it is never a source of
a computed value.

**Standing.** A publication is :class:`Current`, or withdrawn by a maintainer
ruling recorded in an issue, :class:`~orpheus.reference.withdrawal.Withdrawal`
``(reason, issue)`` (a publication found unfit as a whole). An ERRATUM is not a
withdrawal: it is a cited correction to one value, so it is data, a
:class:`~orpheus.reference.reading.Printed` citing the erratum in place of
the value it corrects. The
same sum, :data:`Standing`, is the standing of a reference certificate, so a
withdrawn reference reads as withdrawn wherever it is held. A withdrawn
publication refuses every read, naming its issue.

**Admission.** The specification is a
:data:`~orpheus.specification.Specification`; every printed observable must
be one the specification can pose
(:func:`~orpheus.specification.admit_observable`, the one admission every
reader reuses), so a publication cannot print an eigenvalue for a source
question or a point value on the infinite medium.

A published solution is a corroborating ANCHOR for a reference certificate
(:class:`~orpheus.reference.certificate.Corroboration`), and never a
reference a production solution is verified against (the user's ruling of
2026-09-25).
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any, TypeAlias, final, get_args

from orpheus.numerics.content import ContentIdentity, FrozenMapping, content_digest
from orpheus.numerics.observable import Observable
from orpheus.numerics.scalars import parse_member
from orpheus.reference.reading import Printed
from orpheus.reference.withdrawal import Withdrawal
from orpheus.specification import Specification, admit_observable

__all__ = ["Current", "NotPrinted", "PublishedSolution", "Standing"]


@final
@dataclass(frozen=True, eq=False)
class Current(ContentIdentity):
    """In good standing: not withdrawn."""


Standing: TypeAlias = Current | Withdrawal
"""Whether a publication or a reference certificate stands, or is withdrawn with its reason and issue."""


class NotPrinted(LookupError):
    """The publication does not print the observable asked for."""


@final
@dataclass(frozen=True, eq=False)
class PublishedSolution(ContentIdentity):
    """The values a cited work prints for a specification, per observable."""

    specification: Specification
    printed: Mapping[Observable, Printed]
    standing: Standing

    def __post_init__(self) -> None:
        parse_member(self.specification, get_args(Specification), "PublishedSolution", "the specification", "a specification")
        parse_member(self.standing, get_args(Standing), "PublishedSolution", "the standing", "a standing")
        if not isinstance(self.printed, Mapping) or not self.printed:
            raise ValueError("PublishedSolution: the publication prints nothing; it maps at least one observable to its printed value")
        for observable, value in self.printed.items():
            parse_member(observable, get_args(Observable), "PublishedSolution", "a printed key", "an observable")
            parse_member(value, (Printed,), "PublishedSolution", f"the value printed for {observable!r}", "a printed value")
            admit_observable(observable, self.specification)
        object.__setattr__(self, "printed", FrozenMapping(self.printed.items()))
        content_digest(self)

    def read(self, observable: Observable) -> Printed:
        """The value printed for ``observable``; refused if withdrawn or not printed."""
        match self.standing:
            case Withdrawal(reason=reason, issue=issue):
                raise ValueError(f"PublishedSolution: withdrawn (#{issue}: {reason}); it cannot be read")
            case Current():
                pass
        printed: Any = self.printed.get(observable)
        if printed is None:
            raise NotPrinted(
                f"PublishedSolution: the publication does not print {observable!r}; it prints "
                f"{', '.join(repr(key) for key in self.printed)}"
            )
        return printed
