r"""A reference's reading of an observable: an enclosure, or a printed value (#405 P2 step 1).

A reading is a claim about one observable, and whose claim it is decides its
type (the user's ruling of 2026-10-02, "G3"):

* a **reference's** claim is about the exact answer of its own equation.
  :data:`ReferenceReading` is an :class:`~orpheus.numerics.enclosure.Enclosure`
  (a guarantee the reference derives: the value and a bound on its distance
  to the exact one; an exact value has bound 0) or a :class:`Printed` value
  (a cited author's claim, never recomputed);
* a **production** method's claim is a self-report that verification puts on
  trial: :data:`~orpheus.numerics.outcome.ProductionReading`.

The two sums share no member, so a production reading cannot be passed where
a reference's claim is required. A statistical interval (Monte Carlo's
``Estimated``) is a production claim only, and Monte Carlo is never a
reference (#505).

**A printed value is its decimal text.** The trailing zeros a float drops are
the precision claim (``1.0`` and ``1.00`` are two claims), so the text is the
one source: it is canonicalised by :class:`~decimal.Decimal` (``1.0E0``,
``+1.0`` and ``10E-1`` are one claim), and the value (the correctly rounded
double) and the half unit in the last printed digit are derived from it,
never stored beside it. The citation must name the place the value is
printed (a table, a page): a value is printed at a place, and 8 of the 47
Sood cases print their values in different places (P1 step 2c).
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from decimal import Decimal, InvalidOperation
from typing import TypeAlias

from orpheus.data.citation import Citation
from orpheus.numerics.content import ContentIdentity, content_digest
from orpheus.numerics.enclosure import Enclosure

__all__ = ["Printed", "ReferenceReading"]


@dataclass(frozen=True, eq=False)
class Printed(ContentIdentity):
    """A value as printed in a cited work: its decimal text and the place it is printed."""

    text: str
    citation: Citation

    def __post_init__(self) -> None:
        if not isinstance(self.text, str):
            raise TypeError(f"Printed: a printed value is its decimal text, got a {type(self.text).__name__}")
        try:
            number = Decimal(self.text.strip())
        except InvalidOperation:
            raise ValueError(f"Printed: {self.text!r} is not a printed decimal number") from None
        if not number.is_finite():
            raise ValueError(f"Printed: {self.text!r} is not a finite printed number")
        if not math.isfinite(float(number)):
            raise ValueError(f"Printed: {self.text!r} lies beyond the range of a double")
        object.__setattr__(self, "text", str(number.copy_abs() if number.is_zero() else number))
        if not isinstance(self.citation, Citation):
            raise TypeError(f"Printed: the citation is a Citation, got a {type(self.citation).__name__}")
        if self.citation.locator is None:
            raise ValueError(
                f"Printed: the citation {self.citation.bibkey!r} has no locator; a printed value is printed at a place"
            )
        content_digest(self)

    @property
    def value(self) -> float:
        """The printed number, correctly rounded to a double."""
        return float(Decimal(self.text))

    @property
    def half_unit(self) -> Decimal:
        """Half a unit in the last printed digit, exactly."""
        last_digit = int(Decimal(self.text).as_tuple().exponent)  # finite: admitted at construction
        return Decimal(5).scaleb(last_digit - 1)


ReferenceReading: TypeAlias = Enclosure | Printed
"""A reference's claim about one observable of the exact answer of its equation."""
