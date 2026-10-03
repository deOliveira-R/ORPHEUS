r"""The reference certificate: what a reference claims, per observable, and how it is established (#405 P2 step 5).

A reference's reading is a guarantee about the exact answer of its own
equation (the user's ruling "G3"). The certificate records, per observable,
the :class:`Claim`: the declared target the bound must meet, and how the
enclosure is ESTABLISHED. There are two ways, and only two (the user's
ruling of 2026-10-03, "no ladder certifies"):

* :class:`Exact`: the exact value as a symbolic expression with its
  self-check named; its enclosure is the double nearest it and the rounding
  between them, bounded from above;
* :class:`DerivedBound`: an a posteriori bound with computable constants (the
  Nyström residual bound :math:`\lVert x - x_n\rVert \le \lVert (I -
  K_n)^{-1}\rVert\,\lVert (K - K_n)x\rVert`, Atkinson 1997 Thm 4.1.2; a series
  tail; interval arithmetic; an integrator's own error control), carrying its
  enclosure and naming its method.

Richardson extrapolation is not a way: it estimates the limit of a sequence
and assumes the limit is the right one, which is what verification must
prove (ERR-006: the right rate to the wrong limit). A refinement sequence is
kept as FALSIFYING evidence only (:class:`Refinement`): correct enclosures at
every resolution all contain the exact answer, so they must share a point,
and the claim's enclosure must meet that common part. Its observed orders
are reported and never decide anything.

**Corroboration.** A :class:`PublishedSolution` of the same problem is an
anchor (:class:`Corroboration`): wherever it prints an observable the
certificate claims, the printed interval and the claimed enclosure must
intersect, decided on the outward-rounded ends. Monte Carlo and a fine-mesh
production run are never anchors (#505), and cannot be spelled as one.

**The state is derived, never stored** (:attr:`ReferenceCertificate.state`):
the :class:`~orpheus.reference.withdrawal.Withdrawal` when the standing is
one; otherwise :class:`Invalid` with every reason (a bound above its
target, an anchor disjoint from a claim, a refinement whose enclosures share
no point or miss the claim), or :class:`Valid`. A reference whose family
cannot yet derive a bound has no certificate, and so cannot anchor a
verification (the user's ruling of 2026-10-03: certification is decoupled
from caching).
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from fractions import Fraction
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any, TypeAlias, cast, final, get_args

from orpheus.numerics.content import ContentIdentity, FrozenMapping, content_digest
from orpheus.numerics.enclosure import Enclosure
from orpheus.numerics.mesh_free_function import parse_srepr
from orpheus.numerics.observable import Observable
from orpheus.numerics.scalars import parse_finite_real, parse_member, parse_positive_real
from orpheus.reference.published import Current, PublishedSolution, Standing
from orpheus.reference.withdrawal import Withdrawal

if TYPE_CHECKING:
    import sympy

__all__ = [
    "Claim",
    "Corroboration",
    "DerivedBound",
    "Establishment",
    "Exact",
    "Invalid",
    "ReferenceCertificate",
    "Refinement",
    "State",
    "Valid",
]

#: The working precision (decimal digits) an exact expression is evaluated at.
_EXACT_DIGITS = 60


def _named(text: object, where: str, noun: str) -> str:
    if not isinstance(text, str) or not text.strip():
        raise ValueError(f"{where}: {noun} is a non-empty text, got {text!r}")
    return text


@final
@dataclass(frozen=True, eq=False)
class Exact(ContentIdentity):
    """The exact value, a real constant expression (SymPy ``srepr`` text), and the self-check that proves it."""

    expression: str
    by: str

    def __post_init__(self) -> None:
        import sympy

        if not isinstance(self.expression, str):
            raise TypeError(f"Exact: the expression is srepr text, got a {type(self.expression).__name__}")
        parsed = sympy.sympify(parse_srepr(self.expression, "Exact: the expression"))
        if parsed.free_symbols or not parsed.is_real or not parsed.is_finite:
            raise ValueError(f"Exact: the expression must be a finite real constant, got {parsed}")
        object.__setattr__(self, "expression", sympy.srepr(parsed))
        object.__setattr__(self, "by", _named(self.by, "Exact", "the self-check"))
        content_digest(self)

    def enclosure(self) -> Enclosure:
        """The double nearest the exact value, and the rounding between them, bounded from above."""
        import mpmath
        import sympy

        value = cast("sympy.Expr", sympy.sympify(parse_srepr(self.expression, "Exact: the expression")))
        if isinstance(value, sympy.Rational):
            exact_rational = Fraction(int(value.p), int(value.q))
            centre = float(exact_rational)
            gap_rational = abs(Fraction(centre) - exact_rational)
            bound = float(gap_rational)
            if Fraction(bound) < gap_rational:
                bound = math.nextafter(bound, math.inf)
            return Enclosure(centre, bound)
        with mpmath.workdps(_EXACT_DIGITS):
            exact = mpmath.mpf(value.evalf(_EXACT_DIGITS))
            centre = float(exact)
            # the evaluation's own error is below one unit in the working precision's last digit
            evaluation = abs(exact) * mpmath.mpf(10) ** (1 - _EXACT_DIGITS)
            gap = abs(mpmath.mpf(centre) - exact) + evaluation
            bound = float(gap)
            if mpmath.mpf(bound) < gap:
                bound = math.nextafter(bound, math.inf)
        return Enclosure(centre, bound)


@final
@dataclass(frozen=True, eq=False)
class DerivedBound(ContentIdentity):
    """An enclosure established by an a posteriori bound with computable constants, and the method that derived it."""

    derived: Enclosure
    method: str

    def __post_init__(self) -> None:
        parse_member(self.derived, (Enclosure,), "DerivedBound", "the enclosure", "an enclosure")
        object.__setattr__(self, "method", _named(self.method, "DerivedBound", "the method"))
        content_digest(self)

    def enclosure(self) -> Enclosure:
        """The derived enclosure."""
        return self.derived


Establishment: TypeAlias = Exact | DerivedBound
"""How a claim's enclosure is established: exactly, or by a derived bound (never by extrapolation)."""


@final
@dataclass(frozen=True, eq=False)
class Claim(ContentIdentity):
    """One observable's claim: the target its bound must meet, and how its enclosure is established."""

    target: float
    established: Establishment

    def __post_init__(self) -> None:
        object.__setattr__(self, "target", parse_positive_real(self.target, "Claim", "the target"))
        parse_member(self.established, get_args(Establishment), "Claim", "the establishment", "an establishment")
        content_digest(self)

    def enclosure(self) -> Enclosure:
        """The claimed enclosure."""
        return self.established.enclosure()


@final
@dataclass(frozen=True, eq=False)
class Corroboration(ContentIdentity):
    """A published solution of the same problem, as an anchor, with why it is independent of the reference."""

    anchor: PublishedSolution
    independence: str

    def __post_init__(self) -> None:
        parse_member(self.anchor, (PublishedSolution,), "Corroboration", "the anchor", "a published solution")
        if isinstance(self.anchor.standing, Withdrawal):
            raise ValueError(
                f"Corroboration: the anchor is withdrawn (#{self.anchor.standing.issue}: {self.anchor.standing.reason})"
            )
        object.__setattr__(self, "independence", _named(self.independence, "Corroboration", "the independence note"))
        content_digest(self)


@final
@dataclass(frozen=True, eq=False)
class Refinement(ContentIdentity):
    """One observable's enclosures at a sequence of resolutions: falsifying evidence, never a certifier."""

    observable: Observable
    members: tuple[tuple[float, Enclosure], ...]

    def __post_init__(self) -> None:
        parse_member(self.observable, get_args(Observable), "Refinement", "the observable", "an observable")
        members = tuple(self.members)
        if len(members) < 2:
            raise ValueError("Refinement: a refinement sequence has at least two members")
        admitted = []
        for parameter, enclosure in members:
            admitted.append(
                (
                    parse_finite_real(parameter, "Refinement: a resolution parameter"),
                    parse_member(enclosure, (Enclosure,), "Refinement", "a member's reading", "an enclosure"),
                )
            )
        parameters = [parameter for parameter, _ in admitted]
        if any(later >= earlier for earlier, later in zip(parameters, parameters[1:])):
            raise ValueError(f"Refinement: the resolution parameters must decrease strictly (a refinement), got {parameters}")
        object.__setattr__(self, "members", tuple(admitted))
        content_digest(self)

    def common_part(self) -> tuple[float, float] | None:
        """The interval every member's enclosure contains, or ``None`` if they share no point."""
        low = max(enclosure.ends()[0] for _, enclosure in self.members)
        high = min(enclosure.ends()[1] for _, enclosure in self.members)
        return (low, high) if low <= high else None

    def observed_orders(self) -> tuple[float, ...]:
        """The order each three consecutive members suggest; reported, never decisive."""
        orders = []
        for (h0, a), (h1, b), (h2, c) in zip(self.members, self.members[1:], self.members[2:]):
            far, near = a.value - b.value, b.value - c.value
            if far == 0.0 or near == 0.0 or far / near <= 0.0 or h0 == h1 or h1 == h2:
                orders.append(math.nan)
                continue
            orders.append(math.log(far / near) / math.log(h0 / h1))
        return tuple(orders)


@final
@dataclass(frozen=True, eq=False)
class Valid(ContentIdentity):
    """Every claim meets its target, every anchor and every refinement agrees with it."""


@final
@dataclass(frozen=True, eq=False)
class Invalid(ContentIdentity):
    """At least one check failed; every failed check's reason."""

    reasons: tuple[str, ...]


State: TypeAlias = Valid | Invalid | Withdrawal
"""A reference certificate's state, derived from its claims, its evidence and its standing."""


def _meets(a: Enclosure, b: Enclosure) -> bool:
    low_a, high_a = a.ends()
    low_b, high_b = b.ends()
    return max(low_a, low_b) <= min(high_a, high_b)


@final
@dataclass(frozen=True, eq=False)
class ReferenceCertificate(ContentIdentity):
    """A reference's claims per observable, the evidence that can refute them, and its standing."""

    claims: Mapping[Observable, Claim]
    corroborations: tuple[Corroboration, ...] = ()
    refinements: tuple[Refinement, ...] = ()
    standing: Standing = field(default_factory=Current)

    def __post_init__(self) -> None:
        if not isinstance(self.claims, Mapping) or not self.claims:
            raise ValueError("ReferenceCertificate: a certificate claims at least one observable")
        for observable, claim in self.claims.items():
            parse_member(observable, get_args(Observable), "ReferenceCertificate", "a claimed key", "an observable")
            parse_member(claim, (Claim,), "ReferenceCertificate", f"the claim for {observable!r}", "a claim")
        object.__setattr__(self, "claims", FrozenMapping(self.claims.items()))
        corroborations = tuple(self.corroborations)
        for corroboration in corroborations:
            parse_member(corroboration, (Corroboration,), "ReferenceCertificate", "a corroboration", "a corroboration")
            if not any(observable in self.claims for observable in corroboration.anchor.printed):
                raise ValueError("ReferenceCertificate: an anchor shares no observable with the claims")
        object.__setattr__(self, "corroborations", corroborations)
        refinements = tuple(self.refinements)
        for refinement in refinements:
            parse_member(refinement, (Refinement,), "ReferenceCertificate", "a refinement", "a refinement")
            if refinement.observable not in self.claims:
                raise ValueError(f"ReferenceCertificate: a refinement of {refinement.observable!r}, which is not claimed")
        object.__setattr__(self, "refinements", refinements)
        parse_member(self.standing, get_args(Standing), "ReferenceCertificate", "the standing", "a standing")
        content_digest(self)

    @property
    def state(self) -> State:
        """Derived: the withdrawal, or every failed check, or valid."""
        match self.standing:
            case Withdrawal():
                return self.standing
            case Current():
                pass
        reasons: list[str] = []
        for observable, claim in self.claims.items():
            enclosure: Any = claim.enclosure()
            if enclosure.bound > claim.target:
                reasons.append(f"{observable!r}: the bound {enclosure.bound!r} exceeds the target {claim.target!r}")
        for corroboration in self.corroborations:
            for observable, printed in corroboration.anchor.printed.items():
                claim = self.claims.get(observable)
                if claim is not None and not _meets(claim.enclosure(), printed.enclosure()):
                    reasons.append(
                        f"{observable!r}: the claim {claim.enclosure()!r} and the anchor's printed "
                        f"{printed.text} ({printed.citation.bibkey}, {printed.citation.locator}) are disjoint"
                    )
        for refinement in self.refinements:
            common = refinement.common_part()
            if common is None:
                reasons.append(f"{refinement.observable!r}: the refinement's enclosures share no point")
                continue
            low, high = self.claims[refinement.observable].enclosure().ends()
            if max(low, common[0]) > min(high, common[1]):
                reasons.append(f"{refinement.observable!r}: the claim misses the refinement's common part {common!r}")
        return Invalid(tuple(reasons)) if reasons else Valid()
