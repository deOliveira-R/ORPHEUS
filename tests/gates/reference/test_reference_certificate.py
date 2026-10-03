r"""The reference certificate (#405 P2 step 5, R5c).

Specified by the test-architect (2026-10-03, ``.claude/plans/reference_p2_spec.md``
§1.5), on the user's step-5 ruling ("no ladder certifies"):

* a claim per observable pairs a target :math:`\tau` with how its bound is
  ESTABLISHED: ``Exact`` (an exact expression, stored as SymPy ``srepr``
  text, whose enclosure is its correctly rounded double and the rounding) or
  ``DerivedBound`` (an a-posteriori bound with computable constants, an
  ``Enclosure`` and its method). Each yields an ``Enclosure``;
* corroborations: an anchor (a ``PublishedSolution`` in step 5) with an
  independence note; its printed interval must meet the claim's enclosure;
* refinements: FALSIFYING evidence only. A refinement sequence's enclosures
  must share a common point (every correct bound contains the exact answer,
  so a disjoint pair proves a bound wrong); its observed orders are reported,
  never a certifier;
* the state ``Valid | Invalid(reasons) | Withdrawal`` is DERIVED from the
  evidence and the standing, never stored.

Every interval decision is taken in exact arithmetic (``Fraction``,
``Decimal``) so a boundary row is decidable.
"""

from __future__ import annotations

import ast
import dataclasses
import math
from fractions import Fraction
from pathlib import Path
from typing import Any

import pytest

from orpheus.numerics.content import content_digest
from tests.gates._content_identity_helpers import (
    Entry,
    check_equal_pair,
    check_perturbation,
    check_pickle,
    check_population,
    leg,
    pair_ids,
    param_id,
    perturbation_ids,
    require,
)
from tests.gates.reference import _step5 as s5

pytestmark = pytest.mark.foundation

_HERE = "tests/gates/reference/test_reference_certificate.py"
_ROOT = Path(__file__).resolve().parents[3]


def _interval(e: Any) -> tuple[Fraction, Fraction]:
    return Fraction(e.value) - Fraction(e.bound), Fraction(e.value) + Fraction(e.bound)


# ═════════════════════════════════════════════════════════════════════════════
# R5c.1 — Exact: the correctly rounded double and its rounding
# ═════════════════════════════════════════════════════════════════════════════

_EXACT = (
    ("one-third", "Rational(1, 3)"),
    ("two-sevenths", "Rational(2, 7)"),
    ("three-quarters-a-double", "Rational(3, 4)"),
    ("golden-ratio-quadratic-irrational", "(1 + sqrt(5)) / 2"),
    ("a-2-group-k-like-root", "(Rational(37, 20) + sqrt(Rational(301, 400) + Rational(3, 50))) / 2"),
)


@pytest.mark.parametrize("expression", [e for _, e in _EXACT], ids=[n for n, _ in _EXACT])
def test_r5c_1_exact_encloses_its_expression_with_the_rounding(expression: str) -> None:
    r"""The enclosure's value is the correctly rounded double of the exact
    expression and its bound is :math:`|{\rm fl}(x) - x|` rounded up, checked
    against mpmath at 120 digits: the interval contains x, the bound is at
    most half an ulp of the value rounded up (so at most one ulp), and it is
    0 exactly when x is a double (3/4)."""
    import mpmath as mp
    import sympy

    mp.mp.dps = 120
    expr = sympy.sympify(expression)
    e = s5.exact(expr).enclosure()
    x = mp.mpf(sympy.N(expr, 130))
    require(abs(mp.mpf(e.value) - x) <= mp.mpf(e.bound), f"{expression}: {e!r} misses it")
    require(e.value == float(x), f"{expression}: value {e.value!r} is not the correctly rounded {float(x)!r}")
    require(e.bound <= math.ulp(e.value), f"{expression}: bound {e.bound!r} > one ulp")
    require((e.bound == 0.0) == (mp.mpf(e.value) == x), f"{expression}: bound {e.bound!r} on an exact double is not 0, or the converse")


def test_r5c_1_two_spellings_sympy_canonicalises_alike_are_one_value() -> None:
    """``Rational(2, 6)`` and ``Rational(1, 3)`` parse to one SymPy object, so
    their stored ``srepr`` text is one value. NOT claimed: that every two
    spellings of one real number are one value (``sqrt(4)`` and ``2`` are, by
    SymPy's evaluation; a sum SymPy does not simplify is not); identity is by
    the canonical text, a cache miss never a wrong hit."""
    import sympy

    E = s5.c(s5.CERTIFICATE, "Exact")
    a, b = E("Rational(2, 6)", "x"), E(sympy.srepr(sympy.Rational(1, 3)), "x")
    require(a == b and content_digest(a) == content_digest(b), f"{a.expression!r} vs {b.expression!r}")


_BAD_EXACT = (
    ("not-srepr", ("1/3", "x"), ValueError),
    ("a-symbol", ("Symbol('k')", "x"), ValueError),
    ("not-finite", ("oo", "x"), ValueError),
    ("not-real", ("I", "x"), ValueError),
    ("a-call-outside-sympy", ("__import__('os')", "x"), ValueError),
    ("empty-by", ("Rational(1, 3)", ""), ValueError),
    ("a-float-expression", (0.5, "x"), TypeError),
)


@pytest.mark.parametrize("args, error", [b[1:] for b in _BAD_EXACT], ids=[b[0] for b in _BAD_EXACT])
def test_r5c_1_exact_admission(args: tuple[Any, Any], error: type[BaseException]) -> None:
    """Refused: text that is not ``srepr`` (parsed through a whitelist, never a
    bare ``eval``), a free symbol, a non-finite or non-real value, an empty
    provenance, a bare float."""
    with pytest.raises(error):
        s5.c(s5.CERTIFICATE, "Exact")(*args)


# ═════════════════════════════════════════════════════════════════════════════
# R5c.2 — DerivedBound; R5c.3 — the claim against its target
# ═════════════════════════════════════════════════════════════════════════════


def test_r5c_2_a_derived_bound_is_its_enclosure_and_its_method() -> None:
    d = s5.derived(1.25, 1e-9, "Nystrom residual bound")
    require(d.enclosure() == s5.enclosure(1.25, 1e-9), f"{d.enclosure()!r}")
    with pytest.raises(ValueError, match="method"):
        s5.c(s5.CERTIFICATE, "DerivedBound")(s5.enclosure(1.0, 0.0), " ")
    with pytest.raises(TypeError, match="Enclosure"):
        s5.c(s5.CERTIFICATE, "DerivedBound")((1.0, 1e-9), "m")


_TARGETS = (
    ("bound-below-target", 1e-9, 1e-8, "Valid"),
    ("bound-equal-to-target", 1e-8, 1e-8, "Valid"),
    ("bound-one-ulp-above-target", math.nextafter(1e-8, 1.0), 1e-8, "Invalid"),
    ("bound-far-above-target", 1e-3, 1e-8, "Invalid"),
)


@pytest.mark.parametrize("bound, target, kind", [t[1:] for t in _TARGETS], ids=[t[0] for t in _TARGETS])
def test_r5c_3_a_claim_meets_its_target_or_the_certificate_is_invalid(bound: float, target: float, kind: str) -> None:
    """``Valid`` iff the established bound is at most the target (closed: equal
    is valid); ``Invalid`` names the observable and both numbers."""
    cert = s5.certificate({s5.eigenvalue(): s5.claim(target, s5.derived(1.0, bound))})
    require(s5.state_kind(cert) == kind, f"state {cert.state!r}, expected {kind}")
    if kind == "Invalid":
        require(any("Eigenvalue" in r for r in cert.state.reasons), f"the reasons {cert.state.reasons} do not name the observable")


_BAD_CLAIMS = (
    ("non-positive-target", lambda: s5.claim(0.0, s5.derived()), ValueError),
    ("infinite-target", lambda: s5.claim(math.inf, s5.derived()), ValueError),
    ("established-by-an-enclosure", lambda: s5.claim(1e-8, s5.enclosure(1.0, 0.0)), TypeError),
    ("established-by-a-measured", lambda: s5.claim(1e-8, s5.mod("orpheus.numerics.outcome").Measured(1.0)), TypeError),
    ("no-claims", lambda: s5.certificate({}), ValueError),
    ("claim-key-not-an-observable", lambda: s5.certificate({"k": s5.claim(1e-8, s5.derived())}), TypeError),
)


@pytest.mark.parametrize("build, error", [b[1:] for b in _BAD_CLAIMS], ids=[b[0] for b in _BAD_CLAIMS])
def test_r5c_3_admission(build: Any, error: type[BaseException]) -> None:
    """A target is a positive finite real; an establishment is ``Exact`` or
    ``DerivedBound`` (a bare ``Enclosure`` has no provenance; a production
    ``Measured`` is never a reference's establishment); a certificate claims
    something, about observables."""
    with pytest.raises(error):
        build()


# ═════════════════════════════════════════════════════════════════════════════
# R5c.4 — the anchor law; R5c.5 — admissible anchors
# ═════════════════════════════════════════════════════════════════════════════


def _anchor(text: str, standing: Any = None) -> Any:
    return s5.published(spec=s5.eigen_medium(), printed_map={s5.eigenvalue(): s5.printed(text)}, standing=standing)


_ANCHORS = (
    ("inside", 1.4, 1e-9, "1.40000", "Valid"),
    ("overlapping-at-the-edge-by-1e-12", 1.400005 + 1e-6 - 1e-12, 1e-6, "1.40000", "Valid"),
    ("disjoint-by-one-printed-unit", 1.40002, 1e-9, "1.40000", "Invalid"),
    ("disjoint-below", 1.39990, 1e-6, "1.40000", "Invalid"),
)


@pytest.mark.parametrize("value, bound, text, kind", [a[1:] for a in _ANCHORS], ids=[a[0] for a in _ANCHORS])
def test_r5c_4_the_anchor_must_meet_the_claim(value: float, bound: float, text: str, kind: str) -> None:
    """For every observable both the claim and the anchor hold, the claim's
    enclosure and the printed enclosure (``Printed.enclosure()``) intersect,
    decided exactly; a disjoint pair makes the certificate ``Invalid``,
    naming the observable and the anchor's citation. The edge row overlaps by
    1e-12 (an exactly touching row is not constructible: a printed half unit
    5·10^-k is never a double, so the two ends cannot be made to coincide)."""
    cert = s5.certificate({s5.eigenvalue(): s5.claim(1e-3, s5.derived(value, bound))}, [s5.corroboration(_anchor(text))])
    require(s5.state_kind(cert) == kind, f"state {cert.state!r}, expected {kind}")
    if kind == "Invalid":
        require(any("SoodForsterParsons2003" in r for r in cert.state.reasons), f"the reasons {cert.state.reasons} do not cite the anchor")


_BAD_ANCHORS = (
    ("a-withdrawn-anchor", lambda: s5.corroboration(_anchor("1.40000", s5.withdrawal())), ValueError),
    ("an-empty-independence-note", lambda: s5.corroboration(_anchor("1.40000"), ""), ValueError),
    ("an-enclosure", lambda: s5.corroboration(s5.enclosure(1.4, 1e-6)), TypeError),
    ("a-float", lambda: s5.corroboration(1.4), TypeError),
    ("a-printed-value", lambda: s5.corroboration(s5.printed("1.40000")), TypeError),
    ("a-production-reading", lambda: s5.corroboration(s5.mod("orpheus.numerics.outcome").Measured(1.4)), TypeError),
)


@pytest.mark.parametrize("build, error", [b[1:] for b in _BAD_ANCHORS], ids=[b[0] for b in _BAD_ANCHORS])
def test_r5c_5_only_a_current_publication_anchors(build: Any, error: type[BaseException]) -> None:
    """Step 5's admissible anchor is a CURRENT ``PublishedSolution`` with a
    non-empty independence note (step 6 adds ``ReferenceSolution``, editing
    this row on purpose). A production reading or solution never anchors:
    a fine-mesh production run is not a reference, and Monte Carlo never is
    (#505)."""
    with pytest.raises(error):
        build()


def test_r5c_5_an_anchor_sharing_no_observable_is_refused() -> None:
    """A corroboration that shares no observable with the claims is vacuous
    and reads as evidence (X1): the certificate refuses it at construction."""
    anchor = s5.published(spec=s5.eigen_slab(), printed_map={s5.point(0.5, 0): s5.printed("0.25")})
    with pytest.raises(ValueError, match="no observable"):
        s5.certificate({s5.eigenvalue(): s5.claim(1e-3, s5.derived())}, [s5.corroboration(anchor)])


# ═════════════════════════════════════════════════════════════════════════════
# R5c.6 — refinements falsify, never certify
# ═════════════════════════════════════════════════════════════════════════════


def _members(*pairs: tuple[float, float, float]) -> list[tuple[float, Any]]:
    return [(h, s5.enclosure(v, b)) for h, v, b in pairs]


def test_r5c_6_enclosures_sharing_a_point_leave_the_state_valid() -> None:
    """Three members whose enclosures share the point 1.0 (the exact answer of
    a model ``J_h = 1 + h^2`` with bounds ``2 h^2``): ``Valid``. The observed
    orders are returned (two, both 2 to rounding), and they decide nothing."""
    r = s5.refinement(s5.eigenvalue(), _members((0.4, 1.16, 0.32), (0.2, 1.04, 0.08), (0.1, 1.01, 0.02)))
    cert = s5.certificate({s5.eigenvalue(): s5.claim(0.05, s5.derived(1.01, 0.02))}, refinements=[r])
    require(s5.state_kind(cert) == "Valid", f"state {cert.state!r}")
    orders = tuple(r.observed_orders())
    require(len(orders) == 1 or len(orders) == 2, f"observed orders {orders}")
    require(all(abs(o - 2.0) < 1e-9 for o in orders), f"observed orders {orders}")


def test_r5c_6_a_disjoint_pair_makes_the_certificate_invalid() -> None:
    """The coarsest member claims [1.15, 1.17], the finest [1.0, 1.02]: no common
    point, so one of the bounds is wrong and the certificate is ``Invalid``,
    naming the observable."""
    r = s5.refinement(s5.eigenvalue(), _members((0.4, 1.16, 0.01), (0.2, 1.04, 0.08), (0.1, 1.01, 0.01)))
    cert = s5.certificate({s5.eigenvalue(): s5.claim(0.05, s5.derived(1.01, 0.01))}, refinements=[r])
    require(s5.state_kind(cert) == "Invalid", f"state {cert.state!r}")


def test_r5c_6_a_claim_outside_the_refinements_common_part_is_invalid() -> None:
    """The members agree with each other ([0.96, 1.12] and [0.99, 1.03] share
    [0.99, 1.03]) while the claim [1.09, 1.11] meets its target and misses that
    common part: the exact answer lies in the common part if the members'
    bounds are correct, so the claim's bound is wrong, and the certificate is
    ``Invalid``, naming the observable and the missed common part. Control: the
    same refinement with a claim [1.00, 1.04] inside it is ``Valid`` (the
    members' agreement alone does not make the row red)."""
    r = s5.refinement(s5.eigenvalue(), _members((0.2, 1.04, 0.08), (0.1, 1.01, 0.02)))
    require(r.common_part() is not None, "activation: the members share no point, so the row tests the other law")
    missed = s5.certificate({s5.eigenvalue(): s5.claim(0.05, s5.derived(1.10, 0.01))}, refinements=[r])
    require(s5.state_kind(missed) == "Invalid", f"state {missed.state!r}")
    require(any("Eigenvalue" in reason and "misses the refinement's common part" in reason for reason in missed.state.reasons),
            f"the reasons {missed.state.reasons} do not name the missed common part")
    inside = s5.certificate({s5.eigenvalue(): s5.claim(0.05, s5.derived(1.02, 0.02))}, refinements=[r])
    require(s5.state_kind(inside) == "Valid", f"control: state {inside.state!r}")


def test_r5c_6_a_perfect_order_never_certifies_and_a_wild_one_never_refutes() -> None:
    """The ruling's two sides. (i) A sequence with the theoretical order 2
    exactly cannot make ``Valid`` a claim whose bound exceeds its target.
    (ii) A sequence whose observed order is nonsense (values moving away,
    then back) but whose enclosures share a point leaves a meeting claim
    ``Valid``: the order is information, not a verdict."""
    perfect = s5.refinement(s5.eigenvalue(), _members((0.4, 1.16, 0.32), (0.2, 1.04, 0.08), (0.1, 1.01, 0.02)))
    over = s5.certificate({s5.eigenvalue(): s5.claim(1e-3, s5.derived(1.01, 0.02))}, refinements=[perfect])
    require(s5.state_kind(over) == "Invalid", "a perfect order certified a bound above its target")
    wild = s5.refinement(s5.eigenvalue(), _members((0.4, 1.0, 0.5), (0.2, 1.3, 0.4), (0.1, 1.01, 0.02)))
    meets = s5.certificate({s5.eigenvalue(): s5.claim(0.05, s5.derived(1.01, 0.02))}, refinements=[wild])
    require(s5.state_kind(meets) == "Valid", f"a wild observed order refuted: {meets.state!r}")


_BAD_REFINEMENTS = (
    ("one-member", lambda: s5.refinement(s5.eigenvalue(), _members((0.1, 1.0, 0.1)))),
    ("parameters-not-decreasing", lambda: s5.refinement(s5.eigenvalue(), _members((0.1, 1.0, 0.1), (0.2, 1.0, 0.1)))),
    ("a-member-not-an-enclosure", lambda: s5.refinement(s5.eigenvalue(), [(0.2, 1.0), (0.1, s5.enclosure(1.0, 0.1))])),
)


@pytest.mark.parametrize("build", [b for _, b in _BAD_REFINEMENTS], ids=[n for n, _ in _BAD_REFINEMENTS])
def test_r5c_6_admission(build: Any) -> None:
    with pytest.raises((TypeError, ValueError)):
        build()


def test_r5c_6_a_refinement_of_an_unclaimed_observable_is_refused() -> None:
    r = s5.refinement(s5.point(0.5, 0), _members((0.2, 1.0, 0.1), (0.1, 1.0, 0.05)))
    with pytest.raises(ValueError, match="claim"):
        s5.certificate({s5.eigenvalue(): s5.claim(1e-3, s5.derived())}, refinements=[r])


# ═════════════════════════════════════════════════════════════════════════════
# R5c.7 — the standing; R5c.8 — the state is derived; R5c.9 — no ladder certifies
# ═════════════════════════════════════════════════════════════════════════════


def test_r5c_7_a_withdrawn_certificate_is_withdrawn_whatever_its_evidence() -> None:
    """A ``Withdrawal`` standing is the state, IS the same object, on a
    certificate whose evidence alone would be ``Valid`` and on one that
    would be ``Invalid``."""
    w = s5.withdrawal(issue=516)
    good = s5.certificate({s5.eigenvalue(): s5.claim(1e-3, s5.derived(1.0, 1e-9))}, standing=w)
    bad = s5.certificate({s5.eigenvalue(): s5.claim(1e-12, s5.derived(1.0, 1e-9))}, standing=w)
    require(good.state is w and bad.state is w, f"{good.state!r}, {bad.state!r}")


def test_r5c_8_the_state_is_derived_never_stored() -> None:
    """``dataclasses.fields(ReferenceCertificate)`` is exactly ``(claims,
    corroborations, refinements, standing)``: no stored state, verdict or
    validity flag (a stored verdict drifts from its evidence)."""
    names = tuple(f.name for f in dataclasses.fields(s5.c(s5.CERTIFICATE, "ReferenceCertificate")))
    require(names == ("claims", "corroborations", "refinements", "standing"), f"fields {names}")


_STRUCK = ("ConvergedLadder", "Extrapolated", "RichardsonEstimate", "LadderBound")


def test_r5c_9_no_ladder_certifies() -> None:
    """The user's ruling: no establishment by convergence. By AST over every
    ``.py`` under ``orpheus/`` (input count printed, more than 300): none of
    the struck classes is defined; positive control: ``DerivedBound`` is
    found in ``orpheus/reference/certificate.py``. The establishment sum names
    exactly ``{Exact, DerivedBound}``."""
    import typing

    files = sorted((_ROOT / "orpheus").rglob("*.py"))
    print(f"R5c.9: {len(files)} files parsed")
    require(len(files) > 300, f"activation: {len(files)}")
    found: dict[str, list[str]] = {}
    for path in files:
        for node in ast.walk(ast.parse(path.read_text(), filename=str(path))):
            if isinstance(node, ast.ClassDef):
                found.setdefault(node.name, []).append(str(path.relative_to(_ROOT)))
    require(found.get("DerivedBound") == ["orpheus/reference/certificate.py"], f"control: {found.get('DerivedBound')}")
    struck = {n: found[n] for n in _STRUCK if n in found}
    require(not struck, f"struck establishment classes defined: {struck}")
    members = {k.__name__ for k in typing.get_args(s5.c(s5.CERTIFICATE, "Establishment"))}
    require(members == {"Exact", "DerivedBound"}, f"Establishment names {sorted(members)}")


# ═════════════════════════════════════════════════════════════════════════════
# R5c.10 — content identity
# ═════════════════════════════════════════════════════════════════════════════


def _cert(**changes: Any) -> Any:
    args = dict(claims={s5.eigenvalue(): s5.claim(1e-3, s5.derived(1.4, 1e-9))},
                corroborations=[s5.corroboration(_anchor("1.40000"))], refinements=[], standing=None)
    args.update(changes)
    return s5.certificate(args["claims"], args["corroborations"], args["refinements"], args["standing"])


ROSTER: tuple[Entry, ...] = (
    Entry(cls=s5.c(s5.CERTIFICATE, "Exact"), base=s5.exact, parts=("expression", "by"),
          perturb={"expression": (leg("another rational", lambda: s5.exact("2/7")),),
                   "by": (leg("another identity", lambda: s5.exact(by="derive_other")),)},
          pairs=(("two builds", s5.exact, s5.exact),)),
    Entry(cls=s5.c(s5.CERTIFICATE, "DerivedBound"), base=s5.derived, parts=("derived", "method"),
          perturb={"derived": (leg("bound one ULP", lambda: s5.derived(bound=math.nextafter(1e-9, 1.0))),),
                   "method": (leg("another method", lambda: s5.derived(method="series tail")),)},
          pairs=(("two builds", s5.derived, s5.derived),)),
    Entry(cls=s5.c(s5.CERTIFICATE, "Claim"), base=lambda: s5.claim(1e-3, s5.derived()), parts=("target", "established"),
          perturb={"target": (leg("one ULP", lambda: s5.claim(math.nextafter(1e-3, 1.0), s5.derived())),),
                   "established": (leg("exact instead", lambda: s5.claim(1e-3, s5.exact())),)},
          pairs=(("int-like target", lambda: s5.claim(1, s5.derived()), lambda: s5.claim(1.0, s5.derived())),)),
    Entry(cls=s5.c(s5.CERTIFICATE, "Corroboration"), base=lambda: s5.corroboration(_anchor("1.40000")),
          parts=("anchor", "independence"),
          perturb={"anchor": (leg("another printed digit", lambda: s5.corroboration(_anchor("1.40001"))),),
                   "independence": (leg("another note", lambda: s5.corroboration(_anchor("1.40000"), "another method")),)},
          pairs=(("two builds", lambda: s5.corroboration(_anchor("1.40000")), lambda: s5.corroboration(_anchor("1.40000"))),)),
    Entry(cls=s5.c(s5.CERTIFICATE, "Refinement"),
          base=lambda: s5.refinement(s5.eigenvalue(), _members((0.2, 1.04, 0.08), (0.1, 1.01, 0.02))),
          parts=("observable", "members"),
          perturb={"observable": (leg("a point value", lambda: s5.refinement(s5.point(), _members((0.2, 1.04, 0.08), (0.1, 1.01, 0.02)))),),
                   "members": (leg("one value by one ULP", lambda: s5.refinement(s5.eigenvalue(), _members((0.2, math.nextafter(1.04, 2.0), 0.08), (0.1, 1.01, 0.02)))),)},
          pairs=(("two builds", lambda: s5.refinement(s5.eigenvalue(), _members((0.2, 1.04, 0.08), (0.1, 1.01, 0.02))),
                  lambda: s5.refinement(s5.eigenvalue(), _members((0.2, 1.04, 0.08), (0.1, 1.01, 0.02)))),)),
    Entry(cls=s5.c(s5.CERTIFICATE, "Valid"), base=lambda: s5.c(s5.CERTIFICATE, "Valid")(), parts=(), perturb={},
          pairs=(("two builds", lambda: s5.c(s5.CERTIFICATE, "Valid")(), lambda: s5.c(s5.CERTIFICATE, "Valid")()),)),
    Entry(cls=s5.c(s5.CERTIFICATE, "Invalid"), base=lambda: s5.c(s5.CERTIFICATE, "Invalid")(("a reason",)), parts=("reasons",),
          perturb={"reasons": (leg("another reason", lambda: s5.c(s5.CERTIFICATE, "Invalid")(("another reason",))),)},
          pairs=(("two builds", lambda: s5.c(s5.CERTIFICATE, "Invalid")(("a reason",)), lambda: s5.c(s5.CERTIFICATE, "Invalid")(("a reason",))),)),
    Entry(cls=s5.c(s5.CERTIFICATE, "ReferenceCertificate"), base=_cert,
          parts=("claims", "corroborations", "refinements", "standing"),
          perturb={"claims": (leg("another target", lambda: _cert(claims={s5.eigenvalue(): s5.claim(2e-3, s5.derived(1.4, 1e-9))})),),
                   "corroborations": (leg("none", lambda: _cert(corroborations=[])),),
                   "refinements": (leg("one", lambda: _cert(refinements=[s5.refinement(s5.eigenvalue(), _members((0.2, 1.4, 1e-3), (0.1, 1.4, 1e-6)))])),),
                   "standing": (leg("withdrawn", lambda: _cert(standing=s5.withdrawal())),)},
          pairs=(("two builds", _cert, _cert),)),
)


@pytest.mark.parametrize("entry", ROSTER, ids=[e.id for e in ROSTER])
def test_r5c_10_the_population_is_the_parts(entry: Entry) -> None:
    check_population(entry)


@pytest.mark.parametrize("entry, part, the_leg", perturbation_ids(ROSTER),
                         ids=[param_id(e.id, p, lg[0]) for e, p, lg in perturbation_ids(ROSTER)])
def test_r5c_10_each_part_moves_the_digest(entry: Entry, part: str, the_leg: Any) -> None:
    check_perturbation(entry, part, the_leg)


@pytest.mark.parametrize("entry, pair", pair_ids(ROSTER), ids=[param_id(e.id, p[0]) for e, p in pair_ids(ROSTER)])
def test_r5c_10_equal_content_is_one_value(entry: Entry, pair: Any) -> None:
    check_equal_pair(entry, pair)


@pytest.mark.parametrize("entry", ROSTER, ids=[e.id for e in ROSTER])
def test_r5c_10_pickle(entry: Entry) -> None:
    check_pickle(entry)
