r"""The reference solution and its reading (#405 P2 step 6, R6.1-R6.8; step 7b.1 re-posed R6.5).

Specified by the test-architect (2026-10-03, ``.claude/plans/reference_p2_spec.md``
§1.6), on the user's rulings of 2026-10-03 (a reference's reading is its
natural extension, one answer; no stored representation in P2; certification
decoupled from caching) and qa F2 of step 5:

* ``ReferenceSolution(specification, derivation, certificate)``: the
  derivation (a ``Derivation``, ``evaluate(observable) -> Evaluation``, with
  ``Evaluation = Exact | DerivedBound | Uncertified``: the natural extension
  evaluated with its derived bound where the family has one, on demand, for
  any admissible observable) is a plain reference; the certificate is a
  ``ReferenceCertificate`` or ``None``;
* construction enforces what a certificate cannot: every claimed observable is
  posable on the specification (``admit_observable``), every anchor prints for
  the SAME specification, and no claim is on a ``Ratio`` (a ratio's reading is
  the quotient of its operands' readings, one definition);
* ``read(observable)``: ``admit_observable`` first; a ratio by
  ``Enclosure.__truediv__`` of its operands' readings; an eigenvalue or a
  linear observable is EVALUATED by the derivation: an ``Exact`` or a
  ``DerivedBound`` reads as its enclosure (a derived enclosure is a guarantee
  on its own, G3), an ``Uncertified`` reads as itself; when the certificate
  claims the observable, the claim's enclosure and the established one must
  share a point (refused, "disagrees", otherwise), and a claimed observable
  the derivation leaves uncertified is refused ("uncertified");
* step 7b.1 (the user's ruling of 2026-10-03, R6.5 re-posed, R6.8): a family
  that cannot derive a bound reads ``Uncertified(value)``, a reading with no
  guarantee, which only the explicit uncertified comparison accepts.

The rows posed on the exact infinite medium (R6.6 and the exact family's R6.8
rows) state claims about that derivation, so they live with it, in
``tests/gates/derivations/test_exact_medium_reference.py`` (the layer-DAG
re-filing, 2026-10-26). What stays here is the reference layer's own: the
reader over the test-side ``TableDerivation``, and ``Uncertifiable``'s working
precision on ``Exact`` itself.
"""

from __future__ import annotations

import ast
import os
import subprocess
import sys
from fractions import Fraction
from pathlib import Path
from typing import Any

import pytest

from tests.gates._content_identity_helpers import require
from tests.gates.reference import _step5 as s5
from tests.gates.reference import _step6 as s6

pytestmark = pytest.mark.foundation

_ROOT = Path(__file__).resolve().parents[3]


# ═════════════════════════════════════════════════════════════════════════════
# R6.1 — construction
# ═════════════════════════════════════════════════════════════════════════════


def test_r6_1_the_fields_and_their_types() -> None:
    """The three fields, all required (no default on the certificate: a family
    without a derived bound passes ``None`` and says so at the call); the
    test-side table derivation satisfies the ``Derivation`` protocol; a
    certificate of ``None`` constructs."""
    import dataclasses
    import inspect

    cls = s5.c(s6.SOLUTION, "ReferenceSolution")
    names = tuple(f.name for f in dataclasses.fields(cls))
    require(names == ("specification", "derivation", "certificate"), f"fields {names}")
    defaulted = [p.name for p in inspect.signature(cls).parameters.values() if p.default is not inspect.Parameter.empty]
    require(not defaulted, f"defaulted fields {defaulted}")
    spec, derivation, cert = s6.medium_table()
    require(isinstance(derivation, s6.derivation_protocol()), "the table derivation is not a Derivation")
    s6.reference(spec, derivation, cert)
    s6.reference(spec, derivation, None)


class _RetiredVerb:
    """A derivation spelling step 6's verb ``establish``, retired by step 7b.1 (``evaluate``): not a ``Derivation``."""

    def establish(self, observable: Any) -> Any:
        return s6.exact_of(Fraction(7, 5))


_BAD = (
    ("specification-a-string", lambda d, c: s6.reference("medium", d, c), TypeError, "specification"),
    ("derivation-without-evaluate", lambda d, c: s6.reference(s5.eigen_medium(), object(), c), TypeError, "derivation"),
    ("derivation-spelling-the-retired-verb", lambda d, c: s6.reference(s5.eigen_medium(), _RetiredVerb(), c), TypeError, "derivation"),
    ("certificate-a-mapping", lambda d, c: s6.reference(s5.eigen_medium(), d, {}), TypeError, "certificate"),
)


@pytest.mark.parametrize("build, error, fragment", [b[1:] for b in _BAD], ids=[b[0] for b in _BAD])
def test_r6_1_admission(build: Any, error: type[BaseException], fragment: str) -> None:
    _, derivation, cert = s6.medium_table()
    with pytest.raises(error, match=fragment):
        build(derivation, cert)


def test_r6_1_a_claim_the_specification_cannot_pose_is_refused() -> None:
    """qa F2: a certificate claiming a ``PointValue`` cannot be held by an
    infinite-medium reference (the medium has no position); the refusal is
    ``admit_observable``'s (R6.2's route covers the function)."""
    _, derivation, _ = s6.medium_table()
    cert = s5.certificate({s5.point(0.5, 0): s6.exact_claim(Fraction(1, 2))})
    with pytest.raises(ValueError, match="position"):
        s6.reference(s5.eigen_medium(), derivation, cert)


def test_r6_1_an_anchor_of_another_specification_is_refused() -> None:
    """qa F2, ``[M]`` before this step: a certificate claiming k on the medium,
    anchored by a publication printing for a DIFFERENT specification (the
    slab), read ``Valid``. The holder refuses it, naming the specification;
    the same anchor printed for the holder's own specification is admitted."""
    spec, derivation, _ = s6.medium_table()
    foreign = s5.published(spec=s5.eigen_slab(), printed_map={s5.eigenvalue(): s5.printed("1.40000")})
    own = s5.published(spec=spec, printed_map={s5.eigenvalue(): s5.printed("1.40000")})
    claims = {s5.eigenvalue(): s6.exact_claim(Fraction(7, 5))}
    with pytest.raises(ValueError, match="specification"):
        s6.reference(spec, derivation, s5.certificate(claims, [s5.corroboration(foreign)]))
    s6.reference(spec, derivation, s5.certificate(claims, [s5.corroboration(own)]))


def test_r6_1_a_ratio_claim_is_refused() -> None:
    """PROPOSAL (§1.6): a ratio's reading is the quotient of its operands'
    readings, so a claim on a ``Ratio`` would be a second definition of one
    number; the holder refuses it, naming the quotient."""
    spec, derivation, _ = s6.medium_table()
    ratio = s5.c(s5.OBSERVABLE, "Ratio")(s6.group_flux(0), s6.group_flux(1))
    cert = s5.certificate({ratio: s6.exact_claim(Fraction(390, 77))})
    with pytest.raises(ValueError, match="quotient"):
        s6.reference(spec, derivation, cert)


# ═════════════════════════════════════════════════════════════════════════════
# R6.2 — read admits first, through the one function
# ═════════════════════════════════════════════════════════════════════════════

_UNPOSABLE = (
    ("point-value-on-the-medium", lambda: s5.point(0.5, 0)),
    ("weight-with-two-regions", lambda: s6.group_flux(0, regions=2)),
    ("weight-with-three-groups", lambda: s6.group_flux(0, groups=3)),
)


@pytest.mark.parametrize("make", [m for _, m in _UNPOSABLE], ids=[n for n, _ in _UNPOSABLE])
def test_r6_2_an_unposable_observable_is_refused_before_the_derivation_runs(make: Any) -> None:
    spec, derivation, cert = s6.medium_table()
    ref = s6.reference(spec, derivation, cert)
    derivation.calls.clear()
    with pytest.raises(ValueError):
        ref.read(make())
    require(not derivation.calls, f"the derivation ran on an unposable observable: {derivation.calls}")


def test_r6_2_read_admits_through_the_specifications_one_function(monkeypatch: pytest.MonkeyPatch) -> None:
    """ROUTE gate: with ``admit_observable`` rebound in every binding to a decoy,
    a posable read raises the decoy's error (the reader calls it) and so does
    an unposable one (no second admission in the reader refuses first)."""
    import importlib

    spec, derivation, cert = s6.medium_table()
    ref = s6.reference(spec, derivation, cert)
    module_name, name = s5.ADMISSION
    honest = getattr(importlib.import_module(module_name), name)

    def decoy(*args: Any, **kwargs: Any) -> None:
        raise ValueError("DECOY-ADMISSION")

    rebound = 0
    for mod in list(sys.modules.values()):
        if mod is not None and getattr(mod, name, None) is honest:
            monkeypatch.setattr(mod, name, decoy)
            rebound += 1
    require(rebound >= 1, "activation: nothing bound admit_observable")
    for observable in (s5.eigenvalue(), s5.point(0.5, 0)):
        with pytest.raises(ValueError, match="DECOY-ADMISSION"):
            ref.read(observable)


# ═════════════════════════════════════════════════════════════════════════════
# R6.3 — a read is the derivation's established enclosure, checked against a claim
# ═════════════════════════════════════════════════════════════════════════════


@pytest.mark.parametrize("which", ["eigenvalue", "group-0", "group-1"])
def test_r6_3_a_read_is_the_derivations_established_enclosure(which: str) -> None:
    """The reading equals ``derivation.evaluate(observable).enclosure()``, and
    the derivation was asked for exactly that observable, on demand."""
    spec, derivation, cert = s6.medium_table()
    ref = s6.reference(spec, derivation, cert)
    observable = {"eigenvalue": s5.eigenvalue(), "group-0": s6.group_flux(0), "group-1": s6.group_flux(1)}[which]
    derivation.calls.clear()
    reading = ref.read(observable)
    require(derivation.calls == [observable], f"{which}: the derivation was called with {derivation.calls}")
    require(reading == derivation.evaluate(observable).enclosure(), f"{which}: read {reading!r}")


def test_r6_3_the_reading_is_the_established_enclosure_not_the_claims() -> None:
    """When the claim (exact 7/5) and the established enclosure (a derived bound
    1.4 +- 1e-6) agree, the reading is the ESTABLISHED one: the certificate
    corroborates, it does not replace the derivation's answer."""
    spec, derivation, cert = s6.medium_table()
    derived = s5.derived(1.4, 1e-6, "a residual bound")
    derivation.table[s5.eigenvalue()] = derived
    reading = s6.reference(spec, derivation, cert).read(s5.eigenvalue())
    require(reading == derived.enclosure(), f"read {reading!r}, established {derived.enclosure()!r}")


def test_r6_3_an_unclaimed_observable_reads_its_established_enclosure() -> None:
    """G2's point: an observable the certificate never declared (a weight
    (1, 1), a test mesh's cell average in general) reads the derivation's
    enclosure, with and without a certificate."""
    spec, derivation, cert = s6.medium_table()
    both = s5.flux(((1.0, 1.0),))
    derivation.table[both] = Fraction(30, 7) + Fraction(11, 13)
    for certificate in (cert, None):
        reading = s6.reference(spec, derivation, certificate).read(both)
        require(s6.contains(reading, Fraction(30, 7) + Fraction(11, 13)), f"{reading!r}")


def test_r6_3_a_derivation_that_disagrees_with_its_claim_is_refused() -> None:
    """A derivation establishing k = 7/5 + 1/10^9 exactly, against an exact claim
    7/5: the two enclosures share no point (the step-5 family law), so the
    reference cannot be BUILT, naming the disagreement. The claim and the
    establishment are one quantity in two places (X4), so they are checked once,
    at construction, and ``read`` is a pure evaluation (the elegance review of
    step 6, finding 4)."""
    spec, derivation, cert = s6.medium_table()
    derivation.table[s5.eigenvalue()] = Fraction(7, 5) + Fraction(1, 10**9)
    with pytest.raises(ValueError, match="disagree"):
        s6.reference(spec, derivation, cert)


@pytest.mark.parametrize("standing", ["invalid", "withdrawn"])
def test_r6_3_a_non_valid_certificate_still_reads(standing: str) -> None:
    """Ruled 2026-10-03: ``read`` answers on an ``Invalid`` or ``Withdrawn``
    certificate (only the ``VerificationCertificate`` refuses it)."""
    spec, derivation, _ = s6.medium_table()
    if standing == "invalid":
        cert = s5.certificate({s5.eigenvalue(): s5.claim(1e-30, s6.exact_of(Fraction(7, 5)))})
        require(s5.state_kind(cert) == "Invalid", f"activation: {cert.state!r}")
    else:
        cert = s5.certificate({s5.eigenvalue(): s6.exact_claim(Fraction(7, 5))}, standing=s5.withdrawal())
    reading = s6.reference(spec, derivation, cert).read(s5.eigenvalue())
    require(s6.contains(reading, Fraction(7, 5)), f"{reading!r}")


# ═════════════════════════════════════════════════════════════════════════════
# R6.4 — a ratio is the quotient of its operands' readings
# ═════════════════════════════════════════════════════════════════════════════


def test_r6_4_a_ratio_reads_through_the_one_quotient(monkeypatch: pytest.MonkeyPatch) -> None:
    """ROUTE gate (X4): with ``Enclosure.__truediv__`` and ``Ratio.quotient``
    each wrapped by a counter, a ratio's reading equals ``read(a) / read(b)``,
    the division ran once and it ran inside the one ratio rule,
    ``Ratio.quotient`` (the reader keeps no private quotient: improved
    2026-10-26, when R7.19, the cross-answer witness of that rule, moved to
    the homogeneous tests); the quotient's enclosure contains the exact ratio
    (30/7)/(11/13)."""
    Enclosure = s5.c(s5.ENCLOSURE, "Enclosure")
    Ratio = s5.c(s5.OBSERVABLE, "Ratio")
    honest = Enclosure.__truediv__
    honest_quotient = Ratio.quotient
    calls: list[Any] = []
    quotients: list[Any] = []

    def counting(self: Any, other: Any) -> Any:
        calls.append((self, other))
        return honest(self, other)

    def counting_quotient(self: Any, read: Any) -> Any:
        quotients.append(self)
        return honest_quotient(self, read)

    spec, derivation, cert = s6.medium_table()
    ref = s6.reference(spec, derivation, cert)
    a, b = ref.read(s6.group_flux(0)), ref.read(s6.group_flux(1))
    monkeypatch.setattr(Enclosure, "__truediv__", counting)
    monkeypatch.setattr(Ratio, "quotient", counting_quotient)
    reading = ref.read(Ratio(s6.group_flux(0), s6.group_flux(1)))
    require(len(calls) == 1, f"the quotient ran {len(calls)} times")
    require(len(quotients) == 1, f"the ratio rule Ratio.quotient ran {len(quotients)} times")
    require(reading == honest(a, b), f"{reading!r} != {honest(a, b)!r}")
    require(s6.contains(reading, Fraction(30, 7) / Fraction(11, 13)), f"{reading!r} misses the exact ratio")


# ═════════════════════════════════════════════════════════════════════════════
# R6.5 — no derivable bound: the reading is Uncertified (re-posed, step 7b.1)
# ═════════════════════════════════════════════════════════════════════════════

#: An observable the table derivation leaves uncertified (no derived bound), and a second one.
_U_OBS = ((0.5, 2.0),)
_U_VALUE = 2.5
_V_OBS = ((2.0, 0.5),)
_V_VALUE = 0.75


def _uncertified_table() -> tuple[Any, Any, Any]:
    """``medium_table`` with two observables the derivation leaves uncertified."""
    spec, derivation, cert = s6.medium_table()
    derivation.table[s5.flux(_U_OBS)] = s6.uncertified(_U_VALUE)
    derivation.table[s5.flux(_V_OBS)] = s6.uncertified(_V_VALUE)
    return spec, derivation, cert


@pytest.mark.parametrize("certified", [True, False], ids=["with-a-certificate", "without"])
def test_r6_5_an_observable_the_derivation_cannot_bound_reads_uncertified(certified: bool) -> None:
    """The user's ruling of 2026-10-03 (step 7b.1): a family with no derived
    bound for an observable evaluates it ``Uncertified(value)``, and ``read``
    returns that value AS IT IS: an ``Uncertified`` (never an ``Enclosure``,
    so no float without a bound passes for a guarantee), equal to what
    ``evaluate`` returned, the derivation asked once. With and without a
    certificate (the certificate does not claim the observable)."""
    spec, derivation, cert = _uncertified_table()
    ref = s6.reference(spec, derivation, cert if certified else None)
    derivation.calls.clear()
    reading = ref.read(s5.flux(_U_OBS))
    require(derivation.calls == [s5.flux(_U_OBS)], f"the derivation was called with {derivation.calls}")
    require(type(reading) is s6.uncertified_class(), f"read {reading!r}, not an Uncertified")
    require(reading == s6.uncertified(_U_VALUE) and reading.value == _U_VALUE, f"read {reading!r}")
    require(not hasattr(reading, "enclosure"), "an uncertified reading offers an enclosure")


_RATIOS = (
    ("uncertified-over-enclosed", (_U_OBS, None), lambda a, b: _U_VALUE / b),
    ("enclosed-over-uncertified", (None, _U_OBS), lambda a, b: a / _U_VALUE),
    ("uncertified-over-uncertified", (_U_OBS, _V_OBS), lambda a, b: _U_VALUE / _V_VALUE),
)


@pytest.mark.parametrize("operands, expected", [r[1:] for r in _RATIOS], ids=[r[0] for r in _RATIOS])
def test_r6_5_a_ratio_with_an_uncertified_operand_reads_uncertified(operands: Any, expected: Any) -> None:
    """A ratio reads as the quotient of its operands' readings (R6.4), so a
    ratio with an uncertified operand, in either position, is uncertified:
    ``Uncertified(a.value / b.value)``, where an enclosed operand contributes
    its double (here the group-0 flux 30/7, exact, read as an enclosure).
    ROUTE: exactly one division ran, an ``Uncertified`` dunder."""
    spec, derivation, cert = _uncertified_table()
    ref = s6.reference(spec, derivation, cert)
    enclosed = ref.read(s6.group_flux(0))
    require(type(enclosed) is s5.c(s5.ENCLOSURE, "Enclosure"), f"activation: group 0 reads {enclosed!r}")
    num, den = (s5.flux(o) if o is not None else s6.group_flux(0) for o in operands)
    U = s6.uncertified_class()
    calls: list[str] = []
    honest = {name: getattr(U, name) for name in ("__truediv__", "__rtruediv__")}
    with pytest.MonkeyPatch.context() as mp:
        for name, method in honest.items():
            def counting(self: Any, other: Any, _name: str = name, _method: Any = method) -> Any:
                calls.append(_name)
                return _method(self, other)
            mp.setattr(U, name, counting)
        reading = ref.read(s5.c(s5.OBSERVABLE, "Ratio")(num, den))
    require(len(calls) == 1, f"the Uncertified quotient ran {calls}")
    want = expected(enclosed.value, enclosed.value)
    require(type(reading) is U and reading.value == want, f"read {reading!r}, expected Uncertified({want!r})")


@pytest.mark.parametrize("claimed", ["contains-the-value", "disagrees-with-the-value"])
def test_r6_5_a_certificate_cannot_claim_an_uncertified_observable(claimed: str) -> None:
    """A certificate claim is a target a derived bound is held to; an
    observable the derivation leaves uncertified has no bound to hold, so the
    reference cannot be built: ``ValueError`` naming the uncertified
    evaluation. The check runs BEFORE the claim-versus-establishment agreement:
    a claim that would also DISAGREE (exact 7/5 against the value 2.5) is
    refused as uncertified, not as a disagreement (the wiring order, lessons
    ``L43c``); a claim whose enclosure contains the value is refused too (the
    value would otherwise pass for corroborated)."""
    spec, derivation, _ = _uncertified_table()
    target = s5.derived(_U_VALUE, 1e-3, "a test bound") if claimed == "contains-the-value" else s6.exact_of(Fraction(7, 5))
    cert = s5.certificate({s5.flux(_U_OBS): s5.claim(1.0, target)})
    with pytest.raises(ValueError, match="uncertified") as refusal:
        s6.reference(spec, derivation, cert)
    require("disagree" not in str(refusal.value), f"refused for the agreement, not the uncertified evaluation: {refusal.value}")


_NOT_EVALUATIONS = (
    ("an-enclosure", lambda: s5.enclosure(1.4, 0.0)),
    ("a-printed-value", lambda: s5.printed("1.40000")),
    ("a-measured", lambda: _measured(1.4)),
)


def _measured(value: float) -> Any:
    from orpheus.numerics.outcome import Measured

    return Measured(value)


@pytest.mark.parametrize("make", [m for _, m in _NOT_EVALUATIONS], ids=[n for n, _ in _NOT_EVALUATIONS])
def test_r6_5_what_a_derivation_returns_is_admitted_as_an_evaluation(make: Any) -> None:
    """A ``Derivation`` is an open protocol, so the reader admits what it
    returns: an object outside ``Evaluation = Exact | DerivedBound |
    Uncertified`` (an ``Enclosure``, a ``Printed``, a ``Measured``) is refused,
    ``TypeError`` naming the evaluation, never read as a reading. The sum
    itself names exactly the three members."""
    import typing

    members = set(typing.get_args(s6.evaluation()))
    want = {s5.c(s5.CERTIFICATE, "Exact"), s5.c(s5.CERTIFICATE, "DerivedBound"), s6.uncertified_class()}
    require(members == want, f"Evaluation names {sorted(m.__qualname__ for m in members)}")
    spec, derivation, cert = s6.medium_table()
    derivation.table[s5.flux(_U_OBS)] = make()
    with pytest.raises(TypeError, match="evaluation"):
        s6.reference(spec, derivation, cert).read(s5.flux(_U_OBS))


# ═════════════════════════════════════════════════════════════════════════════
# R6.7 — the layer
# ═════════════════════════════════════════════════════════════════════════════

_LAYER_SCRIPT = """
import sys
import orpheus.reference.solution as m
print(m.__file__)
print(sorted({k.split('.')[1] for k in sys.modules if k.startswith('orpheus.')}))
"""


def test_r6_7_the_reader_imports_no_derivation() -> None:
    """A cold import of ``orpheus.reference.solution`` loads only
    ``{reference, specification, numerics, data, geometry}``: the derivations
    implement the protocol and import the reference package, never the
    reverse."""
    env = {**os.environ, "PYTHONPATH": str(_ROOT)}
    out = subprocess.run([sys.executable, "-O", "-c", _LAYER_SCRIPT], cwd=_ROOT, env=env, capture_output=True, text=True)
    require(out.returncode == 0, out.stderr)
    file, packages = out.stdout.strip().splitlines()
    require(file.startswith(str(_ROOT)), f"the subprocess imported {file} (L22)")
    loaded = set(ast.literal_eval(packages))
    require("reference" in loaded, f"activation: {loaded}")
    allowed = {"reference", "specification", "numerics", "data", "geometry", "values"}
    require(loaded <= allowed, f"loads {sorted(loaded - allowed)}")



# ═════════════════════════════════════════════════════════════════════════════
# R6.8 — Uncertifiable carries the value at the working precision (step 7b.1)
# ═════════════════════════════════════════════════════════════════════════════


#: A rational standing in for the group flux the exact medium supplied until the layer-DAG re-filing
#: (2026-10-26): the cancellation in ``r z + r 10^-k`` depends on ``r``, and 17/23 is a value at which
#: SymPy's default 15 digits are wrong at both offsets below while the working precision is right
#: (``[M]`` 2026-10-26, a scan of p/q < 30: 5 of the fractions scanned had this property; at 17/23 the
#: 15-digit value is 5.0e3 and 8.0e34 relative off at k = 130 and 160, the 60-digit value 3.3e-17 and
#: 7.2e-17), so the working precision is what the rows see.
_SCALE = Fraction(17, 23)


def _offset_approximation(exponent: int) -> Any:
    """``Uncertifiable.approximation`` of ``r (z + 10^-k)``, z the Machin zero and r ``_SCALE``.

    SymPy distributes the rational over the sum, so the expression is ``r z + r 10^-k``,
    the shape a derivation's flux integral of the weight ``z + 10^-k`` has. ``Exact``
    refusing it as ``Uncertifiable`` is the activation (strict ``evalf`` cancels).
    """
    import sympy

    expression = sympy.Rational(_SCALE.numerator, _SCALE.denominator) * (s6.machin_zero() + sympy.Integer(10) ** -exponent)
    with pytest.raises(s5.c(s5.CERTIFICATE, "Uncertifiable")) as refusal:
        s5.c(s5.CERTIFICATE, "Exact")(sympy.srepr(expression), "the working precision")
    return refusal.value.approximation


def _relative_error(approximation: float, exponent: int) -> Fraction:
    truth = _SCALE * Fraction(1, 10**exponent)
    return abs(Fraction(approximation) - truth) / truth


_OFFSETS = (("at-1e-130", 130), ("at-1e-160", 160))


@pytest.mark.parametrize("exponent", [e for _, e in _OFFSETS], ids=[n for n, _ in _OFFSETS])
def test_r6_8_the_uncertified_value_is_evaluated_at_the_working_precision(exponent: int) -> None:
    """The fallback's DEFINITION, not a guarantee (the value stays
    uncertified): the value an uncertified reading carries is the expression's
    non-strict ``evalf`` at the working precision ``_EXACT_DIGITS`` (60),
    carried by ``Uncertifiable.approximation``. On ``r z + r 10^-k`` (r =
    17/23, z the Machin zero) the evaluation at ``evalf``'s default 15 digits
    is all cancellation (``[M]`` 2026-10-26: 5.0e3 times off at k = 130 and
    8.0e34 times at k = 160); at 60 digits both carry the value within 1e-15
    relative of the exact ``r 10^-k`` (``[M]`` 3.3e-17 and 7.2e-17). (Re-posed
    2026-10-26 from the exact medium's group-1 flux onto ``Exact`` itself: the
    precision is the certificate's; that a derivation reads this value is
    ``tests/gates/derivations/test_exact_medium_reference.py``, R6.8.)"""
    relative = _relative_error(_offset_approximation(exponent), exponent)
    require(relative <= Fraction(1, 10**15), f"the uncertified value is {float(relative):.3g} relative off")


@pytest.mark.xfail(strict=True, reason=(
    "a FIXED working precision has an offset beyond which cancellation swamps the value: [M] 2026-10-26 the "
    "expression (17/23)(z + 10^-250) carries an approximation 4.7e58 times off at 60 digits (2026-10-03, on the "
    "exact medium's group-0 flux: 2.8e59). Ruled out of scope for 7b.1 (2026-10-03): no precision schedule earns "
    "a guarantee; the fix is to CERTIFY such values with an absolute radius from evalf's accuracy record, #568 "
    "(P4), when this row XPASSes and is rewritten to assert the value is certified"))
def test_r6_8_the_uncertified_value_survives_a_deep_cancellation() -> None:
    """CHALLENGE: the same definition as the row above at k = 250: the
    approximation is the expression's value within 1e-15 relative (the value
    stays uncertified; the row pins the fallback's definition, not a
    guarantee)."""
    relative = _relative_error(_offset_approximation(250), 250)
    require(relative <= Fraction(1, 10**15), f"the uncertified value is {float(relative):.3g} relative off")


def test_r6_8_uncertifiable_is_the_one_cannot_certify_refusal() -> None:
    """``Uncertifiable`` is a ``ValueError`` (so the step-5 rows that match a
    ``ValueError`` still hold) raised by ``Exact`` on the Machin zero; a
    non-real expression raises a plain ``ValueError`` (the discrimination)."""
    import sympy

    U = s5.c(s5.CERTIFICATE, "Uncertifiable")
    require(issubclass(U, ValueError), "Uncertifiable is not a ValueError")
    with pytest.raises(U, match="certified digits"):
        s5.c(s5.CERTIFICATE, "Exact")(sympy.srepr(s6.machin_zero()), "activation")
    with pytest.raises(ValueError) as plain:
        s5.c(s5.CERTIFICATE, "Exact")(sympy.srepr(sympy.acos(2)), "a non-real constant")
    require(not isinstance(plain.value, U), f"a non-real expression is refused as uncertifiable: {plain.value}")
