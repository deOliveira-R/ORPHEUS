r"""The verification certificate and the comparison verbs (#405 P2 step 7a, R7.1-R7.10; step 7b.1, R7.11-R7.18).

Specified by the test-architect (2026-10-03, ``.claude/plans/reference_p2_spec.md``
§1.7a), on the design "P2, the design as ruled" item 7, the rulings of §0 item 3
(the algebraic-error floor is required for an ORDER verdict only) and the
carried review note S4 (the verbs call ``answer.read(observable)`` themselves:
the reading verb is the only producer of a production reading):

* ``verify_agreement(answer, observable, reference, tolerance, algebraic_error)``
  returns a ``VerificationCertificate``: the reference must be held ``Valid``
  (refused BEFORE anything is read); the floor ``b_ref <= tol / 10``; the
  verdict ``|m - v_ref| + b_ref <= tol``, both decided EXACTLY; the algebraic
  error recorded, never required;
* ``verify_order(answers, observable, reference, tolerance, algebraic_errors,
  order, band)`` returns an ``OrderVerification``: every algebraic error
  ``Measured`` or ``Asserted`` at most ``tol / 10`` (``NotYet`` and
  ``NotApplicable`` refused as unestablished), every error resolved above the
  tolerance, the observed orders returned and held against the caller's order
  and band;
* the first production reading: ``HomogeneousResult.read``.

* step 7b.1 (the user's ruling of 2026-10-03): a reference reading an
  observable ``Uncertified`` (its family derives no bound) is refused by both
  verbs ("uncertified", before production is read); the explicit, weaker
  claim is ``compare_uncertified(answer, observable, reference, tolerance)``,
  an ``UncertifiedComparison`` (``|m - v| <= tol``, decided exactly; not a
  certificate), which refuses a non-``Valid`` certificate before any read
  and a CERTIFIED reading (``ReadingCertified``, naming ``verify_agreement``).

DECLARED LIMIT: nothing pairs the answer with the reference's specification
(production results hold none until P4's projection); the caller pairs them.
"""

from __future__ import annotations

import ast
import dataclasses
import math
import os
import subprocess
import sys
from fractions import Fraction
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from orpheus.numerics.outcome import Asserted, Measured, NotApplicable
from tests.gates._content_identity_helpers import require
from tests.gates.reference import _step5 as s5
from tests.gates.reference import _step6 as s6
from tests.gates.reference import _step7 as s7
from tests.gates.specification import _fixtures as sf

pytestmark = pytest.mark.foundation

_HERE = "tests/gates/reference/test_verification.py"
_ROOT = Path(__file__).resolve().parents[3]
_K = Fraction(7, 5)  # not a double, so its exact enclosure has a positive bound


# ═════════════════════════════════════════════════════════════════════════════
# R7.1 — the reference is held Valid, before anything is read
# ═════════════════════════════════════════════════════════════════════════════

_NOT_VALID = (
    ("no-certificate", "none", "no certificate"),
    ("invalid", "invalid", "invalid"),
    ("withdrawn", "withdrawn", "4321"),
)


@pytest.mark.parametrize("certificate, fragment", [r[1:] for r in _NOT_VALID], ids=[r[0] for r in _NOT_VALID])
def test_r7_1_a_reference_not_held_valid_is_refused_before_any_reading(certificate: str, fragment: str) -> None:
    """``ReferenceNotValid`` (a ``ValueError``), naming why: no certificate; an
    ``Invalid`` one (its reason); a ``Withdrawn`` one (its issue). Neither
    answer is read: a refused reference must not cost a production solve's
    reading, and the order of refusals is the contract."""
    answer = s7.FixedAnswer({s5.eigenvalue(): 1.4})
    reference = s7.reference_of(_K, certificate)
    if certificate == "invalid":
        require(s5.state_kind(reference.certificate) == "Invalid", f"activation: {reference.certificate.state!r}")
    reference.derivation.calls.clear()
    error = s7.not_valid()
    require(issubclass(error, ValueError), "ReferenceNotValid is not a ValueError")
    with pytest.raises(error, match=fragment):
        s7.verify_agreement(answer, s5.eigenvalue(), reference, 1e-3, s7.NOT_YET)
    require(not answer.reads and not reference.derivation.calls, f"read before refusing: {answer.reads}, {reference.derivation.calls}")


def test_r7_1_the_order_verb_holds_the_reference_valid_too() -> None:
    answers = [(0.2, s7.FixedAnswer({s5.eigenvalue(): 1.44})), (0.1, s7.FixedAnswer({s5.eigenvalue(): 1.41}))]
    with pytest.raises(s7.not_valid(), match="no certificate"):
        s7.verify_order(answers, s5.eigenvalue(), s7.reference_of(_K, "none"), 1e-3,
                        [Measured(0.0), Measured(0.0)], 2.0, 0.2)


# ═════════════════════════════════════════════════════════════════════════════
# R7.2 — the roles: the verb reads the answer itself
# ═════════════════════════════════════════════════════════════════════════════


def test_r7_2_the_verb_calls_the_answers_read_exactly_once() -> None:
    """S4: the verb is handed the ANSWER and reads it; the certificate's
    reading is what ``read`` returned (identity), read once."""
    reading = Measured(1.4)
    answer = s7.FixedAnswer({s5.eigenvalue(): reading})
    require(isinstance(answer, s7.answer_protocol()), "the double is not a ProductionAnswer")
    cert = s7.verify_agreement(answer, s5.eigenvalue(), s7.reference_of(_K), 1e-3, s7.NOT_YET)
    require(answer.reads == [s5.eigenvalue()], f"reads {answer.reads}")
    require(cert.reading is reading, f"the certificate holds {cert.reading!r}, not the answer's reading")


_BAD_READINGS = (
    ("an-enclosure", lambda: s5.enclosure(1.4, 0.0)),
    ("a-printed-value", lambda: s5.printed("1.40000")),
    ("a-bare-float", lambda: 1.4),
    ("asserted-evidence", lambda: Asserted(1e-9, "a convergence claim")),
)


@pytest.mark.parametrize("make", [m for _, m in _BAD_READINGS], ids=[n for n, _ in _BAD_READINGS])
def test_r7_2_a_reading_that_is_not_a_production_reading_is_refused(make: Any) -> None:
    """An answer whose ``read`` returns a reference's reading (an ``Enclosure``,
    a ``Printed``), a bare number, or another ``Evidence`` member is refused,
    ``TypeError`` naming the production reading: a reference cannot be
    verified against itself through a production slot (G3)."""
    class Returns:
        def read(self, observable: Any) -> Any:
            return make()

    answer = Returns()
    with pytest.raises(TypeError, match="production reading"):
        s7.verify_agreement(answer, s5.eigenvalue(), s7.reference_of(_K), 1e-3, s7.NOT_YET)


def test_r7_2_an_object_without_read_is_not_an_answer() -> None:
    """A ``Measured`` (a reading, not an answer) is refused as the answer."""
    with pytest.raises(TypeError):
        s7.verify_agreement(Measured(1.4), s5.eigenvalue(), s7.reference_of(_K), 1e-3, s7.NOT_YET)


_BAD_ARGUMENTS = (
    ("algebraic-error-a-float", dict(algebraic_error=1e-9), TypeError, "algebraic error"),
    ("algebraic-error-none", dict(algebraic_error=None), TypeError, "algebraic error"),
    ("tolerance-zero", dict(tolerance=0.0), ValueError, "tolerance"),
    ("tolerance-infinite", dict(tolerance=math.inf), ValueError, "tolerance"),
    ("tolerance-a-string", dict(tolerance="1e-3"), TypeError, "tolerance"),
)


@pytest.mark.parametrize("change, error, fragment", [r[1:] for r in _BAD_ARGUMENTS], ids=[r[0] for r in _BAD_ARGUMENTS])
def test_r7_2_the_arguments(change: dict[str, Any], error: type[BaseException], fragment: str) -> None:
    args = dict(answer=s7.FixedAnswer({s5.eigenvalue(): 1.4}), observable=s5.eigenvalue(),
                reference=s7.reference_of(_K), tolerance=1e-3, algebraic_error=s7.NOT_YET)
    args.update(change)
    with pytest.raises(error, match=fragment):
        s7.verify_agreement(**args)


def test_r7_2_no_argument_has_a_default() -> None:
    """Every call site states its evidence (§0 item 3): no default on either verb."""
    import inspect

    for name in ("verify_agreement", "verify_order"):
        params = inspect.signature(s5.c(s7.VERIFICATION, name)).parameters.values()
        defaulted = [p.name for p in params if p.default is not inspect.Parameter.empty]
        require(not defaulted, f"{name} defaults {defaulted}")


# ═════════════════════════════════════════════════════════════════════════════
# R7.3 — the reference floor; R7.4 — the agreement verdict, exactly
# ═════════════════════════════════════════════════════════════════════════════

_TOL = 2.0**-16
_B = 2.0**-20  # the reference bound, dyadic: every row below is exactly representable


def _cert(reading: float, bound: float = _B, tolerance: float = _TOL) -> Any:
    answer = s7.FixedAnswer({s5.eigenvalue(): reading})
    return s7.verify_agreement(answer, s5.eigenvalue(), s7.reference_with_bound(1.0, bound), tolerance, s7.NOT_YET)


def test_r7_3_the_floor_is_a_tenth_of_the_tolerance_inclusive() -> None:
    """``floor_holds`` iff the reference bound is at most a tenth of the
    tolerance, decided exactly: at ``tol = 10 * b`` it holds (10 * 2^-20 is a
    double), one ULP of tolerance below it does not; a failed floor makes
    ``agrees`` false even for an exact reading, and ``require`` raises for the
    FLOOR, not the reading (the order ``AgreementCertificate`` fixed)."""
    at = _cert(1.0, _B, 10 * _B)
    require(at.floor_holds and at.agrees, f"at the floor: {at!r}")
    below = _cert(1.0, _B, math.nextafter(10 * _B, 0.0))
    require(not below.floor_holds and not below.agrees, f"below the floor: {below!r}")
    with pytest.raises(AssertionError, match="floor"):
        below.require()


_VERDICTS = (
    ("exactly-at-the-tolerance", 1.0 + _TOL - _B, True),
    ("one-ulp-beyond", math.nextafter(1.0 + _TOL - _B, 2.0), False),
    ("inside-but-for-the-reference-bound", 1.0 + _TOL, False),
    ("below-the-reference", 1.0 - _TOL + _B, True),
    ("exact", 1.0, True),
)


@pytest.mark.parametrize("reading, agrees", [v[1:] for v in _VERDICTS], ids=[v[0] for v in _VERDICTS])
@pytest.mark.rests_on(f"{_HERE}::test_r7_3_the_floor_is_a_tenth_of_the_tolerance_inclusive")
def test_r7_4_agreement_is_the_distance_plus_the_reference_bound(reading: float, agrees: bool) -> None:
    """``agrees`` iff ``|m - v_ref| + b_ref <= tol``, exactly: the boundary row
    agrees, one ULP beyond does not, and a reading within ``tol`` of the
    reference's centre but not of every point of its enclosure does not (the
    worst case over the enclosure, the verdict a theorem about the exact
    answer). ``require`` raises "disagrees" exactly when it does not agree."""
    cert = _cert(reading)
    require(cert.floor_holds, "activation: the floor fails")
    require(cert.agrees is agrees, f"agrees = {cert.agrees}, expected {agrees}")
    if agrees:
        cert.require()
    else:
        with pytest.raises(AssertionError, match="disagrees"):
            cert.require()


# ═════════════════════════════════════════════════════════════════════════════
# R7.5 — the algebraic error is recorded under agreement
# ═════════════════════════════════════════════════════════════════════════════


@pytest.mark.parametrize("evidence", [s7.NOT_YET, s7.DIRECT, Measured(1e-14), Asserted(1e-12, "a claim")],
                         ids=["not-yet", "not-applicable", "measured", "asserted"])
def test_r7_5_agreement_records_any_algebraic_evidence(evidence: Any) -> None:
    """The comparison bounds the total error whatever its split (§0 item 3), so
    every ``Evidence`` member gives a verdict and is carried unchanged."""
    answer = s7.FixedAnswer({s5.eigenvalue(): 1.0})
    cert = s7.verify_agreement(answer, s5.eigenvalue(), s7.reference_with_bound(1.0, _B), _TOL, evidence)
    require(cert.algebraic_error is evidence and cert.agrees, f"{cert!r}")


# ═════════════════════════════════════════════════════════════════════════════
# R7.6 — the order verdict
# ═════════════════════════════════════════════════════════════════════════════


def _ladder(order: float, hs: tuple[float, ...] = (0.4, 0.2, 0.1), c: float = 0.5) -> list[tuple[float, Any]]:
    """Production answers ``1 + c h^order`` against an exact reference 1."""
    return [(h, s7.FixedAnswer({s5.eigenvalue(): 1.0 + c * h**order})) for h in hs]


def _order(answers: Any, *, errors: Any = None, order: float = 2.0, band: float = 0.1, tolerance: float = 1e-3) -> Any:
    errors = [Measured(0.0)] * len(answers) if errors is None else errors
    return s7.verify_order(answers, s5.eigenvalue(), s7.reference_with_bound(1.0, 0.0), tolerance, errors, order, band)


def test_r7_6_the_observed_orders_are_returned_and_held_against_the_declared_order() -> None:
    """A second-order ladder holds against order 2 (band 0.1): the observed
    orders are returned, one INTERVAL per consecutive pair (the errors are
    known within the reference bound plus the algebraic error, qa of step
    7a), each containing 2; a first-order ladder does not hold against order 2."""
    second = _order(_ladder(2.0))
    orders = tuple(second.observed_orders)
    require(len(orders) == 2 and all(low <= 2.0 + 1e-9 and 2.0 - 1e-9 <= high for low, high in orders), f"observed orders {orders}")
    require(second.holds, f"{second!r}")
    first = _order(_ladder(1.0))
    require(not first.holds, f"a first-order ladder held against order 2: {tuple(first.observed_orders)}")


_UNESTABLISHED = (
    ("not-yet", [s7.NOT_YET] * 3),
    ("not-applicable", [s7.DIRECT] * 3),
    ("measured-above-the-floor", [Measured(0.0), Measured(2e-4), Measured(0.0)]),
    ("asserted-above-the-floor", [Asserted(2e-4, "c"), Measured(0.0), Measured(0.0)]),
)


@pytest.mark.parametrize("errors", [u[1] for u in _UNESTABLISHED], ids=[u[0] for u in _UNESTABLISHED])
def test_r7_6_the_order_verdict_requires_the_algebraic_floor(errors: Any) -> None:
    """§0 item 3: attributing the error to the discretisation needs every rung's
    algebraic error ``Measured`` or ``Asserted`` at most ``tol / 10``;
    ``NotYet`` and ``NotApplicable`` cannot establish it. Refused,
    ``Unestablished`` (a ``ValueError``) naming the algebraic error. Positive leg:
    ``Measured`` and ``Asserted`` at the floor hold (the row above, and
    ``Asserted(5e-5)`` here)."""
    with pytest.raises(s7.unestablished(), match="algebraic error"):
        _order(_ladder(2.0), errors=errors)
    _order(_ladder(2.0), errors=[Asserted(5e-5, "c"), Measured(5e-5), Measured(0.0)])


def test_r7_6_an_error_below_the_tolerance_is_unresolved() -> None:
    """An order observed on errors below the tolerance measures the floors, not
    the discretisation: a rung whose error is below ``tol`` is refused,
    "unresolved" (``c = 1e-3`` puts the finest rung's error at 1e-5)."""
    with pytest.raises(ValueError, match="unresolved"):
        _order(_ladder(2.0, c=1e-3))


def test_r7_6_the_order_verdict_holds_the_reference_floor() -> None:
    answers = _ladder(2.0)
    with pytest.raises(s7.unestablished(), match="reference's bound"):
        s7.verify_order(answers, s5.eigenvalue(), s7.reference_with_bound(1.0, 2e-4), 1e-3, [Measured(0.0)] * 3, 2.0, 0.1)


def test_r7_6_the_order_verdict_reads_each_answer_once() -> None:
    answers = _ladder(2.0)
    _order(answers)
    require(all(a.reads == [s5.eigenvalue()] for _, a in answers), f"{[a.reads for _, a in answers]}")


# ═════════════════════════════════════════════════════════════════════════════
# R7.7 — returned, never stored
# ═════════════════════════════════════════════════════════════════════════════


def test_r7_7_the_verdicts_are_returned_values_not_stored_ones() -> None:
    """``VerificationCertificate`` and ``OrderVerification`` are frozen
    dataclasses and NOT ``ContentIdentity`` (no cache keys them); the verdicts
    ``floor_holds``, ``agrees`` and ``holds`` are DERIVED properties, not
    fields; two calls return two objects."""
    from orpheus.numerics.content import ContentIdentity

    for name, derived in (("VerificationCertificate", ("floor_holds", "agrees")), ("OrderVerification", ("holds",))):
        cls: type = s5.c(s7.VERIFICATION, name)
        params = getattr(cls, "__dataclass_params__", None)
        require(params is not None and params.frozen, f"{name} is not a frozen dataclass")
        require(not issubclass(cls, ContentIdentity), f"{name} is a content value")
        fields = {f.name for f in dataclasses.fields(cls)}
        require(not fields & set(derived), f"{name} stores {sorted(fields & set(derived))}")
        for prop in derived:
            require(isinstance(getattr(cls, prop, None), property), f"{name}.{prop} is not a property")
    a, b = _cert(1.0), _cert(1.0)
    require(a is not b, "the verb returned a stored object")


# ═════════════════════════════════════════════════════════════════════════════
# R7.8 — the first production reading: HomogeneousResult.read
# ═════════════════════════════════════════════════════════════════════════════


def _flux(weights: Any) -> Any:
    return s5.flux((tuple(weights),))


def test_r7_8_the_homogeneous_result_reads_its_own_answer() -> None:
    # The ONLY guard of production's weight reading against a hand computation:
    # production and the reference read a weight through one shared function,
    # values_without_position, so R7.9's agreement cannot see that function
    # reverse the groups (qa of step 7a, finding 4).
    """``read`` returns ``Measured``: the eigenvalue is ``k_inf``; a group's
    indicator flux integral is that group's flux; a weighted flux integral is
    the weighted sum (per unit volume, the result's gauge); a ratio is the
    quotient of its operands' readings; a constant ``Symbolic`` weight reads
    as its table twin. 2 groups (a 1-group fixture cannot tell a transposed
    weight, lessons ``L86d``)."""
    from orpheus.numerics.mesh_free_function import Symbolic

    result = s7.solve_homogeneous(sf.fuel())
    k = result.read(s5.eigenvalue())
    require(type(k) is Measured and k.value == result.k_inf, f"{k!r}")
    for g in range(2):
        phi = result.read(s6.group_flux(g))
        require(type(phi) is Measured and phi.value == float(result.flux[g]), f"group {g}: {phi!r}")
    w = (0.25, 3.0)
    weighted = result.read(_flux(w)).value
    require(abs(weighted - (w[0] * result.flux[0] + w[1] * result.flux[1])) <= 4 * math.ulp(weighted), f"{weighted!r}")
    ratio = result.read(s5.c(s5.OBSERVABLE, "Ratio")(s6.group_flux(0), s6.group_flux(1))).value
    require(ratio == float(result.flux[0]) / float(result.flux[1]), f"ratio {ratio!r}")
    constant = s5.c(s5.OBSERVABLE, "FluxIntegral")(Symbolic.of(0.25, 3))
    require(abs(result.read(constant).value - weighted) <= 4 * math.ulp(weighted), "a constant Symbolic weight reads differently from its table")


_REFUSED_READS = (
    ("point-value", lambda: s5.point(0.5, 0), "position"),
    ("three-groups", lambda: _flux((1.0, 1.0, 1.0)), "group"),
    ("two-regions", lambda: s6.group_flux(0, regions=2), "region"),
)


@pytest.mark.parametrize("make, fragment", [r[1:] for r in _REFUSED_READS], ids=[r[0] for r in _REFUSED_READS])
def test_r7_8_what_a_zero_dimensional_answer_cannot_read_is_refused(make: Any, fragment: str) -> None:
    result = s7.solve_homogeneous(sf.fuel())
    with pytest.raises(ValueError, match=fragment):
        result.read(make())


# ═════════════════════════════════════════════════════════════════════════════
# R7.9 — the first verification: production homogeneous against the exact medium
# ═════════════════════════════════════════════════════════════════════════════

def _library_a_4g() -> Any:
    from orpheus.derivations.common.xs_library import get_mixture

    return get_mixture("A", "4g")


# fuel(4) has four identical groups, so it cannot see a group reversal; library
# mixture A in 4 groups has distinct group fluxes (qa of step 7a, finding 3)
_MIXTURES = (("fuel-2g", sf.fuel), ("fission-only-2g", sf.fission_only), ("fuel-4g", lambda: sf.fuel(4)),
             ("library-A-4g", _library_a_4g))
#: The tolerance, relative: ``[M]`` 2026-10-03 the production k and group fluxes
#: sit within 1.65e-16 relative of the exact ones over these three mixtures
#: (600x margin), and the exact reference's bound (half an ulp, ~1.1e-16
#: relative) meets the floor tol/10 with 90x margin.
_REL_TOL = 1e-13


def _medium(mixture: Any) -> Any:
    from orpheus.data.cells import CellCoefficient, Channel
    from orpheus.numerics.question import Eigen
    from orpheus.specification import InfiniteMediumSpecification

    # The fission gauge is DECLARED (the user's ruling of 2026-10-08: the default production counts the (n,2n)
    # emission too): these rows judge the homogeneous solver, fission-gauged at the density 100 until #517 moves it onto the declaration.
    fission = CellCoefficient.every(Channel.FISSION_EMISSION)
    return InfiniteMediumSpecification(0, mixture, Eigen(fission, gauge=fission))


@pytest.mark.parametrize("make", [m for _, m in _MIXTURES], ids=[n for n, _ in _MIXTURES])
def test_r7_9_the_homogeneous_solver_agrees_with_the_exact_medium(make: Any) -> None:
    """REFERENCE (exact rationals) at L1: the production 0-D solve of each
    mixture agrees with the exact infinite medium on k, on each group's flux
    and on the first two groups' ratio, at 1e-13 relative, the algebraic error
    ``NotApplicable`` (a direct solve). This is the step's own positive leg;
    the next row is the instrument's negative one."""
    mixture = make()
    reference = s6.exact_medium_reference(_medium(mixture))
    result = s7.solve_homogeneous(mixture)
    groups = len(result.flux)
    ratio = s5.c(s5.OBSERVABLE, "Ratio")(s6.group_flux(0, groups=groups), s6.group_flux(1, groups=groups))
    for observable in [s5.eigenvalue(), ratio] + [s6.group_flux(g, groups=groups) for g in range(groups)]:
        scale = abs(reference.read(observable).value)
        cert = s7.verify_agreement(result, observable, reference, _REL_TOL * scale, s7.DIRECT)
        cert.require()


def test_r7_9_a_perturbed_production_disagrees() -> None:
    """X1, the negative: production solved with ``nu Sigma_f`` scaled by
    ``1 + 1e-11`` (k moves by 1e-11 relative, 100x the tolerance) DISAGREES
    with the exact reference of the unperturbed mixture, on k."""
    mixture = sf.fuel()
    perturbed = dataclasses.replace(mixture, SigP=np.asarray(mixture.SigP) * (1.0 + 1e-11))
    reference = s6.exact_medium_reference(_medium(mixture))
    cert = s7.verify_agreement(s7.solve_homogeneous(perturbed), s5.eigenvalue(), reference,
                               _REL_TOL * abs(reference.read(s5.eigenvalue()).value), s7.DIRECT)
    require(cert.floor_holds and not cert.agrees, f"{cert!r}")


# ═════════════════════════════════════════════════════════════════════════════
# R7.10 — the layer
# ═════════════════════════════════════════════════════════════════════════════

_LAYER_SCRIPT = """
import sys
import orpheus.reference.verification as m
print(m.__file__)
print(sorted({k.split('.')[1] for k in sys.modules if k.startswith('orpheus.')}))
"""


def test_r7_10_the_verbs_import_no_method() -> None:
    """A cold import of ``reference.verification`` loads only
    ``{reference, specification, numerics, data, geometry}``: the verbs are
    generic over a ``ProductionAnswer``; no method package (homogeneous, sn)
    and no derivation is imported."""
    env = {**os.environ, "PYTHONPATH": str(_ROOT)}
    out = subprocess.run([sys.executable, "-O", "-c", _LAYER_SCRIPT], cwd=_ROOT, env=env, capture_output=True, text=True)
    require(out.returncode == 0, out.stderr)
    file, packages = out.stdout.strip().splitlines()
    require(file.startswith(str(_ROOT)), f"the subprocess imported {file} (L22)")
    loaded = set(ast.literal_eval(packages))
    require("reference" in loaded, f"activation: {loaded}")
    allowed = {"reference", "specification", "numerics", "data", "geometry"}
    require(loaded <= allowed, f"loads {sorted(loaded - allowed)}")


# ═════════════════════════════════════════════════════════════════════════════
# The step-7a review round (the elegance review, 2026-10-03)
# ═════════════════════════════════════════════════════════════════════════════


def test_r7_review_a_certificate_cannot_be_built_against_an_invalid_reference() -> None:
    """The certificate holds the answer and the reference and reads both at
    construction, so building one directly against an Invalid reference is
    refused, not read as "agrees" (built directly, it read agrees=True)."""
    certificate_class = s5.c(s7.VERIFICATION, "VerificationCertificate")
    answer = s7.FixedAnswer({s5.eigenvalue(): 1.4})
    with pytest.raises(s7.not_valid(), match="invalid"):
        certificate_class(answer, s5.eigenvalue(), s7.reference_of(Fraction(7, 5), "invalid"), 1e-6, s7.Measured(0.0))


def test_r7_review_the_failures_are_typed() -> None:
    """require() raises ReferenceTooLoose for the floor and Disagreement for
    the comparison, both AssertionErrors."""
    loose = s7.verify_agreement(s7.FixedAnswer({s5.eigenvalue(): 1.4}), s5.eigenvalue(),
                                s7.reference_with_bound(1.4, 1e-3), 1e-3, s7.Measured(0.0))
    with pytest.raises(s5.c(s7.VERIFICATION, "ReferenceTooLoose")):
        loose.require()
    off = s7.verify_agreement(s7.FixedAnswer({s5.eigenvalue(): 1.5}), s5.eigenvalue(),
                              s7.reference_of(Fraction(7, 5)), 1e-6, s7.Measured(0.0))
    with pytest.raises(s5.c(s7.VERIFICATION, "Disagreement")):
        off.require()


@pytest.mark.parametrize("resolutions, errors", [((0.1, 0.1), (Fraction(1), Fraction(1, 4))),
                                                 ((0.2, 0.0), (Fraction(1), Fraction(1, 4))),
                                                 ((0.4, 0.2, 0.1), (Fraction(1),))],
                         ids=["equal-resolutions", "zero-resolution", "lengths-differ"])
def test_r7_review_an_order_verification_admits_its_refinement(resolutions: Any, errors: Any) -> None:
    """An OrderVerification admits its resolutions as the reference's
    Refinement does (one admission: at least two, positive, decreasing
    strictly) and one error per resolution; equal resolutions crashed in log,
    and 3 resolutions with 1 error held vacuously, before the review."""
    order_class = s5.c(s7.VERIFICATION, "OrderVerification")
    with pytest.raises(ValueError):
        order_class(s5.eigenvalue(), resolutions, errors, tuple(Fraction(0) for _ in errors), 2.0, 0.1)


def test_r7_review_the_order_interval_refuses_a_wrong_order_hidden_by_the_floors() -> None:
    """True order 1.7 at h = 0.2, 0.1 with the finest error 1.2 tol and both
    floors near their limit read 1.963 as a point estimate, inside 2 ± 0.1 (qa
    of step 7a, finding 1). With the errors known only within the floors the
    order interval is wide, and the verdict does not hold."""
    order_class = s5.c(s7.VERIFICATION, "OrderVerification")
    tol = Fraction(1, 10**6)
    fine = Fraction(12, 10) * tol  # the true errors of an order-1.7 method: coarse / fine = 2**1.7
    coarse = fine * Fraction(32490, 10000)
    u = Fraction(999, 1000) * tol / 10 * 2  # the reference bound and the algebraic error, each at its floor
    measured = (coarse + u, fine - u)  # the floors shift each error by up to u, in the worst direction
    point = math.log(float(measured[0] / measured[1])) / math.log(2.0)
    require(1.9 <= point <= 2.1, f"activation: the point estimate {point} should sit inside 2 ± 0.1")
    verdict = order_class(s5.eigenvalue(), (0.2, 0.1), measured, (u, u), 2.0, 0.1)
    require(not verdict.holds, f"a true order 1.7 held against 2 ± 0.1: {verdict.observed_orders}")


def test_r7_review_a_sign_change_is_not_convergence() -> None:
    """Errors +0.04 then -0.01 read as a clean order 2, but the answers cross
    the reference: not monotone convergence, so the verdict does not hold."""
    order_class = s5.c(s7.VERIFICATION, "OrderVerification")
    verdict = order_class(s5.eigenvalue(), (0.2, 0.1), (Fraction(4, 100), Fraction(-1, 100)), (Fraction(0), Fraction(0)), 2.0, 0.1)
    require(not verdict.monotone and not verdict.holds, f"{verdict!r}")
    with pytest.raises(s5.c(s7.VERIFICATION, "OrderNotObserved")):
        verdict.require()


def test_r7_review_the_agreement_verdict_is_exact() -> None:
    """v = 1, b = 2^-60, m = 2, tol = 1: the exact total error 1 + 2^-60 exceeds
    the tolerance, so it disagrees; in floats 1 + 2^-60 rounds to 1 and agreed
    (qa of step 7a, finding 2: every R7.4 row was dyadic)."""
    cert = s7.verify_agreement(s7.FixedAnswer({s5.eigenvalue(): 2.0}), s5.eigenvalue(),
                               s7.reference_with_bound(1.0, 2.0**-60), 1.0, s7.Measured(0.0))
    require(cert.floor_holds and not cert.agrees, f"{cert!r}")


# ═════════════════════════════════════════════════════════════════════════════
# Step 7b.1 — the uncertified reading (the user's ruling of 2026-10-03)
# ═════════════════════════════════════════════════════════════════════════════
#
# R7.11: both verbs refuse an uncertified reading of a Valid reference, before
# production is read. R7.12-R7.18: the explicit uncertified comparison.

_U_K = 1.40625  # dyadic: the uncertified value of k in the fixtures below


def _u_ratio() -> Any:
    """The ratio group-0 flux (exact 30/7) over group-1 flux (uncertified): uncertified through its denominator."""
    return s5.c(s5.OBSERVABLE, "Ratio")(s6.group_flux(0), s6.group_flux(1))


def test_r7_11_verify_agreement_refuses_an_uncertified_reading_of_a_valid_reference() -> None:
    """The reference is ``Valid`` (its certificate claims the group-0 flux) and
    reads k ``Uncertified``: an uncertified value cannot anchor a
    verification, so ``ReadingUncertified`` names it ("uncertified"), through
    the verb and through the certificate built directly, and production is
    NEVER read (a counting double). Activation: the standing is ``Valid`` and
    the reference was read for k (the refusal is the reading's, not the
    standing's). Positive leg: the same reference verifies its certified
    observable, the group-0 flux. Since the review round the refusal is its
    own type, ``ReadingUncertified`` (a ``ValueError``), unrelated to
    ``ReferenceNotValid``, which names the standing only."""
    error = s7.reading_uncertified()
    require(issubclass(error, ValueError), "ReadingUncertified is not a ValueError")
    require(not issubclass(error, s7.not_valid()) and not issubclass(s7.not_valid(), error),
            "ReadingUncertified and ReferenceNotValid are related by subclassing")
    reference = s7.uncertified_reference(_U_K, "valid")
    require(s5.state_kind(reference.certificate) == "Valid", f"activation: {reference.certificate.state!r}")
    for build in (lambda a: s7.verify_agreement(a, s5.eigenvalue(), reference, 1e-3, s7.NOT_YET),
                  lambda a: s5.c(s7.VERIFICATION, "VerificationCertificate")(a, s5.eigenvalue(), reference, 1e-3, s7.NOT_YET)):
        answer = s7.FixedAnswer({s5.eigenvalue(): _U_K})
        reference.derivation.calls.clear()
        with pytest.raises(s7.reading_uncertified(), match="uncertified"):
            build(answer)
        require(not answer.reads, f"production was read before the refusal: {answer.reads}")
        require(s5.eigenvalue() in reference.derivation.calls, f"activation: the reference was not read ({reference.derivation.calls})")
    certified = s7.verify_agreement(s7.FixedAnswer({s6.group_flux(0): float(s7.F0)}), s6.group_flux(0), reference, 1e-12, s7.NOT_YET)
    require(certified.agrees, f"{certified!r}")


def test_r7_11_a_ratio_with_an_uncertified_operand_is_refused_too() -> None:
    """A ratio reads as the quotient of its operands' readings (R6.4/R6.5), so
    a ratio over an uncertified group-1 flux is uncertified and cannot be verified:
    refused "uncertified", production unread."""
    reference = s7.uncertified_reference(_U_K, "valid")
    answer = s7.FixedAnswer({_u_ratio(): 3.0})
    with pytest.raises(s7.reading_uncertified(), match="uncertified"):
        s7.verify_agreement(answer, _u_ratio(), reference, 1e-3, s7.NOT_YET)
    require(not answer.reads, f"production was read: {answer.reads}")


def test_r7_11_verify_order_refuses_an_uncertified_reading_of_a_valid_reference() -> None:
    """The order verb refuses the uncertified reading as the agreement verb
    does ("uncertified"), with every algebraic error established (so no other
    refusal fires first), and reads no answer of the refinement."""
    reference = s7.uncertified_reference(_U_K, "valid")
    answers = [(0.2, s7.FixedAnswer({s5.eigenvalue(): 1.44})), (0.1, s7.FixedAnswer({s5.eigenvalue(): 1.41}))]
    with pytest.raises(s7.reading_uncertified(), match="uncertified"):
        s7.verify_order(answers, s5.eigenvalue(), reference, 1e-3, [Measured(0.0), Measured(0.0)], 2.0, 0.2)
    require(all(not a.reads for _, a in answers), f"production was read: {[a.reads for _, a in answers]}")
    require(s5.eigenvalue() in reference.derivation.calls, "activation: the reference was not read")


# ── R7.12: the standing refusals, before anything is read ───────────────────

_STANDING = (
    ("invalid", "invalid", "invalid"),
    ("withdrawn", "withdrawn", "4321"),
)


@pytest.mark.parametrize("certificate, fragment", [r[1:] for r in _STANDING], ids=[r[0] for r in _STANDING])
def test_r7_12_a_reference_whose_certificate_is_not_valid_is_refused_before_any_reading(certificate: str, fragment: str) -> None:
    """``compare_uncertified`` against a reference whose certificate exists and
    is not ``Valid`` (``Invalid``: its reason; ``Withdrawn``: its issue) raises
    ``ReferenceNotValid`` BEFORE anything is read: neither the reference's
    derivation nor the production answer is called."""
    reference = s7.uncertified_reference(_U_K, certificate)
    if certificate == "invalid":
        require(s5.state_kind(reference.certificate) == "Invalid", f"activation: {reference.certificate.state!r}")
    reference.derivation.calls.clear()
    answer = s7.FixedAnswer({s5.eigenvalue(): _U_K})
    with pytest.raises(s7.not_valid(), match=fragment):
        s7.compare_uncertified(answer, s5.eigenvalue(), reference, 1e-3)
    require(not answer.reads and not reference.derivation.calls,
            f"read before refusing: {answer.reads}, {reference.derivation.calls}")


@pytest.mark.parametrize("certificate", ["none", "valid"])
def test_r7_12_no_certificate_or_a_valid_one_is_compared(certificate: str) -> None:
    """The positive legs: an uncertified family has no certificate (its normal
    state, not a refusal: the verb differs from ``verify_agreement`` here,
    R7.1), and a ``Valid`` certificate that does not claim the observable is
    compared too. The comparison holds what the reference read (equal to the
    derivation's ``Uncertified``) and what the answer read (identity), each
    read once."""
    reference = s7.uncertified_reference(_U_K, certificate)
    reading = Measured(_U_K)
    answer = s7.FixedAnswer({s5.eigenvalue(): reading})
    reference.derivation.calls.clear()
    comparison = s7.compare_uncertified(answer, s5.eigenvalue(), reference, 1e-3)
    require(answer.reads == [s5.eigenvalue()], f"reads {answer.reads}")
    require(reference.derivation.calls == [s5.eigenvalue()], f"the reference was read {reference.derivation.calls}")
    require(comparison.reading is reading, f"the comparison holds {comparison.reading!r}, not the answer's reading")
    require(comparison.reference_reading == s6.uncertified(_U_K), f"{comparison.reference_reading!r}")
    require(type(comparison.reference_reading) is s6.uncertified_class(), "activation: the reference reading's type")
    require(comparison.agrees, f"{comparison!r}")


# ── R7.13: a certified reading is refused: a test is forced up ──────────────

_CERTIFIED = (
    ("valid-certificate-exact-k", lambda: s7.reference_of(_K, "valid"), s5.eigenvalue),
    ("derived-bound", lambda: s7.reference_with_bound(1.0, 1e-9), s5.eigenvalue),
    ("enclosed-flux-of-a-certified-uncertified-family", lambda: s7.uncertified_reference(_U_K, "valid"), lambda: s6.group_flux(0)),
)


@pytest.mark.parametrize("make, observable", [r[1:] for r in _CERTIFIED], ids=[r[0] for r in _CERTIFIED])
def test_r7_13_a_reading_the_reference_certifies_is_refused(make: Any, observable: Any) -> None:
    """An ENCLOSED reading of a reference whose certificate is present (so
    ``Valid``, by the standing check): an ``Exact`` k, a derived bound, the
    exact flux of a family whose other readings are uncertified. It is refused, ``ReadingCertified`` (a ``ValueError``) naming
    ``verify_agreement``: the day a family derives a bound, every
    uncertified comparison against it reds and must be re-posed as a
    verification (the weaker claim cannot outlive its reason). Production is
    not read. ``ReadingCertified`` is not ``ReferenceNotValid`` either way,
    so the two refusals are told apart."""
    error = s7.reading_certified()
    require(issubclass(error, ValueError), "ReadingCertified is not a ValueError")
    require(not issubclass(error, s7.not_valid()) and not issubclass(s7.not_valid(), error),
            "ReadingCertified and ReferenceNotValid are related by subclassing")
    answer = s7.FixedAnswer({observable(): 1.0})
    with pytest.raises(error, match="verify_agreement"):
        s7.compare_uncertified(answer, observable(), make(), 1e-3)
    require(not answer.reads, f"production was read: {answer.reads}")


_NO_CERTIFICATE_ENCLOSED = (
    ("exact-k", lambda: s7.reference_of(_K, "none"), s5.eigenvalue, float(_K)),
    ("exact-flux-of-an-uncertified-family", lambda: s7.uncertified_reference(_U_K, "none"), lambda: s6.group_flux(0), float(s7.F0)),
    ("derived-bound-wider-than-the-tolerance", lambda: s6.reference(s5.eigen_medium(), s6.TableDerivation(
        {s5.eigenvalue(): s5.derived(1.0, 0.25, "a test bound")}), None), s5.eigenvalue, 1.0),
)


@pytest.mark.parametrize("make, observable, value", [r[1:] for r in _NO_CERTIFICATE_ENCLOSED], ids=[r[0] for r in _NO_CERTIFICATE_ENCLOSED])
def test_r7_13_an_enclosed_reading_of_a_reference_without_a_certificate_is_compared(make: Any, observable: Any, value: float) -> None:
    """The review round (2026-10-03): ``verify_agreement`` refuses a
    reference with no certificate (R7.1), so the comparison takes ANY reading
    of one, an ``Enclosure`` included: the two verbs together cover every
    state of a reference in good standing. The comparison holds the
    reference's enclosure (equal to ``read``) and uses only its value: a
    derived bound of 0.25 does not stop a reading within 2^-16 of the value
    from agreeing, and twice the tolerance disagrees."""
    reference = make()
    require(reference.certificate is None, "activation: the reference has a certificate")
    enclosure = reference.read(observable())
    require(type(enclosure) is s5.c(s5.ENCLOSURE, "Enclosure"), f"activation: the reading is {enclosure!r}")
    require(enclosure.value == value, f"activation: the enclosure's value {enclosure.value!r} is not {value!r}")
    tolerance = _UT * value
    near = s7.compare_uncertified(s7.FixedAnswer({observable(): value + tolerance / 2}), observable(), reference, tolerance)
    require(near.reference_reading == enclosure, f"{near.reference_reading!r}")
    require(near.agrees, f"{near!r}")
    beyond = s7.compare_uncertified(s7.FixedAnswer({observable(): value + 2 * tolerance}), observable(), reference, tolerance)
    require(not beyond.agrees, f"{beyond!r}")


def test_r7_13_an_uncertified_ratio_is_compared() -> None:
    """The positive leg of the forcing: a ratio over an uncertified group-1 flux
    reads ``Uncertified(F0 / v)`` (the enclosed numerator contributes its double),
    so it is COMPARED, and the comparison holds that reading."""
    reference = s7.uncertified_reference(_U_K, "valid")
    value = float(s7.F0) / _U_K
    comparison = s7.compare_uncertified(s7.FixedAnswer({_u_ratio(): value}), _u_ratio(), reference, 1e-12)
    require(comparison.reference_reading == s6.uncertified(value), f"{comparison.reference_reading!r}")
    require(comparison.agrees, f"{comparison!r}")


# ── R7.14: the verdict, inclusive and exact; R7.15: Disagreement ───────────

_UT = 2.0**-16
_U_VERDICTS = (
    ("exactly-at-the-tolerance-above", 1.0 + _UT, True),
    ("one-ulp-beyond-above", math.nextafter(1.0 + _UT, 2.0), False),
    ("exactly-at-the-tolerance-below", 1.0 - _UT, True),
    ("one-ulp-beyond-below", math.nextafter(1.0 - _UT, 0.0), False),
    ("equal", 1.0, True),
)


def _compare(reading: float, value: float = 1.0, tolerance: float = _UT) -> Any:
    return s7.compare_uncertified(s7.FixedAnswer({s5.eigenvalue(): reading}), s5.eigenvalue(),
                                  s7.uncertified_reference(value, "none"), tolerance)


@pytest.mark.parametrize("reading, agrees", [v[1:] for v in _U_VERDICTS], ids=[v[0] for v in _U_VERDICTS])
def test_r7_14_the_comparison_agrees_within_the_tolerance_inclusive(reading: float, agrees: bool) -> None:
    """``agrees`` iff ``|m - v| <= tol``, inclusive at the boundary (dyadic
    rows, exactly representable), one ULP beyond on either side refused. No
    reference bound enters (there is none): the row that disagreed under
    ``verify_agreement`` for the bound alone (R7.4, "inside but for the
    reference bound") agrees here, so the weaker claim is a different
    verdict, not a looser tolerance."""
    comparison = _compare(reading)
    require(comparison.agrees is agrees, f"agrees = {comparison.agrees}, expected {agrees}")


def test_r7_14_the_verdict_is_decided_exactly() -> None:
    """v = -0.1 (the double), m = 1e17, tol = 1e17: the exact distance
    1e17 + 0.1000000000000000055... exceeds the tolerance, but its float
    rounds to 1e17 == tol, so a float verdict agrees. The exact verdict does
    not; one ULP more tolerance (1e17 + 16) agrees."""
    require(1e17 - (-0.1) == 1e17, "activation: the float distance rounds onto the tolerance")
    require(not _compare(1e17, -0.1, 1e17).agrees, "the verdict was decided in floats")
    require(_compare(1e17, -0.1, math.nextafter(1e17, math.inf)).agrees, "the control leg disagrees")


def test_r7_15_require_raises_disagreement_exactly_when_it_does_not_agree() -> None:
    """``require()`` raises the step-7a ``Disagreement`` (an
    ``AssertionError``) naming the unverified reference value when the comparison does
    not agree, and returns when it does; there is no floor to raise first
    (the reference has no bound), so ``ReferenceTooLoose`` never fires, even
    with a tolerance far below the value's own rounding."""
    error = s7.disagreement()
    require(issubclass(error, AssertionError), "Disagreement is not an AssertionError")
    with pytest.raises(error, match="unverified reference value"):
        _compare(math.nextafter(1.0 + _UT, 2.0)).require()
    _compare(1.0 + _UT).require()
    _compare(1.0, 1.0, 2.0**-1000).require()


# ── R7.16: the production reading; R7.17: the arguments; R7.18: returned ────

_U_BAD_READINGS = (
    ("an-enclosure", lambda: s5.enclosure(_U_K, 0.0)),
    ("a-printed-value", lambda: s5.printed("1.40625")),
    ("an-uncertified-value", lambda: s6.uncertified(_U_K)),
    ("a-bare-float", lambda: _U_K),
    ("asserted-evidence", lambda: Asserted(1e-9, "a convergence claim")),
)


@pytest.mark.parametrize("make", [m for _, m in _U_BAD_READINGS], ids=[n for n, _ in _U_BAD_READINGS])
def test_r7_16_the_production_reading_is_measured_only(make: Any) -> None:
    """An answer whose ``read`` returns a reference's reading (an
    ``Enclosure``, a ``Printed``, an ``Uncertified``: a reference compared
    with itself through the production slot), a bare number, or another
    ``Evidence`` member is refused, ``TypeError`` naming the production
    reading; an object without ``read`` is not an answer."""
    class Returns:
        def read(self, observable: Any) -> Any:
            return make()

    with pytest.raises(TypeError, match="production reading"):
        s7.compare_uncertified(Returns(), s5.eigenvalue(), s7.uncertified_reference(_U_K), 1e-3)


def test_r7_16_a_reading_is_not_an_answer() -> None:
    with pytest.raises(TypeError):
        s7.compare_uncertified(Measured(_U_K), s5.eigenvalue(), s7.uncertified_reference(_U_K), 1e-3)


_U_BAD_ARGUMENTS = (
    ("tolerance-zero", dict(tolerance=0.0), ValueError, "tolerance"),
    ("tolerance-negative", dict(tolerance=-1e-3), ValueError, "tolerance"),
    ("tolerance-infinite", dict(tolerance=math.inf), ValueError, "tolerance"),
    ("tolerance-a-string", dict(tolerance="1e-3"), TypeError, "tolerance"),
    ("reference-a-table", dict(reference={}), TypeError, "reference"),
)


@pytest.mark.parametrize("change, error, fragment", [r[1:] for r in _U_BAD_ARGUMENTS], ids=[r[0] for r in _U_BAD_ARGUMENTS])
def test_r7_17_the_arguments(change: dict[str, Any], error: type[BaseException], fragment: str) -> None:
    args = dict(answer=s7.FixedAnswer({s5.eigenvalue(): _U_K}), observable=s5.eigenvalue(),
                reference=s7.uncertified_reference(_U_K), tolerance=1e-3)
    args.update(change)
    with pytest.raises(error, match=fragment):
        s7.compare_uncertified(**args)


def test_r7_17_no_argument_has_a_default() -> None:
    """The verb and the comparison's constructor default nothing: the
    tolerance a migrated test keeps is stated at its call."""
    import inspect

    for target in (s5.c(s7.VERIFICATION, "compare_uncertified"), s5.c(s7.VERIFICATION, "UncertifiedComparison")):
        params = list(inspect.signature(target).parameters.values())
        names = [p.name for p in params]
        require(names == ["answer", "observable", "reference", "tolerance"], f"{target.__name__} takes {names}")
        defaulted = [p.name for p in params if p.default is not inspect.Parameter.empty]
        require(not defaulted, f"{target.__name__} defaults {defaulted}")


def test_r7_18_the_comparison_is_returned_not_stored_and_is_not_a_certificate() -> None:
    """``UncertifiedComparison`` is a frozen, ``@final`` dataclass, NOT a
    ``ContentIdentity`` and not a ``VerificationCertificate`` (it certifies
    nothing); its fields are the four arguments and the two readings (the
    readings ``init=False``, read at construction); ``agrees`` is a derived
    property, not a field, and it has no floor; two calls return two
    objects."""
    from orpheus.numerics.content import ContentIdentity

    cls: type = s5.c(s7.VERIFICATION, "UncertifiedComparison")
    params = getattr(cls, "__dataclass_params__", None)
    require(params is not None and params.frozen, "UncertifiedComparison is not a frozen dataclass")
    require(getattr(cls, "__final__", False) is True, "UncertifiedComparison is not @final")
    require(not issubclass(cls, ContentIdentity), "UncertifiedComparison is a content value")
    require(not issubclass(cls, s5.c(s7.VERIFICATION, "VerificationCertificate")), "an uncertified comparison is a certificate")
    fields = {f.name: f.init for f in dataclasses.fields(cls)}
    require(fields == {"answer": True, "observable": True, "reference": True, "tolerance": True,
                       "reading": False, "reference_reading": False}, f"fields {fields}")
    require(isinstance(getattr(cls, "agrees", None), property), "agrees is not a property")
    require(not hasattr(cls, "floor_holds"), "an uncertified comparison has a floor")
    a, b = _compare(1.0), _compare(1.0)
    require(a is not b, "the verb returned a stored object")


# ═════════════════════════════════════════════════════════════════════════════
# R7.19-R7.20 — the ratio rule, spelled once (#405 P2 step 7b.2 prerequisite)
# ═════════════════════════════════════════════════════════════════════════════


def test_r7_19_every_answer_reads_a_ratio_through_the_one_quotient(monkeypatch: pytest.MonkeyPatch) -> None:
    """ROUTE (the one-definition witness, X4): ``Ratio.quotient(read)`` is the
    ratio rule for every answer. Rebinding it to a decoy that divides the
    other way moves the readings of BOTH the production homogeneous answer and
    the exact reference at once, each calling it once, so neither keeps a
    private rule. Activation: honestly, both read the group-0/group-1 ratio
    (not its reciprocal), and it is not 1."""
    mixture = sf.fuel()
    result = s7.solve_homogeneous(mixture)
    reference = s6.exact_medium_reference(_medium(mixture))
    Ratio = s5.c(s5.OBSERVABLE, "Ratio")
    ratio = Ratio(s6.group_flux(0), s6.group_flux(1))
    honest = (result.read(ratio).value, reference.read(ratio).value)
    require(abs(honest[0] - float(result.flux[0]) / float(result.flux[1])) <= 4 * math.ulp(honest[0]) and honest[0] != 1.0,
            f"activation: the production ratio {honest[0]!r}")
    calls: list[str] = []

    def reversed_quotient(self: Any, read: Any) -> Any:
        calls.append(type(read.__self__).__name__)
        return read(self.denominator) / read(self.numerator)

    monkeypatch.setattr(Ratio, "quotient", reversed_quotient)
    moved = (result.read(ratio).value, reference.read(ratio).value)
    require(sorted(calls) == ["HomogeneousResult", "ReferenceSolution"], f"the quotient was called by {calls}")
    for h, m, who in zip(honest, moved, ("production", "reference")):
        require(abs(m * h - 1.0) <= 1e-14, f"the {who} reading did not move with the rule: {h!r} -> {m!r}")


_MEASURED_QUOTIENTS = (
    ("dyadic", lambda: Measured(3.0) / Measured(2.0), 1.5),
    ("non-dyadic", lambda: Measured(1.0) / Measured(3.0), 1.0 / 3.0),
    ("negative", lambda: Measured(-3.0) / Measured(2.0), -1.5),
)


@pytest.mark.parametrize("quotient, value", [r[1:] for r in _MEASURED_QUOTIENTS], ids=[r[0] for r in _MEASURED_QUOTIENTS])
def test_r7_20_the_quotient_of_two_measurements_is_a_measurement(quotient: Any, value: float) -> None:
    """A ratio observable's production reading: ``Measured / Measured`` is
    ``Measured(a.value / b.value)``, the float quotient."""
    q = quotient()
    require(type(q) is Measured and q.value == value, f"{q!r}")


def _untyped(value: Any) -> Any:
    """A value the type checker does not narrow (the row tests a RUNTIME refusal)."""
    return value


_NOT_MEASURED_OPERANDS = (
    ("over-an-enclosure", lambda: Measured(1.0) / s5.enclosure(2.0, 0.0)),
    ("an-enclosure-over", lambda: s5.enclosure(2.0, 0.0) / Measured(1.0)),
    ("over-an-uncertified", lambda: Measured(1.0) / s6.uncertified(2.0)),
    ("over-a-bare-float", lambda: Measured(1.0) / _untyped(2.0)),
)


@pytest.mark.parametrize("quotient", [r[1] for r in _NOT_MEASURED_OPERANDS], ids=[r[0] for r in _NOT_MEASURED_OPERANDS])
def test_r7_20_a_measurement_does_not_divide_by_a_reference_reading(quotient: Any) -> None:
    """A production reading and a reference reading never combine (G3, the
    disjoint sums): ``Measured / Enclosure`` and the reverse, ``Measured /
    Uncertified`` and ``Measured / float`` are a ``TypeError``. Division by a
    zero measurement raises ``ZeroDivisionError``."""
    with pytest.raises(TypeError):
        quotient()
    with pytest.raises(ZeroDivisionError):
        returned = Measured(1.0) / Measured(0.0)
        require(False, f"a division by a zero measurement returned {returned!r}")
