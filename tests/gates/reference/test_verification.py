r"""The verification certificate and the comparison verbs (#405 P2 step 7a, R7.1-R7.10).

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
    orders are returned, one per consecutive pair, each 2 to rounding; a
    first-order ladder does not hold against order 2."""
    second = _order(_ladder(2.0))
    orders = tuple(second.observed_orders)
    require(len(orders) == 2 and all(abs(o - 2.0) < 1e-9 for o in orders), f"observed orders {orders}")
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

_MIXTURES = (("fuel-2g", sf.fuel), ("fission-only-2g", sf.fission_only), ("fuel-4g", lambda: sf.fuel(4)))
#: The tolerance, relative: ``[M]`` 2026-10-03 the production k and group fluxes
#: sit within 1.65e-16 relative of the exact ones over these three mixtures
#: (600x margin), and the exact reference's bound (half an ulp, ~1.1e-16
#: relative) meets the floor tol/10 with 90x margin.
_REL_TOL = 1e-13


def _medium(mixture: Any) -> Any:
    from orpheus.data.cells import CellCoefficient, Channel
    from orpheus.numerics.question import Eigen
    from orpheus.specification import InfiniteMediumSpecification

    return InfiniteMediumSpecification(0, mixture, Eigen(CellCoefficient.every(Channel.FISSION_EMISSION)))


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
