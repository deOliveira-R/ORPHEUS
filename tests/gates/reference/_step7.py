"""The step-7a gates' adapter and fixtures (#405 P2, ``.claude/plans/reference_p2_spec.md`` §1.7a).

Every production name the step-7a gates reach is spelled ONCE here, resolved on
its module at call time (a rebinding battery arm reaches every row, lessons
``L100``); a rename is one edit in this file.
"""

from __future__ import annotations

from fractions import Fraction
from typing import Any

from orpheus.numerics.outcome import Asserted, Measured, NotApplicable, NotYet
from tests.gates.reference import _step5 as s5
from tests.gates.reference import _step6 as s6

VERIFICATION = "orpheus.reference.verification"
HOMOGENEOUS = "orpheus.homogeneous.solver"

NOT_YET = NotYet(564, "no stage-A error estimator of a functional yet")
DIRECT = NotApplicable("a direct solve: no iteration, no algebraic error beyond rounding")


def verify_agreement(answer: Any, observable: Any, reference: Any, tolerance: Any, algebraic_error: Any) -> Any:
    return s5.c(VERIFICATION, "verify_agreement")(answer, observable, reference, tolerance, algebraic_error)


def verify_order(answers: Any, observable: Any, reference: Any, tolerance: Any, algebraic_errors: Any, order: Any, band: Any) -> Any:
    return s5.c(VERIFICATION, "verify_order")(answers, observable, reference, tolerance, algebraic_errors, order, band)


def not_valid() -> type[BaseException]:
    return s5.c(VERIFICATION, "ReferenceNotValid")


def unestablished() -> type[BaseException]:
    return s5.c(VERIFICATION, "Unestablished")


def answer_protocol() -> Any:
    return s5.c(VERIFICATION, "ProductionAnswer")


def solve_homogeneous(mixture: Any) -> Any:
    return s5.c(HOMOGENEOUS, "solve_homogeneous_infinite")(mixture)


class FixedAnswer:
    """A production answer double: reads a table of observables, counting its reads."""

    def __init__(self, table: dict[Any, Any]) -> None:
        self.table = dict(table)
        self.reads: list[Any] = []

    def read(self, observable: Any) -> Any:
        self.reads.append(observable)
        value = self.table[observable]
        return Measured(value) if isinstance(value, float) else value


def reference_of(k: Fraction, certificate: Any = "valid") -> Any:
    """The infinite medium with an exact derivation of k (and nothing else), and the asked certificate."""
    spec = s5.eigen_medium()
    derivation = s6.TableDerivation({s5.eigenvalue(): k})
    if certificate == "valid":
        cert = s5.certificate({s5.eigenvalue(): s6.exact_claim(k)})
    elif certificate == "invalid":
        cert = s5.certificate({s5.eigenvalue(): s5.claim(1e-300, s6.exact_of(k))})
    elif certificate == "withdrawn":
        cert = s5.certificate({s5.eigenvalue(): s6.exact_claim(k)}, standing=s5.withdrawal(issue=4321))
    else:
        cert = None
    return s6.reference(spec, derivation, cert)


def reference_with_bound(value: float, bound: float) -> Any:
    """The infinite medium whose derivation establishes k as ``value +- bound`` (a derived bound), certified."""
    spec = s5.eigen_medium()
    derived = s5.derived(value, bound, "a test bound")
    derivation = s6.TableDerivation({s5.eigenvalue(): derived})
    cert = s5.certificate({s5.eigenvalue(): s5.claim(1.0, derived)})
    return s6.reference(spec, derivation, cert)


__all__ = ["Asserted", "Measured"]
