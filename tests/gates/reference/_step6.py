"""The step-6 gates' adapter and fixtures (#405 P2, ``.claude/plans/reference_p2_spec.md`` §1.6).

Every production name the step-6 gates reach is spelled ONCE here, resolved on
its module at call time (a rebinding battery arm reaches every row, lessons
``L100``); a rename is one edit in this file. The step-5 adapter supplies the
observables, readings and certificate builders.
"""

from __future__ import annotations

from fractions import Fraction
from typing import Any

import numpy as np

from orpheus.numerics.mesh_free_function import RegionwiseConstant
from tests.gates.reference import _step5 as s5

SOLUTION = "orpheus.reference.solution"
EXACT_MEDIUM = "orpheus.derivations.common.exact_homogeneous"


def reference(specification: Any, derivation: Any, certificate: Any) -> Any:
    return s5.c(SOLUTION, "ReferenceSolution")(specification, derivation, certificate)


def not_certified() -> type[BaseException]:
    return s5.c(SOLUTION, "NotCertified")


def derivation_protocol() -> Any:
    return s5.c(SOLUTION, "Derivation")


def exact_medium_reference(specification: Any) -> Any:
    return s5.c(EXACT_MEDIUM, "exact_infinite_medium_reference")(specification)


def exact_medium_of(mixture: Any) -> Any:
    return s5.c(EXACT_MEDIUM, "exact_infinite_medium_of")(mixture)


# ── a test-side derivation ──────────────────────────────────────────────────


class TableDerivation:
    """A derivation whose natural extension is a table keyed by observable, counting its calls.

    An entry is an exact ``Fraction`` (established as ``Exact``) or an
    establishment itself; an observable with no entry has no derivable bound,
    so ``establish`` raises ``NotCertified`` (the family-without-a-bound case).
    """

    def __init__(self, table: dict[Any, Any]) -> None:
        self.table = dict(table)
        self.calls: list[Any] = []

    def establish(self, observable: Any) -> Any:
        self.calls.append(observable)
        if observable not in self.table:
            raise not_certified()(f"TableDerivation: no derived bound for {observable!r}, so it cannot be certified")
        entry = self.table[observable]
        return exact_of(entry) if isinstance(entry, Fraction) else entry


def exact_of(value: Fraction) -> Any:
    import sympy

    return s5.exact(sympy.Rational(value.numerator, value.denominator))


def group_flux(group: int, groups: int = 2, regions: int = 1) -> Any:
    """The flux integral whose weight is the indicator of one group (on ``regions`` regions)."""
    w = np.zeros((regions, groups))
    w[:, group] = 1.0
    return s5.c(s5.OBSERVABLE, "FluxIntegral")(RegionwiseConstant(w))


def exact_claim(value: Fraction, target: float = 1e-12) -> Any:
    return s5.claim(target, exact_of(value))


def medium_table() -> tuple[Any, TableDerivation, Any]:
    """The eigen infinite medium with a table derivation and an exact certificate on k and both group fluxes."""
    k, f0, f1 = Fraction(7, 5), Fraction(30, 7), Fraction(11, 13)
    spec = s5.eigen_medium()
    table = {s5.eigenvalue(): k, group_flux(0): f0, group_flux(1): f1}
    cert = s5.certificate({s5.eigenvalue(): exact_claim(k), group_flux(0): exact_claim(f0), group_flux(1): exact_claim(f1)})
    return spec, TableDerivation(table), cert
