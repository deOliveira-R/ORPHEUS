"""The step-6 gates' adapter and fixtures (#405 P2, ``.claude/plans/reference_p2_spec.md`` §1.6 and §1.7b.1).

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


def uncertified(value: Any) -> Any:
    """``Uncertified(value)``: a reference's value with no derived bound (step 7b.1)."""
    return s5.c(s5.READING, "Uncertified")(value)


def uncertified_class() -> type:
    return s5.c(s5.READING, "Uncertified")


def evaluation() -> Any:
    """The ``Evaluation`` sum a derivation returns: ``Exact | DerivedBound | Uncertified``."""
    return s5.c(SOLUTION, "Evaluation")


def derivation_protocol() -> Any:
    return s5.c(SOLUTION, "Derivation")


def exact_medium_reference(specification: Any) -> Any:
    return s5.c(EXACT_MEDIUM, "exact_infinite_medium_reference")(specification)


def exact_medium_of(mixture: Any) -> Any:
    return s5.c(EXACT_MEDIUM, "exact_infinite_medium_of")(mixture)


# ── a test-side derivation ──────────────────────────────────────────────────


class TableDerivation:
    """A derivation whose natural extension is a table keyed by observable, counting its calls.

    An entry is an exact ``Fraction`` (evaluated as ``Exact``), a ``float``
    (evaluated as ``Uncertified``: the family-without-a-bound case, step
    7b.1), or an evaluation itself (an ``Exact``, a ``DerivedBound``, an
    ``Uncertified``, or a deliberately wrong object for an admission row). An
    observable with no entry is a fixture error (``KeyError``), never a
    production case.
    """

    def __init__(self, table: dict[Any, Any]) -> None:
        self.table = dict(table)
        self.calls: list[Any] = []

    def evaluate(self, observable: Any) -> Any:
        self.calls.append(observable)
        entry = self.table[observable]
        if isinstance(entry, Fraction):
            return exact_of(entry)
        if isinstance(entry, float):
            return uncertified(entry)
        return entry


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
