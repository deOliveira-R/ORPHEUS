r"""The first production reading and the first verification, on the homogeneous answer (#405 P2 step 7a, R7.8, R7.9; step 7b.2, R7.19).

Specified by the test-architect (2026-10-03, ``.claude/plans/reference_p2_spec.md``
§1.7a and §1.7b.2). These rows lived in ``tests/gates/reference/test_verification.py``
until the layer-DAG re-filing (Phase 1b, 2026-10-26, ``.claude/plans/layer_dag.md``):
their claims are about the homogeneous solver's answer, which the reference layer
cannot name.

* R7.8: ``HomogeneousResult.read``, the first production reading, reads its own
  answer (k, group fluxes, weighted sums, ratios, a constant ``Symbolic``) and
  refuses what a zero-dimensional answer cannot read;
* R7.9: the production 0-D solve agrees with the exact infinite medium under
  ``verify_agreement``, and a perturbed production disagrees (the X1 negative);
* R7.19: ``Ratio.quotient`` is the ratio rule of every answer: the homogeneous
  answer and a reference solution both read a ratio through it.
"""

from __future__ import annotations

import dataclasses
import math
from typing import Any

import numpy as np
import pytest

from orpheus.derivations.common import exact_homogeneous
from orpheus.homogeneous import solver
from orpheus.numerics.outcome import Measured
from tests.gates._content_identity_helpers import require
from tests.gates.reference import _step5 as s5
from tests.gates.reference import _step6 as s6
from tests.gates.reference import _step7 as s7
from tests.gates.specification import _fixtures as sf

pytestmark = pytest.mark.foundation


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

    result = solver.solve_homogeneous_infinite(sf.fuel())
    k = result.read(s5.eigenvalue())
    require(type(k) is Measured and k.value == result.k_inf, f"{k!r}")
    for g in range(2):
        phi = result.read(s6.group_flux(g))
        require(type(phi) is Measured and phi.value == float(result.flux[g]), f"group {g}: {phi!r}")
    w = (0.25, 3.0)
    weighted = result.read(_flux(w)).value
    require(abs(weighted - (w[0] * result.flux[0] + w[1] * result.flux[1])) <= 4 * math.ulp(weighted), f"{weighted!r}")
    ratio = result.read(s5.OBSERVABLE.Ratio(s6.group_flux(0), s6.group_flux(1))).value
    require(ratio == float(result.flux[0]) / float(result.flux[1]), f"ratio {ratio!r}")
    constant = s5.OBSERVABLE.FluxIntegral(Symbolic.of(0.25, 3))
    require(abs(result.read(constant).value - weighted) <= 4 * math.ulp(weighted), "a constant Symbolic weight reads differently from its table")


_REFUSED_READS = (
    ("point-value", lambda: s5.point(0.5, 0), "position"),
    ("three-groups", lambda: _flux((1.0, 1.0, 1.0)), "group"),
    ("two-regions", lambda: s6.group_flux(0, regions=2), "region"),
)


@pytest.mark.parametrize("make, fragment", [r[1:] for r in _REFUSED_READS], ids=[r[0] for r in _REFUSED_READS])
def test_r7_8_what_a_zero_dimensional_answer_cannot_read_is_refused(make: Any, fragment: str) -> None:
    result = solver.solve_homogeneous_infinite(sf.fuel())
    with pytest.raises(ValueError, match=fragment):
        result.read(make())


# ═════════════════════════════════════════════════════════════════════════════
# R7.9 — the first verification: production homogeneous against the exact medium
# ═════════════════════════════════════════════════════════════════════════════

def _library_a_4g() -> Any:
    from orpheus.data.macro_xs.synthetic import get_mixture

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


@pytest.mark.rests_on(
    "tests/gates/homogeneous/test_kinf_exact_reference.py::test_kinf_is_the_exact_eigenvalue_within_the_derived_bound",
    "tests/gates/homogeneous/test_kinf_exact_reference.py::test_flux_is_the_exact_gauged_eigenvector_within_the_derived_bound",
    "tests/gates/derivations/test_exact_medium_reference.py::test_r6_6_the_exact_infinite_medium_reads_the_exact_answer",
)
@pytest.mark.parametrize("make", [m for _, m in _MIXTURES], ids=[n for n, _ in _MIXTURES])
def test_r7_9_the_homogeneous_solver_agrees_with_the_exact_medium(make: Any) -> None:
    """REFERENCE (exact rationals) at L1: the production 0-D solve of each
    mixture agrees with the exact infinite medium on k, on each group's flux
    and on the first two groups' ratio, at 1e-13 relative, the algebraic error
    ``NotApplicable`` (a direct solve). This is the step's own positive leg;
    the next row is the instrument's negative one."""
    mixture = make()
    reference = exact_homogeneous.exact_infinite_medium_reference(_medium(mixture))
    result = solver.solve_homogeneous_infinite(mixture)
    groups = len(result.flux)
    ratio = s5.OBSERVABLE.Ratio(s6.group_flux(0, groups=groups), s6.group_flux(1, groups=groups))
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
    reference = exact_homogeneous.exact_infinite_medium_reference(_medium(mixture))
    cert = s7.verify_agreement(solver.solve_homogeneous_infinite(perturbed), s5.eigenvalue(), reference,
                               _REL_TOL * abs(reference.read(s5.eigenvalue()).value), s7.DIRECT)
    require(cert.floor_holds and not cert.agrees, f"{cert!r}")


# ═════════════════════════════════════════════════════════════════════════════
# R7.19 — the ratio rule, spelled once (#405 P2 step 7b.2 prerequisite)
# ═════════════════════════════════════════════════════════════════════════════


def test_r7_19_every_answer_reads_a_ratio_through_the_one_quotient(monkeypatch: pytest.MonkeyPatch) -> None:
    """ROUTE (the one-definition witness, X4): ``Ratio.quotient(read)`` is the
    ratio rule for every answer. Rebinding it to a decoy that divides the
    other way moves the readings of BOTH the production homogeneous answer and
    a reference solution at once, each calling it once, so neither keeps a
    private rule. Activation: honestly, both read the group-0/group-1 ratio
    (not its reciprocal), and it is not 1. The reference is the test-side
    table derivation's (``_step6.medium_table``, group fluxes 30/7 and 11/13):
    the rule is the reader's, so any derivation witnesses it (re-based
    2026-10-26 off the exact medium, a derivation the claim does not need)."""
    mixture = sf.fuel()
    result = solver.solve_homogeneous_infinite(mixture)
    reference = s6.reference(*s6.medium_table())
    Ratio = s5.OBSERVABLE.Ratio
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
