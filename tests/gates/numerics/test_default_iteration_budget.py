"""The budget is DERIVED from the tolerance, and it is one law (#340, N3).

Five hardcoded ``max_inner`` constants shipped in SN plus a sixth in
:class:`~orpheus.numerics.iteration.KEigenvalue`, in three
``(budget, tolerance)`` combinations.  A constant cannot track a tolerance it
does not know about, and `[M]` all of them were SHORT on the very
configuration that opened #340.

These gates pin the two properties that make the replacement defensible:

* the budget **moves with the tolerance by the geometric law**, not by an
  arbitrary factor — the shipped pair differed by 5x where the law says 1.5;
* the prior (:func:`default_iteration_budget`, no data) and the posterior
  (:meth:`StoppingCriterion.projected_iterations`, fitted from a trajectory)
  are **one statement** evaluated at different intercepts, so they cannot
  drift apart into two spellings of the budget law.

Runtime note (vv Mode 8): collected module, so bare ``assert`` survives ``-O``.
"""

from __future__ import annotations

import math

import pytest

from orpheus.numerics.convergence import (
    StoppingCriterion,
    default_iteration_budget,
    resolve_iteration_budget,
)

#: `[M]` 2026-08-09, real converged solves of the d=3 all-reflective pure
#: absorber (LS-S4, 2g, 3x4x5 cells) — the configuration whose truncated exit
#: opened #340.  `needed` is OBSERVED (``len(history.flux_residuals)`` on a
#: run with `converged is True`), not fitted.  Reproduce by solving at
#: `max_inner=6000` and reading the length.
_MEASURED_D3_ABSORBER = [
    (1e-9, 1007),
    (1e-12, 1473),
    (1e-13, 1631),
    (1e-15, 2031),
]


@pytest.mark.foundation
class TestTheBudgetTracksItsTolerance:
    """A constant cannot; the law must."""

    @pytest.mark.parametrize("tolerance", [1e-6, 1e-8, 1e-10, 1e-12, 1e-15])
    def test_a_tighter_tolerance_costs_strictly_more(
        self, tolerance: float
    ) -> None:
        assert default_iteration_budget(tolerance) > default_iteration_budget(
            tolerance * 100.0
        )

    def test_the_ratio_between_two_tolerances_IS_the_log_ratio(self) -> None:
        r"""``n(tol) = |ln tol| / |ln rho|``, so the ratio of two budgets is
        the ratio of the log-tolerances — free of ``rho`` entirely.

        `[M]` the tree shipped ``1000/200 = 5.0`` for a pair whose law-ratio
        is ``ln(1e-12)/ln(1e-8) = 1.5``.  That factor-of-3.3 discrepancy is
        the whole argument for deriving rather than choosing.
        """
        ratio = default_iteration_budget(1e-12) / default_iteration_budget(1e-8)
        assert ratio == pytest.approx(math.log(1e-12) / math.log(1e-8), rel=2e-3)

    @pytest.mark.parametrize("tolerance, needed", _MEASURED_D3_ABSORBER)
    def test_it_covers_the_configuration_that_opened_the_issue(
        self, tolerance: float, needed: int
    ) -> None:
        """⭐ The acceptance criterion, against MEASURED counts.

        `[M]` the shipped ``1000`` covers **none** of these four rows; the
        derived budget covers **all** of them.  The first row is the sharp
        one — 1007 needed against 1000 shipped, a shortfall of **seven
        sweeps**, which is the margin the founding defect's gate rode for
        months before a correct quadrature change (#337) moved it out.
        """
        assert needed > 1000, "row no longer discriminates against the old cap"
        assert default_iteration_budget(tolerance) >= needed

    @pytest.mark.parametrize("tolerance", [0.0, -1e-8, 1.0, 2.0])
    def test_a_tolerance_outside_the_unit_interval_is_refused(
        self, tolerance: float
    ) -> None:
        """The law is sized against a unit initial residual, so a tolerance
        at or above 1 is already met before the first iteration and a budget
        for it is meaningless."""
        with pytest.raises(ValueError, match=r"\(0, 1\)"):
            default_iteration_budget(tolerance)

    @pytest.mark.parametrize("served_rate", [0.0, 1.0, 1.5, -0.5])
    def test_a_NON_CONTRACTING_served_rate_is_refused(
        self, served_rate: float
    ) -> None:
        """At ``rho >= 1`` there is no finite budget to derive, so promising
        one would be the same lie ``projected_iterations`` refuses to tell."""
        with pytest.raises(ValueError, match="contract"):
            default_iteration_budget(1e-8, served_rate=served_rate)

    def test_a_slower_served_rate_buys_a_bigger_budget(self) -> None:
        assert default_iteration_budget(
            1e-8, served_rate=0.995
        ) > default_iteration_budget(1e-8, served_rate=0.9)


@pytest.mark.foundation
class TestThePriorAndThePosteriorAreOneLaw:
    """⭐ The single-source gate.

    ``default_iteration_budget`` evaluates the geometric law with no data;
    ``StoppingCriterion.projected_iterations`` fits ``(intercept, rate)`` from
    a trajectory and evaluates the SAME law.  Feeding the projection a decay
    that matches the prior's own assumptions must therefore reproduce the
    prior EXACTLY — off-by-one included.  A second spelling of the law would
    show up here as a one-iteration disagreement, which is precisely the kind
    of drift that goes unnoticed between two hand-written formulas.
    """

    @pytest.mark.parametrize("tolerance", [1e-6, 1e-8, 1e-10, 1e-12])
    @pytest.mark.parametrize("served_rate", [0.9, 0.986, 0.995])
    def test_a_unit_decay_at_the_served_rate_reproduces_the_prior(
        self, tolerance: float, served_rate: float
    ) -> None:
        prior = default_iteration_budget(tolerance, served_rate=served_rate)
        posterior = StoppingCriterion(
            name="residual",
            # r0 = 1 at exactly the served rate: the prior's own assumptions,
            # made into data.  Long enough that the tail fit is clean.
            trajectory=tuple(served_rate**k for k in range(400)),
            tolerance=tolerance,
        ).projected_iterations()
        assert posterior == prior


@pytest.mark.foundation
class TestNoneMeansDeriveAndAnIntIsDeliberate:
    """The resolution verb, spelled once for six call sites."""

    def test_none_derives(self) -> None:
        assert resolve_iteration_budget(None, 1e-8) == default_iteration_budget(
            1e-8
        )

    def test_an_explicit_budget_is_never_second_guessed(self) -> None:
        """Starving a solve on purpose is correct API use — `[M]` #340's audit
        found 7 deliberate truncations in the suite, every one of them passing
        an explicit budget.  A resolution that "helpfully" raised a small
        explicit cap would silently retune all seven."""
        assert resolve_iteration_budget(3, 1e-12) == 3

    def test_an_explicit_budget_ignores_the_tolerance_entirely(self) -> None:
        assert resolve_iteration_budget(50, 1e-6) == resolve_iteration_budget(
            50, 1e-15
        )
