"""The homogeneous population precondition: exactly the eight shipped cases, each producing."""

from __future__ import annotations

import pytest

from ._homogeneous_population import EXPECTED_CASES, mixture_cases

pytestmark = pytest.mark.foundation


def test_the_homogeneous_population_is_the_eight_producing_cases() -> None:
    """The population precondition every consumer of ``_mixture_cases`` rests
    on: exactly the eight shipped cases, each PRODUCING (an eigenvalue entry on
    a non-producing mixture is ``k = 0`` and a dead row)."""
    cases = mixture_cases()
    if tuple(sorted(cases)) != EXPECTED_CASES:
        pytest.fail(f"population changed: {sorted(cases)}")
    barren = [name for name, mix in cases.items() if not mix.is_producing]
    if barren:
        pytest.fail(f"non-producing mixtures in the eigenvalue population: {barren}")
