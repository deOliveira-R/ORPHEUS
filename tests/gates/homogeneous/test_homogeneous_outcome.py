r"""The infinite medium's answer carries its outcome — step 3, U1 (2026-09-17).

``HomogeneousResult.outcome`` is an ``EigenOutcome``: the k-eigen question over
the 0-D pencil, the gauged representative as the posed ``(ng, 1)`` column, λ = k∞
with its one-point trajectory, and the ``ScaleGauge`` that fixed the representative
(νΣf·φ = 100 n/cm³/s — the one shipped deliberately-named target, recorded as the
functional that ran rather than applied anonymously). ``k_inf`` and ``flux`` are
READ off it; the byte-stability gate pins that the arithmetic did not move.

⚠ ≥ 2G on purpose: ``[M]`` ``IntegratedReactionRate.evaluate`` silently accepts an
``(ng,)`` vector and returns the WRONG number (200.0 at 2g where the column reads
100.0); at 1g the two coincide and a law written there would be blind.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.common.xs_library import get_mixture
from orpheus.homogeneous.solver import solve_homogeneous_infinite
from orpheus.numerics.outcome import EigenOutcome
from orpheus.numerics.posing import K_MAP

pytestmark = pytest.mark.foundation


def _require(condition: object, message: str) -> None:
    if not condition:
        raise AssertionError(message)


@pytest.mark.parametrize("ng", ["2g", "4g"])
def test_law_the_outcome_reproduces_k_inf_and_the_gauge_target(ng: str) -> None:
    result = solve_homogeneous_infinite(get_mixture("A", ng))
    out = result.outcome
    _require(isinstance(out, EigenOutcome) and out.posing.spectral_map is K_MAP, "the kind is the eigen one, under the k map")
    _require(out.lam == result.k_inf and out.keff == result.k_inf, "k_inf IS the outcome's λ (one storage)")
    _require(out.trajectory == (result.k_inf,), "a direct solve has a one-point trajectory")
    _require(out.state.shape == (result.flux.size, 1), "the state is the posed (ng, 1) column")
    _require(np.shares_memory(result.flux, out.state), "flux is a VIEW of the state, not a copy")
    _require(abs(out.gauge.functional(out.state) - 100.0) <= 1e-14 * 100.0, "the recorded section lands on its target: νΣf·φ = 100 n/cm³/s")
    _require(out.gauge.target == 100.0, "…and says so")
    _require(abs(out.rayleigh() - out.lam) <= 8 * np.finfo(float).eps * abs(out.lam), "the posing's own quotient reproduces k∞ (``[M]`` 2.2e-16)")


def test_negative_a_moved_target_is_seen(ng: str = "2g") -> None:
    result = solve_homogeneous_infinite(get_mixture("A", ng))
    out = result.outcome
    _require(abs(out.gauge.functional(2.0 * out.state) - 100.0) > 50.0, "a rescaled state no longer sits on the section — the law has teeth")
