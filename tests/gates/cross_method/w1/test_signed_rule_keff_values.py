r"""G04w (d): the LOCAL value leg of the signed-rule gate (N04, spec 7.5).

keff of the vacuum-faced 2 by 2 box on ``lebedev(13)`` equals 0.29309086 (scattering order 0) and 0.27773609 (order 1), and on ``lebedev(11)`` 0.29347604 and 0.27820053, within the solve's own convergence bound read
from its exit report (``Asserted(bound=1e-06, ...)``), not a chosen tolerance. The values were recomputed on the base tip by this gate's probe (``[M]`` 0.29347604426414825, 0.2782005320318106,
0.29309086444039906, 0.27773609108832553). A toleranced comparison of a computation with its recorded value is host-keyed, so the file is not in the CI step. Green at the base tip; the mutation of
``test_signed_rule_is_a_legal_factor.py`` (reinstall the negative-weight refusal in ``Factor._require_constructible``) reddens the two ``lebedev(13)`` rows here by the refusal.
"""

from __future__ import annotations

import pytest

from tests.gates.cross_method.w1.test_signed_rule_is_a_legal_factor import _require, _vacuum_box_solve

pytestmark = pytest.mark.foundation

# the keff of the vacuum-faced box: (quadrature order, scattering order) -> the opening's printed value
LOCAL_VALUES = {(13, 0): 0.29309086, (13, 1): 0.27773609, (11, 0): 0.29347604, (11, 1): 0.27820053}


@pytest.mark.parametrize("key", sorted(LOCAL_VALUES))
def test_local_keff_values_within_the_solves_own_bound(key: tuple[int, int]) -> None:
    solution = _vacuum_box_solve(*key)
    bound = float(solution.exit_report.balance.bound)
    _require(bound == 1e-06, f"the exit report's asserted bound reads {bound}, not 1e-06")
    keff = float(solution.outcome.keff)
    _require(abs(keff - LOCAL_VALUES[key]) <= bound, f"keff {keff!r} on lebedev({key[0]}), scattering order {key[1]}, differs from {LOCAL_VALUES[key]} by more than the solve's bound {bound}")
