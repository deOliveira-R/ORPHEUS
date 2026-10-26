"""The quotient of two measurements is a measurement (#405 P2 step 7b.2, R7.20).

Specified by the test-architect (2026-10-03, ``.claude/plans/reference_p2_spec.md``
§1.7b.2). The row lived in ``tests/gates/reference/test_verification.py`` until
the layer-DAG re-filing (Phase 1b, 2026-10-26, ``.claude/plans/layer_dag.md``):
it states ``Measured.__truediv__``, a numerics law. Its sibling, a measurement
never divides by a reference reading, stays in the reference tests (it names
``Uncertified``, a reference type).
"""

from __future__ import annotations

from typing import Any

import pytest

from orpheus.numerics.outcome import Measured
from tests.gates._content_identity_helpers import require

pytestmark = pytest.mark.foundation


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
