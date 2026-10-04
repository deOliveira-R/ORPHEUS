r"""``Symbolic.without`` preserves the function's values (#405 P2 step 7b.2.1, gate R7b2.10).

``without(*coordinates)`` is the constructive face of "independent of these
coordinates": it may only drop a coordinate the function does not depend on,
and it must return THE SAME FUNCTION. Until the review round of 7b.2.1 it
called ``sympy.simplify`` first, which rewrote ``Piecewise((1, sin(3r) > 0),
(0, True))`` as its first period alone (``[M]`` qa: the value at r = 2.5 went
from 1 to 0), and no gate pinned value preservation (85 of 85 passed across
the fix). Each row lifts the function before and after at sample points.

First red: the ``simplify`` restored (battery arm ``without_simplifies``,
``scratch/reference_architecture/p2/ta/step7b2/read/battery/``). Row 1 is the
catcher of ERR-096 (re-dropped in process under ``-O``: it reddens, row 2 stays
green).
"""
from __future__ import annotations

import pytest
import sympy

from orpheus.numerics.mesh_free_function import Symbolic

pytestmark = pytest.mark.foundation

_R, _MU, _PHI = Symbolic.r, Symbolic.mu, Symbolic.phi
_SAMPLES_R = (0.3, 1.1, 2.5, 4.0, 7.7)
_SAMPLES_PHI = (0.0, 0.9, 2.3, 5.1)


def _lift(expression: sympy.Expr, point: dict) -> float:
    return float(expression.subs(point, simultaneous=True))


@pytest.mark.catches("ERR-096")
def test_r7b2_10_without_preserves_a_periodic_step() -> None:
    """The sin step keeps its value at every sample (1 at r = 2.5, inside its SECOND period: the activation)."""
    step = sympy.Piecewise((1, sympy.sin(3 * _R) > 0), (0, True))
    before = Symbolic.of(step)
    after = before.without(_MU, _PHI)
    assert _lift(step, {_R: 2.5}) == 1.0
    for r in _SAMPLES_R:
        assert _lift(after.expressions[0], {_R: r}) == _lift(before.expressions[0], {_R: r}), r


def test_r7b2_10_without_removes_a_cancelling_direction_dependence_and_keeps_the_values() -> None:
    """``(sin²φ + cos²φ)·(r + 1)`` does not depend on φ: ``without(phi)`` admits it, and the result equals the
    original at every (r, φ) sample (to 4 ulp: the substitution is exact, the float lift rounds)."""
    expression = (sympy.sin(_PHI) ** 2 + sympy.cos(_PHI) ** 2) * (_R + 1)
    before = Symbolic.of(expression, 2 * _R)
    after = before.without(_PHI)
    for g in range(2):
        for r in _SAMPLES_R:
            for phi in _SAMPLES_PHI:
                a = _lift(before.expressions[g], {_R: r, _PHI: phi})
                b = _lift(after.expressions[g], {_R: r, _PHI: phi})
                assert abs(a - b) <= 4 * abs(a) * 2.0**-52, (g, r, phi, a, b)
