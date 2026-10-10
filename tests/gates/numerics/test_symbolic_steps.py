"""``Symbolic.steps``: the one definition of where a symbolic weight is not smooth in r, and the integrator that reads it.

Re-homed in P1 step (e1b) of the characteristic-reference campaign
(``.claude/plans/characteristic_reference_architecture.md``, "Step (e)") from the
old trajectory-resolvent family's ``test_trajectory_resolvent_reference.py``
(rows R7b2.6 and R7b2.2.3), which step (e2) deletes with the family. The rows
never read the family: they pin :meth:`orpheus.numerics.mesh_free_function.Symbolic.steps`
and its exact consumer :meth:`orpheus.mesh.structured.Mesh1D.cell_integrals`, which
outlive it; without them both lose their only direct gates (``[M]`` 2026-10-10,
``git grep -P "\\bsteps\\("`` over ``tests/``: 8 lines in the old file and 3 in its
API helper, nothing else). The bodies are the old rows', with the helper
``_trajectory_resolvent_api.steps`` written out (sorted floats of
``weight.steps((low, high))``) and its refusal pattern copied.
"""
from __future__ import annotations

import math

import numpy as np
import pytest
import sympy

from orpheus.geometry import CoordSystem
from orpheus.numerics.mesh_free_function import Symbolic

pytestmark = pytest.mark.foundation

#: The old API helper's ``UNLOCATED_REFUSAL``: production's scope-edge message for an unlocated non-smooth construct.
_UNLOCATED_REFUSAL = r"(?i)polynomial in r|not located|allow|scope boundary"


def _steps(weight: Symbolic, low: float, high: float) -> tuple[float, ...]:
    """The step locations of ``weight`` inside ``[low, high]``, sorted."""
    return tuple(sorted(float(x) for x in weight.steps((low, high))))


def test_the_steps_are_sorted_inside_the_range_and_a_smooth_weight_has_none() -> None:
    """A box's two steps, sorted; the range clips them; ``r^2 < 2`` steps at sqrt(2); a smooth weight has none; a
    step on a non-polynomial argument is refused with production's scope-edge message. Succeeds
    ``test_r7b2_6_symbolic_steps``. First red (``[M]`` 2026-10-10): the clipping to the range dropped."""
    r = Symbolic.r
    box = Symbolic.of(sympy.Piecewise((1, (r >= sympy.Rational(3, 5)) & (r < sympy.Rational(7, 10))), (0, True)))
    assert _steps(box, 0.0, 2.0) == pytest.approx((0.6, 0.7), abs=0, rel=1e-15)
    assert _steps(box, 0.0, 0.65) == pytest.approx((0.6,), abs=0, rel=1e-15)
    assert _steps(Symbolic.of(sympy.Piecewise((1, r**2 < 2), (0, True))), 0.0, 2.0) == pytest.approx((math.sqrt(2.0),), rel=1e-15)
    assert _steps(Symbolic.of(r**2 + 1), 0.0, 2.0) == ()
    with pytest.raises(ValueError, match="polynomial in r"):
        _steps(Symbolic.of(sympy.Piecewise((1, sympy.sin(3 * r) > 0), (0, True))), 0.0, 2.0)


def test_steps_at_algebraic_and_transcendental_locations() -> None:
    """A step at pi/4 and at sqrt(2)/2 (SymPy's polynomial root finder over QQ[pi] or EX raised on these, qa F3).
    Succeeds ``test_r7b2_2_3_steps_at_algebraic_and_transcendental_locations``."""
    r = Symbolic.r
    quarter_pi = Symbolic.of(sympy.Piecewise((1, r < sympy.pi / 4), (0, True)))
    root_half = Symbolic.of(sympy.Piecewise((1, r < sympy.sqrt(2) / 2), (0, True)))
    assert _steps(quarter_pi, 0.0, 2.0) == pytest.approx((math.pi / 4,), rel=1e-15)
    assert _steps(root_half, 0.0, 2.0) == pytest.approx((math.sqrt(2.0) / 2,), rel=1e-15)


@pytest.mark.parametrize("make", [lambda r: sympy.arg(r - 1), lambda r: sympy.atan2(r - 1, 0)], ids=["arg", "atan2"])
def test_an_unlisted_non_smooth_construct_is_refused(make) -> None:
    """``Symbolic.steps`` is an ALLOW-LIST: a construct whose jumps it does not locate (``arg(r - 1)`` jumps at 1,
    ``atan2(r - 1, 0)`` too) is refused, never integrated across silently. Succeeds
    ``test_r7b2_2_3_an_unlisted_non_smooth_construct_is_refused``."""
    with pytest.raises(ValueError, match=_UNLOCATED_REFUSAL):
        _steps(Symbolic.of(make(Symbolic.r)), 0.0, 2.0)


def test_the_mesh_integrates_a_step_at_pi_over_4() -> None:
    """``Mesh1D.cell_integrals`` of ``r < pi/4`` on the A|B|A sphere's mesh totals 4 pi / 3 (pi / 4)^3, to 64 ulp
    (qa F3, the regression on the mesh's exact integrals). Succeeds
    ``test_r7b2_2_3_the_mesh_integrates_a_step_at_pi_over_4``."""
    from orpheus.mesh import CellsByCount, Mesher
    from tests.gates.sn.verification.analytical._aba_reference import aba_geometry

    mesh = Mesher(aba_geometry(CoordSystem.SPHERICAL)).partition(tuple(CellsByCount.uniform_width(n) for n in (2, 4, 2))).mesh
    weight = Symbolic.of(sympy.Piecewise((1, Symbolic.r < sympy.pi / 4), (0, True)))
    total = float(np.sum(mesh.cell_integrals(weight)))
    expected = 4.0 / 3.0 * math.pi * (math.pi / 4) ** 3
    assert abs(total - expected) <= 64 * math.ulp(expected), (total, expected)
