r"""Production pairs typed members: the formal path is the one a solve takes (the route gate of the space-membership unit U2).

The user's ruling R1 (2026-10-30, plan ``.claude/plans/space_membership_covector_gauge.md``): only the formal machinery,
a typed member handed to a space verb, is guaranteed; a bare array passes only the shape leg of the parse, so an
equal-shape array of another space is indistinguishable from a member. That guarantee covers production only if
production is on the formal path, and nothing in a green solve says so.
This file says so: :meth:`~orpheus.numerics.space.FunctionSpace.coordinates_of`, the exact-shape membership verb of the
scalar verbs (``inner_product``, ``norm``) and of every co-vector pairing and mint (``Covector.evaluate``, ``flat``), is
made to refuse a bare array, and one eigenvalue solve per sharp method must complete.

Only the pairings are held to it. The action verb ``leading_coordinates_of`` (the metric legs, an operator's matvec,
the Riesz legs, ``.H``) legitimately carries a trailing batch, ``A (x) I_P``, and is not patched. The exact verb has one
action caller, the plain binding's admission
(:func:`~orpheus.transport.operators.lift.admit_array`, since the user's ruling D1, 2026-10-30): an operator acts on the
operator algebra's bare iterate (#664), so its calls are exempted, identified by the caller's CODE OBJECT (not a name
match), and not counted. The exemption cannot swallow a pairing only while every caller of ``admit_array`` is an
action verb; :func:`test_the_exempted_door_is_reached_only_by_actions` asserts that over ``orpheus/`` by AST (`[M]`
2026-10-30: 6 of 6 callers, ``apply``/``apply_transpose``/``solve`` of the isotropic transfer and the multiplier).

First red: at the parent of U2c (``64cc3a7a``) every one of the 7 exercised sites passed an ``ndarray`` (`[M]` the spec's
caller census); after U2c's typed callers and before ``FunctionSpace.integral``, the two ``flat(unit)`` sites
(``diffusion/solver.py`` and ``homogeneous/solver.py``) still handed the space its own unit as a bare array
(`[M]` 2026-10-30, a probe recording the raw callers' frames: 1 call each, SN 0). Mutation: re-spell either site
``space.flat(space.unit)``. The liveness count (typed calls > 0) is the positive control: a spy that never fires would
pass every row.
"""
from __future__ import annotations

import ast
import pathlib
import sys

import pytest

from orpheus.numerics import space as space_module
from orpheus.numerics.space import Member
from orpheus.transport.operators import lift

pytestmark = pytest.mark.foundation


@pytest.fixture
def typed_calls(monkeypatch) -> list[int]:
    """Make the exact-shape membership verb refuse a bare array from a pairing; count the typed calls it admits.

    A call whose immediate caller is the plain binding's admission (an action, #664's bare iterate) is passed through
    and not counted."""
    honest = space_module.FunctionSpace.coordinates_of
    action = lift.admit_array.__code__
    count = [0]

    def strict(self, x):
        if sys._getframe(1).f_code is action:
            return honest(self, x)
        if not isinstance(x, Member):
            raise AssertionError(f"a bare {type(x).__name__} reached {self.name!r}.coordinates_of from production")
        count[0] += 1
        return honest(self, x)

    monkeypatch.setattr(space_module.FunctionSpace, "coordinates_of", strict)
    return count


def test_homogeneous(typed_calls):
    from orpheus.data.macro_xs.synthetic import get_mixture
    from orpheus.homogeneous.solver import solve_homogeneous_infinite

    solve_homogeneous_infinite(get_mixture("A", "2g"))
    assert typed_calls[0] > 0, "the spy never fired: the homogeneous solve reached no exact membership verb"


def test_diffusion_eigen(typed_calls):
    from tests.gates.diffusion import _wave2_helpers as H

    H.solve(*H.build("slab", ("vacuum",)))
    assert typed_calls[0] > 0, "the spy never fired: the diffusion solve reached no exact membership verb"


def test_sn_eigen_vacuum(typed_calls):
    """A 2-group two-region vacuum slab: reaches the production rate, ``compute_keff`` and the leakage bridge."""
    from orpheus.sn.solver import solve_sn
    from tests.gates.sn.architecture import _sn_consumers_helpers as S
    from tests.gates.transport import _carry_cases as carry

    mesh, quad, conditions = S.build("slab", ("vacuum",))
    solve_sn(carry.materials(), mesh, quad, boundary_conditions=conditions, **carry.EIGEN)
    assert typed_calls[0] > 0, "the spy never fired: the SN solve reached no exact membership verb"


def test_the_exempted_door_is_reached_only_by_actions():
    """The exemption above is sound only while ``admit_array`` is an action's door: a pairing (an inner product, a
    co-vector's evaluate) re-spelled through it would pass the spy unseen. Every production call of it sits in an
    action verb, by AST over ``orpheus/``; the positive control is the census itself (6 callers known today)."""
    root = pathlib.Path(lift.__file__).resolve().parents[2]
    callers = []
    for f in sorted(root.rglob("*.py")):
        for fn in ast.walk(ast.parse(f.read_text())):
            if isinstance(fn, (ast.FunctionDef, ast.AsyncFunctionDef)):
                for n in ast.walk(fn):
                    if isinstance(n, ast.Call) and getattr(n.func, "id", getattr(n.func, "attr", None)) == "admit_array":
                        callers.append((f.relative_to(root).as_posix(), fn.name))
    assert len(callers) >= 6, f"the census lost its known callers: {callers}"
    stray = [c for c in callers if c[1] not in {"apply", "apply_transpose", "solve"}]
    assert not stray, f"admit_array, the pairing spy's exempted door, is reached from a non-action: {stray}"
