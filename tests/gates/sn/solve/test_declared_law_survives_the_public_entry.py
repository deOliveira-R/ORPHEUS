r"""A boundary condition DECLARED for the faces reaches a public fixed-source solve (wave 2 re-posing).

The declaration channel, landed alongside campaign phase **P4** of
`.claude/plans/archive/affine_boundary_source_channel.md`, and re-posed at the wave-2 carve onto the conditions package.
The channel's claim is unchanged: **a declared condition must change a public solve**, and it is a claim about the
ARRIVAL of the declaration, not about the detail of the delivered ``q``.

What changed is the mechanism, and with it what can go wrong:

1. **A face carries a role, never a law.** A face is glued to a partner (a deck) or faces the exterior (an ``int`` tag); the
   condition a tag stands for is declared at the solve, by the keyword-only ``boundary_conditions`` of every entry, as a
   :class:`~orpheus.transport.boundary_conditions.WallCondition` (an affine pair: a response network and an inflow). The old
   split between a ``BC`` tag, which could not carry a function, and an already-typed law, which could, is gone with the string
   tag: there is ONE shared exterior type and its facades.

2. ⭐ **Every public solver entry point rebuilds the method mesh**
   (:func:`~orpheus.sn.solver.solve_sn_fixed_source` calls ``_as_problem(...)``, which constructs a fresh
   :class:`~orpheus.sn.problem.SNProblem`), so anything installed on a constructed one would be DISCARDED. The conditions
   therefore ride the ENTRY, which is what the public API already takes: ``boundary_conditions=`` is forwarded by
   ``_as_problem`` to the binding of the new problem, and no ``problem=`` parameter is needed anywhere.

⚠ **What this module does NOT claim.** It pins that the declaration ARRIVES and changes the answer, not that the delivered
``q`` is correct in detail. The value claim is
``test_declared_inflow_reaches_the_rhs.py::test_the_declared_boundary_law_holds_on_the_answer`` (``gamma_-psi|_f == q_f``,
ERR-075), which reads the boundary law off the converged answer. Kept separate because they fail for different reasons: that
one for a mis-delivered source, this one for a declaration that never arrives.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry import StructuredGeometry
from orpheus.geometry.deck import MirrorDeck
from orpheus.mesh import CellsByCount, Mesher
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.solver import solve_sn_fixed_source
from orpheus.transport.boundary_conditions import BoundaryRefused, RefusalAxis, WallCondition
from tests.gates.sn._test_helpers import placeholder_materials
from tests.gates.transport._boundary_conditions_helpers import local

pytestmark = pytest.mark.l1

_VALUE = 2.5
_N_ORD = 8
_NG = 2
_NX = 8


def _solve(xmin_condition, inner: str = "source_iteration") -> np.ndarray:
    """Drive the FULLY PUBLIC entry point — no private helper anywhere.

    The left face is the tag 0, bound to ``xmin_condition``; the right face is the tag 1, a vacuum. A mirror DECK on the left
    (``xmin_condition`` a deck) needs no condition.
    """
    left = xmin_condition if isinstance(xmin_condition, MirrorDeck) else 0
    geom = StructuredGeometry.slab((0.0, 2.0), (0,), left=left, right=1)
    mesh = Mesher(geom).partition(CellsByCount.uniform_volume(_NX)).mesh
    quad = Quadrature.gauss_legendre(n_ordinates=_N_ORD)
    source = np.ones((quad.n_ordinates, _NG, len(mesh.centers)))
    conditions = {1: WallCondition.vacuum()}
    if left == 0:
        conditions[0] = xmin_condition
    solution = solve_sn_fixed_source(
        placeholder_materials(ng=_NG), mesh, quad, source,
        inner_solver=inner, boundary_conditions=conditions,
    )
    return np.asarray(solution.scalar_flux.values)


class TestTheDeclarationSurvivesTheRebuild:
    r"""⭐ The keystone: a declared condition changes a public solve."""

    @pytest.mark.parametrize("inner", ["source_iteration", "krylov"])
    def test_a_declared_inflow_changes_the_answer(self, inner: str) -> None:
        r"""`[M]` ``|φ_declared − φ_vacuum|_inf = 1.877290`` on this fixture (HEAD ``5d04c296``).

        The row that fails if the declaration is dropped anywhere along
        entry → ``_as_problem`` → ``bind`` → ``realize_boundary`` → ``from_binding`` → RHS. Before the channel it was dropped
        at the first step and would have been dropped again at ``_as_problem``.

        Both inner solvers, because they consume the composite RHS through
        different bodies — and because Krylov was unreachable with a declared
        prescribed inflow at all until P3 retired the affine operator.
        """
        declared = _solve(WallCondition.uniform(_VALUE), inner)
        vacuum = _solve(WallCondition.vacuum(), inner)
        moved = float(np.max(np.abs(declared - vacuum)))
        assert moved > 1e-6, (
            f"[{inner}] a DECLARED uniform inflow left the answer identical "
            f"to vacuum (|Δφ| = {moved:.3e}) — the declaration was discarded "
            f"somewhere between the entry and the solve, which is exactly "
            f"the defect this channel closes."
        )

    def test_the_spellings_of_a_vacuum_are_one_float_program(self) -> None:
        r"""⭐ The facades are one value, not several.

        ``WallCondition.vacuum()``, ``specular(0.0)`` and ``white(0.0)`` name the same physics and are ONE declared value (one
        digest), realize the same ``ZeroMap``, and so must produce the same float program. BIT-identical, not merely close: the
        facades' only extra work is the construction of the network, so any difference means the realized maps diverged after
        the parse.

        This is the row that would catch a realizer that quietly treated a zero-weight wall differently from the empty network —
        the failure mode a "does it change the answer?" test cannot see.
        """
        reference = _solve(WallCondition.vacuum())
        for spelling in (WallCondition.specular(0.0), WallCondition.white(0.0)):
            np.testing.assert_array_equal(
                reference, _solve(spelling),
                err_msg=(
                    f"declaring {spelling!r} and declaring the vacuum gave different answers — the facades of one "
                    f"value have diverged"
                ),
            )


class TestTheChannelDoesNotWeakenTheOtherDeclarations:
    """The conditions and the decks are both first-class; a declared source is an addition, not a detour."""

    @pytest.mark.parametrize("role", [WallCondition.vacuum(), MirrorDeck()], ids=["vacuum", "mirror"])
    def test_a_condition_or_a_deck_still_resolves(self, role) -> None:
        """A vacuum condition and a mirror deck both drive a solve to a finite answer.

        Cheap, but it is the row that reddens if the declared-source path were placed so that it swallowed the other
        declarations too (e.g. a ``binding.sources`` read broad enough to match every tagged face).
        """
        phi = _solve(role)
        assert np.all(np.isfinite(phi))
        assert float(np.max(np.abs(phi))) > 0.0

    def test_an_unregistered_condition_still_refuses(self) -> None:
        r"""⭐ The declared-source path must not become a hole in the admission table.

        SN's realizer table gates which conditions it admits. A method-only condition SN registers no realizer for (a return
        fraction below zero, as diffusion's zero flux) is refused on the STATE_CONE axis by the table's unregistered-unit
        door, and the refusal rides the same ``boundary_conditions=`` channel: it must not sail through the binding and fail
        later, or worse, silently.
        """
        with pytest.raises(BoundaryRefused) as info:
            _solve(local().NegativeProbe())
        assert info.value.admission.axis is RefusalAxis.STATE_CONE
