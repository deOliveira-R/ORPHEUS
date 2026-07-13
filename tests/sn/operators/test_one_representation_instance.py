r"""S6.5 — the "one representation instance" discriminating tests (#222 capstone).

The structural payoff of the operator/representation re-layering: the
:class:`~orpheus.sn.loss_representation.LossRepresentation` instance that
:meth:`StreamingOperator.apply` consumes for the matvec :math:`(L+C)\psi` MUST
BE (object identity) the instance that :meth:`InvertibleOperator.solve`
consumes for the forward substitution :math:`(L+C)^{-1}q`.  This makes the L21
invariant ("matvec ≡ sweep — two actions of ONE operator") a TYPE FACT
enforced by construction, closing coding-elegance Smell #16 (two doors to one
operator's representation).

NEGATIVE PRE-CONDITION (why these are tripwires, not tautologies): at the
pre-S6.5 HEAD the doors each called ``default_for(sn_mesh)`` INDEPENDENTLY —
``apply`` via the ``loss_representation`` cached_property, ``solve`` via
``transport_sweep`` → ``default_for`` (a FRESH frozen-dataclass instance per
sweep call), and the G-S resolvent via a THIRD ``default_for`` per solve.  The
spies below capture the CALL-TIME ``self`` of each door, so the tests observed
genuinely different objects and FAILED at pre-S6.5 HEAD (recorded in the S6.5
commit body); they pass once both doors consume the operator's ONE instance.

Spying call-time ``self`` (the S6.4 one-walk SPY pattern) rather than
comparing property handles keeps the test non-tautological: a property defined
as delegation would make ``A.loss_representation is L.loss_representation``
true by spelling; what the spy proves is that the SOLVE PATH actually runs on
that object.

NOTE the scope boundary DISSOLVED at the coupled-block campaign step 6
(R-6.1): the module-level ``transport_sweep`` operator-free entry — which
legitimately selected via ``default_for`` per call — is retired; every
production solve now rides an operator (test callers build a fresh one via
``tests.sn._test_helpers.sweep_once``).  The one-instance theorem is about a
LONG-LIVED operator's doors, and is now the only selection discipline left.

``-O``-safe (vv Mode 8): ``pytest.fail`` only — no bare asserts.
``foundation`` — software-structure invariants (no theory ``:label:``).
"""
from __future__ import annotations

import pytest

from orpheus.sn.operators.streaming import StreamingOperator
from orpheus.transport.operators.multiplication_operator import MultiplicationOperator
from tests.sn._test_helpers import cart2d_2g_nonsquare, het_operands

pytestmark = pytest.mark.foundation


def _spy_capture(monkeypatch, cls, method_name: str, seen: dict, key: str):
    """Wrap ``cls.method_name`` to record the call-time ``self`` under ``key``."""
    real = getattr(cls, method_name)

    def spy(self, *args, **kwargs):
        seen[key] = self
        return real(self, *args, **kwargs)

    monkeypatch.setattr(cls, method_name, spy)


def test_apply_and_solve_share_one_representation_instance(monkeypatch):
    """[L0 structural] the representation ``apply`` runs on IS the one
    ``solve`` runs on IS the operator's ``loss_representation`` handle.

    2-D Cartesian is the discriminating config (gate memo
    ``s6_relayering_verification.md`` §4): at d=2 the two doors route through
    genuinely different actions (the loss_action matvec vs the scheduled
    sweep), so a fresh-construction door is observable as a distinct
    frozen-dataclass instance.  The spy class is the PRODUCTION d=2
    representation — ScanMarch since the S6.9 Fork-B2 default flip
    (MovingFrontierWindow through S6.9; re-point the spies if the default
    moves again — the ``seen`` capture failing empty is the tripwire).
    """
    from orpheus.sn.loss_representation import ScanMarch

    seen: dict[str, object] = {}
    _spy_capture(monkeypatch, ScanMarch, "loss_action", seen, "apply")
    _spy_capture(monkeypatch, ScanMarch, "sweep", seen, "solve")

    sn = cart2d_2g_nonsquare()
    sig_t, psi = het_operands(sn)[:2]
    L = StreamingOperator(sn)
    C = MultiplicationOperator.from_mesh(sig_t, sn)
    A = L + C

    _ = L.apply(psi)
    _ = A.solve(psi)

    rep_apply = seen.get("apply")
    rep_solve = seen.get("solve")
    if rep_apply is None or rep_solve is None:
        pytest.fail(
            f"spy capture incomplete (apply={rep_apply!r}, solve={rep_solve!r}) "
            "— a door no longer routes through ScanMarch; the spy "
            "targets need re-pointing at the production representation."
        )
    if rep_apply is not rep_solve:
        pytest.fail(
            "S6.5 NOT satisfied: StreamingOperator.apply and "
            f"InvertibleOperator.solve ran on DIFFERENT LossRepresentation "
            f"instances (apply id={id(rep_apply)}, solve id={id(rep_solve)}). "
            "L21 (matvec ≡ sweep) is a coincidence, not a type fact — the "
            "two-doors smell (#222 S6.0) is open again."
        )
    if rep_apply is not L.loss_representation:
        pytest.fail(
            "the shared call-time instance is NOT the operator's "
            "loss_representation cached_property — a third construction site "
            "is feeding both doors."
        )
    if A.loss_representation is not L.loss_representation:
        pytest.fail(
            "InvertibleOperator.loss_representation does not delegate to its "
            "streaming leaf's instance — the composite's handle is a separate "
            "construction."
        )


def test_scheduled_solve_runs_the_operators_instance(monkeypatch):
    """[L0 structural] the reified ``M = (L+C−B_lower)`` scheduled solve's
    interior kernel is BOUND to the operator's one instance (door 3 of the
    S6.5 inventory; migrated from the dissolved ``_GaussSeidelResolvent``
    at #226 step 2 — retirement-means-test-migration).

    Pre-S6.5 the G-S path passed ``interior=default_for(sn_mesh)
    ._sweep_interior`` — a fresh instance per scheduled solve.  Post-#226
    the scheduled walk routes through the representation's OWN ``sweep``
    door (which injects its own kernel); the spy captures the ``self`` of
    every ``_sweep_interior`` execution during ``M.solve`` and demands it
    IS ``A.loss_representation`` (= ``M.loss_representation``, the
    invertible operand's one instance).
    """
    from orpheus.sn.operators.boundary import SNBoundaryOperator
    from orpheus.sn.loss_representation import ScanMarch
    from orpheus.sn.loss_representation.sweep_schedule import SweepSchedule

    captured: list[object] = []
    real = ScanMarch._sweep_interior

    def spy(self, *args, **kwargs):
        captured.append(self)
        return real(self, *args, **kwargs)

    monkeypatch.setattr(ScanMarch, "_sweep_interior", spy)

    sn = cart2d_2g_nonsquare()
    sig_t, psi = het_operands(sn)[:2]
    L = StreamingOperator(sn)
    C = MultiplicationOperator.from_mesh(sig_t, sn)
    A = L + C
    M = A - SNBoundaryOperator(sn).split(SweepSchedule.gauss_seidel(sn)).lower

    _ = M.solve(psi)

    if not captured:
        pytest.fail(
            "M.solve never executed ScanMarch._sweep_interior — the spy "
            "target needs re-pointing at the production interior kernel."
        )
    rep = A.loss_representation
    if M.loss_representation is not rep:
        pytest.fail(
            "M.loss_representation must BE the invertible operand's one "
            "instance (S6.5)."
        )
    strays = [obj for obj in captured if obj is not rep]
    if strays:
        pytest.fail(
            f"S6.5 NOT satisfied: {len(strays)}/{len(captured)} scheduled "
            "interior-kernel executions ran on a representation instance "
            "that is NOT the operator's loss_representation (door 3 of the "
            "S6.5 inventory is open)."
        )
