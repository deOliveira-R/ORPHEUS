r"""S6.4 — the "one octant walk" discriminating tests (#222).

The structural payoff of the unified kernel-parameterized octant walk: a SWEEP
(:math:`(L+C)^{-1} q`) and a LOSS_ACTION (:math:`(L+C)\psi`, the matvec) on the
SAME 2-D Cartesian mesh MUST exercise the SAME
:meth:`~orpheus.sn.loss_representation._OctantWalk._interior_walk` traversal —
they fork only at the cell kernel (solve vs apply) + the emit policy, NEVER in
a duplicated octant frame.

Staging (the gate memo ``s6_4_unified_walk_verification.md`` §4):

* **sub-step (a)**: ``_OctantWalk`` landed with BOTH matvec variants
  (window + scan-march) routed through it; the SPY test was
  ``xfail(strict=True)`` (the sweep door still used ``sweep_octant_group``'s
  private frame — a recorded gap).
* **sub-step (b)** (this commit): the sweep frames came into the walk
  (``_OctantWalk.sweep_group`` × the representations' ``_sweep_interior``
  kernels; ``sweep_octant_group`` + ``_sweep_2d_scanmarch`` retired) → the
  SPY flipped xfail → PASS and the marker was removed.

``-O``-safe (vv Mode 8): ``pytest.fail`` only — no bare asserts.
``foundation`` — software-structure invariants (no theory ``:label:``).
"""
from __future__ import annotations

import pytest

from orpheus.sn.operator import CollisionOperator, StreamingOperator
from tests.sn._test_helpers import cart2d_2g_nonsquare, het_operands

pytestmark = pytest.mark.foundation


def test_sweep_and_loss_action_hit_one_octant_walk(monkeypatch):
    """[L0 structural] sweep AND loss_action exercise the SAME
    ``_OctantWalk._interior_walk`` (the one-walk claim, observed at runtime).

    NEGATIVE PRE-CONDITION (why this is a tripwire, not a tautology): before
    S6.4(b) the octant orchestration was duplicated in lockstep —
    ``sweep_octant_group`` (sweep) vs the matvec frames — so the sweep door
    did NOT reach the shared walk and this test FAILED (landed
    ``xfail(strict=True)`` at sub-step (a), the recorded gap).  Since (b)
    both doors hit the ONE walk: L21 (matvec ≡ sweep) is a CODE FACT.
    """
    from orpheus.sn.loss_representation import _OctantWalk

    hits: list[str] = []
    real_walk = _OctantWalk._interior_walk

    def spy(self, *args, **kwargs):
        hits.append("walked")
        return real_walk(self, *args, **kwargs)

    monkeypatch.setattr(_OctantWalk, "_interior_walk", spy)

    sn = cart2d_2g_nonsquare()
    sig_t, psi = het_operands(sn)
    L = StreamingOperator(sn)
    C = CollisionOperator(sn, sig_t)
    A = L + C

    # (1) the MATVEC direction routes through the shared walk:
    hits.clear()
    _ = L.apply(psi)
    if not hits:
        pytest.fail(
            "loss_action did NOT route through _OctantWalk._interior_walk — "
            "the matvec still uses a private duplicated octant frame "
            "(the S6.4(a) carve regressed)."
        )

    # (2) the SWEEP direction routes through the SAME shared walk:
    hits.clear()
    _ = A.solve(psi)
    if not hits:
        pytest.fail(
            "sweep did NOT route through _OctantWalk._interior_walk — a "
            "private sweep octant frame has re-appeared (the S6.4(b) "
            "one-walk unification regressed)."
        )
    # Both directions hit the ONE walk → L21 (matvec ≡ sweep) is a CODE FACT.


def test_octant_walk_is_kernel_parameterized_not_boolean():
    """[L0 structural] ``_OctantWalk`` forks on a KERNEL OBJECT + EMIT policy,
    NOT a boolean ``is_solve``/``is_apply`` flag (coding-elegance Smell #3).

    A boolean dispatch (``if solve: cell_kernel_batch else
    residual_kernel_batch``) inside the shared walk would re-introduce the
    twin-path branch the carve eliminates.  The walk's interior MUST receive
    its cell operation as a callable/object parameter.  Source-inspection by
    necessity (the anti-pattern is a CODE shape, not a runtime value) — via
    the AST's identifiers, so docstrings/comments that NAME the anti-pattern
    don't trip it.  ``-O``-safe.
    """
    import ast
    import inspect
    import textwrap

    from orpheus.sn.loss_representation import _OctantWalk

    smells = {"is_solve", "is_apply", "is_matvec"}
    tree = ast.parse(textwrap.dedent(inspect.getsource(_OctantWalk)))
    identifiers = {
        node.id for node in ast.walk(tree) if isinstance(node, ast.Name)
    } | {
        node.attr for node in ast.walk(tree) if isinstance(node, ast.Attribute)
    } | {
        node.arg for node in ast.walk(tree)
        if isinstance(node, (ast.arg, ast.keyword)) and node.arg is not None
    }
    offenders = sorted(identifiers & smells)
    if offenders:
        pytest.fail(
            f"_OctantWalk carries boolean direction flag(s) {offenders} — the "
            "carve degraded into the boolean-flag anti-pattern "
            "(coding-elegance Smell #3).  The cell operation MUST be a "
            "kernel OBJECT/callable parameter, not a flag."
        )


def test_both_matvec_variants_share_the_walk(monkeypatch):
    """[L0 structural] BOTH matvec interior kernels (window + scan-march)
    route through the ONE ``_OctantWalk._interior_walk`` — the (a)-level
    one-frame fact (the former Fork-B1 lockstep twins are gone).

    Passes from sub-step (a) — unlike the sweep SPY above, which waits
    for (b).
    """
    from orpheus.sn.loss_representation import (
        MovingFrontierWindow, ScanMarch, _OctantWalk,
    )

    hits: list[str] = []
    real_walk = _OctantWalk._interior_walk

    def spy(self, *args, **kwargs):
        hits.append("walked")
        return real_walk(self, *args, **kwargs)

    monkeypatch.setattr(_OctantWalk, "_interior_walk", spy)

    sn = cart2d_2g_nonsquare()
    sig_t, psi = het_operands(sn)

    for rep_cls in (MovingFrontierWindow, ScanMarch):
        hits.clear()
        _ = rep_cls(sn).loss_action(sig_t, psi)
        if not hits:
            pytest.fail(
                f"{rep_cls.__name__}.loss_action did NOT route through "
                "_OctantWalk._interior_walk — its octant frame is a private "
                "duplicate again (the S6.4(a) carve regressed)."
            )
