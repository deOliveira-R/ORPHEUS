r"""S6.4 — the "one octant walk" discriminating tests (#222).

The structural payoff of the unified kernel-parameterized octant walk: a SWEEP
(:math:`(L+C)^{-1} q`) and a LOSS_ACTION (:math:`(L+C)\psi`, the matvec) on the
SAME 2-D Cartesian mesh MUST exercise the SAME
:meth:`~orpheus.sn.loss_representation._OctantWalk._interior_walk` traversal —
they fork only at the cell kernel (solve vs apply) + the emit policy, NEVER in
a duplicated octant frame.

Staging (the gate memo ``s6_4_unified_walk_verification.md`` §4):

* **sub-step (a)** (this commit): ``_OctantWalk`` exists and BOTH matvec
  variants (window + scan-march) route through it.  The matvec leg of the SPY
  test passes; the SWEEP leg still fails (``sweep_octant_group`` keeps its
  private frame) → the SPY test is ``xfail(strict=True)`` — a recorded gap.
* **sub-step (b)**: the sweep frames come into the walk → the SPY flips
  xfail → xpass; REMOVE the marker in that commit.

``-O``-safe (vv Mode 8): ``pytest.fail`` only — no bare asserts.
``foundation`` — software-structure invariants (no theory ``:label:``).
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry import BC, CoordSystem, Mesh2D
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.geometry import SNMesh
from orpheus.sn.operator import CollisionOperator, StreamingOperator
from orpheus.transport.fields.angular_flux import AngularFlux
from orpheus.transport.fields.boundary_flux import BoundaryFlux
from orpheus.transport.timed_full_field import TimedFullField
from tests.sn._test_helpers import placeholder_materials

pytestmark = pytest.mark.foundation


def _cart2d_2g_nonsquare(nx: int = 5, ny: int = 7) -> SNMesh:
    """2-D Cartesian, reflective, 2G, NON-SQUARE (the x↔y-swap moat) —
    the discriminating config: the octant frame + pure-z branch + interior
    recurrence all degenerate on a 1G flat square box (vv §H1/§H2)."""
    mesh = Mesh2D(
        edges_x=np.linspace(0.0, 2.0, nx + 1),
        edges_y=np.linspace(0.0, 3.0, ny + 1),
        mat_map=np.zeros((nx, ny), dtype=int),
        coord=CoordSystem.CARTESIAN,
        bc_xmin=BC("reflective"), bc_xmax=BC("reflective"),
        bc_ymin=BC("reflective"), bc_ymax=BC("reflective"),
    )
    return SNMesh(mesh, Quadrature.level_symmetric(4), placeholder_materials(ng=2))


def _het_operands(sn: SNMesh):
    """Heterogeneous σ_t + a non-flat random state (≥2G, non-degenerate)."""
    rng = np.random.default_rng(20260611)
    sig_t = rng.uniform(0.3, 3.0, size=(sn.ng, *sn.spatial_shape))
    psi = TimedFullField.zeros(bulk=AngularFlux, boundary=BoundaryFlux, mesh=sn)
    psi.bulk.values[...] = rng.standard_normal(psi.bulk.values.shape)
    for face in psi.boundary.layout.faces:
        fv = psi.boundary.face_view(face)
        fv[...] = rng.standard_normal(fv.shape)
    return sig_t, psi


@pytest.mark.xfail(
    strict=True,
    reason=(
        "S6.4(b) pending: the SWEEP direction still uses sweep_octant_group's "
        "private octant frame — only the matvec routes through "
        "_OctantWalk._interior_walk at sub-step (a).  Remove this marker in "
        "the (b) commit (the recorded-gap discipline)."
    ),
)
def test_sweep_and_loss_action_hit_one_octant_walk(monkeypatch):
    """[L0 structural] sweep AND loss_action exercise the SAME
    ``_OctantWalk._interior_walk`` (the one-walk claim, observed at runtime).

    NEGATIVE PRE-CONDITION (why this is a tripwire, not a tautology): before
    S6.4(b) the octant orchestration is duplicated in lockstep —
    ``sweep_octant_group`` (sweep) vs the matvec frames — so the sweep door
    does NOT reach the shared walk and this test FAILS (the xfail above
    records the gap).  After (b) both doors hit the ONE walk and it passes.
    """
    from orpheus.sn.loss_representation import _OctantWalk

    hits: list[str] = []
    real_walk = _OctantWalk._interior_walk

    def spy(self, *args, **kwargs):
        hits.append("walked")
        return real_walk(self, *args, **kwargs)

    monkeypatch.setattr(_OctantWalk, "_interior_walk", spy)

    sn = _cart2d_2g_nonsquare()
    sig_t, psi = _het_operands(sn)
    L = StreamingOperator(sn, sig_t)
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
            "sweep did NOT route through _OctantWalk._interior_walk — the "
            "sweep still uses sweep_octant_group's private octant frame "
            "(the one-walk unification, S6.4(b), is open)."
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

    sn = _cart2d_2g_nonsquare()
    sig_t, psi = _het_operands(sn)
    L = StreamingOperator(sn, sig_t)

    for rep_cls in (MovingFrontierWindow, ScanMarch):
        hits.clear()
        _ = rep_cls(sn).loss_action(L, psi)
        if not hits:
            pytest.fail(
                f"{rep_cls.__name__}.loss_action did NOT route through "
                "_OctantWalk._interior_walk — its octant frame is a private "
                "duplicate again (the S6.4(a) carve regressed)."
            )
