r"""The angular closure is ONE object, and every per-cell consumer reaches IT.

The P4.9a §6c witness (``.claude/plans/streaming_path_says_what_it_is.md``,
phase P4.9a; spec: ``scratch/p4_9a_verification_plan.md`` §2).

⛔ RED BEFORE P4.9a, by construction — and that is this gate's whole point.
Before the un-weld the degenerate cylindrical-axis branch of
``_OneDimScanWalk._run`` closed the angular axis by calling
``DiamondDifference.update``, whose Morel-Montry block was an inline
Pattern-2 twin of ``sn/angular/closure.py`` — not a closure OBJECT at all,
so leg 1 (``len(seen) > 0``) could not pass however the gate were written.
[M] 2026-08-28, at 5c4f56d7: 8 test files named ``precompute_psi_state``,
3 named ``outgoing_angular_state``, and the intersection was EMPTY — nothing
in the tree compared the two spellings.  [M] the branch is reached ONLY on a
cylindrical mesh with ``n_phi ≡ 2 (mod 4)`` (one bit-exact ``η = 0`` ordinate
per μ-level); at ``n_phi ≡ 0 (mod 4)`` ``DiamondDifference.update`` is called
ZERO times in a full eigenvalue solve, which is why every frozen snapshot in
the tree was blind to the twin.

Deviation from the plan-memo's leg 3, recorded here so nobody restores it:
the ruled P4.9a design (Q1) keeps the NON-degenerate solve fast path on the
closure-MINTED scan constants — no closure method call in the hot loop — so
"the vectorized branch reaches the same object" is witnessed at the MINTING
site (the ``StreamingCoefficientCache`` field gate in
``tests/sn/sweep/core/test_cache.py``) and, per-cell, at the MATVEC's own
degenerate arm (``cell_contribution``), which this gate asserts instead.
The fp(4, 8) solve is kept as a CONTROL: it must fire the march recorder
ZERO times (the hot loop stays call-free).
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry import BC, CoordSystem, Mesh1D
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.sn.operators.streaming import StreamingOperator
from orpheus.transport.fields.angular_boundary_flux import AngularBoundaryFlux
from orpheus.transport.operators.multiplication_operator import (
    MultiplicationOperator,
)
from orpheus.transport.source_sinks import AngularSourceSink
from tests.sn._test_helpers import het_operands, placeholder_materials, sweep_once

pytestmark = [pytest.mark.foundation]

_NG = 2
_N_CELLS = 8


def _cylinder_mesh(n_phi: int) -> SNMesh:
    """Cylinder fixture; ``n_phi ≡ 2 (mod 4)`` activates the per-cell path."""
    mesh = Mesh1D(
        edges=np.linspace(0.01, 2.0, _N_CELLS + 1),
        mat_ids=np.zeros(_N_CELLS, dtype=int),
        coord=CoordSystem.CYLINDRICAL,
        bc_left=BC("reflective"),
        bc_right=BC("vacuum"),
    )
    quad = Quadrature.folded_product(n_mu=4, n_phi=n_phi)
    return SNMesh(mesh, quad, placeholder_materials(ng=_NG))


def _run_one_sweep(sn_mesh: SNMesh) -> None:
    rng = np.random.default_rng(20260828)
    sig_t = rng.uniform(0.3, 3.0, size=(_NG, *sn_mesh.spatial_shape))
    q_values = rng.standard_normal((sn_mesh.quad.N, _NG, *sn_mesh.spatial_shape))
    source = AngularSourceSink(values=q_values, space=sn_mesh.angular_bulk_space)
    boundary = AngularBoundaryFlux.zeros(sn_mesh.angular_trace)
    sweep_once(source, sig_t, sn_mesh, boundary)


def _run_one_matvec(sn_mesh: SNMesh) -> None:
    sig_t, psi, seed = het_operands(sn_mesh)
    L = StreamingOperator.pose(sn_mesh)
    C = MultiplicationOperator.from_mesh(sig_t, sn_mesh)
    if seed is None:
        (L + C).apply(psi)
    else:  # pragma: no cover - cylinder het_operands carries no seed
        from orpheus.numerics.coupled_system import CoupledField

        from tests.sn._test_helpers import joint_m_grid

        grid, _space = joint_m_grid(sn_mesh, L + C)
        grid.apply(CoupledField(systems=(psi, seed)))


def test_every_per_cell_consumer_reaches_the_mesh_closure(monkeypatch):
    """[foundation] Solve + matvec degenerate arms reach ONE closure object.

    Leg 1 (REACHED) is deliberately the first assert: pre-carve it is the
    only statement that can fail (the recorder stays empty because no
    closure method exists on the degenerate solve route).
    """
    sn = _cylinder_mesh(n_phi=6)
    # Leg 0 (P4.9b) — the POSED operator's slot is the hub's instance:
    # the loop-closer between the operator's field and the identity every
    # later leg asserts. Without it the consumers could reach the hub's
    # object while a posed operator held a different one.
    assert StreamingOperator.pose(sn).angular_closure is sn.angular_closure
    closure_cls = type(sn.angular_closure)

    march_seen: list[object] = []
    contrib_seen: list[object] = []

    # ``advance_psi_half`` does not exist pre-carve: ``raising=False`` lets
    # the spy install anyway, and the walk simply never calls it — leg 1
    # then fails on the empty recorder, which is this gate's red reading.
    orig_march = getattr(closure_cls, "advance_psi_half", None)

    def march_spy(self, *args, **kwargs):
        march_seen.append(self)
        assert orig_march is not None, (
            "advance_psi_half was called but has no production body"
        )
        return orig_march(self, *args, **kwargs)

    orig_contrib = closure_cls.cell_contribution

    def contrib_spy(self, *args, **kwargs):
        contrib_seen.append(self)
        return orig_contrib(self, *args, **kwargs)

    monkeypatch.setattr(closure_cls, "advance_psi_half", march_spy, raising=False)
    monkeypatch.setattr(closure_cls, "cell_contribution", contrib_spy)

    _run_one_sweep(sn)

    # Leg 1 — REACHED: the degenerate solve branch called the closure at all.
    assert len(march_seen) > 0, (
        "the degenerate solve branch never reached the angular closure "
        "(pre-P4.9a it marches through DiamondDifference's inline twin)"
    )
    # Leg 2 — IDENTITY: every march call was on the mesh's own closure.
    assert all(obj is sn.angular_closure for obj in march_seen)

    # Leg 3 — the OTHER per-cell consumer: the matvec's degenerate arm
    # reaches the SAME instance through cell_contribution.
    contrib_seen.clear()
    _run_one_matvec(sn)
    assert len(contrib_seen) > 0
    assert all(obj is sn.angular_closure for obj in contrib_seen)

    # Leg 4 — NON-VACUITY: a second mesh's closure is a DIFFERENT object and
    # its calls land on ITS instance — the recorder discriminates identities.
    sn2 = _cylinder_mesh(n_phi=6)
    assert sn2.angular_closure is not sn.angular_closure
    march_before = len(march_seen)
    _run_one_sweep(sn2)
    new_calls = march_seen[march_before:]
    assert len(new_calls) > 0
    assert all(obj is sn2.angular_closure for obj in new_calls)

    # Control — the NON-degenerate fast path stays call-free (fp(4, 8) has
    # zero degenerate ordinates and consumes MINTED constants, not methods).
    march_seen.clear()
    _run_one_sweep(_cylinder_mesh(n_phi=8))
    assert len(march_seen) == 0, (
        "the fast path must consume the closure's minted scan constants, "
        "never a per-cell method call (Q1's perf half)"
    )


def test_spatial_protocol_carries_no_angular_member():
    """[foundation] The L2 visit family is purely spatial (row 3's half).

    Asserted by ``dataclasses.fields`` — not ``hasattr`` on an instance,
    which a defaulted field or ``getattr`` fallback would still answer
    (coding-standards, the string-form clause).
    """
    import dataclasses

    from orpheus.transport.spatial.scheme import (
        CellResult,
        CellVisit,
        UpstreamState,
    )

    def field_names(cls) -> set[str]:
        return {f.name for f in dataclasses.fields(cls)}

    assert field_names(CellVisit).isdisjoint({"tau", "c_in", "c_out"})
    assert "angular_upstream" not in field_names(UpstreamState)
    assert "outgoing_angular_state" not in field_names(CellResult)
