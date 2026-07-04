r"""Generate the FROZEN pre-carve 1-D walk matvec baselines (Phase 2.5 S0.3).

The #280 apply-loop unification (2.5a) claims a BIT-IDENTICAL relocation of
``_OneDimScanWalk._apply_walk`` (fwd matvec) and ``loss_action_transpose``
(adjoint matvec) into one orientation-parametrized loop frame. The existing
``array_equal`` matvec canaries are SELF-REFERENTIAL (they compare
``op.apply[_transpose]`` against the SAME ``loss_action[_transpose]`` body on
a fresh operator instance — override-not-leak discriminators), so a
relocation that moves BOTH the SUT and the reference together leaves them
green even if values shifted. These snapshots are the structurally
independent anchor: captured BEFORE the relocation, they cannot move with it
(gate-spec ``a3_solve_transpose_verification.md`` §10.2).

Scope + lifetime:

* ``(L+C)`` only — the object the 2.5a carve rebuilds (``−B`` is a separate
  leaf, untouched).
* The curvilinear rows freeze the LAGGED Carlson-seed formulation. Per
  ruling R10 the #282 fix (2.5d) re-poses the curvilinear seed — the sphere
  / cylinder rows are then RE-CAPTURED (principled re-baseline, recorded in
  the roadmap); the Cartesian slab row stays bitwise through the whole
  phase.
* Disposition at 2.5e close-out: keep as permanent walk-value canaries or
  retire once the spy + oracle chain supersedes them — decided EXPLICITLY
  (the fuller-view rule).

Every case's mesh/operator/input construction lives HERE (one source; the
test re-runs the same build). Inputs are seeded ``_random_composite`` fields
(non-flat bulk AND boundary — a zero boundary would null the boundary
in↔out swap the transpose exercises).

Run:  ``.venv/bin/python -m tests.sn.regression._generate_walk_baselines``
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np

from orpheus.sn.operators.streaming import StreamingOperator
from orpheus.transport.operators.multiplication_operator import MultiplicationOperator
from tests.sn.operators.test_g_adjoint_reciprocity import (
    _make_cyl,
    _make_slab,
    _make_sphere,
    _random_composite,
)
from tests.sn.regression._generate_snapshots import SNAPSHOT_DIR


@dataclass(frozen=True)
class WalkMatvecCase:
    """One frozen-baseline row: a named build of (sn, LC, psi, phi)."""

    name: str
    build: Callable[[], tuple]

    @property
    def snapshot_name(self) -> str:
        return f"walk_matvec_{self.name}.npz"


def _per_group_sigma(sn, ng: int = 2) -> np.ndarray:
    """Per-group-varying σ_t (the reciprocity-file recipe) — keeps the group
    axis non-degenerate so a group-swap relocation defect moves values."""
    return np.stack(
        [np.full(sn.spatial_shape, 0.5 * (1.0 + 0.5 * g)) for g in range(ng)],
        axis=0,
    )


def _build(mesh_builder, seed: int):
    sn, _ = mesh_builder(ng=2)
    sig_t = _per_group_sigma(sn)
    lc = StreamingOperator(sn) + MultiplicationOperator.from_mesh(sig_t, sn)
    rng = np.random.default_rng(seed)
    psi = _random_composite(sn, rng)   # forward-matvec input
    phi = _random_composite(sn, rng)   # adjoint-matvec input
    return sn, lc, psi, phi


CASES: tuple[WalkMatvecCase, ...] = (
    WalkMatvecCase("slab_2g", lambda: _build(_make_slab, 20260705)),
    WalkMatvecCase("sphere_2g", lambda: _build(_make_sphere, 20260706)),
    WalkMatvecCase("cyl_2g", lambda: _build(_make_cyl, 20260707)),
)


def run_case(case: WalkMatvecCase) -> dict[str, np.ndarray]:
    """Apply both orientations; return the four value blocks."""
    _, lc, psi, phi = case.build()
    fwd = lc.apply(psi)
    adj = lc.apply_transpose(phi)
    return {
        "fwd_bulk": np.asarray(fwd.bulk.values, dtype=np.float64),
        "fwd_trace": np.asarray(fwd.boundary.values, dtype=np.float64),
        "adj_bulk": np.asarray(adj.bulk.values, dtype=np.float64),
        "adj_trace": np.asarray(adj.boundary.values, dtype=np.float64),
    }


def main() -> None:
    for case in CASES:
        blocks = run_case(case)
        path = SNAPSHOT_DIR / case.snapshot_name
        np.savez(
            path,
            fwd_bulk=blocks["fwd_bulk"],
            fwd_trace=blocks["fwd_trace"],
            adj_bulk=blocks["adj_bulk"],
            adj_trace=blocks["adj_trace"],
        )
        print(f"captured {path.name}: " + ", ".join(
            f"{k}{v.shape}" for k, v in blocks.items()
        ))


if __name__ == "__main__":
    main()
