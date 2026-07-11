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


def _draw_with_seed(sn, rng):
    """One (2-block composite, ψ½ leaf) draw — the SAME rng order (bulk →
    boundary → seed) as the pre-B.2d 3-block builder, so the frozen
    baselines reproduce bit-identically; the leaf rides the walk's explicit
    legs (``None`` on non-carrying meshes)."""
    from orpheus.transport.fields.radial_characteristic_flux import (
        RadialCharacteristicFlux,
    )

    composite = _random_composite(sn, rng)
    seed_leaf = None
    if sn.radial_characteristic_space is not None:
        seed_leaf = RadialCharacteristicFlux.zeros_on(sn)
        seed_leaf.values[...] = rng.standard_normal(seed_leaf.values.shape)
    return composite, seed_leaf


def _build(mesh_builder, seed: int):
    sn, _ = mesh_builder(ng=2)
    sig_t = _per_group_sigma(sn)
    lc = StreamingOperator(sn) + MultiplicationOperator.from_mesh(sig_t, sn)
    rng = np.random.default_rng(seed)
    # The random ψ½ legs freeze the seed-fed matvec rows too (#282 route
    # (a) → B.2d explicit legs); the object-level anchor pins them
    # independent of the metric.
    psi, psi_seed = _draw_with_seed(sn, rng)   # forward-matvec input
    phi, phi_seed = _draw_with_seed(sn, rng)   # adjoint-matvec input
    return sn, lc, (psi, psi_seed), (phi, phi_seed)


CASES: tuple[WalkMatvecCase, ...] = (
    WalkMatvecCase("slab_2g", lambda: _build(_make_slab, 20260705)),
    WalkMatvecCase("sphere_2g", lambda: _build(_make_sphere, 20260706)),
    WalkMatvecCase("cyl_2g", lambda: _build(_make_cyl, 20260707)),
)


def run_case(case: WalkMatvecCase) -> dict[str, np.ndarray]:
    """Apply both orientations; return the four value blocks."""
    sn, lc, (psi, psi_seed), (phi, phi_seed) = case.build()
    if psi_seed is None:
        fwd = lc.apply(psi)
        adj = lc.apply_transpose(phi)
    else:
        from orpheus.transport.source_sinks import (
            RadialCharacteristicSourceSink,
        )

        fwd = lc.apply(
            psi,
            radial_characteristic_flux=psi_seed,
            radial_characteristic_source=(
                RadialCharacteristicSourceSink.zeros_on(sn)
            ),
        )
        # phi's ψ½ leaf rides the role-erased ``seed_cot`` leg (exactly the
        # slot the pre-eviction 3-block adjoint consumed).
        adj = lc.apply_transpose(
            phi,
            seed_cot=phi_seed,
            seed_cot_out=RadialCharacteristicSourceSink.zeros_on(sn),
        )
    return {
        "fwd_bulk": np.asarray(fwd.interior.values, dtype=np.float64),
        "fwd_trace": np.asarray(fwd.boundary.values, dtype=np.float64),
        "adj_bulk": np.asarray(adj.interior.values, dtype=np.float64),
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
