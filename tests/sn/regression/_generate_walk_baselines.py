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

from orpheus.sn.loss_representation import FullFieldWavefront
from orpheus.sn.operators.streaming import StreamingOperator
from orpheus.transport.operators.multiplication_operator import MultiplicationOperator
from tests.sn._test_helpers import (
    cart2d_2g_nonsquare,
    random_radial_characteristic_field,
)
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
    """One (2-block composite, ψ½ composite) draw — the SAME rng order (bulk →
    boundary → seed) as the pre-B.2d 3-block builder, so the frozen
    baselines reproduce bit-identically; the ψ½ composite rides the walk's
    explicit legs (``None`` on non-carrying meshes). Phase C 4e retired the
    unified ψ½ leaf: the seed draw is the split-native composite whose
    per-slot values are bit-faithful to the retired unified single-buffer
    fill (see :func:`~tests.sn._test_helpers.random_radial_characteristic_field`)."""
    composite = _random_composite(sn, rng)
    seed_leaf = random_radial_characteristic_field(sn, rng)
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


class _OracleLC:
    r"""Rep-layer ``(L+C)`` pair for the cart2d row (#310 C3, spec §5.4).

    The multi-D adjoint DID NOT EXIST pre-carve, so this row is a NEW-value
    re-baseline (L17/2.5b), captured POST-carve after the
    ``test_multi_d_reverse_walk`` object oracles verified it.  Both
    orientations are captured at the REP layer — the full-cochain ORACLE
    arm (``FullFieldWavefront.loss_action`` / ``loss_action_transpose``)
    — because that oracle, not any deferral, is the stable structural
    reference: the production ``lc.apply_transpose`` went live at the
    #310 C4 flip (and the windowed PRODUCTION reverse is pinned
    bit-identical against this same oracle object), so the baseline is
    the value canary for the production path through the rep layer.
    """

    def __init__(self, sn, sig_t: np.ndarray) -> None:
        self._rep = FullFieldWavefront(sn)
        self._sig = sig_t

    def apply(self, psi):
        return self._rep.loss_action(self._sig, psi)

    def apply_transpose(self, phi):
        return self._rep.loss_action_transpose(self._sig, phi)


def _build_cart2d(seed: int):
    sn = cart2d_2g_nonsquare()               # rectangular nx≠ny (spec §5.5)
    sig_t = _per_group_sigma(sn)
    lc = _OracleLC(sn, sig_t)
    rng = np.random.default_rng(seed)
    psi, psi_seed = _draw_with_seed(sn, rng)   # seed is None (non-carrying)
    phi, phi_seed = _draw_with_seed(sn, rng)
    return sn, lc, (psi, psi_seed), (phi, phi_seed)


CASES: tuple[WalkMatvecCase, ...] = (
    WalkMatvecCase("slab_2g", lambda: _build(_make_slab, 20260705)),
    WalkMatvecCase("sphere_2g", lambda: _build(_make_sphere, 20260706)),
    WalkMatvecCase("cyl_2g", lambda: _build(_make_cyl, 20260707)),
    WalkMatvecCase("cart2d_2g", lambda: _build_cart2d(20260724)),
)


def run_case(case: WalkMatvecCase) -> dict[str, np.ndarray]:
    """Apply both orientations; return the four value blocks."""
    sn, lc, (psi, psi_seed), (phi, phi_seed) = case.build()
    if psi_seed is None:
        fwd = lc.apply(psi)
        adj = lc.apply_transpose(phi)
    else:
        # step 6: the joint matvec rows are THE GRID's block actions
        # (presence structural — the walk carries only the decoupled
        # block); the pinned blocks stay the System-A bulk/trace views.
        from orpheus.numerics.coupled_system import CoupledField

        from tests.sn._test_helpers import joint_m_grid

        grid, _space = joint_m_grid(sn, lc)
        fwd = grid.apply(CoupledField(systems=(psi, psi_seed))).systems[0]
        adj = grid.apply_transpose(
            CoupledField(systems=(phi, phi_seed)),
        ).systems[0]
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
