---
name: Trajectory Resolvent hindsight architectural review
description: 2026-05-03 hindsight review of trajectory_resolvent (6 geometries × 2 topologies on shared variant_alpha_core). Identifies B1-B2-B3 as top refactor candidates ranked by leverage × low-risk. Plan deliverable; gated on test-architect's cross-method regression net.
type: project
---

# Trajectory Resolvent — hindsight architectural review (review-only dispatch)

**Date**: 2026-05-03
**Branch**: `main` (review only — no code change in this dispatch)
**Plan deliverable**: `.claude/plans/trajectory_resolvent_hindsight_refactor.md`

## Why dispatched

User asked: "look at the trajectory resolvent with the clarity of
hindsight and what are shared concepts that were not unified, and
elegant architecture that has not been implemented or leveraged."
Companion to a parallel `cross-domain-attacker` dispatch (frame-
level perspective) and a parallel `test-architect` dispatch (the
cross-method regression net that gates the implementation phase).

## Inventory verdict

`trajectory_resolvent/` is **8 modules / 4992 LOC** organised as:

- `variant_alpha_core.py` (436 LOC): 4 closure primitives — already
  unified (Class A in the plan).
- 6 geometry modules (sphere/cylinder/slab/slab_asym/hollow_sphere/
  annulus) sharing 12 structurally-identical power-iteration loops,
  12 result dataclasses, and 5 chord/trajectory primitive helpers
  scattered across files.

The boundary between "core" and "per-geometry" lives at the closure.
Phase-2 unification (2026-05-02) shipped that boundary. After 4 more
geometries (Phase-3A through 3C-2), the per-geometry overhead has
grown to ~6× the core's contribution — time to push the boundary
outward.

## Top-3 refactor candidates (ranked leverage × low-risk)

1. **B2 — `power_iterate_variant_alpha` driver** (HIGH leverage,
   LOW risk, S cost). 12 byte-identical outer loops collapse to
   one. Bit-equality preserved by IEEE-754 reproducibility.

2. **B3 — `GreensResult` + `PhaseSpaceMesh` tagged union** (MED
   leverage, LOW risk, S cost). 12 dataclasses → 1 + 6. New
   `k_history`, `residual_history` fields fall out for free.

3. **B1 — `ChordOracle` Protocol + per-geometry implementations**
   (HIGH leverage, MEDIUM risk, M cost). 6 chord-arithmetic blobs
   collapse to one Protocol + 6 concrete implementations. Cross-
   pollinates with peierls_nystrom + MoC + CP. Promotes to
   `orpheus.derivations.common.chord_oracle` after stabilising
   inside trajectory_resolvent.

These three are the "Variant α 2.0 architecture" the cross-domain-
attacker memo flagged. They land in order R1→R2→R3 *after* the
test-architect's cross-method regression net is alive.

## Class C (premature abstraction, do NOT unify)

- C1: phase-space dimensionality is NOT parametrisable — sphere's
  `(n_r, n_mu)` and cylinder's `(n_r, n_mu_axial, n_phi_az)` reflect
  different residual symmetry groups (SO(3) vs SO(2)×R), a load-
  bearing physical distinction.
- C2: rank-N closure (1≤N≤arbitrary) is NOT yet justified — only
  rank-1 and rank-2 have analytical references.
- C3: sphere MR segmentation and hollow-sphere b-partition look
  similar but are mathematically distinct (smooth-σ_t segmentation
  vs trajectory-class change at inner surface). B1 handles them
  with different concrete oracles, not unification.

## Cross-pollination (post-stabilisation)

- `ChordOracle` belongs in `common/`, not in trajectory_resolvent
  (Nyström + MoC + CP all consume it). Promote after R3 stabilises.
- `chord_half_lengths` in `common/kernels.py` is the proto-
  ChordOracle; absorb when promoting.
- `numerics.power_iteration` already exists; consume or upgrade
  after B2 is validated.
- Resolvent T = (I−S)^{-1} is structurally an asset across pillars
  (fn_method has its own `power_iterate_dominant_eigenmode`).
  Unified `numerics.dominant_eigenmode` driver is a "wait" item.

## Cardinal Rule 1 anchors (regression nets per refactor)

1. Bit-equality preservation (IEEE-754 floor) at every test of
   every geometry.
2. Branch-1 SymPy modules in `origins/specular/*.py` (algebra-of-
   record) are UNTOUCHED — every `derive_*()` foundation test must
   pass byte-equal before/after each refactor.
3. Cross-method regression net (test-architect deliverable) MUST
   be alive in CI before R1 begins. Without that net, a refactor
   bug compensating an existing bug would be invisible.
4. Failure-mode-induced bugs catalogued per refactor + their
   catching tests (B1 → variable swap, sign flip, missing factor;
   A4 → ERR-002 transpose convention; B4 → Signature 4 quadrature
   constant; etc.).

## Implementation order (gated on test-architect deliverable)

- R0: wait for cross-method regression net to land
- R1: B2 (driver) — 4-6 hours
- R2: B3 (result type + mesh) — half day
- R3: B1 (ChordOracle Protocol) — 2-3 days
- R4: promote ChordOracle to `common/` — half day
- R5: opportunistic polish (B4, B5, A4, A5)

## Convergence with cross-domain-attacker

The plan converges with the cross-domain-attacker memo
(`variant_alpha_family_hindsight.md` from 2026-05-02). Same top-3
refactors; same Class B vs Class C boundaries; same "wait for
instance N+1" treatment of full G-structure abstraction. The two
memos read together as frame-level (their memo) + code-level (this
memo) views of the same architectural diagnosis.

## Status

Review-only dispatch. Plan deliverable shipped at
`.claude/plans/trajectory_resolvent_hindsight_refactor.md`. No code
change. Implementation gated on the test-architect's cross-method
regression net. No archivist dispatch (no Sphinx narrative is
owed — the plan is internal architectural document, not theory).
