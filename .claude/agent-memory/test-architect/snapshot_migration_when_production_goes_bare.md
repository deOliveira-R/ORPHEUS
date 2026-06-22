---
name: snapshot-migration-when-production-goes-bare
description: Reusable recipe for migrating a per-sweep/per-step regression-snapshot harness when production code goes BARE (drops an intermediate, moves coupling external, changes per-step output while preserving the converged fixed point). Shared-driver SoT, schema=persisted∩compared, vacuum-bit-id correctness gate, snapshot-inheritance-needs-anchor, false-@catches retirement, term-activation re-verify.
metadata:
  type: feedback
---

When a refactor makes a production path BARE — it drops an intermediate, moves
a coupling step external (e.g. intra-sweep reflective BC → external
`reflect_outflow_into_inflow` called once per source iteration), or changes the
PER-STEP output while preserving the converged FIXED POINT — a per-step
regression-snapshot harness built on the old per-step output breaks. The
canonical instance (Wave-O #208 O.4b: `_sweep_2d_wavefront` went bare;
intra-sweep Gauss-Seidel → inter-sweep Jacobi, same fixed point, slower rate,
different per-sweep output). This recipe is reusable for any such migration.

**1. Mirror production in BOTH driver and generator via ONE shared helper
(coding-elegance Pattern 2).** Factor the production-mirroring logic (the
`/sum_w` per-ordinate projection, the reflect-inject + sweep loop) into the
TEST module; the generator IMPORTS it. Generator-vs-test drift becomes
unspellable. Inject the now-external coupling UNCONDITIONALLY (a no-op for the
vacuum case since B=0; verify bit-exact no-op before relying on it).

**2. Snapshot schema = exactly what is PERSISTED AND COMPARED — no more.** The
old schema often stores a richer intermediate (full interior-edge buffers) that
the test only ever SLICED (the boundary edges). When the bare path makes the
intermediate ephemeral/unrecoverable, the new schema stores ONLY the surviving
surface (the four boundary face views). The `.npz` files show as `M` in git
(schema flip) even where vacuum VALUES are bit-identical — expected, not a bug.

**3. The migration's correctness gate = VACUUM BIT-IDENTITY.** Before
regenerating, capture the legacy vacuum compared-arrays to a temp dir. After
regen, assert the vacuum cases reproduce them within nulp (got 0.0 / ~1e-16).
This proves bare ≡ legacy for ZERO inflow — the check that the bare path didn't
change the vacuum path. If vacuum DIFFERS → STOP, the bare refactor changed the
vacuum path = a BUG, do not commit. (The non-vacuum cases LEGITIMATELY change
per-step output — they are gated differently, see 5.)

**4. Snapshot inheritance from NEW (bare) code REQUIRES a structurally-independent
anchor (vv §1).** The reflective/non-vacuum baselines now inherit from the bare
code, so old-vs-new bit-id is not available for them. Add a closed-form anchor
case: all-reflective + uniform Q + uniform Σ_t → φ = (diag Σ_t − Σ_s^T)^{-1} Q
via `numpy.linalg.solve` (NOT a snapshot). That anchor makes the inherited
reflective snapshots correctness evidence rather than self-reference. State the
inheritance source (legacy-bit-id for vacuum, anchor-grounded for reflective)
explicitly in each case's docstring.

**5. A `@catches("ERR-NNN")` that no longer holds is a FALSE coverage claim —
remove it, don't keep it for history.** If the bare path no longer exercises
the code path the ERR lived in (e.g. no intra-sweep BC → ERR-003 intra-sweep
ordering is gone from this harness), DELETE the marker here and confirm the ERR
stays caught elsewhere (the 2G-convergence test still gates the converged
eigenvalue the external-reflect path reaches). A stale `@catches` is worse than
no marker — it asserts coverage that isn't there.

**6. Re-verify term ACTIVATION after adding any inject (vv Mode-7 discipline).**
An inject that silently NULLS an intended term is a silent coverage loss.
Confirm post-migration that the reflect inject does not zero the cases'
intended terms: multigroup asymmetric group ratios stay non-trivial, the
aniso-source branch still shifts the scalar flux, reflective faces stay nonzero.

Bare + external-reflect converges SLOWER (the migration's rate cost), so the SI
driver needs a deeper inner budget (`max_inner`, lower `inner_tol`) to reach
the same `rtol` — the converged VALUE is invariant; only the count grows.

See [[regression-tolerance-design]] (nulp vs SAFETY×conv_tol).
