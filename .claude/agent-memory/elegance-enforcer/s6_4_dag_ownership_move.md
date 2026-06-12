---
name: s6-4-dag-ownership-move
description: S6.4(c) DAG-ownership move + S6.4(d) full-field-oracle fold — both PASS; per-shape lru_cache, the _OctantWalk ×3→×1 consolidation, in-edge-only-seed≡ι_* proof, stale-WavefrontFlux-docstrings concern
metadata:
  type: project
---

S6.4(c) of #222 (worktree `sn-nd-layout`): the per-octant `SweepDependencyGraph`
family moved OFF the mesh (`SNMesh.sweep_graphs` attribute + the curvilinear `= None` illegal
slot, both scrubbed from `geometry.py`) ONTO the `_DAGWavefront` representation family via a
`@property sweep_graphs` that routes to a module-level `@lru_cache(maxsize=8)`
`_graphs_for_shape` behind the `SweepDependencyGraph.for_shape(spatial_shape)` classmethod.

**VERDICT: PASS** (4 CONCERN-level nits, none blocking; all are "documentation/test-shape,"
not bug-habitat). Textbook Pattern-4 illegal-states-unrepresentable win + a
many-meshes→one-family-per-shape SSOT tightening. The mesh is pure geometry again.

**Why each reviewed axis cleared:**
- **`lru_cache(maxsize=8)` bound**: NOT a fragility. Graphs are pure `(shape×octant)` derived
  objects; eviction only re-pays construction (docstring says exactly this). `is`-identity is
  asserted ONLY in the test (`test_same_shape_meshes_share_the_graphs`, a cache-hit witness) —
  production NEVER compares graph identity, it consumes graph CONTENT (`graph.residual`/`.apply`/
  `.levels`). Eviction can never produce a wrong answer, only a re-build. Bound of 8 is generous.
- **property-per-access** (`for_shape` called each octant loop): cache hit is a dict-lookup on a
  small tuple key — negligible vs the per-octant DAG walk. Caching on the representation INSTANCE
  would BREAK same-shape sharing (two reps of one shape → distinct dicts) — correctly rejected.
- **mutable-dict-shared-across-callers**: documented "treat as IMMUTABLE." `MappingProxyType`
  would be strictly stronger (illegal-states-unrepresentable for mutation) — CONCERN not
  VIOLATION: no caller mutates; the inner `SweepDependencyGraph` values are the real shared state
  (frozen dataclass). The standing follow-up if this reopens.
- **test-migration completeness**: COMPLETE. Content pins (4 octants, d=1 chain, hand-`_diag_cache`
  golden, type contract, lebedev-independence, no-rebuild) all carry through `for_shape`; the
  "absent on curvilinear" pins → structural `not hasattr(mesh,"sweep_graphs")` ×4 coord + the
  `geometry.py`-source-grep + the DAG-free-rep grep. `l0`→`foundation` correct (no theory label).
  Pre-existing bare-`assert` Mode-8 gap retired. 18 pass under -O.
- **`_sweep_full_field` free-function reaching `for_shape`** (the documented (d) holdout): PASS.
  Predecessor sweep-body that `FullFieldWavefront.sweep` delegates to; the matvec twin
  `loss_action` is an instance method using `self.sweep_graphs`. Both consume the SAME
  single-source `for_shape` cache → NOT a twin-delivery smell.

**Standing follow-up if S6.4(d) reopens this**: MappingProxyType the returned mapping
(enforcement > documented-immutability). Low priority — no live mutator.

## S6.4(d) — full-field ORACLE folds into the shared walk (reviewed uncommitted, 2026-06-11)

**VERDICT: PASS** with two doc/scope CONCERNS (none block the core carve). The FINAL
consolidation: octant frames ×3 → ×1, O.4b boundary blocks ×4 → ×1. `_sweep_full_field`
(the third sweep frame, the (c) holdout above) DELETED; `FullFieldWavefront` is now
kernels-only — `sweep` routes through `_sweep_2d_scheduled(SweepSchedule.jacobi,
reflect=None, interior=self._sweep_interior)`, `loss_action` through
`_OctantWalk.loss_action(..., self._loss_action_interior)`. Same frame as
`MovingFrontierWindow`/`ScanMarch`; the interior storage policy is the ONLY fork.

LOAD-BEARING ARGUMENTS VERIFIED (do NOT re-litigate):
- **In-edge-only seed ≡ historical whole-trace ι_* seed, BYTE-identical.** `_octant_face_cochain`
  zeros each per-axis buffer (`n_a+1` slots), seeds ONLY the domain in-edge from `inflow[a]`.
  Sound: `graph.apply`/`graph.residual` walk `levels` in upwind topological order — every
  interior + out-edge slot is an upstream out-face (diamond `ψ_out=2ψ̄−ψ_in`), written before
  read; in-edge is the ONLY read-but-never-walk-written slot. OLD `WavefrontFlux.seed` touched
  ONLY the edge slot too (interior already 0) — old `pf[oct_idx].copy()` was already zeros+in-edge.
  "Whole-trace ι_*" = whole-trace across FACES, edge-only spatially.
- **Per-octant `_edge_outflow` shed ≡ old `wavefront.absorb` whole-trace writeback** — pinned
  bit-for-bit by the oracle test's post-sweep boundary-trace assert (ran GREEN).
- in/out-edge indices: in=`0`(+a)/`spatial[a]`(−a), out mirror — consistent w/ `_inflow_faces`.
- **axis binding now positional over `range(ndim)` at 3 producers** (`_octant_face_cochain`,
  `_inflow_faces`, operands `str_axes`) — BETTER SSOT than old `WavefrontFlux.axes` indirection.
- d=1: `_sweep_2d_scheduled` made d-generic (`spatial=sig_t.shape[1:]` vs `ng,nx,ny`); d=1 spine
  pinned by `test_cumprod_1d_equals_full_field_spine` (nulp). 72 tests GREEN under -O.

THE TWO CONCERNS (carry to (f) review):
1. **Stale WavefrontFlux orchestrator-binding docstrings**: `sweep_graph.py:750-752`,
   `cell_update.py:289-293`, `sweep.py:854-855` still claim the per-axis tuple is "born from
   `WavefrontFlux.face(a)` over `WavefrontFlux.axes` at the orchestrator" — FALSE post-(d) (raw
   `range(ndim)` loop). The "rationale-comment misstating a load-bearing invariant" tell. Fix or
   (f)-defer-with-tracked-list.
2. **Bundled scope creep**: the diff un-xfails two `@slow` cylinder twin-path tests
   (`test_l1_standoff_slab_cylinder.py`), unrelated to the carve, justified by a prior D-K heal +
   an XPASS(strict) run I could NOT verify (slow-suite forbidden). Strict-xfail removal = contract
   change → own commit or cited verifying run.

**WavefrontFlux orphaning correctly flagged-not-deleted.** 5b deliberately un-orphaned it as the
typed fuller-view oracle (aggressive-retirement EXCEPTION); after (d) its boundary algebra
DISSOLVED into the frame → exception no longer applies (fuller view now = raw `_octant_face_cochain`
buffers, same `window≡full` pin). Delete is distinct (f) scope (latent consumer
`interior_face_space.py`). Deferral disciplined IFF (f) is a real tracked artifact listing:
WavefrontFlux delete + `_sweep_2d_scheduled` rename (name now lies — admits d=1) + the 3 stale docs.
