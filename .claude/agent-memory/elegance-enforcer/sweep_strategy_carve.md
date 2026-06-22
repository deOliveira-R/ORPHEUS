---
name: sweep-strategy-carve
description: The SweepStrategy carve (C3.4/C3.5 re-scoped) — S1 sweep-side polymorphic dispatch verdict + the locked design's review lens for S2-S5.
metadata:
  type: project
---

The SweepStrategy carve (branch `worktree-sn-nd-layout`, plan
`.claude/plans/sn_sweep_strategy.md`, designed 2026-06-10, LOCKED) replaces the
scattered procedural sweep/matvec dispatch with a first-class polymorphic
`SweepStrategy`. Continuation of [[nd-generic-carve-tuple-axis]] (C3.0-C3.3).

**The locked architecture (do NOT re-litigate — review faithfulness only):**
- Protocol `SweepStrategy(sweep, residual[S2], supports)` + frozen-dataclass base
  `_SweepStrategy` (mesh field + `__post_init__` guard raising `IncompatibleStrategy`)
  + family base `_DAGWavefront` (shared Cartesian `supports`).
- Three leaves: `CumprodScan`→`_sweep_1d_unified` (1-D any geom), `MovingFrontierWindow`
  →`_sweep_2d_wavefront` (prod opt), `FullFieldWavefront`→`_sweep_2d_full_field` (oracle).
- Selection SSOT: `supports(mesh)->Compatibility(ok,reason)` per strategy + priority-ordered
  `SWEEP_STRATEGIES` registry + `default_for(mesh)` (first applicable). Construction guard
  reuses `supports`. New `SNMesh.is_cartesian` (= `curvature is None`) is the geometry SSOT,
  ORTHOGONAL to `is_1d` (= `ndim==1`); selection keys on BOTH axes.
- Governing principle: construct general (capability), select narrow (policy in
  supports/default_for, NOT in code), specialize only on MEASURED cost.

**S1 VERDICT (sweep-side only, working tree): PASS-with-nits.** Reviewed
`sweep_strategy.py` (NEW), `geometry.py:643 is_cartesian`, `sweep.py:104 transport_sweep`
rewire, `test_unified_sweep_dispatch.py` migration.

Durable rulings (reinforce these were RIGHT — don't re-flag):
1. **Protocol + `_SweepStrategy` ABC is NOT a Pattern-2 twin** — legitimate
   interface(runtime_checkable Protocol)/shared-impl(dataclass guard) split. Python's
   only way to share the guard across 3 leaves while keeping the public contract a Protocol.
2. **`_DAGWavefront` is NOT hollow** — hosts the shared `supports` (both wavefront leaves
   report identical 2-D-Cartesian capability). Without it the predicate duplicates → Pattern-7
   drift habitat at S3 widening. Tie broken by REGISTRY ORDER (window before oracle), not by
   supports — correct capability-vs-policy split (locked decision 5).
3. **`supports` declaring the per-phase HONEST `ndim==2` (not eventual `ndim>=2`) is RIGHT,
   not a smell.** Advertising d-general now = Pattern-4 lie in the other direction (3-D mesh
   → shape-crash instead of clean `IncompatibleStrategy`). Widening triggers plan-tracked
   (S3 oracle, S4 window) + documented at declaration site = anti-pattern-11 satisfied.
4. **Moment-guard relocation `transport_sweep`→`CumprodScan.sweep` is RIGHT** (Pattern 4:
   the type that can't produce moments carries its own refusal). VERIFIED a MOVE not a
   duplication (git diff: old `if reduced is not None: if moment_projection` block deleted).
   `FullFieldWavefront`'s moment refusal is a DISTINCT reason (oracle-not-implemented vs
   1-D-algorithmically-impossible) — not a twin.
5. **Wrappers are byte-faithful** — each narrows the uniform Protocol kwargs to its backend's
   actual surface (`CumprodScan` forwards `initial_guess` not `moment_projection`;
   `MovingFrontierWindow` the reverse; `FullFieldWavefront` neither). Verified vs wrapped sigs.

Two CONCERN nits (optional, inline-fixable, NON-blocking — approved to commit):
- The lazy `from .sweep_strategy import default_for` in `transport_sweep` (correct cycle-break;
  consumer-side comment good) lacks a RECIPROCAL comment at `sweep_strategy.py:117`'s
  `from .sweep import ...`. Reciprocal-cross-reference discipline.
- `default_for -> SweepStrategy` returns `cls(mesh)` from `type[_SweepStrategy]` registry;
  "every registered entry is a CONCRETE leaf (not the abstract base whose `sweep` raises
  NotImplementedError)" is convention-only. A mis-registered base fails at call-time-deep, not
  registration. Comment at `SWEEP_STRATEGIES:375` suffices; `@abstractmethod` only if it
  provably doesn't perturb the bit-identity gate (frozen-dataclass/ClassVar trap caution).

**S2 VERDICT (matvec side, working tree): PASS.** Reviewed `operator.py` diff (apply →
`sweep_strategy.residual`, apply_transpose → `residual_transpose`, `_apply_1d`/`_apply_1d_transpose`
verbatim extraction, `sweep_strategy` cached_property) + the +96 `sweep_strategy.py` S2 additions
(`residual`/`residual_transpose` on Protocol + base + 3 leaves + shared `_DAGWavefront` deferral).

Durable rulings (reinforce these were RIGHT — don't re-flag):
1. **NO twin-delivery-plumbing smell — the predicted #1 hazard did NOT materialize.** The matvec is
   single-sourced at the operator (`apply` has ONE return: `strategy.residual(self, psi)`; no second
   fold). Each twin provably consumes ONE discretization: `MovingFrontierWindow.sweep`→
   `sweep_graphs[OctantLabel(σ)].apply_windowed`, `.residual`→`...residual_windowed` (forward/matvec
   faces of ONE `SweepDependencyGraph`, same DD kernel); `FullFieldWavefront` likewise (op comment
   "SAME cell kernel as residual_windowed"); `CumprodScan` shares the 1-D `_compute_LpC` body. Pinned
   by the existing `window≡full` oracle. This is the GOOD shape my inst-knowledge #1 names: routes
   single-sourced at the operator AND consuming the same graph object.
2. **`strategy.residual(operator, psi)` reaching `operator._apply_2d_cartesian` is the RIGHT seam,
   NOT a layering smell.** Strategy = the matvec-VARIANT selector; the matvec BODIES stay
   operator-owned (they need σ_t, M_spatial, sweep_graphs — all operator state). Passing the whole
   `operator` (not destructured pieces) is correct: the strategy's job is to pick WHICH operator
   method, not to re-plumb its args. Weighed against A2D-1: relocating `_apply_2d_cartesian` into the
   strategy would change its source-hash + widen blast radius — WRAP-not-relocate confirmed (diff
   touches only docstring/dispatch lines around `_apply_2d_cartesian`, never its body).
3. **`_apply_1d`/`_apply_1d_transpose` extraction is bit-identical + cleanly homed.** Verbatim move
   from apply/apply_transpose's 1-D branches (only `sn_mesh`→`self.sn_mesh` rename + import relocation
   into the new methods). Geometry asymmetry (1-D in `_apply_1d`, 2-D in `_apply_2d_cartesian`, both
   operator methods selected by the strategy) is HONEST — they ARE different algorithms; symmetry is
   restored at the strategy-selection layer, not forced into one body. Symmetry-in-math = symmetry-in-
   code satisfied.
4. **Shared `_DAGWavefront.residual_transpose` deferral is RIGHT single-sourcing.** ONE
   `NotImplementedError` raise inherited by both 2-D leaves (the adjoint is deferred for the whole
   DAG family, not per-buffer-policy). `_DAGWavefront` now carries real shared content (supports +
   residual_transpose) — confirmed NOT hollow: neither leaf redefines either; leaves override only
   sweep+residual (the genuinely-differing parts).
5. **cached_property `sweep_strategy` mirrors M_spatial/M_angular_redist correctly** (selection fixed
   by mesh, stable for the operator's life). Lazy-import cycle-break is consistent + acceptable:
   `sweep_strategy.py` imports `sweep` at module load; both `transport_sweep` and the cached_property
   reach `default_for` via function-local imports — no new cycle.

DEVIATION RULINGS (both JUSTIFIED — agree with the implementer):
- **D1: keep the 3 `_MSpatial` raise-guards** (not the plan's "remove all 4"). CORRECT. L20 audit
  verified: `_compute_decomposition` has 3 non-strategy callers (transient orchestrator `op.py:264`,
  `M_spatial.apply:1123`, `M_angular_redist.apply:1241`) → its guard is load-bearing. These are
  RAISE-not-route guards (Pattern-4 illegal-state refusal at internal helpers), NOT dispatch points —
  keeping them does not reintroduce scattered dispatch. The plan's "remove 4" was written before the
  caller audit; the evidence overrides the letter. The ONLY removed raise was the live `apply`/
  `apply_transpose` 2-D dispatch branch — exactly the scattered dispatch the carve targets.
- **D2: 2-D adjoint keeps `NotImplementedError` (not `IncompatibleStrategy`).** CORRECT. The mesh IS
  compatible (forward matvec works); only the adjoint FEATURE is deferred. `IncompatibleStrategy` is
  a strategy↔mesh-incompatibility signal — the wrong concept here. Pattern-4 (right exception names
  the right illegal state): a deferred-feature is not an illegal-pairing. The raise message names the
  deferral (O.2b lands 1-D reverse sweep first) — honest, never a silent wrong answer.

NO nits this round (the 2 S1 optional nits — reciprocal lazy-import comment, registry concrete-leaf
invariant — are S1-scope, untouched by S2; not re-raised). S2 is a clean PASS.

**S3 VERDICT (d-generic spine + cumprod≡spine equivalence + adapter retirement): clean PASS, ZERO
nits.** First NON-bit-identical phase. S3a committed `37ce528`; S3b (equivalence-test rewrite) the only
working-tree change. Reviewed `geometry.py:660 streaming(axis)`, `operator.py:1791 _apply_full_field`,
`sweep.py:1149 _sweep_full_field`, `sweep_strategy.py:336 supports-override`, both test files.

Durable rulings (reinforce — don't re-flag):
1. **In-place generalization left NO 2-D residue.** `ng=shape[0]; spatial=shape[1:]`; `*spatial` buffers;
   `"ng...,n->g..."` ellipsis reduction; `str_x`/`str_y` GONE from both spine bodies (survive only in the
   d=2 WINDOW paths `sweep.py:1046`/`operator.py:1697` = S4 scope). One special-case branch was
   GENERALIZED AWAY (`sx/sy/if sx==0 and sy==0` → `signs=label[:ndim]; if not any(signs)`), none added.
2. **The `label[:ndim]` spelling-divergence is HONEST, not Pattern-7 drift.** Sweep uses
   `sweep.label.signs[:ndim]` (`.signs` extracts tuple from an `OctantLabel`, `sweep_schedule.py:76`);
   matvec uses `octant.label[:ndim]` (`quad.octants`→`DiscreteMeasurePartition.label` is a PLAIN TUPLE,
   `measure.py:684`). Different source TYPES (OctantLabel vs tuple) → different spelling, same in-plane
   tuple out. Projecting the S² phantom out-of-plane sign is principled (`OctantLabel` docstring
   `sweep_graph.py:88` allows phantom sign_z the in-plane sweep is invariant under).
3. **Pure-z guard `not any(signs)` provably dormant at d=1** — a real 1-D octant is `(±1,)`, `any` always
   truthy → short-circuit structurally unreachable; the `# (d≥2 only)` comment is verified-correct.
4. **`streaming(axis)` pull-forward from S5 = JUSTIFIED, not scope-creep.** Governing-principle "construct
   general": the spine GENUINELY needs the d-uniform accessor NOW (`tuple(streaming(a) for a in
   range(ndim))` is the line that retires the hand-listed `(streaming_x,streaming_y)` order-drift habitat
   the prior code's OWN comment flagged). Two real consumers at landing (sweep+matvec spine) → clears
   Pattern-6. Body `(self.streaming_x,self.streaming_y)[axis]` is int-keyed (not stringly), 2 boundary
   guards (Cartesian `AttributeError` + `0<=axis<ndim` `IndexError`) with domain-precise messages.
   Correctly does NOT yet retire the d=2 WINDOW's hand-listed tuples (that's S4) — pays its way at landing,
   seeds the S4 retirement without reaching into it.
5. **`FullFieldWavefront.supports` widen to any-d Cartesian (drop `ndim==2`) = HONEST, not over-claim.**
   The HONEST INVERSE of the S1 window ruling: there, advertising `ndim>=2` for a d=2-only body would be a
   Pattern-4 lie; HERE the spine BODY is genuinely d-general (DAG walk + `OctantLabel` 3-tuples), so
   advertising `is_cartesian` MATCHES the impl. d=3 "over-claim" is inert (no 3-D mesh constructor exists);
   the only future gap (`streaming_z`) fails LOUD at the accessor guard, never silent. Factoring is right:
   base keeps d=2 for the window (its body IS d=2-only), spine overrides because it ALONE widened — minimal
   honest delta, NOT a predicate that should live elsewhere.
6. **Two full-field bodies (`_sweep_full_field` oracle vs `_sweep_2d_wavefront` window) = the accepted
   `feedback_aggressive_retirement` fuller-view-oracle exception, DOUBLY pinned now (inst-knowledge #3).**
   Shared `DiamondDifference.cell_kernel_batch` (math can't drift, only storage). Pin #1 d=2 `window≡full`
   `np.array_equal` (angular+scalar+TRACE, unchanged). Pin #2 NEW d=1 `CumprodScan≡FullFieldWavefront` at
   nulp — the oracle is now ALSO the d-generic correctness spine the cumprod opt is validated against.
   Corner-probing satisfied (not just ℓ=0 scalar). Oracle exception, correctly kept.
7. **Equivalence-test rewrite = cleanest possible form + COMPLETE retirement.** 5 adapter scaffolds
   (`_wavefront_1d_sweep` stub, `_cumprod_1d_sweep`, `_SpineNotLanded`, `_spine_landed`,
   `test_diamond_difference_importable` keepalive) + 2 `xfail(not _spine_landed())` guards RETIRED; the
   equivalence is now 2 lines through the real strategy API. Grep `orpheus/ tests/` = ZERO source refs to
   any retired name (only `docs/_build/html` generated artefacts, out of scope, regen on sphinx-build).
   This is `feedback_retirement_means_test_migration` done right (rewired, not delete-only).
8. **`assert_array_almost_equal_nulp` @ `_NULP_BOUND=128` = the RIGHT principled-equiv gate (NOT
   bit-id).** Cumprod parallel-prefix vs spine sequential-forward-sub = same affine recurrence, different
   FP association → ULP drift principled. Bound DERIVED (`depth≈nx≤16 × ×8 safety = 128`) with the
   tighten/re-justify contract in-comment (anti-pattern-17 satisfied). Contrast the d=2 oracle's
   `np.array_equal` (same per-cell eval order there) — right gate chosen per each pair's FP structure.
   Mode-9 config adequate: 2g asymmetric-downscatter + het per-cell Σ_t + non-uniform source + non-zero
   inflow (psi_in≠psi_out everywhere; H1/H2 dodged) + the independent k_inf=1.875 transfer-matrix anchor.

**S4 VERDICT (widen window to `frontier_dim=d−1`): PASS-with-nits.** The core is `sweep_graph.py`
(+584/−352): new mesh-time `_LevelFrontier`(per-level `read`/`write`/`seed`/`shed` selectors) +
`_FrontierPlan`(whole-sweep `spatial_shape`/`free_bbox`/`det`/`levels`/`is_point`) + d-tuple-of-slabs
`_MovingFrontier`(`seed`/`incoming`/`emit`/`shed`); `_build_frontier_plan` REPLACES the d=2-only
`_window_metadata`; `apply_windowed`/`residual_windowed` now d-agnostic walks on per-axis TUPLES
(`inflow`/`str_axes_octant`/`capture`); graph field `window_slots|None`+`window_edges|None` → one
non-optional `window_plan: _FrontierPlan`. Call sites (`sweep.py:938`, `operator.py:1750`) tuple-packed
byte-faithful (x=0,y=1 order kept); A2D-1 hash regenerated w/ honest WRAP-not-relocate changelog.

Durable rulings (reinforce — RIGHT):
1. **Mesh-time-plan factoring is the CORRECT abstraction, not over-structured.** All d-dependent index
   arithmetic in `_build_frontier_plan` (L16 zero per-sweep recompute); walk d-agnostic. The 3-class
   split = data-the-walk-reads (`_FrontierPlan`, frozen, mesh-time) vs mutable-state-advanced (`_MovingFrontier._win`,
   per-sweep) vs one-level-addressing (`_LevelFrontier`, named `read`/`write`/`seed`/`shed` = Pattern-3 vs
   the anonymous-positional-tuple smell). Keep all three — collapsing reintroduces the C3 anon-tuple smell.
2. **NO twin duplication beyond the irreducible solve/apply kernel.** `apply_windowed`/`residual_windowed`
   share the same `_MovingFrontier` + `window_plan` + `_FrontierPlan` levels + `seed→incoming→kernel→emit→shed`
   skeleton; differ ONLY in `cell_kernel_batch` vs `residual_kernel_batch` + output reduction. Pinned by
   `window≡full` oracle now at d=1/d=2/d=3.
3. **`box_contiguous = det <= 1` (`:567`) is clean + load-bearing.** ONE predicate ("free region is interval
   ⟺ d≤2") → `slice` (d=2 bit-id + contiguity view) vs fancy-array (d≥3 simplex) at 2 mesh-time sites;
   the walk stays uniform. Verified `slice(c[0], c[-1]+1)` reproduces retired `win_x[...,i0:i1+1]` exactly.
4. **Complete retirement.** `window_slots`/`window_edges`/`_window_metadata`/`_bounds`/`seed_x`/`seed_y` +
   the d=2-only RuntimeError guards GONE (grep: zero live refs); `|None` correctly dropped (plan exists d≥1).
   Equivalence test rewritten d-generically off the retired `inflow_x=` kwargs.
5. **Naming reads like the cochain trace algebra** — `seed`(ι_*)/`incoming`(gather)/`emit`(advance)/`shed`(ι*),
   each docstring states the correspondence. Master standard PASS.

VERIFIED BYTE-FOR-BYTE vs legacy (the high-stakes d=2 bit-id check): x-seed → ghost col 0 from `inflow[0][jj@i0]`
(matches `seed_x(prev, inflow_x[:,:,jj[0]])`); y-seed → slot `i1`=max-local-i where local_j==0 from `inflow[1][ii@j0]`
(matches `seed_y(prev, i1, inflow_y[:,:,ii[-1]])`); the determined-axis `write.append(read)` is read-only slice reuse
(no alias risk). d=1 path traced end-to-end: `det=0`, `free_bbox=()`, slab `(N_oct,ng,2)`, `read=()`, `is_point`
re-inserts/squeezes the kernel's `n_diag=1` — principled (genuine `frontier_dim==0` geometry).

Two CONCERN nits (optional, NON-blocking, approved to commit):
- **`is_point` triple-branch** (`incoming:290`/`emit:302`/`shed:319`) — anti-pattern-7 singular-point branch,
  but MILD: uniform (one root cause = kernel's fixed `n_diag` contract vs the d=1 no-free-axis slab),
  verified by the d=1 admission, and the alternative (phantom length-1 free axis) is a LESS-honest data
  structure. RECOMMEND LEAVE AS-IS (option a); a `_as_diag`/`_from_diag` helper (option b) is a 3-line collapse
  not worth the indirection. Do NOT pursue the phantom-axis route.
- **`det==d-1` / `free_bbox = shape[:det]` coupling** (`:565`) — `is_point`(count-based) and `det`(prefix-slice)
  are two spellings of the free/determined partition; consistent TODAY (STRENGTH, `det` is the SSOT). Latent IFF
  S5+ ever makes determined-axis a POLICY (non-last) → then `free_bbox` must become the non-`a==det` comprehension,
  not a prefix slice. Optional one-line comment; no action now.

Carve scorecard: S1 PASS-with-nits, S2 PASS, S3 clean PASS, S4 PASS-with-nits. NEXT = S5 (retire the
2-D orchestration: `OctantLabel.sign_x/sign_y/streams_in_2d` shims + `str_x`/`str_y` hand-list in
`sweep_octant_group` + `_sweep_2d_wavefront` → d-generic; `MovingFrontierWindow.supports` widens off `ndim==2`).

**S5.1a VERDICT (ScanMarch strategy, issue #222, working tree): PASS-with-nits.** New
`_scanmarch_row`+`_sweep_2d_scanmarch` (`sweep.py:1186/1248`), `ScanMarch` leaf + 3rd registry slot
(`sweep_strategy.py:484/589`), `test_scan_march_equivalence.py` (G2 oracle). OPT-IN (Fork B1:
registered AFTER window → `default_for` unchanged, A2D-1/affine anchors stay free-green).

**THE DUPLICATION VERDICT — accept (a), Fork-B1-transitional, NOT factor now.** `_sweep_2d_scanmarch`
re-emits the octant scaffold from `sweep_octant_group` (oct_idx/sx/sy projection, pure-z `Q/σ_t`
branch, sx_eff/sy_eff μ=0 ride-up, x/y in/out face strings, inflow read, `face_view[oct_idx]=capture`
shed). This is a REAL Pattern-2 octant-projection twin (verified line-by-line: `:1318-1358` ≡
`:874-924`), the ONLY difference being the interior solve (anti-diagonal `apply_windowed` vs row-march).
WHY (a) is right NOW, not weasel-deferral: (1) `sweep_octant_group` is MONOLITHIC — scaffold + the
single `apply_windowed` call are inlined in ONE loop body with NO interior-solve callback seam; the
factor-now option (b) requires cutting that seam INTO the production window (`_for_each_inplane_octant`
yielding projected octant + a solve callback), which is exactly the bit-identity-anchor surgery Fork B1
forbids until the default flips. (2) The twin is BOUNDED (2 instances, the 3rd never coming — d≥3 march
recurses INSIDE `_sweep_2d_scanmarch`, reusing the SAME scaffold, not adding a 3rd). (3) The collapse
destination is NAMED and plan-tracked (S5.3 default-flip + window retire/unify). HAZARD if it lingers:
the day a 3rd octant-projection consumer or a projection-convention edit (new face, metric-weighted
shed) lands on one-not-the-other. ⚠ DOC GAP (the nit): `_sweep_2d_scanmarch`'s docstring SAYS "Mirrors
`_sweep_2d_wavefront`'s octant projection + boundary-trace I/O exactly; only the interior schedule
differs" — TRUE but does NOT name the duplication as a tracked transitional cost nor point at S5.3. The
mirror-claim reads as a courtesy note, not a Pattern-2 IOU. REQUIRED before commit: one sentence naming
the `sweep_octant_group` scaffold twin + the S5.3 consolidation trigger (anti-pattern-11: a transitional
duplication needs a written removal trigger, not just a "mirrors X" aside). Factoring-now (b) REJECTED:
trades a bounded+tracked+anchor-safe duplication for unbounded bit-identity-anchor risk on the
production path — wrong call under Fork B1.

Durable rulings (reinforce — RIGHT):
1. **`_scanmarch_row` reverse/moveaxis is CLEAN, not procedural.** `ordinate_scan` scans axis-0; the
   octant batch carries x on the last axis → `moveaxis(-1,0)`/`moveaxis(0,-1)` is the honest adapter
   (named `out_x_sweep`/`in_x_sweep`/`psi_avg_sweep` in SWEEP order, then `[...,::-1]` back to MESH
   order for −x). The `x_reverse` reflect-coefficients-scan-reflect-back is the standard chain-reversal,
   no cleaner closed form. ψ̄=½(in+out) + out_y=2ψ̄−ψy_in = the WDD diamond read verbatim — master
   standard PASS (reads like "scan x, close the diamond, thread y").
2. **`_sweep_2d_scanmarch` reads like the algorithm.** α=2s_x/D−1, β=2(Q+s_y·ψy_in)/D, D=σ_t+s_x+s_y
   all named with the docstring math; `capture_x`/`psi_y_in`/`out_y`/`x_outflow` carry domain meaning;
   the y-row march direction `range(ny)` vs `range(ny-1,-1,-1)` is honest per-octant. Scan-along-x +
   march-over-y is legible.
3. **The deferral raises are HONEST (anti-pattern-11 + Pattern-4).** `ScanMarch.residual` d=2 raise
   names S5.1b + the A2D-1-free-green reason; `residual_transpose` d=2 names the O.2b ordering. Mesh is
   compatible, only the FEATURE pends → `NotImplementedError` not `IncompatibleStrategy` (matches the
   D2 ruling). 1-D residual/residual_transpose correctly route to `operator._apply_1d[_transpose]`.
4. **Moment einsum is single-sourced at the contraction.** Scan-march per-row `"nlm,ngi,n->lmgi"` is
   the per-y-row slice of the window's per-diagonal `"nlm,ngd,n->lmgd"` (`sweep_graph.py:926`) — same
   `Y·ψ_avg·w` projection, same `moment_buf[...] +=` accumulation; `(moment_buf, None)` return mirrors
   `_sweep_2d_scheduled` (no aliasing). Principled-equiv (cross-row `+=` reorder), pinned by the G2
   `ScanMarch≡MovingFrontierWindow` moment test (ℓ≥1 anti-degeneracy guard live under -O via pytest.fail).
5. **`CumprodScan.sweep`/`ScanMarch.sweep` 1-D moment-refusal is a DELIBERATE twin, ACCEPTABLE.** Both
   raise the identical L21 message ("moment output is 2-D Cartesian only…"). This is a 2-instance
   verbatim duplication of a refusal string — TODAY acceptable (the refusal is the SAME illegal-state,
   and ScanMarch's 1-D arm genuinely forwards to `_sweep_1d_unified`, the SAME body CumprodScan wraps).
   Borderline Pattern-7 (one convention, two spellings); if a 3rd 1-D-moment-refusing site appears,
   collapse to a shared `_reject_1d_moment_output()` helper. NOT blocking — flag for S5.3.

S5.1a nit (REQUIRED before commit, inline-fixable): the `_sweep_2d_scanmarch` docstring must name the
`sweep_octant_group` scaffold as a tracked Fork-B1 transitional twin + the S5.3 consolidation trigger
(upgrade the "Mirrors … exactly" aside into a Pattern-2 IOU). Optional (S5.3): collapse the 1-D
moment-refusal twin string.

**S5.1b VERDICT (row-march MATVEC twin, issue #222, working tree): PASS — clean, ZERO blocking nits.**
New `StreamingOperator._apply_2d_cartesian_scanmarch` (`operator.py:1791`), the Pattern-2 extraction
`_x_scan_faces` (`sweep.py:1186`) with `_scanmarch_row` REFACTORED to call it (`sweep.py:1269`),
`ScanMarch.residual` 2-D wired (`sweep_strategy.py:549`), G2.c `test_scanmarch_residual_equals_oracle`.

**`_x_scan_faces`-EXTRACTION VERDICT: clean SSOT, correctly factored (NOT over/under).** This IS the
ONE genuinely-shared piece — the x-face recurrence solve (reverse→`ordinate_scan`→`in_x=concat(seed,
out[:-1])`→un-reverse) — between the SWEEP (α=attenuation, β=affine-source) and the MATVEC (α=−1,
β=2ψ̄). The two consumers differ ONLY in (a) the coefficients passed and (b) what they do with the
returned faces: sweep computes `0.5*(in+out)` + `out_y=2ψ̄−ψy_in`; matvec uses `in_x` for the
`−s_x·in_x` residual term and re-derives `out_y` from the KNOWN probe. Return triple `(in_x, out_x,
x_outflow)` = the UNION of both consumers' needs (`out_x` load-bearing for the sweep, `_`-discarded by
the matvec — the honest idiom; trimming would force two return-shape variants, a worse smell). The
refactor is BIT-IDENTICAL for the sweep: old computed `0.5*(in+out)` in sweep-order then reversed,
new returns reversed faces then computes `0.5*(in+out)` in mesh-order — reversal is an element
permutation, `0.5*(a+b)` elementwise → `reverse(½(in+out))≡½(reverse(in)+reverse(out))` bit-for-bit;
`x_outflow=out_x_sweep[...,-1]` unreversed, unchanged. G2 sweep tests staying green confirms. Master
standard PASS (docstring states the α/β correspondence for BOTH directions).

**MATVEC-TWIN VERDICT: SAME accept-don't-factor as the sweep twin (S5.1a). The matvec case does NOT
differ.** `_apply_2d_cartesian_scanmarch` duplicates `_apply_2d_cartesian`'s octant-projection block
(`:1723-1735`), face strings (`:1740-1743`), `LpC−sig_t·probe` subtraction (`:1761-2`), AND the ENTIRE
O.4b boundary-residual block (`:1772-1782`) — VERIFIED byte-identical line-by-line, the ONLY difference
being the interior solve (`graph.residual_windowed` anti-diagonal vs the y-row march loop). Identical
Fork-B1 situation: `_apply_2d_cartesian` routes through `graph.residual_windowed` (the A2D-1-hashed
production window matvec) with NO interior-solve callback seam — factoring a `_for_each_inplane_octant`
helper would re-carve that hashed path, which Fork B1 forbids until S5.3. Twin BOUNDED (2 instances;
d≥3 recurses inside, reusing the scaffold). **The IOU note IS adequate** — the `.. note::` (`:1822-1834`)
names the Pattern-2 duplication, the Fork-B1/A2D-1 reason, "**Edit both in lockstep.**", AND the S5.3
consolidation trigger; it explicitly cross-references the sweep sibling's same-named note. This is the
S5.1a nit (which I required there) applied PROACTIVELY here — the matvec docstring did NOT need the
correction the sweep did. anti-pattern-11 satisfied (transitional duplication WITH a written removal
trigger). Accept (a); factor-now (b) REJECTED (same reasoning).

Durable rulings (reinforce — RIGHT):
1. **Reads-like-the-math PASS.** `α=−1, β=2ψ̄` (reflection-scan reuse) is LEGIBLE not cryptic — docstring
   states "since ψ̄ is KNOWN the closure out_x=2ψ̄−in_x IS a first-order recurrence in the faces (a
   pure-reflection scan)"; `alpha_reflect = np.full(..., -1.0)` is self-naming, `2.0*psi_bar_row` is the
   β at the call site. The body reads "reconstruct x-faces from the probe → evaluate the (Σ_t+s_x+s_y)ψ̄
   −s_x·in_x−s_y·in_y residual → march y (out_y=2ψ̄−ψy_in threaded)" — the WDD diamond verbatim.
2. **Named intermediates domain-meaningful.** `LpC`(=(L+C)ψ̄ accumulator, docstring'd), `D_row`(=σ_t+s_x
   +s_y, the cell denominator), `in_x_row`/`x_outflow`/`cap_x`(domain x-outflow per y-row)/`psi_y_in`/
   `out_y`/`streamed`(per-face shed dict) all carry physical/cochain meaning. `probe`=ψ̄ (consistent w/
   the sibling). NO anonymous reductions.
3. **`ScanMarch.residual` 2-D wire is HONEST + symmetric.** 1-D→`_apply_1d` (geometry-blind), 2-D→the
   new method; the raise replaced by the real dispatch; docstring states L21 + the α=−1 row-march +
   the G2.c oracle pin. `residual_transpose` d=2 STILL the O.2b deferral raise (unchanged — correct,
   the adjoint is genuinely unwired). The strategy now row-marches in BOTH directions = symmetry-in-
   math→symmetry-in-code (the sweep/matvec faces of ONE ScanMarch).
4. **G2.c test = correct conventions.** Reuses `_build_mesh` + the principled-equiv `_RTOL`/`_ATOL`
   (NOT nulp — the docstring at `:54-63` justifies why: near-zero shed elements amplify ULP). Pins
   BOTH bulk AND the O.4b boundary-block residual (the latter named the face-shed-convention-drift
   catcher — exactly the corner the role grid demands). Fresh deterministic seed offset (`13`). Types
   the probe as `AngularFlux`/`BoundaryFlux` (a FLUX is the apply INPUT) → operator mints
   `AngularSourceSink`/`BoundarySourceSink` (Aψ is a SOURCE/SINK) — role grid respected (inst-knowledge
   #4). The matvec is single-sourced at the operator (`ScanMarch.residual` has ONE return) — NO twin-
   delivery-plumbing smell (the inst-knowledge #1 hazard did NOT materialize, same as S2).

NO blocking nits. One OPTIONAL S5.3-deferred observation (NOT raised as a nit): `_apply_2d_cartesian`
uses `sn_mesh.sweep_graphs[OctantLabel((sx_eff,sy_eff))]` for its interior; the scanmarch matvec uses
the row loop — at S5.3 when the window retires, the two interiors collapse to one scan-march body and
the scaffold twin dissolves with it (the named destination). No action now.

Carve scorecard: S1 PASS-with-nits, S2 PASS, S3 clean PASS, S4 PASS-with-nits, S5.1a PASS-with-nits,
S5.1b PASS (zero nits).

**S6.3 VERDICT (operator/representation re-layering — move the walk OFF the operator, issue #222,
working tree): PASS-with-nits.** The four `_apply_*` (`_apply_1d`/`_apply_2d_cartesian`/
`_apply_2d_cartesian_scanmarch`/`_apply_full_field` + the two `_transpose`) DELETED from
`StreamingOperator`; bodies MOVED into the leaves' `loss_action` (CumprodScan/ScanMarch 1-D →
`M_spatial._compute_LpC`; MovingFrontierWindow → the windowed `residual_windowed` walk;
FullFieldWavefront → the WavefrontFlux full-cochain walk). Reviewed `operator.py` (−512 net),
`loss_representation.py` (+415), `solver.py`+test retargets, the NEW `test_loss_action_convention.py`.

Durable rulings (reinforce — RIGHT):
1. **The `−C` glue is single-sourced now, exactly 2×.** It was 5× duplicated across the deleted
   `_apply_*`; now lives in EXACTLY `operator.py:1469` (`apply`) + `:1519` (`apply_transpose`) — one
   per direction, the irreducible minimum (forward + adjoint are genuinely different walks). VERIFIED
   by grep: `- self.sigma_t[None] *` appears at precisely those 2 lines; ZERO σ_t SUBTRACTIONS leaked
   into the moved `loss_action` bodies (the 3 `LpC[oct_idx] = sig_t * probe` are pure-z ASSIGNMENTS
   setting (L+C)ψ̄=Σ_tψ̄ for no-in-plane-streaming ordinates — correct, not a double-count). No leak.
2. **`apply` reads like `L = (L+C) − C` character-for-character.** `lpc = loss_action(...)`;
   `out_bulk = lpc − σ_t·ψ`. Named intermediate `lpc`=(L+C)ψ; the comment states the algebra +
   "applied ONCE here (it was 5× duplicated)". Master standard PASS.
3. **The sweep/apply asymmetry is intentional + documented.** `sweep → (L+C)⁻¹q` DIRECTLY (no `+C`
   add-back — grep: ZERO `+ sig_t`/`+ self.sigma_t` in any sweep body), because the inverse is of the
   FULL composite (`InvertibleOperator.solve` docstring: "Invert (L+C)ψ = rhs"); `apply` is on the
   leaf `L` so it subtracts `−C`. Clean.
4. **Fork-B1 deferral to S6.4 is plan-sanctioned + correctly fenced.** MovingFrontierWindow.loss_action
   + ScanMarch.loss_action share the octant frame + O.4b boundary block (differ only in interior walk:
   `residual_windowed` anti-diag vs `_x_scan_faces` row-march). BOTH carry `.. note::` with
   `**Edit both in lockstep.**` + consolidation trigger S6.4/`_OctantWalk2D` (re-homed from the
   superseded S5.3 — grep: ZERO stale S5.3 refs). Bounded (2 instances), reciprocal, tracked. Pulling
   forward to S6.4 = right (plan stages S6.3=move/S6.4=extract, each bit-id-gated; factoring INTO this
   commit would re-carve the A2D-1-hashed window mid-move). anti-pattern-11 satisfied.
5. **A2D-1 source-hash pin correctly REGENERATED, not source-frozen.** Retargeted
   `inspect.getsource(StreamingOperator._apply_2d_cartesian)` → `MovingFrontierWindow.loss_action`;
   new SHA recorded; changelog documents the relocation with OUTPUT-byte-identity proven by the
   window≡full MATVEC oracle + affine-carve golden (NOT source identity). The A2D-1 "regenerate +
   assert output via oracle" discipline done right. Test green.
6. **Retirement complete + test-migrated.** `grep "def _apply_" operator.py` = ZERO; ZERO live refs in
   `orpheus/`; the tests/ hits are comments/docstrings except the A2D-1 pin (retargeted) +
   `test_two_d_cartesian_loss_action_returns_result` (renamed off `..._routes_through_apply_2d_cartesian`
   — the contract rename reflects the new home). retirement-means-test-migration satisfied.
7. **The NEW `test_loss_action_convention.py` is exemplary — the architectural defense.** Pins
   `loss_action=(L+C)ψ` NON-tautologically: (a) flat-reflective `L·ψ_flat=0` anchors loss_action=σ_t·ψ
   to the independent k_inf analytical fact; (b) the `−C` glue cross-checked against a SEPARATELY-built
   `CollisionOperator` (subtrahend verified =σ_t·ψ, not the production self-check); (c)
   `loss_action − apply == C·ψ` catches a leaf that returns bare `L·ψ`. `-O`-safe, foundation, ≥2G het.
   This is exactly the test that would catch the S6.6 `ExplicitMatrix` author who returns `M_stream@ψ`.

THE NIT (CONCERN, inline-fixable, NON-blocking — should fix before commit). **The Protocol +
`_LossRepresentation` base `loss_action` docstrings are STALE — they still state the PRE-S6.3 `L·ψ`
convention while the convention flipped to `(L+C)ψ`.** `loss_representation.py:208` (Protocol method):
`r"""The forward matvec twin L ψ…this applies L"""` — WRONG post-S6.3. Worse, the SAME file's
class-bullet at `:188` says `(L+C)ψ` — self-contradiction at the contract surface. ALL FOUR concrete
leaves were correctly updated to "FULL within-group loss `(L+C)ψ`"; only the Protocol (the canonical
contract a new-representation author reads FIRST) + base were left behind. Bug habitat: the plan
EXPLICITLY anticipates `ExplicitMatrix.loss_action` at S6.6 (`loss_action=M@ψ`); an author consulting
the Protocol docstring reads "return L·ψ", returns the streaming-only matrix, and the operator's
`apply` subtracts C a SECOND time → double-counted collision. The convention pin guards EXISTING
leaves but a NEW leaf is exactly where the stale contract bites before its pin is written. Required:
update `:188`/`:208-215` (Protocol) + the base `:276` to "(L+C)ψ; operator subtracts C".

**S6.4(a) VERDICT (`_OctantWalk` extraction — collapse the Fork-B1 octant-frame twin, #222, working
tree): PASS — ZERO blocking nits.** THE deferred-debt payment named at S5.1a/S5.1b/S6.3. New
`_OctantWalk` (`loss_representation.py:328-535`: `_ApplyOperands` apply-direction problem-data bundle +
`_inflow_faces`/`_outflow_faces` d-generic helpers over `_AXIS_NAMES` + `_interior_walk` shared frame +
`loss_action` apply frame); `MovingFrontierWindow.loss_action`/`ScanMarch.loss_action` COLLAPSED to
3-line delegators (`return _OctantWalk(self.mesh).loss_action(operator, psi, self._loss_action_interior)`),
each supplying only its `_loss_action_interior` kernel; both `.. note:: Edit both in lockstep` IOUs
DELETED. A2D-1 source-hash pin (`TestT4dApply2DCartesianSourceHashPin`) RETIRED. New
`tests/sn/operators/test_one_octant_walk.py` (AST anti-flag tripwire + both-variants SPY + strict-xfail
sweep SPY). Focused gates green (2 passed, 1 strict-xfail), confirmed by me.

Durable rulings (reinforce — RIGHT):
1. **The Fork-B1 twin is GENUINELY ELIMINATED, not relocated.** Octant projection + pure-z + face strings
   + inflow read + outflow shed + O.4b boundary block — written 2× through S6.3 — now ONE place. Both IOUs
   deleted. Pinned at runtime by `test_both_matvec_variants_share_the_walk` (passes from (a)). This is the
   NAMED collapse destination of inst-knowledge #1/#2; the deferred elegance debt was PAID. (Pattern 2.)
2. **The two-LAYER kernel signature is the load-bearing extensibility win.** `_interior_walk`-level
   `interior` returns CAPTURE ONLY (direction-neutral); the `loss_action`-level kernel returns
   `(LpC_oct, capture)`; the `run_interior` CLOSURE adapts (side-effects `LpC[oct_idx]=LpC_oct` into the
   closed-over accumulator, returns capture). ⟹ sub-step (b) supplies a DIFFERENT closure (solve emits
   angular/moment internally, returns capture) and `_interior_walk` does NOT change. "matvec≡sweep is one
   walk" becomes a forthcoming CODE fact; strict-xfail SPY records the (b) gap. (Pattern 5 + 1.)
3. **Closures (`pure_z`/`run_interior`/`shed`) are EMIT-POLICY collaborators, NOT procedural lambdas.**
   They capture per-call-fresh accumulators (`LpC`/`streamed`); a Protocol/class per policy would be
   Pattern-6 over-abstraction (1-3-line adapters, no reusable 2nd consumer). Direction carried by
   OBJECTS, never a boolean `is_solve` — pinned MECHANICALLY by the AST identifier tripwire
   `test_octant_walk_is_kernel_parameterized_not_boolean` (inspects `ast.Name`/`Attribute`/`arg`, so
   docstrings naming the anti-pattern don't false-trip). Anti-pattern 3 prevented by construction + tested.
4. **`FullFieldWavefront.loss_action` is the DOCUMENTED (d) holdout, verified not accidental.** Still has
   its own inline `for octant in quad.octants` + pure-z + boundary block (`:844-877`); walks `graph.residual`
   with the typed `WavefrontFlux` whole cochain (genuinely different storage → different capture shape →
   correctly deferred to (d)). `_OctantWalk` docstring + gate-memo staging name (d). Removal trigger tracked.
5. **A2D-1 retirement is `feedback_retirement_means_test_migration` done right.** Source-hash on a now-SHARED
   walk trips on every legit refactor with no behavior signal → retired; tripwire job → the
   `assert_array_equal` `window≡full` MATVEC output oracle vs the structurally-distinct full-field oracle
   (NON-SQUARE cases catch x↔y swap = the inst-knowledge corner). Test class deleted, `MovingFrontierWindow`
   import dropped, successor note in oracle docstring, `operator_algebra.rst` Key-Facts + body annotated
   with full S6.3→S6.4(a) provenance. Complete.
6. **In-plane projection SINGLE-SOURCED at `_octant_sweep` (`sweep_schedule.py:173`).** Matvec
   (`_interior_walk` via `SweepSchedule.jacobi(...).groups[0].sweeps`) + sweep-side `sweep_octant_group`
   consume the SAME `OctantSweep` objects → projection cannot drift sweep↔matvec. The legacy inline
   `sx=label[0]; sy=label[1] if len>=2 else 0` retired. The "two-spellings-of-one-partition" tell did NOT fire.

Inst-knowledge tells checked, NONE fired: capture allocs are fresh (`np.empty_like(inflow)`, not views) →
no aliased-slot trap; no SNMesh rebuild-`.copy()` seam (`_OctantWalk(self.mesh)` carries the built mesh);
pure-z/grazing rationale comments state correct invariants. The S6.3 stale-contract nit I raised is FIXED
(Protocol `:210-226` + base `:285-295` now state `(L+C)ψ`) — the (b)/(d) author reads a clean contract.

Carve scorecard: S1 PASS-with-nits, S2 PASS, S3 clean PASS, S4 PASS-with-nits, S5.1a PASS-with-nits,
S5.1b PASS, S6.3 PASS-with-nits (stale-contract-docstring), S6.4(a) PASS (zero nits — the deferred-twin payment).

**S6.4(b) VERDICT (the SWEEP frames come into the shared walk, #222, working tree): PASS — 3 stale-text
CONCERNs, ARCHITECTURE clean, ZERO structural nits.** The deferred-debt payment finished on the SOLVE side:
`sweep_octant_group` (162 lines, was in `__all__`) + `_sweep_2d_scanmarch` (the "edit both in lockstep" IOU)
GENUINELY DELETED, bodies relocated as `MovingFrontierWindow._sweep_interior` / `ScanMarch._sweep_interior`
behind `_OctantWalk.sweep_group` (the solve mirror of `.loss_action`). `_sweep_2d_scheduled` now
kernel-parameterized (REQUIRED `interior` kwarg). The strict-xfail SPY FLIPPED to PASS — L21 (matvec≡sweep)
is a CODE FACT. Confirmed `3 passed` under -O by me.

Durable rulings (reinforce — RIGHT):
1. **`_SweepEmit` guarded type is a real Pattern-4 win, NOT ceremony.** Predecessor took 4 loose optional
   kwargs → a half-wired output (`moment_buf` given, `Y=None`) was REPRESENTABLE + einsum-crashed. The
   `__post_init__` `if angular==moment: raise` makes the angular-XOR-moment illegal state unconstructible.
   The `if emit.moment_buf is None` branch IS the established Phase-5c output-mode idiom (mirrors
   `apply_windowed`'s own DI), NOT the direction-flag anti-pattern: direction carried by OBJECTS
   (`_sweep_interior`/`_loss_action_interior`, `_SolveOperands`/`_ApplyOperands`, `_SweepEmit` vs `LpC`/`streamed`
   closures); `moment_buf is None` tests WHICH OUTPUT BUFFER. AST tripwire confirms no boolean direction leaked.
2. **`_sweep_2d_wavefront`'s `interior=None` window default is single-sourced at the KERNEL BODY, not a twin.**
   `_sweep_interior` is ONE method; `MovingFrontierWindow`/`_DAGWavefront` are frozen mesh-only dataclasses →
   a fresh-default instance ≡ the bound `self` byte-for-byte (kernel reads only `self.mesh`). Default is justified
   (the historical bare-window entry; 2 direct test consumers `test_2d_octant_sweep_equivalence:390` /
   `test_2d_full_field_oracle:90` depend on it). Same shape S1's `default_for→cls(mesh)` approved. Latent habitat
   (future per-instance field the default wouldn't carry) inert today.
3. **`solver.py` G-S `default_for(sn_mesh)._sweep_interior` cross-module private reach = transitional plumbing
   with a NAMED trigger (inst-knowledge #1).** Selection single-sourced through `default_for` (same SSOT as
   `transport_sweep` + the `interior=None` default); comment names S6.5 (operator's own instance) + S6.4(f)
   (cycle dissolution). Kernel single-sourced. CONCERN-not-VIOLATION, consistent w/ S6.4(a)/S2.
4. **Apply/solve symmetry is clean; the 2 asymmetries are PRINCIPLED (inst-knowledge #6).** (a) apply frame
   whole-Jacobi (`.loss_action` builds its own group) vs solve per-group (`sweep_group`, schedule-driven) —
   only solve has a G-S consumer; per-group apply granularity = Pattern-6 speculation correctly avoided. (b)
   apply uses `LpC`/`streamed` closures not `_SweepEmit` — matvec has ONE output mode, solve has the genuine
   angular-XOR-moment bifurcation. Math IS asymmetric → code matching it is correct.
5. **Retirement complete (the deletions are real, not relocated-twins).** grep: zero LIVE `sweep_octant_group`/
   `_sweep_2d_scanmarch` calls; surviving refs are past-tense notes EXCEPT the 3 CONCERNs below. `__all__`
   trimmed. The DD import dropped from sweep.py.

THE 3 CONCERNs (stale text the deletion left behind; NOT architecture; fix before commit):
- **Required (Cardinal Rule 3, Sphinx-breaking)**: `:func:` xrefs to DELETED symbols —
  `tests/sn/sweep/cartesian_2d/test_scan_march_equivalence.py:3` (`_sweep_2d_scanmarch` → repoint
  `ScanMarch._sweep_interior`) + `docs/theory/operator_algebra.rst:5072` (`sweep_octant_group` → `_OctantWalk.sweep_group`).
- **Should-fix prose**: `operator_algebra.rst:4795` "`sweep_octant_group` loops octants" → name `_OctantWalk._interior_walk`
  (the cross-octant `+=` moment call graph that page documents).
- **Should-fix message**: `test_one_octant_walk.py:110-114` live `pytest.fail` string still says "sweep still uses
  `sweep_octant_group`'s private frame … S6.4(b) is OPEN" — misstates current state (a misstated-invariant bug
  habitat). The docstrings were past-tensed; this string was missed.

Carve scorecard: …, S6.4(a) PASS (zero nits), S6.4(b) PASS (3 stale-text CONCERNs, architecture clean).
Pattern for the (c)-(f) reviews: DELETION sub-steps leave dangling `:func:` xrefs + descriptive prose + live
test-failure-message strings naming the deleted symbol; grep `orpheus/ tests/ docs/` for the deleted name BEFORE
verdict and triage each hit (past-tense-note OK vs live-xref/prose/message = CONCERN).
