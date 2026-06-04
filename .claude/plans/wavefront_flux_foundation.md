# WavefrontFlux — typed interior face-flux cochain (Wave O / #205 / #208)

**Branch:** `refactor/field-role-typing` (worktree). **Base HEAD:** `266fcf5`
(B.5.2 landed). **Mode:** main-agent DIRECT authorship, turn-by-turn user
steering (`feedback_no_method_implementer_for_surgical_carves`) — touches the
SN sweep + matvec hot path; commit per phase; bit-identity gates throughout.
**Standing exclusions from commits:** `.claude/agent-memory/`, `.mcp.json`,
hooks, `derivations/diagnostics/`.

**Design basis (READ FIRST — the native-frame validation):**
- cross-domain-attacker: `.claude/agent-memory/cross-domain-attacker/field_role_typing_faceflux_frames.md`
  (the cochain frame; field+views NOT per-face objects; the `(octant×face)`
  G-S graph; the KBA/Pautz anchor; the honest "representation-only, no rate" disclaimer).
- explorer (sweep/cell-update structure): `.claude/agent-memory/explorer/sn_si_reflective_gauss_seidel_recovery_surface.md`.
- B.5.2 role-grid precedent: `.claude/plans/b52_boundary_residual_retype.md` (the
  type-only / bit-identical carve discipline). Wave O umbrella: `wave_o_operator_typing.md`.

---

## STATUS (2026-06-04) — Phase 0 PASSED; verified worktree state

**⚠ READ-THE-WORKTREE DISCIPLINE.** All `file:line` and `grep` MUST target the
worktree (`/Users/rodrigo/git/nuclear/ORPHEUS/.claude/worktrees/field-role-typing`),
NEVER the main checkout (`/Users/rodrigo/git/nuclear/ORPHEUS`, branch
`refactor/sn-operator-algebra` — a DIFFERENT, older branch whose `sweep.py`/`operator.py`
DIFFER from this branch). When in doubt, confirm the live source with
`inspect.getsource(fn)` under `PYTHONPATH=<worktree>`. (A wrong-tree read this
session falsely concluded the matvec was an FD stencil; Phase 0 caught it.)

**Verified worktree state (this branch, via `inspect.getsource`):**
- **O.4b is ALREADY COMPLETE here.** `_apply_2d_cartesian` routes through
  `graph.residual` → `residual_batch` (DD closure); `_compute_gradients` (the FD
  stencil) is **RETIRED**; `matvec ≡ sweep` is a valid gate. The matvec carries
  raw `psi_x`/`psi_y` interior buffers (`operator.py:1402-1407` seed, `:1429`
  `graph.residual`) → it IS a WavefrontFlux typing target (Phase 3 is REAL).
- **Both sweep and matvec are BARE-boundary** (O.4b Phase E done): seed the given
  inflow trace, no in-sweep `bc.apply`, reflection via the external sibling `−B`;
  the matvec emits an active-trace boundary RESIDUAL as `BoundarySourceSink`
  (outflow slot `streamed−given`, inflow slot `given`). The seed/absorb the
  WavefrontFlux `ι_*`/`ι*` types is unchanged by bareness.

**Phase 0 de-risk: PASSED** (`derivations/diagnostics/diag_phase0_wavefront_derisk.py`,
`python -O`, against the bare worktree sweep):
- (a) typed `ι_*`/`ι*` round-trip **bit-identical** to raw seed/absorb; `ι_*∘ι*=id`
  (biproduct law); zero-copy (`shares_memory`); negative control (transposed seed)
  breaks → gate non-vacuous.
- (b) full bare sweep with a flat-buffer `WavefrontFlux` backing (interior
  `FaceLayout` + reshape views) **bit-identical** (`max|Δ|=0.0`) on
  angular/scalar/boundary, ALL reflective + vacuum cases — AND the two-array
  backing equally bit-identical (localizer: no layout drift).
- (c) **no perf regression** (typed/raw median ratio 0.96–1.01, < +5% gate).
- §3a′ axis-parametric smoke: the interior `FaceLayout` builds for `axes=(0,)`,
  `(0,1)`, `(0,1,2)` through ONE path (`FaceLayout.from_named_shapes` is already
  axis-parametric — no new primitive needed).
→ **The flat-buffer `WavefrontFlux` substrate is SOUND. Proceed to Phase 1 mint.**

**Verification corrections (supersede §5 where they conflict):**
- **L16 gate = FOCUSED SUBSET, not full `tests/sn`** (which does NOT finish in
  20 min). Use the carve-path fast tier + the hot-path microbench
  (`diag_l16_wavefront_microbench.py`), with a re-captured **worktree-branch**
  baseline (the test-architect's baseline was measured on the MAIN checkout).
- **Pre-existing reds (worktree, 2026-06-04):** `test_solver_components.py` is RED
  (`NotImplementedError: R-1 Step E 2-D Cartesian SI deferred`) — exclude from the
  tier; cylinder #206/#196; #212 hang. `test_krylov_curvilinear_precond_safety`/
  `b1pp`/`restart` are GREEN here (memory was stale).

---

## 0. What this is (and the native frame)

Type the SN interior cell-face angular fluxes — currently raw ephemeral numpy
arrays `psi_x (N_oct,ng,nx+1,ny)` / `psi_y (N_oct,ng,nx,ny+1)` — as a named
field **`WavefrontFlux`**. This kills the Pattern-3 "unnamed quantity"
anti-pattern (flux-bearing tensors with no type identity) and **names the trace
operator the sweep already applies by hand** (`psi_x[:,:,0,:] =
boundary.face_view("xmin")` is `ι_*`; the writeback is `ι*`).

**Native frame = discrete exterior calculus / cochains** (attacker-validated):
- The per-ordinate angular flux on faces is a primal **1-cochain** `ψ⁽¹⁾_Ω ∈ C¹`
  (a value per oriented face).
- The cell-average ψ̄ is a **0-cochain**; the diamond closure
  `ψ_out = 2ψ̄ − ψ_in` is the averaging map `C¹ → C⁰`.
- The boundary trace is the **pullback `ι*`** under `∂Ω ↪ Ω`; "absorption =
  identity" is `ι_* ∘ ι* = id` on the boundary chain.
- The full face cochain splits as the **biproduct** `C¹ = C¹_int ⊕ C¹_∂`, which
  IS the #208 `V_bulk ⊕ V_boundary` shape: `WavefrontFlux = C¹_int` (interior),
  `BoundaryFlux = C¹_∂` (boundary), coupled by `ι*`/`ι_*` at the domain edges.

**Honest scope (load-bearing — do NOT oversell):** this is a
**representation / elegance win only**. It does NOT change asymptotic cost,
does NOT recover the SI rate, and does NOT enable parallelism on its own (the
attacker stated this as an absence, not a hedge). The seed/absorb stays an
**inherent cheap copy** (`O(boundary faces)`, negligible vs the `O(volume)`
sweep) at the **persistent-boundary / ephemeral-interior lifetime boundary** —
true zero-copy is precluded (and unnecessary) because `BoundaryFlux` persists
across SI iterations while `WavefrontFlux` is rebuilt each sweep. The payoff is
the **type** (named field, typed `ι*`/`ι_*`, reads-like-the-cochain-math,
illegal-states-unrepresentable), and it is the clean substrate the G-S schedule
(`si_gauss_seidel_recovery.md`) then lands on without rework.

---

## 1. The name + the concept hierarchy

- **`WavefrontFlux` (NEW)** — the interior face-flux field the wavefront sweep
  propagates: the flux across all (interior) faces traversed as the wavefront
  advances from the inflow trace to the outflow trace. The interior 1-cochain
  `C¹_int`. **Ephemeral, single-role (flux), vectorized field + views.** The
  name carries the operational signature (the propagating sweep state), NOT a
  static single face — that is why it supersedes the earlier "FaceFlux" working
  name. (The wavefront at topological level `k` (`i+j=k`) is a *slice* of it;
  see §3 storage decision.)
- **`BoundaryFlux` (EXISTS)** — the boundary face-flux trace `C¹_∂`:
  persistent, role-typed (Flux / Source / Residual — the role grid lives here
  because the BC consistency residual `ψ.inflow − B·ψ.outflow − q.inflow` lives
  here), inflow/outflow-partitioned (`TraceSpace`), `FaceLayout`-backed.
- **`ι_*` (seed)** : `BoundaryFlux.inflow → WavefrontFlux` domain-edge slots.
  **`ι*` (absorb)** : `WavefrontFlux` domain-edge slots → `BoundaryFlux.outflow`.
  Typed trace operators replacing the raw edge-slice copies. Clean copies (the
  orderings already agree — proven by the working `psi_x[:,:,0,:] =
  face_view(...)`), not zero-copy (lifetime split, §0).
- **NO standalone per-face `BoundaryFaceFlux` / `FaceFlux` object** — the
  attacker REJECTED per-face objects on three grounds (vectorization/L16; the
  cochain frame is storage-granularity-indifferent; biproduct consistency).
  Per-face access stays a **view** (`face_view`), as `BoundaryFlux` already does.

Role grid placement (the WHY): interior faces are the transient off-diagonal of
a triangular factor → no residual role → `WavefrontFlux` is flux-only. The role
grid is a 0-cochain (cell) concept that the boundary 1-chain inherits only via
its BC residual.

---

## 2. Current state (what exists vs. what is raw) — all `file:line` at HEAD `266fcf5`

- **Raw interior face fluxes (= untyped WavefrontFlux):** `sweep.py:853-854`
  (`psi_x = np.zeros((N,ng,nx+1,ny))  # ephemeral interior cache`; `psi_y`) and
  the matvec mirror `operator.py:1402-1403`.
- **Seed `ι_*`:** `sweep.py:857-860` (`psi_x[:,:,0,:] = boundary_flux.face_view("xmin")`,
  `…[:,:,nx,:]`=xmax, `psi_y[:,:,:,0]`=ymin, `…[:,:,:,ny]`=ymax). Matvec mirror in
  `_apply_2d_cartesian`.
- **Absorb `ι*`:** `sweep.py:940+` (`boundary_flux.face_view("xmin")[:] = psi_x[:,:,0,:]`, …).
- **The carrier:** `SweepCellSlice` (`cell_update.py:181`) holds `psi_x`/`psi_y`
  (`:295-296`) + `psi_avg_probe` (`:304`), mutated in place.
- **The walk:** `SweepDependencyGraph` (`sweep_graph.py:127`); `from_cartesian_2d`
  (`:191`, levels `i+j=k`, reversed per octant sign); shared `_make_slice`
  (`:264`); `apply` (`:300`, solve→`update_batch`) + `residual` (`:386`,
  matvec→`residual_batch`). DD closure: `diamond.py` `update_batch` (`:273`,
  scatter `ψ_out=2ψ̄−ψ_in` at `:354-355`), `residual_batch` (`:361`, `:424-425`).
- **The 2-D sweep:** `_sweep_2d_wavefront` (`sweep.py:756-940`): seed 4 faces
  up front (`:857-860`), `for octant in quad.octants` (`:852`), per-octant
  `graph.apply` (`:903-925`), absorb after loop (`:931-943`).
- **The 2-D matvec:** `_apply_2d_cartesian` (`operator.py:~1334-1477`): same
  shape, `graph.residual`, emits the boundary block (B.5.2: `BoundarySourceSink`).
- **The 1-D sweep:** `_run_1d_sweep` (`sweep.py:382-739`): edge fluxes + two
  direction groups (`:523`). **The 1-D matvec boundary:** `_compute_LpC`
  (`operator.py:~534`) + `_compute_decomposition` (`:~575`) (the twin).
- **`BoundaryFlux` (the model to mirror):** `transport/fields/boundary_flux.py`;
  base machinery in `transport/fields/_bases.py` (`from_mesh` `:426`, `zeros_on`
  `:457`, `face_view` `:398`, `layout` `:416`); `FaceLayout` in
  `numerics/face_layout.py:108` (`faces`, `total_size`, per-face slots).
- **`BoundaryFaceFlux`:** a RETIRED module (`spatial/boundary_face_flux.py`,
  subsumed by the Phase-C sweep-frame matvec; refs in
  `pole_angular_closure.py:78,151`). Do NOT resurrect it as an object.

---

## 3. The design

### 3a. Storage path — A FIRST, THEN B (user directive 2026-06-04)
`WavefrontFlux`'s name is the *frontier*; two backings express it, and we do
**A, secure the win, then go STRAIGHT to B** (not a deferred maybe):
- **(A) Full interior face-field — FIRST.** Per-axis dense arrays (the current
  `psi_x`/`psi_y`), the wavefront at level `k` is a slice. Preserves the proven
  `(N_oct,ng,n_diag)` vectorized walk EXACTLY (the hot path reads
  `WavefrontFlux` views, byte-identical indexing) → minimal risk, type-only,
  bit-identical. Memory `O(N·ng·nx·ny)` (unchanged). This isolates the *typing*
  carve from the *storage* carve so each lands clean.
- **(B) Moving-frontier window — STRAIGHT AFTER A's win.** Store only the active
  diagonal(s) (`O(N·ng·(nx+ny))` — a genuine memory win and literally "the flux
  at a wavefront"). Changes the walk to a rolling 2-diagonal buffer. Lands on
  the A-proven `WavefrontFlux` type (the API already fits), so B is a storage
  swap behind a stable type, gated by the same bit-identity (the *converged*
  solution is unchanged; only the live working-set shrinks). Its own phase, its
  own elegance review.
**Why split:** A makes the type bit-identical (free verification by inheritance);
B then changes only the backing behind that type. Doing both at once would
conflate a type change with a storage change — two bug habitats in one commit.

### 3a′. 3D-READINESS imperative (build the types axis-parametric NOW)
We have **1D + 2D = two working instances**, so `feedback_unify_after_two_instances`
*licenses* a dimension-generic foundation now (3D = the validating third
instance, the future "general foundation" session — see
`.claude/plans/nd_foundation.md`). The cochain frame has no `d` baked in, so the
NEW `WavefrontFlux` / interior `FaceLayout` / `ι*` **types + interfaces MUST be
built axis-parametric** (`feedback_architecture_forward_not_legacy_fit`: extend
to N-D, not N+1):
- A face field over **`axes`** (1→2→3 face orientations), NOT hardcoded
  `psi_x`/`psi_y` attributes — e.g. `WavefrontFlux.face(axis)` over a
  `tuple[axis] of dense arrays` or a single layout-indexed buffer.
- The interior `FaceLayout` keyed by `(axis, position)`, the boundary = the
  `2d` codim-1 edge faces (2 per axis) — the same generalization the trace
  `FaceLayout` needs.
- `ι*`/`ι_*` over the edge subset of EACH axis (not a 4-face literal).
**But** keep the *implementation* targeting 1D+2D (no `Mesh3D`, no 8-octant
code, no `i+j+k` walk yet — those land in the `nd_foundation` session). The
contract: nothing in the new types may *hardcode* "2 axes / 4 faces / 4
octants"; the 2D-specific call sites (`from_cartesian_2d`, `DiamondDifference`'s
explicit x/y, `_sweep_2d_wavefront`) stay 2D for now but are flagged in
`nd_foundation.md` as the generalization targets. **Test the axis-parametricity:
a unit test that the `WavefrontFlux`/`FaceLayout` API accepts `axes=(0,)` (1D)
and `axes=(0,1)` (2D) through the SAME code path** — proving 3D (`axes=(0,1,2)`)
is a parameter, not a fork.

### 3b. The type
- `WavefrontFlux` carries the interior face cochain: the x-normal field
  (`nx+1` × `ny` faces) and the y-normal field (`nx` × `ny+1` faces) per
  `(N_oct, ng)` — i.e. wraps `psi_x` + `psi_y` (2-D); just the x-faces (1-D).
  Backing = a flat buffer + an **interior `FaceLayout`** (mirror `BoundaryField`)
  OR the two dense arrays with an explicit layout descriptor (decide for
  zero-copy vectorized views). `UNITS = ANGULAR_FLUX_UNITS` (cm⁻²·s⁻¹·sr⁻¹ —
  same as `BoundaryFlux`; the trace is all-flux). Field + zero-copy views; NO
  per-face objects.
- Boundary = the domain-edge slots: `WavefrontFlux` x-edges (i=0, i=nx) ≙
  `BoundaryFlux` xmin/xmax; y-edges (j=0, j=ny) ≙ ymin/ymax. The orderings agree
  (current code proves it) → `ι*`/`ι_*` are clean per-face copies.
- Typed `ι_*(BoundaryFlux) → WavefrontFlux` (seed) + `ι*(WavefrontFlux) →
  BoundaryFlux` (absorb), single source of truth for the trace map (replaces the
  4 raw edge-slice assignments at each of the sweep + matvec sites).

### 3c. The biproduct (consistency with #208)
`WavefrontFlux ⊕ BoundaryFlux = C¹` mirrors the cell+trace `TimedFullField`
(`bulk: BulkField ⊕ boundary: BoundaryField`) one locus down (face locus). Keep
the `ι*`/`ι_*` as the projection/injection with `π_∂ ι_* = id` — the same
biproduct laws the operator algebra already uses.

---

## 4. Phasing (de-risk first; turn-by-turn; bit-identity every commit)

**PROACTIVE per the handoff protocol:** this is an operator-algebra carve
crossing the typed↔packed boundary → **dispatch test-architect FIRST** for the
bit-identity verification plan (which existing snapshots pin the raw behaviour;
the typed-`ι*` round-trip pin; the **L16 perf gate** — full `pytest tests/sn -q`
wall-clock must not regress, since this touches the hot path).

- **Phase 0 — DE-RISK (diagnostic, excluded).** Prototype the interior
  `FaceLayout` + `WavefrontFlux` wrapper on a small 2-D mesh. Prove: (a) the
  typed `ι_*`/`ι*` round-trip is BIT-IDENTICAL to the current raw seed/absorb
  (orderings agree); (b) the vectorized walk on `WavefrontFlux` views is
  bit-identical to the raw `psi_x`/`psi_y` walk (same `(N_oct,ng,n_diag)`
  indexing, no copy in the hot path); (c) a micro-bench shows no per-cell/
  per-face Python crept in (L16). STOP if values move or the hot path slows.
- **Phase 1 — mint `WavefrontFlux`** (type + interior `FaceLayout` + `ι*`/`ι_*`
  operators) with unit tests (units; field+views; the trace round-trip;
  zeros/from-buffer factories). No production wiring yet.
- **Phase 2 — wire the 2-D sweep ✅ DONE (`992b0c0`).** `_sweep_2d_wavefront`'s
  raw `psi_x`/`psi_y` → typed `WavefrontFlux` (face(0)/face(1) zero-copy views);
  seed/absorb → typed `wavefront.seed`/`wavefront.absorb` (`ι_*`/`ι*`). Gate:
  bit-identical (octant snapshots + Gate-K k_inf=1.875, 36 passed); L16
  perf-neutral (stash A/B median ratios ~1.00); elegance self-review.
- **Phase 3 — wire the 2-D matvec ✅ DONE (`0e3e16c`).** `_apply_2d_cartesian`'s
  raw `psi_x`/`psi_y` → typed `WavefrontFlux`; `streamed` dict → `wavefront.edge_view`
  (new accessor, removes the face→slot duplication). Gate: bit-identical (126
  passed: Gate-K, matvec≡sweep, BC extraction); **A2D-1 source-hash REFRESHED**
  (`f683f229…`→`12697ab3…`, deliberate-edit tripwire, behavior-neutral).
- **Phase 4 — ⚠ EVALUATED → DEFERRED to `nd_foundation` (2026-06-04).** The 1-D
  sweep/matvec is a parallel-prefix SCAN, NOT a wavefront: its interior fluxes
  are transient `(nx,K,ng)` chain-ordered scan output (`_compute_LpC` has no
  buffer at all), the structural antithesis of the cochain `(N,ng,nx+1)`, and
  the boundary is already typed. Forcing `WavefrontFlux` in = WRONG-FIT (risks
  the L16 cumprod, multiplies concepts). The ONE shared seam (the boundary-trace
  exchange + DD-closure averaging) unifies cleanly only when `nd_foundation`
  re-expresses both folds as one `d`-generic walk — recorded as the
  load-bearing collapse seam in `nd_foundation.md` §2.3 (incl. the HARD
  constraint that the collapse must not regress the d=1 scan's parallel-prefix
  efficiency). Evaluation: explorer memo `wavefront_flux_1d_tightening_verdict.md`.
  The §3a′ 3D-readiness (the type is axis-parametric) is the realized 1-D
  benefit; a future 3-D *Cartesian* sweep is a wavefront and uses WavefrontFlux.
  (Storage path (A) is delivered by Phases 0–3 for the 2-D wavefront — the
  locus where the interior cochain is persistent. 1-D needs no storage-A.)
- **Phase 5 — docs + retirement (A close-out).** archivist: the cochain frame +
  `WavefrontFlux ⊕ BoundaryFlux = C¹` biproduct in `operator_algebra.rst` /
  `discrete_ordinates.rst`; extend the storage×role grid (#205) with the FACE
  locus; retire the bare-array references (the `# ephemeral interior cache`
  comment becomes the typed story). The `BoundaryFaceFlux`-pending note in
  `transport/fields/__init__.py` resolved (no per-face type; views suffice).
- **Phase 6 — storage path (B): the moving-frontier window (STRAIGHT after A).**
  Swap the `WavefrontFlux` backing from the full per-axis face-field to a
  rolling 2-diagonal window (`O(N·ng·(nx+ny))`), behind the A-proven type/API.
  De-risk first (the converged solution is bit-identical; only the live
  working-set shrinks — pin peak memory drops + values unchanged). The walk
  becomes "advance the frontier" (consume incoming diagonal, produce outgoing).
  Its own elegance review + the L16 wall-clock gate (the rolling buffer must not
  add per-cell Python). 1-D may be a no-op (the frontier is a point); the win is
  2-D (and the future 3-D, where the volume↔surface ratio makes it largest).
- **(Later, separate plans):** the **G-S schedule** on the typed substrate via
  the `(octant×face)` reflective graph (`si_gauss_seidel_recovery.md`); then the
  **N-D / 3-D general foundation** (`nd_foundation.md`) — the future "tighten +
  generalize" session this plan's §3a′ makes 3D-ready by construction.

---

## 5. Verification (must-stay-green; the change is TYPE-ONLY → bit-identical)
Like B.5.2: `WavefrontFlux` wraps the SAME buffers; `.values` unchanged; only
the type + the named `ι*`/`ι_*` differ. So:
- **Bit-identity:** 2-D octant snapshots (`test_2d_octant_sweep_equivalence`),
  Gate-K k_inf=1.875 (`test_2d_l2_matvec_correctness`), A2D-1 hash (refresh),
  1-D `test_native_matvec` + `test_streaming_operator_decomposition`
  (Resolution-A twin), `test_bc_extraction_2d`/`_matvec`.
- **Correctness:** k_inf 1.875; SI≡Krylov; flat-balance; MMS L1 (1-D+2-D+
  curvilinear). **Perf (L16):** full `pytest tests/sn -q` wall-clock not
  regressed — THE risk for a hot-path retype.
- **Sentinel 36/36 (no `-O`).** Default `python -O`; env `PYTHONPATH=$PWD
  /Users/rodrigo/git/nuclear/ORPHEUS/.venv/bin/python`. **NO `continuous_get`**
  (#212). Bash: `> file 2>&1; echo exit=$?`.
- **PRE-EXISTING REDS (orthogonal):** cylinder matvec #206/#196; #212 hang;
  `test_krylov_curvilinear_precond_safety`/`test_b1pp`/`test_krylov_restart`
  (#200-adjacent Krylov budget).

## 6. Risk ranking
1. **Hot-path perf regression (HIGH, L16).** A typed wrapper that breaks the
   zero-copy vectorized view → per-cell Python → 10–20× slab regression. The
   walk MUST operate on `WavefrontFlux`'s backing arrays as zero-copy views with
   byte-identical indexing. Phase 0 micro-bench + the full-suite wall-clock gate
   are the guards.
2. **`ι*`/`ι_*` ordering (MEDIUM).** The typed trace must match the raw
   edge-slice ordering. Already proven to agree (the working `=`), but pin it.
3. **1-D / 2-D layout twin (MEDIUM).** `WavefrontFlux` must cover the 1-D
   direction-group edges AND the 2-D per-axis faces without a twin path
   (Cardinal Rule 2). Design the interior `FaceLayout` geometry-agnostically.
4. **Biproduct drift vs #208 (LOW).** `WavefrontFlux ⊕ BoundaryFlux` must
   compose like the cell+trace `TimedFullField`; reuse the projection/injection.
5. **Scope creep into storage-(B) (LOW).** Resist the moving-frontier rewrite in
   this carve — it is a separate, post-typing optimization.

## 7. The G-S forward pointer (why this is the right foundation)
The motivating SI Gauss–Seidel recovery lands on THIS substrate: the
`(octant×face)` reflective graph (attacker) — edges `(producer octant, face) →
(consumer octant, face)`; order-respecting edges fold into the inverted
triangular factor (G-S), cycle edges lag (Jacobi); KBA/Pautz is the literature
home; our `i+j=k` levels ARE the KBA diagonals. With typed `WavefrontFlux` +
`BoundaryFlux` + `ι*`, the G-S schedule is an explicit typed composition
(`sweep octant → ι* absorb → −B reflect → ι_* seed next octant`) rather than an
implicit buffer-timing dance. **Fold the `(octant×face)` construction into
`si_gauss_seidel_recovery.md` regardless of build order** (it sharpens that plan).

## 8. Pickup checklist (fresh post-compaction session)
1. Confirm branch `refactor/field-role-typing`, `git log --oneline -3` near
   `266fcf5`+, tree clean of tracked non-excluded paths.
2. Read the attacker memo (the cochain frame + field+views + the `(octant×face)`
   G-S graph) and this plan §0–§3 (the honest "representation-only" scope; the
   storage (A)-vs-(B) fork).
3. **Dispatch test-architect FIRST** (bit-identity + the L16 perf gate plan).
4. **Phase 0 de-risk** (typed `ι*` round-trip bit-identical + walk bit-identical
   + no hot-path slowdown) — STOP if it fails.
5. Phases 1→4 turn-by-turn, bit-identity + wall-clock gate each commit;
   elegance-enforcer per phase; A2D-1 hash refresh at Phase 3.
6. Phase 5 docs (archivist: cochain frame + `C¹ = C¹_int ⊕ C¹_∂` biproduct +
   storage×role FACE locus). Update `wave_o_operator_typing.md` + auto-memory.
7. Then the G-S schedule (`si_gauss_seidel_recovery.md`) on the typed substrate.
