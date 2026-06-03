# Wave O — Operator-Role Typing (BulkOperator / FullOperator / BoundaryOperator) + Typed Sources / Residuals

**Branch:** `refactor/moment-space-and-layering`
**Worktree:** `.claude/worktrees/moment-space-and-layering/`
**GitHub Issue:** [#208](https://github.com/deOliveira-R/ORPHEUS/issues/208)
**Phase status:** PENDING. Lands AFTER Wave T (`wave_t_tensor_network.md`). Precedes parent-plan steps P3.4 + P3.6.

**Date filed:** 2026-05-28 (during D-H.1b). Plan stub: 2026-05-30 (this file).

**Status:** ACTIVE on `refactor/field-role-typing` — see ⭐ CURRENT EXECUTION below (the §0–§6 stub is historical/stale).

---

## ⭐ CURRENT EXECUTION (2026-06-02) — branch `refactor/field-role-typing`, scope B (BC-extraction INCLUDED)

**This section supersedes §0–§6 below.** The 2026-05-30 stub assumed branch
`refactor/moment-space-and-layering` with Wave T *pending*; both are now stale
(Wave T is DONE here; B.3/B.5.2/B.6 of the field-role-typing refactor already
shipped most of #208's field-side substrate). Authoritative ground-truth map
(every `file:line` anchor): `.claude/agent-memory/explorer/issue_208_operator_algebra_surface.md`.

### Already shipped on this branch (feeds #208 — do NOT re-do)
- **L1 operator algebra** (`numerics/operator.py`): `LinearOperator` is a
  `@runtime_checkable Protocol`; `LinearOperatorMixin` installs dunders +
  `.H`/`adjoint()`. **The adjoint closure laws `(A+B).H`, `(A∘B).H`, `(αA).H`
  are FULLY WIRED** (`OperatorSum`:619 / `OperatorProduct`:691 /
  `ScaledOperator`:748), cap-gated by `CAP_APPLY_TRANSPOSE`. Capability tags +
  `MissingCapability` raised at composition time.
- **B.5.2**: bulk operator `.apply` outputs retyped `AngularFlux →
  AngularSourceSink`; `ZeroOperator(codomain_zero=…)` (operator.py:822;
  docstring already cites #208).
- **B.3/B.5.1**: bulk residual leaves `AngularResidual`/`ScalarResidual` (+
  `from_balance`); `BoundaryResidual` (residuals/boundary_residual.py:82) has
  `from_balance` AND its docstring documents the held mistyped wiring.
- **B.6**: locus-typed `TimedFullField` (`bulk: BulkField`, `boundary:
  BoundaryField`) — this IS the "TimedFullAngularSource/Residual composite" the
  stub imagined; **no new wrapper type needed**.
- `InflowSourceSpec` rename (`c7060f0`) freed the `BoundarySource` name; the L2
  boundary leaf is `BoundarySourceSink` (source_sinks/boundary_source_sink.py:113)
  — minted, EMPTY, no consumer.
- **Boundary architecture is 3-layer** (Issue #186): `BoundaryTraceLaw` ABC
  (geometry/boundary/_base.py:74, pure affine `γ₋ψ=R·G·γ₊ψ+q`, NO `apply`) →
  `SNBoundaryRealizer.realize(law)→LinearOperator` (boundary_realizer.py:123,
  the sole bridge) → realized L1 op.

### Scope decision (2026-06-02): INCLUDE BC-extraction (Option B)
The full architecturally-complete Wave O: `(L_full + C − S − F − B)ψ = q`
canonical, `B` a top-level `BoundaryOperator`. Rationale (user):
architecture-forward (`[[feedback_architecture_forward_not_legacy_fit]]`) — the
cleanest end state (bare bulk-streaming Full leaf + algebraic boundary adjoints)
and it retires the BC-absorbed sweep rather than building an adjoint around it.
Accepted cost: **O.4 restructures the geometry-non-uniform, L1-verified sweep DAG.**

### Operator-algebra framing (LOCKED 2026-06-02 — dagger biproduct category + block metric)
The native structure (cross-domain-attacker, twice-confirmed; memory
`.claude/agent-memory/cross-domain-attacker/issue_208_operator_algebra_frames_full_access.md`):
operators form a **dagger inverse biproduct category with a non-trivial block
metric**. What this MANDATES for Wave O:

- **Biproduct** — `V = V_bulk ⊕ V_boundary` ⟹ operators ARE block matrices
  (a theorem, not a convention); N-way `⊕` (kinetics `flux⊕precursors`) is free.
  *Sanity check during O.1/O.2:* `OperatorSum.apply` must be
  biproduct-by-construction (a 3rd `⊕` slot needs NO new `if`), else the N-way
  story is aspirational.
- **Metric `G` is LOAD-BEARING (correctness, not optional).** The dagger `†` is
  the **G-adjoint** `A† = G⁻¹AᵀG`, NOT the plain transpose; boundary block
  `G_s = |Ω·n|·w_n`. "Adjoint via order-reversal" yields only the transpose → the
  WRONG adjoint of the off-diagonal coupling until `G` is populated. → Design
  decision #2.
- **Dagger functor** — the already-wired `.H` closure laws `(A+B).H` / `(A∘B).H`
  / `(αA).H` ARE its axioms; algebra-level adjoint composition is free once each
  leaf's `.H` is the correct G-adjoint.
- **`B` is a gluing cocycle** (reflective = deck transformation `Ω→R(Ω)`;
  periodic = base quotient), not a generic `A_ss` — `B†` is a free
  inverse-permutation, and O.4 extraction = lifting the cocycle into
  global-section assembly (shares the MoC fiber-bundle frame). Independent
  validation that Option B (extraction) is right.
- **Transport-resolvent backbone** — the four method `solve`s are quadratures of
  one resolvent `(Ω·∇+Σ_t)⁻¹`; diffusion is its asymptotic LIMIT → predicts the
  self-adjoint-elliptic exception. Confirms the cross-method mixin layering
  (universal `apply+solve+adjoint` base + `TriangularSolve`(SN/MoC) /
  `DenseFactorization`(CP) / `SelfAdjointElliptic`(diffusion); causal-ordering
  stays SN-local per `[[feedback_unify_after_two_instances]]`).

**DEFERRED to stay lean (filed #213):** capabilities as a morphism class
`Iso ⊂ PartialIso ⊂ General` (restriction/inverse category — the principled
replacement for the string-tag capability set). Orthogonal axis (invertibility,
not block-touch); Wave O keeps the cap-tags and hand-wires
`StreamingOperator.solve_transpose` (the reverse sweep).

### Reframed substep ledger (status vs current branch)
| Step | Scope | Status (explorer) |
|---|---|---|
| **O.1 ✅** | `BlockRole` enum + `BulkOperator`/`FullOperator` markers (value-based `isinstance` metaclass) @ L1; tag C/S/F=BULK, L/InvertibleOperator=FULL. `797b505` | DONE — bit-identical (numerics 572 / operators 397 / primitives 271 / transport 200 / sweep-core 326, Gate 1.3 still xfail; sentinel 36); +16 foundation tests; elegance-reviewed APPROVE-WITH-NITS. BoundaryOperator marker deferred to O.4. |
| **O.2** | adjoint: `StreamingOperator.apply_transpose` + `OperatorSum` Protocol routing → **Gate 1.3 green** | ~50% (L1 closure wired; `CollisionOperator` adjoint ships; MISSING the streaming adjoint + routing) |
| **O.3a** | boundary-residual typing (apply-output boundary `BoundaryFlux→BoundaryResidual` + C/S/F zeros) | **MERGED into O.4a** — entangled with extraction (the pre-extraction boundary slot is overloaded; see the O.3a-entanglement finding above). The typing lands as a consequence of O.4a, not before it. |
| **O.3b** | typed full source (`BoundarySource`) + `BoundarySourceSink.from_spec` + non-vacuum-BC MMS | UNTOUCHED — sequenced AFTER O.2 (needs `B` / `q_inflow` from O.4). Leaves minted (B.3); genuinely-new compute. |
| **O.4** | **BC-extraction** → sibling `B` BoundaryOperator; bare bulk-streaming `L_full` | UNTOUCHED (highest risk; geometry-non-uniform) |
| **O.5** | retire implicit-zero-boundary shortcuts | UNTOUCHED (gated on O.1+O.3) |

### O.1 close-out — carry-forward acceptance criteria (from the elegance review)
The O.1 elegance review (APPROVE-WITH-NITS) surfaced two forward-composition
gaps that are correct lean boundaries at O.1 (no live `block_role` consumer yet)
but MUST be pinned so they don't calcify:

- **O.2 — adjoint propagates the role.** Block-role is adjoint-invariant
  (`A_bb† : bulk→bulk`, `A_ss† : boundary→boundary`; a full 2×2 transposes to a
  full 2×2). `_AdjointOperator` (operator.py:459) currently reports
  `block_role=None`, so `isinstance(L.H, FullOperator)` is FALSE today. O.2 MUST
  give `_AdjointOperator` a `block_role` property returning the inner operator's
  role. **Acceptance test:** `(L).H` is a `FullOperator`. (Else the adjoint
  composite falls into the unclassified branch and O.2's role-dispatch routes it
  wrong — a path bug that passes shape checks.)
- **O.2 — derive the sum role + RETIRE the hardcoded line (twin-path).** O.1 sets
  `InvertibleOperator.__init__`'s role with a hardcoded `self.block_role =
  BlockRole.FULL`; generic `OperatorSum` stays `None`. When O.2 adds general
  sum-role derivation from operands (FULL if any operand FULL, else BULK if any
  BULK, …), it MUST **retire** the hardcoded `InvertibleOperator` line in the
  SAME change — otherwise the special-case and the general derivation are a twin
  path (anti-pattern #1 / ERR-026 shape: a later operand change updates one and
  the stale copy diverges).
- **O.4 — retire the deprecated boundary alias when the marker lands.**
  Introducing the `BoundaryOperator` block-role marker MUST retire the
  deprecated `geometry.boundary.BoundaryOperator` alias (a misnamed
  `BoundaryTraceLaw`) in the same wave — no dangling two-name overlap
  (`[[feedback_aggressive_retirement]]`).

### O.4 design sketch (the load-bearing carve — refine when it starts)
- `B` = `BoundaryOperator`: reads `ψ.boundary.outflow`, writes
  `ψ.boundary.inflow` via the realized `R·G` law; the affine inflow `q_inflow`
  moves to the RHS `q.boundary`.
- `L_full` = bare bulk-streaming `FullOperator`: given bulk ψ +
  `ψ.boundary.inflow`, produce bulk balance + `ψ.boundary.outflow`. No BC logic
  inside.
- Boundary block of the residual `(…−B)ψ−q`: `ψ.inflow − (R·G·ψ.outflow +
  q_inflow)` — GMRES/SI drives it to zero. **This IS the held boundary residual
  (O.3), generalized** — O.3 and O.4 are the same wiring at two levels.
- Today's reflective fixed-point lives INSIDE the sweep (reflective wrap
  operator.py ~520); extraction moves it OUT to the outer solve (which already
  drives the full residual to zero).
- **Geometry-non-uniform → two sub-carves**: **O.4a 1-D** (slab/sphere/cyl;
  extract from `_compute_LpC`, where `bc_outer/bc_inner.apply` + reflective wrap
  live) then **O.4b 2-D Cartesian** (extract from `_apply_2d_cartesian`; this
  path forces the boundary residual to ZERO today (line 1450) → 2-D must ADD the
  boundary-residual compute, a genuinely new term).
- **Equivalence class for O.4 is CONVERGENCE-equivalence, NOT bit-identity**
  (`vv-principles` §bit-identity): moving the reflective coupling from
  inside-the-sweep to the outer solve changes the iterate trajectory → iterates
  differ, but the converged ψ solves the same linear system, so it matches to
  solver tol. Verify against MMS/analytical converged ψ + preserved L1 MMS
  slopes. **Vacuum-BC cases stay bit-identical** (B=0; given-zero-inflow sweep =
  today's sweep).

### Design decisions (baked in)
1. **Q1 — keep the 3-layer boundary split.** `BoundaryTraceLaw`@L1 → `realize` →
   L1 op stays; O.1 adds a `BoundaryOperator(LinearOperator)` Protocol that
   *classifies the realized op* (the split is already type-enforced; don't
   collapse).
2. **`|Ω·n|` boundary metric in O.2 — LOAD-BEARING (the dagger `G`).** The
   adjoint is the G-adjoint `A† = G⁻¹AᵀG`; the boundary block carries
   `G_s = |Ω·n|·w_n` (partial current `J± = ∫|Ω·n|ψ`). The hook EXISTS:
   `TraceSpace` carries `omega_dot_n` + `FaceLayout` leaf-data, and
   `_AdjointOperator.apply` (operator.py:497) already reads
   `inner_product_weights`. O.2 MUST: (a) populate
   `TraceSpace.inner_product_weights = |omega_dot_n| ⊙ w_n`; (b) make each
   leaf's `.H` the G-adjoint; (c) add a **discriminating Gate-1.3 test —
   reciprocity in the G-metric on a NON-trivial off-diagonal block** (today's
   Euclidean bulk-only reciprocity cannot detect a metric-blind boundary
   adjoint; without (c) the bug is invisible).

### Sequencing (LOCKED 2026-06-02; AMENDED — O.3a MERGED into O.4)
**O.1 ✅ → O.4a (1-D extraction — ALSO lands the boundary typing) → O.4b (2-D
extraction) → O.2 (adjoint on the bare shape + Gate 1.3 green) → O.3b (typed
source / `from_spec` / non-vacuum MMS) → O.5.**

**⚠ O.3a-entanglement finding (2026-06-02 — why O.3a is no longer a standalone
precursor).** Attempting O.3a (retype the apply-output boundary `BoundaryFlux →
BoundaryResidual`) FAILED with a `TimedFullField boundary type mismatch:
BoundaryResidual vs BoundaryFlux` in the source-iteration source sum. Root
cause: the pre-extraction boundary slot is **overloaded**. The Krylov matvec
wants `L.apply(ψ).boundary` = the BC-absorbed defect `γ₊ψ − bc_estimate`
(`BoundaryResidual`), and the `OperatorSum` closure forces C/S/F's boundary
zeros to match (all-or-nothing — retyping only `L` breaks the sum). But the
source-iteration RHS `S.apply(ψ) + F.apply(ψ) + q_ext` is *typed* arithmetic,
and `q_ext.boundary = self._boundary_flux` (solver.py:572) is a `BoundaryFlux`
**carrying the reflective BC inflow trace** — the sweep reads `rhs.boundary` as
its BC seed (a flux). The SAME `S.apply` boundary output is summed in both
contexts, so it cannot be `BoundaryResidual` (matvec) AND the flux-seed type
(SI) at once; forcing uniform `BoundaryResidual` would mislabel the inflow-flux
seed as a residual (a Cardinal-Rule-1 type lie). **There is no honest uniform
type for the pre-extraction boundary slot** — this IS the source/residual
boundary asymmetry #208 exists to fix, and it is load-bearing exactly as the
dagger-framing predicted. O.4's extraction resolves it by construction (inflow →
`B`'s domain = a flux; outflow consistency defect → the residual), so the clean
boundary typing is a *consequence* of extraction, not a precursor. The O.3a
production edits were applied then **reverted** (HEAD = `797b505` O.1; the
revert is green). The boundary-typing retype now lands AS PART of O.4a.

Both the main agent and the test-architect independently recommend
**O.4-before-O.2**: the O.2 adjoint built around the BC-absorbed sweep would be a
THROWAWAY (O.4 deletes that forward shape). O.4 is verified by forward-only
references (MMS / `Q/Σ_t` / homogeneous `k_∞`) so there is NO dependency
inversion; the restructuring anchors (Gate 1.1, MMS slopes, Resolution-A
bit-exact, sentinel) are already green and need no adjoint. Build the adjoint
ONCE on the final bare-bulk-streaming shape. (Rejected alternative
`O.1→O.2→…→O.4→O.5` would flip Gate 1.3 green earlier but discard the BC-absorbed
adjoint.)

### Verification strategy (test-architect plan: `.claude/agent-memory/test-architect/issue_208_wave_o_verification_plan.md`)
Per-substep equivalence class + gate:
- **O.1** Protocols — BIT-IDENTITY (pure tagging).
- **O.2** L adjoint + routing — BIT-IDENTITY forward; the reverse sweep is NEW,
  verified by Gate 1.3 reciprocity vs a dense-probe transpose.
- **O.3a** boundary-residual retype — VALUE-IDENTITY (matvec already computes the
  face defect; only the python type changes).
- **O.3b** from_spec + non-vacuum MMS — GENUINELY-NEW compute (all 3 bit-identity
  criteria pass: named inflow-trace intermediate / MMS-imposed inflow reference /
  exact single step).
- **O.4** BC-extraction — CONVERGENCE-EQUIVALENCE (non-vacuum: iterates differ,
  converged ψ matches an EXTERNAL reference under the `iter×cond×ULP` drift
  bound) / BIT-IDENTITY (vacuum, B=0).
- **O.5** retire shims — BIT-IDENTITY.

New tests: Protocol-conformance with **exclusivity rows as the negative half**
(L11); `OperatorSum` dispatch (bulk leaves → zero boundary, full leaf's residual
survives); **Gate 1.3 reciprocity** — EXTEND to +slab (cheap vacuum anchor) and
+2-D-as-xfail (the 2-D adjoint is a separate harder deliverable, lands with
O.4b); the **`|Ω·n|·w_n` boundary-block adjoint** reciprocity (NEW — Gate 1.3
tests only `.bulk.values` Euclidean today, so the partial-current metric is
unpinned; pair with a negative control proving the weighting is load-bearing);
O.4 convergence-equivalence (reflective/albedo/white/periodic across geometries;
vacuum bit-identity); non-vacuum MMS with the **anti-bias ansatz `(A(x)+μB(x))/W`**
(overrides the inherited isotropic-vanishing `sin(πx/L)` bias).

Must-stay-green every commit: operator-algebra-core gate; L0 streaming-equilibrium
(Gate 1.1); L1 curvilinear + 2-D MMS; Resolution-A bit-exact; sentinel. **No
`continuous_get` eigenvalue gate** (blocked by the registry bug, issue #212 —
use direct MMS / homogeneous-`k_∞`).

**Gate 1.3** flips green at O.2 for sphere+cyl (+slab); 2-D Cartesian reciprocity
is the separate harder adjoint, lands with O.4b.

### Open architectural questions for O.4 (resolve before O.4 test code; test-architect §7)
1. Who populates `TraceSpace.inner_product_weights = |Ω·n|·w_n` (and when — O.2
   needs it for the boundary-block adjoint).
2. 2-D adjoint timing — deferred to O.4b (the 2-D reverse sweep + the
   newly-added 2-D boundary residual).
3. Is `B*` free via the already-wired `(A+B).H` closure once `B` is a
   Protocol-classified `BoundaryOperator` with its own analytic `.H`?
4. Non-vacuum MMS ansatz amplitude (scattering ratio `c<1`, `A(x)>0` for
   physicality).

### QA watch (test-architect flag)
Candidate **Mode-8 "passive-boundary blind spot"**: a zero-by-construction
boundary slot satisfies any bulk-only test vacuously (generalizes H3 to the
trace equation). NOT yet a `vv-principles` row / ERR-NNN — earns one only if
O.4b surfaces a concrete bug. QA to watch during O.4b.

### Cross-refs
Ground-truth map: `.claude/agent-memory/explorer/issue_208_operator_algebra_surface.md`.
Parent field-role plan: `field_role_typing_view_g.md` (#208 sequenced after
B.6). Memory: `[[project_wave_o_operator_algebra]]`.

---

## ⭐⭐ O.4 EXECUTION PLAN (detailed — THE NEXT OPERATION; absorbs O.3a)

**Compaction-resume pointer.** Branch `refactor/field-role-typing`; HEAD = `d7e1316`
(O.4a.2 Commit 1). Before executing, READ: the `⭐⭐⭐ O.4a.2 IN-FLIGHT STATE`
block immediately below (the live working-tree state), then this section +
extraction-surface map
`.claude/agent-memory/explorer/issue_208_o4_extraction_surface.md`, the
**flip-wiring map** `.claude/agent-memory/explorer/issue_208_o4a2_flip_wiring_surface.md`
(the `file:line` for the solver/sweep/−B surgery), the test-architect matvec spec
`.claude/agent-memory/test-architect/issue_208_o4a2_matvec_extraction_spec.md`,
and the verification plan
`.claude/agent-memory/test-architect/issue_208_wave_o_verification_plan.md`.

---

## ⭐⭐⭐ O.4a.2 — ✅ COMPLETE (HEAD `2bdc66d`, 2026-06-03)

**`(L_full + C − S − F − B)ψ = q` is now canonical on BOTH the matvec/Krylov path
(Commit 2 `4c0ff96`) AND the SI/bare-sweep path (Commit 3 `2bdc66d`).** The
bc-in-sweep mechanism is RETIRED for 1-D (2-D wavefront stays bc-in-sweep until
O.4b). Sphinx theory update dispatched to the archivist.

**Commit 3 (`2bdc66d`) — bare sweep + SI −B:** `transport_sweep` 1-D entry reads
the seeded inflow directly (no `bc.apply`); `_solve_timed_full_field` seeds from
`rhs.boundary` (NOT `initial_guess.boundary` — the "leave as-is" guidance below
was WRONG); `_solve_source_iteration` folds `S+B` + retires
`q_ext.boundary=_boundary_flux`; new `_reflect_outflow_into_inflow` helper drives
−B for the DIRECT loops (`_solve_fixed_source_si` + final eigenvalue sweep),
1-D-guarded (`reduced is not None`). **TEST-MIGRATION FINDING (handoff
under-scoped):** the bare sweep broke 5 tests calling `transport_sweep` directly
in reflective loops (relying on bc-in-sweep) — all migrated to drive −B via the
helper + 2 invertible-operator tests to the `rhs.boundary` seed contract.
**#212 (PRE-EXISTING):** the full eigenvalue suite hangs on the `continuous_get`
Peierls registry eager-build (3 tests deselected; solver converges ~0.3s).
**DEFERRED to O.2** (elegance-enforcer): unify the 3 parallel `(L+C)` matvec
impls; honest `L+C−S−F−B` driver composition (retires the `S+B` fold — type-checks
only because `S.domain is None`, the forcing tripwire) + the two −B routes
collapse; a trace-only `B.reflect_into_inflow(boundary)` entry point (retires the
helper). Gates: vacuum bit-identical, reflective convergence-equiv
(streaming-equilibrium both solvers, k∞, si_carve→k∞, keff, invertible Q/Σ,
sentinel no-`O`), 506 + 73 + keff_2d 4.

---

### ✅ COMMIT 2 LANDED (HEAD `4c0ff96`) — the co-land was SPLIT into 2 commits

**The "co-land everything" framing below was SPLIT** (2026-06-03) after grounding
revealed the bare-sweep is a delicate GLOBAL `transport_sweep` change touching 3
callers + the seed-source logic. Splitting keeps every gate green at each commit
(honours the co-land INTENT — no broken intermediate state) while de-risking:

- **Commit 2 (DONE, `4c0ff96`):** the **−B BC-extraction on the matvec/Krylov
  path**. Matvec flip (`_compute_LpC`+twin: keystone deleted, given-inflow read,
  **KEPT** the outflow defect `streamed−ψ.outflow`, ADDED the inflow identity;
  dead `bc_outer`/`bc_inner` retired → explicit `curvature!="cartesian"`
  case-split), `B` **inflow-row projection** (`SNBoundaryOperator._apply_faces`),
  `SNSolver._scattering_with_boundary_op` = `S+B` folded into BOTH Krylov S-args.
  **Krylov path on `−B`; SI path UNCHANGED** (still bc-in-sweep + partner-flux
  seeding — it routes through `(L+C).solve`, a separate WDD sweep, not the flipped
  matvec). Also `0442dce` refreshed a pre-existing stale A2D-1 hash pin
  (orthogonal, from O.4a.1-α). Gates green: vacuum bit-identical (matvec 18
  baselines + regression snapshots), reflective convergence-equiv (curvilinear
  streaming-equilibrium krylov no-`O`, k∞, keff), operator suite 455,
  elegance-enforcer CONCERNS resolved.
- **Commit 3 (DONE, `2bdc66d`):** the **SI bare-sweep extraction** (see the
  COMPLETE banner above). Retired the bc-in-sweep mechanism for 1-D. The
  "REMAINING WIRING" / mechanics below are the historical handoff — accurate
  except the corrected items called out in the banner (seed-from-`rhs.boundary`;
  the 5-test migration; the helper for the direct loops).

**TWO design findings the original handoff MISSED (corrections):**
1. The whole-trace `B` needed **inflow-row PROJECTION**. The full-face realized
   law (specular permutation) leaks `R·ψ.inflow` onto the OUTFLOW slots, which
   `−B` would subtract into the outflow-definition residual. Fixed: `apply` emits
   on inflow slots only, `apply_transpose` on outflow slots only. (Empirically
   confirmed before the fix.)
2. `_solve_timed_full_field` seeds the sweep from `initial_guess.boundary`,
   IGNORING `rhs.boundary` (where the `Bψ` source lands). So the Commit-3 bare
   sweep MUST seed inflow from `rhs.boundary`, NOT the iterate — **the original
   "leave `_solve_timed_full_field` as-is" guidance is WRONG.**

**DEFERRED to O.2** (Cardinal-Rule-2 flags surfaced by the elegance review): the
3 parallel `(L+C)` matvec impls (`_compute_LpC` / `_compute_decomposition` / the
solve-sweep) want unification; the honest `L+C−S−F−B` driver composition (so the
`S+B` fold retires — it type-checks ONLY because `S.domain is None`, the O.2
forcing-function tripwire); realized-law self-projection (so `B`'s mask is
intrinsic).

---

### Original "co-land" handoff (below) — superseded by the split above; mechanics still accurate

**HEAD `d7e1316`. The matvec side of the co-land flip is DONE but UNCOMMITTED in
the working tree** (it persists across compaction — it is on disk, not in
context). The flip commits ATOMICALLY once green; do NOT revert these, do NOT
`git checkout`/`stash` them.

**Uncommitted working-tree manifest (verify with `git status` post-compaction):**
- `orpheus/sn/operator.py` — **DONE: the design-(a) matvec flip** in BOTH
  `_MSpatialOperatorSum._compute_LpC` AND its `_compute_decomposition` twin:
  (1) keystone `inflow_full = bc_outer.apply(outflow_at_boundary.T)` DELETED;
  backward sweep seeds from the GIVEN `face_outer`; (2) slab inner seed +
  `outer_inflow_estimate` read the GIVEN trace (not `bc_*.apply`); (3) boundary
  emission KEEPS the outflow defect `ψ.outflow − streamed` (the `I·ψ.outflow`
  diagonal — UNCHANGED, vacuum bit-identical) and ADDS the inflow identity
  `ψ.inflow` on the inflow slots (the `I·ψ.inflow` diagonal). Curvilinear pole
  Carlson seed KEPT.
- `tests/conftest.py` — **DONE:** added the `--capture-baseline` pytest option.
- `tests/sn/_data/` (untracked) — **DONE:** 18 captured pre-carve vacuum baselines
  (slab/sphere/cyl × seed 0/1/2, bulk + boundary).
- `tests/sn/operators/test_bc_extraction_matvec.py` (untracked) — the
  test-architect's matvec gate file. NEEDS RECONCILIATION (see below).

**VALIDATED (apply level):** vacuum bit-identity 18/18 byte-identical (bulk +
boundary, slab/sphere/cyl); `Q/Σt` no pole spike; `(L+C).apply == M_spatial.apply
+ M_angular_redist.apply` consistent (slab + sphere); 3b (`TestLFullReadsInflow`)
flips XPASS for all 3 geometries (L_full now reads `ψ.boundary.inflow`).

**KEY DESIGN CORRECTION (supersedes the "raw outflow" prose in the mechanics
bullets below — those are WRONG, fix them):** the canonical `(L_full + C − S − F
− B)ψ = q` requires `L_full` to KEEP the outflow self-consistency defect (the
`I·ψ.outflow` diagonal — `ψ.outflow` is a stored unknown `−B` reads as its
input). "Emit raw outflow / drop the subtraction" would make `ψ.outflow`'s row
singular and break `−B`-as-sibling. The ONLY matvec changes are: delete keystone,
read given inflow, ADD the inflow identity. The outflow emission is UNCHANGED
(today's `computed − stored`, sign-free since `q.outflow ≡ 0`).

**TWO LOAD-BEARING CONSTRAINTS (from the flip-wiring map):**
1. `OperatorSum` does NOT propagate `CAP_SOLVE` (`numerics/operator.py:689-692`)
   → `−B` CANNOT join `LC = L+C` (strips `solve`). It FOLDS into the subtracted
   S-argument: `self.scattering_op + SNBoundaryOperator(self.sn_mesh)` (S FIRST so
   the domain-check skips — S.domain is None; `SNBoundaryOperator.domain =
   sn_mesh.trace`). One change serves BOTH paths: Krylov matvec →
   `(L+C).apply − (S+B).apply − F.apply = (L+C−S−F−B)`; SI source →
   `LC.solve((S+B+F)ψ + q)` = the `Bψ` term.
2. SI double-counting → the sweep MUST go BARE NOW (read seeded inflow directly,
   NO in-entry `bc.apply`); else `bc.apply(B·ψ.outflow)` double-applies B. (Krylov
   tolerates the BC-absorbed sweep as a mere preconditioner; SI's sweep IS the
   solution map.)

**REMAINING WIRING (the next session's work — `file:line` from the flip-wiring map):**
1. **`−B` S-arg fold** (`solver.py`): `_solve_krylov` (~733-746), `_solve_source_iteration`
   (~592-598), `_solve_fixed_source_krylov` (~1461-1473): replace the
   `self.scattering_op` S-argument with `self.scattering_op +
   SNBoundaryOperator(self.sn_mesh)`.
2. **Bare sweep** (`sweep.py`): drop entry `bc_*.apply` (slab `bc_left/right_obj.apply`
   ~473-474; curvilinear `bc_outer_obj.apply` ~488) — read the seeded inflow trace
   directly; persist the RAW outflow (~585/587 slab, ~730 curv). (2-D wavefront
   ~889-900 is O.4b — leave.)
3. **Retire the SI partner-flux seeding** `solver.py:572` (`q_ext_composite.boundary
   = self._boundary_flux`) — the inflow becomes a live solved unknown carried in
   `ψ.boundary` (threaded via `to_flat`). `_boundary_flux` as the FINAL-solution
   carrier (`1011/1042/1276/1316`) STAYS.
4. **`_solve_fixed_source_si` direct loop** (`solver.py:1236-1278`, slab/2-D default)
   — no S-arg triple; add the `Bψ` term to its source explicitly.
5. **`_solve_timed_full_field`** (`operator.py` ~1988-2025) — the inverse action;
   leave as-is (consistency enforced by the outer −B).

**GATE MIGRATION + TEST RECONCILIATION (retirement = test migration):**
- `tests/sn/operators/test_bc_extraction_matvec.py`: **unmark 3b**
  (`TestLFullReadsInflow` — flipped XPASS); **rewrite 3c** (`TestLFullEmitsRawOutflow`
  → assert the inflow-IDENTITY emission) + **GATE4** (`TestVacuumBoundaryDefectVsRaw`
  → assert the outflow-DEFECT is KEPT, i.e. output DOES depend on input outflow);
  **delete 3a** (`TestWholeTraceBOperator` + `_realized_B_for_face` — the whole-trace
  `B` is canonically tested in `test_sn_boundary_operator.py`). Fix the naming
  mismatch: the test imports `assemble_boundary_operator` from `orpheus.sn.operator`
  — migrate it to `from orpheus.sn.boundary_operator import SNBoundaryOperator`.
- `tests/sn/sweep/core/test_phase_c_gates.py`: Gate 1.1 per-ordinate flat-flux
  (~279-322) — migrate to `(L+C−B)` / a consistent reflective inflow trace; the
  **`==2 bc_right.apply` call-count assert** in
  `test_bc_trace_contract_capture_and_compare_sphere` (~673) → the flipped matvec
  calls `bc_outer.apply` **0 times** (reads given `face_outer`) → migrate to 0 (or
  redirect to `B.apply`).
- `tests/sn/operators/test_streaming_operator_decomposition.py` Resolution-A
  (~147-204, reflective, zero-input boundary): the reference `_LC_matvec` helper
  must match the flipped emission (incl. the inflow-identity rows);
  `TestResolutionADifferentFromPriorWrong` reaches into `_compute_decomposition`
  (~327).
- `tests/sn/sweep/curvilinear/test_streaming_equilibrium_curvilinear.py` (~137-140,
  full-solver reflective) — the through-solver acceptance gate for the `−B` wiring.

**FULL GATE (commit only when green):** vacuum bit-identical END-TO-END (through
solver) + reflective CONVERGENCE-EQUIVALENCE — converged ψ matches MMS / `Q/Σt` /
homogeneous `k∞` to solver tol (iterates differ, converged ψ matches an EXTERNAL
reference); the old-vs-new drift tripwire `< 2×solver_tol`. **NO `continuous_get`**
(fixture bug #212). Then: fix the "raw outflow" prose in the mechanics bullets,
run the operator-algebra-core + sentinel (no -O) gates, elegance-enforcer on the
full flip, commit, update plan + memory.

---

### The keystone (single load-bearing change)
The matvec line `inflow_full = bc_outer.apply(outflow_at_boundary.T)`
(`operator.py:521` in `_compute_LpC`, `:795` in `_compute_decomposition`) is the
intra-call reflective re-apply that makes `L` couple bulk↔boundary. **Deleting
it decouples `L`**; the `R·G` it performed moves to the sibling `−B`. Everything
else follows from this.

### The algebra (block structure on `V = V_bulk ⊕ V_inflow ⊕ V_outflow`)
```
            ψ_bulk         ψ_inflow      ψ_outflow
r_bulk    [ L_bb+C−S−F     L_{b,in}        0       ]        = q_bulk
r_inflow  [ 0              I              −B        ] · ψ   = q_inflow   (consistency)
r_outflow [ −L_{out,b}     0               I       ]        = 0          (definition)
```
- **`L_full` (FULL):** `L_bb` (bulk balance), `L_{b,in}` (inflow seeds the sweep),
  `−L_{out,b}` (bulk produces the outflow trace). Reads `ψ.inflow` as a GIVEN;
  writes `ψ.bulk` + `ψ.outflow`. NO BC logic.
- **`C, S, F` (BULK):** `A_bb` only.
- **`B` (BOUNDARY):** `ψ.outflow → ψ.inflow` via the realized `R·G`. The affine
  inflow `q_inflow` → `q.boundary` (`BoundarySource`).
- The identities `I·ψ_in`, `I·ψ_out` make the trace values unknowns; the outer
  Krylov/SI drives `r_inflow = ψ.inflow − B·ψ.outflow − q_inflow → 0` and
  `r_outflow = ψ.outflow − streamed(ψ.bulk) → 0`. **Placing the boundary
  identities is the core wiring task of O.4a.2** (map Area 1).

### SI ⇄ Krylov unification (the elegance)
Today the reflective coupling is implicit and DIFFERENT in the two paths (matvec
re-applies intra-call @521; solve-sweep persists outflow across outer iterations
via `self._boundary_flux`). Post-extraction BOTH use `B` explicitly: Krylov
solves `(L_full+C−S−F−B)ψ=q` via `.apply`; SI is the splitting
`(L_full+C).solve((S+F+B)ψ + q)` — the `Bψ` term supplies the inflow each
iterate. One operator `B`, two consumers (apply vs the SI source).

### Boundary typing (the absorbed O.3a — now honest, no overload)
- `ψ.boundary.inflow` / `ψ.boundary.outflow`: `BoundaryFlux` (fluxes; `B`'s
  domain/codomain).
- the consistency defect `r_inflow`: `BoundaryResidual` (the named held residual).
- `q.boundary` (prescribed inflow): `BoundarySource` (zero for vacuum/reflective;
  non-zero only for `PrescribedInflow`).
The pre-extraction overload is GONE: inflow (flux) and defect (residual) are now
different slots, not one slot wearing two hats.

### Sub-steps (turn-by-turn, each its own commit + gate)
- **O.4a.1 — `B` as a `BoundaryOperator` leaf. [✅ DONE 2026-06-03 —
  α `9a7f216` / β1 `a44fac5` / β2 `1069755` / γ `91d3249`].** Tagged the
  `SNBoundaryRealizer.realize` outputs with `block_role = BlockRole.BOUNDARY`
  (producer-site `_as_boundary` helper, 6 linear branches; `PrescribedInflow`
  NOT tagged — it is `q.boundary`), landed the `BoundaryOperator` marker
  (`numerics/operator.py`), `_BoundBoundaryOperator` forwards the role.
  Sub-commits: **α** prod alias→canonical rewire (boundary_realizer + 2-D
  defaults; bit-identical); **β1** ~139 test-ref migration (incl. deleting the
  alias-identity test); **β2** alias-def deletion + full doc-correctness sweep
  (archivist: 6 `:class:` xref fixes + prose + 4 rst pages; the
  BulkOperator/FullOperator/**BoundaryOperator** MARKER refs LEFT — distinct
  from the retired geometry alias); **γ** marker + tagging + tests (foundation:
  all 9 linear realizations advertise `BoundaryOperator` + exclusivity,
  PrescribedInflow negative, mesh `bc_*` forwarding). Verified bit/value-
  identical, Sphinx clean, elegance-enforcer PASS.
  **REFINEMENT vs this plan:** `domain = codomain = mesh.trace` + the whole-trace
  apply lift were **deferred to O.4a.2** — they are coupled to the `−B`
  OperatorSum consumer (the whole-trace `B` is a block-diagonal-over-faces
  assembly whose shape is *defined by* that wiring; build-primitive-not-product,
  `coding-elegance` Pattern 6). O.4a.1's gate (`B.apply(outflow)` ≡ legacy
  `bc.apply`) is per-face and satisfied by construction (apply untouched).
  Realized-op adjoint readiness (for O.2): reflective `PermutationOperator&I`
  (involution, `apply_transpose` ✓), vacuum `IncomingOrdinateMaskTensor&I` ✓,
  periodic self-adjoint ✓, albedo ✓, **white `AngularAverageOperator&I`**
  (`apply_transpose` NOT advertised — self-adjoint only under `|Ω·n|·w`; waits
  on O.2's metric).
**⭐ CO-LAND DECISION (2026-06-03, user-approved).** Grounding revealed the
plan's O.4a.2/O.4a.3/O.4a.4 split is IDEALIZED — the boundary semantics are
**coupled**: the matvec EMITS the boundary slot (defect today), the solver
CONSUMES it, and the boundary *unknown structure* flips (inflow goes from
"recomputed each matvec by the keystone @519" to an independent unknown driven
by the `ψ.inflow − B·ψ.outflow` consistency residual — a solver-level change).
**Also:** the must-stay-green set is NOT vacuum-only as first assumed — Gate 1.1
(`sphere_GL4_reflective`/`cyl_LS4_reflective`, matvec-level per-ordinate
flat-flux), Resolution-A decomposition (`bc=reflective`), and the curvilinear
streaming-equilibrium are all REFLECTIVE + matvec-level; they test the current
*coupled* `(L+C)` and MUST MIGRATE to the extracted `(L+C−B)` (retirement = test
migration), not "stay green unchanged". So O.4a.2/3/4 are **co-landed as one
flip** with the whole-trace `B` as a safe additive precursor:
- **Commit 1 (additive, no behavior change):** `SNBoundaryOperator` (new module
  `orpheus/sn/boundary_operator.py`) — the whole-trace `B` as a BOUNDARY-block
  leaf on `TimedFullField` (zero bulk; `boundary.face_view(face) =
  sn_mesh.bc_<face>.apply(ψ.boundary.face_view(face))` per trace face),
  `block_role = BlockRole.BOUNDARY`, `domain = codomain = sn_mesh.trace`,
  `capabilities = APPLY ∪ (APPLY_TRANSPOSE iff all per-face laws advertise it)`.
  Gate: per-face-restricted `B.apply ≡ legacy bc.apply` (pins the face→law
  wiring), zero bulk, role/domain/capabilities. Nothing consumes `B` yet.
- **Commit 2 (the coupled flip):** bare `L_full` matvec (delete keystone @519/
  twin) + single-pass sweep + solver inflow-unknown + `−B` consistency residual,
  co-landed, + MIGRATE the reflective matvec gates to `(L+C−B)`. Gate: vacuum
  bit-identical (apply + end-to-end) + reflective convergence-equivalence
  (MMS / `Q/Σ_t` / `k_∞`). Subsumes the O.4a.2(matvec)/O.4a.3(sweep)/O.4a.4(solver)
  mechanics below (kept as the mechanics reference).

The mechanics referenced by Commit 2 (formerly the separate O.4a.2/3/4 bullets):
- **[matvec mechanics]** In `_compute_LpC` (386-593) +
  `_compute_decomposition` (595-884): **delete 521/795**; read the backward-sweep
  seed from `ψ.boundary.inflow` (not `bc_outer.apply`); read the inner inflow from
  `ψ.boundary.inflow` (slab, was 456/716); boundary output = RAW outflow (drop the
  `− bc_estimate` at 576-578/859-861); KEEP the curvilinear pole Carlson seed
  (454/714 — regularity, not BC). Wire the boundary identities + `−B` at the
  OperatorSum level. **Gate:** vacuum matvec bit-identical; reflective matvec
  convergence-equiv; Gate 1.1, Resolution-A bit-exact (vacuum), L1 MMS slopes.
- **O.4a.3 — bare `L_full` (solve/sweep).** In `sweep.py:transport_sweep` (99):
  retire the entry `bc_*.apply` (473-474/488); persist-back RAW outflow (585/587,
  no BC); inflow from the seeded `boundary_buf` inflow slots. **Gate:** the sweep
  given a fixed `ψ.boundary.inflow` is a single pass = same bulk+outflow as before
  for a consistent inflow.
- **O.4a.4 — outer-solve rewire + boundary typing (`solver.py` + `iteration.py`).**
  `q.boundary → BoundarySource`; retire the SI `rhs.boundary = self._boundary_flux`
  seeding (572) — inflow becomes a solved unknown in `ψ.boundary.inflow`; SI adds
  the `Bψ` term to the source; the outer loop drives `r_inflow → 0`. Apply the
  honest boundary types. **Gate:** reflective/albedo/white/periodic converged ψ
  matches MMS / `Q/Σ_t` / homogeneous `k_∞` to solver tol (convergence-equiv);
  vacuum bit-identical; the `old_vs_new` drift tripwire (< 2×solver_tol).
- **O.4b — 2-D Cartesian (SEPARATE design).** `_apply_2d_cartesian` (1333-1463)
  fills the cell-centre proxy from BCs + forces the boundary output to zero (1456).
  Extraction requires making `face_view` an ACTIVE trace (read `ψ.boundary.inflow`,
  emit the defect) — NOT a mechanical port (prior active-trace attempt missed `k_∞`
  by ~10%; map Area 3). Highest-risk; pin with 2-D `Q/Σ_t` (L0) + 2-D MMS slope +
  the bulk-ψ drift tripwire + assert the new boundary residual → 0 at convergence.
  **The 2-D adjoint (Gate 1.3 2-D id) lands HERE, not at O.2.**

### Then O.2 (adjoint, on the bare shape)
`StreamingOperator.apply_transpose` (reverse sweep) + populate
`TraceSpace.inner_product_weights = |Ω·n|·w_n` (trace_space.py:276 — `omega_dot_n`
@309 + `quadrature.weights` both in scope; `_AdjointOperator` already reads it
@630-632) → the white-BC adjoint becomes free + the boundary-block reciprocity
becomes correct. Gate 1.3 flips green (sphere+cyl+slab); 2-D id stays xfail until
O.4b.

**O.2 carry-forward — composers must DERIVE block_role from operands (NOT stamp):**
1. `OperatorSum`/`ScaledOperator` derive `block_role` from operands (e.g. BULK if
   all-BULK; BOUNDARY if all-BOUNDARY) — and **RETIRE the hardcoded
   `InvertibleOperator.__init__` FULL stamp** (twin-path avoidance, from the O.1
   review).
2. `_AdjointOperator` must propagate `block_role` (pin `(L).H is FullOperator`,
   `B.H is BoundaryOperator`).
3. **`realize_recursively` composed BCs (`boundary_realize.py:206-215`) carry NO
   role today** (bare `OperatorSum` → `None`) — the γ elegance-enforcer NIT. When
   O.2 lands sum-role derivation, a mixed-BC `B` (`0.3·Reflective + 0.7·White`)
   auto-derives BOUNDARY from its all-BOUNDARY leaves. **Do NOT stamp the composer
   in `realize_recursively`** — derive it, else a twin path with O.2's derivation.
   No live consumer reads composed-BC role until O.2's `OperatorSum.apply`
   role-dispatch, so this is correctly deferred (not a γ bug).

### Open-Q resolutions (test-architect §7, from map Area 6)
1. metric @ `TraceSpace.from_mesh_and_quadrature` (276); both factors available;
   tangentials → 0 (positive-SEMI-definite — exclude tangentials from the pos-def
   test).
2. `B*` free via `(A+B).H` for reflective/vacuum/periodic/albedo; **white** needs
   the `|Ω·n|·w` metric first (couples to O.2) — the one non-free case.
3. 2-D adjoint at O.4b (forward residual must exist first).
4. non-vacuum MMS: new dataclass in `derivations/continuous/mms/sn.py`
   (`(A(x)+μB(x))/W`, `c<1`, `A>0`, derived `external_source`, prescribed inflow
   `γ₋ψ`, `ProblemSpec` BC `prescribed`). This is an **O.3b** item; O.4 verifies
   via `Q/Σ_t` + `k_∞` + existing vacuum MMS, adding the non-vacuum MMS at O.3b.

### Equivalence + gates (test-architect §1)
O.4a.1 = value-identical (`B.apply` ≡ legacy `bc.apply`); O.4a.2/.3/.4 =
convergence-equivalence for non-vacuum (iterates differ; converged ψ matches an
EXTERNAL reference under iter×cond×ULP), bit-identical for vacuum (B=0).
Must-stay-green every commit: operator-algebra-core, Gate 1.1, L1 curvilinear+2-D
MMS slopes, Resolution-A bit-exact, sentinel (no -O). **NO `continuous_get`**
(fixture bug #212 — use direct MMS / `Q/Σ_t` / `k_∞`).

### Risk ranking (do 1-D fully, verify, THEN 2-D)
O.4b 2-D ADD (highest — new compute, the ~10% `k_∞` cautionary tale) >
O.4a.2 reflective-curvilinear matvec extraction (ERR-006/026-prone seed) >
O.4a.4 solver rewire > O.4a.1 `B` leaf + alias (mechanical). A fresh
`test-architect` pass for the O.4a test CODE is optional — the verification
PLAN above + the referenced memo are sufficient to start O.4a.1.

---

## 0. Pickup checklist (read first)  *(historical stub — superseded by ⭐ CURRENT EXECUTION above)*

If you are picking this plan up in a fresh session:

1. **Read this plan top-to-bottom.** No section is optional.
2. **Verify Depth B + Wave T completed**: `git log --oneline -30` should show the D-A through D-K commits AND the Wave T commits (T.1–T.5; see `wave_t_tensor_network.md` §6 for the substep ledger).
3. **Read the GitHub issue**: `gh issue view 208` — the original ask + the post-D-H.1b conversation that filed it.
4. **Read these grand-report sections**:
   - §5.5 — Field hierarchy (the Flux / Source / Residual roles that need typed treatment).
   - §5.7 — Operator hierarchy (the Bulk / Full / Boundary classification at the type system).
   - §16A.10 (lines 3142–3197) — BC as tensor network; the boundary-only operator's algebraic shape.
   - §32.4–32.6 — Space / Field / Operator specs that Wave O instantiates the missing Protocol surface for.
5. **Read the existing L1 operator algebra**:
   - `orpheus/numerics/operator.py` — `LinearOperator` Protocol, `OperatorSum`, `OperatorProduct`, capability dispatch.
   - `orpheus/sn/operator.py` — the SN operator leaves (`StreamingOperator`, `CollisionOperator`, `InvertibleOperator`).
   - `orpheus/geometry/boundary.py` — the boundary-realizer family (`SpecularBoundaryOperator`, `VacuumBoundaryOperator`, etc.).
6. **Read `[[project_wave_o_operator_algebra]]`** memory and `[[project_transport_state_container]]` memory.
7. **Pick up at the leftmost incomplete step in §6** below.

The `[[feedback_no_method_implementer_for_surgical_carves]]` rule applies: this is the main agent's work with turn-by-turn user steering. Do NOT batch via method-implementer.

---

## 1. The principle

**Operators on `(BulkSpace ⊕ BoundarySpace)` have a natural 2×2 block structure that the type system currently erases.** Today, every operator inherits from a single `LinearOperator` Protocol regardless of which blocks it activates. This conflates three distinct algebraic shapes:

- **Bulk-only** operators (`C` = collision, `S` = scattering, `F` = fission): only the `L_bb` block is non-zero. No boundary action. Action: `(bulk, boundary) → (L_bb·bulk, 0)`.
- **Full** operators (`L` = streaming): all four blocks live. Reads upstream boundary, produces downstream boundary. Action: `(bulk, boundary) → (L_bb·bulk + L_bs·boundary, L_sb·bulk + L_ss·boundary)`.
- **Boundary-only** operators (`BoundaryRealizer` family — specular, vacuum, white, albedo, periodic): only the `L_ss` block is non-zero. Action: `(bulk, boundary) → (0, L_ss·boundary)`.

**Asymmetry 2 — Source / Residual are full fields, not bulk-only.** The deeper insight (D-H.1b conversation, 2026-05-28): a source has both a bulk part (volumetric `Q(r, Ω, g)`) AND a boundary part (prescribed incoming flux `ψ|_{Γ_-}`). Vacuum BCs make the boundary part zero, but MMS with non-trivial inflow wants typed boundary sources. Similarly, residuals from `(L+C)ψ − q` have both bulk and boundary parts — the streaming matvec already produces both today but stuffs them into an `AngularFlux + BoundaryFlux` pair without a typed wrapper.

So `AngularSource`, `AngularResidual` should be **TimedFullField-shaped types** (with explicit bulk + boundary members), not bulk-only `Field` types. The current implicit-zero-boundary pattern shipped in D-H.1b is a placeholder; Wave O lifts it to typed.

### 1.1 What Wave O closes

| Current state (post-Depth-B) | Post-Wave-O state |
|---|---|
| All operators inherit `LinearOperator` regardless of block structure | `BulkOperator`, `FullOperator`, `BoundaryOperator` Protocols at L1 distinguish the three algebraic shapes |
| `OperatorSum.apply` calls every leaf with the same `TimedFullField`; bulk-only leaves return implicit-zero boundary | `OperatorSum.apply` dispatches by Protocol: bulk-only leaves see `bulk`; full leaves see `(bulk, boundary)`; boundary-only leaves see `boundary` |
| `AngularSource` and `AngularResidual` don't exist as distinct types — `AngularFlux` is the shape container for all three roles | `TimedFullAngularSource`, `TimedFullAngularResidual` ship as TimedFullField composites; the `(role × storage) = 12-cell` matrix from Issue #205 becomes machine-checkable |
| MMS machinery emits `AngularFlux` as source; the boundary part is implicit-zero | MMS emits typed `TimedFullAngularSource` whose boundary part is non-zero for non-vacuum BCs |
| `(L+C).apply(ψ)` returns `TimedFullField` with both bulk and boundary residual; consumers must know the type erasure | The return type IS `TimedFullAngularResidual`; the operator's codomain explicitly types the algebraic structure |
| Gate 1.3 (apply ↔ apply_transpose reciprocity) is xfail-strict pending Wave O | Gate 1.3 lands: every operator's adjoint is built analytically by Protocol-dispatched `.H` propagation through `OperatorSum` |

### 1.2 The "two architectural asymmetries" framing

The user's diagnosis from the D-H.1b conversation (2026-05-28, recorded in `[[project_wave_o_operator_algebra]]`):

> **Asymmetry 1** — operators have implicit bulk/boundary character.  No type distinguishes their character.
>
> **Asymmetry 2** — Source / Residual are full fields, not bulk-only.  Sources and residuals have bulk AND boundary parts.

Both are dimensional-typing sins (Issue #205 scope). Wave O makes them both machine-checkable.

### 1.3 Dependency on Wave T

Wave O sits on top of Wave T (`wave_t_tensor_network.md`). Wave T factors operators into tensor products `A = A_x ⊗ A_ω ⊗ A_g`. Wave O classifies each TENSOR FACTOR (and each resulting `TensorProductOperator` / `SumOfTensorProductsOperator`) as Bulk / Full / Boundary.

Specifically:
- **Boundary realizers** (T.1: vacuum, periodic, white, albedo as `K ⊗ I ⊗ I` tensor networks) are `BoundaryOperator` instances.
- **Fission** (T.2: `F = χ ⊗ νΣ_f`) is a `BulkOperator`.
- **Scattering** (T.3: `Σ_ℓ Σ_ℓ ⊗ A_ℓ ⊗ G_ℓ`) is a `BulkOperator` (no boundary action; the moment-folding is bulk-only).
- **Streaming** (T.4: `L_spatial + L_angular_redist`) is a `FullOperator` — its `L_spatial` factor reads/writes face traces, its `L_angular_redist` factor is bulk-only.

Wave O cannot land before Wave T because the Protocol classification needs the factored shape to type each factor correctly. Doing Wave O first would force the classification on the un-factored flat-axis legacy operators, which is the wrong shape.

### 1.4 Expected complications

- **`OperatorSum.apply` dispatch** currently consumes a single `TimedFullField` and routes to every leaf. Post-Wave-O, the dispatch must split the composite into `(bulk, boundary)` and route to bulk-only / full / boundary-only leaves per their Protocol — without breaking the algebra closure laws `(A + B).H = A.H + B.H` and `(A ∘ B).H = B.H ∘ A.H`.
- **The codomain space** of a full operator is `BulkSpace ⊕ BoundarySpace`. Wave O may need `DirectSumSpace` (deferred to P3.6 in Depth B) earlier than planned, or a lightweight ad-hoc shape until P3.6 lands. The `[[feedback_unify_after_two_instances]]` rule applies: don't ship `DirectSumSpace` until P3.6 also wants it (flux ⊕ precursors).
- **Backward compatibility with implicit-zero**: D-H.1b's implicit-zero-boundary pattern works because `TimedFullField` defaults boundary to zero. Wave O can keep the implicit-zero shorthand at the leaf level while adding typed `AngularSource` / `AngularResidual` at the composition level. Both should round-trip.
- **`apply_transpose` Protocol completeness**: Gate 1.3's xfail-strict status reflects that the analytical-adjoint propagation through `OperatorSum.H` is not yet wired for the operator algebra. Wave O lands this: each Bulk / Full / Boundary Protocol declares its own `.H` analytic, and `OperatorSum.H` composes by Protocol.

These complications are deferred to "discover during execution" — the architectural commitment to Protocol-classified operators is non-negotiable; the concrete shape of each substep refines when the substep starts.

---

## 2. Dependencies

Wave O cannot start until **both Depth B AND Wave T complete**. Specifically:

- Depth B D-A through D-K: the typed `Field` / `TimedFullField` substrate.
- Wave T T.1 through T.5: every production operator factored as `TensorProductOperator` / `SumOfTensorProductsOperator`. The Wave O Protocol classification names a factor's algebraic role; that role is only well-defined once factoring is in place.

Out-of-scope for Wave O:
- **`DirectSumSpace`** — deferred to P3.6 per Depth B §11.1 invariant #13 and `[[feedback_unify_after_two_instances]]`. Wave O may surface a second consumer (the operator codomain `BulkSpace ⊕ BoundarySpace`); if so, the unification with P3.6's `flux ⊕ precursors` lands at P3.6.
- **Issue #201** — `split AngularFlux from AngularSource / AngularResidual` — the role-typing slice. Wave O lands its bulk-side; the angular-storage slice of the 12-cell matrix lands at #201 close-out.
- **Issue #205** — the full 12-cell matrix. Wave O is a foundational layer; #205 closes when every cell ships its concrete type.
- **MoC and CP adoption** — the operator Protocols ship at L1 and any method can consume them. MoC / CP migrations land in their own waves.

---

## 3. Substep ledger

Each substep has a single named gate and a single commit boundary (1–3 commits per substep).

| Step | Scope | Gate | Status |
|---|---|---|---|
| **O.1** | `BulkOperator` / `FullOperator` / `BoundaryOperator` Protocols at L1 in `orpheus/numerics/operator.py`. Each Protocol declares: `apply(TimedFullField) → TimedFullField` semantics; capability set; analytic `.H` adjoint. | All existing operators retain runtime behaviour; capabilities API gains the Protocol classification. | PENDING |
| **O.2** | `OperatorSum.apply` / `.H` dispatch by Protocol: bulk-only leaves receive `bulk` only; full leaves receive `(bulk, boundary)`; boundary-only leaves receive `boundary` only. Algebra closure laws `(A + B).H = A.H + B.H` and `(A ∘ B).H = B.H ∘ A.H` hold by construction. | `test_phase_c_gates.py` Gate 1.3 (apply ↔ apply_transpose reciprocity) flips from xfail-strict to passing across slab + sphere + cylinder + 2-D Cartesian. | PENDING |
| **O.3** | `TimedFullAngularSource` and `TimedFullAngularResidual` typed wrappers — TimedFullField composites with explicit `bulk + boundary` members. MMS machinery emits typed sources for non-vacuum BCs. | New L0 invariants: source's `boundary` part is non-zero IFF the BC is non-vacuum; residual's `boundary` part equals the streaming matvec's face residual exactly. | PENDING |
| **O.4** | Migrate existing consumers: `solve_sn_fixed_source` accepts `TimedFullAngularSource`; `(L+C).apply` returns `TimedFullAngularResidual`; the implicit-zero-boundary shorthand from D-H.1b retires at the consumer boundary (not at the operator leaf). | All 12-file operator-algebra-core tests stay green; cardinality preserved; new test pin for non-vacuum boundary-source MMS lands. | PENDING |
| **O.5** | Retire the implicit-zero-boundary shortcuts shipped in D-H.1b. Each bulk-only leaf's `apply(TimedFullField) → TimedFullField` returns its result via the Protocol-classified `BulkOperator` API; the `OperatorSum.apply` dispatch fills the implicit-zero boundary at the algebra level, not at the leaf. | Final clean-up: zero `# implicit-zero boundary` comments in production; the typed contract is universal. | PENDING |

Each substep follows the veto pattern (Lessons L12 + L13 + L17 + L18 + L20) established by Depth B Wave D-I / D-J.

---

## 4. Verification strategy

### 4.1 New tests added by Wave O

- `tests/numerics/test_operator_protocols.py` — Protocol conformance (BulkOperator / FullOperator / BoundaryOperator); each existing operator advertises the correct Protocol.
- `tests/numerics/test_operator_sum_protocol_dispatch.py` — `OperatorSum.apply` dispatches by Protocol; bulk + boundary parts route correctly.
- `tests/numerics/test_typed_source_residual.py` — `TimedFullAngularSource` / `TimedFullAngularResidual` invariants (boundary part non-zero iff non-vacuum BC).
- `tests/sn/test_mms_non_vacuum_boundary.py` — MMS with non-vacuum BC; verifies the boundary part of the typed source matches the analytical inflow.

### 4.2 Tests that MUST stay green at every commit

- 12-file operator-algebra-core gate (the post-D-J baseline: 286 passed, 13 xfailed).
- Resolution A bit-exact decomposition gate (`test_streaming_operator_decomposition.py`).
- L1 MMS pins (Test 3.1 in `test_2d_l2_matvec_correctness.py`; curvilinear MMS in `test_mms_curvilinear.py`).
- L0 streaming-equilibrium anchor (`tests/sn/test_phase_c_gates.py` Gate 1.1).
- Gate 1.3 (apply ↔ apply_transpose reciprocity) — currently xfail-strict; FLIPS GREEN at O.2 (the load-bearing Wave O deliverable).

### 4.3 Bit-identity judgment per step

- **O.1**: Protocol classification only — runtime behaviour identical. Bit-identity required.
- **O.2**: Dispatch reorganisation — bit-identity expected on existing leaves; if FP-non-associativity drift surfaces, the three-criteria rule of `vv-principles` §"Bit-identity vs principled-equivalence" applies.
- **O.3 / O.4 / O.5**: Typed-wrapper introduction at the consumer boundary; the bulk part of every result MUST be bit-identical to pre-Wave-O. The boundary part is a NEW component (previously implicit-zero); its non-zero values are NEW invariants pinned by O.3's tests.

---

## 5. Exit Route

When Wave O completes:

- **Gate 1.3 flips green** — operator-algebra reciprocity is gated to round-off across all geometries.
- **Issue #208 closes.**
- **Issue #205 becomes much closer to closeable** — the operator-side typing lands; the field-side (`AngularFlux` / `ScalarFlux` / `HarmonicMomentField` split into Flux / Source / Residual) is the remaining work (#201 + scattered other tickets).
- **P3.4 (Problem/Solver split) becomes implementable on top of the role-typed algebra.** `CriticalityProblem.loss: FullOperator` and `CriticalityProblem.fission: BulkOperator` become machine-checkable type signatures.

The branch state at Wave O completion: `refactor/moment-space-and-layering` ready for P3.4 sequencing.

---

## 6. Cross-references

- **Parent plan:** `.claude/plans/moment_space_and_layering_plan.md` — Phase 3 sequencing (post-2026-05-30 update: `Depth B → Wave T → Wave O → P3.4 → P3.6`).
- **Depth B plan:** `.claude/plans/depth_b_field_on_function_space.md` — the typed Field substrate Wave O extends; §11.3 documents the post-Wave-T → Wave O ordering.
- **Wave T plan:** `.claude/plans/wave_t_tensor_network.md` — the tensor-product factoring Wave O classifies.
- **Project memory:** `[[project_wave_o_operator_algebra]]` — the original 2026-05-28 decision context.
- **Issue #208:** `gh issue view 208` — the ticket.
- **Issue #205:** `gh issue view 205` — the full cross-method field architecture matrix Wave O lands the operator-side of.
- **Issue #201:** `gh issue view 201` — the angular-storage role-typing slice (Wave O's downstream).
- **Grand report:** `.claude/plans/neutron_transport_grand_report_v3.md` §5.5 (Field hierarchy), §5.7 (Operator hierarchy), §32.4–32.6 (specs).
- **Lessons:** L17 (convention crosswalk before carve — applies to the Protocol dispatch carve at O.2), L18 (Pattern 7 producer-side normalisation — applies to the implicit-zero-boundary collapse at O.5).
