# Stencil-assembly + DSA roadmap — k-estimator → spatial reification → DSA (#2)

**Status: IN EXECUTION — Phase 1 opened 2026-07-03 on `refactor/k-estimator-unification`
(off main @ `3a19133`); rulings R1–R6 COLLECTED (see the Rulings section). Campaign
chain amended per R1.**

Campaign chain: **Phase 1** k-estimator unification (#259 + #291) → **Phase 2** spatial
substrate promotion + the ASSEMBLY third mode (#272 + #158 + #253, user ruling below) →
**Phase 2.5 (per R1)** the #280 orientation×kernel walk unification, its own
branch/campaign → **Phase 3** DSA for SN (#2, folding #200's Krylov posture; wires into
the UNIFIED walk). Each phase is a separate branch + ff-merge, in order; later phases
consume earlier ones.

## Provenance (how this plan was decided)

Follows directly from campaign #290 (merged @ `3a19133`; see
`.claude/plans/diffusion_integration_290.md`). The user chose the sequence
warm-up → spatial extraction → DSA, then supplied the architectural key (below). The
recon facts in this file were verified against the tree at `3a19133`.

## THE USER RULING (2026-07-03) — assembly is the missing third consumption mode

> DD and LD carry only the closure for forward substitution and Krylov application
> because we don't yet have a matrix representation for them; a matrix representation
> demands the DD/LD stencil implementation, just as diffusion's does. Extracting the
> spatial discretization therefore ALSO builds the assembly machinery.

Verified against the tree — the claim is true and *conservative*:

- `orpheus/numerics/operator.py:751` — `as_matrix` is **basis-vector probing** (column j
  = apply to e_j; the "functor out of the operator category"; `max_dimension=4096`).
  NO operator in the codebase has stencil assembly. SN's only matrix route is O(N)
  sweeps — the historical 640 s matvec-twin gate (#236) is this cost made visible.
  Diffusion's `MatrixInverseOperator(FlattenedOperator(A, template))` also
  materializes by probing (tolerable only because 1-D N is small).
- The per-cell stencil coefficients ALREADY EXIST as data: `orpheus/sn/spatial/
  cell_balance.py` (`CellBalanceTerms`, `cell_balance_terms`) for DD;
  `linear_discontinuous.py` (`_LDCellTerms`) + `_ubld.py` (`assemble_ubld`,
  `assemble_inflow_axis`, `per_cell_solve` — LD already assembles per-cell block
  systems) for LD. The closure IS the stencil, consumed in-flight.

**The reframe that makes it safe (load-bearing design rule):** assembly is NOT a new
discretization — it is the third consumption mode of the ONE closure algebra:

| Mode | Consumer | Walk |
|---|---|---|
| `solve` | sweep = forward/triangular substitution | cells in walk order |
| `apply` | Krylov matvec = action | cells in any order |
| `assemble` | **NEW** — emit the SAME coefficients as (row, col, value) | cells, scatter into global sparse |

**Twin-path guardrail (Pattern 2; the Phase-F lesson):** `assemble` MUST consume the
same coefficient source (`CellBalanceTerms` / `_LDCellTerms` / closure config) that
`solve`/`apply` consume. NEVER a parallel stencil implementation. Pinned by permanent
equivalence gates: `assembled @ x ≡ apply(x)` (ULP-tight) and per-ordinate
`triangular_solve(assembled, q) ≡ sweep solve` on Cartesian.

**What one assembly layer buys (why this ordering is right):**
1. Sparse matrix representation for SN L(+C) — replaces O(N)-probe twins in tests
   (perf), enables sparse-direct/preconditioner consumers (#200), spectral analysis.
2. Diffusion assembly upgrades from probed-dense to direct-sparse (same machinery;
   optional sparse-LU resolvent).
3. **#284 becomes provable at the object level**: sweep solve ≡ triangular solve of
   the assembled matrix in walk order (Cartesian: expected EXACT). The sphere's #282
   lagged pole seed becomes VISIBLE as non-triangularity of the assembled operator in
   walk order — a characterization gate, not a fix, in this campaign.
4. **Consistent DSA by derivation**: the low-order operator can be DERIVED as the
   discrete moment reduction of the ASSEMBLED SN operator (Alcouffe's consistency,
   computed rather than hand-imposed) — Phase 3's foundation.

## Execution mode + compaction protocol

**Surgical mode throughout** (main agent writes; NO `method-implementer`; other agents
dispatched freely — the #290 pattern). **R6 RULED**: Phase 2a's mechanical relocation
is ALSO full-surgical — the delegation exception is struck.

Compaction points marked ⏸ Cn. Before each: (1) commit everything; (2) append a
`## CHECKPOINT` / per-phase `**STATUS:**` block here (done @ hash / deviations /
rulings); (3) tell the user it is safe to `/compact`. Re-anchor after compaction from
THIS file + `git log` (trust git over any summary — process-discipline rule).

Canonical gates every phase (the #290 set): serial walls
`.venv/bin/python -O -m pytest <paths> -m "not slow" -p no:xdist -p no:cacheprovider
--timeout=300`; full `-m "not slow"` tree before each merge (0 reds); CLI
`npx pyright orpheus/` = 1 (the #288 residual — diffusion + the six folders at 0);
sphinx `-W` exit 0 (generate_rst only if derivations change); `python -m
tests._harness.audit` exit 0; demo k = 1.022173 where diffusion is touched;
mutation-verify every new gate's teeth (in-process monkeypatch, never git-checkout on
uncommitted files; pytest-plugin sentinels for -O bare-assert hazards — Mode 8/11).

---

## Phase 1 — k-estimator unification (#259 root + #291 symptom) [~1 session]

The diffusion P5 discipline, transplanted while hot: derive the k-update denominator
THROUGH the loss operator (`⟨1,(Aψ)_bulk⟩_V` = absorption + leakage by the column-sum
theorem + telescoping of the conservative divergence) so no term CAN be forgotten —
then unify the fragmented copies.

**Inventory (verified @ `3a19133`) — six `compute_keff`/`compute_production_rate`
sites:**
- `orpheus/diffusion/solver.py:296,306` — the CLEAN pattern (P/⟨1,(Aψ)_bulk⟩;
  IntegratedReactionRate #270 arm). The template.
- `orpheus/sn/solver.py:941,966` — **#291: omits leakage** (biased k-update for
  vacuum-bounded eigenproblems).
- `orpheus/cp/solver.py:689`, `orpheus/moc/core.py:211` — hand-rolled copies.
- `orpheus/numerics/eigenvalue.py:123,163` — the protocol declarations.
- `orpheus/numerics/iteration.py:1217,1233` — per #259, the KEigenvalue
  estimator-injection seam is DEAD in production — retirement candidate (3-search
  audit; rewire-not-delete for behavioral tests).

**Steps:**
1. **P0 dispatch — explorer**: per-site map (who calls each estimator in production;
   which tests pin each; whether SN's L1 vacuum anchors are sensitive to the #291
   bias or reflective-dominated). Output: rewire table + expected-shift table.
2. **Characterize #291 FIRST** (before fixing): run a vacuum-bounded SN anchor, log
   reported-k vs `direct_eigenvalue`-style ground truth (or balance-identity defect).
   Decide bit-identical vs principled-re-baseline per vv-principles BEFORE the carve.
3. **The carve**: one estimator discipline — production via the typed
   `IntegratedReactionRate`; denominator derived through each method's loss operator
   (per-method wiring, ONE shared spelling where the algebra permits — extract the
   helper only if ≥2 sites are literally isomorphic; no premature abstraction).
   SN gains the leakage term structurally (#291). CP/MoC rewired to the same
   discipline (mark OPTIONAL — drop to a follow-up if friction; they are outside the
   six-folder focus).
4. **Retire** the dead iteration.py estimator-injection seam (#259) — 3-search audit.
5. **Gates**: per-method balance identity `P/k = absorption + leakage`; SN
   reported-k ≡ eigen-solve k on a vacuum anchor (tolerance from step 2); teeth: a
   leakage-drop mutation must go RED on a vacuum case and stay GREEN on reflective
   (that asymmetry IS the #291 catcher).
6. Close #291 + #259 via trailers; comment #270 (CP/MoC arms status).

⏸ **C1** (small phase — compaction optional; checkpoint block regardless).

---

## Phase 2 — spatial substrate promotion + the assembly mode (#272 + #158 + #253)

### 2-P0 — dispatches (before any move)
- **explorer**: the promotion split map for `orpheus/sn/spatial/` — for each module,
  method-generic vs SN-sweep-specific (see draft split below), with the full 3-search
  blast radius (graph callers + text grep code/tests/docs + direct constructors).
- **test-architect** (MUST trigger — cross-subsystem carve): the bit-identity gate
  plan for the relocation + the equivalence-gate spec for the assembly mode
  (assembled≡apply, triangular≡sweep, probed≡assembled).
- **cross-domain-attacker**: the assembly framing — global DOF numbering as the flat
  functor; block structure (per-ordinate triangular streaming blocks, moment-coupled
  scattering); whether the (row,col,value) emission wants its own algebraic type.

### 2a — Relocation: `sn/spatial/` → `transport/spatial/` (discharges #272)
Draft split (P0 verifies; verified inventory @ `3a19133`):
- **PROMOTE** (method-generic trial-space/closure layer): `scheme.py`
  (`DiscretizationSchemeBase` + registry), `diamond.py` (`DiamondDifference`),
  `linear_discontinuous.py` (`LinearDiscontinuous`, `_LDCellTerms`), `_ubld.py`
  (per-cell mass/assembly kernels), `cell_balance.py` (`CellBalanceTerms`),
  the moment-count policy (#253 — single-source cell `per_axis^d` / face
  `per_axis^{d-1}` in ONE place at promotion time).
- **STAY in sn/** (sweep-walk + angular machinery, not spatial trial space):
  `pole_angular_closure.py`, `psi_half_angle_seed.py`, `sweep_cache.py`, `scan.py`;
  `pairing.py` per P0's call.
- Bit-identical relocation (git mv + import rewires; NO behavior). Full 3-search
  retirement audit on the old paths; Sphinx xref sweep (unresolved :mod:/:class:
  refs render silently — grep the built HTML).
- #158's ask (cell updates as first-class concrete types) is largely ALREADY TRUE
  (registry-keyed classes exist) — sharpen contracts/docstrings at the new home,
  close or comment #158 honestly per what's left.

### 2b — The assembly mode (the reification)
1. **Numerics home**: an assembled-sparse operator (scipy CSR/COO carrier) in
   `orpheus/numerics/` conforming to `LinearOperator` (apply = sparse matvec;
   `as_matrix` = densify). **R2 RULED**: structural assembly is a separate
   `assemble()` → sparse assembled operator; `as_matrix` keeps its dense
   Mat-functor semantics and DELEGATES to densified assembly when available
   (probing remains the fallback for assembly-less operators; #213 not blocked on).
2. **Global DOF numbering** = the existing flat layouts (`FullField` ravel order,
   `FlattenedOperator`'s template contract, face layouts) — do NOT mint a second
   numbering; the flat functor already exists.
3. **Emitters**: DD via `CellBalanceTerms` → per-cell (row,col,value); LD via
   `_LDCellTerms`/`assemble_ubld` blocks → block scatter (moment-blocked DOFs — the
   #253 single-sourcing is a genuine prerequisite here, not a rider). Scope:
   **Cartesian (slab + 2-D) only**; curvilinear assembly is OUT (blocked on the
   #282 pole-seed structure — see the characterization gate below).
4. **Consumers wired**:
   - SN: per-ordinate `L(+C)` blocks assemble triangular-in-walk-order; scattering
     stays the K_iso/kernel form (its dense kernel exists via `full_scatter_kernel`).
   - Diffusion: `LeakageOperator`/family assembly replaces probing for the matrix
     path — gate bit-identical vs the probed matrix, THEN optionally switch the
     resolvent to sparse LU (perf; principled-equivalent gate).
   - Tests: replace O(N)-probe twin gates with assembled gates where they exist
     (keep ONE probed≡assembled pin per family as the permanent oracle — the
     fuller-view-oracle exception, decided EXPLICITLY per the retirement rule).
5. **Gates (the teeth of the whole campaign)**:
   - `assembled @ x ≡ apply(x)` per scheme × geometry × ng (ULP-tight).
   - Per-ordinate `scipy triangular solve(assembled L+C) ≡ sweep solve` (Cartesian
     DD + LD; bit-or-ULP per test-architect's spec) — **this discharges #284's
     contract question with an object-level proof.**
   - **#282 characterization**: assembling the spherical walk order and showing the
     pole-seed row breaks triangularity = the lag made visible. Record on #282 as
     evidence; the FIX stays in #282/#280's scope.
   - Mutation teeth: a coefficient sign-flip in the emitter must red BOTH the
     equivalence gate and (via the shared source) the sweep gates — proving there is
     ONE source. If only the equivalence gate reds, a twin path exists → stop, fix.

### 2c — Docs + close-out
- Theory: `loss_representations.rst` / `discrete_ordinates.rst` — the three-modes
  section (solve/apply/assemble as consumption modes of one closure algebra);
  `operator_algebra.rst` — `as_matrix` probing vs structural assembly; dev-history
  rows. Close/comment #272, #253, #158, #284; comment #282 with the
  characterization; comment #200 (blocks now assemblable).

⏸ **C2**.

---

## Phase 3 — DSA for SN (#2; folds #200's Krylov posture; runs AFTER the #280 walk
unification per R1 — the accelerator wires into the unified orientation×kernel walk)

### 3-P0 — dispatches (MUST, before the plan-of-record for this phase)
- **literature-researcher** — check `/Users/rodrigo/git/nuclear/ORPHEUS/scratch/
  literature/` FIRST (user maintains it; all NSE volumes local): Alcouffe 1977
  (consistent DSA for DD), Adams & Larsen 2002 review (fast iterative methods —
  the ρ = 0.2247c Fourier results + partial-consistency taxonomy), Adams & Martin
  1992 (M4S), Warsa–Wareing–Morel 2004 (Krylov-DSA robustness where SI+DSA
  degrades). If a paper is not local: report back, do NOT pivot to a secondary
  source without user approval.
- **test-architect** — the verification plan: FP-invariance gates, spectral-radius
  measurement protocol, mutation matrix. MUST encode the ⚠#215 / vv Mode-9 trap:
  every FP-invariance gate runs an ANISOTROPIC config (vacuum/heterogeneous) AND a
  diagonal cubature — never the isotropic reflective box alone.
- **cross-domain-attacker** — the restriction/prolongation pair as a Petrov–Galerkin
  frame on the angular axis (moment-0 analysis face, isotropic reconstruction face)
  — reuse the frame vocabulary; do not mint an ad-hoc R/P pair.

### 3a — Consistency by derivation (Phase 2's payoff)
Compute the discrete P1/diffusion moment reduction of the ASSEMBLED Cartesian DD
operator (slab first) and compare — as MATRICES — against the diffusion family's
assembly on the same mesh:
- If ≡ (or ≡ up to a characterized boundary row): the landed `LeakageOperator` IS the
  consistent partner — Alcouffe's result recovered computationally. Gate it, done.
- If ≢: **R4 RULED "decide from 3a data"** — bring the matrix diff to the user with
  its STRUCTURE characterized (boundary row only? interior stencil? a scaling?) and
  collect the derived-stencil vs partial-consistency ruling THEN. Do not pre-commit
  either way; both remain cheap once assembly exists.

### 3b — The accelerator
1. **R** (restriction): moment-0 of the typed SN residual (`evaluate_residual`,
   `sn/solver.py`) → diffusion source on the scalar composite; boundary arm = the
   ℓ=0 half-range moment under the SHARED `|Ω·n|w` metric (the #290 ruling-2 seam).
2. **Low-order solve**: the (3a) operator over `DiffusionMesh.from_material_mesh(
   sn_mesh)` — same axes/materials/BCs by construction. Correction-equation BCs:
   vacuum faces → Marshak (𝒜=0; the error problem has zero inflow — the #290
   vacuum law is EXACTLY right), reflective → reflective.
3. **P** (prolongation): isotropic correction onto the φ moments (+ the PG-frame
   shape from 3-P0).
4. **Wiring**: both postures through the existing iteration layer — SI+DSA
   (`SourceIteration` wrap) and Krylov-preconditioned (`KrylovAcceleration`; folds
   #200's intent). The Protocol trigger FIRES: declare `full_field_space` on
   `TransportMethod` (the recorded P7b deferral — first generic consumer arrives).
5. **Scope**: slab/Cartesian DD = arm 1. **LD arm 2 gated on measured need** —
   **R5 RULED "decide from measured table"**: measure LD's ρ under the arm-1
   low-order operator across the c→1 × optical-thickness sweep, present the table
   to the user, collect the build-LD-consistent ruling then (M4S / DG-diffusion
   territory if triggered) — deliberately demand-driven, per defer-until-
   consumption. Curvilinear DSA: OUT (blocked on #282; documented seam).

### 3c — Gates
- **FP-invariance (Mode 9, #215-proof)**: converged accelerated flux ≡ plain-SI
  fixed point to solver tolerance, on vacuum+heterogeneous anisotropic config AND a
  diagonal cubature; separately from any rate claim.
- **Spectral radius**: measured ρ vs the ~0.2247c Fourier bound on the classic
  homogeneous-infinite/periodic problem; degradation curve vs cell optical thickness
  (the consistency stress axis); Krylov iteration counts vs SI+DSA vs plain SI —
  the c→1 table IS the teaching artifact (demo + theory page).
- **Correction→0 gate** at convergence (the correctness-safe-by-construction pin).
- **R/P object gates** (Mode 12): conservation `⟨1, R r⟩ = ⟨1, r⟩`-family +
  adjoint-consistency of the pair; mutation-verified.
- **No-masking control**: a seeded transport-operator mutation must converge (fast)
  to the SAME wrong answer with DSA on — the accelerator must not launder bugs.
- ⚠ within-group solves must NOT route through the #215 σ_r-fold spelling.
- Full #290-style close-out: theory docs (`discrete_ordinates.rst` DSA section with
  the Fourier analysis + literature eq numbers; cross-refs to `diffusion_1d.rst`'s
  DSA-seam section), dev-history rows, `Closes #2`, dispositions for #200/#215,
  comment #22 (N-D interaction).

⏸ **C3** + campaign close-out.

---

## Rulings — COLLECTED 2026-07-03 (execution start)

- **R1 (#280 sequencing)**: **BETWEEN Phases 2 and 3** — #280's orientation×kernel
  walk unification runs as its own branch (Phase 2.5) after assembly lands and
  before DSA, so the accelerator wires into the UNIFIED walk. (Deviates from the
  recorded recommendation "after the campaign": the user chose the cleaner
  consumption order for DSA over campaign brevity.)
- **R2 (assembly spelling)**: **separate `assemble()`** → sparse assembled operator;
  `as_matrix` keeps its dense Mat-functor semantics and DELEGATES to densified
  assembly when available; probing stays the fallback for assembly-less operators.
  The parked #213 capability-morphism redesign is not blocked on.
- **R3 (Phase-1 CP/MoC arms)**: **include, drop on friction** — attempt the CP/MoC
  rewires after SN+diffusion land; on substrate friction, file a labeled follow-up
  instead and close #259 scoped to sn/numerics/diffusion.
- **R4 (3a outcome)**: **decide from 3a data** — deliberately deferred; bring the
  derived-vs-landed matrix diff (structure characterized) to the user at 3a.
- **R5 (LD arm criterion)**: **decide from measured table** — deliberately deferred;
  present the measured-ρ c→1 × optical-thickness sweep at Phase 3, rule then.
- **R6 (2a delegation)**: **full surgical** — the mechanical relocation is main-agent
  work too; no delegation exception.

## Issue map

| Phase | Drives | Touches / comments | Follow-ups it may file |
|---|---|---|---|
| 1 | #259, #291 | #270 (CP/MoC arms) | CP/MoC estimator rewire (if R3 defers) |
| 2 | #272, #253, (#158) | #284 (discharged), #282 (characterized), #200, #256 | curvilinear assembly (post-#282) |
| 2.5 (R1) | #280 | #282 (structure fix candidate) | — |
| 3 | #2 | #200 (posture folded), #215 (disposition), #22, #293 | LD low-order op (per R5 table), curvilinear DSA |
