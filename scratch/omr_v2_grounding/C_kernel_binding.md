# Grounding C — §I.3 (kernel ≠ operator) + Phase 1 vs HEAD of main

Reconciliation of `.claude/plans/orpheus-operator-machinery-report-v2.md` §I.3
(:71-102), corrections #13/#20 (:265/:272), III-A mis-factorings (:288-296) and
the Phase 1 table (:404-414) against **HEAD `7aae9bf1` (main, 2026-08-19)**.
Plan dated 2026-08-08. All `[M]` = measured this dispatch (command given).
Siblings: space layer = A, operator machinery (splitting/pencil) = B, flux
ontology = D — one-line handoffs only where noted.

## Headline

**The plan's §I.3 was written ≥6 weeks stale.** Every load-bearing §I.3 premise
about S/F/C dates to a tree state that ended at `bbe8a51d` (2026-06-26, #261):
the plan's v1 operator audit was "retained" into v2 unverified while the tree
had already landed (a) declared Optional `(domain, codomain)` on S/F/C, (b) the
`.kernel` exposure + `IntegralKernelOperator` **Protocol** (an explicit
counter-doctrine to correction #13), (c) the production `R∘Λ∘M` binding via
`frame.conjugate(Λ)` (`cc6d022e` 2026-06-24), and (d) the ratified
**frame-eigenbasis-ownership** ruling (`d30d4a68` 2026-06-25) that contradicts
the v2 amendment "the frame and measure live on the space's axes". What
survives of the diagnosis at HEAD: **apply-time carrier dispatch on S/F/C is
real and unchanged** (`[M]` `git log --since=2026-08-08 -- orpheus/transport/operators/`
→ **0 commits**), and no tightness gate exists. Post-plan, the closure/L story
moved substantially (#344 `loss_kernel_gauge`, closure-seam gates) in a shape
the plan does not anticipate.

## 1. Per-operator dispatch-state table

`[M]` sources: class bodies read at HEAD; dispatch sites via
`grep -rn singledispatch orpheus/transport/operators/`.

| Op | Class (file:line) | Declared spaces? | Apply-time parsing? | Binding site / measure |
|---|---|---|---|---|
| **S** | `ScatteringOperator` — `orpheus/transport/operators/scattering.py:356` | YES, Optional: `domain`/`codomain` properties return `full_field_space: FullFieldSpace \| None` (endomorphic; `None` on bare construction ⟹ guard skips). P4.5 W-D pattern, landed ≤ `bbe8a51d` 2026-06-26 | YES: `@singledispatchmethod _apply_impl` (scattering.py:1093), **4 arms** — `FullField→FullField`, `ScalarFlux→ScalarSourceSink`, `AngularFlux→AngularSourceSink`, `HarmonicMomentFlux→AngularSourceSink` — with honest per-carrier `@overload` static surface (#257 S8c) | Frame **owned by the operator**: `frame` `cached_property` (scattering.py:488) = `HarmonicFrame.from_galerkin(quadrature.angular_frame(L))`; `kernel` (scattering.py:599) = `frame.conjugate(LegendreMomentScattering)` = typed `R∘Λ∘M`; `full_scatter_kernel` (:656) = `R∘(Λ_{ℓ≥0}+N₂ₙ)∘M`. Measure enters via the frame's `DiscreteMeasure` |
| **F** | `FissionOperator` — `orpheus/transport/operators/fission.py:147` | YES, Optional: same W-D pattern (`full_field_space \| None`; F never enters a production `OperatorSum`, so the space is the "honest cross-method tag") | YES: `@singledispatchmethod` (fission.py:469), **3 arms** — `FullField→FullField`, `ScalarFlux→ScalarSourceSink`, **bare `np.ndarray`** (KEigenvalue/depletion/diffusion outer boundaries) | `kernel` property (fission.py:258) = rank-1 `TensorProductOperator` `RankOneOperator(χ, νΣf, axis=0) & IdentityOperator()`; adjoint = free dyad swap. Data single-sourced on `mat_xs` (`MaterialXSField`) |
| **C** | `MultiplicationOperator` — `orpheus/transport/operators/multiplication_operator.py:130` (`LinearOperator["FullField"]`) | YES, Optional: `space: FunctionSpace \| None = None` field; `domain`/`codomain` (:341/:357) return it; `#261: folded up from the retired CollisionOperator` | YES: `@singledispatchmethod` (:376), **2 arms** — `FullField` (with an internal `isinstance` family parse Angular/Scalar, #289 "parse loudly at the seam") and **bare `np.ndarray`** (meshless escape hatch, #276) | No integral structure (local/diagonal — NOT an `IntegralKernelOperator`); multiply delegated to a `DiagonalOperator` engine built once in `__post_init__`. **Zero frame imports** |
| **L(+C)** | `StreamingCollisionOperator` — `orpheus/sn/operators/streaming.py:445`, typed `OperatorSum["FullField","FullField",StreamingOperator,MultiplicationOperator]` | YES, fixed via the OperatorSum generic (FullField-endomorphic) | **NO** dispatch: `apply(psi: FullField) -> FullField`; no singledispatch in streaming.py (`[M]` grep) | `.solve` = WDD sweep (`.inverse()` → `SweepOperator`); matvec/adjoint/sweep = three actions of ONE operator (L21), all routed through the strategy from `default_for(self.sn_mesh)` (streaming.py:404); spatial closure bound at **mesh construction** (`SNMesh(scheme=…)`, default `DiamondDifference()`, `augmented_mesh.py:208,265-266`) |

Note the reconciliation the two facts force: the plan's cause-1 (measure
ambiguity, S/F) is **resolved** — the measure is bound at construction through
the operator-owned frame. The surviving apply-time dispatch on all three is the
plan's cause-2 (carrier shape) — and the tree treats it as a deliberate
cross-method surface, not a defect: S's `ScalarFlux` arm carries an explicit
in-code keep ruling ("Deliberately retained — a named-future-consumer surface"
for #205; scattering.py:1172), and C's ndarray arm is the #276 meshless
cross-model hatch.

## 2. Claim-by-claim verdicts

### §I.3 (plan :71-102)

| Claim | Verdict | Evidence |
|---|---|---|
| "`ScatteringOperator`, `FissionOperator`, `CollisionOperator` have no unique (domain, codomain) pair" | **CHANGED (stale at write)** | Declared Optional pair on S/F/C since ≤ 2026-06-26 (`[M]` `git log -S "Operator-algebra space metadata" -- orpheus/transport/operators/scattering.py` → only `bbe8a51d`). Nuance keeping the *spirit* alive: the declared pair describes only the composite arm; the other arms are genuinely different arrows (e.g. `ScalarFlux→ScalarSourceSink`), so "one object, several arrows" is still true |
| "…`CollisionOperator`…" (the name) | **REFUTED at write** | Class did not exist on 2026-08-08 — collapsed into `MultiplicationOperator` at `bbe8a51d` 2026-06-26 (#261) |
| "dispatch on the argument at apply time" | **CONFIRMED** | S :1093 (4 arms), F :469 (3 arms), C :376 (2 arms); unchanged since the plan (`[M]` 0 commits on `transport/operators/` since 2026-08-08) |
| "They are kernels, not operators" | **CONTESTED BY LANDED DOCTRINE** | `integral_kernel_operator.py` (landed `bbd9c5b6`/`bbe8a51d` 2026-06-26): "A Kernel REFINES LinearOperator … NOT disjoint"; `@runtime_checkable` Protocol requiring the LinearOperator surface PLUS `.kernel`; category gate `tests/transport/test_integral_kernel_category.py` |
| "Polymorphism belongs at construction (`.on(V)`), not at apply" | **PARTIALLY-LANDED, differently** | Measure-binding IS at construction (frame `cached_property`); carrier polymorphism stays at apply BY DESIGN (cross-method arms). No `.on(V)` exists: `[M]` `grep -rn "\.on(" orpheus/` → only `GeneratingMeasure.on(a,b)` (interval restriction, unrelated) |
| Two-cause table (S/F = measure ambiguity → bind to frame; C = carrier shape → one carrier) | **CAUSE 1 RESOLVED pre-plan; CAUSE 2 SURVIVES on all three** | S/F bound to the frame since `cc6d022e` 2026-06-24; carrier dispatch spans S (4), F (3), C (2) — the residual is cause-2 everywhere, and the tree's landed answer is typed multi-arm entry points, not "one carrier type" |
| "R·Λ·M is the definition of the binding" | **CONFIRMED + LANDED (pre-plan)** | `FrameBase.conjugate` (frame.py:205) returns `OperatorProduct(R, OperatorProduct(Λ, M))`; production since `cc6d022e`; Λ = `LegendreMomentScattering` (scattering.py:115) carries real spaces so the composability guard validates natively ("NO cast") |
| "F is the rank-1 case \|χ/4π⟩⟨νΣf\|" | **CONFIRMED** | fission.py:258 `kernel` = `RankOneOperator(χ, νΣf, axis=0) & IdentityOperator()`; 0-ULP crosscheck `tests/sn/operators/test_fission_kernel_crosscheck.py` |
| v2 amendment: "frame/measure/mass matrix all live on the **space's axes**; `.on(V)` is only as well-defined as V's axes are" | **REFUTED by a PRE-plan ruling** | The tree ratified the opposite ownership 2026-06-25 (`d30d4a68`, "W-E resolution — scattering owns its eigenbasis frame (Funk–Hecke)"): frame.py:205ff "The frame is then OWNED by that operator …, not by the phase space — see `:ref:frame-eigenbasis-ownership`" (docs/theory/foundations/frame.rst; also operator_algebra.rst:2877). The angular binding works today with NO Phase S |
| bind(K)† = bind(K†) ⟺ R = M* ⟺ tight — gate every binding on it | **NOT LANDED (as a gate)** | `[M]` `grep -rn "tightness\|frame_bounds\|B_over_A\|is_tight" orpheus/` → 1 hit, unrelated prose (rules_product.py:170). The DISCRIMINATION is carried structurally instead: `GalerkinFrame` (test is trial) advertises `M* = R` up to the dual factor; `PetrovGalerkinFrame` explicitly does not (frame.py:334-373 docstrings); `FrameBase.gram` REFUSES `GramStructure.DENSE` trials with `NotInvertible` (frame.py:297) |
| v2 addendum: the tightness FAMILY (4 joints) documented as one family | **NOT LANDED** | No such documentation exists (same grep; no "tightness" in docs/theory beyond quadrature-degree prose) |
| LossRepresentation split rule: answer→operator, cost→strategy; "closure → binding of L; traversal → strategy; ERR-026 closes structurally" | **SPLIT SUBSTANTIVELY REALIZED — but around a LIVE LossRepresentation, and ERR-026 was ALREADY closed** | Answer-side spatial closure = `SNMesh.scheme` (mesh construction, augmented_mesh.py:208/265); cost-side traversal = the strategy Protocol `LossRepresentation` (loss_representation/`__init__.py:251`) selected by `default_for(mesh)` (streaming.py:404); strategies: `CumprodScan`:1359 / `MovingFrontierWindow`:1540 / `FullFieldWavefront`:1822 / `ScanMarch`:2203. ERR-026: **CLOSED 2026-06-12** (error_catalog.rst:1806, ERR-058 closure-seed fix, #98→#99→#168→#195) — pre-plan |
| "the closure is also what gives the discrete space its trace content" (v2, §I.3 last ¶) | **VINDICATED, now MEASURED (post-plan)** | #344: ker A is a TRACE object (bulk share of the null projector `1.1e-28`); DD's `ψout = −ψin` face involution generates it; LD (in-cell slope) has `dim ker A == 0` on the identical box (loss_kernel_gauge.py module doc §1) |

### Corrections #13 and #20

| Correction | Verdict | Evidence |
|---|---|---|
| **#13** — `IntegralKernelOperator` is "a contradiction in terms; split into `IntegralKernel` (data) and the bound arrow produced by `.on(V)`" | **NOT LANDED; tree carries the explicit counter-design, pre-plan** | The Protocol (2026-06-26) inverts the plan's arrow: operator is primary, `.kernel` EXPOSES the integral structure (itself a LinearOperator — R∘Λ∘M *includes* M, so it is NOT measure-free). The measure-free residue that DOES exist: `LegendreMomentScattering` (scattering.py:115) / `N2NMomentOperator` (:300) — moment-space cores holding only `mat_xs` data + SH spaces. `[M]` `grep -rn "class \w*Kernel" orpheus/` → 5 classes, **0 of 5** are data-only integral-kernel objects (Protocol, a Basis, 3 unrelated Protocols) |
| **#20** — Phase 1 gated on Phase S ("the live space layer cannot answer; audit F1–F3") | **GATING MOOTED for the angular binding; Phase S itself not landed** | The operator-owned frame bypasses the V-axes dependency entirely (ruling above, 2026-06-25 — 6 weeks pre-plan). No S1/S2-shaped axis refactor exists (`[M]` `orpheus/numerics/space.py` has 4 lexical "axes" hits; spaces are concrete classes under `orpheus/numerics/spaces/`). Metrics DO already ride spaces in two places the F1–F3 audit framing under-credits: `SphericalHarmonicSpace.inner_product_weights` = g_C (ERR-039 Phase 1, 2026-05-26) and `G_bulk = V_cell·w_n` read via `SNMesh.full_field_space` (augmented_mesh.py:845). → **for sibling A** to adjudicate F1–F3 in detail |

### III-A rows (plan :288-296)

| Row | Verdict | Evidence |
|---|---|---|
| Frame hierarchy `FrameBase → PetrovGalerkinFrame → GalerkinFrame` (#268); scattering documents S = R∘Λ∘M via Funk–Hecke | **CONFIRMED, and stronger than claimed** | frame.py:114/:334/:374; M/R realized as `_FrameAnalysis`:412 / `_FrameReconstruction`:446, exposed at `.analysis`:195 / `.reconstruction`:200; Funk–Hecke at frame.py:223 + scattering.py:26. Not just "documented": `frame.conjugate` IS the production path (cc6d022e) |
| "missing: the kernel/operator split around it and the tightness gate" | **CONFIRMED as written, misleading in frame** | The plan's split (data + bound arrow) is absent — but a kernel/operator DISTINCTION (Protocol + `.kernel` exposure) already existed at write time, resolved in the opposite direction |
| `KEigenvalue` fuses Problem and Solver; split into `CriticalityProblem` + `PowerIteration` | **CONFIRMED at HEAD (not done)** | `KEigenvalue(Generic[V])` at `orpheus/numerics/iteration.py:1192`; `[M]` `grep -rn "class CriticalityProblem" orpheus/` → 0 |
| `power_iteration`'s `EigenvalueSolver` Protocol speaks fission vocabulary | **CONFIRMED at HEAD** | `orpheus/numerics/eigenvalue.py:106` Protocol; `compute_fission_source` at :132, consumed at :420 |
| "`LossRepresentation`: deleted; contents split by the answer/cost rule" | **REFUTED (as deletion)** | Alive as the strategy Protocol, loss_representation/`__init__.py:251` ("One algorithm for the within-group transport solve and its twin", `sweep`/`sweep_transpose`/`loss_action`/`loss_action_transpose`). The answer/cost split happened AROUND it (see §I.3 row above): it became the cost-side owner, not a deleted object |
| "ERR-039: … land the in-flight fix *with* the tightness gate" | **REFUTED (stale ≥2.4 months at write)** | Catalogue (error_catalog.rst:3107-3260): CAUGHT + fixed at introduction 2026-05-10 (Wave 0); re-fixed properly 2026-05-26 (Phase 1, `refactor/moment-space-and-layering`, commits 0eb9cf3..c5be4b0 — branch merged). No ERR-039 branch exists at HEAD (`[M]` `git branch -a`) |
| "`domain=None` escape hatch: verified live; closes as the quotient instance, forcing consumer = homogeneous solver" | **"Live" CONFIRMED; closure NOT LANDED** | The Optional-space pattern is the SHIPPING default on S/F/C (`None` ⟹ guard skips); no `EnergyGroupSpace` (`[M]` grep → 0); the arrow-narrowing campaign tracks the closure differently: per-law `xfail(strict=True)` markers ARE the todo list (test_b3_domain_narrowing.py:21, the seven-law domain gate) |
| Guide drift: `BlockOperator` (v3) vs `CoupledOperator` (code) | **CONFIRMED** | `CoupledOperator` at `orpheus/numerics/coupled_system.py:617`; no `BlockOperator` class (`[M]` grep) |

### Phase 1 table (plan :404-414)

| # | Item | Verdict | Evidence |
|---|---|---|---|
| 1.1 | `.on(V)` binding for S and F; gate: "apply-time singledispatch deleted from S, F" | **NOT LANDED** — gate red on both: singledispatch alive (S :1093, F :469); no `.on(V)`. The R/M-supplying machinery exists but is frame-owned-by-operator, not V-supplied (see #20) |
| 1.2 | C stays a multiplier; gate: "no frame import in multiplication_operator.py" | **GATE SATISFIED (was already true at write)** — `[M]` `grep -n frame multiplication_operator.py` → 1 hit = a memory-file *path in a docstring citation* (:92), zero imports. The work content was #261 (2026-06-26, pre-plan) |
| 1.3 | L's binding carries the closure; ERR-026 closes structurally; gate: apply/solve closure-identity test | **PARTIALLY-LANDED pre-plan, EXTENDED post-plan, premise stale** — ERR-026 closed 2026-06-12; the L21 one-operator design (sweep/matvec/adjoint through ONE strategy, single-sourced σ) shipped under #257/#280; no test literally named "closure-identity" (`[M]` grep → 0) but `tests/sn/sweep/core/test_sweep_vs_apply_consistency.py` is the ERR-026 tripwire and **12** test files carry `catches("ERR-026")` (`[M]` `grep -rl | wc -l`). Post-plan the closure seam got the interrogation surface the plan never designed (see §3) |
| 1.4 | Tightness gate on every binding; "land with ERR-039" | **NOT LANDED; premise refuted** — no tightness machinery (grep above); ERR-039 closed May 2026 with a STRUCTURAL class-closure instead: four operators (S₀, Πᵀ=w·S₀, Π*=g_C·S₀ via the generic `_AdjointOperator` + space metric, R=(2ℓ+1)·S₀) given "separately-typed homes; conflating any pair is structurally prevented" (catalogue :3204-3218). Live catchers: **8** `catches("ERR-039")` in `tests/numerics/test_spherical_harmonic_space.py` (`[M]` grep -c), incl. `test_H_equals_g_C_times_S0`:226, `test_R_equals_2l_plus_1_times_S0`:155. NOTE: the classes the catalogue's fix narrative names (`MomentProjection`, `HarmonicMomentReconstruction`) were themselves retired by the Frame campaign — `orpheus/numerics/projection.py` now holds ONLY the `AnalysisOperator`/`ReconstructionOperator` ABCs (:88/:136); the four-operator separation lives on as frame faces + space metrics. 1.4 is a standalone item if still wanted |
| 1.5 | Split `IntegralKernelOperator`; gate: "no kernel claims operatorhood" | **NOT LANDED — and the tree's landed doctrine claims the exact opposite** ("the Kernel is a *refinement* of" LinearOperator, gated by `test_integral_kernel_category.py`). This is a design CONFLICT to adjudicate, not a pending task |
| 1.6 | ~~EnergyGroupSpace~~ → S2 quotient instance | **NOT LANDED** (S2 absent; no EnergyGroupSpace — consistent with the supersession, nothing to do until S2) |
| 1.7 | Cache kernels, not bound operators; gate: "binding twice from one kernel shares Λ storage" | **NOT LANDED as architected; intent satisfied differently; tree carries a counter-ruling** — S/F bound operators ARE cached (on `SNSolver.__init__`, solver.py:1333/:1339); L+C **deliberately NOT cached**, with the in-code reason (solver.py:1323-1332: a solver-held copy "can silently drift from the one the sweep uses (it did: the former `self.L`/`self.S`/`self.F` triple was production-dead and misnamed)"); kernels rebuilt per access (S.kernel "Built fresh per access"; F.kernel plain property); frame cached per operator; `full_field_space` cached on SNMesh (:1041); `face_transmission_spectrum` cached per (closure, ndim) (scheme.py:372-375). The GATE's intent holds structurally: every kernel/operator is a read-through view onto ONE `MaterialXSField`, so Λ storage is shared by construction, not by a cache layer |

## 3. What the plan does not know (post-2026-08-08 movement)

`[M]` `git log --oneline --since=2026-08-08 -- orpheus/sn/operators/ orpheus/sn/loss_representation/ orpheus/transport/`:

- **#344 loss-kernel-gauge campaign** (`f934ff57` 2026-08-14, `b51bc802` 2026-08-15, `1a2be025` 2026-08-15, plus `5def63b0`):
  - `orpheus/sn/operators/loss_kernel_gauge.py` — `A = L+C−S−B` is **exactly singular** on all-reflective DD Cartesian boxes; kernel built in **closed form** (no eigensolve); `LossKernelBasis(Basis)`:679, `LossKernelGauge(LinearOperator)`:1019 = the G-orthogonal projector, an endomorphism of the **boundary-trace space** (the kernel is a trace object — bulk share `1.1e-28`).
  - **The gauge consumes the plan's own binding machinery**: each block projector is `frame.conjugate(InverseMetricOperator(frame.gram))` — "ONE spelling of 'project onto a span' in the package and this is a consumer of it" (loss_kernel_gauge.py:1035-1043). The frame-conjugation constructor became load-bearing in exactly §I.3's direction, frame-owned.
  - Wiring: `SNMesh.loss_kernel_gauge` property (augmented_mesh.py:974); every solver exit routes `_exit_gauge_trace` (solver.py:589); `Solution.IterationHistory.gauge_correction` (solution.py:232); exhaustiveness gate `tests/sn/solve/test_every_entry_gauges_its_trace.py`.
  - `DiscretizationSchemeBase.face_transmission_spectrum(ndim)` (`orpheus/transport/spatial/scheme.py:911`) — "Does this closure leave any face mode UNDAMPED? — **ASKED, never declared**" (drives the scheme's own `cell_kernel_batch` on a probe cell; ask-don't-tabulate). This is the closure-interrogation seam the plan's 1.3 gestured at ("the closure decision is simultaneously the trace-content decision") realized as a queryable spectral contract on the scheme type.
- **`c33178ef`** 2026-08-12 — "gate the angular-closure seam, and stop crediting Hebert for a tau he never wrote": the ANGULAR (M-M τ) closure seam got its own gate; touches sn/operators + loss_representation.
- **`1689faf4`** 2026-08-08 — THE FLIP: `SNMesh(CYLINDRICAL)` admits exactly the carrying rules (R12a admission at the mesh).
- **`f8ecb4f6`** 2026-08-08 — ERR-078 fix (ψ½ march solve dropped the outflow-row rhs).
- **`a6fd7a08`** 2026-08-10 — "the per-group volume integral belongs to the MESH" (→ sibling A: metric/measure ownership movement).
- **`adc887d6`** 2026-08-09 / **`144cdf51`** 2026-08-17 — falsified-claims fixes + 85 declared `implements` edges (docs).
- **Meanwhile `orpheus/transport/operators/` and `orpheus/numerics/{frame,projection}.py` have 0 commits since 2026-08-08** — the S/F/C/Protocol surfaces the plan targets are exactly as they were at plan-write; the plan's staleness is inherited from v1, not caused by post-plan drift.

## 4. Issue-overlap map

| Issue | State | Maps to |
|---|---|---|
| **#359** — M-M angular recurrence has THREE production spellings + two adjoints (declared, bit-identity-gated twin) | OPEN | Plan 1.3 territory: the ANGULAR closure's algebra is not yet one bound object. A "binding carries the closure" refactor is the natural dissolution site; note the closure seam is now GATED (`c33178ef`) but not unified |
| **#316** — frame-backed `angular_moment(ℓ)` reduction (DSA d₀/d₁ as the frame's ℓ-th analysis-face row) | OPEN | Plan 1.1's consumer-side completion: M (the analysis face) as the single reduction verb across in-sweep accumulation / `integrate_angular` / DSA restriction. Confirms the frame IS the binding home — and that binding-vocabulary adoption is incomplete downstream |
| **#279** — migrate CP/MoC isotropic in-scatter onto shared `K_iso = IsotropicScattering + IsotropicN2N` (`orpheus/transport/operators/isotropic_scattering.py:229/:345`; SN routed 0-ULP at P2) | OPEN | Plan 1.7's intent in cross-method form: ONE kernel (the Σ_s0ᵀ group-transfer) → many per-method arrows differing only in the angular wrapper. The kernel/binding separation already operating as a migration program |

## 5. Refuted premises (one-line structural reasons)

1. **"S/F/C have no unique (domain, codomain) pair"** — declared Optional pair landed ≤ `bbe8a51d` 2026-06-26, six weeks pre-plan; the v1 audit was retained into v2 without the §7 tree reconciliation.
2. **"CollisionOperator"** (named as live) — collapsed into `MultiplicationOperator` 2026-06-26 (#261); the plan names a class that did not exist when it was written.
3. **"ERR-039 … land the in-flight fix"** — nothing was in flight on 2026-08-08: fixed 2026-05-10, re-fixed 2026-05-26, branch merged; the class was closed by typed separation + space metrics, not by a tightness gate.
4. **"ERR-026 closes structurally"** (as a Phase 1 outcome) — already CLOSED 2026-06-12 by the ERR-058 closure-seed fix (catalogue :1806); the surviving r→0 issue is ERR-059, a different root (WONTFIX-for-DD).
5. **"LossRepresentation: deleted"** — alive at loss_representation/`__init__.py:251` as the cost-side strategy Protocol; the answer/cost split was realized by MOVING the closure to the mesh and the selection to `default_for`, not by deleting the object.
6. **v2 amendment "frame/measure live on the space's axes ⟹ Phase S is a prerequisite of the binding"** — the tree ratified operator-owned eigenbasis frames (`d30d4a68` 2026-06-25, `:ref:frame-eigenbasis-ownership`) and the binding has been production since 2026-06-24 with no Phase S; the prerequisite claim re-derived a question the tree had answered.
7. **Correction #13's prescription** (kernel = data, operator = `.on(V)` arrow) — the tree's landed, gated doctrine is the inverse (Kernel REFINES operator; `.kernel` exposes, nothing binds); adjudicate the conflict before any Phase 1.5 work, because landing #13 as written means REVERSING a documented design, not filling a gap.

## Structural summary (durable; line numbers drift, re-derive via Nexus)

At HEAD the tree has **already built the §I.3 binding machinery in its own
dialect**: Λ (moment-space cores holding only XS data) + M/R (frame faces
carrying the measure) + a construction-time binding verb (`frame.conjugate`),
with the frame owned by the operator whose eigenbasis it is (Funk–Hecke), and
the kernel EXPOSED off the operator rather than the operator produced from the
kernel. What the plan's Phase 1 still correctly identifies as absent: the
apply-time carrier dispatch on S/F/C (deliberate, cross-method), a numerical
tightness/B-over-A gate (discrimination is type-structural instead), and any
`.on(V)` spelling. What it gets wrong is the baseline: its diagnosis describes
the pre-#261 tree, and two of its prescriptions (#13, the v2 space-ownership
amendment) now contradict ratified, documented, gated design decisions rather
than fill gaps.
