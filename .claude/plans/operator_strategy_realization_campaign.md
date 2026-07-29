# Operator · Splitting · Realization — the three-layer separation campaign

> ## ✦ COLD-PICKUP MARKER — read this block first, then §2, then §9
>
> **Written 2026-07-28 at a deliberate compaction point. Re-anchor from THIS FILE +
> `git log`, never from a conversation summary.**
>
> ### Where the tree is
> Tier 1 (the naming-honesty carve) **merged to `main` @ `b0a003b4`**; the three
> ruling commits (`9cc69420`, `6946fe30`, `360a8087`) sit on top.
>
> **P0 is IN PROGRESS on branch `refactor/operator-strategy-layers`.**
> **O-1 is DISCHARGED @ `9a546640`** — `tests/sn/architecture/` exists and carries
> AC-b (R7 pinned, `xfail(strict=True)`), AC-b′, and M-5/M-6/M-7 with every number
> MEASURED on that commit. No production change has been made. `WithinGroupSystem`
> and `_select_si_splitting` are now SAFE TO TOUCH — they were not, before that commit.
>
> Any "merged / not merged / in-flight" line in this file is a point-in-time snapshot
> that lies forward. **Verify with `git status -sb` / `git merge-base --is-ancestor`;
> never trust the prose** (`.claude/rules/process-discipline.md`).
>
> ### Read order
> 1. this block · 2. §2 (the acceptance criterion — it was CORRECTED, do not use the
> first version) · 3. `.claude/plans/campaign_verification_plan.md` §0 and §6 (the
> measured red set, and the 12 ordering constraints) · 4. §4 (issue reconciliation —
> this campaign owns no new roadmap) · 5. §5 (six rulings not to regress).
>
> ### Four corrections already absorbed — do NOT re-derive them
> 1. **The obvious acceptance criterion is signature-tautological.** "The posed equation
>    is invariant across `inner_schedule`" measures exactly 0.0 and *cannot* be made red:
>    `build_within_group_system` takes no strategy parameter. A hand-injected falsifier
>    moves it (5.15e-02), so it *looks* non-vacuous. Use the three-legged criterion in §2.
> 2. **The real red is R7 — the splitting is a TWIN.** `gauss_seidel` runs
>    `ScheduledInvertibleOperator` + `SNMaskedBoundaryOperator` while the record still
>    advertises `StreamingCollisionOperator` + `SNBoundaryOperator`. `jacobi` is the
>    control leg.
> 3. **A 1-D slab HIDES R7** — the G-S arm fires only on a seedless multi-D Cartesian
>    mesh, so a slab reports `True` on both rows. The first probe of this used a slab and
>    got the wrong answer. **Use 2-D Cartesian.**
> 4. **`N_AA = S + B_a` is NOT a violation** of what O.2a retired. O.2a retired the fixed
>    *named triple* (`SNSolver._scattering_with_boundary_op`); a sum inside a grid slot is
>    a lowering detail, not a role assignment. See §1.1 and §6/P3.
>
> ### The immediate next action
> ~~P0 item 1 (O-1): write AC-b RED.~~ **DONE @ `9a546640`.** Next: P0 items 3 and 4
> — the leaf gates G1.1/G1.3/G1.4/G1.5 (O-2/O-3, before P1 makes the degradations
> unrepresentable) and the P-4 performance baseline (O-6, a hard blocker for stage 3).
> See §11 for live status.
>
> ### At risk if the working tree is disturbed
> `.claude/skills/vv-principles/{SKILL.md,error_catalog.md}` carry **two new failure-mode
> sharpenings** (Mode-8 *signature-tautological gate*; Mode-9 *the degeneracy kills the
> RATE gate too*, `ρ ≡ 0` on the invariant subspace). Those files are **forbidden to
> commit** by standing policy, so they exist ONLY in the working tree with no git history
> to recover from. A `git checkout` / `restore` / `stash` on those paths destroys them
> irreversibly (lesson L28 — this has already cost the project two catalog entries once).

**Status:** PLAN OF RECORD, awaiting sign-off. Authored 2026-07-28 from a design
dialogue (see §10 for what it supersedes).
**Branch (when it starts):** `refactor/operator-strategy-layers`
**Binding lens:** [[feedback-build-the-machinery-operator-algebra]] — the goal is the
MACHINERY, and the operator-algebra formalism is realized ALWAYS (the four questions:
*what is it · how is it posed · what invariant tests it · how is it solved*). Risk and
effort set the SEQUENCE, never the WHETHER.

---

## 1. Thesis

The code conflates three genuinely independent degrees of freedom:

| | what it answers | today |
|---|---|---|
| **the operator** | what the physics *is* | entangled with the splitting |
| **the splitting** | which part is implicit | committed at operator-construction time |
| **the realization** | how it is computed | selected by `str` flags on one entry point |

Conflating them makes the design space **multiplicative** — visible as
`inner_solver` × `inner_schedule` × `acceleration` × geometry × scheme ×
carrying/seedless, with the illegal combinations hand-managed by refusals scattered
across the solver (`_within_group_si` carries two; the DSA arm two more).

Separating them makes it **additive**: a strategy becomes a *value*, and "is this
combination legal?" stops being a maintained refusal and becomes a schedulability
question the machinery answers.

### The pipeline

```
1. PHYSICS REALIZATION   leaves L, C, S, F, B (+T=1/v) — each MONOMORPHIC
                         (exactly one domain, one codomain), method-owned,
                         built by a realizer
            ↓
2. PROBLEM POSING        role assignment into  A_loss ψ = λ M ψ
                         + the μ → physical-eigenvalue map.  METHOD-AGNOSTIC.
                         k:  A = L+C−S−B,      M = F,    k = μ
                         α:  A = L+C−S−F−B,    M = 1/v,  α = −1/μ   ← F MIGRATES
                         adjoint = a daggered row, not a layer
            ↓
3. STRATEGIC PARTITIONING   a partition of the state space induces a block grid;
                         a POSING assigns each block — or a PIECE of one — to
                         implicit or explicit:   A_loss = M_impl − N_expl
            ↓
4. LOWERING + SOLVE      realize the schedule, form M_impl⁻¹, run a thin algorithm
```

**The compact statement, to be preserved verbatim through every rewrite of this plan:**

> a **partition** induces the block grid · a **posing** assigns each block (or part of
> one) to implicit or explicit · the **schedule** follows from the implicit set's
> dependency structure

### 1.0 The PAIRED-CONSTRUCTION ruling (user, 2026-07-29)

> **The loss and gain operators are built as a PAIR, so nothing can get lost.** A split of
> one operator into two must happen **at a single place**, with both halves **absorbed
> into their targets immediately and unmissably** — never "construct `B_lower` here, fold
> it into the loss; then N lines below construct (or forget) `B_upper` and hand-roll it
> into the gain, hoping it was done right."

**Formalization — make the ASSIGNMENT the primitive and DERIVE `M`, `N`:**

* `A` decomposes into **pieces** (blocks, or sub-blocks after refinement);
* every piece carries exactly **one label** — implicit or explicit;
* `M = Σ(implicit pieces)`, `N = Σ(explicit pieces)`.

Then `A = M − N` holds **by construction of the labelling, not by assertion** — a piece
cannot be lost, because every piece has a label and the labels partition the piece set.
A *split* is a **refinement**: replace `P` by `(P₁, P₂)` with `P₁ + P₂ = P`, gated once at
the refinement, then label each independently. The only mutation the pair admits is a
**transfer**, invariant-preserving by arithmetic:

```
move P from explicit to implicit:   M ← M − P ,   N ← N − P
                                    M′ − N′ = M − P − N + P = A     ✓
```

**Why this is the strongest available control: it makes the σ_r bug UNSPELLABLE.** The bug
wanted `M = LC − Σ_s0·𝕀`, but `Σ_s0·𝕀` **is not a piece of `A`** (the pieces are `LC`,
`S = Σ_s0·P_iso + Σ_{s,ℓ} + …`, `B`). Refining `S` and labelling `Σ_s0·P_iso` implicit
gives the honest `M = LC − Σ_s0·P_iso`, which the capability check then refuses because it
is not sweep-invertible. The bug required **substituting a different operator during the
move**; `transfer` gives that substitution nowhere to happen.

**Current state — the precedent is half-right.** `SNBoundaryOperator.split`
(`sn/operators/boundary.py:431`) IS atomic and returns a *named* pair, and its docstring
already states the defence: *"Returns a named pair so the two construction sites cannot be
swapped silently."* What is missing is the other half: the absorption is two loose
expressions (`return LC - parts.lower, (S, parts.upper)`, `solver.py:774`) — write `(S,)`
and the complement vanishes with nothing to catch it — and it happens in
`_select_si_splitting`, a **different place** from where the record's `(M, N)` was built.
That second fact **is R7**.

**The productive caveat — this discipline is NARROWER than general splittings, and the
narrowing is a feature.** Classical Jacobi (`M = D`, `N = D − A`) invents an `M` that is
not a sum of `A`'s pieces. So there are two categories, and the discipline separates them:

| | how `M` is formed | `A = M − N` | spectral gate |
|---|---|---|---|
| **partition splitting** | assignment over `A`'s own pieces | **by construction** | rate only |
| **approximation / preconditioner** | invent `M`, define `N = M − A` | by definition | **REQUIRED** — nothing structural bounds ρ |

Sweep, angular Jacobi, boundary G-S, block-triangular — all category 1, safe by
construction. **The σ_r fold is category 2**, which is exactly why it needs the §3 spectral
gate and why its home is #200 rather than the splitting machinery. The discipline
classifies it automatically instead of relying on memory.

Two consequences worth stating because they are easy to lose:

* **Posing precedes partitioning.** The posing decides which leaves are even *in*
  `A_loss` — for α-eigenvalue `F` migrates into the loss. Partition first and you
  partition the wrong operator.
* **A leaf can be split ACROSS the implicit/explicit boundary.** `B` does not go
  implicit *or* explicit — it splits into `B_lower` (folded into the sweep) +
  `B_upper` (lagged). The implicit operator is a subset of *pieces*, not of leaves.
  This is why `SNBoundaryOperator.split` exists (`sn/operators/boundary.py:431-448`).

---

### 1.1 The precedent — this codebase already made this move once, and it worked

The campaign is not a new idea imported from outside. **Wave O step O.2a made exactly
this move at one layer**, and recorded the reasoning
(`docs/theory/foundations/boundary_conditions.rst:754-772`):

> **Why variadic — the fixed triple encoded a false posing distinction.** :math:`S`,
> :math:`F` and :math:`B` are *homogeneous* in the driver… The fixed :math:`(A, S, F)`
> triple gave :math:`S` and :math:`F` named slots the *resolvent layer* never uses — it
> was **encoding a posing-layer role assignment at the iteration layer, where it does not
> belong**. Collapsing the triple to a homogeneous `*gains` bag **moves the role
> distinction back to the posing layer (its proper home)** and lets a fourth gain slot in
> as a *data addition* rather than a new named slot.

That is this campaign's thesis, verbatim, applied to one interface. The campaign
generalizes it: O.2a moved *role assignment* out of the iteration layer; P1–P5 move
*carrier choice* out of application (into construction), *splitting choice* out of
equation construction (into strategy), and *strategy choice* out of `str` flags (into
data). Same argument, three more places.

The dual is also on record and is the failure this campaign prevents: the σ_r fold put a
**numerical** decision (which term may be lagged) inside an **operator** construction,
where nothing could gate it — and it shipped 46–56 % silent errors.

---

## 2. Acceptance criterion — corrected

> ### ⚠ The first draft of this section was WRONG, and the correction is instructive.
>
> It proposed *"changing a solve strategy must not touch operator construction"*, asserted
> as RED. **It is GREEN — at exactly 0.0 — and it cannot be made red**, because
> `build_within_group_system(sn_mesh, mat_xs, *, scattering_op, scattering_order)` takes
> **no strategy parameter at all** (verified by `inspect.signature`). The invariance is a
> SIGNATURE fact, not a behavioural one: the knob cannot physically reach the object, so
> the gate is green in every possible run. A hand-injected falsifier *does* move it
> (5.15e-02 when `B` is dropped), which is exactly what makes this trap dangerous —
> the probe looks non-vacuous.
>
> This is now catalogued as a `vv-principles` **Mode-8 sub-case: the SIGNATURE-tautological
> gate**. A whole campaign can be anchored on such a criterion and be unfalsifiable from
> the first commit.

### The real red — R7, the splitting is a TWIN

**MEASURED independently by the main agent** (2-D Cartesian, seedless, heterogeneous 2G,
`level_symmetric` S4, reflective+vacuum):

| `inner_schedule` | driver runs | record advertises | `is` |
|---|---|---|---|
| `jacobi` | `StreamingCollisionOperator` + `SNBoundaryOperator` | *the same objects* | **True** |
| `gauss_seidel` | `ScheduledInvertibleOperator` + `SNMaskedBoundaryOperator` | `StreamingCollisionOperator` + `SNBoundaryOperator` | **False** |

The record (`WithinGroupSystem`) advertises one splitting; the G-S driver **silently
discards it and re-derives another** via `_select_si_splitting` (`solver.py:862`). The
`jacobi` row is the built-in **control leg** proving the asymmetry is real.

*(Note for whoever re-runs this: a 1-D slab shows `True` on BOTH rows — the G-S arm only
fires on a seedless multi-D Cartesian mesh. The main agent's first probe used a slab and
got the wrong answer. Use 2-D.)*

**This is the campaign's motivating defect, stated exactly:** stage 2 (posing) and stage 3
(splitting) are welded into one record, so there is no named boundary at which "strategy
may enter" can be asserted — and because there is no boundary, a second splitting grew
beside the first and the first went silently stale.

### The criterion, three-legged

* **AC-a — the strategy-entry boundary is a SIGNATURE fact.** No function on the
  `mesh + materials → posed pencil` chain accepts any strategy token
  (`inner_solver`, `inner_schedule`, `tol`, `restart`, `corrector`, …).
  **RED today because the chain does not exist** — the function that produces the posed
  `A` also produces the splitting, so the gate is *unwritable* against `main`. That is
  the strongest possible statement that the boundary is missing.
* **AC-b — one splitting per posed equation.** The driver runs the *same objects* the
  record advertises (R7). **RED today on `gauss_seidel`, green on `jacobi`.**
* **AC-c — the posed object is byte-identical across strategy**, fingerprinted on the
  OBJECT (sha256 of its image of a fixed probe basis), never on a spectrum or a `k` —
  Mode 12. **Unwritable today**: no pencil exists, and `A − M` currently *raises*
  (`IncompatibleOperatorComposition` — `A` and its own `M` live on different carriers).

Full specification, mutation register and tolerances: §9.

---

## 3. Dependency readiness — verified against the tree, not remembered

Per [[feedback-plan-states-dependency-readiness]]. All at `main` @ `b0a003b4`.

| capability | state | evidence |
|---|---|---|
| `LinearOperator[Domain, Codomain]` two-param Protocol | **LANDED** | `numerics/operator.py:504`; TypeVars `:117-118` with PEP-696 `Codomain` default |
| `Vector` Protocol (#256 step 1) | **LANDED** | `numerics/vector.py:94` |
| Source-role carriers (`ScalarSourceSink`, `HarmonicMomentSourceSink`, `AngularSourceSink`) | **LANDED** | `transport/source_sinks/` — a 7-module package |
| Typed factored edges `M`, `Λ`, `R` | **LANDED + documented** | `operator_algebra.rst:2941-2953`; realized as `frame.reconstruct(Λ.apply(φ))/W` |
| `CoupledOperator` block grid — matvec, block solve, triangularity, `inverse()`, `assemble()` at offsets | **LANDED** | `numerics/coupled_system.py:617`; `_triangular_orientation` `:850`; `CoupledSubstitutionOperator` `:1212` |
| Declared structural zeros | **LANDED** | `None` block IS the zero map — *"coupling sparsity is structural, no zero arithmetic runs"* (`:627-628`) |
| `BlockRole` = a partition-slot tag on bulk⊕trace | **LANDED** | `operator.py:208-231` — `BULK` (C,S,F) / `BOUNDARY` (B) / `FULL` (L, the only irreducible one) |
| Realizer seam precedent | **LANDED** | `BoundaryTraceLaw → SNBoundaryRealizer → fixed-signature operator` |
| `ρ` as a typed method | **LANDED** | `Displacement.contraction_ratio(previous)`, `transport/displacements/_displacement.py:92`; `true_error_estimate(ρ)` `:120` |
| Four-tier eigenvalue architecture | **DOCUMENTED + realized** | `operator_algebra.rst:3938+`; commits `650032e`/`7603c8e` |
| SCC criterion for a dependency digraph | **LANDED (narrow)** | `derivations/discrete/sn/sweep_acyclicity.py` @ `d9b092c2` — trace digraph only; generalizes |
| DSA consistent low-order system | **LANDED** | `sn/acceleration/dsa.py`; #2 closed 2026-07-27 |

### Open gaps that BITE this campaign

| gap | evidence | which phase it blocks |
|---|---|---|
| `FissionOperator.domain -> "FunctionSpace \| None"` | verified by introspection | P1 — the root |
| `.H` silently Euclidean when spaces are `None` | `numerics/operator.py:1221-1226` applies the metric only when both are non-`None`; `homogeneous/solver.py:193` constructs space-anonymously | P1 |
| `explicit_gains` — presence-dependent arity + positional `B_a`-LAST | `coupled_system.py:320-322`, parsed at `solver.py:855` | P3 |
| `A == M − N` **ungated** | `WithinGroupSystem` has no `__post_init__` check; correctness rests on two constructor arms written correctly side by side | P0 |
| `1/k` applied OUTSIDE `F` | `compute_fission_source`, `solver.py:1126-1132` — *"the 1/k division stays at this level"* | P4 — the pencil has nowhere to put λ |
| No pencil type | 0 hits for `GeneralizedEigenPencil`; `KEigenvalue(A,S,F)` holds the triple but its ONLY production consumer is the adjoint path (`solver.py:2472`) | P4 |
| One instance = 7 arrows (F) / 6 (S); three dispatch mechanisms | explorer map, `scratch/review_map_reaction_operators.md` | P2 |

---

## 4. Issue reconciliation — this plan owns NO new roadmap

The predecessor plan cited **zero** GitHub issues while duplicating ~10. Cardinal Rule 4
makes issues the plan-of-record; this campaign is an *execution vehicle* for issues that
already exist, and it must not mint a parallel backlog.

| issue | relationship | phase |
|---|---|---|
| **#296** block-operator abstraction | **This campaign IS its forward quest.** #296 already maps space→block-triangular, angle→block-Jacobi, energy→block-G-S. Do not re-derive. | P5 |
| **#261** relocate C/F/S cores | **PREREQUISITE.** Also carries the *lone L3 blocker* (`MaterialXSField`'s `SNMesh` type vs `tests/test_layer_imports.py`) the predecessor plan omitted. | P2a |
| **#256** field-typed algebra | Its **open fork 3** ("multi-arm union vs typed siblings") is answered by this campaign: typed siblings, at construction. | P2c |
| **#205** storage × role typing | Becomes the domain/codomain of the monomorphic operators rather than a separate typing project. | P1/P2 |
| **#213** capabilities as a morphism class | The capability-tag → declared-arrow move is the same idea; `StreamingCollisionOperator` (renamed @ `8367346f`) is the first instance. | P2b |
| **#288** cross-class source-sink dunder | Pressure relieved (roles become operator endpoints), not necessarily closed. | P1 |
| **#289** `FullField` generic over leaves | The `π_bulk` / `j_bulk` named injections are its restriction/injection pair. | P5 |
| **#200** block-inverse preconditioner | **Consumer.** Also the honest home of the σ_r fold (§7). | P6 |
| **#277** RQI / bordered matrix | **Consumer** — needs the pencil; already keeps `(A,F)` un-fused on ndarrays. | P4 |
| **#273** energy-group Gauss–Seidel | **Consumer** — becomes a partition choice, and is the natural *second axis* demonstration. | P5 |
| **#260** `SumOfTensorProductsOperator` | Λ = ⊕ℓ Σ_s,ℓ is the reduced operator's natural form. | P2a |
| **#297** `Ravellable` Protocol | Touched by the block flat-dimension work. | P3 |
| **#318** NDA | **Consumer** — a nonlinear posing that re-poses between iterations. | future |
| **#320** wire the sweep-cycle criterion | The schedule layer is its consumer; `sweep_acyclicity.py` generalizes. | P5 |
| **#219** MethodSpace + **registry** | ⚠ **PARTLY STALE** — its F-C proposes a `MethodSpaceBuilder` *registry*, but the realizer registry was **DISSOLVED** at #290 P7b (`44d583e`) in favour of `TransportMethod`. This campaign uses **explicit realizers, no registry**. Update #219 F-C rather than execute it. | P2b |
| #263, #287, #302, #306 | adjacent; untouched | — |

**New issues this campaign may file:** exactly one — the σ_r fold's corrected splitting
(§7), unless it folds into #200.

---

## 5. Rulings this campaign MUST NOT regress

Each was made deliberately, with rationale and gates. Violating one is a session failure.

1. **Pattern M stays until #261 lands.** `operator_algebra.rst:513-524` parks the
   apply-spelling question on #261 *by ruling*: "deciding the spelling before the cores
   move would be premature; the sharing should dictate the form." P2c is that decision,
   arriving in the sanctioned order. Pinned by the 630-line
   `tests/sn/operators/test_operators_apply_typed.py`.
2. **The fused ndarray scattering kernel survives** as the 0-ULP oracle
   (`test_scattering_kernel_crosscheck.py`). Under monomorphism it is not a second
   *arrow* — it is a second *realization* of one arrow, living inside the operator's
   apply, with the typed edges as its oracle.
3. **`HarmonicFrame` stays in `transport/`.** All its consumers are inside `transport/`
   (`scattering.py` ×7, `frames/__init__` ×3, `reaction_rate_functional` ×2); moving it
   to `sn` strands them and violates `tests/test_layer_imports.py`.
4. **No global string-keyed realizer registry** (#290 P7b dissolved it). Realizers are
   constructed and passed explicitly.
5. **`Kernel` keeps its established meaning** — `IntegralKernelOperator`
   (`transport/operators/integral_kernel_operator.py:164`): *a Kernel IS a
   LinearOperator, and adds the `kernel` member*; grand-report §5.6 suffix law: "`Kernel`
   means it is integrated against a measure." Pinned by a 13-test intrinsic-property gate
   with negative cases. **This campaign mints no non-callable "kernel" layer** — the
   reduced operators already exist and are already monomorphic.
6. **The diffusion resolvent stays a direct LU.** `NEVER` route it through the
   structure-keyed `A.inverse()` — the Green splitting diverges for fine-mesh elliptic
   operators (`diffusion/solver.py:48-50`).

---

## 6. Phases

Each phase is independently mergeable and leaves the full suite green. Compatibility
shims live one merge cycle only (`coding-standards.md`).

### ⬢ P0 — Instrument and gate the status quo *(no production change)*

1. The **acceptance-criterion gate** (§2), red today.
2. The **splitting invariant** `A == M − N` on `WithinGroupSystem`, with mutation teeth:
   re-introduce `N = 0` (the exact σ_r bug) and confirm RED. *This alone retro-catches a
   46–56 % silent-error defect.*
3. The **performance harness** (§8 R1) — wall-clock + allocation baselines, captured
   BEFORE any composition lands.
4. Baselines: representative 1-D/2-D SN eigenvalue, fixed-source, anisotropic,
   radial-characteristic, diffusion, homogeneous, adjoint.

**Exit:** every later phase has a gate that can red.

### ⬢ P1 — Monomorphic domains *(the root, first half)*

1. `domain`/`codomain` become **non-Optional** on the reaction operators.
2. Give the meshless/homogeneous path a real space (today it passes
   `basis_shape=(ng,1)` explicitly *because* the operators carry no domain —
   `homogeneous/solver.py:194,202`).
3. The `.H`-without-metric trap becomes **unrepresentable**, not guarded.

**Exit:** no production operator can be constructed without stating its arrow.

### ⬢ P2 — The realizer seam *(the root, second half)*

* **P2a — relocate the C/F/S cores** (#261 steps 1–2). Respect its named L3 blocker.
  Λ's natural form is #260's sum-of-tensor-products.
* **P2b — `SNReactionRealizer`**, on the `SNBoundaryRealizer` precedent. Each method
  builds its own monomorphic realized operator as a composition
  `j_bulk ∘ E₀ ∘ F_G ∘ M₀ ∘ π_bulk`. **Composition is the algebra of record and the
  equivalence oracle; a fused path may follow, pinned bit-identical.**
* **P2c — collapse the dispatch.** Now — and only now — the spelling question is
  answerable (ruling 1). The `@overload` surface, the three dispatch mechanisms, and
  #256 fork 3 all resolve together.

**Exit:** every realized operator has one public `apply` signature and a stated arrow.

> ⏸ **COMPACTION POINT.** Commit, then re-anchor from this file + `git log`.

### ⬢ P3 — Shape symmetry: loss, M and N take the SAME partition

> **User ruling (2026-07-29) — the SHAPE-SYMMETRY principle.** The gain operator composes
> like the loss operator; implicit and explicit are the same *kind* of thing
> operator-wise, so one partition machinery builds both and lowering serves both. The
> machinery that builds System A must build System B and the coupled A⊕B.
>
> **This is required, not aesthetic.** Implicit-vs-explicit is a property of the
> **assignment**, not of the operator — G-S *moves a piece of `B` across the boundary*.
> A block cannot move from `M` to `N` if `M` and `N` have different shapes. Same-shape is
> the precondition that makes the piece-split expressible; its absence is exactly why the
> current G-S must un-name and re-derive, which **is R7**.
>
> **MEASURED (2026-07-29) — the carrying arm already complies; the seedless arm is the
> anomaly, breaking it three ways:**
>
> | | `loss` | `implicit_operator` | `explicit_gains` |
> |---|---|---|---|
> | carrying (sphere) | `CoupledOperator(2×2)` | `CoupledOperator(2×2)` | `tuple[CoupledOperator(2×2)]` |
> | seedless (slab) | `CoupledOperator(1×1)` | `StreamingCollisionOperator` | `tuple[Scattering, SNBoundary]` |
>
> **This RETIRES the "gain-grid partition" open question** posed in an earlier draft
> ("does the seedless gain grid partition by system or by bulk⊕trace?"). The question was
> mis-framed: it treated the gain grid as independently choosable. It is not — the SYSTEM
> carries a partition and `loss` / `M` / `N` all inherit it. Choose once; three follow.
>
> **The symmetry extends to `F`.** `A`, `M`, `N` and `F` are all maps from the state space
> to the source/residual space — which is precisely what makes the pencil contract
> `domain(A) == domain(M_eig)` expressible (§P4). The pencil is this same symmetry one
> level up, and it is unspellable today for the same reason: `Optional` domains and a
> tuple for `N`.
>
> **Where the symmetry stops** (both real, neither a shape asymmetry):
> *capability* — `M` must be invertible, `N` need not be (that is #213's morphism class,
> `Iso ⊂ PartialIso ⊂ General`); and *lowering* — `M`'s inverse realization is
> schedule-dependent while `N` only ever needs `apply`. Stage 4 is allowed to treat them
> differently; that is what a lowering stage is for.
>
> **Composability is already the direction of record.** #296: *"`CoupledField` = N systems,
> each a complete `(interior ⊕ boundary)` composite — today's `FullField` generalizes to a
> carrier-generic `System[Interior, Boundary]`."* Unrealized only because bulk⊕trace and
> systems are two unrelated types today (`Composite`/`FullField` vs `CoupledField`).
>
> **P3 therefore becomes**: make the seedless arm shape-consistent using the partition the
> carrying arm already uses (concrete, small, independently valuable). Making the partition
> a first-class *choosable* object stays P5. The acceptance target is the table above with
> all three columns identical on both rows.

1. `explicit_gains: tuple[...]` → a `CoupledOperator` gain grid on **both** arms. The
   carrying arm already holds exactly this (`gains=(N,)` with
   `N = CoupledOperator([[N_AA, None], [emission, B_b]])`); only the seedless arm is a
   bare pair.
2. The positional `B_a`-LAST convention and the arity-by-mesh parse dissolve.

   ⚠ **This does not re-open the `S+B` fold O.2a retired — but not for the obvious
   reason.** The carrying arm *already* sums them: `N_AA = S + B_a`
   (`coupled_system.py:521`). What O.2a retired was the **fixed named triple** — the
   deleted `SNSolver._scattering_with_boundary_op`, which packed the reflection into the
   middle slot of a `(A, S, F)` signature. A summation inside a *grid slot* is a
   **lowering** detail; a named slot in the driver signature is a **role assignment**.
   The grid also preserves the property O.2a was protecting — extensibility by *data
   addition* (a future B-trace or α-time term is a new row/column, not a new named
   parameter). ~~Open design question: does the seedless gain grid partition by system or
   by bulk⊕trace?~~ **RETIRED as mis-framed** — see the shape-symmetry ruling above: the
   gain grid is not independently choosable; the SYSTEM carries the partition and
   `loss`/`M`/`N` inherit it.
3. `M⁻¹N` becomes formable ⇒ **ρ is probeable as an operator expression.**

**Exit:** the RHS is an operator, not a positionally-parsed tuple.

### ⬢ P4 — Posing + the pencil

1. `GeneralizedEigenPencil(lhs, rhs, spectral_parameter)` with the contract
   `domain(A) == domain(M)`, `codomain(A) == codomain(M)` — **now non-vacuous**, because
   P1 made domains non-Optional.
2. Move `1/k` out of `compute_fission_source` and onto the pencil's spectral parameter.
3. Posing rows: k / α / fixed-source / adjoint. **Constraint (L23a): the pencil is an
   ADDITIONAL Layer-2a object; it must NOT displace `power_iteration`'s late binding**,
   which is strictly more general (it admits CP BiCGSTAB, diffusion's direct inverse, the
   homogeneous dense solve — none of which has an `(L,S,F)` factorization).

**Exit:** the eigenproblem is an object; α and shift-invert are posing rows, not engines.

> ⏸ **COMPACTION POINT.**

### ⬢ P5 — Partition + schedule

1. `Partition` — restriction/injection pairs, with the **reconstruction identity**
   `Σᵢ JᵢRᵢ = I` as its correctness gate.
2. A **second constructor** on the existing block interface: derive `A_ij = Rᵢ A Jⱼ` from
   `(A, partition)`. *Not* a parallel `PartitionedOperatorView` type — one interface,
   two constructors. (This resolves the fork between the predecessor plan and #296.)
3. **Schedule** = SCC of the implicit set; generalize `sweep_acyclicity.py` off the trace
   digraph (#320).
4. Demonstrate **two** axes: the carrier partition (already built) + energy G-S (#273).

**Exit:** a strategy is a value; `_select_si_splitting` and the `str` flags retire.

> **User ruling (2026-07-28): `ScheduledInvertibleOperator` is expected to DISAPPEAR with
> this phase — do not rename it beforehand.** It exists only as what
> `_select_si_splitting` returns when the G-S schedule folds `B_lower` into the implicit
> operator. Once "fold a piece of `B` into the implicit set" is a *partition + assignment*
> rather than a bespoke type, the class has no reason to exist. Its name states a
> capability rather than the object (the defect fixed for `InvertibleOperator` →
> `StreamingCollisionOperator` at `8367346f`), but renaming a symbol slated for retirement
> is churn. **Retire, don't rename.** If P5 slips or the class survives the carve, the
> rename returns to the table.

### ⬢ P6 — `solve auto` + the σ_r fold's honest home

1. A policy object reading `Displacement.contraction_ratio` and branching:
   ρ small → SI; ρ → 1 → Krylov (splitting-invariant) or DSA.
2. The σ_r fold re-enters as a **correct** splitting (§7) with its measured ρ, wired to
   #200 as a preconditioner candidate — never as a stationary iteration.

**Exit:** the system can choose its own strategy, and every splitting carries a spectral
claim.

---

## 7. The σ_r fold — corrected, and its honest role

Measured this session (`/tmp/fold_algebra.py`, to be promoted into an algebra of record):

* The historical bug was `N = 0` — claiming `A = M`. The **consistent** splitting is
  `M = (L+C) − Σ_s0·𝕀`, `N = M − A = −Σ_s0·(𝕀 − P_iso)` — the *anisotropic remainder*.
* As a **stationary iteration it is provably divergent**: `ρ = c/(1−c)` — 0.25 at c=0.2,
  **1.00 at c=0.5**, 9.0 at c=0.9. On a real 40-cell S₈ slab: 6.91 at c=0.9, 17.99 at
  c=0.99. *This is the divergence in the project's memory; it was never a bug.* The
  mechanism: it annihilates the slow isotropic mode (eigenvalue exactly 0) and pays by
  amplifying the anisotropic modes.
* As a **Krylov preconditioner** it is defensible: measured κ(M⁻¹A) 18.99 vs the plain
  sweep's 26.84 at c=0.99; a wash at moderate c. Caveat: the `σ_r > 0` guard degrades it
  exactly where it is most wanted (c→1).

⇒ Keep `foldable_part`/`residual_part`, **fixed**: give them the honest splitting API
(return `(M, N)`, gated by `A == M − N`) rather than retiring or merely renaming them.
The concept is correct; only `N = 0` was wrong.

---

## 8. Risks and controls

**R1 — Performance regression. THE dominant risk, with precedent.**
Lesson L16: a refactor once moved hoisted vectorised work into a per-cell Python fold and
cost **10–20× on slab, ~6× on the suite**, with *every correctness gate green*.
`A_ij = Rᵢ A Jⱼ` as a composition risks exactly that shape.
**Control:** composition is the algebra-of-record and equivalence oracle; the fused path
is the production realization, pinned bit-identical. The P0 performance harness is a
merge gate for P2b and P5 — **neither may merge on correctness alone.**

**R2 — Abstraction that does not pay.** Pattern 6 requires ≥2 instances. Satisfied (the
ψ½ 2×2 is instance one; the hard-wired axes are the rest) — but a `Partition` type with a
single derived-grid consumer would be speculative.
**Control:** root-first; consumers pull. P5 ships **two** axes or it does not ship.

**R3 — The spectral constraint has no structural proxy.** `ρ(M⁻¹N) < 1` is not derivable
from sparsity; the σ_r fold is the worked counterexample.
**Control:** every posing carries a spectral claim with a gate (§9), and Mode-9 discipline
— **never** gate a splitting on the fully-reflective isotropic box, where the wrong
formulation is accidentally exact.

**R4 — A long-lived branch.** **Control:** P0/P1 are independently valuable and
separately mergeable; the campaign can stop after any phase.

**R5 — Structural zeros still need declaring** for performance (not for correctness).
**Control:** a partition declares its sparsity, as `CoupledOperator` already does with
`None`.

---

## 9. Verification plan

**The verification plan is `.claude/plans/campaign_verification_plan.md`.** It is
NORMATIVE for this campaign — where it and this file disagree, it wins (it is the one
grounded in measurement). Headlines:

* **It corrected §2** (above) — the acceptance criterion was signature-tautological.
* **New test package** `tests/sn/architecture/`, `pytestmark = [foundation]` (software
  invariants, no theory `:label:` — the verifies⊥level doctrine).
* **Per-stage intrinsic-property gates** G1.1–G1.6 (leaves), G2.1–G2.5 (posing),
  G3.1–G3.5 (partition), G4.1–G4.4 (lowering) — each with its defining law.
* **A 22-entry mutation register** (M-1…M-22), every one naming the exact mutation and
  the expected RED signature, several with **measured** magnitudes (drop `B` → 5.15e-02;
  sign-flip `N` on the sphere arm → 8.57e-03; the `_GaussSeidelResolvent` falsifier →
  2.667). M-3 is a deliberate **negative** control (changing `restart` must NOT red AC-c).
* **The spectral gate** §3, with the Mode-9 trap stated in its exact algebraic form: an
  isotropic seed puts the σ_r fold's iteration operator in an invariant subspace where
  `N ≡ 0`, so `ρ` reads **identically 0** — the splitting looks *optimal*. The harness
  must project out the isotropic component, and ships the isotropic seed as a permanent
  **control leg** (S4).
* **The performance gate** §5 — three legs, and the *catcher* is the deterministic
  **call-count scaling** leg (leaf-kernel entries must not scale with `n_cells`), not the
  timing leg. Wall-clock is normalised against an in-process calibration rather than an
  absolute ms threshold.
* **12 phase-ordering constraints** (O-1…O-12) from the verification side. Two are hard
  blockers: **O-1** (write AC-b RED *before* touching `WithinGroupSystem` — fix it first
  and the fix is unprovable) and **O-6** (the perf baseline must be captured at the end of
  P2; stage 3 may not merge without it).
* **§4.2 ratifies four PRINCIPLED re-baselines** with their three `vv-principles` criteria
  discharged — two of which are *already* non-zero on `main` (the coupled arm's
  `A == M − N` measures 3.5e-17, not 0). A uniform bit-identity contract is refused.

Two additions to the plan's own scope that came from it:

* **G2.3-α is mandatory and currently has zero coverage.** The α posing is the *reason*
  `F` must migrate, and there is no α reference in the tree. Build the independent
  closed form first (`scipy.linalg.eig` on the hand-posed `G×G` pencil, ~15 lines) — O-4.
  Otherwise the first α number the code produces becomes its own baseline.
* **File ERR-070** for the σ_r fold when S5/S6 first catch a regression — #215 closed
  without an ERR entry, and the catalog ends at ERR-069.

Standing requirements regardless:

* Every math-bearing type ships a test of its **defining laws**
  ([[feedback-test-intrinsic-properties]]): partition ⇒ reconstruction identity;
  splitting ⇒ `A = M − N`; pencil ⇒ the domain/codomain contract; schedule ⇒ acyclic ⇒
  one-pass exact.
* Every mutation gate names its exact mutation and expected RED signature.
* Where a composition replaces a fused reduction, **principled re-baseline** per the
  `vv-principles` three criteria — never a silently loosened tolerance.

---

## 10. What this supersedes

* `.claude/plans/sn_operator_realization_and_posing_plan.md` — **SUPERSEDED.** Its
  diagnosis was largely right; its grounding was against the grand report rather than the
  tree, so ~⅔ re-derived existing architecture, three items would have regressed
  deliberate rulings, and it cited no issues. Retain as archaeology.
* `.claude/plans/sn_operator_realization_plan_REVIEW.md` — the grounding pass that
  produced this plan's §3/§4/§5. Input, now folded in.
* `.claude/plans/operator_inverse_algebra_carve.md` Phases 2–5 (`.inverse()` returns an
  operator; retire `CAP_SOLVE`) — **absorbed** into P2b/P4.
* `.claude/plans/coupled_block_operator_campaign.md` — **DONE**; it is instance #1 of the
  block machinery P5 generalizes. Archaeology, not a twin.

---

## 11. P0 — live status

**Branch `refactor/operator-strategy-layers`. No production change in P0**: everything
here is a gate that reds against `main`, or a number recorded from it.

1. ✅ **AC-b + AC-b′ — DONE @ `9a546640`** (O-1 discharged).
   `tests/sn/architecture/` minted: `_config.py` (the verification plan's §8
   mandatory-configuration table, ONE home — O-11's rule applied to the fixtures) and
   `test_stage_separation.py`. **14 passed, 1 xfailed; pyright CLI 0 errors.**
   * **AC-b** pins R7 as `xfail(strict=True)`, and `--runxfail` confirms it reds for the
     *right* reason ("driver ran `ScheduledInvertibleOperator`, record advertises
     `StreamingCollisionOperator`"), not incidentally — an xfail otherwise hides any
     failure (Mode-11 discipline).
   * Control legs shipped alongside: `jacobi`; both sphere rows (schedule is inert on
     the carrying arm); the **Krylov** path green on both arms (R7 is SI-specific — if
     that row reds, the defect spread); and **the slab trap as a named row**, so the
     wrong answer this campaign already got once cannot be got again.
   * **AC-b′** per-arm and MEASURED: seedless exactly 0 ULP, carrying 4.48e-17 (passes
     at `nulp=8`, fails at 4). A uniform bit-identity contract is refused.
   * **M-5** 1.00e-01 / 1.83e-02 · **M-7** 3.32e-02 / 3.66e-02 · **M-6** the exact σ_r
     fold: honest `N` 2.93e-17 green, `N=0` **5.43e-03 RED** on an anisotropic state and
     **3.57e-18 INVISIBLE** on an angularly-flat one. Leg 3 is Mode 9 in closed form and
     ships as a permanent control; the *mechanism* (`Σ_s0·P_iso ≡ Σ_s0·𝕀` on a flat flux,
     measured 2.045 both) is asserted directly so the control stays falsifiable.
   * ⇒ **`WithinGroupSystem` and `_select_si_splitting` are now safe to touch.**
2. ⏳ **G1.1 / G1.3 / G1.4 / G1.5** — the leaf gates, red against today's
   `FissionOperator` (`domain is None`), `SNBoundaryOperator` (raw `AttributeError`) and
   `F.H` (silent Euclidean). O-2/O-3: they must exist **before** P1 makes any of those
   unrepresentable, else the degradation loses its catcher. Same `xfail(strict=True)`
   convention as AC-b.
3. ⏳ **P-4** — capture the performance baselines (O-6; a baseline measured after a
   regression is worthless). P-1 (call-count scaling) is the real catcher and is
   contention-immune; the P-2 wall-clock constant must be captured on a **quiescent**
   tree or it is not trustworthy.

**The convention this phase establishes:** every currently-red campaign gate ships as
`xfail(strict=True)` naming its red-set ID, the phase that flips it, and the instruction
to delete the marker on XPASS. **The strict-xfail set IS the campaign's todo list**, and
strictness means no phase can silently fix something without acknowledging it.
