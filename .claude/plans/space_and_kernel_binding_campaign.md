# Operators born bound — the space substrate, the cone, and kernel→operator realization

**Campaign 1 of the two-campaign return to the operator-machinery thread.** Chartered
2026-08-19 from the user's rulings (§0) on the grounded
`orpheus-operator-machinery-report-v2.md` (its Part VI carries the tree reconciliation at
HEAD `7aae9bf1`; memos `scratch/omr_v2_grounding/{A,B,C,D}_*.md` carry the `[M]`
evidence). Successor-in-MECHANISM to the operator-strategy campaign's P1/P2 — that plan's
gates and rulings SURVIVE and bind here (see the 2026-08-19 banner in
`operator_strategy_realization_campaign.md`). **Campaign 2** (late lowering: the
LossRepresentation overturn → `GeneralizedEigenPencil`, resolvent machinery with
complexification, partitioning/schedule) is chartered only after this campaign lands.

**The outcome.** (1) The spaces a transport problem lives on EXIST as first-class
objects — Energy first, then Spatial and Angular, with SN's space their composition —
carrying structural identity, space-owned metrics, and per-axis partition hooks.
(2) Flux lives in the **positive cone K of an ordered vector space V**; the affine-torsor
field algebra is overturned. (3) S and F are **kernels** (representation-free data) that a
constructor binds to a frame and a space, returning fully-bound monomorphic OPERATORS;
apply-time dispatch retires.

---

## 0. Rulings of record (user, 2026-08-19) — none may be regressed

1. **Flux lives in the positive cone K ⊂ V.** "What matters is correctness. We should not
   be bound by past mistakes." The torsor field algebra is overturned and the cone
   implemented. Adjudication + consequences: §4 (phase CS3) and
   `orpheus-operator-machinery-report-v2.md` Part VI.5 D1.
2. **Kernel → Operator at construction.** "What today is named `ScatteringOperator` will
   become `ScatteringKernel` and inside a constructor will bind the Kernel to the Frame to
   return the actual `ScatteringOperator`, fully bound." Operators stop being overloaded
   on application; they are specialized at construction. (Correction #13's direction; the
   "Kernel REFINES LinearOperator" Protocol doctrine is overturned — Part VI.5 D3.)
3. **Space first.** "Before we sort this, we will work on space definition, because the
   Homogeneous solver relies on the existence of an Energy space … Then we will add
   Spatial and Angular spaces, and SN will be a composition of Energy, spatial and angular
   spaces (and this will power the partitioning machinery later on)."
   ⚠ Execution-order note (user, 2026-08-19, post-compaction): "first" binds
   spaces-before-the-BINDING (CS4), which the re-ordered sequence preserves. The cone
   carve (CS3) was verified independent of the space work and now executes before the
   space phases — see the §4 sequencing block for the derivation.
4. **Campaign split.** This campaign = Space + Kernel→Operator realization (+ the cone
   carve, which must precede the dispatch collapse — §4 rationale). Campaign 2 = the
   LossRepresentation overturn — "currently an early decision on partitioning" — to
   `GeneralizedEigenPencil`, resolvent machinery (with complexification), and the
   partitioning machinery: "a late decision on how the system is finally lowered."
5. **Inherited and still binding** (operator-strategy campaign): PAIRED CONSTRUCTION,
   SHAPE SYMMETRY, L23a (the pencil is ADDITIVE; `power_iteration`'s late binding stays),
   "`ScheduledInvertibleOperator` disappears with the partition phase — retire, don't
   rename", O-4/O-5 precede any α work (⛔⛔ never resolve α by grepping `alpha`).

---

## 1. Dependency readiness — what the tree provides (verified 2026-08-19 vs `7aae9bf1`; spot-re-certified at `00085baf` post-CS3-R — the P1 ledger row reproduces `[M]` 105 passed / 21 xfailed, and the two CS3-consumer rows below carry their landed state)

| capability | state | evidence |
|---|---|---|
| `FunctionSpace` is carrier-generic (PEP-696 `Carrier` TypeVar) | shipped — the surface CS2 builds beneath | memo A §3.F5 (`space.py:92-100`) |
| Frame machinery `FrameBase → PetrovGalerkinFrame → GalerkinFrame`; `frame.conjugate(Λ)` IS the production R∘Λ∘M | shipped — CS4's binding verb exists | memo C §1 (`frame.py:114/:205`, `scattering.py:488/:599`) |
| Moment-space cores holding only XS data (`LegendreMomentScattering`, `N2NMomentOperator`; F's rank-1 dyad) | shipped — the Λ half of the kernel split already separated | memo C §2 (`scattering.py:115/:300`, `fission.py:258`) |
| Q6 measure metadata (`ReferenceMeasure`, `ExactnessClaim`, `invariance_group`, `half_range_clean`) | shipped — the GENERATOR side of the Angular axis; CS2 FORGETS from it | memo A §4 S2 row |
| `DiscreteMeasure.partition_by` | shipped (2026-05-10) — the measure-level partition primitive the per-axis hooks forget from | memo A §4 (`measure.py:1042`) |
| The P1 ledger: 12 `xfail(strict=True)` markers + `tests/sn/architecture` | shipped; `[M]` 105 passed / 21 xfailed reproduced 2026-08-19 | memo B §8.8 (node ids) |
| #340 record machinery (`IterationRecord`, `StoppingCriterion`, named criterion trajectories, `IterationBudget`) | shipped — ✅ CS3 relocated the iterate diagnostics there (`f9d571b5`: `increment_norms` + derived `contraction_ratios` / `true_error_estimate()`) | memo B H3/§6.1 |
| Cone fragments: `is_positivity_preserving` (DD `False` + numerical witness), ψ≥0 realizer refusal, coefficient cone-as-predicate battery, ray normalization in `power_iteration` | shipped — ✅ CS3 added the element predicate `Field.cone_violations` + the DD witness on top (`f9d571b5`); #390 tracks the flag's first production reader | memo D §3 |
| Gaps that bite: F1 identity aliasing (`[M]` live probes: `Γ₊(GL8) == Γ₊(product(4,2))` True with unequal metrics); F4 densification (`[M]` state-size weight tensor); NO GroupAxis anywhere; 4–5 metric doctrines over 6+ sites (F3); #369 (measure-identity twin); #295/#297 (layout/carrier fragmentation) | live — this campaign's subject | memo A §§2–3 |

⚠ Worktree note: if this campaign runs in a worktree, rebuild Sphinx inside it and
`use_workspace` the Nexus graph (L22). The line numbers above drift — re-derive via Nexus
at pickup.

---

## 2. Phase CS1 — the homogeneous solver poses its problem on a real Energy space

**Goal.** Energy stops being an anonymous integer in a shape tuple: an Energy space
exists as a first-class object, and the homogeneous (infinite-medium) solver is its
forcing consumer — its operators know they live on Energy alone.

**Proposed means** (hypothesis, 2026-08-19): a `GroupAxis`/Energy space minted by
forgetting from the group structure (§I.8's forgetful decomposition), realized as the
maximally-quotiented member of the symmetry-reduction family (§I.11: trivial spatial ⊗
ℓ=0 ⊗ Energy — the "per unit volume" intensive convention IS the normalized quotient
measure). The homogeneous solver constructs it; F/C there bind `domain = EnergySpace`.
Design detail (nodal axis ⟹ coordinate cone; the V/V* collapse hook of tightness-family
(iv) declared even if unexercised) fixed at design time.

**Done when** (checkable): the four model-generic strict-xfails
`test_model_generic_leaf_declares_a_space[C-2g, C-4g, F-2g, F-4g]` are DELETED on XPASS
(`tests/sn/architecture/test_monomorphic_leaves.py:670` at `00085baf`);
`homogeneous/solver.py`'s space-less operator constructions (the
`_assemble_loss_operator(mat_xs)` / `FissionOperator.from_solver_data(mat_xs=…)` pair at
`:193-194`, which leave `.domain` None at runtime — ⚠ there is no literal
`domain is None` string to grep; the predicate is "the constructions pass a space") and
BOTH `basis_shape=(ng, 1)` spellings (`:194` and `:202` — the double-pass) are gone.

**§8 blast-radius note (enabler ≠ neutral).** Giving F/C a real domain flips the
composability guards from *skipped* to *active* on every composition that today rides
the `None`-skip — BOTH of them (re-anchored by SYMBOL at `00085baf`, CS3-R prep; line
numbers drift): `OperatorSum`'s `_agreed_space((a, b), "domain"/"codomain", …)`
(skips per-operand on None) and `OperatorProduct`'s `A.domain == B.codomain` check
(`IncompatibleOperatorComposition`, ≈`operator.py:1629` — skipped when either is None;
the homogeneous path composes BOTH: `C − K_iso` is a sum, `K = M⁻¹ @ F` a product).
And `.H` flips from Euclidean-fallback to metric-applied — `_AdjointOperator`'s
`inner_codomain.apply_metric(y) if inner_codomain is not None else y`
(≈`operator.py:1297-1307`). Enumerate the newly-checked compositions and the `.H`
consumers on the homogeneous path BEFORE landing; gate at the tier the change is
observable.

## 3. Phase CS2 — SN's space is a composition of Energy, Spatial, and Angular

**Goal.** SN's bulk space is literally the composition Energy ⊗ Spatial ⊗ Angular of
first-class per-axis objects — one composition mechanism, structural identity, metrics
owned by the spaces, per-axis partition hooks recorded for Campaign 2.

**Proposed means** (hypothesis): `SpatialAxis = forget(Mesh)` (owns V_cell — volumes move
OUT of operators; sequence carefully, ERR-067-class adjoint checks per move),
`AngularAxis = forget(Quadrature)` (weights + parity closure; directions stay generator
data; the Q6 registry tags are the raw data), Energy from CS1. The identity charter
(report item 0.8): structural `__eq__`/`__hash__` by interning at the generators,
logged-then-flip (S3's intermediate mode); metric layer per-axis (deletes
`_broadcast_metric` and F4's dense outer product — a memory assertion that no state-sized
weight tensor is allocated).

**Done when**: memo A §2's live aliasing probes REFUSE (`Γ₊(GL8) == Γ₊(product(4,2))` →
False; FullField cross-decomposition → False; `SpatialMomentSpace(4,1) == (2,2)` →
False); the three product mechanisms are one (`[M]` census); the live mints route through
it — **the live mints are `augmented_mesh.py:1099` + `_bases.py:381/:503`; the space.py
factories are DEAD (test-only)** — and the LD moment-mass inline product
(`augmented_mesh.py:1095-1098`) becomes a per-axis metric composition; suite green
bit-identical, or principled-equivalence per `vv-principles` where a reduction order
legitimately moves (document the three criteria per relaxation).

**Scope decisions to fix at design time** (open, not ruled): #295 (layout Protocol) in or
explicitly out; #297 (carrier contracts) in or out; `WithTrace` tag now vs with the DFEM
stress case (report S4 — the weakest-confidence design point; the
`_RadialCharacteristicSubSpace` split pair is the concrete stress case in the tree);
whether the Γ ladder re-bases onto the axis machinery now or later (it is the
best-designed part of the layer — memo A praise table — churn only with cause).

## 4. Phase CS3 — flux lives in the cone: the field-algebra overturn

**Goal.** One field type per physical role, living in V; cone membership an element
PREDICATE; cone-preservation a REALIZATION flag; iterate semantics (ρ, ‖Δψ‖/(1−ρ)) on
the iteration layer. The affine gates and the displacement type family retire.

**Design rulings (2026-08-19, pre-carve — user on Q1/Q2, standing doctrine on Q3/Q4;
the measured forks are `scratch/cs3_verification_plan.md` §9.2):**
- **Q1 — ρ is defined on the SPACE norm** (user). Today's spelling already is
  `space.norm(values)`; Euclidean agreement is an accident of `inner_product_weights is
  None`. Consequence accepted: CS2's physical metric legitimately REDS the capture gate
  (`tests/numerics/test_si_diagnostic_trajectory.py`), which then re-derives with a
  regeneration note — principled over bit-identical.
- **Q2 — the diagnostic surface is `IterationRecord`** (user): named trajectories on the
  durable record via the #340 named-criterion channel, not solver attributes (#366's
  defect class). The two pinning tests re-point to thread the record out.
- **Q3 — the cone predicate returns the offending INDICES** (vv anti-#14: return the
  structure; emptiness = membership; `bool` derivable).
- **Q4 — `affine_combination` dissolves** (0 production callers; in V the blend is
  ordinary arithmetic; the relaxation concept's future home is the iteration layer).

**Sequencing — RE-ORDERED 2026-08-19 (user): CS3 executes FIRST. Execution order is
CS3 → CS1 → CS2 → CS4.** The user asked whether the cone carve is independent of the
space work; re-derived against the tree, it is. The original rationale (kept per
plan-authoring §3, refuted half marked):

- CS4's collapsed `apply` signature must be written against the settled field algebra,
  or the new bound operators bake in the flux-vs-displacement refusal (#331's fork) that
  Campaign 2's Krylov/resolvent work would immediately re-open — **STANDS** (CS3 before
  CS4; unaffected by the re-order).
- ~~Running CS3 after CS2 means the re-typed fields land on the final spaces (one pass
  over the field layer, not two).~~ ⛔ REFUTED 2026-08-19 on the independence
  re-derivation: CS3 changes the ROLE algebra (dunders, leaf set); CS2 changes the SPACE
  objects the fields point at — disjoint parts of the field layer. The one shared site is
  the fiber check, and `[M]` it is NOT the torsor machinery's: the mesh-binding guard
  lives in the retained `_check_partner` chain (`transport/fields/_bases.py:199` "Add the
  mesh-binding guard on top of Field's class/space gate"; `:643` "class identity, space
  equality, AND mesh identity"; `:837-838` names the space-aliasing hazard and closes it
  by mesh identity for face fields), and the base dunders route every same-class pair
  through it (`numerics/field.py:277-282` — `__add__`/`__sub__` call `_check_partner`
  first). Retiring the mixin drops `+`/`−` through to exactly that chain, so the fiber
  discipline ("different problems don't mix") survives the carve at today's exposure
  under either order. Cone-first in fact SHRINKS CS2's blast radius: the 7 displacement
  leaves retire before CS2 re-points space wiring, halving the leaf inventory whose
  identity semantics flip.
  ⚠ Carve-time check owed (cheap, one gate): a cross-mesh same-shape `ψ + ψ'` must still
  REFUSE after the flip — the negative control that the fall-through really lands on the
  mesh-binding chain. `_check_torsor_partner`'s own mesh arm (`_flux_role.py:196`)
  retires WITH the mixin: its job was the fiber check on a cross-class partner that
  Layer 1 cannot see, and it ends when the cross-class partner does.

**Proposed means** (the carve; blast radius = memo D §9, [M] 16 production files / 16
test files / 5 doc pages): retire `FluxRole`'s gates, the 7 displacement leaves, the
Rep-keyed registry and per-block displacement minting; KEEP the Layer-1 cross-class +
space/mesh identity gates (fiber discipline — "fluxes of different problems don't add"
was always about the fiber, never about affine structure); relocate
`contraction_ratio`/`true_error_estimate` onto the iteration layer
(`SourceIteration`/`IterationRecord` — the stop already rides flat norms there, memo D
§2); re-point the foundation battery (10 tests), the 7 TypeError-pin files and ~16
affine raise-sites; re-type DSA's displacement-minted returns; rewrite
`field_algebra.rst` → the cone chapter (report D7; the §"Why affine" argument at :242 is
the overturned one; 4 `:eq:`-cited labels are APIs — every citer re-read per
coding-standards); sweep the 5 torsor doc pages + `coupled_system.py:108,239,292` +
`coding-elegance` anti-pattern #18's worked example (the skill must not teach the retired
ontology); close #331 (operators are linear on V — the A/V ceremony dissolves).

**Done when**: `grep -r torsor` returns only the retirement note (report 0.7's real
gate); #331 closed with `S.apply(Δψ)` legal; the SI convergence diagnostics reproduce
from the iteration layer `[M]` on a c→1 fixture (old-vs-new ρ trajectory equivalence
gate BEFORE the type retires); a cone-membership predicate exists and DD's negative-flux
witness exercises it; suite green with re-pointed batteries and migrated
`catches`/`verifies` markers (retirement = marker migration).

**What survives untouched** (memo D's measured invariants): the stopping logic, the
source/sink/residual vector-space algebra (#288 orthogonal), the Krylov flat path,
`is_positivity_preserving` + witness, the realizer's ψ≥0 refusal, quadrature-weight
positivity, ray normalization.

**CS3 execution ledger (2026-08-19, branch `refactor/cone-field-algebra`):**

- ✅ Step 0 — pre-carve instruments `d7737f6d`: `scratch/cs3_verification_plan.md` +
  the ρ-trajectory capture gate (`tests/numerics/test_si_diagnostic_trajectory.py`,
  5/5 at rtol=1e-12; composite-norm control 4.71e-3; five-mutation battery).
- ✅ Step 1 — diagnostics → the record `c3e66b18`: `IterationRecord.increment_norms`
  + derived `contraction_ratios`/`true_error_estimate()` (⛔ NOT in `criteria` —
  `converged` is `all(cleared)`, a ρ entry would flip every verdict);
  `_principal_bulk_leaf` type-agnostic walk; `where_largest` promoted to `Field`;
  the three Displacement methods deleted same commit (twin prevention). `[M]` the
  frozen trajectory reproduces; sweep 3439 passed; walls green.
- ✅ Step 2 — the algebra flip `993fa280`: `FluxRole` DELETED (file + 7 leaf bases),
  `±` fall through to the fiber-guarded `Field` algebra; DSA re-typed in the same
  commit; 11 test files re-spelled per plan §3; the base-point surrogates in
  `test_declared_law_is_linear.py` became direct additivity + the difference law
  (sharper: the direct spelling reds on affine); 2-D matvec gate STRENGTHENED.
  `[M]` 4821 passed / walls green / mutations: M-add reds 12 value legs, M-mesh
  reds exactly the fiber row. Closes #331 at merge.
- ✅ Step 3 — package retirement `5efd2178`: `orpheus/transport/displacements/`
  deleted (zero consumers; marker-migration set measured EMPTY); ~~all production
  prose swept to dated past tense~~ ⛔ REFUTED by CS3-R (2026-08-19, same day):
  the sweep's own filters were partial — ~20 live-tense sites survived (the RC
  family's ⊖ prose incl. 2 error messages, dsa.py's noun, 5 gates' premise,
  2 `_bases.py` docstrings a case-blind grep missed, 3 AGENT.md briefs, 3
  agent memories) — repaired across CS3-R sweeps 1–5; the universal lacked its
  denominator (plan-authoring §2); the ρ≈c anchor renamed
  `tests/sn/solve/test_si_convergence_diagnostics.py`. `[M]` 9859 collected clean.
- ✅ Step 4 — the cone predicate `3b9e8651`: `Field.cone_violations` (offending
  indices, most-negative first; −0.0 member, nan violation; docstring carries the
  four does-NOT-claim negatives); unit legs + the one-parameter DD witness pair
  (`[M]` thick leg min ψ = −6.399383e-01, 2 of 8 negative, report ≡ the set).
- ✅ Step 5 — docs `8a6fc353` + follow-ups `e341e12a`/`c634919f` (archivist +
  main agent): `field_algebra.rst` → the cone chapter (602 → ~1540 lines; the
  six-argument adjudication; the affine era as dated history with claim-by-claim
  falsification); labels: `affine-torsor-algebra` RETIRED → `flux-vector-algebra`,
  diagnostics labels RENAMED `iterate-*`, residual label KEPT+annotated,
  `positive-cone-definition` MINTED, `affine-bc-form` untouched + disambiguated;
  32 dead Python-domain refs → 0; `operator_algebra`'s role-axis argument
  re-derived; coding-elegance #18 REVERSED; cross-domain-frames A.1/Shape-3 and
  numerical-bug-signatures Sig-9 stop teaching the torsor. `[M]` sphinx `-E -W`
  clean, warning set byte-identical to baseline; sentinels 0 violations (539→540);
  15 torsor survivors, each dated-historical.

---

### ⏸ COMPACTION POINT #1 — CS3 ⏹ COMPLETE (2026-08-19)

**Measured baseline at close**: full tree `-m "not slow"` serial = `[M]`
**1 failed / 9542 passed in 1:01:42** — the 1 red INHERITED from main
(`627e64d6` added `tests/sn/__init__.py`, refuting the xref checker's
namespace-package premise fixture; repaired at `88e46240`, re-pointed at
`tests/geometry`). **This branch's contribution: 0 new reds.** Walls: capture
gate 5/5 (rtol=1e-12), bit-identity 3/3 + DD 11/11 under `-W error::DriftWarning`.
`dead_references` → 0 after `c634919f`. Budget note: the gate ran 62 min
uncontended — the ≥90-min budget held with margin.

**Corrections that supersede older text in THIS file**: the §4 sequencing block
(CS3 first — the one-pass argument ⛔ refuted; the fiber guard was never the
torsor's); the four design rulings block (space-norm ρ / record surface /
indices predicate / `affine_combination` dissolves); the verification plan's
G1 ✅ dissolved (#389), G4 ⛔ false measurement, §4.5 leg-c ⛔ refuted
(boundary block nonzero on vacuum), §4.1 leg-d blind to the subtract mutation.

**Durable lessons** (candidates for rules/skills on distillation):
1. A trajectory diagnostic must NOT ride a criteria list whose conjunction IS
   the verdict (`converged = all(cleared)`) — the §8 enabler class, caught at
   design time by reading the consumer.
2. A direct algebraic law can be strictly SHARPER than its constrained-era
   surrogate (direct additivity reds on affine; base-point independence
   structurally cannot) — a carve that legalizes a spelling should re-derive
   which gates STRENGTHEN, not only which re-spell.
3. The mass-delete shadow hazard has a second face: an untracked `__pycache__`
   of a deleted PACKAGE materializes a phantom PEP-420 namespace package that
   still imports (`__file__ is None`, 0 members) — `git clean` the directory as
   part of the delete.
4. `[M]` `tests/sn/acceleration` has no `__init__.py` — `627e64d6`'s "audit
   finds no other broken chain" is REFUTED (coverage attribution for that
   subtree is broken the same way tests/sn's was). Belongs to the nexus
   test-DAG thread, not this campaign.

**Owed micro-edit**: `field_algebra.rst`'s development-history row and
`history.rst`'s ready-to-paste row carry `*(in development)*` — swap to the
merge hash on the next docs touch (CS1's D8 is the natural ride).

**State**: CS3 ⏹ COMPLETE; CS3-R (below) is the next step by user instruction;
CS1 (Energy space) and CS2/CS4 unstarted — after CS3-R, resume from §2 with the
§0 rulings and this file's §1 readiness table.

---

## 4-R. Phase CS3-R — the carve survives a clear-context adversarial review

**Chartered by the user at the CS3 close (2026-08-19, verbatim intent):** one
more pass through the machinery "just to check that all the torsor machinery was
replaced by cone one and find any opportunities to improve expressiveness (for
example, using cone properties when we can, expressing residuals, etc)" —
explicitly a review of the carve's own work **with clear eyes/context**, i.e.
run AFTER the compaction that follows this entry.

**Goal.** Two verdicts, independently reached: (a) **completeness** — no torsor
machinery survives un-replaced, where "machinery" means CONCEPTS and SHAPES, not
the word (vestigial names, torsor-shaped call patterns, base-point spellings,
docs/tests that re-teach the retired design by structure); (b) **expressiveness**
— the cone ontology's new affordances are USED where they pay (cone properties
at consumers, residual expression, gates that can now state laws directly).

**Method** (binding): the adversarial-first review discipline
([[feedback-adversarial-phase-before-balance]]) — Phase 1 harsh (how would I
break this / make it 100× better; reshaping is fair game; "it works" is not a
defence), Phase 2 re-evaluate and WITHDRAW what does not survive, with every
refuted candidate recorded with its structural reason (process-discipline).
Load `coding-elegance` + the nexus smell sweeps (`twin_paths`,
`discriminations`, `native_place`, `dead_functions`, `dead_references` — the
graph at HEAD post-merge is current). Read FIRST: this file §0/§4 + compaction
point #1; `docs/theory/foundations/field_algebra.rst` (the cone chapter);
`scratch/cs3_verification_plan.md`; memo `scratch/omr_v2_grounding/D_*.md`.

**Leads** (recorded at compaction prep by the carve's own author — UNVERIFIED
hypotheses, deliberately not pre-judged; the reviewer re-derives or refutes
each, and a refutation is first-class output):

1. `Field.where_largest` and `Field.cone_violations` share the flat-scan /
   argsort / unravel shape — a possible Pattern-2 twin inside the carve's own
   additions (one parameterized magnitude-map primitive?).
2. Residual expressiveness (the user's own lead): the SI stop and record ride
   flat `_l2_norm` spellings by design (conv_tol semantics; the retired
   doctrine's own note said switching the STOP to `.l2` re-interprets the
   tolerance) — re-weigh under principled-over-bit-identical whether the typed
   residual algebra should now carry more of it, and whether
   `affine-typed-residual-eq`'s kept label matches what ships.
3. Cone properties at consumers: `power_iteration`'s ray normalization, the BC
   realizer's cone refusal, the #390 flag-reader — is there a coherent
   cone-vocabulary surface rather than three idioms?
4. `CoefficientRole` is now an empty marker whose complement-content died —
   re-run type-vs-property on its two remaining justifications (taxonomy
   isinstance; future multiplier home).
5. DSA's `apply(self, displacement: ...)` parameter name and any other
   vestigial `displacement`/`disp`/`⊖`/mint/base-point vocabulary in live
   signatures, locals, and test names (grep the CONCEPT's spellings, not
   "torsor").
6. The Krylov `to_flat` role erasure: under V a typed Krylov iterate is newly
   representable — likely #289/#288 territory (out of CS3-R scope to BUILD;
   in scope to STATE the now-cheaper design).
7. `tests/sn/operators/test_declared_law_is_linear.py`'s long header was
   part-rewritten in place — full-coherence re-read.
8. The composite delegation notes and `_principal_bulk_leaf`'s walk — could the
   composite own a `principal_leaf` property instead of numerics duck-walking
   it (native-place question)?

**Done when** (checkable): a findings memo lands (scratch/, incrementally
written) with every Phase-1 candidate carrying a Phase-2 verdict + structural
reason; in-session-fixable items are FIXED inline (process-discipline: no
open/close issue noise); cross-session items are FILED with module labels; the
completeness verdict states its denominators (which spellings were grepped,
which smell sweeps ran, over which trees); and this section gains its ✅ row in
the ledger.

**CS3-R execution ledger:**

- ✅ CS3-R ⏹ COMPLETE (2026-08-19, branch `refactor/cs3r-cone-review`, six
  sweeps `755f99b5`/`77c7cc68`/`a740d7ba`/`d79adb27`/`8bde369b`/`8112e1b0`).
  Findings: `scratch/cs3r_findings.md` + the independent census
  `scratch/cs3r_census_qa.md`. **Completeness was REFUTED at review start**
  (~20 live-tense survivals: the RC family's ⊖ prose incl. 2 error messages,
  DSA's noun → *increment*, 5 gates' false "illegal" premise, 3 AGENT.md
  briefs — one PRESCRIBING the torsor as the fix — 1 skill clause of
  derivative staleness, 3 test-architect memories, 2 production docstrings my
  case-blind grep missed and the census caught) and is complete against two
  independent filters at close. Expressiveness: 5 gates re-spelled to the
  direct law (`[M]` blend-form 2.2e-16 BLIND vs direct 1.168 LOADED on an
  affine mutant); `principal_bulk_leaf` polymorphic on the carriers (numerics
  walker + `_NormedLeaf` retired); `Field._index_tuples` single-sourced;
  leads 2/3/4/7 WITHDRAWN with structural reasons in the memo. Filed #391
  (member-contract collapse — CS3 dissolved its recorded blocker), #392
  (`face_areas` shim, surfaced by the gate's warning escalation); commented
  #289 (typed Krylov now representable). `[M]` wide slice 4908 passed /
  0 own reds (the 1 red = the escalated DeprecationWarning on the
  pre-existing #392 shim; canonical invocation green); capture gate 5/5
  frozen; pyright `orpheus` = 1 (the ratchet's recorded state, 0 mine);
  sphinx clean ×3; `dead_references` 0/56 at close.

## 5. Phase CS4 — S and F are kernels; operators are born bound

**Goal.** `ScatteringKernel`/`FissionKernel` are representation-free data; a constructor
binds Kernel × Frame → the actual operator, fully bound (one domain, one codomain, ONE
public `apply` signature); apply-time dispatch retires; the binding carries its
correctness gate.

**Proposed means** (hypothesis): the kernels wrap what already exists
(`LegendreMomentScattering` + `mat_xs`; the rank-1 χ⊗νΣf dyad); the constructor mints
the frame from (the kernel's eigenbasis STRUCTURE — Funk–Hecke — × the space's angular
measure) — ruling D2's reconciliation: the frame remains the operator's eigenbasis
(preserving `frame-eigenbasis-ownership`'s content) while the SPACE supplies the measure.
C stays a multiplier binding space-only (no frame — its data object is the coefficient
field itself). The cross-method arms (S's `ScalarFlux` entry, the ndarray hatches —
in-code keep rulings #205/#276, `scattering.py:1172`) become separately-bound arrows or
explicit adapters — design decision at carve time with those rulings re-read, not
silently dropped. `IntegralKernelOperator`'s Protocol + category test re-scope to the new
split (ruling D3). **The tightness gate (report 1.4, the one genuine Phase-1 gap — `[M]`
0 hits at HEAD) lands WITH the binding**: `bind(K)† = bind(K†) ⟺ the frame is tight`,
with a deliberately non-tight frame as the ERR-039-class negative control.

**Done when**: the remaining P1 ledger deletes — R1 annotation ×5 and R2 refusal ×3
(S/F/C at the binding; L/B ride the CS2 SN space) — for the campaign total of **all 12
P1 strict-xfails deleted**; `[M]` no `singledispatchmethod` on S/F `apply`; the non-tight
negative control REDs; #359's three M-M spellings re-checked against the new binding
(unify or re-file with the measured residue).

**Surgical posture.** This is the operator-algebra carve family: the MAIN AGENT writes
with the user steering; `method-implementer` NOT dispatched (`delegation.md`);
`test-architect` dispatched BEFORE the carve (proactive trigger); the L17 convention
crosswalk written to `.claude/plans/` before code.

---

## 6. What this campaign does NOT do (Campaign 2 + parked)

- Everything LOWERING: the LossRepresentation overturn, `GeneralizedEigenPencil`,
  resolvent machinery + complexification (report S5 → Campaign 2 prerequisite),
  partitioning/schedule/BlockDigraph, SCC wiring (#324/#320), the monodromy close (#300,
  blocked on #299), the σ_r honest home. The grounded state of that territory:
  report Part VI.2.
- The Phase-C battery beyond what CS3 itself needs (its value-level subset — CW brackets,
  SI-MONOTONE, MONOTONE-SOLVE — attaches to the #340 seam any time post-ruling; the
  step/SC witness question D6 is decided in that planning).
- The doc workstream beyond D7 (cone chapter, rides CS3) and D8's essentials (spaces
  chapter, rides CS2).
- Carried ordering constraint: **O-4** (independent α reference, ~15 lines
  `scipy.linalg.eig` on a hand-posed G×G pencil, + a velocity FIXTURE — the data layer
  carries no `v_g`) before F migrates in Campaign 2's posing work.

## 7. Risks and controls

- **Enabler blast radius** (plan-authoring §8): every domain that stops being `None`
  activates guards and metric-`.H` paths that were silently skipped — enumerate per
  phase, gate at the observable tier (the ERR-076/ERR-067 defect family).
- **Metric moves** (F3 → spaces): each volume/weight relocation owes an adjoint-path
  check (ERR-067 class) and bit-identity or documented principled-equivalence.
- **CS3's re-pointing** owes the full retirement discipline: 3-search audit per symbol +
  concept-grep, marker migration, message-fragment pins, and the labelled-equation rule
  for `field_algebra.rst`'s 4 labels.
- **Identity flip** (CS2/S3): logged-then-flip intermediate mode; inventory the aliasing
  sites before the flip (the surgical-edit discipline); #369's registry fence is the
  measure-side precedent.
- **Compaction points** at every phase boundary (≥4 phases ⟹ mandatory); commit-then-
  checkpoint; this file is the plan of record and is a LIVING document under
  `plan-authoring.md` — surprises get logged here in place.
