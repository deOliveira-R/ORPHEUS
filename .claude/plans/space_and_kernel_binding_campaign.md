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

## 1. Dependency readiness — what the tree provides (verified 2026-08-19 vs `7aae9bf1`)

| capability | state | evidence |
|---|---|---|
| `FunctionSpace` is carrier-generic (PEP-696 `Carrier` TypeVar) | shipped — the surface CS2 builds beneath | memo A §3.F5 (`space.py:92-100`) |
| Frame machinery `FrameBase → PetrovGalerkinFrame → GalerkinFrame`; `frame.conjugate(Λ)` IS the production R∘Λ∘M | shipped — CS4's binding verb exists | memo C §1 (`frame.py:114/:205`, `scattering.py:488/:599`) |
| Moment-space cores holding only XS data (`LegendreMomentScattering`, `N2NMomentOperator`; F's rank-1 dyad) | shipped — the Λ half of the kernel split already separated | memo C §2 (`scattering.py:115/:300`, `fission.py:258`) |
| Q6 measure metadata (`ReferenceMeasure`, `ExactnessClaim`, `invariance_group`, `half_range_clean`) | shipped — the GENERATOR side of the Angular axis; CS2 FORGETS from it | memo A §4 S2 row |
| `DiscreteMeasure.partition_by` | shipped (2026-05-10) — the measure-level partition primitive the per-axis hooks forget from | memo A §4 (`measure.py:1042`) |
| The P1 ledger: 12 `xfail(strict=True)` markers + `tests/sn/architecture` | shipped; `[M]` 105 passed / 21 xfailed reproduced 2026-08-19 | memo B §8.8 (node ids) |
| #340 record machinery (`IterationRecord`, `StoppingCriterion`, named criterion trajectories, `IterationBudget`) | shipped — where CS3 relocates the displacement diagnostics | memo B H3/§6.1 |
| Cone fragments: `is_positivity_preserving` (DD `False` + numerical witness), ψ≥0 realizer refusal, coefficient cone-as-predicate battery, ray normalization in `power_iteration` | shipped | memo D §3 |
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
`test_model_generic_leaf_declares_a_space[C-2g, C-4g, F-2g, F-4g]` are DELETED on XPASS;
`homogeneous/solver.py`'s `F.domain is None` state and the `basis_shape=(ng,1)`
double-pass are gone (`[M]` grep).

**§8 blast-radius note (enabler ≠ neutral).** Giving F/C a real domain flips
`OperatorSum`'s composability guard from *skipped* to *active* on every composition that
today rides the `None`-skip (`operator.py:582`), and `.H` from Euclidean-fallback to
metric-applied (`operator.py:1221-1226`). Enumerate the newly-checked compositions and
the `.H` consumers on the homogeneous path BEFORE landing; gate at the tier the change is
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

**Why this phase sits between the spaces and the binding** (sequencing rationale, mine —
re-orderable at the first compaction point): CS4's collapsed `apply` signature must be
written against the settled field algebra, or the new bound operators bake in the
flux-vs-displacement refusal (#331's fork) that Campaign 2's Krylov/resolvent work would
immediately re-open; and running CS3 after CS2 means the re-typed fields land on the
final spaces (one pass over the field layer, not two).

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
