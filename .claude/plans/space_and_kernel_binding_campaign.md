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
   space phases — see the §4 sequencing block for the derivation. ⚠ 2026-08-20: CS1.5
   (the Medium un-weld) inserted between CS1 and CS2 — chartered in the CS1 design
   record (its Q11 ruling carries the reasoning).
   ⚠⚠ 2026-08-20 (user, later the same day): **the order is RE-RULED — the kernel
   carve's energy core pulls forward**: CS3 ✓ → CS1 ✓ → **CS4a** (kernel types +
   construction-time binding for the energy family C/IsoS/IsoN2N/F; S-angular's
   frame and the L/B xfail rows still ride CS2 — §5's own done-when already
   carves them out) → **CS1.5′** (medium/assignment, designed against the LANDED
   binding signature) → CS2 → **CS4b** (the apply-time dispatch retirement
   remainder). Reason (the pinned dead-end; user: "the sequencing might just
   have been a bad charter on my part"): the CS1.5 design kept orbiting
   carrier-interface questions — where `mat_map`/`spatial_shape`/`volume_measure`
   live — that are ARTIFACTS of apply-time-polymorphic operators reading a
   carrier; ruling 2's construction binding dissolves all six of the homogeneous
   path's carrier demands (`[M]` the demand→home mapping in
   `.claude/plans/kernel_and_medium_objectives.md`). The charter error was the
   maximalist reading of "spaces before binding" — ALL spaces before ANY
   binding — when a binding needs only ITS space, and the energy-family spaces
   ship since CS1.
   ⚠⚠b 2026-08-20 (post-selection, same day): the adversarial design rounds +
   rulings F1–F4 (synthesis of record `scratch/cs4a_synthesis.md`; operative
   charter = §5's CS4a/CS4b/CS4c subsections) subdivide the arc and re-place
   the medium: **CS4a** (kernel core — ScatteringKernel/N2NKernel/
   FissionKernel + construction bindings on the energy four; NO apply-arm
   deletion — the R-A per-instance feeding census refutes space-keyed arm
   selection today, vv-principles #29) → **CS4b** (fields constructed from
   SPACE — the third weld the rounds exposed; the fabricated path retires
   HERE, where it becomes structurally unneeded) → **CS1.5′** (Medium) → CS2
   → **CS4c** (the dispatch collapse: feeding normalization, #205
   re-litigation, arm deletion, C-mandatory with its 131/43-site migration,
   S→kernel shell). Phases subdivide further as needed (user: "as many as
   necessary to chunk this work into targeted sessions with compaction in
   between").
   ⚠⚠c 2026-08-21 (user, post-CS4a): a **CS4a-R review round precedes CS4b**
   — "before CS4b, let's do a round of review of the kernel and operator
   machinery with cleaner context." The order reads CS3 ✓ → CS1 ✓ → CS4a ✓ →
   **CS4a-R** → CS4b → CS1.5′ → CS2 → CS4c. Charter = §5's CS4a-R
   subsection; precedent = §4-R (CS3-R).
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

**Design record (discussion CONCLUDED 2026-08-20, rounds 1–6 — every fork ruled or
resolved):** `.claude/plans/cs1_energy_space_design.md` is the CS1 design record AND
the post-compaction resume surface — read its "▶ THE PROTOCOL FROM HERE" section
FIRST. Headline rulings: EnergyAxis (identity = ng + edges content; counting-measure
theorem); explicit `(ng, 1)` quotient point (the collapse doctrine, Appendix A =
archivist's source); axes-tuple order = the recorded view convention (#394 filed for
the JAX-era layout freedom); Q8 = CS1 rewrites the space.py identity doctrine in
migration form; **Q11 = CS1.5 chartered — "the Medium un-weld" — a NEW phase between
CS1 and CS2** (Medium minted ABOVE MaterialMesh, which keeps its name as the
medium × mesh PULLBACK; homogeneous re-homes to `Medium.infinite_homogeneous(mix)`;
`bulk_space` migrates carrier→Medium there; the fake `AxisMesh(edges=[0,1])` mint
dies there and ONLY there). **ABSORBED 2026-08-20** — protocol step 2 (clear-context
re-evaluation + §P grounding) executed at `273d431a`; the record's §P results block
carries the evidence; the record stays the full doctrine, this section is now the
executable summary:

**Ruled means (CS1 ships):** `orpheus/numerics/axis.py` — `BasisKind`, frozen `Axis`
(structural eq/hash; label + shape rank ≥ 1 + `weights | None` ≡ counting-deliberate
+ NODAL/MODAL), `EnergyAxis` (`from_grid`/`synthetic`; identity = ng + edges BYTES —
`EnergyGrid` is `eq=False`, so the axis reads content, never grid identity; docstring
carries the counting theorem, the V/V* hook, the faces reading, the descending
fast-first convention). `FunctionSpace` gains `axes` (compare=False), `of_axes`
(deterministic INJECTIVE names — `(name, shape)` is the identity until S3),
the per-axis metric path, `has_coordinate_cone`, `__mul__` axis-threading; the
`space.py:145-150` identity doctrine rewritten in MIGRATION form (Q8).
`MaterialMesh.bulk_space` — the UNIFORM cached mint
`of_axes(energy_axis, Axis("spatial", spatial_shape, weights=volumes, NODAL))`,
honest on EVERY member (degenerate carrier ⟹ quotient point, weight 1.0; meshed
carrier ⟹ the CS2 scalar-bulk seed, volumes as weights); energy arm = the materials'
unanimous content-equal grid, else `synthetic(ng)`. Homogeneous rewiring: C rides
the extended `from_mesh` chain (`bulk_space` fallback AFTER `full_field_space`, so
SN/diffusion short-circuit untouched); IsoS/N2N/F get `space=V`; F's field widens
`FullFieldSpace | None` → `FunctionSpace | None` and renames → `space` (the grounded
§6b batch lives in the record; ⭐ RULED 2026-08-20: name = `space`, scope = **F+S** —
S renames, its type stays narrow);
`Field.cone_violations` consults `has_coordinate_cone`. Execution: steps 0–5 with
step 3 SPLIT — 3a rename/widen (behavior-neutral) then 3b wiring + xfail handling
(§6b-fused: the chain extension alone XPASSes the C rows). Fences: no
Spatial/Quadrature/Harmonic axis classes, no ⊕, no identity flip, no metric
relocation or `_broadcast_metric` deletion, no Optional→mandatory flip, and no
`from_materials` internals (CS1.5 owns the fake-`AxisMesh` retirement).

**Done when** (checkable): the four model-generic strict-xfails
`test_model_generic_leaf_declares_a_space[C-2g, C-4g, F-2g, F-4g]` are DELETED and
replaced by the positive homogeneous floor — ⭐ CORRECTED 2026-08-20: only the **C**
rows XPASS mechanically (forcing deletion); the **F** rows cannot (their arm calls
`from_solver_data` BARE and no-default-derivation is RULED), so they are deleted in
the same §6b commit on the retired-mirror warrant
(`tests/sn/architecture/test_monomorphic_leaves.py:670`, re-verified at `273d431a`;
ledger `[M]` 105 passed / 21 xfailed reproduced 2026-08-20, 1.9 s; homogeneous suite
23 passed, 2.4 s — the byte-stability budget);
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
observable. ✅ ENUMERATED 2026-08-20 (design record, §P results items 2/3): all
construction sites classified (three CS1-edits, the rest threads-already or
documented stays-None); NO existing test asserts homogeneous-op `.H` values, so the
vv#19 loaded/blind pair is new coverage with nothing demoted; the metric flip is
×1.0/÷1.0 or a counting skip — bit-identical by construction, gate still owed.
**Verification battery of record (protocol step 3, 2026-08-20):**
`scratch/cs1_verification_plan.md` — 40 gates, 23-mutation battery with
MUST-STAY-GREEN anti-claims, reviewed + adjudicated in the design record §T-R
(findings F1–F12 folded; the seven open rulings resolved; ⛔ headline: a scalar
metric is a provable `.H` non-catcher, so the vv#19 control is a per-GROUP
weighted toy and the cell measure is guarded by space identity alone).

### ⏸ COMPACTION POINT #2 — CS1 ⏹ EXECUTED (2026-08-20, steps 1–5 on `feature/cs1-energy-space`, unmerged at write time — trust git for merge status)

| step | commit | shipped | teeth run |
|---|---|---|---|
| 1 | `1afff47b` | `orpheus/numerics/axis.py` (BasisKind; frozen Axis, structural eq/hash PER SUBCLASS; EnergyAxis identity = ng + edges CONTENT; Q-T1/Q-T3 canonicalization: all-ones→None, −0.0 killed via `+0.0` which also forces the defensive copy, non-finite refused, NO non-negativity guard; weighted/MODAL EnergyAxis refused — the counting theorem as a construction invariant) + `tests/numerics/test_axis.py` A1–A12 | M11/M12/M20/M1/M2/M10 exact incl. M2's A3-stays-green anti-claim |
| 2 | `f4876354` | `of_axes` (derived name = the S3 identity bridge — blake2b over type-tagged length-prefixed `Axis._structural_bytes()`, injective + process-stable); per-axis metric path (exact interior-axis placement; pairing single-sourced through `apply_metric` so ERR-067 is unspellable); 3-valued `has_coordinate_cone`; `__mul__` threads axes, never fabricates, and BRIDGES axis-borne measures into mixed dense products (value-correctness; CS2 retires the densifier); `dual()` threads axes; two construction guards (one metric source; shape ≡ concat); Q8 migration-form doctrine in `space.py` | M3–M9 + M20′ exact incl. M5's B5/B6-stay-green (densified metric is value-identical) and M4's subprocess-leg-only |
| 3a | `e8769897` + amend `24a991ba` | Q5 executed: `full_field_space` → `space` on F and S; F widened to `FunctionSpace \| None`; **S ALSO widened (user, 2026-08-20, superseding the narrow-S clause — `[M]` C/IsoS/IsoN2N already carry the wide type; the record's Q5 row carries the supersession)**; 8 keyword call sites; mesh properties + `FullFieldSpace` class keep their names | behavior-neutral: 1939+202 passed; pyright ratchet unchanged |
| 3b | `6bd782ab` | `MaterialMesh.bulk_space` (uniform formula; energy arm mat_map-reachable; CS1.5 re-point comment); `from_mesh` chain + D8 re-words; solver rewiring (IsoS/N2N/F get `space=V`; BOTH production `(ng,1)` spellings deleted — the M23 grep obligation discharged in the commit); 4 xfail rows deleted (C forced by XPASS, F on the retired-mirror warrant) → D1 floor + D2/D3 witnesses + D4a/D4b vv#19 pair + D6/D8/D9/D11 (`tests/homogeneous/test_operator_spaces.py`) + B9/B10/B11 + D7 poison path gate + D12 mirror re-points + **D5 byte gate 8/8 BIT-identical across the rewiring (the counting theorem measured; fixture captured at `24a991ba`; migration gate — retire after the merge cycle)**. ⚠ Q-T5 correction: `[M]` only region A is producing (the record's "A/C" enumeration was stale; the `is_producing` screen governs) | M13–M19+M22: M13 control 6 reds; M17 → B9/B10 ONLY with D5/D4a/D1 GREEN (the F2 proof); M19 → D4b ALONE with D4a green (vv#19 measured); M14's homogeneous-stays-green |
| 4 | `6da1b23c` | `Field.cone_violations` consults `has_coordinate_cone` (False ⟹ typed refusal naming space+reason; None ⟹ legacy, documented) + E1/E2/E3 | M21 control: E1 red ALONE |
| 5 | (this row's commit) | `spaces.rst` seed from Appendix A (dialectic included) + toctree; `api/numerics` automodule for axis; `field_algebra.rst` cone-consult micro-edit; this ledger row | Sphinx `-W` clean |

Full sweep at step 4: numerics+homogeneous+sn/architecture+transport+diffusion
**3113 passed / 17 xfailed** (ledger: the 4 R1 value rows gone, annotation family
intact per F12); pyright `orpheus` = 1 (ratchet terminal state) throughout.
**Done-when reconciliation:** (1) ✅ xfails deleted on the two warrants, floor
shipped; (2) ✅ `[M]` grep at `6bd782ab`: zero production `(ng, 1)` spellings
(docstring prose only); (3) ✅ D2/D3 red-capable on the established fragments +
E1 modal witness; (4) ✅ byte-stable 8/8; (5) ✅ this row + the step-5 docs.
**NOT done by CS1 (fences held):** no Spatial/Quadrature/Harmonic axis classes,
no ⊕, no identity flip (S3), no metric relocation (CS2), no Optional→mandatory
(CS4), `from_materials` internals untouched (CS1.5 owns the fake-AxisMesh death).

**⭐ BRANCH-HOLD RULING (user, 2026-08-20):** `feature/cs1-energy-space` is HELD
unmerged — *"We will significantly advance before the next merge to main"* — so
the following phases (CS1.5 → CS2 …) land on THIS branch, and the ≥90-min
pre-merge full gate runs once, before that eventual merge, not per phase. A
session picking this up must NOT merge after its own phase without the user's
go. Per-phase verification stays as CS1 practiced it (scoped suites + battery +
pyright + Sphinx `-W`). Merge-status: trust git, never this sentence.

## 2.5 Phase CS1.5 — the medium is a concept with its own home (the un-weld)

**⏸ RE-SEQUENCED 2026-08-20 (§0 ruling 3's ⚠⚠): CS4a — the kernel core — now
precedes this phase, and this section's grounded design (with the four fork
rulings of 2026-08-20) is DEMOTED to a CONTESTANT in the joint kernel+medium
design rounds.** The carrier-side machinery below (the union's quotient member,
`from_medium`'s quotient arm) was interface design for demands CS4a dissolves;
the physics-side objectives (medium expressible, generator lattice, conformity,
retirement honesty) survive, restated in
`.claude/plans/kernel_and_medium_objectives.md`. The grounding census and the
`[M]` fact base remain the shared evidence. Body kept per plan-authoring §3.

**Goal.** The physical medium (what fills space) exists as its own type ABOVE
the discretization; the homogeneous solver poses on it; the fake
`AxisMesh(edges=[0,1])` mint is unspellable. CS2's axis generators then bind to
the right generators from day one (`SpatialAxis` forgets from the MESH, the
quotient axis from the MEDIUM's symmetry — Q11's own reason).

**Charter**: `cs1_energy_space_design.md` Q11 + §B2 rounds 4–6.
**Grounding + rulings + step decomposition (the §F-equivalent)**:
`.claude/plans/cs15_medium_unweld_design.md` (⏹ GROUNDED + 5 forks RULED
2026-08-20; census of record `scratch/cs15_grounding_census.md`).

**The grounded shape** (headline — the design note is authoritative): `Medium`
= the NAMED PAIRING `(InfiniteHomogeneous marker | StructuredGeometry) × {id →
Mixture}` — the loose pair every solver already passes, REUSING the shipped
`Region`/`StructuredGeometry` (the grounding surprise: §B2's "interface
positions + per-segment materials" already ships at
`geometry/structured_geometry.py`); pullback union
`spatial_structure: MeshBackedRegions | QuotientPoint` inside `MaterialMesh`
(`[M]` all 10 spatial writers already funnel through `_init_data`; type-absence
is production-clean); `from_medium` both arms + conformity guard with its §6c
witness; eg-coherence refusal at Medium construction **with the user's regime
caveat** (right for coarse/GENDF multigroup, NOT universal — fine/ACE data
legitimately unequal, unionization is the modern reconciliation; carried into
docstring + theory page).

**Steps** (design note §5): 1 mint Medium (additive) → 2 union restructure in
place (from_materials keeps its signature, stops faking — §6b clean; byte gate
8/8) → 3 from_medium + re-home + retire from_materials (ONE step by §6b: 1
production + 16 test sites in 7 files + the `infinite_medium.rst:1115` xref) →
4 docs + the census §7 promise ledger (8 items, incl. the unmarked D6/D7
re-points).

**Done when**: the promise ledger is discharged item-by-item; `grep -rn
"from_materials" orpheus/ tests/ docs/` returns only past-tense history; the
conformity guard's first red (hand-built non-conforming `Mesh1D` REFUSED with a
region-naming reason) is a committed witness; D5 byte gate 8/8 bit-identical
across the rewiring; scoped suites + pyright terminal-1 + Sphinx `-W` clean.

**CS4a-R amendment (2026-08-21, XD-4 — attribution in
`scratch/cs4a_r_findings.md`):** the Medium's grid-coherence invariant
should realize the THREE-outcome structure the energy-arm rule currently
collapses to two: *agree* → `from_grid`; *no grid anywhere* → `synthetic`;
*disagree* → REFUSE (or refine — `EnergyGrid.overlap_to` already ships the
refinement maps). Today `EnergyAxis.synthetic` stands for BOTH "the problem
has no grid" and "the materials disagreed" (`[M]` indistinguishable
outputs), the disagree arm is first-wins order-dependent (`[M]`
`[2g,4g]`→(2,), swapped→(4,)), and `[M]` two same-partition grids differing
at the last ulp silently lose their grid (exact content-equality is the
RIGHT relation — transitive, deterministic; the near-equal case is evidence
the disagree arm fires on real multi-provenance data). Unreachable from
production today (the carrier's ng-unanimity gate + singleton callers —
the review's MA-1 withdrawal); the refusal belongs HERE, at the Medium,
not at the classmethod.

**Design input (user question, 2026-08-21) — `SNMesh` is a misnomer, and
the remedy is the DECOMPOSITION this phase + CS2/CS4b already charter, not
a rename.** The class is a conjunction: geometry mesh (spatial measure +
BC declarations) × materials (energy measure) × assignment (`mat_map` —
not a measure; the coupling that makes Σ a function ON the phase space) ×
quadrature (angular measure) × scheme (the within-cell closure that
refines the spatial factor) — `augmented_mesh.py`'s own module name
confesses it. Ruled reading: the phase space PROPER is the object the
class MINTS (post-CS2 literally `of_axes(Angular, Energy, Spatial)` ⊕
trace), so ⛔ do NOT spend the name `SNPhaseSpace` on the container — the
Medium takes materials+assignment (this phase), the minted space takes the
phase-space role (CS2), fields consume the space (CS4b), and the residual
— mesh + quadrature + scheme (+BCs on the geometry) — is the
DISCRETIZATION-choices tuple (`SNDiscretization` if a class survives; it
may dissolve into the solver entry). Decide the residual's name when the
un-welds land, not before. Formal note for the BC constituent: a BC is a
law on the trace pair (Γ₊, Γ₋) whose measure is the SPATIAL measure's
boundary restriction (dV → dA, the divergence-theorem counterpart) × the
`|Ω·n̂|`-weighted ANGULAR measure (which also defines the Γ± split) ×
energy (diagonal; an albedo may couple groups) — already realized as the
trace summand's metric `G_trace = |Ω·n̂|·w_n` in `full_field_space`.

**⭐ Refined by the cross-domain tournament (2026-08-21; memo =
`scratch/sn_posing_ontology_and_names.md` — every claim file:line-grounded;
main-agent spot-verifications below marked `[M]`):**

* **The criterion, operational form:** the residual owns exactly what
  determines the **DOF set and its Gram** — removal test: *did ψ's shape or
  G change?* YES ⟹ carrier; NO ⟹ operator machinery. Its ONE irreducible
  behavior is **mutual admissibility** of the choices. The FEM borrowing
  that re-derives it: *a datum belongs to the carrier iff it CREATES DOFs*
  (an FEM quadrature is an integration rule — no DOFs — while SN's angular
  quadrature is a collocation: the ordinates ARE DOFs).
* **The scheme-vs-frame symmetry INVERTS rather than levels:** frame ↔
  CLOSURE (both operator-realization data; the closure should leave for
  the same reason the frame already did), quadrature ↔ ELEMENT (both
  DOF-creating; the element face — `spatial_basis_per_axis`,
  `moment_mass_diagonal`, θ — stays for the same reason quad does). And
  the scheme does NOT fission into two classes: both faces read ONE θ
  (the closure is a Galerkin statement about the element) — the move is
  mint-and-forget, one scheme object with two READERS (the space mint
  reads the element face; the streaming binding reads the closure face).
* **Forgetfulness-after-binding is uniform, with one limit:** every
  constructed object keeps the arrows its contract needs and forgets its
  minter — valid only when the retained set is CLOSED under the
  operations later demanded ("the arrows plus the laws they satisfy";
  the PG-frame adjoint is the counterexample class).
* **The tree's real inventory exceeds the measures by six objects** (the
  reduced streaming stencil, `_streaming_axes`, the realized-BC operator
  table, `pole_angular_closure` — a bound strategy with a BACK-REFERENCE
  to the mesh — `curvature: str` beside `coord`, the law registry, plus
  minted spaces/gauge/dag_walk): four of the five L-machinery residents
  are the StreamingOperator's, misfiled on "the phase space".
* **Name tournament verdict:** `SNPose` ⛔ KILLED (`[M]` "posing" is spent
  on the `(A_loss, M)` eigenvalue arrangement — `eigenvalue.py:23`);
  `SNRealization` ⛔ KILLED (word spent on operator realization AND
  direction inverted — this object is realization's INPUT);
  `SNProblem`/`SNModel` ⛔ (falsified by the Medium extraction);
  **winner `SNDiscretization`** (names the invariant ROLE; term of art;
  degenerates correctly to `DiffusionDiscretization` — whose residual
  under the criterion is `axes` alone, `[M]` `_init_core` is a documented
  twin). ⚠ CONDITIONAL on the **funnel ruling** (user's, at CS1.5′/CS2
  grounding): if the space mint is the unique funnel every construction
  passes through, the four cross-choice refusals live there and **no
  class** beats every name (the residual dissolves into
  `solve_sn(medium, mesh, quadrature, scheme)`); if not, the class is
  worth ~40 lines and exactly ONE invariant (inadmissible combinations
  unspellable) — there is no middle version that is not the God object
  again. `[M]` the invariant is live: `SNMesh(slab_1D, lebedev)`
  CONSTRUCTS today (110 S² ordinates on a slab) while
  `tests/numerics/test_registry.py:1210` asserts
  `not slab.admits_domain(lebedev)` — filed as an issue.
* **⭐ The §8 mechanism (why the mesh became the God object): IDENTITY
  SCARCITY.** Things accrete to the mesh because `FunctionSpace.__eq__`
  is `(name, shape)` — too weak to key caches or gauges on (the tree
  SAYS so at `augmented_mesh.py:986-994`). ⟹ sequencing constraint for
  CS2, recorded in its amendment block.

## 3. Phase CS2 — SN's space is a composition of Energy, Spatial, and Angular

**Goal.** SN's bulk space is literally the composition Energy ⊗ Spatial ⊗ Angular of
first-class per-axis objects — one composition mechanism, structural identity, metrics
owned by the spaces, per-axis partition hooks recorded for Campaign 2.

**⭐ CS2 design inputs (user + main agent, 2026-08-20 — three adjudications from
the adjoint/binding discussion; each supersedes the older means-sketch where
they conflict):**

1. **The SN adjoint does NOT carry the dyad trap — keep the two-layer
   architecture it proves.** `[M]` verified at `87fde3f7`: leaves ship the
   honest EUCLIDEAN transpose (fission's composite arm: unweighted Σ/W then
   `Kᵀ` then w-broadcast — `fission.py:429-460`), and `.H` conjugates with the
   space-owned metric (`_AdjointOperator.apply` = `G_V⁺ ⊙ apply_transpose(G_W
   ⊙ y)`, `operator.py:1290-1314`; SN interior `G_bulk = V_cell·w_n`,
   `augmented_mesh.py:845/:1102`) — composing to the metric-correct form
   (weighted-integrate/constant-broadcast, factors swapped). Gated at the
   discriminating tier: `test_g_adjoint_reciprocity.py` (⟨Aψ,φ⟩_G duality with
   an INDEPENDENT g_inner + the |Ω·n̂| negative control),
   `test_fission_adjoint.py` (`composite_pairing_identity`,
   `weight_swap_discriminator` — literally "WHICH side carries the quadrature
   weights"). ⚠ The surviving nuance: the leaf transpose arm reads `w` off
   the OPERAND (`bulk.mesh.quad.weights`) — correct today, but the same
   third-weld read vv#30 catalogues; under input 2 it becomes a SPACE read.
   ⚠ **The user's price assessment (2026-08-20): the wrapper is a PAID cost
   of the dyad discipline, not a free lunch — a target for EVENTUAL
   reconsideration.** Its retirement condition, stated now so the target is
   concrete: under the factored composite realization (inputs 2–3 + the F4
   addendum), a channel's `.H` becomes the reversed composite of
   space-supplied adjoint PAIRS — metric-exact with NO G-conjugation,
   because the conjugation is absorbed into the pair construction. Each
   channel that moves to the factored form exits the wrapper's clientele;
   the wrapper retires (aggressive-retirement) when its clientele empties —
   monolithic Euclidean leaves (e.g. the #280/#310 transpose sweeps) are
   what keep it alive until then.
2. **Axes are measure-ACCESSIBLE, never measure-WELDED (user refinement,
   2026-08-20, superseding this input's first spelling "measure-complete").**
   Storing the measure OBJECT (or the Quadrature) on the axis would weld
   space to measure and abandon the forgetful-map concept. The ruled shape —
   the generalization of what `EnergyAxis` already does — is a three-part
   law: **(a) forget the GENERATOR** (`from_grid` consumes the grid's
   surface only — edges, n_groups — under a TYPE_CHECKING import; the axis
   never holds the grid; same for `AngularAxis = forget(Quadrature)`: the
   registry tags, ordering conventions, construction provenance are
   forgotten); **(b) store the IDENTITY CONTENT that distinguishes the
   measure** (energy: edges bytes — implemented, ruled Q2, precisely so two
   axes claiming the same ng gate operations by grid content; angular: the
   NODES by the same argument — two rules with equal weights and different
   directions are different measures; spatial likewise — "a fair
   consideration for the other Axes", user); **(c) expose an ACCESSOR that
   mints the measure object on demand** (`EnergyAxis.grid() → EnergyGrid` —
   rebuilt from stored edges; numerics→data is layer-legal, "any layer →
   input"; `AngularAxis.measure() → DiscreteMeasure` — numerics-native; the
   return-type asymmetry is honest to each axis's nature). The ℓ≥1 moment
   maps need `P_ℓ(Ω_n)` — the nodes — so weights-only is insufficient
   regardless; with (b)+(c) the scattering binding collapses to
   (kernel, space) — the separate `quadrature=` argument (a second source
   beside `space=`, a latent disagreement channel) disappears, and ruling
   2's "frame minted from kernel eigenbasis × the SPACE's angular measure"
   is satisfied through the accessor. ⛔ This supersedes the means-sketch's
   "directions stay generator data" clause below (the DIRECTIONS-as-identity
   move; the generator itself stays forgotten, so the clause's intent —
   don't weld the quadrature object into the axis — survives in (a)).
   **⭐ Deepened (user conjecture verified, 2026-08-20): the general
   structure is ALREADY HALF-BUILT — all three generators natively speak
   `DiscreteMeasure`** (`EnergyGrid.as_measure()` at `energy_grid.py:193`,
   documented symmetric with `Mesh1D.volume_measure`; the Q6 quadrature
   rules RETURN `DiscreteMeasure` natively, `rules_1d.py:41`). ⟹ a generic
   **`Axis.from_measure(label, measure, kind)`** collapses the per-axis
   constructors, and the axis metric = the view's weights UNIFORMLY —
   energy's counting included (`as_measure` carries `w_g = 1` with the
   rationale in place: φ_g is already group-integrated). **The reduction's
   principled stopping point**: `as_measure`'s nodes are integer group
   INDICES with unit weights, so for energy the view is DELIBERATELY
   DEGENERATE as an identity — two grids with equal ng have identical
   views; that degeneracy IS the counting theorem. The law: **identity =
   measure-view content where the view is faithful (spatial, angular) + a
   supplement exactly where the theorem degenerates it (energy: the edges
   bytes, ruled Q2)** — EnergyAxis's specialization reduces to that one
   supplemented field with its stated reason; `from_grid` and the
   weighted-refusal reduce to the generic constructor (+ Pattern-4 guard on
   the direct path). **Prerequisite, `[M]` probed: `DiscreteMeasure.__eq__`
   RAISES today** (default dataclass eq over ndarrays,
   `ValueError: truth value…`) — content-equality must be minted on
   `DiscreteMeasure` (the `Axis._structural_bytes` discipline generalizes
   verbatim: nodes + weights + support bytes) before measure equality can
   gate anything; a latent hazard worth closing regardless.
3. **The streaming/boundary formalization criterion (the user's, verbatim in
   substance): a FULL operator — bulk AND boundary action on
   `FullFieldSpace = interior ⊕ trace` — is what you get by binding to the
   DISCRETIZATION SCHEME; the scheme is exactly what connects bulk to
   boundary.** `StreamingKernel`'s representation-free datum is GEOMETRY
   data (the σ-free advection structure per coordinate system), and its
   binding is (kernel, space, scheme): the scheme (`DiscretizationSchemeBase`,
   already first-class on `SNMesh`) supplies the closure that couples
   interior to trace, and B's laws attach at the trace seam the scheme
   creates. Refinement: the scheme also CO-GENERATES the spatial axis — LD
   carries per-cell linear moments, i.e. a MODAL spatial factor
   (`BasisKind.MODAL` exists for exactly this), so scheme enters both the
   space construction and the closure. This is the L/B formalization the CS4
   sketch lacked ("L/B ride the CS2 SN space" now has its mechanism), and
   the binding-arity table per channel reads: energy (kernel, space);
   angular scattering (kernel, space→frame); streaming (kernel, space,
   scheme); boundary (law, trace-seam-of-scheme).

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

**CS4a-R amendments (2026-08-21, review round — attribution in
`scratch/cs4a_r_findings.md`):**

* **EE-7 — `FunctionSpace.inner_product` gains a SHAPE refusal.** `[M]` on a
  `(3,1)` axis-built space, ⟨x,(3,)⟩ silently returns 18.0 where the correct
  pairing is 6.0 (broadcast outer-product + sum), and a bare scalar operand
  is accepted. The admission half of the ERR-067 repair never landed; the
  homogeneous solver's `.reshape(ng,1)` calls are load-bearing today. Land
  the refusal on the raw surface when the composite overrides are designed
  (the same step that sizes the surface).
* **EE-8 — split the ng guard into its two natures**: the space-owned QUERY
  (`FunctionSpace.energy_extent` — the `has_coordinate_cone` shape; the
  thing the chartered axis-keyed strengthening needs) and the
  binding-owned POLICY (refuse + name the operator; moves into CS4c's
  binding base). Record the guard's two latent limits while touching it:
  it returns on the FIRST `EnergyAxis` (an undeclared choice — a
  domain≠codomain condensation space would carry two), and `(space, ng)`
  structurally assumes an ENDOMORPHISM (a 47g→2g condensation binding has
  no single ng conforming to both ends).
* **XD-3 — the chartered strengthening compares axis IDENTITY, not
  extent.** `[M]` the shipped guard reads `(extent,) = axis.shape` and
  never touches `axis.edges`, so a 2g kernel from grid A binds silently to
  a space on grid B — cardinality where the CS1 doctrine says CONTENT
  decides identity. Unspellable on the 4 live bindings today (space and
  data minted from the same mixture); becomes spellable exactly when
  arbitrary spaces bind. Design direction: the data side carries its own
  `EnergyAxis` through the hoisted rule, and the guard's strengthened form
  is one `==` on the axis value object.
* **Ontology-tournament facts (2026-08-21, memo
  `scratch/sn_posing_ontology_and_names.md`) — the identity-scarcity
  SEQUENCING constraint:** the `of_axes` derived-name flip (measure bytes
  in the digest) must land BEFORE any relocation of caches/gauges off the
  mesh, or they re-accrete for the same correct reason the tree already
  states (`augmented_mesh.py:986-994`: `FullFieldSpace.__eq__` is
  `(name, shape)` — a cache keyed there is keyed on a SIZE). Two red
  witnesses for the flip, both `[C]` at HEAD: two meshes differing only
  in a BC declaration mint `==` composite spaces; two schemes with
  different `moment_mass_diagonal` mint `==` composite spaces (the bulk
  block is hand-named `"sn_bulk"`). And `is_same_phase_space` demands
  `materials is` — problem-identity conflated with space-identity (the
  measure-vs-measurable-function frame's red witness) — re-scope it when
  the axis-built composite lands.
* **XD-7 — reserve `_agreed_space`'s stated expiry.** `operator.py:363-374`
  says in its own words that per-leg bindings make agreement the wrong law
  and "a product-space constructor is what has to be built" — CS2's
  frame-at-binding IS that expiry arriving. `[M]` two genuine per-leg
  bindings today either lose the leg information or red a well-posed sum;
  the resolution law must become BY LEG. No charter line reserved this
  until now.

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
CS3 → CS1 → CS1.5 → CS2 → CS4 (CS1.5 = the Medium un-weld, inserted 2026-08-20 by
the Q11 ruling in the CS1 design record).** The user asked whether the cone carve is independent of the
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

**⚠ SPLIT + PULLED FORWARD 2026-08-20; CHARTERED post-selection the same day —
the CS4a/CS4b/CS4c subsections below are the operative charter (§0 ruling 3's
⚠⚠/⚠⚠b carry the order).** The unified sketch at the end of this section
predates the split; where they conflict, the subsections govern. ⛔ The
sketch's done-when "all 12 P1 strict-xfails deleted" is CORRECTED: `[M]` the
tree holds **16** strict-xfail rows (5 R1 + 3 R2 + 8 R6); 4 R1-VALUE rows were
remedied at CS1-3b (recorded in-file, `test_monomorphic_leaves.py:668-677`)
and R6×8 was never P1's — the campaign total reads **"the remaining 8"**
(R1×5 + R2×3), apportioned CS4a → CS4c below.

### CS4a — the kernel core (CHARTERED 2026-08-20; ⏹ EXECUTED 2026-08-21 — ledger at the end of this subsection)

**Goal.** The interaction kernels exist as representation-free data and the
four energy-family operators are BOUND at construction — space carried and
validated, the homogeneous problem posed on Energy ⊗ point minted from the
mixture — while NO apply arm is deleted (the R-A census forbids space-keyed
selection until CS4c's feeding normalization).

**Selection record**: `scratch/cs4a_synthesis.md` (winner-plus-grafts + the
refutations R-A…R-F, each spot-re-verified); assemblies/attacks
`scratch/cs4a_assembly_*.md` / `scratch/cs4a_attack_*.md`; operator census
`scratch/cs4a_operator_facts.md`; objectives
`.claude/plans/kernel_and_medium_objectives.md` (per its honesty rows, O1
completes at CS4b, O4 at CS4c, O5 at CS2 entry).

**Rulings (user, 2026-08-20, forks F1–F4)**:
1. **F1** — fields are a TWO-STEP: CS4a lands operator→kernel; CS4b reworks
   the field mesh-requirement → SPACE ("fields need mesh to inquire about
   space; fields might be constructed from space alone" — `CrossSectionField`
   under the same criterion). The fabricated path retires at CS4b, where it
   becomes structurally unneeded — NOT at CS4a.
2. **F2** — phases subdivide freely; C's mandatory-space flip lands in one of
   them (CS4c, with its `[M]` 131/43-site test migration).
3. **F3** — the La13511/registry re-expression is POST-CAMPAIGN: ⭐
   verification cases CONFORM to the elegant production design; they are
   never a design driver. The Medium's layered arm ships witness-only,
   knowingly (the C10 tension accepted).
4. **F4** — roster = `ScatteringKernel` + `N2NKernel` + **`FissionKernel`**:
   the user's stress test overrides the today-only type-vs-property reading —
   a moment-based solver gives fission a frame-bound second realization, and
   Campaign 2's `GeneralizedEigenPencil`/resolvent (the α-eigenvalue
   resolvent) needs F rebindable as a first-class datum
   ([[feedback-defer-only-when-architecture-vague]]: concrete architecture +
   imminent consumers ⟹ build). Collision stays a property (ruling 2
   verbatim), the binding arity leaving a later coefficient-kernel admissible
   without re-carve. Re-examine the roster at Campaign 2's chartering.
   ⭐ **F4 addendum — the dyad's adjoint boundary (user question, 2026-08-20,
   design input for CS2/CS4c's angular binding).** The dyad fully determines
   the adjoint exactly where it IS the whole operator: energy-only spaces (+
   the spatially-diagonal extension), where counting measure makes metric
   adjoint ≡ transpose ≡ the χ↔νΣf factor swap BY THEOREM (the verification
   plan's own blindness measurement is this theorem read as a gate). On
   ANGULAR/MOMENT spaces the dyad UNDERDETERMINES the adjoint: F decomposes
   as (isotropic embedding E) ∘ (χ broadcast) ∘ (νΣf contraction) ∘ (angular
   retraction R), and the metric adjoint of the w-weighted retraction is the
   CONSTANT embedding — the Euclidean transpose is the w-WEIGHTED embedding,
   differing by a G-conjugation (the measured G6.3/87 % failure class) —
   with the Σw = 4π placement owned ONCE by the binding or the adjoint
   silently gains/loses it (failure-mode #3). ⟹ the angular/moment binding
   MUST realize F as the composite with R/E supplied by the SPACE's measure
   (or the frame: the ℓ=0 projection/injection pair is adjoint iff the frame
   is TIGHT — which is exactly why the sketch's tightness gate
   `bind(K)† = bind(K†) ⟺ tight` exists), never as a monolithic apply with
   a hand-coded adjoint. `.H` then composes through the existing operator
   algebra and the factor swap is a theorem of the factors. Caveat: metric
   adjoints need nondegenerate per-factor metrics — bulk angular weights are
   positive (safe); degeneracy lives on trace metrics, which F never
   touches.

**Steps**:

- **K0** — the per-row xfail marker split (`pytest.param`, the `_G13_ROWS`
  pattern at `test_monomorphic_leaves.py:741-747`; the function-level
  `_R1_XFAIL`/`_R2_XFAIL` marks cannot flip partially — `[M]` probed
  `2 failed [XPASS(strict)]`).
- **K1** — mint `orpheus/transport/kernels.py`: `ScatteringKernel` (the
  per-material Legendre stack {Σ_sℓ}, `p0`/`truncated(L)`/`ng` — the ℓ-index
  IS the Funk–Hecke eigenbasis index, the ground for CS2's frame-at-binding
  and the twin heal), `N2NKernel` (Σ_2n, multiplicity 2 — the loss-side
  channel ruling intact), `FissionKernel` (the (χ, νΣf) FACTORS, never the
  outer product; transpose = the factor swap). Views over `Mixture`;
  `MaterialXSField`'s dense caches absorb/delegate. The energy-arm rule
  (content-equal edges ⟹ `from_grid` else `synthetic`) HOISTS to one shared
  function; `bulk_space` re-points to it (no second spelling).
- **K2** — the bindings: the four energy constructors take (kernel/data,
  space[, assignment — by ARITY, absent on quotient bindings; ⚠ the arity
  arm has NO CS4a witness — its first meshed binding is CS4b/CS4c's, stated
  here so the gap reads as scoped, not missed]); the ng-conformity typed
  refusal keys on **presence + the ENERGY AXIS where the space carries axes**
  (⛔ verification-plan finding 2 CORRECTED this clause 2026-08-20: keying on
  "the space's shape" is UNRUNNABLE — `[M]` `space.shape[0]` is NOT ng on
  the SN composite (flat `(64,)`; interior ordinate-first `(4,2,6)`), a
  plugin probe destroyed 250 of 845 rows; the axis-keyed guard is live on
  `[M]` 192 of 1022 constructions and the fraction goes in the guard's
  docstring per vv#28; universal strengthening arrives with CS2's axes —
  main-agent ruling on the plan's blocking Q1 = its option B);
  `solve_homogeneous_infinite` mints its space from the mixture
  (`of_axes(energy, point)`) instead of reading `carrier.bulk_space`; the
  two homogeneous `IntegratedReactionRate` sites re-pose to the space's
  inner product (⚠ this INVERTS the D5 sensitivity partition — volumes stop
  entering the homogeneous path, so the byte gate's sensitivity note is
  re-dated in the same step); F space-mandatory (`[M]` 4 production + 9 test
  sites in 5 files). Kernel views are READ-ONLY (`[M]` the shipped
  `sig_s_legendre` returns a writable cache alias — mutating it moves the
  loss matrix; the mint closes this). ⛔ SUPERSEDED-in-scope 2026-08-21
  (CS4a-R EE-4): the mint closed the hazard only for KERNEL consumers —
  then a set of size zero — while the producer alias stayed live for the
  ~20 carrier consumers (`[M]` re-probed: a caller mutation reached the
  loss matrix at review time). The review froze the producer's dense
  caches at build (`_build_dense_caches`), closing it tree-wide; the
  kernel copy remains the guarantee that survives carrier rework. The degenerate carrier STAYS as the
  coefficient/mat_xs source until CS4b (F1). D5 byte gate 8/8 = exit
  criterion. The `.. implements::` directive left pointing at a symbol the
  re-posed path stops calling re-points in the same step (plan finding h).

**Verification plan of record**: `scratch/cs4a_verification_plan.md`
(2026-08-20; 24 gates K0 3 / K1 9 / K2 12, 25 battery arms incl. 3 positive
controls; scope 1373 rows / 26.19 s per arm; the DriftWarning wall carries
`[M]` 2 named pre-existing reds — baseline recorded before anything is
added). Its three charter refutations are folded above and in the done-when.

**Done when**: kernels exist with the slice-consistency 0-ULP crosscheck
re-pointed at the kernel object; all four homogeneous operators bound to the
mixture-minted E⊗pt space (the grep obligation is TWO predicates — no
`mat_xs.mesh.bulk_space` read left in `homogeneous/solver.py` AND
`from_mesh`'s chain not resolving the carrier's — per the plan's (d)); the
adjoint check ships as the plan's (c) re-scope (⛔ the chartered
"counting-measure adjoint THEOREM gate" CANNOT RED — `[M]` 0.0 under defect
and fix, ≤2.2e-16 at 56 000× volume spread, `[G,Aᵀ]=0` by theorem — so the
red-capable PREMISE leg lands plus one corollary row labelled THEOREM
carrying the blindness table and naming D4b as the only loaded partner; the
docstring may NOT claim flip coverage); R1-F/R2-F flip via K0's split (`[M]`
R1-C is NOT free — it stays; the ledger reads 16 → 14 at CS4a's close); K0
adds the permanent strict-introspection gate (`[M]` no `xfail_strict` in
`pyproject.toml`, so a lost `strict=` is silent); NO apply arm deleted (a
grep/AST obligation + the plan's 12-cell behavioural matrix, not a
mutation); D5 8/8; pyright terminal-1; Sphinx `-W` clean.

**Fences**: no apply-arm deletion; no field-class edits beyond what K2's rate
re-pose strictly needs (fields are CS4b's); no S/angular changes; no Medium;
C8 (frames/CS2) intact; solver entries unchanged.


#### ⏹ CS4a EXECUTED 2026-08-21 — the ledger (K0 → K1 → K2a → K2b on `feature/cs1-energy-space`, pushed; BRANCH-HOLD stands)

| step | commit | content |
|---|---|---|
| K0 | ✅ LANDED `069e2caa` | per-row `_R1_ROWS`/`_R2_ROWS` split (`_G13_ROWS` shape) + the PERMANENT strict-introspection gate `test_ledger_xfail_marks_are_strict`. `[M]` G0.1 98 node-ids byte-identical; G0.2 16 XFAIL lines identical; battery M0.1/M0.2/M0.3 per prediction (M0.2's zeros recorded: G0.1/G0.2 both blind, G0.3 the only instrument) |
| K1 | ✅ LANDED `15bbf935` | `orpheus/transport/kernels.py` minted (ScatteringKernel / N2NKernel / FissionKernel, frozen + read-only views over ONE `Mixture`; χ law through `enforce_emission_spectrum`; truncation refuses beyond-order) + `EnergyAxis.from_materials` = the ONE energy-arm rule, `bulk_space` re-pointed. G1.1–G1.9 (35 rows). Battery 9 arms; 2 test-file defects the battery itself caught were repaired in-flight (the symmetric control row's carrier made all-symmetric; G1.5's leg order made M1.6/M1.7 distinguishable) |
| K2a | ✅ LANDED `9f1d4190` | `_pose_space(mix)` mints Energy ⊗ point; `_assemble_loss_operator(mat_xs, space)`; rates → `space.inner_product` (`[M]` value no-op: D5 8/8, frozen rates 0 ULP); the shared `_energy_conformity` guard wired into all four constructors; G2.1–G2.10 (minus G2.11); the three §7.4 doc edits incl. the `.. implements:: normalisation :by:` re-point |
| K2b | ✅ LANDED `49b29391` | F space-MANDATORY (field + both annotations + `from_solver_data`); the 9-test-site §6b unit; D11 DELETED; R1-F/R2-F markers deleted — **ledger 16 → 14** (R1: L,C,S,B; R2: C,S; R6: 8); G2.11. Q3 ruled: hand-built survey space (no carrier exists behind the synthetic mat_xs) |

**Findings that SUPERSEDE plan/charter claims (edited here per plan-authoring §3):**

- ⛔ REFUTED 2026-08-21 (by the K2a battery, M2.2 re-run against the landed
  gates): *"the space swap has NO runtime witness; §7.1's grep is the ONLY
  evidence"* — TRUE when measured (the plan's tree had no G2.4/G2.5), FALSE
  after K2a: `[M]` re-wiring the solve to `mat_xs.mesh.bulk_space` reds
  **G2.4 + G2.5** (each wiggles one side of the measure identity, so they
  distinguish the mint from the carrier read). The greps remain obligations;
  they are no longer the only evidence.
- `[M]` the D5 sensitivity INVERSION (F3) measured on BOTH sides: pre-K2
  `volumes ×2` moved `flux 397.946 → 198.973` and the space weight was
  bit-identically inert; post-K2 `volumes ×2` moves NOTHING on the
  homogeneous path (G2.4 + D5 green; 24 off-path reds are the wider tree's
  legitimate volume-measuring tests) and the space weight ×2 reds **D5 8/8**
  + G2.5/G2.6/G2.1/D1 (19 rows).
- `[M]` M2.10c (the charter's original shape-keyed ng guard, run as a
  mutation): **275 rows destroyed** (207 failed + 68 errors) — the
  executable confirmation of verification-plan F2 (250/845 at plan time).
- `[M]` M2.13 (iso pair space-mandatory): **84 reds** — the F2/R-C scoping's
  measurement; the iso flip belongs to CS4c with its migration batch.
- `[M]` M2.12: the registry-keyset gate is structurally BLIND to the iso
  pair's `isinstance` arms (predicted; the behavioural matrix G2.8 is the
  instrument that catches it).
- ⚠ M2.6's value-motion catcher is **D5**, not G2.3 — G2.3 reads the pairing
  directly (its own spec), so solver-side contraction errors are caught by
  the byte gate's 6 non-1g rows.
- The campaign total *"the remaining 8"* now reads **"the remaining 6"**
  (R1-C/R2-C/R1-S/R2-S → CS4c; R1-L/R1-B → CS2). `[M]` R1-C confirmed NOT
  free (F5).

**Deferred consciously:** the arity arm ships with NO CS4a witness (§2(h).6 —
its first meshed binding is CS4b/CS4c's); `FissionKernel` has one gate and no
production consumer until CS4c rebinds F (Q2 — stated in its docstring);
IntegratedReactionRate survives for its meshed diffusion/SN consumers.

**The wide run (`tests/sn` whole, off-critical-path close-out):** `[M]`
2026-08-21 @ `49b29391`: **1 failed / 3275 passed / 3 skipped / 51 xfailed
in 2h06**. The one red
(`test_phase_e_trajectory_resolvent_flux_shape_crosscheck[cyl_2g_3reg_folded_4x8_dd_n40]`,
0.1268 vs the 0.12 gate, group 0 only) is **NOT CS4a's**: `[M]` deterministic
(bit-identical in isolation) AND fails at `32d2f548` (pre-CS4a, worktree
run); the §4.4 diff audit concurs (no SN-path value line). Tracked as
**#397** (leading lead: the 6.3 ω-fold re-baselined this snapshot's k_eff
while the MR reference is fold-independent; the 0.12 bound is pre-fold).
The 1-ULP DriftWarning wall's two pre-existing cyl reds are **#396**.

### CS4a-R — the kernel/operator machinery survives a clear-context adversarial review (user-chartered 2026-08-21; ⏹ EXECUTED 2026-08-21 — ledger at the end of this subsection)

**Goal.** The machinery CS4a landed — and the operator layer it now sits in —
holds up under a reviewer who did not build it: every Phase-1 attack is
either repaired, filed, or withdrawn in Phase 2, and whatever survives
RESHAPES the CS4b/CS4c charters before CS4b spends a line. (Precedent: §4-R,
where the same round on CS3 refuted its completeness claim and strengthened
5 gates — the review earns its session.)

**Why clear context is the method, not a convenience.** The builder's context
defends the build (it re-reads its own rationale as evidence). The reviewer
reads the TREE first — code, then gates, then docs — and only THEN the
charter/ledger (§5 CS4a's ⏹ block) to learn which oddities were *ruled*
rather than accidental. Same discipline as CS3-R, whose census caught what
the author's own grep missed (case-blind filter; the two-filter
countermeasure is now standing: any completeness claim gets a second,
independently-vocabularied filter).

**Scope — the surfaces (all landed `069e2caa`…`49b29391`, pushed):**

* production: `orpheus/transport/kernels.py` (the three kernels);
  `orpheus/transport/operators/_energy_conformity.py` (the shared guard);
  the four energy operators' binding surfaces
  (`fission.py` — space now MANDATORY; `multiplication_operator.py`,
  `isotropic_scattering.py` — `__post_init__` guards);
  `orpheus/numerics/axis.py` (`EnergyAxis.from_materials`) +
  `orpheus/transport/mesh/material_mesh.py` (`bulk_space` re-point);
  `orpheus/homogeneous/solver.py` (`_pose_space`, the re-posed rates).
* the WIDER operator machinery these now live in (`numerics/operator.py`'s
  algebra, the adjoint two-layer architecture, the apply-arm dispatch the
  R-A census froze) — Phase 1 may attack the layer, not only the diff;
  reshaping is fair game, including the CS4b/CS4c charters themselves.
* gates: `tests/transport/test_kernels.py` (50 rows),
  `tests/homogeneous/test_operator_spaces.py` (the D-suite + G2.*),
  `tests/sn/architecture/test_monomorphic_leaves.py` (the 14-row ledger +
  the strict gate); the two doc pages
  (`docs/theory/foundations/{spaces,infinite_medium}.rst`).

**Method.** Two phases, strictly ordered ([[feedback-adversarial-phase-
before-balance]]): Phase 1 harsh — *how would I break this / how would I
make it 100× better* — with "it works"/"it landed yesterday" inadmissible as
defences; Phase 2 re-evaluates and WITHDRAWS what does not survive (a
"don't touch" verdict belongs there, as a retracted attack). Dispatch
freely (elegance-enforcer, qa, cross-domain-attacker; an independent census
agent for any completeness claim); `method-implementer` stays out
(surgical posture). Verify any claim against the tree before adopting it
— the ⏹ ledger's `[M]` rows carry their configurations.

**Known deferrals are NOT findings** (each is chartered, with a home —
re-flagging them costs the round its signal):
F10's IsoS/IsoN2N untyped ScalarFlux→ndarray fall-through (recorded in
G2.8, repaired at CS4c); the arity arm's missing witness (§2(h).6 —
CS4b/CS4c's first meshed binding); `FissionKernel`'s one-gate/no-consumer
status (Q2 — stated in its docstring, consumer = CS4c's rebind);
`IntegratedReactionRate` surviving for meshed consumers (CS4b moves the
mesh dependency); the iso pair + C keeping Optional space (CS4c, with the
131/43-site migration); the ng guard's 4-of-13 reach (vv#28 — CS2's axes
strengthen it); the ledger's remaining 6 xfails (apportioned CS4c/CS2).
Open issues already tracking review-adjacent state: #396 (DriftWarning
wall reds), #397 (the MR flux-shape red — pre-existing, measured), #393,
#394, #395.

**Fences that still bind any in-review fix:** C8 (no S/frames reach from
the kernel module); R-A (NO apply-arm deletion before CS4c's feeding
census); BRANCH-HOLD (everything lands on `feature/cs1-energy-space`; no
merge without the user's go); the canonical runner + copy-aside battery
discipline for any gate touched.

**Deliverables / done-when:** a findings file
(`scratch/cs4a_r_findings.md`) with every finding carrying its Phase-2
verdict (CONFIRMED-fixed-inline / CONFIRMED-filed / WITHDRAWN, each with
the structural reason — a refuted attack is first-class output); inline
fixes committed under the standing verification discipline (battery scope
green, D5 8/8, pyright terminal-1, Sphinx `-W` if docs move); any CS4b/
CS4c charter amendments edited into their subsections IN PLACE with
attribution; a §-R-style ledger appended here on completion. Done when
Phase 2 has adjudicated every Phase-1 attack and the plan's next section
(CS4b) reflects whatever survived.

#### ⏹ CS4a-R EXECUTED 2026-08-21 — the review ledger

**Reviewers**: main agent (adjudication + own tree read) + elegance-enforcer +
qa + cross-domain-attacker + an independent census agent (two-filter
discipline). **Findings file**: `scratch/cs4a_r_findings.md` — every finding
carries its Phase-2 verdict; probes `scratch/cs4a_r_probe_*.py` +
`scratch/cs4a_r_census_*.py`.

| commit | content |
|---|---|
| ✅ LANDED `c7f9fa8d` | production repairs: XD-6 intensive condensed XS (⟨Σ,φ⟩/⟨1,φ⟩ — `[M]` shipped spelling scaled with the point weight, D5-bit-identical fix); EE-4 producer caches frozen (`[M]` the alias was live; the ledger's "mint closes this" superseded in place); EE-3 revised-not-swapped (`[M]` `CrossSectionField.ng` is a mesh read-through — the reviewer's fix was refuted at the tree); 7 docstring truth repairs; EE-5 api-page re-point + `.. implements::` re-point |
| ✅ LANDED `d61e097b` | gate strengthening: QA-F1 per-site guard witnesses (`[M]` 1 of 4 sites had one; teeth re-proved); QA-F4 introspected strict-gate population (`[M]` the `_R7_ROWS` evader now reds); QA-F5 both-properties R1; QA-F2 licence re-scopes + the hand-authored G1.4b convention pin; G1.5 frozen-producer upgrade; G2.5 rate leg → invariance witness; QA-F3/F6/F7/F8/F9/F10/F11/F12 |
| (this commit) | charter amendments in place (CS1.5′ ×1, CS2 ×4, CS4b ×3, CS4c ×4); the vv-principles #17 hoisted-guard sharpening (qa's, with its `[M]` table); fleet agent memories |

**Scoreboard** (details in the findings file): ~34 Phase-1 attacks across 5
reviewers → 19 CONFIRMED-fixed-inline, 12 CONFIRMED-charter-amendment (all
edited into their subsections above), 1 revised-at-fix-time (EE-3 — the
review process catching its own reviewer), 13 WITHDRAWN with structural
reasons (first-class output: MA-1/XD-4a heterogeneous-ng — user + the
carrier's eager unanimity gate; XD-5 two-denominators — correct layering;
qa's W1–W7 incl. three of its own probe-refuted attacks; elegance's 6 incl.
the ClassVar trap its own memory asserted, now corrected).

**The round's highest-value findings**: XD-1 (`[M]` the CS4c tightness gate
as originally chartered paired the adjoint law with the tightness
hypothesis — the negative control would have PASSED; re-specified as
Galerkin leg + multiplicativity leg + ℓ=0 blindness record); QA-F1 (the
hoisted guard's witness table — now a vv-principles #17 sharpening); XD-6
(a physics-semantics convention that had been recorded as an observation,
now decided: condensed XS are intensive); EE-4 (a founding-purpose closure
claim that was true and vacuous).

**Verification at close**: battery scope **1436 / 1 / 15** (+2 = the new
gates); whole `tests/sn/operators` **1230 / 1 / 5** baseline-exact; D5
**8/8** through both commits; pyright ratchet `{"transport": 1}`; Sphinx
`-W` clean; `dead_references` 0. BRANCH-HOLD stands; fences C8 / R-A
respected (no apply arm touched; no S/frame reach from the kernel module).

### CS4b — fields are space elements (goal chartered; grounding at its own session)

**Goal.** A field is constructible from its SPACE alone — the mesh
requirement retires across the field layer ("fields need mesh to inquire
about space", user F1; `[M]` sizing inputs: `BulkField.mesh` mandatory on all
10 leaves; ≥11 composite-arm output-mint sites; `CrossSectionField` probed
unconstructible without a mesh) — and THEN the fabricated path retires:
`from_materials` + the fake axis + the invented node, `MaterialXSField`'s
meshless admission (the mesh→assignment un-weld), the `mesh is None` sentinel
honesty, SN's typed promotion refusal (the bare assert `-O` strips today),
the `.areas` wrong-message, the `infinite_medium.rst:1115` xref +
`dead_references`, the `[M]` 16 test-site migration. **O1's
no-fabricated-data tell completes HERE**, re-measured by the byte gate.
⛔ **The "and THEN the fabricated path retires" half + the O1 tell are
SUPERSEDED (user F5 ruling, 2026-08-21)**: the SN path is the FIRST kernel
consumer — designing the binding against the homogeneous path is greedy/local
optimization. The fabricated path's last consumer (`homogeneous/solver.py:229`)
dissolves at CS4c's homogeneous CODA (after the SN rebind), and the O1 tell
completes THERE. CS4b's scope = the field/space machinery collapse
(constructor flip, axis-built cached spaces, retract/embed/restrict verbs,
space-equality identity, factory decomposition, the honest-refusal +
sentinel + `.areas` + doc repairs, the test migration). Design record:
`.claude/plans/cs4b_fields_design.md` (rounds 0-2; grounding census +
taxonomy census + kernel stress memo in scratch/). The CLASS-MERGE question
(one Flux class vs per-family leaves) is a NAMED decision point immediately
after CS4b lands — decided on static-typing + units-decomposition evidence,
not silently dropped (user F1 ruling).

**CS4a-R amendments (2026-08-21, review round — attribution in
`scratch/cs4a_r_findings.md`):**

* **EE-1 — the typed reaction-rate co-vector lands here.** The homogeneous
  rates read the space's pairing on bare arrays (recorded ruling: deliberate
  at CS4a; the tree's only raw two-ndarray `inner_product` production site).
  When fields become space elements, mint the space-bound reaction-rate
  co-vector (a functional, NOT a field method — `CrossSectionField
  .inner_product` requires same-class partners, which is exactly why a
  functional is the type that pairs across classes), re-point the three
  solver lines and the `.. implements:: normalisation` at it.
* **XD-10 — charter language precision:** the honest sentence is *"a field's
  SIZING is derivable from its space; its ROLE is not"* (`FunctionSpace`
  explicitly disclaims units/role, `space.py:150-156`). And
  `has_coordinate_cone` is `None` on axes-less spaces — retiring the mesh
  requirement does not by itself retire the unanswerable cone predicate;
  say which fields inherit `None` and why that is acceptable.
* **EE-5 — the grep obligation for ANY data-flow removal in this phase (and
  CS4c) is `orpheus/` + `docs/` + `tests/`, by CONCEPT as well as symbol.**
  CS4a's done-when spelled two `.py` predicates and
  `docs/api/homogeneous.rst:88` shipped present-tense-false for a day
  (found by the review; "the carrier's axis-built bulk_space" — greppable
  by symbol AND by the concept "threaded domain").

**⏹ CS4b DESIGN COMPLETE (rounds 0–3c, 2026-08-21/22) — execution NOT
started.** Surfaces, in pickup order: the design record
`.claude/plans/cs4b_fields_design.md` (⚠ its own header: later rounds
supersede earlier fork text IN PLACE — read the round sections first) → the
verification plan `scratch/cs4b_verification_plan.md` (7 steps S1–S7; the
22-arm identity table §11.1; gate specs with named witnesses; §16 done-when
predicates; §0 baselines) → the census memos (`cs4b_grounding_census.md`,
`cs4b_field_taxonomy_census.md`, `cs4b_kernel_as_frame_stress.md`, all in
scratch/). Ruled at design: effort is NEVER a decision criterion (round 1);
the machinery collapse is CS4b's committed scope and the CLASS-MERGE is a
named decision point after; F2 = space CONTENT equality (`[M]` every space is
law-blind, correctly — laws are operator data); SN-first kernel consumer
(F5 reversal); Q3 no-gate (`[M]` no concrete site); Q5 = **complete the frame
square** — the frame gains `dual_reconstruction`, the field-layer /Σw
spellings are its hand-rolled ℓ=0 shadows, the factor is G₀⁻¹ DERIVED; Q6
refusal ratified with the exact-modal-cone extension filed as **#400** (and
the moment-tail bug as **#399**). Step sequence S1 carrier mints → S2 the
METRIC step (`[M]` bulk `l2` moves 41 % — CS4b owns the
`test_si_diagnostic_trajectory` re-derivation, ρ≈c anchor licence) → S3
consumers + the 22-arm re-key → S4 the mesh field deletes → S5 sugar
migration → S6 verbs + #399 → S7 repairs/EE-1/docs. The full ruling table
(Q1–Q8, epistemic markers) = the design record §Round 3c. Execution posture:
surgical (main agent writes, user steers); it OPENS by writing the L17
convention crosswalk from the two surfaces, then S1. Session commits:
`466e6756`/`57270819`/`455d5f9a`/`34df88cb`/`3f358998`/`44a4a4b8` + the fold.

**⏸ CS4b EXECUTION — COMPACTION POINT (2026-08-22): S1–S3 ✅ LANDED, S4
opens next.** Commit table (all pushed; `git merge-base --is-ancestor`
each before trusting this row):

| unit | commits |
|---|---|
| L17 crosswalk | `2b4bfed3` (+ resolutions fold `e8ecdd80`) |
| S1 carrier mints + G1 gates + 6-arm battery | `4069155b` |
| Q7 ruled **option 3** (user): `_synthetic_for_tests` retired, real-carrier fixtures | `be608c01` |
| S2a the metric step (space-source flip; `_SPACE_NAME` retires; G-M1..3 re-derived under the ρ≈c anchor, move `[M]` 1.117e-3) | `07e0fe77` |
| S2b Q2 (composite interior IS the mint; the dense oracle moves in-test) | `8a205cbf` |
| S3 (four commits): partner family + digests → operator arms via `space_on` → 8 manufactured witnesses + permission legs → F4/read-throughs + the 22-arm battery | `9138b3c3` `8b3d0cd4` `5da03d05` `a82d31e4` |

**Rulings landed IN execution (supersede the pre-execution text where they
touch it):** ⭐ the B5 refusal machinery SUPERSEDED (user walk-through — the
corners are unspellable by CURRENCY; the int is a retiring surface; the
scheme-binding doctrine recorded at crosswalk B5); Q7 = option 3;
`moment_axis` aligned to the `is_multi_moment`/`NotImplementedError` idiom;
the trace/ray spaces carry CONTENT-DIGEST names (the F2 doctrine's
face-family leg — one mint site each); the FaceField layout sub-arm RETIRED
tautological-under-the-digest (`[M]` battery-measured blind); the S3/S5
consumer boundary refined (factory-ARG `.mesh` reads migrate at S5 —
crosswalk "S3 execution refinements"). ⚠ #401: the strict DriftWarning wall
carries 2 INHERITED cylinder reds (`[M]` bit-identical at the session base;
CS4b contributes zero, bisected three ways) — adjudicate before the next
step that must read the wall.

**✅ S4 LANDED 2026-08-22 — S4a `554ff10b` (F4) + S4 `1333135e` (F5 + Q1
+ the deletion, one §6b unit).** The rulings below were executed the same
session they were made; the landing details are in the two commit
messages (the honest record: the space-route consumer re-points, the
harmonic factors[1] rebuilds, the widened-member S6 defers, the ~116
test re-keys, the retired `as_per_ordinate`). Landing-time discoveries
that later steps must know:

- ⭐ ~~**The frame square's structural asymmetry**~~ ⛔ MISCHARACTERIZED
  AS STATED (user diagnosis, 2026-08-22; corrected by the S4-AMENDMENT
  below): the `space=` need on `reconstruct` was an UNBOUND OPERATOR's
  missing codomain, not a structural fact — the frame was never bound
  to its two field spaces at construction. The honest residual fact is
  the CONSTRUCTOR-INPUT direction: the moment space is derivable from
  the angular space + `L`, never the reverse (the angular target
  carries the quadrature axis and, widened, the scheme's mass-bearing
  moment axis) — so the constructor takes the ANGULAR space and
  derives the moment space. The production caller passes its posed
  composite's `interior_space` — at CONSTRUCTION, not per verb call.
- ⚠ **`Axis.weights is None` IS the all-ones measure** (canonicalized
  at construction — `[M]` GL2's w=[1,1] reads back None). Every weight
  reader carries the branch.
- **pyright ratchet 1 → 0, re-baselined `{}`** (§8.4) — the #226
  terminal state; S4's honest typing also retired the pre-existing
  transport error.
- The widened (`spatial_moments > 1`) harmonic members
  (`scalar_flux`/`truncate` — self-derivation blocked by the scheme
  axis) are LOUD defers naming S6/G6.7; `scalar_flux` takes `space=`
  for the widened caller (the windowed moment arm passes its composite
  interior's marginal).
- Suites at the landing: transport 536 / sn core 1690 / architecture
  114 / solve 202 / fast set 2666 (D5 8/8) / cp+mc 180 (slow
  deselected). `sn/verification`, `derivations`, and the slow marks
  belong to the ONE ≥90-min pre-merge gate (branch-hold ruling).
- S7-owed corpus staleness (chartered, same branch):
  `field_algebra.rst`'s "mesh identity" layer rows;
  `operator_algebra.rst`'s `_phase_space_shape` xref; the frame-square
  + naming-rule pages. #402 filed (the LD-certificate un-skip rider).

**⛔ S4-AMENDMENT CHARTERED FIRST (user ruling, 2026-08-22, post-S4
review) — an operator is not an operator without its two spaces.**
The user diagnosed the S4 `space=` parameter on `reconstruct` as the
symptom, and ruled the repair at TWO levels:

1. **Bind the typed frame (the local repair).** `HarmonicFrame` binds
   its two full field spaces at CONSTRUCTION: the constructor takes the
   ANGULAR field space (the one production site,
   `ScatteringOperator.frame`, holds the posed composite's
   `interior_space` — iterate-width, LD included) and derives + stores
   the moment space once (`_moment_space_for` moves from
   apply-time-per-operand to binding-time; the derivation direction is
   a CONSTRUCTOR-input fact — moment = f(angular, L), never the
   reverse). The verbs lose `space=`; admission = content equality
   against the bound domain; outputs ride the bound codomain. ⛔ The
   "frame square's structural asymmetry" recorded at the S4 landing is
   MISCHARACTERIZED as stated (it read an unbound operator's missing
   codomain as a structural fact) — correct the S4 bullet above, the
   topic memory, and the harmonic_frame docstrings IN this amendment.
   `[M]` the numerics-level faces (`frame.analysis`/`.reconstruction`)
   are already fully bound with tested `.H` (test_frame.py's
   adjoint-for-free rows) — the defect is ONLY the typed lift
   `M ⊗ I_cells`, which was never minted as a bound operator.
2. **⭐⭐ THE STRONGER FORM (the user's root-cause ruling, verbatim in
   substance): the base class/protocol must DEMAND domain and
   codomain.** "The reason this happened in general is because the base
   class or protocol of operator is not demanding a domain and
   co-domain. If we define this at the base, and then go fixing
   everything that uses the operator machinery, the problem will
   disappear — Analysis, Reconstruct, Synthesis, Scattering, Fission,
   etc. everything will need to have its domain and co-domain defined."
   `[M]` the base's own docstring names the hazard it grants:
   `operator.py:658-676` defaults BOTH to `None` and "when either
   operand of a composition has None ... the composability check is
   SKIPPED — preserving the legacy duck-typed behaviour". Blast
   surface, measured 2026-08-22 at `1333135e`: **77 LinearOperator
   subclasses** outside operator.py; **44 `def domain` overrides**
   (⟹ ~33 inherit the None default); 5 overrides return a literal
   None; plus Optional-returners (`ScatteringOperator.space:
   FunctionSpace | None` et al.). ⚠ RECONCILE WITH CS4c's charter
   before designing: CS4c already owns "C space-mandatory (the
   131/43-site migration)", the iso pair's Optional space, and the
   binding-base `BoundOperator(datum, space)` — the amendment must
   decide what pulls FORWARD (the base-level DEMAND + the .H/R2
   hazard closure) vs what stays CS4c's migration (per §6b, the
   call-site sets decide, not the tidy narrative). The R2 hazard the
   mandate closes is already catalogued: `test_monomorphic_leaves` —
   "an anonymous leaf's `.H` is a bare Euclidean transpose wearing
   the Hilbert adjoint's name".

**⏹ S4-AMENDMENT EXECUTED 2026-08-22 (design ruled + all four steps
landed the same session: `33950d81` → `6e04a749` → `6fc247fb` →
`aa508d3f`; ledger below; [M] full fast set 7790/0 no-cache, pyright
0). The ruled order continues at S5.**

*Corrected census (supersedes the charter's "77 subclasses / ~33
inherit None" — that was a text-grep over-count, and a first AST pass
mis-attributed inheritance by DFS order instead of reporting all
sources).* `[M]` AST census over `orpheus/`
(scratchpad `census2.py`): **56 transitive LinearOperator subclasses;
the whole inverse family (`InverseOperator`, `MatrixInverseOperator`,
`GreenOperator`, `SweepOperator`, `CoupledSubstitutionOperator`)
already derives via `InverseWrapMixin`'s domain↔codomain swap;
composites derive; `_CarrierMatvecOperator` is scipy's class (out of
scope). True silent set: 4 concrete leaves + 2 ABCs** —
`IdentityOperator`, `DiagonalOperator`, `RankOneOperator`,
`DSACorrection`, plus `AnalysisOperator`/`ReconstructionOperator`
(whose docstrings already state the V/W pair "lives in the runtime
domain/codomain"). `FissionOperator` is the in-tree precedent (space
MANDATORY since CS4a K2, non-Optional properties).

**The six rulings:**

1. **Base demand = "abstract now, narrow later."** The root's
   None-returning default bodies are DELETED (`@property
   @abstractmethod`); every class must answer — with spaces, by
   derivation, by the pointwise law, or an explicit documented
   Optional override naming its owning campaign (S/C/iso → CS4c,
   L/B → CS2). The annotation stays `FunctionSpace | None`; the
   terminal narrowing + composability-skip retirement land when those
   campaigns finish. No CS4c/CS2 call-site set is pulled.
2. **`.adjoint()`/`.H` refuses unbound** — an adjointable-but-unbound
   leaf's `.H` today silently returns the Euclidean transpose wearing
   the Hilbert adjoint's name (the R2 hazard). Non-pointwise operators
   must declare both spaces to take `.H`.
3. **The natural family = `PointwiseOperator`** (ruled name): the
   space-polymorphic stratum of the multiplier algebra the tree
   already documents (`MultiplicationOperator`'s own law-suite
   `M[1]=I, M[0]=0`). Law: endomorphic at the operand's space, no
   stored pair (`domain`/`codomain` None BY LAW, discriminated by
   type); **`.H = self`** — real multipliers commute with every
   diagonal metric, so the adjoint is metric-free. Members: `Identity`
   (×1), `Diagonal` (×f), the endo `Zero` (×0). Boundary: Permutation
   is endomorphic but basis-coupling — excluded, its `.H` genuinely
   needs the metric. Bound multipliers (`MultiplicationOperator`,
   `InverseMetricOperator`) stay bound classes. `[M]` the user's
   albedo reading confirmed from physics: law = (scalar multiplier
   α∈[0,1]) ∘ (geometric link Γ₊→Γ₋); 0 = Zero, 1 = Identity,
   between = Scaled; the natural piece is the multiplier, the bound
   piece the link. `MultiplicationOperator`-as-collision is NOT a
   misnomer — #261 ruled mechanism-at-class/role-at-binding
   (`collision = MultiplicationOperator(...)` at all 3 sites).
4. **The hetero-Zero weld: FULL un-weld now.** `[M]` three layers:
   (a) the fission (B,B) hooked Zero exists because `CoupledOperator`
   refuses all-None columns ("drop the system") — impossible for a
   pencil member sharing the domain — and the refusal is load-bearing
   (`apply_transpose`'s empty `reduce`); (b) the missing operator is
   the RESTRICTION: `F_posed = [[F],[E]]-stack ∘ restrict_A`
   (rectangular stacks are already blessed; 1-member coupling "is the
   legitimate degenerate"); (c) vacuum is NOT the same weld — B3.2's
   zero morphism stands (the α→0 member of the law family, realized
   structurally, link never built) — but its `_zero_rows` closures
   duplicate what its already-bound codomain knows. ⟹ the
   `codomain_zero`/`transpose_zero` hooks lose both production
   consumers and RETIRE; `ZeroOperator` splits into the natural endo
   zero (stateless ×0 echo, joins `PointwiseOperator`) and
   **`ZeroMorphism`** (born-bound: both spaces REQUIRED, apply/
   transpose mint zeros from the bound pair, `.H` = the swapped
   `ZeroMorphism` — metric-free trivially).
5. **Restrictions: one concept, mint only the member projection.**
   All the tree's restrictions are ONE split pair (embedding e,
   retraction r; r∘e = id; r = e† in the right metrics; P = e∘r the
   orthogonal projector), in two realizations: index-subset
   (gather/scatter — `TraceRestrictionOperator`, Γ± half-splits,
   block projections, member projection) and measure-weighted
   (constant-broadcast / w-average — the kernel's pair,
   `integrate_angular`, `from_isotropic`). Many wear cloaks
   (`.interior` reads, `y.systems[i]`, index arrays). This amendment
   mints ONLY **`SystemRestrictionOperator`** (coupled → 1-member
   coupling; adjoint = extension-by-zero via the space's wired
   `zeros()` seam; `TraceRestrictionOperator`'s sibling, carrying the
   split-pair law-suite as its intrinsic-property test). The family
   ABC waits for its first generic consumer (CS4c's verb) per the
   type-minting rule; the concept goes into the operator-algebra
   corpus at S7.
6. **Level 1 (the frame) as chartered**: `HarmonicFrame(basis,
   measure, *, angular_space)` derives + stores `moment_space` once
   (named `angular_space`/`moment_space` — absolute names; the verbs
   run opposite ways so domain/codomain are direction-relative);
   verbs lose `space=`, admission = content equality, outputs ride the
   bound pair; `ScatteringOperator.frame` demands the posed space (the
   windowed-arm refusal moves onto the property — the kernel path
   demands posedness too; production always poses, `[M]`
   solver.py:1403 + both derived-sibling ctors).

**Execution order (§6b — the hook retirement's call-site set is
{vacuum site, fission site} + the hook params, so the split lands only
after BOTH consumers move):**
- **A1** ✅ LANDED `6e04a749` (+ the S4 sweep-straggler batch
  `33950d81` its full-fast sweep surfaced) — level 1: the frame
  binding + scattering re-point + test migrations + TestBinding gates
  + the "structural asymmetry" record corrections. The per-axis
  spatial-moment read hoisted to the static
  `BulkField._spatial_moments_per_axis_of(space)` (one source for
  fields AND the binding-time derivation).
- **A2** ✅ LANDED `6fc247fb` — `SystemRestrictionOperator` (born
  bound; split-pair law suite ×8) + the fission posing re-spelled as
  stack ∘ restrict. ⭐ In-execution discovery: the extension-by-zero's
  member CLASS is load-bearing (`[M]` the primal mint put a
  flux-classed ray zero into the daggered chain — the cross-class gate
  refused), so `CoupledSpace` gained the materialization seam's
  COTANGENT twin (`dual_zeros=` / `dual_zeros()`, #276 A4 duality
  typing; the SN builder wires both at both sites via
  `_zero_full_field_dual`).
- **A3+A4** ✅ LANDED `aa508d3f` (fused in execution — the `.H` refusal needed the
  metric-free predicate, which needed the family): `PointwiseOperator`
  (Id ×1 / endo-Zero ×0 / Diagonal ×f; `.H = self`;
  `is_metric_free_adjoint` the recursive predicate, derived on
  Sum/Product/Scaled/TensorProduct/SumOfTP; the bound multiplier
  `MultiplicationOperator` DECLARES it True — bound-ness and
  metric-freeness are orthogonal); the Zero split (`ZeroOperator`
  natural stateless; `ZeroMorphism` born-bound, zeros from the
  declared pair + the payload-tail convention — `[M]` the trace seam's
  spaces are structural-axes-only in the `face_method_space` fixtures,
  full-shape on `SNMesh.angular_trace`); the vacuum site migrated
  (`_zero_rows` + the `codomain_zero`/`transpose_zero` hooks RETIRED);
  the root's None defaults DELETED (`@abstractmethod`, four legitimate
  answers documented at the property); `_AdjointOperator.__init__`
  refuses an unbound non-metric-free inner (R2 closed —
  `test_dropping_the_binding_BREAKS_the_weighted_law` upgraded from
  measuring the silent Euclidean degradation to pinning the refusal);
  `OperatorProduct`'s domain/codomain gained the pointwise-conforming
  arms; `DSACorrection` derives from its mesh; `RankOneOperator`
  declares the documented CS4c debt; 10 test doubles across 6 files
  declare the unbound state; the faceless-vacuum escape CLOSED
  (refusal, not silent realization).

THEN the ruled order continues: S5
(the sugar-FACTORY retirement — the 909-site `zeros_on`/`from_mesh`
spelling migration; the boundary re-sharpening above already moved the
consumer-frame reads into S4) → S6 (verbs + the frame square + #399 —
its G6.3/G6.8 adjoint gates PRESUPPOSE the bound frame) → S7 (repairs
+ EE-1 + docs).

**⛔ FRAME-FACTORY RE-CARVE CHARTERED (2026-08-23, post-amendment user
review) — the typed frame faces become BOUND OPERATORS minted by the
frame; the frame reverts to the shared (basis, measure) factory.**

**✅ RULED 2026-08-24 — the design session CONCLUDED; execution plan
installed.** All leans ratified after five pressure-test rounds + the
desktop master-report reconciliation, and the frame square MEASURED
(`[M]` Parseval probes: the stored moment metric is WRONG-SIDE, ratio
118.7 vs 1.000 with the inverse; `M.H = R/W` to 5.6e-17 with the right
one; the Parseval metric is the inverse of the frame's DISCRETE Gram —
frame-owned, not a basis constant).
**▶ RESUMES AT: execute `.claude/plans/frame_square_recarve.md` —
F-0 (metric truth) → F-1 (the mint).** That file is the execution
authority (rulings, probe numbers, §6b call-site sets, the CS4c /
Phase-S debt ledger); the charter text below stays as the DESIGN
RECORD.

**Goal (outcome).** Analysis/Reconstruction at the transport level are
first-class BOUND OPERATORS (domain/codomain = the two full field
spaces — the S4-amendment's demand landing on the right objects), and
the frame is the SHARED per-(quadrature, L) factory that mints them —
one projection table, many bindings.

**The user's diagnosis (confirmed by census).** A1 bound the field
spaces TO THE FRAME — a proxy: it bound the factory when the ruling
("Analysis and Reconstruct are operators") was about the products. The
frame is not an operator; it is an operator factory, and it is SHARED:
`[M]` four consumers beyond the scattering operator's typed verbs, none
wanting A1's field-space binding —
- `loss_representation/__init__.py:4821` — the windowed accumulation is
  principled-equivalent to "the post-sweep flat `FrameBase.analysis`
  projection" (the user's remembered moment projection — exact);
- `sn/operators/windowing.py:71` — `BulkAnalysisOperator(frame,
  sn_mesh)`: frame IN, bound operator OUT — the factory pattern already
  shipping; its docstring says "``frame`` must be the SCATTERING
  operator's own angular frame … the same object, single source" AND
  "mint the synthesis-side lift when a consumer arrives"; its
  `codomain` returns None with "no minted FunctionSpace yet … the
  honest 'not yet typed'" — a debt A1's `_derive_moment_space` can fill
  and A1's placement could not reach;
- `sn/acceleration/dsa.py:569` — reads the frame's ℓ=1 table row (the
  single source of the μ row);
- `sn/operators/loss_kernel_gauge.py:1171` — builds its own
  `GalerkinFrame(basis, measure)` for the Gram machinery.

**A1's three symptoms** (why the fused shape is wrong): (1) the frame
stopped being shareable — one consumer's field-space context welded to
the per-(quadrature, L) table object; (2) the typed faces are METHODS —
no `@`-composition, no `.H`, not operators under the base demand,
breaking the level symmetry (numerics `FrameBase` correctly mints
`_FrameAnalysis`/`_FrameReconstruction` as bound OPERATOR instances,
adjoint-for-free pinned in test_frame.py); (3) the kernel ndarray faces
(`conjugate` — no field spaces needed) became pose-demanding for no
reason.

**Proposed means (2026-08-23, MY assessment, NOT ruled — the user has
pressure tests pending):**
- `HarmonicFrame` reverts to (basis, measure) + carrier knowledge; its
  identity returns to the table's identity (frames equal ⟺ shared
  projection).
- The frame MINTS bound typed operators on request — proposed spelling
  `frame.analysis_on(angular_space)` / `frame.reconstruction_on(...)`
  (⚠ spelling + class names e.g. `HarmonicAnalysisOperator` are
  HYPOTHESES; `[M]` no such symbols exist at HEAD). Domain = the
  angular field space, codomain = the derived moment space (A1's
  `_derive_moment_space` + admission checks MOVE to the mint —
  derivation direction fact unchanged: moment = f(angular, L), never
  reverse). Role-polymorphic `apply` (the ScatteringOperator overload
  precedent; spaces are role-blind — S4c's "role transition = same
  space, new class"). They join the `AnalysisOperator`/
  `ReconstructionOperator` ABC family — algebra citizens WITH `.H`.
- `ScatteringOperator` mints its pair at binding and keeps the frame
  ONLY as an accessor — CS4c's mint-and-FORGET contract arriving early
  ("auxiliary structure minted at binding, retained only as an
  accessor; the frame accessor + a declared analysis-verb field"). ⚠
  Honest caveat: the hot aniso kernel is the frame-FUSED `conjugate`
  chain (the 0-ULP canary), so S genuinely keeps the frame accessor
  until CS4c re-expresses the kernel through the binding base. XD-1's
  "the binding must DECLARE which analysis verb it realizes" is
  satisfied structurally — the minted operator IS the declaration.
- `BulkAnalysisOperator` fills its codomain debt from the same mint
  (or is later re-expressed as minted-analysis-on-bulk ⊕ trace
  identity).

**The strongest argument**: S6's G6.3/G6.8 adjoint gates need the
typed lift's Hilbert adjoint (metric sandwich over the FIELD spaces).
Method-verbs have no route to `.H`; minted bound operators get it from
the exact machinery the amendment hardened (and `_AdjointOperator` now
refuses to fake).

**What carries over from A1** (this is a completion, NOT a revert):
the derivation logic, the content-equality admission, the
"structural asymmetry" record corrections, and the TestBinding gates
(re-keyed frame → minted operators). Sequencing: this is S6's OPENING
MOVE ("verbs + the frame square"); independent of S5's mechanical
migration — either order works; my lean (unruled): re-carve before S5
while A1's context is hot and before anything grows a dependency on
the frame-level binding.

The original proposal text, for the record (all three ruled above):

1. **F4's replacement — the composite becomes an element of
   `FullFieldSpace`** (the ruled elegant form): the slot gate compares per
   BLOCK (`space.interior_space is/== interior.space` + trace — R2: the
   composite's own `==` is name+shape, block-blind), and `Composite.mesh`
   retires. Open design bit: WHERE the composite's space comes from at
   construction — `TimedFullField.zeros` has the carrier
   (`mesh.full_field_space`); direct `FullField(interior=, boundary=)`
   builds do not. ⚠ LD nuance: the cached composite's interior is the
   width-widened product (`==`-comparable, not `is`, to a per-call widened
   field space) — the slot gate needs the `is`-fast-path + `==` form.
   The transitional `# type: ignore[attr-defined]` lines in
   `full_field.py` die here with the reads.
2. **F5's pose home** — "does this problem carry System B" is the
   OPERATOR/pose's property (the field's space cannot carry it; the user's
   scheme-binding doctrine is the lens: augmentation binds what the method
   needs). Candidates: the `evaluate_residual` signature gains the pose/
   carrier; or the guard moves into the coupled-assembly seam where the
   pose is decided. Its first-ever witness:
   `test_space_content_witnesses.py::TestF5ResidualPoseGuard`.
3. **Q1 — `MomentResidual`**: the S4 `_from_balance` flip
   (`cls(values=lhs.values − rhs.values, space=lhs.space)`) dissolves the
   2-arg blocker (`L` recoverable via
   `space.find_factor(SphericalHarmonicSpace).L`). Standing
   recommendation (round 3c, unvetoed): record-as-choice at
   `_from_balance`, MINT only with a consumer. The flip retires the two
   `type: ignore` in `numerics/field.py` — a ratchet DECREASE is itself a
   re-baseline obligation (§8.4).

**✅ ALL THREE RULED (user, 2026-08-22, S4 design proposal session):**

1. **F4 = DERIVED + digest names.** `Composite.space` is a cached property
   `FullFieldSpace.from_blocks(interior.space, boundary.space)`;
   `from_blocks` adopts the of_axes rule (name derived from member
   `(name, shape)` content, `full_field#<digest>`; the `name=` role param
   RETIRES — G2.3: role is class identity, so `"radial_characteristic"`
   leaves the SPACE name). The `__post_init__` cross-slot mesh gate
   retires with NO construction-time replacement (coherence is not
   spellable from block content without a carrier reference; the
   S3-landed seam arms are the law — the currency doctrine). `[M]` 51
   direct production build sites made stored-`space=` boilerplate;
   composite `==` becomes exactly as content-keyed as its members.
   Consequences: `test_full_field_space.py` name pin → prefix pin;
   `test_typed_residual_evaluation.py:228` prefix pin survives per R4;
   `test_radial_characteristic_field.py:267` re-keys to space objects;
   `test_mixed_mesh_composite_rejected`'s TWIN fixture flips to the F2
   positive leg. Named non-goal: `CoupledSpace` untouched.
2. **F5 = SYSTEM-BOUND.** `evaluate_residual(system: WithinGroupSystem,
   psi, q_ext)`; the pose guard becomes the loss grid's ARITY (typed
   refusal names the arity, never the mesh); `[M]` all 4 production
   certificate sites' `system.loss if coupled else _bare_loss_arm(system)`
   ternaries die and `_bare_loss_arm` folds into the seedless arm;
   `TestF5ResidualPoseGuard` re-keys to the arity refusal. 2 prod + 6
   test call sites migrate.
3. **Q1 = RECORD-AS-CHOICE.** The flip `cls(values=lhs.values −
   rhs.values, space=lhs.space)` is `[M]` exact for all six shipped
   residual leaves (uniformly `(values, space, mesh)` today). Refinement
   over the recorded design input: a future `MomentResidual` needs a
   small `from_balance` override threading the operand's stored
   `(L, spatial_moments)` fields — NOT the `find_factor` route (the
   fields exist on `MomentField`). The stale "lands on cls's OWN space"
   docstring paragraph (false since S2a — role-blind shared mints) is
   fixed in the same edit. Rider: the flip may also dissolve
   `_residual_is_expressible`'s LD skip — un-skipping is a behaviour
   change needing witnesses (vv #26), FILED as a follow-up issue, not
   ridden.

Then S4 proper: mesh retires from the 2 declaration roots; the §8.3
test-side `.mesh` reads (~116) re-point; `_phase_space_shape` collapses to
space-vs-values; `Field.zeros`'s `mesh=` extra retires. After S4: S5 (the
909-site sugar migration, NOW INCLUDING the factory-arg `.mesh` reads per
the refined boundary) → S6 (verbs + the frame square + #399) → S7
(repairs + EE-1 + docs). Standing: BRANCH-HOLD; C8/R-A fences.

⚠ **S4/S5 boundary re-sharpened at execution (2026-08-22, §6b):** the
"factory-arg `.mesh` reads migrate at S5" refinement and "S4 deletes a
field nothing reads" were in tension — `[M]` AST census at the S4b
commit: **21 consumer-side field-receiver reads survive in production**
(scattering 7, fission 6, multiplication 4, harmonic_frame 2, isotropic 1
— plus ~7 under field-vars `Q`/`moment`/`lhs`/`seed`/`chi`) and the
fields' own internal machinery still reads `self.mesh` (~20 sites:
`_bases` 11, `harmonic_moment_flux` 6, `angular_flux` 1,
`scalar_source_sink` 1). §6b rules: the deletion's unit of work is the
call-site set, so **the S4 deletion batch absorbs re-pointing every
surviving field-receiver read** (to the operator's own carrier member or
the space route — mechanical, in the old factory spelling where one
survives), and **S5 keeps only the sugar-FACTORY retirement** (the
909-site `zeros_on`/`from_mesh` spelling migration). The `self` reads on
OPERATOR classes (134 raw, dominated by `loss_representation`) are
operator-side declared members — the #226 F2 pattern — and survive by
design.

### CS4c — the dispatch collapse (goal chartered; after CS2 per the ruled order)

**Goal.** Apply-time overloading retires: the per-instance
feeding-normalization census first (the R-A `[M]` 6-of-12 fact; the method is
vv-principles #29), the #205 ndarray keep-ruling re-litigated at the call
sites that feed raw vectors (the k-eigenvalue path), arm deletion per
binding, C space-mandatory (the 131/43-site migration), S →
`ScatteringKernel` shell with the iso pair + `LegendreMomentScattering` as
ℓ=0/moment bindings of the one datum (the O7 twin heal; R1-S/R2-S flip). May
pull earlier than post-CS2 at its chartering if frame-independence holds (C8
reserves only the frame MINT and the L/B + R6 rows).

**⭐ SN-FIRST consumer ordering (user F5 ruling, 2026-08-21):** the SN path
drives the binding design — it exercises every axis, retract/embed/restrict,
and all kernels; the homogeneous path re-points LAST as the degenerate coda
(where `from_materials`' one production consumer dissolves and O1's tell
completes). Designing the binding against the homogeneous path was rejected
as greedy/local optimization.

**Binding-base design inputs (round-2 kernel re-engagement, 2026-08-21 —
full reconciliation in `cs4b_fields_design.md` §Round 2):** the common-part
abstraction is the binding base as a dataclass ABC (EE-6 lifted to
`BoundOperator(datum, space)`) — space admission, domain/codomain, and the
mint-and-FORGET contract (auxiliary structure minted at binding, retained
only as an accessor; the frame accessor + a declared analysis-verb field is
where XD-1's declaration lives). Per-channel recipes live on the bound
operator's constructor (operators import kernels — C8 direction), kernels
stay array-verb data (the Basis analogue). Datum KINDS stay three under the
one base: integral kernel (S/N2N/F), multiplier (C), differential-stencil
(L — shares the binding shape, is NOT a kernel; `[M]`
`IntegralKernelOperator` is the strict nonlocal-integral discriminator and
L/C fail it). The RESTRICTION verb joins retract/embed (trace measure = the
restricted bulk measure dV→dA × |Ω·n̂|w — already how `angular_trace` is
built; composite block projections and Γ± half-splits are the true
subselections; restriction† = extension-by-zero; the R∘G rewiring stays on
the boundary thread #367).

**The original unified sketch (predates the split; the subsections govern):**

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

**CS4a-R amendments (2026-08-21, review round — attribution in
`scratch/cs4a_r_findings.md`):**

* **EE-6 ⭐ — lift the BINDING BASE first, before the 131/43-site migration.**
  The doctrine's central noun has no type: four hand-wired `__post_init__`
  bodies each remember the guard, and `[M]` (QA-F1) three of four sites had
  no witness until the review added per-site rows. A base/mixin carrying
  `space`, an abstract `data_ng`, and ONE `__post_init__` makes the
  migration inheritance instead of 13 hand-edits, dissolves the
  `operator="Literal"` strings into `type(self).__name__`, gives the C
  binding a semantic `data_ng` operand (the review kept `values.shape[0]`
  — `CrossSectionField.ng` is a MESH read-through CS4b retires), unifies
  the family's incoherent dataclass equality/mutability matrix (MA-4: F
  mutable + space-excluded eq; C `eq=False`; iso pair frozen full-eq), and
  guarantees the S family cannot fall out of the guard wiring at the
  rebind (CEN-4: `ScatteringOperator` carries `mat_xs` + `space` and no ng
  guard TODAY).
* **XD-1 ⭐⭐ — the tightness gate is RE-SPECIFIED; as originally chartered
  its negative control PASSES.** `[M]` (probe `scratch/cs4a_r_probe_
  binding_laws.py`, L=3): under the shipped un-normalized analysis verb
  the adjoint law `bind(K)† = bind(K†)` holds at 2.2e-16 for a
  deliberately NON-tight rule (`M† = R` unconditionally), while
  MULTIPLICATIVITY `bind(K₁K₂) = bind(K₁)bind(K₂)` is what tightness
  certifies (breaks at 1.58 non-tight, 1e-14 tight). The gate lands as
  TWO legs + a blindness control: (i) Galerkin leg — the wrong-EMBEDDING
  negative control (w-weighted vs constant; the measured G6.3 87 % class);
  (ii) tightness leg — multiplicativity on a genuinely non-tight rule at
  ℓ≥1; (iii) record that ℓ=0 discriminates NOTHING (`[M]` both laws hold
  for every rule — §6c). The binding must DECLARE which analysis verb
  (`analysis` vs `project`/canonical-dual) it realizes — the chartered
  equivalence has no truth value until it does. The Σw = 4π placement is
  a VALUE contract, not the adjoint law's hypothesis (`[M]` cancels under
  c ∈ {1, 2π, 137}). ⚠ The truncation-naturality square (Borrowing 3) is
  NOT chartered as a gate: `[M]` the probe's square commuted at 0.0 for
  every rule including non-tight — no demonstrated red, no gate (§6c).
  ⭐ **XD-1 sharpening (2026-08-21, kernel stress test + executed probes):
  the adjoint leg's operand does not EXIST as chartered.** `[M]` (P-1) no
  kernel carries any of {T, transpose, adjoint, dagger, H, conj} —
  `bind(K†)` is unspellable at HEAD. For fission it lacks a TYPE, not just
  an API: `[M]` (P-10) `FissionKernel(chi=νΣf, nu_sig_f=χ)` is REFUSED by
  the kernel's own χ-simplex guard ("EmissionSpectrum is not normalized:
  0.75") — the factor swap leaves the invariant, so the adjoint kernel is
  a DIFFERENT type. The gate lands only after the design decides: kernel-
  level dagger machinery with a typed image, or the adjoint leg re-
  expressed operator-side (`bind(K).H` vs an independently-assembled
  adjoint). Also `[M]` (P-2) fused vs split channel binding
  (`conjugate(Λ+N2N)` vs the operator-sum of conjugates) differ at
  5.6e-17 — FP association only, but the 0-ULP crosscheck pin moves under
  any split; `[M]` (P-3) all four ℓ=0 transfer spellings agree bit-exact
  today (four AGREEING copies — the Pattern-2 debt without a live bug);
  `[M]` (P-4) an over-order `scattering_order` is refused by an accidental
  bare `IndexError`, not a typed guard — the kernel's DERIVED `order`
  subsumes the operator's independent int at the rebind. Full verdict:
  `scratch/cs4b_kernel_as_frame_stress.md` — FrameBase is a 2-field
  zero-math BINDER (the rich object is Basis, 6 array-returning verbs) ⟹
  kernels gain representation-free ARRAY-returning verbs; binding stays
  the external third object = EE-6 lifted to `BoundOperator(datum,
  space)` (deleting the abstract `data_ng`); C stays a deliberate 3+1.
* **XD-2 — the (n,2n) multiplicity single-source obligation, with its
  discriminator.** `[M]` 12 production literal `2`s across
  material_xs_field/iso/cp/moc/mc; the kernel ClassVar is the 13th home
  and unread. The rebind lands a count gate (no production literal outside
  the channel datum — reds at 12 today; a rebind of IsoN2N alone leaves
  11, which is the difference between "single source" and "thirteenth
  home"). Design input, cleared for ≥2 instances: the (Σ_r reaction-XS,
  ν multiplicity, P emission-law) channel triple — `[M]` MC already
  decomposes exactly this way (`mc/solver.py:443-449`: row-normalized
  sampling + w×2), and fission already stores a per-group ν implicitly
  (SigP/SigF). ⚠ `Mixture.Sig2` is an MT=16-only slot (gendf), so
  (n,3n)/(n,4n) have no home — note for the data layer, not this carve.
* **XD-9 — the condensation dual-pair gate is FissionKernel's first
  non-self-referential gate**: `bind(condense(K)) == T·bind(K)·T⁺` iff χ
  is marginalized and νΣf averaged (the ruled pair, `mixture.py:442-448`)
  — red-capable in the average/average direction. Until it lands, the
  kernel's dyad gate compares the kernel against itself, and the latent
  twin (`FissionKernel.dyad()` vs `FissionOperator.kernel`'s
  `outer(chi, production_rate)`) heals at this rebind.

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
