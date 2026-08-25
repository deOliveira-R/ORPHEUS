# The collapse pair is frame-induced — the axis marginal rides the stage-2 generator

**Status: ✅ EXECUTED 2026-08-24 (post-compaction, as deferred). Naming RATIFIED
by the user the same day: retraction / section ("I gladly accept the names
retraction and section"). §7's reconciliation ran clean (ancestors confirmed,
census unchanged at zero consumers, IndicatorBasis's node contract verified
friendly to synthetic index nodes, no energy consumers). The §5 commit landed —
see the campaign plan §5 ledger for the hash and suite verdict.**

*(Original status, superseded: RULED 2026-08-24 — execution deferred past a
context compaction; this file was the resume surface for S6.0b.)* Branch `feature/cs1-energy-space` (BRANCH-HOLD
stands). Posture: surgical (main agent writes, user steers). Canonical runner:
`.venv/bin/python -O -m pytest -p no:randomly` SERIAL.

▶ The executing session: read this file whole, run §7's reconciliation, then
execute §5's one commit (S6.0b). The campaign plan
(`space_and_kernel_binding_campaign.md` §5) holds the surrounding S6 order
(S6.0b → S6.1 #399 → S6.2 re-homes); `scratch/s6_design_brief.md` holds the
#399/re-home reconnaissance (its verb-machinery sections are SUPERSEDED by this
file).

---

## 1. Goal, separately from means

**Goal.** The axis-collapse operators (the angular-flux reduction and its
section — and any future single-mode axis collapse) are induced by the Frame
machinery — ONE generator discipline across the tree — instead of being a
hand-derived parallel spelling; the space-level verbs survive as the mint
surface; the operators are the forgetful retention of the generator's output.

**Why (the architecture fact that forced it).** S6.0 (`048144db`, shipped
2026-08-24) built `AxisRetractionOperator`/`AxisEmbeddingOperator` with
hand-spelled kernels. The design dialogue established they are a concept-level
twin of the frame at rank one — `[M]` the exact correspondence is measured (§3)
— which is precisely Cardinal Rule 2's stop condition. The user's summary
ruling: *"a retraction is a special case of an Operator constructed by a Frame
use... the retraction operator is a specialized form of the Frame operator
output."*

## 2. The ruled floor (all ✅ user, 2026-08-24, this dialogue)

1. **The verbs are legit and live on `FunctionSpace`** — the space is "a neat
   place" for them; only the REALIZATION changes.
2. **Frame is the engine.** The retraction is the frame construction with every
   argument determined by "collapse this whole axis": basis = single-region
   `IndicatorBasis` over the axis's index set; measure = the axis's measure.
   The faces specialize: R = the analysis face's content; `R.H` = the
   reconstruction face's content; E = reconstruction ∘ G⁻¹ with **G = the
   frame's 1×1 `discrete_gram`** — the section's divisor is the F-0 Parseval
   metric of the rank-one frame, induced, never a hand convention.
3. **The forgetful-map discipline** (the user's correction that reshaped the
   design): build the frame eagerly AT THE MINT SITE, extract the operators,
   **discard the frame**. Neither the axis nor the operator keeps the
   generator. Retaining the frame — or adding an `Axis`→measure accessor —
   needs its OWN justification (e.g. a second consumer needing the identical
   frame instance for measure-consistency/anti-aliasing), which today does not
   exist: consistency is carried by content-determinism + gates (the F2
   doctrine), not instance sharing.
4. **The stage-2 generator discipline** (user, verbatim — recorded in memory as
   `feedback_stage2_generator_discipline`): *"A stage-2 generator induces
   structure on both the space and the operator, and the two inductions must be
   minted together, at one site. Frame: induces the HarmonicAxis metric (space
   side), mints Analysis/Synthesis (operator side) — consistency is the
   tightness gate. Scheme: induces the trace descriptor and basis kind (space
   side), mints the closure (operator side) — consistency is one closure
   serving both apply and solve, which is ERR-026's structural closing. Mesh
   and Quadrature are the degenerate cases (space side only). Forgetting =
   retaining the induced parts; accessors are provenance."*
5. **Eager-then-available**: the reason the frame-machinery cost objection
   dissolves is WHEN the frame is built — *"if it's such a basic operation, why
   didn't we build the frame before?"* Realized here as lazy-memoized minting
   on the SPACE (one mint per space per axis; carriers cache spaces, so
   `sn.angular_bulk_space.retraction(...)` is warm after first touch — no new
   carrier property needed). `AngularFlux.integrate_angular` becomes
   re-keyable through it later (G6.5 pins that re-key bit-identical).
6. **Cross-method genericity**: the machinery is numerics-level — any
   axis-built `FunctionSpace` from any method's carrier (SN, diffusion, CP…).
   No mesh type appears anywhere in it. Which axes ADMIT is the collapse
   doctrine's business (§4), not the method's.
7. **Names must be canonical** (user: *"only accept the name I say if they are
   the canonically correct names"*). Analysis + proposal in §6.

## 3. The measured floor `[M]` (configurations attached; re-verify per §7 on pickup)

- **The four-square closes under `.H`** — `[M]` 2026-08-24, synthetic 3-axis
  product (angular w=[.3,.7,.5,.5] ⊗ energy counting ⊗ spatial V=linspace),
  inline probe (transcript; re-runnable in 10 lines): `B == R.H` True;
  `w-mean == E.H` True; `mean∘B == id` True; `E∘R == the w-mean projector`
  True. Prose names: R = π\_\* (fiber integration), R.H = π\* (pullback — the
  (π\_\*, π\*) adjunction), P = E∘R = the conditional expectation onto
  axis-constant functions.
- **R ≡ the shipped angular reduction, bit-identical** — gated
  (`test_axis_marginal.py::test_g65…`, np.array_equal on the real SN carrier;
  the kernel deliberately spells the same einsum).
- **E ≡ the `from_isotropic` kernel, bit-identical** — gated (`…test_g66…`,
  np.array_equal; ÷Σw-then-broadcast, same float ops).
- **E is the iso column of the SH frame's physical adjoint** — `[M]`
  `scratch/probe_s6_q5_dissolution.py`: sphere GL L=1 (DIAGONAL Gram, Parseval
  dressed): `face.H(e₀φ) == E(φ)` to 2.2e-16; `reconstruction(e₀φ)/W == E(φ)`
  bit-exact (0.0) on BOTH geometries. Slab L=2 (DENSE Gram): `face.H`
  un-physical by the continuum-metric factor (ℓ=0 ratio 8π = 4π·W_slab) — the
  RECORDED F-0 limitation (CS4c matrix-metric debt), not a blocker: the
  collapse pair rides the measure, not the SH metric.
- **`R.codomain` content-equals the sibling scalar mint** — gated
  (`…test_the_marginal_of_the_angular_bulk_is_the_scalar_bulk`): the angular
  marginal IS `bulk_space` (S1 shares the axes verbatim). This row IS this
  generator's clause-3 family-boundary consistency — half of the tightness
  gate.
- **Frame faces are views, not retentions** — `[M]` read
  `orpheus/numerics/frame.py:526-600`: `_FrameAnalysis`/`_FrameReconstruction`
  hold `frame: FrameBase`. ⟹ a minted operator that kept a face would keep the
  frame; true forgetting = copying the induced data out (§5.2).
- **ULP honesty bounds** (do not over-pin): G6.4 (`R.H == Σw·E`) is nulp-tier
  in general (`[M]` 1.1e-16 synthetic; exact on the §9 production fixture);
  multi-dim-axis `R∘E=id` is 1-ULP (`[M]` 2.2e-16 — flattened 2-D measure
  reassociation); the 1-D angular `R∘E=id` is bit-exact. The gate file already
  pins each at its measured tier.
- **Zero production consumers of the S6.0 pair** — `[M]` grep 2026-08-24: the
  symbols appear only in `operator.py`, `space.py`,
  `tests/numerics/test_axis_marginal.py` (the `test_symmetry.py` hit is an
  unrelated group-theory test name). ⟹ §5's §6b set is exactly those three
  files; the re-carve is one clean commit.
- **F-0's indicator claim** (cited, not yet re-measured here): the recarve
  ledger records "indicator frames get 1/m_R (same theorem)" — for ONE region,
  `discrete_gram` is 1×1 with entry Σw. **Execution must pin this** (§5's
  gram-derivation gate); until then treat as F-0's claim, not this plan's
  measurement.

## 4. Clause gating — the collapse doctrine applied (and the shipped defect)

The collapse doctrine (`docs/theory/foundations/spaces.rst` §"which axes
survive a degeneracy") decides each axis's codomain:

| clause | axis instance | verdict | S6.0b behavior |
|---|---|---|---|
| 3 (compact canonical orbit) | **angle** | axis DROPS; rebroadcast lives on the arrows | admit — the shipped drop-form is correct here |
| 2 (partition-integration, L¹) | **energy** | axis PERSISTS at its one-cell member (Bateman: `⟨σ̄,φ⟩` consumes the partition) | **REFUSE `EnergyAxis`** with the condensation pointer — the clause-2 machinery EXISTS (`EnergyGrid.overlap_to`, the PG condensation frames); minting a drop here would twin it (Cardinal Rule 2) |
| 2 | spatial | persists (one-cell keeps V_cell) | untyped today — see below |
| 1 (non-compact orbit quotient) | homogeneous spatial | NORMALIZE; quotient point persists | lives in the homogeneous carrier; not this machinery |

✅ REMEDIED 2026-08-24 by S6.0b (this plan's §5 commit): a typed `EnergyAxis`
is now REFUSED at the mint with the condensation pointer
(`TestFrameInduction::test_typed_energy_axis_is_refused_with_the_condensation_pointer`);
the untyped-generic admission is the retitled
`TestAxisGeneric::test_untyped_axis_is_admitted_whatever_its_label`.
*(Original hazard, superseded: the S6.0 `retraction("energy")` DROPPED the
energy axis — clause-2-wrong — and the old test pinned the wrong behavior.)*

Untyped generic NODAL axes (synthetic test spaces; the untyped spatial axis)
stay ADMITTED with clause-3 semantics documented — the clause gate becomes
fully structural when CS2's typed axes land (**recorded CS2 hand-off**: the
clause verdict becomes axis-family polymorphism then). Refusing by label
string was rejected (stringly).

## 5. The means (proposed 2026-08-24; mechanics to verify at execution per §7)

### 5.1 The mint — one site, both inductions, forgetful

A module-level constructor at the frame layer (`orpheus/numerics/frame.py` —
the generator's home), shape:

    _collapse_pair(space, axis_index) -> (AxisRetractionOperator, AxisSectionOperator)

Internally, at mint time only: build `DiscreteMeasure` from the axis (synthetic
index nodes, the axis weights — counting when `None`; support tag from the
label; a LOCAL helper, **not** a public `Axis` accessor — ruled floor #3);
build the single-region `IndicatorBasis`; construct the literal
`FrameBase(basis, measure)`; read the induced parts off it — the kernel
weights from `frame.measure.weights`, **the section divisor from
`frame.discrete_gram` (1×1)**; construct both operators together (one site,
both inductions); let the frame go out of scope. If constructing the literal
frame at mint time proves impractical for a corner (e.g. `IndicatorBasis`'s
node/points contract fights synthetic index nodes — UNVERIFIED, check first),
the fallback is: derive the induced data directly + pin the frame equivalence
in the tightness gate ONLY — but prefer the literal construction (consistency
by construction beats consistency by gate).

### 5.2 The operators — the forgetful retention (mostly as shipped)

`AxisRetractionOperator` / `AxisEmbeddingOperator→AxisSectionOperator` (§6)
keep their S6.0 shells — admission, dims bookkeeping, bound product spaces,
shape guards, the einsum/broadcast kernels (G6.5/G6.6 pin them bit-identical
to the shipped consumers) — their retained state IS already the forgetful form
(`_flat_weights`, `_total_weight`, spaces, dims). What changes: construction
routes through §5.1 (the public path — direct `__init__` becomes the mint's
internal); `_total_weight` is set from the frame's `discrete_gram` entry; the
docstrings re-derive from the generator story (+ the π\_\*/π\*/section/
conditional-expectation vocabulary of §3/§6); the class docstring records the
discipline (ruled floor #4).

### 5.3 Caching — on the space, lazily

`FunctionSpace.retraction(axis_label)` / `.section(axis_label)` memoize the
minted pair per axis label in the frozen dataclass's `__dict__` (the F-0
`basis_space` pre-seed precedent proves the pattern on frozen dataclasses).
One mint per space per axis; carriers cache spaces ⟹ warm for SN/diffusion/…
with zero new carrier surface. Both verbs share one mint (both artifacts
minted together — the one-site clause).

### 5.4 Gates (the §6c witnesses exist today)

- All shipped G6 rows survive re-keyed (names per §6).
- **The energy-refusal flip** (§4): the clause-2 test becomes
  `pytest.raises(..., match=...condensation...)` — its witness is the typed
  `EnergyAxis` in the tree today.
- **The tightness gate** (the generator's two-inductions-consistent gate,
  promoted from the transcript probes): the minted pair ≡ the literal
  single-region frame's faces on random inputs (R vs analysis content; R.H vs
  reconstruction content; E vs reconstruction∘G⁻¹), + the shipped
  `R.codomain == bulk_space` row as the clause-3 boundary half.
- **The gram-derivation pin**: the section's divisor equals the 1×1
  `discrete_gram` entry of the literal frame (turns §3's cited F-0 claim into
  this plan's `[M]`).

### 5.5 §6b / sizing

Call-site set: `orpheus/numerics/operator.py`, `orpheus/numerics/space.py`,
`tests/numerics/test_axis_marginal.py` — `[M]` zero production consumers (§3).
**One commit.** Suite gate: the gate file + tests/numerics + neighbor trees
(transport, sn/operators, sn/mesh — the S6.0 battery); pyright 0.

**Done when:** the gate file is green with the energy row a REFUSAL and the two
new gates present; `grep -rn "AxisEmbeddingOperator" orpheus/ tests/` returns
only history/prose (if §6 ratified); the mint function is the only construction
path (`grep` for direct `AxisRetractionOperator(` outside frame.py/tests
returns nothing); pyright 0.

## 6. Naming — the canonical-name analysis (ruled floor #7)

**`retraction` — canonical, KEEP.** The categorical name for the left inverse
of a split pair (split epimorphism; Mac Lane CWM §I.5); the collapse doctrine's
own prose uses "the retract rule". Content name for docstrings: fiber
integration / pushforward π\_\* along the projection.

**`embedding` — NOT canonical for this object → `section` (✅ RATIFIED by the
user 2026-08-24 post-compaction; the renames below are EXECUTED).** Two failures: (a) any injective
structure-preserving map is an embedding — `R.H` (the pullback π\*) is one
too, so the name cannot discriminate the pair our two-types design exists to
discriminate (the ERR-051 weld was hiding in the name); (b) E is definable
only relative to R (`R∘E = id` IS its definition), and the canonical name for
the right inverse of a retraction is **section** (split monomorphism; the
section–retraction pair). Renames if ratified: verb `embedding → section`;
class `AxisEmbeddingOperator → AxisSectionOperator`; "embedding" survives in
prose only as the generic adjective. Docstring/theory vocabulary: R = π\_\*,
R.H = π\* (the (π\_\*, π\*) adjunction), E = the measure-normalized section,
E∘R = the conditional expectation onto axis-constant functions.

~~**If the user vetoes**: keep `embedding` …~~ — no veto: ratified, executed
(`AxisEmbeddingOperator → AxisSectionOperator`; verb `embedding → section`;
zero hits of either old spelling remain in orpheus/, tests/, docs/).

## 7. Pickup reconciliation (plan-authoring §7)

1. `git merge-base --is-ancestor 048144db HEAD` (S6.0 shipped) and confirm
   `2690a434` (S5.4) is an ancestor; branch `feature/cs1-energy-space`.
2. Re-run the consumer census (§3 last bullet) — a consumer may have appeared;
   if so the §5.5 set widens (§6b).
3. Read `IndicatorBasis.__init__`/`evaluate`'s node contract BEFORE writing
   §5.1 (the one UNVERIFIED mechanic — synthetic index nodes).
4. Check whether the naming ratification arrived (§6); if not, ask before the
   rename half.
5. The shipped energy-drop hazard (§4) — confirm no energy-marginal consumer
   appeared in the interim.

## 8. After S6.0b (the campaign order, unchanged) — ✅ ALL EXECUTED 2026-08-24

*(Terminal note: S6.1 `ffb8f286`, S6.2 `78925753`+`53e7d207`, and S7
`1f8e0323`/`2e054bfc`/`26699740` all landed the same day this plan
executed; the campaign plan §5 is the ledger. The section below is the
charter as written, kept per plan-authoring §3. ⚠ One §3 correction
inherited from the S7 audit: this plan's §3 "gram ≡ weights.sum, 8 of 8"
row's universal consequence was REFUTED at GL8 — 1 ULP; see the
campaign ledger's bannered row and correction commit `6734bf15`.)*

S6.1 = #399 (space-derived moment members; witnesses in
`scratch/cs4b_verification_plan.md` §5 — three already-raising LD rows flip
positive). S6.2 = the math-bearing factory re-homes, where the
`from_isotropic` retirement decision RETURNS with E frame-induced (the earlier
deferral: user 2026-08-24, "to discuss the special case, we need the general
case"), plus the `integrate_angular` re-key through the cached pair (G6.5
licenses bit-identically). The χ-class general profile (WeightedIndicatorBasis
single-region — the Petrov-Galerkin test side; its docstring's own second
example is the χ emission collapse) is the recorded LATENT generalization —
same mint, different basis — built only when its consumer lands (CS4c/kernel
adjacency).
