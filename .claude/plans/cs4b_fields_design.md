# CS4b design record — fields are space elements

**Status: ROUND 1 RECORDED** (round 0 = grounding; round 1 = user rulings +
redirections, 2026-08-21; two investigations in flight: field-taxonomy member
census → `scratch/cs4b_field_taxonomy_census.md`, Kernel-as-Frame stress test
→ `scratch/cs4b_kernel_as_frame_stress.md`). Charter:
`space_and_kernel_binding_campaign.md` §5 "### CS4b" (+ the CS4a-R amendments
EE-1/XD-10/EE-5 recorded there). Grounding census (every number below with its
predicate + file:line): `scratch/cs4b_grounding_census.md` @ `466e6756` — cite
it, do not re-copy it (plan-authoring §9).

**Goal (domain terms).** A field is a pair (values, space) whose space answers
every structural question the field has — sizing, measure, membership identity.
The mesh is the space's *provenance*, not the field's *attribute*. (XD-10
bound: SIZING and — per the census, MEASURE — derive from the space; ROLE does
not, role is class identity.)

## ⭐ Meta-ruling (user, 2026-08-21, round 1) — effort is NOT a criterion

Verbatim: *"we never decide on something because it is easier. The difficulty
of something is merely a question of how many sessions it will take. … our
objective is to create the best code. Correct, elegant, ergonomic, efficient.
It will take as long as it takes."* Migration size (the 632-vs-22 framing) is
retracted as a justification everywhere in this file; a convenience surface is
justified only on the four criteria, or accepted solely as a **labeled
intermediate state before full migration**. Recorded durably in
[[feedback-build-the-machinery-operator-algebra]].

## Grounding verdict (corrections to the charter, census §9)

- "10 leaves" → **2 declaration roots** (`BulkField._bases.py:170`,
  `FaceField:782`), 10 ABCs, **20 concrete leaves**. Migration is written
  against the 2 roots.
- "≥11 output-mint sites" → **~107 production construction sites** (~65
  operator-arm / ~28 solver / ~14 promotion; §2). The dunder algebra rides
  free (`replace`-based); the factories and `_from_balance` do not.
- "16 test sites" → **15 direct-ctor sites; 632 factory calls in 86 test
  files** if factory signatures change (§5) — the F4 fork.
- The fabrication's ENERGY half already landed (CS4a K2 `_pose_space`);
  `from_materials` survives as XS-data supplier with ONE production consumer
  (`homogeneous/solver.py:229`); the rebind dissolving it is chartered CS4c —
  the F5 fork.
- EE-1's co-vector partially exists (`ReactionRateFunctional` fiberwise,
  `InnerProductFunctional` generic axis-contraction) — extend, don't mint
  (existence-check done, census §9.6).
- `cross_section_field.py:89-91` false `mesh : SNMesh` docstring — **✅ FIXED
  on sight 2026-08-21** (this session, pre-design).

## The forks

### F1 — WHICH space family do bulk fields get?

Census facts: three disjoint bulk-space mint families coexist (§6c: per-leaf
name+shape TAGS / axis-built `MaterialMesh.bulk_space` / `full_field_space`'s
"sn_bulk" Gram interior); the bulk MEASURE (`mesh.quad.weights` reads at
`_bases.py:409`, `angular_source_sink.py:186`) is not derivable from the tag
family; `has_coordinate_cone`: axis-built → True, tags → None (§10); face
families already read CACHED rich spaces off the mesh — the pattern that works.

- **(A) Axis-built per-role cached spaces on the carrier** *(proposed,
  2026-08-21, unratified)*: `ScalarField.space := mesh.bulk_space` (exists,
  cached); `AngularField.space := mesh.angular_bulk_space` (new
  cached_property: `of_axes(angular axis from quad ⊗ EnergyAxis ⊗ spatial
  axis)` [+ LD moment tail]); the moment family's cell-group factor re-points
  at the cached scalar space. Unifies mint family (i)→(ii); (iii)'s interior
  can re-point at the same object (Gram equivalence V·w — bit-vs-principled
  adjudicated by the byte gates at test-architect time). Measure/ng/N become
  axis reads; cone answerable; the `_space_for_mesh` twin mint RETIRES.
- **(B) Re-point at the existing rich spaces without unification**: angular →
  `full_field_space.interior`, scalar → `bulk_space`. Less construction; keeps
  three families; ng stays positional; measure stays fused inside G.

Why (A): it is the campaign's own direction (CS2 = axes on composites — CS4b
lays the field-layer slice of that substrate), and it makes the census's
sharpened XD-10 (measure not derivable today) actually false afterward.

**Round 1 (user) — F1 WIDENED to a review of the field taxonomy itself.**
Verbatim: *"with space formalized, there is the possibility that we don't need
an angular flux field and a scalar flux field as separate classes, or at least
there is some opportunity for simplification. This is especially true if the
angular space has all the information for a retract operation leading to its
collapse."* I.e. the FAMILY axis (Angular/Scalar/Moment/…) may collapse into
space structure (which axes the space has), leaving classes to carry ROLE only.
Member census dispatched (T1–T5). The adjudication tension already visible:

- **For collapse**: the moment family is the in-tree precedent (L lives in the
  SPACE); the dunder algebra is family-uniform; face families' distinctive
  methods already read `self.space`; per-family `_phase_space_shape` becomes a
  space-axes predicate; units may decompose along axes (the /sr is the angular
  factor's density — units compositional along the tensor structure).
- **Against full class-collapse**: Pattern 4 / the coding-standards decision
  lattice — *"axis changes the SHAPE ⇒ class"*: ψ+φ is today a STATIC type
  error; one class with shape-as-data demotes it to runtime-only. The census's
  T3 (static-typing surface) measures what that costs.
- **The likely synthesis to pressure-test**: family classes survive as thin
  TYPED VIEWS (static ψ/φ distinction kept), while the family MACHINERY
  (validation, integration, factories) collapses into space-driven generic
  bodies + the F3 retract/embed operators. To be adjudicated against the
  census, with the user.

**F1-sub (surfaced by `_from_balance`): per-FAMILY spaces, not per-LEAF.**
Today each leaf tags its own `_SPACE_NAME`, so space identity duplicates ROLE
— exactly what XD-10 says the space must not carry. Under per-family cached
spaces, `_from_balance` collapses to `cls(values=lhs.values − rhs.values,
space=lhs.space)` — no factory, no mesh, no navigation; role transition is
purely class transition (Layer 1 already owns it). *(proposed)*

### F2 — identity re-homing (the §8 surface)

Census §7, measured: space-`==` is a DEMOTION at every family (tags blind to
volumes+BCs; axis-hash BC-blind); two `from_mesh` calls today give `space is`
→ False (uncached mint). ⟹ the honest replacement for the `mesh is` gate is
**cached-space `is`** — which REQUIRES F1's move (the mint lives on the
carrier, cached; factories source it, never mint). *(proposed; near-forced)*

**Round 1 sharpening (mine, after the F3 axis reads — for the user's eye):**
on axis-built spaces, `==` is CONTENT equality (axes' structural bytes incl.
measure — `space.py:161-175`), so the "demotion" verdict was measured on the
TAG family and does not transfer. The BC-blindness of bulk-space `==` may be
mathematically CORRECT rather than a defect: by the DOF-set+Gram criterion the
bulk space genuinely does not depend on BCs — BCs enter the TRACE spaces
(which ARE BC-sensitive), so trace-field partner gates carry that
discrimination. ⟹ the candidate doctrine: partner identity = space CONTENT
equality (axis-built `==`), with cached-`is` as the fast path — provenance
(which mesh instance) stops being an arithmetic gate, content does the work.
⚠ This shifts "different problems don't mix" from provenance-identity to
content-identity — a doctrinal call the user must ratify (CS3 ruling 1 kept
fiber discipline as "space/mesh identity"; this says the SPACE half suffices).

Same-step obligations (§6b — the gate's full call-site set):
- `BulkField._check_partner` mesh-`is` → space-`is` (`_bases.py:196`); the
  `FaceField` copy (`:832`); `ScalarSourceSink.__add__`'s private spelling
  (`scalar_source_sink.py:155`).
- The **17 operator/solver-side `field.mesh is not …` gates** (census §7 list)
  re-keyed per site.
- `FullField.__post_init__`'s `getattr(x,"mesh",None)`-tolerant gate
  (`full_field.py:265-274`) — ⚠ becomes a SILENT NO-OP for migrated leaves;
  re-key in the same step. Elegant form: **FullField becomes an element of
  `FullFieldSpace`** — the composite carries the composite space and the slot
  gate reads `space.interior is interior.space` (+ trace) — the composite-level
  spelling of "fields are space elements". Also re-homes the `FullField.mesh`
  property (`:279`).
- Test pins of the gate messages (~10 files, census §7) migrate with it.

### F3 — sibling-space navigation for cross-family derivations

The problem: `AngularFlux.integrate_angular` mints a ScalarFlux **from
self.mesh** (`angular_flux.py:125`); same shape at
`HarmonicMomentFlux.scalar_flux` (:240), `ScalarSourceSink.as_per_ordinate`
(:198), `truncate/extend` (L→L′). Once fields hold no mesh, the method cannot
reach its sibling space. No third way exists: either spaces are navigable or
the caller/owner supplies the codomain.

- **(a)** The carrier mints a navigable space FAMILY (spaces know their
  marginals) — re-introduces the weld one level down; rejected unless the
  funnel ruling later creates a natural family object.
- **(b)** The maps become bound OPERATORS at the carrier (angular integral,
  moment truncation) and field methods retire/delegate — realizes the algebra
  (the standing lens), but pulls operator-binding machinery into a field-layer
  phase; the binding BASE is chartered CS4c (EE-6).
- **(c)** ~~*(proposed bridge)*: the method takes its codomain space as an
  argument~~ ⛔ SUPERSEDED same day by the user's reframe + the measured
  answer below.

**Round 1 (user reframe + `[M]` resolution): retract/embed, owned by the
PRODUCT space.** User: *"this seems like a retract and embed question. I think
(but I'm not sure) that the space has all the information to retract. I'm not
sure about all the information to embed. But even if it doesn't, if we left an
accessor on space to access the original Discrete Measure, then the
information exists."* Measured against the axis layer (2026-08-21):

- **Retract: YES, fully.** An axis-built space stores `axes: tuple[Axis, ...]`
  (`space.py:196`), each axis carrying its FACTOR MEASURE (`axis.py:16,114-121`;
  `None` ≡ counting, canonicalized; per-factor storage, never the outer
  product). Dropping the angular axis: the remaining axes ARE the marginal
  space (`of_axes(*rest)` — name derived injectively from content, so the
  reconstruction is canonical-by-`==`), and the dropped axis's `weights` IS
  the integration kernel. Conditional on F1-(A): TAG spaces have `axes=None`
  and cannot retract — a further argument for (A).
- **Embed: owned by the CODOMAIN product, so the marginal never navigates.**
  The product knows its factors: which axes the operand's space must match,
  and which axis to broadcast along (with `weights` available for the
  isotropic ÷Σw normalization — today's `from_isotropic` math). Callers
  (solvers/operators) hold the richer space via the carrier.
- **The categorical structure is real**: with the isotropic convention,
  `R ∘ E = id` on the marginal (a genuine retraction pair); `E ∘ R` is the
  isotropic projector on the product — the K_iso family (#276). Condensation
  (XD-9's `T·bind(K)·T⁺`) is retract/embed along the ENERGY axis;
  homogenization along the SPATIAL axes — one primitive, three campaign
  consumers ("build primitives, not products").
- **Realization *(proposed)*: as OPERATORS minted by the product space** —
  `space.retraction(axis_label) → LinearOperator` /
  `space.embedding(axis_label) → LinearOperator` — so `integrate_angular`
  becomes (sugar over) `space.retraction("angular") @ psi`, realizing the
  algebra instead of the welded einsum (`_bases.py:409`). The DiscreteMeasure
  accessor (user's fallback) is noted as the escape hatch if a genuinely
  non-axis measure ever needs to ride along — not needed for the shipped
  families.

### F4 — factory survival

~~632 factory calls vs 15 direct-ctor sites ⟹ factories survive as
conveniences *(near-forced)*~~ ⛔ RETRACTED same day — the framing was
effort-based, which the meta-ruling forbids ("can be accepted only as an
intermediate state before full migration", user).

**Re-derived on the four criteria (round 1, proposed):** the factories
DECOMPOSE — they are not one kind of thing.

- **Pure sugar** (`from_mesh`, `zeros_on`, `from_ndarray`): each saves exactly
  one property read (`space=mesh.<role>_space`) over the primary constructor /
  `Field.zeros(space)`; a second construction idiom beside the primary is a
  Pattern-2 seam with no ergonomic payload ⟹ **retire, full migration**; any
  mesh-delegating stage they pass through en route is a LABELED intermediate
  state, never the destination.
- **Math-bearing** (`from_isotropic` = the isotropic EMBED ÷Σw;
  `from_mesh_and_L` = SH(L) space construction; `from_face_arrays` =
  slot-layout assembly): these are retract/embed/space-construction machinery
  wearing factory clothes ⟹ **re-home to the space/operator layer** (F3's
  operators; space builders), then retire the factory spelling.
- `_from_balance` stops needing any factory under F1-sub (per-family spaces):
  `cls(values=lhs.values − rhs.values, space=lhs.space)`.

### F5 — the homogeneous rebind slice (the "and THEN" tension)

Census §9.5: CS4a K1's kernels have ZERO production consumers; the rebind that
dissolves `from_materials`' last consumer is chartered CS4c; the CS4b charter
says "and THEN the fabricated path retires". Either:

- **(pull, proposed)**: a homogeneous-only rebind slice lands in CS4b — the
  O9 ~10 operator-construction sites in `homogeneous/solver.py` re-point at
  the kernels, `from_materials`' last production consumer dissolves, the
  fabricated path retires as chartered, O1's byte-gate tell completes (D5 8/8
  is the wall). Kernels get their first production consumers (closes the
  unconsumed-machinery gap early).
- **(defer)**: fabricated path survives into CS4c; CS4b's done-when must NOT
  claim the O1 tell (charter edit required either way — §3, edit in place).

**Round 1 (user): REDIRECTED — neither branch taken yet.** The kernel design
itself is re-opened first: *"right now, kernel look like a VERY thin class.
But I think there is some inspiration to be taken from Frame. Frame assembles
the frame object and generated the analysis and synthesis operators. It seems
to me like Kernel could be a more robust class that assembles the linear
operators. You might want to stress test that perspective."* Stress test
dispatched (cross-domain-attacker, P1–P8: generation set, layering inversion,
cross-method reuse, the XD-1 analysis-verb home, the C asymmetry, EE-6
compatibility, the thinness diagnosis, foreign frames) →
`scratch/cs4b_kernel_as_frame_stress.md`.

**Stress-test VERDICT (2026-08-21, sharpest claims probe-verified):** the
Frame analogy is off by one layer — `[M]` `FrameBase` is a 2-field frozen
dataclass (basis, measure) implementing ZERO math (`frame.py:113-132`; both
apply faces one-line delegations); the RICH object is `Basis` (6
representation-free verbs, zero runtime imports). `FrameBase(basis, measure)`
IS `bind(kernel, space)` — **the precedent argues FOR the chartered external
binder**, and kernel-generates-operators (design b) is refuted on two
structural grounds (import direction — a kernels↔scattering cycle once CS4c
re-points; binding is BINARY, neither operand owns it — the third object is
EE-6, lifted to `BoundOperator(datum, space)`). The deciding rule: **a data
object's verbs return ARRAYS; only the binder returns OPERATORS.**

**The user's thinness diagnosis SURVIVES, re-aimed**: kernels are thin
against `Basis` (the 6-verb data analogue), not against Frame — so the
surviving design (c) = thin-data kernels ENRICHED with representation-free
array/kernel-returning verbs (`truncated` ships; dagger-where-typed,
`condensed` per XD-9's ruled pair, channel algebra) + the external binder +
C's deliberate 3+1 (its `Id`-frame absence is `IntegralKernelOperator`'s sole
discriminator — uniformising blinds a working gate).

**Probe results (executed this session, `scratchpad/probe_kernel_stress.py`):**
P-1 no adjoint-family member on any kernel — the chartered `bind(K)† =
bind(K†)` gate's operand DOES NOT EXIST; P-10 the fission factor swap is
REFUSED by the χ-simplex guard — K† of a FissionKernel is a DIFFERENT TYPE
(recorded as the XD-1 sharpening in the campaign plan CS4c block); P-2 fused
vs split channel binding differ at 5.6e-17 (FP association — the 0-ULP pin
moves under any split); P-3 all four ℓ=0 spellings agree bit-exact (agreeing
Pattern-2 copies, no live bug); P-7 both kernel docstring equalities hold;
P-4 over-order refused via accidental bare `IndexError` (untyped — the
kernel's derived `order` subsumes the operator's int at the rebind).

**F5 resolution *(proposed, now safe)*:** with (b) refuted, the chartered
external-binder design STANDS, so the homogeneous rebind slice can pull into
CS4b without double-bind risk: the O9 ~10 sites re-point at kernels through
the binder, `from_materials`' last consumer dissolves, O1's tell completes,
D5 8/8 walls it. The kernel VERB enrichment (the re-aimed thinness) lands at
CS4c with the binding base. Open probe for the test-architect: P-6 (is any
SHIPPED rule a non-tight witness, or is the XD-1 gate's negative leg
custom-rule-only — a §6c question).

### F6 — EE-1's integrated reaction-rate co-vector (obligation, not a fork)

Extend against the shipped pair (`ReactionRateFunctional` fiberwise /
`InnerProductFunctional` axis-contraction) — the integrated pairing ⟨Σ,φ⟩_G
with codomain scalar; re-point the three homogeneous rate lines + the
`.. implements:: normalisation`. Existence-check done (census §9.6); the
extend-vs-new adjudication happens at the execution step with both files open.

## Non-fork obligations (carried from charter + census, all same-phase)

1. The bare assert → typed refusal: `sn/mesh/augmented_mesh.py:322` (probed:
   messageless plain / deep `AttributeError` under `-O`); model on
   `diffusion/augmented_mesh.py:211-218` (the honest O8 half).
2. `.areas` wrong-message (`material_mesh.py:517-523` — 3 arms share
   `_areas = None`, message true for 1).
3. The `mesh is None` two-meanings sentinel (`material_mesh.py:207,492`) —
   the un-weld names the states apart.
4. Docs: `infinite_medium.rst` 4-step narrative re-write (HALF-STALE already
   — `_pose_space` undocumented there; sites census §4f), `spaces.rst:1037`;
   `dead_references` baseline **0** at HEAD — exit must hold it.
5. `_from_balance` flips WITH the factories (verified weld,
   `numerics/field.py:356-358`).
6. Charter number corrections edited in place at the rulings commit
   (plan-authoring §3): the 10/≥11/16 rows.
7. EE-5 grep obligation: every data-flow removal swept over `orpheus/` +
   `docs/` + `tests/`, by CONCEPT as well as symbol.

## Protocol tail (per the ratified per-phase protocol)

Rulings on F1–F5 (user) → fold into the campaign plan §5 CS4b → dispatch
**test-architect** (proactive trigger: the carve crosses
numerics/transport/sn/homogeneous; brief = this file + census) → compact →
execute, surgical posture (main agent writes, user steers).
