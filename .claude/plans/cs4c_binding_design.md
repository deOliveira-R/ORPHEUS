# CS4c design record — every operator receives its two spaces and its minted data

**STATUS: LIVING — design round CLOSED 2026-08-30; second-pass review (§12) folded in; F6 ruled option (i). EXECUTION OPENED 2026-08-30: test-architect + step-0 feeding census dispatched in parallel on the user's go. Remaining deferred forks: F-D (awaits step-0), F-E (step 3/4).**
This is the design-round record for the CS-ladder remainder (Campaign 2's
opening act; plan of record `.claude/plans/cs_ladder_remainder.md`; charter
`.claude/plans/space_and_kernel_binding_campaign.md` §CS4c/§CS2). It is a
LIVING document under `plan-authoring.md`: rulings are edited in as they land,
refuted text stays with its refutation, every claim carries its epistemic
marker. Nothing here is executed until the round closes and test-architect
runs (the proactive trigger).

Ground memos (both at `fc60ea64`, 2026-08-30, read-only explorer dispatches):
- `scratch/cs4c_opener_structural_ground.md` — ledger rows, frame-independence
  verdict, leaf surfaces, CS2 residue, #359 remains, the 10-item charter
  corrections list.
- `scratch/cs4c_opener_count_census.md` — the four count re-censuses with
  predicates (S4/S5 DISSOLVED; C-migration 7-vs-133 nuance; R-A roster intact).

Marker key: **[R]** = ruled by the user (dated, quoted where load-bearing);
**[P]** = main-agent proposal, presented and awaiting ruling; **[M]** =
measured (with source); *hypothesis* = unverified means-sketch.

---

## 0. Ground facts the design stands on (all [M], 2026-08-30)

1. **The pull-forward's license HOLDS — via PROVENANCE.** The frame mint needs
   only what the space can answer through the CS5 generator channel
   (`axis.generator_as(Quadrature, consumer=…)`); the S² nodes live on the
   generator, not on the axis; the production bulk space is generator-wired
   axis-built (`augmented_mesh.py:1183`). The live S-side mint to re-route is
   `ScatteringOperator.frame` (`transport/operators/scattering.py:485-512`),
   which today reads the operator's own mandatory `quadrature` field.
2. **The acceptance instrument**: 3 marker constants / 14 strict-xfail
   param-rows in `tests/sn/architecture/test_monomorphic_leaves.py` —
   R1{L,C,S,B}, R2{C,S}, R6{B×8}. CS4c owns R1-C/R2-C/R1-S/R2-S; the CS2
   residue owns R1-L/R1-B/R6.
3. **Dissolved charter inputs**: S5 (909-site sugar) landed `2690a434`;
   S4 (21 field-`.mesh` reads) at 0/0/0; #359's three-spellings premise
   dissolved by P4.9a (one canonical body, `sn/angular/closure.py:1405-1435`);
   CS2 witness 2 (moment mass) remedied by construction.
4. **Binding-relevant landed substrate**: collapse pair
   `FunctionSpace.retraction`/`section` (frame-induced, born bound, S6.0/b);
   `DualSpace` with axes+metric threaded (P7 S2); `HilbertMetric` family
   (P7) whose docstring reserves the Riesz-leg composition; `ScatteringKernel`
   + `FissionKernel` minted with **0 production consumers** — this ladder is
   their first (⚠ second-pass correction 2026-08-30: the family is THREE —
   `N2NKernel` also ships, `kernels.py:222-246`, already carrying
   `multiplicity: ClassVar[int] = 2` as the datum home XD-2 names); L's binding ruled + P4.9b landed (transitional
   `(sn_mesh, spatial_closure, angular_closure)` ctor; end state = the
   4-tuple `(domain, codomain, spatial_closure, angular_closure)` riding O-3).
5. **Dispatch surface**: exactly 3 `singledispatchmethod` sites —
   S (4 arms, no ndarray), C (FullField/ndarray), F (FullField/ScalarFlux/
   ndarray). ⚠ Step-0 census correction (2026-08-30): S's `FullField` arm is
   a RE-DISPATCHER (`scattering.py:1189` re-enters the dispatcher on
   `psi.interior`), so "4 arms" over-counts bodies and under-counts
   branching — construction-time selection relocates that branch one frame
   in; the assembly-owned tier-3 lift is what removes it. The ndarray k-path is alive and protocol-wide
   (`numerics/eigenvalue.py:420` → ndarray-typed realizers in SN/diffusion/
   cp/moc). The 8 `MaterialXSField.apply_*` arms (:743-1023) are consumed
   ONLY by the S family (14 sites) — O-6's move rides the S rebind.
6. **Posing-filtration rulings bind this design** (R13–R23, 2026-08-25):
   CS1.5′ superseded in shape (`Materials` landed; `InfiniteMedium` is
   post-ladder, R23); `MaterialXSField` dissolution ratified (R13);
   B's reshape waits for the crystallized operator concept (R18).

---

## 1. The binding base — [R] corrected: EVERY bound operator carries its TWO spaces; `(datum, space)` was a guide, not the signature

**[R] (user, 2026-08-30, verbatim in substance):** *"space is not enough to
determine domain and codomain. It would be enough naively… but once we have
Riesz operators implemented, even the naive version will die. At the very
least every BoundOperator will need its domain and codomain, and different
amounts of data (the collision here is the easy one, but the scattering will
require kernel and frame for example)."*

**Main-agent corroboration — the tree already proves this four ways:**
- **EE-8's own recorded limit** (charter, CS4a-R): the `(space, ng)` guard
  "structurally assumes an ENDOMORPHISM — a 47g→2g condensation binding has
  no single ng conforming to both ends." The charter knew the one-space
  spelling was the degenerate case.
- **ERR-076** is the founding measurement: an operator that derived NO spaces
  from its factors silently degraded `(K_ω ⊗ I).H` from the Hilbert adjoint
  to the Euclidean transpose (87 % on a production law).
- **XD-7**: `_agreed_space`'s stated expiry — per-leg bindings make agreement
  the wrong law; resolution must become BY LEG. Per-leg = each leg knows its
  own two spaces.
- ⭐ **A SHIPPED witness (step-0 census, [M] 2026-08-30)**: on 2-D
  Cartesian, S at `sn/solver.py:1406` receives a composite whose interior is
  the MOMENT space while its bound `space` names the ANGULAR-interior
  composite — `domain ≠ codomain` live in production, the single `space`
  field naming only the codomain (`scratch/cs4c_feeding_census.md` §3).
  ERR-076 is the historical version of exactly this.
- **L's ruled end state is already the corrected shape**: the O-3 4-tuple
  `(domain, codomain, spatial_closure, angular_closure)`. The user's
  correction makes S/C/F converge to the SAME discipline — one shape across
  the 3+1+1 datum kinds. And the moment route's internal factors (analysis
  `V → moments`, reconstruction `moments → V`) are non-endomorphisms even
  when the composite is one; with `DualSpace` shipped, `riesz_lower: V → V*`
  is non-endo by nature.

**Consequence for the base [P]:** `BoundOperator` is a dataclass ABC carrying
`domain` + `codomain` (both mandatory) + the shared behavior: per-END
admission (the energy-conformity guard runs against the domain's AND the
codomain's energy axis — one check collapses the two for endomorphisms), the
mint-and-forget storage contract, and the declared analysis-verb field
(XD-1's home). Each channel declares its own datum fields — C: coefficient;
S: kernel + minted faces; F: kernel (+ faces on angular spaces); L (later,
O-3): the closure pair. The abstract `data_ng` stays deleted (XD-1 stress
verdict). Endomorphism sugar (`codomain` defaulting to `domain`) may live on
the CLASSMETHOD tier, never on the base.

---

## 2. The frame is CONSTRUCTED OUTSIDE and handed in — [R]; the operator mints its faces from it and forgets it (accessor stays for provenance during prototyping)

**[R] (user, 2026-08-30):** minting the frame inside `ScatteringOperator`
construction *"prevents us from using the same Frame for the Fission Operator
and for the Windowing solution method. In principle the Frame needs to be
constructed and given as an argument, together with domain and codomain. The
Frame should act as an Analysis and Reconstruct Operator factory: during
ScatteringOperator construction it uses the frame factory to produce the
operators it needs, then forgets the frame (an accessor can stay behind for
now for provenance and during prototyping, until we're sure it can be
completely gone)."* Same concept applies to `FissionOperator`.

**This is the posing-filtration §5c interning doctrine, now applied to S/F
[M]:** *"sharing justifies interning — one frame per axis pair, many arrows
per consumer (Windowing minting its own differently-bound arrows from the
same frame is the use case). The cache triad: kernel per cross-section set,
frame per axis pair, arrows per space."* It is also the stage-2 generator
discipline verbatim (Frame/Scheme induce; forgetting = retaining the induced
parts; accessors are provenance). F-1's commit subject already says it:
*"the frame is the shared factory; faces are the bound operators."*

**[P] — the synthesis with the ctor/classmethod discipline (§3):** the EXACT
constructor takes the minted faces (what the operator retains); the
extract-and-mint CLASSMETHOD takes `(kernel, frame, domain, codomain)`,
mints the faces, and feeds the ctor, attaching the frame reference as a
`compare=False` provenance field (the accessor, retirement-tracked). Two
reasons the faces — not the frame — sit on the exact ctor:
(a) the XD-1 gate's negative controls REQUIRE constructing
deliberately-wrong pairs (the wrong-embedding control), which a
frame-welded ctor forbids without doctoring the frame — the exact defect
the readiness assessment flagged ("two shipped gates require a doctored
closure to get its own values");
(b) vv#28's fixture trap closes by gating classmethod-built ≡ ctor-built.
✅ **[R] CONFIRMED (user, 2026-08-30):** *"Your read and my intentions
match."* The frame arrives at the classmethod tier; the faces sit on the
exact ctor; the provenance field rides until proven removable.

✅ **[R] RULED (user, 2026-08-30): the HUB holds the shared frame** — by
symmetry with the P4.9b scheme ruling (the hub keeps the generator; the
operator holds the minted products). One frame per (axis content, L),
shared by S, F, and the Windowing method; consumers that hold only the
space still reach the generator through `generator_as` (the CS5 channel).

---

## 3. The construction discipline for ALL operator families — ✅ [R] RULED (user, 2026-08-30): the three-tier discipline is ENDORSED

**The question [R-posed]:** establish ONE discipline in wide context — do
Bulk/Full/Boundary operators each get exactly their space (bulk/full/trace),
or does everyone get the full space and extract? The user's stated habit:
*"the class constructor takes just what it needs, and a class method takes
more than it strictly needs and extracts exactly what it needs to feed the
construction."*

**[P] The three-tier discipline** (the habit, promoted to law, with the
extraction placed once):

- **Tier 1 — the exact constructor.** Takes precisely what the instance
  RETAINS: its two spaces + its operative data (coefficient / kernel +
  minted faces / closure pair / realized law table). Nothing richer. This is
  Pattern 4 (illegal states unrepresentable), the cheap-fixture surface, and
  what the Riesz legs read (`A* = domain.riesz_raise ∘ A.dual() ∘
  codomain.riesz_lower` needs both ends' metrics — exact spaces make it a
  theorem).
- **Tier 2 — extract-and-mint classmethods.** Take richer objects (full
  space, frame, scheme, kernel, mesh, mixture) and extract/mint exactly,
  then call the ctor. `MultiplicationOperator.from_mesh` is the shipped
  precedent; P4.9b's `.pose(sn_mesh)` is another. R1 (posing arc) is
  honoured: operators stop REQUIRING the mesh (ctor tier) while transitive
  reach at the factory tier stays legal. **Every classmethod owes an
  equivalence gate** (classmethod-built ≡ ctor-built on the same inputs) so
  the exact-ctor fixtures provably represent production (closes vv#28's
  simple-ctor-vs-composite-factory blindness).
- **Tier 3 — lifts between space kinds are explicit verbs, never arms.**
  A bulk-born operator enters a full-space composition through the
  restriction/extension family (`embed ∘ A_bulk ∘ restrict` — extension by
  zero on the trace), owned by the ASSEMBLY at one site, not by a FullField
  arm on every class. The repeated "which block do I act on" conditional is
  the missing type; the verb is the type.

**The per-family answer this yields:** which space a binding takes is a
property of the BINDING, not of the class — recorded per row in the
binding-arity table. The within-group members bind the composite (that is
the space the sum lives on); the k-outer F binds the bulk space (its operand
is a flux distribution, never a trace — ⛔ ASPIRATIONAL: [M] step-0 census
2026-08-30, the tree binds `sn_mesh.full_field_space` at `sn/solver.py:1412`
while feeding bare `(ng,*spatial)`; step 4 makes the bulk binding true); law objects bind trace spaces
exactly (Γ₊ → Γ₋), with the boundary realization applying the crystallized
concept later (R18). At construction, the declared domain selects the ONE
action body — construction-time selection is exactly what the R-A census
forbade **until feeding is normalized**, which is step 0/step 5's job; after
normalization it is the legal mechanism O2 names.

✅ **[R] (user, 2026-08-30), the load-bearing articulation:** *"As long as we
get the essence right (the ctor), then we have a lot of flexibility to change
class methods going forward as the code evolves. Maybe we change our mind on
how the class method should work based on experience enforcing V&V discipline
… but if we get the design right, the core survives."* ⟹ the exact ctor is
the DESIGN COMMITMENT; the classmethod tier is deliberately evolvable without
touching the core. Endorsed as proposed, tier 3's lift-verbs included.

**Angles weighed** (for the record): testing (exact ctors = minimal
fixtures + hand-built negative controls); illegal states (a bulk S cannot
hold a trace it never touches); sharing/interning (factories at tier 2 keep
one mint site per consumer); sum agreement (lifting keeps the within-group
sum single-space today, deferring XD-7's per-leg law to the consumer
campaign); vv#28 (guards keyed on what the object always carries — the
exact spaces); #394's JAX lowering (exact retained arrays are the static
boundary).

---

## 4. Kernels COLLAPSE to one datum per interaction — [R]; operators ask for the truncation they need at construction

**[R] (user, 2026-08-30):** *"The kernels at least should collapse. When the
operators are constructed, they ask for what they need — SN asks for
anisotropy up to whatever order it wants; diffusion asks for just the
isotropic part at construction of its operator."*

[M] `ScatteringKernel` already carries exactly this surface (`p0`,
`truncated(L)`, the ℓ-stack over one `Mixture`). The O7 twin heal follows:
`LegendreMomentScattering`, the iso pair, and `N2NMomentOperator` become
bindings/truncations of the one datum; the 8 `MaterialXSField.apply_*` arms
(S-family-only consumers, 14 sites) become the kernel's array verbs (O-6,
R13's ratified dissolution advancing as a by-product).

⭐ **Step-0 census addenda ([M] 2026-08-30, `scratch/cs4c_feeding_census.md`):**
(a) the satellite mint rate is 1-per-APPLY, not 1-per-solve (up to 911 LMS
instances in one Krylov k-solve) — §4's collapse converts it to
1-per-construction, a measurable cost win with its baseline now recorded;
(b) forward and adjoint S take structurally DIFFERENT internal routes
(`N2NMomentOperator.apply` has zero forward traffic; the forward path builds
the (n,2n) source through `add_n2n_source`) — a twin-path fact the O-6 arm
absorption must price; (c) the iso pair is a measured ASYMMETRIC arrow
(`ScalarFlux → bare ndarray`, 4816 applies, chart-dependent bare-ndarray leg
on sphere/cylinder).

**Open residue [P]:** the operator-CLASS fate. The kernel collapse is ruled;
whether `IsotropicScattering`/`IsotropicN2N` survive as named ℓ=0 binding
recipes (diffusion's spelling) or dissolve into truncated
`ScatteringOperator` bindings is a step-3/4 design decision — decided with
the census's per-binding carrier table in hand.

---

## 5. The solution path is DEFERRED to its own campaign — [R]

**[R] (user, 2026-08-30, verbatim in substance):** *"I don't want to put much
focus on the solution path yet (which consumes what we're doing) — the
solution path will have a detailed campaign, and it will drastically reshape
the machinery of the consumers, because now we understand the ontology of the
consumers much better than when we first designed it."*

Consequences for this ladder: step 5 (arm deletion + #205/#276) takes the
**solver-side adapter** route — operators end whole, the ndarray hatches
close via explicit adapters at the solver seams; typing the k-outer iterate
and every consumer reshape rides the consumers campaign. The within-group
assembly keeps its current shape (tier-3 lift at the existing single site);
nothing in this ladder redesigns iteration/posing machinery.

---

## 6. The tightness gate's adjoint leg — [P] operator-side, explained

The binding's chartered correctness gate exists to catch the ERR-039 family
(a claimed `Π* = R` that wasn't). CS4a-R's XD-1 probe REFUTED the original
spelling — `bind(K)† = bind(K†)` holds at 2.2e-16 even for a deliberately
non-tight rule under the shipped un-normalized analysis verb, so its negative
control would have PASSED. The re-specified gate: (i) the Galerkin leg with
the wrong-EMBEDDING negative control (w-weighted vs constant — the measured
G6.3 87 % class); (ii) the MULTIPLICATIVITY leg
(`bind(K₁K₂) = bind(K₁)·bind(K₂)` — breaks at 1.58 non-tight, 1e-14 tight)
at ℓ≥1; (iii) the recorded ℓ=0 blindness (both laws hold for every rule).

The unresolved half: any adjoint-LAW leg needs the operand `bind(K†)`, and
**the adjoint kernel does not exist as a type** — [M] no kernel carries any
dagger API, and `FissionKernel(chi=νΣf, nu_sig_f=χ)` is REFUSED by its own
χ-simplex guard (the swap breaks the normalization invariant), so F's
adjoint kernel is a DIFFERENT type, not a re-parametrization. Two routes:

- **(a) kernel-side dagger machinery**: mint typed adjoint images
  (`AdjointFissionKernel` etc.). Cost: a type family with no other consumer
  until Campaign 2's pencil possibly wants F† as a first-class datum.
- **(b) operator-side [P, lean]**: realize S/F as composites with
  space-supplied faces (the F4 addendum's mandate), land the Riesz legs
  (accepted suggestion 4), and then `.H` is the reversed composite of the
  factors' adjoints — a THEOREM of the factors. ~~The gate compares `bind(K).H`
  against an independently-ASSEMBLED adjoint (built from daggered factors,
  never through `.H`)~~, plus legs (i)–(iii). Kernel daggers are minted at the
  pencil, with their consumer, if Campaign 2 demands the datum.

  ⛔ **REFUTED 2026-08-30 (test-architect pre-carve round,
  `scratch/cs4c_verification_plan.md` §1.1): the struck comparison is a
  TAUTOLOGY and cannot red under ANY embedding.** [M]
  `(RKM)† = M†K†R†` is an algebraic identity for any three maps and any
  nondegenerate metrics — a wrong embedding enters BOTH sides and cancels
  (measured ≤ 2.24e-16 under correct, constant, AND unweighted embeddings).
  The RULING's substance stands untouched (no kernel dagger; adjointness is
  operator-side, by theorem — the theorem-ness is exactly what the probe
  confirmed). The ERR-039 CATCHER is leg (i), the Galerkin defect on the
  FACES (`‖M† − R‖/‖R‖`: [M] 0.0 correct vs 2.0e-1…3.2 wrong). The identity
  is kept as a documented structural theorem in the gate module's docstring,
  not as an assertion (vv#19: a reading that cannot change is not evidence).

  ⛔ **And leg (ii)'s negative control must be chosen by MEASURED redness,
  not by non-tightness magnitude** (plan §1.2): [M] `gauss_legendre(L)` is
  maximally non-tight (`‖MR−I‖ = 1.0`) yet bit-clean on multiplicativity for
  zonal operands (≤ 5.9e-16, 200 seeds × L∈{1,2,3}); `equispaced_equal` reds
  at ≥ 1.2e-2 and is the control that ships.

✅ **[R] RULED (user, 2026-08-30): route (b) — and the kernel-dagger
question is RESOLVED: no kernel dagger, in this ladder or later, absent a
consumer that can state the missing hypotheses.** The user's argument,
confirmed with one refinement:

- **The user (correct):** an adjoint only makes sense once the kernel is
  realized into an operator — domain and codomain with their Hilbert
  metrics are necessary to define it (⟨Aψ,φ⟩_{G_W} = ⟨ψ,A†φ⟩_{G_V} reads
  both), which is why the operator level pins it by theorem; *"if it cannot
  be pinned by theorem the math is probably incomplete to define the
  meaningful dagger."*
- **The refinement:** a metric-free DATA involution does exist (argument
  swap / matrix transpose / factor swap) — well-defined as pure data. What
  it is NOT is an adjoint: `bind(swap(K)) = bind(K)†` is a CONDITIONAL
  theorem of the binding, with hypotheses only the binding can check
  (counting measure on energy — the CS4a gate that "CANNOT RED, [M] 0.0
  under defect and fix" is this theorem's triviality; tight frame on
  angular — XD-1's multiplicativity finding; nondegenerate metrics). A
  kernel-level `.dagger` claiming adjointness would place the claim where
  its hypotheses are unknowable — vv#28's shape, the F4 addendum's
  Σw-placement warning, and ERR-002's transpose-convention drift habitat
  reopened. Where the swapped data is genuinely needed (inside the
  independently-assembled adjoint), it is array-verb consumption of the one
  datum, never a type.
- **Type-level corroboration [M]:** `FissionKernel(chi=νΣf, nu_sig_f=χ)`
  is refused by its own χ-simplex guard — the type CORRECTLY says the
  adjoint image is not a forward fission datum (the simplex invariant is
  directional physics). And L23(b) already rules "the adjoint is a daggered
  POSING row" — even Campaign 2's pencil consumes bound adjoints, not
  kernel daggers.

---

## 7. Accepted suggestions — [R] all six (user, 2026-08-30)

1. **CS2 witness 1 (BC) re-ruled correct-by-doctrine**: spaces are law-blind
   (F2's DOF-set+Gram criterion); the done-when's BC aliasing probe becomes a
   POSITIVE gate (law-blindness by design); BC-keyed state stays on the mesh.
   Charter edit owed at the CS2 section (in place, per plan-authoring §3).
2. **The S3 flip is a bridge-retirement**: digest names already carry axis
   content through `(name, shape)` equality; the flip makes identity
   structural directly and retires the bridge. Not blocking steps 0–5.
3. **The densifier retirement goes native**: `_tensor_product_factored_metric`
   (P7's new consumer, `space.py:901-905`) reads per-axis diagonals natively;
   `_dense_axes_weights` retires with both clients re-pointed.
4. **The Riesz-leg debt rides this ladder**: mint `riesz_lower`/`riesz_raise`
   on the space wrapping the `HilbertMetric` faces; retire `_AdjointOperator`
   into `A* = domain.riesz_raise ∘ A.dual() ∘ codomain.riesz_lower`
   (`.H` consumers unchanged); #375 expected to dissolve; Λ's factor
   collection rides step 3's frame work. The RESTRICTION verb stays parked
   for the boundary thread unless step 4's composite genuinely needs it.

   ⚠ **Three pre-carve measurements bound this step (plan §1.5/§1.6/§2,
   2026-08-30):** (a) the honest round-trip law is
   `riesz_raise ∘ riesz_lower = P_range(G)` — identity on the bulk,
   the tangential-zeroing projector on the trace ([M] a legal
   `product(4,4)` 2-D mesh has 32/64 tangential slots, trace defect 2.87;
   all four ledger fixtures carry 0 tangential slots and are blind) — a
   `== id` gate is both blind on the corpus and a false red in production;
   a singular-metric fixture SHIPS with the legs (R3 resolved: yes).
   (b) the two verbs live on the PRIMAL space only — `DualSpace.of(V)`
   deliberately carries the primal's metric, so a generic
   `apply_metric`-as-lower on a DualSpace yields G² ([M] `[0.25,4,16]` for
   `w=[.5,2,4]`); never write the double-Riesz involution gate.
   (c) ⭐ the retirement is a measured COVERAGE UPGRADE — splitting the
   paired metric mutation into single legs takes the ledger 9/20 → 20/20
   rows red — with the dual obligation that the flat-metric BLINDNESS
   control keeps the PAIRED mutation ([M] a single leg reads `|1−c|`
   exactly on the flat slab: honest arithmetic, false red); a stale
   `_METRIC_CONSTRAINED` after the split is silent coverage loss.
5. **#359 hygiene now**: comment posted re-pointing the site map (canonical
   body + deliberate gated ULP-twin; premise dissolved; residue =
   close-verify). No code.
6. **`_interior_space`'s unwitnessed refusal**: skipped iff step 3 lands this
   campaign; otherwise owes one witness.

---

## 8. The step ladder (updated hypothesis — re-shaped by §§1–7; step design per phase)

| step | content | ledger rows | notes |
|---|---|---|---|
| 0 | ✅ DONE 2026-08-30 — the feeding census (vv#29), `scratch/cs4c_feeding_census.md` (11 entries × 13 sites × 23 verbs; per-arm activation controls; 11/11 headline numbers non-perturbed) | — | HOMO traffic did NOT move (only the space's provenance did); SN C is minted per-outer and NEVER applied (fused override reads its data — vv#29 mode (d)); #205 ScalarFlux arm corroborated at 0 traffic with an AST-closed reference set |
| 1 | `BoundOperator(domain, codomain)` base (§1) + the Riesz legs + `_AdjointOperator` retirement (§7.4) | — | base admits L later; per-END admission |
| 2 | C exact-ctor mandatory space; `from_mesh` = tier-2 classmethod | R1-C, R2-C | 7 truly-anonymous sites + pin re-keys; #276 arm adjudicated at step 5 |
| 3 | S rebind: kernel + frame-handed-in faces (§2), quadrature field dispenses via generator channel, O-6 arm absorption, iso/LMS/N2N as truncations (§4), XD-1 gate (§6), Λ collection | R1-S, R2-S | extend the pinned `generator_as` AST gate; re-key the `domain is None` pin |
| 4 | F rebind: `FissionKernel` first consumer; composite realization with space faces (F4 addendum; collapse pair ships); XD-9 condensation gate; XD-2 (n,2n) count gate | — | |
| 5 | arm deletion + #205/#276 re-litigation with census evidence; solver-side adapters (§5) | — | construction-time body selection becomes legal here |
| 6 | CS2 residue: S3 bridge-retirement + densifier-native (§7.2/.3); L/B minimal annotation flips; R6 carrier guard at `boundary.py:714` | R1-L, R1-B, R6 | minimal by design — O-3's 4-tuple and R18's B reshape NOT pre-empted; ⚠ three boundary classes carry Optional |
| coda | homogeneous path re-points last (F5); `from_materials` consumer dissolves; O1 completes | — | `InfiniteMedium` (R23) stays post-ladder |

Standing: surgical posture; test-architect BEFORE step 1's first edit; L17
crosswalk first; per-step batteries; predicted-then-measured deltas vs
9950/0 (13 trees rc=0); branch + ff-merge.

⚠ **Second-pass review 2026-08-30 (§12) amends this table**: steps 1+2 fuse
to ONE merge unit (F2); step 2's §6b population is the direct-ctor re-key
set, not the 7 anonymous sites (F3); the R1 reader must be extended in the
same commit as the first field-tier flip (F4); step 3's test-side set is
[M] 99 calls / 18 files (F5); XD-2's count gate moves to step 3 and carries
a user fork on its denominator (F6); XD-1's multiplicativity operand needs a
spelling decision (F7).

---

## 9. Open forks

| id | fork | status |
|---|---|---|
| F-A | the three-tier construction discipline (§3) | ✅ [R] 2026-08-30 — endorsed; the ctor is the surviving core, classmethods deliberately evolvable |
| F-B | frame at ctor vs classmethod tier (§2 synthesis) | ✅ [R] 2026-08-30 — confirmed (classmethod + provenance field) |
| F-C | adjoint-leg route (§6) | ✅ [R] 2026-08-30 — operator-side (b); NO kernel dagger (resolution recorded in §6) |
| F-D | per-binding space table (within-group composite / k-outer bulk / …) | census LANDED 2026-08-30 — observed table at `scratch/cs4c_feeding_census.md` §3 ([M]: bound space is a faithful domain in HOMO only, 1 of 4 families; DIFF needs {composite, scalar-bulk} per binding; SN C's space is a label on data). Ruling at step design. |
| F-E | iso-pair operator-CLASS fate (§4 residue) | ✅ [R] 2026-08-30 — decided at step 3/4 *"with the advantage of hindsight … a sharper understanding of the ontology"* |
| F-F | CS2-residue placement | standing lean (late, step 6) — presented, unobjected |
| F-G | frame interning site (§2) | ✅ [R] 2026-08-30 — the hub |

## 10. Supersessions this round creates (owed edits, executed when the round closes)

- ✅ Charter §CS4c: the `(datum, space)` EE-6 spelling annotated as
  guide-superseded by §1's two-space correction (in place, 2026-08-30).
- ✅ Charter §CS2: witness 1/2 banner at the ontology-tournament bullet
  (in place, 2026-08-30).
- ✅ `.claude/plans/cs_ladder_remainder.md`: §3 staleness list resolved from
  the two censuses; §4 opening-protocol status stamped (2026-08-30).
- ✅ #359: comment posted 2026-08-30
  (github.com/deOliveira-R/ORPHEUS/issues/359#issuecomment-5469489901).

## 11. ▶ NEXT ACT (ruled 2026-08-30): a fresh-context SECOND-PASS REVIEW of this record, then the dispatch

**[R] (user, 2026-08-30):** *"I will give you the opportunity to compact
context and review the plan with fresh context to see if there is something
you'd tighten on a second pass. After you review it, we will do the proper
agent dispatch."*

The review's method (the CS3-R/CS4a-R clear-context precedent, scaled to a
PLAN): read the tree-facing claims fresh — this record top-to-bottom, the two
ground memos, the charter banners — and attack before balancing:
(a) **plan-authoring defects in the §8 ladder**: §6b step boundaries vs
call-site sets (does any step leave a call site speaking the old signature?);
§6c gates landing with their witnesses (what shipped input does each new gate
REJECT on landing day?); §10 done-when reachability (run every tell's grep at
its stated denominator, intersect with the declared untouched set);
(b) **contradictions between the recorded rulings and the tree** (re-derive,
don't trust this file's tense — §7 discipline);
(c) **unstated denominators/predicates** in any ladder row or ruling;
(d) anything the design round missed that a builder would trip on
(the review may propose additions, not only cuts).
Findings tighten THIS document in place ([M]-marked, dated). THEN:
test-architect dispatch (the pre-carve MUST trigger) + step 0, per §8's
standing block — both on the user's go after the review is presented.

## 12. Second-pass review findings (2026-08-30, fresh context per §11) — all [M] at `14201b48`

Method: §11's attack list, run against the tree with fresh context. Every
claim below was measured this session (commands: AST censuses via
`.venv/bin/python`, greps stated inline). Findings are edits-in-place
obligations on the steps they name; none reopens a ruled fork except F6
(one small user fork).

### F1 — step 1's `_AdjointOperator` retirement is underpriced, and its §6b set includes the ACCEPTANCE INSTRUMENT itself

- **(a) `A.dual()` does not exist on operators.** [M] `grep -rn "def dual"`
  → only `FunctionSpace.dual` (`space.py:804`) + `CoupledField.dual_zeros`.
  The §7.4 formula (verbatim reserved at `metric.py:84` ✓) therefore mints
  THREE new objects: an operator-level dual wrapper (apply =
  `inner.apply_transpose`, spaces = the duals) + the two riesz-leg
  operators. The retirement is a re-expression, not a re-plumbing.
- **(b) Four test files import/monkeypatch the private class directly**
  ([M] grep `_AdjointOperator` over `tests/`):
  `test_operator_capability_predicates.py` (:25/:177/:186/:280 — direct
  construction of the wrapper), `test_coupled_operator.py` (:78/:866 —
  monkeypatches `.apply`), `test_inverse_adjoint_coherence.py`
  (:92/:340/:377 — the M-ADJ mutation legs), and ⭐
  **`test_monomorphic_leaves.py` — THE LEDGER — at :608 module-scope
  `_TRUE_ADJOINT_APPLY = _operator_module._AdjointOperator.apply` + three
  monkeypatch sites (:939/:987/:1030, the M-10 family)**. Deleting the
  class kills the ledger file at COLLECTION (vv#17's third pipeline
  failure — reads as 0 failed). ⟹ the retirement commit migrates the
  mutation instruments to mutate the riesz legs / dual wrapper instead
  (coding-standards: retirement = test+marker migration), applied to the
  campaign's own acceptance instrument FIRST.
  ⚠ Pre-carve round extensions (plan §1.7/§1.8, 2026-08-30): [M] there are
  TWO collection-killers, not one — `test_operator_capability_predicates.py`
  constructs the wrapper inside a module-level parametrize list (`:280`);
  and the retirement retires FOUR surfaces, not one: `apply`, `inverse()`
  (the #280 swap law), `is_invertible`, and `apply_transpose`'s refusal —
  [M] the last has ZERO witnesses tree-wide, so its successor spelling owes
  one in the same commit.
- **(c) ~40 docstring references** to `_AdjointOperator` across `orpheus/`
  (`:class:` Python-domain refs — silent at every Sphinx severity) —
  `dead_references` sweep owed at the step, magnitude now stated.
- Day-1 witnesses exist ✓: the adjoint-coherence + factored-adjoint-identity
  gates exercise the new composition on landing.
- Nuance recorded: the unbound-`.H` refusal ALREADY ships
  (`operator.py:1313`, S4-amendment ruling 2026-08-22, with the pointwise
  self-adjoint exemption at `multiplication_operator.py:264-265`) — the
  retirement changes the MECHANISM of `.H`, not its admission.

### F2 — steps 1+2 are one MERGE UNIT (§6c)

The base ABC's per-END admission has no concrete subclass until step 2 —
no shipped input it can reject on landing day. Fuse: land base + Riesz
legs + C rebind on one branch (commits may stay separate inside it); the
R1-C/R2-C flips are the base's first witnesses.

### F3 — step 2's "7 sites" is the WRONG-PREDICATE count (§2 filter)

The 7 counts space-ANONYMOUS direct ctors. The step as ruled (§1: base
carries `domain`+`codomain` as mandatory fields; endomorphism sugar
classmethod-only) re-keys EVERY direct ctor call: [M] **26 direct C sites**
(3 production + 23 test, census 1's own numbers re-bucketed), while the
127 `from_mesh` spellings survive ONLY IF tier-2 keeps a one-space
`space=` parameter — a design decision the record did not state. ⟹ stated
now: **tier-2 classmethods keep the endomorphism one-space spelling** (the
sugar §1 licenses); the §6b population per step is the DIRECT-ctor re-key
set + attribute-read set (`.space` reads: [M] 161 Load-context sites in
`orpheus/`, mixed operator/field receivers — the per-step census
discriminates receivers at design time).

### F4 — §10: the R1 acceptance reader CANNOT SEE the success shape

[M] `_domain_annotation` (`test_monomorphic_leaves.py:680-705`) walks the
MRO for a PROPERTY (`vars(klass).get(prop_name)` → `fget.__annotations__`).
A mandatory dataclass FIELD never appears in `vars(klass)` → the reader
returns `("<not found>", "<not found>")`, which contains neither "None"
nor "Optional" → **the strict xfail xpasses VACUOUSLY** — identically for
a correct flip, for a field still annotated `FunctionSpace | None`, and
for `domain` being deleted outright. ⟹ the first flip commit EXTENDS the
reader (also read dataclass `__annotations__` through the MRO; treat
`<not found>` as FAIL), in the same commit as the SUT change, so the
instrument keeps discriminating. (§10's question — "what does the
instrument print on success?" — answered: PASS, but also PASS on two
failure shapes.) Pre-carve round: [M] independently confirmed and worse —
the reader also passes on `field(default=None)`, kw-only fields, and a
deleted attribute (plan §3.1). R6 RESOLVED from §1's own ruling: the base
is a dataclass ABC carrying domain+codomain as FIELDS, so the reader
extension (plan §3.2 G-B1 + meta-test G-B2) is MANDATORY.

### F5 — step 3's test-side denominator, now measured

[M] AST census over `tests/` (call nodes, direct + `from_*`):
**99 construction calls in 18 files** — ScatteringOperator 7 direct +
8 `from_solver_data`, LegendreMomentScattering 19, IsotropicScattering 23,
IsotropicN2N 18, N2NMomentOperator 4, FissionOperator 2 + 18
`from_solver_data`. Step 3 (≈71 S-family calls) and step 4 (≈20 F calls)
now carry their populations. ⚠ census 1's variable-call/registry-loop
check (the 2026-08-29 §6b spelling inventory) was run for C ONLY — step
3's opener census owes the same check for the S family.

### F6 — XD-2's count gate is §10-third-shape DESIGNED-RED at its stated denominator ⟶ one small user fork

[M] the multiplicity-2 literals at HEAD: `material_xs_field.py`
:1020/:1036/:1067, `isotropic_scattering.py:448`, `cp/solver.py`
:568/:629/:697/:784, `moc/core.py` :184/:316/:366 (+1 MC site per the
charter's own census) — **7–8 of ~12 live in cp/moc/mc**, solver families
in this ladder's untouched set. As chartered ("no production literal
outside the channel datum") the gate cannot reach green from the
SN/diffusion rebind alone. Fork:
- **(i) [lean]** the step includes the mechanical ~8-line re-point of
  cp/moc/mc literals to `N2NKernel.multiplicity` — datum CONSUMPTION
  (recycling the vocabulary, not improving those methods on their own
  terms — consistent with the sharpening-order law);
- **(ii)** the gate's denominator shrinks to SN+diffusion with the
  cp/moc/mc residue recorded in the gate's own docstring.

✅ **[R] RULED (user, 2026-08-30): option (i)** — the step includes the
mechanical re-point of every multiplicity literal to
`N2NKernel.multiplicity`; the gate's done-when is 0 production literals
outside `kernels.py`, at the full ~12-site denominator.

⚠ Pre-carve round corrections to this finding (plan §1.3/§1.4/§6.3):
- [M] the denominator is **14**, not ~12 — `material_xs_field.py:809/:862`
  were missing from every prior list; and two members evade a `2.0 *`
  filter (`moc/core.py:316` is an INTEGER `2 *`; `mc/solver.py:447` is
  `w *= 2.0`, an AugAssign) — the count gate's predicate is validated
  against all four spellings (plan §6.1).
- [M] the XD-9 pair's path: `orpheus/data/macro_xs/mixture.py:437-450`
  (the charter's `transport/mixture.py` path never existed at HEAD), and
  `condense` has a second (bilinear) branch obeying a different law.
- [M] §6d run on the re-point: cp/moc/mc → transport is layer-LEGAL
  (L3 → L2) but all three edges are 0 today and cost ~254 ms cold import
  each (eager `transport/__init__`); MC's spelling hoists to a module
  constant (single source — a literal with an exclusion would be the
  thirteenth home again) unless the user objects.
- R5 RESOLVED: the two new sites sit inside the arms O-6 absorbs at step
  3 — they die by absorption; the other 12 re-point; the count gate lands
  at END of step 3 and is green from landing day.
Also: XD-2's SUBJECT (the N2N truncations) lands at step 3; the record
placed the gate at step 4 — ⟹ **the count gate moves to step 3** so gate
and witness land together (§6c).

### F7 — XD-1's multiplicativity operand has NO SPELLING

[M] `kernels.py` API: `p0` + `truncated` only — no product/compose/
`__matmul__`. The gate's `bind(K₁K₂)` operand needs a decision:
a kernel composition verb (Funk–Hecke: composition of zonal kernels =
the ℓ-stacks' elementwise product — domain-true) vs gate-side array
construction from the moment stacks. Lean (defer-until-2): the GATE
builds the product kernel from the stacks; mint the verb only when a
second consumer appears. Named here as a test-architect design input.

### F8 — small corrections and step-note additions

- The kernel family is **three** (§0.4 corrected in place): `N2NKernel`
  ships with `multiplicity: ClassVar[int] = 2` (`kernels.py:224`) — the
  n2n datum is ALREADY collapsed; §4's O7 heal consumes it rather than
  minting it.
- Step 5's notes gain the [M] monkeypatch pair
  `tests/homogeneous/test_homogeneous.py:166/:174` (patches
  `MultiplicationOperator.apply` — an arm-deletion §6b member, census 1's
  own warning).
- Step 6's notes gain: (a) the `_R6_XFAIL` marker text cites
  `boundary.py:343` verbatim — re-point at the flip; (b) the L/B
  annotation flips owe a constructor-site census (who passes
  `space=None` to `StreamingOperator` + the three boundary classes)
  before any annotation turns non-Optional.
- §7.4's "(`.H` consumers unchanged)" — TRUE for the public spelling
  ([M] 149 production / 381 test `.H` mentions, no signature change);
  FALSE for the four private-importing test files — covered by F1(b).

### §12-bis — the pre-carve verification round (2026-08-30): plan, census, and ruling status

Both mandated pre-carve artifacts landed 2026-08-30 (both at HEAD
`2f44ed4e`, numbers re-measured there — cite the artifacts, never copy):
**`scratch/cs4c_verification_plan.md`** (test-architect; 17 gate rows,
5 batteries all ≤ 24 s, cumulative prediction ≈ 10 000/0/19sk/227des/56xf
at the coda) and **`scratch/cs4c_feeding_census.md`** (step 0; folded into
§§0–4, 8–9 above).

Further pre-carve facts the step designs consume (plan § refs):
- **Ledger coverage gap (plan §1.10):** the R1 instrument covers [M]
  4 of 8 classes carrying an Optional space annotation, in two spellings —
  `IsotropicScattering`, `IsotropicN2N`,
  `RadialCharacteristicBoundaryOperator`, `SNMaskedBoundaryOperator` flip
  UN-GATED; steps 3/6 either extend the ledger rows or record the gap in
  the flip commit.
- **`_R2_XFAIL`'s stated mechanism is present-tense-false for S** (plan
  §1.9): an anonymous `ScatteringOperator.H` now RAISES `MissingAdjoint`;
  only C degrades (via `is_metric_free_adjoint=True`) — the marker's
  evidence text rides the flip.
- **§6b members no AST call census can see** (plan §8.4): 3
  variable-mediated constructions (`test_isotropic_scattering.py`
  :124/:183/:220) + 2 monkeypatch surrogates (`test_homogeneous.py`
  :166/:174) + the set-literal registry pins (`test_kernels.py:660-671`).
- **The tier-2 equivalence-gate roster is larger than three** (plan §4):
  `StreamingOperator.pose` + the three kernel `from_mixture`s also owe
  G-C1 rows. R7 RESOLVED: yes — §3's ruled "every classmethod" discipline
  admits no datum-tier exemption.
- **dead_references baseline: [M] 0 dead / 52 checked** — every step's
  exit re-runs it against this zero.

**Ruling status:** R3 ✅ (singular-metric fixture ships — §7.4 rider),
R4 ✅ (legal; MC hoists to a module constant), R5 ✅ (absorption),
R6 ✅ (fields ⟹ reader extension), R7 ✅ (from_mixture owes gates).
**R1 ✅ [R] CONFIRMED (user, 2026-08-30: "Corrections proven by math are
always welcome")** — the tautological leg is dropped; Galerkin-on-faces is
the catcher; the identity stays as documented theorem.

**R2 ✅ [R] RULED (user, 2026-08-30): the NAMED-COMPOSITE direction — the
AdjointOperator survives as a FIRST-CLASS operator realizing the dagger
arrow.** ~~My own earlier proposal — `.H` returns the composed product,
the swap law a theorem of reversed factor inverses~~ ⛔ REFUTED at design
time: [M] `OperatorProduct.inverse()` is a wrap-delegate
(`InverseOperator`), demoting the #280 swap law from object identity; and
role passthrough dies (`test_psi_half_coupling:2895` pins
`a_ab.H.system_role`). The ruled design, with the user's own articulation:
*"during the construction of the AdjointOperator we literally construct an
operator with domain and codomain swapped — the domain of A becomes the
codomain of A†, and the codomain of A the domain of A† — and this exchange
is realized in code."* Maximum-leverage consequences, all step 1:
- apply routes through three first-class legs built at construction
  (`codomain.riesz_lower`, `A.dual()`, `domain.riesz_raise`) —
  bit-identical (G-A1);
- `A.H.H is A` becomes an OBJECT IDENTITY (adjoint() on the adjoint
  returns `inner`);
- `apply_transpose` becomes a THEOREM of the legs
  (`(A*)ᵀ = ♭_W ∘ A ∘ ♯_V`, metrics symmetric by admission) — closing
  **#375** (both defects: the dead-end AND the 0-witness raising stub);
- the swap law stays #280's object identity; `dual(dual(A)) = A` object
  identity too;
- the legs are public space verbs (`V.riesz_lower : V → V*`,
  `V.riesz_raise : V* → V`, PRIMAL-only — `DualSpace` REFUSES them, making
  the G² trap unspellable);
- the class is PROMOTED to public `AdjointOperator` in its own mechanical
  follow-up commit (3-search audit + dead_references before/after).

⚠ §7.4's "retire `_AdjointOperator`" is hereby SUPERSEDED-IN-MECHANISM:
the retirement target was the INLINE metric recipe in its apply, never the
arrow object; the chartered outcome (adjoint = algebra, legs
single-sourced and individually mutable, DualSpace consumed) lands whole,
plus the dagger-law identities a deletion could never have carried.
