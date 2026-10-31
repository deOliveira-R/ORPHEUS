# Space membership, the co-vector and the kernel gauge — a living plan

Opened 2026-10-30, from layer-DAG Phase 2 (`layer_dag.md`). Ontology being searched: the objects are proposed here and refined with the user before any code beyond the bug fixes below.

## Why this exists

The Phase 2 spec (`scratch/layer_dag/refiling/specs_2/numerics_laws.md` §5c) found a wrong answer. `IntegratedReactionRate.evaluate` on the 2g Mixture-A production rate reads 0.225 for a flat `(2,)` flux, where the exact ⟨νΣf, φ⟩ is 0.095. The cause is in numerics: `InnerProductFunctional.evaluate` lets numpy broadcast a flat carrier against an `(ng, 1)` weight into an outer product.
- The reviews: qa `scratch/layer_dag/ipf_fix/qa.md`, elegance `scratch/layer_dag/ipf_fix/elegance.md`.
- The map: explorer `scratch/layer_dag/ipf_fix/explorer_map.md`.
- What they found: the same hazard sits in `HilbertMetric.pairing` and `FunctionSpace.inner_product`. No space checks that an element is its member.

The user's rulings, 2026-10-30:
- **Membership:** "it seems like the right answer, but it also seems to be elegance debt and should be marked as so … we need a central mechanism that operates correctly by construction to stop these bad constructions. Maybe it's a strong definition of inner product or a development of spaces."
- **Co-vector:** "make your best and principled proposal."
- **KernelGauge:** "Would this be useful for the other methods? It seems generally useful. If so, numerics is probably the right place for it."

## Measured facts the proposal rests on

All `[M]` from the explorer, 2026-10-30, unless marked otherwise. The runtime plugin wraps the space verbs; its control recorded both planted non-members.

**Membership**
1. **No space or metric verb checks membership.** `DiagonalMetric.validate_for` returns None.
   - Membership is checked only on the element side: `Field.__post_init__`, and two `require_member` parsers (`FullField`, tagged #457; `RadialCharacteristicField`).
   - On the operator side, 10 of 67 `apply` bodies check their input. That is a floor.
2. **"Member" has two readings.**
   - Exact shape.
   - Prefix plus payload: the space describes the structural axes, and trailing axes ride through. This is coded in `ZeroMorphism` (2 sites) and stated in prose for `TraceRestrictionOperator`.
   - **Runtime census:**
     - Production traffic: 0 payload-tailed elements.
     - Tests: 7 functions in 4 files pass payload-tailed elements through `inner_product` and the Riesz legs (they test frame-face `.H` and the moment-head metric). One deliberate SN mutant does the same.
3. **The carrier's shape against the weight's.**
   - `InnerProductFunctional.evaluate`: 35 260 of 35 260 SN-run calls have them equal.
   - Leading alignment (`DiagonalMetric`'s `_broadcast_leading`) and numpy's trailing alignment coexist in numerics (X4).
4. **Callers pass bare arrays.** All 6 production `inner_product` callers do. Shape cannot tell apart two spaces of equal shape and different weights; only a typed member compared by content identity can.

**Co-vector**

5. **No dual element exists.**
   - 0 fields on `space.dual()`. There is no Covector class.
   - The Riesz legs ♭/♯ return bare arrays.
   - The theory (`spaces.rst`, `spaces-dual-space`): "in coordinates a functional is an array of V's shape … what makes it a different space is its metric."
6. **⟨1,·⟩_V is spelled 4 ways:** inside `IntegratedReactionRate`, diffusion's `_volume_integral`, homogeneous `space.inner_product`, and SN's `self.volume`.
   - `IntegratedReactionRate` agrees with the space's full pairing in 20 085 of 20 085 calls.
   - But a dense metric overlay on the energy block is accepted (the same overlay on the spatial block is refused). Under it the full pairing reads 83 and `IntegratedReactionRate` reads 46.

**Gauge and singular losses**

7. **`KernelGauge` has one realization:** SN's `LossKernelGauge`, consumed only by `SourceOutcome`, built 3 times, all in SN.
   - The projector itself is already generic numerics (`frame_faces.conjugate` with `CrossGramInverse` over a `GalerkinFrame`). SN supplies only the kernel basis.
   - Three defects:
     - the Protocol's `gauge(trace)` assumes the kernel lives on the boundary trace;
     - `sn/solver.py:804` reads `.blocks`, which is off the Protocol;
     - the metric is read from the legacy `inner_product_weights` slot, not `_resolved_metric`.
8. **Singular losses elsewhere.**
   - Diffusion: all-reflective, Σa = 0, a fixed source gives max|A⁻¹q| = 1.7e17, silently. The control Σa = 0.1 gives 13.2. The exact-zero-pivot refusal never fires.
   - Homogeneous: refuses with `LinAlgError`, but only by pivot luck.
   - CP: rescales `phi *= 1/max(phi)` without recording it.

## The proposals

### P1. Membership: the space owns it, the typed member makes it unspellable

- **What a member is.** A member of a space V is an element whose `.space` IS V by content identity, so a `Field`. This is the by-construction tier. It is the only check that catches two spaces of equal shape and different metric.
  - Every space verb that takes elements (`inner_product`, `pairing`, `norm`, `apply_metric`, `apply_inverse_metric`, the Riesz legs) asks the space one verb, `V.require_member(x)`.
  - That verb is where a `Field`'s space is compared with V.
- **The debt tier.** A bare array is admitted only with `x.shape == V.shape`, exactly.
  - That one check is tagged `ELEGANCE-DEBT[guard]`. Its retirement trigger: the 6 production callers pass Fields.
  - Strict equality closes both the outer-product hazard and the leading-against-trailing alignment hole, so the X4 pair collapses to one convention. [Corrected 2026-10-30, after R1: the raw tier's check runs only under `__debug__`, so under `-O` a raw outer-product call still returns the wrong number (`[M]` 1.4 for 0.68). The closure under `-O` comes from production passing typed members, #664.]
- **The payload convention is a tensor product, named.** An element with a trailing payload is a member of V ⊗ P, acted on by op ⊗ I_P. The proposal: spell it, with `TensorProductSpace` and the identity on the payload factor, at its 2 coded sites and in the 7 test functions, so the space verbs never see a non-member.
  - Each of the 7 tests is re-posed on the explicit V ⊗ P, per the ladder.
- **Passing `cone-membership-is-a-predicate`.** Half 1 holds: every admitted value is legal. Half 2 holds only once the payload is spelled: every legal value is admitted.

### P2. The co-vector: an element of the dual space, paired canonically

- **What it is.** A `Covector` is a member of V* (a field on `V.dual()`), with `evaluate(x) -> float`.
  - `evaluate` is the canonical duality pairing: coordinates contracted, NO metric, x required to be a member of V.
  - The metric enters only in how a co-vector is minted, never in how it pairs.
- **Two mints, kept apart because they differ.**
  - `V.flat(w)`, the Riesz map, reads V's metric.
  - The measure mint: a density σ against V's integration measure, σ·μ. This is how a reaction rate R = ∫ σ φ dμ is born.
  - They agree exactly when V's metric is its integration measure: 20 085 of 20 085 physical-flux calls. They differ under a metric overlay (83 against 46).
  - The reaction rate is the measure mint. `IntegratedReactionRate` becomes the measure mint of its cross section, and ⟨1,·⟩_V's four spellings become one.
- **Consumers (Pattern 6, more than two):** `IntegratedReactionRate`, `ScaleGauge`'s functional, the Phase 2 fixture's `partial(space.inner_product, w)` (one line when this lands), and ⟨1,·⟩_V.
- **What it is not.** `InnerProductFunctional` is a different object: a fibre-wise contraction V ⊗ W → W, the row of the dyad `outer(v, f)`. It stays, under its own name. Whether it should be renamed to say what it is waits for the user.
- **A finding to list:** the energy-block overlay is accepted where the spatial-block overlay is refused, so the two blocks follow different rules.

### P3. The kernel gauge: numerics owns the projector; each method rules gauge-or-refuse

- **The numerics realization.** A `KernelProjectionGauge(space, kernel_basis)`: the G-orthogonal projector Π onto the span of an explicit kernel basis, with σ = I − Π.
  - It is assembled from the generic pieces SN already uses.
  - SN's `LossKernelGauge` keeps only what is SN's: `gauge_freedom`, `_reflection_orbits`, `_build_block_table`, the counting laws. It hands its basis to the numerics projector, so there is one source.
- **Fixed in the same move:**
  - the Protocol parameter is renamed from `trace` to the element;
  - `.blocks` either joins the Protocol or leaves the solver;
  - the metric is read from `_resolved_metric`.
- **Gauge or refuse is the method's ruling, decided by whether its kernel is physical.**
  - SN's kernel is a discretization artefact (diamond-difference face transmission and mirror-odd modes, independent of Σa). Projecting it out is correct.
  - Diffusion's and homogeneous's kernel is physical: Σa = 0 in a closed domain is a steady source problem with NO solution (neutrons accumulate). A gauge there would print a finite wrong answer. The right response is a structural refusal (closed walls and zero absorption everywhere), the analogue of SN's `gauge_freedom` predicate, with an ERR entry.
- **Theory.** The general page (`operator_algebra.rst`, `kernel-gauge-section`) must say that canonicality is the method's claim. The minimum-norm argument is SN's (`cartesian_multid.rst`, `sn-loss-kernel-gauge`).

## Rejected candidates

- **An exact-shape check inside the space verbs, alone:** fails half 2. It refuses the 7 legitimate payload tests, because the payload convention is unnamed.
- **`partial(space.inner_product, w)` as the co-vector:** inherits the outer-product hazard (elegance, `[M]` probe1: 0.285 and 0.225 admitted). It also reads the metric where a reaction rate needs the measure.
- **Gauging diffusion's singular loss:** the posing has no solution, so a gauge answers an ill-posed question with a number.
- **Keeping trailing broadcast as a convention:** it is the hazard itself. 0 of 35 260 calls use it.

## Order

1. **The bug fixes, before Phase 2's moves** (wrong answers first, `process-discipline`):
   - **The co-vector fix**, completed on `fix/inner-product-carrier-shape`:
     - strict `w.shape == x.shape`;
     - `IntegratedReactionRate` checks its carrier against `cross_section.space`, tagged debt;
     - the docstrings, and the stale `homogeneous/solver.py:508-513` prose and its no-op reshape;
     - an ERR entry with its `catches` row;
     - the commit body names the 8473 memo entries it stales.
   - **The diffusion and homogeneous structural refusal**, with an ERR entry.
2. **Phase 2's test moves** (no production code) proceed in parallel with the design discussion. The 4 kernel-gauge rows still go to `sn/`: they test SN's object's conformance.
3. **P1, P2 and P3** after the user's rulings, as W1/W3 builds with test-architect specs.

## ⏸ Open questions for the user

- **P1:** the payload convention as V ⊗ P, made explicit. Is that the right reading?
- **P2:** the canonical pairing plus two mints, with the reaction rate as the measure mint. Should `InnerProductFunctional` be renamed?
- **P3:** confirm the physical-kernel argument: diffusion and homogeneous refuse, SN gauges.

## Rulings and refinements, 2026-10-30 (after the first proposal)

- **R1, the user: guarantees hold only through the formal machinery.** Well-posed mathematics is guaranteed structurally only on the formal path: specialized constructors taking the packed formal objects, and typed members.
  - A raw constructor is legal, but it must be handed everything the specialized one would be.
  - Checks for misuse of a raw path sit under `if __debug__:`, which `-O` strips.
  - Each raw/specialized pair gets one test showing both build the same object.
  - Memory: `feedback_guarantees_through_formal_machinery.md`.
  - Consequence: the canonical `-O` suite never runs a `__debug__` refusal, so its row needs a non-`-O` run or a skip under `-O`, stated as such.
- **R2, the user: diffusion and homogeneous refuse a singular loss** (P3's gauge-or-refuse split, confirmed).
- **P1 refined, after the strongest counter-argument (the payload is a batch, not a factor).** The two readings split by kind of verb.
  - Actions (`apply`, ♭, ♯, `.H`) accept V plus a trailing batch, lifted as A ⊗ I_P. The rule is leading shape = domain shape, and nothing is added to the hot path.
  - Scalar verbs (`inner_product`, `pairing`, `norm`) take exact members only. A batch has no inner product. Today they sum over the payload with counting weights, which drops the payload's own measure (the 0.225 class).
  - Pairing over the payload requires an explicit `TensorProductSpace(V, P)`.
  - Of the 7 payload tests, only the direct `inner_product` callers (`test_axis_metric_is_a_derived_object.py`) re-pose on V ⊗ P. The Riesz-leg and `.H` callers are actions and stay.
  - Status: RULED, the user, 2026-10-30: "I agree with your readings. Let's do those."
- **P2 refined: one mint, not two (the Riesz-uniqueness counter-argument wins).**
  - The `Covector` is a member of V*, evaluated by the canonical pairing, minted only as `V.flat(w)` with V's own metric.
  - On every space with a measure, the metric IS the measure. The leak is closed by refusing an overlay on a counting-measure block (counting is a measure; today the energy block admits a dense overlay: 83 against 46).
  - Spaces whose metric is a frame Gram (moments, frame faces) get reaction-rate co-vectors by pullback along the reconstruction (`E.dual`), never by a second flat.
  - `IntegratedReactionRate` becomes `space.flat(σ)`, and ⟨1,·⟩_V's four spellings become `space.flat(1)`. `InnerProductFunctional` stays: it is a different object, a fibre-wise contraction.
  - Owed before the build: a census of energy-block (counting-measure) overlays in all spaces, production and tests. If one is a legitimate second inner product, the two-mint design returns.
  - Status: RULED, the user, 2026-10-30 (same message).
- **R3, the user, 2026-10-30: "Fredholm as the decision rule, kernel measurement as witness."**
  - **The decision.** Gauge-or-refuse is decided by the Fredholm alternative on the posing (A, q), recorded in the posing's admissibility:
    - trivial kernel: solve;
    - a kernel with q orthogonal to ker A^H: a solvable family, so gauge (the minimum-norm member);
    - q not orthogonal: no solution, so refuse.
  - This replaces the "physical vs discretization kernel" split, which was a heuristic. It agrees with R2 on every case discussed, and gauges a closed pure scatterer with a balanced source (∫q = 0).
  - **Predictors.** The structural predicates (the absorbing-chain predicate for the loss; SN's `gauge_freedom` for the DD kernel) predict ker A a priori and serve production.
  - **The witness.** The measured near-kernel is the witness: σ_min, the right and left singular vectors and the gap, a numerics primitive. A gate pins that predicted dim = measured dim. A posing no predictor covers is measured.
  - **Not refused:** the near-singular continuum (σ_min tracks Σa). Only round-off-level σ_min decides; above it, the answer is returned with its condition number.
  - `[M]` 2026-10-30, dense SVD, diffusion 8-cell mirror slab: Σa = 0 1g gives σ_min/σ_max 9.9e-18, a gap to 2.85e-2, |⟨u_min, q⟩|/|q| = 0.91, and max|x| 1.73e17 ≈ ⟨u, q⟩/σ_min. The control Σa = 0.1 gives 4.7e-2 and 13.2.
- **Corrections and orchestrator decisions, 2026-10-30 (the test-architect's spec, `scratch/membership_covector/spec/spec.md`).**
  - [REFUTED 2026-10-30 for "which spellings are ⟨1,·⟩_V"]: "the four spellings become `space.flat(1)`". The fact:
    - only homogeneous and diffusion are ⟨1_V,·⟩, with the space's UNIT, not `np.ones` (on a frame-dressed head the unit is not all ones);
    - `IntegratedReactionRate` is ⟨σ,·⟩_V = `flat(σ)`;
    - SN's `self.volume` is a per-group spatial integral, which belongs to the spatial axis marginal.
  - **A dense form is admitted only on a MODAL (coefficient) axis**, read from `Axis.kind`. A Gram is a form on coefficients; a nodal axis carries a measure, counting included.
    - `[M]` 22 of 137 positioned forms are reclassified, 0 of them in production (all 108 production forms are on frame-dressed modal heads).
  - **The gauge correction is reported in the G-norm** (U3). Today it uses the Euclidean `np.linalg.norm`.
  - **Production callers are converted to typed members** (U1, U2). `[M]` 6 of 6 scalar-verb sites can be typed; 8 reaction-rate sites take `ScalarFlux(values, bulk_space)`.
    - The residue: the `InnerProductFunctional` inside the rank-one fission operator, unbound until that operator gets its spaces (its own tagged debt).
- **The R3 specification and the build order, 2026-10-30** (`scratch/membership_covector/spec/spec.md` §R3; draft `drafts/test_fredholm_admissibility.py`).
  - **The four cases, `[M]`:**
    - diffusion closed pure scatterer with a positive source: ⟨w, q⟩/|q| = 0.91, refuse;
    - the same scatterer with a balanced source: 1.5e-17, gauge (max|x| 219, agreeing with a bordered solve to 2.6e-14);
    - SN DD box at c = 0.9: the left null vectors are trace-only, so a source with no trace part pairs at about 2e-15, gauge;
    - SN DD box at Σa = 0: union kernel 3 = 2 + 1, and a uniform source pairs at 0.52, refuse.
  - **The left (conserved) basis is not the class indicator.** The groups that feed the class carry the probability of reaching it (`[M]` 0.5 = Σs(1→2)/(Σs(1→2)+Σa₁)).
  - **The U3 projector admits summands on G-orthogonality, not only disjoint supports.** On the Σa = 0 DD box the two kernels overlap on the trace.
  - **The witness is `NearKernel`**, a new unit U3b. Its predictor rows are gated against numpy's dense SVD, an independent oracle.
  - **The near-singular range is solved.** The LU's reciprocal-condition estimate goes into `ExitReport.admissibility`.
  - **Order: U1, U2, U3, U0, U3b.** U0 needs U2's co-vector and U3's projector.
  - The raw-tier debt issue is #664.
- **U1 landing notes, 2026-10-30** (the full `-O` suite on the worktree, `scratch/membership_covector/u1full/driver.log`: 12 failures in 2 targets, every other target green).
  - **The member type lives in `space.py`.** A lazy `Field` import inside the membership check made a `field → space` edge, refused by `test_concept_order` and the layer gate's cycle row. Now `space.py` declares a nominal `Member` (coordinates and space) and `Field` subclasses it.
  - **Three more rows pinned the metric-blind identity**: `test_dense_metric` A9, `test_space::…independent_of_weights` and `test_space_identity_is_structural` G1.4. They are restated: equal-shape spaces with different metrics are unequal and the comparison does not raise.
  - **[REFUTED 2026-10-30 for U1] qa F1, "the raising leg asks V* through `V.dual().apply_metric`".** The fact: a space that types its own carrier (it overrides the metric legs to act on `x.values`, as in `test_coupled_operator._ToySpace`) has a plain `DualSpace` that knows nothing of the carrier, so 7 coupled-operator rows crashed. ♯ keeps asking V (`apply_inverse_metric`). That is a known defect: it admits a typed member of V and refuses a typed covector. **It moves to U2's typed dual elements (#633), with the deleted test `TestTheRaisingLegActsOnTheDual` restored there.**
- **U2 progress, 2026-10-30.**
  - **U2a LANDED `138f7d02`**: a positioned form lives only on a Modal axis (P_modal); 218 dense-form fixture rows re-posed, new gate `test_overlay_is_a_form_on_a_modal_block.py`.
  - **U2b LANDED `64cc3a7a`**: `Member` is a frozen dataclass (`values`, `space`, the shape check) that `Field` and `Covector` inherit, one spelling of the pair. `Covector` lives in `space.py`, not `functional.py` as the spec proposed: `V.flat` mints it and the concept order runs space before functional. `V.flat(w) = Covector(G w, V.dual())`, `V.unit = ⊗ axis.unit`. ♯ asks V* again (F1 restored); the F1 crash was the test double `_ToySpace`, whose plain dual did not type its carrier, fixed in the fixture (`_ToyDualSpace`), no production guard. IPF guard: exact shape under `__debug__`, tagged #664.
  - **U2c LANDED `36dbac9d`**: IRR is `<flat(sigma), phi>` (a cached `covector`), its `density` field retired; SN, diffusion and homogeneous hand it typed `ScalarFlux`es; `FunctionSpace.integral` (flat of the constant function as a typed member) replaces diffusion's `_volume_integral` and the homogeneous <1, phi>. Re-baseline as a DELTA overlay, `tests/gates/transport/carry_rebaselines.json` (event U2c: 36 of 46 eigenvalue-path keys, at most 18 ULP; 0 of 18 fixed-source fluxes), read through `_carry_cases.expected` and gated by `test_the_rebaselines_are_declared_deltas`; `carry_head.npz` is never rewritten. `[M]` per call on identical flux: at most 3 ULP for diffusion (12 cases), 2 for SN (4 bodies). Route gate `tests/gates/cross_method/test_production_pairs_typed_members.py` (the exact verb refuses a bare array; one solve per sharp method). Surprise: the first route-gate probe found the only raw callers were the two `flat(unit)` sites, the space handed its own unit as a bare array, which is what minted `integral`.
  - **Not typed by U2 (recorded, not a defect of U2):** `homogeneous/solver.py` `fission.production_rate.evaluate(phi)` (the dyad's IPF on a raw array) stays on the IPF raw tier (#664).
  - **U2 docs and review LANDED `6aa26d00`**: the archivist's pages plus ERR-121 (energy-block overlay leak) and ERR-122 (the 0.225 read), each with a catches marker verified by an in-process arm. The inline review findings are fixed: the layer gate (rows moved to `tests/gates/homogeneous/test_reaction_rate_on_the_pose.py`), IRR's adjoint through `coordinates_of`, the overlay gate's 64-ULP ceiling and event chain, the #666 tag, the docstrings. Filed: #666 (composite spaces bypass membership), #667 (two volume-integral routes), #668 (identity keyed on the metric's spelling), #669 (the debug rows never run in CI).
  - **Open for the user (qa F1):** under `-O` the raw path `HomogeneousProblem(...).production_rate.evaluate(np.array([0.6, 0.4]))` reads 0.225, as R1's raw tier allows; the formal path refuses it on every run. qa's guard-free spelling `np.sum(np.multiply(c, x, out=np.empty_like(c)))` refuses an upward broadcast on every run and is bit-identical on 2000 of 2000 draws, but still admits a `(1,)` carrier. The default is R1 as ruled.
  - **Next:** U3 (KernelProjectionGauge in numerics; SN's LossKernelGauge rebuilt on it), a surgical carve with the user steering.
  - **U2d is NOT built in U2.** #633 is the flux/source ROLE typing of dual arrows; the issue places it in S6c (ruling 77, moments as covectors). The composite `.values` path in `full_field_space.py` stays, its comment citing #633.

## U1b. One parse at the domain; every extension named (elegance debt, ruled 2026-10-30)

**Ruling (the user, 2026-10-30).** ERR-122 is one symptom of an operation that is not well defined at a deeper level. Fixing it in `Covector.evaluate` alone (qa F1's `np.multiply(out=...)`) is a symptom patch and is not taken. This is elegance debt and is paid now. A root mechanism that protects raw and typed operands alike is not a revision of R1. R1 limits what is GUARANTEED and forbids gates scattered per site; it never forbade robustness at the root. Order: U1b lands before U3, because U3's gauge would be built on the exposed verbs.

**The measurement** (`[M]` the explorer, 2026-10-30, HEAD `6aa26d00`, `scratch/membership_covector/broadcast/census.md`, with the AST scripts, `triage.json`, `probe.py` and the `-O`/debug logs beside it):
- **History.** 6 of 12 recorded incidents share the mechanism: ERR-062, ERR-087, ERR-119, ERR-122, `e42fd498` and `34f465cc`. ERR-055 is adjacent.
- **Census.** 48 of 61 in-scope public verbs (numerics, transport, sn, diffusion, homogeneous) combine an operand with a stored array through broadcasting and have no exact shape check that runs on every run:
  - 9 are debug-only (#664);
  - 31 have no check;
  - 8 have a partial check.
  
  A further 29 internal per-cell kernels are out of scope.
- **Probes.** 32 rows on production objects under `-O`: 20 silent and 12 raising. The silent rows split three ways:
  - 5 outer products, ERR-122's 0.225 among them;
  - 11 implicit pullbacks x⊗1, a different member than the one passed;
  - 4 non-members out.
  
  8 of the 20 stay silent with asserts on.
- **Admission is spelled at least 8 ways.** Each name:
  - `FunctionSpace.coordinates_of` and `leading_coordinates_of`, which take a typed member first, with the raw tier checked in debug only;
  - transport `admit_array`, which is raw only, refuses typed members, and is exact on every run;
  - `_AxisMarginalBase._require`, exact on every run;
  - `require_member` (twice);
  - `FissionMomentOperator._admit`;
  - `harmonic_frame._admit`;
  - `DiagonalOperator._check_shape`, which checks rank only.
  
  A concept shared across this many sites is a Cardinal Rule 2 stop signal.
- **Production holds bare arrays.**
  - `power_iteration` keeps its iterate as a bare array (`eigenvalue.py:392`, `:421`).
  - SN unwraps the field at `sn/solver.py:2043` and never re-wraps it for `compute_group_*_rate`.
  - Homogeneous solves on a raw array at `homogeneous/solver.py:486-487`.

**What is missing, mathematically** `[R]`. NumPy's `*` is the Hadamard product on the broadcast closure of shapes. Applied to two shapes that differ, it silently performs one of three maps. Two of them are legitimate when NAMED; the third never is.
1. **The tensor action A ⊗ I along axes the operator does not act on.**
   - A cross section acts on (energy, space). It acts as M_f ⊗ I_rest on the ordinate axis, and on the LD moment axis, where (fψ)_k = f ψ_k.
   - A metric block or a frame face acts on its own block and carries the payload: A ⊗ I_P.
   - Placement is read BY LABEL at construction (`_broadcast_axes`, ERR-119). For an unlabelled trailing payload, placement is by leading POSITION (`leading_coordinates_of`, `_broadcast_leading`).
   - This is the operand-side identity, not a pullback of the coefficient. On a MODAL axis the pullback u⊗f would multiply only the average slot, which is wrong.
2. **The pullback π* = u ⊗ id of a member onto a larger product.** It exists, named: `AxisPullbackOperator`, the Hilbert adjoint of the axis retraction (`numerics/axis_marginal.py`). It uses the axis UNIT u, not ones.
3. **Never legitimate:** numpy extending the OPERAND to the stored array's shape (x⊗ones, with the wrong unit on a modal axis), or both operands at once (the outer product).

**The design** (to be steered step by step; W3).
1. **One parse.** `coordinates_of` (exact) and `leading_coordinates_of` (leading shape, for A ⊗ I_P) check a bare array's shape on EVERY run; a typed member is compared by content as today. This makes one tier. The ELEGANCE-DEBT[guard] #664 tag retires: the check is the parse, not a guard. Re-scope or close #664 in the same commit. The other admission spellings delegate their shape leg to it and keep only what is theirs (the composite/BulkLift refusal, roles).
2. **Every object that acts on elements parses at its own domain through that one verb.**
   - A bound operator parses at its domain space.
   - An object below the space level owns an axis and parses through the one-axis space that axis induces: `DiscreteMeasure.integrate`, the basis `reconstruct`/`analyze` einsums, `FrameAnalysis`/`FrameReconstruction`, and `InnerProductFunctional`, which binds to its space.
   - An engine is minted by its bound owner with the full domain shape (stage-2 discipline, so `DiagonalOperator`'s check becomes exact). Engines are not checked independently.
   - Per object, decided once in the spec.
3. **The tensor action is placed by label at construction.** `_broadcast_axes` moves into numerics as the one placement mechanism (an operator's axes against its domain's axes). The engine stores the placed coefficient, an `expand_dims` view, so the elementwise multiply and the reduction tree are unchanged. `_broadcast_leading` stays as the honest spelling of A ⊗ I_P on an unlabelled payload, after a leading parse. Verify it is not a weld before touching it.
4. **The production bare-array sites are converted to typed members:** the SN `:2043` unwrap feeding `compute_group_*_rate`, the homogeneous `:486`, `Solution.reaction_rate_density`, and whether `power_iteration`'s iterate is typed. Extend `typed_calls` in `test_production_pairs_typed_members.py` to count before and after.

**Verification.**
- The 20 silent probe rows are the first reds of the new gates, run under `-O`. Seed: `scratch/membership_covector/broadcast/probe.py`.
- The debug-only rows (the ERR-122 catcher, "flat carrier refused") become `-O` coverage. This retires #669's debug-row exposure for them.
- Bit identity is expected. Any new overlay event in `carry_rebaselines.json` means a reduction tree changed, and is investigated, not absorbed.
- Docs that go stale: ERR-122's status line, the two-tier paragraph in `spaces.rst`, and the `coordinates_of` and `leading_coordinates_of` docstrings. Run `dead_references` after any helper retires.

**Ruled (the user, 2026-10-30, AskUserQuestion):** D1 = admit typed: the plain bindings admit a typed member of their end by content; the composite-FullField refusal (BulkLift) stays. D2 = the two-extension reading stands; build on it. Next: the test-architect's spec (`scratch/membership_covector/u1b/spec.md`). **Orchestrator decision, 2026-10-30 (the spec's first question):** design step 2 as written was unspellable. `metric`, `basis` and `measure` sit below `axis` and `space` in `test_concept_order.py`'s TIERS, which counts lazy imports. So the parse splits into its two legs:
- the SHAPE leg, membership of the coordinate space R^shape (exact, or leading for A⊗I_P), is one primitive at the lowest tier;
- the IDENTITY leg, content comparison of a typed member, stays in `FunctionSpace`.

`coordinates_of` and `leading_coordinates_of` compose the two legs. An object that owns its domain shape (a measure's n_points, a basis layout, a metric block) parses with the primitive. An engine that does not own its full shape (`DiagonalOperator`) is minted by its owner with that shape. Option (b), where engines never parse and owners parse one frame up, was refuted for this question: it left `DiscreteMeasure.integrate` silent at its own public verb.

**Asked (kept for the record).**
- (D1) Under the unified parse, should transport's plain bindings (`admit_array`) admit a typed member of their end, compared by content? Today they refuse it with "pass `.values`".
- (D2) Is the two-extension reading right: the tensor action placed by label, and the explicit pullback `AxisPullbackOperator`?

Sizing `[R]`: about 48 verbs plus 3 production sites, in numerics (most), transport, SN, diffusion and homogeneous; 2 to 3 sessions.

- **U1b spec LANDED** (`scratch/membership_covector/u1b/spec.md`, the test-architect, 2026-10-30). The parse sites of the 48 verbs:

  | Parse site | Verbs |
  |---|---|
  | a bound space | 20 |
  | a composite end | 4 |
  | an owned shape, through the primitive | 11 |
  | the table's node axis | 4 |
  | an engine leg minted with its owner's shape | 2 |
  | delegates to its owner | 4 |
  | none of its own | 3 |

  - **New module:** the primitive is `numerics/coordinate_space.py`, in a new first TIERS row.
  - **First reds:** 20 of 20 silent probe rows red today.
  - **Cost:** 0.21 µs per parse, about 0.3 % of an SN solve predicted; not measurable end to end.
  - **Bit identity:** checked three ways: the carry overlay, a 32-key head capture with a positive control, and the SN drift corpus as a delta.
  - **Commit order:** spec §8 (primitive, space parse, below-space objects, faces, engine, D1, observables, docs).
  - **Typed conversion:**
    - homogeneous `:486-487` converts, bit-equal on 3 of 3 mixtures;
    - the group rates and `reaction_rate_density` have no production caller, so the verbs are typed instead;
    - `power_iteration`'s bare iterate is #664's re-scoped residue.
  - **Flags:**
    - F1: the bases' node legs.
    - F2: `RadialCharacteristicReconstruction` has `domain=None`.
    - F3: `TransferMaterialField`.
    - F4: legacy axes-less broadcast weights.
    - F5: the `DiagonalOperator` engine.
    - F6: the Lambertian pair is bound wider than the factor it acts on.
- **Ruled (the user, 2026-10-30, AskUserQuestion), on the spec's NEEDS.**
  - **N1 = one placement.** `DiagonalOperator` takes a coefficient, a shape, and the axes the coefficient spans, all fixed at construction. `axis=k` is no longer a separate mode; it becomes the placement `(k,)`. The 123 test constructions are re-posed with their shape. No rank-agnostic operator remains.
  - **Error type = `NotInCoordinateSpace(ValueError)` only.** A wrong shape is a wrong value of the right type. `admit_array`'s non-array refusal stays a `TypeError`, and the tests that pin a `TypeError` for a wrong shape are re-pinned.
- **U1b progress, 2026-10-30.**
  - Step 1 LANDED `53736b2b`: the primitive.
  - Step 2 LANDED `ab37d078`: the space's parse runs on every run. `[M]` `-O`: numerics 5643, transport 1451, homogeneous 138, diffusion 148, cross_method 207 passed.
  - **Deviation from the spec, §6:** `FunctionSpace.flat` keeps calling `apply_metric`. Composite spaces override `apply_metric`, so routing `flat` to `_diagonal_apply_metric` would bypass their override. The re-parse is harmless, because the parse returns its operand. Spec §6's `bare_actions == 0` counter is therefore dropped; the route gate counts bare actions excluding the `flat` frame.
  - Step 3 LANDED `123166fb`: the five basis contractions are concrete and `@final` on `Basis`; each parses its operand at the table's modes or nodes; the 7 families implement `_hooks`. The hook bodies are 5 einsums repeated across the families, a twin path filed as #671. ⚠ The step shipped with a red Sphinx `-W`: 5 `implements :by:` directives named the now-final family verbs. The archivist re-targeted them to the `_hooks` in step 8.
  - Step 4 LANDED `ef113fdd`: the frame faces and `CrossGramInverse` parse at their bound ends.
  - Step 5 LANDED `112106eb`: `DiagonalOperator(coefficient, *, shape, spans)`, its placement read by label (`FunctionSpace.positions_spanned_by`); `_check_shape`, `_broadcast` and `_broadcast_axes` retired. `[M]` The energy-extent guard must run before the engine mint, or the placement `ValueError` preempts it.
  - Step 6 LANDED `6b4f342e` + `f04c96de`: D1. `admit_array` delegates to `space.coordinates_of`. The pairing route gate exempts `admit_array` by code object.
  - Step 7 LANDED `4882293a`: every production silent admission parses through the primitive. The gate is `tests/gates/cross_method/test_one_parse_on_production.py`: 20 of 20 silent rows red at a `6aa26d00` archive; routing over 25 sites; composite `require_member` strict xfail (#666). `[REFUTED 2026-10-30]` The homogeneous k∞ conversion (`:486-487`): the hub's `production_rate` is measure-carrying, so a bare parse moved k∞ under a weight-2 measure (`test_operator_spaces.py`). The spec's probe ran at unit weight only. It was reverted, so there is no `bare_ipf` counter. `[M]` Full suite green over 33 targets; head capture 32 keys, 0 moved.
  - **Review (qa and the elegance-enforcer, 2026-10-30, on `6aa26d00..4882293a`): no production correctness defect.** Fixed in U1b:
    - The axis marginal held a space but parsed shape-only, refusing a typed `Member` ("values.shape ()") against D1. It now parses at the ACTION verb, `leading_coordinates_of`: it is R ⊗ I_P. `[REFUTED 2026-10-30]` First tried with `coordinates_of`: the pairing spy fired, because the exact verb is the pairings' and an operator's `apply` is an action.
    - The pairing gate's prose claimed the exemption was counted, but the code did not count it. The exemption is now made sound by an AST gate: every `admit_array` caller is an action verb, 6 of 6.
    - Refusal rows added for the two sites that only the route row had guarded: the production rate and `SourcePosing`.
    - Owner tags of the form `repr(space)` prove "at this space", not "at this site". This was accepted and recorded in the gate's docstring.
  - **Filed:**
    - #670: one parse per operand at the operator base. `[M]` 4 parses per `apply_metric`, +1.7–2.7 % per solve. It folds in the multiplier's double parse and the per-call owner repr; blocked by F2 and #666.
    - #671: the basis-hook twin einsums.
    - #672: `FissionMomentOperator`'s domain is (1,) but it acts on (1, ng, nx).
    - #673: `DiagonalOperator.on` plus the ERR-119 residue (the spans are checked by extent, not by label).
    - #674: F4, the legacy weights slot is never validated.
    - Comments posted: #664 re-scoped to the bare iterate (its door is `admit_array`); F6 on #633.

## U3. The kernel gauge in numerics (started 2026-10-30)

- **Ruled by the user (2026-10-30, AskUserQuestion, at the U3 start):**
  - **Summand input = a spanning set.** `KernelProjectionGauge` takes the spanning columns of each summand. Numerics runs the √G-weighted SVD, which does the rank cut and the G-orthonormalisation in one step. SN supplies only its pair generators. The rank tolerance moves to numerics and owes a numerics-side justification.
  - **`LossKernelGauge` is retired.** SN keeps a builder function, `gauge_freedom`, `_reflection_orbits`, the pair generators and the counting laws; the builder returns the numerics object.
  - **The `KernelGauge` Protocol stays,** with its parameter renamed off `trace`. `KernelProjectionGauge` conforms, and test doubles stay legal.
  - **The DD left basis (R3.2's design point) is deferred to U0.** It is designed there together with the absorbing-chain predictor's left basis; U3 is the projector only.
- **Already ruled in the spec (§U3.2 and R3), not re-asked:**
  - overlapping summands are admitted iff mutually G-orthogonal, and a non-orthogonal overlap is refused;
  - `_exit_gauge_trace` reports `‖x − σ(x)‖_G/‖x‖_G`;
  - the metric is read from `_resolved_metric`, and a non-diagonal metric is a declared scope boundary;
  - the solver reads `dimension == 0`, not `.blocks`.
- **Constraints carried from U1b:**
  - the new operator's admission parses at its domain's action verb and gets a row in the routing gate;
  - the summand API is flat indices plus a diagonal weight, so that U0 can pose it on diffusion's full field.
- **U3 code LANDED `2431e2f0`.** SN's gauge is bit-identical on 18 of 18 saved keys.
- **U3 review (qa and the elegance-enforcer, 2026-10-31): no correctness defect.** Ruled by the user (AskUserQuestion, 2026-10-31):
  - (a) the metric verb `HilbertMetric.diagonal(shape)` and `FunctionSpace.metric_diagonal()` land in U3, with composites concatenating in flat order. A dual space is now admitted; it had been refused as "non-diagonal".
  - (b) overlapping summands that are not G-orthogonal stay refused.
  - (c) Π is spelled in the operator algebra (`Σ γ_b.H @ Π_b @ γ_b`). SN's outputs are re-baselined as a declared delta: `[M]` at most 1.2e-16 of max|x|, from the summation order.
  - Also fixed: `KernelSummand` parses its fields; duplicate labels, out-of-range supports and non-1-D spaces are refused; a declared orthogonality to the constants is checked before the zero is stored.
  - Filed: #675 (a flat-index verb on the trace layout) and #676 (`apply` means the section on one gauge and the projector on the other). Commented on #671 (the new basis joins the twin-hook family).
- **For U0, owed by (c):** a `KernelProjectionGauge` on a `FullFieldSpace` (and probably on a `CoupledSpace`) crashes in `apply`. `gather.H` = ♯∘dual∘♭ calls the composite's `apply_metric`, which takes a `CompositeField`, not flat coordinates (#666). It is landed as a strict xfail in `test_kernel_projection_gauge.py`. U0 poses the gauge on diffusion's full field, so U0 fixes it together with #666's flat-coordinate metric.

## U0. The Fredholm decision on the posing (started 2026-10-31)

**The tree against the spec's R3** (`[M]` the explorer, 2026-10-31, HEAD `25d3febe`, grep and Read, `scratch/membership_covector/u0/explorer_map.md`). Three of R3's assumptions do not hold:
- **Diffusion has no source door.** `solve_fixed_source` (`diffusion/solver.py:328`) is a protocol method reached only inside `power_iteration`. Diffusion builds no posing (0 hits for `SourcePosing(`, `EigenPosing(`, `OperatorPencil(` under `orpheus/diffusion/`) and no exit report (`DiffusionResult`, `:389`). A pure scatterer has no eigen posing (F = 0). So R3's balanced-source case and the draft's diffusion refusal rows have no production entry to fire at.
- **R3.4 cannot be written as spelled.** `Asserted.bound` is a required float. `ExitReport(` is built at 1 site in orpheus (`sn/solver.py:683`). `MatrixInverseOperator` exposes no condition estimate, and `lu_factor(overwrite_a=True)` destroys A before ‖A‖₁ could be read.
- **Eigen posings are checkable.** Four of SN's five entries record `NotApplicable("an eigen question has no admissibility")`. Under R3 an eigen posing on a singular A is admissible iff the fission emission is orthogonal to ker Aᴴ.

**The SN DD left kernel is not predicted by the right basis** (`[M]` `scratch/membership_covector/u0/left_basis_probe.py`, dense SVD, worktree `25d3febe`, -O):
- on (3,4) LS2 ng=1 c=0.9 and LS4 ng=2 c=0.9, the left null vectors are trace-only (bulk mass 1e-29) and lie entirely OUTSIDE span(Φ), span(GΦ) and span(G⁻¹Φ) (projection residual 1.00 on all three); the right null vectors lie in span(Φ) (residual 1e-14);
- Lᵀ R = 0: the zero eigenvalue is not semisimple (a Jordan chain), so ker A ⊂ range A. Consequence: a bordered or deflated factorisation A + u (Gu)ᵀ is singular for the SN DD kernel (it needs yᵀu ≠ 0); the projector gauge is unaffected.
- The left basis has two consumers: the Fredholm test of a source with a trace component on a d ≥ 2 all-reflective box, and the adjoint entries' gauge. `[M]` the adjoint entries hand the FORWARD projector to the daggered posing (`sn/solver.py:2918-2920`), which is not ker Aᴴ; this is inert today (the transpose solve is 1-D only, `:2722-2728`), a latent defect for #280.

**Orchestrator decisions (the ontology is known):**
- The verdict (Solve | Gauge | Refuse) is a numerics object beside `gauge.py`: it is the Fredholm alternative on (left co-vectors, projector, q) and knows no neutrons. The absorbing-chain predictor reads `Mixture` and `Binding` and lives in transport; `gauge_freedom` stays SN's.
- Faces conserve iff `BoundaryCondition.return_fraction == 1.0` (the abstract verb, not an `isinstance(WallCondition)` check); decks conserve.
- The class block's pinned solve and the left solve A_TTᴴ w_T = −A_CTᴴ w_C reuse U3's `gather` restriction.
- Order: #666 (flat-coordinate composite metric), the numerics verdict, the transport predictor, SN, homogeneous, diffusion, then the test-architect re-poses the drafts.

**Ruled by the user (2026-10-31, AskUserQuestion), all four as recommended:**
- **Diffusion source door: build it.** A public diffusion fixed-source entry, matching `solve_sn_fixed_source`: it poses `SourcePosing`, returns a `SourceOutcome` with an `ExitReport`, and runs the Fredholm verdict there.
- **Conditioning: in U0.** `MatrixInverseOperator` gains LAPACK's reciprocal-condition estimate (gecon, with ‖A‖₁ read before the factorisation overwrites A). It travels on the Solve verdict into every exit report; the homogeneous result carries the verdict.
- **Eigen entries: every entry.** An eigen posing on a singular loss pairs the emission spectrum with the left co-vectors: refuse if fission emission reaches an absorption-free closed class, gauge if it is orthogonal to the left kernel, solve on a trivial kernel.
- **SN DD left basis: derived in U0.** Its generators are found and gated against a dense SVD, and the adjoint posings get the left projector (the latent defect at `sn/solver.py:2918-2920`).
- **Ruled by the user (2026-10-31, AskUserQuestion, second checkpoint), both as recommended:**
  - **The verdict slot.** `ExitReport` gains `solvability: Solve | Gauge`, the numerics Fredholm verdict with its reason, the source's pairing and the conditioning. `admissibility` keeps the posing's own predicate (the multiplying source's subcriticality) and reads NotApplicable elsewhere. A Refuse raises a `ValueError` subclass at the entry, so no outcome holds one.
  - **The singular factorisation.** Diffusion and homogeneous invert A + Π, where Π is the U3 G-orthogonal kernel projector. `[M]` `scratch/membership_covector/u0/deflated_probe.py`, 3 singular diffusion fixtures: cond(A + Π) 35–42 against 1e17 for A; relative residual ≤ 5.8e-16; agreement with the bordered dense oracle ≤ 9.9e-16. It needs a semisimple kernel (LᵀR nonsingular; `[M]` |LᵀR| 0.68–0.96 on the 3). SN's DD kernel is not semisimple (Lᵀ R = 0), so SN keeps iterate-then-project.
- **The SN DD left basis is derived** (`[M]` the numerics-investigator, `scratch/membership_covector/u0/left_basis/MEMO.md`): ker Aᵀ = span(JΦ), J = diag(−(Ω·n̂_f)·A_f) on the trace, 0 on the bulk. Per block it is `_block_generators` with the column `−sign(Ω·n̂)·sawtooth·Y·checker`. Containment both ways ≤ 1.2e-14 on 9 of 9 cases (d = 2 and 3, graded, unequal weights); the right-generator mutation reads 1.00. Lᵀ R = 0 follows: each mode's in and out values have equal magnitude, so the net boundary current of φ² vanishes.
- **#666's flat half LANDED in the working tree** (not committed yet): `act_on_direct_sum_coordinates` in `space.py`, used by `FullFieldSpace` and `CoupledSpace`'s metric legs, so a composite's metric acts on its flat coordinates as well as its typed carrier. The two U3 strict xfails passed (removed), and #614's two leg-1 rows passed with reciprocity (removed from the blocked set). `[M]` -O: numerics 5739 passed and 3 failed: the 2 #614 XPASSes and 1 `test_traced_memo_findings::q10`, which passed on rerun. sn/operators 1375 passed, transport 1451 passed, layer and concept-order gates passed.
- **U0 step 2 LANDED `54aef2ca`:** `numerics/fredholm.py` (PredictedKernel, Solve, Gauge, NoSolution, fredholm_alternative, deflate, `.H`), `numerics/evidence.py` (the evidence sum left `outcome.py`), `transport/absorbing_chain.py`. The verdict's statistic is the coordinate cosine on each annihilator's DECLARED support: the test-architect measured the cancellation ratio first and found it divides noise by noise when a co-vector carries round-off where it should be zero (`[M]` a dense-SVD DD basis read 0.70 on a bulk-only source); the cosine reads 2.4e-15 there.
- **SN's physical kernel** (`[M]` the numerics-investigator, `scratch/membership_covector/u0/physical/MEMO.md`): right = the class's energy null vector e (K e = 0, K = diag σ_t − S0ᵀ − 2 S2ᵀ), constant and isotropic on every bulk cell and trace slot; left = u_g · w_n · (V_cell on the bulk ⊕ J on the trace), J = −(Ω·n̂)·A_f, u = 1 on the class and u_F = Σ_{g∈C} σ_s(f→g)/(σ_t,f − σ_s(f→f)) on a feeder. Both hold on 9 of 9 cases (d = 1–3, mirror and white(1.0), graded, n2n inside the class) when every present material shares e and u_F. With different materials in the class they fail (residual 0.27): the right vector is a null vector of the class block, the feeder co-vector a transposed block solve.
- **Ruled by the user (2026-10-31): SN's heterogeneous absorption-free closed class is a SCOPE-BOUNDARY.** SN predicts the closed forms; a class whose present materials do not share e and u_F is refused at one declared door. machinery: the class-block null-vector solve through SN's own sweep (U3b's iterative NearKernel route is its natural home). The census found 0 such posings in the test tree.
- **U0 step 3 LANDED `c2a88dfc`:** `MatrixInverseOperator.reciprocal_condition` (gecon on the getrf factors and the 1-norm read before `overwrite_a`). `[M]` the Hager-Higham theorem row is blind to a wrong norm inside its band (ratio 1.87); the definitional row reddens both norm arms.
- **Ruled by the user (2026-10-31): SN's closed-form scope is ONE declared door** (`_require_closed_form_scope`), four entries: a non-Cartesian mesh, a moment closure, a tangential ordinate on a closed face, and materials with different energy vectors. machinery: the class-block null-vector solve through the sweep (U3b). The census found 0 singular test posings outside the scope.
- **Orchestrator decisions in the U0 wiring (working tree, 2026-10-31):**
  - The SN DD left basis is `_balance_weights[support] × right generators`: one diagonal W = w_n V_cell on the bulk ⊕ −w_n (Ω·n̂) A_f on the trace, shared with the physical left co-vector (W(u ⊗ 1)). No second generator loop. `[M]` `scratch/membership_covector/u0/sn_predictor_check.py`: the predicted right and left spans equal the dense ones both ways on 4 of 4 boxes, ≤ 7.4e-14, union included.
  - SN's gauge is on the full field (`_exit_gauge` replaces `_exit_gauge_trace`); the correction is in the full field's G-norm. The adjoint entries gauge with `predicted_kernel.H`.
  - The eigen verdicts probe the range of F (and of F^H on the adjoint) with the emission of a uniform unit flux: every factor is non-negative, so the probe pairs to zero iff the whole range does.
  - `KernelProjectionGauge` admits a multi-axis leaf space on its C-order flat chart (`_flat_chart`): the homogeneous pose is (ng, 1). The 1-D half of the U3 boundary is lifted; the dense-metric half stays.
  - Diffusion and homogeneous mint their bases with `PredictedKernel.of_closed_blocks` (an SVD of each closed block, one solve on its feeders) from the dense loss; they invert `kernel.deflate(loss)`.
  - The homogeneous H4 spy row (0 LU calls on a non-fissile singular medium) is relaxed: the LU of A + Π is legitimate, and the ScaleGauge refusal still fires.
- **U0 review (qa and the elegance-enforcer, 2026-10-31), on the per-entry wiring.** qa found two HIGH defects: the adjoint eigen verdict could never refuse (`[M]` it gauged FIS_FEEDS and returned λ ≈ 1 with a Rayleigh gap of 0.77 where `solve_sn` refused), and a `Solve` record on a loss with a 112-dimensional kernel (product(4,4), 2-D all-mirror). The elegance review found the verdict decided at 9 sites and the outcome's gauge stored apart from the posing's kernel.
- **Ruled by the user (2026-10-31): the posing derives the verdict.** The pencil carries the `PredictedKernel` (`kernel=`, rebased onto the pencil's space; `.H` daggers it); `EigenPosing.solvability` refuses unless the kernel is the WHOLE pencil's (M annihilates the modes and the left kernel annihilates range M; symmetric under `.H`); `SourcePosing.of_pencil(pencil, σ, q)` carries the kernel and `solvability` runs the Fredholm alternative on q; `SourceOutcome.gauge` is read off the posing (the slot is retired). The emission probes and `_solvability` retired; diffusion gains a flat pencil and eigen posing on the solver.
- **Orchestrator decisions with the redesign:** the DD predictor's tangential orbits are a sixth entry of the one scope door (with the void-region entry); `ClosedClass(members, feeders).energy_vectors(mixture)` and `removal_matrix` in transport are the one spelling of the class's energy kernel (SN reads it; the `2.0` literal became `N2N_MULTIPLICITY`); `MaterialMesh.present_materials`; `KernelProjectionGauge.displacement` is the one gauge reading; `FlattenedOperator.apply_transpose`; `of_closed_blocks` asserts nullity one; `.H` uses the metric's pseudo-inverse; `DiffusionSourceSolution.record`; the diffusion source door measures its balance.
- **Filed:** #677 (one coordinate-to-group map), #678 (the flat chart as a FunctionSpace verb), #679 (diffusion's loss and kernel on the hub), #680 (the SN trace metric has no face area). Commented on #665 (a fissile class at round-off removal: λ = −2.3e16 from `solve_sn_adjoint`, unflagged).
- **U0 LANDED (2026-10-31).** Commits:

  | commit | step |
  |---|---|
  | `5ce66815` | a composite space's metric acts on its flat coordinates (#666's flat half; #614 leg 1) |
  | `54aef2ca` | the Fredholm verdict (`numerics/fredholm.py`), `numerics/evidence.py`, `transport/absorbing_chain.py` |
  | `c2a88dfc` | the LU's reciprocal-condition estimate |
  | `ef5eea29` | every posing derives its verdict from the pencil's predicted kernel; SN, diffusion, homogeneous wired; ERR-123 |
  | `119383f9` | qa's re-review: the tangential-only DD kernel's prose; the M ≥ 0 guarantee cited |
  | `f74d7847` | the theory pages, the `verifies` markers, the retirement residue; closes #460 |

  `[M]` -O at `ef5eea29`: sn 4146, numerics 5763, transport 1466, cross_method 262, diffusion 168, homogeneous 149, reference 352 passed; the layer, test-layer, docstring-xref, elegance-ledger and vv-audit gates passed; pyright 0 errors on the changed files. At `f74d7847`: sphinx -E -W exit 0; `dead_references` 3, all pre-existing test helpers.
  - qa's re-review of the redesign: no blocking finding. Both eigen entries refuse together. The Rayleigh gap is ≤ 4.5e-10 on every Gauge return. The multiplying source at σ = 1 matches the dense minimum-G-norm member to 1.5e-13. `fredholm_alternative` has one caller (`posing.py`). qa's F2 is withdrawn: the 112 is the T component, in ker G.
  - Owed to U3b: the iterative NearKernel witness (the SN scope door's machinery, the class-block null-vector solve through the sweep; SN's conditioning reading, now `NotApplicable`). Open elsewhere: #666 (the scalar-verb membership half), #614 (leg 2), #677–#680.
- **Next:** U3b (the NearKernel witness), then one memo regeneration, the merge and CI.
