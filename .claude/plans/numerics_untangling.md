# The `numerics` layer's concepts in one acyclic order — a living plan

Opened 2026-10-18 from the verification-tree campaign's sub-cone measurement (`verification_tree_refinement.md`, finding S13; ruling 27 there). Status: **the ontology is being searched** (`plan-authoring` §0). Nothing here is implemented. Implementation starts when the user rules the plan polished. Execution is to be delegable (the verification-tree plan's ruling 18), so the plan's final form cuts each step into packages with an acceptance check.

## Goal, in the domain's terms

The `numerics` layer holds the mathematics every method rests on: axes, discrete measures, bases, the frame that pairs a measure with a basis, quadratures, function spaces, symmetry groups and invariance, operators. Each concept is defined from the concepts before it and never reaches a concept after it, so its claims can be verified using only what comes before. The concepts form one acyclic order with joins: for example, the Discrete Measure and Basis concepts both feed into the Frame (the verification-tree plan's ruling 25). That order is the base of the whole verification cone, so every sub-cone above rests on it being right.

## The user's instructions (verbatim, 2026-10-18)

> "2. settle it first." (the numerics knot, before the census)

> "3. It's fine that the SN is essentially a big blob depending on SNProblem *for now*. This is a fairly strong indication that we have not been diligent enough with our rule of "__init__ constructor being the raw constructor for testing and class methods the specialized constructors that receive convenient arguments that needs unravelling", so we need to revit it in its own campaign. Untangling the numerics layer is a pretty good start."

## The rulings by dependence (read this first; the ledger below is the chronological record)

The rulings are numbered in the order they were given. They depend on each other in the concept order, which is not that order, so this index regroups them. Each group is one concept of the order (G0's `TIERS`, `tests/gates/numerics/test_concept_order.py`), upstream first. Inside a group a ruling comes after the rulings it rests on. A ruling's text lives once, in the ledger or in its wave's section; this index holds a one-line gloss, what the ruling rests on (`rests on`), and its state:
- `[LANDED <hash>]`: the code is on `main`;
- **in flight**: being built now;
- **owed**: ruled, not yet built;
- **open**: a question still to be ruled.

Glosses are summaries. When a ruling's exact wording matters, read it at its source.

**How to use it mid-implementation:** before changing a concept, read its group and every group above it that its rulings rest on. A ruling that rests on an **owed** one cannot be built before it, unless its own text says how it stands in meanwhile (for example ruling 53 on degenerate metrics before S6).

### P. How the campaign works (no concept dependence; they govern every group)
- **3.** The test DAG's sub-cones are the code's ontology, so the concept order and the base of the verification cone are one object.
- **14.** The code is never the oracle. What SHOULD be is argued from the mathematics, the ontology and elegance.
- **19.** The order is derived from what is joined to what (rests on 3), never ruled edge by edge. The import cycles were the signal.
- **1.** The knot is settled before the verification census. **2.** The raw-`__init__` rule is its own campaign (#611), applied here where an edge breaks it.
- **28.** Split execution: M packages are delegated; for an S package the main agent writes the design, the user rules, and the main agent writes the code with checkpoints. **29.** The census interleaves, upstream first. **30.** No compatibility re-exports.
- **5.** The operator module is cut, not the concepts. This is the cut's charter, `[LANDED e60330d2]` with M7.

### A. Symmetry and the G-space (the support of a measure)
- **9.** Test both readings and think broadly. It produced 13.
- **13.** The support carries the group action; the measure contributes only equality up to relabelling; no `DiscreteMeasure` subtypes (rests on 9).
- **43.** The action is typed by an abstract `GSpace(Manifold)` (rests on 13). **48** amends it: point identity (`coordinates`) is on every `Manifold`.
- **47.** A motion acts on M iff g(M) = M, one law, no special case (rests on 43).
- **46.** A total G-space verb `ambient_points` (rests on 43).
- **42.** The invariance verbs stay on the measure as compositions, and `invariance.py` is deleted (rests on 13, 43, 47).
- **44.** The bare-interval tests are re-posed on S²/O2_x (rests on 42, 43).
- **49.** The measure's two support guards are elegance debt, #624: the machinery owed is a measure typed by its support (rests on 42).
- State: `[LANDED 7621fc5f]` (S5). Owed from here: #624.

### B. The measures: recurrence, discrete measure, quadrature
- E′ finding "Gauss is not a join" (no ruling number): the recurrence is upstream of both kinds of measure. `[LANDED 2071c98e]` (M1).
- **31.** The orbit space's reference measure is `pushforward_reference`, in `generating_measure` (rests on A's quotients). `[LANDED 659a4656]` (M2).
- **21.** A quadrature is a discrete measure with a fibration; its cached frames are a cache, not part of the concept. It sits before the axis.
- **8** then **12.** The quadrature part of the measure family is cleaned (rests on 21). **32** is the surviving rule list, Quadrature's classmethods; amended: `folded_product` is not registered. **33.** One name per level. `[LANDED e60330d2]` (M6).

### C. The frame (the join of basis and measure)
- **16.** A frame is the sampling map of a basis on a measure, with its dual. The choice of measure and basis fixes the representation; their fineness fixes the resolution.
- **22.** The Gram datum is the frame, at the join of measure and basis, upstream of the axis (rests on 16, 19). Its faces are chain steps after the operator.
- **17.** Every axis carries the pullback of its factor's own measure; "integrates order L" is tightness against the continuum metric (rests on 16, 22). Standard S_N stays; Galerkin quadrature is open as #613.
- **34.** The frame's faces are downstream constructors, `FrameAnalysis.of(frame)` and its siblings (rests on 22). `[LANDED 1ecbc64a]` (M3).
- **58.** The Galerkin test-space predicate is the identity `frame.test is frame.basis`, not `isinstance(frame, GalerkinFrame)`: a `PetrovGalerkinFrame(b, m, b)` whose test family is its trial family IS a Galerkin frame and gets the dressed test space (the user, 2026-10-20: "what the mathematics favor"). Lands in S7; the tree has the isinstance behaviour until then (M3).
- **59.** `Axis.generator` stays as provenance for now (the user, 2026-10-20: it "proves itself quite useful"); the user asked for a better proposal (S7's design section).
- **63. S7's spaces are minted on the space side, through their axis** (the user endorsed the orchestrator's proposal, 2026-10-20). The sampled space of a measure is `FunctionSpace.of_measure(measure, label)`, the space over `Axis.of_measure(measure, label)`; it replaces `DiscreteMeasure.space`, which today builds L²(μ) inline as a bare `FunctionSpace(name="L2[…]", inner_product_weights=w)`, a second spelling of the axis-built head. The coefficient space of a basis is the space over `Axis.of_basis(basis, label)` (the modal axis carrying the unit, the continuum Gram and the index shape, read off the basis with no `space` import); the typed spaces keep their `for_basis` mints, and the indicator families gain one. Generic consumers (`frame_faces.basis_space(frame)`) reach the right space through ONE dispatch, `coefficient_space(basis)` in `spaces/`, where each typed space registers the basis it spans and an unregistered basis gets the plain axis-built space. `basis.py` and `measure.py` then import no space module, closing G0's three S7 `PENDING` edges. **Owed before code:** (a) a census of `.space` reads on a statically unknown basis, which decides whether the registry is needed at all; (b) what `SphericalHarmonicSpace` and `LegendreSpace` add beyond `from_L` and `isotropic_slot`; if nothing structural, they collapse into the axis-built space and the registry disappears; (c) whether the axis-built L²(μ) has the same content identity as today's inline one (a re-baseline question if not).
- **64. `Axis.generator` is split by what consumers read** (the user endorsed the orchestrator's proposal, 2026-10-20; refines 59). `[M]` 2026-10-20, an AST pass over `orpheus/` (excluding the homonym in `derivations/common/withdrawal.py`): 6 production reads. Four recover a quadrature's directions (`mu_x`) to build streaming or a harmonic frame (`sn/angular/closure.py:1675`, `:2236`; `sn/mesh/reduced_operator.py:516`; `transport/frames/harmonic_frame.py:441`). One is a presence guard (`harmonic_frame.py:483`). One re-truncates a moment head through the basis or frame (`spaces/moment_head.py:119`). (1) **Directions are content:** a nodal axis carries its nodes beside its weights, and the four direction reads read `axis.nodes`. `[R]` Today axis equality excludes the generator, so two rules with equal weights and different directions (a rule and its reflection) mint equal axes and equal spaces, and a field bound to one passes a space check against an operator built on the other; a probe is owed to confirm this before the change. (2) **Provenance becomes a typed closed sum** `Origin = FromMeasure | FromQuadrature | FromBasis | FromFrame | Counting`, read by an exhaustive `match`, replacing the four-way union and `generator_as(kind, consumer=…)`. The axis `kind` is derived from the origin (NODAL from a measure or a quadrature, `Modal(unit)` from a basis or a frame), so `Axis(kind=NODAL, generator=basis)` becomes unspellable. What remains of provenance is the genuine section of the forgetful map: re-generation at another order (`at_order`).

### D. The axis
- **4.** The axis comes after what defines it; its generator is kept as provenance. **22** names the generator: an axis is a forgetful image of its measure or its frame (rests on B and C).
- **open:** the forgetful-map concept (ruling 4's "until a forgetful-map concept is chosen"). No issue yet.
- **7** then **11.** `mass_matrix` retires onto the frame's Gram; the axis receives its generator; mints move to the axis side (rests on 22). **35** splits it: M4a and M4b `[LANDED cac6d401, 4e1f2182]`; the coefficient space of a basis and the sampled space of a measure are **S7, owed** (G0 PENDING: `("space","basis")`, `("space","measure")`, `("spaces","basis")`).
- **6** then **10.** The space holds the operators extracted from the frame, not the frame; the axis marginalises itself, modal and nodal alike (rests on 17, 22). **36** refines 10: the axis owns its unit u, not the covector. **37.** The unit u is minted with the axis. **38.** The section divides by the correctly rounded mass. **39.** The unit is content on modal axes only. **40.** A modal axis is minted with its unit. **41.** The axis kind is the closed sum `NODAL | Modal(unit)` (rests on 39, 40). `[LANDED 659a4656]` (S1).
- **15** then **18.** The curvilinear LD refusal is machinery not yet built; curvilinear LD gets a Gauss-nodal spatial axis (rests on 17). **Owed** as #612.

### E. The space, the metric, the dual and the adjoint's posedness
- **25.** Whether the adjoint belongs with the algebra waits on whether it is always posed (rests on 22's metric chain). Settled by 26.
- **26.** A degenerate metric is replaced by its live subspace, and the dual space carries G⁻¹, so the adjoint is posed everywhere; the adjoint, the dual and the Riesz operators stay with the algebra (rests on 25). **Owed** as S6 (#617).
- **27.** The lying adjointability flags are fixed, pinned by one honesty gate (rests on 25's census). `[LANDED 31f14c93]`.

### F. The operator algebra (the sum, the product, the zero)
- **45.** The zero belongs to the algebra: `ZeroOperator` and `ZeroMorphism` live in `operator.py`, and `0·A` refuses until S2. `[LANDED 126cf514]`.
- **23.** A sum is a formal linear combination of terms, unordered; a product is a flattened ordered word (rests on 19 and on L(V, W) being a vector space).
- **50.** A sum's inverse is named by its split: `GreenOperator(M, N)` for A = M − N, and a sum of two or more terms has no inverse of its own (rests on 23, and on 20 below, whose first half it brings forward).
- **51.** `ScaledOperator` retires into the combination; `{A: 1}` collapses to A (rests on 23, 45).
- **52.** `+` never fuses; the fused sweep composites become atoms, built explicitly by their owner (rests on 23, and on 20's "the realisation is chosen by the Strategy").
- **53.** The adjoint distributes over the composites (rests on 23, 26, 27). Its precise form, the orchestrator's reading of 53 from the S2 gate spec's measurements (`[M]` probes `p1_pointwise_dense.py`, `p2_word_adjoint_singular_middle.py`): the DUAL distributes over every sum and word, exactly and with no metric; the ADJOINT is ♯ ∘ dual ∘ ♭, taken once at the composite's ends. `Σ cᵢ wᵢ.H` and `B.H @ A.H` are theorems only where every term and every middle space has a non-degenerate metric, so they hold everywhere only after S6 (26). Today's whole-composite adjoint is exactly this form, so its bits do not move.
- The orchestrator's calls inside 51–53 (S2 gate spec, 2026-10-20): terms collect by object identity; the formal arrows (adjoint, dual, inverse) collect by (kind, the inner's key, the ends); float coefficients are not exact, and that is recorded, not hidden by a tolerance; every zero is absorbed, including the scalar 0 for `sum()`; `(αA)⁻¹` seeds its inner solve with α·x₀.
- **54.** A pin moved only by the flattening of a named sum is re-baselined, that row alone, under the three conditions; the named sum keeps its name (rests on 23, 51). First instance: the coda record `homo_2eg_n2n`.
- **55.** A block or system role lives on the LEAVES: the site that builds an atom stamps it, and every composite derives its role by the join; the zero map a combination mints inherits the join of the cancelled terms' roles (rests on 23; the S2 review's F1).
- **56.** A multiplier refuses a space without labelled axes rather than guess which axes its coefficient spans (rests on 9's space-first reading; the S2 review).
- **57.** The two disagreeing seed paths of a scaled inverse (αA)⁻¹ move to S4, which owns how the Strategy realises an inverse (rests on 20).
- State: `[LANDED 4e5ee36b]` (S2). Owed from here: 57 and the seed paths to S4; the fused atoms' assembly note to S3.

### G. After the algebra: assembly, the inverse family, the Strategy
- **24.** Assembly is a homomorphism out of the algebra: each leaf returns its raw matrix, and `[Σ cᵢTᵢ] = Σ cᵢ[Tᵢ]` folds over the terms (rests on 23). **Owed** as S3. It removes G0's PENDING edge `("assembled_operator","operator")`.
- **60.** `as_matrix` is the operator's matrix representation, fixed by its spaces' coordinates (the axes in order), so it stays on `LinearOperator` and is realised by probing alone. Structural assembly is a faster algorithm for the same matrix and lives downstream as `dense(op)` in `assembled_operator`; `MatrixInverseOperator` calls `dense` (the user, 2026-10-20, S3 checkpoint: "if space knows its axis, and LinearOperator knows its space, LinearOperators should have no problem assembling the associated matrix"; resolved as: the spaces fix where each entry sits, the algorithm fixes how it is computed).
- **61.** The domain composites that emit an assembly (`CoupledOperator`'s block law, `BulkLift`'s embedding, `FlattenedOperator`'s pass-through) are atoms to the free algebra: each keeps its law in its own module and recurses through the downstream `assemble` functor (the user, 2026-10-20, S3 checkpoint).
- **62.** A leaf's emission is a raw sparse matrix and only the functor wraps it; the three `FlattenedOperator` identity pins are re-posed as value pins, and the registry row that passes the result on as an operator wraps it through the functor (the user, 2026-10-20, S3 checkpoint).
- **20.** The inverse is formal: a sum's inverse is `InverseOperator(sum)`, realised by the Strategy's splitting and driver; `InverseOperator` keeps its name (rests on 23). Its first half is in S2 (ruling 50). Its second half is **owed** as S4: the Strategy's driver realises `InverseOperator(sum)`, and Green's adjoint becomes Green of the inner adjoint. That second half rests on 53. The S2 commit removes G0's PENDING edge `("green_operator","operator")`.

- **65. The Strategy mints the Green** (the user, 2026-10-20, S4 checkpoint; refines 20's second half). `[M]` the S4 census (`scratch/numerics_untangling/packages/S4_census.md`): `GreenOperator(` has 28 call sites, 1 in production (its own `dual`), 27 in tests; the Strategy realises (M−N)⁻¹ by hand at 6 driver sites (`SourceIteration(seeded_inverse(M), N)`), so the tree holds two realisations of one inverse and only the hand-wired one runs. The Strategy (its `Splitting` A = M − N and its driver choice) realises the formal arrow `InverseOperator(A)` into a `GreenOperator`, and the 6 sites consume that Green. `A.inverse()` on a multi-term sum keeps refusing, naming the Strategy route: an inverse with no Strategy in scope is unspellable. The plan row's premises were stale: the lazy Green import is gone and G0 holds no S4 edge.
- **66. A Green's adjoint stays the base ♯∘dual∘♭ until S6** (the user, 2026-10-20; amends 20's adjoint clause). Its `dual` is already the exact `Green(Mᵀ, Nᵀ)`. `GreenOperator(self.inner.H)` cannot be built (the split cannot be read back from `(M−N).H`), and the distributed `Green(M.H, N.H)` is the honesty gate's independent witness (`test_adjointability_is_honest.py`, `_sn_green`); an override of that form was removed on 2026-10-18 because the witness then agreed with it to exactly 0.0. Revisited at S6, where ruling 53's distributed form holds everywhere.
- **67. One seed law for a scaled inverse** (the user, 2026-10-20; settles 57). `(αA)⁻¹` owns its seed law once (forward α·x₀, transpose likewise); `seeded_inverse` honours the seed it is handed instead of dropping it (`_SeededExactApply`, `iteration.py:291`); gates pin the forward, the transpose and the driver-slot paths, each with a reddening mutation.

### What this order says to build next
Every ruling in A to D has landed except S7 (35's second half, with C's open predicate) and #612, #613, #624. What is owed is in E, F and G, and the dependences fix the sequence:
1. **S2** (F: 23, 50–53). It rests only on landed rulings, plus 26 where its own text says how it stands until S6.
2. **S3** (24) and **S4** (20's second half). Both rest on S2's combination. S4 also rests on 53.
3. **S6** (26). After it, 53's distributed forms hold everywhere.
4. **S7** (35's second half, and C's open predicate, which needs the user's ruling first).

## Rulings ledger

1. **The knot is settled before the verification census** (the user, 2026-10-18). The census goes by layer from the base (verification-tree ruling 5), and `numerics` is the base.
2. **The raw-`__init__` rule is its own campaign, after this one** (the user, 2026-10-18; #611). This plan applies the rule wherever a `numerics` edge it untangles is an instance of the rule being broken.

3. **The test DAG's sub-cone structure is a statement of the code's ontology** (the user, 2026-10-18): "Notice that in broad strokes, the DAG of the tests, and more specifically, the way the sub-cones connect is really a statement of the code ontology (what machinery we have, what machinery depends on what other machinery and at what level and how they stack up)." So the concept order this plan writes is the same object as the verification cone's base. A knot in one is a knot in the other.
4. **F1: the axis comes after what defines it** (the user, 2026-10-18): "The Axis subcone has to come after what is required to define the axis, and an axis needs either a DiscreteMeasure (a nodal axis) or a Basis (a modal axis). The axis keeps it's generator for provenance, at least for now, because we have not decided yet on the proper way to implement forgetful maps concept." So the four typing edges into `axis` are the correct direction: axis is downstream of measure and basis. Each typing import is then an import of an upstream module, and the knot is that `measure`, `basis`, `frame` and `quadrature` import `axis` at runtime (the census's paths). The generator stays as provenance until a forgetful-map concept is chosen (no issue yet; candidate for one).
5. **F2: the operator module is cut, not the concepts** (the user, 2026-10-18, agreeing): the abstract `LinearOperator` is split from the composers, and the assembled and Green operators sit between them.
6. **F3 is an ontological question** (the user, 2026-10-18): "In principle, a frame connecting a discrete measure to a basis seems clean. If the space should have a frame and how it retracts or sections is a reasonable question (and if the space should have the frame of just extract the operators from the frame (this makes more sense, and seems to be essentially a marginalizing operator, so an axis carries the method to marginalize itself for example)." The user's leaning: the space does not hold the frame; it holds the operators extracted from it; the collapse pair is a marginalisation, and an axis carries the method to marginalise itself. Related: #525 (posing unit 2, "every collapse, projection and section is one verb of one frame hierarchy", plan `posing_u2_frame_collapse.md`).
7. **F4: `mass_matrix` is likely legacy; `basis.space` is likely the wrong arrow and a malformed concept** (the user, 2026-10-18): "the mass_matrix(measure) might be a legacy leftover from when we first developed it. basis.space might be the wrong arrow. Also, I think space typically is a composite of multiple axes, each one with either a DiscreteMeasure or a Basis, so this might be a malformed concept."
8. **F5: the quadrature part of `DiscreteMeasure` needs organising and cleaning** (the user, 2026-10-18): "this seems to indicate that there is some organizing and cleaning up to do in the Quadrature part of DiscreteMeasure."
9. **F6: test both readings, and think broadly** (the user, 2026-10-18): "Test both readings (and aditional perspectives, for example, does every DiscreteMeasure needs this invariance method? If this something related specifically to quadratures, so is this a justification to have AngularMeasure as a specialization of DiscreteMeasures? Think broadly."

10. **F3: the axis marginalises itself; the space maps to the space without it** (the user, 2026-10-18): "Did I understand well that your reading is "The axis knowing the method to marginalize itselt, and the space containing the axis being able to point to point/map to another space where the axis is marginilized" is the appropriate construction? It seems good to me. It also seems like a model axis should be allowed to marginalize (why would nodal be, and modal not?)". Ruled, with the orchestrator's precision:
    - The axis owns one covector: the integrals of its basis functions against its measure, `m_i = ∫ u_i dμ`, so that the marginal of a field `c` along it is `⟨m, c⟩`. On a nodal axis this is the weights (the all-ones vector in the metric). On a modal axis it is the moments of the basis (for Legendre, `e_0` scaled by the measure's mass). The attacker's "unit", the coordinates of the constant 1, is `G⁻¹ m` with `G` the axis's metric.
    - The space owns the lift into the product. The retraction is `R = ⟨m, ·⟩ ⊗ id_rest`, and the section is its Moore–Penrose right inverse `E = R^H (R R^H)⁻¹`.
    - Every axis whose basis functions are integrable against its measure marginalises, modal and nodal alike. The present modal refusal is retired, not kept as a boundary. An axis that cannot (an energy axis, whose collapse is the condensation's) refuses with its own reason.
11. **F4: retire `mass_matrix`; the axis receives its generator** (the user, 2026-10-18): "retire mass_matrix. Regarding giving basis a method to build an axis, is this not an inversion of giving axis a basis? Also, it there any advantage in basis creating an axis, instead of an axis receiving the basis (I suppose this is the instantiation of a modal axis, right?)". Ruled:
    - `mass_matrix` retires onto the frame's `discrete_gram` under `retirement-audit`. Its closed-form gates and their `verifies` and `catches` markers migrate to the Gram, along with the two docs `implements` directives. The abstract method and the two test stubs it forces go.
    - A modal axis is instantiated by receiving its basis, in its raw `__init__`. A convenience constructor lives on the axis side as a classmethod (#611's rule). The one advantage of `basis.space` (the basis knows which axis subclass, Legendre or harmonic) is kept by that classmethod on each axis kind. Leaving it on the basis is the backward `basis -> axis` runtime import the census found.
    - The same inversion on the nodal side: `DiscreteMeasure.axis()` (`measure.py:391`; 2 production callers, plus `Quadrature.axis` with 5) builds a nodal axis from the measure. By this ruling the nodal axis receives the measure, and the verb moves to the axis side. `[R]` the orchestrator's reading of the ruling, to confirm.
    - A space is built from its axes by `FunctionSpace.of_axes`. `for_basis`'s second constructor goes.
12. **F5: the quadrature part of the measure family is cleaned** (rulings 8 and 11): the two rule lists become one, the second names of one verb become one name, and the production-dead factories and the never-used re-exports retire. Each retirement goes under `retirement-audit`.
13. **F6: the support carries the group action** (the user, 2026-10-18: "If you see this as the correct path forward, I endorse it."). The manifold carries its group's action (a G-space). The measure contributes one group-free verb, equality up to a relabelling of its points (#403's equality). Invariance is that equality after the pushforward along the action. A support with no action refuses the question rather than zero-padding an answer. The 2026-08-29 ruling stands: no `DiscreteMeasure` subtypes.

14. **The code is never the oracle; (a) and (b) are investigated** (the user, 2026-10-18): "(a) and (b) I'm genuinely not sure. It seems like we need to investigate, gather more data, reason about it from multiple perspectives and adversarially attack perspectives to find the correct way to think about it. The important part is that we should never take the current state of code as a correct oracle. It is by definition under development, which means we should always be willing to refine our understanding and change the code. The development is guide by math, ontology (what machinery, how they connect, direction of connection, and the stackup) and elegance." A reading of the tree establishes what IS; what SHOULD be is argued from the mathematics, the ontology and elegance, and the tree changes to match.
15. **(c): the curvilinear LD refusal is machinery not yet built** (the user, 2026-10-18): "Is it the right way for the LD to behave as it is? ... Is the reason why it doesn't work for curvilinear a side-effect of the modal description of it being under-developed? ... a refusal here is just something not yet build. there is no principled reason that I can think of for why we couldn't analyse a cylinder or a sphere with LD scheme." Filed as #612. The user's attack is one of the investigation's questions (below).

16. **(a) A frame is the sampling map of a basis on a measure, with its dual** (the user, 2026-10-18): "I think that is a fair reading, because we choose the Discrete Measure and Basis, which fundamentally chooses the representation, but how fine they are is also a resolution aspect." The choice of the measure and the basis fixes the representation, and how fine each is fixes the resolution. One construction serves both: the angular frame between representations, the homogenisation and condensation frames between resolutions.
17. **(b) Every axis carries the pullback of its factor's own measure** (the user, 2026-10-18: "I adopt the pulled back carry."). A nodal axis carries its weights. A modal axis carries its basis's Gram under the factor's measure: the discrete Gram for SN's angle, the continuum mass for LD's cells.
    - The user then asked: "Regarding the orders a quadrature can or cannot integrate, is this a tightness statement?" Yes, `[R]` the orchestrator, in `frame.rst`'s own definition (a frame is tight with constant c when its frame operator is c·I, in named metrics).
    - A frame is tight iff its canonical dual is a scalar multiple of the frame. Production's reconstruction is the canonical dual for the continuum measure (the addition-theorem factor). So the square `M.H = R/W` closing with one scalar is the tightness of the sampling frame, with the quadrature weights on the ordinates and the continuum Gram on the coefficients.
    - Its condition: the rule integrates every product of two live basis functions exactly. For harmonics up to L, the rule is a weighted 2L-design on them.
    - Three regimes, `[M]` `frame_bridge/attack_probes/p6_tree_square.py` and `scratch/adjoint_underresolved/`:
      - tight (34 of 37 probed rows);
      - a frame but not tight (injective, dual not proportional: LS4 L=3, 0.067; LS8 L=5, 0.039);
      - not a frame (more live modes than the rule separates, lower bound zero: LS4 L=5, 66 modes on 24 ordinates).
    - Under the pullback rule the frame is Parseval in its own pulled-back metric always, by construction. So "the rule integrates order L" is tightness against the continuum metric: whether the discrete measure reproduces the continuum one on the products.
    - The remaining method choice: standard S_N (the continuum dual; exact in the tight regime, correct and self-adjoint but not tight elsewhere, today's code) against Galerkin quadrature (a dual that closes every identity in all three regimes, a different method). Ruled (the user, 2026-10-18): "We keep Sn for now and file Galerkin as an issue (because it will require us to investigate that deeper)." Standard S_N stays the method; Galerkin quadrature is #613, an investigation before any build.
18. **(c) Curvilinear LD gets a Gauss-nodal spatial axis** (the user, 2026-10-18): "The nodal seems correct for today. Maybe the modal one has a different use (finite element?) which is not currently in scope." Within each cell, the Lagrange basis at the 2-point Gauss rule of the weight `r^{max(d−2,0)}`: one nodal spatial axis from a discrete measure, with both curvilinear LD forms exactly diagonal. The modal within-cell axis is not retired by this ruling: it may serve a finite-element method outside today's scope. #612 carries the ruling.

19. **Order the concepts by what is joined to what; the cycles were the signal** (the user, 2026-10-18, on step 3's residual edges):
    - On the operator split: "Option 1 seems correct, but before that simply making a ruling on it, try to check how the ontology stacks up. For example, Riesz operators arise when operator and what other concept get joined? In other words, don't simply rely on a ruling. Instead, try to make a minimal organization of the ontology by sub-cones to see exactly where the things join. Like we did with frame."
    - On the Green operator: "does all inverse of an operator sum returns a GreenOperator mathematically? I think the inverse of the sum of the streaming and collision operators ends up being a GreenOperator. Is that the only GreenOperator we can make, or other sums also result in Green operator mathematically?"
    - On the lift: "An operator maps a space onto another space (this is implicit in that an operator has domain and co-domain). Therefore, the operator sub-cone must come after space for the ontology to be well-formed. This is a strong indication that option 1 is at least directionally correct from principled thinking, and trying to break the cycles was the signal."
    - On the reference measure: "Option 1 looks like it takes my other answers directionally. You might want to rethink this anyway and see if the ontology is well-formed with option 1. The question is "what are we joining and when?""
    - So the order is derived, not ruled edge by edge. Candidate E below is that derivation.

20. **The inverse is formal, and a sum's is realised by the Strategy's splitting; the name goes to the stronger match** (the user, 2026-10-18): "Option 1 seems good. InverseOperator being taken is not a reason… it depends on why inverse operator is taken, and who holds the stronger match to the inverse operator concept."
    - `[M]` by reading `operator.py:2697`: `InverseOperator` is the generic member of the #226 inverse family, "the name is earned by exactly the universal contract and nothing more". Its `apply` delegates to the forward's own `solve`, and the round trip holds. Its siblings earn their names by a distinguishing invariant: `SweepOperator` (a triangular, seed-independent direct solve), `MatrixInverseOperator` (a dense factorisation), `GreenOperator` (the Neumann series of a splitting).
    - So `InverseOperator` IS the formal inverse, and the stronger match. What a sum lacks is a `solve` of its own, which is why `OperatorSum.inverse` reaches for `GreenOperator`.
    - The orchestrator's reading `[R]`: a sum's inverse is `InverseOperator(sum)`, realised by a solve the Strategy supplies (its splitting `A = M − N` and a driver). The Neumann property becomes a property of that realisation, not a type chosen by the sum's spelling order. The fused `L + C` inverse stays structural (`SweepOperator`). No new name is needed.
21. **Quadrature is a discrete measure with a fibration** (the user, 2026-10-18, adopting the recommendation): a discrete measure whose levels are the fibres of the quotient map. It sits before the axis; its cached angular frames are a cache on the frame side, not part of the concept.

22. **The Gram datum is the frame, placed at the join of measure and basis** (the user, 2026-10-18: "Yes, that looks correct. The question is "is there a place that this should fall onto, or should this be its own place. (I.e. is this concept in the wrong cone or is this a cone of its own)."). The orchestrator's answer `[R]`:
    - Not a cone of its own. It is a real join: both a basis and a measure are needed, and neither is an ancestor of the other. So it cannot fall onto either parent's cone.
    - It falls onto the frame. `FrameBase` is "a Basis bound to a DiscreteMeasure", and its Gram is the metric that binding induces; with a continuous measure it is the continuum pairing, the same concept with exact integrals.
    - E's error was placing the whole frame after the operators. The frame sits where the verification tree's ruling 25 put it, at the join of measure and basis, upstream of the axis.
    - Its faces (analysis and reconstruction) and the Riesz operators are not joins: the frame and the space are ancestors of the operator, so they are chain steps after it (the antichain test).
    - It also answers F1's open question. An axis is a **forgetful image**: a nodal axis keeps a measure's weights and forgets its nodes; a modal axis keeps a frame's coefficient index set and Gram, and forgets the nodes and the basis functions. `Axis.generator` points at the domain of that forgetful map, which is why it looked like provenance pointing downstream (ruling 4's "until a forgetful-map concept is chosen").
23. **A sum is a formal linear combination of operator terms** (the user, 2026-10-18: "Agreed. Your reading is endorsed."): an unordered collection of signed, scaled terms, so the abelian-group laws (associativity, commutativity) and scalar distributivity hold as equality. A product stays an ordered sequence, flattened, since composition is associative but not commutative.

24. **Assembly is a homomorphism, moved out of the algebra** (the user, 2026-10-18: "I agree with your recommendations. Let's go ahead."). Assembly maps the operator algebra into the matrix algebra; with a sum a formal linear combination it is a fold over the terms, `[Σ cᵢTᵢ] = Σ cᵢ[Tᵢ]`, downstream of the algebra; each leaf returns its raw matrix.
25. **The adjoint's placement waits on whether the adjoint is always posed** (the user, 2026-10-18, interrupting before the adjoint ruling): "I just realized there is one question we need before ruling the adjoint machinery as part of the algebra. I think this might depend on if the adjoint is always well-formed and can always be posed. Because if the answer is "no" then that means if the machinery is all together, we might need guards of various types in code, and this seems like the place that encodes this question "is the adjoint always available"."
    - The orchestrator's mathematical split `[R]`. In finite dimensions, a LINEAR map between spaces with NON-DEGENERATE metrics always has a unique adjoint (Riesz representation). It fails in four ways:
      - (1) a space with no metric: only the metric-free transpose exists;
      - (2) a degenerate metric: the pseudo-inverse "adjoint" is a true adjoint only if A maps ker G_V into ker G_W. This case is live here, since the dense Grams are singular;
      - (3) an affine, not linear, map: only its linear part has an adjoint (ERR-075's class);
      - (4) no transpose action implemented: the adjoint exists, but cannot be realised.
    - Cases 1 to 3 are questions of posing, which types can decide. Case 4 is realisation, a capability, as with the inverse (ruling 20).
    - `[M]` the orchestrator: 52 classes define `is_adjointable`. `.H` and `.dual()` both refuse when it is false (`operator.py:980`, `:1027`), so today one flag answers "is a transpose coded" and gates even the metric-free transpose.
    - A census of each false flag's reason is dispatched (`scratch/numerics_untangling/adjoint_posed/`).

    - **The census** (`[M]` the explorer, `scratch/numerics_untangling/adjoint_posed/report.md`; saved by the orchestrator). The adjoint is not always posed, and every failure mode gets past `is_adjointable`, which means only "a transpose is coded":
      - case (3), affine maps: already solved by type. Every affine law splits into a linear factor and a typed source channel, and 0 affine operators remain.
      - case (1), a missing or wrong metric: `RankOneOperator` has `None` ends, its flag reads `True`, and its `.H` raises. `DualSpace` carries the primal's G instead of G⁻¹, so the Riesz legs' `.H` raises `TypeError` while their flag reads `True`.
      - case (2), a degenerate metric: live. On `folded_product(4,4)` at L = 4, the frame's reconstruction and the moment transfer alone map the domain metric's kernel outside the codomain's, so their `.H` is a true adjoint only on the range: reciprocity over all inputs misses by 1.00 and 0.495, and by ≤ 1.6e-14 on the range. Composites are exact, and no production `.H` daggers those two alone.
      - case (4), no transpose coded: the flag lies both ways. `SweepOperator` reads `True` on 2-D poses while `S.H.apply` raises, so `solve_sn_adjoint` in 2-D fails inside the solve. `LeakageOperator` and `DiffusionBoundaryOperator` read `False`, while their assembled matrix transposes. `_transpose_substitution_ready` tests `callable(...)`, which a raising method passes.
    - The orchestrator's reading `[R]`: with cases 1 and 2 made structural, every linear operator between non-degenerate metric spaces has an adjoint, so the operators form a dagger category and the adjoint term belongs with the algebra. Posing then needs no guard, and one realisation capability remains (a transpose action), shared by the dual and the adjoint, as with the inverse (ruling 20).
26. **Degenerate metrics are replaced by the live subspace, so the adjoint is exact everywhere** (the user, 2026-10-18: "replace to get exact adjoint everywhere."). A space whose metric is singular is replaced by its live subspace, the metric's range, on which the metric is non-degenerate. The pseudo-inverse doctrine then has no range restriction left to carry. With case (1) closed by type (every operator end a space; the dual carrying G⁻¹) and case (3) already closed, the adjoint is always posed, which settles ruling 25: the adjoint, the dual and the Riesz operators stay with the algebra as dagger-category structure.
27. **The lying adjointability flags are fixed now** (the user, 2026-10-18: "Let's fix it."). The fix is first pinned by one honesty gate: for every operator the registries can build, `is_adjointable` is true exactly when `.H` constructs and applies. The gate is red today on every lying case. The `PrescribedInflow` docstring's false example is corrected with it.

**Answered from the same exchange, now ruled (22 and 23):**
- *"I think the Riesz operators also need Gram no? Where do they get that from?"* Yes. ♭ maps v to Gv, with G the space's metric. The space's metric is assembled from its axes', and each axis carries its factor's pulled-back measure (ruling 17): weights on a nodal axis, a Gram on a modal one. So the chain is: Gram datum (basis × measure), then axis metric, then space metric, then Riesz, the frame's faces and the adjoint. This is the case for the split: the Gram datum becomes the one source for every downstream metric. Today the same quantity is computed in four places: the frame's `discrete_gram`, the retiring `mass_matrix`, the LD scheme's moment mass, and the volume measure's cell weights.
- *"That's basically asking of a sum of operators should be associative, right?"* Both associativity and commutativity, and scalar distributivity. L(V, W) is a vector space under addition, so sums obey the abelian-group laws and distribute over scalars. Today's `OperatorSum` honours none of them as equality: it is a binary tree with no `__eq__`, so `(A + B) + C` and `A + (B + C)` are different trees, and `A + B` differs from `B + A`. A sum as an unordered collection of signed terms (a formal linear combination) is the normal form that honours all three. Composition is associative but not commutative, so a product stays an ordered sequence, flattened.

## Measured state

All `[M]` 2026-10-18, the test-architect's scripts under `scratch/verification_tree/subcones/`. The tree was at `main` `c5c21b25`, and `orpheus/` was unchanged during the run.
- The module graph is built from each module's imports. A module's **top-level** imports form an acyclic graph: 44 modules, 91 edges, depth 7, with 8 joins. `frame <- basis, measure, projection, space` is one of the joins.
- **Function-scope (lazy) imports** add 19 edges (`lazy_edges.py`). 7 of them close a cycle with the top-level graph. In each line below, the module after the arrow imports the module before it inside a function, while the module before it already imports the one after it at top level:
  - `assembled_operator -> operator`
  - `frame -> space`: verified by the orchestrator. `FunctionSpace._axis_collapse_pair` imports `frame._collapse_pair` (`orpheus/numerics/space.py:558`). Nexus `native_place` independently places `_collapse_pair` on `FunctionSpace`.
  - `generating_measure -> manifold`
  - `generating_measure -> measure`
  - `green_operator -> operator`
  - `spaces -> basis`
  - `traced_memo -> _traced_memo_boot`
- **`TYPE_CHECKING` imports** add 13 edges beyond every runtime import (`tc_edges.py`). 10 of them close a cycle with the runtime graph. In each line, the module after the arrow names the module before it for typing only, while the module before it depends on the one after it at runtime:
  - `basis -> axis`, `frame -> axis`, `measure -> axis`, `quadrature -> axis`
  - `measure -> basis`
  - `measure -> invariance`
  - `exactness -> manifold`
  - `quadrature -> spaces`
  - `posing -> iteration`
  - `space -> assembled_operator`
- The cycles span 21 modules (runtime plus lazy) and 14 modules (runtime plus typing). The concepts involved are axis, basis, measure, frame, quadrature, the function spaces, symmetry and invariance.
- **Not yet measured:**
  - what each of the 17 edges is used for (which symbol, in which def, for what);
  - whether a lazy import is deferred for import order alone or hides a real dependence;
  - the edges inside one module, between concepts that share it.

## First readings of the knot (`[R]` by the orchestrator, from the edge list only; each is a question, not a finding)

- **The axis names four of its downstream types** (basis, frame, measure, quadrature) for typing. CLAUDE.md says "on an axis: a discrete measure, a basis and the frame pairing them". So either the axis *holds* those objects, or it *is* a coordinate interval that they are defined over. If it holds them, `Axis` welds two concepts: the coordinate (upstream of everything) and the discretisation laid on it (downstream of measure and basis).
- **The measure and the basis each refer to the other.** The measure depends on the basis at runtime, and the basis names the measure for typing. Under ruling 25 they are parallel sub-cones that meet only at the frame. A runtime dependence of measure on basis is then either a leak or a concept the frame should own.
- **The space reaches its frame** (`frame -> space`, lazy). The frame is minted on a space, and the space caches a pair the frame mints. The stage-2 generator rule (memory: "Frame/Scheme induce space + operator together at one site") says the frame induces the space's data. So the cache may belong to the frame's output and not to the space.
- **The generating measure is reached lazily from the measure and from the manifold.** The generating measure is downstream of both, so this is either a factory on the wrong type or a classmethod that should sit on the generating measure (#611's rule).
- **`traced_memo -> _traced_memo_boot`** is the memo's own boot split, probably deliberate (a boot module must import before the memo). It needs its reason checked, not a redesign.

## The edge-use census (step 1, landed 2026-10-18)

`[M]` the explorer, `scratch/numerics_untangling/edge_uses.md`; scripts beside it. The tree is `main` `c5a4c15d` (`orpheus/` identical to the capture's `c5c21b25`). Acceptance was checked by the explorer and re-read by the orchestrator: 17 rows, with the edge set equal to the scripts' output; 26 import sites, each an AST import statement; 101 cited test ids, each a context in `cov_numerics.sqlite`. The positive control `space.py:558` was found first.

**The import graph does not say which imports can move.** A probe executed each site's import at module top, in memory, in fresh interpreters importing either module first. 6 of the 7 lazy imports fail with ImportError either way; the boot split works only with `traced_memo` imported first. Of the 10 typing imports, only 2 load at runtime in both orders: `posing -> iteration` and `space -> assembled_operator`.

**The 17 edges, grouped by what the importing site does** (the orchestrator's grouping, read from the table):

| family | edges | what the sites do (`[M]`, the table) | reading (`[R]`, the orchestrator; for the user) |
|---|---|---|---|
| F1. The axis keeps its generator | 8, 10, 11, 15 | One field: `Axis.generator: DiscreteMeasure \| Basis \| Quadrature \| FrameBase \| None` (`axis.py:186`), excluded from identity (`compare=False`). It is written at 6 sites, and `Axis.generator_as(kind)` returns it after an `isinstance` check. All 5 production `generator_as` calls narrow to `Quadrature` (`harmonic_frame.py:451`, `:493`; `sn/mesh/reduced_operator.py:522`; `sn/angular/closure.py:1675`, `:2236`). One bare read, `spaces/moment_head.py:118`, branches on `GalerkinFrame` and `TruncatedBasis` and calls `.at_order(L_new)`. | CLAUDE.md: the stage-2 generators are "read for their induced data and then forgotten". The axis does not forget its generator, so the upstream object carries a reference to four downstream ones. The 5 readers need the quadrature back from the axis. Either what they need is induced data the axis should carry by name, or the quadrature, and not the axis, is what they should have been handed. |
| F2. Base-class module builds its subclass | 1, 5, 17 | `OperatorSum`, `OperatorProduct` and `ScaledOperator.assemble` in `operator.py` construct `SparseAssembledOperator`, a subclass of `LinearOperator`, which is defined in `operator.py`. `OperatorSum.inverse` constructs `GreenOperator(self)`, a class parameterised by `OperatorSum`. `assembled_operator.py` names `FunctionSpace` for typing only, and that import loads in both orders. | The algebra really is closed: the assembly of a sum is an assembled operator, and the inverse of a sum is a Green operator. So the dependence is mutual by definition. `[R]` The knot is in the module cut: `operator.py` holds both the abstract `LinearOperator` and the composers. Splitting the abstract type from the composers puts the assembled and Green operators between them, and no concept changes. |
| F3. The space calls its frame | 2 | `FunctionSpace._axis_collapse_pair` calls `frame._collapse_pair(self, label)` and caches the (retraction, section) pair on the space (`space.py:558-560`). Its callers are `FunctionSpace.retraction` and `.section`. Nexus `native_place` places `_collapse_pair` on `FunctionSpace`. | The pair is the frame's to mint ("frame-induced at ONE site", the docstring). The stage-2 rule says the frame induces the space and its operators together. Either the frame mints the pair when it mints the space, or retraction and section are the space's own operations and `native_place` is right. |
| F4. The basis does the frame's work, or reaches up to its space | 6, 12 | 5 bases each define `mass_matrix(measure)`, which reads `measure.nodes` and `.weights` and returns `∫ u_i u_j` by quadrature. It has 0 production references, and `base.py` calls it "a one-off diagnostic". `LegendreBasis.space` and `SphericalHarmonicBasis.space` call `<Space>.for_basis(self)`, a classmethod of the downstream space. | A mass matrix pairs a basis with a measure, which is the frame's Gram (ruling 25's join), so it sits on the wrong concept. It is also dead in production. `basis.space` is a reverse accessor to a classmethod that already exists, so it is a convenience pointing downstream (#611's rule, read in reverse). |
| F5. The measure reaches the generating measure | 3, 4, 9 | `measure.gauss_legendre(n)` and `gauss_chebyshev(n)` call `LEGENDRE.gauss(n)` and `CHEBYSHEV_T.gauss(n)`: 0 production call sites, re-exported from `orpheus/numerics/__init__.py`, called by tests only. `DiscreteMeasure.generating_measure` has 0 references. `manifold._sphere_mod_o2` stores `LEGENDRE` in `Quotient.reference`, typed by the protocol `ReferenceMeasure`; it is read in `quadrature/registry.py:1043`. | A generating measure (a weight function whose Gauss rules are discrete measures) is upstream of the rules it generates. The measure module's factories and property are then classmethods or methods that belong on the generating measure, and 3 of the 4 sites have no production consumer. The quotient's reference measure is a real use, already typed by a protocol, so it may be an artefact of the protocol's home module. |
| F6. Functions annotated with a downstream type | 13, 16 | 9 `invariance` functions take `measure: DiscreteMeasure`, read only `.nodes`, `.support` and `.weights`, and are called from `DiscreteMeasure`'s own verbs (`measure.py:1082-1123`) and from quadrature. `native_place` puts 4 of them on `DiscreteMeasure`. 3 defs in `spaces/angular_trace_space.py` take a `Quadrature` and read `.N`, `.nodes`, `.weights` and `.mean_axis_cosine`. | Invariance is a question asked OF a measure ("invariance is the measure's question", CLAUDE.md's vocabulary), so these are the measure's methods living in a kernel module, or the kernel is a lower concept that needs only a weighted node set (a protocol). The angular trace space consumes a quadrature as its data, and that is a real dependence on the downstream quadrature. |
| F7. Deliberate or harmless | 7, 14 | `_traced_memo_boot` imports `traced_memo` as the last step of a child process's ordered boot. `KEigenvalue` stores an `EigenPosing` and reads `.pencil.rhs` and `.rayleigh`; the import loads in both orders. | The boot split is the boot's design, kept with its reason. `posing -> iteration` can be a runtime import, and it is an artefact. |

**Production-dead defs found on the way** (`[M]` AST over `orpheus/` by name, tests excluded): `DiscreteMeasure.generating_measure` (0 references); `measure.gauss_legendre` and `gauss_chebyshev` (0 calls; the 8 bare `gauss_legendre(...)` calls in `orpheus/derivations/` bind a different function of the same name); `mass_matrix` on 5 bases (0 references). Each is a retirement candidate under `retirement-audit`. Each also needs the `coding-elegance` Pattern 6 check: is it zero-consumer AND self-conceding ("a one-off diagnostic"), or is it a fuller-view oracle that a permanent test consumes?

**The census's own limits:** coverage covers `tests.gates.numerics` only; `_boot` runs in a child process the capture does not measure; Nexus `native_place` scores only module-level functions and its replies stopped at 22 of 50 rows; method and property callers are found by name, without resolving the receiver.

## The F3 to F6 investigations (landed 2026-10-18)

Reports: `scratch/numerics_untangling/f4/report.md`, `f5/report.md`, `f6/report.md` (explorers, facts only), and `ontology_f3_f6.md` (the cross-domain attacker), each with its scripts beside it. Each was read by the orchestrator with session context. The explorers' reports were saved by the orchestrator, because the harness refuses a sub-agent's write of a report file.

**F3, the collapse pair.**
- `[M]` (the attacker; `ontology_f3_f6/probe_f3.py`): no frame mints a collapsed space. All 6 production collapses are along `"angular"`, on products built by `FunctionSpace.of_axes`.
- The pair is integration along one factor of a product space (the pushforward on densities), with its Moore–Penrose right inverse `E = R^H (R R^H)^-1` as the section. `R R^H` is the measure's total mass `Σw`, so ERR-051's divisor is derived, not chosen. `E∘R` is the conditional expectation onto fields that are constant along the axis.
- The pair needs one datum per axis: the coordinates of the constant function 1 on that axis. The attacker calls it the axis's **unit**, an invented name. Together with the axis's metric it builds the pair.
  - On a nodal axis the unit is the all-ones vector: `[M]` it matches today's pair to ≤ 9e-16.
  - On a modal Legendre axis it is `e_0`: `[M]` it matches the nodal result to 0.0.
- So a modal marginal exists. The docstring that refuses one names the wrong covector, which `[M]` misses by 9.1 %.
- The frame's `_collapse_pair` reads back exactly what it is given: weights equal the axis weights, and the Gram entry equals `Σw` to 0 to 6 ULP. Its proposed home is a gate oracle, not the mint.
- **Proposal** `[HYPOTHESIS]`: the axis owns its unit; the space owns the lift into the product (`R = ⟨unit, ·⟩ ⊗ id_rest`, `E` its pseudo-inverse). This fits ruling 6's leaning ("an axis carries the method to marginalise itself") for the one-axis datum. The whole operator needs the product, which only the space has. It agrees with #525's mathematics (its scope item 9), but not with its title, since no frame object owns the pair.
- **Open:** whether `E.apply` stays bit-identical. That needs `Σw` summed in the Gram's order, or the gram-derivation gate re-pinned at its nulp tier.

**F4, `mass_matrix` and `basis.space`** (`[M]` the explorer).
- `mass_matrix` meets both of `coding-elegance` Pattern 6's trim conditions:
  - 0 production calls, now and in its history;
  - docstrings conceding it ("a one-off diagnostic", "Not on any hot path", "no consumer");
  - bit-equal to the frame's `discrete_gram` on 11 of 11 fixtures (one formula under two spellings, X4), and no test compares the two.
  - It has 6 production definitions, the sixth in `sn/` (`LossKernelBasis`).
- What still points at it, and so must migrate (`retirement-audit`):
  - closed-form test gates with `verifies` and `catches` markers (ERR-039, `real-sh-discrete-orthogonality`);
  - `implements` directives in `frame.rst:516-536` and `spherical_harmonics.rst:328`;
  - the `@abstractmethod`, which forces two test stubs.
- Its history: it landed as `discrete_mass_matrix` (`0eb9cf39`, 2026-05-26) and was made abstract in the Frame carve (`c65154a0`, 2026-06-23). The frame's `discrete_gram` came two months later (`0317373d`). The user's "legacy leftover" (ruling 7) is consistent with this history.
- A stale role the dead-reference tool misses: `spherical_harmonic_basis.py:153` and `:198` still say `:meth:discrete_mass_matrix`, a name removed in `28d92758`.
- `basis.space` is not dead. It has production consumers on every frame family (`FrameBase.basis_space`, `test_space`, `moment_head.truncated_head`).
- For Legendre and spherical harmonics, `for_basis` mints a single-axis space, with the basis as that axis's generator. It calls the raw constructor and bypasses `FunctionSpace.of_axes`, "the ONE composition mechanism". Its result is `==` and hash-equal to `of_axes` on its own axis; only the subclass and the name differ.
- Indicator and loss-kernel bases return a plain `FunctionSpace(name, shape)` with no axes.
- So ruling 7's "malformed concept" reads `[R]`: a basis builds one axis, and a space is a product of axes (`of_axes`). `basis.space` is that one-axis product under a second constructor.

**F5, the quadrature machinery** (`[M]` the explorer).
- `Quadrature` holds a `DiscreteMeasure`; it is not a subclass. The two differ in identity: the measure is unhashable and raises on `==`, while the quadrature has content equality.
- Every numerics Gauss–Legendre rule comes from `LEGENDRE.gauss(n)` by delegation: a pure alias, one wrapper that adds σ_x, one that adds the polar orbit, and the `Quadrature` classmethod and the registry entry on top. They are bit-equal.
- The derivations code computes GL separately (numpy, mpmath, its own Golub–Welsch), within 1 ulp in nodes and up to about 1e4 ulp in weights.
- Production solver code never constructs a rule: all 13 production `Quadrature` factory calls are in `derivations/continuous/mms/sn.py`.
- The exactness machinery decides nothing in production: its two deciding sites sit behind `select_quadrature`, which has 0 production callers.
- Cleaning candidates:
  - two lists of the same rules (`Quadrature`'s 5 classmethods and the registry's 4 factories; the registry builds a bare `DiscreteMeasure`, and `folded_product` is unregistered);
  - pass-throughs under second names (`n_points` against `n_ordinates`/`N`; `permutation_under` against `ordinate_permutation`, with two constants of the same value 1e-13);
  - 8 re-exported rule names never imported through the package;
  - 3 production-dead factories (`measure.gauss_legendre`, `gauss_chebyshev`, `equispaced`).
- The stale theory-page sentence and code comment it found are fixed [REMEDIED 2026-10-18 @cf1ec9bc].

**F6, invariance** (`[M]` the explorer's spy over `numerics`, `geometry`, `transport` and `sn`, serial).
- 0 of 20 168 production-entered invariance calls are on a non-angular measure. Every one is on S² or a quotient of it.
- The other production invariance question goes around the kernel: `rules_product._derived_product_group` asks the circle factor and the μ factor through `RigidMotion.preserves`, at their own dimension.
- **A latent wrong answer.** On a support that is not a quotient, the kernel zero-pads the nodes into R³ and answers (`invariance.py:142-160`):
  - an energy-group measure reports σ_y, σ_z, SO2_x and O2_x invariance True;
  - `equispaced(0,1,4)` and `equispaced(-1,1,4)` disagree on σ_x, since the answer depends on the origin of an embedding the support never names;
  - the chart μ-rule reports SO2_x True and SO2_y False, the "x" being the padding's column 0.
  - Production never asks these questions, so no production answer is wrong today.
- The group acts on the support (`induced_action` on a quotient), not on the measure. A primitive over a bare weighted point set already exists (`permutation_preserving`, `transformation.py:334`).
- **A prior ruling bears on the user's `AngularMeasure` question:** "RULED 2026-08-29 … `DiscreteMeasure` does NOT branch into Spatial/Angular/Energy subtypes … `Quadrature` IS the angular measure-plus-fibration" (`space_and_kernel_binding_campaign.md:2248`). The name `AngularMeasure` is taken by a Protocol in `sn/angular/redistribution.py:224`, which `Quadrature` implements.
- **Proposal** `[HYPOTHESIS]` (the attacker):
  - invariance is a property of the pair (measure, a group's action on its support); the support carries the action (a *G-space*: a manifold with an action of a subgroup of O(3));
  - the measure contributes one group-free verb, equality up to a relabelling of its points (the equality #403 owes);
  - invariance is that equality after the pushforward along `support.action(g)`;
  - `[M]` the composition agrees with today's kernel on 75 of 75 rule × finite-group rows, 45 of them invariant, and a perturbed weight breaks both;
  - a support with no action then refuses the question instead of answering it;
  - no `AngularMeasure`: the chart μ-rule is "angular" with no action, and a spatial R² point set asks a meaningful D_4h question. If anything is specialised, it is the manifold;
  - the concept order becomes groups, then the G-space manifold, then the measure, then invariance as a composition; census edge 13 vanishes.

**Side finding, under diagnosis.** Two traced-memo gates failed when many test processes shared the host: `test_m3_15_two_processes_racing_on_one_miss_leave_one_valid_entry` (the entry read `Absent`, where `Hit` was expected) and `test_b13_a_miss_indexes_only_bytes_it_has_not_indexed`. Both pass serially: `[M]` 56 of 56. The numerics-investigator is reproducing both under controlled load (`scratch/memo_race/`). The first may be a real race in the memo's write.

## What an axis needs to be constructed (the user's question, 2026-10-18: "Evaluate it. Check if an axis only needs a measure or a basis to be constructed.")

`[M]` by reading, the orchestrator, 2026-10-18, at `main` `e6354b31`: `orpheus/numerics/axis.py:117-466`; `spaces/spherical_harmonic_space.py:193-215` (`for_basis`); `frame.py`, `_dressed` (the frame's re-dressing of a coefficient space); `measure.py:391-420` (`DiscreteMeasure.axis`).

**An axis's identity** is its class plus `(label, shape, weights, kind)`, plus each family's own fields: `EnergyAxis.edges`, `LegendreAxis.spent_axis`. `generator` is provenance only (`compare=False`). Where each part comes from:

| part | nodal axis (from a measure) | modal axis (from a basis) | energy axis |
|---|---|---|---|
| `label` | the role in the posing (`"angular"`, `"spatial"`); not in the measure | the role (`"harmonic"`); not in the basis | `"energy"`, fixed |
| `shape` | `(n_points,)` | the basis's coefficient layout (`(L+1,)`, `(L+1, 2L+1)`) | `(n_groups,)` |
| `weights`, the metric [the energy column below REFUTED 2026-10-18, see the exercise] | the measure's weights | **a Gram, which the basis alone does not fix.** Minted from the basis, it is the continuum Gram `basis.metric_per_ell`, the basis against its continuous reference measure. Re-dressed by a frame, it is the Parseval inverse of the discrete Gram over a quadrature. When that Gram is dense, the axis's `weights` are emptied and the metric is positioned on the space (`FactoredMetric`), because the axis can hold only a diagonal | none: the counting measure, as a theorem (group integrals pair with group averages without widths) |
| `kind` | NODAL, implied by the measure | MODAL | NODAL, enforced |
| family content | — | `spent_axis`, from the basis's domain `S²/O(2)_a` | `edges`, from the energy grid |

So:
- **A nodal axis needs a measure, plus its role label.** Everything else is the measure's.
- **A modal axis needs a basis AND a measure.** The basis gives the shape and the kind. The metric needs a measure to integrate against: the continuum reference measure, or a quadrature. The two give different axes, since the weights are in the identity. The tree carries one moment space per carrier, the frame's (ruling R-6.2c-1, cited in `HarmonicAxis`'s docstring). That is why `FrameBase.basis_space` re-dresses the basis-minted space, and why every frame family reads `basis.space` only to replace its metric (F4).
- **An energy axis needs neither:** a partition of energy (the group edges), with its counting measure and the group-indicator basis implicit.
- **Ruling 10's marginal covector** `m_i = ∫ u_i dμ` also needs both the basis and the measure.

**The reading this suggests** `[R]`, the orchestrator, for the user: **every axis is what a frame induces**, a measure paired with a basis.
- A nodal axis is the degenerate frame: a measure with its own point-evaluation basis, whose Gram is `diag(w)`. `DiscreteMeasure.axis` mints that frame's axis without naming the frame.
- A modal axis is a Galerkin frame's axis. Its metric is that frame's Gram. The basis-only mint is the frame of the basis with its continuum reference measure, a second frame the tree does not name.
- An energy axis is the frame of a partition with its indicator basis under the counting measure.

That is CLAUDE.md's vocabulary taken literally ("the frame pairing them, which induces the space's metric"), and the stage-2 generator rule ("Frame/Scheme induce space + operator together at one site; read for their induced data and then forgotten"). Its consequences, each to rule:
1. **The concept order.** `Axis` and `FunctionSpace` are value types with raw constructors (label, shape, metric, kind), and import nothing above them. The frame, the join of measure and basis, mints them as their generator. So the order is: measure and basis in parallel, then the frame, which mints axes and spaces. The types themselves sit below the frame. This refines ruling 4: an axis comes after what defines it, and what defines a modal axis is a frame, not a basis.
2. **F4, refined.** A modal axis receives its metric from a frame. So the convenience constructor is the frame's (`FrameBase.basis_space` today), not the basis's (`basis.space`). The continuum mint becomes a named frame (the basis with its reference measure), or it retires if nothing but tests needs a continuum-metric head.
3. **The axis's metric is a diagonal slot.** A general modal metric (a dense Gram) cannot sit on the axis today; it is positioned on the space instead. Under "the metric is an object" (CLAUDE.md), the axis would carry a metric object, diagonal or dense.

### The exercise: what induces an axis, grounded in the spatial schemes (2026-10-18)

**The user's framing** (verbatim): "I think technically every axis does end up needing either at least a discrete measure or a basis ... even the energy axis, which looks like it doesn't need any, actually gets an energy grid, which itself it a discrete measure ... When we're generating the spatial axis, what induces the axis is the discretization scheme that we choose for the discretization, for example, diamond difference induces a nodal axis and linear discontinuous induces a modal axis. In other words, nodal or modal is an axis representation that is induced by the discretization choice for that axis. Now you can sharpen that, refute that, improve that ... when we pick a linear discontinuous spatial discretization and the spatial axis is modal, why would it need also a discrete measure? or it wouldn't? Does it require a frame? Why? and do both require a frame, even the diamond difference one? why?"

**Facts read** (`[M]` the orchestrator, `main` `37382ee2`):
- `EnergyGrid.as_measure()` (`data/energy_grid.py:194`) returns a `DiscreteMeasure`: nodes are the group indices, weights all 1 (the counting measure), and the support is `EnergyGroups(n)`. The grid also has `as_basis()`, the partition. Its own section header calls the pair "The dual frame views (source = measure, target = basis)". So the user is right: the energy axis comes from a discrete measure, and the table's energy column ("needs neither") is refuted.
- The spatial bulk space (`transport/mesh/material_mesh.py:560-590`) is `of_axes(energy, spatial)`. In 1-D the spatial axis is `self.volume_measure.axis("spatial")`: a NODAL axis generated by the mesh's volume measure (nodes at the cell centres, weights the cell volumes). That holds under every scheme.
- Linear discontinuous (LD) does not turn that axis modal. It ADDS a second factor: `SpatialScheme.moment_axis` (`transport/spatial/scheme.py:1693`) mints a MODAL axis `"spatial_moment"` of shape `(2^d,)`, the within-cell tensor-Legendre moments. Its weights are `moment_mass_diagonal`, the unit-volume Legendre mass (`[1, θ]` at d = 1). It has no generator, and its docstring says "the mass IS the scheme's". Diamond difference (DD) is the slopeless closure: width 1, so the moment-tail policy appends no factor (#240's byte-identity invariant). The bulk metric is `G = V·w ⊗ moment_mass_diagonal`.
- That factorisation (cell volume times a reference-cell mass) is admitted only for a constant volume element (`_assert_moment_mass_is_expressible`, the guard reading `has_constant_volume_element`). On a curvilinear cell, `∫ P_k P_l r^{d-1} dr` depends on the cell, and the product form does not hold.

**The sharpening** `[R]`, the orchestrator:
1. **A discretisation chooses a trial space; nodal or modal is the choice of basis within it.** LD's space, discontinuous piecewise-linear functions, factors as (which cell) ⊗ (the polynomial within the cell). The first factor is nodal: the cell partition, its indicator basis and its volume measure. The second is modal in ORPHEUS: Legendre moments. The same LD space could be carried nodally, with two Lagrange points per cell; ORPHEUS chose moments. So DD does not "induce a nodal axis and LD a modal one". Both induce the nodal cell axis. LD adds a within-cell factor, which DD has at width 1 and drops.
2. **Each axis needs a pair, a basis and a measure, because its metric is a Gram; one member of the pair is explicit and the other is canonical.**
   - A nodal axis has an explicit discrete measure (a point set, or a partition with its volumes). Its basis is the canonical one of that set (point evaluation, or the cell indicators), so the Gram is `diag(w)`. The energy axis, the DD cell axis and the SN ordinate axis are all this case. The energy grid states the pair itself (`as_measure` / `as_basis`).
   - A modal axis has an explicit basis. Its measure is the factor's own continuous measure (Lebesgue on the reference cell, or the uniform measure on S²), and its Gram is exact, in closed form. So **the LD moment axis needs no discrete measure**: its mass is a closed-form integral of Legendre polynomials over the reference cell, and nothing is sampled. That is why it carries no generator.
   - So the user's "either a discrete measure or a basis" holds, with the other member canonical. Ruling 4 stands as written.
3. **The frame is not what induces an axis; it is the bridge between two representations of one factor.** In this tree a frame is "a `Basis` bound to a `DiscreteMeasure`" (`FrameBase`'s docstring), with analysis and reconstruction faces between the nodal samples and the modal coefficients. It is needed exactly when one factor is carried in both representations at once, with conversions between them. The angular factor in SN is that case: the sweep runs on ordinates (nodal, the quadrature's axis), and scattering on moments (modal, the harmonic basis's axis). Then the modal axis's metric must be the discrete Gram under the quadrature, not the continuum one, so that analysis and reconstruction stay adjoint to each other (ERR-080's slab frame). That is why the frame re-dresses the head (F4), and why the tree carries one moment space per carrier, the frame's (R-6.2c-1).
   - **LD needs no frame today.** Nothing in the solve converts the spatial moments to point samples: the cross sections are piecewise constant, and their product with a linear field projects exactly. A spatial frame would appear the day a coefficient within a cell is not polynomial (a continuous σ(x), a manufactured source), so that its integrals need a quadrature.
   - **DD needs no frame either.** It carries one representation, the cell averages, whose metric is the volume measure.
4. **So "every axis is what a frame induces"** is [REFUTED 2026-10-18, by the exercise] for the question "what generates an axis". The fact it rested on survives in a narrower form: a modal axis's metric is a Gram, and on a factor bridged by a frame, that Gram is the frame's.

**The concept order this gives** `[R]`, for the user's ruling:
- the discrete measure, and the basis, in parallel, each with its canonical partner;
- the axis, which receives either one: nodal from a measure, modal from a basis;
- the frame, the join of a measure and a basis on one factor, which bridges the factor's two axes and re-dresses the modal one's metric;
- the space, a product of axes.
- The scheme is the spatial factor's generator (CLAUDE.md: "the scheme does the same for the spatial axis"). It chooses the trial space, and so whether a within-cell modal factor exists.

**Questions this leaves for the user:**
- (a) Is the frame a bridge between two axes of one factor, downstream of both, rather than the generator of an axis?
- (b) On a bridged factor, does the modal axis carry the frame's discrete Gram (today's R-6.2c-1)? Then a modal axis has two legitimate metrics: continuum where the factor is purely modal (LD's spatial moments), discrete where it is bridged (SN's angular moments).
- (c) Is the curvilinear LD guard (no product form when the volume element is not constant) the declared edge of machinery not yet built, a cell-dependent moment metric, and so an issue?

## Step 3: the target order, checked mechanically (2026-10-18)

`[M]` the explorer, `scratch/numerics_untangling/step3/report.md` (saved by the orchestrator), scripts and `target_*.json` beside it. Positive control: its builder reproduces `subdag.static_graph()` exactly (91, 110 and 123 edges) and the census's 17 edges.

**The result.** With rulings 5, 10 (three alternatives), 11, 12 and 13 applied, the full graph (top plus lazy plus typing) is still not acyclic. Every residual cycle passes through one of 7 edges; remove those and alternative B2 for the lift is acyclic, in 17 layers. Of today's 17 cycle-closing edges, 5 are gone, 4 remain on no cycle, and 8 remain on cycles.

**The 7 edges, with the orchestrator's reading** (`[R]`, for the user):
1. **`LinearOperator`'s own dunders build the composers** (`__add__` → `OperatorSum`, `__mul__` → `ScaledOperator`, `__matmul__` → `OperatorProduct`). An algebra is closed under its operations, so the abstract type and its composites are one concept, and ruling 5's split cuts through it. The split that ruling 5 needs instead moves the leaves out (the Riesz operators, the axis marginal operators, the leaf operators) and keeps the algebra, its type and its composites, together.
2. **`OperatorSum.inverse` builds a `GreenOperator`,** whose class is `InverseWrapMixin[OperatorSum]` and which wraps the iteration drivers. Inversion is an operation of the algebra, but this one chooses a solution strategy (an iteration) inside the algebra. CLAUDE.md: "the Strategy inverts only M: `M.inverse()` is the sweep". So the Problem/Solution boundary may be crossed here. That is a question to investigate before the split, not a module move.
3. **`SupportsAssembly` annotates `SparseAssembledOperator`** for typing only. The assembly is a homomorphism into the matrix representation, and naming its target type in a protocol is a forward reference: an artefact.
4. and 5. **`Axis.generator`'s union names `FrameBase` and `Quadrature`** for typing. Ruling 4 keeps the generator as provenance until a forgetful-map concept is designed, so these two are owed to that design. No issue tracks it yet.
6. **`manifold`'s catalogue stores `LEGENDRE` as the orbit space S²/O(2)'s reference measure.** The reference measure of an orbit space is the pushforward of the sphere's uniform measure along the quotient map, which is a measure on the manifold, downstream of it. The manifold knowing a downstream measure is a leak. Its one reader is the quadrature registry (`registry.py:1043`), so the knowledge can live there, or be derived as that pushforward.
7. **The memo's boot split**, deliberate: kept with its reason.

**The lift (ruling 10), measured three ways.**
- A: the space mints the lift (a 2-cycle between `space` and the operators).
- B: the operator type builds it, `AxisRetractionOperator.of(space, axis)`. The cycle survives through the Riesz operators.
- B2: B, plus the Riesz operators built on the operator side. The space imports no operator, and the cycle is gone.
- B2 refines ruling 10's "the space owns the lift" into "the space owns the product structure the lift is built from; the operator layer builds the lift". That is #611's rule (raw types upstream, convenient constructors downstream). It re-points 6 production callers of `.retraction(` and `.section(`.

**The small choices the check surfaced:**
- where the fold `DiscreteMeasure.quotient` goes, since it calls `certificate_under` and must leave `measure` with the five invariance verbs: `invariance`, or `quadrature` with the certificate passed in;
- whether the indicator bases' `.space` retires too: it must, or the cycle `{axis, basis, space}` remains (`[M]`, `KEEP_INDICATOR_SPACE=1`);
- `_padded_metric_tensor` moves with `HarmonicAxis`'s new classmethod;
- `ReferenceMeasure` and `OrthogonalSystem` move into `manifold.py`, or upstream of it.

## Candidate ontologies

None yet. The edges' uses are measured first (below). Then each knot is read as one of four cases:
- a leak, where a concept is reached from the wrong side;
- a weld, where two concepts are held in one type, so the type sits on both sides of the order;
- a mutual definition the mathematics really has: a duality, or an adjoint pair defined together;
- an artefact: a deferred import for import order only, or a type hint that could name a protocol.

**E. The `numerics` concepts as sub-cones with their joins** (drafted 2026-10-18 by the orchestrator, for ruling 19; `[R]` from the mathematics, each join with the reason it is a join; the comparison with the code is `[M]` from step 3's report and the edge census).

Read each line as "X: its upstream concepts (the join), and why". A join is a concept that needs two or more upstream concepts at once, the way the frame needs a measure and a basis (ruling 25 of the verification tree).

1. **Primitives:** scalars, units, content identity, small enums, and the metric forms as plain bilinear forms on arrays (`DiagonalMetric`, `DenseMetric`). They depend on nothing.
2. **Symmetry groups:** `SubgroupOfO3`, rigid motions. Group theory, upstream of anything they act on.
3. **Manifold:** a set with charts.
   - **G-space** = join(manifold, group action): the manifold carries its group's action (ruling 13).
   - **Orbit space** `M/H` = join(G-space, subgroup): the quotient, with its `induced_action`.
4. **Continuous measure on a manifold** (`GeneratingMeasure`, `UniformMeasure`, `ProductMeasure`): after the manifold.
   - **The pushforward of a continuous measure** along a map of manifolds = join(continuous measure, manifold map).
   - The **reference measure of an orbit space** is the pushforward of the base's uniform measure along the quotient map, so it is a join(continuous measure, orbit space).
   - Its identification with `LEGENDRE` is a theorem through the Archimedes chart: the uniform measure on S² projects to the uniform measure on the axis.
   - So it sits after both the manifold and the continuous measures. The manifold storing `LEGENDRE` (step 3's edge 6) is this join written into the upstream concept. That answers the user's "what are we joining and when?" for the reference measure, and it is the user's option 1.
5. **Discrete measure** (point masses on a manifold, with the group-free relabelling equality): after the manifold, parallel to the continuous measures.
   - **Exactness claim** = join(discrete measure, continuous reference measure, orthogonal system): a statement relating the two measures.
   - **Gauss construction** (`gauss(n)`) = join(continuous measure, discrete measure): it makes the second from the first.
   - **Invariance** = join(discrete measure, G-space): the relabelling equality after the pushforward along the action (ruling 13).
6. **Basis** (functions on a domain): after the manifold. Parallel to both kinds of measure.
   - **Continuum Gram** (the mass of a basis against a continuous measure) = join(basis, continuous measure).
7. **Axis.** A nodal axis receives a discrete measure; a modal axis receives a basis and carries its factor's pulled-back measure (ruling 17). Its metric's join therefore depends on the factor:
   - join(basis, continuous measure) gives the continuum Gram (LD's within-cell moments);
   - join(basis, discrete measure), the frame's data, gives the discrete Gram (SN's angular moments).
8. **Space** = a product of axes, its metric derived from them: after the axis.
9. **Operator** (a map from a space to a space, with domain and codomain): after the space (the user's reading). Together with it, its **algebra**: sum, product, scaling, tensor product.
   - L(V, W) is itself a vector space under sum and scaling, and composition makes the algebra. These are properties of the set of operators, not of any single operator, so the element type and its operations are ONE concept. That is the user's option 1 for the split: the abstract type and its composites stay together.
   - The dunders on the base class are where Python spells the algebra's operations (Pattern 1), not a dependence of "operator" on "sum".
10. **The metric-bearing constructions** = join(operator, space with a metric):
    - the dual space and `_DualOperator`;
    - the **Riesz operators** ♭ and ♯ (the answer to the user's question: an operator joined with a space's metric);
    - the **adjoint** `.H = ♯ ∘ dual ∘ ♭` (after the Riesz operators and the product).
    - Step 3's Riesz SCC (`{field, functional, linear_operator, space}`) is this join written backwards: `space.riesz_lower` and `riesz_raise` (`space.py:856`, `:876`) mint downstream objects. Alternative B2 (built on the operator side) is the join written forwards.
11. **Field** (an element of a space) after the space. **Functional** (an element of the dual) after the dual space. **Rank-one operator** = join(field, functional, operator).
12. **Frame** = join(basis, discrete measure), as data: the sampling map Φ and the weights. Its **faces** (analysis and reconstruction) are operators between two coordinate spaces, so the frame-with-faces is also after the space and the operator.
    - The **marginalisation lift** (retraction and section) = join(space, axis covector, operator): built on the operator side from the space's product structure and the axis's covector (ruling 10 refined by B2).
    - The **assembled operator** (the homomorphism into the matrix representation) is after the operator algebra.
13. **Inverse family.**
    - A **formal inverse** T⁻¹ is after the operator: a typed object meaning "the inverse of T".
    - Its **realisations** are a separate axis: direct (structural: a triangular solve, the fused `L + C` sweep), or through a **splitting** (the Neumann series of M − N with M invertible). The splitting is the Strategy's `Splitting` (CLAUDE.md, the Solution half), so it sits after the operator algebra AND after the iteration drivers.
    - **Today's `GreenOperator`** is the splitting realisation. Its docstring says "the name is EARNED by the distinguishing G-Neumann invariant". It derives the splitting from the sum's SPELLING ORDER ("the left-spine head is the DESIGNATED preconditioner"). So the algebra's sum encodes which term is implicit, which is the Strategy's decision written into the Problem half: a weld, `[R]`.
    - The answer to the user's question: mathematically every invertible operator has a Green's operator, its inverse, whose kernel is the Green's function. Green-ness is not a property of a sum, nor of `L + C` alone:
      - `(L + C)⁻¹` is the uncollided transport Green's operator;
      - `(L + C − S)⁻¹` is the full transport Green's operator;
      - a sum with no streaming (`C − S`) has an inverse with a local kernel.
    - What the tree's class names is a realisation.
    - **Proposed:** `OperatorSum.inverse` returns the formal inverse. The Strategy realises it, through the splitting it already owns. Edge 2 of step 3 (`green_operator -> composers`, the algebra reaching the drivers) then disappears.
14. **Iteration drivers** (Richardson, power iteration, Krylov), **posing** (the pencil, `EigenPosing`, `SourcePosing`) and **outcomes**: after the operator algebra, posing before the drivers that consume it.
15. **Quadrature:** a composite of a discrete measure, a level structure and its cached angular frames, so after the frame. Concrete spaces (`spaces/`) are after the axis, the space and their bases.

**What E says about step 3's 7 residual edges:**
- edge 1 (the dunders): not a cycle between concepts, the algebra is one concept (item 9);
- edge 2 (Green): the Strategy's splitting written into the algebra (item 13);
- edge 3 (`SupportsAssembly`'s forward reference): the assembled operator after the algebra (item 12); a typing artefact;
- edges 4 and 5 (`Axis.generator`): provenance pointing downstream, owed to a forgetful-map design (ruling 4);
- edge 6 (`LEGENDRE` in the manifold): the reference measure's join written into the manifold (item 4);
- edge 7 (the boot split): deliberate.

**Where E's order differs from the code** (`[M]` step 3): the axis moves after the measure and the basis; the space and the operators swap (today `space` imports operators, while E puts the operator after the space); invariance moves after the measure; the operator layer builds the Riesz operators and the lift; the reference measure leaves the manifold.

## Refuted candidates

None yet.

## Next steps (`[HYPOTHESIS]`, 2026-10-18)

1. **Measure the 17 edges' uses.** For each edge: the importing def, the symbol, what the def does with it, and the claim of any test that executes that def. The census is by AST, with a positive control on `space.py:558`. This step is delegable as one package with a fixed table, and its acceptance check is that the edge set equals the scripts' 17.
2. **Read each knot with the user**, and classify it under the four cases above. This is the ontology search, and it is not delegable.
3. Write the ordered concept list (the `numerics` sub-cone DAG at concept grain), with each concept's home module. Its acyclicity is the done-when.

**The memo-race diagnosis** [REMEDIED 2026-10-18 @829b79d8]: neither red was a memo defect (`scratch/memo_race/report.md`). M3.15 launched its racing processes with a hardcoded `-O` while the memo key carries the optimisation level. B13 counted every file of the shared `.cache/sources`. Both tests now read only their own run.

### The investigation of (a), (b) and the user's attack (dispatched 2026-10-18)

Three perspectives in parallel, each adversarial in its first pass (W5). A separate second pass re-evaluates them; it never bounds the first. Reports go under `scratch/numerics_untangling/frame_bridge/`.
1. **The facts** (explorer). The places where LD, the axis or the space assume a product metric or a diagonal one. The deferred `scheme.mint(chart)` and its CS5 gate. `BundleMeasure`, with its 5 references. The frame's re-dressing. How curvilinear DD carries `r^{d-1}`.
2. **The literature** (literature-researcher):
   - LD and discontinuous finite elements in spherical and cylindrical SN: the moment against the nodal form, the mass with `r^{d-1}`, positivity, and the cell at r = 0;
   - Galerkin quadrature (Morel 1989) and the discrete against the continuum moment metric, which bears on (b);
   - any treatment of a broken polynomial space as a fibred space.
3. **The frames** (cross-domain attacker):
   - (a): is the frame a bridge between two axes of one factor, or a generator?
   - (b): which Gram a bridged modal axis carries, and whether a discrete Gram is a variational crime or the consistent choice (Petrov–Galerkin, Galerkin quadrature);
   - the user's attack: is LD's within-cell factor a fibre whose metric varies over the cells (a vector bundle with a fibre metric), and is the diagonal-only axis metric the gap that makes curvilinear LD unspellable?

### The three reports, and the second pass (2026-10-18)

Reports (each read with session context; the explorer's saved by the orchestrator):
- facts: `frame_bridge/facts/report.md`;
- literature: `frame_bridge/literature.md`;
- frames: `frame_bridge/attack.md`, with probes in `frame_bridge/attack_probes/`;
- one wrong-answer diagnosis it triggered: `scratch/adjoint_underresolved/report.md`.

**The adjoint check** (`[M]`, the numerics-investigator). The attacker found the angular frame's two faces disagreeing on 3 of 37 rule × order rows, all where the order outruns the rule: `folded_product(2,4)` L = 2, LS4 L = 3, LS8 L = 5 (`attack_probes/p6_tree_square.py`, reproduced by the orchestrator). Production's SN adjoint is nonetheless correct.
- Slab forward/adjoint reciprocity holds to round-off on every row.
- A positive control (dropping the group transpose) breaks it to 0.92.
- The production scattering operator is standard S_N, `(1/W) R Λ M`, self-adjoint in the quadrature weights on every rule (Morel 1989, Eqs. 39a, 39b). It never uses the frame's square.
- What was wrong: the `frame.py` docstring claiming the square closes "across every shipped sphere family" [REMEDIED 2026-10-18, this pass's commit], and the absence of any reciprocity gate above L = 1 (the test-architect is adding the rung).

**Second pass: the re-evaluation** (`[R]` the orchestrator, unless marked). It weighs the three first-pass reports against each other; it does not bound them.

(a) **What a frame is.**
- The first-pass reading ("a bridge between two representations of one factor, needed only when both are carried") is refuted twice by the attacker (`[M]`):
  - 7 of the 10 production frame constructions are homogenisation and condensation frames, mapping fine cells to coarse regions: two resolutions, not two representations;
  - on 30 of 33 angular rows the frame is an embedding, not a change of basis.
- The reading that survives both is the attacker's. A frame is the **sampling map** (the basis evaluated on a measure's points), taking a basis's coefficients to a measure's samples, together with its canonical-dual inverse.
  - On the angle, the basis is the harmonics and the measure the quadrature.
  - On homogenisation, the basis is the coarse-region indicators and the measure the fine volume measure.
  - It is a morphism between two coordinate spaces of possibly different dimension. Representations and resolutions are two uses of one construction.

(b) **Which metric a modal axis carries.** The attacker's rule unifies what looked like two legitimate metrics: **every axis carries the pullback of its factor's own measure**.
- A nodal axis pulls back along the identity: its weights.
- A modal axis pulls back along its basis's synthesis: its Gram, `ΦᵀWΦ` when the factor's measure is discrete, the continuum Gram when it is continuous.
- So SN's angular moments carry the discrete Gram, because the quadrature is the angle's measure. LD's spatial moments carry the continuum mass, because the cell's Lebesgue measure is integrated exactly. That is one rule applied to two factors whose measures differ. It agrees with today's tree in both places.
- `[M]` (the attacker) the alternative that ruling R-6.2c-1 rejected was mis-specified: it put the continuum `4π/(2ℓ+1)` on covariant moments (a variance error, failing Parseval on 33 of 33 rows; the variance-correct dual passes on 32 of 33). The ruling's outcome stands under the pullback rule; its recorded evidence does not.
- What is left is one genuine choice, on the rows where the order outruns the rule. There the frame's reconstruction uses the continuum addition-theorem factor while its metric uses the discrete Gram, so the frame mixes two measures.
  - Make one datum single-sourced on the frame (Galerkin quadrature: GQ2-like, Morel 1989 and Sanchez–Ragusa 2011). Every identity then closes, but it changes the scattering on those rows: it is a different method.
  - Or keep standard S_N, which is correct and self-adjoint, and state that the frame's square does not close there. That is what the remedied docstring and the coming rung do.
  - The literature takes both (Shands et al. 2025 builds the moments orthogonal under the quadrature, so the weights survive). This is a method choice for the user, not a defect.

(c) **Curvilinear LD.** The user's attack holds and sharpens.
- The carrier stays a tensor product: a bundle of constant rank over a finite set of cells is trivial. What fails is the metric's Kronecker form, `G = ⊕_c V_c ĝ(τ_c)`, `τ = r_in/h`. It is Kronecker exactly when the volume element is constant, the guard's predicate.
- `[M]` (the facts): no space today admits a metric spanning two axes (the probe was refused on an axis-built space and admitted on an axes-less one). The mass is spelled twice (space side and kernel side). The `scheme.mint(chart)` that would build mass and redistribution as two products of one form was designed and never built, blocked on #423 (no polynomial `Basis`).
- Two repairs:
  1. **A Gauss-nodal spatial axis.** Within each cell, take the Lagrange basis at the 2-point Gauss rule of the weight `r^{max(d−2,0)}`. Both curvilinear LD forms are then exactly diagonal. `[M]` the attacker, ≤ 1.3e-13; the orchestrator re-derived it by exactness of the Gauss rules:
     - sphere: the mass `∫bᵢbⱼ r² dr` is degree 3 under weight r, and the redistribution `2∫bᵢbⱼ r dr` degree 2;
     - cylinder: the mass is degree 3 and the redistribution degree 2 under weight 1.
     The redistribution becomes the multiplier `(d−1)/r` at the nodes. The spatial factor is then ONE nodal axis from a discrete measure (two points per cell, with the volume-element weights), machinery that exists. No moment axis and no cross-axis metric are needed. It is the same Galerkin scheme in another basis, not lumping: the solution is identical in exact arithmetic, and positivity and lumping stay #158's separate choice. The literature agrees: the robust curvilinear schemes are nodal, and lumping is defined only in the nodal basis (Palmer and Adams 1993).
  2. **A modal axis with a fibre metric varying over the cells,** which needs the cross-axis metric machinery the space refuses today.
  - The second pass favours (1): it removes a concept instead of adding one, and it makes the curvilinear guard unnecessary rather than satisfied. The moment reading (cell average and slope) stays available through a frame, the sampling map of the Legendre moments on the Gauss nodes, wherever an output needs it.

### Candidate E attacked, and its second pass (2026-10-18)

The first pass is the cross-domain attacker's, `scratch/numerics_untangling/ontology_e_attack.md`, with probes in `ontology_e/`. The second pass is the orchestrator's, `[R]` unless marked.

**The antichain test** (the attacker's name). A declared join is real only if no parent is an ancestor of another; otherwise it is a step along a chain. `[M]` 9 of E's 14 declared joins pass, all in the measure layer: G-space, reference measure, invariance, exactness, continuum Gram, frame data, modal axis, rank-one operator, quadrature. 5 fail: the orbit space, item 10, item 12's lift, the assembled operator, and the inverse. So the operator layer is a chain. Its cycles come from later concepts written into earlier ones, not from missed joins. The second pass accepts the test as the join criterion, and it becomes one of the mechanical gates below.

**The corrections E takes:**
1. **The Gram datum** (the attacker's name). Item 7 needs the frame's data for SN's modal angular axis, while item 12 places the frame after the space and the operator; E contradicted itself. The concept both need is one join, (basis, measure), producing the Gram: `ΦᵀWΦ` against a discrete measure, `∫uᵢuⱼdμ` against a continuous one. It sits before the axis. The frame's faces (analysis and reconstruction) are operators and sit after the operator. So E's item 12 splits in two: the data before the axis, the faces after the operator.
2. **Exactness and tightness are one concept, with a measured caveat.** The attacker's claim: tightness at order L is equivalent to exactness to degree 2L. `[M]` the orchestrator, `scratch/numerics_untangling/ontology_e/tight_vs_exactness.py` (37 rows):
   - the two agree on 30;
   - 6 cannot be decided, because the folded product rules carry no exactness claim, a gap of their own;
   - one is a counterexample: `gauss_legendre(2)` at L = 2 is tight while exact only to degree 3, because that row is not injective (3 modes on 2 nodes) and the square is read on the analysis range.
   So the equivalence holds for the live, distinguishable modes, which ruling 17's tightness statement already says.
3. **Gauss is not a join.** A Gauss rule is a truncation of the Jacobi matrix (the orthogonal polynomials' three-term recurrence), which also generates the continuous measure (Favard's theorem; Golub and Welsch 1969). So the recurrence sits upstream of both measures, and E's item 5 "Gauss construction" moves there. `GeneratingMeasure.exact_recurrence` is that object today.
4. **Quadrature is a discrete measure whose levels are the fibres of the quotient map,** so it sits before the axis. Its cached angular frames are a cache, not part of the concept. That revises E's item 15.
5. **Item 9's conclusion stands, its reason does not.** Hom-sets being vector spaces needs no composite types (plain matrices have the structure). The composites exist because the algebra is kept as unevaluated TERMS, whose type and constructors define each other. That is why the element type and its composites are one concept.
   - `[M]` `operator.py` defines no `__eq__`, so a sum is an ordered term with no equality. That is where the spelling-order behaviour of the Green operator comes from.
   - Proposed: a sum is an unordered set of signed terms, so `C + L` fuses like `L + C`, and `A − S` equals `(−S) + A`.
6. **Item 13: the formal inverse, with the splitting explicit.**
   - `[M]` (the attacker's spy over 7 gate tiers): 35 `GreenOperator` constructions, 0 from production code outside `numerics`; the diffusion and homogeneous solvers document that they avoid `A.inverse()` on a sum. The fix is small.
   - The docstring's "distinguishing invariant" is the definition of the inverse once the series converges.
   - The formal inverse is a term of item 9's algebra; its realisation takes the splitting as an explicit argument.
   - The `numerics` concept is the abstract `A = M − N`; `Splitting.from_schedule` is `sn/`'s instance.
   - The fused `L + C` inverse is STRUCTURE, because its sweep order is canonical. A Gauss–Seidel schedule or a cycle cut is a Strategy choice.
   - The name `InverseOperator` is taken today (an inverse for a value leaf).
7. **Item 4 holds more widely, item 4's identification less.** `[M]` the pushforward and the descended degree filtration hold for all 3 quotients in the catalogue (O(2), σ and trivial; n = 2 to 5). "Uniform projects to uniform", and the identification with `LEGENDRE`, hold for O(2) only. The join also needs the base's filtration, and a base with an invariant finite measure.
8. **Item 10: the dual space is metric-free.** It belongs with the space, and the dual operator with the operator. Only the Riesz operators, `G⁻¹` and the adjoint join the metric.
   - **A finding, `[M]` by reading:** `DualSpace.of` gives the dual the primal's metric G (`space.py:1376-1385`, "the dual carries the SAME metric as the primal (L²-Riesz)").
   - The metric induced on V* is G⁻¹, and the operator side knows it: `RieszLowerOperator` refuses a dual space ("♭ on V* would apply G where the dual-side map needs G⁻¹ (the G² trap)", `operator.py:1362`), and the dual operator's adjoint routes through the primal (`:1510`).
   - So the tree carries the wrong metric object on the dual, kept from biting by guards. That is guard debt (`coding-standards`, "A guard is elegance debt"). No wrong answer is known; who reads a dual space's inner product is not yet measured.

**Delegation** (verification-tree ruling 18): E is checkable mechanically by 4 gates, each red on today's tree (the attacker's proposal; the second pass accepts it):
- the antichain test on the declared joins;
- imports grouped by concept at symbol grain, typing imports included, forming an acyclic graph;
- each join's constructor takes exactly one type from each parent;
- a fresh interpreter's `sys.modules` after importing a concept contains none of its descendants, which catches lazy imports.

No gate decides whether a join is mathematically real; that stays a reading.

### E′: candidate E revised by rulings 20 to 23 and the second pass (2026-10-18)

The order, by concept. "After" means "imports from"; a join names all its parents. Each item marks what changed from E.
1. Primitives: scalars, units, content identity, enums, and the metric forms as bilinear forms on arrays.
2. Symmetry groups.
3. Manifold; G-space = join(manifold, group action); orbit space = the quotient, a chain step after the G-space (changed: E called it a join).
4. **The orthogonal-polynomial recurrence** (the Jacobi matrix; today `GeneratingMeasure.exact_recurrence`), upstream of both kinds of measure (new: Gauss is not a join).
5. Continuous measure on a manifold, and discrete measure (point masses, with the relabelling equality), in parallel; the Gauss rule is the recurrence's truncation (changed). Exactness = join(discrete measure, continuous reference measure). Invariance = join(discrete measure, G-space). The orbit space's reference measure = join(continuous measure, orbit space), a pushforward; the identification with `LEGENDRE` is O(2)'s theorem only (changed).
6. **Quadrature** = a discrete measure with a fibration, its levels the fibres of the quotient map (moved up from E's 15; ruling 21).
7. Basis, parallel to the measures.
8. **Frame** = join(basis, measure): the binding, its sampling map and its Gram, discrete or continuum (moved up from E's 12; ruling 22). Tightness against the continuum metric is exactness on the live modes (ruling 17, measured caveat).
9. **Axis** = a forgetful image: nodal from a discrete measure (its weights), modal from a frame (its index set and Gram) (changed; ruling 22).
10. Space = a product of axes, with its metric assembled from theirs; the dual space, metric-free, with it (changed: E put the dual after the operator).
11. Operator and its algebra, one concept kept as terms: a sum a formal linear combination, a product an ordered sequence (ruling 23).
12. Chain steps after the operator: the dual operator; the Riesz operators and `G⁻¹` (the space's metric); the adjoint; the frame's faces; the marginalisation lift (built on the operator side, B2); the assembled operator; the formal inverse `InverseOperator`, with realisations (structural, as the fused sweep; or a splitting supplied by the Strategy, ruling 20).
13. Field (after the space), functional (after the dual space), rank-one operator = join(field, functional, operator).
14. Iteration drivers, posing and outcomes, after the operator algebra.
15. Concrete spaces (`spaces/`), after the axis, the space and their frames.

**Findings E′ carries, each owed an issue or a step:**
- the dual space carries the primal's metric, kept from biting by guards (the second pass, item 8);
- the folded product rules carry no exactness claim (the tightness check);
- `Axis.generator` becomes the forgetful map's source (ruling 22).

**E′ checked mechanically** (`[M]` the explorer, `scratch/numerics_untangling/step3/eprime_report.md`; the control still reproduces today's graphs, 91, 110 and 123 edges, and the 17 cycle-closing edges).
- The target graph has 47 nodes, and its top-level part is acyclic. 19 layers, counting each cycle as one node: primitives and the recurrence; content and symmetry; manifold; basis and exactness; measure; the frame's data, the generating measure and invariance; quadrature; axis; space; field; functional; the operator cycle; projection; the frame's faces and posing; concrete spaces; convergence; eigenvalue; iteration; the Green operator.
- Two cycles remain:
  1. `{assembled_operator, operator, operator_steps}`. It is closed by the three composite `assemble` laws, which build `SparseAssembledOperator` (`operator.py:1948`, `:2161`, `:2300`; typed by `SupportsAssembly`, `:90`), and by `LinearOperator.adjoint` and `.dual`, which build `AdjointOperator` and `_DualOperator`.
  2. The memo's boot split, which is deliberate.
- With the verbs and the laws resolved either way, the graph is acyclic in 19 layers.
- Of step 3's 7 residual edges: 3 are gone, 2 now point upstream (the generator edges, ruling 22), and 2 remain on a cycle. Of today's 17: 7 are gone, 8 remain on no cycle, and 2 are on a cycle.
- The explorer's choices:
  - `IdentityOperator`, `PointwiseOperator` and `TensorProductOperator` stay with the algebra, because the algebra builds them; 9 leaves move to `operator_steps`.
  - The frame's faces are the face classes plus 10 methods, `FrameBase.project` among them; the frame's data then imports only `basis`, `manifold` and `measure`.
  - Quadrature is freed by moving `angular_frame` (with its cache and helper) and `Quadrature.axis`.
  - `DualSpace` already sits in `space.py`.
  - The recurrence separates cleanly.
  - A migration is owed for the tests that pin `Quotient.reference` (`tests/gates/numerics/test_manifold.py:1170-1185`).
- **The orchestrator's readings for the two remaining rulings** (`[R]`, put to the user 2026-10-18):
  1. The adjoint, the dual operator and the Riesz operators stay WITH the algebra. With metric-bearing spaces (ruling 22) the operators form a dagger category, so the adjoint is a term constructor of the algebra's structure, like `+` and `@`, and the transpose is the dual functor's action. The Riesz operators are leaves the algebra builds (the explorer's own criterion). Their parents (the space's metric, the dual space) are upstream. And `.H` is the syntax of the Master Standard.
  2. The assembly moves OUT, as a homomorphism from the operator algebra into the matrix algebra. With a sum a formal linear combination (ruling 23) it is a fold over the terms, `[Σ cᵢTᵢ] = Σ cᵢ[Tᵢ]`. Each leaf returns its raw matrix, and the homomorphism, downstream, wraps it.

## ⏸ COMPACTION POINT — 2026-10-18, ruling 27's fix in flight, NOT committed

**Where the work is.** Branch `fix/adjointability-honest`, cut from `main` at `37dc85e0`. The code is UNCOMMITTED in the working tree:
- `orpheus/numerics/{operator,green_operator,coupled_system}.py`;
- `orpheus/sn/loss_representation/__init__.py`;
- `orpheus/sn/operators/{sweep_operator,streaming,radial_characteristic,boundary,windowing}.py`;
- `orpheus/transport/operators/multiplication_operator.py`;
- `tests/gates/numerics/test_coupled_operator.py`, `tests/gates/sn/operators/test_inverse_adjoint_coherence.py`;
- the untracked new gate `tests/gates/numerics/test_adjointability_is_honest.py`.

Pristine copies of every touched `orpheus/` file are under the session scratchpad `pristine/`. That path is lost at compaction; `git diff` against `37dc85e0` is the durable record.

**What the fix does** (ruling 27, the user's rulings "A + B now, C tracked" and "Add it now"):
- A new forward capability, `LinearOperator.is_transpose_invertible`, with the reader `transpose_invertible(op)`: the transposed solve `(Aᵀ)⁻¹` exists and runs. One law: an inverse is adjointable iff its forward is transpose-invertible.
- Implemented by: `StreamingCollisionOperator` (matvec factors AND the new representation capability `has_transpose_sweep`, True only on the 1-D `CumprodScan`); `RadialCharacteristicOperator`; `DiagonalOperator` and `MultiplicationOperator` (`solve_transpose = solve`); `OperatorProduct` (A⁻ᵀB⁻ᵀ); `CoupledOperator`.
- Readers: `SweepOperator.is_adjointable`; `InverseOperator` (new `apply_transpose = inner.solve_transpose`); `CoupledOperator._transpose_substitution_ready`.
- `GreenOperator.adjoint()` = Green of Σ tᵢᴴ, with an honest flag.
- `SNMaskedBoundaryOperator.apply_transpose` = the existing rows branch; `BulkAnalysisOperator` gets Mᵀ ⊕ Id.
- Riesz legs and `RankOneOperator` set to False. **That last one is a regression; see below.**
- The honesty gate: 44 passed, 6 strict-xfail against #614 (the diffusion carrier defect, group C), split by leg (`_CARRIER_BLOCKED_LEG1/2`). Leg 1 threads a cold start for the Green row (`adjoint_act`).
- `test_predicate_flip_is_surgical` is re-posed as `test_the_inverse_family_reads_one_transposed_solve_predicate`. The coupled-operator stubs declare the capability.

**The elegance review: NOT READY** (`scratch/adjoint_honesty/elegance.md`, read by the orchestrator). Owed before commit:
1. **REGRESSION.** `RankOneOperator.is_adjointable = False` breaks `IsotropicFission.apply_transpose`: its kernel `outer(χ, νΣf) & IdentityOperator()` refuses through `TensorProductOperator`'s factor check. `[M]` `tests/gates/sn/operators/test_isotropic_fission.py`: 4 failed (`TestTranspose`). Cause: one flag answers two questions. **Fix:**
   - keep `is_adjointable` = "a transpose is coded";
   - extract posing from `AdjointOperator.__init__` as one `posed(op)`, with the DualSpace-end case tagged for ruling 26;
   - leg 1 splits: `.dual()` iff the flag; `.H` iff flag ∧ `posed(op)`;
   - revert `RankOneOperator` and the Riesz legs to True.
   This is the transposable/adjointable split, owed NOW.
2. **Single source.** Hoist the inverse law onto `InverseWrapMixin`, deleting the copies in `InverseOperator`, `SweepOperator` and `CoupledSubstitutionOperator`, and `SweepOperator.apply_transpose`'s `isinstance` backstop. `MatrixInverseOperator` and `GreenOperator` override. Correct the mixin docstring's "adjoint axis deferred", now false.
3. `transpose_invertible` returns a `TypeGuard` via a `SupportsTransposeSolve` Protocol, removing two `# type: ignore[attr-defined]`.
4. Move `solve_transpose = solve` to `PointwiseOperator` (it is stated twice).
5. **Gate repairs:**
   - fix Green's default cold start in Green itself (it takes the right-hand side's class, not the iterate's), then delete `adjoint_act`;
   - the module docstring's "not marked xfail" is now false;
   - add `TensorProductOperator` and `IsotropicFission` rows;
   - `test_inverse_adjoint_contract_keystone_v2`'s `[InverseOperator]` and `[GreenOperator]` rows still expect False: re-baseline them in the same commit;
   - `SweepOperator.is_adjointable`'s docstring first sentence still names the removed `isinstance`.
6. **Tracked to E′:** the flag is primitive at a leaf (`invertible(A.dual())` needs a `_DualOperator` swap law); `(A+B)* = A* + B*` belongs on `OperatorSum` (ruling 23), after which Green's `adjoint` is `GreenOperator(self.inner.H)`.

**Also pending:**
- the qa review (dispatched; report to `scratch/adjoint_honesty/qa/report.md`);
- the wide serial run (`scratch/adjoint_honesty/wide/summary.txt`, `numerics`, `geometry`, `transport`, `diffusion`, `sn`, `-m "not slow"`). It was started before the regression was known, so expect the fission and keystone reds;
- the `PrescribedInflow` docstring example (`geometry/boundary/prescribed_inflow.py:93`, "ignores psi_out; returns q") is false and is corrected in this commit;
- the finding that the harmonic analysis face's `.H` refuses its own typed codomain element (TypeError in 2-D, ValueError on the slab) gets an issue.

**Issues filed this session:** #610 (memo early cutoff), #611 (the raw-`__init__` constructor-rule campaign), #612 (curvilinear LD, ruled Gauss-nodal), #613 (Galerkin quadrature, investigate), #614 (diffusion composite carriers).

**Rulings 1 to 27 stand** (above). Candidate E′ is the target order, with rulings 24 (assembly out), 25 and 26 (the adjoint stays with the algebra as dagger structure; degenerate metrics replaced by live subspaces) and 27.

⏸ **Next, in order:**
1. Apply the elegance findings 1 to 5 on the branch. Do NOT edit `orpheus/` while the wide run is in flight; check `summary.txt` for DONE first.
2. Read qa's report and fold it in.
3. Re-run the honesty gate, the fission, keystone, coherence, coupled and adjoint files, then the wide run.
4. qa re-check; commit (code plus the gate), Sphinx `-W` with the regenerated `matrix.rst` staged; merge `--ff-only`; push; watch CI.
5. Then E′ cut into packages (verification-tree ruling 18), now with ruling 26's live-subspace move, case (1)'s typed ends and the dual metric among them. Then E′ cut into packages (verification-tree ruling 18). Then E revised in place as E′ [done above], then the step-3 instrument re-run on E′'s moves. The earlier step-3 sketch follows: the ordered concept list: measure and basis in parallel; the axis (nodal from a measure, modal from a basis, each carrying its factor's pulled-back measure); the frame, the sampling map with its dual; the space, a product of axes; quadrature; invariance as a composition; and the operator module cut into the abstract operator, then the assembled and Green operators, then the composers. Each concept gets its home module and the edges its ruling removes; acyclicity is the done-when. Then the packages (verification-tree ruling 18). "Every axis is what a frame induces" is refuted by the exercise. Then step 3: the ordered concept list of `numerics` (groups, the G-space manifold, measure and basis in parallel, the axis that receives either, the frame as their join, the space as a product of axes, quadrature, invariance as a composition, and the operator module cut abstract-then-composers), with each concept's home module and the edges each ruling removes. Its acyclicity is the done-when. Then the work cut into packages (verification-tree ruling 18).

## ⏸ STATE — 2026-10-18 (later), ruling 27's fix revised, under second review, NOT committed

This supersedes the compaction point above for the fix. Its rulings and step-5 list stand.

**The revision** (branch `fix/adjointability-honest`, still uncommitted; the record is `git diff 37dc85e0 -- orpheus tests`, plus the untracked gate):
- **The split.** Two questions, two predicates. `is_adjointable` means a Euclidean transpose is coded, and gates `.dual()`. The new `posed(op)` means the ends give Riesz legs: declared and not a `DualSpace`, or metric-free. `.H` needs both, and `AdjointOperator.__init__` reads `posed`. Rank-one and Riesz are back to True. qa and the elegance review found the same root cause independently (`scratch/adjoint_honesty/qa/report.md` F1).
- **One inverse law on `InverseWrapMixin`**: transpose = forward's transposed solve and dually, plus metric-free inheritance. The `InverseOperator` and `SweepOperator` copies are deleted. `SupportsTransposeSolve` makes `transpose_invertible` a `TypeGuard`.
- **`_DualOperator` swap law** `(Aᵀ)⁻¹ = (A⁻¹)ᵀ` (pulled forward from item 6 above).
- **Green**: `dual()` = Green of the dual sum and `adjoint()` = Green of the adjoint sum, through one `_transposed` helper; `apply_transpose = dual().apply`. Cold start = the first Neumann partial sum, `[M]` bit-identical to zero-seeded.
- **Elegance finding 4 refuted for the question "one home for the self-transpose law"**: `IdentityOperator` has no `solve`, and `MultiplicationOperator` is not a `PointwiseOperator`.

**Measured** `[M]` 2026-10-18:
- the gate: 51 passed, 6 strict xfail (#614);
- the formerly red and adjoint-touching files: 516 passed;
- the mutation battery (`scratch/adjoint_honesty/split/mut_battery.py`): 7 of 8 arms red;
- the radial-flag arm reddens 8 tests outside the gate;
- pyright: 0 errors on the touched files.

**Filed:** #615, the analysis face's `.H` refusing its typed moment.

**In flight:** the wide run (`scratch/adjoint_honesty/wide2/summary.txt`, every gate directory, `-m "not slow"`), the qa re-check (`qa2/`), the elegance re-check (`elegance2.md`), and the archivist (docs, ERR-087, Sphinx `-W`).

⏸ **Next:**
1. Read the three reports and the wide-run summary; fix what they find.
2. Add the `catches("ERR-087")` markers the archivist lists.
3. Commit code, gate and docs together, with `matrix.rst` and `error_index.md` regenerated and staged; merge `--ff-only`; push; watch CI.
4. E′ cut into packages.

## ⏹ RULING 27's FIX — LANDED 2026-10-19, `main` `31f14c93` (+ `ef50f730`, the history hash)

This closes the two state sections above for the fix; their rulings stand.

- `[M]` 2026-10-19, before the commit: the full serial suite with `-O` over every gate directory plus derivations and the root files ran with 0 failed (numerics 4539 passed and 6 strict xfail on #614; sn 3834; derivations 3098; `scratch/adjoint_honesty/wide4/summary.txt`). The xref and catalogue gates passed 46, pyright reported 0 errors on the touched files, and CI `gates` run 39215087307 succeeded.
- The catalogue entry is **ERR-117**; ERR-087 was already taken.
- Green's cold start: the first Neumann partial sum is one counted step, checked before it is driven, and the driver is built per call with the remaining total budget (`[M]` 30 of 30 configurations use the same number of applies; 20 of 20 converged ones are bit-identical).
- A deferred edge recorded in the code: `ScaledOperator` has no native transposed solve. One is owed when a scaled block first sits on a coupled grid's diagonal (the `is_transpose_invertible` base docstring, `orpheus/numerics/operator.py`). No issue is filed, because no consumer exists.
- Filed during the fix: #615 (the analysis face's `.H` refuses its typed moment), #616 (the multi-D transposed sweep), #617 (`DualSpace` carries G, not G⁻¹), #618 (the driver cold start), #619 (`C.apply` accepts a source).

⏸ **Next:** cut E′ into delegable packages, each with an acceptance check (the verification-tree plan's ruling 18). The executor is a Sonnet fan-out or an ultracode session.

## E′ cut into packages (draft, 2026-10-19, `[R]` the orchestrator; open to the user's rulings)

**Purpose.** Execute E′ (the target order above) so that `orpheus/numerics/` imports in one acyclic order, the import order being the concept order. The census of the numerics tests (`verification_tree_refinement.md`, phase 2) then reads each concept against code that already sits at its place. Each package below is written for an agent with no session context (verification-tree ruling 18).

**One correction to E′'s mechanical report, from rulings 25 and 26.** `AdjointOperator`, `_DualOperator` and the Riesz operators stay in the algebra module, because the adjoint is dagger-category structure (ruling 26) and `AdjointOperator.__init__` builds `♯ ∘ dual ∘ ♭`. The report's "verbs kept" variant moved them to the chain-step module, which is what left the cycle `{operator, operator_steps}`. With them kept in the algebra module, that cycle does not arise. The chain-step module then holds only the nine concrete leaves and the marginalisation lift. `[R]`, to be measured by package G0's instrument once the target order is written as data.

**Two kinds of package.**
- **M, mechanical: delegable.** A move, a re-home or a retirement whose target is fully specified here. Behaviour is unchanged. Callers are re-pointed in the same package, with no compatibility re-exports (`coding-standards`, "Retire as you go"). An M package ends bit-identical on every test.
- **S, semantic: not delegable as written.** A package that changes what an object is or what it returns: a ruling realised as code (the sum's normal form, the live subspace, the lift's covector). Each one needs a design pass before its code. By the workflows rule's surgical-carve exception, an operator-algebra carve is written by the main agent with the user steering. This conflicts with verification-tree ruling 18 (delegate), so the executor is a question for the user (below).

**The acceptance instrument every package shares (G0).**
- A new gate, `tests/gates/numerics/test_concept_order.py`, in the idiom of `tests/gates/test_layer_imports.py`. It holds E′'s order as data: each `orpheus/numerics/` module is assigned a rank. It asserts that every import edge (top-level, lazy and `TYPE_CHECKING`, all three kinds) goes from a lower rank to a higher one.
- Its allow-list `PENDING` names every edge that violates the order today. Each entry names the package that removes it. A package's acceptance includes its entries leaving the list, with the gate green. The list only shrinks; an entry added to it is a review blocker.
- The gate's first red, which it must show when it lands (X1): with `PENDING` emptied, it fails on exactly today's violating edges. Its builder is checked against `scratch/numerics_untangling/step3/control.py` (91, 110 and 123 edges at `bcad95f5`).
- The memo pair `{traced_memo, _traced_memo_boot}` is a deliberate boot split (F7), and is a permanent named exemption, not a pending entry.

**Common acceptance for every package**, run by the orchestrator on the return, serially (`reference_test_execution_env`):
1. G0 is green, and the package's `PENDING` entries are gone.
2. Behaviour: `python -O -m pytest -p no:cacheprovider -m "not slow"` over `tests/gates/{numerics,transport,sn,diffusion,homogeneous,geometry,cross_method}` gives 0 failed. For an M package, the passed count equals the pre-package count, except for tests the package lists as moved or retired.
3. `npx pyright` shows 0 new errors on the touched files.
4. `sphinx-build -W` is clean. `dead_references` shows no new row, with a positive control (a deliberately dead `:class:` on a scratch page).
5. The return lists every caller it re-pointed. The orchestrator re-runs the caller census (an AST pass over `orpheus/`, `tests/` and `docs/`, with a known member of each shape) and requires 0 imports of the moved names from their old modules.

### The packages

| id | kind | ruling | what | files (the conflict set) | depends on |
|---|---|---|---|---|---|
| G0 | M | E′ | the concept-order gate with `PENDING` (above) | new test file | none |
| M1 | M | E′ item 4 (E7) | extract the orthogonal-polynomial recurrence into `numerics/recurrence.py`: `_at_digits`, the four closed forms, and the Jacobi and Laguerre recurrences, lifted from closures inside `_jacobi` and `_laguerre` into builder functions of (a, b). The Gauss path stays in `generating_measure`. | `generating_measure.py`, new `recurrence.py` | G0 |
| M2 | M | E′ item 5 (E7) | `Quotient.reference` leaves `manifold`: remove the lazy `LEGENDRE` import (`manifold.py:1659`) and the typing `ReferenceMeasure` import (`:86`); the orbit space's reference measure becomes a function of the measure layer. Its one production reader is `quadrature/registry.py:1043`, and `tests/gates/numerics/test_manifold.py:1170-1185` migrate. **The home is a judgement (owed below).** | `manifold.py`, `generating_measure.py` or `measure.py`, `quadrature/registry.py`, the test | G0, M1 |
| M3 | M | 22 (E3) | split `frame.py` into the frame's data (binding, sampling map, Gram; imports only `basis`, `manifold`, `measure`) and `frame_faces.py` (`CrossGramInverse`, `_dressed`, `_FrameAnalysis`, `_FrameReconstruction`, and the 10 face methods listed under E3, `project` among them). `_collapse_pair` and `_AxisCollapsePair` stay put until S1 deletes them. | `frame.py`, new `frame_faces.py`, 24 + 11 + 13 + 10 production face-method callers (E3's counts, receivers unresolved) | G0 |
| M4 | M | 11 | the axis receives its generator. `DiscreteMeasure.axis` / `.space`, `Quadrature.axis` (5 callers: `sn/problem.py:461`, `:1337`; `sn/mesh/reduced_operator.py:605`, `:670`, `:775`), `basis.space` and `FunctionSpace.for_basis` become classmethods on the axis kinds (#611's rule). Retire `mass_matrix` onto the frame's `discrete_gram` under `retirement-audit`, migrating its closed-form gates with their `verifies`/`catches` markers and the two docs `implements` directives. | `measure.py`, `basis.py`, `axis.py`, `space.py`, `quadrature/directional.py`, the callers | G0 |
| M5 | M | 21 (E4) | free the quadrature: `angular_frame`, its cache and `_harmonic_basis` move to `frame_faces`. Callers: `sn/problem.py:1457`, `transport/frames/harmonic_frame.py:454`. | `quadrature/directional.py`, `frame_faces.py`, 2 callers | M3, M4 |
| M6 | M | 12 | the quadrature cleanup: one rule list, one name per verb, retire the production-dead factories and the unused re-exports, each under `retirement-audit` with its dependency table. **The list of names is owed by the orchestrator** (from the F5 report) before dispatch. | `quadrature/`, `measure.py` | M5 |
| S1 | S | 10 (B2) | the axis owns its marginal covector `m_i = ∫ u_i dμ`; the space owns the product structure; the operator side builds the lift `R = ⟨m, ·⟩ ⊗ id` and `E = Rᴴ(RRᴴ)⁻¹`. A modal axis marginalises, and the energy axis refuses with its own reason. Deletes `_collapse_pair` and `_AxisCollapsePair`, and removes `space.py`'s operator imports (`:86`, `:87`, `:558`, `:871`, `:887`). Re-points 6 production callers of `.retraction(` and `.section(`. | `axis.py`, `space.py`, `frame.py`, `operator.py` | M3, M4 |
| M7 | M | 5 (E1) | the operator module is cut: the nine concrete leaves (`ZeroOperator`, `ZeroMorphism`, `PermutationOperator`, `checked_space_extent`, `TraceRestrictionOperator`, `SumOfTensorProductsOperator`, `DiagonalOperator`, `RankOneOperator`, `outer`) and the lift move to a chain-step module. `[M]` 2026-10-19 the importers by AST: 15 import sites in 9 production files and 72 in 40 test files. | `operator.py`, the new module(s), 49 files | S1 |
| S2 | S | 23 | `OperatorSum` becomes a formal linear combination (an unordered collection of signed, scaled terms, with `__eq__`); `OperatorProduct` a flattened ordered sequence. `(A+B)ᴴ = Aᴴ + Bᴴ` on the sum (the item tracked to E′ by ruling 27's review). | `operator.py`, the sum's consumers | M7 |
| S3 | S | 24 | assembly as a homomorphism downstream of the algebra: the three composite `assemble` laws (`operator.py:1948`, `:2161`, `:2300`) and `SupportsAssembly` move to `assembled_operator`, and each leaf returns its raw matrix. Removes the cycle `{assembled_operator, operator}`. | `operator.py`, `assembled_operator.py`, every leaf's `assemble` | S2 |
| S4 | S | 20 (E2) | a sum's inverse is `InverseOperator(sum)`, realised by a solve the Strategy supplies (its splitting and driver); removes the lazy `GreenOperator` import (`operator.py:1916`). Green's `adjoint` then becomes `GreenOperator(self.inner.H)`. | `operator.py`, `green_operator.py`, `iteration`, the Strategy (`sn/`), every `.inverse()` on a sum | S2 |
| S5 | S | 13 | invariance as a composition: the manifold carries its group's action (a G-space); the measure contributes equality up to relabelling (#403); a support with no action refuses rather than zero-padding (`invariance.py:150`). | `manifold.py`, `measure.py`, `invariance.py`, `rules_product.py` | M2 |
| S6 | S | 26, #617 | degenerate metrics replaced by the live subspace; `DualSpace` carries G⁻¹. Closes case (1) and case (2) of ruling 25. | `space.py`, `metric.py`, `frame.py`, `operator.py`'s Riesz | S1 |

**The waves this order gives** (packages in one wave touch disjoint files):
- wave 1: G0;
- wave 2: M1, M3, M4 (M1 touches the generating measure; M3 the frame; M4 the measure, basis, axis and space; disjoint);
- wave 3: M2, M5, S1;
- wave 4: M6, M7, S5;
- then S2, S3, S4, S6 in series, because each touches `operator.py`.

**Done-when for the whole cut:** G0's `PENDING` is empty; the graph is acyclic with the memo pair as the one named exemption; the full suite (every directory plus derivations and the root files, serial, `-O`) shows 0 failed; and the orchestrator's re-run of the E′ instrument on the real tree reproduces the target layers.

**What stays out of the packages (judgement, owed before dispatch):**
1. M2's home for the orbit space's reference measure.
2. M6's list of names, from the F5 report.
3. The names of the new modules. The session names `frame_faces` and `operator_steps` are placeholders. Proposals: the frame's faces in `frame_faces.py`; the concrete leaves in `elementary_operators.py`; the lift in `axis_marginal.py`.
4. Each S package's design.

**Rulings on the cut** (the user, 2026-10-19, choosing the recommended option each time):
- **28. Split execution.** M packages go to a Sonnet fan-out or an ultracode session. For each S package, the main agent writes the design, the user rules on it, and the main agent writes the code with checkpoints (the workflows rule's surgical-carve exception, read together with verification-tree ruling 18).
- **29. The census interleaves.** Each concept is censused as soon as its packages land, upstream first (verification-tree ruling 22: the cone is built as it is sorted).
- **30. No compatibility re-exports.** Every caller is re-pointed in the package that moves the name, and the old import path stops working.

⏸ **Next:**
1. Settle the three judgement items above: M2's home, M6's name list, the module names.
2. Build G0, and measure the correction above (no `{operator, operator_steps}` cycle with the adjoint kept in the algebra) on it.
3. Dispatch wave 2 (M1, M3, M4) in worktrees.
4. Write S1's design for the user.

**Rulings on the judgement items** (the user, 2026-10-19, the recommended option each time):
- **31. M2's home.** The orbit space's reference measure is one function, `pushforward_reference(orbit_space)`, in `generating_measure.py`, keyed on the quotient's group. The per-entry data (`LEGENDRE` for the O(2) quotients, `None` elsewhere) leaves the manifold's catalogue.
- **32. M6's surviving rule list is `Quadrature`'s classmethods.** Each registry entry references a classmethod and returns a `Quadrature`, and `folded_product` is registered. **Amended 2026-10-19 (wave 4, the orchestrator):** `folded_product` is NOT registered. Its built rule carries no exactness claim, because S²/σ_y has no reference measure, and registering it would mean typing its degree. See the wave-4 section.
- **33. One name per level.** The measure keeps the generic names (`n_points`, `permutation_under`), and the quadrature keeps the angular ones (`n_ordinates`, `ordinate_permutation`). `N` retires, and the two 1e-13 constants become one.
- **Module names** (the orchestrator's call): the frame's faces go in `frame_faces.py`, the nine concrete leaves in `elementary_operators.py`, and the lift in `axis_marginal.py`.

**G0 LANDED** (`tests/gates/numerics/test_concept_order.py`). `[M]` 2026-10-19: its edge builder equals the step-3 control's (108 edges between modules = 123 less the 15 into the package root); a planted backward import is detected; with `PENDING` emptied the order test reddens on exactly 17 edges. **The 17 pending edges, by package:**
- M2: 2;
- M3: 2, plus `frame` importing `space`, shared with S1;
- M4: 6;
- M6: 1;
- S1: 2;
- S3: 1;
- S4: 1;
- S5: 1.

M1, M5 and M7 remove no pending edge. M1 adds the `recurrence` module and M7 two new modules, already ranked; they are moves the order needs for the concepts, not repairs. An edge from `measure` into `invariance` is pending for S5 (the measure's five invariance wrappers).

⏸ **Next:** dispatch wave 2 (M1, M3, M4), each in its own worktree with a brief written to ruling 18's form. Then write S1's design for the user.

**Rulings on M3 and M4** (the user, 2026-10-19, the recommended option each time). Reading the bodies showed that neither package was mechanical as first written.
- **34. The frame's faces are downstream constructors** (ruling 22 read literally, and #611's rule). Each operator-valued face is built from the frame by a classmethod on its own operator class: `FrameAnalysis.of(frame)`, `FrameReconstruction.of(frame)`, `CrossGramInverse.of(frame)`. These replace the private `_FrameAnalysis` and `_FrameReconstruction`, and the `gram_inverse` property. The composing faces (`conjugate`, `reconstruct_after`, `project`) and the space-valued faces (`basis_space`, `test_space`, `measure_space`) become functions in `frame_faces.py` that take the frame first; S7 may re-home the space-valued ones.
  - The faces are `cached_property` values today, and `GalerkinFrame.test_space` returns the cached `basis_space` to keep `is` identity. A downstream constructor keeps both properties by memoising onto the frame instance exactly as `cached_property` does: one helper in `frame_faces.py` that writes into the instance's `__dict__`. A cache is not a concept (ruling 21), so this is legal upstream state written by a downstream module.
  - The Galerkin override of `test_space` becomes a dispatch on the frame's test basis being its trial basis (`frame.test is frame.basis`), the identity the override expresses. That is a one-site essential branch, not a type switch over frame classes.
- **35. M4 splits.** M4a (retire `mass_matrix`) and M4b (the axis mint moves to the axis side: `Axis.of_measure`, `Axis.of_quadrature`) are mechanical. The coefficient space of a basis (`basis.space`, `for_basis`) and the sampled space of a measure (`DiscreteMeasure.space`) become **S7**, designed by the main agent with the user after wave 2. G0's `PENDING` labels are updated to match.

Wave 2 = M1, M3, M4a and M4b. Baseline at dispatch: `main` `e7ad3127` plus this commit. M1's trusted checker is `scratch/numerics_untangling/packages/m1_snapshot.py` (`[M]` on `main`: 234 rules and recurrence rows, digest `e2d578f3…`).

**36. The axis owns its unit, not its covector** (the user, 2026-10-19, refining ruling 10). The axis owns u, the coordinates of the constant function 1: ones on a nodal axis, the expansion of 1 on a modal one. The covector of integration is m = G_k u, from the space's metric block. Then R = mᵀ ⊗ id, Rᴴ = π* = u ⊗ id (exact under the product metric, with no inverse of G), and E = u/⟨m,u⟩ ⊗ id. The nodal case is today's collapse pair. An axis whose basis cannot represent 1 refuses with that reason. The modal and positioned-form refusals retire. Design: `scratch/numerics_untangling/packages/S1_design.md`. The lift's four classes move to `axis_marginal.py` in S1, not in M7.

**S1's gate spec corrected the design** (`[M]` the test-architect, `scratch/numerics_untangling/packages/S1_gate_spec.md`, on `cd9e7e18`):
- Every production `SNProblem.moment_space(L)` head is frame-dressed. It holds covariant moments under the pseudo-inverse of the discrete Gram B (to 1.9e-15, 9 of 9 heads), so the constant's coordinates there are B e₀, not e₀. With u = e₀ the marginal is the integral divided by Σw on 8 of 9 heads, and not proportional on `folded_product(4,4)` at L = 4.
- Reciprocity and R∘E hold for any nonzero u: 38 of 38 rows stay green with the wrong unit. Only the unit and integration gates see it.
- The product-metric condition fails on 0 of 112 (space, axis) pairs over 39 multi-axis spaces: `_require_positioned_over` makes it structural on axis-built spaces.
- The LD spatial-moment axis is MODAL with no generator.
- `dsa.py:595-597`'s "bit-exactly on every rule" is false: R∘E is not bit-exact on 118 of 120 rules (≤ N/3 ULP).
- Across the gates, 179 of 293 rows are red on today's tree, each for its named reason.

**37. The unit is minted with the axis** (the user, 2026-10-19). u is an axis field set by whatever mints the axis, like the weights: ones for a measure; the coefficients of 1 for a basis; the analysis of 1, Φᵀw from the frame's data, when a frame dresses a head; e₀ at the LD scheme's moment-axis mint. This is the stage-2 rule: the forgetful image keeps the induced data, and the generator stays provenance.
**38. The section's divisor is the correctly rounded mass** ⟨m,u⟩ (`math.fsum`) (the user, 2026-10-19). Moved pins are re-baselined under the three-condition rule and listed in the S1 commit (61 to 79 of 120 rules differ from today's bits; the LS4 2-D regression cases by 2 ULP).
- The orchestrator's calls: the "constant not in the span" refusal is dropped (F5: nothing in the tree can trigger it, so it would be an unwitnessed guard); `dsa.py`'s false comment is corrected in the S1 commit; the archivist gives the modal marginal a theory label, so G2 can carry `verifies`.

## ⏹ WAVE 2 — LANDED 2026-10-19 (M1 `2071c98e`, M3 `1ecbc64a`, M4a `cac6d401`, M4b `4e1f2182`, review fixes `9e2dd251`)

Integrated on `refactor/eprime-wave2` and fast-forwarded to `main`.
- **M1.** `orpheus/numerics/recurrence.py` holds `RecurrenceCoefficients`, the private lift `_at_digits`, the four closed forms, `jacobi_recurrence` and `laguerre_recurrence`. `[M]` the trusted checker `m1_snapshot.py` prints the same digest (234 rows, `e2d578f3…`) after the move.
- **M3.** `orpheus/numerics/frame_faces.py`: `FrameAnalysis.of`, `FrameReconstruction.of`, `CrossGramInverse.of`, and the functions `basis_space`, `test_space`, `measure_space`, `conjugate`, `reconstruct_after`, `project`. Two departures from ruling 34's text:
  - The cache is not the `cached_property` slot under the old names. The elegance review found that writing `frame.__dict__["analysis"]` resurrects `frame.analysis` depending on call order. The faces live in one private slot, `frame.__dict__["_frame_faces"]`, behind one helper `_face`.
  - The Galerkin override is `isinstance(frame, GalerkinFrame)`, not `frame.test is frame.basis`. The two differ on `PetrovGalerkinFrame(b, m, b)`: the identity predicate gives it the dressed test space, while the isinstance predicate gives it `test.space`, undressed. The tree had the isinstance behaviour, so M3 kept it. The mathematics favours the identity predicate. **This is a user ruling owed at S7.**
- **M4a.** `Basis.mass_matrix` and its 5 overrides are retired; tests read `GalerkinFrame(...).discrete_gram`. `frame-discrete-gram` lost its `documented` sentinel. It is now verified by `test_frame_discrete_gram_against_lebedev`, which qa's perturb arm reddens.
- **M4b.** `Axis.of_measure` and `Axis.of_quadrature` are classmethods on the axis. The `AngularMeasure` Protocol's `axis` member is deleted, and the reduced-operator factories are annotated `Quadrature`. Recorded, not taken: `Axis.of_measure` ignores `cls` (S7 revisits the axis constructors).
- **Review.** qa passed all five checks: 56 of 56 identities; 9 of 9 markers migrated; the mutation arms redden every migrated test; 0 silent monkeypatches. The elegance review raised two blockers (the cache and three false sentences), both fixed in `9e2dd251`.
- `[M]` wide run on the integrated branch: 17 of 17 test directories exit 0 (numerics 4542 passed, sn 3834, derivations 3098, the root gates 611, …; `scratch/numerics_untangling/packages/wave2_wide/summary.txt`). After the fixes, numerics, transport, data and sn reran with the same counts. Sphinx -W exit 0. dead_references 0 of 66 (positive control: the new nodes are present).
- **One review fix was wrong** (`[M]` CI run 39241651367, red; remedied `84cb71b4`). The elegance NIT "make `mixture.py`'s import form uniform" changed it to `from orpheus.numerics import frame_faces`. An input-layer module must import numerics by submodule (`tests/gates/test_layer_imports.py`), so the old form was load-bearing. The post-fix re-run had left out the root gates. For later waves: the post-fix re-run includes `tests/gates/test_*.py`.
- G0's `PENDING` now holds 13 edges: M2 2, S7 3, M6 1, S1 4, S3 1, S4 1, S5 1.

⏸ **Next:** wave 3, which is M2, M5 and S1. The main agent writes S1 from `S1_design.md` (its REVISED section) and `S1_gate_spec.md`, under rulings 36 to 38. M2 and M5 are delegable, one worktree each.

## Wave 3 (in flight, 2026-10-19): M2, M5, S1

**Rulings at S1's first checkpoint** (the user, 2026-10-19, the recommended option each time):
- **39. The unit is content on MODAL axes only.** `Axis.content_parts` drops `unit` on a NODAL axis, where it is all ones by the kind, so every nodal digest is the one it had before the unit existed (`[M]` a nodal axis's digest, before and after: `26b7af9d…`). On a MODAL axis the unit records the variance of the coordinates (coefficients against moments), which the weights do not fix. The case that decided it: two dense frame-dressed heads of the same order (weights `None`, the form positioned on the space) were EQUAL axes with different units, B e₀ on different rules.
- **40. A MODAL axis is minted with its unit**, and a unit-less MODAL mint is refused at construction. A modal factor that cannot integrate is then unrepresentable, and no guard at the marginal is needed. Cost: 10 test fixtures gained a unit (an AST census, its positive control the new refusal row).

**State.**
- M2 is back on `worktree-agent-acb35fb3873af4fb3`, at `ee80ff47`: `pushforward_reference` lives in `generating_measure`, numerics has 4542 passed, and its 2 PENDING edges are gone. Its review is owed at integration.
- M5 is running.
- S1 is committed as `5f850a69` on `refactor/s1-axis-marginal`. Its wide run is in progress, as are qa (with the spec's §6 mutation arms on the landed module), the elegance review and the archivist.
- **Integration note:** M5 generalises `frame_faces._face`. S1 adds `instance_cache.cached_on`, which is meant as the one home of that memo. At merge, `_face` is re-pointed onto it.

**Ruling at S1's review** (the user, 2026-10-19, the recommended option):
- **41. The axis kind is a closed sum, `BasisKind.NODAL | Modal(unit)`.** The elegance review measured two states the field-based design allowed:
  - `replace(nodal_axis, kind=MODAL)` yielded a modal axis whose unit was all ones, which bypassed ruling 40;
  - a nodal unit of the wrong shape was accepted.
  qa found a third: a test built `Axis(kind="nodal")` with a string, and the old enum read it as modal.
  - The modal arm is the value class `Modal(unit)`. `Axis.unit` is read from the kind, so ruling 39 holds by the type.
  - Construction refuses three things: a kind that is neither arm, a zero unit, and a unit over another shape.
  - It retired the `content_parts` override, the nodal branch of the unit check and the narrowing assert. Landed at `71896877`.

**S1 review outcome.**
- qa: every §6 mutation arm reddens as predicted on the landed module. Two coverage gaps were closed with rows that redden under their own arm:
  - the LD unit row;
  - the closed-form adjoint against the sandwich on the 19 modal rows.
- The elegance review's C1, C3 and C4 were taken:
  - `FunctionSpace._axis_block(k)`;
  - one memo, `instance_cache.cached_on`;
  - `total_weight` renamed `mass`.
- **Owed to S6** (qa F2, which predates S1): `basis_space(frame)` holds two variances. The analysis face writes moments there, while `project` and the reconstruction read coefficients. Φᵀw is the right unit for the moment reading only. No production path reaches the Petrov–Galerkin case.
- **Deferred NIT:** e₀ is hand-minted at 3 sites. It becomes a basis verb ("the coefficients of 1") when a second basis normalisation appears.

**Integration** (`refactor/eprime-wave3` at `7baf444c`): S1 with M2 (`ee80ff47`) and M5 (`e4c3160c`). Its wide run is in progress, and the archivist, qa and the elegance enforcer are confirming.

## Wave 3 LANDED (2026-10-19, `main` `659a4656`)

- `[LANDED 659a4656]` M2, M5 and S1, with the review fixes. The integrated tree was `7baf444c`; the fix commit is `659a4656`, preceded by the plan commit `ee9b9cc0`.
- `[M]` The wide run on `7baf444c` (`-O`, serial, `-m "not slow"`): 17 of 17 test directories exit 0 (numerics 4864 passed, sn 3834, derivations 3098, the root gates 613; `scratch/numerics_untangling/packages/wave3_wide/summary.txt`).
- `[M]` After the fixes, the touched files plus `tests/gates/test_*.py` gave 1364 passed. pyright found 0 errors in the touched files. Sphinx -W exit 0. dead_references 0 of 66.
- The review fixes applied after the wide run:
  - the alias `AxisKind = BasisKind | Modal`;
  - `frame_faces._SLOT`;
  - `_harmonic_basis(support: Manifold, ...)`;
  - the `_overlay_forms` docstring;
  - the modal-kind refusal test split into three parametrized legs.
- The branches and worktrees of waves 1–3 are deleted.

⏸ **Next:** wave 4, which is M6, M7 and S5. Then S2, S3, S4 and S6 in series. S7's design is still to be settled with the user; the PG-degenerate ruling is owed there. qa F2 goes to S6. The numerics test census interleaves with the waves.

## Wave 4 (in flight, 2026-10-19): M6, M7, S5

**A correction to the wave table.** M6 and S5 both edit `measure.py`, so they are not disjoint. M6 and M7 ran in parallel. S5's code waits for M6.

**M6 and M7 returned** (integrated on `refactor/eprime-wave4` at `c27830f5`; its wide run, qa and the elegance review are in progress).
- M7 (`c6f164e2`): the nine leaves are in `elementary_operators.py`. `[M]` numerics 4866 passed, the baseline exactly.
- M6 (`fcd2828e`, `e834dca5`, `1c42a5e0`):
  - the aliases `measure.gauss_legendre` and `gauss_chebyshev` and the property `DiscreteMeasure.generating_measure` are retired, and `PENDING[("generating_measure", "measure")]` is gone;
  - 8 root re-exports are retired;
  - `Quadrature.N` is retired (509 reads rewritten; the 65 kept non-quadrature reads are listed in `M6_callers.txt`), and `_REFLECTION_ATOL` with it;
  - the registry's factories are the classmethods, `build` and `select_quadrature` return `Quadrature`, and the new gate is `test_registry_rule_list.py` (13 rows).
  - `[M]` numerics 4879 passed = 4866 + 13.
- **Orchestrator rulings on M6's open points:**
  - `folded_product` is NOT registered, against ruling 32's letter. `[M]` its built rule carries no exactness claim (25 of 25 sampled), because `quotient` drops the claim: S²/σ_y has no reference measure. Registering it would mean typing a degree, which "values are computed, never typed" forbids. The missing reference measure on S²/σ_y is the prerequisite.
  - The Gauss–Legendre registry key `n` became `n_ordinates`, to match the classmethod (ruling 33).
  - `AngularField.N` (`transport/fields/_bases.py`) is left for a later sweep.
- `measure.equispaced` is kept, by the orchestrator's judgement: it is the only non-Gauss rule on an interval, and the measure-algebra tests use it as their fixture.

**S5 design:** `scratch/numerics_untangling/packages/S5_design.md`. **Rulings** (the user, 2026-10-19, the recommended option each time):
- **42. The invariance verbs stay on the measure, as compositions.**
  - `permutation_under(g) = relabelling_to(pushforward(support.action(g)))`; the certificate, the invariance test and the fold are built on it.
  - `invariance.py` is deleted: `OrbitCertificate` and the lattice walk (`maximal_satisfying`, `maximal_subgroups`, `candidate_groups(points)`) go to `symmetry`, and the windows go to `measure`.
  - This supersedes G0's comment "a composition the measure does not own", which was written before the caller count (about 130 test calls) and the fold's dependence on the certificate were measured.
- **43. The action is typed by an abstract `GSpace(Manifold)`.**
  - It carries `action(motion)` (`Quotient.induced_action` is renamed to this) and `coordinates(points)`, the point identity of a quotient.
  - Sphere, Ball, RealSpace, Circle, Interval and Quotient are G-spaces. On an Interval, x → −x acts iff the interval is symmetric.
  - EnergyGroups, IndexSet, Product and FundamentalDomain are not, and asking them raises `TypeError`.
- **44. The 93 test calls on a bare `Interval[-1,1]`** are re-posed on S²/O2_x through `on_orbit_space`. A test that pinned only the padding's verdict is deleted, with its reason.

**The S5 gate spec returned** (the test-architect, `scratch/numerics_untangling/packages/S5_gate_spec.md`). The baseline is frozen at `c27830f5`: `.claude/worktrees/s5/tests/gates/numerics/data/s5_invariance_baseline.json`, 104 rules, 3242 verdicts, 2306 finite rows. The draft gate module has 1062 rows; `[M]` 38 mutation arms on a model of the design, none blind. Its findings, and what became of each:
- F1, accepted (orchestrator): `relabelling_to` is `permutation_between` plus a weight leg that compares the TWO measures' weights. `permutation_preserving` compares one array, which is right only for a pushforward.
- F2, accepted: walk-versus-bruteforce agreement is the same definition evaluated twice, not evidence of downward closure. That refutes F6's Frame-2 first test for this question. Downward closure keeps its own law gate.
- F5, accepted: an action refuses a translation, and `permutation_under` turns a refused action into None, never a raise, because the three callers of `ordinate_permutation` branch on None.
- F6, accepted: `spherical_product` refuses a polar factor whose support is not `COSINE_INTERVAL` at its top, as the azimuthal factor already does. The pinned row becomes a discrimination row.
- F8: the "93 calls" are 10 tests. 9 are re-posed and 1 is deleted (`test_a_BARE_interval_keeps_the_column_zero_convention`, which pins the padding as the contract).
- F9: `_derived_product_group`'s docstring loses its "independent Hasse walk" credit. The surviving independent pin is the scipy-assignment certificate check.
- **46. (F3, the user, 2026-10-19) A total G-space verb `ambient_points(points)`.** On a Quotient it returns the orbit barycentres in the base's ambient space; on every other G-space it is the identity. It sits beside `ambient_dim`. `is_invariant_under` refuses when `ambient_dim` is not 3. It replaces the private `invariance._embedded_nodes`, which has 15 test reads and is imported at module scope by `tests/_harness/references.py`.
- **47. (F4, the user, 2026-10-19) A motion acts on M iff g(M) = M, with no special case.** `GSpace.action` is ONE law: each family answers whether g preserves the data that defines it, and the action is then `x ↦ g(x)` from M to M.
  - Sphere, Ball, Circle: g fixes the origin.
  - RealSpace(d): E(d).
  - Interval [a, b]: g permutes {a, b}, which leaves the identity and x ↦ a+b−x.
  - Quotient: g normalises H.
  The design table's "x → −x iff a = −b" is that rule's answer on [−1, 1], not a branch. The user asked whether the special case had any advantage; it has none. What is special about [−1, 1] is that it is the chart of S²/O2_x, where μ → −μ is the induced action of σ_x. That is a commuting-square gate to add: `Interval(-1, 1).action(mirror)` equals the quotient's `action(σ_x)` on the chart.

**The wave-4 reviews** (`[M]` qa: ready; elegance: one blocker, which predates the wave). The fixes are on `fix/wave4-review` at `126cf514`, cut from `c27830f5`.
- **ERR-118 (the elegance review's B1).** `positive_weights` was typed on the Lebedev spec, and the selector returned order 13, which has a weight of −0.372. Fixed by reading `positive_weights` and `level_structured` off the built rule.
- **Ruling 45 (the user, 2026-10-19): the zero belongs to the algebra.**
  - `ZeroOperator` and `ZeroMorphism` return to `operator.py`.
  - `0·A` refuses under a SCOPE-BOUNDARY tag until S2, which makes it `ZeroMorphism(A.domain, A.codomain)`.
- **Owed to S2:** `SumOfTensorProductsOperator` is a composite, not a leaf. It imports `operator._agreed_space` and sums by hand. It becomes a refinement of the formal linear combination, or it retires with #260.
- **Kept or filed:**
  - `checked_space_extent` stays in `elementary_operators`: it is a transitional guard that #330 retires.
  - Elegance S4: the `AngularMeasure` Protocol is a second spelling of `Quadrature`'s surface. Filed as #620.
  - `AngularField.N`. Filed as #621.
  - qa F4: `examples/` is broken at import. Already #540, the same defect.
  - `folded_product`'s registration waits on the S²/σ_y reference measure. #370 carries it, with a wave-4 status comment.
- Elegance S2 was done instead of filed: `QuadratureSpec.parameters`, a second copy of the classmethod's signature, retired at `52c9253d`. That closes #371.
- `[M]` The wide run on `c27830f5` passed all 17 suites (`scratch/numerics_untangling/packages/wave4_wide/summary.txt`; numerics 4879, sn 3834, derivations 3098 passed). The fixes were ff-merged onto `refactor/eprime-wave4` at `52c9253d`, and the touched suites were re-run (`wave4_fix/summary.txt`).

**M6 and M7 landed on `main` at `e60330d2`** (2026-10-19). `[M]` The wide run passed all 17 suites, and CI succeeded.

**S5 written** (2026-10-19; branch `refactor/s5-invariance-composition`, code `553ecea8`, docs uncommitted when this was written). The main agent wrote it with two checkpoints, under rulings 42-47.
- **The review rulings** (the user, 2026-10-19, the recommended option each time):
  - **48.** Point identity (`coordinates`) is on `Manifold`, not only `GSpace`, which amends ruling 43. `relabelling_to` is group-free and works on every support.
  - **49.** The measure's two support guards (`_gspace`, `_acted_on_in`) are tagged `ELEGANCE-DEBT[guard] #624`. The machinery is a measure typed by its support.
  - Checkpoint 1: `quotient` moved from `Manifold` to `GSpace`, a Quotient's base is a `GSpace`, and `candidate_groups(points, *, atol)` takes its window as a required argument.
- **The elegance blockers, both fixed, each with a witness that reddens under its defect** (`[M]` an in-process mutation of each, one red each):
  - B1 is ruling 48.
  - B2: the refusal is a typed closed sum, `OutOfDimension | DoesNotMap | DoesNotNormalise`, behind one door, `GSpace.refusal`. It replaces an `except ValueError`, which read an internal bug as "not a symmetry".
- **The elegance concerns, done:** `RigidMotion.preserves` and `permutation_preserving` are retired; `Quotient.acting_dimension` reads its base's; the refusal contract is documented per kind.
- `[M]` **qa: ready.** The compositions match the deleted kernel term by term. 9 of 9 mutation arms red, the positive control 653 reds. The 104-rule record is blind to the weight leg; its witnesses are the synthetic relabelling rows and the LS8 1e-9 control.
- `[M]` **Tests:** `test_invariance_composition.py` has 1079 rows and 43 mutation arms, none blind. The touched suites passed 2436 of 2436 after the fixes.
- **Performance** (`[M]` min of 15 runs): the composition is faster than the old kernel, for example `ordinate_permutation` 0.28 → 0.16 ms and `symmetry_groups` on Lebedev-17 102 → 56 ms.
- **Issues filed:** #620 (the `AngularMeasure` Protocol), #621 (`AngularField.N`), #622 (the chart μ-rule's phase rests on an O(3) tag that S5 made uncheckable), #623 (`Interval(-inf, inf)` against `RealSpace(1)`), #624 (a measure typed by its support).

- `[M]` **The wide run on `553ecea8` passed all 17 suites** (`scratch/numerics_untangling/packages/s5_wide/summary.txt`; numerics 5955, sn 3834, derivations 3098, root 612 passed). Sphinx `-W` is clean; `dead_references` reads 0, and reads 1 with a planted retired name (`induced_action`) as its positive control.

**S5 merged @ `7621fc5f`** (2026-10-20). CI succeeded; the worktree and branch are deleted.

**S2 designed** (2026-10-20; design `scratch/numerics_untangling/packages/S2_design.md`, census `S2_census.md` beside it: an AST pass over 1185 files and a runtime spy over 15586 tests). A sum is a formal linear combination, `{word ↦ coefficient}` with no order; a word is one operator or a flattened product; composition is not distributed over a sum. The census found what stood in the way: the sum's head contract (`is_invertible = a.is_invertible`, used by `GreenOperator`'s left-spine split, pinned by 5 tests and used by no production path), three fused composites chosen by type tests inside dunders (so `L + C` fused and `C + L` did not), the float association of a left-nested sum, and an adjoint that wrapped the whole sum.
- **The design rulings** (the user, 2026-10-20, the recommended option each time):
  - **50. A sum's inverse is named by its split.** S4's first half moves into S2: a sum of two or more terms has no inverse of its own, and `GreenOperator` is built from an explicit split `A = M − N`, as `GreenOperator(M, N)`. The left-spine walk and the lazy import leave. S4 keeps the Strategy driver and Green's adjoint.
  - **51. `ScaledOperator` retires into the combination.** `αA` is the one-term combination; `{A: 1}` collapses to `A`; `(αA)⁻¹ = α⁻¹A⁻¹` stays algebra-closed.
  - **52. `+` never fuses.** The three type-dispatching dunders are deleted. The fused composites (`StreamingCollisionOperator`, `ScheduledInvertibleOperator`, `WindowedSweep`) become atoms with named legs and an `expansion()`, built explicitly by the owner that needs the sweep (the Splitting's implicit part, the loss representation).
  - **53. The adjoint distributes.** `(Σ cᵢwᵢ).H = Σ cᵢwᵢ.H`, and a word's adjoint is reversed. `A.H.H` on a sum is a new object equal to `A`; a leaf keeps `is`.
    - **Its precise form** (the orchestrator, 2026-10-20, from the S2 gate spec's measurements): the dual distributes over every sum and word, and the adjoint is ♯ ∘ dual ∘ ♭ taken once at the composite's ends. The literal `Σ cᵢ wᵢ.H` misses reciprocity by 1.26e-3 when a term has no spaces (a pointwise term on a dense metric), and `B.H @ A.H` misses by 6.4e-2 when the middle space's metric is singular (`[M]` `scratch/numerics_untangling/packages/S2_gates/probes/p1_pointwise_dense.py`, `p2_word_adjoint_singular_middle.py`). Those two forms are theorems where every metric involved is non-degenerate, which is everywhere after S6.
  - Taken as recommended, with no objection: terms collect by object identity (D2), the zero is the empty combination (D6), a product keeps identity factors (D7), and `SumOfTensorProductsOperator` retires with its 8 tests re-pointed (D8).

**S2 gate spec written** (the test-architect, 2026-10-20; `scratch/numerics_untangling/packages/S2_gate_spec.md`, with the instruments under `S2_gates/` beside it). `[M]` 40 rows in `tests/gates/numerics/test_linear_combination.py`: 29 red and 11 green on today's tree. On the spec's prototype of the ruled design, 36 are green; the 4 the prototype does not model are red. A 33-arm mutation battery reddens every row the prototype models, and its positive control reddens 15. The fusion blast set has 3 production owners: `sn/coupled_system.py:553`, `sn/splitting.py:448` and `sn/solver.py:1007`. It also has 78 test sites, 75 of them `L + C`. About 95 re-baselined sites in 30 files is `[R]`; the landed code's red loop decides the number. #625 is filed: `PointwiseOperator.H` on a dense metric; production does not reach it, and its home is S6.
- **54. A pin moved only by flattening a named sum is re-baselined** (the user, 2026-10-20). The homogeneous loss `C − K_iso`, with `K_iso = S + 2Σ₂ᵀ`, now folds as `(C − S) − 2Σ₂ᵀ`. That moves the coda record `homo_2eg_n2n` at A[0,0] by 1 ulp, and the library mixtures do not move. That one row is regenerated under the three conditions: named intermediates, the raw-mixture reference as the independent anchor, and drift of one more level of addition. `K_iso` keeps its name, and the record's docstring is annotated with the date. Rejected: keeping `K_iso` an atom so the bytes stay, because that re-welds what 52 removed; and spelling the loss to reproduce the old order of addition, because the code would then serve a floating-point order.
- The orchestrator's calls from the spec, made 2026-10-20:
  - `0·I` (all cancelled terms pointwise) is `ZeroOperator()`, and a combination refuses only when a non-pointwise term has no spaces.
  - The fold subtracts at −1 and accumulates in place only into a buffer it owns. `[M]` peak 3.0086× against the 4.0× ceiling of `test_p3_peak_allocation_per_dof`.
  - `ScheduledInvertibleOperator.expansion()` is the flattened `L + C − B_lower`.

**S2 merged @ `4e5ee36b`** (2026-10-20). The main agent wrote the production carve with checkpoints. The test-architect re-baselined 46 test files and one fixture (110 log lines: RE-POSE 77, INVERTED 12, DEMOTED 11, DEAD 4; `scratch/numerics_untangling/packages/S2_rebaseline_log.md`). The archivist wrote the theory pages and minted 11 labels, which 35 gate rows now `verifies`.
- `[M]` The mutation battery on the landed code: 39 of 39 rows reddened by at least one arm, none blind (`S2_rebase/mutation_battery_landed.md`). Rows c40–c43 were added at review, each with its own reddening mutation.
- `[M]` The bit probe, run by qa and by the test-architect independently: 55 SAME, 5 refused alike, 10 MOVED, all `homogeneous_A2g_n2n`, at most 8.8e-16 relative (ruling 54). The S_N loss `A_AA` and every pencil `at(σ)`, forward and adjoint, are SAME.
- **A latent defect, fixed in S2 by the user's ruling:** on an LD carrier the collision multiplier broadcast σ over the ordinate axis only, so `C.apply` refused alone; the `L + C` fusion had hidden it. `_broadcast_axes` reads the axes σ does not span by label (ruling 56 for the label-less case).
- **The S_N realizer's α-branches collapsed** into `α·(K_ω ⊗ I)` (the user's ruling): `[M]` 1619 of 1619 calls gave the dedicated zero's exact spaces.
- **Review fixes:**
  - the `FusedAtom` base, whose non-realised faces come from one cached expansion;
  - `_Ends` retired;
  - `composition` keeps a one-term combination's ends;
  - one refusal builder;
  - roles on the leaves (ruling 55);
  - the fold accumulates in place only where exact (qa F2, c43);
  - the transpose seed scaled like `apply` (qa F3).
- **Declined, with the reason:** building `M` at the labelling site (elegance finding 3). The labelled terms are the single source; a stored `M` would go stale under a re-labelling.
- **Owed to S4:** ruling 57, the scaled-inverse seed paths; and `GreenOperator._gains` keeps N as a tuple beside the combination.
- **Owed to S3:** the elegance review's note that a fused atom's `assemble` now follows its expansion, so assemblable leaves reach the atoms with no edit.
- **The gates before the merge** (2026-10-20):
  - `[M]` the wide run `-m "not slow"`: 18 739 passed, 0 failed, 261 skipped, 59 xfailed (78 min); the three files edited after it started, plus the 11 root gate files, re-ran green (700 passed).
  - `[M]` pyright on `orpheus/` (the six core folders' configuration): 0 errors; the planted probe was reported, so the zero is read.
  - `[M]` Sphinx `-W` exit 0; `dead_references` 0 of 68, with a planted `ScaledOperator` reference read as 1.
  - The LD defect is catalogued as ERR-119, caught by 3 decorated LD tests; the DD rows are its controls. The re-dropped defect reddens exactly the LD rows and the refusal row.

**S3 written, reviewed** (2026-10-20). Assembly is a homomorphism out of the algebra (rulings 24, 60–62). `operator.py` keeps only the flag `is_assemblable`; `assembled_operator.py` owns the functor `assemble(op)` (the sum, word and fused-atom laws, then the leaf arm, which canonicalises each emission to CSR), `dense(op)`, `SupportsAssembly`, `assemblable` and `MissingAssembly`. `as_matrix` probes only; `MatrixInverseOperator` densifies through `dense`. Eight classes implement `assembled_matrix(emit)`; the three domain composites recurse through `emit`. G0's S3 edge is removed (planted-control checked) and `matrix_inverse_operator` has its own tier.
- `[M]` The census (`scratch/numerics_untangling/packages/S3_census.md`): of 1 245 base `as_matrix` calls over five suites, 300 took the assembly route, and the only production site among them was `MatrixInverseOperator.__init__` (159).
- `[M]` Bit identity: 176 of 176 LU input matrices over 139 tests identical before and after (the test-architect); 0 of 112 CSR and dense records differ over 56 operators against a HEAD tree (qa, with a one-byte corruption as control).
- `[M]` Gates: `tests/gates/numerics/test_assembly_functor.py`, 25 rows plus the review's G1h, G10 and the size row; the battery's 22 arms each redden their target rows (`S3_gate_spec.md` §7–8). 14 test files re-baselined, with marker counts unchanged.
- **Review fixes:** the refusal names the deepest part that cannot emit (qa F1); the fused atom's expansion is built once and read through `FusedAtom.expanded`, and the atom's flag is final (the user's ruling on elegance F1); the size refusal is worded for both algorithms; the fold's dead `in_place` flag is removed (elegance F6); five stale docstrings.
- **The user's ruling on elegance F2:** the composites' derived flags stay beside the functor, like the other capability axes, kept equal by gate G10 (25 registry operators: `is_assemblable` iff `assemble` succeeds).
- **Declined, with the reason:** naming the `emit` callback type once (`flat_operator` shares `assembled_operator`'s tier and cannot import it); `assemble_ordinate_blocks` building its own carrier (a consumer wrapping matrices it holds); public names for the private helpers `_left_fold`, `_resolve_basis_shape`, `_check_dimension` (a cross-module import of a private helper already had a precedent in `matrix_inverse_operator`).
- **Follow-ups filed:** #626, threading the full law path into `MissingAssembly` (elegance F5); #627, writing the domain composites as atoms whose expansions use restriction and extension maps, so the functor owns every law (elegance (d)).
- `[M]` The wide run `-m "not slow"` on the settled tree (2026-10-20, `-O`, serial): 18 792 passed, 0 failed, 261 skipped, 59 xfailed (60 min). A first run had 1 failure, a traced-memo "skeleton changed" error caused by editing `operator.py` while that run was in progress; the test passed in isolation, and the second run is the one recorded here. pyright: 0 errors besides a planted probe (the control). Sphinx `-W` exits 0, regenerating `matrix.rst`; `error_index.md` is unchanged. `dead_references` reads 0 of 68, and 1 with a planted `orpheus.numerics.operator.MissingAssembly` docstring reference as its positive control.

⏸ **Next:** merge S3; then S4 under rulings 65–67, then S6, then S7 under rulings 58, 63 and 64.
