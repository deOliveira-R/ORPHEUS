# The `numerics` layer's concepts in one acyclic order — a living plan

Opened 2026-10-18 from the verification-tree campaign's sub-cone measurement (`verification_tree_refinement.md`, finding S13; ruling 27 there). Status: **the ontology is being searched** (`plan-authoring` §0). Nothing here is implemented. Implementation starts when the user rules the plan polished. Execution is to be delegable (the verification-tree plan's ruling 18), so the plan's final form cuts each step into packages with an acceptance check.

## Goal, in the domain's terms

The `numerics` layer holds the mathematics every method rests on: axes, discrete measures, bases, the frame that pairs a measure with a basis, quadratures, function spaces, symmetry groups and invariance, operators. Each concept is defined from the concepts before it and never reaches a concept after it, so its claims can be verified using only what comes before. The concepts form one acyclic order with joins: for example, the Discrete Measure and Basis concepts both feed into the Frame (the verification-tree plan's ruling 25). That order is the base of the whole verification cone, so every sub-cone above rests on it being right.

## The user's instructions (verbatim, 2026-10-18)

> "2. settle it first." (the numerics knot, before the census)

> "3. It's fine that the SN is essentially a big blob depending on SNProblem *for now*. This is a fairly strong indication that we have not been diligent enough with our rule of "__init__ constructor being the raw constructor for testing and class methods the specialized constructors that receive convenient arguments that needs unravelling", so we need to revit it in its own campaign. Untangling the numerics layer is a pretty good start."

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

## Candidate ontologies

None yet. The edges' uses are measured first (below). Then each knot is read as one of four cases:
- a leak, where a concept is reached from the wrong side;
- a weld, where two concepts are held in one type, so the type sits on both sides of the order;
- a mutual definition the mathematics really has: a duality, or an adjoint pair defined together;
- an artefact: a deferred import for import order only, or a type hint that could name a protocol.

## Refuted candidates

None yet.

## Next steps (`[HYPOTHESIS]`, 2026-10-18)

1. **Measure the 17 edges' uses.** For each edge: the importing def, the symbol, what the def does with it, and the claim of any test that executes that def. The census is by AST, with a positive control on `space.py:558`. This step is delegable as one package with a fixed table, and its acceptance check is that the edge set equals the scripts' 17.
2. **Read each knot with the user**, and classify it under the four cases above. This is the ontology search, and it is not delegable.
3. Write the ordered concept list (the `numerics` sub-cone DAG at concept grain), with each concept's home module. Its acyclicity is the done-when.

⏸ **Next:** the user's readings of the F3 to F6 findings (2026-10-18); the memo-race diagnosis in parallel.
