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

⏸ **Next:** F1 and F2 ruled; F3 to F6 investigated in parallel (dispatched 2026-10-18; reports under `scratch/numerics_untangling/`): F3, the space, the frame and marginalisation; F4, the consumers of `mass_matrix` and `basis.space`; F5, the quadrature API of the measure family; F6, which measures ask invariance questions, and whether that marks a specialisation. Then the readings with the user.
