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

⏸ **Next:** rulings 16 to 18 recorded, and ruling 17's method choice closed (standard S_N; Galerkin quadrature #613). Step 3, the ordered concept list: measure and basis in parallel; the axis (nodal from a measure, modal from a basis, each carrying its factor's pulled-back measure); the frame, the sampling map with its dual; the space, a product of axes; quadrature; invariance as a composition; and the operator module cut into the abstract operator, then the assembled and Green operators, then the composers. Each concept gets its home module and the edges its ruling removes; acyclicity is the done-when. Then the packages (verification-tree ruling 18). "Every axis is what a frame induces" is refuted by the exercise. Then step 3: the ordered concept list of `numerics` (groups, the G-space manifold, measure and basis in parallel, the axis that receives either, the frame as their join, the space as a product of axes, quadrature, invariance as a composition, and the operator module cut abstract-then-composers), with each concept's home module and the edges each ruling removes. Its acyclicity is the done-when. Then the work cut into packages (verification-tree ruling 18).
