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

⏸ **Next:** step 1, the edge-use census.
