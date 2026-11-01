# Edges as tested contracts — each import edge of the production graph states what crosses it, and a test proves the provider serves it

Living plan (`plan-authoring` §0), opened 2026-10-31. The ontology is being searched; implementation waits for the
user's ruling that this plan is polished. Sequencing: after `layer_dag.md` Phase 2 (the user's ruling of 2026-10-31,
memory `feedback_testing_graph_before_features`): the census below is measurement and may run earlier; the build may not.

## 1. Goal, in the domain's terms

The layer DAG (`layer_dag.md`) makes the production modules nodes and their imports edges. A node is established when
its taxonomy is settled; the tool for a taxonomy is the class hierarchy (nominal: what a thing IS). An edge is a
dependence of one node on another; the goal is that every edge that carries a SERVICE states its contract as one named
object, and that a test proves each provider honours the contract, so that:

1. a change that alters the MEANING of what crosses an edge (its law) reddens a test at the edge, at the lowest node
   that can see it, before any downstream pin;
2. the surface a consumer depends on is explicit and countable (today it is whatever attributes it happens to read);
3. the cached test graph (`verification_tree_refinement.md`, the cone store) gains named contract nodes.

Proposed means `[HYPOTHESIS 2026-10-31]`: `typing.Protocol` for the static half of a contract, a parametrised law suite
for its behavioural half.

## 2. What a contract is — two halves

- **Static half (spelling).** A `Protocol`: method names and signatures. Checked by pyright where an object is passed
  to a parameter annotated with the Protocol, or at a class that subclasses the Protocol explicitly.
  `runtime_checkable` adds an attribute-presence check only. Neither checks meaning.
- **Behavioural half (law).** The meaning of what crosses: ⟨Ax, y⟩ = ⟨x, Aᵀy⟩; R∘E = id; flux in the positive cone; a
  quadrature exact to its declared degree. A law suite parametrised over every conformer is the only proof. `[R]` In
  ORPHEUS the law half carries most of the value, because pyright is not a hard gate today (`[M]` 226 errors on the
  branch against 232 on `main`, 2026-10-31, `npx pyright` with the worktree config).

**Limit, stated once.** A contract localises WHICH promise broke. It does not cut numerical dependence: a provider can
keep every law and change its values (another quadrature, another scheme order), and every downstream pin reading
those values stays in its cone. `[R]` So contract nodes sharpen the cone store; they do not replace its numeric edges.

## 3. The ownership question — where the Protocol lives

### 3.1 Disambiguation: two directions on one edge

"Upstream" is ambiguous on an edge, and the two readings below were arguing about different edges until separated.

- **Import direction.** The higher-knowledge module imports the lower (the layer criterion). Call the lower end the
  LOWER node.
- **Service direction.** Who calls whom with what. Two kinds of edge, both present in the tree:
  - **(i) Service edge.** The higher node imports the lower and calls it; the LOWER node is the PROVIDER. Example:
    `sn` reads a `Quadrature` from `numerics`.
  - **(ii) Framework edge (inversion of control).** The lower node defines an algorithm that consumes an object it
    cannot import; the higher node implements the contract and passes the object down; the LOWER node is the
    CONSUMER. `[M]` 2026-10-31, `git grep`: `numerics/eigenvalue.py` `EigenvalueSolver` is implemented by
    `diffusion/solver.py:196` and `cp/solver.py:1132`; `reference/verification.py:122` `ProductionAnswer` is consumed
    by the reference layer and satisfied by L3 answers; `numerics/operator.py` `LinearOperator`.

### 3.2 Reading D — the consumer owns the Protocol (Go / hexagonal style)

Arguments for:
- D1 Interface segregation: the consumer lists only what it calls, in its own vocabulary; an owner-side Protocol tends
  toward the provider's whole surface (a "header interface"), which restates the class (a Pattern-2 duplicate).
- D2 It surfaces the consumer's implicit assumptions (an ordering, a sign, a normalisation the provider never promised).
- D3 Dependency inversion: the abstract side should not depend on the concrete side.

Stress test:
- On a framework edge (ii) the consumer IS the lower node, so D and Reading U coincide: no conflict.
- On a service edge (i) the consumer is the higher node. Then the provider cannot import the Protocol (layer gate), so
  it can neither declare conformance (`class Quadrature(Proto)`, pyright-checked at definition) nor host the law test.
  The law test lands in the consumer's test tree, so a change of meaning at the provider reddens later and wider: in
  every consumer's tree, once per consumer. `[R]`
- D3 inverts in this codebase. In Go the consumer (business logic) is the abstract side and the provider (a database) the
  detail. Here layers are ordered by knowledge: the lower node is the MORE abstract one (`numerics` is mathematics).
  The Stable Abstractions Principle (Martin) says dependencies point toward stable, abstract packages, which is the
  lower node. `[R]`
- D2 is real, but Pattern 7 (normalise at the definition site) says a convention a consumer relies on belongs to the
  producer; the repair is to PROMOTE the assumption into the provider's contract, not to keep it consumer-side.
- D1 is real and is the main risk of Reading U; it is answered by measurement (§4, M2), not by ownership.

### 3.3 Reading U — the provider (upstream) owns the Protocol (the user's reading)

Arguments for (the user's, 2026-10-31):
- U1 The consumer already imports the provider, so the Protocol adds no edge.
- U2 The provider can test that it serves the Protocol, so a change of the meaning (the law) of what it serves reddens
  at the provider.

Stress test:
- U1 holds, and holds for Reading D too (a structural Protocol adds no import either way). The discriminating point is
  where the CONFORMANCE TEST can live: under U it lives in the provider's own test tree, imports nothing above it, and
  sits at the root of the provider's cone, so it reddens first. That is the property the cached test graph wants.
- U2 holds and is the stronger argument: a law is a property of the produced quantity (Pattern 7), so its owner is the
  producer. U also allows explicit declaration (`class P(Proto)`), a definition-site check D cannot have.
- Risk U-a, header interfaces: the Protocol grows to the provider's full surface and constrains parts no consumer
  uses. Mitigation: the Protocol's surface is the union of what consumers use, measured (M2).
- Risk U-b, several providers in different lower nodes serve one consumer: "the provider's node" is then not one node.
  Resolution: the Protocol lives at the meet, the lowest layer every provider imports.
- Risk U-c, a framework edge: the provider is the HIGHER node and cannot host a Protocol the lower consumer must import.
  Resolution: the Protocol lives at the lower node, which is the consumer there.

### 3.4 Synthesis `[HYPOTHESIS 2026-10-31]`, awaiting the user

Neither "upstream" nor "downstream" is the rule; the layer criterion already is: **a Protocol is a module-level concept,
so its home is the lowest-knowledge layer whose vocabulary suffices to define it.** On both kinds of edge that is the
LOWER node (or the meet of several), because only the lower node can be imported by both ends.

- Service edge (i): the lower node is the provider, so this is Reading U, and the provider's tests prove conformance
  and law.
- Framework edge (ii): the lower node is the consumer, so this is Reading D, and the conformance tests for each higher
  implementer live in the implementer's tree (it imports the Protocol), with the law suite itself defined once at the
  lower node and parametrised by the implementers (as the operator-law gates are today; to be verified, M4).
- A Protocol in a HIGHER node is legitimate only when its vocabulary is the higher node's own, i.e. it names a role
  that exists only there. `[M]` 2026-10-31 an example candidate: `mc/solver.py:76` `MCGeometry` (to be read before
  classifying, M1).

`[M]` 2026-10-31 the tree already leans this way: of the 42 `Protocol` classes in `orpheus/` (AST: a `ClassDef` whose
base name ends in `Protocol`; script `scratchpad/protocols.py`, not tracked), 23 are in `numerics/`, 8 `transport/`,
4 `sn/`, 3 `reference/`, 2 `mesh/`, 1 `data/`, 1 `mc/`. There are 20 direct `ABC` subclasses (6 numerics, 6 sn,
5 transport, 2 mesh, 1 geometry). Whether each Protocol is service-edge or framework-edge is not yet measured (M1).

## 4. Measurements owed (each with its predicate; none run yet unless marked)

- **M1** Classify the 42 Protocols: (i) service, (ii) framework, or (iii) higher-node role; for each, its implementers
  (explicit subclass vs structural only) and whether a law test exists. Instrument: AST plus `protocol_conformers`
  (Nexus; needs `scripts/setup.sh --no-docs` first).
- **M2** For every cross-node import edge of the production graph (the `layer_dag.md` node set), the surface the
  consumer uses from the provider: names imported, attributes read on the provider's types, parsed by AST with a type
  census (`retirement-audit`, members spelled without the symbol). Classify each edge: (a) value type, the type is the
  contract; (b) service with ≥ 2 providers or a narrow slice of a wide class, a Protocol candidate; (c) single provider
  used in full, stays concrete; (w) weld, the consumer reads a provider's private names.
- **M3** For each (b) edge, the laws the consumer relies on, read from the consumer's code and docs, and whether a test
  pins each today.
- **M4** How the existing law gates are organised (the operator adjoint gates, the frame and quadrature gates): defined
  once and parametrised over conformers, or re-written per conformer.

## 5. Candidate ontologies

- **C1 Protocol per edge.** Every edge gets a Protocol. Unspellable: an undeclared dependence. Welded: value-type edges
  gain a duplicate of the type. `[R]` over-builds; expected to fail Pattern 6's concept-count test on (a) and (c) edges.
- **C2 Protocol per service role, at the lower node, with a law suite** (§3.4). Unspellable: a provider whose law
  changed with no red at its own node. Welded: nothing known yet.
- **C3 ABCs instead of Protocols.** Nominal; conformance enforced at instantiation for abstract methods. Welds taxonomy
  with contract: a class serving two edges needs two bases, and test doubles must subclass. `[R]` kept for taxonomy only.

## 6. Refuted candidates

- **"The consumer owns the Protocol" as a general rule** [REFUTED 2026-10-31 for this codebase's service edges, by the
  stress test in §3.2]: on a service edge the provider could not import it, so neither declare conformance nor host
  the law test; and the dependency-inversion argument behind it inverts where the lower layer is the abstract one. The
  FACT it establishes: the consumer's usage is the right MEASURE of a Protocol's surface (M2), not the right home.

## 7. Open questions for the user

- Q1 Is the synthesis (home = lowest layer whose vocabulary suffices) the rule you want, or do you want a stricter
  "always the provider" even on framework edges (impossible under the layer gate, as far as `[R]` reasoning goes)?
- Q2 Does "prove it serves the Protocol" mean the law suite (behavioural), or also a static gate (pyright on the
  conformance assignments), given pyright is not yet a hard gate?
- Q3 Scope of the first census: `numerics`, `transport`, `sn` only (the six core folders' busiest edges), or the whole
  production graph?

## 8. Rulings ledger

- 2026-10-31 The user: open a living plan; discuss and measure before planning the build.

## 9. Start condition

Implementation starts when the user rules this plan polished, and not before `layer_dag.md` Phase 2 lands.
