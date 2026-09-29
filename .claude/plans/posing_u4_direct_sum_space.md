# A system's unknowns live in one labelled direct-sum space, and every block structure is one partition (posing unit 4)

Status: SCHEDULED, not opened (2026-09-29). Issue #527; charter issue #522. Charter: `.claude/plans/posing_sequence.md`, "The work, ordered by dependence"; ruled content: "The ontology as it stands", "Layer 1: the direct-sum space and the trace" and the `TransportMethod` paragraph. Evidence: the prototype and its two attacks, `scratch/posing_sequence/open_items/coproduct/` (`memo.md`, `proto/`, `attack_structure/`, `attack_elegance/`, `proto2/`, `proto3/`). Counts `[M]` 2026-09-29 by `git grep -n '\b<symbol>\b' -- orpheus tests derivations tools`. Depends on units 0 and 2 only through names.

## Goal, in the domain's terms

Bulk and trace, neutrons and precursors, and the nodes of a discrete axis are all summands of a direct sum; a system's operator is a block operator over it. Today two carriers (`FullFieldSpace`, `CoupledSpace`) spell extension by zero two ways and two enums (`BlockRole`, `SystemRole`) declare block structure that nothing reads. After this unit there is one space, one block operator, and the block pattern is derived from the operators' own structure.

## Scope (ruled 2026-09-28/29)

1. `DirectSumSpace` of labelled summands; `FullFieldSpace` and `CoupledSpace` become it (the prototype reproduces both bit for bit: 77 gates, 48 of 48 mis-ordering controls red). The points live on the measure: `DiscreteMeasure.__add__` widened to differing supports, with a support type (not built by the prototype). Labels are part of the identity (a label gate for isomorphic summands).
2. `CoupledOperator` keyed by (row, column) label is the one block operator; `from_blocks` returns it; block substitution lives here once.
3. The injections carry `assemble()`; an alien field is refused (`x.space == self.domain`); carrier access through the #391 `parts`/`with_parts` protocol; `Composite` stays the field; `CoupledField` gains labels.
4. `SpacePartition` from `FunctionSpace.partition_by`; block patterns DERIVED from tensor factors (every bulk leaf) and DECLARED per axis on the 4 leaves without them (streaming, leakage, the two boundary operators); `BlockRole` (117 lines, 29 files) and `SystemRole` (53, 10) retire.
5. The trace is the trace summand; `angular_trace` and `scalar_trace` become one role. The docs carry the loss-row reading (bulk → trace is `[A_tb | A_tt(out ← in)]`).
6. `TransportMethod` (45 lines, 19 files) retires: `axes` to `MaterialMesh`; `BOUNDARY_OPERATOR_REGISTRY` to the realizer's admission table; `realize_boundary_law` into the trace frame plus the realizer; `bc` to the boundary term.
7. `[OPEN]`, not blocking: whether flux and source grade the space (`V → V*`), which would make a role the space an arrow is bound on (unit 7's investigation).
8. Related issues: #296 (reify the block-operator abstraction: absorbed), #391, #287, #403.

## Gates

The prototype's 77 equivalence gates against production on the SN slab, the diffusion slab and the SN sphere, each with its mis-ordering control; the label gate for isomorphic summands; each derived block pattern against the dense matrix with a mis-declared control (orientation needs ≥ 3 groups with one-way coupling); `couple(split(A)) == A`.

## Done-when (hypothesis)

`FullFieldSpace`/`CoupledSpace` are one class; `BlockRole`, `SystemRole`, `TransportMethod` grep to 0; the gates green with their first reds; the full `-O` suite and `-W` build clean.

## Opening obligations

Re-run the consumer census of both carriers and both enums (AST); read the prototype memos, `full_field.py`, `coupled_system.py`, `operator.py`'s role machinery; W3 with test-architect.

## Sizing

3–5 sessions `[R]`.
