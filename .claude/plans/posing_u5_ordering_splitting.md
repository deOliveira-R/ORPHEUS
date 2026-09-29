# A solve is an ordering of a system and the splitting it induces, with its contraction rate predicted (posing unit 5)

Status: SCHEDULED, not opened (2026-09-29). Issue #528; charter issue #522. Charter: `.claude/plans/posing_sequence.md`, "The work, ordered by dependence"; ruled content: "The ontology as it stands", "Layer 3: the solution — the structure level". Evidence: `scratch/posing_sequence/open_items/coproduct/proto2/`, `proto3/`, `proto4/` and `scratch/posing_sequence/open_items/ordering_attack/{structure,elegance}/`. Depends on unit 4.

## Goal, in the domain's terms

Every iteration the tree runs (the sweep, source iteration, DSA, boundary Gauss–Seidel, group ordering) is the same construction: order the system's unknowns (a partition and a partial order, with declared cuts), split at the order into members solved in sequence with the cut couplings lagged, and read the contraction rate from the small operator the cut leaves on its interface. The mistakes a driver can make today (a split without an order, an undeclared lag, a cut that leaves a cycle, a cold inexact inner solve, stopping on the increment) become unspellable or refused.

## Scope (ruled 2026-09-29)

1. **`Ordering`**: a `SpacePartition` plus a PARTIAL order of its blocks (the SCC condensation DAG after declared cuts, with levels; Jacobi the antichain, Gauss–Seidel a chain), DERIVED or DECLARED (a sequence with its induced cut); a cut that leaves its SCC alive is refused; the adjoint of an ordering is the reversed ordering with the transposed cut; the graph step is `order_along` without a cut; the directed/reciprocal kind reported per axis, the DAG stored only where directed.
2. **`Splitting`** (the output of `split()`): the members in order, each a `System`, and the lagged coupling with its INTERFACE; `(M, N)` derived, never stored or exposed, `M` never assembled; the interface operator `K`; the predicted contraction rate `ρ(K)` with its Collatz–Wielandt bracket when `K ≥ 0`, reported otherwise; Krylov on the interface; singular iff `1 ∈ spec K` (the trace-elimination test).
3. **The inner contract**: warm correction from the previous outer iterate with an inner count; a residual gate on every solve; an `IterationRecord` with children (#340's contract).
4. **`augment` and `eliminate`** as far as the structure level needs them (eliminate carries its lift and back-substitution and keeps its terms); source iteration as augment (the scattering term's image) then an ordering over (flux, source) with the back cut; DSA as the Galerkin-correction augment (a coarse summand, ordered (flux, coarse, source)); the SN sweep as the triangular ordering over (inflow, bulk in flight order, outflow).
5. **Energy Gauss–Seidel** restricted to the energy SCCs in production SN (exact on downscatter-only problems); boundary Gauss–Seidel as a declared ordering; a 2-D boundary Gauss–Seidel gate BEFORE `Splitting.from_schedule` retires.
6. **Predecessors re-spelled or retired**: `Splitting.from_schedule`'s public halves, the octant schedule's hand ordering, `SourceIteration` and `KrylovAcceleration` (their recipe to the ordering, their instance at a point to the binding), `power_iteration` stays the binding's; the elegance attack's predecessor map (`ordering_attack/elegance/memo.md`) is the census.
7. **Out of this unit by design**: Chebyshev/Anderson (`numerics/iteration.py`), the eigenvalue loop, nonlinear feedback, overlapping Schwarz.
8. **Names**: `Ordering`, `Splitting`, `Member`, `augment`, `eliminate`, the tree's axis names (never string tags); typed inner strategies per member are a later item (today an `inner=` callback).
9. Related issues: #273 (energy-group Gauss–Seidel: absorbed), #324 (SCC trace digraph into the `B_lower`/`B_upper` splitting: absorbed), #343 (the octant order as a rate lever: answered by the predicted `ρ(K)`), #315, #321, #469, #296.

## Gates (the charter's structure-level list)

`couple(split(A)) == A`; `M − N == A` for every splitting; the SCC singularity prediction against the SVD nullity; `K`'s spectrum against the full error-propagation operator's; the Collatz–Wielandt bracket where `K ≥ 0`, reported where not (a diamond-difference control); the reversed-ordering adjoint; the refusal of a cycle-leaving cut; the residual gate with the missing-block control; the Gauss–Seidel direction by block-triangularity; the DSA augment's spectrum against `(I − P A_c⁻¹ R A) A_r⁻¹ T`; the 2-D boundary Gauss–Seidel gate; production's source-iteration iterate reproduced (proto4: 1.7e-16).

## Done-when (hypothesis)

Every production iteration is spelled as an ordering and its splitting; `Splitting.from_schedule`'s halves are private or retired; the gates green with their first reds; #273 and #324 closed.

## Opening obligations

Re-read proto4 and the ordering attack; census the production drivers (AST); read the SN iteration theory pages and #340's contract; W3 with test-architect.

## Sizing

5–8 sessions `[R]`.
