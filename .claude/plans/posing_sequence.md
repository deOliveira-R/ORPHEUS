# The posing sequence after the material mesh — a living plan (ontology being searched)

Opened 2026-09-25 from #405 step 2, phase P1 (`.claude/plans/reference_cache.md`, "P1, fourth exchange" and "P1, fifth exchange"), where the reference specification's stored question (`Eigen(k)`, `Eigen(c)`, ...) raised what an eigenvalue question IS in the operator algebra. The ontology is not settled; implementation waits for the user's ruling that this plan is polished (`plan-authoring` §0). Markers: `[M]` measured, `[R]` reasoned, `[HYPOTHESIS]`, `[REFUTED YYYY-MM-DD]`.

## The goal, in the domain's terms

From the operators a method has posed on its space to a computed answer, name every step as an object in its right place, so that: the question asked (k, c, α, noise, a fixed source, a subcritical source) is data, never a code path; every legitimate way of computing an answer to one question (a shift, a resolvent, a splitting along any axis, Jacobi or Gauss–Seidel, a schedule, a lowering) is a value the machinery can derive or check from the operators' own declared structure; and no step carries data that belongs to another.

## What is already ruled (do not regress)

- The loss and gain operators are built as a PAIR: the assignment of pieces to sides is the primitive, the sides are derived, the only mutation is a transfer, a split is a refinement (`operator_strategy_realization_campaign.md` §1.0, the user, 2026-07-29).
- Loss, `M` and `N` take the same partition (shape symmetry, the same plan's P3, the user, 2026-07-29).
- Building the pencil is the Problem's last step; manipulating it (resolvent, partitioning) is the Solution side; the source-driven problem builds something else as its last step (`consumers_step2_design.md` §1, the user, 2026-09-13). Shape (B) landed: `OperatorPencil` (layer 1), `EigenPosing` and `SourcePosing` (layer 2), four cells, the `(M and q)` cell a composition (`consumers_step2_design.md` §3.5).
- The partition is per axis and Space-owned; the block digraph and its order are computed from operator and partition; the schedule is the splitting's (`orpheus-operator-machinery-report-v2.md` §I.10).
- The resolvent is not a new type: `R(σ) = pencil.at(σ).inverse()` (the same report, §I.4).
- The Solution side changes the trajectory of solving a problem, never the problem (restated by the user, 2026-09-25, attack 3 below).

## The first-pass synthesis under attack `[R]` (the orchestrator, 2026-09-25; full text in `reference_cache.md` "P1, fifth exchange")

1. The pencil and the splitting are two labellings of one structure: the pieces of one operator sum `T = L + C − S − N_2n − F − B`. A pencil labels each piece "left" or "carries the parameter" (k: fission carries; c: all collision emission carries); a splitting labels each piece of `at(σ)` implicit or explicit. Both admit `transfer` and `refine`.
2. The owner of a transform follows from what it preserves: the splitting's transfer preserves the operator (`M − N = A`), so it is the Strategy's; the pencil's transfer preserves only `at(1)`, so it changes the question and is the Problem's. Everything that keeps the question (a shift, a partition, the implicit/explicit labelling, the schedule, the lowering) consumes the pencil and is never one of its methods.
3. Five kinds of object: the factors; one generic labelled piece set with `transfer` and `refine`; the pencil built from a labelling; a per-axis partition; the Strategy values composing them.
4. The specification names the carrier in data-layer vocabulary (`FISSION`, `EMISSION`); `CriticalParameter` keeps the non-affine parameters.

## The user's attacks (2026-09-25, verbatim)

*"(1) When we designed the EigenPencil we thought that the Pencil was the most convenient way to express the system at the end of Problem. Maybe that is true, of maybe what we call today Problem should just define the Operators and that should be fed to a Problem Object that defines the shape of the Pencil. Or maybe we we right and the Pencil is the right way to end the Problem (since there all operators are defined and we have recognized the problem as an eigenvalue problem of SOME sort, leaving the question of "which sort of eigenproblem this is to the next step). (2) We might need a step between Problem and Solution to perform SOME transfers, more specifically the ones that change the problem in a fundamental way, like k eigenvalue, alpha, neutron noise, etc. (3) We have recognized the Solution side as only changing the trajectory of solving the problem, not the problem itself, and I think this reading was ontologically correct. But this means, among other things, that from the terms on RHS we should be capable of knowing HOW the problem can be split. For example, maybe BoundaryOperator or ScatteringOperator break lower triangularity. We need both the capability of Splitting the Operators into upper-lower pair and decide which one is implicit and explicit. Finally, from what is on the explicit side, we should be capable of figuring out HOW we can split. Jacobi? G-S? etc. It's possible that the ontology of the sequence is not complete yet. (Like the difficulty we had to reach the ontology of layering we need to emerge regarding the sequence of steps to reach the MaterialMesh and then turn that into Problem."*

## The review (W5, dispatched 2026-09-25)

The cross-domain-attacker and the elegance-enforcer attack the synthesis and the user's three attacks, first pass adversarial, second pass a separate re-evaluation. Their memos land in `scratch/posing_sequence/`.

## Candidate ontologies

(Filled from the review.)

## Refuted candidates

(Each with its structural reason.)

## Rulings ledger

- 2026-09-25: the review is ordered before any ruling on the synthesis (the user: "Let's take option 2 and attack it first").

## The P1 dependency

P1 of `reference_cache.md` needs only the question's vocabulary in the specification (step 7 of its carve). Steps 1 to 6 of P1 do not depend on this plan and may proceed.
