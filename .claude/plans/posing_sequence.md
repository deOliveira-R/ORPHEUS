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

## The review's findings (2026-09-25; memos `scratch/posing_sequence/attacker_memo.md` and `elegance_memo.md`, probes beside them)

Measured facts `[M]`:
- The k choice is fixed in the record, one level below the pencil: `WithinGroupSystem.loss` and `.production` are built for k (`sn/coupled_system.py:420-480`) and `hub.pencil` is made from them (`sn/problem.py:1181`; the homogeneous hub the same, `homogeneous/solver.py:336-353`). No c or α pencil can be built from that record.
- The Strategy does not read the pencil: `Splitting.from_schedule` reads `system.factors` (`sn/splitting.py:294`).
- The loss's signed terms are written twice: the expression `LC − S − N2N − B_a` and the carrying grid (`sn/coupled_system.py:647, 730`), and the `LossTerm(op, ±1)` tuples (`sn/splitting.py:292-325`), kept equal only by the law test.
- The subcritical multiplying source bypasses both the pencil and the splitting: fission is appended by hand as an extra gain (`sn/solver.py:3613`), and the splitting certifies itself against `at(0)` while the posing is `at(1)`.
- `ALPHA_MAP` is wrong: with `at(σ) = A − σM` and `M = 1/v`, μ = −α, but the map returns −1/μ (a 1×1 probe: true α −1.5, the map −0.667). No hub mints it; four pages cite it. Confirmed by the orchestrator from `numerics/posing.py:80-84`.
- Uniform dilation of the geometry is affine: the critical dilation factor is the eigenvalue of the pencil `(L, −(C − S − N_2n − F))` (probe: 4.111390152 both ways). Varying only the outer extent is not affine.

The converged answer to the user's attacks `[R]`:
1. **Attack (1) lands.** The Problem should end at its anchored term set: the signed operators `T_p` of one balance, and the physical point (the anchor, where `T = L + C − S − N_2n − F − B`; a term present in the data but absent at steady state, such as `1/v`, has coefficient 0 there), plus the source when one is posed. The pencil is not the Problem's last step; it is a view the QUESTION derives (`OperatorPencil.along(terms, carriers)`), and `hub.pencil` survives as the k question's view.
2. **Attack (2) lands, and the middle stage is the question itself.** A question is a curve through the anchor (which terms carry the parameter) plus what is asked of the resolvent on it: a pole (an eigenvalue), the value at a point (a fixed source, noise at ω, a Laplace variable), a path, a contour, or a root over a family of Problems (a non-affine critical parameter). It is a VALUE, which is why the reference specification can store it. Its laws: the family passes through the physical operator; the eigenvalue scale is derived from the carrier set; the adjoint of a pole needs no datum, of a source one detector; a rational family (delayed α) becomes affine on an enlarged Problem (neutrons plus precursors); one question's answer can anchor another (noise at the k-normalised point). The line between question and Strategy is a gate: a Strategy choice may not change the answer. This reverses two landed rulings: "α and noise are one Problem, two Strategies" (`consumers_step2_design.md` §3.5) and F11 (the fission-suppressed `at(0)` as a Strategy point): those points determine the answer, so they belong to the question.
3. **Attack (3) lands with one correction.** Triangularity is a property of a combination of operators under an order, never of one operator (the retired `creates_sweep_cycle` flag failed for this reason, `sweep_acyclicity.py:24-31`). Three declarations suffice: D1, the coupling pattern of each operator per axis (on the kernel, the scheme or the boundary law); D2, the rank of the closing edges per strongly connected component (close or split, report v2 §I.6); D3, the sign structure that lets a comparison theorem rank two splittings (Varga's regular splittings; the step scheme satisfies it on 384 of 384 probe configurations, diamond difference on 12 of 384). The lower and upper parts then follow from a partition and an order (generalising `SNBoundaryOperator.split`); Jacobi and Gauss–Seidel are not read off the explicit side, they are the two extremes of the implicit set (empty order and full topological order). `BlockRole` is a lossy seed of D1 (its join labels `C + B` as full although the sum is block-diagonal).
4. **The sequence:** MaterialMesh; the method's discretisation (spaces, per-axis partitions); the factors with their D1 to D3 declarations; the anchored term set (the Problem's end); the question (the family, the pencil as a view, the posing); the Strategy recipe; instantiation at the question's points; lowering and drivers; the Solution. Report v2 Part IV had this order; `consumers_step2_design.md` §3.1 merged the fourth and fifth steps, which is the inversion that landed.

## Candidate ontologies

- **The anchored term set, then the question** (the converged answer above) `[R]`.
- The landed shape: the pencil as the Problem's terminal object, the question a spectral map over it (`consumers_step2_design.md` shape (B)). Kept as the k question's view.

## Refuted candidates

- **The k pencil as the Problem's end:** it fixes the k choice in the record, and the Strategy reaches around it to the factors.
- **A `transfer` verb that changes the question:** `1/v` is not a term of the physical operator, so no transfer produces α; and the pencil's law `at(1) = T` holds for every transfer, so a c question that forgets a term (the ray emission, `sn/coupled_system.py:690`) passes silently. A carrier SET checked against the term set is complete by construction.
- **A separate stage that performs transfers:** its only content is the question value.
- **A new generic piece-set type in numerics:** it would be a third spelling of the terms; the operator sum should expose its terms.
- **A closed enum of carrier bundles:** `EMISSION` names a bundle; the building block is a set of channels.
- **Every point of the parameter is a Strategy choice:** refuted for the points that determine the answer (noise, fission suppression, Laplace).
- **"Operator X breaks lower triangularity":** triangularity belongs to a combination under an order.
- **The iteration type read from the explicit side:** Jacobi and Gauss–Seidel are the extremes of the implicit set.
- **A generic `Partition` type built now:** one bespoke instance (octant groups); energy Gauss–Seidel (#273) would be the second.
- **`linearize()` on `OperatorPencil`:** a degree-1 object has nothing to linearise; the enlarged Problem does.
- **"A critical size is always non-affine":** uniform dilation is affine; the outer extent alone is not.
- **"Fission is dense in space":** it is local to each cell; the fission matrix `A⁻¹F` is dense.

## The P1 seed, as the review proposes it `[R]`

- `Eigen(carriers)`: a set of channels (`FISSION`, `SCATTERING`, `N2N`; `TIME` when α is posed), with `EMISSION` a named constant for the c set; the eigenvalue scale derived from the set.
- `FixedSource`: the source (or, adjoint, the detector) and a set of SUPPRESSED channels, empty by default; today's σ = 0 is `{FISSION}` suppressed. It stores no σ.
- `CriticalParameter` for the outer extent stays separate (non-affine).
None of it needs operators, and none of it is undone by the campaign that builds the rest.

## Rulings ledger

- 2026-09-25: the review is ordered before any ruling on the synthesis (the user: "Let's take option 2 and attack it first").

## The P1 dependency

P1 of `reference_cache.md` needs only the question's vocabulary in the specification (step 7 of its carve). Steps 1 to 6 of P1 do not depend on this plan and may proceed.
