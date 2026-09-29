---
name: ordering-cut-interface-augmentation-frames
description: Verdicts of the 2026-09-29 structural attack on "Ordering as the primitive" (proto4 System/Ordering/Split/expose/reduce) and on the pipeline-of-morphisms second design — the spectrum of a cut lives on its INTERFACE (rank N), K ≥ 0 is the CW predicate (not regularity), Krylov runs on the interface, DSA is a THIRD-summand augmentation (not a hole), cold-started inexact inner solves converge to a wrong fixed point, the stored order is a chain, topology relocates cycles (1 ∈ spec K is the guard), expose is a section of reduce.
metadata:
  type: project
---

# Ordering / cut / interface / augmentation (attack memo
`scratch/posing_sequence/open_items/ordering_attack/structure/memo.md`, 8 probes beside it;
re-ground every number against the live tree before acting)

**Why:** the user ruled proto4's architecture (`system.order_along(partition, cut).split().solve(q)`,
`Piece` a `System`, `(M, N)` derived, `expose` as REDUCE's inverse) and owed it an attack before any
carve; the coordinator added a second design (System → System morphisms with a topological step).

**How to apply (durable verdicts, each with its number):**
- **The interface frame.** For any cut, `K = (M⁻¹N)[c, c]` (`c` = columns of `N`) or `K_r = W M⁻¹ U`
  (`N = U W`, `rank N`) carries the nonzero spectrum of `M⁻¹N` (`[M]` ≤ 3e-9 on 8 cuts; `rank N` = 8 /
  8 / 24 / 88 against `n` = 160 … 800). `spec(M⁻¹N) = λ/(1+λ)` over `spec(A⁻¹N)` (Möbius, ≤ 1e-11).
  **CW predicate = `K ≥ 0`, not regularity**: task 4 and task 2 have `M⁻¹ ≱ 0` (60 / 120 negatives)
  and `K ≥ 0` (brackets hold). Signed `K` (DD τ > 2): report the interface eigenvalues or the observed
  ratio; `ρ(|K|) = 1.041` on a converging mesh is the false alarm. **GMRES on the interface**: 5 / 5 /
  17 iterations vs 12 / 8 / 255 passes. **Cut ranking** = `ρ(K)` at `rank N` substitutions; equal `|B| =
  1` weights give ρ 0.964776 / 0.965731 (weighted FAS refuted). Over-cutting a simple cycle at `k`
  places: `ρ = g^{1/k}` (√g measured). `expose` INFLATES the interface to `n` (160 vs 24; 40 vs 14):
  expose `range T`, which is the CP moment space (SI's interface matrix IS the CP matrix of the resolvent).
- **The stored order is a chain.** `Ordering.components` is a linear extension; the DAG lives in
  `forward`. Task 4: depth 20 × width 8 (8 independent reflective chains); Jacobi = one level of 10 yet a
  sequence of 10; 2-D step scheme matches the grid-poset (depth `nx+ny−1`, width `4·min`). Production's
  `SweepSchedule.groups` grades at OCTANT granularity. Fix: store the DAG + levels; `sequence=` = a linear
  extension.
- **Cold-started inexact inner solves converge to the wrong fixed point** `(B⁻¹ + A − D) x = q`:
  m = 1 cold → relative error 1.0 accepted by the increment stop; tol 1e-2/1e-4/1e-8 → error = the
  tolerance. Warm (correction from the previous outer iterate) = the two-stage iteration, right fixed
  point, ρ predicted to 6 digits (0.9098 / 0.8280 / 0.6873 at m = 1 / 2 / 4). proto3's `nest_matrix` was
  warm; proto4's `inner=(piece, rhs)` callback is cold — the contract has no slot for state.
- **DSA is INSIDE the object**: third summand δ with `A_c δ − R T ψ + R s = 0` (`A_c = R A P`), source
  `s = T(ψ + Pδ)`, order (flux, coarse, source), cuts {flux←source, coarse←source}; the split's spectrum =
  `(I − P A_c⁻¹ R A) A_r⁻¹ T` to ≤ 2e-12 on 3 fixtures; ℓ≤1 Galerkin pair converges (ρ 0.502 / 0.185 /
  0.748 vs SI 0.841 / 0.849 / 0.918), P0 Galerkin diverges (no leakage term: ρ 13 – 39). Not regular
  (`N` has −0.36). The 2-summand fibre of `reduce` spells only SI (the flux iteration is `A_r⁻¹T` at
  every fibre point: `ρ(cut n←c) = ρ(I − A_nn⁻¹ S_nn)` to 12 digits).
- **`expose`/`reduce`** is section/retraction: `reduce∘expose = id` by theorem (`E = I`), `expose∘reduce`
  undefined (reduce forgets terms and 3 blocks: `KeyError`). `.H` commutes under ANY source metric
  (1.1e-17 for `G`, `G⁻¹`, random diagonal, random dense SPD; wrong-sign control 8e-2): the source
  metric cancels in the Schur complement. Regularity transfers as weak-regularity.
- **Mixed nest** (term split AND axis split at one level) = `expose` + cut `{(source_g ← flux_g') : g' ≥ g}`
  on `group × member`: ρ and step equal proto3's `of_rule` nest (m = 1) to 7e-16.
- **Declared-pattern blindness**: with graph AND classification read from a declaration missing the
  (thermal ← fast) block (512 entries), sequence (1,0) is accepted, 1 pass, error 0.37, residual 6e-2, no
  exception. A residual gate on `passes == 1` is mandatory in matrix-free production.
- **Topology (periodic / torus / cover)**: cycles RELOCATE (4 SCCs of 12 → 4 of 10; reflective 2 of 24 →
  cover 2 of 20); interface spectrum = loop gains (torus `Π 1/(1+τ)`; cover `g²`); closure `(I − K)y = b`
  exact to 2e-16 (`|c| + 2` substitutions; Sherman–Morrison at `|c| = 1`); images to 4.4e-16 (random
  source; control 0.32); with scattering the SCC is not simple and the interface is 14 = 10 moments + 4
  wrap, still exact; **singular iff `1 ∈ spec K`** (c → 1: 0.7113 / 0.99618 / 1.00000000; void → 1):
  proto2's gain-1 guard is the `|c| = 1` case. REFUTED (mine): "topology restores regularity" — the
  step-scheme trace form is regular; proto2's non-regularity is DD's sign chain (a SCHEME trait).
- **Second design (pipeline)**: the morphisms are typed by LEVEL (1 manifold: quotient/cover = pullback
  along the covering map; 2 summands: reduce ⇄ its sections (`augment` family: expose, Galerkin
  correction), couple ⇄ split; 3 partition: condense = `order_along` with no cut, cut/order, split →
  (pieces, forward, LAGGED COUPLING WITH ITS RANGE = the interface); 4 drivers). Composition order is
  FORCED: a level-k morphism changes what level k+1 reads. `condense` is not a third noun (the recursion
  covers fine-cut-then-coarse-order). Names: `Split` → `Splitting` (perfect match, the tree's word);
  `expose` invented → `augment` (0 hits, unspent).
- **Layers by design**: Krylov/Chebyshev/Anderson (`numerics/iteration.py`, over the interface
  iteration); eigenvalue loop (`numerics/eigenvalue.py` + `pencil.at`); nonlinear feedback (Newton driver).
  Outside: overlapping Schwarz / multisplittings (not partitions; `[HYPOTHESIS]` overlap = an augmentation
  with identity constraints), non-triangular `M`.

Related: [[coproduct-direct-sum-two-level-frames]] (REDUCE as a rank guard: now `1 ∈ spec K`),
[[power-iteration-vs-keigenvalue-morphism]] (the opaque resolvent consumes `Split.solve`),
[[dsa-saddle-point-mixed-fem-frames]] (`R A P` = Schur inheritance; here it is the third summand),
[[pullback-pair-frobenius-perron-frames]] (the topological step is the deck/pullback),
[[shift-ontology-taxonomy-frames]] (the Möbius chart `λ/(1+λ)` of a cut).
