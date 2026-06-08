---
name: issue-208-operator-algebra-frames
description: Frame attack on the Wave O / issue #208 operator-algebra framing (V=V_bulk⊕V_boundary, 2×2 block operators, capability bundle, adjoint-for-free, cross-method SN/MoC/CP/diffusion). CONFIRMED most-native = dagger inverse biproduct category with a CO-EQUAL non-trivial boundary metric G + method-specific solve-refinement mixins. Strong payoff frames: biproduct (block matrices are a THEOREM), dagger+G-metric (two asymmetries are ONE structure; †=G⁻¹AᵀG; Gate 1.3 is metric-blind), inverse-category (solve is a partial-inverse morphism with solve∘apply=id oracle, not a string tag), triangular-solve-as-refinement (cross-method layering), transport-resolvent (predicts diffusion exception). REFUTED: operad/PROP (no extra payoff over dagger biproduct), homology (B∘B≠0), spectral-on-typing. NEW smell candidate: metric-blindness. Merged the degraded full-access re-run 2026-06-08.
metadata:
  type: project
---

# Issue #208 (Wave O) operator-algebra — foreign-frame attack

Branch `refactor/field-role-typing`, worktree `.claude/worktrees/field-role-typing`.
Candidate framing: `V = V_bulk ⊕ V_boundary`; operators = 2×2 block matrices
`[[A_bb,A_bs],[A_sb,A_ss]]`; operator = typed domain/codomain + capability bundle
`{apply, solve, apply_transpose, solve_transpose}`; only streaming `L` is full,
C/S/F bulk-only, B boundary-only; adjoint "for free" via wired closure laws + sweep
order-reversal; boundary inner product carries `|Ω·n|·w_n`.

Ground truth read: `.claude/agent-memory/explorer/issue_208_operator_algebra_surface.md`
(Area 1 = L1 algebra, closure laws fully wired, gap is only the streaming-leaf adjoint;
Area 5 = Gate 1.3 tests `.bulk.values` Euclidean only; Area 6 = TraceSpace
inner_product_weights=None today but omega_dot_n + layout already stored).

> MERGE NOTE (2026-06-08 hygiene): absorbs the former
> `issue_208_operator_algebra_frames_full_access.md` (a degraded re-run — Read/Grep/
> Nexus all DENIED, code facts UNVERIFIED). That run added two durable refinements,
> now folded in here against the verified code: (a) the block-metric G promoted to
> CO-EQUAL with the biproduct (#1b below), and (b) the **metric-blindness** smell
> (Part C candidate, bottom). The full-access run's code-level claims were
> superseded by THIS file's verified file:line facts. Standing permission-degradation
> flag also folded to the cross-agent note, not the library.

## CONFIRMED most-native framing (RANKED)

**#1 — Dagger inverse biproduct category** (the unifying structure).
The candidate's TWO architectural asymmetries are ONE named structure:
- `V_bulk ⊕ V_boundary` is a **biproduct** ⟹ morphisms between biproducts ARE
  matrices of morphisms (Mac Lane VIII.2 additive-category matrix calculus). The
  2×2 block algebra is a THEOREM, not a modeling choice. The implicit-zero boundary
  is the projection identity `π_s ι_b = 0` — deletes the `codomain_zero` hook
  (operator.py:822) + the per-leaf `BoundaryFlux.zeros_on` returns. N-way ⊕
  (kinetics flux⊕precursors, multiphysics) is free: Mat of an additive category
  is closed under finite biproducts. Protocol classification (Bulk/Full/Boundary)
  = which ι/π pair the morphism factors through (a factorization, not a tag).
- `.H` with `(A+B).H=A.H+B.H`, `(A∘B).H=B.H∘A.H`, `(αA).H=ᾱA.H` IS a **dagger
  functor** (the three wired closure laws ARE the functor axioms; "adjoint for free"
  = "define † on generators, extend by functoriality"). The `|Ω·n|·w_n` metric is
  what makes † the HILBERT adjoint not the naive transpose.
- A **dagger biproduct category** carries both at once. The `solve` partial inverse
  lives in the **inverse-category** layer (Cockett–Lack restriction category): `solve`
  is a partial-isomorphism morphism `f°` with `f°f=id` on domain-of-definition, NOT a
  peer string tag. The right base = **dagger inverse biproduct category** with a
  non-trivial metric on V_boundary.

**#1b — Non-trivial block metric G (CO-EQUAL with the biproduct, not subordinate).**
`G = diag(I_bulk, |Ω·n|·w_n on V_boundary)`. The dagger `†` MUST be the G-adjoint
`G⁻¹AᵀG`, NOT the bare transpose. Without G the dagger is ambiguous and
"adjoint-for-free" is unverifiable — this is the ROOT of the Gate 1.3
metric-blindness trap (deliverable #1 below). The `solve∘apply=id` partial-inverse
law (`f°∘f=f̄` on domain-of-definition) is the free correctness oracle of the
inverse-category layer: a `General` morphism (C/S/F, no inverse) hitting `.solve`
is a TYPE error, not a runtime `MissingCapability` sentinel.

**#5 — Restriction/inverse category** answers Q2 (where solve lives): operator
declares a MORPHISM CLASS `Iso`/`PartialIso`/`General`; the class IMPLIES the
capability set. Illegal-states-unrepresentable replaces the 4-tag powerset +
MissingCapability-at-composition. For a partial isometry `solve_transpose=(solve).H`
— a relation the flat-tag design cannot express.

**Cross-method layering verdict (Q3):** apply/solve/adjoint + biproduct + dagger =
METHOD-UNIVERSAL base (covers SN/MoC/CP/diffusion). Causal-ordering/triangularity is
an SN/MoC-ONLY **refinement mixin** (`TriangularSolveMixin`: solve=forward-subst along
method-supplied causal ordering; adjoint=back-subst on reversed ordering). CP gets
`DenseFactorizationSolveMixin` (solve = dense LU/Cholesky inverse), diffusion gets
`SelfAdjointSolveMixin`/`EllipticSolveMixin` (solve=CG/MG, solve==solve_transpose).
The candidate framing's RISK: it lists causal-ordering at the SAME level as
apply/solve/adjoint — if the base bakes in ordering, CP (dense) and diffusion
(elliptic) cannot inhabit it. Triangularity MUST be a refinement.

**Transport-resolvent (Feynman–Kac) backs the layering analytically:** SN/MoC/CP/MC
`solve` = four QUADRATURES of one object, the resolvent `(Ω·∇+Σ_t)⁻¹` (Peierls kernel);
adjoint solve = backward semigroup (`Ω→−Ω` = path reversal). Diffusion is the EXCEPTION
— P1/asymptotic LIMIT of the resolvent, not a quadrature — which is exactly why its
solve is elliptic-self-adjoint while the others are characteristic-triangular. One
principle predicts BOTH the layering split AND the diffusion exception.

**Sheaf/fiber-bundle (MoC axis):** B (R·G boundary law) is a GLUING COCYCLE, not a
generic A_ss block — periodic BC = base-manifold QUOTIENT (true S¹ ring), reflective =
deck transformation `Ω→R(Ω)`. O.4 BC-extraction = "lift the gluing cocycle out of the
local solve into global-section assembly." Shares the gluing abstraction with the
committed MoC fiber-bundle frame (`[[project_moc_structure]]`) even though SN/MoC do
NOT share the SweepGraph.

## HIGH-VALUE FREE DELIVERABLES (first tests that discriminate)

1. **Gate 1.3 is metric-blind (verification gap, follows from dagger frame).** It tests
   `.bulk.values` Euclidean only (explorer Area 5). The moment an off-diagonal block
   A_bs is non-trivial, the correct adjoint is the Hilbert adjoint w.r.t. the
   `|Ω·n|·w_n` boundary metric. Test: run reciprocity twice (weights=None vs
   `|omega_dot_n|⊙w_n`) with a non-trivial boundary block; pass-Euclidean/fail-metric
   proves the metric is load-bearing and current Gate 1.3 insufficient. PAIR WITH the
   negative control the test-architect already flagged (wave_o_operator_typing.md §verif).
2. **Involution round-trip `L.H.H == L`** (dagger axiom f††=f) — free invariant test
   the moment StreamingOperator.apply_transpose ships.
3. **Third-summand test** for the biproduct frame: add a trivial 3rd ⊕ slot and check
   OperatorSum.apply needs NO new branch (projection relations handle it). If it needs an
   `if`, the code is block-matrix-by-convention not biproduct. Discriminates whether the
   biproduct is latent or must be built — directly informs N-way ⊕ future-proofing.
4. **Capability-set lattice check:** enumerate observed capability sets ({APPLY},
   {APPLY,SOLVE}, {APPLY,SOLVE,APPLY_TRANSPOSE}, {APPLY,APPLY_TRANSPOSE} per explorer
   Area 1) and check they map 1:1 onto `{General⊂PartialIso⊂Iso}×{non-dagger⊂dagger}`.
   1:1 ⟹ inverse-category frame fits, tags derivable from class. An orphan combo
   ({SOLVE} w/o {APPLY}) would refute it.
5. **k/adjoint-eigenvalue overlay:** forward generalized EV `(L+C)ψ=(1/k)Fψ` vs adjoint
   EV via wired `.H` must give same k to tol — strong argument the algebra MUST carry
   the metric for the adjoint flux to be correct; sensitivity `δk=⟨φ,δLψ⟩/⟨φ,Fψ⟩`
   becomes one inner product.

## REFUTED / low-signal for this problem class
- **Operad / PROP** — dagger biproduct categories ARE the relevant compact-structure
  PROP; naming it "PROP" adds zero concrete payoff over "dagger biproduct category".
- **Homology / chain complex** — tempting via the word "boundary" but `B∘B ≠ 0`
  (two reflections compose to a non-trivial map); trace γ + extension are a dagger
  adjoint pair, not a differential. No `∂²=0`.
- **Spectral / multiplication-operator ON THE TYPING QUESTION** — fired strongly for
  the µ-kernel (`[[phase5_continuous_mu_frames]]`) but the block-typing question is
  "right type for a block-structured morphism," not "diagonalize this operator." Spectral
  re-enters only at the eigenvalue overlay (pollination), not the core framing.

## ELEGANCE SMELLS observed (candidate #16 promotion material)
- "Natural 2×2 structure" / "adjoint for free" = analogy-language for a THEOREM
  (biproduct + dagger functor). Structure described as "natural" when forced ⟹ native
  frame not yet named. (Same shape as Smell #15: missing variational principle.)
- `codomain_zero` callable + implicit-zero shims = correction term standing in for the
  biproduct projection identity `π_s ι_b = 0`.
- Capability string-tag set checked at composition (MissingCapability) = stringly-typed
  dispatch standing in for a typed morphism-class statement.
- TWO geometry-non-uniform BC-absorption paths (1-D seed/reflect vs 2-D cell-fill, 2-D
  drops the boundary residual entirely) = structurally distinct paths to one operator +
  a correctness asymmetry (Smell #16-candidate, already tracked in
  `[[issue_168_phase_c_sweep_frame]]`).
- **Metric-blindness** (NEW, from the full-access re-run — Part C candidate): a
  metric/inner-product-dependent property (adjointness, orthogonality, sensitivity)
  checked under the WRONG inner product, structurally unable to detect a wrong
  weighted block. TELL: the adjoint/reciprocity test reads `.bulk.values` Euclidean
  only (Gate 1.3, explorer Area 5) while a weighted block (V_boundary, `|Ω·n|·w_n`)
  carries a non-trivial measure. Discriminating test: build `L†` wrong in the
  boundary block (use `Aᵀ_ss` not `G⁻¹Aᵀ_ss G`); the Euclidean gate PASSES, the
  metric gate FAILS. Promote to Part C after the G-adjoint Gate 1.3 fix ships +
  the discriminating test is demonstrated.

## Cross-refs
`[[trajectory_resolvent_foreign_frames]]` (rank-N=MPO bond dim; ⊗ is orthogonal to this ⊕),
`[[issue_168_phase_c_sweep_frame]]` (sweep-as-matvec + Smell #16),
`[[phase5_continuous_mu_frames]]` (resolvent/multiplication-operator precedent),
`[[project_moc_structure]]` (MoC fiber-bundle, shares the gluing abstraction).
