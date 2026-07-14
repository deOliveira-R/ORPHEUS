---
name: assembly-mode-functor-dag-frames
description: Frame verdict on the SN/transport ASSEMBLY mode (the 3rd consumption of the ONE closure algebra — solve/apply/assemble; campaign Phase 2b, Cartesian slab+2D, curvilinear OUT). THREE verdicts. Q1 assembled≡probed IS FUNCTORIALITY (not a numeric coincidence): as_matrix is the additive-monoidal functor F:Op→Mat (F(A+B)=F(A)+F(B), F(A@B)=F(A)F(B), F(αA)=αF(A), F(A⊗B)=F(A)⊗F(B)); assemble = the SAME functor into scipy.sparse; the flat C-ravel order IS the monoidal/Kronecker-factor structure (ordinate⊗group⊗space for SN, group⊗space for diffusion — per-carrier, pinned by to_flat). SumOfTensorProductsOperator (L=Σ_a D_a⊗Ω_a⊗I_g) already names it. Q2 per-octant (L+C) block IS the weighted incidence matrix of the sweep DAG: triangular ⟺ acyclic walk order (levels = the permutation P; Pᵀ(L+C)P strictly lower-tri); #282 pole-seed lag = a back-edge, detectable as an above-diagonal entry on the permuted sparsity. Assemble per-ordinate blocks + lift (block-diag over ordinates + separate scattering I_space⊗Λ_moment), NOT monolith; DSA R·A·P = a clean triple product on the ANGLE factor. Q3 emission earns NO new type (scipy COO/CSR is the carrier; the sum-concat/scale/compose laws are INHERITED from the operator algebra via the sparse-Mat functor — a COO-builder-with-laws is a twin of OperatorSum). Local-to-global map = the EXISTING carrier ravel (FullField.to_flat / FlattenedOperator template) — SHARE it (Smell #16 shape-2), do NOT re-derive; reify as a type only when an UNSTRUCTURED (MoC/DG connectivity-table) consumer lands. Branch-verified on refactor/pyright-burndown worktree 2026-07-04.
metadata:
  type: project
---

# Assembly mode — structural frame verdict (attacker Phase 2b)

Read DIRECTLY on the worktree (L-005). Ground facts: `as_matrix`
(`numerics/operator.py:751`) probes column j = ravel(A·e_j) in C-order,
docstring "functor OUT of the operator category (Op→Mat)"; `FlattenedOperator`
(`numerics/flat_operator.py`) = `flat∘A∘unflat`, template pins the ravel;
`FullField.to_flat` = `concat(bulk.values.ravel(), boundary.values)` (direct
sum, C-order bulk); `assemble_ubld` (`sn/spatial/_ubld.py`) ALREADY assembles
per-cell `2^d×2^d` blocks via batched Kronecker (`M=⊗ mass_1d`, `_GRAD_1D`,
`_FOUT_1D`, `_FIN_TRACE`); `SumOfTensorProductsOperator` + `TensorProductOperator`
(`&`) already carry the Kronecker algebra AND the doc `L=Σ_a D_a⊗Ω_a⊗I_g`;
`SweepDependencyGraph.levels` = the anti-hyperplane topological order (cells
`i+j=const` independent). Diffusion `LeakageOperator` (`diffusion/operators.py`)
= tridiagonal FD, resolvent built by PROBING (`MatrixInverseOperator(FlattenedOperator(A,template))`).

## The dominant smell — #16 shape-4 (fires BEFORE the code)
`assemble` is the THIRD consumption of the ONE closure algebra already shared by
`solve` (sweep) and `apply` (matvec). Both already single-source the per-cell
coefficients (`CellBalanceTerms`, `assemble_ubld`, `D1ClosedForm`). The FIX is
mandatory-by-construction: the assembler MUST consume the SAME `assemble_ubld` /
`cell_balance_for_streaming` blocks — re-deriving the stencil in the emitter is
the twin. Also shape-1: probe (`as_matrix`) vs structural assemble are two paths
to one operator; the `assembled≡probed` gate turns the correctness CLAIM into a
theorem (Q1).

## Q1 — assembled ≡ probed IS FUNCTORIALITY (math-demands)
`F: Op→Mat` is an ADDITIVE-MONOIDAL functor: objects→flat dim (via `.shape`/ravel),
operators→matrices, with `F(A+B)=F(A)+F(B)`, `F(A@B)=F(A)·F(B)`, `F(αA)=αF(A)`,
`F(A⊗B)=F(A)⊗F(B)` (Kronecker). `as_matrix`=F into dense; `assemble`=F into sparse.
`assembled≡probed` is functoriality — a THEOREM — iff (1) each LEAF's assemble ==
its probe, (2) composites recurse via the homomorphism laws. The **flat C-ravel
order IS the monoidal structure** (the lexicographic tensor order of the axis
factors): SN `(N,ng,*spatial)` = ordinate⊗group⊗space; diffusion `(ng,nx)` =
group⊗space. Per-carrier, pinned by `to_flat` — NOT a hidden bug, but assembly
MUST emit in the carrier's own order or `F(A⊗B)=F(A)⊗F(B)` fails. RECOMMENDATION:
`assemble()` recurses over Sum/Product/Scaled/TensorProduct via scipy.sparse
add/matmul/scale/kron (leaves override); `as_matrix` on an assemble-able op
delegates to `assemble().toarray()`. Kronecker-per-factor emitters (Q1b) are
clean ONLY on uniform mesh+const-XS (translation-invariant `_GRAD_1D`); general
per-cell coefficients vary ⟹ per-cell-block COO scatter is the general leaf.

## Q2 — per-octant (L+C) = weighted incidence matrix of the sweep DAG (math-demands)
`(L+C)_oct` is lower-triangular in walk order (C diagonal; L strictly-lower +
self-term). `levels` IS the permutation P; `Pᵀ(L+C)P` strictly-lower-tri ⟺ the
walk order is acyclic. The whole within-group op is NOT triangular (octants have
different orders, scattering couples them) — only PER-OCTANT. So assemble
per-ordinate blocks + lift: block-diag over ordinates (streaming+collision, no
ordinate coupling) + scattering as a SEPARATE Kronecker term `I_space⊗Λ_moment`
(Funk-Hecke: dense-in-angle, diagonal-in-space). NOT the monolith. On uniform
mesh the per-octant block is block-Toeplitz-triangular. **#282 pole-seed lag =
a back-edge (feedback arc)**: the Carlson seed reads the most-inward ordinate's
outer value = a row referencing a not-yet-emitted column ⟹ an above-diagonal
entry on the permuted sparsity. Curvilinear angular redistribution ALSO couples
ordinates within a cell (Morel-Montry) ⟹ the block-diag-over-ordinate assumption
breaks — WHY curvilinear is OUT (only characterized). **DSA payoff**: because the
assembled A carries the space⊗angle factor structure, the moment reduction R·A·P
(R=moment analysis M, P=reconstruct, both on the ANGLE factor) commutes through
the spatial factor ⟹ per-block assembly makes the diffusion synthetic operator a
CLEAN TRIPLE PRODUCT on the angle block.

## Q3 — the emission earns NO new type (structure-demands)
Apply L-004. The (row,col,value) triple stream = scipy COO internals. The
candidate morphisms (sum-concat, scalar-scale, compose-with-diagonal, walk-order
permutation) are ALREADY the operator algebra (`OperatorSum`/`ScaledOperator`/
`OperatorProduct`/`PermutationOperator`). A "COO-builder with composition laws"
re-implements OperatorSum one layer down = a TWIN. So: scipy.sparse (COO→CSR) is
the carrier; the algebra is INHERITED via the sparse-Mat functor (Q1). The ONE
thin new leaf is a serialization carrier `SparseAssembledOperator(LinearOperator)`
wrapping the scipy matrix so `MatrixInverseOperator`/scipy consumers work — the
EXACT parallel of `FlattenedOperator`+`MatrixInverseOperator` (earned by the
scipy-consumer boundary, NOT a new algebra). The **local-to-global map** = the
carrier's EXISTING ravel (`FullField.to_flat` / `FlattenedOperator` template /
`SweepDependencyGraph.levels` global indices) — SHARE it (Smell #16 shape-2: one
DOF numbering re-derived in the assembler is the bridge = the un-named map).
Reify a `LocalToGlobalMap` TYPE only when an UNSTRUCTURED consumer (MoC per-ray,
DG connectivity table) lands — structured-grid ravel is a closed-form index
formula, one realization, zero applied non-id morphism (L-004 clause a fails).

## Discriminating first tests (all can RED)
- Q1 functoriality: monkeypatch every leaf `apply` to raise; `A.assemble()` MUST
  still succeed (reads coefficients, not applies) — a probe-fake REDs. THEN
  `assemble(A).toarray()` allclose `A.as_matrix()` @1e-13; mutate a leaf assemble
  to DROP the `F_in` inflow term ⟹ gate REDs.
- Q2 triangularity: `Pᵀ(L+C)_oct P` strictly-lower-tri (`triu(k=1)==0`) on slab;
  REVERSE the up/down face coupling ⟹ above-diagonal entries RED. Sphere: assert
  the per-octant block is NOT strictly-lower (positive #282 back-edge detection).
- Q3 no-new-type: `assemble(A+B)` equals `assemble(A) + assemble(B)` under scipy
  `+` AND `assemble(A).toarray()` allclose `FlattenedOperator(A,template).as_matrix()`
  (shared ravel); a group-major assembler on an ordinate-major carrier ⟹ index
  permutation mismatch REDs. Assert `type(A.assemble())` is scipy CSR / the thin
  wrapper, NOT a bespoke `StencilEmission`.

## Cross-method pollination
- **Diffusion FD = the SECOND assembly consumer** (unify-after-2 trigger): the
  `LeakageOperator` tridiagonal resolvent is built by PROBING today (O(n²) applies
  + dense LU) — assembly emits ~3n nonzeros direct (O(n)); first test = assemble
  the diffusion `A=L+C−S−B`, compare to the `FlattenedOperator.as_matrix()` probe.
- **FEM/DG**: assembly = pushforward along the local-to-global map; `assemble_ubld`
  ALREADY emits the local element matrices (`2^d×2^d`), the scatter (COO along the
  ravel) is the missing half. The scatter IS `coo_matrix((block.ravel(),(rows,cols)))`
  with scipy summing duplicates — no ORPHEUS type.

## Refuted (durable UNEXPLORED for this problem class)
- **Homology / chain complex (DEC)** — the face cochain `C¹=C¹_int⊕C¹_∂` is real
  (sweep_graph `_MovingFrontier` ι_*/ι* ALREADY names it) but assembly is a
  MATRIX-realization question, not a `∂²=0` one; `B∘B≠0` in transport (L-001). The
  DEC frame belongs to the trace algebra, not the assembler.
- **Block-circulant / FFT-diagonalization** — the sweep is TRIANGULAR not
  circulant (causal, not periodic); block-Toeplitz on uniform mesh but the pure
  Toeplitz→circulant embedding needs periodic BC AND buys nothing over the O(n)
  triangular solve. Fires only if a fully-periodic within-group operator (no
  sweep) ever needs a direct solve.
- **MPO / tensor-network truncation** — `SumOfTensorProductsOperator` is a FIXED
  rank-d sum (d spatial axes), no bond-dimension truncation knob (L-001).
- **Category theory (beyond the plain functor)** — the additive-monoidal functor
  F:Op→Mat IS the concrete categorical content and it is NAMED (as_matrix
  docstring); no further adjunction/Kan lever produces a test.

Cross-refs: [[flux_to_sourcesink_operator_contract_frames]] (the Op[Din,Cout]
contract F assembles), [[streaming_apply_transpose_frame]] (the (I−N) triangular
factorization Q2 assembles), [[field_role_typing_faceflux_frames]] (the C¹
cochain the refuted-DEC frame belongs to), [[issue_208_operator_algebra_frames]]
(the biproduct block structure = the block-diag-over-ordinate + boundary lift).
