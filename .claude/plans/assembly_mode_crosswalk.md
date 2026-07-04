# Assembly-mode crosswalk — 2b design contract (stencil-assembly campaign)

**Status: IN EXECUTION (2026-07-04, branch `refactor/spatial-promotion-assembly`).**
The binding inputs: roadmap 2b (`stencil_assembly_dsa_roadmap.md`) + R2 + the two
2-P0 crux blocks + test-architect L16. This file records the *realized* design —
the convention contract every emitter and gate cites. Verified against the tree
@ `f4a6749`.

## The one law

`assemble()` is the THIRD consumption mode of the ONE closure algebra
(solve = sweep substitution, apply = matvec, assemble = emit the SAME
coefficients as (row, col, value) into a scipy-sparse carrier). It is the same
additive-monoidal functor `Op → Mat` that `as_matrix` already documents
(`operator.py` "the functor OUT of the operator category"), landing in a sparse
carrier: leaves override `assemble()`, composites recurse via the homomorphism
laws (Sum→`+`, Product→`@`, Scaled→`*`; TensorProduct→kron DEFERRED — no 2b
consumer, the diffusion loss and SN L(+C) trees contain none).

## The numerics home (commit 1)

- **`orpheus/numerics/assembled_operator.py`** (new leaf module — the
  `flat_operator.py` / `matrix_inverse_operator.py` placement family):
  `SparseAssembledOperator(LinearOperator)` — thin wrapper over a
  `scipy.sparse` CSR matrix (constructor accepts COO; `.tocsr()` performs the
  FEM duplicate-summing). `apply(x)` = CSR matvec on **flat 1-D vectors** (the
  `FlattenedOperator` serialization convention — an assembled operator is a
  serialization, not new structure); `apply_transpose` = `matrix.T @ x`
  (`is_adjointable` True); `as_matrix` override = densify (honoring the
  basis-shape-consistency + `max_dimension` contract, the
  `MatrixInverseOperator.as_matrix` precedent); `is_assemblable` True and
  `assemble()` returns `self` (idempotent — the functor is closed).
- **`operator.py` three-layer surface** (the SupportsInverse idiom, exactly):
  1. predicate `is_assemblable` (base default `False`; composites recurse:
     Sum/Product = both legs, Scaled = operand);
  2. narrowing target `SupportsAssembly(LinearOperator, Protocol)` declaring
     `assemble() -> SparseAssembledOperator` (TYPE_CHECKING import — no
     runtime cycle);
  3. checked bridge `assemblable()` (`TypeGuard`).
  Refusal exception: **`MissingAssembly(TypeError)`** — the assembly-axis
  sibling of `NotInvertible` / `MissingAdjoint`, raised by the composer
  `assemble()` bodies when an operand cannot emit.
- **R2 realized in `as_matrix`**: the probing loop is EXTRACTED byte-identical
  to `_as_matrix_by_probing(shape)` (the RETAINED fuller-view pathway G3
  forces); `as_matrix` keeps its gate resolution, then delegates to
  `assemble().matrix.toarray()` when `assemblable(self)` (with a
  dimension-consistency ValueError), else probes. No operator is assemblable
  until an emitter lands ⟹ the delegation is a no-op for every existing call
  site at commit 1 (bit-safety by construction).
- **`FlattenedOperator`**: `is_assemblable` = inner's; `assemble()` = inner's
  (flat-size-checked). The flat serialization is transparent to assembly —
  the typed operator's emission is ALREADY in the flat layout.

## Global DOF numbering (never re-derive)

The flat functor is `FullField.to_flat()` = `concat(bulk.values.ravel(),
boundary.values)` — verified at `full_field.py:367-382`:

- **bulk block** = C-ravel of `bulk.values`; diffusion scalar bulk `(ng, nx)`
  ⟹ flat index `g·nx + i` (group-major).
- **trace block** at offset `n_bulk`; the scalar trace's per-face slot is
  `(2, ng)` C-ordered at `FaceLayout.faces[f].offset` with row 0 = J⁺
  (`ScalarTraceSpace.OUTFLOW_ROW`), row 1 = J⁻:
  `col(J±, f, g) = n_bulk + offset(f) + row·ng + g`.
- `n_bulk` / `n_trace` read off the operator's own `FullFieldSpace`
  (`bulk_space.shape` / `trace_space.shape`) — the `OperatorSum` equal-domain
  guard is what makes summing assembled matrices sound (all leaves of one loss
  carry the SAME composite space by construction).

## Diffusion emitters (commit 2 — the FIRST consumer)

One-source table (the coefficient source each `assemble()` consumes — the
SAME attributes/kernels `apply` consumes; L16's mutation teeth bite here):

| Leaf | Coefficient source | Emission |
|---|---|---|
| `LeakageOperator` (L) | `self._conductance` / `_areas` / `_volumes` / `_face_closures` (precomputed from `_interior_conductance` / `_boundary_closure`) | direct structural: interior-face ±A·g/V scatter (per-face contributions duplicate-summed on the bulk diagonal — nulp vs apply's diff-then-divide grouping, per L16), edge-row trace coupling ±A/V on (J⁺,J⁻), outflow-defect row (1, −c_φ, −c_J⁻), inflow identity row |
| `DiffusionBoundaryOperator` (B) | `mesh.bc` realized laws | per-face `law.as_matrix(basis_shape=(ng,))` onto the inflow row — ng probes THROUGH the law's own apply (one-source; the laws are Zero/Identity/Scaled·I) |
| `MultiplicationOperator` (C) | `self.coefficient.values` (the same array the engine multiplies) | bulk diagonal `broadcast_to(values, bulk_shape).ravel()` — family-blind (prepend-broadcast covers both the scalar `(ng,nx)` and a future angular `(N,ng,…)` bulk); bit-exact vs apply (one multiply per entry) |
| `IsotropicScattering` / `IsotropicN2N` (S) | the production einsum kernels `mat_xs.apply_p0_in_scatter` / `apply_n2n` | **group-impulse probing**: ng bare-kernel calls on `e_g' ⊗ 1_cells` extract every cell's `ng×ng` block exactly (unit inputs ⟹ exact coefficients), scattered block-diagonal over cells. Deliberately NOT `dense_per_material()` — that is the storage-side transpose-convention ORACLE (vv L11) and must stay realization-independent of production |

`is_assemblable` on C/S leaves: True iff `space` is a block-bearing
`FullFieldSpace` (space-anonymous instances honestly refuse — no layout to
emit into). L/B always (they carry the mesh).

Consumer switch: `MatrixInverseOperator(FlattenedOperator(loss, template))`
picks up the assembled matrix AUTOMATICALLY through the as_matrix delegation
once every loss leaf assembles — the existing diffusion suite becomes the
regression net (all its k/trace/balance gates run the assembled path; nulp
matrix movement ⟹ k moves ≪ every gate tolerance). Sparse-LU resolvent =
optional perf follow-up, NOT taken in 2b.

## SN emitters (commits 3–4) — Cartesian ONLY

Per-ordinate-per-group spatial block M over the bulk C-ravel, by a SYMBOLIC
walk of the sweep DAG in walk order (the same `SweepDependencyGraph` order the
sweep and matvec use): carry per-face the affine row (sparse coefficient row
over upstream cell DOFs ⊕ boundary-inflow constant); per cell emit
`denom·e_cell − Σ_a coupling_a·face_row_a`, then update
`face_out_row_a = (1/w)·e_cell + α·face_row_a` (α = −(1−w)/w).

- DD: denom + couplings from `_cartesian_streaming_diagonal` (the ONE
  Cartesian source `cell_kernel_batch` / `residual_kernel_batch` /
  `cartesian_scan_coefficients` all consume); w = ½ ⟹ α = −1, 1/w = 2 — the
  face chain has MEMORY: the assembled block is triangular-in-walk-order with
  dense chains along each axis (1-D: dense lower triangular — the honest
  matrix of the marching reconstruction; 2-D: per-axis chains, O(nx+ny)/row).
- LD: cell blocks from `assemble_ubld` (the SAME batched systems
  `per_cell_solve`/the LD apply einsum consume), inflow scatter from
  `assemble_inflow_axis`, upwind face trace from the cell moments — block
  scatter; the #253 `cell_moment_count`/`face_moment_count` policy is the DOF
  stride source.
- Curvilinear: OUT (the M-M angular chain couples ordinates; #282's lagged
  pole seed is the walk-order back edge). The #282 characterization gate
  asserts the defect POSITIVELY (L16: `np.any(triu≠0)` loud-flip, never
  xfail).

## Gates (L16 verbatim, the teeth)

- G1 `assembled@x ≡ apply(x)` — nulp/rtol≈1e-11, never array_equal, never a
  scalar functional; het + asymmetric SigS + non-uniform h; non-flat
  fixed-seed random x; ≥2G.
- G2 `scipy solve_triangular(PᵀMP) ≡ sweep solve` — rtol; LAPACK independence
  EARNS the #284 discharge; triangularity leg `triu(PᵀMP, 1) == 0` exact.
- G3 `_as_matrix_by_probing ≡ assemble().to_dense()` — exactly ONE pin per
  family (diffusion loss + one DD slab block); the delegation tautology trap.
- Teeth: monkeypatch sign-flip in the SHARED source (`_interior_conductance` /
  `_cartesian_streaming_diagonal` / `assemble_ubld`) must red BOTH the new
  gates AND the existing sweep/stencil suites; only-new-reds ⟹ twin ⟹ STOP +
  ERR-NNN. Negative test: every leaf `apply` monkeypatched to raise ⟹
  `assemble()` still succeeds.
- New tests under `tests/transport/spatial/` (assembly core + SN emitters) and
  `tests/diffusion/` (the diffusion equivalence gates, beside the stencil
  gate they extend).
