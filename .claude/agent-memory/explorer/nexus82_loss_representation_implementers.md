---
name: nexus82-loss-representation-implementers
description: Ground truth for `.. implements::` on the 17 equations of theory/methods/sn/loss_representation — 12 have implementers, 4 are laws with NONE, 1 is doc-drifted; plus the finding that the true LpC implementers were absent from the guess pool entirely.
metadata:
  type: reference
---

# nexus#82 — who implements the `loss_representation` page's 17 equations

Determined 2026-08-17 against `main` @ `a1c90aac` (graph build `8bb695d3`), by reading
the equation TEXT first and then the code — never the existing `implements` edges (all
466 are token guesses on `loss` / `inverse` / `adjoint`).

⚠ **`file:line` below is re-derive-via-Nexus, not the headline.** The durable claims are
§2 (which equations are LAWS, not computations) and §3 (the Resolution-A mechanism drift).

## 1. The verdict table

| label | implementer (dotted path) | file:line | why |
|---|---|---|---|
| `loss-rep-LpC` | `orpheus.sn.operators.streaming.StreamingCollisionOperator` | `orpheus/sn/operators/streaming.py:445` | the object that IS `(L+C)` — an `OperatorSum(StreamingOperator, MultiplicationOperator)` whose `.apply` is the LHS and `.solve` the system's solution |
| | `orpheus.sn.operators.streaming.StreamingOperator` | `orpheus/sn/operators/streaming.py:169` | the discretised `L = Ω·∇\|_WDD` half (streaming + curvilinear angular redistribution, σ-free since #257 S8b) |
| `loss-rep-resolution-a` | `orpheus.sn.loss_representation._LossRepresentation.streaming_action` | `orpheus/sn/loss_representation/__init__.py:485` | ⚠ see §3. The one symbol whose value IS `Lψ` derived from the full loss. Realised as `loss_action(σ=0)`, NOT as a subtraction |
| | `orpheus.sn.operators.streaming.StreamingOperator.apply` | `orpheus/sn/operators/streaming.py:312` | the operator-level door returning `Lψ`; the `-C` glue is its documented contract and is pinned numerically by `test_loss_action_convention.py` against an independent `M[σ_t]` |
| `loss-rep-affine-cell` | `orpheus.transport.spatial.diamond.DiamondDifference.residual_kernel_batch` | `orpheus/transport/spatial/diamond.py:444` | literally `residual = denom * psi_bar - numer` with `denom = Σ_t + Σ_a 2g_a` |
| | *(secondary, curvilinear arm)* `orpheus.sn.loss_representation._OneDimScanWalk._apply_walk` | `orpheus/sn/loss_representation/__init__.py:3039` (formula at `:3264`) | `m_full = (denom*psi_cell - numer_upstream)/V` — same form, curvilinear population, via `orpheus.transport.spatial.cell_balance.cell_balance_for_streaming` (`cell_balance.py:123`), which builds `denom`/`numer_upstream`. Declare only if the equation is meant to cover curvilinear; the written `c_a = 2\|μ_a\|/Δa` is the Cartesian spelling |
| `loss-rep-affine` | `orpheus.sn.loss_representation._LossRepresentation.streaming_action` | `orpheus/sn/loss_representation/__init__.py:485` | the equation NAMES `streaming_action`; the method exists *because* of this identity (`return self.loss_action(self._zero_sigma_for(psi), psi)`) |
| | `orpheus.sn.loss_representation._LossRepresentation.streaming_action_transpose` | `orpheus/sn/loss_representation/__init__.py:492` | the transpose sibling, same identity at σ=0 |
| `loss-rep-leaf-sum` | **NONE** | — | derivation identity about a superseded path; see §2 |
| `loss-rep-removal-sigma` | **NONE** | — | notation/definition with no production site; see §2 |
| `loss-rep-affine-kernel-maps` | `orpheus.sn.loss_representation.assembly._probe_coefficient_blocks` | `orpheus/sn/loss_representation/assembly.py:187` | returns `(diag, inflow, trace, memory)` = `(Ā, C_b, T_a, α_ab)` and its unit-probe reading is valid *only* because the maps are affine — it is the code the equation licenses |
| | *(the kernels the affine form describes)* `DiamondDifference.residual_kernel_batch` / `orpheus.transport.spatial.linear_discontinuous.LinearDiscontinuous.residual_kernel_batch` | `diamond.py:444` / `linear_discontinuous.py:628` | DD gives the 1×1 blocks, LD the `2^d×2^d` UBLD blocks; declare if the equation is read as a statement about the kernels rather than about the extraction |
| `loss-rep-walk-order-rows` | `orpheus.sn.loss_representation.assembly.assemble_ordinate_blocks` | `orpheus/sn/loss_representation/assembly.py:276` | the walk body IS the equation: `cell_rows -= inflow[b] @ row_in`; `row_out = trace[a]@E_c + Σ_b memory[a][b]@row_in`. (`ordinate_walk_order`, `assembly.py:162`, produces the permutation the triangularity claim is stated in — a companion, not the rows) |
| `loss-rep-sweep-global-conjugation` | `orpheus.sn.loss_representation.assembly.assemble_ordinate_blocks` | `orpheus/sn/loss_representation/assembly.py:276` (nested `_to_global_frame` at `:371`) | `vals * dof_signs[rows] * dof_signs[cols]` = `Φ M Φ`. The nested closure is not a graph node, so the enclosing function is the addressable implementer. Φ comes from `octant_moment_frame_signs` (that is `ld-ubld-octant-moment-frame-signs`'s equation, not this one) |
| `loss-rep-scanmarch` | `orpheus.sn.loss_representation.ScanMarch` | `orpheus/sn/loss_representation/__init__.py:2172` | the class IS the schedule; the piecewise d=1/d=2 split lives in `.sweep` (`:2253`) and `.loss_action` (`:2408`) |
| | `orpheus.sn.loss_representation.ScanMarch._sweep_interior` | `orpheus/sn/loss_representation/__init__.py:2302` | the `scan(x)∘march(y)` body, SOLVE direction |
| | `orpheus.sn.loss_representation.ScanMarch._loss_action_interior` | `orpheus/sn/loss_representation/__init__.py:2444` | the `scan(x)∘march(y)` body, APPLY direction |
| `loss-rep-scanmarch-solve` | `orpheus.transport.spatial.diamond.DiamondDifference._cartesian_streaming_diagonal` | `orpheus/transport/spatial/diamond.py:327` | the single source of `couplings[a] = 2 g_a` and the explicit left fold `S = ((Σ_t + 2g_0) + 2g_1) + …`; consumed by all three Cartesian producers |
| `loss-rep-scanmarch-solve-affine` | `orpheus.transport.spatial.diamond.DiamondDifference.cartesian_scan_coefficients` | `orpheus/transport/spatial/diamond.py:694` | `a = 2.0*scan_diag*inverse_denom - 1.0`, `inverse_denom = 1/denom`, `w = _DD_W` |
| | `orpheus.transport.spatial.scheme.DiscretizationSchemeBase.source_emission` | `orpheus/transport/spatial/scheme.py:1382` | `return QV * inverse_denom / w` — the β half |
| | `orpheus.sn.loss_representation.ScanMarch._sweep_interior` | `orpheus/sn/loss_representation/__init__.py:2302` | forms `QV_eff = Q + c_y·ψ_{y,in}` (`:2394`) — the transverse fold inside β's numerator lives here, nowhere else |
| `loss-rep-scanmarch-apply` | `orpheus.transport.spatial.diamond.DiamondDifference._reflection_coeffs` | `orpheus/transport/spatial/diamond.py:368` | `alpha = -(1.0-w)/w`, `beta = psi_bar/w` — the w-generic arithmetic |
| | `orpheus.transport.spatial.diamond.DiamondDifference.reflect_scan_coefficients` | `orpheus/transport/spatial/diamond.py:750` | DD's instantiation at `w = _DD_W` giving the stated `(−1, 2ψ̄)`. ⚠ `_x_scan_faces` (`sn/sweep/scan.py:302`) only RUNS the recurrence — a consumer, not an implementer |
| `loss-rep-scanmarch-apply-residual` | `orpheus.transport.spatial.diamond.DiamondDifference.residual_kernel_batch` | `orpheus/transport/spatial/diamond.py:444` | the ÷V matvec kernel; `ScanMarch._loss_action_interior` calls it (`:2497`) — the caller is a consumer |
| `loss-rep-facewise-separable` | **NONE** | — | definitional iff-criterion; see §2 |
| `loss-rep-adjoint-inverse-swap` | `orpheus.numerics.operator._AdjointOperator.inverse` | `orpheus/numerics/operator.py:1335` | `return inner_inverse.H` — the RHS spelled literally; this is what makes the law an object identity rather than a numerical coincidence |
| `loss-rep-metric-adjoint-solve` | `orpheus.numerics.operator._AdjointOperator.apply` | `orpheus/numerics/operator.py:1288` | `G_V⁺ ⊙ apply_transpose(G_W ⊙ y)` — the conjugation, verbatim. There is deliberately no `_AdjointOperator.solve` |
| | `orpheus.sn.operators.sweep_operator.SweepOperator.apply_transpose` | `orpheus/sn/operators/sweep_operator.py:158` | supplies the second equality `(A⁻¹)ᵀ = solve_transpose`, by returning `self.inner.solve_transpose(b)` |

## 2. The NONE set — what the inference should never have guessed at

Four equations have **no** implementer. Three are LAWS or DEFINITIONS; one is a
derivation step about a path the tree deliberately overrides.

- **`loss-rep-leaf-sum`** (`L.apply(ψ) + C.apply(ψ) = M(σ_r)ψ`) — a *derivation identity*
  proving the pre-carve leaf sum was value-correct by coincidence. Doubly historical: (a)
  `StreamingCollisionOperator.apply` **overrides** `OperatorSum.apply`, so the leaf sum is
  unreachable for `(L+C)` in production; (b) since #257 S8b `L.apply` returns
  `streaming_action(ψ)` and never sources `σ_t` at all, so the equation's own underbrace
  (`M(σ_t)ψ − σ_tψ`) describes a path that no longer exists. `OperatorSum.apply` computes
  the generic LHS but is not accountable to *this* claim.
- **`loss-rep-removal-sigma`** (`σ_r = σ_t − Σ_{s,0}^{g→g}`) — a notation line. The page
  says it itself: *"There is no production caller of the removal form yet."* `[M]`
  `MaterialXSField.foldable_sigma` (`material_xs_field.py:1213`) returns only the
  **subtrahend** `Σ_{s,0}^{g→g}`, never `σ_r`. The only sites forming `σ_t − σ_s0` are in
  `orpheus/derivations/discrete/sn/dsa.py` (`:632`, `:1023`) and that is the **DSA
  low-order removal** `σ̂_R` — a different operator in a different physics context.
  Declaring it would be a false attribution.
- **`loss-rep-facewise-separable`** (`M_d = ⊕_a M^(1)_a ⟺ …`) — a definitional criterion
  naming a property. The property is carried as a **declared `ClassVar[bool]`**
  (`DiscretizationSchemeBase.transverse_coupling_is_facewise`: `scheme.py:793` default
  `False`, `diamond.py:147` `True`), never computed. `ScanMarch.supports` READS the tag;
  `ScanMarch._sweep_interior` is exact *because of* the criterion but does not compute it.
- **(borderline, resolved to a declaration)** `loss-rep-adjoint-inverse-swap` reads like a
  law but is NOT NONE: `_AdjointOperator.inverse` is written to make it true by
  construction, so the code genuinely IS the identity.

## 3. ⚠ Doc drift found — `loss-rep-resolution-a`'s stated mechanism is stale

The page's §"Resolution A — the operator's only glue" says `StreamingOperator.apply`
*subtracts* `σ_t ⊙ ψ.bulk` from the full loss. `[M]` at HEAD it does not:

- `StreamingOperator.apply` (`streaming.py:349`) → `self.loss_representation.streaming_action(psi)`
- `_LossRepresentation.streaming_action` (`__init__.py:490`) → `self.loss_action(self._zero_sigma_for(psi), psi)`

i.e. `L` is realised as `(L+C)|_{σ=0}`, not as `(L+C) − C`. The two are equal **by
`loss-rep-affine`**, so the equation is still true and still gated — but as the
*justifying identity*, not as transcribed arithmetic. `grep` for the subtraction returns
nothing in `orpheus/`. The change is #257 S8b; the gate that pins the new mechanism is
`tests/sn/operators/test_streaming_operator_decomposition.py::TestPureLIsLossActionAtZeroSigma`,
whose docstring records that this class once pinned *the opposite*.

Stale prose that inherits it (all present-tense-false): `loss_representation/__init__.py`
docstrings at `:333`, `:653`, `:1096`, `:1215`, `:3017`, `:3431`, and
`streaming.py:547-548` (*"Resolution A subtractive form: L.apply(ψ) = M(ψ; σ_t) − σ_t ⊙ ψ_cell"*).

Second observation: the page **overloads** `loss-rep-resolution-a` with a second claim —
"the composite owns its matvec, single-sourcing σ" (`StreamingCollisionOperator.apply`,
`streaming.py:659`), which is what `test_removal_form_matvec_sweep.py` actually verifies
under that label. The equation TEXT states only `L = (L+C) − C`, so the table above
declares to the text. If the label is meant to carry both claims,
`StreamingCollisionOperator.apply` / `.apply_transpose` (`:659` / `:680`) belong too.

## 4. The guess pool was not merely noisy — for one equation it was disjoint from the truth

`loss-rep-LpC`'s 23 guesses are all token matches on `loss` **inside
`orpheus.sn.loss_representation`** (plus the two `.loss_representation` property
accessors). Neither `StreamingCollisionOperator` nor `StreamingOperator` — the actual
`(L+C)` and `L` — appears in the pool at all. For `loss-rep-adjoint-inverse-swap` the
opposite: the one true implementer (`_AdjointOperator.inverse`) IS in the pool, ranked
31st of 35 by degree. So a "take the top guess" heuristic is wrong in both directions.

Related: [[lessons-L1]] (blast radius is graph + grep + constructors + doc nodes),
[[affine-operator-split-convention]].
