# A3 — `solve_transpose` + wire `.H` to carry `CAP_SOLVE` (adjoint SN transport, #276)

> Implementation plan (companion to the verification spec `.claude/plans/a3_solve_transpose_verification.md`).
> Branch `feature/sn-adjoint-transport`. Surgical, **main-agent-direct, user-steered, NO method-implementer**;
> per-gate cadence. Campaign home: `.claude/plans/kiso_cross_model_extraction.md`. Phase spec:
> `.claude/plans/adjoint_sn_transport.md:130`.

## Goal

The SN loss `(L+C)` gains a **direct reverse-DAG forward-substitution `solve_transpose`** (solving
`(L+C)ᵀ x = b` exactly, single-pass — NOT Krylov-on-transpose, NOT μ-reversal), and
`_AdjointOperator.solve` wires `A.H.solve(b) = G_W⁻¹ · inner.solve_transpose(G_V · b)`, closing the
`operator.py:648-651` deferral so `(L+C).H` advertises `CAP_SOLVE`. This is the one genuine new
primitive of the adjoint chain (per-leaf table `adjoint_sn_transport.md:89`); A4 (daggered eigenvalue,
φ*) then feeds the UNCHANGED `power_iteration` the daggered triple whose inner solve is this.

## The structure — the fourth face of a 2×2 (the elegant framing)

|                       | **solve** (invert, triangular) | **apply** (matvec)             |
|-----------------------|--------------------------------|--------------------------------|
| **forward** (fwd DAG) | `sweep`                        | `loss_action`                  |
| **transpose** (rev DAG)| **`sweep_transpose`** ← NEW   | `loss_action_transpose` ← exists |

`sweep_transpose : sweep :: loss_action_transpose : loss_action`. The reverse-walk skeleton is ALREADY
built + validated in `_OneDimScanWalk.loss_action_transpose` (`loss_representation/__init__.py:2602`,
pinned by `test_g_adjoint_reciprocity`): reversed cell traversal, the boundary in↔out swap
(reads OUTFLOW cotangents / writes INFLOW), the `mirror = quad.reflection_index("x")` Carlson coupled-pole
seed permutation, `closure.angular_adjoint`, and the **ψ-independent coefficient reuse** (the SAME
`cell_balance_for_streaming` denom + `closure.cell_contribution` the forward uses — Pattern 2, no twin
algebra). `sweep_transpose` reuses that exact skeleton but with **SOLVE cell semantics** (solve for
`x_cell` directly in reverse order) instead of matvec cotangent-accumulation.

**Math (verified):** `A : V→W`, the G-adjoint `A* = G_V⁻¹ Aᵀ G_W`, so `(A*)⁻¹ = G_W⁻¹ (Aᵀ)⁻¹ G_V`.
Hence `A.H.solve(b) = G_W⁻¹ · solve_transpose(G_V · b)` — the metric usage **swaps** relative to
`_AdjointOperator.apply` (which is `G_V⁻¹ · apply_transpose(G_W · ·)`) because solve inverts. `solve_transpose(c)`
solves `Aᵀ x = c` (no metric — the bare Euclidean transpose inverse; the metric lives in the wrapper).

**Direct, not Krylov (test-architect):** the forward `sweep` is one forward-substitution pass (`sweep`→`_run`,
no outer loop), so `sweep_transpose` is exact single-pass. Slab tolerance is FP-only (`rtol≈1e-11`), NOT
`inner_tol`-scaled. On sphere/cyl the Carlson coupled-pole seed reads the iterate (`__init__.py:1005`) →
direct machinery iterate-threaded for the seed (`rtol≈1e-9`), still not Krylov. The forward matvec
`loss_action(sigma, psi)` has NO `initial_guess`, so the dense `M` (G2 oracle) is exact even on the sphere.

## Scope (follows `loss_action_transpose`'s scope — not new debt)

- **1-D `_OneDimScanWalk`** (cartesian slab + curvilinear sphere/cyl): IMPLEMENT. This is what A3's gates
  (slab + sphere) and A4 (∞-medium + slab + sphere) need.
- **Multi-D Cartesian**: DEFER (`raise NotImplementedError`, mirroring `loss_action_transpose:2647`).
- **2-D wavefront `_DAGWavefront`**: defer unless a gate needs it (A3 gates are 1-D). Decide at A3.2.
- `ScanMarch` multi-D: defers (consistent with its `loss_action_transpose:1859` deferral).

## Phases (each: green gate + forward-safety → present → commit when the user asks)

- **A3.0 — verification scaffold.** Land the test-architect harness `tests/sn/operators/test_loss_transpose_solve.py`:
  the **G2 dense `(L+C)ᵀ⁻¹` oracle** (build `M` from the FORWARD `apply` on flat composite basis vectors via
  `to_flat`/`from_flat` `full_field.py:337/355`; `solve_transpose(b) == np.linalg.solve(M.T, b)` — shares NO
  code with the reverse-walk, the keystone) + G1 round-trip; config builders (heterogeneous ≥2G per-group-σ_t
  slab AND a **sphere** leg — MANDATORY, the μ-reversal mirror only lives under `curvature != cartesian`;
  vacuum BC; `_random_composite` non-zero-boundary input). Gates RED until A3.2 (or co-develop).
- **A3.1 — the transpose-solve cell relation.** The reverse-order cell solve: invert the transposed WDD cell
  recurrence using the SAME ψ-independent coefficients `loss_action_transpose` reuses. Forward cell:
  `denom·ψ_cell = |μ|A_total·ψ_in + angular_numer + V·q`, then `ψ_out = 2ψ_cell − ψ_in`. The transpose-solve
  inverts the transposed recurrence walking reverse DAG order. Mirror `loss_action_transpose`'s exact algebra
  (the reverse recurrences at `:2741-2762`) — but SOLVE for the cell unknown, do not accumulate cotangents.
- **A3.2 — `_OneDimScanWalk.sweep_transpose(sigma, b) → FullField`.** The reverse-DAG forward-substitution,
  reusing the `loss_action_transpose` walk frame (reverse cells, boundary swap, `mirror`, `angular_adjoint`)
  with the A3.1 solve cell. Add `sweep_transpose` to the `LossRepresentation` Protocol (`__init__.py:229`) +
  the 1-D rep; multi-D/2-D `raise NotImplementedError`. GATES **G1** (round-trip `solve_transpose ∘ apply_transpose = I`)
  + **G2** (dense oracle), slab + sphere; mutations: forward-DAG order, +μ/−μ mirror → RED.
- **A3.3 — `InvertibleOperator.solve_transpose` + `CAP_SOLVE_TRANSPOSE`.** Route `solve_transpose(b) →
  representation.sweep_transpose(self.sigma, b)` (mirror `apply_transpose → loss_action_transpose` at
  `streaming.py:757`); add `CAP_SOLVE_TRANSPOSE` to `InvertibleOperator.capabilities`.
- **A3.4 — `_AdjointOperator.solve` + wire `CAP_SOLVE`.** Add `CAP_SOLVE_TRANSPOSE = "solve_transpose"`
  constant (`operator.py`); implement `_AdjointOperator.solve(b) = inner_codomain.apply_inverse_metric(
  inner.solve_transpose(inner_domain.apply_metric(b)))` (the metric-swap vs `apply`); advertise `CAP_SOLVE`
  iff inner has `CAP_SOLVE_TRANSPOSE`; retire the 648-651 "solve does NOT propagate" deferral comment.
  GATES **G4** (Mode-11 sentinel: `(L+C).H.solve` executes `solve_transpose`, counter>0, NOT forward `solve`,
  counter==0), **G5** (`(L+C).H.apply((L+C).H.solve(b)) == b` value round-trip), **G6** (cap flip: net-new
  `CAP_SOLVE` on `(L+C).H` ONLY; `S.H`/`F.H`/`L.H` stay solve-less; named 0-ULP forward canaries).
- **A3.5 — G3 full-loss `(L+C−S−B)` reciprocity (the apply-side unblock).** Extend
  `test_g_adjoint_reciprocity.py` from `(L+C−B)` to the full `A_loss = L+C−S−B` now that S† exists
  (per-group via one-hot-per-group φ, vv L27; asymmetric SigS ≥2G; S lifts to FullField). **NOTE:** G3 is an
  adjoint-MATVEC gate (the full `L+C−S` is NOT sweep-invertible — S couples groups), not a `solve_transpose`
  gate; A3 *unblocks* it and it feeds A4's daggered eigenvalue posing. Keep per the brief's Deliverable 2.

## Capability-propagation decision

`CAP_SOLVE_TRANSPOSE` does NOT need to propagate through `OperatorSum`/`OperatorProduct`/`ScaledOperator`
for A3 — the loss solve lives on the fused `InvertibleOperator`, not `OperatorSum` (sum-solve never
propagates). `_AdjointOperator.solve` only needs the inner (the `InvertibleOperator`) to carry
`CAP_SOLVE_TRANSPOSE`. Revisit ONLY if A4's daggered posing composes operators that need it (the daggered
`A_loss` feeds `power_iteration` via the `InvertibleOperator`'s solve — expected fine).

## Critical files

- **NEW primitive:** `orpheus/sn/loss_representation/__init__.py` (`sweep_transpose` on the Protocol `:229` +
  `_OneDimScanWalk` near `loss_action_transpose` `:2602`; the cell solve A3.1).
- **cell relation:** `orpheus/sn/spatial/cell_balance.py` (`cell_balance_for_streaming` — reuse, don't twin)
  + `orpheus/sn/spatial/diamond.py` (the closure `cell_contribution`/`angular_adjoint`).
- **operator wiring:** `orpheus/sn/operators/streaming.py` (`InvertibleOperator.solve_transpose` + cap, near
  `apply_transpose` `:757`); `orpheus/numerics/operator.py` (`CAP_SOLVE_TRANSPOSE` const, `_AdjointOperator.solve`
  `:642-702`, the deferral `:648-651`).
- **gates:** `tests/sn/operators/test_loss_transpose_solve.py` (NEW, G1/G2/G4/G5/G6),
  `tests/sn/operators/test_g_adjoint_reciprocity.py` (extend, G3).

## Discipline

Surgical, main-agent-direct, user-steered, NO method-implementer. Per-phase: green canonical gate
(`.venv/bin/python -O -m pytest <paths> -m "not slow" -q -rfE --timeout=300 -p no:xdist -p no:cacheprovider`)
+ CLI-pyright `orpheus/` ≤ **412** (should net DOWN — retiring the solve-deferral confession). Commit/push/merge
ONLY when the user asks; stage explicitly; trailer `Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>`;
no `# type: ignore` (`cast` OK); NEVER `git checkout/restore/stash` on uncommitted files. **Re-baseline, NOT
bit-identity** for the new transpose-solve values; the G2 dense oracle + reciprocity are the structural
completeness proofs. The G2 dense oracle (independent of the reverse-walk code) is the keystone safety net
against a bug copied into both `apply_transpose` and `sweep_transpose`.
