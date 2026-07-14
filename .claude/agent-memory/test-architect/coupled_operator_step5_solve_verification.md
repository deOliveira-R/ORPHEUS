---
name: coupled-operator-step5-solve-verification
description: The Phase-C step-5 (#41) "THE BLOCK SOLVE" gate spec for the Coupled Block Operator campaign (branch refactor/sn-walk-unification @ d881243d, step 4 + 4e COMPLETE, wall 6340/0 · 37 xfailed). Step 5 poses the resolvent block-triangular DIRECT (M re-poses as the honest upper-triangular CoupledOperator grid [[LC,Seeding],[None,A_BB]]; M.solve = block back-substitution) + the scattering ITERATIVE (SI over M⁻¹N unchanged) + the ρ-honest block residual r=Aψ−q via evaluate_residual as the STOP test, not ‖Δψ‖. THE headline architectural finding: grid.is_invertible refers to the LOSS grid A=system.loss (the full 2×2 CoupledOperator), invertible via the LU/materialize branch (MatrixInverseOperator) — NOT M; step 5 flips CoupledOperator.is_invertible → the dormant xfail test_a2a_grid_swap_law_inverse_arm (coherence :659) goes LIVE, requiring the NET-NEW MatrixInverseOperator.apply_transpose (LUᵀ free). Load-bearing gates: row-6 dense-LU oracle (M.solve ≡ np.linalg.solve(dense(M)), ~5.5e-16, psi_half :391 — the durable floor the substitution re-points onto), the ρ-honest stop unit gate (r_running ≡ N·Δψ AND ≡ evaluate_residual when M exact), the in-M-lag CERTIFICATE (the free-identity is blind to in-M lag by construction — #282's seed-lag was IN M — the end-of-solve evaluate_residual closes exactly that hole), φ=Q/Σ_t + k_eff anchors, singular-grid loud-fail (MatrixInverseOperator LinAlgError at construction). THE stop-test re-baseline ledger: frozen k-snapshots (test_dd_regression :100) move ~inner_tol-scale (solver.py:835-837 — a stopping-norm change shifts converged k, keff_tol-tight snapshots can't absorb); the ρ=c convergence-rate ratios (test_si_convergence_rate) are RELATIVE→robust. D2a/D3/D5 user-ruling branches spec'd. Sibling of step4 / 4e-unweave memos (L7/L13/L16/L18/L19).
metadata:
  type: reference
---

# Verifying step 5 (#41) — THE BLOCK SOLVE (resolvent direct + scattering iterative + ρ-honest stop)

**When this recipe applies.** A coupled-block campaign has landed the N-general
`CoupledField`/`CoupledOperator` machinery (step 4), re-typed the off-diagonal
blocks + built `build_within_group_system` (the co-producing builder), evicted the
ray from FullField, and UN-WOVEN the walk (4e: the walk routes System B through the
extracted `RadialCharacteristicOperator`/`…Seeding` engines). What REMAINS: the
resolvent `M` is still spelled as the `CoupledInvertibleOperator` fused-walk facade
(its `solve` calls `self.fused.solve(rhs_a, radial_characteristic_source=rhs_b,
radial_characteristic_flux=flux_buf)` with explicit leaf kwargs — coupled_system.py
:383-400), and the LOSS grid `A` (the full 2×2 `CoupledOperator`) has NO solve
(`is_invertible` base `False`, coupled_system.py:580-582). **Step 5 is the solve:**
(i) re-pose `M` as the HONEST upper-triangular grid `[[LC, Seeding], [None, A_BB]]`
whose `solve` is block back-substitution (the fused delegation dissolves); (ii) give
the numerics `CoupledOperator` a `solve`/`is_invertible`/`inverse` mode; (iii) flip
`grid.is_invertible` on the LOSS grid via the materialize/LU branch, converting the
dormant swap-law xfail to a live gate; (iv) wire the ρ-honest block residual
`r = Aψ − q` (`evaluate_residual`, solver.py:240) as the STOP test. Governing memos
(same dir): `coupled_operator_step4_verification.md` (the M1-M5 machinery + the
assemble≡probe floor + A2 reciprocity), `coupled_operator_4e_unweave_verification.md`
(the survivors ledger + row-6 oracle + the L18 cold-residual acceptance number).

**Verified against HEAD `d881243d`.** Every `file:line` below is read off that tree
(⚠ best-estimate for production internals — grep the symbol, not the line). Canonical
`.venv/bin/python -O -m pytest -p no:xdist --timeout=300 -p no:cacheprovider` SERIAL;
every tooth in-process `monkeypatch` (NEVER `git checkout` an uncommitted file — L4);
every gate fires under `-O` (`np.testing`/`pytest.fail`, never bare `assert` outside a
pytest-rewritten test body — vv Mode 8).

## THE headline architectural finding (resolve BEFORE writing gates)

**`grid.is_invertible` in the dormant xfail refers to the LOSS grid `A = system.loss`
(the full 2×2 `CoupledOperator`), NOT the resolvent `M`.** Verified: the xfail
`test_a2a_grid_swap_law_inverse_arm` (coherence :659-682) reads `grid, _space =
_within_group_grid(sn)` where `_within_group_grid` (:593) returns `system.loss` — the
`CoupledOperator` grid `[[A_AA, A_AB], [-emission, A_BB]]`. So there are TWO distinct
solve objects, and step 5 realizes both:

| object | shape | invertibility route | solve mechanism | consumer |
|---|---|---|---|---|
| **M** (`WithinGroupSystem.resolvent`) | upper-triangular `[[LC, Seeding], [None, A_BB]]` | **triangular** (None-pattern strictly-upper ∧ every diagonal block invertible) | block back-substitution: `ψ_B=A_BB.solve(q_B)`; `ψ_A=LC.solve(q_A−Seeding·ψ_B)` | the per-iteration SI step (`M.inverse()` feeds `SourceIteration`) |
| **A** (`WithinGroupSystem.loss`) | full 2×2 `[[A_AA, A_AB], [-emission, A_BB]]` | **materialize/LU** (square ∧ below size gate → `MatrixInverseOperator(self)`) | one LU of the assembled/probed matrix + backsolve | the EXTRACT dense-LU oracle + the swap-law inverse arm (flips the xfail) |

`M` is upper-triangular because the emission (B←A) is LAGGED into `N` (the gains) —
`M`'s `(B,A)` slot is structurally `None`, so System B (the ray march) does not depend
on System A, and back-substitution (System B first → feeds the seed to System A) is
valid AND is exactly the curvilinear sweep order. `A` is a FULL 2×2 (both off-diagonals
present) — NOT triangular — so its inverse is the whole-system direct solve via
materialization ("if we assembled the matrix it must work"). Do NOT conflate them: the
triangular-substitution `M.solve` is the row-6-oracle-pinned per-iteration step; the
materialize-LU `A.inverse()` is the swap-law/EXTRACT-oracle whole-system solve.

## The FIVE landings — the equivalence bar SPLITS by sub-commit

| sub-commit | what lands | equivalence bar | oracle | RED-if |
|---|---|---|---|---|
| **5a (D1) numerics solve** | `CoupledOperator.solve`/`solve_transpose`/`is_invertible`/`inverse` (triangular detect + materialize branch) + `MatrixInverseOperator.apply_transpose` (LUᵀ) | **structural (constructability) + LUᵀ value** | synthetic upper/lower-triangular toy grids; `M.T` dense cross-check | non-triangular claims triangular / substitution order wrong / `LUᵀ ≠ transpose-then-LU` |
| **5b (D2) M re-pose** | `resolvent` = the honest triangular grid; `M.solve` = substitution; `M.solve_transpose` = forward-substitution on Mᵀ | **principled-equiv `< 1e-11` (~5.5e-16), NOT bit-id** (M-M row) | **row-6 dense-LU** `test_extract_to_dense…` (psi_half :391) re-points onto the substitution | value off the 5.5e-16 floor / substitution drops the Seeding subtraction |
| **5c (D5) ρ-honest stop** | free-identity `r=N·Δψ` in `SourceIteration` (all arms) + ONE end-of-solve `evaluate_residual` certificate (solver) | **unit-identity (r_running ≡ N·Δψ) + certificate ≡ evaluate_residual** | explicit `N.apply(Δψ)` on a small config; `evaluate_residual(A,ψ,q)` when M exact | free-identity normalized by ‖ψ‖ not ‖q‖ / a gain dropped from rhs bookkeeping / in-M lag NOT caught by the certificate |
| **5d (D2a/D3) retirement** | `CoupledInvertibleOperator`/`CoupledSweepOperator` retire (or stay facades); the :712 fused-matvec twin (route-through-grid or keep) | **branch-dependent (bit-id keep / principled-equiv route)** | `test_si_single_primitive_contract` isinstance rows re-point | a retired symbol's tooth goes silently toothless (grep the deletions) |
| **5e xfail→live + anchors** | `grid.is_invertible` flips True; `test_a2a_grid_swap_law_inverse_arm` xfail→live; `Q/Σ_t` + k_eff anchors; the stop-test re-baseline ledger | **swap-law `array_equal` (object identity) + closed-form 1e-10** | the `MatrixInverseOperator(grid).H` graph both spellings run; φ=Q/Σ_t + k_inf | swap-law not bit-id / anchor drift beyond re-baseline band |

**The DURABLE load-bearing gate is 5b's row-6 oracle** (`M.solve ≡ np.linalg.solve(
dense(M), b)` at ~5.5e-16) — it was DESIGNED (step 0) as the anchor every EXTRACT/re-pose
step pins against, and the M-substitution re-pose IS such a step. **The highest-risk NEW
teeth are 5c's certificate** (the in-M-lag catcher, which the free-identity is blind to
by construction) and **5a's `MatrixInverseOperator.apply_transpose`** (with ZERO negative
tests today, and load-bearing for the swap-law inverse arm).

---

## 5a (D1) — the numerics `CoupledOperator` solve mode (semantics-agnostic)

Home: EXTEND `tests/numerics/test_coupled_operator.py` (the step-4 M1-M5 suite —
synthetic toy blocks, N-general, NOT SN) with a `TestBlockSolve` section + a NEW
`tests/numerics/test_matrix_inverse_operator.py` transpose section (or extend the
existing one). Build synthetic 2×2 grids over two DISTINCT toy System types (asymmetric
dense/diagonal blocks) so triangularity, substitution order, and type-safety are all
hand-verifiable.

- **S1 — triangular detection is STRUCTURAL (None-pattern).** `CoupledOperator.
  is_invertible` on an upper-triangular grid `[[D11, D12], [None, D22]]` (D11, D22
  invertible) is True via the triangular branch; a FULL grid `[[D11, D12], [D21, D22]]`
  is True via the materialize branch (square ∧ below size gate); a grid with a
  NON-invertible diagonal block on the triangular branch refuses. Tooth: a detector that
  reads a full grid as triangular (ignores the present `(B,A)` block) → the substitution
  gives the WRONG answer (it drops the `D21` coupling) — RED against the materialize-LU
  reference `np.linalg.solve(assemble(), b)`.
- **S2 — the substitution VALUE + order (the core new arithmetic).** For an upper-
  triangular `[[D11, D12], [None, D22]]`: `solve([q1, q2])` = `x2 = D22.solve(q2); x1 =
  D11.solve(q1 − D12.apply(x2))`. Assert `array_equal` (toy dense blocks — the block
  substitution adds no reduction) vs a hand-summed reference. **Teeth (each a distinct
  RED):** (i) solve `D11` before `D22` (wrong order — reads an unsolved `x2`); (ii) `q1`
  instead of `q1 − D12·x2` (drop the off-diagonal subtraction); (iii) `+D12·x2` (sign
  flip). Config: **asymmetric blocks D12 ≠ 0, D22 ≠ D11** (a symmetric grid is
  order-transpose-blind, Mode-12) + x = non-flat fixed-seed random.
- **S3 — `solve_transpose` = forward-substitution on the transposed grid.** `Mᵀ` of an
  upper-triangular `M` is LOWER-triangular `[[D11ᵀ, None], [D12ᵀ, D22ᵀ]]`; `solve_
  transpose([b1, b2])` = `x1 = D11.solve_transpose(b1); x2 = D22.solve_transpose(b2 −
  D12.apply_transpose(x1))`. Assert `solve_transpose(b) ≡ np.linalg.solve(dense(M).T, b)`
  (`array_equal` toy dense). Tooth: forward-vs-back substitution swapped on the transpose
  → RED. **THE Mode-12 note:** `solve_transpose` on a SYMMETRIC toy grid is
  transpose-blind — mandate asymmetric D12 ≠ D21-position blocks.
- **S4 — `inverse()` is the STRUCTURE-KEYED factory (taxonomy §13).** On the triangular
  branch `inverse()` returns the substitution wrap (`InverseWrapMixin` — `apply` runs the
  substitution, `solve` = the forward matvec, `inverse()→self`); on the materialize
  branch it returns `MatrixInverseOperator(self)`. Assert `type(grid.inverse())` per branch
  + `grid.inverse().apply(grid.apply(x)) ≈ x` (round-trip, principled-equiv for the
  materialize branch, `array_equal` for the toy triangular).
- **S5 — `MatrixInverseOperator.apply_transpose` (NET-NEW, LUᵀ free).** THE new numerics
  primitive: the transpose of a materialized matrix is the transpose backsolve
  `lu_solve(lu, b, trans=1)` (LAPACK's `trans` flag on the SAME LU factors — no
  re-factorization). Gate: on an **ASYMMETRIC** small dense operator, `MatrixInverse
  Operator(A).apply_transpose(b) ≡ np.linalg.solve(A.as_matrix().T, b)` (`array_equal`/nulp)
  AND `is_adjointable` is now True. **Teeth:** (i) forget the `trans=1` flag (returns
  `A⁻¹b` not `A⁻ᵀb`) → RED on the asymmetric matrix, GREEN on a symmetric one (Mode-12 —
  MANDATE asymmetric); (ii) `is_adjointable` stays False → the swap-law inverse arm (5e)
  cannot construct `grid.H.inverse()`. **This is the arm the whole swap-law rides
  (`_AdjointOperator.is_invertible = invertible(inner) ∧ adjointable(inner.inverse())`,
  operator.py:1251 — `grid.inverse()` MUST be adjointable).** ⚠ verify `InverseWrapMixin`
  does NOT supply `is_adjointable`/`apply_transpose` (operator.py:2058 — "the ADJOINT axis
  is NOT part of the back-half"), so `MatrixInverseOperator` must declare BOTH explicitly.

---

## 5b (D2) — the M re-pose to the honest triangular grid (principled-equiv)

`build_within_group_system` (coupled_system.py:570) currently sets `resolvent =
CoupledInvertibleOperator(LC, space, sn_mesh)` (:700 — the fused-walk facade). D2 re-poses
`M` as the genuine upper-triangular `CoupledOperator [[LC, Seeding], [None, A_BB_march]]`
(the M docstring's block form, :311), whose `solve` is the 5a block back-substitution.
Home: `test_psi_half_coupling.py` (re-point the M-fixtures + row-6 oracle).

- **B1 — the row-6 oracle re-points onto the substitution (THE durable floor, 5.5e-16).**
  `test_extract_to_dense_is_principled_equivalent_not_bit_identical` (:391) compares
  `M_op.solve(q)` vs `np.linalg.solve(_dense(M_op.apply), q)` at `< 1e-11`. After D2,
  `M_op.solve` IS the block substitution (`LC.solve` + `A_BB.solve` + `Seeding` subtract)
  — a DIFFERENT reduction tree from the fused joint walk it replaces → the row-6 bar
  (~5.5e-16 principled-equiv, NOT bit-id) is EXACTLY the M-M row it was designed for. ⚠
  the oracle's `_dense(M_op.apply)` self-probe is TAUTOLOGICAL for the VALUE claim (both
  sides densify the same M) — the row proves solve ≡ LU-of-M (reduction-tree independence,
  LAPACK ⟂ substitution), NOT that M is the right operator; the RIGHT-operator anchor is
  B4 (φ=Q/Σ_t). Keep the row-6 pin; ADD B4.
- **B2 — `test_welded_sweep_is_exact_direct_inverse` (:369) STAYS, doc reword.**
  `M.solve(M.apply(ψ)) ≈ ψ < 1e-12` — after D2 the "welded sweep" is the block
  substitution; the round-trip holds (principled-equiv ~1e-15 ≪ 1e-12). Reword the
  docstring "welded sweep IS the exact inverse" → "block-triangular substitution"; NO
  value re-baseline.
- **B3 — the substitution transpose ≡ the reverse-scan (bit-id INHERITANCE vs the pre-D2
  M).** `M.solve_transpose` after D2 = forward-substitution `x_A = LC.solve_transpose(b_A);
  x_B = A_BB.solve_transpose(b_B − Seedingᵀ·x_A)`. Anchor to a FROZEN pre-D2 baseline
  (`--capture-baseline` on `CoupledInvertibleOperator.solve_transpose` BEFORE the re-pose)
  at the principled-equiv `< 1e-11` bar (both spell the reverse scan; reassociation drift).
  The structurally-independent cross-check: `np.linalg.solve(_dense(M).T, b)` (dense Mᵀ ⟂
  the substitution). Tooth: swap forward-vs-back substitution on the transpose → RED vs the
  dense Mᵀ.
- **B4 — φ=Q/Σ_t END-TO-END on the carrying sphere (the redistribution catcher, the RIGHT-
  operator anchor).** G-d3.1 `test_g_d3_1_flat_flux_equilibrium_on_the_carrying_sphere`
  (psi_half :3701 per the 4e memo — verify the line) STAYS: non-fissile reflective sphere,
  uniform `Q`, uniform `Σ_t` → `φ = Q/Σ_t` everywhere, rtol 1e-10, BOTH drivers. This is a
  PHYSICS invariant the substitution preserves; a dropped/sign-flipped `Seeding` subtraction
  or a wrong A_BB march moves φ off Q/Σ_t O(1). **THE primary 5b correctness anchor** (the
  row-6 oracle is reduction-tree-independence; B4 is the operator-is-correct claim). ⚠
  §0.6: flat-flux `Q/Σ_t` is the single most powerful curvilinear diagnostic AND it nulls
  redistribution — pair it with a HETEROGENEOUS spatial-convergence leg (keff differences
  shrink under refinement) so the redistribution/α-recursion terms activate.
- **B5 — the ρ(M⁻¹N) splitting rate STAYS (`test_outer_si_splitting_rate…` :326).** After
  D2, `_dense(system.resolvent.apply)` probes the block-matvec grid instead of the fused
  facade; ρ = 0.371 holds (~5.5e-16 drift ≪ the O(1) `< c+1e-6` bound). RE-POINT the
  accessor, value HOLDS.

---

## 5c (D5) — the ρ-honest stop (free-identity + certificate)

**The mechanism split (D5 pending ruling — spec BOTH).** The stop test lives in
`SourceIteration.solve` (iteration.py:648-653, today `res = ‖Δψ‖/max(‖ψ‖,1e-30)`; the
docstring :645-647 explicitly warns "switching the stopping norm would re-interpret tol").

- **D5(i) honest-per-iteration** — pass `A` (or a residual callback) to `SourceIteration`;
  evaluate `r = A·ψ − q` via a real forward apply each iteration (`+1 M matvec/iter ≈ 2×
  walk cost`). Gate: the stop value ≡ `‖evaluate_residual(A, ψ, q)‖` each iteration.
- **D5(ii) free algebraic identity + end-certificate (RECOMMENDED)** — inside
  `SourceIteration`, the loop already forms `rhs = q_ext + Σ_g g.apply(ψ)` (= `q + N·ψ_n`).
  The residual of `A = M − N` when `M⁻¹` is exact is `r_{n+1} = N·(ψ_n − ψ_{n+1}) = −N·Δψ`
  (derived: `M·ψ_{n+1} = rhs_n ⟹ A·ψ_{n+1} − q = rhs_n − N·ψ_{n+1} − q = N·(ψ_n − ψ_{n+1})`).
  ZERO marginal cost (retain `rhs_prev`; `r = rhs_prev − rhs`), uniform across ALL arms
  (windowed in moment space). PLUS the honest end-of-solve certificate: ONE driver-level
  `evaluate_residual(system.loss, ψ, q_ext)` after convergence — because the free identity
  is BLIND to in-M lag by construction (#282's seed-lag was IN M, not in N), the certificate
  is the lag-death classifier that closes exactly that hole.

**Scope (D5 pending ruling — spec BOTH).** `coupled arm only` vs `all SourceIteration arms`.
The certificate is stop-AGNOSTIC (it evaluates on the converged flux), so it SHOULD cover
the Krylov path too (GMRES is already residual-stopped, but the certificate is the honest
cross-check). The free-identity stop scope is the pending axis.

**Normalization (D5 pending ruling).** `‖r‖/max(‖q_ext‖, tiny)` (equation-relative) vs
today's `‖Δψ‖/max(‖ψ‖, 1e-30)` (iterate-relative). Equation-relative is the honest choice
(the residual is a rate-density, its natural scale is the source); q_ext≈0 falls back to
absolute (a zero source has a zero solution — refutation r5).

Home: `tests/numerics/test_iteration.py` (the free-identity unit gate — semantics-agnostic)
+ a solve-layer certificate test in `test_psi_half_coupling.py` / a new
`tests/sn/solve/test_coupled_solve_certificate.py`.

- **C1 — the free-identity unit gate (r_running ≡ N·Δψ, semantics-agnostic).** On a SMALL
  synthetic `SourceIteration` (toy `M⁻¹`, toy gain `N`, a few iterations), assert the
  running stop value equals an EXPLICITLY-computed `N.apply(ψ_prev) − N.apply(ψ_new)` at
  each iteration (`array_equal`/nulp on the toy dense). **Teeth:** (i) normalize by `‖ψ‖`
  instead of `‖q_ext‖` → the stop value diverges from the explicit residual O(1); (ii) drop
  a gain from the rhs bookkeeping → `r ≠ N·Δψ`; (iii) sign flip. This is the Mode-12-safe
  gate (the invariance group of `‖·‖` contains only global scaling — an equation residual
  is OUTSIDE it; assert the OBJECT `r` element-wise, not just its norm).
- **C2 — the certificate ≡ `evaluate_residual` when M is exact (the exact-M bridge).** On a
  DIRECT-solve config (M exact, e.g. the carrying sphere with the substitution `M.solve`),
  after convergence assert `‖free_identity_r‖ ≈ ‖evaluate_residual(A, ψ, q)‖` to solver
  tolerance — i.e. the free identity IS the true residual when M⁻¹ is exact. This is the
  bridge that licenses the free-identity as the per-iteration stop.
- **C3 — the in-M-lag CERTIFICATE catcher (THE highest-value NEW row).** Monkeypatch the
  substitution `M.solve` to consume a STALE `ψ_B` (an in-M lag — the #282-family defect: the
  seed march reads a lagged ray). The free-identity stop STILL reports "converged" (it only
  sees N·Δψ, which shrinks) — but the END-CERTIFICATE `evaluate_residual(A, ψ, q)` reports
  `‖r_cold‖ ≫ tol` and raises/flags LOUDLY. **The gate: with the lag mutation active, the
  certificate RED (loud) while the free-identity is GREEN (fooled) — the asymmetry IS the
  lag-death classifier's proof.** Without this, the ρ-honest stop is a Mode-11 vacuous claim
  (the free identity measures the UNCHANGED N sibling; the in-M defect is off its call
  graph). PAIR with the `test_282` C(i) cold-residual (`‖A·solve(b)−b‖<1e-11`, the L18
  acceptance number) which certifies the DIRECT solve consistency this certificate rides on.
- **C4 — the certificate covers the windowed arm via reconstruction (refutation r3).** The
  windowed `M = P @ (L+C).inverse()` (`P` a coisometry — iteration.py:498-500) is NOT an
  exact inverse, so the free-identity `r = N·Δψ` in MOMENT space is an increment-class
  indicator, NOT the true full-system residual. The certificate `evaluate_residual` evaluates
  on the RECONSTRUCTED full-angular flux (`Solution.angular_flux`, which the windowed finalize
  already produces) — so it IS the honest check for the windowed arm. Gate: on a 2-D windowed
  solve, the end-certificate `evaluate_residual` on the reconstructed flux is `< SAFETY×tol`;
  the moment-space free-identity is documented as a rate indicator (NOT credited as the true
  residual for the windowed arm). ⚠ decide at design time: exempt the windowed arm's
  free-identity from the "≡ true residual" claim, keep only its stop role.
- **C5 — zero-gain honest exit (refutation r4).** A zero-gain `SourceIteration` (N=0, `rhs =
  q` constant) → `ψ_1 = M⁻¹q = A⁻¹q` exact; the free identity gives `r_1 = N·Δψ = 0` (exit at
  ≤2 iterations, one-lagged). Gate: zero-gain SI stops in ≤2 iterations AND the certificate
  `evaluate_residual(A, ψ_1, q)` is machine-zero (honest — one M⁻¹ apply IS the exact solve
  when A=M).

---

## 5e — the xfail→live conversion + the swap-law inverse arm

Step 5 flips `CoupledOperator.is_invertible` on the LOSS grid `A` (via 5a's materialize
branch). This XPASSes `test_a2a_grid_swap_law_inverse_arm` (coherence :653-682, currently
`xfail(strict=False)`).

- **E1 — the exact xfail→live edit (vv xfail discipline).** DELETE the `@pytest.mark.xfail(
  strict=False, reason="grid.is_invertible is False until the step-5 block solve lands…")`
  decorator (:653-658). The test body ALREADY asserts `grid.H.inverse().apply(b) ≡
  grid.inverse().H.apply(b)` via `np.testing.assert_array_equal` on BOTH systems (:673-682)
  — a bit-identical OBJECT-IDENTITY gate (both spellings run the SAME `MatrixInverseOperator(
  grid).H` graph through `_AdjointOperator.inverse() = inner.inverse().H`, operator.py:1265;
  L19). Keep `array_equal` (NOT rtol — object identity, bit-identical BY CONSTRUCTION). ⚠ the
  `strict=False` was correct while dormant (a raise → xfail, not a false green); after the
  flip it XPASSes and MUST become a plain live gate (the vv discipline: an xfail for a
  not-yet-landed feature flips to live when the feature lands).
- **E2 — the NET-NEW swap-law preconditions (the arm the whole thing rides).** `grid.H.
  inverse()` requires `grid.H.is_invertible = invertible(grid) ∧ adjointable(grid.inverse())`
  (operator.py:1251). So: (a) `grid.is_invertible` True (5a materialize branch); (b) `grid.
  inverse()` (= `MatrixInverseOperator(grid)`) ADJOINTABLE (5a's `apply_transpose` + `is_
  adjointable=True`). Add a POSITIVE control (`grid.H.inverse().apply(b)` succeeds, correct
  value vs the dense-LU of `dense(grid).T`) BEFORE the bit-identity object gate — a bare
  `array_equal` of two failing constructions is a false green.
- **E3 — the substitution-inverse adjointability chain (M's arm, already LIVE — STAYS).**
  `test_carrying_fused_wrap_refuses_and_the_coupled_sibling_carries` (:479), G1c/G2c/G4c
  (:502/:521/:544), `test_a2a_grid_forward_h_reciprocity` (:618, the M-ADJ-metric tooth) —
  all STAY (they live on M, not the loss grid). Under D2a-retire they RE-POINT (§ retirement).
- **E4 — the two end-to-end value anchors (L2, closed-form structural).** All equivalence
  gates (5a-5b bit-id/principled-equiv) are inheritance = necessary-NOT-sufficient. Anchor
  the SOLVE value structurally-independently: (a) non-fissile reflective sphere fixed-source
  → `φ=Q/Σ_t` (B4 — the redistribution catcher, the single most powerful curvilinear
  diagnostic); (b) fissile A/2g+χ → `k_eff` (≥2G, Cardinal Rule) vs the derivations closed
  form (`test_kinf_homogeneous`) at rtol 1e-10. ⚠ **non-fissile → k_inf is a nan dead gate**
  (refutation r6) — φ=Q/Σ_t for the non-fissile, SEPARATE fissile for k. Mode-12 lens: k
  anchors the FORWARD only (`eig(K)`); the swap-law/A2 reciprocity is the adjoint complement
  (zero adjoint coverage in k, L19). ⚠⚠ the k probes are STOP-TEST-SENSITIVE — state the
  honest tolerance NOW (see the re-baseline ledger); the closed-form anchors at 1e-10 are the
  truth, the 15-digit 4e k-snapshots move and must widen/re-baseline.

---

## THE stop-test-change re-baseline ledger (a deliverable — grep-driven)

The D5 stop-test change re-interprets `tol` (‖Δψ‖-iterate-relative → ‖r‖-equation-relative),
which shifts the converged fixed point by ~`inner_tol`-scale. **The load-bearing evidence:
solver.py:835-837 documents EXACTLY this** — "a schedule change shifts the converged k_eff by
~inner_tol (1e-10-scale — same fixed point, vv Mode 9; only the inner SI stopping differs),
which the keff_tol-tight regression snapshots cannot absorb." The D5 change IS an inner-SI-
stopping change, so the same ledger applies. Partition the existing green gates:

**STOP-TEST-SENSITIVE (must re-baseline or widen):**
- `tests/sn/regression/test_dd_regression.py:100-102` — `result.keff` vs frozen `snap["keff"]`
  via `_regression_assert` (replaced the magic `rtol=1e-12`/`atol=1e-13` floors, :3-8). The
  frozen k moves ~inner_tol-scale → RE-BASELINE the snapshots in the D5 commit (against the
  closed-form anchor, NOT old-vs-new k — Signature 10 / stale-snapshot discipline).
- The 15-digit k probes from 4e (`test_psi_half_coupling` `TestWithinGroupSystemAnchors`
  :3652, and any `assert_allclose(k, <frozen>, rtol=1e-12)` in `tests/sn/solve/`) — widen to
  `SAFETY × inner_tol` (~1e-6) or re-baseline. Grep `tests/sn/solve/ tests/sn/regression/` for
  `keff` + `1e-1[012]`.
- ANY `history.n_inner ==`-frozen-count pin (grep `n_inner ==`) — a stop change moves the
  count. NONE found in `test_si_convergence_rate` (see below), but grep the tree.

**STOP-TEST-ROBUST (assert UNCHANGED — a red here = the change leaked somewhere it shouldn't):**
- `tests/sn/verification/analytical/test_si_convergence_rate.py` — measures `history.n_inner`
  but via RELATIVE ratios (no frozen count): `n_gs < n_jac` (:219), `n_gs >= n_kry` (:224),
  `n_si >= 1.5·n_krylov` (:305), and `n_inner/n_predicted ∈ [0.6, 1.2]` (:280-282, the analytic
  ρ=c anchor). Wide-band relative → ROBUST to a uniform stop change. ⚠ ONE caveat: `test_si_
  jacobi_rate_matches_scattering_ratio` (:266) computes `n_predicted` from the analytic ρ=c
  (`log(tol)/log(c)`) — if the stop re-interprets `tol`, recompute `n_predicted` CONSISTENTLY
  with the new stop semantics (the [0.6, 1.2] band absorbs ±1 iteration on counts of ~10-40,
  but the CENTER must track the new stop). This is a convergence-RATE claim (flux-shape-
  INDEPENDENT, 1G-legitimate per AGENT.md §5) — DISTINCT from the value claims; do NOT
  re-baseline it, re-derive `n_predicted`.
- The frozen `walk_matvec_*` baselines, the affine-carve bit-id, the row-6 oracle, `test_282`
  C(i)-C(iv), the coherence G2c `array_equal`, `test_recovery_preserves_kinf_2g` (:315,
  closed-form-anchored) — all evaluate on FIXED inputs / closed-form references, NOT on the
  iteration's stopping. STAY unchanged; a red = a real regression.

**The ledger IS a deliverable:** the D5 commit must (a) re-baseline the frozen k-snapshots
against the closed-form anchor, (b) re-derive `n_predicted` in the rate test, (c) leave the
robust set untouched and run it as the tripwire (any robust red = the stop change leaked into
a value/structural path).

---

## The migration / retirement ledger (D2a / D3 branch-spec'd)

**D2a (PENDING) — retire `CoupledInvertibleOperator`/`CoupledSweepOperator` vs keep as SN
facades.** Inventory (`grep -rn "CoupledInvertibleOperator\|CoupledSweepOperator" tests/`,
verified at tip): 7 files.

| test file | refs | D2a-RETIRE fate | D2a-KEEP fate |
|---|---|---|---|
| `test_si_single_primitive_contract.py` | 7 | RE-POINT: `isinstance(L_eig, CoupledSweepOperator)` / `isinstance(L_eig.inner, CoupledInvertibleOperator)` (:168-170) → the new triangular-`CoupledOperator` inverse type (an `InverseWrapMixin` substitution wrap) | STAYS (facades persist) |
| `test_inverse_adjoint_coherence.py` | 8 | RE-POINT `_coupled_M` (:449) to build the triangular `CoupledOperator` grid; G4c's Mode-11 wrap (:544, `monkeypatch.setattr(CoupledInvertibleOperator, "solve"…)`) RE-AIMS onto the grid's `solve` — else silently toothless | STAYS |
| `test_psi_half_coupling.py` | 11 | RE-POINT `_joint_M` (:206) → the triangular grid; row-6 oracle (:391) re-points onto the grid `.solve` | STAYS (`_joint_M` builds the facade) |
| `test_282_direct_seed_fixed_point.py` | 2 | verify the Mode-11 sentinels (`test_mode11_direct_solver_executes_on_sphere_not_slab`) re-aim onto the successor solve — a deleted symbol → silent no-op | STAYS |
| `test_radial_characteristic_metric.py` | 2 | audit | STAYS |
| `test_krylov_curvilinear_precond_safety.py` | 2 | audit | STAYS |
| `test_phase_c_gates.py` | 2 | audit | STAYS |

**RETIRE is the aggressive-retirement-rule default** (coding-standards: superseded code is
noise; M IS a numerics `CoupledOperator`, so the SN facade is a same-job-done-worse
predecessor). The RE-POINT is mechanical (the fixtures build the successor grid). ⚠ the
retirement audit's THREE searches: (1) graph `callers` on `CoupledInvertibleOperator`/
`CoupledSweepOperator`; (2) text-grep code + tests + `docs/` (an unresolved `:class:` renders
as plain text with NO `-W` warning — grep `docs/theory/` for the retired names); (3) direct
constructors (`CoupledInvertibleOperator(` / `CoupledSweepOperator(` call sites — production
`_within_group_si` :739 `isinstance(system.resolvent, CoupledInvertibleOperator)`, and
`build_within_group_system` :700). RECOMMENDATION: **D2a-RETIRE** (the algebra says M IS a
numerics CoupledOperator; a named SN facade over it is a twin-path habitat) — but FLAG the
`_within_group_si` :739 isinstance dispatch must re-key on the triangular-grid resolvent type.

**D3 (PENDING) — the :712 fused-matvec twin.** `orpheus/sn/operators/radial_characteristic.py
:712` (the tracked-transient-twin note — the A_AB seed rows in the FORWARD matvec).
- **D3-ROUTE** (route `M.apply` through the grid block matvec too): under D2a-RETIRE this is
  FORCED (a `CoupledOperator` grid's `apply` IS the block matvec). The Krylov forward-matvec
  FP stream moves ~5.5e-16-class → RE-BASELINE the frozen `walk_matvec_sphere_2g.npz`
  fwd_bulk (anchored to `_dense(grid.apply)`, NOT old-vs-new ULP — L18). Retires the
  `_seed_rows_forward/_transpose` inline placements' production callers.
- **D3-KEEP** (keep the fused matvec): the twin persists behaviour-pinned; `M.apply` stays
  the fused walk, only `M.solve` re-poses. More bit-identity survives (the sphere fwd_bulk
  STAYS bit-exact). RECOMMENDATION: **D3-KEEP under H1-narrow** (extract the genuinely-
  separable A_BB solve; leave the forward matvec fused — it is already single-sourced via the
  shared kernel, and forcing the block matvec buys reassociation drift for no architectural
  gain), UNLESS D2a-RETIRE forces D3-ROUTE (then re-baseline the sphere fwd_bulk).

**`TestCoupledSolve` home:** `test_psi_half_coupling.py` (the ψ½ instance's suite — where
`TestRegressionFloor`/`TestCoupledBuilder`/`TestWithinGroupSystem` live) for the SN-bound
rows (B1-B5, C2-C5, E3-E4); `tests/numerics/test_coupled_operator.py` for the semantics-
agnostic 5a (S1-S4); `tests/numerics/test_matrix_inverse_operator.py` for S5; `tests/numerics/
test_iteration.py` for C1; `test_inverse_adjoint_coherence.py` for E1-E2.

---

## The singular-grid loud-failure leg (deliverable #4 — I AGREE with the numpy convention)

`is_invertible` advertising the ROUTE while singularity fails loudly at `inverse()`
construction IS the right design (numpy.linalg.solve convention), and `MatrixInverseOperator`
ALREADY implements it: (a) exact singularity → `LinAlgError` at CONSTRUCTION (zero LU pivot,
matrix_inverse_operator.py:159-166); (b) rectangular → `ValueError` (:140-146); (c) too-large
→ `MatrixTooLarge` (eager, :137-138). The `is_invertible` predicate is STRUCTURAL (reads the
operand tree — square ∧ materializable), and singularity is a VALUE property caught eagerly at
`inverse()` (the docstring's `(−S)+(L+C)` witness, :38-52). **Gate:** a square-assemblable but
SINGULAR grid (e.g. a gain grid `N` config with a zero-row diagonal block, or a degenerate
2×2 where a diagonal block is singular) → `grid.inverse()` raises `LinAlgError` LOUDLY at
construction, NOT a garbage backsolve. Assert `pytest.raises(LinAlgError, match="exactly
singular")`. Tooth: a route that returns an `inf`/`nan`-bearing backsolve instead of raising →
RED. I do NOT disagree with the convention — advertise the route structurally, fail-loud on
the value eagerly; it matches the existing `MatrixInverseOperator` eager-guard design and
Cardinal Rule 1 (never return a non-inverse).

---

## Mutation ledger (monkeypatch-only; NEVER `git checkout` an uncommitted file — L4)

| gate | tooth (in-process) | target `file:line` (⚠ best-estimate) | expected RED |
|---|---|---|---|
| S1 | read a full grid as triangular (ignore the `(B,A)` block) | NEW `CoupledOperator._is_triangular` | wrong answer vs materialize-LU ref |
| S2 (i/ii/iii) | swap solve order / drop `−D12·x2` / sign-flip `+D12·x2` | NEW `CoupledOperator.solve` substitution | `array_equal` fail vs hand-sum (asymmetric blocks) |
| S3 | forward↔back substitution swap on Mᵀ | NEW `CoupledOperator.solve_transpose` | fail vs `solve(dense(M).T, b)` (asymmetric) |
| S5 (i) | forget `trans=1` (return `A⁻¹b` not `A⁻ᵀb`) | NEW `MatrixInverseOperator.apply_transpose` | RED asymmetric / GREEN symmetric (Mode-12) |
| S5 (ii) | `is_adjointable` stays False | NEW `MatrixInverseOperator.is_adjointable` | swap-law inverse arm cannot construct |
| B1 | naive dense `M⁻¹` ignoring inflow/seed rows | (the row-6 gate's own reference) | `> 1e-11` |
| B4 | drop / sign-flip the `Seeding` subtraction in the substitution | NEW `M.solve` (`build_within_group_system` re-pose) | φ moves off `Q/Σ_t` O(1) |
| C1 | normalize by `‖ψ‖` not `‖q‖` / drop a gain | NEW `SourceIteration` free-identity | stop-value ≠ `N·Δψ` O(1) |
| C3 (HEADLINE) | `M.solve` consumes a STALE `ψ_B` (in-M lag) | NEW `M.solve` (or `A_BB.solve` staleness) | certificate `‖r_cold‖ ≫ tol` LOUD while free-identity GREEN |
| singular | route returns `inf/nan` backsolve not raise | `matrix_inverse_operator.py:159` | no `LinAlgError` → RED |
| E2 | `grid.inverse()` not adjointable | `MatrixInverseOperator.is_adjointable` | `grid.H.inverse()` construction raises |

---

## Refutations (fire before the ink dries)

1. **(r1) Is the substitution REALLY principled-equiv on the REFLECTIVE corner config?**
   `B_b` is in `N` (the gains), so the corner datum arrives via `rhs_B`'s corner slot in BOTH
   routes (fused threads it via `radial_characteristic_source`'s corner slot, solver.py:2077-
   2084; the substitution's `A_BB.solve(rhs_B)` reads `rhs_B`'s corner). During SI the corner
   reflection is the LAGGED `B_b ∈ N` — `rhs_B.corner = q_B.corner + (B_b·ψ_prev).corner`.
   VERIFY no hidden data-flow difference (the d3-era product-cyl ig-hazard analog): a gate that
   drives BOTH `M.solve` (substitution) and the pre-D2 fused `M.solve` on a reflective-sphere
   `rhs_B` with a NONZERO corner, asserting `< 1e-11` on ψ_B's corner block. If the fused route
   needs the manual corner-threading choreography (:2077-2084) but the substitution's `A_BB.
   solve` reads the corner from `rhs_B` directly, CONFIRM they agree (they should — the reflect
   is lagged into N in both). This is a MANDATORY 5b gate, not an assumption.
2. **(r2) Does the free-identity stop move a keff_tol-tight snapshot beyond its band?** ρ≈0.371
   → `1/(1−ρ)≈1.59`, so the increment understates the true error ~1.59× — the residual stop
   shifts the converged k by ~inner_tol-scale (solver.py:835-837). QUANTIFY: the shift is
   ~1e-8-to-1e-10, which is BELOW `keff_tol=1e-7` (outer convergence unaffected) but ABOVE the
   frozen 15-digit k-snapshots' band → the snapshots (`test_dd_regression:100`) RE-BASELINE; the
   closed-form anchors at rtol 1e-10 ABSORB it (borderline — a 1e-10 shift is at the 1e-10 band
   edge; state the honest anchor tol as `SAFETY×inner_tol ≈ 1e-8` for the value probes, 1e-10
   for the analytic `Q/Σ_t`/k_inf).
3. **(r3) The windowed-moment arm under scope=all — is `‖rhs_n − rhs_{n+1}‖` a legitimate
   equation residual there?** NO — the windowed `M = P @ (L+C).inverse()` (`P` a coisometry,
   iteration.py:498-500) is NOT an exact inverse, so the exact-M premise FAILS and `r = N·Δψ`
   in moment space is an increment-class indicator, not the true residual. RESOLUTION (C4):
   exempt the windowed arm's free-identity from the "≡ true residual" claim (keep only its stop
   role); the end-CERTIFICATE `evaluate_residual` on the RECONSTRUCTED full-angular flux (which
   the windowed finalize already produces) IS the honest check and covers the windowed arm.
4. **(r4) Zero-gain SI — is the exit honest?** YES (C5): `ψ_1 = M⁻¹q = A⁻¹q` exact (A=M when
   N=0); the free identity gives `r_1 = N·Δψ = 0` (exit ≤2 iterations); the certificate confirms
   `evaluate_residual(A, ψ_1, q)` machine-zero. Gate it (≤2 iterations + machine-zero
   certificate) — a zero-gain path that iterates or exits non-zero is a bookkeeping bug.
5. **(r5) The q_ext≈0 edge (the normalization guard).** `‖r‖/max(‖q_ext‖, tiny)` — for q_ext≈0
   the guard caps the relative residual at absolute (a zero source ⟹ zero solution; the same
   guard family as today's `max(‖ψ‖, 1e-30)`). Gate: q_ext exactly zero → the solve returns ψ=0
   with a machine-zero absolute certificate (no div-by-zero, no spurious non-convergence).
   Document the fallback (equation-relative degrades to absolute as q_ext→0 — honest).
6. **(r6) non-fissile → k_inf is a nan dead gate** — E4 uses fixed-source `φ=Q/Σ_t` for the
   non-fissile sphere + a SEPARATE fissile A/2g for k, ≥2G (Cardinal Rule).
7. **(mine — the headline) `grid.is_invertible` is the LOSS grid A, not M.** M is ALREADY
   invertible (`CoupledInvertibleOperator.is_invertible=True`, :444). The dormant xfail flips on
   the LOSS grid's materialize-LU route, requiring `MatrixInverseOperator.apply_transpose`
   (NET-NEW, zero negative tests) for `grid.inverse()` to be adjointable
   (`_AdjointOperator.is_invertible` clause 2, operator.py:1251). Do NOT spec the flip as an M
   property — it is a `CoupledOperator` (numerics) property realized via `MatrixInverseOperator`.
8. **(mine) the in-M-lag certificate is NOT redundant with the row-6 oracle.** The row-6 oracle
   (`M.solve ≡ dense-LU of M`) proves the substitution is a faithful inverse of the ASSEMBLED M
   — but if the PRODUCTION M carries an in-M lag (the seed march reads a stale ray, #282
   family), the assembled `_dense(M.apply)` carries the SAME lag on both sides (tautology-blind,
   L18 snapshot-inheritance). The certificate `evaluate_residual(A, ψ_converged, q)` evaluates
   the TRUE `A·ψ − q` (A includes the un-lagged coupling) — it is the ONLY gate that catches an
   in-M lag. C3 is mandatory; the row-6 oracle does NOT subsume it.
9. **(mine) the certificate needs the COUPLED q_ext.** `evaluate_residual(system.loss, ψ, q_ext)`
   requires `q_ext` as a `CoupledField [q_A, q_B]` on a carrying mesh (`_build_fixed_source_rhs`
   :2138 builds it; eigenvalue q_ext = coupled fission source). The certificate lives at the
   SOLVER layer (where `system.loss` + the coupled q_ext are in scope), NOT in `SourceIteration`
   (which is semantics-agnostic and holds only M⁻¹ + N). This keeps the free-identity in
   `SourceIteration` (needs only gains) and the certificate in the solver (needs A + coupled q).

## Result contract

Land in FIVE sub-commits (5a numerics-solve D1 / 5b M-re-pose D2 / 5c ρ-honest-stop D5 / 5d
retirement D2a+D3 / 5e xfail→live + anchors + re-baseline), sequencing verification-first: **5a
FIRST (semantics-agnostic, synthetic — the substitution + `MatrixInverseOperator.apply_transpose`
proven in isolation), THEN 5b (M re-pose, row-6 principled-equiv re-point), THEN 5c (the stop,
with C3 the in-M-lag certificate as the headline NEW row), THEN 5d (retirement, D2a/D3 ruled),
THEN 5e (flip the xfail + re-baseline the k-snapshots)**. NEW: `tests/numerics/test_coupled_
operator.py::TestBlockSolve` (S1-S4) + `test_matrix_inverse_operator.py` transpose (S5) +
`tests/numerics/test_iteration.py` free-identity (C1) + a solve-layer certificate test (C2-C5) +
`test_psi_half_coupling.py::TestCoupledSolve` (B-rows, E3-E4). CONVERT: `test_a2a_grid_swap_law_
inverse_arm` (coherence :653) xfail→live + the E2 positive control. RE-POINT (D2a-RETIRE):
`_joint_M`/`_coupled_M`/the row-6 oracle/`test_si_single_primitive_contract` isinstance rows/
`_within_group_si` :739 dispatch. RE-BASELINE: the frozen k-snapshots (`test_dd_regression:100`,
against the closed-form anchor); RE-DERIVE `n_predicted` (`test_si_convergence_rate:266`, NOT
re-baseline). STAYS (do NOT touch): `test_si_convergence_rate` relative ratios, the frozen
`walk_matvec_*`/affine-carve/row-6/`test_282` C(i)-C(iv)/G2c `array_equal` (run as the leak
tripwire — a robust red = the stop change leaked into a value/structural path). **Load-bearing
deliverables:** the **row-6 dense-LU oracle** (M.solve ≡ LU-of-M, the 5.5e-16 floor the
substitution re-points onto) + **the in-M-lag CERTIFICATE C3** (the free-identity is blind to
in-M lag BY CONSTRUCTION — the certificate is the ONLY catcher, #282's lag was IN M) + **the
φ=Q/Σ_t + k_eff anchors E4** (the operator-is-correct claim, not just reduction-tree-independent)
+ **`MatrixInverseOperator.apply_transpose` S5** (the arm the swap-law inverse rides — NET-NEW,
Mode-12-asymmetric) + **the stop-test re-baseline ledger** (the frozen k-snapshots move; the ρ=c
ratios are robust). End-to-end acceptance per sub-commit: full `tests/sn -m "not slow"` +
`tests/numerics` + `tests/diffusion` + `tests/cp` + ratchet `transport:1` + sphinx -W; through 5a
the wall STAYS 6340/0, and after 5e only the frozen k-snapshots moved (re-baselined) + one xfail
went live (37→36 xfailed) — every count else holds.
