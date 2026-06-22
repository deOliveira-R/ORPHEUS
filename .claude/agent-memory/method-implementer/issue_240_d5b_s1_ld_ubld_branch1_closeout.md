# D5b-S1 Branch 1 — d-generic UBLD SymPy algebra-of-record (closeout)

**Issue:** #240 / #38 / #37, sub-step D5b-S1 **Branch 1 only**.
**Branch:** `feature/sn-space-angle-tier2`, MAIN checkout (NOT a worktree).
**Date:** 2026-06-16. **Status:** Branch 1 COMPLETE; NOT committed (user commits).
**Ships nothing to production** — symbolic reference + proofs only. Branch 2
(numpy production primitive) is a SEPARATE follow-on dispatch.

## What this is

The canonical Branch-1 (SymPy, State 1A closed-form) algebra-of-record for
the dimension-generic **Unlumped BiLinear Discontinuous** (UBLD) per-cell
Galerkin system of the SN equations on a Cartesian cell — the tensor-product
bilinear/trilinear DG-P1 object (`{1,x,y,xy}`, 2^d moments), NOT the
simplex-P1 `{1,x,y}` (the Adams-2001 thick-diffusion-limit trap).

## Symbolic structure (the synthesis)

Built via Kronecker products of the verified 1-D LD factor operators in the
Legendre moment basis `{1, P₁}` on width `h` (reverse-engineered from the
production slab 2×2 + confirmed by the d=1 reduction oracle):

- **Mass** `M_1d = diag(h, θh)` → `M = M₁ ⊗ … ⊗ M_d`.
- **Streaming** `G_1d = |μ|·[[0,0],[-2,0]]` (∫B_i ∂_x B_j) →
  `G = Σ_a μ_a·(M₁ ⊗ … ⊗ G_1d[axis a] ⊗ … ⊗ M_d)` — gradient on the active
  axis, **mass on every transverse axis** (the volume-integral
  factorisation; this is the load-bearing build choice).
- **Downstream-face** `F_out_1d = |μ|·outer(B(+1),B(+1)) = |μ|·[[1,1],[1,1]]`
  → `F_out = Σ_a μ_a·(M₁ ⊗ … ⊗ F_out_1d[axis a] ⊗ … ⊗ M_d)`.
- **Upstream-face** (inflow → RHS): `B(-1)=[1,-1]` test-weighting on the
  active axis, mass on transverse axes, times `|μ_axis|`; the upstream face
  is a `2^{d-1}`-moment transverse object.
- Assembled: `A = G + F_out + Σ_t·M` (2^d×2^d dense non-symmetric);
  `R = M·S⃗ + Σ_axes F_in·ψ_in_traces`.

The d=1 case is the Kronecker-with-one-factor identity → reduces EXACTLY to
the production slab 2×2. The `xy` cross-moment coupling EMERGES from the
algebra at d≥2 (the `θ²` diagonal weight on the xy moment, `θ³` on xyz);
NO 4×4/8×8 entry hand-transcribed. Confirmed: d=1 → production 2×2 to
symbolic zero; d=2 → 4×4 with xy coupling; d=3 → 8×8 with θ³.

## Deliverables (file paths)

1. **SymPy assembler:**
   `orpheus/derivations/discrete/sn/ld_ubld.py` — `assemble_ubld` (Kronecker
   build), `assemble_inflow_axis`, `per_cell_solve`, `downstream_face_trace`,
   + 5 `derive_*()` verification functions. Home matches the sibling
   `balance.py` (discrete-SN symbolic discretisation the production solver
   must satisfy).
2. **Foundation gate:**
   `tests/sn/spatial/test_ld_ubld_symbolic.py` — 6 `@pytest.mark.foundation`
   tests (no `verifies(...)`), `-O`-safe (function-call assertions via
   `_require`/`_require_zero`/`_require_zero_matrix` + `np.testing`, NOT bare
   assert — Mode 8). One test per claim + an extra anchoring the symbolic
   primitive to the LIVE production `LinearDiscontinuous.update`.
3. **Sphinx stub:** `docs/theory/discrete_ordinates.rst` §`ld-ubld-multidim`
   (after `sn-affine-outgoing-face-reconstruction`) — labels
   `ld-ubld-cell-system`, `ld-ubld-d1-reduction`, `ld-ubld-exact-on-bilinear`
   + `:mod:` cross-refs + a `.. todo:: Archivist expansion` marker. Sphinx
   builds clean (exit 0; the only 6 warnings are PRE-EXISTING SyntaxWarnings
   in unrelated test files).
4. **ERR-060** logged in `.claude/skills/vv-principles/error_catalog.md`
   (the caught bug — see below). `@pytest.mark.catches("ERR-060")` on the d2
   test.
5. This closeout memo.

## The two oracles (all proven `sympy.simplify(diff) == 0`)

**Oracle (i) — d=1 reduction.** `derive_d1_reduction_to_production`:
`diff_A = [[0,0],[0,0]]`, `diff_R = [[0],[0]]`, `diff_face = 0`,
`diff_psi_bar = 0`, `diff_psi_hat = 0`. `S` and `D₂'` recovered as the
production closed forms:
`S = (μ² + (Σ_t h + μ)(Σ_t h θ + μ))/(Σ_t h θ + μ)`, `D₂' = Σ_t h θ + μ`.
The ÷V view (`derive_d1_kernel_view_equals` — production `_kernel_terms`)
and ×V view (`derive_d1_scan_view_equals` — production
`affine_scan_coefficients`, transmission `a` source-independent) both reduce
to the same d=1 with `diff_psi_bar = diff_psi_out = 0` — the "single-source
the math" proof for Branch 2.

**Oracle (ii) — exact-on-bilinear (d=2).** `derive_d2_exact_on_bilinear`:
ψ = a+bx+cy+dxy fed via DG-exact upstream-x/upstream-y face moments +
projected source moments → solved 4 moments − exact projections =
`Matrix([[0],[0],[0],[0]])`. The xy cross moment is exercised (d≠0 symbolic).

## VERIFICATION PASTE-BACK (L12)

Final clean gate, canonical `python -O` invocation:

```
$ .venv/bin/python -O -m pytest tests/sn/spatial/test_ld_ubld_symbolic.py -v
collected 6 items
tests/sn/spatial/test_ld_ubld_symbolic.py::TestOracleId1Reduction::test_d1_reduction_to_production_schur PASSED [ 16%]
tests/sn/spatial/test_ld_ubld_symbolic.py::TestOracleId1Reduction::test_d1_divV_kernel_view_equals_reduction PASSED [ 33%]
tests/sn/spatial/test_ld_ubld_symbolic.py::TestOracleId1Reduction::test_d1_timesV_scan_view_equals_reduction PASSED [ 50%]
tests/sn/spatial/test_ld_ubld_symbolic.py::TestOracleId1Reduction::test_d1_symbolic_primitive_matches_production_update PASSED [ 66%]
tests/sn/spatial/test_ld_ubld_symbolic.py::TestOracleIIBilinearExactness::test_d2_exact_on_bilinear PASSED [ 83%]
tests/sn/spatial/test_ld_ubld_symbolic.py::TestD3StructuralReadiness::test_d3_assembles_8x8_with_theta_cubed PASSED [100%]
======================== 6 passed, 1 warning in 20.61s =========================
```

LD floor + new gate combined: `25 passed, 1 warning` (19 existing
`test_linear_discontinuous.py` + 6 new) — no regression.

Key `sympy.simplify(...) == 0` outputs (from the module CLI):
```
[PASS] V_d1: diff_A = Matrix([[0,0],[0,0]]); diff_R = Matrix([[0],[0]]);
       diff_face = 0; diff_psi_bar = 0; diff_psi_hat = 0;
       S = (mu**2 + (Sigma_t*h + mu)*(Sigma_t*h*theta + mu))/(Sigma_t*h*theta + mu);
       D2_prime = Sigma_t*h*theta + mu
[PASS] V_d1_kernel: diff_psi_bar = 0; diff_psi_out = 0
[PASS] V_d1_scan:   diff_psi_bar = 0; diff_psi_out = 0; a_source_independent = True
[PASS] V_d2_bilinear: diff = Matrix([[0],[0],[0],[0]])
[PASS] V_d3: size = (8,8); xyz_collision_weight = theta**3
```

## ERR-060 — the caught bug (the oracle did its job)

First draft of `assemble_inflow_axis` dropped the `|μ_axis|` streaming factor
on the inflow RHS (Mode #3 missing factor). INVISIBLE to all three d=1
oracles (the d=1 RHS is built inline in `derive_d1_reduction_to_production`,
not routed through `assemble_inflow_axis`); CAUGHT by Oracle (ii) — the d=2
exact-on-bilinear gate is the FIRST consumer of the multi-axis inflow
assembler. **Mutation-verified `-O`-safe:** re-dropping the factor → d2 test
FAILS (returncode 1, via `pytest.fail`), d1 tests stay GREEN (proof they're
blind to the bug class). Fixed: `return mu_axis * sp.Matrix(...)`. This is
the algebra-of-record discipline working as designed — bug never reached
production.

## Notes / decisions

- **State 1A** chosen (closed-form): the per-cell solve is a finite dense
  symbolic LUsolve; all claims close as `simplify(diff)==0`; the slowest is
  the d2 exact-on-bilinear (~15-18 s, the closed-form double integration of
  the manufactured source/trace moments) — within the 30 s State-1A budget.
- **Home:** `orpheus/derivations/discrete/sn/` (the discrete-SN symbolic
  discretisation the production solver must satisfy) — sibling to
  `balance.py`/`contamination.py`. NOT `continuous/` (that's references).
- **`assemble_inflow_axis`** supports `axis ∈ {0, d-1}` (all the oracles
  need; d=2 uses axis 0 and 1). The general interior-axis interleave
  (d≥3, 0<axis<d-1) raises NotImplementedError — Branch-2 territory.
- **Structural-frame check** (cross-domain-frames, before discretization):
  the Kronecker/tensor-product algebra IS the native frame (the contract
  prescribed it); elegance detector fired NO smells; no cross-domain-attacker
  dispatch needed.
- **NOT touched:** any `orpheus/sn/` production file, `affine_scan_coefficients`
  / `_kernel_terms` / `_schur_terms`. NO new runtime types. NO git add.

## Owed (downstream)

- **Branch 2** (the numpy production primitive for the d≥2 UBLD
  `cell_kernel_batch`/`residual_kernel_batch`) — separate dispatch; the
  face-payload contract widening (`2^{d-1}`-moment face cochain, literature
  memo §3) is its load-bearing architecture item.
- **Archivist** DISPATCH_REQUEST emitted (the rich narrative for the
  `ld-ubld-multidim` stub) — `followup: false`.
