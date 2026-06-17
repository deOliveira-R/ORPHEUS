# D5b-S1 Branch 2 — production numpy d-generic UBLD primitive + Rule-of-Three collapse (closeout)

**Issue:** #240 / #38 / #37, sub-step D5b-S1 **Branch 2**.
**Branch:** `feature/sn-space-angle-tier2`, MAIN checkout (NOT a worktree).
**Date:** 2026-06-16. **Status:** Branch 2 COMPLETE; NOT committed (user commits + runs review agents).
**Bit-identity-sensitive refactor of SHIPPED production code, done in ISOLATION (no d≥2 wiring).**

## What this is

The numpy Branch-2 counterpart of the committed Branch-1 SymPy
algebra-of-record (`orpheus/derivations/discrete/sn/ld_ubld.py`). Two layers,
ONE source of truth, plus the collapse of the three production 1-D LD views
(the Rule-of-Three) onto a single shared d=1 closed-form helper.

## Deliverables (file paths)

1. **NEW production primitive:** `orpheus/sn/spatial/_ubld.py`
   - `assemble_ubld` / `assemble_inflow_axis` / `per_cell_solve` — the
     d-generic dense Kronecker primitive, BATCHED over a cell/anti-diagonal
     stack (`np.linalg.solve` over the trailing 2 axes). Numpy mirror of the
     symbolic assembler one-for-one (factors `_GRAD_1D=[[0,0],[-2,0]]`,
     `_FOUT_1D=[[1,1],[1,1]]`, `_FIN_TRACE=[1,-1]`, `mass_1d=diag(h,θh)`;
     `_batched_kron` is the numpy batched analog of `sympy.kronecker_product`).
     CANONICAL d-generic source: d=1 (now) AND d≥2 (S2 wires onto it). NOT the
     production d=1 path (that would be the L16 per-cell-solve regression).
   - `d1_closed_form(g, sig_t, theta) -> D1ClosedForm` — the analytic Schur of
     the primitive's d=1 2×2, VECTORIZED (no dense solve) = the production fast
     path. `D1ClosedForm` carries the SCALE-FREE invariants
     (`g`, `g_over_theta`, `d2`, `k`, `w`, `eff_denom`, `theta`) + three view
     methods: `kernel_rhs` (÷V), `scan_xV` (×V transmission `(a,inv_denom,w)`),
     `schur_xV` (×V per-cell `(S, eff_source, eff_numer, θ·ŝ, |μ|A_down, D₂')`).
2. **COLLAPSED production views:** `orpheus/sn/spatial/linear_discontinuous.py`
   — all THREE views now single-source through `d1_closed_form` (imports ONLY
   `d1_closed_form`, NOT `per_cell_solve`/`np.linalg.solve` — L16 audit clean):
   - `_schur_terms` (×V per-cell) → `cf.schur_xV(h, s_bar, s_hat, psi_in)`.
   - `_kernel_terms` (÷V DAG kernel) → `cf.eff_denom`, `cf.kernel_rhs`, `cf.w`.
   - `affine_scan_coefficients` (×V scan) → `cf.scan_xV(V_full)`.
3. **NEW Branch-2 gate:** `tests/sn/spatial/test_ld_ubld_primitive.py` — 10
   `@pytest.mark.foundation` tests, `-O`-safe (`np.testing` / `pytest.fail`,
   the `_require` helper; NO bare assert — Mode 8). Three groups:
   - `TestPrimitiveMatchesSymbolic` (4): numpy d=1 A/M/G/F_out == symbolic
     (small input); d=1 solved moments == symbolic dense; batched-over-stack
     shape contract; **d=2 exact-on-bilinear** (`@catches("ERR-060")`, the
     xy-coupling structural-independence gate — the numpy analog of the
     symbolic Oracle ii).
   - `TestClosedFormEqualsDenseReduction` (3): the shared closed form == the
     dense d=1 `per_cell_solve` in ALL THREE views (÷V kernel, ×V scan, ×V
     per-cell Schur incl. slope) to machine ε. **This closes the
     elegance-enforcer's Branch-1 CONCERN in code.**
   - `TestProductionViewsAnchoredToPrimitive` (3): the LIVE production scheme
     (`update` / `cell_kernel_batch` / `affine_scan_coefficients`) ==
     the dense primitive's d=1 solve, 2G heterogeneous (the LINK PROOF —
     anchors the production views to the d-generic primitive, not only to each
     other).
4. **Sphinx stub:** `docs/theory/discrete_ordinates.rst` §`ld-ubld-branch2-primitive`
   (after the Branch-1 `ld-ubld-multidim` todo) — labels
   `ld-ubld-scale-free-invariants` + `ld-ubld-rule-of-three-collapse`, `:mod:`
   cross-refs to `_ubld`, a `.. todo:: Archivist expansion (D5b-S1 Branch 2)`.
   Sphinx builds clean (exit 0; only PRE-EXISTING warnings: 6 SyntaxWarnings in
   unrelated test files + the pre-existing `ld-cartesian-1d`/`ld-slab`
   verifies-marker info messages from `test_mms_ld_slab.py`, untouched by me).
5. This closeout memo.

## The single-source synthesis (Rule-of-Three collapse)

The d=1 Schur rides two SCALE-FREE (×V≡÷V) dimensionless invariants:
`k = (g/θ)/(g/θ + Σ_t)`, `w = 1/(1+k)`, with `g = |μ|A_down/V` the ÷V
streaming-over-volume. Every view's coefficients are an algebraic function of
`(g, Σ_t, k, w)` × a power of `V`. The crosswalk that made it work:
- ÷V Schur diagonal `eff_denom = (g+Σ_t) + g·k`; ×V `S = V·eff_denom`.
- ÷V slope denom `d2 = g/θ + Σ_t`; ×V production `D₂' = θ·V·d2 = Σ_t·θ·V + |μ|A_down`
  (THE bug-prone scaling — `D₂'` carries the extra `θ` the ÷V `d2` does NOT;
  caught at first write by the `schur_xV` link test, fixed `V·d2`→`θ·V·d2`).
- ×V transmission `a = m(1+k)²/S − k` with `m = |μ|A_down = g·V` (source-indep).

## Bit-identity decision (the load-bearing V&V call)

Single-sourcing means ONE reduction tree; the three views used DIFFERENT
reduction trees before, so re-baselining ≥2 of them is INHERENT. The decision:
- **The strict gate `tests/sn/sweep/core tests/sn/solve` carries NO LD
  numerical golden** — its LD items are STRUCTURAL only (routing / trait probes
  / "refuses 2d LD"); its numerical snapshots are DD-only = the BIT-IDENTICAL
  NEGATIVE CONTROL. It stayed `513 passed / 1 skipped / 4 xfailed` pre==post
  (if it had moved, I'd have broken the 1-D DD math — it did not).
- **LD numerical behavior lives in the LD-specific gates** at `rtol=1e-12`
  (`test_ld_two_paths_scan_equals_dag_oracle` rtol 1e-12; MMS conv rtol 1e-9),
  which ABSORB the ~1-ULP re-baseline by construction. `30 passed / 1 xfailed`
  pre==post. The re-baseline is the helper's `g·k` / `m·k` associations vs the
  legacy inline `(g·g/θ)/d2` / `m·p/d2` — pure FP-non-associativity, bounded by
  (reduction depth)×ULP. Sanctioned by vv-principles 3-criteria (named
  intermediates `k`/`w`/`eff_denom`; verified vs the structurally-independent
  dense primitive; drift = FP association). No golden file moved (no LD `.npy`).

## VERIFICATION PASTE-BACK (L12) — literal stdout in the report message.

Strict gate `513 passed, 1 skipped, 4 xfailed` (pre==post, bit-identical
negative control held). LD gates `30 passed, 1 xfailed` (pre==post). New
primitive+link gate `10 passed`. Sphinx exit 0.

## ERR-060 mutation re-verified (the d=2 gate has teeth, -O-safe)

Re-dropped the `|μ_axis|` factor in `_ubld.assemble_inflow_axis` (the same
Mode-3 missing-factor as Branch 1) under canonical `python -O`:
`test_d2_exact_on_bilinear` FAILED (rc 1, via `pytest.fail` — fires under -O),
`test_divV_kernel_view_equals_dense` (d=1) stayed GREEN. The d=1 paths are
structurally blind to the multi-axis inflow factor (their RHS is built inline /
single-face); the d=2 exact-on-bilinear gate is the FIRST multi-axis consumer
and the load-bearing teeth. Same epistemic structure as Branch 1's ERR-060.
ERR-060 already logged (Branch-1); the Branch-2 d2 test carries `@catches`.

## Notes / decisions

- **Home:** `orpheus/sn/spatial/_ubld.py` (private, near `linear_discontinuous.py`
  + `scheme.py`/`scan.py`/`diamond.py` — the per-cell discretization layer).
  NOT `orpheus/derivations/` (that's symbolic references). NOT exported from
  `orpheus/sn/spatial/__init__.py` (private primitive; the scheme imports it
  directly; S2 will decide on any public surface).
- **structural independence (above the trusted-library line):** Branch 2 is
  hand-written numpy; Branch 1 is `lambdify`/symbolic SymPy. They share ONLY
  numpy/sympy primitives (below the line). The primitive unit tests cross-check
  numpy assembly against the symbolic matrices — a SEPARATE structural angle.
- **NOT wired d≥2:** the `_kernel_terms` `len(s_axes)!=1` raise stays IN PLACE
  (S2 closes it). `θ=1/3`, the slope-row sign convention (`_LDCellTerms.slope`),
  and every external API/signature are UNCHANGED.
- **NOT touched:** `loss_representation.py`, `sweep_graph.py`, `scattering`, the
  iterate, `scheme.py` (the base reconstruction staticmethods are CONSUMED, not
  modified). NO git add. The pre-existing uncommitted `error_catalog.md`
  (Branch-1's ERR-060) + auto-regen `matrix.rst` (foundation +10) left as-is.

## Owed (downstream)

- **S2** — wire the d≥2 bilinear cell-batch kernel onto the dense primitive +
  the `2^{d-1}`-moment face-cochain payload widening (literature memo §3); close
  the `_kernel_terms` `len(s_axes)!=1` raise. The primitive is built and
  d=2-exact-on-bilinear-verified — S2 is wiring + face-payload, not new algebra.
- **Archivist** — DISPATCH_REQUEST emitted for the rich narrative of the
  `ld-ubld-branch2-primitive` stub (followup:false).

## LESSON (proposed skill note)

Single-sourcing a Rule-of-Three across DIFFERENT reduction conventions (×V /
÷V / ×V-scan) is a Pattern-7 (normalise-at-definition) win, but the convention
crosswalk is load-bearing: the ×V `D₂' = θ·V·d2` (extra θ) vs the ÷V `d2`
scaling is exactly the bug habitat. Write the ×V↔÷V↔scan scaling table FIRST
(which quantity carries the θ, which carries the V), then the helper. The
`schur_xV` link test (vs the dense primitive) caught my first-draft `V·d2`
(missing θ) at write time — the d=1 link proof is not ceremony, it is the teeth.
