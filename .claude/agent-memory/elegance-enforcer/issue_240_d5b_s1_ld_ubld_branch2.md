---
name: issue-240-d5b-s1-ld-ubld-branch2
description: D5b-S1 Branch-2 PRODUCTION numpy UBLD primitive + Rule-of-Three collapse of the 3 slab-LD views onto one d1_closed_form helper — PASS-WITH-NITS
metadata:
  type: project
---

# #240 D5b-S1 Branch 2 — production numpy UBLD + Rule-of-Three collapse (PASS-WITH-NITS)

Sibling of [[issue-240-d5b-s1-ld-ubld-branch1]] (the SymPy algebra-of-record). Branch 2 is the
PRODUCTION numpy counterpart: NEW `orpheus/sn/spatial/_ubld.py` (d-generic dense Kronecker primitive
`assemble_ubld`/`assemble_inflow_axis`/`per_cell_solve` + the vectorized `d1_closed_form`/`D1ClosedForm`
fast path) and the MODIFIED `linear_discontinuous.py` collapsing the THREE slab-LD views
(`_schur_terms` ×V, `_kernel_terms` ÷V, `affine_scan_coefficients` ×V-scan) onto the ONE helper.
Branch `feature/sn-space-angle-tier2`, MAIN checkout. Date 2026-06-16. Verdict **PASS-WITH-NITS**.

## The 4 load-bearing axes ALL PASS PROVABLY (not by inspection — by SymPy)

1. **Rule-of-Three genuinely collapsed (Pattern 2 / Cardinal Rule 2).** All three views derive
   coefficients from `d1_closed_form`; only the ×V/÷V/×V-scan SCALING is applied at the call site.
   I verified ALL SEVEN coefficient crosswalks reduce symbolically to 0 (`sp.simplify(diff)==0`):
   D₂', Schur S, eff_source, eff_numer, scan `a`, kernel_rhs, `w`. So the ONLY difference from the
   pre-carve inline code is FP-association (the documented ~1-ULP re-baseline). No view re-spells the
   2×2 independently — grep-confirmed each view body is now ~3 lines around `cf = d1_closed_form(...)`.

2. **d=1 STAYS on the fast path (L16).** Grep-clean: `linear_discontinuous.py` imports ONLY
   `d1_closed_form` (line 175); `per_cell_solve`/`np.linalg.solve`/`assemble_ubld` appear NOWHERE in
   production LD. The dense Kronecker path is the d-generic REFERENCE only (test + S2 d≥2). Confirmed
   no per-cell dense solve leaked into the d=1 production path.

3. **The primitive reads like the math.** `assemble_ubld` mirrors the symbolic Kronecker build
   one-for-one (factors `_GRAD_1D`/`_FOUT_1D`/`_FIN_TRACE`/`mass_1d`; `_batched_kron` = batched
   `sympy.kronecker_product`); the active-axis-grad / transverse-mass factorisation is identical to
   Branch 1. `D1ClosedForm` carries the named scale-free invariants `k`,`w`,`g`,`d2`,`eff_denom`
   (Pattern 3). Dead-weight audit: every dataclass field is consumed (`g`/`g_over_theta`/`d2`/`k`/
   `theta` internally, `w`/`eff_denom` internally AND externally). Zero unused fields.

4. **The ×V↔÷V crosswalk (D₂' = θ·V·d2, the extra-θ bug habitat) cleanly single-sourced.** The extra
   θ that the ×V slope denom carries (and the ÷V `d2` does NOT) lives as ONE named intermediate
   `d2p = self.theta * V * self.d2` in `schur_xV` (`_ubld.py:343`), proven `==` legacy `Σ_t·θ·h+|μ|`.
   The closeout records a first-draft `V·d2` (missing θ) caught at write time by the `schur_xV` link
   test vs the dense primitive — the d=1 link proof is the teeth, not ceremony. Matches the team's
   Pattern-7-crosswalk-table-FIRST discipline.

## Gates RE-RAN (all green)
- New gate `test_ld_ubld_primitive.py` 10 passed (-O). Includes the d=2-exact-on-bilinear ERR-060
  catcher (xy coupling) + the 3-view closed-form==dense reduction (the elegance Branch-1 CONCERN
  closed in code) + the LINK proof (LIVE production update/kernel/scan == dense d=1).
- `tests/sn/spatial tests/sn/sweep/core` 506 passed / 1 skip / 4 xfail (the bit-id DD negative
  control held — strict gate carries NO LD golden, its LD items are structural-only).
- LD two-paths scan≡DAG oracle + MMS-LD 5 passed / 1 xfail. Single-sourcing holds end-to-end.

## NITS (none commit-blocking; flagged for do-now where cheap)
- **Test :284** `np.all(...)` returns `numpy.bool_` → `_require(condition: bool)` — genuine-but-trivial
  typing nit (runtime fine). Cheap fix `bool(np.all(...))`. Same class as :111/:112 `.subs`-on-list.
- **Pyright :175 `_ubld` import-unresolved + :263 `key=` registry kwarg + :527 `A_total` not accessed**
  = ALL benign. :175/:263 are #226 D3-rename pyright-config staleness ([[reference-pyright-lsp-rooting]]);
  runtime import works (tests pass). :527 `A_total` is the established CONTRACT-UNIFORMITY-not-smell:
  DD's `affine_scan_coefficients` consumes it (`diamond.py:587 a=2|μ|·A_total/denom−1`); LD leaf
  doesn't need it but the polymorphic protocol (scheme.py:905) requires the uniform signature. Same
  ruling as [[issue-158-ld-spatial-occupant]].
- **`assemble_inflow_axis` `axis∈{0,d-1}` else-raise** is NOT anti-pattern-#7 boundary-special-case —
  it is an honest-per-phase capability boundary mirroring Branch 1 exactly (d≤2 needs only {0,d-1};
  interior-axis interleave deferred to S2 per Pattern 6). Fails loudly, doesn't mis-assemble. PASS.

## Carry-forward for S2 (when the d≥2 bilinear kernel wires onto the primitive)
- `_kernel_terms` `len(s_axes)!=1` raise stays IN PLACE (S2 closes it). The d2 inflow factor is the
  Mode-3 missing-|μ_axis| habitat both Branch 1 (ERR-060) and Branch 2's re-mutation caught — the
  d=1 paths are STRUCTURALLY BLIND to it; the d=2-exact-on-bilinear gate is the first multi-axis
  consumer with teeth. Any S2 face-cochain widening (2^{d-1}-moment payload) must re-probe it.
