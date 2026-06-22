---
name: issue-240-d5b-s1-ld-ubld-branch1
description: D5b-S1 Branch-1 SymPy d-generic UBLD reference + foundation gate — PASS-with-nits; the canonical algebra-of-record production reduces TO (not a copy), Kronecker SSOT, L11 thinness, pyright float|None narrows.
metadata:
  type: project
---

# #240/#38/#37 D5b-S1 Branch 1 — d-generic UBLD SymPy reference (PASS-WITH-NITS)

Branch `feature/sn-space-angle-tier2`, MAIN checkout. Files: `orpheus/derivations/discrete/sn/ld_ubld.py`
(SymPy d-generic Kronecker assembler, SHIPS NOTHING) + `tests/sn/spatial/test_ld_ubld_symbolic.py`
(6 `@foundation` oracles, all green; CLI 5 derivations all symbolic-zero). The multi-D analog of LD on a
Cartesian cell = BILINEAR/TRILINEAR DG-P1 (UBLD), 2^d moments, NOT simplex-P1 (Adams-2001: simplex FAILS
thick diffusion limit on quads, bilinear PASSES; the `xy` cross moment is the load-bearing term).

## The two load-bearing axes — BOTH PASS provably

1. **Kronecker SSOT (Pattern 5+2)**: 4 one-D factors (`mass_1d`/`grad_1d`/`fout_1d`/`fin_trace_weight`,
   :131-171) defined ONCE; `assemble_ubld` (:179-230) composes d-generically (`for a in range(d)` /
   `kronecker_product`) with ZERO per-d special-casing — d=1/2/3 run the SAME assembly. `A=G+F_out+Σ_t·M`
   (:229) is character-for-character MRM-2016 Eq.12. `xy` coupling EMERGES from `M₁⊗M₂`; no 4×4/8×8 entry
   hand-transcribed. General-first (anti-pattern-15 AVOIDED): d-generic ships first, production = rank-1 collapse.

2. **Canonical-source-NOT-copy (anti-pattern-1)**: symbolic module imports ZERO production (`import sympy`
   only). Production `_LDCellTerms`/`_kernel_terms`/`affine_scan_coefficients` proven to REDUCE TO the
   symbolic via the d=1 Kronecker-with-one-factor identity (`diff_A==zeros` CLI-confirmed). Correct
   relationship: symbolic=general algebra-of-record, production=specialization. A copy would be the smell;
   this is the inverse.

## NITS (all cheap, none architectural; approval conditions)

- **L11 THINNESS (CONCERN)**: the 3 SYMBOLIC d=1 oracles (`derive_d1_reduction`/`_kernel_view`/`_scan_view`)
  compare `assemble_ubld` against HAND-TRANSCRIBED SymPy copies of production formulas (comments cite
  prod LINE NUMBERS :342/:393/:435) — two re-typings agreeing, weaker than brief item-4 ideal. The GENUINE
  structural independence rides on ONE numeric test `test_d1_symbolic_primitive_matches_production_update`
  (:117-179, calls LIVE `LinearDiscontinuous().update`). Mitigation real but thin (Q̂=0 flat, 1 ordinate,
  1 cell). REMEDY (rec, not blocking): extend that numeric test to also tie ÷V `_kernel_terms` + ×V
  `affine_scan_coefficients` to live prod, so all 3 view-equals oracles each get a live witness.
  NOTE the d=2 bilinear oracle (:469-550) IS fully independent MMS (manufactures source/traces by
  independent SymPy integration of `ψ=a+bx+cy+dxy`) — and ERR-060 (`error_catalog.md:4469`, dropped
  |μ_axis| streaming factor) proves the d=2 oracle ALREADY caught a bug d=1 was blind to → L11 working.
- **prod LINE-NUMBER citations (CONCERN, Smell-7/anti-11)**: :342/:393/:435 cite `linear_discontinuous.py`
  line ranges (currently accurate: _schur_terms 332-335, _kernel_terms 443-453, affine_scan 564-571) →
  fragile, no integrity check. Fix → SYMBOL citations (Nexus/grep rename-tracked).
- **pyright ld_ubld.py (genuine, fix)**: :275/:277 `Expr`-not-`Matrix` return lie (`mu_axis * sp.Matrix(kron)`
  types as Expr; lift the `sp.Matrix(...)` wrap outside the scalar product — the return contract is
  "2^d-vector", Branch-2 consumers index it). :337 `simplify` overload = cosmetic `# type: ignore`.
- **DEAD-WEIGHT param (genuine, trivial)**: `downstream_face_trace(..., axis=0)` :290 NEVER reads `axis`
  (body branches on `d` only) — "for future use" smell. Drop until d≥2 face lands (Branch 2). Brief named
  "290 unused axis" — pyright doesn't flag (it's a ruff ARG check) but it's real dead weight.
- **pyright test float|None (genuine, fix)**: :145/:146 `float(st.abs_mu/.volume)` + :177/:178
  `assert_allclose(res.outgoing_spatial_flux|None)` — pyright does NOT narrow through the `_require()`
  helper (not assert/TypeGuard). Mirror prod `linear_discontinuous.py:321` direct-narrow; keep -O-safe
  `_require` for the runtime fire.
- **pyright test :138 = NOISE (do NOT fix)**: `slab_streaming(mesh, quad)` Quadrature-vs-AngularMeasure =
  Protocol invariance (`mu_x` mutable). Codebase-wide (test_diamond/test_linear_discontinuous + 9 files);
  `slab_streaming` documented to take "any AngularMeasure-shaped object". Same family as #226 config noise.
  The `orpheus.sn.spatial.scheme` import-unresolved (:50/:51) = known #226 D3-rename config staleness, also NOISE.

## Exemplary
-O-safe guards (`_require`/`_require_zero`/`_require_zero_matrix`, :54-77, function-call not bare-assert);
no `verifies(...)` correctly justified (algebra invariants = foundation level not L0/L1); `catches(ERR-060)`
grounded; one test per derive_* claim; honest scope-bounding (`assemble_inflow_axis` raises on interior-axis
d≥3 interleave = Pattern-6 deferral named, not special-case smell).

## Branch-2 forward watch (no issue needed, in closeout/plan)
Multi-D face cochain must carry 2^{d-1}-moment object (memo §3); affine-scannable likely does NOT survive
multi-D (memo §7.2). `downstream_face_trace`/`assemble_inflow_axis` are the widening seams — Branch-2
reviewer: ensure face-object widening single-sources through these two, no parallel face-assembly path.
