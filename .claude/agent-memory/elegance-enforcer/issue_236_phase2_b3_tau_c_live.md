---
name: issue-236-phase2-b3-tau-c-live
description: "#236 Ph2 B3 (the LIVE τ fold) elegance verdict — PASS-WITH-NITS, no must-fix, uncommitted on feature/sn-spatial-angular-product"
metadata:
  type: project
---

# #236 Phase 2 Step B3 — "the live τ fold" — PASS-WITH-NITS (no blocking nits)

**Branch** `feature/sn-spatial-angular-product` (MAIN checkout, NOT a worktree), UNCOMMITTED at review.
This is the LIVE-consumer step after Steps A/B1/B2 (closure PRODUCES τ + c, stamps `CellVisit.c_in/c_out`).
B3 makes the live sweep + scan + matvec CONSUME the closure-owned τ/c so geometry-τ goes dead → Step C deletes it.

**6 production seams** (all PASS): `pole_angular_closure.py` (`tau_per_ordinate` accessor+cache; `_build_c_per_ordinate_cache`→`_build_per_ordinate_cache` now gathers 3 constants) · `scheme.py` (`CellVisit.tau: float = 1.0`) · `geometry.py` (`_make_cell_visit` stamps τ) · `cell_balance.py` (`cell_balance_terms(..., *, c_in, c_out)` stops rebuilding c) · `diamond.py` (`DD.update` reads `visit.tau`/`visit.c_*`) · `sweep_cache.py` (`GeometryCoefficients` sources τ from `closure.tau_per_ordinate`, derives `tau_inv`/`mm_a_in_coeff`).

## The 4 load-bearing rulings (all PASS, grep-verified)
1. **τ SINGLE-SOURCED + Step-C precondition HOLDS.** c-formula now in EXACTLY ONE production site (`pole_angular_closure.py:1108-1109`, the closure ctor); `cell_balance.py` rebuild DELETED. Grep: every `.tau_mm` in the 3 consumers is COMMENT-ONLY; the only LIVE `.tau_mm` reads are `reduced_operator.py` (geometry producer Step C deletes) + the closure ctor's `reduced.tau_mm` (and even that is superseded — closure produces τ via `morel_montry_tau_per_level(quad,coord):1070`, NOT `reduced.tau_mm`). "Nothing live reads tau_mm" = MET.
2. **Primitive/product split RIGHT (Pattern 5).** Closure exposes PRIMITIVE τ only; `tau_inv=1/τ`/`mm_a_in_coeff=(1−τ)/τ` are consumer-local L16 perf-hoists (scan recurrence coeffs, `loss_representation.py:3270-3271`) — pushing them onto the closure would put a consumer's layout opt into the producer contract (wrong direction). The Q2 "DD raw-τ `(ψ̄−(1−τ)ψ_in)/τ` vs scan pre-split `tau_inv·ψ̄−mm_a_in_coeff·ψ_in`" = L21 ONE-operator-two-applications (verified algebraically identical); BOTH now read the single τ source → NOT a collapse target. Standing Pattern-2 twin, acceptable, hazard reduced not introduced.
3. **Pattern 4 + separation.** DD geometry-blind AND closure-blind (reads `visit.tau`/`visit.c_*` as DATA, never the closure type). `CellVisit.tau=1.0` default CORRECT (neutral M-M weight; 0.0 = ÷0 landmine in the recurrence; docstring states it). Unbound legacy sets `_tau_per_ordinate_cache=None`+accessor raises, base annotations `_tau_per_level: |None`/`_tau_per_ordinate_cache: |None` both declared — mirrors c-cache bit-for-bit.
4. **c-fold signature reads like domain.** `cell_balance_terms(..., *, c_in, c_out)` kw-only — c are INPUTS not derived; prevents the same-typed-float positional swap (anti-#13/#31). Improvement over old `st.tau_mm` reach-in.

## Q5 (deferred-debt) — NOT worsened, NO new dup
`cell_balance_for_streaming` (`cell_balance.py:120`, the assembly twin / `diamond.py:289-303` Phase-6 target) NOT in B3 diff (grep -c = 0). DD.residual path already sourced c from `visit.c_*` (B2); B3 leaves it. Residual path reads NO τ (apply direction needs only c-constants in denom/numer) → already consistent.

## NITS (all FOLLOW-UP, none block commit)
- **NIT-1 (Pattern-11 record):** theory `.. todo:: Archivist expansion needed.` (`discrete_ordinates.rst:1171`) — established file convention (5×), `todo_include_todos=True` so renders + build-safe, carries full B3 content inline + names Step C as trigger → anti-#11 EXCEPTION (tracked). Only follow-up: ensure a tracked issue/plan step owns the archivist expansion.
- **NIT-2 (standing tripwire, at threshold):** raw-τ recurrence test-reference now in ≥2 test bodies (`_c_surrogate.py` covers c; `test_diamond.py TestBitIdenticalCurvilinear` bodies cover the raw-τ recurrence the surrogate does NOT supply — legit independent, NOT collapse). N=2 = the unify trigger; on a THIRD raw-τ test ref, extract `tau_recurrence_reference(psi_avg,psi_in,tau)` beside `c_from_streaming_terms`.
- **NIT-3 (record):** "0-ULP / verified-in-process |visit.tau−st.tau_mm|≡0" asserted in docstring but the probe is described not checked-in as a named test (your recurring "verified-in-process no artifact" tell). Empirically sound — 141 touched-module tests pass (DriftWarning-escalating snapshots + surrogate round-trip). Cosmetic: docstring should cite the gating snapshots rather than imply a dedicated identity assert.

## Gate ran
`pytest tests/sn/sweep/core/{test_diamond,test_cell_balance_for_streaming,test_ordinate_scan,test_sweep_cache}.py` → 141 passed, 1 skipped. Bit-identity claim holds operationally.

## Series continuity
B3 = the LIVE step (B2 was zero-live-risk `DD.residual`). After B3 the geometry-τ producer is dead → **Step C** = retire `spherical_streaming`/`cylindrical_streaming` τ + `StreamingTerms.tau_mm`/`alpha_*`. #248 owns the orphaned-Protocol question (B3 mints `tau_per_ordinate` on BOTH Protocol+ABC, no regression). Exemplary phased-carve consumer step: predecessor inline-c REMOVED not duplicated.
