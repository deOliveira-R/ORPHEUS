---
name: sotp-separability-verification
description: Durable architectural-verification lesson from Wave T (LANDED). The SOTP master condition — a sum-of-tensor-products separable operator requires a Cartesian-product per-axis decomposition; coupled physics (per-material XS, per-ordinate M-M curvilinear recurrence) is NOT SOTP-able and MUST fall back to OperatorSum. How to gate a "separable" claim and which drift bound applies.
metadata:
  type: feedback
---

Wave T (T.3 scattering, T.4 streaming) carved SN operators toward a
sum-of-tensor-products (SOTP) separable form. LANDED on `main`
(`AngularRedistributionOperator`, `M_spatial`/`M_angular_redist` in
`orpheus/sn/operator.py`). The durable test-design lesson:

**THE SOTP MASTER CONDITION.** An operator is SOTP-able (expressible as
`Σ_k A_k ⊗ B_k ⊗ ...` per-axis) ONLY when it admits a CARTESIAN-PRODUCT per-axis
decomposition — each factor acts on ONE axis independently. Coupled physics
breaks this:
- **Per-material cross-sections** couple the spatial axis to the energy/group
  axis (Σ varies per cell AND per group) → not a clean tensor product.
- **The Morel-Montry curvilinear sequential half-grid recurrence**
  (per-ordinate angular redistribution) couples the angular axis to itself
  across ordinates → NOT SOTP-able. It MUST fall back to a bespoke
  `OperatorSum` leaf (`AngularRedistributionOperator`), not a separable factor.

When verifying a "this operator is separable / SOTP" claim, the gate is an
`assert_separable` foundation test that FAILS on a heterogeneous OperatorSum —
NEVER assert separability on a config where per-material XS or per-ordinate
coupling is present; the claim is FALSE there and the test must catch it.

**ROUTE A vs ROUTE B — which bit-identity bound applies.** Two implementation
routes for a separable rewrite: Route A (composition-shaped, preserves the FP
reduction tree → strict `array_equal` bit-id) and Route B (sum-of-products,
re-associates the reduction → `nulp` relaxation). Choose per the vv-principles
3-criteria gate. T.4's expanded drift bound for the streaming rewrite was
`iter_count × nx × ULP ≈ 4e-13` (ABOVE the naive 1e-13 but BELOW 1e-12) →
Route A default when the bound is tight enough to stay sub-conv-tol. The
scattering T.3 default was Route B (`nulp = 4·(L+1)`, the per-Legendre-order
einsum budget). The principle: a separable rewrite that re-associates a
reduction is principled-equivalence (named per-axis intermediates), not bit-id
— gate it at nulp(reduction_depth), anchor it at k_∞ + MMS, NOT at old-vs-new
ULP alone.

**The slab is the degenerate SOTP case (structural-preservation).** In slab
geometry the angular redistribution is a ZeroOperator (no curvilinear
α-redistribution) → the operator IS trivially separable there. So a slab-only
separability test PROVES NOTHING about the coupled curvilinear case (vv §H7 /
ERR-006 family — slab hides redistribution bugs). The curvilinear companion is
MANDATORY, exactly as for any redistribution-bearing claim.

See [[regression-tolerance-design]] (Route A array_equal vs Route B nulp), the
`coding-elegance` skill (Pattern 5 build-the-primitive; the OperatorSum
fallback for non-separable physics).
