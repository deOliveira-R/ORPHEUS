---
name: phase5a-5c-angular-windowing-lessons
description: Durable verification lessons from the SN 2-D angular-windowing carves (5a representation-only → moment-space SI iterate; 5c output-dropping → in-sweep moment accumulation). LANDED, tests live. THE TRAP (existing 2-D snapshots are all P0-isotropic), the moment-tensor fuller-view oracle that catches ℓ/m drift, nulp(n_diag) relaxation rationale.
metadata:
  type: feedback
---

The angular-windowing carves are LANDED on `main`. 5a: the 2-D Cartesian SI
iterate is `HarmonicMomentField` moments not full `AngularFlux`
(`_MomentWindowedResolvent`, gated `sn_mesh.reduced is None`). 5c: the
per-anti-diagonal sweep OUTPUT is dropped — moments accumulate in-sweep inside
`apply_windowed` (`einsum("nlm,ngd,n->lmgd", Y_oct, ψ_avg, w_oct)`), the
post-sweep flat `MomentProjection.apply` kept ONLY as a verification oracle.
Tests live at `tests/sn/solve/test_2d_anisotropic_windowing.py`. This note
keeps the WHY.

**1. THE TRAP (the canonical Cardinal-1 / Mode-7 instance for this codebase).**
EVERY existing 2-D snapshot was P0 ISOTROPIC (`2d_1g_LS4_dd_15x15`,
`2d_2g_LS4_dd_8x4_het_si`); the only P1-anisotropic snapshots were 1-D
(slab/sphere); the `2d_octant_05_qaniso` case is anisotropic-SOURCE (it
BYPASSES the moment projection). **Consequence: a windowing that silently
drops φ_ℓ≥1 PASSES THE ENTIRE EXISTING 2-D SUITE.** The fix was a NEW snapshot
`2d_2g_p1_aniso_dd_8x4_het_si` (fuel|moderator 8×4, B-mix μ̄=0.6,
vacuum-x/reflective-y, scattering_order=1, SI) + a NEW live anisotropic
windowing file. MUTATION-CONFIRMED: a φ₀-only mutant (P0-SI vs P1-Krylov-ref)
diverges 22% → fails allclose by 5 orders. **General rule: before carving any
moment-reduction path, audit whether the existing snapshots EXERCISE the
moments being reduced. If they are all isotropic, the carve cannot be gated by
them — manufacture an anisotropic case FIRST.**

**2. The Krylov path is a GENUINE cross-check, not a twin.** 5a windows ONLY
the SI iterate; Krylov stays full-angular. So `SI ≡ Krylov` (rtol 1e-6) is a
real structural cross-check for ℓ=0 (different angular representation, same
fixed point) — NOT self-consistency. Keep the cross-check leg pointed at the
UN-windowed path.

**3. The fuller-view moment-tensor oracle catches ℓ/m-slot drift (Mode-5) that
the scalar cross-check is BLIND to.** 5c's `_MomentWindowedResolvent.solve`
collapsed to forward `base.solve_moments(...)`; the full flat
`MomentProjection.apply` was kept as the oracle, pinned by
`test_2d_windowed_moments_in_sweep_equal_post_projection` — full moment tensor
SUT-vs-oracle, INCLUDING the ℓ≥1 block. This is load-bearing because the
scalar-flux cross-check (G1(a) SI≡Krylov, G1(b) Σwψ=φ) sees ONLY ℓ=0 — a
ℓ/m-slot swap or a dropped ℓ≥1 moment is INVISIBLE to it. The oracle leg is
PROCEDURALLY-not-structurally-independent (shares Y / weights / kernel) so it
MUST lean on the scalar cross-check + k_inf for the independent ground for ℓ=0,
and the new test pins the ℓ≥1 block those grounds can't see. PAIR it with an
ℓ≥1-non-vacuity precondition (`max|m[1]|/max|m[0]|>1e-3`) so it cannot pass as
a latent dud on an accidentally-isotropic config.

**4. nulp(n_diag) is the PRINCIPLED relaxation when per-anti-diagonal `+=`
reorders the sum.** Cross-octant accumulation reorders the N-ordinate sum vs a
flat `MomentProjection.apply` reduce → bit-id breaks at FP-non-associativity.
Accept per the 3-criteria: (1) per-octant partial moment IS the named principled
intermediate φ_ℓ^m; (2) anchored by closed-form k_inf + the SI≡Krylov-FULL
ground; (3) drift bounded `reduction_depth = N·eps` (measured 2.74e-16 / 4 ULP
≤ 4·N·eps). K=32 = next-pow-2 above the N=24 reduction depth. The existing DD
regression snapshot ALREADY drifts ~6920 ULP at its `SAFETY×conv_tol≈1e-11`
gate (5b inheritance) → ~10× headroom → 5c needed NO new relaxation there, and
strict `-W error::DriftWarning` is already moot for that case (do NOT
regenerate the snapshot, do NOT gate 5c on strict bit-id of the het-aniso
case). 5a representation-only PREFERS same-site reduce → bit-id `array_equal`
holds; only the OUTPUT-dropping 5c reorders.

**5. Anti-scope: keep the one-shot reconstructions full-angular.** Both
final-recon paths (eigen plain `transport_sweep`, fixed-source
`base_resolvent.solve`) stay full sweeps so user-facing `Solution.angular_flux`
is bit-id → the recon pin STAYS `np.array_equal`. 1-D / curvilinear / Krylov
stay full (gate `sn_mesh.reduced is None`). A windowing carve is 2-D-SI-ONLY;
any drift in 1-D/sphere/cyl/matvec/Krylov is a BUG, gated by the all-geometry
DD regression `test_dd_regression`.

**6. MEASURED win = 3.06× peak memory** (S8/4g/24×24: 2.26→0.74 MB) — the
transient full-angular per-sweep machinery is what 5c eliminates; the win grows
with angular order N (transient ∝ N, moment tensor fixed (L+1)(2L+1)). This is
why 5a alone gave 18× smaller PERSISTENT iterate but only ~1.2× PEAK (the
per-sweep transient dominated until 5c dropped it).

See [[wavefront-flux-carve-lessons]] (the orthogonal face-window carve, which
keeps outputs bit-id where windowing drops them), [[regression-tolerance-design]]
(nulp vs SAFETY×conv_tol), [[eigen-on-nonfissile-mixture-is-malformed]] (the
mixture-B fixed-source reformulation used to build the aniso snapshot without
fission).
