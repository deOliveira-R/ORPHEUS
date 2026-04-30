---
name: Issue #129 planar-limit resolution
description: Hollow cyl r_0→R does NOT reduce to slab at fixed thin L; the clean limit is Σ_t·L → ∞ (optically thick), and Cauchy's mean-chord theorem gives <chord>=2L identically.
type: project
---

Issue #129 (Phase G.4) is **resolved with redefined comparison**.

**The wrong claim** (original plan): hollow cyl `r_0 = 0.999 R, R = 1`
≈ slab `L = 0.001` at 1e-8. Empirical baseline showed 22.5% rel_diff.

**The correct picture**:
- Cauchy's mean-chord theorem proves `<chord>_slab(L) = <chord>_cyl-annulus(r_0, R) = 2L` EXACTLY at every R. Leading-order P_esc agrees identically.
- At fixed L=1e-3, varying R from 1 → 1000: k_eff gap PLATEAUS around 3.4% (does NOT vanish). This is the curvature-cleanest scan and it shows there is no clean limit at fixed thin L.
- At fixed R=1000, varying L: gap collapses to 1e-5 for L ≥ 0.1. Optical thickness drives kernel agreement, not curvature.
- The kernels disagree on the SECOND chord moment: cylinder has tangential chords ~√(2RL) that slab does not; slab has 1/|µ| tail regularised by E₃.

**Why**: At fixed thin L the structural mismatch persists — cylinder
`Ki_1` integrates over axial direction analytically, leaving an
in-plane chord distribution that supports tangential rays of length
~√(2RL). Slab `E_1` integrates over the in-plane y-z plane,
leaving the µ-direction with chord L/|µ|. These two distributions
share their first moment (Cauchy) but differ at second order.

**The proposed replacement test (G.4-redefined)**:
Slab L=1.0 vs hollow cyl r_0=999, R=1000, same XS. Σ_t·L=1
(optically marginal). Stage 1c shows rel_diff = 5e-6 already at N=3
dps=20.

**Why it matters for future work**: any future "limit comparison"
between curvilinear and Cartesian Peierls kernels MUST go through
the optically-thick regime, not the thin-shell regime. The
"geometric limit" intuition is broken by the Bickley axial
pre-integration — the kernels see different axial physics that
geometric thinning does not erase.

Diagnostics:
- `derivations/diagnostics/diag_issue129_planar_limit_stage1_rscan.py` (R-scan, L/R fixed)
- `derivations/diagnostics/diag_issue129_planar_limit_stage1b_fixed_L.py` (fixed L, vary R)
- `derivations/diagnostics/diag_issue129_planar_limit_stage1c_L_scan.py` (fixed R, vary L)
- `derivations/diagnostics/diag_issue129_planar_limit_stage1d_pesc.py` (Cauchy mean-chord identity)

Resolution memo: `.claude/plans/issue-129-planar-limit-resolution.md`
