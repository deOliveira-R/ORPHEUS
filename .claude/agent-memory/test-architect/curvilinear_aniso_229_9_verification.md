---
name: curvilinear-aniso-229-9-verification
description: Verification spec for #229 (path-I angular-floor retune) + #9 (path-II P1 Legendre scattering, NEW curvilinear coverage). The TWO-unrelated-anisotropic-paths split; sphere has a pre-floor O(h²) window but CYLINDER DOES NOT at any practical quadrature (floor dominates — abandon the cyl spatial claim); the L0 operator-admits P1-source-without-MMS-source trick; the P1-lowers-keff leakage-monotone direction.
metadata:
  type: project
---

Pre-implementation spec for closing #229 + #9 (anisotropic scattering in
curvilinear 1-D sphere/cylinder SN). Ladders RUN on `main` 2026-06-13. The
committed tests will be the live record; this keeps the WHY + the measured
windows so a fresh session doesn't re-run the ladders.

**THE LOAD-BEARING FRAMING (two unrelated "anisotropic" paths).**
- **Path (I)** = geometric angular redistribution `(1−μ²)/r ∂_μψ` (sphere) /
  `ξ²B/r` (cyl), threaded by the M-M Carlson α-dome. **P0-only**; the
  "anisotropy" is in the angular-flux ANSATZ. The existing curvilinear aniso
  MMS cases (`SN{Spherical,Cylindrical}AnisotropicMMSCase`, ansatz
  `(A+ζB)/W`) exercise ONLY this. #229 is a path-(I) test-design floor.
- **Path (II)** = P1+ Legendre SCATTERING moments `R·Λ·M` (`scattering.py
  build_aniso_source`, scattering_order≥1), geometry-AGNOSTIC, wired
  identically for all geometries via `_within_group_triple` (the `S` in the
  `(L+C),S,B` triple carries P1 when scattering_order=1). NO curvilinear test
  exercised path (II) — #9 is NEW coverage. `sigma_s1` is "reserved, not
  wired" on the aniso MMS cases — KEEP it that way (see #9 L0 below).

**Phase 0 prereq RED** = pure 4-D→3-D test-index migration. `external_source`
returns `(N,ng=1,nx)` (rank-d carve `6ae3da8`); the two `@foundation`
symbolic guards still index `Q_numerical[:,0,:,0]`. Fix → `[:,0,:]` (lines
259+313 of `test_sn_mms_anisotropic_symbolic.py`). No deeper contract issue
(producer correct; `Q_sympy_grid/sum_w` already `(N,nx)`).

---

## #229 — path-(I) angular-floor retune (Option-2 two-claim split)

⭐ **THE KEY ASYMMETRY: sphere has a pre-floor O(h²) window, CYLINDER DOES
NOT — at ANY practical quadrature.** Measured (SI inner):
- SPHERE S16 nx{10,20,40}: orders 1.9794, 1.9780 (deterministic, ~4% margin
  > 1.9); floor(nx=160)=7.29e-4 sits 5× below the segment's finest (3.82e-3).
  Assert `min(orders[:2]) > 1.9` on the COARSE segment — NOT `orders[-2:]`
  (the issue's current code includes the floor-degraded 1.715/0.675).
  S32 extends the clean window to nx=80; floor drops to 2.89e-4 (2.52× per
  Sxx doubling).
- CYLINDER: even n_mu=16 (N=512) reaches only 1.803 on the COARSEST `5,10,20`
  step, then 1.195. The angular floor DOMINATES the spatial error before
  O(h²) establishes. **ABANDON the cyl spatial-O(h²) claim** (do NOT assert
  a cyl spatial rate gate). Runtime is NOT the blocker (n_mu=16×nx=40=2.5s);
  the MATHEMATICS is — the floor is structural. vv anti-pattern #5/#17: a
  claim that cannot hold MUST NOT be asserted; pin what IS true instead.

**Option-2 second claim = VERIFIED FLOOR (angular convergence at fixed-fine
nx).** Make the floor a verified quantity by asserting it SCALES:
- sphere nx=160: `err(S32) < err(S16)/2.0` (measured 2.52×; gate 2.0× margin).
- cyl nx=80: `err(n_mu16) < err(n_mu8)/2.0` (measured 2.38×; gate 2.0×).
Cyl floor at fixed nx=80: n_mu 4→8→16 = 1.90e-2, 7.39e-3, 3.11e-3 (2.57/2.38×).
Sphere floor at fixed nx=160: S4→S8→S16→S32 = 4.69e-3,2.82e-3,7.29e-4,2.89e-4.

**ERR-058-class re-floor catcher (vv Mode 7, question d).** Keep one aniso
curvilinear case with a WIDE absolute band: cyl `psi_qext_floor_band`
`1e-3 < err < 5e-2` (catches a re-floored ~5e-2 wrong-FP) + `@catches
("ERR-026")`. Sphere `psi_qext_band` loosen the current `1e-3` upper to
`5e-3` (the S16 floor at nx=40 is 3.82e-3 — the existing 1e-3 is what the
xfail rides on; 5e-3 passes the legit floor, catches the wrong-FP class).

**ALL 5 xfails REMOVED.** 6 equation labels (`sn-mms-{sph,cyl}-aniso-
spatial-convergence` + 4 `{sph,cyl}-aniso-{psi,qext}`) all migrate to GREEN
tests (question e). The cyl spatial-convergence LABEL re-attaches to the
floor-scaling test (semantics shift "spatial O(h²)"→"angular floor scales" —
ARCHIVIST must update `discrete_ordinates.rst:3355` narrative). The 5th
marker `test_mms_prescribed_inflow_sphere_activates_redistribution` rides the
SAME floor → split its band to `1e-8<err<5e-3`, drop strict-xfail; its
converged-VALUE `assert_allclose(phi,ref,2e-2)` STAYS (the dropped-q.boundary
catcher, floor-independent).

---

## #9 — path-(II) P1 Legendre scattering in curvilinear (NEW)

⭐ **L0 trick: operator-admits with a KNOWN aniso angular-flux input —
NO σ_s1-wired MMS source needed (avoids the costly symbolic P1-source
derivation).** Feed ψ_ref,n=(A+ζB)/W to the within-group `S` at
scattering_order=1, isolate P1 as `S₁.apply(ψ)−S₀.apply(ψ)`, assert PER-
ORDINATE (NOT weight-summed — α-dome telescopes, anti-pattern #8) vs a
structurally-independent hand-ref:
- SPHERE (fully SH-table-INDEPENDENT, strongest): `q_n^P1 =
  (1/W)·3·μ_n·Σ_s1·φ_1`, `φ_1 = B(r)·Σ(w·μ²)/W`. **VALIDATED rel 3.9e-15.**
- CYL (explicit SH moment-sum, indep of `R·Λ·M` einsum): `q_n^P1 =
  (1/W)·3·Σ_s1·Σ_m Y₁ᵐ(Ω_n)φ₁ᵐ`, `φ₁ᵐ=Σ_n w_n Y₁ᵐψ_n`. **VALIDATED 6.2e-15.**
1g is LEGITIMATE here — this is a flux-shape/OPERATOR claim (the per-ord P1
source reads φ_1, flux-shape-dependent by construction), NOT an eigenvalue
claim. Cardinal Rule bars 1G EIGENVALUE only. Add explicit `_require(
max|S₁−S₀|>1e-6)` negative control (a dropped-P1 makes S₁−S₀≡0 → fails the
non-zero hand-ref match anyway, but pin it; sphere max=3.7e-3, cyl 5.6e-4).

⭐ **L1 directional EIGENVALUE: P1 forward-peaked (μ̄>0) LOWERS k_eff vs P0.**
Physics: positive μ̄ preserves forward direction → in a finite vacuum-bounded
sphere, forward-preserved scattered neutrons more likely cross the outer
boundary → ENHANCED LEAKAGE → lower k_eff. VALIDATED robust:
- homog sphere R=4/10/25: Δ=P1−P0 = −3.76e-3 / −1.32e-3 / −2.88e-4 (sign
  always neg; |Δ| GROWS as sphere shrinks = leakage-monotone signature = the
  structural negative control a sign-flipped/absorption-mimicking P1 violates).
- HET fuel-core(r<5)+mod-shell R=10: Δ=−1.40e-2 (140× the 1e-3 bar).
Materials: fuel=`get_mixture("A","2g")` (ONLY fissile 2g; SigF=[.01,.08],
asymmetric downscatter-only P0 → avoids ERR-002+1G degeneracy; P1
SigS[1]=[[.019,.005],[0,.045]] forward-peaked), mod=`get_mixture("C","2g")`.
Two L1 rows: (1) het-sphere `keff_P1<keff_P0` AND `1e-3<(P0−P1)<5e-2`;
(2) leakage-monotone `(P0−P1)|R=4 > (P0−P1)|R=25 > 0` (the mechanism pin).
Mirrors 2-D `test_keff_2d.py:256 test_p1_changes_heterogeneous_keff` but
asserts SIGN not just presence. Reachable via public `solve_sn(...,
scattering_order=L)`.

**L2 SUBSUMED/DEFERRED** — L0 (source threads at machine precision) + L1
(het+2G+curv+aniso eigenvalue) already cover it; a P1-convergence L2 needs
the σ_s1-MMS-source (costly) AND rides the SAME #229 floor → marginal.
Cheapest future L2 = self-convergence (P1 het sphere keff vs mesh), weaker
than L0+L1. Defer with a one-line #9 closeout note, don't build now.

**#9 LANDED 2026-06-13** (branch `fix/curvilinear-aniso-pole-and-clamp`,
parallel to W1 `b2d8a6d`). Files: L0 =
`tests/sn/verification/mms/test_curvilinear_aniso_scattering_p1.py` (2 `@l0`
+ `verifies("pn-scatter","flux-moments")` — sphere SH-table-indep + cyl
explicit-Y₁ᵐ-sum, both per-ord vs hand-ref + `peak>1e-6` neg-control); L1 =
`TestSphereP1DirectionalEigenvalue` in
`tests/sn/eigenvalue/test_keff_curvilinear.py` (2 `@l1`+`verifies
("pn-scatter")` — het-vac-sphere sign+band, homog-vac leakage-monotone
R4>R25). ⭐ ONE SPEC CORRECTION: the L1 leakage rows REQUIRE `BC.vacuum`
outer — the `curvilinear_homogeneous_mesh`/`_two_region_mesh` helpers
DEFAULT to `BC.reflective`, and a reflective sphere has NO leakage → P0≡P1
(k=k_inf=1.875, Δ≈1e-10). Pass `bc=BC.vacuum` (helpers accept it; `from
orpheus.geometry import BC`). Measured GL8: HET-vac R10 Δ=1.40e-2; homog-vac
R4 Δ=3.75e-3 > R25 Δ=2.89e-4 > 0 — matches spec targets. L0 relerr sphere
4.7e-15 / cyl 5.6e-15 (machine precision, spec said 3.9e-15/6.2e-15 — same
order; my probe A(r)/B(r) differs from validation run, claim is
machine-precision agreement). `pn-scatter`/`flux-moments` are the REAL
geometry-agnostic labels (existed; prior tests only 2-D Cartesian — #9 is
first curvilinear exercise). NO production code touched; path-II works.
GREEN under `python -O`. **L2 still deferred per above** (no σ_s1-MMS built;
self-convergence is the cheapest future row but weaker than L0+L1).

**xs_library**: `from orpheus.derivations.common.xs_library import
get_mixture`; `m.SigS` = len-`nLeg` list of SCIPY SPARSE (`.todense()`),
SigS[0]=P0, SigS[1]=P1; `m.SigT`, `m.SigF`, `m.chi`. Only A is fissile in 2g.

See [[lessons]] (L1 homog-blind-to-curvilinear), [[phase4-46-nonvacuum-mms-
ansatz]] (the prescribed-inflow MMS this builds beside),
[[regression-tolerance-design]].
