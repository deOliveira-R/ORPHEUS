# DSA for SN (#2) — verification specification (PRE-implementation)

**Status:** design-only, PRE-carve. This spec SHAPES the build (3a/3b/3c);
it does not describe landed code. It refines the roadmap Phase-3 gate list
(`stencil_assembly_dsa_roadmap.md` §3c) into a full battery with configs,
tolerances, stabiliser analyses, and a mutation matrix.

**Prerequisites already landed:** `A_diff` (the in-algebra diffusion
operator, #290 @ `3a19133`) with `assemble()` / `as_matrix()` / the
`MatrixInverseOperator(FlattenedOperator(A, template))` exact-LU resolvent;
`evaluate_residual` (`sn/solver.py`) returning the typed
`FullField(interior=AngularResidual, boundary=AngularBoundaryResidual)`;
the `ScatteringOperator.foldable_part/residual_part` σ_r split; the
`SourceIteration` / `KrylovAcceleration` drivers with the ρ-honest STOP and
`residual_history` / `contraction_ratios` surfaced.

**Cross-dispatch dependencies (roadmap 3-P0, the MAIN agent owns these):**
`literature-researcher` supplies the consistent low-order stencil
coefficients (Larsen 1982 modified-4-step, Alcouffe 1977, Adams & Larsen
2002 §II for ρ=0.2247c) that **D1/D2** compare against;
`cross-domain-attacker` lands the R/P Petrov–Galerkin frame vocabulary
(moment-0 analysis face, isotropic reconstruction face) that **D7/D8** are
spelled in. Where a gate references a theory `:label:`, it is a
`TBD-but-shaped` placeholder the P8 doc pass fixes.

---

## 0. The organizing principle — the correction→0 partition (READ FIRST)

DSA's defining property is that its low-order **correction → 0 at
convergence** (the source of the correction is the scattering residual
`Σ_s(φ^{l+½} − φ^l)`, which vanishes at the fixed point). This is not a
detail — it is the axis the entire battery is organized on, because it
**partitions the failure modes into two disjoint classes with disjoint
catchers**:

| Bug lives in… | Effect on the converged FP | Effect on the rate | Caught by |
|---|---|---|---|
| **The within-group transport operator** `A = L+C−S−B` (the σ_r-fold #215; a sweep sign flip; a wrong closure fed to BOTH sweep and low-order) | **CHANGES the FP** (value-wrong) | — | **FP-invariance** D3–D5 + routing guard D10 |
| **The accelerator machinery** {R, P, `A_diff`, correction sign/scale, low-order D/removal, boundary rows} | **NONE** — correction→0 makes it value-safe *by construction* | degrades / diverges | **object gates** D1/D2/D7/D8 + **rate/stability** D11–D13 (+ D12 divergence) |

**Consequence that the plan-of-record MUST internalize:** the
FP-invariance gates are **structurally BLIND** to every R/P/`A_diff`/
correction-sign error — those leave the fixed point *identically* unchanged
(not sub-floor, *identically*: the correction is zero at the FP regardless
of how wrong `A_diff` is). A verification plan that gates DSA correctness
**only** on "converged flux ≡ SI fixed point" is Mode-9-complete for the
#215 class and Mode-9-**empty** for the entire accelerator-quality class. The
object gates (D1/D2/D7/D8) and the rate gates (D11–D13) are therefore
**load-bearing, not supplementary** — they are the *only* catchers for the
majority of the plausible-error surface. The mutation matrix (§4) makes this
concrete: of the eight canonical implementation errors, exactly **one** (the
σ_r-fold) reds the FP gates; the other seven are caught only by object/rate
gates.

The correction→0 gate (D6) is the *proof of the safety property itself* — and
it carries its own blind spot (a **dead R** that returns 0 makes δφ≡0, passing
"correction→0" trivially), closed by the non-trivial-first-iterate pairing.

---

## 1. Claim-layer / pillar table (per AGENT.md §1.5 — gated before any row is written)

| Gate | Claim layer | Pillar | Structurally-independent ground | 1G legit? |
|---|---|---|---|---|
| D1 SN-DD reduction ≡ stencil | term (L0) / object | closed-form (hand-posed P1 reduction) | analytical moment reduction on a 3-cell uniform slab, computed by hand | n/a (object) |
| D2 reduction ↔ `A_diff` consistency | derivation-consistency / object | — (two ORPHEUS assemblies) | independence supplied by **D1's** hand-posed anchor, NOT by D2 itself | n/a (object) |
| D3–D5 FP-invariance | flux-shape **invariance** (L2) | plain-SI FP (un-accelerated), correctness backing = existing SN L1 anchors | the plain-SI FP is the un-accelerated reference; its *correctness* rides the existing homogeneous k_inf / `Q/Σ_t` anchors | **NO — ≥2G** |
| D6 correction→0 | correctness-safety property (L1) | closed-form (the residual→0 identity) | the scattering residual vanishes at the FP by construction | ≥2G for the pairing |
| D7 R conservation | object / balance identity (foundation) | closed-form (discrete `⟨1,Rr⟩=⟨1,r⟩`) | hand-written moment-0 identity with explicit `w_n`, `V_i` | n/a (object) |
| D8 R/P adjoint-consistency | object / frame law (foundation) | closed-form (`R₀ ≡ M₀^adjoint` under the frame Gram) | the PG-frame `R.H` adjoint identity (frame.rst) | n/a (object) |
| D9 no-masking | flux-shape **co-variance** (L2) | plain-SI FP on a *mutated* operator (both must move together) | the deliberately-seeded transport bug | ≥2G |
| D10 routing guard | structural (foundation/sentinel) | — (AST / call-graph) | the σ_r-fold call signature | n/a |
| D11 ρ measurement | **rate** claim | closed-form / semi-analytical (ρ=0.2247c, Adams & Larsen 2002 §II) | the Fourier dispersion bound | **YES — rate is flux-shape-independent (declared)** |
| D12 reflective stability | rate + convergence (boolean) | closed-form (ρ<1 ⟹ converges) | plain-SI converges on the same slab (the spike control) | ≥2G for the value pairing |
| D13 c→1 × thickness table | rate characterization | semi-analytical (ρ=0.2247c band) | the Fourier bound | YES for ρ; ≥2G for the count table's value guard |

**MMS is absent by design** — DSA makes no eigenvalue claim (it is an
accelerator; the eigenvalue is the SN solver's, verified elsewhere), and MMS
cannot prove the invariance/rate claims that are DSA's actual content.

---

## 2. Config catalog (the fixtures — build them ONCE, reference by name)

Config-blindness (AGENT.md §0.6, lessons L1) is the dominant test-design
hazard here; each config below is chosen to BREAK a specific blindness. The
`make_mixture` / library-A/B/C/D trap (lessons L1: **`make_mixture` hardcodes
`SigL=0` and every A/B/C/D fixture has `Sig2=0`; there is no `sig_l`
parameter**) means every anisotropic (ℓ≥1) config MUST build the `Mixture`
constructor **directly** with a non-zero `SigS[1]` P1 moment — a `make_mixture`
call silently nulls the ℓ≥1 physics D3 exists to protect.

| Name | Geometry / dims | Groups | Materials | BC | Quadrature | Scatter | Purpose / breaks which blindness |
|---|---|---|---|---|---|---|---|
| **CFG-ANISO-VAC** | 1-D slab | ≥2G | het (fuel A + moderator B zones) | **vacuum** | S8 GL (non-uniform w) | **ℓ≥1 (P1/P3, `SigS[1]≠0` built directly)** | primary FP-invariance; breaks flat-flux + isotropic-snapshot + 1G blindness; the σ_r-fold is 56% wrong here |
| **CFG-ANISO-HET** | 1-D slab | ≥2G | strong contrast (absorber C + scatterer D) | vacuum | S8 GL | ℓ≥1 | redistribution/spatial-distribution stress; the σ_r-fold is 46% wrong here |
| **CFG-DIAG-2D** | 2-D Cartesian | ≥2G | het | vacuum | **level-symmetric / lebedev (NOT axis-aligned product)** | ℓ≥1 | angular-schedule interaction; breaks the ERR-056 shared-face blindness (vv Mode 9 (b)) |
| **CFG-HOMOG-INF** | 1-D slab, periodic/infinite | **1G legit (rate)** | homogeneous | reflective/periodic | S8 GL | isotropic | the classic ρ=0.2247c Fourier problem; rate claim ⟹ 1G legitimate (declared, AGENT.md §5) |
| **CFG-REFLECT** | 1-D slab **and** 2-D box | ≥2G | homogeneous + het variant | **fully reflective** | S8 GL (1-D), level-sym (2-D) | isotropic + ℓ≥1 variant | the historical divergence regime (the naive spike DIVERGED here); c up to 0.99 |
| **CFG-TINY** | 1-D slab, 3–4 cells | 2G | het (2 zones) | vacuum + reflective variants | S4 GL | ℓ≥1 | object-level matrix gates D1/D2/D7/D8 — small enough to hand-pose |
| **CFG-THICK** | 1-D slab | ≥2G | homogeneous, cell optical thickness Σ_t·Δx ∈ {0.1, 1, 10, 100} | vacuum + reflective | S8 GL | isotropic | the consistency stress axis (thick cells expose partial-consistency degradation) |

**The isotropic fully-reflective box is DELIBERATELY ABSENT as a
correctness config** — it is the designed-green degeneracy (vv Mode 9): the
σ_r-fold is *exact* on it, so a FP-invariance gate there is a false green.
It appears only as CFG-REFLECT's *stability* discriminator (D12), never as a
value gate.

---

## 3. The gates

Each gate: **purpose · config · measured quantity · tolerance+justification ·
stabiliser analysis · mutation rows · marker/level · runtime tier**. Every
mutation is validated under the canonical `python -O` (lessons L4 / vv Mode
8: bare `assert` in production/helper modules is stripped — gates use
`np.testing.assert_*` / `pytest.raises` / explicit `raise`; the D10 routing
guard uses an AST/counter tooth, never a bare assert).

---

### Phase 3a — consistency by derivation

#### D1 — the SN-DD moment reduction ≡ hand-posed P1 stencil (+ derivation-side sanity pins)

- **Purpose.** The consistent low-order operator is the P1/moment reduction
  of the SAME discrete DD transport equations. D1 verifies the *reduction
  itself* is correct — term by term — against a structurally-independent
  hand-posed stencil, BEFORE it is compared to `A_diff` (D2). Without D1,
  D2 is a comparison of two possibly-both-wrong ORPHEUS assemblies.
- **Config.** CFG-TINY (3-cell uniform slab first — where the 3-point
  stencil is hand-writable; then a 4-cell 2-zone variant for the
  material-interface face coefficient).
- **Measured quantity.** The **matrix** `reduce_P1(assemble(L_sn + C_sn −
  S_sn))` — the moment-0/moment-1 reduction of the assembled DD transport
  operator, current eliminated to a scalar stencil — element-wise against
  the analytical reduction `T_hand` computed by hand from the moment
  equations + the DD closure `ψ_face = ½(ψ_in+ψ_out)`.
- **Tolerance.** `assert_allclose(atol=1e-12, rtol=1e-12)` — NOT
  `array_equal` (lessons L16: the reduction's sum order ≠ the hand-posed
  order; principled-equivalent, reduction-depth ULP, vv §bit-identity
  crit-3). The **derivation-side sanity pins** are separate assertions on
  the *same* matrix (reuse the diffusion suite spellings verbatim —
  `tests/diffusion/test_operators.py::test_column_sum_conservation`,
  `::test_m_matrix_sign_pattern`): `1ᵀ A_lo = Σ_a` per group (column-sum
  conservation, exact `atol=1e-12`); M-matrix sign pattern (positive
  diagonal, non-positive off-diagonal, diagonal-dominant ⟹ non-negative
  inverse ⟹ positive flux); `(−∇·D∇)@1 = 0` (constant-flux null space of
  the leakage part).
- **Stabiliser analysis.** Measured functional = the **matrix** (object-
  level, per vv Mode 12 remedy: "pin the OBJECT, not a functional of it").
  Invariance group = **trivial (identity)** — element-wise matrix equality
  has no stabiliser, so it is Mode-12-safe by construction. Contrast the
  *wrong* design: gating on `eig(A_lo)` (invariant under similarity +
  transpose) or on a k-functional (the #290 lesson: `k=λ_max(A⁻¹F)` is
  invariant under factor-swap AND transpose) would be blind to the
  reduction errors this gate exists to catch. The column-sum pin is a
  *balance* functional (invariant under per-term cancellation, anti-#8) —
  hence it is a SANITY pin, NOT the primary gate; the element-wise matrix
  equality is the primary.
- **Mutation rows** (each reddens under `-O`):
  - M-D1-DFACE — D face-interpolation swap (arithmetic vs harmonic mean) →
    element-wise reds (interior stencil off).
  - M-D1-REMOVAL — Σ_r vs Σ_t confusion in the removal coefficient →
    element-wise + column-sum pin red.
  - M-D1-CLOSURE — wrong DD closure weight (½ → some other) → element-wise
    reds; **AND** (the L16 one-source proof) the mutation is applied to the
    **shared** face-closure source, so it MUST also red the existing SN
    sweep/matvec suites. If ONLY D1 reds and the sweep suite stays green, a
    **twin-path stencil exists → STOP, fix, file ERR-NNN** (the whole point
    of consistent DSA is one shared closure).
- **Marker/level.** `foundation` (algebraic-reduction invariant; the
  sanity pins are software invariants) + `verifies("dsa-consistent-loworder"
  TBD)` on the hand-posed anchor (closed-form → the physics claim).
- **Runtime.** Fast (tiny matrices) — **inner loop**.

#### D2 — the SN reduction ↔ landed `A_diff` matrix-consistency characterization

- **Purpose.** Answer the roadmap R4 question *from data*: is the landed
  #290 `A_diff` the consistent partner (Alcouffe recovered computationally),
  or does it differ, and *how*? This gate is a **characterization harness**
  that emits the *structured* matrix diff — it is NOT finalizable as a
  pass/fail CI gate until the user rules (R4), because the RT0/harmonic-mean
  `A_diff` may legitimately differ from the DD moment reduction in the
  interior.
- **Config.** CFG-TINY (uniform + 2-zone).
- **Measured quantity.** The structured difference `Δ = reduce_P1(assemble(
  SN_DD)) − assemble(A_diff)` on `DiffusionMesh.from_material_mesh(sn_mesh)`
  (same axes/materials/BCs by construction). Decomposed into three blocks:
  **interior stencil** `Δ[1:-1, 1:-1]`, **boundary rows** `Δ[{0,N-1}, :]`,
  and a **global scaling** probe (`reduce/A_diff` element-wise ratio on the
  interior).
- **Tolerance / operational "≡ up to a boundary row".**
  - *"≡"* verdict: `assert_allclose(Δ, 0, atol=1e-12)` everywhere.
  - *"≡ up to a boundary row"* verdict (operational): the interior block
    matches — `assert_allclose(Δ[1:-1,1:-1], 0, atol=1e-12)` — **AND** the
    support of the nonzero diff is confined to boundary rows:
    `set(nonzero_rows(Δ)) ⊆ {0, N-1}`. This is the "consistent interior,
    boundary-closure difference" case.
  - *"≢"* verdict: neither holds → the harness **prints** the interior-diff
    structure (is it a constant multiple? a `1/(3Σ_t)` vs consistent-D
    factor? a per-face pattern?) and the gate is left as an
    `xfail(strict=False, reason="R4 ruling pending — bring diff to user")`
    until the user rules derived-stencil-vs-partial-consistency.
- **Stabiliser analysis.** Object-level matrix diff → trivial invariance
  group (Mode-12-safe). The real hazard is NOT a stabiliser but
  **structural independence** (lessons L16): if the SN reduction and
  `A_diff` share the face-closure primitive, `Δ=0` could hold because both
  inherit the same closure bug. **D1's hand-posed anchor is the independence
  guarantee** — D2 says "they agree", D1 says "and the agreed-on value is
  right". Never credit D2 alone.
- **Mutation rows.**
  - M-D2-LOWORDER-D — inconsistent D in the low-order operator (arithmetic
    face mean where consistent-DD wants harmonic, or `1/(3Σ_t)` where the
    reduction wants the DD-consistent coefficient) → the interior block of
    Δ becomes nonzero; the harness reports it as an interior-stencil
    difference (this is the R4 "≢" signal, not necessarily a bug — but it
    changes the verdict).
  - M-D2-BOUNDARY — one face's boundary row wrong → Δ support extends
    beyond `{0,N-1}` OR the boundary block magnitude jumps; caught by the
    boundary-row decomposition. (This same error is a **rate/stability** bug
    downstream, caught by D12 — the boundary is where the naive spike
    diverged.)
- **Marker/level.** `foundation` (object consistency) — no `verifies()`
  until the R4 ruling names the low-order operator's provenance.
- **Runtime.** Fast — **inner loop** (but the gate's *verdict* is gated on
  the R4 user ruling; ships as characterization + xfail until then).

---

### Phase 3b — accelerator correctness (the FP-invariance battery + object gates + controls)

#### D3 — FP-invariance: SI+DSA on an ANISOTROPIC config (the #215-proof, ℓ≥1)

- **Purpose.** The load-bearing correctness gate for the *value-changing*
  bug class: the DSA-accelerated converged flux MUST equal the plain-SI
  fixed point. Run on a config that **breaks the isotropic-reflective
  degeneracy** where the σ_r-fold is accidentally exact (vv Mode 9; the
  #215 trap shipped 46–56% anisotropic errors while isotropic tests passed).
  Includes **anisotropic scattering (ℓ≥1)** to prove the moment DSA does not
  touch — DSA accelerates only φ_0, so the ℓ≥1 moments must be *identical*.
- **Config.** CFG-ANISO-VAC (primary) and CFG-ANISO-HET (redundant
  material-contrast stress). Both vacuum, both ≥2G, both `SigS[1]≠0` built
  on the `Mixture` constructor directly.
- **Measured quantity.** The **full angular flux** `ψ_DSA(n,i,g)` vs
  `ψ_SI(n,i,g)` — NOT merely the scalar φ_0. (Asserting only φ_0 would be
  blind to a P-prolongation that leaks onto ℓ≥1.) Both solved to a tight
  ρ-honest residual `tol=1e-10` (the driver's `SourceIteration` STOP, the
  *equation residual* — NEVER the increment; Signature 9 / ρ-blind stopping,
  else near-critical configs pass a non-converged iterate silently).
- **Tolerance.** `assert_allclose(ψ_DSA, ψ_SI, rtol=SAFETY × tol)` with
  `SAFETY=10`, `tol=1e-10` ⟹ `rtol≈1e-9`. Justification (lessons L7): both
  paths converge the SAME operator to residual `< tol`; the flux agreement
  is bounded by the residual→error amplification `1/(1−ρ)`, absorbed into
  `SAFETY` (the accelerated ρ≈0.22 is small; the plain-SI ρ≤0.9 is the
  binding one, ×10 headroom). Read `tol` off the run config (single source
  of truth, `assert_regression` style), never hardcode. **Two-anchor
  (lessons L2):** the invariance assertion (a) is necessary-not-sufficient
  (it cannot tell you the shared FP is *correct*); pair with (b) — on a
  degenerate-to-homogeneous variant, `ψ_SI` matches the existing SN
  closed-form `Q/Σ_t` (fixed-source) / k_inf (eigenvalue) anchor. Do NOT
  re-mint the anchor; reference the existing SN L1 gate.
- **Stabiliser analysis.** The functional (`ψ_DSA ≡ ψ_SI`, element-wise)
  has no algebraic stabiliser. The Mode-9 blind spot is a **config-induced
  coincidence**, not an algebraic invariance: on isotropic-flux regimes the
  σ_r-fold's error (`Σ_s0·I` vs `Σ_s0·P_iso`) is annihilated, so the wrong
  operator has the same FP as the right one — the config *is* the
  stabiliser. CFG-ANISO-VAC/HET exit that stabiliser (anisotropic flux ⟹
  `P_iso ψ ≠ ψ`). **The disjoint-blind-spot pairing:** D3 is blind to
  R/P/`A_diff` errors (correction→0 makes them FP-invisible) — those are
  caught by D7/D8/D11–D13. D3 and the rate gates have DISJOINT blind spots;
  that disjointness is *why both exist* (the task's explicit requirement).
- **Mutation rows.**
  - **M-D3-SIGMAR** (the reason this gate exists) — wire the within-group
    solve as the σ_r-fold `A_wg.solve(S_residual·ψ + q)` (`S_foldable =
    Σ_s0·P_iso` realized as the diagonal σ_r-sweep) → `ψ_DSA` diverges
    46–56% from `ψ_SI` on CFG-ANISO-VAC/HET → RED. On the isotropic-
    reflective box it would be GREEN (the designed-green degeneracy — the
    proof the config choice matters). **#215 has no ERR entry** (catalog
    ends ERR-069); D3 is the FIRST catcher — `catches("ERR-070" TBD)` and
    file the ERR when the σ_r-fold is formally re-caught by this gate.
  - **M-D3-P-ONTO-L1** — P prolongs the isotropic correction onto ℓ≥1
    moments (not just φ_0) → the ℓ≥1 angular flux moves → RED *only* because
    the measured quantity is the full ψ, not φ_0 (the config-blindness
    closure for the isotropic-snapshot trap, lessons L1).
- **Marker/level.** `l2` (integration: multi-group, heterogeneous,
  accelerator ∘ solver) + `verifies("dsa-fixed-point-invariance" TBD)` +
  `catches("ERR-070" TBD)`.
- **Runtime.** Moderate (two converged solves). CFG-ANISO-VAC is the
  **inner-loop** row; CFG-ANISO-HET is `slow`.

#### D4 — FP-invariance: Krylov-preconditioned-by-DSA (the second posture)

- **Purpose.** The Krylov posture (DSA as a right/left preconditioner inside
  `KrylovAcceleration`, folding #200's intent) has a DIFFERENT plug point
  (the scipy `_as_scipy_linop` preconditioner slot, generalized per
  operator_algebra.rst:1594) and a different failure surface (a broken
  preconditioner slows GMRES but GMRES is splitting-invariant, so it can
  MASK a preconditioner bug by still converging to the right FP). D4 proves
  the Krylov+DSA converged flux equals the plain-Krylov FP.
- **Config.** CFG-ANISO-VAC + CFG-ANISO-HET (same as D3).
- **Measured quantity.** Full angular flux `ψ_{Krylov+DSA}` vs
  `ψ_{Krylov}` (plain GMRES on `(L+C−S−B)`). GMRES `tol=1e-10`.
- **Tolerance.** `assert_allclose(rtol=1e-9)`. GMRES converges to the true
  residual (not increment), so no Signature-9 caveat; but assert on the
  reconstructed angular flux, not the moment iterate.
- **Stabiliser analysis.** GMRES is splitting-invariant ⟹ the plain-Krylov
  FP is preconditioner-independent, so D4's *invariance* is nearly
  tautological for a CORRECT operator — which is exactly why D4 alone cannot
  validate the preconditioner (a broken preconditioner still yields the
  right FP, just slowly). **D4's teeth are in the paired rate gate D13**:
  the Krylov+DSA iteration count must DROP vs plain-Krylov; a broken
  preconditioner reds D13 (no speedup / GMRES stagnation), not D4. Stated
  plainly: **D4 proves safety, D13 proves the preconditioner works.**
- **Mutation rows.**
  - M-D4-PRECOND-SIGN — sign-flip the preconditioner correction → GMRES
    still converges to the right FP (D4 GREEN) but iteration count balloons
    → D13 RED. (Documents the D4/D13 division of labor — a rate bug, not a
    value bug.)
  - M-D4-PRECOND-DEAD — preconditioner returns zero (identity precond) →
    D4 GREEN, D13 RED (no speedup over plain Krylov). Confirms D4 is
    value-safe-blind by design.
- **Marker/level.** `l2` + `verifies("dsa-fixed-point-invariance" TBD)`.
- **Runtime.** Moderate; CFG-ANISO-VAC inner-loop, HET `slow`.

#### D5 — FP-invariance on a DIAGONAL cubature (angular-schedule interaction)

- **Purpose.** Any angular-schedule dependence (a Gauss-Seidel within-group
  angular sweep; a moment-0 restriction that hardcodes the product-
  quadrature weight layout) is invisible on an axis-aligned `product`
  quadrature where each face sees one octant. The ERR-056 shared-face class
  (vv Mode 9 (b)) requires a **diagonal/level-symmetric** cubature to
  expose it. D5 runs FP-invariance on a 2-D level-symmetric config.
- **Config.** CFG-DIAG-2D (2-D Cartesian, ≥2G, het, level-symmetric or
  lebedev, vacuum, ℓ≥1).
- **Measured quantity / tolerance.** As D3 (full angular flux, `rtol=1e-9`).
- **Stabiliser analysis.** The axis-aligned `product` quadrature is the
  stabiliser (each face = one octant ⟹ no shared-face coupling). Level-
  symmetric breaks it (shared faces couple octants). A DSA restriction that
  assumed product-weight structure moves the FP here but not on product.
- **Mutation rows.**
  - M-D5-PRODUCT-WEIGHTS — R's moment-0 hardcodes the product-quadrature
    per-ordinate weight indexing → on level-symmetric the restricted source
    is mis-weighted → (because correction→0, the FP is *still* unchanged!) —
    so this actually reds **D7** (conservation) not D5. *Correction:* D5's
    genuine mutation is a within-group **angular Gauss-Seidel schedule** bug
    (if the build uses one) that changes the FP on shared-face cubatures →
    D5 RED, product GREEN. If the build uses no angular schedule, D5
    degenerates to a redundant 2-D FP-invariance check — keep it as the 2-D
    coverage row and note the schedule-mutation is inactive.
- **Marker/level.** `l2` + `verifies("dsa-fixed-point-invariance" TBD)`.
- **Runtime.** `slow` (2-D converged solves).

#### D6 — correction→0 at convergence (the correctness-safety property, + the dead-R closure)

- **Purpose.** Prove the property the whole partition (§0) rests on: the
  low-order correction δφ → 0 at the fixed point. This is *why* DSA is
  correctness-safe by construction.
- **Config.** CFG-ANISO-VAC (≥2G, ℓ≥1) — the correction must vanish even
  where the flux is anisotropic.
- **Measured quantity.** `‖δφ^{(final)}‖ / ‖φ^{(final)}‖` at convergence
  (the last-iterate correction norm, scalar).
- **Tolerance.** `< SAFETY × tol` (the correction scales with the
  scattering residual, which is `< tol` at the stop). **PAIRED (the blind-
  spot closure):** `‖δφ^{(first)}‖ / ‖φ‖ > 1e-2` (the correction is
  NON-trivial early). Without the pairing, a **dead R that returns 0**
  makes δφ≡0 every iterate and passes "correction→0" *trivially* — the
  gate's Mode-12 blind spot (the measured norm's stabiliser contains
  "δφ identically zero"). The first-iterate non-triviality assertion exits
  that stabiliser.
- **Stabiliser analysis.** `‖δφ‖→0` is a scalar-norm functional whose
  stabiliser contains **every δφ that is small for the wrong reason**
  (dead R, dead P, dead `A_diff.solve` returning 0). The first-iterate
  lower bound is the object-level companion that pins δφ is *alive*; the
  conservation gate D7 pins δφ is *right*.
- **Mutation rows.**
  - M-D6-DEAD-R — R returns 0 → δφ≡0 → "correction→0" GREEN (trivially) but
    first-iterate lower bound RED. Proves the pairing is load-bearing.
  - M-D6-NO-VANISH — feed the correction from the *external* source `q`
    (which does NOT vanish) instead of the scattering residual → δφ does not
    →0 → "correction→0" RED. Proves the gate catches a correction that
    breaks the safety property (would make DSA value-changing).
- **Marker/level.** `l1` + `verifies("dsa-correction-vanishes" TBD)`.
- **Runtime.** Fast — **inner loop**.

#### D7 — R conservation object gate `⟨1, R r⟩_V = ⟨1, r⟩_{V,Ω}` (+ the anti-#8 per-cell closure)

- **Purpose.** The restriction R (moment-0 of the angular residual → scalar
  diffusion source) must conserve the integrated residual. This is the
  *object-level* catcher for the "R missing quadrature weights" error that
  the FP gates are blind to (correction→0).
- **Config.** CFG-TINY (small enough to hand-verify; **S4 GL — non-uniform
  weights load-bearing**: an equal-weight quadrature nulls the
  missing-weights mutation).
- **Measured quantity + the discrete identity spelled out.** For an angular
  residual `r_{n,i,g}` (ordinate n, cell i, group g), the restriction is the
  moment-0 with quadrature weights: `(R r)_{i,g} = Σ_n w_n r_{n,i,g}`.
  Conservation:

  ```
  ⟨1, R r⟩_V  =  Σ_i V_i (R r)_{i,g}  =  Σ_i V_i Σ_n w_n r_{n,i,g}
              =  ⟨1, r⟩_{V,Ω}
  ```

  The gate asserts LHS == RHS with the SAME `V_i` and `w_n` on both sides
  (a random het `r` per group).
- **Tolerance.** `assert_allclose(rtol=1e-12)` (a telescoping-exact
  identity; near machine precision). **PAIRED (anti-#8 closure):** a
  **per-cell delta-source distribution** check — feed `r = δ_{n0,i0}` (one
  ordinate, one cell) and assert `R r == w_{n0}` at cell `i0`, `0`
  elsewhere. The conservation *sum* is a balance functional (invariant under
  a weight-swap-between-cells that preserves the total — Mode-12/anti-#8);
  the delta-source pins the *distribution*, catching the compensating error
  the sum cannot see.
- **Stabiliser analysis.** `⟨1, R r⟩` is a balance/telescoping functional;
  its stabiliser = every per-cell/per-ordinate error that cancels in the
  global sum (the canonical anti-#8 blind spot — "particle balance holds by
  construction"). The delta-source distribution check is the OUTSIDE-the-
  stabiliser companion. The boundary arm (the ℓ=0 half-range moment) is
  checked separately in D8 under the `|Ω·n|w` metric.
- **Mutation rows.**
  - M-D7-NO-WEIGHTS — `R r = Σ_n r_{n,i}` (drop `w_n`) → conservation sum
    RED on the non-uniform GL quadrature (GREEN on equal-weight — the
    config-blindness note). AND the delta-source check reds (`R δ = 1 ≠
    w_{n0}`).
  - M-D7-WEIGHT-SWAP — swap weights between two cells (preserves the total)
    → conservation sum GREEN (compensating), delta-source distribution RED.
    Proves the pairing is load-bearing.
  - M-D7-WRONG-VOLUME — R uses `V_i` inconsistent with the transport mesh's
    → conservation RED (the `V_i` on the two sides differ).
- **Marker/level.** `foundation` (balance identity — software/math
  invariant) + `verifies("dsa-restriction-conservation" TBD)`.
- **Runtime.** Fast — **inner loop**.

#### D8 — R/P adjoint-consistency object gate (the PG-frame law)

- **Purpose.** The restriction R and prolongation P must be an adjoint pair
  under the angular measure (R = moment-0 analysis face, P = isotropic
  reconstruction face of the Petrov–Galerkin frame — the cross-domain-
  attacker's 3-P0 deliverable). This is the object-level catcher for the
  "P missing 1/(4π)-family normalization" error the FP gates are blind to.
- **Config.** CFG-TINY (small matrices; **bulk** arm) + a boundary-face
  fixture (the **boundary** arm).
- **Measured quantity.** Two arms:
  - *Bulk:* `P ≡ R^adjoint` under the angular Gram — spelled as the frame
    identity `R₀ = M₀.H` (frame.rst `R.H` adjoint-for-free machinery).
    Element-wise on the small dense realizations `P.as_matrix()` vs
    `R.as_matrix().T` (or the frame-Gram-weighted transpose, per the frame's
    declared metric — the convention the cross-domain-attacker fixes).
  - *Boundary:* the ℓ=0 **half-range** moment of the angular trace under the
    **shared `|Ω·n| w` metric** (the #290 ruling-2 partial-current seam) —
    `P_boundary ≡ R_boundary^adjoint` under `⟨·,·⟩_{|Ω·n|w}`.
- **Tolerance.** `assert_allclose(atol=1e-12)` element-wise (object-level).
- **Stabiliser analysis.** Spelled as the **matrix** `P − R^T` (object) →
  trivial invariance group (Mode-12-safe). The wrong design — checking a
  single inner product `⟨P e, ψ⟩ = ⟨e, R ψ⟩` for one `(e, ψ)` — is a scalar
  functional whose stabiliser contains normalization errors that happen to
  cancel for that pair; the full-matrix identity has no such escape. Spell
  it as the matrix, not a sampled bilinear form.
- **Mutation rows.**
  - M-D8-P-NORM — P missing the `1/(4π)`-family normalization (or carrying
    `1/W` where the metric expects unweighted) → `P ≠ R^T` under the
    declared measure → RED. (This same error is a **rate** bug downstream —
    a mis-normalized P is a bad accelerator — caught also by D13; the object
    gate localizes it.)
  - M-D8-BOUNDARY-METRIC — the boundary arm uses the wrong half-range metric
    (Euclidean instead of `|Ω·n|w`) → boundary adjoint RED. Guards the #290
    ruling-2 seam.
- **Marker/level.** `foundation` (frame law) + `verifies("dsa-rp-frame-
  adjoint" TBD)`. **Depends on the cross-domain-attacker's PG-frame
  landing** — if the R/P pair is spelled ad-hoc rather than as a frame,
  this gate has no `R.H` to assert against; flag as a blocking dependency.
- **Runtime.** Fast — **inner loop**.

#### D9 — no-masking control (DSA converges fast to the SAME wrong answer)

- **Purpose.** The converse of D3–D5. Seed a bug in the **within-group
  transport operator** and confirm DSA-accelerated SI converges (fast) to
  the SAME wrong flux as plain SI. The accelerator must not *launder* a real
  solver bug into a plausible-but-different answer — DSA changes the RATE,
  never the VALUE, in BOTH directions (D3–D5: don't move the right FP; D9:
  faithfully reproduce the wrong FP).
- **Config.** CFG-ANISO-VAC with a **seeded transport mutation** (e.g., a
  scaled Σ_s, or a perturbed streaming coefficient — a bug in `A = L+C−S−B`,
  NOT in the DSA machinery). ≥2G.
- **Measured quantity.** `ψ_{DSA on buggy A}` vs `ψ_{SI on buggy A}` —
  both converged; assert they AGREE (both wrong, identically).
- **Tolerance.** `assert_allclose(rtol=1e-9)` (same FP, different rate) +
  a `!=` guard that both DIFFER from the un-mutated correct FP by ≫ tol
  (proves the seeded bug is active — a Mode-10 activation check; without it
  a null mutation passes vacuously).
- **Stabiliser analysis.** The functional (`ψ_DSA ≡ ψ_SI` on the buggy
  operator) is element-wise → no algebraic stabiliser. The activation guard
  (both differ from the correct FP) is the anti-Mode-10 companion — the
  seeded mutation must be in a term CFG-ANISO-VAC activates (a Σ_s or
  streaming perturbation is active on anisotropic het flux).
- **Mutation rows.**
  - M-D9-LAUNDER — if DSA's correction were (wrongly) built from the
    *external* source rather than the residual, DSA on the buggy operator
    would converge to a DIFFERENT flux than SI (it would "correct toward"
    the diffusion answer) → D9 RED. Confirms D9 catches a laundering
    accelerator.
  - Control: the un-seeded run (both agree AND match the correct FP) — the
    no-op leg that pins the asymmetry.
- **Marker/level.** `l2` + `verifies("dsa-fixed-point-invariance" TBD)`.
- **Runtime.** Moderate — `slow` (two converged solves + control).

#### D10 — the #215 routing structural guard

- **Purpose.** A structural tripwire: the within-group solve inside the DSA
  loop must NOT route through the σ_r-fold spelling (`A_wg.solve(S_residual)`
  with the σ_r-sweep realizing `S_foldable = Σ_s0·P_iso`). This is the
  *design-time* guard that reddens if anyone wires the trap the docstrings
  at `scattering.py:968-971` warn against — complementing D3 (which catches
  it *numerically* after the fact).
- **Config.** N/A (structural / call-graph, not numerical).
- **Measured quantity.** A Mode-11-style **in-process wrap** (autouse
  fixture / monkeypatch counter) on `CollisionOperator`'s σ_r-sweep path:
  assert the DSA within-group solve does NOT invoke the σ_r-fold realization
  (counter == 0), and/or an **AST tripwire** forbidding the
  `A_wg.solve(S_residual...)` pattern in the DSA module. The within-group
  solve MUST route through the honest `(L+C−S−B)` resolvent (SI or Krylov),
  never the diagonal-in-angle σ_r removal.
- **Tolerance.** Exact (counter == 0 / AST match absent). Uses an explicit
  `pytest.fail` / `raise`, NOT a bare assert (Mode 8, `-O`).
- **Stabiliser analysis.** Structural gate — no numerical functional, no
  stabiliser. The tooth is "does the σ_r-fold path execute?" A green D10
  with a fired wrap-counter on the honest resolvent (counter > 0 there)
  proves the DSA loop exercises the correct path (the Mode-11 "gate executes
  the intended line" discipline).
- **Mutation rows.**
  - M-D10-WIRE-FOLD — wire the within-group solve as the σ_r-fold →
    counter > 0 (or AST match) → RED. (Pairs with D3's M-D3-SIGMAR: D10
    catches it structurally at design time, D3 catches the 46–56% error
    numerically.)
- **Marker/level.** `foundation` + `sentinel` (cheap always-on structural
  tripwire) + `catches("ERR-070" TBD)`.
- **Runtime.** Fast — **inner loop / sentinel** (run without `-O` per the
  sentinel convention).

---

### Phase 3c — rate + stability

#### D11 — ρ measurement protocol vs the 0.2247c Fourier bound

- **Purpose.** Verify the *acceleration works as theory predicts*: the
  measured spectral radius of SI+DSA is the c-independent ρ ≈ 0.2247c
  (Adams & Larsen 2002 §II), on the classic homogeneous problem.
- **Config.** CFG-HOMOG-INF (homogeneous, periodic/infinite slab, S8 GL).
  **1-group is LEGITIMATE here** — this is a RATE claim (the iteration-map
  spectral radius), which is flux-shape-independent BY DESIGN (AGENT.md §5;
  the Cardinal Rule bars only 1G *eigenvalue* claims — the claim layer is
  declared "rate", not "eigenvalue"). The Fourier analysis is itself
  1-group.
- **Measured quantity + estimator + asymptotic-regime detection.**
  - **Estimator:** ρ_est = geometric mean of the tail of successive
    **residual** ratios `r_{n+1}/r_n` (read off `SourceIteration`'s
    `residual_history`; or equivalently the driver's `contraction_ratios`
    from the typed `FluxDisplacement`). **Use the residual, NOT the
    increment** — Signature 9: the increment ratio is ρ-honest only away
    from c→1; the residual ratio is the ρ-honest estimator (the driver's
    STOP already uses the equation residual).
  - **Asymptotic-regime detection:** discard the transient — take the last
    K ratios where `std(tail) / mean(tail) < 0.05` (the ratio has
    stabilized). If it never stabilizes within `max_iter`, the estimate is
    invalid → fail loudly (do not report a transient ρ).
- **Tolerance.** `|ρ_est − 0.2247·c| < 0.03` in the asymptotic regime, at
  c ∈ {0.5, 0.9}. Justification: the 0.2247 is the continuous-Fourier bound;
  discrete DD + finite S_N introduces an O(few %) correction, so a 0.03
  absolute band (not a tight rtol) is the honest tolerance — a tighter band
  would be reference contamination (the bound is asymptotic, the measurement
  is discrete). This is a **precision claim on ρ**; pair with D13's
  c-independence (ρ flat across c).
- **Stabiliser analysis.** ρ measures the iteration operator's spectral
  radius — invariant under similarity/transpose of the iteration operator,
  which is fine (ρ IS the quantity). **ρ's blind spot: it does NOT see the
  converged VALUE** — a wrong-but-fast DSA has a good ρ and a wrong flux. So
  D11 is PAIRED with D3 (FP-invariance): D11 says "fast", D3 says "and
  right". Disjoint blind spots (the task's requirement) — D11 cannot catch a
  value bug, D3 cannot catch a rate bug.
- **Mutation rows.**
  - M-D11-INCREMENT-STOP — measure ρ from the increment ratio instead of
    the residual → at c=0.99 the estimator under-reports (Signature 9) →
    the c-sweep (D13) shows a spurious ρ collapse near c→1 → RED. Proves the
    estimator choice is load-bearing.
  - M-D11-NO-ACCEL — replace DSA with identity (no correction) → ρ_est ≈ c
    (plain SI) ≠ 0.2247c → RED. Sanity that the gate measures acceleration.
- **Marker/level.** `l1` (closed-form/semi-analytical Fourier bound) +
  `verifies("dsa-spectral-radius" TBD)`.
- **Runtime.** `slow` (needs many iterations to reach the asymptotic
  regime); the c ∈ {0.9} single point is the inner-loop subset.

#### D12 — reflective-BC stability gate (the historical failure mode, positive)

- **Purpose.** A *positive* gate on exactly the config where the naive
  (inconsistent FD) DSA spike DIVERGED: fully-reflective, c up to 0.99. With
  *consistent* DSA this MUST converge with ρ in band. This is the load-
  bearing stability discriminator — the reflective box has no leakage to
  damp the diffusion↔transport consistency error, so an inconsistent
  low-order operator diverges here (and ONLY here — vacuum is forgiving).
- **Config.** CFG-REFLECT (fully-reflective slab AND 2-D box), ≥2G, c ∈
  {0.9, 0.95, 0.99}.
- **Measured quantity.** (a) **Converges** — SI+DSA reaches the residual
  tol within `max_iter` (the boolean the spike failed); (b) ρ_est in band
  (< 0.3, one-sided — reflective ρ may exceed the 0.2247c vacuum value
  slightly but must stay bounded well below 1); (c) **value pairing** — the
  converged flux matches the plain-SI reflective FP (`rtol=1e-9`), so
  "converges" is not "converges to the wrong stable point" (a stable-but-
  wrong FP is the Mode-12 blind spot of a boolean convergence gate).
- **Tolerance.** Convergence within `max_iter` (assert `len(
  residual_history) < max_iter` AND final residual `< tol`); ρ_est `< 0.3`;
  value `rtol=1e-9` vs plain SI.
- **Stabiliser analysis.** "Converges" is a boolean whose stabiliser
  contains "converges to a wrong stable FP" — closed by the value pairing
  (c). ρ_est < band is a rate functional blind to the value — also closed
  by (c). Together: converges + fast + to the right answer.
- **Mutation rows.**
  - M-D12-INCONSISTENT-LO — use an inconsistent low-order operator (a plain
    FD Laplacian `D=1/(3Σ_t)` + 5-point stencil, NOT the DD-consistent
    reduction — the exact naive-spike bug) → **DIVERGES** on CFG-REFLECT
    (residual grows) → RED. This is the single highest-value stability
    mutation — it reproduces the historical failure. (The SAME inconsistent
    operator is GREEN on vacuum — the reason the gate is reflective.)
  - M-D12-BOUNDARY-ROW — a wrong reflective boundary row in the low-order
    operator → diverges or converges to a wrong FP on reflective → RED
    (caught by the convergence + value pairing). Pairs with D2's
    M-D2-BOUNDARY (object-level) — D12 is the runtime consequence.
- **Marker/level.** `l1` (the stability property + closed-form ρ<1
  convergence) + `verifies("dsa-reflective-stability" TBD)`.
- **Runtime.** `slow` (c=0.99 reflective SI is thousands of plain-SI
  iterations for the plain-SI reference; the DSA path is fast but the
  reference is slow). Inner-loop subset = c=0.9 slab only.

#### D13 — c→1 × optical-thickness sweep + the iteration-count table

- **Purpose.** The characterization artifact (the teaching table for the
  theory page + demo) AND the CI rate-regression gate. Demonstrates
  c-independence (the ρ≈0.2247c headline) and bounds the partial-consistency
  degradation across cell optical thickness (the consistency stress axis
  where partially-consistent schemes fail).
- **Config.** CFG-HOMOG-INF (c-sweep) + CFG-THICK (thickness sweep) +
  CFG-REFLECT (BC discriminant). c ∈ {0.9, 0.95, 0.99, 0.999};
  Σ_t·Δx ∈ {0.1, 1, 10, 100}; BC ∈ {vacuum, reflective}. Four solver
  columns: SI, SI+DSA, Krylov, Krylov+DSA.
- **Measured quantity.** The iteration-count table `n(solver, c, thickness,
  BC)` = `len(residual_history)` per cell of the table.
- **Tolerance / what is a HARD CI gate vs a characterization table:**
  - **HARD CI gate (the inner-loop subset — one c, one thickness, both
    BCs):** `n_{SI+DSA} < 0.5 · n_SI` at c=0.9 (conservative speedup floor;
    the measured spike was 0.05–0.13×, so 0.5× is a loose-but-real
    regression tripwire — a rate regression that nothing measures is a
    silent un-gated claim). **c-independence:** `n_{SI+DSA}(c=0.99) < 1.5 ·
    n_{SI+DSA}(c=0.9)` (the count stays ~flat as c→1 — the whole point).
    **Krylov+DSA:** `n_{Krylov+DSA} < n_{Krylov}` (the preconditioner helps
    — the D4 teeth).
  - **Characterization table (the full sweep, `slow`, docs/demo):** the full
    `n(solver, c, thickness, BC)` grid, emitted to the theory page. No hard
    per-cell tolerance — it is the artifact, not a gate. (But the reflective
    column MUST show finite counts everywhere — a `inf`/`max_iter` entry
    there is a D12 stability failure surfacing in the table.)
- **Stabiliser analysis.** Iteration count is a rate functional, blind to
  the value — paired with D3/D12 (value) throughout. The conservative 0.5×
  floor is chosen so FP-non-associativity / machine-jitter in the count
  (±1 iteration) cannot false-red the gate, while a real regression (DSA
  stops accelerating) reds hard.
- **Mutation rows.**
  - M-D13-NO-SPEEDUP — any accelerator-quality bug (M-D4-PRECOND-SIGN,
    M-D8-P-NORM, M-D2-LOWORDER-D) → the count does not drop → CI gate RED.
    This is where the R/P/`A_diff` error class (FP-invisible) is finally
    caught. **D13 is the catcher of last resort for the entire rate-only
    bug class.**
  - M-D13-C-DEPENDENT — a partially-consistent low-order operator → the
    count climbs with c (loses c-independence) → the c-independence gate
    RED. Distinguishes consistent (flat) from partial (c-climbing).
- **Marker/level.** CI subset: `l1` + `verifies("dsa-spectral-radius"
  TBD)`. Full sweep: characterization (NO `verifies()`, `slow`) — its
  FAILURE is the signal to update the table (lessons L5: a negative-
  regression characterization gate).
- **Runtime.** `slow` (full grid); the one-c/one-thickness/both-BC subset is
  the inner-loop CI gate.

---

## 4. The mutation matrix — which gate reds each canonical implementation error

For every plausible error, the gate that reds AND the confirmation the claim
is OUTSIDE that gate's stabiliser. **The partition (§0) is visible in the
last column: exactly ONE error (the σ_r-fold) reds the FP gates; the rest
are object/rate-caught.**

| # | Implementation error | Primary catcher | Stabiliser check (why it's caught, not designed-green) | FP gates blind? |
|---|---|---|---|---|
| 1 | **R missing quadrature weights** | **D7** (conservation sum) + **D13** (rate degrades) | D7 uses non-uniform GL weights (equal-weight nulls it); the delta-source distribution check catches the weight-swap variant | **YES** (correction→0) |
| 2 | **P missing 1/(4π)-family normalization** | **D8** (P ≠ R^T, object matrix) + **D13** (rate) | D8 spells the full matrix `P−R^T`, not a sampled bilinear form (which could cancel) | **YES** |
| 3 | **sign flip in the correction update** | **D11/D13** (rate degrades) + **D12** (may diverge on reflective) | if it still converges, FP unchanged (correction→0) → rate is the only signal; the count balloons or diverges | **YES** |
| 4 | **correction applied to the wrong iterate** | **D13** (rate) | correction→0 leaves the FP unchanged; only the rate moves | **YES** |
| 5 | **low-order operator: inconsistent D / removal** | **D1/D2** (matrix structure) + **D12** (diverges on reflective) + **D13** (c-dependent count) | D1's element-wise object gate has no stabiliser; D12's reflective divergence is the runtime consequence (vacuum is forgiving — the config choice matters) | **YES** |
| 6 | **boundary rows wrong on one face** | **D2** (boundary-row block) + **D8** (boundary adjoint) + **D12** (reflective) | D2 decomposes the diff into interior/boundary blocks; D8's `|Ω·n|w` metric arm; D12 is where boundary bugs diverge | **YES** |
| 7 | **the σ_r-fold spelling (#215)** | **D3/D4/D5** (FP-invariance, 46–56% wrong) + **D10** (routing guard, structural) | D3 runs the ANISOTROPIC config that exits the isotropic-flux stabiliser where the fold is exact; D10 catches it at design time | **NO — this is the one that reds them** |
| 8 | **dead R / dead P / dead A_diff (returns 0)** | **D6** (first-iterate non-triviality) + **D7/D8** (object) | D6's paired lower bound exits the "δφ small for the wrong reason" stabiliser | **YES** (correction→0 GREEN trivially) |

**The load-bearing reading:** if the plan-of-record shipped only the
FP-invariance gates (D3–D5), errors 1–6 and 8 — *seven of eight* — would
ship silently. The object gates (D1/D2/D7/D8) and rate gates (D11–D13) are
the battery's actual coverage of the accelerator; the FP gates cover only
the σ_r-fold class.

---

## 5. Gate → build-phase mapping

| Build phase | Gates | Rationale |
|---|---|---|
| **3a** (consistency by derivation) | **D1** (reduction ≡ stencil + sanity pins), **D2** (reduction ↔ `A_diff` characterization) | Land with the assembled DD operator + the moment reduction; D2's *verdict* is gated on the R4 user ruling (ships as characterization + xfail until then). |
| **3b** (the accelerator) | **D3, D4, D5** (FP-invariance, both postures + diagonal cubature), **D6** (correction→0), **D7, D8** (R/P object gates), **D9** (no-masking), **D10** (routing guard) | Land with the R/P frame + the low-order solve + the SI+DSA / Krylov+DSA wiring. D6/D7/D8/D10 are the object/property gates that land WITH the machinery; D3/D4/D5/D9 need the full loop. D8 is BLOCKED on the cross-domain-attacker's PG-frame; D1/D2 are BLOCKED on the literature-researcher's consistent-stencil coefficients. |
| **3c** (rate + stability) | **D11** (ρ vs 0.2247c), **D12** (reflective stability), **D13** (c→1 × thickness table) | Land with the tuned accelerator; D13's full sweep is the theory-page/demo artifact. The LD arm-2 decision (roadmap R5) reads off D13's measured table. |

**Inner-loop fast subset** (the pre-merge quick gate, `-m "not slow"`): D1,
D2, D6, D7, D8, D10 (all fast object/property gates) + D3 CFG-ANISO-VAC + D11
c=0.9 + D12 c=0.9 slab + D13 one-c/one-thickness/both-BC. **Full sweep**
(`slow`): D3 HET, D4/D5, D9, D11 full, D12 c=0.99, D13 full grid.

---

## 6. Honest scope — what this battery does NOT cover

- **Curvilinear DSA — OUT (blocked on #282).** DSA here is slab/Cartesian DD
  (arm 1). The curvilinear R/P (spherical/cylindrical moment restriction
  through the angular-redistribution closure) is a documented seam, blocked
  on the #282 curvilinear-transpose work. No gate here touches curvilinear;
  a `.. note::` on the theory page must state it.
- **The LD arm-2 — demand-driven (roadmap R5), NOT gated here.** LD's ρ
  under the arm-1 low-order operator is MEASURED by D13's sweep and the
  table is brought to the user; the build-LD-consistent decision (M4S /
  DG-diffusion) is deferred. If arm-2 lands, it needs its own D1/D2/D3
  analogs (the LD moment reduction ≠ the DD one — a distinct consistency
  gate).
- **Eigenvalue (k-effective) DSA — the outer iteration, NOT covered.** This
  battery gates the *within-group* (inner) acceleration. DSA inside a
  k-eigenvalue outer iteration adds the `n_inner=None` gap (the SN
  eigenvalue path does not surface the inner count — `si_convergence_rate_
  verification.md`); the rate gates D11/D13 run on FIXED-SOURCE problems
  precisely to sidestep that gap. A k-eigenvalue DSA rate claim needs the
  inner-count surfaced first (a separate issue).
- **TSA (transport-synthetic acceleration, #5) and Chebyshev/Wielandt
  outer acceleration (#75) — sibling accelerators, out of scope.** DSA is
  the diffusion-synthetic inner accelerator only.
- **Validation (L3) — not applicable.** DSA changes the iteration, not the
  equation; there is no experimental datum for an accelerator. The battery
  is L0/L1/L2 + foundation only.
- **The R4 low-order-operator provenance ruling.** D2 CHARACTERIZES the
  reduction↔`A_diff` diff but does not DECIDE derived-stencil-vs-partial-
  consistency — that is a user ruling the D2 harness informs. Until ruled,
  D2 is a characterization gate, not a pass/fail CI gate.

---

## 7. Design tensions for the plan-of-record to decide

1. **The 3a consistency gate is a characterization, not a verdict, until
   R4.** D2 can only become a hard equality CI gate AFTER the user rules
   (from D2's structured diff) whether the landed RT0/harmonic-mean `A_diff`
   IS the consistent DD partner, or whether a NEW consistent stencil must be
   derived. The spec ships D2 as a characterization harness + xfail; the
   plan-of-record must sequence the R4 ruling BEFORE 3b (the low-order solve
   in 3b consumes whichever operator R4 blesses).

2. **The shared-face-closure structural-independence question (L16).**
   Consistent DSA REQUIRES the SN sweep and the low-order operator to share
   the face-closure — which makes D2's `Δ=0` partly a bit-identity
   *inheritance*, not an independent check. The independence MUST come from
   D1's hand-posed anchor. The plan must confirm the closure is genuinely
   single-sourced (the L16 one-source proof: a closure sign-flip reds BOTH
   the SN sweep suite AND D1) — if a twin stencil exists, STOP and fix
   before D2 means anything.

3. **D8 depends on the R/P being a Petrov–Galerkin FRAME, not an ad-hoc
   pair.** If the cross-domain-attacker's 3-P0 frame does not land (R/P
   spelled as bare functions), D8 has no `R.H` adjoint-for-free to assert
   against, and the adjoint-consistency gate degrades to a hand-written
   inner-product check (weaker, sample-dependent). The plan should treat the
   PG-frame as a HARD prerequisite for D8, not optional polish.

4. **The measurand for the rate gates assumes fixed-source.** D11/D13 read
   `len(residual_history)` from `SourceIteration` on FIXED-SOURCE problems.
   The eigenvalue path leaves `n_inner=None` (open gap). If the plan wants a
   k-eigenvalue DSA rate claim, it must first surface the inner count — flag
   as a dependency, do not silently gate on an unmeasured count.

5. **ρ estimator: residual vs increment (Signature 9).** The spec mandates
   the RESIDUAL ratio (the driver's ρ-honest STOP), not the increment. If
   the build's ρ-measurement helper reads the increment (the `contraction_
   ratios` from `FluxDisplacement`), it under-reports near c→1 and the
   c-independence gate (D13) sees a spurious ρ collapse. The plan must wire
   the estimator off the residual history, not the displacement.

6. **Is there an ERR entry for the σ_r-fold?** No (catalog ends ERR-069).
   D3/D10 are the FIRST catchers. The plan should file **ERR-070** (Mode 9,
   "the σ_r-fold ships 46–56% anisotropic error while isotropic tests pass")
   when D3 first reddens on the seeded σ_r-fold, and attach `catches(
   "ERR-070")` to D3 and D10 — mutation-verified, not topic-tagged
   (vv "a catches() marker is a coverage claim").

7. **Tolerance single-source.** Every FP-invariance tolerance is `SAFETY ×
   tol` off the run-config `tol` (regression-tolerance-design memory), never
   hardcoded. The plan must expose the DSA loop's `tol` as the shared source
   the generator and the gate both read.
