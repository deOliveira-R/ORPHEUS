# R-1 Step 4 Step G — verification plan

**Branch**: `refactor/sn-operator-algebra` @ `6cd8355`
**Date**: 2026-05-22
**Author**: test-architect (dispatched per `r1_step4_session2_followup.md`
§2.3; `subagent-handoff-protocol` proactive-trigger row "operator-
algebra carve crossing subsystem boundaries").

**Companion documents** (READ FIRST — this plan does not duplicate
their content):

* `.claude/plans/r1_step4_session2_followup.md` (the master plan;
  §2.3-§4 set Step G's scope).
* `.claude/plans/r1_step4_g_dependency_audit.md` — callgraph audit
  (CORRECTION section + SURPRISE-1..8 are the source of truth for
  what retires).
* `.claude/plans/r1_step4_g_convention_crosswalk.md` — the API
  contract this plan defends.
* `.claude/lessons.md` L17-L21 (durable lessons from session 1).
* `.claude/skills/vv-principles/error_catalog.md` ERR-049 + ERR-050
  (the bug classes this plan structurally prevents from re-
  introduction).

---

## Scope and inputs

Step G's structural carves: **G0 + G1 + G2**.

| Step | Action | Verification weight |
|---|---|---|
| **G0** | NEW native-shape unified matvec (consumes typed `AngularFlux` + `BoundaryFlux` natively; no `eq_map` slot lookup) | HIGH — gates G3b/G3c |
| **G1** | `_solve_fixed_source_krylov` → `KrylovAcceleration` + typed AngularFlux (1-D only) | HIGH — public MMS API |
| **G2** | `_solve_fixed_source_si` → `SourceIteration` + typed AngularFlux (1-D + 2-D Cartesian) | HIGH — public MMS API |
| **G3a..G3g** | Mechanical retirements (delete helpers, classes, decoders, monkey-patch fixtures) | LOW — sentinel: existing green tests after structural carves still pass |
| **G3h** | Cartesian build_equation_map retirement — DEFERRED to Phase A | — |

This plan covers G0/G1/G2 deeply; G3a-G3g get a tests-retire-with-
symbol table + a sentinel "all extant green tests still green" gate.

---

## V&V hierarchy applied to Step G

Per `vv-principles` §"Hierarchical claim taxonomy" — every claim in
this plan is annotated with its layer + pillar + structural reference.

| Claim level | What gets pinned for Step G | Structural reference (the pillar) |
|---|---|---|
| **L0 — term-level** | Per-ordinate flat-flux residual (catches ERR-026 if Carlson seed wires regress); /sum_w producer-side normalisation; eq_map slot-lookup ≡ quad-derived-mask equivalence | Closed-form (uniform ψ ⇒ streaming=0 + uniform Q/Σ_t fixed-source identity); structurally independent of any iteration scheme |
| **L0 — software invariant** | Native-matvec ≡ legacy packed-face-slot matvec on the same (typed) input ψ; AngularFlux round-trip identity; capability advertisement (`L`, `C`, `L+C` advertise apply/solve correctly) | Closed-form (matvec is purely linear; double-application equivalence is an algebraic identity) |
| **L1 — equation verification** | Multi-group homogeneous reflective ⇒ `k_inf = νΣ_f/Σ_a` (closed-form A⁻¹F dominant eigenvalue); slab MMS O(h²) convergence; curvilinear MMS O(h²) convergence (under the assumption ERR-026's `xfail` markers retire structurally with Step G — see §3.4 below); 2-D Cartesian MMS O(h²) | Closed-form k_inf (homogeneous-reflective limit); MMS for spatial convergence (NOT eigenvalue per §"Three pillars") |
| **L1 — semi-analytical** | Case singular-eigenfunction reference (1G 2-region slab, S8); trajectory_resolvent Variant α Green's function (cylinder MR 2G) | Semi-analytical pillar — both references already shipped in `tests/sn/test_l1_standoff_slab_cylinder.py` |
| **L2 — integration** | 11 SN regression snapshots (`tests/sn/regression/snapshots/*.npz`) | Self-convergence regression: bit-identity (slab) OR principled-equivalence (curvilinear, FP-noise bounded) |
| **foundation** | InvertibleOperator (L+C).solve round-trip identity; (L+C).apply = L.apply + C.apply; OperatorSum reduction invariants | Closed-form algebraic identity |

**L4 (code-to-code) appears nowhere in this plan as a correctness
claim.** The §"Equivalence gate" rows that compare native-matvec
output to legacy packed-face-slot output are SANITY checks; their
backing L0/L1 anchors are named explicitly on each row.

**1-group degeneracy guard** (Cardinal Rule 6): every k_inf /
eigenvalue gate here uses ≥2 groups (2eg or 4eg). 1G appears ONLY in
the existing snapshot suite (e.g. `cyl_1g_homogeneous_LS4_dd_n20`)
where 1G is the regression — not the verification claim.

---

## G0 — native-shape unified matvec verification

**What ships**: a new function (working name
`transport_operator_matvec_native(psi: AngularFlux, sn_mesh, sig_t,
*, boundary_flux: BoundaryFlux) → AngularFlux`) that consumes typed
`AngularFlux.values` `(N, ng, nx, ny)` + `BoundaryFlux.xmax_face`
`(N, ng)` + `BoundaryFlux.xmin_face` (slab) directly, derives
inflow-ordinate masks from `sn_mesh.quad`, and returns native-shaped
cell + face residuals. Replaces `transport_operator_matvec_unified`'s
packed-face-slot signature.

### G0.1 — Existing L0/foundation tests that pin today's packed-face-slot interface

These tests CURRENTLY exercise the packed-face-slot interface. Each
row notes whether the test migrates (rewires to the new native
matvec) or retires (pinned a contract that dissolves).

| Test file:line | Class/test | Tier | Pins | Step G action |
|---|---|---|---|---|
| `tests/sn/test_unified_matvec_slab.py:80` | `test_unified_slab_zero_psi_gives_zero` | L0/foundation | Linearity: matvec(0)=0 | **MIGRATE** — rewire call site to native matvec; assert still holds |
| `tests/sn/test_unified_matvec_slab.py:94` | `test_unified_slab_constant_psi_gives_sigma_t` | L0 | Σ_t·ψ subtraction on uniform ψ (slab) | **MIGRATE** — load-bearing L0; structural |
| `tests/sn/test_unified_matvec_slab.py:184` | `test_unified_slab_l1_homogeneous_kinf_2g[nx]` | L1 | 2G homogeneous-reflective slab k_inf via unified matvec | **MIGRATE** — load-bearing L1 anchor |
| `tests/sn/test_unified_matvec_slab.py:243` | `test_unified_slab_differs_from_legacy_fd_O_h` | L1 | Unified WDD ≠ legacy 1st-order FD (the order-of-accuracy delta) | **RETIRE at G3c** — no legacy FD path to compare against post-G3c |
| `tests/sn/test_unified_matvec_sphere.py:126` | `test_unified_constant_psi_gives_sigma_t` | L0 | Σ_t·ψ subtraction on uniform ψ (sphere) | **MIGRATE** |
| `tests/sn/test_unified_matvec_sphere.py:145` | `test_unified_zero_psi_gives_zero` | L0 | Linearity (sphere) | **MIGRATE** |
| `tests/sn/test_unified_matvec_cylinder.py:240` | `test_unified_cylinder_matches_hand_reference` | L1 | Hand-derived 1-cell cylinder matvec | **MIGRATE** |
| `tests/sn/test_unified_matvec_cylinder.py:280` | `test_unified_cylinder_zero_psi_gives_zero` | L0 | Linearity (cylinder) | **MIGRATE** |
| `tests/sn/test_unified_matvec_cylinder.py:295` | `test_unified_cylinder_constant_psi_gives_sigma_t` | L0 | Σ_t·ψ on uniform ψ (cylinder) | **MIGRATE** |
| `tests/sn/test_unified_matvec_cylinder.py:437` | `test_unified_cylinder_l1_mr_2g_trajectory_resolvent` | L1 | 3-region 2G ABA cylinder vs trajectory-resolvent Green's function | **MIGRATE** — load-bearing L1 (semi-analytical pillar) |
| `tests/sn/test_unified_matvec_cylinder.py:489` | `test_unified_cylinder_l1_homogeneous_kinf_2g` | L1 | 2G cylinder homogeneous k_inf | **MIGRATE** |
| `tests/sn/test_b1pp_verification.py:114` | `test_b1pp_lplusc_is_full_rank` | L0 | Typed (L+C) full rank on slab/sphere/cylinder | **STAYS** (uses typed leaves; not packed-face-slot) |
| `tests/sn/test_b1pp_verification.py:167` | `test_b1pp_constant_flux_collapses_to_collision` | L0 | Native B1'' face-aware contract | **STAYS** — already native |
| `tests/sn/test_b1pp_verification.py:224` | `test_b1pp_lplusc_gmres_converges_fp_noise` | L1 | GMRES on (L+C) at FP-noise on all 3 geometries | **STAYS** |
| `tests/sn/test_b1pp_verification.py:279` | `test_b1pp_decode_encode_roundtrip` | L0 | `_with_traces` round trip | **STAYS through G3d**; retires when `_with_traces` retires (G3d) |

### G0.2 — NEW tests for the native-shape interface

The crosswalk's "Path-forward convention contract" (one-page table,
§Axes 3 + 6) defines the API. New L0 tests pin that contract
directly.

| New test (file:: name) | Tier | Pins | Catches |
|---|---|---|---|
| `tests/sn/test_native_matvec.py::test_inflow_ordinate_mask_derived_from_quad_matches_eq_map_slot_lookup[slab/sphere/cylinder]` | L0 | The `quad.mu < 0` mask equals the legacy `eq_map.face_outer_ordinate` slot ordering (modulo permutation). **Sentinel for the structural fix.** | Convention drift between quad-derived and slot-derived inflow masks |
| `tests/sn/test_native_matvec.py::test_native_matvec_consumes_typed_angular_flux[slab/sphere/cylinder]` | L0 | Shape contract: input `AngularFlux(values=(N,ng,nx,ny))` + `BoundaryFlux(xmax_face=(N,ng),xmin_face=(N,ng))` → output `AngularFlux(values=(N,ng,nx,ny))` with face residual in `boundary.xmax_face` | Typed/packed shape swap |
| `tests/sn/test_native_matvec.py::test_native_matvec_zero_psi_gives_zero[slab/sphere/cylinder]` | L0 / foundation | Linearity at zero | Sign-flip / variable-swap |
| `tests/sn/test_native_matvec.py::test_native_matvec_uniform_psi_gives_minus_sigma_t_psi[slab/sphere/cylinder]` | L0 | Uniform ψ ⇒ streaming term = 0 ⇒ matvec(ψ) = -Σ_t·ψ at every cell, every ordinate | ERR-006 / ERR-007 / ERR-026 re-introduction (per-ordinate flat-flux residual mode) |
| `tests/sn/test_native_matvec.py::test_native_matvec_per_ordinate_flat_flux_residual[sphere/cylinder]` | L0 / structural-fix sentinel | Per-ordinate flat-flux balance (Signature 1 in numerical-bug-signatures) | ERR-026 family — invisible to global-balance tests; the canonical curvilinear diagnostic |
| `tests/sn/test_native_matvec.py::test_native_matvec_face_residual_appears_only_on_boundary_face[slab/sphere/cylinder]` | L0 | Native matvec emits non-zero face residual at outer (slab: also inner) boundary face, zero elsewhere | B1'' face-state convention drift |
| `tests/sn/test_native_matvec.py::test_native_matvec_apply_is_linear[slab/sphere/cylinder]` | L0 / foundation | `matvec(a·ψ₁ + b·ψ₂) = a·matvec(ψ₁) + b·matvec(ψ₂)` | Non-linearity (e.g. accidental conditional branch) |

**Tag conventions**: every new test carries `@pytest.mark.l0` or
`@pytest.mark.foundation`, the appropriate `@pytest.mark.verifies(...)`
linking to `docs/theory/sn_operator_algebra.rst` labels (which Phase 6
must add), and where applicable `@pytest.mark.catches("ERR-049",
"ERR-050", "ERR-026")` for regression-catch traceability.

### G0.3 — Equivalence gate: legacy packed-face-slot vs native shape

| Geometry | Bit-identity expected? | Reason | Test |
|---|---|---|---|
| **Slab** | **NO** — principled-equivalent | The legacy `SNStreamingOperator.apply` slab path uses the Cartesian decoder (no face) + cell-centre proxy; the native matvec uses B1'' face-aware (the test_b1pp_verification.py contract). The legacy path is O(h) at the boundary; the native path is O(h²). **The difference IS the intended improvement** (per audit augmentation A Difference 1). | New `tests/sn/test_g0_equivalence.py::test_g0_native_vs_legacy_slab_principled_equivalent` — relative diff scales as O(h) at boundary cell; assert pre-bridge ULP-bounded match interior, structured O(h) at outer face. Three-pillar attestation required (see §4 below) |
| **Sphere** | **YES** — bit-identical (audit A Difference 3) | The inward hemisphere at the outer face IS the M-M coupled-pole seed | New `tests/sn/test_g0_equivalence.py::test_g0_native_vs_legacy_sphere_bit_identical` — `np.array_equal` |
| **Cylinder** | **NO** — principled-equivalent (audit A Difference 2; ~4e-3 at nx=40) | Legacy cell-centre proxy ≠ B1'' face-aware. The difference IS the cylinder twin-path bug (`tests/sn/test_l1_standoff_slab_cylinder.py::test_cylinder_l1_sweep_vs_krylov_twin_path` xfail). G0 fixes manifestation #6 by construction (Pattern 2). | New `tests/sn/test_g0_equivalence.py::test_g0_native_vs_legacy_cylinder_principled_equivalent` — relative diff bounded by O(h) at the outer cell; matches B1'' L1 convergence. Three-pillar attestation required |
| **2-D Cartesian** | **N/A** — defer to Phase A (G3h) | 2-D Cartesian native matvec is out of scope (SURPRISE-2). Legacy path retains | (no test) |

**Acceptance**: green G0.1 + G0.2 + G0.3. G0 ships ONLY after this
suite is green; G1/G2 depend on G0.

---

## G1 — `_solve_fixed_source_krylov` onto `KrylovAcceleration` verification

**What ships**: `_solve_fixed_source_krylov` body collapses to:

```
(L+C-S) @ ψ = q_ext_typed   # KrylovAcceleration.solve
```

2-D Cartesian raises `NotImplementedError` (mirrors `_solve_krylov`'s
guard at the eigenvalue inner; per SURPRISE-4/5 the principled
landing zone is `_solve_fixed_source_si` via G2).

### G1.1 — Existing L1 MMS gates that MUST continue passing

The L1 MMS suite is the public-API contract for `solve_sn_fixed_source`.
Step G must NOT change the convergence rate the MMS measures.

| Test file:line | Test | Tier | Current state | Step G expectation |
|---|---|---|---|---|
| `tests/sn/test_mms.py:34` | `test_sn_1d_slab_mms_converges_second_order` | L1 | PASSING — `orders > 1.9`, `1e-8 < errors[-1] < 1e-4` | **PASS unchanged.** Slab uses SI (no curvilinear default flip) — G2 covers this; G1 does not affect slab path |
| `tests/sn/test_mms_curvilinear.py:67` | `test_sn_spherical_mms_converges_second_order` | L1 | **xfail strict** (ERR-026 PARTIAL — pre-asymptotic L2 magnitude) | **CONTINUE xfail strict** — Step G does not close ERR-026's pre-asymptotic-magnitude branch; the convergence rate must remain O(h²). New L1 sentinel below pins the rate independent of the absolute-magnitude assertion |
| `tests/sn/test_mms_curvilinear.py:124` | `test_sn_cylindrical_mms_converges_second_order` | L1 | **xfail strict** (same as sphere) | **CONTINUE xfail strict** with same rate sentinel below |
| `tests/sn/l1_analytical/test_mms_curvilinear_aniso_dd_convergence.py:95` | `test_sn_spherical_aniso_mms_converges_second_order` | L1 | **xfail strict** (ERR-026 PARTIAL — anisotropic ansatz, pre-asymptotic) | **CONTINUE xfail strict** + rate sentinel |
| `tests/sn/l1_analytical/test_mms_curvilinear_aniso_dd_convergence.py:151` | `test_sn_cylindrical_aniso_mms_converges_second_order` | L1 | **xfail strict** | **CONTINUE xfail strict** + rate sentinel |
| `tests/sn/test_mms_aniso.py:38` | `test_sn_p1_aniso_mms_converges_second_order` | L1 | PASSING | **PASS unchanged** (slab path; SI default; G2 covers) |
| `tests/sn/test_mms_heterogeneous.py:56` | `test_sn_heterogeneous_mms_converges_second_order` | L1 | PASSING | **PASS unchanged** (slab path) |
| `tests/sn/test_mms_2d.py:45` | `test_sn_2d_cartesian_mms_converges_second_order` | L1 | PASSING | **PASS unchanged** (2-D Cartesian → SI per SURPRISE-5; covered by G2) |
| `tests/sn/test_mms_2d.py:90` | `test_sn_2d_cartesian_2g_heterogeneous_mms_converges_second_order` | L1 | PASSING | **PASS unchanged** (G2 2-D path) |

**The xfail-strict curvilinear MMS gates are the augmentation A
answer**: Step G does NOT unblock ERR-026's pre-asymptotic-magnitude
branch (the absolute-error assertion `1e-8 < errors[-1] < 1e-3`).
Issue #195 is the tracking artifact. Step G's responsibility is to
**preserve the convergence rate** (the `orders > 1.9` assertion) —
which the rate-only sentinel below pins independently.

### G1.2 — Existing L1 analytical k_inf gates (multi-group, heterogeneous)

These gates use `solve_sn` (eigenvalue), NOT `solve_sn_fixed_source`,
BUT they exercise the same `_solve_krylov` carve and the same
operator-algebra path. G1's structural changes must not degrade them.

| Test file:line | Test | Tier | Current state |
|---|---|---|---|
| `tests/sn/l1_analytical/test_kinf_homogeneous.py:125` | `test_kinf_homogeneous[ng_key×coord×inner_solver]` | L1 | 35 pass + 5 xfailed (sphere-4eg-krylov per #200) — closed-form `k_inf = νΣ_f/Σ_a` |
| `tests/sn/l1_analytical/test_kinf_homogeneous.py:184` | `test_kinf_homogeneous_spectrum[ng_key×coord×inner_solver]` | L1 | Same — closed-form A⁻¹F dominant eigenvector |
| `tests/sn/test_invertible_operator.py:730` | `test_si_carve_recovers_analytical_kinf[coord×ng_key]` | L1 | PASSING — 6 cases (slab/sphere/cylinder × 2eg/4eg). **Already the L1 sentinel for ERR-049** (`@pytest.mark.catches("ERR-049")`) |
| `tests/sn/test_invertible_operator.py:657` | `test_fixed_source_homogeneous_reflective_recovers_q_over_sigma[coord]` | L1 | PASSING — `ψ_n = q_n / Σ_t` on every geometry |
| `tests/sn/test_krylov_curvilinear_precond_safety.py:187` | `test_identity_preconditioner_recovers_kinf[coord]` | L1 | PASSING — **Already the L1 sentinel for ERR-050** (`@pytest.mark.catches("ERR-050")`) |
| `tests/numerics/test_iteration.py:506` | `test_keigenvalue_matches_solve_sn_2g_slab` | L1 | PASSING — compares typed `KEigenvalue` (typed adapters around (L+C)) against legacy `solve_sn` |

**`test_keigenvalue_matches_solve_sn_2g_slab` — Augmentation D
answer**: this is the load-bearing L1 anchor that compares typed
Krylov against legacy `SNStreamingOperator.apply`. Per audit
augmentation A, slab is **principled-NON-equivalent** between the
two paths (the legacy uses Cartesian decoder + cell-centre proxy,
not B1'' face-aware). However, this test runs on 2G slab where the
algebraic delta is small and the SI inner converges to the analytical
fixed point; current tolerance is `abs(keff - expected) < 1e-9`.
**Audit augmentation A's "cylinder ≈ 4e-3 at nx=40" delta does NOT
apply here** (this test is slab, not cylinder; slab's algebraic
delta is dominated by boundary-cell O(h) and reaches `< 1e-9` at the
n=10 mesh used). Therefore:

- **At G1 ship**: test continues passing (G1 doesn't change the
  iteration's fixed point on slab; SI still converges the SAME way).
- **At G3f (SNStreamingOperator retirement)**: test is **MIGRATED**
  to compare typed `KEigenvalue` against typed `(L+C)-driven
  KrylovAcceleration`. The legacy bundle is gone; the comparison
  becomes typed-vs-typed. Tolerance can tighten (both sides use the
  same operator algebra at this point).

### G1.3 — NEW L0/L1 pins for G1 (convention catchers, structural-fix sentinels)

| New test | Tier | Pins | Catches | Structural-fix sentinel for |
|---|---|---|---|---|
| `tests/sn/test_fixed_source_g1.py::test_g1_krylov_no_eq_map_in_callstack[slab/sphere/cylinder]` | L0 / foundation | Stack trace of `_solve_fixed_source_krylov` for 1-D geometries contains NO call to `build_equation_map*` or `solution_to_angular_flux*` or `_build_rhs_*` | Re-introduction of any retired packed-face-slot symbol | The structural-fix sentinel for the retirement itself |
| `tests/sn/test_fixed_source_g1.py::test_g1_krylov_2d_cartesian_raises_notimplemented` | L0 | 2-D Cartesian input to `solve_sn_fixed_source(inner_solver="krylov")` raises `NotImplementedError` with explicit deferral message | Silent regression of 2-D Cartesian Krylov (audit SURPRISE-4) | The structural defer for Phase A |
| `tests/sn/test_fixed_source_g1.py::test_g1_krylov_external_source_per_ord_density_not_iso[slab/sphere/cylinder]` | L0 | When called with `external_source = np.full((N,ng,nx,ny), q_iso/sum_w)` (per-ord density convention), recovered ψ = q_iso/(W·Σ_t)/... matches `PerOrdinateSource.from_isotropic(q_iso, sn_mesh)` exactly. Producer-side normalisation only | Re-introduction of consumer-side `/sum_w` bridge | ERR-049 sentinel (re-affirm) |
| `tests/sn/test_fixed_source_g1.py::test_g1_krylov_curvilinear_mms_rate_only[sphere/cylinder]` | L1 / structural-fix sentinel | At `nx ∈ {20, 40, 80}` the L2 error ratios give `orders > 1.9` (rate-only; absolute magnitude tolerance is the xfail branch tracked by #195) | Step G regression of curvilinear MMS rate, INDEPENDENT of ERR-026 magnitude branch | The rate-vs-magnitude split per Augmentation A |
| `tests/sn/test_fixed_source_g1.py::test_g1_krylov_2g_het_slab_keff_matches_g1_si[2eg]` | L1 / structural-fix sentinel | Heterogeneous slab 3-region 2G keff via Krylov inner ≡ via SI inner (both routes through G1 + G2 typed carves) to ≤ 1e-8 rel | Twin-path divergence post-carve (L14 manifestation #6 sentinel for the carve) | Twin-path agreement post-G1 + G2 |
| `tests/sn/test_fixed_source_g1.py::test_g1_krylov_solver_attrs_l_is_gone[slab/sphere/cylinder]` | L0 / foundation | After G3f the `SNSolver` instance does NOT have a `solver.L` attribute that holds `SNStreamingOperator`. (Pre-G3f: skip or xfail.) | Mid-retirement regression where `self.L` lingers | The "self.L removed entirely" pre-condition for G3f (audit) |

### G1.4 — Snapshot regeneration policy (three-pillar attestation)

The 11 SN regression snapshots live at
`tests/sn/regression/snapshots/*.npz`. Per `vv-principles` §"Bit-
identity vs principled-equivalence", any snapshot whose bit-value
changes under G1 must be regenerated with a three-pillar attestation.

**Inventory of snapshots vs Step G impact**:

| Snapshot | Solver path | G1 impact | G2 impact | G0 impact | Regenerate? |
|---|---|---|---|---|---|
| `slab_2g_homogeneous_dd_n20` | `solve_sn` (eigenvalue, SI inner) | None (eigenvalue inner unchanged) | None | None | **NO** — bit-identical |
| `slab_2g_3reg_dd_n40` | `solve_sn` (eigenvalue, SI inner) | None | None | None | **NO** |
| `sphere_2g_homogeneous_dd_n20` | `solve_sn` (eigenvalue, default→Krylov post-Step-D) | None | None | None | **NO** — unchanged |
| `sphere_2g_3reg_dd_n40` | `solve_sn` (eigenvalue, Krylov inner) | None | None | None | **NO** — unchanged |
| `cyl_1g_homogeneous_LS4_dd_n20` | `solve_sn` (eigenvalue, Krylov inner) | None | None | None | **NO** |
| `cyl_1g_homogeneous_product_dd_n20` | `solve_sn` (eigenvalue) | None | None | None | **NO** |
| `cyl_2g_3reg_LS4_dd_n40` | `solve_sn` (eigenvalue, Krylov inner) | None | None | None | **NO** |
| `slab_2g_p1_aniso_dd_n20` | `solve_sn` (eigenvalue, SI inner; P1 scattering) | None | None | None | **NO** |
| `sphere_2g_p1_aniso_dd_n20` | `solve_sn` (eigenvalue, Krylov inner; P1) | None | None | None | **NO** (currently a known pre-existing issue per `r1_step4_session2_followup.md` §"Reality check" — 5 failing on baseline) |
| `2d_1g_LS4_dd_15x15` | `solve_sn` (eigenvalue, 2-D Cartesian wavefront sweep, SI) | None | None | None | **NO** — 2-D wavefront sweep untouched |
| `slab_fixed_source_dd_n20` | `solve_sn_fixed_source` (default SI on slab) | None (slab default is SI) | **G2 path** — same fixed point, FP-noise drift | None | **Bit-identical expected**; if not, three-pillar attest |

**Conclusion**: Step G changes the *path* through which fixed-source
problems are solved, but the discrete fixed point of `(L+C-S)·ψ =
q_ext` is unchanged — the typed `KrylovAcceleration` and the typed
`SourceIteration` converge to the SAME ψ that the legacy code converged
to (per `vv-principles` §"Bit-identity vs principled-equivalence"
criterion 2: structurally-independent reference is the analytical
homogeneous-reflective `k_inf` already pinned by
`test_kinf_homogeneous` AND the closed-form `ψ = Q/Σ_t` already
pinned by `test_fixed_source_homogeneous_reflective_recovers_q_over_sigma`).

Drift, if any, is FP-non-associativity at `iteration_count × ULP`
scale — well below the curvilinear `rtol=5e-6` and the slab
`rtol=1e-12` bounds.

**If any snapshot's bit value DOES change** (FP-non-associativity
exceeds the regression tolerance), the three-pillar attestation:

1. **Pillar 1 — L0 streaming-equilibrium / flat-flux invariant**: at
   the touched snapshot's configuration, the per-ordinate flat-flux
   residual on a uniform-medium variant of the same problem MUST be
   < 1e-13 (the `test_g0_per_ordinate_flat_flux_residual` family on
   the matching geometry).
2. **Pillar 2 — Pomraning pole isotropy** (where applicable —
   sphere/cylinder snapshots only): the pole/μ=−1 face flux MUST
   match the M-M Carlson coupled-pole seed identity to FP-noise.
3. **Pillar 3 — Structurally independent cross-check**: the snapshot's
   k_inf (homogeneous-reflective rows) MUST match the closed-form
   `νΣ_f / Σ_a` to `< 1e-10` (the `test_kinf_homogeneous` row at the
   same `ng_key × coord`).

**NEVER** regenerate a snapshot whose value drifted because "the new
number is closer to the old one." Drift WITHOUT three-pillar
attestation is rejected.

### G1.5 — NotImplementedError gate for 2-D Cartesian Krylov

**Augmentation C answer**: Today there is NO existing regression test
that exercises 2-D Cartesian fixed-source via the Krylov inner. The
2-D MMS suite (`tests/sn/test_mms_2d.py`) calls
`solve_sn_fixed_source(...)` without specifying `inner_solver`, so it
gets the default which routes to SI for 2-D Cartesian (the curvilinear
default flip in `solve_sn_fixed_source` does not fire). The 2-D
Cartesian path through Krylov inner has been unreachable from any
test since at least Wave-E.

**Implication**: G1's `NotImplementedError` is a clean defer. NO
existing test regresses; the SURPRISE-5 routing (2-D Cartesian → SI
via G2) preserves all 2-D MMS tests. The principled-equivalence
contract is trivially satisfied: nothing to be principled-equivalent
to.

The new test `test_g1_krylov_2d_cartesian_raises_notimplemented`
pins this contract explicitly so a future regression that silently
re-enables Krylov-on-2D-Cartesian fails loudly.

---

## G2 — `_solve_fixed_source_si` onto `SourceIteration` verification

**What ships**: `_solve_fixed_source_si` body collapses to:

```
SourceIteration(L+C, S, ZeroOperator).solve(q_ext_typed)
```

Geometry-agnostic — slab + sphere + cylinder + 2-D Cartesian all
land here. The 2-D Cartesian path is the SURPRISE-5 landing zone.

### G2.1 — Existing tests that pin today's SI behaviour

Same family as G1.1/G1.2 but the SI-side rows:

| Test | Tier | Step G action |
|---|---|---|
| `tests/sn/test_mms.py::test_sn_1d_slab_mms_converges_second_order` | L1 | **PASS unchanged** (1-D slab routes to G2 via SI default) |
| `tests/sn/test_mms_aniso.py::test_sn_p1_aniso_mms_converges_second_order` | L1 | **PASS unchanged** |
| `tests/sn/test_mms_heterogeneous.py::test_sn_heterogeneous_mms_converges_second_order` | L1 | **PASS unchanged** — heterogeneous slab, the load-bearing redistribution catcher |
| `tests/sn/test_mms_2d.py::test_sn_2d_cartesian_mms_converges_second_order` | L1 | **PASS unchanged** (2-D Cartesian → SI per G2 SURPRISE-5) |
| `tests/sn/test_mms_2d.py::test_sn_2d_cartesian_2g_heterogeneous_mms_converges_second_order` | L1 | **PASS unchanged** — heterogeneous 2-D, multi-group, redistribution + group coupling stress |
| `tests/sn/test_invertible_operator.py::test_si_carve_recovers_analytical_kinf[coord×ng_key]` | L1 | **PASS unchanged** — the existing L1 sentinel for ERR-049 stays the gate |

### G2.2 — NEW L0/L1 pins for G2

| New test | Tier | Pins | Catches |
|---|---|---|---|
| `tests/sn/test_fixed_source_g2.py::test_g2_si_geometry_agnostic_dispatch[slab/sphere/cylinder/cart_2d]` | L0 / foundation | `_solve_fixed_source_si` for each geometry dispatches through `transport_sweep` natively (no per-geometry branch); the SI inner converges to the same fixed point on uniform-medium input | Geometry-specific helper code re-introduction |
| `tests/sn/test_fixed_source_g2.py::test_g2_si_2d_cartesian_landing_zone_kinf` | L0 | 2-D Cartesian fixed-source via G2 SI route recovers `Q/Σ_t` on homogeneous reflective (the `test_fixed_source_homogeneous_reflective_recovers_q_over_sigma` analog for 2-D) | 2-D Cartesian regression at the landing zone for SURPRISE-5 |
| `tests/sn/test_fixed_source_g2.py::test_g2_si_external_source_per_ord_density_not_iso[slab/sphere/cylinder]` | L0 | Same A1-convention pin as G1.3 row, for the SI route | ERR-049 sentinel (SI route) |
| `tests/sn/test_fixed_source_g2.py::test_g2_si_initial_guess_threads_through_explicit_kwarg[slab/sphere/cylinder]` | L0 / foundation | SI's `_solve_with_seed` passes previous iterate via `initial_guess=` to `InvertibleOperator.solve` (not via `rhs(1)` history) — Phase 1.2 contract | Re-introduction of the silent `rhs(1)` history coupling (ERR-050 sentinel for SI) |

### G2.3 — 2-D Cartesian landing-zone tests (SURPRISE-5 from the audit)

The 2-D Cartesian Krylov path raises `NotImplementedError` (G1.5);
2-D Cartesian SI is the answer. The existing 2-D MMS tests
(`test_sn_2d_cartesian_mms_converges_second_order`,
`test_sn_2d_cartesian_2g_heterogeneous_mms_converges_second_order`)
already validate this route — they will continue to pass after G2
because the underlying `transport_sweep` dispatch is unchanged.

The NEW `test_g2_si_2d_cartesian_landing_zone_kinf` adds an L0
homogeneous-reflective sentinel so the landing zone is explicitly
pinned (it's not an MMS, it's a closed-form Q/Σ_t check).

---

## Cross-cutting — `test_phase_c_gates.py` 11 sites (SURPRISE-7)

`tests/sn/test_phase_c_gates.py` has 11 instantiation sites of
`SNStreamingOperator(sn_mesh, sig_t)`. Per audit, they retire at
G3f. Each must either migrate (typed equivalent exists) or retire
(pins a soon-retired contract).

| Test (file line) | Pins | G3f action | Capability-equivalent in path forward |
|---|---|---|---|
| `:115 test_apply_linearity_under_sweep_frame` (×3 geom) | L+C linearity in sweep frame | **MIGRATE** → `StreamingOperator + CollisionOperator` typed sum | `test_b1pp_verification.py::test_b1pp_lplusc_is_full_rank` + new typed linearity test |
| `:201 test_apply_curvilinear_per_ordinate_flat_flux_residual` | ERR-026 catcher (per-ordinate flat-flux on sphere/cylinder) | **MIGRATE** → equivalent on `(L+C).apply` typed | New `tests/sn/test_native_matvec.py::test_native_matvec_per_ordinate_flat_flux_residual` (G0.2 row) |
| `:248 test_apply_apply_transpose_reciprocity_under_sweep_frame` | (L+C) self-adjoint / reciprocity | **MIGRATE** → typed `.H` reciprocity | New `tests/sn/test_operators_apply_typed.py` row for adjoint reciprocity (the `.H` propagation contract from coding-elegance Pattern 1) |
| `:277 test_apply_face_fluxes_match_sweep_recurrence_spherical` | Sphere apply ≡ sphere sweep at face | **MIGRATE** → typed sweep ≡ typed apply via shared SNCellOperator (Phase 1.2 contract) | New `tests/sn/spatial/test_sweep_vs_apply_consistency.py::test_typed_sweep_apply_face_match[sphere]` |
| `:318, :361, :411, :450, :471 BC trace contract suite` (×5) | BC trace at vacuum/reflective/mixed sphere | **MIGRATE** → on typed L (the BoundaryRealizer carries the BC convention; the legacy `op.apply` is a wrapper) | New typed BC-trace tests (a 5-row sub-suite on typed `StreamingOperator`) |
| `:586 test_sweep_curvilinear_per_ordinate_flat_flux_residual` | ERR-026 catcher (sweep side) | **STAYS** (already on typed sweep) | — |

**Net effect**: `test_phase_c_gates.py` shrinks to the 1 sweep-side
test that uses typed primitives. The 10 `SNStreamingOperator`-bearing
tests migrate to the typed equivalents listed in the right column,
which are either new (G0.2 family) or existing (`test_b1pp_verification.py`,
`test_operators_apply_typed.py`).

**Issue #199 body amendment**: add this row inventory + the explicit
"migrate or retire" table BEFORE G3f starts.

---

## Tests retiring at G3f and their replacements

The audit identifies the test files whose existence depends on
`SNStreamingOperator`. At G3f the class deletes; these tests retire
WITH it.

| Test file | Retire / migrate at G3f? | Capability-equivalent path-forward test |
|---|---|---|
| `tests/sn/test_snstreamingoperator.py` (22+ tests) | **DELETE** entirely | Each invariant is covered by typed equivalents: linearity → `test_operators_apply_typed.py`; sweep-frame matvec → `test_b1pp_verification.py` + new G0.2 family; eq_map slot consistency → retires with packed-face-slot layout (no replacement needed) |
| `tests/sn/test_streaming_operator_decomposition.py` (8 sites using `_with_traces`) | **MIGRATE** at G3d (when `_with_traces` retires) or **DELETE** if its pins are subsumed | The 8 `_with_traces` sites test the decoder family. Audit says they pin the typed↔legacy bridge — that contract dissolves entirely. Replace with typed-only round-trip tests in `test_operators_apply_typed.py` if any invariant doesn't already exist there |
| `tests/sn/test_streaming_operator.py::TestCompositionEquivalence` (xfail block at :286-373) | **DELETE the xfail block** at G3f | The comparison was `SNStreamingOperator.apply ≟ (L+C).apply`. No `SNStreamingOperator` exists post-G3f; the block has nothing to compare. **The fact that the xfails flip green during G3c → G3f is the intended improvement, NOT a regression** (audit augmentation A) |
| `tests/sn/test_l1_standoff_slab_cylinder.py` (monkey-patch fixtures at :92, :128, :132) | **DELETE the `_patch_apply_to_unified` context manager + its uses at every test** | Per audit + Augmentation B: the monkey-patches were a pre-Step-G stand-in. Post-G3f, `solve_sn(inner_solver="krylov")` natively routes through typed `(L+C)`, so no patch is needed. The 4-leg L14 standoff tests SURVIVE and use unpatched `solve_sn` directly — pin **structural invariants** (sweep ≡ Krylov on the same operator algebra), NOT the legacy bundle's cell-centre-proxy contract. **The `xfail strict` cylinder twin-path test flips green** — also the intended improvement |
| `tests/sn/test_unified_matvec_slab.py` (monkey-patch sites at :149, :175, :179) | **DELETE the monkey-patch fixtures**; the unified-matvec tests themselves migrate to G0.2 family OR stay if they exercise the typed path | Per augmentation B: these tests pin **structural invariants** (zero-psi, constant-psi, kinf-homogeneous, FD-vs-WDD O(h) divergence). The structural invariants migrate to native matvec via G0.2; the FD-vs-WDD test retires (no legacy FD post-G3c) |
| `tests/sn/test_unified_matvec_cylinder.py` (monkey-patch sites at :399, :427, :431) | Same as slab | Migrate or retire per same rationale |
| `tests/sn/test_unified_matvec_sphere.py` | **STAY** (no monkey-patch; uses typed primitives) | — |
| `tests/sn/spatial/test_apply_matvec_cylinder_invariants.py` (uses `build_equation_map_cylindrical` at :56, :104) | **MIGRATE at G3e** to typed cylinder invariants on `StreamingOperator.apply`-typed | New typed-invariant tests under `tests/sn/spatial/` |
| `tests/numerics/test_iteration.py::test_keigenvalue_matches_solve_sn_2g_slab` | **MIGRATE at G3f** to compare `KEigenvalue(L+C, S, F).solve` against typed `(L+C)-driven KrylovAcceleration`-based `solve_sn` (post-G3f `solve_sn` routes through typed) | Same test name, typed implementation |

**Augmentation B answer**: per-file recommendation —
`test_l1_standoff_slab_cylinder.py`: structural invariants worth
keeping (sweep ≡ Krylov twin-path agreement IS the L14 contract;
that survives the carve and becomes a typed-twin-path test);
`test_unified_matvec_{slab,cylinder}.py`: mostly structural
invariants (linearity, uniform-ψ); the monkey-patch sites are
exclusively pre-Step-G plumbing and delete cleanly. **Do NOT retire
these test files entirely — migrate their invariants.** The
`SNStreamingOperator`-only tests (`test_snstreamingoperator.py`,
`test_streaming_operator_decomposition.py`'s `_with_traces` cluster,
`test_phase_c_gates.py`'s 10 SN-bundle rows) are the ones that
retire wholesale.

---

## Verification gates between steps — the checkpoint sequence

Each step has an explicit green-bar checkpoint that MUST hold before
the next step starts.

| Step | Gate before next step starts | What confirms it |
|---|---|---|
| **G0** | Native matvec ≡ legacy packed-face-slot matvec (bit-identical sphere, principled-equivalent slab+cylinder per §G0.3); ALL G0.1 migrated tests + G0.2 new tests pass | `pytest tests/sn/test_native_matvec.py tests/sn/test_g0_equivalence.py tests/sn/test_unified_matvec_*.py` |
| **G1** | `_solve_fixed_source_krylov` for 1-D geometries calls NO retired symbol; 2-D Cartesian Krylov raises `NotImplementedError`; ALL G1.1 + G1.2 existing L1 gates green; ALL G1.3 new tests pass; G1.4 snapshot policy holds (no regen needed) | `pytest tests/sn/l1_analytical/ tests/sn/test_mms*.py tests/sn/test_invertible_operator.py tests/sn/test_krylov_curvilinear_precond_safety.py tests/sn/test_fixed_source_g1.py tests/sn/regression/` |
| **G2** | `_solve_fixed_source_si` for all geometries routes via `SourceIteration`; 2-D Cartesian lands here; ALL G2.1 existing + G2.2 new tests pass | `pytest tests/sn/test_mms*.py tests/sn/test_fixed_source_g2.py tests/sn/test_invertible_operator.py` |
| **G3a** | `_make_sweep_preconditioner` + `_build_rhs_*` deleted; zero production callers verified by `grep` | `grep -rn "_build_rhs_\|_make_sweep_preconditioner" orpheus/` returns ONLY docstring refs |
| **G3b** | typed `StreamingOperator._apply_typed` + `CollisionOperator._apply_typed` use native matvec; `_ensure_eq_map`, `_eq_map`, `n_unknowns` removed | grep check + `pytest tests/sn/test_operators_apply_typed.py tests/sn/test_invertible_operator.py` |
| **G3c** | `transport_operator_matvec_unified` (packed-face-slot signature) deleted | grep check |
| **G3d** | `pack_with_traces` / `solution_to_angular_flux_with_traces` / `build_equation_map_with_traces` / `EquationMap` deleted | grep check; `tests/sn/test_streaming_operator_decomposition.py` migrated or deleted |
| **G3e** | `solution_to_angular_flux_{spherical,cylindrical}` + `build_equation_map_{spherical,cylindrical}` deleted | grep check; `tests/sn/spatial/test_apply_matvec_cylinder_invariants.py` migrated |
| **G3f** | `SNStreamingOperator` deleted; 11 sites in `test_phase_c_gates.py` migrated (per §"Cross-cutting" table); monkey-patch fixtures in `test_l1_standoff_slab_cylinder.py` + `test_unified_matvec_{slab,cylinder}.py` removed; xfails in `test_streaming_operator.py::TestCompositionEquivalence` deleted; `test_keigenvalue_matches_solve_sn_2g_slab` migrated; **ALL existing L1 MMS gates + L1 kinf gates STILL green** | full `pytest tests/sn tests/numerics tests/sn/l1_analytical tests/sn/regression`; grep check `grep -rn SNStreamingOperator orpheus/` returns nothing |
| **G3g** | `AngularFlux.to_flat_with_traces` + `from_flat_with_traces` deleted | grep check; no production callers |
| **Phase 6 / H** | Sphinx narrative shipped; `sphinx-build docs docs/_build/html -W` clean | sphinx command |

**Cylinder twin-path xfail flip**: between G3c and G3f the audit
augmentation A says the cylinder twin-path xfails in
`test_l1_standoff_slab_cylinder.py::test_cylinder_l1_sweep_vs_krylov_twin_path`
+ `test_cylinder_l1_refinement_both_paths` go from `xfail strict` to
PASSING. This is the **intended improvement**, NOT a regression. The
gate at G3f-1 (right before SNStreamingOperator deletes) MUST detect
this transition and remove the `xfail` markers in the same commit.

---

## Bug-class containment — which ERR-NNN this plan structurally prevents

The crosswalk's bug-class containment table is the canonical view;
this plan adds the test-level enforcement:

| Bug class | Sentinel test in this plan | Step that ships the sentinel |
|---|---|---|
| **ERR-049 — convention drift /sum_w** | `test_g1_krylov_external_source_per_ord_density_not_iso[slab/sphere/cylinder]` + `test_g2_si_external_source_per_ord_density_not_iso[slab/sphere/cylinder]` + existing `test_si_carve_recovers_analytical_kinf` (`@catches("ERR-049")`) | G1.3 + G2.2 |
| **ERR-050 — silent precond fallback** | Existing `test_identity_preconditioner_recovers_kinf` + `test_default_sweep_preconditioner_recovers_kinf_on_slab` (`@catches("ERR-050")` after ERR-050 logged); + new `test_g2_si_initial_guess_threads_through_explicit_kwarg` | G2.2 |
| **ERR-026 — curvilinear sweep WDD wrong fixed point** | Existing `test_per_ordinate_flat_flux_consistency` + `test_spherical_sweep_vs_bicgstab_flat_flux` (Signature 1); new `test_native_matvec_per_ordinate_flat_flux_residual[sphere/cylinder]` | G0.2 |
| **ERR-006/007 — α recursion + ΔA/w** | Same as ERR-026 sentinels (the per-ordinate flat-flux is the canonical diagnostic) | G0.2 |
| **Twin-path divergence (L14 manifestation #6)** | New `test_g1_krylov_2g_het_slab_keff_matches_g1_si[2eg]` (twin-path agreement at L1, multi-group, heterogeneous); existing `test_cylinder_l1_sweep_vs_krylov_twin_path` flips green at G3f | G1.3 + G3f |
| **Convention drift between SI and Krylov producers** | Same A1-convention pins (G1.3 + G2.2) on both inner solvers; they share the SAME `ScatteringOperator.apply` post-A1 | G1 + G2 |
| **2-D Cartesian silent Krylov regression** | New `test_g1_krylov_2d_cartesian_raises_notimplemented` | G1.3 |

---

## Open questions / decisions still to make at G0/G1/G2 design time

These are explicit deferrals — the implementer (main agent or
method-implementer sub-agent) must resolve them during G0/G1/G2 code
design.

1. **`solve_sn_fixed_source` public-API signature** — crosswalk
   §"Public API stays unchanged" defers the decision: keep
   `external_source: np.ndarray` (wrap to `AngularFlux` internally
   at entry) OR tighten to `AngularFlux` over a deprecation cycle.
   **Recommended resolution at G1 design time**: preserve
   bare-ndarray entry since `solve_sn_fixed_source` is L1-MMS-facing
   and external callers exist. Tightening can come in a later wave.
2. **Heterogeneous slab snapshot policy post-G1** — if FP-non-
   associativity drift on `slab_2g_3reg_dd_n40` exceeds the
   `rtol=1e-12` bound (it shouldn't — both paths use SI with the
   same producer-side normalisation post-A1), regenerate with the
   three-pillar attestation per §G1.4. **Recommended**: re-run the
   regression suite immediately after G2 to detect any drift; do
   NOT preemptively regenerate.
3. **Whether `solver.L` removal in `SNSolver.__init__` (G3f) is a
   one-shot delete or a deprecation cycle** — audit notes that
   `orpheus.sn.SNStreamingOperator` is re-exported as PUBLIC API
   from `orpheus/sn/__init__.py`. A clean public-API-breaking
   removal is preferred (Cardinal Rule 1 — principled over lazy);
   but a one-cycle deprecation shim might be necessary for external
   callers per the project's general convention. **Recommended**:
   verify no external callers via `grep` across `derivations/` +
   `scratch/` + `docs/` (the audit already did this — only doc
   refs); proceed with clean deletion + clear migration note in the
   commit message + a Sphinx narrative update.
4. **What to do with the 2-D Cartesian Krylov path** — G1 raises
   `NotImplementedError`; G2 catches 2-D Cartesian via SI. **Open
   question**: should `solve_sn_fixed_source` (the dispatcher at
   `solver.py:1262`) intercept the 2-D Cartesian + Krylov request
   and silently route to SI (warning the user), OR propagate the
   `NotImplementedError` from G1? **Recommended**: propagate;
   silent routing violates the principle of least surprise. Document
   the 2-D Cartesian + Krylov deferral in the docstring + raise a
   clear error.
5. **Whether `test_phase_c_gates.py` becomes one of the test files
   that retires entirely** — after the 10 migration rows + the 1
   keeper, the file has only 1 test. **Recommended**: relocate the
   1 keeper to `tests/sn/spatial/test_streaming_equilibrium_curvilinear.py`
   (which already has the matching cluster), then delete
   `test_phase_c_gates.py`. Discuss with user before deciding.
6. **G3f flip of `xfail strict` markers in `test_streaming_operator.py`
   + `test_l1_standoff_slab_cylinder.py`** — when these flip green,
   they MUST be removed in the same commit that lands G3f, not in a
   follow-up. Audit augmentation A flags this; the verification gate
   at G3f-1 (full L1 sweep) DETECTS the transition. **Decision**:
   the G3f commit message includes the xfail-removal as an explicit
   line item; the test changes are part of the G3f atomic commit.

---

## Summary table — verification artifact inventory

What this plan ships across G0/G1/G2:

| Artifact | Count | New / Existing | File path |
|---|---|---|---|
| L0 native matvec contract pins | 7 (× 1-3 geometries = ~15-20 cases) | NEW | `tests/sn/test_native_matvec.py` |
| L0/L1 G0 equivalence gate | 3 (slab/sphere/cylinder) | NEW | `tests/sn/test_g0_equivalence.py` |
| L0/L1 G1 contract pins | 6 | NEW | `tests/sn/test_fixed_source_g1.py` |
| L0/L1 G2 contract pins | 4 | NEW | `tests/sn/test_fixed_source_g2.py` |
| Existing L1 MMS gates preserved | 9 (across mms*.py) | EXISTING | tests/sn/test_mms*.py |
| Existing L1 kinf / homogeneous gates preserved | 6 | EXISTING | tests/sn/l1_analytical/, tests/sn/test_invertible_operator.py |
| Existing L1 standoff gates preserved + flipped | 4-leg L14 family on slab + cylinder | EXISTING | tests/sn/test_l1_standoff_slab_cylinder.py |
| Regression snapshot policy (no regen expected) | 11 | EXISTING | tests/sn/regression/snapshots/ |
| Tests retiring at G3a-G3g | ~22 (`test_snstreamingoperator.py`) + 10 (`test_phase_c_gates.py` migrations) + monkey-patch sites | EXISTING (delete) | listed in §"Tests retiring at G3f" |
| Tests with `catches=ERR-NNN` markers | Existing ERR-049 (3), ERR-050 (2), ERR-026 (3); + new ERR-049 (3), ERR-050 (1), ERR-026 (2) | mixed | various |

---

## Cross-references

* `.claude/plans/r1_step4_g_dependency_audit.md` — what retires.
* `.claude/plans/r1_step4_g_convention_crosswalk.md` — API contract.
* `.claude/plans/r1_step4_session2_followup.md` §2.3 + §3 + §4 —
  scope.
* `.claude/lessons.md` L17 (convention crosswalk), L18 (Pattern 7),
  L19 (None-default), L20 (retirement audit), L21 (sweep/matvec
  single strategy).
* `.claude/skills/vv-principles/SKILL.md` + `error_catalog.md`
  ERR-049 + ERR-050 + ERR-026 + ERR-006 + ERR-007.
* `.claude/skills/coding-elegance/SKILL.md` Pattern 2 (composition
  over duplication; the manifestation #6 / #7 cure) + Pattern 7
  (definition-site convention).
* `.claude/skills/numerical-bug-signatures/SKILL.md` Signature 1
  (curvilinear sweep divergence under refinement).
* GitHub issues: #199 (Step G omnibus — body needs §"Cross-cutting"
  table amendment), #174 (`_build_rhs_cartesian` — subsumed by G3a),
  #195 (curvilinear MMS pre-asymptotic magnitude — NOT closed by
  Step G), #200 (block-inverse face preconditioner — NOT closed by
  Step G), #160 (`SNStreamingOperator` — closes at G3f).
