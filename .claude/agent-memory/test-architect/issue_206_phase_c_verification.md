---
name: issue-206-phase-c-verification
description: Phase-C-specific verification plan for #206 — moving the 1-D SN matvec (_compute_LpC / _compute_LpC_transpose) OFF _MSpatialOperatorSum INTO _OneDimScanWalk as the apply-direction kernel. Per-leg verdict on the four-leg standoff (which legs are blind to the relocation + the fix), the bit-id-vs-nULP golden decision keyed on the density-vs-scan denom grouping choice, the decomposition-split twin gate, the transpose/adjoint + degenerate-cylinder gates, the Mode-9 convention-drift gate, and the per-sub-step (slab→sphere→cyl) gate-run recipe with route-arounds. Supersedes the PHASE C section of issue_206_cellupdate_seam_verification.md (line anchors re-verified live 2026-06-14 @ HEAD 03d5c02).
metadata:
  type: project
---

# Issue #206 Phase C — verification plan (test-architect, 2026-06-14)

Branch `refactor/sn-cellupdate-seam-slab` @ HEAD `03d5c02` (Phase B LANDED:
`_OneDimScanWalk` extracted, sweep relocated; matvec STILL on the operator).
This note is the Phase-C-specific deliverable; the broader 3-phase spec +
the L17 crosswalk live in [[issue-206-cellupdate-seam-verification]] (its line
anchors have drifted — the anchors below are re-verified LIVE this session).

## What Phase C does (verified against live code)

Move the 1-D matvec OFF `_MSpatialOperatorSum` (`orpheus/sn/operator.py`) INTO
`_OneDimScanWalk` (`orpheus/sn/loss_representation.py:1720`) as the
apply-direction kernel, mirroring the 2-D `_OctantWalk` carve (S6.3b `3a79ab3`).

LIVE anchors (re-grepped 2026-06-14, the audit's were stale):
- `_compute_LpC` — `operator.py:343` (the forward matvec; density form
  `m_full = (denom·psi_cell − numer_upstream)/V[i]` at **`:467`**, degenerate-cyl
  copy at **`:543`**; boundary block `m_boundary` **`:565-588`**).
- `_compute_LpC_transpose` — `operator.py:597` (the adjoint; the curvilinear
  `closure.angular_adjoint` reverse factor rides here).
- `_compute_decomposition` — `operator.py:774` (the **byte-twin dual-emission
  walk** of `_compute_LpC`; 3 LIVE production consumers — see §3).
- The current 1-D matvec twin = `CumprodScan.loss_action` (`loss_representation.py:730`)
  which TODAY forwards to `operator.M_spatial._compute_LpC(psi)` (`:741`); the
  transpose twin = `CumprodScan.loss_action_transpose` (`:743`) → `:755`. **Phase C
  re-targets these two methods to call into `_OneDimScanWalk`.**

## ⭐⭐ THE CRUX — the density-vs-scan denom grouping (the gate-tolerance decision)

The bit-identity boundary is precisely locatable in the live code:

- **Matvec (density form)** — `cell_balance_for_streaming`
  (`spatial/cell_balance.py:120`) builds `denom = (streaming_denom_term +
  angular_denom_term)[None,:] + collision_denom_term[:,None]` (`:247-250`,
  `streaming = 2·abs_mu·A_down`, `collision = total_xs·V`), then
  `_compute_LpC:467` does `(denom·psi_cell − numer_upstream)/V`. **`denom` is
  built and CONSUMED directly; no reciprocal.**
- **Sweep (scan form)** — `affine_scan_coefficients` (`spatial/diamond.py:379`)
  builds `denom = geometric_streaming_term[:,None,:] + collision_volume_term`
  (`:460`, SAME grouping order: streaming+curvature, then +collision) then
  `inverse_denom = 1.0/denom` (`:461`). The sweep body then does `b =
  2·(QV_chain + ang_contrib)·inverse_denom` (`loss_representation.py:2170`).
  **`denom` is reciprocated ONCE and cached; the scan multiplies by it.**

⟹ The two `denom` constructions use the **same FP grouping order**
(verified line-by-line). So Phase C has a genuine DESIGN CHOICE:

**Option DENSITY (RECOMMENDED for the gate)** — the relocated matvec kernel
reconstructs `denom` directly (keep consuming `cell_balance_for_streaming`'s
`denom`, or `denom = 1.0/inverse_denom` ONLY if it then immediately
`denom·ψ` — but the reciprocal round-trip `1/(1/denom)` is NOT bit-exact, so
prefer building `denom` from the same primitives). This stays **BIT-IDENTICAL**
to the current matvec → 0 ULP → the A-NEW matvec leg passes strict under
`-W error::DriftWarning`.

**Option SCAN** — the relocated matvec consumes the cached `inverse_denom` and
re-expresses `(denom·ψ − numer)/V` as the scan's `2·QV·inv_denom` grouping
(the elegance pull: ONE coefficient cache for sweep AND matvec). This costs a
**reciprocal-of-reciprocal per cell** → principled-equivalence at nULP, NOT
bit-identity. **nULP bound: ≤ nx ULP** (one extra reciprocal round-trip per
cell in the nx-cell WDD recurrence; the per-cell relative perturbation is ≤2
ULP and the recurrence is length nx, so the reduction-depth bound `reduction_depth
= nx` the A-NEW matvec leg already uses (`test_affine_carve_baseline.py:320`)
ABSORBS it). Structural-independence justification per `vv-principles`
§"Bit-identity vs principled-equivalence": (1) `inverse_denom`/`denom` are
NAMED, inspectable intermediates (the DD attenuation cache); (2) the value is
verified against the structurally-independent ground in leg 1
(`kinf_homogeneous` closed-form + trajectory_resolvent semi-analytical), NOT
against the old value alone; (3) the drift is FP-non-associativity bounded by
`reduction_depth × ULP = nx × ULP` — exactly the `kind="direct"` bound.

**GATE DECISION — state which the gate enforces:**
- The A-NEW matvec leg (`test_affine_carve_baseline.py::TestAffineCarveMatvecBaseline`,
  `kind="direct", reduction_depth=nx`) **already admits the Option-SCAN drift**
  (Layer 1 = `assert_array_almost_equal_nulp(nulp=nx)` passes ≤nx ULP) AND
  **surfaces it** (Layer 2 `DriftWarning` fires on ANY n_ulp>0).
- **If the main agent picks DENSITY:** run the A-NEW matvec leg under
  `-W error::tests.sn.regression._regression_assert.DriftWarning` in `sweep/core`
  (where the conftest does NOT override — FINDING 1) → STRICT bit-id, MUST be
  0 ULP. This is the recommended gate (cheapest, strongest).
- **If the main agent picks SCAN:** DROP the `-W error` escalation on the matvec
  leg ONLY; the bare `kind="direct", reduction_depth=nx` nULP gate is the
  contract. Document the relaxation in the test (it already has the docstring
  hook at `:55-63`: "principled-equivalence at nULP"). The SWEEP leg stays
  strict (it is untouched by the matvec relocation). DO NOT regenerate the
  golden — the nULP gate compares to the pre-carve `eab05ab`/`be4a57b` baseline,
  which is the structurally-independent inheritance anchor.

The recommendation is **DENSITY → strict bit-id**: it is the smaller blast
radius (the matvec is a pure relocation, not a re-expression), keeps the A-NEW
matvec leg as a sharp strict gate, and the elegance win of sharing the cache is
marginal (the matvec already shares `cell_balance_for_streaming` with `DD.residual`
via Pattern 2 — `diamond.py:267`). If the main agent insists on SCAN for
elegance, the nULP gate is sound but the strict A-NEW matvec arm is forfeited.

---

## 1. PER-LEG VERDICT on the four-leg standoff (L14)

| Leg | Test (LIVE path) | Pins relocation? | Tol | Blind-spot / fix |
|-----|------------------|------------------|-----|------------------|
| 1 (sweep ≡ indep ref) | `tests/sn/verification/analytical/test_phase_c_crosscheck.py` | **SWEEP only — NOT the matvec** | value (kinf exact; rel<2% MR) | ⚠ GAP: leg 1 runs `solve_sn` (the SWEEP path); it NEVER calls `(L+C).apply`. A matvec-only relocation bug (Krylov path) is invisible here. FIX: leg 1 stays the leg-1 ground for the SWEEP/value; the MATVEC value-ground is leg 2 + the Krylov-vs-SI consistency in `test_sweep_vs_apply_consistency.py` (which DOES drive the matvec via Krylov). |
| 2 (matvec ≡ ref) | `test_bc_extraction_matvec.py::TestStreamingEquilibriumValue::test_flat_flux_per_ordinate_balance_no_pole_spike` + `TestVacuumMatvecBitIdentity` | **partial** | bit-id (vacuum) / spike-factor 2× (flat) | ⚠ **DOUBLE GAP**: (a) the spike gate is FLAT-ψ → vv §H2 NULLS curvilinear redistribution → blind to a `dA_w·c_out` relocation bug; AND (b) it is **BARE assert** (`:477`) → INERT under -O (Mode-8). FIX: the A-NEW matvec leg (het σ_t + NON-flat ψ, `np.testing`-based, -O-safe) is the real leg-2 matvec ground; the spike gate is a complementary L0 anchor that MUST run WITHOUT -O. **The transpose dense-probe oracle the prior note cited (`diag_p42_adjoint_oracle.py`) NO LONGER EXISTS** — leg-2-transpose is now ONLY the reciprocity gate (see §4). |
| 3 (sweep ≡ matvec TWIN) | `test_loss_action_convention.py` + `test_sweep_vs_apply_consistency.py` + `test_streaming_operator_decomposition.py` | **YES — the load-bearing relocation leg** | mixed (see below) | The non-tautological flat-reflective anchor (`test_loss_action_is_full_loss_LpC_flat_reflective`, -O-safe `np.testing`) + the ≥2G het `apply == loss_action − C·ψ` (`test_apply_equals_loss_action_minus_independent_collision_het`, -O-safe) DIRECTLY pin "the relocated matvec is still (L+C)". ⚠ `test_streaming_operator_decomposition.py::TestResolutionADecomposition::test_bit_exact_uniform_sigma_t` (`:180/:197` bare `assert rel<1e-14`) is INERT under -O; its sibling `TestSubtractiveDefinition` (`:257` `assert_array_equal`) IS -O-safe. |
| 4 (refinement) | `test_phase_c_crosscheck.py` flux-shape/MR rows | SWEEP only | value | Convergence RATE necessary not sufficient (vv §5); leg 1 supplies the converged VALUE ground. Same matvec-blindness as leg 1. |

**Headline per-leg finding:** **legs 1 + 4 are GREEN-but-blind to the matvec
relocation** (they exercise the SWEEP path only — `solve_sn`/SI). The matvec
relocation is pinned by **leg 3** (the twin gate — directly `(L+C).apply` ≡
`loss_action − C`) and **leg 2** (matvec ≡ reference value). The Krylov path
(which drives the matvec end-to-end) is pinned by `test_sweep_vs_apply_consistency.py`
(`test_solve_sn_si_vs_krylov_consistency_homogeneous_sphere` etc. — SI sweep ≡
Krylov matvec keff). **Do NOT rely on leg 1/4 to catch a matvec-only sign/index
drift.**

---

## 2. THE MATVEC-RELOCATION GOLDEN DECISION

**The existing A-NEW gate SUFFICES; no Phase-C-specific golden is needed** —
BUT the runtime mode + tolerance per leg must be set explicitly.

`tests/sn/sweep/core/test_affine_carve_baseline.py` already captures (at the
pre-carve `eab05ab`/`be4a57b` baseline, committed) BOTH:
- **Sweep leg** (`TestAffineCarveSweepBaseline`, `:221`) — one `transport_sweep`
  on het σ_t + non-flat random Q, ≥2G, slab/sphere/cyl; `reduction_depth=nx`
  (angular) / `nx+N` (scalar). UNTOUCHED by the matvec relocation → stays STRICT.
- **Matvec leg** (`TestAffineCarveMatvecBaseline`, `:283`) — one `(L+C).apply(ψ).bulk`
  on `het_operands` (het σ_t + non-flat random ψ + random boundary trace, every
  redistribution term ACTIVE), ≥2G, slab/sphere/cyl; `reduction_depth=nx`.

Per-leg tolerance for Phase C:
- **Sweep leg:** STRICT bit-id (`-W error::...DriftWarning` in sweep/core).
  MUST be 0 ULP (the sweep is not touched by Phase C).
- **Matvec leg:** DENSITY option → STRICT bit-id (0 ULP, `-W error`). SCAN
  option → nULP, bound `≤ nx ULP` (the `kind="direct"` bound at `reduction_depth=nx`),
  drop the `-W error` escalation on this leg only, keep the bare nULP assert.

Structural-independence justification (the nULP path, if chosen): the matvec
VALUE is grounded in leg 1/leg 2's structurally-independent references
(`kinf_homogeneous` closed-form, trajectory_resolvent semi-analytical, Q/Σ_t
closed-form), NOT in the old matvec value. Old-vs-new ULP is the inheritance
inertness check; the references are the correctness ground. All three
`vv-principles` criteria hold (named intermediate, independent reference,
FP-non-associativity bounded by reduction depth).

---

## 3. THE DECOMPOSITION-SPLIT GATE (the twin Phase C must NOT break)

⚠ **CRITICAL ARCHITECTURAL HAZARD specific to Phase C.** `_compute_decomposition`
(`operator.py:774`) is the **byte-twin dual-emission walk** of `_compute_LpC`
(same cell loop, same closure routing, but emits `(M_spatial, M_angular_redist)`
into two buffers instead of the fused `(L+C)ψ`). It has **3 LIVE production
consumers** (re-verified):
- `operator.py:269` — transient orchestrator (`_SpatialSweepDirection.apply` →
  `m_spat, m_ang`).
- `operator.py:1158` — `_MSpatialOperatorSum.apply` → `m_spat, _`.
- `operator.py:1276` — `AngularRedistributionOperator.apply` → `_, m_ang`.

**If Phase C moves `_compute_LpC` into `_OneDimScanWalk` but leaves
`_compute_decomposition` on the operator, the two walks DIVERGE — re-introducing
the exact twin-path the carve exists to kill (Pattern 2 violation, the Phase-F
shape).** The carve MUST either (a) also route `_compute_decomposition`'s walk
through the SAME `_OneDimScanWalk` frame (the right answer — single source), or
(b) explicitly document why the dual-emission split stays a separate frame and
pin that the two STILL agree.

**Gate that proves the split stays correct through the carve:**
`tests/sn/operators/test_streaming_operator.py::TestT4cAlgebraDecompositionInvariantCurvilinear`
(`:854`) — `test_sphere_LpC_equals_M_spatial_plus_M_angular_redist` (`:882`) +
`test_cylinder_...` (`:908`) + slab `TestT4bAlgebraDecompositionInvariantSlab`
(`:714`). These assert `(L+C)ψ ≡ M_spatial·ψ + M_angular_redist·ψ` PER-ORDINATE-BULK
via `np.testing.assert_allclose`/`assert_array_equal` (**-O-SAFE** — verified).
**This is the gate that catches a divergent decomposition twin.** MUST stay
green AND MUST be re-run after the carve to confirm the 3 consumers still get a
correct split. Add a NEW assertion if Phase C unifies the frames: that
`_compute_decomposition`'s two outputs SUM to `_compute_LpC`'s one output
bit-identically (a single-source tripwire) — `M_spatial.apply + M_angular_redist.apply
== (L+C).apply` at `array_equal` (this is the Pattern-2 single-source pin; if the
frames are unified it is bit-id by construction, if not it is the divergence
catcher).

---

## 4. THE TRANSPOSE/ADJOINT GATE + the degenerate-cylinder gate

### Adjoint (`_compute_LpC_transpose`, `operator.py:597`)

**The ONLY -O-firing gate that exercises the curvilinear angular adjoint
(`closure.angular_adjoint`):**
`tests/sn/operators/test_g_adjoint_reciprocity.py::test_g_adjoint_reciprocity_full_block`
(`:199`) — `⟨Aψ,φ⟩_G = ⟨ψ, A.Hφ⟩_G` for `A=L+C` (reflective BC so the A_ss
block is live), random NON-flat ψ/φ on bulk+trace, over
`slab/slab_2g/sphere_2g/cyl` (the `_BUILDERS` dict `:122-126`), all `pytest.fail`
(**-O-SAFE**). PAIRED with the L11 negative control
`test_wrong_trace_metric_breaks_reciprocity` (`:267`, slab+sphere) — dropping
|Ω·n| MUST break reciprocity (proves §2 is metric-aware not blind, per ERR-051).
**MUST stay green through Phase C** — it walks `_compute_LpC_transpose` via
`A.H.apply`. This is the load-bearing adjoint gate.

⚠ **GAP — the standalone transpose dense-probe oracle is GONE.** The prior note
cited `derivations/diagnostics/diag_p42_adjoint_oracle.py` as the leg-2-transpose
ground; it NO LONGER EXISTS (verified absent). The reciprocity gate is now the
SOLE adjoint cover. Reciprocity is a STRONG structural gate (it pins
`Aᵀ` against `A` via the metric, structurally independent of the forward
implementation) — but it is a CONSISTENCY check, not a value-against-reference
check. **Recommendation:** this is sufficient (reciprocity + L11 control is the
canonical adjoint verification per `test_g_adjoint_reciprocity`'s own design); a
dense-probe `Aᵀ` oracle would be additive but is NOT required. If the main agent
wants belt-and-suspenders, a `kind="direct"` snapshot of `(L+C).apply_transpose(φ).bulk`
on `het_operands` φ (mirror of the A-NEW matvec leg) is the cheapest add — it
inherits bit-id/nULP from the same baseline mechanism.

### Degenerate-cylinder branch (THE HARD PART — no scan/slab analogue)

The cylinder pure-azimuthal ordinate (`|η| < 1e-15`) has NO `ordinate_scan`
analogue: the matvec handles it as a direct per-cell `(denom·ψ − numer)/V` with
`A_downstream=0.0`, `psi_face_in=zero_face` (`operator.py:504-544`); the sweep
handles it via the SLOW per-cell `cell_update.update` path (`loss_representation.py:2133-2157`).
These are **structurally DIFFERENT branches in the two paths** — Phase C moving
the matvec to `_OneDimScanWalk` MUST preserve the matvec's degenerate branch
(it CANNOT route through the scan). Pins:
- `test_g_adjoint_reciprocity[cyl]` (drives the degenerate branch via `A.apply`/`A.H.apply`).
- The A-NEW matvec leg `[CYL]` (`test_affine_carve_baseline.py`, het σ_t + non-flat ψ
  — the degenerate ordinate is in the quad, so the branch is exercised; bit-id/nULP).
- The decomposition invariant `test_cylinder_LpC_equals_M_spatial_plus_M_angular_redist`
  (`:908`).
- ⚠ The cylinder build in the A-NEW gate uses `level_symmetric` cubature
  (`_build_sn_mesh:134` — "the curvilinear cubature the matvec serves") which
  POPULATES degenerate ordinates; confirm the relocated matvec's degenerate
  branch fires under this quad (a `level_symmetric` cyl has the `|η|<eps`
  ordinate). **If the degenerate branch is mishandled, the A-NEW matvec[CYL] leg
  goes red (NaN or wrong value at the degenerate column).** This is the sharpest
  degenerate-branch catcher.

---

## 5. MODE-9 CONVENTION-DRIFT GATE (per-ordinate, het, aniso, ≥2G, non-flat)

The matvec touches the curvilinear M-M angular closure + the Carlson
coupled-pole seed (ERR-058). A wrong angular relocation is INVISIBLE on the flat
box (vv §H2 nulls redistribution; vv §H3 per-ordinate balance telescopes). The
gates that run on STRESSING configs (NOT the degenerate box):

- **Existing, -O-safe, sufficient for the spatial+collision matvec:**
  `tests/sn/solve/test_affine_carve_bit_identity.py` — `si_slab_2g_het` (P1
  anisotropic, ≥2G, het σ_t, non-flat; sha256, `raise AssertionError` so -O-safe)
  + the 2-D row. The slab arm is the end-to-end ≥2G-het-aniso anchor.
- **Existing, -O-safe, the per-ordinate twin:** `test_loss_action_convention.py::
  test_apply_equals_loss_action_minus_independent_collision_het` (≥2G het non-flat,
  `np.testing`).
- **GAP — no SPHERE/CYL anisotropic per-ordinate matvec golden.** The A-NEW
  matvec leg uses `het_operands` (het σ_t + non-flat ψ — activates `dA_w·c_out`
  redistribution) but **isotropic-only** (no σ_s1 / P1 source in the matvec
  probe; the matvec is `(L+C)` which has no scattering — so P1 enters only via
  the *moments of ψ*, and `het_operands` ψ is non-flat-in-μ so the M-M
  redistribution IS exercised). Verdict: the curvilinear M-M redistribution path
  IS covered by A-NEW[SPH/CYL] + reciprocity[SPH/CYL] on non-flat ψ — the
  redistribution `dA_w·c_out` term is ACTIVE on any non-flat-in-μ ψ regardless
  of scattering order. **A dedicated sphere-aniso companion is NOT required for
  the matvec** (the matvec is `(L+C)`, scattering-free; the aniso angular
  coupling that would need a P1 companion lives in the SCATTERING operator S, not
  in this carve). The Mode-2 (μ swap) / Mode-4 (recurrence drift) catchers are
  the per-ordinate (NOT weight-summed) A-NEW matvec leg + reciprocity — both
  per-ordinate by construction. **CONFIRM the A-NEW matvec leg asserts on
  `out.bulk.values` (the full `(N,ng,nx)` per-ordinate array — it does,
  `:319`), NEVER on weight-summed scalar flux.** A μ_x sign / mirror-index drift
  preserves Σw·ψ but corrupts per-ordinate ψ (Mode-2) — the per-ordinate array
  assertion is the catcher.

**Carlson coupled-pole seed (ERR-058) — the invisible-on-flat-ψ hazard.** The
matvec reads `outflow_at_inner.T[quad.reflection_index("x")]` (`operator.py:495`);
the sweep reads `pole_outflow[mirror[global_n]]` (`loss_representation.py:2130`).
Both consume the SAME mirror permutation + the inward ordinate's pole-face
outflow (Pattern 2). Phase C MUST keep the matvec's seed routing.
Catchers: `test_g_adjoint_reciprocity[sphere_2g]` (non-flat) + the sphere
`test_phase_c_crosscheck` flux-shape rows + the A-NEW matvec[SPH] leg (non-flat ψ
drives the O(h)-wrong pole-cell-centre read out of cancellation). The
historical pole-CELL-centre read was exact on flat ψ — so a flat-ψ gate would
pass while this hides (vv §H2).

---

## 6. PER-SUB-STEP GATE-RUN RECIPE (slab → sphere → cylinder)

Env: HOST → `.venv/bin/python`. Canonical: `python -O -m pytest`. The carve
proceeds geometry-by-geometry; run the GEOMETRY-SCOPED gate at each sub-step,
then the full set at the end. Route-arounds baked in (FINDING + pre-existing reds).

### Standing route-arounds (apply to EVERY command)
- `-k "not (vacuum_bulk_bit_identical_1d and SPH)"` — 3 stale post-ERR-058 SPH
  vacuum-matvec snapshots (CONFIRMED red at clean HEAD this session).
- `-k "not (sphere_1g_apply_bit_identical or sphere_2g_apply_bit_identical)"` —
  2 stale post-ERR-058 sphere apply snapshots (CONFIRMED red this session; fire
  under -O via `assert_allclose`).
- `--deselect tests/sn/eigenvalue/test_keff_slab.py::test_heterogeneous_absolute_keff` —
  #212 `continuous_get` hang (orthogonal).
- `-p no:cacheprovider` (clean collection).

### THE STRICT BIT-ID GATE (where `-W error` actually escalates — FINDING 1)
`-W error::DriftWarning` is INERT in `tests/sn/regression/` (its conftest forces
`always`). It ONLY escalates in `tests/sn/sweep/core/` + `tests/sn/solve/`. Use
the QUALIFIED path (bare `error::DriftWarning` → AttributeError):
`-W "error::tests.sn.regression._regression_assert.DriftWarning"`.

**STRICT gate (DENSITY option — recommended):**
```
.venv/bin/python -O -m pytest \
  tests/sn/sweep/core/test_affine_carve_baseline.py \
  tests/sn/solve/test_affine_carve_bit_identity.py \
  -W "error::tests.sn.regression._regression_assert.DriftWarning" \
  -p no:cacheprovider -q
```
Per sub-step scope: `-k "SLB or slab"` (slab), `-k SPH` (sphere), `-k CYL`
(cylinder). GREEN = the geometry's sweep+matvec legs 0-ULP-identical to the
pre-carve baseline. (SCAN option: drop `-W error` on the matvec leg, keep the
bare `kind="direct" nulp=nx` — see §2.)

### THE VALUE/TOLERANCE GATE (the four-leg standoff + decomposition + adjoint)
```
.venv/bin/python -O -m pytest \
  tests/sn/operators/test_bc_extraction_matvec.py \
  tests/sn/operators/test_loss_action_convention.py \
  tests/sn/operators/test_streaming_operator_decomposition.py \
  tests/sn/operators/test_g_adjoint_reciprocity.py \
  tests/sn/sweep/core/test_sweep_vs_apply_consistency.py \
  tests/sn/solve/test_affine_carve_bit_identity.py \
  tests/sn/operators/test_streaming_operator.py \
  -k "not (vacuum_bulk_bit_identical_1d and SPH) and not (sphere_1g_apply_bit_identical or sphere_2g_apply_bit_identical)" \
  -p no:cacheprovider -q
```
**Live baseline this session (HEAD 03d5c02, fast subset minus crosscheck):
122 passed / 3 skipped / 3 deselected.** Adding `test_streaming_operator.py`
brings the decomposition invariant + capability gates (the 2 sphere-apply reds
are -k-excluded). GREEN = legs 2/3/4(value) + decomposition + adjoint hold.

### THE LEG-1 STRUCTURALLY-INDEPENDENT GROUND (slow — run once per geometry milestone)
```
.venv/bin/python -O -m pytest \
  tests/sn/verification/analytical/test_phase_c_crosscheck.py \
  -p no:cacheprovider -q
```
SLOW (trajectory_resolvent cylinder-MR GL quadrature is heavy — minutes). The
FAST leg-1 closed-form arm is `-k test_sn_spherical_homogeneous_kinf_recovery_2g`
(verified 1 passed, 0.5s this session). Run the full file at the slab→sphere→cyl
milestones, not every micro-commit.

### THE Mode-8 BARE-ASSERT GATE (run WITHOUT -O — these are INERT under -O)
The following Phase-C-load-bearing checks use BARE assert (stripped by -O,
verified via the pytest "assertions not in test modules" warning + the explicit
`^\s+assert ` audit this session):
- `test_bc_extraction_matvec.py:477` — `test_flat_flux_per_ordinate_balance_no_pole_spike`
  (the Mode-3 missing-ΔA/w pole-spike L0 catcher) — **bare assert, INERT under -O**.
- `test_streaming_operator_decomposition.py:180/197` — `TestResolutionADecomposition::
  test_bit_exact_uniform_sigma_t` (leg-3 bit-exact twin) — **bare, INERT**.
  (sibling `TestSubtractiveDefinition:257` IS -O-safe.)
- `test_streaming_operator.py` — 36 bare-assert lines (capability,
  `assert isinstance`, decomposition structural). The `assert_*` rows fire under
  -O; the bare-`assert` rows do not.
- `test_phase_c_gates.py:442/446` (Gate 1.2 bit-id) + `:396` (Gate 1.3
  reciprocity) — bare assert. **NOTE: `test_phase_c_gates.py` is Issue #168's
  Phase C (the 2-D sweep-frame apply rewrite), NOT #206's 1-D matvec relocation
  — it is adjacent coverage, not a #206 gate. Run it but do not treat its name
  as a #206 gate.**
```
.venv/bin/python -m pytest \
  tests/sn/operators/test_bc_extraction_matvec.py \
  tests/sn/operators/test_streaming_operator_decomposition.py \
  tests/sn/operators/test_streaming_operator.py \
  -k "not (vacuum_bulk_bit_identical_1d and SPH) and not (sphere_1g_apply_bit_identical or sphere_2g_apply_bit_identical)" \
  -p no:cacheprovider -q
```
(NO `-O`.) GREEN here = the pole-spike Mode-3 catcher + the bit-exact leg-3 twin
+ the decomposition structural asserts actually FIRED. **Mode-8 fix-at-touch
recommendation:** at the carve commit, migrate `test_flat_flux_per_ordinate_balance_no_pole_spike`'s
`:477` bare `assert` → `np.testing.assert_array_less` / `pytest.fail` so the
Mode-3 catcher fires under the canonical -O. (It is the sharpest curvilinear
missing-ΔA/w gate; leaving it -O-inert is a false-green hazard per vv Mode 8.)

### "GREEN means" per sub-step
- **Slab sub-step:** STRICT gate `-k "SLB or slab"` 0-ULP; VALUE gate slab rows
  green; `si_slab_2g_het` sha256 green. (Slab has no curvilinear redistribution —
  the simplest milestone; the density form is trivially bit-id.)
- **Sphere sub-step:** STRICT gate `-k SPH` 0-ULP (DENSITY) or ≤nx ULP (SCAN);
  `test_g_adjoint_reciprocity[sphere_2g]` + `test_sphere_LpC_equals_M_spatial_plus_M_angular_redist`
  green; A-NEW matvec[SPH] green; pole-spike[SPH] green (run WITHOUT -O). The
  3 stale-SPH reds REMAIN red (route around — they are not ours).
- **Cylinder sub-step:** STRICT gate `-k CYL` 0-ULP/≤nx ULP; A-NEW matvec[CYL]
  green (the degenerate-ordinate branch fires under `level_symmetric` cubature —
  the sharpest degenerate catcher); `test_g_adjoint_reciprocity[cyl]` +
  `test_cylinder_LpC_equals_M_spatial_plus_M_angular_redist` green;
  `test_sweep_vs_apply_consistency` cylinder rows green. #206 cylinder-matvec
  history (#206 itself) — confirm no NEW divergence vs the baseline.

---

## Self-improvement note (no new failure mode)

Phase C is the canonical Pattern-2 twin-path unification (same shape as
Phase F / S6.3b). The ONE genuinely Phase-C-specific architectural hazard is
**§3 — the `_compute_decomposition` dual-emission twin that must NOT diverge
when `_compute_LpC` relocates**. This is NOT a new vv failure mode (it is
Pattern-2 / anti-pattern #1 "twin path", already catalogued); the gate is the
decomposition invariant (`test_streaming_operator.py::TestT4c...`) plus the
recommended single-source tripwire (`M_spatial + M_angular_redist == (L+C)` at
`array_equal`). The Mode-8 inertness of the pole-spike Mode-3 catcher
(`:477` bare assert) is the fix-at-touch item. No skill-table append needed.
