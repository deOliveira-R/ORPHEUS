---
name: issue-240-phase2-step-b-sigma-leak
description: Verification spec for #240 Phase 2 Step B — the "Σ matvec leak" (InvertibleOperator inherits OperatorSum.apply = L.apply+C.apply instead of owning loss_action(σ_C)). ⭐ THE PREMISE CORRECTION established this session by probe — the leak is NOT a numerical bug (M(σ)ψ is AFFINE in σ in the forward direction → leaky apply ≡ M(σ_r)ψ to ≤2 ULP, the round-trip op.apply(op.solve(q))≈q PASSES under the leak). It is a principled-equivalence/elegance carve. Gate file + xfail-teeth design.
metadata:
  type: project
---

# #240 Phase 2 Step B verification spec — the Σ matvec leak

**Status:** PRE-IMPLEMENTATION. Branch `feature/sn-space-angle-tier2`,
HEAD `28d76c9` (Step A `d717d4d` LANDED). Host env, canonical `python -O`.
Spec stub file WRITTEN + verified GREEN: `tests/sn/operators/test_removal_form_matvec_sweep.py`
(12 passed / 8 xfailed under `-O`, fully Mode-8-clean).

## ⭐⭐ THE PREMISE CORRECTION (the headline — established by probe this session)

The driving brief framed the leak as a NUMERICAL defect: "matvec ≠ sweep⁻¹
for σ_r ≠ σ_t — silent L21 violation". **This is FALSE as a numerical claim.**

Probe (`/tmp/probe_affine_robust.py`, `/tmp/probe_roundtrip_teeth.py`) across
slab/sphere/cyl/2D-Cartesian × vacuum/reflective × ≥2G het non-flat ψ:

- `leaky_apply(σ_t,σ_r) ≡ M(σ_r)ψ` to **≤2 ULP (rel ≤2.5e-16)** — EVERY geom/BC.
- `op.apply(op.solve(q)) ≈ q` and `op.solve(op.apply(ψ)) ≈ ψ` **PASS at ≤1.3e-15
  under the CURRENT leaky code** (slab/cyl/2D).

**Mechanism:** in the FORWARD (apply) direction the WDD matvec is
`M(σ)ψ = (denom·ψ_cell − numer_upstream)/V` with `denom = streaming_term + σ·V`
⟹ `M(σ)ψ = streaming_action(ψ) + σ·ψ_cell` — **AFFINE in σ**. So the inherited
`OperatorSum.apply` leak `L.apply(σ_t)+C.apply(σ_r) = (M(σ_t)ψ−σ_t·ψ)+σ_r·ψ =
streaming_action(ψ)+σ_r·ψ = M(σ_r)ψ`. The non-affinity of M in σ (σ in `1/denom`
of the cell-AVERAGE) lives in the SWEEP/SOLVE direction, NOT apply. The
forward action never inverts the denominator.

⟹ **The proposed round-trip gate has NO TEETH** against the leak (cannot
distinguish leaky from override — a vv Mode-9 degeneracy IN THE GATE DESIGN
ITSELF). The "negative control demonstrating the gate fails under pre-fix
leaky code" CANNOT be a numerical round-trip — the pre-fix code is value-correct.

**What the carve actually fixes:** a principled-equivalence/ELEGANCE defect
(coding-elegance Pattern 1+2): the composite `(L+C)` OWNS a single matvec
`loss_action(σ_C)`; realising it as `L.apply+C.apply` re-derives the streaming
action through σ_t then cancels it (redundant round-trip, twin-path realisation
of ONE operator's action). It becomes a genuine LATENT trap only if a future
refactor makes `L.apply` non-affine in σ. NOT a new ERR-NNN (no wrong value
ever shipped). Verification problem = the #240-Phase-2-Step-A shape
(value-preserving re-association): same-value nULP + structural-independence +
the teeth must be a STRUCTURAL gate.

## The fix being verified

1. `loss_action(operator, psi)` → `loss_action(sigma, psi)` (3 bodies read
   `operator.sigma_t`: `loss_representation.py:626` `_OctantWalk`, `:1839`
   `_OneDimScanWalk._apply_walk` arg, `:2288` transpose; forwarders 749/763/
   933/1132/1384/1460 thread `operator`→`sigma`). Symmetric with `sweep(Q,sig_t,…)`.
2. `StreamingOperator.apply/apply_transpose` (`operator.py:805/856`) pass
   `self.sigma_t` (leaf keeps own σ; Resolution-A `−self.sigma_t` unchanged).
3. NEW `InvertibleOperator.apply/apply_transpose` overrides:
   `return self.loss_representation.loss_action(self.sigma, psi)` (the FULL
   `M(σ_C)ψ` directly, no leaf decomposition). InvertibleOperator (`operator.py:1130`)
   currently has NO apply override → inherits `OperatorSum.apply` (`numerics/
   operator.py:740` = `a.apply(x)+b.apply(x)`). `CAP_APPLY_TRANSPOSE` already
   advertised (both L,C have it) → override keeps it; `.H` wiring transitive.

## Probe-established value facts (production σ_C==σ_t, override vs leaf-sum)

- slab (`_OneDimScanWalk`): **BIT-IDENTICAL 0/32** → `array_equal`.
- sphere (`_OneDimScanWalk`): **BIT-IDENTICAL 0/32** → `array_equal`.
- cylinder: 2/192, ≤1 ULP → nULP. 2-D-Cartesian (`_OctantWalk`): 26/768, ≤2 ULP → nULP.
- boundary block: bit-identical ALL arms.

## ⭐ SPHERE ROUND-TRIP DIVERGES (the geometry-fragility finding)

Probe `/tmp/probe_sphere_rt.py`: sphere bare-operator `op.apply(op.solve(q))`
DIVERGES (9.8e5 single, 2e99 iterated-40) — the curvilinear M-M sweep reads the
PREVIOUS ITERATE for the Carlson coupled-pole closure (ERR-058), so a single
`op.solve` is NOT the one-shot inverse of `op.apply`. **Cylinder round-trips
CLEANLY (3.1e-15)** — its degenerate pure-azimuthal ordinate does not break the
one-shot inverse. ⟹ the round-trip gate is slab/cyl/2D ONLY; sphere's matvec≡
sweep claim is carried by the STRUCTURAL teeth gate (a) `apply==M(σ_r)` (does
NOT round-trip) + the existing `TestInvertibleSolveBridgeRegression` fixed-point
bridge (production σ).

## THE GATE FILE: `tests/sn/operators/test_removal_form_matvec_sweep.py`

`@pytest.mark.foundation` + `verifies("loss-rep-resolution-a")`. NO `catches`.
4 groups:

- **(a) STRUCTURAL TEETH** (the ONLY discriminator) — `xfail(strict=True,
  raises=(AssertionError,AttributeError))` until override:
  - `test_invertible_apply_is_M_of_C_sigma_bit_identical` (slab/sph/cyl/2D):
    `op.apply.bulk == M(σ_r)ψ` via `array_equal`, where `M(σ_r)` = SEPARATE
    `StreamingOperator(σ_r).loss_action(L_ref,ψ)` (leaf whose OWN σ_t IS σ_r —
    unambiguous; built with CURRENT signature so reference doesn't error pre-fix).
    Override→bit-id (xpass); leak→leaf-sum ≤2 ULP off→AssertionError (teeth fire).
  - `test_invertible_apply_transpose_is_M_transpose_of_C_sigma_bit_identical`
    (slab/sph/cyl — **NOT cart2d**: `ScanMarch.loss_action_transpose` is a
    deferral raise, 2-D adjoint not impl).
- **(b) VALUE GROUND** (pass under leak AND override — NOT teeth):
  - `test_removal_form_matvec_sweep_roundtrip` (slab/cyl/2D — sphere EXCLUDED):
    `op.apply(op.solve(q))≈q` + `op.solve(op.apply(ψ))≈ψ` at rtol 1e-10.
  - `test_removal_form_apply_value_equals_M_of_sigma_r` (all 4): `op.apply ==
    M(σ_r)` at atol 1e-12 (affine-in-σ characterisation; the SI structural ref).
  - `test_removal_form_kinf_independent_reference_2g`: `pytest.xfail` — needs
    #200 solver entry (fold within-grp self-scatter into diagonal). Closed-form
    eigenvalue ground lands with #200; operator-level round-trip is Step-B's ground.
- **(c) PRODUCTION-σ INVARIANT** `test_production_sigma_apply_value_preserved`
  (σ_C==σ_t): override == leaf-sum; slab/sph `array_equal` (True), cyl/2D
  `nulp=4`+`pytest.fail` drift bound 1e-14 (Mode-8). Pins value-preservation.
- **(d) NEGATIVE CONTROL** `test_removal_form_nonpositive_sigma_r_rejected`:
  `InvertibleOperator(L,C(σ_r≤0))` raises "strictly positive" (`operator.py:1244`).
  The vv-L11 positive(removal cases)+negative(here) pair for the σ_r>0 contract.

σ_r built as numpy array directly (`σ_t − Σ_s0`, bounded draw → >0) — independent
of mesh.materials (which carry no scatter) → removal form representable on EVERY
geometry. 2-D = NON-SQUARE 4×5 `level_symmetric` (genuine mu_y, #214-safe; x↔y catcher).

## EXISTING-TEST RULING TABLE (each ruled live this session)

| Test | file:line | call | Verdict |
|------|-----------|------|---------|
| `test_loss_action_is_full_loss_LpC_flat_reflective` | `test_loss_action_convention.py:110` | `rep.loss_action(L,psi)` | MIGRATE → `loss_action(sig_t,psi)` |
| `test_apply_equals_loss_action_minus_independent_collision_het` | `:151` | `rep.loss_action(L,psi)` | MIGRATE → `loss_action(sig_t,psi)` |
| `test_apply_vs_sweep_2d_residual_cancellation` (ERR-026) | `test_2d_l2_matvec_correctness.py:174` | `A.apply==L.apply+C.apply` atol1e-12 | **KEEP-as-allclose** (≤2 ULP « 1e-12, survives). It is FALSE-by-construction for σ_r≠σ_t → NO LONGER the load-bearing ERR-026 gate. RE-POINT the ERR-026 anchor to the new teeth gate (a) `apply==M(σ_r)` + keep this as the σ_t consistency check. Update docstring. |
| `test_apply_equals_l_plus_c_on_typed_flux` | `test_invertible_operator.py:259` | slab `array_equal` | **KEEP** (slab bit-id 0/32 survives strict). Stale docstring claim "matches inherited OperatorSum" → update: now matches loss_action(σ_C), bit-id at σ_C==σ_t on slab. |
| `test_both_matvec_variants_share_the_walk` | `test_one_octant_walk.py:149` | `rep_cls(sn).loss_action(L,psi)` | MIGRATE → `loss_action(sig_t,psi)` |
| scanmarch/full-field oracles | `test_scan_march_equivalence.py:205-206`, `test_2d_full_field_oracle.py:136-137` | `<rep>(sn).loss_action(L,state)` | MIGRATE → `loss_action(sig_t,state)` (sig_t in scope) |
| `test_dag_free_representations_never_reference_the_substrate` | `test_dag_ownership.py:130` | NAME-strings `"loss_action"` + `inspect.getsource` for `"sweep_graphs"` | **KEEP unchanged** (no signature assertion). |
| `test_g_adjoint_reciprocity_full_block` | `test_g_adjoint_reciprocity.py:199` | `A=L+C−B`, `A.H.apply(phi)` | **KEEP green** (σ_C==σ_t → ≤2 ULP, rel<1e-12 survives). The `(L+C)` sub-term routes through NEW `InvertibleOperator.apply_transpose` via `.H` (OperatorDifference distributes); transitively exercises the override. Capability `apply_transpose` stays advertised. |
| `_LC_matvec` helper | `_test_helpers.py:343` | `(L+C).apply(psi)` | becomes `loss_action(σ_C)` under override. `test_streaming_operator_decomposition.py:167/180` compares `_LC_matvec(σ_t)` vs `L.apply+C.apply` at `rel<1e-14` (NOT array_equal) → survives ≤2 ULP. ⚠ bare-assert `:180/:197` INERT under -O (Mode-8 pre-existing). |
| `test_b1pp_verification.py:137/232/312` | free expr `L.apply+C.apply` | **KEEP unchanged** (does NOT route through override; full-rank claim, σ_C==σ_t). |
| `test_LC_apply_composite` | `test_streaming_operator.py:494` | `A.apply` type/shape only | KEEP (stale docstring). |

## RE-BASELINE / STRICT GATE

The strict DriftWarning gate: `python -O -m pytest tests/sn/sweep/core
tests/sn/solve -W "error::tests.sn.regression._regression_assert.DriftWarning"`
(the `tests.sn.regression...` path — `orpheus.sn...` is WRONG, silently fails to
escalate; bare `error::DriftWarning`→AttributeError). DD SWEEP/SOLVE snapshots
STAY STRICT (apply re-association doesn't touch solve). APPLY snapshots that pin
`(L+C).apply` (slab/cart2d via `TestT4bPreT4RegressionSnapshot`) → already
`assert_regression(kind=direct,reduction_depth=nx)` (Phase-1) → ≤2 ULP absorbed.
NOTE: the production path is σ_C==σ_t where slab/sphere are BIT-ID, cyl/2D ≤2 ULP
→ existing apply snapshots that touch cyl/2D may need the kind=direct nULP if
they were strict; slab stays strict. Confirm at impl.

## ROUTE-AROUNDS (7 pre-existing reds; NEVER all tests/sn — #212 hang)
`-k "not (vacuum_bulk_bit_identical_1d and SPH) and not (sphere_1g_apply_bit_identical
or sphere_2g_apply_bit_identical) and not test_2d_mesh_resolution and not
two_d_cartesian_loss_action"`. Green floor: tests/sn/operators+spatial = 507p/4skip
(+ the 2 #214 reds route around). New file: 12p/8xf under -O.

## Cross-links
- Extends [[issue_240_phase2_s_axes_scanmarch_spec]] (Step A; same value-preserving
  re-association + DriftWarning pattern) + [[issue_206_phase_c_verification]]
  (density-vs-scan; the loss_action carve precedent) + [[feedback_eigen_on_nonfissile_mixture]]
  (removal form as fixed-source). vv §"Bit-identity vs principled-equivalence"
  (3-criteria); vv Mode-9 (degeneracy IN the gate design — the round-trip is blind).
- coding-elegance Pattern 1+2 (composite owns its matvec; twin-path = leaf sum).
