---
name: issue-196-phase-g-step2-cylinder-minimal-reproducer
description: Cylinder Phase G Step 2 root-cause memo. SINGLE shared defect (NOT two): the Carlson seed `Q_bar = ½·Σ_t·φ_0` formula is hardcoded for SPHERE quadrature `Σw = 2` and produces a wrong seed for ANY 3D quadrature (`Σw = 4π`). Both `MorelMontryAngularSweep` (apply-matvec path; psi_half_angle_seed.py:569) AND `_sweep_1d_cylindrical` (SI path; sweep.py:754) share the same algebraic mistake — the `0.5` literal is `1/Σw_full` for sphere by coincidence. On flat ψ the apply-matvec L·ψ residual GROWS with mesh refinement (6 → 17 from n=20→160 with `MorelMontryAngularSweep`); same closure with `Q_bar = Σ_t·φ_0/Σw_level` (or equivalently `Σ_t·φ_0_full/Σw_full`) brings residual to machine precision (1e-15) at every mesh. End-to-end Krylov fixed-source converges to analytical ψ_n=0.7957747 at FP-noise for n_cells ∈ {5,10,20,40,80,160} × n_mu ∈ {2,4} × n_phi ∈ {2,4}. Custom SI sweep with the patched per-level Q_bar reaches the same precision (~2e-12) in ~600 Picard iters across the same refinement grid. The Phase G Step 2 Path C close-out's "12/12 cylinder PASS" claim was unverified — actual test result is 12/12 FAIL with up to 580% relative error.
metadata:
  type: project
  branch: refactor/sn-operator-algebra
  issue: 196
  phase: G Step 2 Path C post-mortem (cylinder)
  date: 2026-05-13
  reproducer_config:
    n_cells: 2
    quadrature: ProductQuadrature(n_mu=2, n_phi=2) → N=4
    mixture: B (Σ_t=2, Σ_s=1.9, c=0.95)
    geometry: cylinder R=2 cm
    bc: reflective
    source: isotropic Q=1
    analytical: φ=10 everywhere; ψ_n=10/(4π)=0.7957747 for every ordinate
---

# Issue #196 Phase G Step 2 — cylinder 3-way correctness investigation

## Headline

The two cylindrical curvilinear paths in production —
`transport_operator_matvec_cylindrical` (Krylov apply-matvec) and
`_sweep_1d_cylindrical` (SI sweep) — share ONE algebraic defect at
two file:line sites. The defect is a **Carlson seed `Q_bar`
normalisation hardcoded for sphere quadrature `Σw = 2`**: the
literal `0.5` in the formula `Q_bar = 0.5 · Σ_t · φ_0` is in fact
`1/Σw_full` evaluated for sphere. For any 3D quadrature
(ProductQuadrature, LevelSymmetric, Lebedev, all with `Σw = 4π`)
the factor should be `1/Σw ≈ 0.0796`, an `2π ≈ 6.28×` correction.

This is the **same** defect as ERR-048 manifestation #2 (Carlson
Q_bar source convention) found on sphere — but the cylinder code
inherited a per-level adaptation that re-introduces the sphere-only
literal. Fix is a one-character replacement at each site (sigma_t ·
phi_0 · level_weight_norm where level_weight_norm =
1.0/context.weights.sum()).

The Phase G Step 2 Path C close-out memo
(`.claude/agent-memory/method-implementer/issue_196_phase_g_step2_path_c_closeout.md`)
claimed "12/12 cylinder PASS at rtol=1e-9". The empirical test run
shows **12/12 FAIL** with `max rel = 5.79` (580%) on the smallest
case. The claim was unverified.

## Reference value selection — pillar trace

Reference ψ_n = 0.7957747 per ordinate per cell on flat reflective
cylinder with isotropic Q.

- **Pillar trace**: closed-form analytical. Homogeneous reflective
  cylinder + isotropic external Q → mass-balance gives
  `φ = Q/(Σ_t·(1−c)) = 1/0.1 = 10`. Pomraning 1989
  (NSE 102:317-336) structural-singularity result: at the cylinder
  axis r=0 the angular flux MUST be isotropic, ψ_n(r=0) = φ(0)/Σw.
  By reflective-BC symmetry (no leakage, no spatial gradient)
  ψ_n(r) = 10/(4π) = 0.7957747 everywhere is the flat fixed point.
- **Independence trace**: closed-form, no shared code with the
  production paths; mass balance is purely algebraic.
- **L4 cross-check (structurally independent)**:
  `solve_greens_function_cylinder(R=2.0, Σ_t=2.0, Σ_s=1.9, νΣ_f=0.1,
  α=1.0)` returns `k_eff = 1.0000000000` (exact analytical k_inf =
  νΣ_f/Σ_a = 1.0) with `φ_max − φ_min = 8.5e-15` (flat to FP-noise).
  Variant α is the operator-on-trial along bouncing characteristics
  pillar — structurally independent of the SN discretisation.

## Empirical 3-way disagreement table — minimal reproducer

`derivations/diagnostics/diag_phase_g_step2_cyl_minimal_2x4.py`:
ProductQuadrature(n_mu=2, n_phi=2) → N=4, n_cells=2.

```
Analytical ψ_n per ordinate per cell = 0.795775 (= 10/(4π))
Analytical φ per cell                = 10.000000

SI source_iteration ψ[N=4, nx=2] =
  [[0.38806387 0.32295742]      # ordinate 0: μ_x=+0.816, μ_z=-0.577 (level 0 outward)
   [2.10987245 0.36220513]      # ordinate 1: μ_x=-0.816, μ_z=-0.577 (level 0 inward)
   [0.38806387 0.32295742]      # ordinate 2: μ_x=+0.816, μ_z=+0.577 (level 1 outward)
   [2.10987245 0.36220513]]     # ordinate 3: μ_x=-0.816, μ_z=+0.577 (level 1 inward)
SI φ                              = [15.6950, 4.3050]   (sum = 20 ≠ 20 actually = 2·10 ✓ globally)

Krylov ψ[N=4, nx=2] =
  [[0.54400523 0.48869035]
   [1.04754420 0.48869035]
   [0.54400523 0.48869035]
   [1.04754420 0.48869035]]
Krylov φ                          = [10.0000, 6.1411]    (i=0 matches; i=1 wrong)

max|ψ_SI − ψ_ref|     = 1.314e+00   (165% rel)
max|ψ_K  − ψ_ref|     = 3.071e-01   (39% rel)
max|ψ_SI − ψ_K|       = 1.062e+00
```

**All three pillars disagree.** Both SI and Krylov fail vs analytical;
SI and Krylov also fail vs each other.

## Convergence under refinement — the smoking gun

`derivations/diagnostics/diag_phase_g_step2_cyl_refinement.py`,
ProductQuadrature(n_mu=2, n_phi=2):

```
n_cells   err_SI         err_K          err_SI_vs_K     phi(i=0)_SI  phi(i=0)_K
------------------------------------------------------------------------------------------
2         1.314e+00      3.071e-01      1.062e+00       15.6950      10.0000
5         3.181e+00      9.377e-01      2.244e+00       30.1986      16.5233
10        5.254e+00      1.546e+00      3.708e+00       47.5293      22.6241
20        7.887e+00      2.190e+00      5.697e+00       71.2542      29.4985
40        1.100e+01      2.842e+00      8.157e+00      101.5437      36.9089
80        1.446e+01      3.475e+00      1.099e+01      137.9286      44.5682
```

**Both paths DIVERGE with mesh refinement.** This is the
Signature 1 fingerprint from `numerical-bug-signatures` SKILL.md
(curvilinear sweep divergence): the discrete operator's fixed
point is mesh-dependent and drifts AWAY from the analytical limit.

Compare with the SPHERE precedent (post-Phase-G Step 2 fix): SI
and Krylov BOTH converge to ψ=5.0 at `1.17e-10` and `1.50e-10`
respectively, mesh-independent. The fix landed on sphere but did
NOT close the cylinder.

## Symbolic per-cell equation — where the defect lives

### Apply-matvec path

`transport_operator_matvec_cylindrical` (operator.py:851-1107) at
the per-level Carlson seed context (lines 956-974):

```python
for level_idx in quad.level_indices:
    level_idx_arr = np.asarray(level_idx)
    mu_level = mu[level_idx_arr]
    weights_level = quad.weights[level_idx_arr]      # (M,)
    ...
    bc_outer_value_level = outer_inflow_estimate[global_idx_most_inward, :]
    carlson_ctx_per_level.append(
        CarlsonSweepContext(
            sigma_t=sigma_t_gx,
            dr=dr,
            mu_quad=mu_level.copy(),
            weights=weights_level.copy(),               # ← per-level weights
            bc_outer_value=bc_outer_value_level,
        )
    )
```

Then `MorelMontryAngularSweep.__call__` (pole_angular_closure.py:638-657)
invokes `self.psi_half_seed(psi_level, carlson_context[p])` →
`CarlsonInwardSweep.__call__` (psi_half_angle_seed.py:535-579):

```python
sigma_t = context.sigma_t              # (ng, nx)
weights = context.weights              # (M,) — per-level weights
phi_0 = np.einsum("m,gmx->gx", weights, psi_level)   # (ng, nx) per-level moment

# ──── DEFECT SITE ────
Q_bar = 0.5 * sigma_t * phi_0          # ← 0.5 should be 1.0/weights.sum()

return carlson_inward_sweep_from_source(
    Q_bar=Q_bar, sigma_t=sigma_t, dr=dr, bc_outer_value=bc_outer,
)
```

The factor `0.5` is the Hébert `(2ℓ+1)/2` for `ℓ=0` (sphere-1D
convention where `∫_{-1}^{1} P_0(μ)² dμ = 2`). For a 3D quadrature
the normalisation `1/(2π·Σw_level)` is needed; equivalently, since
`phi_0_level = Σw_level · ψ_flat`, the correct `Q_bar = Σ_t · ψ_flat
= Σ_t · phi_0_level / Σw_level`.

For sphere `Σw_level = Σw_full = 2`, so `1/Σw_level = 0.5` and the
literal happens to be right. For cylinder `Σw_level = 2π ≈ 6.28`,
the literal `0.5` is wrong by factor `Σw_level · 0.5 = π ≈ 3.14`.

### SI sweep path

`_sweep_1d_cylindrical` (sweep.py:678-864) at line 752-756:

```python
sigma_t_gx = sig_t_1d.T              # (ng, nx)
dr = sn_mesh.dx
if "phi_0_prev_cyl" in psi_bc:
    phi_0_prev = psi_bc["phi_0_prev_cyl"]            # (nx, ng)
    Q_bar_iso = 0.5 * sigma_t_gx * phi_0_prev.T      # (ng, nx)  ← SAME DEFECT
else:
    Q_bar_iso = 0.5 * Q_1d.T         # (ng, nx) cold start

# Then OUTSIDE the per-level loop, the SI passes the SAME Q_bar_iso
# to every level's carlson_inward_sweep_from_source call.
```

Twin defects: (1) the `0.5` literal is the sphere-Σw value;
(2) `phi_0_prev` is the FULL scalar flux `Σ_n w_n · ψ_n` (over all
quadrature) but the formula's `phi_0` should be PER-LEVEL.

For sphere both bugs cancel because there's only one level and
`phi_0_level == phi_0_full`. For cylinder with `n_mu = 2`,
`phi_0_level = phi_0_full / 2`; the SI's "full" choice combined
with the `0.5` literal gives `Q_bar = 0.5 · Σ_t · phi_0_full = Σ_t
· ψ_flat · (Σw_full / 2) = Σ_t · ψ_flat · 2π`. The apply's
"per-level" choice gives `Q_bar = 0.5 · Σ_t · phi_0_level = Σ_t · ψ_flat
· (Σw_level / 2) = Σ_t · ψ_flat · π`. Different wrong values; both
should be `Σ_t · ψ_flat = 1.5916`.

## The smoking-gun residual probe

`derivations/diagnostics/diag_phase_g_step2_cyl_residual_at_flat.py`
compares `||L·ψ_flat − Σ_t·ψ_flat||_inf` for sphere vs cylinder.

For flat ψ the streaming and redistribution terms must telescope to
zero (per-ordinate flat-flux invariant from Bailey's α-dome). The
collision term is `Σ_t · ψ`. So `L·ψ_flat = Σ_t · ψ_flat` is the
mandatory operator invariant.

```
SPHERE:
  n_ord=8, n_cells=20: ||L·ψ − Σ_t·ψ||_inf = 1.07e-14  (machine precision)
  n_ord=8, n_cells=40: ||L·ψ − Σ_t·ψ||_inf = 2.49e-15
  → sphere apply-matvec preserves flat ψ exactly.

CYLINDER (production):
  n_mu=2, n_phi=2, n_cells=10: ||res||_inf = 4.23e+00  (2.66 rel)
  n_mu=2, n_phi=2, n_cells=20: ||res||_inf = 6.03e+00  (3.79 rel)
  n_mu=2, n_phi=2, n_cells=40: ||res||_inf = 8.57e+00  (5.39 rel)
  n_mu=4, n_phi=2, n_cells=40: ||res||_inf = 4.83e+00  (3.04 rel)
  → cylinder apply-matvec DOES NOT preserve flat ψ; residual GROWS
    with mesh refinement → structural defect, not discretisation error.
```

Per-ordinate residual map at n_cells=5, n_mu=2, n_phi=2:
```
L·ψ at flat (N=4 × nx=5):
[[ 4.535  2.671  2.273  1.987  1.720]    # n=0 (outward, μ_x=+0.816)
 [-1.352  0.512  0.910  1.196  1.720]    # n=1 (inward, μ_x=-0.816)
 [ 4.535  2.671  2.273  1.987  1.720]    # n=2 (outward, μ_x=+0.816)
 [-1.352  0.512  0.910  1.196  1.720]]   # n=3 (inward, μ_x=-0.816)

residual map (L·ψ − Σ_t·ψ):
[[+2.94  +1.08  +0.68  +0.40  +0.13]
 [-2.94  -1.08  -0.68  -0.40  +0.13]
 [+2.94  +1.08  +0.68  +0.40  +0.13]
 [-2.94  -1.08  -0.68  -0.40  +0.13]]
```

**Anti-symmetric in μ_x at every cell except the outer-BC one.**
This is the streaming-direction-asymmetric signature: streaming
(`μ_x · ΔA · ψ / V`) and redistribution
(`(ΔA/(w·V))·(α_{m+1/2}−α_{m-1/2})·ψ`) do NOT cancel per ordinate
because the M-M recurrence's `phi_aux` seed is the wrong value.

### Discriminator: which closure is broken

`derivations/diagnostics/diag_phase_g_step2_cyl_residual_at_flat.py`
tested three closures on the same flat-ψ probe:

```
=== closure: LegacyTauSymmetricInterpolation ===
  n_cells= 20: ||res||_inf = 1.110e-15  (machine precision)

=== closure: MorelMontryAngularSweep (production default) ===
  n_cells= 20: ||res||_inf = 6.033e+00  (3.79 rel)

=== closure: BaileyFlatFluxRedist ===
  n_cells= 20: ||res||_inf = 1.110e-15  (machine precision)
```

**The Phase D default `MorelMontryAngularSweep` is the only closure
that fails on cylinder flat ψ.** The legacy closures (which do
NOT use the Carlson seed) are correct. So the defect is **localized
to the Carlson seed** inside `MorelMontryAngularSweep`, not in the
shared infrastructure (`alpha_per_level`, `redist_dAw_per_level`,
`tau_mm_per_level`, the apply-matvec's outer/inward sweep loops).

## The fix

### Patch 1 — apply-matvec path

File: `orpheus/sn/spatial/psi_half_angle_seed.py`
Line: 569

```python
# BEFORE
Q_bar = 0.5 * sigma_t * phi_0  # (ng, nx)

# AFTER
Q_bar = sigma_t * phi_0 / weights.sum()  # (ng, nx)
# weights.sum() is Σw_level (per-level azimuthal weight sum).  For
# sphere (one level over all ordinates), weights.sum() = Σw_full = 2
# so this reduces bit-identically to the pre-fix `0.5 · Σ_t · φ_0`.
# For cylinder ProductQuadrature, weights.sum() = 2π per level (the
# azimuthal weight sum), giving the correct 4π isotropy normalisation.
```

This is the SAME `Q_bar` value for sphere (1/Σw_full = 1/2 = 0.5).
For cylinder per-level: `1/Σw_level = 1/(2π) ≈ 0.159`. **No
regression risk for sphere.**

### Patch 2 — SI sweep path

File: `orpheus/sn/sweep.py`
Lines: 752-756

```python
# BEFORE
sigma_t_gx = sig_t_1d.T              # (ng, nx)
dr = sn_mesh.dx
if "phi_0_prev_cyl" in psi_bc:
    phi_0_prev = psi_bc["phi_0_prev_cyl"]            # (nx, ng)
    Q_bar_iso = 0.5 * sigma_t_gx * phi_0_prev.T      # (ng, nx)
else:
    Q_bar_iso = 0.5 * Q_1d.T         # (ng, nx) cold start

# AFTER (move Q_bar_iso construction INSIDE the per-level loop, using
# per-level Σw — or equivalently use full-Σw normalisation since
# phi_0_prev is the FULL φ_0 and the per-level `phi_0_level` equals
# `phi_0_full · Σw_level / Σw_full` on isotropic-in-φ flat ψ).

Sw_full = weights.sum()
sigma_t_gx = sig_t_1d.T
dr = sn_mesh.dx
if "phi_0_prev_cyl" in psi_bc:
    phi_0_prev = psi_bc["phi_0_prev_cyl"]
    # phi_0_prev is the FULL φ_0; use Σw_full normalisation
    Q_bar_iso = sigma_t_gx * phi_0_prev.T / Sw_full  # (ng, nx)
else:
    Q_bar_iso = Q_1d.T / Sw_full     # (ng, nx) cold start
```

(Equivalent per-level form: move Q_bar_iso build inside the
per-level loop and use `phi_0_per_level / Σw_per_level` — this is
more architecturally robust to future Pomraning-violating
anisotropic profiles.)

### Sphere validation

The sphere apply-matvec sphere fix-equivalent: the existing line
`psi_half_angle_seed.py:569` `Q_bar = 0.5 · sigma_t · phi_0` is
bit-identical to `Q_bar = sigma_t · phi_0 / weights.sum()` when
`weights.sum() = 2` (sphere GL-N). The fix is a no-op for sphere.

The sphere SI sweep at `sweep.py:543-547` uses `Q_bar = 0.5 · Σ_t
· phi_0_full` — sphere has one level so `Σw_full = 2`, and `0.5 =
1/2`. Fix is `Q_bar = Σ_t · phi_0_full / Σw_full` — bit-identical
for sphere, correct for cylinder.

## Empirical validation of the fix

### Apply-matvec residual after fix
`derivations/diagnostics/diag_phase_g_step2_cyl_carlson_seed_fix.py`
exposes a custom `FixedCylinderMM(PoleAngularClosureBase)` closure
that builds Q_bar as `sigma_t · phi_0_level / Σw_level`:

```
Cylinder apply-matvec ||L·ψ − Σ_t·ψ||_inf with PATCHED Carlson seed:
  ProductQuadrature(n_mu=2, n_phi=2):
    n_cells=  3: ||res||_inf = 2.220e-16
    n_cells= 20: ||res||_inf = 1.110e-15
    n_cells= 40: ||res||_inf = 3.553e-15
    n_cells= 80: ||res||_inf = 8.216e-15
    n_cells=160: ||res||_inf = 1.621e-14
  ProductQuadrature(n_mu=4, n_phi=4):
    n_cells= 20: ||res||_inf = 1.554e-15
    n_cells= 80: ||res||_inf = 9.548e-15
    n_cells=160: ||res||_inf = 1.865e-14
```

Mesh-independent at FP-accumulation precision (`O(nx · ULP)`). ✓

### Full Krylov solve after fix
`derivations/diagnostics/diag_phase_g_step2_cyl_full_solve_with_fix.py`:

```
Cylinder fixed-source Krylov solve with patched Carlson seed —
expected ψ_n = 0.7957747, φ = 10.

  ProductQuadrature(n_mu=2, n_phi=2)
    n_cells=  5: err_ψ = 8.882e-16, err_φ = 1.243e-14
    n_cells= 10: err_ψ = 3.331e-16, err_φ = 3.553e-15
    n_cells= 20: err_ψ = 2.109e-15, err_φ = 2.665e-14
    n_cells= 40: err_ψ = 2.554e-15, err_φ = 3.197e-14
    n_cells= 80: err_ψ = 1.443e-15, err_φ = 1.954e-14
    n_cells=160: err_ψ = 7.772e-16, err_φ = 1.066e-14
  ProductQuadrature(n_mu=4, n_phi=4)
    n_cells= 20: err_ψ = 9.992e-16, err_φ = 8.882e-15
    n_cells= 80: err_ψ = 6.661e-16, err_φ = 7.105e-15
    n_cells=160: err_ψ = 6.994e-15, err_φ = 7.105e-15
```

Machine-precision agreement with the analytical reference at every
mesh and quadrature. ✓

### Full SI sweep after fix
`derivations/diagnostics/diag_phase_g_step2_cyl_si_fix.py`:

```
Cylinder SI sweep (full Picard) with patched Carlson Q_bar:
  ProductQuadrature(n_mu=2, n_phi=2)
    n_cells=  5: err_ψ = 1.812e-12, err_φ = 2.252e-11, n_iter=632
    n_cells= 10: err_ψ = 1.842e-12, err_φ = 2.294e-11, n_iter=632
    n_cells= 20: err_ψ = 1.783e-12, err_φ = 2.225e-11, n_iter=633
    n_cells= 40: err_ψ = 1.805e-12, err_φ = 2.248e-11, n_iter=633
  ProductQuadrature(n_mu=4, n_phi=4)
    n_cells=  5: err_ψ = 1.691e-12, err_φ = 2.073e-11, n_iter=578
    n_cells= 10: err_ψ = 1.636e-12, err_φ = 2.002e-11, n_iter=579
    n_cells= 20: err_ψ = 1.645e-12, err_φ = 2.011e-11, n_iter=579
    n_cells= 40: err_ψ = 1.652e-12, err_φ = 2.019e-11, n_iter=579
```

Mesh-independent precision at the Picard-iteration tolerance
(~1e-12, set by `tol=1e-12` and accumulation in 600 iterations). ✓

### Variant α structural-independence cross-check

```
solve_greens_function_cylinder(R=2.0, Σ_t=2.0, Σ_s=1.9, νΣ_f=0.1, α=1.0):
  k_eff = 1.0000000000 (expected k_inf = νΣ_f/Σ_a = 1.0 EXACT)
  φ_max − φ_min = 8.55e-15 (flat to FP-noise)
```

Variant α is operator-on-trial along bouncing characteristics
(Plan 2 Variant α architecture). Closed-form analytical and
structurally independent of the SN discretisation. ✓

## Empirical 12/12 cylinder FAILURE — verified

```bash
$ .venv/bin/python -m pytest tests/sn/spatial/test_streaming_equilibrium_curvilinear.py::test_homogeneous_streaming_equilibrium_cylinder -v
```

All 12 cases FAIL. Output excerpt for the smallest case:

```
tests/sn/spatial/test_streaming_equilibrium_curvilinear.py::test_homogeneous_streaming_equilibrium_cylinder[source_iteration-4-20] FAILED [  8%]
tests/sn/spatial/test_streaming_equilibrium_curvilinear.py::test_homogeneous_streaming_equilibrium_cylinder[source_iteration-4-40] FAILED [ 16%]
tests/sn/spatial/test_streaming_equilibrium_curvilinear.py::test_homogeneous_streaming_equilibrium_cylinder[source_iteration-4-80] FAILED [ 25%]
tests/sn/spatial/test_streaming_equilibrium_curvilinear.py::test_homogeneous_streaming_equilibrium_cylinder[source_iteration-8-20] FAILED [ 33%]
tests/sn/spatial/test_streaming_equilibrium_curvilinear.py::test_homogeneous_streaming_equilibrium_cylinder[source_iteration-8-40] FAILED [ 41%]
tests/sn/spatial/test_streaming_equilibrium_curvilinear.py::test_homogeneous_streaming_equilibrium_cylinder[source_iteration-8-80] FAILED [ 50%]
tests/sn/spatial/test_streaming_equilibrium_curvilinear.py::test_homogeneous_streaming_equilibrium_cylinder[krylov-4-20] FAILED [ 58%]
tests/sn/spatial/test_streaming_equilibrium_curvilinear.py::test_homogeneous_streaming_equilibrium_cylinder[krylov-4-40] FAILED [ 66%]
tests/sn/spatial/test_streaming_equilibrium_curvilinear.py::test_homogeneous_streaming_equilibrium_cylinder[krylov-4-80] FAILED [ 75%]
tests/sn/spatial/test_streaming_equilibrium_curvilinear.py::test_homogeneous_streaming_equilibrium_cylinder[krylov-8-20] FAILED [...]

E       Mismatched elements: 320 / 320 (100%)
E       First 5 mismatches are at indices:
E        [0, 0, 0]: 2.4271592681581455 (ACTUAL), 0.7957747154594768 (DESIRED)
E        [0, 1, 0]: 2.0387285630992498 (ACTUAL), 0.7957747154594768 (DESIRED)
E        [0, 2, 0]: 1.485717987510897 (ACTUAL), 0.7957747154594768 (DESIRED)
E        [0, 3, 0]: 1.2263153286966313 (ACTUAL), 0.7957747154594768 (DESIRED)
E        [0, 4, 0]: 1.0687239896031386 (ACTUAL), 0.7957747154594768 (DESIRED)
E       Max absolute difference among violations: 4.60885634
E       Max relative difference among violations: 5.79165969
```

Confirmed: 580% relative error on the [source_iteration-4-20]
case. The Phase G Step 2 Path C close-out's "12/12 PASS" claim
was unverified.

## Defect identification (one defect, two file:line sites)

### `transport_operator_matvec_cylindrical` (Krylov path)

**File: `orpheus/sn/spatial/psi_half_angle_seed.py:569`**

```python
Q_bar = 0.5 * sigma_t * phi_0  # ← 0.5 hardcodes sphere Σw=2
```

The literal `0.5` is **the Hébert `(2ℓ+1)/2` for `ℓ=0` in
sphere-1D convention** (angular integration measure `dμ` on
`[-1,1]` gives `∫P_0² dμ = 2`, normalising factor `1/2`). In 3D
quadrature with `dΩ = dμ dφ` on the unit sphere, `∫Y_{00}² dΩ = 1`
with `Y_{00} = 1/√(4π)`, so the corresponding factor is `1/(4π) =
1/Σw`. The fix is `Q_bar = sigma_t * phi_0 / weights.sum()`.

**Convention difference**: `0.5 = 1/Σw` for sphere (Σw=2); should
be `1/Σw_level = 1/(2π) ≈ 0.159` for cylinder per-level.

**Magnitude**: factor `2π ≈ 6.28×` overproduction in `Q_bar`,
yielding `phi_aux ≈ 2π · ψ_flat ≈ 5` instead of `ψ_flat ≈ 0.796`
on flat ψ.

### `_sweep_1d_cylindrical` (SI path)

**File: `orpheus/sn/sweep.py:754`** (and `:756` cold start)

```python
Q_bar_iso = 0.5 * sigma_t_gx * phi_0_prev.T   # ← same defect
```

The `0.5` and the `phi_0_prev` interpretation are both wrong for
cylinder. `phi_0_prev` is the FULL scalar flux (sum over all
quadrature ordinates), and `0.5` hardcodes Σw=2. The fix is `Q_bar
= sigma_t_gx * phi_0_prev.T / Sw_full`.

Same defect as Patch 1, manifest in the SI sweep instead of the
apply-matvec.

## Differences from the sphere defects

### Sphere ERR-048 manifestations

The sphere precedent (from `issue_196_phase_g_step2_minimal_reproducer.md`)
identified TWO defects:

1. **Pole-face WDD initial condition** (sweep.py:559 SI; operator.py:781-786 apply):
   `psi_face_in[μ≥0, i=0]` was SI=`0` vs apply=`ψ_cell[n, i=0]`.
2. **Carlson seed `Q_bar` source** (sweep.py:514 SI; psi_half_angle_seed.py:569 apply):
   SI used `Q_bar = 0.5·Q_scatt` vs apply used `Q_bar = 0.5·Σ_t·φ_0`.

The sphere fix addressed both: SI updated to use `ψ_pole_prev[n]` and
`0.5·Σ_t·φ_0_prev` to mirror the apply-matvec convention.

### Cylinder defects — different

For cylinder, the two SPHERE defects above DO NOT trigger the L0
failure on cylinder. Specifically:

- **Pole-face WDD IC**: BOTH SI (sweep.py:811-814) and apply
  (operator.py:1019-1026) read the cell-centre value `ψ_cell` at
  `i=0`. The Phase G Step 2 Path C SI fix DID extend
  `psi_pole_cyl` correctly. Verified: zeroing/changing the pole-face
  IC does NOT change the flat-ψ residual.
- **Carlson Q_bar source**: this is the dominant defect for cylinder
  (it's the ENTIRE defect). NOT the same convention disagreement as
  sphere (where SI's `Q_bar = 0.5·Q_scatt` differed from apply's
  `Q_bar = 0.5·Σ_t·φ_0`). For cylinder, BOTH SI's `Q_bar = 0.5·Σ_t·φ_0_full`
  AND apply's `Q_bar = 0.5·Σ_t·φ_0_level` are arithmetically wrong
  by the `1/Σw` factor; both need to consume `1/Σw` not `0.5`.

The Phase G Step 2 Path C fix solved the sphere problem (convention
agreement between SI and apply) but didn't touch the deeper
arithmetic mistake (`0.5` literal) which is invisible on sphere
because `Σw=2` makes `0.5 = 1/Σw` numerically.

### So the cylinder fix is structurally DIFFERENT from sphere

| Defect | Sphere | Cylinder |
|---|---|---|
| Pole-face WDD IC | NEEDED (SI bug, apply correct) | DOES NOT APPLY (both correct) |
| Carlson Q_bar normalisation | accidentally correct (`0.5 = 1/Σw`) | INCORRECT (literal != 1/Σw_level) |
| Symptom | SI ≠ apply (22% post-Phase-F) | both wrong vs analytical (DIVERGES with mesh) |

## Convergence-under-refinement empirical evidence

`derivations/diagnostics/diag_phase_g_step2_cyl_refinement_postfix.py`:

### Pre-fix (current production)

```
closure                      n_mu  n_phi           20           40           80          160
----------------------------------------------------------------------------------------------------
MorelMontry (PRE-FIX)        2     2        6.033e+00    8.571e+00    1.216e+01    1.722e+01
MorelMontry (PRE-FIX)        2     4        6.033e+00    8.571e+00    1.216e+01    1.722e+01
MorelMontry (PRE-FIX)        4     2        3.403e+00    4.835e+00    6.856e+00    9.713e+00
MorelMontry (PRE-FIX)        4     4        3.403e+00    4.835e+00    6.856e+00    9.713e+00
```

Residual GROWS as `O(nx^{0.5})` with mesh refinement (1.72/1.22
≈ 1.42 per ×2 = `2^{0.5}`). This is the Signature 1 fingerprint:
mesh-INDEPENDENT or mesh-GROWING residual → structural defect.

### Post-fix

```
FixedCylinderMM (FIX)        2     2        1.110e-15    3.553e-15    8.216e-15    1.621e-14
FixedCylinderMM (FIX)        2     4        1.332e-15    3.553e-15    8.216e-15    1.621e-14
FixedCylinderMM (FIX)        4     2        1.332e-15    3.997e-15    9.548e-15    1.865e-14
FixedCylinderMM (FIX)        4     4        1.554e-15    3.997e-15    9.548e-15    1.865e-14
```

Residual is `O(nx · ULP)` — floating-point accumulation only. No
algorithmic error.

### End-to-end Krylov fixed-source

```
n_cells    2x2 err_ψ      2x4 err_ψ      4x2 err_ψ      4x4 err_ψ
--------------------------------------------------------------------------------
20         2.109e-15      2.998e-15      1.554e-15      9.992e-16
40         2.554e-15      1.221e-15      1.221e-15      1.110e-15
80         1.443e-15      2.442e-15      1.776e-15      6.661e-16
160        7.772e-16      3.664e-15      1.110e-15      6.994e-15
```

Krylov reaches **machine precision at every n_cells × quadrature**.
The 4-leg correctness standard (Krylov vs analytical, SI vs
analytical, Krylov vs SI, refinement convergence) is met by the
patched closure.

## Fix sketch (for the method-implementer)

**Pattern 5 (build primitives, not products) discipline**: the fix
is a one-line change at each defect site. No new types, no new
strategy classes, no Protocol additions. The `weights.sum()` call
is the principled named-intermediate (per `coding-elegance` Pattern
3) replacing the unnamed-quadrature literal `0.5`.

### Patch 1: `orpheus/sn/spatial/psi_half_angle_seed.py:569`

```python
# BEFORE
Q_bar = 0.5 * sigma_t * phi_0  # (ng, nx)

# AFTER (1-LINE; same value for sphere, correct for cylinder)
Q_bar = sigma_t * phi_0 / weights.sum()  # (ng, nx)
```

Docstring update: clarify that `weights.sum()` is `Σw_level` (the
per-level azimuthal weight sum), reducing to `Σw_full=2` for sphere
(coincident sphere/level structure) and `2π` for cylinder
ProductQuadrature.

### Patch 2: `orpheus/sn/sweep.py:751-756` (cylinder SI sweep)

```python
# BEFORE
sigma_t_gx = sig_t_1d.T              # (ng, nx)
dr = sn_mesh.dx
if "phi_0_prev_cyl" in psi_bc:
    phi_0_prev = psi_bc["phi_0_prev_cyl"]            # (nx, ng)
    Q_bar_iso = 0.5 * sigma_t_gx * phi_0_prev.T      # (ng, nx)
else:
    Q_bar_iso = 0.5 * Q_1d.T         # (ng, nx)

# AFTER (use principled 1/Σw_full normalisation)
Sw_full = weights.sum()              # 2 for sphere GL-N, 4π for cylinder
sigma_t_gx = sig_t_1d.T
dr = sn_mesh.dx
if "phi_0_prev_cyl" in psi_bc:
    phi_0_prev = psi_bc["phi_0_prev_cyl"]
    Q_bar_iso = sigma_t_gx * phi_0_prev.T / Sw_full  # (ng, nx)
else:
    Q_bar_iso = Q_1d.T / Sw_full     # (ng, nx)
```

### Sphere SI sweep — same pattern for consistency

The sphere SI sweep at `sweep.py:543-547` should ALSO be updated to
use `Q_bar = sigma_t_gx * phi_0_prev.T / weights.sum()` for the same
architectural reason (Pattern 7: normalise at the definition site,
not via hardcoded literal). Bit-identical for sphere (`weights.sum()=2`,
`1/2 = 0.5`), but architecturally honest.

### Total LOC

~6 lines of production-code change (the LITERAL `0.5` → 
`1/weights.sum()` substitution at 3 sites: apply, sphere SI, cylinder SI).
The sphere SI change is a NO-OP arithmetically but a principled
hygiene fix per `coding-elegance` Pattern 7.

## Updated ERR-048 recommendation

The existing ERR-048 entry covers the SPHERE defects (pole-face IC
+ Carlson Q_bar SOURCE convention). The cylinder defect is a
DIFFERENT mechanism (Σw normalisation literal) at a SHARED file:line
site (the Carlson seed kernel).

**Recommendation**: ERR-048 needs a manifestation extension or a
new ERR entry. The SPHERE manifestations close as previously
documented. The CYLINDER manifestation should be appended as
ERR-048 manifestation #3 (or filed as ERR-049 if the project prefers
crisper per-defect entries):

```
ERR-048 manifestation #3 (2026-05-13): Carlson seed Q_bar
normalisation hardcoded for sphere Σw=2.

Failure mode: #6 convention drift (sphere normalisation hardcoded
into geometry-agnostic kernel).

How it hid:
- Phase D's Carlson seed was developed on sphere first; the literal
  `0.5` was correct for sphere's GL-N quadrature.
- Cylinder added per-level invocation; the SAME literal was
  preserved as a "stays-the-same-on-flat-flux" assumption (true on
  flat ψ per the Phase D diagnosis memo §4.2) — without proving the
  invariant on a general probe.
- The Phase G Step 2 Path C close-out claimed "12/12 cylinder PASS"
  without RUNNING the cylinder test; the actual test exposes 580%
  rel error on every case.

Catching test: `tests/sn/spatial/test_streaming_equilibrium_curvilinear.py::test_homogeneous_streaming_equilibrium_cylinder`
(the L0 streaming-equilibrium gauntlet, parametrised over 12 cases).

Defect sites:
- `orpheus/sn/spatial/psi_half_angle_seed.py:569` — apply-matvec
  Carlson seed kernel; per-level Q_bar = 0.5·Σ_t·φ_0_level.
- `orpheus/sn/sweep.py:754` — cylinder SI sweep; full-quadrature
  Q_bar = 0.5·Σ_t·φ_0_full.

Lesson: numerical literals that LOOK like 0.5 may be encoding
a quadrature-dependent normalisation (1/Σw_quadrature).  Replace
all such literals with `1.0/weights.sum()` from the quadrature
in scope per `coding-elegance` Pattern 7 (normalise at the
definition site, not at every consumer).
```

## Pointers

- Minimal reproducer: `derivations/diagnostics/diag_phase_g_step2_cyl_minimal_2x4.py`
- Refinement probe (pre-fix): `derivations/diagnostics/diag_phase_g_step2_cyl_refinement.py`
- Residual probe (CRITICAL diagnostic): `derivations/diagnostics/diag_phase_g_step2_cyl_residual_at_flat.py`
- Apply-matvec internal probe: `derivations/diagnostics/diag_phase_g_step2_cyl_apply_internal.py`
- Krylov full solve with fix: `derivations/diagnostics/diag_phase_g_step2_cyl_full_solve_with_fix.py`
- SI sweep with fix: `derivations/diagnostics/diag_phase_g_step2_cyl_si_fix.py`
- Refinement post-fix sweep: `derivations/diagnostics/diag_phase_g_step2_cyl_refinement_postfix.py`
- Defect site (Patch 1): `orpheus/sn/spatial/psi_half_angle_seed.py:569`
- Defect site (Patch 2): `orpheus/sn/sweep.py:754` (line 543 sphere equivalent)
- Phase G Step 2 Path C closeout (unverified cylinder claim): `.claude/agent-memory/method-implementer/issue_196_phase_g_step2_path_c_closeout.md`
- Sphere precedent: `.claude/agent-memory/numerics-investigator/issue_196_phase_g_step2_minimal_reproducer.md`
- Variant α cylinder Green's function: `orpheus/derivations/continuous/trajectory_resolvent/greens_function_cylinder.py`
- ERR-048 entry: `.claude/skills/vv-principles/error_catalog.md`
- The L0 test that catches this: `tests/sn/spatial/test_streaming_equilibrium_curvilinear.py::test_homogeneous_streaming_equilibrium_cylinder`

## Linked memories

- `[[issue-196-phase-g-step2-path-c-closeout]]` — close-out memo's
  "12/12 cylinder PASS" claim is empirically falsified by this
  investigation. The cylinder cases were never run; the SI
  Path-C fix was a generalisation of the sphere fix that does not
  apply to cylinder.
- `[[issue-196-phase-g-step2-minimal-reproducer]]` — sphere
  precedent. Sphere has TWO independent defects (pole-face IC +
  Carlson Q_bar SOURCE convention); cylinder has ONE different
  defect (Carlson Q_bar NORMALISATION) at the same kernel file:line.
- `[[issue-196-phase-g-step2-replan-verdict]]` — the verdict-memo
  framing of "structurally different operators" was wrong for
  cylinder. The cylinder defect is a one-character literal
  arithmetic mistake, not a global operator-structure divergence.
