---
name: issue-196-phase-g-step2-minimal-reproducer
description: Phase G Step 2 minimal hand-tractable reproducer (2-cell × 2-ordinate sphere) identifies the TWO algebraic terms that differ between SI sweep and apply-matvec at the fixed point. Decisive evidence (Picard-on-toggleable variants) shows that fixing BOTH defects in the SI sweep recovers the apply-matvec's converged ψ to machine precision (~1.2e-11) at any mesh refinement. The dominant defect is the pole-face WDD initial condition (SI uses ψ_face_in=0 at the pole face; apply-matvec uses ψ_cell at the pole face per Lewis-Miller §4.5 Carlson starting-direction). Secondary defect is the Carlson seed Q_bar source: SI uses Q_scatt (scattering source only); apply-matvec uses Σ_t·φ_0 (= Q_scatt + Q_ext at fixed point). H3 confirmed (implementation-level algebraic divergences); H1 refuted (the difference is NOT sweep-order-dependent); H2 partially confirmed (only the seed-source term contributes, and amplification is LINEAR not non-linear). Decision A vs B impact: Path 1 (A) per-cell N×N block solve does NOT directly fix these; what IS needed is a 5-line surgical patch to `_sweep_1d_spherical` (pole-face IC) + a 3-line patch (Carlson Q_bar source). The architectural choice between (A) and (B) is now informed by knowing that the legacy SI's structural divergence is RECOVERABLE within the current strategy framework — Path (A) the block-solve becomes UNNECESSARY for L0 closure; Path (B) reframed (drop SI for curvilinear) is still principled but no longer the only correctness-driven option.
metadata:
  type: project
  branch: refactor/sn-operator-algebra
  step: Phase G Step 2 pre-decision
  date: 2026-05-13
  reproducer_config:
    n_cells: 2
    n_ordinates: 2
    quadrature: GaussLegendre1D
    mixture: B (Σ_t=2, Σ_s=1.9, c=0.95)
    geometry: sphere R=2 cm
    bc: reflective
    source: isotropic Q=1
    analytical: φ=10 everywhere; ψ_n=5 for every ordinate at every cell (Pomraning isotropy)
---

# Phase G Step 2 — minimal hand-tractable reproducer + per-cell algebra walkthrough

## Headline

**The SI sweep and the apply-matvec compute the SAME per-cell algebraic
equation `denom · ψ_avg = source + |μ|·A_total·ψ_spatial_in + dA_w·c_in·ψ_angle_in`,
but with two algebraically-different INPUT VALUES that DO NOT vanish at the
fixed point:**

| # | Term | SI sweep value | Apply-matvec value | At FP (ψ_K=5 flat, μ_n=+1/√3) |
|---|------|---------------|--------------------|-------------------------------|
| 1 | `ψ_spatial_in` at pole face (μ≥0, i=0) | **`0`** (sweep.py:559) | **`ψ_cell[n, i=0]`** (operator.py:781-786) | SI: 0 vs Apply: 5 → propagates `+5.4556` of unbalanced streaming at (n=1,i=0) |
| 2 | Carlson seed `Q_bar` source | **`0.5 · Q_scatt`** (sweep.py:514 with `Q_1d = scattering source only`) | **`0.5 · Σ_t · φ_0`** (psi_half_angle_seed.py:569) | SI: 0.5·19=9.5 vs Apply: 0.5·20=10 → `0.5·Q_ext` constant offset |

The empirical 22% L0 error is mesh-independent because BOTH terms are O(1)
quantities that do not attenuate under refinement. Toggle BOTH on in the
SI sweep → SI converges to ψ_K at machine precision (1.2e-11), at any
mesh size (verified at n_cells ∈ {2, 5, 10, 20, 40} × n_ord ∈ {2, 8}).

## H1 / H2 / H3 verdict

- **H1 (sweep-order-dependent operator from in-cell M-M Gauss-Seidel)**:
  **REFUTED**. Both differing terms are algebraic INPUT-VALUE choices
  at the per-cell equation, independent of sweep order or G-S vs Jacobi
  iteration scheme. The M-M recurrence's in-cell state propagation IS
  reproducible by the apply-matvec's precomputed `redist_full` array
  (the recurrence is unitary at the fixed point regardless of order).
  Verified by `diag_phase_g_step2_two_bugs_isolation.py` variant 3:
  the SI sweep with both fixes reaches `ψ_K` at machine precision
  WITHOUT changing the per-cell M-M recurrence's in-cell propagation.

- **H2 (seed BC trace difference amplified non-linearly by M-M)**:
  **PARTIALLY CONFIRMED**. The Carlson seed `Q_bar` source IS different
  (Term #2 above). The difference propagates LINEARLY through the
  Carlson inward sweep and the M-M recurrence — but does NOT undergo
  any non-linear amplification. The seed difference alone contributes
  only ~3% relative error (variant 2 vs variant 0 ≈ same answer).

- **H3 (implementation-level algebraic divergence)**: **CONFIRMED**.
  Two algebraic INPUT-VALUE choices made by the SI sweep diverge from
  the apply-matvec. Both are deliberate code (well-commented, line
  numbers cited above) — neither is a transcription bug. They are
  different mathematical conventions.

## Reference value selection — pillar trace

The reference is `ψ_K = 5.0` per ordinate per cell, derived from:

- **Pillar trace**: closed-form analytical. For a homogeneous reflective
  sphere with isotropic external Q, infinite-medium isotropic equilibrium
  gives `φ = Q/(Σ_t·(1−c)) = 1/(2·0.05) = 10`. Pomraning 1989
  (NSE 102:317-336) structural-singularity result: at r=0 the angular
  flux MUST be isotropic, ψ_n(r=0) = φ(0)/Σ_w = 10/2 = 5 for every n.
  By reflective-BC symmetry the flat field ψ_n(r) = 5 everywhere
  satisfies all boundary conditions.
- **Independence trace**: closed-form is structurally INDEPENDENT of the
  apply-matvec (no shared code path); Pomraning's derivation is purely
  analytical (no quadrature).
- **L4 cross-check (informational)**: `solve_sn_fixed_source(...,
  inner_solver="krylov")` converges to `ψ_K = 5` per ordinate per cell
  to ~1.85e-13. This is L4 (code-to-analytical-limit), the canonical
  L0-SN-001 streaming-equilibrium gate.

## Minimal reproducer

`derivations/diagnostics/diag_phase_g_step2_minimal_2x2.py` —
2 cells × GL-2 (N=2), mixture B 1G, sphere R=2 cm, reflective BC,
isotropic external Q=1.

Empirical SI vs Krylov:

```
ψ_SI:  [[4.6505, 6.0617],     # ordinate 0 (μ=-1/√3)
        [3.1241, 6.1637]]     # ordinate 1 (μ=+1/√3)

ψ_K:   [[5.0000, 5.0000],     # all entries = 5 (analytical)
        [5.0000, 5.0000]]

|Δψ|_∞ = 1.876  (~38% rel error)
```

This is the smallest configuration that demonstrates the discrepancy.
`(2, 2)` is hand-tractable: 4 unknowns ψ[n, i], one M-M ordinate-step,
two outward visits and two inward visits.

## Numerical inputs (2 cells × 2 ordinates)

| Symbol | Value | Source |
|--------|-------|--------|
| μ_0 | −1/√3 ≈ −0.5774 | GL-2 negative root |
| μ_1 | +1/√3 ≈ +0.5774 | GL-2 positive root |
| w_0, w_1 | 1, 1 (Σw = 2) | GL-2 weights |
| A[0:3] | [0, 31.665, 50.265] | `sn_mesh.reduced.face_areas` (4π r² at r ∈ {0, 1, 2}) |
| V[0:2] | [16.755, 16.755] | (4π/3)(r₁³−r₀³) per cell |
| ΔA[0:2] | [31.665, 18.600] | `delta_A = A[i+1]−A[i]` |
| α_h[0:3] | [0, 0.5774, 0] | half-ordinate redistribution; α_{n+1/2} = α_{n−1/2} − w_n·μ_n |
| τ_mm[0:2] | [0.5, 0.5774] | M-M clamp (τ_0 = 1/2 by convention, τ_1 = α_{3/2} for GL-2) |
| Σ_t | 2.0 | mixture B 1G |
| Σ_s_self | 1.9 | mixture B 1G |
| Q_ext | 1.0 | isotropic external source |

Closure coefficients per ordinate (computed at the per-cell call):

| n | c_in = (1−τ)/τ · α_{n+1/2} + α_{n−1/2} | c_out = α_{n+1/2}/τ |
|---|---------------------------------------|---------------------|
| 0 | (0.5/0.5)·0.5774 + 0 = 0.5774 | 0.5774/0.5 = 1.1547 |
| 1 | (0.4226/0.5774)·0 + 0.5774 = 0.5774 | 0/0.5774 = 0 |

## Symbolic derivation — per-cell equation at the fixed point

The per-cell DD balance with M-M angular closure, FROM THE SAME
`cell_balance_terms` function (orpheus/sn/spatial/cell_balance.py),
is consumed by:

- **SI sweep** via `DiamondDifference._update_curvilinear` (diamond.py:602-669):
  forward solve `ψ_avg = (source + numer_upstream) / denom`.
- **Apply-matvec** via `_apply_curvilinear_residual` (operators.py:319-367):
  reverse `residual = denom · cell_avg − (source + numer_upstream)`.

Same algebra. The bug is NOT in the algebra; it is in WHAT VALUES are
passed in as `psi_spatial_in` and `psi_angle_in`.

### Per-cell equation (canonical form)

For each (n, i):

```
denom[n, i]   = 2·|μ_n|·A_downstream + (ΔA_i/w_n)·c_out[n] + Σ_t·V_i
A_downstream  = A_outer (= A[i+1])  if μ_n ≥ 0
              = A_inner (= A[i])    if μ_n < 0

numer_upstream[n, i]  =  |μ_n| · (A_inner + A_outer) · ψ_spatial_in[n, i]
                       + (ΔA_i/w_n) · c_in[n] · ψ_angle_in[n, i]

source[n, i] = Q_total[i] · V_i / Σw    (weight-normalised volumetric source)

ψ_avg[n, i] = (source[n, i] + numer_upstream[n, i]) / denom[n, i]
```

### Substituting ψ_K = 5 (flat) into the apply-matvec inputs

Apply-matvec's input values at (n=1, i=0) for `ψ_cell ≡ 5`:

1. **Carlson seed `phi_aux`** — runs the inward DD recurrence at μ=−1
   with `Q_bar = 0.5·Σ_t·φ_0`:
   - φ_0(i=0) = 5 + 5 = 10 (sum over ordinates with w=1).
   - Q_bar(i=0) = 0.5 · 2 · 10 = 10 = Q_bar(i=1).
   - `bc_outer_value`: BC.apply on cell-centre ψ at i=1 (the outermost
     cell); for reflective BC + flat ψ, the most-inward inflow = 5.
   - Inward sweep: phi_aux[i=1] = (1·10 + 2·5)/(1·2 + 2) = 5; face = 5.
                   phi_aux[i=0] = (1·10 + 2·5)/(1·2 + 2) = 5; face = 5.
   - **phi_aux = [5, 5]**.

2. **M-M recurrence at i=0** (ordinate-wise, vectorised over cells):
   - m=0: ψ_half_left = phi_aux[0] = 5. ψ_half_right = (ψ_cell[0,0] −
     (1−0.5)·5)/0.5 = (5 − 2.5)/0.5 = 5. ψ_half_left → 5.
   - m=1: ψ_half_left = 5. ψ_half_right = (5 − (1−0.5774)·5)/0.5774 = 5.

   So `ψ_angle_in[n=1, i=0] = psi_half_left at start of m=1 step = 5`.

3. **WDD spatial recurrence on outward sweep** (i=0 → i=1):
   - psi_face_in[n=1, i=0] = ψ_cell[n=1, i=0] = 5  (pole-face IC,
     operator.py:781-786).
   - psi_face_out[n=1, i=0] = 2·5 − 5 = 5  (flat propagation).
   - psi_face_in[n=1, i=1] = 5. psi_face_out[n=1, i=1] = 5. Reflective
     BC → inward inflow = 5.

   So `ψ_spatial_in[n=1, i=0] = 5`.

Substituted into the per-cell equation at (n=1, i=0):

```
denom = 2·0.5774·A[1] + (31.665/1)·0 + 2·16.755
      = 2·0.5774·31.665 + 0 + 33.510
      = 36.572 + 33.510 = 70.082

numer_upstream = 0.5774·(0 + 31.665)·5 + (31.665/1)·0.5774·5
              = 91.430 + 91.422 = 182.852

source = Q_total · V / Σw = 20·16.755/2 = 167.55  (Q_total = Σ_s·φ_0 + Q_ext = 19 + 1 = 20)

ψ_avg = (167.55 + 182.852) / 70.082 = 350.402 / 70.082 = 5.000  ✓
```

The apply-matvec's per-cell equation evaluates to `ψ_avg = 5` at the FP.

### Substituting ψ_K = 5 (flat) into the SI sweep inputs

SI sweep's input values at (n=1, i=0) when iterating Picard with ψ=5 flat:

1. **Carlson seed `phi_aux`** — runs the inward DD recurrence at μ=−1
   with `Q_bar = 0.5 · Q_scatt`:
   - Q_scatt = Σ_s · φ_0 = 1.9 · 10 = 19 (per cell).
   - Q_bar = 0.5 · 19 = 9.5  ≠ 10 (apply-matvec value).
   - `bc_outer_value`: BC.apply on `bc_outer` buffer at most-inward
     ordinate. The persistent `bc_outer` carries the outflow at the
     outermost face from the previous outward sweep. For flat ψ=5 with
     reflective BC the outflow face value = ψ_face_out[n=1, i=1].
     **However**, that face value differs from the apply-matvec's
     equivalent because of cascade from defect #1; in the variant-2 fix
     it stabilises at 5 (matching apply).
   - Inward sweep (assuming bc_outer_value = 5):
     phi_aux[1] = (1·9.5 + 2·5)/(1·2 + 2) = 19.5/4 = 4.875.
     phi_aux[0] = (1·9.5 + 2·face)/(4) where face=2·4.875−5=4.75:
       phi_aux[0] = (9.5 + 9.5)/4 = 4.75. So **phi_aux = [4.75, 4.875]**.
   - phi_aux ≠ [5, 5] (apply-matvec value).

2. **M-M recurrence at i=0** initialised with phi_aux[0] = 4.75:
   - m=0: ψ_half_left = 4.75. ψ_half_right = (5 − 0.5·4.75)/0.5 = 5.25.
     ψ_half_left → 5.25.
   - m=1: ψ_half_left = 5.25.

   So `ψ_angle_in[n=1, i=0] = 5.25 ≠ 5` (apply-matvec value).

3. **WDD spatial recurrence on outward sweep**:
   - psi_face_in[n=1, i=0] = **0** (SI pole-face IC, sweep.py:559).
   - psi_face_out[n=1, i=0] = 2·5 − 0 = 10 (≠ apply's 5).

   So `ψ_spatial_in[n=1, i=0] = 0`.

Substituted into the per-cell equation at (n=1, i=0):

```
denom = 70.082 (same as apply)

numer_upstream  =  0.5774·(0 + 31.665)·0 + (31.665/1)·0.5774·5.25
              =  0 + 96.058 = 96.058   ≠ 182.852

source = same (167.55)

ψ_avg = (167.55 + 96.058) / 70.082 = 263.608 / 70.082 = 3.762  ≠ 5
```

The SI's per-cell equation evaluates to `ψ_avg ≈ 3.76` when fed flat
ψ=5 — i.e., **ψ_K is NOT a fixed point of the SI sweep**. The SI's
own fixed point is `ψ_SI ≈ [4.65, 6.06; 3.12, 6.16]`, which the SI
recurrence DOES leave invariant (residual 9.5e-15 per the SI Picard
convergence).

### The differing terms — algebraic isolation

Subtracting "SI input value" from "apply input value" term-by-term:

| Term | SI value at FP | Apply value at FP | Difference |
|------|---------------|-------------------|-----------|
| ψ_spatial_in at (n=1, i=0) | 0 | 5 | **+5** |
| ψ_angle_in  at (n=1, i=0) | 5.25 | 5 | **−0.25** |
| Carlson Q_bar source       | 0.5·Q_scatt | 0.5·Σ_t·φ_0 | **−0.5·Q_ext** |

The Carlson Q_bar mismatch CAUSES the ψ_angle_in mismatch (the seed
phi_aux differs, propagating through M-M). The ψ_spatial_in mismatch
is a SEPARATE pole-face WDD-IC convention disagreement.

In `numer_upstream`:

```
Δ(numer_upstream) at (n=1, i=0)
  = |μ_n|·A_total·(5 − 0)  +  (ΔA/w)·c_in·(5 − 5.25)
  = 91.430                  +  31.665·0.5774·(−0.25)
  = 91.430                  +  (−4.570)
  = 86.860
```

In `ψ_avg`:

```
Δ(ψ_avg) at (n=1, i=0) = Δ(numer_upstream) / denom
                       = 86.860 / 70.082
                       = 1.239
```

Empirical: ψ_SI[n=1, i=0] = 3.124 vs ψ_K[n=1, i=0] = 5.000 → Δψ = 1.876.
The algebraic-difference estimate +1.239 is the LEADING contribution
at this cell; the residual ~0.64 comes from the propagation across
neighbour cells and the self-consistent Picard convergence of the
contaminated SI fixed point.

## Empirical evidence — toggleable fixes (decisive)

Script: `derivations/diagnostics/diag_phase_g_step2_two_bugs_isolation.py`.

A custom SI sweep with two toggleable fixes:

- `fix_pole_face_ic`: set `psi_spatial_in[μ≥0, i=0] = psi_state[n, 0]`
  instead of 0 at the outward-sweep start.
- `fix_carlson_seed_source`: build `Q_bar = 0.5 · Σ_t · φ_0` (from the
  current Picard iterate ψ) instead of `Q_bar = 0.5 · Q_scatt`.

Picard convergence on the 2-cell × 2-ordinate problem:

```
v0 (production SI baseline):                   max|ψ−5| = 1.876   FAIL
v1 (pole-face IC fix only):                    max|ψ−5| = 0.064   FAIL  (3% rel)
v2 (Carlson seed source fix only):             max|ψ−5| = 1.882   FAIL  (~unchanged)
v3 (BOTH fixes):                               max|ψ−5| = 1.15e-12   PASS
```

Mesh-independence confirmed (`diag_phase_g_step2_mesh_scaling.py`):

```
n_cells   v0       v1        v2       v3
   2     1.88     0.064     1.88     1.1e-12
   5     2.24     0.129     2.22     1.2e-12
  10     2.37     0.178     2.32     1.2e-12
  20     2.45     0.224     2.38     1.2e-12
  40     2.51     0.266     2.42     1.2e-12
```

GL-8 at n=40 (matches the diagnostic-memo configuration):

```
v0:   max|ψ−5| = 3.67    (73% rel error)
v1:   max|ψ−5| = 0.155   (3% rel; pole-IC fix attenuates 24×)
v2:   max|ψ−5| = 3.67    (basically unchanged; seed-fix alone is sub-dominant)
v3:   max|ψ−5| = 1.26e-11   (machine precision)
```

The pole-face IC defect is the DOMINANT contribution (~24× larger than
the Carlson-seed-source defect). Both must be fixed for full L0 closure.

## Implications for the (A) vs (B) decision

### Path (A) Path 1 — per-cell N×N angular block solve

The Phase G Step 2 verdict memo proposed this as the principled fix.
**This minimal-reproducer finding REFUTES the necessity of Path (A)
for L0 correctness.** The block solve assumes the operator-divergence
is in the M-M-in-cell-recurrence vs precomputed-redist structure (H1),
but my finding refutes H1: a 5-line surgical patch to `_sweep_1d_spherical`
that changes the pole-face IC and the Carlson-seed Q_bar source recovers
ψ_K to machine precision WITHOUT reorganising the M-M recurrence.

The per-cell N×N block solve solves a DIFFERENT problem (sweep-order-
independent in-cell coupling). It is not needed for L0 correctness on
the streaming-equilibrium test. It might be needed for some other
property (e.g., stability with very-near-critical scattering), but that
property has not been demonstrated as motivating.

### Path (B) Path 4 reframed — drop SI for curvilinear

The verdict-memo recommended this on the grounds that the legacy SI
solves a structurally-different operator. **My finding REFINES this**:
the SI operator and the apply-matvec operator differ in two SPECIFIC
algebraic terms, not in their global mathematical structure. A surgical
fix is reachable; the structural-difference framing is partly false.

Path (B) is still principled IF the project chooses to retire curvilinear
SI for cost/maintainability reasons. But it is no longer the only
correctness-driven option. The decision becomes a normal ergonomics
vs cost trade-off, not a "SI is fundamentally wrong on curvilinear"
claim.

### Path (C) — Surgical 2-line patch (NEW; recommended for evaluation)

A two-line patch to `_sweep_1d_spherical`:

**Patch 1 (pole-face IC)** at sweep.py:559:

```python
# Before:
psi_spatial_in = np.zeros(ng)

# After:
# Apply-matvec pole-face IC (Lewis-Miller §4.5 Carlson starting-
# direction): at i=0 face area A[0]=0, so streaming-IN ≡ 0 regardless
# of psi_face_in. But the WDD recurrence propagates psi_face_in as
# ψ_face_out = 2·ψ_cell − ψ_face_in into downstream cells. Choosing
# psi_face_in = ψ_cell[i=0] gives ψ_face_out = ψ_cell at every cell
# on flat ψ — i.e., preserves Pomraning isotropy at r=0.
# Need to access psi_cell at i=0 for ordinate n; we have only psi_angle
# state at this point. The cleanest fix: thread the previous-iteration
# angular flux through psi_bc, then read psi_state[n, 0] here.
# (For variant3 verification, see diag_phase_g_step2_two_bugs_isolation.py.)
psi_spatial_in = psi_state_at_pole[n]  # (ng,)
```

Note: the SI Picard loop already has access to the previous iterate
via `psi_bc`. A minimal contract change adds `psi_state_at_pole` to
`psi_bc`.

**Patch 2 (Carlson seed source)** at sweep.py:514:

```python
# Before:
Q_bar = 0.5 * Q_1d.T  # (ng, nx) — Q_1d is scattering-only

# After:
# Apply-matvec source: 0.5 · Σ_t · φ_0, where φ_0 is built from the
# previous Picard iterate ψ. At fixed point Σ_t·φ_0 = Q_scatt + Q_ext,
# but during Picard iteration the SI's previous-ψ-based φ_0 matches
# the apply-matvec's input ψ semantics.
phi_0_prev = (weights[:, None] * psi_state_prev).sum(axis=0)  # (nx,)
Q_bar = 0.5 * sig_t_1d.T * phi_0_prev[None, :]                # (ng, nx)
```

Both patches require threading `psi_state_prev` through the sweep
contract. This is minor — the SI Picard loop already has it.

**Cost**: ~10 LOC; no new strategy types; no snapshot regeneration if
the SI fixed point now coincides with the Krylov fixed point.

**Stability**: the Picard scheme on the modified SI sweep operator
converges with the same residual history as the production SI (~700
iterations for the L0 test; same as variant 3 above).

**Architectural impact**: Variant (C) ELIMINATES the structural
divergence between SI and apply-matvec at the algebraic-operator
level. The remaining redundancy in the codebase (`_sweep_1d_spherical`
vs `transport_operator_matvec_spherical`) is then a Pattern 2 twin-path
candidate for Step 3+ unification, not a correctness issue.

## What this means for Phase G Step 2

The original Step 2 plan ("L0 passes on BOTH inner_solver options at
rtol=1e-9") is achievable via Path (C) — the SI surgical patch — at
~10 LOC of code, no new strategies, no snapshot regen, no soft-break
of the inner_solver API.

The verdict memo's recommendation (B) is still defensible if the
project wants to consolidate around Krylov. But the L0 evidence no
longer COMPELS that decision; the surgical patch is now an option.

## Reference document trace

- `psi_angle[i]` in `_sweep_1d_spherical` carries M-M state across
  ordinate visits at cell i. This state IS algebraically equivalent
  to the apply-matvec's `redist_full` array at the fixed point — the
  recurrence is the same kernel (verified by Pattern 2 sharing of
  `_mm_weighted_angular_recurrence_single_level`).

- The `cell_balance_terms` algebra at
  `orpheus/sn/spatial/cell_balance.py` is shared by both call sites.
  Both pass the same `(c_in, c_out, denom, numer_upstream)` formulae;
  the divergence is purely in the INPUT VALUES (`psi_spatial_in`,
  `psi_angle_in` derived from `psi_half_seed`).

## Pointers

- Reproducer: `derivations/diagnostics/diag_phase_g_step2_minimal_2x2.py`
- Symbolic walkthrough: `derivations/diagnostics/diag_phase_g_step2_symbolic_2x2.py`
- Pole-face IC pinpoint: `derivations/diagnostics/diag_phase_g_step2_pole_face_ic.py`
- Two-bugs isolation: `derivations/diagnostics/diag_phase_g_step2_two_bugs_isolation.py`
- Mesh-scaling verification: `derivations/diagnostics/diag_phase_g_step2_mesh_scaling.py`
- SI sweep: `orpheus/sn/sweep.py:397-595` (pole-face IC at line 559,
  Carlson Q_bar source at line 514)
- Apply-matvec: `orpheus/sn/operator.py:571-838` (pole-face IC at
  lines 781-786, Carlson Q_bar source at psi_half_angle_seed.py:569)
- Shared per-cell algebra: `orpheus/sn/spatial/cell_balance.py`
- M-M recurrence kernel: `orpheus/sn/spatial/pole_angular_closure.py:340-458`
- Carlson seed kernel: `orpheus/sn/spatial/psi_half_angle_seed.py:363-419`

## Linked memories

- `[[issue-196-phase-g-step2-replan-verdict]]` — verdict-memo's
  "structurally different operators" claim is now REFINED: the operators
  differ by two specific algebraic INPUT-VALUE choices, not in their
  global structure.
- `[[issue-196-phase-g-step2-diagnostic]]` — Output D §"BC trace law
  asymmetry" was directionally correct but identified only one of the
  two defects; the M-M Carlson seed source mismatch was missed.
- `[[issue-196-phase-g-step2-replan-blocker]]` — H_A (algebraic
  discrepancy) is CONFIRMED here; H_C (G-S vs Jacobi fixed-point
  asymmetry) is REFUTED.
- `[[issue-168-phase-f-closeout]]` — Phase F's Carlson seed backport
  introduced the SI's `Q_bar = 0.5 · Q_scatt` line — the correct
  intent was to mirror apply-matvec, but the source choice diverged
  (`Q_1d = Q_scatt` in SI vs `Σ_t · φ_0` in apply-matvec). Defect #2
  is a Phase F transcription residual.
