---
name: Phase F Step 2 mesh-refinement convergence study — SN pole-cell defect diverges, outer-cell defect O(h)
description: Numerics-investigator Step 2 of Issue #168 Phase F. SN sphere_2g_3reg cell-0/cell-1 ratio diverges away from 1 with mesh refinement (0.522 at n=40 → 0.473 at n=320); cell-(N-1)/cell-(N-2) converges toward 1 at ~O(h). Variant α at matching refinement has innermost-node ratio → 1 cleanly (1.00002 at n_r=96) and outermost-node ratio → 1 (1.001 at n_r=96). SN is unambiguously the outlier. Fires Step 2 decision branch 3 (divergent → high urgency) for the pole cell and effectively branch 2 (structural defect) for the band 0.473-0.522; outer cell is branch 1 (O(h^p)). Proceed to Step 3 with two distinct suspect sites.
type: project
---

# Phase F Step 2 — SN curvilinear mesh-refinement convergence study

**Branch**: `refactor/sn-operator-algebra` 2026-05-12. Phase F Step 2
of `.claude/plans/issue_168_phase_f_curvilinear_boundary_eigenvector.md`.

**Predecessor**: Phase E (`6708a4a`); Step 1 confirmed the stored
`sphere_2g_3reg_dd_n40` snapshot reproduces bit-identically under the
current code, so the 9.09× rise (0.132 → 1.202 in g=0) is the current
SN solver's eigenvector, not a stale artefact.

## Headline finding

**The two boundary anomalies in the SN sphere_2g_3reg eigenvector have
DIFFERENT convergence behaviours, and BOTH are pathological relative
to the Variant α reference:**

- **Cell 0 (pole, r ≈ h/2)** — ratio cell-0/cell-1 **DIVERGES AWAY
  from 1** under refinement: 0.522 at n=40 → 0.473 at n=320. The
  |ratio − 1| grows from 0.478 to 0.527 across an 8× refinement.
  This is **Step 2 decision-table branch 3** (DIVERGENT, high
  urgency). The asymptote (linear fit in h) lands at ~0.473 — a
  fixed structural value, NOT the expected ~1.0 for a smooth flux.
- **Cell N-1 (outer reflective face, r = R - h/2)** — ratio
  cell-(N-1)/cell-(N-2) converges toward 1 at ~O(h¹): 0.887 →
  0.961, with log-log convergence rate ~0.76 across the last
  refinement step. This is **branch 1** (O(h^p), file follow-up).

- **Variant α at matching effective refinement** (n_r ∈ {24, 36, 48,
  72, 96} composite per-region GL nodes) has both innermost-region-0
  and outermost-region-2 ratios → 1 cleanly: innermost ratio drops
  from 1.0019 (n_r=24) to 1.00001 (n_r=96); outermost ratio from
  1.028 to 1.001. **Variant α has no boundary-cell pathology**.

The pole-cell SN behaviour combined with Variant α's clean
convergence makes the **SN spherical-pole closure** unambiguously
the responsible structural defect. The outer cell is a separate,
milder, O(h) artefact at the reflective face.

## 1. Probe scripts

### 1.1 SN probe — `/tmp/phase_f_sn_mesh_refinement.py`

Self-contained Python script. Calls `_sphere_3region("2g", n)` from
`tests/sn/regression/_generate_snapshots.py` (verified divisibility-
by-4 constraint) → `run_case(cfg)` → reads
`result.scalar_flux[:, 0, :]` of shape `(N, 2)`. Saves results to
`/tmp/phase_f_step2_sn_results.npz`.

Run with `.venv/bin/python /tmp/phase_f_sn_mesh_refinement.py`. Total
wall time ~55 s for n ∈ {40, 80, 160, 320}.

### 1.2 Variant α probe — `/tmp/phase_f_variant_alpha_mesh_refinement.py`

Self-contained Python script. Calls `solve_greens_function_sphere_mr`
from `orpheus.derivations.continuous.trajectory_resolvent.greens_function`
on the same A|B|A radii (`_MR_RADII = [0.5, 1.5, 2.0]`) and 2G XS as
the Phase E test harness `_run_sphere_2g_3reg_full`. The composite
per-region GL allocates `n_r/4 : n_r/2 : n_r/4` per region. Saves to
`/tmp/phase_f_step2_var_alpha_results.npz`.

Run with `.venv/bin/python /tmp/phase_f_variant_alpha_mesh_refinement.py`.
Total wall time ~12 min for n_r ∈ {24, 36, 48, 72, 96} (n_r=96 alone
takes 6 min — power-iteration is dominant).

Output API: `phi_g.shape = (G, n_r)` — transpose to `(n_r, G)` for
indexing parity with SN's `(N, ng)`.

## 2. SN convergence table

Problem: sphere_2g_3reg (A|B|A reflective, 2G, GL-8, R=2.0 cm).
Cells distributed per `n_per_region = (n//4, n//2, n//4)`, uniform
spacing within each region (h = 2.0/n).

```
SN
n_total  k_eff           c0/c1 g=0   cN-1/N-2 g=0  iters  conv
   40    1.3578153066    0.522335    0.887070      n/a    True
   80    1.3576649296    0.505775    0.902471      n/a    True
  160    1.3576295736    0.488757    0.934138      n/a    True
  320    1.3576226569    0.473246    0.961195      n/a    True

(iters: solve_sn doesn't expose outer_iterations on SNResult, no flag returned)

Group g=1 ratios:
n_total  c0/c1 g=1   cN-1/N-2 g=1
   40    0.513435    0.847001
   80    0.493841    0.864501
  160    0.478590    0.905047
  320    0.465600    0.942414

Profile end-points (g=0):
n_total  sf[0]        sf[1]        sf[-2]       sf[-1]
   40    0.13223      0.25314      1.20242      1.06663
   80    0.10191      0.20149      1.13270      1.02223
  160    0.07672      0.15696      1.06092      0.99104
  320    0.05638      0.11913      1.01178      0.97252
```

**k_eff** converges nicely to 1.35762266 (relative drift n=40→320 ≈
1.4e-4, expected for an O(h²) DD scheme on a 9× concentrated flux).

**Cell-0 absolute value sf[0,g=0] DECREASES under refinement** while
the rest of the eigenvector also redistributes downward as the peak
(sf[-2]) drops 1.20 → 1.01 (the L∞-normalised mode reshapes). The
relative shape — cell-0/cell-1 — drifts FURTHER from unity.

### 2.1 Convergence-rate analysis

```
|ratio - 1| under refinement (h = R/n_total):

n_total    h          |c0/c1 - 1|     |cN-1/N-2 - 1|
   40    0.050000    4.776648e-01     1.129303e-01
   80    0.025000    4.942250e-01     9.752923e-02
  160    0.012500    5.112429e-01     6.586177e-02
  320    0.006250    5.267538e-01     3.880501e-02

Log-log convergence rates (slope of log|ratio-1| vs log h):
  40 → 80   : -0.0492   (pole DIVERGES)     +0.2115   (outer)
  80 → 160  : -0.0488   (pole DIVERGES)     +0.5664   (outer)
 160 → 320  : -0.0431   (pole DIVERGES)     +0.7632   (outer ~h^0.76)

Linear extrapolation of c0/c1 in h:
  ratio = 0.473 + 1.06 * h    →  asymptote = 0.473 at h → 0
```

The negative log-log slope for the pole means **|ratio − 1| GROWS**
proportional to h^(-0.05) — slow growth but the wrong sign. The
asymptote of the cell-0/cell-1 ratio is **0.473**, a fixed structural
defect value, NOT 1.

The outer ratio converges with a healthy log-log slope of ~0.5-0.76
(roughly O(h^¾)) toward 1, consistent with a first-order BC-trace
error.

## 3. Variant α convergence table

Same A|B|A 2G XS data, composite per-region GL. Region 0 and region 2
each get n_r/4 nodes; region 1 gets n_r/2. Innermost-region-0 node
r-coordinate scales like the smallest GL abscissa within the region:
~h * (1 - cos(π/4n)) / 2 ≈ const/n² in the limit; e.g. n_r=96 →
24 nodes in region 0 → innermost at r=0.0012.

```
Variant α
n_r  k_eff           ratio_inner g=0   ratio_outer g=0   iters  conv
 24  1.3580855287    1.001949          1.027508          182    True
 36  1.3767592454    1.000432          1.008354          182    True
 48  1.3708983556    1.000148          1.004619          182    True
 72  1.3794476273    1.000031          1.001390          181    True
 96  1.3789344438    1.000010          1.001004          181    True

Group g=1:
n_r  ratio_inner g=1   ratio_outer g=1
 24  0.999356          0.989452
 36  0.999921          0.997624
 48  0.999972          0.997776
 72  0.999994          1.000405
 96  0.999998          1.000037
```

**ratio_inner → 1 monotonically and exponentially fast** (1e-2 →
1e-5). **ratio_outer → 1 also**, slightly slower but cleanly.

Note: Variant α k_eff oscillates non-monotonically with n_r (Phase E
finding: composite-GL has a quadrature floor, k_eff hops above the
SN n=40 value of 1.35782 at n_r ≥ 36). At n_r=96 the k_eff is
**1.37893** vs SN n=320's **1.35762** — Variant α is 1.6% high
relative to SN's mesh-refined k_eff. That gap is a different issue
(Phase E quadrature floor); irrelevant to the shape question.

What matters here: **at every Variant α refinement, the analog
boundary-cell ratios are within ~3% of 1 and converging toward 1**.

## 4. Step 2 decision

The plan's decision matrix:

| Branch | SN cell-0/cell-1 trend                                  | Action |
| ------ | ------------------------------------------------------- | ------ |
| 1      | Approaches 1 with rate ~h¹ or h²                        | O(h^p), file follow-up, skip Step 3 |
| 2      | Stays fixed ~0.5 across refinements                     | Structural, proceed to Step 3 |
| 3      | Goes the WRONG way (worse with refinement)              | DIVERGENT, high urgency, file ERR-NNN |
| Other  | Variant α also shows anomalies                          | Need 3rd reference, escalate |

**Pole cell (cell 0) fires BRANCH 3.** |ratio − 1| grows under
refinement (slope ~−0.05 in log-log against h, so monotonically
divergent). The asymptotic value 0.473 is also Branch-2-like (a
structural plateau), but the sign of the drift is what triggers
branch 3 specifically.

**Outer cell (cell N-1) fires BRANCH 1.** Converges toward 1 at rate
roughly O(h^¾) — slower than the O(h²) DD interior, consistent with
a first-order boundary BC-trace error. Not a high-priority defect
relative to the pole-cell finding.

**Variant α does NOT show analogous anomalies.** Both innermost and
outermost ratios converge cleanly to 1. SN is unambiguously the
outlier.

### Recommendation

**Proceed to Step 3** with the diagnosis pre-narrowed:

- **Primary suspect**: SN spherical-pole closure (cell 0). The
  cell-0/cell-1 ratio plateaus at ~0.47, indicating a structural
  defect in the **pole-cell scalar-flux integration** OR the
  **per-ordinate ψ at i=0**. Step 3 should immediately run the
  per-ordinate-ψ-vs-ordinate diagnostic at i=0 (per the plan's
  Step 3b) on the heterogeneous eigenvector. Phase D's Gate 1.1
  pinned per-ordinate identity on FLAT ψ; the failure-mode here
  is on NON-FLAT ψ — the Carlson seed is necessary but not
  sufficient, as the plan hypothesised.

- **Secondary suspect**: SN reflective-BC handling at cell N-1. The
  O(h^¾) trend toward 1 is consistent with a first-order BC-trace
  truncation; not divergent, but slower than the interior DD's
  O(h²). Lower priority than the pole defect; may resolve once the
  pole closure is fixed (the eigenmode reshapes globally).

- **Likely fix sites** (from the plan + present diagnosis):
  - `orpheus/sn/operator.py::transport_operator_matvec_spherical`
    pole-face IC at line ~738 — the WDD spatial pole-face seed for
    outgoing μ > 0 ordinates.
  - `orpheus/sn/spatial/pole_angular_closure.py::_mm_weighted_angular_recurrence_single_level`
    — the M-M angular half-angle recurrence at i=0, which Phase D
    pinned for flat ψ but not for non-flat eigenvectors.
  - `orpheus/sn/spatial/psi_half_angle_seed.py::CarlsonInwardSweep`
    — the Phase D Carlson seed has an L=0 isotropic-only assumption
    flagged in its docstring; that simplification may be the source
    of the structural pole-cell defect on non-flat eigenvectors.
  - `orpheus/sn/solver.py::_scalar_flux_from_angular` line 655 —
    the angle-integration `Σ w_n ψ_n` itself is geometry-agnostic,
    SO if per-ordinate ψ at i=0 is structurally wrong, the bug is
    upstream in the sweep, not in this helper.

## 5. Convergence flag

All four SN runs report `converged = True` from `SNResult`. The
`outer_iterations` attribute is NOT exposed on the SNResult dataclass
under the current API (the probe attempted both `outer_iterations`
and `iterations`; both returned -1, the default fallback). The
solver internally calls `solve_sn(..., max_outer=500, keff_tol=1e-12,
flux_tol=1e-10, max_inner=300, inner_tol=1e-10)` — these are the
same tight tolerances used by the snapshot generator. No
`converged=False` observed; no max-iteration saturation evident
from the k_eff values (which converge monotonically to 5 sig figs
across refinements).

All five Variant α runs report `converged = True` with
`iterations = 181-182` (well under `max_iter=500`).

The plan's "non-monotone power iteration at high refinement" risk
did NOT materialise for the SN side at n=320; k_eff converged
monotonically and tightly. Variant α k_eff is non-monotone but the
plan flagged that as a pre-existing Phase E phenomenon (composite-
GL quadrature floor), distinct from this investigation's scope.

## 6. Hypothesis check

The plan's hypothesis: *"a residual SN integration / cell-averaging
defect at curvilinear boundary cells where the per-ordinate angular
flux ψ_{n, i=0, g} is correct-on-flat-ψ (Phase D's Carlson seed
pinned this via Gate 1.1) but the angle-integrated scalar flux
φ_{i=0, g} = Σ w_n ψ_{n, i=0, g} on non-flat eigenvectors has a
residual bias. The same defect at the outer reflective face (i =
nx-1) explains why cell 39 sits below cell 38."*

**Empirical evidence SUPPORTS this hypothesis WITH AN AMENDMENT:**
the pole defect is not a bias that converges out as h → 0, it is a
**STRUCTURAL** defect that converges to a fixed ratio 0.473 (or
something near it). The "bias" wording in the hypothesis is too
weak — a bias-that-vanishes would be branch 1. The empirical
behaviour is branch 3 (divergent away from unity).

This means Phase D's Carlson seed indeed fixed the FLAT-ψ
per-ordinate identity (Gate 1.1) but the underlying spherical-pole
closure has a SECOND defect that activates only on non-flat
eigenvectors (which heterogeneous MR produces but homogeneous
problems do not — a homogeneous reflective sphere has a flat
eigenmode and would not exercise this defect).

The outer-cell defect is consistent with the plan's hypothesis: an
O(h) BC-trace truncation that converges out at refinement, and may
require attention but is **not the urgent problem**.

### Action: open a NEW ERR-NNN entry?

Per the diagnostic-cascade close-out protocol (vv-principles
"Log every caught bug"), this would warrant a new ERR-NNN if the
defect is a fresh class. Since the defect is a **manifestation of
ERR-026** (already catalogued; the curvilinear sweep WDD wrong
fixed-point) and is tracked at row #6 of its manifestation table
in the existing ERR-026 entry, **no new ERR-NNN is needed** — the
Phase F closeout (after Step 3 lands the fix) updates the existing
ERR-026 row #6 from OPEN to either CLOSED (if Step 3 succeeds) or
to a more specific subcategory.

If Step 3 surfaces a NEW failure mode distinct from ERR-026 (e.g.
a Carlson L=0 isotropic-only-assumption defect that is structurally
different from the WDD-fixed-point issue), that would justify a new
ERR-NNN at Phase F close.

## 7. Falsified alternatives

- **"Variant α is the outlier, SN is right"**: falsified. Variant α
  has independent quadrature (Sanchez-Pomraning style), independent
  spatial discretization (composite per-region GL), independent
  iteration (Schur+power-iteration). It converges to ratio 1.000 at
  the analog boundary points; the structural defect cannot be on
  the Variant α side without contradicting its convergence pattern.

- **"Both are wrong about the BC"**: falsified by Variant α's clean
  ratio_outer → 1.0010 at n_r=96. The reflective BC is correctly
  enforced by both methods at the boundary; SN's slow O(h^¾) outer
  convergence is a discretisation-truncation issue, not a BC
  interpretation issue.

- **"O(h²) discretisation error that vanishes"**: falsified for the
  pole cell. The asymptote of the ratio is 0.473, not 1. The
  |ratio − 1| GROWS with refinement, not shrinks.

- **"Stale snapshot artefact"** (Step 1 hypothesis): falsified by
  Step 1's bit-identical regeneration AND now by the fresh n ∈ {80,
  160, 320} runs all reproducing the same structural ratio band.

## 8. Files referenced

### Created (not committed per plan)

- `/tmp/phase_f_sn_mesh_refinement.py` — SN probe script
- `/tmp/phase_f_variant_alpha_mesh_refinement.py` — Variant α probe
- `/tmp/phase_f_step2_sn_results.npz` — saved SN results
- `/tmp/phase_f_step2_var_alpha_results.npz` — saved Variant α results
- `/tmp/phase_f_step2_sn_run.log` — SN run log (53 s wall)
- `/tmp/phase_f_step2_var_alpha_run.log` — Variant α run log (~12 min wall)

### Read

- `tests/sn/regression/_generate_snapshots.py:121,368` — `_sphere_3region`,
  `run_case`
- `tests/sn/test_phase_c_crosscheck.py:187,534` — Variant α XS data,
  `_run_sphere_2g_3reg_full` (reused at multiple n_r)
- `orpheus/derivations/continuous/trajectory_resolvent/greens_function.py:836,841`
  — `GreensFunctionMRResult` with `phi_g.shape=(G, n_r)`

## 9. Step 3 dispatch readiness

Per the plan's sub-agent dispatch chain:

> 3. Step 3a–3c (diagnostic deep dive): continue with numerics-
>    investigator (resume the same instance if Step 2 just
>    completed). They produce the per-ordinate diagnostic + identify
>    the fix site.

Step 3 has all the information it needs to proceed:

1. **Locked target**: pole cell of the sphere_2g_3reg MR eigenvector.
   The outer-cell defect is secondary and may resolve with the pole
   fix.
2. **Diagnostic strategy** is Step 3b in the plan: run SN on the MR
   problem; extract `result.angular_flux[g, n, i=0, 0]` per
   ordinate; check symmetry of ψ(r=0, μ) across ordinates. Pomraning
   1989 says ψ should be approximately ISOTROPIC at r=0 — if it
   isn't, the per-ordinate values are structurally wrong upstream of
   the Σ w_n ψ_n integration.
3. **Fix-site candidates** narrowed to four files (see §4
   recommendation above). The Carlson seed file's docstring WARNING
   about L=0 isotropic-only-assumption is a particularly strong
   candidate for the non-flat-eigenvector activation pattern
   observed here.

The Step 3 diagnostic dispatch should re-use the
`_sphere_3region("2g", 40)` reproducer (smallest stable case),
extract `result.angular_flux`, and produce a per-ordinate plot at
i=0 vs i=1 for both g=0 and g=1. The expected (Pomraning) signature
is near-isotropic ψ at r=0; deviation reveals the bug class.
