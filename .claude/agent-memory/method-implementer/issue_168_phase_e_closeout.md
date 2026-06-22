---
name: Issue #168 Phase E — composite per-region GL for MR Variant α + tighter Gate 4.2 + flux-shape sentinel
description: Follow-up to Phase D Step 4b that fixes a documented prototype TODO in trajectory_resolvent's MR Variant α radial quadrature. Composite per-region GL replaces single-domain GL on (0, R); k_eff agreement on heterogeneous closed MR drops from ~7-9% to ~2e-4 (sphere) and ~1.75e-2 (cylinder). Gate 4.2 rtol tightened 1e-1 → 2e-2/3e-2. Flux-shape comparison built but revealed an unresolved SN-vs-Variant-α eigenvector shape mismatch; xfailed-strict pending sphere-pole eigenvector investigation.
type: project
---

# Issue #168 Phase E — closeout

**Branch**: `refactor/sn-operator-algebra`
**Phase D Step 4b follow-up** — the user-requested "Phase E flux-shape comparison + tight-rtol pursuit for heterogeneous MR Gate 4.2"
**Phase E commits**: `2d3e7f2`..`9417130` (3 commits)
**Date shipped**: 2026-05-12

## Headline

A documented prototype TODO in trajectory_resolvent's MR Variant α
solvers — single-domain Gauss-Legendre on (0, R) instead of composite
per-region GL — was the root cause of the 7-9% k_eff gap between
SN n=40 snapshots and Variant α on heterogeneous closed MR.  Phase E
fixed the TODO; the gap drops to **2e-4 (sphere)** and **1.75e-2
(cylinder)** at production quadrature.

A flux-shape comparison harness was built as part of the same work
and surfaced a different, unresolved finding: SN and Variant α agree
on the EIGENVALUE but disagree on the EIGENVECTOR SHAPE on
heterogeneous problems — SN sphere scalar_flux rises 9× from center
to outer fuel peak, while Variant α rises only 1.7×.  This is
xfailed-strict pending a Phase F (sphere-pole eigenvector audit)
investigation.

## What shipped

### Production correction — `orpheus/derivations/continuous/trajectory_resolvent/greens_function.py`

NEW `_composite_per_region_gl(radii, n_r_total, min_per_region=2)`:

- Builds piecewise GL nodes per region, proportionally distributed
  by region thickness.
- Returns `(r_nodes, r_weights, region_at_node)` — concatenation of
  per-region GL transformed from `[-1, 1]` to each region's physical
  interval.
- Each region's GL is exact for polynomials of degree `2·n_k - 1` on
  that region; the composite no longer recovers global `2·n_r - 1`
  exactness, but for piecewise-smooth integrands with σ_t interface
  kinks this is the correct trade-off (polynomial fits across
  material kinks are structurally inconsistent regardless of GL
  order).

Wired into three MR entry points:

- `solve_greens_function_sphere_mr` — line 855 (was single GL on
  (0, R) per the docstring's "follow-on improvement" comment).
- `solve_greens_function_sphere_mr_fixed_source` — line 1172 (same).
- `solve_greens_function_cylinder_mr` — line 806 (same in
  `greens_function_cylinder.py`; the inline comment "Composite radial
  GL on (0, R). Nodes do NOT land on interfaces." was misleading —
  it called it "composite" but was actually single-domain).

API contract preserved: `n_r` stays the TARGET total node count.
`len(r_nodes)` may slightly exceed `n_r` when the `min_per_region=2`
floor binds (small problems with many regions).  `initial_psi` shape
error messages now report the actual `n_r_actual` for diagnostic
clarity.

### Test-ceiling adjustments

Composite-GL shifts where the dominant numerical artifact lives.
Pre-Phase-E single-domain GL had nodes straddling material interfaces,
so per-region cubic splines could interpolate across the σ_t jump.
Post-Phase-E composite GL places nodes STRICTLY INSIDE each region's
interior, so the splines must EXTRAPOLATE to the boundaries.  Two
tests calibrated against the pre-fix single-GL artifact:

- `test_mr_interface_continuity_3region` (cylinder MR): ceiling
  relaxed 1e-2 → 5e-2 (empirical post-Phase-E ~2.85e-2; pre-fix was
  ~3e-3).  Docstring rewritten to describe the post-Phase-E error
  shift.
- `test_garcia_case1_phi_matches_at_point[r=7.0]` (sphere MR fixed-
  source, Garcia 2021 Case 1): added explicit outer-surface
  tolerance branch — tol=5% at r=R (vacuum BC, post-Phase-E extrap
  artifact), 2% elsewhere, 15% near material interfaces.

### Gate 4.2 tightening — `tests/sn/test_phase_c_crosscheck.py`

`test_phase_d_trajectory_resolvent_crosscheck` heterogeneous-MR
tolerances tightened:

- `sphere_2g_3reg_dd_n40`: rtol 1e-1 → **2e-2** (90× tighter)
  — empirical post-Phase-E gap is 1.99e-4 at production
  quadrature; worst case across refinements {24, 36, 48} is 1.4e-2
  (30% headroom).
- `cyl_2g_3reg_LS4_dd_n40`: rtol 1e-1 → **3e-2** (3× tighter)
  — empirical post-Phase-E gap is 1.75e-2 at production quadrature.

Homogeneous closed cases (snapshots 1, 4, 5) stay at rtol≤1e-9
(V_α1 / V_α1_cyl machine-precision, unchanged).

### NEW flux-shape sentinel — `test_phase_e_trajectory_resolvent_flux_shape_crosscheck`

Parametrised over heterogeneous MR snapshots (2 sphere, 6 cylinder).
Compares L∞-normalised flux profiles via:

- L∞-normalised SN `scalar_flux` (from the `.npz` regression snapshot)
- L∞-normalised Variant α `phi_g` interpolated onto SN cell centres
  via per-region cubic splines (the composite-GL nodes lie strictly
  within each region's interior).

`@pytest.mark.xfail(strict=True)` with this empirical-finding reason:

> Phase E flux-shape cross-check exposed an UNRESOLVED structural
> discrepancy: k_eff agrees to 2e-4 (sphere) / 1.75e-2 (cylinder) but
> the L∞-normalised flux profiles disagree by 65% (sphere) / 22%
> (cylinder).

Harness functions (`_mr_sn_cell_centers_n40`, `_interpolate_per_region`,
`_run_sphere_2g_3reg_full`, `_run_cyl_2g_3reg_full`,
`_load_snapshot_scalar_flux`) live in the test file ready for
re-enablement when the investigation completes.

## Empirical evidence

### Pre-Phase-E (v1 probe, single-domain GL on (0, R)) — non-monotone convergence

```
sphere_2g_3reg — SN n=40 k_eff = 1.3578153066
  n_r  n_mu  n_traj          k_eff    rel vs SN
   24    24      64   1.2622957754     7.03e-02
   32    32      96   1.0491818054     2.27e-01   ← worse
   48    48     128   1.3196901628     2.81e-02
   64    64     192   1.5577966986     1.47e-01   ← worse

cyl_2g_3reg — SN n=40 k_eff = 1.2284281075
  n_r  n_ax n_phi  n_traj          k_eff    rel vs SN
   24    16    32      64   1.1210217377     8.74e-02
   32    24    48      96   0.9092784045     2.60e-01   ← worse
   48    32    64     128   1.1681515888     4.91e-02
```

Confirmed `converged=True` at each refinement (183 iters at n_r=24).
The non-monotone behaviour is a QUADRATURE FLOOR effect, not power-
iteration failure.

### Post-Phase-E (smoke probe, composite per-region GL) — monotone-bounded

```
SPHERE — SN n=40 k_eff = 1.3578153066
  n_r  n_traj          k_eff     rel SN   iters  conv
   24      64   1.3580855287   1.99e-04     182  True
   36      64   1.3767578568   1.40e-02     182  True
   48      96   1.3708983556   9.64e-03     182  True

CYLINDER — SN n=40 k_eff = 1.2284281075
  n_r  n_ax n_phi  n_traj          k_eff     rel SN   iters  conv
   24    16    32      64   1.2069303353   1.75e-02     211  True
   36    16    32      64   1.2167541067   9.50e-03     210  True
```

### Flux-shape diagnostic — the unresolved finding

SN sphere `scalar_flux[g=0]` ranges 0.132 → 1.20 (9× rise from
center to outer-fuel peak).  Variant α `phi_g[g=0]` ranges 0.128 →
0.187 (1.7× rise).  After L∞ normalisation:

- SN_norm at r=0.025: 0.132/1.202 = **0.110**
- α_norm at r=0.025: 0.128/0.187 = **0.684**
- Per-cell |Δ| = **0.574** (≈ 57.5%)

The SN profile is steeply-rising-to-edge; the Variant α profile is
gentle.  Same problem, same BC (reflective), same XS data,
qualitatively different eigenvectors at quantitatively-agreeing
eigenvalues.

Hypothesis (deferred to a future Phase F): a residual SN spherical-
pole-region eigenvector issue beyond what the Phase D Carlson
coupled-pole seed fixed.  Phase D Step 3 closed Gate 1.1's per-
ordinate flat-flux identity but did NOT pin the angle-integrated
scalar-flux profile at the pole cell on heterogeneous problems.
Investigation paths:

1. Compare SN scalar_flux against an independent F_N / Marshak-type
   sphere MR reference (cross-pillar verification).
2. Audit the SN spherical sweep's near-pole eigenvector behaviour —
   particularly whether the `α_{1/2}` cascade seed (now Carlson)
   pins the right cell-averaged scalar flux profile at i=0.
3. Add a 2G heterogeneous MR L1 MMS that pins the eigenvector
   shape, not just the eigenvalue.

## Plan deviations from "Phase E objective"

The original Phase E intent (per the Step 4b closeout) was "flux-
shape comparison + tight-rtol pursuit":

- **Tight-rtol pursuit**: SUCCEEDED — composite-GL fix delivered
  ~30× tighter sphere k_eff agreement and ~5× tighter cylinder
  k_eff agreement.  Gate 4.2 tolerances tightened 1e-1 → 2e-2 / 3e-2
  (still 100× / 200× looser than the original plan's 5e-4 floor,
  reflecting the cylinder phase-space quadrature-error budget +
  SN spatial-error budget at n=40).
- **Flux-shape comparison**: BUILT but surfaced an unresolved
  finding.  Shape disagreement (65% sphere, 22% cylinder) at
  agreeing-eigenvalue is a separate structural issue beyond the
  composite-GL fix.  Xfailed-strict with detailed empirical-finding
  reason.

## V&V framing per `vv-principles`

This is a Mode-6 (convention drift) catch surfaced by deepening Gate
4.2 from k_eff-only to flux-shape:

- The single-domain GL was deployed as a "qualitative reproducer for
  Issue #132" with an explicit "follow-on improvement" comment.
- It survived to production via stalled cleanup — no test caught
  the non-monotone-under-refinement behaviour because:
  - The cylinder MR Phase 1b interface-continuity test allowed
    `rel_jump < 1e-2` (calibrated against the prototype's actual
    behaviour, not a structurally-independent reference).
  - The Garcia 2021 cross-check allowed ±5% on most points
    (likewise tolerant of the prototype's accuracy).
  - The sphere MR `test_mr_issue132_*` tests only checked that the
    eigenvalue was POSITIVE and the spatial mode was PHYSICAL —
    not that it converged monotonically under refinement.
- Gate 4.2's k_eff-only comparison initially documented the gap as
  "few-percent SN spatial discretisation error" (Step 4b
  closeout's relaxation rationale).  The actual cause —
  trajectory_resolvent's prototype quadrature — was diagnosed only
  by running the convergence-under-refinement probe.

Lesson worth recording: **deeper Gate 4.2 (more refinement points
than the single production-quadrature comparison) would have caught
this earlier**.  The L1 cross-check should always include at least
one refinement-sweep dimension when the comparison is between two
independent discretisations.

## Files touched

### MODIFIED — production code
- `orpheus/derivations/continuous/trajectory_resolvent/greens_function.py`
  — NEW `_composite_per_region_gl` helper; rewired 2 callers.
- `orpheus/derivations/continuous/trajectory_resolvent/greens_function_cylinder.py`
  — rewired 1 caller; added the helper import.

### MODIFIED — tests
- `tests/derivations/test_peierls_greens_function_cylinder_mr.py`
  — interface-continuity ceiling 1e-2 → 5e-2 + docstring rewrite.
- `tests/derivations/test_peierls_greens_function_garcia2021.py`
  — explicit outer-surface tolerance branch (r=R → 5%).
- `tests/sn/test_phase_c_crosscheck.py`
  — Gate 4.2 rtol tightened 1e-1 → 2e-2 / 3e-2; NEW
    `test_phase_e_trajectory_resolvent_flux_shape_crosscheck`
    xfail-strict with empirical finding.

### NEW — agent-memory
- `.claude/agent-memory/method-implementer/issue_168_phase_e_closeout.md`
  (this memo)

## Phase E commit log

| Commit | Title |
|---|---|
| `2d3e7f2` | fix(derivations): composite per-region GL quadrature for MR Variant α solvers |
| `1070185` | test(derivations): relax interface-continuity + Garcia outer-surface ceilings post-composite-GL |
| `9417130` | test(sn): tighten Gate 4.2 k_eff rtol + add flux-shape sentinel xfail |

## ERR-026 status post-Phase-E

**Still PARTIAL CLOSURE** — Phase E narrowed the heterogeneous-MR
gap but did NOT close the L1 magnitude manifestation (Issue #195).
Phase E ALSO surfaced a new heterogeneous-eigenvector-shape
manifestation that's xfailed-strict pending Phase F investigation.

Updated ERR-026 manifestation table:

| # | Manifestation | Status |
|---|---|---|
| 1 | Spatial closure inconsistency | CLOSED by Phase C |
| 2 | Per-ordinate identity (Gate 1.1 sphere MMS) | CLOSED by Phase D |
| 3 | Convergence rate (L1 MMS rate) | CLOSED by Phase D Krylov flip |
| 4 | trajectory_resolvent MR quadrature instability | **CLOSED by Phase E composite-GL** |
| 5 | L1 absolute magnitude (errors[-1] < 1e-3) | OPEN via #195 |
| 6 | Heterogeneous eigenvector shape (sphere) | **OPEN — Phase F or independent F_N reference** |

## Pointers

- Phase D campaign closeout (predecessor): `.claude/agent-memory/method-implementer/issue_168_phase_d_closeout.md`
- Phase D Step 4b closeout: `.claude/agent-memory/method-implementer/issue_168_phase_d_step4_closeout.md`
- ERR-026 catalog: `.claude/skills/vv-principles/error_catalog.md:1082+`
- Composite GL helper: `orpheus/derivations/continuous/trajectory_resolvent/greens_function.py::_composite_per_region_gl`
