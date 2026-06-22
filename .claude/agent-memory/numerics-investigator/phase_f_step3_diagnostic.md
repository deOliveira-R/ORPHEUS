---
name: Phase F Step 3 diagnostic — SI/sweep path is missing the Phase D Carlson seed; apply-matvec already has it
description: Numerics-investigator Step 3a-3c of Issue #168 Phase F. Diagnosis pins the SN sphere boundary-cell defect to `_sweep_1d_spherical`'s hardcoded `psi_angle = np.zeros((nx, ng))` at sweep.py:474. The Phase D Carlson coupled-pole seed fix was applied only to the apply-matvec (operator.py path); the SI/sweep path is missing the backport. Krylov inner-solver path (which goes through the apply-matvec) converges cleanly under refinement (r01: 1.029 → 1.001 at n=40→160) while SI diverges (0.522 → 0.489). Fix is mechanical: backport CarlsonInwardSweep into the spherical sweep's angular-recurrence seed.
type: project
---

# Phase F Step 3 — Diagnostic deep-dive closeout

**Branch**: `refactor/sn-operator-algebra` 2026-05-12. Phase F Step 3
of `.claude/plans/issue_168_phase_f_curvilinear_boundary_eigenvector.md`.

**Predecessor**: Step 2 memo
`.claude/agent-memory/numerics-investigator/phase_f_step2_mesh_refinement.md`
(SN pole ratio DIVERGES under refinement; Variant α converges cleanly).

**Sister memo (Phase D analog)**:
`.claude/agent-memory/numerics-investigator/phase_d_gate_1_1_sphere_mms_diagnosis.md`
(the apply-matvec-path twin of this Step 3 finding; closed by Phase D).

## Headline finding

**The exact same bug Phase D diagnosed and fixed in the apply-matvec
path (`orpheus/sn/operator.py`) is STILL PRESENT in the SI/sweep path
(`orpheus/sn/sweep.py::_sweep_1d_spherical:474`).** Phase D fixed the
apply-matvec by introducing the Phase D Carlson coupled-pole seed
(`CarlsonInwardSweep`) consumed via `MorelMontryAngularSweep`. The
SI/sweep path was never updated; it still initialises the half-angle
face flux to zero — the same legacy "Phase B ZeroSeed" that Phase D
explicitly diagnosed as the wrong-fixed-point bug.

`solve_sn` defaults to `inner_solver="source_iteration"` for all
geometries (per `orpheus/sn/solver.py:828` and design note at
`solver.py:863-875`). The default-flip was intentional in Phase E
(keeps the 6 curvilinear regression snapshots bit-identical) but
its consequence — that production heterogeneous MR eigenvalue runs
go through the buggy SI path — was not surfaced.

This is the structural defect responsible for the Step 2 divergent
pole-cell behaviour. The fix is mechanical: backport the Phase D
Carlson seed into `_sweep_1d_spherical`.

## 1. Per-ordinate ψ at the sphere pole (Step 3b probe)

Script: `/tmp/phase_f_step3a_per_ordinate.py`. Runs the
`sphere_2g_3reg` n=40 snapshot through `solve_sn` (default SI inner
solver) and extracts `result.angular_flux` shape `(N_ord, nx, ny=1,
ng)`. Inspects per-ordinate ψ(μ) at the pole cell (i=0,
r_centre=0.025).

**Pomraning 1989** structural-singularity result: ψ(r=0, μ) should be
near-isotropic (the streaming operator's μ-dependence vanishes at
r=0). Deviation reveals per-ordinate sweep pathology.

Empirical (n=40, g=0, GL-8):

| Cell | r_centre | mean(ψ) | std(ψ)/mean | (max-min)/mean |
| --- | --- | --- | --- | --- |
| **i=0 (pole)** | 0.025 | 5.49e-2 | **0.520** | **1.389** |
| i=1            | 0.075 | 1.07e-1 | 0.460 | 1.253 |
| i=20 (mid)     | 1.025 | 3.70e-1 | 0.310 | 0.912 |
| i=N-2          | 1.925 | 5.39e-1 | 0.310 | 0.835 |
| **i=N-1 (outer)** | 1.975 | 4.95e-1 | 0.244 | 0.630 |

The pole cell has by far the SHARPEST μ-dependence in the whole
mesh — exactly the opposite of the Pomraning isotropy prediction.
ψ peaks at small |μ| (μ≈±0.18, ψ≈9.0e-02) and minimises at large
|μ| (μ=±0.96, ψ≈1.4e-02) — a 6.4× spread.

The Σ w_n ψ_n integration to `sf[0, 0]` is correct to machine
precision (-1.5e-12 absolute drift); the bug is unambiguously
upstream of `_scalar_flux_from_angular`. The per-ordinate ψ at i=0
is structurally wrong.

## 2. Isolation: SI vs Krylov inner solver (Step 3c probe)

Script: `/tmp/phase_f_step3c_si_vs_krylov.py`. Runs `solve_sn` twice
on the same `sphere_2g_3reg` n=40 problem — once with
`inner_solver="source_iteration"` (the default; sweep path), once
with `inner_solver="krylov"` (apply-matvec path).

Result:

| Probe | k_eff | sf[0]/sf[1] | sf[N-1]/sf[N-2] | ψ@i=0 cv |
| --- | --- | --- | --- | --- |
| SI (sweep) | 1.35781531 | **0.5223** | 0.8871 | 0.520 |
| Krylov (apply-matvec) | 1.38464040 | **1.0288** | 0.9745 | 0.445 |

The Krylov path **eliminates the pole anomaly**: ratio is 1.029
(consistent with a smooth eigenvector). The outer-cell ratio
likewise improves from 0.887 to 0.975. The k_eff also changes by
2.0% — confirming the global eigenmode reshapes when the pole
closure is correct.

Mesh-refinement under Krylov (Step 3c sanity check, script
`/tmp/phase_f_step3c_krylov_refinement.py`):

```
   n    k_eff(SI)  r01(SI)   k_eff(Kr)  r01(Kr)   rN-1/N-2(Kr)
  40   1.35781531   0.5223   1.38464040  1.0288    0.9745
  80   1.35766493   0.5058   1.38261730  1.0125    0.9906
 160   1.35762957   0.4888   1.38167934  1.0018    0.9958
```

Krylov path:
- r01: 1.0288 → 1.0018 (factor 16× under 4× refinement → **O(h²)**)
- rN-1/N-2: 0.9745 → 0.9958 (factor 6× under 4× refinement → ~O(h²·⁵))

SI path: r01 worsens monotonically as Step 2 reported.

**Krylov-path SN now matches Variant α's clean O(h^p) convergence**.
The two paths give different k_eff values (1.38168 Kr at n=160 vs
1.37893 Var-α at n_r=96), but both converge cleanly toward their
own fixed point in the right pattern. Variant α's residual 0.2%
gap is a separate Phase E composite-GL quadrature floor (see Step 2
§3 footnote) — not relevant to this Phase F structural fix.

## 3. Code-walk: WHERE the bug lives

### 3.1 SI/sweep dispatch chain

`solve_sn` (`orpheus/sn/solver.py:824`, default `inner_solver=
"source_iteration"`) →
`power_iteration` (`solver.py:899`) →
`SNSolver._solve_source_iteration` (`solver.py:374,409`) →
`transport_sweep` (`sweep.py:99-151`) →
**`_sweep_1d_spherical`** (`sweep.py:396-545`) ← bug lives here.

### 3.2 The bug line

`orpheus/sn/sweep.py:474`:

```python
# Angular "face flux" between successive ordinates: ψ_{n-1/2, i}
# Shape (nx, ng). Initialised to zero for the first ordinate (α_{1/2}=0).
psi_angle = np.zeros((nx, ng))
```

This is the **EXACT TWIN** of the Phase D legacy bug at
`pole_angular_closure.py:411`:

```python
psi_half_left = np.zeros((ng, nx), dtype=psi_level.dtype)
```

Both lines initialise the M-M angular recurrence's half-angle face
flux ψ_{1/2,i,g} to zero — the **wrong term initialization** Phase D
diagnosed (Mode 3 failure mode per
`numerical-bug-signatures/SKILL.md`).

The Phase D diagnostic memo
(`phase_d_gate_1_1_sphere_mms_diagnosis.md` §4) explicitly stated:

> Production file: `orpheus/sn/spatial/pole_angular_closure.py`
> Function: `_mm_weighted_angular_recurrence_single_level`
> Line: 411
> [...] The literal `psi_half_left = np.zeros(...)` is the wrong
> Carlson seed. Per Hébert §3.9.4 Eqs. (3.432)-(3.435), the correct
> seed is the inward μ = −1 spatial sweep's cell-centred output
> `φ̄_{1/2,i}` — which IS a function of (Σ_t, source moments, BC),
> NOT identically zero.

The identical statement applies, verbatim, to
`sweep.py:474`. The math is the same; the production code path is
different.

### 3.3 Why Phase D didn't fix this site

Phase D was scoped to the apply-matvec path because Gate 1.1 MMS
runs through the apply-matvec via `SNStreamingOperator.apply` (the
xfail-strict marker on
`tests/sn/test_phase_c_gates.py::test_apply_curvilinear_per_ordinate_flat_flux_residual`
parametrised over the M-M closure). Phase D's
`MorelMontryAngularSweep` Carlson-seed default (`geometry.py:188-191`)
closes Gate 1.1 perfectly for the apply-matvec — confirmed by
n=40 SI vs Krylov split here (Krylov route consumes the apply-
matvec with the Carlson default; pole closure clean).

The SI/sweep path uses the `cell_update` strategy
(`DiamondDifference`), NOT the `pole_angular_closure` strategy. The
two strategies are separate fields on `SNMesh`
(`geometry.py:142,188`). The sweep's per-cell recurrence (in
`diamond.py::_update_curvilinear`) does the same M-M weighted-DD
math as the apply-matvec's M-M recurrence, but consumes
`upstream_state.angular_upstream` from the sweep orchestrator's
`psi_angle[i]` buffer (which is the hardcoded zero seed at sweep
start).

Phase D's bit-identical contract intentionally preserved the SI/
sweep path's output (the 6 curvilinear regression snapshots are
SI-generated). The intentional preservation meant the structural
defect on heterogeneous MR eigenvectors stayed live in the
production default path.

## 4. Fix site — concrete recommendation for `method-implementer`

### 4.1 Primary fix: `orpheus/sn/sweep.py:474` (the spherical sweep)

The fix is mechanical and mirror-images Phase D's apply-matvec fix:

**Replace** the hardcoded zero initialisation at line 474:

```python
# CURRENT (BUG):
psi_angle = np.zeros((nx, ng))
```

**with** a Carlson coupled-pole inward sweep that produces the
cell-centred ψ_{1/2,i,g} per Hébert §3.9.4 Eqs. (3.432)-(3.435).
The simplest realisation reuses `CarlsonInwardSweep` from
`orpheus/sn/spatial/psi_half_angle_seed.py` — but that strategy's
`__call__` consumes a `psi_level: (ng, M, nx)` array (the current
input ψ at all ordinates of the level) which is NOT yet available
at sweep start before the first ordinate iteration.

Two implementation options for `method-implementer`:

**Option A — call Carlson at iteration k = 0** (cleanest, matches
Hébert §3.9.4 structurally):

The Carlson inward μ = −1 sweep is **source-driven**, not
ψ-current-iterate-driven, when read literally per Hébert Eq.
(3.432). It needs:
- The cell-averaged source `Q̄_i` evaluated at μ = −1, built from
  the source's Legendre moments via `(2ℓ+1)/2 · Q_ℓ · P_ℓ(−1)`.
- The outer-face inflow at μ = −1 (from BC trace).
- Σ_t and Δr per cell.

For an L=0 isotropic operator (current ORPHEUS), this reduces to
`Q̄_i = Q_iso(r_i) / 2` (the `weight_norm = 1/Σw = 1/2` factor at
`sweep.py:481` divides by 2 to convert isotropic-source-per-cell to
angular-source-density at each ordinate). For the SI inner-solver
loop, `Q_1d[i,g]` is the within-group source which already includes
scattering + fission + fixed contributions — so the Q̄_i for the
Carlson sweep is `Q_1d[i,g] / Σ w` (matches the existing
`weight_norm` machinery).

Concrete intervention at `sweep.py:474`:

```python
# Phase F fix: Carlson coupled-pole seed for ψ_{1/2,i,g}.
# Replaces the legacy hardcoded zero (Phase D analog for the
# SI/sweep path).
psi_angle = _carlson_inward_sweep_seed(
    Q_1d, sig_t_1d, sn_mesh.dx,
    bc_outer, bc_outer_obj, weights, mu, ng, nx,
)
```

with `_carlson_inward_sweep_seed` a local helper running Hébert
Eqs. (3.434)-(3.435):

```python
def _carlson_inward_sweep_seed(
    Q_1d, sig_t_1d, dr, bc_outer, bc_outer_obj, weights, mu, ng, nx,
):
    """Phase F: cell-centred ψ_{1/2,i,g} via Hébert §3.9.4 inward
    μ=-1 sweep. Backport of Phase D's CarlsonInwardSweep for the
    SI/sweep path.
    """
    # Outer-face inflow at μ=-1: apply BC to current outflow
    # estimate (initially zero on first SI iteration; the loop
    # updates bc_outer in-place after each outward sweep).
    inflow_full = bc_outer_obj.apply(bc_outer)  # (N, ng)
    most_inward_idx = int(np.argmin(mu))
    bc_outer_value = inflow_full[most_inward_idx, :]  # (ng,)

    # Cell-averaged source at μ=-1, L=0 (isotropic): Q̄_i = Q_1d[i]/Σw.
    # (Matches sweep.py:481-482 weight_norm; P_0(-1)=1; (2·0+1)/2=1/2;
    # factor 1/2 inherited from the weight_norm = 1/Σw = 1/2 for GL.)
    Q_bar = Q_1d * (1.0 / weights.sum())  # (nx, ng)

    # Hébert (3.434)-(3.435) inward DD recurrence, vectorised across g.
    phi_aux = np.zeros((nx, ng), dtype=Q_1d.dtype)
    phi_face = bc_outer_value.copy()  # (ng,)
    for k in range(nx - 1, -1, -1):
        denom = dr[k] * sig_t_1d[k] + 2.0     # (ng,)
        phi_cell = (dr[k] * Q_bar[k] + 2.0 * phi_face) / denom
        phi_aux[k] = phi_cell
        phi_face = 2.0 * phi_cell - phi_face
    return phi_aux
```

**Option B — reuse `CarlsonInwardSweep`** directly by adapting its
input contract:

Since `CarlsonInwardSweep.__call__(psi_level, context)` consumes a
ψ-current array shaped `(ng, M, nx)` to build the L=0 moment, and
the SI/sweep path has `Q_1d` (the source) and not a current ψ at
sweep start, the contract mismatch makes direct reuse awkward.
However, after sweep iteration k > 0 the previous iteration's
angular_flux exists — `CarlsonInwardSweep` could be invoked at the
END of iteration k to feed iteration k+1. This is more invasive
(requires looping at the iteration level, not the sweep level) and
is NOT recommended.

**Method-implementer's choice**: Option A is the correct fix.
Option B's source-vs-ψ contract mismatch is exactly why Phase D's
sister apply-matvec path uses `psi_level` (the operator's input ψ),
while the SI/sweep path naturally uses `Q_1d` (the within-group
source). Both are mathematically equivalent on the fixed-point
solution where `L · ψ = Σ_t · ψ + (scattering + fission)` ⇒
`Σ_t · φ_0(r) = Q̄(r)` and the inward sweep gives the same answer
either way.

### 4.2 Cylindrical analogue: `sweep.py::_sweep_1d_cylindrical`

The cylindrical sweep is at `sweep.py:552-...` and has its own
per-μ-level loop. Per the Phase D diagnostic §7.3 architectural
observation:

> The M-M recurrence currently uses `psi_half_left = 0` for the
> CYLINDRICAL case too, but cylindrical Gate 1.1 passes empirically
> [...because] each μ-level's α-dome ends at α=0 by antisymmetry,
> so the wrong `psi_half_left = 0` cancels at the level boundary.

The cylindrical sweep is the same structural class as the spherical
one. Phase F's fix should also apply the per-level Carlson seed to
the cylindrical sweep for architectural alignment — and verify the
6 cylindrical regression snapshots still pass under the fix. The
per-level cylindrical context already exists in the apply-matvec
side (`operator.py::transport_operator_matvec_cylindrical` builds
`carlson_context: list[CarlsonSweepContext]`); the SI/sweep path
needs the analog.

Cylindrical regression behaviour is currently O(h^p) clean (Step 2
did NOT identify a divergent ratio for cylinder) — but architectural
alignment is the Cardinal Rule 2 imperative. The cylindrical fix is
LOWER URGENCY than the spherical fix (the bug there is invisible),
but should ship in the same Phase F change.

### 4.3 Bit-identity break — expected and unavoidable

The 6 curvilinear regression snapshots (`tests/sn/regression/snapshots/`)
were SI-generated. The Carlson fix in the sweep path will change
their values. The regression contract per `vv-principles` §"Bit-
identity vs principled-equivalence" requires:

1. **Principled new formulation**: yes — Hébert §3.9.4 Carlson seed
   is the canonical mathematical form, replacing a hardcoded zero
   that Phase D explicitly diagnosed as wrong. Each intermediate
   is named (`phi_aux[i,g]` ↔ Hébert `φ̄_{1/2,i}`).
2. **Structurally-independent reference**: yes — Variant α
   (`solve_greens_function_sphere_mr`) is a Schur + power-iteration
   solver on the integral form, structurally independent of SN's
   per-cell sweep. Step 2 empirical evidence shows Variant α
   converges to clean O(h^p) ratios at all probes; the fixed SN
   will match Variant α's qualitative shape.
3. **FP-non-associativity bound**: this is NOT a non-associativity
   change — the Carlson seed is a different VALUE, not a different
   reduction order. The drift is algorithmic (intended) and
   well-defined (Hébert §3.9.4 is the closed-form ground).

The snapshot regeneration follows the Phase D campaign precedent
(Phase D regenerated curvilinear snapshots in
`issue_168_phase_d_step3_closeout.md`).

### 4.4 Important: do NOT apply this fix to `_sweep_1d_cumprod` or `_sweep_2d_*`

The fix is curvilinear-specific. The Cartesian sweep paths
(`_sweep_1d_cumprod`, `_sweep_2d_*`) have NO angular redistribution
term, NO α-cascade, NO half-angle face flux ψ_{1/2,i,g}. They
must NOT be touched.

## 5. Confidence: NEAR-CERTAIN

The diagnosis is closed on three independent pieces of empirical
evidence:

1. **Step 2 mesh-refinement DIVERGENCE** (`phase_f_step2_mesh_refinement.md`):
   SI ratio 0.522 → 0.473 with refinement vs Variant α's clean
   convergence to 1. Structural defect proven.
2. **Step 3a per-ordinate isotropy violation**: ψ at i=0 has cv=0.52
   (max/min=6.4×). Pomraning 1989 predicts near-isotropic. Per-
   ordinate sweep is structurally wrong upstream of `_scalar_flux_from_angular`.
3. **Step 3c SI-vs-Krylov split**: switching to `inner_solver="krylov"`
   (which goes through the apply-matvec with the Phase D Carlson
   seed) ELIMINATES the pole anomaly (r01: 1.029 → 1.001 under
   refinement). Same problem, same solver, same materials, same
   quadrature, same mesh — ONLY the inner-solver dispatch changes.
   That isolates the bug to the SI/sweep path's M-M recurrence
   seed.

The line `sweep.py:474 psi_angle = np.zeros((nx, ng))` is the
twin of the Phase D legacy bug at `pole_angular_closure.py:411`
(which Phase D explicitly diagnosed and fixed via `CarlsonInwardSweep`).
The fix is mechanical: backport the Carlson coupled-pole sweep into
the SI/sweep path.

The remaining ~ 0.18% k_eff gap between Krylov SN (n=160) and
Variant α (n_r=96) at 1.38168 vs 1.37893 is the Phase E composite-
GL quadrature floor identified in Step 2 §3 — a separate, milder
issue NOT part of Phase F's scope.

## 6. Verification path for `method-implementer`

After implementing the fix in `_sweep_1d_spherical`:

1. **Re-run the Step 3c SI-vs-Krylov probe**: SI must now produce
   the same eigenvalue and r01 as Krylov did pre-fix (1.029 at
   n=40, 1.0018 at n=160). Diagnostic script:
   `/tmp/phase_f_step3c_si_vs_krylov.py`.

2. **Re-run Step 3a per-ordinate probe**: ψ at i=0 cv must drop
   from 0.520 to ~0.445 (Krylov post-fix value). Pomraning
   isotropy improves but does NOT reach machine zero at finite h
   — the discrete pole-cell value reflects ψ at r=0.025, not at
   r=0. Diagnostic script:
   `/tmp/phase_f_step3a_per_ordinate.py`.

3. **Phase D Gate 1.1 (apply-matvec) regression**: must STAY GREEN.
   The fix is in the SI/sweep path; the apply-matvec was not
   touched. Run:
   `pytest tests/sn/test_phase_c_gates.py::test_apply_curvilinear_per_ordinate_flat_flux_residual -v`.

4. **NEW Gate 1.6 (sweep-path per-ordinate flat-flux residual)**:
   the dual of Gate 1.1 — flat-ψ probe through the sweep path
   (not the apply-matvec) must satisfy `L·ψ = Σ_t·ψ` per ordinate
   to ~1e-12. Promote this from a diagnostic script to a
   permanent test at `tests/sn/spatial/test_sweep_carlson_seed.py`
   or similar.

5. **Phase E flux-shape sentinel**:
   `tests/sn/test_phase_c_crosscheck.py::test_phase_e_trajectory_resolvent_flux_shape_crosscheck`
   should xpass after the fix; remove the xfail marker.

6. **Curvilinear regression snapshots**: 6 curvilinear snapshots
   will need regeneration:
   - sphere_homogeneous_1g_n40
   - sphere_homogeneous_2g_n40
   - sphere_2g_3reg_dd_n40
   - cylinder_homogeneous_1g_ls4
   - cylinder_homogeneous_2g_ls4
   - cylinder_2g_3reg_ls4
   Generate via `python -m tests.sn.regression._generate_snapshots
   --case ...` (or the equivalent CLI; the framework already
   regenerates on `pytest tests/sn/regression -v --regen`).
   Verify each regenerated snapshot passes Gate 4.2 k_eff
   tolerance AND the new flux-shape sentinel.

7. **Apply-vs-sweep consistency test**: Phase C explicitly aligned
   the apply matvec and sweep math. After the Phase F fix, the
   two paths should produce the SAME eigenvalue and the SAME
   flux profile (within iteration tolerance). The current SI-vs-
   Krylov gap (k_eff 1.358 vs 1.385 at n=40) is the *manifestation*
   of the bug — post-fix, the SI value must match Krylov's.
   This becomes a new structural invariant test:
   `tests/sn/spatial/test_sweep_vs_apply_consistency.py`.

## 7. Files referenced

### Created (not committed)

- `/tmp/phase_f_step3a_per_ordinate.py` — Step 3b per-ordinate ψ probe
- `/tmp/phase_f_step3c_si_vs_krylov.py` — Step 3c SI vs Krylov split
- `/tmp/phase_f_step3c_krylov_refinement.py` — Krylov refinement sanity
- `/tmp/phase_f_step3a_psi_n40.npz` — saved ψ for inspection

### Read

- `orpheus/sn/sweep.py:393-545` — `_sweep_1d_spherical` (bug site at 474)
- `orpheus/sn/sweep.py:99-151` — `transport_sweep` dispatch
- `orpheus/sn/solver.py:824-920` — `solve_sn` (defaults to SI)
- `orpheus/sn/operator.py:680-840` — `transport_operator_matvec_spherical` (has Carlson seed)
- `orpheus/sn/spatial/pole_angular_closure.py` — the Phase D `MorelMontryAngularSweep` + `_mm_weighted_angular_recurrence_single_level`
- `orpheus/sn/spatial/psi_half_angle_seed.py` — the `CarlsonInwardSweep` / `ZeroSeed` strategies
- `orpheus/sn/spatial/diamond.py:550-624` — `_update_curvilinear` (per-cell recurrence the sweep dispatches into)
- `orpheus/sn/geometry.py:120-200` — `SNMesh.__init__` (defaults `pole_angular_closure`)

### Sister artefacts

- `.claude/agent-memory/numerics-investigator/phase_f_step2_mesh_refinement.md` — Phase F Step 2
- `.claude/agent-memory/numerics-investigator/phase_d_gate_1_1_sphere_mms_diagnosis.md` — Phase D analog
- `.claude/agent-memory/literature-researcher/phase_d_carlson_coupled_pole.md` — Hébert §3.9.4 math
- `.claude/plans/issue_168_phase_f_curvilinear_boundary_eigenvector.md` — Phase F plan

## 8. Failure-mode characterisation per `vv-principles`

This bug is **Mode 3 — Missing factor / wrong term initialization**
(identical to Phase D's ERR-026 ancestor). The hardcoded `psi_angle
= 0` is a wrong term initialization with 1-group/homogeneous-flat
invisibility.

Why it survived Phase D's regression suite:
- Phase D fixed the apply-matvec path; the SI/sweep path's
  `psi_angle = 0` was untouched.
- The 6 curvilinear regression snapshots were SI-generated under
  the wrong seed; they "pass" by tautology (the bug is baked into
  the contract).
- Phase D's Gate 1.1 MMS xfail-strict marker is on the apply-matvec
  path; the SI/sweep path doesn't run that probe.
- 1-group homogeneous reflective passes by degeneracy
  (Sigma_a-only k_eff; flat eigenmode masks any seed bug).
- Multi-region heterogeneous WAS the only stressor — and Phase E's
  Variant α comparison was xfailed for flux shape, not eigenvalue.
  The k_eff stayed within 2e-4 of Variant α at n=40 (Phase E
  achievement), masking the structural shape defect.

This is a **double Mode 3** — same anti-pattern, same fix, two
production-code paths. Phase D closed one; Phase F closes the
other.

## 9. Recommended ERR-026 manifestation update

After Phase F ships:

- ERR-026 manifestation table row #6 ("heterogeneous eigenvector
  shape on SI/sweep path") flips OPEN → CLOSED.
- Add a NEW ERR-NNN entry for "Phase D fix not backported to
  SI/sweep path until Phase F" — OR fold this into the existing
  ERR-026 entry's lessons (the latter is simpler; the manifestation
  is the same bug class).
- Add a Phase F-discovered anti-pattern to `vv-principles`
  `SKILL.md` §Anti-patterns:
  > **Whenever a fix is applied to one of two structurally-mirrored
  > production paths (apply-matvec vs SI/sweep, etc.), MUST check
  > the OTHER path for the same defect.** Mode 3 wrong-term-
  > initialization defects often appear in pairs; fixing one path
  > without auditing its sister is a Cardinal Rule 2 (architecture)
  > violation that ERR-026 instantiated twice.

## Pointers

- **Diagnostic scripts** (this Step 3): `/tmp/phase_f_step3a_per_ordinate.py`,
  `/tmp/phase_f_step3c_si_vs_krylov.py`,
  `/tmp/phase_f_step3c_krylov_refinement.py`.
- **Step 2 memo**: same directory `phase_f_step2_mesh_refinement.md`.
- **Phase D analog memo**: `phase_d_gate_1_1_sphere_mms_diagnosis.md`.
- **Phase F plan**: `.claude/plans/issue_168_phase_f_curvilinear_boundary_eigenvector.md`.
- **Bug site**: `orpheus/sn/sweep.py:474` (and cylindrical analog
  in `_sweep_1d_cylindrical`).
- **Existing canonical fix to backport**:
  `orpheus/sn/spatial/psi_half_angle_seed.py::CarlsonInwardSweep`
  (Hébert §3.9.4 Eqs. 3.432-3.435).
- **`vv-principles` reference**: Mode 3 (missing factor / wrong
  term initialization) + Anti-pattern: "twin-path fix
  incompleteness".
