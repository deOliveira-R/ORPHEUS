# Phase F — SN curvilinear boundary-cell eigenvector defect

**Tracking issue**: Continuation of Issue [#168](https://github.com/deOliveira-R/ORPHEUS/issues/168) (ERR-026). Closes the heterogeneous-eigenvector-shape manifestation that Phase E surfaced.
**Branch**: `refactor/sn-operator-algebra` (continue from Phase E tip `6708a4a`)
**Predecessor closeouts**:
- `.claude/agent-memory/method-implementer/issue_168_phase_e_closeout.md` (immediate prior, names the empirical evidence)
- `.claude/agent-memory/method-implementer/issue_168_phase_d_closeout.md` (Carlson coupled-pole seed)
- `.claude/agent-memory/numerics-investigator/phase_d_gate_1_1_sphere_mms_diagnosis.md` (the prior numerics-investigator probe pattern this work will extend)

---

## Context (load into a fresh-context Claude)

### The problem the prior phases left open

Phase E (commits `2d3e7f2..6708a4a`, 2026-05-12) shipped the composite per-region GL correction to trajectory_resolvent's MR Variant α solvers, dropping the SN-vs-Variant-α heterogeneous-MR k_eff gap from 7–9% to **2e-4 (sphere) / 1.75e-2 (cylinder)**. As a by-product Phase E built a flux-shape comparison harness and discovered a separate, unresolved issue:

**The two solvers AGREE on the eigenVALUE but DISAGREE on the eigenVECTOR shape.**

Empirical evidence captured in
[`tests/sn/test_phase_c_crosscheck.py::test_phase_e_trajectory_resolvent_flux_shape_crosscheck`](../../tests/sn/test_phase_c_crosscheck.py) (xfail-strict):

| | SN scalar_flux range | Variant α phi_g range | k_eff rel |
|---|---|---|---|
| Sphere 2g 3reg | 0.132 → 1.20 (9.1×) | 0.128 → 0.187 (1.7×) | 2e-4 |
| Cyl 2g 3reg LS4 | (similar pattern, smaller magnitude) | … | 1.75e-2 |

### The diagnostic pattern that points at the root cause

When I inspected the SN sphere `scalar_flux[g=0]` for the `sphere_2g_3reg_dd_n40` snapshot:

```
0.132 0.253 0.292 0.324 0.349 0.370 0.388 0.403 0.415 0.424   ← region 0 fuel A
0.474 0.561 0.629 0.678 0.716 0.746 0.769 0.787 0.800 0.810
0.817 0.822 0.825 0.828 0.829 0.831 0.833 0.836 0.840 0.848   ← region 1 mod B
0.861 0.878 0.899 0.934 0.987 1.058 1.135 1.198 1.202 1.067   ← region 2 fuel A
```

The middle of the profile (cells 1–38) looks like a plausible fuel-mod-fuel reflective eigenmode for the asymmetric A|B|A layout. The anomaly lives **only at the two boundary cells**:

- **Cell 0** (r=0.025, the spherical-pole cell): 0.132 — drops **1.9×** to its neighbor 0.253. Every other adjacent-cell ratio is ≤ 1.15×.
- **Cell 39** (r=1.975, the outer-reflective cell): 1.067 — drops FROM cell 38's peak 1.202, even though ∂φ/∂r=0 at r=R should make r=R the local MAX.

Variant α's `phi_g[g=0]` on its composite-GL nodes shows no analogous anomaly — the innermost and outermost GL nodes' values flow smoothly with their neighbors.

### Best-guess hypothesis

A residual SN integration / cell-averaging defect at curvilinear **boundary cells** where:

- The **per-ordinate angular flux** ψ_{n, i=0, g} is `correct-on-flat-ψ` (Phase D's Carlson coupled-pole seed pinned this via Gate 1.1) but…
- The **angle-integrated scalar flux** `φ_{i=0, g} = Σ w_n ψ_{n, i=0, g}` on **non-flat eigenvectors** has a residual bias.
- The same defect at the outer reflective face (i = nx-1) explains why cell 39 sits below cell 38 despite the reflective BC.

Phase D's Carlson seed fix is a NECESSARY but not SUFFICIENT condition for the boundary-cell scalar flux to be correct. The seed pins per-ordinate identity on FLAT ψ probes (Gate 1.1); it does not directly pin the **angle-integrated** scalar flux on **non-flat** eigenvectors.

This hypothesis is the Phase F STARTING POINT, not its conclusion. Steps 1–3 below empirically falsify or confirm it.

---

## Goal

After Phase F ships, ALL of these hold simultaneously:

1. `tests/sn/test_phase_c_crosscheck.py::test_phase_e_trajectory_resolvent_flux_shape_crosscheck` xpasses → marker removed (xfail-strict will self-enforce removal).
2. SN and Variant α agree on the EIGENVECTOR SHAPE on heterogeneous MR to within the chosen Gate-4.2-shape tolerance (target: ≤ 5% per-cell after L∞ normalisation; defended empirically).
3. `error_catalog.md` ERR-026 manifestation table: #6 (heterogeneous eigenvector shape) flips OPEN → CLOSED. ERR-026 remains PARTIAL only on #5 (L1 absolute magnitude per #195).
4. The boundary-cell defect (if real) is **fixed in production code**, not patched around in tests. Tests pin the fix; the fix lives in `orpheus/sn/`.

If after Step 2 the empirical evidence shows the anomaly converges out under mesh refinement (O(h^p)), the deliverable shifts: file a follow-up issue for the convergence-rate of the boundary-cell error, and tighten the existing flux-shape comparison to a refinement-dependent bound. ERR-026 #6 stays OPEN in that case (it would become a CONVERGENCE manifestation, similar to #5).

---

## Step 1 — Snapshot freshness check (~10 min)

**Goal**: confirm the stored `sphere_2g_3reg_dd_n40` snapshot was actually generated under the post-Phase-D defaults. If pre-Phase-D, regenerate and check whether the stored 9× profile is an obsolete artefact.

### What to do

1. **Inspect snapshot metadata**:

   ```bash
   .venv/bin/python -c "
   import numpy as np
   d = np.load('tests/sn/regression/snapshots/sphere_2g_3reg_dd_n40.npz')
   for k in d.keys():
       v = d[k]
       if v.shape == ():
           print(f'  {k}: {v.item()}')
   "
   ```

   Read the `generator_commit` field. Compare to:
   - Phase D commits `9512459..7d1cdd8` (Phase D shipping window — the snapshot SHOULD have been regen'd at or after `c44fe9b` if Phase D's regeneration trigger fired).
   - Phase E commits `2d3e7f2..6708a4a` (Phase E did NOT touch the SN side, so didn't regen).

2. **If pre-Phase-D**: the snapshot is stale. Regenerate:

   ```bash
   .venv/bin/python -m tests.sn.regression._generate_snapshots --case sphere_2g_3reg_dd_n40
   ```

   (Check `_generate_snapshots.py` for the actual CLI — may need to call a specific function. If the module has no CLI, write a one-off script that imports `_sphere_3region` and calls the snapshot-writer.)

   Save the regenerated `scalar_flux` separately, then compare to the stored values:

   ```bash
   .venv/bin/python -c "
   import numpy as np
   stored = np.load('tests/sn/regression/snapshots/sphere_2g_3reg_dd_n40.npz')['scalar_flux']
   # ... regenerate locally with current code, call fresh_sf
   # Compare:
   print(f'max abs diff = {np.max(np.abs(stored - fresh_sf))}')
   print(f'max rel diff = {np.max(np.abs(stored - fresh_sf) / np.abs(stored))}')
   "
   ```

3. **If post-Phase-D** (most likely — Phase D Step 3 closeout claims snapshots are bit-identical because the SI/sweep path generates them, not the apply matvec): proceed to Step 2. The stored profile IS what the current SN solver produces.

### Empirical decision

- **Stored snapshot is post-Phase-D AND regenerated values match bit-identically**: confirms the 9× profile is the current SN solver's output. Proceed to Step 2.
- **Stored snapshot is pre-Phase-D AND regen produces a DIFFERENT (smoother) profile**: the issue was a stale snapshot. Update the snapshot, re-run the Phase E flux-shape sentinel — if it now passes, remove the xfail marker and ERR-026 #6 closes. If it still fails, proceed to Step 2.
- **Stored snapshot is post-Phase-D AND regen produces a different profile** (unexpected): something non-deterministic in the SN solver. File a follow-up issue + investigate before continuing.

### Deliverable

A short note (maybe inline in this plan or in a step-1 memo) recording the `generator_commit` and the regen comparison. No production-code changes expected in Step 1.

---

## Step 2 — Mesh-refinement convergence study (~30 min wall-time + analysis)

**Goal**: determine whether the cell-0 / cell-N-1 anomalies are O(h^p) discretisation errors that converge out, or structural defects that persist at all refinements.

### What to do

1. **Run SN at multiple mesh refinements** on the same heterogeneous MR sphere problem. Build a probe at `/tmp/phase_f_sn_mesh_refinement.py`:

   ```python
   import numpy as np
   from orpheus.derivations.common.xs_library import get_xs, get_mixture
   from orpheus.geometry.structured_geometry import StructuredGeometry, Region, BC
   from orpheus.geometry.mesh import Mesh1D, RegionMesh
   from orpheus.sn.quadrature import GaussLegendre1D
   from orpheus.sn.solver import solve_sn_eigenvalue   # or solve_sn — check the API

   # Replicate _sphere_3region for n in {40, 80, 160, 320}
   # Use the SAME XS data and quadrature as the snapshot generator
   fuel = get_mixture("A", "2g")
   mod = get_mixture("B", "2g")

   for n_total in [40, 80, 160, 320]:
       n_per_region = (n_total // 4, n_total // 2, n_total // 4)
       geom = StructuredGeometry(
           geometry="SPH",
           regions=(
               Region(mat_id=0, outer_thickness_cm=0.5),
               Region(mat_id=1, outer_thickness_cm=1.0),
               Region(mat_id=0, outer_thickness_cm=0.5),
           ),
           bcs=(BC.reflective,),
       )
       mesh = Mesh1D.from_geometry(
           geom,
           region_meshes=tuple(RegionMesh(n_cells=n) for n in n_per_region),
       )
       quadrature = GaussLegendre1D.create(n_ordinates=8)
       result = solve_sn_eigenvalue(
           materials={0: fuel, 1: mod}, mesh=mesh, quadrature=quadrature,
           scattering_order=0,
       )
       sf = result.scalar_flux[:, 0, :]   # (n_total, 2)
       print(f"n={n_total}  k_eff={result.k_eff:.10f}")
       print(f"  cell 0:    sf[g=0]={sf[0, 0]:.6e}  sf[g=1]={sf[0, 1]:.6e}")
       print(f"  cell 1:    sf[g=0]={sf[1, 0]:.6e}  sf[g=1]={sf[1, 1]:.6e}")
       print(f"  ratio 0/1: {sf[0, 0] / sf[1, 0]:.4f}  (expected ~1)")
       print(f"  cell N-1:  sf[g=0]={sf[-1, 0]:.6e}")
       print(f"  cell N-2:  sf[g=0]={sf[-2, 0]:.6e}")
       print(f"  ratio N-1/N-2: {sf[-1, 0] / sf[-2, 0]:.4f}  (expected ~1 for reflective)")
   ```

   (Sanity check API names against `orpheus/sn/solver.py` first — `solve_sn_eigenvalue` may be `solve_sn` or `solve_sn_keff` or similar.)

2. **Tabulate the convergence**: collect `cell-0/cell-1` and `cell-(N-1)/cell-(N-2)` ratios as a function of `n_total`. If the ratios approach 1 as `n_total → ∞`, the anomaly is O(h^p). If they stay fixed at ~0.5 (cell 0/1) and ~0.89 (cell N-1/N-2), the anomaly is structural.

3. **Also run the Variant α at matching refinement** — `_run_sphere_2g_3reg_full()` already exists in the test file; just adapt and call at production quadrature. Compute the same boundary-cell ratios on Variant α's `phi_g`. If Variant α stays smooth (ratios ~1) at all quadrature levels, the SN is the outlier.

### Empirical decision

| SN cell-0/cell-1 ratio under refinement | Interpretation | Next step |
|---|---|---|
| Approaches 1 as `n → ∞`, with rate ~h¹ or h² | O(h^p) discretisation error at boundary cells. Bug is real but converges out. | File a follow-up issue tracking the convergence rate (not blocking); update Phase E flux-shape sentinel to use a refinement-aware tolerance; **skip Step 3**. |
| Stays fixed at ~0.5 across refinements | Structural defect that does NOT converge out. The 9× profile shape is a real bug. | Proceed to Step 3 with the diagnosis narrowed: scalar-flux integration at `r_{i-1/2} = 0` cells is the suspect. |
| Goes the WRONG way (gets worse with refinement) | Divergent behaviour at boundary cells. Very likely a structural bug worse than expected. | Proceed to Step 3 with **high urgency** — file a fresh ERR-NNN entry in `error_catalog.md`. |
| Variant α also shows anomalies at boundary cells at high refinement | Both methods have boundary issues; one of them is wrong about the model (BC interpretation). | Need a 3rd reference (F_N or independent literature value) — escalate to user; **don't proceed to Step 3** without it. |

### Deliverable

- Probe script at `/tmp/phase_f_sn_mesh_refinement.py` (one-off, not committed; results pasted into the closeout memo).
- A short summary table in the Step 2 closeout note (file: `.claude/agent-memory/numerics-investigator/phase_f_step2_mesh_refinement.md`).
- Decision recorded explicitly: which branch of the table fires, and why.

---

## Step 3 — Deep audit of SN spherical pole / outer-face scalar flux

**Goal**: identify the specific code site where the boundary-cell scalar flux acquires its anomalous value, propose a fix, validate.

### Triggered only if Step 2 confirms structural defect

If Step 2 shows the anomaly converges out, **skip this step**.

### Sub-steps

#### 3a. Trace `scalar_flux[i=0, g]` from per-ordinate ψ

Per `orpheus/sn/solver.py:655` `_scalar_flux_from_angular`:

```python
sf[ix, iy, :] = np.sum(fi[:, :, ix, iy] * quad.weights[None, :], axis=1)
```

This is `Σ w_n ψ_n` — should be correct for any geometry. But for sphere, ψ_n at i=0 is computed by the **M-M angular recurrence + per-cell balance** which has special handling at the pole. Trace:

- Where is `fi[g, n, i=0, 0]` (ψ_n at the i=0 cell-centre) computed in the spherical curvilinear sweep? Likely `orpheus/sn/sweep.py::_sweep_1d_spherical` or `orpheus/sn/operator.py::transport_operator_matvec_spherical`.
- For each outgoing ordinate `n` with μ_n > 0, the sweep evaluates ψ_{n, i=0} using the WDD spatial diamond relation `ψ_{n, face_out} = 2·ψ_{n, cell} - ψ_{n, face_in}` AND the M-M angular recurrence `ψ_half_left → ψ_half_right`. Verify both formulas at i=0 produce the right cell-centred value.
- For each incoming ordinate `n` with μ_n < 0, the sweep at the spherical pole uses the reflected partner. Is the reflection correct?

#### 3b. Compare per-ordinate flux at i=0 vs i=1 on the eigenvector

Build a probe that runs SN on the heterogeneous MR sphere AND extracts per-cell, per-ordinate ψ:

```python
# After solve_sn_eigenvalue, instead of using the convenience scalar_flux,
# dig into the result for the angular flux:
psi = result.angular_flux  # (ng, N_ord, nx, 1) — check actual attribute name
print("Per-ordinate ψ at i=0 (pole):")
for n in range(quad.N):
    print(f"  n={n}  μ_n={quad.mu_x[n]:+.4f}  ψ={psi[0, n, 0, 0]:.6e}")
print("Per-ordinate ψ at i=1:")
for n in range(quad.N):
    print(f"  n={n}  μ_n={quad.mu_x[n]:+.4f}  ψ={psi[0, n, 1, 0]:.6e}")
```

Symmetry check at the pole: ψ(r=0, μ) should be approximately isotropic (Pomraning 1989 — the streaming operator's μ-dependence vanishes at r=0). So at i=0 (r=0.025 in n=40 mesh), ψ should be APPROXIMATELY constant across all ordinates.

If ψ at i=0 is strongly μ-dependent → the per-ordinate values have a structural issue. The Σ w_n ψ_n integration is then "correct integration of a wrong angular shape".

#### 3c. Audit the reflective BC handling at i=N-1

For the outermost cell with reflective BC at r=R:
- The BC trace law (Phase C / §16A.3) maps outgoing ordinates to incoming ordinates by reflecting μ → -μ. The Phase D Carlson seed provides bc_outer_value for the auxiliary inward ordinate; the actual BC apply at the boundary edge handles all incoming ordinates from the outflowing WDD-propagated face values.
- Audit: at i=N-1, the cell-centred ψ_n for incoming ordinates (μ_n < 0) is computed from the WDD spatial relation using the reflected outflow. Does this match what we'd expect for a reflective BC at i=N-1?
- Specific check: for a flat eigenvector ψ_n(r) = constant, the reflective BC at i=N-1 should produce a self-consistent solution. Phase D's Gate 1.1 already pins this for flat-ψ. But for NON-FLAT ψ, is the cell-N-1 scalar flux consistent with the surrounding cells?

#### 3d. Propose the fix

Likely fix sites:
- The cell-balance equation at i=0 for the M-M angular recurrence: verify the cell volume `V_0 = (4/3)π·r_{1/2}^3` is correctly applied. If the sweep uses a cell-INTEGRATED relation (with V_0 multiplying both sides), check that the eventual cell-AVERAGED value is correctly derived.
- The pole-face IC at i=0 for the WDD outward sweep (`psi_face_in = fi[:, outgoing_mask, i0, 0].copy()`, Lewis-Miller §4.5 Carlson seed per Phase C). This is currently the cell-centred value of the pole cell; verify it's appropriate for non-flat eigenvectors.
- The outer-face BC trace at i=N-1: verify the reflective mapping produces the right incoming ψ_n for non-flat eigenvectors.

#### 3e. Implement fix + verify

Once a fix site is identified:
1. Implement the fix in `orpheus/sn/`. Likely in `operator.py` or `sweep.py`.
2. Run the Gate 1.1 test suite to confirm no regression (the per-ordinate flat-flux identity must still hold).
3. Re-run the heterogeneous MR snapshot to confirm the new SN profile has smooth boundary cells.
4. Re-enable `test_phase_e_trajectory_resolvent_flux_shape_crosscheck` (remove the xfail marker). It should now xpass.
5. Re-run the SN regression suite — the 6 curvilinear snapshots will need to be regenerated under the fix. Verify the regenerated snapshots also pass Gate 4.2 + the flux-shape check.

### Deliverable

- Production-code fix in `orpheus/sn/`.
- New foundation test pinning the boundary-cell scalar flux invariant (the dual of Gate 1.1's per-ordinate identity but for the angle-integrated scalar flux).
- Regenerated curvilinear regression snapshots.
- Phase F closeout memo at `.claude/agent-memory/method-implementer/issue_168_phase_f_closeout.md`.

---

## Reference materials (load these into the post-compaction agent)

### Closeout memos (predecessor phases)

- `.claude/agent-memory/method-implementer/issue_168_phase_e_closeout.md` — Phase E (composite-GL + flux-shape sentinel)
- `.claude/agent-memory/method-implementer/issue_168_phase_d_closeout.md` — Phase D campaign (Carlson seed + default flips)
- `.claude/agent-memory/method-implementer/issue_168_phase_d_step3_closeout.md` — Phase D Step 3 (Protocol + matvec wiring)
- `.claude/agent-memory/method-implementer/issue_168_phase_d_step4_closeout.md` — Phase D Step 4b (Gate 4.2 implementation)
- `.claude/agent-memory/method-implementer/issue_168_phase_c_closeout.md` — Phase C (sweep-frame matvec)

### Numerics-investigator memos

- `.claude/agent-memory/numerics-investigator/phase_d_gate_1_1_sphere_mms_diagnosis.md` — the methodology / template for the per-ordinate residual probe Step 3b will extend

### Literature memos

- `.claude/agent-memory/literature-researcher/phase_d_carlson_coupled_pole.md` — Hébert §3.9.4 math; the L=0 isotropic-only assumption in `CarlsonInwardSweep` flagged here may be relevant
- `.claude/agent-memory/literature-researcher/sphere_sn_pole_closure_canonical.md` — Pomraning 1989 structural-singularity reference

### Catalog

- `.claude/skills/vv-principles/error_catalog.md:1082+` ERR-026 entry with all phase narratives

### Test artefacts to extend

- `tests/sn/test_phase_c_crosscheck.py::test_phase_e_trajectory_resolvent_flux_shape_crosscheck` — xfailed sentinel, harness functions ready for re-enablement
- `tests/sn/test_phase_c_gates.py::test_apply_curvilinear_per_ordinate_flat_flux_residual` — Phase D Gate 1.1, the model for a parallel "Gate 1.6: per-cell scalar-flux symmetry at boundary cells" we may add in Phase F
- `tests/sn/diagnostics/gate_1_1_sphere_mms_failure.py` — template for the Step 3b per-ordinate diagnostic

### Production code to audit (Step 3)

- `orpheus/sn/operator.py` — `transport_operator_matvec_spherical` at the sweep-frame body (the pole-face IC at line ~734, the M-M angular recurrence call at the per-cell loop, the BC trace law at the boundary edge)
- `orpheus/sn/sweep.py` — `_sweep_1d_spherical` (the sweep-based path used by snapshot generation; SI not apply)
- `orpheus/sn/spatial/pole_angular_closure.py` — `_mm_weighted_angular_recurrence_single_level` (the Phase D Carlson seed consumer)
- `orpheus/sn/spatial/psi_half_angle_seed.py` — `CarlsonInwardSweep` (the Phase D seed implementation, watch out for the L=0 isotropic-only assumption flagged in the docstring WARNING block)
- `orpheus/sn/solver.py` — `_scalar_flux_from_angular` (line 655, the angle-integration helper)

### Phase E commit log

| Commit | Title |
|---|---|
| `2d3e7f2` | fix(derivations): composite per-region GL quadrature for MR Variant α solvers |
| `1070185` | test(derivations): relax interface-continuity + Garcia outer-surface ceilings post-composite-GL |
| `9417130` | test(sn): tighten Gate 4.2 k_eff rtol + add flux-shape sentinel xfail |
| `6708a4a` | chore(agent-memory): index Issue #168 Phase E closeout |

---

## Verification

### Acceptance criteria for Phase F closure

1. **Step 1 outcome documented**: snapshot freshness verified; either confirms post-Phase-D AND the 9× profile is reproducible, OR regenerates and yields a different result.
2. **Step 2 outcome documented**: SN cell-0/cell-1 ratio measured at n ∈ {40, 80, 160, 320}; either converges to 1 (O(h^p) error) or stays fixed (structural defect).
3. **If structural defect (Step 3 fires)**:
   - Production code fix in `orpheus/sn/`.
   - Phase D's Gate 1.1 (per-ordinate flat-flux residual) stays at residual ≤ 1e-12 (no regression).
   - NEW Gate 1.6 (per-cell scalar-flux at boundary cells under a non-flat eigenvector probe) PASSES.
   - `test_phase_e_trajectory_resolvent_flux_shape_crosscheck` xpasses; xfail marker removed.
   - 6 curvilinear regression snapshots regenerated (and they pass Gate 4.2 + flux-shape).
4. **If O(h^p) convergence (Step 2 alone)**:
   - Probe results documented.
   - Flux-shape sentinel updated to use a refinement-aware bound (or replaced with a `test_phase_f_boundary_cell_convergence_rate` test).
   - Follow-up issue filed for the boundary-cell convergence-rate.
   - ERR-026 #6 reclassified: "OPEN — heterogeneous eigenvector shape" → "OPEN — boundary-cell convergence rate at curvilinear boundaries" with the documented O(h^p).

### Test sequence (end-to-end gate)

```bash
# Step 2: mesh-refinement probe (~30 min)
.venv/bin/python /tmp/phase_f_sn_mesh_refinement.py

# Step 3: foundation tests stay green (~5 min)
.venv/bin/python -m pytest tests/sn/test_phase_c_gates.py tests/sn/spatial/ -q

# Step 3 + 4: flux-shape sentinel (~20 min when re-enabled)
.venv/bin/python -m pytest tests/sn/test_phase_c_crosscheck.py -v

# Phase F end-to-end gate (after fix lands; ~25 min)
.venv/bin/python -m pytest tests/sn/regression/ tests/sn/test_phase_c_*.py -v
```

### Sphinx build

After Phase F lands the production-code fix, the Sphinx narrative should describe:
- The boundary-cell defect (what it was, how it manifested)
- The fix (what changed in `orpheus/sn/`)
- The new Gate 1.6 invariant test

Dispatch the archivist for this expansion (Step 6-equivalent of Phase D's archivist dispatch pattern).

---

## Sub-agent dispatch chain (proposed)

1. **Step 1 (cheap, do it inline)**: main agent runs the npz inspection + regen probe directly.
2. **Step 2 (mesh-refinement study)**: dispatch **`numerics-investigator`** with a brief that points at this plan + the Phase E closeout. They build the probe, run it, write the Step 2 summary memo.
3. **Step 3a–3c (diagnostic deep dive)**: continue with `numerics-investigator` (resume the same instance if Step 2 just completed). They produce the per-ordinate diagnostic + identify the fix site.
4. **Step 3d–3e (production-code fix + verification)**: dispatch **`method-implementer`** with a brief that points at the numerics-investigator's fix-site identification. They implement, verify, regenerate snapshots.
5. **Sphinx narrative**: dispatch **`Archivist`** after the fix lands.
6. **Final closeout** (this is by main agent): closeout memo + sn_reshape.md update + ERR-026 manifestation table refresh + close out remaining markers if applicable.

---

## Risks + mitigations

- **Risk**: Step 2 reveals the convergence is non-monotone (similar to the pre-Phase-E Variant α MR finding). Power iteration's residual at higher refinement may produce false eigenvalues.
  - **Mitigation**: check `result.converged` and `result.iterations` for each n_total. Tighten `inner_tol` if needed; allow `max_inner` to grow as `n_total` grows.

- **Risk**: Step 3 identifies the fix site but the fix has wider side effects (e.g. changes the Phase C / D bit-identity of regenerated snapshots). The 6 curvilinear regression snapshots would need regeneration.
  - **Mitigation**: expected and acceptable. The Phase D campaign already regenerated curvilinear snapshots; Phase F may do it again. Verify the regenerated values reproduce Phase E's Gate 4.2 tightened tolerances + the new shape-comparison sentinel.

- **Risk**: The boundary-cell defect is in the SI/sweep path AND the apply matvec — fixing one but not the other re-introduces apply ↔ sweep inconsistency (which Phase C explicitly aligned).
  - **Mitigation**: any production fix must update BOTH paths AND be pinned by a test that exercises BOTH. Phase C's apply-vs-sweep consistency invariant is the gate.

- **Risk**: After Step 2 + Step 3, the SN profile converges to Variant α's profile — but Variant α may itself have a residual error (its 1.7× shape might be off, just less off than SN's 9× shape). The "right" answer requires a 3rd reference.
  - **Mitigation**: the cylinder F_N (Issue #170) or sphere F_N (Issue #171) is the structurally-distinct cross-pillar reference. If after Phase F the SN and Variant α agree but we still need ground truth, F_N becomes the next priority.

- **Risk**: The boundary-cell defect is structural in the M-M angular closure design itself (e.g. Hébert §3.9.4 has known limitations on non-flat eigenvectors that the Carlson seed alone doesn't address).
  - **Mitigation**: this would be the **literature-researcher** path — dispatch them to find references that discuss SN spherical-pole eigenvector accuracy on heterogeneous problems. Bailey-Morel-Chang 2010 NSE 165, Pomraning 1989, and Lewis-Miller 1993 are the candidate sources.

---

## Out of scope (tracked separately)

- **#195** L1 absolute magnitude on curvilinear MMS — Phase F does NOT close this manifestation (it's a different bug class).
- **#170** Cylinder F_N pillar — Phase F may benefit from F_N as a 3rd reference but doesn't directly implement it.
- **#171** Sphere F_N pillar — same.
- **Wave C-extension** LD / EC / Step cell updates — separate refactor; Phase F's boundary-cell fix should compose cleanly with future cell-update strategies via the same Protocol/composition pattern Phase D established.

---

## Final note for the post-compaction agent

When you pick this up:

1. **Re-read this plan first** before any code action.
2. **Re-load** `.claude/agent-memory/method-implementer/issue_168_phase_e_closeout.md` for the immediate prior context.
3. **Start with Step 1** even if it seems trivial — the snapshot freshness check is a cheap eliminator that may short-circuit the whole investigation.
4. **At each step's empirical-decision point**, write a short note (memo or inline) recording which branch fired and why, BEFORE proceeding. This makes the investigation auditable.
5. **Don't skip Step 2's "Variant α boundary check"** (sub-step under Step 2's bullet 3). Confirming Variant α has smooth boundaries at all refinements is what makes SN "the outlier"; if Variant α also has anomalies, the conclusion shifts.
6. **If at any point the empirical evidence contradicts the hypothesis, STOP** and re-scope. Don't improvise math fixes.
7. **Phase E test ceilings adjusted in commit `1070185` are NOT what the bug is** — those are accommodating the composite-GL extrapolation trade-off, unrelated to the SN boundary-cell defect. Don't tighten them as part of Phase F.
