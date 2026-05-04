# Plan — Issue #138 full collapse (cyl/sph façades retired; slab native = V-of-V)

**Date:** 2026-04-29
**Branch (current):** `refactor/peierls-facade-narrow` (carries the prior narrow commit `a2b6bcd`)
**Branch (new):** `refactor/peierls-138-full-collapse` (cut from `a2b6bcd`)
**Issue:** #138
**Predecessor commit:** `a2b6bcd` (the narrow-scope refactor)

---

## Decision context

User invoked Cardinal Rule 2 (Infrastructure is critical) to override the
`peierls_nystrom.rst §theory-peierls-api-posture` "permanent wrappers"
designation. Carve-out: **keep `peierls_slab.py` (the native E₁ Nyström)
as verification-of-verification** until a future session implements the
parallel Green's function approach and uses it as the new oracle. At
that point the slab native path will be retired too.

This means the original Issue #138 Phase B (lift `_basis_kernel_weights`
into the unified slab branch and deprecate `solve_peierls_eigenvalue`)
is **explicitly out of scope** — fusing the kernel-weight machinery
into the unified path would destroy the cross-check (both routes would
share the same numerics).

Phase A (cyl/sph façade collapse) is **explicitly in scope** and goes
all the way: wrappers, dataclasses, adapters all retired. Façades
become registry-only modules.

---

## End-state

### `orpheus/derivations/peierls_cylinder.py` — registry-only (~150 LoC)

Keep:
- Module docstring (rewritten)
- `from . import peierls_geometry as _pg`
- `GEOMETRY = _pg.CYLINDER_1D` singleton
- `_MAT_IDS_CYL = {1: [2]}` constant
- `_F4_CYL_TOL = {0.1: "1.4 %", 0.2: "5.4 %", 0.3: "13 %"}` constant
- `_build_peierls_cylinder_case(...)` (calls `_pg.solve_peierls_1g(GEOMETRY, ...)` directly, returns `_pg.PeierlsSolution`)
- `_build_peierls_cylinder_hollow_f4_case(...)` (calls `_pg.solve_peierls_mg(geometry, ...)` directly, returns `_pg.PeierlsSolution`)

Drop:
- `composite_gl_r`, `_lagrange_basis_on_panels` re-exports — already dropped in narrow commit
- `_rho_max`, `_optical_depth_along_ray`, `_which_annulus`, `build_volume_kernel`, `compute_P_esc`, `compute_G_bc`, `build_white_bc_correction` — already dropped in narrow commit
- `continuous_cases() -> []` stub — already dropped in narrow commit
- `PeierlsCylinderSolution` dataclass + `phi()` method
- `_soln_to_cylinder` adapter
- `solve_peierls_cylinder_1g` wrapper
- `solve_peierls_cylinder_mg` wrapper

### `orpheus/derivations/peierls_sphere.py` — registry-only (~150 LoC)

Mirror of cylinder. Keep `_MAT_IDS_SPH`, `_F4_SPH_TOL`, `_build_peierls_sphere_case`, `_build_peierls_sphere_hollow_f4_case`. Drop `PeierlsSphereSolution`, `_soln_to_sphere`, `solve_peierls_sphere_{1g,mg}`.

### `orpheus/derivations/peierls_slab.py` — UNCHANGED

V-of-V status preserved. No code touched. Sphinx docstring may add
a forward-reference to the Green's function plan as the retirement gate.

### `orpheus/derivations/peierls_geometry.py` — UNCHANGED

Already canonical. No additions needed.

### `orpheus/derivations/peierls_cases.py` — UNCHANGED logic, lazy-import targets stable

The lazy imports at lines 128, 139, 157 still hit
`_build_peierls_*_case` / `_build_peierls_*_hollow_f4_case` functions
that survive in the registry-only modules. The `_F4_CYL_TOL` / `_F4_SPH_TOL`
imports at lines 426–427 also survive.

### Sphinx docs

`peierls_nystrom.rst`:
- §theory-peierls-api-posture: rewrite. Move `solve_peierls_cylinder_*`,
  `solve_peierls_sphere_*` rows from "Permanent" to a new "Retired
  (2026-04-29, Issue #138)" subsection. Keep `solve_peierls_mg`,
  `solve_peierls_1g` as the canonical permanent surface. Document the
  L104a re-application: "Wrappers are temporary unless really
  necessary; the cyl/sph wrappers' parameter renames and shape-specific
  dataclasses turned out not to clear that bar after Cardinal Rule 2
  scrutiny".
- §theory-peierls-slab-polar-retirement: append a "Retirement gate"
  subsubsection citing the Green's function plan as the trigger for
  full slab retirement. Reference
  `.claude/plans/peierls-greens-function-approach.md` (existing
  scratch plan).
- Lines 950–998 (Permanent table), 988–998 (V-of-V table): update to
  reflect new posture.
- Migration notes: add a `:rubric:` listing canonical replacements
  (`solve_peierls_mg(_pg.CYLINDER_1D, ...)`, `PeierlsSolution`).

`collision_probability.rst`:
- Lines 3081, 3223, 3227, 4060, 4067, 4280, 4283, 4285 — references to
  `solve_peierls_cylinder_1g`, `PeierlsCylinderSolution`,
  `solve_peierls_sphere_1g`, `PeierlsSphereSolution`. Redirect each to
  the canonical `_pg` symbol (with `CYLINDER_1D` / `SPHERE_1D` binding).

---

## LoC delta projection

| File | Before narrow | After narrow | After this | Δ from narrow |
|---|---|---|---|---|
| `peierls_cylinder.py` | 716 | 558 | ~150 | −408 |
| `peierls_sphere.py`   | 694 | 556 | ~150 | −406 |
| `peierls_slab.py`     | 718 | 705 | 705  |  0 |

Plus one Sphinx posture rewrite, one cross-ref sweep, ~14 test sites
rewritten.

---

## Phase order — sequenced for safety

The order minimises broken-baseline windows. Each phase ends with a
green test run before proceeding.

### Phase 0 — Cut the branch and lock baseline

1. `git checkout -b refactor/peierls-138-full-collapse` (off `a2b6bcd`)
2. Run the same 200-test green-baseline subset that passed under the
   narrow commit:
   - batch 1 (migrated cyl/sph + flux): 80 tests
   - batch 2 (foundational + rank-n): 116 tests
   - test_peierls_convergence (L0 slab): 4 tests
3. Confirm Sphinx -W still builds clean.
4. Snapshot the L0 error-catalog ERR coverage (30/31).

### Phase 1 — Sphinx posture flip (docs lead the code)

This phase is **docs-only** so any Sphinx warning surfaces immediately.

1. Rewrite §theory-peierls-api-posture to flip the posture:
   - "Permanent" table: only `solve_peierls_mg`, `solve_peierls_1g`.
   - New "Retired (2026-04-29, Issue #138)" table: cyl/sph wrappers
     + dataclasses, with `_pg.solve_peierls_*` + `PeierlsSolution` as
     the canonical replacement and a brief migration recipe.
   - "V-of-V" table: unchanged for now (slab native still listed).
   - Add a `Why retired` paragraph citing Cardinal Rule 2 + L104a
     re-application + the user's 2026-04-29 decision.
2. Add a `Retirement gate` subsubsection under
   §theory-peierls-slab-polar-retirement specifying that the slab
   native path will be retired once the Green's function approach
   ships and demonstrates parity from a third independent route.
3. Run `sphinx-build -W` — expected to produce **broken cross-ref
   warnings** for the not-yet-deleted `peierls_cylinder.solve_*` etc.
   Use `nitpick_ignore` or temporary `:noindex:` only as a *last*
   resort — preferred fix is to update each cross-ref site to the
   canonical symbol now (the deletion in Phase 4 then introduces no
   new warnings).
4. Commit: `docs(peierls): #138 — flip api-posture; cyl/sph wrappers retired`.

### Phase 2 — Migrate test files off cyl/sph wrappers + dataclasses

For each test that imports a soon-to-be-deleted symbol, replace with
the canonical `_pg` call. The migration is mechanical:

| Old call | New call |
|---|---|
| `solve_peierls_cylinder_1g(radii, sig_t, sig_s, nu_sig_f, ..., n_beta=N, ..., n_phi=M, ...)` | `_pg.solve_peierls_1g(GEOMETRY, radii, sig_t, sig_s, nu_sig_f, ..., n_angular=N, ..., n_surf_quad=M, ...)` |
| `solve_peierls_sphere_1g(radii, sig_t, sig_s, nu_sig_f, ..., n_theta=N, ..., n_phi=M, ...)` | same pattern |
| `solve_peierls_cylinder_mg(...)` | `_pg.solve_peierls_mg(GEOMETRY, ...)` |
| `solve_peierls_sphere_mg(...)` | `_pg.solve_peierls_mg(GEOMETRY, ...)` |
| `PeierlsCylinderSolution(r_nodes=..., n_quad_y=...)` | `PeierlsSolution(r_nodes=..., geometry_kind="cylinder-1d", n_quad_angular=..., ...)` |
| `PeierlsSphereSolution(r_nodes=..., n_quad_theta=...)` | `PeierlsSolution(r_nodes=..., geometry_kind="sphere-1d", n_quad_angular=..., ...)` |

**Inner-radius variants:** when the old wrapper had `inner_radius != 0`,
it constructed a fresh `_pg.CurvilinearGeometry(kind="cylinder-1d", inner_radius=...)`.
The migrated callers must construct the same geometry inline and pass
it to `_pg.solve_peierls_*`.

Files to edit (concrete site count):

- `tests/derivations/test_peierls_cylinder_eigenvalue.py` —
  `solve_peierls_cylinder_1g` ×8 (lines 70, 90, 96, 115, 132, 153, 177, 197, 221)
- `tests/derivations/test_peierls_cylinder_white_bc.py` —
  `solve_peierls_cylinder_1g` ×4 (lines 114, 131, 138, 155)
- `tests/derivations/test_peierls_sphere_eigenvalue.py` —
  `solve_peierls_sphere_1g` ×4 (lines 50, 68, 99, 129)
- `tests/derivations/test_peierls_sphere_white_bc.py` —
  `solve_peierls_sphere_1g` ×1 (line 56)
- `tests/cp/test_peierls_cylinder_flux.py` —
  `PeierlsCylinderSolution(...)` constructor ×1 (line 94) + import
- `tests/cp/test_peierls_sphere_flux.py` —
  `PeierlsSphereSolution(...)` constructor ×1 (line 96) + import
- `tests/derivations/test_peierls_multigroup.py` lines 499 and 520 —
  `_build_peierls_cylinder_hollow_f4_case` and
  `_build_peierls_sphere_hollow_f4_case` are KEPT, no migration needed

After every file edit, run that file's tests in isolation to confirm
parity. (Same chunk-by-chunk discipline used in narrow commit to dodge
the OOM ceiling.) **The migrated tests use the CURRENT wrapper code as
their oracle** — they should pass before Phase 3 starts.

Commit cadence: one commit per geometry pair, e.g.
- `test(peierls): #138 — migrate cylinder eigenvalue + white_bc tests off façade wrappers`
- `test(peierls): #138 — migrate sphere eigenvalue + white_bc tests off façade wrappers`
- `test(peierls): #138 — migrate flux dataclass constructors to canonical PeierlsSolution`

### Phase 3 — Migrate registry constructors off cyl/sph wrappers

Inside `peierls_cylinder.py` and `peierls_sphere.py`, the
`_build_peierls_*_case` and `_build_peierls_*_hollow_f4_case`
functions currently call the soon-to-be-deleted
`solve_peierls_cylinder_*g` / `solve_peierls_sphere_*g` wrappers.
Rewrite each to call `_pg.solve_peierls_1g(GEOMETRY, ...)` /
`_pg.solve_peierls_mg(geometry, ...)` directly.

Specifically, `_build_peierls_cylinder_case` constructs a
`PeierlsCylinderSolution` for the normalised flux. Replace that
construction with `dataclasses.replace(sol, phi_values=phi_normed[:, np.newaxis])`
on the canonical `PeierlsSolution`. Same for sphere and the hollow_f4
variants.

The `phi_fn` closure that wraps `sol.phi(...)` works unchanged — the
canonical `PeierlsSolution.phi()` method is identical to the
shape-specific ones (verified in `peierls_geometry.py:5480-5493`).

Commit: `refactor(peierls): #138 — registry constructors call _pg.solve_peierls_* directly`.

### Phase 4 — Drop the wrapper code from cyl/sph façades

With all consumers migrated, delete:
- `solve_peierls_cylinder_1g`, `solve_peierls_cylinder_mg` (lines 263–396 of cylinder.py post-narrow)
- `_soln_to_cylinder` adapter
- `PeierlsCylinderSolution` dataclass
- (mirror for sphere)

Update each module's docstring to reflect the registry-only role.
Verify `_pg`, `numpy`, `dataclasses` (no longer needed if no `@dataclass` remains in the file — verify), `_reference`, `_xs_library` imports are still all used; remove anything orphaned.

Run the test batches again. Run sphinx -W — should be clean (Phase 1
already redirected the cross-refs).

Commit: `refactor(peierls): #138 — delete cyl/sph wrapper code; façades are registry-only`.

### Phase 5 — Verification + audit

1. Run all peierls test files chunk-by-chunk (200-test green
   baseline as Phase 0). All should pass.
2. Run the L0 convergence gate for slab (`test_peierls_convergence.py`).
   Slab path is unchanged, but verify nothing collateral-broke.
3. Run `python -m tests._harness.audit` — expect ERR coverage 30/31
   identical to baseline; no new orphan equations (`peierls-*` labels
   should still resolve to surviving symbols).
4. Run `sphinx-build -W` — expect clean.
5. (Optional) Spot-check Nexus: `mcp__nexus__provenance_chain` from
   `peierls-cylinder-equation` / `peierls-sphere-equation` to a code
   node — should land on `_pg.solve_peierls_mg` (not a façade wrapper).

### Phase 6 — Merge, close, log

1. Fast-forward merge: `git checkout main && git merge --ff-only refactor/peierls-138-full-collapse`.
2. Push: `git push origin main`.
3. `gh issue close 138` with a comment summarising scope (cyl/sph fully
   collapsed; slab kept as V-of-V pending Green's function in a future
   session) + LoC delta + commit list.
4. Log to `.claude/lessons.md` if any new lesson emerged from the
   migration (e.g., a non-obvious test discovery).
5. Delete branch local + remote.

---

## Risks & mitigations

### R1 — Breaking the L0 convergence gate (Cardinal Rule 5)

Slab native path is untouched in this plan. L0 risk only arises if a
test that *does* call slab somehow imports a deleted symbol via a
chain. **Mitigation:** Phase 0 baseline establishes the green-state
fingerprint; Phase 5 re-verifies. Phase 2 commits are per-geometry so
bisection is trivial if a test surprises us.

### R2 — Sphinx cross-ref churn

44 cross-refs to façade public APIs across `peierls_nystrom.rst` and
`collision_probability.rst` (per Issue #138 audit). **Mitigation:**
Phase 1 lands the Sphinx churn first, so subsequent code deletions
introduce no new warnings.

### R3 — Test files I missed during the narrow migration

Sanity check before Phase 4: `grep -rn "from orpheus.derivations.peierls_\(cylinder\|sphere\) import" tests/`
should show only `_build_peierls_*` and `GEOMETRY` and (during the
phase) the still-being-migrated `solve_peierls_*g` lines that Phase 2
has not yet touched.

### R4 — `dataclasses.replace` on canonical PeierlsSolution

The canonical dataclass is `frozen=True`. `dataclasses.replace(sol, ...)`
is the documented pattern for frozen dataclass updates. **Verification:**
write the replacement once in `_build_peierls_cylinder_case` (Phase 3),
run the test that exercises the normalised-flux path (the
`test_peierls_cylinder_eigenvalue.py` cases that go through the
registry path), confirm bit-equivalence vs the baseline.

### R5 — `inner_radius` geometry construction churn

The wrappers had a `geometry = _pg.CurvilinearGeometry(kind="...", inner_radius=...)
if inner_radius != 0.0 else GEOMETRY` pattern. Tests calling the
hollow path will need this geometry construction inline. To DRY,
consider a tiny module-level helper `_geometry_for(inner_radius)` in
each migrated test file. **Decision deferred to Phase 2** — only
introduce the helper if more than one test site needs it.

### R6 — OOM ceiling on large pytest runs

The container's ~4GB RAM ceiling kills suites > ~150 mpmath-heavy
tests in one process. **Mitigation:** Phase 0 and Phase 5 use the
same chunk-by-chunk batching that worked for the narrow commit
(batch 1 cyl/sph, batch 2 foundational, test_peierls_convergence
single-file).

### R7 — `peierls_cases.py` lazy-import targets

`peierls_cases.py` lazy-imports `_build_peierls_cylinder_hollow_f4_case`,
`_build_peierls_sphere_hollow_f4_case`, `_F4_CYL_TOL`, `_F4_SPH_TOL`.
All survive in the registry-only modules. **Verification:** Phase 5
runs the registry path end-to-end via
`reference_values.continuous_get("peierls_cyl1D_hollow_1eg_1rg_r0_10")`
or equivalent.

---

## Estimated session count

- Phase 0–1: 1 session (~30 min, mostly Sphinx writing)
- Phase 2: 1 session (mechanical, but careful — 6 test files × 1–8 sites each)
- Phase 3–4: 1 session (the actual code deletion)
- Phase 5–6: 1 session (verification + merge)

Total: 3–4 sessions for the full collapse. The narrow commit took ~1 session.

## What this plan does NOT do

- Lift `_basis_kernel_weights` into the unified slab branch (Issue #138
  Phase B.1 original) — would destroy V-of-V cross-check
- Delete `peierls_slab._build_kernel_matrix` / `_build_system_matrices`
  (Issue #138 Phase B.2 original) — those are the V-of-V machinery
- Deprecate `solve_peierls_eigenvalue` (Issue #138 Phase B.3 original)
  — that IS the V-of-V entry point
- Update `tests/l0_error_catalog.md` ERR-NNN entries — slab native
  path unchanged

These are all deferred to the post-Green's-function-implementation
session, when the third independent route is available to assume the
oracle role.
