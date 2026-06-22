---
name: issue-240-d5b-s3-a-inc-c-closeout
description: #240/#38/#37 D5b-S3-A — fold Increment C (spatial-moment iterate φ̂ + Σ_s·φ̂ scattering + source seams) into the SN SI/forward path. IN-FLIGHT (NOT complete). LANDED + VERIFIED the load-bearing Σ_s⊗I_spatial scattering-lift bridge (byte-identical at per_axis==1; GATE 4 holds 513/1/4). BLOCKER surfaced + characterized: carrying φ̂ in the iterate requires WIDENING the typed-field space-shape contract (ScalarField/AngularField/MomentField validate shape==(ng,*spatial) — no room for a trailing 2^d axis) — a foundational typed-field change that is the make-or-break design decision and was NOT in the crosswalk's seam map. The scattering lift + the field-space widening are the two halves; lift done, widening owed.
metadata:
  type: project
---

# #240 D5b-S3-A — Inc C moment-iterate fold: IN-FLIGHT closeout

**Branch** `feature/sn-space-angle-tier2`, base HEAD `f45a219` (S2 done).
**NOT committed** (main agent reviews + commits). Host env, `.venv/bin/python -O`.
Canonical `python -O -m pytest`; NEVER all `tests/sn` (#212).

## STATUS — what is DONE vs OWED

**DONE + VERIFIED (safely shippable now):**
- The `Σ_s ⊗ I_{spatial-moment}` scattering lift (the LOAD-BEARING bridge,
  crosswalk §"The Bridge rows" #1). Three einsum subscripts in
  `orpheus/sn/material_xs_field.py` gained a trailing `...` spectator on the
  cell axis:
  - `apply_p0_in_scatter` (`:468`): `"fg,fc->gc"` → `"fg,fc...->gc..."`.
  - `apply_n2n` (`:494`): `"fg,fc->gc"` → `"fg,fc...->gc..."`.
  - `apply_legendre_scattering_moments` (`:517`): `"mfc,fg->mgc"` →
    `"mfc...,fg->mgc..."`.
  This makes `Σ_s` carry NO spatial-moment index → a trailing `2^d` axis rides
  through every per-material group-matmul as a SPECTATOR broadcast. **The
  projection primitives** (`MomentProjection.apply`/`apply_transpose`,
  `HarmonicMomentReconstruction.apply` in `orpheus/numerics/projection.py`)
  ALREADY use `...` for trailing axes — verified, no change needed. The
  scattering `apply` arms + `_aniso_source_from_moment_values` +
  `_assemble_per_ordinate_source` thread the trailing axis through unchanged
  (the `/sum_w` is a scalar divide; field `+` is value-wise).

**THE `apply_p0_in_scatter` BROADCAST VERDICT (test-architect's flagged
assumption) — RESOLVED, with a CORRECTION to the crosswalk:** the crosswalk
§"The Bridge rows" #1 (and the brief) ASSUMED `apply_p0_in_scatter` "broadcasts
over the trailing axis" and that no reshape was needed. **That assumption is
FALSE as the code stood:** the einsum `"fg,fc->gc"` hard-codes the cell axis as a
SINGLE index `c`, so a rank-3 `phi (ng, n_cells, 2^d)` RAISES
(`operand has more dimensions than subscripts given ... no '...'`). The fix is
NOT a reshape — it is the one-character subscript change `fc->gc` ⇒
`fc...->gc...` (add the `...` spectator). **Verified rank-2-exact**
(`np.array_equal(einsum("fg,fc->gc",...), einsum("fg,fc...->gc...",...)) == True`
when no trailing axis) → BYTE-IDENTICAL at DD/Step (`per_axis==1`, the trailing
axis is absent → `...` matches nothing). This is the cleanest possible Pattern-7
producer-side lift. **The brief's "if it does not broadcast you must re-check DD
bit-id" branch is satisfied: the `...` change IS the fix AND DD bit-id is proven
(GATE 4 = 513/1/4 unchanged).**

**GATE 4 (DD/Step bit-identity, the negative control) — HOLDS:**
`python -O -m pytest tests/sn/sweep/core tests/sn/solve -W "error::tests.sn.regression._regression_assert.DriftWarning"`
= **513 passed / 1 skipped / 4 xfailed**, IDENTICAL pre/post the scattering lift.
(The test-architect memo said "562" but the live count at HEAD `f45a219` is
513/1/4 — the memo flagged "RE-CONFIRM at pickup"; 513 is the live S2 baseline,
matching the prior strict-gate count.) `tests/sn/operators tests/sn/spatial`:
the scattering lift introduced ZERO new failures (561 pass / 7 fail, the 7 being
the documented PRE-EXISTING reds — sphere 1-D matvec SPH ×3, `Face 'ymin' mu_y`
×2, sphere curvilinear apply ×2; confirmed via `git stash` that they fail at
clean HEAD too).

**OWED (the rest of S3-A — NOT done; see §BLOCKER for why):**
1. The φ̂ iterate carrier widening (BOTH carriers — `AngularFlux` un-windowed +
   `HarmonicMomentField` windowed). **BLOCKED on the typed-field space widening.**
2. The cell-emit φ̂ accumulation (`_CellSolve.cell` sweep_graph.py:883-895; the
   moment-mode emit must accumulate the full `2^d` spatial axis as `...lmgdp`).
3. The two source seams (d≥2 wavefront `_ubld_system` genuine `(2^d,ng)` Q;
   d=1 scan `D1ClosedForm.kernel_rhs`/`schur_xV` slope thread + the affine `b`).
4. The 1-D matvec slope iterate (the #37 krylov-slab tripwire SUT).
5. GATE 1 (flip the 1-D tripwire + the 2-D analog), GATE 2 re-run, GATE 5.

## ⚠ THE BLOCKER (the make-or-break design decision — NOT in the crosswalk)

**Carrying φ̂ in the iterate requires WIDENING the typed-field space-shape
contract.** The crosswalk's seam map (file:line) is accurate for the SWEEP /
SCATTERING / SOURCE-CARRIER plumbing, but it does NOT address the TYPED-FIELD
layer, which is where the spatial-moment axis must LIVE between sweeps:

- The iterate is `TimedFullField(bulk=AngularFlux | HarmonicMomentField)`.
- `Field.__post_init__` (`orpheus/numerics/field.py:197`) GATES
  `self.values.shape != self.space.shape → raise`. `AngularField` /
  `ScalarField` / `MomentField` / `AngularSourceSink` / `ScalarSourceSink`
  spaces are `(ng, *spatial_shape)` (or `(N, ng, *spatial)` / `(L+1, 2L+1, ng,
  *spatial)`), derived in `from_mesh` / `zeros_on` from `mesh.spatial_shape`.
  **There is NO slot for a trailing `2^d` spatial-moment axis.** A bulk field
  carrying `(N, ng, *spatial, 2^d)` violates the space-shape gate at
  construction (Pattern 4 — illegal-states-unrepresentable — firing CORRECTLY;
  the trailing axis is currently an illegal state).
- `ScalarSourceSink.zeros_on(mesh)` allocates `(ng, *spatial)` — the scattering
  output accumulator `iso` has NO trailing axis, so the moment-aware
  `add_iso_source(iso, phi_with_2^d)` would broadcast-FAIL on the in-place `+=`
  even though `apply_p0_in_scatter`'s einsum now accepts it (the ACCUMULATOR is
  the wrong shape, not the einsum).

**What is NOT a blocker:** `TimedFullField.to_flat`/`from_flat`
(`timed_full_field.py:439`) use `template.bulk.values.shape` + `.ravel()` →
shape-agnostic, the Krylov flat vector grows automatically with the bulk shape.
`AngularFlux.integrate_angular` (`einsum("n,ng...->g...")`) is moment-agnostic.

**The design choice (owed — a focused design pass FIRST):** the typed-field
spaces need an OPTIONAL trailing spatial-moment axis (length `per_axis^d`),
gated so `per_axis==1` keeps the EXACT current shape (no trailing axis → DD
byte-identical). Candidate: a `spatial_moments: int = 1` parameter on the
field-space factories (`from_mesh` / `zeros_on` / `from_mesh_and_L`), threaded
from `sn_mesh.scheme.spatial_basis_per_axis ** ndim`, that appends `(n,)` to
`space.shape` iff `n > 1` (reuse `_ubld.face_moment_tail`'s policy = the single
source for "append iff > 1"). This is a Pattern-4 widening of the make-illegal-
states-unrepresentable gate: the spatial-moment axis becomes a FIRST-CLASS,
typed part of the field space, not a bare-ndarray trailing axis. It ripples to
every `from_mesh`/`zeros_on` consumer but is byte-identical at `per_axis==1`.

⭐ **RECOMMEND** before resuming: dispatch a design pass (or proceed with the
`spatial_moments` factory parameter) to widen the field spaces FIRST, then the
iterate/emit/source seams become straightforward (the einsum lift is already
moment-agnostic). The scattering lift landed here is the half that needed NO
field-space change (it operates on raw `.values` ndarrays inside the per-
material loop, below the space gate); the field-space widening is the other half.

## Files changed (this session)
- `orpheus/sn/material_xs_field.py` — 3 einsum subscripts gained `...`
  (the `Σ_s ⊗ I_spatial` lift). Docstrings updated (the lift + the rank-2-exact
  no-op invariant). NOTHING ELSE touched.

## Gates run (this session)
- GATE 4: `tests/sn/sweep/core tests/sn/solve -W error::DriftWarning` →
  **513P/1skip/4xf**, IDENTICAL pre/post (byte-identical by construction).
- `tests/sn/operators tests/sn/spatial` → 561P/7F (7 PRE-EXISTING, stash-confirmed).
- LD MMS floor (pre-change baseline, for GATE 2 reference):
  `tests/sn/verification/mms/test_mms_ld_slab.py + test_mms_ld_2d.py` →
  **8P/1xf** (the 1xf is `test_ld_thick_diffusive_limit_xfail`, the #37 tripwire
  to flip; all gate-2 transport-MMS green).

## DID NOT gate on S3-FP==S2-FP
Confirmed: no gate asserts the S3 fixed point equals the S2 flat-source FP. Inc C
is a PHYSICS-COMPLETION (the FP CHANGES — diffusion-limit-consistent). The
landed scattering lift is byte-identical ONLY at `per_axis==1` (DD/Step), where
`S_full ≡ S_flat` by construction (the `...` matches nothing); at LD the lift
will move the FP (correctly) once the iterate carries φ̂.

## Out-of-scope boundary hit
- The d≥2 Krylov 2-D LD matvec raise (`_CellResidual.cell` sweep_graph.py:929)
  is S3-B — UNTOUCHED.
- The typed-field space widening (the BLOCKER above) was NOT in scope as written
  but is a HARD PREREQUISITE for the iterate carrier — surfaced + characterized
  here for the resume.

## LESSON (candidate for the crosswalk / coding-elegance)
A Pattern-7 producer-side lift that lives in a per-material loop on RAW `.values`
ndarrays (below the typed-field space gate) is cheap + byte-identical (the
`fc->gc` ⇒ `fc...->gc...` change). But the iterate that FEEDS it lives ABOVE the
space gate (typed `Field` with a rigid `(ng,*spatial)` space). A crosswalk that
maps the sweep/scattering/source plumbing but skips the TYPED-FIELD space layer
will under-scope the carve: the einsum accepts the wider axis, but the typed
accumulator/iterate cannot HOLD it without a space-contract widening. **Map the
field-space layer (the Pattern-4 shape gate) as an explicit crosswalk row
whenever a carve adds an axis to a quantity that travels between solver steps.**
