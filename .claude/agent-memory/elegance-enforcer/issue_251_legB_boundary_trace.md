---
name: issue-251-legB-boundary-trace
description: Elegance review of #251 Leg B (moment-resolved boundary trace — transverse face-slope carry for 2-D Cartesian LD) on refactor/sn-foundation-cleanup — PASS-WITH-NITS
metadata:
  type: project
---

# Elegance Review: #251 Leg B — moment-resolved BOUNDARY trace (transverse face-slope carry)

Branch `refactor/sn-foundation-cleanup`, UNCOMMITTED. The boundary twin of #247 Leg A
(committed `d9396a2`). Scope: `geometry.boundary_face_layout` (storage lever),
`loss_representation._inflow_to_moments` (windowed widen) + `_octant_face_cochain` (FFW
oracle seed) + the 4 outflow capture-collapse DROPS (MFW solve/apply + FFW solve/apply),
`boundary_source_sink.prescribed_inflow` (two-arm producer), and the #251 gate block in
`test_mms_ld_2d.py`.

## Summary Verdict
**PASS-WITH-NITS** — architecturally sound; the change reads like the domain (a face
cochain carrying transverse Legendre moments) and the Leg-A NIT-1 lesson WAS internalized
(`is_moment_valued_by_flat_rank` minted in `moment_layout.py`, `is_moment_valued_by_rank`
delegates to it, Leg B reuses the primitive — zero open-coded rank arithmetic this time).
Three nits, all do-now, none block the math. Gates verified green: 11 #251 foundation,
the L1 consumption/mutation, sweep 672✓, spatial+regression 88✓ (DD/Step byte-identity
held; the 6920-ULP aniso-DD DriftWarning is pre-existing, within tol).

## Adjudications the brief asked for (all resolved)

- **Single-source reuse of rank/tail/average**: PASS for the RANK + TAIL + AVERAGE
  concepts. `is_moment_valued_by_flat_rank` (the Leg-A NIT-1(a) fix) is the single rank
  core; both `_inflow_to_moments:401` and FFW `_octant_face_cochain:1263` call it; the
  windowed/FFW discriminators are NOT a duplicated predicate — they call the SAME
  primitive (oracle twin, justified: FFW is the verification spine, never prod default).
  `face_moment_tail` single-sources the tail; `AVERAGE_MOMENT` single-sources slot 0.
  ⭐ EXCEPTION — the face-COUNT `per_axis^{d-1}` is NOT single-sourced (NIT-1 below).
- **Producer shape-discrimination** (`arr.shape == view.shape` vs `== view.shape[:-1]`):
  CLEAN (analogous to Leg A's `_build_fixed_source_rhs` shape-equality). NOT a twin of
  the consumer's rank discrimination — two DISTINCT caller contracts at two DIFFERENT
  boundaries: the producer parses CALLER INTENT (full slot vs scalar-to-broadcast) with
  the known `view.shape` in hand (shape-equality, no S4 collision hazard); the consumer
  parses TRACE SHAPE with no view (rank, S4-safe). Pattern 4 "parse at the boundary"
  done right; the fall-through `else` raise is a clean illegal-state rejection.
- **4 outflow-collapse drops**: COMPLETE — grep confirms zero remaining
  `tuple(c[..., AVERAGE_MOMENT] for c in capture)` collapses at the 4 sites; the only
  surviving `[..., AVERAGE_MOMENT]` reads are the seed-widen (`:413`), the Leg-A bulk lift
  (`:2796`, out of scope) and `s_bar` scalar-flux extraction (`:2977`, different concept).
  DD/Step byte-identity preserved BY CONSTRUCTION: `n_face_moments==1 → face_moment_tail
  == () →` no trailing axis; the OLD `if n_face_moments > 1:` branch was already skipped
  for DD, so removing it is a no-op for DD (verified: sweep 672✓ + DD regression green).
- **trace_space.py untouched = RIGHT layering, not a missed abstraction**: CONFIRMED. The
  capture→trace seam is automatic: the per-octant capture is allocated
  `np.zeros_like(boundary.face_view(face))` (`loss_representation.py:815`) /
  `np.empty_like(face)` from the moment-resolved inflow (`:1071`/`:1157`), so it inherits
  the widened `boundary_face_layout` shape; `_MovingFrontier.shed`
  (`sweep_graph.py:360-365`) writes `out_faces[axis][:,:,out_mask]` into
  `capture[axis][(slice,slice)+capture_target]` — basic slicing over axes 0,1,2 leaves the
  trailing moment axis spanned, so the cochain shed broadcasts over it for free. The
  explorer's "moment-ready by accident" prediction is correct; the trace metric/views
  broadcast over trailing axes. No TraceSpace edit needed.

## Findings

### PASS — Leg-A NIT-1 lesson INTERNALIZED: rank primitive reused, not re-spelled
`moment_layout.py` now ships `is_moment_valued_by_flat_rank(array, flat_ndim)` as the rank
CORE; `is_moment_valued_by_rank(array, reference)` delegates
(`return is_moment_valued_by_flat_rank(array, reference.ndim + 1)`). This is EXACTLY the
Leg-A NIT-1(a) preferred fix landed, and it let the SOURCE/SHAPE-only consumers (the lift,
now the boundary `_inflow_to_moments` + FFW seed) join the single source — they know only
the flat rank, not a reference array. The standing Leg-A tell ("a shape-only consumer
open-codes a 3rd spelling because the reference-array signature is friction") is RESOLVED.

### PASS — master standard: reads like the domain
A face cochain carrying transverse Legendre moments. `_inflow_to_moments` is a 3-arm
rank-discriminated widen (identity DD/Step / widen-slot-0 scalar / pass-through+validate
moment-resolved); the FFW seed mirrors it; the producer is a 2-arm shape-discriminated
slot assignment. Each arm is named with its physical contract (the bar/slope split, the
"scalar default is blind to the transverse slope — the Leg-B asymmetry"). Pattern 6 / L11
honored: production ACCEPTS the moment-resolved face, never PROJECTS it (the projector
lives at the MMS call site, `leggauss` only).

### NIT-1 (do-now, Pattern 2 / Cardinal Rule 2) — the FACE-moment count `per_axis^{d-1}` now has TWO spellings
**Location**: `geometry.py:1070` (`n_face_moments = self.scheme.spatial_basis_per_axis **
(self.ndim - 1)`) vs `loss_representation.py:323-324` (`_n_face_moments`:
`per_axis ** (self.mesh.ndim - 1)`).
**Skill reference**: Pattern 2 (single source) + anti-pattern #1 (two spellings of one
predicate) + the "unify after two instances" amendment to Pattern 6.
**Problem**: #251 introduces the SECOND spelling of the per-FACE transverse-moment count.
The cell count `per_axis^d` already has several open-coded spellings (PRE-EXISTING; see
Arch-Opp), but the FACE count was single-sourced at `_n_face_moments` until this diff
open-coded it in `geometry`. The diff's OWN docstring even cites `_n_face_moments` ("the
same … policy the cell-cochain `_n_face_moments` … key on") — the recurring Leg-A NIT-2
tell: the docstring names the primitive the code does not call.
**Bug-habitat argument**: the storage PRODUCER (`geometry.boundary_face_layout`, sets the
trace slot width) and the cochain CONSUMER (`loss_representation`, allocates/sheds into
that slot) must agree on the face width or the capture-shed seam silently mis-shapes. They
are coextensive TODAY by the identical formula; the day the face-moment policy changes (a
non-tensor face basis, a scheme-dependent face count that diverges from `per_axis^{d-1}`,
or a d=3 face that is `per_axis^2`), the change lands in one module and the other silently
disagrees — a wrong-shaped trace the rank discriminator may even accept.
**Required change**: single-source the face count. Mechanically feasible homes (need both
`per_axis` from the scheme AND `ndim` from the mesh):
  (a) **Preferred** — a free function `face_moment_count(per_axis, ndim) -> int` in
      `moment_layout.py` (the home of all moment-layout policy, beside `face_moment_tail`),
      called by BOTH `geometry.boundary_face_layout` and `_n_face_moments`. Mirrors the
      `is_moment_valued_by_flat_rank` single-source move.
  (b) Or a mesh property `SNMesh.n_face_moments` that `_LossRepresentation._n_face_moments`
      also delegates to (precedent: `scheme.is_multi_moment` exists "rather than
      re-spelling `spatial_basis_per_axis > 1`", scheme.py:868). geometry then reads
      `self.n_face_moments`.
Either collapses the two spellings; (a) keeps the policy in the layout-policy leaf.

### CONCERN (oracle-only, do-now-if-cheap) — FFW seed lacks the width-validation the windowed twin has
**Location**: FFW `_octant_face_cochain:1263-1266` (`elif
is_moment_valued_by_flat_rank(inflow[a], flat_face_ndim): buf[tuple(in_edge)] = inflow[a]`)
vs the windowed `_inflow_to_moments:401-409` (which raises `ValueError` on
`face.shape[-1] != n`).
**Skill reference**: Pattern 4 (illegal states) + symmetry-in-math ⇒ symmetry-in-code.
**Problem**: the windowed widen validates the trailing moment width and raises a clear
"expected 2^(d-1)" ValueError; the FFW oracle seed discriminates with the SAME primitive
but does NOT validate — it assigns straight into the `(…, n)`-tailed buffer view. Probed
numpy behavior: a width-`n` (correct) or width-3 (wrong) input behaves as expected
(correct OK; 3 raises a GENERIC broadcast ValueError, no "2^{d-1}" guidance), BUT a
**width-1 trailing axis silently BROADCASTS** the single moment across all `n` slots — a
wrong-physics seed that does not raise. A moment-resolved `(N,ng,*transverse,1)` face is
rank-classified as moment-valued and then silently fanned out.
**Bug-habitat argument**: small and oracle-only (FFW is the verification spine, never the
production default; the windowed prod path IS guarded), so this is a CONCERN not a
VIOLATION — but the FFW oracle's job is to be the bit-exact reference, and a silently
mis-broadcast seed would make the oracle agree with a wrong window for the wrong reason.
The asymmetry (one twin validates, the other does not) is itself the smell: a future
maintainer reading the FFW path will assume the same guard the windowed path advertises.
**Required change**: factor the moment-resolved seed validation into the shared primitive
both paths call. Cleanest: have `_inflow_to_moments` BE the single seed-normalizer and
have the FFW `_octant_face_cochain` call it (the FFW currently re-implements the
discrimination inline) — collapsing the oracle-twin discrimination into ONE method that
validates once. If that is judged too large for Leg B, at minimum lift the
`face.shape[-1] != n` check into a tiny shared helper
(`assert_face_moment_width(face, n)`) called from both the windowed widen and the FFW seed.

### NIT-2 (do-now, anti-pattern #11) — geometry docstring names `_n_face_moments`/`_spatial_moment_tail` it does not call
**Location**: `geometry.py` docstring (the "Spatial-moment tail" block, ~:1028-1044) +
the inline comment `:1068-1069`.
**Skill reference**: anti-pattern #11 (a claim the code does not honor is a bug habitat).
**Problem**: the docstring asserts the lever uses "the same … policy the cell-cochain
`_n_face_moments` / `_spatial_moment_tail` key on" and the same `face_moment_tail` single
source — true for the TAIL (it does call `face_moment_tail`), FALSE for the COUNT (it
open-codes `per_axis^{d-1}`, NIT-1). A future maintainer greps for the shared count
primitive, finds none, and is misled about the single source.
**Required change**: auto-resolved by applying NIT-1(a)/(b) (the code then matches the
docstring). If NIT-1 is deferred, soften the docstring to state the count is open-coded
here, matching `_n_face_moments` by hand, with a removal trigger.

## Architectural Opportunities
The CELL-moment count `per_axis^d` has MULTIPLE pre-existing open-coded spellings
(`_bases.py:196`, `spatial_moment_space.py:206/256`, `loss_representation.py:802`,
`_spatial_moment_tail:358`) — the cell analog of NIT-1, but PRE-EXISTING (not introduced
by #251). If NIT-1(a) lands a `face_moment_count` helper, a sibling `cell_moment_count`
(or one `moment_count(per_axis, n_axes)` taking the axis count) would let all the cell-side
spellings join too. Worth a `module:sn` issue (`type:improvement`) so a fresh session can
sweep both counts onto the moment-layout leaf — do NOT fold the cell sweep into Leg B
(scope), but the face count IS in Leg B's scope (NIT-1).

## Approval Conditions
1. **NIT-1** — single-source the face-moment count `per_axis^{d-1}` (preferred:
   `face_moment_count(per_axis, ndim)` in `moment_layout.py`, called by BOTH
   `geometry.boundary_face_layout` and `_n_face_moments`; or a `SNMesh.n_face_moments`
   property both delegate to). The two spellings must collapse to one.
2. **CONCERN** — give the FFW oracle seed the same moment-width validation the windowed
   widen has (preferred: route the FFW seed through `_inflow_to_moments` so the
   discrimination + validation live in ONE method; else a shared
   `assert_face_moment_width` helper). Closes the silent width-1 broadcast in the oracle.
3. **NIT-2** — auto-resolved by NIT-1; else soften the geometry docstring's single-source
   claim for the count.

The rank/tail/average single-source (PASS), the no-projection-in-production (PASS), the
4 collapse drops (COMPLETE, DD byte-identical), the trace-untouched layering (CORRECT),
and the producer two-arm contract (CLEAN) are all sound. With NIT-1 + the CONCERN applied
the face-moment COUNT joins the single source and the oracle twin is guarded — then the
carve is clean to commit.
