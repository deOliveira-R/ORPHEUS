---
name: issue-246-moment-axis-predicate
description: Design-scoping map for #246 (4 value-based moment-axis shape probes → one typed SpatialMomentSpace predicate); in-scope data at each site, plumbing verdict, predicate home
metadata:
  type: project
---

# #246 — typed `has_spatial_moment_axis` predicate scoping (branch `refactor/sn-foundation-cleanup`)

**Why:** #240 D5b-S3 (ERR-061) left a 4th value-based "does this array carry the
`2^d` spatial-moment axis?" shape probe (`_reframe`). The elegance-enforcer flagged
it as the 4th of a family; #246 wants ONE typed predicate keyed on the
`SpatialMomentSpace` factor (#64 / homed `orpheus/numerics/spaces/spatial_moment_space.py`).

**How to apply:** when implementing #246, this is the in-scope map. The CRUX finding
governs whether it is a clean predicate-swap (it is NOT, fully — see below).

## ⭐ THE CRUX FINDING (governs the whole issue)

The 4 probes split into **two kinds**, and only ONE kind is a clean typed-predicate swap:

- **Kind A — pure predicate** ("is the moment axis present?", boolean): the `_reframe`
  guard (Site 1) and the `if n_face_moments > 1:` capture-reduction (4 sites inside Site 2's
  consumers) and `_moment_broadcast_sigma`'s rank test (Site 3). These ask a yes/no question.
- **Kind B — width/count** (the actual integer `2^{d-1}` or `2^d`, needed for buffer
  ALLOCATION and to thread `n_face_moments=` into `walk_windowed`/`walk_full`/
  `_octant_face_cochain`/`_inflow_to_moments`): `_n_face_moments` (the property body),
  `_spatial_moment_tail`. These CANNOT become a boolean — they need the number.

#246 should NOT delete `_n_face_moments`/`_spatial_moment_tail` — they are width
derivations. #246 should mint ONE typed predicate and route the Kind-A yes/no probes
through it. The width derivations stay; they are not "shape probes standing in for a
predicate", they are honest counts.

## ⭐ THE PLUMBING REALITY (the load-bearing constraint)

`SpatialMomentSpace` is composed onto a bulk-field space ONLY by
`BulkField._compose_spatial_moments` (`orpheus/transport/fields/_bases.py:157-202`), and
that composition is **gated on an EXPLICIT `spatial_moments` parameter** (default `1`
everywhere). It is NOT auto-read from `mesh.scheme.spatial_basis_per_axis` (deliberate
construct-general/select-narrow discipline, `_bases.py:183-194`: "Reading the scheme by
default would silently widen LD field shapes and break LD byte-identity before the
consumers that fill the axis exist — a Pattern-4 violation").

Consequence: a field carries the `SpatialMomentSpace` factor **iff some producer passed
`spatial_moments=per_axis>1`**. Today the ONLY producers that do are the SWEEP/scattering
OUTPUT wraps (`operator.py:1622-1631`; `solver.py:446/642`; `scattering.py:786/1163`;
`operator.py:1620`; `loss_representation.py:2085-2089`). So:
- The matvec/sweep **probe field `psi.bulk`** (`TimedFullField`) at a `loss_action`/`apply`
  entry DOES carry the factor when it came from an LD sweep output → `field.has_spatial_moment_axis`
  is reachable+correct AT THE FIELD BOUNDARY (the `loss_action`/`apply`/`sweep` method entry).
- BUT all 4 raw-probe sites live INSIDE the per-octant walk, operating on octant-sliced
  `np.ndarray` buffers (`probe[oct_idx]`, `Q_cells`, `psi_avg`, `sig`) — **NO typed field
  in scope**, only `self.mesh` (an `SNMesh`, which has `.scheme.spatial_basis_per_axis` +
  `.ndim` but NO field-space factory and NO `spatial_moment_space()` helper).

So the predicate that is actually reachable at the 4 sites is a **mesh/scheme-keyed**
predicate, NOT a field-keyed one. `field.has_spatial_moment_axis` is the RIGHT API at the
method-boundary (`loss_action`, `sweep`, `apply`); a `mesh`/`scheme`-derived predicate
(or the already-existing `self._spatial_moment_tail != ()`) is what the inner-walk sites can use.

## Site-by-site

### Site 1 — `_reframe` guard (`orpheus/sn/sweep_graph.py:100`)
- Code: `if frame_signs is None or arr.shape[-1] != frame_signs.shape[0]: return arr`
- Question: "does `arr` carry the trailing `2^d` moment axis?" (`frame_signs` is the
  `(2^d,)` involution; `frame_signs is None` already encodes DD/Step).
- In scope: a raw `np.ndarray` (`Q_cells`, `psi_avg`, `residual`, the d=1 probe
  `(K,ng,1[,2^d])`) + `frame_signs` (`(2^d,) | None`). NO field, NO mesh, NO scheme —
  `_reframe` is a free function in `sweep_graph.py`, called from `_CellSolve.cell` /
  `_CellResidual.cell` (`sweep_graph.py:901/903/970/973/975`) and the 1-D scan walk
  (`loss_representation.py:2341/2353/2999/3042`). The cell-ops carry `moment_frame_signs:
  (2^d,)|None` but NOT the scheme/mesh.
- **S4 hazard (CONFIRMED real):** the probe `arr.shape[-1] == frame_signs.shape[0] == 2^d`
  collides if a genuine NON-moment array happens to have a trailing axis of length `2^d`.
  Today safe because (a) no flat-source d≥2 LD path exists and (b) the one non-moment array
  `Q_cells` has trailing `n_diag==1 != 2` at d=1. At d=2, `2^d==4`; a flat source whose
  trailing axis is 4 (e.g. collides with an anti-diagonal/level structure) mis-fires a sign
  flip. The typed predicate ELIMINATES this because it keys on the *space's factor presence*
  (a structural fact set by the producer), not on a coincidental ndarray axis length. What
  the predicate needs that `arr.shape[-1]` lacks: the array's TYPED COMPANION (its `Field`
  or, at minimum, the closure's `per_axis>1` fact from the scheme) — the *intent* that this
  axis is the moment axis, not its size.
- Cleanest replacement: this is the HARDEST site. `_reframe` operates on a bare ndarray with
  no companion. Two honest options: (1) thread a boolean `is_moment_valued` (derived ONCE at
  the walk-entry from the field/scheme) into `_reframe` alongside `frame_signs` — then the
  guard is `if not is_moment_valued or frame_signs is None: return arr`; OR (2) keep
  `frame_signs is None` as the DD/Step gate (it already is the typed-ish gate — `None`
  means per_axis==1) and replace ONLY the fragile `arr.shape[-1] != frame_signs.shape[0]`
  half with the threaded boolean. The `frame_signs is None` half is ALREADY a clean typed
  gate (it comes from `octant_moment_frame_signs` which returns `None` for per_axis==1);
  the SECOND half (the ndarray size compare) is the S4-fragile one to replace.

### Site 2 — `_n_face_moments` property (`orpheus/sn/loss_representation.py:311`)
- Code: `per_axis = self.mesh.scheme.spatial_basis_per_axis; return per_axis ** (self.mesh.ndim - 1)`
- Question (the PROPERTY itself): NOT a predicate — it is a COUNT (`2^{d-1}`), Kind B.
- Consumers (read at): `369` (`_inflow_to_moments`: `n = self._n_face_moments; if n == 1: return inflow` — Kind-A predicate `+` width), `1032`, `1060`, `1068` (`if n_face_moments > 1:` capture reduce — Kind A), `1118`, `1134`, `1140` (Kind A), `1292`, `1295`, `1319` (Kind A), `1366`, `1368`, `1386` (Kind A). So `_n_face_moments` feeds BOTH the width (passed as `n_face_moments=` and via `face_moment_tail`) AND the `> 1` boolean.
- In scope: `self` is a `_LossRepresentation` → `self.mesh` (SNMesh) reachable; `self.mesh.scheme.spatial_basis_per_axis` + `self.mesh.ndim` in scope. NO field at these inner sites (the field is unwrapped to `.bulk.values` at the `loss_action`/`sweep` entry above).
- Cleanest replacement: KEEP `_n_face_moments` (it is an honest count). Introduce a sibling
  predicate `_has_spatial_moment_axis` (or reuse `self._spatial_moment_tail != ()`) and
  replace the FOUR `if n_face_moments > 1:` (1068/1140/1319/1386) + the `if n == 1:`
  (369) with it for READABILITY. NOTE: this is Pattern-6 territory — `n_face_moments > 1`
  and `_spatial_moment_tail != ()` are already trivially correct (both derive from the same
  `mesh.scheme.spatial_basis_per_axis`); swapping them buys NAMING clarity, not a safety gain.
  No plumbing needed (mesh+scheme already in scope). Low-risk, low-reward.

### Site 3 — `_moment_broadcast_sigma` (`orpheus/sn/loss_representation.py:493-515`)
- Code: `return sig[..., None] if moment_valued.ndim > sig.ndim + 1 else sig`
- Question: "does `moment_valued` (the probe/Q) carry a trailing moment axis that `sig`
  (`(ng,*spatial)`) lacks?" Kind A (boolean) but answered by a RANK compare, not a size
  compare — so it is S4-SAFER than Site 1 (rank, not coincidental length).
- In scope: two raw ndarrays (`sig`, `moment_valued`). The `+ 1` absorbs the leading
  ordinate axis (`probe`/`Q` is `(N…,ng,*spatial[,2^d])`, one more leading axis than `σ_t`).
  Callers: `698` (`emit.pure_z(oct_idx, q / _moment_broadcast_sigma(operands.sig_t, q))`),
  `789` (`LpC[oct_idx] = _moment_broadcast_sigma(sigma, probe_oct) * probe_oct`). Both are
  inside `_OctantWalk.sweep`/`loss_action`, where `self.mesh` (SNMesh) IS in scope (`sn_mesh
  = self.mesh` at 753) — and the probe-field `psi` enters `loss_action` typed at the method
  boundary (756: `probe = psi.bulk.values` discards the type one line before).
- Cleanest replacement: pass the `is_moment_valued` boolean (derived ONCE at the
  `loss_action`/`sweep` entry from `psi.bulk.has_spatial_moment_axis` or
  `mesh.scheme.spatial_basis_per_axis > 1`) instead of re-deriving via the rank compare.
  Then `_moment_broadcast_sigma(sig, *, is_moment_valued)` → `sig[..., None] if
  is_moment_valued else sig`. Mild plumbing (one bool threaded into the helper). LOW priority
  — the rank compare is already S4-safe (it is NOT the `2^d`-collision hazard; that is Site 1).

### Site 4 — `len(s_axes) > 1` — CONFIRMED RETIRED
- Grep for `len(s_axes)` / `s_axes) > 1` across `orpheus/`: ZERO hits. Retired by S3 (the
  per-cell gate that used it is gone). The `s_axes` tuple is now consumed positionally
  (`for a in range(ndim)`), not length-probed. No residual anywhere. Nothing to do for #246.

## Recommendation: (B) needs plumbing — but it is SMALL, and partly Pattern-6

This is **NOT** a clean "the space is already reachable at all 4 sites" swap. The
`SpatialMomentSpace` factor lives on the TYPED FIELD, and the field is unwrapped to
`.bulk.values` (a bare ndarray) at the `loss_action`/`sweep`/`apply` method ENTRY — every
one of the 4 inner-walk probe sites sees only ndarrays + `self.mesh`. So:

- **At the method boundary** (`loss_action`, `sweep`, `apply`, `loss_action_transpose`): the
  field-keyed `field.has_spatial_moment_axis` (NEW property on `BulkField`, mirroring the
  existing `spatial_moments_per_axis` at `_bases.py:231-253`) is the honest typed query.
  Derive ONE boolean there.
- **Thread that one boolean inward**: into `_reframe` (Site 1 — replaces the S4-fragile
  `arr.shape[-1] != frame_signs.shape[0]` half), into `_moment_broadcast_sigma` (Site 3 —
  replaces the rank compare). The `_CellSolve`/`_CellResidual` cell-ops would carry it
  alongside `moment_frame_signs` (or, since `frame_signs is None ⟺ not moment-valued`, the
  boolean is REDUNDANT with `frame_signs is not None` — see below).
- **Site 2's `> 1` booleans**: replace with `self._has_spatial_moment_axis` (mesh-keyed,
  no plumbing — `self.mesh` already in scope). Pattern-6 naming win only.

⭐ **Elegant shortcut for Site 1 (the actual S4 fix):** `frame_signs` ALREADY encodes the
predicate — `octant_moment_frame_signs(signs, per_axis)` returns `None` iff `per_axis==1`
(DD/Step). So `frame_signs is not None ⟺ this closure carries the moment axis`. The
S4-fragile half (`arr.shape[-1] != frame_signs.shape[0]`) exists ONLY to skip the
involution for a FLAT array (no moment axis) at a MOMENT closure — i.e. the `Q_cells`
matvec-zero / flat-external case. The robust fix: gate on whether the ARRAY is moment-valued
via a threaded boolean (set when the caller KNOWS the array's typed origin), not via its
coincidental last-axis length. The cleanest signature keeps `frame_signs is None` as the
closure-level DD/Step gate (already typed) and replaces the array-level size compare with the
threaded `is_moment_valued`. This is the ONE change that eliminates the S4 hazard.

## Predicate signature + home

- **Field-level (the honest typed query, at method boundaries):**
  `BulkField.has_spatial_moment_axis` — a `@property -> bool` on
  `orpheus/transport/fields/_bases.py` (next to `spatial_moments_per_axis`, 231-253):
  `return self.spatial_moments_per_axis > 1` (which already does the `find_factor(
  SpatialMomentSpace)` query). Trivial, typed, single-sourced. This is the home #246's
  body asks for ("a method on the typed field").
- **Mesh/scheme-level (for the inner-walk sites with no field):** EITHER reuse the existing
  `_LossRepresentation._spatial_moment_tail != ()` / `_n_face_moments > 1`, OR add a thin
  `DiscretizationSchemeBase.is_multi_moment` (`spatial_basis_per_axis > 1`) so the scheme
  owns its own "do I carry slopes?" fact. The scheme-level home is the more honest one (the
  property IS a scheme attribute), but it is OPTIONAL — `mesh.scheme.spatial_basis_per_axis
  > 1` is already a one-liner in scope everywhere it's needed.
- **DO NOT** put a `mesh.spatial_moment_space()` factory on SNMesh just for this — SNMesh has
  no field-space factory today, and minting one solely to query factor-presence is plumbing
  for no safety gain (Pattern 6). The scheme already holds `spatial_basis_per_axis`; that is
  the single source.

## Precise file:line list a method-implementer would touch

- `orpheus/transport/fields/_bases.py:231` — ADD `has_spatial_moment_axis` property
  (delegates to the existing `spatial_moments_per_axis > 1`).
- `orpheus/sn/sweep_graph.py:87-102` — `_reframe`: replace `arr.shape[-1] !=
  frame_signs.shape[0]` with a threaded `is_moment_valued` boolean (signature change:
  add the param). The S4 fix.
- `orpheus/sn/sweep_graph.py:838-926` (`_CellSolve`) + `929-976` (`_CellResidual`) — carry
  the `is_moment_valued` boolean alongside `moment_frame_signs` (or note it is redundant with
  `moment_frame_signs is not None` for the array-IS-the-iterate case, and the only genuinely-
  flat array is `Q_cells` in the matvec → pass `is_moment_valued=False` for that arg).
- `orpheus/sn/loss_representation.py:493-515` (`_moment_broadcast_sigma`) — replace the rank
  compare with a threaded boolean; callers `698`, `789`.
- `orpheus/sn/loss_representation.py:311-320` (`_n_face_moments`) — KEEP (a count). Optionally
  add `_has_spatial_moment_axis` sibling.
- `orpheus/sn/loss_representation.py:369, 1068, 1140, 1319, 1386` — the `== 1`/`> 1`
  predicates: optionally swap to the named predicate (Pattern-6 naming, no plumbing).
- `orpheus/sn/loss_representation.py:2341, 2353, 2999, 3042` — 1-D scan `_reframe` sites
  (same shared `_reframe`); they inherit the signature change. `frame_signs_for(scheme,...)`
  in scope, scheme reachable.
- (optional) `orpheus/sn/spatial/scheme.py` — `DiscretizationSchemeBase.is_multi_moment`
  (`spatial_basis_per_axis > 1`) if the scheme-level home is preferred for the inner sites.

## Anti-surprise notes
- `_n_face_moments` and `_spatial_moment_tail` are WIDTHS, not predicates — do not "collapse"
  them to a boolean; they feed `np.zeros((...,*tail))` allocations and `n_face_moments=` kwargs.
- The `frame_signs is None` half of Site 1's guard is ALREADY typed (None ⟺ per_axis==1).
  Only the second half (the size compare) is the smell. Don't replace the half that's fine.
- The field carries the factor ONLY because the SWEEP OUTPUT WRAP passes `spatial_moments=
  per_axis` (operator.py:1622). The INPUT-side construct-general default is `1` — so a
  hand-built `AngularFlux.from_mesh(values, mesh)` (no `spatial_moments=`) does NOT carry the
  factor even on an LD mesh. The field-level predicate is only trustworthy on fields that came
  through a moment-aware producer; the scheme-level predicate (`spatial_basis_per_axis > 1`)
  is the unconditional truth. Prefer the SCHEME-level predicate for the inner walk; reserve
  the field-level one for boundary readability where the field provenance is known.
