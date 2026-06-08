---
name: composite-fixed-source-api
description: solve_sn_fixed_source composite API review — durable "what good looks like" examples (Union parse-at-boundary, prescribed_inflow Pattern-4) + the SNMesh re-home seam that defeats the mesh-identity guard. LANDED 7c6624d/d87f3b0/bedc394.
metadata:
  type: project
---

`solve_sn_fixed_source(external_source: np.ndarray | TimedFullField)` +
`BoundarySourceSink.prescribed_inflow` + MMS migration. LANDED (`7c6624d`/`d87f3b0`/`bedc394`).
Kept for the durable reference examples + the recurring SNMesh re-home seam.

## "What good looks like" — reinforce these shapes
- **`np.ndarray | TimedFullField` Union dispatch is a CLEAN parse-at-boundary, NOT a twin path:**
  both arms converge on ONE `return TimedFullField(...)`; the array arm just substitutes
  `BoundarySourceSink.zeros_on` for the boundary (vacuum = a VALUE, not a branch). The legacy array
  reads as "bulk-only special case of the composite." isinstance-on-TYPE (not stringly-typed) so the
  type-checker validates both arms. Template for adding a legacy-input arm without forking a path.
- **`prescribed_inflow` masking `view[inflow,:] = arr[inflow,:]` is Pattern-4 textbook:** outflow
  source slots are physically meaningless → unconstructable. Justified as a 3rd constructor (vs
  `from_face_arrays` = faces-optional, vs deferred `from_spec` = recipe; `from_spec` correctly STILL
  deferred per unify-after-two).
- **Reuse-TimedFullField-not-mint-FixedSource is RIGHT:** the public RHS IS the same object SI/Krylov
  consume internally → no twin concept. New constructors rely on dataclass defaults for
  `_history`/`history_depth` = SSOT gain over the retired inner path that passed them at 2 sites.
- Duplication collapse real: the `q_ext_composite` build (was duplicated in `_solve_fixed_source_si`
  + `_solve_fixed_source_krylov`) is now one `_build_fixed_source_rhs`. Retirement = test migration
  honored (`test_mms_prescribed_inflow.py` moved off the `_within_group_triple`+SI bypass).

## The durable seam (promoted to AGENT.md #7 — "SN carves rebuild SNMesh internally")
`solve_sn_fixed_source` builds its OWN `sn_mesh` then re-homes a caller's composite onto it via
`from_mesh(boundary.values.copy(), sn_mesh)`. Load-bearing TODAY (caller's SNMesh ≠ solver's
instance; `TimedFullField.__post_init__` + `BoundaryField._check_partner` enforce mesh-identity).
BUT the `.copy()` re-wrap DEFEATS the mesh-identity guard before it runs — the day SNMesh
construction gets args-sensitive/cached state, a real mismatch is papered over silently. The honest
entry is `from_setup(sn_mesh, composite)` that honors the guard. This shape recurs across SN solver
carves — grep for solvers rebuilding SNMesh from `(materials, mesh, quad)` when a typed object
already carries a mesh.

(The MMS `fixed_source` byte-identical-bundler twin flagged at review was fixed in `bedc394` —
deduped into a shared bundler. The future history-free `FullField` base, when a 2nd timeless consumer
appears, is the deferred factoring; aesthetic-only today.)

Related: [[wave-o-operator-algebra]] (the BoundarySourceSink this builds), [[phase5-windowing-carves]].
