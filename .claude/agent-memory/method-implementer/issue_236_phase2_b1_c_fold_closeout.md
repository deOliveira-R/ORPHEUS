# Issue #236 Phase 2 — τ carve STEP B1 (c_in/c_out fold onto the closure, site P3)

Branch `feature/sn-space-angle-tier2` — 2026-06-17. NOT committed (staged/committed nothing per brief). BIT-IDENTICAL.

## Scope (narrow, one site)

Folded the ONE free seam of the four-site `c_in`/`c_out` single-source-of-truth
violation: the CollisionCache geometry-coefficient populator
`GeometryCoefficients.from_mesh_and_quad` (site **P3**) now sources `c_in`/`c_out`
from `sn_mesh.pole_angular_closure` instead of rebuilding them inline. The other
THREE c-rebuild sites (need CellVisit threading) are LATER dispatches (B2/B3/C) —
out of scope.

## Changed files

### `orpheus/sn/spatial/pole_angular_closure.py` (the accessor)

PUBLIC polymorphic accessor pair added on BOTH:
- The `@runtime_checkable PoleAngularClosure` **Protocol** (lines ~277-300) —
  `c_in_per_ordinate` / `c_out_per_ordinate` property declarations, so consumers
  TYPED against the Protocol (the cache reads `sn_mesh.pole_angular_closure` which
  is typed `PoleAngularClosure`, the Protocol — NOT the ABC) see the members.
  **This was the load-bearing pyright lesson** (mirrors the prior
  `transverse_coupling_is_facewise` closeout): `PoleAngularClosure` is a Protocol,
  `PoleAngularClosureBase` is the concrete ABC; an accessor on the ABC alone is
  invisible to Protocol-typed call sites → "Cannot access attribute" pyright error.
  Mint on BOTH.
- The ABC `PoleAngularClosureBase` (the implementation): `c_in_per_ordinate` /
  `c_out_per_ordinate` `@property` returning the `(N,)` global-ordinate-gathered c
  via the shared `_gather_per_ordinate(per_level)` helper (pure permutation keyed
  on `self.level_indices`, no arithmetic). Both properties + the helper RAISE
  `RuntimeError` on the legacy unbound (`sn_mesh=None`) mode (no precomputed c /
  level partition) — the populator only ever reaches them mesh-bound.

Did NOT reach into `_c_in_per_level`/`_c_out_per_level` (private) from P3 — the
brief's hard constraint. Both concrete subclasses ALREADY store these privates:
- `MorelMontryAngularSweep` (`pole_angular_closure.py:805-816`): per-μ-level
  `(M_p,)` arrays `alpha_out/tau` and `(1-tau)/tau*alpha_out+alpha_in` — the EXACT
  formula P3 inlined.
- `IdentityAngularClosure` (slab/Cartesian, `:1413-1418`): neutral
  `(np.zeros(N),)` for both (α=0, τ=1 ⇒ c=0). Returns the same zeros P3 inlined
  for slab.

Type-annotation hygiene (net-NEGATIVE pyright): added base annotations
`level_indices` / `_c_in_per_level` / `_c_out_per_level` as `tuple | None` (the
broadest contract — M-M unbound stores None). Widened `IdentityAngularClosure`'s
two `_c_*_per_level` + `level_indices` declared annotations from `tuple` to
`tuple | None` to keep the mutable override compatible (runtime values UNCHANGED —
Identity always binds the non-None tuple). NO subclass c VALUES touched.

### `orpheus/sn/spatial/sweep_cache.py` (the P3 rewire)

- Dropped the now-orphaned `alpha_in`/`alpha_out` `(N,)` allocations (`:236-237`)
  and their `st0.alpha_in`/`st0.alpha_out` reads (`:299-300`). `tau` STAYS (still
  feeds `tau_inv` + `mm_a_in_coeff`). `StreamingTerms` fields LEFT alone (their
  retirement is Step C).
- Replaced the inline `c_out = alpha_out/tau` / `c_in = (1-tau)/tau*alpha_out+alpha_in`
  (`:309-310`) with `closure = sn_mesh.pole_angular_closure; c_out =
  closure.c_out_per_ordinate; c_in = closure.c_in_per_ordinate`. Dispatch is by
  closure TYPE (no `if coord ==` branch). 12-line WHY comment citing Cardinal
  Rule 2 + Pattern 7 + the bit-identity argument.

### `docs/theory/discrete_ordinates.rst` (Sphinx stub)

NEW section `.. _sn-closure-c-constants-owned:` "c_in / c_out are angular-closure
constants — Step B1 (one site folded)" placed right after the Step-A
`sn-tau-closure-owned` todo (which already forward-references "Step B...
consolidates the four-site c_in/c_out duplication"). Contains: a `.. math::` for
both c constants, `:mod:`/`:meth:`/`:attr:` cross-refs to the accessor + the
folded consumer + the anchor gate, `:eq:`mm-weights`/`dd-mm-closure-constants`
cross-refs, and the archivist `.. todo::` marker. Build rc=0; label + heading +
math rendered into HTML (grep-confirmed). Nexus rebuilt during build.

## Bit-identity verdict — PROVEN 0-ULP

The argument: closure computes c from closure-τ (0-ULP == geometry-τ, pinned by
Step-A Leg-1 gate `test_tau_producer_equivalence.py`) and the SAME α the populator
read → closure per-level c == inline c bit-for-bit; per-level→(N,) gather is a
pure permutation. Confirmed three ways:

1. **Anchor gate** `test_cache_populator_matches_cell_balance_terms` (sphere) PASS
   at rtol=1e-14 — pins cache `denom` (carries `(ΔA/w)·c_out`) to
   `cell_balance_terms`.
2. **In-process c-array equality** (`np.array_equal`, max diff EXACTLY 0.000e+00):
   - sphere (nx,N) ∈ {(8,4),(20,8),(12,6)} — c_out AND c_in 0-ULP.
   - cylinder at the EXACT regression configs (`product_2x4` n_levels=2 N=8;
     `LS4` n_levels=2 N=24) — c_out AND c_in 0-ULP (the multi-level gather is
     byte-identical).
3. **DriftWarning-escalating run** `sweep/core + solve` = 520P/1skip/4xf, ZERO
   DriftWarning, no golden moved. The `cyl_1g_homogeneous_product_dd_n20`
   snapshot drift stays at the pre-existing 3.90e-11 (the brief's named acceptable
   boundary) — NOT grown.

## Verification (all `-O`, paste-back in the response)

- V1 cache-equivalence: `27 passed, 1 skipped` (anchor test isolated: `1 passed`).
- V2 `sweep/core + solve` (DriftWarning escalates): `520 passed, 1 skipped, 4 xfailed`, no drift.
- V3 structural anchors (g-adjoint + curvilinear MMS + tau-producer): `23 passed`.
- V4 curvilinear sweep set (route around unified_matvec): `168 passed, 33 deselected`.
- Conformance: `test_pole_angular_closure.py` 13P; `test_discretization_scheme_protocol.py` 17P
  (no `@runtime_checkable` breakage from the new Protocol members); `tests/sn/spatial` 76P.
- Operators set: `7 failed, 505 passed` — the 7 are PRE-EXISTING (3 sphere
  vacuum_bulk_1d[*-SPH] + 2 `Face 'ymin' mu_y` + 2 sphere apply snapshots; the
  documented branch set across prior D5b sessions). PROVEN independent of my
  change: the sphere matvec receives byte-identical c (proof #2 above).
- pyright on both changed files: 31 errors (DOWN from 36 pre-edit baseline) — all
  31 are pre-existing #226 `| None` noise on `_alpha_per_level`/`_dr`/`_V` +
  sweep_cache:280 `face_area_inner+face_area_outer`; ZERO reference my accessors,
  the P3 rewire, the override, or attribute-access. Net-NEGATIVE (removed the 5
  errors my first-draft base annotation introduced).

## Owed

- B2/B3/C: the remaining 3 inline c-rebuild sites (need CellVisit threading) +
  `StreamingTerms.alpha_in`/`alpha_out`/`tau_mm` field retirement (Step C) +
  geometry-side τ producer retirement (Step B).
- archivist: rich-narrative expansion of the `sn-closure-c-constants-owned` stub
  (DISPATCH emitted, followup:false).

## LESSON

When a consumer is typed against a Protocol (not the concrete ABC), a new
polymorphic accessor must be minted on BOTH the Protocol AND the ABC — the
Protocol declares the contract the consumer's type sees, the ABC supplies the one
shared implementation. Minting on the ABC alone leaves the Protocol-typed call
site with a pyright "unknown attribute" (the SAME mechanism as the prior
`transverse_coupling_is_facewise` Protocol-trait mint). For a per-level→global
gather fold, the airtight bit-identity proof is the in-process `np.array_equal`
on the gathered `(N,)` arrays at the EXACT snapshot configs (sphere single-level
AND cylinder multi-level) — it isolates the producer-side fold from the
pre-existing matvec-snapshot reds that a whole-suite run conflates.
