---
name: issue-236-phase2-b3-tau-c-live-closeout
description: "#236 Phase 2 B3 — the LIVE τ fold: DD.update + scan + cell_balance_terms consume closure-owned τ/c off CellVisit instead of StreamingTerms.tau_mm; BIT-IDENTICAL 0-ULP; the precondition for Step C (geometry-τ retirement)"
metadata:
  type: project
---

# #236 Phase 2 B3 — the LIVE τ fold (closeout)

Branch `feature/sn-spatial-angular-product`, MAIN checkout, 2026-06-18.
NOT committed — left staged-or-unstaged for the main agent.

**Why:** B3 is the step that makes the LIVE sweep + scan + matvec-solve
paths CONSUME the angular-closure-owned τ (and the derived c_in/c_out)
instead of the geometry-owned `StreamingTerms.tau_mm`. After B3, NOTHING
in the live sweep / scan / matvec reads `StreamingTerms.tau_mm` — the
precondition Step C verifies before DELETING the geometry-side τ producer.

**How to apply:** B3 is the 3rd-of-4 c-fold sites + the 5th τ consumer.
The mechanism: closure mints a public `tau_per_ordinate` (mirroring B1's
c-accessors) on BOTH the Protocol AND the ABC; `CellVisit` gains a `tau`
field (DEFAULT 1.0, NOT 0.0 — divide-by-zero landmine); `_make_cell_visit`
stamps it; `DD.update` reads `visit.tau` + passes `visit.c_in/c_out` to
`cell_balance_terms` (now a pure consumer, no τ read); `sweep_cache` sources
τ from the closure for the scan split. Step C is the next dispatch.

## The 6 production seams (all bit-identical 0-ULP)

1. **`pole_angular_closure.py`** — `tau_per_ordinate` accessor + cache.
   - Base ABC annotation `_tau_per_level: tuple|None` (was unannotated).
   - `_tau_per_ordinate_cache: np.ndarray|None` field.
   - RENAMED `_build_c_per_ordinate_cache` → `_build_per_ordinate_cache`
     (now gathers 3 constants: c_in, c_out, τ; precondition checks all 3).
   - `tau_per_ordinate` property on the ABC (raises when unbound).
   - `tau_per_ordinate` property declared on the `@runtime_checkable`
     Protocol too (contract completeness — SAME both-sites mint as B1's c).
   - M-M unbound branch: `_tau_per_ordinate_cache = None`.
   - Both bound `__init__`s: renamed `_build_per_ordinate_cache()` call.
   - **OVERRIDE-CONTRAVARIANCE FIX (the pyright catch):** widened Identity's
     `_tau_per_level` annotation `tuple → tuple|None` to match the new base
     annotation (EXACTLY the widening B1 applied to `_c_*_per_level`). The
     bare `tuple` would be a `reportIncompatibleVariableOverride` (mutable
     override invariant). Runtime VALUE unchanged (`np.ones(N)`).
2. **`scheme.py`** — `CellVisit.tau: float = 1.0` (after c_out). DEFAULT 1.0
   = neutral M-M weight; 0.0 would divide-by-zero in DD.update's recurrence.
3. **`geometry.py`** — `_make_cell_visit` stamps
   `tau=float(closure.tau_per_ordinate[global_ordinate])`.
4. **`cell_balance.py`** — `cell_balance_terms` signature gains
   `*, c_in: float, c_out: float`; DELETED `tau = st.tau_mm` + the inline
   c rebuild; consumes passed-in c. No longer reads τ. (Return still carries
   c_in/c_out fields, now = passed-in values; `test_ordinate_scan` reads them.)
5. **`diamond.py`** — `DD.update`: passes `c_in=visit.c_in, c_out=visit.c_out`
   to `cell_balance_terms`; angular recurrence `tau = visit.tau` (was
   `visit.streaming_terms.tau_mm`).
6. **`sweep_cache.py`** — DELETED `tau = np.empty(N)` + `tau[gn] = st0.tau_mm`;
   sourced `tau = closure.tau_per_ordinate`; KEPT `tau_inv = 1/tau`,
   `mm_a_in_coeff = (1-tau)/tau` (scan split, L16 perf-cache). `st0` kept for
   `mu_start`. `loss_representation.py:3269` (scan recurrence) UNTOUCHED —
   re-sourced the data behind `geom.tau_inv`/`mm_a_in_coeff`.

## Caller audit for `cell_balance_terms` (7 sites)

1 production (`diamond.py:205` → `c_in=visit.c_in, c_out=visit.c_out`).
6 test sites retargeted:
- `test_diamond.py:1305,1347` (visits stamped by helpers) → `c_in=visit.c_in,
  c_out=visit.c_out`.
- `test_ordinate_scan.py:90` (`_affine_coefficients_from_visits`, unstamped
  visits) → `c_from_streaming_terms(visit.streaming_terms)`.
- `test_cell_balance_for_streaming.py:158,198` (raw `st`) →
  `c_from_streaming_terms(st)`.
- `test_sweep_cache.py:477` (dag_walk visit, stamped) → `visit.c_in/c_out`.

## TEST-FIXTURE RIPPLE (the failures that surfaced — NOT a prod break)

`DD.update` now reads `visit.c_in/c_out` + `visit.tau` off the visit.
Test fixtures that built `CellVisit` UNSTAMPED (defaulting c=0.0, τ=1.0) and
relied on `update` rebuilding c internally DIVERGED. Root cause = the
fixtures, NOT production (in-process probe proved prod is 0-ULP). Fixed:
- `test_ordinate_scan.py`: the 3 visit-builders (slab/sphere/cylinder) now
  stamp `c_in,c_out = c_from_streaming_terms(st)` (slab → 0.0). 31 dual-view
  failures → green.
- `test_diamond.py`: `TestBitIdenticalCurvilinear` (2) + `TestCylindrical
  Degenerate` (2) synthetic visits stamp `c_in=ref_c_in, c_out=ref_c_out,
  tau=st.tau_mm` (those tests compute ref c inline + the angular ref uses τ).
- The 4 SLAB visits in test_diamond.py (`face_area_downstream=1.0`) LEFT
  unstamped — slab's neutral c=0/τ=1 defaults ARE correct.

## Gates (verbatim summary lines)

Gate A — live-sweep + scan bit-id snapshots (DriftWarning-escalating):
```
588 passed, 2 skipped, 4 xfailed, 2 warnings in 84.81s
```
(`tests/sn/sweep/core tests/sn/spatial tests/sn/sweep/cartesian_2d
-W error::...DriftWarning` — ZERO drift escalation.)

Gate B — solve-level snapshots (eigenvalue + fixed-source, #212 deselected):
```
60 passed, 1 warning in 77.68s
```

Gate C — Leg-1 τ producer-equivalence + B2 c-stamp catcher:
```
12 passed, 1 warning in 0.17s
```

Gate D — structural anchors (matvec twin, curvilinear MMS, g-adjoint):
```
7 failed, 507 passed, 4 skipped, 1 xfailed, 15 warnings in 5.31s
```
The 7 FAILURES ARE PRE-EXISTING — proven by re-running them at HEAD
(`5418e58`) in a read-only git worktree: identical 7 fail, same magnitudes
(3× `test_vacuum_bulk_bit_identical_1d[*-SPH]` + 2× `Face 'ymin' mu_y` +
2× `test_sphere_{1g,2g}_apply_bit_identical`). NONE of the 7 test files are
in my changeset. The mismatch (O(0.8), 100% elements, 3.4 relative) is
structurally incompatible with a 0-ULP change. The B3-relevant anchors
(curvilinear MMS + g-adjoint) pass in isolation: `14 passed`.

Gate E — P3 cache-populator gate (`tau_inv`/`mm_a_in_coeff` match):
```
27 passed, 1 skipped, 1 warning in 0.31s
```
(`test_cache_populator_matches_cell_balance_terms` pins closure-sourced
`tau_inv`/`mm_a_in_coeff`; green.)

Retargeted-file confirmation under `-O`:
```
141 passed, 1 skipped, 1 warning in 0.43s
```

## In-process bit-identity proof (0-ULP, 3 quantities × 4 configs)

`abs(visit.tau − st.tau_mm)`, `abs(closure.tau_per_ordinate − st.tau_mm)`,
and `np.array_equal(cache.tau_inv/mm_a_in_coeff, geometry-τ-derived)` all
EXACTLY 0.000e+00 at sphere(nx=8,GL8)+(12,GL6) single-level AND
cylinder(nx=8,product 4×4)+(12,4×8) multi-level.

## Coverage check (vv Mode-11)

- (i) curvilinear `DD.update` slow path (reads `visit.tau`/`visit.c`):
  EXERCISED by `test_diamond.py::TestBitIdenticalCurvilinear::
  test_spherical_{outward,inward}_bit_identical` + `TestCylindrical
  Degenerate` + `test_ordinate_scan.py::TestDualViewContracts` (sphere/cyl)
  — all GREEN, all call `DD.update` with non-flat angular upstream.
- (ii) CumprodScan fast path with non-flat-in-angle flux: EXERCISED by the
  sphere/cyl solves in gate B. INSTRUMENTED `GeometryCoefficients.
  from_mesh_and_quad` during a 2-group heterogeneous-scatter sphere
  `solve_sn_fixed_source` — cache built with NON-TRIVIAL closure-sourced
  `tau_inv` (max|tau_inv−1|=1.55) + `mm_a_in_coeff` (1.55); the scan
  recurrence `geom.tau_inv·ψ_avg − geom.mm_a_in_coeff·ψ_in` consumes them.

## pyright (touched files)

ZERO new errors from B3. HEAD baseline (worktree) = 62; working tree = 54.
The one error my BASE annotation would have introduced
(`pole_angular_closure.py:1712 _tau_per_level reportIncompatible
VariableOverride`) I FIXED by widening Identity's annotation (net −1). All
54 remaining are pre-existing #226 `StreamingTerms.<field>: float|None`
wrong-rooting noise on `alpha_*`/`abs_mu`/`_dr`/`_V` operand chains — none
mention `tau_per_ordinate`/`tau_per_level`/`CellVisit.tau`/
`_build_per_ordinate`/override.

## Sphinx

NEW stub `sn-tau-c-on-cellvisit-live` in `docs/theory/discrete_ordinates.rst`
(after B2's `sn-closure-c-on-cellvisit`, before "Substituting the WDD
Closure"). `:eq:`/`:meth:`/`:attr:`/`:mod:` xrefs + archivist TODO. Also
fixed the B2 section's now-stale `:meth:_build_c_per_ordinate_cache` →
`_build_per_ordinate_cache`. Build exit 0, no new warnings; label+section in
HTML + Nexus.

## Deviations from brief

- Brief §1.2 said `CellVisit` has `c_in/c_out` after `face_area_downstream`;
  field default for `face_area_downstream` is actually 0.0 (brief's "1.0" was
  the SLAB docstring value, not the field default) — no impact, I added
  `tau: float = 1.0` after `c_out` as specified.
- Brief §2 SEAM 1 said rename + add τ cache "the only change" in both bound
  `__init__`s — TRUE, but pyright surfaced the Identity `_tau_per_level`
  override widening as a REQUIRED extra edit (the base-annotation side-effect
  the brief's SEAM-1 list implied via "add base annotation" but didn't spell
  out). Applied the B1-precedent widening.
- The 35 test failures in gate A/D's first run were FIXTURE inconsistencies
  (unstamped visits), repaired in-place per the brief's SEAM-4 caller-audit
  spirit — not in the original 6-test-site list because they construct
  `CellVisit` directly (not via `cell_balance_terms`).

## OWED

- **Step C** (next dispatch): retire the geometry-side τ producer
  (`spherical_streaming:681-688` / `cylindrical_streaming:798-815` /
  `slab_streaming:495`) + `StreamingTerms.tau_mm` / `alpha_in`/`alpha_out` /
  `ReducedStreamingOperator.tau_mm` / `tau_mm_per_level`. B3 left them in
  place (bit-id precondition); Step C verifies NOTHING live reads them.
- archivist DISPATCH emitted (followup: false) for the rich narrative.
- #248 (orphaned Protocol) — track when C lands.

## LESSON

Moving a CONSUMER from "rebuilds X internally" to "reads X off a data
packet" makes EVERY fixture that built the packet UNSTAMPED a latent red —
even when the production fold is provably 0-ULP. The B2 closeout already
documented this for `residual`; B3 hit it again for `update` because
`update`'s consumers (the dual-view + bit-identical-curvilinear fixtures)
built visits without c/τ and relied on the OLD internal rebuild. Stamp the
fixture (the surrogate / the inline ref) the moment a scheme method starts
reading a field off the visit. Separately: a "new base annotation" on an
ABC whose concrete subclass binds the same name with a NARROWER type
(non-`|None`) is a `reportIncompatibleVariableOverride` — widen the
subclass annotation to match (value unchanged), the exact B1-precedent the
brief's SEAM-1 list implied. The airtight bit-id proof remains in-process
`abs(diff)==0` on the stamped `visit.tau` + `np.array_equal` on the cache
split at the EXACT sphere(single-level)+cyl(multi-level) configs, isolating
the fold from the pre-existing operators-matvec xfails a whole-suite run
conflates.
