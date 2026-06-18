# #236 Phase 2 B2 review-fixes closeout — perf cache (L16) + production-stamp catcher (vv L11 Mode 11) + test-surrogate dedup (Pattern 2)

Branch `feature/sn-spatial-angular-product` 2026-06-18 (NOT committed; finishing pass on the
already-uncommitted B2 carve). THREE review-found fixes, all BIT-IDENTICAL (0-ULP). Did NOT
touch B3 (`DD.update`/`cell_balance_terms`/the 5th τ consumer) or the geometry-τ producers.

## Context (what B2 already did, unchanged here)
The M-M weighted-diamond constants `c_out = α_out/τ`, `c_in = (1-τ)/τ·α_out + α_in` (an
angular-closure property) are single-sourced onto the angular closure and reach the STATELESS
`DiamondDifference.residual` as DATA on `CellVisit.c_in`/`c_out`, stamped by the single
production site `SNMesh._make_cell_visit` (all 4 `dag_walk` yield paths funnel through it) which
reads the closure's `(N,)` `c_{in,out}_per_ordinate` accessors (the B1 per-global-ordinate gather).
`DD.residual` has NO live production callers (the live matvec reads `cell_contribution`); the
CellVisit.c mechanism is the foundation for the LIVE B3 (`DD.update`/`cell_balance_terms`).

## FIX 1 — perf cache (L16), `orpheus/sn/spatial/pole_angular_closure.py`
The accessors re-ran the full `(N,)` per-level→global `_gather_per_ordinate` on EVERY access →
the per-visit `_make_cell_visit` stamp made the visit-producing loop O(N²·nx). The gather is a
pure permutation of the immutable `_c_*_per_level`, so it is now computed ONCE.
- NEW `PoleAngularClosureBase._build_c_per_ordinate_cache()` (shared on the ABC) gathers BOTH
  constants to `(N,)`, marks them READ-ONLY (`setflags(write=False)` — B1's P3 populator holds a
  reference to the shared view; read-only guards corruption), stores `_c_{in,out}_per_ordinate_cache`.
  Carries a precondition guard (`_c_*_per_level is None → RuntimeError`) — Pattern 4 + narrows the
  type so pyright stays clean (the old `_gather(self._c_*_per_level)` call shape tripped
  `reportArgumentType` on the `tuple|None` annotation; the guard clears both).
- BOTH concrete `__init__`s call it after binding `_c_*_per_level`: `MorelMontryAngularSweep`
  (sphere/cylinder), `IdentityAngularClosure` (slab — neutral zeros). Unbound legacy M-M
  (`sn_mesh=None`) sets the cache to `None` directly (accessors raise the existing error).
- The two properties now just `return self._c_*_per_ordinate_cache` (O(1); same object on
  re-access — verified `is`-identical, no realloc).
- **VALUE-IDENTICAL proven:** in-process `np.array_equal(cache, fresh_regather)` True for c_in AND
  c_out at sphere(GL8,nx10) + multi-level cyl(product 2×4, n_levels=2, nx8) + slab(GL4) — the
  cache is byte-for-byte the old re-gather. Cache arrays read-only (`flags.writeable False`).
- **MEASURED (L16):** in-process microbench, `sphere N=32 nx=200` (6400 visits/sweep), 30 reps:
  BEFORE (re-gather per access) 32.08 ms/sweep → AFTER (cached) 21.91 ms/sweep ≈ **1.46× speedup**.
  (Method: temporarily reverted the two properties to `_gather(self._c_*_per_level)` with a
  `# PERF-PROBE-BEFORE` marker, measured, restored by re-edit — L28, never git checkout.)
  The suggested pytest wall-clock (`test_w1_clamp_silent_on_flat.py`) is fixture-overhead-dominated
  (~6 s either way) so the visit-loop microbench is the load-bearing L16 measurement.

## FIX 2 — committed production-stamp catcher (vv L11 Mode 11)
NEW `tests/sn/sweep/core/test_cell_visit_c_stamp.py` (3 `@pytest.mark.foundation` tests). The
original carve had NO committed test exercising `_make_cell_visit`'s c-stamp: the matvec twin reads
`cell_contribution` directly (never `DD.residual`, the only `CellVisit.c` reader), and the diamond /
cell-balance fixtures stamp visits with a SURROGATE → both structurally blind to the PRODUCER (a wrong
global-ordinate map would ship silent). The catcher walks REAL production `dag_walk`:
- `test_sphere_production_stamp_matches_inline` — sphere (GL8, nx10), all N ordinates, single-level
  (`direction_idx` IS the global ordinate).
- `test_multilevel_cylinder_production_stamp_matches_inline` — cylinder `product(n_mu=2,n_phi=4)`
  → 2 μ-levels, all within-level m per level (exercises the `level_indices[p][m]` per-level→global
  permutation across DISTINCT level blocks).
- `test_slab_production_stamp_is_neutral_zero` — Cartesian (GL4) visits have `c_in==c_out==0.0`.
For EVERY visit asserts `visit.c_*` == constants recomputed INLINE from that visit's OWN
`streaming_terms.{alpha_in,alpha_out,tau_mm}` at 0-ULP via the hand-transcribed surrogate (NOT the
closure's own c — vv L11, no tautology). Mode-8 SAFE: `np.testing.assert_array_equal` + `pytest.fail`
(function calls, fire under `-O`; NO bare assert). The visit-count guards also `pytest.fail`.
**0-ULP confirmed:** in-process `max|visit.c − inline|` EXACTLY 0.000e+00 (c_in + c_out) at all 3
geometries → `assert_array_equal` is the right bar (not `allclose`).
**MUTATION-SELF-CHECK (re-verified):** swapped `c_in`↔`c_out` in `_make_cell_visit` → sphere +
cylinder cases RED (`Max relative difference 0.5`); slab stayed GREEN (neutral c==0 makes the swap
invisible there — correct). Restored by RE-EDIT (L28). Catcher GREEN post-restore.

## FIX 3 — test-surrogate dedup (Pattern 2, unify-after-two)
`test_diamond.py::_c_from_streaming_terms` and `test_cell_balance_for_streaming.py::_visit_c` were
two byte-identical hand-recomputes of `c_out=α_out/τ`, `c_in=(1-τ)/τ·α_out+α_in`. Extracted ONE
shared `tests/sn/sweep/core/_c_surrogate.py::c_from_streaming_terms` (directory-local module — the
surrogate is sweep-core's reference; NOT `tests/sn/_test_helpers` which is SN-wide; NOT conftest
which is hooks/fixtures, not importable helpers). Both files + the new catcher import it. Kept
HAND-TRANSCRIBED (does NOT import the production closure — it is the structurally-independent
reference that keeps the apply-vs-solve round-trip honest, vv L11). Note: the SEPARATE inline
`c_out_scalar`/`c_in_scalar` in `_scalar_to_vectorized_inputs` was LEFT — it builds the vectorized
`angular_denom_term`/`angular_numer_upstream` contributions (a different purpose, not a `(c_in,c_out)`
visit-stamp), so it is not a duplicate of the surrogate.

## GATES (all `-O`, HOST `.venv/bin/python`)
- New catcher: `test_cell_visit_c_stamp.py` 3 PASSED.
- sweep/core: 463 passed / 1 skipped / 4 xfailed (was 460 → +3 catcher tests; dedup'd
  test_diamond + test_cell_balance_for_streaming green).
- Bit-identity (DriftWarning ESCALATES via `-W error::tests.sn.regression._regression_assert.DriftWarning`):
  `tests/sn/sweep/core tests/sn/solve` 523 passed / 1 skip / 4 xf, ZERO DriftWarning. (NOTE: the
  qualified `DriftWarning` filter path is `tests.sn.regression._regression_assert.DriftWarning` —
  the bare `error::DriftWarning` raises `AttributeError: module 'builtins' has no attribute`.)
- `cyl_1g_homogeneous_product_dd_n20` 3.9e-11 NOT grown: the phase-C crosscheck
  `test_phase_c_crosscheck.py -k cyl_1g_homogeneous_product_dd_n20` PASSED (k_inf pinned 1.5).
- Structural anchors: g-adjoint + curvilinear-MMS + tau-producer 23 passed.
- B1 accessor consumer (cache must not change values): `test_sweep_cache.py` 27 passed / 1 skip.
- pyright `pole_angular_closure.py`: 30 errors (was 32 pre-fix; the precondition guard cleared the
  two new `_gather` arg-type errors) — ZERO in my edited regions (cache annotations / helper /
  properties); the 30 are pre-existing #226 `|None` noise on `_alpha_per_level`/`_dr`/`_V` etc.
- Sphinx: build exit 0, zero warnings naming `discrete_ordinates.rst`; B2-fix content rendered
  (4 HTML hits). Pre-existing warnings only (LD `verifies` orphan-eq noise + test_fission SyntaxWarning).

## Changed files
- `orpheus/sn/spatial/pole_angular_closure.py` — cache fields + `_build_c_per_ordinate_cache` ABC
  helper (read-only, precondition guard) + both `__init__`s call it + accessors return cache.
- `tests/sn/sweep/core/test_cell_visit_c_stamp.py` — NEW catcher (3 foundation tests).
- `tests/sn/sweep/core/_c_surrogate.py` — NEW shared hand-transcribed surrogate.
- `tests/sn/sweep/core/test_diamond.py` — import shared surrogate, drop local `_c_from_streaming_terms`.
- `tests/sn/sweep/core/test_cell_balance_for_streaming.py` — import shared surrogate, drop `_visit_c`.
- `docs/theory/discrete_ordinates.rst` — extended the existing `sn-closure-c-on-cellvisit` stub
  TODO with the three review fixes (NO new label — no new verifiable claim; bit-identical opt + a
  catcher + a dedup).

## Sphinx stub
NO new label minted (the finishing fixes add no new verifiable claim). Extended the existing B2
`.. _sn-closure-c-on-cellvisit:` TODO with: the O(1) cache (+ `:meth:` xref to the new helper +
measured 1.46×), the production-stamp catcher (closes the vv L11 Mode-11 gap the carve had), and
the surrogate dedup. archivist DISPATCH for the rich-narrative expansion emitted (followup:false).

## LESSON (vv L11 Mode 11 + L16 + Pattern 2)
A carve that moves a quantity onto a NEW production reader (`_make_cell_visit` stamping
`CellVisit.c`) is NOT proven by a green slow twin if that twin's call graph never reaches the new
reader (the matvec reads `cell_contribution`, never `DD.residual`) AND the only other catchers build
the input packet with a SURROGATE (blind to the producer — mutate the stamp, they stay green because
the surrogate carries the same wrong value on both sides). The committed catcher MUST walk the REAL
production producer (`dag_walk` → `_make_cell_visit`) and recompute the retired formula INLINE from
the visit's own inputs, with a mutation-self-check (swap → red) to prove it is a catcher not a tag.
Separately: a property that re-derives an immutable permutation on every access is an O(1)-cacheable
L16 win — cache once at construction, return read-only so shared-reference consumers can't corrupt
it; the airtight value-identity proof is in-process `np.array_equal(cache, fresh_regather)` at the
exact sphere(single-level)+multi-level-cylinder+slab configs. When caching behind a `tuple|None`
private, a precondition guard (Pattern 4) both documents the post-bind invariant AND narrows the
type so pyright stays clean — preferable to leaving the `reportArgumentType` noise.
