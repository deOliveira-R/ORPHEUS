# Issue #236 Phase 2 Step C — geometry-τ retirement: L20 dependency audit

**Branch**: `feature/sn-spatial-angular-product` · **HEAD**: `002844a` ·
**Checkout**: MAIN `/Users/rodrigo/git/nuclear/ORPHEUS` (not a worktree) ·
**Date**: 2026-06-18 · **Scope**: read-only L20 audit (this IS the Step-C plan).

Steps A/B1/B2/B3 relocated the Morel–Montry τ and the derived c_in/c_out
constants onto the angular CLOSURE. The live sweep/scan/matvec path now reads
τ/c from the closure (`CellVisit.tau`/`.c_in`/`.c_out`, stamped in
`_make_cell_visit`). The geometry-side τ producer + the
`StreamingTerms.tau_mm`/`alpha_in`/`alpha_out` fields are the retirement
targets. This audit proves which are orphaned in PRODUCTION and which still
have a load-bearing reader (notably: the test ORACLES still read geometry-τ).

---

## GO / NO-GO — "is geometry-τ now fully orphaned in production?"

**GO (production).** Of the 13 `.tau_mm` lines in `orpheus/`, exactly **2 are
code**, both in `reduced_operator.py` — `:505` (producer assert) and `:519`
(the dead-bake into `StreamingTerms`). The other 11 are docstrings/comments.
There are **ZERO live reads** of `StreamingTerms.tau_mm` / `.alpha_in` /
`.alpha_out` anywhere in the sweep / scan / matvec / scheme path. The live τ/c
flow is closure-sourced (`geometry.py:1330-1332` reads
`closure.{c_in,c_out,tau}_per_ordinate`; `diamond.py:210/211/245` reads
`visit.c_in`/`.c_out`/`.tau`). **The geometry-side τ bake is dead production
code.**

**NO-GO (as a blind delete) — geometry-τ is still load-bearing as a TEST
ORACLE.** Both regression-floor tests use the GEOMETRY producer as their
structurally-independent (vv L11) reference:
- `test_cell_visit_c_stamp.py::_assert_visit_stamp_matches_inline` (`:94-96`)
  cross-checks the stamped `visit.tau` against `st.tau_mm` (the geometry bake)
  — deliberately distinct from the closure producer to avoid a tautology
  (comment `:42-45`).
- `test_tau_producer_equivalence.py` `test_sphere_closure_tau_equals_geometry_factory_0ulp`
  (`:87`) and `test_cyl_closure_tau_equals_geometry_factory_0ulp` (`:148`) read
  `spherical_streaming(...).tau_mm` / `cylindrical_streaming(...).tau_mm_per_level`
  as the reference half of a 0-ULP equality.

⇒ **Step C is a surgical FIELD-excision plus a TEST-ORACLE migration, not a
blind delete.** The oracles must be re-pointed at the structurally-independent
`morel_montry_weights` (contamination.py:156), which the OTHER floor legs
(`test_sphere_tau_matches_independent_reference` `:109`,
`test_cyl_tau_clamp_is_the_only_difference_from_reference` `:175`) already use,
BEFORE the geometry-τ producer/fields are removed. Migrate-then-delete.

---

## Per-target three-surface table

Legend: (a) production callers · (b) test callers · (c) internal-to-orphan.

### T1 — geometry τ PRODUCERS (`reduced_operator.py`)

| Symbol | (a) prod | (b) tests | (c) internal-to-orphan |
|---|---|---|---|
| sphere τ block `:681-688` (`mu_edge` loop + `tau_mm` array) | feeds `:699 tau_mm=tau_mm` field only | `test_reduced_operator.py::test_tau_mm_bit_identical` `:115`, `test_alpha_dome_bit_identical_across_N` `:129`; `test_tau_producer_equivalence.py:87` (oracle) | `mu_edge` array is SHARED — `:700 mu_start=float(mu_edge[0])` reads `mu_edge[0]` (= −1.0). **Surgical**: keep `mu_edge[0]` for `mu_start`, drop the `mu_edge[1:]` fill + `tau_mm` loop. |
| cylinder τ block `:808-815` (`tau`/`eta_edge` loop) | feeds `:826 tau_mm_per_level=…` field only | `test_reduced_operator.py::test_tau_mm_per_level_bit_identical` `:174`, `test_full_per_level_hash_equality` `:184`; `test_tau_producer_equivalence.py:148` (oracle) | shares the per-level loop with `mu_start_per_level.append(-sin_theta)` (`:802`) — `sin_theta` (`:801`) is ALSO consumed by `eta_edge[0]`/`[M]` (τ-only). **Surgical**: keep the `mu_start_per_level` append; drop `eta_edge` + `tau` + the `tau_mm_per_level.append`. |
| slab synthetic `tau_mm=1.0` `:495` | (the dead-bake, slab arm) | `test_reduced_operator.py::test_slab_streaming_terms_neutral_curvature:288`; `test_discretization_scheme_protocol.py:401` | none — single kwarg in the slab `StreamingTerms(...)`. |
| the bake accessor `streaming_terms()` `:519` (sphere) / `:554` (cyl) / `:495` (slab) + assert `:505`/`:534` | **DEAD-BAKE** (no downstream reader of the field) | (consumed indirectly via T3 test reads) | the `assert self.tau_mm is not None` (`:505`) + `assert self.tau_mm_per_level is not None` (`:534`) retire WITH the field. |

**Surgical-excision hazard (T1)**: the τ producers are INTERLEAVED with
still-needed outputs. `spherical_streaming` returns `alpha_half`, `redist_dAw`,
`mu_start`, `face_areas`, `delta_A` — all live. `mu_start=float(mu_edge[0])`
shares the `mu_edge` array with the τ loop. `cylindrical_streaming` shares the
per-level loop body (`sin_theta`, `mu_start_per_level`) with the τ computation.
⇒ **Whole-function deletion is WRONG.** Excise only the τ statements; keep the
α-dome / redist / mu_start / face-area machinery. (α-dome, redist_dAw, mu_start
are consumed by the LIVE Carlson seed `CarlsonSweepContext` and by
`streaming_terms()`'s still-live fields — do not touch them.)

### T2 — `ReducedStreamingOperator.tau_mm` (`:415`) / `.tau_mm_per_level` (`:426`)

| Surface | Finding |
|---|---|
| (a) prod | written by the producers (T1); read at `:505/:519` (sphere bake) and `:534/:545/:554` (cyl bake) — the dead-bake only. No other production reader. |
| (b) tests | `test_reduced_operator.py` `:117/:118/:129` (tau_mm), `:176/:178/:199/:229/:232` (tau_mm_per_level + None-on-slab); `test_tau_producer_equivalence.py:87/:148` (oracle). |
| (c) internal | the producer τ blocks (T1) exist only to populate these; the bake reads (T1 accessor) exist only to forward them into `StreamingTerms.tau_mm`. Retire together. |

### T3 — `StreamingTerms.tau_mm` (`:295`), `.alpha_in` (`:283`), `.alpha_out` (`:292`)

| Surface | Finding |
|---|---|
| (a) prod | **ZERO live reads.** Verbatim grep below: only `reduced_operator.py:519/:554/:495` WRITE `tau_mm=`; `:517/:518/:552/:553/:493/:494` WRITE `alpha_in/out=`. No `st.tau_mm` / `st.alpha_*` READ in any sweep/scan/scheme code. The live consumer `cell_balance_terms` (`cell_balance.py:337/352`) takes `c_in`/`c_out` as ARGUMENTS — it explicitly "does not read here from `st.alpha_*` / `st.tau_mm`" (`:286`). `diamond.py` reads `visit.c_*`/`visit.tau`, NOT `st.*`. |
| (b) tests | reads of `st.tau_mm`/`st.alpha_in`/`st.alpha_out`: `test_reduced_operator.py` `:286-288/:303-305/:326-328/:428-431`; `test_discretization_scheme_protocol.py` `:156/:160/:164/:184/:398-401/:411-414`; `test_diamond.py` `:381/:402-404/:480-481/:495/:505/:567-569/:587-588/:593/:605/:654-656/:682/:970-972`; `test_cell_balance_for_streaming.py` `:119-121` (+constructor `:77-79/:93-95`); `test_cell_visit_c_stamp.py:95` (the ORACLE). |
| (c) internal | these three are the bake's destination; they exist only to carry geometry-τ/α into a packet nobody reads in production. Retire together with T1/T2. The OTHER `StreamingTerms` fields (`chord_length`, `mu`, `face_area_inner/outer`, `delta_A_over_w`, `mu_start`, `volume`, `abs_mu`) are LIVE — `streaming_terms()` and `CellVisit.streaming_terms` stay. |

### T4 — `StreamingTerms.alpha_in` / `.alpha_out` orphan status (VERIFY post-B3)

**CONFIRMED orphaned in production.** The B-audit claim holds: the SOLE live
consumers of α_in/α_out were the c_in/c_out rebuild sites, now retired by
A/B1/B2/B3. The live c_in/c_out are computed in the CLOSURE from the closure's
OWN per-level α (`pole_angular_closure.py:1106-1109`:
`alpha_in = alpha[:M_p]` / `alpha_out = alpha[1:M_p+1]` — these are LOCAL names
bound from the closure's `alpha` array, NOT `StreamingTerms.alpha_*`). Zero
`.alpha_in`/`.alpha_out` attribute reads on any `StreamingTerms`/`st` instance
in `orpheus/` (the `derivations/` `alpha_in`/`alpha_out` hits are an UNRELATED
trajectory-resolvent Green's-function API — different namespace, not a target).

### T5 — legacy `__call__`-arg `tau_mm` on the closure (SEPARATE surface — SURVIVES Step C)

| Site | Finding |
|---|---|
| Protocol `pole_angular_closure.py:305` (sig `:320`, doc `:236/:310`) | the `__call__(..., tau_mm=...)` Protocol signature |
| ABC abstractmethod `:502` (sig `:577`) | abstract `__call__` |
| spherical `__call__` body `:1391`/`:1468`/`:1512-1526` | consumes `tau_mm` ARG (asserts ndarray, threads into the per-cell recurrence) |
| cyl `__call__` body `:1435-1475`/`:1533-1552` | consumes per-level `tau_mm[p]` ARG |
| Identity `__call__` `:1777/:1783` | `del ... tau_mm ...` (ignores) |

**Who CALLS `closure.__call__(..., tau_mm=...)`?** Search for the callers:
this is the UNBOUND-legacy `MorelMontryAngularSweep(sn_mesh=None)` invocation
path where τ is passed positionally as a runtime arg (the closure is NOT
mesh-bound, so `tau_per_ordinate` would raise `:565`). **This is a DIFFERENT
surface from the geometry-τ producer.** The `__call__`-arg `tau_mm` is the
closure's own runtime parameter; the geometry producer's `StreamingTerms.tau_mm`
is what Step C deletes. **T5 does NOT fold into Step C** — it survives unless a
separate decision retires the unbound-legacy closure call path. (Recommend:
file/track as a follow-up; out of Step-C scope. Confirm no PRODUCTION caller
passes `tau_mm=` to `__call__` before any future T5 retirement — the bound path
uses `tau_per_ordinate`, not the arg.)

---

## Load-bearing verbatim grep (L12) — `\.tau_mm\b` in `orpheus/`

```
$ grep -rn '\.tau_mm\b' orpheus/
orpheus/sn/geometry.py:1317:        ``st.alpha_*`` / ``st.tau_mm`` (the inline formula the former   # COMMENT
orpheus/sn/spatial/diamond.py:204:        # them from st.alpha_* / st.tau_mm.                              # COMMENT
orpheus/sn/spatial/diamond.py:240:        # the geometry-owned streaming_terms.tau_mm (which Step C retires); # COMMENT
orpheus/sn/spatial/diamond.py:315:        # them from ``st.alpha_*`` / ``st.tau_mm`` — that inline formula    # COMMENT
orpheus/sn/spatial/pole_angular_closure.py:310:  ... StreamingTerms.tau_mm`` ...                                       # DOCSTRING
orpheus/sn/spatial/pole_angular_closure.py:399:  # ``StreamingTerms.tau_mm``. ...                                      # COMMENT
orpheus/sn/spatial/pole_angular_closure.py:485:  ... StreamingTerms.tau_mm``) via :attr:`tau_per_ordinate`.            # DOCSTRING
orpheus/sn/spatial/pole_angular_closure.py:560:  ... StreamingTerms.tau_mm`` (bit-identical ...                       # DOCSTRING
orpheus/sn/spatial/pole_angular_closure.py:1065: # ... (``reduced.tau_mm`` /                                          # COMMENT
orpheus/sn/spatial/cell_balance.py:286:    here from ``st.alpha_*`` / ``st.tau_mm`` — this helper does not read  # DOCSTRING
orpheus/sn/spatial/scheme.py:273:  ... streaming_terms``\ ``.tau_mm`` (which Step C will retire); ...     # DOCSTRING
orpheus/geometry/reduced_operator.py:505:            assert self.tau_mm is not None                              # CODE (producer assert)
orpheus/geometry/reduced_operator.py:519:                tau_mm=float(self.tau_mm[direction_idx]),                # CODE (dead-bake into StreamingTerms)
```

Counts: 13 total lines · 6 files. **2 CODE** (`reduced_operator.py:505`,
`:519`) · 11 docstring/comment. Cyl analogue `:534`/`:554` and slab `:495` use
`.tau_mm_per_level`/literal `1.0` so don't match `\.tau_mm\b` exactly but are
the same dead-bake (see T1). **Zero live sweep/scan/matvec readers.**

---

## StreamingTerms constructor-call-site ripple (field deletion)

Dropping `tau_mm`/`alpha_in`/`alpha_out` as dataclass fields ripples to every
`StreamingTerms(...)` keyword call-site that passes them:

```
$ grep -rn 'StreamingTerms(' orpheus/ tests/ derivations/
orpheus/geometry/reduced_operator.py:487  (slab)   — drop tau_mm=1.0/alpha_in=0.0/alpha_out=0.0
orpheus/geometry/reduced_operator.py:509  (sphere) — drop tau_mm=/alpha_in=/alpha_out=
orpheus/geometry/reduced_operator.py:546  (cyl)    — drop tau_mm=/alpha_in=/alpha_out=
tests/sn/sweep/core/test_diamond.py:561/648/964    — pass tau_mm=st_real.tau_mm/alpha_in=/alpha_out=
tests/sn/sweep/core/test_cell_balance_for_streaming.py:71/87 — pass tau_mm=/alpha_in=/alpha_out=
```

The fields are keyword-only-with-defaults (`tau_mm: float | None = None`,
etc.), so deletion is clean for any call-site that OMITS them; but the 3
production constructors + 5 test constructors PASS them and must drop the
kwargs. Field deletion is mechanical (defaults, no positional reliance). The
test constructors at `test_diamond.py` / `test_cell_balance_for_streaming.py`
must be migrated/retired in the same step (they construct a `StreamingTerms`
with τ/α to drive an INLINE c-recompute the production no longer does — those
inline-recompute tests pin RETIRED math and should be re-pointed at the closure
c/τ or deleted).

---

## Retirement ORDER (leaves-first, L20)

1. **Migrate test ORACLES off geometry-τ FIRST** (de-risk; keeps the floor
   green throughout):
   - `test_cell_visit_c_stamp.py::_assert_visit_stamp_matches_inline` `:94-96`:
     replace the `st.tau_mm` reference with `morel_montry_weights(quad, coord)`
     (contamination.py:156) gathered to the visit's global ordinate (sphere
     direct; cyl via `level_indices`). Keep it structurally independent of the
     closure stamp (do NOT reference `tau_per_ordinate`).
   - `test_tau_producer_equivalence.py::test_sphere_closure_tau_equals_geometry_factory_0ulp`
     `:77` and `::test_cyl_closure_tau_equals_geometry_factory_0ulp` `:137`:
     these "closure == geometry-factory" legs become VACUOUS once the factory
     τ is gone. Retire them; their structural-independence is preserved by the
     surviving `test_sphere_tau_matches_independent_reference` `:109` /
     `test_cyl_tau_clamp_is_the_only_difference_from_reference` `:175` legs
     (closure vs contamination.py — already geometry-τ-free).
2. **Retire the inline-recompute test constructors** (`test_diamond.py:561/648/964`,
   `test_cell_balance_for_streaming.py:71/87` + the `c_from_streaming_terms`
   reads `:119-121`): re-point at closure c/τ or delete (they pin the RETIRED
   "rebuild c from st.α/st.τ" math).
3. **Delete the geometry τ PRODUCER statements** (T1, SURGICAL):
   - sphere `:681-688` τ block + `:699 tau_mm=tau_mm` — KEEP `mu_edge[0]`/`mu_start`.
   - cyl `:808-815` τ block + `:826 tau_mm_per_level=…` — KEEP `mu_start_per_level`/`sin_theta`.
   - slab `:495 tau_mm=1.0`.
4. **Delete the bake accessor reads** (T1 accessor): `streaming_terms()` lines
   `:505/:519` (sphere), `:534/:545/:554` (cyl) for τ; `:493/:494/:517/:518/:552/:553`
   for α. KEEP the rest of each `StreamingTerms(...)` constructor.
5. **Drop the top-level fields**: `ReducedStreamingOperator.tau_mm` (`:415`),
   `.tau_mm_per_level` (`:426`) [T2]; `StreamingTerms.tau_mm` (`:295`),
   `.alpha_in` (`:283`), `.alpha_out` (`:292`) [T3/T4]. Retire the
   `test_reduced_operator.py` tau/alpha-bit-identical tests (`:115/:122/:174/:184`)
   and the `st.alpha_*`/`st.tau_mm` field-match tests (`:286-288/:303-305/:326-328/:428-431`).
6. **Update Sphinx** (archivist): `discrete_ordinates.rst` (`:205/:209/:953/:954/:1177/:1233/:1244/:1478/:8771/:9336`),
   `structured_geometry.rst` (`:313/:420`). Doc-only; not a code dependency.

**NOT in Step C**: T5 (the closure `__call__`-arg `tau_mm` / unbound-legacy
`MorelMontryAngularSweep(sn_mesh=None)` path) — separate surface, survives;
track as a follow-up if the unbound-legacy call path is to be retired.

---

## Surgical-excision hazards (summary)

- **Do NOT whole-function-delete** `spherical_streaming` / `cylindrical_streaming`
  / `streaming_terms()` — they emit live `alpha_half`/`redist_dAw`/`mu_start`/
  `face_areas`/`delta_A` (and the live `StreamingTerms` fields chord/mu/faces/
  delta_A_over_w/mu_start/volume/abs_mu). Excise only the τ/α statements.
- **`mu_edge` (sphere) and `sin_theta`/`mu_start_per_level` (cyl) are SHARED**
  between the τ block and the still-live `mu_start` outputs. Surgically keep the
  shared scalars.
- **α_in/α_out fields delete cleanly** (T4 confirmed orphaned) but their
  producer arrays `alpha_half`/`alpha_per_level` STAY (live Carlson seed + the
  closure's own c-rebuild reads `alpha`).

## Regression FLOOR — does it depend on geometry-τ surviving?

- `test_tau_producer_equivalence.py`: **MIXED.** The two `*_equals_geometry_factory_0ulp`
  legs (`:77`, `:137`) DO depend on geometry-τ (delete them, step 1). The
  independent-reference legs (`:109`, `:175`) + the Cartesian-raise (`:231`)
  do NOT — they read the closure producer + contamination.py oracle, survive.
- `test_cell_visit_c_stamp.py` τ arm: **DEPENDS** on `st.tau_mm` as the L11
  oracle (`:94-96`). Re-point at `morel_montry_weights` (step 1) BEFORE deleting
  the bake. Until re-pointed, deleting `StreamingTerms.tau_mm` REDS this test.

Both files pass clean today (12 passed, 0.17s) — they are the floor to keep
green via the migrate-then-delete ordering above.
