# B0.3 — coverage repairs (working notes)

Branch `refactor/operator-strategy-layers`. Host `.venv` (Py 3.14), canonical
`.venv/bin/python -O -m pytest`, SERIAL.

Evidence key: **[M]** measured this session · **[R]** read (quoted + `file:line`)
· **[G]** grep/exhaustive.

---

## Baseline, before any edit [M]

```
.venv/bin/python -O -m pytest tests/geometry/ tests/sn/operators/ tests/diffusion/ -q -p no:randomly
→ 1224 passed, 4 skipped, 2 xfailed, 32 warnings in 15.68s
```

Skip reasons (`-rs`) [M]:

* `[3] tests/sn/operators/test_bc_extraction_matvec.py:445` —
  `2-D mesh construction not available here: tuple index out of range`
  (**the three inert sentinels — item 2**; all three are the `seed=0,1,2`
  parametrisation of ONE test, exactly as the audit §9.4 records).
* `[1] tests/sn/operators/test_solver_components.py:662` — 421-group HDF5 data
  absent. Out of scope (an environment precondition, a legitimate skip).

---

## Instrument — the auditor's own 31-mutation harness, REUSED

`$CLAUDE_JOB_DIR/tmp/` still holds the boundary review's harness, so the
before/after is apples-to-apples rather than a re-implementation:

* `mutplugin.py`  — 12 leaf-action mutations (`ORPHEUS_MUT`)
* `mutplugin2.py` — 14 guard-disabling mutations (`ORPHEUS_GUARD`)
* `mutplugin3.py` — 5 adjoint mutations (`ORPHEUS_ADJ`)
* `mutplugin4.py` — the ERR-052 `drop_renorm` mutation (`ORPHEUS_ERR052`)
* `run_mutations.sh` / `run_guards.sh` / `run_adj.sh` — the runners

Every mutation is an in-process monkeypatch with a bite check; **no file under
`orpheus/` or `tests/` is written by the harness**, and no `git checkout` is
used anywhere in this job.

---

## Repair log

### Item 3 — `tests/geometry/test_bc_errors.py`: 9 tautological legs → production legs

`tests/geometry/test_bc_errors.py` (whole file rewritten).

**Before [M]:** every leg was `with pytest.raises(X): raise err` where `err` was
built two lines above AS an `X`. **0 of 14** guard-disabling mutations reddened
the file.

**After [M]** — same 14 mutations, same harness, file run in isolation:

| guard neutered (`ORPHEUS_GUARD=`) | before | after | test that now reds |
|---|:--:|:--:|---|
| `involutive` | green | **RED** | `test_reflection_not_involutive_error` |
| `measure_preserving` | green | **RED** | `test_boundary_geometry_map_not_measure_preserving_error` |
| `inflow_to_outflow` | green | **RED** | `test_reflection_did_not_map_inflow_to_outflow_error` |
| `white_submarkov` | green | **RED** | `test_submarkov_violation_error` |
| `white_positive` | green | **RED** | `test_boundary_response_not_positive_error` |
| `albedo_submarkov` | green | **RED** | `test_submarkov_violation_error` |
| `albedo_positive` | green | **RED** | `test_boundary_response_not_positive_error` |
| `source_on_trace` | green | **RED** | `test_boundary_source_not_on_incoming_trace_error` |
| `sn_zeroflux_refusal` | green | **RED** | `test_boundary_error_base_class_contract` |
| `sn_vacuum_orientation` | green | **RED** | `test_vacuum_applied_to_outgoing_trace_error` |
| `diff_periodic_refusal` | green | **RED** | `test_boundary_error_base_class_contract` |
| `diff_prescribed_refusal` | green | **RED** | `test_boundary_error_base_class_contract` |
| `walk_drop_scalar` | green | green | *(walker, not an error guard — `test_law_composition.py` catches it, 5 reds)* |
| `walk_sum_first_only` | green | green | *(idem, 4 reds)* |

**12 / 14 — up from 0 / 14.** The two remaining are composition-walker
mutations, which are not error-raise guards at all; pointing an *error-type*
file at them would be the same category error the repair removes.

Unmutated control: `11 passed` [M]. Every red above is a `DID NOT RAISE`
(or a redirected raise for `measure_preserving`, where neutering the measure
check lets the mispaired GL-8 table fall through to the ERR-045 check —
still a red, still attributable).

**Two findings surfaced by the repair** (production; NOT fixed here — I may
not edit `orpheus/`):

1. **[G] `IncomingOutgoingTraceClassificationError` (ERR-040) has ZERO
   production raise sites.** Its documented trigger
   `BoundaryTraceLaw.assert_inflow_outflow_classification`
   (`orpheus/geometry/boundary/_base.py:251`) is a no-op ABC default that no
   concrete law overrides. The error is exported, documented in the
   `__init__.py` catalog (`:260`) and type-tested — and unreachable. A fifth
   instance of the review's declared-capability pattern, this one inside the
   error module. Its test is now type-contract-only with the absence recorded
   in the docstring; it deliberately carries NO fake production leg.
2. **[M] `law=` tag convention drift, exactly one site.** Every guard tags
   `law=` with the lowercase registry key — measured verbatim for `"vacuum"`,
   `"reflective"`, `"white"`, `"albedo"`, `"zero_flux"`, `"periodic"` — except
   `assert_source_lives_on_incoming_trace`, which tags
   `law=type(self).__name__` (`_base.py:342`) and so emits
   `"PrescribedInflow"`. The test asserts case/underscore-insensitively so it
   pins "the tag identifies the law" (dropping `law=` still reds it) without
   calcifying the drift.

### Item 4 — the two blind "hand-computed" white tests

`tests/geometry/test_boundary.py:224` (`test_white_bc_4_point_quadrature_hand_computed`)
and `:242` (`test_white_bc_axis_z_on_product_quadrature`), plus the new
textbook-table constants block above them.

**Root cause (Mode 12).** Both fed a **hemisphere-CONSTANT** ψ_out (1.0/7.0 and
2.0/9.0) and asserted the output equalled that constant. A normalised weighted
average of a constant IS that constant *for any weights*, so the measured
functional's invariance group contained the ENTIRE `w·|Ω·n̂|` formula. Neither
contained a hand-computed number: the expected `1.0` fell out of normalisation
alone, so the docstring "explicit hand calculation" over-claimed.

**Repair.** Break the invariance with a non-constant outgoing field, and make
the reference genuinely independent — the published Gauss-Legendre tables
(Abramowitz & Stegun 25.4) typed in as literals, never read off `quad.weights`
(the existing catcher `test_white_bc_returns_cosine_weighted_average`
re-transcribes the production formula — vv L11 procedural, not structural,
independence). Each gate asserts the quadrature MATCHES the published table
first, so a `Quadrature` change reds loudly instead of silently rebasing.

* GL-4: ψ = 1 on the inner outgoing ordinate, 4 on the outer, 7 incoming
  ⟹ `P = (w₁μ₁·1 + w₂μ₂·4)/(w₁μ₁ + w₂μ₂) = 2.723973656470134`.
* product(8,4), z-axis: ψ = μ_z outgoing, 9 incoming
  ⟹ `Σw_k μ_k² / Σw_k μ_k = (1/3)/0.50576403… = 0.659068878836894`, gated at
  `rtol=1e-13` — PLUS the **closed-form continuum anchor the audit found the
  suite lacked** ([G] §6 caveat 1):
  `(∫_{2π}|Ω·n̂|μ dΩ)/(∫_{2π}|Ω·n̂| dΩ) = (2π/3)/π = 2/3`, at `rtol=2e-2`
  (the numerator is exact under GL-8; the denominator carries a 1.14 %
  half-range kink error — an explained, honest gap, still 24 % away from the
  dropped-cosine value 0.5058).

**Mutation evidence [M]** — the harness's three white mutations, run against
the three white value tests in isolation:

| mutation | before | after |
|---|---|---|
| `white_nocos` (drop `\|Ω·n̂\|`) | 1 failed / 2 passed (**both blind**) | **3 failed** |
| `white_badnorm` (norm = Σw) | *(blind by the same mechanism)* | **3 failed** |
| `white_fullsphere` (no hemisphere mask) | *(idem)* | **3 failed** |

Unmutated control: `6 passed` on `-k white`.

---

### Item 5 — `test_bc_equivalence_snapshot.py`: honest label + hard-fail on a missing snapshot

`tests/geometry/test_bc_equivalence_snapshot.py:94` (marker) and `:104`
(`_load_or_skip` → `_load_snapshot`).

**Relabel.** `pytestmark` was `[l1, regression]` while the file's own header
records that the cross-implementation half was deleted at #186 and *"the
snapshots themselves now record realiser-path outputs"* — a self-generated
baseline wearing an L1 label. Now `[foundation, regression]`: `foundation` is
the honest level for "software invariant, no theory-page `:label:`", and
`regression` is what the project's own marker table calls this shape ("NOT a
verification reference"). **The file keeps its value** — it is still the widest
mutation net in the subsystem (9 of 12 leaf mutations); only its self-claim
changed.

**Hard fail.** `_load_or_skip` called `pytest.skip` on a missing baseline, on a
decoupled-roll-out rationale that expired years ago (all 8 `.npz` are tracked).
A deleted or renamed snapshot silently disabled the widest net.

**Mutation evidence [M]** — repoint `_SNAPSHOT_DIR` at an empty directory
(in-process, bite-checked):

| | pre-repair | post-repair |
|---|---|---|
| all 8 baselines missing | `8 skipped` (invisible) | **`8 errors`**, each naming the file |

Unmutated control: `8 passed`.

---

### Item 6 — four conflicting V&V level markers

`tests/diffusion/test_boundary_realizer.py` — module-level
`pytestmark = [foundation]` stacked with `@pytest.mark.l0` on
`TestRulingThreeSemantics`, so its four tests resolved with
`['foundation','l0']`.

**Measured before [M]** — under `-W error::pytest.PytestUnknownMarkWarning`:

```
INTERNALERROR> pytest.PytestUnknownMarkWarning:
  …::TestRulingThreeSemantics::test_vacuum_means_zero_incoming_current
  has conflicting V&V level markers ['foundation', 'l0']; using 'l0'
```

**Repair.** Drop the module blanket; each class states its own level —
`foundation` for the registration / albedo-table / refusal / composition /
method-space classes, `l0` for the ruling-3 semantics class (whose four tests
ARE term-level physics claims against the P1 partial-current dictionary).

**Measured after [M]:** `31 passed` under the same `-W error` invocation.
Partition unchanged and complementary: `-m l0` collects **4/31**,
`-m foundation` collects **27/31** — no test lost or changed its level.

---

### Item 7 — production-registry pollution from the trace-law tests

`tests/geometry/test_boundary_trace_law.py:53` (`_StubLaw`) + the new
import-time eviction and module-scoped autouse fixture.

`__init_subclass__(key=…)` fires at class-creation time, so merely COLLECTING
the file inserted `_stub_for_test` into the production
`BoundaryTraceLaw.registry` for the rest of the session, with no teardown.

**Measured before [M]:**

```
before import: ['albedo','periodic','prescribed_inflow','reflective','vacuum','white','zero_flux']
after  import: [ + '_stub_for_test' ]
```

**Repair.** Evict at import; re-install for exactly this module's tests via an
autouse module-scoped fixture with a `finally` teardown (so an error inside a
test cannot leave the singleton polluted). The two tests that legitimately
require the key still see it.

**Mutation evidence [M]** — an ordering probe appended after the module:

| module under test | probe `'_stub_for_test' not in registry` |
|---|---|
| pre-repair file (`git show HEAD:…`) | **FAILED** |
| post-repair file | **PASSED** (`15 passed`) |

Import-time leak after repair: `set()`.

---

### Item 2 — the three permanently-inert sentinels

`tests/sn/operators/test_bc_extraction_matvec.py:406`
(`test_vacuum_2d_cartesian_bulk_bit_identical`, seeds 0/1/2) — the three skips
were the three parametrisations of this ONE row, exactly as the audit records.

**Root cause.** The body built a 1-D `Mesh1D` and then read
`sn_mesh.spatial_shape[1]`, raising `IndexError` on a 1-tuple; the enclosing
`except Exception as exc: pytest.skip(f"…{exc}")` swallowed it. The
self-described SENTINEL had **never executed an assertion**.

**It is a test-side fixture issue, not a production one** — checked as the
brief asks: `Mesh2D` exists and `SNMesh` accepts it (`tests/sn/_test_helpers.py:284`
already builds one). **No `orpheus/` change was needed.**

**Repair.** (a) build a genuine `Mesh2D`, deliberately NON-SQUARE (nx=3, ny=4)
so an x↔y transposition cannot hide behind a square-box symmetry (vv §H2);
(b) delete the blanket `except` so a construction failure is a FAILURE;
(c) keep the `ny > 1` check as a hard `assert`, never a skip;
(d) capture the three baselines.

**Honest-scope re-pose.** The docstring called this "the SENTINEL that the 1-D
O.4a.2 carve did not perturb the 2-D path". That claim is **not recoverable** —
O.4a.2 landed long ago and a baseline captured today cannot testify about it.
The row is now labelled what it is: a **regression floor captured 2026-07-30**,
whose correctness anchor is NOT the snapshot but the structurally-independent
2-D gates that already exist and were re-run green before capture
(`tests/sn/sweep/cartesian_2d/test_2d_l2_matvec_correctness.py`,
`test_scan_march_equivalence.py`, `test_2d_full_field_oracle.py` → `23 passed`).
Captured arrays sanity-checked: shape `(24,1,3,4)`, all finite, 288/288
non-zero, `|a|max ≈ 2.2–2.5e+01`, seeds mutually distinct.

**Mutation evidence [M]** — two mutations, each bite-checked, run against the
pre-repair file (`git show HEAD:…`) and the repaired one:

| mutation | pre-repair | post-repair |
|---|---|---|
| `drop_sy` — zero the TRANSVERSE (y) streaming coefficient in `ScanMarch._loss_action_interior` (2-D-only: the 1-D arm returns earlier) | **3 SKIPPED** | **3 FAILED** |
| `degenerate` — `SNMesh.__init__` raises `IndexError` on a `Mesh2D` (the pre-repair failure verbatim) | **3 SKIPPED** ("*2-D mesh construction not available here: tuple index out of range*") | **3 FAILED** |

File total: `30 passed, 3 skipped` → **`33 passed, 0 skipped`**.

---

### Item 1 — `catches("ERR-052")` re-posed onto the mechanism

`tests/sn/operators/test_boundary_conditions.py` — the marker moved off
`test_vacuum_keff_lower_than_reflective` onto the new
`test_power_iteration_renormalises_to_unit_production_rate`, plus two shared
helpers (`_err052_fixture`, `_hand_production_rate`).

**Phantom re-confirmed on THIS tree [M]** (the audit measured it at
`73627b71`; the tree has since moved): with ERR-052 re-introduced via the
auditor's own `mutplugin4` (`ORPHEUS_ERR052=drop_renorm` — the
`ProductionRateSolver` narrowing made never to match, i.e. the pre-fix
un-normalised trajectory verbatim), the pre-repair file reports **`11 passed`**.

**The mechanism.** The production fix divides the iterate by its production
rate each outer step, whose documented consequence is the output convention
`P(φ) = ∫_V Σ_g νΣ_{f,g} φ_g dV = 1`. That is the instrumented quantity the
bug moves, so it is what the new test asserts — at `rtol=1e-12`, **no margin**,
and (the point) **no regime dependence**: `P(φ)=1` holds at every outer count,
so this marker cannot go blind the way its predecessor did when the fixture
drifted to 6 outers. The reference is hand-assembled from
`result.scalar_flux.values`, the mesh cell widths and the mixture's `SigP` —
never via `SNSolver.compute_production_rate`, which is the routine the
normaliser divides by (vv L11; a reference computed by the normaliser would be
a tautology). The `(n,2n)`-free-ness of the fixture is asserted, so the hand
formula cannot silently go wrong if the reference case gains a Σ₂ channel.

**Mutation evidence [M]** — measured production rate:

| leg | unmutated | ERR-052 re-introduced |
|---|---|---|
| reflective (k = 1.875) | 0.9999999999999998 | **0.84375** (16 % low) |
| vacuum (k = 0.164) | 1.0 | **0.0739524** (**13.5× low**) |

| file | control | ERR-052 re-introduced |
|---|---|---|
| pre-repair | 11 passed | **11 passed** ← the phantom |
| post-repair | **13 passed** | **2 failed**, 11 passed |

Both tolerance legs (`1e-7` → 6 outers, `1e-12` → 32 outers) red.

**FINDING recorded in the test's docstring — the catalogue's failure
mechanism is NOT reachable in this fixture.** ERR-052 documents
"denormalised-FP underflow after ~30-60 outers". Measured on the vacuum leg
with the bug re-introduced, `|φ|max` sits at **5.6552e-01 after 6, 24 AND 32
outers** — the un-normalised iterate reaches a *stable* magnitude rather than
decaying without bound, because once `k` has converged the `F·φ/k` step is
scale-neutral. So no tolerance tightening could ever have rescued the old
ordering assertion; only moving to the mechanism does. (The `tol`
parametrisation is retained as the regime-independence proof, not as a route
to the underflow.) The ERR-052 catalogue entry's "Test reference" line still
names the OLD test — it needs updating to the new one.

---

## The 31-mutation sweep, re-run — **31 / 31 CAUGHT** (was 30 / 31)

Same instrument, same mutation names, same file sets (`sweep31.sh` drives
`mutplugin` / `mutplugin2`; the 5 adjoint rows run `mutplugin3` on the named
catcher files). Every mutation carries its original bite check.

### 12 leaf-action mutations (`ORPHEUS_MUT`)

| mutation | audit reds | now | verdict |
|---|---:|---:|---|
| `refl_identity` (M1) | 13 | 13 | CAUGHT |
| `refl_shift` (M2) | 12 | 12 | CAUGHT |
| `white_nocos` (M3) | 11 | **13** | CAUGHT (+2 — the repaired white gates) |
| `white_badnorm` (M4) | 16 | 16 | CAUGHT |
| `white_fullsphere` (M5) | 16 | 16 | CAUGHT |
| `vacuum_complement` (M6) | 8 | 8 | CAUGHT |
| `psrc_nomask` (M7) | 3 | 3 | CAUGHT |
| `reflect_nosel` (M8) | 5 | 5 | CAUGHT |
| `albedo_1minus` (M9) | 7 | 7 | CAUGHT |
| `albedo_1minus_diff` (M10) | 3 | 3 | CAUGHT |
| `diff_zeroflux_sign` (M11) | 4 | 4 | CAUGHT |
| `diff_vacuum_one` (M12) | 3 | 3 | CAUGHT |

Only M3 moved, by exactly the two white tests item 4 repaired. M4/M5 were
already caught by the old constant-input form (the algebra: a hemisphere-
constant ψ cancels the cosine factor ONLY when the normaliser stays
consistent, which is M3 alone) — so the +2 lands precisely where the
blindness was.

### 14 guard-disabling mutations (`ORPHEUS_GUARD`)

| mutation | audit reds | now | verdict |
|---|---:|---:|---|
| `involutive` | 1 | **2** | CAUGHT |
| `measure_preserving` | 2 | **3** | CAUGHT |
| `inflow_to_outflow` | 2 | **3** | CAUGHT |
| `white_submarkov` | 1 | **2** | CAUGHT |
| `white_positive` | 1 | **2** | CAUGHT |
| `albedo_submarkov` | 1 | **2** | CAUGHT |
| `albedo_positive` | 1 | **2** | CAUGHT |
| `source_on_trace` | 2 | **3** | CAUGHT |
| `sn_zeroflux_refusal` | 1 | **2** | CAUGHT |
| `sn_vacuum_orientation` | 2 | **3** | CAUGHT |
| `diff_periodic_refusal` | 2 | **3** | CAUGHT |
| `diff_prescribed_refusal` | 1 | **2** | CAUGHT |
| `walk_drop_scalar` | 7 | 7 | CAUGHT |
| `walk_sum_first_only` | 5 | 5 | CAUGHT |

**Exactly +1 on each of the twelve error guards** — the `test_bc_errors.py`
repair, landing one-for-one where it was designed to. The two walker rows are
unchanged, as intended.

### 5 adjoint mutations (`ORPHEUS_ADJ`)

| mutation | reds | verdict |
|---|---:|---|
| `adjointable_any` (A1) | 1 | CAUGHT |
| `white_advertises_transpose` (A2) | 1 | CAUGHT |
| `transpose_nomask` (A3) | 12 | CAUGHT |
| `corner_transpose_swap` (A4) | 3 | CAUGHT |
| `ruled_corner_widened` (A5) | 1 | CAUGHT |

(A3's 12 vs the audit's 9 is a different file set, not a coverage change — I
ran the five named catcher files rather than `run_adj.sh`'s twenty, because the
full set costs ~2 min per mutation. The verdict is what the gate asks for.)

### Plus the 32nd — the one the audit recorded as MISSED

| mutation | audit | now |
|---|---|---|
| `ORPHEUS_ERR052=drop_renorm` | **MISSED** (11 passed) | **CAUGHT** (2 failed) |

**Catch rate: 31/31 on the auditor's set (≥ 30/31 required), and the single
historical miss is closed.**

---

## Suite results

| run | before | after |
|---|---|---|
| `tests/geometry/ tests/sn/operators/ tests/diffusion/` (the required gate) | 1224 passed, **4 skipped**, 2 xfailed | **1229 passed, 1 skipped, 2 xfailed** (20.2 s) |
| the 18-file boundary harness | 303 passed, **3 skipped** | **328 passed, 0 skipped** |
| `tests/test_vv_harness_audit.py tests/numerics/ tests/transport/fields/ tests/sn/primitives/` | — | 1392 passed |

The one remaining skip in the gate is the legitimate environment precondition
(`test_solver_components.py:662`, 421-group HDF5 data absent) — the only skip
in the boundary surface that names a *precondition* rather than an exception.

---

## Docs / V&V-matrix impact

**[X] RETRACTION — an earlier draft of these notes claimed the committed
matrix had "always" recorded the snapshot file as `foundation` while the code
said `l1`, i.e. a pre-existing code↔doc inconsistency. That was WRONG, and it
was my own measurement error: I diffed against a working-tree copy that had
already been regenerated mid-session. There is no such inconsistency.**

The truth, read from git (`git show HEAD:docs/theory/verification/matrix.rst`):

```
Total tests collected: **7051**        L1, 1192, 16.9%   foundation, 4592, 65.1%
diffusion/test_boundary_realizer,       4, 0, 0, 0, 27, 0
geometry/test_bc_equivalence_snapshot,  0, 8, 0, 0,  0, 0     <- L1, consistent with the pre-repair marker
operators/test_boundary_conditions,     0, 0, 0, 0, 11, 0
```

The committed matrix agrees with the pre-repair code exactly. What my repairs
change, and what the next regeneration must pick up:

| row | HEAD | after B0.3 | cause |
|---|---|---|---|
| `geometry/test_bc_equivalence_snapshot` | `0, 8, 0, 0, 0, 0` | `0, 0, 0, 0, 8, 0` | item 5 relabel |
| `operators/test_boundary_conditions` | `0,0,0,0,11,0` | `0,0,0,0,13,0` | item 1 (+2 legs) |
| `Total tests collected` | 7051 | 7053 | item 1 |
| `L1` / `foundation` | 1192 / 4592 | 1184 / 4602 | items 5 + 1 |
| `` `ERR-052` `` coverage count | 1 | 2 | item 1 |
| `diffusion/test_boundary_realizer` | `4,0,0,0,27,0` | **unchanged** | item 6 is level-NEUTRAL |

Item 6 moves no counts at all — `conftest.py`'s documented tiebreak already
resolved `['foundation','l0']` to L0, so the repair removes only the
`PytestUnknownMarkWarning`. Verified by running the audit with a verbatim
pre-repair copy alongside the repaired file in ONE collection:

```
diffusion/test_zzprerepair_diffbr_TMP    4    0    0    0   27    0
diffusion/test_boundary_realizer         4    0    0    0   27    0
geometry/test_zzprerepair_snap_TMP       0    8    0    0    0    0   <- L1
geometry/test_bc_equivalence_snapshot    0    0    0    0    8    0   <- foundation
                                        L0   L1   L2   L3   FD   ??
```

**Working-tree state, for the main agent.** `docs/theory/verification/matrix.rst`
is currently MODIFIED in the working tree — regenerated at 22:08:31 by a
concurrent Sphinx build (`docs/conf.py:125`'s `builder-inited` hook), between
my item-5 edit (22:07:04) and my item-1 edit (22:18:23). It therefore carries
the item-5 relabel but NOT the item-1 `+2` rows. **I did not write it and have
not touched `docs/`.** The next Sphinx build regenerates it completely; no
manual action is needed beyond letting that happen before the commit.

## Findings for the plan (production; not fixed here — `orpheus/` is off-limits)

1. **[G] `IncomingOutgoingTraceClassificationError` (ERR-040) is unreachable** —
   zero `raise` sites in `orpheus/`; its documented trigger
   `assert_inflow_outflow_classification` (`_base.py:251`) is a no-op ABC
   default nobody overrides. A fifth declared-capability instance, inside the
   error module. Fold into the B0.1 "retire the false promises" list, or give
   it a raiser.
2. **[M] `law=` tag drift at exactly one site** — `_base.py:342` tags
   `type(self).__name__` (`"PrescribedInflow"`) where every other guard tags
   the lowercase registry key. One-line fix; the test is written so it will not
   need changing when it lands.
3. **[M] ERR-052's documented failure mechanism is unreachable in its fixture**
   — the un-normalised iterate stabilises at `|φ|max = 5.6552e-01` (6, 24 and
   32 outers alike) instead of decaying to denormal, because `F·φ/k` is
   scale-neutral once `k` converges. The catalogue's "~30-60 outers →
   denormalised FP" is not what this configuration does. Worth a note on the
   ERR-052 entry.
4. **The ERR-052 catalogue "Test reference" line is now stale** — it names
   `test_vacuum_keff_lower_than_reflective`; the marker has moved to
   `test_power_iteration_renormalises_to_unit_production_rate`.
5. **The 2-D vacuum sentinel's original claim is unrecoverable** (item 2) — its
   three baselines are a floor captured today, not pre-carve evidence. If the
   plan wants genuine "the 1-D carve did not perturb 2-D" evidence, that has to
   come from the structurally-independent 2-D gates, not from this snapshot.
6. **`docs/theory/verification/matrix.rst` needs one regeneration** for the
   item-1 `+2` rows (the item-5 half is already in the working tree from a
   concurrent Sphinx build). The next `sphinx-build` does it automatically.
7. **`_load_or_skip` had no external consumers** (three-search retirement audit:
   graph, text-grep across code/tests/docs, direct callers) — the rename to
   `_load_snapshot` is contained. Same for `_stub_for_test`: **zero** consumers
   outside `test_boundary_trace_law.py`, so the registry fixture is safe.
