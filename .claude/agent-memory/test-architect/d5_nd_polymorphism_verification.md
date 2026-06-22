---
name: d5-nd-polymorphism-verification
description: Verification spec for #240 Phase 2 Step D5 — completing spatial-closure polymorphism (DD/Step/LD × CumprodScan/ScanMarch/MovingFrontierWindow/FullFieldWavefront/_OneDimScanWalk × d=1/2/3). Three threads — D5-0 routing-honesty (the CONFIRMED-LIVE 2-D LD→ScanMarch misroute), D5a N-D DD scan-march fold onto the coefficient model (#239, principled re-baseline), D5b N-dim LD on the DAG wavefront (#158 Inc D / #38, the NEW bilinear math). Per-cell coverage matrix = definition-of-done. Multi-D LD MMS (μ-non-trivial, Q̂≠0 slope source). Orphan-label list for D6.
metadata:
  type: project
---

# #240 Phase 2 Step D5 — N-dim complete-polymorphism verification spec

**Status:** PRE-IMPLEMENTATION. Branch `feature/sn-space-angle-tier2`,
HEAD `0cc0cbf`. Host env, canonical `python -O -m pytest`. This spec
SHAPES the carve (the proactive L17 cross-strategy/cross-dim dispatch).

## ⭐ LIVE STATE CORRECTIONS (probed this session — the brief reflects a
## slightly earlier snapshot; these supersede)

1. **The `s_axes` convention change (Phase-2 Step-A / #239 coupled) HAS
   ALREADY LANDED.** Both `diamond.py:cell_kernel_batch/residual_kernel_batch`
   (`:348-352`, `:390-395`) AND the LD kernel (`linear_discontinuous.py:437`
   `g = s_axes[0]  # was 0.5*s`) consume the RAW down-face streaming `g`,
   with the diamond `2` internalised in the DD kernel. The ScanMarch inline
   DD at `loss_representation.py:1320-1321` ALSO already applies `sx2 =
   2.0*s_x` / `sy2 = 2.0*s_y` with the `NOTE(#239)` marker. **So D5a is NOT
   "move the factor-2" (done) — D5a is "fold the inline DD onto
   `affine_scan_coefficients` so LD can ride the 2-D row-march and Step
   shares the path".** The convention groundwork is in; D5a is the
   coefficient-model lift on top of it.
2. **LD Increment B has fully landed.** `LinearDiscontinuous.is_affine_scannable
   = True` (`:278`); LD has its 1-D `affine_scan_coefficients` (`:508`); a 1-D
   slab LD mesh routes to `CumprodScan`; the two-paths gate
   `test_ld_two_paths_scan_equals_dag_oracle` (`test_mms_ld_slab.py:151`)
   exists and passes. So the D5b "LD on the DAG wavefront" is specifically
   the **multi-D (d≥2) bilinear** LD kernel — the 1-D DAG path already works.
3. **`Step` does NOT exist in production.** It appears ONLY as a docstring
   example (`scheme.py:538`, `class Step(...key="step")`) and in the
   advection-scheme narrative (`scheme.py:566` "Step ↔ first-order upwind",
   `is_positivity_preserving=True`). There is NO registered `Step` scheme,
   no `Step` tests. **The matrix's `Step` rows are "NOT-YET (scheme
   unimplemented)" everywhere** — Step is the SLOPELESS sibling of DD
   (affine-scannable, `w` weight differs: upwind `w=1` vs DD `w=½`), so when
   Step lands it rides EXACTLY the DD path set (D5a's fold is generic over
   DD+Step by construction). Flag: if the campaign wants Step verified, it
   needs its own implementation task — NOT in D5's scope, but the D5a fold
   makes Step polymorphic-READY.
4. **The D5-0 misroute is CONFIRMED LIVE (probe):**
   `ScanMarch.supports(2-D Cartesian LD mesh).ok == True` →
   `default_for(2-D LD) == ScanMarch`. ScanMarch's `_sweep_interior`
   (`:1312-1347`) runs inline DD (`alpha = 2·sx2/D_row − 1`) with NO scheme
   dispatch → a 2-D LD mesh **silently computes DD, dropping LD's bilinear
   slope**. The `FullFieldWavefront` 2-D LD path RAISES
   (`cell_kernel_batch(d=2)` → `NotImplementedError "d=1 only"`), so today a
   2-D LD mesh EITHER misroutes silently (ScanMarch, the default) OR raises
   (if the wavefront is forced). D5-0 makes the silent path raise-or-route
   honestly.

## Baselines (live @ HEAD `0cc0cbf`, `-O`)

```
tests/sn/spatial/test_linear_discontinuous.py tests/sn/verification/mms/test_mms_ld_slab.py
  → 24 passed, 1 xfailed   (the 1 xfail = test_ld_thick_diffusive_limit_xfail)

tests/sn/sweep/core tests/sn/solve -W "error::tests.sn.regression._regression_assert.DriftWarning"
  → 505 passed, 1 skipped, 4 xfailed   (the DD strict bit-id gate = D5's negative control)
```

Route-around the 7 pre-existing reds; NEVER run all of `tests/sn` (#212 hang):
```
-k "not (vacuum_bulk_bit_identical_1d and SPH) and not (sphere_1g_apply_bit_identical \
    or sphere_2g_apply_bit_identical) and not test_2d_mesh_resolution \
    and not two_d_cartesian_loss_action"
```
Target dirs: `tests/sn/operators tests/sn/spatial tests/sn/sweep/core
tests/sn/sweep/cartesian_2d tests/sn/solve`.

---

# 1. THE PER-(scheme × strategy × dim) COVERAGE MATRIX — definition-of-done

Cell legend:
- **T:X** = VALID and tested by test X (file::test).
- **GAP** = VALID path, UNtested (a real coverage hole to close).
- **EXCL(R)** = STRUCTURALLY excluded via `supports()`/raise with reason R
  (a mathematical impossibility — the ONLY permitted exclusion).
- **NOT-YET** = capability not implemented (the campaign's open work).
- **n/a** = the strategy is dimension-restricted by spine design (e.g.
  `CumprodScan`/`_OneDimScanWalk` are 1-D-only by construction).

Strategy abbreviations: CS=`CumprodScan`, SM=`ScanMarch`, MFW=`MovingFrontierWindow`,
FFW=`FullFieldWavefront`, ODW=`_OneDimScanWalk` (the internal 1-D scan engine CS/SM
delegate to in d=1).

## DD (Diamond Difference) — affine-scannable, w=½

| Strategy | d=1 | d=2 | d=3 |
|----------|-----|-----|-----|
| CS  | T: `test_mms_ld_slab::routes_to_cumprod_scan`(DD half) + `test_sweep_regression`(slab DD scan) | n/a (CS 1-D-only) | n/a |
| SM  | T: dispatch(DD slab→CS, but SM admits 1-D); covered by `test_scan_march_equivalence` | **D5a TARGET → T: `test_scanmarch_sweep_equals_oracle`** (already DD-verified; D5a keeps it green through the fold) | EXCL("row-march kernels unpack d=2; d≥3 → FFW spine", `ScanMarch.supports`, pinned `TestD3SupportsMatrix::test_scan_march_refuses_d3`) |
| MFW | n/a (MFW Cartesian d≥2) | T: `test_2d_full_field_oracle`(window≡full) | GAP→NOT-YET (#227 d≥3 kernels; FFW serves d≥3) |
| FFW | T: `test_mms_ld_slab`(DD via FFW slab oracle) | T: `test_2d_full_field_oracle` + `test_mms_2d` | T: `TestD3SupportsMatrix::test_d3_cartesian_falls_through_to_full_field_spine` + k_inf-3D (#225) |
| ODW | T: `test_sweep_regression`(1-D DD scan) + slab apply snapshots | n/a | n/a |

## Step (first-order upwind) — affine-scannable, w=1 — **scheme UNIMPLEMENTED**

| Strategy | d=1 | d=2 | d=3 |
|----------|-----|-----|-----|
| CS / SM / MFW / FFW / ODW | **NOT-YET (Step scheme not built — docstring stub only `scheme.py:538`)** | NOT-YET | NOT-YET |

Step is slopeless → rides the SAME path set as DD once implemented. D5a's
coefficient-model fold makes Step polymorphic-READY (generic over DD+Step).
Building Step is a SEPARATE task; D5 only guarantees the path EXISTS for it.
A `xfail(reason="Step scheme #NNN")` placeholder is OPTIONAL (see §5).

## LD (Linear-Discontinuous) — affine-scannable (1-D), bilinear (multi-D)

| Strategy | d=1 | d=2 | d=3 |
|----------|-----|-----|-----|
| CS  | T: `test_mms_ld_slab::routes_to_cumprod_scan` + `test_ld_two_paths_scan_equals_dag_oracle` + `test_sn_1d_slab_ld_mms_converges_second_order` | n/a (CS 1-D-only) | n/a |
| SM  | T: 1-D LD rides SM via ODW (SM.sweep 1-D → ODW) — covered by Krylov≡SI | **EXCL("bilinear slope coupling needs the wavefront", `ScanMarch.supports` NEW predicate — D5-0+D5b) — the structural exclusion** | EXCL(same + d≥3 row-march deferred) |
| MFW | n/a | **NOT-YET → D5b TARGET (multi-D bilinear LD kernel)** | NOT-YET (#227) |
| FFW | T: `test_sn_1d_slab_ld_mms_converges_second_order`(LD via FFW oracle) | **NOT-YET → D5b TARGET (the multi-D LD MMS + two-paths)** | NOT-YET (d≥3 LD bilinear) |
| ODW | T: `test_ld_two_paths` + Krylov≡SI(`test_sn_1d_slab_ld_mms_krylov_matches_si`) | n/a | n/a |
| (curvilinear 1-D) | EXCL("curvilinear LD unpublished", `affine_scan_coefficients` raise on non-neutral curvature + `_require_slab` + `test_ld_curvilinear_scan_rejected`) | — | — |

**Definition-of-done after D5 (the matrix's terminal state):**
- D5-0: the LD `SM d=2/d=3` cell flips from a SILENT MISROUTE to a declared
  `EXCL("bilinear slope coupling needs the wavefront")`.
- D5b: the LD `FFW d=2` and `MFW d=2` cells flip `NOT-YET → T` (the multi-D
  LD MMS + two-paths gate). d=3 LD stays `NOT-YET` (out of D5b scope — flag).
- D5a: the DD `SM d=2` cell stays `T` (the fold is value-preserving; the
  existing oracle gate is the catcher) AND the path becomes scheme-generic
  so the LD-2-D-via-SM remains correctly EXCL (LD never reaches the inline
  DD), while Step-2-D-via-SM becomes polymorphic-ready.
- Every remaining `NOT-YET` carries a tracked issue (#227 d≥3 kernels, Step
  scheme, d=3 LD).

**The acceptance lens:** after D5, EVERY non-`NOT-YET` cell is either `T`
(tested) or `EXCL(reason)` (a `supports()`/raise with a mathematical reason).
There is ZERO `GAP` (untested-but-valid) and ZERO silent code gap. The ONLY
`NOT-YET` cells are explicitly-tracked unimplemented capabilities.

---

# 2. PER-THREAD GATE SPECS

The three threads have DIFFERENT risk profiles and bit-identity postures:

| Thread | Posture | Risk | Sequencing |
|--------|---------|------|-----------|
| **D5-0** routing honesty | bit-identical / structural | LATENT-WRONG (live misroute) | **FIRST** (close the silent bug before any math moves) |
| **D5a** N-D DD scan-march fold | principled ~1-ULP re-baseline | value-preserving (LOUD if broken) | **SECOND** (stabilise SM before D5b new math) |
| **D5b** N-dim LD bilinear | NEW capability (nothing to be identical to) | new-math (MMS-gated) | **THIRD** (the genuinely new closure) |

## ════ D5-0 — ROUTING HONESTY (the latent-misroute fix) ════

**SUT:** a dimension-aware `ScanMarch.supports` predicate. Today
(`loss_representation.py:1226-1233`):
```python
(mesh.is_1d or (mesh.is_cartesian and mesh.ndim == 2)) and mesh.scheme.is_affine_scannable
```
conflates "1-D prefix-scannable" (LD: yes) with "N-D scan-march-compatible"
(LD: NO — bilinear). The fix asks the precise question: a scheme is
scan-march-compatible **in d≥2** only if its multi-D closure is the inline-DD
row-march shape (slopeless: DD, Step) — NOT if it merely has a 1-D affine
recurrence. The cleanest predicate (recommend a scheme trait, NOT an
isinstance — coding-elegance anti-pattern #4): a NEW
`is_scan_march_compatible` ClassVar (or reuse the existing `is_affine_scannable`
for d=1 and gate d≥2 on a NEW trait). The 1-D arm keeps
`is_affine_scannable` (LD stays true there); the d≥2 arm gates on the new
trait (LD false, DD/Step true).

**⚠ Design discipline (flag to method-implementer):** do NOT special-case LD
by name. The predicate must be a SCHEME TRAIT the scheme declares (DD/Step
declare scan-march-compatible-in-N-D=True; LD declares False because its
multi-D closure is bilinear). This keeps the `supports()` honest and
extensible. The trait's MEANING: "my multi-D cell closure is the slopeless
inline-DD row-march `α/β/D_row` shape" — equivalently "I have no transverse
slope coupling".

### GATE D5-0.1 — NEGATIVE routing gate (the headline D5-0 gate)

**Claim:** selection contract (structural). A 2-D Cartesian LD mesh does
NOT select `ScanMarch`; it EITHER selects the wavefront (`MovingFrontierWindow`
or `FullFieldWavefront`) OR construction raises with a structural reason.

**File:** EXTEND `tests/sn/sweep/core/test_unified_sweep_dispatch.py::TestD3SupportsMatrix`
(confirmed live @ `:197`). **Crucial:** the existing `_fake` helper
(`:215-225`) ALWAYS supplies `scheme=SimpleNamespace(is_affine_scannable=True)`
— it CANNOT exercise the LD-multi-D case. **ADD a NEW fake variant** with
the LD-shaped trait (`is_affine_scannable=True` for the 1-D scan trait but
the NEW `is_scan_march_compatible=False`), OR build a REAL 2-D LD `SNMesh`
(the probe pattern verified live this session: `SNMesh(Mesh2D(...), q, mats,
scheme=LinearDiscontinuous())`).

```python
@pytest.mark.foundation
def test_scan_march_refuses_2d_non_slopeless_scheme(self):
    """D5-0: a 2-D Cartesian LD mesh (bilinear slope coupling) must NOT
    select ScanMarch — its row-march runs inline DD, dropping LD's slope.
    This is the CONFIRMED-LIVE misroute (default_for(2D LD) == ScanMarch
    pre-D5-0). After D5-0 the predicate refuses it."""
    sn = _build_2d_ld_mesh()                      # real SNMesh, LD scheme
    if ScanMarch.supports(sn).ok:
        pytest.fail(
            "ScanMarch.supports admitted a 2-D LD mesh — its row-march "
            "kernel runs inline DD (loss_representation.py:1320), silently "
            "dropping LD's bilinear slope. The misroute computes DD, not LD."
        )

@pytest.mark.foundation
def test_2d_ld_default_for_routes_to_wavefront_or_raises(self):
    """D5-0 selection: the LIVE default_for on a 2-D LD mesh lands on a
    WAVEFRONT (MFW/FFW) — never ScanMarch. (Before D5b lands the multi-D
    LD kernel, the wavefront RAISES on the LD multi-d cell_kernel_batch;
    after D5b it runs. EITHER outcome is honest — what is forbidden is the
    silent ScanMarch inline-DD path.)"""
    sn = _build_2d_ld_mesh()
    sel = default_for(sn)                          # MUST be MFW or FFW
    if isinstance(sel, ScanMarch):
        pytest.fail("2-D LD routed to ScanMarch (inline DD) — misroute")
```

- **Tolerance/assertion:** structural (`supports().ok` boolean / `isinstance`
  + `pytest.fail` — `-O`-safe, Mode-8-clean).
- **Config:** real 2-D Cartesian LD mesh, `level_symmetric` quad (genuine
  `mu_y`, #214-safe), NON-SQUARE (`nx=4, ny=3` — x↔y blindness defense).
  `n_groups ∈ {1,2}` not load-bearing for a SELECTION gate (no flux solved),
  but a 2G mesh exercises the `is_scan_march_compatible` read on a realistic
  material map.
- **V&V level / markers:** `@pytest.mark.foundation` (no `verifies` — pins a
  selection contract, not an equation).
- **Failure mode caught:** Mode #6 convention drift (the trait conflation) +
  the AI-adjacent "trait answers the wrong question" — a SELECTION-honesty
  bug. Catches the live D5-0 misroute.

### GATE D5-0.2 — DD/Step routing NO-REGRESSION (the predicate narrowing lost nothing)

**Claim:** the narrower `ScanMarch.supports` still admits DD-in-2-D (and
1-D LD via the 1-D scan trait). `default_for(2-D DD) == ScanMarch` (the
S6.9 Fork-B2 production default); `default_for(1-D LD) == CumprodScan`.

**File:** the existing `TestD3SupportsMatrix::test_scan_march_keeps_d2_and_1d_coverage`
(`:237`) + `test_mms_ld_slab::test_ld_slab_mesh_routes_to_cumprod_scan` (`:58`)
already pin these. **D5-0.2 is the assertion they STAY green** — re-run, do
NOT add new. **⚠ the `_fake` at `:240-244` must keep a 2-D fake that the
narrowed predicate still admits (DD-shaped: `is_scan_march_compatible=True`)
— UPDATE the fake to carry the new trait so the d=2 admission pin exercises
the trait, not just `is_affine_scannable`.**

- **Failure mode caught:** Mode #5 over-narrowing (a too-aggressive predicate
  that also refuses DD-in-2-D would silently drop the production default to
  MFW — a perf regression + a routing surprise). This is the positive half
  of the anti-pattern-#11 contract pair (must-select DD + must-NOT-select LD).

---

## ════ D5a — N-D DD SCAN-MARCH ONTO THE COEFFICIENT MODEL (#239) ════

**SUT:** rewrite ScanMarch's two inline-DD interior kernels off the hard-coded
`alpha/beta/D_row` math and onto the scheme's affine coefficients. Two blocks:
- `_sweep_interior` (`:1312-1347`): `sx2 = 2.0*s_x`, `D_row = sigt + sx2 + sy2`,
  `alpha = 2.0*sx2/D_row − 1`, `beta = 2.0*(Q + sy2*ψy)/D_row`.
- `_loss_action_interior` (`:1426-1460`): `D_row`, `LpC_oct = D_row·ψ̄ −
  sx2·in_x − sy2·ψy`.

Both carry `NOTE(#239)` markers. The fold: the x-scan recurrence consumes
the scheme's `affine_scan_coefficients` (`a, inverse_denom, w`); the
transverse `sy2·ψy_in` is a DIRECT (known-upstream) term folding into the
scan's effective source via `source_emission`; the closure uses
`cell_average`/`outgoing_face_from_average`. The DD reference kernels
(`diamond.py:293-401`, the explicit `zip(s_axes, psi_in)` left-fold) are the
d-generic structural model.

**⭐ This is a PRINCIPLED ~1-ULP RE-BASELINE of existing DD math** (precedent:
#158-B1, the s_axes Step-A re-baseline). The converged VALUES are unchanged;
only the arithmetic path (the `α/β` recurrence vs the coefficient recurrence)
re-associates. Generic over DD AND Step (both slopeless → both fold).

### GATE D5a.1 — DD 2-D ScanMarch ≡ FullFieldWavefront oracle (the headline)

**Claim:** software invariant (foundation). The genericised ScanMarch row-march
≡ the d-generic full-field DD kernel, on the SAME 2-D config.

**File:** the EXISTING `tests/sn/sweep/cartesian_2d/test_scan_march_equivalence.py::
test_scanmarch_sweep_equals_oracle` (`:104`) + `test_scanmarch_residual_equals_oracle`
(matvec leg) ALREADY pin this and are DD-verified. **D5a.1 is the assertion
they STAY green through the fold** (the rewrite must not change the value).
NON-bit-id across schedules (`assert_allclose` nulp-tolerant) — SURVIVES by
design.

- **Config (Mode-9 degeneracy-break — MANDATORY):** confirm the existing
  parametrize covers het / 2G-asymmetric-downscatter / non-flat / NON-SQUARE
  (`nx≠ny`) / `level_symmetric` (genuine `mu_y`). Read the `:104` decorator;
  if no tuple is `(nx≠ny, ng=2, bc=specular)`, ADD one. A flat-flux square
  box NULLS the streaming-redistribution and HIDES x↔y swaps — a coefficient-
  fold bug (e.g. transverse `sy` folded into the wrong scan source) is
  invisible there.
- **Tolerance:** `assert_allclose(rtol=_RTOL, atol=_ATOL)` (nulp-class — the
  two reduction trees differ by FP-association). A FOLD bug (dropped factor,
  wrong transverse coupling) is a clean factor → fires ~1e15 ULP, LOUD.
- **Markers:** `@pytest.mark.foundation` (existing).
- **Failure mode caught:** Mode #1 sign flip (transverse `−sy·ψy` sign in the
  refactored source), #3 missing factor (the `2` lost in the coefficient
  fold), #6 convention drift (`α/β` vs coefficient-model grouping).

### GATE D5a.2 — DD 2-D matvec vs FROZEN snapshot (the value-preserving linchpin)

**Claim:** the re-associated 2-D matvec produces the value FROZEN under the
pre-fold convention to nULP. **A two-paths gate (D5a.1) is BLIND to a bug that
moves BOTH paths identically** (e.g. both ScanMarch AND the FFW oracle, if
the oracle also re-routed — though here only ScanMarch changes, so D5a.1 IS a
true cross-check). Still, PAIR with a frozen-reference gate so a fold bug that
preserves the ScanMarch≡oracle relation but shifts the value is caught.

**File:** EXTEND `tests/sn/operators/test_streaming_operator.py::
TestT4bPreT4RegressionSnapshot` with a `cart2d` apply arm (the `cart2d_*_apply_*`
snapshots EXIST in `pre_t4_snapshots.npz`, frozen pre-#240, currently
UNCONSUMED — per the Step-A note §2a). bulk →
`assert_regression(kind="direct", reduction_depth=mesh.nx)` (the ONLY
admissible drift is ~1 ULP; a value shift blows the nulp bound by ~1e15 ULP).
Boundary trace → STRICT `assert_array_equal` (0 ULP — the outflow defect
reconstructs from the same `2ψ̄−in` faces).

- **NOTE:** the slab apply snapshots (the `TestT4bPreT4RegressionSnapshot`
  slab arms) are ALREADY in place and pin the 1-D path; the DD 1-D scan-march
  does NOT change under D5a (1-D rides CumprodScan/ODW, not the 2-D row-march),
  so the slab arms stay green untouched. The NEW cart2d arm is the 2-D-specific
  frozen reference.
- **Tolerance:** `assert_regression` `kind=direct reduction_depth=mesh.nx`
  (iterative? NO — single matvec → direct; `reduction_depth` = the cell-chain
  depth per [[feedback_regression_tolerance_design]]).
- **Markers:** the regression-snapshot machinery is `@pytest.mark.foundation`
  (existing class).
- **Failure mode caught:** the "both paths move together" blind spot of
  D5a.1 — a frozen-old-convention reference makes a uniform value shift
  unmissable (precedent: the Step-A §3 positive gate).

### GATE D5a.3 — DD strict bit-identity (the SCAN/SWEEP negative control)

**Claim:** the SOLVE-path scan snapshots stay byte-identical (DD scan is
power-of-2-exact: `0.5·in+0.5·out ≡ 0.5·(in+out)`, `2·QV·inv ≡ QV·inv/0.5`).
The fold must NOT leak into the solve path.

**Pin:** `python -O -m pytest tests/sn/sweep/core tests/sn/solve -W
"error::tests.sn.regression._regression_assert.DriftWarning"` MUST stay
**505 passed / 1 skipped / 4 xfailed** (confirmed live this session). The
`tests.sn.regression...` path is load-bearing (`orpheus.sn...` silently fails
to escalate — FINDING-1, [[issue_206_phase_c_verification]]). Re-confirm the
4 xfailed ids are the SAME pre-existing reds (#206 cyl-matvec etc.).

- **If RED on a SWEEP snapshot** → the fold leaked into solve (a bug — DD
  scan must stay byte-identical). **If RED on an APPLY snapshot** → that
  snapshot needed D5a.2 migration (re-association on the matvec). The gate
  DISCRIMINATES the two failure classes.
- **Failure mode caught:** #3 missing factor / #6 convention drift that
  changes the DD scan reduction tree (would shift the byte-identical solve).

### GATE D5a.4 — Step polymorphic-READINESS placeholder (OPTIONAL, tracked)

**Claim:** the genericised ScanMarch interior consumes `affine_scan_coefficients`
GENERICALLY (no inline DD, no scheme name) → when Step lands, Step-2-D rides
the row-march with no further ScanMarch change.

**File:** an `xfail(reason="Step scheme not implemented — #NNN")` placeholder
in `test_scan_march_equivalence.py` asserting `ScanMarch.sweep` with a Step
`scheme` ≡ FullFieldWavefront. **OPTIONAL** — Step is out of D5's scope; this
just TRACKS the readiness so it flips to xpass when Step lands. RECOMMEND only
if the campaign intends Step soon; otherwise the §6 audit (zero inline-DD in
ScanMarch) is sufficient evidence of polymorphic-readiness.

---

## ════ D5b — N-DIM LD ON THE DAG WAVEFRONT (#158 Inc D / #38) — THE NEW MATH ════

**SUT:** supply the multi-axis bilinear LD Schur kernel so
`LinearDiscontinuous.cell_kernel_batch`/`residual_kernel_batch` work for
`len(s_axes) ≥ 2` (close the `linear_discontinuous.py:430-436` raise). LD's
multi-D closure is BILINEAR (an independent slope per streaming axis) — it
does NOT fit the scan-march's x-scan + transverse-direct-y, so LD-in-N-D rides
the DAG wavefront (`FullFieldWavefront`/`MovingFrontierWindow`, already pure
scheme-delegators via `sweep_graph.py` → `scheme.cell_kernel_batch`/
`residual_kernel_batch`). The 1-D Schur intermediates (`_LDCellTerms` `:223`,
`_schur_terms` `:295`, `_kernel_terms` `:411`) are the structural model; the
multi-D kernel generalises the 2×2 to a (1+d)×(1+d) per-cell system (average +
one slope per axis).

**⚠ The multi-D LD discretization FORMULATION (the bilinear closure equations)
comes from a SEPARATE literature-researcher dispatch** — D5b's MMS does NOT
need it (MMS is solution-driven: pick ψ, derive the source from the KNOWN
within-group transport PDE, structurally independent of the LD kernel — L11).

**Bit-identity posture:** NEW capability — nothing to be identical to. The
gates are MMS (convergence-order/flux-shape) + two-paths (software invariant).

### GATE D5b.0 — INVERT the N-D raise pin (the test that must be RETIRED)

**Claim:** the existing pin `test_linear_discontinuous.py::TestLDKernel::
test_cell_kernel_batch_rejects_multi_d` (`:361-367`) asserts
`pytest.raises(NotImplementedError, match="d=1")` for `s_axes=(two, two)`.
When D5b lands, this raise is GONE — **the test must be RETIRED/INVERTED**.

**Action:** REPLACE it (do NOT leave it asserting the old raise — a FALSE
gate, anti-pattern: asserting a contract the carve deliberately reverses):
```python
def test_cell_kernel_batch_admits_multi_d(self):
    """D5b: the multi-D bilinear LD kernel runs (was: raised d=1-only).
    A d=2 s_axes now solves the (1+2)×(1+2) per-cell bilinear system."""
    # 2-D LD cell — runs, returns (psi_avg, (out_x, out_y)) with the
    # bilinear slope eliminated per axis.
    psi_avg, faces = ld.cell_kernel_batch(psi_in=(in_x, in_y),
                                          s_axes=(s_x, s_y), ...)
    assert len(faces) == 2          # one outgoing face per streaming axis
```
- **Markers:** `@pytest.mark.foundation`.
- **NOTE the partner raise** `affine_scan_coefficients` is NOT inverted — LD's
  2-D path is the WAVEFRONT, not the scan; LD stays NON-scan-march-compatible
  in d≥2 (the D5-0 EXCL). The scan-coefficients raise on non-neutral curvature
  stays (curvilinear LD unpublished). Only the DAG-kernel `len(s_axes)!=1`
  raise inverts.

### GATE D5b.1 — multi-D LD round-trip (lockstep, per-cell, foundation)

**Claim:** software invariant (foundation). The multi-D `cell_kernel_batch`
(solve) and `residual_kernel_batch` (apply) describe the SAME bilinear per-cell
system: `residual_kernel_batch` at the `psi_bar` `cell_kernel_batch` solves
for vanishes to FP noise.

**File:** `tests/sn/spatial/test_linear_discontinuous.py::TestLDKernel` — add
`test_residual_zero_at_solved_cell_avg_2d` mirroring the existing
`test_group1_equals_group2_flat`/round-trip structure but for `s_axes=(s_x, s_y)`.
- **Config:** `n_groups ∈ {1, 2}` HETEROGENEOUS per-group `sigt_cells`
  (1G degenerate control + 2G real gate, vv H1). NON-flat per-axis `psi_in`
  (the bilinear slope per axis must be ACTIVE — flat in_x/in_y nulls the
  slope, hiding the multi-D coupling). The round-trip must feed back the FULL
  solved state (if the multi-D LD returns extra slope moments, ALL must round-
  trip — a partial feed passes spuriously).
- **Tolerance:** `atol=1e-12` (the (1+d)×(1+d) solve is a few division ULP
  deep — bump from DD's `1e-13` by the reduction depth; DOCUMENT the bound).
- **Term-activation declaration (Mode 7):** the round-trip ACTIVATES streaming
  (both axes) + collision + the per-axis slope; NULLS nothing.
- **Failure mode caught:** #4 wrong recursion (slope-elimination Schur index
  drift), #2 variable swap (x-slope ↔ y-slope).

### GATE D5b.2 — multi-D LD MMS O(h²) (leg 1, the structurally-indep ref)

**Claim:** convergence-order (math, L1) + flux-shape. The multi-D LD kernel,
run via `FullFieldWavefront` (the d-generic DAG oracle), matches the
manufactured 2-D solution at O(h²). **MMS = the math/flux-shape pillar; it
does NOT and CANNOT prove an eigenvalue (vv hierarchical taxonomy) — there is
NO LD eigenvalue gate.**

**File:** `tests/sn/verification/mms/test_mms_ld_slab.py` is 1-D; add a NEW
file `tests/sn/verification/mms/test_mms_ld_2d.py` (the 2-D LD home).
- **Driver:** `solve_sn_fixed_source(materials, Mesh2D(...), quadrature, Q,
  scheme=LinearDiscontinuous())` — the `scheme=` kwarg threads end-to-end
  (confirmed `solver.py:1972` → `_as_sn_mesh(... scheme=scheme)`). The 2-D LD
  mesh routes to `FullFieldWavefront` (after D5-0 + D5b: NOT ScanMarch).
- **Ladder:** `nx=ny ∈ {20, 40, 80}` (or non-square `(nx,ny)` per pair) →
  assert `orders[-1] > 1.95 ∧ all > 1.85` + the VALUE band against
  `phi_exact` (vv §5: rate ≠ correctness; the manufactured solution IS the
  structurally-independent reference).
- **The MMS DESIGN is §3 below** (μ-non-trivial, Q̂≠0 slope source — the
  Mode-7 override is MANDATORY here, NOT optional).
- **Markers:** `@pytest.mark.l1 @pytest.mark.slow @pytest.mark.verifies(
  "ld-cartesian-2d", "transport-cartesian-2d")` — **`ld-cartesian-2d` is a
  NEW label D6 must mint** (see §6).
- **Failure mode caught:** #1 sign flip (the multi-D slope-row sign — LM-1989
  trap, now per-axis) diverges under refinement; #3 missing factor breaks the
  value band; #5 index error (non-uniform mesh detectably wrong).

### GATE D5b.3 — multi-D LD two-paths: FFW ≡ MFW (leg 3, the headline software invariant)

**Claim:** software invariant (foundation). The multi-D LD discretization run
via `FullFieldWavefront` (anti-diagonal wavefront) ≡ via `MovingFrontierWindow`
(rolling frontier) — the SAME LD via two different DAG SCHEDULES. (Both are
pure scheme-delegators → they call the SAME `cell_kernel_batch`; the agreement
proves the WALK storage policy doesn't perturb LD's cell math.)

**File:** `test_mms_ld_2d.py::test_ld_2d_two_paths_ffw_equals_mfw`.
- **Config (Mode-9 degeneracy-break — MANDATORY):** the §3 stress ansatz
  (het Σ_t(x,y), 2G asymmetric downscatter, μ-non-trivial, NON-SQUARE
  `nx≠ny ≈ 40`). A flat/square/1G config makes the two reduction trees
  accidentally FP-coincident (the schedule-equivalence Mode-9 sub-case —
  [[issue_158_ld_spatial_verification]] B6).
- **Tolerance:** end-to-end converged scalar flux →
  `assert_allclose(rtol=SAFETY × inner_tol, atol=...)` ≈ `rtol=1e-9`
  (iterative → `SAFETY×conv_tol`, [[feedback_regression_tolerance_design]]).
  NOT `array_equal` (FFW and MFW gather/scatter faces in different orders →
  different reduction trees). RECOMMEND ALSO a single-sweep
  `assert_array_almost_equal_nulp(nulp=n_cells)` (fixed-depth, sharper).
- **PAIR with D5b.2** (both paths independently match the MMS ref — the
  structural ground; the two-paths nULP alone is necessary-not-sufficient).
- **VERIFY each leg ran its rep** (assert `default_for(ld_2d_mesh)` is the
  expected wavefront; force MFW explicitly for the oracle leg) — a two-paths
  gate where both legs secretly ran the same rep is a silent false green
  (Mode-8-adjacent).
- **Markers:** `@pytest.mark.l1 @pytest.mark.verifies("ld-cartesian-2d")`.
- **Failure mode caught:** #4 wrong recursion (walk-order-dependent slope
  drift), the schedule-equivalence Mode-9 sub-case.

### GATE D5b.4 — multi-D LD scan-matvec ≡ scan-sweep (leg 3, the matvec twin)

**Claim:** software invariant (foundation). On the 2-D LD wavefront, LD's
forward matvec (`FullFieldWavefront.loss_action` → `residual_kernel_batch`)
≡ LD's SI sweep (`FullFieldWavefront.sweep` → `cell_kernel_batch`). I.e.
Krylov ≡ SI on the 2-D LD path.

**File:** `test_mms_ld_2d.py::test_ld_2d_krylov_matches_si` (mirror the 1-D
`test_sn_1d_slab_ld_mms_krylov_matches_si`).
- **⚠ HAZARD:** the matvec reconstructs the outgoing face via
  `outgoing_face_from_average` — the multi-D LD must supply a working
  `outgoing_face_from_average` in BOTH the sweep AND the matvec direction
  (per axis). A scan/wavefront SWEEP without a verified MATVEC is half a carve
  (the L14 leg-3 standoff demands both directions). If LD's multi-D
  `residual_kernel_batch` reconstruction is wrong, Krylov can't converge or
  diverges from SI.
- **Tolerance:** `rtol=1e-9, atol=1e-11` (the 1-D LD Krylov≡SI tolerance).
- **Config:** the §3 stress config (≥2G het).
- **Markers:** `@pytest.mark.l1 @pytest.mark.verifies("ld-cartesian-2d")`.
- **Failure mode caught:** #6 convention drift (`Lψ` vs `(L+C)ψ` in the
  matvec), #4 wrong recursion in the apply-direction face reconstruction.

### GATE D5b.5 — DD-vs-LD routing-flip (the #158-B3 analogue — different converged values)

**Claim:** DD and LD on the SAME 2-D heterogeneous 2G config converge to
DIFFERENT values, each anchored to its OWN absolute reference. (A test that
asserts DD≡LD would be WRONG — they are different discretizations; a test
that only checks "LD runs" is blind to LD silently computing DD, the exact
D5-0 misroute symptom.)

**File:** `test_mms_ld_2d.py::test_dd_and_ld_2d_converge_to_different_values`.
- **Construction:** drive `solve_sn_fixed_source` TWICE on the SAME mesh/
  materials/source — once `scheme=DiamondDifference()`, once
  `scheme=LinearDiscontinuous()`. Assert the converged scalar fluxes DIFFER
  by MORE than the solver tolerance at COARSE mesh (where the O(h²) closures
  diverge): `assert not np.allclose(phi_dd, phi_ld, rtol=1e-3)`. Anchor each:
  LD → its O(h²) MMS (D5b.2); DD → its existing 2-D MMS / k_inf anchor.
- **Why:** this is the catcher that LD is GENUINELY computing LD, not DD. If
  the D5-0 misroute regressed (LD silently ran ScanMarch inline DD), this
  test would see DD≡LD and FAIL. It is the cross-thread guard tying D5-0
  honesty to D5b correctness.
- **Markers:** `@pytest.mark.l1 @pytest.mark.foundation` (a discrimination
  invariant, not an equation claim) — actually `@pytest.mark.foundation`
  alone (no `verifies` — it's a "these are different schemes" contract).
- **Failure mode caught:** the D5-0 misroute REGRESSION (LD→DD silent
  collapse), Mode #2 variable swap (scheme dispatch swapped).

---

# 3. THE MULTI-D LD MMS DESIGN (the Mode-7 override — MANDATORY)

**CRITICAL (vv Mode 7 + MMS operational rules):** the existing 2-D MMS
(`SN2DCartesianMMSCase`, `phi_exact = sin(πx/Lx)·sin(πy/Ly)`) is
ISOTROPIC-IN-μ by construction (the manufactured ψ has no μ dependence; the
source `external_source` builds `streaming = mu_x·dA/dx + mu_y·dA/dy` but the
flux itself is angularly flat). **An isotropic/flat ansatz NULLS the per-axis
slope rows of the bilinear LD closure — it tests NOTHING about the multi-D
slope coupling, the exact thing D5b introduces.** FORBID it. The simplification
bias must be overridden at write-time (AI has no SymPy derivation cost — pick
the stress-test ansatz).

## The trial ψ(x, y, μ) — with the activated/nulled-term declaration

```
ψ_{n,g}(x,y) = [ A_g(x,y)  +  μ_{x,n}·B_g(x,y)  +  μ_{y,n}·C_g(x,y) ] / W
```
with, per group g:
```
A_g(x,y) = a0_g + a1_g·sin(πx/Lx)·sin(πy/Ly) + a2_g·cos(2πx/Lx)·cos(3πy/Ly)
B_g(x,y) = b0_g + b1_g·sin(2πx/Lx)·sin(πy/Ly)        # x-direction slope driver
C_g(x,y) = c0_g + c1_g·sin(πx/Lx)·sin(2πy/Ly)        # y-direction slope driver
```

**Activated / nulled declaration (vv Mode 7 — the ansatz MUST declare this):**

| Term | Activated by | Why load-bearing |
|------|-------------|------------------|
| **x-axis LD slope** | `μ_x·B_g(x,y)` with `b1≠0` + mixed x-harmonic in A | the bilinear closure's x-slope row sees a genuinely x-varying, μ_x-weighted field — the NEW multi-D coupling |
| **y-axis LD slope** | `μ_y·C_g(x,y)` with `c1≠0` + mixed y-harmonic | the y-slope row, INDEPENDENT of x (the bilinearity: the two slopes do not collapse) |
| **cross-coupling** | `cos(2πx)·cos(3πy)` in A (MIXED scales, different x/y harmonics) | breaks any x↔y symmetry that a square-symmetric ansatz would hide (x↔y swap defense) |
| **boundary closure** | `a0_g > 0` (NON-vanishing at all 4 edges) | real vacuum-inflow BC the LD boundary closure must satisfy with non-trivial interior ([[phase4_46_nonvacuum_mms_ansatz]]: prior MMS all VANISH at boundary → BC untested) |
| **group coupling** | per-group `a/b/c` coeffs + 2G asymmetric downscatter | mode #6 convention drift in the group axis |
| NULLED: nothing | — | the ansatz activates streaming(both axes) + collision + both slopes + BC + group coupling. NO term covered by an active ERR-NNN is nulled. |

**Concrete parameter discipline:**
- `a0_g > 0` (boundary non-vanishing — LOAD-BEARING).
- MIXED spatial scales (`k=1` and `k=2,3` harmonics, different per axis) → the
  coarse-mesh O(h²) error is high-frequency-dominated → STRESSES the closure
  + gives the ladder real dynamic range (NOT floor-limited).
- `μ_x·B + μ_y·C` makes the per-ordinate field genuinely μ-dependent (the
  bilinear LD slopes per axis see a NON-flat field — the isotropic sin·sin
  ansatz is BLIND to this).
- 2G, HETEROGENEOUS Σ_t(x,y) (reuse `_default_hetero_2d_xs_functions`,
  `derivations/continuous/mms/sn.py:950`, Σ_a>0 guaranteed) + non-trivial
  per-group coeffs (downscatter coupling, mode #6).
- NON-SQUARE domain (`Lx ≠ Ly` or `nx ≠ ny`) — x↔y swap defense.

## Derivation approach (SymPy — algebra-of-record Branch 1; L11; NEVER
## hand-transcribe)

A NEW `SN2DCartesianLDStressMMSCase` (or `..._aniso`) in
`derivations/continuous/mms/sn.py` (mirror `SN2DCartesianMMSCase` `:642`):
substitute the ansatz into the 2-D within-group transport PDE
```
μ_x ∂_x ψ + μ_y ∂_y ψ + Σ_t ψ = (1/W)(Σ_s^T φ + Q^ext)
```
and solve SYMBOLICALLY for `Q^ext_n(x,y)`. The scalar flux
`φ_g = ∫ψ dμ` → the `μ_x·B + μ_y·C` terms integrate to zero over a symmetric
quadrature → `φ_g = A_g(x,y)` (the isotropic moment), while the streaming
derivative carries the FULL ansatz (the μ-weighted slope drivers survive in
the per-ordinate source). Pin the symbolic source with a
`@pytest.mark.foundation` derive-test BEFORE consuming it. The source is
**structurally independent of the LD kernel** (derived from the continuous
PDE, not from any cell-update call — L11).

## Expected convergence order, slope-source Q̂≠0, strategy

- **Order:** O(h²) (LD is second-order). Assert `orders[-1] > 1.95 ∧ all > 1.85`.
- **⭐ Q̂≠0 REQUIREMENT (the LM-1989 slope-sign trap):** the manufactured ψ has
  a NON-ZERO slope moment per axis (`∂ψ/∂x`, `∂ψ/∂y` non-constant per cell).
  Unlike the 1-D Increment-B flat-source posture (Q̂=0, slope UNKNOWN solved),
  the multi-D LD MMS SHOULD supply the SLOPE-MOMENT source so the slope-SOURCE
  sign — the OTHER half of the LM-1989 §1.4/§6 trap, untested at Q̂=0 — is
  exercised. **This depends on whether D5b threads the moment source** (the
  1-D #158 Increment C). If D5b is flat-source (Q̂=0, like 1-D Inc B), the MMS
  runs flat-source and the slope-SOURCE sign stays untested (document the
  scope: "multi-D slope-UNKNOWN verified, slope-SOURCE deferred to the
  moment-source increment"). If D5b threads the moment source, the MMS MUST
  supply Q̂≠0 with a non-trivial per-axis slope (the brief's requirement).
  **Flag to method-implementer: declare the moment-source posture; the MMS
  Q̂ follows it.** The non-vanishing per-axis slope drivers (`B`, `C` with
  spatial variation) make Q̂≠0 naturally — the MMS is slope-source-READY.
- **Strategy:** runs on the WAVEFRONT (`FullFieldWavefront` oracle for the
  O(h²) gate; the two-paths gate D5b.3 adds `MovingFrontierWindow`). NOT
  ScanMarch (LD-in-2-D is the structural EXCL).

## The new theory label (flag to D6)

The MMS `@verifies("ld-cartesian-2d")` — **`ld-cartesian-2d` does NOT exist**
as a `:label:` block (see §6). D6 (archivist) MUST mint it in
`docs/theory/discrete_ordinates.rst` (the LD theory section is currently a
`.. todo:` stub `:1498`).

---

# 4. IMPLEMENTATION-SEQUENCE RECOMMENDATION

**The ordering is risk-driven (the §2 thread table):**

### Phase 0 — literature-researcher dispatch (BEFORE D5b code)
The multi-D bilinear LD FORMULATION (the closure equations: the (1+d)×(1+d)
per-cell system, the per-axis slope rows, the upwind-discontinuous multi-D
face closure). This is the one piece D5 does NOT have. Dispatch BEFORE D5b
implementation (NOT before D5-0/D5a, which need no new formulation). The MMS
(§3) does NOT depend on it (solution-driven) — so the MMS case + the D5b
gates can be WRITTEN in parallel with the literature dispatch.

### Phase 1 — cross-domain-attacker dispatch (the new L1 reference frame check)
The multi-D LD MMS is a NEW L1 reference. Per the proactive-dispatch protocol
(subagent-handoff-protocol §"New L0/L1 reference solver derivation"), dispatch
cross-domain-attacker to check the MMS's structural frame BEFORE the SymPy ink
dries: does the 2-D Cartesian ansatz lift cleanly? Is the μ_x·B + μ_y·C slope
structure the right frame for a BILINEAR (not biquadratic) closure? Cheap
insurance against shipping a reference that under- or over-stresses the
bilinear coupling. Dispatch in PARALLEL with Phase 0 (both feed D5b).

### Phase 2 — D5-0 routing honesty (LAND FIRST of the code)
Close the CONFIRMED-LIVE silent misroute BEFORE any math moves. After D5-0, a
2-D LD mesh routes to the wavefront (which RAISES pre-D5b, HONEST) instead of
silently running ScanMarch inline DD. This is the cheapest, highest-value fix
(a one-trait predicate change + the negative routing gate). **Rationale:**
shipping D5a/D5b on top of a live misroute means every 2-D LD test would have
to fight the misroute; closing it first makes the subsequent math land on an
honest routing substrate.

### Phase 3 — D5a DD scan-march fold (stabilise SM before D5b)
The bounded, value-preserving fold. It re-bases the ScanMarch interior onto
the coefficient model (DD+Step generic). LAND before D5b's new math because:
(a) it's value-preserving (low risk — the §2-gates catch a 2× error loud);
(b) it makes ScanMarch scheme-generic, so the D5-0 EXCL (LD never reaches
inline DD) is enforced by the ABSENCE of inline DD, not just by `supports()`;
(c) D5b's two-paths gates run on a stabilised wavefront/scan substrate.

### Phase 4 — D5b multi-D LD bilinear (the new capability)
The genuinely new math (needs Phase 0's formulation). Lands the multi-D LD
kernel + inverts the raise pin + the MMS + two-paths + Krylov≡SI + routing-flip.

### Phase 5 — D6 archivist dispatch (theory labels + narrative)
Mint the orphan + new labels (§6); expand the LD theory `.. todo:` stub into
the rich derivation (the multi-D bilinear closure, the Step/DD/LD advection-
scheme narrative). Per algebra-of-record, the archivist expands AFTER the
SymPy MMS module + the D5b code land.

**Summary sequence:** `[lit-researcher ∥ cross-domain-attacker] → D5-0 →
D5a → D5b → D6`. The two dispatches are parallel and feed D5b; D5-0 and D5a
need no new formulation and can start immediately.

---

# 5. EXISTING TESTS — RETIRE / INVERT / MIGRATE / STAY-GREEN

| Test | file:line | Verdict | Reason |
|------|-----------|---------|--------|
| `test_cell_kernel_batch_rejects_multi_d` | `test_linear_discontinuous.py:361` | **INVERT** (D5b.0) | asserts the `d=1`-only raise; D5b removes it → rewrite to `test_cell_kernel_batch_admits_multi_d`. Leaving it is a FALSE gate. |
| `TestD3SupportsMatrix._fake` | `test_unified_sweep_dispatch.py:215` | **MIGRATE** (D5-0) | always supplies `is_affine_scannable=True` — add the NEW `is_scan_march_compatible` trait so the d=2 admission pin exercises the right trait; add an LD-shaped fake (or real 2-D LD mesh) for the negative gate. |
| `test_scan_march_keeps_d2_and_1d_coverage` | `test_unified_sweep_dispatch.py:237` | **STAY GREEN** (D5-0.2) | the narrowing must NOT lose DD-2-D / 1-D coverage. Re-run; update the fake's trait. |
| `test_ld_slab_mesh_routes_to_cumprod_scan` | `test_mms_ld_slab.py:58` | **STAY GREEN** | 1-D LD routing unchanged by D5 (1-D scan trait kept). |
| `test_scanmarch_sweep_equals_oracle` / `_residual_equals_oracle` | `test_scan_march_equivalence.py:104,179` | **STAY GREEN** (D5a.1) | the DD 2-D fold must not change the value; these ARE the headline D5a catchers. Confirm het/2G/non-square parametrize. |
| `TestT4bPreT4RegressionSnapshot` (slab arms) | `test_streaming_operator.py` | **STAY GREEN** + EXTEND (D5a.2) | slab arms pin 1-D (untouched); ADD a cart2d apply arm (frozen pre-#240 snapshots). |
| DD strict gate (`sweep/core + solve -W error::DriftWarning`) | (invocation) | **STAY 505/1/4** (D5a.3) | the SCAN/SWEEP byte-id negative control. |
| `test_ld_two_paths_scan_equals_dag_oracle` | `test_mms_ld_slab.py:151` | **STAY GREEN** | 1-D LD two-paths unchanged. |
| `test_ld_curvilinear_scan_rejected` | `test_mms_ld_slab.py:202` | **STAY GREEN** | curvilinear LD EXCL unchanged (D5 does not touch curvilinear). |
| `test_ld_thick_diffusive_limit_xfail` | `test_mms_ld_slab.py:235` | **STAY xfail** | the diffusion-limit tripwire; D5 (flat-source posture) does not fix it. Do NOT claim D5 closes it. |
| `test_mms_2d` (DD) | `test_mms_2d.py:46,92` | **STAY GREEN** | DD 2-D MMS — the structurally-indep value reference D5a must preserve. |
| `test_kinf_homogeneous` (≥2G) | `analytical/test_kinf_homogeneous.py` | **STAY GREEN** | the ONLY structurally-indep EIGENVALUE anchor; 1G DEGENERATE (Cardinal Rule). D5a value-preservation rides on it. |

**Anchors that MUST stay green throughout (the green floor):**
`test_kinf_homogeneous`(≥2G) + `test_mms_2d` + `test_mms_ld_slab` + the 2-D
oracles (this session: 56p/3xf) + the DD strict gate (505/1/4) + the LD suite
(24p/1xf).

---

# 6. THE ORPHAN-LABEL PROBLEM (what D6 owes)

Confirmed live this session (`grep ":label:" docs/theory/discrete_ordinates.rst`):

**Labels that EXIST** (the green substrate):
- `transport-cartesian` (`:275`), `transport-cartesian-2d` (`:292`)
- `dd-cartesian-1d` (`:523`), `dd-cartesian-2d` (`:562`)

**ORPHAN labels** (`@verifies`'d but NO `:label:` block — D6 must mint):
- **`ld-cartesian-1d`** — `@verifies`'d at `test_mms_ld_slab.py:81,122,150`;
  NO backing label. The LD theory section is a `.. todo:` stub (`:1498`).
- **`ld-slab`** — `@verifies`'d at `test_mms_ld_slab.py:81`; NO backing label.

**NEW labels D5 tests need** (D6 must mint, do NOT silently depend on them):
- **`ld-cartesian-2d`** — the multi-D LD MMS (`test_mms_ld_2d.py`, D5b.2/.3/.4)
  `@verifies` it. Brand new (no LD 2-D theory exists). D6 mints it WITH the
  multi-D bilinear closure derivation.

**Design discipline:** the D5b LD 2-D tests `@verifies("ld-cartesian-2d",
"transport-cartesian-2d")` — `transport-cartesian-2d` EXISTS (the generic
2-D PDE), so even before D6 mints `ld-cartesian-2d` the test links to ONE real
label (no silent total-orphan). Per the V&V harness, a `verifies` on a
non-existent label is a SOFT gap (the audit flags it) — NOT a test failure —
so the D5b tests collect and pass; the label gap shows in
`python -m tests._harness.audit` and is D6's close-out. **Flag to D6 (the
archivist owes 3 labels):** `ld-cartesian-1d`, `ld-slab` (retroactive, Inc-A
orphans), `ld-cartesian-2d` (new, D5b). The Step rows, if ever built, owe
`step-cartesian-1d`/`step-cartesian-2d` — NOT D5's concern.

---

# 7. FAILURE-MODE COVERAGE SUMMARY (vv §6 + numerical-bug-signatures)

| Failure mode | D5 gate that catches it |
|--------------|------------------------|
| #1 sign flip (LD per-axis slope row; D5a transverse `−sy` sign) | D5b.2 MMS diverges under refinement; D5a.1 oracle |
| #2 variable swap (x-slope↔y-slope; scheme dispatch swap) | D5b.1 round-trip (non-flat per-axis); D5b.5 routing-flip; NON-SQUARE config throughout |
| #3 missing/extra factor (the `2` lost in D5a fold; LD slope weight) | D5a.1 oracle (clean factor → loud); D5a.2 frozen snapshot; D5a.3 strict DD gate |
| #4 wrong recursion (LD Schur slope-elim index; D5a coefficient recurrence) | D5b.1 round-trip; D5b.3 two-paths; D5b.4 Krylov≡SI |
| #5 index error (multi-D cell-chain reorder) | D5b.2 MMS (non-uniform mesh); D5b.3 |
| #6 convention drift (D5a `α/β` vs coeff-model; LD `Lψ` vs `(L+C)ψ` matvec) | D5a.3 strict gate; D5b.4 matvec; D5-0 trait-honesty |
| **Mode 7 MMS simplification bias** | §3 stress ansatz override (μ-non-trivial `μ_x·B+μ_y·C`, a0>0, mixed scales, 2G het, NON-SQUARE) — MANDATORY, the existing sin·sin 2-D MMS is FORBIDDEN for LD |
| Mode 8 `-O` strips bare assert | ALL D5 gates use `np.testing.assert_*` / `assert_regression` / `pytest.fail` — NEVER bare `assert`. Canonical `-O`. (The §5 dispatch tests use `pytest.fail` — `-O`-safe.) |
| **Mode 9 splitting/SCHEDULE/representation equivalence in degenerate regime** | D5b.3 (FFW≡MFW schedule-equiv) + D5a.1 (ScanMarch≡FFW) run on het/2G-asym/non-flat/NON-SQUARE — breaks the FP-coincidence. THE D5-specific Mode-9 hazard: a two-paths/two-schedules gate on a flat-square-1G box is degenerate-blind. |

## ⭐ Mode-9 GENERALIZATION (self-improvement trigger — flag at delivery)

D5b.3 (FFW≡MFW) is a SCHEDULE-equivalence FP claim (two DAG walk policies,
same cell math), NOT a splitting. Per [[issue_158_ld_spatial_verification]]
the Mode-9 skill row's scope ("a splitting") is NARROWER than the phenomenon
— the SAME degenerate-regime-blindness applies to splitting OR acceleration OR
SCHEDULE OR REPRESENTATION equivalence. This is the SECOND D5-area instance
(the 1-D Inc-B B6 was the first). **RECOMMEND sharpening the vv-principles
Mode-9 row: "a splitting/acceleration" → "a splitting / acceleration /
schedule / representation equivalence".** Already in PRACTICE (the #206
Phase-C scan-vs-DAG nULP gates, the Inc-B two-paths) — this just names the
principle. NO new ERR (no caught bug yet) — a skill-row scope sharpen, not an
error-catalog entry.

## NEW failure-mode-table row consideration (self-improvement)

**A SELECTION-honesty bug (the D5-0 misroute) is a NEW test-design failure
shape not in the vv table:** a `supports()` predicate whose trait answers a
DIFFERENT question than the one the strategy needs ("1-D affine-scannable" vs
"N-D scan-march-compatible") → a scheme is admitted to a strategy that runs
the WRONG cell math (inline DD on an LD mesh) SILENTLY. This is distinct from
Mode 7 (MMS bias — the test can't see the bug) and Mode 9 (degenerate regime
— the gate is blind): here the SOLVER is wrong (LD silently computes DD) AND
no test exists to catch it (the existing `_fake` always supplies the
admitting trait). It is a SELECTION/DISPATCH-honesty failure: a polymorphic
dispatch that conflates two predicates. **This is caught by D5-0.1 (the
negative routing gate) + D5b.5 (the DD≠LD discrimination gate).** If a real
solver bug is later shown to have hidden behind this conflation in production,
log an ERR-NNN (failure mode #6 convention-drift variant: trait/predicate
conflation). For now: documented here as a test-design pattern. **Flag at
delivery: this is a candidate Mode-10 ("predicate-conflation silent
misroute") IF it recurs — currently a single instance, so it stays a
documented row here, not a skill promotion.** Per the self-improvement
directive, promote to the vv-principles failure-mode table only on the
SECOND instance or a caught production bug.

---

# 8. DELIVERY CHECKLIST (the implementer's gate sequence)

1. **Phase 0/1 dispatches** — lit-researcher (multi-D LD formulation) ∥
   cross-domain-attacker (MMS frame). Both feed D5b; start immediately.
2. **D5-0** — add `is_scan_march_compatible` trait (DD/Step True, LD False);
   narrow `ScanMarch.supports` d≥2 arm. Gate D5-0.1 (negative: 2-D LD ∉
   ScanMarch) + D5-0.2 (positive: 2-D DD ∈ ScanMarch, 1-D LD ∈ CumprodScan
   stay green). VERIFY the live misroute is closed (re-run the probe).
3. **D5a** — fold ScanMarch `_sweep_interior` + `_loss_action_interior` onto
   `affine_scan_coefficients`/`source_emission`/`cell_average`/
   `outgoing_face_from_average`. Gates D5a.1 (oracle stays green, het/2G/
   non-square) + D5a.2 (cart2d frozen snapshot, nulp) + D5a.3 (strict DD
   505/1/4). Audit §6 grep: zero inline `2*s`/`α=`/`β=` in ScanMarch bodies.
4. **D5b** — multi-D LD bilinear kernel (close `:430` raise). Gate D5b.0
   (INVERT the raise pin) + D5b.1 (round-trip, non-flat per-axis) + D5b.2
   (MMS O(h²), the §3 stress ansatz) + D5b.3 (FFW≡MFW two-paths) + D5b.4
   (Krylov≡SI matvec) + D5b.5 (DD≠LD discrimination).
5. **Green floor** — `test_kinf_homogeneous`(≥2G) + `test_mms_2d` +
   `test_mms_ld_slab` + 2-D oracles (56p/3xf) + DD strict gate (505/1/4) +
   LD suite (24p/1xf) all stay green.
6. **D6 dispatch** — mint `ld-cartesian-1d`, `ld-slab`, `ld-cartesian-2d`;
   expand the LD theory `.. todo:` stub.
7. Route around the 7 pre-existing reds; NEVER run all `tests/sn` (#212).

---

## Cross-links
- Extends [[issue_158_ld_spatial_verification]] (Inc B = the 1-D LD scan; D5b
  = the deferred Inc D multi-D bilinear; the Mode-9 schedule-equiv
  generalization first flagged there) + [[issue_240_phase2_s_axes_scanmarch_spec]]
  (the s_axes convention — LANDED; D5a is the coefficient-model lift on top;
  the §2a cart2d frozen-snapshot precedent) + [[issue_240_phase2_step_b_sigma_leak]]
  (the loss_action carve precedent) + [[issue_206_phase_c_verification]]
  (density-vs-scan nULP) + [[feedback_regression_tolerance_design]]
  (iterative→SAFETY×tol; direct→reduction_depth) + [[phase4_46_nonvacuum_mms_ansatz]]
  (a0>0 boundary-non-vanishing).
- coding-elegance Pattern 4 (the `is_scan_march_compatible` TRAIT — make the
  illegal "LD on the row-march" unrepresentable; NOT an isinstance) + Pattern 7
  (the affine coefficients single-source the fold) + anti-pattern #4 (no
  stringly/isinstance scheme dispatch in `supports`).
- vv-principles Mode 7 (the §3 MMS override) + Mode 9 (the schedule-equiv
  sub-case, generalization flagged) + hierarchical taxonomy (LD = convergence-
  order/flux-shape; NO eigenvalue gate; MMS can't prove eigenvalues) +
  Bit-identity-vs-principled (D5a = principled re-baseline, 3 criteria).
