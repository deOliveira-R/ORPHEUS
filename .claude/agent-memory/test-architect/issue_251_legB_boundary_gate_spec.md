---
name: issue-251-legB-boundary-gate-spec
description: "#251 (Leg B of #247) L1 verification spec — close the BOUNDARY half of the LM-1989 trap in 2-D Cartesian LD SN: the manufactured prescribed inflow varies transversely along each face but _inflow_to_moments zeros the 2^{d-1} transverse face-slope. Structural lift-threading teeth + consumption mutation + face-projection foundation sub-gate + bit-identity guard. ⚠ improves-on-flat is NOT achievable (sub-floor). Branch refactor/sn-foundation-cleanup. PRE-IMPL."
metadata:
  type: project
---

# #251 — 2-D Cartesian LD BOUNDARY transverse-face-slope gate (Leg B of #247)

Issue #251 (split out of #247), L1 VERIFICATION. Branch
`refactor/sn-foundation-cleanup`. Host `.venv/bin/python -O`. I design the gate;
the production change (a moment-resolved boundary trace — widen
`_inflow_to_moments` to carry the transverse face-slope) is owned by the
method-implementer. Extends [[issue-247-moment-source-gate-spec]] (Leg A — the
BULK external+scattering slope source, COMMITTED `d9396a2`) +
[[issue-247-legA-moment-source-closeout]] (the proven structural-teeth
template) + [[feedback_vv_tagging]] + [[feedback_regression_tolerance_design]].

## 0. THE PROBED RESULT THAT SHAPES EVERYTHING (read first)

⭐⭐ Two probed facts dominate the design, BOTH confirming the brief's CRUX and
ONE going further than the brief stated:

### 0a. The transverse face-slope is real but small (~12% of face-average, O(h))

Projecting the manufactured inflow `ψ_{n,g}(x_face, t)/W` onto the transverse
face Legendre basis `{1, ξ}` (bare per-transverse coefficient), face `xmin`
(transverse = y), inflow ordinates only:

```
nc=8 : |face_slope| max=1.91e-2 mean=8.72e-3; |face_bar| max=1.17e-1; slope/bar median=0.163 max=0.464
nc=16: |face_slope| max=1.24e-2 mean=4.59e-3; |face_bar| max=1.25e-1; slope/bar median=0.094 max=0.273
nc=32: |face_slope| max=6.64e-3 mean=2.31e-3; |face_bar| max=1.27e-1; slope/bar median=0.047 max=0.139
```

The slope/bar ratio HALVES per refinement (O(h)) — exactly the signature of an
intra-cell slope (it captures the cell-width variation of the BC trace, so it
is O(h) relative to the cell-average). The "~12% of the face-average" the brief
quoted is the nc=8/16 median; it is NOT a fixed fraction — it decays O(h).

### 0b. ⚠⚠ "IMPROVES ON FLAT" IS NOT ACHIEVABLE FOR LEG B (the sharper-than-Leg-A trap)

PROBED LIVE (faithful surrogate: monkeypatch `_inflow_to_moments` to seed slot 1
with the REAL projected manufactured transverse slope, full public solve, A-error
near the boundary):

```
nc=16 near-bdy A-err: flat(avg/centre)=1.0151e-2  real-slope=1.0303e-2  flipped=1.0091e-2  improves? NO
nc=32 near-bdy A-err: flat=1.9427e-3            real-slope=1.9658e-3  flipped=1.9294e-3   improves? NO
# Also tried slot0=AVERAGE (not centre) + slope: 1.1192e-2 vs avgbar-only 1.1055e-2 — STILL worse.
```

Seeding the TRUE transverse slope makes the converged near-boundary A-error
SLIGHTLY WORSE, and the FLIPPED slope is slightly BETTER — i.e. the boundary
correction's sign relative to A is below the discretization floor and its
contribution to the converged A-error is essentially noise. This is a SHARPER
Mode-10 instance than Leg A: for Leg A the bulk slope source carries real
sub-cell info across the WHOLE domain (improves-on-flat WAS achievable, probed
3.40e-3 < 5.99e-3); for Leg B the boundary signal is too LOCALIZED and too SMALL
(O(h) intra-cell at the edge only) → it sits below the bulk O(h²) A-error.

CONSEQUENCE: **the brief's Deliverable 5 ("improves on flat") must be DROPPED /
REFRAMED.** The positive verification is the STRUCTURAL lift-threading at machine
precision (the analog of Leg A's `array_equal` lift-pass-through), NOT a
converged-flux improvement. The teeth come from TWO O(1) structural places:
1. **The widened trace carries the projected transverse slope at machine
   precision** (Deliverable 1 — the production-change proof; a regression to
   re-zeroing the slope is caught at O(1) here).
2. **A CONSUMED transverse-face-slope sign flip moves the converged
   near-boundary flux ≫ `_CONSUMPTION_TOL`** (Deliverable 2 — the consumption
   proof; probed below), paired with the SCALAR-inflow no-op leg pinning the
   Leg-B asymmetry.

### 0c. The consumption signal IS detectable (the teeth are sound)

PROBED (real manufactured slope, +slope vs −slope flip, full public solve):

```
nc=16: real-vs-flip |Δφ|/|φ| near-bdy = 4.106e-3  global = 3.273e-3
nc=32: real-vs-flip |Δφ|/|φ| near-bdy = 8.383e-4  global = 8.990e-4
```

The flip moves the converged near-boundary flux 4.1e-3 (nc=16) — clears
`_CONSUMPTION_TOL=1e-8` by ~5.6 orders of magnitude. The signal halves under
refinement (O(h), boundary-localized), so the gate is a FIXED-MESH consumption
test, NOT a convergence leg. LINEARITY confirmed (seed slope = k·avg, k∈{.05,.1,.2}
→ flip |Δφ|/|φ| near-bdy = {2.4e-3, 4.8e-3, 9.5e-3}, exactly linear in k) — proves
the cochain GENUINELY consumes the transverse slope, not a numerical artifact.

### 0d. The no-op control is BYTE-IDENTICAL (the bit-identity guard holds)

PROBED: with the widened `_inflow_to_moments` seeding slot 1 = 0 (a SCALAR
inflow → transverse slopes stay zero), the converged flux is `array_equal` to
today's avg-only result (maxdiff 0.0). This is the Deliverable 4 negative control
by construction.

## 1. THE FACE-MOMENT NORMALIZATION CRUX (locked, probed — coordinate with method-implementer)

The production trace today carries one scalar per (ordinate, group,
transverse-CELL): `face_view(face)` shape `(N, ng, n_transverse)`, holding the
**cell-CENTRE** eval of the manufactured inflow `ψ_{n,g}(x_face, t_centre)/W`
(PROBED: face_view ≡ centre-eval, maxdiff 0.0; ONLY inflow ordinates nonzero).
`_inflow_to_moments` (`loss_representation.py:357-378`) widens each face to
`(N, ng, n_transverse, n_face_moments)` with `n_face_moments = per_axis^{d-1} = 2`
for d=2 LD, seeding ONLY `AVERAGE_MOMENT` (slot 0) ← scalar, slope ZERO.

The cochain consumes the per-face moments via `assemble_inflow_axis`
(`_ubld.py:270-323`): the upstream-face moment vector `(..., 2^{d-1})` is
weighted by the **transverse mass** `mass_1d(h_t, θ) = diag(h_t, θ·h_t)` (PROBED:
bar-only seed → `|μ|·h_t·bar` on the `[1,−1]` active trace; slope-only seed →
`|μ|·θh_t·slope`). So — EXACTLY like Leg A's cell mass `M=diag(h,θh)` —

⭐ **the production trace must carry the BARE per-transverse Legendre
coefficients** (NO θ, NO h_t, NO transverse-cell-length): the cochain's
transverse mass adds them. The two face moments, per (ordinate, group,
transverse-cell):

- **slot 0 = face-bar** = the transverse AVERAGE of `ψ_face/W` over the cell
  (`= ⟨ψ_face, P₀⟩ / ⟨P₀,P₀⟩`), the analog of Leg A's `q̄`.
- **slot 1 = transverse face-SLOPE** = `⟨ψ_face, P₁⟩ / ⟨P₁,P₁⟩` along the
  transverse coordinate ξ∈[−1,1] (for a linear `ψ_face = a + b·t`: `slope =
  b·h_t/2`). The analog of Leg A's `q̂_a`. THIS is what `_inflow_to_moments`
  zeros today.

⚠ **THE APPLES-TO-APPLES TRAP (the same crux that bit Leg A — flag for
method-implementer):** the projection's slot-1 normalization MUST match the
cochain's transverse-mass `θ`-convention. The test-side projector computes the
BARE coefficient (`⟨ψ,P₁⟩/⟨P₁,P₁⟩`, NO θ/h_t); the production trace MUST store
the BARE coefficient too (NOT pre-multiplied by `θh_t`). If production stores a
θ-weighted or h_t-weighted slope, the structural lift-threading gate (Deliverable
1) will fail by a constant `θh_t` factor — that is a TRUE bug (the cochain would
double-apply the mass), and the gate catching it is the design intent.

⚠ **THE SLOT-0 CENTRE-vs-AVERAGE CAVEAT (harder than Leg A stated, NEW for Leg
B):** today's scalar trace carries the cell-CENTRE value (PROBED), but the
PROJECTION's slot-0 is the cell-AVERAGE (they differ by O(h²)). The
method-implementer must decide whether the widened slot-0 stays CENTRE
(backward-compat: the scalar-inflow path must remain byte-identical — Deliverable
4) or switches to AVERAGE. RECOMMENDATION: keep slot-0 = the EXISTING scalar
trace (centre) for the scalar-inflow path (so DD/Step + the existing 1-D
prescribed-inflow MMS stay byte-identical), and have the MOMENT-resolved trace
carry (average-bar, slope) only when a moment-resolved inflow is explicitly
supplied. The Deliverable-1 structural gate must therefore compare the SLOPE slot
(slot 1) — NOT slot 0 — against the projected reference (slot 0 may legitimately
differ centre-vs-average). The face-projection foundation sub-gate (Deliverable
3) pins the projector's slot-0 IS the average separately.

L11 structural independence: the test-side projector uses ONLY
`numpy.polynomial.legendre.leggauss` + the closed-form manufactured `ψ_face` — it
NEVER calls `_inflow_to_moments` nor any `_ubld`/`LinearDiscontinuous` /
`assemble_inflow_axis` method. It descends only from `case._drivers` (the
manufactured A/B/C) + `numpy`.

## 2. WHAT THE PRODUCTION CHANGE LOOKS LIKE (for the gate to target)

A DIFFERENT production path than Leg A (which widened the BULK lift in
`solver.py`). Leg B widens the BOUNDARY trace:

- **`_LossRepresentation._inflow_to_moments`** (`loss_representation.py:357`):
  today seeds ONLY slot 0 (the average) and zeros the transverse slopes. The
  widening must seed slot 1 (the transverse face-slope) when a moment-resolved
  inflow is supplied. The single derivation site (shared by FFW + MFW, solve +
  apply directions — `_sweep_interior:1033` and `_loss_action_interior:1120`).
- **`BoundarySourceSink` / `TraceSpace`** (the trace must be able to CARRY the
  `2^{d-1}` per-face moment axis): today `face_view` is `(N, ng, n_transverse)`
  scalar-per-cell. A moment-resolved trace needs `(N, ng, n_transverse,
  2^{d-1})` OR a parallel moment-trace object. The method-implementer + explorer
  decide the TraceSpace widening shape; coordinate the face-moment Kronecker
  order (transverse axes, `[bar, slope]`) with `face_moment_tail`.
- **`case.prescribed_inflow`** (`mms/sn.py:1542`): today emits a per-transverse-
  CELL scalar (centre eval). A moment-resolved boundary test needs the case (or
  the TEST) to also emit the transverse slope. Per Leg A's Pattern-6 precedent
  (NO production projector — the test owns the projection, L11), RECOMMEND the
  TEST projects the manufactured inflow onto the face moments and supplies them;
  production accepts the moment-resolved trace, does NOT compute it.

⚠ The exact TraceSpace/trace-widening API is an EXPLORER task (the brief flags
"coordinate with the explorer's TraceSpace map"). This spec pins the
NORMALIZATION (§1) and the TEETH (§3-§5); the method-implementer + explorer pin
the trace-carrier shape. The gates below are written against the OBSERVABLE
contract: (a) the widened trace/cochain receives the projected slope at machine
precision; (b) a consumed slope flip moves the converged flux.

### ⭐ EXPLORER ALIGNMENT (`[[issue-251-trace-space-widening]]`, live @ d9396a2)

The explorer's deep file:line map CONFIRMS and CONCRETIZES this spec. Key
alignments (fold into the gate targeting):

1. **The trace storage lever is `geometry.py:boundary_face_layout` (1044-1050)** —
   append `face_moment_tail(per_axis**(ndim-1))` to each face slot shape. The
   `TraceSpace` itself is UNCHANGED (its metric + omega_dot_n are per-ordinate and
   broadcast over trailing axes — moment-ready by construction). NOT a separate
   trace object. DD/Step → `face_moment_tail(1) == ()` → byte-identical (the §6
   bit-identity guard holds at the layout level).
2. **`_inflow_to_moments` becomes a PASS-THROUGH (identity)** when the trace
   already carries the moment axis — the inflow tuple arrives `(*face, 2^{d-1})`-
   valued, NOT a scalar to zero-fill. ⭐ This is EXACTLY what Deliverable 1's
   threading gate tests: feed a moment-resolved face → assert slot-1 preserved
   (the pass-through contract). The gate is correctly targeted at the stable
   single producer site (`loss_representation.py:357`).
3. **The producer edit** `boundary_source_sink.py:274` (`view[inflow, :] =
   arr[inflow, :]` → `view[inflow] = arr[inflow]`, span ALL trailing axes) + the
   **4 outflow capture-collapse DROP sites** (`loss_representation.py:1068/1140/
   1319/1386`, guarded `if n_face_moments > 1`) are the rest of the production
   change. The MMS (test) projects the manufactured inflow onto the transverse
   face-Legendre basis and supplies the moment-resolved `face_values` — production
   just packs them (Pattern 6, no production projector, L11 — exactly Leg A).
4. **⚠ Leg B is 2-D-ONLY by construction** — a 1-D slab face is a point
   (`face_shape == ()` → `_n_face_moments = per_axis**0 = 1`, no transverse axis),
   so the 1-D prescribed-inflow MMS is byte-identical (the §6 D4 guarantee is
   structural, not coincidental).
5. **⚠ REFLECTIVE-COUPLING transverse-slope SIGN (a NEW verification concern,
   not covered by my vacuum-BC gates):** the explorer flags that the reflective
   `B` operator (`boundary_operator.py:_reflect_trace`) passes the moment axis
   through the ordinate permutation for STORAGE, but the transverse-slope SIGN
   under reflection across a face must be verified correct (the tangent-plane
   coordinate is preserved under a normal-flip reflection → PROBABLY no slope sign
   flip, but UNVERIFIED). The `SN2DCartesianLDStressMMSCase` is VACUUM-BC, so my
   gates do NOT exercise this. RECOMMEND a SEPARATE reflective-LD MMS leg (or an
   `op.H` adjoint check on the transverse-slope reflection) — flag for a
   numerics/qa review when production lands. This is a Leg-B follow-up, NOT in
   scope for the vacuum-BC gates here, but it is a genuine sign trap (Mode-1) the
   vacuum gates cannot see (H2: vacuum nulls the reflective coupling).

## 3. DELIVERABLE 1 (PRIMARY) — the STRUCTURAL teeth: trace threads the projected slope

The analog of Leg A's `test_ld_2d_external_slope_source_threaded_through_lift`
(`array_equal` lift pass-through). The production producer threads the projected
transverse face-slope through to the cochain at MACHINE PRECISION; the slot is no
longer zeroed.

**Two legs**, both foundation (`@pytest.mark.foundation`, NO `verifies` — a
software invariant on the producer threading, not a solver flux-shape claim):

1. **POSITIVE (the production-change proof):** build a moment-resolved boundary
   trace from the test-side projector (`_project_inflow_to_face_moments` →
   `leggauss` only), feed it through the widened `_inflow_to_moments` (or the
   widened trace constructor), and `np.testing.assert_array_equal` the SLOPE slot
   (slot 1) of the widened object equals the projected reference. (NOT slot 0 —
   the §1 centre-vs-average caveat.) A regression that re-zeroes the slope is
   caught at O(1) here (the converged flux is sub-floor sensitive — §0).

2. **NEGATIVE control (scalar-inflow no-op):** feed a SCALAR inflow (the existing
   `case.prescribed_inflow` — no transverse moment) through the widened
   `_inflow_to_moments` and assert the transverse slopes stay EXACTLY ZERO
   (`array_equal` to a zeros array). The scalar path must be byte-identical to
   today (PROBED: array_equal=True, maxdiff 0.0). THIS is the Leg-B asymmetry —
   the scalar default is correctly blind to the transverse slope.

⚠ Because the trace-carrier API is not yet fixed (§2), write Deliverable 1 to
target `_inflow_to_moments` directly (the single derivation site that IS the
producer of the widened face object): construct the input scalar/moment inflow,
call `rep._inflow_to_moments(inflow)`, inspect the returned tuple's slot 1. This
is the SHARPEST, most stable target (it does not depend on the eventual
TraceSpace shape). `xfail(strict=True)` until production widens it (today it
returns the average-only widening → slot 1 is zero → the positive leg fails →
xfail; once production lands the positive leg xpasses → flip to a plain
foundation gate per [[feedback_vv_tagging]]).

## 4. DELIVERABLE 2 — the consumption-proof mutation control (anti-pattern #11)

The PRIMARY sign-catcher (the consumption proof). Flip the transverse-face-slope
sign → the near-boundary converged flux moves ≫ `_CONSUMPTION_TOL=1e-8` (PROBED
4.1e-3 at nc=16); the SCALAR-inflow gate (no transverse moment) stays GREEN (the
Leg-B asymmetry — the no-op). Runtime monkeypatch reverted in `finally`.

### Mutation table

| # | Path | Mutation | NEW gate must | EXISTING scalar gate (`test_ld_2d_stress_converges_second_order`) |
|---|------|----------|---------------|-----------|
| B1 | BOUNDARY transverse face-slope | seed the moment-resolved trace's slot-1 with the REAL projected slope, then a SECOND solve with slot-1 NEGATED | converged NEAR-BOUNDARY `|Δφ|/|φ|` between +slope and −slope > `_CONSUMPTION_TOL` (PROBED 4.1e-3 ≫ 1e-8) | stays GREEN — it uses the SCALAR `prescribed_inflow` (slot-1 ≡ 0 → flipping zero is a no-op). THIS asymmetry IS the Leg-B Mode-10 gap being closed. |
| B2 | scalar no-op control | seed slot-1 = 0 (scalar inflow) | converged flux `array_equal` to today's avg-only solve (PROBED maxdiff 0.0) | n/a (this IS the today path) |

The teeth are the CONSUMED-flip converged-flux change (O(1) above the 1e-12 fixed
point, ~5.6 orders above `_CONSUMPTION_TOL`), NOT a convergence order and NOT a
value-band against A (the §0b trap — the A-error is sub-floor insensitive to the
slope sign). MEASURE near-boundary specifically (the edge cell mask): the signal
is boundary-localized (near-bdy 4.1e-3 > global 3.3e-3 at nc=16).

⚠ Because the trace-carrier API is not yet fixed, the mutation control is
written as a faithful SURROGATE: monkeypatch `_inflow_to_moments` to seed slot-1
from a test-supplied per-face slope dict (the REAL projected manufactured slope),
revert in `finally`. The +slope and −slope solves differ only in the slope sign.
`xfail(strict=True)` until production lands a moment-resolved trace (today the
surrogate IS the only way to exercise the path; once production accepts a
moment-resolved trace natively, RE-TARGET this gate onto the public API per
[[feedback_vv_tagging]] and drop the monkeypatch — the surrogate pins the
CONSUMER threading but is structurally blind to the PRODUCER/trace stamp, vv
Mode-11; the public-API re-target closes that).

NOTE the Mode-11 hazard explicitly: the monkeypatch surrogate recomputes the
production formula on BOTH sides, so it pins "the cochain consumes slot-1" (flip
it → flux moves) but is BLIND to whether the PRODUCTION trace stamp puts the
projected slope in slot-1 (a stamp that drops the slope would leave the surrogate
green because the surrogate carries the value itself). Deliverable 1 (the
structural lift-threading on the real producer) is the stamp catcher; Deliverable
2 (the surrogate consumption) is the consumer catcher. BOTH are required — name
them as the two halves per the Leg-A closeout's two-O(1)-teeth recipe.

## 5. DELIVERABLE 3 — the face-projection-correctness foundation sub-gate (L11)

A `@pytest.mark.foundation` test (NO `verifies` — a software invariant on the
face projector). Project a KNOWN 1-D polynomial `ψ_face(t) = c0 + c1·t` along the
transverse coordinate of a face cell `[tL, tR]` onto the face Legendre basis
`{1, ξ}` and assert the two coefficients EXACTLY (machine precision, `leggauss`
q_nodes≥2 integrates a linear exactly):

```
face_bar   = c0 + c1·tc                 (slot 0 — transverse cell AVERAGE)
face_slope = (h_t/2)·c1                  (slot 1 — bare transverse P₁ coeff)
```

where `tc = (tL+tR)/2`, `h_t = tR−tL`. The hand-derived coefficients are the
structurally-independent reference. The convention `(h_t/2)·∂_t` on `{1,ξ}`
MATCHES `mass_1d`'s `diag(h_t, θh_t)` Legendre basis — this is what makes the
projection apples-to-apples with the cochain's transverse mass (§1).

**SECOND leg:** project the manufactured inflow on a real mesh and assert slot-0
== a transverse cell-average computed by an INDEPENDENT fine quadrature
(q_nodes=12) to ~1e-8 (pins that the projector's bar IS the cell average, not the
cell-centre value — the §1 caveat: do NOT cross-check slot-0 against the existing
centre-eval `prescribed_inflow`, which differs by O(h²)).

**HOW IT AVOIDS TAUTOLOGY (L11):** the face projector uses only `leggauss` + the
hand-derived integral, NEVER `_inflow_to_moments` / `assemble_inflow_axis` / any
LD cell op; the reference is hand-laid 1-D polynomial algebra, not a production
echo. Branch 1 (`_2d_cartesian_ld_stress_symbolic`) emits only the per-ordinate
`Q_closed`, NOT a symbolic FACE moment — so the "numeric==symbolic face moment"
leg is NOT available (Branch 1 has no symbolic face-moment emitter; the polynomial
hand-check is the L11 ground). This sub-gate is GREEN NOW (no production change
needed — it tests only the test-side projector). It is the 1-D-transverse analog
of Leg A's `test_tensor_legendre_projection_matches_hand_polynomial` (which the
face projector can REUSE — the face projection is the per-axis factor of the
tensor projection).

## 6. DELIVERABLE 4 — the bit-identity guard (DriftWarning + scalar no-op)

The DEFAULT scalar-inflow path (transverse slopes zeroed) MUST stay
bit-identical after the trace widening. Three guards:

1. **The strict DriftWarning regression gate** (the canonical strict-clean gate):
   ```
   pytest tests/sn/sweep/core tests/sn/solve \
     -W "error::tests.sn.regression._regression_assert.DriftWarning"
   ```
   Baseline (Leg A closeout, live): 520 passed, 1 skipped, 4 xfailed. DD/Step
   (`n_face_moments == 1` → `_inflow_to_moments` is the identity) is byte-
   identical; the 2-D-LD scalar-inflow path (slot-1 still zeroed unless a
   moment-resolved trace is supplied) must not move.

2. **The existing 1-D prescribed-inflow MMS** `tests/sn/verification/analytical/
   test_mms_prescribed_inflow.py` (4 tests, slab + sphere — DD single-moment
   closures, `n_face_moments == 1` → identity widening) MUST stay GREEN
   (PROBED green now). The widening only triggers for `n_face_moments > 1` =
   2-D LD, so this file is byte-identical by construction.

3. **The existing LD two-paths gates** `test_ld_2d_two_paths_ffw_equals_mfw` +
   `test_ld_2d_stress_two_paths_ffw_equals_mfw` (the FFW≡MFW bit-id) MUST stay
   GREEN — both legs route through the SAME widened `_inflow_to_moments`, so the
   schedule-equivalence is preserved.

4. **NEW explicit scalar no-op leg** (the Deliverable-1 negative control,
   already counted in §3 leg 2): feeding a scalar inflow through the widened
   `_inflow_to_moments` gives slot-1 ≡ 0 (`array_equal`), and the full solve is
   `array_equal` to today's (PROBED maxdiff 0.0).

## 7. DELIVERABLE 5 — ⚠ DROPPED/REFRAMED (the near-boundary improvement leg)

**The brief's Deliverable 5 ("the moment-resolved boundary trace lands the
near-boundary flux strictly closer to A than the scalar trace") is NOT
ACHIEVABLE** — PROBED (§0b): seeding the real transverse slope makes the
converged near-boundary A-error SLIGHTLY WORSE in every config tested
(centre-bar, average-bar, average-bar+slope), and the FLIPPED slope is slightly
BETTER. The boundary correction is sub-floor: the converged A-error is dominated
by the bulk O(h²) discretization, and the localized O(h)-small boundary-trace
slope sits below it.

REFRAME: the POSITIVE verification of Leg B is the STRUCTURAL lift-threading
(Deliverable 1 leg 1 — the trace carries the projected slope at machine
precision), NOT a converged-flux improvement. This is the §0b lesson: Leg B is a
SHARPER Mode-10 instance than Leg A (Leg A's bulk slope DID improve-on-flat;
Leg B's boundary slope cannot). Document this in the honest-scope note: the Leg-B
transverse face-slope is now CONSUMED (Deliverable 2 reddens on a flip) and
THREADED at machine precision (Deliverable 1), but its converged-flux
contribution is below the discretization floor (so there is NO converged-value
improvement leg — the teeth are structural). This is the canonical resolution of
a Mode-10 gap where the activated term's converged-value contribution is
genuinely sub-floor.

⚠ SELF-IMPROVEMENT FLAG (for delivery): this is a NEW SHAPE within Mode-10 — an
activated-and-consumed term whose converged-value contribution is sub-floor not
just for the SIGN (the canonical Mode-10) but for ANY value claim, so the
"add a companion gate that isolates the term so its error is O(1)" half of the
Mode-10 recipe FAILS (no fixed-source problem makes a boundary-trace slope the
DOMINANT forcing — it is intrinsically a boundary perturbation). The resolution
is STRUCTURAL teeth ONLY (producer-threading at machine precision + consumed-flip
≫ tol), with NO value-improvement leg. The existing vv Mode-10 row already
covers "activated-but-unconstrained → mutation-verify + O(1)-isolating
companion"; THIS instance shows the companion may be UNAVAILABLE (the term has no
O(1)-dominant regime), in which case the structural producer-threading +
consumed-flip teeth are the complete resolution. RECOMMEND a one-line addition to
the vv Mode-10 row at delivery: "if no regime makes the term the dominant forcing
(e.g. a boundary-trace slope), the producer-threading-at-machine-precision +
consumed-flip-≫-tol structural pair is the complete resolution — there is NO
value-improvement leg to add." NO new Mode number (it is a sub-case of Mode-10).
NO ERR until a real production bug (next free ERR-063).

## 8. DELIVERABLE 6 — the negative pin (reject a wrong transverse-moment width)

The trace widening must STILL reject a genuinely-wrong moment width. By analogy
to Leg A's `test_moment_resolved_bulk_still_rejects_wrong_trailing_axis`
(`per_axis**ndim = 4` for the bulk), the BOUNDARY moment width is
`per_axis**(ndim-1) = 2^{d-1} = 2` for d=2 LD. The widened trace constructor /
`_inflow_to_moments` must reject a moment-resolved inflow whose trailing
transverse-moment axis ≠ `2^{d-1}` (a clear ValueError naming the expected
`2^{d-1}`), AND DD/Step (`n_face_moments == 1`, no transverse-moment axis) must
reject any moment-resolved inflow outright.

⚠ The exact rejection site depends on the trace-carrier API (§2). Write the
negative pin against the eventual moment-resolved-trace constructor (coordinate
with explorer). `xfail(strict=True)` until production lands the moment trace.
If the eventual API does not have a single shape-validating entry, the negative
pin lives on `_inflow_to_moments` (assert it raises on a `(..., 3)` transverse
moment width on a 2-D LD mesh, and on any moment width on a DD mesh).

## 9. -O safety (Mode 8) + file location

The canonical invocation is `python -O` (CONFIRMED live: pytest warns "assert
statements are not executed"). ALL new gates MUST use `np.testing.*` /
`pytest.fail` / `pytest.raises` only — NO bare `assert` (false green under -O).
Mirror the existing `test_mms_ld_2d.py` style (every gate there is -O-safe).

**File**: REPLACE the `test_ld_2d_boundary_transverse_face_slope` skip stub at
the END of `tests/sn/verification/mms/test_mms_ld_2d.py` (lines ~1022-1037) with
the real gates. The face projector can REUSE the existing
`_project_scalar_to_tensor_legendre` helper's per-axis factor (a face is the 1-D
transverse projection — the tensor projector's per-axis Legendre coefficient).
`volume_weighted_l2` does NOT exist on this branch — use the local `_l2_2d`
(line 60).

## 10. THE GATES (xfail-strict until production lands the moment trace)

In `tests/sn/verification/mms/test_mms_ld_2d.py`, replacing the #251 skip stub:

| Gate | Mark | State | What it proves |
|------|------|-------|----------------|
| `test_face_transverse_legendre_projection_matches_hand_polynomial` | foundation | GREEN NOW | D3 leg 1 — the face projector reproduces the hand-derived `[c0+c1·tc, (h_t/2)·c1]` to machine precision (L11) |
| `test_face_projection_slot0_is_transverse_cell_average` | foundation | GREEN NOW | D3 leg 2 — slot-0 == fine-quad transverse cell average (not centre) |
| `test_ld_2d_boundary_slope_threaded_through_inflow_to_moments` | foundation, xfail(strict=True, reason=#251 prod) | xfail → flip on prod | D1 leg 1 — `_inflow_to_moments` threads the projected slope into slot-1 at machine precision; + leg 2 scalar no-op (slot-1 ≡ 0) |
| `test_ld_2d_boundary_slope_sign_mutation_reddens` | l1, verifies("ld-cartesian-2d"), xfail(strict=True, reason=#251 prod) | xfail → flip on prod | D2 B1 — consumed transverse-slope flip moves near-bdy converged flux ≫ `_CONSUMPTION_TOL`; + the scalar-inflow no-op leg (B2) pinning the Leg-B asymmetry |
| `test_ld_2d_boundary_scalar_inflow_byte_identical` | foundation, xfail(strict=True, reason=#251 prod) | xfail → flip on prod | D4 — scalar inflow through the widened path == today's solve (`array_equal`) |
| `test_moment_resolved_trace_rejects_wrong_transverse_width` | (no level), xfail(strict=True, reason=#251 prod) | xfail → flip on prod | D6 — reject transverse-moment width ≠ `2^{d-1}`; DD rejects any moment trace |

xfail rationale ([[feedback_vv_tagging]]): the xfail `reason=` names the unlocking
production change ("the moment-resolved boundary trace / widened
`_inflow_to_moments` is not yet implemented — #251"). `strict=True` so the gate
FLIPS to a hard fail (alerting "production landed, un-xfail me") the moment the
slope is threaded. The two foundation projection sub-gates are NOT xfail (they
test only the test-side projector — GREEN now).

⚠ The 4 xfail gates depend on the trace-carrier API the explorer +
method-implementer fix. They are written against `_inflow_to_moments` (the stable
single producer site) via a faithful surrogate where the public API is not yet
moment-aware. Once production lands a moment-resolved trace natively, RE-TARGET
the consumption + threading gates onto the public `solve_sn_fixed_source`
(dropping the monkeypatch surrogate — it pins the consumer but is Mode-11-blind
to the producer stamp; the public re-target closes that) and flip xfail→plain.

## 11. Baseline (live @ HEAD d9396a2, branch refactor/sn-foundation-cleanup)

```
pytest tests/sn/verification/mms/test_mms_ld_2d.py -q -m "not slow"
  → 14 passed, 1 skipped, 3 deselected  (70s)   [1 skip = the #251 stub I replace]
pytest tests/sn/verification/analytical/test_mms_prescribed_inflow.py -q
  → 4 passed  (0.5s)   [the bit-identity-by-construction baseline for D4]
```
`per_axis = 2`, `ndim = 2`, `n_face_moments = per_axis^(d-1) = 2`,
`AVERAGE_MOMENT = 0`, transverse face-moment order `[bar=0, slope=1]`. NEVER run
all `tests/sn` (#212 hang at `test_keff_slab::test_heterogeneous_absolute_keff`);
scope to `tests/sn/verification/mms` + `tests/sn/verification/analytical` +
`tests/derivations`.

## 12. WHERE THE FACE-NORMALIZATION CRUX IS HARDER THAN THE BRIEF STATED

1. **Slot-0 centre-vs-average (NEW, §1):** the brief said "match the production
   face-moment normalization". The harder reality: today's scalar trace carries
   the cell-CENTRE, not the average, so the widened slot-0 has a backward-compat
   decision (keep centre for the scalar path = byte-identical; the moment path's
   slot-0 = average). The Deliverable-1 structural gate must compare SLOT-1 only
   (slot-0 may legitimately differ centre-vs-average). The face-projection
   sub-gate pins slot-0 == average SEPARATELY (against fine quad, NOT against the
   centre-eval `prescribed_inflow`).

2. **The improves-on-flat leg is impossible (§0b, §7):** the brief's Deliverable
   5 assumed the boundary slope improves the converged near-boundary flux toward
   A. PROBED FALSE — the boundary signal is sub-floor (worse, not better, and the
   flip is slightly better). This is a SHARPER Mode-10 than Leg A; Deliverable 5
   is DROPPED and the positive verification is the structural lift-threading.

3. **The trace-carrier API is not yet fixed (§2):** Leg A widened a single
   `solver.py` lift on an EXISTING moment-resolved bulk path (a free ride). Leg B
   needs a NEW moment-resolved TRACE carrier (TraceSpace widening) — an explorer
   task. The gates are written against the stable `_inflow_to_moments` producer
   site (via surrogate) so they do not depend on the eventual trace shape; they
   re-target onto the public API once the trace lands.

## 13. Self-improvement (Mode-10 sub-case flagged)

This is a THIRD Mode-10 instance (after #240 D5b-S4 and #247 Leg A) and the FIRST
where the converged-value contribution is sub-floor for ANY value claim (not just
the sign) → the "O(1)-isolating companion" half of the Mode-10 recipe is
UNAVAILABLE (a boundary-trace slope has no regime where it is the dominant
forcing). The resolution is STRUCTURAL teeth ONLY. Recommend a one-line addition
to the existing vv Mode-10 row at delivery (NO new mode, NO skill rewrite — the
row already names the recipe; this sharpens the "when the companion is
unavailable" branch). NO ERR until a caught production bug (next free ERR-063).
