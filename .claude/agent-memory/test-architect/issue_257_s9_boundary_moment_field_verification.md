---
name: issue-257-s9-boundary-moment-field-verification
description: "#257 S9 PRE-IMPL gate spec — BoundaryMomentField (typed moment-tail boundary trace) + close the moment-state boundary drop in 2-D Cartesian LD SN. ⭐⭐ CONVERGENCE-REGIME VERDICT: the user's coherent-promise improves-on-flat-AT-THE-BOUNDARY value gate is NOT ACHIEVABLE — PROBED across global/edge/first-row regions + 3 transverse frequencies, flat≡mom≡flip all O(h²), boundary slope is sub-floor (REFUTES the brief's value-gate ask, confirms #251 wall is FUNDAMENTAL not regime-specific). Teeth = STRUCTURAL only (type-intrinsic + threading + consumed-flip), which #251 Leg B ALREADY landed live. Branch feature/field-typed-operator-algebra, HEAD 8e0a2cf. PRE-IMPL."
metadata:
  type: project
---

# #257 S9 — `BoundaryMomentField` + close the moment-state boundary drop (gate spec)

PRE-IMPL test-architect spec. Branch `feature/field-typed-operator-algebra`,
HEAD `8e0a2cf`. Host `.venv/bin/python -O`. I design the gate; the production
(`BoundaryMomentField` type + the MMS `prescribed_inflow` moment-slot build +
the moment-state closure) is owned by method-implementer. Extends
[[issue-251-legB-boundary-gate-spec]] (the LIVE #251 Leg B boundary gates —
the structural-teeth precedent S9 INHERITS) + [[issue-247-moment-source-gate-spec]]
(the bulk improves-on-flat precedent) + [[issue-257-s5-functional-category-verification]]
(the S-series type-intrinsic-gate discipline) + [[feedback-test-intrinsic-properties]]
+ [[feedback_vv_tagging]] + [[feedback_regression_tolerance_design]].

## 0. THE STATE RECONCILIATION (read FIRST — the brief is half-stale)

⭐⭐ The brief describes a gap that is LARGELY ALREADY CLOSED. `e5f2b1c`
(#251 Leg B) is an ANCESTOR of HEAD (the foundation-cleanup cluster merged to
`main` @ `05fa1ef`, this branch is off it). So the #251 Leg B PRODUCTION
change is LIVE on this branch:

- `_inflow_to_moments` (`loss_representation.py:448-498`) is the PASS-THROUGH
  for a moment-resolved inflow (threads slot-1 unchanged; rank-discriminated).
- `SNMesh.boundary_face_layout` appends `face_moment_tail(per_axis^d)` (the
  moment tail is in `layout.total_size` — the flat boundary buffer ALREADY
  holds it).
- `BoundarySourceSink.prescribed_inflow` (`boundary_source_sink.py:284-294`)
  accepts a FULL moment slot (`view.shape`) OR a scalar (`view.shape[:-1]` →
  seeds AVERAGE_MOMENT, slope zero — the Leg-B asymmetry).
- The capture/outflow STORES the moment whole (`loss_representation.py:1199`
  comment: "no longer collapsed to the average").
- **All 6 #251 Leg B boundary gates PASS** (no longer xfail — PROBED live: 6
  passed). `_solve_with_boundary_slope` (`test_mms_ld_2d.py:1211`) drives the
  REAL public moment-resolved path END-TO-END (producer stamp + cochain
  consumer); the consumed-flip mutation + the scalar no-op control are GREEN.

So S9 is NOT re-deriving the #251 boundary value/structural gates. S9's ACTUAL
deliverables (per #257 plan §S9 + #256 fork #1):
1. **The typed `BoundaryMomentField`** — today the moment trace rides as an
   untyped trailing axis on `BoundarySourceSink`/`BoundaryFlux`'s flat buffer
   (the buffer FITS it via the layout tail, but no type NAMES "this trace is
   moment-resolved"). The gap is a TYPE, not storage.
2. **The MMS case builds the moment slot directly.**
   `SN2DCartesianLDStressMMSCase.prescribed_inflow` (`mms/sn.py:1542-1598`)
   STILL builds a SCALAR per-face `vals=(N,ng,n_t)` (cell-CENTRE eval) → hits
   the scalar branch → slope dropped. Only the TEST helper
   `_face_transverse_buffers` (`test_mms_ld_2d.py:1161`) projects the slope.
   S9 moves that projection into the production MMS path (single-sourced with
   the bulk #247 projector).
3. **Close the moment-state boundary drop.** Disambiguate: the `(moment_buf,
   None)` at `loss_representation.py:3395` is the HARMONIC-ANGULAR moment
   windowing output (Phase 5c) — the `None` is the boundary trace slot of the
   moment-PROJECTION mode (NOT yet returned). This is a DIFFERENT axis from the
   #251 TRANSVERSE-SPATIAL face-moment. #256 fork #1 scopes S9 to the
   TRANSVERSE-SPATIAL moment (the #251 lever). Coordinate with explorer which
   `None` the user means — likely the boundary trace being a TYPED
   moment-resolved field, not the angular-moment windowing `None`.

## 1. ⭐⭐ THE CONVERGENCE-REGIME VERDICT (the crux — REFUTES the brief)

**The user's coherent-promise value gate — "moment-resolved boundary trace
lands the near-boundary flux strictly closer to A / recovers clean 2nd order
AT THE BOUNDARY where flat does not" — is NOT ACHIEVABLE.** Not for the
existing MMS, and not for any non-trivial-boundary-slope regime I could
construct. This CONFIRMS #251's wall as FUNDAMENTAL, and goes one step further
(I tested the user's exact "make the boundary slope non-trivial" hypothesis and
it still fails). The vv-honesty rule (L11 / Mode-10) requires I report this, not
manufacture a passing gate.

### The probes (all live @ HEAD 8e0a2cf, host -O)

**Probe A — improves-on-flat near the boundary, existing MMS:** seeding the
REAL projected transverse slope makes the converged near-bdy A-error SLIGHTLY
WORSE; the flipped slope is slightly BETTER (sub-floor noise):
```
nc=16 near-bdy A-err: flat=9.40e-3 mom=9.54e-3 flip=9.34e-3  improves? NO
nc=32 near-bdy A-err: flat=1.80e-3 mom=1.82e-3 flip=1.79e-3  improves? NO
```

**Probe B — the user's exact hypothesis (sharp transverse boundary
variation):** added a K× high-frequency transverse term so the boundary trace
oscillates sharply. STILL no improvement (K=3 and K=6):
```
K=3 nc=16: flat=1.01e-2 mom=1.01e-2 improves? NO (ratio 0.997)
K=6 nc=16: flat=3.96e-2 mom=3.97e-2 improves? NO (ratio 0.999)
```

**Probe C — convergence ORDER across regions (global / edge-ring /
first-row), flat vs mom vs flip, nc=16/32/64:** the boundary treatment changes
NOTHING above the floor — every region O(h²) or better, orders identical:
```
global   : flat 2.001/2.005  mom 1.998/2.003  flip 1.999/2.004
edge-ring: flat 2.385/2.461  mom 2.389/2.462  flip 2.386/2.460
first-row: flat 2.499/2.510  mom 2.507/2.516  flip 2.499/2.508
```
The near-boundary flux is even CLEANER than O(h²) (the prescribed inflow is
exact at the face) — there is no order to recover.

### The MECHANISM (why — the load-bearing physics)

**Probe D — face REPRESENTATION error (bar-only vs bar+slope inflow trace):**
THIS is above-floor — the boundary slope genuinely lifts the inflow
representation O(h)→O(h²):
```
nc=16: bar-only=4.04e-3  bar+slope=4.69e-4  ratio= 8.6
nc=32: bar-only=2.04e-3  bar+slope=1.18e-4  ratio=17.2 (bar-only halves=O(h); bar+slope quarters=O(h²))
nc=64: bar-only=1.02e-3  bar+slope=2.96e-5  ratio=34.5
```

So the slope improves the *inflow representation* by a full order — but that
improvement does NOT propagate above-floor into the converged *flux*. The
reason: the bar-only face inflow's O(h) error is a BOUNDARY-LOCALIZED forcing.
Integrated over the codimension-1 boundary into the 2-D-domain L2 norm it
contributes O(h^1.5), which is SUB-DOMINANT to the bulk O(h²) once the
LD interior closure (already O(h²)) carries it inward. The boundary
transverse-slope moment is a SECOND-ORDER correction to an ALREADY-second-order
face representation — it cannot move the converged flux above the bulk floor.

**This is the canonical Mode-10 resolution where the companion-isolating gate
is UNAVAILABLE** (a boundary-trace slope has no regime where it is the dominant
forcing — the boundary is measure-zero in the refinement limit). Exactly the
#251 §0b/§7 finding, now RE-CONFIRMED on this branch AND against the user's
non-trivial-slope hypothesis. **NO value gate. Structural teeth only.**

### The honest framing of the coherent promise (for the user + the docstring)

The coherent promise "LD gives 2nd order EVERYWHERE including the boundary" is
TRUE and ALREADY DELIVERED — but NOT because the transverse boundary-slope
moment is what delivers it. The LD face closure achieves O(h²) at the boundary
from the AVERAGE moment alone (the prescribed inflow is exact at the face
cells; the bulk LD closure carries it inward at O(h²)). The transverse
boundary-slope moment is a representation refinement (O(h)→O(h²) on the
inflow trace) that is genuinely consumed (a flip moves the flux 4.1e-3
near-bdy ≫ tol) and genuinely threaded (machine-precision), but its
converged-flux contribution is sub-floor. So S9 does NOT "remove a
not-at-the-boundary asterisk" on the convergence ORDER (there was none — it is
already O(h²) at the boundary); S9 closes the TYPE/representation gap and the
consumption/threading is the verifiable content. **Recommend the user reframe
the S9 motivation from "recover 2nd order at the boundary" (already true) to
"type the moment-resolved boundary trace + thread the inflow-representation
refinement that the LD closure now consumes" (the real, verifiable content).**

## 2. THE GATE LIST (what S9 SHOULD ship)

Most positive verification of the boundary CONSUMPTION/threading is ALREADY
LIVE (the 6 #251 gates). S9 ADDS: (a) the type-intrinsic gate for
`BoundaryMomentField`; (b) re-target the #251 surrogate-built boundary onto the
production MMS `prescribed_inflow` once it builds the moment slot; (c) the
single-source projection gate; (d) the bulk-regression bit-id guard. NO value
gate (§1).

### GATE A — `BoundaryMomentField` type-intrinsic (anti-pattern #11, foundation)

The S-series intrinsic-property gate (the user directive
[[feedback-test-intrinsic-properties]]; mirrors S1 cone / S5 Functional≠Operator).
`BoundaryMomentField` IS-A `BoundaryField`; its DEFINING property is the
moment-axis-matching `_check_partner` + the moment-resolved layout. Both
directions + discriminator foils, `-O`-safe (`np.testing`/`pytest.raises`):
- POSITIVE: a `BoundaryMomentField` on a moment-resolved trace (`per_axis^d>1`)
  constructs; its arithmetic with a sibling `BoundaryMomentField` closes
  (add/sub/scalar-mul preserve the moment axis); IS a `BoundaryField`
  (`isinstance`), IS a `Vector`, IS a `FullField`-compatible boundary leaf.
- NEGATIVE (the discriminator): `BoundaryMomentField + BoundaryFlux` (or a
  scalar `BoundarySourceSink`) → raises (the `_check_partner` cross-class gate
  — a moment-resolved trace must not silently mix with a scalar trace). A
  `BoundaryMomentField` on a DD/Step (`face_moment_count==1`, no moment axis)
  layout → reject (or be the degenerate scalar — pin which by construction).
- DISCRIMINATOR foil: a same-shape `BoundaryFlux` that happens to have the same
  `total_size` is NOT a `BoundaryMomentField` (class identity, not shape) — the
  type is the moment-ROLE marker, not the buffer width.
⚠ Mode-11: if `BoundaryMomentField` is a thin leaf (storage inherited, like
`BoundaryFlux`), its `_check_partner` is inherited too — VERIFY the
moment-discrimination teeth actually fire (a wrong-class add MUST raise). If the
type adds NO discriminating behaviour beyond class identity, say so (it is then
a naming leaf — still legitimate per #256, but the intrinsic gate is just the
class-identity guard + the Vector/BoundaryField conformance).

### GATE B — re-target the boundary threading + consumption onto the production MMS

The #251 gates `_solve_with_boundary_slope` build the moment slot in the TEST
(`_face_transverse_buffers`) and feed it via the public
`BoundarySourceSink.prescribed_inflow`. When S9 makes
`SN2DCartesianLDStressMMSCase.prescribed_inflow` emit the moment slot DIRECTLY,
RE-TARGET:
- `test_ld_2d_boundary_slope_threaded_through_inflow_to_moments` (foundation):
  KEEP — it already pins `_inflow_to_moments` threads slot-1 at machine
  precision against the leggauss-only `_face_transverse_buffers` reference
  (L11). Add a leg: the PRODUCTION `case.prescribed_inflow(sn)` now produces a
  moment slot whose slot-1 == the leggauss reference (the producer-stamp
  catcher — closes the Mode-11 producer-blindness the surrogate had).
- `test_ld_2d_boundary_slope_sign_mutation_reddens` (l1, verifies
  ld-cartesian-2d): KEEP — the consumed-flip ≫ `_CONSUMPTION_TOL` teeth (PROBED
  4.1e-3 near-bdy at nc=16, linear in slope magnitude). When the MMS builds the
  moment slot, the flip can be applied to the production-emitted slot (not just
  the test surrogate).
- `test_ld_2d_boundary_scalar_inflow_no_op_negative_control` (foundation):
  KEEP — a scalar inflow → slope zero → byte-identical (the Leg-B asymmetry).

### GATE C — single-source projection (Cardinal Rule 2)

The MMS boundary projection MUST be the SAME primitive as the bulk #247
projection. Today: bulk = `_project_scalar_to_tensor_legendre`
(`test_mms_ld_2d.py:509`, the full tensor projector); face = `_face_transverse_legendre`
(`:1050`, the per-AXIS factor). The face projection IS the per-axis factor of
the tensor projection (the #251 spec §5 notes this). S9's production MMS
boundary build must call the SAME shared face-Legendre projector the bulk uses
(NOT a third copy). These projectors are deliberately TEST-side (L11) — the
production MMS case may grow a `_project_inflow_to_face_moments` that descends
ONLY from `case._drivers` + numpy `leggauss` (NEVER `_inflow_to_moments` /
`_ubld` / any LD op). Gate (foundation): the production MMS face projector ==
the leggauss hand-derived `[c0+c1·tc, (h_t/2)·c1]` on a known polynomial
(machine precision — the existing
`test_face_transverse_legendre_projection_matches_hand_polynomial` covers the
test-side projector; add a leg pinning the production MMS projector agrees).
⚠ FACE-NORMALIZATION CRUX (locked, from #251 §1): the trace carries the BARE
per-transverse Legendre coeff `[bar=⟨ψ,P₀⟩/⟨P₀,P₀⟩, slope=⟨ψ,P₁⟩/⟨P₁,P₁⟩]` —
NO θ/h_t (the cochain's transverse mass `mass_1d(h_t,θ)=diag(h_t,θh_t)` adds
them). A θ- or h_t-weighted slope fails Gate B's threading by a constant factor
(a TRUE bug — double-applied mass). SLOT-0 centre-vs-average caveat: today's
scalar trace carries cell-CENTRE; the projection's slot-0 is cell-AVERAGE
(differ O(h²)). The threading gate compares SLOT-1 only.

### GATE D — bulk regression bit-identity (S9 must not perturb anything)

The DD/Step (`per_axis==1` → `face_moment_count==1` → no moment axis) +
the bulk moment tensor MUST stay byte-identical. The canonical strict gate:
```
.venv/bin/python -O -m pytest tests/sn/sweep/core tests/sn/solve \
  -W "error::tests.sn.regression._regression_assert.DriftWarning"
```
Plus: the 6 #251 boundary gates stay GREEN; the bulk #247 gates
(`test_ld_2d_external_slope_source_*`) stay GREEN; the 1-D prescribed-inflow MMS
(`test_mms_prescribed_inflow.py`, 4 tests, `face_moment_count==1` → identity
widening) byte-identical; the keff/MMS regression subset (the plan's subset).
A new TYPE (`BoundaryMomentField`) that does not change the buffer or the
arithmetic is bit-identical by construction — the gate confirms it.

### GATE E — the moment-state boundary-drop closure (IF S9 closes the `None`)

⚠ ONLY if explorer confirms S9 closes the `(moment_buf, None)` angular-moment
windowing boundary slot (NOT just the transverse-spatial type). If so: the
new boundary block in the moment-projection-mode return is provably `==` the old
`None`'s reconstruction (the `solve_moments` resolvent currently reconstructs
the boundary from the scalar path — gate that the typed boundary block byte-
matches it). This is the #256-plan "new boundary block provably == old None"
gate. Coordinate the EXACT seam with explorer — the brief conflates the
transverse-spatial moment (#251 lever, fork #1) with the angular-moment
windowing `None`; pin which before writing this gate.

## 3. MMS CASE DESIGN (the brief's Deliverable 1)

VERDICT: **EXTEND `SN2DCartesianLDStressMMSCase` (no sibling case).** Reasons:
- The existing case ALREADY activates the boundary closure at the FACE-AVERAGE
  moment (`a0>0` non-vanishing inflow on all 4 edges) and the transverse face
  slope IS present (PROBED |slope|max~1.2e-2, O(h)). A sibling case would
  duplicate the driver machinery (Cardinal Rule 2).
- The ONLY change the case needs is `prescribed_inflow` building the moment slot
  (slot-1 = projected transverse slope) instead of the scalar — a PRODUCTION
  COMPLETION of the existing case, not a new ansatz.
- Probes A/B showed a sharper transverse boundary frequency does NOT unlock a
  value gate — so there is NO motivation for a new high-frequency sibling. The
  existing ansatz's boundary slope is non-trivial ENOUGH for the structural
  teeth (the consumed-flip moves the flux 4.1e-3 ≫ tol).

### Mode-7 ACTIVATES / NULLS declaration (for the boundary-slope term)
- ACTIVATES: the bilinear UBLD per-axis SPATIAL slope rows (the bulk, #247);
  the FACE-AVERAGE boundary closure (`a0>0`, #251 already); the TRANSVERSE
  face-SLOPE boundary moment (slot-1, the S9 target — ACTIVATED: PROBED the
  consumed-flip moves the converged flux 4.1e-3, linear in slope magnitude →
  GENUINELY consumed, not carried-and-ignored).
- NULLS: nothing the boundary slope needs. The vacuum-mesh-BC + prescribed-
  inflow-via-`q.boundary` posture NULLS the REFLECTIVE coupling (H2) — so the
  reflective-BC transverse-slope sign trap (#252, the #251 §2.5 explorer flag)
  is NOT exercised. Declare it: S9's vacuum-BC MMS is BLIND to the reflective
  transverse-slope sign (a Leg-B follow-up, #252 filed, NOT in S9 scope).
- The boundary-slope term is ACTIVATED but per §1 its converged-value
  contribution is SUB-FLOOR — Mode-10 with the companion-isolating gate
  UNAVAILABLE. CONSTRAINED by the consumed-flip mutation (Gate B), NOT a value
  band. ≥2G (the case is 2G asymmetric downscatter, L2); heterogeneous
  (`_default_hetero_2d_xs`, Σ_a>0); NON-SQUARE (Lx≠Ly, x↔y-swap defence).

## 4. METHOD-IMPLEMENTER PRE-READ + ANTI-RECOMMENDATIONS

Pre-read (file:line):
- `orpheus/transport/fields/boundary_flux.py:84-127` — the thin-leaf template
  (`BoundaryMomentField` is its sibling: storage/algebra/face-access inherited
  from `BoundaryField`, leaf adds UNITS + `_DISPLACEMENT_CLS` + class identity).
- `orpheus/transport/fields/_bases.py:480-562` — `BoundaryField` base (flat
  buffer keyed on `layout.total_size` — ALREADY holds the moment tail;
  `_check_partner` mesh+layout guard; `face_view` slice).
- `orpheus/transport/source_sinks/boundary_source_sink.py:188-295` — the
  moment-aware `prescribed_inflow` setter (scalar vs full-slot discrimination
  LIVE).
- `orpheus/sn/loss_representation.py:448-498` (`_inflow_to_moments` pass-through)
  + `:1170-1200` (capture stores the moment whole) + `:3390-3429` (the
  `(moment_buf, None)` angular-moment windowing — the OTHER moment axis, do not
  conflate).
- `orpheus/numerics/moment_layout.py:64-78` (`face_moment_count` =
  `per_axis^(d-1)` — the single-source width the type + projection must agree on).
- `orpheus/derivations/continuous/mms/sn.py:1542-1598`
  (`prescribed_inflow` — STILL builds scalar; the production completion target)
  + `:1387-1432` (`_drivers` — the projection source).
- `tests/sn/verification/mms/test_mms_ld_2d.py:1050-1450` (the LIVE #251 gates +
  `_face_transverse_buffers`/`_face_transverse_legendre`/`_solve_with_boundary_slope`
  — re-target onto production, don't rebuild) + `:509` (`_project_scalar_to_tensor_legendre`
  — the bulk projector to single-source the face projection from).

Anti-recommendations (≥5, what NOT to build):
1. **Do NOT retype `BoundaryFlux`** — the per-ordinate reflective trace `B`
   needs `BoundaryFlux` as-is (`boundary_operator.py:249-262`).
   `BoundaryMomentField` is a SIBLING leaf, not a replacement.
2. **Do NOT mirror `BulkField._compose_spatial_moments`** — boundary storage is
   a FLAT `FaceLayout` buffer (the moment tail is in `layout.total_size`), NOT a
   space-factor like the bulk `HarmonicMomentField`. The moment axis is a
   trailing tail on the flat buffer, not a separate space dimension.
3. **Do NOT duplicate the SymPy/leggauss projection** — the production MMS face
   projection MUST be the per-axis factor of the bulk
   `_project_scalar_to_tensor_legendre` (Cardinal Rule 2 / Gate C). These
   projectors are TEST-side (L11) — production accepts the moment slot, does NOT
   compute it (Pattern 6, exactly Leg A/B).
4. **The moment carrier stays the ONE `TimedFullField`** — `MomentFullField` is
   a FACTORY, not a subclass (#256 plan). Do NOT add a `TimedFullField`
   subclass for the moment-resolved boundary.
5. **Do NOT manufacture an improves-on-flat / order-recovery value gate** —
   §1 PROVES it is sub-floor (would falsely RED a correctly-consumed term, or
   pass only by tuning, neither honest). The teeth are STRUCTURAL.
6. **Do NOT conflate the two moment axes** — the #251 transverse-spatial
   face-slope (`face_moment_count`, the S9 fork-#1 target) is ORTHOGONAL to the
   Phase-5c harmonic-angular `moment_buf` windowing `None`
   (`loss_representation.py:3395`). Resolve which the user means with explorer
   before touching the `None`.
7. **Do NOT add a per-cell simplex/cone invariant to the boundary moment** —
   the boundary moment is a Legendre coefficient `[bar, slope]`, NOT a
   probability/coefficient field. It is a `BoundaryField` (flux-trace) leaf, not
   a `CoefficientField`.

## 5. MODE DISCIPLINE
- **Mode-8** (`-O`): all gates `np.testing`/`pytest.fail`/`pytest.raises`; NO
  bare assert. The existing `test_mms_ld_2d.py` is all `-O`-safe — mirror it.
- **Mode-11**: the threading gate (Gate B/C) must EXECUTE the production
  `case.prescribed_inflow` (the producer stamp), not just the test surrogate —
  the #251 surrogate was producer-blind; the production-MMS re-target closes it.
  Sentinel-confirm the moment-resolved `prescribed_inflow`→LD boundary closure
  is on the call graph (Nexus `callers`/grep the moment-slot path).
- **Mode-10**: the boundary-slope term is ACTIVATED + CONSTRAINED by the
  consumed-flip mutation (Gate B). The companion-isolating value gate is
  UNAVAILABLE (§1) — the structural pair (machine-precision threading +
  consumed-flip ≫ tol + scalar no-op control) is the COMPLETE resolution. This
  is the SAME sub-case the #251 spec already added to the vv Mode-10 row.

## 6. SELF-IMPROVEMENT — NO new vv mode, NO ERR

This is the THIRD Mode-10 sub-floor-companion-unavailable instance (#240 D5b-S4
slope-source → #247 Leg A → #251 Leg B → S9). The vv Mode-10 row ALREADY
carries the "if no regime makes the term the dominant forcing … structural pair
is the complete resolution" branch (added by #251). S9 REINFORCES it; no skill
edit needed. NO new failure mode. NO ERR until a caught production bug
(next free ERR-063). The NOVEL finding worth recording (in THIS note, not the
skill): I tested the user's EXACT non-trivial-boundary-slope hypothesis (sharp
transverse frequency, Probe B) and it STILL did not unlock a value gate — so the
sub-floor wall is FUNDAMENTAL to a boundary-trace moment (codimension-1,
measure-zero in the limit), not regime-specific. The face-REPRESENTATION error
DOES improve O(h)→O(h²) (Probe D) but does not propagate above-floor into the
converged flux. This is the sharpest statement of the Mode-10 boundary-trace
sub-case to date.
