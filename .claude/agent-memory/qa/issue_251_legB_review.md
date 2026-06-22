---
name: issue-251-legB-review
description: "#251 Leg B QA review (BOUNDARY transverse face-slope — the boundary twin of #247 Leg A's bulk slope source, for 2-D Cartesian LD). VERDICT SUPPORTED, NO ERR, NO blocker. Branch refactor/sn-foundation-cleanup, UNCOMMITTED. Mode-11 closed (gates drive the production moment-resolved arm 688× on the live public solve + go RED against emulated pre-#251 zero-fill). DD/Step + scalar-inflow byte-identical (sweep/core 460/1/4 strict; 1-D prescribed-inflow 4 green). improves-on-flat genuinely sub-floor (probed flat 2.131e-2 < real-slope 2.163e-2 → drop is honest). Reflective storage pass-through corruption-free; SIGN correctly deferred to #252. 1 minor scope-note: DD identity at n==1 ≠ spec D6 'DD rejects any moment trace' (not a defect — DD never receives a moment inflow)."
metadata:
  type: project
---

# #251 Leg B — BOUNDARY transverse-face-slope (the LM-1989 trap boundary half) — QA review

Branch `refactor/sn-foundation-cleanup`, UNCOMMITTED. Host `.venv/bin/python -O`.
The BOUNDARY twin of [[issue-247-legA-review]] (the BULK slope source, committed
`d9396a2`). 4 files: `geometry.py` (`boundary_face_layout` moment-tail widening),
`loss_representation.py` (`_inflow_to_moments` rank-discriminated pass-through +
`_octant_face_cochain` seed + 4 dropped outflow capture-collapse sites),
`boundary_source_sink.py` (`prescribed_inflow` scalar-or-moment slot assign),
`test_mms_ld_2d.py` (#251 block — 2 xfail-strict gates flipped to live + 4 green).

## OVERALL VERDICT: SUPPORTED. NO ERR-063. NO blocker. Commit OK.

All 4 claims SUPPORTED with mutation evidence. ONE minor scope-note (DD identity
vs spec-D6 wording) — a documentation/spec nit, NOT a defect, NOT commit-gating.

## Claim-by-claim

### Claim 1 (boundary Mode-10 gap closed + sign constrained) — SUPPORTED
The transverse face-slope is now CARRIED (no longer zeroed) AND its sign is
CONSTRAINED (a flip is catchable). Proven by TWO in-process mutations (NO prod edit,
L28-safe), both via the canonical `-O` invocation:
- **Re-zero the slope rows** (the EXACT #251 bug `f[..., 1:] = 0.0` in the
  moment-resolved arm): threading gate RED (AssertionError), consumption gate RED
  (`|Δφ|/|φ|=0.000e+00 ≤ 1e-08` — flipping a zeroed slope is a no-op = the exact
  Mode-10 signature), scalar no-op control GREEN (Leg-B asymmetry), width-reject GREEN.
- **Emulate the PRE-#251 unconditional zero-fill** (no rank discrimination →
  treats a `(...,2)` moment face as scalar → appends a spurious `(...,2,2)` axis):
  threading gate RED ("did not RECOGNISE the moment-resolved inflow"), width-reject
  gate RED ("DID NOT RAISE ValueError"), scalar no-op GREEN. Both previously-xfail
  gates genuinely REQUIRE the #251 change to pass — not vacuously green.

### Claim 2 (scalar-inflow + DD/Step byte-identical) — SUPPORTED
- LAYOUT level: DD trace slots `(24,2,6)`/`(24,2,8)` — NO trailing moment axis
  (`face_moment_tail(1)==()`); LD slots `(24,2,6,2)`/`(24,2,8,2)` (trailing 2 =
  `2^{d-1}`). Confirmed live.
- SNAPSHOT level: `tests/sn/sweep/core` 460 passed / 1 skipped / 4 xfailed under
  `-W error::DriftWarning` — ZERO drift escalation (DD/Step trace untouched). (Main
  agent runs the full `sweep/core`+`solve` 520/1/4 baseline in background.)
- SCALAR-INFLOW path: `test_mms_prescribed_inflow.py` 4 passed (1-D slab face is a
  point → `n_face_moments==1` → identity widening). Full-solve scalar no-op
  (`slope_sign=0.0`) vs today (`slope_sign=None`) maxdiff = **0.000e+00** (reproduced
  independently). The `prescribed_inflow` producer discriminates by EXACT shape: full
  slot (DD's scalar slot IS the full slot) → byte-identical write; scalar-onto-moment
  → seed slot 0 only.
- ORACLE: the two-paths FFW≡MFW gates + Krylov≡SI + external-slope (Leg A) +
  scattering all PASS (9/9) — the 4 dropped capture-collapse sites preserve
  schedule-equivalence (both legs route through the same widened cochain).

### Claim 3 (teeth are STRUCTURAL; dropping "improves-on-flat" is correct) — SUPPORTED, sound & non-circular
INDEPENDENTLY re-probed against the analytic `phi_exact` reference at nc=16:
```
near-bdy A-err: flat=2.1311e-02  real-slope=2.1630e-02  flipped=2.1185e-02
real-slope improves on flat? False
consumption flip moves near-bdy flux: 4.101e-3  (TOL=1e-8)   noise floor (no-op maxdiff) = 0.0
```
The real transverse slope makes the converged near-boundary error SLIGHTLY WORSE
(2.163e-2 > 2.131e-2) and the flipped slope is slightly BETTER — the boundary
correction is genuinely sub-floor below the bulk O(h²) error. So "improves-on-flat"
is NOT achievable; dropping it is HONEST, not hiding a problem (keeping it would
falsely FAIL while the slope is correctly consumed). The positive verification is
TWO O(1) structural signals, neither depending on the converged value matching A:
(a) machine-precision threading (`assert_array_equal` slot-1, reference built by
`leggauss` only — L11, NON-circular: production's arm is a pure pass-through
`widened.append(face)`, so the gate proves pass-through ≠ the old zero-fill); (b)
consumed-flip 4.101e-3 ≫ TOL by 5.6 orders, well above the deterministic noise floor
(0.0). The 4.101e-3 figure triple-agrees across three structurally-distinct
computations (my independent re-derivation, the test-architect surrogate 4.106e-3,
the method-implementer public-path 4.1012e-3). This is the canonical Mode-10
resolution where the term has no O(1)-dominant regime (the companion-gate half of the
recipe is UNAVAILABLE for a boundary-trace slope). NEW Mode-10 sub-case → the
test-architect's one-line skill addition is warranted (see Self-improvement below).

### Claim 4 (no latent reflective bug) — SUPPORTED, #252 is the right disposition
(a) **No storage corruption from the widening.** Built a reflective-xmin/xmax LD-2D
mesh, seeded a random moment-shaped trace `(24,2,6,2)`, ran `_reflect_trace("apply")`:
no crash, moment axis preserved (`(24,2,6,2)` out), and slot-1 follows the SAME
ordinate permutation as slot-0 for all 12 matched inflow ordinates → **cross-moment
corruption count = 0**. The reflective `B` is `PermutationOperator(perm, axis=0)`
broadcasting over trailing axes; the moment axis passes through for STORAGE without a
hard-coded trailing-axis assumption. NO latent storage bug introduced.
(b) **#252 is the right disposition.** The transverse-slope SIGN under reflection is
genuinely UNVERIFIED (the Leg-B MMS is vacuum-BC → H2 nulls the reflective coupling).
This is a Mode-1 sign trap the vacuum gates STRUCTURALLY cannot see — it MUST be a
follow-up, not blocked into #251 (whose shipped scope is the vacuum prescribed-inflow
path, fully verified). #252 is filed OPEN with correct labels (level:L1, module:sn,
module:tests, type:improvement) and a clear physics-expectation description.

## Specific V&V checks (all PASS)
- **L11 structural independence**: `_face_transverse_legendre` uses ONLY
  `numpy.polynomial.legendre.leggauss` + hand-laid algebra — NEVER
  `_inflow_to_moments`/`assemble_inflow_axis`/any LD cell op. The projection-correctness
  foundation sub-gate's reference is the hand-derived `[c0+c1·tc, (h_t/2)·c1]` (a
  linear is exact at q_nodes≥2). Non-tautology confirmed.
- **Anti-pattern #11 (positive AND negative)**: every gate has both legs. The
  consumption gate pairs the consumed-flip (RED on the bug) with the scalar-inflow
  no-op (the Leg-B asymmetry — GREEN under the bug because slot-1≡0). Sound.
- **Mode 11 (gate executes the rewired path)**: CLOSED. Instrumented the public
  `+slope` solve → `_inflow_to_moments` called 344×, the moment-resolved arm fired
  688× (0 scalar, 0 identity). The live gate drives the production moment-resolved
  arm END-TO-END (NOT a surrogate; the monkeypatch is fully removed). The 2
  previously-xfail gates now genuinely pass via production (proven RED against the
  emulated pre-#251 zero-fill).
- **Mode 8 (-O safe)**: ZERO bare asserts in the #251 test block (all
  `np.testing.*`/`pytest.fail`/`pytest.raises`) AND zero bare asserts added to the 3
  production files. Ran under `-O` (the "asserts not executed" warning fired); gates
  fire correctly.
- **Negative pin**: a wrong transverse-moment width (3 ≠ `2^{d-1}=2`) RAISES a
  ValueError naming `2^(d-1)` AND the bad width 3 (Pattern 4 — the relaxation does
  not swallow shape bugs). Correct width (2) passes through with slot-1 preserved.
- **ERR adjudication**: NO ERR warranted. The slope was UNVERIFIED, not WRONG — no
  sign/magnitude bug in the now-consumed path. Evidence: the consumed-flip magnitude
  triple-agrees across three independent computations, the threading is at machine
  precision, and the bulk O(h²) convergence gate still passes (exit 0). The cochain's
  transverse mass `diag(h_t, θh_t)` and the BARE-coefficient projection are
  apples-to-apples (the §1 crux held).
- **L-007 tagging**: 5 gates `foundation`-only (no `verifies`); 1 gate
  (`...sign_mutation_reddens`) `l1 + verifies("ld-cartesian-2d")` — the one
  full-public-solve consumption claim, correctly tied to the eq. No `catches`
  (correct — no caught bug).
- **Single source (Pattern 7)**: `face_moment_tail` is the ONE tail policy (the
  geometry lever reuses it); `is_moment_valued_by_flat_rank` is the ONE rank
  discriminator (windowed :401 + FFW oracle :1263, no open-coded 3rd `.ndim==` spelling).

## Mutation evidence (verbatim, the teeth)
Under the #251 re-zero-slope mutation (in-process plugin, `-O`):
```
FAILED ...::test_ld_2d_boundary_slope_threaded_through_inflow_to_moments - AssertionError
FAILED ...::test_ld_2d_boundary_slope_sign_mutation_reddens - Failed: ... |Δφ|/|φ|=0.000e+00 ≤ 1e-08 ...
2 failed, 2 passed (scalar no-op + width-reject GREEN = the Leg-B asymmetry)
```
Under the emulated PRE-#251 unconditional zero-fill (`-O`):
```
FAILED ...::test_ld_2d_boundary_slope_threaded_through_inflow_to_moments - Failed: ... did not RECOGNISE the moment-resolved inflow ...
FAILED ...::test_ld_2d_boundary_trace_rejects_wrong_transverse_width - Failed: DID NOT RAISE <class 'ValueError'>
2 failed, 1 passed (scalar no-op GREEN)
```
Mode-11 instrumentation (public `+slope` solve): `_inflow_to_moments` total=344,
moment_resolved=688, scalar=0, identity=0 → "Mode-11 CLOSED".

## Minor scope-note (NOT a defect, NOT commit-gating)
Spec §8/D6 said "DD/Step (`n_face_moments==1`) must REJECT any moment-resolved
inflow outright." The implementation early-returns IDENTITY at `n==1` BEFORE the
rank discrimination, so a `(...,2)` array passed to a DD `_inflow_to_moments` returns
unchanged (NO raise). The shipped negative-pin gate only tests the LD wrong-width
case, so it matches the shipped behavior. This is NOT a correctness defect: DD's
trace layout is scalar-only by construction (`face_moment_tail(1)==()`), so a DD
inflow tuple is ALWAYS rank `d+1` — a moment-resolved inflow never reaches a DD
`_inflow_to_moments` in production. The `n==1` identity IS the correct DD contract.
Flag for the spec author only (the D6 "DD rejects any moment trace" half was neither
implemented nor tested, by design); no action required for the commit.

## Self-improvement check
Two pushback-adjacent observations: (1) the Mode-10 "companion-gate-unavailable"
sub-case (a boundary-trace slope has no O(1)-dominant regime → structural teeth ONLY
are the complete resolution) is a NEW shape within Mode-10. The test-architect's §13
recommends a one-line addition to the vv Mode-10 row. This is genuinely novel
(neither #240 D5b-S4 nor #247 Leg A had a term with NO dominant regime — both could
improve-on-flat). RECOMMEND the main agent apply the one-line Mode-10 row addition
to vv-principles SKILL.md at delivery (the test-architect drafted the exact wording).
(2) The "verify a re-zero mutation reddens the PUBLIC-solve consumption gate, AND
emulate the pre-change behavior to prove the xfail→live flip is non-vacuous" recipe
is the Mode-11 closure pattern, already covered by L-027/L-031/lessons on
sentinel-instrumenting + mutating the documented bug — no new lesson needed beyond
folding this Leg-B instance into the existing Mode-10/Mode-11 entries.
