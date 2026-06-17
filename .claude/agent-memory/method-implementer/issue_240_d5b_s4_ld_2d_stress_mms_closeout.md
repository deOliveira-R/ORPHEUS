# #240 D5b-S4 — 2-D Cartesian LD stress MMS (closeout)

Branch `feature/sn-space-angle-tier2` (NOT committed — main agent commits after
review). Host env, `.venv/bin/python`, canonical `python -O -m pytest`.

## What this delivers

The L1 flux-shape verification of the multi-D bilinear (UBLD) Linear-
Discontinuous closure (the S3-landed `(L+C−S_full)` moment matvec): a vv-Mode-7
strengthened MMS whose trial flux is μ-bilinear (activates the per-axis SPATIAL
slope rows) with a NON-vanishing boundary trace (`a0>0` → prescribed inflow).
Closes the **slope-UNKNOWN** half of the LM-1989 slope-row sign trap in d≥2.

## Files changed (manifest)

- `orpheus/derivations/continuous/mms/sn.py` (+592):
  - `_LD2D_STRESS_COEFFS` — single-sourced spatial harmonic amplitudes
    (numerator/denominator pairs) shared by Branch 1 (SymPy `Rational`) and
    Branch 2 (numpy float).
  - `_ld2d_stress_amplitudes()` — Branch-2 float reader of the constant.
  - `_2d_cartesian_ld_stress_symbolic()` — Branch-1 SymPy builder (the
    algebra-of-record; returns ψ, φ=A, Q_closed for the strengthened ansatz).
  - `derive_2d_cartesian_ld_stress_mms()` — Branch-1 substitution-identity
    derive function (residual `simplify == 0`).
  - `SN2DCartesianLDStressMMSCase` — Branch-2 dataclass: `_drivers` (A/B/C +
    derivatives), `phi_exact`, `build_mesh` (non-square vacuum), `build_materials`
    (per-cell het 2G via `_default_hetero_2d_xs_functions`), `external_source`
    (FLAT per-ordinate `(N,ng,nx,ny)`), `prescribed_inflow` (non-vanishing
    face-average trace via `BoundarySourceSink.prescribed_inflow`).
  - `build_2d_cartesian_ld_stress_mms_case(...)` — factory (non-square 1.3×0.9,
    2G c=(1.0,0.4), default level-symmetric S4 / no pure-z; pass Lebedev for the
    pure-z matvec gate).
- `tests/derivations/test_sn_mms_ld_2d_stress_symbolic.py` (NEW, 8 foundation
  tests): SymPy substitution identity (`@verifies("ld-cartesian-2d",
  "transport-cartesian-2d")`) + INDEPENDENT finite-difference residual check
  (structural independence vs SymPy's own `diff`, L11) + φ=A + slope-driver
  per-axis activation + B,C-break-x↔y-reflection + Branch2==Branch1 source
  cross-check (1.5e-16) + prescribed-inflow non-vanishing + structural-
  independence-of-LD-kernel guard.
- `tests/sn/verification/mms/test_mms_ld_2d.py` (+177): D5b.2 headline
  (`test_ld_2d_stress_converges_second_order`, `@l1 @slow
  @verifies("ld-cartesian-2d","transport-cartesian-2d")`) + D5b.4 matvec twin
  (`test_ld_2d_stress_krylov_equals_si`, `@l1 @verifies("ld-cartesian-2d")`) +
  D5b.3 stress upgrade (`test_ld_2d_stress_two_paths_ffw_equals_mfw`, `@l1
  @verifies("ld-cartesian-2d")`, ADDED alongside the existing random-source
  vacuum two-paths pin — did NOT replace it; the random-source pin is an
  orthogonal coverage). Imports `build_2d_cartesian_ld_stress_mms_case`,
  `build_nonvacuum_fixed_source`.
- `docs/theory/discrete_ordinates.rst` (+73): the D5b-S4 stub — `:label:
  ld-cartesian-2d` (NEW, minted here) + the strengthened-ansatz math + `:mod:`/
  `:func:`/`:class:` cross-refs + an archivist TODO + an HONEST-SCOPE note
  (slope-SOURCE half deferred). `docs/verification/matrix.rst` auto-regenerated
  by the Sphinx build (records `ld-cartesian-2d, 4`).

## THE moment-source posture finding (integration point a) — the crux

The production does NOT consume an external slope-moment source Q̂ in d≥2.
Probe-verified (`/tmp/probe_ld_source_moment.py`):
- `solver.py:_lift_external_source_to_moments` lifts a flat `(N,ng,*spatial)`
  external source onto slot 0 (average), **slope rows ZEROED**. The code itself
  documents this: "the strengthened Q̂≠0 ansatz is S4" (`solver.py:1884-1885`).
- The public `solve_sn_fixed_source` entry HARD-VALIDATES the bulk against the
  flat `(N,ng,*spatial)` shape (`solver.py:1876`); a moment-resolved
  `(N,ng,nx,ny,2^d)` external source RAISES `ValueError`.
- The ONLY moment-valued source the production consumes is the iterate-driven
  SCATTERING source `Σ_s·φ̂` (S3's `S_full`), NOT an external Q̂.

So the LM-1989 trap is **NOT fully closed** in d≥2: the slope-UNKNOWN half IS
verified (B,C drive non-trivial per-ordinate fields → the bilinear closure
solves ψ̂_x, ψ̂_y from the average + scattering, always active). The
slope-SOURCE sign half is DEFERRED. This is a SCOPE finding, NOT an S3
integration gap (the production correctly consumes the flat external source +
the scattering slope source — it simply has no entry point for an external Q̂).
Frame 2 §232 anticipated exactly this and instructed: document honestly, do NOT
claim the trap is closed. Done (case docstring + theory note + the gate
docstrings).

## THE boundary-trace finding (integration point b)

The non-vanishing prescribed inflow IS supplied (via
`BoundarySourceSink.prescribed_inflow` + `build_nonvacuum_fixed_source`, the
Phase-4.6 1-D precedent lifted to 2-D) — the boundary closure is stressed at
the FACE-AVERAGE moment (a0>0). But the boundary trace `mesh.trace` carries one
SCALAR per face per ordinate per group (slot shapes `(N,ng,n_transverse)`, NO
`2^{d-1}` transverse-moment axis); `loss_representation._inflow_to_moments`
seeds the scalar onto the face-average moment and zeros the transverse face
slope (the loss-rep docstring flags this as "the #240 D5b-S4 boundary
widening"). So the transverse face-SLOPE inflow moment is DROPPED — same
deferred status as the slope-source. NOT a hack-around: the average-moment
boundary closure is genuinely exercised; the transverse slope is a production-
wiring increment (a moment-resolved trace) beyond this test's scope.

## Mutation verification (Frame 2 §316-327) — the strengthening is load-bearing

Mutated the LD slope-row SIGN on BOTH axes identically (runtime monkeypatch,
reverted in `finally`; `_ubld.py` git-clean confirmed). Three faithful
mutations of the multi-D UBLD slope-row sign:

| Mutation | strengthened result | verdict |
|----------|--------------------|---------|
| `_GRAD_1D` full sign flip | NaN | CAUGHT |
| `_FIN_TRACE` slope `[1,-1]→[1,1]` | order −4.62, finest 20.3 | CAUGHT (weak also caught) |
| `_GRAD_1D[1,0]: -2→+2` (surgical slope-row) | inf | CAUGHT |

Baseline strengthened: order 2.00–2.14, finest 3.5e-3..6.0e-3. The UBLD
bilinear closure is tightly coupled (the slope feeds the face cochain that
propagates through the wavefront), so a slope-row sign error DIVERGES the
iteration rather than producing a subtle wrong-limit — a STRONGER guarantee
than Frame 2's subtle-cancellation scenario (there is no false-green regime to
hide in for this closure). The `_GRAD_1D[1,0]` flip is the most faithful
"same-sign slope-row sign on both axes" (the Kronecker build applies `_GRAD_1D`
per-axis identically). Strengthened MMS catches all three definitively.

## Gate results

- `python -O -m pytest tests/derivations/test_sn_mms_ld_2d_stress_symbolic.py`
  → **8 passed** (the foundation/L11 derive gate).
- `python -O -m pytest tests/sn/verification/mms tests/derivations/test_sn_mms_ld_2d_stress_symbolic.py`
  → **50 passed** (459s; includes D5b.2 headline O(h²) order=2.00,
  D5b.4 Krylov≡SI, D5b.3-stress FFW≡MFW + all pre-existing MMS).
- No-regression: `python -O -m pytest tests/sn/spatial tests/sn/sweep/core
  tests/sn/sweep/cartesian_2d tests/sn/solve -W
  "error::tests.sn.regression._regression_assert.DriftWarning"` →
  **632 passed, 2 skipped, 4 xfailed** (DriftWarning gate clean).
- Sphinx `-b html` → build succeeded (only pre-existing test-file
  SyntaxWarnings; `ld-cartesian-2d` label present in HTML, 3 hits; no
  duplicate/undefined-label warning).
- V&V audit → `ld-cartesian-2d  4 test(s)`; no new orphan; audit exit 0.
- Branch2==Branch1 source cross-check: max |Δ| = 1.5e-16 (machine precision).

## ERR-NNN

None caught (no production bug; the moment-source/boundary-trace posture is a
documented deferred-scope decision, not a bug). The mutation verification is a
throwaway — reverted, not committed.

## Honest scope (what is verified vs deferred)

VERIFIED: multi-D LD bilinear slope-UNKNOWN sign (O(h²) flux-shape) + the
AVERAGE-moment prescribed-inflow boundary closure + the matvec twin (Krylov≡SI)
+ the two-DAG-schedule invariant, all on the strengthened μ-bilinear het-2G
non-square non-vacuum config; mutation-verified non-vacuous.

DEFERRED (the slope-SOURCE half of the LM-1989 trap + the transverse face-slope
inflow moment): needs (1) a moment-resolved EXTERNAL-source entry —
`solve_sn_fixed_source` must accept `(N,ng,*spatial,2^d)` and
`_lift_external_source_to_moments` must thread the slopes instead of zeroing
them; and (2) a moment-resolved boundary trace — `mesh.trace` /
`_inflow_to_moments` must carry the `2^{d-1}` transverse face-moment axis. Both
are production-wiring increments beyond the S3 moment matvec. The manufactured
Q̂ + transverse face-slope ARE derived (Branch 1 is slope-source-ready) — only
the production consumption path is missing. This is a candidate next sub-step
(call it D5b-S5 / "moment-source consumption"); flagged for the main agent.

## Lessons (for the skill)

- **algebra-of-record / single-source**: when Branch 1 (SymPy) and Branch 2
  (numpy) share spatial constants, single-source them as `(num, den)` pairs
  (`Rational` for SymPy, exact float for numpy) — the Branch2==Branch1 cross-
  check then pins the two EVALUATORS agree, not just the symbolic identity. The
  amplitudes cannot drift.
- **vv Mode-7 + the production posture is a load-bearing pre-design probe**: an
  MMS that supplies a slope-moment source is only as good as the production's
  CONSUMPTION path. Probe `_lift_external_source_to_moments` /
  `_inflow_to_moments` BEFORE designing the ansatz's Q̂ — the production may
  zero the very moment the ansatz manufactures. Map "does the entry point
  consume the moment I'm manufacturing?" as a crosswalk row for any MMS that
  targets a multi-moment closure.
- **mutation on a tightly-coupled DG closure diverges, doesn't subtly-wrong**:
  the Frame-2 subtle-cancellation false-green scenario does not arise for the
  UBLD wavefront (the slope feeds the propagating face cochain). The catch is
  catastrophic (NaN/inf), which is a clean mutation signal but means the
  weak-ansatz-false-green contrast is moot — both catch a slope-row sign flip.
  The strengthening is still load-bearing for the OTHER plausible-wrong class
  (a finite x↔y-symmetric slope error confined to the slope channel), which the
  available mutation sites couldn't isolate; document the catch evidence and
  the limit.
