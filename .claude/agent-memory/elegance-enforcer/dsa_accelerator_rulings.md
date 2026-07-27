---
name: dsa-accelerator-rulings
description: Design rulings from the consistent-DSA (#2) Phase-3b review — reusable for 3c (theory page) and any future accelerator/inverse-operator sibling
metadata:
  type: project
---

# Consistent-DSA (#2) accelerator surface — durable rulings

Reviewed the Phase-3b DSA carve (`feature/sn-dsa`, uncommitted) 2026-07-26:
`orpheus/sn/acceleration/dsa.py` (`DSALowOrderSystem` + `DSACorrection`),
the `SourceIteration(corrector=...)` hook, the Krylov DSA preconditioner, and the
`integrate_angular` single-source hoist. Verdict: NO code BLOCKER — the code is
genuinely excellent; the one live finding was a tree-wide L16 stale-doc-contract.

**Why:** these rulings recur in the 3c theory-page review and when reviewing the
Green/Matrix inverse siblings, a P1-DSA extension, 2-D corner moments, or a second
accelerator. **How to apply:** lead with these when a new accelerator/low-order
carve lands; they are what "elegant" looks like for this family.

## What the carve NAILED (reinforce; these are the bar)

- **Derivation-pin = the L-005 gold standard.** Production `_build` is pinned
  entry-for-entry (atol 1e-15, `scattering_order=1` so P1-D is exercised) against a
  SYMBOLIC derivation-of-record (`orpheus/derivations/discrete/sn/dsa.py`
  `build_consistent_dd_system`, which PROVES Larsen (23a-f)/(27)/(38)/(39) rather
  than transcribing). Independent structural reference, not a self-snapshot — verify
  the pin exists (`tests/sn/acceleration/test_dsa_low_order.py::TestProductionTie`).
- **The 3-P0-frame anti-mint verdict** (R = `integrate_angular` = the ℓ=0 analysis
  face of `Quadrature.angular_frame(0)`; P = the `/Σw` normalized isotropic
  injection = the ℓ=0 reconstruction). MINT NOTHING — a third moment-0 spelling is
  the Smell-16 twin. Pinned 0-ULP by the D8 gates (`integrate_angular` ≡ frame row ≡
  displacement tangent-map, `assert_array_equal`). This is the correct resolution of
  "the correction needs a moment-0 reduction."
- **The `integrate_angular` single-source hoist** (`_bases.py`
  `_integrate_angular_values` — the ONE `einsum("n,ng...->g...")` body; `AngularFlux`
  wraps it byte-identically, the NEW `AngularDisplacement.integrate_angular` wraps the
  SAME body as its tangent map → `ScalarDisplacement`). Textbook Pattern 2 + "a linear
  reduction is its own tangent map." Exemplary.
- **The torsor discipline (#208 affine algebra).** `DSACorrection.apply` returns a
  DISPLACEMENT-typed composite (interior `AngularDisplacement`, boundary
  `AngularBoundaryDisplacement`) regardless of input role; the update is the torsor
  `ψ ⊕ Δψ` via `Composite.__add__`'s leaf-level `flux + displacement → flux` (single-
  sourced in `full_field.py` `_check_partner`; `flux + flux` stays a TypeError). The
  corrector is a `LinearOperator[V]` optional param on `SourceIteration` — the no-op
  clean-before-extend extension, NOT a sibling driver class (which would twin the loop).
- **Teeth that bite** (`test_dsa_acceleration.py::TestTeeth`): in-process monkeypatch
  (never a file mutation — matches process-discipline), sign-flip + trace-arm-zero, and
  crucially **residual-VALUE witnesses not iteration-count bars** (a diverging run
  reports max_inner−1 → off-by-one-prone; the value witness ~1e+35 vs <1e-6 is robust).
  D6 pins correction→0 at machine zero (`assert_array_equal`, no tol).

## Ratified DESIGN rulings (do not relitigate; carry into 3c + siblings)

- **The correction consumes the sweep DISPLACEMENT, never the typed residual.**
  `Δψ = ψ^{l+1/2} − ψ^l` → `d0 = interior.integrate_angular().values`. (This is what
  the old forward-promises got wrong — see the stale-contract finding below.)
- **The low-order operator is the DERIVED edge-centered SN-side system, NOT the
  in-algebra `A_diff = L+C−S−B`.** R4 (2026-07-26): the standalone cell-centered
  diffusion loss (`orpheus.diffusion`, RT0/harmonic-mean) measurably DIVERGES as an
  accelerator (ρ up to 54.7, `.claude/plans/dsa_d2_characterization.md`). Two defining
  laws (standalone-accuracy vs consistency-by-derivation), NOT a twin path. Build home
  is SN-side because the coefficients are SN-discretization properties (discrete γ_N,
  W2, scheme ρ).
- **`scattering_order` threads into the build; σ_s1 enters ONLY when the sweep retains
  ℓ≥1.** Consistency is with the ITERATED truncated operator, not the mixture data (a
  P0-truncated sweep's consistent partner carries bare-P0 `D=1/(3σ_t)` even if the
  mixture has P1 rows). `residual_sig_s()`/`foldable_sigma()`'s first production
  consumer — legitimate HERE (correction→0-safe) and only here (the #215 fold-into-
  sweep trap changes the anisotropic fixed point).
- **The trace arm is load-bearing for reflective walls** (discovered in flight, pinned
  by the trace-zero tooth): the production splitting lags reflection as the B gain
  reading the ITERATE's outgoing trace, while Larsen's reflecting row (f1=0) models an
  error equation tracking the CORRECTED state — leaving the trace uncorrected diverges.
  Inject the wall-EDGE solutions `f0[0]/f0[K]` into the boundary-trace displacement.
- **Accept the Krylov dual-interior** (`DSACorrection.apply` takes `AngularDisplacement`
  SI OR `AngularFlux` Krylov): GMRES vectors are role-erased at the scipy `from_flat`
  boundary and rebuilt flux-typed; the swept vector IS the displacement-from-zero. Both
  route through the ONE `integrate_angular`. Honest, not a two-hat throw.

## The one live finding — a clean L16 (stale deferral contract, tree-wide)

The carve LANDED the #2 capability in a DIFFERENT shape than every forward-promise
predicted, and 3 MAIN-TREE doc spots still assert the OLD shape (residual-as-substrate
+ in-algebra `A_diff` as the correction):
- `docs/theory/foundations/field_algebra.rst:522-529` — the strongest: names
  `A_diff = L+C−S−B (now built #290 P4; :doc:/theory/methods/diffusion_1d)` as the
  correction mechanism → points a fresh session AT the divergent operator R4 rejected.
  Also "typed residual is the DSA substrate" + "`as_dsa_source` lands WITH #2" (NOT
  minted — correctly, a third moment-0 spelling).
- `field_algebra.rst:101-102` — "the consistent-DSA substrate" (residual, wrong).
- `operator_tensor_network.rst:565-571` — "DSA consumes the fused residual + a separate
  in-algebra `A_diff`" (both wrong).

Author's deferral plan named only 1 of 3 (the `as_dsa_source` line) → the L16
half-cleanup tell. I DISAGREED with deferring the false-MECHANISM content to 3c: the
full theory page defers, but a confidently-false claim pointing at the divergent
operator is Cardinal-Rule-3 doc-brain contamination — neutralize all 3 now (forward-
pointers) while the R4 context is hot. (`as_dsa_source` NOT minted is a CORRECT
deviation. A sibling worktree `nexus-workspace-wiring` carries the same stale line at
`operator_algebra.rst:3155` — that branch owns it, out of this review's scope.)

## NITs (recorded, non-blocking)

- `dsa.py:327` boundary g0 weight spelled `-sgn*2.0*half_gs*0.5` (= `-sgn*half_gs`,
  bit-identical power-of-2 scaling); interior rows + the reference builder spell the
  same ½g0 as plain `half_gs`. One quantity, two spellings (Pattern 7 readability).
- Trace face keys hardcoded `{"xmin","xmax"}` (`dsa.py:552-559`) — ACCEPTABLE in the
  1-D-admitted arm (`from_face_arrays` validates keys vs the layout, raises loud; the
  wall-edge→face map is genuine 1-D physics). Do NOT lift the FaceField ABC now
  (premature abstract-over-the-difference). Collapse trigger: 2-D corner moments.
- `acceleration: str|None` (`solver.py`) — house-consistent with `inner_solver`/
  `inner_schedule`; anti-#4 applies to the whole family (tracked dispatch-spelling
  debt #259-262), don't type in isolation.

## 3c rate/stability TIER review (tests-only, `test_dsa_rate.py` 65 tests, 2026-07-26)

Verdict PASS / CONCERNS — NO violation, NO rework. Tests-only on the 3b I passed;
every twin is coextensive-today. Durable rulings for 3d / any gate-tier review:

- **lru_cache-over-raw is a MUTATION FOOTGUN → the elegant guard is Pattern-4, not a
  docstring.** `_uniform_solve = lru_cache(_uniform_solve_raw)`; teeth MUST use `_raw`
  or a monkeypatched call served a healthy cached result is a VACUOUS GREEN (a gate
  falsely reporting teeth — the worst L-005 failure). Verified the colliding key is
  live (`test_si_dsa_speedup_floor` caches `(0.9,0.5,vac,dsa)`, exactly what a future
  DSA-mutation tooth wants). CONCERN not VIOLATION (convention holds today). The
  structural fix that makes the illegal state unrepresentable: an **autouse fixture
  keyed on `"monkeypatch" in request.fixturenames`** that `cache_clear()`s in BOTH
  setup and teardown — closes vacuous-green (setup clear ⟹ recompute under patch) AND
  forward-poisoning (teardown clear) regardless of which wrapper the author calls.
  ~free here (only 1 monkeypatch test, uses `_raw` anyway). Reusable pattern for any
  cached-healthy + raw-mutation test split.
- **A negative-control / spectral instrument that produces DELIBERATELY-WRONG pairings
  belongs TEST-LOCAL, never in `orpheus/derivations/` (the algebra of RECORD).** The
  WD sweep `_wd_sweep_matrix` (parametrized by `a`, a family production doesn't have)
  is anchored SUFFICIENTLY by (a) the a=½ S2 machine-zero composite anchor
  (`test_dd_member_anchor`, ρ<1e-12) + (b) the production-chain S2 exactness
  (`TestS2Exactness` via the real solver) jointly pinning sweep⊗low-order⊗cell-update
  ⊗trace. Key discipline the control gets RIGHT: `_consistent_low_order` calls the
  PRODUCTION `DSALowOrderSystem._build` (not a re-transcription) — only the SWEEP and
  the matrix cell-average (`c2e`) are test-local, and both are S2-pinned. NIT on such
  instruments: a docstring "reproduces the DD sweep to 1e-17" that only the COMPOSITE
  ρ checks (L-001: prose names a stronger/direct identity than the code asserts); and
  a dead BC arm (the `_wd_sweep_matrix` reflective branch — all 4 callers pass vac/vac
  ⟹ YAGNI-drop, anti-#10).
- **A hand-rebuilt mutation is a PARTIAL TWIN of the production builder — verify it is
  faithful TODAY, then demand cross-ref + collapse-trigger (L-002).** The lumped-mass
  tooth rebuilds ~40 lines of `_build`'s tridiagonal + Marshak/reflective rows,
  differing ONLY in the mass distribution (verified line-by-line: identical `dh`
  leakage, identical G reused from production, row-sum-conserving lump). Legitimate
  mutation TODAY; CONCERN = a future `_build` leakage/Marshak/γ_N change desyncs the
  twin silently (the 3× count baseline becomes apples-to-oranges). Ideal fix = build
  the mutation as a **mass-delta ON TOP of production `a_low`** (faithful-by-
  construction) — clean for the INTERIOR (mass adds directly) but the boundary lump
  entangles `c1`/orient, so the pragmatic do-now is the reciprocal cross-ref +
  collapse trigger.
- **AST routing-sentinel design (reusable fence pattern).** `_fold_accessor_callers`
  scans `orpheus/` for `ast.Attribute` reads OR `ast.FunctionDef` defs of the fenced
  names, diffs vs an allowlist. Confirmed allowlist == live tree (3 files). Two
  exemplary legs: (1) the SELF-CHECK — `dsa.py must be in callers` reddens if the
  scanner sees nothing (a broken scanner can't pass vacuously, L-005 structural); (2)
  the planted-consumer tooth proves the tripwire bites. Docstring claim "docstrings
  don't register (string constants)" VERIFIED true (contrast L-001). NOTE: the
  FunctionDef arm mixes definition-sites with consumer-sites in the allowlist — the
  pure threat surface is attribute READS; including defs is defensible defense-in-depth
  (catches a definition MOVE to a new module). PASS.
- **NIT recurring in physics gates: a scattering ratio `c` spelled as a bare literal.**
  `RHO_BOUND * 0.9` / `RHO_BOUND * 0.5` (rate.py:230,243) MEANS the A&L bound
  `0.2247·c`; writing the literal (not `RHO_BOUND * c`) obscures that the band edge IS
  the Fourier bound at this c and reads like a magic safety factor (Pattern 3/7).
- **Docstrings-as-interim-record are exemplary (Cardinal-Rule-3 brain) but carry a
  migration trigger.** The file doubles as the 3c record until the theory page lands —
  physics-dense, every band pinned to a measured number, failure modes named
  (Signature 9, Mode-9, #215/ERR-070). NOTE: when the theory page lands, the measured-
  facts prose MIGRATES (not duplicates → else a twin doc-source).
## 3d DELTA review — P1-DSA d₁ arm + ERR-071 sweep-inverse fix (2026-07-27)

Two production deltas on `feature/sn-dsa`, both green + pyright 0/0/0. Verdict
PASS with 1 must-fix (a DOC sign contradiction, not code) + 2 NITs.

- **The anti-mint ruling EXTENDS to ℓ=1: the ℓ=1 analysis row must route through
  the frame single-source exactly as d0 does through `integrate_angular`.** The P1
  arm open-codes `self._w_mu = w*mu` and `d1 = einsum(w·μ, interior.values)` while
  the docstring CLAIMS "the frame's ℓ=1 analysis row … already exist[s]" and the
  frame DOES provide it (`angular_frame(1).table[:,1,1] == μ`, test-pinned 0-ULP in
  `test_moment_pair_injection_object_pins`). This is L-001 + a symmetry break: d0 is
  a FIELD reduction (`interior.integrate_angular()`, encapsulated, frame-backed), d1
  reaches into `.values` + open-codes the contraction. Test-pinned coextensive ⟹ NIT
  (do-now: source the coeff from the frame / a field-level ℓ=1 reduction); the
  build-the-machinery destination is the ℓ=1 sibling of `integrate_angular` so BOTH
  moments are frame-backed field reductions. Same family: `_inv_w2` recomputes W2
  that `_build` already computes (guard-pinned ≈1/3 ⟹ very minor). The (33)
  RECONSTRUCTION open-coding IS accepted (per my 3b bar: reconstruction may be
  open-coded IFF R∘P=I is pinned — and it is). REINFORCE: `_inv_w2` is COMPUTED from
  the quadrature not hardcoded 3.0 (Pattern 14 ✓); `_a_coef` single-sources the (23f)
  weight across g1+‌(28b); the P0/P1 branch is ESSENTIAL (verified `moment1_update`
  ≠ 0 at d1=0 via the −(D/h)Δf₀ term, so P0 byte-identity genuinely needs the branch,
  NOT an anti-#7 special-case).
- **A sweep-inverse "restore the dropped rhs rows" fix belongs in the SOLVE method,
  not "deeper" in the generic march.** ERR-071: `(L+C)⁻¹` seeded the rhs then the
  march clobbered the outflow-trace rows; fix = one post-march
  `boundary_buf[outflow] -= seed[outflow]` via the typed partition
  (`angular_trace.outflow_indices_for_face`). PLACEMENT verdict PASS:
  `_solve_timed_full_field` is the single home that assembles the inverse's output
  across ALL row-classes (inflow=seed, interior=march, outflow=restore,
  tangential=untouched) — the restore is its native outflow-row step. Pushing it into
  the generic sweep/march would entangle the generic interior march with THIS
  operator's boundary-row convention (worse coupling). Forward/inverse consistency is
  LOCKED by the new round-trip gate's pure-outflow leg (single-sources the sign). The
  restore uses the SAME typed partition the forward `loss_action` uses (single-sourced
  "which rows"). This is symmetry-in-math→symmetry-in-code: forward row `y = streamed
  − x` ⟹ inverse `x = streamed − y`, both spelled `streamed − (·)_out`.
- **The MUST-FIX was a SIGN doc contradiction (code correct, docs disagree) — the #1
  failure mode in newly-committed V&V scaffolding.** Authoritative forward =
  `streamed − ψ_out` (`loss_representation/__init__.py:1164` `streamed[out]-given[out]`,
  confirmed by its own :1085/:1195 docs). The fix code (`-= rhs_out` → `streamed −
  rhs_out`) is SIGN-CORRECT. But `streaming.py:58` (class docstring, pre-existing) says
  `psi.outflow - streamed` (negated), and the NEW identity-gate docstring
  SELF-CONTRADICTS ("DEFECT `ψ_out − streamed`" vs parenthetical "`streamed − ψ_out`").
  A green round-trip gate is SIGN-BLIND (only proves apply/solve AGREE) — so I read the
  forward code to settle it, never the gate. Ruling: harmonize ALL doc mentions to
  `streamed − ψ_out`; the new gate's self-contradiction is in-scope must-fix, the
  pre-existing :58 is L-004 (the delta made the convention acute → reconcile now).
- **A cache-guard's docstring must not overclaim its coverage.** My applied IMPROVE-a
  autouse `_no_cache_across_mutation` clears ONLY `_uniform_solve` but claims
  "regardless of which solve wrapper the tooth calls." `_p1_solve` is UNCACHED (safe
  today) so complete NOW — but a future `@lru_cache` on `_p1_solve` + a P1 mutation
  tooth reopens the vacuous-green footgun the guard exists to close, and its own
  docstring would have lied. NIT: generalize the clear to all module caches, or scope
  the claim. (Un-mentioned 3rd delta `transport/method.py` = docstring-only, retires
  the R4-falsified "DSA driver #2 first full_field_space consumer" forward-promise —
  good proactive L-016 hygiene, PASS.)

## Helper dup (3c) — deferred

- **Helper dup (`_phi`,`_inners`,`_TOL`,`_FP_RTOL`) across rate+acceleration files:**
  rule-of-two MET, but author already applied the L-002 cross-ref ("same rationale as
  the 3b file") ⟹ DEFER is defensible; extract-now marginally better scoped to a
  `_helpers.py` (NOT `conftest.py` — these are plain fns/consts) for the introspection
  layer + bands ONLY. The config `_solve`/`_uniform_solve` helpers stay per-tier
  (different mesh/S-order/mixture — L-003 different-consumers, NOT a share target).
