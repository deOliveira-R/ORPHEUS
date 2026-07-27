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
