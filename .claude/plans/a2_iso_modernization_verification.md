# Verification plan — campaign #276 phase A2 (SN scattering iso-modernization + S†)

**Branch:** `feature/sn-adjoint-transport` · **Carve:** A2a (forward iso-through-frame
cleanup) + A2b (S† extend). Proactive `test-architect` dispatch per the operator-algebra-
crossing-subsystems trigger (per-ordinate ↔ moment, normal ↔ adjoint, fast-path ↔ composed).

This plan is governed by the standing discipline (`AGENT.md` §0.5/§0.6/§1.5): every gate is
**provably able to red** (named mutation, fires under `python -O`), its reference is
**structurally independent** of the SUT, and its regime **activates** the term the bug lives
in. Each gate row carries `config / reference / tolerance / claim-layer / mutation-that-reddens`.

---

## 0. The carve, restated structurally (so the gates target the right seam)

The SN scattering operator `S` (`orpheus/transport/operators/scattering.py`) decomposes,
TODAY, into:

| Channel | Math | Today's code path |
|---|---|---|
| ℓ=0 P0 in-scatter | `Σ_{s0}ᵀ @ φ` (`φ = ∫ψ dΩ`) | `add_iso_source` → `MaterialXSField.apply_p0_in_scatter` (per-material matmul) |
| (n,2n) | `2·Σ_{2n}ᵀ @ φ` | `add_n2n_source` → `MaterialXSField.apply_n2n` |
| ℓ≥1 aniso | `(1/W)·R∘Λ_{ℓ≥1}∘M ψ` | `build_aniso_source` = `(1/W)·kernel`, `kernel = frame.conjugate(Λ_{skip_l0=True})` |

`_assemble_per_ordinate_source` combines them: `(iso / Σw) + aniso`.

**A2a re-expresses the iso ℓ=0 P0 channel as the ℓ=0 frame conjugation.** Because ORPHEUS uses
*unnormalized* real harmonics (`Y₀⁰ = 1`): `integrate_angular(ψ)` IS the ℓ=0 analysis `M₀`, and
the iso `/W` broadcast IS the ℓ=0 reconstruction `R₀`. So the iso path equals
`(1/W)·R₀∘Σ_{s0}∘M₀` — exactly `frame.conjugate(Λ_{skip_l0=False})` restricted to ℓ=0. The full
Legendre scattering becomes `(1/W)·frame.conjugate(Λ_{ℓ≥0})` (`Λ` with `skip_l0=False`); (n,2n)
stays a **parallel ℓ=0 term** (separate XS `Σ_2n`, factor 2 — NOT folded into the frame).
The bespoke `add_iso_source`/`add_n2n_source` fast-path assembly retires (internal-only callers;
see §A2a-RETIRE). Forward-only; **no transpose yet**.

**A2b extends:** `LegendreMomentScattering.apply_transpose` (per-ℓ group-flip `Σ_{s,ℓ}ᵀ`) +
`CAP_APPLY_TRANSPOSE`; then `ScatteringOperator.apply_transpose = (1/W)·kernel.apply_transpose`
+ `CAP_APPLY_TRANSPOSE`. Because the iso ℓ=0 now lives IN the frame conjugation, `kernel`'s
Euclidean transpose (propagated by `OperatorProduct.apply_transpose`, operator.py:862-864)
covers **iso + aniso in ONE Sᵀ**. (n,2n) is the one channel OUTSIDE the frame → its transpose
is a separate gated term.

### Structural facts that pin the architecture (verified in-tree, not from memory)

1. **`OperatorProduct.apply_transpose` propagates iff BOTH operands have it**
   (operator.py:841-842, 862-864: `b.apply_transpose(a.apply_transpose(x))`). The frame faces
   `M`/`R` ALREADY advertise `CAP_APPLY_TRANSPOSE` (`frame.py:441/475`, `_AdjointOperator` gives
   `R.H`/`M.H` for free). **So the ONLY missing leaf is `LegendreMomentScattering.apply_transpose`**
   — `Λ` today is `{CAP_APPLY}` ONLY (scattering.py:283-285, the docstring at :262-263 explicitly
   defers it). The whole "S† falls out free" claim rests on this single leaf. ⇒ **A2b's load-
   bearing correctness gate is on `Λ.apply_transpose`** (everything downstream is propagation
   that the operator-algebra tests already cover).

2. **Euclidean transpose, NOT the metric Hilbert adjoint `.H`.** The carve spec wires
   `apply_transpose = (1/W)·kernel.apply_transpose` — the plain matvec transpose `Sᵀ`. The
   fission template (`test_fission_adjoint.py`) verifies exactly this: `op.apply_transpose` +
   **Euclidean** reciprocity `⟨Sφ,ψ*⟩=⟨φ,Sᵀψ*⟩` (Euclidean `.sum()`, no metric). `.H` is a
   DIFFERENT object (`_AdjointOperator.apply` = `G_V⁺ ⊙ apply_transpose(G_W ⊙ y)`, operator.py:676,
   carries the SH `(2ℓ+1)` Gram + `1/W` weighting). **A2b's reciprocity gate MUST use `apply_transpose`
   and a Euclidean inner product** — mixing in `.H`/a metric is a different identity and would
   either false-green or false-red. (Contrast the SEPARATE `test_g_adjoint_reciprocity.py`, which
   deliberately tests the `(L+C-B)` *G-adjoint* `.H` — that is the loss operator's metric adjoint,
   a different surface; do NOT conflate.)

3. **The iso fast-path helpers have NO external production caller** (grep: `add_iso_source` /
   `add_n2n_source` / `build_aniso_source` / `_assemble_per_ordinate_source` are referenced ONLY
   inside `scattering.py` + docstrings + tests). The within-group SI/Krylov driver calls
   **`S.apply(psi)`** (the `@singledispatchmethod` arms). ⇒ **A2a's behavioral change flows into
   the solve via `S.apply`** — the solve-level forward-preservation gates (keff/SI/MMS) genuinely
   reach the re-routed code (no Mode-11 "gate never executes the changed path" risk *at the solve
   level*). The retirement of the bespoke helpers is an INTERNAL refactor (see §A2a-RETIRE).

---

## A2a — forward gates (preservation + equivalence + perf)

### Gate A2a-EQ-1 — THE LOAD-BEARING EQUIVALENCE (the §0.5 named risk for A2a)

**Claim:** `(1/W)·frame.conjugate(Λ_{skip_l0=False}).apply(ψ.values)`'s **ℓ=0 contribution**
reproduces the legacy iso fast-path `(1/W)·add_iso_source(integrate_angular(ψ))`. Equivalently
(the cleaner, more direct probe): the NEW full-`S.apply(ψ)` equals the LEGACY-assembled
`(iso/Σw) + aniso` built from the retired-fast-path math.

- **Claim layer:** flux-shape (a per-ordinate source field, NOT an eigenvalue). Pillar:
  bit-identity / principled-equivalence (`vv §bit-identity`, criteria below).
- **Config:** ≥2G (anti-#3), **anisotropic** `Σ_{s,ℓ}` with ℓ≥1 non-zero so the FULL ℓ≥0
  tensor is exercised (NOT a P0-only fixture — that is the Mode-7/Mode-10 trap flagged in §RISKS),
  **heterogeneous** (≥2 materials with different asymmetric `Σ_{s0}` AND `Σ_{s1}`; L1 — flat-flux/
  homogeneous nulls cross-group + redistribution), `Sig2 ≠ 0` (so the n2n ℓ=0 term is live, #269).
  **Reuse `test_scattering_kernel_crosscheck.py::solver_p1_het`** (it already builds exactly this:
  P0+P1 asymmetric two-material 2G with `Sig2`). Drive with a non-isotropic `AngularFlux` (so
  ℓ≥1 moments are non-zero — `_aniso_psi`).
- **Reference (structurally independent of the NEW frame path):** the LEGACY assembly computed by
  a path that does NOT go through `frame.conjugate(Λ_{skip_l0=False})`:
  `legacy = (apply_p0_in_scatter(φ) + 2·apply_n2n(φ))/Σw + build_aniso_source(ψ)` where
  `φ = ψ.integrate_angular()`. The aniso half is the UNCHANGED `kernel` (skip_l0=True). The iso
  half is the per-material matmul `Σ_{s0}ᵀ@φ` — a DIFFERENT einsum order (per-material matmul vs
  the frame's `M₀→Λ₀→R₀` chain), so this is genuine procedural independence at the iso seam.

  *Pin the equivalence capture pre-carve:* snapshot `legacy` for `solver_p1_het` + `_aniso_psi`
  into a committed `.npy` (or compute it in-test from the pre-carve helpers if they survive as
  internal). Post-carve, `S.apply(ψ)` must reproduce it.
- **Tolerance — principled-equivalence, NOT 0-ULP.** `Y₀⁰=1` makes the math identical but the
  reduction tree differs (per-material `einsum("fg,fc...->gc...")` vs `R₀(Λ₀(M₀·ψ))` = three
  chained einsums + the `(2ℓ+1)` addition-theorem broadcast at ℓ=0 where `2·0+1=1`). Drift is
  FP-non-associativity, bounded by `(reduction_depth)·ULP` (`vv §bit-identity` criterion 3).
  Gate: `np.testing.assert_allclose(new, legacy, rtol=1e-14, atol=1e-300)` — i.e. ~tens of ULP,
  NOT `assert_array_equal`. Document all three §bit-identity criteria in the test docstring:
  (1) principled (the named intermediate is `M₀ψ = φ₀`, the scalar flux — a reactor-physics
  quantity); (2) structurally-independent reference (the per-material-matmul legacy half + the
  closed-form anchor in A2a-EQ-2); (3) drift = reduction-depth ULP.
- **Mutation that reddens:** (a) flip the ℓ=0 block sign in `Λ_{skip_l0=False}` → RED; (b) drop
  the ℓ=0 block (regress to `skip_l0=True` for the iso path) → RED (the iso source vanishes);
  (c) double-count ℓ=0 (leave the old `add_iso_source` call AND add ℓ=0 to the frame) → RED at
  2× iso. Confirm each fires under `-O` (function-call assertions only).

### Gate A2a-EQ-2 — closed-form iso anchor (the structural-independence ground for EQ-1)

EQ-1 alone proves "new == legacy" — necessary but NOT sufficient (`vv §bit-identity` criterion 2;
L2: old-vs-new ULP can't tell you the pre-carve value was *right*). Anchor the iso channel to a
closed form.

- **Claim layer:** flux-shape, closed-form pillar.
- **Config:** flat infinite-medium — **all-reflective, uniform `Σ_t`, uniform `Σ_{s0}` (cross-group
  ≥2G), uniform source `Q`** (no fission needed → use a non-fissile mixture, L7: an eigen solve on
  a moderator is `k=nan`; this is FIXED-SOURCE). Geometry: slab is fine here (the iso P0 channel
  has no curvilinear redistribution; this gate isolates the *energy-group* iso matmul, not angular).
  ≥2G mandatory.
- **Reference:** the closed-form flat-flux infinite-medium solution
  `φ = (diag(Σ_t) − Σ_{s0}ᵀ)⁻¹ Q` (L7; pure linear algebra on the XS matrices, ZERO transport
  discretization — structurally orthogonal to the sweep). At convergence the SN fixed-source
  solve must reproduce this (flat flux ⇒ the sweep is exact for the iso channel).
- **Tolerance:** `rtol = SAFETY(10) × conv_tol` read OFF the run config (L7 `assert_regression`
  pattern — the gate IS the claim, not a magic floor).
- **Mutation that reddens:** transpose `Σ_{s0}` the wrong way in the ℓ=0 frame block (`Σ_{s0}`
  vs `Σ_{s0}ᵀ`) → RED with the wrong group ratio (Mode-6 / Signature-3; the asymmetric cross-group
  `Σ_{s0}` makes the transpose detectable — verify the fixture's `Σ_{s0}` is asymmetric, else this
  gate is blind to the transpose, same precondition as the fission role-swap discriminator).

### Gate A2a-EQ-3 — RE-BASELINE audit on the EXISTING 0-ULP `test_scattering_kernel_crosscheck.py`

**The aniso ℓ≥1 `kernel` (skip_l0=True) is UNCHANGED by A2a** (A2a adds a *separate* full-ℓ frame
conjugation for the iso path; it does NOT touch `kernel`). Therefore:

- **`test_scattering_kernel_crosscheck.py` STAYS 0-ULP** (`assert_array_equal`). It pins
  `S.kernel.apply == _aniso_source_from_moment_values(M·ψ)` for the ℓ≥1 part — both sides
  unchanged. **Do NOT relax it to principled-equiv.** If A2a's implementation accidentally routes
  the aniso path through the new full-ℓ conjugation (changing the reduction order), THIS gate
  reds — which is the CORRECT signal that the implementer broke the unchanged-aniso invariant.
- **Action:** confirm post-carve this file is still 0-ULP green. If it reds, the carve changed the
  aniso path (a scope violation) — STOP. The NEW iso-through-frame equivalence is the SEPARATE
  principled-equiv gate A2a-EQ-1, NOT a relaxation of this one.
- **Sharpening (Mode-7 audit):** `solver_p1_het` is P1 + `Sig2≠0` + 2G het — it ALREADY activates
  ℓ≥1 and n2n. After A2a, ADD one assertion to this file (or A2a-EQ-1): the new full-`S.apply`
  ℓ=0 contribution is NON-ZERO for this fixture (`bool(np.any(iso_part != 0))`), so the equivalence
  gate is not vacuously comparing the iso channel as zeros (the Mode-10 "activated-but-unconstrained"
  guard, mirroring the existing `test_aniso_path_is_non_degenerate`).

### Gate A2a-FWD-1 — eigenvalue forward-preservation (the keff contract)

**Claim:** the converged keff is UNCHANGED to solver tolerance after the iso re-expression (the
iso channel is re-wired, not re-physics'd).

- **Claim layer:** EIGENVALUE (solver claim). Pillar: closed-form (homogeneous k_inf) +
  semi-analytical (heterogeneous Case reference).
- **Existing gates to re-run (named, do NOT re-baseline — they must stay GREEN unchanged):**
  - `tests/sn/verification/analytical/test_kinf_homogeneous.py` (l1) — **the canonical structurally-
    independent eigenvalue benchmark**: dominant eigenvector + eigenvalue of `A⁻¹F` by pure linear
    algebra, NO transport discretization (module docstring). 2eg + 4eg (anti-#3; the 4eg eigenvector
    catches the flux-shape/eigenvector-flip class a 1G k_inf is blind to). **This run exercises the
    iso P0 + cross-group P0 channel through the full SI solve** — exactly the A2a-touched term.
    *Caveat (L1):* homogeneous ⇒ flat flux ⇒ the *angular redistribution* is null here; this gate
    covers the ENERGY-iso channel, NOT angular. Angular coverage is FWD-3 (aniso MMS) + the
    heterogeneous keff below.
  - `tests/sn/eigenvalue/test_heterogeneous_transport.py` (l1, `@catches("ERR-025")`) — two-region
    reflective slab, piecewise `Σ`, mesh-refinement convergence to the Case singular-eigenfunction
    reference (mesh-independent, structurally independent). **Heterogeneous ⇒ activates cross-group
    redistribution** (L1/H2). The vv L2 multi-group heterogeneous mandatory row.
  - `tests/sn/eigenvalue/test_keff_slab.py`, `test_keff_2d.py`, `test_keff_curvilinear.py` — the
    broader keff regression set across geometries.
- **Tolerance:** the existing gates' own tolerances (they encode the claim). Pass = bit-stable to
  their committed thresholds. **A drift beyond FP-non-associativity here is a real regression** — the
  iso channel changed the physics, STOP.
- **Mutation that reddens:** (a) drop the ℓ=0 block in the new conjugation → keff moves O(1) (no
  in-scatter); (b) sign-flip ℓ=0 → keff diverges/moves O(1). Confirm `test_kinf_homogeneous` 4eg
  reds under each. (These mutations are the same as A2a-EQ-1 (a)/(b) but observed at the *eigenvalue*
  layer — the two layers cross-validate.)

### Gate A2a-FWD-2 — fixed-source SI/Krylov forward-preservation

**Claim:** the converged fixed-source flux is unchanged; the SI and Krylov fixed points still agree
(they solve the same operator — the iso re-expression must not split them).

- **Claim layer:** flux-shape (converged field). Pillar: closed-form (`Q/Σ_t` flat-flux) +
  cross-representation (SI ≡ Krylov).
- **Existing gates:**
  - `tests/sn/solve/test_fixed_source_g1.py` + `test_fixed_source_2d_equivalence.py` +
    `test_si_single_primitive_contract.py` — the SI fixed-source path.
  - **The single most powerful curvilinear iso diagnostic (L1, AGENT.md §0.6):** a fixed-source
    flat-flux `Q/Σ_t` case — uniform `Q`, uniform `Σ_t`, all-reflective → exact `φ = Q/Σ_t`
    everywhere, in **sphere AND cylinder** (the curvilinear redistribution is then activated by the
    angular part while the iso channel must still reproduce `Q/Σ_t`). Check this exists in
    `tests/sn/sweep/curvilinear/test_streaming_equilibrium_curvilinear.py`; if it does, re-run it;
    if the iso modernization touches the per-ordinate combine, this gate is the catch.
  - SI ≡ Krylov fixed-point agreement: `tests/sn/sweep/core/test_sweep_vs_apply_consistency.py`
    (the matvec-vs-sweep round-trip — both must still hit the same fixed point post-carve).
- **Config:** ≥2G heterogeneous for the redistribution-bearing row (vv L2); the `Q/Σ_t` anchor is
  flat by construction (its job is the iso-channel closed-form, not redistribution).
- **Tolerance:** `SAFETY × conv_tol` (L7); the existing gates' thresholds.
- **Mutation that reddens:** drop/sign-flip ℓ=0 in the conjugation → the SI fixed point moves O(1)
  AND (if the change is non-symmetric in the SI vs Krylov assembly) SI and Krylov split → both the
  `Q/Σ_t` and the consistency gate red.

### Gate A2a-FWD-3 — anisotropic ℓ≥1 redistribution preservation (the physics L1 reference)

**Claim:** the ℓ≥1 angular redistribution is unchanged (A2a does not touch `kernel`, but the
*combine* in `_assemble_per_ordinate_source` and the apply arms are touched — confirm the aniso
half still threads correctly).

- **Existing gate (structurally-independent physics reference — do NOT re-baseline):**
  `tests/sn/verification/mms/test_curvilinear_aniso_scattering_p1.py` (l1,
  `@verifies("pn-scatter","flux-moments")`) — sphere + cylinder P1 source vs a hand-derived
  reference. This is THE structural-independence ground for the ℓ≥1 scattering physics (the 0-ULP
  crosscheck explicitly cites it as its L1 backing). Plus `test_mms_aniso.py` /
  `test_curvilinear_aniso_convergence.py`.
- **Tolerance:** the gate's own (it stays green unchanged).
- **Mutation that reddens:** break the aniso threading in `_assemble_per_ordinate_source` (e.g.
  apply `/W` twice to the aniso half, or drop the aniso term) → the MMS source disagrees O(1).

### Gate A2a-PERF-1 — hot-path non-regression

**Claim:** the iso modernization does not slow the within-group hot path.

- **Config / what to measure (L16 wall-clock):** full `tests/sn/` (`python -O -m pytest tests/sn -q`,
  `-m "not slow"`) wall-clock, pre-carve vs post-carve, same machine, warm cache. A regression > ~10%
  is a flag.
- **The specific perf risks to confirm (the carve's own §3 questions):**
  1. **P0-only case (`scattering_order=0`, L=0):** the ℓ=0 moment tensor must be *just the scalar
     flux* — `M₀ψ = ∫ψ dΩ`, shape `(1,1,ng,*spatial)` — NOT a full `(L+1,2L+1,…)` allocation. Confirm
     the new conjugation at L=0 builds no larger tensor than the old `integrate_angular`. **Gate:**
     assert the ℓ=0-only moment field shape is `(1, 1, ng, *spatial)` (no wasted angular axes), and
     micro-time `S.apply(ScalarFlux)` / a P0 fixed-source solve pre/post — must not regress.
  2. **ℓ≥1 case:** the new path adds only the ℓ=0 BLOCK to an already-built `(L+1,2L+1,…)` tensor
     (the aniso path already allocates it). Confirm no SECOND full-tensor allocation (e.g. the
     implementer building iso and aniso tensors separately then summing). **Gate:** assert one
     moment-tensor allocation per `apply` (a peak-array or `tracemalloc` probe, or simply confirm
     the aniso and iso share the single `M.apply(ψ)` moment field).
- **Flag:** if the P0-only path regresses (the most common case — most production runs are P0), a
  **microbench is mandatory** before merge (the full-suite wall-clock can mask a P0 regression if
  the suite is aniso-heavy). Recommend `derivations/diagnostics/` microbench timing `S.apply` at
  L=0 across cell counts pre/post; promote to a perf characterization gate only if a regression
  appears (L5 — bound it one-sided, no upper calcification).
- **Mutation that reddens:** N/A (perf gate); the teeth are the pre/post comparison itself. Document
  the baseline number in the plan/commit so a future regression is detectable.

### §A2a-RETIRE — retirement audit (coding-standards: retire-as-you-go + test migration)

`add_iso_source` / `add_n2n_source` retire from the assembly path. The THREE-search blast radius:

1. **Graph callers** (`nexus impact`): internal to `scattering.py` only (verified by grep — no
   external production caller).
2. **Text-grep across code + tests + docs:** `test_scattering_operator.py::TestBitIdenticalExtractionP0`
   (4 tests) pins `add_iso_source`/`add_n2n_source` against a per-cell reference at `rtol=1e-13`.
   **Test-migration decision (per coding-standards):**
   - These are **behavioral** tests (the iso P0 in-scatter math `Σ_{s0}ᵀ@φ` is a correctness
     contract) → **rewire to the successor**, NOT delete. The successor surface is EITHER the
     surviving `MaterialXSField.apply_p0_in_scatter` kernel (which the frame's ℓ=0 reconstruction
     still calls per-material) OR the typed `S.apply(ScalarFlux)` arm (which today routes through
     `add_iso_source`+`add_n2n_source`, lines 1368-1369 — confirm A2a keeps `S.apply(ScalarFlux)`
     producing the iso scalar source; the W-F docstring marks it a kept cross-method surface).
   - The `*_delegators_match_operator_directly` tests pin `SNSolver._add_scattering_source`
     delegation — if those solver delegators also retire, migrate or delete per their nature
     (API-smoke → delete; behavioral → rewire).
   - **Verify migration is non-vacuous:** the migrated assertion must still red on a transposed/
     dropped `Σ_{s0}` (re-confirm the mutation post-migration; a migrated test that lost its teeth
     is a silent coverage loss).
3. **Direct constructors / docs:** the `material_xs_field.py` docstrings (:10-11, :623, :659) cite
   `add_iso_source`/`add_n2n_source` by name — update the prose pointers (an unresolved `:meth:`
   cross-ref renders as plain text with NO `-W` warning, so the Sphinx gate won't catch it). The
   `kernel`/`build_aniso_source` docstrings at :601, :680, :690-691, :719, :732 describe the iso as
   the "separate `add_iso_source` fast path" — these become **stale** (iso is now in the frame).
   Update them or the next session mints a wrong mental model (Cardinal Rule 3).

**The (n,2n) channel does NOT retire** — it stays a parallel ℓ=0 term (`apply_n2n` / the `Σ_2n`,
factor-2 math). Only the *fast-path assembly framing* of the P0 ℓ=0 iso retires.

---

## A2b — transpose gates (S† correctness)

Mirror the A1 fission template `tests/sn/operators/test_fission_adjoint.py` verbatim in structure
(foundation-tagged, `-O`-safe `require`/`np.testing.*`, ≥2G AND 4G). New file:
`tests/sn/operators/test_scattering_adjoint.py`.

### Gate A2b-Λ-1 — `LegendreMomentScattering.apply_transpose` correctness (THE load-bearing leaf)

**Claim:** `Λᵀ` applies the per-ℓ group-TRANSPOSE `Σ_{s,ℓ}ᵀ` on moment space (the forward `Λ`
applies `Σ_{s,ℓ}`; the transpose flips the group axis per ℓ block).

- **Claim layer:** flux-shape (a moment field). Pillar: structurally-independent hand loop.
- **Config:** ≥2G AND 4G; **asymmetric `Σ_{s,ℓ}` per ℓ** (per-group non-proportional, so the
  `[g_from,g_to]` transpose is detectable — the same asymmetry precondition as the fission role-swap).
  P1 minimum (ℓ=0 + ℓ=1 blocks both transposed; with `skip_l0` BOTH settings tested). Heterogeneous
  (≥2 materials).
- **Reference (structurally independent):** an explicit Python double-loop computing
  `(Λᵀφ)_ℓ^m|_g = Σ_{g'} Σ_{s,ℓ}[g, g'] φ_ℓ^m|_{g'}` (note the index order — forward is
  `Σ_{s,ℓ}[g',g]`, transpose is `Σ_{s,ℓ}[g,g']`), per ℓ block, per cell — sharing NO numpy reduction
  with `MaterialXSField.apply_legendre_scattering_moments`. (Add a `hand_derived_legendre_scattering`
  helper to `tests/transport/_integral_kernel_helpers.py`, mirroring `hand_derived_fission_emission`.)
- **Tolerance:** `rtol=1e-13` (procedural-independence, not 0-ULP — the hand loop ≠ einsum order).
- **Mutation that reddens (THE group-transpose mutation):** in `Λ.apply_transpose`, use the
  forward index `Σ_{s,ℓ}[g',g]` instead of `Σ_{s,ℓ}[g,g']` (i.e. forget to transpose) → RED. With
  asymmetric `Σ_{s,ℓ}` the two differ; confirm via a discriminator row (mirror
  `TestRoleSwapDiscriminator`): `hand(transposed) != hand(forward)` for this fixture, else the gate
  is blind to the flip.

### Gate A2b-Λ-2 — `Λ` capability advertises `apply_transpose`

- **Claim:** `CAP_APPLY_TRANSPOSE in LegendreMomentScattering(...).capabilities`. Update the EXISTING
  `test_legendre_moment_scattering.py::test_apply_only` (lines 80-83) which currently ASSERTS
  `CAP_APPLY_TRANSPOSE not in capabilities` — **this assertion FLIPS** (it is the explicit "deferred"
  pin; A2b lands the feature, so the negative assertion becomes a positive one). This is a planned
  re-baseline of a structural pin (NOT a value gate).
- **Mutation that reddens:** ship `Λ.apply_transpose` but forget to add `CAP_APPLY_TRANSPOSE` to the
  `capabilities` frozenset → the propagation through `OperatorProduct.apply_transpose` SILENTLY
  drops the transpose capability (operator.py:841: `_has(a, CAP_APPLY_TRANSPOSE) and _has(b, …)`),
  so `S.apply_transpose` / `S.kernel.apply_transpose` would raise `MissingCapability` → A2b-S-2 reds.
  This gate catches the capability-set/implementation mismatch directly.

### Gate A2b-S-1 — `S†` correctness vs a structurally-independent dense oracle

**Claim:** `S.apply_transpose(ψ*) = (1/W)·kernel.apply_transpose(ψ*)` is the genuine Euclidean
transpose of the full forward `S.apply` (now covering iso + aniso, since iso is in the frame).

- **Claim layer:** flux-shape. Pillar: structurally-independent dense oracle.
- **Config:** ≥2G AND 4G, **anisotropic (P1) + heterogeneous + `Sig2≠0`** (reuse `solver_p1_het`).
  This config exercises iso ℓ=0 + aniso ℓ≥1 + n2n in ONE Sᵀ (the carve's claim that iso-in-frame
  unifies the transpose coverage).
- **Reference (dense `Sᵀ` oracle, structurally independent of the matvec):** materialize the forward
  `S` as a dense matrix by applying `S.apply` to each Euclidean basis vector `e_i` (column-by-column),
  transpose the matrix, and apply to `ψ*`. This is `Sᵀψ*` computed WITHOUT calling
  `apply_transpose` — it shares NO transpose primitive with the SUT (it's the *definition* of the
  matrix transpose via forward applies). For a small fixture (4×3 cells, 2G, ~lebedev N) the dense
  matrix is tractable. (Mirror the `diag_p42_adjoint_oracle.py::validate_composite_adjoint` dense-
  probe pattern cited in `test_g_adjoint_reciprocity.py`.)
- **Tolerance:** `rtol=1e-12` (the dense-build accumulates a few more FP ops).
- **Mutation that reddens:** (a) `Λ` not transposed (A2b-Λ-1's mutation, observed at the S level);
  (b) `S.apply_transpose` forgets the `1/W` → RED by the `Σw` factor; (c) transpose only aniso,
  drop the iso ℓ=0 from the transpose (regress `kernel` to `skip_l0=True` for the transpose path) →
  RED on the iso block.

### Gate A2b-S-2 — Euclidean reciprocity `⟨Sφ,ψ*⟩ = ⟨φ,Sᵀψ*⟩` (the transpose-DEFINING identity)

**Claim:** the defining adjoint property, **per-group, Euclidean** (NOT weight-summed — L27).

- **Claim layer:** flux-shape (an identity). Pillar: the transpose definition itself.
- **Config:** ≥2G AND 4G, **asymmetric `SigS` per ℓ** (a symmetric `SigS` makes `S=Sᵀ` and the gate
  vacuous — Signature-3/L27; the discriminator below pins this), **anisotropic P1 + heterogeneous +
  `Sig2≠0`** (reuse `solver_p1_het`; confirm its `SigS`/`Sig2` are asymmetric — `solver_p1_het`'s
  `p0_a=[[.38,.10],[.05,.90]]` is asymmetric ✓, `Sig2=[[0,.03],[.01,0]]` asymmetric ✓). **Now that
  iso is in the frame, this ONE reciprocity gate covers iso + aniso + n2n in a single `Sᵀ`** (the
  carve's headline simplification — state it in the docstring).
- **Reference:** the identity is self-referential between `apply` and `apply_transpose` — make it
  non-vacuous via the discriminator (next gate). `lhs = (S.apply(φ).values * ψ*).sum()` per the
  Euclidean inner product; `rhs = (φ * S.apply_transpose(ψ*)).sum()`. (Both φ and ψ* are
  asymmetric random fields per group AND cell — `_asymmetric_field`.)
- **Tolerance:** `rtol=1e-12` (the fission template's reciprocity threshold).
- **Mutation that reddens:** any wrong `Sᵀ` (wrong axis, missing `1/W`, un-transposed `Λ`, dropped
  iso/n2n block) breaks the identity → RED. (Reciprocity is the broadest single transpose gate —
  it catches the whole class; A2b-S-1's dense oracle pins the *value*, this pins the *adjoint law*.)

### Gate A2b-S-3 — reciprocity discriminator (the gate genuinely constrains `Sᵀ ≠ S`)

**Claim (Mode-7/Mode-2 precondition):** with asymmetric `SigS`, `S ≠ Sᵀ`, so the reciprocity gate
is not vacuously satisfied by `S` being symmetric.

- **Reference:** `require(not np.allclose(S.apply(x).values, dense_S_transpose @ x), …)` — i.e. the
  forward and the (independently-built dense) transpose genuinely DIFFER for this fixture. If they
  agree, the fixture lost its asymmetry and A2b-S-2 is blind to a `Sᵀ`-that-equals-`S` bug.
- **Mutation that reddens:** symmetrize the fixture's `SigS` → this gate fires (telling you the
  reciprocity gate above would be vacuous). Mirrors `TestRoleSwapDiscriminator` exactly.

### Gate A2b-S-4 — capability + Mode-11 routing sentinel

**Claim:** `S` advertises `apply_transpose`, and `S.apply_transpose` ROUTES THROUGH
`kernel.apply_transpose` (→ `Λ.apply_transpose` via the `OperatorProduct` chain), NOT an inline
reduction.

- **Capability:** `require(CAP_APPLY_TRANSPOSE in S.capabilities, …)` — and confirm `S` (the
  `@dataclass`) gains `CAP_APPLY_TRANSPOSE` in its `capabilities` (today scattering.py advertises
  `{CAP_APPLY}` only — A2b must add it).
- **Mode-11 routing (the strictly-stronger in-process wrap, per `vv` Mode-11 sharpening):**
  monkeypatch `LegendreMomentScattering.apply_transpose` (the NEW leaf) with a counting wrapper;
  call `S.apply_transpose(ψ*)`; assert the counter > 0. A `S.apply_transpose` that inlines the
  group-flip instead of reusing `Λ.apply_transpose` leaves the counter at 0 → RED. (Mirror
  `test_apply_transpose_routes_through_kernel_transpose` which wraps `TensorProductOperator.apply_transpose`
  for fission; here the wrapped target is `LegendreMomentScattering.apply_transpose`.)
- **Mutation that reddens:** implement `S.apply_transpose` as a bespoke per-material `Σ_{s,ℓ}ᵀ` loop
  bypassing `kernel` → the sentinel counter stays 0 → RED.

### Gate A2b-N2N-1 — (n,2n) transpose (the one channel OUTSIDE the frame)

**Claim:** the (n,2n) ℓ=0 term's transpose is `2·Σ_{2n}ᵀ` (the `Σ_2n` group-flip + factor 2),
gated separately because n2n is NOT in the frame conjugation (it's a parallel ℓ=0 term, scattering.py
comment :1066-1068).

- **Claim layer:** flux-shape. Pillar: structurally-independent hand loop.
- **Config:** **non-zero `Sig2`** (per #269 the n2n fixtures MUST be non-vacuous — `make_mixture`
  defaults `Sig2=0` AND every A/B/C/D library mixture has `Sig2=0`, per lessons L1; build `Sig2`
  DIRECTLY on the mixture as `solver_p1_het` does, or reuse `_make_mixture_with_n2n`). Asymmetric
  `Σ_2n` (so the transpose is detectable). ≥2G.
- **Reference:** explicit hand loop `(n2nᵀ source)_g = 2·Σ_{g'} Σ_{2n}[g,g'] φ_{g'}` (transposed
  index), structurally independent of `apply_n2n`.
- **Tolerance:** `rtol=1e-13`.
- **Mutation that reddens:** (a) drop the factor 2 in the n2n transpose → RED at ½; (b) forget to
  transpose `Σ_2n` → RED (asymmetric `Σ_2n` makes it detectable); (c) **omit the n2n channel from
  `S.apply_transpose` entirely** (only transpose the frame conjugation, forgetting n2n is outside it)
  → RED. (c) is the highest-value mutation — it directly tests the carve's "n2n is the one separate
  ℓ=0 channel" claim; confirm A2b-S-2's reciprocity gate (with `Sig2≠0`) ALSO catches it (the
  reciprocity identity must include the n2n contribution on both sides).

---

## §RISKS — Mode-7 / Mode-10 / Mode-11 hazards specific to the iso-modernization

These are the §0.5/§0.6 traps where a green gate would NOT catch the bug. The implementer/qa MUST
confirm each is closed.

- **R-Mode-7 (the P0-only fixture nulls the unified ℓ≥0 path) — THE central risk.** A2a's whole
  point is the FULL ℓ≥0 conjugation. A gate that runs on `scattering_order=0` (P0-only) NEVER builds
  the `ℓ≥1` blocks, so it cannot see a bug in "adding the ℓ=0 block to the already-built tensor"
  (PERF-1.2) NOR a bug where the ℓ=0 frame block is wrong ONLY when ℓ≥1 is present. **Mitigation:**
  EVERY equivalence/correctness gate (A2a-EQ-1, A2a-EQ-2 except its flat closed-form anchor, all
  A2b) runs on the **anisotropic P1 + heterogeneous + Sig2≠0** `solver_p1_het` fixture. The P0-only
  fixture (`solver_2g_p0`) is used ONLY for the perf P0-path gate (PERF-1.1) and the
  `S.apply(ScalarFlux)` iso-scalar arm — NEVER as the sole equivalence reference. Flagged because
  `test_scattering_operator.py` has many `solver_2g_p0` rows; a careless implementer might add the
  A2a equivalence row there and it would pass vacuously.

- **R-Mode-10 (the iso ℓ=0 contribution is sub-dominant in a redistribution-heavy config).** If the
  equivalence gate measures only the FULL `S.apply` and the iso ℓ=0 is a small fraction of the aniso-
  dominated source, a sign error in ℓ=0 could sit below the `rtol`. **Mitigation:** A2a-EQ-3's
  non-degeneracy assertion (`iso_part != 0` AND comparable in magnitude); the A2a-EQ-2 closed-form
  anchor isolates the iso channel with NO aniso present (iso is O(1) there); and the A2a-FWD-1 keff
  mutation observes the ℓ=0 drop at the eigenvalue layer where it is O(1).

- **R-Mode-11 (a transpose gate that never executes `Λ.apply_transpose`).** The dense-oracle gate
  (A2b-S-1) builds `Sᵀ` from FORWARD applies — it does NOT call `apply_transpose` at all (that's its
  structural-independence virtue) — so it does NOT prove the SUT's `apply_transpose` is even reached.
  **Mitigation:** A2b-S-4's in-process wrap of `LegendreMomentScattering.apply_transpose` is the
  Mode-11 catcher — it proves the production `S.apply_transpose` enters the new leaf. The dense oracle
  proves the VALUE; the wrap proves the PATH. Both are required.

- **R-double-`/W` / convention drift (Pattern-7 crosswalk).** The iso path's `/W` lived in
  `_assemble_per_ordinate_source` (`(iso/Σw) + aniso`). The frame conjugation's `R₀` is the iso
  broadcast — does it ALSO carry `/W`, risking a double-`/W`? **Mitigation:** A2a-EQ-1's principled-
  equiv gate (`new == legacy` at 1e-14) catches a double-`/W` directly (it would be off by `Σw`),
  and the `Q/Σ_t` flat-flux anchor (FWD-2) catches a `/W` convention error O(1). The implementer MUST
  write the convention crosswalk (per-ordinate ↔ moment, where `/W` lives) into the carve plan BEFORE
  coding (coding-elegance §crosswalk; the carve crosses the per-ordinate ↔ moment ↔ fast-path seams).

---

## Summary table — re-baseline disposition (EXPLICIT per §1.5)

| Existing gate | A2a/A2b disposition | Tolerance after carve |
|---|---|---|
| `test_scattering_kernel_crosscheck.py` (ℓ≥1 `kernel`) | **STAYS 0-ULP** (aniso unchanged) | `assert_array_equal` (unchanged) |
| `test_scattering_operator.py::TestBitIdenticalExtractionP0` (iso fast-path) | **MIGRATE** to successor surface (retirement) | `rtol=1e-13`, teeth re-confirmed |
| `test_legendre_moment_scattering.py::test_apply_only` (`CAP_APPLY_TRANSPOSE not in`) | **FLIP** (negative pin → positive; A2b lands feature) | structural (capability set) |
| `test_kinf_homogeneous.py` (4eg) | **STAYS GREEN unchanged** (forward-preservation) | gate's own threshold |
| `test_heterogeneous_transport.py` (`@catches ERR-025`) | **STAYS GREEN unchanged** | gate's own threshold |
| `test_curvilinear_aniso_scattering_p1.py` (`@verifies pn-scatter`) | **STAYS GREEN unchanged** | gate's own threshold |
| **NEW** `test_scattering_adjoint.py` (A2b-Λ/S/N2N) | **NET-NEW** (mirror fission template) | `rtol=1e-12..1e-13` |
| **NEW** A2a-EQ-1/EQ-2 equivalence rows | **NET-NEW** (principled-equiv + closed-form) | `rtol=1e-14` / `SAFETY×conv_tol` |

**Net-new files:** `tests/sn/operators/test_scattering_adjoint.py` (A2b), the A2a equivalence rows
(extend `test_scattering_kernel_crosscheck.py` or a new `test_scattering_iso_frame_equivalence.py`),
a `hand_derived_legendre_scattering` + `hand_derived_n2n` helper in
`tests/transport/_integral_kernel_helpers.py`.

**Catalog:** if A2a/A2b catches a real bug during implementation (e.g. the iso ℓ=0 frame block
transpose, the n2n-outside-frame omission), log an ERR-NNN per the `vv-principles` "Log every
caught bug" directive — the iso-through-frame double-`/W` and the n2n-omission are new failure-mode
candidates not yet in the catalog.
