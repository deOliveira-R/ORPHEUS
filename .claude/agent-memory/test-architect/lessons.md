# Test Architect — Lessons

Read at the START of every invocation. Behavioral corrections only:
"what test-design mistake did I make / catch, and what verification
discipline fixed it?" The failure-mode TAXONOMY (Modes 1–11, the three
pillars, the anti-patterns) lives in the preloaded `vv-principles`
skill — NEVER restate it here; cite it (`Mode N`, `anti-#N`, `§pillar`).
The reference inventory + cross-section mixtures live in `AGENT.md`.
Campaign play-by-play is gone (merged archaeology).

The spine behind nearly all of these: **a verification plan is not done
when "the tests pass" — it is done when (a) each gate is provably ABLE
to red (a mutation reddens it, under the canonical `-O`), (b) the
reference is structurally INDEPENDENT of the SUT, and (c) the test's
regime ACTIVATES the term the bug lives in.** Every lesson is one face
of that standard. → This spine is now a STANDING directive in `AGENT.md`
§0.5 (gate-liveness) + §0.6 (config-blindness) + §1.5 (structural
independence) — it governs every plan by default. The L1/L4 entries
below are kept for the forensic war-stories that grounded it.

---

## L1 — Homogeneous / 1-group / flat-flux tests are blind to the hardest math

→ **Promoted to `AGENT.md` §0.6 as a standing directive.** Kept here for
the forensic war-story (the cyl DD suite that documented the `Q/Σ_t`
diagnostic and STILL omitted it, so the bug hid).

The recurring root failure: the convenient test config nulls the exact
term the solver is most likely to get wrong. The blindnesses compound —
know all of them and pick a config that breaks every one that matters:

- **Flat flux** nulls every redistribution / α-recursion / weight-cancel
  term (curvilinear closure, angular redistribution). `vv §H2`.
- **1-group** makes `k = νΣ_f/Σ_a` flux-shape-independent — degenerate.
  Always ≥2G for an eigenvalue claim (`vv §1-group`, anti-#3).
- **Homogeneous** nulls redistribution AND spatial-distribution bugs.
- **Slab geometry** is the degenerate SOTP case (angular redistribution
  is a ZeroOperator) — a slab-only separability/closure test PROVES
  NOTHING about the coupled curvilinear case.
- **Isotropic-source 2-D snapshots** are blind to a dropped φ_ℓ≥1: a
  windowing that silently drops higher moments PASSES an all-isotropic
  suite. Before carving a moment-reduction path, AUDIT whether the
  existing snapshots exercise the moments being reduced; if all are
  isotropic, manufacture an anisotropic case FIRST.
- **A SYNTHETIC fixture can null a property the REAL data exercises — making
  a synthetic-only assertion a brittle trap that FALSE-REDs on real data.**
  (#274 F7.) An old leg condensed a *synthetic* 421-group χ (a literal step
  function in the fastest 10 % of groups) → WIMS-69 and asserted
  `argmax(χ_coarse)==0`. The REAL `pwr_like_mix()` χ peaks at group 61
  (~1.15 MeV) — the 2e7-eV grid ceiling carries ~60 groups of
  small-but-nonzero high-energy tail ABOVE the fission peak, so the real
  χ-peak lands in coarse group ~4, NOT 0; `argmax==0` would FALSE-RED on the
  production library. The synthetic fixture had nulled the high-energy tail
  the real data has, so the brittle assertion only "worked" on the synthetic
  case. **Fix = assert the physically-correct CUMULATIVE-MASS property
  (`χ[:n//2].sum() > 0.5`, the bulk is in the fast half), never a brittle
  argmax/exact-index on a value whose real distribution differs from the
  fixture's.** Trap signature: a `np.argmax(...)==K` / exact-index assertion
  whose K was read off a *synthetic* fixture. Before pinning an index/peak on
  synthetic data, ask "does the REAL data put the peak HERE?" — if the real
  distribution is broader/shifted, pin a cumulative/inequality property
  instead. (Same family as the config-blindness rule, one level up: the
  *fixture*, not the *config*, is unrepresentative — verify the gate against
  REAL data when it exists, which is exactly why F7 replaced the
  synthesize-and-skip leg with `pwr_like_mix()`.)
- **`make_mixture` SILENTLY NULLS two channels — `Sig2` AND `SigL`.**
  (a) `sig_2` defaults to `csr_matrix((ng,ng))` = all-zero, AND it is zero
  on EVERY `xs_library` A/B/C/D fixture; (b) **`make_mixture` has NO
  `sig_l` parameter — it hardcodes `SigL = np.zeros(ng)`
  (`xs_library.py:56`)**, so any `sig_l` you "pass" is dropped on the
  floor. So any (n,2n) term (SN `compute_group_production_rate`'s
  `+2·∫Σ₂·φ`, the CP/MoC/homog net-removal `−2·Σ₂`) and any (n,α)/SigL term
  of the balance identity are identically NULLED (Mode-10) on a
  `make_mixture` fixture; regression tests that touch them go VACUOUSLY
  green. **Trap signature**: a "balanced" fixture built as `sig_t = sig_c +
  sig_l + sig_f + rowsum` THROUGH `make_mixture` is IMBALANCED by exactly
  `sig_l` (the SigT carries it, the Mixture zeros it → `assert_balanced`
  raises). To exercise SigL/Sig2 you MUST build the `Mixture(...)`
  constructor DIRECTLY (the `test_mixture_xs_balance.py` / P5.0
  `test_mixture_condense.py::_balanced_fissile_4g` hand-built pattern), or
  `make_mixture(..., sig_2≠0)` for n2n only (reuse
  `tests/cp/test_verification.py::_make_mixture_with_n2n`). Same shape as
  the anisotropic-case rule above: the convenient builder nulls the
  channel; manufacture the activating fixture FIRST.

How to apply: minimum catch for a curvilinear solver = heterogeneous
spatial-convergence (keff differences shrink under refinement) +
fixed-source flat-flux `Q/Σ_t` (the single most powerful curvilinear
diagnostic — it was documented and STILL omitted from the cyl DD suite,
and the bug hid). Print the analytical-references checklist from
`AGENT.md` and check each item off against a test; do not rely on memory.

---

## L2 — Match the reference to the CLAIM LAYER, and terminate in a structurally-independent ground

The two recurring reference mistakes: using a reference that cannot
prove the claim's layer, and using one that shares structure with the SUT.

- **MMS does NOT prove eigenvalues** (source-driven by construction).
  An eigenvalue claim needs a closed-form (homogeneous k_inf, transfer
  matrix) or semi-analytical reference. MMS reaches only flux-shape +
  convergence-order. Declare each gate's claim layer (convergence-order
  / flux-shape / eigenvalue) and confirm the pillar can prove it
  (`vv §pillars`).
- **Convergence rate is necessary, NEVER sufficient.** O(h²) to the
  wrong limit is still O(h²). A self-referencing Richardson ladder (ref
  built from the same buggy code) converges cleanly to the wrong value.
  ALWAYS pair a rate assertion with a converged-VALUE assertion against
  a structurally-independent reference (anti-#5).
- **No mesh-independent transport eigenvalue reference exists** in this
  codebase (heterogeneous refs are diffusion-based ~0.3% gap or
  self-referencing). Use diffusion-eigenvalue as a cross-check with an
  explicit tolerance, NEVER a precision target. Tracked: Issue #8.
- **Two-anchor template for a pure-refactor matvec carve**: a committed
  snapshot ("didn't move" = bit-id inheritance) is necessary-NOT-
  sufficient (ULP-distance can't tell you the pre-carve value was right)
  — pair it with a closed-form value anchor (`Q/Σ_t`, k_inf). The Krylov
  path is a GENUINE cross-check of the SI path (different representation,
  same fixed point) only when one side is the UN-windowed/UN-carved path.

How to apply: write the (claim-layer, pillar, truth-source) triple per
gate BEFORE drafting it. Force the structural-independence cross-check
from a DIFFERENT angle (a kernel/row-sum check AND a closed-form check),
never two derivations of the same closed form (`vv §structural-indep`).

---

## L3 — Override the MMS simplification bias at write-time; declare what the ansatz ACTIVATES vs NULLS

The "simplest trig that satisfies the BCs" is the default human choice
because hand-derivation is error-prone — and it is exactly the ansatz
that tests nothing (`vv Mode 7`). Two concrete instances that bit:

- **Every prior SN MMS ansatz VANISHES at the boundary** → the
  prescribed-inflow `q.boundary≠0` path was NEVER exercised. The fix
  needed a NON-vanishing-at-face ansatz: `ψ=(A(x)+μ·B(x))/W` with
  `A=a0+a1·sin(kx)`, **a0>0 load-bearing** (makes the incoming current
  γ₋ψ≠0). A sin-only ansatz re-creates the very gap the MMS exists to
  close.
- **An isotropic-in-μ ansatz nulls the angular redistribution** — the
  sweep's hardest math. A curvilinear/Pℓ MMS MUST pair with an
  angularly-non-trivial companion (`ψ=(A(r)+B(r)μ)/W`, B≠0 so ∂ψ/∂μ≠0).
  NEVER ship slab-only or isotropic-only.

For every MMS gate, declare in the docstring which terms the ansatz
**activates** AND which it **nulls** (`vv Mode 7`). If the nulled set
includes a term covered by an active ERR-NNN, redesign the ansatz. When
the test-architect spec mandates a stress ansatz (mixed-scale k, ≥2G
het, a0>0), that is a BINDING CONTRACT — a shipped test that swaps in
the canonical sin/1G/homogeneous case is a gate downgrade, not a style
choice (caught downstream by qa).

How to apply: with SymPy the source is derived programmatically, so
there is NO derivation-cost excuse for the simple ansatz. Pick ψ_chosen
for stress value (high frequency, mixed scales, group coupling,
non-vanishing boundary) and write the activate/null declaration first.

---

## L4 — A gate that cannot red is worse than no gate: prove every gate's teeth bite

→ **Promoted to `AGENT.md` §0.5 as a standing directive** (the gate must
be provably able to red, under `-O`). Kept here for the three concrete
mechanisms + the forensic detail.

The deepest recurring test-design failure: shipping a green gate that is
structurally incapable of catching the bug it claims to catch. Three
distinct mechanisms, all requiring an explicit liveness proof:

- **Mode-8 (`-O` strips bare `assert`)**: the canonical invocation is
  `python -O`. A bare `assert` in a PRODUCTION/helper module (or any
  non-collected module) is a NO-OP under `-O` — an always-on canary that
  cannot die is a FALSE green. Use `np.testing.assert_*` / `pytest.fail`
  / explicit `raise` (function calls, fire under `-O`). Bare asserts in
  COLLECTED `tests/` modules ARE rewritten by pytest and DO fire under
  `-O` (verified) — reserve the Mode-8 flag for `orpheus/` paths and
  always-on sentinels. DEMONSTRATE: a deliberately-wrong value must still
  FAIL under `-O`.
- **The xfail-strict false-positive**: a `strict=True` xfail is
  "satisfied" by failing for ANY reason — a stale array index makes it a
  FALSE xfail (green suite, wrong reason). For features-not-yet-landed
  use `strict=False` + a `reason=` pointing at the unlocking API/label,
  so it naturally un-xfails (flips to xpass) when the feature lands.
  Prove an xfail→live FLIP is non-vacuous: confirm it is RED pre-fix for
  the DOCUMENTED reason and GREEN post-fix.
- **Mutation-validate the catcher**: a `catches(ERR-NNN)` is a COVERAGE
  CLAIM — re-introduce the EXACT documented bug and confirm THIS test
  (not merely some test in the run) reddens. Same-area/same-module is
  NOT coverage. For a self-consistency round-trip (`residual(kernel(q))≈0`,
  both arms share the kernel), instrument a call-COUNT to prove the path
  is even exercised end-to-end — a round-trip pins neither correctness
  nor convention (a matvec runs only under `inner_solver='krylov'`; SI
  sweeps never touch it).

How to apply: for every gate named as evidence, ask "what mutation
reddens this, and does it fire under `-O`?" before crediting it.
Mutation-by-monkeypatch in-process (NEVER `git checkout` a file you have
uncommitted edits in). See qa L-006/L-010 (Mode-8 classification),
`vv Mode 8/11`.

---

## L5 — A characterization gate pins what is TRUE + the regression floor WITHOUT calcifying the limitation

When pinning a documented WONTFIX limitation (a first-order-at-pole
closure, an interpolation floor, an approximation plateau), the gate
must protect the floor without blocking a legitimate future fix:

- **GUARANTEE tests** carry `verifies(...)` and assert the property that
  IS correct (global L2 O(h²)).
- **CHARACTERIZATION tests** carry NO `verifies(...)` and bound the
  limitation ONE-SIDED (pole L∞ order `> 0.8` = "≥ first-order, does not
  regress", NO upper bound → a future O(h²) fix keeps the gate GREEN).
- When a fix claims to "remove a floor", do NOT assert removal — measure
  the floor's SCALING with the OTHER axis (quadrature). If it scales the
  same pre+post, the fix did NOT touch that floor; pin the floor as a
  verified scaling quantity (`err(S32) < err(S16)/2` is FALSIFIABLE — a
  fixed closure-bug floor is quadrature-independent → ratio ~1 → fails).
  `vv anti-#5/#17`.

How to apply: separate the guarantee from the characterization; bound a
limitation below-only; pin a floor by its scaling, never assert its
absence. Promote a diagnostic as a NEGATIVE-regression foundation gate
(no `verifies`) where its FAILURE is the signal to update the gate.

---

## L6 — A Mode-10 sub-floor term is closed by STRUCTURAL teeth, not a tightened value band

When a term's code path RUNS but its error is O(h²)-small (below the
convergence floor — a slope source, a boundary-trace transverse slope),
tightening the converged-flux value band does NOT pin it (the error is
sub-floor by construction). The complete closure is TWO O(1) structural
teeth + a control:

1. The production producer threads the projected moment through at
   MACHINE PRECISION (catches a regression to zeroing), with a
   structurally-independent (leggauss-only) reference.
2. A CONSUMED source-row sign flip moves the converged answer ≫ solver
   tol (catches sign-blindness), PAIRED with the no-op control (scalar /
   zeroed input → byte-identical) that pins the asymmetry.

A SHARPER case has NO value-improvement leg AT ALL: a boundary-trace
slope is sub-floor for ANY value claim everywhere (no regime makes it
the dominant forcing — a correctly-consumed slope can make the converged
value slightly WORSE). Do NOT manufacture a value-improvement gate there
— it would falsely RED a correct term.

How to apply: when a term is localized / higher-order-small, ship
producer-threading-at-machine-precision + consumed-flip-≫-tol + a no-op
control; accept the absence of a value-band gate as the CORRECT
signature. Mode-7's activated/nulled declaration is necessary but not
sufficient — mutation-verify an activated term is also CONSTRAINED.
`vv Mode 10`; mirrors method-implementer L-007 / qa L-037.

---

## L7 — Regression-snapshot tolerance is a CLAIM about what the solver promised, not a magic floor

Hand-picked floors (`rtol=1e-12`, `rtol=5e-6`) encode no claim — they
were either too tight for FP-non-associativity or papered over a real
bug. Replace with a single-source-of-truth `assert_regression`:

- **Iterative** result (k_eff/flux from power/source iteration) →
  `SAFETY × conv_tol`, where conv_tol is the solver's OWN stopping
  criterion read OFF the run config (the SoT shared by generator AND
  test — never hardcoded), and `SAFETY=10` is the `ρ/(1−ρ)`
  amplification headroom (principled, not a fudge).
- **Direct** result (single sweep, no outer iteration) →
  `assert_array_almost_equal_nulp(nulp=reduction_depth)`.
- **Bit-identity** for a pure-refactor is enforced by `-W
  error::DriftWarning` LAYERED on top — "the gate passes" ≠
  "bit-identical"; verify WHICH invocation ran.

Two traps that recur:

- **An ITERATED end-to-end snapshot CANNOT be the bit-identity gate** for
  a zero-numerical-change refactor — committed iterated DD snapshots
  already drift 1000s–100000s ULP under `-W error::DriftWarning`
  (cross-run/cross-machine FP jitter). The bit-id PROOF must descend to a
  single-step DIRECT snapshot on a fixed-seed RANDOM het ≥2G ψ with
  non-zero inflow (flat ψ nulls redistribution), captured pre-carve via
  the root-conftest `--capture-baseline` flag.
- **A non-fissile mixture has NO eigenvalue.** A `solve_sn` (eigen)
  snapshot on a moderator mixture is `k=0/abs→nan` — a silent dead test.
  To exercise a scattering/streaming operator path on a non-fissile
  mixture, reformulate as FIXED-SOURCE (`solve_sn_fixed_source` + uniform
  Q); the same path runs and is well-posed, corroborated by the
  closed-form `φ=(diagΣ_t−Σ_s0ᵀ)⁻¹Q` flat infinite-medium solution.

How to apply: the gate IS the claim — map quantity→config-key explicitly,
make the helper `-O`-safe, and demonstrate both DriftWarning modes
(informational + escalated) with a deliberate ±ULP perturbation. The
principled-equivalence ladder (bit-id ↔ nulp ↔ SAFETY×conv_tol ↔
sha256-golden) is chosen per the `vv` 3-criteria; a re-association that
names per-axis intermediates is principled-equivalent, gated at
`nulp(reduction_depth)`, anchored at k_inf + MMS — never at old-vs-new
ULP alone. Distinguish exact-cancellation from IEEE-754: a "clamp silent
on flat" / "w=½ byte-identical" premise from a τ-bearing or
`a+b`-summed reduction is exact-arithmetic-only — verify the IEEE
micro-fact at the prompt, not the docstring.

---

## L8 — Reuse the registry schema for cross-method infra; tag pairwise agreement L1, gate it max(tol)

When building cross-solver regression infra over an existing case
registry, a thin wrapper (pointer back to the registry case + per-solver
tolerances + pillar/claim-layer/truth-source tags) beats a richer
parallel schema (which loads the migration cost and becomes a second
drifting source of truth). Cross-method agreement tolerance MUST be
`max(tol_a, tol_b)` — tighter is reference contamination (`vv §6`). Tag
both truth gates and pairwise-agreement gates **L1** (the codebase
convention — each method's individual L1 backing makes the agreement
L1-strength when structural independence is genuine), NOT L4. Truth
values MUST trace to PRIMARY citations — transcribing from memory made up
two values the tests caught; verify against the literature memo first.

How to apply: reuse, don't reinvent; declare one-sided coverage honestly
(don't pretend coverage that doesn't exist); the harness knows
L0–L3+foundation only.

---

## L9 — V&V tagging idioms (the registry contract)

The audit is driven by `tests/_harness/registry.py` parsing markers —
mis-tagged tests show as orphans in every `session_briefing`:

- Module `pytestmark=[pytest.mark.verifies("label")]` for a default
  Sphinx label; override per-test for an extra equation. (No module
  `pytestmark` when it would clobber a foundation tag — tag explicitly.)
- `foundation` = a SOFTWARE invariant (data-structure / factory /
  algebraic-reduction); it MUST NOT carry `verifies(<physics-eq>)` —
  stacking them is silent level conflation (the audit credits the
  physics eq with a flat-flux foundation test). `math-origin` promotes
  carry no `verifies` (no `:label:` to point at).
- `@pytest.mark.slow` for >5s; the standard fast gate is `-m "l1 and not
  slow"`. Mark just the slow PARAMS (`pytest.param(..., marks=slow)`),
  not the whole function, when only some params are slow.

How to apply: see qa L-007 (foundation+verifies conflation) for the
removal recipe. `vv-status` rationale comments use (parens) not
[brackets] (docutils reads [xxx] as a citation → duplicate-citation
warnings under `-W`).

---

## L10 — The proactive-test-architect payoff: refute the plan's premises at design time

Dispatched proactively before an operator-algebra carve, the highest-value
output is REFUTING the plan's optimistic premises before the ink dries —
the plan repeatedly claims "bit-identical" / "clean O(h²)" / "improves on
flat" / "this arm is dead" / "N errors clear" and these are measurably
false:

- "iso-MMS bit-identical" was exact-arithmetic-only (≤1 ULP in IEEE).
- "aniso clean O(h²) at S16" was a #229-interpolation floor, untouched by
  the carve (needs S32).
- "improves-on-flat at the boundary" was unachievable (the slope is
  sub-floor; the correctly-consumed slope made the value slightly WORSE).
- **"the bare-`np.ndarray` arm is DEAD" — REFUTED by a runtime trace** (a
  signature-only carve). A production call site whose ARGUMENT is annotated
  `np.ndarray` is NOT evidence the typed arm superseded it — the annotation
  is the strongest reason to suspect the bare arm is LIVE. Trace the actual
  runtime object, do not read the type hint: `SNSolver.initial_flux_
  distribution()->np.ndarray` flowed bare into `F.apply(...)` (dispatched to
  the ndarray arm, returning `ndarray` not `ScalarSourceSink`), and `solve_sn`
  fed a bare `Q` to `add_iso/n2n_source` (in-place mutate, `None` return) —
  both load-bearing, neither dead. The retirement-audit's "stale annotation
  lies forward" cousin: a `T`-typed arm can be live precisely because the
  caller's `T` annotation is honest. The `assert_type(F.apply(arr),
  np.ndarray)` row is the static gate that REDS (→ `reportCallIssue`) if you
  wrongly retire it.
- **"N CLI-pyright errors clear" — never trust the count; assert it.** The
  N errors rarely map 1:1 onto the carve parts; some belong to an adjacent
  under-typing the carve doesn't touch (a `FullField`-arm `psi.bulk.
  integrate_angular()` gap is NOT cleared by retiring the *ndarray* arm).
  Gate with a verbatim post-carve `npx pyright <module>` and assert the
  EXACT residual.
- **"the same fold applies UNIFORMLY across the N solvers" — REFUTED by
  reading all N bodies.** A carve that routes a shared quantity (keff =
  prod/absn) "through one functional across 5 solvers" assumes the 5 sites
  compute the same thing; they rarely do. Reading every `compute_keff`:
  SN denom = `∫Σa·φ` (Σa already folds Σ₂.sum); **CP denom = net-removal
  `∫Σt·φ − ∫Σs·φ − 2·n2n`, NOT absorption** (equal value, different FP path
  AND n2n on the DENOMINATOR not numerator); diffusion denom = `∫Σa·φ +
  leakage`; homogeneous is **0-D (no `Σ_cells` — the volume fold is a
  category error, leave it)**. The "uniform `IRR(νΣf)/IRR(Σa)` fold" is
  clean only at the NUMERATOR; denominators are a different-physics family
  per solver. Gate each family separately; flag the CP net-removal rewire
  and diffusion leakage as per-solver MATH decisions (re-baseline per
  vv §bit-identity), not a uniform mechanical fold.

How to apply: measure the premise (temp-edit + revert, matched-quadrature
ladders, probe the regime the term dominates, RUNTIME-trace the consumer,
verbatim-count the pyright delta) rather than encoding the plan's hope into
a gate. A gate built on a false premise either calcifies a limitation,
falsely reds a correct term, or BREAKS a live path a "signature-only" carve
swore it wouldn't. For a type-contract deliverable the gate is a pyright-
checked `assert_type` block (extend the existing one — one home per
apply-contract); prove its teeth by mutating an overload return and
confirming `reportAssertTypeFailure` fires under CLI-pyright. State the
refutation in the plan so the implementer ships the achievable carve, not
the hoped-for one.

---

## L11 — An "axis-transpose of an existing reduction" reuses the machinery but ONE collapse rule flips — gate the DIFFERING rule, not the shared one

When a new operation is the axis-transpose of one that already shipped
(energy condensation = the energy-axis transpose of spatial homogenization;
both reuse the SAME `PetrovGalerkinFrame.project` = `G⁻¹M`), the implementer
will (correctly) copy the template — and the ONE place the two genuinely
DIFFER is exactly where the copy goes wrong. The gate's whole job is to
isolate the differing rule with a structurally-independent oracle; the
shared machinery is already verified by the original's gate.

- **The differing rule for condensation (vs homogenization): the MATRIX
  2-axis collapse.** Homogenize flux-weights BOTH axes of `Σ_s[g_from,g_to]`
  (`material_xs_field.py:403`, `project(_gather_legendre(l))`). Condense
  SUMS the sink (`Σ @ T`) then flux-AVERAGES the source (`project`). The
  vector channels, the `G⁻¹M` mechanics, the Gram pseudo-inverse, the chi
  birth-group-sum are all IDENTICAL or trivially analogous — only the matrix
  sink-vs-source asymmetry is new. So the load-bearing gate is the
  scattering 2-axis probe with three value-RED mutations: (i) swap
  sink/source, (ii) sum BOTH (drop source flux-weight), (iii) **project
  BOTH** — and (iii) IS the homogenize behavior, so it directly catches
  "implementer copied homogenize verbatim". Verified: only the
  sink-summed/source-averaged collapse preserves the hand-summed in-scatter
  rate `Σ_{g∈G}Σ_{g'∈G'} φ_g Σ_s[g,g']`; project-both BREAKS it (the rate
  gate reds on the copy-bug). Matrix collapse is ONE-ULP (the `@T` matmul
  reduction tree ≠ explicit sum) → `assert_allclose(atol=1e-12)`, never `==`
  (vv §bit-identity); vector channels ARE 0-ULP.
- **chi is the OTHER differing rule.** Homogenize production-weights χ ACROSS
  CELLS (a real average via `emission_frame`); condense is single-material,
  so χ collapses by PURE birth-group SUM `χ @ T` — NOT frame-projected (a
  projection would average → `Σχ ≈ 0.38 ≠ 1`, destroying the simplex). Gate
  the sum + `Σχ` preservation SEPARATELY from the projected channels.
- **The orientation trap is silent and ORDER-not-partition.** The canonical
  `eg` is DESCENDING; `IndicatorBasis.evaluate` buckets with
  `searchsorted(side="right")` = ASCENDING-only. Feeding descending coarse
  edges does NOT break the partition (row-sums stay 1) — it REVERSES the
  coarse-group ORDER (fast fine groups → thermal coarse column). So the
  intrinsic-property test must pin the coarse-group ORDER (fastest fine →
  coarse 0), not just partition completeness; demonstrate the reversal
  directly against the live `IndicatorBasis` (needs no SUT). Mode-6
  convention-drift with a green-partition disguise.

How to apply: for any "transpose/mirror of a shipped reduction" gate, (1)
diff the two production bodies and enumerate every rule that FLIPS; (2) make
the load-bearing gate a mutation-probe whose mutation IS the un-flipped
(copied-from-original) rule, and confirm it REDS against a hand-summed
rate/invariant oracle; (3) keep the shared-machinery checks light (the
original's gate covers them); (4) re-audit the convenient builder for a
silently-nulled channel (L1: `make_mixture` zeros SigL AND Sig2) before
trusting any "balance survives" leg. The skeleton goes GREEN-on-contact
(`@_needs` skipif + oracle math in-test); the convention traps + oracles run
UNCONDITIONALLY so they bite at design time.

---

## L12 — Folding a "fast path" INTO a composed form = principled-equiv, NOT 0-ULP; the UNCHANGED sibling stays 0-ULP; "transpose falls out free" is gated on the missing leaf(es) + a Mode-11 wrap. SPLIT variant (perf-forced) adds: no-tensor sentinel + OperatorSum heterogeneous-space guard + LD-monkeypatch re-pin

When a clean-before-extend carve re-expresses a bespoke "fast path" as an
instance of an already-composed operator — justified by a NORMALIZATION
IDENTITY (`Y₀⁰=1` makes `integrate_angular`=ℓ=0 analysis `M₀`, the iso `/W`
broadcast=ℓ=0 reconstruction `R₀`, so the iso P0 in-scatter =
`(1/W)·frame.conjugate(Λ)` restricted to ℓ=0) — three test-design facts
that recur and are easy to get wrong:

- **The fold is principled-EQUIV, never 0-ULP.** Same math (the identity
  guarantees it), but a DIFFERENT reduction tree (per-material matmul
  `einsum("fg,fc...->gc...")` vs the chained `R₀(Λ₀(M₀·ψ))` + the
  addition-theorem `(2ℓ+1)` broadcast). Gate the equivalence at
  `assert_allclose(rtol≈1e-14)` (reduction-depth ULP, `vv §bit-identity`
  crit-3), NEVER `array_equal`. Pair with a CLOSED-FORM iso anchor
  (`φ=(diagΣ_t−Σ_{s0}ᵀ)⁻¹Q` flat infinite-medium) — new-vs-legacy ULP is
  necessary-NOT-sufficient (L2/crit-2; can't tell you the pre-carve value
  was right).
- **The UNCHANGED sibling kernel MUST STAY 0-ULP.** If the fold adds a
  SEPARATE full-ℓ conjugation (skip_l0=False) for the iso path and leaves
  the ℓ≥1 `kernel` (skip_l0=True) untouched, the existing 0-ULP crosscheck
  on that kernel STAYS `array_equal` — do NOT relax it. Its RED is the
  CORRECT signal that the implementer accidentally re-routed the aniso
  path through the new fold (a scope violation). Audit disposition
  EXPLICITLY: which existing pin stays 0-ULP, which re-baselines to
  principled-equiv, which is net-new (the §1.5 re-baseline table).
- **"S† falls out free" rests on the ONE missing leaf, not the
  propagation.** `OperatorProduct.apply_transpose` propagates iff BOTH
  operands have `CAP_APPLY_TRANSPOSE`; the frame faces M/R already do, so
  the ENTIRE free-transpose claim reduces to adding `apply_transpose` (the
  per-ℓ group-flip `Σ_{s,ℓ}ᵀ`, index `[g,g']` vs forward `[g',g]`) +
  `CAP_APPLY_TRANSPOSE` to the leaf (`LegendreMomentScattering`, which was
  `{CAP_APPLY}`-only). Load-bearing gate = the LEAF transpose vs a hand
  loop (≥2G+4G asymmetric `Σ_{s,ℓ}`, the group-flip mutation reds). A
  forgotten `CAP_APPLY_TRANSPOSE` on the frozenset SILENTLY drops the
  capability through the product (no transpose, raises downstream) — gate
  the capability set directly (and FLIP the old "deferred" negative pin
  `CAP_APPLY_TRANSPOSE not in caps` → positive). The capability-set/impl
  mismatch is a distinct catcher from the value gate.
- **Euclidean transpose ≠ the metric Hilbert adjoint `.H`.** The carve
  wires `apply_transpose = (1/W)·kernel.apply_transpose` (plain matvec Sᵀ);
  the reciprocity gate is Euclidean `⟨Sφ,ψ*⟩=⟨φ,Sᵀψ*⟩` (`.sum()`, per-group,
  L27). Do NOT use `.H` (= `G_V⁺⊙apply_transpose(G_W⊙·)`, carries the SH
  Gram + `1/W` — a DIFFERENT identity that false-greens/reds). The dense
  `Sᵀ` oracle (build the matrix from FORWARD applies, transpose, apply)
  pins the VALUE structurally-independently but does NOT call
  `apply_transpose` → it is Mode-11-blind to the path; PAIR it with an
  in-process WRAP of the leaf transpose (counter>0) that pins the PATH.
- **A channel OUTSIDE the fold needs its own transpose gate.** When (n,2n)
  stays a parallel ℓ=0 term (separate `Σ_2n`, factor 2 — NOT folded into
  the frame), the highest-value mutation is "omit n2n from `S.apply_transpose`
  entirely" (transpose only the frame conjugation). Gate it with a `Sig2≠0`
  fixture (lessons L1: `make_mixture` AND every A/B/C/D library mixture zero
  `Sig2` → build it on the Mixture directly; #269 non-vacuity) + reciprocity
  that includes the n2n term on both sides.

- **The SPLIT variant (#276 A2-Option2): when the full-frame fold REGRESSES
  PERF (materializes the ℓ=0 moment tensor), the production keeps the iso
  SCALAR fast-path and SPLITS `S = (1/W)·[S_iso + S_aniso]` as an
  OperatorSum of two composed operators (S_iso = `R₀∘K_iso∘M₀` scalar faces,
  NO tensor; S_aniso = the UNCHANGED `kernel`). The transpose still falls out
  free — but now via `(A+B)ᵀ=Aᵀ+Bᵀ` (`OperatorSum.apply_transpose` propagates
  iff BOTH summands carry `CAP_APPLY_TRANSPOSE`), so the "free transpose"
  reduces to S_iso's leaves (M₀ᵀ/K_isoᵀ/R₀ᵀ), not just Λᵀ. THREE Option-2-only
  gates the full-frame fold did NOT need: (1) **the no-moment-tensor SENTINEL**
  — an in-process wrap on the tensor verbs asserting **0 calls** at L=0 is the
  LOAD-BEARING PERF invariant (it IS the reason the full-frame reverted); route
  P0 through the frame → counter>0 → RED. (2) **the OperatorSum heterogeneous-
  space guard** — `OperatorSum.__init__` raises iff `a.domain!=b.domain` (by
  `FunctionSpace (name,shape)`), SKIPPED when a side is `None`. S_aniso's outer
  space is the per-ordinate `L2[S^2]/(N,)` (conjugate wraps M/R); S_iso's outer
  must EQUAL it (reuse `frame.measure.space`, Pattern-2) OR be `None`. A
  freshly-MINTED per-ordinate space name (same `(N,)` shape, diff name) RAISES;
  a `None`-skip makes the sum's shape-guard VACUOUS (then the LD-multi-moment
  equivalence gate is the SOLE shape catch). Gate it with a positive-construct
  + negative-mismatch-RAISES pair (anti-#11). (3) **the LD-monkeypatch RE-PIN**
  — a pre-existing slope-source mutation gate that monkeypatches the OLD shared
  combine (`_assemble_per_ordinate_source`) goes VACUOUS when the split retires/
  reshapes it (Mode-11); re-pin to the NEW production line (`K_iso.apply` /
  `apply_p0_in_scatter`, which carries the `...` spatial-moment spectator) +
  sentinel-wrap that the gate EXECUTES it. The LD config is mandatory for the
  oracle gate (the full-frame regression was LD-only). The oracle
  (`full_scatter_kernel` = the moment-tensor form, KEPT per fuller-view-oracle)
  pins forward AND transpose VALUE structurally-independently.

How to apply: for an iso/fast-path→composed fold, the §0.5 named risk is
THE equivalence (principled-equiv 1e-14 + closed-form anchor + UNCHANGED-
sibling-stays-0-ULP). For the "transpose falls out free" extend, gate the
ONE leaf transpose (hand loop, group-flip mutation) + capability-set FLIP +
Euclidean reciprocity (asymmetric SigS, per-group, +discriminator) + the
Mode-11 leaf-WRAP (dense oracle pins value, wrap pins path) + the outside-
the-fold channel's own transpose. For the SPLIT variant add the no-tensor
sentinel + the OperatorSum space-guard construct/mismatch pair + the
LD-monkeypatch re-pin (above). Re-audit the Mode-7 trap: a P0-only
fixture NEVER builds the ℓ≥1 blocks and is BLIND to the unified ℓ≥0 fold (or
the iso/aniso SPLIT) — run every equivalence/transpose gate on an
ANISOTROPIC+het+Sig2≠0 fixture, reserve P0-only for the perf P0-path gate.
`vv §bit-identity / Mode-7 / Mode-11`; mirrors qa/method-implementer
fission-adjoint (A1) template.

- **The CROSS-MODEL energy-operator EXTRACTION variant (#276 #119, the 3rd
  design — SUPERSEDES the SPLIT):** the iso/aniso SPLIT is dropped; the iso
  ℓ=0 projection becomes a standalone `K_iso = IsotropicScattering +
  IsotropicN2N` OperatorSum on the SCALAR flux, routing through the EXISTING
  `mat_xs` scalar verbs (`apply_p0_in_scatter`/`apply_n2n`,
  material_xs_field.py:613/652) + NEW scalar transpose twins
  (`apply_p0_in_scatter_transpose`/`apply_n2n_transpose`, mirror the in-tree
  moment twins at :751/:840), extracted to `transport/operators/` (sn-free),
  consumed at THREE sites incl. a CROSS-MODEL homogeneous LHS-fold. This is
  SIMPLER to gate than the SPLIT: same `mat_xs` verbs + same `+=` order ⇒ the
  forward `K_iso≡fast-path` AND the whole-`_assemble_per_ordinate_source`
  bit-id are **0-ULP `array_equal`** (NO reduction-tree change — the SPLIT's
  frame faces are gone, so it is NOT principled-equiv). The SPLIT's two
  headline risks DISSOLVE: OperatorSum-space is LOW (both summands share
  `mat_xs`+scalar space) and the no-tensor sentinel is moot (scalar path never
  touches the tensor). The risk MOVES to ONE new place: **the cross-model
  `as_dense` consumer holds a bare `Mixture`, NOT a `MaterialXSField`** —
  `homogeneous/solver.py:93` builds `A=diags(SigT)−SigS0_T−2·Sig2_T` from
  `mix.SigS[0].T` with NO mesh / NO `cells_by_material`, so `K_iso.as_dense()`
  REQUIRES the operator be constructible from bare ENERGY data (a `from_mixture`
  ctor → `{0: mix.SigS[0].toarray()}`); else homogeneous must import an SN mesh
  and the "sn-free/model-independent" claim is FALSE (DEFER P4 to an
  AskUserQuestion if not). VERIFY-IN-TREE payoff (L10): the brief's
  "diffusion/homogeneous build A=diags−K_iso" is FALSE for diffusion
  (`diffusion/solver.py:240` is matrix-FREE + 2G-hardcoded
  `(sig_a+sig_s)*fi−sig_s[::-1]*fi[::-1]`, builds NO dense A) ⇒ `as_dense` has
  ONE real consumer, shape = per-material `dict[int,(ng,ng)]` (mirror
  `n2n_matrix`), NOT a block-diagonal `(N·ng)²` (YAGNI). Add a VACUUM bit-id
  correctness gate (Sig2=0 ⇒ `as_dense()≡` legacy `SigS0_T` exactly — differ ⇒
  STOP it's a bug; `snapshot_migration_when_production_goes_bare.md`) + the
  homogeneous keff invariance vs the closed-form analytical k_inf (2eg+4eg,
  `<1e-12`, NOT old-vs-new). Production S† rides the `full_scatter_kernel`
  ORACLE (frame form, not K_iso) ⇒ the S†-oracle-equiv leg is a near-TAUTOLOGY
  (prod IS that expression); the LOAD-BEARING S† pin is Euclidean reciprocity
  `⟨S.apply,χ⟩=⟨ψ,S.apply_transpose⟩` (fwd fast-path vs adjoint frame-form =
  two structurally-different reps). The LD-monkeypatch on
  `_assemble_per_ordinate_source` SURVIVES (it wraps the OUTER combiner the
  carve keeps) — but ADD a Mode-11 wrap on `IsotropicScattering.apply`
  (counter>0) AND a direct-on-`IsotropicScattering` mutation (the old
  monkeypatch is now one layer ABOVE the changed reader, so a stamp bug it
  faithfully passes through could be masked). Plan
  `.claude/plans/a2_kiso_verification.md`; SUPERSEDES `a2_option2_verification.md`.

- **PHASE C continuation (#276 — homogeneous USES `C=MultiplicationOperator`
  + `F=FissionOperator` as operators via APPLY-TO-BASIS, dropping
  `np.diag`/`np.outer`):** A=`C−K_iso` (OperatorSum), F=fission dyad; the dense
  `(ng,ng)` realized by applying each op to the ng single-cell `(ng,*sp)` columns
  `e_i`. The GAP: the 206-caller FullField-only `MultiplicationOperator.apply`
  gains a `singledispatchmethod` bare-`ndarray` arm (`coeff.values*x` = the
  per-ordinate engine action; FullField arm UNCHANGED = the conversion's
  bit-id-inheritance guard ⇒ existing `test_multiplication_operator.py` +
  SN forward/adjoint suite stay 0-ULP green). **The LOAD-BEARING gate is
  CROSS-ARM CONSISTENCY, NOT bare-vs-its-own-`coeff*x` (tautological):**
  `apply(φ_bare (ng,*sp))` ≡ `apply(φ BROADCAST to (N,ng,*sp) FullField)
  .bulk.values[n]` ∀ ordinate n — the FullField-engine path
  (`DiagonalOperator(coeff,broadcast_axes=(0,))` = `coeff[None]*x_full`) is the
  structurally-independent reference, 0-ULP `assert_array_equal` (same two floats
  multiplied). Build the broadcast input: `np.broadcast_to(φ[None],(N,ng,*sp))
  .copy()` into `TimedFullField.zeros(...).bulk.values`. ADD an EXPLICIT
  `_require(bare.shape==(ng,*sp))` — `assert_array_equal` BROADCASTS so it ALONE
  misses the extra-leading-axis mutation (`coeff[None]*x`→`(1,ng,*sp)`). The
  single-cell meshless homogeneous consumer is AXIS-BLIND (spatial=`(1,)` nulls a
  group↔spatial transpose) ⇒ the intrinsic gate MUST use the existing
  `_cartesian_2d_mesh(nx=5,ny=3,ng=2)` + asymmetric σ (config-blindness L1);
  homogeneous k_inf catches only the COMPOSITION (sign-flip `C+K_iso` /
  C-omission `−K_iso` / K_iso-omission `=C` all move k O(1) on the all-channel
  `homo_2eg_n2n` — the existing `test_kinf_exact` covers it, NO new homogeneous
  CORRECTNESS test needed). NEW homogeneous gate worth adding = a SHARP A-level
  procedural pin: expose `_assemble_loss_matrix`→ `assert_allclose(A_applytobasis,
  diag(σt)−(σs0+2σ2)ᵀ, atol=1e-12)` (localizes a sign/term bug faster than the
  end-to-end eig; PROCEDURAL — shares `mat_xs` data ⇒ PAIR with the structural
  `case.k_inf` anchor, L2 two-anchor). Mutations: (i) axis-transpose + (ii)
  group-drop → intrinsic only (homog axis-blind); (iii) extra-leading-axis →
  intrinsic shape-assert (downstream crash/mis-assemble is unreliable).
  `apply_transpose` delegates to `apply` (self-adjoint real coeff ⇒ free bare
  arm; gate pins the delegation survives); `solve` STAYS FullField-only (no
  consumer). Dispatch gate: an unregistered type (ScalarFlux / `object`) → base
  `TypeError` (mirror F). Mode-11 wrap on the ndarray arm (counter>0 during
  `solve_homogeneous_infinite`) confirms the k_inf gate EXECUTES the new reader.
  **Scope = ndarray-ONLY** (F has +ScalarFlux ∵ a live SN scalar-fission
  consumer; C has NO scalar-flux collision consumer yet — Pattern-6 defer;
  trigger = the first diffusion/CP/MoC refounding on C). K_iso's input surface
  is duck-typed `_values_of` (model-portable primary), C's is FullField-primary
  + bare-secondary — the three operators' surfaces differ ∵ their consumer sets
  differ (honest, not inconsistent). `vv §bit-identity / Mode-11`; mirrors the
  K_iso A2 plan above.

---

## L13 — Verifying a "capability-string → typed-operator" carve (`.solve`/`CAP_SOLVE` → `.inverse()`; `CAP_APPLY_TRANSPOSE` → `.H`): the equivalence keystone is bit-id INHERITANCE; the WHOLE frozenset (BOTH axes) retires SAFELY because the "iff both" law lives in the METHOD BODY — licensed by a FAITHFULNESS keystone (derived-predicate ≡ old-frozenset ∀ operator, before deletion); a half-modernized axis is itself the twin-path HAZARD

When a carve replaces a stringly-typed capability + gated method
(`CAP_SOLVE` + `.solve`) with a typed operator (`.inverse()` returns a
`SweepInverseOperator`/`KrylovInverseOperator`), three reusable
test-design facts (#226 B4 Phases 2–5, plan
`issue_226_inverse_operator_verification.md`):

- **The keystone equivalence gate is bit-identity INHERITANCE, not a
  fresh value claim.** `A.inverse().apply(b) == A.solve(b)`
  `assert_array_equal` is bit-exact BY CONSTRUCTION (the inverse
  operator WRAPS the legacy method) — it proves the wrapping changed
  NOTHING, necessary-NOT-sufficient (L2). It RIDES the EXISTING
  closed-form anchors (`test_removal_form_kinf_independent_reference_2g`,
  keff k_inf) for "was right"; do NOT mint a new reference. The
  curvilinear seed (`initial_guess`) MUST be threaded on a
  sphere/cyl seeded variant (slab is BLIND — §0.6); PAIR the value gate
  with a Mode-11 in-process wrap capturing `initial_guess is seed` (the
  converged solve may be seed-INSENSITIVE, nulling the value gate — only
  the wrap pins the threading unconditionally, L4/L12).
  - **The DELIVER-ALL nine advertisers (`is_invertible` → a working
    `.inverse()`): do NOT trust a brief's assumed inverse REALIZATION —
    read the code.** The elegant choice is DELEGATION-to-`solve` (a
    generic `InverseOperator` whose `apply` IS `leaf.solve`, the division
    `b/c`), NOT a reciprocal-coefficient twin — because (a) `(1/c)·b ≠
    b/c` drifts and (b) `1/Σ` is a units-dishonest mean-free-path
    (coding-elegance Pattern 3). This SUPERSEDES the "reciprocal Diagonal
    ~1-ULP" premise and splits into TWO DISTINCT tolerance claims:
    (i) `inverse().apply == solve` is **BIT-IDENTICAL** for ALL nine
    (delegation + bit-id-preserving structural composites — `Perm→take(π⁻¹)`,
    `(AB)⁻¹=B⁻¹A⁻¹`, `(αL)⁻¹=(1/α)L⁻¹`, `Id→x`); (ii) the round-trip
    `inv.apply(A.apply(x))≈x` is still per-type **nulp** (the two
    directions use ×-then-÷, DIFFERENT arithmetic — nulp even when
    equivalence-to-solve is bit-id). Tolerance = the round-trip's FP
    reduction depth (count ops: `array_equal` for gather/passthrough
    Perm/Id; `nulp(2)` for Diagonal ×÷; sweep `rtol=1e-10`, **sphere
    excluded** #282) — NEVER a blanket loose tol. And `(A⁻¹)⁻¹ = A`
    holds by **OBJECT IDENTITY** for the `InverseOperator`-backed leaves
    (`inverse()` returns `self.inner`) — assert `is A`, NOT the brief's
    "action-equality at ~ULP" (that degrades to action-nulp ONLY for the
    recomputed-scalar `ScaledOperator`, `1/(1/α)`; NOT testable on
    `InvertibleOperator` — `SweepOperator.inverse` is #280-deferred).
    The universal round-trip (I1) is the STRONGEST Step-1 gate: it needs
    NO legacy `.solve` (survives the frozenset/solve retirement) and pins
    the inverse against its OWN forward + a closed-form value anchor.
    Every mutation names its ACTIVATING config (§0.6 applies to
    MUTATIONS): non-commuting `Diag@Perm` for the `(AB)⁻¹` order-swap,
    NON-involution 3-cycle for the perm-inverse, α≠1 for the scale,
    non-unit `c` for the delegation. Negative controls RAISE eagerly at
    `.inverse()` (not `assert` — `-O`); α=0 is a ctor `ValueError`, a
    DIFFERENT edge from the value-dependent singular (`min|f|=0`).

- **A LAW-gating capability CAN be retired — the law is in the METHOD
  BODY, not the frozenset — and leaving it as a string while a sibling
  axis modernizes is ITSELF the hazard (twin-path, Cardinal Rule 2).**
  My first pass recommended KEEP `CAP_APPLY_TRANSPOSE` (retire `CAP_SOLVE`
  only) because it gates `OperatorSum/Product.apply_transpose` "propagates
  iff BOTH" (the #276 "S† falls out free" law). The user OVERRULED it: a
  modernized solve axis (`.inverse()`/`SupportsInverse`) beside a
  stringly-typed transpose axis (`CAP_APPLY_TRANSPOSE`) is TWO mechanisms
  for one structural question — forbidden. The resolution: retire the
  WHOLE `capabilities` frozenset (`CAP_APPLY`+`CAP_SOLVE`+
  `CAP_APPLY_TRANSPOSE`) via ONE mechanism per axis (`.inverse()`/`.H`
  operator-returning methods + `SupportsInverse`/`SupportsAdjoint` STATIC
  Protocols + recursive derived `@property is_invertible`/`is_adjointable`
  + eager `NotInvertible`/`MissingAdjoint`). SAFE because the "iff both"
  LAW already lives in the method body (`OperatorSum.apply_transpose =
  a.apply_transpose(x)+b.apply_transpose(x)`; `OperatorProduct.solve =
  b.solve(a.solve(·))`) — the frozenset is only a redundant
  ADVERTISEMENT. The retirement does NOT touch S† (the reciprocity LAW
  `a.H+b.H` is untouched); only the advertisement moves.
  **The FAITHFULNESS keystone is what licenses the deletion:** a SCAFFOLD
  gate asserting `is_invertible == (CAP_SOLVE in caps)` AND `is_adjointable
  == (CAP_APPLY_TRANSPOSE in caps)` for EVERY operator (leaves + every
  composite, the value-dep zero-coeff + half-adjointable asymmetry
  fixtures MANDATORY — else blind to the axis separation), proved GREEN
  during coexistence, THEN deleted WITH the frozenset. Its PERMANENT
  successors are the recursive composition pins (`(a+b).is_adjointable ==
  a.is_adjointable and b.is_adjointable`, etc. — they never reference the
  frozenset) + the migrated closure-law tests + the S† method canaries.
  Mutations M-ADJ-PROP (break the recursive `and`) and M-ADJ-FORGE (force
  `is_adjointable=True` on a half-adjointable sum) both RED the keystone;
  M-ADJ-FORGE downstream is a RAISE not a wrong-value (a half-adjointable
  summand has no transpose to forge — the forge bypasses the eager `.H`
  raise → AttributeError).
- **Runtime-query form = derived `@property`, NEVER `isinstance(Supports*)`.**
  `runtime_checkable` Protocol `isinstance` is CLASS-uniform (checks method
  PRESENCE) → blind on composites (every `OperatorSum` HAS
  `apply_transpose` on the class, so a half-adjointable sum passes
  `isinstance(SupportsAdjoint)`) AND blind to the value-dep singular edge
  (`isinstance(M_zero, SupportsInverse)` True even at `min|f|==0`). Only
  the recursive/data-reading `@property` is instance-accurate. The
  Protocols earn their place as the STATIC contract ONLY (pyright +
  annotation target). The genuine pre-query consumer forcing the property
  is KrylovAccel preconditioner-SELECTION (solve axis); the ADJOINT axis
  has NO external pre-query consumer (read only by composers + the eager
  `.H` raise) → smaller migration than solve.
- **The LITERAL-STRING precondition trap.** Two landed S† canaries
  (`test_g_adjoint_reciprocity:210`, `test_removal_form_matvec_sweep:308`)
  guard with `if "apply_transpose" not in A.capabilities` — a BARE string,
  so a `CAP_APPLY_TRANSPOSE`-CONSTANT grep MISSES them, and they
  `AttributeError` the instant the frozenset deletes. "S† stays green" is
  NOT automatic — the deletion REQUIRES rewiring these preconditions to
  `not A.is_adjointable`. Before any capability retirement, grep BOTH the
  constant AND the literal string `"solve"`/`"apply_transpose"` membership
  reads; the consumer surface is WIDER than the carve plan's sketch
  (CAP_SOLVE wove through DiagonalOperator's value-dep spectrum law, the
  sn/solver resolvents, every composer, KrylovAccel's preconditioner
  SELECTION — not just SourceIteration; the adjoint cap is read by every
  composer + the boundary-op caps-derivation).

- **The adjoint-of-inverse `(A.H)⁻¹ = (A⁻¹).H` is verified by the
  RECIPROCITY identity, not a round-trip tautology.** For `x =
  A.H.inverse().apply(b)`, assert `⟨A.apply(ψ), x⟩_G = ⟨ψ, b⟩_G` for
  random ψ — the FORWARD matvec + G-metric is structurally independent
  of how `x` was computed (transpose-solve via `loss_action_transpose`),
  and reuses the EXACT metric `test_g_adjoint_reciprocity` validates.
  The round-trip `A.H(A.H.inverse(b)) ≈ b` is necessary-not-sufficient
  (both arms share the metric) — keep it but PAIR with reciprocity. The
  metric-vs-Euclidean discriminator (L12: `.H` ≠ Euclidean `Aᵀ`)
  ACTIVATES on sphere/cyl only (non-trivial trace metric) — the
  drop-the-G-metric mutation reds curvilinear, stays green on slab.
  `A.H.inverse()` rides the LANDED transpose matvec ⇒ this gate
  pre-builds #280's substrate; if the carve defers the transpose-solve,
  the gate is `xfail(strict=False, reason="#280 …")` that flips to
  xpass when #280 lands (NEVER strict=True — a stale failure satisfies
  it for the wrong reason, L4).

How to apply: for a capability→operator carve, (1) keystone =
array_equal inheritance + existing closed-form anchor + seed Mode-11
wrap; (2) classify EVERY retired capability string structural-vs-law-
gating and mutation-prove the law-gating one's retirement ≠ no-op;
(3) the value-dependent singular edge (DiagonalOperator zero-coeff,
`CAP_SOLVE iff min|f|>0`) survives as a runtime RAISE inside `.inverse()`
(positive+negative pair, anti-#11; RAISE not bare-assert — `-O`), NOT a
frozenset; the STATIC `ScaledOperator(0.0)` edge is already a ctor
`ValueError`. The `# type: ignore` deletions are pinned by EXTENDING the
`assert_type` C6 block (`test_operators_apply_typed.py`) + CLI-pyright
`reportAssertTypeFailure` teeth (L10). `vv §bit-identity / Mode-8/11 /
anti-#2,#11,#19`; mirrors L12 (the scattering-fold sibling) + L2/L4/L10.

---

## L14 — Verifying the FIRST *iterative* inverse of a family (`GreenOperator` = `OperatorSum.inverse()` wrapping the SI driver) + the `is_invertible=a.is_invertible` contract that FLIPS frozen pins + a wrap-delegate MIXIN extraction

When a carve adds an ITERATIVE inverse-as-operator (the inverse WRAPS a
convergent driver, not a legacy `.solve`) and changes a composite's
invertibility contract to route to it, six reusable test-design facts
(#226 taxonomy §12 step 4, spec `issue_226_inverse_operator_verification.md`
PART III):

- **The first iterative inverse has NO legacy `.solve` to inherit — so its
  keystone is STRUCTURAL-INDEPENDENCE, not bit-id inheritance (contrast
  L13).** Its correctness rides a closed-form dense-LU anchor
  (`np.linalg.solve` of the materialized sum) + the distinguishing
  name-earning invariant, NEVER old-vs-new ULP. And once the composite
  gains `solve(b) := self.inverse().apply(b)` (the conditional-CAP_SOLVE
  backing method), the `inverse().apply ≡ solve` gate is a **TAUTOLOGY**
  (solve IS defined as inverse().apply) → it proves NOTHING; do NOT add it
  as coverage (it was bit-id evidence in L13 only because the inverse
  wrapped a PRE-EXISTING solve). Round-trip tol is DRIVER-TOL
  (`SAFETY×tol`, L7 iterative), NOT nulp — the one row that breaks the
  §12.0 "count arithmetic ops → nulp" rule.
- **A raise-on-non-converge gate (`ConvergenceFailure`) MUST test the TRUE
  residual `‖A·ψ−q‖/‖q‖`, NOT the driver's iterate-INCREMENT.** Signature 9
  (ρ-blind stopping): the increment under-reports the true error by
  `1/(1−ρ)`, so an increment-check passes a non-converged iterate SILENTLY
  as `c→1` — defeating the exception's Cardinal-Rule-1 purpose. Falsifier =
  a `c=0.99` near-critical split where increment<tol but true residual≫tol.
  This is a FIRM design objection to raise if the brief says "final
  residual" ambiguously (the driver returns `(psi, history)` and does NOT
  raise on max_iter — so the wrapping op owns the check). The divergent
  edge PAIRS with a convergent CONTROL (a bare `pytest.raises` with no
  no-raise control is Mode-10-blind); the raise is a `raise` not a bare
  assert (`-O`).
- **`is_invertible = self.a.is_invertible` (left-spine leading term) FLIPS
  frozen `is_invertible is False` pins False→True BY DESIGN** — a plain sum
  with an invertible LEADING term now HAS a general inverse (the Green), so
  "no general (A+B)⁻¹" is retired. Retirement-audit: `grep -rn
  "is_invertible is False" tests/` WITHOUT an adjoint filter (a `grep -iv
  adjoint` HIDES the `is_adjointable…is_invertible is False` compound-line
  hit) and classify EACH: leaf / product (`a and b`) / scaled STAY False;
  plain-sum-with-invertible-leading FLIPS. The predicate is
  "leading-term-preconditionable AT THIS SPELLING", NOT "mathematically
  invertible" (`(−S)+A` is math-invertible yet reports False) — acceptable
  ONLY given the #261 canonical-ordering precedent + zero current consumers
  + a LOUD refusal; doc it + pin the refusal message.
- **The faithfulness keystone (`is_invertible ≡ CAP_SOLVE∈caps`, L13) stays
  green ONLY if the conditional `CAP_SOLVE` add is LOCKSTEP** — both read
  `self.a.is_invertible` (single source). The `a`-recursion needs BOTH
  asymmetric fixtures `inv+applyonly` AND `applyonly+inv`: a both-invertible
  sum is BLIND to a/b/and/or (§0.6); only the `(True,False)` pattern
  uniquely IDs the `a`-rule (b→(F,T), and→(F,F), or→(T,T)). Mutation
  M-SUM-CAPDRIFT (add CAP_SOLVE under `_has(a) and _has(b)` ≠ the property)
  reds faithfulness.
- **The order-dependent-strategy trap has THREE edges, each needing teeth**
  (#261): `L+C` → `SweepOperator` (MRO shadows the generic sum inverse —
  assert `type(...) is SweepOperator`, NOT Green); `C+L` → a Green whose
  collision-preconditioned Richardson DIVERGES → loud `ConvergenceFailure`
  at apply (not silent wrong answer); `(−S)+A` (leading non-invertible) →
  REFUSES at construction with the canonical-ordering message. The Green
  auto-split must FLATTEN the left spine of EXACT-`OperatorSum` nodes ONLY
  (subclasses like `InvertibleOperator` are leaf-terms) — mutation
  M-GRN-FLATTEN (flatten THROUGH the subclass) makes `A_loss.inverse()`
  refuse (leading `L` non-invertible) where it should build a Green.
  Green-vs-hand-`SourceIteration` equivalence (auto-split ≡ canonical SI) is
  bit-id-inheritance (necessary-not-sufficient) → PAIR with the dense-LU
  anchor.
- **A NEW iterative inverse's seed-threading needs its OWN Mode-11 pin — a
  landed per-iterate spy ONE LEVEL BELOW is BLIND to it.** `Green.apply(q,
  initial_guess=x0)` threads x0 as the driver's START; the landed spy
  (wraps the inner `InvertibleOperator.solve`, pins driver→inner per-iterate
  threading) stays GREEN if Green drops x0 (each iterate still seeds the
  next), AND the drop is value-invisible (warm-start changes RATE not the
  converged value) → wrap the INTERNAL `SourceIteration.solve`, assert the
  first call receives x0 (values-equal, not `is` — a defensive copy
  false-REDs `is`). Name-earning invariant triage (taxonomy §13 bar):
  G-Neumann (partial sums `Σ_k(P⁻¹N)^k P⁻¹q`, ratio≈c) + forward dense-LU
  anchor = LOAD-BEARING; G-reciprocity = INCLUDE-lighter (the ONLY thing
  pinning the transposed-operand split without a 2nd oracle — needs a
  NON-symmetric A; Euclidean transpose, NOT `.H` = #280-adjacent); G-kernel
  (`G·δ_j`=col j of A⁻¹) = FOLD into the anchor (δ_j inputs). G-Neumann-L1
  Mode-9 config = het ≥2G VACUUM slab, NOT isotropic reflective box (flat
  flux nulls the multiple-scattering the series expands).
- **A wrap-delegate MIXIN extraction (shared back-half of ≥2 inverse twins)
  is behaviour-preserving → the no-op proof is "existing suites green
  UNCHANGED"; the NEW gate is the pyright ABSTRACT-signature enforcement**
  (`apply(x,/,*,initial_guess=None)` on the mixin → `reportIncompatible
  MethodOverride` if a sibling drops the kwarg — the #285 structural-vs-
  convention decision resolves to STRUCTURAL) + a mixin mutation that reds
  ALL siblings (M-MIXIN-CODOM/INV — single-source proof). `SupportsSeeded
  Apply` conformance is a STATIC (annotation) gate, not runtime (not
  runtime_checkable).

How to apply: for an iterative-inverse-as-operator carve, (1) keystone =
dense-LU anchor + G-Neumann, NOT bit-id (no legacy solve) — and flag the
`inverse().apply≡solve` tautology; (2) demand the raise checks the TRUE
residual (Signature 9 falsifier); (3) grep-audit the `is_invertible is
False` flip WITHOUT an adjoint filter, classify leaf/composite, migrate the
premise ("no general (A+B)⁻¹" retires); (4) faithfulness lockstep +
both-asymmetric-fixtures; (5) three ordering edges (MRO-shadow / divergent-
loud / leading-non-invertible-refuse) + flatten-exactness; (6) the new op's
OWN seed Mode-11 pin; (7) mixin no-op = existing-suites-unchanged + pyright
abstract-sig teeth. `vv §bit-identity / Mode-8/9/11 / anti-#2,#11 /
Signature-9`; the iterative sibling of L13 (which handled the DIRECT/exact
inverses); mirrors L2/L4/L7/L10/L12.

---

## L15 — Verifying the RETIREMENT of a coexisting mechanism (the frozenset→typed-predicate TERMINAL step): keystone-v2 is a structural/value CONTRACT (not a comparison), migration is a mechanical-rule + completeness-regrep (not a table), land-order is ATOMIC, and the pruned-body re-route rides inheritance+independent-anchor+sentinel with a direct/iterative TOLERANCE SPLIT

When a carve RETIRES a coexisting advertisement (the `capabilities`
frozenset, BOTH axes) that a typed successor already replaced (the ADDITIVE
steps = L13/L14 — `is_invertible`/`is_adjointable` + `.inverse()`/`.H` +
`SupportsInverse`/`SupportsAdjoint`), four terminal-step disciplines (#226
§12 step 6, spec PART V §36–§44):

- **The coexistence-era faithfulness scaffold (`is_X ≡ CAP_X∈caps`) DELETES
  with the frozenset — its PERMANENT successor is a structural/value-split
  CONTRACT, not a comparison (keystone v2).** Per operator: `is_invertible
  True ⟹ .inverse() returns`; False ⟹ EITHER `_STRUCTURAL_ABSENT`
  (`not hasattr` — Zero/masks/RankOne/streaming-L/energy S·F leaves/B) OR
  `_VALUE_RAISE` (declares `.inverse()`, raises the eager `NotInvertible` —
  zero-coeff Diag/Mult, non-invertible-head sum/product/scaled). PIN WHICH
  PER TYPE — the split IS the contract. Mirror the adjoint axis
  (`is_adjointable ⟺ .H returns vs raises MissingAdjoint EAGERLY`; uniform,
  no structural-absent — `.H` is on the base). Add bridge-consistency
  (`invertible(op)==op.is_invertible` pins the one-line TypeGuard body vs
  drift — M-BRIDGE `return True` reds it). The recursive-composition pins
  survive VERBATIM (never read the frozenset); only the scaffold CALLS strip
  out. Rewrite the ONE shared helper (Cardinal Rule 2), `-O`-safe
  explicit-`raise` (un-collected module). The two asymmetry fixtures (Zero =
  adjointable-not-invertible; TRUE-zero-coeff Mult; half-adjointable
  `full+apply_only`) stay MANDATORY (blind to the axis SEPARATION otherwise).
- **Migration is a MECHANICAL RULE + a completeness RE-GREP gate, NEVER a
  fixed table.** The test surface is WIDER than any map (127 reflective
  `.capabilities` reads / 33 files; steps keep ADDING files a stale map
  predates — RE-GREP yourself, don't trust the plan's inventory). Rule:
  `CAP_SOLVE∈caps→is_invertible`, `CAP_APPLY_TRANSPOSE∈caps→is_adjointable`,
  strict `==frozenset({...})`→the axis-predicate CONJUNCTION,
  `pytest.raises(MissingCapability)`→by-AXIS (`NotInvertible`/`MissingAdjoint`/
  plain `TypeError match="apply"`). Bare-STRING `"apply_transpose" not in
  op.capabilities` PRECONDITIONS (§1d.2) are MISSED by a CAP-constant grep —
  grep `.capabilities` too; leaving ONE unrewired `AttributeError`s
  post-deletion (a false RED masking whether S† held — M-LITERAL-STRING).
  DONE only when `grep -rn ".capabilities\|MissingCapability\|CAP_*\|_has("
  tests/ orpheus/` returns ZERO.
- **Land-order is W1-additive / W2-ATOMIC / W3-teeth.** The frozenset
  deletion breaks ALL reflective reads at ONCE → W2 is atomic (delete +
  full test-migration in one commit). The exception REPLACE
  (`MissingCapability`→`NotInvertible`/`MissingAdjoint`) must ALSO be
  W2-atomic — a W1 flip breaks every `pytest.raises(MissingCapability)`
  mid-flight (sibling TypeError subclasses don't cross-catch). The
  eager-`.H` flip (raise site `.apply`→`.H`, a REAL behavior change) is the
  ONE W1 change, with its own `.H`-raise-site test migration (M-ADJ-EAGER;
  lazy `wrapper=A.H; raises: wrapper.apply` → `raises(MissingAdjoint): A.H`).
  Re-baseline the pyright ratchet DOWN in-commit (deleting N `# type:
  ignore` + `cast(SupportsInverse)` clears reds; assert the EXACT new floor,
  never trust the count — L10); guard-narrow deletions are CLI-pyright teeth
  (`reportAttributeAccessIssue` — M-GUARD-DELETE-PYRIGHT).
- **A "solve pruned to native realizations" body-rewrite (the composed
  `OperatorProduct.solve` re-route `b.inverse().apply(a.inverse().apply(q))`,
  which FIXES the deleted-solve factors) needs bit-id-INHERITANCE (pre-carve
  `--capture-baseline`) + a structurally-INDEPENDENT dense anchor
  (`np.linalg.solve(as_matrix,q)`) + a Mode-11 SENTINEL.** The re-route
  EXECUTES (`(A@B).inverse().apply` → `InverseOperator.apply=inner.solve` →
  the re-routed body) — but a value-ref spelled the SAME way is tautological,
  so wrap `OperatorProduct.solve` (counter>0). The factor-kind tolerance
  SPLITS: DIRECT factors (Diag/Perm/Scaled/Identity — integer gather /
  identical float path) are `array_equal`; the ITERATIVE sum-factor row
  (Green-vs-Green — math-bit-id but driver-iterated) is
  `assert_regression(SAFETY×driver_tol)` NOT `array_equal` (L7 iterated-drift
  trap). Config-blindness: NON-COMMUTING factors (Diag@3-cycle-Perm) — the
  commuting-D1@D2 designed-green control proves the order-swap mutation is
  BLIND without them.

How to apply: for a coexisting-mechanism retirement, the scaffold gate is
TEMPORARY (deletes with the mechanism) — design its PERMANENT
structural-contract successor FIRST; migrate by rule+regrep not a table
(re-grep the surface yourself — it is wider than any map); retire
atomically (frozenset AND its exception class together); a pruned-body
rewrite rides inheritance + independent-anchor + Mode-11-sentinel with a
direct/iterative tolerance split. The terminal sibling of L13 (the ADDITIVE
design this EXECUTES) / L14 (the iterative Green) / L7 (regression
tolerance) / L10 (assert the pyright delta). `vv §Mode-8/11 / bit-identity /
anti-#2,#11,#19`.
