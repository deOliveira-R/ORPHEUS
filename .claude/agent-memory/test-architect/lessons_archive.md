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
- **Retiring a runtime guard that had NO negative test → its
  replacement's teeth are NET-NEW, not migrated.** Before crediting a
  guard→invariant mechanism-swap as behavior-identical, grep whether ANY
  test asserts the OLD guard's raise (`pytest.raises(match=<guard msg>)`).
  If none — common for internal "belt-and-suspenders" consistency
  asserts (the SN `_require`/`_refuse_starting_direction` biconditional
  had ZERO negative tests across 7 sites) — the retirement silently
  drops a REAL-but-unpinned safety net. The plan MUST WRITE the negative
  test the guard never had, pinned to the NEW consolidated check, or the
  swap loses teeth invisibly. A `pytest.raises` teeth-gate for such a
  check MUST `match=` the SPECIFIC consolidated-check message: a
  downstream crash on the same malformed input (e.g. `precompute_psi_
  state(None)` → AttributeError) satisfies a bare `pytest.raises(...)`
  and gives a FALSE-green teeth claim (the xfail-strict cousin). Pair
  with the positive CONTROL (the correctly-built composite solves) so
  the raise is attributable to the mismatch, not a broken path.

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

---

## L16 — Verifying a NEW *consumption mode* of a shipped per-cell closure algebra (`assemble()` = the 3rd mode beside `solve`/`apply`) + the relocation-is-behavior-free gate

When a carve adds a THIRD consumption mode of the ONE closure algebra
(SN/diffusion `solve`=sweep + `apply`=matvec already share the per-cell
coefficient source; the new `assemble()` emits the SAME coefficients as
(row,col,value) into a global scipy-sparse `LinearOperator`), and relocates
its module DOWN a layer (`sn/spatial`→`transport/spatial`), six reusable
test-design facts (Phase-2 spatial-substrate + assembly campaign,
`.claude/plans/assembly_mode_crosswalk.md`):

- **The ONE-SOURCE PROOF is the campaign's whole point (the Phase-F twin-path
  guardrail).** `assemble` MUST consume the same coefficient source `solve`/
  `apply` do (DD: `cell_balance.cell_balance_for_streaming`; LD:
  `_ubld.assemble_ubld`+`d1_closed_form`; diffusion:
  `operators._interior_conductance`/`_boundary_closure`). The proof =
  a monkeypatch **sign-flip in the SHARED source** must red BOTH the new
  equivalence gates AND the EXISTING sweep/matvec suites. If ONLY the new gate
  reds → a parallel stencil (twin path) exists → **STOP, fix, log ERR-NNN**.
  Classify every mutation shared-source (reds both) vs assembly-local (reds
  only new: DOF-transposition, wrong-neighbor, walk-order-P). **DD is the
  HIGHER twin-path risk than LD**: DD has NO dense block (its matvec is the
  fused scalar `residual_kernel_batch` → a fresh stencil emission is tempting);
  LD's apply already does `einsum(assemble_ubld["A"], psi)` so LD-assembly is
  scatter-of-EXISTING-blocks. The M1 sign-flip on the DD shared source is the
  single highest-value mutation.
- **Sparse-order ≠ apply-order ⇒ NO gate is 0-ULP; the independence is a
  FEATURE.** `assembled@x` sums each row in CSR order; `apply(x)` in the
  native `np.diff`/`einsum`/left-fold order → G1 (`assembled@x≡apply(x)`) is
  **nulp/`rtol=1e-11`** (bound = row-bandwidth×eps), NEVER `array_equal`
  (vv §bit-identity crit-3). `solve_triangular` (scipy→LAPACK `dtrtrs`) is a
  structurally-INDEPENDENT forward substitution vs the ORPHEUS scalar sweep
  recurrence → G2 (`solve_triangular(P·L·Pᵀ,q)≡sweep.solve(q)`) is **rtol**,
  and that independence EARNS G2 its L2-cross-check status + **discharges the
  sweep-inverse contract question object-level** (here #284). G2's
  triangularity leg (`triu(P·L·Pᵀ,k=1)==0`) is a structural EXACT zero.
- **The `as_matrix`-DELEGATES-to-assembly TAUTOLOGY trap (proactive refute,
  L10).** Once `as_matrix` delegates to densified assembly (the R2 design), the
  naive `op.as_matrix()≡op.assemble().to_dense()` is VACUOUS (assembly vs
  itself). The fuller-view-oracle gate MUST force the **RETAINED probing
  fallback** (`_as_matrix_by_probing`, the pre-carve apply-to-basis loop) ≡
  `assemble().to_dense()`. Keep EXACTLY ONE probed≡assembled pin per family as
  the permanent oracle (an explicit keep decision, coding-standards fuller-view
  exception) — the diffusion loss (already has the `as_matrix`≡hand-posed-FD
  stencil l0 gate) + one DD slab block. Probing `A[i,j]=apply(e_j)[i]` extracts
  ONE coefficient (no summation) but via apply-arithmetic ≠ direct emission →
  still nulp, not 0-ULP.
- **Config: the DOF-transposition + face-pairing mutations are transpose-BLIND
  on a SYMMETRIC operator (Mode-12).** A row/col swap on `A=Aᵀ` is invisible →
  the fixture MUST be het + **asymmetric SigS** + **non-uniform h** so `A≠Aᵀ`
  and the transposition/face-swap is observable (reuse the exact
  `tests/diffusion/test_operators.py` fixture — it already pulls every lever
  with the comment "asymmetric Σ_s so a transpose is observable"). G1's `x`
  MUST be **non-flat fixed-seed random** with non-zero inflow (flat x nulls the
  DD streaming coupling, §0.6/H2). ≥2G always (anti-#3); exclude 1G explicitly,
  ship NO 1G "structural smoke".
- **G1 is object-level ONLY as a matvec, NEVER a derived scalar.** `Δx=0` for a
  FIXED implementation bug `Δ=A_asm−A_apply` is caught by a single GENERIC
  random x w.p.1 (kernel measure-zero); K≈8 fixed-seed random = redundancy;
  G3 (probed≡assembly on the full matrix) is the exhaustive object pin on the
  anchor member. NEVER gate `keff(assembled)≡keff(apply)` — a similarity/
  transpose Δ is spectrally invisible (Mode-12, the #226 step-5b overclaim).
- **The out-of-scope-defect CHARACTERIZATION gate flips LOUD, opposite to
  xfail.** Curvilinear assembly is OUT (blocked on the #282 lagged-ψ½ pole
  seed); characterize it by assembling the spherical walk and POSITIVELY
  asserting `np.any(triu(P·L_sph·Pᵀ,k=1)≠0)` (the defect EXISTS) with an
  actionable message. NOT xfail: xfail-`strict=False` flips silently to xpass
  for a not-yet-landed FEATURE; here we want a LOUD RED forcing a rewrite when
  the defect is FIXED (promote to a G2 triangular gate). Check `error_catalog`
  for a #282 ERR → `catches()`; else cite the issue.
- **RELOCATION is behavior-free BY CONSTRUCTION → lean on existing gates, argue
  AGAINST a new snapshot/hash.** A pure `git mv`+import-rewire changes no
  numerical-code byte, so a NEW output-hash/snapshot is redundant with the
  EXISTING `-W error::DriftWarning` snapshot suites (which already pin byte-id;
  an iterated snapshot is the WRONG bit-id gate anyway, L7). The walls are:
  existing DriftWarning suites green + Sphinx `-W` clean (a `:mod:`/`:class:`
  xref at the old path renders as plain text with NO `-W` warning → ALSO)
  grep-completeness `grep -rn "<old.path>" orpheus/ tests/ docs/` returns ZERO
  (L15 regrep) + importer-suite nodeid-set identical modulo the moved files'
  own paths. Tests travel WITH their modules (relocation ≠ retirement — move
  behavioral AND import-path smokes; markers key on nodeid, `registry.clear()`
  per-collection so no surgery).

How to apply: for a new-consumption-mode-of-a-shipped-algebra carve, the
§0.5 named risk is THE one-source proof (shared-source sign-flip reds BOTH new
AND existing; DD is the twin-path-risk cell). Tolerances: G1 nulp / G2
rtol-via-independent-scipy / structural-EXACT triangularity; NEVER `array_equal`
(sparse-order), NEVER a scalar functional (Mode-12). Force the retained probing
path for the fuller-view-oracle (as_matrix-delegates tautology). Mandatory
asymmetric-SigS+non-uniform-h (transpose observability) + non-flat random x.
Characterize an out-of-scope defect with a positive assert-defect + loud-flip
message, not xfail. Relocation: existing DriftWarning + Sphinx-W + regrep-ZERO,
no new hash. The consumption-mode sibling of L11 (axis-transpose of a shipped
reduction) / L12 (fast-path→composed fold); `vv §bit-identity / Mode-8/12 /
anti-#2,#3,#4,#11`; mirrors L7/L10/L15.

## L17 — Verifying the reverse-SCAN transpose-solve (#280) + its apply-loop unification: the assembled-Mᵀ oracle, the self-referential-relocation trap, and the retired-vocab reconciliation

The #280 carve builds the 4th face of the `(L+C)` 2×2 — the reverse-DAG
transpose-solve — as a reverse-SCAN coherent with the forward Blelloch
`_run` (NOT the reverse-LOOP the a3 spec assumed), and first unifies the
apply-loop family `{_apply_walk (fwd matvec), loss_action_transpose (adj
matvec)}`. Full recipe `a3_reverse_scan_transpose_verification.md`; the
durable test-design facts:

- **The operator VOCAB is stale — reconcile against the #226 typed-predicate
  world FIRST (L10 refute-the-premise).** The a3 spec's
  `InvertibleOperator.solve_transpose`/`CAP_SOLVE_TRANSPOSE`/`_AdjointOperator.solve`/
  `CAP_SOLVE` is RETIRED (L13/L15 frozenset→predicate landed). CURRENT
  primitive = **`SweepOperator.apply_transpose`** (`SweepOperator=(L+C).inverse()`;
  the reverse-scan IS the inverse operator's transpose matvec), gated by
  **`SweepOperator.is_adjointable`→True**. The metric transpose-solve
  `A.inverse().H.apply(b)=G⁻¹·apply_transpose(G·b)` falls out of the EXISTING
  `_AdjointOperator.apply` FOR FREE — NO `_AdjointOperator.solve` is minted;
  `_AdjointOperator.inverse()`+`is_invertible` close the coherence
  `A.H.inverse()≡A.inverse().H`. Design gates against the LIVE surface, never
  the retired string tags (`InverseWrapMixin` docstring already names #280).
- **The existing 0-ULP `array_equal` matvec canaries are SELF-REFERENTIAL →
  necessary-not-sufficient for a RELOCATION.** `test_removal_form_matvec_sweep`'s
  fwd/adj `array_equal` gates compare `op.apply[_transpose]` against the SAME
  `loss_action[_transpose]` on a FRESH `StreamingOperator(σ_r)` — a
  reduction-tree "override-owns-it vs leaf-sum-leak" discriminator, NOT a value
  oracle. A relocation moving BOTH the SUT path AND the reference path (same
  relocated code) leaves them GREEN even if values shifted (the twin/Mode-11
  hazard). The genuine 0-ULP relocation proof is a **FROZEN pre-carve baseline
  snapshot** (`assert_regression --capture-baseline`, both orientations,
  slab/sphere/cyl 2G) — structurally independent because captured BEFORE the code.
- **The new shared 1-D loop frame needs its OWN structural pin (there is NONE
  today).** Multi-D has `test_one_octant_walk` (runtime `_interior_walk` spy +
  AST tripwire on `{is_solve,is_apply,is_matvec}`); the 1-D `_OneDimScanWalk`
  has no spy (its `_apply_walk`/`loss_action_transpose` were verbatim
  relocations, three separate skeletons). 2.5a CREATES the shared frame → mint
  the ORIENTATION sibling: a Mode-11-WRAP runtime spy (both orientations route
  through the ONE method, count>0) + an AST tripwire forbidding
  `is_adjoint`/`is_forward`/`is_transpose`/`is_reverse` as identifiers, demanding
  an orientation-carrying OBJECT (like `_ApplyOperands`/`_SolveOperands`).
- **The assembly mode gives a NEW structurally-independent transpose oracle —
  Cartesian ONLY.** `permuted=block.as_matrix()[ix_(order,order)]` is
  walk-order lower-tri; the reverse-solve oracle is
  `solve_triangular(permuted.T, b[order], lower=False)` (LAPACK back-sub,
  shares NO code with the ORPHEUS scan). CATCHES the wrong TRANSPOSED-SCAN-
  COEFFICIENT (a'/b' of the reversed affine recurrence — scan≠loop) + gives the
  exact triangularity certificate + discharges "reverse-substitution IS
  sweep_transpose" object-level (transpose analog of #284). But assembly is
  Cartesian-only (`SNMesh.streaming` gate) → **the SPHERE μ-reversal/mirror bug
  is INVISIBLE to it; the dense-apply G2 `np.linalg.solve(M.T,b)` (M from fwd
  apply) stays the KEYSTONE on the mandatory sphere leg.** LD caveat: if
  `loss_action_transpose` is DD-only (no moment tail), the LD assembled-Mᵀ
  oracle is moot — verify LD-adjoint scope first.
- **Reverse-SCAN's NEW modes over the reverse-LOOP:** (i) wrong transposed scan
  coeffs → G2 dense-Mᵀ + assembled-Mᵀ; (ii) two-denom seam (÷V apply vs ×V
  `_run` scan, #242) → G1 round-trip + G2; (iii) curvilinear μ-level-coupled scan
  reversal + `angular_adjoint` + Carlson mirror → G1/G2 SPHERE leg (slab NULLS);
  all object-level — NEVER credit `sweep_transpose` on `k*`/norm (`eig(Kᵀ)=eig(K)`
  is adjoint-blind, Mode-12; memory "NEVER keff(asm)≡keff(apply)").
- **Scaffold: G3 full-loss `(L+C−S−B)` reciprocity lands PRE-carve** (S† at
  `15185e5`; extend `test_g_adjoint_reciprocity` L+C−B→L+C−S−B, asymmetric SigS
  ≥2G, per-group one-hot φ vv L27) — it hardens the `apply_transpose` surface
  2.5a rebuilds, giving the adjoint bit-identity a composite-level canary.
  #282-fix (CONDITIONAL, unruled): if the augmented ψ(·,μ=−1) state lands in
  2.5, the `test_282_characterization` loud-flip fires (L16) and its successor
  is a spherical G2 gate on the augmented system (seed rows first, back edge →
  forward edge; teeth = swap the seed coupling direction).

How to apply: reverse-DAG transpose-solve or apply-loop-unification carve →
FIRST reconcile the retired `CAP_*`/`solve_transpose` vocab to the live
`SweepOperator.apply_transpose`/`is_adjointable`; frozen-baseline the
relocation (the array_equal canaries are self-referential); mint the
orientation-OBJECT AST tripwire + 1-D spy (none exists); the assembled-Mᵀ
oracle is a real Cartesian L2 cross-check but the SPHERE dense-apply G2 is the
keystone (assembly is Cartesian-only); object-level gates ONLY (Mode-12);
G3 reciprocity is pre-carve. The transpose/inverse-axis sibling of L13/L14/L15
(the operator-taxonomy family) + L16 (the assembly consumption mode);
`vv §bit-identity / Mode-8/11/12 / anti-#2,#3`.

## L18 — Verifying a CARRIER AUGMENTATION's ψ½ block (#282 route (a): the direct μ=−1 starting-direction seed as a typed per-level block). ⚠ SHARPENED 2026-07-06: the "zero-G-weight DOF → reciprocity IDENTICALLY blind → encode the blindness as a positive pin" claim was REFUTED by the derivation-of-record — the STATE metric is G_sd=V_cell (SPD), NOT the ghost zero, so reciprocity CLOSES Mode-12 (a seed-row flip REDS). Still-valid: closed-form-direct-solver pillar, prove-the-block-is-CONSUMED, cold-residual acceptance number, "cylinder unmoved" is a TESTED assertion

The #282 fix adds a typed per-starting-direction-level ψ½ block to the
curvilinear composite (sphere 1 level, cyl n_polar; slab/cart ABSENT) and makes
`solve` march the seed rows FIRST via `carlson_inward_sweep_from_source` on the
TRUE source's q½ block (the lag dies; the sphere solve becomes a genuine
triangular inverse). Full spec `a3_solve_transpose_verification.md` §16. Durable
test-design facts (the sub-phase sibling of L16/L17):

- **⚠ CORRECTED: a carrier DOF's Hilbert STATE metric is set by its OPERATOR
  ROLE, NOT by its angular-INTEGRATION weight — conflating them gives the WRONG
  (ghost) metric and a Mode-12 FALSE-GREEN** (derivation-of-record
  `starting_direction_metric_gauge_derivation.md`, 2026-07-06; this REPLACES my
  earlier "zero-weight DOF is reciprocity-invisible / encode the blindness as a
  positive pin" claim). The Carlson α₁ᐟ₂=0 is the angular THROUGH-FLUX weight
  (correctly excludes ψ½ from `φ=Σw_n ψ_n` — an OPERATOR coefficient); the STATE
  metric is set by ψ½'s role as a first-class RADIAL field (its self-block A_ss
  is a BANDED radial transport operator `μ∂_r+σ`, not a grazing angular face —
  diag_gsd_05). So **G_sd = V_cell (SPD)**, mirroring `G_bulk=V·w`; since
  `apply_transpose≡Aᵀ` EXACTLY, reciprocity is GAUGE-FREE among ALL SPD G_sd.
  **G_sd=0 is the single FORBIDDEN value**: it puts the seed in null(G), severing
  the seed→bulk A_bs coupling from `A.H=G⁻¹AᵀG` — a WRONG adjoint the instant the
  seed carries data (production random-seed defect 1.3e-2), green ONLY because the
  old gate fed a present-but-ZERO seed. **With SPD G_sd, reciprocity DOES catch a
  seed-row transpose error (Mode-12 CLOSED)** — INVERT any gate built on the old
  blindness (rename `_blind_to_` → `_catches_`; assert flip REDS not stays-green).
- **A Mode-12-CLOSURE reciprocity gate needs BOTH legs — control-holds AND
  mutated-reds; the mutated leg ALONE is fooled by a broken baseline.** Flip
  `_seed_rows_forward` (forward A_ss) but NOT `_seed_rows_transpose` (A.H's
  independent reverse mode) → an apply-vs-transpose inconsistency the SPD metric
  sees. Assert (1) UNMUTATED nonzero-seed reciprocity HOLDS `<1e-12` (baseline is
  the honest V_cell adjoint) AND (2) the flip REDS `>1e-6`. WITHOUT leg (1),
  reverting G_sd→0 gives flip-defect 0.107>tol TOO (baseline broken, flip
  invisible — a broken baseline MIMICKING "caught") → false green. Measured dense
  sweep {0/id/V_cell/V·w}×{zero/random/flip}: g_sd=0 → random 0.107 / flip 0.107
  (BLIND, unchanged); V_cell → random 1e-16 / flip 1.000 (CLOSED); faithfulness
  `‖G⁺TG−A.H‖`=2.8e-14. The gauge is INERT (control leg c): corner V[-1]→2·V[-1]
  keeps reciprocity green AND the bulk solve bit-identical (observable-level).
  Feed the independent `_g_inner` seed term from PRODUCTION's
  `starting_direction_space.inner_product_weights` (the gauge must MATCH A.H's,
  else false-red — the seed metric verifies no VALUE, only SPD-ness, so this is
  NOT reference contamination; the load-bearing content is the TRANSPOSE, caught
  because the error is in Aᵀ not G). Canonical Mode-12 gate = the promoted
  `diag_gsd_02` dense sweep; production-path corroboration = the inverted
  `test_282::test_mode12_g_reciprocity_catches_a_seed_row_flip`.
- **The direct-solver structural pillar is the closed-form exponential ODE, NOT
  the recurrence hand-trace.** The μ=−1 equation `−dφ/dr+σφ=q` has
  `φ(r)=q/σ+(φ_R−q/σ)e^{−σ(R−r)}`; DD converges O(Δr²) to it. The EXISTING
  Carlson pins hand-trace the SAME DD recurrence (procedural, NOT structural —
  the ERR-032 shared-upstream hazard) and route through the PROXY `Q̄=σφ₀/Σw`,
  never the true-source direct entry. Add the convergence pin on an ARBITRARY
  `Q_bar`; the GRADED-mesh row is MANDATORY (uniform is blind to a `dr[k]` index
  drift, Mode 5). The q½ source fold needs the anisotropic 2-term
  `P_ℓ(−1)=(−1)^ℓ` hand-check manufactured FIRST (an all-isotropic suite is
  silently blind to a dropped/mis-signed φ_ℓ≥1, §0.6).
- **Prove the new carrier block is CONSUMED before crediting any gate (term
  activation, Mode-11-adjacent).** Zero the source's q½ block → the solve MUST
  move (`solve(b_with_q½) ≠ solve(b, q½=0)`); if equal, the whole augmentation
  is INERT (silent coverage loss). Pair with a Mode-11 WRAP sentinel on
  `carlson_inward_sweep_from_source` (count>0 on sphere solve, ==0 on slab).
- **The cold-residual is THE acceptance number, promoted from a diagnostic —
  a LAG-DEATH classifier, distinct from seed correctness.** `‖A·solve(b)−b‖/‖b‖`
  over the FULL augmented field, COLD: sphere 5.18e5 → `<1e-11` (slab-level, no
  iterate-threading, no sphere slack — the fix EARNS the tight tol); slab/cyl
  controls already ~7e-16 must STAY. It certifies solve↔apply CONSISTENCY (lag
  dead), NOT coefficients (the direct-solver pin + Euclidean `Mᵀ` own those).
  Companions: seed-INSENSITIVITY (Δ 4.57e-2 → bitwise), the decisive classifier
  (`b=A.apply(psi0)`; cold solve recovers `psi0` → the fixed point did not move),
  the coarse `level_symmetric` S4 end-to-end (SI-NaN/Krylov-negative → finite+
  positive; a physicality gate, NOT precision NOT keff), the c=0 no-outer-
  iteration degenerate.
- **"Curvilinear re-baselines" is too coarse for the CYLINDER — assert-unmoved
  FIRST.** The cyl seed is weightless (α=0 AND the per-level α-dome telescopes,
  the #282 0.0-bit row) → the frozen `walk_matvec_cyl_2g` MUST stay `array_equal`
  (a DIFFERENT ψ½ carrier value gives the SAME physical output iff the seed is
  truly annihilated). A cyl move HALTS the carve (FP-reorder → re-capture WITH a
  3-criteria record + "weightless was too strong"; or a leak-bug). Only the
  SPHERE re-captures (principled, [[feedback-principled-over-bit-identical]]);
  the slab stays bitwise. Sharper than the roadmap's blanket "curvilinear
  re-captures" — it pins the seed-weightless INVARIANT a silent re-capture loses.
- **Sequencing (R10: 2.5d BEFORE 2.5b) sets the gate split.** The fix lands
  first → the 2.5b sphere G1/G2 `Mᵀ` chain is written EXACT single-pass (rtol
  1e-11, NO iterate-threading — SUPERSEDES the a3 §0.3 loose-1e-9). 2.5d's own
  successor to the loud-flipped `test_282_characterization` is triangularity
  ONLY (`triu(PᵀM_augP,1)==0` on the matvec-probed AUGMENTED matrix, seed-rows-
  first order — NO assembly-emitter extension); the LAPACK-≡-sweep leg waits for
  the 2.5b reverse-scan.
- **Framing-agnostic by construction.** The ψ½ block's storage (F1 third block /
  F2 angular-boundary-trace / F3 zero-weight extra ordinate) is under concurrent
  design — write gates against the carrier's PUBLIC contract (presence keyed by
  curvature both-ways, to_flat round-trip, add/scalar closure, zero-metric-weight,
  FULL-field residual) so they hold in ALL three; `# FRAMING:` sub-notes ONLY on
  storage-layout assertions. All three give the seed zero G-weight → the
  reciprocity-blindness note is framing-independent.
- **The ℓ=0-ONLY q½ source fold SILENTLY BREAKS the MMS fixed-source — the
  eigenvalue wall is BLIND (migration wave 3, 2026-07-04, concrete war-story of
  the anisotropic-`P_ℓ(−1)`-first warning above).** An MMS fixed-source is
  INHERENTLY anisotropic: the streaming operator `Ω·∇ψ` manufactures a ℓ=1
  source even for an ISOTROPIC trial `ψ=A(r)` (`Q = μA'(r) + σ_t A(r)`). A
  `from_angular_source` fold that folds only ℓ=0 (`q½=½q₀`, its documented
  "honest scope: production sources are isotropic") seeds the pole march with a
  WRONG q½ → the WHOLE inward march corrupts (interior AND pole wrong; sphere L2
  stuck ~6.7e-2 vs the frozen O(h²) 3.7e-3; NEGATIVE convergence order). The
  EIGENVALUE path is BLIND (isotropic fission/scatter → ℓ=0 fold EXACT), so a
  green k_inf / SI≡Krylov wall MISSES it; ONLY the anisotropic MMS FIXED-SOURCE
  pole-cell gate reds. Two corollaries for the migration test-architect: (1)
  before trusting an MMS fixed-source as an ABSORB-gate under a ψ½-source carve,
  verify the fold's MOMENT REACH ≥ the MMS source's anisotropy — a "production
  sources are isotropic" scope note is FALSIFIED by the MMS's own streaming
  source; (2) a pole-cell CONVERGENCE gate cannot be re-baselined (no stored
  value) — a negative order under the carve is a STOP-and-report, never a
  re-capture. The isotropic-source solve (reflective `q/Σ_t`) stays green, which
  is exactly why the break hides until the anisotropic gate runs.

How to apply: a carve that ADDS a zero-measure/zero-weight DOF to a
metric-carrying composite (a starting-direction seed, a ghost cell, any
α=0 phase-space DOF) → the G-adjoint reciprocity gate is IDENTICALLY blind
(Mode 12); constrain the DOF by the EUCLIDEAN transpose oracle + a closed-form
direct-solver pin + the full-field round-trip, encode the blindness as a
positive pin, and PROVE the DOF is consumed (zero its source → the solve moves)
before crediting anything. Promote the seed-lag diagnostics (cold-residual,
seed-insensitivity, decisive classifier) to acceptance gates. Assert-unmoved
FIRST on the config where the DOF is weightless (a silent re-capture discards the
invariance). The direct-solver's structural reference is the closed-form ODE, not
the recurrence hand-trace (ERR-032). Sub-phase sibling of L16/L17;
`vv §Mode-10/11/12 / anti-#8 / §bit-identity`.

## L19 — Verifying the adjoint-inverse SWAP-LAW wiring (#280 2.5c, `A.H.inverse() ≡ A.inverse().H`): Gate 1 is a FORWARD-only reciprocity structurally independent of the reverse-scan, `b` MUST be bulk-only source-carried, Gate 2 is a bit-identical OBJECT identity, the `.H`≠Euclidean "slab stays GREEN" discriminator is FALSE, and the predicate flip MUST propagate to the capability-survival CONTRACT

The 2.5c terminal of the #280 reverse-scan sub-campaign wires the swap law
`(A*)⁻¹ = (A⁻¹)*` as an OBJECT IDENTITY: `_AdjointOperator.inverse() =
inner.inverse().H`, `SweepOperator.apply_transpose = inner.solve_transpose` (the
2.5b reverse-scan), `SweepOperator.is_adjointable` flips True over the
`InvertibleOperator` arm; the metric adjoint-solve `A.H.inverse().apply(b) =
G⁺·solve_transpose(G·b)` falls out of the EXISTING `_AdjointOperator.apply` FOR
FREE. Gate file `tests/sn/operators/test_inverse_adjoint_coherence.py` (helpers
reused from `test_loss_transpose_solve.py`). Durable gate-design facts (the
wiring-proof sibling of L13–L18):

- **Gate 1 (forward-matvec G-reciprocity `⟨Aψ,x⟩_G=⟨ψ,b⟩_G`, `x=A.H.inverse().apply(b)`)
  is structurally independent of the reverse-scan BY CONSTRUCTION.** The
  verification arithmetic is FORWARD `A.apply` + the metric inner product ONLY —
  it NEVER calls `apply_transpose`/`solve_transpose`; the reverse-scan lives
  solely inside the SUT `x`. `⟨Aψ,(A*)⁻¹b⟩_G = ⟨ψ,A*(A*)⁻¹b⟩_G = ⟨ψ,b⟩_G` holds
  iff `x` genuinely solves `A*x=b`, paired against the independently-verified
  forward operator (the strongest available ground: forward apply is the SoT, the
  metric is pinned by `test_g_adjoint_reciprocity`). This is the elegant way to
  gate a reverse solve WITHOUT the reference sharing its code.
- **`b` MUST be bulk-only source-carried; a random FULL `b` falsely reds Gate 1/3
  (~1e-1) even UNMUTATED.** solve/apply (and transposes) are genuine inverses only
  OFF the boundary/seed outflow FREE-DOF slots (#284) AND off the metric null
  space (tangential trace |Ω·n|=0 ⊕ the all-zero-ghost-metric seed). A bulk-only
  `b` (zero boundary/seed) lies in the range of `A*` → Gate 1/3 clean to ~1e-16 on
  slab/sphere/cyl. ψ (the test vector) stays FULL (bulk+trace random) so the trace
  coupling + trace metric ARE exercised in the pairing; seed present-but-ZERO
  (Mode-12 regime). The well-posed round-trip is `(A*)⁻¹∘A*` (=I on bulk);
  Gate 3's `A*∘(A*)⁻¹` is =I only on the source-carried bulk — compare BULK, never
  the full field.
- **Gate 2 (swap law) is a BIT-IDENTICAL object identity, NOT an rtol check.**
  `A.H.inverse()` returns LITERALLY `A.inverse().H`, so both spellings build+run
  the SAME object graph on the same `b` → `assert_array_equal` (0 ULP), stronger
  than the spec's allowed rtol=1e-12. Compare bulk ⊕ boundary ⊕ seed.
- **The M-ADJ-metric "slab stays GREEN" discriminator is FALSE — flag it, don't
  encode it.** The spec predicted skipping the `G⁺`/`G` wrap in
  `_AdjointOperator.apply` reds Gate 1 only on sphere/cyl (slab "metric trivial").
  EMPIRICALLY the composite metric `G = V·w_n` (bulk) ⊕ `|Ω·n|·w_n` (trace) is
  non-trivial on EVERY geometry (non-uniform mesh + the bulk-vs-trace weight
  mismatch that streaming couples), so `A*=G⁻¹AᵀG≠Aᵀ` everywhere and Gate 1 reds on
  slab too (slab 0.33 / sphere 0.19 / cyl 0.13). The `.H`≠Euclidean claim is TRUE
  universally (the metric wrap is load-bearing on every geometry — STRONGER
  coverage than the intended discriminator), just NOT a slab-vs-curvilinear split.
  NEVER design a `.H`≠Euclidean gate as a slab-green discriminator on a composite
  bulk⊕trace metric — the trace weight `|Ω·n|` ≠ the bulk weight `V` on ANY
  geometry (sibling of L1: verify the CONCRETE metric before assuming slab nulls a
  term). The mutation still has teeth (reds all geoms); assert RED on all three.
- **`A.H.inverse()` needs the Design-C `invertible(A.H)` TypeGuard bridge for
  pyright.** `.inverse()` lives on `SupportsInverse`; `.H` returns the base
  `LinearOperator` (which erases it). Wrap: `ah=A.H; if invertible(ah): return
  ah.inverse()` (`invertible` is a TypeGuard — POSITIVE branch only; a trailing
  `pytest.fail` is NoReturn). The bridge doubles as a live
  `_AdjointOperator.is_invertible` pin and reads predicates only (no solve) so it
  does NOT perturb the Mode-11 counters. `A.inverse().H` needs NO bridge (`A` is a
  concrete `InvertibleOperator`).
- **Mode-11 wrap sentinel: wrap `InvertibleOperator.solve_transpose` (>0) AND
  `.solve` (==0).** `A.H.inverse().apply(b)` MUST enter the reverse-scan and NEVER
  the forward solve; the routing mutation (`SweepOperator.apply_transpose →
  inner.solve` forward) reds BOTH counters. All geometries give st=1/solve=0
  (sphere seed is single-call exact post-2.5d — no iterate-threading).
- **The predicate flip is SURGICAL and MUST propagate to the capability-survival
  CONTRACT in the SAME landing.** `SweepOperator.is_adjointable` flips True over
  `InvertibleOperator`; the sibling wrap-delegates `InverseOperator`/`GreenOperator`
  STAY False (no override — keyed on `isinstance(inner, InvertibleOperator) and
  inner.is_adjointable`). **The Part-A landing that flipped it BROKE TWO existing
  forward canaries (both the same flip; blast radius WIDER than any one contract
  table — the L14/L15 grep-audit lesson):** (a)
  `test_capability_survival.py::test_inverse_adjoint_contract_keystone_v2_sn` row
  `((L+C).inverse(), True, False, INVERTIBLE)` (`adj=False`→`True`); (b)
  `test_inverse_operator_equivalence.py::test_inverse_returns_sweep_operator_surface`
  `assert inv.is_invertible is True and inv.is_adjointable is False` (→`is True`).
  A predicate-flip carve's retirement audit MUST grep EVERY `is_adjointable is
  False`/`adj=False` assertion on the flipped type (NOT just the contract-table
  rows) and update them in the same commit; a single-file gate-authoring session
  downstream cannot fix them (scope) — report as required Part-A completion items.

How to apply: verifying an adjoint-inverse / swap-law wiring (`A.H.inverse()` /
`A.inverse().H`) → Gate 1 is the FORWARD-only G-reciprocity (structurally
independent of the reverse-scan, the keystone) on a bulk-only source-carried `b`
with a FULL ψ; Gate 2 is `array_equal` (object identity); route `A.H.inverse()`
through the `invertible()` TypeGuard bridge; Mode-11-wrap the reverse+forward
inner readers. NEVER design a `.H`≠Euclidean gate as a slab-green discriminator —
the composite bulk⊕trace metric is non-trivial on every geometry. A predicate flip
MUST propagate to the capability-survival contract table. Terminal wiring-proof of
the #280 reverse-scan chain (L13→L17 vocab, L18 carrier);
`vv §Mode-11/12 / L11 structural-independence / L1 config-blindness`.

---

## L20 — Verifying a STRUCTURE-ONLY ABC re-parent + stringly-key→typed-key `[K]` generalization (intended bit-identical): the crux is a MULTI-index synthetic fixture (the shared layout is exercised in production at ONE key-index → the generalized index term is never activated), the metric is a REGRESSION-GUARD (not a new gate) with an ERR-067 "NO metric on the ABC" tripwire, and the load-bearing MRO check is sibling-NOT-child

The C2 `FaceField(Field)` carve (design `facefield_codim1_design.md` §5): introduce a
parent ABC that owns a shared flat-buffer discipline ONCE, re-parent three sibling
field types under it (`AngularBoundaryField`/`ScalarBoundaryField`/
`StartingDirectionField`), and generalize `FaceLayout`/`FaceSlot` from `str` keys to a
typed key `[K]` so the pole space (`(level,sign)` keys) stops re-implementing
`_leg_offset`/`cells_slice`/`corner_slice` and shares the ONE impl with the trace
(`str` keys). Intended bit-identical; nothing numeric moves. Five reusable facts:

- **The crux is a MULTI-INDEX SYNTHETIC fixture — the shared layout is exercised in
  PRODUCTION at exactly ONE key-index, so the generalized index term is never
  activated (config-blindness, sharpening of L1/L11).** `StartingDirectionSpace`'s
  offset formula is `(2·pos + sign_pos)·per_sign` where `pos = levels.index(level)`,
  but production (sphere-GL) carries a SINGLE seed level → `pos ≡ 0` always → the
  level-index term is DEAD in every existing pin (`test_starting_direction_carrier.py`
  builds only `_sphere`). A typed-key `FaceLayout[(level,sign)]` bug in the
  level-position mapping is INVISIBLE to the whole suite. Fix = a direct
  `StartingDirectionSpace.for_levels((0,2,5), ng, nx, cell_volumes=…)` synthetic
  (valid — the ctor guard only rejects empty/non-ascending), forcing `pos ∈ {0,1,2}`.
  The LOAD-BEARING mutation is `2*pos → 1*pos`: it REDs the multi-level golden at
  `pos≥1` while the single-level sphere STAYS green — the asymmetry IS the
  config-blindness evidence (same shape as L11's product-RED/LS-GREEN). The golden is
  HAND-COMPUTED (`(2·pos+sign_pos)·per_sign`), NOT SUT-captured, so it holds
  independently of the carve; pin `.offset`+`cells_slice`/`corner_slice` start/stop AND
  slice CONTENTS via an `np.arange` fingerprint fill (a transposed/mis-offset slot
  diverges; a bare shape check is blind to a transpose). The `str` regime is already
  golden'd by `test_face_layout.py::TestFromNamedShapes` (hardcoded `.offset==0/16/32/56`).

- **The metric is a REGRESSION-GUARD, not a new gate — with an ERR-067 tripwire that
  the ABC stays metric-free.** "Structure-only, metric descends per-leaf" means the
  carve MUST NOT touch the per-leaf metrics (trace `|Ω·n|·w`, pole `V_cell`). Lean on
  the EXISTING object-level metric gates kept green + unmodified
  (`test_starting_direction_metric.py::test_shipped_metric_block_values` `atol=0.0`;
  `TestA4SeedStateMetricVcell`; the trace 0-ULP `test_one_kernel_reproduces_trace_metric_at_0ulp`).
  ADD ONE new tripwire: a mutation that gives the ABC a UNIFORM face metric
  (`(1−µ²)·w ≡ 0` at the pole — the retired ghost!) must RED the byte-exact `V_cell`
  gate. That is the ERR-067 guardrail: "NO metric on the `FaceField` ABC" is precisely
  what stops `G_sd=0` from returning. Pole gates are Mode-12-SAFE only because they
  assert the weight ARRAY object (`assert_array_equal` on `cells_view(weights)`), never
  a scalar functional; every pole gate feeds a NONZERO seed + a control leg (the twin
  bulk/trace-only composite) — a zero seed sits in the ghost's invariance group (L18).

- **The load-bearing MRO check is SIBLING-NOT-CHILD.** After `X` and `Y` both become
  children of the new parent `P`, the dispatch that discriminates them
  (`isinstance(self.boundary, BoundaryField)` / `isinstance(self.starting_direction,
  StartingDirectionField)` in `full_field.py:274/280`) survives ONLY if the re-parent
  keeps them SIBLINGS: `isinstance(sd_field, BoundaryField) is False` and
  `isinstance(boundary_field, StartingDirectionField) is False` MUST stay — assert both
  (the cross-slot contamination guard). The break mode to gate against: collapsing the
  INTERMEDIATE name (`BoundaryField`) INTO the new parent `P` breaks
  `isinstance(x, BoundaryField)` at the consumer. Also gate `P in cls.__mro__` (between
  child and `Field`), `P` exported from `__all__`, `isinstance(x, Field)` stays True for
  every leaf (the `functional.py` consumer), and — the subtle seam — the flat-buffer
  validation (moving UP from `BoundaryField.__post_init__` onto `P`) fires EXACTLY ONCE
  through the `super().__post_init__()` chain (a moved-up check must not double-run or
  skip; pin the layout-less-space raise still fires with the right leaf type-name).

- **`to_flat`/`from_flat` round-trip is bit-identical BY CONSTRUCTION iff the flat
  offsets are unchanged** — it reads `.values` (the flat buffer), never `FaceLayout`
  directly. So the offset gate DISCHARGES the round-trip; no new serialization
  mechanism, just keep the existing `TestA2FlatRoundTrip` green (incl. its dropped-seed
  `_two_block_to_flat` mutation tooth) and extend the flat-total assertion to the
  multi-level synthetic.

- **Retest = the offset/metric/MRO seams, not the whole solver.** A structure-only
  relocation is behavior-free by construction; the proof is (offsets byte-identical
  both regimes) ⊕ (metric regression-guards green + the ERR-067 tripwire) ⊕ (MRO
  sibling-not-child + post_init-once) ⊕ (round-trip inherited from offsets). Deliver
  TWO new files (`test_face_layout_typed_key.py`, `test_facefield_hierarchy.py`); keep
  the commit-1 metric suite + `test_face_layout.py` + `test_face_streaming_normal.py`
  UNMODIFIED as guards.

How to apply: for any "introduce a parent ABC + re-parent siblings + generalize a
stringly-keyed layout to `[K]`" bit-identical carve — (1) find the production
key-index cardinality; if the shared layout is exercised at ONE index, MANUFACTURE a
multi-index synthetic and make the load-bearing mutation the generalized index term
(`n*pos→pos`), demonstrated RED on multi-index / GREEN on single-index; (2) hand-compute
the offset golden so it is carve-independent, pin slice CONTENTS not just shape; (3)
treat per-leaf metrics as regression-guards + add ONE "ABC stays metric-free" tripwire
(the ERR-067/`G_sd=0` guardrail); (4) the MRO gate is sibling-not-child +
`post_init`-fires-once, and the break mode is collapsing an intermediate class name;
(5) round-trip is inherited from the offset gate. `vv §Mode-12 / L11 / L1 config-blindness`.

---

## L21 — Verifying a WRAP resolvent operator that DELEGATES to an already-verified free-function engine (`A_BB.solve` WRAPS `carlson_inward_sweep_from_source`; coupled-block campaign step 1c): the equivalence keystone is bit-id INHERITANCE + a Mode-11 call-COUNT spy; structural independence MUST come from the ENGINE's OWN closed-form (NOT its recurrence); and an operator-internal transform line you CANNOT monkeypatch (a `[:, ::-1]` data reversal) is pinned by a spy on the call ARGS + a graded-mesh non-vacuity check

When an operator's `.solve`/`.apply` is a thin WRAP of a standalone engine
already carrying its own L0/L1 gates (the #282 `carlson_inward_sweep_from_source`
had a full closed-form convergence gate BEFORE `RadialCharacteristicOperator`
wrapped it), the operator is the SUT — its gates add the WRAP layer (space
cells/corner views, source-block reading, the two-leg composition) the engine
gate never touches. Five reusable test-design facts, recurring across every
"WRAP→bit-id" step of a block-operator campaign:

- **The equivalence keystone is bit-id INHERITANCE + a Mode-11 call-COUNT
  sentinel, not a fresh value claim.** Monkeypatch the engine in the OPERATOR's
  module namespace (`_rc_mod.carlson_inward_sweep_from_source`) with a spy that
  DELEGATES to the real engine and records `(args, result)`; assert
  `len(calls) == 2·n_levels` (2 legs/level) AND `solve.values ==`
  (`array_equal`) an independent two-leg reference built with the REAL engine
  (imported at test top, UNPATCHED). The count proves the path is EXECUTED (a
  divergent inlined copy leaves the counter at 0 — Mode-11); the bit-id proves
  the WRAP changed nothing. Necessary-NOT-sufficient (L2) → pair with a
  structurally-independent value anchor.
- **Structural independence is from the ENGINE's OWN gate, not just the SUT
  (ERR-032).** The operator's convergence anchor MUST be the analytic
  closed-form ODE (`φ = q/σ + (φ_R − q/σ)e^{−σ(R−r)}` for the μ=−1 ray — a pure
  exponential, no marching), the SAME reference the engine gate uses — NOT a
  hand-re-execution of the DD recurrence (procedural, not structural). Anchor
  the INWARD leg's cells (its entry face φ_R = corner_in is FREELY settable, so
  it maps 1:1 to the closed form); the outward leg's entry is the derived
  pole_face (harder to anchor cleanly → cover it via the spy, below). R=1.0
  keeps σR≤1.3 so the O(Δr²) regime is reached by nx=64 (R=4.0 gave ratios
  3.16→3.75, still short of the [3.4,4.6] band at nx=64); DISTINCT per-group
  (σ,q,φ_R) + max-over-groups error makes it genuinely ≥2G and catches a
  group-axis view transpose.
- **An operator-INTERNAL transform line you CANNOT monkeypatch (the `[:, ::-1]`
  outward-leg reversal + `[::-1]` un-reversal) is pinned by the spy on the call
  ARGS + a graded-mesh non-vacuity check.** You may not `git checkout` the
  committed operator, and monkeypatching a bare slice inside `solve` is
  impossible — so assert the outward call received `dr[::-1]`, `q_plus[:,::-1]`,
  `σ[:,::-1]` (the WRAP's data-carried-orientation contract, the "2.5a
  discipline"). The `dr`-WIDTH reversal is vv Mode-5 BLIND on a uniform mesh
  (`dr[::-1]==dr`), so the gate MUST run on a GRADED mesh and assert
  `dr[::-1] ≠ dr` (non-vacuity) — else a dropped width-reversal hides. (The
  `q`/`σ` reversals are caught on ANY spatially-varying-data mesh; only the
  width reversal needs graded. Measured dropped-reversal delta: outward cells
  move O(1) — 5.3e-1.)
- **Pole continuation / leg-chaining is pinned by the SAME spy:** outward call's
  `bc_outer_value` == inward call's returned `phi_face_final` (the exit IS the
  entry, R13) + assert the exit face is non-trivial (≠0) so the equality is not
  vacuously satisfied by zeros.
- **The ISOLATED resolvent adjoint (`solve_transpose`)'s consistency partner is
  the EUCLIDEAN inner product (plain dot), NOT the `V_cell` Hilbert metric.**
  The pure ray-block transpose is the Euclidean reverse-mode adjoint; the
  metric `.H` is realized ONCE at the composite (L19). Do NOT reuse the sibling
  `_g_recip` (it carries `G_sd`). Gate `⟨solve u,v⟩=⟨u,solve_transpose v⟩`
  (< 1e-11, het σ, ≥2 draws) + assert the source μ=+1 corner cotangent is
  EXACTLY 0 (the q½ fold never writes that slot). Teeth: a sign flip in the
  engine transpose (`−f_bar → +f_bar`) → defect 5.9e-17 → 1.0.

How to apply: for a WRAP operator over a verified engine, ship (1) the count-spy
+ bit-id-reference keystone; (2) an operator-level convergence gate re-anchored
on the engine's OWN closed form (graded MANDATORY, ≥2G distinct-group,
φ_R≠q/σ), teethed by monkeypatching the engine with DD-coefficient mutants
(closure/denom sign, diamond factor RED both meshes; `dr[k]→dr[k−1]` uniform-
BLIND/graded-RED = the Mode-5 keystone); (3) the spy-on-args reversal +
pole-continuation gates on a GRADED mesh with non-vacuity checks (the escape
hatch for un-monkeypatchable internal transforms); (4) Euclidean (not metric)
adjoint reciprocity + the unused-slot-is-zero pin; (5) constructor negative
gates for the non-carrying CONTROL (cylinder needs a level-structured
quadrature — GL is slab-only; `match=` the specific message, + a positive
build control — L4 net-new teeth). `vv §Mode-11/8 / L1 config-blindness / L2
structural-independence / L4 gate-liveness`.

## L22 — Verifying a WELDED-fold UN-WELD (a per-cell/per-level closure hand-rolled at N sites → ONE single source; `A_BA` the ψ½ Schur fold, Coupled Block Operator Step 2): the single-source routing proof is a Mode-11 wrap-COUNTER on the object the sites construct (EXACT `2·n_levels`, not `>0`), bit-identity INHERITS from the manufactured-anisotropic fold contract, the transpose must single-source too (a hand-rolled `0.5` = `coeff[0]`), and the S/F arms feed ℓ=0 ONLY so the contract MUST manufacture ℓ≥1

Rule: to gate an un-weld that extracts an inlined closure to one source, ship
SEVEN gate types — and design against the fold CONTRACT (shape-agnostic), never
a specific call surface (isolate "how the extracted source is invoked" behind a
one-line `# BIND:` helper the implementer flips; run the correctness gates
THROUGH the real consumer operators, which are shape-stable).

Why: (1) the CENTERPIECE is single-source ROUTING — a Mode-11 wrap-counter on
the SAME object the sites construct, asserting BOTH consumers (S-fwd AND F-fwd)
ENTER it; a green-but-unrouted arm = a divergent inlined copy = Mode 11. Count
EXACT `2·n_levels` (not `>0`). The seed arms `local-import` the fold INSIDE
`apply`, so `monkeypatch.setattr(source_module, fold, spy)` IS picked up. (2)
Bit-identity `np.array_equal(out, old_loop_ref)` INHERITS verification from the
fold contract — but PAIR it with a structurally-INDEPENDENT closed form
(`½·emission`, independent of the fold fn) so it pins the coefficient, not just
the routing. (3) The convenient config nulls the hardest term: the production
arms feed **ℓ=0 ONLY**, so an S/F-only gate is STRUCTURALLY blind to `P_ℓ(±1)`
for ℓ≥1 — MANUFACTURE ℓ=0+ℓ=1 (≥2G distinct-group) and assert
`Q̄(±1)=½Q₀±(3/2)Q₁`; the tooth (drop `sign^ℓ`) reds anisotropic (`3·|Q₁|`) but
the SAME mutation on an iso-only input stays 0.0 — the iso-null asymmetry IS the
necessity proof (L1/L11). (4) A TRANSPOSE hand-coded as a bare constant (the
S-adjoint's `0.5` = `coeff[0]`, bypassing the fold) needs its OWN single-source:
gate the helper-level Euclidean contract `⟨fold(m),y⟩=⟨m,fold_T(y)⟩` (tooth: a
`0.6` ℓ=0 coeff reds) AND the operator-level fwd↔adj CONSISTENCY
`⟨A_BA φ,χ̄⟩=⟨φ,A_BAᵀ χ̄⟩` (tooth: patch the FORWARD fold→0.6 → the adjoint's
still-`0.5` DISAGREES → reds). The consistency gate is BLIND to a SHARED
coefficient (both legs scale together) — it catches the fwd↔adj MISMATCH the
hand-rolled constant risks; the VALUE is pinned by the fwd contract. If the
chosen shape is factory-only (no transpose surface), those two carry the whole
contract — do NOT force a `.apply_transpose` gate; flag the binding.

How to apply: (a) refutation #4 — split configs: the fissile arm (F: `χ·νΣf·φ`)
is IDENTICALLY zero on a non-fissile mixture → a VACUOUS 0==0 gate; build a
fissile mixture + a `max|emission|>1e-6` NON-VACUITY guard before the assertion.
(b) refutation #6 — cylinder+slab (`space is None`) are the non-carrying CONTROL
(emit `None` ray + counter 0), NOT "other geometries". (c) exercise the arm on a
HAND-BUILT seed-carrying composite even when the production DRIVER is "dormant"
(the arm fires on any non-None-ray composite). (d) BIND-isolate the direct
extracted surface behind `# BIND:`; those rows are `xfail(strict=False,reason=)`
→ flip to xpass on landing (L4). `vv §Mode-11/8 / L1 config-blindness (ℓ=0
iso-null) / L2 structural-independence / anti-#3 ≥2G / anti-#4 fissile`. Full
recipe → §3 [aba-schur-fold-unweld-verification].

## L23 — Verifying a SYNTHETIC ACCELERATOR (DSA #2, pre-carve spec): the `correction→0` property PARTITIONS the failure surface, and FP-invariance is STRUCTURALLY BLIND to the majority of it — the object+rate gates are the load-bearing catchers, not supplements

The deepest test-design fact for any accelerator whose correction vanishes
at the fixed point (DSA, and the TSA/#5 sibling): the correction→0 property
makes the accelerator correctness-safe BY CONSTRUCTION for EVERY error in
{R (restriction), P (prolongation), the low-order operator `A_diff`,
correction sign/scale} — those leave the converged FP **identically**
unchanged (not sub-floor; identically, because the correction source is zero
at the FP regardless of how wrong the machinery is). This PARTITIONS the
plausible-error surface into two disjoint classes with disjoint catchers:

- **Value-changing class = corruption of the WITHIN-GROUP transport operator
  itself** (`A=L+C−S−B`): the #215 σ_r-fold (`A_wg.solve(S_residual)` with
  the σ_r-sweep = `Σ_s0·P_iso` only for isotropic flux → 46–56% wrong on
  anisotropic), a sweep sign flip. Caught by **FP-invariance** (converged
  flux ≡ plain-SI FP) — but ONLY on an ANISOTROPIC config (vv Mode 9; the
  isotropic-reflective box is the designed-green degeneracy where the fold
  is exact). This is the class FP-invariance exists for — and it is exactly
  ONE of the eight canonical errors.
- **Rate-only class = the accelerator machinery** {R, P, `A_diff`,
  correction sign, low-order D/removal, boundary rows}: FP-invariance is
  **structurally BLIND** to all of it (correction→0). Caught ONLY by
  **object gates** (R conservation `⟨1,Rr⟩=⟨1,r⟩` with NON-uniform weights +
  the anti-#8 per-cell delta-source distribution pairing; P≡Rᵀ adjoint as
  the full MATRIX not a sampled bilinear form; the low-order matrix ≡
  hand-posed stencil, element-wise, Mode-12-safe object-level) + **rate
  gates** (ρ vs 0.2247c Fourier; the c→1 × optical-thickness iteration-count
  table; reflective-BC divergence — the inconsistent-low-order operator
  diverges on reflective, is forgiving on vacuum).

The battery-completeness corollary: a plan that gates a correction→0
accelerator's correctness ONLY on FP-invariance — even a *correctly*
anisotropic FP-invariance — ships SEVEN of eight canonical errors silently.
FP-invariance and the rate gates have DISJOINT blind spots (ρ can't see the
value, the flux can't see the rate) — which is WHY both exist. Three
gate-design closures recur: (1) the **correction→0 gate's** own blind spot
is a DEAD R (δφ≡0 passes "→0" trivially) → pair with a first-iterate
NON-triviality lower bound. (2) the **no-masking** converse: mutate the
transport OPERATOR (not the machinery) → DSA must converge to the SAME wrong
answer as SI (proves it changes the RATE not the VALUE in both directions).
(3) the **#215 routing guard** is a Mode-11 wrap-counter / AST tripwire (the
within-group solve must NOT route through the σ_r-fold), catching at
design-time what FP-invariance catches numerically. Assemble-consistency
(3a) is object-level matrix diff (interior/boundary/scaling decomposition —
"≡ up to a boundary row" = interior block element-wise-equal AND the
nonzero-diff support ⊆ {0,N−1}); its independence rides a hand-posed anchor,
never the two ORPHEUS assemblies agreeing (L16 one-source proof). ρ estimator
= the RESIDUAL ratio (Signature 9: the increment under-reports as c→1);
1G is LEGITIMATE for the ρ claim (rate is flux-shape-independent — declare
the claim layer, AGENT.md §5). Anisotropic FP configs MUST build the
`Mixture` constructor DIRECTLY (`make_mixture` nulls SigL AND every A/B/C/D
mixture has Sig2=0 — lessons L1; a `make_mixture` ℓ≥1 config silently nulls
the moment DSA must leave untouched).

How to apply: for any correction→0 accelerator, FIRST draw the value/rate
partition table (§0 of the spec) — for each of the 8 canonical errors name
the gate that reds and confirm FP-invariance is blind. The object+rate gates
are LOAD-BEARING (the only catchers for 7/8), not supplements; a
FP-invariance-only plan is Mode-9-complete for the σ_r-fold and Mode-9-EMPTY
for the machinery. Sharpens vv Mode 9 (the config degeneracy) with the
gate-completeness corollary. Full spec → `.claude/plans/dsa_verification_spec.md`.
`vv §Mode-9/12 / anti-#8 / Signature-9 / bit-identity`.

---

## L24 — Verifying a THREE-DOF SEPARATION campaign (operator ∥ splitting ∥ realization): MEASURE the proposed acceptance criterion BEFORE the plan is anchored on it — a signature-invariance AC is permanently GREEN, and the real red is the WELD (one record carrying two stages, its stale half silently re-derived downstream)

A campaign that separates degrees of freedom the code conflates arrives
with an acceptance criterion of the shape *"changing X must not touch Y"*.
That shape is a trap: it is usually **already true by SIGNATURE**, so the
gate is unfalsifiable from the first commit. The highest-value output of a
proactive dispatch is REFUTING it with a measurement (L10), then supplying
the criterion that IS red.

- **Run the AC as a probe before writing any gate.** Measured on `main`
  @ `b0a003b4`: `build_within_group_system(sn_mesh, mat_xs, *,
  scattering_op, scattering_order)` takes **no** strategy parameter ⟹ the
  posed `A` is strategy-invariant by construction; and `A x == (M−ΣN) x`
  reads **EXACTLY 0.0** (2-D Cartesian seedless, both schedules) /
  **3.5e-17** (sphere carrying — the coupled grid re-associates, so
  `array_equal` is the WRONG contract there from day one). Falsifiers
  (drop `B` → 5.15e-2; sign-flip `N` → 8.57e-3) prove the probe
  non-vacuous — it is simply GREEN. → the new **vv Mode-8 SIGNATURE-
  tautological** sub-case: a falsifier check PASSES on such a gate, so
  "does it red?" gives FALSE confidence; the honest gate is on the
  SIGNATURE (adding the knob must red it), and the value row demotes to a
  regression floor.
- **The red is the WELD, and it is findable by object IDENTITY, not
  value.** `WithinGroupSystem` carries stage-2 (`loss`) AND stage-3
  (`implicit_operator`/`explicit_gains`) in one record, so there is no
  boundary at which "strategy may enter" is assertable — and because there
  is none, a SECOND splitting grew beside the first: measured
  `_within_group_si(..., inner_schedule='gauss_seidel')` returns
  `base is record.implicit_operator` → **False** (driver runs
  `ScheduledInvertibleOperator`+`SNMaskedBoundaryOperator`; the record
  still advertises `StreamingCollisionOperator`+`SNBoundaryOperator`),
  while `jacobi` → **True**. The `is`-identity gate across the strategy
  ladder catches this; no value gate can (both splittings reconstruct the
  same `A` — that is what makes them both *valid* and the twin *silent*).
  Ship the green arm as the CONTROL leg so the asymmetry is pinned.
- **`A − M` may not even be SPELLABLE.** Measured:
  `IncompatibleOperatorComposition — OperatorSum requires equal domains;
  got CoupledSpace('coupled(full_field)', shape=(128,)) and
  FullFieldSpace('full_field', shape=(128,))`. Same shape, different NAME
  ⟹ `FunctionSpace.__eq__` refuses. So the splitting identity can only be
  checked by APPLYING to a state with a manual carrier unwrap — a
  structural red worth naming in the plan, because it is invisible to
  every value probe.
- **Monomorphism reds are cheap to measure and belong in the plan, not
  the implementation.** `F.domain is None` on the production
  `from_solver_data` build (every other leaf reports a concrete space);
  `F.H` then BUILDS and silently returns a bare Euclidean transpose
  (`_AdjointOperator.apply` applies the metric only when both inner spaces
  are non-`None`). Also: the carrier guard is NON-UNIFORM — one leaf
  raises a typed `TypeError` with a remediation message, another a raw
  `AttributeError`. Gate the CONTRACT (one carrier accepted, all others
  refused with a typed message), never the dispatch SPELLING when a
  corpus ruling parks it (#261) — the plan must respect the parked
  ordering or it will be rejected on process, not on math.
- **The rate DOF needs its own closed form, and the degenerate regime
  kills the RATE gate too.** ρ(M⁻¹N) is the one constraint structure
  cannot supply. Its anchor is hand algebra on the spatially-flat /
  angularly-anisotropic mode (`L=0` ⟹ `ρ = c` standard SI, `ρ = c/(1−c)`
  for the σ_r fold, diverging at `c ≥ ½`); an ISOTROPIC power-iteration
  seed makes `N ≡ 0` and ρ reads **identically 0** — ship that as a
  permanent CONTROL leg, not a one-off mutation. 1G is LEGITIMATE for a ρ
  claim (declare the claim layer). The measurand already exists and is
  read by nothing: `SourceIteration.contraction_ratios`
  (`numerics/iteration.py:742`) — cross-check the operator-side power
  iteration against it (`rtol=5e-2`); disagreement means the driver is not
  running the advertised splitting (the weld regression, above).
- **A composition-over-fusion carve owes a PERF gate whose catcher is
  STRUCTURAL, not wall-clock.** `A_ij = Rᵢ A Jⱼ` is the exact shape that
  cost 10–20× on slab / ~6× on the suite with every correctness gate green.
  Wall-clock is machine-dependent and CI-flaky (normalise against an
  in-process same-FLOP calibration if used at all); the deterministic
  catcher is a **leaf-kernel call-COUNT that must not scale with `n_cells`**
  (refine nx 20→160, assert the count is unchanged) plus a `tracemalloc`
  peak-per-DOF ceiling. Capture the baseline on the commit BEFORE the
  composition phase — a perf gate baselined after the regression is worthless
  — and hard-block that phase's merge on it.

How to apply: for any N-DOF separation campaign, (1) EXECUTE the proposed
AC as a probe and report GREEN loudly if it is — a plan anchored on a false
red ships unverifiable phases; (2) hunt the WELD (one record spanning two
stages) and gate it by object `is`-identity across the strategy ladder,
with the already-green arm as the control; (3) enumerate the *unspellable*
states (cross-carrier composition, absent pencil type) separately from the
red ones — unspellable is worse than red and needs a phase-ordering
constraint, not a tolerance; (4) order the phases so no gate is unwritable
at its phase's merge (write the red gate BEFORE the fix that flips it);
(5) the rate DOF gets a closed form + a degenerate-seed control; (6) the
composition DOF gets a call-count perf gate baselined pre-carve. Plan of
record → `.claude/plans/campaign_verification_plan.md`. Sharpens vv Mode 8
(signature-tautological) + Mode 9 (rate-blind degeneracy); mirrors L10
(refute the premise) and L16 (one-source proof / object-not-scalar gates).

---

## L26 — A metric-adjoint (`.H`) reciprocity gate is blind exactly when `[G, A] = 0`: the decidable criterion is the COMMUTATOR, not "uniform h" — and a leaf whose action commutes with the metric ALGEBRAICALLY needs a second, metric-agnostic mutation or its row is a dead gate

The received G1.4-style recipe — *"make the metric non-degenerate and
non-uniform (non-uniform `h`, curvilinear so `V_cell` spans an order of
magnitude); ship the uniform-`h` leg as the config-blindness control"* —
is **two-thirds wrong, and both errors are measurable in minutes.** The
mutation is the ERR-067/Mode-12 one: drop the metric inside
`_AdjointOperator.apply` so `.H` degrades to the bare Euclidean
transpose. Its residual is the **commutator** `[G, A]`: with diagonal
`G`, `⟨Ax,y⟩_G = yᵀGAx` and `⟨x,Aᵀy⟩_G = yᵀAGx`, so the mutation is
**identically silent iff `GA = AG`**. Three consequences, all MEASURED
on the SN leaf set (P0, `main` @ `b0a003b4`):

- **"Uniform `h`" is NOT the blind config.** A uniform-`h` slab under
  `gauss_legendre(4)` still REDs at 1.3e-1 (`L`) / 4.0e-1 (`S`) / 2.7e-1
  (`F`) — because `G = V_cell·w_n` and the *quadrature weights* vary even
  when `V_cell` does not. Shipping it as the "blind control" would assert
  a falsehood and the control would RED. The honest blind config is a
  **globally CONSTANT** metric: `G = c·𝟙 ⟹ G⁻¹AᵀG = Aᵀ` exactly. Build it
  deliberately — Cartesian slab (constant `V`) × `gauss_legendre(2)` (the
  only SN rule with equal weights, both exactly 1) × `h = 1/√3` so the
  bulk constant `V·w` EQUALS the trace constant `|Ω·n|·w = 1/√3` across
  the bulk⊕trace composite. MEASURED: one unique metric entry, every leaf
  blind at ≤ 4.3e-16 (vs 1.40/1.01/1.46 on the sphere). ALWAYS pin the
  constancy as a precondition (`np.unique(...).size == 1`) or the leg
  silently stops being the blindness proof — verified reddenable by
  swapping the fixture for the sphere.
- **The CURVILINEAR lever only constrains the leaves the metric fails to
  commute with, and the quadrature FAMILY decides which.** A
  `level_symmetric` rule has **constant weights** (MEASURED: one unique
  value), so on a non-uniform sphere `G = V(space) ⊗ c` — which commutes
  with every SPACE-DIAGONAL, angle-mixing operator. Result: `S` and `F`
  measure **exactly 0.0** under the mutation; only `L` (bulk↔trace
  coupled) reds. Swap to `gauss_legendre(4)` and `L`/`S`/`F` all red at
  O(1). So "curvilinear + non-uniform `h`" is NOT sufficient — the metric
  must vary along **both** the spatial and the angular axis. Check the
  weight vector's unique count before trusting the config.
- **Some leaves are metric-blind by ALGEBRA, in every reachable config —
  say so, and close them with a DIFFERENT mutation.** A
  `MultiplicationOperator` (`C`) is diagonal, so `Cᵀ = C` and `G⁻¹CG = C`
  for ANY diagonal `G`: no configuration exists. A specular boundary
  (`B`) is a signed PERMUTATION mapping `μ→−μ`, preserving both `|Ω·n|`
  and `w_n`, so it commutes with the trace metric; the only
  non-permutation face laws (white/albedo) are **refused by `SNMesh`**
  (MEASURED: `ValueError`, only `reflective`/`vacuum` accepted), so there
  is no reachable config either. Those rows are Mode-10
  *exercised-but-unconstrained*. Do NOT assert their silence (it
  calcifies a fact a legitimate change could move) and do NOT claim teeth
  they lack. Close them with a **second, metric-agnostic mutation** —
  scale the adjoint by 2, which reds EVERY leaf at exactly 0.5 relative —
  so each row still provably pins `apply_transpose`'s structure, and
  point at the gate that DOES pin their metric (for `B`: the existing
  L11 drop-`|Ω·n|`-from-the-REFERENCE control in
  `tests/sn/operators/test_g_adjoint_reciprocity.py`).

**Two companion traps caught in the same file, both worth the habit:**

- **The xfail-strict FALSE POSITIVE bit on first write (L4, live).** The
  G1.5 row "an anonymous leaf must refuse construction" is written as
  *construct → if it succeeded, demonstrate the degradation → `pytest.fail`*.
  The demonstration (`.H` vs `apply_transpose` on a meshless `(ng,1)`
  probe) blew up inside `np.einsum` for `ScatteringOperator` — so that
  row "xfailed" while asserting **nothing** about the claim. Fix: the
  demonstration is BEST-EFFORT, wrapped, and its outcome (including its
  own failure) is reported as **evidence text inside the `pytest.fail`
  message, never as the verdict**. Rule: in a strict-xfail body, there
  must be exactly ONE statement that can fail, and it must be the
  documented one — run `--runxfail` and READ every message before
  crediting the gate.
- **Prove the XPASS flip, don't assume it.** Simulate the landing phase
  with a throwaway pytest **plugin** (`-p mysim` on `PYTHONPATH`) that
  monkeypatches the production surface the phase will change
  (`from_solver_data` raises without a space; `_apply_faces` raises a
  typed `TypeError`) — the strict xfails then turn into
  `[XPASS(strict)] …` hard failures, which is the proof the marker
  forces acknowledgment. Cheap, in-process, and NEVER a `git checkout`
  (this tree carries uncommitted-by-policy skill files with no history).

Also durable from the same pass: **a "reuse the existing fixture" brief is
a hypothesis, not an instruction.** The neighbouring SN operator fixtures
carry `placeholder_materials` (SigS/χ/νΣf all zero) — MEASURED: `F` is
then the ZERO operator and its reciprocity row is the tautology `0 == 0`
(the non-vacuity guard reds with `|⟨Ax,y⟩_G| = 0.0` on all four
geometries). Build the fissile `Mixture` DIRECTLY (L1: `make_mixture`
nulls `SigL` AND `Sig2` and has no P1 channel), and state in the file
WHY each neighbour was rejected — with the measurement.

How to apply: for any `.H`/adjoint/metric gate, (1) write `G` explicitly
and ask *which leaves commute with it* BEFORE choosing the fixture —
that, not "non-uniform mesh", is the config criterion; (2) partition the
leaves into metric-CONSTRAINED and metric-INVARIANT by running the drop
mutation, and give the invariant ones a metric-agnostic non-vacuity
mutation; (3) build the globally-constant-metric leg deliberately and
guard its constancy; (4) `--runxfail` every strict xfail and read the
message; (5) prove the flip with a plugin-based landing simulation.
Sharpens vv Mode 12 (the ERR-067 metric-repair closure) with the
commutator criterion, and vv Mode 10 (the third state:
exercised-but-unconstrained *by algebra*, where no isolating regime can
exist). Gate file → `tests/sn/architecture/test_monomorphic_leaves.py`.

---

## L25 — Building the call-count PERF gate (P0 of the three-DOF campaign, EXECUTED): "the count must not scale with `n_cells`" is TOO COARSE — the invariant is PER-AXIS, the regression it catches is EXACTLY value-identical, and a perf baseline is a (number, FIXTURE) pair

L24 §6 specified the composition-cost gate as "a leaf-kernel call-COUNT
that must not scale with `n_cells` (refine nx 20→160, assert the count is
unchanged)". Building it against the tree corrected that in four ways.
Every number below is MEASURED (branch `refactor/operator-strategy-layers`,
`orpheus/` @ `360a8087` ≡ `361471ec`, host `.venv` Py 3.14.3, `python -O`).

**1. The spec's leaf name was STALE — prove the wrap fires BEFORE asserting.**
`DiamondDifference.update_batch` does not exist (renamed at S6.4(e) to
`cell_kernel_batch` / `residual_kernel_batch`). The genuine apply-direction
leaf on `A.apply(x)` is **`residual_kernel_batch`**. The non-vacuity proof
is two-part and the SECOND part is what makes it a proof: (a) the counted
method fires (>0) on a real apply on every geometry; (b) the plausible
SIBLINGS (`cell_kernel_batch`, per-cell `residual`, per-cell `update`)
stay at **exactly 0** through the SAME harness — which simultaneously
demonstrates the harness reports 0 for an off-path name, i.e. that the
vacuity branch is real (vv Mode 11). Never assert a count before both.

**2. The invariant is PER-AXIS. Measure which axes are hoisted FIRST.**
MEASURED arity of ONE `A.apply(x)`:

| path | count | scales with |
| --- | --- | --- |
| 2-D Cartesian `ScanMarch` (production default) | `8 · ny` = 80/80/80/80 at nx=8/16/32/64 | **ny only** — invariant in nx, in `sn_order` (S2/S4/S6 → 80/80/80 at N=8/24/48) and in ng (1/2/4 → 80/80/80) |
| 1-D slab `CumprodScan` **apply** | `2 · nx` = 40/80/160/320 at nx=20/40/80/160 | **n_cells — already a per-cell Python fold** |
| 1-D slab `CumprodScan` **solve** (`ordinate_scan`) | `2` at every nx | nothing — fully hoisted |

So a blanket "count must not scale with `n_cells`" would have been FALSE
on `main` for 1-D and would have shipped either a permanently-red gate or
(worse) a quietly-weakened one. The honest gate is four legs:
(a) strict `==` on the **scanned** axis (the catcher); (b) calls/cell must
DECAY under isotropic refinement (1.00→0.50→0.25 MEASURED — admits the
honest per-row linearity of the marched axis while still forbidding the
per-cell regime); (c) an **exact-arity regression FLOOR** on the
already-folded path (pin `2`, red on BOTH sides — an improvement must
re-baseline consciously, L7); (d) the hoisted sibling as a **CONTROL leg**
proving the instrument can see invariance, so leg (c)'s linearity is a
property of the code and not of the harness.

**3. The regression class is EXACTLY value-identical — that IS the argument.**
Re-dispatching the batched kernel once per column of its level slice
(the faithful L16 fold) gives `np.array_equal` **True**, `max|diff|`
**0.0e+00** — the DD residual kernel is cell-local, so the split changes
the arity and nothing else. **No value gate, at any tolerance, against any
reference, can ever see this.** State that measurement in the module; it is
what demotes wall-clock to corroboration and promotes the count to catcher.
Corollary for the mutation battery: **check the batch SHAPE the leaf
receives before choosing the fold axis.** A cell-axis fold is INERT in 1-D
(`psi_bar` is already `(4, 2, 1)`; 2-D gets `(3, 2, nx)`) — mutate the
ordinate axis there instead (`{20:40,…}` → `{20:200,…}` → RED). A mutation
that fails to red because it is inapplicable is NOT "the gate has no teeth".

**4. A perf baseline is a (number, FIXTURE) PAIR — own the fixture, gate it.**
Do NOT source a perf fixture from the package's shared `_config`: the
correctness gates are free to retune it (a Mode-9 fix bumping S4→S6, a
third region) and every committed constant silently starts describing a
different problem. Import the shared *conventions* (the direct `Mixture`
builder, the one posing site, the carrier-generic probe state), own the
*sizes*, and ship a **fixture-fingerprint gate** (`n_cells`, `n_ordinates`,
`n_groups`, `n_regions`, `n_dof`) whose failure message says RE-CAPTURE.
Read `n_regions` off `material_xs_field().cells_by_material` (typed,
dimension-agnostic, and it counts materials as the OPERATORS resolve them)
— never off `mesh.mat_map`, which is 2-D-only and `Optional`.

**Tolerance discipline, split by contention-immunity.** MEASURED: P-2
(wall clock, ratio to an in-process dense contraction of the same nominal
FLOP count) = 22.45–24.41 across runs → the constant is PROVISIONAL on a
contended host; set it at ~3× and say so. P-3 (`tracemalloc` peak /
`n_dof·8`) = 3.0299–3.0328, spread **0.03 %** → contention-IMMUNE, so
TIGHTEN it rather than padding it: MEASURED `+1/+2/+3` held full-field
temporaries → `4.035 / 5.033 / 6.033`, so `4.0` catches the FIRST extra
field where the plan's suggested `6.0` needs three. Note a PEAK metric is
blind to transients below the existing peak (a 0.6 MB per-call scratch
moved the reading by 0.0006); the mutation that proves its teeth must
HOLD the temporary, which is also the mechanism the gate exists to catch.

How to apply: for any composition-over-fusion carve, (1) find the real
batched leaf by wrapping candidates and reading the counts — never trust a
spec's illustrative method name; (2) tabulate the arity against EVERY axis
(cell, ordinate, group, and each spatial axis separately) before writing a
single assertion; (3) prove the fold mutation is value-identical and put
that number in the docstring; (4) split the legs by contention-immunity —
pad the wall clock, tighten the allocation; (5) fingerprint the fixture.
Delivered: `tests/sn/architecture/test_composition_cost.py` (9 gates,
0.9 s, pyright-clean). Refines L24 §6; sharpens L16 (the perf-regression
precedent) and L4 (prove every gate's teeth).

---

## L28 — REPAIRING blind gates: re-pose onto a regime-INDEPENDENT mechanism, and measure "before" with `git show`, never `git checkout`

Context: boundary-machinery review phase B0.3 (2026-07-30), seven measured
coverage defects in a suite that was already strong (30/31 mutations caught).
Repairing blind gates is a different design problem from writing them, and the
four rules below are what the job actually turned on. Result: **31/31** on the
auditor's own set, plus the historical miss closed.

**1. A DECAYED `catches` marker is repaired by moving to a REGIME-INDEPENDENT
mechanism invariant — NOT by driving the fixture back into the regime.** The
reflex is "the config converges in 6 outers but the bug needs 30–60, so tighten
the tolerance". MEASURED: that is impossible here — with ERR-052 re-introduced
the un-normalised iterate stabilises at `|φ|max = 5.6552e-01` at 6, 24 AND 32
outers, because `F·φ/k` is scale-neutral once `k` converges, so the catalogued
denormal catastrophe is **unreachable in its own fixture at any depth**. The
repair is to assert the *output convention the fix establishes* — here
`P(φ) = ∫νΣ_f φ dV = 1`, which holds at EVERY outer count, so the marker cannot
decay again. Mutation: 1.0 → 0.0739 (13.5×) on the subcritical leg, 1.0 →
0.84375 on the supercritical one; `rtol=1e-12`, no margin. **Always check
whether the catalogued mechanism is reachable before trying to reach it** — and
record the answer, because the next reader will otherwise chase it.
**Corollary:** the mechanism reference must NOT be computed by the routine that
ESTABLISHES it (`compute_production_rate` is what the normaliser divides by) —
hand-assemble it from the mixture + mesh + solution, or the gate is a tautology.

**2. A Mode-12-blind "hand-computed" test needs BOTH halves: break the
invariance AND source the reference independently.** Breaking the invariance
alone (non-constant input) leaves you re-transcribing the production weight
formula — which is what the *existing* catcher already did (vv L11 procedural,
not structural). Type the PUBLISHED table (Abramowitz & Stegun GL nodes/weights)
as literals, assert the quadrature matches it as a precondition, and build the
expected number from the literals. Bonus available cheaply: the closed-form
continuum anchor (`∫|Ω·n̂|μ dΩ / ∫|Ω·n̂| dΩ = 2/3`) at an honestly-stated
discretization tolerance — the audit had flagged its absence.

**3. A tautological-`raises` file is repaired guard-by-guard with the guard
named at `file:line` in each docstring, and the acceptance measurement is the
PER-GUARD red table.** 0/14 → 12/14. Every residual miss must be
*categorically* out of scope and said so (the 2 here are composition-walker
mutations, not error guards) — an unexplained residual is indistinguishable
from an incomplete repair. Two production defects fell out of the exercise that
no amount of reading had surfaced: an error class with **zero** raise sites
(its trigger is a no-op ABC default nobody overrides) and a one-site `law=` tag
convention drift. **Pointing tests at production entry points is itself an
audit of production.** When a drift is found and you may not fix production,
assert the *property* (case-insensitively) rather than the drifted *value* —
pins the contract without calcifying the bug.

**4. Measure "before" with `git show HEAD:<file> > <tmp copy beside it>`.** The
tree carried uncommitted work in `.claude/` and `orpheus/`, so `git checkout` /
`restore` / `stash` were forbidden. Writing the pre-repair content to a
throwaway sibling file runs it under the SAME mutation, in the SAME collection,
against the SAME fixtures — a true A/B, then `rm`. This is how every
before/after row in the report was obtained (phantom re-confirmed as `11
passed`; sentinels re-confirmed as `3 SKIPPED`; registry probe RED).

**5. HYGIENE — for the "committed" side of any doc/generated-artifact diff, use
`git show HEAD:<path>`, never the working-tree file.** A concurrent Sphinx
build regenerated `docs/theory/verification/matrix.rst` mid-session (its
`builder-inited` hook), between two of my edits. Diffing a regeneration against
the already-regenerated working-tree copy manufactured a confident,
circumstantial, and entirely WRONG "the committed V&V matrix has always
disagreed with the code" finding, complete with a git-log narrative. Caught
only by checking `git show HEAD:` and the file mtimes. **On a shared tree,
`ls -la`/mtime is part of the measurement.** Retract loudly when it happens —
a wrong finding costs more than no finding.

**6. Never let a repair reduce the catch rate — prove it, don't assume it.**
Re-run the auditor's own harness rather than a re-implementation
(`$CLAUDE_JOB_DIR/tmp/mutplugin{,2,3,4}.py` survived the session), so the
before/after is apples-to-apples. The per-mutation red COUNTS are the real
signal: `white_nocos` 11→13 and exactly +1 on each of the twelve error guards
proved the repairs landed where designed, and every other count being unchanged
proved nothing was disturbed.

---

## L29 — Verifying a DOMAIN NARROWING (an operator's domain shrinks from a superset to the honest half; SN boundary B3): the old gate does NOT go tautological, it stops RUNNING — and the narrowing's own new index remap is activated only in a regime the existing gate does not use

The campaign shape: production hands a realized law the WHOLE slot and discards
the unwanted rows on the way out (`full = law.apply(face_in); out[sel] =
full[sel]`); the honest physics is `Γ₊ → Γ₋`, so the carve narrows the law's
domain AND codomain. A sibling arm (diffusion) already does it right, so the
carve is a domain narrowing against a worked in-repo precedent. Five reusable
test-design facts, all MEASURED (boundary review §B3, 2026-07-31):

- **The "will this gate go tautological?" question usually has a THIRD answer:
  the gate BREAKS.** The reflex worry (`vv` Mode-8 signature-tautological) is
  that "the discarded rows are zero" becomes unfalsifiable once the codomain no
  longer contains them. MEASURED: it stays reddenable (a pass-through mutation
  reds it with its original message) — but the gate never gets there, because
  its *reference expression* (`bc.apply(whole_slot)`) feeds the narrowed
  operator the wrong shape. **Always run the gate under a full-fidelity
  simulation of the carve** (patch the REALIZER, not only the call site) — a
  call-site-only simulation leaves the reference working and hides the break.
  Diagnostic value: `realizer` mode reds 16 items, `narrow` mode reds 0.
- **The gate's TEETH change even when the assertion survives.** After the
  narrowing the original bug (the law's discarded image leaking) is
  *unspellable* — strictly better than red (L24) — and what remains reddenable
  is only the WRITE-TARGET family (wrong target / extra write). So the honest
  re-pose is on the mechanism: `got[complement(codomain_rows)] == 0` plus a
  SHAPE gate on the realized operator (`apply` maps `(n_dom,…)→(n_cod,…)`) plus
  a `pytest.raises` NEGATIVE for the domain guard. Ship the shape gates as the
  post-carve contract: they are RED pre-carve and GREEN post-carve, which is
  the flip-proof the carve is measuring the right thing.
- **A three-way index split hides in a two-way vocabulary.** The face slot is
  `inflow ⊔ outflow ⊔ TANGENT` (`|Ω·n| ≤ ε` belongs to neither) — 4 of 8
  ordinates on the product-quadrature fixture. The shipped gate asserted
  `got[outflow] == 0` and was MEASURED BLIND to a tangent leak that the
  widened `got[complement(inflow)] == 0` catches. Whenever a gate says "the
  other rows are zero", check whether "other" is really the complement.
- **The narrowing can introduce a NEW index remap whose activating regime the
  existing gate does not use — and the SAME remap appears at ≥2 sites whose
  discriminating fixtures are COMPLEMENTARY.** `out[sel] = full[sel]` needs no
  remap; after narrowing, the image is indexed by position WITHIN the codomain,
  so both (a) the row subset of a schedule split and (b) the law's own table
  need `searchsorted(codomain_rows, ·)`. MEASURED at landing, one mutation per
  site over the whole boundary suite: the **split** remap's naive
  `arange(sel.size)` is EXACTLY correct in 1-D and reds exactly ONE gate (the
  new 2-D one) — the sphere-based split gate stays green; the **law's** remap
  is wrong on the SLAB (the mirror reverses order: `perm[inflow] = [3,2]` →
  local `[1,0]`) and reds 16 items, but the CYLINDER and the 2-D
  level-symmetric fixtures stay green (there `perm[inflow]` is already the
  codomain set in order). So 1-D covers one site and 2-D the other, and
  **neither fixture covers both**. Ship both gates, each with an ACTIVATION
  guard asserting its fixture still discriminates (`local != arange`), or one
  silently decays into the other's blind spot (Mode-8 class 7). **Before the
  carve, enumerate the arithmetic the narrowing ADDS, then find a
  discriminating fixture PER SITE** (§0.6, one level in from config-blindness:
  the blind axis is the *dimension*, and it points opposite ways at the two
  sites).
- **The re-posed WIRING gate is structurally blind to the law's MATH — say so
  and name the other catcher.** After the re-pose, the composite gate's
  reference is `law.apply(γ₊ψ)`, i.e. the SUT's own law object; it therefore
  pins which domain / which rows and takes the law's table as given. MEASURED:
  collapsing the law's permutation to the identity, or rolling it by one,
  reddens NOTHING in the re-posed composite file while changing production's
  output hash. The catcher is the bit-identity gate whose reference is the
  RETIRED expression materialised in numpy from the law DESCRIPTOR (there:
  6 and 10 items red). Two gate families, two threat models — record the split
  in the file's docstring or a future reader reads the wiring gate's green
  under a math mutation as coverage.
- **Bit-identity here is `array_equal`, and a SECOND, independently-written
  narrowed implementation is the strongest evidence.** The narrowing removes
  rows from a gather — no reduction tree changes — so any nulp tolerance would
  hide the bug class. Prove it twice: once with a wrapper
  (`full.apply(embed(γ₊))[cod]`) and once with a genuinely reduced kernel (a
  reduced permutation). sha256-identical on both, forward AND transpose, over
  every reachable fixture. Also measure the laws that DON'T narrow neutrally
  (SN's albedo/periodic realize as full-face identities → their narrowed action
  is a DECISION, not a refactor) and ship those as `xfail(strict=True)` so the
  decision cannot be forgotten when a later phase admits them.
- **The "which laws are still un-narrowed?" gate needs an anti-Mode-12 leg, and
  the answer is usually BIGGER than the plan says.** A shape assertion
  (`apply(|dom|) → |cod|`) CANNOT distinguish `Γ₊→Γ₋` from `Γ₊→Γ₊` whenever
  `|Γ₊| == |Γ₋|` — MEASURED true for EVERY quadrature × face pair in the tree,
  so the error class sits inside the shape functional's invariance group. The
  discriminator that escapes it: feed the law a FULL-face probe and require it
  NOT to emit `N` rows (a narrowed law structurally cannot; an endomorphism
  always does). That single leg is what caught the laws the plan had NOT listed
  — the brief predicted "3 un-narrowed laws"; measurement found **6 rows across
  4 kinds**, because the α=0 / α=1 *fast paths* (a bare `ZeroOperator` /
  `IdentityOperator`) are endomorphisms too, and a rank-0 source law raised on
  the narrowed domain. Enumerate from the REGISTRY and add a
  completeness test (`covered == set(registry)`) so a law added later cannot
  escape the gate.
- **A narrowed law cannot compose with an un-narrowed sibling — budget for the
  mixed-composition gates going xfail, and RE-POSE their probe first.** Every
  `α·narrowed + β·full_face` gate breaks the moment the first law narrows.
  Marking them `xfail(strict=True)` is right, but a probe left on the OLD
  domain fails at the *reference* line with a broadcast error, not at the
  documented composition refusal — a misattributed xfail (Mode-8 class 4) that
  will XPASS for the wrong reason. Move the probe onto the narrowed domain
  FIRST, then `--runxfail` and read the message.
- **Prove the strict-xfail set FLIPS, and watch for a Mode-8 class-4 trap
  inside that very proof.** A landing-simulation plugin that narrows ONE law
  must turn exactly that law's rows `XPASS(strict)` and leave the others
  `XFAIL`. First attempt simulated the α=0 arm as `ScaledOperator(0.0, …)`; the
  constructor REFUSES a zero scalar ("degenerate; use ZeroOperator
  explicitly"), the xfail swallowed the `ValueError`, and the row read as "did
  not flip". The simulation must reproduce the shape production would actually
  emit — otherwise the class-4 defence is itself class-4.
- **Method note that paid for itself twice.** (a) The reduced operator does NOT
  validate its domain — `PermutationOperator(reduced_perm).apply(full_slot)`
  silently truncates and returns wrong numbers (MEASURED), so "hand it the full
  slot" stays spellable unless the carve adds a guard; the negative test for
  that guard is a deliverable, not a nicety. (b) A first tangent-leak mutation
  was VACUOUS (a permutation's image on the tangent rows is structurally zero
  once the domain is `Γ₊`) — caught by the sha256 positive-control table before
  any colour was read. Keep a per-variant hash table; it is cheap and it is the
  only thing standing between you and a confident false "no gate catches this".

---

## L30 — MIGRATING the inherited test surface of a domain narrowing (B3.4a, the phase AFTER L29's carve): a rectangular operator's self-adjointness gate is not lost, it is RECIPROCITY against a sibling the tree already builds; and the guard the narrowing dissolves is replaced by a guard on the CLASSIFIER, whose discriminating fixture is usually ONE quadrature

L29 is the carve; this is the test-migration phase that follows it — the 28
inherited reds of the white/prescribed-inflow narrowing (`AngularAverageOperator`
+ `IncomingSourceOperator` typed `Γ₊ → Γ₋`). Every finding below is MEASURED
(2026-08-01, `scratch/b3_4a_test_migration.md`, harness `scratch/b3_4a_mutations.py`).

- **A "cannot be posed on the narrowed operator" gate is usually a RECIPROCITY
  gate in disguise — look for the mirror object the tree already builds before
  reaching for `xfail`.** The brief pre-authorized a strict xfail for
  `⟨Ax,y⟩_w == ⟨x,Ay⟩_w` on a now-rectangular `A`. It was NOT needed. Read WHY
  the old endomorphic version held: `A = 1_N ⊗ (cos_w/norm)` with the SAME
  masked `cos_w` on both the broadcast and the contraction, so `W·A` was
  symmetric. Narrowed, the weighted adjoint `A* = W₊⁻¹AᵀW₋ = 1_{Γ₊} ⊗
  (cos_w₋/norm₊)` is *itself a Lambertian average*, i.e. **the SAME production
  constructor at the OPPOSITE face** (`from_quadrature(q, axis, -sign)`) —
  because that face's `Γ₊` IS this face's `Γ₋`. So the gate becomes
  `⟨Aψ,φ⟩_{W₋} == ⟨ψ,A*φ⟩_{W₊}` with BOTH sides production objects and nothing
  hand-built. It rides on `Σ_{Γ₊} w|Ω·n| == Σ_{Γ₋} w|Ω·n|`, which is
  BIT-IDENTICAL on every quadrature × face in the tree — assert it INLINE as an
  activation guard (a future asymmetric rule must red HERE, loudly, not be
  absorbed into a tolerance). Teeth are REAL and distinct from conservation: a
  face-sign-dependent `norm` mutation reds both reciprocity rows while the
  `+1`-face conservation rows stay green.
- **When a narrowing retires a private CLASSIFIER, the replacement gate is on
  the classifier, and exactly ONE fixture in the tree discriminates it.** The
  operator used to carry its own `(outward_sign*mu) > 0.0` outflow test, a twin
  of the trace space's `> TANGENTIAL_EPS`. MEASURED across 8 quadrature × face
  pairs, the two sets agree on SEVEN — `product(n_mu=2, n_phi=4)` on **xmax**
  is the only disagreement (|Γ₊|=2 vs strict-set 4; 4 tangential ordinates),
  and the SAME quadrature on **ymin** agrees. So the gate must name that exact
  fixture AND carry an activation guard (`assert not array_equal(gamma_plus,
  retired_strict)`), or it silently decays into a restatement of its siblings.
  Mutation W1 (restore `> 0.0`) reds that ONE test and nothing else in 110.
- **`|Γ₊| == |Γ₋|` on EVERY reachable face makes the codomain claim
  Mode-12 blind — and part of it is IRREDUCIBLE.** Sizing the codomain from
  `Γ₊` inside `from_quadrature` reddened **0 of 110** tests and CANNOT be
  caught: no quadrature × face pair in the tree has unequal half-traces. The
  *operator-level* twin (an `apply` that broadcasts over the INPUT's leading
  axis) IS catchable — but only via a HAND-CONSTRUCTED operator with unequal
  sizes (`n_outflow=3, n_inflow=5`), never from `from_quadrature`. So: for any
  narrowed operator, ship one hand-built unequal-size row as the anti-Mode-12
  leg, and report the constructor-level half as an honest residual blind spot
  (closing it is a TYPING job, not a test job).
- **A `catches(ERR-NNN)` on a re-posed leg silently decays through the SAME
  Mode-12 hole — re-run the mutation, don't reason about it.** The ERR-047
  delivered-`q` leg, re-posed to "the output has exactly |Γ₋| rows", was fed a
  `Γ₊`-sized probe and stayed GREEN under the exact catalogued bug (`q` sized
  from its input). Adding a **FULL-FACE probe** leg (N rows in ⟹ |Γ₋| rows out)
  restored the red. The `catches` marker was a phantom for one edit cycle and
  only the mutation run showed it.
- **A guard the narrowing NEWLY reaches ships with no negative test — write
  it, and write the positive control with it.** Three guards were net-new at
  this phase (white's `_outflow_restriction`, white's orientation cross-check,
  prescribed-inflow's `_outflow_restriction`); each got a `pytest.raises` with
  a SPECIFIC `match=` plus, for the orientation guard, a four-face positive
  control — without which a realize that raised for ANY reason satisfies the
  negatives (L4's attributability rule).
- **The orientation cross-check FINDS the mismatches the old silence hid —
  expect them in the FIXTURES, and check the production tag-parser too.**
  MEASURED: a sibling file built `face="xmin"` while its law declared
  `axis="x", outward_sign=+1` (= `xmax`) — pre-narrowing it averaged the
  installation face's INFLOW hemisphere and reported nothing. And
  `_law_from_tag` builds every parameter-free law as `law_cls()`, so
  `WhiteBoundary` always gets the default `xmax` orientation whatever face is
  being resolved (reflective correctly threads `AXIS_NAMES[label.axis_index]`)
  — LATENT only because the method's admission registry does not yet list
  `white`. When a law carries its own copy of a datum the installation site
  also names, grep the tag-parser for the parameter-free fallback.
- **Method: the migration's own regression proof is the DELTA, not the green.**
  Run the wide suite before and after: `40 failed` → `12 failed` on
  `tests/{geometry,sn/operators}`, and `12 failed, 4225 passed` on
  `tests/{sn,geometry,numerics,transport} -m "not slow"` — the same 12
  out-of-scope rows, so nothing else moved. Pair it with `-rs` (no
  skip-swallowed sentinel: 110 passed, 0 skipped, 0 xfailed).

---

## L31 — Verifying an `A ≡ B` theorem that holds BY SHARED BODY (B3.4b): the equivalence gate is designed-GREEN under a body bug, and two sibling routes through one assembly have DISJOINT discriminating fixtures

B3.4b re-posed `AlbedoBoundary`'s specular pairing into `R` and routed the two
completions through the bodies that already realize `ReflectiveBoundary` /
`WhiteBoundary`, making `albedo(α, SpecularReturn(a)) ≡ reflective(a, α)` and
`albedo(α, IsotropicReturn(a,s)) ≡ white(a,s,α)` **theorems by construction**.
Every finding MEASURED 2026-08-01 (`scratch/b34b_verification_plan.md`; new gate
`tests/geometry/test_reemission_closure.py`; harness selected by `ORPHEUS_B34B`).

- **The design's own justification is the reason the gate cannot verify.** "The
  two routes execute the same code, so the theorem holds" ⟹ a shared body
  carries a SHARED bug and moves both sides identically. MEASURED: writing the
  `to_local` remap as a bare `arange` left `TestEquivalenceTheorems` at
  **60 passed / 0 failed** while the operator was wrong on 3 of 5 quadratures.
  The `≡` gate catches ARGUMENT drift only (wrong axis / dropped α / wrong
  `law_key`); the catcher is an **independent-expression anchor** written from
  raw quadrature data (`reflection_index`, `weights × omega_dot_n`) importing
  nothing from the realizer above the `np.take`/`np.tensordot` line. Ship BOTH
  and say in the closeout which one is evidence.
- **⭐ Two sibling routes through one assembly need COMPLEMENTARY fixtures —
  neither list covers the other, and a "representative quadrature" serves
  neither.** MEASURED: the specular remap mutation (`to_local`→`arange`)
  DISCRIMINATES on `gauss_legendre(4/8)` + `lebedev(17)` and is BLIND on
  `product(2,4)` + `level_symmetric(6)` (there `perm[sorted(inflow)] ==
  sorted(outflow)`, so the local permutation IS `arange`); the diffuse
  classification mutation (`TANGENTIAL_EPS`→`> 0.0`) DISCRIMINATES on
  **`product(2,4)` ONLY** (mis-admits 2 ordinates) and is blind on all three
  others. Parametrise both routes over the FULL list and keep the blind rows as
  documented controls with an asserted fixture-table test (`if _id in
  discriminating: assert not array_equal(local, arange)`) — so the day a
  quadrature change makes the anchor blind, THAT reds instead of the anchor
  silently becoming theatre.
- **Folding N routes into one assembly EXPOSES latent crashes that were never
  gated — the fold's own gate is net-new, not migrated.** MEASURED:
  `ScaledOperator` refuses a zero scalar, and the pre-fold reflective/white arms
  had only an `α == 1` fast path ⟹ `ReflectiveBoundary(axis, 0.0)` and
  `WhiteBoundary(..., 0.0)` were LEGAL laws (α=0 passes every invariant) that
  died in the numerics layer. Gate the fixed behaviour AND pin the asymmetry
  with a **negative control** (`ScaledOperator(0.0, …)` still raises), else the
  gate would also pass if someone had widened the numerics primitive instead.
- **A refusal added to a law with per-amplitude fast paths MUST be parametrised
  over α.** MEASURED: guarding the refusal with `and law.albedo not in (0.0,
  1.0)` reds ONLY the α∈{0,1} rows — the exact hole that kept albedo un-narrowed
  at both fast paths through two campaign phases. Same shape as L29's
  "narrowing is per-branch".
- **A `pytest.raises` on a refusal is teeth-less without MESSAGE legs.**
  MEASURED: replacing the specific message with a generic one keeps
  `exc.value.law == "albedo"` TRUE, so a bare `raises(BoundaryError)` + `law`
  check stays green. Assert the substrings that name the COMPLETIONS and the
  DEFECT; pair with a positive control (the completed law realizes), else an arm
  that refused everything also passes (`vv` Mode-8 class 5).
- **An "exactly one of X, Y is non-trivial" invariant: pose it over the whole
  registry; carry a violator with a LANDING PHASE as `xfail(strict=True)` and a
  PERMANENT one as a named positive assertion.** MEASURED two violators, not the
  one the design named: `ReflectiveBoundary(α<1)` = 2 non-trivial (→ B5, xfail —
  the XPASS forces the marker's deletion) and `AlbedoBoundary(1.0)` bare = ZERO
  non-trivial (permanent — the degenerate factorisation a SCALAR trace forces;
  an xfail with no unlock is a red that never flips, so assert it positively and
  gate the union with an inventory-completeness leg). NEVER scope such a gate to
  "tag-reachable laws": measured, `_law_from_tag` hard-codes `albedo=1.0`, so
  that scoping silently drops the exact row the design wanted flagged.
- **When a method REFUSES a spelling another method realizes, the second arm
  needs an INERTNESS gate.** The sibling claim "on a scalar trace the closure has
  nothing to fix" is a claim about the OTHER arm and nothing pins it: gate
  `diffusion.realize(law(α, C)) ≡ diffusion.realize(law(α))` ∀C. MEASURED to red
  under a `type(law.response_kernel)` dispatch — the one mutation that would
  otherwise silently change every `BC("albedo", …)` answer.
- **Method traps.** (1) A mutation can red a gate by RAISING rather than
  comparing (a shape guard fires first) — that red is Mode-8-class-4-shaped
  *inside the mutation*; name the attributable catcher separately. (2) NEVER
  baseline against a tree being written: the same command read `32 failed`
  mid-edit and `9 failed` once settled. (3) A snapshot GENERATOR silently rots
  one phase behind its harness — measured, `_generate_bc_equivalence_snapshots`
  raises on 6 of 8 cases while the harness's own failure message instructs the
  reader to run it.
- **⭐ The `≡` must ALSO hold where production reads the FACTORS, and THERE the
  gate has teeth.** A realized-output `≡` gate is blind to a shared-body bug;
  but the same two laws are also interrogated by predicates that read their
  FACTORS, and those two interrogations share NOTHING. MEASURED: with a
  specular pairing legal in `R`, four SN sites spelling
  `law.geometry_map.permutes_ordinates` inline returned different answers for
  two laws that realize to the SAME matrix (sweep-schedule `B_lower`/`B_upper`
  split, DSA admission, ruled-corner, `_reflect_corner`) — a wrong fixed point,
  not a slow one. Gate it as the EQUIVALENCE (not a hard-coded expected table,
  which drifts with the design) + a NON-VACUITY leg pinning both answers
  (`False`-for-everything is exactly the bug's shape) + a leg asserting the
  RETIRED single-tier read DIVERGES, so a "simplification" back to it reds.
- **Constructing a break-exactly-ONE-invariant mutant is a design problem, and
  the obvious candidate usually breaks several.** For a specular-pairing table:
  `np.roll(arange(N),1)` breaks measure AND sign AND involution, so the
  earliest check fires and the row proves nothing about independence. And **no
  ODD cycle can ever isolate the involution**: `a→b→c→a` forces
  `sign(b)=−sign(a)`, `sign(c)=sign(a)`, `sign(a)=−sign(a)`. The working
  construction is `π ∘ σ` (true mirror ∘ a swap of two SAME-sign SAME-measure
  ordinates) — and it needs a DEGENERATE measure class, which
  `gauss_legendre(8)` does not have. Carry the QUADRATURE per row, not shared.
- **A refusal that stands on two different arguments must have its message
  assertions KEYED to the argument.** `AlbedoBoundary(α)` bare is refused at
  every α, but for α≠0 the defect is a pairing-by-array-position and for α=0 it
  is a `VacuumInflow` twin (nothing returns ⟹ no pairing is missing). A blanket
  "the message names both completions and the under-determination" pins a FALSE
  reason on the α=0 row. Key each leg, and give the negative ("α=0 must NOT say
  under-determined") a CONTROL that α≠0 does — else it is vacuously true.
- **A generator can rot one phase behind its harness, and "fix the generator"
  can be a SCHEMA question wearing a bug's clothes.** MEASURED: the BC-snapshot
  generator raised on 6 of 8 cases (still `SNMethodSpace.minimal` for narrowed
  laws) while the harness's own failure message told the reader to run it.
  Giving it face-ful spaces only moves the failure to the SCHEMA — a narrowed
  law cannot emit the FULL-FACE `psi_in` the frozen artefacts store, and
  migrating the schema DESTROYS the property that makes those artefacts
  valuable (frozen pre-narrowing ⟹ old-vs-new is independent). STOP and hand
  the decision back; fix the false instruction in the error path meanwhile.
- **Snapshot re-pose when a completion supersedes a retired spelling: inherit a
  FROZEN artefact generated by the SIBLING law, do not regenerate.** MEASURED
  `0.5 × specular_x_lebedev17["psi_in"][inflow]` is bit-identical (maxdiff 0.0)
  to the albedo-with-specular-closure image, so the `≡` theorem gets pinned
  against an artefact predating every line under test — criterion 2 of the
  re-baseline ladder satisfied by construction, and no self-generated baseline.
  The retired law's OWN artefact is unusable (its Γ₋ rows read ordinates the
  narrowed input never sees) ⇒ delete artefact + generator case together.
  (SUPERSEDED for this harness by **L32** — the user ruled the whole file over
  to derived references. The inherit-from-a-sibling recipe stays correct
  whenever a derived expression is NOT available.)

---

## L32 — INVERT a snapshot generator from RECORDER to REFERENCE: the artefact stops being production's echo and becomes the math; the frozen FILE is what stops generator+production drifting together

A snapshot generator that calls production and freezes its output is
**self-referential** — the gate says `production == a recording of production`,
which detects change and certifies nothing, and it BREAKS on every production
signature change (this one broke twice: B3.2, then B3.4a). Restricting a
pre-change artefact rescues *procedural* independence only (`vv` L11): it
certifies "the new path agrees with the old path", worth exactly what the old
path was right — and a narrowing campaign's whole premise is that the old path
was wrong. The repair is to **compute the reference from the law's equation and
freeze THAT**; the generator then imports nothing from the realization layer and
cannot be broken by it again.

- **The precondition on the recipe: the expression must be TOTAL.** A recording
  pins every bit of the output; a derived reference pins only what the
  expression determines. Check per case that the expression fixes the whole
  array (all 7 did: `0`, a gather, a broadcast average, their sum, a partner
  gather). If an expression is partial ("the image is isotropic", magnitude
  unfixed), inverting TRADES a total drift lock for a partial correctness
  claim — then keep the recording, or pair the two.
- **Derive from the EQUATION, never by transcribing the implementation.** The
  specular pairing came from the isometry `Ω ↦ Ω − 2(Ω·n̂)n̂`, realized by
  matching flipped direction cosines back onto the ordinate table with EXACT
  float equality (MEASURED exact + involutive on gl4/gl8/lebedev17/ls4/ls6/
  product24) — never `quadrature.reflection_index`, which is production's own
  table. Demanding an exact unique hit is what keeps it independent: a tolerant
  nearest-match would agree with a WRONG table. The Lambertian came from
  `J⁻ = αJ⁺` + isotropy (normalising over Γ₋), where production normalises over
  Γ₊; assert at generation time that `Σ_Γ₊ w|Ω·n̂| == Σ_Γ₋ w|Ω·n̂|` within
  reduction-order ULP so the choice cannot silently matter — and if it ever
  fails, the derived one is the correct one.
- **⭐ Once the reference is computable, the FROZEN FILE is the only thing left
  standing between a wrong expression and a green gate.** If the harness
  recomputed it, generator and production would drift TOGETHER through one
  shared expression — a different self-reference. So the harness reads
  `psi_in` off disk and imports ONLY the case registry; make that
  **structural, not documentary**: an AST gate asserting (a) the generator
  imports nothing matching `realizer|method_space|orpheus.sn`, (b) the harness
  pulls only `{CASES, case_by_id, BCEquivalenceCase}`, (c) artefacts on disk ==
  registered cases (no orphan). All three proven-red by feeding the predicates
  a mutated source copy — never by editing a tracked file.
- **Fold the duplicated case identity into the registry while you are there.**
  Pre-migration the law + quadrature of every case were spelled TWICE
  (generator and harness) with nothing keeping them in step. The clean field is
  `compose: Callable[[Realize], Operator]` — a function OF the harness's
  realization step, so `lambda realize: realize(Vacuum())` covers rank-1 and
  `0.3*realize(a) + 0.7*realize(b)` covers rank-N, and the generator still
  imports no realization code. NOT a stringly-typed `("WhiteBoundary", {...})`.
- **Re-derive every tolerance from structure; the migration MEASURES its own
  quality.** Gathers and α-folds are reduction-depth 0 ⟹ `array_equal` (a
  tolerance there would admit the bug); an `n`-term positive-summand contraction
  vs a `tensordot` is `κ=1`, so `|Γ₊|·ε` (the probe being `U(0,2)`, hence
  non-negative, is what makes κ=1 — say so). Retire inherited `nulp=4`/`nulp=64`
  folklore. Then compare the DERIVED reference against the retired recording:
  MEASURED **bit-identical on 3 of 7**, **1–2 ULP on 3**, **98 % on 1** — the
  outlier being exactly the case whose law production gets wrong. That table is
  the migration's own evidence that six claims did not move, only their source.
- **Keep the PROBE bit-identical across the migration** (draw full-face with the
  same seed, then restrict to Γ₊). Free traceability: the stored rows are a
  literal subset of what the retired artefacts probed with, so the change is
  provably schema+reference and never the input.
- **A two-face law cannot be probed with a one-face draw.** Periodic's reference
  is the PARTNER face's outflow (`Γ₋(xmin) == Γ₊(xmax)` as ordinate sets —
  assert it, do not assume). With ONE draw shared by both faces the defect is
  unobservable/arguable (the faces ARE identified under the quotient, so the
  identity looks defensible); use a SECOND independent seed and ship a LIVE test
  asserting the two probes differ. Config-blindness, one level up: the
  convenient single-draw fixture nulls the term.
- **A strict xfail for a not-yet-built phase needs a FLIP-PROOF, because it will
  not flip by itself when the fix changes the CALL rather than the body.**
  Production's periodic op is an angular identity fed the wrong half-trace, so
  after B3.4c the test's `apply` line changes and the marker is deleted by hand.
  Ship a LIVE companion showing the body already reproduces the reference when
  fed the partner's Γ₊ — it is simultaneously the flip-proof (it names exactly
  what the phase changes) and a body pin (a scaling/permuting/averaging operator
  fails it), and it is NOT `x == x`. Keep every premise (space cross-check,
  pairing, probe activation) in LIVE tests so the xfail body has exactly ONE
  failable statement; verify with `--runxfail` that it reds on THAT line.
- **Do not add cases that duplicate an existing anchor.** The two B3.4b closures
  were rejected as new cases: `test_reemission_closure.py` already anchors both
  against hand-written expressions over 5 quadratures × 3 α, and the only thing
  a case adds is the frozen artefact — already delivered for the specular
  closure by one extra method on an existing case.
- **Store the face and BOTH index sets on EVERY artefact** (pre-migration only
  the vacuum case carried one). It makes each file self-describing AND turns the
  Layer-1 classification into a frozen independent statement the harness
  cross-checks. Honest scope: MEASURED, that cross-check is BLIND to the
  `> 0.0` tangential twin on this case set (lebedev/LS/GL carry EXACT-zero
  tangential cosines; only `product(2,4)` has round-off ones) — name the blind
  spot and point at the file that does cover it, rather than adding a case.
- **Mutation table (8 in-process monkeypatches, control 43 passed/1 xfailed):**
  local-remap→`arange` (2 red), α-fold dropped (3), cosine weight dropped (3),
  zero-map→echo (1), domain built from Γ₋ (6), local-perm reversed (4),
  normal flipped (16), `TANGENTIAL_EPS→0` (0 — the documented blind). The
  `arange` mutation is quadrature-coincidence-blind on ls4/ls6; a REVERSED
  local perm is the coincidence-proof specular mutation — use it to give a
  rank-N row's specular leaf its own red.

---

## L33 — When the phase's PRODUCTION lands while you plan: the deferred gate whose FLIP IS A NO-OP, and the green gate whose REASON decayed (B3.4c, periodic → the partner face)

Two failure shapes, both measured, both new. The context: B3.4c narrows the SN
periodic law so its domain is the PARTNER face's Γ₊ (`γ₋ψ|_f = γ₊ψ|_{f'}`), and
its plan-of-record shipped a strict-xfail "todo list" naming the phase.

### 1. ⭐ The MARKER WHOSE FLIP IS A NO-OP — a fifth Mode-8 class-4 shape

L32 taught "an xfail for a not-yet-built phase needs a FLIP-PROOF because it
will not flip by itself." That is right and it is not enough. The prescribed
flip — *"when it lands, hand the operator the PARTNER's Γ₊ and delete the
marker"* — turned the xfail into a **character-for-character duplicate of the
LIVE flip-proof sibling next to it**, a row green BEFORE and AFTER the phase.
Nothing about the edit is sensitive to whether the production change landed, so
deleting the marker would signal "phase complete" while asserting nothing new.
**MEASURED, decisively:** the production step landed mid-session and the row
still red byte-identically — same 735/735, same first-mismatch values. So a
"B3.4c gate" survived the entire B3.4c production change untouched.

**The discriminator, and it is cheap:** *diff the xfail body against the live
companion you shipped as its flip-proof.* If applying the documented edit makes
them textually equal, the flip is ceremony. **Rule: an xfail's flip-edit must
touch a statement whose VALUE the production change determines.** The repair
here was to state the claim against the production ANSWER rather than a
hard-coded face — `source = law.geometry_map.domain_face(face)`, then index the
probe dict BY that answer, so a regression makes the gate select the wrong probe
and red. And the composite-level claim (the thing the phase actually builds) is
a NEW test in the operator's own file, never an edit of the law-level row: a
harness that only ever realizes a law against a hand-built method space cannot
state a claim about the COMPOSITION.

### 2. ⭐ The gate that stays GREEN while its REASON becomes false

`B` was documented block-DIAGONAL over faces; a wrap makes it block-STRUCTURED.
Three test rows assert it — and `[M]` **all three stay green**, because every
one is posed on a reflective-only fixture. What decays is the *justification*
("`B` is block-diagonal over faces, **so** the face-subset reflect MUST be the
exact restriction" — the Gauss-Seidel schedule's exactness argument). Nothing
in the test, the tag, or CI notices; it is Mode-8 class 7 applied to a REASON
rather than a fixture. **Rule: when a phase falsifies a structural claim, grep
the claim's WORDS in tests/ as well as its symbols — the rows that assert it on
a now-special-case fixture must be re-scoped from "`B`" to "these laws" IN THE
SAME CHANGE, and the new structure gets its own positive gate.**

### 3. The plan's OWN premises had already moved — steps 1/2/3/5/6 were on disk

Five of seven plan steps were implemented, uncommitted, when the plan began, and
step 4 landed mid-session. **Re-run `git status --porcelain -- orpheus/` at the
START and again before every claim; stamp the plan with the tree state.** A
verification plan written against a stale premise gates the wrong thing.

### 4. Recipes worth keeping

- **Ship the plan's gates as a RUNNABLE dry-run module**, not only as fenced
  code. Transcribing them exposed two errors no amount of reading would: the
  perturbation form `B(base+e) − B(base)` is NOT bit-exact (`[M]` 1 ULP, the
  `(x+1)−x` cancellation — exploit LINEARITY with a zero base instead), and one
  negative's refusal had moved from the realizer to the composite.
- **When a tag registry refuses the law, the gate still exists one layer down.**
  `BC("periodic")` raises (`SNMesh.BOUNDARY_OPERATOR_REGISTRY == {"vacuum",
  "reflective"}`), but `sn.realize_boundary_law(law, face)` installs it — that
  is the ONLY layer at which an end-to-end claim is stateable pre-#189, and it
  must be named as such in the fixture's docstring.
- **A predicate flip on a FACTOR may have zero downstream observable.** `[M]`
  `SpatialWrap.is_adjointable: False → True` changed nothing: the composite
  reads each REALIZED operator (identity body ⇒ already `True`). Gate the
  factor, state the scope, and put the teeth on a transpose SUPPORT gate.
- **Reciprocity `⟨Bx,y⟩ = ⟨x,Bᵀy⟩` is a CONSISTENCY check, not a correctness
  one.** `[M]` the pre-fix operator reciprocated at rel 1.15e-16 because forward
  and transpose were wrong the SAME way; reverting the fix does not red it. It
  catches only the forward/transpose MISMATCH (`[M]` rel 1.5e-2…2.3e-1). Pair it
  with an object-level SUPPORT gate (feed a state on Γ₋(f), assert the deposit
  lands on Γ₊(f′) and the own face stays zero).
- **A slab is the degenerate two-face case for any partner map** — "the
  opposite face" and "the other face" coincide, so a 2-D companion is mandatory.
  `[M]` a cross-axis partner is SHAPE-LEGAL (`|Γ₋(xmin)| = |Γ₊(ymax)| = 12` on
  LS4) while the index sets differ, and a face-NAME permutation certification
  passes it.
- **Mutation battery (7 monkeypatches, control 95/95 green):** `domain_face →
  face` (25 red), `face_opposite → identity` (18), suffix↔sign swapped (12),
  revert the composition's source face (7), wrap body ×2 (6), factor
  `is_adjointable → False` (2). **Three predicted blindnesses CONFIRMED:** the
  revert does NOT red reciprocity; `face_opposite → identity` does NOT red the
  involution row (the identity IS an involution — the NO-FIXED-POINT row is the
  catcher); the suffix swap does NOT red either round-trip row (a bijection's
  round-trip is invariant under relabelling both directions — a transcription of
  the RETIRED expression is the catcher). A convention's round-trip is never
  sufficient; pin the ABSOLUTE table.
- **A same-sized-wrong-SET probe is the only way to gate "compare sets, not
  sizes".** `[M]` `replace(space, inflow_indices=space.outflow_indices)` → the
  real guard raises, a size comparison passes. No mesh produces that input, so
  the discipline is otherwise an unverified claim however well its docstring
  argues for it.

---

## L34 — Making a SYMMETRY EXACT is not a ≤1-ULP change: exactness manufactures TIES, and a downstream sort that was accidentally total becomes under-determined (#325, roots_of_unity)

The framing that arrives with such a change is always *"the node values move
by ≤1 ULP, so anything pinned bit-exactly re-baselines."* That sentence is
true about the **values** and can be catastrophically false about the
**results**, because the whole point of the change is to make previously-
distinct floats **bit-equal** — and anything downstream that ORDERS, GROUPS,
DEDUPES, or takes an `argmin`/`argsort` over them was relying on the noise to
break ties.

`[M]` #325: `roots_of_unity` makes azimuthal mirror partners share `η`
bit-exactly. `rules_product.py` sorts each μ-level by `η` with a default
`np.argsort`. Distinct-η per level at `n_phi=64` drops **60 → 33**; the
per-level ordering changed in **36 of 36** `(n_mu, n_phi)` configurations,
including those with no axis node. End-to-end on a heterogeneous 2-G
cylindrical fixed-source at `n_phi=32`, with the node values held
**bit-identical** and only the ordering varied: **1.008 % flux change**; and
`kind="quicksort"` vs `kind="stable"` — a numpy implementation detail —
differ by **1.775 %**. The node-value drift in the same run was **1.06e-14**.
Twelve orders apart, in the same commit, under one justification sentence.

**The checklist when a change makes something exact:**
1. Grep the consumers for `argsort` / `argmin` / `sort` / `unique` /
   `set(` / dict-keying **on the quantity being made exact**. Each is a
   tie-break that was being decided by the noise you just removed.
2. **Is the sort key INJECTIVE?** `η = sinθ·cos φ` is 2-to-1 over
   `φ ∈ [0,2π)`. A non-injective key means the ordering was never determined
   *by the physics* — exactness only reveals it. That is a live latent defect
   to be ruled on, not a side effect of your change.
3. **Does the ambiguity converge away?** `[M]` it did not: the gap sat flat at
   ~1.0–1.8 % across `n_phi ∈ {32, 64, 128}`. A flat-in-refinement gap is a
   defect; a shrinking one is discretization.
4. **Split the commits.** The ordering ruling lands FIRST, alone, with its own
   reference. Otherwise the value-change commit's DriftWarning run is
   uninterpretable and its justification record is a lie.

**The Mode-12 trap specific to ordering:** *"the level is sorted by η"* is the
wrong functional — sortedness is invariant under permuting equal elements, so
the two orderings that differ by 1.8 % are BOTH sorted. Gate the **full index
tuple** against an independently-constructed one, plus a
**`kind=`-invariance** row (`quicksort`/`stable`/`heapsort`/`mergesort` must
give bit-identical tuples), which is the operational proof the key is
injective. And the canonical curvilinear diagnostic is no help: `[M]` the
homogeneous flat-flux `Q/Σ_t` leg gave the SAME value on all four ordering
legs — the config that AGENT.md §0.6 calls the single most powerful
curvilinear diagnostic is designed-green here, because flat flux nulls exactly
the redistribution the ordering perturbs.

### L34b — the two blind involutions, and the "is it a live bug?" discipline

`[M]` Two measurements from the same session, both worth reusing:

- **A nearest-neighbour partner search is not automatically a bug.** The
  hypothesis "NN matching over a noisy set must be mis-pairing" was
  **REFUTED**: the map is a valid permutation AND involution up to
  `n_phi=1024`, because the mis-pair margin is the node separation
  (`5.0e-3`) versus a `1e-16` perturbation — 13 orders. Measure the margin
  before filing it as a bug; and **KEEP the search** rather than replacing it
  with an index formula, because the sibling family (`lebedev`) has no such
  formula and the replacement would mint a twin path AND delete the only
  detector of the next item.
- **`ref[ref] == id` passes on a residual-0.94 garbage map.** `[M]` for odd
  `n_phi` the x-mirror partner does not EXIST; the unguarded `argmin` returns
  the nearest node anyway and the result is still a permutation and still an
  involution. The only functional outside that stabiliser is the
  **RESIDUAL** `max‖R·x_n − x_{ref[n]}‖`. Same shape as the `face_opposite`
  lesson (L33): an involution law alone is blind because the identity is an
  involution — here, because *any* self-inverse pairing is.

### L34c — "≤1 ULP agreement with `np.cos`" is the wrong reference, not just a wrong number

`[M]` The shipped construction is `3.75 ulp(1.0)` from `np.cos(2πp/q)` — the
stated acceptance criterion FAILS. But `np.cos(2πp/q)` is not the true value:
`fl(2πp/q)` already carries argument-reduction error. Against 100-digit
mpmath the new code is `0.57 ulp(1.0)` and the legacy is `3.72 ulp(1.0)`.
**Gating a more-accurate implementation against the less-accurate one it
replaces gates it against the error it exists to remove.** State the criterion
against the arbitrary-precision value, and state it honestly: `[M]` the legacy
is strictly closer on `30/1040 = 2.9 %` of components (by ≤ `1.11e-16`), so
the claim is *"within 0.57 ulp everywhere"*, never *"always closer"*.

### L34d — the harness that says "N gates are blind" is the one to distrust first

`[M]` My first teeth harness reported **18 of 18 gates never red**. The
mutation had never bound — `pytest.main` imports its own copy of the test
module, so patching a privately-loaded copy does nothing. The corrected
harness patches the SOURCE module, purges the test module from `sys.modules`,
and carries a **positive-control probe** ("did the freshly-imported test
module bind the mutant?") printed on every row. With it: 14 of 17 gates red
from the body battery, the remaining 3 from targeted mutations. This is the
`vv-principles` Mode-8 METHOD WARNING and it cost a full false verdict — an
all-blind result is far more likely to be a broken harness than a broken
suite.

Corollary, and it caught me: the same harness flagged a gate **I had just
written** as never-reddening — an "activation guard" asserting odd `q` has no
node at `π` (`not any(2p == q)`), which is a **theorem about parity**, true
for every input. It survived authoring and a full green run. Run the teeth
harness over your OWN new module before delivering it.

### L34e — a "bit-identical to the legacy adapter" pin whose legacy adapter is gone

`[M]` `test_product_bit_identical_to_legacy_adapter` compares
`product_mu_phi(...)` to `Quadrature.product(...)` — and `Quadrature.product`
CALLS `product_mu_phi`. The named legacy adapter
(`orpheus.sn.quadrature.ProductQuadrature`) no longer exists. Replacing
`product_mu_phi` with **random nodes** leaves the test **green**, while its
docstring claims *"Pins the cylindrical regression snapshots"*. A snapshot
whose only "independent" pin is that test has **no anchor at all**, and
re-baselining it against its own old value is the contaminated re-baseline
L32 is about. **Whenever a pin names a "legacy"/"reference"/"adapter"
counterpart, grep that the counterpart still EXISTS and is not the SUT
reached by another name** — a delegation added later silently converts the
pin into `X == X`.

## L35 — verifying a pure-math PRIMITIVE (a group/algebra type) against pure math

From the G2 gate spec for `orpheus/geometry/transformation.py` (rigid motions
`E(d)`), written PRE-carve — and the module LANDED mid-plan, which is half the
lesson. Full plan: `scratch/g2_verification_plan.md`.

### L35a — the reference pillars for an ALGEBRA are different from a solver's

There is no MMS row and no semi-analytical row: a group-theory type has no
differential operator to manufacture a solution for and nothing reduces to a
quadrature. Every row is **closed-form**, and the structurally-independent
grounds available are: (1) **SymPy under an EXPLICIT unit parameterisation**
(`n̂ = (sinθcosφ, sinθsinφ, cosθ)` — `[M]` imposing `Σnᵢ²=1` by `subs` after
expansion **does not fire** and leaves `4n₀²(Σnᵢ²−1)` residuals); (2) an
**external implementation with a different ALGORITHM** (`scipy…Rotation.
from_rotvec` is quaternion-based, so it is not a re-spelling of Rodrigues);
(3) the **Lie definition** `R = expm(θ(vuᵀ−uvᵀ))` — Padé, structurally unlike
every closed form, and the only **dimension-generic** reference (`[M]` 3.9e-14
vs the shipped constructor, so gate at 1e-12, not 1e-15); (4) **published
tables** (group orders; the cube's Wyckoff site symmetries); (5) **exact
integer arithmetic** (a permutation's bijectivity, `π(gh)=π(g)∘π(h)`,
`ord(g^k)=n/gcd(k,n)`) — this last one needs no reference at all and is the
strongest class available.

### L35b — the group-action HOMOMORPHISM is the deepest cheap gate

`π(g∘h) = π(g)∘π(h)` on an induced permutation cross-checks **two independent
layers** — the element algebra and the nearest-neighbour action on points —
using only integers. It simultaneously pins the composition order, the
row-vs-column action convention, and `π` vs `π⁻¹`. `[M]` the wrong composition
order violates it on **102 of 144** cube/`O_h` pairs. It is **VACUOUS on an
abelian group** (`C_n`: 0 violations either way) — the fixture must be
non-abelian. A checker that computes `π` and returns a `bool` (ERR-073's
shape) makes this law unaskable; returning the permutation makes it free.

### L35c — an involution/order law is Mode-12 BLIND to the AFFINE part

`(σ, t)² = (I, σt + t)`, and `σt = −t` for **every** `t ∈ span(n̂)`. `[M]` all
of `t ∈ {d n̂, 2d n̂, 4d n̂, −2d n̂}` are involutions, while the true fixed
plane moves by `0.37 / 0.74 / 1.48`. Likewise a rotation's ORDER is blind to
its centre. **The only catcher is "fixes its named fixed set POINTWISE"** —
gate the FIXED SET, never the order, for anything affine. Generalisation: for
any type carrying a decomposition `(linear, affine)`, enumerate which laws
factor through the linear part alone — those are designed-green on the affine
half.

### L35d — the SEAT is a theorem, not a convention (and it is one gate)

A `G`-preserved weighted point set has a `G`-FIXED centroid (3-line proof from
`w_{π(i)} = w_i` + linearity). `[M]` a cube shifted to `(2.5,−1.25,0.75)`:
**48/48** `O_h` elements *seated at the centroid* preserve it; **1/48**
unseated do. That single row converts "where do we put the origin?" from a
modelling choice into a computed fact, and it is the abstract form of a live
defect (`Mirror('x').is_invariant(mesh)` reads False for production meshes
solely because they start at `origin = 0.0`).

### L35e — bijectivity and the match WINDOW are INDEPENDENT failure modes

The ERR-073 injectivity guard does **not** make a nearest-neighbour certifier
sound. `[M]` a two-point set off-symmetry by **1e-9** certifies under a
**1e-7** window with a *perfectly injective* π. So the window is a
first-class correctness parameter: it must be an **explicit argument** (a
module constant makes the "the window bites" gate SIGNATURE-tautological), and
it should default relative to the point set's **minimum pairwise separation**
— a computable, set-intrinsic quantity — rather than to an absolute constant.

### L35f — a `-> T | None` collapses N guards into one value; isolate from the INPUT side

`permutes` returns `None` for "image not in the set" OR "not injective" (and
`preserves` adds "weights disagree"). A negative gate asserting `is None`
proves only that *some* guard fired. Where the API cannot be changed, prove
the isolation **from the input**, in the test, with plain numpy: `[M]` the
ERR-073 duplicate-node witness has guard1 passing **48/48**, guard3 **48/48**,
guard2 **0/48** — a clean isolator; and an unequal-weight antipodal pair makes
`permutes` SUCCEED while `preserves` FAILS, isolating the weight guard (which
is **VACUOUS on any equal-weight fixture**, i.e. on every shipped quadrature).

### L35g — ⛔ do NOT assert TIGHTER than the type's own invariant

`RigidMotion.__post_init__` accepts `max|QᵀQ−I| ≤ 1e-10`, so an element
carrying a 1e-11 shear is a **legal value of the type**. A gate asserting
`1e-14` on *any* element asserts a property the type does not promise and
would red on a legal instance. Split it: one row for **the type's invariant**
(and its rejection threshold, via the production constructor), a separate,
stronger row for **the named constructors' actual quality** (`[M]` ≤1e-14, and
`signed_permutation` **exactly 0** — integer entries). Conflating them leaves
the constructors ungated at their real quality while over-claiming the type's.

### L35h — bit-exactness is EARNED PER LAW; measure before choosing the assertion

`[M]` on the same type: identity `e∘g=g∘e=g` **500/500 bit-exact**;
associativity **500/500** on signed permutations but **0/500** on general
rotations; `σ²=I` bit-exact for a COORDINATE normal but **0/300** for a
general one; `g∘g⁻¹` **0/500**; seating fixes `c` only **58/500**. And a
seated reflection of a cell lattice lands its images **5.6e-17–2.2e-16** off
their partners while the induced PERMUTATION is exactly `arange(n)[::-1]` —
integers exact, coordinates not. Pick `array_equal` vs `atol` per row from
measurement; a uniform choice is either a false red or a thrown-away gate.
The coordinate-vs-general normal split is also a **size-blindness** trap in the
argsort-n family: a coordinate-only fixture tests the degenerate diagonal
sign-flip and is structurally blind to the general Householder path.

### L35i — the module can LAND while you are planning it (and that IS the deliverable)

`[M]` the brief said `transformation.py` did not exist; it was written during
the session (untracked, mtime mid-plan). Re-checking the plan against the
shipped code re-scored **10 of 15 API findings as DISCHARGED** and surfaced
**3 new gates the landed design earns** (two spelled actions
`on_points`/`on_directions`; a circle-point rotation primitive that preserves
exactness where an angle cannot; `element_order() -> int | None` where `None`
is the glide/screw law). **Before delivering any pre-carve plan, `ls` the
target path** — the reconciliation is higher-value than the original matrix,
and an unreconciled plan silently specifies a dead API. Same family as the
B3.4c finding (L33) where all seven production steps landed mid-plan.

### L35j — WRITING the gates found three defects, and every one was in the FIXTURE

G2 shipped as `tests/geometry/test_transformation.py` (42 gates / 96 cases,
8.9 s, 32/32 mutations caught, 0 blind). The three reds during authoring were
all mine, and each is a reusable trap:

1. **A single Gram–Schmidt pass is not orthonormal enough for a 1e-12 gate.**
   `[M]` `max|GGᵀ − I| = 2.1e-12` on near-parallel random draws, so the
   production constructor *correctly refused my fixture*. When a shipped guard
   rejects your test input, suspect the fixture first. Build orthonormal
   frames with QR.
2. **An ABSOLUTE tolerance is the wrong contract for a translation residual.**
   It scales `O(ops × ‖t‖ × eps)`, so a fixed `atol` reds on correct code for
   large-‖t‖ draws (`[M]` 2.1e-14 against a 1e-14 bound) — while being too
   loose for small ones. Normalise by `max(1, ‖desired‖_∞)` and MEASURE the
   normalised worst case (`[M]` ≤6.9e-15 over 3000–4000 draws per law).
3. **State the law in the direction that is a FLOAT theorem.**
   `on_points − on_directions == t` is NOT exact (`fl(a+t) − a ≠ t`; `[M]`
   2.2e-16) but `on_points == on_directions + t` IS bit-exact `[M]` 6000/6000,
   because it recomputes the same expression. The true direction is also the
   stronger assertion — check both before settling for a tolerance.

### L35k — the SELF-CONSISTENT gate is blind to a transposed action; the homomorphism is not

`[M]` The mutation `x @ Q` for `x @ Qᵀ` reddens **O2 (the homomorphism) but
NOT O1 (the positive `permutes` gate)**. O1 asserts `g.on_points(P) == P[π]`
where `π` was *built by the same mutated action* — and a transposed action
still maps a group-invariant set onto itself (`Qᵀ` is in the group), so O1 is
self-consistent and green while returning the WRONG permutation. Only
`π(g∘h) = π(g)∘π(h)` catches it, because a transposed action is an
ANTI-homomorphism. Generalisation: **a gate that builds its own reference
through the mutated path can only see errors that break self-consistency** —
pair every "the map is correct" row with a COMPOSITION law it must satisfy.

### L35l — the harness lied first, again (L34d recurrence, and the control caught it)

`[M]` The G2 battery's first run reported **32/32 BLIND** while the summary
lines plainly read `23 failed` / `63 failed` / `38 failed`: the parser wanted
`FAILED` lines that `-q --tb=no` never emits (needs `-rf`), and ANSI colour
codes broke the match. The POSITIVE CONTROL is what exposed it — a mutation
making `reflection` return `+I` cannot leave 42 gates green, so "0 gates" was
a *harness* verdict. Cost: one run. **Never ship a battery without a mutation
whose expected blast radius is large enough that a zero reading is absurd.**

---

## L36 — TYPE-COLLAPSE carve (N types → 1 parameterised type): the class-level gates decay silently

Source: G5 verification plan (`scratch/g5_verification_plan.md`), PRE-carve,
`refactor/operator-strategy-layers` @ `af99f064`. The carve: collapse
`IdentityMap` + `SpecularMirror` (`orpheus/geometry/boundary/_factors.py:443,474`)
into one `SelfPairedDeck(motion: RigidMotion)`.

### L36a — the archetype's signature failure: class-level assertions become tautologies

A type-collapse moves the discriminating information **from the TYPE to a
FIELD**. Every existing gate that asserts on the type keeps passing and stops
discriminating. Two measured instances in ONE file
(`tests/geometry/test_boundary_factors.py`):

* `test_every_production_law_states_both_factors` asserts
  `isinstance(law.geometry_map, geom_cls)` per law. Today: `IdentityMap × 5`,
  `SpecularMirror × 2`, `SpatialWrap × 1`. Post-collapse **7 of 8 rows expect
  the same class**, so "a vacuum law that mistakenly declared a mirror" is no
  longer distinguishable.
* `test_specular_mirror_is_the_only_ordinate_permuting_geometry` asserts
  `permuting_geoms == {SpecularMirror}`. Post-collapse the permuting set and the
  non-permuting set are the SAME class — the gate reduces to "the class is among
  the classes that permute", satisfied by any non-empty subset.

**RULE: for any N→1 type collapse, inventory every `isinstance` / `type(...)` /
class-set assertion over the collapsing family FIRST, and re-pose each onto the
PARAMETER in the same commit that lands the collapse.** A commit that lands the
collapse alone has deleted those gates while the suite stayed green. The
re-posed form asserts the parameter's VALUE (here: the exact signed-diagonal
matrix per law, `array_equal` — bit-exactness is earned, `[M]` 12/12 exact).

Cheap companion the archetype earns for free: `type(A().factor) is
type(B().factor)` — `is`, not `isinstance`. It is the collapse's headline claim
as one line, and it is the ONLY gate that catches "collapsed" into a base class
plus two subclasses.

### L36b — the naming invariant is necessary-not-sufficient, and the missing clause has a concrete inhabitant that does the DEFERRED type's job

The brief's insight — "a face paired with ITSELF has an involutive pairing" — is
true forward and was used as the guard `element_order() in (1,2)`. The CONVERSE
is false. `[M]` (probe p1) the involutions of `E(3)` are FOUR families:

```
identity        order=1 det=+1 fix=3      self-paired  ✔
reflection      order=2 det=-1 fix=2      self-paired  ✔
half-turn       order=2 det=+1 fix=1      maps a face to its OPPOSITE  ✘
inversion (-I)  order=2 det=-1 fix=0      maps a face to its OPPOSITE  ✘
```

So the sketched guard admits exactly the elements that do the **deferred**
type's job (`SpatialWrap`), re-opening the illegal state the carve exists to
close. Correct guard is a conjunction; the second clause (`fixed_subspace_
dimension ≥ dimension − 1`) is the load-bearing one, and it is dimension-generic
(`[M]` d=1: reflection fix=0=d−1 admitted, inversion(1) IS the mirror; d=2:
inversion(2) is a half-turn, det=+1 fix=0, rejected).

**RULE: when a carve's guard is ONE algebraic property named after the concept,
enumerate that property's full conjugacy classes and check each against the
SEMANTIC claim.** Then build the two-way witness pair — here a glide
(`T[0,1,0] ∘ refl_x`: `fix=2` ✔, `order=None` ✘) and `inversion(3)`
(`order=2` ✔, `fix=0` ✘) — so each clause has a mutant that only IT catches,
with DISTINCT `match=` strings. A single generic message lets one clause be
deleted with the other's message still matching.

### L36c — a `RigidMotion`-carrying type's TRANSLATION is bit-identically invisible to every angular functional

`RigidMotion.on_directions` **drops the translation by construction**
(`transformation.py:515-518`). `[M]` `reflection(normal=e_x, offset=off)` gives
an image bit-identical to `offset=0` at `off ∈ {0, 2.5, −17.0}` on
`lebedev(17)`. So a wrong mirror-plane POSITION is designed-green at the
permutation level, the realized-image level, the frozen-snapshot level and any
scalar — no tolerance, refinement or regime can expose it.

Closure is the ERR-067 pattern: **repair the TYPE, not the gate.** Requiring
`motion.is_linear` at construction makes the offset unspellable and turns the
blindness into a real gate. `[M]` the stronger guard `is_linear ∧ fix ≥ d−1`
also *implies* involution (a linear orthogonal `Q` fixing a hyperplane pointwise
is `I` or a reflection ⟹ `Q²=I`), so it needs no `element_order` and carries no
`atol` — the involution property demotes to an asserted theorem. Prefer it.

Sibling permanent blindness worth recording: for ANY involution
`perm == argsort(perm)` identically, so forward-vs-inverse permutation is
unfalsifiable by theorem. Good news for `is_adjointable`; it means the day a
NON-involutive deck map lands (#178), the transpose direction becomes falsifiable
for the first time and needs its own gate.

### L36d — the level lattice: bool / matrix / permutation / certified-table are pairwise incomparable

The recurring question "is a permutation-level gate blind to what a matrix-level
gate catches, and vice versa?" has a clean answer for this shape — **neither
subsumes the other, so all four levels are required**:

| level | sees | blind to |
|---|---|---|
| bool (`permutes_ordinates`) | that *some* relabeling happens | which axis; the offset; whether the quadrature admits it |
| matrix (`motion.linear`) | which mirror, bit-for-bit | whether the quadrature is CLOSED under it |
| permutation (`reflection_index`) | the right relabeling of THIS quadrature | the offset; forward-vs-inverse; a local remap unless the fixture REVERSES order |
| certified table (ERR-042/044/045) | involutive ∧ measure-preserving ∧ inflow→outflow | which mirror was declared (its canonical mutant IS the identity table) |

Two operational corollaries. (i) `domain_face`-style accessors that don't read
the parameter have the WHOLE group as their stabiliser — such a gate can never
be evidence about the parameter; its value is entirely in the cross-type
partition (self-paired ⟺ not-a-wrap), so it MUST ship with the sibling type's
no-fixed-point control. (ii) The certified-table row whose canonical mutant is
the identity table (`ERR-045`,
`tests/geometry/test_bc_universal_invariants.py::TestReflectiveInflowToOutflowInvariant`)
is the row a merge-identity-with-mirror carve must not let decay — it is the gate
saying the two are still distinguishable where it matters.

### L36e — measure the brief's cost estimate before rationing the battery

`[M]` The brief budgeted "`tests/numerics` + `tests/geometry` ≈ 5.5 min". The
subset the carve can actually reach —
`tests/geometry tests/numerics/test_quadrature_directional.py tests/numerics/test_face_layout.py`
— runs in **9.40 s** (`1 failed, 743 passed`, the 1 being the pre-declared
task-#33 red). Off by ~35×, in the helpful direction: a 16-mutation battery goes
from "ration it" to ≈6 min. **Measure the reachable subset, not the directory
the brief names.** Same discipline as refuting an optimistic premise — an
over-stated cost silently shrinks the battery, which is the same loss as a blind
gate but harder to see.

### L36f — a `pytest.raises(ValueError)` written against a `ValueError` SUBCLASS is a false gate that reads correct

`class BoundaryError(ValueError)` (`orpheus/geometry/boundary/_errors.py:44`).
The defect being fixed is that a certification path raises a *bare* `ValueError`
instead of `BoundaryError` — so `[M]` `except ValueError` catches it today and
`except BoundaryError` does not. A gate spelled `pytest.raises(ValueError)` is
therefore GREEN BEFORE THE FIX and stays green after: total loss, one word, and
it looks right in review. **When the fix is "raise the project's typed error
instead of the builtin", the gate MUST name the subclass, and the pre-fix state
of that gate MUST be RED** — ship it `xfail(strict=True)` so the fix's landing
is what deletes the marker.

---

## L37 — OPTIONAL-BINDING carve (a `None` default that silently disables a derivation): the gate you cannot write, and the sentinel that means two things

**Campaign:** ORPHEUS `geometric_transformation_consolidation` G6 — plan
`scratch/g6_verification_plan.md`, written PRE-carve 2026-08-04 at `2ca7dd86`.
Scope: `LinearOperator.domain`/`.codomain` are `Optional[FunctionSpace] = None`;
when `None`, `_AdjointOperator.apply` (`numerics/operator.py:1223-1226`) **skips
both metric applications**, so `A† = G_V⁻¹AᵀG_W` degrades silently to the bare
Euclidean transpose. G6 mints `Γ±(face)` spaces, binds the boundary tier, makes
binding mandatory, moves traversal onto the space.

### a. The carve archetype and where its keystone lives

An **optional-binding** carve is not a value carve. Nothing numerical changes on
the forward path (`apply` never reads `domain`), so a bit-identity keystone is
available and CHEAP — and it is worthless alone, because the whole point of the
carve is the ONE path where numbers DO change (`.H`). The keystone is a
**metric-sensitive reciprocity gate**, and it must be written BEFORE the metric
is installed so it goes RED→GREEN. Written after, it can only ever be green, and
a green metric gate proves nothing about whether the metric does anything.

Corollary the plan states as non-negotiable ordering: **the red→green transition
IS the evidence the metric step is not decorative.** Nothing else produces it.

### b. ⭐ The Mode-12 stabiliser was FOUR cases wide — measure the whole tier, not the one operator

`vv` Mode 12 says compute the commutator `[G, Aᵀ]` at design time. Measured on
the ORPHEUS SN boundary tier, EVERY shipped law except one is in the stabiliser:

| realized law | operator | `[G, Aᵀ] = 0`? |
|---|---|---|
| vacuum / any re-emission at α=0 | `ZeroOperator` | yes — `G⁻¹0ᵀG = 0` |
| reflective, albedo+specular | `PermutationOperator` | yes — ERR-042 weight preservation, **`G₋ == G₊[local]` bit-identical** on GL(4)/product(4,4)/lebedev(17) |
| periodic | index-identity wrap | yes — opposite normals ⟹ same `\|Ω·n\|` |
| prescribed inflow | `IncomingSourceOperator` | rank-0, apply-only |
| **white / albedo+isotropic** | `AngularAverageOperator` | **NO** — `‖A†−Aᵀ‖/‖Aᵀ‖` = 0.209 / 0.684 / 0.612 |

⟹ the adjoint is "right by accident of the one case anyone checked", and the
accident is four cases wide. **A design-time stabiliser audit must enumerate the
WHOLE family the gate could be posed on, not just the operator in front of you** —
otherwise "test the reflective BC harder" reads like a cheaper alternative when
it is provably a non-catcher.

Measured margins for the one catcher (`probe_rank1.py`, lebedev 17, spaces bound,
`outer(1_{Γ₋}, InnerProductFunctional(cos_w/norm))` carrying PRODUCTION `cos_w`):

```
⟨Ax,y⟩_G−  vs  ⟨x, A.H y⟩_G+   rel err 2.2e-16   <- the LAW
⟨Ax,y⟩_G−  vs  ⟨x, Aᵀ  y⟩_G+   rel err 5.5e-1    <- the DEGRADATION
mutation domain:=None  ⟹  A.H y  array_equal  Aᵀ y   (exactly)
activation: trace metric max/min = 1.351 / 3.468 / 5.601  (GL4 / product / lebedev)
```

### c. ⭐ The MUTATION found what code-reading could only assert

Positive control PC1 (`FunctionSpace.apply_inverse_metric := identity`, a
linear/shape-preserving mutation — anti-#18 clean) over
`tests/sn/operators + tests/geometry`:

```
baseline    3 failed, 1597 passed   23.2 s
PC1        33 failed, 1567 passed   23.7 s      ⟹ 30 NET REDS
  19 test_g_adjoint_reciprocity.py · 7 test_inverse_adjoint_coherence.py
   3 test_radial_characteristic_metric.py · 1 test_psi_half_coupling.py
```

**All 30 are bulk / SN-leaf / radial-characteristic gates. NOT ONE narrowed
boundary-law operator reddened** — because `domain is None` short-circuits the
call, so the metric mutation is *unreachable* there. That is empirical proof of
the blindness, obtained by mutation rather than by reading `if x is not None`.

**Reusable technique: run the positive control BEFORE writing the gate, and read
WHICH files it misses.** The absent reds map the blind region exactly, and they
are far more persuasive in a plan than a code citation.

### d. ⭐ A sentinel that encodes TWO states makes the discriminating gate
UNWRITABLE — the type must fix it, not a test

`domain = None` means BOTH "this operator is space-GENERIC by mathematics"
(`IdentityOperator`, `ZeroOperator`, `DiagonalOperator`, `PermutationOperator` —
an identity is the identity on every space) AND "nobody bothered to bind it".
**No gate can distinguish them**, because the two states are the same value.
That is *why* the degradation is silent, and it is why "an unbound operator
cannot be constructed" is the wrong acceptance line for the generic leaves.

The verification verdict is a TYPE demand, not a test: mint an explicit
space-generic marker so the escape is DECLARED. Generalisation — **when an
acceptance criterion asks a gate to separate two states that share one runtime
value, the honest deliverable is "this gate cannot exist; here is the type that
makes it exist", not a cleverer assertion.**

### e. The survey that sizes a "turn on a check that never ran" risk

Instrument the composer constructors and COUNT the skips before binding
anything. Plugin `skipcount2.py`, `tests/geometry + tests/sn/operators`, 23 s:

```
OperatorSum 172 SKIP / 1694 CHECK (9.2 %) · OperatorProduct 6 / 1568 (0.4 %)
178 skipping compositions, 10 distinct operand shapes
  148  Sum(IsotropicScattering, IsotropicN2N)          both sides fully anonymous
    5  Product(BulkAnalysisOperator, SweepOperator)    LEFT CODOMAIN ONLY  <- the dangerous one
```

**The half-bound composition is the highest-risk row**, not the biggest one: the
bound side already commits to a space and nothing has ever compared it to the
other, so a genuine mismatch can be hiding there. A fully-anonymous pair usually
just needs the same space on both sides.

### f. Two acceptance lines that were arithmetically / structurally impossible

1. **"permute the layout and assert BIT-IDENTICAL output."** Measured: a
   49-term reduction under 200 random permutations is bit-identical **50/200
   (25 %)**, worst drift 4 ULP; even `sum(cos_w)` alone is bit-stable 146/200.
   FP addition is not associative — **no architecture can make a REDUCTION
   permutation bit-identical.** The honest replacement is a permutation that
   reorders no addition: the **face packing order** in the flat-buffer layout.
   Measured bit-identical 4/4 per face (slices + index sets) while the flat
   buffer moves — a stronger discriminator, because it reds exactly the
   consumer that reads a self-computed offset instead of a name.
   **Rule: before accepting a bit-identity acceptance line, ask which
   REDUCTIONS the proposed change reorders. Zero ⟹ `array_equal` is honest.
   Any ⟹ the line is impossible and must be re-scoped, not loosened.**
2. **"ZERO diff in any operator file."** Not a pytest gate — a permanent test
   asserting a property of a diff that no longer exists is signature-tautological
   (`vv` Mode 8 class 3). Split into (a) a one-shot commit-scoped
   `git diff --name-only` check quoted in the commit message, and (b) a
   PERMANENT **AST vocabulary gate** forbidding the book-keeping the carve
   relocates (`axis=<int>` signature defaults, `.to_local(`, bare
   `n_restricted`/`n_total`) inside operator modules — with a negative leg
   feeding the walker a fixture string that MUST be flagged, else it is a grep
   that passes because it matched nothing.

### g. The decay list has THREE flavours, not one (the G5 lesson, sharpened)

G5 taught: a type-collapse makes class-level gates DECAY to tautologies, and the
re-pose must land in the same commit. An **optional→mandatory** carve is sharper,
because it changes what is CONSTRUCTIBLE:

- **DIE** — the gate can no longer build its subject (`test_space_anonymous_by_default`
  asserting `M.domain is None`; `test_space_anonymous_leaves_refuse` constructing
  three unbound operators). It reds; it must be deleted or re-posed.
- **DECAY** — still constructs, still passes, assertion is now a theorem of the
  type system (`test_leaf_declares_both_function_spaces`). **Invisible in a run.**
  Re-pose onto what stays contingent: not "a space is bound" but "the space is the
  RIGHT one" — `op.domain is mesh.full_field_space`, **identity not equality**,
  since `FunctionSpace.__eq__` is `(name, shape)` and accepts a look-alike.
- **INVERT** — the gate now pins the degradation as the contract. Canonical
  instance: `test_self_adjoint_M_H_equals_M`, whose docstring literally says
  *"the metric-blind Euclidean `.H` (domain None) reduces to the representation
  transpose"* and asserts `array_equal(M.H.apply(ψ), forward)`. Both halves break:
  the docstring's REASON goes present-tense-false, and the assertion itself reds —
  **measured, binding a metric changes a self-adjoint diagonal's `.H` by 2 nulp**
  (`G⁻¹(f·(G x))` vs `f·x` — `(g·x)/g` is not an IEEE identity).
  **Corollary: `assert_array_equal` on ANY `.H` result of a newly-bound operator
  breaks at ~2 nulp. Grep them all before the carve.**

### h. The strict-xfail set was ALREADY the todo list — and it must be read by ARM

`tests/sn/architecture/test_monomorphic_leaves.py` already carried G6.4's
acceptance gates as `strict=True` xfails whose reasons say *"WHEN THIS XPASSES:
P1 has landed — delete this marker."* Measured `105 passed, 21 xfailed` in 1.31 s,
of which **exactly 12** flip (4 `test_model_generic_leaf_declares_a_space` + 5
`test_leaf_space_annotation_is_not_optional` + 3
`test_leaf_without_a_space_refuses_construction`). The other **8** (`_R6_XFAIL`,
the carrier-guard non-uniformity) are a DIFFERENT campaign phase and must NOT be
deleted. **Deleting markers because "the xfails flipped" is exactly the Mode-8
class-4 misattribution.** Read the xfail set BY ARM, never by file.

### i. The scope alarm a decay-list sweep surfaces for free

The committed xfail `test_model_generic_leaf_declares_a_space` mirrors
`orpheus/homogeneous/solver.py` line for line and MEASURES `C.domain is None` /
`F.domain is None` there — a PRODUCTION path. ⟹ **a "boundary-tier" carve that
retires the `None` default on the BASE class silently reaches the homogeneous
solver, 12 test-local `LinearOperator` subclasses across 6 files, and a
`_multiplier` helper feeding 21 call sites.** The decay list is ~20 tests if the
mandate is boundary-scoped and ~150 if it is tree-wide, and *the plan must choose
before the first line of code*, because everything downstream is sized by it.

**Generalisable: when a campaign step says "make X mandatory", locate WHERE the
optional default is declared. If it is on a shared base, the step is not scoped
to the campaign's tier no matter what the plan's heading says.**

### j. A harness lie reproduced, in the exact anti-#17 shape

`grep "^FAILED"` on a `-q --tb=no` run returned **zero rows** on a run whose own
summary line read `33 failed` — the `FAILED` lines carry ANSI colour codes.
Identical to the 2026-08-03 battery that reported 32/32 BLIND. **Always
`sed -e 's/\x1b\[[0-9;]*m//g'` before parsing, and always cross-check the
extracted count against the summary line.** Cost here: one confusing minute,
because a positive control was already in hand to contradict it.

---

## L38 — BINDING carve, step 5 (bind the deck permutation `Γ₊→Γ₋`): the binding is ungated by construction, and it evaporates at the tensor product

**Context.** ORPHEUS `#330` G6.3 step 5, `refactor/operator-strategy-layers`,
2026-08-04. `PermutationOperator` gains optional `domain`/`codomain`;
`_specular_kernel` binds them; `is_involution` refuses across two different
spaces; `inverse()` swaps them. I was dispatched to design the gates.

### 1. The tree moved TWICE mid-dispatch, and it demoted one of my measurements

The brief said "I am writing the code directly; you design the gates" — i.e.
pre-carve. `[M]` `git status` at 22:39: `operator.py` modified (items 1/3/4
already landed). At 22:41: `realizer.py` modified (item 2 landed). The plan
went from PRE-carve to POST-carve while I was reading context.

The damage: my first bit-identity probe compared the production kernel against a
hand-bound copy of itself. Between probe run 1 (22:39, production UNBOUND) and
run 2 (22:41, production BOUND) the implementer shipped the binding, so both
sides of my `array_equal` were bound and the leg was true for free — the
lesson-#4 demotion arriving through the clock rather than through a rewire. The
tell was a *sign inversion in the data*: run 1 reported
`is_involution=True` on `gauss_legendre(4)`, run 2 reported `False` for the same
fixture. **A measurement that contradicts your own earlier measurement of the
same quantity is a tree-moved signal, not a flaky test.**

Defences adopted: build controls EXPLICITLY (`np.asarray(P.perm).copy()`),
assert `control.domain is None` INSIDE the gate, re-`git diff` before every
claim, and read "I am about to write X" as "X may already be on disk".

### 2. The binding is gated by NOTHING — measured, not argued

In-process pytest plugin rebinding `realizer._specular_kernel`; the plugin
raises if the rebind does not take and prints its entry count in the terminal
summary. Baseline `tests/geometry tests/sn/operators -m "not slow"` =
`3 failed / 1668 passed` in 24 s, `_specular_kernel entered 1252 times`.

| mutation | new reds |
|---|---|
| **POSITIVE CONTROL** — return the IDENTITY permutation, binding intact | **+23** (26 failed) |
| drop the binding (`domain=None, codomain=None`) | **0** |
| ⭐ **SWAP** the binding (`domain=Γ₋, codomain=Γ₊`) | **0** |
| bind both ends to `Γ₊` | **0** |

The control is inside the algebraic class (the identity IS a permutation:
linear, orthogonal, right shape, right spaces), so anti-#18 holds and its 23
reds mean "wrong bijection", not "stopped being an operator".

⭐ **The SWAP is the mutation to design for.** It is a one-word implementer slip;
it survives the extent guard because `|Γ₊| == |Γ₋|` on every shipped quadrature;
it leaves `is_involution` correctly `False` (the two spaces still differ), so
the refusal gate cannot see it; it changes no arithmetic. The ONLY catcher that
can exist is an `is`-identity row asserting WHICH space is WHICH end. And it is
not benign: once the composition is routed through `@`, a swapped binding makes
the LEGAL composition raise and the ILLEGAL one pass.

### 3. The binding evaporates at the tensor product — the step's own goal unmet

`[M]` through `SNBoundaryRealizer().realize(...)`:

| law | returned type | `.domain` |
|---|---|---|
| reflective α=1 | `TensorProductOperator` | **None** |
| reflective α=0.5 | `ScaledOperator` (forwards a `None`) | **None** |
| white | `TensorProductOperator` | **None** |
| — its inner `OperatorProduct` (the committed step-3 chain) | `OperatorProduct` | `Γ₊(xmin)` ✅ |

`OperatorProduct` derives its spaces from its factors; `TensorProductOperator`
does not. So a leaf bound `Γ₊→Γ₋` reaches the realizer's OUTPUT as `None`, and
the campaign's next step ("route `_reflect_trace` through `@` so the
composability check FIRES") **cannot fire** — `[M]` one `None` short-circuits
the check (`bound @ unbound` is ACCEPTED). The already-COMMITTED step 3 has the
identical hole.

⟹ **Measure the binding at the tier the CONSUMER sees.** Rule: before crediting
"the check now fires", compose the object production HANDS OUT, not the leaf you
just bound. And when a sibling step landed the same pattern earlier, check it too
— then ship the gap as a `strict=True` xfail naming the later step, rather than
widening scope.

### 4. `is_involution` — the fixture-dependent attribute, and its four legs

One physical law (mirror about `x` on `xmin`), narrowed to `Γ₊`-local indices:

| quadrature | raw `perm[perm]==arange` | UNBOUND | SAME-SPACE | HALF-BOUND | bound `Γ₊→Γ₋` |
|---|---|---|---|---|---|
| gl(4) / gl(8) / product(4,4) / ls(6) | True | True | True | True | **False** |
| lebedev(17) | **False** | False | False | False | **False** |

Three gate-design consequences:

* **A refusal gate written only on `lebedev(17)` is green on a NO-OP** — its
  `False` comes from the index test, not the space test. Every refusal row needs
  an ACTIVATION GUARD asserting `perm[perm] == arange` is True, and lebedev
  belongs in a separately-named control row.
* **Two positive legs, not one, and they are not interchangeable.** UNBOUND
  (`None`/`None`) and SAME-SPACE (`Γ₊→Γ₊`) both report True; only SAME-SPACE
  proves the clause discriminates on *space inequality* rather than on *presence
  of binding*. An implementation spelled `not (domain is not None)` passes
  UNBOUND and fails SAME-SPACE.
* **HALF-BOUND reports True** (the conjunction needs both ends) — a deliberate
  asymmetry pinned by nothing; it is the sole catcher for an `and`→`or` slip.

### 5. `.H` is metric-blind, and the criterion is EXACT (not a tolerance question)

`[M]` `G_{Γ₋} == G_{Γ₊} ∘ π` **bit-exactly on all five quadratures** — a specular
mirror preserves `|Ω·n|·w_n`, so `G₊⁻¹PᵀG₋ = Pᵀ` in exact arithmetic. This is
`vv` Mode 12's own named example. `.H` vs `apply_transpose`: `array_equal` True
on gl(4) and product(4,4), **1 nulp** on gl(8), ls(6), lebedev(17) — the
`(g·y)/g` round-trip, and note the 0-vs-1 split lands on DIFFERENT fixtures than
step 1's did. Fixtures agree by luck; `assert_array_equal` is a false red on 3
of 5. Tolerance 2 nulp from round-trip DEPTH.

⭐ **Pin the blindness CRITERION as its own gate**, not just the near-equality:
`wm == wp[perm]` is a statement about the physics, it explains why the `.H` row
is allowed to be approximate, and it reds FIRST if a future quadrature breaks
it. Negative legs, measured and quadrature-independent: scaling `G_{Γ₊}` by 3
moves `.H` by exactly `2/3`; scaling `G_{Γ₋}` by 3 moves it by exactly `2`.
Assert the NUMBER, not "it moved".

### 6. Cost — the reachable-subset rule again

The brief's cost ladder said `+ tests/numerics` ≈5 m 45 s. `[M]` the positive
control reddens **zero** files under `tests/numerics`, and no `tests/numerics`
module imports `sn.boundary.realizer`. The whole realizer-side battery lives in
a **24 s** slice. Budgeting off the directory rather than the reachable subset
would have made a 7-mutation battery 40 minutes instead of 3.

⚠ **The ANSI trap fired again**: `grep "^FAILED"` returned nothing against a run
that had plainly failed. `--color=no` when parsing. Third instance.

### 7. A retirement question the plan surfaced but did not decide

`is_involution` has **zero production consumers**, its one would-be consumer now
reports `False` by design, and the new docstring concedes *"nothing in the tree
does"*. Zero consumers + a self-conceding docstring is `coding-elegance`
Pattern-6's trim signal verbatim. Raised as an explicit user ruling (refine vs
retire) rather than decided in a test plan — because if it is retired, three of
the ten gates collapse into "ask the algebra" (`P @ P` raises), which is the
better claim anyway.

---

## L39 — P3 affine-boundary-source carve: the CONSUMPTION-TIER readout, and a "prophylactic" carve that was a live bug

Dispatch 2026-08-05, `refactor/operator-strategy-layers` @ `ef4c3537`. Brief:
plan the verification for retiring `IncomingSourceOperator` and collapsing
`realize(PrescribedInflow)` onto the zero morphism. Brief already carried a
measured double-delivery finding; my job was the gate design.

### 1. ⭐⭐ The keystone gate was NOT one of the three candidates offered

The brief offered (a) `B(0)==0` linearity, (b) two-user-path flux agreement,
(c) superposition in `q`. The winner was a FOURTH: **evaluate the boundary law
on the converged answer.** For prescribed inflow `L = 0`, so
`γ₋ψ|_f == q_f` — read from `sol.angular_flux.boundary.face_view(f)[inflow]`.

`[M]` het-2G scattering slab, identical composite bulk on every row:

| config | SI `γ₋(xmin)` | Krylov |
|---|---|---|
| declared, HEAD (both channels) | `5.000000000000` | RAISES |
| declared, pre-P2′ (operator channel only) | `2.500000000000` | RAISES |
| all-vacuum + composite `q_∂` (source channel only) | `2.500000000000` | `2.500000000000` |
| vacuum control | `0.000000000000` | `0.000000000000` |

**Three configurations, three bit-exactly distinguishable readings, ONE
assertion, NO reference solver.** Double = `2q`, loss = `0`, wrong `L` =
`q + Lγ₊ψ`. Closed-form pillar (the BC *is* the definition), discretization-,
mesh-, quadrature- and materials-independent. `assert_array_equal` is EARNED
(`5.0` is `2.5+2.5`, exact; `2.5` is a copy).

⭐ **Generalisable rule: when a carve is about "is this term delivered, and how
many times?", look for a tier where the governing EQUATION can be evaluated on
the converged answer.** A trace/boundary/interface DOF is usually such a tier,
because the posed system's rows there ARE the constitutive law. That beats a
snapshot (self-referential), beats a flux comparison (needs tolerance + a
non-vacuity leg), and beats linearity (blind to the source channel).

### 2. ⛔ Candidate (c) — superposition in `q` — is a Mode-12 NON-CATCHER

The double delivery IS `q → 2q`, **which is still exactly affine in `q`**.
`φ(a) = φ(0) + a·s` holds for every `s`, including `s = 2 s_true`. The scale
factor is in the measured functional's invariance group ⟹ designed-green on the
whole defect class at every tolerance and refinement. Only the *coefficient*
pinned against an independent reference has teeth. **A "linearity in the
parameter" gate can never detect a wrong CONSTANT multiplying that parameter** —
check this before adopting any superposition row.

### 3. ⭐⭐ The "prophylactic architecture" framing was FALSE — measure the fence

The plan's §1 asserted TWO fences and concluded the carve was prophylactic. Both
failed: (i) the `block_role=None` stamp fences nothing because
`SNBoundaryOperator._face_laws` collects EVERY face's law with no role filter
(`|B(0)|_inf = 2.5`); (ii) the consequence is a **hard raise** on Krylov —
`ConvergenceCertificateError: ‖Aψ−q‖/‖q‖ = 1.718` — because an affine
`A(x) = A_lin(x) − c` breaks GMRES's Arnoldi relation.

⭐ And the raise fired on **BOTH sides of P2′**, so the Krylov path was already
dead pre-P2′: P2′ was a double-delivery regression on **SI only**. **A "not a
live bug, fenced twice" claim is a HYPOTHESIS — apply `B(0)` (or the fence's own
predicate) to the object the consumer actually holds.** A stamp only fences if
somebody reads it; grep the consumer for the filter before believing the fence.

### 4. ⚠ anti-#18 in its sharpest form: the OBVIOUS mutation is out-of-class

The natural mutation ("put the affine operator back") **breaks linearity**, so
every solve-level red comes from the operator ceasing to be linear, not from the
delivery count. Reds under it must NOT be credited as single-delivery coverage.

⭐ **The in-class mutation for a delivery-COUNT claim is `q + q` in the SOURCE
channel** — everything stays linear, so every red is attributable to the count
alone. Test: *a gate that reds under the affine reinstatement but NOT under
`q + q` is not a single-delivery catcher.* Second in-class mutation, for the
carve's own content: `L := IdentityOperator` (perfectly linear) — gate A
(`B(0)==0`) is BLIND to it; only the trace gate sees it. That asymmetry is why
one gate cannot replace the other.

### 5. ⭐ The tree was being written UNDER me — three times, and each time the
### reconciliation outvalued the forecast

`realizer.py` had lost the `IncomingSourceOperator` import while the call site
still referenced it (`NameError`) — so I measured "before" from
`git worktree add <tmp> <HEAD-hash>` (NEVER `git checkout`; the tree carries
irrecoverable uncommitted state). Then, mid-plan: four test files migrated, a
brand-new `tests/numerics/test_zero_operator_spaces.py` appeared, and the 8-red
debt I had measured went to **`555 passed`** before I finished. So §2 became an
AUDIT of what landed and §8 the residual-gap list.

⭐ **Re-measure the red count at the END of a planning dispatch, not just the
start.** A migration table whose status column is stale by minutes is worse than
no table — add an explicit "READ §8 FIRST" banner rather than silently shipping
it. And the deliverable that survived was **the gap list**: measured absent by
grep at the end — Gate D, Gate A, the two-channel gate, the Krylov row, the
production SWAP row, the `checked_space_extent` negative, the `is_adjointable`
widening, the un-mutation-verified `catches` re-home, `DEAD TARGETS 2`, a
label-less `:ref:`, and the whole battery.

### 6. ⭐ A measured number recorded in a COMMENT is not a gate

`|B(0)| = 2.5` — the campaign's central measurement — appears in the tree ONLY
as prose at `tests/sn/operators/test_operator_block_role.py:203`. Excellent
documentation; zero teeth. **When auditing a landed carve, grep for the
measurement's NUMBER and ask whether any `assert` consumes it.** The
"honest metadata that gates nothing" shape.

### 7. Did the retired oracle's independence survive? Audit, don't assume

`test_boundary_source_from_specs.py` used `IncomingSourceOperator.apply` as its
shape oracle. **Verdict: no demotion** — that operator's entire body was
`source.evaluate((n_inflow,) + psi_out.shape[1:])`, a ONE-LINE ADAPTER, so
inlining it to `spec.evaluate((|Γ₋|, ng))` preserves the independence, which
always lived in *who produces each side* (the bridge vs the user's own spec
object), never in the adapter.

⭐ **The discriminator for "did this rewire demote the gate?": is the retired
symbol a SOURCE of the expected value, or a FORWARDER of it?** A forwarder can
be inlined for free. Corollary trap recorded as a refuted candidate: making the
new oracle dimension-generic by reading `slot_shape[1:]` — the SUT's own
derivation — WOULD be the demotion.

### 8. The `|Γ₊| ≠ |Γ₋|` Mode-12 fixture migrated STRICTLY STRONGER

`TestIncomingSourceOperator` was the only place unequal half-trace extents were
reachable (equal on every production face), which is what discriminates "fills
the codomain" from "echoes the input". `ZeroOperator` is equally hand-buildable
AND adds a transpose direction and two space identities the apply-only affine
operator never had. ⭐ **When a Mode-12 fixture's subject is retired, check
whether the successor is also hand-constructible — if it is, the claim migrates
and usually GAINS legs.** Two production docstrings had independently written
the same argument (`_narrowed_zero_operator`: *"wrong in principle and merely
lucky in practice"*), which is how I knew the claim was load-bearing.

### 9. The production SWAP hole that a same-object gate cannot see

The landed `test_prescribed_inflow_realizes_the_same_object_vacuum_does` asserts
`prescribed.domain is vacuum.domain` — that the two AGREE, never that **either
is the RIGHT space**. Both bound swapped ⟹ green, and `|Γ₊| == |Γ₋|` on every
production face ⟹ `checked_space_extent` passes too. `L37` measured that three
wrong bindings each produced ZERO new reds over 1668 tests. ⭐ **An
agreement-between-two-siblings row is NOT an is-it-correct row** — a collapse
carve naturally produces the former and it feels like coverage.

### 10. Doc-debt instrumentation: two instruments, and a hole between them

`[M]` `tools/check_docstring_xrefs.py orpheus tests docs` → `DEAD TARGETS 2
across 9 sites`, with file:line and role — far better than a hand grep. BUT it
gates only **fully-qualified Python-domain** targets (7073 of 14541 roles), and
`sphinx -W`/`-n` cannot see a docstring in an un-`automodule`'d module (verified:
`orpheus.sn.boundary.*` and `orpheus.numerics.operator` have none). ⟹ a
`:ref:`bc-affine-source-channel`` cited twice by the new code with **no target
anywhere** is caught by NEITHER. ⭐ **After a carve adds a `:ref:` from a
docstring, grep the label in `docs/` by hand — no instrument reports it.** And
raw grep found 13 files vs the checker's 9 sites: prose mentions are invisible to
the checker, so run BOTH and triage BY TENSE (past-tense provenance STAYS).

### 11. Fixture facts (add to the config-blindness inventory)

* `_build_fixed_source_rhs`'s two arms take **different bulk sources** —
  per-ordinate `(N, ng, *spatial)` array vs an `AngularSourceSink`. Feeding one
  to each leg of a two-channel comparison gave `φ[0] = 3.083` vs `2.480`: a bulk
  difference misread as a channel difference. **Use ONE spelling on every leg.**
* Two channels reaching the same fixed point are **NOT bit-identical**:
  `|φ_D0 − φ_C|_inf = 1.998e-13` at `inner_tol = 1e-13` (rel ≈ 0.3 × tol),
  `array_equal = False`. Gate at `SAFETY(10) × inner_tol`.
* `tests/sn/solve` naming trap: `test_affine_carve_bit_identity.py` is **#208's**
  `FluxDisplacement` carve, NOT the affine BOUNDARY carve — all-vacuum configs,
  `sha256`-frozen, and 3 of its params are pre-existing reds. `[M]` its signature
  is unchanged by P3, which is the useful fact: **a `sha256` gate that MOVES
  during a carve that should not touch its path is a real signal.**
* Costs `[M]`: the P3 blast-radius slice **24.5 s** (555 passed);
  `tests/sn/solve + tests/numerics` **7 m 27 s** (2016 passed, the 3 known reds).
