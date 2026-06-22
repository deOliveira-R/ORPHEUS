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
of that standard.

---

## L1 — Homogeneous / 1-group / flat-flux tests are blind to the hardest math

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
twice the plan claimed "bit-identical" / "clean O(h²)" / "improves on
flat" and all three were measurably false:

- "iso-MMS bit-identical" was exact-arithmetic-only (≤1 ULP in IEEE).
- "aniso clean O(h²) at S16" was a #229-interpolation floor, untouched by
  the carve (needs S32).
- "improves-on-flat at the boundary" was unachievable (the slope is
  sub-floor; the correctly-consumed slope made the value slightly WORSE).

How to apply: measure the premise (temp-edit + revert, matched-quadrature
ladders, probe the regime the term dominates) rather than encoding the
plan's hope into a gate. A gate built on a false premise either calcifies
a limitation or falsely reds a correct term. State the refutation in the
plan so the implementer ships the achievable gate, not the hoped-for one.
