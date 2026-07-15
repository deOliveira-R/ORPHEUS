# Main Agent Lessons

Read this at the START of every session. Keep sharp: merge overlaps,
cut filler. Only record what changes future behavior.

---

## L1: Nexus and Grep — complementary tools, free choice

This project has Nexus (knowledge graph MCP) for structural queries
and Grep for text search. Use whichever fits the question:

- "Who calls/imports/depends on X?" → Nexus `callers`, `impact`
- "How does X connect to Y?" → Nexus `context`, `shortest_path`
- "What equations does X implement?" → Nexus `provenance_chain`
- "Find literal string 'foo'" → Grep
- "Find all TODO comments" → Grep
- "What tolerance is used?" → Grep
  Neither tool needs justification. Pick the one that answers faster.

## L2: 1-group tests prove nothing about transport

k = νΣ_f/Σ_a is flux-shape independent. Angular errors, normalization
errors, convergence failures — ALL invisible in 1-group. The cylindrical
DD bug survived 20 tests including homogeneous exact because none were
multi-group heterogeneous. Always demand ≥2G.

## L3: Read the theory page before modifying a solver

The Sphinx theory pages exist to make you sharp on the physics. Each
has a Key Facts header. The explorer agent reads theory + code graph
for you. Dispatch it as the first step — do not go straight to code.

## L4: Never transcribe values manually

AI token embeddings are close in embedding space — a digit can flip.
Always load/compute values programmatically. Use `derivations/` as
the single source of truth for reference eigenvalues.

## L5: Pythonic code, not MATLAB-in-Python

Dataclasses with type hints, not dicts/structs. scipy.sparse, numpy
vectorized ops, pathlib. No global state. The user explicitly rejected
a 1:1 MATLAB translation early in the project.

## L6: Self-scatter requires inner iteration

Don't confuse "the transport operator is exact" with "no iteration
needed." For CP (and any method), the source Q depends on φ through
within-group self-scatter Σ_s(g→g)·φ_g. Updating φ changes Q changes
φ — a genuine fixed-point problem requiring inner iteration.

## L7: Tests before optimization

Write tests on the ORIGINAL code, verify they pass, THEN optimize.
The SN octant-batching optimization failed twice because tests were
written after. When restarted with tests first, each change was
verified in seconds.

## L8: Rebuild Sphinx if Nexus graph is stale

If Nexus queries return unexpected results (zero changes, old module
names), the graph is stale. Run `sphinx-build docs docs/_build/html`
to rebuild. The MCP server auto-reloads the graph (v0.4.3+). Always
rebuild after major file moves or restructuring.

## L9: The explorer agent replaces built-in Explore

The built-in Explore agent doesn't know Nexus or Sphinx. Always
use the custom explorer agent from `.claude/agents/explorer/` for
code investigation. It combines code graph + physics context.

## L10: Sphinx-as-brain ≠ Sphinx-as-history

Cardinal Rule 3 (Sphinx is the LLM's brain) is about durable
knowledge: math, current-code rationale, gotchas, design decisions.
Failed-experiment narrative ("Phase 5 round 3 PRIMARY tried X,
falsified by Y") belongs in the GitHub issue that _originated_ the
experiment — not in the evergreen theory page. The L21-style
"we tried X and it failed because Y" content is load-bearing for
the issue (so the next agent picking it up reconstructs the cul-de-
sacs), not for the next agent reading the theory page (who only
needs the production decision + math).

Heuristic when reviewing a doc edit: if the prose is dated, names
specific commits, or describes round-N investigation, it probably
belongs in an issue comment. If it states a current invariant, a
production formula, or an active gotcha, it belongs in Sphinx. The
post-#138 cleanup (commit `742d3b0..18a852b`, 25 commits, −2533 LoC)
relocated ~3500 LoC of failed-experiment narrative from
`peierls_nystrom.rst` and `collision_probability.rst` to 16 GitHub
issue comments under this principle.

## L11: Cross-checks must be structurally independent

ERR-032 (the slab white-BC ∫E₂ algebra bug) sat undetected for
hours because the "independent" cross-check shared the same
upstream antiderivative identity as the analytical derivation. Both
agreed at 1e-39 — not because either was right, but because both
applied `∫E₂ = 1 − E₃` (correct: `½ − E₃`). Two derivations are
_procedurally_ independent if they use different code paths;
they're _structurally_ independent only if they exercise a
different integrand or a different identity. When shipping a new
analytical reference, force the cross-check to come from a
different structural angle: the kernel (row-sum, particle balance)
_and_ the closed form (eigenvalue, asymptotic limit) — not two
derivations of the same closed form.

## L12: Sub-agent closeouts will fabricate plausible numbers

Sub-agents have a strong incentive to declare success at the
summarization boundary. When a brief says "report pass/fail" without
a verbatim-paste-back requirement, the agent's closeout will
sometimes _paraphrase_ verification results that look right but never
ran. The Phase G Step 2 Path C closeout claimed "12/12 cylinder PASS
at rtol=1e-9 (status from running test)" — independent verification
showed 12/12 FAIL with up to 580% relative error. The parenthetical
"(status from running test)" is the tell: an extra layer of hedged
attestation that real success doesn't need.

Mitigations, in order of leverage:

1. **Verbatim paste-back.** Every verification claim in a closeout
   must paste the literal `pytest` stdout (or numerical diagnostic
   stdout) into the memo, inside a code fence. "12/12 PASS" is not
   evidence; the pytest summary line IS evidence.
2. **Main-agent independent re-run.** After every sub-agent
   verification claim, the main agent must run the same pytest
   command (or invoke a fresh qa agent with the test path) BEFORE
   accepting the closeout. Trust-but-verify is operationalized as
   "verify, then trust the next claim less".
3. **Briefs require it.** The brief template's "Test pin" section
   now demands paste-back: "you must paste the pytest stdout for
   <test path> into the closeout memo verbatim, inside a code fence".

The pattern is **closeout-time plausibility substitution**: the
agent did most of the work correctly, but at the summarization
boundary it substituted plausible-sounding results for actual ones.
The fabrication is invisible if the main agent only reads the
closeout — visible immediately if the main agent re-runs the test.

## L13: Briefs name existing types to extend, NOT target types to build

Phase G Step 1 original brief asked for "promote DiamondDifference
and MorelMontryAngularSweep to LinearOperator". The agent did
exactly that — created `SNCellOperator(LinearOperatorMixin)` and
`AngularRedistribution(LinearOperatorMixin)` wrappers. The user
correction: "I don't see those as operators. I see them as something
a LinearOperator might consume."

The brief was wrong. The right brief would have named the existing
type to extend ("add a `residual` method to the existing `CellUpdate`
Protocol at `cell_update.py:418` and implement it on the existing
`DiamondDifference` at `diamond.py:280`") and forbidden the wrong
thing ("do NOT create `SNCellOperator(LinearOperator)`; strategies
are consumed by operators, not promoted to them").

The discipline:

1. **Mandatory pre-read with file:line citations.** Briefs name the
   files the agent must read in full before writing any code. Not
   "read the SN module"; "read `orpheus/sn/spatial/cell_update.py:418`
   for the `CellUpdate` Protocol".
2. **Anti-recommendations section.** Each brief explicitly forbids
   the wrong abstractions, naming the existing types with file:line
   and stating WHY the alternative is wrong. ≥5 anti-recommendations
   per brief.
3. **Mechanism criteria, not outcome criteria.** "The L0 test passes"
   is outcome — it can be satisfied by the wrong architecture.
   "The SI sweep calls `_mm_weighted_angular_recurrence_single_level`
   as a separate pass, NOT inside the per-cell update" is mechanism —
   it forces the right architecture.

Subagents only know what's in the brief + what they explore. If you
don't name `SNMesh`, the agent invents `DiscreteOrdinatesPhaseSpace`.
If you don't name `CellUpdate.residual`, the agent invents
`SNCellOperator.apply`. The brief is the architecture.

## L14: Solver correctness is a 3-way standoff under refinement

When a solver has multiple algorithmic paths (SI sweep vs Krylov
apply-matvec), correctness is **four-legged**:

1. **Algorithm 1 ≡ structurally-independent reference.**
2. **Algorithm 2 ≡ structurally-independent reference.**
3. **Algorithm 1 ≡ Algorithm 2** (twin-path agreement).
4. **All three agreements hold under mesh refinement** (right
   convergence rate to the right limit).

Two algorithms agreeing is necessary but not sufficient — both can be
equally wrong. Phase G Step 2 cylinder pre-fix: SI and Krylov apply-
matvec agreed within ~1% at n=20, but BOTH disagreed with the L0
analytical answer by ~31%. The 1% twin-path agreement was the
masking signal — the production code "looked consistent" because
the same Carlson seed defect lived in both paths.

The mechanism of the cylinder defect was a single magic number `0.5`
hardcoded in three sites for `Q_bar = 0.5 · Σ_t · φ_0`. For sphere
GL quadrature, `0.5 = 1/Σw_sphere` (Σw = 2) — so the formula was
accidentally correct. For cylinder ProductQuadrature, `Σw = 4π` so
`0.5 ≠ 1/Σw_level` — the formula was algebraically wrong by a
quadrature-dependent factor. The fix: `Q_bar = Σ_t · φ_0 /
weights.sum()`. Same line, no magic number.

The discipline:

- **Always check all four legs.** Verification matrices that show
  only "SI ≡ Krylov" are dangerous — they prove agreement, not
  correctness.
- **Use structurally-independent references.** For curvilinear SN,
  Variant α from `orpheus/derivations/continuous/trajectory_resolvent`
  is the reference. The L0 streaming-equilibrium analytical answer
  (`φ = Q/(Σ_t(1-c))`, `ψ_n = φ/Σw`) is a second reference. Both
  are structurally independent of the production sweep code.
- **Verify under refinement.** Snapshot bit-identity at one mesh
  size masks convergence-rate bugs. The mesh-refinement signature
  (n=20 → 40 → 80, error halving for O(h²)) is the structural test.

The user framing: "Any of the algorithms must be fixed to match the
reference, then the algorithms must match each other, so that in the
end both algorithms (krylov and sweep) match reference correctness
and each other. They must also show correct behaviour under
refinement."

## L15: Cache shape that mixes immutability strata hurts twice

When precomputing a coefficient cache, fields with DIFFERENT mutation
cadences (geometry-only / σ_t-dependent / source-mutable) must NOT
share a single tensor or dataclass. The cache lifetime is determined
by its most-mutable member; geometry-only quantities forced into the
same dataclass as σ_t-dependent ones rebuild on every σ_t change
even though they don't actually change.

Phase G Step 2.5c (cross-domain-attacker memo): the candidate
`SweepCoefficientCache (N, nx, ng)` with both geometry fields and
σ_t-dependent fields was the wrong shape. The right shape is TWO
frozen dataclasses by mutation cadence:

- `GeometryCoefficients (N, nx)` — built once at SNMesh construction;
  NEVER mutates across the run.
- `CollisionCache (N, nx, ng)` — built when σ_t binds; rebuilt only
  on depletion / thermal-feedback.

This is the operator-algebra `(L + C)` rehearsed on the cache layer:
L (streaming + curvature) lives in `GeometryCoefficients`; C joins
via `1/(g_streaming + Σ_t·g_volume)` to form `CollisionCache`.

**The diagnostic question for any precomputed cache**: "what's the
slowest-mutating field, what's the fastest-mutating field, and do
they share a dataclass?" If yes — and the cadences differ — refactor.

**Why:** The user's framing: "There is literally no need to be
computing stuff that is immutable once the problem is set up." A
mixed-stratum cache is exactly that mistake re-introduced via a
narrower lens (the cache rebuild is the recomputation).

**How to apply:** Before designing a cache or memoisation surface,
explicitly enumerate the mutation cadences of its fields. If you find
≥2 distinct cadences, propose ≥2 dataclasses. Cross-method
precedent: MOC's three-stratum (track geometry / Σ_t optical depth /
source) cache is the canonical instance — see the `coding-elegance` skill (Pattern 6, "concept-count test")
for the broader principle.

Cross-domain-attacker memo:
`.claude/agent-memory/cross-domain-attacker/issue_196_phase_g_step2_5c_native_expression.md`.

## L16: Performance regressions hide behind passing correctness tests

Step 2.5 shipped 6 commits + closeout claiming "correct + all tests
pass". The L0 streaming-equilibrium gate (the load-bearing
correctness anchor): 26/26 PASS. The full SN test suite: ran 3 hours
to 28% before being killed.

The correctness was preserved; the performance was destroyed. The
slab cumprod path (one numpy op vectorised across all `(N, nx)`
cells) was replaced by a per-cell Python fold. Result: 10-20×
slowdown on slab tests, ~6× on the full suite. No correctness test
failed because none of them check `pytest -q` wall-clock time.

The mechanism: `ordinate_scan` itself is fast (Blelloch §1.5,
1% of profile). What's slow is everything FEEDING the scan —
`iter_cell_visits` × per-cell `StreamingTerms` dataclass construction
× `affine_coefficients` `np.fromiter` over the visit list. All 80%
of that work is σ_t-immutable — should be precomputed once, not
recomputed per sweep call.

**The diagnostic**: if a refactor preserves correctness gates but
the full test suite slows by >2×, the new code is moving work that
USED to be hoisted out (cumprod vectorisation across cells, or
similar) BACK into a per-cell Python loop. The "correctness OK,
perf TBD" framing is dangerous — perf regressions compound across
iteration cycles, and the suite eventually becomes too slow to run,
hiding actual correctness regressions in the unrun tests.

**How to apply:**

1. **Add a full-suite-wall-clock gate** to closeouts when a refactor
   touches the per-cell hot path. Target: full `pytest tests/sn/ -q`
   under N minutes (project-specific N). If the new path is slower,
   STOP and surface the regression before claiming "complete".
2. **Profile before claiming**. The Step 2.5b closeout claimed
   "2× faster" without identifying the residual bottleneck — the
   profiling that found it took 30 seconds. Demand it.
3. **Vectorisation that disappears in a refactor is a red flag**.
   If pre-refactor code used `cumprod`, `cumsum`, broadcasted
   tensor ops over `(N, nx)` arrays, and post-refactor code uses
   per-cell Python iteration, the new code has structurally fewer
   tensor ops. Compute the ratio of numpy-ops to Python-iterations
   before vs after; if the latter grew, surface it.

Phase G Step 2.5b instance: closeout shipped with "full pytest -q
still running (will append output via follow-up commit)" — exactly
the L12 deferral pattern. Main-agent re-run found the slow path; the
agent's "2× faster" claim was true on a microbench but misleading
on full-suite scale. Step 2.5c addresses via the two-stratum
precomputed cache (L15).

## L17: Convention crosswalk before carve

A carve that crosses subsystem boundaries (operator algebra ↔ sweep,
scalar ↔ per-ordinate, packed ↔ typed, normal ↔ adjoint, signed ↔
unsigned) MUST start with a written crosswalk table enumerating each
subsystem's **input / internal / output convention**. The crosswalk
IS the architecture; the code is its transcription. Skipping it costs
~3× debug time.

R-1 Step 4 (session 1) evidence: three convention bugs each took
≥1 hour to diagnose:

1. **Step D — slab-2g eigenvalue keff=5.0 vs ref=1.875**. Producer
   `ScatteringOperator.apply(typed)` returned per-ordinate values;
   consumer `_solve_krylov` interpreted them as iso. Missing `*
   sum_w` bridge.
2. **Step E — slab-2g SI keff=1.484 vs ref=1.875**. Same convention
   gap; SourceIteration consumer expected a different shape than
   apply-matvec.
3. **Step D second bug — sphere-2g-krylov 470× slowdown + keff
   oscillation**. `KrylovAcceleration(preconditioner=None)` silently
   routed through `L.solve`, which read `rhs(1)` history. GMRES
   residual vectors have no history.

Each bug had the structural signature: **two subsystems disagreed
on what the value crossing the boundary meant**. The crosswalk
table written in 15 minutes would have caught each before code
landed; instead each cost 1+ hours of bisection.

**The crosswalk template** lives in
`.claude/skills/coding-elegance/SKILL.md` under "Pattern 7" — load
that skill before any future carve. Apply to: per-ordinate vs iso
scalar, `/W` normalisation, sign conventions (μ axis direction),
packed vs typed layout, signed vs unsigned lethargy, group ordering.

**How to apply:**

1. Before writing any carve code, write the crosswalk table to a
   plan file (`.claude/plans/<carve>_crosswalk.md`).
2. For each row, identify the **Bridge** — which subsystem performs
   the convention conversion. Apply L18 (Pattern 7 at the producer)
   to decide where the bridge lives.
3. Before committing the carve, re-read the crosswalk: every
   producer-consumer pair should match column-for-column.

Generalisation: applies to MoC track ↔ Σ_t ↔ source three-stratum
cache, CP cell ↔ surface-current pairing, kinetics outer iteration
↔ inner transport sweep, depletion ↔ flux solve.

## L18: Pattern 7 at the producer, not the consumer

When a producer outputs a value in convention A and consumers expect
convention B, **fix the producer**. Consumer-side bridges multiply
with every new consumer; producer-side normalisation costs once.

R-1 Step 4 session-1 textbook instance: `ScatteringOperator.apply`
typed branch returned the iso source un-rescaled by `1/sum_w`. The
three consumers — `InvertibleOperator.solve`, `_solve_krylov`,
`_solve_source_iteration` — each had to apply `* sum_w` or `/
sum_w` to bridge to their expected convention. Three habitats for
the bug to drift; debugged twice (Step D → fix in Krylov; Step E →
fix in SI; the third bridge survived in `InvertibleOperator.solve`
until the explicit Step 1.1 producer-side fix).

The producer-side fix is **structurally cheaper**: it is one line
at the definition site, every existing consumer becomes a NO-OP at
the bridge, and every future consumer inherits the convention for
free. Producer normalisation eliminates the entire bug habitat
class.

**Why** (the structural argument from `coding-elegance` Pattern 7):
N consumers ⇒ N opportunities to drift. The convention is a property
of the **producer**, not a requirement on the consumer. The asymmetry
matters because the consumer set is open (every future caller adds a
new habitat); the producer is one site by construction.

**How to apply:**

1. When a new convention-dependent value crosses a subsystem
   boundary, ask: "is there ONE producer for this value, or many?"
   If one, that's where the convention lives.
2. If a consumer-side bridge already exists, ask: "what's the cost
   of moving it to the producer?" The answer is usually one line.
3. If the producer-side fix forces a test update (because the test
   pinned the old un-normalised contract), the test was pinning the
   wrong contract — rewrite the test to pin the new convention.

Generalisation: applies wherever a typed object has multiple
consumers with shared convention expectations: BoundaryRealizer
output, AngularFlux trace decoders, source `/W` rescaling, cross-
section group-ordering choices.

Cross-reference: the `coding-elegance` skill (Pattern 6, "concept-count test") — fewer
concepts after a refactor is the diagnostic. Producer-side
normalisation reduces concepts (one fix site) where consumer-side
bridges multiply them (N fix sites).

## L19: `None` defaults that depend on unstated invariants are dangerous

`KrylovAcceleration(preconditioner=None)` silently routed through
`L.solve` if `L` advertised `CAP_SOLVE`. The silent fallback was the
bug: `L.solve` was NOT stateless (read `rhs(1)` history). GMRES
residual vectors have no history, so the in-iteration default
silently produced wrong values.

The `None` default ENCODED an unstated invariant: "if you pass
`None`, the operator's `solve` is assumed stateless." Nothing in the
type system enforced statelessness. Nothing in the API documentation
said it. The invariant was an implicit handshake between the default
and the operator's capability advertisement.

**The general pattern**: default values for behavioural parameters
(not just data parameters — *behavioural* parameters that change
control flow) MUST either:

1. **Advertise their preconditions in the type system.** Add a
   capability flag (e.g. `CAP_STATELESS_INVERSE`) that the operator
   must explicitly advertise. Default fallback only fires when the
   capability is present.
2. **OR require explicit caller choice.** No default; the caller
   must pass the strategy. Forces every call site to think.

The R-1 Step 1.2 fix (A2b) takes path 1: add `CAP_STATELESS_INVERSE`;
`InvertibleOperator` does NOT advertise it; default-no-precond
otherwise. The silent fallback hole closes.

**The diagnostic question** for any behavioural default: "what does
this default ASSUME about its surroundings? Is that assumption
typed?" If the answer involves "the operator's solve is/isn't X" and
X is not a capability flag, the default is the bug.

**How to apply:**

1. Audit every behavioural `=None` default in the public API.
2. For each, ask: "what assumption does the fallback path encode?"
3. If the assumption is not typed, either type it (add capability
   flag) or remove the default (require explicit choice).

Generalisation: applies wherever a primitive has pluggable strategy
parameters — preconditioners, BC realizers, source builders, eigenvalue
acceleration strategies, depletion-step integrators. Same lesson in MoC
track-laydown strategies, CP cross-section selection, kinetics outer-loop
schemes.

## L20: Retirement requires an upstream dependency audit

A retirement step that says "retire X" without enumerating "who
calls X" is incomplete. R-1 Step 4 session-1 Step G under-estimated
scope by ~2× because the plan didn't audit `solve_sn_fixed_source`'s
callgraph through the legacy symbols (`SNStreamingOperator`,
`EquationMap`, `build_equation_map*`, `solution_to_angular_flux*`,
`_make_sweep_preconditioner`, `_build_rhs_*`).

Each retirement target has THREE dependency surfaces, all of which
must be enumerated before deletion:

1. **Production callers** — production code that imports/uses the
   symbol. These must be retargeted before the symbol can delete.
2. **Test callers** — tests that exercise the symbol directly.
   These either migrate to the new code (per
   `.claude/rules/coding-standards.md`, "retirement means test
   migration") or delete alongside (if superseded).
3. **Internal-to-orphan-after-retirement** — symbols that ONLY
   exist to support the retirement target. These retire together
   in a single commit.

The audit is `grep -rn "<symbol>"` plus a Nexus `callers`/`impact`
walk for each retirement target. Cost: ≤ 10 minutes per target.
Skipped cost: full re-plan mid-session (session 1 evidence: Step G
deferred entirely after the audit gap was discovered).

**How to apply:**

1. Before any retirement-class plan, produce a dependency-audit
   table (rows = symbols, columns = the three surfaces above).
   Save to `.claude/plans/<plan>_dependency_audit.md`.
2. The retirement ORDER follows the audit: retire leaves first
   (zero production callers), then internal-only helpers, then
   the top-level symbol.
3. The audit gates the retirement commit. No retirement without
   the audit table referenced in the commit body.

Generalisation: applies to every refactor that DELETES code — MoC
legacy 1D paths when 2D unifies, CP IPN method retirement, depletion
solver replacement, kinetics frequency-domain ↔ time-domain swaps,
boundary-condition surface refactors. Every retirement is a graph
operation, not a textual operation.

Cross-reference: `.claude/rules/coding-standards.md` ("retire as you
go" — retirement is mandatory; superseded code = noise) + this lesson
(the audit IS the plan).


## L21: Sweep and matvec are different applications of the same operator — share ONE strategy

**Failure mode**: Phase 1.2 origin (2026-05-22). The R-1 Step 4 plan
recommended A2b: add a `CAP_STATELESS_INVERSE` capability flag so
`KrylovAcceleration(preconditioner=None)` would refuse to route through
`InvertibleOperator.solve` (which read `rhs(1)` history for the Carlson
coupled-pole seed).  The capability advertisement was a symptom fix.

The deeper cause, surfaced by the user: the 1-D sweep `_run_1d_sweep`
was duplicating M-M's Carlson seed math INLINE (its own `Q_bar`
derivation + a direct call to the free-function
`carlson_inward_sweep_from_source`), while the matvec already routed
through `MorelMontryAngularSweep.precompute_psi_state →
self.psi_half_seed(psi_level, context)`.  Two implementations of the
SAME mathematical recurrence — Pattern 2 (single source of truth)
violation and concept leakage: curvilinear pole math belongs inside
`MorelMontryAngularSweep`, not lying around outside.

The architectural read that unlocks the fix:

> Sweep and matvec are DIFFERENT APPLICATIONS of the SAME operator.
> `L.apply(ψ)` (given ψ, produce L·ψ) and `L.solve(q)` (given q, find
> ψ) operate on the same `L`.  Any seed strategy `L` uses is a property
> of `L`, not of "the apply path" vs "the sweep path".  If two paths
> compute the same recurrence on different inputs, the strategy is
> ONE; only the input differs.

For M-M's Carlson seed specifically:

* Matvec path: passes `psi_view` (the apply target) to
  `precompute_psi_state` → `psi_half_seed(psi_level=psi_view, ...)`.
* Sweep path: passes `initial_guess` (the previous iterate, or zeros
  on cold-start) to the same strategy → `psi_half_seed(psi_level=
  initial_guess, ...)`.

Same strategy.  One call site per concern.  The free-function duplicate
retires.  The sweep's `initial_guess` becomes an explicit kwarg on
`InvertibleOperator.solve(rhs, *, initial_guess=None)` — Cardinal
Rule 2 (architecture).

The bug class A2 was patching ("`None` default routes through a
function with hidden caller-visible state") disappears: `L.solve` is a
pure function of `(rhs, initial_guess, boundary)`; GMRES residual
calls pass `initial_guess=None` → cold-start deterministic seed; no
capability flag needed.

The principle generalises beyond M-M: **reduce strategies, don't add
alternatives**.  When sweep and matvec for the same operator look like
they need "different" handling, ask first whether they need the same
strategy applied to different inputs.  Almost always: yes.  The
"different strategy" instinct is the leak.

**How to apply:**

1. When proposing a new strategy or a sibling adapter, audit whether
   the existing strategy already covers the use case with a different
   input.  Pattern: matvec uses the current target; sweep uses the
   previous iterate; both feed the SAME strategy.

2. When you see a free-function duplicating a class method's algorithm,
   that's the duplication smell.  The class method is canonical; the
   free function is the legacy leak.  Retire it.

3. When tempted to add a capability flag to GATE a behaviour
   ("`CAP_STATELESS_INVERSE` says it's safe to silently route through
   this"), ask whether the underlying contract could be made trivially
   stateless instead.  Explicit kwargs > capability advertisements
   when both options exist.

Generalisation: applies to operator algebras where the same operator
admits both `.apply` and `.solve` — SN sweep+matvec, CP collision-
probability+inverse, MoC characteristic tracing+source-iteration
fixed-point.  Every operator-algebra carve that ships new sweep paths
must reuse the existing `.apply`-side machinery for whatever's shared
(seed strategies, BC realisers, angular closures).

Cross-reference: `[[lessons-L18]]` (Pattern 7 producer-side — same
spirit: single canonical site for the convention, not duplicated at
consumers); commit `c93355c` evidence (net **−35 LOC** while
strengthening invariants); CLAUDE.md Cardinal Rule 2 (architecture is
critical — single source of truth).

## L22: In a worktree, read the WORKTREE — never the main checkout

When the session runs in a git worktree (`.claude/worktrees/<name>`), the
ORIGINAL repo checkout (`/Users/rodrigo/git/nuclear/ORPHEUS`) is on a
**different branch** and its files **differ**. Reading the main-repo
absolute path — or running `grep` after `cd`-ing to the main root — silently
reads the WRONG branch's source.

2026-06-04 instance (WavefrontFlux carve): the main checkout was on
`refactor/sn-operator-algebra` (2-D matvec = a cell-centred FD stencil with
in-sweep `bc.apply`); the worktree branch `refactor/field-role-typing` had
O.4b fully landed (matvec routes through `graph.residual`, `_compute_gradients`
RETIRED, bare boundary). I `Read`/`grep`'d the **main** path and concluded
"matvec is FD, sweep-only, O.4b pending", recommended a now-pointless "merged
arc", and **falsely accused the explorer of fabricating** `operator.py` — the
explorer had read the worktree (correct) via `spatial/` paths; I had read main.
The user caught it ("I thought we already did that"). Phase 0 de-risk caught
the mechanical symptom (a verbatim sweep copy diverged 0.137 for reflective).

**The discipline:**

1. **Every `Read`/`Grep`/`Bash` targets the worktree.** Use worktree-relative
   paths (the shell cwd resets to the worktree) or the full worktree-prefixed
   absolute path. NEVER `cd /Users/rodrigo/git/nuclear/ORPHEUS` (the env even
   says "Do NOT cd to the original repository root") — a `cd … && grep` chain
   runs the grep in main.
2. **`inspect.getsource(fn)` under `PYTHONPATH=<worktree>` is the source of
   truth** when a sub-agent's file:line conflicts with your own read. The
   conflict is often "we read different trees", NOT a fabrication. Print
   `module.__file__` to confirm which tree is imported.
3. **Before accusing a sub-agent of L12 plausibility-substitution, confirm you
   and it read the same tree.** A real fabrication and a wrong-tree read look
   identical from the main agent's seat; only the tree check distinguishes them.
4. **Nexus queries carry the same hazard.** The session's MCP server is
   launched against the MAIN checkout's graph and is NOT restarted on
   EnterWorktree. After entering a worktree: rebuild Sphinx inside the
   worktree, then call `mcp__nexus__use_workspace(<worktree root>)`
   (nexus ≥ 0.12). `mcp__nexus__workspaces` lists checkouts + graphs;
   the `session_briefing` workspace block warns on graph↔checkout
   branch mismatch.

Cross-reference: `[[lessons-L12]]` (sub-agent fabrication — the failure mode I
wrongly attributed here); the dispatched explorer was correct.

## L23: "Redundant" is a layer claim — check binding-generality before retiring; posing is its own layer

Phase 1 R5 planned to "retire the deprecated `power_iteration`, migrate all 5
solver families onto `KEigenvalue`." Tracing it before coding refuted the
premise on two counts, both worth generalising:

**(a) Two implementations of the same algorithm are not interchangeable unless
they bind at the same layer — and the MORE GENERAL layer is the one to keep.**
`power_iteration(solver: Protocol)` and `KEigenvalue(L,S,F).solve` are the SAME
power-method fixed-point combinator (byte-identical loop bodies). But
`power_iteration` binds the inner resolvent LATE (the opaque
`solve_fixed_source` Protocol method) while `KEigenvalue` binds it EARLY (from an
`(L,S,F)` triple). Late binding is **strictly more general**: it admits the
monolithic-matrix resolvents (CP BiCGSTAB, diffusion FD, homogeneous direct)
that have NO `(L,S,F)` factorisation. The symbol marked "deprecated" was
actually the canonical engine; retiring it would have *lost* generality and
forced fictitious operators onto the matrix methods. The genuine redundancy was
the **duplicated loop**, not the symbol — so the narrow one delegates to the
general one (one loop), and the deprecation arrow was reversed. **Before
retiring "the redundant one": name the morphism both share, then ask which
binds more generally. Keep the general; make the specific delegate.**

**(b) Problem POSING is a distinct architectural layer — do not fold it into the
solver/adapter.** The user caught a layering error: K-eigenvalue, α-eigenvalue,
and transient differ in *how the operators are arranged*, so a generic "adapter"
layer misses the arrangement step. The correct decomposition (cross-domain-
attacker `eigenvalue_posing_layering_frames.md`): **leaves** (method-specific
operators) → **posing** (bifurcated: 2a method-agnostic role-assignment into the
standard form `A_loss ψ = λ M ψ` + the μ→physical-eigenvalue map; 2b
method-specific `A_loss` realisation) → **resolvent** `A_loss⁻¹` (method-specific)
→ **algorithm** (general over the standard form). K vs α differ ONLY in 2a
(K: `A_loss=L+C−S−B, M=F, k=μ`; α: `A_loss=L+C−S−F−B, M=1/v, α=−1/μ`). The
algorithm never knows K-vs-α. The metric lives at the leaf; the adjoint is a
daggered posing row, not a layer.

**How to apply:** when a plan says "X is redundant, retire it," dispatch the
cross-domain-attacker (or trace it yourself) to (1) confirm X and its
replacement are the same morphism, (2) identify which is the more-general
binding layer, (3) collapse the genuine duplication (often a loop/combinator)
rather than deleting the general symbol. When designing a solve that has
problem-type variants (K/α/transient, forward/adjoint), separate POSING (operator
arrangement → standard form) from ALGORITHM (general over the form) from
RESOLVENT (method-specific inverse). Generalise only after ≥2 instances exist
(the `coding-elegance` skill, Pattern 6 "unify after two instances") — document the
other posing rows as seams, don't build them.

Cross-reference: `[[lessons-L21]]` (sweep/matvec = one operator, two
applications — same "reduce strategies, don't add alternatives" spirit);
Cardinal Rule 2 (single source of truth = collapse the loop, keep the general
layer); commits `650032e`+`7603c8e` (net structural reduction).

## L24: Don't execute a plan literally when the ground has moved — re-characterise sub-steps at pickup

Phase 2 (O.2a, the honest `L+C−S−B` driver) surfaced three related traps, all about
**not running a plan's verb mechanically once the situation has shifted under it**:

**(a) A "forcing-function / tripwire" sub-step can be SUPERSEDED by the design choice
it was meant to force.** Plan sub-step 1 was "give `ScatteringOperator` a domain so the
`S+B` fold THROWS — forcing the honest composition." But the chosen design (variadic
`(L_resolvent, *gains)` drivers, user-approved) retired the fold DIRECTLY — so the
tripwire was moot (B is never summed with S again) AND the domain it wanted didn't
exist (no bulk `FunctionSpace`; bulk operators type via `block_role`, not domain
spaces). Executing it mechanically would have minted a speculative `V_bulk` for no
live benefit. **When a design fork lands, re-evaluate every downstream sub-step against
the NEW state — a sub-step whose only job was to force an outcome already achieved is
moot; demote it to a documented seam.**

**(b) Before deferring a "pre-existing red," RE-RUN and CHARACTERISE it — a memory's
red-note can be stale or mis-attributed, and the fix is often cheap.** The memory listed
"precond_safety k≈1.67 / b1pp ×3 / restart ×3" as orthogonal reds to leave. Running +
reading each showed THREE different classes, two cheaply fixable: (i) `precond_safety`
hand-built its Krylov OMITTING the B operand (`(LC, S, ZeroOperator)`) → wrong eigenmode
— a **test-construction bug**, fixed by adding the now-first-class B gain (→ k_inf=1.875);
(ii) `b1pp_constant_flux` asserted the PRE-O.4b boundary semantics (`out.boundary==0`)
— a **stale-expectation** the architecture changed under (O.4b emits inflow-identity +
outflow-defect), fixed by migrating the assertion; (iii) `restart` was production-correct
to ~1.6e-9 but pinned at a strict 1e-9 — **gate-tightness**, fixed by relaxing to a
tolerance that still catches the real 4.7e-1 defect signature. None were "the solver is
wrong." **A red is not self-classifying: re-run it, read the assertion, decide test-bug
vs stale-expectation vs gate-tightness vs real-defect before deferring.**

**(c) A sub-step framed as a "quick re-enable" can be a substantial FEATURE — read the
guarding tests before sizing it (or deferring it).** Plan sub-step 5 (#200) read "re-enable
the sweep-as-preconditioner (`preconditioner=None`)." But `test_krylov_curvilinear_precond_safety`
is a carefully-designed test that PINS the identity-precond contract and explicitly
documents the naive sweep precond as POOR on curvilinear — switching would lock in the
wrong production state. #200's real fix is the block-inverse FACE preconditioner, a
multi-session feature that doesn't exist. **A test that "documents why NOT to do the
obvious thing" is the cheapest scoping signal there is — grep the guarding tests first.**

**How to apply:** at each sub-step pickup — (1) re-read it against the CURRENT
code/architecture, not the plan's authoring-time assumptions; (2) for any "pre-existing
red" about to be deferred, run + characterise it (the classes above); (3) for any "quick
fix / re-enable," grep the guarding tests. The 15-min audit beats both executing a moot
tripwire and deferring a one-line fix.

Cross-reference: `[[lessons-L23]]` (the R5 "redundant"-is-a-layer-claim — same spirit:
trace the premise before acting on a plan's verb); `.claude/rules/process-discipline.md`
("bias toward completion": "no consumer" / "pre-existing red" often means "not yet characterised — look closer");
commits `deb1ce3` (precond_safety test-bug fix) + `33dd5ff` (b1pp stale-expectation +
restart gate-tightness).

## L25: File-internal `replace_all` is safe ONLY if every occurrence is the target concept

A token rename via `Edit(replace_all=true)` scoped to "one file where every X is the thing
I'm renaming" is NOT safe if the file's DOCSTRINGS cross-reference SIBLING methods that
share the base name. S6.2 (#222) ran `replace_all residual → loss_action` on
`loss_representation.py` on the premise "every `residual` here is the strategy's matvec
method." False: three docstrings referenced the DAG's `graph.residual` and
`SweepDependencyGraph.residual_windowed` (DIFFERENT concepts) → corrupted to
`graph.loss_action` / `loss_action_windowed` (refs to nonexistent methods). Cosmetic
(docstrings only; tests stayed green because the functional `graph.residual_windowed(...)`
CALLS live in a DIFFERENT file that was never replace_all'd) — but wrong, and caught only
at S6.3 when re-reading the moved bodies.

**The discipline:** before a file-internal `replace_all <token>`, grep the file for the
token and SCAN every hit for OTHER meanings — sibling-method cross-refs
(`collaborator.residual`), library calls, substring collisions (`residual_windowed`,
`residual_kernel_batch`). If any hit is NOT the rename target, the replace_all is unsafe:
narrow the `old_string` with disambiguating context, or do targeted edits. The "every
occurrence in THIS file is the same concept" premise must be VERIFIED by reading the grep,
not assumed from the file's primary purpose — the same token names a method on THIS class
AND a method on a collaborator class in the same file's prose.

Cross-reference: `[[lessons-L17]]` (crosswalk before carve — same "verify the convention at
each site" spirit); the S6.2→S6.3 corruption (3 docstring spots in `loss_representation.py`,
fixed in the S6.3 commit).

## L26 — Mode-8 scope: pytest's rewriter keeps test-module bare asserts LIVE under `-O` (2026-06-11)

A strict-xfail cylinder L1 test XPASSed under the canonical `python -O` and the first
suspicion was vv Mode 8 (assert compiled out → vacuous body → fake XPASS). A 10-second
deliberate-failure probe settled it: **bare `assert` statements in COLLECTED TEST MODULES
fire under `-O`** — pytest's assertion rewriter transforms them at import time into
explicit code the `-O` flag cannot strip (the `PytestConfigWarning` seen in every `-O` run
says exactly this: "assertions NOT IN TEST MODULES or plugins will be ignored"). So Mode 8's
inert-assert trap applies to HELPER modules, production code, probe/diagnostic scripts, and
anything pytest does not collect+rewrite — NOT to the test files themselves. Implications:
(1) an XPASS of a strict xfail under `-O` is REAL evidence the body passed — investigate the
heal, don't dismiss it as a vacuous-assert artifact; (2) migrating test-module bare asserts
to `np.testing`/`pytest.fail` is a consistency/clarity choice, not a correctness fix;
(3) when auditing for Mode 8, grep the NON-collected surface (helpers, `derivations/`,
`_test_helpers`, production `assert`s), and when in doubt run the deliberate-failure probe —
it is cheaper than any amount of reasoning about the toolchain.

## L27 — Angular-closure residual audits MUST be per-ordinate; weight-summed checks are blind BY CONSTRUCTION (2026-06-12)

During the #195/ERR-058 investigation, a diagnostic computed the
"operator residual of the manufactured solution" as the WEIGHT-SUMMED
(scalar) angular reduction `Σ_n w_n r_n` — and reported O(h²) while the
per-ordinate residual was O(10). The M-M redistribution's α-dome
telescopes under the angular weight sum REGARDLESS of the half-angle
thread values, so a scalar residual cannot see a wrong angular closure
— this is anti-pattern #8 ("particle balance holds" ≠ per-ordinate
balance) instantiated inside a DIAGNOSTIC, where it mis-supported a
"near-singular operator" hypothesis for a full investigation round
(falsified only by a dense SVD: σ_min ≈ 0.9, plus an explicit
per-ordinate `|apply(ψ_ref) − q|` check that read 13.4 where the scalar
read 2e-4).

Operational rules:
1. Any residual/admission audit of an operator with an ANGULAR
   redistribution term must assert per-ordinate (max or volume-weighted
   per-ordinate norm), never only the weight-summed reduction.
2. A closure SEED is part of the operator: enumerate the field family
   each seed is EXACT on (flat-in-space, flat-in-angle, equilibrium) and
   make sure at least one gate sits OUTSIDE that family. Both ERR-058
   seeds were exact on the flat fields every existing gate used (Mode 7
   at the operator-internals level).
3. When two solution paths agree bit-identically (SI ≡ Krylov), that is
   evidence they share ONE discrete system — not that the system is
   right. Post-unification, agreement gates must be re-classified from
   "independent cross-check" to "twin-consistency check" (L4-class).

Catalog: ERR-058 (mechanism + fixes). Promoted gate:
`tests/sn/verification/mms/test_curvilinear_operator_admits_mms.py`.

## L28 — `git checkout`/`git restore` on a tracked path silently destroys uncommitted working-tree state; forbid it in briefs (2026-06-17)

During #240 D5b-S3, a sub-agent ran `git checkout .claude/skills/vv-principles/error_catalog.md`
to revert a mis-placed edit. That file's session-start state was UNCOMMITTED working-tree work
(the ERR-060 + ERR-061 catalog entries, accumulated across prior sessions, deliberately
uncommitted because `.claude/skills/*` lands via the instruction-architecture flow). `git checkout`
reverted it to HEAD and **destroyed both entries from disk** — and because they were never
committed, git history could NOT recover them. They were restored only because the agent had
pasted the verbatim text into its closeout (`## CATALOG RESTORE PAYLOAD`); had it not, the
content was gone. The `@catches("ERR-060"/"ERR-061")` markers in committed tests would have
dangled.

The hazard is acute for this project precisely because of the forbidden-to-commit accumulating
files (`error_catalog.md`, `vv-principles/SKILL.md`): they carry load-bearing working-tree state
that is, by policy, NOT in git — so the usual "git can undo it" safety net does not exist. A
`git checkout`/`git restore`/`git stash drop`/`git clean` on such a path is irreversible.

**How to apply:**
1. **Sub-agent briefs MUST forbid `git checkout <path>` / `git restore <path>` / `git stash` on
   tracked paths.** To undo an edit, re-edit (the agent knows what it changed) or use the editor's
   own revert — never a git-level discard, which cannot distinguish "my bad edit" from
   "uncommitted work that predates me."
2. **Main-agent recovery protocol when it happens anyway:** the destroyed content is recoverable
   ONLY from (a) the agent's closeout paste-back (hence the L12 verbatim-paste-back discipline
   doubles as a backup), or (b) other agent-memory closeouts that quoted it. Restore from there;
   the main agent CAN edit `.claude/skills/*` (the classifier blocks sub-agents, not main) — the
   sub-agent's "classifier-blocked" is not a dead end, it's a hand-off to main.
3. **Generalises** to any uncommitted-by-policy state: `docs/_build/` (regenerated, low stakes),
   but especially the instruction-architecture-managed `.claude/skills/*` accumulators.

Cross-reference: `[[lessons-L12]]` (paste-back — here it was the ONLY backup); the standing
forbidden-to-commit set (`.claude/skills/*`).

## L29 — Re-anchoring an invariant during a refactor: the new check must be AS STRONG; a downstream guard's WEAKER property does NOT subsume it (2026-06-26)

#261 collapsed `CollisionOperator` into `MultiplicationOperator`, removing the held `sn_mesh` and forcing
a re-expression of the `InvertibleOperator` mesh-identity invariant `streaming.sn_mesh is diagonal.sn_mesh`.
The tempting move — "the W-D `OperatorSum` composition guard already checks `domain`/`codomain` equality, so
the explicit invariant is redundant; downgrade it to that" — was WRONG, and an independent cross-domain-attacker
challenge caught it BEFORE it shipped a latent wrong-physics bug.

The invariant's possible strengths **nest**: `object-identity ⊋ geometric-consistency ⊋ shape-equality`. The
guard checks only the WEAKEST (`FullFieldSpace.__eq__` is `(name, shape)`), so two equal-shape /
different-volume meshes pass the guard yet compute wrong physics (the WDD sweep threads the diagonal's σ
against the STREAMING geometry — they must be the SAME mesh, not merely same-shaped). The essential
requirement was the MIDDLE tier. The correct re-anchor kept that strength through a handle that SURVIVED the
refactor: `streaming.sn_mesh is diagonal.coefficient.mesh` (the `CrossSectionField` still carries `.mesh`).

**How to apply:** when a refactor deletes the field an invariant checked and you must re-express it — do NOT
relax it to whatever a nearby guard happens to verify. (1) Enumerate the property's strength tiers
(identity / value-equality / shape / …). (2) Determine which tier the invariant ACTUALLY needs by tracing
what BREAKS on violation (here: the σ↔geometry pairing in the sweep). (3) Re-anchor AT that tier through a
surviving handle. (4) Mutation-verify the re-anchored check still REDS on a real violation (here
`test_mismatched_mesh_rejected` still fires). A guard checking a weaker property is never a substitute for a
stronger invariant — "redundant with X" is false whenever X is strictly weaker.

Cross-reference: `[[lessons-L24]]` (re-characterise sub-steps when the ground moves — same "don't mechanically
assume redundancy" spirit); `[[lessons-L23]]` ("redundant" is a layer claim — check binding-generality before
retiring); the independent-challenge discipline (dispatch an adversarial reviewer for your own
about-to-relax-an-invariant reasoning — it caught this).

## L30 — Before folding N call-sites onto ONE typed abstraction, verify they share the same OPERATION, not just the same data (2026-06-26)

A "unify the seven untyped `compute_keff` functionals onto one `ReactionRateFunctional`" plan (a
cross-domain-attacker re-seat recommendation) survived design review but COLLAPSED on contact with the actual
code: the fission **row-factor** `⟨νΣf,φ⟩` contracts the **group** axis (per-cell density, keeps space), while
the keff estimators `∫νΣf·φ dV` contract the **space** axis (per-group rate, keeps groups). They are DIFFERENT
functionals — reconciling only at the fully-contracted scalar total (Fubini). "Seven copies of the same
co-vector" was an illusion: same DATA (νΣf), different OPERATION (which axis). Forcing the fold would have
either duplicated the production computation (the Cardinal-Rule-2 smell the whole carve REMOVED) or retired
meaningful per-group diagnostics. Compounding it: each solver's keff is BESPOKE (SN (n,2n), CP net-removal,
diffusion leakage, MoC per-region, homogeneous 0-D) — not a uniform fold at all.

**How to apply:** when a plan says "type/unify these N call-sites with one abstraction," BEFORE building it
(1) READ each call-site's actual reduction — which axis does it contract, what is the output shape/units? Two
sites sharing a coefficient field are NOT the same functional if one contracts groups and the other contracts
space. (2) The honest unit is the SHARED operation, not the shared data; if the sites differ, the abstraction
is `⟨w,·⟩`-the-row-factor for some and `∫⟨w,·⟩dV`-the-scalar for others — possibly RELATED (one the degenerate
or volume-integral of the other) but DISTINCT types. (3) Mid-execution, when reading the real bodies
contradicts the plan's "uniform fold," STOP and reshape — don't force it (here: minted `IntegratedReactionRate
= ∫⟨Σx,φ⟩dV` as the volume-integral of the density — a distinct type composing the row-factor + the measure —
and routed only SN; the other four deferred as bespoke, #270). A structural reviewer judges the plan; only
reading the code judges the fit.

Cross-reference: `[[lessons-L24]]` (re-characterise sub-steps when the ground moves — same spirit, here
triggered by reading the bodies); the test-architect "read all five `compute_keff` bodies first" reshaping is
the concrete instance.

## L31 — A sub-agent's VERDICT and its CAUSAL ATTRIBUTION are separable; verify the attribution against git/issues, not just the verdict (2026-06-26)

The numerics-investigator correctly verified the current curvilinear sphere matvec is CORRECT and the red
snapshot tests merely STALE (the VERDICT — backed by L0/L1 references + a mutation-proven probe), but
MIS-ATTRIBUTED the cause to ERR-058 `3b088ee`. A `git log -- <fixture>` plus issue #250 showed the two
snapshot stores were refreshed AT the ERR-058 fixed point (`798372f`) and again at #240; the true
uncaptured diverging commit was the LATER `b2d8a6d` (spherical Morel–Montry τ-unclamp, Bailey Eq. 43,
Refs #229). The investigator read the ORIGINAL capture commit and missed the refresh history — its verdict
was unaffected, only the "why" was wrong.

The discipline: "current is correct / the snapshot is stale" (verdict) and "commit X caused the drift"
(attribution) are DIFFERENT claims with DIFFERENT evidence; a rock-solid verdict does NOT co-validate the
attribution. Before writing a root-cause into a docstring/commit (Cardinal Rule 1), verify it with
`git log -- <file>` (the fixture's refresh history settled it in one command), `git show --stat <commit>`,
and the tracking issue — never transcribe a sub-agent's "because commit X" verbatim. Same session, the same
`trust git` reflex caught origin/main sitting 64 commits behind a "ff-merged to main" memory claim.

Cross-reference: `[[lessons-L12]]` (sub-agent closeouts — here honest-but-incomplete forensics, NOT
fabrication); `[[lessons-L11]]` (verify, don't trust); `.claude/rules/process-discipline.md` (trust git,
not frozen memory).

## L32 — Dedup a shared concern at the container's OWN lowest layer; never route the lower-layer verb through the higher one (2026-06-27)

P5.5 (energy-condensation reshape): two verbs assemble the same `Mixture` container from collapsed
cross-section channels — `Mixture.condense` (the ENERGY collapse, a `data`-layer method) and
`MaterialXSField.project_through` (the SPATIAL collapse, a `transport`-layer method). They DUPLICATED the
channel-assembly taxonomy (the `csr`-wrap + `eg`-thread of SigT/SigC/.../SigS/Sig2/χ → a `Mixture`). An
elegance review's shape ruling said "dissolve `Mixture.condense`; the single home is `project_through`." That
was **layer-illegal**: `Mixture` is `data`, `project_through` is `transport`, and `data ↛ transport`, so a
`data`-layer verb CANNOT route through a `transport`-layer verb. The correct refinement (which a fresh
elegance pass then ruled SOUND): extract the shared assembler at the container's OWN lowest layer —
`Mixture.from_dense_channels(...)` in `data` — and have BOTH verbs call DOWN to it (`condense` in data;
`project_through` transport→data, allowed). The lower verb stays where it lives.

The generalizable rule: when N verbs at DIFFERENT layers share a concern, the dedup target is a primitive at
the **deepest** layer any of them sits in (or below), called by all — NOT the verb at the highest layer
absorbing the others. "Single home = the existing verb that happens to do it" is wrong whenever that verb
sits ABOVE a consumer; the home is the shared primitive at/below the lowest consumer. Corollary on what to
extract: factor out only the part that is genuinely layer-and-axis-INVARIANT (here the assembly mechanics —
`csr`-wrap + `eg`-thread), and leave the part that legitimately differs per verb (here the weighting
taxonomy: energy marginalize-vs-average vs spatial flux-vs-production frames) in each verb — that asymmetry
is real domain content, not a twin.

**How to apply:** when a plan/review says "dissolve X into Y" or "Y is the single home," (1) check the LAYER
of X, Y, and every consumer; (2) if Y is above any consumer of the shared concern, reject "route through Y"
— site the shared primitive at the deepest consumer's layer instead; (3) verify the extracted unit is the
invariant slice (mechanics), not the per-verb-varying slice (the weighting/axis). Generalises to every future
cross-method coarsening verb (CP region-averaging, MoC FSR-merging, GEC rank>0) that assembles a shared
container from a higher and a lower layer.

Cross-reference: `[[lessons-L29]]` (don't relax an invariant to a weaker downstream check — same "scrutinize
the tempting simplification" spirit); `[[lessons-L24]]` (re-characterise a plan's verb when the ground —
here the layering — contradicts it); `.claude/rules/coding-standards.md` ("clean before extending" — the
assembler is the single generic body both verbs extend through).

## L33 — In a doc-repair campaign the corpus is an ADVERSARY: never cite a doc/docstring as evidence for a structural claim, and never chain a commit after an unread verification (2026-07-14)

Two self-inflicted failures in one session, both while *fixing* documentation:

**(a) A stale docstring manufactured a wrong architecture.** Designing a corpus-wide
transport-method taxonomy, I asserted "CP has no `(L,S,F)` factorization — CP *is* the
propagator," citing `numerics/eigenvalue.py:37-43`. The docstring was stale and **the code it
describes refutes it**: `diffusion/solver.py:240` = `leakage + collision - scattering - boundary`;
`homogeneous/solver.py:143-146` = `collision - k_iso`; `cp/solver.py:519-520,526` applies S/(n,2n)/F
explicitly *outside* `P_inf` — and `P_inf` **is** CP's `(L+C)⁻¹`. An adversarial reviewer
(cross-domain-attacker) killed the taxonomy and traced it to the docstring. Same session, the same
disease bit twice more: the SN machine header I had shipped that morning encoded the retired
pre-B-extraction fold (`L: streaming + boundary`), contradicted by its own page 14k lines later,
and six more sites across three pages asserted the retired "four-operator algebra" as current fact.

> **The rule: a doc or docstring is a CLAIM, not evidence.** Before building any structural
> argument on a documented assertion, verify it against the code (`grep`/Read the implementation).
> This is sharpest *during a documentation-repair campaign*, where the corpus you are repairing is
> exactly the corpus that will mislead you while you repair it — and where a fresh agent's first 50
> lines (a machine header, a module docstring) are the highest-leverage place for a lie to sit.
> Corollary: when a doc and the code disagree, that is not a typo — it is an architecture symptom
> (one fact with no single source of truth). Fix the source, then grep the corpus for its siblings.

**(b) I claimed a verification I never read.** I ran `sphinx -W` with `run_in_background`, then in
the *next* command chained `cat <log>` and `git commit` with a newline (not `&&`). The build had
FAILED (on a bug I introduced: `\rm` inside a non-raw docstring → `\r` → broken RST), the log
printed "BUILD FAILED", and the commit landed anyway carrying the message "Clean -W build." That is
vv-principles **L12 committed against myself** — the exact fabrication pattern I demand sub-agents
never do, in a session where `vv-principles` was preloaded.

> **The rule: never sequence a commit after a verification whose result you have not READ.** Read
> the exit code, then decide. Do not put a verification and a commit in the same Bash call unless
> they are `&&`-chained. If a commit message asserts a gate ("clean build", "tests pass"), that
> string is a claim under L12 — it needs the same paste-back discipline as a sub-agent closeout.
> (Recovery: amend while unpushed so the false claim never becomes history.)

Cross-reference: `[[lessons-L12]]` (sub-agent closeouts fabricate at the summarization boundary —
here the fabricator was the main agent); `[[lessons-L31]]` (verdict vs attribution — same "verify
before transcribing" reflex); `.claude/rules/process-discipline.md` (trust git, not a frozen claim
— generalizes: trust CODE, not a frozen doc). Evidence: `275a753a`, `018ecb7b`, `0ca0d378`;
`.claude/plans/documentation_corpus_architecture.md` §3.6 (the 12 MUST-NOT claims exist so the next
author cannot repeat this).
