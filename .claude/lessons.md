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

- `StreamingCoefficientCache (N, nx)` — built once at SNMesh
  construction; NEVER mutates across the run. *(Shipped 2026-05 as
  `GeometryCoefficients`; renamed 2026-08-26 — see the two banners
  below.)*
- `CollisionCache (N, nx, ng)` — built when σ_t binds; rebuilt only
  on depletion / thermal-feedback.

This is the operator-algebra `(L + C)` rehearsed on the cache layer:
L (streaming + curvature) lives in `StreamingCoefficientCache`; C joins
via `1/(g_streaming + Σ_t·g_volume)` to form `CollisionCache`.

⛔ **`SweepCoefficientCache` is the REJECTED name above and must never
be reused.** It is recorded here as the wrong shape, so a class wearing
it would make this lesson read as condemning the shipped design — and
this file is loaded at every session start. `[M]` 2026-08-26: no class
by that name ever existed (`git log -S "class SweepCoefficientCache"`
is empty); it was a proposal in the cross-domain-attacker memo. It was
proposed *again* as the rename target for `GeometryCoefficients` and
caught only because the sweep's residual grep surfaced this entry.
⟹ this is `plan-authoring` §3's ambiguous-name hazard: a refutation
attached to a NAME condemns every later design that shares it, so the
disambiguation belongs here, at the banner, where a summariser reads.

⚠⚠ **And L15 applies one level down, to the class it holds up as
correct.** `[M]` 2026-08-26, three meshes against one quadrature
(uniform nx=6, uniform nx=20, GRADED nx=6): the "right shape" class
mixes **three** invalidation strata, not one — 7 fields bit-identical
across all three meshes (`abs_mu`, `c_in`, `c_out`, `tau_inv`,
`mm_a_in_coeff`, `is_degenerate`, `level_ordinates`), 4 that differ on
every re-mesh (`A_down`, `A_total`, `dA_w`, `V`), and 2 that turn on
ordinate sign and cell COUNT only (`chain_idx`, `chain_idx_inv`). So
the whole object rebuilds on any re-mesh including the 7 that provably
cannot change. **The worked example of a lesson is not exempt from
it** — and nothing prompts you to re-apply a rule to its own exhibit,
because the exhibit is what taught you the rule. The split is
scheduled at the un-weld plan §4; the name was made stratum-agnostic
in the meantime so it does not bless the weld.

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
three consumers — `StreamingCollisionOperator.solve`, `_solve_krylov`,
`_solve_source_iteration` — each had to apply `* sum_w` or `/
sum_w` to bridge to their expected convention. Three habitats for
the bug to drift; debugged twice (Step D → fix in Krylov; Step E →
fix in SI; the third bridge survived in `StreamingCollisionOperator.solve`
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
`StreamingCollisionOperator` does NOT advertise it; default-no-precond
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
`StreamingCollisionOperator.solve` (which read `rhs(1)` history for the Carlson
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
`StreamingCollisionOperator.solve(rhs, *, initial_guess=None)` — Cardinal
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

> **⚠ MECHANISM SUPERSEDED (2026-07-28) — the PRINCIPLE stands, the plumbing
> above is archaeology.** The `initial_guess` kwarg described here no longer
> carries the sweep's Carlson seed. **#282 route (a)** retired the
> extrapolated-iterate channel entirely: the curvilinear seed is now the
> composite's **first-class ψ½ STATE** (`radial_characteristic.cells(p, -1)`,
> read in `MorelMontryAngularSweep.precompute_psi_state`), i.e. upstream state
> in the augmented walk order rather than a previous iterate threaded inward.
> `StreamingCollisionOperator.solve` (renamed from `InvertibleOperator`,
> 2026-07-28) now **accepts and DROPS** `initial_guess` — it is an exact direct
> inverse with nothing to seed. The kwarg survives ONLY as
> `SupportsSeededApply` protocol conformance (`numerics/iteration.py`), which
> `SourceIteration` type-checks its `A_inv` against and threads unconditionally;
> six other exact inverses accept-and-drop identically. A warm start today is
> purely an **iteration-layer** concept — it enters through the RHS
> (`rhs = q + Σ gains·ψ_prev`), never through the sweep's interior recurrence.
> The one operator that genuinely consumes a seed is the *iterative*
> `GreenOperator`.
>
> **What still holds** (and is why this lesson stays): sweep and matvec are two
> applications of ONE operator, and a seed strategy is a property of the
> operator, not of the path. #282 route (a) is that principle applied *harder* —
> it removed the strategy parameter instead of duplicating it.

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
a re-expression of the `StreamingCollisionOperator` mesh-identity invariant `streaming.sn_mesh is diagonal.sn_mesh`.
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
`.claude/plans/archive/documentation_corpus_architecture.md` §3.6 (the 12 MUST-NOT claims exist so the next
author cannot repeat this).

## L34 — A path built from SEGMENTS is invisible to a path-grep; and know which reference class your build actually gates (2026-07-15)

Corpus Phase A moved 33 documentation pages. The move's blast-radius audit was three greps
(`:doc:` roles, toctree entries, raw `docs/theory/...` strings) and it found 258 real sites. It
**could not see the single biggest hazard**, because that hazard contains no greppable path:

```python
DOCS_THEORY_DIR = REPO_ROOT / "docs" / "theory"     # no "docs/theory/" substring exists
```

`tools/verification/generate_capability_matrices.py` assembles its output directory from
*segments*, fires on Sphinx's `builder-inited`, wipes with a **non-recursive** glob, and rewrites
its `.inc.rst` outputs. Post-move it would have wiped nothing (the real files had moved out of its
glob), written **duplicates** at the dead flat path, and orphan-errored — while every grep in the
audit reported clean. A dispatched **explorer** found it by *reading the tool* rather than
searching for the string. Compounding it, the pinning test carried a **twin copy** of the same
constant (`tests/.../test_capability_matrices.py:29`) — two sources of truth for one path, which
is what made the drift possible at all (Pattern 2).

> **Rule 1: when auditing a move/rename, grep the LAST SEGMENT (`"theory"`, `"peierls_nystrom"`),
> not the joined path** — and read any tool that *writes into* the moved tree. `Path` /
> `os.path.join` / f-string assembly, and any config value built the same way, are structurally
> invisible to the obvious search. Corollary: a path constant gets **one** home; consumers
> (including tests) **import** it.

The same phase corrected a load-bearing false belief about the build gate. `.claude/rules/
coding-standards.md` claimed an unresolved `:func:`/`:class:`/**`:ref:`** "renders as plain text
with no `-W` warning." A 60-second probe (a scratch Sphinx project with one deliberately-broken
ref of each species) measured otherwise on Sphinx 9.1.0:

| broken ref | warns? | gated by `-W`? |
|---|---|---|
| `:doc:` → missing page | **yes** (`ref.doc`) | ✅ |
| `:ref:` → missing label | **yes** (`ref.ref`) | ✅ |
| `:func:` / `:class:` / `:mod:` | **no** (needs `-n`) | ❌ |
| raw path string in prose/docstring | **no** — it is just text | ❌ |

> **Rule 2: know which reference class your gate covers, by measurement, before you rely on it.**
> Here it *inverted* the plan: page moves and label retirements ARE build-gated (so `:doc:`/toctree
> churn needs no hand-audit — let `-W` find it, and it found all 7), while the classes that need a
> hand-built gate are the Python-domain roles and **raw text**. Phase A's raw-path gate was a
> filesystem check (57 live / 0 dead) because no build can do it. It found 3 dead pointers, all
> **pre-existing**: `axis.py` promised `docs/theory/sn_dim_agnostic.rst` — a page git shows was
> **never written** — citing "D9 of `<a plan that no longer exists>`", which is *precisely* the
> untraceable reference #231 exists to eliminate.

**The hazard is not hypothetical: it had already fired 3× in this repo**, including on this very
branch — the `galerkin_projection.rst`→`frame.rst` rename two days earlier orphaned 2 raw paths;
`peierls_greens.rst`→`trajectory_resolvent.rst` orphaned 16. **Every rename owes a raw-path
sweep**; the build will not do it for you.

> **Rule 3 (learned the same day, by re-running the gate at compaction): "archaeology" is a
> per-FILE judgment, not a per-DIRECTORY one.** Phase A excluded `.claude/plans/` from the sweep
> wholesale, reasoning that a plan's `← discrete_ordinates.rst` names a move's *source* and
> rewriting it corrupts the record. True — but true of **exactly one file**
> (`documentation_corpus_architecture.md`, which describes the moves). Every *other* plan carried
> **live pointers meant to be followed** ("the ruling is archived in X", "Doc drift for close-out:
> `X:281`" — the latter in an *active* campaign roadmap). The over-broad exclusion left **9 dead
> pointers across 7 plans**, including the 2 the explorer had explicitly warned me about. Before
> excluding a directory, ask *what the text DOES* in each file: a pointer to be followed gets
> swept; a record OF the move does not. Also: a two-hop rename (`galerkin_projection` → `frame` →
> `foundations/frame`) needs the **historical alias** in the map — the current-name map alone
> silently misses every pointer that predates the first hop.

Cross-reference: `[[lessons-L20]]` (retirement requires a dependency audit — this is the audit's
blind spot); `[[lessons-L33]]` (a doc is a claim, not evidence — here the *rule* was the wrong
claim, and measurement settled it); `.claude/rules/coding-standards.md` (corrected in `08e58ee6`);
`.claude/plans/archive/documentation_corpus_architecture.md` §7.1.

## L35 — Moving a labelled block: `-W` gates the LINK, never the PROSE around it — two silent-falsehood classes need a bidirectional grep (2026-07-16)

L34 established that `-W` gates broken `:ref:`/`:doc:` (so Phase A trusted the build to find
page-move/label churn — and it did). Phase B moved four labelled blocks *between* pages and
found the limit of that trust: because `:ref:`/`:eq:` are **path-immune** (a global label
resolves from anywhere), a move breaks **no link** — the build stays green — while the PROSE
the author wrote *around* those links silently goes false. Two classes, neither `-W`-visible:

1. **Directional words.** "`:ref:`windowing-retyped` **above**`" was true when the anchor sat
   above on the same page; after it moves to another page it is nowhere near "above". The
   `:ref:` still resolves (green build); the word "above" is now a lie.
2. **Page-qualifiers.** "`:ref:`bc-extraction` **in :doc:`…/operator_algebra`**`" named the old
   home; after the block moves to `boundary_conditions` the qualifier sends readers to the
   wrong page. The `:doc:` resolves (the old page still exists — green build); the claim is
   false.

Both are the D5 disease the campaign exists to kill (a page that MISINFORMS), and both are
invisible to every build. The catch is a **grep run THREE ways per move**:

- **staying-content → moved-label** (on the SOURCE page): refs to the departed labels carrying
  "above/below/here/this section" (now cross-page) OR a `in :doc:`<source>`` self-qualifier.
- **moved-content → staying-label** (on the moved block): refs BACK into content that stays,
  same two classes.
- **bystander page-qualifiers** (any THIRD page): a tree-wide grep for `<moved-label>` co-located
  with `:doc:`<old-page>``. Block 3's two move-files were clean; **3 stale qualifiers hid in
  `sn/index.rst`, a bystander that was neither source nor dest.**

> **Rule: after moving a labelled block, grep the moved-label set three ways (source-page refs,
> moved-block refs, bystander qualifiers), filtering for directional words and `:doc:`
> qualifiers. The green `-W` build proves the LINKS resolve; it says NOTHING about whether the
> surrounding prose still tells the truth.** Separate a real falsehood from a false positive: an
> *intra-block* "below" survives (the whole block moves together, internal order preserved); a
> demonstrative "that section" bound to the immediately-preceding ref survives cross-page; a
> *categorical* "one locus down" (cell↔face hierarchy) is not document-position. Fix only what
> the move actually falsified.

Operational note for **Phase C** (which moves ~63k lines): the sub-agent's default scan is
`:doc:`-forward and dest-file-focused; it MISSES the reverse direction and the bystanders. I
found the Block-2 falsehoods post-hoc, then folded the bidirectional+bystander grep into the
Block-3 brief and the archivist executed it cleanly — **put the three-way grep in the brief, or
run it yourself before committing.**

Cross-reference: `[[lessons-L34]]` (which reference class the build gates — this is its prose
corollary: even a GATED role leaves ungated prose); `[[lessons-L33]]` (the corpus is an
adversary that MISINFORMS — a move is a mechanical way to mint fresh misinformation);
`.claude/plans/archive/documentation_corpus_architecture.md` §7.2.

## L36 — Moving a label between files: the incremental Sphinx build raises a PHANTOM `duplicate label` — force `-E` for the gate (2026-07-16)

[[lessons-L35]] is about the PROSE a label move falsifies; this is about the BUILD TOOL itself
lying *during* the move. When a `.. _label:` moves from page A to page B, an INCREMENTAL
`sphinx-build -W` (reusing the saved `docs/_build/doctrees` environment) can raise
`WARNING: duplicate label <x>, other instance in <A>` — even though A no longer contains it.
Sphinx re-reads the changed/new files, but the label registry from the cached environment still
holds A's stale registration when B's fresh copy is registered; the two collide transiently. A
plain `-W` *rerun* then flips to exit 0 once the environment settles — so on a label move
NEITHER incremental result is trustworthy (one false-RED, the next false-GREEN).

> **Rule: the Phase-C / any-label-move build gate is `sphinx-build -E -W` (`-E` = discard the
> saved environment, full re-read), NOT plain `-W`.** `-E` makes the label registry consistent
> in one pass. Confirm independently that each moved label is single-homed:
> `grep -c '<label>' <source>` (must be 0) + tree-wide `grep -rn '^\.\. _<label>:'` (exactly 1).
> Do NOT `rm -rf docs/_build` to force it — that trips the host permission gate; `-E` achieves
> the same without deletion.

Worked (ch3 `angular_quadrature` extraction): the `--keep-going` incremental build showed the
phantom `duplicate label quadrature-types`; `grep -c` proved it 0-in-source / 1-tree-wide;
`-E -W` came back "build succeeded", 0 warnings.

Cross-reference: `[[lessons-L35]]` (the prose corollary), `[[lessons-L34]]` (which reference
class the build gates); `.claude/plans/archive/sn_split_catalog.md` (the Phase-C build-gate section).

## L37 — Editing Python source under a RUNNING full-tree gate false-reds every `inspect.getsource`/AST self-inspection test (2026-07-24)

The pre-merge full-tree gate runs ~51 min serial. Editing a Python source file on disk
mid-run — even a pure-DOCSTRING edit — false-reds any test that source-inspects the
edited module: `inspect.getsource(obj)` resolves the object's line RANGE from the code
object loaded at collection time but reads the BYTES from the file *now on disk*, so an
insertion above the object shifts the window and `ast.parse` dies with the tell-tale
`File "<unknown>", line 1` SyntaxError. (#310 C3: the C3-b docstring fixes landed at
~minute 40 of the C3-a full tree → `test_one_dim_walk_is_orientation_object_not_boolean`
+ `test_frame_is_kernel_parameterized_not_boolean` red with exactly that signature;
both green on the frozen tree, 33/33.)

> **Rule: while a full-tree (or any long) pytest gate is running, Python sources are
> FROZEN — stage doc/plan/markdown work instead, and land code edits either before
> launch or after the verdict.** If a long run reds an AST/source-inspection test with
> `File "<unknown>", line 1`, FIRST check whether the file was edited during the run
> (`git log --since` vs the run's start time) and re-run the failing files on the
> quiescent tree before treating it as a regression. A gate run concurrent with edits
> gates NEITHER the old nor the new tree.

Cross-reference: the nohup/E-core pacing note (`reference_test_execution_env`) — long
runs invite parallel work; this lesson bounds WHICH work is safe (non-imported files
only).

## L38 — A deliberate mutation in the SHARED working tree is indistinguishable from a production bug to a parallel agent (2026-08-10)

[[lessons-L37]] says Python sources are FROZEN while a long gate runs, and bounds
which work is safe. This is its multi-agent twin, and it bites in the opposite
direction: not "my edit false-reds a running test" but **"my edit makes a parallel
agent report a defect that does not exist."**

#340 N6. I flipped one production line to measure a mutation's blast radius —
`if history.converged:` → `if history.fully_converged:  # M8 PROBE — REVERT ME` —
and ran a 4-minute slice. A `test-architect` was concurrently writing gates against
the same tree. It read the guard, then watched a solve with `converged=True` emit a
warning that the guard it had just read says is impossible, and was one step from
reporting **"the landed guard is broken."** It was saved only by running
`inspect.getsource` on the anomaly, which showed it the probe comment and an mtime
of 05:59:29.

The failure is *mine*, and the near-miss is instructive about which half I got
wrong. I DID warn it — "I am running your M8 pre-measurement myself right now …
Do not run M8 yourself — you would collide with my run." That warns about
**duplicated work**. It says nothing about **a mutated file**, and the agent had
no reason to read "collide" as "the source you are about to read is lying."

> **Rule: a mutation to a tracked file is a change to every concurrent reader's
> ground truth. Before mutating, either (a) serialise — no agent dispatches live;
> or (b) state the mutation explicitly in every live agent's brief: the file, the
> line, the exact replacement text, and when it reverts.** "I am running an
> experiment" is not that statement. A comment marker in the mutated line
> (`# M8 PROBE — REVERT ME`) is necessary and was what actually rescued this, but
> it only helps an agent that already suspects the source and goes looking.

Corollary for the reader, which generalises past mutations: **when an observation
contradicts source you have read, re-read the source with `inspect.getsource`
before concluding the code is wrong.** Your file read and the running process can
disagree — because another agent edited between them, because an untracked shadow
copy serves the import (the `coding-standards` mass-delete hazard), or because the
module was imported before the edit. All three look identical from the seat of the
agent forming a verdict, and all three are cheap to discriminate and expensive to
get wrong.

⭐ The collision also handed the campaign M8 for free — the agent's bracketed
re-measurement of the flipped guard (`130 passed`) independently reproduced mine.
That is luck, not a mitigation.

Cross-reference: `[[lessons-L37]]` (the single-agent form — edits under a running
gate); `[[lessons-L22]]` (you and a sub-agent read different trees — same
"confirm you share ground truth before accusing" reflex); `[[lessons-L12]]`
(sub-agent fabrication — here the agent was RIGHT to doubt itself, and the
discipline that saved it was verification, not caution);
`.claude/rules/process-discipline.md` (mutation-testing an uncommitted file:
never `git checkout` to revert — this adds *never mutate silently beside a live
agent*).

---

## L39 — The COMPLEMENT of a guard reaches the states its partner never visits (2026-08-10)

`[M]` #340 N6b step 2. `_certify_within_group_exit` had been calling
`evaluate_residual` on every mesh in the suite for months, guarded on
`record.converged`. The new exit diagnostic was written as its deliberate
complement — same operator, same iterate, same rhs shape, guarded on
`not record.fully_converged` — on the reasoning that a call already proven safe
under one guard is safe under the other.

It is not. `evaluate_residual` REFUSES a bare System-A residual on a mesh that
carries starting-direction levels (it would silently omit `r_B` — the Mode-12
blindness the split-residual mint exists to prevent). The certificate had never
met that refusal, because on curvilinear meshes the truncated solves are exactly
the ones it returns early from. **The complement inherits the partner's code
path and NONE of its field coverage.** The slice went **9 → 25** reds, all 16
new ones curvilinear.

> **Rule: when you add a consumer as the complement of an existing guard, its
> preconditions are UNVERIFIED for your call, however long the partner has run.
> Enumerate the callee's own refusals — read its guard clauses, not just its
> signature — and ask which of them the partner's guard has been hiding.**

Three things make this worth its own entry rather than a footnote:

1. **It is invisible to review by construction.** "This exact call already
   exists twenty lines up and has run for months" is a *true* statement that
   points the wrong way. The diff shows a duplicated call; what changed is the
   set of states it runs in, which the diff cannot show.
2. **The tell is available for free, before any run** — the callee raises, and
   its `raise` sites enumerate its preconditions. One read of
   `evaluate_residual`'s guards would have surfaced both the LD refusal (which
   an agent DID find, because it is documented in the partner's docstring) and
   the carrying-mesh refusal (which nobody found, because it is not).
   ⭐ A precondition documented only in the *callee's* raise is exactly the one a
   partner's docstring will not warn you about.
3. **Guard-complementarity is otherwise a GOOD pattern** and this is not an
   argument against it. The certificate asserts when convergence was claimed,
   the diagnostic reports when it was not; no solve pays for both applies. Keep
   the shape — just do not let "complementary guards" read as "shared coverage".

Cross-reference: `vv-principles` #17 (verify the instrument on a known positive
before trusting a negative — here the instrument was the *reasoning*, and the
known positive was a mesh class the partner never reached); `[[lessons-L30]]`
(same data, different operation — its sibling: same call, different states).

## L40 — A STOPPING CRITERION can be Mode-12 blind to the very mode its loop exists to converge (2026-08-11)

**Context.** #340 N4 gave MoC's inner sweep loop a real stopping rule. The
obvious criterion was `‖Δφ‖_F/‖φ‖_F` — the quantity MoC's OUTER already
measures, already implemented, free to reuse.

**`[M]` It is nearly blind, and on the SIMPLEST fixture it fails hardest.**
`moc_cyl1D_1eg_1rg`, cold boundary flux:

| sweep | `‖Δφ‖/‖φ‖` | `‖Δψ_b‖/‖ψ_b‖` |
|---|---|---|
| 4 | **`0.000000e+00`** | `3.488e-02` |
| 12 | `0.000000e+00` | `3.104e-04` |
| first `< 1e-5` | **2** | **18** |

A break on `‖Δφ‖` alone stops at sweep 2 with the boundary angular flux four
orders from converged — and the loop's own docstring says it exists "to
converge the boundary angular fluxes".

**The mechanism, and why it generalises.** φ is a **volume moment**; the
cyclic-track closure's slow mode lies in that functional's stabiliser. A
1-group problem has no group coupling to iterate, so Δφ collapses in one pass
while the geometric feedback around reflective tracks is untouched. This is
`vv-principles` **Mode 12** applied to a *stopping criterion* instead of a
*test*: the standard question is "can this gate SEE the error?"; the same
question belongs to every convergence test — **"can this criterion see what
the loop is FOR?"**

**The discriminator, and it is cheap.** Compare the loop's measured quantity
against the loop's **carried STATE**. If the state is bigger than the
measurement — MoC carries `_fwd_bflux`/`_bwd_bflux` across outers while
measuring only a volume moment of φ — then some component of the state is
unmeasured, and it is exactly where a slow mode hides. `IterationRecord.criteria`
is a TUPLE precisely so a level can declare more than one; declaring both makes
*which mode binds* readable instead of assumed.

⚠ The counter-intuitive part worth remembering: the fixture where the blindness
was total was the **1-group** one, i.e. the case a reviewer would reach for as
the easy sanity check. Simplicity removed the coupling that would have made the
proxy work.

## L41 — A probe that FREEZES one half of a coupled pair can be measuring its own freeze (2026-08-11)

**Context.** After N4, MoC's outer 2 cost as much as outer 1 (`[172, 176, 3]`)
even though the boundary angular flux is carried across outers. To isolate the
cause from the outer's physics, I probed with a **fixed** fission source:
repeated `solve_fixed_source` on one instance.

**`[M]` The probe read 176 sweeps per call — identically BEFORE and AFTER the
real fix.** First reading: "hypothesis refuted." Correct reading: *this probe
cannot see it.* `solve_fixed_source` renormalises the flux it returns, so
holding `q` fixed while φ is rescaled is itself an inconsistency the probe
introduced — and it dominated the one under study. The genuine defect (the
normalisation rescaled φ but not the carried ψ_b, splitting one solution
vector's scale in half) shows up only end-to-end, where the outer recomputes
`q` from the normalised flux: **64 → 32, 351 → 177, 259 → 128** sweeps.

**The rule.** Freezing a variable is the standard way to isolate a term, and it
is safe only when the frozen variable is genuinely INDEPENDENT of the one under
study. When the two are coupled by a shared normalisation, an update rule, or a
conservation constraint, the freeze introduces a NEW inconsistency whose size is
unknown — and a null result then says nothing about the hypothesis. Before
believing a refutation from a frozen-variable probe, ask: **what does the freeze
itself break, and is that bigger than the effect I am looking for?**

Sibling of `vv-principles` #17 (verify the instrument on a known-positive before
trusting a negative): here the instrument was not broken, it was
*mis-configured* — and the failure direction was the flattering one, "your
hypothesis is wrong" rather than "your probe is blind".

## L42 — `python -O script.py` strips the probe's OWN assertions; the canonical `-O` is safe only inside pytest (2026-08-11)

**Context.** #340 N4.7 moved the ConvergenceWarning emitter out of `sn/solver.py`.
The bit-identity gate was a standalone probe: lift the OLD function's source out
of `git show HEAD:` with `ast`, exec it, and compare its message against the new
one across every advice arm. I ran it the way this project runs everything —
`.venv/bin/python -O probe.py`.

**`[M]` It printed `7/7 character-identical` and `both versions SILENT on a
fully-converged tree`. Three of its own checks had not executed.** `-O` sets
`__debug__ = False` and the compiler **omits every `assert` statement**, so
`assert len(caught) == 1`, `assert old_msg != new_msg` (the check that the
intended clause changed AT ALL) and `assert not caught` (the entire
converged-silence claim) were absent from the bytecode. The output was
byte-identical to a genuine pass. Re-run without `-O`: all three fire, all three
green — but that was luck, not evidence.

**Why this is a trap specific to THIS project.** `.claude/rules/vv-testing.md`
makes `python -O -m pytest` canonical *and states why it is safe*: pytest
AST-rewrites assertions inside test modules, so `-O` strips only production
(`orpheus/`) asserts. That reasoning is correct and it **does not transfer one
inch outside a pytest-collected test module**. A scratch probe, a `tools/`
script, a `__main__` block, and `tests/_harness/` helpers (which is exactly why
they already use explicit `raise`) all get silently disarmed.

**The rule.** A verification instrument that is not a collected pytest test must
either (a) run WITHOUT `-O`, or (b) use `if not cond: raise AssertionError(...)`.
Never `assert` in a standalone probe you intend to trust. And when a probe
reports a clean pass, ask whether its checks could have run at all — the
`vv-principles` #17 discipline (verify the instrument before trusting its
verdict) applies to the probe's *execution model*, not only to its logic.

⚠ The failure direction is the flattering one, as always: a disarmed probe never
reports a failure, so it looks like confirmation.

## L43 — A library message that names the CALLER's variable has crossed a boundary it cannot see across — and parameterising the name is the WRONG fix (2026-08-11)

**Context.** SN's truncation warning ended: *"... Or read
`solution.history.fully_converged` and handle it."* Perfectly helpful, shipped
for months, and **wrong in two independent ways**: it presumes the caller named
their variable `solution` (a fact no library can know), and #340 N4.7 was about
to make the same emitter serve CP, MoC and diffusion — whose entries return a
`*Result` and have no `solution` anywhere.

**The tempting fix, and why it is a trap.** Add a `verdict=` parameter so each
public entry passes the path on its own return type. It is one line per call
site and it reads as obviously correct. It is the **exact construction that step
N6a had just retired**: `budget_name` / `budget` / `tol` were also per-entry
facts passed in by the call site, and they were removed *because* a fact
asserted at a call site drifts from the object it describes — on a tree they
described the wrong level, every number real and every pairing wrong. A verdict
spelling is the same object with a different name, and it would have been
re-introduced into the very function that was cleaned of them.

**The rule.** When a diagnostic wants to tell the reader *how to interrogate the
result*, name the **attribute and its TYPE**, never a variable path — the type
is a fact the library owns, the variable name belongs to the caller. Here:
"read `fully_converged` on the `IterationRecord` this solve returned". True in
every family, no parameter, and each result type's own `record` docstring
already carries the last hop.

**Generalisation, and the grep for it.** Any message, docstring or error string
containing a dotted path whose FIRST segment is a plausible local-variable name
(`solution.`, `result.`, `solver.`, `mesh.`, `cfg.`) is a candidate. The
discriminator: *could this library have chosen that name?* If not, it is
guessing. Sibling of `feedback-lossy-return-type-is-the-root-cause` — both are
about a producer asserting something only the consumer can know.

## L44 — Three ways a verification instrument silently disarms itself, all in one step, all flattering (2026-08-11)

**Context.** #340 N4.7. Landing the carve needed three instruments: a
message-parity probe, an eight-site attribution probe, and a mutation battery.
**All three reported success while not measuring what they claimed**, and none
of the three failures announced itself.

| # | mechanism | what it printed | what was true |
|---|---|---|---|
| 1 | ran the standalone probe under `python -O` | `7/7 identical` + `both versions SILENT` | three `assert`s absent from the bytecode, incl. the entire silence claim (→ [[lessons-L42]]) |
| 2 | `except Exception: print(...); continue` in the site loop | `all sites: exactly one warning, attributed OUTSIDE orpheus/` | **3 of 7 sites raised in their fixture and were never measured** — the skip did not reach the failure list |
| 3 | mutation anchored on the bare token `"stacklevel=3"` | `41 passed` for both attribution mutants ⟹ "the gates are blind" | `.replace(…, 1)` hit the **docstring**, which documents `stacklevel=3` four times; the code was never mutated |

**The common shape, and it is not "be careful".** Each failure made the
instrument report the *reassuring* verdict: (1) a check that cannot fail, (2) a
gap that reads as coverage, (3) a blind gate that reads as a blind gate — which
is the nastiest, because #3 would have been *acted on*: the honest response to
"41 passed under a mutation" is to strengthen or delete the gate, and I would
have weakened a gate that was already correct.

**Three rules, each grep-checkable:**

1. **A standalone probe must not use `assert`** — `if not cond: raise
   AssertionError(...)`, or run without `-O`. (L42 has the project-specific why:
   the canonical `python -O -m pytest` is safe only because pytest rewrites
   assertions *inside collected test modules*.)
2. **A skipped case is a FAILURE, never a `continue`.** Any loop that measures
   N things and can bail on one must append to the failure list on the bail
   path. "Not measured" and "measured OK" must never print the same.
3. **A mutation anchor must match SYNTAX, not a token** — anchor on
   `"ConvergenceWarning,\n        stacklevel=3"`, not `"stacklevel=3"`. And note
   the aggravator: **a well-documented function is MORE likely to defeat a
   token anchor**, because the prose discusses the very constant the mutation
   targets. Asserting "the target string was present" does not help — it was
   present, in the docstring.

**The generalisation worth carrying past these three.** `vv-principles` #17
says verify the instrument on a known-positive before trusting a negative.
This step says the same thing about the instrument's *plumbing* rather than its
logic: **ask what this probe would print if it measured nothing at all.** If
that output is indistinguishable from success, the probe has no negative leg
and its green is worth nothing. In all three cases above the answer was "it
would print exactly what it just printed".

Sibling of [[lessons-L38]] (a mutation in a shared tree is indistinguishable
from a production bug) — both are about the *evidence pipeline* failing while
the code under test is fine.

## L45 — A DIAGNOSTIC can be annihilated by the very symmetry the design introduced (2026-08-11)

The instrument was right, the code was right, and the reading was worthless.

`[M]` Q5.6.4. BMC's contamination factor **β** is *the* diffusion-limit
diagnostic for an SN angular differencing scheme — zeroing it is the entire
reason the Morel–Montry τ exists. The in-tree evaluator existed and ran
clean. It reported **round-off for every candidate cell-edge convention**,
including one that provably **diverges the solve** (NaN from `n_φ ≥ 16`, τ
outside `[0,1]`). That contradiction is what forced a control:

| edge set fed to β | β |
|---|---|
| production | `+6.94e-18` |
| garbage — edges scaled 0.5× | `+3.47e-18` |
| garbage — edges **cubed** | `+1.73e-18` |
| garbage — **random**, antisymmetrised | `−3.47e-18` |
| one edge nudged (breaks antisymmetry) | **`−3.53e-03`** |

**β on a σ_y-folded arc is a symmetry identity, not a measurement.** The fold
makes the nodes antisymmetric and the α dome symmetric, so
`term_{M−1−m} = −term_m` and the sum cancels pairwise for ANY antisymmetric
edge set. β could only ever see antisymmetry.

⭐ **The generalisable shape, and it is nastier than ordinary Mode-12
blindness:** the blinding structure was introduced *by the campaign's own
central achievement*. The σ_y fold is the thing we built; it silently
annihilated the functional that would have judged what we built next. So the
question to ask after ANY symmetry-introducing carve is not "does my
diagnostic still run?" but **"is my diagnostic still a function of the thing I
am about to vary?"** — and the cheapest way to answer it is to feed the
diagnostic deliberate garbage in the varied slot (`vv-principles` #17,
extended from harnesses to *analytic* instruments).

⚠ Corollary for gates: a β gate on a folded cylinder rule would be green
forever while wearing an authoritative name. Either gate the UNFOLDED parent,
or gate the antisymmetry and say that is what you tested.

The replacement that DOES discriminate was in the same literature and needed
no solve — **ν-closure**: the cell march implied *by* τ must land on the
level's far endpoint. `[M]` `1.000000` for any derived τ, `1.016389` for the
retired clamp, `1.164784` for `τ ≡ ½`. Both of the latter correspond to no
partition of the level at all.

## L46 — TWO OBJECTS SHARING A LETTER, and the summary that contained both without noticing (2026-08-11)

`[M]` Q5.6.4. A literature memo returned, in the same report, both:

* *"τ was chosen to zero β"* (BMC Table II: β at round-off for the correct τ), and
* *"δ = 0 … would not give the correct diffusion limit"* where `δ = 2τ − 1`,
  and Lathrop's own recursion makes `δ ≡ 0 ⟹ β ≡ 0`.

Read together: the correct scheme both has and has not β = 0. I built a design
on the reading that fit my hypothesis, and it was the wrong one.

**They are two different quantities wearing one letter, and they are
near-OPPOSITES:** BMC's β is ONE SCALAR (the `J⁽²⁾` contamination), zero iff τ
is Morel–Montry; Lathrop's β is A SEQUENCE (α's pointwise defect), zero iff
`τ ≡ ½` — the *diamond* scheme. They also live at different orders: leading
(diffusion limit) vs first (truncation).

⭐ **The rule: when a summary of a source contains an internal tension, the
tension IS the finding — resolve it at the source before designing.** Do not
pick the branch that fits the hypothesis in flight; that is exactly when the
cost is highest, because the design is already leaning on it. The resolution
here also came with independent verification: reproducing BMC's Table I to
every printed digit proved which definition was live.

⚠ Same family as the three-way `tau` overload in this codebase (closure
weight / optical depth / critical half-thickness in mfp). When a symbol is
overloaded, **say which sense you mean at every use site**, and give the
nomenclature its own documented home — `derivations/discrete/sn/angular_differencing.py`
now carries it.

## L47 — A FIXTURE can be structurally incapable of seeing what its gate is credited with (2026-08-11)

`[M]` Q5.6.4. After re-posing the cylinder cell partition from η-midpoint to
ω-midpoint, `test_cell_visit_c_stamp`'s cylinder row stayed **green under
`assert_array_equal`** — bit-exact, no tolerance — while its oracle's τ
disagreed with production by `4.46e-2`.

Not a tolerance leak: its fixture is `folded_product(n_mu=2, n_phi=4)`, i.e.
**M = 2 ordinates per level**, and at M = 2 the interior chord midpoint
`(η₀+η₁)/2 = 0` **IS** the arc edge at `ω = π/2`. The two partitions are
**bit-identical** there, diverging only from M = 3 (`3.17e-02`, `4.46e-02`,
`6.82e-02` at M = 3/4/6).

⭐ **A green bit-exact gate is not evidence that its fixture can see the thing
you changed.** Before crediting (or blaming) a gate for a change, ask what the
SMALLEST fixture it runs is, and whether the change is even expressible there.
Degenerate small cases collapse distinctions — M = 2 has one interior edge, and
one interior edge is where every reasonable partition agrees.

⚠ The dual is the trap that bit me first: the same row went RED at an
intermediate state (oracle clamped, production unclamped) and GREEN after I
removed the clamp — which *reads* as "the fix worked" when the truth was "the
disagreement shrank below what this fixture can express". A red that turns
green is not proof of correctness; check whether the discriminating power was
ever there.

## L48 — Take the PREDICATE, not the recipe: how to port a literature scheme whose mechanism does not fit (2026-08-11)

User instruction, Q5.6.4, on being told our azimuthal rule is incompatible
with the literature's angular closure: *"Let's try to fix our own before
changing to any other implementation. We should be able to get this right, and
even if the literature suggests another method or closure, what we need might
not be their closure or method but the CONCEPT they used to make it accurate.
That is probably transferable."*

`[M]` It was exactly right, and the shape is reusable. The literature's
cell edges (cumulative quadrature weight in the radial cosine) are
**inapplicable** to our rule: ordinates land outside their own cells, worsening
with refinement (0/4 → 4/8 → 12/16 → 28/32), and the solve diverges. The
reason is one line: their construction presumes *the weight equals the cell's
measure in the marching variable*, and our equispaced-ω rule has equal weights
while an arc cell's η-measure ∝ `sin ω`.

⟹ **their edge recursion is one SOLUTION to a predicate system, not the
system.** Extracting the predicates (P0 α-closure, P1 `c = Σwη² = ⅔`
diffusion, P2 barycentric τ, P3 ordinate-inside-its-own-cell, P4 conservation)
made the port trivial: choose the partition that satisfies **P3** for OUR
nodes (the ω-midpoint), and **P2 then DETERMINES τ** with no freedom left —
`τ_m = ½ + ½·cot(ω_m)·tan(Δω/4)`, verified to `1.67e-16`.

⭐ **The procedure, when a published mechanism does not fit your
discretisation:** (1) find what CONDITION the mechanism was built to satisfy;
(2) check whether the condition is even stated separately from the mechanism —
usually it is, a few equations earlier; (3) evaluate the condition on your own
construction; (4) satisfy it your own way. Copying the mechanism is what
produces a "faithful" port that diverges; copying the condition is what
transfers.

⚠ And a bonus that fell out: **P3 became a THEOREM** on the new partition (on
a monotone arc the ω-midpoint edges bracket their own node, so `τ ∈ (0,1)` is
forced — `[M]` 4000 random arcs), whose only equality case is a node on Σ. So
cylinder-P3 reduces to the fold's own well-posedness criterion `Σ = ∅`. When a
port is right, the predicates tend to collapse into each other rather than
multiply.

---

## L49 — An instrument can be REFERENCE-limited, and the tell is that refining the SUT makes agreement WORSE (2026-08-11)

A cross-method comparison has TWO error sources, and the usual mental model —
"the reference is right, the SUT converges toward it" — silently assumes one of
them is zero. When it is not, the comparison has a **floor**, and past that
floor the number grades the reference, not the code.

**The tell is counter-intuitive and unmistakable once you look for it: refine
the SUT and the agreement gets WORSE.** That is the SUT converging *past* the
reference. Monotone improvement means you are still above the floor; monotone
degradation means you are below it and measuring the reference's own error.

> `[M]` 2026-08-11, Q5.6.4. The trajectory-resolvent flux-shape cross-check was
> the memo's capital-letters *"THE DECISIVE MEASUREMENT"* and the sole evidence
> for a `1.92× worse` verdict that nearly reverted a correct carve. Laddered in
> `n_φ`, all four τ conventions collapsed to `≈2.8–3.0e-2` by `n_φ = 32`, within
> 10 % of each other. Holding `n_φ = 32` and sweeping `nx` 40 → 80 → 160 then
> decided whose error that was: chord+absorber `3.03 → 3.23 → 3.37e-2`, arc
> `2.76 → 2.92 → 2.92e-2`, ω `2.76 → 2.97 → 3.10e-2` — **rising, monotonically,
> in every column.** The floor is the reference's own discretisation, and NO
> refinement on the SN side can push past it. Its dynamic range for an angular
> claim was exhausted by `n_φ = 16`; the only rung where it could rank anything
> was the coarsest, where its own error is 20–40 % of the signal it was grading.

**The procedure, and it is cheap.** Before crediting a cross-method comparison:
sweep the SUT's own resolution against the FIXED reference. Falls ⟹ still
above the floor, the number is about your code. Flat ⟹ at the floor. **Rises ⟹
below it, and every ranking you read there is noise plus reference error.**

⚠ **The consequence for the GATE, not just the investigation.** A
reference-limited gate cannot be tightened by improving the SUT — only by
refining the *reference* — and its tolerance is therefore a statement about the
reference that its docstring almost never makes. Ours carried
`tol_per_cell = 1.2e-1` against a `≈3e-2` floor while claiming to pin "the
resulting shape agreement". Say the bound, or the next person will read a pass
as evidence about the code.

Sibling of L47 (a fixture structurally unable to see what its gate is credited
with) and of the Mode-12 family: this is the same blindness located in the
REFERENCE rather than in the fixture or the functional.

---

## L50 — Two parallel sub-agents: one reads the other's UNCOMMITTED output as established fact (2026-08-12)

Dispatching agents in parallel is the standing posture (`delegation`), and it is
right. But agents share one working tree, and **an untracked file written by
agent A is indistinguishable, to agent B, from a landed artefact** — same path,
same content, same apparent authority. B has no way to ask "is this committed?"
unless told to, and the failure is silent and confident.

> `[M]` 2026-08-12, Q5.6.4. An `archivist` and a `test-architect` ran
> concurrently. The archivist had drafted a correct line — *"no ψ̂ positivity
> gate exists"* — then found `tests/sn/sweep/curvilinear/test_psi_half_positivity.py`,
> concluded *"the gate landed 2026-08-11"*, and **withdrew its own correct
> draft**. `git ls-files` says the path *"did not match any file(s) known to
> git"*: it was the test-architect's in-flight work, created minutes earlier in
> response to the main agent's brief, which had itself asserted the gate was
> missing. A closed loop of three, none of it checked.

⭐ The substance survived — the numbers in that module were right — which is
exactly what makes this dangerous: **a provenance error with correct content
leaves no symptom.** Had the test-architect's module been wrong or abandoned,
the archivist's report would have cited a file that never landed.

**The rules that follow:**

1. **A brief that states a NEGATIVE about the tree** ("X has no gate", "nothing
   asserts Y", "there is no consumer") must be re-verified by the recipient, and
   the brief should say so. It is the class of claim most likely to be stale by
   the time the agent reads it — including because of another agent you
   dispatched in the same turn.
2. **Before an agent credits a file as established, it checks tracked status.**
   `git ls-files --error-unmatch <path>` is one call. Anything untracked is a
   LIVE RESULT, never a landed gate, and must be reported as such.
3. **The main agent owns deconfliction**, because only it knows the full
   dispatch set. Scope parallel agents to disjoint trees (this session: the
   archivist was explicitly fenced out of `tests/`, which is the only reason the
   collision was a misreading rather than a merge conflict) — and when two
   agents *must* touch overlapping ground, tell each what the other is doing.

⚠ Generated artefacts are the second-order case: the `-E` Sphinx build
regenerated the V&V matrix and absorbed rows from the other agent's untracked
tests, so committing the docs pass alone would have pinned a matrix referencing
files git did not know about. **Regenerate once, after all parallel work lands.**

---

## L51 — An instrument has a TYPE: CONSTRAINT, RANKER, or DIAGNOSTIC. Using one as another is the error (2026-08-12)

**Context.** The SN #235 angular-closure campaign had killed five τ instruments
for blindness and needed a sixth. Two candidates arrived the same day and
*both* were misclassified — in opposite directions — by people (me included)
who had each verified the instrument was measuring *something* real.

**The taxonomy, and it is the whole lesson:**

| type | answers | may it rank designs? |
|---|---|---|
| **CONSTRAINT** | "does this scheme satisfy property P?" pass/fail | ⛔ never |
| **RANKER** | "which of these schemes is more accurate?" ordering | ✅ that is its job |
| **DIAGNOSTIC** | "how large is residual R?" a number | ⛔ never, until correlated |

**The two founding misclassifications:**

* `[M]` **#319's diffusion-limit flux-dip test was proposed as a RANKER; it is
  a CONSTRAINT, and one *structurally biased toward the incumbent*.** The
  diffusion limit's angular content is `span{1, μ}`, which is exactly the
  closure's kernel — the barycentric τ is exact there *by construction*, so
  the test scores the shipped scheme zero and every alternative positive, at
  every order, on every material. It is a perfectly good pass/fail question and
  can never order two schemes.
* `[M]` **the endpoint defect `D` was proposed as a RANKER; it is a
  DIAGNOSTIC.** Reference-free, pointwise, genuinely τ-loaded (ranks garbage
  2.6–45× above production) — and **uncorrelated with accuracy**: Pearson r on
  log = `+0.75 → +0.26 → +0.06` at `n_φ = 8/16/32` against an analytic MMS,
  with 0/4 ranks agreeing at the two finer orders. "It moves when τ moves" is
  necessary and nowhere near sufficient.

**Why the two errors look identical from inside.** In both cases the
instrument *is* sensitive to the parameter, the numbers *are* real, and the
ordering *is* stable. What fails is the inference from "sensitive to τ" to
"ranks τ by accuracy". Sensitivity is about the instrument; ranking is a claim
about the world.

⟹ **Declare the type BEFORE the first number, in the instrument's own
docstring.** It costs one line, it is unfalsifiable-by-drift (unlike a
tolerance), and it is the only thing standing between a constraint and a
design decision built on it. Both promoted gates from this campaign now carry
an explicit "what this CANNOT do" paragraph naming the type.

Companion to `vv-principles` #24 (which supplies the four checks an instrument
must pass to be a ranker at all) and to [[lessons-L49]] (reference-limited).

---

## L52 — A truncation-order statement about the EQUATION is not a ranking of solution ACCURACY (2026-08-12)

**The failure, mine, and it survived three hours and one confident report.**
Lathrop (2000) NSE 134 Eq. (30) proves that with `δ = 2τ−1` the curvilinear
angular truncation is `O(δΔμ + Δμ²)`, so *"only with `μ_m = μ̄` (`δ = 0`) is the
truncation order `O(Δμ²)`"* — every weighted diamond is first-order in angle.
`[M]` our shipped τ has `(τ−½)/w` unbounded (`= M/4`), so it is squarely in the
first-order class. I inferred: **therefore the shipped scheme is worse than
plain diamond**, and reported it with a measured 8× to back it.

⛔ **The inference is invalid, and the measurement that "confirmed" it came
from a fixture in the scheme's own kernel.** On a fixture built to excite the
closure, the shipped scheme **wins** by 34 % / 12 % at `n_φ = 16 / 32`.

**Why the two are not the same claim.** A truncation order bounds the *local
residual* the discrete equation leaves when the exact solution is substituted.
Solution accuracy is that residual propagated through the operator's inverse
and measured in some norm on some problem. Between them sit: the amplification
of the inverse, the other error channels the scheme also controls (here, the
diffusion limit), the norm, and the fixture. Any of the four can dominate, and
in this case two did.

⟹ **A published order statement is a hypothesis about accuracy, never a
measurement of it.** Cite it to *predict*, then measure — and per
`vv-principles` #24(d)/(e), measure on a fixture that is neither in the
scheme's kernel nor outside the regime where its advantage lives.

⚠ The generalisation that makes this worth a lesson rather than a footnote:
the literature statement was **true**, the numerical measurement was
**correct**, and the conclusion was still **false**. Two sound inputs, one
invalid join. When a literature result and a measurement agree, check that
they are about the same object before treating the agreement as corroboration.

## L53 — A stale INVENTORY is a silent denominator: 11 directories in neither column (2026-08-12)

Task #51's red inventory was carefully written, `[M]`-marked, and explicitly
honest about its limits: it named the two slices it HAD measured and listed the
seven directories it had NOT. Both lists were true. The defect was the **gap
between them** — `find tests -mindepth 1 -type d` returns 11 directories that
appear in *neither*, so they were not "unmeasured", they were **uncounted**.

`[M]` The uncounted set: `_harness`, `_mutation`, `cross_method`, `derivations`,
`homogeneous`, `numerics`, `transport`, `sn/acceleration`, `sn/architecture`,
`sn/l1_analytical`, `sn/regression`, `sn/solve`. One of them held **3 further
reds** (`sn/solve`, the #333 sha256 gate — and OPEN issue #333 had said so since
2026-08-06, in writing). The reported denominator was 20; the true one was 23.

⚠ **Why the honest caveat made it worse, not better.** A note saying "NOT
measured: A, B, C" reads as a complete statement of ignorance, so the next
session budgets for A, B, C and treats everything else as covered. An inventory
with no caveat at all would have prompted "what about the rest?"; this one
answered that question wrongly and authoritatively. Same family as
`plan-authoring` §2's quantifier rule (a universal claim carries its
denominator), one level up: here the *partition itself* was incomplete, and
neither list was individually false.

⟹ **The check is mechanical and costs one command: enumerate the universe, then
subtract both lists.** For a test inventory that is `find tests -type d`; for a
call-site audit it is the full `grep`, not the union of the files you happened
to open. Write the leftover set explicitly, even (especially) when it is empty —
"measured A+B, not-measured C, leftover ∅" is a *different and stronger* claim
than "measured A+B, not-measured C", and only the first one can be checked.

**Corollary — reconcile against the issue tracker, not only against the tree.**
#333 named its own red arms and its own blocker ("Blocked on #327", long since
closed at `414f2cb6`). One `gh issue list` cross-read would have surfaced the
missing rows before any measurement. Cardinal Rule 4 says issues ARE the log;
that cuts both ways — an inventory that does not consult them is re-deriving,
badly, what the log already knows.

`[M]` Closed out: 3403 + 2194 + 3446 passed across the previously-uncounted and
re-verified slices, 0 reds. Every one of the 23 was PRE-EXISTING (byte-identical
failure sets at detached worktree `adb73fd5`); none was a regression.

## L54 — An instrument that never RAN reports in the safe-looking direction; verify execution before reading output (2026-08-15)

L49/L51 and `vv-principles` #24 all concern an instrument that **runs** and
measures the wrong thing. This is the prior failure: the command **never
executed at all**, and what the harness printed *read like a result*. Three
instances in one session, three different mechanisms, all pointing the
flattering way:

1. **zsh does not word-split an unquoted `$VAR`.** `T="a.py b.py c.py"; pytest
   $T` passes ONE argument — a path that does not exist. pytest printed
   `1 warning in 0.01s` and exited 0. That reads as a fast clean run; **zero
   tests were collected.** (Bash *would* have split it, so the idiom is
   correct-looking and shell-dependent. This is the same defect the SN #344
   campaign had already logged once, where a quoted `"$T"` collapsed three test
   paths into one.)
2. **`pytest -p <path>` silently loads nothing.** `-p` takes an importable
   MODULE NAME; given a filesystem path it does not error. The mutation plugin
   never installed, the grep matched nothing, and the empty output was
   indistinguishable from "the mutation reddened nothing" — i.e. from a real
   verdict of *blind*. Fix: `PYTHONPATH=<dir> ... -p <module_name>`, and have
   the plugin **print a banner** so its absence is visible.
3. **A `pgrep -f <pattern>` guard inside a command whose own command line
   contains `<pattern>` matches ITSELF**, so the condition never becomes false
   and an `until` loop hangs forever. A hang reads as "still running", which is
   the one status nobody investigates.

**The unifying tell: the output was EMPTY or trivially short, and empty was
interpreted as a measurement.** Absence of failure is not evidence of success;
absence of *anything* is not evidence at all.

**How to apply — three checks, each one line:**

- **Assert the denominator.** Every pytest invocation whose result you will
  quote must have its **collected count** read, not just its exit status. `N
  passed` where you expected ~83 is a red flag; `1 warning in 0.01s` is a dead
  run. (`vv-principles` #17's "check the collected count", generalised past
  mutation batteries.)
- **Make the instrument announce itself.** A mutation plugin, a monitor, a
  census hook — each prints one banner line on activation. Then "no banner" and
  "no findings" are different observations.
- **A guard must be able to become false.** Before an `until`/`while` loop,
  ask what event terminates it, and check the guard cannot be satisfied by the
  loop's own process. Prefer waiting on a **PID** (`ps -p $PID`) over a
  `pgrep -f` pattern match.

Cross-reference: `vv-principles` #17 (the harness lies before the code does, and
in the safe-looking direction) — this is that principle applied one step
earlier, to whether the harness ran; `[[lessons-L12]]` (paste-back — a verbatim
summary line makes a zero-collection run visible where a paraphrase hides it).

## L55 — A census is only as wide as its SPELLING; for a code construct, use a parser (2026-08-15)

[[lessons-L54]] is about an instrument that never RAN. This is the next
failure along: the instrument runs, matches, and reports a confident
number — while structurally unable to see one of the forms it is
censusing. **Three instances in one session, all mine, all reporting the
reassuring answer.**

| # | the census | what it could not see | cost |
|---|---|---|---|
| 1 | `grep 'app\.config\.nexus_'` to prove "all config reads are routed" | `getattr(app.config, "nexus_infer_implements", True)` — the name is a **string** there, not an attribute access | 2 settings stayed unrouted while the census printed clean; a config key shipped **inert** |
| 2 | exact-string `replace_all` of a multi-line `"--db", type=Path, default=Path(...),` | six **single-line** `add_argument("--db", type=Path, default=Path(...))` — no trailing comma, closing paren instead | every `runtime-*` verb ignored the config and returned `[]`, which reads exactly like "nothing ingested" |
| 3 | `ls -d ~/a ~/b /x/*/c` to locate a checkout | zsh `nomatch` **aborts the whole command line** on the failed glob, so the `ls` never ran | reported "no local checkout" of a repo that was right there |

⭐ **The rule, and it is mechanical: to census a CODE CONSTRUCT, walk the
AST, not the text.** An `ast.walk` for `add_argument("--db", …)` sees the
call however it is line-broken, spaced, or aliased; a grep sees only the
spellings someone thought of. The gate that now pins case 2 is exactly
that, and it is mutation-verified — re-introducing the missed single-line
form reddens it.

⚠ And the reason this is not merely "grep harder": **the failure direction
is always flattering.** A census that misses a form reports FEWER
offenders, i.e. "you are done". Same shape as
[[lessons-L44]]'s three disarmed instruments and `vv-principles` #17 — ask
of any census *"what would this print if it matched nothing, and is that
distinguishable from success?"*

⭐ Two more from the same session, worth carrying because both return
EMPTY rather than erroring:

- **SQLite double quotes bind to COLUMNS, not strings.** `where a.key="source"`
  joined against a table that HAS a `source` column resolves to that column
  and returns zero rows. Single-quote every string literal. (Cost: one false
  "this edge type has no provenance" reading.)
- **`grep -c 'X'` proving a retirement is clean** is the same defect as case 1
  whenever `X` has a second spelling — `.claude/rules/coding-standards.md`
  already carries the three-spellings-of-a-math-symbol rule; this is its
  code-construct twin.

### Two more, 2026-08-16 — the census was *right* and still incomplete

Cases 1–3 all missed a SPELLING. These two missed nothing of the sort; the
pattern was correct and the answer was still partial, which is why they need
separate naming.

- ⭐⭐ **A census keyed to one member of a FAMILY misses the family.** Moving
  the nexus graph, I censused `grep -rln "graph\.db"` and adjudicated every
  hit. What moved was the **directory**, which holds `graph.db`, `graph.json`
  AND `graph.html` — and `docs/index.rst` linked the explorer at
  `<_nexus/graph.html>`, so the shipped homepage rendered a **404**. Nothing
  catches it: a raw relative hyperlink is checked by Sphinx at *no* severity
  (unlike `:doc:`/`:ref:`). ⟹ when a **container** moves, census the
  container, not the one file you were thinking about: anchor on the
  separator — `grep -rnE "_nexus/|graph\.(db|json|html)"`. (⚠ and the bare
  container name was useless here: `_nexus` matches every `mcp__nexus__*`
  tool name, 559 KB of output — which is *why* the narrower `graph.db`
  looked like the sensible choice.)
- ⭐ **`head -N` on your own census turns a complete result into a confident
  partial one.** I piped a `grep -rn` through `head -20`, read the 20 lines
  as the answer, and concluded a file was clean that has **8 hits** — it sat
  below the cut, in traversal order, not alphabetical. ⟹ census with
  `grep -rln` (file list, no bodies) or pipe through `wc -l` first; reserve
  `head` for *reading* hits, never for *counting* them.

Cross-reference: [[lessons-L54]] (the instrument never ran — this is the
instrument running half-blind); [[lessons-L25]] (`replace_all` is safe only
if every occurrence is the target concept — case 2 is its dual: safe on what
it matched, and it did not match everything); [[lessons-L56]] (what the
half-blind census then FAILS to notice, because the broken consumer went
quiet instead of loud).

---

## L56 — "nothing found" and "I looked in the wrong place" must not print the same thing (2026-08-16)

[[lessons-L55]] is about a search that cannot see part of what it looks
for. This is its consumer-side twin: the search is fine, the **reader of
its result** cannot tell an empty answer from a broken one — so a
misconfiguration is delivered as a fact.

**Four instances in one session, three of them independent code paths,
all failing in the reassuring direction.**

| the consumer | what it printed when MISCONFIGURED | what that is indistinguishable from |
|---|---|---|
| `runtime-ingest --kind coverage` | `nodes: 0 / edges: 0 / unresolved: 0`, **exit 0** | a workload that genuinely touched nothing indexed |
| `nexus-brief.sh`, `nexus-dead-refs.sh` | nothing at all (`[ -f "$db" ] \|\| exit 0`) | a project with no graph built |
| `/doc-health` | `(nexus or graph not found — build the docs)` | docs never built |
| `nexus status` from outside the project | `does not exist / run 'nexus analyze' first` | a graph that was never built |

The shared mechanic: each has a **legitimate empty case**, and the
failure path was routed into it. A hook whose contract is "quiet exit 0
when there is no graph" cannot express "I looked in the wrong place" —
so when `[graph].output` moved, two of three hooks went **dark rather
than red**, and the breakage presented as an absence for as long as
nobody wondered why the briefs had stopped.

⭐ **The rule, and it is one line of output: a failure must name the
thing it looked for.** Every message above became actionable by adding
the resolved path — `no graph at <path> (the path .nexus/config.toml
declares)`. That is what makes a stale path *visible* instead of
inferred; the old messages named no path at all, which is precisely why
a wrong one could survive a campaign.

⭐⭐ **And the design rule underneath it: a count of drops is not enough —
tally them BY REASON.** `unresolved: 0` was true and useless, because the
coverage records were discarded *upstream* of that counter. The repair
was a `JoinLedger` with `outside_scope` / `unindexed_file` /
`no_enclosing_node`, since those three have three different remedies and
one number collapses them. `[M]` the same artifact went 0 → **1691** bound
nodes; the ledger then reads `bound: 121 / unindexed: 37 / outside: 0`
over 158 files, and the wrong-root case reads `unindexed: 158` — a
signature you can diagnose on sight.

⚠ **The asymmetry that hides this class, and the reason a green
neighbour proves nothing:** cProfile's `co_filename` and viztracer's
event names are **absolute**, so neither backend ever hit the relative-key
path. A clean cProfile ingest was not evidence that coverage ingest
worked — exactly `vv-principles` #19's "a positive reading cannot
discriminate loaded from blind", applied to a sibling code path instead
of to a gate.

⟹ **The question to ask of any consumer that can return empty:** *what
does this print when its input is misconfigured, and is that
distinguishable from success?* If not, it is not an instrument. This is
also the standing form of [[lessons-L54]] (the instrument that never ran)
— there the run was missing, here the *reason* is.

Cross-reference: nexus #59 (the same defect surveyed across `callers` and
the whole `runtime-*` family), nexus #56 (fixed, `5796c6d`).

## L57 — Of two twins, the bug lives in the one whose wrong answer is INVISIBLE; so a defect report naming a helper is first a question about DUPLICATION (2026-08-16)

A review reported nexus's `_unparse_attribute` as fabricating `calls` edges
and priced the repair at **one line** (make it return `None` on an
unresolvable root). The mechanism was right and the framing was wrong: the
function should not have existed. `_dotted_name` — **360 lines above, in the
same file** — is the same computation and *already* returned `None` on an
unresolvable root. Landing the reported one-line fix would have left two
implementations of "reconstruct a dotted name", which is the condition that
produced the bug. The correct repair was a **retirement**, and it removed 28
production lines while dissolving the `Name`-vs-`Attribute` branch at the call
site as a by-product.

⭐ **The transferable part is WHICH twin carried the defect, because it was
predictable.** The two had different consumers, and the consumers had very
different feedback loops:

| twin | consumer | a truncated answer is… |
|---|---|---|
| `_dotted_name` | decorator/marker parsing | **visible** — `@pytest.mark.l0` stops being recognised |
| `_unparse_attribute` | call-target resolution | **invisible** — it mints a plausible edge nobody checks |

The duplicate serving the *low-feedback* consumer is where the bug survives,
and it survives *for as long as the duplication does*. So the twins do not
merely risk drifting — the drift is **biased**, toward whichever copy nobody
can see is wrong.

⟹ Two rules, both cheap:
1. **Before fixing a named helper, grep its own file for a function that does
   the same job.** One grep. A defect in a helper is evidence about the
   helper's *uniqueness* before it is evidence about its logic.
2. **When you do find twins, suspect the one with the weakest feedback loop
   first**, and fix by retiring toward the one whose consumer would have
   screamed — it is the one with the better-tested contract.

⚠ **And the direction of failure decides the priority.** This one failed
false-**ALIVE**: it invented callers. `[M]` on ORPHEUS's graph — 1741 `calls`
pairs removed and **0 added**, 291 of them pointing at a real indexed symbol,
57 → 5 fabricated self-loops (the graph claimed
`SumOfTensorProductsOperator.apply` calls itself), 94 real symbols whose
incoming call edges were **all** fabricated, `dead_functions` candidates
2919 → **2967**. A false-alive defect is strictly worse than a gap: a gap
announces itself when you go looking, while fabricated evidence *answers the
question* — and answers it plausibly. Note those 94 are **not** thereby dead;
most are protocol-typed methods reached by dynamic dispatch, which the static
graph cannot see either way. What changed is that it stopped presenting
invention as evidence.

⭐ **The same session produced the same shape a second time**, which is why
this is a lesson and not a note: the role→objtype map (`func` → `function`)
was a **local dict** inside the docstring scanner, so the doctree walker built
ids from the raw role. `[M]` 316 short-prefix nodes and **265 symbols carrying
both spellings at once**, each holding part of the symbol's edges. Same cause
(one concept, two homes), same invisibility (a split id resolves to *something*),
same repair (hoist to one home — `_mappings.REFTYPE_OBJTYPE_MAP` — and let both
producers read it). Two instances in one file family in one afternoon is a
statement about the module, not about luck.

⭐⭐ **A THIRD and FOURTH instance the same day, and the fourth carries the
sharper rule.** #68's equation-namespace leak was again two sibling branches in
ONE function with unequal discipline: the `:math:` branch guarded by a
**blocklist** (reject `\`, `{`, `}`) three lines above a Python branch asking
the opposite, stronger question (`_is_dotted_identifier` — *is this a
well-formed name?*). Inline math clears a blocklist trivially, so `math:equation:*`
was **51 % LaTeX** (`[M]` 956 of 1860 unresolved). And its newline half was the
doctree walker missing the wrap-normalisation the docstring scanner had always
had — the same one-of-two-producers split as the objtype map.

⚠ **The fourth instance also shows the cost of fixing one by SHARING**, which is
the right repair and is not free: hoisting `_normalize_wrapped_target` so both
producers read it exposed that the function was *already wrong*. It collapsed
ALL whitespace whenever the result matched the dotted shape — and a sentence of
letters plus a full stop matches once its spaces are gone. Downstream,
`_classify_unresolved` refuses a non-identifier as napoleon noise, so prose was
already being dropped correctly; normalising first **disguised prose as an
identifier** and walked it through that gate. `[M]` 48 junk classes minted
(`py:class:allkeyvariables.`, `py:class:dictmappingmaterialIDtoMixture.`).

⟹ Two rules fall out, both cheap:
3. **Never normalise before classifying.** A normaliser makes strings *look
   more like* what a classifier accepts, so running it first can only ever
   widen the classifier's acceptance — never narrow it. If a gate downstream
   rejects malformed input, normalisation belongs *after* it.
4. **When a pass is meant to be purely subtractive, check that its node delta
   IS purely subtractive.** This one was found by arithmetic, not by a test:
   957 nodes dropped, net **−910** — so 48 appeared, and a pass that only
   removes things cannot add any. One query. The suite was green throughout,
   and no gate I had written would have caught it, because I had not thought to
   assert a *direction*.

⭐⭐ **A FIFTH instance, and it sharpens the principle past twins: the SAME
code path survives a defect for as long as its output is unobservable — so
turning a READ into a WRITE is a diagnostic move.** The four cases above are
two copies with unequal feedback. This one has no copy. `find_project_root`
walked up looking for `.nexus/`, and `~/.nexus/` exists on every machine that
has run the MCP server, because that is where `usage.jsonl` lives. So **every
project without a config of its own resolved its root to `$HOME`** and read
the wrong project's settings — for months, silently, because a wrong answer
about "which project am I in" still *looks* like an answer. Moving the graph
store into `.nexus/` turned that same resolution into a file write, and
`graph.db` appeared in the home directory within one test run. `[M]`
`find_project_root("<repo>/tests/roots")` → `/Users/rodrigo`.

⟹ **When a resolver has been trusted for a long time and never checked, ask
what would make its answer VISIBLE**, and arrange that once. A read is
checked only by someone who already suspects it; a write checks itself.

⚠ And the root cause is worth naming on its own, because it is a naming
defect and it recurs: **one directory name carried two scopes** — `.nexus/`
as a *project's* settings and `.nexus/` as the *machine's* state. No walk can
tell those apart, so the search had to be BOUNDED (at `$HOME`, and at a
checkout) rather than made smarter. Sibling of [[lessons-L46]] (two objects
sharing a letter): when one name spans two scopes, the fix is a boundary, not
a better lookup.

Cross-reference: `coding-elegance` Pattern 2 (single source of truth) — this is
its *diagnostic* direction, running from a bug report back to the duplication
rather than from duplication forward to a predicted bug; `[[lessons-L56]]`
(a defect that reports in the reassuring direction); `.claude/rules/coding-standards.md`
(retirement as a first-class deliverable). Landed nexus `3e137ff`, the
follow-up spelling commit, and `c51672c` (the store move + the two bounds).

---

## L58 — A promise that code did NOT do something cannot be gated by what the code returns (2026-08-16)

Nexus Track 1.1 added a staleness pass to the MCP tool boundary, promising it
would be **absent to the byte** when the graph is fresh: check a cached set,
and if nothing changed hand back the *original* payload string rather than a
re-serialisation. Two gates asserted that promise — one by object identity,
one by the flag's absence — and both were green.

`[M]` The mutation battery came back **GREEN on the two arms those gates
exist to cover**:

| mutation | why the gate could not see it |
|---|---|
| delete the fresh-graph early return | the walk then runs, marks nothing, and returns the same object anyway — identity still holds |
| re-serialise even when nothing was marked | unreachable from a fresh graph: the early return preempts it |

Both gates were about a **route**, and both asserted an **output**. A
correct implementation and one that does the work and throws it away are
indistinguishable in the return value, so the gate is green under either —
wearing a name that says the promise is covered.

⭐ **The general shape, which is broader than caching.** Any promise of the
form *"this does not do X when Y"* — a short-circuit, a memo, a fast lane, a
no-op guard, "costs nothing when there is nothing to do" — is a claim about
what the code SKIPPED. Skipping leaves no trace in the result. That is the
whole point of skipping.

⟹ **Two moves, both cheap:**

1. **Instrument the route, not the output.** Monkeypatch the expensive call
   *in the module's own namespace* and assert the call list is empty. Pick
   the cheapest observable on the skipped path — the parse, the subprocess,
   the `stat` — never a timing measurement, which is a flaky proxy for the
   same question. Here: a `json.loads` counter, and `parses == []`.
2. **Gate BOTH sides of the branch.** The case that exercises a fast lane's
   *fallback* is not the case the promise is about. The re-serialisation
   mutation only became visible from a **dirty tree with an unaffected
   file** — the one configuration that reaches the walk and marks nothing.
   A single gate on the happy side can never reach the other arm.

⚠ **And assert identity, not equality, whenever the promise is "unchanged".**
`to_json(json.loads(x)) == x` for anything `to_json` produced, so equality
cannot see a needless round-trip. `is` can.

Sibling of `[[lessons-L54]]` (an instrument that never RAN reports in the safe
direction) — there the instrument was dead, here it was alive and pointed at
the wrong thing, and both fail green. Promoted to `vv-principles` **#26**,
which is where the review-time form lives; it is #19's discrimination rule
(a positive reading cannot separate loaded from blind) moved from a gate's
*sensitivity* to a function's *route*. Landed nexus `8de24ec`.

## L59 — A fixture is blind on the axis where its author had only ONE example in front of them, and that is predictable BEFORE the mutation (2026-08-16)

`vv-principles` #24(d)/(e) says the fixture can be the lie, and #17's
per-arm battery is how you find out. This is the sharpening the graph
campaign paid for **four times in three days**: you do not have to wait
for the battery, because *which* axis a fixture cannot see is
predictable. It is the axis the author never had to think about — the one
where every example available to them looked the same.

| the blind fixture | the axis it could not see | what its author had in front of them |
|---|---|---|
| `rich_graph` spelling `file_path` as `"alpha.py"` | absolute vs relative stored paths | one graph, hand-written, relative |
| the byte-identity done-when | the fresh path (an edited file cannot reach it) | a dirty tree, because that is what the feature is for |
| every graph fixture in `test_runtime_ingest` | absolute vs relative stored paths, **again** | Sphinx builds, which are always absolute |
| the decorator gate's defs at 10 / 30 / 52 | siblings within `DECORATOR_WINDOW = 8` | a fixture written for readability, where you space things out |

`[M]` the last two were arms **6 of 6** of a battery whose other four
reddened cleanly; both mutations were **no-ops against the fixture**, so
both reported GREEN and would have read as "well covered".

⭐ **The rule, and it is one question asked at fixture-writing time:**
*what did I have only one of when I wrote this?* One path spelling, one
graph state, one parity, one ordering, one file layout, one seed. That
is the axis the gate is about to certify without testing. Either add the
second example or write in the docstring that the axis is unexercised —
what must not happen is a green battery arm being read as coverage.

⚠ The tell that makes it findable in review, before any mutation: **the
fixture is more regular than the world.** Real `@property` blocks sit 5
lines apart; the fixture spaced them 20. Real graphs mix absolute and
relative; the fixture picked one. Regularity is what makes a fixture
readable and it is exactly what makes it blind — so the readability that
recommends a fixture is a *warning sign* about it, not a defence.

Cross-reference: `vv-principles` #24(d)/(e) (configuration includes what
the fixture cannot see), #17 (per-arm granularity, which is how these
surfaced), #13's refinement-ladder disguise — `8/16/32/64` is the same
defect wearing arithmetic, a single congruence class the author never
had a reason to break. Sibling of [[lessons-L57]]: there the *defect*
hides where the wrong answer is invisible, here the *gap* hides where the
fixture is uniform. Landed nexus `7db466d` / `c497ddb`.

---

## L60 — A type SIGNATURE is evidence about its author's assumptions, never about a theorem (2026-08-26)

**The setup.** An architecture decision hung on one question: does the
Morel–Montry angular weight `τ` become cell-dependent under a
linear-discontinuous spatial scheme? If yes, the shipped
`tau_per_ordinate` accessor's arity had to widen before the
`AngularClosure` ABC could be frozen.

**What I argued.** I looked at the three functions that compute the
flux-dip contamination coefficient and its weight —
`contamination_beta(quad, geometry, *, edges)`,
`morel_montry_tau_per_level(quad, coord)`,
`angular_cell_edges_per_level(quad, coord)` — observed that **not one of
them can accept a spatial argument**, and offered that as evidence that
`β` (and therefore `τ`) is spatially independent. I also measured `τ`
bit-identical across five meshes including a graded one, which is real
evidence about the *shipped implementation*.

**Why the argument is invalid** (caught by the derivation agent I had
briefed to attack it; the *conclusion* turned out true, via a different
route). The signature has no spatial argument **because Morel–Montry and
Bailey–Morel–Chang derived `β` with space held CONTINUOUS, and say so**.
The absence is a faithful record of the derivation's scope. Had the
joint (angular × LD-spatial) problem introduced spatial dependence, the
answer would have been a *new function carrying a spatial argument* —
the existing one's arity says exactly nothing about that. The argument
assumed the corpus already contained the answer to a question nobody had
yet asked.

⟹ **`vv-principles` Mode 8, in a new dress**: the reading cannot come
out any other way, so it carries no information. A signature is a
*tautological* instrument for any question about mathematics — it can
only ever report what its author already believed.

**What actually settled it** (and the shape to reach for first): a
theorem with checked hypotheses. A scalar convex combination commutes
with every linear map, so `τ`'s two defining conditions — cone
membership `τ ∈ [0,1]` and the barycentric condition — are the *same
scalar statement in every moment component*; a per-(ordinate, cell) `τ`
is an overdetermined system whose every row returns the same value.
Hypotheses: `τ` is `r`-independent by construction, the moment
projection is linear, the flux cone is convex. **None of them mentions a
basis** — which is the whole content. Corroborated by the redistribution
operator factoring as a TENSOR PRODUCT `R_spatial ⊗ A_angular`, whose
angular factor's free symbols are `{μ, w, τ}`.

**The generalisable rule.** When a question is about *mathematics*, rank
the evidence: a theorem with checked hypotheses ≫ a measurement of the
shipped code ≫ the shape of an API. The last of those feels like
evidence because it is cheap, checkable, and phrased in the same
vocabulary as the question — which is precisely what makes it seductive.
Before citing an interface as evidence, ask: **"what would this
interface look like if the answer were the opposite?"** If the honest
answer is "exactly the same, until somebody changed it", you have a
tautology.

⚠ And the aggravator specific to this case: the invalid argument sat
*beside* a genuine measurement (`τ` bit-identical across five meshes),
and the measurement's authority carried the argument. That is
`plan-authoring` §2's `[M]`-scope defect wearing an API instead of a
conjunction — a right number under one proposition lending its badge to
a neighbouring claim it never tested.

Cross-reference: `vv-principles` #17 (verify the instrument on a known
positive), #19 (a positive reading cannot discriminate loaded from
blind), Mode 12 (ask what the functional annihilates — here, what the
signature *could* record). Sibling of [[lessons-L58]]: there a promise
about what code did NOT do could not be gated by what it returns; here a
claim about what mathematics does NOT depend on could not be gated by
what an API accepts. Charter §5d.8; memo
`scratch/tau_under_ld_dip_analysis.md` §5.

---

## L61 — An unvalidated FILTER and a clean tree print the same thing; and in a commit message, the shell's own quoting silently deletes evidence (2026-08-26)

**Six** incidents in one session, **two distinct mechanisms**, one shared
shape: **a pattern that never matched what it was written to match, and output
that looked like a clean result.** Mechanism (c) below is the worse one — it
prints no error at all.

**(a) The census that returned a false zero.** Auditing P3's rename targets for
string-form references — `getattr(x, "name", default)`, the form that fails
*silently in the default's direction* and had already broken two tests during
P1's `curvature` retirement — I wrote

```zsh
grep -rnE "['\"]$sym['\"]" --include='*.py' orpheus/ tests/ | wc -l
```

zsh could not parse the embedded quotes, printed `(eval): bad math expression`,
and `wc -l` faithfully reported **`0`**. Twice. A broken filter and a
consumer-free symbol are **indistinguishable in the output** — and the failure
direction is flattering: it says *"nothing to do"*.

⟹ **Re-run any nontrivial pattern in Python, and prove the filter against a
POSITIVE CONTROL before trusting a negative.** The third attempt asserted a
known member first (`MorelMontryAngularSweep` → 2 known `__all__` hits) and
refused to report if the control came back empty. That one line is the whole
difference between a measurement and a ritual.

**(b) The same defect eating a COMMIT MESSAGE — and this one is worse, because
it damages a permanent artifact and nothing warns.** A commit written as
`git commit -m "… ChartConnection's \`angular\` field needs a type …"` has its
backticks **command-substituted** by the shell. zsh printed
`command not found: angular` and committed anyway; the recorded message read
*"ChartConnection's  field needs a type"* — the identifying word **deleted**,
the sentence still grammatical, the commit still green.

This project's prose style puts backticks in nearly every commit message, so
the exposure is not occasional — it is every commit.

⟹ **Never pass a commit message through `-m "…"`. Use `-F -` with a QUOTED
heredoc** (`<<'MSG'` — the quotes on the delimiter are what disable
substitution). And after any message that was not written that way, run
`git log -1 --format=%B` and read it: the damage is invisible at commit time and
permanent afterwards.

**(c) — and the one that emits NO diagnostic whatsoever: `grep` here is
`ugrep`, and an anchor inside an alternation group silently matches nothing.**
`[M]` `grep` is a shell function wrapping **ugrep 7.5.0**, not GNU/BSD grep.
Isolated on a fixture containing `square Gram, while`:

| pattern | matches |
|---|---|
| `grep -E 'Gram'` | 1 ✅ |
| `grep -E '(^\|[^a-z_])Gram'` | **0** ⛔ exit 1, **no error** |
| `grep -E '([^a-z_])Gram'` | 1 ✅ |
| `grep -P '(?<![A-Za-z_])Gram'` | 1 ✅ |
| `grep -E '\bGram\b'` | 1 ✅ |

⚠ **That construction is precisely the retirement-audit idiom** — *"the symbol,
but not preceded by a letter"* is how you separate `gram` from `pro`**`gram`**`s`,
or a retired module name from a surviving attribute spelled identically. It hit
me twice: once checking whether the `gram → pairing` keeps had survived (shell
said 0, Python found the 6 that were correctly kept), and once on a residual
sweep. Both read as *"clean, nothing left"*.

⟹ Recorded in `.claude/rules/nexus-tools.md` with the safe forms
(`\b…\b`, `-P` lookbehind, or no anchor). **For any COMPLETENESS claim, use
Python.**

⭐ **What unifies (a), (b) and (c), and why this is one lesson rather than three:** in
all three, the failure surfaced *on a channel nobody was reading, or on none at
all* — an `(eval)` warning above a plausible number, a `command not found` above
a successful commit, and in (c) **nothing whatsoever** — while the artifact
looked correct. This is
`vv-principles` #17's broken-harness failure (*the instrument lies before the
code does, and it lies in the safe-looking direction*) applied to the shell
itself, and `plan-authoring` §2's FILTER clause sharpened once more: writing the
pattern down is not validating it, and **a filter's own error output is not
where you will look.**

⚠ The aggravator specific to this session: I had *already* written the FILTER
clause's newest sharpening into `plan-authoring` hours earlier, in this same
session, about this same campaign. The habit did not transfer to my own shell
commands — which is the same non-transfer recorded at `coding-standards`'
string-form clause (P1 ran the check for `mu_start`, then failed to run it for
`curvature` one commit later, the same afternoon).

Cross-reference: `plan-authoring` §2 (the FILTER and VIEWPORT clauses — never
`head` an enumeration, never trust an unvalidated pattern), `vv-principles` #17
(validate the instrument on a known positive), `coding-standards` (a symbol grep
cannot see a name inside a string). Sibling of [[lessons-L56]] — there,
"nothing found" and "looked in the wrong place" printed the same thing; here,
"pattern matched nothing" and "pattern never ran" do.

## L62 — A positive control drawn from the TREE validates a filter only for shapes the tree already exhibits (2026-08-27)

`vv-principles` #17 and [[lessons-L61]] establish the rule: validate a filter
against a positive control before trusting a negative. This is the boundary of
that rule, and it was found by a census that ran the control correctly and was
still wrong twice.

**Context.** Censusing whether any non-SN solver computes `delta_A`
(`= A_{i+1/2} − A_{i−1/2}`, the contracted connection coefficient) under any
disguise. Ten filter families over 334 files, each validated before its zero was
read. Tree control: family 1 must find the two known sites in
`geometry/reduced_operator.py`. It did.

**Two defects survived anyway, and they fail differently:**

| # | the filter | why the control did not catch it |
|---|---|---|
| 1 | the slice-difference `X[1:] − X[:-1]` was **single-axis only** | the known site is 1-D, so the control had nothing multi-axis in it. It silently could not see `phi[:, 1:] − phi[:, :-1]`. |
| 2 | the closed-form family required an **explicit leading coefficient** | ⭐ **no control was possible**: the known site is a *slice*-difference, so the closed-form family had **no tree member to validate against at all**. |

Defect 2 is the lesson. It missed `np.pi * (r_outer**2 - r_inner**2)` in
`moc/geometry.py:284` — the single most important near-miss in the whole census,
literally a difference of two areas and numerically `ΔA_sph/4`. It was found by
**reading the file**, not by any filter.

> **The rule: a positive control drawn from the tree can only validate a filter
> for shapes the tree already exhibits. Where a filter family has no exemplar in
> the tree, its control is vacuous and reading is the only instrument.**

⟹ Two checks, both cheap:
1. **Validate against a LITERAL of the shape, not against a tree member.** An
   in-memory fixture of the exact form you claim to catch works even when the
   tree contains no instance. (The census did this for eight of ten families;
   the two that failed are the two that used the tree as their control.)
2. **Enumerate your filter families and ask which have no tree exemplar.** Those
   are exactly the ones whose zero carries no information, and they must be
   closed by reading rather than by grepping.

⚠ The failure direction is the usual flattering one: a family with no exemplar
reports zero everywhere, which reads as "the concept does not occur" rather than
"I cannot see this shape".

⭐ And the adjudication that made the near-miss decidable is worth carrying on
its own: when two quantities have the same arithmetic form, discriminate by the
**dimension of the measure being differenced**. `delta_A` differences a
*(d−1)*-dimensional SURFACE measure (`2πr`, `4πr²`) — a connection coefficient.
`πr²` is a *d*-dimensional VOLUME measure (in MoC's 2-D plane the volume measure
IS an area) — a cell volume. The decisive witness was a shipped gate asserting
`region_areas.sum() == pitch**2`: the regions PARTITION the domain, and a sum of
surface differences telescopes to `A_out − A_in`, never to a domain measure.

Cross-reference: [[lessons-L61]] (an unvalidated filter and a clean tree print
the same thing — this is its boundary case, where the filter IS validated and
the validation is vacuous); [[lessons-L55]] (a census is only as wide as its
spelling); `vv-principles` #17.
