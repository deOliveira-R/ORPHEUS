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
source) cache is the canonical instance — see `[[feedback-elegance-causes-collapse]]`
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
