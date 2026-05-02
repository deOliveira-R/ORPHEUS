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
