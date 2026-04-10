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
