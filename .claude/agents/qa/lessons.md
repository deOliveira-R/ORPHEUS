# QA Agent Lessons Learned

Read at START of every invocation. Keep sharp: merge overlaps, cut
filler. AGENT.md has the full error catalog, anti-patterns, and
verification hierarchy — never duplicate that here. Only record
what changes future behavior beyond what AGENT.md already encodes.

---

## L-001 — Test count is not coverage (2026-04-05, cylindrical DD)

20 passing tests (homogeneous exact, conservation, balance,
non-negativity) missed a fundamental 2-term bug. Failure signature:
keff diverging under mesh refinement (1.15 → 0.90 → 0.52).

**Behavioral rule**: When reviewing any "all tests pass" claim for a
new geometry or solver, first ask: "Is there a heterogeneous
mesh-refinement convergence test?" If not, the suite proves nothing
about the transport operator's correctness.
