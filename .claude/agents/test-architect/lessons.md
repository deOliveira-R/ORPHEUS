# Test Architect Lessons

Read this at the START of every invocation.

---

## L1: Curvilinear geometry needs heterogeneous convergence tests

Homogeneous tests (exact k_inf, particle balance, conservation,
non-negativity) are ALL blind to redistribution bugs in curvilinear
coordinates because flat flux makes geometry terms vanish. The
**minimum** test that catches these bugs: heterogeneous spatial
convergence (keff differences must shrink with mesh refinement).
Always include one for any curvilinear solver.

## L2: Follow your own procedure — especially fixed-source diagnostics

The fixed-source flat-flux test (Q/Σ_t) was already in AGENT.md
as "the single most powerful diagnostic for curvilinear bugs" but was
not included in the cylindrical DD test suite. It would have caught
the bug immediately (flux spike at r=0). When designing a test plan,
explicitly check every item in the analytical references list and
confirm each has a corresponding test.

## L3: No mesh-independent transport eigenvalue reference exists yet

Heterogeneous eigenvalue references are either diffusion-based (~0.3%
gap) or self-referencing (Richardson extrapolation). Case's method
would close this gap. Tracked: GitHub Issue #8.
