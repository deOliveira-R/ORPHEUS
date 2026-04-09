# Test Architect Lessons

Read at START of every invocation.

---

## L1: Homogeneous tests are blind to curvilinear redistribution bugs

Flat flux makes geometry terms (alpha terms, area ratios) vanish.
Homogeneous k_inf, particle balance, conservation, non-negativity —
all pass even with broken curvilinear discretization. **Minimum catch**:
heterogeneous spatial convergence (keff differences shrink with mesh
refinement). Mandatory for any curvilinear solver.

## L2: Walk the analytical-references checklist explicitly

The fixed-source flat-flux test (Q/Sigma_t) was already documented as
"the single most powerful diagnostic for curvilinear bugs" but was
omitted from the cylindrical DD test suite — and the bug hid. When
designing a test plan, enumerate every item in the analytical references
list from AGENT.md and confirm each maps to at least one test. Do not
rely on memory; print the checklist and check it off.

## L3: No mesh-independent transport eigenvalue reference exists

Heterogeneous eigenvalue references are diffusion-based (~0.3% gap) or
self-referencing (Richardson extrapolation). Case's method would close
the gap. Tracked: GitHub Issue #8. Use diffusion eigenvalue as a
cross-check with explicit tolerance, never as a precision target.
