# QA Agent Lessons Learned

This file is the QA agent's memory. Each entry records bugs found,
bugs missed, verification gaps, and anti-patterns discovered.

**Read this file at the START of every invocation.**

---

## 2026-04-05 — Initial lesson from cylindrical DD investigation

- **Bug missed by 20 tests**: wrong α recursion + missing ΔA/w factor
  in cylindrical sweep. Passed: homogeneous exact (1G/2G/4G), particle
  balance, conservation, flux non-negativity, single sweep finite.
- **Verification level gap**: no L1 heterogeneous convergence test with
  mesh refinement. The divergent keff (1.15→0.90→0.52) was the ONLY
  observable that caught the bug.
- **Anti-pattern**: "homogeneous exact implies correct" — false for
  curvilinear geometries where redistribution errors cancel for flat flux.
- **Test to add (added)**: `test_heterogeneous_1g_spatial_convergence` —
  keff differences must decrease with mesh refinement.
- **Key insight**: per-ordinate flat-flux consistency is the fundamental
  correctness criterion for curvilinear SN. Test it with a fixed-source
  diagnostic (uniform Q, check volume-avg and flux range).
