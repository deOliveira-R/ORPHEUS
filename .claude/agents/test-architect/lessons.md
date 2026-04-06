# Test Architect Lessons

Read this at the START of every invocation.

---

## 2026-04-05 — Cylindrical DD verification gaps

- **20 tests passed, bug survived**: homogeneous exact, particle
  balance, conservation, flux non-negativity — none caught the
  heterogeneous divergence.
- **The test that caught it**: heterogeneous 1G spatial convergence
  (keff must converge monotonically with mesh refinement).
- **Missing test type**: fixed-source flat-flux diagnostic — would
  have shown the flux spike at r=0 immediately.
- **Rule confirmed**: 1-group is degenerate. The 2G resolution
  independence test (4×8 vs 8×8) was the xfail that motivated
  the investigation.
- **Analytical reference gap**: no mesh-independent transport
  eigenvalue (Case's method) exists yet. Current references are
  either diffusion-based (~0.3% gap) or self-referencing (Richardson
  extrapolation). See IMPROVEMENTS.md DO-00000000-005.
