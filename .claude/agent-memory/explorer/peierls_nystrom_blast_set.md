---
name: peierls-nystrom-blast-set
description: Shape of the peierls_nystrom package's consumers (solve vs closed-form primitive split, registry reach, sole-verifier markers) measured 2026-09-24 for the Q-R2 quarantine
metadata:
  type: project
---

Measured 2026-09-24 at HEAD 901f64ca (#405 step 2, Q-R2: the user ruled the Peierls Nyström references not research-grade; they must neither be used nor run until improved).

- `peierls_nystrom/geometry.py` mixes two halves in one module: Nyström SOLVE machinery (build_volume_kernel, closure operators, solve_peierls_1g/mg) and closed-form primitives (compute_P_ss/T_specular/P_esc/G_bc, reflection_*, CurvilinearGeometry). A quarantine keyed on "imports peierls_nystrom" over-reaches: 12 files import only primitives.
- 0 production modules outside the package import it; the registry reaches it only by `pkgutil.walk_packages` + `continuous_case_builders`. Nexus impact shows flat_source_cp/cp.solver upstream: those are docstring `references` edges, not imports.
- Sole-verifier casualties: `flat-source` (a CP equation, only via cp/test_peierls_flux.py), 6 ERR entries (027-030, 032, 063).
- Per-test attribution needs a module-local call closure (helpers/fixtures) plus the subprocess-string import in cp/test_peierls_rank_n_protocol.py.

**Why:** the quarantine and later the "research-grade" rebuild will need this cone again.
**How to apply:** re-run the censuses (scratch scripts are session-local; rebuild from this description) before trusting counts; see [[sn-k-solve-walls-and-build-cadence]] for the other slow-tier shape.
