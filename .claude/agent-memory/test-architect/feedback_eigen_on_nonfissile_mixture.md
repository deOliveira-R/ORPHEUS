---
name: eigen-on-nonfissile-mixture-is-malformed
description: A k-eigenvalue snapshot on a non-fissile mixture is malformed (k=production/absorption=0/abs → 0/0 nan); reformulate as fixed-source to exercise the same operator path
metadata:
  type: feedback
---

An eigenvalue (`solve_sn`) snapshot/test on a **non-multiplying** mixture
(`Σ_f = νΣ_f = 0`, e.g. xs_library mixture B "moderator-like") is
MALFORMED: `k = production/absorption = 0/abs` and the
production-rate-normalised power iteration divides by zero. The frozen
snapshot is all-NaN and gates nothing — a silent dead test.

**Why:** ORPHEUS SN regression redesign 2026-06-01. The two
`slab/sphere_2g_p1_aniso` regression cases existed to pin the Pℓ
Galerkin assembly path, but were built as `solve_sn` (eigen) on mixture
B. The legacy snapshots stored `keff=nan` + all-nan flux; regenerating
them raised `ZeroDivisionError` in `compute_keff`.

**How to apply:** when a verification/regression case wants to exercise
a SCATTERING/STREAMING operator path (Pℓ Galerkin, curvilinear closure,
angular redistribution) but the mixture has no fission, reformulate as a
**fixed-source** problem (`solve_sn_fixed_source` with a uniform
external source), NOT an eigenvalue problem. The same operator path runs
and the problem is well-posed.

Independent corroboration for the reformulated fixed-source case:
a uniform isotropic source + reflective BCs gives the flat
infinite-medium scalar flux `φ = (diag(Σ_t) − Σ_{s0}^T)^{-1} Q`
(closed-form linear algebra on the P0 scattering matrix — structurally
independent of the transport discretisation; the P1 moment vanishes for
flat flux since the current is zero). For mixture B 2G with Q=1/group it
is `[5, 38]`. High scattering ratio (c≈0.97) → SI needs a deep inner
budget (`max_inner≈3000`) to reach the convergence floor; use explicit
`inner_solver="source_iteration"` (the flat reflective solution is exact
for the WDD sweep, so SI is bit-stable and ~100× faster than the
curvilinear krylov auto-flip here).

General test-design rule it sharpens: before snapshotting/asserting a
k-eigenvalue, confirm the mixture is multiplying. A non-multiplying
medium has NO eigenvalue. (Related: 1-group eigenvalue degeneracy is a
DIFFERENT problem — there k exists but is shape-independent.)

See [[regression-tolerance-design]] for the gate that runs these cases.
