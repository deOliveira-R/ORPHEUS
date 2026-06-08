---
name: si-convergence-rate-verification
description: Durable lesson — how to VERIFY a solver convergence-RATE (iterations-to-converge) claim against a structurally-independent target (analytic SI spectral radius ρ_J=c). LANDED test_si_convergence_rate.py. The measurand (history.n_inner), the still-OPEN eigenvalue-n_inner-is-None API gap, the rate-not-eigenvalue 1G-OK exception.
metadata:
  type: feedback
---

Verifying a convergence-RATE claim (iterations-to-converge, an acceleration
gain) is a distinct test-design problem from verifying a VALUE. The SN SI
Gauss-Seidel rate-recovery gate is LANDED:
`tests/sn/verification/analytical/test_si_convergence_rate.py` (6 passed,
4 xfailed). This note keeps the durable method.

**1. The measurand is `history.n_inner`, and there is a STILL-OPEN gap.**
A rate claim needs the iteration COUNT, not the converged value.
`solve_sn_fixed_source(...).history.n_inner` carries it. **OPEN VERIFICATION
GAP (preserve this breadcrumb):** the EIGENVALUE path `solve_sn(...)` leaves
`n_inner = None` — `KEigenvalue`/`power_iteration` drops `_inner_residuals` at
`iteration.py`. A rate gate on the eigenvalue inner is therefore UNWRITABLE
until a recovery PR adds an `IterationHistory.total_inner_iterations` seam
(`IterationHistory` lives in `orpheus/sn/solution.py`). The G-S→Jacobi SI rate
regression in Wave-O O.4a was INVISIBLE precisely because nothing measured
n_inner — the only assertion was `test_solution.py` checking `n_inner is None`.
A rate claim that nothing measures is a silent un-gated claim.

**2. The structurally-independent TARGET is the analytic SI spectral radius
ρ_J = c (scattering ratio).** For a homogeneous reflective medium the
source-iteration error contracts at ρ = c (Adams-Larsen / Hageman-Young
consistently-ordered theory). The iteration count then scales as
log(tol)/log(c). This is the closed-form anchor that makes a rate measurement
correctness evidence rather than self-consistency. Measured: mixture-B 2g SI
= 655 iters vs log(tol)/log(c)=728 (ratio 0.90); the model-problem G-S radius
ρ_GS ~ c² → ~halves the count.

**3. Gate the rate CONSERVATIVELY and keep the value-correctness guards
separate.** The rate gate is `n_inner < 0.75 × Jacobi_baseline` (xfail-strict
RED today → flips green on the recovery). 0.75 (not the theoretical 0.5) is
self-sufficient and robust to FP-reassociation drift in the baseline — DEFER
tightening to 0.55 (the literature-researcher Adams-Larsen pull) until needed.
PAIR the rate gate with independent VALUE guards: G-1 k_inf=1.875 rtol 1e-10,
G-2 SI≡Krylov rel-Linf<1e-8, G-3 foundation flat-balance limit, G-4 explicit
vacuum=128±2 (with `BC.vacuum` set on the MESH, not a kwarg — same G-4 dud as
[[phase4-46-nonvacuum-mms-ansatz]]).

**4. The 1G-OK EXCEPTION (rate claims, unlike eigenvalue claims).** A rate
claim ρ=c is FLUX-SHAPE-INDEPENDENT BY DESIGN — it is a property of the
iteration map, not the solution. So a 1-group homogeneous slab is a LEGITIMATE
rate-claim test (contrast: a 1G EIGENVALUE test is degenerate per Cardinal-6).
State the claim layer explicitly: rate-not-eigenvalue → 1G is fine here.
Regime choice: mixture B (c→1, non-fissile) is the right stress — mixture D was
REJECTED (Σ_t≈0.04 thin breaks the log-model, ratio 3.87).

**5. Robustness regime + the `Solution.compare` SAME-MESH gotcha.** Use a
non-fissile high-c mixture (mixture B, c→1) for the c→1 stress. `Solution.compare`
is SAME-MESH-ONLY (`self.mesh is other.mesh`, raises "meshes differ" on any two
separate solver calls) → compare `scalar_flux.values` arrays DIRECTLY for the
SI≡Krylov leg, not via `.compare`. NO module `pytestmark` (it clobbers the
foundation tag on G-3 with a conflict warning) — tag each test explicitly.

vv Mode-9 (verify splittings to the same FP on anisotropic/diagonal-stressing
configs) and ERR-056 (shared-face premature reflect) are the failure modes this
work surfaced; both are now canon in the vv-principles skill — NOT re-logged
here. See [[regression-tolerance-design]] (SAFETY×conv_tol),
[[eigen-on-nonfissile-mixture-is-malformed]].
