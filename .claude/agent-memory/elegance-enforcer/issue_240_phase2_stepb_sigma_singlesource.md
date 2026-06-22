---
name: issue-240-phase2-stepb-sigma-singlesource
description: #240 Phase 2 Step B — InvertibleOperator owns its matvec; σ single-sourced from C; loss_action(operator)→loss_action(sigma) across 9 classes. PASS-with-nits.
metadata:
  type: project
---

#240 Phase 2 Step B (`feature/sn-space-angle-tier2`, uncommitted at review). **PASS-WITH-NITS** (no BLOCKER).

**The carve** closes the two-sources-for-σ leak: `InvertibleOperator=L+C` previously inherited
`OperatorSum.apply = L.apply+C.apply` where `L.apply` read σ from `L.sigma_t` while `.solve` read it
from `C` — value-correct only by the affine-in-σ coincidence (`M(σ)ψ = streaming_action(ψ)+σ·ψ`).
Now: (1) protocol `loss_action(self, operator, psi)` → `loss_action(self, sigma, psi)` across all
**9 classes** (Protocol, `_LossRepresentation`, `_OctantWalk`, `CumprodScan`, `_DAGWavefront`,
`MovingFrontierWindow`, `FullFieldWavefront`, `ScanMarch`, `_OneDimScanWalk`) + `_transpose`;
(2) `InvertibleOperator.apply/apply_transpose` overrides = `loss_action(self.sigma, psi)` (σ from C);
(3) `_require_typed_composite` helper single-sources the typed+mesh-identity guard for both operators.

**RULINGS (verified, durable):**
- **Signature uniformity = COMPLETE.** All 14 `def loss_action*` typed `sigma: "np.ndarray"`; zero
  production caller passes `operator`; zero representation method reads `.sigma_t` off an op handle
  (grep-clean). Internal helpers (`loss_action_decomposed`/`_apply_walk`/`sweep`) keep their plain-array
  `sigma_t` param — that's the ESTABLISHED #206 convention, NOT a twin.
- **The override is NOT "sometimes distributes, sometimes doesn't" confusion — it REMOVES asymmetry.**
  `OperatorSum` (numerics/operator.py:660) docstring: apply/apply_transpose distribute, **solve does
  NOT propagate** (no generic `(A+B)⁻¹`). `InvertibleOperator` ALREADY overrode `solve` (WDD sweep).
  Pre-carve state was the genuine asymmetry (solve fused, apply distributed). Carve makes apply/transpose
  ALSO fused → the composite owns its whole apply/solve/transpose triad coherently. Correct direction.
- **Resolution-A two-conventions web is COHERENT.** `StreamingOperator.apply` (L leaf) calls
  `loss_action(self.sigma_t)` then SUBTRACTS `σ_t·ψ` → bare `Lψ`. `InvertibleOperator.apply` (composite)
  calls `loss_action(self.sigma)` and does NOT subtract → full `M(σ)ψ`. Both consume the SAME `loss_action`
  returning `(L+C)ψ`. Documented accurately at both sites + theory `:label:loss-rep-removal-form-matvec`.
- **`_require_typed_composite` = good Pattern-2 extraction** (was duplicated in StreamingOperator.apply +
  apply_transpose; now 4 consumers). Naming/placement fine.
- **Tests EXEMPLARY** (`test_removal_form_matvec_sweep.py`): constructs σ_r≠σ_t (the removal form
  production CANNOT yet build) → proves override correct-by-construction in the regime the carve unlocks;
  pins apply+apply_transpose bit-id vs `M(σ_r)`, matvec≡sweep round-trip, + structurally-indep k_inf.
- **Re-baseline HONEST** (5 .npy: cyl+2D matvec ≤5 ULP, FP-non-assoc — override drops redundant
  `(x−σψ)+σψ` round-trip). Theory todo-block is a TRACKED artifact (anti-pattern-11 exception).

**NITS (do-now, approval conditions):**
- NIT-1 (Pattern-3/7 latent): `loss_action(sigma)` forwards `sigma` into `_apply_walk(sigma_t=...)` →
  the diamond cell-balance kernel param literally named `sigma_t` RECEIVES σ_r in the removal form.
  Latent only (prod = σ_t today; solver.py:217/894/969 all `total_cross_section`); but this carve is
  WHAT makes σ≠σ_t cleanly reachable, so the lie survives one layer deeper. Rename internal helper
  param `sigma_t`→`sigma` (or `sigma_diag`) along the loss_action→_apply_walk thread to finish the job.
- NIT-2 (anti-pattern-11, stale doc on live code): `InvertibleOperator` CLASS docstring still says
  "**inherits** the OperatorSum apply (the sum of the operand actions)" (operator.py:~1167-1172) — the
  carve OVERRODE it; and "apply_transpose is **reserved for Phase H**" (~1200-1201) — the carve ADDED it.
  Both refuted by the same commit. Fix both.

Discriminator carry-fwd: when a carve makes σ≠σ_t reachable, grep every `sigma_t`-named param on the
forward thread — the name is a convention assertion the value can now violate.
