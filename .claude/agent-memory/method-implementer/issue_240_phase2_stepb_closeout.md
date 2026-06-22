---
name: issue-240-phase2-stepb-closeout
description: #240 Phase 2 Step B — InvertibleOperator owns its matvec (loss_action(self.sigma)); loss_action protocol signature (operator,psi)→(sigma,psi); removal-form C(σ_r) correct-by-construction. NOT pushed; NOT committed (user commits after review).
metadata:
  type: project
---

# #240 Phase 2 Step B — composite owns its matvec (correct-by-construction)

Branch `feature/sn-space-angle-tier2`. 2026-06-15. NOT committed (user reviews
+ commits). An ARCHITECTURE carve, not a bug fix: the matvec value never shipped
wrong (the inherited leaf sum was value-correct by the affine-in-σ coincidence);
the carve makes `(L+C)` own ONE `loss_action(σ_C)` single-sourcing σ from C —
removing the latent affine-in-σ trap and the `operator`-as-σ-carrier smell.

**Why:** `.solve` sources σ from C (`self.sigma`), but the inherited
`OperatorSum.apply = L.apply + C.apply` sourced σ from L's OWN `sigma_t`. Two
sources that agreed only because production builds L and C from the same σ_t. The
forward WDD matvec is AFFINE in σ (`M(σ)ψ = streaming_action(ψ) + σ·ψ`), so the
leaf sum `L.apply(σ_t) + C.apply(σ_r) = streaming_action + σ_r·ψ = M(σ_r)ψ` is
the right value but the wrong source. The override `loss_action(self.sigma)`
single-sources σ from the diagonal.

**How to apply:** when refactoring an `OperatorSum` subclass whose `.solve`
sources a coefficient from one operand, audit whether the inherited `.apply`
sources it from the OTHER operand — a latent twin-source. If the action is affine
in that coefficient, the leaf sum is value-correct-by-coincidence; override
`.apply` to single-source from the same operand `.solve` uses.

## The two coordinated changes (production)

1. **Protocol signature `(operator, psi)` → `(sigma, psi)`** on `loss_action` /
   `loss_action_transpose` across ALL of `orpheus/sn/loss_representation.py`:
   `LossRepresentation` Protocol, `_LossRepresentation` base, `_OctantWalk`
   (body `sig_t = operator.sigma_t` → `sig_t = sigma`; keeps its 3rd `interior`
   param), `CumprodScan`, `_DAGWavefront` (transpose deferral), `MovingFrontierWindow`,
   `FullFieldWavefront`, `ScanMarch`, `_OneDimScanWalk` (forward → `_apply_walk(sigma,...)`;
   transpose body `sgx = operator.sigma_t` → `sgx = sigma`). `loss_action_decomposed`
   / `_apply_walk` were ALREADY taking a plain `sigma_t` — UNTOUCHED. Docstring at
   `_apply_walk` (was ~1943) updated. Caller sites in `orpheus/sn/operator.py`:
   `StreamingOperator.apply` (`loss_action(self.sigma_t, psi)`) +
   `apply_transpose` (`loss_action_transpose(self.sigma_t, phi)`), Resolution-A
   `−σ_t·ψ` subtraction KEPT.

2. **NEW `InvertibleOperator.apply` / `apply_transpose` overrides** (near `solve`,
   `orpheus/sn/operator.py`): `return self.loss_representation.loss_action(self.sigma, psi)`
   / `...loss_action_transpose(self.sigma, phi)`. Input guards factored into a NEW
   module-level helper `_require_typed_composite(method, sn_mesh, field)` (Pattern 2
   single-source — `StreamingOperator.apply` + `apply_transpose` now consume it too,
   replacing the duplicated inline TypeError + mesh-identity guards). `apply_transpose`
   propagates into `InvertibleOperator.capabilities` automatically (OperatorSum adds
   CAP_APPLY_TRANSPOSE iff both L and C have it — they do). Multi-D Cartesian transpose
   raises the representation's deferral (never silent wrong).

## Verification result (HOST .venv/bin/python -O)

- Teeth gate `tests/sn/operators/test_removal_form_matvec_sweep.py`: **19 passed, 1
  xfailed** (the 7 strict-xfail-until-override teeth — 4× `apply_is_M_of_C_sigma` +
  3× `apply_transpose_is_M_transpose` — FLIPPED to plain pass; `_XFAIL_UNTIL_OVERRIDE`
  marker REMOVED; the 1 xfail is the `#200` k_inf solver-entry stub that stays xfail).
- Comprehensive: `tests/sn/operators spatial sweep/core sweep/cartesian_2d solve` with
  route-arounds → **1080 passed, 6 skipped, 7 deselected, 5 xfailed, 0 failed** (147 s).
- Strict DriftWarning gate `tests/sn/sweep/core tests/sn/solve -W error::...DriftWarning`
  → **505 passed, 1 skipped, 4 xfailed** (after the 3 re-baselines below).
- `tests/sn/regression` → 13 passed (pre-existing 6921-ULP SI 2-D DriftWarning within
  tol, NOT mine — SI rides solve not apply).
- G-adjoint reciprocity `test_g_adjoint_reciprocity.py` → 25 passed, 4 xfailed (`.H`
  wiring intact: `A=L+C−B` routes `(L+C).apply`/`.apply_transpose` through the override;
  value-equal at production σ).

## Re-baselines (3 — all principled, ≤5 ULP / rel ≤4e-17, vv-principles 3-criteria)

The override drops the `(x − σψ) + σψ` round-trip the leaf sum carried → the MORE
accurate path, so APPLY snapshots captured pre-override re-associate. SWEEP/SOLVE
snapshots UNTOUCHED (they ride `solve`, not `apply` — slab/sphere stay bit-identical).

1. `tests/sn/_data/affine_carve_baseline/matvec_bulk_{SPH,CYL}.npy` re-captured
   (`(L+C).apply` matvec leg: SPH 1 ULP, CYL 4 ULP). SLB bit-identical (not re-captured).
2. `tests/sn/_data/bc_extraction_baseline/vacuum_bulk_CYL_seed{0,1,2}.npy` re-captured
   (vacuum `(L+C).apply` bulk: seed0 46 ULP / rel 3.93e-17 on a near-zero output element,
   seed1 2, seed2 6). SLB bit-identical, SPH stays the pre-existing structural red
   #195/#209 (~1e15 ULP) — NOT re-captured.
3. `tests/sn/solve/test_affine_carve_bit_identity.py` GOLDEN: the TWO
   `krylov_2d_p1_aniso_het` sha256 hashes regenerated (2-D Krylov: the ≤5-ULP 2-D apply
   accumulated through GMRES inner_tol=1e-12 → converged bytes shift at ~1e-12). Verified
   structurally-independent: Krylov-2D converged φ ≡ SI-2D φ (which is bit-identical to
   ITS golden, doesn't use the override) to 3.9e-12 rel. SI 2-D + slab hashes UNTOUCHED.

## Re-baseline of the gate's OWN tolerance (test-architect spec correction)

`test_production_sigma_apply_value_preserved` was specced `bitexact=True` for
slab/sphere from a probe that compared two override-style `loss_action` paths (0/32).
The ACTUAL gate compares override-vs-leaf-sum: re-probed slab 2 ULP, sphere 2, cyl 4,
2-D 5. Dropped the `bitexact` flag; uniform `nulp=6` (headroom above 5) + the rel<1e-14
FP-re-association bound. The bit-identity premise is NOT a property of the override-vs-
leaf-sum pair (the override removes the round-trip). Documented inline + in the section
(c) comment.

## Test migrations (signature `(L,…)`/`(operator,…)` → `(σ,…)`)

- `test_removal_form_matvec_sweep.py`:283/323/441 (the teeth + value-ground references
  `loss_action(L_ref, psi)` → `loss_action(sig_r, psi)`).
- `test_loss_action_convention.py`:110/151 (`rep.loss_action(L, psi)` → `(sig_t, psi)`).
- `test_one_octant_walk.py`:149 (`rep_cls(sn).loss_action(L, psi)` → `(sig_t, psi)`;
  removed now-unused local `L`; `StreamingOperator` import stays, used by the other test).
- `test_scan_march_equivalence.py`:205-206 + `test_2d_full_field_oracle.py`:136-137
  (`.loss_action(L, state)` → `(sig_t, state)`; removed now-unused `L` + `StreamingOperator`
  import in BOTH).
- `test_scan_march_end_to_end.py`:94 — `MovingFrontierWindow.loss_action` is a spy-wrapped
  method REFERENCE (`*args,**kwargs`) — signature-agnostic, NOT migrated.

## Docstring re-anchors (item 3)

- `test_2d_l2_matvec_correctness.py` `test_apply_vs_sweep_2d_residual_cancellation`:
  KEPT (allclose at 1e-12 survives ≤2 ULP), docstring noted it is now an
  allclose-consistency check; the ERR-026 structural anchor MOVED to the teeth gate.
- `test_invertible_operator.py` `test_apply_equals_l_plus_c_on_typed_flux` (slab,
  array_equal): docstring noted override-owns-matvec but stays bit-identical on slab.

## Sphinx stub

`docs/theory/loss_representations.rst` NEW `.. _loss-rep-removal-form-matvec:`
subsection with a `.. todo::` archivist-expansion marker (NOT rich narrative). No
new `:label:` — the carve sharpens the EXISTING `loss-rep-resolution-a` equation; the
teeth gate already `verifies("loss-rep-resolution-a")`. DISPATCH archivist for the
narrative expansion.

## Lessons / skill-evidence

- **affine-in-σ leaf sum** = a value-correct-by-coincidence twin source. New
  `coding-elegance` Pattern-2 instance: an `OperatorSum` subclass whose `.solve` and
  inherited `.apply` source the same coefficient from DIFFERENT operands. Candidate
  anti-pattern #19 addition if it recurs.
- A test-architect's "bit-identical" probe can measure the WRONG pair (here: two
  override-style paths, 0/32) and miss that the production gate compares override-vs-
  leaf-sum (1-2 ULP even on slab/sphere). When a `bitexact=True` gate fires post-impl,
  re-probe the ACTUAL pair before assuming a regression.
- APPLY-path snapshots re-baseline under an apply-direction carve; SWEEP/SOLVE snapshots
  do NOT (they ride `solve`). The blast-radius split (SI bit-identical, Krylov drifts)
  is the structural-independence cross-check that proves the apply re-baseline principled.
