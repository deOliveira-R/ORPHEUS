# CS4c step 4 — pre-registered per-tree gate predictions (DRAFT; finalize before the exit gate)

Baseline (§15): **10079 passed / 0 failed / 19 skipped / 66 xfailed**, 13 trees rc=0
(driver `scratch/_cs4c3_gate_driver.sh`, log `scratch/_cs4c3_fast_gate.log`).

## Test-count arithmetic (running tally; [M] where run)

| change | tree | delta |
|---|---|---|
| test_isotropic_fission.py NEW (4a) | sn/operators | **+12** |
| test_material_field.py fission tiers NEW (4a) | transport | **+7** [M] (35 − 28 pre) |
| test_kernels arm matrix F-ScalarFlux cell re-keyed (same count) | transport | 0 |
| test_fission_adjoint/crosscheck/operator/typed/capability/psi_half/residual/codomain re-keys (same row counts) | sn/operators | 0 |
| capability survey +1 row member (F energy + F angular in one list — list rows not parametrized) | sn/operators | 0 |
| homogeneous/diffusion re-keys (same counts) | homogeneous, diffusion | 0 |
| 4c N2N harmonization gate re-keys | sn/operators | TBD |
| 4d G-F1 (law + 3 controls + activation precondition) | transport or data | +5 (plan §7) |
| 4d G-F2 | transport | +1 |
| 4d G-C1 rows (IsotropicFission classmethod≡ctor landed IN test_isotropic_fission at 4a) | — | counted in +12 |

## Behavioural notes for the readout
- F composite transpose: divide-order ULP change (product route divides by W AFTER
  the dyad; old divided before). [M] max rel 1.3e-16 on the 4b probe. Every pinning
  gate is tolerance-based (test_fission_adjoint rtol 1e-12 — [M] green).
- Forward arms bit-identical ([M] 4a probe: gathers + scalar arms array_equal;
  composite forward ≤ 1 ULP vs hand tensordot, same integrate_angular route).
- keff paths: compute_fission_source arithmetic unchanged (same dyad, same buffer
  values) ⟹ snapshot suites must NOT move. Any keff snapshot delta = a defect.
- test_kernels row-3 conformity: reach WIDENED (wrong-ng on axes-less SN composite
  now refuses through the derived interior) — re-keyed, same row count.
- Mode-11 sentinel (crosscheck): caught the k-outer re-route as designed; re-keyed
  onto IsotropicFission.apply in the same commit.
- tests/sn/solve/test_sn_adjoint_certification.py patches FissionOperator.apply_transpose
  (3 sites) — the adjoint solve's F routing changed homes; EXPECT these to need
  re-keying to the live surface (check in the sn tree run).
