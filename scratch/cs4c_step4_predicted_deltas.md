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

## Battery B-4 — RUN 2026-08-31, in-process plugin arms (crash-safe by construction: subprocess-only mutation, no disk writes)

Scope: test_kernels.py + test_isotropic_fission.py + test_fission_adjoint.py + test_fission_kernel_crosscheck.py (91 rows).

| arm | mutation | measured verdict |
|---|---|---|
| B4.1 positive control | `FissionKernel.dyad` → zeros | **2 reds** — G1.7 (dyad-direction theorem) + G-F1 law row. ⚠ NOT "broad" as the plan predicted: `dyad()` is the oracle VIEW; the production route consumes the gathered FACTORS directly (kernel.chi/nu_sig_f), by design (`kernels.py` dyad docstring). The harness-alive role is served; the broad control is B4.6. |
| B4.2/3/4 wrong-morphism | in-gate CTRL rows | [M] separations restated in the gate (6.421e-1 / 1.685e0 / 7.087e-2); red-capable at the landing commit with no production change (§6c). |
| B4.5 degenerate target | in-gate | 2 precondition-red rows (1-fine-per-coarse + flat-φ) green as refusal tests. |
| B4.6 factor swap | `gather_chi` ↔ `gather_nu_sig_f` | **11 reds** across the four modules (G-F2, hand-rolled forward + dual-dyad transpose, the flip control, adjoint rows) — the broad catcher, multiple independent routes. |

## Test-count arithmetic additions (4d)
| change | tree | delta |
|---|---|---|
| G-F1 (law + 3 CTRL params + 2 B4.5 rows) | transport | +6 |
| G-F2 | transport | +1 |

## FINAL per-tree predictions (pre-registered 2026-08-31, before the exit gate; [M] baselines from scratch/_cs4c3_fast_gate.log, module deltas measured by --collect-only at HEAD vs 4c613ce7)

| tree | baseline | predicted | delta source |
|---|---|---|---|
| numerics | 2538 | 2538 | test_iteration re-keys are body-swaps, count-neutral |
| transport | 630 / 1sk | **645 / 1sk** | +7 material_field ([M] 28→35) +7 kernels ([M] 53→60: G-F1 ×4 + B4.5 ×2 + G-F2) +1 tier-2 G-C1 angular-F row (added post-registration, pre-gate — [M] 7→8 in test_tier2_equivalence_s_family.py) |
| geometry | 727 / 4sk / 1xf | unchanged | — |
| data | 227 | unchanged | — |
| homogeneous | 50 | unchanged | re-keys count-neutral |
| diffusion | 113 | unchanged | re-key count-neutral |
| cp | 141 | unchanged | — |
| moc | 121 | unchanged | — |
| mc | 39 / 2xf | unchanged | — |
| cross_method | 81 | unchanged | — |
| sn | 3361 / 1sk / 47xf | **3373** | +12 test_isotropic_fission (NEW); n2n [M] 10→10; every other re-key body-swap |
| derivations | 1637 / 13sk / 11xf | unchanged | — |
| root+harness | 414 / 5xf | 414 **± label delta** | the P7 lesson: the documented-label registry gate gains a param per new equation label — RESOLVE from the archivist's return BEFORE the readout, not after |

**Total: 10106 / 0 / 19sk / 66xf (+ label-gate delta), 13 trees rc=0.**
