---
name: characteristic-reference-sn-consumers
description: The 13 SN rows reading the old trajectory-resolvent reference (P1 step d), their single re-point seam, and the traps (label collisions, pinned names, scattering order)
metadata:
  type: project
---

Stamped at HEAD `39ebc4b1` (2026-10-09). Delete this once the characteristic-reference campaign's step (e) merges.

**Shape:**
- 13 SN test functions (17 ids) read the old reference. An AST census over 1132 tracked files, with 7 positive controls, finds no 14th (`scratch/characteristic_architecture/p1_step_d/explorer/census.py`).
- 11 of the 13 go through ONE seam, `_aba_reference.aba_reference` / `aba_reference_at`.
- The 2 partial-reflector rows and the 3 Gate-4.2 edge ids call the array solvers directly.

**Why:** step (d) re-points them; (e) retires the family.

**How to apply:**
- The spec's "D9" (self-convergence) and the rung-5 spec's "[D9]" (the Resolution refusal, `test_characteristic_reference.py`) are different rows. A grep of "D9" finds the wrong one.
- `test_crosscheck_harness.py` pins the xfail rows by NAME and by the substring "#516" in the reason, so renaming or rewording reds a foundation row.
- "SN Garcia rows" do not exist. The Garcia consumers are the old family's own rows.
- The SN partial-reflector rows run at scattering order 0 (verified at runtime, with a P1 control), so `isotropic_mixture("A")` poses the same problem.
- The new reference has no error-estimate callable (only `Uncertified`). The cylinder has one-axis-down steps only (degree 3.4e-8, transport 6.9e-11).

Related: [[reference-producer-landscape]], [[trajectory-resolvent-family-inventory]].
