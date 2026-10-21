---
name: s6bc-census-facts
description: S6b/S6c census facts — PG frame is a superclass of Galerkin; SN leaves annihilate tangential trace slots; a dual head fails MomentHead
metadata:
  type: project
---

Written by the orchestrator at the explorer's request (2026-10-21); the explorer's write-scope hook refused its own path. Content is the explorer's, uncurated.

- `GalerkinFrame` and `HarmonicFrame` SUBCLASS `PetrovGalerkinFrame`, so `isinstance(f, PetrovGalerkinFrame)` admits every SN frame. Use the exact type to count the PG case (a first table of 1 904 reads was wrong for this reason; exact-type count 0).
- Every production SN leaf with a degenerate end maps ker G (tangential trace slots) to exactly 0 (87 of 87 rows; control 0.70). Today's `.H` is the quotient adjoint; the `A.H.H` gap lives on the ker→ker block only.
- `isinstance(h.dual(), MomentHead)` is False for both head families (the dual lacks `L`, `isotropic_slot`, `degree_block`, `truncated`), while `h.dual() == h` is True today.
- Report: `scratch/numerics_untangling/packages/S6bc_census/report.md`.
