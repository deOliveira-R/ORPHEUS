---
name: adjointability-flag-semantics
description: What is_adjointable means at edf55489 (2026-10-18). One flag means "transpose coded" and gates both .dual() and .H; it is True on multi-D sweeps that raise, and no flag checks the kernel condition (folded_product(4,4) L=4 R and Λ fail it)
metadata:
  type: project
---

Measured 2026-10-18 at `edf55489`. Re-verify before repeating.

- 52 `is_adjointable` definitions. None reads a metric. The base default is `False`, inherited by 15 of 76 subclasses (8 concrete, all reason "no transpose coded"). `.dual()` and `.H` both gate on it. `.H` adds two gates the flag does not report: the UNBOUND refusal (`RankOneOperator` is True and its `.H` raises) and the Riesz legs' `DualSpace` `TypeError` (`riesz_lower.H` raises).
- Refused with a known matrix: `LeakageOperator` and `DiffusionBoundaryOperator`. `op.assemble()` (`SparseAssembledOperator`, True) daggers the same matrix.
- True while the transpose raises: `SweepOperator` on every 2-D pose. `S.H` builds and `.apply` raises `NotImplementedError` from `sweep_transpose` (1-D-scan-only, #280). The docstring's "multi-D defer via inner.is_adjointable" is stale.
- Degenerate metrics: dense heads (`G⁺`), dead slots, and the `|Ω·n|w` trace. The ψ½ space refuses a zero weight. The kernel condition `G_W A ker G_V = 0` FAILS for frame `R` and for `Λ` on `folded_product(4,4)` L=4 (`REC_all` O(1); `A.H.H is A` false, `HH = 0.38`). It holds on L, C, B and the sweep over the trace. No production `.H` daggers a frame face or `Λ` alone; composites dagger by their ends.
- Instruments: `scratch/numerics_untangling/adjoint_posed/` (`adjoint_census_plugin.py` wraps `AdjointOperator.__init__`; it reddens a forward-apply-counting gate, so run without `-x`).

See [[harmonic-frame-folded-quadrature-plumbing]], [[affine-operator-split-convention]], [[moment-metric-slots-and-cross-axis-forms]].
