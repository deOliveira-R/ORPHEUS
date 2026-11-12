---
name: v-star-density-chart-sites
description: Where ORPHEUS condenses V* (sources/currents/responses) into densities, as of refactor/true-ends 5d40c6d9 (2026-11-11, #707 R6) — two V* charts coexist, two sweep conventions, one /W home, trace metric lacks area.
metadata:
  type: project
---

Shape of the data side V* at `5d40c6d9` (orpheus identical to `a00539f4`). Full tables:
`scratch/true_ends/census/r6_density_sites.md`; runtime spy `scratch/true_ends/census/r6/plug/src_spy.py`
(wraps the 7 SourceSink leaves' dataclass `__init__` + a `values` data descriptor for test reads).

- **Two V* charts already coexist.** `FunctionSpace.flat` → `Covector` on `V.dual()` stores G·w (integrated);
  the 7 `*SourceSink` leaves store densities. R6 = SourceSink joins Covector.
- **The sweep has two conventions.** 1-D scan walk multiplies the density by V and solves the integrated cell
  balance (`QV_per_ord = Q_per_ord * V_b`); the multi-D DAG kernel is per-volume (`numer = Q_cells`, s_axes = |μ|/Δa).
  Neither ever sees w_n: the producer divides by W.
- **/W has one home**: `AngularLift._combine` (`iso / self.total_weight`), ~3 M mints per not-slow run; plus
  `TransferOperator._redistribute_*` and `AngularSourceSink.from_isotropic`.
- **Trace metric = |Ω·n|·w_n with NO face area** (`_build_trace_metric_weights`), while the SN left-kernel
  `_trace_flux_weights` spells (Ω·n)·A_f. The boundary inflow and the reflection B emit trace "sources" in ψ units.
- **Gauges**: forward SN k and gauge read flux-side `IntegratedReactionRate` (chart-proof); the adjoint
  `KEigenvalue` gauge/rayleigh(w=1) and diffusion `compute_keff` (`_integration.evaluate(loss.apply(...))`) read
  V* values and move/raise under an integrated chart.

**Why:** design-v2 of #707 sequences the carve on these seams. **How to apply:** re-verify the anchors by grep
before citing (the carve will move them); the spy plugin is reusable for any "who mints class X" question.
Census lesson: the word "reference" in a test body is NOT a frozen snapshot — classify FROZEN by
`np.load|.npz|snapshot|digest` only (147 → 51 on this census). See [[census-predicates-bound-reference-and-activation-traceback]].
