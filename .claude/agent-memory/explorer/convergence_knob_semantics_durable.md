---
name: convergence-knob-semantics-durable
description: The five iterating families' convergence knobs — which are the SAME quantity (dk, max_outer), which differ only cosmetically (dphi), and which are NOT the same quantity at all (inner_tol = residual in SN vs increment in CP/MoC, differing by 1/(1-rho)).
metadata:
  type: project
---

Measured 2026-08-13 at `0f5ca91c` for issue **#364** (single-sourcing the
convergence-knob defaults). Line numbers drift; the SEMANTICS below do not.

**The fact.** Five families iterate (`sn`, `cp`, `moc`, `diffusion`,
`numerics.KEigenvalue`); `homogeneous` does not (direct `dominant_eigenpair`).
They share the ONE `power_iteration` outer loop, and each realizes
`measure_stopping_criteria` itself — deliberately, so the loop never has to pick
a metric.

**Why:** #364 asks whether the per-family tolerance differences are decisions or
drift. The answer turns entirely on whether the knobs denote the same quantity.
Three of the five knobs do; one emphatically does not, and that one is where a
naive single-source would silently change every family's meaning.

**How to apply:** before ANY proposal to share a convergence literal across
families, classify the knob:

- **`dk` / `max_outer` — SAME quantity.** All five: absolute `|k - k_old|`, no
  normalisation, `k ~ 1`; and one outer = one power iteration. A shared literal
  is dimensionally defensible.
- **`dphi` — same SHAPE, four cosmetic differences, all small.** CP judges in
  relative `l-inf`, the other four in relative `l2` (stated at
  `numerics/eigenvalue.py` `measure_stopping_criteria` Protocol docstring and
  `cp/solver.py`'s own). `[M]` the ratio `l-inf/l2` measures **1.13–3.4**, i.e.
  far under the 10x tolerance gap the families ship — **the norm choice does NOT
  justify a different literal.** Also: only SN / diffusion / KEigenvalue
  implement `compute_production_rate`, so only they get `power_iteration`'s
  unit-production renormalisation; CP opts out BY DESIGN (it conditions with
  `phi *= 1/max(phi)`), MoC never adopted it — so CP/MoC's `dphi` carries scale
  drift where SN/diffusion's is pure shape. And diffusion's carrier is the FLAT
  COMPOSITE (bulk ⊕ boundary trace), `[M]` trace = 0.46 % of the norm on a
  50-cell 2G vacuum slab — small, but it grows with the face/cell ratio.
- ⛔ **`inner_tol` — NOT THE SAME QUANTITY. This is the load-bearing finding.**
  - SN source-iteration: the **ρ-honest EQUATION RESIDUAL** `||A psi - q|| /
    ||q_ext||` — explicitly "a DELIBERATE re-interpretation of `tol`"
    (`SourceIteration`'s class docstring derives it).
  - SN Krylov: scipy GMRES `rtol` on the **preconditioned** residual vs `||b||`.
  - CP Gauss-Seidel inner: a relative `l2` **INCREMENT**, per group.
  - MoC inner: a relative Frobenius **INCREMENT**, on the scalar flux AND the
    boundary angular flux.
  - diffusion: no inner tolerance (exact LU resolvent).

  `[M]` residual/increment `= 1/(1-rho)` **exactly** (the docstring predicts it
  from Banach). Measured on a 40-state contraction: `c=0.5` → the two agree
  (1 extra pass); `c=0.9` → increment is **10x looser** (+22 passes); `c=0.99` →
  **100x** (+458); `c=0.999` → **1000x** (+6904). ⟹ the literal `1e-8` in SN and
  the literal `1e-8` in CP/MoC are true-error claims differing by an UNBOUNDED
  factor. **A single source must name the METRIC, not share the number.**
- **Inner BUDGET — different units per family.** SN/KEigenvalue derive it from
  the tolerance (`resolve_iteration_budget` / `default_iteration_budget` at
  `_SERVED_RATE = 0.986`); CP ships a literal `100` **per group and
  Gauss-Seidel-only** — and `solver_mode` defaults to `"jacobi"`, whose
  `_solve_fixed_source_jacobi` is a SINGLE PASS with no loop, so **CP's inner
  knobs are dead on the production default path**; MoC ships a literal `200`
  transport sweeps; diffusion has no inner level. Krylov additionally counts
  restart CYCLES where the record counts Arnoldi steps (the #349 /
  `IterationBudget` exchange rate).

**The `flux_tol = 10 x keff_tol` regularity is real but unwritten.** `[M]` at the
converged tail `dphi/dk` measures **9.2–12** (CP), **3.4** (SN 1-D slab het 2G
GL8), **4.4–12** (SN 2-D LS4) — so 10x is approximately what makes the two
criteria CO-BIND (clear at the same outer). Nothing in the tree states the rule.
It is worth writing down; it is not a per-family calibration worth preserving.

See also [[lessons-L24]] (a solver's nesting shape is per-ENTRY-POINT, and CP/MoC
/diffusion each nest differently) and [[sn-multigroup-axis-structure]].
