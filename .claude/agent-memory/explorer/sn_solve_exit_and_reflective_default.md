---
name: sn-solve-exit-and-reflective-default
description: The SN solve exit is THREE Solution-construction sites (not the one its docstring claims); an unset-BC Mesh2D/axis-tuple resolves to ALL-REFLECTIVE at SNMesh level; and the tangential-ordinate bucket is where a trace-space kernel becomes visible.
metadata:
  type: reference
---

Three durable structural facts about the SN return path and its boundary
conditions. Measured 2026-08-14 on `refactor/track-b-remainder` @ `5c39226b`;
line numbers are re-derive-via-Nexus, the SHAPE is what survives.

## 1. The solve exit is THREE construction sites, not one

`orpheus/sn/solver.py` builds a `Solution`/`AdjointSolution` at **four** return
statements funnelling into **three** distinct paths:

| path | site | who reaches it |
|---|---|---|
| `_package_solution` | `return _package_solution(` ×2 | `solve_sn` (eigenvalue) AND `_package_adjoint_solution` ← BOTH adjoint entries |
| direct `Solution(...)` | `return Solution(` in `_solve_fixed_source_si` | `solve_sn_fixed_source` (SI arm) |
| direct `Solution(...)` | `return Solution(` in `_solve_fixed_source_krylov` | `solve_sn_fixed_source` (Krylov arm) |

⛔ `_package_solution`'s own docstring says *"The ONE `SolutionBase`
construction convention … The single boundary where converged iterates become
the typed return"*. **Present-tense FALSE** — the forward fixed-source family
bypasses it entirely. Any change to "what the solver returns" that is installed
at `_package_solution` silently misses `solve_sn_fixed_source`.

Carrier asymmetry across the three, also load-bearing:
- `_package_solution` wraps `_cell_average_angular(...)` — the AVERAGE_MOMENT
  slot only — into `TimedFullField(_history=(), history_depth=2)`.
- the SI arm returns the raw iterate composite (`angular_out`), moment tail and all.
- the Krylov arm returns `psi_full`, whose boundary block is documented as
  *"the matvec's B1'' face residual"* — a defect, not a flux trace.

## 2. Unset BC ⟹ ALL-REFLECTIVE, and it is the SNMesh's default, not the entry's

`resolve_boundary_conditions` (`orpheus/transport/method.py`) carries
`default = BC("reflective")` — the infinite-lattice/eigenvalue convention — and
it runs unconditionally inside `SNMesh.__init__`. Consequences, all `[M]`:

- `Mesh2D(edges_x, edges_y, mat_map)` with **no** `bc_*` arguments →
  `{xmin,xmax,ymin,ymax} = ReflectiveBoundary` on **every** route: the
  eigenvalue entries (`solve_sn` / `solve_sn_adjoint` pass
  `boundary_condition=None` to `_as_sn_mesh`), the axis-tuple d=3 entry, AND
  every direct `SNMesh(mesh, quad, mats)` in a test. The canonical k_inf
  lattice fixture is written with zero BC declarations.
- the two fixed-source entries default `boundary_condition="vacuum"`, but
  `_apply_default_bcs` fills faces **only when ALL are None**. ⚠ So declaring
  ONE face and leaving the rest unset gives the *reflective* default on the
  rest, not the entry's: `Mesh2D(..., bc_xmin=BC("vacuum"))` +
  `boundary_condition="vacuum"` measures
  `{xmin: VacuumInflow, xmax/ymin/ymax: ReflectiveBoundary}`.
- default scheme is `DiamondDifference`. ⟹ "unset BCs + default scheme + Mesh2D"
  IS the all-reflective-Cartesian-DD configuration.

## 3. The tangential ordinate bucket is the third bucket, and it is where a trace kernel shows

`AngularTraceSpace.inflow_indices_for_face` selects `< -TANGENTIAL_EPS`,
`outflow_indices_for_face` selects `> +TANGENTIAL_EPS`; the module comment says
tangential ordinates (`|Ω·n| <= eps`) *"belong to `"full"` alone"*. So any
functional summed over a **signed half-range** (half-range φ±, J±) misses them,
while any functional summed over the **full face** (`Σ_n w_n ψ_n` over all
ordinates, `Σ|trace|`) includes them. That is exactly why an `|Ω·n|^0` face
moment can see a trace-space null direction that every half-range functional
annihilates — and why the split only bites on quadratures that HAVE tangential
ordinates (`product`, `lebedev`), not on `level_symmetric`. (L-010's
non-exhaustive-predicate lesson, with a physical consequence attached.)

## 4. Zero production consumers of the returned trace

`grep "\.boundary_flux\|angular_flux\.boundary" orpheus/ examples/` returns
only `orpheus/sn/solution.py` (the property definition) plus two docstring
mentions of a retired module path. Nothing inside the package reads the trace
a solve returns. `Solution.compare()` diffs `angular_flux.interior.values` —
**bulk only** — so it is trace-blind too.
