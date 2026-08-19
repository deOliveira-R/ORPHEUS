---
name: nexus82-cp-implementers
description: CP equation→code truth (nexus#82) — the CP method has exactly TWO implementations, the three `_*_cp_matrix` are zero-arithmetic facades, and the page's `vv-status rationale` comments did NOT answer any of the 15 target equations (the TEST pytestmark comments did)
metadata:
  type: project
---

Third in the equation→implementer series after [[nexus82-operator-algebra-implementers]]
and [[nexus82-loss-representation-implementers]]. Derived 2026-08-18 from the tree,
not from `.claude/inventories/implements_declaration_inventory.md` (refuted on this
population). Full per-equation output: `.claude/inventories/declare_cp.md`.

**Why:** ORPHEUS is building a licensing-grade V&V ledger; an equation with no declared
implementer makes every coverage claim on it unadjudicable. The 15 CP equations carried
**207** claims between them.

## ⭐ The method correction — the brief's stated method did NOT work here

The dispatch brief said the authored `.. (vv-status rationale)` comment above/below each
`.. math::` "is usually the whole answer" (9 of 11 classified that way on a prior page).
`[M]` on `docs/theory/methods/collision_probability.rst`: **6** rationale comments exist
on 2296 lines, and **0 of 15** target equations has one. The page's rationale density is a
page-by-page property, not a corpus property — do not budget on the prior page's rate.

**What answered instead, in descending yield:**

1. **Inline prose naming the method with a `:meth:` role + a quoted code line.** This page
   does it constantly and precisely: *"In the code: ``P_out = 1 - P_cell.sum(axis=1)``
   (:meth:`CPMesh._apply_white_bc`)"*, *"character-for-character, the solver line
   ``rcp[i, i] += sti * t[i] - (0.5 - _e3(tau_i))``"*. Grep the page for `:meth:` and for
   `` :: `` literal blocks before anything else.
2. ⭐ **The test module's own `pytestmark = pytest.mark.verifies(...)` COMMENTS.** These are
   authored, greppable, and frequently name the exact symbol whose breakage reds the gate —
   e.g. `tests/cp/test_properties.py:36-38`: *"test_reciprocity also pins the rearranged
   lower-triangle form the code uses to fill P_ji (**a factor error in `_normalize_rcp`
   breaks the identity**)"*. That single comment settled `reciprocity-lower-triangle`.
   **This surface is not in the brief's method and is the one I would reach for first next
   time.** Find it with `grep -rn -- "\"<label>\"" tests/ orpheus/`.
3. **A neighbouring equation's own label** as a *scoping* answer — see the geometry
   partition below.

## ⭐ Measured: all 15 had ZERO `implements` edges — not even guesses

`[M]` 2026-08-18 against `.nexus/graph.db` (13 206 `implements` edges total): every one of
the 15 had **no** incoming `implements` edge of any provenance. Their labels (`p-inf`,
`tau-m`, `self-double-integral`) share no name token with any symbol, so the heuristic
never fired.

⟹ **The "declaring one stands the guess down for the whole equation" hazard is
population-dependent — measure it before letting it shape the answer.** One graph query
(`g.in_edges(<eq node>)` filtered to `type=='implements'`) tells you whether omission
costs a guess or costs nothing. Where it costs nothing, the risk calculus flips: breadth
is cheap and the only real cost of omission is a claim left unadjudicable.

## The durable CP shape (survives line-number churn)

- **Exactly TWO implementations of the flat-source CP construction**, structurally
  parallel line for line:
  1. **Production** — `cp.solver.CPMesh`, a 4-step pipeline: `_compute_slab_rcp` (E₃,
     scalar) **or** `_compute_radial_rcp` (cyl+sph, y-quadrature) → `_normalize_rcp` →
     `BC_REGISTRY`-selected closure (`_apply_white_bc` / the vacuum lambda).
     `compute_pinf_group` is the orchestrator running those in order.
  2. **Derivations** — `derivations.continuous.flat_source_cp.geometry.build_cp_matrix`,
     ONE geometry-polymorphic body for all three geometries (slab = the degenerate 1-point
     y-rule at `y=0, w=0.25`), with rcp assembly + normalisation + white-BC closure all
     inlined.
- ⛔ **`slab._slab_cp_matrix` / `cylinder._cylinder_cp_matrix` / `sphere._sphere_cp_matrix`
  are post-Phase-B.2 thin facades with ZERO arithmetic** — each pre-binds a
  `FlatSourceCPGeometry` singleton and returns `build_cp_matrix(...)`. The theory page
  still cites them as carrying the formulas (`P_inf_g[:,:,g] = P_cell + np.outer(...)`,
  `P_in = sig_t_g * t_arr * P_out`, `S_cell = 2.0 * np.pi * r_cell`) and says the formula
  appears "independently in all three derivation scripts". `[M]` grep for those lines in
  the three files → **0 hits**. There is ONE derivations implementation, not three; the
  independence the page advertises is gone.
- **The reduced CP `rcp[i,j] = Σ_t,i V_i P_ij` is a materialized array in three functions
  and is symmetric BY CONSTRUCTION** — that is what makes `reciprocity` declarable rather
  than a bare law. `_normalize_rcp` (row-divide by `Σ_t,i V_i`) is the site whose factor
  can break it, which the test comment names.
- ⛔ **`reciprocity-lower-triangle`'s "This halves the computation cost" is
  present-tense-false.** No lower-triangle fill-in exists; all three implementations run
  the full N² double loop.
- **`P_inf[:, :, g]` appears 4× in `cp/solver.py`**: 3 transposed contractions
  (`_solve_fixed_source_jacobi`, `_solve_fixed_source_gs`, `_compute_balance_residual`) +
  1 assembly assignment in `solve_cp`. The GS arm and the residual diagnostic are the two
  the theory page does NOT name — a page-driven declaration misses both.
- **`CPSolver._compute_balance_residual` forms BOTH sides of `neutron-balance`** and
  returns their norm; it exists for no other purpose. When hunting an equation's
  implementer, grep for a `*_residual` method before anything else — a residual routine is
  the equation written out twice.

## Page-scoping rules this page taught (generalize to any multi-geometry method page)

- **A chapter's opening sentence scopes its labels.** "Optical Path Construction Along a
  Chord" opens *"the most geometrically intricate part of the **cylindrical and spherical**
  CP implementations"* — so `tau-m`/`tau-p` are radial-only, and the slab's structurally
  identical `gap_d`/`gap_c` belongs to its OWN labels (`cumulative-optical-path-slab`,
  `dd-slab`, `dc-slab`). Declaring the slab site would double-cover a construction the page
  deliberately partitions.
- **Before declaring a general primitive, grep the corpus for a label that already owns
  it.** `kernels.e_n`/`e_n_mp` look like `e3-def` implementers; the general `E_n` has its
  own labels on another page (`en-definition` + 3 siblings,
  `docs/theory/verification/reference_solutions.rst:232-255`) and its consumers are
  Peierls-Nyström, not CP. Same move saved `_second_difference` from being rehomed — its
  native label `cp-second-difference-operator` is real and `tested`
  (`references/peierls_nystrom.rst:6071`).

## The directive form (read from the shipped precedent, not recalled)

`docs/theory/methods/sn/loss_representation.rst:219` — one block per implementer:

```rst
.. implements:: <equation-label>
   :by: orpheus.dotted.path.To.Symbol

   Body prose saying why this IS the implementation.
```

Negative form (landed at `bb075c93`, 11 equations, all on `operator_algebra.rst` +
`sn/loss_representation.rst`): `.. no-implementation:: <label>` with
`:kind: identity|law|canonical-form|definition`.
