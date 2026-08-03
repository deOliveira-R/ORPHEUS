---
name: nonsn-geometric-transform-census
description: Durable census of geometric transformations in ORPHEUS's non-SN solver families (moc/mc/cp/diffusion/homogeneous/derivations) — who has what, and the two algebraic identities that survive line drift.
metadata:
  type: project
---

Census taken 2026-08-03 (`refactor/operator-strategy-layers`, `b0a003b4`→`bfedc621`);
full audit at `scratch/geom_transforms_audit_cross_solver.md`. Line numbers there
are perishable — the claims below are algebraic/structural and are not.

**Why:** the point-symmetry machinery (`numerics/symmetry.py`,
`numerics/roots_of_unity.py`, and now `numerics/quadrature/rules_circle.py`) exists;
the recurring question is "which solver families re-implement it by hand?"

**How to apply:** treat this as the routing map for any cross-method symmetry
question. Re-derive `file:line` via Nexus/grep; do NOT re-derive the shape.

## The shape (durable)

- **Zero hand-built rotation/reflection MATRICES exist outside
  `orpheus/numerics/symmetry.py`.** Measured twice (open + close) over all of
  `moc mc cp diffusion homogeneous derivations` — every matrix-shape pattern
  (`[[c,-s`, `eye(3)-2`, `-2*np.outer`, `Matrix([[...cos`) returns ZERO hits.
  The transformations are all angle arithmetic, arithmetic sign, or modulo.
- **Trig is FIVE lines** in all non-SN production: `moc/quadrature.py` (the
  azimuthal set), `moc/geometry.py` (its cos/sin), `mc/solver.py` (direction
  sampling). Nothing else.
- **MoC owns the only two hand-rolled `_orbit_closure` clones** in the codebase:
  `_reflected_azi_index` (angular partner by `argmin`) and `_find_link` (spatial
  partner by `argmin`). Neither has `_orbit_closure`'s two guards (distance
  window; bijectivity/ERR-073). `_find_link` computes `best_dist` and discards it.
- **MC has NO reflection at all** — its only BC is periodic (`% pitch`), the
  torus deck map. The Householder `d − 2(d·n)n` appears exactly ONCE in the whole
  scope, as *prose* in `peierls_nystrom/origins/specular/r_matrix.py`, never coded.
- **CP's images are an ADDITION.** `gap_c = bnd_pos[i] + bnd_pos[j]` is the
  mirror `s → −s`; marked only by a comment. Its folds (`chord_half_lengths` +
  `2.0*`, `y ∈ [0,R]`) are done by *choosing a half domain up front*, so they are
  invisible to a grep for `abs(`.
- **diffusion/ and homogeneous/ have no geometric transformation.** Clean negative.
- **Four spellings of one Neumann-sum-over-the-reflection-orbit**: CP's
  `1/(1−P_inout)`, Peierls' `(I − TR)⁻¹`, the continuous-µ `1/(1−e^{−2aµ})`, and
  `billiard.py`'s named `T = (I−S)⁻¹`. Peierls' own docstring proves it equals
  CP's at rank 1. The unifying concept is named; the facade owns none of them.

## The two algebraic identities (will outlive every line number)

1. **MoC's azimuthal set is the odd `4n`-th roots of unity.**
   `φ_k = π(k+½)/n_azi = 2π(2k+1)/(4·n_azi)`, so it equals
   `roots_of_unity(2k+1, 4·n_azi)` = the **upper half of
   `periodic_trapezoid(2·n_azi, shift=STAGGERED)`**. Measured identical to
   ≤5e-16; only the weight normalisation differs (MoC `Σω=1` vs circle-rule
   `Σw=2π`, ratio exactly `π`).
2. **MoC's reflection partner has a closed form: `k ↦ n_azi − 1 − k`** (i.e.
   `[::-1]`). `_reflected_azi_index` rediscovers it by `argmin` on every call,
   and is invoked twice per track with identical arguments. The `numerics`
   route makes that mirror **bit-exact at every `n_azi`**; MoC's `linspace+cos`
   makes it bit-exact at **none**.

## The gate that cannot bite

`tests/moc/test_verification.py::test_reflective_links_form_cycles` asserts
"following `fwd_link` returns to a visited track". For ANY total map on a finite
set the orbit repeats within `n+1 < max_steps`, so the cycle assertion is
satisfiable by construction; only the separate `assert current >= 0` (dangling
link) can fail. It also ignores `fwd_link_fwd`. It cannot see a non-bijection, a
wrong partner, or a large link gap. Measured: the link map IS a bijection today,
but the geometric match is loose — the reflected track's entry point sits up to
**2.5 % of the pitch** from the exit point, unreported.

Related: [[quadrature-landscape-durable]] (MoC's unnamed `[0,π)` quotient +
Σω=1 vs SN's `[0,2π)` + Σw=4π — corroborated and now located).
