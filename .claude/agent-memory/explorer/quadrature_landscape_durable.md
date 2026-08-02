---
name: quadrature-landscape-durable
description: The durable SHAPE of ORPHEUS's quadrature landscape — which axes (range/spacing/rule/exactness-space/generation) have >=2 realizations, where each family lives, and the measured LS_N exactness over-claim.
metadata:
  type: project
---

Surveyed 2026-08-01 on `refactor/operator-strategy-layers`; full evidence in
`scratch/quadrature_landscape_survey.md` (that file is the transient artifact —
this note is the part that survives churn). Re-derive any `file:line` via Nexus.

**Why:** the campaign question is "which quadrature choices deserve first-class
abstraction?", decided by "does the choice have >=2 genuine realizations, so an
invariant can be tested ACROSS them?" See [[lessons-L16]] for the measurement
technique.

**How to apply:** when any quadrature/angular-closure question comes up, start
from this shape instead of re-deriving it.

## The five axes and their realization counts

The one expression `np.linspace(0, 2π, n_phi, endpoint=False)` + `w=2π/n_phi`
(the SN azimuth) fuses **four** choices, not three:

| axis | in-tree realizations | verdict |
|---|---|---|
| **RANGE** | 8 (`[-1,1]`, `[0,1]`, `[a,b]`, `[0,2π)`, **`[0,π)`**, `[0,π]`, `[0,∞)`, `S²`) | earns it — strongest case |
| **SPACING** | Gauss / Chebyshev-Gauss / Laguerre / equispaced (3 offsets) / tabulated TY / adaptive | earns it |
| **RULE on the CIRCLE** | **ONE** (periodic trapezoid) and mathematically should be — it IS the circle's Gauss rule | stays a property |
| **RULE on an INTERVAL** | 5, + Gauss-Lobatto studied-and-declined | earns it |
| **EXACTNESS SPACE** | **9 distinct spaces collapsed into one unlabelled `int`** | earns it; the axis with the most existing realizations AND the only measured falsehood |
| **NODE GENERATION** (algebraic vs trig-evaluated) | see the 6-way kind census below | earns it (#325) |

## Generation-KIND census (measured 2026-08-02, Q3 symmetry-exactness sweep)

The earlier "2 exact vs 2 evaluated" split was too coarse. The real taxonomy —
GROUP-ACTION (representative + exact ops) / IMPOSED-SYMMETRY / HYBRID / FORMULA /
TABULATED / EXTERNAL:

- **GROUP-ACTION: exactly ONE** — `numerics/roots_of_unity.py` (integer quarter-split +
  octant fold + sign flips; ZERO tolerances, both fixed points decided by integer
  arithmetic). **Still zero production consumers, still not exported** from
  `numerics/__init__.py`; the 3 intended repoints (`rules_product`, `moc/quadrature`,
  `derivations/discrete/sn/balance.py`) are all un-done and BLOCKED on the argsort-tie
  ruling (exactness manufactures the ties).
- **IMPOSED-SYMMETRY: two** — `GeneratingMeasure.gauss` (`(x−x[::-1])/2`, gated by the
  COMPUTED `is_symmetric`) and numpy's `leggauss` inside `rules_1d`. Note these are two
  DIFFERENT code paths for the same rule, deliberately 1-4 ULP apart (snapshots pinned
  to `leggauss`).
- **HYBRID: `level_symmetric_sn` is NOT "exact"** — the 8-fold sign replication IS an
  exact group action (so D_2h is bit-exact), but the orbit REPRESENTATIVE is FORMULA:
  `η` comes from the level table while `ξ = sqrt(sinθ²−η²)`, two arithmetic paths to the
  same value ⟹ the index-6 coordinate-PERMUTATION half of its O_h claim is only
  round-off-approximate.
- **`lebedev` is EXTERNAL/unverified in-tree**, not "exact": scipy's tables, no ORPHEUS
  group action applied, invariance checked only in the test suite.
- **FORMULA:** `rules_product` (azimuth), MoC (azimuth + per-track cos/sin), MC flights.
- **TABULATED:** the MoC Tabuchi-Yamamoto polar table only.

**The CHECKER has the same split** (`numerics/symmetry.py`): `_octahedral_ops` (48
signed-permutation matrices) and `_reflections`/`_inversion_op` are EXACT; but
`_cyclic_ops` (= the n-th roots of unity!), `_vertical_mirrors`, `_rotation_z`,
`_rotation_about_axis`, `_icosahedral_ops` are `cos`/`sin` FORMULA. ⟹ the `Dnh(n_phi)`
product-rule claim is the ONE where BOTH sides are trig evaluations, which is what the
`_NODE_WINDOW_FACTOR = 100` widening (match window `1e-11` vs weight window `1e-13`)
absorbs. Lebedev/LS-vs-O_h is the cleanest pairing.

**The invariance machinery's production reach is ~nil.** `is_invariant`'s only
production caller is `registry.AngularSymmetry.admits_symmetry`, reachable only through
`select_quadrature`, which has NO production caller (meshes are built by direct
`Quadrature.<factory>` + `SNMesh(...)`). `maximal_invariance_groups`, `orbit_certificate`
and `singular_set` are **test-only**. So every shipped `invariance_group` is a DECLARED
tag; the computed check runs in `tests/`, not at construction. (Two tags have shipped
false: the product rule's `SO2`, and `cartesian2d`'s `O_h` residual — both retired
2026-08-02; product is now `Dnh(n_phi)`.)

**Partner maps: 3 SEARCH vs 1 FORMULA.** `_find_reflections` (`quadrature/directional.py`,
O(N²) `argmin`, **no distance guard at all**), `moc/geometry._reflected_azi_index`
(`argmin` over φ, no guard — the answer is exactly `n_azi−1−a` by index arithmetic), and
`MOCMesh._find_link` (nearest track endpoint, `best_dist` never thresholded). The one
FORMULA is `Quadrature.gauss_legendre`'s `identity[::-1]`, legitimate ONLY because
numpy imposes the ±μ mirror. The checker's `_orbit_closure` is a 4th `argmin` search
(window `atol*100`). The tree already owns the FORMULA endpoint: `PermutationOperator`
(`numerics/operator.py`) — `argsort` inverse + `perm[perm]==arange` involution, exact.

## Durable structural facts

- **The RANGE is ALREADY a named polymorphic property — but only in derivations.**
  `derivations/continuous/peierls_nystrom/geometry.py` `angular_range` returns
  `(-1,1)` for slab-polar and `(0,π)` otherwise. `numerics/quadrature` has nothing
  equivalent; every range there is a literal.
- **MoC has already taken a symmetry quotient and never named it.** Its azimuth
  is `[0,π)` MIDPOINT with weights summing to **1**; SN's is `[0,2π)`
  LEFT-ENDPOINT with weights summing to `4π`. Two normalisation conventions, no
  crosswalk. MoC's rule is the periodic trapezoid on the QUOTIENT circle `R/πZ` —
  spectral only on π-periodic integrands (true only AFTER the forward/backward
  track pairing). Measured, stated nowhere.
- **The OFFSET does not change exactness on the circle** (measured: left-endpoint
  and midpoint both exact for `1<=|k|<=n-1`). It is a pure symmetry /
  degeneracy-avoidance choice, currently made 3 different ways with 3 different
  (or absent) justifications.
- **MoC is entirely outside the tagged system**: `MOCQuadrature` is a plain frozen
  dataclass — no `DiscreteMeasure`, no `invariance_group`, no
  `degree_of_exactness`, no registry entry. Its Tabuchi-Yamamoto table is the one
  exactness claim in the tree that is NOT a polynomial degree (it is a minimax fit
  on the Bickley `Ki₃` family) — the case that an `int` field structurally cannot
  carry.
- **`degree_of_exactness` has exactly ONE reader that computes with it**
  (`tests/numerics/test_spherical_harmonic_basis.py`, `L = deg // 2`); every other
  reference is a tag-pinning assert. The selector does NOT read the attribute — it
  calls a per-family INVERTED formula in `registry.py`, so the formula is
  duplicated (rule-side tag + registry-side inverse) with no consistency check.
- **`level_symmetric_sn` is an EQUAL-weight construction, not Carlson-Lathrop.**
  One distinct weight per rule for every N. Measured degree of exactness is
  **3 for every N** while the tag says `N-1`; degree 3 is free for ANY `O_h`
  orbit with `Σw=4π`, so the tag adds nothing even when right (N=4). Consequence:
  the discrete SH Gram is non-orthogonal at 25-45% for `L>=2`, and
  `quadrature.angular_frame(L).gram` assumes the ANALYTIC Gram `4π/(2ℓ+1)`.
  Lebedev and product are clean to 1e-14.
- **`min(2n_mu-1, n_phi-1)` for the product rule is CORRECT and SHARP** (measured
  on both branches). The `min` over "polynomials" and "trig polynomials" is
  legitimate because for `x^a y^b z^c` on `S²` the max azimuthal frequency is
  exactly `a+b <= d`. The defect is that the reconciling mapping is unstated.
- **Shipped-but-unconsumed second realizations already exist**: `measure.equispaced`
  (midpoint), `measure.gauss_chebyshev`, `derivations.gauss_laguerre` — all 0
  production consumers. They are the >=2 the criterion asks for; nothing tests an
  invariant across them.
- **Gauss-Lobatto**: fully characterised, measured (~1.2x error penalty, same
  `N→∞` limit), and **architecturally DECLINED** — because `τ_raw=0` breaks the
  Morel-Montry recurrence AND the R12a `τ_raw ∈ (0,1)` predicate. That leak is the
  proof the RULE choice is not encapsulated.
- **`#123` is the project's own ratification of the >=2 criterion**: it mandates a
  ">= 2-quadrature signed-error stability protocol" for CP rank-N closure claims.

## Live tracked items on these axes

`#325` (L0, generation/spacing), `#326` (L1 bug, RANGE — the cylindrical level
double-covers `[0,2π)`), `#128` (QMC/Sobol' — a 4th kind of exactness space:
Hardy-Krause variation), `#109` (τ-coordinate + Gauss-Laguerre + Gauss-Jacobi),
`#123` (the >=2 protocol), `#265` (quadrature laws as invariant value objects),
`#235` (2-D (η,φ) cylindrical closure). `D_6h`/hex is a DOCUMENTED deferral with
the group already implemented (`numerics/symmetry.py` `Dnh(6)`, "forward-looking
— not yet consumed") and **no** `D_6h`-invariant rule in tree.

## Doc tension to expect

`docs/theory/methods/sn/curvilinear_one_group.rst` `sn-direct-seed-circle-vs-interval`
argues the FULL circle is load-bearing ("a periodic redistribution axis gives
edge-inclusion for free"). #326 says that full circle is a double cover. That
section needs rewriting whichever way #326 lands.
