# Issue #325 — blast-radius audit

**Scope**: replacing trig-evaluated circular node sets with an algebraically
generated (group-action) primitive.

**Audit run**: 2026-08-01, branch `refactor/operator-strategy-layers`,
HEAD `1659d756`, working tree clean of production edits. Nexus graph active
workspace `/Users/rodrigo/git/nuclear/ORPHEUS`, built from `cea73f8d` (2
commits behind HEAD; both are boundary/plan commits, no quadrature touch).
numpy 2.4.4, Python 3.14 (`.venv`).

`[M]` = measured, verbatim output pasted.

---

## §0 HEADLINE

### 0.0 The one-paragraph verdict

**The VALUE blast radius is tiny and fully measured; the ANALYSIS blast radius
is larger than the issue describes.** Running the whole product-consuming
surface (three batches: 3024 tests counted plus a third full sweep over
`tests/sn`+`tests/derivations`+`tests/transport`) with the algebraic primitive
swapped in produces exactly **one test failure** — and it is not a value re-baseline, it is
a test whose fixture depends on the defect existing — plus exactly **one**
stored-snapshot value movement (`cyl_1g_homogeneous_product_dd_n20`, 1–2 ULP,
inside its own correctness gate, `DriftWarning` only). What IS larger than the
issue says: (a) the intra-level `argsort` order changes in 24/24 configurations
(§0.1) — though §9 proves this is FP-only, not a change of recurrence; (b) the
tie-resolution non-determinism it exposes **already exists today on
`level_symmetric`** (§6), which is a live latent reproducibility hole
independent of #325; (c) the MoC consumer row is incomplete — the reflection
site the issue does not list needs a different fix (§3.c); (d) the doc blast
radius is real, including a **provably vacuous test guard** (§5.1) and a
**measured doc table that #325 falsifies** (§5.3).

**Recommended framing:** #325 is a principled re-baseline whose numerical risk
is negligible and whose *design* work (tie-break policy, MoC reflection, doc
re-framing, the vacuous eps guard) is the real content. Fourteen additions to
the acceptance criteria in §7.

### Three claims in the issue body that do not survive measurement

### 0.1 ⛔ "Nothing changes structurally: same count, same ordering, same weights" — FALSE

`rules_product.py:139` sorts each μ-level by the direction cosine `η = μ_x`:

```python
order = np.argsort(mu_x[level_arr])
```

Today the azimuthal mirror pair `(k, n−k)` produces `cos φ` values that differ
in the last bits, so `μ_x` has a strict total order within a level and the
sort is unambiguous. Under #325 those values become **bit-equal** — an exact
TIE — and tie resolution changes the emitted order. `[M]`:

```
 n_mu  n_phi |  OLD ties  NEW ties |  order CHANGED?
    2      3 |         0         2 |            True
    2      4 |         0         2 |            True
    2      6 |         0         4 |            True
    2      8 |         0         6 |            True
    2     12 |         0        10 |            True
    2     16 |         0        14 |            True
    2     24 |         0        22 |            True
    2     32 |         6        30 |            True
    4      3 |         0         4 |            True
    4      4 |         0         4 |            True
    4      6 |         0         8 |            True
    4      8 |         0        12 |            True
    4     12 |         0        20 |            True
    4     16 |         0        28 |            True
    4     24 |         0        44 |            True
    4     32 |         8        60 |            True
    8      3 |         0         8 |            True
    8      4 |         0         8 |            True
    8      6 |         0        16 |            True
    8      8 |         0        24 |            True
    8     12 |         0        40 |            True
    8     16 |         4        56 |            True
    8     24 |         0        88 |            True
    8     32 |        22       120 |            True

  level order CHANGED in 24/24 configurations
```

The `n_mu=2, n_phi=8` detail, level 0, at 20 significant digits `[M]`:

```
  idx  OLD mu_x                     NEW mu_x
    0  8.16496580927726034460e-01  8.16496580927726034460e-01
    1  5.77350269189625842081e-01  5.77350269189625842081e-01
    2  4.99959962173948737480e-17  -0.00000000000000000000e+00
    3  -5.77350269189625731059e-01  -5.77350269189625842081e-01
    4  -8.16496580927726034460e-01  -8.16496580927726034460e-01
    5  -5.77350269189625842081e-01  -5.77350269189625842081e-01
    6  -1.49987988652184615081e-16  0.00000000000000000000e+00
    7  5.77350269189625620037e-01  5.77350269189625842081e-01

  OLD level_indices[0] = [4, 5, 3, 6, 2, 7, 1, 0]
  NEW level_indices[0] = [4, 3, 5, 2, 6, 1, 7, 0]
```

Independently reproduced. The coordinator's finding stands. **Note the extra
detail:** `n_phi=32` ALREADY has ties today (6 / 8 / 22 depending on `n_mu`) —
so the "unambiguous total order" premise is already false at `n_phi=32` on the
CURRENT tree. #325 does not create the hazard; it universalizes it.

### 0.2 ⛔ `np.argsort`'s default tie resolution is NOT deterministic across kinds

`[M]`, on the post-#325 (tied) `μ_x` array:

```
  n_phi=  4: all-kinds-agree=True  quicksort==stable=True
  n_phi=  8: all-kinds-agree=True  quicksort==stable=True
  n_phi= 12: all-kinds-agree=True  quicksort==stable=True
  n_phi= 16: all-kinds-agree=True  quicksort==stable=True
  n_phi= 24: all-kinds-agree=False  quicksort==stable=False
       quicksort: [12, 11, 13, 14, 10, 9, 15, 8, 16, 7, 17, 6, 18, 5, 19, 4, 20, 3, 21, 22, 2, 23, 1, 0]
       mergesort: [12, 11, 13, 10, 14, 9, 15, 8, 16, 7, 17, 6, 18, 5, 19, 4, 20, 3, 21, 2, 22, 1, 23, 0]
          stable: [12, 11, 13, 10, 14, 9, 15, 8, 16, 7, 17, 6, 18, 5, 19, 4, 20, 3, 21, 2, 22, 1, 23, 0]
        heapsort: [12, 11, 13, 14, 10, 9, 15, 8, 16, 7, 17, 6, 18, 5, 19, 4, 20, 3, 21, 22, 2, 23, 1, 0]
  n_phi= 32: all-kinds-agree=False  quicksort==stable=False
       quicksort: [16, 15, 17, 18, 14, 13, 19, 20, 12, 21, 11, 22, 10, 9, 23, 8, 24, 25, 7, 26, 6, 27, 5, 28, 4, 29, 3, 2, 30, 31, 1, 0]
       mergesort: [16, 15, 17, 14, 18, 13, 19, 12, 20, 11, 21, 10, 22, 9, 23, 8, 24, 7, 25, 6, 26, 5, 27, 4, 28, 3, 29, 2, 30, 1, 31, 0]
          stable: [16, 15, 17, 14, 18, 13, 19, 12, 20, 11, 21, 10, 22, 9, 23, 8, 24, 7, 25, 6, 26, 5, 27, 4, 28, 3, 29, 2, 30, 1, 31, 0]
        heapsort: [16, 15, 17, 18, 14, 13, 19, 20, 12, 21, 11, 22, 10, 9, 23, 8, 24, 25, 7, 26, 6, 27, 5, 28, 4, 29, 3, 2, 30, 31, 1, 0]
```

`np.argsort` default `kind="quicksort"` is an *introsort* — its tie resolution
is an implementation detail with no cross-version / cross-platform guarantee.
**Verdict:** this is a latent reproducibility bug that already exists at
`n_phi = 32` (where ties already occur today) and that #325 would extend to
every `n_phi`. It is fixable in the same change — `kind="stable"` plus a
tie-breaker that is *physically* meaningful (see §0.3) — but the fix is not
in the issue's acceptance criteria and must be added.

### 0.3 ⛔ The sort KEY itself is defective: `η = sinθ·cos φ` is 2-to-1 over `φ ∈ [0, 2π)`

`η` cannot order a full azimuthal circle: `φ` and `2π − φ` share `η`. `[M]`,
post-#325 (`n_mu = 2`, level 0), showing `η` in emitted order and the `ξ = μ_y`
that follows it:

```
  n_phi=4:   eta = [-0.8165 -0.      0.      0.8165]
             xi  = [-0.      0.8165 -0.8165  0.    ]
             n_unique_eta = 3 / 4
  n_phi=6:   eta = [-0.8165 -0.4082 -0.4082  0.4082  0.4082  0.8165]
             xi  = [-0.      0.7071 -0.7071  0.7071 -0.7071  0.    ]
             n_unique_eta = 4 / 6
  n_phi=8:   eta = [-0.8165 -0.5774 -0.5774 -0.      0.      0.5774  0.5774  0.8165]
             xi  = [-0.      0.5774 -0.5774  0.8165 -0.8165  0.5774 -0.5774  0.    ]
             n_unique_eta = 5 / 8
  n_phi=12:  eta = [-0.8165 -0.7071 -0.7071 -0.4082 -0.4082 -0.      0.      0.4082  0.4082  0.7071  0.7071  0.8165]
             xi  = [-0.      0.4082 -0.4082  0.7071 -0.7071  0.8165 -0.8165  0.7071 -0.7071  0.4082 -0.4082  0.    ]
             n_unique_eta = 7 / 12
```

Sorting the full circle by `η` **interleaves the two half-circles** (`ξ > 0`
and `ξ < 0` alternate). The order is `η`-monotone, which is what the Bailey
α-recursion needs (see §A.3), but the ±ξ interleave within each `η` tie is
arbitrary in BOTH the old and the new code. Today rounding noise picks it;
after #325 the sort algorithm picks it. Neither is a decision the physics made.

### 0.4 The "≤ 1 ULP" claim is true in ABSOLUTE terms, false in ULP-DISTANCE terms

`[M]` — comparing the current `np.cos(2πk/n)` against the proposed primitive:

```
 n_phi  max_ulp_cos  max_ulp_sin   max_abs_delta
     3            5            3    4.996004e-16
     4 4371439675322042890 4368955796522032135    1.836970e-16
     5            3            1    1.665335e-16
     6           11 4368955796522032135    6.106227e-16
     8 4371439675322042890 4368955796522032135    2.220446e-16
    10            3 4368955796522032135    2.220446e-16
    12 4371439675322042890 4368955796522032135    6.106227e-16
    16 4371439675322042890 4368955796522032135    6.106227e-16
    20 4371439675322042890 4368955796522032135    2.775558e-16
    24 4371439675322042890 4368955796522032135    8.326673e-16
    32 4371439675322042890 4368955796522032135    6.106227e-16
    48 4371439675322042890 4368955796522032135    9.992007e-16
    64 4371439675322042890 4368955796522032135    6.106227e-16
```

The 4.37e18 figures are the on-axis nodes crossing to **exactly `0.0`**: there
really are ~4.4e18 representable doubles between `6.12e-17` and `0.0`. The
issue's "≤ 1 ULP" is the *absolute* delta (≤ 9.99e-16 measured, slightly above
the 8.3e-16 the issue quotes — I see 9.99e-16 at `n_phi = 48`), i.e. ~1 ULP of
a value near 1.

**This distinction bites a real gate.** `tests/sn/regression/_regression_assert.py`
implements the ULP-distance metric with exactly the monotone-int mapping used
above, and `kind="direct"` calls `np.testing.assert_array_almost_equal_nulp`.
Any regression output element that crosses to (or away from) exactly zero
blows the nulp budget catastrophically, not marginally. See §2.3.

### 0.5 Correction to the issue's blast-radius count

The issue body says *"76 test files and 23 production sites touch `product`."*
`[M]` That number is wrong — it counts the token `product`, which is
massively overloaded in this tree (`OperatorProduct`, `TensorProductOperator`,
`itertools.product`, `product_state`, "production", "dot product", "inner
product"). Measured properly:

| population | measured |
|---|---|
| test files calling `.product(` | **50** |
| test call sites `.product(` | **100** |
| production `Quadrature.product(` call sites | **3** (all in `orpheus/derivations/`) |
| production modules importing `product_mu_phi` | **4** (`numerics/__init__`, `quadrature/__init__`, `quadrature/directional`, `quadrature/registry`) |

The three production `Quadrature.product(` sites `[M]`:

```
orpheus/derivations/discrete/sn/contamination.py:224:        q = Quadrature.product(n_mu=4, n_phi=n_phi)
orpheus/derivations/continuous/mms/sn.py:2094:    quadrature = Quadrature.product(n_mu=n_mu, n_phi=n_phi)
orpheus/derivations/continuous/mms/sn.py:3847:    quadrature = Quadrature.product(n_mu=n_mu, n_phi=n_phi)
```

There is **no** production consumer outside `orpheus/derivations/` — the
solver never selects `product` for itself. `select_quadrature` (the registry
route that CAN return `ProductQuadrature`) has **zero production callers**;
its only callers are `tests/numerics/test_registry.py`. A false positive worth
naming: `tests/numerics/data/step6_product_reroute_baseline.npz` is about
`OperatorProduct.solve`, not the product quadrature — it holds
`DiagonalOperator`/`PermutationOperator` results on a seeded random vector and
is completely immune.

---

## §A — The widened question: IS the α-recursion order-sensitive in fact?

Measured, not reasoned. Four variants of the SAME cylinder SN solve, isolating
VALUES from ORDER — `A` = old nodes/old order (production today), `B` = new
nodes/new order (#325 as proposed), `C` = **old nodes / NEW order** (isolates
the pure ORDER effect), `D` = new nodes / OLD order (isolates the pure VALUE
effect). `[M]`:

```
==================================================================================
CYL 2G 3-region reflective n=40  product(2x4)
==================================================================================
  variant      k_eff                       dk vs A      rel_dk    rel dphi    ULP(phi)
  A old/old    1.2281570659809602        0.000e+00   0.000e+00   0.000e+00           0
  B new/new    1.2281570659809606        4.441e-16   3.616e-16   1.413e-15          11
  C old/NEW    1.2281570659809604        2.220e-16   1.808e-16   8.977e-16           8
  D new/OLD    1.2281570659809609        6.661e-16   5.424e-16   1.321e-15          10

==================================================================================
CYL 2G 3-region reflective n=40  product(4x12)
==================================================================================
  variant      k_eff                       dk vs A      rel_dk    rel dphi    ULP(phi)
  A old/old    1.229241194323012         0.000e+00   0.000e+00   0.000e+00           0
  B new/new    1.2292411943230113       -6.661e-16  -5.419e-16   5.764e-15          49
  C old/NEW    1.2292411943230115       -4.441e-16  -3.613e-16   3.756e-15          27
  D new/OLD    1.2292411943230113       -6.661e-16  -5.419e-16   5.704e-15          50

==================================================================================
CYL 1G homogeneous VACUUM R=8cm n=40  product(2x4)
==================================================================================
  variant      k_eff                       dk vs A      rel_dk    rel dphi    ULP(phi)
  A old/old    1.43283764876626          0.000e+00   0.000e+00   0.000e+00           0
  B new/new    1.43283764876626          0.000e+00   0.000e+00   1.744e-15          13
  C old/NEW    1.4328376487662597       -2.220e-16  -1.550e-16   4.278e-16           3
  D new/OLD    1.43283764876626          0.000e+00   0.000e+00   1.283e-15           9

==================================================================================
CYL 1G homogeneous VACUUM R=8cm n=40  product(4x8)
==================================================================================
  variant      k_eff                       dk vs A      rel_dk    rel dphi    ULP(phi)
  A old/old    1.4304706910361304        0.000e+00   0.000e+00   0.000e+00           0
  B new/new    1.43047069103613         -4.441e-16  -3.104e-16   7.767e-16           5
  C old/NEW    1.4304706910361302       -2.220e-16  -1.552e-16   9.215e-16           7
  D new/OLD    1.4304706910361302       -2.220e-16  -1.552e-16   6.214e-16           4
```

(A `product(4x8)` 3-region row and the `4x12` vacuum row, omitted for length,
sit in the same band; the homogeneous-reflective k_inf cases are degenerate —
they return 1.5 to 4e-16 regardless of variant, so they carry no signal.)

### VERDICT on §A — reassuring, but with a caveat that must be stated

**The pure ORDER permutation (variant C) moves the converged answer by at most
`4.4e-16` relative in `k_eff` and `5.8e-15` relative (49 ULP) in the scalar
flux.** That is round-off scale. The α-recursion is order-sensitive only
through FP non-associativity, not structurally: the level permutation does NOT
attach a different α to a given ordinate in exact arithmetic, it only changes
the summation order of the running sum.

So #325 does **not** reclassify to "changes a production recurrence to a
different answer". It stays a re-baseline. **But three caveats are load-bearing:**

1. The variant-C column is not zero. The permutation is a real perturbation at
   the ULP scale, and it is **the same size as the value shift** (compare C and
   D columns — they are indistinguishable). Any re-baseline narrative that says
   "the values moved ≤1 ULP" is telling half the story; the *order* moved too,
   and both contribute comparably.
2. The measurement covers CYL 1-D only (the geometry that consumes
   `level_indices`). It does not cover a future d-generic or 2-D consumer of
   `LevelStructure`.
3. The tie-resolution non-determinism (§0.2) means the ORDER itself is not a
   reproducible function of the inputs at `n_phi ∈ {24, 32}` unless
   `kind="stable"` is pinned. That is the item to fix, and it is what makes the
   ≤50-ULP order effect *unreproducible* rather than merely small.

### Who reads the per-level ordinate ORDER

`[M]` grep + Nexus. Production consumers of `level_indices` / `level_mu` /
`pole_angular_closure.level_indices`:

| site | what it does with the order |
|---|---|
| `orpheus/transport/radial_characteristic_field.py:289,292,313` | per-μ-level ordinate slice for the ψ½ carrier |
| `orpheus/sn/operators/radial_characteristic.py:780,789,852,859` | `A_BB` — the level loop feeding the ψ½ march |
| `orpheus/sn/sweep/pole_angular_closure.py:595,893,904,926,1019,1070,1178,1184,1278` | **the α-recursion itself** (`numer_bar` / `alpha_half` built per level, in level order) |
| `orpheus/sn/loss_representation/__init__.py:2827,2902,3518,4089-4091,4155,4609-4611,4621` | the curvilinear walk + matvec level loops |
| `orpheus/sn/mesh/augmented_mesh.py:1261,1317,1468` | `(mu_level, ordinate_idx) → global_n` index translation |
| `orpheus/sn/sweep/cache.py:247-252` | per-level ordinate cache |
| `orpheus/geometry/reduced_operator.py:208,498,576-582,779,791,807` | `alpha_half` / `redist_dAw` construction (the Bailey recursion's storage) |
| `orpheus/derivations/discrete/sn/contamination.py:135,189` | derivation-side level loop |

`pole_angular_closure.py` is the α-recursion's home and the order-sensitive one.
Note `rules_sphere.py:215` applies the **same** `argsort` pattern to the
level-symmetric family — so any fix to the sort must be applied there too, or
the two families diverge in convention.

---

## §1 — The `product` consumer inventory, split by sensitivity

### 1.a Bit-exact pins — measured EMPIRICALLY, not by grep

A grep for `assert_array_equal` in the 50 product-consuming test files returns
~200 hits, and **almost all of them are immune**: they are *route-equivalence*
assertions (route A vs route B computed in the same process from the same
nodes), where both sides move together. Classifying them by reading would be
error-prone, so I measured instead.

**Method `[M]`:** a pytest plugin swapped `product_mu_phi` for the algebraic
`circle_nodes` implementation from the issue body (patching the definition
module, `directional`, `registry`, and both `__init__` re-exports plus the
registry spec's `factory`), then ran all 50 product-consuming test files:

```
1 failed, 1283 passed, 13 deselected, 61 xfailed, 10 warnings in 265.17s (0:04:25)
```

**Exactly ONE test breaks:**

```
FAILED tests/sn/operators/test_angular_average_operator.py::TestOutflowHalfTraceIsTheDomain::test_domain_is_the_tangential_band_gamma_plus
  - AssertionError: fixture no longer discriminates: Γ₊ [0, 4] and the retired '> 0.0' set [0, 4] agree,
    so this gate cannot see the classifier the narrowing retired.
```

This is **not** a value re-baseline. It is a test whose fixture depends on the
defect existing: `product(2, 4)` was chosen precisely because it carries
small-but-nonzero tangential cosines, so the retired `> 0.0` classifier and the
`> TANGENTIAL_EPS` one disagree. #325 removes the discriminator and the test's
own ACTIVATION guard fires — exactly as designed (`vv` Mode-8 class 7). See §6
for the migration options; this needs a `test-architect` decision, not a
re-baseline.

**No bit-exact pin against a frozen constant moved.** Zero `assert_array_equal`
failures, zero snapshot failures in that batch.

### 1.b Tight-tolerance pins — none bit

`[M]` No `rtol`/`atol` gate at or below 1e-14 in the 50 files failed. The
tightest gates in the batch that *do* consume product nodes and stayed green
include `tests/sn/primitives/test_quadrature.py:375` (`rtol=1e-14` on `delta_A`)
and `tests/sn/operators/test_angular_average_operator.py:157`
(`rtol=1e-14` on the conservation row).

### 1.c Loose — the remainder

1283 tests in the batch, 1282 green under the swap.

### 1.d Correction to the issue's population claim

See §0.5. **50** test files (not 76) and **3** production `Quadrature.product(`
call sites (not 23), all three inside `orpheus/derivations/`.


---

## §2 — The snapshot stores

`[M]` Full census: 33 `.npz` + 21 `.npy` + 2 `.json` stored baselines under
`tests/`. I traced the quadrature each one is built from, via its generator /
consumer.

| store | quadrature | product-dependent? | regenerates from | regeneration routes through the CHANGED code? |
|---|---|---|---|---|
| `tests/sn/regression/snapshots/cyl_1g_homogeneous_product_dd_n20.npz` | **`product(2,4)`** | ⚠ **YES — the only one** | `tests/sn/regression/_generate_snapshots.py` (`--case cyl_1g_homogeneous_product_dd_n20`) | ⚠ **YES** — the generator calls `Quadrature.product` + `solve_sn` |
| the other 12 `tests/sn/regression/snapshots/*_dd_*.npz` | GL-8 / LS4 / LS6 / Lebedev-5 | no | same generator + `_generate_2d_octant_snapshots.py` | n/a |
| `tests/sn/regression/snapshots/walk_matvec_{slab,sphere,cyl,cart2d}_2g.npz` | `_make_cyl` uses **`level_symmetric(4)`** (NOT product) | **no** | `_generate_walk_baselines.py` | n/a |
| `tests/geometry/snapshots/bc_equivalence_*.npz` (7) | Lebedev-17 ×3, LS4 ×2, LS6, GL-8 | **no** | `tests/geometry/_generate_bc_equivalence_snapshots.py` | n/a |
| `tests/sn/_data/bc_extraction_baseline/*.npy` (18) | GL (SLB/SPH), **LS4** (CYL) | **no** | self-regenerating on first run (writes if missing) | n/a |
| `tests/sn/_data/bc_extraction_2d_baseline/*.npy` (3) | LS4 | **no** | ditto | n/a |
| `tests/sn/_fixtures/wave_t_t3/pre_t3_snapshots.npz` | Lebedev-17 | **no** | `_capture_pre_t3_snapshots.py` | n/a |
| `tests/sn/_fixtures/wave_t_t4/pre_t4_snapshots.npz` | GL + LS | **no** | `_capture_pre_t4_snapshots.py` | n/a |
| `tests/sn/sweep_ref_2g.npy` | Lebedev-17 / GL-8 / LS4 | **no** | hand-captured | n/a |
| `tests/numerics/data/step6_product_reroute_baseline.npz` | **none** — it is `OperatorProduct`, not the product QUADRATURE | **no** (false positive on the token) | pre-carve capture script | n/a |
| `tests/_harness/pyright_baseline.json` | n/a | no | pyright | n/a |

### 2.1 ⚠ The one contaminated-re-baseline risk

`cyl_1g_homogeneous_product_dd_n20.npz` is the single product-derived stored
baseline, and its regeneration path is **exactly** the pattern the crosswalk
§16.1 warns about: `_generate_snapshots.py` builds the case with
`Quadrature.product(n_mu=2, n_phi=4)` and calls `solve_sn` — i.e. the generator
is a **recorder of production output**, not an independent reference. Re-running
the generator after #325 records whatever the new code produces; the resulting
assertion is `production == a recording of production`.

**However — and this materially reduces the risk — this generator already
carries the compensating discipline the BC-equivalence generator got at B3.4b.**
Its module docstring states the protocol explicitly (`_generate_snapshots.py:26-31`):

> When this script is run on a different solver state than the snapshots
> were generated against, the reproductions will mismatch — that is the
> intended behaviour. The protocol for legitimate updates: (1) audit why
> the new output is correct (with V&V evidence), (2) re-run the
> generator, (3) commit both the new snapshot AND the audit narrative
> in the same commit.

and `test_dd_regression.py:26-38` records the **independent corroboration** each
value carries: the `cyl_1g_homogeneous_*` rows are homogeneous-reflective, so
`k_eff` is pinned by the closed form `k_inf = νΣ_f/Σ_a = 1.5`, independent of
the quadrature. That closed-form anchor is what makes a regeneration honest here
— it is NOT a bare recording.

**Verified `[M]`:** all four variants of the cylinder homogeneous-reflective
solve return `k_eff = 1.5 ± 4.4e-16` (see §A), so the closed-form anchor holds
under the swap.

### 2.2 The gate that actually fires

`tests/sn/regression/_regression_assert.py` is a two-layer gate:

* Layer 1 **correctness** — `kind="iterative"` → `SAFETY(=10) × conv_tol`.
  `test_dd_regression.py` uses `kind="iterative"` for *every* row, so a ≤50-ULP
  drift passes with ~10 orders of margin.
* Layer 2 **drift tripwire** — a `DriftWarning` on ANY non-bit-identical value,
  escalatable with `-W error::DriftWarning`.

**So `cyl_1g_homogeneous_product_dd_n20` will NOT hard-fail on a plain
`pytest` run; it will emit a `DriftWarning`, and it WILL hard-fail under
`-W error::DriftWarning`** — which is the flag the project reserves for
"pure-refactor, zero-numerical-change" PRs. #325 is not such a PR, so the
correct handling is: regenerate that ONE snapshot, record the justification
(node values now exact-symmetric; k_eff still 1.5 against the closed form), and
say so in the commit.

### 2.3 ⚠ The `kind="direct"` / `nulp` gates — check these before regenerating

Four call sites use `assert_regression(kind="direct")`, whose gate is
`np.testing.assert_array_almost_equal_nulp` — the metric that explodes across
a zero crossing (§0.4):

```
tests/sn/regression/test_walk_matvec_baselines.py:45-50   reduction_depth=1   ← the TIGHTEST gate in the tree
tests/sn/operators/test_bc_extraction_matvec.py:357-359   reduction_depth=sn_mesh.nx
tests/sn/operators/test_streaming_operator.py:731,811     reduction_depth=self._PURE_L_NULP
tests/sn/sweep/core/test_affine_carve_baseline.py:208-213
```

`[M]` **None of the four is product-derived** — `_make_cyl` (the walk-baseline
cylinder builder, `tests/sn/operators/test_g_adjoint_reciprocity.py:122`) uses
`Quadrature.level_symmetric(4)`; `test_bc_extraction_matvec.py:206` uses
`level_symmetric`; the others are GL/LS. So the nulp family is currently
immune. **But this is a fixture accident, not a structural guarantee** — the
same file defines `_make_cyl_product` (line 134) on `product(2,4)` explicitly
"the DEGENERATE-class row". If any nulp-gated baseline is ever repointed to the
product builder, a zero-crossing element makes the gate fail by ~4e18 nulp.

---

## §3 — The MoC consumer inventory — and the issue's MoC row is INCOMPLETE

### 3.a/b Bit-exact and tight-tolerance pins: NONE

`[M]` `tests/moc/test_quadrature.py` is the only file that pins `phi`, and only
three of its tests touch it:

| test | assertion | verdict |
|---|---|---|
| `test_azimuthal_angles_in_range` | `phi >= 0`, `phi < pi` | immune |
| `test_azimuthal_angles_uniform_spacing` | `np.allclose(np.diff(phi), pi/n, atol=1e-14)` | immune (1e-14 ≫ 4.4e-16) |
| `test_shapes` | `phi.shape == (16,)` | immune |

There is **no** test anywhere in `tests/moc/` that pins an azimuthal angle
VALUE, a `cos φ` / `sin φ`, or the reflected-track partner index. `[M]` The
whole `tests/moc/` tree is 6 files. **Honest answer: MoC azimuthal angles are
essentially unpinned.** A #325 change there is unobservable by the current
suite — which is a coverage gap, not a safety margin.

### 3.c ⚠ The structural finding: `circle_nodes(2k+1, 4n)` does not fit MoC's data shape

`MOCQuadrature` stores **`phi` (the angle)**, not `(cos φ, sin φ)`. The issue's
consumer table says "repoint `moc/quadrature.py:89` to `circle_nodes(2k+1, 4n)`",
but that returns cosines/sines. Repointing therefore requires a **dataclass
shape change** (`phi` → `cos_phi` / `sin_phi`, or both), and both current `phi`
consumers must move:

| site | what it does |
|---|---|
| `orpheus/moc/geometry.py:316-320` | `cos_phi = np.cos(phi[a])`, `sin_phi = np.sin(phi[a])` — the track-generation direction vector. **Direct drop-in for `circle_nodes`.** |
| `orpheus/moc/geometry.py:222-229` `_reflected_azi_index` | **the site the issue does not list** |

`_reflected_azi_index` is MoC's twin of `Quadrature.reflection_index` — a
nearest-neighbour `argmin` over the **angle**:

```python
def _reflected_azi_index(phi: np.ndarray, azi_index: int) -> int:
    """Both vertical and horizontal reflections map phi -> pi - phi."""
    phi_refl = np.pi - phi[azi_index]
    if phi_refl < 0:
        phi_refl += np.pi
    if phi_refl >= np.pi:
        phi_refl -= np.pi
    return int(np.argmin(np.abs(phi - phi_refl)))
```

`[M]` It currently returns the **correct** answer for every `n_azi` measured —
because the spacing `π/n` dwarfs the ~4.4e-16 residual:

```
MoC _reflected_azi_index vs the exact index-arithmetic partner (n-1-k)
 n_azi  agrees?  search result                            exact n-1-k
     2     True  [1, 0]                                   [1, 0]
     3     True  [2, 1, 0]                                [2, 1, 0]
     4     True  [3, 2, 1, 0]                             [3, 2, 1, 0]
     6     True  [5, 4, 3, 2, 1, 0]                       [5, 4, 3, 2, 1, 0]
     8     True  [7, 6, 5, 4, 3, 2, 1, 0]                 [7, 6, 5, 4, 3, 2, 1, 0]
    16     True  [15, 14, ...]                            [15, 14, ...]
    32     True  [31, 30, ...]                            [31, 30, ...]
    64     True  [63, 62, ...]                            [63, 62, ...]

Is pi - phi[k] BIT-EQUAL to phi[n-1-k] today?
  n_azi=  2: bit-equal 2/2   max|diff| = 0.000e+00
  n_azi=  3: bit-equal 1/3   max|diff| = 4.441e-16
  n_azi=  4: bit-equal 4/4   max|diff| = 0.000e+00
  n_azi=  6: bit-equal 3/6   max|diff| = 2.220e-16
  n_azi=  8: bit-equal 5/8   max|diff| = 2.220e-16
  n_azi= 16: bit-equal 10/16   max|diff| = 4.441e-16
  n_azi= 32: bit-equal 13/32   max|diff| = 4.441e-16
```

The correct #325 treatment of this site is **not** a `circle_nodes` repoint —
it is deleting the search: `φ_k = π(2k+1)/(2n)` gives
`π − φ_k = π(2(n−1−k)+1)/(2n) = φ_{n−1−k}` **by integer arithmetic**, so the
partner is `n_azi − 1 − azi_index`, exactly, with no search and no tolerance.
That is the same "make the mistake unspellable" move the issue makes for
`reflection_index` (§4), and it belongs in the same change. **Add it to the
consumer table.**


---

## §4 — Where `reflection_partners` is BUILT

### The matching code (durable shape)

`reflection_index` (`orpheus/numerics/quadrature/directional.py:373`) is a pure
dict lookup into `Quadrature.reflection_partners`. The dict is built eagerly at
construction by `_compute_sphere_reflection_partners`
(`directional.py:170-187`), which is **SHARED, not per-family** — it is called
identically from all three sphere factories:

```
directional.py:568   Quadrature.lebedev          → _compute_sphere_reflection_partners
directional.py:588   Quadrature.level_symmetric  → _compute_sphere_reflection_partners
directional.py:614   Quadrature.product          → _compute_sphere_reflection_partners
```

(Nexus `callers` on `_find_reflections` confirms exactly these three plus the
composer; `Quadrature.gauss_legendre` takes a different, non-search route —
the `i ↔ N−1−i` GL symmetry.)

The matcher itself is the O(N²) nearest-neighbour the docstring concedes
(`directional.py:125-149`):

```python
def _find_reflections(tx, ty, tz, rx, ry, rz) -> np.ndarray:
    n = len(tx)
    ref = np.empty(n, dtype=int)
    for i in range(n):
        dist = (rx - tx[i]) ** 2 + (ry - ty[i]) ** 2 + (rz - tz[i]) ** 2
        ref[i] = np.argmin(dist)
    return ref
```

`_compute_sphere_reflection_partners` calls it three times, once per axis,
with the axis component negated.

### How inexact is the match TODAY

`[M]` — `exact` counts rows where the partner node is **bit-equal** to the
reflected node:

```
    product 2x4 axis=0: exact=  0/8   max_node_resid=2.999760e-16  max_||mu|-|mu_p||=0.000000e+00
    product 2x4 axis=1: exact=  2/8   max_node_resid=1.999840e-16  max_||mu|-|mu_p||=0.000000e+00
    product 2x4 axis=2: exact=  8/8   max_node_resid=0.000000e+00  max_||mu|-|mu_p||=0.000000e+00
    product 4x8 axis=0: exact=  0/32  max_node_resid=3.455092e-16  max_||mu|-|mu_p||=3.330669e-16
    product 4x8 axis=1: exact=  4/32  max_node_resid=2.303395e-16  max_||mu|-|mu_p||=2.220446e-16
    product 4x8 axis=2: exact= 32/32  max_node_resid=0.000000e+00  max_||mu|-|mu_p||=0.000000e+00
   product 8x16 axis=0: exact=  0/128 max_node_resid=7.216450e-16  max_||mu|-|mu_p||=3.330669e-16
   product 8x16 axis=1: exact=  8/128 max_node_resid=6.106227e-16  max_||mu|-|mu_p||=6.106227e-16
   product 8x16 axis=2: exact=128/128 max_node_resid=0.000000e+00  max_||mu|-|mu_p||=0.000000e+00
   product 4x12 axis=0: exact=  0/48  max_node_resid=1.054712e-15  max_||mu|-|mu_p||=1.054712e-15
   product 4x12 axis=1: exact=  4/48  max_node_resid=7.216450e-16  max_||mu|-|mu_p||=5.551115e-16
   product 4x12 axis=2: exact= 48/48  max_node_resid=0.000000e+00  max_||mu|-|mu_p||=0.000000e+00
            LS4 axis=0: exact= 24/24  max_node_resid=0.000000e+00  max_||mu|-|mu_p||=0.000000e+00
            LS4 axis=1: exact= 24/24  max_node_resid=0.000000e+00  max_||mu|-|mu_p||=0.000000e+00
            LS4 axis=2: exact= 24/24  max_node_resid=0.000000e+00  max_||mu|-|mu_p||=0.000000e+00
          LEB17 axis=0: exact=110/110 max_node_resid=0.000000e+00  max_||mu|-|mu_p||=0.000000e+00
          LEB17 axis=1: exact=110/110 max_node_resid=0.000000e+00  max_||mu|-|mu_p||=0.000000e+00
          LEB17 axis=2: exact=110/110 max_node_resid=0.000000e+00  max_||mu|-|mu_p||=0.000000e+00
            GL8 axis=0: exact=  8/8   max_node_resid=0.000000e+00  max_||mu|-|mu_p||=0.000000e+00
            GL8 axis=1: exact=  8/8   max_node_resid=0.000000e+00  max_||mu|-|mu_p||=0.000000e+00
            GL8 axis=2: exact=  8/8   max_node_resid=0.000000e+00  max_||mu|-|mu_p||=0.000000e+00
```

Three things this says that the issue does not:

1. **The issue's claim is confirmed and is worse on the x-axis than stated.**
   For `product`, `reflection_index(0)` (`"x"`) has **zero** exact matches at
   every size, and `reflection_index(1)` (`"y"`) has only the `sin(0)=0` rows.
   Axis 2 (`"z"`) is already exact — `μ_z` comes from `leggauss`, which is
   bit-symmetric.
2. **The partner's `|Ω·n̂|` genuinely differs** — up to `1.05e-15` at
   `product(4,12)`. This is the ERR-042 measure-preservation invariant being
   satisfied only to round-off, exactly as the issue says.
3. **Every other family is already exact.** LS4, Lebedev-17 and GL-8 are
   `exact = N/N` on all three axes with `max_resid = 0.0`. Only `product`
   needs this fix — which supports "one primitive, targeted at the violators".

### Production consumers of `reflection_index`

`[M]` (grep + Nexus; `reflection_index` is reached through several
BC-realization layers that `callers` alone under-reports):

| file:line | context |
|---|---|
| `orpheus/geometry/boundary/_specular.py:145` | `certify_measure_preservation` (ERR-042) |
| `orpheus/geometry/boundary/_specular.py:184` | `certify_involution` (ERR-044) |
| `orpheus/geometry/boundary/_specular.py:220` | `certify_inflow_outflow` (ERR-045) |
| `orpheus/sn/boundary/realizer.py:383` | builds the realized `PermutationOperator` table |
| `orpheus/sn/loss_representation/__init__.py:3294` | pole-face seed gather |
| `orpheus/sn/loss_representation/__init__.py:3470` | `mirror` for the curvilinear walk |
| `orpheus/sn/loss_representation/__init__.py:4189` | matvec `mirror` |
| `orpheus/sn/loss_representation/__init__.py:4593` | transpose-matvec `mirror` |
| `orpheus/derivations/discrete/sn/sweep_acyclicity.py:287-399` | the acyclicity proof takes it as a parameter |

Plus doc-level references in `geometry/boundary/{__init__,reflective,albedo,_factors}.py`
and `numerics/operator.py:2190` — those are prose, but `numerics/operator.py:2190`
carries a **dangling** `:func:`~orpheus.sn.quadrature.Quadrature.reflection_index``
target (see §5.3).

### What #325 buys here, and the one thing to check

`[M]` Post-#325 the reflected node IS bit-equal to a real node, so
`np.argmin(dist)` finds a `dist == 0.0` entry. **Check before assuming it is
unique:** `argmin` returns the FIRST minimum, so the exact-lookup claim needs
"no two nodes coincide" (true for all shipped rules) — this is worth an assert,
not an assumption. The issue's acceptance criterion ("replace the search or add
a gate that the match is exact") is the right shape; prefer the **gate**
(`dist.min() == 0.0` for every row) over silently keeping the search.

---

## §5 — The doc blast radius

Per `.claude/rules/coding-standards.md`: Python-domain xrefs are **not** gated
by `-W` unless the build runs `-n`, so the searches below are the only thing
that catches this class.

### 5.1 ⛔ CONFIRMED: `angular_trace_space.py` module docstring over-claims

The brief's suspicion is **correct and measured**.
`orpheus/numerics/spaces/angular_trace_space.py:79-85` says:

> Empirically (``eps_probe.py``, Gauss-Legendre ``N=2..64`` + Lebedev orders
> ``3..53``): nominally-tangential cosines are **exactly** ``0.0``

`[M]` That is TRUE for the two families it cites and **FALSE for `product`**:

```
M1 — does `product` today put cosines in the (0, TANGENTIAL_EPS] band?
n_mu n_phi  axis n_in_band n_exact_0   max_in_band
   2     4  mu_x         4         0  1.499880e-16
   2     4  mu_y         2         2  9.999199e-17
   4     8  mu_x         8         0  1.727546e-16
   4     8  mu_y         4         4  1.151697e-16
   6     8  mu_x        12         0  1.783906e-16
   6     8  mu_y         6         6  1.189271e-16
   8    16  mu_x        16         0  1.805800e-16
   8    16  mu_y         8         8  1.203867e-16
   4    12  mu_x         8         0  1.727546e-16
   4    12  mu_y         4         4  1.151697e-16

M2 — the SAME probe on the families the docstring cites
   GL8  mu_x in_band=  0 exact0=  0
   GL8  mu_y in_band=  0 exact0=  8
   GL8  mu_z in_band=  0 exact0=  8
   GL7  mu_x in_band=  0 exact0=  1
 LEB17  mu_x in_band=  0 exact0= 12
 LEB17  mu_y in_band=  0 exact0= 12
 LEB17  mu_z in_band=  0 exact0= 12
   LS4  mu_x in_band=  0 exact0=  0
```

`product` puts **`n_mu` to `2·n_mu` cosines per axis strictly inside
`(0, TANGENTIAL_EPS]`** — the exact band the docstring says is empty.

**Worse: the sentence claims the property is what makes `TANGENTIAL_EPS`
principled, and the guard test that "guards the gap" cannot see the
violation.** `tests/numerics/test_angular_trace_space.py:333-337` builds its
quadrature list from `gauss_legendre` and `lebedev` ONLY:

```python
quads = [Quadrature.gauss_legendre(n) for n in (2, 3, 4, 5, 7, 8, 16, 32, 64)]
for order in (3, 5, 7, 11, 17, 29, 53):
    quads.append(Quadrature.lebedev(order))
```

and its assertion is exactly the property `product` violates:

```python
assert not np.any((amu > 0.0) & (amu <= TANGENTIAL_EPS)), (
    f"{type(q).__name__} {axis} has a cosine in the round-off "
    f"band (0, eps] — the exactly-zero assumption is violated")
```

**So the guard is a vv Mode-7 vacuity: it asserts a property over a family set
chosen (accidentally) to exclude the only member that breaks it.** Adding
`product` to that list on the CURRENT tree makes it RED. That is the perfect
pre-#325 `xfail(strict=True)` — it becomes green exactly when #325 lands, and
it *is* the issue's acceptance criterion "a gate that every production
quadrature has EXACTLY zero tangential cosines". `[M]` post-#325 measurement:

```
M7 — post-fix: does `product` then have EXACT tangential zeros?
  n_mu=2 n_phi=4: mu_x in_band=0 exact0=4 | mu_y in_band=0 exact0=4
  n_mu=4 n_phi=8: mu_x in_band=0 exact0=8 | mu_y in_band=0 exact0=8
  n_mu=8 n_phi=16: mu_x in_band=0 exact0=16 | mu_y in_band=0 exact0=16
  n_mu=4 n_phi=12: mu_x in_band=0 exact0=8 | mu_y in_band=0 exact0=8
  n_mu=2 n_phi=24: mu_x in_band=0 exact0=4 | mu_y in_band=0 exact0=4
```

Also note the docstring cites `eps_probe.py` as its evidence and **that file
does not exist anywhere in the tree** (`[M]` `find . -name 'eps_probe*'` →
nothing; the only two mentions are the docstring itself and the test's prose).
A citation to a deleted probe script is an unfalsifiable provenance claim — fix
it in the same pass.

### 5.2 Doc prose asserting the equispaced-φ convention

| file:line | claim | #325 effect |
|---|---|---|
| `docs/theory/methods/sn/curvilinear_one_group.rst:3896-3897` | names ``np.linspace(0, 2π, n_φ, endpoint=False)`` verbatim | must be repointed to the new primitive |
| `docs/theory/methods/sn/curvilinear_one_group.rst:3903-3906` | "for **even** `n_φ` the grid hits `φ = π` exactly … `μ_y = sinθ sin π = 0`" | `[M]` **`φ[n/2] == np.pi` is True and `cos == -1.0` is True today, but `sin(π) = 1.2246467991473532e-16`, so `μ_y = 0` is exact-arithmetic true and FLOAT-FALSE today.** #325 makes it float-true (`-0.0`). The τ_raw,0 = 0 / no-seed argument is unaffected — it rides on `cos`, which is already exact. |
| `docs/theory/methods/sn/curvilinear_one_group.rst:3942` | the "equispaced (trapezoidal, spectral)" table row | unchanged (the NODES are the same points; only their float representation improves) |
| `docs/theory/methods/sn/angular_quadrature.rst:51-70` | the `product` section; "Within each level, ordinates are sorted by increasing `η = sinθ cos φ` to match the α recursion convention from `BaileyMorelChang2010`" | **the §0.3 question lives here.** The doc says the α-recursion wants **η-ordering**, and this is the only place the convention is stated. It does NOT say what to do when η ties — which is precisely what #325 creates. This page must gain the tie-break rule. |
| `docs/theory/foundations/discrete_measures.rst:352-357` | `μ_S2 = gauss_legendre(n_mu) * equispaced(0, 2π, n_phi)` — names the `orpheus.numerics.measure.equispaced` primitive | check whether `equispaced` itself is a second trig site (it is a *measure* factory, not a cosine factory — verify before repointing) |
| `orpheus/numerics/quadrature/rules_product.py:6-10, 113-114` | "the bit-identical contract enforced by the regression snapshots requires this exact convention" | **must be rewritten** — this sentence is the reason a future session would refuse to change the line #325 changes. |

### 5.3 Measured doc claims that #325 makes stale

`docs/theory/foundations/boundary_conditions.rst` carries a **measured table**
about `product`'s classifier disagreement that #325 zeroes out. `[M]` the
before/after numbers a doc update needs:

```
product(2, 4)   N = 8
   face | OLD tang OLD misadmit | NEW tang NEW misadmit
   xmin |        4            2 |        4            0
   xmax |        4            2 |        4            0
   ymin |        4            0 |        4            0
   ymax |        4            2 |        4            0
   zmin |        0            0 |        0            0
   zmax |        0            0 |        0            0
  xmax tangential Omega.n  OLD: [-1.4998798865218462e-16, -1.4998798865218462e-16,
                                  4.9995996217394874e-17,  4.9995996217394874e-17]
  xmax tangential Omega.n  NEW: [-0.0, -0.0, 0.0, 0.0]

product(4, 8)   N = 32
   xmin |        8            4 |        8            0
   xmax |        8            4 |        8            0
   ymin |        8            0 |        8            0
   ymax |        8            4 |        8            0
```

**The tangential COUNT is unchanged** (4 of 8 stays 4 of 8), so these claims
survive:

* `boundary_conditions.rst:145` — "a cylinder under `product(n_mu=2, n_phi=4)`
  carries **4 of 8** ordinates tangential at `xmax`" ✅ still true
* `boundary_conditions.rst:4389` — same ✅ still true
* the three-way partition / "never spell 'not inflow' as 'outflow'" doctrine ✅
  still true and still load-bearing

**These become stale (present-tense-false):**

* `boundary_conditions.rst:1000-1032` — the classifier-divergence table. The
  `product(2,4)` row's "rows the `> 0.0` test claims that Γ₊ does not" goes
  `2 → 0`; the whole "face-asymmetric within a single quadrature … `ymin`
  carries the same tangential count with zero mis-admissions, because the sign
  flip moves the round-off across zero" narrative dissolves.
* `boundary_conditions.rst:1123-1137` — "the four tangential ordinates carry
  `Ω·n̂` of `+5.0e-17` (twice) and `-1.5e-16` (twice) — round-off, not exact
  zero". After #325 they are exact `±0.0`.
* `boundary_conditions.rst:1080` — the "rows the `> 0.0` test mis-admits"
  column header.

Per the user's articulation standard (present-tense-false is the bug; past-tense
history stays), these are a **re-framing** job, not a deletion: the B3.4a
narrative is genuine history and the numbers are its evidence. The fix is to
tense them ("**before #325** the four tangential ordinates carried …") and add
the post-#325 row. This is an `archivist` hand-off, and it is not optional —
`boundary_conditions.rst` is the page that teaches this bug class.

### 5.4 Dangling Python-domain xrefs in the files #325 touches

`[M]` `orpheus/sn/quadrature` **does not exist** (`ls orpheus/sn/` → no
`quadrature.py`, no `quadrature/`). These `:class:` / `:func:` roles therefore
resolve to nothing and render as plain text with **no `-W` warning**:

```
orpheus/numerics/quadrature/rules_product.py:8    :class:`orpheus.sn.quadrature.ProductQuadrature`
orpheus/numerics/quadrature/rules_product.py:63   :class:`orpheus.sn.quadrature.ProductQuadrature`
orpheus/numerics/quadrature/rules_product.py:102  :class:`orpheus.sn.quadrature.ProductQuadrature`
orpheus/numerics/measure.py:1102                  :class:`orpheus.sn.quadrature.ProductQuadrature`
orpheus/geometry/reduced_operator.py:742          :class:`~orpheus.sn.quadrature.ProductQuadrature`
orpheus/numerics/operator.py:2190                 :func:`~orpheus.sn.quadrature.Quadrature.reflection_index`
tests/numerics/test_rules_product.py:6, 61        :class:`~orpheus.sn.quadrature.ProductQuadrature[.create]`
```

Three of the seven are in `rules_product.py` — the file #325 rewrites — so they
are free to fix in the same commit. (`docs/theory/foundations/spherical_harmonics.rst:579`
and `docs/.../discrete_measures.rst:734` mention `orpheus.sn.quadrature` as
literal text, not roles; those are prose-stale but harmless.)


---

## §6 — ⭐ The finding that reframes §0.2: the tie hazard ALREADY EXISTS on `level_symmetric`

`rules_sphere.py:211-215` sorts level-symmetric levels with the **same**
`np.argsort(mu_x[idx])` pattern. `level_symmetric` is the family that is
already algebraic and already exact — so it already has **massive exact ties**.
`[M]`:

```
Does the LEVEL_SYMMETRIC family ALREADY have intra-level eta ties today?
(rules_sphere.py:214 uses the same np.argsort(mu_x[idx]) pattern)
            family   ties  levels   argsort kinds agree?
level_symmetric(4)     18       2                   True
level_symmetric(6)     36       3                  False
level_symmetric(8)     60       4                  False
level_symmetric(16)   216       8                  False
```

**Consequences:**

1. #325 does **not introduce** the tie / non-determinism hazard. It exists
   today on `level_symmetric` — the family the cylinder DD snapshot, the 2-D
   octant snapshots, and the walk baselines all use — and `LS6+` already shows
   sort-kind disagreement. #325 merely extends the hazard to `product`.
2. That makes the fix cheaper to justify and **independent of #325's merits**:
   pinning `kind="stable"` (plus a physically meaningful secondary key) closes a
   latent reproducibility hole that is live right now. A future numpy introsort
   change would silently permute LS6+ level order and move every derived value
   at the ULP scale — invisible under the correctness gates, RED under
   `-W error::DriftWarning`.
3. It also means the "exact algebraic generation ⟹ ties" relationship is not a
   #325 side-effect but the *expected* consequence of exact symmetry. The
   already-exact families demonstrate it. Framing #325 as "it creates ties" is
   backwards: exact symmetry creates ties, and ties are the *evidence* the
   symmetry is exact.

The secondary key the physics wants is stated at
`docs/theory/methods/sn/angular_quadrature.rst:64-67` only as far as "sorted by
increasing η to match the α recursion convention (Bailey)". It is silent on
ties. **Deciding the tie-break is a physics question that must be answered from
the Bailey Eq. 50 convention before the sort is pinned** — my read of §0.3 is
that η-monotone is what the recursion needs and the ±ξ interleave is free, but
that should be confirmed against the reference, not assumed. This is a
`literature-researcher` / `archivist` item, not a coding one.

Also noted in passing at `rules_sphere.py:212`: `tol = 1e-12` decides μ-level
membership by float comparison — a third magic tolerance in the same
neighbourhood, exactly-zeroable by the same argument (`|μ_z|` values come from a
table, so equality is exact). Out of #325's scope but the same family of defect.

---

## §7 — What the change must include (beyond the issue's acceptance criteria)

Additions the measurements above force. Everything already in the issue's
acceptance list still stands.

| # | item | why | owner |
|---|---|---|---|
| A1 | **Pin `np.argsort(..., kind="stable")`** in BOTH `rules_product.py:139` and `rules_sphere.py:214` | §0.2 + §6: tie resolution is implementation-defined and already non-deterministic on LS6+ | code |
| A2 | **Decide and DOCUMENT the tie-break** the α-recursion needs (η-monotone + what within a tie?) against Bailey Eq. 50, in `angular_quadrature.rst` | §0.3: `η` is 2-to-1 over `φ ∈ [0,2π)`; the current order interleaves the two half-circles arbitrarily | archivist / literature |
| A3 | **MoC: replace `_reflected_azi_index` with `n_azi − 1 − k`** (integer arithmetic, no search) | §3.c: the issue's MoC row is incomplete; `circle_nodes` alone does not reach this site | code |
| A4 | **MoC: decide the `MOCQuadrature` data shape** (`phi` vs `cos_phi`/`sin_phi`) before repointing | §3.c: `circle_nodes` returns cosines; MoC stores the angle | code |
| A5 | **Add `product` to `test_eps_sits_in_the_round_off_to_genuine_gap`'s quadrature list** — as `xfail(strict=True)` BEFORE the fix, flipping green when it lands | §5.1: the guard is currently vacuous over the only violating family (vv Mode 7). This is the issue's "gate that every production quadrature has EXACTLY zero tangential cosines" | test-architect |
| A6 | **Migrate `test_domain_is_the_tangential_band_gamma_plus`** | §1.a: its ACTIVATION guard fires because #325 removes the discriminator. Re-anchor on a synthetic quadrature carrying a small-but-nonzero cosine (the classifier's contract is about the CODE, not about which shipped rule happens to trip it) — or retire it explicitly | test-architect |
| A7 | **Fix the `angular_trace_space.py:73-91` docstring** — the "exactly `0.0`" claim, the family list it was measured over, and the citation to the non-existent `eps_probe.py` | §5.1 | archivist |
| A8 | **Re-frame `boundary_conditions.rst` §1000-1137** to past tense + add the post-#325 numbers | §5.3 | archivist |
| A9 | **Repoint `curvilinear_one_group.rst:3896-3906`** (the `np.linspace` spelling; the `μ_y = 0` claim becomes float-true) | §5.2 | archivist |
| A10 | **Rewrite `rules_product.py:6-10, 113-114`** — the "bit-identical contract enforced by the regression snapshots requires this exact convention" sentence is the standing objection to the change itself | §5.2 | code |
| A11 | **Fix the 7 dangling `orpheus.sn.quadrature.*` Python-domain xrefs** (3 of them in `rules_product.py`) | §5.4 — silent under `-W` without `-n` | code |
| A12 | **Regenerate `cyl_1g_homogeneous_product_dd_n20.npz`** with the justification recorded, per the generator's own protocol | §2.1/§2.2 — the ONE snapshot that moves | code |
| A13 | **Guard the exact lookup** in `_find_reflections`: assert `dist.min() == 0.0` per row rather than silently keeping the argmin search | §4 | code |
| A14 | **Note the signed-zero** the primitive produces (`-0.0` on half the axis nodes, `[M]` `xmax` → `[-0.0, -0.0, 0.0, 0.0]`) | any consumer doing `np.sign`, `1/x`, or `copysign` on a direction cosine sees `-0.0`. `==` / `<` / `>` comparisons are unaffected | code review |

### Not needed

* **The degenerate-ordinate branch is safe.** `[M]`
  `orpheus/sn/mesh/augmented_mesh.py:1079` `_DEGENERATE_ABS_ETA_THRESHOLD = 1e-15`
  and the predicate is `abs(eta) < threshold` (`:1270`, `:1473`) — `6e-17` and
  `0.0` both satisfy it, so the pure-azimuthal branch fires identically before
  and after. `_octant_sign_predicate`'s `_OCTANT_SIGN_EPS = 1e-15`
  (`directional.py:95`) likewise. No branch flips.
* **No `select_quadrature` production caller** — the registry route to
  `ProductQuadrature` is test-only, so the registry needs no migration beyond
  the factory swap.

---

## §8 — Measured summary of the empirical test sweep

Method: a pytest plugin (`/Users/rodrigo/.claude/jobs/c30e4f25/tmp/plugin_325.py`)
swapped `product_mu_phi` for the algebraic implementation from the issue body,
patching the definition module, `directional`, `registry`, both `__init__`
re-exports, and the registry spec's `factory`. Canonical invocation
`python -O -m pytest ... -m "not slow"`, serial.

| batch | scope | result |
|---|---|---|
| 1 | the 50 test files that call `.product(` | `[M]` `1 failed, 1283 passed, 13 deselected, 61 xfailed in 265.17s` |
| 2 | `tests/sn/regression tests/moc tests/cross_method tests/numerics tests/geometry` (the INDIRECT consumers grep misses — `test_dd_regression.py` reaches `product` only through `_generate_snapshots.CASES`) | `[M]` `1741 passed, 4 skipped, 11 deselected, 1 xfailed, 7 warnings in 252.59s` |

The only failure, batch 1:

```
FAILED tests/sn/operators/test_angular_average_operator.py::TestOutflowHalfTraceIsTheDomain::test_domain_is_the_tangential_band_gamma_plus
  - AssertionError: fixture no longer discriminates: Γ₊ [0, 4] and the retired '> 0.0' set [0, 4] agree
```

The only value movement, batch 2 (both PASS the correctness gate; both raise
the informational tripwire):

```
tests/sn/regression/test_dd_regression.py:101: DriftWarning: cyl_1g_homogeneous_product_dd_n20: k_eff drifted 1 ULP / 1.48e-16 rel (within tol 1.0e-11)
tests/sn/regression/test_dd_regression.py:110: DriftWarning: cyl_1g_homogeneous_product_dd_n20: scalar_flux drifted 2 ULP / 2.62e-16 rel (within tol 1.0e-09)
```


---

## §9 — ⭐ Why the order effect is FP-only: a structural proof, not just a measurement

§A measures the order effect at round-off scale. Here is *why*, which is the
part that makes the conclusion durable rather than case-dependent.

The cylindrical α-recursion (`docs/theory/methods/sn/curvilinear_one_group.rst:229-241`,
`:eq:`alpha-cylindrical``, Bailey–Morel–Chang 2010) is

  α_{p,m+½} = α_{p,m−½} − w_m η_m,   α_{p,½} = 0

with the ordinates on level `p` sorted by increasing `η`. Consider two
adjacent ordinates `a`, `b` with **`η_a == η_b` exactly** (a #325 tie). Swapping
them leaves:

* every α **before** the pair unchanged (same prefix),
* every α **after** the pair unchanged (the pair contributes
  `−w_a η_a − w_b η_b`, symmetric in the swap),
* the **intermediate** α between them changing by `w_a η_a − w_b η_b
  = η(w_a − w_b)`.

So the tie-break is a genuine mathematical choice **iff the two tied ordinates
carry different weights.** `[M]` They never do — the weights are constant within
a μ-level for every production family:

```
          product(2,4): weights constant per level = True  max weight spread WITHIN a tie = 0.000e+00
          product(4,8): weights constant per level = True  max weight spread WITHIN a tie = 0.000e+00
         product(4,12): weights constant per level = True  max weight spread WITHIN a tie = 0.000e+00
    level_symmetric(4): weights constant per level = True  max weight spread WITHIN a tie = 0.000e+00
    level_symmetric(6): weights constant per level = True  max weight spread WITHIN a tie = 0.000e+00
    level_symmetric(8): weights constant per level = True  max weight spread WITHIN a tie = 0.000e+00
```

(For `product` this is structural: `w_{p,m} = w_GL(μ_p)·2π/n_φ` is independent
of `m` by construction — `rules_product.py:116,134`.)

**Conclusion: the intra-level tie-break is provably α-invariant in exact
arithmetic on every shipped quadrature. The ≤50-ULP movement §A measures is
pure FP non-associativity, nothing more.** #325 is a re-baseline, not a change
of recurrence.

**This also gives A1/A2 their honest justification.** Pinning `kind="stable"` is
not needed for *correctness* — it is needed for *reproducibility*: without it,
the bit-values of every cylinder snapshot are a function of numpy's introsort
implementation. And the tie-break rule A2 should be documented as **"any; the
recursion is invariant because weights are constant within a level"**, with the
invariance stated as the reason — which is a much stronger doc sentence than
picking an arbitrary secondary key. If a future quadrature ever carries
level-varying weights, that invariance breaks, and the doc sentence is exactly
what makes the breakage visible.


### §8.1 Batch 3 (broadest sweep) — COMPLETE

A third sweep over `tests/sn tests/derivations tests/transport` (the MMS /
verification / derivation surface, which reaches `product` through
`orpheus/derivations/continuous/mms/sn.py` and
`orpheus/derivations/discrete/sn/contamination.py`, and which covers the large
majority of `tests/sn` that batch 1's 50-file list did NOT).

`[M]` The run emitted its **complete per-test outcome stream to `[100%]`**, then
the harness killed the process during the warnings-summary teardown — so the
`N passed` summary line is absent, but every per-test marker was written:

```
progress stream terminal state : [100%]
F markers                      : 1
E markers                      : 0
```

**Identification of the single `F`, by elimination `[M]`:** the known
fixture-discrimination test is inside batch 3's collection (`tests/sn` is in
scope), and it fails under the patch when run alone —

```
$ PYTHONPATH=.../tmp .venv/bin/python -O -m pytest \
    "tests/sn/operators/test_angular_average_operator.py::TestOutflowHalfTraceIsTheDomain::test_domain_is_the_tangential_band_gamma_plus" \
    -p plugin_325 -q --no-header
FAILED tests/sn/operators/test_angular_average_operator.py::TestOutflowHalfTraceIsTheDomain::test_domain_is_the_tangential_band_gamma_plus
  - AssertionError: fixture no longer discriminates: Γ₊ [0, 4] and the retired ...
1 failed, 1 warning in 0.09s
```

— so with batch 3's total `F` count equal to 1, that test **is** the one
failure. **No other test in `tests/sn`, `tests/derivations` or
`tests/transport` fails under the swap, and there are no errors.**

Combined across all three batches: the entire non-slow SN / geometry / numerics
/ MoC / cross-method / derivations / transport surface breaks in **exactly one
place**, and that place is a guard whose fixture depends on the defect.

To reproduce the summary line (the run is cheap to redo, ~40 min serial):

```
PYTHONPATH=/Users/rodrigo/.claude/jobs/c30e4f25/tmp \
  .venv/bin/python -O -m pytest tests/sn tests/derivations tests/transport \
  -p plugin_325 -q --no-header -m "not slow" -rf
```

(The plugin is at `/Users/rodrigo/.claude/jobs/c30e4f25/tmp/plugin_325.py`; it
patches live objects only and touches no tracked file.)

---

## §10 — What I did NOT cover

Stated so the gaps are not mistaken for clean results.

1. **`-m slow` tests were deselected** in all three batches (11–13 per batch).
   The heavy curvilinear L1 rows live there and they are the ones most likely to
   sit near a convergence-tolerance edge. Run them before merging. This is the
   single largest remaining gap in the empirical evidence.
2. **`-W error::DriftWarning` was not run.** Under that flag
   `cyl_1g_homogeneous_product_dd_n20` WILL fail (§2.2) — by design.
3. **Batch 3's numeric summary line is absent** — the process was killed during
   teardown after emitting the complete `[100%]` outcome stream. The failure
   count (1) and error count (0) are read off the per-test markers and the
   identity of the failure is established by elimination (§8.1), not from a
   summary line. Re-run if you want the count in the record.
4. **The MoC half was not swap-tested.** Repointing MoC requires the dataclass
   shape decision (§3.c), so there was nothing to swap in. Given §3.a/b, the
   current suite would not detect the change anyway.
5. **`derivations/discrete/sn/balance.py:367`** (the algebra of record) was read
   only for its call shape, not swap-tested. Its symbolic tests
   (`tests/derivations/`) are inside batch 3.
6. **Cross-platform / cross-numpy tie resolution** was measured on numpy 2.4.4,
   macOS arm64 only. The §0.2 claim is that the guarantee is absent, not that a
   specific other platform differs.
7. **`orpheus/numerics/measure.py::equispaced`** (a midpoint-rule measure
   factory with the same convention lineage) was inspected but not audited as a
   fourth #325 consumer — it produces angles, not cosines, so it is not a trig
   site, but the design of record should say whether it is in or out of scope.
