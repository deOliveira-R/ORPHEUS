---
name: angular-layer-hidden-transformations
description: Where the ANGULAR/quadrature/moment layer performs rotations, reflections and coordinate changes WITHOUT naming them as group elements — the durable map (SH frame rotation, octant=D_2h stratification, twin partner engines, the six open #325 sites).
metadata:
  type: project
---

Audited 2026-08-03 on `refactor/operator-strategy-layers` @ `bfedc621`; evidence
(with `file:line` + measurements) in `scratch/geom_transforms_audit_angular.md` —
that file is the transient artifact, this note is the part that survives churn.
Re-derive any line number via Nexus. Companion to
[[quadrature-landscape-durable]].

**Why:** the campaign question is "which hand-rolled transformations could share
the `symmetry.py` / `roots_of_unity.py` machinery?" These are the answers that a
grep cannot re-find, because most of them are not spelled as transformations.

## The four structural facts

1. **The SH basis works in a DIFFERENT frame from everything else, and applies
   the change by variable choice.** `numerics/basis/spherical_harmonic_basis.py`
   declares the polar axis to be `μ_x` with azimuth in the `(μ_y, μ_z)` plane.
   Measured: the implied `(x,y,z)_std = (μ_y, μ_z, μ_x)` is the **120° rotation
   about the `(1,1,1)` body diagonal** — a proper `O_h` element, present in
   `_octahedral_ops()`, and **NOT expressible as a `SubgroupOfO3` tag** (`Cn(3)`
   is realized about **z** and measurably does not contain it). Every other
   angular consumer (`_cyclic_ops`, `_vertical_mirrors`, `Dnh(n_φ)`,
   `LevelStructure.azimuth`, `hemisphere = sign(μ_z)`) measures azimuth about
   **z**. No crosswalk exists. Because the rotation is applied by *which array
   is called `cos θ`*, there is no matrix for any check to test.

2. **The "octant partition" is the orbit-type stratification of `(Z₂)³ = D_2h`,
   and "octant" is a misnomer.** Measured part counts: `lebedev(17)` → **26**
   (8 chambers + 18 walls), `product(4,8)` → 16, `level_symmetric(8)` → 8,
   `gauss_legendre(8)` → 2. Labels with a `0` component are exactly `Fix(σ_a)`
   — the set `OrbitCertificate.singular_set` already decides EXACTLY (integer
   identity `π_M(i)==i`). The group is already expressible and already named
   elsewhere: `registry.GEOMETRY_ANGULAR_SYMMETRY` assigns `Dnh(2)` to
   cylinder/cartesian2d, and `Dnh(2)` realizes exactly the 8 diagonal sign
   matrices (measured `|G|=8`). Three views of one object, mutually unaware.

3. **`_orbit_closure` is NOT the only permutation-computing engine.** Five
   others: (a) `symmetry._is_reflection_invariant_1d` — `argsort`+`searchsorted`,
   window `atol*10` (vs `_orbit_closure`'s `argmin`, window `atol*100`), and it
   **discards** the permutation it builds; (b) exact FORMULA endpoints that never
   route through it (`arange(N)[::-1]` for GL; `j = n_half−1−p−k` for LS_N;
   `[::-1]` as the Reynolds/group-average in `GeneratingMeasure.gauss`);
   (c) unguarded `argmin` searches (`moc/geometry._reflected_azi_index`,
   `MOCMesh._find_link`, `radial_characteristic_field`'s most-inward ordinate,
   `pole_angular_closure._edge_seed_stencil`). The tree already owns the typed
   endpoint the angular layer doesn't use: `numerics.operator.PermutationOperator`.

4. **The curvilinear fiber is a CIRCLE that the code only ever treats as an
   ordered interval.** The α-recursion (`geometry/reduced_operator.py`) marches
   `argsort(η)` order; the M-M τ builds edges `[−sinθ, midpoints…, +sinθ]` —
   the circle cut at `φ=π` with the two endpoints never identified. No rotation,
   no reflection, no wrap is spelled. The one reflection that IS there is the
   pole continuation `ψ(0,+μ)=ψ(0,−μ)`, correctly routed through
   `quad.reflection_index("x")`. `LevelStructure.fiber` holds the CORRECT
   `(hemisphere, azimuth)` ordering and the α-recursion does **not** consume it
   (issue #326).

## Two live defect leads (measured, not verdicts)

- **`Quadrature.spherical_harmonics`'s slab claim is FALSE at `ℓ≥2`.** It
  promises "only `m=0` is non-zero" for GL1D; measured, the `m>0` slots carry
  ~0.83 at `ℓ=2,3`. Mechanism: for a slab `sinθ≈0.9` so the `on_axis` guard never
  fires, `(cos φ, sin φ)` becomes `(0,0)` (not on `S¹`), `arctan2(0,0)=0`, so
  `cos(mφ)≡1`. `ℓ≤1` is clean only because it is hardcoded. Reconstruction with
  the full table differs from the `m=0`-only route by **4.4×**. Unreachable by
  any test: the only `scattering_order≥2` test uses `lebedev(17)`.
- **The moment layer destroys the exactness the quadrature just bought.**
  `spherical_harmonic_basis` collapses `(cos φ, sin φ)` through `arctan2` and
  re-evaluates `cos(mφ)`. At `m=1` that is a pure identity round trip. Measured
  on `product(4,8)`: `arctan2` route keeps **0** exact zeros, the de-Moivre route
  (`T_m(cos φ)`, `sin φ·U_{m−1}(cos φ)`) keeps **8**.

## Epsilon family (all idle, one justified by a ghost)

`_OCTANT_SIGN_EPS = 1e-15` (`quadrature/directional.py`), `TANGENTIAL_EPS =
8.88e-16` (`spaces/angular_trace_space.py`), and a bare `1e-14` in
`pole_angular_closure._edge_seed_stencil` — three different numbers for one
question. Measured min nonzero `|cos|` over the shipped rules: `1.57e-1`. The
first is defended by "keep in lockstep with `_DEGENERATE_ABS_MU_THRESHOLD`",
which **exists nowhere** in `orpheus/` or `tests/`.

## Where hand-built 3×3 matrices live

`numerics/symmetry.py` (7 constructors: `_rotation_z`, `_reflections`,
`_inversion_op`, `_vertical_mirrors`, `_octahedral_ops`, `_icosahedral_ops`,
`_rotation_about_axis`) and ONE re-spelling in
`quadrature/directional._compute_sphere_reflection_partners`
(`eye(3)` with a `−1`, i.e. a local copy of `_reflections`). Nowhere else in
production — `basis/`, `transport/`, `sn/`, `moc/` apply transformations only as
index permutations or sign flips. There is **no Wigner-D anywhere** (`grep
wigner` returns only Wigner–Seitz pin cells); the only group action on a MOMENT
vector in the tree is `fold_moments_to_radial_characteristic`'s `(±1)^ℓ`.
`_inversion_op` is a group element with no tag (no `C_i`/`S_n` entry).
