# Buried geometric transformations — the ANGULAR / QUADRATURE / MOMENT layer

Audit on `refactor/operator-strategy-layers` @ `bfedc621`. Tree clean at open and
at close for every file cited (only `.claude/skills/vv-principles/*` and untracked
`scratch/` differ). Read-only: nothing in `orpheus/` or `tests/` was modified.

Scope: `numerics/quadrature/`, `numerics/symmetry.py`, `numerics/roots_of_unity.py`,
`numerics/moment_layout.py`, `numerics/basis/`, `transport/fields/`,
`transport/operators/scattering.py`, and the ANGULAR parts of `sn/sweep/`. Two
adjacent files carry angular transformations and are included, flagged:
`geometry/reduced_operator.py` (the α-recursion) and
`numerics/spaces/angular_trace_space.py` (the Ω·n̂ table). MoC is listed only
where it is a #325 site.

Spelling classes: **(M)** explicit matrix · **(P)** index permutation ·
**(S)** sign flip / coordinate negation · **(A)** angle arithmetic ·
**(O)** other.

---

## The table

### A. `numerics/symmetry.py` — the machinery itself (the only place elements are matrices)

| file:line | the transformation, in math | spelled TODAY | group element? | EXACT? |
|---|---|---|---|---|
| `symmetry.py:1017` `_rotation_z(θ)` | `R_z(θ)` | **(M)** built from `np.cos/np.sin(θ)` | yes — `C_n` generator; expressible | **NO.** Only ever called at `θ = 2πk/n` — an exact rational multiple of 2π. Measured vs `roots_of_unity`: max Δcos `5.0e-16`; `C_4`/`C_8` get **0** exact zeros where roots-of-unity gets 2. #325 site |
| `symmetry.py:1027` `_reflections(axis)` | `σ_a = I − 2 ê_a ê_aᵀ` | **(M)** `eye(3)` with one `-1.0` | yes — `Mirror(a)`; expressible | **bit-exact** (signed permutation) |
| `symmetry.py:1041` `_inversion_op()` | `−I` | **(M)** | yes, but **NOT expressible as a tag** — there is no `Inversion`/`C_i`/`S_n` entry. Reachable only inside `I_h` | bit-exact |
| `symmetry.py:1092` `_cyclic_ops(n)` | `{R_z(2πk/n)}` | **(M)** via `_rotation_z` | `Cn(n)` | **NO** — inherits `_rotation_z` (above) |
| `symmetry.py:1097` `_vertical_mirrors(n)` | Householder `I − 2 n̂n̂ᵀ`, `n̂` at `kπ/n + π/2` | **(M)** + **(A)** (`np.cos/np.sin` on `kπ/n + π/2`) | `D_nh`'s σ_v set; expressible only as part of `Dnh(n)` | **NO** — `kπ/n + π/2 = 2π(2k+n)/(4n)`, again an exact rational multiple of 2π. #325 site |
| `symmetry.py:1132` `_octahedral_ops()` | `(±1)³ ⋊ S₃` — 48 signed permutation matrices | **(M)** assembled from `itertools.product`/`permutations` | `O_h` | **bit-exact** (this is why Lebedev's `O_h` claim is exact on both sides) |
| `symmetry.py:1157` `_icosahedral_ops()` | `I_h` = ⟨R₅, R₃⟩ ∪ `−I·⟨…⟩` | **(M)** via Rodrigues, closed by BFS with `np.allclose(atol=1e-10)` | `I_h` | **NO** — `√5` golden-ratio vertices, `cos/sin(2π/5)`, `cos/sin(2π/3)`; deduped by tolerance |
| `symmetry.py:1217` `_rotation_about_axis(a, θ)` | Rodrigues `I cosθ + sinθ K + (1−cosθ) aaᵀ` | **(M)** | the ONE arbitrary-axis rotation constructor in the tree. **NOT expressible as a tag**: `Cn(n)` is realized about **z** only | **NO** |
| `symmetry.py:911` `_is_reflection_invariant_1d` | find `π` with `−x_i = x_{π(i)}` | **(P)** built by `argsort` + `np.searchsorted`, window `atol*10` | `Mirror("x")` on the `(μ,0,0)` embedding | evaluated; **and the permutation is DISCARDED** — returns `bool` |
| `symmetry.py:1346` `_orbit_closure` | `π_M` with `M x_i = x_{π_M(i)}` | **(P)** built by full `(n,n)` distance matrix + `np.argmin`, window `atol*_NODE_WINDOW_FACTOR (=×100)` | any `ops` list | evaluated; **returns the permutation** (`OrbitCertificate`) |
| `symmetry.py:1296` `OrbitCertificate.singular_set` | `Σ = {x : Stab(x) ≠ {e}}` | **(P)** integer identity `perm == arange` | the orbifold singular set | **exact** — no tolerance. The one place in the tree where a fixed-locus is decided exactly |
| `symmetry.py:1441` `_distinct_azimuths` | count of distinct `φ` among off-axis nodes | **(A)** `arctan2` + `mod 2π` + a `1e-9` gap test + a wrap-around fix-up | bounds `C_n`/`D_nh` by divisors | evaluated, hard-coded `1e-9` (NOT `atol`) |

### B. `numerics/roots_of_unity.py` — the generating face (the exact endpoint)

| file:line | the transformation, in math | spelled TODAY | group element? | EXACT? |
|---|---|---|---|---|
| `roots_of_unity.py:232` | `4p = quad·q + r` quarter-turn split | **(O)** integer `np.divmod` | the `C_4` coset decomposition of `C_q` | exact (integer) |
| `roots_of_unity.py:238` | octant fold `r ↦ min(r, q−r)` | **(O)** integer `np.minimum` | the reflection `σ` of the first quadrant | exact (integer) |
| `roots_of_unity.py:243-245` | undo the fold: `(c,s) ↦ (s,c)` | **(P)** `np.where` coordinate **swap** | the 45° mirror | exact — a swap, not arithmetic |
| `roots_of_unity.py:251-253` | 45° diagonal fixed point `2r == q` | **(O)** integer test → `np.sqrt(0.5)` on both components | the swap's fixed point | exact by construction |
| `roots_of_unity.py:256-258` | quarter turns as `(±c, ±s)` reassignment | **(S)** `np.select` over `quad ∈ {0,1,2,3}` | `C_4` acting on `S¹` | exact |
| `roots_of_unity.py:262` | `−0.0 → +0.0` | **(O)** `+ 0.0` | signed-zero canonicalisation | exact |
| `rules_circle.py:269-274` | `shift ∈ ℚ/ℤ` reduction, then `ζ_{bn}^{bm+a}` | **(O)** `Fraction % 1` + `roots_of_unity(b·arange(n)+a, b·n)` | the rotation of the node lattice by `s` steps | **bit-exact** |

**`roots_of_unity` now has exactly ONE production consumer** — `rules_circle.periodic_trapezoid`
(`rules_circle.py:274`), reached from `rules_product.product_mu_phi`
(`rules_product.py:330`). It is still not exported from `numerics/__init__.py`.
(This supersedes the "zero production consumers" state recorded 2026-08-02.)

### C. `numerics/quadrature/` — the rules

| file:line | the transformation, in math | spelled TODAY | group element? | EXACT? |
|---|---|---|---|---|
| `directional.py:110` `_octant_sign_predicate` + `:102` `_OCTANT_SIGN_EPS=1e-15` | `ℓ(Ω) = (sign μ_x, sign μ_y, sign μ_z)` — the orbit-type stratification under the coordinate-mirror group | **(S)** three thresholded sign tests | **yes, and unnamed**: the label set is exactly the chambers ⊔ walls of `(Z₂)³ = D_2h`. `symmetry.py` CAN express the group (`Dnh(2)`, measured `\|G\|=8` = the 8 diagonal sign matrices) — the partition just doesn't say so | evaluated by `eps`; the eps is **provably idle** — measured min nonzero `\|cos\|` is `1.57e-1` (product 4×5) … `2.24e-1` (LS8), 14 orders above `1e-15`. The docstring's "keep in lockstep with `_DEGENERATE_ABS_MU_THRESHOLD`" (`:101`) is **dangling** — that symbol no longer exists anywhere in `orpheus/` or `tests/`; the surviving sibling `TANGENTIAL_EPS = 8.88e-16` (`angular_trace_space.py:164`) is a **different number** |
| `directional.py:159`, body `:209-215` `_compute_sphere_reflection_partners` | `σ_a` partner map per axis | **(M)** builds `eye(3)` with `[a,a] = -1`, then **(P)** via `_orbit_closure` | `Mirror(a)`; expressible — and this IS the one site that routes a partner map through the group machinery | evaluated (`atol=1e-13`, node window `1e-11`), but **certified** (bijection + weights) |
| `directional.py:574` `Quadrature.gauss_legendre` | slab `μ → −μ`, `partner(i) = N−1−i` | **(P)** `np.arange(N)[::-1]` | `Mirror("x")` on `(μ,0,0)` | **bit-exact** — index arithmetic, no search. Legitimate only because `GeneratingMeasure.gauss` **imposes** the mirror (below) |
| `directional.py:575-576` | `μ_y ≡ μ_z ≡ 0 ⟹ σ_y, σ_z fix every node pointwise` | **(P)** `identity` | `Mirror("y")`, `Mirror("z")` | exact, and correct: this is the same derivation `_check_invariance_1d` (`symmetry.py:876-891`) writes out |
| `directional.py:302` `axis_cosines` | the embedding `μ ↦ (μ,0,0)`; axes past `dim` → `zeros` | **(O)** column index + zero fallback | the slab-into-`S²` inclusion | exact (`np.zeros`) |
| `generating_measure.py:328-352` | **Reynolds/group-average**: `x ← (x − x[::-1])/2`, `w ← (w + w[::-1])/2` — the projection onto the `σ`-invariant subspace | **(P)** array reversal `[::-1]` | `Mirror("x")` — but the reversal is only the reflection permutation *because* `eigh` returns ascending eigenvalues. Never named | makes the result **bit-exactly** mirror-symmetric (measured defect `3.3e-16 → 0.0`) |
| `rules_sphere.py:289-294` | third level index `j = n_half − 1 − p − k` | **(O)** integer arithmetic | closure under the `S₃` coordinate permutation of `O_h` | **bit-exact** — this is what makes all 48 `O_h` ops exact (was `sqrt(1−μ²−η²)`, which realized only `D_2h`) |
| `rules_sphere.py:303-310` | 8-fold sign replication `(±η, ±ξ, ±μ)` | **(S)** triple `for s in (-1,1)` | `(Z₂)³ = D_2h` — the octant group again, generated here and classified in `directional.py:110`, with no shared name | **bit-exact** |
| `rules_sphere.py:327` | level membership `\|μ_z\| == level_mu[p]` | **(O)** float **equality** (was `tol=1e-12`) | the fiber of the invariant `\|μ_z\|` | **exact** — the values are copied out of `mu_levels` |
| `rules_sphere.py:328` | intra-level order by `η` | **(P)** `np.argsort(mu_x[idx])`, default `quicksort`, **no `kind=`** | none — `η = sinθ cos φ` is **2-to-1 on the fiber** (documented at `rules_sphere.py:164-175`; issue #326) | evaluated; unstable sort over genuine ties |
| `rules_sphere.py:335` | recover the fiber angle `φ = arctan2(μ_y, μ_x) mod 2π` | **(A)** | the chart on the level circle | evaluated; branch cut at `φ=0` |
| `rules_sphere.py:204-218` `LevelStructure.fiber` | order the fiber by `(hemisphere, azimuth)` | **(P)** `np.lexsort` | an ordering of the circle `S¹` — the *correct* one, unlike `level_indices` | evaluated |
| `rules_product.py:356-362` | `(μ, φ) ↦ (sinθ cos φ, sinθ sin φ, μ)` — the embedding into `S²` | **(A)/(O)** `sqrt(1−μ²)` × the roots-of-unity components | the parametrisation of `S²`, NOT a group element. This is the pushforward the theory page (`discrete_measures.rst:362-386`) flags as missing from the composition algebra | `cos φ`, `sin φ` exact; `sinθ` evaluated |
| `rules_product.py:382` | intra-level order by `η` | **(P)** `np.argsort(..., kind="stable")` | none (see above) | evaluated, but the tie-break is now **named** |
| `rules_product.py:425-428` | `φ = arctan2(sin φ, cos φ) mod 2π`, tiled per level | **(A)** round trip back to the angle chart | the chart | evaluated (`≤8.9e-16`); order preserved |
| `rules_product.py:428` | `hemisphere = sign(μ_z)` | **(S)** `np.sign` | the `σ_z` orbit label | exact (`μ_z` is a GL node, never 0 for even `n_mu`) |

### D. `numerics/basis/` + the moment layer

| file:line | the transformation, in math | spelled TODAY | group element? | EXACT? |
|---|---|---|---|---|
| `spherical_harmonic_basis.py:44-51` (prose) + `:421-427` (code) | **the polar axis is `μ_x`, azimuth in the `(μ_y, μ_z)` plane** — i.e. every `Y_ℓ^m` is the standard harmonic composed with a **frame rotation** `R : (x,y,z) ↦ (y,z,x)` | **(O)** a coordinate *convention*, realized by which array is assigned to `cos_theta` | **YES, and this is the big one.** Measured: `R` has `det=+1`, order 3, axis `(1,1,1)/√3`, angle 120° — the 3-fold rotation of `O_h`. It **IS** in `_octahedral_ops()` (measured `True`). It is **NOT expressible** as a `symmetry.py` tag: `Cn(3)` is realized about **z** and does **not** contain `R` (measured `False`); only `OctahedralOh` contains it, as 1 of 48 unaddressable elements | the matrix is a bit-exact signed permutation — but it is never applied *as a matrix*; it is applied by choosing which column to call `cos θ`, so nothing can check it |
| `spherical_harmonic_basis.py:415-417` | `Y_1^{-1}=μ_z`, `Y_1^0=μ_x`, `Y_1^{+1}=μ_y` | **(O)** hardcoded | the same `R`, spelled a *second* time and consistently | exact |
| `spherical_harmonic_basis.py:427` + `:435-436` | `φ = arctan2(sin φ, cos φ)`, then `cos(mφ)`, `sin(mφ)` | **(A)** | de Moivre: `cos(mφ) = T_m(cos φ)`, `sin(mφ) = sin φ·U_{m−1}(cos φ)` — the `C_m` action on the circle | **NO, and needlessly.** The code already holds `(cos φ, sin φ)` and throws them away through the angle chart. Measured on `product(4,8)`: `\|cos(m·arctan2) − T_m(cos φ)\| = 1.1e-16 / 3.9e-16 / 5.0e-16` at `m=1/2/3`, and the Chebyshev route keeps **8 exact zeros** at `m=1` and `m=3` where the `arctan2` route keeps **0**. At `m=1` it is a pure identity round trip that loses bits |
| `spherical_harmonic_basis.py:423-426` | on-axis guard `sinθ < 1e-15` → `(cos φ, sin φ) := (1, 0)` | **(O)** gauge fixing at the chart singularity | choosing a representative of the `SO(2)` fiber over a pole | evaluated; **and it does not fire for a slab** — see Finding 3 |
| `numerics/moment_layout.py` (whole file) | — | — | — | **GAP**: this module is *spatial* tensor-Legendre layout only. There is **no angular moment-layout module**; the `(ℓ, m)` slot convention (`slot = l+m`, entries outside `\|m\|≤ℓ` zero) is restated in ≥3 places: `spherical_harmonic_basis.py:196-199`, `harmonic_moment_flux.py:13`, `material_xs_field.py:964`, `scattering.py:198` |
| `radial_characteristic_space.py:522-583` `fold_moments_to_radial_characteristic` | `Q̄(μ=±1) = Σ_ℓ ((2ℓ+1)/2) Q_ℓ (±1)^ℓ` | **(S)** `sign**l` | the **antipodal map's action on Legendre moments** — the 1-D Wigner-D at the inversion element | exact (`(±1)^ℓ`), and the **only** place in the tree where a group element's action on a MOMENT vector is written down |
| `transport/operators/scattering.py:911-955` | `Q^aniso = (1/W)·R Λ M ψ` | **(O)** operator composition | the addition theorem `Σ_m Y_ℓ^m(Ω)Y_ℓ^m(Ω') = P_ℓ(Ω·Ω')` — i.e. the `SO(3)`-invariance of the scattering kernel — is the *reason* moments diagonalise `Λ`. Named in prose (`:922`, and `spherical_harmonic_basis.py:24-27`), never as a transformation | n/a |
| **anywhere** | a **Wigner-D matrix** / a rotation applied to a moment vector | — | — | **DOES NOT EXIST.** `grep -rni wigner orpheus/` returns only "Wigner–Seitz" (a pin-cell geometry, unrelated). If a moment rotation were needed today, the only ingredients are `_rotation_about_axis` (private, `symmetry.py:1217`) and re-tabulating `Y` at rotated directions — there is no `D^ℓ_{m'm}` anywhere |

### E. `sn/sweep/` + the curvilinear fiber march (the angular half)

| file:line | the transformation, in math | spelled TODAY | group element? | EXACT? |
|---|---|---|---|---|
| `geometry/reduced_operator.py:688-690` (sphere) | `α_{n+1/2} = α_{n−1/2} − w_n μ_n` | **(O)** a Python `for` loop cumulative sum over the ordinate ORDER | the march around the polar interval. Permutation-invariant as a telescoping sum (my L-014) | evaluated; guarded by `assert \|α_N\| < 1e-12` |
| `geometry/reduced_operator.py:778-786` (cylinder) | same, per level, over `η` sorted ascending | **(O)** loop over `level_indices` | **the march around the fiber CIRCLE** — and no rotation, reflection or wrap is spelled anywhere. The circle is cut open into an ordered list by `argsort(η)`, which is 2-to-1 on it (#326) | evaluated |
| `geometry/reduced_operator.py:806-809` | level start `η_{1/2} = −sinθ_p = −√(1−μ_z²)` | **(A)** | the `φ = π` point of the level circle | evaluated |
| `pole_angular_closure.py:576-586` (sphere) | `μ_edge` from `−1` by cumulative `w`; `τ = (μ−μ⁻)/(μ⁺−μ⁻)` | **(O)** loops | the angular cell edges on `[−1,1]` | evaluated; `1e-15` degenerate fallback |
| `pole_angular_closure.py:593-611` (cylinder) | `η_edge`: `−sinθ`, midpoints, `+sinθ` | **(A)/(O)** | the fiber circle cut at `φ=π` and laid on an interval `[−sinθ, +sinθ]`; **the circle's periodicity is discarded here** | evaluated |
| `pole_angular_closure.py:1012-1030` `_edge_seed_stencil` | find the 2 most-inward **distinct-μ** ordinates | **(P)** `np.argsort(mu)` + a `\|Δμ\| > 1e-14` distinctness scan | selecting an orbit representative on the `η`-tie | evaluated; a **third** distinctness epsilon (`1e-14`), different from `1e-15` and `8.88e-16` |
| `pole_angular_closure.py:253-271` `_gather_per_ordinate` | per-level `(M_p,)` → global `(N,)` | **(P)** `out[level] = values` | the level partition as a permutation | exact |
| `transport/radial_characteristic_field.py:314` | level's most-inward ordinate | **(P)** `ords[argmin(mu[ords])]` | orbit-representative selection | evaluated (no guard); the same quantity `_edge_seed_stencil` finds a different way |
| `loss_representation/__init__.py:3294` | pole continuation `ψ(0, +μ) = ψ(0, −μ)` | **(P)** `outflow_at_inner.T[quad.reflection_index("x")]` | `Mirror("x")` — correctly routed through the certified table | exact index gather |
| `loss_representation/__init__.py:3470, 4189, 4593` | the same, in the transpose / matvec arms | **(P)** `mirror = quad.reflection_index("x")` | `Mirror("x")` | exact gather |
| `sn/sweep/psi_half_angle_seed.py:125, 303` | the inward (`μ=−1`) leg marched by reversing the radial arrays | **(P)** `[:, ::-1]` on `Q`, `σ`, `dr` | the radial orientation reversal, not an angular group element (it pairs with the antipodal `sign ∈ {−1,+1}` leg index) | exact |
| `numerics/spaces/angular_trace_space.py:202-251` `build_omega_dot_n` | `Ω·n̂_f` | **(S)** `sign * mu[axis]` — the outward normal is an `(axis, sign)` pair, never a vector | the face normal is `±ê_a`; the projection is a coordinate **selection**, and the ± IS the `min/max` half of the mirror convention | exact given the cosines |
| `numerics/spaces/angular_trace_space.py:431, 440` | `Γ₋ / Γ₊ / tangential` three-way split | **(S)** `row < −TANGENTIAL_EPS`, `row > +TANGENTIAL_EPS` | the `σ_a`-fixed locus is exactly the tangential set — the same `Σ` `OrbitCertificate.singular_set` decides **exactly** | evaluated (`8.88e-16`) |
| `transport/fields/angular_boundary_flux.py:159` | `J = Σ_m sign(Ω·n̂)·\|Ω·n̂\|w_m ψ_m` | **(S)** `np.sign(row) * metric` | reconstructing the signed projection from `\|·\|·w` and a sign | exact |

### F. Out-of-scope-but-adjacent (MoC) — listed only because they are #325 sites

| file:line | transformation | spelled | exact? |
|---|---|---|---|
| `moc/quadrature.py:89` | `φ_a = π(2a+1)/(2 n_azi)` on the QUOTIENT circle `[0,π)` | **(A)** `np.linspace + π/(2 n_azi)` | **NO** — `= 2π(2a+1)/(4 n_azi)`, exactly `roots_of_unity(2a+1, 4·n_azi)` |
| `moc/geometry.py:319-320` | `cos φ_a`, `sin φ_a` from the stored angle | **(A)** | **NO** — same |
| `moc/geometry.py:222-229` `_reflected_azi_index` | `φ ↦ π − φ` on `[0,π)` | **(A)** + **(P)** `argmin(\|φ − φ_refl\|)`, **no distance guard** | the answer is exactly `n_azi−1−a` by index arithmetic; computed by search instead |
| `moc/geometry.py:412-437` `_find_link` | nearest track endpoint | **(P)** `best_dist` never thresholded | evaluated, unguarded |

---

## The 7 questions

### 1. Is there ANY rotation of moments (Wigner-D)? What exists if one were needed?

**No Wigner-D anywhere.** `grep -rni 'wigner' orpheus/` returns only *Wigner–Seitz*
pin-cell geometry. No `D^ℓ_{m'm}`, no Euler angles, no moment-space rotation
operator. The moment layer is rotation-free: `SphericalHarmonicBasis`
(`spherical_harmonic_basis.py`) exposes `evaluate`/`synthesize`/`analyze`/
`reconstruct`/`mass_matrix` and nothing that acts on the `m` index.

But a rotation **is** already applied to the SH basis, silently: the project's
polar axis is `μ_x`, not `μ_z` (`spherical_harmonic_basis.py:44-51`, code at
`:421-427`). Every ORPHEUS `Y_ℓ^m` is the textbook harmonic composed with the
frame change `(x,y,z)_std = (μ_y, μ_z, μ_x)`. Measured: that map is a proper
rotation of order 3 about the `(1,1,1)` body diagonal (120°) — a genuine `O_h`
element, present in `_octahedral_ops()`, and **not expressible as a
`symmetry.py` tag** (`Cn(3)` is about `z` and does not contain it). Because it is
applied by *which column is called `cos θ`* rather than by a matrix, no
invariance check, no adjoint, and no test can see it.

If a moment rotation were needed today the ingredients are: `_rotation_about_axis`
(private, `symmetry.py:1217`, Rodrigues, trig-evaluated) to build the matrix, and
re-tabulating `Y` at `nodes @ R.T` — i.e. rebuild-and-reproject, never a
closed-form `D^ℓ`. The nearest existing thing to a moment-space group action is
`fold_moments_to_radial_characteristic` (`radial_characteristic_space.py:522`),
which applies `P_ℓ(±1) = (±1)^ℓ` — the antipodal map on Legendre moments.

### 2. Is the OCTANT partition a group-theoretic object or ad-hoc sign classification?

It **is** the orbit-type stratification of `(Z₂)³ = D_2h`, and it is spelled
ad-hoc. `directional.py:110` `_octant_sign_predicate` labels each node by
`(sign μ_x, sign μ_y, sign μ_z)` against `_OCTANT_SIGN_EPS = 1e-15`
(`directional.py:102`); `Quadrature.octants` (`directional.py:493`) hands it to
`DiscreteMeasure.partition_by`.

Measured, the label set is exactly "chambers ⊔ walls" of the coordinate-mirror
group:

| rule | full-sign chambers | strata with a `0` component (= walls) | total parts |
|---|---|---|---|
| `lebedev(17)` | 8 | 18 | 26 |
| `product(4,8)` | 8 | 8 | 16 |
| `level_symmetric(8)` | 8 | 0 | 8 |
| `gauss_legendre(8)` | 2 | 0 | 2 |

So "octants" is a misnomer for every rule except LS_N: the partition returns 26
entries for Lebedev. The zero components are precisely membership in
`Fix(σ_a)` — the same set `OrbitCertificate.singular_set` (`symmetry.py:1296`)
decides by the **exact integer identity** `π_M(i) == i`. The group is
expressible today: `SubgroupOfO3.Dnh(2)` realizes exactly the 8 diagonal sign
matrices (measured `|G| = 8`), and `registry.py:686-693` already assigns
`Dnh(2)` as the owed residual for `cylinder`/`cartesian2d`. The partition just
does not name it, and asks the membership question with a float threshold that
is provably idle (min nonzero `|cos|` across the shipped rules: `1.57e-1`, i.e.
14 orders above the eps).

The docstring pointer at `directional.py:100-101` — "*Matches the pure-z
degenerate-ordinate threshold in `orpheus.sn.sweep` and
`orpheus.transport.spatial.diamond` (`_DEGENERATE_ABS_MU_THRESHOLD`); keep in
lockstep*" — is **dangling**: `_DEGENERATE_ABS_MU_THRESHOLD` exists nowhere in
`orpheus/` or `tests/`. The surviving sibling `TANGENTIAL_EPS`
(`angular_trace_space.py:164`) is `8.88e-16` ≠ `1e-15`.

### 3. Does anything BESIDES `_orbit_closure` compute "the permutation induced by a transformation on a node set"?

Yes — **five** other sites, in three kinds.

**(a) A second matching engine inside `symmetry.py` itself.**
`_is_reflection_invariant_1d` (`symmetry.py:911-940`) finds the `−x` partner by
`np.argsort` + `np.searchsorted` + a two-candidate `min`, with window `atol*10`
— against `_orbit_closure`'s `(n,n)` distance matrix + `np.argmin` with window
`atol*100` (`symmetry.py:1393-1395`). Two engines, two windows, and the 1-D one
**throws the permutation away** (returns `bool`). Both measured working on their
own inputs; nothing compares them.

**(b) Formula endpoints (correct, exact, unshared).**
`Quadrature.gauss_legendre` (`directional.py:574`) uses `arange(N)[::-1]` —
`partner(i) = N−1−i` by index arithmetic. `rules_sphere.py:289-294` gets the
third level index by `j = n_half − 1 − p − k`. `generating_measure.py:351-352`
uses `[::-1]` as the reflection permutation (implicitly, via `eigh`'s ascending
order). None of these route through, or is checked against, `_orbit_closure`.

**(c) Unguarded nearest-neighbour searches.**
`moc/geometry.py:229` `_reflected_azi_index` — `argmin(|φ − φ_refl|)`, no
distance threshold (the answer is `n_azi−1−a` in closed form).
`moc/geometry.py:412-437` `_find_link` — nearest track endpoint, `best_dist`
never compared to anything.
`transport/radial_characteristic_field.py:314` and
`pole_angular_closure.py:1019-1030` both find "the level's most-inward ordinate"
— one by bare `argmin(mu)`, the other by `argsort` plus a `1e-14` distinctness
scan. Same quantity, two spellings, two (or zero) tolerances.

The tree also already owns the *typed* endpoint nobody in the angular layer
uses: `PermutationOperator` (`numerics/operator.py:2172`) — inverse via
`argsort`, involution detected by `perm[perm] == arange`, exact.

### 4. Curvilinear angular redistribution: is any rotation/reflection of the fiber circle spelled explicitly?

**No.** The fiber over a polar level is a circle and the march is around it, but
every step treats it as an **ordered interval**:

- `reduced_operator.py:778-786` — `α_{m+1/2} = α_{m−1/2} − w_m η_m` accumulated
  in a Python loop over `level_indices[p]`, whose order came from
  `argsort(mu_x)` (`rules_sphere.py:328` / `rules_product.py:382`).
- `pole_angular_closure.py:593-611` — the level's angular edges are built as
  `η_edge = [−sinθ_p, midpoints…, +sinθ_p]`, i.e. the circle **cut at
  `φ = π`** and laid flat. Neither endpoint is identified with the other; the
  periodicity is simply dropped.
- `reduced_operator.py:806-809` — the starting direction is `−sinθ_p`, the
  `φ = π` point, computed as `−√(1−μ_z²)`.

There is no `φ → φ + 2π/n` rotation, no `φ → −φ` reflection of the level, and no
wrap. The only reflection that *does* appear on the fiber is the **pole
continuation** `ψ(0,+μ) = ψ(0,−μ)`, spelled as an index gather through the
certified table (`loss_representation/__init__.py:3294, 3470, 4189, 4593`:
`quad.reflection_index("x")`) — that one is exact and correctly named `Mirror("x")`.

Two consequences the ordering choice carries, both already documented in-tree:
the sort key `η = sinθ cos φ` is **even in φ**, hence 2-to-1 on the fiber
(`rules_sphere.py:164-175`; issue #326), and `LevelStructure.fiber`
(`rules_sphere.py:204-218`) exists as the *correct* `(hemisphere, azimuth)`
lexsort but is **not** what the α-recursion consumes.

### 5. Where are `cos`/`sin` evaluated on exact rational multiples of 2π? (#325's remaining surface)

Complete enumeration for production `orpheus/` (excluding `derivations/`,
`plotting.py`, and the docstring mentions). Three of the six are **inside the
symmetry checker itself** — which is exactly the "remaining half" the
`rules_product.py:50-56` caution names.

| # | site | angle | exactly `roots_of_unity`-able? |
|---|---|---|---|
| 1 | `symmetry.py:1019` `_rotation_z` ← `_cyclic_ops` (`:1094`) | `2πk/n` | **YES, directly** — `roots_of_unity(arange(n), n)` |
| 2 | `symmetry.py:1126` `_vertical_mirrors` | `kπ/n + π/2 = 2π(2k+n)/(4n)` | **YES** — `roots_of_unity(2k+n, 4n)` |
| 3 | `symmetry.py:1220` `_rotation_about_axis` ← `_icosahedral_ops` (`:1182, :1188`) | `2π/5`, `2π/3` | **partially** — the `cos θ`/`sin θ` yes; the axis carries `√5` and a normalisation |
| 4 | `spherical_harmonic_basis.py:435-436` | `mφ`, `φ` recovered from the direction cosines | **partially, and a better fix exists**: use `T_m(cos φ)` / `sin φ·U_{m−1}(cos φ)` from the components the code already holds — no angle, no `arctan2` (`:427`), no branch cut |
| 5 | `moc/quadrature.py:89` | `π(2a+1)/(2 n_azi)` | **YES** — `roots_of_unity(2a+1, 4·n_azi)` |
| 6 | `moc/geometry.py:319-320` | `cos φ_a`, `sin φ_a` from the stored angle | **YES** — same numerators |

Not #325 sites: `mc/solver.py:412-413` (random directions, not rational);
`fuel/`, `kinetics/`, `thermal_hydraulics/`, `geometry/factories.py`,
`geometry/mesh.py` (all `linspace` over radii/lengths, no trig).

Already closed: `rules_product.py` (the azimuthal factor, `:330`, via
`rules_circle.py:274`) — the one repoint that has landed.

Measured impact of the open sites: `_cyclic_ops(4)` and `_cyclic_ops(8)` produce
**zero** exact zeros where `roots_of_unity` produces two; the `Dnh(n_φ)` check
therefore has trig on **both** sides, which is what
`symmetry._NODE_WINDOW_FACTOR = 100` (`symmetry.py:1241`) absorbs. Item 4
measured on `product(4,8)`: the `arctan2`→`cos(mφ)` route keeps **0** exact
zeros where the Chebyshev-from-components route keeps **8** (at `m=1` and `m=3`).

### 6. Is there any place that builds a 3×3 rotation matrix by hand?

Every hand-built 3×3 in production `orpheus/` lives in **two files**
(`grep -rn 'np\.eye(3)'`):

- `symmetry.py:1017` `_rotation_z` — `R_z(θ)`, trig.
- `symmetry.py:1029` `_reflections` — `eye(3)` with one `−1.0`.
- `symmetry.py:1043` `_inversion_op` — `−eye(3)`.
- `symmetry.py:1127` `_vertical_mirrors` — Householder `eye(3) − 2 n̂n̂ᵀ`, trig normal.
- `symmetry.py:1150-1153` `_octahedral_ops` — 48 signed permutation matrices assembled slot-by-slot.
- `symmetry.py:1221-1226` `_rotation_about_axis` — Rodrigues; the ONLY arbitrary-axis constructor.
- `directional.py:209-210` — `reflection = np.eye(3); reflection[axis, axis] = -1.0`, i.e. a **local re-spelling of `_reflections(axis)`** two modules away (it imports `_orbit_closure` from `symmetry` but not `_reflections`).

Nowhere else. In particular there is **no** hand-built rotation in `basis/`,
`transport/`, `sn/`, or `moc/` — those layers apply transformations as index
permutations or sign flips only.

### 7. `angular_frame` / `axis_cosines`: what IS the embedding convention, and is it stated once?

The convention is: **a quadrature's node array is a `(N, d)` block of direction
cosines whose column `i` is `Ω·ê_i`; a 1-D polar marginal embeds as `(μ, 0, 0)`;
missing columns are zero.**

It is stated in **six** places, and they agree — but a seventh, the SH basis,
uses a *different* frame and that disagreement is invisible:

1. `directional.py:302-324` `axis_cosines` — the canonical accessor; zero-pads
   past `nodes.shape[1]`.
2. `directional.py:330-354` — `mu_x/mu_y/mu_z` as views on axes 0/1/2.
3. `directional.py:356-396` — the cylindrical aliases `η/ξ/μ` = columns 0/1/2,
   with the explicit warning that `mu_x` is *misleading* in cylindrical context.
4. `directional.py:446-487` `angular_frame` — re-states it operationally by
   `np.column_stack([axis_cosines(0), axis_cosines(1), axis_cosines(2)])` with
   `support=SPACE_SPHERE`, i.e. **a slab measure on `[-1,1]` is re-declared as
   living on `S²`** at frame-construction time.
5. `symmetry.py:876-891` `_check_invariance_1d` — *derives* the 1-D mirror arm
   from the same embedding ("`(μ,0,0)` ⟹ `σ_y`, `σ_z` fix every node
   pointwise, `σ_x` is the one real test"), citing `axis_cosines`.
6. `Mirror`'s docstring (`symmetry.py:198-213`) and
   `discrete_measures.rst:528-536` restate it a third and fourth time in prose.

The disagreement: `spherical_harmonic_basis.py:44-51` declares "*The polar axis
is `μ_x` … azimuth is measured in the `(μ_y, μ_z)` plane*". That is a
**different frame** from the one every other consumer assumes (`μ_z` axial —
see `directional.py:347-354`, `rules_product.py:20-27`, `registry.py`'s
`Dnh(n_φ)` about `z`, and `symmetry.py`'s entire `C_n`/`D_nh`/`σ_h` realization,
all about **z**). It is self-consistent and it composes correctly with the
`(μ,0,0)` slab embedding — which is presumably why it was chosen — but it means
the SH basis and the symmetry checker measure azimuth about **perpendicular
axes**. There is no crosswalk between the two statements anywhere.

---

## The 3 findings that most surprised me

### 1. The spherical-harmonic basis silently applies an `O_h` element that `symmetry.py` cannot name — and it is applied by *variable choice*, not by a matrix

`spherical_harmonic_basis.py:421-427` sets `cos_theta = mu_x` and measures
azimuth in the `(μ_y, μ_z)` plane. Measured, the implied frame change
`(x,y,z)_std = (μ_y, μ_z, μ_x)` is the **120° rotation about the `(1,1,1)` body
diagonal** — proper (`det = +1`), order 3, and a member of `_octahedral_ops()`.
`SubgroupOfO3` can express `Cn(3)`, but that is realized about **z** and
measurably does *not* contain this element; only `OctahedralOh` does, as one of
48 elements with no individual handle.

What makes this the most surprising row in the audit is not the rotation — it is
that the whole rest of the angular layer measures azimuth about **z**
(`_cyclic_ops`, `_vertical_mirrors`, `Dnh(n_φ)`, `LevelStructure.azimuth`,
`hemisphere = sign(μ_z)`), so the basis and the checker are working in frames
related by a group element neither of them mentions. And because the rotation is
applied by choosing which array to call `cos θ`, there is no matrix for
`is_invariant` to test, no adjoint to check, and no place a rename would break.

### 2. The "octant partition" is a 26-entry orbifold stratification wearing an 8-chamber name, decided by an epsilon that is 14 orders of magnitude too loose to matter — and whose only justification is a dangling pointer

`Quadrature.octants` (`directional.py:493`) is documented as the eight-way
sign decomposition. Measured: `lebedev(17)` → **26** parts (8 chambers + 18
walls), `product(4,8)` → **16**, `gauss_legendre(8)` → **2**. The label's zero
components mark exactly `Fix(σ_a)` — the singular set — and `symmetry.py` already
computes that set **exactly** (`OrbitCertificate.singular_set:1296`, integer
identity `π_M(i) == i`, no tolerance), while `registry.py:686-693` already names
the group (`Dnh(2)`, measured `|G| = 8` = the coordinate sign matrices). Three
pieces of the same object, none aware of the others.

And the threshold: `_OCTANT_SIGN_EPS = 1e-15` is justified by
"*Matches … `_DEGENERATE_ABS_MU_THRESHOLD`; keep in lockstep*"
(`directional.py:100-101`) — a symbol that **exists nowhere** in `orpheus/` or
`tests/`. The surviving sibling `TANGENTIAL_EPS = 8.88e-16` is a different
number, and a third, `1e-14`, sits in `pole_angular_closure.py:1023`. Measured
minimum nonzero `|cos|` over the shipped rules: `1.57e-1`. All three epsilons are
idle, and the one comment defending them points at a ghost.
(This is the delegation-shaped-claim pattern: prose that hands a responsibility
to a named mechanism elsewhere, where the mechanism is what's missing.)

### 3. `#325`'s remaining surface is mostly *inside the symmetry checker*, and the moment layer throws away exactness the quadrature just bought

I expected the open `roots_of_unity` repoints to be rule-side (MoC, the SN
azimuth). The SN azimuth is **done** (`rules_circle.py:274` ← `rules_product.py:330`),
and three of the six remaining sites are `symmetry.py`'s own operators:
`_rotation_z` (`:1019`), `_vertical_mirrors` (`:1126`), `_rotation_about_axis`
(`:1220`). Measured, `_cyclic_ops(4)`/`(8)` produce **zero** exact zeros where
`roots_of_unity` produces two — so the `Dnh(n_φ)` invariance claim on the product
rule now has an exact node set checked by an inexact operator set, and
`_NODE_WINDOW_FACTOR = 100` absorbs the difference from one side only.

The sharper one is in the moment layer. `spherical_harmonic_basis.py:427` takes
the two components `(cos φ, sin φ)` it has just computed, collapses them through
`arctan2` into an angle, and then re-evaluates `np.cos(m*phi)` / `np.sin(m*phi)`
(`:435-436`). At `m = 1` that is a pure identity round trip. Measured on
`product(4,8)`: the `arctan2` route keeps **0** exact zeros; the de-Moivre
route from the same components (`T_m(cos φ)`, `sin φ·U_{m−1}(cos φ)`) keeps
**8** at `m=1` and **8** at `m=3`, differing by `1.1e-16`…`5.0e-16`. The exact
azimuths bought by the `roots_of_unity` carve are destroyed one layer above the
quadrature, by a chart the code did not need to enter.

### Bonus — a falsified docstring found while measuring the above (worth a separate look)

`Quadrature.spherical_harmonics` (`directional.py:434-438`) claims: "*For slab
GL1D quadratures only the `m = 0` harmonics `Y_l^0(μ_x)` carry non-zero values;
the other slots are filled with zeros.*"

Measured on `Quadrature.gauss_legendre(8).spherical_harmonics(3)`:

```
 l=0:  m=0: 1.00e+00
 l=1:  m=-1: 0.00e+00   m=0: 9.60e-01   m=+1: 0.00e+00        <- true
 l=2:  m=-2: 0.0  m=-1: 0.0  m=0: 8.83e-01  m=+1: 8.34e-01  m=+2: 8.37e-01
 l=3:  m=-3..-1: 0.0     m=0: 7.73e-01  m=+1: 8.04e-01  m=+2: 7.37e-01  m=+3: 7.51e-01
```

Mechanism: for a slab, `axis_cosines(1) = axis_cosines(2) = 0` but
`sinθ = √(1−μ_x²) ≈ 0.9`, so the `on_axis` guard at
`spherical_harmonic_basis.py:423` **never fires**; `(cos φ, sin φ)` becomes
`(0, 0)` — not a point of `S¹` — and `arctan2(0.0, 0.0) = 0.0`, so
`cos(mφ) ≡ 1` for every `m`. Only the `m < 0` slots (which carry `sin(mφ)`) are
zero. `ℓ ≤ 1` is clean because it is hardcoded (`:411-418`); the defect begins
exactly where the chart is entered.

Consequence, measured with the module's own einsums (`analyze` then
`reconstruct`) on a random axisymmetric slab flux at `L=3`: the reconstruction
with the full table differs from the `m = 0`-only reconstruction by a factor of
**4.4** (`max|ΔR| = 6.40`, relative `3.44`). The discrete Gram over the same
table shows the spurious slots carrying twice the mass of the true one
(`⟨Y_2^0,Y_2^0⟩ = 0.4 = 2/5` correct; `⟨Y_2^{+1},Y_2^{+1}⟩ = ⟨Y_2^{+2},Y_2^{+2}⟩ = 0.8`).

Why nothing caught it: the only `scattering_order ≥ 2` test in the tree
(`tests/sn/operators/test_scattering_operator.py:1584`) uses
`Quadrature.lebedev(17)`, a genuine 3-D cubature where the chart is
well-defined. The slab `P_{≥2}` path is fixture-unreachable. I am reporting the
measurement, not a verdict on the end-to-end physics — that needs a
`numerics-investigator` with the slab `P_ℓ` normalisation chain in hand.
