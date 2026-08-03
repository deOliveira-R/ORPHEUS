# Geometric transformation machinery — consolidation campaign

> **STATUS: MAPPING.** Audits in flight. Nothing implemented. This file is the
> plan of record; it is written to survive a compaction and be picked up cold.
>
> **Branch** `refactor/operator-strategy-layers`. Base commit at authoring:
> `bfedc621` (Q5.0.2 — `Z2` retired, `Mirror(axis)` minted).

---

## 0. The user's directives, verbatim in substance (2026-08-03)

1. **Consolidate the geometric-transformation machinery in ONE appropriate
   place, and deploy it from that single place.**
2. **Sequence is fixed and non-negotiable:**
   1. consolidate the machinery,
   2. **verify it mathematically**,
   3. **only then** clear the tests that need transformations to lean on it.
3. **This runs BEFORE Q5.1.** Rationale, verbatim in substance: *"it will never
   be as cheap as today to do this"* — the campaign currently holds the whole
   symmetry/quadrature context, and that advantage decays.
4. After mapping + sequencing, **compact**, then start implementation with clean
   context.

### Why step 2 must precede step 3 — the principle, already established

This is the project's **verify-a-shared-primitive-against-pure-math** ruling
applied to transformations: verify the shared primitive against
domain-independent ground truth FIRST; after that, a system-under-test and its
oracle BOTH using it is *not* contamination, because **independence lives in
the INPUT ASSEMBLY, not in the shared routine.**

Inverted — tests migrated onto UNVERIFIED machinery — the tree loses its
independent check silently: every gate would then agree with production because
both call the same unproven code. That is the failure this ordering exists to
prevent, and it is why "migrate the tests" is a *later* step and not a
bookkeeping tail.

---

## 1. What the machinery IS today (read directly, not inferred)

Three faces, all in `orpheus/numerics/`:

| face | where | what |
|---|---|---|
| **generate** | `roots_of_unity.py` | circle points BY the group action; mirrors bit-exact |
| **realize** | `symmetry.py::_realized_ops` | tag → finite matrix set; `_close_group` closes a generating set |
| **act & certify** | `symmetry.py::_orbit_closure` | applies ops to nodes, returns the induced PERMUTATION in an `OrbitCertificate` |

Primitives, all **private**, all trading in bare `np.ndarray`:

| symbol | `symmetry.py:` | what |
|---|---|---|
| `_rotation_z(theta)` | 1017 | 3×3 about z |
| `_reflections(axis)` | 1027 | the three COORDINATE mirrors, as a diagonal sign flip |
| `_cyclic_ops(n)` | 1092 | `C_n` about z |
| `_vertical_mirrors(n)` | 1097 | `D_nh`'s vertical planes — **Householder `I − 2 n̂ n̂ᵀ`, built inline at :1127** |
| `_octahedral_ops()` | 1132 | `O_h` |
| `_icosahedral_ops()` | 1157 | `I_h` |
| `_rotation_about_axis(axis, theta)` | 1217 | **a GENERAL Rodrigues rotation** |
| `_close_group(ops)` | 617 | closure under multiplication |
| `_orbit_closure(...)` | 1346 | → `OrbitCertificate` (`:1245`) carrying `permutations` |

`_rotation_x` / `_rotation_y`: **DEAD** — measured zero consumers tree-wide.

### The two observations that motivate the campaign

1. **The general reflection already exists, unnamed.** `_vertical_mirrors`
   builds `I − 2 n̂ n̂ᵀ` inline (`:1127`) — reflection about an ARBITRARY plane —
   while `_reflections(axis)` builds the three coordinate mirrors as a diagonal
   sign flip, and `Mirror(axis)` can name ONLY those three. The tree constructs
   the object its type system cannot express. Same shape for rotation:
   `_rotation_about_axis` is general, private, untyped.
2. **Every transformation crosses a boundary as a bare `np.ndarray`.** An
   element of O(3) has composition, an inverse, a determinant, an order, an
   axis-or-plane, and an action on points; today each is implicit or recomputed
   per site. This is anti-pattern #13 landing on the most algebraic object in
   the codebase.

---

## 2. Audit findings

### 2a. Boundary layer — LANDED, `scratch/geom_transforms_audit_boundary.md`

**⛔ It corrected a claim I had made to the user.** The 3-D specular path does
NOT reflect by hand: `_compute_sphere_reflection_partners`
(`directional.py:207-216`) builds the matrix and delegates to
`symmetry._orbit_closure`, bit-identical to `orbit_certificate(...).permutations`
by measurement. That unification landed in Q5.0.1 (`a7695148`) the same day.

The real findings:

* **FOUR type vocabularies for σ_a that do not know about each other** —
  `SpecularMirror` (`_factors.py:474`), `SpecularReturn`, `SpecularReemission`
  (`_factors.py:724`), and `symmetry.Mirror`. `_factors.py` **never imports
  `symmetry`**.
* **One remaining hand-rolled site**, and it is the counter-intuitive one:
  the **1-D** path (`directional.py:570-577`, `identity[::-1].copy()`) is the
  ONLY BIT-EXACT reflection in the layer, and it is the one that BYPASSES the
  machinery. The principled 3-D route decides by `argmin` at `atol=1e-13`
  (measured residual `0.0`, runner-up margin 0.26–1.16 — the map is not at
  risk; the tolerance lives in the CERTIFICATION, not the permutation).
* **A live defect with the wrong error type.** `product(4,5)` (odd `n_φ`) fails
  the certificate on axis 0, so `reflection_index` (`directional.py:416-424`)
  raises a bare **`ValueError`**, not a `BoundaryError`. Same odd-`n_φ` family
  as ERR-074.
* **A third spelling as an integer ±1 key swap**: `_reflect_corner`
  (`sn/operators/boundary.py:928-933`) applies the mirror to the
  off-quadrature `μ=±1` ray by hand, and **unscaled** (documented defect at
  `:897-907`).
* ⭐ **A REAL BOUNDARY ON THE WHOLE IDEA.** The periodic BC's `SpatialWrap`
  (`_factors.py:523-627`) is a **TRANSLATION** `x ↦ x + L`. `SubgroupOfO3` is
  origin-fixing, so it **structurally cannot express it**. *Point group ≠ space
  group.* Any "one transformation machinery" must decide whether it models
  O(3) or the full Euclidean group — a DESIGN decision, not a refactor.
* **Latent consumer:** `Cn`/`Dnh` are fully realized with **zero** BC-layer
  consumers; the nine rotation mentions there are all prose pointing at issue
  **#178 `SymmetryBoundary`** (the rotational/sector BC that does not exist).
* Partial duplication of `_orbit_closure`'s weight guard: `_specular.py:114-163`
  vs `symmetry.py:1403`. And `_face_domains` (`sn/operators/boundary.py:361`)
  is `_orbit_closure`'s permutation-or-refuse SHAPE on the FACE set.
* `is_involution` (`numerics/operator.py:2243`) reports `False` on
  `lebedev(17)` — mathematically right (a `Γ₊→Γ₋` bijection between different
  sets) while `operator.py:2188-2191` asserts the opposite in prose. **Zero
  production consumers**, so not a live bug.

### 2b. Angular / quadrature — LANDED, `scratch/geom_transforms_audit_angular.md`

**`symmetry.py` is the ONLY place in the whole tree where a transformation is a
matrix.** Seven constructors live there; exactly one re-spelling exists outside
(`directional.py:209-210`, a local `eye(3)`+`−1` copy of `_reflections`, in a
module that already imports `_orbit_closure`).

* ⭐⭐ **THE SH BASIS APPLIES AN `O_h` ELEMENT BY VARIABLE CHOICE, AND THE
  MACHINERY CANNOT NAME IT.** `spherical_harmonic_basis.py:421` sets
  `cos_theta = mu_x`, so every project `Y_ℓ^m` is the textbook harmonic composed
  with a rotation. `[M]` That rotation has `det=+1`, **order 3**, axis
  `(1,1,1)/√3`, 120° — an element of `_octahedral_ops()` — and `Cn(3)`
  (realized about z) measurably does NOT contain it; only `OctahedralOh` does,
  as 1 of 48 unaddressable elements. **Every other angular consumer is about
  z** (`_cyclic_ops`, `Dnh(n_φ)`, `hemisphere = sign(μ_z)`, `rules_product`).
  Six places state the `(μ,0,0)`/z-polar convention and agree; this is the
  seventh and it disagrees. **No crosswalk exists.** Because the rotation is
  applied by *which array gets called* `cos θ`, there is no matrix to test, no
  adjoint to check, and a rename breaks nothing.
* ⭐⭐ **A LIVE FALSIFIED DOCSTRING with a measured consequence → GitHub #328.**
  `Quadrature.spherical_harmonics` (`directional.py:434-438`) claims that for a
  slab GL rule only `m=0` slots are non-zero. `[M]` **Independently reproduced
  by the main agent**: `m≠0` slots carry **0.74–0.88**. Mechanism: for a slab
  `sinθ≈0.9` so the `on_axis` guard (`spherical_harmonic_basis.py:423`) never
  fires, `(cos φ, sin φ)` becomes `(0,0)` — not a point of `S¹` — and
  `arctan2(0,0)=0`, so `cos(mφ) ≡ 1`. Audit measurement (not re-verified):
  full-table vs `m=0`-only reconstruction differ by **4.4×**, and the discrete
  Gram shows the spurious slots at **twice** the true mass. Fixture-unreachable:
  the only `scattering_order ≥ 2` test uses `lebedev(17)`.
* ⭐ **"Octants" is a 26-part orbifold stratification wearing an 8-chamber
  name.** `[M]` lebedev(17) → **26** parts (8 chambers + 18 walls),
  product(4,8) → 16, LS(8) → 8, GL(8) → 2. It IS the `(Z₂)³ = D_2h` orbit-type
  stratification (`directional.py:102,:110`), spelled as three thresholded
  signs. The group is ALREADY named 
  (`registry.py:686-693`, `Dnh(2)`) and the walls are ALREADY computed exactly
  by `OrbitCertificate.singular_set` — **three views, mutually unaware**.
  Its `_OCTANT_SIGN_EPS = 1e-15` is justified by a **GHOST SYMBOL**
  (`_DEGENERATE_ABS_MU_THRESHOLD`, which exists nowhere in `orpheus/` or
  `tests/`); `[M]` min genuine `|cos|` across shipped rules is `1.57e-1`, so all
  three epsilons in this family (`1e-15`, `TANGENTIAL_EPS`, `1e-14`) are idle.
* ⭐ **#325's remaining surface is mostly INSIDE THE CHECKER.** Six open sites;
  three are `symmetry.py`'s own operators. `[M]` `_cyclic_ops(4)`/`(8)` produce
  **zero** exact zeros where `roots_of_unity` produces two — so the product
  rule's `Dnh(n_φ)` claim now has an **exact node set checked by inexact
  operators**, absorbed by `_NODE_WINDOW_FACTOR = 100` from one side only.
  Full list: `_rotation_z` ← `_cyclic_ops`; `_vertical_mirrors`;
  `_icosahedral_ops` (partial, axis carries `√5`);
  `spherical_harmonic_basis.py:435-436`; `moc/quadrature.py:89`;
  `moc/geometry.py:319-320`.
* **FIVE permutation engines besides `_orbit_closure`**, incl. a second matching
  engine *inside the same file* (`_is_reflection_invariant_1d`, window `atol×10`
  vs `×100`, and it **discards π** to return a bool).
* **`_inversion_op` (−I) is NOT expressible** — no `C_i`/`S_n` tag; reachable
  only inside `I_h`. `_rotation_about_axis` is likewise **not taggable** (`Cn`
  is about z only).
* **The Reynolds/group average is already in the tree, unnamed**:
  `generating_measure.py:328-352` symmetrises via `(x − x[::-1])/2`, which is
  `Mirror("x")` averaging — `[M]` it is what turns `3.3e-16` into `0.0`.
* **`roots_of_unity` now has exactly ONE production consumer** and is still not
  exported from `numerics/__init__.py`.
* No Wigner-D anywhere; the only group action on a moment vector is
  `radial_characteristic_space.py:522-583`'s antipodal `(±1)^ℓ`. There is **no
  angular moment-layout module** — the `(ℓ,m)` slot convention is restated in
  ≥4 places.


### 2c. Spatial / mesh — LANDED, `scratch/geom_transforms_audit_spatial.md`

**Verdict: the spatial layer speaks PERMUTATIONS and SIGN VECTORS, not 3×3
matrices — which is precisely why none of it can be spelled in `symmetry.py`'s
vocabulary today.** Zero hand-built transformation matrices in the entire
spatial scope.

* ⭐⭐ **THE CATEGORY BOUNDARY IS THE TRANSLATION, NOT THE REFLECTION.** The
  expected explanation ("`symmetry.py` speaks `(N,3)` nodes, the mesh speaks
  index lattices") was **MEASURED FALSE**: `_orbit_closure` fed a **centred**
  1-D cell lattice under `_reflections("x")` returns the permutation
  `[5 4 3 2 1 0]` — exactly the `np.arange(n)[::-1]` that
  `sweep_graph.py:577-580` writes by hand. *The machinery already produces the
  spatial answer.* What kills it is `Mesh1D.from_geometry` starting at
  `origin = 0.0` (`geometry/mesh.py:405,:490,:516`), so
  `Mirror('x').is_invariant(mesh.volume_measure)` is **False** for the
  production slab AND sphere, **True** only for a hand-centred `[-1,1]`.
  **The spatial group is `E(d) = O(d) ⋉ ℝ^d`, and the affine part has NO TYPE
  ANYWHERE** (`origin: float`, an inline `- pitch/2` at `factories.py:117`, and
  the half-cell's mirror plane is "wherever x=0 happens to be").
  Sharper: for CYL/SPH, `r ∈ [0,R]` is the **fundamental domain**, never
  G-invariant — `False` is the CORRECT answer, and it says the layer needs a
  **quotient vocabulary**, not an invariance predicate.
* ⭐ **`octant_moment_frame_signs` (`transport/spatial/_ubld.py:96-148`) is a
  full GROUP REPRESENTATION and nothing says so** — the 1-D character rep of
  `(Z₂)^d` on the tensor-Legendre basis, verified a homomorphism
  (`χ(g)χ(h) == χ(gh)` over all `2^d × 2^d` pairs at d=1,2,3, exactly). The
  docstring and the theory page say only "involution" — true per element, a
  whole structural level short. **It was born from a bug** (ERR-061); naming the
  rep would have made that bug *unspellable* rather than findable-by-debugging.
* **The sweep reversal is spelled ELEVEN times in four idioms**, and the one
  done RIGHT is the newest: `_reverse_octant_traversal`
  (`loss_representation/__init__.py:750-787`) + `_CellResidualTranspose`
  (`sweep_graph.py:1046-1090`) map `o ↦ −o` and look up the mirror octant's
  graph, leaving the walk untouched — "the whole reversal is DATA". Meanwhile
  four methods within ~100 lines of `augmented_mesh.py` (`:1275,:1398,:1426,
  :1499`) each re-write the same `range` ternary, and `scan.py:302-363` takes a
  literal `x_reverse: bool`.
* The r=0 pole continuity `ψ(0,+μ)=ψ(0,−μ)` IS a named, permutation-realized
  deck transformation (`loss_representation/__init__.py:3294` et al.) — but via
  `quad.reflection_index("x")`, a hard-coded axis literal.
* The chart map is **never written down** — only the integrated Jacobian
  (`geometry/coord.py:46-135`), whose docstring claims to be "the single point
  where coordinate-system dependence lives"; three exceptions listed in the
  audit. Coordinate identity is carried in **four** parallel encodings.
* Dangling Python-domain xrefs to `SNMesh._setup_spherical`/`_setup_cylindrical`
  (deleted): `geometry/reduced_operator.py:16-17,:653,:747` and
  `docs/theory/foundations/structured_geometry.rst:395-396`. Sphinx `-W` cannot
  see these.

### 2d. Cross-solver — LANDED, `scratch/geom_transforms_audit_cross_solver.md`

**Headline: no non-SN module builds a rotation or reflection matrix at all.**
Trig appears on five lines in all non-SN production.

* ⭐⭐ **MoC HAND-ROLLS `_orbit_closure`, TWICE, GUARD-FREE.**
  `_reflected_azi_index` (`moc/geometry.py:222-229`) applies `σ: φ ↦ π−φ` then
  `np.argmin(np.abs(phi - phi_refl))`; `_find_link` (`:412-435`) does the
  spatial half by squared-distance argmin. `symmetry.py:1383-1407` does the
  identical job **with** a distance window, a bijectivity check (ERR-073) and a
  weight check. MoC has none — and `_find_link` computes `best_dist`
  (`:422,431-432`) and **discards it**.
* ⭐⭐ **MoC's azimuthal set IS the rule `numerics` shipped hours ago.**
  `φ_k = π(k+½)/n_azi = 2π(2k+1)/(4·n_azi)` is exactly the upper half of
  `periodic_trapezoid(2·n_azi, shift=STAGGERED)`. `[M]` over
  `n_azi ∈ {2..64}`: identical node sets to ≤5.0e-16, weights differing by
  exactly π (normalisation only). The `numerics` route makes the mirror
  `k ↦ n−1−k` **bit-exact at every `n_azi`**; MoC's `linspace+cos` at **none**.
* **CP's method of images is an addition**: `gap_c = bnd_pos[i] + bnd_pos[j]`
  (`cp/solver.py:276,:354`) is the reflection `s ↦ −s`, marked only by a
  three-word comment.
* **MC has no reflection at all** — `BC_REGISTRY` holds only `"periodic"`; and
  no lab-frame rotation exists because direction is re-sampled absolutely at
  every collision (`mc/solver.py:409-413`).
* `d − 2(d·n)n` appears **once in the whole scope**, as PROSE at
  `derivations/continuous/peierls_nystrom/origins/specular/r_matrix.py:48-52`,
  never coded — while its matrix form already exists at `symmetry.py:1127`.
* **Measured negatives**: the MoC link map IS a bijection today on all four
  configurations tested, so the missing guard does not currently bite — but the
  geometric match is loose, the reflected track's entry point sitting up to
  **8.96e-02 cm** from the exit point (**~2.5 % of the 3.6 cm pitch**),
  unreported. And `tests/moc/test_verification.py:556-588` **cannot see any of
  it**: for any total map on a finite set the orbit repeats within
  `n+1 < max_steps`, so its cycle assertion is satisfiable by construction.
* `diffusion/` and `homogeneous/` contain no geometric transformation.

---

## 3. SYNTHESIS — what the four audits say together

**The machinery is not missing. It is UNREACHABLE.** Every constructor a
consumer could want already exists in `symmetry.py` — the general Householder
(`:1127`), the general Rodrigues (`:1217`), `O_h`'s 48 signed permutations
(`:1132`), inversion (`:1041`). What is missing is a **vocabulary** the
consumers can speak.

### The vocabulary mismatch, measured

`symmetry.py` trades in **3×3 matrices acting on `(N,3)` nodes**. Its consumers
speak four *different* languages, and none of them is matrices:

| consumer | speaks | example |
|---|---|---|
| BC / specular | **permutations** of ordinates | `reflection_index(axis)` |
| spatial / sweep | **permutations** of a cell lattice + **sign vectors** | `np.arange(n)[::-1]`, `Π_a s_a^{o_a}` |
| MoC | **angle arithmetic** + argmin search | `φ ↦ π − φ` |
| moment layer | **sign powers** on `ℓ` | `(±1)^ℓ` |

`_orbit_closure` is the ONE bridge (matrix → permutation) and it is already
correct — three guards, incl. the ERR-073 bijectivity check. **The consolidation
is therefore mostly about making that bridge REACHABLE, not about writing new
math.**

### The three findings that shape the design

1. **The category boundary is the TRANSLATION.** `[M]` `_orbit_closure` on a
   *centred* cell lattice already returns the spatial answer; what defeats it is
   `origin = 0.0`. The spatial group is `E(d) = O(d) ⋉ ℝ^d` and **the affine
   part has no type anywhere**. *Point group ≠ space group* — this is a scope
   decision, not a refactor. (Periodic's `SpatialWrap` is the same boundary from
   the BC side.)
2. **A transformation can hide as a CONVENTION.** The SH basis applies a real
   `O_h` element by choosing which array is called `cos θ`. There is no matrix,
   so there is nothing to test — and it has already produced a live falsified
   docstring (#328). *A transformation applied by naming rather than by
   multiplication is invisible to every gate.*
3. **A representation is not an element.** `octant_moment_frame_signs` is the
   character rep of `(Z₂)^d` — `[M]` verified a homomorphism — and is documented
   only as "an involution". It was **born from ERR-061**. `symmetry.py` today
   has no way to say "this is a representation of that group".

### What is genuinely duplicated (Cardinal Rule 2 targets)

| the operation | spellings | worst offender |
|---|---|---|
| `σ_a` on ordinates | **7+** | 4 BC vocabularies + `_reflect_corner`'s ±1 key swap |
| `_orbit_closure`'s job | **5 clones** | MoC's two guard-free `argmin`s |
| sweep-sense reversal | **11**, in 4 idioms | 4 `range` ternaries in ~100 lines of one file |
| `(Z₂)³` octant classification | **3 mutually-unaware views** | generated, classified, and named in 3 modules |
| the equispaced circle | **2** | MoC's azimuths ARE `periodic_trapezoid(2n, STAGGERED)` |

---

## 4. Open design questions — ANSWER BEFORE WRITING CODE

1. **Scope: `O(3)` or `E(3)`?** The translation is outside `O(3)`. Decide
   explicitly. A defensible middle: keep the consolidated core point-group, and
   give the affine part its own small type so `origin`/`pitch/2` stop being bare
   floats — without claiming a space-group implementation.
2. **Where is "one appropriate place"?** `symmetry.py` mixes THREE concerns:
   the transformation primitives, the subgroup lattice, and the invariance
   checker. Consolidation likely means **extracting the primitives into their
   own module** that `symmetry.py` then consumes — which also gives
   `roots_of_unity` a natural home next to them.
3. **What is the TYPE?** Candidate: a value type carrying composition, inverse,
   determinant, order, and axis-or-plane, with `Reflection(normal)` /
   `Rotation(axis, angle)` constructors and the coordinate mirrors as the exact
   special case. Must answer: does it also carry its **action** (matrix) vs its
   **realization on a given node set** (permutation)?
4. **Does `Mirror(axis)` widen to `Reflection(normal)`?** Minted hours ago with
   a deliberately narrow `{x,y,z}`. Widening is cheapest NOW.
5. **Exactness tiering.** Coordinate mirrors and signed permutations are
   bit-exact; general Householder/Rodrigues are not. The type should make the
   tier LEGIBLE, not uniform. Closing the three checker-side `roots_of_unity`
   sites is in scope here (#325's remaining half).
6. **Do representations get a name?** (finding 3 above). Possibly out of scope —
   but decide, do not drift.
7. **What is the mathematical verification (the step-2 gate)?** Against PURE
   MATH, never another ORPHEUS path: group axioms (closure / associativity /
   identity / inverse), `MᵀM = I`, `det = ±1`, involution for reflections,
   order-`n` for `C_n`, Householder and Rodrigues against independent
   constructions, conjugation `R σ R⁻¹`, and the orbit–stabiliser count on a
   known orbit. Every law mutation-verified.

---

## 5. Sequencing

| step | what | gate |
|---|---|---|
| **G0** | answer §4 with the user | a ruling recorded here |
| **G1** | extract/mint the primitives in ONE place; widen `Mirror`; export `roots_of_unity` | tree green; realizations bit-identical |
| **G2** | ⭐ **mathematically verify** them against pure math | the law suite, every law mutation-verified |
| **G3** | route `symmetry.py`'s seven constructors through the consolidated core | realizations bit-identical |
| **G4** | close the checker-side `roots_of_unity` sites (`_cyclic_ops`, `_vertical_mirrors`) | `Dnh(n_φ)` check exact on BOTH sides |
| **G5** | route the BC layer's 4 σ_a vocabularies + `_reflect_corner` | tree green; the `product(4,5)` `ValueError` becomes a `BoundaryError` |
| **G6** | MoC: replace the 2 guard-free `argmin`s with the certificate; adopt `periodic_trapezoid` | the 8.96e-2 cm link gap becomes visible/asserted |
| **G7** | ⭐ **migrate the tests** onto the verified machinery | coverage preserved; no gate weakened; re-run each migrated gate's justifying mutation |
| **G8** | retire superseded spellings + dead code (`_rotation_x/_y`, the `directional.py:209` re-spelling, the ghost-justified epsilons) | grep clean |

**Deliberately NOT in scope** (file/track, do not fix here): #328 (the SH
`m≠0` defect), the SH polar-axis crosswalk, `#178 SymmetryBoundary`, the
11-fold sweep-reversal cleanup, the dangling `_setup_spherical` xrefs.

---
## 6. Standing constraints for the implementation session

* `main` is always green; never commit to `main`; `--ff-only`; no squash-merge.
* **NEVER** `git checkout <path>` / `git restore <path>` / `git stash` /
  `git clean` on a tracked path — `.claude/skills/vv-principles/{SKILL.md,error_catalog.md}`
  carry uncommitted IRRECOVERABLE state and are **forbidden to commit**. Commits
  must be path-scoped to exclude them. To revert a mutation, copy-to-tmp and
  restore, or monkeypatch in-process.
* Host env: `.venv/bin/python`; canonical `python -O -m pytest`; **xdist is
  UNSTABLE — run SERIAL**.
* Commit trailer: `Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>`.
* Temp files under `$CLAUDE_JOB_DIR/tmp`, not `/tmp`.
* **L37**: while a long pytest gate runs, Python sources are FROZEN.
* **Mutation testing MUST reach the test process.** Monkeypatching the parent
  and running pytest in a SUBPROCESS reads GREEN for every mutation — it cost a
  full round in Q5.0.2. Mutate the SOURCE (copy to tmp, mutate, restore), and
  check the `-k` filter is not deselecting the very tests that should catch it.

### The red-set baseline — SCOPED, and the scoping bit me once

* `tests/sn -m "not slow"`: **6 deliberate reds**, documented in `81689a58`.
* `tests/geometry`: **a 7th**, `test_bc_equivalence_snapshot.py::
  TestWhiteXminPartial03GLSnapshot::test_matches_the_frozen_scaled_lambertian`
  — `[M]` 8/60 elements off by 2.2e-16 abs / 1.1e-15 rel vs `rtol` 8.9e-16.
  Bisected: fails at `292a1ba5` too, so it predates E3/E4/Q5.0.2; likely origin
  is Q3's GL re-baseline (`bc89b62e`). **Every "the same six" report in the
  quadrature campaign is scoped to `tests/sn` — run `tests/geometry` too.**
* pyright: **1 error**, the accepted #288 residual in
  `transport/operators/scattering.py:757`.
