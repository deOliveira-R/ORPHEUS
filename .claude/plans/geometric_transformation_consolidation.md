# Geometric transformation machinery — consolidation campaign

> **STATUS: G1–G5 · G6.0 · G6.1 LANDED. G6.3 steps 0–5 DONE; step 6 DONE (absorbed by the
> affine campaign's P3 — bind `ZeroOperator`'s spaces); step 8.0 LANDED `72f5ce97`
> (see §7h.4 — it turned out to carry a live `.H` repair, not just plumbing);
> **step 8 LANDED** (§7h.5 — the check now runs on the production path, and the
> transpose leg is DERIVED rather than written).
> ✅ **§7d.3 LANDED 2026-08-07** (`c7ca338e` test tier → successor;
> `b5ac130e` the table DELETED + docs/label/markers migrated — see the
> §7d.3 landing record for the audit trail). **STEP 7 IS COMPLETE.**
> Steps 7a–7c AND 7d.0–7d.2 landed the same day
> (`e13313a8` … `f11a129f`): `PairedDeck` retired `SpatialWrap`; ONE
> kernel derives every deck arrow from the motion via
> `Quadrature.ordinate_permutation` (⛔ the earlier `RigidMotion.permutes`
> spelling here was the AFFINE action — ordinates are directions); the
> certifications and the sweep pole seed read the same source; the
> reflection table is GONE — `ordinate_permutation` is structurally the
> only path.
> ▶ **RESUMES AT G6.3b — the flagship metric gate on the factored
> Lambertian** (task #45; ⛔ its original "absorb
> `AngularAverageOperator.apply_transpose`" half is a STALE MECHANISM
> TITLE — that operator was retired at step 3b, 2026-08-04, factored
> into `IsotropicEmissionOperator @ PartialCurrentOperator` where each
> link has one honest transpose, so the deferral DISSOLVED; re-derive
> the remaining scope from the tree per §7.2 before designing), then
> G6.5 — and track G is done, unblocking Q5.1.
> **Read ⏸ COMPACTION POINT #6 (END of file) FIRST — it carries the
> G6.3b hazards, the owed wide slice, and the six-commit table — then
> §7d.3's landing record and the G6 design rows :1496 / :1905-1945.**
>
> ⏸ #6 → #5 → #4 → #3 → #2 (each additive; #5's 7d.3 items and #4's
> step-7 items are HISTORY but
> its gate costs and lessons stand). #3 covers what changed under track G
> *while G was paused* (the Γ ladder grew a `directions` field; the
> `level_symmetric` weights changed and the family now REFUSES above S12) —
> ⚠ #3's "Step 8.0 — the next step" section is HISTORY, it landed. Then ⏸ #2 (end of
> §7h) for the red set and the five-times-refuted-scope lesson; step 5's own findings
> are §7h.3. **#3 does not supersede #2 — it is additive.**
> This file is the plan of record; it is written to survive a compaction and be
> picked up cold. Verify every hash below with
> `git merge-base --is-ancestor <hash> HEAD` before trusting this header — it is
> a snapshot, and the `NEXT =` field is the one that rots (see the main memory
> index's three index disciplines, 2026-08-03).
>
> | phase | commit(s) | what landed |
> |---|---|---|
> | **G1** | `6acb6a8a` | `orpheus/geometry/transformation.py` — the dimension-generic rigid-motion core |
> | **G1b** | `cc5703ed` | the `E(d)→E(D)` embedding, the `Permutation` type, `NotAFinitePointGroupError` |
> | **G2** | `387fe625` | verified against PURE MATH — 42 gates / 96 cases, 32 mutations / 32 caught / 0 blind |
> | **G3** | `3bf029de` + `0c3d4e65` | `numerics/symmetry.py` speaks `RigidMotion`; the 1-D/3-D arm split deleted; the 1-D invariance tag corrected to `Mirror("x")` |
> | **G4** | `cfe0a448` (§7d) | the checker's `C_n` / `σ_v` re-based on `roots_of_unity`; mirrors become the coset `C_n·σ₀`; `_rotation_z` retired |
> | **G4.5** | `af99f064` (§7d) | the index-permutation route to an exact `0.0` — proposed, measured, **abandoned** with its reason |
> | **G5** | `3cce383a` + `241097b3` (§7e) | `SelfPairedDeck` minted and verified, then `IdentityMap` + `SpecularMirror` retired onto it; two decaying gates re-posed in the same commit |
>
> **G6 (in progress)** — ⭐ **every operator knows its domain, codomain and space;
> the SPACE owns shape and traversal.** The adjoint then falls out of
> well-posedness instead of being hand-rolled, and a layout change happens in ONE
> place with zero edits to any mathematical operation. `[M]` half of this is
> ALREADY built and correct (`_AdjointOperator` derives `A† = G_V⁻¹AᵀG_W` from the
> spaces alone) — the defect is that the binding is OPTIONAL, so it silently
> degrades to the Euclidean transpose.
>
> **G6.0 is DONE**; read **§7g** for what it measured and the two user rulings that
> re-scoped the step. Headlines: the Γ ladder is **three** tiers, not two; **G** is
> metric-blind in every shipped case while **R** is where the metric bites (and R
> is an endomorphism, so it can be self-adjoint); binding comes from the
> **realizer per law, never the class**; **no operator is silently wrong today**,
> so G6 closes a trap rather than fixing a live defect; and the tree-wide
> **mandate moved to #330**. §7f is the original scope with its falsified claims
> corrected in place. Steps G7/G8/G9 are the former G6/G7/G8, renumbered.
>
> **Branch** `refactor/operator-strategy-layers`. Base commit at authoring:
> `bfedc621` (Q5.0.2 — `Z2` retired, `Mirror(axis)` minted).
> **Merge-status: trust git, not this header** — `git merge-base --is-ancestor`.

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

`symmetry.py` trades in **3×3 matrices on `(N,3)` nodes**. Its consumers speak
four other languages: **permutations** of ordinates (BC) and of cell lattices
(spatial), **sign vectors** (the octant rep), **angle arithmetic** (MoC), and
**sign powers on ℓ** (moments). `_orbit_closure` is the ONE bridge
(matrix → permutation) and it is already correct — three guards, incl. the
ERR-073 bijectivity check. **So the carve is mostly about making that bridge
reachable, not about writing new math.**

### Duplication (Cardinal Rule 2 targets)

| operation | spellings | worst offender |
|---|---|---|
| `σ_a` on ordinates | **7+** | 4 BC vocabularies + `_reflect_corner`'s ±1 key swap |
| `_orbit_closure`'s job | **5 clones** | MoC's two guard-free `argmin`s |
| sweep-sense reversal | **11**, 4 idioms | 4 `range` ternaries in ~100 lines of one file |
| `(Z₂)³` octants | **3 mutually-unaware views** | generated, classified, named in 3 modules |
| the equispaced circle | **2** | MoC's azimuths ARE `periodic_trapezoid(2n, STAGGERED)` |

---

## 4. RULINGS — settled with the user, 2026-08-03

### R1. The core is `RigidMotion`, and the affine part is FORCED

Not a preference — a proof. **You cannot describe a reflection or a rotation
*in space* without the affine part.** A hyperplane is `{x : n̂·x = d}`; only
`d = 0` passes through the origin. A rotation needs a centre. Only the
origin-fixing special case is linear.

**The relation between affine and non-affine is ONE operation** — conjugation
by a translation:

```
Q seated at c  ≡  translate(c) ∘ Q ∘ translate(−c)  =  (Q, (I − Q)c)
```

which reproduces every case:

| element | geometric description | `Q` | `t` |
|---|---|---|---|
| reflection | normal `n̂`, offset `d` | `I − 2n̂n̂ᵀ` | `2d·n̂` |
| rotation | plane, angle `θ`, centre `c` | `R` | `(I − R)c` |
| inversion | centre `c` | `−I` | `2c` |
| translation | vector `v` | `I` | `v` |

Composition `(Q₁,t₁)∘(Q₂,t₂) = (Q₁Q₂, Q₁t₂+t₁)`; inverse `(Qᵀ, −Qᵀt)`.

⭐ **This is load-bearing for a consumer that is broken TODAY.** `[M]`
`Mirror('x').is_invariant(mesh)` reads **False** for the production slab and
sphere, because meshes start at `origin = 0.0` while the mirror plane sits at
`x = a`. That is `d ≠ 0` — unexpressible today, one field away from
expressible.

### R1b. Two axes, not one — and the scope argument had to be replaced

⛔ **Corrected 2026-08-03, on the user's challenge.** The original wording ran
two independent questions together as "rigid vs affine". They are:

| | question | answer | argued by |
|---|---|---|---|
| **Axis 1** | is the **translation** in the core? | **YES — forced** | the proof above: you cannot say *where* a mirror is without it |
| **Axis 2** | is a **non-isometric linear part** in the core? | **NO** | see below |

`[M]` Axis 1 is implemented and working: `reflection(normal=x, offset=2.5)`
carries `t = [5,0,0]`, is not linear, and fixes its own plane; the production
cell lattice reads `None` under a mirror at `x=0` and
`Permutation([5,4,3,2,1,0])` under the same mirror at `x=½`. `[M]` Axis 2 is
enforced: shear, uniform scale and anisotropic scale are all rejected.

**The Axis-2 justification was resting on an argument this plan had already
disavowed.** It read *"No consumer needs shear or scaling"* — the very
"weak argument" the affine proof was introduced to replace. And the affine proof
does not cover it: it argues for **including translation** and says nothing
about **excluding shear**. Replaced with the structural argument:

**Orthogonality is not a convenience — it is what makes "symmetry" a well-posed
question.**

* A symmetry **preserves the structure**; for a weighted point set that
  structure is distance and volume. A non-isometry carries the node set to a
  *different* node set, so `permutes` / `preserves` would not be answering a
  hard question — they would be answering a malformed one.
* A point group **is by definition** a subgroup of `O(d)` seated at a point,
  not of `GL(d)`.
* `close_group` of a non-isometry is generically **infinite** (scaling by 2
  generates ℤ), so the finite-group machinery loses its meaning and
  `NotAFinitePointGroupError` would fire on ordinary input.
* `inverse()` is `Qᵀ` — **exact**. A general linear part makes it a solve
  carrying conditioning, and every G2 law that closes at ~1e-15 today
  (`g∘g⁻¹ = e`; adjoint-of-uplift = uplift-of-inverse) degrades to
  tolerance-bound.

The consumer count survives only as a footnote, never as the argument.

**Where a genuinely non-rigid affine map would go.** Not here. The real
candidate in this tree is the **reference-cell → physical-cell map** behind the
LD/UBLD tensor-Legendre basis, which carries an anisotropic scaling. That is a
**chart** — a different object with a different job — and it belongs beside the
curvilinear charts `(r,θ) ↔ (x,y)`, which are **non-linear** and stay outside as
a different concept, not a degenerate case.

### R2. Elements are DIMENSION-GENERIC, parameterised by the complement of what they fix

| | reflection | rotation |
|---|---|---|
| fixed set | dim `n−1` (hyperplane) | dim `n−2` |
| **parameterised by the complement** | 1-D → a **normal** | 2-D → a **rotation plane** + angle |
| det | −1 | +1 |

ℝ²: a rotation fixes the origin (dim 0). ℝ³: fixes a line — which is why we say
"axis", the normal of the rotation plane. For `n ≤ 3` every rotation is
*simple*, so this covers us completely.

⭐ **This DISSOLVES the 1-D/3-D arm split rather than reconciling it.**
`symmetry.py` is 3-D-locked (`np.eye(3)`, `SubgroupOfO3`), which is why
`_check_invariance` has two arms — 1-D data must be forced through the
`(μ,0,0)` embedding. Q5.0.2 made the arms *agree* by deriving the 1-D case from
that embedding; a dimension-generic `Reflection` makes ℝ¹ first-class (normal
`[1]` → `[−1]` → `x ↦ −x`) and the split is deleted, not reconciled.

⚠ **Naming trap, live in the tag shipped 2026-08-02.** `Mirror("z")` reads in
English as "through the z axis" = `diag(−1,−1,1)`, `det = +1` — which is
`Cn(2)`, a rotation, and is *also* expressible. Only the docstring
disambiguates. Fix: `normal=` keyword-only, so the wrong reading is unspellable.

### R3. Placement — `geometry`, and the layer test PERMITS it

⛔ **A claim I made and had to retract.** I asserted `numerics → geometry` was
forbidden, from a 0-import count plus the layer *description*. The rule says
otherwise: `FORBIDDEN_EDGES["numerics"] = L2 | L3` — transport and the solver
families only. The docstring is explicit: *"Imports flow only L3 → L2 → L1 …
**and any layer → input**."* `geometry` is an input layer; **every** layer may
import it, `numerics` included. The 0 count was current usage, not a rule.

And the dependency runs the way the user argued: **geometry determines its point
group → that drives quadrature selection**, which `GEOMETRY_ANGULAR_SYMMETRY`
already encodes. The inverse ("choose a quadrature, then decide the geometry")
is absurd. Putting the walk in numerics puts the upstream concept downstream of
its own consumer.

```
geometry/transformation.py   dimension-generic elements + closure
                             + the point-set orbit primitive   (numpy only)
geometry/  (symmetry)        the 3-D point-group lattice + candidates + the walk
numerics/                    the WEIGHT-aware wrapper only
```

`SubgroupOfO3` keeps its name — the named entries (`O_h`, `I_h`, `D_∞h`,
`SO(3)`) genuinely are 3-D, even though the elements below them are not.

### R4. `permutes` vs `preserves` — two predicates, not an optional flag

`_orbit_closure`'s three guards separate cleanly. Guards 1–2 (image-is-a-node,
**bijection**) need only points; guard 3 (weights match) is the measure
specialisation, and it is ONE LINE.

```
g.permutes(points)    -> permutation  | None    "does g map this set onto itself?"
g.preserves(measure)  -> certificate  | None    "...and match the weights?"
```

`preserves` composes `permutes` + one guard. An optional `weights=None` would be
a boolean flag changing *which question is asked* — anti-pattern #3.

**And geometry wants `preserves`, not `permutes`**: a mesh has weights too — the
**cell volumes**. A mesh whose centres are mirror-symmetric but whose volumes
are not is not symmetric for any integration purpose.

### R5. The uplift is BY CONSTRUCTION, not by method

An adjoint is undetermined until the domain and codomain are. So the uplift
cannot be `element.as_operator()` — that would also make geometry name an
operator type. Instead a **law type** binds the element to a space pair:

```
ReflectionLaw(element, domain, codomain)     # built where the spaces live
```

taking its operation from `geometry/transformation.py`. Geometry never names an
operator; no cycle. Verification law worth building around: for orthogonal
`Q⁻¹ = Qᵀ` and permutation `P⁻¹ = Pᵀ`, so

```
T.realize_on(μ).as_operator().H  ≡  T.inverse().realize_on(μ).as_operator()
```

**the adjoint of the uplift is the uplift of the inverse element.**

### R6. Representations — build the permutation family, DOCUMENT the character family

The five collapse into two structurally different families.

**A — permutation reps on finite point sets** (ordinates, faces-as-normals, cell
lattices): 3 of 5. The bridge exists and is correct. **BUILD.**

**B — character/sign reps on moment bases** (`(±1)^ℓ` on Legendre moments;
`octant_moment_frame_signs`' `(Z₂)^d`): 2 of 5, and *not permutations* —
diagonal sign operators indexed by moment order. **DEFER the type**, by the
project's own type-minting criterion (b): no non-identity morphism is applied
(nothing composes, restricts, tensors, or takes characters of them).

**But fix the prose NOW** — `octant_moment_frame_signs` should say "the
character rep of `(Z₂)^d` on the tensor-Legendre basis", not "an involution".
One line, and ERR-061 is the evidence for what not saying it costs.

**Revisit trigger:** when the **angular moment-layout module** is built (there
is none; the `(ℓ,m)` convention is restated in ≥4 places). That is where the
character reps get a home *and* a second morphism (restriction to the octant
subgroup) — the ≥2-instances-plus-a-morphism moment.

### R7. `Mirror` is the SURFACE; `Reflection` is what it produces; the ACTION is separate

User's vocabulary, adopted: **Mirror** = the surface (normal + offset);
**Reflection** = the transformation it produces; **Specular** = one *action* it
can carry.

⛔ **A MAP I GOT WRONG IN CONVERSATION — corrected here so it does not
propagate.** I placed white/albedo as `G = reflection, R = Lambertian/αI`. FALSE.

| BC | G (deck transformation) | R (constitutive response) |
|---|---|---|
| specular / reflective | the reflection | identity |
| **white** | **identity** | Lambertian (rank-one) |
| **albedo** | **identity** | `αI` |
| specular-albedo | the reflection | `αI` |
| periodic | **translation** | identity |
| sector (#178) | **rotation** | identity |

The code was already right (`white.py:114-132`, `albedo.py:142-158`:
`IdentityMap` unconditionally), and the campaign's own theorem names the trap I
fell into: *"G is unobservable exactly when R is rank-one — the theorem that hid
the Lambertian in the geometry slot."*

⭐ **The user's physical framing, which sharpens the membership rule.** Albedo
and white are **real mirrors**, not conceptual ones — their deck transformation
is the identity and the mirror physics lives in the RESPONSE. Specular is the
*ideal* limit whose response is a delta in the mirror direction; Lambertian is
the perfectly-diffuse limit; **a real Ni/Ti supermirror neutron guide is a
peaked-but-not-delta kernel in between**, which needs a genuine `R` — not a
scalar, not a rank-one.

⟹ **The same reflection belongs in G or in R depending on whether the domain is
actually quotiented.** Specular at a *symmetry plane* IS a deck transformation
(the half-domain really is a quotient). Specular at a *physical mirror surface*
is NOT — nothing is quotiented; it is constitutive. Same math, different slot.
This is the campaign's stated sufficient test applied to the guide case, and it
is why neutron-guide modelling lands in `R`.

⭐ **And G need not be a reflection — the BC family is ONE law parameterised by
a rigid motion:** `ψ_in(x,Ω) = ψ_out(g⁻¹x, Q_g⁻¹Ω)`, with `g` a reflection
(specular), a translation (periodic), or a rotation (sector, #178 — e.g. the
two cut faces of a 60° hex wedge, related by `C₆`). That explains BOTH the four
unrelated BC vocabularies AND why `Cn`/`Dnh` are fully realized with **zero**
consumers: #178 is the consumer, blocked on exactly this machinery.

### R8. The phase-space quotient is METHOD-specific, not geometry's

The `S²` in `S²/G⁰` is the **direction** sphere of the transport equation.
Mechanism: the spatial symmetry group's **linear part** acts on directions, so
`ψ(x,Ω) = ψ(gx, Q_gΩ)` and ψ depends on direction only through `S²/G⁰`.

But **phase space is only determined once a method is picked**: diffusion has no
`S²` (a moment closure), CP has integrated angle out. So the quotient is a
concern of methods carrying a discrete measure on `S²`, NOT a universal one.

⟹ The 1-D arm's residual content — *"`SO(2)` acts trivially on a polar
marginal"* — is a statement about a **polar marginal**, which only SN has. It
does **not** migrate to geometry; it stays with the method. (This resolves an
earlier open doubt.)

### R9. It is the MESH, not the raw shape, that computes the point group

Four steelman arguments were put against "geometry computes its own group".
Two fell, two stand — and neither survivor blocks:

| # | argument | verdict |
|---|---|---|
| 1 | a shape has no point set for `_orbit_closure` | **FALLS** — CSG vertices are the intersections of 3 surfaces, a canonical finite set. Caveat: *curved* primitives have no vertices and must declare their group analytically (a cylinder is `D_∞h`); both feed the same machinery |
| 2 | the mesh can BREAK the shape's symmetry | **STANDS** — and becomes a feature, see below |
| 3 | materials break it (a centred vs off-centre pin: `D_4h` vs `C_1`) | **STANDS**, but is *resolvable* — see the centre-of-mass ruling |
| 4 | BCs break it | **FALLS** — geometry already knows the BCs; only the *Law* is unrealised, and that is method-specific. Knowing which faces are reflective/vacuum/periodic is enough |

⟹ The claim sharpens from "geometry computes its point group" to **"a MESH
computes its point group"** — the mesh has the point set, the volumes, AND the
materials. `Mesh1D` (`geometry/mesh.py`) carries `edges + mat_ids + coord +
volumes`: exactly the four things `Sym(domain) ∩ Sym(materials)` needs. Package
placement is unaffected.

The problem's group is `Sym(domain) ∩ Sym(materials) ∩ Sym(BCs) ∩ Sym(source)`.

### R10. ⭐ The ORIGIN is not a convention — it is the CENTRE OF MASS

User's insight, and it closes R1's loose end. **Every symmetry of a body
preserves its mass distribution, hence fixes its centre of mass.** So the centre
of mass lies in the fixed-point set of *every* element of the point group, and
is therefore the canonical seat.

Computable without cross-sections: from material density directly, or from
element/nuclide densities under the per-material homogeneity assumption that
already holds.

⟹ The affine part is not merely *representable* — it is **determinable**. This
is what turns "`origin = 0.0` defeats `Mirror('x')`" from a modelling choice
into a computed fact.

### R11. ⭐ `Sym(geometry) ⊇ Sym(mesh)` is a MESH-QUALITY METRIC

Fell out of R9 unintentionally and is worth keeping. The ideal shape's group is
an **upper bound**; the mesh realises some subgroup of it. **The gap is a
meshing objective** — a good mesh preserves as much of the domain's symmetry as
it can. A future heuristic-meshing criterion, recorded here so it is not lost.

---

## 5. Open questions — NOT blocking, to settle during implementation

1. ⛔ **RESOLVED BY MEASUREMENT, AND THE ANSWER INVERTED THE QUESTION
   (2026-08-03).** I had this recorded as *"`_distinct_azimuths` /
   `candidate_groups` may be quadrature-shaped wearing geometry clothes — a
   square lattice's symmetry is not bounded by its azimuth count."* **False.**
   The bound's reasoning — a `C_n` rotation (`n>1`) fixes no azimuth, so the
   off-axis azimuths fall into FREE orbits of size `n`, hence `n | count` — is
   a fact about ℝ³, not about the sphere. `[M]` sound on every set tried:

   | point set | `n_az` | divisors | true `C_n` |
   |---|---|---|---|
   | 3×3 / 4×4 / 5×5 square lattice | 8 / 12 / 16 | … | 4 / 4 / 4 |
   | 2×4 rectangle | 8 | 1,2,4,8 | 2 |
   | hex 1-ring / 2-ring | 6 / 12 | … | 6 / 6 |
   | cube vertices | 4 | 1,2,4 | 4 |
   | triangle / pentagon | 3 / 5 | … | 3 / 5 |
   | 1-D slab cell centres | 2 | 1,2 | 2 |

   **The real limitation is that the candidate set is AXIS-LOCKED**, which is a
   different and more tractable problem: `Cn`/`Dnh` are realized about **z**
   only and `Mirror` offers **coordinate normals** only. `[M]` A 3×3 lattice
   standing in the x–z plane has a genuine `C_4` about **y**; the z-locked walk
   reads `n_az = 2` and finds nothing. `[M]` A lattice rotated 30° about z
   **keeps** its `C_4` (rotations about z commute with the re-orientation) while
   **both** `σ_x` and `σ_y` fail and only the rotated normal permutes.

   ⟹ That asymmetry — *the rotation family survives re-orientation about its own
   axis, the mirror family does not* — is the campaign's own ruling
   **"containment ≠ subconjugacy: literal for the gate, the aligning ROTATION
   for re-orientation"** landing on a concrete case. The fix is not a polymorphic
   candidate generator; it is to let the walk carry an **orientation** (a
   `RigidMotion` conjugating the realization), which G1's `conjugated_by` already
   supplies. Bears directly on "the mesh computes its own group", since a mesh
   has no obligation to be axis-aligned.
2. **Is `Mesh1D` misplaced?** There is a `transport/mesh/` package (`axis.py`,
   `material_mesh.py`, `material_xs_field.py`) and discretisation arguably
   belongs there, while `Mesh1D` sits in `geometry/`. Noted by the user; a
   SEPARATE carve, not this one. (For symmetry purposes `Mesh1D` is
   well-placed — it has geometry AND materials.)
3. **Naming collision** `Mirror` (the group `{e,σ}` / the surface) vs
   `Reflection` (the element). Lean: `Reflection` = the element in
   `transformation.py`; `Mirror` = the surface/tag; `Mirror`'s realization is
   `[Reflection.through_hyperplane(...)]`.
4. **Does `symmetry.py` move whole, or leave the lattice behind?** Whole-move is
   ~1660 lines + 105 tests, mostly mechanical. "Never as cheap as today" argues
   whole.

---

## 6. Sequencing

| step | what | gate |
|---|---|---|
| **G1** ✅ `6acb6a8a` | mint `geometry/transformation.py`: dimension-generic elements, `RigidMotion`, closure, the point-set orbit primitive; `permutes`/`preserves` | tree green; realizations bit-identical |
| **G2** ✅ | ⭐ **mathematically verify** against PURE MATH — group axioms, `QᵀQ=I`, `det=±1`, involution, order-`n`, Householder/Rodrigues vs independent constructions, conjugation `RσR⁻¹`, orbit–stabiliser | every law mutation-verified |
| **G3** ◐ | route `symmetry.py`'s seven constructors through the core (**DONE**); delete the 1-D/3-D arm split (**BLOCKED — needs a ruling, see §7c**) | realizations bit-identical |
| **G4** ✅ §7d | close the checker-side `roots_of_unity` sites (`_cyclic_ops`, `_vertical_mirrors`) | ⛔ "`Dnh(n_φ)` exact on BOTH sides" is **unachievable** — angle-addition is not a float theorem. Landed criterion: exact at `n_φ ∈ {1,2,4}`, and mirror residual **≡** rotation residual everywhere |
| **G5** | the BC layer's 4 σ_a vocabularies + `_reflect_corner`; `ReflectionLaw` binding | `product(4,5)`'s `ValueError` → `BoundaryError` |
| **G6** ⭐ NEW | **bind the boundary operators' domain/codomain to the trace spaces** (§7f) | a Γ₋→Γ₊ composition is REFUSED, not silently skipped |
| **G7** (was G6) | MoC: replace the 2 guard-free `argmin`s with the certificate; adopt `periodic_trapezoid` | the 8.96e-2 cm link gap becomes visible |
| **G8** (was G7) | ⭐ **migrate the tests** onto the verified machinery | coverage preserved; re-run each gate's justifying mutation |
| **G9** (was G8) | retire superseded spellings + dead code (`_rotation_x/_y`, the `directional.py:209` re-spelling, the ghost-justified epsilons) | grep clean |

**Bit-identity — two sets, two different gates.** *Preservation* (must stay
bit-identical): `_reflections`, `_octahedral_ops`, `directional.py:574`
(`arange(N)[::-1]`), `rules_sphere.py:289`, all `roots_of_unity` output.
*Improvement* (deliberately re-baselined): `_cyclic_ops`, `_vertical_mirrors`,
where `[M]` `_cyclic_ops(4)/(8)` produce **zero** exact zeros against
`roots_of_unity`'s two. Separate gates, or the second silently licenses drift in
the first.

**Deliberately NOT in scope** (filed/tracked, not fixed here): #328 (the SH
`m≠0` defect), the SH polar-axis crosswalk, #178 `SymmetryBoundary`, the
11-fold sweep-reversal cleanup, the dangling `_setup_spherical` xrefs, the
adaptive/exactness-based quadrature selection (see below), `Mesh1D`'s placement.

**Parked with its seam identified:** quadrature selection when the mesh is
irregular and the answer is "no particular affinity". `select_quadrature`
already runs staged gates (stage 0 domain, stage 1 symmetry), so an
affinity/exactness score is a natural **stage 2**, reached only when stage 1
returns nothing. The user's two variants — (a) assume a uniform flux and
pre-estimate exactness, (b) an ML initial-guess framework predicting flux on its
(much smaller) manifold — differ only in what prior seeds the estimate, so one
seam serves both. Monte Carlo is attractive as the estimator precisely because
it needs no quadrature to evaluate the integrand — no circularity. A standing
interest of the user's; revisit deliberately.

---

## 7. G1 — what landed, and what it measured (`6acb6a8a`)

`orpheus/geometry/transformation.py`, ~700 lines, `numpy`-only, exported from
`geometry/__init__`. **Zero production consumers** — G3 wires it. `pyright`
clean on both touched files; `tests/numerics` + `tests/geometry` = 2338 passed,
1 failed (the known `WhiteXminPartial03GL` red, task #33, bisected to
`292a1ba5` — unchanged).

**The surface.** `RigidMotion(linear, translation)` with orthogonality enforced
at construction and only there; `@` composition, `.inverse()`,
`.conjugated_by()`, `.seated_at()`; `on_points` / `on_directions` as two named
actions with **no `__call__`** (a direction has no position — translating one is
the bug the split prevents); `.determinant`, `.fixed_subspace_dimension`,
`.element_order()` (returns the order, not "is it n?" — `gⁿ = e` is satisfied by
every element whose order *divides* n); constructors `identity`,
`translation_by`, `inversion`, `reflection(normal=, offset=)`,
`rotation_from_circle_point`, `rotation`, `rotation_about_axis`,
`signed_permutation`; free function `close_group`.

**R1 and R2 demonstrated, not asserted.** `[M]` A 6-cell production lattice
under a mirror at `x=0` → `None`; the SAME lattice under
`reflection(normal=x, offset=½)` → `[5 4 3 2 1 0]`, exactly the
`arange(n)[::-1]` at `sweep_graph.py:577-580`. `[M]` A bare `(4,1)` GL rule
under `reflection(normal=[1.0])` → `[3 2 1 0]`, identical to the `(μ,0,0)`
3-D embedding — the arm split dissolved.

**Bit-identity, measured against the legacy primitives:**

| | vs | max │Δ│ |
|---|---|---|
| `reflection(normal=e_i)` | `_reflections(axis)` | **0.0** (all three) |
| `signed_permutation` family | `_octahedral_ops()` | **0.0** (48/48) |
| `preserves(...)` permutations | `_orbit_closure(...)` | **0.0** on lebedev(17), LS(8), product D₅ₕ, D₈ₕ |
| planar rotation | `_rotation_z` | 1.11e-16 ⟵ improvement set |
| `rotation_about_axis` | `_rotation_about_axis` | 7.77e-16 ⟵ improvement set |

`[M]` Closure reproduces every order: `C_n`/`D_nh` for n = 2,3,4,6,8; `O_h` = 48
**and the closed set IS `_octahedral_ops()`'s 48 matrices**, not merely a group
of the right size; `I_h` = 120 with 60 proper.

**`rotation_from_circle_point` is the primitive, `rotation(angle=)` the
convenience** — an angle is a lossy parameterisation of a point on S¹.
`[M]` it returns exactly `0.0` where `np.cos(np.pi/2)` = 6.1e-17. **This is the
hook G4 needs** to re-base `_cyclic_ops`/`_vertical_mirrors` on
`roots_of_unity`.

### ⚠ A G3 PRECONDITION discovered here, currently held by accident

`numerics/__init__.py:40` imports `symmetry`, so G3's back-edge
(`symmetry` → `geometry.transformation`) lands inside a live package cycle.
`[M]` It imports clean at every entry point (`orpheus`, `orpheus.numerics`,
`orpheus.numerics.symmetry`, `orpheus.sn.solver`) — **because every one of
geometry's 16 numerics imports is SUBMODULE-level** (`from
orpheus.numerics.X import Y`), never `from orpheus.numerics import Y`. A
partially-initialised package serves the former and fails the latter. Nothing
enforces this today. **G3 must add the gate** alongside the import, in
`tests/test_layer_imports.py`'s family.

---

## 7b. G2 — the gates, and what they measured

`tests/geometry/test_transformation.py` — **42 gates, 96 parametrised cases**,
8.9 s serial under `-O`, module-level `@pytest.mark.foundation`, SymPy rows
in-file and deliberately NOT `slow`-marked (a `slow`-marked foundation gate is
deselected by the canonical `-m "not slow"` sweep and stops guarding). Plan +
closeout: `scratch/g2_verification_plan.md`.

**Mutation battery: 32 mutations, 32 caught, 0 blind.** Twelve redden EXACTLY
one gate, so coverage isolates rather than smears. Main-agent independently
re-ran three and reproduced the agent's counts row-for-row (M21 → 1 failed,
M13 → 13, M01 → 23), with the source restored byte-identically.

**The findings worth keeping:**

- ⭐⭐ **`σ² = e` is Mode-12 designed-green on the ENTIRE affine-offset family.**
  `Qt = (I−2n̂n̂ᵀ)(αn̂) = −αn̂ = −t`, so `g² = (Q², Qt+t) = (I,0)` for *every* α.
  `[M]` `t ∈ {d, 2d, 4d, −2d}·n̂` are **all** involutions while the true mirror
  plane moves by `0.37 / 0.00 / 0.74 / 1.48`. The offset is this campaign's
  motivating defect and the involution law **cannot see it** — its only catcher
  is "the element fixes its named fixed set POINTWISE".
- ⭐ **The seat is a THEOREM, not a convention** — a `G`-preserved weighted point
  set has a `G`-fixed centroid. `[M]` a cube shifted to `(2.5,−1.25,0.75)`:
  **48/48** seated `O_h` elements preserve it, **1/48** unseated. This is what
  converts **R10** from a modelling choice into a computed fact.
- ⭐ **The permutation homomorphism is the only gate that pins three
  conventions at once** (composition order, row-vs-column action,
  `π` vs `π⁻¹`) — none is checkable alone. `[M]` a transposed point action
  reddens it but does NOT redden the direct `g(P) == P[π]` check, because that
  check builds `π` with the same mutated action and a transposed action still
  maps a `G`-invariant set onto itself.
- ⚠ **VACUITY, measured not argued:** on the **abelian** `C_4` the reversed
  composition order agrees **16/16**; on `O_h`, **42/144**. A `C_n` fixture
  would have made the gate completely vacuous — the same shape as the
  stable-sort tie-break gate this campaign already shipped vacuous at
  n ∈ {4,8,12}. Four gates carried this risk; all four now have non-abelian /
  mixed-parity / unequal-weight fixtures.
- ⭐ **NEVER ASSERT TIGHTER THAN THE TYPE'S OWN INVARIANT.** `__post_init__`
  admits `max|QᵀQ−I| ≤ 1e-12`, so a `1e-13` shear is a **legal** element and a
  gate demanding `1e-14` orthogonality of an arbitrary element asserts a
  property the type does not promise. Split: one gate for the type's invariant
  and its rejection threshold, another for the constructors' (far better)
  actual quality — `signed_permutation` is exactly `0`.
- **An absolute `atol` is the wrong contract for a translation residual** — it
  scales `O(ops × ‖t‖ × ε)`. Gates scale by `max(1, ‖desired‖_∞)`; worst over
  3000–4000 draws per law is **6.9e-15**.
- **`on_points − on_directions == t` is NOT a float theorem** (2.2e-16) while
  **`on_points == on_directions + t` is bit-exact 6000/6000.** Gated in the
  true — and stronger — direction.
- **The harness lied before the code did, again.** The battery's first run
  reported 32/32 BLIND while the summaries plainly read `23 failed` / `63
  failed`: the parser wanted `FAILED` lines that `-q --tb=no` never emits.
  **The positive control is what exposed it** — a mutation making `reflection`
  return `+I` cannot leave 42 gates green. Never run a mutation battery without
  one.

**Laws with no pure-math reference** (recorded so none is smuggled in wearing
one): the `D_nh` mirror-plane placement is a literature *setting*, not a theorem
(`[M]` orthogonality, det, closure and group order are all preserved by rotating
the mirror set, so no intrinsic law can distinguish them); the `(n̂,d)`/`(−n̂,−d)`
gauge; the numeric VALUE of the match window (lower-bounded by measurement,
upper-bounded only by the consumer's own minimum point separation — so default
it *relative to* that intrinsic quantity, which also retires
`_NODE_WINDOW_FACTOR` honestly); and the physical centroid computation, which is
G8's, not G2's (the test-migration step, renumbered when G6 was inserted).

---

## 7c. G3 — routed; and the arm split turned out to hide a CONVENTION CONFLICT

**Done.** All seven constructors now build `RigidMotion`s; `_close_group`
delegates to the core (the hand-rolled I_h BFS with its own dedup and its own
cap is **retired**); `_orbit_closure` is now a ~10-line delegation to
`preserves` and keeps only the *measure*-level question; `OrbitCertificate`
carries `RigidMotion` + `Permutation` and is off its `np.eye(3)` lock;
`_finite_contains` compares elements, not matrices. The one production consumer
outside the file — `directional.py`'s hand-built `np.eye(3)` with a flipped
sign, a **fourth** re-spelling of a reflection in a module that already imported
the checker — is retargeted.

**Gate: `tests/numerics` 1762 passed, 0 failed.** That includes
`test_lebedev_is_EXACTLY_octahedral` and
`test_level_symmetric_is_EXACTLY_octahedral`, which assert `landing == 0.0` for
all **48** `O_h` operators — the strongest available bit-identity evidence.

`[M]` Preservation set exactly `0.0`: `_reflections` (all three),
`_octahedral_ops` (48, in order), `_inversion_op`. Improvement set:
`_rotation_z` `1.1e-16`, `_rotation_about_axis` `1.1e-15`, `_vertical_mirrors`
`4.4e-16` — the last because the core normalises the normal and the old code
did not. `[M]` The raw normal was *already* unit (0 or `1.1e-16` off), so this
is a pure ULP re-rounding, neither better nor worse (`max|QᵀQ−I|` is `4.4e-16`
both ways, marginally worse at n=5). Exactness arrives in **G4** — but only at
`n_φ ∈ {1,2,4}`, and G4 measured that no construction can do better at the rest
(angle-addition is not a float theorem). See §7d.

Two tests migrated with the retirement (a retirement includes its tests):
`_oh_exactness`'s `nodes @ np.asarray(g).T` → `g.on_points(nodes)`, and
`np.linalg.det(M)` → `M.determinant`.

### The G3 precondition is now GATED, and the hazard is REAL

Two new gates in `tests/test_layer_imports.py`: a structural one (an input-layer
package may not import the `numerics` **package**, only its submodules) and an
end-to-end one (six entry points must import in a **fresh interpreter** —
in-process is worthless, everything is already in `sys.modules`).

`[M]` **Mutation-verified.** Inserting `from orpheus.numerics import
SubgroupOfO3` at true module level in `geometry/mesh.py` (AST-anchored, so it
parses) gives exactly the predicted
`ImportError: cannot import name 'SubgroupOfO3' from partially initialized
module 'orpheus.numerics'`. The structural gate reds on `[geometry]`; the
end-to-end gate reds on `orpheus.numerics`, `orpheus.numerics.symmetry` and
`orpheus.sn.solver` — while `orpheus.geometry` stays **green**, correctly, since
importing geometry directly never triggers the cycle.

### The arm split — RULED (user, 2026-08-03): make it honest about dimension

**Done.** `_check_invariance_1d` and `_is_reflection_invariant_1d` are
**deleted** (the latter was one of the five duplicate permutation engines the
audit found — a second matching loop *inside the same file*, with a different
window, that discarded π to return a bool). `_check_invariance_3d` is renamed
`_invariance_on_points`; there is now ONE arm, fed by `_embedded_nodes`, which
lifts the DATA (`μ ↦ (μ,0,0)`, `(x,y) ↦ (x,y,0)`) rather than the group —
`O_h` and `I_h` genuinely are 3-D and there is nothing to restrict them to.

| | before | after |
|---|---|---|
| `Sym(GL(8))` | **`O3`** | the three mirrors |
| asymmetric: `SO2`, `Cn(4)` | **True** ⚠ | **False** ✓ |
| live path `Mirror('x')` | True / False | **True / False — preserved** |

**Two answers to the user's questions, both measured.**

*Why did a 1-D measure answer `O(3)`?* `[M]` `GL(8)` was certified invariant
under **15/15** candidate groups — FALSE for nothing. `_check_invariance_1d` had
exactly ONE discriminating branch (the `μ→−μ` test) and waved the rest through.
Its output was a **one-bit function of the input** wearing a nineteen-group
lattice walk: pass the one test and every group reads True, so the maximal
element is necessarily the TOP of the lattice. It could not have answered
anything else.

*What else over-promises?* `[M]` **Nothing — it was confined to the 1-D arm.**
Monotonicity (`A ⊆ B ∧ P(B) ⟹ P(A)`) over 9 measures × 15 groups: **0
violations** (the 2026-08-02 run found 68; Q5.0.x closed them). Generator-set vs
full closed group: **0 disagreements** — unlike ERR-072's sampling, these
generating sets genuinely generate. Continuous groups on 3-D rules: **all
False**, exact criteria properly applied.

⭐ **AND THE "0 VIOLATIONS" CARRIES ITS OWN LESSON.** Monotonicity read clean
even for the 1-D case — because when *everything* reads True the implication is
vacuously satisfied. **A consistency law is blind to uniform
over-certification.** Anti-pattern #15 catches *inconsistency*, never *uniform
falsehood*; only comparing against a COMPUTED answer catches that.

⭐⭐ **The root cause was a FIELD DOING TWO JOBS.** `gauss_legendre_on_mu`
declared `invariance_group=SubgroupOfO3.SO2`, defended in-comment as *"the tag
names the group that was INTEGRATED OUT to produce the μ-marginal"*. True and
useful — and a **different fact from the one the field's name promises**. Its
own defence gave the tell: *"SO(2) acts trivially on μ, so EVERY measure on
[-1,1] satisfies it"* — **a tag satisfied by every possible value carries no
information and cannot be wrong, so it could never have been checked.** Now
`Mirror("x")`. The spent half is not lost: it already lived in
`AngularSymmetry.continuous_isotropy`, where it derives the support.

**A gate whose witness was DISSOLVED by the fix.**
`test_selector_asks_the_nodes_not_the_declared_tag` pinned "a lattice query
against the declared tag would reject GL for a slab" — a trap that existed
*only because the declaration was wrong*. Fixing the tag made its assertion
false. **A gate whose witness can be removed by fixing an unrelated bug is
pinning a coincidence, not a mechanism.** Re-posed to inject a
true-but-not-maximal declaration (`Trivial`) and show stage 1 is unmoved — which
cannot be dissolved, because the injected tag IS true.

Also re-posed: `test_gauss_legendre_1d_so2_invariant` → its **inversion**, and
`test_so3_on_a_polar_marginal_requires_reflection_symmetry`, whose own premise
("`SO(2)` and `C_n` leave the polar cosine alone") was the `(0,0,μ)` reading.

### ⛔ THE ARM SPLIT WAS NOT A REFACTOR — IT HID A CONVENTION CONFLICT

`_check_invariance_1d` uses **two incompatible embeddings of the polar cosine**,
and deleting the split would silently pick one:

* the **mirror** arm is derived from `(μ, 0, 0)` — μ on **x**, so `σ_x` is the
  one real test (this is the LIVE selection path, and it is correct);
* the **rotational** arm asserts "`C_n`/`SO(2)` rotate about z, which does not
  move the polar cosine" — true only if μ is on **z**. Under `(μ,0,0)` a
  rotation about z *does* move the node.

`[M]` measured on `gauss_legendre_on_mu(8)`:

| | shipped 1-D arm | the `(μ,0,0)` embedding |
|---|---|---|
| `σ_x`, `σ_y`, `σ_z` | True | True ✓ |
| **`C_4` about z** | **True** | **False** |
| `Sym(GL(8))` | **`O3`** | would be `{e, σ_x}` |

And `[M]` on a deliberately **asymmetric** 4-point measure: `Mirror('x')` reads
`False` (correct) while **`SO2` and `Cn(4)` read `True`**, and the walk reports
`Sym = [SO2, Mirror('z'), Mirror('y')]`. That is **ERR-072's shape** — a
continuous group certified without being tested — and it disagrees with the 3-D
arm's own exact criterion (`_is_axis_supported`), which would say `False`.

**Not currently a live defect:** `GEOMETRY_ANGULAR_SYMMETRY`'s slab/sphere rows
are `Mirror("x")`, which is correct on both symmetric and asymmetric input. The
unsound answers sit on non-consumed rows.

**Why it needs a RULING, not a refactor.** R8 already holds that *"`SO(2)` acts
trivially on a polar marginal"* is a statement about a **polar marginal**, which
only SN has — so it does not migrate to geometry. But that means deleting the
arm changes what `is_invariant` MEANS for a 1-D measure, from "the 3-D rule this
marginal came from is G-invariant" (unknowable here) to "this point set in ℝ¹ is
carried onto itself" (checkable, and `Sym = {e, σ_x}`). The second is the honest
one and the dimension-generic core makes it free — but it changes many answers
and is a semantics decision, not a mechanical one. **This is the same family as
the SH polar-axis crosswalk (six sites say z, the seventh says x) and #328.**

---

## 7d. G4 — the checker re-based, and the acceptance line CORRECTED

**Done.** `_cyclic_ops` builds `C_n` from `roots_of_unity(arange(n), n)` as exact
circle points through `rotation_from_circle_point` — the hook G1 minted for
exactly this. `_vertical_mirrors` is now the **coset** `C_n·σ₀` (σ₀ = the
`φ = 0` mirror, normal `ê_y`), i.e. the group fact `D_n = C_n ⊔ C_nσ` spelled
directly. `_rotation_z` had exactly one consumer and **retires** with it.

### ⛔ The acceptance line as written was NOT achievable, and the reason is structural

The header said "`Dnh(n_φ)` **exact on BOTH sides**". Carrying the node at
azimuth `2πj/n` to the node at `2π(j+k)/n` requires
`cos a cos b − sin a sin b = cos(a+b)` to hold in IEEE-754, and the
angle-addition formula **is not a floating-point theorem**. So exactness exists
only where every root involved is an axis point. `[M]` sweeping `n_φ = 1…24`:
the landing residual is exactly `0.0` at **n_φ ∈ {1, 2, 4}** and `1.1e-16` to
`3.3e-16` everywhere else — and no construction can do better at the rest.

This is *not* a shortfall to be fixed in a later phase; it is a proof that the
line asked for something impossible. The honest criterion, and what landed:
**exact where exactness exists, and the mirror half no longer worse than the
rotation half anywhere.**

### What it measured

| | before | after |
|---|---|---|
| `C_4` matrix entries exactly `0`/`±1` | 30/36 | **36/36** |
| `σ_v(4)` entries exactly `0`/`±1` | 26/36 | **36/36** |
| `C_8` / `σ_v(8)` exactly representable | 54/72 · 47/72 | **68/72 · 68/72** |
| landing on `product(4,4)` (rot · mir) | `2.1e-16` · `4.7e-16` | **`0.0` · `0.0`** |
| landing on `product(8,16)` (rot · mir) | `7.1e-16` · `1.3e-15` | `2.2e-16` · **`2.2e-16`** |

⭐ **The mirror residual is now EXACTLY EQUAL to the rotation residual at every
one of `n_φ = 1…16`** — because σ₀ is a coordinate mirror, hence a bit-exact
signed diagonal, so composing with it is a column sign flip that adds no
round-off. Previously the mirrors were 1.7–4.8× *worse* than the rotations (`[M]`
over the six rules swept), because
each was built from a normal at the half-angle `kπ/n + π/2` — a root of unity of
order **4n**, a different generator from the rule's own — and then normalised.
That equality is the falsifiable half of the change and is gated as such.

### Verdict-neutrality, proven by VALUE not by proxy

The re-base only *tightens* the landing, and a tighter landing can only make a
match easier — so the hazard is one-directional and real: a rule that read
`False` could start reading `True`. `[M]` **3234 verdicts** (77 shipped rules ×
42 groups, product/level-symmetric/Lebedev), old construction vs new:
**0 changed** (767 True / 2467 False both ways). Green tests are the proxy; this
is the value comparison (vv-principles anti-pattern #12).

### Gates — `tests/numerics/test_symmetry_exactness.py`, 30 cases

A separate file from `test_symmetry.py` on purpose: that one gates *group
structure*, this one gates the *floating-point realization*. Orthogonality,
determinant, closure and group order all survive a less-accurate construction of
the same group, so nothing in the structural file can see the operator set drift
into a wider tolerance — and the drift would be absorbed by the very node-match
window the checker uses to decide invariance.

`[M]` **Mutation battery, 4 mutations + control, all caught**, and the isolation
is clean:

| mutation | reds | reads |
|---|---|---|
| control (none) | **0/30** | the harness does not red spuriously |
| M1 revert `_cyclic_ops` to `np.cos(angle)` | 9 | the rotation-side gates |
| M2 revert `_vertical_mirrors` to the half-angle normal | 29 | every mirror gate; the rotation-only control correctly stays green |
| M3 `_landing_residual` returns `0.0` always | **exactly 1** | the instrument-liveness control, and nothing else |
| M4 mirror planes rotated by π/2 (the historical convention bug) | 11 | residual blows to `1.97e-01` |

⭐ **M4 reds on ODD `n_φ` only** (1,3,5,7,9,11,13,15) — a direct measurement of
the claim the `_vertical_mirrors` docstring has carried since G3, that a π/2
plane rotation maps the normal set onto itself for even `n` and is therefore
invisible there. The designated convention gate
(`test_vertical_mirror_planes_follow_the_dnh_setting`) also reds, so the division
of labour holds: that gate owns the convention, these own the accuracy.

### Gate

`tests/numerics/test_symmetry.py` + `roots_of_unity` + `registry` +
`rules_product`: **433 passed, 0 failed** (5:16). `tests/numerics` +
`tests/geometry`: see the commit body. `pyright` CLI clean on both touched
files; `check_docstring_xrefs` **0 dead targets across 813 files**.

### ⛔ G4.5 PROPOSED AND ABANDONED — the index-permutation route to a true `0.0`

The obvious follow-on: since the residual is a floating-point artefact of the
matvec, answer "does `g` permute the nodes?" by **integer index arithmetic**
(`j ↦ (j+k) mod n`) instead of geometrically — which `[M]` reproduces the node
set bit-identically, `0.0` at every `n`. Proposed as a G4.5, **measured, and
abandoned.** Four independent findings, any one of which is disqualifying:

| | measured | verdict |
|---|---|---|
| **margin** | window `1e-11` vs residual `1.1e-16`–`4.1e-16` | 24 000×–88 000×; buys **no** correctness |
| **growth** | `n_φ` 8→256 moves the residual only `1.1e-16`→`4.1e-16` | sub-linear; never approaches the window |
| **cost** | `Dnh(64).is_invariant` = 154 ms / 256 nodes | not a bottleneck at production sizes |
| **coverage** | `level_symmetric(4/8/16)` and `lebedev(11/17)` have **NO** azimuthal index structure | product rules only |

⭐ **The performance case and the exactness case point at DIFFERENT rules.** The
slow path is the brute-force candidate scan on a 110-node Lebedev — and Lebedev
is exactly what the index route cannot serve. An optimisation that misses the
only expensive case is not an optimisation.

**The architectural kill, and it is the decisive one.** Because the route serves
only product rules, the geometric path **must stay** — so this can never retire
its predecessor and is structural debt by construction (`.claude/rules/
coding-standards.md`, "retire as you go"). It would be a permanent SECOND
implementation of "does `g` permute this node set", in the module whose own
durable finding is that `_orbit_closure` is *the ONE bridge* matrix→permutation,
and where this campaign already retired **five** clones of that very job.

**The codebase had already ruled**, in `permutes`'s own body:

> *"one match rule cannot disagree with itself the way a fast path plus a
> fallback can."*

That comment defends the single distance-matrix rule against precisely this
shape. G4.5's version is strictly worse than the one it rejects: the two paths
would answer by **different mathematics**, so a disagreement is a silent wrong
verdict rather than a slow one.

Also required, and not free: `LevelStructure.azimuth` is `float64` — the integer
numerator is not carried, so the route needs `(numerator, denominator)` threaded
from `periodic_trapezoid` through `product_mu_phi`, i.e. a change to the
generating measure while Q5 is in flight on that surface.

**What survives is the insight, not the code:** the residual lives in the
MATVEC, not the operator, so no factorisation of the operator reaches it — `[M]`
Cartan–Dieudonné (rotation = two reflections) is 2–4× WORSE at every `n`,
because each reflection is `I − 2n̂n̂ᵀ` plus a normalisation and then you multiply
two of them. And `{1,2,4}` is exactly the **signed-permutation** set: a rotation
whose entries are all `0`/`±1` is one of `I`, `−I`, `±90°`, and only there does
the action become a relabel with no arithmetic. That is a closed form for the
empirical exact set, not a coincidence.

(`[M]` a related non-finding: `rotation_from_circle_point` computes
`Q[0,0] = 1 + (c−1)`, which differs from the raw root in 32 entries — and is
sometimes BETTER, landing on exactly `−0.5` at `n = 3` where the raw root is
`−0.49999999999999994`. Writing entries directly is worse at `n = 3/5/6`, same
at `8/16`, better at `12`. Pure ULP noise; leave it alone.)

### Refuted candidates — recorded with their structural reason

Two trig sites survive in `symmetry.py` and **both are correctly out of scope**,
measured rather than assumed:

- **`_icosahedral_ops` buys NOTHING from the re-base.** Feeding
  `roots_of_unity(1, 5)` instead of `np.cos(2π/5)` produces a **bit-identical**
  matrix (`[M]` max|Δ| = `0.0`; 0/9 entries exactly representable either way).
  The icosahedral axes are golden-ratio irrationals, so the AXIS sets the
  exactness ceiling and the angle cannot move it — `[M]` 34/1080 entries exact
  over the closed group, landing `1.2e-15` on the 12 vertices. Do not re-attack
  this; it would also require minting a `rotation_about_axis_from_circle_point`
  sibling for no measured gain.
- **`_distinct_azimuths` runs point→angle**, not angle→point: `arctan2` on the
  off-axis nodes to bound which cyclic families are candidates. Roots of unity
  go the other direction and do not apply.

### Retired

`_rotation_z` (one consumer, replaced). Its docstring's point — that the
rotation happens in the `(x,y)` **plane** and "about an axis" is a 3-D-only
convenience — is already carried by `rotation_about_axis` in the core, which is
where it belongs. `_Z_AXIS` is now unused, but `[M]` it was **already unused at
HEAD** — pre-existing, and left alone rather than folded into this change.

---

## 7e. G5 — the self-paired deck half, and a refuted premise in the plan's own G5 line

**Scope, ruled by the user (2026-08-04):** do **Specular + Identity first, to try
out the design**. Periodic and rotational are deferred because they need
**surface pairs**: by Poincaré, a face may be paired with ITSELF (an involution,
domain face = codomain face — the reflection) or with a DISTINCT face (the
codomain of one is the domain of the other — periodic, rotational). Only the
self-paired half is nameable by one face, so only it is buildable today.

The code already carried the distinction in Thurston's language: `SpecularMirror`
documented a reflecting face as a quotient **with fixed points** (an orbifold
reflector boundary), explicitly contrasted with `SpatialWrap`, *"whose action is
free"*. Fixed points ⟺ self-paired; free ⟺ genuinely paired.

### ⛔ The plan's own G5 line was WRONG about the four σ_a vocabularies

This plan said "FOUR type vocabularies for σ_a that do not know about each
other", implying they should merge. A re-audit refuted it: they are **four
different categories** — `SpecularMirror` (deck `Γ₊→Γ₋`), `SpecularReemission`
(constitutive kernel `Γ₋→Γ₋`), `SpecularReturn` (a curried factory), and
`symmetry.Mirror` (a subgroup tag). Three of the four are a deliberate,
documented, sweep-schedule-load-bearing split, and merging them re-opens exactly
the conflation **B3.0** removed — *unobservably*, because a rank-one `R`
annihilates `G`. The user's scope is the correct one and touches `R` not at all.

(Also refuted: my re-audit brief claimed B3.4a/b/c postdated the boundary audit.
Backwards — they landed 2026-08-01 and the audit *cites them by name*. Only G3
and G4 postdate it, and both touched only `numerics/symmetry.py`.)

### Stage 1 — `SelfPairedDeck` minted and verified (`3cce383a`)

Guard: **`is_linear ∧ dim Fix ≥ d − 1`**. ⛔ My first design guarded on
`element_order() in (1, 2)` — the CONVERSE of the insight, and wrong. A
proactive `test-architect` dispatch caught it before any code: `E(3)` has FOUR
involution families, and the half-turn (`order 2, det +1, fix 1`) and inversion
(`order 2, det −1, fix 0`) map a face to its **opposite** — precisely
`SpatialWrap`'s deferred job. An order-only guard admits the elements the type
exists to exclude.

The shipped guard is strictly stronger, carries **no tolerance**, and is a
**theorem-carrier**: a linear orthogonal `Q` fixing a hyperplane pointwise is
`I` or a reflection, so `Q² = I` FOLLOWS. Involution became a derived property
to assert, not a premise to trust. Both clauses are load-bearing — a glide
passes the fixed-set clause and fails linearity; an inversion is the reverse.

⭐ **The linearity clause closes a class no gate could ever close.**
`on_directions` DROPS the translation, so a mirror plane at the wrong POSITION
is bit-identical in the permutation, the realized image, the snapshots and any
scalar — Mode-12 designed-green at every tolerance. Unspellable is the only
closure that exists. (`[M]` identical at `offset ∈ {0, 2.5, −17.0}`.)

`[M]` 23 gates; 7 mutations + control, all caught. M4 first produced **no
output**, and reading that rather than the summary found a real defect: two
`parametrize` lists built the SUT at COLLECTION time, so a guard regression
became a collection ERROR that took the whole file down. Made lazy.

### Stage 2 — migration, and the two gates that would have decayed silently

Six production sites moved (`white`/`vacuum`/`albedo`/`zero_flux`/
`prescribed_inflow` → `SelfPairedDeck.identity()`, `reflective` →
`.mirror(axis=…)`); both retired classes deleted from `_factors.py` and the
package surface.

⭐ **`realizer.py:726` passes `axis=law.axis`, not `geometry_map.axis`** — the
realized permutation comes from the LAW, so the G slot is purely declarative and
realization is untouched. That is why the migration is contained to
`geometry/boundary/` + tests, and why bit-identity of the permutation is
structural rather than something to gate.

**`test-architect`'s highest-value finding: two gates decay to TAUTOLOGIES on
the collapse, silently green.** `[M]` measured after the fact —
`SelfPairedDeck` appears on BOTH the permuting side (`reflective-a1`,
`reflective-partial`) and the non-permuting side (vacuum, white, albedo,
prescribed_inflow, zero_flux), so the sets are **not disjoint** and a type-level
assertion cannot separate them. Both re-posed in the SAME commit as the
collapse:

- `test_every_production_law_states_both_factors` — the spec table's *type*
  column became the expected **value**, compared by equality. Strictly stronger:
  the type column could not distinguish `reflective-a1` (axis x) from
  `reflective-partial` (axis **y**) — both were just `SpecularMirror` — and the
  value column pins the axis.
- `test_specular_mirror_is_the_only_ordinate_permuting_geometry` — re-posed onto
  the **answer set** (which law IDs permute) plus the mechanism one level down
  (a self-paired element permutes angle iff `det = −1`). Immune to any future
  re-typing, which the type-level form was not.

Two gates got *stronger* for free: `test_factors_are_frozen` had a vacuous
branch for the field-less `IdentityMap` and now runs for every geometry factor.

---

## ⏸ COMPACTION POINT #1 — before G6 (HISTORICAL; superseded by #2 at the end of §7h)

Everything above is LANDED and committed; everything below is SCOPED and not
started. A session picking up cold re-anchors from **this file + `git log`**,
never from a conversation summary.

- **HEAD at the point:** `a6878a82` on `refactor/operator-strategy-layers`.
  Verify every hash in the header table with
  `git merge-base --is-ancestor <hash> HEAD` — all nine were ancestors when this
  marker was written, and the header is a snapshot that can rot.
- **Tree clean.** No uncommitted production, plan, docs or agent-memory state.
- **Known reds — 4, none of them this campaign's, and that is MEASURED not
  assumed** (a `git worktree` at HEAD and again at `95fa693c` reproduced the
  three SN ones with BYTE-IDENTICAL signatures):
  `TestWhiteXminPartial03GLSnapshot::test_matches_the_frozen_scaled_lambertian`
  (~1 ULP, task #33); `test_cart2d_1g_vacuum_apply_principled_equiv` (1152 ULP);
  `test_cart2d_2g_specular_apply_principled_equiv` (296 ULP);
  `TestBitIdenticalCurvilinear::test_spherical_inward_bit_identical`
  (`assert False`). The last three are the quadrature campaign's deliberate set.
- **Gate costs, measured:** `tests/geometry` ≈ 10 s · `tests/geometry +
  tests/sn/operators + tests/sn/sweep/core` ≈ 100 s · full tree at
  `-m "not slow"` ≈ 52 min. Canonical invocation `python -O -m pytest`, SERIAL.
- **Static gates:** `tools/check_docstring_xrefs.py orpheus tests docs` must
  read `DEAD TARGETS : 0`; `sphinx -E -W` (use `-E`, never `rm -rf docs/_build`)
  builds with 0 warnings; `npx pyright` is the trusted checker (the streamed
  LSP diagnostics are the #226 artifact — advisory only).
- ⛔ **Never `git checkout <path>` / `git restore` / `git stash` / `git clean`
  on a tracked path** — this tree carries uncommitted-by-policy state elsewhere
  and a git-level discard is irrecoverable (lesson L28). To compare against
  another commit, use `git worktree add` (as the red-set triage above did).

---

## 7f. G6 (SCOPED, not started) — every operator knows its spaces; the SPACE owns shape and traversal

### The principle (user, 2026-08-04) — this is the scope, the boundary case is only its first consumer

**Tracked tree-wide as #330**; G6 below is its BOUNDARY-TIER slice. The general
statement (every operator, mandatory binding, the space owning traversal) lives
in the issue so it survives this campaign; #295 (`LayoutBearingSpace`) is a
concrete sub-piece of its traversal half.

> Every operator, bulk or boundary, **must know its domain, codomain and
> space**. Without that binding the adjoint has to be hand-rolled instead of
> falling out for free from mathematical well-posedness. Whatever **shape**
> something has is an algorithmic / book-keeping detail, **not math**. Whoever
> imposes a shape — *the* thing that determines the shape of everything
> downstream of it — must contain the abstraction used to **traverse, iterate
> and operate** on everything inheriting that shape. Any shape change happens
> at that one place, and everyone downstream gets it automatically **without a
> single change to the mathematical operation**. Book-keeping is not math; it is
> a detail necessary to correctly apply mathematical operations.

**Origin.** A review question after G5 (*"does the realizer produce a
LinearOperator with bound domain/codomain on the trace space, with adjoint?"*),
whose measured answer generalised into the statement above. G6 was first scoped
as "bind the BOUNDARY operators"; that was the symptom. This is the scope.

### ⭐ Half of this is ALREADY BUILT, and built correctly

`_AdjointOperator.apply` (`numerics/operator.py:1204-1227`) computes the Hilbert
adjoint **entirely from the spaces**:

```python
z      = inner_codomain.apply_metric(y)        if inner_codomain is not None else y
result = self.inner.apply_transpose(z)
result = inner_domain.apply_inverse_metric(result) if inner_domain is not None else result
```

i.e. `A† = G_V⁻¹ Aᵀ G_W`, with the module's own comment: *"The space owns the
metric; the adjoint wrapper is metric-representation-agnostic."* `FunctionSpace`
owns `inner_product` / `norm` / `apply_metric` / `apply_inverse_metric`, plus
the space algebra (`__mul__` → `TensorProductSpace`, `dual()`). **That IS the
principle, already implemented.** G6 is a COMPLETION, not an invention.

### ⛔ The defect: the binding is OPTIONAL, so the derivation SILENTLY NO-OPS

`domain` / `codomain` are `Optional[FunctionSpace]` defaulting to `None`, and
when `None` **both metric applications are skipped** — `.H` degrades to the bare
Euclidean transpose with no error and no warning.

⛔ **The count below was wrong and G6.0 corrected it.** This line read *"28 of 54
override; **26 inherit the `None` default**"* until 2026-08-04. `54` and `28` are
right *if* you count only `cls.__dict__` overrides — but **8** of the remaining
26 inherit a genuine, working binding from a non-`LinearOperator` ancestor
(`InverseWrapMixin` ×5, `OperatorSum` ×2, `OperatorProduct` ×1) and are NOT
degraded. `[M]` the true ledger is **54 total / 36 bound / 18 unbound**, of which
**6 are abstract** ⟹ **12 concrete unbound classes.** A static override census
over-states the surface; see §7g for the runtime census that settles it.

The bulk tier (fission, scattering, leakage, multiplication) binds; the boundary
and utility tier (`PermutationOperator`, `TraceRestrictionOperator`) does not.

⚠ **The degradation is PER-END, and half-bound is worse than either extreme.** An
operator binding `codomain` but not `domain` gets `G_W` applied and `G_V⁻¹`
skipped — neither the Hilbert adjoint nor the Euclidean transpose. `[M]` two
shipped classes are half-bound at runtime **100 % of the time** and are hot:
`RadialCharacteristicReconstruction` (3089 applies) and `WindowedSweep` (1000).

That `None` default is [[lessons-L19]] exactly — *a default that encodes an
unstated invariant*. It silently means "assume Euclidean", and nothing in the
type system says so.

**And it is a live CORRECTNESS hazard, not hygiene.** `[M]` the SN trace metric
`|Ω·n|·w` restricted to `Γ₋` is **non-constant** — max/min = `1.35`
(gauss_legendre 4), `3.47` (product 4,4), `5.6` (lebedev 17) — so it genuinely
matters. The specular permutation **preserves** it (ERR-042's
weight-preservation), so `[G, Aᵀ] = 0` and dropping the metric is invisible
*there*. The boundary adjoint is therefore **right by accident of the one case
anyone checked** — the Mode-12 commutator criterion. Any boundary operator that
does NOT preserve the trace weight gets a **silently wrong** `.H` today.

⛔ **Two corrections G6.0 forced here.**

1. **No operator is silently wrong TODAY.** `[M]` `.H` is read on exactly **7**
   AST-level expressions in all of `orpheus/` (2 internal to `operator.py`, 5 in
   `sn/solver.py`'s adjoint), every one on a **bound composite**; and
   `SNMesh.BOUNDARY_OPERATOR_REGISTRY` admits only `{reflective, vacuum}`, both
   weight-preserving. So G6 **closes a trap before its consumer arrives** rather
   than fixing a live defect. The architectural argument is untouched — realization
   IS binding to a space — but the urgency claim was overstated, and the honest
   version is the one that goes in the commit message.
2. **The Ni/Ti supermirror `R` of §R7 DOES NOT EXIST in the tree.** Do not plan
   against it. `[M]` the *only* non-weight-preserving boundary operator shipped is
   the Lambertian `AngularAverageOperator` (`orpheus/sn/boundary/angular.py:41`),
   reached by the white / albedo+`IsotropicReturn` arm — measured
   `‖A†−Aᵀ‖/‖Aᵀ‖` = **0.209 / 0.684 / 0.612** (GL4 / product / lebedev), with the
   weighted law holding at `2.2e-16` and the metric-dropped answer wrong by
   `5.5e-1`. Four of five shipped laws are Mode-12 blind (vacuum → `ZeroOperator`;
   reflective and albedo-specular → weight-preserving permutation, `G₋ == G₊[local]`
   **bit-identical** on all three quadratures; periodic → opposite normals;
   prescribed-inflow → rank-0).

### The half that is NOT built: shape and traversal ownership

`FunctionSpace` carries `shape` as a **passive tuple** and has no iteration /
indexing / traversal abstraction. So book-keeping leaks into the math:

- `PermutationOperator(perm, axis=0)` — an **axis** inside a mathematical object.
- `gamma_out.to_local(perm[inflow])` (`realizer.py`) — the local↔global index
  map computed at the **call site**, not owned by the space.
- `TraceRestrictionOperator` carries `n_restricted`, a **length**, where a
  codomain space belongs.

⟹ a layout change today touches every operator. That is the precise thing the
principle forbids.

### What is measured true today

| question | measured |
|---|---|
| realized reflective BC is a `LinearOperator` | ✔ `TensorProductOperator`, `isinstance` True |
| `.H` exists | ✔ returns `_AdjointOperator` |
| `domain` / `codomain` | ✘ **`None`** on `gauss_legendre(4)`, `product(4,4)`, `lebedev(17)` |
| a trace `FunctionSpace` exists to bind to | ✔ `AngularTraceSpace(FunctionSpace)`, `spaces/angular_trace_space.py:299` |
| the half-traces Γ₊(f) / Γ₋(f) exist as SPACES | ✘ they are **index arrays** + a `TraceRestrictionOperator` |
| `TraceRestrictionOperator` binds spaces | ✘ it carries `n_restricted`, a **length** |

**So the gap is layer-wide and consistent, not an inconsistency**: the whole
boundary operator family is *shape*-typed, never *space*-typed. And the cost is
named in `LinearOperator.domain`'s own docstring — *"when either operand has
`None` … the composability check is **SKIPPED**"*. B3.2's narrowing to
`G : Γ₊ → Γ₋` is real in the guards and the prose and **absent from the type
system**, so composing a boundary law the wrong way round is not refused.

### ⛔ The design constraint, and it is MEASURED not hypothetical

`FunctionSpace.__eq__` is `(name, shape)` — the weakest of the nesting tiers
([[lessons-L29]]). Measured, for every shipped quadrature:

```
gauss_legendre(4)  |Γ₊(xmin)|=2   |Γ₊(xmax)|=2   same size=True  same INDICES=False
product(4,4)       |Γ₊(xmin)|=4   |Γ₊(xmax)|=4   same size=True  same INDICES=False
lebedev(17)        |Γ₊(xmin)|=49  |Γ₊(xmax)|=49  same size=True  same INDICES=False
```

⟹ **a half-trace space whose `name` does not encode the FACE (and the sign)
compares EQUAL to its opposite face's**, so a cross-face composition would
type-check while being wrong — the exact class G6 exists to refuse, re-admitted
by the mechanism meant to close it. The name is load-bearing; `shape` alone is
never sufficient here.

### Substeps — in dependency order, each independently landable

| step | what | note |
|---|---|---|
⛔ **RE-SCOPED 2026-08-04 after G6.0 (user ruling).** G6 does **the boundary tier,
complete**. The tree-wide **mandate moves to #330** — because G6.0 measured that
retiring the `None` default at the `LinearOperator` base breaks
`orpheus/homogeneous/solver.py` (a production path) until the **energy/bulk** tier
gets real spaces, which is a different campaign's work. That is a *principled*
boundary, not a risk-dodge: the mandate is blocked on work G6 does not own.
Decay cost of the two scopes, `[M]`: **~20** test sites boundary-scoped vs **~150**
tree-wide.

| step | what | note |
|---|---|---|
| **G6.0** | ✅ **DONE** — survey + triage of the unbound operators, and the runtime binding census | `scratch/g6_0_operator_binding_survey.md`; findings folded into §7g |
| **G6.1** | mint the **three-tier Γ ladder**: `Γ(f)` (per-face slot), `Γ₊(f)`, `Γ₋(f)`, each with face **and** sign in the `name` | `[M]` the plan previously scoped only the halves — the **middle tier has no type either**, and is NOT recoverable as `Γ₊ ⊕ Γ₋` (see §7g) |
| **G6.2** | carry the **restricted metric** onto each (`partial_current_metric` restricted to the tier) | so a half-trace pairing is PHYSICAL, not Euclidean — the ERR-067 family |
| **G6.3** | bind the boundary tier: `γ± : Γ(f) → Γ±(f)`, `G : Γ₊(f) → Γ₋(f)`, `R : Γ₋(f) → Γ₋(f)`, and periodic's cross-face `G : Γ₊(f) → Γ₋(f_opp)` | the wrap is where a face-blind name silently passes. **Bind from the REALIZER, per law — never from the class** (§7g) |
| **G6.3b** | **absorb `AngularAverageOperator.apply_transpose`** (the step its own docstring assigns to boundary phase B5) | user ruling: without it the flagship gate can only run on a stand-in, which is the proxy-evidence pattern `vv-principles` #12 forbids |
| **G6.5** | move **traversal** onto the trace spaces: they own iteration, the axis, and the local↔global index map; boundary operators stop carrying `axis=` and the realizer stops doing `to_local` | the second half of the principle, **boundary-scoped**; the tree-wide sweep is #330 |
| ~~G6.4~~ | ~~make the binding MANDATORY~~ | **MOVED TO #330** — see the re-scoping note above |

### ⭐ The acceptance test that actually falsifies the principle

The principle's own words are testable, and this is the gate that decides
whether G6 succeeded rather than merely compiled:

> **Change a layout in ONE place — the space — and assert (a) ZERO diff in any
> operator file, and (b) BIT-IDENTICAL numerical output.**

⛔ **The concrete instruction here was IMPOSSIBLE as written, and G6.0 measured
the ceiling.** This paragraph said *"permute the ordinate storage order (or move
the group axis) … every gate must stay green with `np.array_equal`."* `[M]`
floating-point summation is not associative, so a 49-term reduction is
bit-identical under reordering only about **25 %** of the time — the criterion
demanded something arithmetic cannot deliver. Same failure class as G4's
acceptance line (§7d): *compute the ceiling before implementing.*

**Honest replacement, and it is a STRONGER discriminator, not a weaker one:**
permute the **face packing order** in the layout. `[M]` bit-identical **4/4 per
face** while the flat buffer genuinely moves — so it tests exactly the claim
(book-keeping changed, mathematics did not) without smuggling in a
reduction-order change that no implementation could survive. Every operator body
must be untouched and every gate green under `np.array_equal`. If any operator
needs an edit, book-keeping is still living in the math and G6.5 is not done.

Supporting acceptance:

- A deliberately reversed composition (`Γ₋→Γ₊` fed where `Γ₊→Γ₋` is required)
  **raises**, naming both spaces. Mutation-verified, with a positive control.
- `Γ₊(xmin) != Γ₊(xmax)` — the measured collision below, now refused.
- **An unbound operator cannot be constructed** (G6.4). The mutation is
  re-introducing the `None` default: it must red.
- ⭐ **A metric-sensitive adjoint gate**: build a boundary operator that does NOT
  preserve the trace weight, and assert `⟨Ax, y⟩_W = ⟨x, A†y⟩_V` in the WEIGHTED
  pairing. `[M]` this gate is IMPOSSIBLE to write today (the metric is dropped),
  and it is the one that would have caught ERR-067. Its control leg: the
  specular case, where the commutator vanishes and weighted ≡ Euclidean.
- Bit-identical everywhere binding is added: spaces change no numbers.
- The pre-existing red set is unchanged (task #33 + the 3 quadrature-campaign
  SN reds, signatures 1152 ULP / 296 ULP / `assert False`).

### Risks, in the order they are likely to bite

1. **Turning on a check that has never run.** Every `None` today means a skipped
   composability test; binding spaces makes those tests fire tree-wide for the
   first time. Some will be REAL bugs and some will be false positives from
   legitimately-untyped legacy paths. **Survey before binding** — count how many
   compositions currently skip, and triage — rather than discovering it as a
   wall of reds.
2. **The periodic wrap crosses faces** (`domain_face → face_opposite`), so it is
   the one law whose domain and codomain belong to DIFFERENT faces. It is the
   natural first witness that the face-encoded name works, and the natural first
   casualty if it does not.
3. **A half-trace space is per-(face, quadrature)**, so naive construction in a
   hot path could allocate per call. Build it on the trace space (which is
   already cached) rather than at each realization.

### Explicitly OUT of scope, with the reason

**Do NOT give `G` a self-adjointness property.** It is honestly false for the
narrowed operator: Γ₊ and Γ₋ are DISJOINT index sets (measured on all three
quadratures), and self-adjointness requires domain = codomain. The local-
coordinate matrix *looks* symmetric on `gauss_legendre(4)` and `product(4,4)`
and is **NOT** on `lebedev(17)` (`M ≠ Mᵀ`, `M² ≠ I`) — so the apparent symmetry
is an artefact of local index ordering and is quadrature-dependent. The
involution belongs to the **un-narrowed full-trace reflection**, where it is a
theorem; narrowing is exactly what destroys it as stated.

⚠ **One live documentation falsehood to fix while here.**
`numerics/operator.py:2188-2191` (`PermutationOperator`) states *"SN specular
reflection through `reflection_index` is an involution"* and exposes
`is_involution` "for downstream consumers that benefit from knowing
self-adjointness". True of the full-trace permutation, **false of what the
realizer produces** — and it sits beside a flag inviting precisely the wrong
inference. The boundary audit flagged the symptom (`is_involution` reads `False`
on `lebedev(17)`) as "zero production consumers, not a live bug"; it is not a
code bug, it IS a prose bug. Its own caveat — self-adjointness *"in the
**unweighted** inner product"* — is the second half: the trace pairing is
weighted by `|Ω·n|·w`, so the physical-metric claim is separate and needs
ERR-042's weight-preservation.

---

## 7g. G6.0 — the survey, and the two rulings that re-scoped G6 (2026-08-04)

Deliverables: `scratch/g6_0_operator_binding_survey.md` (inventory + triage +
runtime census) and `scratch/g6_verification_plan.md` (7 gates, 8-mutation
battery + 2 controls). Design measurements: `scratch/g6_design_measurements.md`.
Every `[M]` below is a probe, not an inference. The §7f claims these falsify
have been corrected **in place** above (present-tense-false is the bug).

### ⭐ R1 — the space ladder is THREE tiers; the middle one was never scoped

| tier | shape | space today |
|---|---|---|
| `Γ` whole boundary | `(n_faces·n_ordinates,)` | ✔ `AngularTraceSpace` |
| **`Γ(f)` per-face slot** | `(n_ordinates,)` | ✘ **MISSING — and unscoped until now** |
| `Γ₊(f)` / `Γ₋(f)` | `(|Γ±|,)` | ✘ missing |

`[M]` `_face_restrictions` builds `TraceRestrictionOperator(..., n_total=n_ordinates)`
— the **per-face** count, not `layout.total_size`. And the middle tier is not
recoverable as `Γ₊ ⊕ Γ₋`, because the tangential class is non-empty:

| quadrature | `\|Γ₊\|` | `\|Γ₋\|` | sum | `n_total` | **tangential** |
|---|---|---|---|---|---|
| `gauss_legendre(4)` | 2 | 2 | 4 | 4 | **0** |
| `product(4,4)` | 4 | 4 | 8 | 16 | **8** |
| `lebedev(17)` | 49 | 49 | 98 | 110 | **12** |

⚠ **`gauss_legendre(4)` has ZERO tangential ordinates** — a gate written only on
it is blind to the entire tier (`vv-principles` Mode 7: the ansatz nulls the term
it was meant to exercise). Every partition-touching G6 gate MUST include
`product(4,4)` or `lebedev(17)`.

Corollary: the face-collision constraint applies **one tier higher** than §7f
says. `Γ(xmin)` and `Γ(xmax)` share `shape=(n_ordinates,)`, so with `__eq__` on
`(name, shape)` the face must be in the name for the middle tier too.

### ⭐ R2 — G and R are different KINDS of thing, and only R can be gated numerically

User ruling, 2026-08-04, confirmed by the code's own Protocol
(`geometry/boundary/_factors.py:399`): **the Lambertian is `R`, not `G`.** A
Lambertian reflective wall is `G = identity` with `R = Lambertian`. The
distinction is *geometry is a conceptually enforceable quotient; response is an
assumption about a real surface* — which is exactly why `_factors.py` files
`SpecularReemission` under **R**: *"a polished wall's specular return, which is
constitutive rather than geometric because a wall is not a quotient."*

| factor | typing | what a gate can see |
|---|---|---|
| **`G : Γ₊ → Γ₋`** | cross-space, NOT an endomorphism | `[M]` **metric-blind in every shipped case** — no numerical gate can ever see it. Its whole value is the **type-level refusal** of a reversed composition. |
| **`R : Γ₋ → Γ₋`** | **endomorphism** — domain == codomain | where the metric bites (Lambertian, ‖A†−Aᵀ‖/‖Aᵀ‖ = 0.209/0.684/0.612). The flagship gate lives here, and needs only `Γ₋(f)`. |

Two consequences worth carrying: the flagship gate is **more reachable** than
scoped (one space, not the Γ₊/Γ₋ pair); and because `R` is an endomorphism it
**can** be self-adjoint — a free, sharp invariant to gate — while §7f's
"OUT of scope: do not give `G` a self-adjointness property" stays exactly right,
since `Γ₊` and `Γ₋` are disjoint.

### ⭐ R3 — bind from the REALIZER, per law; NEVER from the class

`[M]` `TensorProductOperator` and `ScaledOperator` each land in "right by
accident" **or** "would be silently wrong" **depending on which law is realized
into them**. A class-level binding would therefore be a lie for half its
instances. This is the single most important implementation constraint in G6.3.

Related, and it constrains ordering: `[M]` `OperatorSum.domain` is a
**first-non-`None` fallback**, so **both legs of a sum must be bound in one
change** — otherwise the composite advertises a space it never verified.

### R4 — the composability check's blast radius, as a number

`[M]` **178 skipping compositions** over `tests/geometry + tests/sn/operators`
(`OperatorSum` 172/1694 = 9.2 %; `OperatorProduct` 6/1568 = 0.4 %), in 10
distinct operand shapes. Concentrated in **test fixtures, not production**:
`Sum(IsotropicScattering, IsotropicN2N)` alone is 148 of them, and production
`solve_sn_adjoint` skips **0 of 595**. Highest-risk shape is the 5×
`Product(BulkAnalysisOperator, SweepOperator)` — **left codomain unbound only**,
i.e. half-bound, where a genuine mismatch can hide because the bound side already
committed to a space and nothing ever compared them. **Bind, run, then triage —
and if it raises, record it as an ERR-NNN candidate before "fixing" it by
widening the space.**

### R5 — the decay list under the RULED scope is empty where it matters

The verification plan's decay list (12 strict xfails + ~8 gates) was written for
the **tree-wide mandate**. Under the ruled boundary scope, `[M]`:

- The 12 `strict=True` xfails in `tests/sn/architecture/test_monomorphic_leaves.py`
  are about the **SN bulk ladder** (L/C/S/F/B) and the meshless homogeneous path
  — they do **NOT** flip. (They are #330's todo list, and its acceptance gates
  are already committed there — the `xfail(strict=True)`-set-IS-the-todo-list
  technique from the operator-strategy campaign.)
- `[M]` **zero** gates in `tests/geometry`, `tests/sn/boundary`, `tests/sn/operators`
  assert the boundary tier is unbound (`domain is None` / `space_anonymous`) —
  so the **DECAY class, the dangerous one that stays green while discriminating
  nothing, is EMPTY here.**
- What remains is the **BREAK** class, and it is visible as reds: any
  `assert_array_equal` on a `.H` of a newly-bound operator breaks at **1–2 nulp**,
  because `G⁻¹(f·(G x))` vs `f·x` is not bit-identical. Re-derive each tolerance
  from the metric round-trip depth — **do NOT relax to `allclose`**; the
  tolerance is the claim.

### R6 — the positive control, MEASURED (not predicted)

`apply_inverse_metric := identity` — linear and shape-preserving, so it is clean
under `vv-principles` anti-#18 (it does not break a structural law and inflate
the red count). `[M]` **3 failed → 33 failed = 30 net reds**, all in the
metric-aware `.H` family, and **zero boundary-law operators reddened** — because
`domain is None` short-circuits the call. That is *empirical proof of the
blindness* G6 exists to remove, and it is the control every G6 mutation run must
carry (anti-#17).

The test-architect also **reproduced anti-#17's exact failure mode** in passing:
`grep "^FAILED"` returned zero rows on a run whose own summary read `33 failed`
(ANSI codes). The instrument lies before the code does.

### Refuted / corrected, with the structural reason

- ⛔ **"Mint a generic (`EuclideanSpace` / anonymous) space first"** — proposed by
  both agents as Tier 0, **REJECTED by the user, and the code agrees.**
  `homogeneous/solver.py:136-141` hand-passes `basis_shape=(ng,1)` *"because the
  meshless operators carry no `domain` space to derive it from"* — it **has** an
  `EnergyGrid` and **needs** an energy space. So that `None` is the defect wearing
  the sentinel's other meaning, not a mathematical category. A generic space would
  have **institutionalized the ambiguity under a new name**. Correct fix: every
  production `None` is a defect and gets its real space (#330); test fixtures get
  an explicit small named space (`FunctionSpace("R3", (3,))`), one line each.
- ⛔ **The Ni/Ti supermirror `R` (§R7)** — does not exist in the tree. Removed as
  a motivating case above.
- ⛔ **"26 unbound operators"** — a static-override artifact; 8 inherit real
  bindings. True surface is 12 concrete classes.
- ⛔ **"No operator is half-bound"** — the explorer wrote this first and corrected
  itself by runtime census. Two classes are half-bound 100 % of the time, and
  half-bound is worse than either extreme.

---

## 7h. G6.3 — DESIGN OF RECORD: materialize G and R, then bind (user ruling, 2026-08-04)

**G6.1 LANDED** (`34f465cc`): `AngularFaceTraceSpace` + the three accessors
`trace.face_space/outflow_space/inflow_space(face)`, cached, metric-carrying,
face+role in the name. 56/56 gates, mutation battery with a live positive
control, 2294 numerics+transport passing, wide gate at the 4 known reds.
It also fixed a **pre-existing base-class defect** (see §7h.0 below).

### ⛔ G6.3's ORIGINAL SCOPE RESTED ON A FALSE PREMISE — measured

§7f scoped G6.3 as "bind `γ±`, `G`, `R`, and periodic's cross-face `G`". `[M]`
**`G` and `R` are not operators at all.** `law.geometry_map` returns a
`SelfPairedDeck` / `SpatialWrap` and `law.response_kernel` a
`SpecularReemission` / `LambertianReemission` / `ScalarResponse` — **descriptors**,
with no `.apply`, no `.domain`, and `isinstance(_, LinearOperator) is False`.
`SNBoundaryRealizer.realize` emits **ONE** operator per law with `R∘G` already
collapsed. So two of the four scoped rows have **zero construction sites**, and
`R : Γ₋→Γ₋` is un-bindable — nothing in the SN realizer is an endomorphism of
`Γ₋`.

**The sharpest tell** (`[R]` `realizer.py:770-778` vs `:349`): `albedo` +
`SpecularReemission` calls **the same `_specular_kernel` helper `reflective`
uses**, funnelling into the same `_attenuated_kernel_operator` body — while
declaring the mirror lives in **`G`** for one and **`R`** for the other. One
construction site would need two different codomains.

⟹ **the G/R split is precise at the DECLARATION tier and collapsed at the
REALIZATION tier.** `_factors.py:399` states it exactly (*"a polished wall's
specular return … is constitutive rather than geometric because a wall is not a
quotient"*), and the realizer then welds the two factors into a single un-named
product. That is the *welded un-named operation* the standing
build-the-machinery ruling names as **a failure to realize the algebra**.

### ⚠ And binding the collapsed operator would be INERT, not safe

`[M]` the binding was installed as a pytest plugin and run: **4941 bindings,
~5100 tests, ZERO new failures** (the pre-existing reds identical to baseline),
16 periodic cross-face bindings firing. Read that as **inert**: `_reflect_trace`
composes with three raw `.apply` calls (`sn/operators/boundary.py:527`, `:570`),
never `@` or `+`, so **no composability check exists on the path G6.3 types**.
Two hand-rolled size checks already stand in for the absent gate
(`realizer.py:385`, `:502`). Binding alone ships honest metadata that gates
nothing.

### ⭐ THE RULING (user, 2026-08-04): materialize, then bind

> Make `G` and `R` real `LinearOperator`s so the law is spelled **`R @ G`**,
> bind each to its space, and route the composition through the operator
> algebra so the check actually fires.

Payoff beyond tidiness: **`R` is an endomorphism `Γ₋→Γ₋` and CAN be
self-adjoint; `G : Γ₊→Γ₋` cannot** (disjoint index sets). Today neither claim is
expressible — which is precisely why `PermutationOperator`'s `is_involution`
docstring is already flagged present-tense-false for the narrowed operator
(§7f). Materializing makes one claim true and the other unspellable.

### ⛔⛔ CORRECTION 2026-08-04 — the `R @ G` common-intermediate design above is ALSO WRONG

The paragraphs that stood here took the two Protocols at face value and derived a
uniform `R∘G` split with `Γ₋` as the intermediate. **A review of the math with the
user refuted it, and the tree's own realization was the witness all along.** This
was the SECOND false premise in G6.3; the lesson is recorded in §7h.1.

`[R]` `AngularAverageOperator` (`sn/boundary/angular.py`) — the realization of
`LambertianReemission`, i.e. of a **response** — types itself
***"Lambertian re-emission, `Γ₊ → Γ₋`"*** and its body is

```python
psi_avg = (cos_w * psi).sum(axis=0) / self._norm     # collapses the ANGLE axis
return broadcast_to(psi_avg[None, ...], (n_inflow,) + psi.shape[1:])
```

`[M]` **No realized response is an endomorphism of `Γ₋`.** `LambertianReemission`
→ `Γ₊→Γ₋`; `SpecularReemission` → `Γ₊→Γ₋`; `ScalarResponse` → a commuting
`ScaledOperator`. So `R : Γ₋ → Γ₋` is not merely un-bindable — it does not exist,
and *the response carries the crossing* for every law where the response is the
non-trivial factor.

⟹ **TWO Protocol docstrings are present-tense-FALSE and must be corrected**
(the code is right; the declaration tier is wrong):

| `geometry/boundary/_factors.py` | false claim | truth |
|---|---|---|
| `BoundaryResponseKernel` | "`R : \Gamma_- \to \Gamma_-`" | `R : Γ₊ → Γ₋` |
| `BoundaryGeometryMap` | "the crossing … is part of **this factor**, not of the response" | whichever factor is NON-TRIVIAL carries the crossing |

### ⭐ The CORRECTED algebra (user, 2026-08-04) — the factorization is PER LAW KIND

Not one uniform split. The campaign theorem *"exactly one of `G`/`R` is
non-trivial"* sharpens into: **a law is either a deck transformation (atomic) or
a response (composable), and the non-trivial factor is the one that crosses.**

| law kind | structure | intermediate | adjoint |
|---|---|---|---|
| **deck transformation** | ⭐ **ATOMIC** — a measure-preserving bijection does not factor into two meaningful pieces. `Γ₊(f) → Γ₋(f)` self-paired, or `→ Γ₋(f_partner)` through paired faces. **Pure geometry: a law imposed by theorem, transport-method-agnostic** — it needs the method only to know WHICH SPACE to act on. | **none** | **guaranteed, a THEOREM**: the composition operator of a bijection `g` has `G⁻¹ = G_{g⁻¹}`, and measure-preservation makes that inverse the transpose. Purely geometrical. |
| **response** | **`N` COMPOSED operations**, `Γ₊(f) → … → Γ₋(f)`. Lambertian = an **outflow angle contraction** `C` then an **isotropic broadcast** `B`. Constitutive: an assumption about a real surface. | ⭐ **`S(f)`** — the angle-integrated per-face scalar current. Already exists as a TYPE: `ScalarTraceSpace` ("per-face `(J⁺, J⁻)` partial-current pairs… already angle-integrated"). It is `psi_avg` in the body above — **a real physical quantity, currently an anonymous local**. | **conditional** — see the theorem below |

### ⭐⭐ SUBSUMED BY THE CHAIN VIEW (user, 2026-08-04) — one structure, not two

The two-row table above is **correct but over-articulated**: it makes "atomic"
and "composed" two *kinds*, hence two code arms. The user's reframing collapses
them:

> A boundary law is a **sequence from outflow to inflow**. The first link's
> domain is :math:`\Gamma_+(f)`, the last link's codomain is
> :math:`\Gamma_-(f)`, and the interior is whatever the physics needs. A deck
> transformation is simply a sequence of **length 1**.

⟹ **"atomic" is not a kind, it is a degenerate length.** One code path.

**And the user's diagnosis of WHY the endomorphism trouble arose is correct:** it
was an artifact of forcing a fixed **two**-factor decomposition with
**pre-declared** types. In a chain the interior types are *determined by the
chain*, never declared. Checked against every shipped law:

| law | chain | endomorphism? |
|---|---|---|
| specular | `Γ₊ → Γ₋` | none |
| periodic | `Γ₊(f') → Γ₋(f)` — ⛔ **ASPIRATIONAL, not descriptive.** `[M]` 2026-08-07 its link is a bare `IdentityOperator`, the ONLY one of six that is not a typed arrow between the half-traces. Step 7 makes the row true | none |
| white / Lambertian | `Γ₊ → S(f) → Γ₋` | none |
| albedo + specular | `Γ₊ → Γ₋ → Γ₋` (or `Γ₊ → Γ₊ → Γ₋`) | exactly where one BELONGS — the scalar amplitude |
| vacuum | the zero chain | none |

An endomorphism appears where one naturally occurs and is **never required**.
The `R : Γ₋ → Γ₋` declaration was the whole problem, self-inflicted.

### ⭐ The chain adjoint TELESCOPES at any length — `[M]`

Probe: a `7→3→5→4` chain (3 links, **2** interior spaces), interior metrics
scaled across 15 orders of magnitude.

```
max|chained - direct|                       = 3.553e-15
<Ax,y>_V3 - <x,A*y>_V0                      = 8.882e-16
interior metrics × 1e-8 and × 4e7           = 7.105e-15
```

⟹ `(A₁…Aₙ)* = Aₙ*…A₁* = G₀⁻¹(A₁…Aₙ)ᵀGₙ` — **every** interior metric cancels,
not just one. So the theorem in the next subsection is the **chain** theorem, and
the 2-link response is its special case. **No interior space needs a declared
type; it needs only a non-degenerate metric.** Length 1 is the same formula with
no special arm.

### ⭐ NO NEW TYPE — the machinery was already right, the USAGE was wrong

One refinement to the chain: a physically real wall can be *partly* polished,
`α_spec·P + α_diff·L` — two channels in **parallel**. So the structure is
compositions **and sums** — which is exactly the operator algebra ORPHEUS
already has. `OperatorProduct` and `OperatorSum` exist; `@` already checks
composability at every link.

⟹ **Do NOT mint a `BoundaryChain` / tensor-network type.** It would be a new
abstraction where the algebra suffices, against defer-until-≥2-instances. Bind
the endpoints to `Γ₊(f)` / `Γ₋(f)` and **the chain constraint enforces itself**,
link by link. A boundary law is an **expression in the operator algebra**, not a
bespoke two-factor struct.

⚠ **Two structures, easy to conflate — both already live in the code:**

| spelling | meaning | example |
|---|---|---|
| `@` | **sequential** composition — the chain | `B @ C` (contract, then broadcast) |
| `&` | **tensor product across axes** | `k_omega & IdentityOperator()` — angle ⊗ group |

A boundary law is a **`@`-chain whose links may themselves be `&`-products**.
Different directions; do not merge them.

### ⭐ THE THEOREM that settles whether a response adjoint is well-defined

`[M]` probe `$CLAUDE_JOB_DIR/tmp/g6_response_adjoint.py`. For `R = B∘C`:

```
R* = C*B* = (G₊⁻¹ Cᵀ G_S)(G_S⁻¹ Bᵀ G₋) = G₊⁻¹ Cᵀ Bᵀ G₋ = G₊⁻¹ Rᵀ G₋
```

**The intermediate metric CANCELS.** Measured over `G_S` spanning ELEVEN orders
of magnitude — Euclidean `0.0`, `1e-6` → `1.11e-16`, `3.7e5` → `1.11e-16` — with
the weighted adjoint law `⟨Rx,y⟩_{G₋} = ⟨x,R*y⟩_{G₊}` holding at **exactly 0.0**.
At `G_S = 0` the cancellation **BREAKS** (error `7.6e-1`).

⟹ the requirement on `S(f)` is exactly **ONE binary condition: a NON-DEGENERATE
metric.** Not "the physically correct metric".

Two consequences to carry:

1. ⚠ **Only the COMPOSITE is metric-free.** `C*` and `B*` *individually* depend
   on `G_S`, and factoring exists precisely so the factors become usable — so
   `S(f)`'s metric must still be chosen DELIBERATELY. It just cannot make the
   composite wrong.
2. ⭐ **This retroactively gives a G6.1 gate its real reason.**
   `test_the_half_trace_metric_is_strictly_positive` reads as hygiene; it is in
   fact **the precondition for the factored adjoint to be exact**. The `Γ(f)`
   full tier has a DEGENERATE metric (zeros on the tangential rows), so it can
   never serve as an intermediate; the halves can, because they exclude them.

### ⭐ The payoff: the missing transpose falls out for FREE

`AngularAverageOperator` reports `is_adjointable=False` and defers its transpose
to boundary phase **B5**. Factored, there is nothing to defer:

* `Cᵀ(s) = cos_w · s / norm` — the outer product
* `Bᵀ(φ) = Σ_{Γ₋} φ` — the sum over inflow
* ⟹ `Rᵀ(φ) = (cos_w / norm) · Σφ`

That is the campaign's whole thesis in one line: **the adjoint falls out of
well-posedness instead of being hand-rolled.** G6.3b's "absorb the B5 step"
stops being scope creep and becomes a consequence.

⛔ **`|Γ₊| == |Γ₋|` MUST be a GUARD wherever a local-index pairing is used.**
`_narrowed_zero_operator` (`realizer.py:248-256`) already says why:
*"`|Γ₊| == |Γ₋|` on every reachable fixture … an accident, not a contract."*
Measured 2/2, 4/4, 49/49. Note the corrected design needs this in FEWER places —
a factored response never pairs by index at all (it goes through `S(f)`), so only
the deck-transformation arm carries the constraint, and there it is intrinsic to
the bijection rather than assumed.

### Per-law realization map `[M]` — the ONE operator today, typed `Γ₊(f) → Γ₋(f)`

| law | site | class today | cross-face |
|---|---|---|---|
| vacuum | `realizer.py:274-290` | `ZeroOperator` | no |
| reflective α=1 / 0<α<1 / α=0 | `:349` / `:352` / `:345` | `TensorProductOperator(Perm@:421, Id)` / `ScaledOperator(TP)` / `ZeroOperator` | no |
| white (3 regimes) | `:349` / `:352` / `:345` | `TP(AngularAverage@:559, Id)` / Scaled / Zero | no |
| albedo+Specular | `:770-778` | **identical to reflective** | no |
| albedo+Isotropic | `:780-789` | **identical to white** | no |
| periodic | `:859` | `TensorProductOperator(Id, Id)` | **YES** |
| prescribed | `:885-888` | `IncomingSourceOperator` | no |

Three different *outermost* classes for one law (α=1 / 0<α<1 / α=0) ⟹ **bind at
`_attenuated_kernel_operator`'s single exit, not per-arm** (B3.4b already
funnelled four laws through it).

### Constructor readiness `[M]`

- **Nothing to change (2):** `ScaledOperator`, `_BoundBoundaryOperator` — pure
  passthrough properties. This dissolves §7g R3's worry for `ScaledOperator`: it
  *forwards*, so the per-law decision is made once, on its operand.
- **Additive constructor change (6):** `TraceRestrictionOperator`,
  `PermutationOperator`, `TensorProductOperator`, `ZeroOperator`,
  `AngularAverageOperator`, `IncomingSourceOperator` — each has its own
  `__init__`, no `__slots__` blocker, `domain`/`codomain` plain overridable
  properties.
- **No `__init__` at all (1):** `IdentityOperator` — **do not touch** (11
  production + 61 test sites).
- `[M]` **nothing computes a codomain from its domain**, so explicit binding
  conflicts with nothing. Only forwarding (`Scaled`, `_Bound`), projection
  (`OperatorProduct`), and `OperatorSum`'s first-non-`None` fallback.

### Order of work

⛔ **The order below is REVISED for the corrected algebra.** The superseded
version sequenced "materialize `R : Γ₋→Γ₋`" as step 2; no such operator exists.

| # | step | why here |
|---|---|---|
| **0** | ✅ **DONE** `ab7c3711` — correct the false Protocol docstrings | ⭐ **THREE claims were live, not two.** The unanticipated one: the Lambertian's *"IS self-adjoint under the cosine-weighted inner product"*, asserted in FOUR sites while the same docstring 8 lines up says the operator *"maps BETWEEN two half-traces"*. True of the pre-B3.4a full-face endomorphism (`[M]` `‖R*−R‖ = 2.8e-17`); **type-incoherent** after the narrowing — the same defect class as `PermutationOperator.is_involution`, so it is a PATTERN in this tier, not a slip. 15 candidate sites, **9 corrected, 3 deliberately left** (they scope the crossing claim to deck maps, where it is true). Also gave the adjoint theorem a permanent home (`tests/numerics/test_factored_adjoint_identity.py`) instead of leaving it in a probe |
| **1** | ✅ **DONE** — bind the space-side `γ±` (`angular_trace_space._face_restrictions`) | `γ± : Γ(f) → Γ±(f)` in the TYPE system. ⭐ **The check now FIRES**: wrong-half AND cross-face compositions raise `IncompatibleOperatorComposition` naming both spaces (shape could never catch cross-face — `\|Γ₊(xmin)\| == \|Γ₊(xmax)\|` everywhere). Constructor cross-checks the spaces against `n_total`/`len(indices)`, because a mis-binding is invisible at apply-time. Binding stays OPTIONAL here — the realizer's own producer is unbound until step 3, and a mandatory-binding gate would be false today (a test pins that transitional state). ⚠ **`γ±`'s `.H` EQUALS its transpose by construction** — `Γ₊(f)`'s metric IS `Γ(f)`'s restricted, so the weights cancel; the value is the type-level refusal, not numerics. Needs its negative leg (a ×3 metric shifts by O(1) relative vs ≤2 ulp) or the gate would pin "the metric is ignored". `[M]` the round-trip `(g·y)/g` costs **1** nulp on tangential-bearing quadratures, **0** on `gauss_legendre(4)`; tolerance set to 2 from the round-trip DEPTH, not fitted |
| **2** | mint `S(f)`, the angle-integrated scalar current space, on `ScalarTraceSpace`, with a **non-degenerate** metric | the theorem's one requirement; also names `psi_avg`, today an anonymous local |
| **3** | factor the Lambertian into `C : Γ₊(f) → S(f)` and `B : S(f) → Γ₋(f)`, each bound; the law is `B @ C` | the transpose falls out here (`Rᵀ = CᵀBᵀ`), closing the deferred B5 step |
| **4** | ✅ **DONE — but NOT as scoped; the premise was refuted, for the FOURTH time in this step** | see §7h.2 |
| **5** | ✅ **DONE — and it refuted step 8, the FIFTH refutation in this step** | see §7h.3. The deck arm is bound `Γ₊(f) → Γ₋(f)`, `γ₊` with it, `is_involution` RETIRED; `TensorProductOperator` turns out to drop the binding, so **step 8 gains a substep** |
| **6** | ✅ **DONE — landed inside the affine campaign's P3, `8d552395`** | ⭐ vacuum is NOT "the zero chain" as a distinct structure; it is a **length-1 chain whose single link is the zero morphism** (user, 2026-08-05) — same shape as step 5's deck arm, different kernel. `ZeroOperator` gained `domain`/`codomain`, bound in `_narrowed_zero_operator` via `checked_space_extent` against the half-trace index arrays, and gated by the NEW `tests/numerics/test_zero_operator_spaces.py` (which also inherited the `\|Γ₊\| ≠ \|Γ₋\|` Mode-12 discrimination from a retired operator — no face fixture can make it, since the two half-traces are equal-sized on every reachable face). Everything else about prescribed inflow lives in **`.claude/plans/affine_boundary_source_channel.md`** — read its ⏸ COMPACTION POINT #2. ⚠ **Do not re-plan vacuum/prescribed here** — that plan is the authority |
| **7** | ✅ **LANDED 2026-08-07 (realization phase): the deck transformation is UPLIFTED once, and periodic is the translational case** | `e13313a8` + `c3bb7341` + `635c22e5` — `PairedDeck` minted (SpatialWrap retired); the arrow derived from the MOTION via `Quadrature.ordinate_permutation` → `RigidMotion.linear_part.preserves` (⛔ the plan's `RigidMotion.permutes` spelling was the AFFINE action — ordinates are directions; see §7h.6's landing note); `_assert_wrap_identification` dissolved INTO the kernel (its partner return is now consumed, closing the old threading claim). ▶ **REMAINING: §7d — the `reflection_index` retirement audit** (the user's gate; `[M]` 33 live refs, sweep pole seed included) |
| **8.0** | ✅ **LANDED `72f5ce97` — and it was NOT the plumbing step this row said it was** | see §7h.4. The anchor flipped `XPASS(strict)` exactly as the flip-proof predicted. ⛔ **This row's framing was wrong**: it scoped 8.0 as an enabler for 8, so a session trusting it would have shipped a **live `.H` repair with zero gates** — `[M]` binding the product moved the realized *white* law's Hilbert adjoint by **87 % relative** (specular: **0.0**, its metric cancels) |
| **8** | ✅ **LANDED `84d7861f` — the check now runs on the production path** | see §7h.5. `face_action = law @ γ₊`, ONE operator serving BOTH legs (`(law∘γ₊)ᵀ = γ₊ᵀ∘lawᵀ`). `[M]` bit-identical on all four shipped law kinds; flip-proof measured (reverting reddens the B3.4c re-injection gate, and **only** that one of 8) |

### §7h.6 — step 7 RE-SCOPED: uplift the deck transformation, once

> ✅ **LANDED 2026-08-07 in three commits** — the realization phase is
> COMPLETE; only the retirement substep (§7d below, appended at landing)
> remains. 7a `e13313a8` (PairedDeck retires SpatialWrap; the guards are
> complements of the fixed-subspace criterion; `RigidMotion.is_translation`);
> 7b `c3bb7341` (`_deck_kernel` — one body from the motion for
> reflective/albedo-specular/periodic; `Quadrature.ordinate_permutation`;
> `_specular_kernel` + `_assert_wrap_identification` retired onto it;
> periodic's link is the bound `PermutationOperator(arange)`
> Γ₊(f')→Γ₋(f); both transitional pins flipped strictly ON LANDING and
> were re-scoped); 7c `635c22e5` (142 gate rows; the π⁻¹ convention
> pinned on `product(4,8)` — NOT `product(4,4)`, whose local permutation
> is degenerately `arange`; the periodic metric row). `[M]` wide slice
> `tests/{numerics,sn,transport,geometry,diffusion} -m "not slow"` =
> **7 failed / 6159 passed**, 17:30 — the SAME pre-existing seven as
> ⏸ #4's baseline (2× cart2d principled_equiv · 3× affine-carve sha256
> (#333) · 1× spherical diamond · 1× WhiteXminPartial03GL (#33)); the
> passing count grew 6017 → 6159 with the new gates. Complement exposure
> covered by an import smoke (the only structural delta reaching
> cp/moc/mc/derivations is directional.py's new module-level imports —
> all packages import clean).
>
> ⭐ **The landing's own surprise — ordinates are DIRECTIONS.** `[M]`
> `RigidMotion.permutes/preserves` match through `on_points` (the AFFINE
> action), so the unit wrap translation matched to **None** — a false
> "not a symmetry" — where the deck truth is the identity permutation.
> Found twice independently (main-agent probe; test-architect's battery,
> which pins it on the OBSERVABLE: motions sharing a linear part yield
> bit-identical kernels, 26/78 rows red under the affine spelling). Fixed
> by minting `RigidMotion.linear_part` (the E(d)→O(d) projection) and
> matching through it. The plan's own sketch (`permutes() → Permutation`
> below) carried the conflation — a proposed means, not a measured one,
> and the labels held.
>
> ⚠ **Two behavior deltas, deliberate, no pins broken:** periodic no
> longer demands the method space's OWN-face `outflow_indices` (its
> domain is the PARTNER's Γ₊, derived from the face name); and a specular
> law on a rule without mirror closure (odd-n_φ product) now raises the
> kernel's `BoundaryError` at realization instead of `reflection_index`'s
> `ValueError` lookup miss (`BoundaryError` subclasses `ValueError`, so
> type-level pins held; the `test_quadrature_directional.py:294` message
> pin is on the TABLE api, untouched until §7d).

**USER RULING, 2026-08-07** (verbatim): *"for the periodic (translation rigid
motion) and rotational symmetry (rotational rigid motion), we need a
`PairedDeck`, not a `SelfPairedDeck` (because only reflection is involution).
We're going to uplift, do all the checks and retire the old one when we're
successful."* Arrived from the question **"why does the periodic case need the
identity operator? Why not just the translational deck transformation — a law
uplift of the rigid-motion machinery in geometry?"**

**Goal (the OUTCOME).** A deck transformation's action on the trace is DERIVED
from the rigid motion, by one body, for every law that has one. Specular,
periodic and (latent) rotational sector symmetry then differ only in **which
motion** they carry — never in which code path realizes them.

#### Why the current tree is not that — `[M]` two measured gaps

1. **The deck slot is typed asymmetrically.** `SelfPairedDeck` holds
   `motion: RigidMotion`; `SpatialWrap` (periodic) holds `axis: str`. G5's own
   argument sits unapplied in `SelfPairedDeck`'s error message — *"the deck
   element is the transformation itself, not a name for one — an axis LETTER
   cannot say which plane"*. A letter equally cannot say which **lattice
   vector**, nor which **rotation angle**.
2. ⛔ **The specular arm does NOT go through the rigid-motion core.** `[M]`
   `_specular_kernel` builds its permutation from
   `quadrature.reflection_index(axis)` → `_find_reflections`, a
   nearest-neighbour `argmin` table — **ERR-074's own site** (*"never checked
   that the reflection PERMUTES the nodes"*). Meanwhile
   `RigidMotion.permutes` exists, is G2-verified, and **returns a
   bijectivity-checked `Permutation`** — the exact guard whose absence was
   ERR-073. ⟹ **two paths for one concept, and the boundary tier is on the
   unhardened one.** Cardinal Rule 2.

   ⛔ **REFUTED 2026-08-07 (reconciled against the tree at pick-up, §7
   discipline) — the HAZARD half of this claim was STALE.** `[M]`
   `_find_reflections` was retired at `a7695148` (2026-08-02, Q5.0.1):
   the sphere table has been CERTIFIED since — built through
   `_orbit_closure` → `RigidMotion.preserves` (bijection + weights, the
   full ERR-073 discipline), and GL1D's is a closed form. What SURVIVED,
   and what the step actually fixed: the **architectural** half —
   realization read the axis LETTER and never consulted the law's own
   motion (G5's deck slot unread), and only axis-aligned mirrors were
   realizable at all (the table's `range(3)` loop). The fix landed in a
   DIFFERENT campaign between this section's writing and its reading;
   the ⛔ was true when written. The re-check cost minutes; trusting it
   would have mis-framed the commit as a hardening instead of a
   consolidation.

#### The taxonomy, and the criterion that decides it

`[M]` 2026-08-07, against `SelfPairedDeck`'s live guard:

| motion | order | fixed_dim (d=3) | `SelfPairedDeck` |
|---|---|---|---|
| reflection | 2 | **2 = d−1** | accepted — this IS the self-paired case |
| rotation 90° | 4 | 1 | refused |
| **rotation 180° (half-turn)** | **2 — an involution** | 1 | refused |
| translation | ∞ | 3 | refused |

⭐ **Involution is NECESSARY but not SUFFICIENT for self-pairing.** The
half-turn is an involution and still maps a face to its *opposite*. The
deciding criterion is the **fixed subspace**: a self-paired face must LIE in
the motion's fixed set (`fixed_dim ≥ d−1`), which among the deck families only
the reflection satisfies. The existing guard already encodes this, and its
docstring already names the half-turn and the inversion as *"precisely
`SpatialWrap`'s deferred job"* — so the tree anticipated the sibling type and
then under-typed it.

⟹ **`PairedDeck(motion: RigidMotion)`** — the genuine-pair deck element,
covering translation (periodic), rotation (sector symmetry), half-turn and
inversion. `SpatialWrap` retires onto it.

⚠ **Rotational sector symmetry is a LATENT consumer, not speculation** — 1/4,
1/6, 1/8 core symmetry is ordinary reactor physics and the reason to build the
paired type generically rather than special-casing the wrap
(`process-discipline`: *"no current consumer" often means "consumer not yet
wired"*).

#### Why `IdentityOperator` was the wrong answer, stated once

A rigid motion is a linear part (acting on DIRECTIONS) plus a translation
(acting on POSITIONS). A torus wrap's linear part **is** the identity, so
`on_directions` returns its input — which is why the arm reached for
`IdentityOperator`. But the trace arrow is `Γ₊(f') → Γ₋(f)` and for a
translation **the face crossing is the entire content**. `IdentityOperator` is
an endomorphism `V → V` that names no spaces; it can never be an arrow between
two. The linear part being trivial is a **value of the motion**, not a
different kind of object.

#### Proposed shape (2026-08-07, NOT yet implemented)

```
PairedDeck(motion)      SelfPairedDeck(motion)
        \                      /
         one uplift: RigidMotion + trace.directions
                     |
             permutes()  ->  Permutation      (bijectivity-checked, ERR-073)
                     |
      PermutationOperator(π, domain=Γ₊(f'), codomain=Γ₋(f))
```

`[M]` the input is already there: **P6 gave every trace tier a `directions`
field** (`16835c68`, added for the inflow source) — exactly `permutes`'s
argument. The uplift is cheaper now than when G5 ran.

**Done when:** ✅ one body realizes every deck law (`_deck_kernel`,
`c3bb7341`); ✅ the wrap deck carries a motion (`PairedDeck`, `e13313a8` —
`SpatialWrap` retired outright rather than gaining a field); ✅ the boundary
tier no longer calls `reflection_index` (`[M]` zero refs in
`sn/boundary/realizer.py`); ⛔ **PENDING — the user's gate**:
`reflection_index` RETIRED across its consumers — **§7d below.**

#### §7d — the `reflection_index` retirement audit (the remaining substep)

*Goal:* one path for "ordinate permutation induced by a rigid motion" —
`Quadrature.ordinate_permutation` (landed) is the source; the certified
axis-mirror TABLE either becomes its private cache or its consumers
re-derive. *The user's sequencing ruling:* retire only after the checks
pass — they now have (142 rows green; wide slice at ⏸ #5).

`[M]` **The explorer's full census landed 2026-08-07**
(`scratch/audit_reflection_index_retirement.md` — the earlier "10 call
sites" and "33 live refs" were both partial): **7 production call sites,
44 test lines / 17 files, 5 `monkeypatch.setitem(quad.reflection_partners,
…)` injection sites, 2 duck-typed `_TableQuad` fakes, 21 doc mentions /
8 pages**, plus the Nexus equation node `quadrature-reflection-index`
(`angular_quadrature.rst:239`, 23 `tests` edges via the file-level
verifies list at `tests/sn/primitives/test_quadrature.py:35`) — the
marker-migration surface.

**The cache-vs-rederive question is ANSWERED by frequency class:**

- **HOT — 4 sites, `sn/loss_representation/__init__.py`**
  (`_apply_walk`:3294, `loss_action_transpose`:3470, `_run`:4189,
  `_run_transpose`:4593), all static ``"x"``, all the r=0 coupled-pole
  seed, sitting inside the `(L+C)` apply/solve/transpose bodies —
  10²–10⁴ lookups per solve; the TABLE is their cache. Migration =
  **hoist ONE derivation to walk/mesh construction** and let the 4 sites
  read the stored array — which also moves the σ_x-closure refusal from
  mid-iteration to construction, a strict improvement.
- **WARM — 3 sites, `geometry/boundary/_specular.py`** (:145/:184/:220,
  the ERR-042/044/045 certifications via `assert_realizable`,
  data-driven axis). O(N²) rederive is affordable — **but the collapse
  onto `preserves` is NOT lossless per-check**: ⚠ ERR-044's involution
  has NO motion-path carrier (`preserves` certifies bijection, not
  involution — assert it ON the derived π, `π∘π == id`, where the
  catalog names it), and ⚠ ERR-042's check is against the TRACE measure
  ``w·|μ_a|`` where `preserves` checks bare ``w`` — keep each
  certification's own claim; let them consume the derived π, never
  substitute the weaker certificate for the stronger check.
- **COLD — 1 site**: `derivations/…/sweep_acyclicity` takes the ARRAY
  as a parameter; only the test harness calls the accessor. Zero
  ceremony.

Also owed by the substep: **no gate pins the realize()-level refusal**
of an uncertified mirror (the deck-kernel refusal rows call the kernel
directly; `BoundaryError ⊂ ValueError` makes a lazy re-pin vacuous —
mint an explicit realize-path pin); **delete**
`test_deck_kernel.py::TestTheSpecularArmInheritsTheRetiredTable` (its
docstring carries its own deletion order); re-point the
`test_quadrature_directional.py:294` message pin;
`docs/api/discrete_ordinates.rst:67` is the ONE Sphinx-rendered xref
that dangles visibly (the other 20 doc mentions are grep-only); full
three-search audit + shortest-fragment message grep at the end.

✅ **7d.0 landed** `ae018628`: the audit's four found falsehoods
fixed on sight — two docstrings step 7 itself falsified
(`albedo.py`/`reflective.py` "realizes through ``reflection_index``" →
the motion spelling) and two `_find_reflections` citations dead since
Q5.0.1 (`curvilinear_numerics.rst:2228`,
`test_coupled_pole_mu_level_invariant.py:47`).

✅ **7d.1 landed** `a17771db` — the WARM tier: the three certifications
derive their array from the mirror motion (`_specular_pairing` →
`ordinate_permutation`), each keeping ITS OWN claim (the audit's
carrier findings honoured: the involution asserted on the derived π,
the TRACE-measure check unchanged); the five `reflection_partners`
injection sites re-posed onto `ordinate_permutation`
(instance-scoped `setattr`), the two `_TableQuad` fakes converted, and
the missing realize()-level refusal pin minted (message-fragment
matched — `BoundaryError ⊂ ValueError`).

✅ **7d.2 landed** `f11a129f` — the HOT tier: `_ensure_pole_mirror()`
on `_OneDimScanWalk` (the `_ensure_geom_cache` mesh-stash idiom)
derives the coupled-pole pairing ONCE from
`SelfPairedDeck.mirror(axis="x").motion`; the four per-sweep sites
read the stash. A rule with no x-mirror closure refuses at first use
with the pairing's diagnosis.

⟹ **EVERY production consumer of `reflection_index` is migrated.** The
table is production-DEAD but API-alive: its definition tier
(`directional.py`), its own gates, and the test/doc mentions remain —
a deliberate, coherent intermediate state.

#### ✅ 7d.3 — the remaining piece: DELETE the table and migrate its tier

`[M]` remaining census (2026-08-07, post-7d.2): **27 files, ~140
lines** — `tests/numerics/test_quadrature_directional.py` (20, the
table's own gates → re-pose onto `ordinate_permutation` or delete),
`directional.py` (18, the definition tier: `reflection_index`,
`reflection_partners` field ⚠ a DATACLASS field — direct-constructor
audit required, `_compute_sphere_reflection_partners`, the GL1D
closed-form partners), the two curvilinear sweep test modules (14+14,
mostly docstring mentions — tense-discriminate),
`sweep_acyclicity.py` (8, the COLD parameter surface + its harness
call), and ~12 smaller test/doc files. Plus:
`TestTheSpecularArmInheritsTheRetiredTable` deletes per its own
docstring; the `:294` message pin re-points; the Nexus equation node
`quadrature-reflection-index` (`angular_quadrature.rst:239`, 23
`tests` edges via `tests/sn/primitives/test_quadrature.py:35`)
re-poses; `docs/api/discrete_ordinates.rst:67` is the one RENDERED
xref; ~21 doc mentions across 8 pages tense-discriminate. End with
the three-search audit + shortest-fragment message grep + a full wide
slice. The worklist is `scratch/audit_reflection_index_retirement.md`.

✅ **LANDED 2026-08-07 in two green-tree commits** (§6b: the call-site
set, not the description, drew the boundary — consumers first, then
deletion+docs together because the xref checker gates red between):

- **`c7ca338e` (7d.3a)** — every SURVIVING test consumer off the table,
  discriminated by claim class: REFERENCE gates → the new independent
  helper `tests/_harness/references.py::mirror_partner_indices` (pure
  numpy negate-column + exact-match + injectivity; re-pointing them at
  production's source would have been the rewire-demotion);
  PRODUCTION-DATUM gates (coupled-pole contract, ξ-mirror battery,
  phase-C reference sweep, MMS fixtures, the ERR-044 mutant base) →
  the SAME expression production evaluates. The `:294` refusal pin's
  successor minted at the sweep tier:
  `test_a_rule_without_the_x_mirror_refuses_before_any_march`
  (product(4,5) + cylinder solve → ValueError, fragment
  `"cannot seed the r = 0 pole"` — that refusal was UNPINNED until
  here).
- **`b5ac130e` (7d.3b)** — definition tier deleted (`reflection_index`,
  the `reflection_partners` FIELD (constructor audit: 1 omission-style
  direct construction, safe; 4 factory passes died with the population),
  `_compute_sphere_reflection_partners`, GL1D closed form,
  `_resolve_axis_to_index` — its "Unknown axis label" refusal died with
  it; the vocabulary lives on `_mirror_motion`, already gated in
  `test_paired_deck`). Q4 gates re-posed keeping each claim (Q4.5 grew
  per-row `expected_axes` so availability is explicit and the loop
  cannot silently vacate; the witness now runs on the `(N,3)` embedding
  so degenerate axes are exercised). The oracle class deleted per its
  own docstring, tombstone comment naming the real surviving carriers.
  Equation label `quadrature-reflection-index` →
  `quadrature-ordinate-permutation` — ⭐ the old section still
  documented the PRE-Q5.0.1 bare-argmin formula (ERR-074's own body)
  and cited `_find_reflections`; the re-pose made the equation the
  certified definition `Ω_π(n) = QΩ_n ∧ w_π(n) = w_n`, node degree 29
  post-build. 8 pages + `_errors.py` tense-discriminated — ⭐ the three
  `_errors.py` "reflection-index table" diagnostics were caught only by
  the CONCEPT grep (hyphenated prose is invisible to the symbol grep).

`[M]` audit at `b5ac130e`: symbol grep = past-tense history + the
`sweep_acyclicity` PARAMETER name + negative mentions only;
message-fragment greps ("precomputed reflection", "Unknown axis label")
return only past-tense provenance; Nexus post-build: the resolved
method node is GONE, 3 residual unresolved edges all in UNTRACKED
`scratch/derivations/diagnostics/diag_issue168_*` (archaeology, not
consumers); sphinx `-W` clean, matrix regenerated; pyright ratchet
green; xref DEAD TARGETS 0. Touched batteries 4231 passed / 3 failed =
exactly the pre-existing baseline (2× `cart2d_*_principled_equiv` +
`WhiteXminPartial03GLSnapshot` #33).

**(a) The wide slice at this HEAD was NOT measured** — launched twice
in background (2026-08-07), externally stopped both times before
completing; per the session's stated protocol the second stop was
taken as deliberate and no third launch was made. ⚠ The full-slice
measurement is OWED and is assigned to the NEXT wide run (G6.3b's
landing gate, or task #51's green-the-baseline run — whichever comes
first). Evidence in lieu, all `[M]`: the touched-directory batteries
at `b5ac130e` (tests/numerics + geometry + sn/operators + sn/primitives
+ acyclicity + coupled-pole) = **4231 passed / 3 failed — exactly the
pre-existing baseline reds**; the sweep/mms trio 29 passed / 7 xfailed
(the #326 battery); every one of the census's 17 table-consuming test
files was individually migrated AND re-run green; production has been
table-free since `f11a129f`, whose wide slice measured 7/6163; and
`fc5c7a67` on top is prose-only. The un-re-run remainder of tests/
reaches 7d.3 only through production behavior, which 7d.3 did not
change — an EXPECTATION, not a measurement.

**(b) qa review returned 7/7 checklist PASS** (`[M]` demotion audit:
21/21 in-class reds on the reference gates under a to_local mutation,
none production-vs-itself; the refusal pin reds when the refusal is
removed; the ERR-074 catcher live; no consumer reads a deleted symbol)
with **eight prose-truth findings, all fixed at `fc5c7a67`** — four
present-tense-false "precomputes the reflection-partner map" claims
inside the gutted package (caught only by the FIELD-hyphenation concept
grep; no build severity renders those modules), one half-done docstring
correction, one false relative clause, one vocabulary-stale doc line,
and ⭐⭐ the measured independence-narrowing: `[M]` an axis-letter swap
inside `_mirror_motion` left `test_specular_response_pins_to_geometry`
all-green while 78 sibling gates red — both its sides now consume ONE
`SelfPairedDeck.mirror(axis).motion`, so it cross-checks the DERIVATION
but no longer the letter→normal convention (that lives in the reference
helper's deliberately-local `_AXIS_INDEX`). qa minted **vv-principles
#22**: independence has TWO axes — derivation and INPUT resolution
(the shared-constructed-input sibling of coding-standards'
rewire-demotion clause). qa also confirmed item-7b's carrier split:
the deleted oracle's mirror codomain-by-`is` leg rides on
`TestPeriodicIsBoundToThePartner` ALONE (`[M]` M5: dropping the
kernel's `codomain=` reds 5 periodic rows, split gate stays green —
its codomain assertion is None/None-satisfiable); the tombstone
comment says so.

✅ **The transitional pins flipped as designed** (2026-08-07): the
boundary-operator pin reddened with its own re-scope message and is now
the live wrong-face catcher; `test_b3_domain_narrowing`'s periodic xfail
XPASSed(strict) and the todo-list mechanism retired with its last row —
empty set, mechanism deleted, completion note in place. A THIRD flip the
plan had not tracked: `TestSpecMutationsPropagate`'s wholesale
`geometry_map → identity` mutations now redden at REALIZATION (the kernel
refuses an incoherent mirror — itself the proof the realizer stopped
reading the letter), so those probes moved one property up
(`permutes_ordinates`, the flag the consumers actually read — vv #18's
mutate-inside-the-legal-class).

### §7h.5 — step 8: the check fires, and what it does NOT reach

✅ **LANDED `84d7861f`.**

**Shipped.** `_reflect_trace` builds `face_action = law @ γ₊` per face and uses
that ONE operator for both legs, because `(law ∘ γ₊)ᵀ = γ₊ᵀ ∘ lawᵀ` is
`OperatorProduct`'s own transpose law. `[M]` forward and transpose
bit-identical on all four shipped law kinds; the wide slice is **7 failed /
6017 passed**, the same pre-existing seven.

⭐ **The transpose leg is now DERIVED, not written**, and that is the real
prize. The ⚠ trap the docstring had been carrying since B3.2 — *"the transpose
must scatter over Γ₊, never Γ₋"*, whose violation is bit-identical for
off-diagonal permutation laws and was caught historically only by one
grid-reciprocity arm on a het-VACUUM sphere — is no longer a choice the code
makes. γ₋ is not a factor of the thing being transposed.

⚠ **But "unspellable" was an OVERCLAIM, caught in review of my own prose.** A
hand-written `γ₋ᵀ(lawᵀ(·))` still runs silently, because `|Γ₊| = |Γ₋|` makes
the shapes agree. The true, narrower claim: re-opening the trap now requires
ABANDONING the composition — a visible structural edit — instead of choosing
the wrong one of two adjacent names. Both docstrings say the narrow thing.

**⛔ THE FINDING — composing NARROWS the contract from duck-type to morphism.**
`[M]` `_reflect_trace` previously only ever called `law.apply`, so a bare class
with an `apply` was a sufficient surrogate; `law @ γ₊` raises
`TypeError: unsupported operand type(s) for @` on one. Production is
unaffected (every `_face_laws` entry is realizer output, already an operator),
so the fix was to make the test's `_NoTransposeLaw` subclass `LinearOperator`
— i.e. honour the contract it stands in for — **not** to add a guard for a
case the type system now covers. Generalisable: *routing a call site through
the algebra silently raises its operand requirement from "has the verb" to "is
the type", and the only consumers that notice are duck-typed test doubles.*

**⭐ The flip-proof, and the honesty it forced.** Reverting `_reflect_trace` to
the sequential spelling reddens **1 of the 8** new rows —
`test_a_wrong_face_domain_map_now_RAISES_on_apply`, which re-injects B3.4c by
swapping `_face_domains` on a NON-periodic mesh. The other seven gate the
*algebra* the catcher rests on and cannot see the production change at all.
That split is now stated in the class docstring, because eight green rows
otherwise read as eight guards (`vv-principles` #18: *by what mechanism does
THIS gate see THIS property?*).

**`[M]` The cost, measured before anyone "optimizes" it away.** Composing per
call is **1.50×** on the face action itself (2.79 µs vs 1.86 µs; 0.94 µs of
pure `OperatorProduct` construction). At the tier that matters it is **0.42 %**
of a representative solve — a 40-cell 2-group slab with one reflective face
builds **530** face actions in 118 ms. ⚠ The step-8 wide run took **19:59**
against 18:15 for step 8.0; that is **NOT this change** — even a suite
consisting of nothing but that solve for 20 minutes would pay ≈5 s, not 104 s.
Machine variance. Do not add a `face_action` cache for it.

⚠ **Scope the check does NOT reach: periodic.** It realizes to an unbound
`I & I`, so one `None` short-circuits — and periodic is precisely the law
whose domain face differs from its installation face, i.e. the one the check
was designed for. **Step 7** binds it;
`test_PERIODIC_is_the_one_law_the_check_cannot_police_yet` is the transitional
pin and its failure message tells the next session to re-scope the ⚠ note in
`_reflect_trace`'s docstring alongside it.

### §7h.4 — step 8.0: the "plumbing" step that was carrying a live repair

✅ **LANDED `72f5ce97`.** ▶ **Track G now resumes at step 8** — its precondition
is discharged: `[M]` `_reflect_trace:553` spells the face action as
`law.apply(gamma_out.apply(face_in))`, and `law.domain == gamma_out.codomain`
now holds for every arm except periodic (still `I & I`, unbound — step 7's
job), so routing it through `@` will FIRE rather than skip.

**Shipped.** `TensorProductOperator` and `SumOfTensorProductsOperator` derive
their spaces; `OperatorSum`'s two inline transcriptions of the same check are
retired onto the one shared law, `_agreed_space`. The strict xfail
`test_the_realized_operator_carries_the_binding` lost its mark.

**The law, and why it is agreement rather than position.** A composite answers
`domain`/`codomain` from its operands in exactly one of two ways. `A @ B`
answers **by position** — correct, because swapping the operands makes a
different operator. Apply that to `A ⊗ B` and it self-contradicts: `A & I`
would be bound at the domain and unbound at the codomain while `I & A` is the
reverse, so an **order-independent operator would carry order-dependent
spaces**. The commutative composites (`+`, `⊗`) therefore resolve by
**agreement** — every operand that speaks must say the same thing, silence
contributes nothing, disagreement RAISES. `OperatorSum` already worked this
way; the tensor product was the outlier, and the fix made the two one body.

`[M]` agreement is sufficient *because the bindings here are WHOLE-space*: a
`PermutationOperator` tagged `axis=0` on a `(4,3)` trace declares
`domain.shape == (4,3)` — both axes. The factors are not separate legs; each
bound one describes the whole space and at most one can be non-trivial. The
refusal message names the product-space constructor that genuine per-leg
bindings would require, instead of electing a leg.

**⛔ THE SURPRISE — and it is about this plan's own prose, not the tree.**
The order-table row above sold 8.0 as an enabler: *"the binding is metadata,
not enforcement."* `[M]` **it is not metadata.** `_AdjointOperator.apply`
BRANCHES on whether the spaces are `None` to decide whether to apply the
metrics, so filling them in changes `.H` from the Euclidean transpose to the
partial-current Hilbert adjoint. On a 3-group `gauss_legendre(8)` xmin face:

| law | ‖bound − unbound‖ / ‖bound‖ |
|---|---|
| specular (mirror) | **0.0** — `G_{Γ₋} = G_{Γ₊}∘π` bit-exactly; the metrics cancel |
| Lambertian (white) | **0.87** |

So today's realized white boundary `.H` is *wrong*, and step 8.0 repairs it.
A session designing to the row would have shipped that repair silently: the
mutation battery confirms the boundary tier reddened **nothing** until the
consequence gates were added (`TestTheRealizedLawIsMETRICCorrect`, which
parametrizes BOTH laws and asserts the Lambertian fixture can still tell a
weighted adjoint from a bare one — a reflective-only fixture certifies a
Lambertian defect).

⟹ the transferable clause is now `plan-authoring` §8: **a field a consumer
branches on is an INPUT, not metadata** — a step that populates one is never
behaviour-neutral, whatever the plan calls it. Grep `is None` on the field
before believing the word.

**⭐ The mutation battery.** `[M]` all four on ONE slice — TP module + specular
chain + Lambertian factored + realizer + typed-residual (the `OperatorSum` arm)
— baseline **212 passed**, restored to 212 after:

| mutation | reds |
|---|---|
| drop the derivation (regress to pre-8.0) | **18** |
| swap `domain` ↔ `codomain` (positive control) | **8** |
| position rule instead of agreement | **17** — incl. the order-invariance row |
| remove the disagreement refusal | **5** |

⚠ **Two of these numbers were WRONG in an earlier draft of this table** (7 and
3), because the first three mutations were measured on a 189-test slice and the
fourth on 199, before the realizer consequence gates and the typed-residual
module were in it. Same defect `plan-authoring` §4 names: a number without its
configuration is not merely unreusable, it is *usable wrongly* — and here it
would have UNDER-stated the position rule's blast radius by 2.4×, which is the
direction that makes a gate set look thinner than it is. Re-measured on one
slice before landing.

**⛔ THE SECOND SURPRISE — a retirement audit that grepped its own wording.**
Folding `OperatorSum`'s inline check into the shared helper reworded its
exception message. `[M]` the audit grepped `requires equal domains`, found only
the definition site, and concluded no consumer. **Two gates
(`test_typed_residual_evaluation.py`) pin `"equal domains"` — the two-word
fragment — and both went red in the wide run** (9 failed vs the baseline 7; the
other 7 identical). A grep for the longer phrase cannot match a shorter pin, so
the audit was *strictly less sensitive* than the thing it was auditing, and it
returned an empty answer confidently.

Fixed by KEEPING the vocabulary (`f"{owner} requires equal {role}s"`) rather
than updating the gates to a new phrase — the phrase is load-bearing
provenance, and prefixing `owner` makes it say WHICH composite refused now that
the law is shared. ⭐ And the retirement moved the raise **one frame out**: the
gates asserted `co_qualname == "OperatorSum.__init__"`, which is now
`_agreed_space`. Re-pointing the frame alone would have DEMOTED them (the
helper is reachable on behalf of the tensor product too), so they now assert
the innermost frame **and** the caller. Clause added to `coding-standards`
retirement.

⚠ **Harness note — THREE instrument failures in this one step**, all of the
anti-pattern-#17 family (*the instrument lies before the code does, and it lies
in the safe-looking direction*), but by three different mechanisms:

| # | the instrument | what it reported | why |
|---|---|---|---|
| 1 | `tail -2 \| head -1` on pytest output | the pytest **docs footer**, six times | the summary line is the LAST line, not the second-to-last |
| 2 | its rewrite | `1 warning in 0.01s` = "clean" | the run collected **0 tests** — an empty run and a green run look identical |
| 3 | ⭐ `until ! pgrep -f "pytest tests/cp …"` | **"still running" forever** | `pgrep -f` matched the WAITER'S OWN command line, which contains the pattern — the liveness check detected itself |

⭐ **(3) is the new one and the nastiest**, because a never-terminating wait is
indistinguishable from a slow suite: two monitors sat in it until their
timeouts, and the run they were watching was killed unobserved at 58 %
(0 failures at that point). **A liveness check must not be able to match
itself** — filter by log content or let the job's own exit notify, never
`pgrep -f` on a string your waiter also carries. Belongs in `vv-principles`
**#17** alongside the fail-open parser.

The battery was only trusted after the baseline printed **199 passed** — #17's
positive control, applied to the instrument rather than the code.

### §7h.3 — step 5: bound, a flag retired, and step 8 refuted from below

**Shipped.** `_specular_kernel` returns its `PermutationOperator` bound
`Γ₊(f) → Γ₋(f)`; `_outflow_restriction`'s `γ₊` bound `Γ(f) → Γ₊(f)` with it, so
`P @ γ₊` now types the WHOLE face action end-to-end. `inverse()` **inverts** the
binding (it was dropping it). `PermutationOperator.is_involution` RETIRED.
Gates: `tests/sn/operators/test_specular_deck_chain.py` (79 + 1 strict xfail),
plus 2 rows in `test_permutation_operator.py`.

**⛔ THE FIFTH REFUTATION — and this one refutes a FUTURE step, not this one.**
`[M]` `TensorProductOperator` derives NOTHING from its factors, so
`SNBoundaryRealizer.realize(...)` returns `domain=None`. The binding is real at
the inner permutation and invisible at the object the realizer hands out. Step 8
("route `_reflect_trace` through `@` so the check FIRES") composes exactly that
object, and one `None` short-circuits the check — **so step 8 as written would
have landed, passed, and gated nothing.** The already-committed step-3 Lambertian
chain has the identical hole. New substep 8.0; committed as a strict xfail whose
flip-proof is measured (`XPASS(strict)` when `TP.domain → ops[0].domain`).

⭐ **The transferable form: a step can be refuted by a LATER step's precondition,
and the campaign's own order table is what hides it.** Steps 1–7 each verified
their own tier and were each correct there. The hole lives in the seam that step
8 was scheduled to exercise — so it was invisible to every step that came before
it and would have been *invisible to step 8 too*, which would have gone green.
Refutations 1–4 came from reading a claim against the realization; this one came
from asking **"what will the step AFTER this one actually compose?"** Ask it one
step ahead, not at the step.

**⭐ The `is_involution` retirement — a flag whose VALUE tracked the quadrature.**
`[M]` for ONE law (mirror about `x` on `xmin`) the narrowed permutation satisfies
`perm[perm] == arange` on `gauss_legendre(4/8)`, `product(4,4)`,
`level_symmetric(6)` — and NOT on `lebedev(17)`. The physics does not vary with
the quadrature; the value tracked how `to_local`'s `searchsorted` orders the
locals. The user ruled RETIRE over REFINE, and the reason generalizes: a refined
flag still answers, and an answer can be wrong; `P @ P` **cannot be formed**, so
routing the question through the algebra replaces a value that can drift with a
composition that cannot exist. Rewire, not delete — the two tests became "the
square is the identity", asked of `@`.

⚠ Blast radius the retirement actually had: TWO sentences reading *"a permutation
is a bijection (invertible, **involution-detectable**, …)"* — one in
`operator.py`, one in the theory page — became present-tense-false with the
attribute. Neither contains the symbol `is_involution`. **Grep the CONCEPT.**

**⭐ The mutation battery, before and after.** Baseline `3 failed / 1747 passed`:

| mutation | reds BEFORE the gates | reds AFTER |
|---|---|---|
| identity permutation (positive control) | +23 | **+28** |
| ⭐ **swap `domain` ↔ `codomain`** | **0** | **+10** |
| bind both ends to `Γ₊` | **0** | **+18** |
| drop the binding | **0** | **+31** |

The swap is the one that matters: a one-word slip, invisible to
`checked_space_extent` (`|Γ₊| == |Γ₋|` everywhere), zero arithmetic change — and
after step 8 it would make the LEGAL composition raise and the ILLEGAL one pass.
Its 10 reds land **entirely** in `TestTheBindingIsTheRightWayRound`; nothing else
in the slice can see it. `⟹` **an `is`-identity assertion naming which space is
which end is the only instrument that can exist for this class.**

⚠ And the swap does NOT red the `.H` rows — correctly. `G_{Γ₋} = G_{Γ₊}∘π`
bit-exactly on all five quadratures (a mirror preserves `|Ω·n|·w_n`), so the
swapped adjoint is numerically the unswapped one. **The metric is blind to the
swap for the same reason `.H` is blind to the metric here** — which is why the
blindness criterion is pinned as its own gate rather than left implicit under a
`.H ≈ transpose` row that would silently be pinning a coincidence.

⚠ **Harness note, twice in one session and both from the anti-pattern-#17 list:**
a `FAILED`-line scan returned nothing under `-q --tb=no` (no such lines emitted),
and then returned nothing again under `-rf` because **ANSI colour codes break
`^FAILED`**. Use `--color=no -rf`. Both times the empty result read as "no reds".

### §7h.2 — step 4: the arm collapse is IMPOSSIBLE, and what shipped instead

Step 4 was scoped as *"put `α` on `S(f)`; this collapses `_attenuated_kernel_operator`'s
THREE arms into chain manipulations — the chain as-is / prepend a scalar link /
**the zero chain**."* `[M]` **both halves are wrong.**

⛔ **The zero arm cannot collapse.** `ScaledOperator` *refuses* a zero scalar:
`ValueError: ScaledOperator with zero scalar is degenerate; use ZeroOperator
explicitly.` The guard is deliberate and the α=0 arm's own docstring already
said so — *"**Not** a convenience: the general path cannot express it."* The
three arms are not structural debt awaiting a chain; they are three genuinely
different objects, and *a surface that returns nothing IS a vacuum*.

⚠ **I measured `0.0 * chain` in NUMPY and reported it as free.** The operator
algebra refuses it. Same error as §7h.1's, one level down: *measuring the wrong
tier and generalising*. The chain being correctly typed makes the ZEROS right;
it does not make the SPELLING legal.

⛔ **And moving `α` onto `S(f)` is not uniform.** The specular chain has length
1 and no interior, so `α`-on-`S(f)` is a per-law-kind placement — the arm split
returning. `[M]` it also costs a re-baseline: `α·(J/n)` vs `(αJ)/n` differ
**34.8 %** of the time (≤2 nulp) over random floats, and 35 % at α=0.3
specifically — though **0 %** at α=0.5, because a power of two multiplies
exactly. Two fixtures said "bit-identical"; they were lucky. The affected
snapshots include `TestWhiteXminPartial03GLSnapshot`, which is the α=0.3 white
case **already red as task #33** — so a re-baseline would let the carve
re-capture the very anchor meant to pin it.

⟹ **α stays on the composite.** Its position is unobservable except in ULPs
(scalars commute), so moving the multiply buys interpretation at the price of a
re-baseline and a reintroduced branch. The `J⁻ = αJ⁺` reading is documented
where it belongs instead.

### ⭐ What step 4 actually delivered — the user's pin (2026-08-04)

> A deck transformation admits **no attenuation** — and that is the *proof* it
> is geometry. The specular reflection as a RESPONSE is a stand-in for a real
> mirror, which we may attenuate because real mirrors are not ideal. That we
> can check the response against the geometric reflection is a **special case**
> where the two can be pinned against each other. Even a scalar attenuation is
> already an imposed response — real mirrors rarely attenuate uniformly.

`_factors.py` already flagged `ReflectiveBoundary(axis, albedo≠1)` as a
deliberate violation, justified by *taxonomy*. This is the **consequence** form,
which is stronger: **if it can be attenuated, it was never a quotient.**

`tests/geometry/test_specular_response_pins_to_geometry.py` (15 gates) makes the
special case a gate. The two sides are independently derived — the geometric one
applies `SelfPairedDeck.mirror`'s `RigidMotion` (the G1–G5 core, verified against
pure math) to the direction cosines and reads off the induced permutation; the
response side is `quadrature.reflection_index(axis)`, the rule's own table.
Neither consults the other. `[M]` they agree EXACTLY on `gauss_legendre(4)`/`(8)`,
`product(4,4)` (x, y), `lebedev(17)` (x, y) and `level_symmetric(6)` (z).

⭐ **A finding from the teeth check worth keeping: the EXACTNESS leg catches a
class the equality leg cannot.** A scaled mirror `diag(-1,1,1)·1.1` induces the
**correct permutation** while its images are not nodes — `matches_table=True`,
`exact=False`. Equality alone passes it. The residual is what enforces the
**measure-preserving** half of the definition: a scaling is not an isometry, so
it is not a deck element however right its combinatorics look.

### §7h.1 — the transferable lesson: TWO false premises, both from a docstring

G6.3 was scoped twice on claims read from prose, and refuted twice by the
realization:

1. "`G` and `R` are operators" — they are **descriptors**; the realizer welds
   `R∘G`. Caught by the mapping dispatch.
2. "`R : Γ₋ → Γ₋`, so the split is `R @ G` with intermediate `Γ₋`" — **no realized
   response is an endomorphism of `Γ₋`**; the response *carries the crossing*.
   Caught by the USER asking what space sits between the two operators.

⭐ **The lesson: a Protocol docstring is a DECLARATION, not a measurement — and
when a tier declares a typing, CHECK IT AGAINST THE REALIZATION before designing
on it.** Both refutations were sitting in one `[R]` read of
`AngularAverageOperator`'s first line, which types itself `Γ₊ → Γ₋` in direct
contradiction of the Protocol it implements. The declaration/realization gap is
exactly what this campaign exists to close, so it should have been the FIRST
thing checked, not the last.

Corollary for the fleet: the sub-agent DID report "`R : Γ₋ → Γ₋` is un-bindable —
nothing in the SN realizer is an endomorphism of `Γ₋`" and the main agent
under-weighted it, treating it as a scoping detail rather than a refutation of
the algebra. **A sub-agent's negative structural finding deserves the same
scrutiny as a positive one.**

⛔ **Do NOT descend into the `PermutationOperator` / `AngularAverageOperator`
factors.** `[M]` `Γ₊(f)` is already the product space (`(12,2,7)` on 2-D
2-group), and `is_involution` becomes a type-incoherent claim on a
non-endomorphism.

### §7h.0 — the base-class defect G6.1 exposed (fixed in `34f465cc`)

`FunctionSpace._diagonal_inner_product` used numpy's default **TRAILING**
broadcast while `_diagonal_apply_metric` used the documented **LEADING**
convention — *the same metric applied along different axes by two methods of one
space*. Invisible whenever `w.ndim >= x.ndim` (every case the tree exercised),
so it shipped. `[M]` non-square shapes RAISED from `inner_product` while
`apply_metric` worked (`SphericalHarmonicSpace.from_L(3)` at its production
`(L+1, 2L+1, ng, *spatial)` layout); square shapes **silently disagreed** (456 vs
552 on a `(3,3)` probe). Either way `⟨Ax,y⟩ = ⟨x,A†y⟩` is false by construction,
since `_AdjointOperator` builds `A†` from `apply_metric` while the pairing
judging it came from `inner_product`. **The ERR-067 family one layer down.**
Fixed by routing both through `_broadcast_metric`; bit-identical over 2294 tests,
as its no-op argument predicts.

⚠ **Wants an `ERR-NNN` catalog entry — NOT written**, because
`.claude/skills/vv-principles/error_catalog.md` carries uncommitted-by-policy
state forbidden to commit. Flagged for the user; recorded meanwhile in
`scratch/g6_design_measurements.md`.

---

## ⏸ COMPACTION POINT #2 — G6.3 steps 0–4 done, steps 5–8 remain

Supersedes #1. A session picking up cold re-anchors from **this file + `git log`**,
never from a conversation summary.

- **HEAD `78e6b289`** on `refactor/operator-strategy-layers`. **Tree clean.**
  18 commits since #1; all verified ancestors. Verify any hash you rely on with
  `git merge-base --is-ancestor <hash> HEAD`.
- **What landed** (G6.3 step → commit):
  | step | commit | what |
  |---|---|---|
  | scope | `979a6bb9` `e8885d02` `df503d3b` | G6.0 survey; the design of record; the CHAIN reframing |
  | G6.1 | `34f465cc` | the Γ ladder + a base-class metric-convention fix |
  | 1 | `bff3bd96` | `γ±` bound; the composability check FIRES |
  | 2 | `38466932` `c17deea3` | `S(f)`; and `ker(G) == tangential` — one theorem, not two facts |
  | 3 | `b4290873` `47c33d7b` `b4f0f5c9` `7d70e5b5` | the Lambertian factored; realizer spells the chain; `AngularAverageOperator` retired; corpus repointed |
  | 4 | `78e6b289` | the specular RESPONSE pinned against the geometric deck transformation |
- **REMAINING: steps 5–8** (§7h order table). **Step 8 is the one that matters
  most** — until `_reflect_trace` composes through `@`, the binding is honest
  metadata that gates nothing (`[M]` 4941 bindings, ZERO new failures, because
  the path uses raw `.apply`).
- **Known reds — still exactly 4, unchanged all session**, and none this
  campaign's: `TestWhiteXminPartial03GLSnapshot::test_matches_the_frozen_scaled_lambertian`
  (~1 ULP, task #33); `test_cart2d_1g_vacuum_apply_principled_equiv` (1152 ULP);
  `test_cart2d_2g_specular_apply_principled_equiv` (296 ULP);
  `TestBitIdenticalCurvilinear::test_spherical_inward_bit_identical`
  (`assert False`).
- **Gate costs, measured THIS session:** `tests/geometry` ≈10 s ·
  `tests/geometry + tests/sn/operators` ≈25 s · `+ tests/numerics` ≈5 m 45 s ·
  `+ tests/sn/sweep` ≈4 m 50 s · the wide slice
  (`geometry + numerics + sn/{operators,sweep,architecture} + diffusion`)
  ≈9 m 45 s. ⚠ **Always pass `-m "not slow"`** — without it a `tests/sn` run
  exceeded 26 min and is not comparable to the red baseline, which was
  established WITH the deselection.
- ⚠ **Verify a test path exists before running it.** A run that collects
  nothing exits 0 in 0.01 s and looks identical to a green one. This bit twice
  this session (`tests/sn/boundary` does not exist).
- **Static gates:** `tools/check_docstring_xrefs.py orpheus tests docs` →
  `DEAD TARGETS : 0` · `sphinx -E -W` → 0 warnings (use `-E`, never
  `rm -rf docs/_build`) · the pyright **ratchet test** is the gate — bare
  `npx pyright` reports ~2000 errors because it scans `scratch/`.
- ⛔ **Never `git checkout <path>` / `git restore` / `git stash` / `git clean`
  on a tracked path** (lesson L28). Compare via `git worktree add`.

### ⭐⭐ The meta-lesson, and it is the session's most transferable output

**G6.3's scope was refuted FOUR times, always the same way: a claim read from
PROSE, contradicted by the REALIZATION.**

| # | the claim | the refutation |
|---|---|---|
| 1 | "`G` and `R` are operators" | they are **descriptors**; the realizer welds `R∘G` |
| 2 | "`R : Γ₋→Γ₋`, so split `R @ G` through `Γ₋`" | **no realized response is an endomorphism of `Γ₋`** — the response CARRIES the crossing |
| 3 | "the factorization is per-law-KIND" | it is a **chain**; atomic is a degenerate LENGTH, and the machinery (`@`, `+`) already existed |
| 4 | "α on `S(f)` collapses the three arms into … the zero chain" | `ScaledOperator` **refuses** a zero scalar; the arms are three different objects |

Three were caught by the USER asking a question, not by a gate. ⟹ **Before
designing on a Protocol/docstring typing, CHECK IT AGAINST THE REALIZATION** —
one `[R]` read of the implementing class's first line would have caught #1 and
#2 together.

**And its sharpened form, which bit twice more:** *measure at the tier you will
build on.* Refutation #4's first measurement was `0.0 * chain` in **numpy**,
reported as free; the **operator algebra** refuses it. Same shape as measuring
a doc claim instead of the code. A measurement at the wrong tier is not weak
evidence — it is evidence about a different question.

**Corollary for the fleet:** a sub-agent DID report "`R : Γ₋→Γ₋` is un-bindable"
and the main agent under-weighted it as a scoping detail. **A sub-agent's
negative structural finding deserves the same scrutiny as a positive one.**

### Other lessons this stretch, ranked by transferability

1. ⭐ **A correct composite is NOT evidence its factors are right.** The theory
   page derived the correct `Rᵀ` from two WRONG per-factor formulas; the
   composite measured exactly 0.0 and no gate could see it. Dual of the
   metric-cancellation theorem — interior errors cancel just as interior metrics
   do. Gate each factor separately.
2. ⭐ **Fixtures agree by luck; check whether it is a contract.** `α·(J/n)` vs
   `(αJ)/n` was bit-identical on two fixtures and differs **34.8 %** of the time
   in general (0 % at α=0.5 — powers of two multiply exactly).
3. ⭐ **A mutation the property cannot SEE proves nothing either way.** A
   teeth-check that came back bit-identical read as "the reference is
   toothless"; the mutation was inside the stabiliser (symmetric hemispheres).
4. ⭐ **A retirement's rewire can DEMOTE a gate silently** — re-pointing a
   bit-identity reference at the successor compares a value with itself. Hit
   TWICE in one retirement, from different directions.
5. **Grep the CONCEPT, not the symbol.** A claim inversion survived 70 lines
   from its own cause because it named the OTHER class, not the property.
6. **The exactness leg catches what the equality leg cannot** — a scaled mirror
   induces the right permutation with images that are not nodes.

---

## 8. Standing constraints for the implementation session

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

---

## ⏸ COMPACTION POINT #3 — track G resumes at G6.3 step 8.0

A session picking up cold re-anchors from **this file + `git log`**, never from a
conversation summary. This point does not supersede #2 — #2's red set, gate
costs and scope lessons all stand. It records what changed **while track G was
paused**, which is the part a cold reader cannot reconstruct.

### State

**HEAD `c3475228`** on `refactor/operator-strategy-layers`. Tree clean.
Verify any hash with `git merge-base --is-ancestor <hash> HEAD`.

### ⛔⛔ READ FIRST — track G's OWN objects changed while G was paused

Two campaigns ran on this branch between #2 and here, and both reached into
track-G surface. A cold session that assumes G's world is as it left it will be
wrong about three things:

1. **The Γ ladder grew a field.** `AngularTraceSpace` and
   `AngularFaceTraceSpace` now carry **`directions`** — the ordinate
   :math:`\Omega_n`, restricted per tier by the SAME expression the metric uses
   (`16835c68`, affine-campaign P6 step 1). G6.1 minted that ladder; it is not
   the object G6.1 left. `[M]` dimension-generic: `(m,)` in 1-D, `(m, 3)` in
   3-D; all three tiers carry it; four gates in
   `tests/numerics/test_angular_face_trace_space.py`.
2. **`level_symmetric` weights CHANGED at S6 and above** (`df33913d`, #327 Fork
   A). Per-:math:`O_h`-orbit moment-matched, degrees now 5/7/9/11 at
   S6/S8/S10/S12. **S2 and S4 are bit-identical** (single orbit ⟹ forced), which
   is why most of the tree did not move. ⛔ **The family now REFUSES above
   S12** — `tests/numerics/test_symmetry.py` and `test_rules_sphere.py` had
   their S16/S20 sweeps moved to S12, and
   `tests/sn/mesh/test_radial_characteristic_slot_coordination.py` too. A
   G-step that reaches for `level_symmetric(16)` will raise at COLLECTION time.
3. **A new always-on rule exists**: `.claude/rules/plan-authoring.md`
   (`f0c03984`), a LIVING document whose target is surprises-per-campaign → 0.
   It governs how this file should be written from now on, and it already has a
   clause (§6b) that bears on step 8.0 — see below.

### Step 8.0 — the next step, verified ready

**Goal.** `TensorProductOperator` derives `domain`/`codomain` from its factors,
so the binding survives to the object the realizer hands out.

**Why it blocks step 8.** `[M]` TP drops the binding, so
`SNBoundaryRealizer.realize(...)` returns `domain=None`, and **one `None`
short-circuits the composability check** — step 8 as written would land, pass,
and gate nothing. `[M]` 4941 bindings measured, zero failures: a green that
means nothing.

**The anchor, re-verified at `c3475228`:**
`tests/sn/operators/test_specular_deck_chain.py::test_the_realized_operator_carries_the_binding`
is `@pytest.mark.xfail(strict=True)` and still **XFAIL** (not XPASS). Flip-proof
already recorded in #2: patching `TensorProductOperator.domain -> ops[0].domain`
turns it XPASS(strict).

⚠ **TP is a SHARED numerics primitive** with consumers well beyond the boundary
tier. It needs its own gate set and a **wide** re-run, not a boundary-scoped
one. Enumerate its call sites before designing — and per
`plan-authoring` §6b, **check whether the step boundary cuts across them**: if
step 8.0 and step 8 turn out to share a call-site set, they are one step, not
two. (That exact trap cost a cycle in the affine campaign's P6.)

### The red baseline — UNCHANGED at 7, and now better understood

`[M]` `tests/{numerics,sn,transport,geometry,diffusion} -m "not slow"` =
**7 failed / 5988 passed**, ≈17–20 min. Same seven as #2:

* 2× `test_streaming_operator` `cart2d_*_principled_equiv` (1152 / 296 ULP vs a
  256 budget)
* 3× `test_affine_carve_bit_identity` (**sha256** — no magnitude available; the
  instrument defect is **#333**)
* 1× `test_diamond` `spherical_inward_bit_identical`
* 1× `TestWhiteXminPartial03GLSnapshot` (task #33)

⭐ **They are now fully adjudicated and UNBLOCKED** — see
`.claude/plans/green_the_seven_quadrature_reds.md`. Verification criterion (1)
is discharged (L0 nodes vs `numpy.leggauss` 1–3 ULP; L0 exactness both at ε with
**neither implementation dominating**; **462 L1/L2 SN tests green**). The
"defer behind #327" rationale in that file is **REFUTED in place** — all four
`level_symmetric`-riding reds use `sn_order=4`, where Fork A is bit-identical.
Only the re-capture and #333 remain, not investigation.

### Gate costs, measured this session

| gate | cost |
|---|---|
| wide slice `{numerics,sn,transport,geometry,diffusion} -m "not slow"` | 17–20 min |
| `tests/numerics` alone | ≈5.5 min |
| `tests/sn -m "(l1 or l2) and not slow"` | 6.6 min, 462 passed |
| `sphinx -E -W` | ≈4 min, 0 warnings |
| pyright ratchet (`tests/test_pyright_ratchet.py`, NOT bare `npx pyright`) | ≈17 s |
| `tools/check_docstring_xrefs.py orpheus tests docs` | DEAD TARGETS 0 |

### ⭐⭐ Durable lessons from the paused interval

1. ⭐⭐ **Ask what algebra the CARRIER has before asserting an algebraic law.**
   P5's plan demanded `B(x+y) = B(x)+B(y)`; it RAISES — flux states are an
   affine torsor with no origin. The honest law is base-point independence of
   the increment. A spec written in the wrong algebra fails loudly if you are
   lucky and silently if the carrier is permissive.
2. ⭐⭐ **A guard whose subject can DECLINE it is not a guard.** `[M]` a source
   returning zeros at the ERR-047 rank-1 probe and `7.0` at the delivery shape
   skipped certification and delivered anyway. Retired by making the claim
   structural, not by hardening the probe.
3. ⭐⭐ **Count the ORBITS before costing a "re-implement this primitive" fix.**
   Six lines of `np.unique` over sorted `|node|` triples refuted a committed
   deferral rationale: S2/S4 are a single orbit, so the fix is *forced* to be a
   no-op at the most-used order. Blast radius went from "~200 sites" to ONE
   baseline.
4. ⭐⭐ **The imperative `pytest.xfail()` can never XPASS.** It aborts
   unconditionally, so a fixed rule with a stale tag keeps reporting XFAIL. Use
   `pytest.mark.xfail(strict=True)` — and FLIP-PROOF it. #327's nine strict
   xfails turned into nine XPASS failures the moment the weights landed, which
   is what forced the tag correction into the same commit. **Step 8.0's anchor
   is a mark, not imperative — confirmed at `c3475228`.**
5. ⭐ **Capture full, filter at READ time.** A pytest log filtered at write time
   showed "zero occurrences of the refusal message" and nearly produced a wrong
   diagnosis; the evidence had been destroyed at capture.
6. ⭐ **A symbol grep misses a value that is a LOOP VARIABLE.** `level_symmetric(16)`
   found nothing at two sites spelled `for sn_order in (4, 8, 16)`. Same
   subject-and-value-on-different-lines mechanic as `vv-principles` #21.
7. ⭐ **Both-directions wrongness identifies a formula for a DIFFERENT object.**
   `level_symmetric` over-claimed from S6 and *under*-claimed at S2. A number
   wrong both ways is not a mis-calibrated measurement of this thing.

---

## ⏸ COMPACTION POINT #4 — track G resumes at G6.3 step 7 (periodic)

Additive to #3 (which is additive to #2). Nothing in #2's red set or #3's
"track G's own objects changed while paused" section is superseded — but #3's
*"Step 8.0 — the next step"* block is now HISTORY, and so is step 8. A session
picking up cold re-anchors from **this file + `git log`**, never a summary.

### State

**HEAD `811b6407`** on `refactor/operator-strategy-layers`. **Tree clean.**
`[M]` all five commits of this stretch verified ancestors:

| commit | what |
|---|---|
| `72f5ce97` | step 8.0 — the tensor product's spaces, by AGREEMENT |
| `b89ccd2d` | its plan/rules/ERR-076 record |
| `84d7861f` | step 8 — the face action is COMPOSED |
| `4a87cc34` | its plan/rules record |
| `811b6407` | V&V matrix regen |

### ⛔⛔ READ FIRST — three things a cold session would get wrong about step 7

1. ⛔ **The order-table row's line numbers were STALE and are now corrected.**
   It cited `_assert_wrap_identification` at `:512` and its call site at
   `:851`; `[M]` 2026-08-07 they are **`:506`** and **`:970`**. The *claims*
   held (it does return the partner; the call site does discard it) — only the
   anchors rotted. ⟹ **cite the SYMBOL, never the line.** A line number is the
   one field in a plan that rots without anyone editing it.
2. ⛔⛔ **The goal is NOT "bind periodic" — that framing was mine and it was
   backwards** (corrected 2026-08-07 after the user asked whether this session
   contradicted the operation-chain concept). Under the concept a boundary law
   is a **`@`-chain** from :math:`\Gamma_+` to :math:`\Gamma_-` and a deck
   transformation is the degenerate **length-1** case
   (corpus: `bc-deck-length-one-chain`). Periodic is such a case. So the
   OUTCOME is: **periodic's single link becomes a typed arrow, the same kind of
   object every other law's link already is.** The binding is what you GET from
   that, not the objective.

   `[M]` the survey that shows the gap — angular factor of each realized law:

   | law | link | typed arrow? |
   |---|---|---|
   | specular, albedo-specular | `PermutationOperator` | ✅ |
   | white / Lambertian | `OperatorProduct` (2 links) | ✅ |
   | vacuum, prescribed | `ZeroOperator` | ✅ |
   | **periodic** | ~~`IdentityOperator`~~ → `PermutationOperator(arange)`, bound Γ₊(f')→Γ₋(f) | ✅ **since `c3bb7341` — `[M]` all five arrows typed; the survey's gap is CLOSED** |

3. ⛔⛔ **AND THE STEP WAS RE-SCOPED AGAIN, 2026-08-07 — read §7h.6 BEFORE
   anything below.** The user's question *"why does periodic need the identity
   operator at all — why not the translational deck transformation, a law
   uplift of the rigid-motion machinery?"* re-scoped step 7 from a
   periodic-only fix to **uplifting the deck transformation once**: mint
   `PairedDeck` (translation, rotation, half-turn, inversion — only the
   REFLECTION is self-paired), derive the trace arrow from the MOTION via the
   bijectivity-checked `RigidMotion.permutes`, and retire `reflection_index`
   (ERR-074's site, 10 call sites incl. the sweep) once the checks pass. The
   identity-permutation sketch in item 4 below is now an IMPLEMENTATION DETAIL
   of that uplift — the arrow for the translation case — not the step's goal.

4. ⭐⭐ **The reframing DISSOLVES the "blocker" I first recorded here.**
   The original note read *"`IdentityOperator` has no `__init__` and cannot be
   bound, and G6.0 ruled do not touch it"* — true, and beside the point.
   `IdentityOperator` is an **endomorphism** :math:`V \to V`; periodic needs an
   **isomorphism between two DIFFERENT spaces**,
   :math:`\Gamma_+(f') \to \Gamma_-(f)`. The identity is polymorphic
   precisely *because* it names no spaces, so it can never be that arrow. You
   do not bind it — you **replace it with the arrow**. ⟹ the unbindability is
   not an obstacle to route around; it is the diagnosis. Conflating *"the
   action is the identity"* with *"the operator is `IdentityOperator`"* is what
   made this arm the odd one out.

5. ⭐ **The arrow, measured 2026-08-07** — a permutation is the trivial
   relabeling between two DISTINCT index sets, which is exactly what a torus
   wrap induces (`_assert_wrap_identification` proves the local index needs no
   relabelling, so the permutation is `arange`):

   ```python
   PermutationOperator(np.arange(n_inflow), axis=0,
                       domain=Γ₊(partner), codomain=Γ₋(face)) & IdentityOperator()
   ```

   `[M]` on a GL-8 2-group slab: `Γ₊(xmax)` and `Γ₋(xmin)` are both `(4, 2)`
   so the CROSS-FACE arrow is well-posed; it is **value-identical** to
   `IdentityOperator` on `apply` AND `apply_transpose`; and the product carries
   `Γ₊(xmax) → Γ₋(xmin)` through step 8.0's agreement law. It also puts
   periodic in the SAME SHAPE as specular — which is not a bonus, it is the
   point: one body for every deck transformation.

   ⚠ **The one measured difference:** `IdentityOperator().apply(x) is x` is
   **True**; the permutation returns a **copy**. Anything relying on
   identity-aliasing of the periodic face would change. Check before landing.

   ⚠ Threading: `_assert_wrap_identification` RETURNS the partner face and the
   call site discards it — that part of the old row still holds.

*(Proposed means, 2026-08-07, NOT yet implemented — the measurements above are
real, the design decision is not made.)*

### The gate that will flip, and what it obliges

`tests/sn/operators/test_sn_boundary_operator.py::TestTheFaceActionIsCOMPOSED::
test_PERIODIC_is_the_one_law_the_check_cannot_police_yet` asserts periodic is
UNBOUND. It is the transitional pin and it **must redden when step 7 lands**;
its failure message says to re-scope, alongside it, the ⚠ note in
`_reflect_trace`'s docstring that says the same thing. Two places, one claim.

⭐ **Why periodic matters more than "the last arm".** It is the ONLY law whose
domain face differs from its installation face — the only off-diagonal block of
`B` — which makes it *exactly* the law step 8's composability check was
designed to police, and the only one it currently cannot (one `None`
short-circuits). Until step 7, the check is live on the four diagonal arms and
absent on the one that historically carried the defect (B3.4c, `[M]` 98 %
relative error when live).

### Red baseline — UNCHANGED at 7, all pre-existing

`[M]` `tests/{numerics,sn,transport,geometry,diffusion} -m "not slow"` =
**7 failed / 6017 passed**, 19:59. Same seven as #2/#3: 2× `test_streaming_operator`
`cart2d_*_principled_equiv` · 3× `test_affine_carve_bit_identity` (sha256, the
instrument defect is **#333**) · 1× `test_diamond` `spherical_inward_bit_identical`
· 1× `TestWhiteXminPartial03GLSnapshot` (task #33). Adjudicated and unblocked —
`.claude/plans/green_the_seven_quadrature_reds.md`.

### Gate costs, measured this stretch

| gate | cost |
|---|---|
| wide slice `{numerics,sn,transport,geometry,diffusion} -m "not slow"` | **18–20 min** (18:15 / 19:59 — the spread is VARIANCE, see §7h.5) |
| `tests/sn/operators + tests/geometry + tests/transport` | ≈50 s |
| `tests/sn/operators + tests/geometry` | ≈25 s |
| `sphinx -W` (incremental) | ≈1 min, 0 warnings |
| pyright ratchet | ≈17 s |
| `tools/check_docstring_xrefs.py orpheus tests docs` | DEAD TARGETS 0 |

⛔ **The complement of the wide slice (`cp`, `moc`, `mc`, `derivations`, …)
cannot be run whole**: `[M]` long background jobs are killed around the 50-min
mark, and split into chunks the symbolic suites individually exceed the
foreground cap. The proportionate substitute, used at step 8.0: **grep the
reach**. `[M]` zero files in `cp/moc/mc/derivations` import `numerics.operator`,
and only THREE import `orpheus.sn|diffusion|transport` at all — those three run
in **5.7 s**. Chunks that DO fit: light dirs + root test files **684 passed /
55 s**; `mc + moc` **140 passed / 8:23**; `derivations` files 1–35
**538 passed / 7:12**.

### ⭐⭐ Durable lessons from this stretch

1. ⭐⭐ **A field a CONSUMER BRANCHES ON is an input, not metadata.** Now
   `plan-authoring` §8. Step 8.0 was scoped as an enabler; `_AdjointOperator`
   branches on `domain is None` to decide whether to apply metrics, so
   populating it moved the realized white boundary's `.H` by **87 %** —
   and **0.0** on the specular mirror, whose metric cancels. An unmeasured
   change ships **ungated**: the boundary tier reddened nothing until
   consequence gates were written for it. ERR-076.
2. ⭐⭐ **Retiring a MESSAGE STRING: grep the SHORTEST distinctive fragment.**
   `[M]` the audit grepped `requires equal domains`; two gates pin
   `"equal domains"`. A longer pattern is *strictly less sensitive* than the
   consumer's pin, so it returns a confident empty answer. Now in
   `coding-standards`, with its twin: a retirement onto a SHARED helper moves
   the raise's provenance one frame out, silently widening any
   innermost-frame assertion.
3. ⭐⭐ **Routing a call site through the ALGEBRA raises its operand
   requirement from "has the verb" to "IS the type"** — and only duck-typed
   test doubles notice, so it surfaces as a broken test rather than a contract
   change. Fix the surrogate; never add a guard the type system now covers.
4. ⭐ **Count the CATCHERS, not the green rows.** `[M]` reverting step 8
   reddens **1 of its 8** new rows; the other seven gate the algebra beneath
   it. Eight green rows read as eight guards unless the docstring says
   otherwise (`vv-principles` #18).
5. ⭐ **A number without its configuration mis-states its own claim.** Two
   battery figures were measured on a 189-test slice and reported beside two
   from a 199-test slice; re-measured on ONE 212-test slice they were **17 and
   5**, not 7 and 3 — a 2.4× UNDER-statement, in the flattering direction.
6. ⭐ **A liveness check must not be able to match itself.** `[M]`
   `until ! pgrep -f "pytest tests/cp …"` matched the WAITER'S own command
   line, so two monitors never terminated and the run they watched was killed
   unobserved. Third instrument failure of the stretch, all anti-pattern #17,
   all three by different mechanisms (see §7h.4's table). Belongs in
   `vv-principles` #17 — **not yet added there**, because `SKILL.md` was on the
   forbidden-to-commit list; it is recorded here instead.
7. ⭐ **Measure the leaf, then measure it AT THE TIER THAT PAYS.** Composing
   the face action is **1.50×** on the leaf and **0.42 %** of a solve. The
   wide run's +1:44 was variance; a cache would have been machinery bought for
   nothing.

---

## ⏸ COMPACTION POINT #5 — track G resumes at §7d.3 (delete the table, migrate its tier)

Additive to #4 (→ #3 → #2). #4's step-7 items are HISTORY — the chain-table
gap is closed in place, the pins it said must flip HAVE flipped and were
re-scoped — but its gate costs and durable lessons stand. A session picking
up cold re-anchors from **§7h.6's landing banner + its §7d block + this
file + `git log`**, never a summary.

### State

**Branch `refactor/operator-strategy-layers`, tree clean** (only untracked
`scratch/`). The stretch after #4 landed THIRTEEN commits — verify each with
`git merge-base --is-ancestor`:

| commit | what |
|---|---|
| `e13313a8` | 7a — `PairedDeck` retires `SpatialWrap` (motion-typed deck tier; `RigidMotion.is_translation`) |
| `c3bb7341` | 7b — `_deck_kernel`: one body for every deck arrow; `Quadrature.ordinate_permutation`; `RigidMotion.linear_part`; periodic BOUND |
| `635c22e5` | 7c — 142 gate rows (test-architect + periodic metric row); V&V matrix |
| `6d00e917` | plan: realization-phase landing record + the ERR-074 stale-hazard refutation in place |
| `ae018628` | 7d.0 — four audit-found falsehoods fixed on sight; census recorded |
| `a17771db` | 7d.1 — certifications derive from the motion; injection sites re-posed; realize-level refusal pin |
| `f11a129f` | 7d.2 — pole seed hoisted to `_ensure_pole_mirror()` (mesh-stash) |
| `ba775c4d` | plan: 7d.0–7d.2 record + 7d.3 census |
| *(+ this checkpoint commit; earlier in the stretch: `72f5ce97` `b89ccd2d` `84d7861f` `4a87cc34` `811b6407` were #4's steps 8.0/8, already in its table)* | |

### ⛔⛔ READ FIRST — what a cold session would get wrong about 7d.3

1. ⛔⛔ **The table being alive is DELIBERATE, not an oversight.**
   `reflection_index` / `reflection_partners` /
   `_compute_sphere_reflection_partners` are production-DEAD (every
   production consumer migrated at 7b/7d.1/7d.2) but API-alive with their
   own green gates. Do NOT "fix" anything by rewiring production back onto
   the table; 7d.3 is deletion + test/doc migration ONLY.
2. ⛔ **The worklist file is UNTRACKED**:
   `scratch/audit_reflection_index_retirement.md`. If it has evaporated,
   the §7d block above carries the census summary and the re-grep is one
   command (`grep -rn "reflection_index\|reflection_partners" --include="*.py" --include="*.rst" orpheus/ tests/ docs/ derivations/`).
   `[M]` post-7d.2 census: **27 files / ~140 lines**, all test/doc tier +
   the definition tier.
3. ⛔ **`reflection_partners` is a DATACLASS FIELD** — deleting it changes
   `Quadrature`'s constructor. The factories PASS it (breaks loudly); the
   field has a default, so omission sites are safe; grep direct
   `Quadrature(` constructions (the docstring blesses them for
   restrict/pushforward wrappers) before deleting.
4. ⛔ **The `:294` pin trap**: `test_quadrature_directional.py` pins
   `pytest.raises(ValueError, match="no precomputed reflection partner")`.
   `BoundaryError ⊂ ValueError`, so a lazy re-point stays green against
   the WRONG refusal — pin the surviving refusal's own fragment (the
   7d-era wording is `"no specular pairing"` at the certification tier,
   `"cannot seed the r = 0 pole"` at the sweep tier).
5. **`TestTheSpecularArmInheritsTheRetiredTable`** (in
   `test_deck_kernel.py`) deletes outright per its own docstring — a
   one-time equivalence claim against the retired path. Do not weaken it
   into something else.
6. **The equation node** `quadrature-reflection-index`
   (`angular_quadrature.rst:239`) carries **23 `tests` edges** via the
   file-level verifies list at `tests/sn/primitives/test_quadrature.py:35`
   — the marker-migration surface. Re-pose the equation (the concept "the
   ordinate permutation a mirror induces" survives as
   `ordinate_permutation`) rather than silently dropping the edges.
7. **Rendered-vs-grep-only docs**: `docs/api/discrete_ordinates.rst:67` is
   the ONE Sphinx-rendered xref that dangles visibly on deletion; the
   other ~20 doc mentions are invisible to every build severity — **grep
   is the only gate** (the corpus-reference-integrity lesson).
8. **Tense-discriminate the ~140 lines**: many are docstring/prose
   mentions. Past-tense history ("read the table until 7d") STAYS;
   present-tense claims re-point; imperative instructions are MUST-FIX.
   Two intentional past-tense mentions were minted by 7d itself
   (`_ensure_pole_mirror`'s and `_specular_pairing`'s docstrings).
9. End of 7d.3 = the three-search audit + the shortest-fragment message
   grep + a full wide slice + sphinx `-W` + xref checker + pyright.

### Red baseline

`[M]` at 7c (`635c22e5`): wide slice = **7 failed / 6159 passed**, 17:30 —
the SAME pre-existing seven as #2/#3/#4 (2× `cart2d_*_principled_equiv` ·
3× `test_affine_carve_bit_identity` sha256 (#333) · 1×
`test_diamond spherical_inward` · 1× `WhiteXminPartial03GLSnapshot` (#33)).
After 7d.0–7d.2 the focused batteries reproduce only baseline reds
(geometry+sn/operators 1937/3 · sweep+solve+eigenvalue 938/4); the wide
slice at THIS checkpoint's HEAD was launched with this checkpoint — its
result lands in the line below.

`[M]` wide slice at checkpoint HEAD (`ba775c4d` + this checkpoint's
worktree): **7 failed / 6163 passed**, 21:41 — the SAME pre-existing
seven; the passing count grew 6159 → 6163 with 7d.1's realize-level pin
and the re-posed rows. The baseline has now been reproduced at FOUR
consecutive states of this stretch (#4's, 7c's, the focused 7d
batteries', and this one).

### Gate costs (unchanged from #4 where re-measured)

| gate | cost |
|---|---|
| wide slice | 17:30–19:59 (variance) |
| `tests/sn/{sweep,solve,eigenvalue}` `-m "not slow"` | ≈7:50 |
| `tests/geometry + tests/sn/operators` | ≈25 s |
| sphinx `-W` incremental | ≈1 min |
| pyright ratchet | ≈16 s |
| `tools/check_docstring_xrefs.py orpheus tests docs` | DEAD TARGETS 0 |

### ⭐⭐ Durable lessons of the 7a–7d.2 stretch

1. ⭐⭐ **Ordinates are DIRECTIONS — feed the matcher the element that acts
   on them.** `RigidMotion.permutes/preserves` match through `on_points`
   (the affine action); the wrap translation matched to `None` — a false
   "not a symmetry". Minted `RigidMotion.linear_part` (the E(d)→O(d)
   projection) and match through it. General form: when a point-set
   matcher is fed a set living in a quotient vocabulary, the acting
   element is the PROJECTED one. Found independently twice (probe +
   architect battery, 26/78 rows red under the affine spelling).
2. ⭐⭐ **A certification and its realization must read ONE object — and
   the injection surface must move WITH the read.** The certifications
   certified the precomputed table while realization derived from the
   motion: bit-identical by gate, structurally twinned. After the swap,
   five `monkeypatch.setitem(reflection_partners, …)` sites poisoned a
   table nothing read — the negatives would have gone silently vacuous
   (vv #17's harness-lies-first, at the fixture tier). Re-posed onto
   `ordinate_permutation`; the duck-typed fakes gained the new contract
   door.
3. ⭐⭐ **A mutation that makes a law INCOHERENT reds at realization —
   probe consumers by mutating the PROPERTY they read** (vv #18,
   mutate-inside-the-legal-class). `TestSpecMutationsPropagate`'s
   wholesale `geometry_map → identity` mutations began reddening at mesh
   construction (the kernel refuses an incoherent mirror — itself proof
   the realizer stopped reading the letter); the probes moved one
   property up to `permutes_ordinates`, the flag the consumers actually
   read.
4. ⭐ **A convention pin needs a fixture where the convention is
   OBSERVABLE.** The π-vs-π⁻¹ keystone on `product(4,4)` was degenerate
   (its local permutation is `arange`, same as the wrap's); moved to
   `product(4,8)`. And a mutation that reds only by RAISING is attributed
   to nothing — pair it with an in-range mutation that reds by VALUE.
5. ⭐ **`BoundaryError ⊂ ValueError` makes type-only pins vacuous across
   refusal migrations** — a re-pointed refusal gate must pin the
   surviving message's distinctive fragment.
6. ⭐ **The epistemic labels PAID**: the plan's own proposed means carried
   the affine/linear conflation, and because it was labelled *proposed,
   NOT verified*, the reader measured before building. Likewise §7
   reconciliation caught a hazard claim (ERR-074 "unhardened") whose fix
   had landed in a DIFFERENT campaign between writing and reading —
   refuted in place, no new clause needed (the existing §3/§7.3 caught
   it; logging a one-off would be the loop's failure mode (a)).
7. ⭐ **A count's SCOPE is part of its configuration** (§4 applied to
   censuses): "10 call sites" → "33 live refs" → "7 production + 44 test
   lines + 21 doc mentions" were three DIFFERENT greps, each correct for
   its pattern and scope. State what was counted, or the next reader
   plans against the wrong magnitude.

## ⏸ COMPACTION POINT #6 — step 7 CLOSED; track G resumes at G6.3b (the flagship metric gate)

Additive to #5 (→ #4 → #3 → #2). #5's 7d.3 items are HISTORY — the
retirement it scoped is fully landed — but its cold-session hazards about
the DELETED table now read as "do not resurrect", and its gate costs and
lessons stand. A session picking up cold re-anchors from the **status
header → this checkpoint → §7d.3's landing record → the plan's G6 rows
:1496 and :1905-1945** (the G6.3b design inputs), never a summary.

### State

**Branch `refactor/operator-strategy-layers`, HEAD `93e71172`, tree clean
except three deliberately-uncommitted files** (qa agent-memory ×2 and ⛔
`.claude/skills/vv-principles/SKILL.md` — the FORBIDDEN-to-commit file,
whose working delta now contains BOTH the user's edits AND qa's 35-line
anti-pattern **#22**; never commit it). The stretch after #5 landed SIX
commits — verify each with `git merge-base --is-ancestor`:

| commit | what |
|---|---|
| `c7ca338e` | 7d.3a — surviving test tier off the table, BY CLAIM CLASS (reference gates → `tests/_harness/references.py::mirror_partner_indices`; production-datum gates → production's exact expression); the sweep-tier refusal pin minted (`"cannot seed the r = 0 pole"`, previously UNPINNED) |
| `b5ac130e` | 7d.3b — the table DELETED (`reflection_index`, the `reflection_partners` FIELD, `_compute_sphere_reflection_partners`, GL1D closed form, `_resolve_axis_to_index`); Q4 gates re-posed; oracle class deleted per its own docstring; equation node re-posed `quadrature-ordinate-permutation`; 8 doc pages + `_errors.py` migrated. Net −97 lines |
| `fc5c7a67` | qa's eight prose-truth findings fixed; the specular-response independence claim honestly NARROWED with the measurement |
| `1b6beb86` | plan: the §7d.3 landing record, incl. the owed-wide-slice statement |
| `93e71172` | the G6.3b resume pointer's STALE MECHANISM TITLE corrected (plan-authoring §1 widened + surprise-log row) |
| *(+ this checkpoint commit)* | |

Issue **#334** filed (file-level `verifies` lists over-credit equation
coverage — qa informational, pre-existing).

### ⛔⛔ READ FIRST — what a cold session would get wrong about G6.3b

1. ⛔⛔ **G6.3b is the flagship metric gate ONLY.** Its original "absorb
   `AngularAverageOperator.apply_transpose`" half is DISCHARGED: that
   operator was retired at step 3b (2026-08-04, #330), factored into
   `IsotropicEmissionOperator @ PartialCurrentOperator`
   (`orpheus/sn/boundary/angular.py`) where each link has one honest
   transpose and the composite adjoint falls out by composition. Grep
   `AngularAverageOperator` — every hit is past-tense. Re-derive the
   gate's scope from the TREE (§7.2): read `angular.py`'s module
   docstring, plan rows :1496 (the vv #12 rationale — "without it the
   flagship gate can only run on a stand-in") and :1905-1945 (the
   factored-adjoint payoff + the G6.1 strictly-positive-metric gate as
   its precondition). Design anchors that EXIST: the periodic metric
   row (`test_sn_boundary_realizer.py::test_the_PERIODIC_metric_agrees_across_the_pair`
   — weights bit-exact, `.H`-vs-transpose at `nulp=4`, the FP-non-
   associativity mechanism documented) is the deck-tier sibling of what
   the Lambertian tier needs; ERR-067 ("the G-metric is a free parameter
   no eigenvalue gate can see") is the failure class the gate exists for.
2. ⛔⛔ **The wide slice is OWED as G6.3b's landing gate.** The last
   MEASURED wide slice is **7 failed / 6163 passed at `39b29b42`**
   (pre-7d.3, 21:41). At `b5ac130e`..HEAD it was launched twice and
   externally stopped twice — NOT measured; §7d.3's landing record
   carries the in-lieu evidence (touched batteries 4231/3 — exactly the
   three of the seven baseline reds reachable in those dirs; the other
   four (3× affine sha256 #333 · diamond spherical_inward) were NOT
   re-run this stretch) labelled as EXPECTATION. Run the full slice
   before or with G6.3b's landing; expect the same pre-existing seven.
3. ⛔ **The reflection table is GONE — do not resurrect.** Production and
   tests read `Quadrature.ordinate_permutation(motion)` (structurally
   the only source). The independent TEST reference is
   `tests/_harness/references.py::mirror_partner_indices` — its local
   `_AXIS_INDEX` is DELIBERATE (vv #22: it is the only surviving
   cross-check of the axis-letter→normal convention; importing
   `AXIS_NAMES` there would close it silently).
4. ⛔ **The killed-slice protocol** (this stretch's precedent): a
   background gate externally stopped twice is closed out HONESTLY —
   record not-measured, cite in-lieu evidence as expectation, assign the
   measurement forward — never relaunch a third time against a possible
   deliberate stop, and NEVER write a `[M]` line from an output that
   does not exist (a forward-written slice count was caught and
   stripped before commit this stretch).
5. **Memory surfaces are current and consistent**: MEMORY.md index line
   + topic file both say "resumes at G6.3b — the flagship metric gate"
   with the stale-title ⛔ and the owed-slice ⚠. The index is at
   ~19.5 KB nearing its 24.4 KB cap — a fuller per-owner distillation
   pass is QUEUED follow-up (only the track-G line was trimmed this
   stretch; other owners' lines need lossless-cut verification first).
6. **After G6.3b**: G6.5 (traversal onto the trace spaces), then track G
   is DONE and Q5.1 unblocks (the carve-out condition). The seven
   baseline reds are tasks #51/#33/#333's business — unchanged all
   stretch.

### Red baseline

`[M]` last measured wide slice: **7 failed / 6163 passed** at `39b29b42`
(21:41) — 2× `cart2d_*_principled_equiv` · 3× `test_affine_carve_bit_identity`
sha256 (#333) · 1× `test_diamond spherical_inward` · 1×
`WhiteXminPartial03GLSnapshot` (#33). This stretch: NOT re-measured wide
(two external stops — hazard 2 above); `[M]` touched batteries at
`b5ac130e` = 4231 passed / 3 failed (the three reachable baseline reds),
sweep/mms trio 29 passed / 7 xfailed (#326), qa's independent baseline
331 passed / 1 xfailed on the seven reference files.

### Gate costs

| gate | cost |
|---|---|
| wide slice | 17:30–21:41 measured historically; UNCONTENDED only — foreground batteries alongside stretch it past 60 min (this stretch, killed before completing) |
| tests/numerics + geometry + sn/operators + sn/primitives + acyclicity + coupled-pole | 5:33 (measured UNDER background-slice contention; expect less clean) |
| sphinx `-W` incremental | ≈1 min |
| pyright ratchet | ≈16 s |
| `tools/check_docstring_xrefs.py orpheus tests docs` | DEAD TARGETS 0 |

### ⭐⭐ Durable lessons of the 7d.3 stretch

1. ⭐⭐ **Independence has TWO axes — derivation and INPUT resolution**
   (qa-minted vv-principles **#22**). The specular-response cross-check
   kept independent derivations, but both sides came to consume ONE
   `SelfPairedDeck.mirror(axis).motion`: `[M]` an axis-letter swap left
   it all-green while 78 sibling gates red. Single-sourcing an input
   closes a convention cross-check SILENTLY; the fix was narrowing the
   claim + keeping the reference helper's axis map deliberately local.
2. ⭐⭐ **A retirement's test migration is claim-class TRIAGE, one
   retirement, two opposite correct answers**: reference gates must get
   an INDEPENDENT datum (re-pointing them at production's source is the
   rewire-demotion — every one becomes production-vs-itself), while
   production-datum gates must get production's EXACT expression (an
   independent lookalike would gate nothing production reads). `[M]` qa's
   in-class mutation: 21/21 reference gates red, none self-comparing.
3. ⭐ **The CONCEPT grep out-catches the symbol grep — grep every
   SPELLING of the concept, not the symbol's**: two rounds in one
   retirement ("reflection-index table" hyphenated prose in
   `_errors.py`; "reflection-partner map" ×3 INSIDE the gutted package).
   No build severity renders those modules; grep was the only gate, and
   the symbol grep alone would have shipped four present-tense
   falsehoods.
4. ⭐ **§1 governs every RESUME POINTER** (plan-authoring widened + the
   surprise-log row): a close-out copies the next step's title verbatim
   from the frozen task list — reconcile the pointer's named objects
   against the tree (one grep per symbol) before it ships. The
   propagation was same-day into THREE surfaces; a user question, not
   the protocol, caught it.
5. ⭐ **A `[M]` line is a TRANSCRIPTION, never a forecast**: results are
   written only from outputs that exist. The checkpoint protocol's
   "PENDING AT WRITE TIME" idiom is the honest placeholder; this stretch
   proved both halves (a forward-written slice count stripped before
   commit; the killed slice recorded as not-measured and assigned
   forward rather than papered over).
6. ⭐ **A symbol retirement's two-commit shape**: consumers first (tree
   green, the old path still self-gated), then deletion + its own gates
   + docs + label TOGETHER (the xref/verifies gates red between —
   §6b's call-site-set discipline applied to code, Sphinx labels, and
   `verifies` strings alike).
