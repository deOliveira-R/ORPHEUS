# Geometric-transformation audit — SPATIAL layer

**Scope**: `orpheus/geometry/` (minus `boundary/`), `orpheus/transport/mesh/`,
`orpheus/sn/mesh/`, `orpheus/sn/sweep/` (spatial/topological parts),
`orpheus/transport/spatial/`, plus the two shared L1 primitives the spatial
layer routes through (`orpheus/numerics/face_layout.py`,
`orpheus/sn/loss_representation/sweep_graph.py` — the sweep DAG).

**Branch** `refactor/operator-strategy-layers`, HEAD `b0a003b4` at open.
Tree clean except `.claude/skills/vv-principles/*` + untracked `scratch/`.
Line numbers current at that HEAD; re-derive via Nexus if drifted.

**Verdict up front (the honest category report).** The spatial layer's
transformations are almost entirely **NOT point-group elements in the
`symmetry.py` sense** — they are (a) *index-lattice order reversals*,
(b) *a face-role involution on a 2-element set*, (c) *sparse
axis-aligned normals encoded as a name→(axis, sign) bijection*, and
(d) *integrated Jacobians of a coordinate chart whose chart map is
never written down*. There are exactly **two** constructs in the
spatial layer that ARE genuine group representations, and only one of
them is used as a transformation. Details in the table + Q&A.

---

## The table

| file:line | the transformation, in math | how it is spelled TODAY | what KIND | EXACT? |
|---|---|---|---|---|
| `orpheus/numerics/face_layout.py:130-136` | `σ: {min,max} → {−1,+1}`, and its derived inverse | **(O)** two dicts, the second a comprehension-inverse of the first | not a transformation — an **encoding convention** (the outward normal `n̂_f = sign·ê_axis` stored sparsely) | bit-exact (int) |
| `orpheus/numerics/face_layout.py:209-235` `face_opposite` | `f ↦ f'` with `n̂_{f'} = −n̂_f`; an involution on face names | **(S)** sign negation inside the name bijection (`face_name(axis, -sign)`) | **point-group element** — this is `σ_axis` restricted to the face set, i.e. `Mirror(axis)` acting on codim-1 faces. NOT expressible in `symmetry.py` today (that machinery acts on `(N,3)` node arrays, not on named faces) | bit-exact (int sign flip + string lookup) |
| `orpheus/numerics/spaces/angular_trace_space.py:229-241` (`build_omega_dot_n`) | `Ω·n̂_f = sign_f · μ_{axis(f)}` | **(S)** scalar sign multiply against a column selected by axis index | coordinate-**projection**, not a transformation (the normal is never materialised as a vector) | evaluated (the μ column is a quadrature cosine); the sign is exact |
| `orpheus/sn/loss_representation/sweep_graph.py:572-573` | per-axis face-role swap `(in,out) = (0,1) ↦ (1,0)` under `s_a → −s_a` | **(P)** ternary on the sign, producing an int pair | **point-group element restricted to a 2-set** — the axis mirror `σ_a` acting on the cell's two `a`-faces. Same object as `face_opposite`, spelled independently | bit-exact (int) |
| `orpheus/sn/loss_representation/sweep_graph.py:577-580` | `x_a ↦ −x_a` on the cell index lattice: `local→global` map is `arange(n)` or `arange(n)[::-1]` | **(P)** `np.arange(n)[::-1]` | **point-group element, genuinely** — the reflection `σ_a` of the Cartesian cell lattice, realized as the index permutation `i ↦ n−1−i`. It IS the reduction of `Mirror(a)` to the discrete cell set; not expressible in `symmetry.py` (no lattice/permutation-on-cells realization there) | bit-exact (integer permutation) |
| `orpheus/sn/loss_representation/sweep_graph.py:1046-1090` (`_CellResidualTranspose` docstring + its use) | the adjoint runs the **mirror octant's** DAG (`−signs_eff`); `face_in(−o) = face_out(o)` and `levels(−o) = reverse(levels(o))` | **(O)** *the transformation is USED, not recomputed* — the caller hands the walk the `−signs` graph from the same cached family, and the unchanged walk realizes the transposed data flow | **point-group element used as a transformation** — the full inversion `−I` on the octant sign lattice `{±1}^d`. This is the ONE place in the spatial layer where a group element does real work rather than being open-coded | bit-exact (the addressing; the per-cell algebra is the kernel VJP, which is *not* μ-reversal) |
| `orpheus/transport/spatial/_ubld.py:96-148` (`octant_moment_frame_signs`) | `sign[o₀..o_{d−1}] = Π_a (s_a)^{o_a}` — the sweep⇄global spatial-moment frame map | **(S)** explicit sign vector, built by a nested Python loop over flat moment index | **point-group element, as a REPRESENTATION** — the diagonal (character) representation of `(Z₂)^d = ⟨σ_x,…,σ_d⟩` acting on the tensor-Legendre moment basis. Docstring names it "an INVOLUTION"; it is never named as a group element. NOT expressible in `symmetry.py` (that realizes `O(3)` on ℝ³ nodes, not a rep on a polynomial basis) | bit-exact (±1) |
| `orpheus/sn/mesh/augmented_mesh.py:1275-1278` | 1-D sweep sense: cell order `0…nx−1` vs `nx−1…0` | **(O)** `range()` vs `range(..., -1, -1)` | **ordering reversal** (the same `σ_x` as above, but as a loop-order flag — no permutation object exists) | bit-exact |
| `orpheus/sn/mesh/augmented_mesh.py:1398-1400` (`_iter_cartesian_visits`) | same | **(O)** ditto — 2nd independent spelling | ordering reversal | bit-exact |
| `orpheus/sn/mesh/augmented_mesh.py:1426-1431` (`_iter_spherical_visits`) | same **+** the inner/outer face-role swap | **(O)** `range` flip **+ (P)** `select_outer: bool` picking `face_area_outer` vs `face_area_inner` | ordering reversal + the 2-set face involution as a **boolean flag** | bit-exact |
| `orpheus/sn/mesh/augmented_mesh.py:1499-1504` (`_iter_cylindrical_visits`) | same | **(O)** + **(P)** — 4th independent spelling | ordering reversal + face involution as a bool | bit-exact |
| `orpheus/sn/sweep/scan.py:285-288` (`ordinate_scan_transpose`) | the transposed first-order recurrence = the same scan on the **reversed** sequence with shifted multipliers | **(P)** `a[1:][::-1]`, `v[::-1]`, result `[::-1]` | **ordering reversal** used as an honest algebraic identity (the transpose of a lower-bidiagonal solve is the upper one) — it is `σ_x` on the chain index, exact | bit-exact (a permutation of a float array) |
| `orpheus/sn/sweep/psi_half_angle_seed.py:125`, `:303` | inward radial march = outward march on reversed arrays | **(P)** `[:, ::-1]` / `[::-1]` on `Q`, `σ_t`, `dr` | **ordering reversal** (`σ_r` on the radial cell index), reused so the inward march shares the outward kernel | bit-exact |
| `orpheus/geometry/coord.py:46-95` | the 1-D chart Jacobian, integrated: `V = ∫J dr`, `A = J(r)` | **(C)** `match coord:` over three closed-form expressions | **coordinate chart change — but only its INTEGRATED Jacobian.** The chart map itself (`x = r cosθ`) is never written | evaluated (π, powers) |
| `orpheus/geometry/coord.py:100-135` | 2-D `(r,z)` volume `π Δ(r²) Δz` | **(C)** ditto, 2 cases | coordinate chart change (Jacobian only) | evaluated |
| `orpheus/geometry/factories.py:60-74` (`_subdivide_zone`) | inverse-Jacobian: `r_k = √(r₀²+ (k/n)Δ(r²))`, `∛(...)` | **(C)** closed form per coord system | coordinate chart change (the **inverse** of the `coord.py` Jacobian, derived independently — see Q2) | evaluated (`sqrt`/`cbrt`); volumes returned exactly to dodge the round-trip |
| `orpheus/geometry/factories.py:113-121` (`pwr_pin_2d`) | `r = √((x−p/2)² + (y−p/2)²)` — the polar radius of a translated origin | **(C)** + **affine translation** written inline as `- pitch/2` | **coordinate chart change + affine/translation**, both hand-written; the ONLY place in the spatial layer that evaluates a Cartesian→polar map | evaluated |
| `orpheus/geometry/structured_geometry.py:419` | Wigner–Seitz `r_cell = pitch/√π` (equal-area square→disc) | **(C)** one-line formula | **not a transformation** — an area-preserving *model substitution* (a different domain, not a map of the same one) | evaluated |
| `orpheus/geometry/reduced_operator.py:688-690` (sphere) | `α_{n+½} = α_{n−½} − w_n μ_n` | **(A)** cumulative Python loop | **connection coefficient of the chart** (module docstring L1-8 names it: "the same connection-coefficient operator on SO(3) viewed in two coordinate charts"). NOT a point-group element — it is a *Christoffel-like* term, i.e. the derivative of the frame rotation, not a finite group element | evaluated (products of quadrature data); pinned by `assert abs(alpha[N]) < 1e-12` |
| `orpheus/geometry/reduced_operator.py:783-785` (cylinder) | `α^{(p)}_{m+½} = α^{(p)}_{m−½} − w_m η_m` | **(A)** the SAME loop body, restated per level | ditto — a **verbatim twin** of the sphere recursion 95 lines above | evaluated |
| `orpheus/geometry/reduced_operator.py:808` and `orpheus/sn/sweep/pole_angular_closure.py:598` | `sinθ_p = √(1−ξ_p²)` — the level's azimuthal extent | **(C)** the same closed form, **twice, in two modules** | coordinate chart relation (polar↔azimuthal), restated | evaluated |
| `orpheus/sn/mesh/augmented_mesh.py:1552-1556` (`_setup_cartesian`) | `g_a = |μ_a|/Δa` | **(S)** `np.abs(...)` on the axis cosine | **not a transformation** — the `abs` is the *magnitude* of the streaming coefficient (the direction is carried by the octant label), plain arithmetic | evaluated |
| `orpheus/transport/mesh/material_xs_field.py:684-687` | `(N_cells, ng) → (ng, *spatial)` | **(O)** `.T.reshape` | **not a transformation** — a memory-layout transpose (axis reorder of a data buffer), no geometry | bit-exact |

### Table, continued (found after the first pass)

| file:line | the transformation, in math | how it is spelled TODAY | what KIND | EXACT? |
|---|---|---|---|---|
| `orpheus/sn/loss_representation/__init__.py:3294` | the r=0 pole continuity `ψ(0,+μ) = ψ(0,−μ)`: the outward sweep's centreline seed IS the inward sweep's centreline outflow at the **mirror ordinate** | **(P)** `outflow_at_inner.T[quad.reflection_index("x")]` — an ordinate PERMUTATION array indexed straight into the buffer | **point-group element, named, permutation-realized** — `Mirror("x")` under the tree's canonical `(μ,0,0)` embedding. This is the **deck transformation of the spatial quotient**: the 1-D radial mesh is the 3-D ball modulo the rotation group, and this is the identification at its singular point | bit-exact (integer index array); the permutation itself was built by nearest-neighbour matching at quadrature construction (evaluated once) |
| `orpheus/sn/loss_representation/__init__.py:3470`, `:4189`, `:4593` | same map, used in the two adjoint walks + the per-ordinate curvilinear march | **(P)** `mirror = quad.reflection_index("x")` — three more independent fetch sites, all with the axis literal `"x"` hard-coded | ditto | bit-exact |
| `orpheus/sn/loss_representation/__init__.py:721-733` (`_inflow_faces`) | octant sign → the domain face it enters: `s_a ≥ 0 ↦ a`-min | **(S)** `face_name(a, -1 if s >= 0 else +1)` | point-group element on the face 2-set again (5th spelling of the same σ_a), but routed through the ONE bijection | bit-exact |
| `orpheus/sn/loss_representation/__init__.py:736-747` (`_outflow_faces`) | `inflow ∘ σ` | **(O)** *the transformation is USED*: `face_opposite(face) for face in _inflow_faces(...)` | the involution applied as a map instead of a second transcription — the docstring at `:740-746` says exactly this ("that sentence is the implementation rather than a comment beside a second transcription") | bit-exact |
| `orpheus/sn/loss_representation/__init__.py:750-787` (`_reverse_octant_traversal`) | octant label `o ↦ −o` on `{±1}^d` (grazing `0` first promoted to `+1`) | **(S)** `-(+1 if s == 0 else s)` inside a generator | **point-group element as DATA** — the inversion `−I` on the octant sign lattice; the reverse walk then just *selects a different graph*. NOT an involution on grazing labels (`0 ↦ −1 ↦ +1`), by design | bit-exact |
| `orpheus/sn/loss_representation/__init__.py:999` | grazing map `0 ↦ +1` | **(S)** `+1 if s == 0 else s` | **not a transformation** — a choice of orbit representative for a degenerate (zero-streaming) octant | bit-exact |
| `orpheus/sn/loss_representation/__init__.py:1875` / `:1914` | the domain in-edge / out-edge face index swaps `0 ↔ n_a` under `s_a → −s_a` | **(P)** two ternaries on the sign, mirrored between the two functions | the σ_a face involution on the `{0, n_a}` endpoint pair — 6th/7th independent spelling | bit-exact |
| `orpheus/sn/loss_representation/__init__.py:2758-2770` (`_reverse_traversal`) | reverse-mode program order: legs reversed AND each leg's cell chain reversed | **(P)** `leg.cells[::-1]` inside `reversed(legs)` | **ordering reversal** used as the honest adjoint of a program (docstring names the pole-handoff edge flip it induces) | bit-exact |
| `orpheus/sn/sweep/scan.py:283-289` (`ordinate_scan_transpose`) | transposed bidiagonal solve = forward scan on the reversed, shifted chain | **(P)** `a[1:][::-1]`, `v[::-1]`, `[...][::-1]` | ordering reversal, algebraically exact (transpose of a triangular solve) | bit-exact |
| `orpheus/sn/sweep/scan.py:345-363` (`_x_scan_faces`) | reverse the chain, scan, reverse back | **(P)** `alpha[..., ::-1]` etc., driven by a **`x_reverse: bool` parameter** | ordering reversal, and the one place it is a **boolean flag parameter** rather than data or an object | bit-exact |
| `orpheus/sn/sweep/psi_half_angle_seed.py:125`, `:303` | inward radial march via the outward kernel on reversed arrays | **(P)** `[:, ::-1]` on `Q`/`σ_t`/`dr` | ordering reversal | bit-exact |
| `orpheus/geometry/structured_geometry.py:127-134` + `:44-61` | CYL/SPH carry ONE endpoint because the domain is an **orbit space** — the centreline is the singular point of the quotient | **(O)** an int lookup table `{"SLB": 2, "CYL": 1, "SPH": 1}` + a module docstring that says "orbit-space classification" and "orbit-space rank" | **quotient by a group** — and the QUOTIENT is named (in orbit-space language) while the GROUP never is (no `SO(3)`, no `O(2)`, no `Mirror`) | bit-exact (a table) |
| `orpheus/geometry/structured_geometry.py:430-488` (`pwr_slab_half_cell`) | the half-cell is the unit cell modulo the mirror through the fuel centreline | **(O)** `origin = 0` + a DEFAULT `BC("reflective")` at both ends; the docstring at `:441-443` says "Encodes the conceptual symmetry … about the fuel centreline" | **quotient by `Mirror`** — realized entirely as a boundary condition; the group is not named and no map exists | n/a |
| `orpheus/transport/mesh/axis.py:118-147` (`FaceLabel.face_name`) | `(axis_index, endpoint) ↦ "{axis}{min|max}"`, with `"outer" ↦ "max"` | **(O)** a `_ENDPOINT_SUFFIX` dict lookup + f-string | **not a transformation** — a naming crosswalk; but note it silently identifies the radial `outer` endpoint with a Cartesian `max` face, which is what lets the curvilinear mesh reuse the Cartesian sign convention | bit-exact |
| `orpheus/sn/sweep/pole_angular_closure.py:1012-1030` (`_edge_seed_stencil`) | linear extrapolation in μ to the level's starting-direction edge | **(A)** `t = (μ_start − μ_{m0})/(μ_{m1} − μ_{m0})`, with an `argsort` + a `1e-14` distinctness scan | **not a transformation** — an interpolation stencil. Listed because it is the *alternative* to the `reflection_index` pole handoff on non-carrying levels: the same physical continuity, one route a group element, the other a 2-point extrapolation | evaluated (`argsort` order is data-dependent; the `1e-14` guard is a tolerance) |
| `orpheus/transport/spatial/diamond.py:723-738` (`reflect_scan_coefficients`) | `α = −1` in the face recurrence (`ψ_out = 2ψ̄ − ψ_in`) | **(S)** a `−1` multiplier | **NOT a transformation** — a false positive from the word "reflect". This is DD's diamond *closure*, a scalar affine relation between two face values; nothing geometric moves | bit-exact (`−(1−½)/½ == −1` exactly) |
| `orpheus/transport/spatial/_ubld.py:89-91` (`_GRAD_1D`, `_FOUT_1D`) | literal 2×2 matrices | **(M)** `np.array([[0,0],[-2,0]])`, `np.array([[1,1],[1,1]])` | **not transformations** — LD element matrices (a gradient and a face outer product) on the tensor-Legendre basis. The ONLY hand-written small matrices in the whole spatial scope; see Q6 | bit-exact (integers as floats) |

---

## Runtime verification of the structural claims

Probe script (read-only, written for this audit):
`/Users/rodrigo/git/nuclear/ORPHEUS/scratch/_geom_probe_spatial.py`.
Run with `.venv/bin/python -O`. Results:

**P1 — `face_opposite` IS `Mirror(axis)` on the face set.** Involution
over all six `FACE_NAMES`: `True`. And building the six outward normals,
applying `symmetry._reflections(axis)`, and mapping back reproduces
`face_opposite` on that axis + identity on the other two — **exactly**,
for `σ_x`, `σ_y`, `σ_z`. So the codebase's face involution and
`symmetry.py`'s reflection matrices are literally the same group element
on two different realizations.

**P2 — the mirror-octant DAG claim, precisely.** For every octant at
shapes `(5,)`, `(4,3)`, `(3,3,2)`: `face_in(−o) == face_out(o)` and
`face_out(−o) == face_in(o)` — **True everywhere**. `levels(−o) ==
reverse(levels(o))` — **True as SETS at every d**, and **elementwise only
at d = 1**. At d ≥ 2 the within-level cell ORDER differs (both graphs
order C-major over their own local lattice). This does not weaken the
transpose: cells within an anti-hyperplane level are mutually
independent, and `walk_full` uses the SAME `cell_idx` arrays for the
gather and the scatter, so within-level order is immaterial. Worth
saying explicitly, because the docstring's "the mirror graph's levels
ARE the forward's levels reversed" is exact at level granularity and
*not* elementwise past 1-D.

**P3 — `octant_moment_frame_signs` IS the character representation of
`(Z₂)^d`.** For `d ∈ {1,2,3}` and every one of the `2^d` octants, the
vector equals `χ_o(g) = Π_a s_a^{o_a}` computed independently; the map
`g ↦ χ(g)` is a **group homomorphism** (`χ(g)·χ(h) == χ(gh)` for all
pairs, exactly); every element is an **involution**. `per_axis == 1`
returns `None` (DD/Step). So this is not "an involution" — it is the
full 1-dimensional-character rep of the coordinate-mirror group,
restricted to the tensor-Legendre moment basis.

**P4 — `reflection_index("x")` IS `Mirror("x")`, exactly.** On
`gauss_legendre(8)`, `product(4,4)`, `level_symmetric(6)`, `lebedev(9)`:
`max|σ_x(nodes) − nodes[perm]| = 0.000e+00` in every case; the index map
is a bijection, an involution, and weight-preserving. So the spatial
pole handoff is a genuine, exactly-realized group action on the ordinate
set.

**P5/P6 — the category boundary, measured.** `symmetry._orbit_closure`
on a *centred* 1-D cell lattice (`x ∈ [−1,1]`, 6 cells, weights =
volumes) returns a certificate whose permutation is exactly
`[5 4 3 2 1 0]` — **the `np.arange(n)[::-1]` the sweep graph writes by
hand**. On the *production* geometry (`x ∈ [0,2]`, the mesh
`Mesh1D.from_geometry` actually builds, `origin=0`) the certificate is
`None`, and `SubgroupOfO3.Mirror('x').is_invariant(mesh.volume_measure)`
is `False` for both the slab-at-origin-0 and the sphere. **The missing
ingredient is the translation, not the reflection.** This is the single
most decidable finding of the audit — see Q&A below.


---

## The seven questions

### Q1 — sweep-direction reversal: where, and is it ever a transformation?

It is spelled **eleven** times, in four distinct idioms, and it is a
transformation-as-an-object in exactly **two** of them.

*Loop-order flag (never an object):*

1. `orpheus/sn/mesh/augmented_mesh.py:1275-1278` — `dag_walk_cell_indices`
2. `orpheus/sn/mesh/augmented_mesh.py:1398-1400` — `_iter_cartesian_visits`
3. `orpheus/sn/mesh/augmented_mesh.py:1426-1431` — `_iter_spherical_visits`
   (+ `select_outer: bool` face swap)
4. `orpheus/sn/mesh/augmented_mesh.py:1499-1504` — `_iter_cylindrical_visits`
   (+ `select_outer: bool`)
   — all four are literally `range(nx) if μ ≥ 0 else range(nx-1, -1, -1)`.
5. `orpheus/sn/sweep/scan.py:345-363` — `_x_scan_faces(..., x_reverse: bool)`:
   reverse in, scan, reverse out. The one **boolean flag parameter**.

*Face-role / edge-index swap under the sign (a 2-set permutation, written
as a ternary):*

6. `orpheus/sn/loss_representation/sweep_graph.py:572-573` — `face_in`/`face_out`
7. `orpheus/sn/loss_representation/__init__.py:1875` — in-edge index `0 ↔ n_a`
8. `orpheus/sn/loss_representation/__init__.py:1914` — out-edge index `n_a ↔ 0`
9. `orpheus/sn/loss_representation/__init__.py:721-733` — `_inflow_faces`
   (this one at least routes through the `face_name` bijection)

*An actual index-permutation object:*

10. `orpheus/sn/loss_representation/sweep_graph.py:577-580` —
    `axis_map = [np.arange(n) if s >= 0 else np.arange(n)[::-1] ...]`.
    This is the reversal AS DATA: a per-axis `local → global` index array
    that every downstream level index goes through. It is the reflection
    `σ_a` of the cell lattice, materialized as a permutation.
11. `orpheus/sn/sweep/scan.py:283-289`, `orpheus/sn/loss_representation/__init__.py:2758-2770`,
    `orpheus/sn/sweep/psi_half_angle_seed.py:125,:303` — `[::-1]` slices used
    as the honest algebraic reversal of a recurrence/program.

**Represented as a transformation rather than a flag — yes, twice, and
deliberately.**

- `orpheus/sn/loss_representation/__init__.py:750-787` (`_reverse_octant_traversal`)
  + `orpheus/sn/loss_representation/sweep_graph.py:1046-1090`
  (`_CellResidualTranspose`): the adjoint walk does not reverse anything.
  It maps the octant label `o ↦ −o`, **looks up the mirror octant's graph
  in the same cached family**, and runs the UNCHANGED forward walk. The
  docstring is explicit: "the whole reversal is DATA". This is the one
  place the reversal is a group element doing work.
- `orpheus/sn/loss_representation/__init__.py:736-747` (`_outflow_faces`)
  is `face_opposite ∘ _inflow_faces` — the involution applied instead of
  a second transcription (the twin it replaced is described in its own
  docstring at `:740-746`).

Neither is named as a group element, and neither is connected to
`symmetry.py`.

### Q2 — curvilinear (r, θ) ↔ Cartesian: one place or restated? Jacobian derived or independent?

**The chart map itself is never written down anywhere in the spatial
layer.** There is no `x = r cos θ`. `Mesh1D` (`orpheus/geometry/mesh.py:210-534`)
stores `edges` that "are radii for cylindrical/spherical" and a
`CoordSystem` tag; nothing converts. The only place a Cartesian→polar
map is *evaluated* is
`orpheus/geometry/factories.py:117` (`r = √((x−p/2)² + (y−p/2)²)`) — the
2-D pin-cell material painter, and even there only the radial component,
with the origin shift `−pitch/2` hand-written.

What DOES exist is the chart's **integrated Jacobian**, in one canonical
home: `orpheus/geometry/coord.py` — `compute_volumes_1d` (`:46-69`),
`compute_areas_1d` (`:72-95`), `compute_volumes_2d` (`:100-135`). The
module docstring claims it is "the **single point** where
coordinate-system dependence lives" (`coord.py:3-4`). That claim is
**almost** true. Three exceptions:

1. `orpheus/geometry/factories.py:60-74` (`_subdivide_zone`) writes the
   **inverse** Jacobian independently (`r_k = √(r₀² + f·Δ(r²))`,
   `cbrt(...)`) and returns exact volumes rather than round-tripping. Its
   docstring explains why (`sqrt(x)**2 != x`, ERR-020 drift) — a
   *deliberate* second derivation with a stated reason, not drift.
2. `orpheus/geometry/structured_geometry.py:419` — `pitch/√π`, an
   area-equivalence, in the geometry layer rather than `coord.py`.
3. `sinθ_p = √(1 − ξ_p²)` (the level's azimuthal half-extent) is written
   **twice**: `orpheus/geometry/reduced_operator.py:808` and
   `orpheus/sn/sweep/pole_angular_closure.py:598`.

The **connection coefficients** (the part of the chart change that is not
a volume factor) live in `orpheus/geometry/reduced_operator.py`, whose
module docstring at `:1-12` gives the sharpest framing anywhere in the
tree: the spherical `(1−μ²)/r ∂_μ` and cylindrical `−(1/r)∂_φ(ξ·)` terms
are "**the same connection-coefficient operator on SO(3) viewed in two
coordinate charts**". The concept is named. The *implementation* is two
verbatim-twin accumulation loops (`:688-690` sphere, `:783-785`
cylinder) with the `assert abs(alpha[N]) < 1e-12` closure check on the
sphere arm only.

Coordinate-system identity is carried in **three** parallel encodings:
`CoordSystem` enum (`geometry/coord.py:36`), `AxisCoord` enum
(`transport/mesh/axis.py`, crosswalked at `:474-479`), and a bare string
`SNMesh.curvature ∈ {None, "spherical", "cylindrical"}`
(`sn/mesh/augmented_mesh.py:327,:334,:1559`) that the loss
representation re-normalizes to a FOURTH vocabulary including
`"cartesian"` (`sn/loss_representation/__init__.py:3088-3091`).

### Q3 — is there a symmetry-reduced mesh, and is its group named?

**Yes, three of them; the quotient is named once (in orbit-space
language), the GROUP is never named, and no deck transformation object
exists.**

1. **CYL / SPH — the radial mesh is an orbit space.**
   `orpheus/geometry/structured_geometry.py:127-134` gives CYL and SPH
   **one** BC endpoint, and `:44-61` explains why: "The endpoint count
   comes from the **orbit-space classification** of the geometry's
   billiard table … CYL and SPH have one outer surface plus an implicit
   centreline reflection (**orbit-space rank 1**)." Restated verbatim in
   `docs/theory/foundations/structured_geometry.rst:118-121`. The word
   "orbit-space" appears; `SO(3)`, `O(2)`, `Z₂` never do.
   `orpheus/transport/mesh/axis.py:164-168` completes the picture from
   the other side: "the pole at r = 0 is intentionally NOT an endpoint —
   it carries an angular-closure regularity condition, not a BC trace
   law." That is the singular set of the quotient, described without the
   word.
   The **deck transformation of this quotient exists and is exact**:
   `ψ(0,+μ) = ψ(0,−μ)`, realized as
   `outflow_at_inner.T[quad.reflection_index("x")]`
   (`orpheus/sn/loss_representation/__init__.py:3294`, and three more
   fetch sites). Measured to BE `Mirror("x")` to 0 ULP (P4).
2. **`pwr_slab_half_cell`** (`orpheus/geometry/structured_geometry.py:430-488`):
   "Encodes the conceptual symmetry of a square PWR unit cell about the
   fuel centreline. The geometry starts at x = 0 (the symmetry plane…)"
   (`:441-443`). The quotient by `Mirror` is realized purely as
   `origin = 0` plus a default `BC("reflective")`. No group, no map, no
   fundamental-domain type.
3. **`wigner_seitz_pin_cell`** (`:366-428`) is NOT a symmetry reduction
   and should not be mistaken for one — it replaces a square by an
   equal-area disc (a different domain), and its `BC("white")` models a
   lattice, not a mirror. Its docstring calls it "the unit-cell symmetry
   assumption" (`:394`), which over-claims: the periodic-lattice
   *translation* quotient is real, the square→disc step is a
   substitution.

There is **no** octant/wedge/sector mesh anywhere in the spatial layer —
`grep -i 'wedge|sector|quarter|one-eighth'` returns nothing in scope.

### Q4 — face normals: computed geometrically, or assigned by convention?

**Assigned by convention, from the face NAME — but single-sourced
through a bijection, and never materialized as a vector.**

- The convention is one dict: `_FACE_SUFFIX_SIGN = {"min": -1, "max": +1}`
  at `orpheus/numerics/face_layout.py:130`, with its inverse **derived**
  at `:134-136` ("the two directions of one bijection must not be able to
  disagree").
- `face_normal(face) -> (axis, sign)` (`:179-206`) returns the outward
  normal `n̂_f = sign · ê_axis` **stored sparsely**; `face_name` (`:139-176`)
  is the forward direction. `FACE_NAMES` (`:241-245`) is derived from the
  product, not listed.
- The one consumer that needs `Ω·n̂` computes it as `sign * mu[axis]`,
  never as a dot product:
  `orpheus/numerics/spaces/angular_trace_space.py:229-241`.
- The comment block at `face_layout.py:106-123` records that before
  campaign phase **B3.4c** this convention was transcribed at **five**
  sites (a `_FACE_NORMALS` table in the trace space, two verbatim twins
  in the SN boundary layer, one in the sweep schedule, and a
  `face_name.endswith("max")` reverse-parse in the method registry). The
  collapse is recent and, from the greps here, complete inside the
  spatial scope.
- Nothing derives a normal from mesh **geometry** (edges, vertices,
  cross products). The mesh is axis-aligned by construction, so the name
  carries all the orientation there is. `FaceLabel.face_name`
  (`orpheus/transport/mesh/axis.py:118-147`) renders the radial `"outer"`
  endpoint as `"max"`, which is what lets a curvilinear mesh reuse the
  Cartesian sign convention unchanged.

### Q5 — does the cell ordering / DAG ever encode a reflection or rotation?

**Reflection: yes, twice, and one of them is the audit's best find.
Rotation: never.**

- **Inside one graph**, `orpheus/sn/loss_representation/sweep_graph.py:577-580`
  builds the per-axis reversal as an index permutation and `:572-573`
  swaps the face roles — i.e. the octant `o`'s graph IS the `(+1,…,+1)`
  graph composed with `Π_a σ_a^{[s_a<0]}`. That composition is not
  written; each graph is rebuilt from scratch per octant
  (`_graphs_for_shape`, `:1127-1143`).
- **Across graphs**, the adjoint DOES reuse one octant's order for
  another by symmetry: `_CellResidualTranspose`
  (`sweep_graph.py:1046-1090`) + `_reverse_octant_traversal`
  (`loss_representation/__init__.py:750-787`) run the reverse walk on the
  `−o` graph, relying on `face_in(−o) == face_out(o)` and
  `levels(−o) = reverse(levels(o))`. Verified (P2): the face-role
  identity holds exactly at d = 1, 2, 3; the level reversal holds as
  SETS at all d and elementwise only at d = 1 (immaterial — within-level
  cells are independent and the gather/scatter share one index array).
- The Gauss-Seidel schedule (`sweep_schedule.py:113-169`) orders octant
  GROUPS but derives each group's outgoing faces from the label's own
  signs (`_outgoing_faces`, `:242-255`) — no symmetry reuse, and
  deliberately so (`:141-150` explains that a diagonal cubature makes a
  face shared by ≥2 groups).
- No rotation appears anywhere: the mesh lattice is axis-aligned, so the
  only lattice automorphisms the code exercises are the `2^d` sign
  flips. The axis **permutations** that would complete `O_h` (the `S₃`
  half of `_octahedral_ops`) are never applied to a mesh.

### Q6 — any hand-built 2×2 / 3×3 transformation matrix?

**No — not one, in the entire spatial scope.** The only literal small
matrices are `orpheus/transport/spatial/_ubld.py:89-91`
(`_GRAD_1D = [[0,0],[-2,0]]`, `_FOUT_1D = [[1,1],[1,1]]`) and they are LD
**element** matrices (a gradient operator and a face outer product on the
tensor-Legendre basis), built up to higher `d` with `_batched_kron`
(`:164-302`) — the tensor-product structure of the basis, not a change of
frame. `np.diag` appears only in `material_xs_field.py` (group-transfer
matrices).

Every hand-built transformation matrix in the repo lives in the ANGULAR
module: `symmetry._reflections` / `_rotation_z` / `_vertical_mirrors` /
`_octahedral_ops` / `_icosahedral_ops` / `_rotation_about_axis`.
The spatial layer's transformations are all permutations and sign
vectors — which is exactly why none of them can currently be spelled in
`symmetry.py`'s `(3,3)`-matrix vocabulary.

### Q7 — hard-coded axis / index literals where the axis is a parameter?

Sorted by how much the literal is actually doing:

**Load-bearing and genuinely parametric (the strongest cases):**

- `quad.reflection_index("x")` at
  `orpheus/sn/loss_representation/__init__.py:3294`, `:3470`, `:4189`,
  `:4593` — four independent sites, all with `"x"` inline. The axis here
  is "the radial axis", which is always axis 0 on a 1-D curvilinear mesh,
  so the literal is *true*; it is nonetheless the mesh's radial axis, not
  the letter x, and the method itself accepts an int (`directional.py:405-408`
  calls the string form "the legacy SN slab tag … back-compat for the
  unmigrated sweep paths").
- The 1-D curvilinear face names `"xmin"` / `"xmax"` are inline at
  `loss_representation/__init__.py:3124-3126`, `:3356-3376`,
  `:3475-3535`, `:3729-3731`, `:3907-3908`; `sn/operators/boundary.py:808,842,915`;
  `sn/acceleration/dsa.py:232,687-690`. The 2-D/3-D paths derive their
  faces (`_inflow_faces`, `_outgoing_faces`); only the 1-D radial paths
  hard-code. Explicitly acknowledged at
  `sn/mesh/augmented_mesh.py:166-169` ("keeping `sn_mesh.bc["xmin"] ==
  "vacuum"` style comparisons … a solid sphere/cylinder has only ONE
  entry (`"xmax"`)").
- `sn/solver.py:1463` documents a paired literal/`face == "xmin"` compare.

**Positional-by-axis literals that are honest sugar** (documented as
such, not defects): `spatial_shape[0]` at
`transport/mesh/material_mesh.py:210-211` (`nx` sugar) and
`sweep_graph.py:455-459` (`nx`/`ny` 2-D compat accessors, marked "retire
with the d-generic orchestration"); `axes[0].coord` at
`transport/mesh/axis.py:474-479` (guarded — the function refuses a
multi-axis tuple containing a curvilinear axis).

**Where the literal has already been eliminated** (worth recording as
the counter-examples): `AXIS_NAMES` (`face_layout.py:100`) is the single
axis↔name crosswalk; `_setup_cartesian` builds the streaming tuple over
`range(ndim)` with "NO phantom axis"
(`augmented_mesh.py:1545-1556`); `det = d − 1` is called out as "the
SINGLE SOURCE OF TRUTH for the free/determined axis partition … change
`det` only" (`sweep_graph.py:637-642`).

---

## Gaps / adjacencies (things expected but not found)

- **No affine/translation type anywhere.** `origin` is a bare `float`
  parameter on `Mesh1D.from_geometry` (`geometry/mesh.py:405`), the pin
  centre is `pitch/2` inline (`factories.py:117`), and the half-cell's
  symmetry plane is "wherever x = 0 happens to be". This is the gap that
  makes P5/P6 fail.
- **Dangling doc xrefs**: `orpheus/geometry/reduced_operator.py:16-17`,
  `:653`, `:747` and `docs/theory/foundations/structured_geometry.rst:395-396`
  carry `:class:` / `:meth:` references to `SNMesh._setup_spherical` /
  `SNMesh._setup_cylindrical`, which **no longer exist** (grep: only
  `CPMesh._setup_*` in `orpheus/cp/solver.py:174,187` survives, an
  unrelated namesake). `structured_geometry.rst:414` itself says they
  "no longer exist" — the prose is right and the xrefs above it dangle.
  Python-domain xrefs render as plain text with no `-W` warning, so the
  Sphinx gate cannot see this.
- **Adjacency, out of scope**: `chord_half_lengths(radii, y)` in the CP
  solver (`orpheus/cp/solver.py:170-172`) is the same `√(r² − y²)` chart
  relation the spatial layer never writes. Another agent's territory;
  flagged only so it is not missed at the seam.

---

## The 3 findings that most surprised me

### 1. The category boundary is the TRANSLATION, not the reflection

I expected the verdict "spatial mirrors are not expressible in
`symmetry.py` because that module speaks `(N,3)` node arrays while the
mesh speaks index lattices." **Measured, that is false.**
`symmetry._orbit_closure` fed a 1-D cell lattice (cell centres embedded
as `(x,0,0)`, weights = cell volumes) under `_reflections("x")` returns a
valid `OrbitCertificate` whose permutation is **exactly `[5 4 3 2 1 0]`**
— the very `np.arange(n)[::-1]` that `sweep_graph.py:577-580` writes by
hand. `SubgroupOfO3.Mirror('x').is_invariant(mesh.volume_measure)`
dispatches cleanly through the 1-D arm and answers. The machinery
already produces the spatial layer's permutation.

What kills it is one thing: **production meshes are not centred on the
origin.** `Mesh1D.from_geometry` starts at `origin = 0.0`
(`geometry/mesh.py:405,:490,:516`), so the production slab lives on
`[0, L]` and the sphere on `[0, R]`. Measured:
`Mirror('x').is_invariant` = `False` for both, `True` only for a
hand-centred `[-1, 1]` slab. The spatial layer's mirrors are reflections
about the **domain's** mid-plane or about `r = 0` in a *radial*
coordinate — i.e. elements of the Euclidean group
`E(d) = O(d) ⋉ ℝ^d`, not of `O(3)`. The affine part is the whole gap,
and the codebase has no type for it anywhere: `origin` is a bare
`float`, the pin centre is `- pitch/2` inline (`factories.py:117`), and
the half-cell's symmetry plane is "wherever x = 0 happens to be"
(`structured_geometry.py:441-443`).

(The sphere/cylinder case has a second, sharper reading: `r ∈ [0,R]` is
the *fundamental domain* of the quotient, and a fundamental domain is
never `G`-invariant. `is_invariant == False` is the CORRECT answer, and
it says the spatial layer needs a **quotient/fundamental-domain**
vocabulary, not an invariance predicate.)

### 2. `octant_moment_frame_signs` is a full group representation, and nothing says so

I went in expecting a sign-flip hack. It measures as the exact
one-dimensional **character representation** of `(Z₂)^d` on the
tensor-Legendre spatial-moment basis: `χ_o(g) = Π_a s_a^{o_a}`, verified
against an independent construction and verified to be a **group
homomorphism** (`χ(g)·χ(h) == χ(g·h)` for all `2^d × 2^d` pairs at
d = 1, 2, 3, exactly) with every element an involution
(`orpheus/transport/spatial/_ubld.py:96-148`).

Both the docstring (`:124` "The map is an INVOLUTION") and the theory
page (`docs/theory/methods/sn/cartesian_multid.rst:2328,:2344`) say
"involution" — true of each element, but it undersells the object by a
whole structural level: it is not *an* involution, it is the rep of the
coordinate-mirror group. And this object was born out of a **bug**
(ERR-061, the diffusion-limit root cause, narrated at
`cartesian_multid.rst:2285-2330`: "the forward ordinate had
`ψ̂ = +0.048` but the backward had `−0.028` — opposite signs, the smoking
gun"). A group-theoretic framing — "the moment basis carries a
non-trivial rep of the octant group; ordinate sums must be taken in a
single frame" — would have made that bug unspellable rather than
findable-by-debugging.

### 3. The sweep-direction reversal is spelled eleven times, and the one place it is done RIGHT is the newest

I expected either "a flag everywhere" or "an object everywhere". It is
bimodal, and the split runs by age.

The good instance is the adjoint (#310 C3): `_reverse_octant_traversal`
(`loss_representation/__init__.py:750-787`) maps the octant label
`o ↦ −o` and then **looks up a different graph in the same cached
family**; the walk body is byte-for-byte untouched, and the transposed
data flow falls out of `face_in(−o) == face_out(o)` (verified exactly at
d = 1, 2, 3). Its sibling `_outflow_faces` (`:736-747`) is
`face_opposite ∘ _inflow_faces` — a map applied, not a twin transcribed,
with a docstring that says precisely that.

Meanwhile, four methods sitting inside ~100 lines of ONE file
(`sn/mesh/augmented_mesh.py:1275`, `:1398`, `:1426`, `:1499`) each
re-write `range(nx) if μ ≥ 0 else range(nx-1, -1, -1)`, two of them
alongside a hand-rolled `select_outer: bool` face swap; the domain-edge
index flip is written twice more with the ternary mirrored between the
two functions (`:1875` / `:1914`); and `_x_scan_faces`
(`sn/sweep/scan.py:302-363`) takes a literal **`x_reverse: bool`
parameter** — the boolean-flag-parameter anti-pattern the project's own
`coding-elegance` skill flags. The codebase already knows the elegant
form. It just has not been propagated backwards into the 1-D paths.

**Runner-up.** The CYL/SPH radial mesh's quotient IS named — in
orbit-space language, twice (`structured_geometry.py:44-61` and
`docs/theory/foundations/structured_geometry.rst:118-121`), including
the phrase "orbit-space rank". It is the only place in the whole spatial
layer where a quotient is named at all, and even there the *group* never
is. Its deck transformation is separately alive, exact, and named
"mirror" 500 lines away in a different package
(`loss_representation/__init__.py:3294`).

---

## Verification note (L-012 close-out)

HEAD moved `b0a003b4 → bfedc621` DURING this dispatch (the `Z2 → Mirror`
retirement + the quadrature-campaign E-series landed). Touched in-scope
files: `numerics/face_layout.py`, `numerics/spaces/angular_trace_space.py`,
`sn/loss_representation/__init__.py`, `sn/loss_representation/sweep_schedule.py`,
`sn/mesh/augmented_mesh.py`, `sn/mesh/method_space.py`,
`sn/operators/boundary.py`, `sn/acceleration/dsa.py`, `transport/method.py`.
**Every load-bearing `file:line` in this report was re-verified by
`sed -n Np` against `bfedc621` at close** — all 24 spot-checked citations
still land on the quoted construct. `orpheus/geometry/*.py` (minus
`boundary/`), `orpheus/transport/spatial/`, `orpheus/transport/mesh/`,
`orpheus/sn/sweep/` and `sweep_graph.py` were NOT touched in that range.
