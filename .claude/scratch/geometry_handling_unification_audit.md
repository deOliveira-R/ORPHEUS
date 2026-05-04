# Geometry-Handling Unification Audit

**Author**: explorer
**Date**: 2026-05-03
**Trigger**: Pre-R1/R2/R3 architectural reconnaissance for the
`trajectory_resolvent` hindsight refactor. User directive:
*"see if we can leverage something from /orpheus/geometry/. It's
totally fine if not."*
**Read scope**: `orpheus/geometry/` (3 files, ~470 LOC),
`orpheus/derivations/continuous/trajectory_resolvent/` (8 modules,
~5000 LOC), `orpheus/derivations/continuous/fn_method/` (slab,
sphere, core, multi_group), `tests/cross_method/` (protocol,
adapters, cases, test_eigenvalue), the in-flight refactor plan
`.claude/plans/trajectory_resolvent_hindsight_refactor.md`,
`.claude/agent-memory/cross-domain-attacker/trajectory_resolvent_foreign_frames.md`,
plus the existing usages of `Mesh1D`/`Zone` in
`orpheus/derivations/continuous/sood_registry/{la13511,builders}.py`.

This is a **research + assessment deliverable**. No code changes.

---

## Executive summary (TL;DR)

The three subsystems represent geometry on **three structurally
different levels of the stack** and the comparison is genuinely
asymmetric:

1. `orpheus/geometry/` is a **discretised mesh** abstraction
   (cell edges + per-cell material IDs + BC tags). It exists to
   feed sweep-based / collision-probability / diffusion solvers
   that produce per-cell scalar/angular flux on a fixed grid.

2. `trajectory_resolvent` represents geometry as **a continuous
   shape descriptor + an internal trajectory parameterisation**.
   The "mesh" the solver builds (Gauss-Legendre nodes on
   `(0, R)`, `(-1, 1)`, `(0, 2π)`) is a **quadrature** in
   phase-space, not a material mesh. The shape descriptor is
   today a few scalars (`R`, or `(R_in, R_out)`, or `L`) plus a
   pair of α reflectivities.

3. `fn_method` does not have any geometry data structure at all.
   It is parametrised on a single dimensionless number `c`
   (multiplication ratio), takes the geometry as a **string tag
   chosen by the file the function lives in** (`slab/`,
   `sphere/`, `cylinder/` directories), and returns a
   dimensionless critical dimension in mean free paths.

**The Mesh1D abstraction is load-bearing for production
solvers but a poor fit for either reference family.** The two
reference families don't *need* a discretised mesh — they need a
shape descriptor + cross sections. There is, however, a real
opportunity at a layer above geometry: the **CrossMethodCase**
protocol that already exists (`tests/cross_method/protocol.py`)
already plays the unification role I would otherwise nominate.

**Bottom line recommendation**: do NOT try to make `trajectory_resolvent`
or `fn_method` consume `Mesh1D` directly. Instead, recognise the
**already-existing** unification layer — the Sood registry's
`MeshTemplate` (currently a method-agnostic `(geometry_kind,
critical_dimension, n_groups, BCs)` recipe that already produces a
production `Mesh1D` *on demand*) — and **promote it into a small
shared "shape descriptor" type that all three families can agree on
at their public interfaces**. The call sites where geometry data
crosses a method boundary (cross-method tests, sood_registry case
construction) are the only places this matters; the per-method
internals stay as they are.

This is not a missed unification — it's a **half-complete one** that
already works for the production-side path (sood_registry →
`MeshTemplate.build()` → `Mesh1D` → `solve_cp` / `solve_sn`) and just
needs the reference-side path filled in. R1/R2/R3 plan should be
informed by this finding but **not blocked by it**: the geometry
unification is independent of the ChordOracle / power_iterate work
and can land in any order.

---

## Q1 — Per-subsystem geometry inventory

### 1.1 Production: `orpheus/geometry/`

**Files**:

- `coord.py`: `CoordSystem` enum {CARTESIAN, CYLINDRICAL, SPHERICAL};
  pure functions for volume / surface formulas.
- `mesh.py`: `BC` (kind tag + params dict), `Mesh1D` (frozen
  dataclass), `Mesh2D` (frozen dataclass).
- `factories.py`: `Zone`, `mesh1d_from_zones`, plus PWR-convenience
  factories.

**Primary types**:

| Type | Role | Public attributes |
|------|------|-------------------|
| `CoordSystem` | Coord-system selector | `CARTESIAN`, `CYLINDRICAL`, `SPHERICAL` |
| `BC` | Solver-agnostic boundary tag | `kind: str`, `params: dict[str, float]`; convenience instances `BC.vacuum`, `BC.reflective`, `BC.white` |
| `Mesh1D` | 1-D discretised mesh | `edges: ndarray (N+1,)`, `mat_ids: ndarray (N,)`, `coord: CoordSystem`, `precomputed_volumes: ndarray | None`, `bc_left: BC | None`, `bc_right: BC | None`; derived properties `widths`, `centers`, `volumes`, `surfaces`, `total_width`, `N` |
| `Mesh2D` | 2-D discretised mesh | `edges_x`, `edges_y`, `mat_map`, `coord`, four `bc_*` faces |
| `Zone` | Zone-based mesh recipe | `outer_edge`, `mat_id`, `n_cells` |

**Geometry contract** (what a function accepts when it talks
about geometry):

```python
def solve_sn(materials: dict[int, Mixture], mesh: Mesh1D | Mesh2D, ...): ...
def solve_cp(mesh: Mesh1D, materials: dict[int, Mixture], ...): ...
def solve_moc(mesh: Mesh1D, ...): ...
```

The contract is `mesh: Mesh1D` (or `Mesh2D`) PLUS `materials:
dict[int, Mixture]`. The mesh carries **edges + material IDs +
BC tags**; the materials dict carries the cross-section data
keyed by the integer IDs the mesh references.

**BC representation**: `BC` is a frozen dataclass with `kind: str`
and a `params: dict[str, float]` for things like albedo. Each
solver registers a `BC_REGISTRY` (e.g.,
`orpheus/sn/geometry.py:78`) that maps the `kind` string to a
solver-specific factory at construction time. Reflective and
vacuum are universal; white/Marshak/etc. are solver-extensible.

**Mesh / discretisation representation**: arbitrary cell edges
(ascending). Per-cell material ID (`mat_ids` shape `(N,)`).
Equal-volume subdivision is supported via the `Zone` factory
which precomputes volumes as a scalar per zone (avoids the
`cbrt(x)**3 != x` ULP-loss issue, see `factories.py` docstring).

### 1.2 `trajectory_resolvent`

**Files** (all under
`orpheus/derivations/continuous/trajectory_resolvent/`):

- `variant_alpha_core.py`: pure functions `compute_resolvent_T{,_rank2}`,
  `apply_variant_alpha_closure{,_rank2}` on already-computed
  `(F, B, τ_first, τ_period)` scalar/array tuples.
- `greens_function.py` (sphere — 1G + MG + multi-region + fixed-source).
- `greens_function_cylinder.py` (cylinder — 1G + MG).
- `greens_function_slab.py` (slab symmetric — thin facade over asym).
- `greens_function_slab_asymmetric.py` (rank-2 — 1G + MG).
- `greens_function_hollow_sphere.py` (rank-1 + rank-2 b-partition).
- `greens_function_annulus.py` (3D analog of hollow sphere).

**Primary "geometry" types**: NONE. There is no `Geometry`
dataclass anywhere in this package. The geometry is encoded in:

1. **The function name itself** —
   `solve_greens_function_sphere`, `_cylinder`, `_slab`,
   `_hollow_sphere`, `_annulus`, `_slab_asymmetric`. The
   geometric kind is a *static dispatch key* implemented as
   per-file functions.
2. **Positional / keyword arguments**: `R`, `(R_in, R_out)`,
   `L`, plus per-surface reflectivity `α` or `(α_left,
   α_right)` / `(α_in, α_out)`.
3. **Internal quadrature grids** built inside each `solve_*` —
   `r_nodes` from `np.polynomial.legendre.leggauss(n_r)` mapped
   to `[0, R]` (sphere) or `[R_in, R_out]` (hollow), `mu_nodes`
   on `[-1, 1]`, `phi_az_nodes` on `[0, 2π)` (cylinder/annulus).

**Geometry contract** (what each `solve_*` accepts):

| Solver | Geometry args | BC args | Quadrature args |
|--------|---------------|---------|-----------------|
| sphere 1G | `R: float` | `alpha: float = 1.0` | `n_r, n_mu, n_traj_quad` |
| sphere MG / MR | `R: float` (MG) or `radii: ndarray` (MR) | `alpha: float` | same |
| cylinder | `R: float` | `alpha: float` | `n_r, n_mu_axial, n_phi_az, n_traj_quad` |
| slab sym | `L: float` | `alpha: float` | `n_x, n_mu, n_traj_quad` |
| slab asym | `L: float` | `alpha_left, alpha_right` | same |
| hollow sphere | `R_in, R_out: float` | `alpha_in, alpha_out` | `n_r, n_mu, n_traj_quad` |
| annulus | `R_in, R_out: float` | `alpha_in, alpha_out` | `n_r, n_mu_axial, n_phi_az, n_traj_quad` |

**BC representation**: a single (or pair of) `float` named
`alpha`/`alpha_left`/`alpha_right`/`alpha_in`/`alpha_out`,
constrained to `[0, 1]`. There is **no `BC` instance and no
`kind` string**. The continuous family `α=0` is vacuum, `α=1`
is perfect specular, `α∈(0,1)` is partial-reflection albedo.
This is structurally rich (Variant α absorbs three production
BC kinds into one parametric family) but vocabulary-wise
disjoint from production's `BC("vacuum")`/`BC("reflective")`
tag system. White / Marshak diffuse BCs are NOT supported by
trajectory_resolvent (different closure structure — see
`trajectory_resolvent_hindsight_refactor.md` §A1).

**Mesh / discretisation representation**: there is **no
material mesh**. Single-region solvers use the homogeneous
medium implied by the scalar `(σ_t, σ_s, νσ_f)` arguments. The
multi-region sphere solver (`solve_greens_function_sphere_mr`)
does carry a region-edge concept — `radii: ndarray` (outer
radii, ascending, ending in `R`) plus per-region cross-section
arrays of shape `(n_regions, G)` / `(n_regions, G, G)`. There
are NO cell edges; the radii vector is just a region-boundary
list. The internal radial quadrature `r_nodes` is a single
Gauss-Legendre on `(0, R)` annotated with `region_at_node[i]`
saying which region each quadrature node falls in. Multi-region
trajectory machinery (`_trajectory_segments`, `_chord_segments`)
splits each chord at the `radii` values to compose
piecewise-σ_t attenuation along a single trajectory.

The internal `r_nodes` are a **quadrature in the integral
operator**, not a material mesh. They're not exposed as a "mesh"
type and they have no analogue of `mat_ids` / `volumes` /
`edges`.

### 1.3 `fn_method`

**Files** (under
`orpheus/derivations/continuous/fn_method/`):

- `core/{dispersion.py, fn_matrix.py, moments.py}`: shared F_N
  primitives (Case ν₀ root, F_N matrix assembly, moment
  recursions). Pure scalar / matrix algebra in `c` and the
  collocation grid `ξ`.
- `slab/{one_group.py, reflected.py, flux_reconstruction.py}`:
  bare-critical slab F_N (Siewert-Benoist + Grandjean-Siewert
  1979); reflected slab F_N (Neshat-Maiorino 1980).
- `sphere/{one_group.py, flux_reconstruction.py}`: bare-critical
  sphere F_N (Siewert-Thomas 1986).
- `cylinder/__init__.py`: stub (Mitsis WH not convergent for
  cylinder).
- `multi_group/k_inf.py`: pure rational algebra in cross
  sections; no spatial discretisation.
- `cases.py`: capability_rows() metadata.

**Primary "geometry" types**: NONE.

**Geometry contract** (what each `solve_fn_*` accepts):

| Solver | Geometry args | BC args | Internal params |
|--------|---------------|---------|-----------------|
| `solve_fn_slab_bare_critical(c, ...)` | none | implicit (vacuum both ends) | `n_modes, a_min, a_max, n_bracket, bisect_tol` |
| `solve_fn_sphere_bare_critical(c, ...)` | none | implicit (vacuum at R) | `n_modes, R_min, R_max, n_bracket` |
| `solve_fn_slab_reflected_critical(c_core, c_reflector, reflector_half_thickness, ...)` | `Δ: float` | implicit (vacuum at outer reflector) | `n_modes, max_outer_iter, tol` |

The geometric kind (slab vs sphere) is a **directory-level
dispatch key** — `fn_method.slab` vs `fn_method.sphere`. The
reflected slab adds **one geometric scalar** (`reflector_half_thickness`
in mfp). No BC kwargs at all; vacuum is implicit in the F_N
formulation. The output is a **single scalar critical
dimension** (`a_critical_mfp` for slab, `R_critical_mfp` for
sphere, `tau_critical_mfp` for reflected slab) plus the F_N
expansion coefficients.

**BC representation**: NONE — the F_N method is
formulated against vacuum only; the reflected slab buys
non-vacuum-on-the-outer-vacuum-still-on-the-edge by adding a
finite reflector region (still vacuum at the global outer face).
There is no `BC` object, no string tag, no albedo parameter.
Different boundary conditions would mean a different solver
(e.g. albedo BCs in F_N exist in literature but are not
implemented here).

**Mesh / discretisation representation**: NONE — F_N is a
collocation method on **moments of the angular flux at the
slab/sphere boundaries**. The "discretisation parameter" is
`n_modes` (the F_N order N). Internal "grid" is the collocation
points `ξ` on `[0, 1]` ∪ `{ν₀}` — a fixed N+1-vector of
moments, NOT a spatial mesh. There are no cells, no edges, no
material IDs.

### 1.4 Production-side existing bridge: `MeshTemplate` in
`sood_registry/la13511.py`

This is the most important finding for the audit. There is
ALREADY a method-agnostic geometry descriptor in the codebase
that consumes both production and reference paths:

```python
@dataclass(frozen=True)
class MeshTemplate:
    geometry: str            # "slab" / "sphere" / "cylinder" / "infinite" / "ISLC"
    critical_dimension_mfp: float | None
    critical_dimension_cm: float | None
    n_groups: int
    mat_id: int = 0
    bc_left: BC = BC.vacuum
    bc_right: BC = BC.vacuum

    @property
    def domain_extent_cm(self) -> float: ...

    def build(self, n_cells: int = 64) -> Mesh1D:
        # Builds a production Mesh1D via mesh1d_from_zones
```

This type already exists in `orpheus.derivations.continuous.sood_registry.la13511`
and:

- consumes `BC` from `orpheus.geometry` (production type),
- produces a `Mesh1D` from `orpheus.geometry` (production type),
- carries the geometric data needed to feed reference solvers
  (`critical_dimension_mfp` for fn_method, `domain_extent_cm`
  for trajectory_resolvent — both derivable from the same
  scalars).

The `cross_method/adapters.py` adapter classes today bypass this
by extracting cross sections directly from `case.registry_case.materials`
and reading the truth dimension off `case.truth_value`. That works
but **the geometry contract for each adapter is method-bespoke**:

- `FNSlabAdapter.solve(case)` reads `c` via `_extract_c(case)` →
  computes `(σ_t, σ_s, νσ_f)`, then ratios — passes only `c`.
- `TrajectoryResolventSlabAdapter.solve(case)` reads
  `(σ_t, σ_s, νσ_f)` plus the case's truth half-thickness in
  mfp, converts to `L_full_cm = 2 * a_truth_mfp / σ_t` —
  passes `(L, σ_t, σ_s, νσ_f, alpha=0)`.

The adapters do today exactly what a shared geometry descriptor
would do. They could be cleanly factored if both reference
families took a `(geometry_kind, critical_dimension, materials,
bc)` tuple, but they currently don't.

---

## Q2 — Concept-by-concept comparison

| Concept | Production (`orpheus/geometry/`) | trajectory_resolvent | fn_method | Convergeable? |
|---|---|---|---|---|
| **1-D radial mesh (sphere)** | `Mesh1D(coord=SPHERICAL)` w/ N+1 edges + N material IDs | NONE; `R: float` + internal `r_nodes` quadrature on `(0, R)` | NONE; `R_critical_mfp: float` is the OUTPUT, never input | **No** for reference solvers — they don't need a mesh, they need a shape descriptor. The `R: float` input IS the shape descriptor. |
| **1-D axial mesh (slab)** | `Mesh1D(coord=CARTESIAN)` over `[0, 2a]` | NONE; `L: float` + `x_nodes` quadrature on `(0, L)` | `a: float` is OUTPUT for bare; for reflected, `Δ: float` is INPUT | **No** for the same reason. |
| **Multi-region radial mesh** | `Mesh1D` w/ varying `mat_ids[i]` | `radii: ndarray` (region outer radii) + per-region XS arrays of shape `(n_regions, ...)`. NO edges/cells. | `(c_core, c_reflector, Δ)` triple — NO ndarray | **Partial**. trajectory_resolvent's `radii` could derive from `Mesh1D.edges[unique_breakpoints_in_mat_ids]`, but the per-region XS arrays would still need a separate `materials: dict[int, Mixture]` lookup. Not a free win. |
| **Hollow / annular topology** | NOT SUPPORTED by Mesh1D; would need `edges` to start at R_in not 0 — Mesh1D allows this (just ascending) but no factory currently produces it | `(R_in: float, R_out: float)` pair | NOT SUPPORTED | **Partial** — production Mesh1D would need a `mesh1d_from_zones(origin=R_in)` call, which is already supported. trajectory_resolvent could *consume* a Mesh1D over `[R_in, R_out]` and read `R_in`, `R_out` off `mesh.edges[0], mesh.edges[-1]`. But the solver doesn't iterate over cells; it just needs the two scalars. Pure overhead unless multi-region annulus is added. |
| **Boundary condition** | `BC(kind: str, params: dict)` w/ per-solver `BC_REGISTRY` resolving the kind | `alpha: float ∈ [0, 1]`; per-surface `(alpha_in, alpha_out)` etc. | NONE (vacuum implicit) | **Yes, with a thin adapter.** `BC("vacuum") ↔ alpha=0`, `BC("reflective") ↔ alpha=1`, `BC("partial", {"albedo": x}) ↔ alpha=x`. White/Marshak in production map to NOTHING in trajectory_resolvent (different closure structure). The mapping is one-way: production → trajectory_resolvent works; trajectory_resolvent's continuous α covers a richer space than production's tag system can express, but the cross-method tests only need the two corner cases (`α=0`, `α=1`). |
| **Material per region** | `mat_ids: ndarray (N,)` + separate `materials: dict[int, Mixture]` lookup | scalar `(σ_t, σ_s, νσ_f)` for 1G; arrays of shape `(G,)` / `(G, G)` for MG; arrays of shape `(n_regions, ...)` for MR | scalar `c` (1G) or scalar matrix (MG `k_inf`); F_N never sees a region-keyed mixture | **Partial.** trajectory_resolvent could derive its inputs from `(materials_dict, mat_ids, regions_in_order_of_increasing_mat_id)`. For MG/MR cases this would be a real ergonomics win in `cross_method/adapters.py`. For F_N it remains a c-scalar derivation. |
| **Coordinate system** | `CoordSystem` enum {CARTESIAN, CYLINDRICAL, SPHERICAL} | implicit in the function name (`_sphere`, `_cylinder`, `_slab`) | implicit in the directory (`slab/`, `sphere/`, `cylinder/`) | **Yes.** The `CoordSystem` enum is the natural method-agnostic dispatch key. `MeshTemplate` already has a `_GEOMETRY_TO_COORD` dict that maps the string tag to a `CoordSystem`. |
| **Volume / fission-rate Jacobian** | `Mesh1D.volumes` derived from `coord` + `edges` | inlined in each `solve_*`: sphere `4π · ∫ φ(r) r² · w_r`, cylinder `2π · ∫ φ(r) r · w_r`, slab `∫ φ(x) · w_x` | NONE — `c` is dimensionless, no spatial integral | **Partial — high leverage.** The geometry-specific `fission_rate(phi)` closures inside trajectory_resolvent's 12 power-iteration loops (see `trajectory_resolvent_hindsight_refactor.md` Phase 1) compute exactly the Jacobian that `Mesh1D.volumes + np.sum` already encapsulates for production solvers. **This is the one place where adopting `Mesh1D` would meaningfully reduce code.** Refactor candidate B5 (`mesh.volume_integral`) in the hindsight plan. |
| **Quadrature** | NOT IN GEOMETRY MODULE — lives in `sn/quadrature.py` etc., per-solver | inline `np.polynomial.legendre.leggauss(n_*)` in each `solve_*`. No abstraction. | F_N collocation grid `ξ` (Grandjean-Siewert 3+evenly-spaced) or Siewert-Thomas Chebyshev | **No.** Quadrature is intentionally not in `orpheus/geometry/`. The hindsight plan's `PhaseSpaceMesh` (B3) is the right home. |
| **Bounce-trajectory parameterisation** | NOT APPLICABLE | `_first_leg_2d_chord`, `_impact_parameter`, `_bounce_period_2d_chord`, `_trajectory_segments`, `_chord_segments` (per geometry) | NOT APPLICABLE — F_N never traces a trajectory; the chord physics is collapsed into the moment-recursion formulae `A_α(ξ)`, `B_α(ξ)` | **No — true asymmetry.** The chord oracle is a trajectory_resolvent-native concept (and a peierls_nystrom-native one too — see hindsight plan §P1). It has zero analog in production geometry or fn_method. The right home is `orpheus/derivations/common/chord_oracle.py` (per the hindsight plan), NOT `orpheus/geometry/`. |
| **Critical-dimension OUTPUT representation** | NOT APPLICABLE — production solvers solve k-eigenvalue at fixed mesh, output is `k_eff` | NOT APPLICABLE — same as production. trajectory_resolvent solves at fixed `R/L` and outputs `k_eff`. To find a critical dimension you wrap with a root-finder. | INTRINSIC OUTPUT — `solve_fn_*_bare_critical(c, ...)` returns `a_critical_mfp` / `R_critical_mfp` directly. The dispersion-root finder is internal. | **No, by design.** F_N is a critical-dimension-INPUT-→-c-OUTPUT (or c-INPUT-→-critical-dim-OUTPUT) solver; trajectory_resolvent is a fixed-config-→-k_eff solver. The `cross_method/protocol.py` `ScalarTag` union (`"k_eff"` vs `"a_critical_mfp"` vs `"R_critical_mfp"` vs `"tau_critical_mfp"`) handles this correctly today by tagging which scalar each adapter returns. |

---

## Q3 — Architectural opportunities for R1/R2/R3

### a) Adopt-as-is — where reference solvers should switch to production geometry types

**Single nominee**: `BC` for the cross-method adapter layer.

The mapping `BC.vacuum ↔ alpha=0` and `BC.reflective ↔ alpha=1`
is unambiguous. The `BC.params` dict already accommodates an
`{"albedo": x}` extension for partial reflection; defining
`BC("partial", {"albedo": x}) ↔ alpha=x` is a one-line BC factory
in trajectory_resolvent's analog of `SNMesh.BC_REGISTRY`. This
would let cross-method tests express "vacuum" once via the
production tag system and have both solvers consume it.

**Cost**: trivial. Define a thin
`alpha_from_bc(bc: BC) -> float` helper in the cross-method
adapters layer (NOT in `orpheus/geometry/`). The adapters today
hardcode `alpha = 0.0` for bare-critical cases; replacing this
with a `BC`-driven mapping makes new BC kinds (white, etc.)
automatically signal "trajectory_resolvent unsupported, raise" at
the protocol level rather than failing inside the solver.

**Risk**: zero. Pure ergonomics, no math change.

**Does NOT apply to** `Mesh1D` itself. trajectory_resolvent
genuinely doesn't need a mesh data structure — its quadrature is
phase-space, not material-space.

### b) Promote-and-unify — concepts that should be promoted to a shared layer that both reference AND production consume

**Top candidate**: `MeshTemplate` — the method-agnostic geometry
recipe.

Currently lives in `orpheus.derivations.continuous.sood_registry.la13511`
and is used to build `Mesh1D` instances for production-solver
verification tests. Its fields
`(geometry_kind, critical_dimension_{mfp,cm}, n_groups, mat_id,
bc_left, bc_right)` are method-agnostic — they describe the
domain shape and BC, not the solver. If promoted to
`orpheus.derivations.common.geometry_template` (or similar), both:

- production solvers consume `template.build(n_cells)` → `Mesh1D`
- trajectory_resolvent consumes `template.domain_extent_cm` →
  `R` or `L` directly, without needing the discretised mesh
- fn_method consumes `template.critical_dimension_mfp` and the
  geometry kind dispatch
- cross_method adapters consume the template once, dispatch all
  three families from the same case description

**Cost**: small (~half day). The class moves verbatim; its three
existing call sites get import updates. The `cross_method/cases.py`
populated case sets get refactored to use the shared template
instead of inline `geometry: str` strings.

**Risk**: low. The class is already well-tested; its consumers
(`builders.py`) are minimal. The main risk is that adoption by
trajectory_resolvent's adapter creates a circular import:
`orpheus.derivations.continuous.fn_method` and
`orpheus.derivations.continuous.trajectory_resolvent` would both
import from `orpheus.derivations.common.geometry_template`, which
itself imports `BC`/`Mesh1D` from `orpheus.geometry`. The
dependency direction is acyclic so this is fine.

**Naming caveat**: don't name it `Geometry` or `GeometrySpec` —
the project already has `CoreGeometry`
(`orpheus/diffusion/solver.py`), `CurvilinearGeometry`
(peierls_nystrom), `FlatSourceCPGeometry`. The token "Geometry"
is overloaded. **`MeshTemplate` is exactly right** — it's a
recipe, not the discretised object. (The cross-domain-attacker
memo `trajectory_resolvent_foreign_frames.md` ranks "Bundle
refactor (Frame 4)" as a wait-for-instance-N+1 candidate; the
`MeshTemplate` is the *light* version of bundle's BaseAtlas — a
shape descriptor without the differential-geometry machinery.)

**Second candidate (lower priority)**: a `materials_at_region(template,
materials_dict, region_idx) -> Mixture` adapter. This factors out
the `case.registry_case.materials[0]` extraction that
`adapters._extract_1g_xs` does today; applies cleanly to MG and
MR cases. Pure ergonomics; lands as a corollary of the
`MeshTemplate` promotion.

### c) Keep separate — true asymmetries that justify separate types

These are the asymmetries the user must accept; trying to unify
them creates either premature abstraction or false unification.

1. **Bounce-trajectory parameterisation (chord oracle)** —
   trajectory_resolvent only. F_N's chord physics is collapsed
   into the analytic moment recursions; production solvers'
   chord physics is collapsed into per-cell streaming
   coefficients. The `ChordOracle` Protocol (hindsight plan §B1)
   belongs in
   `orpheus/derivations/common/chord_oracle.py` because peierls_nystrom
   needs it too — but it does NOT belong in
   `orpheus/geometry/` and it does NOT speak to fn_method. Cross-domain-attacker
   memo ranks this as the "fiber bundle" pattern at the highest
   leverage; this audit confirms it.

2. **Phase-space quadrature (`PhaseSpaceMesh`)** — trajectory_resolvent
   only (per the hindsight plan §B3). The
   `(n_r, n_mu, n_phi_az, n_traj_quad)` parameter family is
   genuinely native to bouncing-trajectory operators; `Mesh1D`
   has no analog. Should remain inside trajectory_resolvent.

3. **F_N order `n_modes` and collocation grid** — fn_method only.
   Has no analog in trajectory_resolvent or production geometry.

4. **F_N's c-only parametrisation** — fn_method-internal. The
   reduction `(σ_t, σ_s, νσ_f) → c = (σ_s + νσ_f)/σ_t` is the
   F_N method's normal form for the dispersion relation.
   Production solvers don't reduce; trajectory_resolvent doesn't
   reduce. The reduction lives in `_extract_c` in the adapters
   layer, exactly where it belongs.

5. **Continuous `α ∈ [0, 1]` parametrisation** —
   trajectory_resolvent only. Richer than production's tag
   system; cannot be folded into `BC` without losing the
   parametric continuity that makes Variant α work. The ONLY
   bridge is the corner-case mapping in (a) above.

6. **Critical-dimension-as-output vs k_eff-at-fixed-config** —
   structural difference between fn_method (root-find solver)
   and trajectory_resolvent / production (forward solver). The
   `CrossMethodCase.truth_tag` Literal handles this at the test
   protocol level; it should NOT be unified at the solver level.

7. **Reflected-slab geometry** — fn_method native (Neshat-Maiorino
   1980 implementation lives in `slab/reflected.py`). The
   reflected-slab F_N adapter notes (in
   `tests/cross_method/adapters.py:120-167`) that there is
   currently no trajectory_resolvent counterpart. The
   trajectory_resolvent slab-asymmetric module has the algebraic
   machinery (rank-2 closure with `α_left ≠ α_right`) but does
   NOT today host a finite reflector — that would need the slab
   asymmetric module to gain a multi-region capability with a
   per-region σ_t, mirroring sphere's MR extension. **Until
   that work lands, the reflected slab is one-sided coverage,
   period. The geometry abstraction won't fix this.**

### d) Adapter pattern — where direct unification is too invasive but a thin adapter could let cross-method tests use one geometry spec for everything

**This is the recommended path.**

The cross-method test infrastructure (`tests/cross_method/`)
ALREADY has the right shape:

- `CrossMethodCase` carries `geometry: str`, `truth_tag`,
  `truth_value`, `tolerances`, `notes`.
- `SolverAdapter` Protocol with per-adapter `solve(case) ->
  ScalarResult`.
- `ADAPTERS_BY_NAME` registry.
- Per-adapter helpers `_extract_1g_xs`, `_extract_c`,
  `_slab_a_truth_mfp`, etc., that bridge the gap.

The adapter layer is **already** the unification point. What's
missing is that the case construction is currently one of:

- `case.registry_case = LA13511_CASES["Ua-1-0-SL"]` (registry
  case carries a `MeshTemplate` already), OR
- `case.registry_case = None` plus parameters parsed from
  `case.notes` (k=v string).

If the cross-method case held a `mesh_template: MeshTemplate |
None` field directly (rather than relying on
`case.registry_case.mesh_template`), all four adapter helpers
collapse to:

- `template.domain_extent_cm` for the shape scalar.
- `template.bc_left`, `template.bc_right` for the BC.
- `mixture_to_fn_arrays(case.materials[mat_id])` for the XS.

The `_parse_notes_kv` ad-hoc string parsing in
`adapters.py:441-454` (used today for reflected slab and
closed-sphere cases without a registry entry) goes away —
those cases would just construct an inline `MeshTemplate`
with their non-standard parameters.

**Cost**: small to medium (~1 day). Touches
`tests/cross_method/cases.py` (4 case sets, ~12 cases),
`tests/cross_method/adapters.py` (helpers consolidate), and
adds one helper in
`orpheus/derivations/continuous/sood_registry/extractors.py`.

**Risk**: low. The adapter helpers are well-tested by the
existing cross-method test suite (`test_eigenvalue.py`, ~20
tests). Bit-equality of converted XS arrays is preserved (no
arithmetic changes).

**Verifying-tests-become-trivial benefit (the user's stated goal)**:
right now to add a new (geometry × method) cross-check test, the
case author writes a `CrossMethodCase` with the right
`tolerances` map AND the case authors must remember which
adapter helper handles their notes-vs-registry split. After this
change, the case author writes a `MeshTemplate` (one line) and
opts in adapters via `tolerances` (one line each). For
`(homogeneous-sphere, fn-method × trajectory-resolvent ×
peierls-nystrom-when-it-lands)` cross-checks this is a 3-line
case. The user's "implementation of tests between those for
formal verification becomes trivial" is realised at the case-
construction level — exactly the layer where it pays.

---

## Implications for R1/R2/R3 refactor

The hindsight plan currently is:

- **R1 (B2)**: `power_iterate_variant_alpha` driver — eliminate
  12 power-iteration loops. LOW risk, HIGH leverage.
- **R2 (B3)**: `GreensResult` + `PhaseSpaceMesh` — collapse 12
  result dataclasses. LOW risk, MEDIUM leverage.
- **R3 (B1)**: `ChordOracle` Protocol — extract per-geometry
  chord/trajectory primitives. MEDIUM risk, HIGHEST leverage.

This audit suggests **two adjustments**:

### Adjustment 1 — `PhaseSpaceMesh` (R2 / B3) should NOT consume `Mesh1D`

The hindsight plan tentatively considers `PhaseSpaceMesh` as a
tagged-union dataclass with sphere/cylinder/slab/etc.
specialisations. This audit confirms it should remain
trajectory_resolvent-internal — the `Mesh1D` abstraction
genuinely doesn't fit (it's a discretised cell-edge mesh, not a
phase-space quadrature). The two share zero load-bearing
attributes.

What `PhaseSpaceMesh` **should** absorb (per B5 in the hindsight
plan): the `volume_integral(phi)` method that today is inlined
as `fission_rate` closure in 12 power-iteration loops. This is
the `Mesh1D.volumes`-equivalent for trajectory_resolvent, but
reimplemented on Gauss-Legendre weights rather than cell edges.
It should NOT delegate to `Mesh1D`.

### Adjustment 2 — INSERT a new candidate "R0.5": promote `MeshTemplate` to `common/`

This is the only finding in this audit that suggests a new R-phase
slot. The promotion is independent of R1/R2/R3 — it's a code-move
+ adapter-rewire, NOT a Variant α refactor.

Slot it BEFORE the cross-method regression net is extended (the
cross-method net just landed at commit `b215896`; further extension
is upcoming). Doing the promotion first means:

- New cross-method cases use `MeshTemplate` from day one.
- The existing 12 cross-method cases get a one-time refactor
  (mechanical, ~1 hour) to consume the new shared type.
- When the third method (peierls_nystrom or singular_eigenfunction
  via fn_method-cylinder pillar) adds its adapter, it consumes
  `MeshTemplate` directly with no special-case parsing.

**Cost vs payoff**: ~half a day of work; saves ~1 day of glue
when each new method-adapter pair joins the cross-method net.

### What does NOT change

R1, R2, R3 stand. None of the conclusions in this audit alter
the priority ordering or risk classification of the hindsight
plan's main candidates. The geometry audit lives in a layer
ABOVE Variant α (cross-method case construction) and is
orthogonal to the chord-oracle / power-iterate work.

The user's stated discipline ("build cylinder first, test it,
then unify and see if the tests still hold") still applies: the
`MeshTemplate` promotion is itself a small unification, and per
`feedback_unify_after_two_instances.md` it is justified because
the `MeshTemplate` already has TWO instances (the production-side
sood_registry usage and the to-be-built reference-side adapter
usage). The fact that `MeshTemplate` exists today as
production-side-only is the precondition that makes promotion
safe rather than premature.

---

## Recommendation

**Concrete next step**: BEFORE starting R1, do a small scout
commit:

1. Move `MeshTemplate` from
   `orpheus.derivations.continuous.sood_registry.la13511` to
   `orpheus.derivations.common.geometry_template` (new module).
2. Add a `from_critical_dimension_cm()` classmethod for
   constructing without a Sood registry case (so reflected-slab
   and closed-sphere inline cases can use it).
3. Refactor `tests/cross_method/cases.py` to carry `mesh_template:
   MeshTemplate` directly (or leave it on `case.registry_case` for
   registry-backed cases — both paths supported).
4. Refactor `tests/cross_method/adapters.py` helpers to consume
   `MeshTemplate` instead of `_parse_notes_kv` / `_extract_1g_xs`.
5. Run full `tests/cross_method/` suite + `tests/derivations/`
   suite — bit-equality at every test.

This is a ~half-day commit. It clears the path for the
upcoming peierls_nystrom adapter and the singular_eigenfunction
cylinder adapter without touching trajectory_resolvent's
internals. R1/R2/R3 then proceeds unchanged.

**Do NOT**:

- Adopt `Mesh1D` inside trajectory_resolvent's `solve_*`
  functions. The "mesh" inside those solvers is phase-space
  quadrature, not material-space, and forcing it through a cell-
  edge abstraction is a textbook premature unification.
- Attempt to express `α ∈ [0, 1]` partial-reflection in
  production's `BC` system right now. The continuous Variant α
  parametrisation is structurally richer than the tag-system; a
  bridge exists for the corner cases (the cross-method adapter
  helper above) but a full unification would force production to
  carry albedo-on-vacuum semantics it doesn't need.
- Block R1/R2/R3 on this audit. The geometry handling
  unification is a pre-requisite for the NEXT cross-method
  adapter to land elegantly, NOT for the chord-oracle /
  power-iterate refactor inside trajectory_resolvent.

**Honest divergence note** (per
`feedback_unify_after_two_instances.md` and the user's "It's
totally fine if not"): the bouncing-trajectory chord oracle, the
phase-space quadrature, and the F_N collocation grid are
genuinely native to their own methods. Pretending they share a
geometry abstraction with `Mesh1D` would be false unification.
The `MeshTemplate` promotion is the ONE place where unification
pays off — at the case-description boundary where method-agnostic
reasoning is real.

---

## Key files (absolute paths)

**Production geometry (read first if implementing)**:
- `/workspaces/ORPHEUS/orpheus/geometry/__init__.py`
- `/workspaces/ORPHEUS/orpheus/geometry/mesh.py` (`Mesh1D`,
  `Mesh2D`, `BC`)
- `/workspaces/ORPHEUS/orpheus/geometry/factories.py` (`Zone`,
  `mesh1d_from_zones`)
- `/workspaces/ORPHEUS/orpheus/geometry/coord.py` (`CoordSystem`)

**Existing production-side bridge (the half-built unification)**:
- `/workspaces/ORPHEUS/orpheus/derivations/continuous/sood_registry/la13511.py:76-218`
  (`MeshTemplate` class — promotion candidate)
- `/workspaces/ORPHEUS/orpheus/derivations/continuous/sood_registry/builders.py`
  (consumers: `build_materials`, `build_mesh`)
- `/workspaces/ORPHEUS/orpheus/derivations/continuous/sood_registry/extractors.py`
  (`mixture_to_fn_arrays` — XS bridge already used by adapters)

**trajectory_resolvent geometry handling**:
- `/workspaces/ORPHEUS/orpheus/derivations/continuous/trajectory_resolvent/variant_alpha_core.py`
  (4 closure primitives — pure functions, no geometry)
- `/workspaces/ORPHEUS/orpheus/derivations/continuous/trajectory_resolvent/greens_function.py:307-461`
  (`solve_greens_function_sphere` — canonical signature: `R, σ_t,
  σ_s, νσ_f, *, alpha=1.0, n_r, n_mu, n_traj_quad, max_iter, tol,
  initial_psi`)
- `/workspaces/ORPHEUS/orpheus/derivations/continuous/trajectory_resolvent/greens_function.py:913-1322`
  (multi-region sphere — `radii: ndarray` is the closest analog
  to a mesh)
- `/workspaces/ORPHEUS/orpheus/derivations/continuous/trajectory_resolvent/greens_function_cylinder.py:389-547`
  (`solve_greens_function_cylinder` — adds
  `n_mu_axial, n_phi_az`)
- `/workspaces/ORPHEUS/orpheus/derivations/continuous/trajectory_resolvent/greens_function_slab_asymmetric.py:327-450`
  (`solve_greens_function_slab_asymmetric` — adds
  `alpha_left, alpha_right`)
- `/workspaces/ORPHEUS/orpheus/derivations/continuous/trajectory_resolvent/greens_function_hollow_sphere.py:428-575`
  (adds `R_in, R_out`, `alpha_in, alpha_out`)
- `/workspaces/ORPHEUS/orpheus/derivations/continuous/trajectory_resolvent/greens_function_annulus.py:468-634`
  (3D analog of hollow_sphere)

**fn_method geometry handling**:
- `/workspaces/ORPHEUS/orpheus/derivations/continuous/fn_method/__init__.py`
- `/workspaces/ORPHEUS/orpheus/derivations/continuous/fn_method/cases.py`
  (capability matrix)
- `/workspaces/ORPHEUS/orpheus/derivations/continuous/fn_method/slab/one_group.py:121-260`
  (`solve_fn_slab_bare_critical(c, ...)` — geometry IS the file)
- `/workspaces/ORPHEUS/orpheus/derivations/continuous/fn_method/sphere/one_group.py:174-383`
  (`solve_fn_sphere_bare_critical(c, ...)`)
- `/workspaces/ORPHEUS/orpheus/derivations/continuous/fn_method/slab/reflected.py`
  (reflected-slab — adds one geometry scalar)
- `/workspaces/ORPHEUS/orpheus/derivations/continuous/fn_method/core/fn_matrix.py`
  (geometry sign ±1 dispatch — slab=+1, sphere=-1)

**Cross-method test infrastructure (the unification-payoff layer)**:
- `/workspaces/ORPHEUS/tests/cross_method/protocol.py`
  (`CrossMethodCase`, `SolverAdapter`, `ScalarResult`, `ScalarTag`)
- `/workspaces/ORPHEUS/tests/cross_method/adapters.py`
  (existing 6 adapters; helpers `_extract_c`, `_extract_1g_xs`,
  `_parse_notes_kv` — the eventual targets for `MeshTemplate`-based
  refactor)
- `/workspaces/ORPHEUS/tests/cross_method/cases.py`
  (5 case sets, 13 cases total — refactor candidate)
- `/workspaces/ORPHEUS/tests/cross_method/test_eigenvalue.py`
  (truth gates + cross-method agreement gates)

**Refactor plan this audit informs**:
- `/workspaces/ORPHEUS/.claude/plans/trajectory_resolvent_hindsight_refactor.md`
  (R1/R2/R3 main plan)
- `/workspaces/ORPHEUS/.claude/agent-memory/cross-domain-attacker/trajectory_resolvent_foreign_frames.md`
  (Frame 4 / Bundle refactor — wait-for-instance-N+1; this audit
  agrees, but the `MeshTemplate` promotion is the LIGHT version
  that does pay today)
- `/workspaces/ORPHEUS/.claude/scratch/cross_method_test_protocol_assessment.md`
  (Phase-1 assessment that built the protocol — this audit's
  recommendations are direct extensions)

---

## Q5 — Input-level unification (follow-up audit, 2026-05-03)

**Trigger**: User directive after the dataclass-level audit landed:
*"the input that geometries get might still be possibly the same.
So they might not be unified at the dataclass level, but at the
level of input to dataclass."*

The Q1–Q4 conclusions stand: `Mesh1D` (production) and the
trajectory_resolvent / fn_method input shapes cannot be unified at
the dataclass level. The new question is whether the **input** that
constructs each can come from a single canonical case-spec. Below
I trace each pipeline upward, identify the convergence point, and
check whether `MeshTemplate` already covers it.

### 5.1 Per-pipeline input-shape traces

#### Production `Mesh1D` — what calls its constructor in real use

Nexus call graph (`callers(py:class:orpheus.geometry.mesh.Mesh1D)`,
81 results) splits cleanly into three classes:

1. **`mesh1d_from_zones`** (factory) — constructed from a
   `list[Zone]` + `coord: CoordSystem`. The zone is itself
   `(outer_edge, mat_id, n_cells)`. Upstream callers of
   `mesh1d_from_zones` (33 distinct, mostly tests + the
   `factories.homogeneous_1d/pwr_*` convenience builders).
2. **`MeshTemplate.build`** — wraps `mesh1d_from_zones` with a
   single inline `Zone(outer_edge=domain_extent_cm, mat_id, n_cells)`.
3. **MMS test cases** — `SNSlabMMSCase.build_mesh(n_cells)`,
   `SNSphericalMMSCase.build_mesh(n_cells)`, etc., each constructs
   a `Mesh1D` from method-specific inline arrays (manufactured-
   solution domains have known L / R / nx by construction).

The example demo (`examples.discrete_ordinates.demo_discrete_ordinates.main`)
constructs `Mesh1D` directly with literal edge arrays.

**Source of inputs**:
- For the **Sood / sood_registry path** (the only path that
  cross-method tests touch): the case-spec layer is
  `La13511Case.mesh_template: MeshTemplate`. `MeshTemplate` carries
  `(geometry_kind, critical_dimension_mfp, critical_dimension_cm,
  n_groups, mat_id, bc_left, bc_right)`. Production solvers consume
  `case.mesh_template.build(n_cells=64)` → `Mesh1D`. **The case-spec
  layer is `MeshTemplate`. Period.** This is what `build_materials`
  / `build_mesh` / `build_cp_params` in
  `sood_registry/builders.py` formalise.
- For the **MMS path**: each case is its own dataclass
  (`SNSlabMMSCase`, `SNSphericalMMSCase`, ...) which carries the
  manufactured-solution domain inline. There is NO shared case-spec
  layer here; MMS cases are not consumed by reference solvers.
  Cross-method coverage of MMS is not on the table — MMS's reference
  is the manufactured solution itself, not another solver. So MMS
  is OUT of the unification question.
- For the **PWR factories** path: `pwr_pin_equivalent()`,
  `pwr_slab_half_cell()`, etc., are zone-based recipes. They are
  used in tests + the diffusion-solver default. Not a candidate for
  reference-solver consumption (they're heterogeneous; reference
  solvers handle homogeneous + a few multi-region cases).

So for the question that matters — **does the cross-method case
spec map to a single canonical input shape across all three
families?** — the production-side answer is unambiguously
`MeshTemplate`.

#### trajectory_resolvent — what calls its `solve_*` in real use

Nexus `callers(py:function:...trajectory_resolvent.greens_function.solve_greens_function_sphere)`,
20 results, splits as:

1. **`tests.cross_method.adapters.TrajectoryResolventSphereAdapter.solve`**
   and `TrajectoryResolventSphereClosedAdapter.solve` —
   `tests/cross_method/adapters.py:273-307`, 334-365.
2. **`tests.derivations.test_fn_la13511_sphere_xverif.test_fn_sphere_vs_variant_alpha_sphere_at_sood_ua_1_0_sp`**
   — pre-protocol cross-verification test.
3. **`tests.derivations.test_peierls_greens_function_*`** (multiple)
   — internal V&V tests calling `solve_greens_function_sphere`
   directly with literal cross sections.

The 1G / MG slab and cylinder solvers have similar caller
distributions (cross-method adapters + per-method V&V tests).

**Source of inputs**:
- The **cross-method adapter path**: input comes from
  `CrossMethodCase` via the helpers `_extract_1g_xs(case)` and
  `_sphere_R_truth_mfp(case) / sigma_t`. Concretely the adapter
  reads:
    - `case.registry_case.materials[0]` → 1G XS via
      `mixture_to_fn_arrays`.
    - `case.truth_value` (with `case.truth_tag == "R_critical_mfp"`)
      → R_truth_mfp, then `R_cm = R_truth_mfp / sigma_t`.
    - `self.alpha` (adapter constant 0 or 1) — NOT from case.
    - `self.n_r, self.n_mu, self.n_traj_quad` — adapter constants,
      NOT from case.
- The **per-method V&V test path**: each test inlines its own
  literal `(R, sigma_t, sigma_s, nu_sigma_f, alpha, ...)` —
  no shared case-spec involved.

**Critical observation**: `case.registry_case.mesh_template` is
ALREADY there. The adapter just doesn't use it. It reads
`case.truth_value` (which IS `mesh_template.critical_dimension_mfp`
for these cases — they're the same number) and re-derives
`R_cm` via `R_truth_mfp / sigma_t`. Meanwhile
`case.registry_case.mesh_template.critical_dimension_cm`
already holds the right answer in cm.

This is exactly the symptom Q4 anticipated. The geometry input
that trajectory_resolvent's adapter actually needs is:
1. `R: float` in cm — present in
   `case.registry_case.mesh_template.critical_dimension_cm`.
2. `(sigma_t, sigma_s, nu_sigma_f)` — present in
   `case.registry_case.materials[mat_id]` via
   `mixture_to_fn_arrays` (with `mat_id =
   case.registry_case.mesh_template.mat_id`).
3. `alpha: float` — **NOT in MeshTemplate**, derived from
   `bc_left/bc_right` via `BC.vacuum → 0`, `BC.reflective → 1`
   (and a future `BC.partial(albedo) → albedo` extension).

**Closed-sphere k_inf adapter** (`TrajectoryResolventSphereClosedAdapter`):
- `case.registry_case = None` for this case.
- Pulls everything from `case.notes` via `_parse_notes_kv`:
  `sigma_t, sigma_s, nu_sigma_f, R_cm`.
- `alpha = 1` from adapter constant (corresponds to
  `MeshTemplate(geometry="sphere", bc_left=BC.reflective,
  bc_right=BC.reflective)`).

#### fn_method — what calls its `solve_fn_*` in real use

Nexus `callers(py:function:...fn_method.slab.one_group.solve_fn_slab_bare_critical)`,
25 results — and `solve_fn_sphere_bare_critical`, 12 results.
Both split as:

1. **`tests.cross_method.adapters.{FNSlabAdapter, FNSphereAdapter,
   FNReflectedSlabAdapter}.solve`**.
2. **Per-method V&V tests** (`tests.derivations.test_fn_la13511_*`,
   `test_atkinson_product_nystrom`, `test_carlvik_galerkin_xverif_fn`).

**Source of inputs**:
- The **cross-method adapter path**: the adapter reads only `c =
  _extract_c(case) = (sigma_s + nu_sigma_f) / sigma_t` from
  `case.registry_case.materials[0]`. The geometry (slab vs sphere)
  is encoded in WHICH adapter is chosen — `FNSlabAdapter` calls
  the slab solver, `FNSphereAdapter` the sphere solver. There is
  no `geometry` field passed to the solver itself. The geometry
  comes in via the directory dispatch (`fn_method.slab` vs
  `fn_method.sphere`) and `geometry_sign = ±1` inside the matrix
  assembler.
- **Reflected slab adapter** (`FNReflectedSlabAdapter`):
  `case.registry_case = None`. Pulls
  `(c_core, c_reflector, reflector_half_thickness_mfp)` from
  `case.notes` via `_parse_notes_kv`. Three case-specific scalars
  — none of them in `MeshTemplate` today (more on this below).

### 5.2 Shared prefix in the input pipeline

Walking up each pipeline to the level at which all three families
read from the **same conceptual structure**, the convergence point is:

```
CrossMethodCase
  ├── case.registry_case: La13511Case | None
  │     ├── materials: dict[int, Mixture]   ← shared
  │     └── mesh_template: MeshTemplate     ← shared
  │           ├── geometry: str             ← shared
  │           ├── critical_dimension_cm     ← shared
  │           ├── critical_dimension_mfp    ← shared
  │           ├── mat_id: int               ← shared
  │           ├── bc_left: BC               ← shared
  │           └── bc_right: BC              ← shared
  ├── case.geometry: str   (duplicates mesh_template.geometry)
  ├── case.truth_tag       (method-output kind, not input)
  └── case.tolerances      (computational config, per-adapter)
```

For **Sood-backed cases**, every input each adapter needs is
already in `case.registry_case.{materials, mesh_template}`. The
adapters today do these conversions inline and inconsistently:

| Adapter need        | Source today                                          | Could come from                                              |
|---------------------|-------------------------------------------------------|--------------------------------------------------------------|
| `c` (FN)            | `_extract_c` ← materials[0]                           | (same) materials[mat_id]                                     |
| `(σ_t, σ_s, νσ_f)`  | `_extract_1g_xs` ← materials[0]                       | (same) materials[mat_id]                                     |
| `R_cm` / `L_cm`     | `case.truth_value / sigma_t` (re-derives from truth!) | `mesh_template.critical_dimension_cm` directly               |
| Geometry kind       | adapter constant                                      | `mesh_template.geometry` (≡ `case.geometry`)                 |
| `alpha`             | adapter constant                                      | from `mesh_template.bc_left/bc_right` via a one-line mapper  |

For **inline (non-registry) cases** — reflected slab and
closed-sphere k_inf — `case.registry_case = None` and the adapter
falls back to `_parse_notes_kv` to pull values from `case.notes`.
This is the EXACT symptom of input-level non-unification: the
moment a case can't fit in `MeshTemplate`, the adapters open a
side-channel through free-form text. Q4's anticipated answer is
confirmed by the code.

### 5.3 Divergence point

Each pipeline diverges from the shared prefix at a different layer
that is **genuinely method-specific**:

- **trajectory_resolvent** diverges at the *quadrature* layer:
  `(n_r, n_mu, n_traj_quad, n_phi_az, n_mu_axial, max_iter, tol,
  initial_psi)`. None of these are in `MeshTemplate` and none
  belong there — they are computational parameters, not physical
  ones. The adapter encodes them as adapter-class constants.
- **fn_method bare** diverges at *F_N order + bracket*:
  `(n_modes, a_min, a_max, n_bracket, bisect_tol, max_bisect)`.
  Same nature: computational, not physical.
- **fn_method reflected** diverges at the *reflector geometry*:
  `(c_core, c_reflector, reflector_half_thickness_mfp)`. The first
  two are physical (per-region `c` ratios); the third is genuinely
  geometric (reflector thickness in mfp). The reflector geometry
  here is structurally a multi-region `MeshTemplate` (core + two
  reflectors) but with a c-only parametrisation — see 5.4.
- **trajectory_resolvent multi-region sphere**
  (`solve_greens_function_sphere_mr`) diverges at the
  *region-boundary list*: `radii: ndarray` plus per-region
  `(σ_t, σ_s, νσ_f)` arrays. Structurally this is a multi-region
  `MeshTemplate` + `materials: dict[int, Mixture]` tuple — exactly
  the production shape with `mesh_template.geometry = "sphere"`,
  edges = the radii vector, and `mat_ids[i]` = which region each
  inter-radius shell belongs to. **This is the only divergence
  point where the input shape genuinely IS the production
  shape, just with a different consumer.**

The clean separation: **physical input** (geometry kind,
dimensions, materials per region, BC) is unifiable; **computational
input** (quadrature orders, F_N order, bracket, tolerances) is not
and shouldn't be.

This matches the cross_method protocol's existing partition between
`CrossMethodCase.tolerances` (which the case author opts adapters
into, with adapter-side numerical-parameter selection) and the case
itself (which describes the physics).

### 5.4 Q3 — what does the existing `MeshTemplate` already cover?

Verdict: **`MeshTemplate` covers ≈ 80 % of the unifiable input
layer for single-region cases. The reflected slab and the
multi-region sphere don't fit, but the gap is structural (multi-
region) not method-specific.**

Coverage today (from the current `MeshTemplate` definition at
`orpheus/derivations/continuous/sood_registry/la13511.py:76-214`):

| Field                       | Covered? | Comment                                                              |
|-----------------------------|----------|----------------------------------------------------------------------|
| `geometry: str`             | YES      | ``"slab" / "sphere" / "cylinder" / "infinite" / "ISLC"``             |
| `critical_dimension_mfp`    | YES      | F_N consumes this                                                    |
| `critical_dimension_cm`     | YES      | trajectory_resolvent + production consume this                       |
| `n_groups: int`             | YES      | adapter routing key                                                  |
| `mat_id: int`               | YES      | indexes into materials dict                                          |
| `bc_left: BC`, `bc_right`   | YES      | maps to α via BC↔α adapter (single line)                             |
| `domain_extent_cm` property | YES      | slab→2a, sphere/cylinder→R                                           |
| `build(n_cells)` method     | YES      | produces production `Mesh1D`                                         |

NOT covered by the current `MeshTemplate`:

| Concept                              | What it would need              | Belongs in MeshTemplate?  |
|--------------------------------------|---------------------------------|---------------------------|
| Multi-region radial topology         | `regions: list[(r_outer, mat_id)]`  | YES — it's geometry       |
| Hollow / annular topology            | `inner_edge: float \| None`     | YES — it's geometry       |
| Reflector kind for fn_method (c_refl) | per-region `c` (or per-region Mixture) | NO — already in `materials[mat_id]` |
| F_N order `n_modes`                  | computational                   | NO — adapter-side         |
| trajectory n_quad triple             | computational                   | NO — adapter-side         |
| `alpha` continuous albedo            | physical BC parameter           | EXTEND BC.params dict     |

The two **physical** gaps (multi-region, hollow) are exactly what
the user already expressed in the
`feedback_unify_after_two_instances.md` discipline: don't unify
prematurely, build the second instance first. Hollow sphere /
annulus already exist in trajectory_resolvent as standalone
solvers; multi-region sphere already exists in trajectory_resolvent
+ fn_method's reflected slab. So the *second instance* test is
satisfied for both. Extending `MeshTemplate` to carry an optional
`regions: tuple[(outer, mat_id), ...]` (or equivalent) is therefore
justified. **But this extension is MORE WORK than the ≈ 80 %
coverage that already exists**, and the prior audit's
half-day commit (move + use it from cross-method adapters) does
not require it. Multi-region MeshTemplate is a follow-on, not a
prerequisite.

The continuous-α / `BC.params` extension is the missing **one
line** in the BC layer. `BC.partial = BC("partial", {"albedo":
x})` plus a `alpha_from_bc(bc) -> float` helper covers every α
value trajectory_resolvent needs, without growing `MeshTemplate`.

### 5.5 Cross-method harness implication

Q4 anticipated this and the trace confirms it:

The `_parse_notes_kv` in `tests/cross_method/adapters.py:441-454`
exists ENTIRELY because two case sets (reflected slab,
closed-sphere k_inf) cannot fit through `MeshTemplate` today:

- **Closed-sphere k_inf** — fits `MeshTemplate` perfectly, just
  with `bc_left = bc_right = BC.reflective` (which is what α=1
  means). The `_closed_sphere_params` notes parsing is pure
  symptom — reading
  `(sigma_t, sigma_s, nu_sigma_f, R_cm, alpha=1.0)` from a string
  when the same data could come from
  `(case.materials[0], case.mesh_template.{critical_dimension_cm,
  bc_*})`. This case set ALREADY has all the data structures
  available (`materials[0]` is a fully-populated `Mixture`); it
  just bypasses them because `case.registry_case = None` blocks
  the `_extract_1g_xs` path.

  **Fix is one line**: give the `CrossMethodCase` an inline
  `materials: dict[int, Mixture] | None` and `mesh_template:
  MeshTemplate | None` field, set them on the inline cases
  (alongside the existing `registry_case = None`), and have the
  helpers consult `case.registry_case.{materials, mesh_template}`
  OR `case.{materials, mesh_template}` as appropriate.

- **Reflected slab** — does NOT fit current `MeshTemplate`. It
  needs `(c_core, c_reflector, reflector_half_thickness_mfp)`
  which is a **multi-region geometry + per-region c**. Today
  the F_N reflected solver takes those three scalars directly;
  there is no production analog yet (no ORPHEUS production solver
  for multi-region critical-eigenvalue with that parametrisation).
  Multi-region `MeshTemplate` extension would carry it (3 zones:
  reflector ↔ core ↔ reflector with `mat_ids = (1, 0, 1)` and
  per-mat-id `c` derivable from `materials`). But this is the
  multi-region extension, not the half-day promotion.

  Until then, the reflected-slab adapter pulls
  `(c_core, c_reflector, Δ)` from notes — and that's the ONLY
  place this gap matters in the current cross-method test suite.

So the cross-method adapter simplification path is:

| Step | What it does                                                              | Δ adapter helpers     |
|------|---------------------------------------------------------------------------|------------------------|
| 1    | Promote `MeshTemplate` to `common/`                                       | (no behavior change)   |
| 2    | Add `materials`, `mesh_template` fields to `CrossMethodCase` for inline cases (closed-sphere k_inf) | drop `_closed_sphere_params` |
| 3    | Use `mesh_template.critical_dimension_cm` directly in trajectory_resolvent adapters instead of `truth_value / sigma_t` | drop `_sphere_R_truth_mfp`, `_slab_a_truth_mfp` (or simplify) |
| 4    | Add `BC ↔ α` mapper                                                       | drop adapter `alpha` constants in favour of `mesh_template.bc_*` |
| 5    | (later) Extend `MeshTemplate` for multi-region — covers reflected slab    | drop `_reflected_slab_params` |

Steps 1-4 land in the half-day commit the prior audit recommends.
Step 5 is the multi-region extension and waits for the
cylinder-Variant-α work or another prompt to add it.

After steps 1-4, the only `_parse_notes_kv` user is the reflected
slab. After step 5, no adapter helper parses notes; every input is
a typed field on `MeshTemplate` or `Mixture`.

### 5.6 Verdict

**Yes — the input layer IS unifiable, and `MeshTemplate` already
sits at the unification point**. For single-region, single-physics
cases (the vast majority of the cross-method corpus today), every
adapter's needed input is derivable from
`(case.materials, case.mesh_template)`. The fact that adapters
today re-derive `R_cm` from `truth_value / sigma_t` (instead of
reading `mesh_template.critical_dimension_cm`) is a SYMPTOM, not
a structural barrier — the field exists, it's just unused.

The Q5 finding does NOT introduce a new architectural slot; it
**reinforces** the prior audit's recommendation to promote
`MeshTemplate` to `common/`. The promotion is the same half-day
commit. The new colour the input-level analysis adds is:

- **Step 2 (inline-case fields on `CrossMethodCase`)** is the
  highest-leverage individual change. Today the closed-sphere
  k_inf case is forced through a notes string for no structural
  reason; promoting it to typed `materials` + `mesh_template`
  fields kills `_closed_sphere_params` in one stroke.
- **Steps 3-4 (consume `MeshTemplate.critical_dimension_cm` and
  `BC ↔ α` directly)** make the trajectory_resolvent adapters
  trivial — they become "extract XS, extract domain extent,
  extract α from BC, run the solver, return scalar". No special
  helpers, no truth-value re-derivation.

Step 5 (multi-region MeshTemplate) is a follow-on once the
cylinder-Variant-α work or another multi-region instance lands.

### 5.7 Honest divergence note (what is NOT unifiable at the input level)

Three classes of input genuinely cannot live in a shared case-spec:

1. **Computational parameters** (`n_modes`, `n_r`, `n_mu`,
   `n_traj_quad`, `max_iter`, `tol`, ...). These belong on the
   adapter, not the case. The cross-method protocol already
   handles this correctly — the adapter declares its defaults; the
   case declares the *required tolerance* (`tolerances` map) and
   the adapter scales its computational parameters accordingly.
   Forcing them onto `MeshTemplate` would be a category error
   (physics-spec ≠ numerics-spec).

2. **Method-specific intermediates** (F_N's
   collocation grid `ξ`, trajectory_resolvent's `r_nodes` /
   `mu_nodes`, fn_method reflected's iteration scheme). These are
   *outputs* of solver-internal setup, not inputs. They have no
   home in a shared spec.

3. **Method-specific inputs without a physical-spec analog**
   (trajectory_resolvent's `initial_psi` for warm-starts, F_N's
   `bracket = (a_min, a_max)` for the dispersion root-find).
   These are computational hints that sit alongside computational
   parameters. Adapter constants today; could become adapter
   `__init__` arguments tomorrow; do NOT belong on the case.

The asymmetry between **dataclass** and **input** is therefore:

- **Dataclass**: unification is rejected (Q1-Q3 finding stands —
  Mesh1D ≠ trajectory_resolvent geometry args ≠ fn_method directory
  dispatch).
- **Input**: unification is accepted at the *physical-spec* layer
  (geometry kind, dimensions, materials, BC). The
  *computational-spec* layer stays per-adapter, by design.

This is a cleaner partition than the dataclass-only audit
suggested. The user's instinct — "the input might still be the
same even when the dataclass isn't" — is correct, with the
qualifier that "input" splits into physical and computational
sub-layers and only the physical sub-layer unifies.

### 5.8 No new refactor slot — confirms prior recommendation

The prior audit's "R0.5: promote `MeshTemplate` to `common/`" is
the ONE change that captures the input-level unification too. No
new pre-R1 step is needed beyond that promotion + the `materials
| mesh_template` fields on `CrossMethodCase` for inline cases (a
few extra lines in the same commit).

Verification that this doesn't fight any other architectural
direction:

- **Does NOT block R1/R2/R3** (chord oracle, power iterate, phase-
  space mesh). Orthogonal layers.
- **Does NOT block multi-region MeshTemplate** (Step 5 above) —
  that extension is additive.
- **Does NOT force production solvers to change**. They already
  consume `Mesh1D` via `mesh_template.build(n_cells)`. They
  remain unchanged.
- **Does NOT force trajectory_resolvent or fn_method internals to
  change**. They keep their per-method input signatures. ONLY the
  cross-method **adapter** layer (and the few inline cases'
  representation) shifts.
- **Does NOT close off the partial-α extension**. `BC.params =
  {"albedo": x}` is the natural slot; the BC system already
  supports it.

### 5.9 Files touched by this audit (read scope, no edits)

In addition to the prior-audit files, this follow-up read:
- `/workspaces/ORPHEUS/orpheus/cp/solver.py:704-749`
  (`solve_cp` signature — confirms Mesh1D + materials dict input)
- `/workspaces/ORPHEUS/orpheus/sn/solver.py:508-568`
  (`solve_sn` signature — same)
- `/workspaces/ORPHEUS/orpheus/derivations/continuous/sood_registry/builders.py`
  (the four-line build_materials/build_mesh/build_cp_params
  functions — confirms how the production-side bridge actually works)
- `/workspaces/ORPHEUS/orpheus/derivations/continuous/sood_registry/extractors.py`
  (`mixture_to_fn_arrays` — the only XS-shape adapter that already
  lives at the right layer)
- `/workspaces/ORPHEUS/tests/cross_method/cases.py:281-446`
  (REFLECTED_SLAB_CASES + CLOSED_SPHERE_KINF_CASES — the two
  inline case sets that bypass `MeshTemplate` today)
- `/workspaces/ORPHEUS/orpheus/derivations/continuous/trajectory_resolvent/greens_function.py:307-461`
  (`solve_greens_function_sphere` — confirms the input shape
  `(R, σ_t, σ_s, νσ_f, alpha, n_r, n_mu, n_traj_quad, max_iter,
  tol, initial_psi)`)
- `/workspaces/ORPHEUS/orpheus/derivations/continuous/fn_method/slab/one_group.py:121-260`
  (`solve_fn_slab_bare_critical` — confirms the input shape
  `(c, n_modes, a_min, a_max, n_bracket, bisect_tol, max_bisect)`)
- Nexus `callers` results for `Mesh1D`,
  `mesh1d_from_zones`, `solve_greens_function_sphere`,
  `solve_fn_slab_bare_critical`, `solve_fn_sphere_bare_critical`,
  `MeshTemplate.build` — confirms the upstream caller distribution
  used in 5.1.
