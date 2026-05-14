---
name: issue-196-phase-g-step2-5c-immutability-audit
description: Step 2.5c immutability audit — layered enumeration of immutable quantities, BC composition analysis, CellUpdate API recommendation post-cache, and a SweepCoefficientCache design proposal. Profiling target: 80% time in iter_cell_visits + per-cell StreamingTerms construction, all of which is σ_t-immutable once the problem is bound.
metadata:
  type: project
---

# Issue #196 Phase G Step 2.5c — Immutability Audit

User directive: "There is literally no need to be computing stuff that is
immutable once the problem is set up." Profiling: scan is 1% of sweep
time; 80% is `iter_cell_visits` + per-cell `StreamingTerms` dataclass
construction — every byte of which is invariant under σ_t. The audit
below enumerates what is invariant at each layer and proposes a
single-tensor `SweepCoefficientCache` whose population reduces the
inner-iteration sweep to ~3 numpy ops per ordinate.

---

## Q1 — Layered immutability table

Each row lists ONE quantity, its current home (file:line where it is
computed), and how it is currently surfaced to the per-cell update.

### Layer A — Geometry-immutable (invariant across the entire run)

| Quantity | Shape | Home (file:line) | How produced today |
|---|---|---|---|
| `mesh.edges` | (N+1,) | `orpheus/geometry/mesh.py:237` (frozen field) | Set at construction. |
| `mesh.widths` (=`dx`, `dr`) | (N,) | `orpheus/geometry/mesh.py:293` `@property` → `np.diff(edges)` | Recomputed on every access (no `cached_property`). |
| `mesh.centers` | (N,) | `orpheus/geometry/mesh.py:299` | Same — recomputed each access. |
| `mesh.volumes` | (N,) | `orpheus/geometry/mesh.py:314` — `precomputed_volumes` OR `compute_volumes_1d` | One-shot via `Mesh1D.from_geometry`. |
| `mesh.surfaces` (`face_areas`) | (N+1,) | `orpheus/geometry/mesh.py:321` | Computed each access via `compute_surfaces_1d`. |
| `reduced.face_areas` | (N+1,) | `orpheus/geometry/reduced_operator.py:585,679` | Computed ONCE in factory, stored on dataclass. |
| `reduced.delta_A` | (N,) | `orpheus/geometry/reduced_operator.py:588,680` | Stored. |
| `reduced.alpha_half` (sphere) | (N_ord+1,) | `orpheus/geometry/reduced_operator.py:595-602` — Bailey-Hébert α-dome recursion | Stored. |
| `reduced.redist_dAw` (sphere) | (nx, N_ord) | `orpheus/geometry/reduced_operator.py:605` — `delta_A[:, None] / w[None, :]` | Stored. |
| `reduced.tau_mm` (sphere) | (N_ord,) | `orpheus/geometry/reduced_operator.py:614-618` — M-M clamp | Stored. |
| `reduced.alpha_per_level` (cylinder) | list[(M+1,)] | `orpheus/geometry/reduced_operator.py:685-693` | Stored. |
| `reduced.redist_dAw_per_level` (cylinder) | list[(nx, M)] | `orpheus/geometry/reduced_operator.py:697-702` | Stored. |
| `reduced.tau_mm_per_level` (cylinder) | list[(M,)] | `orpheus/geometry/reduced_operator.py:713-730` | Stored. |
| `sn_mesh.dx`, `sn_mesh.dy` | (nx,), (ny,) | `orpheus/sn/geometry.py:197-205` | Set at `__init__`. |
| `sn_mesh._volumes` | (nx, ny) | `orpheus/sn/geometry.py:200` | Set at `__init__`. |
| `sn_mesh.streaming_x[n,i] = 2|μ_x|/dx` | (N_ord, nx) | `orpheus/sn/geometry.py:845-847` — `_setup_cartesian` | Stored. **Only 2-D Cartesian uses it.** Slab/sphere/cylinder leave it as `None`. |
| `sn_mesh.streaming_y[n,j]` | (N_ord, ny) | `orpheus/sn/geometry.py:849-851` | Stored. |
| `sn_mesh.sweep_graphs` (2-D) | dict[OctantLabel, SweepDependencyGraph] | `orpheus/sn/geometry.py:865-871` | Stored. |

### Layer B — Quadrature-immutable

| Quantity | Shape | Home (file:line) | Notes |
|---|---|---|---|
| `quad.mu_x, mu_y, mu_z` | (N,) | `orpheus/sn/quadrature.py:230,283,349,437` (cached in `__init__`) | View of `measure.nodes`. |
| `quad.weights` | (N,) | same | |
| `quad.N` | int | same | |
| `quad.level_indices` (LS/PQ) | list[(M,) int] | `orpheus/sn/quadrature.py:363,451` | Cylindrical level structure. |
| `quad.octants` | tuple | `_OctantsMixin.octants` `@cached_property` `orpheus/sn/quadrature.py:113-142` | First-access cached. **The only quadrature quantity using `cached_property`.** |
| `quad._ref_x, _ref_y, _ref_z` | (N,) int | `orpheus/sn/quadrature.py:287-289,354-356,442-444` | Reflection-partner indices, set at `__init__`. |
| `weight_norm = 1 / sum(w)` | float | `orpheus/sn/solver.py:227` | Set at `SNSolver.__init__`. **Also recomputed inside both sweep bodies at `sweep.py:226,411`.** Twin-path smell. |

### Layer C — `StreamingTerms` per-(cell, direction) — fully geometry-immutable

`StreamingTerms` carries: `chord_length`, `mu`, `face_area_inner`,
`face_area_outer`, `delta_A_over_w`, `alpha_in`, `alpha_out`, `tau_mm`,
`volume`, `abs_mu`. **All 10 fields are functions of `(mesh, quadrature)`
alone — nothing depends on σ_t, scattering, or BC.**

Current home: `orpheus/geometry/reduced_operator.py:385-505`
`ReducedStreamingOperator.streaming_terms(cell_idx, direction_idx, mu_level_idx)`.

**Current invocation pattern — the inefficiency**:
- `SNMesh._iter_cartesian_visits` (`orpheus/sn/geometry.py:573-598`)
  calls `self.reduced.streaming_terms(cell_idx=i, direction_idx=ordinate_idx)`
  **once per (cell, ordinate)**.
- Per call: builds a fresh frozen-slotted `StreamingTerms` dataclass
  with 10 float fields wrapped from already-stored arrays.
- The slab branch lookup `reduced_operator.py:419-446` is a 30-line
  `if` arm that re-extracts `mesh.widths[i]`, `mesh.volumes[i]`,
  `quadrature.mu_x[direction_idx]`, hard-codes 5 neutral curvature
  literals, and constructs the dataclass.
- Sphere/cylinder branches re-extract `face_areas[i]`,
  `face_areas[i+1]`, `redist_dAw[i, n]`, `alpha_half[n]`,
  `alpha_half[n+1]`, `tau_mm[n]` per call.
- `iter_cell_visits` then re-resolves `face_area_downstream` per visit
  by reading `st.face_area_outer` or `st.face_area_inner` based on the
  sweep direction sign (`geometry.py:600-629`, `_iter_spherical_visits`).

Every numerical value referenced is already in storage in `reduced.*`
or `quad.*`. The work is index gather + dataclass allocation, executed
once per (cell × ordinate × outer iteration × inner iteration). At
`N_ord = 16, nx = 160, ng = 2` that's `16 × 160 × 4096` = 10.5M
allocations per converged eigenvalue problem.

`visits = list(sn_mesh.iter_cell_visits(...))` at `sweep.py:274,497`
materialises one such list per `(ordinate, sweep call)`, then the
strategy's `affine_coefficients` (`diamond.py:208-247`) **re-extracts
every scalar field** via 9 `np.fromiter` calls. The `StreamingTerms`
dataclass is built only to be torn apart by the next layer.

### Layer D — σ_t-immutable (invariant across outer iterations and BCs)

| Quantity | Shape | Home | Currently computed |
|---|---|---|---|
| `total_xs * V` (collision contribution) | (nx, ng) | inside `diamond.affine_coefficients` `:261` | **Per ordinate per sweep call**. Same value for every direction. |
| `2|μ| · A_down` | (nx,) — depends on direction-sign class | inside `diamond.affine_coefficients` `:258` | **Per ordinate per sweep call**. Depends only on (|μ_n|, A_face_seq), not on σ_t or source. |
| `dA_w · c_out` (curvature redistribution part of denom) | (nx, ng) | inside `diamond.affine_coefficients` `:259` | **Per ordinate per sweep call**. `c_out = α_out/τ` purely from quadrature. |
| `denom` (full) | (nx, ng) | inside `diamond.affine_coefficients` `:257-261` | Recomputed inside every `affine_coefficients` call. **σ_t × V is the only iteration-variant piece** when σ_t actually changes; if σ_t is fixed across inner iterations (it is — `self.sig_t` is bound at `solver.py:186`), the **entire denominator is iteration-invariant**. |
| `inverse_denom = 1 / denom` | (nx, ng) | NEVER precomputed | Implicit in `a = 2|μ|·A_total/denom - 1` and `b = 2·(...)/denom`. Recomputed per call. |
| `a_attenuation[n,i,g] = 2|μ_n|·A_total[i]/denom[n,i,g] - 1` | (N_ord, nx, ng) | inside `affine_coefficients` `:276` | Recomputed per ordinate per sweep call. Pure function of (μ, geometry, σ_t). |
| `face_area_total[i] = A_inner[i] + A_outer[i]` | (nx,) | inside `affine_coefficients` `:275` | Recomputed each call. |
| `c_in = (1-τ)/τ · α_out + α_in` | (N_ord,) per geometry | inside `affine_coefficients` `:252` | Recomputed each call. |
| `chain_cell_idx` per direction sign (slab) | (nx,) int | `sweep.py:274` — `list(sn_mesh.iter_cell_visits(...))` | Re-materialised once per (ordinate, sweep call). Slab depends only on direction sign. |
| `bc_left_face / bc_right_face` shape (N, ng) | (N, ng) | `psi_bc` dict `sweep.py:236-237` | Persisted across sweeps. Inflow VALUE is per-sweep; the BUFFER + the reflective permutation map are σ_t-immutable. |

### Layer E — per-outer-iteration

Quantities that change across outer (power) iterations but are invariant
across inner Krylov / SI iterations within one outer step.

| Quantity | Shape | Home | When invalidated |
|---|---|---|---|
| `fission_source = χ · (νΣ_f·φ) / k` | (nx, ny, ng) | `solver.py:268-279` | Each outer power step. |
| `Q_bar_iso` (Carlson seed source) | (ng, nx) | `sweep.py:426-430` | Read from `phi_0_prev` cache; updated at end of each sweep (`sweep.py:586`). |
| `psi_bc["psi_pole"]` (Carlson pole-face IC) | (N, ng) | `sweep.py:585` | Updated end of each sweep. |
| `phi_aux` (Carlson seed) | (ng, nx) per level | `sweep.py:457-462` | Recomputed per sweep, but a pure function of (Q_bar, σ_t, dr, bc_outer). |

### Layer F — per-inner-iteration

| Quantity | Shape | Home |
|---|---|---|
| isotropic source `Q + S·φ + F·φ/k + (n,2n)` | (nx, ny, ng) | `solver.py:401-403` / `_solve_fixed_source_si:1095` |
| anisotropic source `Q_aniso = build_aniso(ψ_prev)` | (N, nx, ny, ng) or None | `solver.py:406` |
| `b = 2(source + dA_w·c_in·ψ_a_in) / denom` | (nx, ng) per ordinate | `diamond.affine_coefficients:277` |

### Layer G — per-sweep / per-direction (BC inflow)

| Quantity | Shape | Home |
|---|---|---|
| `inflow_left = bc_left.apply(bc_left_face)` | (N, ng) | `sweep.py:252` |
| `inflow_right = bc_right.apply(bc_right_face)` | (N, ng) | `sweep.py:253` |
| 2-D `full_face_x = sn_mesh.bc_xmin.apply(psi_x[:, 0, :, :])` | (N, ny, ng) | `sweep.py:729` |

---

## Q2 — Precomputation opportunity matrix

The σ_t-immutable rows of Layer D are the precomputation prize. For
the canonical `N_ord = 16, nx = 160, ng = 2` problem at float64:

| Quantity | Natural home | Shape | Storage | Population cost (one-shot) | Invalidation trigger |
|---|---|---|---|---|---|
| `chord_lengths` | `SweepCoefficientCache` | (nx,) | 1.3 kB | Already in `mesh.widths`. | Mesh change. |
| `volumes` | cache | (nx,) | 1.3 kB | Already in `mesh.volumes`. | Mesh change. |
| `face_areas` (curvilinear) | cache | (nx+1,) | 1.3 kB | Already in `reduced.face_areas`. | Mesh change. |
| `face_area_total[i]` | cache | (nx,) | 1.3 kB | one `+` over (nx,). | Mesh change. |
| `abs_mu[n]`, `signed_mu[n]` | cache | (N,) | 130 B | Already in `quad.mu_x`. | Quad change. |
| `face_area_downstream[n, i]` (slab=1.0, sphere/cyl=A_in or A_out by direction sign) | cache | (N, nx) | 21 kB | (N,) gather of `A_in/A_out`. | Mesh+quad change. |
| `c_in[n], c_out[n]` (slab=0.0) | cache | (N,) | 260 B | Pure-quadrature scalars. | Quad change. |
| `delta_A_over_w[n, i]` | cache | (N, nx) | 21 kB | Already in `reduced.redist_dAw`. | Mesh+quad change. |
| `streaming_face[n, i] = 2|μ_n|·A_down[n,i]` | cache | (N, nx) | 21 kB | One outer-product. | Mesh+quad change. |
| `curvature_redist[n, i] = dA_w[n,i] · c_out[n]` | cache | (N, nx) | 21 kB | Two broadcast mults. | Mesh+quad change. |
| **`denom[n, i, g] = streaming_face[n,i] + curvature_redist[n,i] + σ_t[i,g]·V[i]`** | cache | (N, nx, ng) | **41 kB** | broadcast sum over (N, nx, ng). | **σ_t change ONLY.** |
| **`inverse_denom[n, i, g] = 1/denom`** | cache | (N, nx, ng) | 41 kB | reciprocal. | σ_t change. |
| **`a_attenuation[n, i, g] = 2|μ_n|·A_total[i] · inverse_denom[n,i,g] − 1`** | cache | (N, nx, ng) | 41 kB | one fma. | σ_t change. |
| **`cumprod_a[n, i, g] = prod_{k≤i} a[n,k,g]`** (in chain order, per direction sign) | cache | (N, nx, ng) | 41 kB | one `np.cumprod`. | σ_t change. |
| `chain_idx[sign]` (slab direction-sign cell order) | cache | (nx,) int, two copies | 2.6 kB | reversed `arange`. | Mesh change. |
| `chain_idx[n]` (curvilinear, per ordinate) | cache | (N, nx) int | 21 kB | populated from current `iter_cell_visits` once. | Mesh+quad change. |

**Total σ_t-cache size at (N=16, nx=160, ng=2): ~280 kB.**

Population cost: dominated by one `np.cumprod` over (N, nx, ng) plus
a handful of broadcasts — single-digit milliseconds, paid once when
σ_t is bound. The per-inner-iteration savings are the entire content
of `_iter_cartesian_visits`, `_iter_spherical_visits`,
`_iter_cylindrical_visits`, `ReducedStreamingOperator.streaming_terms`,
and `DiamondDifference.affine_coefficients` — i.e. the 80% identified
by profiling.

### What about an N-dependent (per-outer-iteration) cache?

`Q_bar_iso` (Carlson seed source) and `psi_pole` IC are functions of
the **previous-iteration scalar flux**. They can live on a separate
per-outer-iteration cache (a `CarlsonSeedCache` with one method
`refresh(phi)`), populated at the start of each sweep. The Carlson
inward sweep itself is already a 6-line tensor recurrence
(`psi_half_angle_seed.py:411-419`) and is fast enough to leave alone.

### What about the per-inner cache?

`b = 2·(source + ang_contrib) · inverse_denom` reduces to two
broadcasts + one multiply once `inverse_denom` and `dA_w·c_in[n,i]`
are precomputed. The `b` tensor is per-inner-iteration but evaluated
in O(N·nx·ng) tensor time, not per-cell Python time. No caching
needed beyond Layer D.

---

## Q3 — BC composition analysis

### Anatomy of `bc.apply(face_buffer)` today

The realized BC operators after Issue #188 / #186 are all
**single-arg, axis-tagged LinearOperators**:

- **Vacuum** → `IncomingOrdinateMaskTensor(inflow_indices, n_ordinates, axis=0)`
  (`orpheus/sn/boundary_realizer.py` mapping; impl
  `orpheus/numerics/operator.py:843` family). `apply` is one
  `np.take`-style mask write — one allocation, O(N·ng) work for a
  1-D face buffer.

- **Reflective** → `PermutationOperator(reflection_index('x'),
  axis=0)`. `apply` is **`np.take(x, self.perm, axis=self.axis)`**
  (`operator.py:823-824`) — one numpy gather. The permutation array
  is built ONCE in `quadrature.create()`
  (`quadrature.py:287-289,354-356,442-444`) and never recomputed.

- **Mixed faces** → `OperatorSum` etc.; all O(N·ng) single-pass tensor
  ops.

**`bc.apply(...)` is already a fully-vectorised tensor op.** There is
nothing per-cell or per-ordinate-Python inside it. The whole BC
composition costs O(N·ng) per face per sweep — negligible vs the
solver's other work.

### Could BC be folded into the σ_t cache?

For **vacuum** + **reflective**: yes, partially.

- Vacuum BC reads no state. `apply(zeros) → zeros` at the inflow
  ordinates. The inflow value at each (n, sweep) is a constant zero;
  this is structurally absorbed by initialising `psi_face_in = 0` for
  vacuum-inflow ordinates and skipping the `bc.apply(...)` call.
- Reflective BC reads the **PREVIOUS sweep's outgoing flux** at the
  partner ordinate (`sweep.py:236-239,311-313`). This is a
  per-sweep dependency by physical necessity — the Picard / Gauss-Seidel
  fixed-point that source-iteration converges to. It cannot be
  precomputed.

**Verdict**: BC apply is fundamentally a per-sweep operation
(reflective dependency). It is also already a single tensor op. No
gain from caching; no friction from leaving it as-is. The cleanest
sweep skeleton calls `bc.apply(face_buffer)` once per face per sweep
direction and feeds the result into the per-ordinate scan.

The 2-D Cartesian path (`sweep.py:728-740`) already collapsed the
N-fold per-ordinate BC apply down to one per-octant. Step 2.5c
inherits that pattern.

---

## Q4 — CellUpdate API surface post-cache

The user's question is whether `CellUpdate.affine_coefficients(visits, ...)`
collapses to a cache populator. Yes, with one subtlety.

### What survives

**`CellUpdate.update(visit, ...)`** is the human-legible reference
implementation of the per-cell algebra and the partner of the
`residual(...)` dual-view. The Pattern 2 contract
(`diamond.py:283-314`, `cell_balance.py`) is "update and residual share
the algebra via `cell_balance_terms`". This Pattern 2 architecture is
what keeps the L+C operator-algebra correctness gate intact. **Keep it.**
The slow Python `update(...)` becomes the L1 dual-view validator: for
any cache-driven sweep result `ψ_avg[n, i, g]`, the per-cell residual
evaluated through `update(...).residual(...)` must be zero to ULP.

### What goes

**`CellUpdate.affine_coefficients(visits, total_xs, source, angular_state)`**
in its current form is doing TWO orthogonal jobs:

1. **σ_t-immutable coefficient assembly** (the `denom`, the
   `2|μ|·A_total/denom − 1` ratio that becomes `a`). This is what the
   cache subsumes.
2. **σ_t- and source-VARIANT b-vector assembly** (`b = 2·(source +
   dA_w·c_in·ψ_a_in)·inverse_denom`). This stays per-inner-iteration.

After the cache, the call collapses to:

```python
# inside the sweep, per-ordinate, per-inner-iteration:
b = 2.0 * (source_chain + dA_w_c_in[n][:, None] * psi_a_in_chain) * cache.inverse_denom[n_chain]
psi_face = cache.cumprod_a[n_chain] * (psi_in + np.cumsum(b / cache.cumprod_a[n_chain], axis=0))
```

Two natural API designs:

**(A) Cache populator method on the strategy**
```python
class DiamondDifference(...):
    def populate_cache(self, sn_mesh: SNMesh, sig_t: np.ndarray) -> "SweepCoefficientCache":
        ...
```
Strategy-specific: LD / EC / Step would each implement their own
populator returning a strategy-specific cache shape. Locality of
implementation is good (the math lives next to `update`/`residual`).

**(B) Tensor-form coefficient assembly + free `ordinate_scan`**
```python
class DiamondDifference(...):
    def coefficient_tensors(
        self, sn_mesh: SNMesh, sig_t: np.ndarray,
    ) -> tuple[ndarray, ndarray, ndarray]:
        """Return (a_attenuation, inverse_denom, cumprod_a)."""
```
Lower-level, three named tensors instead of a dataclass.

**Recommendation: (A).** `SweepCoefficientCache` as a frozen dataclass
returned by `cell_update.populate_cache(sn_mesh, sig_t)`. Its fields
ARE the named tensors; consumers read them by name; the strategy is
the single source of truth for cache shape (Pattern 2 dual-view stays
clean — `update` and `populate_cache` derive from the same
`cell_balance_terms` algebra).

The current `affine_coefficients(visits, ...)` then becomes a thin
adapter on the cache (vectorised gather along the chain axis) — kept
for the dual-view test, slower-than-cache but algebraically identical.

---

## Q5 — Curvilinear extras

The curvilinear sweep has three angular pieces beyond the slab
spatial-only structure. Each has a different invariance class.

### Carlson coupled-pole seed

`Q_bar = σ_t · φ_0_prev / Σw_full` (`sweep.py:425-430`) is
per-outer-iteration. The inward μ=−1 sweep itself
(`psi_half_angle_seed.py:411-419`) is a vectorised-over-groups
sequential-in-cells loop of ~6 lines. Total cost: O(nx·ng) per sweep.
**Already tensor-shaped along the group axis; sequential along the
spatial axis by design.** Recommend: leave alone, move to a
per-outer-iteration `CarlsonSeedCache` only if profiling demands.

### M-M angular redistribution

The Morel-Montry angular closure `ψ^a_out = (ψ_avg − (1−τ)·ψ^a_in)/τ`
appears two places:

- In the per-cell `update(...)` (`diamond.py:166-170`).
- In the vectorised sweep skeleton, where it's already pulled out as
  a per-chain numpy expression (`sweep.py:555-562`):
  ```python
  tau_chain = np.fromiter(...)[:, None]
  psi_angle_out_chain = (psi_avg_chain - (1 - tau_chain) * psi_angle_in_chain) / tau_chain
  ```

`(1-τ)/τ` and `1/τ` are pure-quadrature constants. **They belong on
the cache as `mm_a_in_coeff[n] = (1-τ[n])/τ[n]` and `mm_inv_tau[n] =
1/τ[n]`.** Then the chain step is `psi_angle_out = mm_inv_tau[n] · psi_avg − mm_a_in_coeff[n] · psi_angle_in`.

### Angular-state thread across ordinates within a level

`psi_angle[i]` is threaded ordinate-by-ordinate within one μ-level
(`sweep.py:469,540,566`). Each ordinate's chain reads
`psi_angle[chain_cell_idx]` (from the previous ordinate's chain write)
and writes back the M-M outgoing value. This is the **angular DAG**
that runs orthogonally to the spatial DAG.

Vectorisability: this is the second affine recurrence in the system
(the first is the spatial one). Per ordinate, given the cell-average
chain `ψ_avg_chain`, the angular update IS a free function of
`ψ_angle_in_chain`. The full level loop, viewed as a 2-D problem with
axes (m ∈ {0..M-1} ordinates, i ∈ {0..nx-1} cells), is a sequence of
M coupled spatial chains. The level **cannot be done as a single
ordinate_scan** because each ordinate's `a, b` depends on the previous
ordinate's `ψ_a_out` (the angular dependence couples the chains).

But ψ_a_in enters the chain only through `b` (additively, via
`ang_contrib = dA_w · c_in · ψ_a_in`). So the per-level structure
splits into two interleaved scans:

1. **Spatial scan** for ordinate m given ψ_a_in_chain[m]:
   `ψ_face_chain[m] = cumprod_a[m] · (ψ_in + cumsum((b_source_chain + dA_w·c_in[m]·ψ_a_in_chain[m]) / cumprod_a[m]))`
2. **Angular update** for ordinate m+1:
   `ψ_a_in_chain[m+1] = (ψ_avg_chain[m] − mm_a_in_coeff[m]·ψ_a_in_chain[m]) / τ[m]`

This is a **rank-1 augmented recurrence** — exactly the structure
`ordinate_scan` decomposes for spatial. The full level update is the
pair-monoid composition along the M axis of (spatial-scan, M-M-update)
pairs. The work-efficient implementation is sequential along M with
each M-step being O(nx·ng) tensor work. **Same fold-with-accumulator
pattern Step 2.5 already established for the angular axis** — no new
mathematical structure, just precompute the per-level `(τ, 1−τ)/τ`
constants.

---

## Q6 — Mathematical notation for the native expression

### Slab — the canonical Blelloch form

The recurrence:
```
ψ^s[n,i+1,g] = a[n,i,g] · ψ^s[n,i,g] + b[n,i,g]
```
with
```
a[n,i,g] = (2|μ_n| − σ_t[i,g]·V[i]) / (2|μ_n| + σ_t[i,g]·V[i])
b[n,i,g] = 2·q[n,i,g] / (2|μ_n| + σ_t[i,g]·V[i])
```
where for slab V[i] = dx[i]. (`A_total = 2`, `A_down = 1`, so
`2|μ|·A_total/denom − 1 = (2|μ|·2 − denom)/denom = (4|μ| − 2|μ| − σ_t·dx)/(2|μ| + σ_t·dx)`,
which simplifies to the form above after writing `denom = 2|μ| + σ_t·dx`
for slab.)

Closed-form solution (Blelloch 1990 §1.5, already implemented as
`ordinate_scan`):
```
ψ^s[n, i, g] = cumprod_a[n,i,g] · (ψ_0[n,g] + cumsum_{k≤i}(b[n,k,g] / cumprod_a[n,k,g]))
```

### Transport-literature naming

The pair `(a, b)` is canonically called the **transmission–emission**
pair (Lewis & Miller 1984 §5.3 in the context of DD; Larsen-Adams
SAILOR framework calls them **per-cell sweep coefficients**). MOC
calls the analog the **track-segment characteristic exponential** `e^(-τ)`
and the **emission term** `(1 − e^(-τ))·Q/Σ_t`. The slab DD pair is the
finite-volume analog of MOC's `(e^(-τ), Q/Σ_t·(1−e^(-τ)))` — same
algebraic structure, different closure.

For curvilinear with the curvature redistribution, the denominator
acquires the M-M term and the additive term acquires `dA_w·c_in·ψ_a_in`.
The notation generalises naturally:
```
a[n,i,g] = 2|μ_n|·A_total[i] / denom[n,i,g] − 1
b[n,i,g] = 2·(q[n,i,g] + dA_w[n,i]·c_in[n]·ψ_a_in[n,i,g]) / denom[n,i,g]
denom[n,i,g] = 2|μ_n|·A_down[n,i] + dA_w[n,i]·c_out[n] + σ_t[i,g]·V[i]
```
Slab degenerates to the canonical form by `A_total=2, A_down=1, dA_w=0`.

### Does the code read like the math?

Today: no. The slab path reads `a, b = cell_update.affine_coefficients(visits, sig_t_chain, QV_chain, None); psi_face_chain = ordinate_scan(a, b, psi_spatial_in)` — which is close, but `affine_coefficients` is a 70-line dataclass-tear-down that builds `a` and `b` from per-cell scalars via 9 `np.fromiter` calls (`diamond.py:212-247`).

After the cache: yes. The slab inner loop is exactly the math:
```python
b = 2.0 * source_chain * cache.inverse_denom[n_chain]
psi_face = cache.cumprod_a[n_chain] * (psi_in + np.cumsum(b / cache.cumprod_a[n_chain], axis=0))
psi_avg  = 0.5 * (psi_in_chain + psi_face)  # DD cell-average closure
```
Three numpy expressions, one per mathematical operation.

---

## Recommendation — `SweepCoefficientCache`

### Class design

```python
# new module: orpheus/sn/spatial/sweep_cache.py

@dataclass(frozen=True, slots=True)
class SweepCoefficientCache:
    """Per-(ordinate, cell, group) precomputed sweep coefficients.

    Built once when σ_t is bound to the solver. Invariant across all
    BC changes, source evaluations, inner Krylov steps, and outer
    power iterations (until σ_t changes — material rebuilds invalidate).
    """

    # Geometry-blind tensors — populated for every geometry, with
    # slab carrying neutral curvature values (cf. Step 2.5 StreamingTerms).

    # Chain (sweep-direction-resolved) cell orderings.
    chain_idx: np.ndarray            # (N_ord, nx) int — cell index in chain order per ordinate
    chain_idx_inverse: np.ndarray    # (N_ord, nx) int — inverse mapping

    # Per-(ordinate, cell) geometry quantities IN CHAIN ORDER.
    abs_mu: np.ndarray               # (N_ord,)
    A_down: np.ndarray               # (N_ord, nx)  (1.0 slab, 0.0 cyl-deg, face area curvilinear)
    A_total: np.ndarray              # (N_ord, nx)  (2.0 slab, A_in+A_out curvilinear)
    dA_w: np.ndarray                 # (N_ord, nx)  (0.0 slab, ΔA/w_n curvilinear)
    V: np.ndarray                    # (N_ord, nx)  cell volume in chain order
    tau_inv: np.ndarray              # (N_ord,)     1/τ_n   (slab: 1.0)
    mm_a_in_coeff: np.ndarray        # (N_ord,)     (1-τ_n)/τ_n  (slab: 0.0)
    c_in: np.ndarray                 # (N_ord,)     (slab: 0.0)
    c_out: np.ndarray                # (N_ord,)     (slab: 0.0)

    # σ_t-dependent tensors — invalidated only when σ_t rebinds.
    denom: np.ndarray                # (N_ord, nx, ng)
    inverse_denom: np.ndarray        # (N_ord, nx, ng)
    a_attenuation: np.ndarray        # (N_ord, nx, ng)
    cumprod_a: np.ndarray            # (N_ord, nx, ng)  precomputed along chain axis

    # Diagnostics / dual-view validation.
    is_degenerate: np.ndarray        # (N_ord,) bool — cyl-deg ordinates (skip scan)

    @classmethod
    def from_diamond_difference(
        cls,
        sn_mesh: SNMesh,
        sig_t: np.ndarray,           # (nx, ng) — slab/curvilinear 1-D
    ) -> "SweepCoefficientCache":
        ...
```

### Home

`orpheus/sn/spatial/sweep_cache.py`. Sits alongside `diamond.py`,
`cell_balance.py`, `scan.py`, `cell_update.py`. The cache populator
delegates per-cell to `cell_balance_terms(...)` then stacks results
into the tensor layout — Pattern 2 single source of truth preserved.

### Lifecycle

- **Populate**: `SNSolver.__init__` builds the cache once after
  binding `self.sig_t` (insert at `solver.py:260` near the L operator
  construction).
- **Invalidate**: any `replace(self.sig_t, new_value)` must rebuild
  the cache. Tag with a `sig_t_id = id(sig_t)` cross-check OR more
  robustly bind the cache on the solver and assert
  `solver._cache.sig_t_id == id(solver.sig_t)` in dev mode.
- **Live for**: every sweep call within one solver instance lifetime.

### Slab sweep skeleton, post-cache (≤ 30 LOC)

```python
def _sweep_1d_cartesian(Q, sig_t, sn_mesh, psi_bc, Q_aniso=None):
    cache = sn_mesh.cell_update.cache_for(sn_mesh, sig_t)  # cached if σ_t unchanged
    quad, N, nx = sn_mesh.quad, sn_mesh.quad.N, sn_mesh.nx
    ng = Q.shape[2]

    weight_norm = 1.0 / quad.weights.sum()
    QV_iso = (Q[:, 0, :] * sn_mesh.volumes[:, 0, None]) * weight_norm   # (nx, ng)
    has_aniso = Q_aniso is not None

    # BC inflow at both faces — one tensor op each.
    bc_left_face  = psi_bc.setdefault("bc_1d_left_face",  np.zeros((N, ng)))
    bc_right_face = psi_bc.setdefault("bc_1d_right_face", np.zeros((N, ng)))
    inflow_left   = sn_mesh.bc_left.apply(bc_left_face)    # (N, ng)
    inflow_right  = sn_mesh.bc_right.apply(bc_right_face)  # (N, ng)

    angular_flux = np.zeros((N, nx, 1, ng))
    scalar_flux  = np.zeros((nx, ng))

    for n in range(N):
        mu_n = quad.mu_x[n]
        chain = cache.chain_idx[n]
        QV_chain = (QV_iso[chain] if not has_aniso
                    else QV_iso[chain] + Q_aniso[n, chain, 0, :] * sn_mesh.volumes[chain, 0, None] * weight_norm)
        psi_in = inflow_left[n] if mu_n >= 0 else inflow_right[n]

        # ── The three native tensor ops ────────────────────────────
        b = 2.0 * QV_chain * cache.inverse_denom[n]                                       # (nx, ng)
        psi_face = cache.cumprod_a[n] * (psi_in + np.cumsum(b / cache.cumprod_a[n], axis=0))  # (nx, ng)
        psi_in_chain = np.concatenate(([psi_in], psi_face[:-1]), axis=0)
        psi_avg_chain = 0.5 * (psi_in_chain + psi_face)
        # ───────────────────────────────────────────────────────────

        psi_avg = psi_avg_chain[cache.chain_idx_inverse[n]]
        angular_flux[n, :, 0, :] = psi_avg
        scalar_flux += quad.weights[n] * psi_avg
        (bc_right_face if mu_n >= 0 else bc_left_face)[n] = psi_face[-1]

    return angular_flux, scalar_flux[:, None, :]
```

### Curvilinear sweep skeleton, post-cache

Same pattern; the only difference is the per-level loop and the
angular-state thread. The cache carries `mm_a_in_coeff[n]` and
`tau_inv[n]`. Per-level body becomes:

```python
phi_aux = carlson_inward_sweep_from_source(Q_bar, σ_t, dr, bc_outer_value)
psi_angle = phi_aux.T.copy()                                    # (nx, ng)

for m_local, global_n in enumerate(ordinates_in_level):
    chain = cache.chain_idx[global_n]
    if cache.is_degenerate[global_n]:
        ... # rare path; keep the per-cell Python loop with cell_update.update
        continue

    QV_chain = QV_iso[chain] + (Q_aniso[global_n, chain, 0, :] * V_chain * weight_norm if has_aniso else 0)
    psi_a_in_chain = psi_angle[chain]
    psi_in = inflow_full[global_n] if μ[global_n] < 0 else psi_pole[global_n]

    # ── Native tensor ops ─────────────────────────────────────
    ang_contrib = cache.dA_w[global_n] * cache.c_in[global_n] * psi_a_in_chain
    b = 2.0 * (QV_chain + ang_contrib) * cache.inverse_denom[global_n]
    psi_face = cache.cumprod_a[global_n] * (psi_in + np.cumsum(b / cache.cumprod_a[global_n], axis=0))
    psi_in_chain = np.concatenate(([psi_in], psi_face[:-1]), axis=0)
    psi_avg_chain = 0.5 * (psi_in_chain + psi_face)
    psi_a_out_chain = cache.tau_inv[global_n] * psi_avg_chain - cache.mm_a_in_coeff[global_n] * psi_a_in_chain
    # ──────────────────────────────────────────────────────────

    psi_angle[chain] = psi_a_out_chain
    angular_flux[global_n, chain, 0, :] = psi_avg_chain
    scalar_flux[chain] += w[global_n] * psi_avg_chain
    if μ[global_n] >= 0 and abs(μ[global_n]) > ε:
        bc_outer[global_n] = psi_face[-1]
```

### What this kills

- `_iter_cartesian_visits`, `_iter_spherical_visits`, `_iter_cylindrical_visits`
  → `cache.chain_idx[n]` (one array gather replaces a Python generator
  yielding `nx` dataclass instances per ordinate).
- `ReducedStreamingOperator.streaming_terms(...)` from the hot path →
  the cache absorbs everything; the method survives as the L1
  reference for cache population only.
- `DiamondDifference.affine_coefficients(...)` from the hot path →
  replaced by `cache.cumprod_a` lookup + per-inner b assembly.
- `iter_cell_visits(...).list(...)` allocations →
  the cache is pre-allocated tensors; zero per-sweep allocation in
  the spatial-coefficient path.

### What this preserves

- `cell_update.update(visit, ...)` as the dual-view reference for
  `residual` and as the per-cell L1 oracle (degenerate cylindrical
  cells still take this path).
- The Pattern 2 dual-view test
  (`tests/sn/spatial/test_cell_balance.py` + scan tests) — the cache
  is derived from `cell_balance_terms` so the dual-view assertion
  becomes `cache_driven_sweep(...) == cell_update.update(...).cell_average_flux`
  by construction.
- BC operator surface — already tensor, kept as-is.
- The Carlson seed and M-M angular thread — already tensor-shaped;
  augmented with cache constants but architecturally unchanged.

---

## Quick references

| Concern | Citation |
|---|---|
| `StreamingTerms` per-call allocation site | `orpheus/geometry/reduced_operator.py:385-505` |
| `affine_coefficients` per-call extraction | `orpheus/sn/spatial/diamond.py:208-247` |
| `cell_balance_terms` single-source algebra | `orpheus/sn/spatial/cell_balance.py:120-199` |
| Slab sweep current body | `orpheus/sn/sweep.py:172-315` |
| Curvilinear sweep current body | `orpheus/sn/sweep.py:322-588` |
| 2-D wavefront body | `orpheus/sn/sweep.py:605-773` |
| BC apply contract (1-arg LinearOperator) | `orpheus/geometry/boundary/_bound_compat.py:112-116` + `orpheus/numerics/operator.py:823-824` (`PermutationOperator.apply`) |
| Quadrature reflection-index precomputation | `orpheus/sn/quadrature.py:287-289` (cached at `create()`) |
| `ordinate_scan` Blelloch §1.5 form | `orpheus/sn/spatial/scan.py:136-137` |
| `SNSolver` σ_t binding | `orpheus/sn/solver.py:186` |
