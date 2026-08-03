# #326 — the ORDERING CONTRACT on `LevelStructure.level_indices`

Map of every production consumer of the per-level ordinate order, what each
one assumes, and what a half-range (quotient) level would have to satisfy.

- Branch `refactor/operator-strategy-layers`, HEAD `fadb827a`.
- Working tree at open AND at close: only `.claude/skills/vv-principles/*`
  modified + `scratch/_from_instruction_arch_worktree/` untracked. **No
  production file is mid-edit**, so every `file:line` below is at HEAD.
- All numeric claims below are MEASURED (probe transcripts inline), not read.

---

## 0. Premise check — the tree has moved AHEAD of the brief

Two corrections to the brief's coordinates, both from commit `3afb52c2`
("a level is a FIBER of an invariant"), **one commit before HEAD**:

1. **`LevelStructure` already carries the fiber coordinate and the successor
   ordering.** `orpheus/numerics/quadrature/rules_sphere.py:137-211` — the
   dataclass now has `polar_invariant: PolarInvariant` (`SIGNED_MU_Z` for
   `product_mu_phi`, `ABS_MU_Z` for `level_symmetric_sn`), `azimuth`,
   `hemisphere`, and a **`fiber(level)` accessor** returning
   `np.lexsort((azimuth[members], hemisphere[members]))`. Its docstring names
   #326 and states why it was NOT made the stored order: *"the stored order is
   what the cylindrical sweep consumes today, and changing it moves results
   (issue #326). This gives the correct ordering a name and a home first."*
   `fiber()` has **zero production callers** — `grep -rn '\.fiber(' orpheus/
   tests/` returns only `tests/numerics/test_rules_sphere.py:340,362`. It is a
   built-and-tested but unwired successor.
   The `level_indices` docstring already carries the `.. warning::` that the
   sort key is 2-to-1 on the fiber.

2. **`orpheus/sn/quadrature.py` DOES NOT EXIST.** The adapter surface named in
   the brief dissolved into `orpheus/numerics/quadrature/directional.py`.
   Stale `:class:` refs to `orpheus.sn.quadrature.ProductQuadrature` /
   `.LevelSymmetricSN` survive at `rules_product.py:8,91,130`,
   `rules_sphere.py:372`, `reduced_operator.py:741-742`, `measure.py:1259` —
   dangling Python-domain xrefs (the silent class: no `-W` warning without
   `-n`).

The two line numbers the brief cites are still correct:
`rules_product.py:172` and `rules_sphere.py:312`.

---

## 1. Consumer table

Every production site (`orpheus/`, non-test) that reads `level_indices` /
`LevelStructure`. "Order-sensitive" = would give a DIFFERENT answer under an
arbitrary permutation of a level's members.

| # | `file:line` | what it reads | order-sensitive? | why |
|---|---|---|---|---|
| **P1** | `numerics/quadrature/rules_product.py:170-173` | PRODUCES: `order = np.argsort(mu_x[level_arr])` | **PRODUCER** | key `η` is 2-to-1 on the fiber; `argsort` default `kind="quicksort"` is **not stable** ⇒ tie order unspecified |
| **P2** | `numerics/quadrature/rules_sphere.py:311-313` | PRODUCES: `idx = where(abs(mu_z)==level_mu[p])`; `argsort(mu_x[idx])` | **PRODUCER** | same key. Level is **4-to-1** here (both `mu_z` signs AND both `ξ` signs) |
| **P3** | `numerics/quadrature/rules_sphere.py:197-211` | `LevelStructure.fiber(level)` — the lexsort successor | n/a | **0 production callers** |
| **P4** | `numerics/quadrature/directional.py:501-511` | `Quadrature.level_indices` → `self.level_structure.level_indices` | pass-through, **not cached** | see §5 |
| **A1** | `geometry/reduced_operator.py:778-786` | **α recursion** `alpha[m+1] = alpha[m] - w[m]*eta[m]` | **YES — core** | §2 |
| **A2** | `geometry/reduced_operator.py:790-795` | `redist_dAw = ΔA[:,None]/w[level_idx][None,:]` → `(nx, M)` | **grouping only** | pure gather; position-covariant with everything else that reads `[p][m]` |
| **A3** | `geometry/reduced_operator.py:805-809` | `sin_theta = sqrt(1 - mu_z[level_idx[0]]**2)`; `mu_start = -sin_theta` | **grouping only (but see note)** | `mu_z` is CONSTANT on a level, so `[0]` is a level LABEL, permutation-immune. The **sign** (`-`) presumes the march STARTS at the inward edge — that is A1/A4's assumption, restated |
| **A4** | `sn/sweep/pole_angular_closure.py:594-611` | **`eta_edge` construction** | **YES — core** | §2 |
| **A5** | `sn/sweep/pole_angular_closure.py:664-675` | `morel_montry_tau_per_level` — clamps A4's raw τ to `[0.5, 1]` | **YES (it HIDES A4's damage)** | `max(0.5, min(1.0, τ_raw))` silently absorbs out-of-range τ. §2 |
| **A6** | `sn/mesh/augmented_mesh.py:811-816` | **R12a predicate**: `eps < float(tau_level[0]) < 1.0 - eps` | **YES — first element** | reads `τ_raw` at **position 0 of the level**. Which physical ordinate is at position 0 decides whether the mesh carries a ψ½ state block at all |
| **C1** | `sn/sweep/pole_angular_closure.py:904-906` | `self.level_indices = tuple(np.asarray(lvl) for lvl in quad.level_indices)` | **DESTRUCTURING BOUNDARY** | the closure SNAPSHOTS the order; §5 |
| **C2** | `sn/sweep/pole_angular_closure.py:926-935` | `alpha_in = alpha[:M_p]`, `alpha_out = alpha[1:M_p+1]`; `c_out = α_out/τ`, `c_in = (1-τ)/τ·α_out + α_in` | **YES (inherited)** | reads α at positions `m` and `m+1` as the two faces of one angular cell |
| **C3** | `sn/sweep/pole_angular_closure.py:262-271` | `_gather_per_ordinate`: `N = sum(lvl.size ...)`; `out[level] = values` | **grouping + COVERAGE** | assumes the partition covers `[0, N)`. **IndexError** on a non-covering partition — §4 |
| **C4** | `sn/sweep/pole_angular_closure.py:717-721` | **the M-M recurrence** `psi_half[m+1] = (psi_level[m] - (1-τ_m)·psi_half[m])/τ_m`, seeded at `psi_half[0]` | **YES — core** | marches strictly in POSITION order; position 0 is the seed |
| **C5** | `sn/sweep/pole_angular_closure.py:1019-1030` | `_edge_seed_stencil`: `order = np.argsort(mu)`; `m0 = order[0]`; first `cand` with `\|mu-mu[m0]\|>1e-14` | **order-INVARIANT by re-derivation** | it RE-SORTS. Returns *positions*, so the returned `m0/m1` track the permutation. The `1e-14` skip is the mirror-pair degeneracy hard-coded |
| **C6** | `sn/sweep/pole_angular_closure.py:1070-1084` | `precompute_psi_state`: `psi_level = psi_g_first[:, level_idx, :]` | **YES (inherited)** | the gather that defines "position m" for C4 |
| **C7** | `sn/sweep/pole_angular_closure.py:1184-1213` | `angular_adjoint`: reverse loop `for m in range(M-1,-1,-1)`, `psi_bar[:, level_idx[m], :] += phb/tau_m` | **YES (inherited)** | exact transpose of C4 — same contract, reversed |
| **C8** | `sn/sweep/pole_angular_closure.py:1111-1119` | `cell_contribution`: `self._dAw_per_level[p][cell, within_positions]`, `_c_in_per_level[p][within_positions]` | **positional API** | `within_positions` are within-level POSITIONS. Every caller must use the same order |
| **W1** | `sn/sweep/cache.py:246-254` | `level_visits_iter.append((p_idx, m_local, int(global_n)))` | **grouping + COVERAGE** | `chain_idx`/`A_down`/… are `np.empty((N, nx))`; only visited `global_n` rows get written — §4 |
| **W2** | `sn/mesh/augmented_mesh.py:1259-1262` | `global_n = int(level_indices[mu_level][ordinate_idx])` | positional decode | `ordinate_idx` is a within-level position — consistent as long as ONE order is used everywhere |
| **W3** | `sn/mesh/augmented_mesh.py:1315-1330` | `_representative_ordinate`: `cand = where(eta_at_level > ±eps)[0]`; `return int(cand[0])` | **grouping only** | any ordinate of the right sign gives the same cell chain (documented). Needs ≥1 member of each sign |
| **W4** | `sn/mesh/augmented_mesh.py:1468-1469` | `_make_cell_visit` global-ordinate lookup | positional decode | same as W2 |
| **L1** | `sn/loss_representation/__init__.py:2827-2852` | `_dag_legs`: `within_mask = mu_level > +eps` / `< -eps`; `ordinates = level_idx_arr[within_mask]`; `within = np.where(within_mask)[0]` | **grouping only** | mask + `np.where` are permutation-covariant. BUT the leg's `ordinates`/`within` inherit the level's order and feed C8 |
| **L2** | `sn/loss_representation/__init__.py:2901-2914` | `_degenerate_positions`: resolves each `\|μ_x\|<eps` global ordinate to `(level, position)` | **grouping + COVERAGE** | `deg_level`/`deg_within` are appended ONLY on a hit; a level that no longer contains a degenerate ordinate silently shortens the lists while `n_deg = global_deg.size` stays — §4 |
| **L3** | `sn/loss_representation/__init__.py:3518` | `numer_bar` tuple built `for li in closure.level_indices` | inherited | shape only |
| **L4** | `sn/loss_representation/__init__.py:4089-4091, 4210` | `_run`: `level_ordinates_list = [list(li) for li in level_indices]`; `for m_local, global_n in enumerate(...)` | **YES — the live sweep** | this loop IS the M-M angular thread order |
| **L5** | `sn/loss_representation/__init__.py:4154-4155, 4620-4621` | seed fold: `m0_local,_,t = closure._edge_seed_stencil(p)`; `m0_global = level_indices[p][m0_local]` | order-INVARIANT (via C5) | decodes C5's position through the SAME array |
| **L6** | `sn/loss_representation/__init__.py:4222-4238, 4303-4307` | **coupled-pole seed**: `psi_in = pole_outflow[mirror[global_n]]` for `mu_n>0`; captured at 4307 for `mu_n<0` | **YES — SIGN-PARTITION order** | §2b. Comment at :4186-4188 states the assumption in words |
| **L7** | `sn/loss_representation/__init__.py:4652` | transpose: `for global_n in reversed(ordinates_in_level)` | **YES (inherited)** | the exact reverse of L4 |
| **R1** | `sn/operators/radial_characteristic.py:780,788-795` | `RadialCharacteristicSeeding.apply`: `ordinates = level_indices[p]`; `within = arange(size)`; `out[:, ordinates, i] = -upstream_numer/V` | **positional** | scatter back through the same order |
| **R2** | `sn/operators/radial_characteristic.py:852-860` | `apply_transpose`: `numer_bar = tuple(-ob[:, level_idx, :]/V ...)` | **positional** | transpose of R1 |
| **R3** | `transport/radial_characteristic_field.py:289-304` | `legvander(mu_p, ords.size-1)`; `einsum("n,nl,ngx->lgx", w_p, legendre, q_p)` | **order-INVARIANT; CARDINALITY-sensitive** | a sum over `n`. But the Legendre DEGREE is `M_p - 1` — halving the level halves the moment order |
| **R4** | `transport/radial_characteristic_field.py:312-315` | `most_inward = ords[int(np.argmin(mu[ords]))]` | **order-INVARIANT by re-derivation** | the second, independent tie-break #326's body flags. `argmin` returns the FIRST minimum, so it is `ords[0]` today but does not depend on it |
| **D1** | `orpheus/derivations/discrete/sn/contamination.py:26-91, 135-149, 189-199` | a **second, independent** `_alpha_dome` + `_cell_edge_cosines` + `morel_montry_weights` over `quad.level_indices` | **YES** | a TWIN PATH inside the shipped package — §6 |
| **Doc** | `docs/theory/foundations/structured_geometry.rst:323,386`; `methods/sn/curvilinear_one_group.rst:626,854`; `methods/sn/curvilinear_numerics.rst:2232`; `methods/sn/angular_quadrature.rst:41`; `methods/sn/index.rst:539`; `foundations/boundary_conditions.rst:5149` | prose describing the `level_indices[p][m]` contract | — | doc blast radius of any change |

**Protocol surface:** `orpheus/geometry/reduced_operator.py:207-212` declares
`AngularMeasure.level_indices` as a `@property` on the structural Protocol —
this is the typed contract `cylindrical_streaming` requires
(`reduced_operator.py:763-770` raises if `level_structure is None`).

---

## 1b. Order-sensitive vs grouping-sensitive — the short answer

**Genuinely needs "ascending in η" (would be WRONG under an arbitrary
permutation):**

- A4 `eta_edge[m+1] = 0.5*(eta[m] + eta[m+1])` with `eta_edge[0] = -sinθ`,
  `eta_edge[M] = +sinθ` — **the only place that reads `[m]` and `[m+1]` as
  geometric neighbours on the angle axis.** Everything else inherits from it.
- A1 α recursion — the α at position `m` is `-Σ_{k<m} w_k η_k`, a *cumulative
  integral in ω*. Permutation-invariant at its ENDPOINTS (telescoping ⇒
  `α[M] = -Σ w η = 0` always) but every INTERIOR value moves.
- A6 R12a `τ_raw[0]` — reads position 0.
- L6 coupled-pole — needs all `η<0` members BEFORE all `η>0` members (§2b).
- A3's minus sign / `mu_start = -sinθ` — the march starts at the most-inward
  edge.

**Grouping-only (correct under ANY permutation of a level's members):**
A2, C3 (given coverage), C8 (positional but self-consistent), W1–W4, L1, L3,
R1, R2 — every "gather by `level_idx`, scatter by `level_idx`" pair.

**Order-invariant by explicit re-derivation** (these RE-SORT or RE-ARGMIN, so
they survive a reorder for free): C5 `_edge_seed_stencil`
(`order = np.argsort(mu)`), R4 `argmin(mu[ords])`.

---

## 2. The α / `eta_edge` monotonicity finding

### 2a. The exact expressions

`orpheus/sn/sweep/pole_angular_closure.py:588-611`
(`morel_montry_tau_raw_per_level`, CYLINDRICAL branch):

```python
for level_idx in quad.level_indices:
    eta = quad.mu_x[level_idx]
    M = len(level_idx)
    sin_theta = np.sqrt(1.0 - mu_z[level_idx[0]] ** 2)
    eta_edge = np.zeros(M + 1)
    eta_edge[0] = -sin_theta                              # <-- assumes eta[0] is the MINIMUM
    for m in range(M - 1):
        eta_edge[m + 1] = 0.5 * (eta[m] + eta[m + 1])     # <-- [m], [m+1] as NEIGHBOURS
    eta_edge[M] = sin_theta                               # <-- assumes eta[M-1] is the MAXIMUM
    tau = np.empty(M)
    for m in range(M):
        deta = eta_edge[m + 1] - eta_edge[m]
        tau[m] = (
            (eta[m] - eta_edge[m]) / deta if abs(deta) > 1e-15 else 0.5
        )
```

`orpheus/geometry/reduced_operator.py:778-786` (`cylindrical_streaming`) —
comment at :777 states the contract in words, *"Ordinates are ordered by
increasing η within each level."*:

```python
for level_idx in angular_measure.level_indices:
    eta = angular_measure.mu_x[level_idx]
    w   = angular_measure.weights[level_idx]
    alpha = np.zeros(M + 1)
    for m in range(M):
        alpha[m + 1] = alpha[m] - w[m] * eta[m]
```

The monotonicity assumed is **`eta` non-decreasing over the level**, so that
(i) `eta_edge` is non-decreasing, (ii) `deta ≥ 0`, and (iii) each `eta[m]` lies
inside `[eta_edge[m], eta_edge[m+1]]` ⇒ `τ_raw ∈ [0, 1]`.

### 2b. The SECOND, coarser monotonicity: the sign partition

`orpheus/sn/loss_representation/__init__.py:4186-4188` says it outright:

> *"Mirror partners share a level, and the M-M thread sweeps inward ordinates
> first, so the captured value is always data."*

Enforced at `:4222-4238` (`if mu_n < 0: psi_in = inflow_full[...]` else
`psi_in = pole_outflow[mirror[global_n]]`) and `:4303-4307` (the capture).
`pole_outflow` is a zero-initialised `(N, ng)` buffer, so if an outward
ordinate is visited before its `reflection_index("x")` partner it silently
reads **zero** — no exception. This is a *weaker* requirement than full
ascending-η: it needs only **all `η<0` before all `η>0`**, which ascending η
supplies.

### 2c. What actually happens under a NON-ascending order — MEASURED

Probe: `Quadrature.product(n_mu=2, n_phi=8)`, cylindrical, level 0 reversed.

```
BASE (ascending η)     τ_raw[0] = [0.0, 1.0, 0.0, 1.0, 3.46e-16, 1.0, 9.28e-16, 1.0]
                       τ       = [0.5, 1.0, 0.5, 1.0, 0.5,       1.0, 0.5,      1.0]

REVERSED (descending)  τ_raw[0] = [ 1.0790,  1.0, 3.85e-16, 1.0, 3.46e-16, 1.0, 9.28e-16, -0.0790]
                       τ       = [ 1.0,     1.0, 0.5,      1.0, 0.5,      1.0, 0.5,       0.5]
any NaN? False
```

**Verdict: silent wrong answer.** No NaN, no raise.

- `deta` goes negative at the ends, so `τ_raw` leaves `[0,1]`: `+1.079` and
  `-0.079`.
- The `abs(deta) > 1e-15` guard only catches an *exact* collapse, never a
  sign flip.
- `morel_montry_tau_per_level` (`:669-675`) then applies
  `max(0.5, min(1.0, τ_raw))` — the structural clamp **absorbs the damage**,
  turning `1.079 → 1.0` and `-0.079 → 0.5`. The out-of-range value never
  reaches a consumer that could notice it.
- A6 (R12a) reads `τ_raw[0]`. Base gives `0.0` (non-carrying); reversed gives
  `1.079` (still non-carrying, by luck — it is outside `(0,1)` on the other
  side). **An intermediate ordering that lands `τ_raw[0] ∈ (0,1)` would flip
  the mesh to CARRYING**, switching on the entire ψ½ state block, the
  `RadialCharacteristic*` spaces, and the System-B operator — a structural,
  not numerical, change.

For contrast, the already-built `LevelStructure.fiber(0)` (lexsort in
`(hemisphere, azimuth)`) gives:

```
fiber(0) = [0 1 2 3 4 5 6 7]   η = [0.8165, 0.5774, 0.0, -0.5774, -0.8165, -0.5774, -0.0, 0.5774]
τ_raw    = [1.0790, 0.2929, 0.5, 0.7071, 0.5, 0.2929, 0.5, 0.5469]
```

i.e. the fiber ordering is a *single monotone traverse of the circle*, but
starting at `φ=0` — so it is monotone in ω, NOT in η, and A4's η-midpoint edge
construction is the wrong edge rule for it. `τ_raw[0] = 1.079` is again
out-of-range; the `{0, 1}` trichotomy fingerprint disappears entirely.
**Wiring `fiber()` in behind the existing `eta_edge` code would be a silent
re-discretisation, not a re-ordering.**

---

## 3. The half-range change-list

Two proposals must be separated before scoping (they read as one and have
opposite blast radii):

- **(a1) fold the ALGORITHM** — march `M_p/2` positions, lift the mirror
  partners back into the full `(N, ng, nx)` buffer. `level_indices` still
  covers `[0, N)`; only the *march* is halved.
- **(a2) fold the STATE** — `level_indices` holds `M_p/2` members and the
  partition **no longer covers** `[0, N)`. This is the "2× fewer unknowns"
  headline.

Everything in **§3.1 is common**; **§3.2 fires only under (a2)**.

### 3.1 Common to both folds (the level's `M_p` and weights change)

| site | what changes |
|---|---|
| `rules_product.py:156-173` / `rules_sphere.py:303-313` | the producer: emit the folded index list + a doubled-weight vector (or a per-level weight side-channel). Note the DiscreteMeasure's weights must NOT change — it is a full-S² cubature consumed by 2-D Cartesian, `symmetry.py`, and the registry's `degree_of_exactness` |
| `reduced_operator.py:783` `alpha = np.zeros(M+1)` | derived from `M` — no code change, but α now runs over the half range and its interior values change |
| `reduced_operator.py:794` `ΔA[:,None]/w_level[None,:]` → `(nx, M)` | shape derived. **Invariant under uniform weight doubling** (α scales with `w`, `ΔA/w` as `1/w`) — measured previously |
| `pole_angular_closure.py:599` `eta_edge = np.zeros(M+1)` and `:603` `eta_edge[M] = sin_theta` | the half range's upper endpoint is `η = ±sinθ` only if the fold's endpoint is `ω = 0` or `π`. Hébert's midpoint half-range puts NO node on the endpoint ⇒ `τ_raw[0]` becomes strictly interior ⇒ **every cylindrical level becomes R12a-CARRYING** |
| `augmented_mesh.py:811-816` R12a predicate | see above — the fold's construction (endpoint-inclusive vs midpoint) decides the mesh's carrying set. This is the single highest-consequence flip in the list |
| `loss_representation:4130-4181` + `4612-4645` the `#280 2.5b` direct-seed fold | its precondition is *"non-carrying, product t=0 / dead level-symmetric"*. It raises `NotImplementedError` loudly if a non-carrying level has a live 2-point stencil. If every level becomes carrying, this whole block goes dead and route (a) takes over |
| `radial_characteristic_field.py:299` `legvander(mu_p, ords.size-1)` | moment order is `M_p - 1`; halving `M_p` halves the reconstructed angular order of the ψ½ source. Only fires on CARRYING levels — which is exactly what the midpoint fold creates |
| `radial_characteristic_space.py` `for_levels(levels, ...)` | takes level INDICES, not member counts — unaffected |
| SH moment / scattering functionals | the ONE genuine break: a `ξ`-ODD moment (`φ_1^ξ`) goes from `-1.3e-16` to `+2.94` under naive halve-and-double. Even-parity moments reproduce to 5e-16. Live consumer: cylindrical P1 anisotropic scattering. Under **(a1)** this break VANISHES (moments are taken on the lifted, exactly-even ψ) |

### 3.2 Fires ONLY under (a2) — the non-covering partition

These are the sites that assume `⋃_p level_indices[p] == range(N)`. **Three of
the four fail silently or with a confusing error, not with a clear one.**

| site | failure mode | evidence |
|---|---|---|
| `pole_angular_closure.py:262-271` `_gather_per_ordinate` | **`IndexError`**. It computes `N = sum(int(lvl.size) for lvl in level_indices)` — the *partition* size, not the ordinate count — then does `out[level] = values` with GLOBAL indices. Probe: `level_indices = ([0,2,4,6],[8,10,12,14])` → `IndexError: index 8 is out of bounds for axis 0 with size 8`. Its own docstring at :264-267 anticipates *"a future non-covering partition then yields a deterministic 0"* — but the `N` computation makes that unreachable | measured |
| `sweep/cache.py:230-254` `GeometryCoefficients` | **silent garbage**. `chain_idx`/`A_down`/`A_total`/`dA_w`/`V` are `np.empty((N, nx))` with `N = quad.N`; only rows for visited `global_n` are written. Unfolded partners keep uninitialised memory | read |
| `loss_representation:2901-2914` `_degenerate_positions` | **length mismatch → `IndexError` at a distant site**. `global_deg = where(abs(mu_x) < eps)` runs over ALL N ordinates; `deg_level`/`deg_within` are appended only on a hit. Consumers at `:3306-3317` and `:3677-3705` do `n_deg = global_deg.size` then index `deg_level[col_idx]` | read |
| `loss_representation:4192-4210` `_run` / `:4647-4652` `_run_transpose` | ordinates absent from any level are **never swept**: `angular_flux[global_n]` stays zero and `scalar_flux` is under-counted unless the doubled weights exactly compensate | read |
| `_dag_legs` `:2828-2853` | half the legs disappear; downstream `_loop_walk` just does less work — no error | read |
| `reflection_index("x")` (`directional.py:373-394`) | the x-mirror maps `η → -η`, `ξ → ξ`, i.e. `φ → π-φ`, which **is closed on `ω ∈ [0,π]`** ⇒ L6's coupled-pole map SURVIVES the fold. The `y`-mirror (`ξ → -ξ`) is exactly what gets quotiented, and the previously-flagged `(-η,+ξ)` vs `(-η,-ξ)` pole-map ambiguity DISSOLVES | measured previously |
| `Quadrature.octants` (`directional.py:482`) | `partition_by(_octant_sign_predicate)` over the full measure — 8 octants → 4 under a state fold |
| Krylov `n_dof` / `restart` (ERR-053 family) | the composite's `to_flat` sizes from `(N, ng, nx)`; a state fold halves `N` |
| `PoleAngularClosureBase.c_*_per_ordinate` `(N,)` caches | built by C3; same coverage assumption |

### 3.3 `level_symmetric` on a cylinder needs its own ruling

`P2` fibers over `|μ_z|` (`PolarInvariant.ABS_MU_Z`), so a level is
**4-to-1** over η (both `μ_z` signs AND both `ξ` signs). Both redundancies are
physical for an infinite 1-D cylinder, but there is no LS analogue of "the
`ω = π` endpoint is already a node", so a 4-to-1 fold makes **every** LS level
R12a-carrying. `select_quadrature("cylinder", ...)` rejects LS
(`SO2 ⊄ O_h`) but that selector has no production caller;
`cylindrical_streaming` accepts any level-structured quadrature, so LS-on-cyl
is reachable by direct `SNMesh(...)` and is pinned by ~6 test modules.

---

## 4. `Quadrature` adapter surface — can a quotient type hide behind it?

**`LevelStructure` is passed through WHOLE, not destructured — at the
`Quadrature` layer.**

`orpheus/numerics/quadrature/directional.py`:
- `:237` — `level_structure: LevelStructure | None = None` is a stored
  dataclass field.
- `:488-523` — `n_levels`, `level_indices`, `level_mu` are read-through
  `@property`s, **not `cached_property`**, each with a `None` fallback
  (`1`, `[arange(N)]`, `array([0.0])`).
- `:576-618` — the two factories (`.level_symmetric`, `.product`) pass the
  producer's `structure` straight in; `.gauss_legendre` and `.lebedev` pass
  `level_structure=None`.

So a new quotient type that satisfies the same three-attribute surface (plus
whatever the SN side reads) drops in **at the `Quadrature` layer with zero
changes there**.

**But `Quadrature` is not the real boundary.** The destructuring happens one
layer down, in three places that SNAPSHOT the list:

| destructuring site | what it snapshots |
|---|---|
| `sn/sweep/pole_angular_closure.py:904-906` (`MorelMontryAngularSweep.__init__`) | `self.level_indices = tuple(np.asarray(lvl) for lvl in quad.level_indices)` — and the base class **re-declares** `level_indices: tuple[np.ndarray, ...]` at `:175` as part of the closure's own contract |
| `sn/sweep/cache.py:246-254` (`GeometryCoefficients.from_mesh_and_quad`) | `level_ordinates` + a flattened `level_visits_iter` |
| `geometry/reduced_operator.py:779, 791, 807` (`cylindrical_streaming`) | three separate re-reads, one per derived array |

Downstream, most consumers read **`mesh.pole_angular_closure.level_indices`**
(the closure's copy), not `quad.level_indices`: `loss_representation:2827,
2902, 3518`, `radial_characteristic_field:289`,
`radial_characteristic.py:780, 852`. A minority still read the quadrature
directly: `cache.py:247`, `augmented_mesh.py:1261, 1317, 1468`,
`loss_representation:4089, 4155, 4609, 4621`.

**Consequence for the design question:** a quotient level type is
transparently insertable behind `Quadrature`, but the *ordering contract*
itself lives on `PoleAngularClosureBase.level_indices` — that is the type the
SN sweep actually programs against, and it is currently an untyped
`tuple[np.ndarray, ...]` with the contract carried only in prose.

---

## 5. Things that surprised me

1. **The clamp is what makes a mis-ordering invisible.** `τ_raw` going to
   `+1.079` / `-0.079` is a loud, checkable violation — and
   `morel_montry_tau_per_level:669-675` immediately squashes it into
   `[0.5, 1]`. The raw producer was split out for the R12a predicate, so a
   `0 ≤ τ_raw ≤ 1` assertion is available *for free* at exactly the right
   place and is not made.

2. **The R12a predicate reads position 0 of an order that has a tie there.**
   `augmented_mesh.py:815` is `eps < float(tau_level[0]) < 1.0-eps`. The
   whole ψ½/System-B machinery switches on a number produced from whichever
   of two ULP-separated ordinates `np.argsort` happened to place first. The
   trichotomy is documented as "structural" (`τ_raw = 0` on product,
   `1` on level-symmetric, interior on sphere) — measured here, the
   product level's `τ_raw` is `[0, 1, 0, 1, 3e-16, 1, 9e-16, 1]`, i.e. the
   *whole vector* alternates. It is the doubled-covering fingerprint, and
   position 0's `0` is the only entry that gets consulted.

3. **A twin path lives inside the shipped package.**
   `orpheus/derivations/discrete/sn/contamination.py:26-91` re-implements
   `_alpha_dome` and `_cell_edge_cosines` (cylindrical branch **line-for-line
   identical** to `pole_angular_closure.py:599-603`) plus a third
   `morel_montry_weights` τ producer at `:186-200`. Its docstrings even
   spell the contract (*"sorted in increasing order"*, *"correct for η-sorted
   ordinates with potential duplicate η values (from paired ±ξ ordinates)"* —
   the codebase already knew about the 2-to-1 degeneracy when this was
   written). Any change to the edge rule must move BOTH.

4. **Two consumers already defend themselves by re-deriving.**
   `_edge_seed_stencil` re-`argsort`s and `most_inward` re-`argmin`s. #326's
   body calls the second one "a second, independent tie-break on the same
   degeneracy" — true, but it also means those two sites are the only ones
   that would survive a reorder unchanged. They are accidental evidence of
   what a de-coupled design looks like.

5. **`fiber()` is not a drop-in for `level_indices`.** It is monotone in ω
   (the circle), not in η. The existing `eta_edge` midpoint rule is an
   η-space construction; feeding it a ω-monotone order changes the
   discretisation (τ_raw becomes `[1.079, 0.293, 0.5, 0.707, 0.5, 0.293,
   0.5, 0.547]` — a completely different angular cell partition), it does not
   merely re-index it.

6. **`_gather_per_ordinate` computes `N` from the partition, not the
   quadrature** — a one-line latent bug (`N = sum(lvl.size ...)` should be
   `quad.N`) that its own docstring's "future non-covering partition" comment
   shows was thought about but not closed.

7. **The `hemisphere`/`azimuth` arrays are per-ORDINATE (global), while
   `level_indices` is per-LEVEL.** Under an (a2) state fold the fiber
   coordinates keep full length while the index lists halve — the two halves
   of `LevelStructure` would carry different notions of "the ordinate set".

---

## Appendix — probe transcripts

```
# Quadrature.product(n_mu=2, n_phi=8): N=16, n_levels=2
level 0: idx=[4 5 3 6 2 7 1 0]
         eta=[-0.816497 -0.57735 -0.57735 -0. 0. 0.57735 0.57735 0.816497]
         xi =[ 0.       -0.57735  0.57735 -0.816497 0.816497 -0.57735 0.57735 0.]
tau_raw[0] = [0.0, 1.0, 0.0, 1.0, 3.46e-16, 1.0, 9.28e-16, 1.0]
tau[0]     = [0.5, 1.0, 0.5, 1.0, 0.5,      1.0, 0.5,      1.0]

# same level, order REVERSED
REVERSED tau_raw[0] = [1.0790, 1.0, 3.85e-16, 1.0, 3.46e-16, 1.0, 9.28e-16, -0.0790]
REVERSED tau[0]     = [1.0,    1.0, 0.5,      1.0, 0.5,      1.0, 0.5,       0.5]
any NaN? False

# same level, LevelStructure.fiber(0) ordering (lexsort hemisphere, azimuth)
fiber(0) = [0 1 2 3 4 5 6 7]
eta      = [0.816497 0.57735 0. -0.57735 -0.816497 -0.57735 -0. 0.57735]
tau_raw  = [1.0790, 0.2929, 0.5, 0.7071, 0.5, 0.2929, 0.5, 0.5469]

# _gather_per_ordinate under a non-covering partition
full cover  -> (16,)
half cover  -> IndexError: index 8 is out of bounds for axis 0 with size 8
```

Test modules touching `level_indices` / `LevelStructure` / `.fiber(` (23):
`tests/numerics/test_rules_product.py`, `test_quadrature_directional.py`,
`test_rules_sphere.py`; `tests/sn/_test_helpers.py`;
`tests/sn/mesh/test_radial_characteristic_slot_coordination.py`;
`tests/sn/operators/test_loss_transpose_solve.py`, `test_psi_half_coupling.py`;
`tests/sn/primitives/test_snmesh_consumes_reduced.py`, `test_quadrature.py`;
`tests/sn/verification/analytical/test_kinf_homogeneous.py`;
`tests/sn/verification/mms/test_mms_ordering_blindness.py`;
`tests/sn/sweep/test_cyl_direct_seed_fold.py`;
`tests/sn/sweep/core/{test_ordinate_scan,test_cell_visit_c_stamp,test_dag_walk}.py`;
`tests/sn/sweep/curvilinear/{test_unified_matvec_cylinder,test_azimuthal_mirror_symmetry,test_alpha_closed_form,test_pole_angular_closure,test_coupled_pole_mu_level_invariant,test_cyl_sweep_regression}.py`;
`tests/sn/sweep/slab/test_unified_matvec_slab.py`;
`tests/geometry/test_reduced_operator.py`.
