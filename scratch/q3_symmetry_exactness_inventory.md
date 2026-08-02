# Q3 — Symmetry-exactness inventory of node/direction/angle generators

**Scope:** production tree `orpheus/`, excluding `orpheus/derivations/**` (reported
separately in §7). Branch `refactor/operator-strategy-layers`, surveyed 2026-08-02.

**Question:** where does ORPHEUS GENERATE a set of nodes/directions/angles that is
*supposed* to carry a symmetry, and is that symmetry EXACT (an integer permutation of
node indices) or only approximate (to ~1e-16)?

**Inventory only — no fixes proposed.** Every entry reports:
(a) the generating expression, (b) the claimed symmetry group, (c) the generation
KIND — `TABULATED` / `FORMULA` (cos/sin of a grid) / `GROUP-ACTION` (a representative
acted on by exact operators), (d) the nearby tolerance constant used to test it.

Legend for (c):

| kind | meaning |
|---|---|
| **GROUP-ACTION** | representative + exact operators (sign flips, swaps, index arithmetic). Symmetry is bit-exact by construction. |
| **FORMULA** | `cos`/`sin` (or `sqrt`) evaluated independently per node on a parameter grid. Symmetry holds only to round-off. |
| **TABULATED** | literal float table. Symmetry is whatever the digits happen to encode. |
| **EXTERNAL** | delegated to numpy/scipy; kind is whatever that library does. |

---

## 1. `orpheus/numerics/quadrature/` — the SN angular rules

### 1.1 `rules_product.py` — GL(μ) × equispaced(φ) product rule on S²  ⚠ **PRIMARY FINDING**

**File:** `/Users/rodrigo/git/nuclear/ORPHEUS/orpheus/numerics/quadrature/rules_product.py`

(a) **Generation expression** — `rules_product.py:129-146`:

```python
phi_pts = np.linspace(0, 2 * np.pi, n_phi, endpoint=False)   # :129
w_phi   = 2.0 * np.pi / n_phi                                # :130
...
mu_gl, w_gl = np.polynomial.legendre.leggauss(n_mu)          # :125
...
sin_theta = np.sqrt(1.0 - mu_val**2)                         # :142
mu_x[idx] = sin_theta * np.cos(phi_pts[m])                   # :145
mu_y[idx] = sin_theta * np.sin(phi_pts[m])                   # :146
mu_z[idx] = mu_val                                           # :147
```

(b) **Claimed symmetry** — `invariance_group=SubgroupOfO3.Dnh(n_phi)` (`:168`),
i.e. **D_{n_φ h}, order 4·n_φ**. The module docstring (`:24-32`) states the claim
explicitly and — importantly — says it is *computed*, not asserted:

> "`maximal_invariance_groups` **computes** exactly that from the nodes, so the tag
> above is checked rather than asserted."

The docstring `:34-49` records that this rule carried a **false `SO(2)` tag until
2026-08-02**, and that the false tag was load-bearing (it was the only reason the
registry's geometry gate admitted this rule for a cylinder).

(c) **Kind: FORMULA.** The azimuth is `np.cos`/`np.sin` of a `linspace` — the exact
construction `roots_of_unity.py` was written to replace. The polar factor is
`EXTERNAL` (numpy `leggauss`, which *does* impose the ±μ mirror — see §1.5).
`sin_theta = sqrt(1 - μ²)` is FORMULA but is applied per-LEVEL (one value per μ), so
it does not itself break the azimuthal symmetry.

Consequences, measured previously and recorded in the module's own prose and in
`tests/sn/_test_helpers.py:994-1010`:
- `cos(π/2) = 6.12e-17`, `sin(π) = 1.22e-16` — the on-axis ordinates are not exactly
  on-axis, and the defect is **asymmetric between components** (`sin(0)` IS exactly 0).
- The x-mirror `k -> (n-k) mod n` is bit-exact for **no** `q` under `linspace+cos`
  (`roots_of_unity.py:107-118`).

(d) **Nearby tolerances:** none inside `rules_product.py` itself. The rule has no
`atol`/`eps` of its own — the tolerance lives downstream in
`symmetry.py`'s `is_invariant(..., atol=...)` (§4) and in
`directional.py`'s `argmin` partner search (§5.1), which has no tolerance at all.

**Additional symmetry-adjacent site in the same function — the per-level sort:**

`rules_product.py:151-154`:
```python
# Sort by increasing η (mu_x) for cylindrical sweep convention.
level_arr = np.array(level_idx)
order = np.argsort(mu_x[level_arr])
level_indices.append(level_arr[order])
```
This sort key `η = sinθ·cos φ` is **2-to-1 over φ** (the ξ-mirror pair shares η). Today
the two partners differ by ~1 ULP of `cos φ`, so the (unstable, `quicksort`) `argsort`
is accidentally total. Making the azimuth exact **manufactures exact ties** here — this
is the documented `#325`-blocked-by-`#326` coupling (test-architect memory L34; the
issue-326 memory in `.claude/agent-memory/explorer/`). Report as a symmetry-exactness
*consequence* site, not a generator.

**Doc/citation note (pre-existing, unrelated to exactness):** `rules_product.py:54-58`
and `:82` cite "Bailey, Adams, Yang & Zika (2009) … Eq. 50" for the η-ordering
convention; that citation is retracted in
`docs/theory/methods/sn/curvilinear_one_group.rst` ("Citation correction").

### 1.2 `rules_sphere.py::lebedev_sphere` — Lebedev O_h rule

**File:** `/Users/rodrigo/git/nuclear/ORPHEUS/orpheus/numerics/quadrature/rules_sphere.py`

(a) `rules_sphere.py:91-96`:
```python
from scipy.integrate import lebedev_rule
pts, w = lebedev_rule(order)
nodes = np.ascontiguousarray(pts.T)  # (N, 3)
```

(b) **Claimed symmetry:** `invariance_group=SubgroupOfO3.OctahedralOh` (`:101`).
Docstring `:64-67`: "O_h (full octahedral group, 48 elements). The construction is
O_h-invariant by design (Lebedev 1976) — every node sits in an O_h-orbit and the
corresponding weights are equal across the orbit."

(c) **Kind: EXTERNAL** (scipy `lebedev_rule`). Not inspected here; scipy's Lebedev
tables are generated from orbit representatives, so the O_h action is *likely* exact
up to the tabulated representative's own rounding, but ORPHEUS neither generates nor
verifies that at construction. **No ORPHEUS-side group action is applied.**

(d) **Nearby tolerances:** none in the function.

### 1.3 `rules_sphere.py::_build_level_symmetric_arrays` — Carlson-Lathrop LS_N

**File:** `/Users/rodrigo/git/nuclear/ORPHEUS/orpheus/numerics/quadrature/rules_sphere.py`

(a) **Generation** — `rules_sphere.py:165-206`. Level cosines by formula:
```python
mu1_sq = 1.0 / (sn_order * (sn_order + 2) / 4)            # :168
delta  = 2.0 * (1.0 - 3.0 * mu1_sq) / (sn_order - 2)      # :169
mu2_levels = mu1_sq + np.arange(n_half) * delta           # :170
mu_levels  = np.sqrt(mu2_levels)                          # :172
...
xi = np.sqrt(max(xi_sq, 0.0))                             # :184
```
then the octant is replicated by an **explicit sign product**, `:194-201`:
```python
for eta, xi, mu_z in octant_dirs:
    for s_eta in (-1, 1):
        for s_xi in (-1, 1):
            for s_mu in (-1, 1):
                all_eta.append(s_eta * eta)   # ...
```

(b) **Claimed symmetry:** `invariance_group=SubgroupOfO3.OctahedralOh` (`:277`).
Docstring `:232-234`: "O_h — invariant under all 48 rotation / reflection elements of
the octahedral group." Module docstring `:8-15` repeats it.

(c) **Kind: HYBRID — GROUP-ACTION for the sign subgroup, FORMULA for the orbit
representative.** The 8-fold sign replication `:195-201` IS an exact group action
(multiplication by ±1 is exact on floats), so the three coordinate reflections σ_x,
σ_y, σ_z are **bit-exact**. But the coordinate *permutation* half of O_h (the
`(η,ξ,μ) -> (ξ,μ,η)` 3-cycles and transpositions, 6 of the 48 cosets) is only
approximate: `η` is drawn from `mu_levels[k]` (an exact table entry) while `ξ` is
`sqrt(sinθ² − η²)` (`:184`) — a *different* arithmetic path to what should be the same
value. So LS_N's O_h claim is exact on the sign subgroup D_2h (order 8) and
round-off-approximate on the index-6 permutation part.

(d) **Nearby tolerances — TWO, both bare literals:**
- `rules_sphere.py:182`: `if xi_sq < -1e-14: continue` — the negative-radicand guard
  that decides whether a candidate direction is admitted to the orbit at all.
- `rules_sphere.py:212-213`: `tol = 1e-12` then
  `idx = np.where(np.abs(np.abs(mu_z_arr) - level_mu_vals[p]) < tol)[0]` — the level
  membership test. **This is a symmetry test executed with a tolerance**: "which
  ordinates lie on level p" is decided by a float comparison rather than by index
  arithmetic, even though the generation loop `:174-201` knows the answer exactly
  (level `p` contributed a known contiguous slice of `octant_dirs`).

Same `argsort(mu_x[idx])` per-level ordering as §1.1 at `rules_sphere.py:214`.

### 1.4 `rules_sphere.py` — the exactness tag

`level_symmetric_sn` records `degree_of_exactness=sn_order-1` (`:278`) with the
docstring hedging "Conservatively … consumers that need a tighter guarantee should
refer to Lewis & Miller Table 4-2" (`:236-244`). Prior measurement (explorer memory
`quadrature_landscape_durable.md`) found the ACTUAL degree is **3 for every N** and
this is an EQUAL-weight construction (`w_octant = 4π/(8·n_octant)`, `:188`), not
Carlson-Lathrop moment-fitted. Recorded here because the exactness tag and the
invariance tag are the two "claims about the node set" the module makes, and one of
them is known-false. Tracked as **#327**.

### 1.5 `rules_1d.py::gauss_legendre_on_mu` — GL on [-1,1]

**File:** `/Users/rodrigo/git/nuclear/ORPHEUS/orpheus/numerics/quadrature/rules_1d.py`

(a) `rules_1d.py:80`: `nodes, weights = np.polynomial.legendre.leggauss(n)`

(b) **Claimed symmetry:** `invariance_group=SubgroupOfO3.SO2` (`:94`). The docstring
`:51-55` is explicit that the ±μ mirror is *deliberately not* the tagged property:

> "Stronger symmetry would require nodes paired (μ_i, −μ_i) symmetrically; **this
> happens to hold for Gauss-Legendre but is not the maximal property we tag here**"

and the inline comment `:85-93` defends `SO2` as naming the group that was *integrated
out* to give the μ-marginal (the measure lives on the quotient S²/SO(2)=[-1,1]).

(c) **Kind: EXTERNAL — and this is the ONE site where the reference implementation
already imposes the symmetry.** `numpy.polynomial.legendre.leggauss` finishes with
`x = (x - x[::-1])/2` and `w = (w + w[::-1])/2`, so the ±μ mirror IS bit-exact here.
That bit-exactness is what makes `Quadrature.gauss_legendre`'s
`identity[::-1].copy()` partner map (§5.2) a legitimate FORMULA rather than a
coincidence.

(d) **Nearby tolerances:** none. Comment `:100-107` records that `leggauss`
Newton-refines and differs from the tree's generic Golub-Welsch body by **1-4 ULP** in
the nodes, and that every SN slab snapshot is pinned to `leggauss`'s bits.

### 1.6 `registry.py` — the geometry gate that CONSUMES the invariance tag

Read separately in §4.3; it is not a node generator.

---

## 2. `orpheus/numerics/roots_of_unity.py` — the exact generator, **NOT WIRED**

**File:** `/Users/rodrigo/git/nuclear/ORPHEUS/orpheus/numerics/roots_of_unity.py`

(a) **Generation** — `roots_of_unity.py:232-262`, an integer quarter-split + octant
fold + sign-flip table:
```python
quad, r = np.divmod(4 * p, q)                       # :232  exact integer split
lo = np.minimum(r, q - r)                           # :238  fold to first OCTANT
theta = np.pi * lo / (2.0 * q)                      # :239
cos_oct, sin_oct = np.cos(theta), np.sin(theta)     # :240  trig on octant ONLY
swap = 2 * r > q                                    # :243  undo fold = coord SWAP
diagonal = 2 * r == q                               # :251  integer, not a tolerance
cos_quad = np.where(diagonal, _SQRT_HALF, cos_quad) # :252
cos = np.select(quadrants, [cos_quad, -sin_quad, -cos_quad, sin_quad])   # :257
return cos + 0.0, sin + 0.0                         # :262  canonicalise signed zero
```

(b) **Claimed symmetry:** the dihedral action on the q-th roots of unity — C_q
rotation `p -> p+1`, x-mirror `p -> -p`, y-mirror `p -> q/2 - p`, and the 45°
diagonal `|cos| == |sin|`. Docstring `:107-118` tabulates measured bit-exactness:
x-mirror exact for **every** q; y-mirror exact for every **even** q; diagonal exact
for every `q % 8 == 0`. Against `linspace+cos`: **no** q for any of the three.

(c) **Kind: GROUP-ACTION.** The only true group-action node generator in the tree.

(d) **Tolerances: ZERO — by design.** Both fixed points are decided by *integer*
arithmetic: the axis by `r == 0` (falls out — `:81-84`), the 45° diagonal by
`2*r == q` (`:251`). `_SQRT_HALF` (`:153`) is `float(np.sqrt(0.5))`, deliberately NOT
`np.cos(np.pi/4)` because that differs from `np.sin(np.pi/4)` by 1 ULP.

### 2.1 CONFIRMED: zero production consumers

Verified two ways.

**Graph** (`mcp__nexus__context` on
`py:function:orpheus.numerics.roots_of_unity.roots_of_unity`): every incoming `calls`
edge originates in `tests/` — `tests/numerics/test_roots_of_unity.py` (the intrinsic-
property gate module) and `tests/sn/_test_helpers.py::_product_rule_with_ordering`.

**Text grep** (whole repo, `*.py`/`*.rst`/`*.md`): the only non-test, non-scratch,
non-plan hits are inside `roots_of_unity.py` itself. It is **not re-exported** —
absent from `orpheus/numerics/__init__.py`.

**Intended consumers**, per `scratch/issue325_verification_plan.md:779` (step 6,
"export `roots_of_unity`; repoint the three consumers"):

| # | intended consumer | current expression | status |
|---|---|---|---|
| 1 | `orpheus/numerics/quadrature/rules_product.py` (plan cites `:115`; live line is **`:129`**) | `np.linspace(0, 2π, n_phi, endpoint=False)` + `np.cos`/`np.sin` at `:145-146` | **NOT repointed** |
| 2 | `orpheus/moc/quadrature.py:89` | `np.linspace(0, π, n_azi, endpoint=False) + π/(2·n_azi)` | **NOT repointed** |
| 3 | `orpheus/derivations/discrete/sn/balance.py:367` | (oracle tree — §7) | **NOT repointed** |

The campaign plan `.claude/plans/quadrature_machinery_campaign.md:57` states the
same in its own words: *"Not exported, no consumer repointed — the repoint is Q5's
business."*

The repoint is **blocked**, not merely pending: test-architect memory records that
`[M]` exactness makes the azimuthal mirror partners share `η` bit-exactly, so
`rules_product.py`'s per-level `np.argsort` (§1.1) acquires real ties, and the
ordering alone moves a 2-G cylindrical fixed-source flux by **1.008 %** at `n_phi=32`
(vs a **1.06e-14** node-value drift). The sort-key ruling is a prerequisite commit.

**Also present but unwired:** the module is the *generating* face of a
generate/check pair whose *checking* face is `orpheus/numerics/symmetry.py`
(`roots_of_unity.py:11-14`, `:122-125`). The check face IS wired (§4); the generate
face is not.

---

## 3. `orpheus/numerics/generating_measure.py` + `measure.py` — the 1-D primitives

### 3.1 `GeneratingMeasure.gauss` — **the symmetry IS imposed here** ✅

**File:** `/Users/rodrigo/git/nuclear/ORPHEUS/orpheus/numerics/generating_measure.py`

(a) `generating_measure.py:302-334` — Golub-Welsch, then the imposition:
```python
nodes, eigenvectors = np.linalg.eigh(jacobi_matrix)      # :302
weights = beta[0] * eigenvectors[0, :] ** 2              # :303
if self.is_symmetric:                                    # :305
    nodes   = (nodes - nodes[::-1]) / 2.0                # :327
    weights = (weights + weights[::-1]) / 2.0            # :328
weights = weights * (beta[0] / weights.sum())            # :334
```

(b) **Claimed symmetry:** the reflection `x_i = -x_{n-1-i}`, `w_i = w_{n-1-i}` for an
EVEN weight. The predicate that gates it is `is_symmetric` (`:220-246`), a **computed
property** — "is the weight even, `w(-x) = w(x)`?" — with the docstring noting
(`:234-240`) that "scipy carries this as a hand-set `symmetrize` boolean".

(c) **Kind: EXTERNAL (`eigh`) + IMPOSED SYMMETRY.** This is the pattern the Q3 hunt is
looking for, already adopted. The comment `:321-326` states the transport rationale
verbatim:

> "an angular quadrature's invariance group is load-bearing (a reflecting boundary is
> exactly representable only if the node set is closed under the face reflection), and
> an EXACTLY symmetric rule makes that a matter of **integer index arithmetic rather
> than of tolerance**."

(d) **Tolerances:** none in the imposition. Measured defect `3.3e-16 (n=8)` /
`8.6e-16 (n=32)` → **exactly 0.0** (`:313-319`), with accuracy improving as a bonus.

**Reach:** `measure.gauss_legendre` (`measure.py:1191-1193`) and
`measure.gauss_chebyshev` (`:1237-1239`) route through `.gauss`. **BUT** — see §1.5 —
the SN production path `rules_1d.gauss_legendre_on_mu` does **NOT**: it calls
`np.polynomial.legendre.leggauss` directly (`rules_1d.py:80`) and only *tags*
`generating_measure=LEGENDRE` (`:108`). The comment `:100-107` says the divergence is
deliberate (1-4 ULP; snapshots pinned). numpy's `leggauss` performs the same
`(x - x[::-1])/2` imposition, so both routes are ±μ-exact — but by two different
code paths that are 1-4 ULP apart.

### 3.2 `measure.equispaced` — the midpoint circle rule

(a) `measure.py:1283`: `nodes = a + (np.arange(n) + 0.5) * h`

(b) **Claimed symmetry:** the docstring `:1257-1264` names it explicitly — midpoints
"integrate constants exactly **while preserving symmetry under reflection through the
centre of the interval**", in contrast to the left-endpoints the SN product rule uses.

(c) **Kind: FORMULA** on the *parameter* (angle) axis, not on the circle. The measure
never leaves parameter space, so no `cos`/`sin` is applied here — the symmetry break
would arrive at whatever maps `[a,b]` onto the circle.

(d) **Tolerances:** none.

**Production consumers: ZERO.** (Confirmed in explorer memory
`quadrature_landscape_durable.md` and re-grepped: `equispaced` and `gauss_chebyshev`
are exported from `orpheus/numerics/__init__.py` but have no production call site.)

---

## 4. Consumers of `orpheus/numerics/symmetry.py`'s invariance machinery

### 4.1 Complete consumer list

Grepped `is_invariant|singular_set|orbit_certificate|maximal_invariance_groups` across
`orpheus/` (excl. derivations) and `tests/`:

| consumer | site | tolerance passed |
|---|---|---|
| **`AngularSymmetry.admits_symmetry`** (`orpheus/numerics/quadrature/registry.py:648`) | `return self.discrete_residual.is_invariant(measure)` | **none — takes the default `atol=1e-13`** |
| `tests/numerics/test_symmetry.py` (~40 call sites) | `is_invariant`, `maximal_invariance_groups`, `singular_set`, `orbit_certificate` | default `atol` everywhere except `test_singular_set_membership_is_exact_not_thresholded:1019`, which **sweeps** `atol` to prove Σ does not move with it |
| `tests/numerics/test_registry.py:623-673` | `admits_symmetry` / `admits_domain` | default |

**That is the whole list.** `admits_symmetry` is the *only* production consumer of the
invariance machinery. `maximal_invariance_groups`, `orbit_certificate` and
`singular_set` have **zero production consumers** — they are test-only today, despite
`rules_product.py:30-32` citing `maximal_invariance_groups` as the thing that "computes
… so the tag above is checked rather than asserted" (the check runs in the test suite,
not at rule-construction time).

**And `admits_symmetry` itself is only reachable through `select_quadrature`, which has
no production caller** — grepped: every `select_quadrature` call site is in
`tests/numerics/test_registry.py`. Meshes are built by direct `Quadrature.<factory>` +
`SNMesh(...)` construction. So the computed-invariance gate is, in practice, an
opt-in test-suite instrument.

### 4.2 The tolerances inside the machinery

| constant / parameter | site | value | what it decides |
|---|---|---|---|
| `is_invariant(atol=...)` | `symmetry.py:422` | **`1e-13`** default | weight equality between matched orbit partners. Docstring `:469-472`: "matches the floating-point noise on a 1-D Gauss-Legendre weight computation" |
| `candidate_groups(atol=...)` | `symmetry.py:1327` | `1e-13` | off-axis cut `rho > atol` for the azimuth count |
| `maximal_invariance_groups(atol=...)` | `symmetry.py:1376` | `1e-13` | forwarded to `is_invariant` |
| `orbit_certificate(atol=...)` | `symmetry.py:1447` | `1e-13` | forwarded to `_orbit_closure` |
| `singular_set(atol=...)` | `symmetry.py:1470` | `1e-13` | forwarded; the Σ membership itself is tolerance-FREE |
| `_NODE_WINDOW_FACTOR` | `symmetry.py:1102` | **`100.0`** | node-match window `= atol * 100 = 1e-11`. Comment `:1094-1101` justifies the asymmetry vs the weight window |
| `_distinct_azimuths` | `symmetry.py:1318`, `:1321` | **`1e-9`** (bare literal, twice) | how far apart two azimuths must be to count as distinct — the coarsest tolerance in the module |
| `_is_reflection_invariant_1d` | `symmetry.py:800`, `:802` | `atol*10` (nodes) / `atol` (weights) | 1-D ±x closure |
| `_close_group` dedup | `symmetry.py:1064` | `atol=1e-10` (bare literal) | matrix-level dedup of generated group elements |
| `SubgroupOfO3.contains` | `symmetry.py:603` | `1e-9` (bare literal) | matrix-norm containment test; comment `:539` notes "the closest pair in C_n differs by more than the 1e-9 granularity" |

### 4.3 ⚠ The checking face is itself a SEARCH

`symmetry.py:1244-1268` (`_orbit_closure`) — the core of every invariance answer:

```python
moved = nodes @ M.T                                                   # :1245
dist = np.linalg.norm(moved[:, None, :] - nodes[None, :, :], axis=2)  # :1254
pi = np.argmin(dist, axis=1)                                          # :1255
if np.any(dist[index, pi] > atol * _NODE_WINDOW_FACTOR):              # :1256
    return None
if np.unique(pi).size != n:                                           # :1262  (ERR-073)
    return None
if np.any(np.abs(weights[pi] - weights) > atol):                      # :1264
    return None
```

The permutation π that certifies invariance is obtained by **O(N²) nearest-neighbour
`argmin` with a `1e-11` acceptance window** — not by index arithmetic. `OrbitCertificate`
is honest about this (`:1120-1124`): "the only place a tolerance enters is matching nodes
while BUILDING π, which is the one place the question is honestly numerical." If the
node sets were group-generated, that "one place" would collapse to an integer map.

Note `_orbit_closure` *does* verify bijectivity (`:1262`), which a plain nearest-neighbour
match does not — added under **ERR-073**.

---

## 5. Reflection-partner / mirror-index maps: SEARCH vs FORMULA

This is the section the brief calls out specifically. Five production sites compute a
partner map. **Three are SEARCH, two are FORMULA.**

| # | site | map | kind | tolerance |
|---|---|---|---|---|
| 5.1 | `numerics/quadrature/directional.py:125-147` `_find_reflections` | ordinate → reflected ordinate, all 3 axes | **SEARCH** (`np.argmin`, O(N²)) | **none — unconditional `argmin`, no distance check at all** |
| 5.2 | `numerics/quadrature/directional.py:543-548` `Quadrature.gauss_legendre` | slab GL, `i ↔ N-1-i` | **FORMULA** (`identity[::-1]`) | none needed |
| 5.3 | `moc/geometry.py:222-229` `_reflected_azi_index` | azimuthal index → its `π − φ` partner | **SEARCH** (`np.argmin`) | **none** |
| 5.4 | `moc/geometry.py:412-435` `MOCMesh._find_link` | track → reflected track by endpoint proximity | **SEARCH** (nearest point, Python loop) | **none — `best_dist` is never compared to a threshold** |
| 5.5 | `numerics/symmetry.py:1254-1255` `_orbit_closure` | node → group image | **SEARCH** (`argmin`) | `atol * 100 = 1e-11` (§4.3) |

### 5.1 `_find_reflections` — the SN specular table

`orpheus/numerics/quadrature/directional.py:144-147`:
```python
for i in range(n):
    dist = (rx - tx[i]) ** 2 + (ry - ty[i]) ** 2 + (rz - tz[i]) ** 2
    ref[i] = np.argmin(dist)
```
Called three times per sphere quadrature at construction
(`_compute_sphere_reflection_partners`, `:183-187`) for Lebedev, level-symmetric AND
product (`:568`, `:588`, `:614`).

The docstring `:139-141` names the alternative explicitly:

> "Conceptually the same operation as a `DiscreteMeasure.pushforward` followed by a
> **snapping** step; the explicit nearest-neighbour search lives here … because the SN
> consumer requires *integer* indices."

**There is no distance guard.** `argmin` always returns *something*; a node set NOT
closed under the reflection yields a silently-wrong table. The guard is downstream and
post-hoc — the three `_specular.py` invariants (§5.6).

This table is the load-bearing one: it is what every reflective/albedo/specular boundary
realizes through. Consumers (`reflection_index`): `orpheus/sn/boundary/realizer.py:383`,
`orpheus/geometry/boundary/_specular.py:145,184,220`,
`orpheus/sn/loss_representation/__init__.py:3294,3470,4189,4593`.

### 5.2 The slab GL partner — the one FORMULA in the SN path

`orpheus/numerics/quadrature/directional.py:543-548`:
```python
identity = np.arange(N)
partners = {
    0: identity[::-1].copy(),  # GL x-reflection: i ↔ N-1-i
    1: identity,               # 1-D: mu_y == 0
    2: identity,               # 1-D: mu_z == 0
}
```
Pure index arithmetic — exactly the "search becomes index arithmetic" endpoint. It is
*legitimate* only because numpy's `leggauss` imposes the ±μ mirror bit-exactly (§1.5,
§3.1). Docstring `:383-385` states the reason: "For slab GL1D the partners are derived
from the `i ↔ N-1-i` symmetry of the GL nodes."

### 5.3 / 5.4 The MoC partner maps

`orpheus/moc/geometry.py:222-229`:
```python
def _reflected_azi_index(phi: np.ndarray, azi_index: int) -> int:
    """Both vertical and horizontal reflections map phi -> pi - phi."""
    phi_refl = np.pi - phi[azi_index]
    if phi_refl < 0:   phi_refl += np.pi
    if phi_refl >= np.pi: phi_refl -= np.pi
    return int(np.argmin(np.abs(phi - phi_refl)))
```
The MoC azimuth (§6.1) is the MIDPOINT rule on `[0, π)`, whose `π − φ` partner is
exactly `n_azi − 1 − a` — pure index arithmetic. The code searches for it instead, and
the wrap-fix branches (`:225-228`) are float comparisons standing in for `mod n`
arithmetic. No tolerance; `argmin` cannot fail loudly.

`orpheus/moc/geometry.py:412-435` (`_find_link`) is a second, spatial, tolerance-free
nearest-neighbour search: it returns `best_idx` with no comparison of `best_dist`
against anything, so a track with no genuine reflected partner silently links to the
least-bad one.

### 5.5 See §4.3.

### 5.6 The three post-hoc guards on the SEARCH result

`orpheus/geometry/boundary/_specular.py` certifies the `_find_reflections` output
after the fact. These are checks, not generators, but they are where the tolerance for
"is the symmetry exact enough?" actually lives:

| invariant | site | test | tolerance |
|---|---|---|---|
| measure-preserving (ERR-042) | `_specular.py:145-163` | `np.allclose(m[perm], m, rtol=_MEASURE_RTOL, atol=0.0)` | **`_MEASURE_RTOL = 1e-12`** (`:100`) |
| involution (ERR-044) | `_specular.py:184-191` | `np.array_equal(ref[ref], np.arange(N))` | **EXACT — integer identity, no tolerance** |
| inflow→outflow (ERR-045) | `_specular.py:220-241` | `np.sign(partner_mu) != np.sign(mu_axis)`, gated by `np.abs(mu) > TANGENTIAL_EPS` | **`TANGENTIAL_EPS = 4·eps ≈ 8.9e-16`** (`orpheus/numerics/spaces/angular_trace_space.py:164`) |

`TANGENTIAL_EPS` is precisely the constant that an exact node set would retire: with
`np.linspace`+`cos`, the `φ = π/2` ordinate has `μ_x = 6.12e-17` (nonzero); with
`roots_of_unity` it is exactly `0.0` and the tangential class is decided by `== 0.0`.
`symmetry.py:1482-1486` already names this family:

> "The three ad-hoc float comparisons the tree grew for this question
> (`_OCTANT_SIGN_EPS`, `_MU_DIRECTION_EPS`, `TANGENTIAL_EPS`) were all asking it
> numerically; measured across 29 production rules, the separation between 'zero' and
> 'nonzero' cosines is a factor of 2.7e13, so the tolerance was never doing real work."

The three siblings, all live:

| constant | site | value |
|---|---|---|
| `TANGENTIAL_EPS` | `orpheus/numerics/spaces/angular_trace_space.py:164` | `4.0 * np.finfo(np.float64).eps` |
| `_OCTANT_SIGN_EPS` | `orpheus/numerics/quadrature/directional.py:95` | `1e-15` |
| `_MU_DIRECTION_EPS` | `orpheus/sn/loss_representation/__init__.py:2726` | `1e-15` |

`directional.py:91-94` explicitly asks that `_OCTANT_SIGN_EPS` be "kept in lockstep"
with `_DEGENERATE_ABS_MU_THRESHOLD` in `orpheus.sn.sweep` /
`orpheus.transport.spatial.diamond` — a lockstep maintained by comment, not by code
(the named threshold no longer exists under that name; the live sibling is
`_MU_DIRECTION_EPS`).

---

## 6. `orpheus/moc/`, `orpheus/cp/`, `orpheus/mc/`

### 6.1 `moc/quadrature.py` — MoC azimuthal + Tabuchi-Yamamoto polar

**File:** `/Users/rodrigo/git/nuclear/ORPHEUS/orpheus/moc/quadrature.py`

**Azimuthal:**

(a) `moc/quadrature.py:89`:
```python
phi = np.linspace(0, np.pi, n_azi, endpoint=False) + np.pi / (2 * n_azi)
omega_azi = np.full(n_azi, 1.0 / n_azi)     # :90
```
i.e. **midpoints** of an even partition of `[0, π)`, weights summing to **1**.

(b) **Claimed symmetry:** never named as a group. Two symmetry facts are stated in
prose only:
- `MOCQuadrature` docstring `:41-43`: "the supplementary angle (φ + π) is the same
  physical track traversed in the opposite direction" — i.e. the rule lives on the
  **quotient circle R/πZ**, an unnamed Z₂ quotient of the SN convention.
- `moc/geometry.py:223` / `:377`: "Both vertical and horizontal reflections map
  phi -> pi - phi" — the D₂-type reflection the reflective BC consumes.

The rule carries **no `DiscreteMeasure`, no `invariance_group`, no
`degree_of_exactness`, no registry entry** — `MOCQuadrature` is a plain frozen
dataclass entirely outside the tagged quadrature system.

(c) **Kind: FORMULA** on the angle; the direction cosines are evaluated LATER and
PER-TRACK-GENERATION, at `moc/geometry.py:319-320`:
```python
cos_phi = np.cos(phi_arr[a_idx])
sin_phi = np.sin(phi_arr[a_idx])
```
Note the midpoint offset means `φ = π/(2·n_azi)·(2a+1)` — a rational multiple of `2π`
with denominator `4·n_azi`, so it IS a root of unity and `roots_of_unity(2a+1, 4*n_azi)`
would generate it exactly. Its `π − φ` mirror is exactly the index `n_azi − 1 − a`.

(d) **Nearby tolerances:** `moc/geometry.py:115` and `:125` — `abs(cos_phi) > 1e-15` /
`abs(sin_phi) > 1e-15`, the axis-aligned-ray guards in `_ray_box_intersections`. These
are the MoC analogue of `TANGENTIAL_EPS`: with `n_azi` even the midpoint grid never
lands exactly on an axis, so they are latent; with a different offset they would fire on
values that "should" be exactly 0.

**Polar (Tabuchi-Yamamoto):**

(a) `moc/quadrature.py:20-33` — `_TY_TABLES`, a literal dict of `(sin_theta, weight)`
pairs for `n_polar ∈ {1,2,3}`, e.g. `np.array([0.166648, 0.537707, 0.932954])`.

(b) **Claimed symmetry:** none, correctly. The comment `:16-18` states "Weights sum to
0.5 (one hemisphere); full-sphere sum is 1.0" — the ±μ hemisphere mirror is implicit in
the half-space convention, never realized as a node-set closure.

(c) **Kind: TABULATED** (6 significant digits, `Yamamoto et al. 2007 Table 2`).

(d) **Tolerances:** none. Note the exactness claim this rule *does* make is not a
polynomial degree at all — it is a minimax fit on the Bickley `Ki₃` family, the one
in-tree exactness claim an `int degree_of_exactness` structurally cannot carry.

### 6.2 `orpheus/cp/` — **no directional node set**

`orpheus/cp/solver.py` performs the angular integration **analytically**: the `E₃` /
`Ki₃` / `exp(−τ)` kernels (`:94-95`, `:184`, `:191`) already contain
`∫dμ` / `∫dφ`. The only quadrature is **spatial** — the `y` (chord-offset) rule at
`cp/solver.py:170`:
```python
q = composite_gauss_legendre(self.mesh.edges.tolist(), p.n_quad_y)
```
imported from `orpheus.derivations.common.quadrature` (`:44`). No symmetry group is
claimed for it and none is needed: it is a composite rule over cell edges, whose
panel structure is intentionally asymmetric. **No finding.**

### 6.3 `orpheus/mc/` — direction sampling

(a) `orpheus/mc/solver.py:410-413`:
```python
theta = np.pi * rng.random()
phi   = 2.0 * np.pi * rng.random()
dir_x = np.sin(theta) * np.cos(phi)
dir_y = np.sin(theta) * np.sin(phi)
```

(b) **Claimed symmetry:** **none is stated anywhere** — no isotropy claim in the
function, in `_random_walk`'s docstring, or in the module docstring
(`mc/solver.py:1-20`, which describes Woodcock delta tracking + analog absorption and
notes it is a "Port of MATLAB `monteCarloPWR.m`"). There is no `dir_z`: the walk is 2-D
in `(x, y)` on a pin-cell, and `(dir_x, dir_y)` is **not a unit vector** (its norm is
`|sin θ|`), so the sampled `theta` also scales the in-plane step.

(c) **Kind: FORMULA** on pseudo-random draws. There is no node set and no group action,
so the exactness question does not arise in the same form — but recorded here because
`theta = π·U` is uniform in the **angle**, not in `cos θ`; the isotropic measure on
`S²` is `cos θ ~ U(−1,1)`. Whether the `sin θ` factor is intended as the 2-D projection
of a 3-D isotropic flight is not stated. **Flagged as an unclaimed/undocumented
convention, not as a symmetry-exactness defect** — a physics question outside this
inventory's scope.

(d) **Tolerances:** none.

---

## 7. `orpheus/derivations/**` — the oracle tree (reported separately)

Two interesting cases; the rest of the derivations quadrature is geometric/adaptive
(kink-aware panel subdivision at tangent angles), where no symmetry group is claimed.

### 7.1 `derivations/discrete/sn/balance.py:360-374` — the product-rule TWIN

```python
mu_gl, w_gl = np.polynomial.legendre.leggauss(N)          # :360
alpha_sum = np.sum(w_gl * mu_gl)
assert abs(alpha_sum) < 1e-14                             # :362
...
phi = np.linspace(0, 2 * np.pi, n_phi, endpoint=False)    # :367
sin_theta = np.sqrt(1 - mu_z**2)                          # :370
eta = sin_theta * np.cos(phi)                             # :371
alpha_sum = np.sum(w_phi * eta)
assert abs(alpha_sum) < 1e-14                             # :374
```

(b) **Claimed symmetry** — stated in the docstring `:345-350` verbatim: "For GL
quadrature: `Σ w_n μ_n = 0` (**antisymmetry**)" and "For cylindrical (per level):
`Σ w_m η_m = 0` (**symmetry of η = sinθ cos φ over equally-spaced φ**)".

(c) **Kind: FORMULA** — a byte-level twin of the production `rules_product.py:129-146`.
This is the **third** intended `roots_of_unity` consumer named in
`scratch/issue325_verification_plan.md:779`.

(d) **Tolerance:** `1e-14`, twice (`:362`, `:374`). This is the clearest instance of the
principle in the brief: a claimed symmetry (`Σ w·η = 0` is an antisymmetry statement
about the node set) verified with an `assert abs(...) < 1e-14` where an exactly
generated node set would make it `== 0.0`. With `roots_of_unity` the mirror pairs
`cos(φ_k) = −cos(φ_{k+n/2})` are bit-exact, so `Σ w·η` would cancel in exact float
arithmetic.

### 7.2 `derivations/common/quadrature.py:227-246` — `_leggauss`, the dual-precision route

```python
if dps <= 53:
    nodes, wts = np.polynomial.legendre.leggauss(n)   # :239 — numpy DOES impose the mirror
    return nodes.astype(float), wts.astype(float)
with mpmath.workdps(dps):
    mp_nodes, mp_wts = mpmath.gauss_quadrature(n, "legendre")   # :242
return (np.fromiter(...), np.fromiter(...))                    # :243-246
```

(b) **Claimed symmetry:** none stated. (c) **Kind: EXTERNAL, two different externals.**
The `dps <= 53` branch inherits numpy's imposed `(x - x[::-1])/2`; the `dps > 53` branch
does **not** — `mpmath.gauss_quadrature` is not symmetrized here and the result is cast
to float64 node-by-node. Whether the two mirror partners round to exactly mirrored
float64s is **not measured** by anything in the tree. (d) **Tolerances:** none.

Consumers of the high-precision branch: `gauss_legendre(..., dps=...)`,
`gauss_legendre_visibility_cone`, `composite_gauss_legendre`,
`observer_angular_quadrature`, `surface_centred_angular_quadrature`. The last two
subdivide at tangent angles `arcsin(r_k/r_obs)` and their "backward mirrors `π − ω_k`"
(`quadrature_recipes.py:203-204`) — the mirrored panel breakpoints are computed by
`π − ω`, a FORMULA that is exact only up to the round-off in `π`.

### 7.3 `derivations/continuous/peierls_nystrom/geometry.py` — the named RANGE

`angular_range` returns `(-1, 1)` for slab-polar and `(0, π)` otherwise. This is the
only place in the tree where the angular RANGE is a named polymorphic property rather
than a literal; `numerics/quadrature` has no equivalent. Not a symmetry-exactness site,
but it is the existing prior art for the quotient-vs-cover distinction that §1.1 and
§6.1 both run into.

---

## 8. Downstream sites where the (in)exactness has teeth

Not generators — but these are where an approximate symmetry currently costs a
tolerance, an `argsort` tie-break, or a search. Listed because they are the blast
radius of any change at a §1 site.

| site | expression | what it stands in for |
|---|---|---|
| `sn/sweep/pole_angular_closure.py:598` | `sin_theta = np.sqrt(1.0 - mu_z[level_idx[0]] ** 2)` | **twin** of `rules_product.py:142` — the same quantity recomputed by a second path at consumption time |
| `sn/sweep/pole_angular_closure.py:602` | `eta_edge[m+1] = 0.5 * (eta[m] + eta[m+1])` | the ξ-mirror pair's midpoint edge **collapses onto the node** — the {τ=1, τ→0} split of #326 |
| `sn/sweep/pole_angular_closure.py:608` | `(eta[m]-eta_edge[m])/deta if abs(deta) > 1e-15 else 0.5` | the tolerance that catches that collapse; exact nodes make `deta` exactly `0.0` |
| `sn/sweep/pole_angular_closure.py:584` | same shape on the sphere's `dmu` | ditto |
| `sn/sweep/pole_angular_closure.py:1019-1023` | `order = np.argsort(mu)` then `if abs(mu[cand]-mu[m0]) > 1e-14` | the **distinct-η search** for the edge-seed stencil — a tie-break threshold that goes live the moment the nodes become exact |
| `rules_product.py:153` / `rules_sphere.py:214` | `np.argsort(mu_x[...])` (unstable quicksort) | per-level ordering; §1.1 |
| `numerics/basis/spherical_harmonic_basis.py:421-436` | `phi = np.arctan2(sin_phi, cos_phi)` then `np.cos(m*phi)`, `np.sin(m*phi)` | a **round-trip through the angle**: any exactness in the node's `(μ_x, μ_y, μ_z)` is destroyed by `arctan2` → `cos`. Guarded by `on_axis = sin_theta < 1e-15` (`:423`) |
| `transport/radial_characteristic_field.py:314` | `most_inward = ords[int(np.argmin(mu[ords]))]` | "most inward ordinate per level" by SEARCH over μ; with an ordered exact level it is index `0` |
| `sn/loss_representation/__init__.py:3294, 3470, 4189, 4593` | `quad.reflection_index("x")` | the r=0 pole continuity `ψ(0,+μ)=ψ(0,−μ)` — four consumers of the §5.1 SEARCH-built table |
| `sn/boundary/realizer.py:383`, `geometry/boundary/_specular.py:145,184,220` | `quadrature.reflection_index(axis)` | every reflective / albedo-specular BC |

**The FORMULA endpoint already exists in the tree:**
`orpheus/numerics/operator.py:2243` `PermutationOperator` — `self.inverse_perm =
np.argsort(perm)`, `is_involution = np.array_equal(perm[perm], np.arange(n))`
(`:2247-2249`). Its docstring `:2274` states the property the whole Q3 question is
after: *"`argsort` of a permutation is exactly involutive in integer math, so
`(P⁻¹)⁻¹` reproduces `perm` EXACTLY."* Once a symmetry is an integer permutation,
every downstream operation is exact for free.

---

## 9. Summary table — every generator, by kind

| # | site | what it generates | kind | claimed group | nearby tolerance |
|---|---|---|---|---|---|
| 1.1 | `numerics/quadrature/rules_product.py:129,145-146` | S² product ordinates | **FORMULA** | `Dnh(n_phi)`, computed-and-checked | none in-file |
| 1.2 | `numerics/quadrature/rules_sphere.py:91-96` | Lebedev ordinates | EXTERNAL (scipy) | `OctahedralOh`, by construction | none |
| 1.3 | `numerics/quadrature/rules_sphere.py:165-206` | LS_N ordinates | **HYBRID** — GROUP-ACTION on signs (`:195-201`), FORMULA on the orbit rep (`:184`) | `OctahedralOh` | `1e-14` (`:182`), `1e-12` (`:212`) |
| 1.5 | `numerics/quadrature/rules_1d.py:80` | GL μ nodes | EXTERNAL — **numpy imposes the mirror** | `SO2` (deliberately not the ±μ mirror) | none |
| 2 | `numerics/roots_of_unity.py:232-262` | q-th roots of unity | **GROUP-ACTION** ✅ | dihedral `C_q` + mirrors + diagonal | **none — integer arithmetic** |
| 3.1 | `numerics/generating_measure.py:327-328` | any Gauss rule for an even weight | EXTERNAL + **IMPOSED SYMMETRY** ✅ | reflection, gated by computed `is_symmetric` | none |
| 3.2 | `numerics/measure.py:1283` | midpoint rule on `[a,b]` | FORMULA (parameter space) | centre reflection, in prose | none |
| 6.1a | `moc/quadrature.py:89` (+ `moc/geometry.py:319-320`) | MoC azimuths / track directions | **FORMULA** | unnamed `R/πZ` quotient + `φ→π−φ` | `1e-15` (`geometry.py:115,125`) |
| 6.1b | `moc/quadrature.py:20-33` | TY polar `sinθ`/weights | **TABULATED** | none (correctly) | none |
| 6.3 | `mc/solver.py:410-413` | MC flight directions | FORMULA on RNG | **none stated** | none |
| 7.1 | `derivations/discrete/sn/balance.py:367,371` | product-rule twin | **FORMULA** | "antisymmetry" / "symmetry over equally-spaced φ" | **`1e-14` ×2** |
| 7.2 | `derivations/common/quadrature.py:239` vs `:242` | GL nodes, two precisions | EXTERNAL ×2 — only the numpy branch imposes | none stated | none |

**Kind census (production tree):**
- **GROUP-ACTION:** 1 (`roots_of_unity`) — and it has **zero production consumers**.
- **IMPOSED-SYMMETRY:** 2 (`GeneratingMeasure.gauss`; numpy's `leggauss` inside
  `rules_1d`).
- **HYBRID:** 1 (`level_symmetric_sn` — exact on the sign subgroup, approximate on the
  coordinate-permutation part).
- **FORMULA:** 3 production (`rules_product`, `moc/quadrature`+`moc/geometry`,
  `mc/solver`) + 1 derivations twin.
- **TABULATED:** 1 (TY polar).
- **EXTERNAL, unverified:** 1 (scipy Lebedev).

**Partner-map census:** 3 SEARCH (`_find_reflections`, `_reflected_azi_index`,
`_find_link`) + 1 SEARCH in the checker (`_orbit_closure`) vs **1 FORMULA**
(`Quadrature.gauss_legendre`'s `identity[::-1]`). Only the FORMULA one sits on a node
set whose symmetry is imposed.

---

## 10. Gaps and caveats

- **`orpheus/numerics/quadrature/registry.py`** was read only for its consumer role
  (`admits_symmetry`, `GEOMETRY_ANGULAR_SYMMETRY`); it generates no nodes. Its Stage-1
  gate is the only production reader of the invariance machinery, and it is unreachable
  in production because `select_quadrature` has no production caller.
- **scipy's `lebedev_rule`** was NOT inspected. The `O_h` claim at `rules_sphere.py:101`
  is inherited from the library, unverified in-tree at construction time (only in
  `tests/numerics/test_symmetry.py`).
- **`mpmath.gauss_quadrature`'s** symmetry (derivations `dps > 53` branch) is
  unmeasured.
- **No production site re-verifies a declared `invariance_group` at construction.**
  Every `DiscreteMeasure` is *tagged*; the computed check runs only in the test suite.
  Two false tags have shipped (the `SO(2)` product tag, retired 2026-08-02
  `rules_product.py:34-49`; the `cartesian2d` `O_h` residual, retired same day
  `registry.py:664-666`), which is what `symmetry.py:1297-1300` is referring to when it
  says "such claims have already shipped false here twice".
- **Doc blast radius noticed in passing (NOT a symmetry finding):** 18 references to the
  module `orpheus.sn.quadrature` survive in `orpheus/` and `docs/`, but that module does
  not exist (the SN adapter hierarchy was retired into
  `orpheus/numerics/quadrature/directional.py::Quadrature`). Sites:
  `numerics/operator.py:2190`; `numerics/measure.py:1184,1259`;
  `numerics/quadrature/rules_sphere.py:21,87,145,238,266`;
  `numerics/quadrature/rules_product.py:8,77,116`;
  `numerics/quadrature/rules_1d.py:11,74`;
  `geometry/reduced_operator.py:168,741,742`;
  `docs/theory/foundations/spherical_harmonics.rst:579`;
  `docs/theory/foundations/discrete_measures.rst:836`. These are Python-domain roles, so
  they render as plain text and produce **no `-W` warning** without `-n`.

---

## 11. ADDENDUM — the GROUP OPERATORS are generated the same way

A class of site the brief's list does not name but that belongs to the same question:
`symmetry.py` generates the matrices `M` that the invariance check applies. **Half of
them are exact sign/permutation matrices; half are `cos`/`sin` of an angle grid.** So a
`Dnh`/`Cn` invariance answer compares a FORMULA-generated node set against a
FORMULA-generated operator — two independent round-off sources, which is exactly what
the `_NODE_WINDOW_FACTOR = 100` widening (`:1102`) is absorbing.

| operator family | site | expression | kind |
|---|---|---|---|
| `_reflections(axis)` | `symmetry.py:888-899` | `M = np.eye(3); M[i,i] = -1.0` | **GROUP-ACTION / exact** |
| `_inversion_op()` | `symmetry.py:902-904` | `-np.eye(3)` | **exact** |
| `_octahedral_ops()` (all 48) | `symmetry.py:1007-1015` | `M[i, p] = signs[i]` over `product([-1.,1.],3) × permutations(3)` | **GROUP-ACTION / exact** — pure ±1 signed permutation matrices |
| `_rotation_z(theta)` | `symmetry.py:878-885` | `c, s = np.cos(theta), np.sin(theta)` | **FORMULA** |
| `_cyclic_ops(n)` — the `C_n` elements | `symmetry.py:955` | `[_rotation_z(2.0 * np.pi * k / n) for k in range(n)]` | **FORMULA** — these ARE the `n`-th roots of unity; `roots_of_unity(k, n)` generates them bit-exactly |
| `_vertical_mirrors(n)` — the `σ_v` of `D_nh` | `symmetry.py:984-990` | `theta = k*np.pi/n + np.pi/2.0`; `n_vec = np.array([np.cos(theta), np.sin(theta), 0.0])`; `M = I − 2·n̂n̂ᵀ` | **FORMULA** — also a root of unity (`2π(2k+n)/(4n)`) |
| `_rotation_about_axis(axis, theta)` | `symmetry.py:1080-1087` | Rodrigues, `c, s = np.cos(theta), np.sin(theta)` | **FORMULA** |
| `_icosahedral_ops()` | `symmetry.py:1032-1072` | golden-ratio vertex table + BFS closure under `R5`, `R3` | **FORMULA + tabulated seed**, with a `np.allclose(..., atol=1e-10)` dedup at `:1064` and an `iteration_cap = 200` safety bound at `:1058` |

Consequence worth stating plainly: the two node families that ARE (partly)
group-generated — Lebedev and level-symmetric — are checked against the **exact**
`_octahedral_ops`, so that pairing is the tree's cleanest. The `Dnh(n_phi)` claim on
the product rule (§1.1) is the *only* one where **both** sides are `cos`/`sin`
evaluations, and it is also the newest tag (2026-08-02).

`_vertical_mirrors`' docstring `:965-982` carries a directly relevant war story: the
mirror-plane *convention* was wrong (normals placed at `kπ/n` instead of `kπ/n + π/2`),
`Dnh(n).is_invariant(product(4, n))` read `False` for odd `n = 3,5,7`, and — the key
sentence — **"Orthogonality, determinant, closure and group order are all preserved by
a rotated mirror set, so none of those checks can see this; only comparing the angles
against the convention can."** A generator-property audit cannot see a convention error.

The docstring also states the setting the whole system assumes (`:965-971`): "the
principal axis along z with a vertex on the **x-axis** … This is the setting every
azimuthal rule in the tree is built in (`np.linspace(0, 2*pi, n, endpoint=False)` puts
a node at `φ = 0`)." The MoC azimuth (§6.1) is the **midpoint** grid and therefore does
NOT put a node at `φ = 0` — it is in the half-step-rotated setting. MoC is outside the
`DiscreteMeasure`/`SubgroupOfO3` system entirely, so nothing currently notices.

---

*End of inventory. No fixes proposed; every claim above is `file:line`-anchored against
the working tree on `refactor/operator-strategy-layers` as of 2026-08-02.*
