# Audit — geometric transformations buried in the BC layer

**Branch** `refactor/operator-strategy-layers` · repo `/Users/rodrigo/git/nuclear/ORPHEUS`
· audit date 2026-08-03. Working tree at open AND close: production files **CLEAN**
(only `.claude/skills/vv-principles/*` modified, `scratch/` untracked) — no
concurrent carve, so every `file:line` below is at HEAD.

**Scope:** `orpheus/geometry/boundary/`, `orpheus/sn/operators/boundary.py`,
`orpheus/sn/boundary/`, `orpheus/numerics/spaces/angular_trace_space.py`,
`orpheus/numerics/face_layout.py`. (+ `numerics/quadrature/directional.py` and
`numerics/operator.py` where the BC layer's transformations are actually realized.)

Spelling classes: **(M)** explicit matrix · **(P)** index permutation ·
**(S)** sign flip / coordinate negation · **(A)** angle arithmetic ·
**(O)** other.

---

## 0. The reference machinery (pre-read findings)

`orpheus/numerics/symmetry.py`:

- `Mirror(axis)` — `symmetry.py:179-241`. Parameterised by the mirror's NORMAL,
  restricted to `("x","y","z")` (`__post_init__`, `symmetry.py:234-241`).
  **Only the three COORDINATE mirrors are expressible.** Arbitrary-plane: NO.
- `_reflections(axis)` — `symmetry.py:1027-1038`: `I` with `M[i,i] = -1.0`.
  **(M)**, bit-exact.
- `_cyclic_ops(n)` — `symmetry.py:1092-1094` → `_rotation_z` `symmetry.py:1017-1024`:
  `cos`/`sin`. **NOT bit-exact.**
- `_vertical_mirrors(n)` — `symmetry.py:1097-1129`: Householder `I − 2 n nᵀ` from
  `cos/sin`. **NOT bit-exact.**
- `_realized_ops(tag)` — `symmetry.py:576-611`; `None` for continuous groups.
- `_orbit_closure(nodes, weights, ops, atol)` — `symmetry.py:1346-1407`: `(N,N)`
  distance matrix → `argmin` → **bijectivity** (`np.unique(pi).size != n`, ERR-073)
  → weight match. Returns `OrbitCertificate.permutations` — *the node permutation
  induced by each group element*.
- `OrbitCertificate` — `symmetry.py:1244-1343`: `singular_set` (integer identity
  `perm == idx`, exact), `stabilizer_order`, `orbits()`.

**The load-bearing upstream fact** (landed 2026-08-02, ~1 day old):
`directional.py:207-216` (`_compute_sphere_reflection_partners`) **already routes
through `_orbit_closure`**. It builds `reflection = np.eye(3); reflection[axis,axis]
= -1.0` **(M)** and asks the symmetry machinery for the permutation. So for every
**3-D** cubature the specular partner map IS `symmetry.py`'s group action.
**MEASURED equal**: `orbit_certificate(q.measure, SubgroupOfO3.Mirror(a)).permutations`
== `q.reflection_index(axis)`, bit-identically, on `product(2,4)`, `product(4,5)`,
`lebedev(17)` (probe §E). The **1-D Gauss-Legendre** path is NOT
(`directional.py:570-577`): hardcoded `identity[::-1]` **(S/P)**.

---

## 1. Findings table

| # | file:line | the transformation, in math | spelled TODAY | is it a group element? | EXACT? |
|---|---|---|---|---|---|
| **F1** | `orpheus/numerics/quadrature/directional.py:207-216` | `σ_a : Ω ↦ Ω − 2(Ω·ê_a)ê_a`, a ∈ {x,y,z}; the induced node permutation | **(M)** 3×3 `np.eye(3)` with `[axis,axis] = -1.0`, handed to `_orbit_closure` → `certificate.permutations[0]` **(P)** | **YES — `Mirror(axis)`**, literally `symmetry._reflections(axis)`. **ALREADY UNIFIED**: this IS the symmetry machinery, measured bit-identical. | Matrix exact. Permutation found by `argmin` over an (N,N) distance matrix at `atol=1e-13` ⟹ **tolerance-matched in MECHANISM**. **Measured exact in FACT**: match residual is `0.0e+00` for every node of `product(2,4)`, `level_symmetric(6)`, `lebedev(17)` on all 3 axes, with runner-up margin 0.26–1.16 (§B). |
| **F2** | `orpheus/numerics/quadrature/directional.py:570-577` | `σ_x` on a 1-D polar marginal: `μ ↦ −μ`, embedded as `(μ,0,0)` | **(S/P)** `identity[::-1].copy()` — hardcoded order reversal; axes 1,2 get bare `identity` | **YES — `Mirror("x")`** (`symmetry.py:876-891` derives exactly this from the `(μ,0,0)` embedding). **NOT routed through the machinery** — the one specular map in the tree that bypasses `_orbit_closure`. The `1: identity, 2: identity` entries ARE `Mirror("y")`/`Mirror("z")` fixing every node POINTWISE — the same fact `_check_invariance_1d` derives at `symmetry.py:888-891`, transcribed. | **BIT-EXACT** — pure integer index arithmetic. Measured: `nodes[perm] == -nodes` **exactly** (max abs of `μ[p]+μ` = `0.0`). The only structurally-exact partner map. |
| **F3** | `orpheus/geometry/boundary/_factors.py:474-520` (`SpecularMirror`) | `G` = mirror about `axis`: `Ω ↦ Ω − 2(Ω·n̂)n̂` | **(O)** pure DESCRIPTOR: one field `axis: str`, no matrix, no `as_operator`. Realizes (in the realizer) to `quadrature.reflection_index(axis)`. | **YES — `Mirror(axis)`**, same `{x,y,z}` parameter set, same math, **different type**. Two independent spellings of one concept with **no shared type and no cross-reference**: `_factors.py` never imports `symmetry`. | n/a (spec). Its realization inherits F1/F2. |
| **F4** | `orpheus/geometry/boundary/_factors.py:523-627` (`SpatialWrap`) | `G` = translation `x ↦ x + L` along `axis` — the deck transformation of `ℝ^d/Λ` | **(O)** descriptor `axis: str`; the realized action is a FACE-NAME map `face_opposite(face)` at `_factors.py:627` **(O)** — a translation VECTOR is never built | **YES in principle — a TRANSLATION.** `symmetry.py` **CANNOT express it**: `SubgroupOfO3` is a subgroup of `O(3)` (linear, origin-fixing); a translation is not in `O(3)`. Would need the Euclidean group `E(3)` / a space group. **Explicitly out of reach today.** | **EXACT** — string parse + integer sign negation, no float. |
| **F5** | `orpheus/numerics/face_layout.py:209-235` (`face_opposite`) | the spatial half of F4, realized: `n̂_f ↦ −n̂_f` | **(S)** `face_name(axis, -sign)` — a literal sign negation on the outward-normal sign, via the `_FACE_SUFFIX_SIGN` bijection (`face_layout.py:130-136`) | The **antipodal map on the 6-element face index set**. Not in `O(3)`'s node action, so not expressible in `symmetry.py`. It IS a `ℤ₂` action — on faces, not on ordinates. | **EXACT** — integer sign flip + dict lookup. Docstring calls it an involution (`face_layout.py:214`) — **MEASURED TRUE over all 6 faces** (§D), **but never asserted in production code**. |
| **F6** | `orpheus/sn/operators/boundary.py:361-382` (`SNBoundaryOperator._face_domains`) | the induced permutation of the FACE set: `f ↦ G_f.domain_face(f)` | **(P)** a `dict[str,str]` built by comprehension, then certified by `sorted(domains.values()) != sorted(faces)` | This is **exactly `_orbit_closure`'s job, one level up and on a different set**: "find the permutation a transformation induces on this finite set, and refuse if it is not a bijection." Faces are not `O(3)` nodes, so the machinery does not literally apply — but the *shape* is identical (permutation-or-refuse). | **EXACT** — string sorting/comparison. |
| **F7** | `orpheus/sn/boundary/realizer.py:382-420` (`_specular_kernel`) | the mirror NARROWED to `Γ₊ → Γ₋`: `j ↦ position of π(inflow[j]) inside Γ₊` | **(P)** `gamma_out.to_local(perm[inflow])` → `np.searchsorted` (`operator.py:2451-2477`), then `PermutationOperator(local_perm, axis=0)` | Still `Mirror(axis)`, **restricted to a half-trace**. NOTE the restriction is a bijection between **two different index sets**, so it is a *coset* representative, not a group element of a group acting on Γ₊. | **EXACT** — `searchsorted` on integer arrays + a membership guard. |
| **F8** | `orpheus/sn/operators/boundary.py:928-933` (`RadialCharacteristicBoundaryOperator._reflect_corner`) | the SAME specular mirror at the **off-quadrature** `μ = ±1` ray corner: `corner(level, −1) ← corner(level, +1)` | **(S)** a literal `±1` sign key swap in a Python loop — `out.corner(level, -1)[...] = seed.corner(level, +1)` | **YES — the same `Mirror(axis)`**, applied to a direction that is **not in the quadrature**. Gated by `law_permutes_ordinates(law)` (`_base.py:71-130`) rather than by any matrix. This is the mirror spelled a **third** way, in a third vocabulary. | **EXACT** — an integer key swap, no arithmetic at all. |
| **F9** | `orpheus/numerics/spaces/angular_trace_space.py:229-251` (`build_omega_dot_n`) | the face's outward normal `n̂_f = sign·ê_axis` and the signed projection `Ω·n̂_f` | **(S)** `omega_dot_n[f_idx] = sign * mu_axis` — the normal is stored SPARSELY as `(axis_index, sign)` (from `face_normal`), never as a 3-vector | The `±ê_a` normal set IS the `O_h` axis-vector orbit; the face-side flip `sign = ±1` is the `Mirror(axis)` action on that normal. Expressible only indirectly. | `sign * mu_axis` is **EXACT** (multiplication by ±1.0 is exact in IEEE-754). The *classification* built on it (`< −TANGENTIAL_EPS` / `> +TANGENTIAL_EPS`, lines 431/440) is tolerance-based, but measured bit-identical to a strict compare on every rule but `product` (§I: `product(2,4)` has 4 tangential of 8; `lebedev(17)` 12 of 110). |
| **F10** | `orpheus/geometry/boundary/_specular.py:220-241` (`assert_specular_pairing_maps_inflow_to_outflow`) | the mirror's SIGN-CLASS action: `sign(μ_{a,π(n)}) = −sign(μ_{a,n})` | **(S)** `np.sign(partner_mu) != np.sign(mu_axis)` — a sign comparison, plus a `TANGENTIAL_EPS` exemption | Not a group element itself — it is the **certification** that F1/F2's permutation realizes `Mirror(axis)`'s hemisphere exchange. `_orbit_closure` does NOT check this (it checks position + weight + bijectivity, not the sign class). | Tolerance-based via `TANGENTIAL_EPS` and `np.sign`. Fires at realization for reflective (`reflective.py:184-186`) and for albedo-with-specular-closure (`albedo.py:243-246`). |
| **F11** | `orpheus/geometry/boundary/_specular.py:114-163` (`assert_specular_pairing_measure_preserving`) | the mirror preserves the trace measure `m_n = w_n|μ_{a,n}|` | **(P)+(O)** `cosine_measure[perm]` vs `cosine_measure`, `np.allclose(rtol=1e-12)` | The measure-preservation half of "`Mirror(axis)` is an isometry". **`_orbit_closure` already checks the WEIGHT half** (`symmetry.py:1403`, `|w[pi] − w| > atol`) — this adds the `\|μ_a\|` factor, which `_orbit_closure` gets for free from the position match. **Partial duplication** (see Q2). | `rtol=1e-12`, `atol=0.0`. |
| **F12** | `orpheus/geometry/boundary/_specular.py:184-191` (`assert_specular_pairing_involutive`) | `π ∘ π = id` — that `σ_a² = e` | **(P)** `np.array_equal(ref[ref], np.arange(N))` | The order-2 statement for `Mirror(axis)`. `_orbit_closure` **does not** prove it (a bijection need not be involutive) — but it IS implied by the matrix `M² = I` **plus an exact position match**. So on today's data it is a theorem, and it is checked anyway. **MEASURED True on every certified axis of every production rule** (§A). | **EXACT** — integer array comparison. |
| **F13** | `orpheus/numerics/operator.py:2243-2249` (`PermutationOperator`) | `π⁻¹` via `np.argsort(perm)`; `is_involution` via `perm[perm] == arange` | **(P)** | The inverse group element `σ_a⁻¹` (= `σ_a` for a mirror). | **EXACT** (integer argsort). ⚠ See "Surprise 3" — the `is_involution` flag is **False** for the narrowed lebedev kernel, and the docstring at `operator.py:2188-2191` claims otherwise. |
| **F14** | `orpheus/geometry/boundary/_factors.py:672-721` (`LambertianReemission`) + `orpheus/sn/boundary/angular.py:308-312` | `R = α·(1/norm)·1_{Γ₋} ⊗ (w·\|Ω·n̂\|)\|_{Γ₊}` — a rank-one cosine-weighted average | **(O)** `cos_w = weights[outflow] * omega_dot_n[outflow]`; apply = contract + `np.broadcast_to` | **NO geometric transformation hides here — CONFIRMED.** The `axis`/`outward_sign` fields select a *hemisphere* via `face_name(...)` + `build_omega_dot_n`; that is a FACE-ORIENTATION lookup (F9), not a transformation of directions. The Lambertian is not multiplicative and not a bijection — the `_factors.py:34-61` membership test rejects it from `G`. Its `G` is `IdentityMap` (`white.py:114-132`, `albedo.py:142-158`). | The hemisphere selection is `TANGENTIAL_EPS`-classified (F9). |
| **F15** | `orpheus/geometry/boundary/_factors.py:635-669` (`ScalarResponse`) | `R = α·I` | **(O)** a bare float. | **NO transformation.** Vacuum/prescribed/zero-flux/bare-albedo all carry `IdentityMap` (`vacuum.py:64`, `prescribed_inflow.py:129`, `zero_flux.py:112`, `albedo.py:143`). | n/a. |
| **F16** | `orpheus/geometry/boundary/_factors.py:724-795` (`SpecularReemission`) | `R = α·P_mirror` — the SAME mirror as F3, in the RESPONSE tier | **(O)** descriptor `alpha, axis`; realizes through the same `reflection_index(axis)` (`realizer.py:767-776` → `_specular_kernel`) | **YES — `Mirror(axis)` again.** Deliberately a different TYPE from `SpecularMirror` (documented at `_factors.py:733-762` as the point, not a smell). So `Mirror(axis)` has **four** spellings in the BC layer: `SpecularMirror`, `SpecularReturn`, `SpecularReemission`, and `symmetry.Mirror`. | inherits F1/F2. |
| **F17** | `orpheus/transport/method.py:315-318, 341-344, 353-...` (`_law_from_tag`) | which mirror plane / which wrap axis a declared BC means | **(S/O)** `AXIS_NAMES[label.axis_index]` for reflective; `face_normal(label.face_name) → (axis_index, outward_sign)` for white | The **choice of group element from the face**. This is where the `Mirror(axis)` parameter is actually bound, and it is derived from the face, not hardcoded — correctly. | **EXACT** (string parse). |

**Nothing else in the BC layer performs a geometric transformation.**
`zero_flux.py`, `prescribed_inflow.py`, `_source.py`, `_composition.py`,
`_bound_compat.py`, `_realizer.py`, `_errors.py` and `FaceLayout`/`FaceSlot`
(`face_layout.py:248-400`) carry none — checked by full read / targeted grep for
`[::-1] argsort np.sign outward_sign 2*np.pi np.pi- np.eye @M .T flip reverse
mirror opposite partner conjugate permut`.

---

## 2. The six specific questions

### Q1 — How does the specular BC build its partner map? Is the involution/permutation guarded or asserted?

**By MATRIX ACTION + nearest-neighbour search, with a bijectivity guard, for
3-D; by INDEX ARITHMETIC for 1-D. Both, then, get three MORE asserted invariants
at realization.** Precisely:

1. **Construction (3-D)** — `directional.py:207-216`: builds the explicit
   reflection matrix, calls `_orbit_closure`. Inside (`symmetry.py:1383-1404`):
   `moved = nodes @ M.T`; full `(N,N)` distance matrix; `pi = np.argmin(dist, axis=1)`;
   then **three guards, all inside `_orbit_closure`** —
   (a) `dist[index, pi] > atol * 100` ⟹ reject (some image is not a node),
   (b) **`np.unique(pi).size != n` ⟹ reject — the BIJECTIVITY guard (ERR-073)**,
   (c) `|w[pi] − w| > atol` ⟹ reject.
   A rejected axis is **OMITTED from the dict**, so `reflection_index(axis)` raises
   rather than returning a wrong map (`directional.py:416-424`).
2. **Construction (1-D)** — `directional.py:570-577`: `identity[::-1]`, **no guard
   at all**, correct by GL-node symmetry.
3. **Certification (both)** — `_specular.py:244-255` `assert_specular_pairing_valid`
   fires **three independent asserted invariants**: measure-preservation (F11),
   **involution (F12, asserted — `np.array_equal(ref[ref], arange)`)**, and
   inflow→outflow (F10). Wired from `reflective.py:210-213` and `albedo.py:240-246`,
   both reached from `SNBoundaryRealizer.realize` at `realizer.py:636-639`.

So: **the involution is ASSERTED, not guaranteed by construction.** The
bijectivity IS guaranteed (guard (b)). The module docstring at `_specular.py:41-57`
tabulates why all three are needed.

**Measured**: on every production rule × every certified axis, involution holds and
the match residual is exactly `0.0` (§A, §B).

**Live consequence found**: `product(4, 5)` (odd `n_φ`) is **not** `σ_x`-closed, so
axis 0 is omitted, and `ReflectiveBoundary(axis="x").assert_realizable(...)` raises
a bare **`ValueError`** — *not* a `BoundaryError` — because `reflection_index` raises
before any BC-layer code runs (`directional.py:420`; measured §F). The refusal is
correct; the error type escapes the BC layer's own taxonomy.

### Q2 — Does anything in the BC layer duplicate `_orbit_closure`'s job?

**Partially — and the biggest duplicate was already retired one day ago.**

- **NO** hand-written "find the permutation induced by a transformation" loop
  survives in `orpheus/geometry/boundary/` or `orpheus/sn/boundary/`. Grep for
  `for … in range / enumerate / np.argmin / np.isin / np.searchsorted /
  np.flatnonzero` over the whole BC layer returns 6 hits, all inflow/outflow
  `flatnonzero` classification (`realizer.py:484,545`, `angular.py:291-292`,
  `_specular.py:231`) — none of them build a partner map.
- **The predecessor DID duplicate it.** `directional.py:168-190` records that until
  2026-08-02 this ran "a bare `argmin` per node — no distance threshold, no
  injectivity check, no weight comparison", i.e. a hand-rolled, *unguarded*
  `_orbit_closure`. That is now gone.
- **F11 partially duplicates guard (c).** `_orbit_closure` already asserts
  `w[π] = w` (`symmetry.py:1403`); `assert_specular_pairing_measure_preserving`
  asserts `w[π]·|μ_a[π]| = w·|μ_a|`. Given an exact position match the second
  factor is implied by the first, so on any map that came out of `_orbit_closure`
  F11 is redundant. It is NOT redundant for the **1-D** path (F2), which never
  touches `_orbit_closure`.
- **F6 (`_face_domains`) is `_orbit_closure`'s shape on a different set** — "compute
  the induced permutation of a finite set and refuse if not a bijection"
  (`sn/operators/boundary.py:361-381`). Different set (faces, not `S²` nodes), so
  not literally the same machinery, but the same certificate concept.

### Q3 — Periodic: is the translation ever represented as a transformation?

**No. It is a face-partner lookup end to end.** There is no translation vector, no
lattice `Λ`, no domain length anywhere in the periodic path:

- `PeriodicBoundary` carries one field, `axis: str` (`periodic.py:88`).
- `SpatialWrap.domain_face(face)` (`_factors.py:577-627`) parses the face name,
  checks the axis matches, and returns `face_opposite(face)` — a **string**.
- `face_opposite` (`face_layout.py:209-235`) is `face_name(axis, -sign)` — an
  integer sign flip.
- The realized operator is the bare `IdentityOperator() & IdentityOperator()`
  (`realizer.py:856`); the whole content of the wrap lives in **which face's `Γ₊`
  the composition feeds it** (`sn/operators/boundary.py:523-525`).
- The identification `Γ₊(f') ≡ Γ₋(f)` is a **guard**, not a construction:
  `_assert_wrap_identification` (`realizer.py:423-511`) compares the two index sets
  with `np.array_equal` and raises otherwise.

`SpatialWrap.permutes_ordinates` is `False` (`_factors.py:560-563`) — correct, the
translation acts on the *spatial* coordinate only, so every angular discretization
is trivially equivariant under it (`docs/theory/foundations/boundary_conditions.rst:665-673`).

⚠ **The consequence for `symmetry.py`**: the periodic transformation is the ONE
in the BC layer that the point-symmetry machinery structurally cannot host. A
translation is not an element of `O(3)`. This is worth naming explicitly because
the `_factors.py:83` sentence "Mirrors, **translations** and rotations are `G`"
reads as if all three were the same kind of object; two of the three are in
`O(3)`, one is not.

### Q4 — Albedo / white: is a geometric transformation hiding inside them?

**No — CONFIRMED from the code, on both laws, both tiers.**

- `WhiteBoundary.geometry_map → IdentityMap()` (`white.py:114-132`);
  `AlbedoBoundary.geometry_map → IdentityMap()` **unconditionally**
  (`albedo.py:142-158`), *including* when the closure is `SpecularReturn`.
- `IdentityMap` (`_factors.py:442-470`) has `permutes_ordinates = False`,
  `domain_face(face) → face`. It carries nothing.
- `AngularAverageOperator` (`angular.py:41-352`) stores only `cos_w` (a `(|Γ₊|,)`
  float array), `norm`, `n_inflow`. Its `apply` is a contract + `np.broadcast_to`.
  **No index map, no permutation, no coordinate transform.**
- Its `axis`/`outward_sign` fields **look** geometric but are consumed at
  `angular.py:289-290` as `face_name(AXIS_NAMES.index(axis), outward_sign)` →
  `build_omega_dot_n` — i.e. purely a *face-orientation lookup* to select which
  hemisphere to integrate. That is F9, not a transformation of directions.
- The **one** place a transformation genuinely lives inside an albedo law is when
  its closure is `SpecularReturn` (`_factors.py:847-861`) → `SpecularReemission`
  (F16), and the code is explicit that this is the mirror in the **R** tier
  (`_factors.py:733-762`, `realizer.py:767-776`). That is by design, not hiding.

The theory page states the reason this misassignment was historically possible:
a rank-one `R` annihilates `G` (`R∘G = R` for any measure-preserving `G` when
`v = |Ω·n̂|` is `G`-invariant), so the `G` slot had no observable consequence —
`_factors.py:148-167`, `docs/theory/foundations/boundary_conditions.rst:568`.
`HemisphericalAverage` DID sit in the geometry slot until phase B3.0; it does not
now (`white.py:122-130`, `angular.py:50-56`).

### Q5 — Hard-coded axis literals where the axis is really a parameter?

**Yes, in exactly one shape: `axis: str = "x"` DEFAULTS on nine types.** Measured
(§G):

| type | file:line | default |
|---|---|---|
| `SpecularMirror` | `_factors.py:495` | `axis="x"` |
| `SpatialWrap` | `_factors.py:558` | `axis="x"` |
| `SpecularReturn` | `_factors.py:858` | `axis="x"` |
| `IsotropicReturn` | `_factors.py:881-882` | `axis="x", outward_sign=+1` |
| `LambertianReemission` | `_factors.py:699-701` | `axis="x", outward_sign=+1` |
| `SpecularReemission` | `_factors.py:774-775` | `axis="x"` |
| `ReflectiveBoundary` | `reflective.py:88` | `axis="x"` |
| `PeriodicBoundary` | `periodic.py:88` | `axis="x"` |
| `WhiteBoundary` | `white.py:108-109` | `axis="x", outward_sign=+1` |

These defaults are **a live hazard the codebase has already been bitten by twice**,
and it says so: `transport/method.py:296-303` — *"A law that declares its own
orientation MUST be listed here. The `law_cls()` fall-through hands such a law its
dataclass DEFAULTS, which are an orientation — and one that is right for exactly
one face."* White fell through until B3.4a (every face got `axis="x", sign=+1`);
periodic fell through until B3.4c (every face got `axis="x"`) —
`transport/method.py:319-344, 353-364`.

Other axis literals, all correct/derived:
- `symmetry.py:1030-1037` `_reflections`: `if axis == "x"` … — a dispatch, not a
  hardcode.
- `_specular.py:103-111` `_axis_cosines`: uses `AXIS_NAMES.index(axis)` with an
  explicit comment against a local `{"x":0,…}` literal. **Correct.**
- `angular_trace_space.py:184-199` `_quadrature_axis`: `if axis == 0/1/2` — a
  dispatch over the three columns. Correct but three-way hardcoded.
- `sn/operators/boundary.py:842, 915` — `self.sn_mesh.bc["xmax"]`, a **literal
  face name**. Justified in-place ("a seed-carrying mesh is 1-D curvilinear:
  exactly ONE boundary face"), so it is a documented 1-D-curvilinear invariant,
  not a parameter.
- `vacuum.py:42` — `np.flatnonzero(quad.mu_x < 0)` inside a **docstring example**.
  Hardcodes axis 0 AND uses a strict `< 0` rather than `TANGENTIAL_EPS`. Doc drift
  only (the same docstring also still claims the SN realization is an
  `IncomingOrdinateMaskTensor`, `vacuum.py:29-32` — retired at B3.2, see
  `realizer.py:690-708`).

### Q6 — Does anything in the BC layer construct a rotation?

**No. Not one.** Grep for `rotation|Rotation|rotational|SymmetryBoundary|sector`
over the whole BC layer returns **9 hits, all prose**:

- `geometry/boundary/__init__.py:25` and `_factors.py:83` — "a mirror, a spatial
  wrap, **a rotation**" in the `G`-membership description.
- `_factors.py:191-192` — the quotient table's "rotational (⅛-core) / by a finite
  rotation / cone points" row. Documented, unimplemented.
- `_factors.py:552-555` — *"Non-opposite gluing — a hex partner, a rotational
  quotient — is genuinely a different object and needs an explicit partner map.
  That is issue **#178** (`SymmetryBoundary`), deliberately NOT this type."*
- `realizer.py:498` and `_errors.py:125` — error messages pointing at #178.

So the rotational case is a **named, deliberately-deferred gap**, tracked at #178.
Note the asymmetry with `symmetry.py`: `Cn(n)` and `Dnh(n)` ARE fully realized
there (`symmetry.py:1092-1129`) with `_orbit_closure` able to produce the induced
node permutation — **the machinery for a rotational BC's angular action already
exists and has zero BC-layer consumers.** (Its spatial half — which sector maps
to which — would still need the non-`O(3)` partner map #178 describes.)

---

## 3. The three findings that most surprised me

### Surprise 1 — the unification the hypothesis predicts has ALREADY HAPPENED for the 3-D path, and it happened one day ago

`_compute_sphere_reflection_partners` (`directional.py:207-216`) does **not** do
its own nearest-neighbour search anymore. It builds the explicit reflection matrix
and hands it to `symmetry._orbit_closure`. I measured the two routes bit-identical:
`orbit_certificate(q.measure, SubgroupOfO3.Mirror(a)).permutations[σ]` ==
`q.reflection_index(axis)` on `product(2,4)`, `product(4,5)`, `lebedev(17)`, all
axes (§E).

The import comment at `directional.py:85-89` states the intent verbatim: *"The
reflection-partner map is a group action on the measure, so it is certified by the
SAME machinery that proves invariance — one source of truth."*

**So the audit's real answer is not "the codebase reflects by hand" — it is "the
codebase reflects by hand in exactly ONE remaining place (the 1-D GL path,
`directional.py:570-577`) and names the mirror in FOUR type vocabularies that do
not know about each other."** `SpecularMirror` / `SpecularReturn` /
`SpecularReemission` (`geometry/boundary/_factors.py`) and `symmetry.Mirror`
(`numerics/symmetry.py`) all carry the same `axis ∈ {x,y,z}` parameter, mean the
same `σ_a`, and `_factors.py` never imports `symmetry`.

### Surprise 2 — the 1-D path is the ONLY bit-exact one, and it is the one that bypasses the machinery

I expected the opposite. `identity[::-1]` at `directional.py:574` is a hardcoded
transcription with no guard whatsoever — and it is *structurally* exact:
`nodes[perm] == -nodes` holds bit-for-bit (§C). The 3-D path, which *is* routed
through the principled machinery, decides its permutation with
`np.argmin(dist, axis=1)` at `atol=1e-13` — tolerance in mechanism. It happens to
be exact in fact today (residual `0.0e+00`, runner-up margin 0.26–1.16, §B), so
the tolerance is doing no work — the same measured story as `TANGENTIAL_EPS`
(`angular_trace_space.py:76-91`) and as `singular_set`'s note that the three ad-hoc
epsilons "were never doing real work" (`symmetry.py:1624-1631`).

**This is the #325 connection.** If node generation moves from `cos`/`sin`
evaluation to algebraic construction, the 3-D partner map's residual stays `0.0`
and its runner-up margin stays O(0.3) — the map is not at risk. What IS
tolerance-sensitive is the *certification* layer built on top: F10's `np.sign`
comparison at `TANGENTIAL_EPS`, and F11's `rtol=1e-12`. And the exactness has a
second-order consequence already visible: `product(4,5)` loses axis 0 entirely
(`is_invariant(Mirror("x")) = False`, §E) — the odd-`n_φ` mirror-plane parity
documented at `directional.py:175-183`.

### Surprise 3 — the mirror is spelled a THIRD time as an integer `±1` key swap, and `PermutationOperator.is_involution` is silently false on the narrowed operator

`RadialCharacteristicBoundaryOperator._reflect_corner`
(`sn/operators/boundary.py:928-933`) applies the specular mirror to the
**off-quadrature** `μ = ±1` ray by writing

```python
out.corner(level, -1)[...] = seed.corner(level, +1)
```

— a literal sign-key swap in a Python loop, gated on `law_permutes_ordinates(law)`.
The docstring is explicit that this is the mirror ("the mirror of `μ = +1` is
exactly `μ = −1`, so the pairing is exact off-quadrature",
`sn/operators/boundary.py:882-884`), and that the OPERATOR route structurally
cannot reach it ("an (ordinate ⊗ group) operator over the QUADRATURE rows",
`:878-880`). So the same `Mirror(axis)` has one realization as a certified orbit
permutation and a second, disjoint realization as a hand-written `+1 → −1`
assignment for the one direction the quadrature omits. It also carries a
documented physics defect: the swap is **unscaled** — it does not multiply by
`R` (`sn/operators/boundary.py:897-907`), correct only at `α = 1`.

Second half of the surprise, found by probe (§H): the **narrowed** specular kernel
`PermutationOperator(local_perm)` built at `realizer.py:420` reports
`is_involution = False` for `lebedev(17)` (49×49 local perm) while reporting
`True` for `gauss_legendre(8)`, `product(2,4)`, `level_symmetric(6)`. That is
mathematically *correct* — since B3.2 the operator is a bijection `Γ₊ → Γ₋`
between two **different** index sets, so "involution" is a category error on it —
but `PermutationOperator`'s own docstring still asserts the opposite:

> "In particular, SN specular reflection through `Quadrature.reflection_index` is
> an involution; periodic shifts and rotational reorderings are not."
> — `orpheus/numerics/operator.py:2188-2191`

The claim is true of the GLOBAL `reflection_index` (measured True everywhere, §A)
and false of the NARROWED operator the SN realizer actually constructs. **Not a
live bug** — `is_involution` has zero production consumers (grep: only
`tests/numerics/test_permutation_operator.py:51-68`; the `derivations` hits are an
unrelated `derive_octant_frame_sign_is_involution`) — but it is a computed,
exported, docstring-blessed flag whose value silently depends on the quadrature.

---

## 4. Two smaller items worth recording

- **F5's involution is documented but never asserted.** `face_opposite`'s docstring
  calls it "an involution" (`face_layout.py:214`); I measured it true over all six
  faces (§D). The BC layer asserts the *angular* mirror's involution three ways
  (F12) and the *spatial* one zero ways. `_face_domains` (F6) asserts the weaker
  permutation property on the derived face map, which is what actually catches a
  half-declared periodic pair.
- **The odd-`n_φ` refusal escapes the BC error taxonomy.** `ReflectiveBoundary`
  and `AlbedoBoundary(·, SpecularReturn)` on `product(4,5)` raise a bare
  `ValueError` from `directional.py:420`, not a `BoundaryError` (§F), because
  `reflection_index` is consulted inside `_specular.py` before any BC-layer
  guard runs. Every other realization refusal in the layer is a `BoundaryError`
  carrying `law=`.

---

## Appendix — probe scripts

`scratch/_probe_bc_geom.py`, `scratch/_probe_bc_geom2.py`,
`scratch/_probe_bc_geom3.py` (untracked). Run with `.venv/bin/python`.
Sections referenced above: §A involution/fixed-point census · §B match-residual
and runner-up margin · §C 1-D exactness · §D `face_opposite` involution ·
§E `symmetry.Mirror` ≡ `reflection_index` · §F odd-`n_φ` refusal · §G axis
defaults · §H narrowed-kernel involutivity · §I tangential-band census.
