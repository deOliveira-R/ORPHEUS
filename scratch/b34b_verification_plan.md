# B3.4b verification plan — albedo's re-emission closure

> **STATUS 2026-08-01 — §11 (G) EXECUTED.** Every row of the migration table
> below is done except Row 12 (the generator), which is **STOPPED and
> documented** — see §11.2. Affected suite `tests/{geometry,sn/operators,
> diffusion}`: **9 failed → 0 failed** (`1543 passed, 5 skipped, 10 xfailed`).
> The §9.2 pairing-certification gate is **built and M13-verified RED** (5
> reds; control 7 passed). Mutation table re-run end-to-end: **9 of 9 red,
> control `239 passed, 1 xfailed`**. Two NET-NEW gates were added for
> production that moved mid-run — `law_permutes_ordinates` (§11.3) and the
> α-keyed refusal message (§11.4).

**Phase:** boundary-machinery B3.4b (`refactor/operator-strategy-layers`).
**Design of record:** `.claude/plans/b3_domain_narrowing_crosswalk.md` §11.1, §12.
**Status of the tree when this plan was measured:** production carve **already
written and uncommitted** by the main agent — `orpheus/geometry/boundary/`
{`_factors.py`, `albedo.py`, `reflective.py`, `__init__.py`, `_specular.py` (NEW)}
and `orpheus/sn/boundary/realizer.py`, on top of `7897c25c`. Every `[M]` below was
measured against THAT tree with `.venv/bin/python -O`, serial, 2026-08-01.
**This plan touches no file under `orpheus/`.**

Legend: `[M]` measured this session · `[R]` read from the tree · `[D]` design
judgement.

---

## 0. Executive summary — the four things that matter

### 0.1 The positional-pairing claim is **CONFIRMED**, and it is sharper than stated

`[M]` A bare `AlbedoBoundary(α)` inside `ι₋ ∘ law ∘ γ₊` produces exactly

```
out[sorted(inflow)[j]] = α · face_in[sorted(outflow)[j]]
```

verified `np.array_equal`-true on `gauss_legendre(4)` xmax, `gauss_legendre(8)`
xmin, `product(2,4)` xmax, `lebedev(17)` xmax, `level_symmetric(6)` xmax. The
pairing is by **array position** in two independently-sorted index sets, with no
geometric content. Use it in the docs.

**The sharpening — the coincidence is quadrature-dependent, and that is the
better argument.** `[M]` On `product(2,4)` xmax and `level_symmetric(6)` xmax the
positional pairing happens to **equal the specular pairing**
(`perm[sorted(inflow)] == sorted(outflow)`); on `gauss_legendre(4/8)` (slab — the
mirror *reverses* order) and on `lebedev(17)` it does **not**. So the old bare
albedo was not merely "unphysical": it silently agreed with a mirror on some
quadratures and disagreed on others, with nothing in the code or the tests
recording which. A configuration-dependent accident is a stronger justification
for the refusal than "no geometry", because it shows the answer was already
inconsistent across the tree.

**This is also a config-blindness trap for the NEW gates** — see §3.

### 0.2 (B) recommendation — pose the invariant over the **whole registry**, and SPLIT the two exceptions

Recommendation: **registry-wide parametrised gate, with `ReflectiveBoundary(axis,
α<1)` carried as `xfail(strict=True)` naming B5, and `AlbedoBoundary(1.0)` bare
carried as a NAMED, ASSERTED exception with its own positive test — not an
xfail, and not a silent exclusion.** Reasoning in §5; the short version:

- `[M]` the invariant has **two** violators today, not one. `§12.4` names
  `ReflectiveBoundary(axis, α<1)` (**two** non-trivial factors). It does **not**
  name `AlbedoBoundary(1.0)` bare, which has **zero** non-trivial factors
  (`G = IdentityMap`, `R = ScalarResponse(1.0) = I`). That row survives B3.4b
  because the diffusion arm still realizes it.
- `ReflectiveBoundary(α<1)` has a **landing phase** (B5), so `xfail(strict=True)`
  is exactly right — it is the campaign's own todo-list technique, and the
  XPASS-failure forces the marker's deletion when B5 lands. An *exclusion* would
  be invisible and B5 could land without anyone re-covering the row.
- `AlbedoBoundary(1.0)` bare has **no** landing phase — it is a permanent,
  correct member of the law algebra whose factorisation is degenerate *because*
  the scalar trace forces `G = id` by dimension. A permanent xfail is the wrong
  tool (lessons L4: an xfail with no unlock is a red-that-never-flips). Make it a
  positive assertion.
- Do **not** scope the gate to "tag-reachable laws". `[M]` `_law_from_tag`
  hard-codes `albedo=1.0` for reflective, so `ReflectiveBoundary(α<1)` is not
  tag-reachable — scoping that way *silently drops the very row §12.4 wants
  flagged*. And "tag-reachable" is method-dependent (`AlbedoBoundary` is absent
  from `SNMesh.BOUNDARY_OPERATOR_REGISTRY` but present in the diffusion one), so
  the invariant's scope would depend on which registry you consulted. The
  invariant is a property of the **law algebra**, not of a method's admission
  list.

### 0.3 (D) recommendation — RE-POSE by **inheriting the frozen `specular_x_lebedev17` artefact**; regenerate NOTHING

`[M]` `0.5 × snapshot["specular_x_lebedev17"]["psi_in"][inflow]` is
**bit-identical (maxdiff 0.0, non-vacuous)** to
`AlbedoBoundary(0.5, SpecularReturn("x"))` realized on `lebedev(17)` `xmax`.

So the albedo row does not need a regenerated baseline at all: it inherits a
**pre-B3.2 frozen artefact that was generated by a different law's realization**,
scaled by the α the test declares. That is strictly stronger than a regenerated
snapshot AND it encodes the `≡` theorem directly at the artefact level. Full
reasoning + the three re-baseline criteria applied explicitly in §7.

`[M]` The `albedo_05_lebedev17` artefact itself is **unusable** by any completion
(its `psi_in` is `0.5 × psi_out` on the full face, so its `Γ₋` rows read the
inflow ordinates the narrowed input never sees) — confirmed by direct comparison,
which returns `False`. It is retired, not re-posed.

### 0.4 What I think is WRONG or missing in §12 — five items

Detail in §12 of this document. Headlines:

1. **`§12.3`'s refusal criterion and `§11.1`'s exactly-one invariant are NOT the
   same statement, and the plan reads as if they were.** `[M]`
   `AlbedoBoundary(0.5)` bare *satisfies* "exactly one of G, R non-trivial"
   (`R = 0.5·I` is non-trivial) while being exactly the law SN refuses. The
   refusal criterion is about `R`'s **type** (scalar amplitude vs angular
   kernel), not about the non-triviality count. Gate them separately.
2. **The α = 0 edge of the refusal is undecided.** `AlbedoBoundary(0.0)` bare is
   physically complete on an angular trace (nothing returns; no pairing needed),
   yet the carve refuses it. I think refusing is *right* — but for a reason §12
   does not give, and the reason must be written down or the next reader will
   "fix" it. See §12.2.
3. **`§12.4` under-counts the invariant's violators** (see §0.2).
4. **The snapshot GENERATOR is already broken and nobody noticed.** `[M]`
   `python -m tests.geometry._generate_bc_equivalence_snapshots` raises on
   **6 of 8** cases today (B3.2/B3.4a narrowed the laws; the generator still
   calls `SNMethodSpace.minimal`). After B3.4b only `periodic` survives; after
   B3.4c, none. The harness's own failure message *instructs the reader to run
   it*. This is an inherited retirement gap that B3.4b's own commit surfaces.
   §9.4.
5. **The specular-closure orientation guard is weaker and less attributable than
   white's.** `[M]` `AlbedoBoundary(0.5, SpecularReturn("y"))` installed on
   `xmax` raises a raw `ValueError` from `TraceRestrictionOperator.to_local`
   (no `law=` attribution, not a `BoundaryError`), whereas the diffuse route
   raises an attributed `BoundaryError`. Two sibling routes, two error
   contracts. §12.5.

---

## 1. What is being verified — the term/claim enumeration (§1)

The carve introduces or changes exactly these claims:

| # | Claim | Where it lives `[R]` |
|---|---|---|
| C1 | `SpecularReemission(α, axis)` satisfies the `BoundaryResponseKernel` Protocol with `amplitude → α`, `is_zero → α == 0`, `is_adjointable → True` | `_factors.py` |
| C2 | `SpecularReturn(axis).kernel(α) == SpecularReemission(α, axis)`; `IsotropicReturn(axis, s).kernel(α) == LambertianReemission(α, axis, s)` — the closure is amplitude-FREE, so α has ONE home | `_factors.py` |
| C3 | `AlbedoBoundary.geometry_map is IdentityMap()` **unconditionally** (independent of the closure) | `albedo.py:124-140` |
| C4 | `AlbedoBoundary.response_kernel` = `ScalarResponse(α)` iff `reemission is None`, else `reemission.kernel(α)` | `albedo.py:142-154` |
| C5 | **T1**: `albedo(α, SpecularReturn(a)) ≡ reflective(a, α)` as realized operators | `realizer.py` albedo arm → `_specular_kernel` |
| C6 | **T2**: `albedo(α, IsotropicReturn(a, s)) ≡ white(a, s, α)` as realized operators | `realizer.py` albedo arm → `_checked_angular_average` |
| C7 | SN **REFUSES** `AlbedoBoundary(α, reemission=None)` at every α, naming the two completions | `realizer.py` albedo arm |
| C8 | Diffusion realizes `AlbedoBoundary(α)` **exactly as before**, and is INERT to the closure | `diffusion/boundary_realizer.py` (unchanged) |
| C9 | The realized albedo law's domain is `Γ₊` (the seven-law domain gate flips) | `realizer.py` + `_outflow_restriction` |
| C10 | The specular pairing's ERR-042/044/045 certification fires on BOTH carriers | `_specular.py` + `albedo.assert_realizable:224-227` |
| **C11** | **(net-new, unannounced)** `α = 0` on ANY kernel route now realizes to the narrowed **zero map** instead of dying in `ScaledOperator` | `realizer.py::_attenuated_kernel_operator` |
| **C12** | **(net-new, unannounced)** `reflective` and `white` at α=0 are now realizable — previously they RAISED | same |
| **C13** | **(net-new)** an unknown closure type raises a distinct `BoundaryError` naming the available closures | `realizer.py` albedo arm, final raise |

`[M]` C11/C12 are real behaviour changes the §12 design text does not mention.
`ScaledOperator(0.0, IdentityOperator())` still raises
`ValueError: ScaledOperator with zero scalar is degenerate; use ZeroOperator
explicitly.` — so before this carve, `ReflectiveBoundary(axis, 0.0)` and
`WhiteBoundary(..., 0.0)` were **legal laws that could not be realized**. They
now realize. That is a fix, and it needs its own gate (§9.1) because nothing in
the inherited suite covers it.

---

## 2. Claim layer, pillar, structural independence (§1.5 gate)

Every row here is a **structural / algebraic** claim about operator
construction — none is a physics claim, so none carries `verifies(<equation>)`.

| Gate family | Claim layer | Pillar | Truth source | Structurally independent of the SUT? |
|---|---|---|---|---|
| A1 `≡` theorems (T1/T2) | *neither* — an **identity between two construction routes** | — (regression / bit-identity inheritance) | the sibling route | **NO — and that is admitted.** Both routes execute one body, so this gate can only catch ARGUMENT drift. It is *necessary, not sufficient*; it must be PAIRED with A2/A3. |
| A2/A3 independent-expression anchors | flux-shape (of the boundary operator) | closed-form | hand-written numpy from `quad.reflection_index` / `quad.weights` × `omega_dot_n` | **YES** — nothing imported from the realizer above the `np.take` / `np.einsum` line. |
| A4 snapshot inheritance | regression | frozen artefact | `bc_equivalence_specular_x_lebedev17.npz` (generated PRE-B3.2, by a *different law*) | **YES** — the artefact predates every line under test. |
| B G/R invariant | algebraic invariant over the law table | closed-form (the §11.1 table) | the table, transcribed | **YES** — the table is prose the code must satisfy. |
| C refusal | contract | — (negative gate) | the production entry point + `match=` on the specific message | **YES** provided the raise comes from `SNBoundaryRealizer.realize` (vv Mode-8 class 5). |
| D domain gate | typing | — | `\|Γ₊\|`-row probe + the full-face non-endomorphism leg | **YES** (Leg B leaves the shape functional — vv Mode 12). |

**Nothing here is an eigenvalue claim, so no MMS/eigenvalue mismatch is
possible.** The `≥2G` cardinal rule does not bite (there is no eigenvalue), but
every probe below carries a **trailing group axis** so a group↔ordinate axis
transpose cannot hide.

**Tolerance ruling — `np.array_equal`, everywhere, for A1/A2/A4.** Justified
structurally, not measured-and-rounded-up:

- T1/T2 route through **one construction body** with identical arguments, so the
  two operators hold the *same* `local_perm` / `AngularAverageOperator` object
  data and their `apply` executes the *identical* reduction tree. The predicted
  FP drift is **exactly 0** (vv §bit-identity criterion 3: the bound is
  `reduction_depth × ULP`, and the reduction depth *difference* is zero). Any
  non-zero drift means the routes did **not** execute the same body — which is
  precisely the bug this gate exists to catch, so `assert_allclose` at any rtol
  would be a gate downgrade that admits the failure mode.
- A2 compares a gather / scaled-gather against the same operation written by
  hand in the same order → predicted **0 ULP**, so `array_equal`.
- **A3 is the one exception.** `[R]` `AngularAverageOperator.apply` computes
  `(cos_w[:,None,…] * psi).sum(axis=0) / norm` with
  `cos_w = w[Γ₊] * (Ω·n)[Γ₊]` and `norm = cos_w.sum()`; a hand-written
  `np.tensordot(cos_w, psi, axes=(0,0))` is the same value by a **different
  reduction order**. Bound = `|Γ₊| × ULP`; on `lebedev(17)` that is
  `49 × 2.2e-16 ≈ 1.1e-14`, so **`assert_allclose(rtol=1e-14, atol=0)` IS the
  bound, not a rounded-up measurement**. `[M]` measured max relative error:
  `5.1e-16` (GL4), `7.1e-16` (product 2×4), `2.4e-16` (LS6), **`3.2e-15`**
  (lebedev 17) — inside the bound, and the largest is on the deepest reduction,
  which is the predicted shape. `[D]` Do not promote A3 to `array_equal`; do
  not loosen it past `1e-14` without recomputing the depth.
- A4 is `array_equal` for the same reason the existing specular snapshot rows
  are (`test_bc_equivalence_snapshot.py:268`, `:307`): a scale-then-gather with
  no re-association. `[M]` verified maxdiff 0.0.

---

## 3. Config-blindness declaration (§0.6) — the MEASURED fixture matrix

**This is the load-bearing section.** `[M]` For the specular route, the
"forgot the `to_local` remap / used `arange`" mutation is **INVISIBLE** on two of
the five candidate fixtures:

| fixture | `γ₊.to_local(perm[inflow])` | equals `arange`? | verdict |
|---|---|---|---|
| `gauss_legendre(4)` xmax | `[1, 0]` | **False** | **DISCRIMINATES** (slab mirror reverses order) |
| `gauss_legendre(8)` xmin | `[3, 2, 1, 0]` | **False** | **DISCRIMINATES** |
| `product(2,4)` xmax | `[0, 1]` | **True** | **BLIND** — keep as a documented CONTROL |
| `lebedev(17)` xmax | `[0,1,2,3,4,5,6,7,…]` but reorders at positions 9–12 | **False** | **DISCRIMINATES** |
| `level_symmetric(6)` ymax | `[0,1,…,23]` | **True** | **BLIND** — keep as a documented CONTROL |

**Binding requirement (specular):** every specular gate that claims to catch a
remap bug MUST run on `gauss_legendre` **and** `lebedev(17)`. Running only on
`level_symmetric` or `product` — the two most "natural" multi-D choices — is
designed-green.

### ⭐ The two routes need COMPLEMENTARY fixtures — neither list covers the other

`[M]` For the **diffuse** route the discriminating fixture is the exact
complement. Re-introducing the B3.4a twin (classify `Γ₊` with `> 0.0` instead of
`> TANGENTIAL_EPS`) admits extra ordinates on:

| fixture | `\|Γ₊\|` | `\|{Ω·n > 0}\|` | mis-admitted | verdict for the DIFFUSE gate |
|---|---|---|---|---|
| `gauss_legendre(4)` xmax | 2 | 2 | — | **BLIND** |
| `product(2,4)` xmax | 2 | 4 | ordinates 1, 5 | **DISCRIMINATES** |
| `lebedev(17)` xmax | 49 | 49 | — | **BLIND** |
| `level_symmetric(6)` ymax | 24 | 24 | — | **BLIND** |

So:

| mutation class | discriminates on | blind on |
|---|---|---|
| specular remap (`to_local` → `arange`) | `gauss_legendre`, `lebedev(17)` | `product(2,4)`, `level_symmetric(6)` |
| diffuse classification (`TANGENTIAL_EPS` → `> 0.0`) | **`product(2,4)`** | `gauss_legendre`, `lebedev(17)`, `level_symmetric(6)` |

**Neither gate's fixture set can be dropped in favour of the other's, and a
single "representative quadrature" cannot serve both.** This is exactly the
crosswalk §9.1 finding ("`to_local`, two sites, complementary fixtures — neither
covers the other") reproduced one phase later from the albedo side. Parametrise
BOTH routes over the full four-fixture list and keep the blind rows as
documented controls; do not prune.

Other blindnesses to break, per gate:

- **Trailing axes.** Every probe is `(n_out, ng, nx)` with `ng != nx` (e.g.
  `(n_out, 4, 2)`), so a group↔spatial transpose reds.
- **α coverage.** `α ∈ {0.0, 0.3, 1.0}` on EVERY equivalence row — the three
  values are three DIFFERENT production branches (`_narrowed_zero_operator` /
  `ScaledOperator` / bare TP). `[M]` `α=0` was the branch that could not be
  realized at all before this carve.
- **Both axes.** At least one row on `axis="y"` (`level_symmetric(6)` ymax), so
  an `axis` hard-code to `"x"` reds.
- **Both faces.** At least one `min` face (`gauss_legendre(8)` xmin,
  `outward_sign=-1`), so an `outward_sign` hard-code to `+1` reds.
- **Random, non-flat probe.** A flat probe makes the Lambertian average equal to
  the input and nulls the whole redistribution; use `default_rng(seed).normal`.
- **Non-vacuity.** Every bit-identity row asserts `np.count_nonzero(expected) > 0`
  — an all-zero comparison passes for any broken implementation. (α=0 rows are
  the exception BY DESIGN and must assert zero-ness + shape instead.)

---

## 4. (A) The `≡` theorems — gate design

### A1 — the route-equivalence gate (necessary, NOT sufficient)

New file `tests/geometry/test_reemission_closure_equivalence.py`
(`pytestmark = [pytest.mark.foundation]`).

Parametrised over `(quadrature, face, axis, outward_sign) × α ∈ {0.0, 0.3, 1.0}`
using the §3 fixture list, three legs per row:

1. **forward** `np.array_equal(albedo_op.apply(x), sibling_op.apply(x))`
   with `x = rng.normal(size=(n_out, 4, 2))`;
2. **transpose** — `adjointable(albedo_op) == adjointable(sibling_op)` FIRST
   (a capability drop is silent and gates the whole composite through
   `SNBoundaryOperator.is_adjointable`), then, when both are adjointable,
   `np.array_equal(albedo_op.apply_transpose(y), sibling_op.apply_transpose(y))`;
3. **non-vacuity** — `count_nonzero(sibling_op.apply(x)) > 0` for α ≠ 0.

`[M]` Current tree: all 15 (quad × α) rows pass both legs for T1 and T2,
**including the Python type** (`ZeroOperator` / `ScaledOperator` /
`TensorProductOperator` on both sides). `[M]` Adjointability: the specular route
is adjointable at every α; the diffuse route is adjointable **only at α=0**
(`LambertianReemission.is_adjointable is False`) — and *both sides agree*, which
is the assertion. Do **not** assert `adjointable(...) is True` unconditionally;
assert the two routes AGREE, plus pin the known value in a separate row so a
future B5 flip is visible.

`[D]` **Do NOT assert Python-type equality** as the primary claim. It is true
today `[M]` but it over-constrains: a legitimate future fast-path difference on
one route would red a gate that claims to be about the *action*. Pin the type
only in the dedicated `test_operator_block_role` rows (§11), which is where type
claims already live.

### A2 — the specular INDEPENDENT-EXPRESSION anchor (the real catcher)

Same file. The reference is written from `quad.reflection_index(axis)` and the
raw index sets — it imports nothing from `orpheus.sn.boundary.realizer`:

```python
perm = quad.reflection_index(axis)
inflow = np.sort(space.inflow_indices)
outflow = np.sort(space.outflow_indices)
full = np.zeros((quad.N,) + x.shape[1:]); full[outflow] = x
expected = alpha * full[perm[inflow]]          # row j reads the mirror of inflow[j]
```

`np.array_equal(op.apply(x), expected)`, over the §3 fixture list. **MUST include
`gauss_legendre` and `lebedev(17)`** (§3). The `product(2,4)` /
`level_symmetric(6)` rows stay in the parametrisation as **documented control
legs** with a docstring note that the remap mutation is invisible there — a
control that pins the coincidence, not filler.

This is the gate that catches a shared-body bug (A1 cannot: both routes carry the
same wrong permutation).

### A3 — the diffuse INDEPENDENT-EXPRESSION anchor

```python
odn  = build_omega_dot_n(quad, (face,))[0]          # face-name -> signed projection
out_w = (quad.weights * np.abs(odn))[outflow]       # NB: Γ₊ only, TANGENTIAL_EPS-classified
current = np.tensordot(out_w, x, axes=(0, 0))       # cosine-weighted outflow current
expected = np.broadcast_to(alpha * current / out_w.sum(), (inflow.size,) + current.shape)
```

`assert_allclose(op.apply(x), expected, rtol=1e-14, atol=0)` — reduction-order
difference, bounded and measured in §2.

`[R]` **The normalisation is confirmed against the source**
(`orpheus/sn/boundary/angular.py:308-314, 343-347`): `norm = cos_w.sum()` over
`Γ₊` only, i.e. the conservative "total returned current = α × incident current"
convention, **not** the textbook `1/π`. `WhiteBoundary`'s own docstring says the
same (`white.py:44-51`). Use `Σ_{Γ₊} w|Ω·n|`.

Mode-7 declaration for A3: the ansatz **activates** the cosine weighting, the
Γ₊-restricted sum, the α fold and the Γ₋ broadcast; it **nulls** nothing (a
random non-flat probe with distinct trailing axes). A *flat* probe would null the
whole averaging term — do not use one.

---

## 5. (B) The exactly-one-of-`G`,`R` invariant — gate design

### 5.1 The decidable predicate

```python
def _g_trivial(G) -> bool:  return isinstance(G, IdentityMap)
def _r_trivial(R) -> bool:  return isinstance(R, ScalarResponse) and R.alpha == 1.0
```

`[D]` `R` trivial means `R = I`, i.e. "no physics asserted" — that is §11.1's own
wording (*"`R = I` exactly when the BC is a pure symmetry statement"*).
`ScalarResponse(0.0)` (vacuum) is **non**-trivial: returning nothing IS a
constitutive statement.

### 5.2 MEASURED state of the whole registry `[M]`

| law | `G` | `R` | # non-trivial |
|---|---|---|---|
| `VacuumInflow()` | `IdentityMap` | `ScalarResponse(0.0)` | 1 ✓ |
| `ReflectiveBoundary(x, 1.0)` | `SpecularMirror(x)` | `ScalarResponse(1.0)` | 1 ✓ |
| **`ReflectiveBoundary(x, 0.7)`** | `SpecularMirror(x)` | `ScalarResponse(0.7)` | **2 ✗** |
| `WhiteBoundary(x, +1, α)` | `IdentityMap` | `LambertianReemission(α,…)` | 1 ✓ |
| `AlbedoBoundary(0.0)` | `IdentityMap` | `ScalarResponse(0.0)` | 1 ✓ |
| `AlbedoBoundary(0.5)` | `IdentityMap` | `ScalarResponse(0.5)` | 1 ✓ |
| **`AlbedoBoundary(1.0)`** | `IdentityMap` | `ScalarResponse(1.0)` | **0 ✗** |
| `AlbedoBoundary(α, SpecularReturn(x))` | `IdentityMap` | `SpecularReemission(α,x)` | 1 ✓ |
| `AlbedoBoundary(α, IsotropicReturn(x,+1))` | `IdentityMap` | `LambertianReemission(α,…)` | 1 ✓ |
| `PeriodicBoundary()` | `SpatialWrap(x)` | `ScalarResponse(1.0)` | 1 ✓ |
| `PrescribedInflow(src)` | `IdentityMap` | `ScalarResponse(0.0)` | 1 ✓ |
| `ZeroFluxBoundary()` | `IdentityMap` | `ScalarResponse(-1.0)` | 1 ✓ |

### 5.3 The gate

New class `TestExactlyOneFactorIsNonTrivial` in
`tests/geometry/test_boundary_factors.py` (its native home — that file already
owns the factor tier).

- **Row set** = the table above (every registry law × representative amplitudes,
  including both closures and both edge α).
- `ReflectiveBoundary(x, 0.7)` → `pytest.mark.xfail(strict=True, reason=<B5>)`.
  Reason text names: the §11.1 discriminator, that this object *is*
  `AlbedoBoundary(0.7, SpecularReturn("x"))` wearing the geometry costume, that
  it is unreachable from a tag (`_law_from_tag` hard-codes `albedo=1.0`), and
  that **B5 retires `ReflectiveBoundary.albedo`** — so deleting the marker is
  part of B5's definition of done.
- `AlbedoBoundary(1.0)` bare → **NOT in the parametrisation.** It gets its own
  positive test:

  ```python
  def test_bare_albedo_at_unit_amplitude_carries_NEITHER_factor(self):
      """G = I and R = I: the one registry member that asserts nothing.

      Not a defect and not deferred — it is the degenerate factorisation a
      SCALAR trace forces (G is the identity by dimension there), which is
      exactly why an ANGULAR method cannot realize it and refuses. Pinned
      positively so the §11.1 invariant's scope is honest rather than
      quietly excluding its own counter-example.
      """
  ```

  and that test cross-links to the SN refusal gate (§6) in its docstring.
- **Inventory-completeness leg** (mirroring
  `test_b3_domain_narrowing.py::test_the_law_inventory_is_complete`): assert
  `{covered rows} ∪ {named exceptions} ⊇ {type(BoundaryTraceLaw.create(k)) for k
  in BoundaryTraceLaw.registry}`. Without it, a law added tomorrow escapes the
  invariant silently — and the exclusion of `AlbedoBoundary(1.0)` would be
  indistinguishable from an oversight.

**Teeth (§F):** mutate `AlbedoBoundary.geometry_map` to return
`SpecularMirror(self.reemission.axis)` when the closure is specular — i.e. put
the mirror back in `G`, the exact thing §11.1 forbids. The specular-closure rows
go to 2 non-trivial factors and RED. Positive control: the `IsotropicReturn` and
`WhiteBoundary` rows stay green, proving the mutation was scoped and did bite.

---

## 6. (C) The refusal gate — both directions, with teeth

### 6.1 SN refuses — `tests/sn/operators/test_sn_boundary_realizer.py`

Replaces the retired `TestRealizeAlbedo` (§11). Requirements, each mapping to a
named Mode-8 class:

```python
@pytest.mark.parametrize("alpha", [0.0, 0.5, 1.0])
def test_bare_albedo_is_refused_naming_both_completions(self, alpha):
    quad = Quadrature.gauss_legendre(8)
    space = face_method_space(quad, face="xmax")
    with pytest.raises(BoundaryError) as exc:                 # PRODUCTION entry point
        SNBoundaryRealizer().realize(AlbedoBoundary(alpha), space)
    assert exc.value.law == "albedo"                          # attribution
    msg = str(exc.value)
    assert "SpecularReturn" in msg and "IsotropicReturn" in msg   # names BOTH completions
    assert "ARRAY POSITION" in msg.upper()                    # names the actual defect
```

- **Mode-8 class 5 (self-satisfied `raises`)**: the raising call is
  `SNBoundaryRealizer().realize`, never a locally-built exception. Non-negotiable.
- **The `match=`/message legs are load-bearing**, not decoration (lessons L4): a
  bare `pytest.raises(BoundaryError)` is satisfied by the *`_outflow_restriction`*
  guard, by `assert_realizable`, or by any other refusal in the arm — a false
  green. The three message assertions pin THIS raise. `[M]` the current
  production message contains all three substrings.
- **α coverage is mandatory** — `[M]` all three α were separate code paths
  pre-carve (`ZeroOperator` / `IdentityOperator` / `ScaledOperator`), and a
  refusal added to only the general branch would leave the two fast paths live.
- **Positive CONTROL, same test class**: the SAME α with each closure realizes
  and produces a `Γ₋`-shaped image. Without it, an albedo arm that raises for
  *every* input also passes.
- **Ordering hazard** `[D]`: `ZeroFluxBoundary` also satisfies "G is identity and
  R is a non-zero scalar". Its refusal branch runs FIRST and is pinned by
  `test_b3_domain_narrowing.py::test_zero_flux_is_refused_outright`
  (`exc.value.law == "zero_flux"`). Keep that test; if the implementer ever
  generalises the albedo refusal into a factor-derived predicate, that test is
  the tripwire that catches the branch reorder.
- **Unknown-closure leg (C13)**: a duck-typed object with a `kernel` method that
  is neither `SpecularReturn` nor `IsotropicReturn` must raise the *other*
  `BoundaryError` and the message must name both available closures. This is a
  real reachable path (`ReemissionClosure` is a `runtime_checkable` Protocol, so
  a third-party shape satisfies it structurally).

### 6.2 Diffusion does NOT refuse — `tests/diffusion/test_boundary_realizer.py`

`[M]` The diffusion arm today returns, for `DiffusionMethodSpace()` and probe
`[1.0, 2.0]`:

| law | operator | `apply` |
|---|---|---|
| `AlbedoBoundary(0.5)` | `ScaledOperator` | `[0.5, 1.0]` |
| `AlbedoBoundary(0.0)` | `ZeroOperator` | `[0.0, 0.0]` |
| `AlbedoBoundary(1.0)` | `IdentityOperator` | `[1.0, 2.0]` |
| `AlbedoBoundary(0.5, SpecularReturn("x"))` | `ScaledOperator` | `[0.5, 1.0]` |
| `AlbedoBoundary(0.5, IsotropicReturn("x", +1))` | `ScaledOperator` | `[0.5, 1.0]` |

Two NET-NEW gates:

- **CD-1 (the arm does not move)** — pin the first three rows' operator TYPE and
  `apply` value bit-exactly. This is the "must not move" claim, and it currently
  has no explicit α∈{0,1} coverage of the *bare* law tied to the B3.4b change.
- **CD-2 (the closure is INERT on a scalar trace)** — `[D]` **this is the gate
  §12 does not ask for and most needs.** Assert
  `diffusion.realize(AlbedoBoundary(α, C)).apply(v)` is `array_equal` to
  `diffusion.realize(AlbedoBoundary(α)).apply(v)` for `C ∈ {None, SpecularReturn,
  IsotropicReturn}` and α ∈ {0, 0.5, 1}. Rationale: `§12.3` says the scalar trace
  has one dof so `α·I` IS the whole story — that is a *claim about the diffusion
  arm* and nothing pins it. A future "read `response_kernel` instead of
  `amplitude`" tidy-up would silently change diffusion answers; CD-2 is the only
  thing that would catch it.
  - Teeth: monkeypatch the diffusion realizer's albedo branch to dispatch on
    `type(law.response_kernel)` → CD-2 reds.

---

## 7. (D) The snapshot re-pose — recommendation and reasoning

### 7.1 Recommendation

1. **DELETE** `TestAlbedo05Lebedev17Snapshot`
   (`tests/geometry/test_bc_equivalence_snapshot.py:207-233`). The claim it pins
   ("realizing a bare albedo on a FACELESS method space yields a full-face
   `α·I`") ceases to be a production behaviour — it is refused. There is no
   successor for THAT claim.
2. **ADD**, inside the existing `TestSpecularXLebedev17Snapshot`, a second
   method that inherits the **same frozen artefact**:

   ```python
   def test_albedo_with_specular_closure_matches_the_same_frozen_image(self, snapshot):
       r"""The ≡ theorem, pinned at the FROZEN artefact.

       ``AlbedoBoundary(α, SpecularReturn("x"))`` and ``ReflectiveBoundary("x", α)``
       realize to the same matrix (B3.4b), so the albedo route's image is the
       SAME pre-B3.2 frozen ``psi_in``, restricted to Γ₋ and scaled by α. The
       reference is therefore a snapshot generated BEFORE the narrowing, by a
       DIFFERENT law — nothing is regenerated, and the theorem is asserted
       against an artefact that predates every line under test.
       """
       quad  = Quadrature.lebedev(17)
       space = face_method_space(quad, face="xmax")
       op = SNBoundaryRealizer().realize(
           AlbedoBoundary(0.5, SpecularReturn(axis="x")), space,
       )
       actual   = op.apply(snapshot["psi_out"][space.outflow_indices])
       expected = 0.5 * snapshot["psi_in"][space.inflow_indices]
       assert np.count_nonzero(expected), "vacuous comparison"
       np.testing.assert_array_equal(actual, expected)
   ```

   `[M]` bit-identical, maxdiff **0.0**, non-vacuous, on the current tree.
3. **RETIRE** `tests/geometry/snapshots/bc_equivalence_albedo_05_lebedev17.npz`
   **and** its generator entry
   (`_generate_bc_equivalence_snapshots.py:166-171`) together. An orphaned
   baseline that `generate_all` keeps re-creating is exactly the half-alive
   artefact `coding-standards.md` forbids.
   `[D]` This is the one item I would put to the user as a checkpoint: deleting a
   committed baseline is irreversible-in-spirit even if git remembers it. The
   alternative (keep the `.npz`, keep the generator case, mark the case
   "diffusion-only") is defensible but leaves a snapshot no test reads.

### 7.2 The three re-baseline criteria — applied explicitly

The recommendation **regenerates nothing**, so the criteria are formally not
engaged. Applied anyway, to show why the *alternative* (regenerate `psi_in` from
the new albedo code) is inferior:

| criterion | inherit-the-specular-artefact (recommended) | regenerate `albedo_05` |
|---|---|---|
| 1. principled at every step | ✓ — every intermediate is named (`psi_out`, the frozen specular image, the α fold) | ✓ (same body) |
| 2. verified against a **structurally-independent** reference | ✓ **the artefact IS the independent reference** — frozen pre-B3.2, produced by `ReflectiveBoundary`, restricted and scaled here | ✗ — a baseline regenerated from the very code it gates is not a reference at any level; the file's own B0.3 RELABEL paragraph (`test_bc_equivalence_snapshot.py:45-66`) already says so |
| 3. drift is FP-non-associativity, dimensionally explained | ✓ — predicted drift 0, `[M]` measured 0 | n/a (no old value to drift from) |

**Is this "re-baselining to make a test pass"?** No, and the discriminator is
sharp: nothing is being *silenced*. The old claim was **retired by refusal**, not
found failing; the new claim is a *different* claim (the `≡` theorem) whose
reference is an artefact that existed before the phase started. Re-baselining
would be taking today's albedo output and freezing it — which is precisely what
this avoids.

### 7.3 What is lost, honestly

The `albedo_05_lebedev17` row was one of the eight in "the widest mutation net in
the boundary subsystem (9 of 12 leaf mutations)". `[D]` Its contribution to that
net was the α-fold on a pass-through operator — which the new row still carries
(the `0.5 ×` factor), on the same quadrature, plus a permutation the old row did
not exercise. So the net does not shrink; it gains the specular pairing on a row
that previously tested only a scalar multiply.

**Optional add-on** (`[D]`, user's call): a NEW
`albedo_isotropic_05_lebedev17` case would be genuinely new coverage — `lebedev(17)`
carries **12 tangential ordinates per axis** (the most of any production
quadrature `[M]` per crosswalk §13.2) and no existing white snapshot runs on it.
That case *would* need a regenerated baseline, and it is admissible under the
three criteria **only if it ships in the same commit as the A3 independent
hand-numpy anchor**, which is what supplies criterion 2. I would not spend it in
B3.4b — it belongs with B4's factor-dispatch work.

---

## 8. (E) The domain gate flip — `tests/sn/operators/test_b3_domain_narrowing.py`

### 8.1 The finding that makes this NOT optional bookkeeping

`[M]` **Under the refusal alone, with `_LAWS` unchanged, the three `albedo_*`
rows keep xfailing — for the WRONG REASON.** Measured with `--runxfail`: they
fail on `BoundaryError` out of `realize`, never reaching the documented
"emits N rows / is an endomorphism" assertion. That is vv **Mode-8 class 4**, the
MISATTRIBUTED strict xfail: the rows look like committed coverage of the
narrowing, assert nothing about it, and — because they can never XPASS — the
"the xfail set IS the todo list" mechanism silently stops working for exactly
those rows. B3.4b could land "complete" with three permanently-green dead gates.

### 8.2 The edit

`_LAWS` (lines 276-290) — replace the three albedo rows with **four** completed
rows, all `deferred=False`:

```python
"albedo_specular_0":   (AlbedoBoundary(0.0, SpecularReturn(axis="x")), False),
"albedo_specular_1":   (AlbedoBoundary(1.0, SpecularReturn(axis="x")), False),
"albedo_specular_05":  (AlbedoBoundary(0.5, SpecularReturn(axis="x")), False),
"albedo_isotropic_05": (AlbedoBoundary(0.5, IsotropicReturn(axis="x", outward_sign=+1)), False),
```

The fixture is `gauss_legendre(8)` face `xmax` (line 335-336), so
`axis="x", outward_sign=+1` matches and the diffuse orientation cross-check is
green by construction. **Keep all three α** — the α=0 row exercises the NEW
`_narrowed_zero_operator` route (C11) and is the only place in this file that
does.

`_B34_XFAIL` (lines 256-270) — rewrite so it describes **periodic only**. Leaving
the albedo clause in makes the surviving marker a misattributed *reason*
(documents albedo, applies to periodic) — the same class of defect the crosswalk
§13.4 already found in `_MIXED_LAW_XFAIL`.

### 8.3 Is a "bare albedo" row still meaningful in this gate? — YES, with a changed claim

It stops being a **domain** claim and becomes a **refusal** claim, and the file
already has the template. Add, immediately after
`test_zero_flux_is_refused_outright`:

```python
def test_bare_albedo_is_refused_outright(self) -> None:
    """A bare ``AlbedoBoundary`` is not a deferred narrowing — it is
    under-determined on an angular trace, and SN says so.

    Pinned HERE, next to zero_flux, so the parametrised gate's law
    inventory stays provably complete: registry = the rows above ⊎
    {zero_flux, bare albedo}. Without it the disappearance of the three
    ``albedo_*`` xfails would be indistinguishable from an oversight.
    """
```

asserting `exc.value.law == "albedo"` and the two-completion message substrings
(§6.1). `test_the_law_inventory_is_complete` (line 379) needs **no change** —
`AlbedoBoundary` is still covered by the four completed rows — but its docstring
should note that the bare spelling is covered by the refusal test, mirroring how
`ZeroFluxBoundary` is handled at line 389.

---

## 9. Net-new gates the carve OPENED (not in the §12 brief)

### 9.1 C11/C12 — `α = 0` now realizes to the narrowed zero map

`[M]` `ScaledOperator(0.0, …)` raises `ValueError: ScaledOperator with zero
scalar is degenerate`. Before this carve, `ReflectiveBoundary(axis, 0.0)` and
`WhiteBoundary(axis, sign, 0.0)` were legal laws (α=0 satisfies every invariant
incl. sub-Markov) that **died in the numerics layer**. `[M]` They now realize to
`ZeroOperator` with `Γ₋`-shaped, all-zero output and a working transpose.

New gate in `tests/sn/operators/test_sn_boundary_realizer.py`, parametrised over
`{reflective, white, albedo+specular, albedo+isotropic} × α = 0`:

- returns a `ZeroOperator`;
- `apply` emits `|Γ₋|` rows, all zero;
- `apply_transpose` emits `|Γ₊|` rows, all zero — the asymmetric-space leg. On a
  fixture with `|Γ₊| == |Γ₋|` this is Mode-12 blind on shape alone, so ALSO
  assert `adjointable(op) is True` and that the transpose of a **non-zero** Γ₋
  probe is zero (the zero map's defining law, not its shape);
- `[D]` and — the leg that actually pins the fix — a **negative control**:
  `ScaledOperator(0.0, IdentityOperator())` still raises, so the gate is
  demonstrably about the realizer's routing, not about a numerics change.

**Teeth:** delete the `if alpha == 0.0:` early return in
`_attenuated_kernel_operator` → all four rows red with the `ScaledOperator`
`ValueError`. Positive control: the α=0.3 and α=1.0 rows stay green.

### 9.2 C10 — the specular pairing is certified on BOTH carriers

`[R]` `albedo.assert_realizable:224-227` fires `assert_specular_pairing_valid`
only when `isinstance(self.reemission, SpecularReturn)`. This closes a hole the
phase would otherwise open (a wrong `reflection_index` table caught on one route,
silently realized on the other). It needs a gate, because "both carriers fire it"
is a claim with no test today.

Gate in `tests/geometry/test_bc_universal_invariants.py` (its native home): a
monkeypatched `Quadrature.reflection_index` returning a table that breaks each of
the three invariants in turn (mispaired weight class → ERR-042; self-map →
ERR-045; 3-cycle → ERR-044) must raise the SAME error type from BOTH
`ReflectiveBoundary(axis)` **and** `AlbedoBoundary(α, SpecularReturn(axis))`,
each with its own `law=` attribution (`"reflective"` vs `"albedo"`).

**Teeth:** delete the `if isinstance(self.reemission, SpecularReturn):` clause in
`albedo.assert_realizable` → the three albedo rows go green-when-they-should-red
(i.e. the gate reds). Positive control: the reflective rows still raise.

`[D]` A cheaper structural companion worth having: assert
`AlbedoBoundary(α, SpecularReturn(a)).assert_geometry_map_measure_preserving(quad)`
is a **no-op** even on a broken table — because this law's `G` really is the
identity. That pins the *tier* claim (`albedo.py:214-220`), which is the thing
§11.1 is about, and it would red if someone "fixed" the hole by overriding the
`G` hook instead of extending `assert_realizable`.

### 9.3 The `_checked_angular_average` generalisation

`[R]` The helper lost its `law: WhiteBoundary` parameter and now takes
`axis` / `outward_sign` / `law_key` keywords. The ERR-041-pattern orientation
guard must be shown to still bite **from the albedo side**:

`AlbedoBoundary(0.5, IsotropicReturn(axis="y", outward_sign=+1))` installed on
`xmax` must raise `BoundaryError` with `law == "albedo"` (not `"white"` — the
attribution is now parametrised and a hard-coded `law="white"` would be a silent
mis-attribution). Positive control: the matching orientation realizes.

**Teeth:** hard-code `law_key="white"` in the raise → the albedo row reds on the
`exc.value.law` assertion.

### 9.4 The snapshot GENERATOR — an inherited retirement gap this phase surfaces

`[M]` `python -m tests.geometry._generate_bc_equivalence_snapshots` **raises on
6 of 8 cases today** (`vacuum`, `specular_x`, `specular_y`, `white_xmax`,
`white_xmin`, `mixed` — all `BoundaryError: … without outflow_indices`), because
it still builds `SNMethodSpace.minimal(quad)` for laws B3.2/B3.4a narrowed. Only
`albedo_05` and `periodic` run — exactly the two un-narrowed laws. **After B3.4b
only `periodic` runs; after B3.4c, none.**

Meanwhile `_load_snapshot`'s failure message
(`test_bc_equivalence_snapshot.py:136-143`) tells the reader to run it. That is a
present-tense-false instruction in a test-facing error path — a doc bug under the
standing "a falsified doc is fixed in the SAME change" directive.

`[D]` **Recommendation:** B3.4b fixes the generator's `_build_payload` to use a
face-ful space (the same `face_method_space` shape the harness uses) for every
narrowed case, and drops the `albedo_05` case. This is ~15 lines and it restores
a claim the harness makes about itself. If the main agent judges it out of scope,
it MUST be filed as an issue in the same commit — a generator that cannot
generate is exactly the half-alive artefact the retirement doctrine targets.
`[D]` I lean **fix it now**: the phase is already touching the harness, and
B3.4c will hit the same wall.

Add a gate either way: a `foundation` test that imports the generator and calls
`_build_payload(case)` for every case in `CASES`, asserting no exception. That
makes the harness's own instruction verifiable instead of aspirational.

Also `[R]`: `_generate_bc_equivalence_snapshots.py` annotates
`quad: AngularQuadrature` in three signatures — a name that **does not exist**
(the module imports `Quadrature`). Harmless at runtime under
`from __future__ import annotations`, but it is a dangling reference a `pyright`
pass on `tests/` would flag. And `_build_payload:339` uses a bare `assert` in a
**non-collected** module — vv Mode-8's genuine domain. Both are one-line fixes to
make while the file is open.

---

## 10. (F) The mutation table — every NEW gate, its RED mutation, its positive control

Mutate **in-process by monkeypatch** (never `git checkout` — two tracked files
carry uncommitted state). Every row must be demonstrated under
`.venv/bin/python -O -m pytest`, serial.

| # | Gate | Mutation that MUST turn it RED | Positive control proving the mutation bit |
|---|---|---|---|
| M1 | A1 (T1 forward) | in the albedo arm, pass `axis="y"` to `_specular_kernel` instead of `law.reemission.axis` | the `axis="y"` (LS6 ymax) row stays green ⇒ the mutation is an argument swap, not a body break |
| M2 | A1 (T1) | pass `1.0` instead of `law.albedo` to `_attenuated_kernel_operator` on the albedo route | the α=1.0 rows stay green; α∈{0.0,0.3} red — the signature of an amplitude drop |
| M3 | A1 (capability leg) | force `is_adjointable → False` on `SpecularReemission` and have the realizer honour it | the diffuse rows (already `False`) stay green ⇒ scoped |
| M4 | **A2 (specular anchor)** | `_specular_kernel`: `gamma_out.to_local(perm[inflow])` → `np.arange(inflow.size)` | `[M]` `product(2,4)` and `level_symmetric(6)` rows stay **GREEN** while `gauss_legendre(4/8)` and `lebedev(17)` red — this two-sided outcome IS the proof the fixture choice is load-bearing |
| M5 | A2 | drop the `α` fold (`float(alpha) * base` → `base`) | the α=1.0 row stays green, α=0.3 reds |
| M6 | A3 (diffuse anchor) | in `_checked_angular_average`, pass `-outward_sign` | the orientation guard fires first ⇒ expect a RAISE, not a value red; if it does NOT raise, the guard is the thing that broke |
| M7 | A3 | re-introduce the B3.4a twin: classify Γ₊ with `> 0.0` instead of `> TANGENTIAL_EPS` inside `AngularAverageOperator.from_quadrature` | `[M]` reds on **`product(2,4)` xmax ONLY** (2 extra rows admitted: ordinates 1 and 5); `gauss_legendre(4)`, `lebedev(17)` and `level_symmetric(6)` all stay GREEN (`extra = []`). That asymmetry IS the control, and it is why the diffuse fixture list must contain a `product`-family member — crosswalk §13.2 measured the same thing from the other side |
| M8 | B (G/R invariant) | `AlbedoBoundary.geometry_map` → `SpecularMirror(reemission.axis)` for a specular closure | the isotropic + white rows stay green |
| M9 | C (SN refusal) | delete the `if law.reemission is None: raise` block (fall through to the unknown-closure raise) | the message-substring assertions red while `exc.value.law == "albedo"` still passes — proving the `match` legs, not the bare `raises`, carry the teeth |
| M10 | C (SN refusal, α coverage) | guard the refusal with `if law.reemission is None and law.albedo not in (0.0, 1.0)` | ONLY the α=0.0 and α=1.0 rows red — the fast-path hole, reproduced |
| M11 | CD-2 (diffusion inertness) | diffusion albedo branch dispatches on `type(law.response_kernel)` | CD-1 (bare-law rows) stays green ⇒ the mutation is closure-specific |
| M12 | 9.1 (α=0 zero map) | delete the `if alpha == 0.0` early return in `_attenuated_kernel_operator` | α=0.3 / α=1.0 rows stay green; the four α=0 rows red with `ValueError` |
| M13 | 9.2 (pairing certification) | delete the `isinstance(self.reemission, SpecularReturn)` clause in `albedo.assert_realizable` | the `ReflectiveBoundary` rows still raise ⇒ scoped to the albedo carrier |
| M14 | 9.3 (orientation attribution) | hard-code `law="white"` in `_checked_angular_average`'s raise | the white row stays green; only the albedo row reds |
| M15 | E (domain gate) | revert `_LAWS` to `deferred=True` for the completed albedo rows | rows become `XPASS(strict)` = FAILED — proving the flip is what makes the todo-list mechanism work |
| M16 | D (snapshot inheritance) | scale by `0.7` instead of `0.5` in the albedo route | the sibling `ReflectiveBoundary` snapshot row stays green ⇒ the new row is independently sensitive |

**METHOD WARNING (vv):** before believing any "the gate is blind" verdict,
prove the mutation RAN and changed numbers — print the mutated operator's
`apply` on a fixed probe and confirm it differs from the unmutated one. M4 in
particular has a *legitimately* green half (the two blind quadratures); do not
read that half as "the mutation didn't bite" without checking the other half red.

### 10.1 MEASURED outcomes — run against `tests/geometry/test_reemission_closure.py`

Harness: `/Users/rodrigo/.claude/jobs/c30e4f25/tmp/mut.py`, selected by
`ORPHEUS_B34B=M4|M7|M8|M9|M10|M11|M12|M13`, in-process monkeypatch only.
Unmutated baseline: **`178 passed, 1 xfailed`** (0.4 s, `python -O`, serial).

| # | result | reds | scoped-green control |
|---|---|---|---|
| M4 | ✅ RED | **6** — `gl4_xmax`, `gl8_xmin`, `lebedev17_xmax` × α ∈ {0.3, 1.0}, on the specular ANCHOR | `product24_xmax` + `ls6_ymax` **stay green** (the predicted blindness) — and, crucially, `TestEquivalenceTheorems` is **60 passed / 0 failed** under this mutation |
| M7 | ✅ RED | **6** — all on **`product24_xmax` ONLY** | `gl4`, `gl8`, `lebedev17`, `ls6` stay green — the exact complement of M4 |
| M8 | ✅ RED | **3** — the two `albedo_specular_*` G/R rows + `test_geometry_map_is_the_identity_for_every_closure` | the isotropic + white + bare rows stay green ⇒ scoped |
| M9 | ✅ RED | **3** — `test_refused_at_every_amplitude[0.0/0.5/1.0]` | `exc.value.law == "albedo"` still passes under the mutation ⇒ the reds are attributable to the **message** legs, not to the bare `raises` |
| M10 | ✅ RED | **2** — α=0.0 and α=1.0 ONLY | α=0.5 stays green ⇒ the fast-path hole reproduced exactly |
| M11 | ✅ RED | **6** — only the `isotropic_*` closure rows of CD-2 | the `bare` and `specular_x` rows stay green ⇒ closure-specific |
| M12 | ✅ RED | **8 of 12** in `TestZeroAmplitudeRealizesTheNarrowedZeroMap` | the 4 `negative_control` legs stay green (`ScaledOperator` still refuses) ⇒ the gate measures the realizer's routing, not the numerics layer |
| M13 | ❌ **NOT CAUGHT — 178 passed** | — | **HONEST GAP, see below** |

**⭐ The M4 result is the plan's most important measurement.** Under the
`to_local → arange` mutation, `TestEquivalenceTheorems` is **60 passed, 0
failed** while the operator is wrong on three of five quadratures. That is the
§12.4 disagreement demonstrated rather than argued: the `≡` theorems are
green under a real body bug, and the independent-expression anchor is the only
thing that reds. Do not credit T1/T2 with coverage they cannot provide.

**⚠ M13 is an OPEN gate, not a passed one.** Removing the
`assert_specular_pairing_valid` call from `AlbedoBoundary.assert_realizable`
leaves `test_reemission_closure.py` **entirely green**. That is correct — the
file has no broken-reflection-table fixture — and it means **§9.2's gate is
still owed**. Until it is written, `_specular.py`'s central claim ("one
certification, two laws") is asserted only by its docstring. The gate belongs in
`tests/geometry/test_bc_universal_invariants.py`, which already owns the
ERR-042/044/045 rows; it needs a monkeypatched `Quadrature.reflection_index`
producing three separately-broken tables and must show that BOTH carriers raise,
each with its own `law=` attribution.

**M7's secondary reds are incidental, and the plan says so rather than
claiming them.** M7 also reds two `TestEquivalenceTheorems[*-product24_xmax]`
rows — but by RAISING (the mutated operator expects 4 rows and is fed 2), not by
comparing. Those are Mode-8-class-4-shaped reds *inside the mutation*, not
evidence the ≡ gate caught anything. The attributable catchers are
`test_matches_the_hand_written_cosine_weighted_average` and
`test_the_image_is_isotropic_over_gamma_minus`, both `product24_xmax`-only.
Same caveat applies to M12's α=0 ≡ rows.

---

## 11. (G) Every existing test — MIGRATE / RE-POSE / DELETE

`[M]` The complete inherited blast radius, measured twice (once by simulating the
refusal as a pytest plugin against the pre-carve tree, once against the live
carve): **9 failures across 4 files**, plus 3 xfail rows that decay silently.
`tests/geometry/test_boundary.py` (2 rows) was **not** in the briefed list.

| # | `file:line` | Test | Verdict | One-line reason |
|---|---|---|---|---|
| 1 | `tests/sn/operators/test_sn_boundary_realizer.py` — `TestRealizeAlbedo::test_albedo_zero_realizes_to_zero_operator` | behavioural | **RE-POSE** | the claim "α=0 → `ZeroOperator`" survives, but on a completed law and a face-ful space; fold into the §9.1 α=0 gate |
| 2 | same — `test_albedo_one_realizes_to_identity_operator` | API-shape | **DELETE** | `IdentityOperator` at α=1 is gone by design (the bare TP of the closure's kernel replaces it); the surviving claim is T1/T2 (§4) |
| 3 | same — `test_albedo_half_realizes_to_scaled_tensor_product` | behavioural | **RE-POSE** | "α ∉ {0,1} → `ScaledOperator(α, TP)`" still holds — re-pose onto `AlbedoBoundary(0.5, SpecularReturn("x"))` on a face-ful space |
| 4 | same — `TestRealizeAlbedo` (class) | — | **RENAME + ADD** | becomes the refusal + completion home; add §6.1's refusal rows and §9.1's α=0 rows here |
| 5 | `tests/geometry/test_boundary.py:528-537` — `test_albedo_bc_scales_outgoing` | behavioural | **RE-POSE** | "`α·ψ_out`" is exactly the T1 statement on a slab where the mirror reverses order; re-pose onto `AlbedoBoundary(0.5, SpecularReturn("x"))` + `_realize_narrowed_for_face_right` and assert against the **specular** expectation, not `0.5*psi_out` |
| 6 | `tests/geometry/test_boundary.py:540-570` — `test_albedo_zero_and_vacuum_agree_on_inflow_rows` | behavioural | **RE-POSE** | the α=0 ≡ vacuum claim is now *structurally* true (both route to `_narrowed_zero_operator`); re-pose to compare the two realized operators directly on Γ₊ and assert both are `ZeroOperator` — it becomes a §9.1 row |
| 7 | `tests/geometry/test_boundary.py:54-69` — `_realize_for_sn` helper | — | **RETIRE-OR-NARROW** | its docstring says "for the laws that are STILL full-face — **albedo and periodic only**"; after B3.4b that is periodic only. Update the docstring in the same change (falsified prose = a bug) and plan its deletion at B3.4c |
| 8 | `tests/sn/operators/test_operator_block_role.py:168-170` — `_LINEAR_LAWS` rows `albedo_0/1/03` | behavioural | **MIGRATE** | rewire to completed laws: `AlbedoBoundary(0.0, SpecularReturn(_FACE_AXIS))`, `(1.0, SpecularReturn)`, `(0.3, IsotropicReturn(_FACE_AXIS, _FACE_OUTWARD_SIGN))`. Keep all three α — this file is where the realized **type** claim lives, and α=0 now yields `ZeroOperator` |
| 9 | `tests/geometry/test_bc_equivalence_snapshot.py:207-233` — `TestAlbedo05Lebedev17Snapshot` | characterization of a retired behaviour | **DELETE** | its frozen `psi_in` is `0.5 × psi_out` on the FULL face; no `Γ₊ → Γ₋` map can reproduce the Γ₋ rows (`[M]` direct comparison returns `False`). Successor = §7.1 item 2 |
| 10 | `tests/geometry/snapshots/bc_equivalence_albedo_05_lebedev17.npz` | artefact | **DELETE** | orphaned once (9) goes; `generate_all` would keep re-creating it |
| 11 | `tests/geometry/_generate_bc_equivalence_snapshots.py:166-171` — the `albedo_05_lebedev17` case | — | **DELETE** | pairs with (10) |
| 12 | `tests/geometry/_generate_bc_equivalence_snapshots.py:348-357` — `_build_payload`'s `SNMethodSpace.minimal` | broken infra | **FIX** (§9.4) | `[M]` already raises on 6 of 8 cases; B3.4b makes it 7 of 8 |
| 13 | `tests/sn/operators/test_b3_domain_narrowing.py:286-288` — `_LAWS` albedo rows | **decaying** xfail | **MIGRATE (mandatory)** | `[M]` under the refusal these xfail for the WRONG reason and can never XPASS — Mode-8 class 4. See §8 |
| 14 | `tests/sn/operators/test_b3_domain_narrowing.py:256-270` — `_B34_XFAIL` | reason text | **RE-WRITE** | must lose the albedo clause; periodic's stays for B3.4c |
| 15 | `orpheus/sn/boundary/realizer.py` module docstring bullet for albedo | doc | already **DONE** `[R]` | the main agent rewrote it in this carve |
| 16 | `tests/geometry/test_law_composition.py:177,199` | descriptor algebra only | **NO CHANGE** | `[M]` confirmed green under the refusal — never realizes albedo |
| 17 | `tests/sn/acceleration/test_dsa_low_order.py` | admission guard | **NO CHANGE** | `[M]` confirmed green — refuses before realization |
| 18 | `tests/geometry/{test_bc_universal_invariants,test_bc_errors,test_bound_compat,test_boundary_factors,test_boundary_factor_consumers}.py` | descriptor-level | **NO CHANGE** | `[M]` all confirmed green under the refusal |
| 19 | `tests/diffusion/test_boundary_realizer.py` | diffusion arm | **NO CHANGE** | `[M]` green; CD-1/CD-2 are SHIPPED in the new file (below) |

### 11.1 What is already written, and the ONE Pattern-2 hazard in the split

**SHIPPED** — `/Users/rodrigo/git/nuclear/ORPHEUS/tests/geometry/test_reemission_closure.py`
(net-new, `foundation`, **178 passed / 1 xfailed**, 0.85 s under `python -O`).
It owns the gates with no older home: the ≡ theorems (§4 A1), both
independent-expression anchors (A2/A3), the G/R invariant + its two named
exceptions (§5), the SN refusal message contract (§6.1), the diffusion
CD-1/CD-2 pair (§6.2), and the α=0 narrowed-zero-map family (§9.1).

**Pattern-2 hazard, flagged deliberately:** the migration of `TestRealizeAlbedo`
will naturally want to re-mint a refusal test. **Keep ONE.** The shipped version
is the *message-contract* one (attribution + both completion names + the
"array position" defect, across all three α) — the part most easily
under-gated. If `test_sn_boundary_realizer.py` grows a second, it should assert
something the shipped one does not (e.g. that a FACELESS space refuses for the
*outflow_indices* reason before the closure question is even reached), or be
deleted.

**§9.2 — DELIVERED 2026-08-01.**
`tests/geometry/test_bc_universal_invariants.py::TestSpecularPairingCertifiedOnBothCarriers`
(7 tests). Three broken tables, each wrong in exactly ONE invariant, each
asserted to red on BOTH carriers with its OWN `law=` attribution; plus the
tier leg (albedo's `G` hook stays a no-op even on a catastrophically broken
table — it reds if someone "closes the hole" by overriding the `G` hook), the
diffuse/bare SCOPE control, the positive control over three real quadratures,
and a realize-time seam leg. `[M]` **M13 now reds 5 of the 7**; control 7
passed. Previously M13 left the whole B3.4b gate file green.

`[M]` **One design finding en route, and it changed a row.** The obvious
ERR-044-only mutant — `np.roll(arange(N), 1)` — breaks *all three* invariants,
so ERR-042 fires first and the row proves nothing about independence. Worse,
**no odd cycle can ever be the ERR-044-only mutant**: `a→b→c→a` forces
`sign(a) = −sign(a)`. The construction that works is `π ∘ σ` with `σ` swapping
two SAME-sign, SAME-measure ordinates — and it needs a **degenerate measure
class**, which `gauss_legendre(8)` does not have (its four per-side measures
are all distinct). The quadrature is therefore per-row, with
`level_symmetric(6)` supplying the degeneracy.

### 11.2 Row 12 (the generator) — STOPPED, and why it is not a repair

`[M]` The three SAFE repairs landed: the four `AngularQuadrature` annotations
naming a nonexistent symbol, the bare `assert` in a non-collected module (→
explicit `raise`; vv Mode 8's genuine domain), and the retirement of the
`albedo_05_lebedev17` case.

**The `SNMethodSpace.minimal` breakage is NOT fixable as a fix.** Giving the
generator face-ful spaces gets past the method space and then fails on the
**schema**: a narrowed law consumes `Γ₊` and emits `Γ₋`, so
`op.apply(psi_out)` on the full face is a shape error (`[M]` `white_xmax_LS4`:
`psi.shape[0] = 24, expected |Γ₊| = 12`). The natural next step — store a
`Γ₋`-shaped `psi_in` — would **invalidate every frozen artefact**, because the
harness deliberately consumes the FULL-FACE pre-B3.2 `psi_in` and restricts it
(cases 3–6 each carry a note saying so). That restriction is the artefacts'
entire value: frozen BEFORE the narrowing, they make old-vs-new an independent
statement instead of the code checking itself.

So it is left **broken-and-documented**: a `.. warning::` block in the
generator states the measurement and enumerates the three real options
(reconstruct the retired full-face semantics / migrate the schema and re-anchor
every case against a structurally-independent reference / retire the generator
and declare the artefacts permanently frozen). **That is a user decision.**
The harness's `_GENERATOR_HINT` — which told the reader to run a generator that
raises — was corrected in the same change (falsified prose is a bug).

### 11.3 NET-NEW — the ≡ theorem at the FACTOR-READING consumers

Production gained `law_permutes_ordinates` mid-run, because the ≡ theorem was
**false at four factor-reading sites** (sweep schedule, DSA admission,
`_has_ruled_corner_action`, `_reflect_corner`). `TestEquivalenceTheorems` is
structurally blind to that class — it reads realized OUTPUT — so this needed
its own gate, and unlike A1 **it has teeth: the two routes share no body
here.** `TestEquivalenceHoldsAtTheFactorReadingConsumers` (5 tests): the
equivalence itself (stated as an equivalence, not a hard-coded table, so it
cannot drift with the design), a non-vacuity leg pinning BOTH answers
(`False`-for-everything is precisely the pre-B3.4b bug's shape), a scope
control over the four non-pairing laws, and a leg asserting the retired
single-tier read **diverges** — so a "simplification" back to
`geometry_map.permutes_ordinates` reds. New mutation **M17** does exactly that:
`[M]` 3 red, the isotropic row green (scoped).

### 11.4 NET-NEW — the refusal message is keyed to α

`[M]` The α=0 message repair landed mid-run. The refusal now says two
different things for two different reasons, and the gate says so: at α ≠ 0 it
must name **both completions** and the **array-position** defect; at α = 0 it
must name the **`VacuumInflow` twin** and must **NOT** claim
under-determination (nothing returns, so no pairing is missing — reciting the
α≠0 reason there sends the reader after a defect that is not present). The α=0
leg carries a CONTROL asserting α=0.5 *does* cite under-determination, so the
negative is discriminating rather than vacuously true of every message.

**Nothing in this list is a pure delete-and-lose-coverage.** Rows 2 and 9 are the
only deletions of *claims*, and both are claims the phase deliberately retires
(an `IdentityOperator` at α=1, and a full-face endomorphism's frozen image);
their behavioural content is carried forward by the T1/T2 gates and the
snapshot-inheritance row.

---

## 12. Where I disagree with §12's design

### 12.1 The refusal criterion and the exactly-one invariant are DIFFERENT statements

§12.3 argues the refusal from "`α·I` is a `Γ₊ → Γ₊` endomorphism and
`IdentityMap` supplies no crossing", and §11.1 states "EXACTLY ONE of `G`, `R` is
non-trivial" as the invariant that "wants a gate". Read together they suggest one
criterion. They are not.

`[M]` `AlbedoBoundary(0.5)` bare **satisfies** the exactly-one invariant
(`G = I` trivial, `R = 0.5·I` non-trivial) and is **refused**. `[M]`
`AlbedoBoundary(1.0)` bare **violates** it (zero non-trivial factors) and is
refused for the same reason as the 0.5 case. So the count does not decide the
refusal.

The decidable refusal predicate is about `R`'s **type**, not the count:

```
under-determined on an angular trace  ⟺  G supplies no crossing
                                          AND R is a non-zero SCALAR amplitude
                                      ⟺  isinstance(G, IdentityMap)
                                          and isinstance(R, ScalarResponse)
                                          and not R.is_zero
```

`[D]` Ask: write both down separately in the docs and gate them separately (§5
and §6). Presenting them as one statement will produce a future "simplification"
that merges them and breaks one.

### 12.2 The α = 0 edge of the refusal is undecided — and I think the carve is right for an unwritten reason

By the predicate above, `AlbedoBoundary(0.0)` bare is **complete** on an angular
trace: nothing returns, so no pairing is needed — the zero map is fully
determined. Yet the carve refuses it. `[R]` §12.3 says "SN refuses it" flat and
never addresses α=0.

I think refusing is right, for two reasons §12 does not give:

1. **Continuity.** `ScalarResponse.is_zero` is documented as an *exact* compare
   with the rationale *"a near-zero albedo is a weak reflector, not a vacuum"*
   (`_factors.py:450-454`). If a weak reflector needs a closure, so does its
   α → 0 limit; otherwise realizability is discontinuous in a parameter at
   exactly one float value.
2. **Twin path.** `VacuumInflow` already IS the α=0 law. Letting a bare
   `AlbedoBoundary(0.0)` realize keeps a second spelling of vacuum alive on the
   SN arm — Pattern 2, in the very campaign that exists to remove twins.

`[D]` Write that reasoning into the refusal's docstring, and keep the α=0 row in
the refusal gate (§6.1) so the decision is pinned either way. Without it, the
next reader sees an "obviously wrong" refusal of a well-determined law and
"fixes" it.

### 12.3 §12.4 under-counts the invariant's violators

Detailed in §0.2 / §5.2. `AlbedoBoundary(1.0)` bare carries **zero** non-trivial
factors and is not mentioned. It is not a defect — it is the degenerate
factorisation a scalar trace forces — but a gate posed over the registry will hit
it, and it must be a named, asserted exception rather than a surprise.

### 12.4 §12.2's "the two routes execute the same body" argument is stated as if it were the verification

§12.2: *"The `≡` theorems then hold because the two routes literally execute the
same code, not because two transcriptions agree — and they get pinned as tests
anyway, since a shared body can still be reached with different arguments."*

The second clause is right and the plan should lead with it, because the first
clause is exactly what makes the `≡` gate **structurally non-independent**: a
shared body carries a shared bug, so T1/T2 are green under *any* body error.
`[M]` Demonstrated: mutation M4 (the `arange` remap) leaves T1 **green on every
quadrature** while the operator is wrong on three of five. The `≡` gates are
necessary-not-sufficient and the independent-expression anchors A2/A3 are the
actual catchers. Say so in the commit message, or a future reader will credit
T1/T2 with coverage they cannot provide.

### 12.5 The specular closure's orientation guard is weaker and less attributable than the diffuse one

`[M]` On `level_symmetric(6)` `xmax`:

- `AlbedoBoundary(0.5, IsotropicReturn("y", +1))` → attributed `BoundaryError`
  with `law="albedo"` and a message explaining the hemisphere mismatch;
- `AlbedoBoundary(0.5, SpecularReturn("y"))` → raw
  `ValueError: TraceRestrictionOperator.to_local: row(s) [...]` — a numerics-layer
  error with **no `law=` attribution** and no explanation that the closure's axis
  disagrees with the installation face.

Both are refusals, so nothing is silently wrong. But two sibling routes with two
error contracts is the kind of asymmetry that becomes a bug the day someone
catches `BoundaryError` around a realization loop. `[D]` Either wrap the
`to_local` failure in an attributed `BoundaryError` in `_specular_kernel`, or
give the specular route the same declared-vs-installed cross-check the diffuse
route has. I lean **the cross-check**, because it is the ERR-041 pattern and it
would also catch the case where the wrong axis happens *not* to cross index sets.
Gate it either way (§9.3), and if the decision is "leave the `ValueError`", the
gate should assert the `ValueError` explicitly so the asymmetry is *documented*
rather than accidental.

---

## 13. Run order, and what "done" means

1. Write the new/changed tests (§4, §5, §6, §7, §8, §9). They red immediately
   where the production is not yet right — that is the point.
2. Run, serially, under the canonical invocation:
   ```
   .venv/bin/python -O -m pytest tests/geometry tests/sn/operators tests/diffusion -q -p no:randomly
   ```
   `[M]` **Pre-migration baseline on the live carve WITH
   `tests/geometry/test_reemission_closure.py` present: `9 failed, 1489 passed,
   5 skipped, 13 xfailed`** — the nine rows of §11 and nothing else; the new file
   contributes **178 passed, 1 xfailed** and zero collateral. Target: **0
   failed**, with the 3 albedo xfails GONE (not converted to passes-by-xfail).

   ⚠ A first run of this command during the carve reported `32 failed`,
   including `test_psi_half_coupling.py` and `test_typed_residual_evaluation.py`
   — a TRANSIENT caused by reading `orpheus/sn/operators/boundary.py` mid-edit.
   It reproduced as 9 once the tree settled. **Do not baseline against a tree
   being written**; re-run before believing any count.
3. `--runxfail` and READ every remaining xfail reason — confirm each reds for ITS
   documented reason (this is what caught the §8.1 misattribution).
4. Run the §10 mutation table. Every row must be demonstrated RED with its
   positive control GREEN, in-process by monkeypatch.
5. Wide gate:
   ```
   .venv/bin/python -O -m pytest tests/{geometry,sn,diffusion,numerics,transport} -m "not slow" -q -p no:randomly
   ```
   `[M]` Reference point on the pre-carve tree with the refusal simulated:
   `9 failed, 4342 passed, 6 skipped, 111 deselected, 67 xfailed` in 11m27s
   (serial, contended host). Target: `0 failed`, xfail count down by 3.
6. `npx pyright orpheus/` = 1 (the #288 floor); `sphinx -E -W` exit 0.

**A gate is not done when it passes. It is done when its §10 mutation has been
run and observed RED under `python -O`.**
