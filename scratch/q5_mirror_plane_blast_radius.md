# Q5 — `Z2` → parameterized mirror: blast-radius map

Audited 2026-08-02 on `refactor/operator-strategy-layers` @ `b54c1066`
(working tree clean except `.claude/skills/vv-principles/*` — **the carve is
NOT underway**; `orpheus/numerics/symmetry.py` and
`orpheus/numerics/quadrature/registry.py` are both clean at HEAD, git-verified
via `git status --short` + `git diff --stat`).

Line numbers current at that HEAD; re-derive via Nexus if drifted.
Probes: `scratch/_q5_probe.py`, `scratch/_q5_probe2.py`.

> **DRIFT NOTICE (L-007/L-012).** `orpheus/numerics/quadrature/directional.py`
> went clean → **modified (+68/−41) DURING this audit** (re-checked at close).
> Every `directional.py:NNN` citation below is **working-tree**, NOT at HEAD.
> Two things the in-flight diff does, both load-bearing for this question:
> (i) it **retires `_find_reflections`** (the guard-free O(N²) `argmin`, HEAD
> `:125`) and (ii) it rewrites `_compute_sphere_reflection_partners` to take a
> `DiscreteMeasure` and route through `_orbit_closure` — i.e. **the certified
> per-axis mirror check in §7.2 is UNCOMMITTED work-in-progress, not landed.**
> `symmetry.py` and `registry.py` are clean at HEAD; all their citations are
> HEAD-true.

---

## 0. The one-paragraph shape

`_NamedSubgroup.Z2` is a **named, parameter-free** enum entry realized as
**exactly one matrix**, `diag(1, 1, −1)` = σ_z (`symmetry.py:510-511` →
`_reflections("z")` at `:888-899`; measured, §3.0 below). Three machines consume
the tag and they do **not** all reach that realization:

| machine | entry point | reaches σ_z? |
|---|---|---|
| **containment** (`contains`) | `_contains` `symmetry.py:608` → `_finite_contains` `:598` | **YES** — literal matrix-set containment over `_close_group(_realized_ops(tag))` |
| **invariance, 3-D nodes** | `_check_invariance_3d` `:815` | **YES** — `_orbit_closure(..., _reflections("z"), ...)` |
| **invariance, 1-D nodes** | `_check_invariance_1d` `:730`, `Z2` arm `:759-767` | **NO** — short-circuits to `_is_reflection_invariant_1d` (`x → −x` on the single coordinate). `_realized_ops` is **never called** (measured: 0 calls) |

So `Z2` means **two different things** depending on the node shape, and the 1-D
meaning ("`x → −x` on whatever the one coordinate is") is **axis-free**, while
the 3-D meaning is **σ_z specifically**. The geometry table's slab/sphere entries
hit only the 1-D meaning.

---

## 1. Every `Z2` use site

### 1a. Production (`orpheus/`) — 4 sites, 2 files

| `file:line` | what the caller MEANS | breaks under parameterization? |
|---|---|---|
| `orpheus/numerics/symmetry.py:109` | **the enum member itself.** Comment `{e, σ} — single reflection / 180° rotation` is the ambiguity: it reads as "any order-2 group" but the realization is σ_z | **it IS the change** |
| `orpheus/numerics/symmetry.py:205,221-224` | **lattice edges.** `(Trivial,Z2)`, `(Z2,Dinfh)`, `(Z2,Oh)`, `(Z2,Ih)`, `(Z2,O3)`. Semantically "the improper order-2 group sits under everything except the proper-rotation groups" | **3 of 5 edges are ALREADY DEAD CODE** — see §5 |
| `orpheus/numerics/symmetry.py:510-511` | **the realization.** `Z2 → _reflections("z")` = σ_z | the axis moves from a literal to a field |
| `orpheus/numerics/symmetry.py:759-767` | **1-D arm.** `x → −x` on the single coordinate — plane-free, "the inversion-like closure" of the marginal | **the load-bearing semantic mismatch** (see §3) |
| `orpheus/numerics/symmetry.py:815-819` | **3-D arm.** σ_z closure, docstring says *"pick z → −z as a canonical representative. (Any single reflection works; the choice is convention.)"* — **that comment is false**, measured §3.6/§3.7 | the "convention" claim must go |
| `orpheus/numerics/symmetry.py:1278` | singleton install | needs a sibling install or a constructor |
| `orpheus/numerics/symmetry.py:1336` | `candidate_groups` named list — `Z2` is offered to the lattice walk for **every** measure | a family needs an axis enumeration here (like `Cn`/`Dnh` get divisors at `:1345-1347`) |
| `orpheus/numerics/quadrature/registry.py:671` (`slab`) | **`μ → −μ`** — "the forward/backward reflection that pairs the two sweep senses and that a reflecting end face consumes" (`registry.py:655-658`). NOT σ_z, NOT "any order-2 group" | **YES — this is the conflict** (§3) |
| `orpheus/numerics/quadrature/registry.py:675` (`sphere`) | same: `μ_r → −μ_r` on the radial cosine | **YES** |

Also relevant, **not** a `Z2` site but the same concept unnamed:

| `file:line` | what it is |
|---|---|
| `directional.py:157-215` `_compute_sphere_reflection_partners` *(UNCOMMITTED — see drift notice)* | builds a **per-axis** mirror certificate inline: `reflection = eye(3); reflection[axis,axis] = -1` for `axis in range(3)`, then `_orbit_closure`. This IS `Mirror(axis).is_invariant` + its `OrbitCertificate`, spelled by hand. **The latent consumer.** |
| `directional.py:400-422` `reflection_index(axis: int \| str)` | the mirror plane is ALREADY a runtime parameter at this layer — `0/1/2` or `"x"/"y"/"z"` |
| `orpheus/geometry/boundary/_specular.py:103-111` | `_axis_cosines(quad, axis)` → `AXIS_NAMES.index(axis)` — specular BCs are already per-plane |

### 1b. Tests — 21 sites, 3 files

| `file:line` | what it means | breaks under parameterization? |
|---|---|---|
| `tests/numerics/test_symmetry.py:67` | `_NAMED` list — "the 8 named entries" | **YES if `Z2` leaves the enum** (list is parametrized into 3+ tests) |
| `test_symmetry.py:96` | `O3.contains(G)` for every named G | needs `O3.contains(Mirror(a))` — **measured FALSE today for a new tag**, §5D |
| `test_symmetry.py:121-123` | `test_z2_chain_via_octahedral`: `Oh ⊃ Z2 ⊂ O3` — "the abstract improper order-2 group" | `Oh ⊃ Mirror(any)` computes True; `O3 ⊃ Mirror` needs a new arm |
| `test_symmetry.py:208` `test_dnh_reflection_in_dnh` | `Dnh(n).contains(Z2)` for n ∈ {1,2,3,4,6} — docstring says *"`Z_2` (a single reflection) sits inside every `D_nh`"* | **THE SHARPEST BREAK.** True only for σ_z. Measured: `Dnh(1).contains(Mirror('x')) = False`, `Dnh(3).contains(Mirror('x')) = False` (§5C). The test's own prose ("a single reflection") is the false generalization |
| `test_symmetry.py:240` | `Z2.is_subgroup_of(Oh)` readability synonym | computed, survives |
| `test_symmetry.py:329` `test_gauss_legendre_1d_z2_reflective` | **`x → −x` on a 1-D GL measure** — explicitly the polar mirror, plane-free | survives *iff* the new member routes to the 1-D branch |
| `test_symmetry.py:378,383` | repr/equality smoke: `Z2 == Z2`, `"Z2" in repr(Z2)` | **repr and `.name` must be extended** (§6) |
| `test_symmetry.py:679` | `test_realized_group_orders_are_textbook`: `"Z2": (Z2, 2)` — group ORDER 2 | survives (one generator ⇒ order 2) |
| `test_symmetry.py:730` | `test_computed_containment_obeys_the_order_relation_laws` — `Z2` in the `finite` list; calls `_close_group(_realized_ops(g._tag))` | survives; **and this is the test that will catch a wrong new lattice arm** |
| `test_symmetry.py:1160-1175` `test_z2_is_improper_so_it_is_not_inside_the_rotation_groups` | **`det = −1`, so `Z2 ⊄ SO(2)/SO(3)`, `Z2 ⊂ O(3)`.** Prose: *"realized as σ_z"* — the ONE test that already names the plane | survives for the realization half; the `O3.contains` half needs the new arm |
| `test_symmetry.py:1221` `test_invariance_is_downward_closed_on_polar_marginals` | monotonicity on the **1-D path**; docstring: *"the sphere-side gate cannot see 1-D bugs: that path never calls `_orbit_closure` and has its own dispatch"* | **the project's own written admission of the split in §0.** A Mirror family must be added to `tags` here or the new member is untested on the 1-D path |
| `tests/numerics/test_registry.py:105,109` | pins `GEOMETRY_ANGULAR_SYMMETRY` by value equality — slab/sphere ⇒ `Z2` | **YES — a table-value change breaks this by construction** |
| `test_registry.py:582,595` | `AngularSymmetry(SO2, Z2).support == "[-1,1]"`; the `NotImplementedError` guard | survives (only `continuous_isotropy` is read) |
| `test_registry.py:668` | `not slab.discrete_residual.is_subgroup_of(gl.invariance_group)` — i.e. `not SO2.contains(Z2)`, the "declared tag is a trap" demo | **fragile**: needs `not SO2.contains(Mirror(a))`, which is accidentally True today only via the `return False` fallthrough (§5D) |
| `tests/numerics/test_measure.py:596,599` | `"SubgroupOfO3.Z2" in repr(...)` — **a literal string assertion on repr** | **YES if the repr text changes** |
| `test_measure.py:661,666` | `invariance_group=Z2` round-trips through `consolidate()` | survives (opaque value) |

### 1c. Docs — `docs/theory/foundations/discrete_measures.rst`

| `file:line` | what it says | status |
|---|---|---|
| `:425` | slab `→ SO(2) × Z_2` | consistent with the table; would need the plane named |
| `:440` | `Trivial ⊂ Z_2 ⊂ O_h ⊂ O(3)` | fine |
| `:461` | *"the eight named entries … `Trivial`, `Z_2`, `SO(2)`, **`O(2)`**, `O_h`, `I_h`, `SO(3)`, `O(3)`"* | **ALREADY STALE** — `O2` was renamed `Dinfh` on 2026-08-02 (`symmetry.py:117-129`). Pre-existing drift, not caused by this carve |
| `:489` | *"the `Z_2` reflection `x → −x` is the only non-trivial 1-D check"* | correct, and it is the **plane-free** reading |
| `:629,:638` | the geometry table rows: slab/sphere owe `Z_2`, prose says **`μ → −μ`** | must be re-spelled with the plane |

No `:func:`/`:class:` Python-domain xref to `Z2` exists (it is an enum member, referenced only as inline literal), so the Sphinx `-W` gate will **not** catch a doc/code divergence here. Grep `discrete_measures.rst` by hand.

---

## 2. The 1-D path (question 2)

`symmetry.py:706-713` picks the branch:

```python
is_1d = measure.dim == 1 or measure.support in (SPACE_INTERVAL_M11, "[0,1]", "R")
if is_1d:
    return _check_invariance_1d(tag, measure, atol)
```

For a slab GL measure (`nodes.shape == (8,)`, `dim == 1`, `support == "[-1,1]"`)
the 1-D branch is taken. The `Z2` arm (`:759-767`) calls
`_is_reflection_invariant_1d` (`:775-804`), which:

- takes `nodes[:, 0]` if 2-D else the flat array (`:777`),
- for each node `x`, binary-searches for `−x` (`:787-790`),
- requires a partner within `atol*10` **and** matching weight within `atol`.

So **yes, it is literally `μ → −μ`** on whichever single coordinate the measure
carries — and **no, it never touches `_realized_ops` or `_orbit_closure`.**

Measured (probe §3):

```
slab GL(8) + Z2 -> _realized_ops calls = 0   _orbit_closure calls = 0
lebedev(9)  + Z2 -> _realized_ops calls = 0   _orbit_closure calls = 1
```

(`_realized_ops` reads 0 for Lebedev too because `_check_invariance_3d` inlines
`_reflections("z")` at `:819` rather than routing through `_realized_ops` — a
second, independent spelling of the same realization.)

The 1-D branch is **not vacuous** — an asymmetric μ-set is correctly rejected:

```
Z2.is_invariant([0.1, 0.5, 0.9]) = False
```

Note the arm is a **set membership over 6 named tags** (`:759-766`) plus a
`Dnh` isinstance (`:769`); a new parameterized mirror type gets **`False` by
fallthrough at `:772`** unless explicitly added.

---

## 3. The geometry-table conflict (question 3) — MEASURED

### The answer: **(c) — something else, and it is more interesting than (a) or (b).**

**Today it is (b)-harmless, for TWO independent reasons — and (a) is a real
latent inconsistency that is one embedding-change away from firing.**

### 3.0 What `Z2` actually is

```
_realized_ops(Z2) = [[1,0,0],[0,1,0],[0,0,-1]]     -> diag [ 1.  1. -1.]  == sigma_z
_reflections('x') = [[-1,0,0],[0,1,0],[0,0,1]]
_reflections('y') = [[1,0,0],[0,-1,0],[0,0,1]]
```

### 3.1 Why it is harmless TODAY — reason 1: the 1-D branch

Slab and sphere both declare `continuous_isotropy = SO2`, so
`AngularSymmetry.support` derives `"[-1,1]"` (`registry.py:626-631`). Their
rules are 1-D, so `Z2` routes to `_check_invariance_1d` and the σ_z matrix is
never built. Measured:

```
slab         G0=SO2      Gamma=Z2       support='[-1,1]'
sphere       G0=SO2      Gamma=Z2       support='[-1,1]'
cylinder     G0=Trivial  Gamma=D_2h     support='S^2'
cartesian2d  G0=Trivial  Gamma=D_2h     support='S^2'

SLAB rule GaussLegendre1D:  nodes.shape=(8,)  dim=1  support='[-1,1]'
   Z2.is_invariant  = True
   -> _check_invariance_1d(Z2, mu, 1e-13) directly: True
```

### 3.2 Why it is harmless TODAY — reason 2: stage 0 fires first

`select_quadrature` runs **stage 0 (domain) before stage 1 (symmetry)**
(`registry.py:882-908`). Every S² rule is rejected for slab/sphere before `Z2`
is ever asked:

```
GaussLegendre1D      support='[-1,1]'  stage0(domain)=True  stage1(Z2)=True
LebedevSphere        support='S^2'     stage0(domain)=False  -> REJECTED before Z2 is ever asked
LevelSymmetricSN     support='S^2'     stage0(domain)=False  -> REJECTED before Z2 is ever asked
ProductQuadrature    support='S^2'     stage0(domain)=False  -> REJECTED before Z2 is ever asked
```

So **no production path can currently reach "slab μ + σ_z".**

### 3.3 A SYMMETRIC μ-set cannot see the difference at all

This is why the inconsistency has never surfaced. Embedding a **symmetric** GL
μ-set into 3-D leaves two all-zero columns, which every coordinate mirror fixes:

```
(mu,0,0)  X-embed: sigma_x=True  sigma_y=True  sigma_z=True   -> Z2 (3-D route) = True
(0,0,mu)  Z-embed: sigma_x=True  sigma_y=True  sigma_z=True   -> Z2 (3-D route) = True
```

The σ_x/σ_z distinction is **invisible on any μ-set that already satisfies the
owed residual**. The gate agrees with itself only because the answer is True
either way.

### 3.4 The DISCRIMINATING measurement — an ASYMMETRIC μ-set

Break `μ → −μ` and the two embeddings diverge:

```
asym (mu,0,0) X-embed: sigma_x=False sigma_y=True  sigma_z=True    Z2(3-D)=True
asym (0,0,mu) Z-embed: sigma_x=True  sigma_y=True  sigma_z=False   Z2(3-D)=False
```

**`Z2(3-D) = True` on a μ-set that violates `μ → −μ`.** That is a **false
certification in the dangerous direction** (the ERR-072 family: a gate that
CERTIFIES a non-invariant rule), and it fires the moment a slab μ reaches the
3-D branch in the tree's own canonical embedding.

### 3.5 The embedding really is `(μ, 0, 0)`

Not an inference — it is written down twice:

- `directional.py:300-322` `axis_cosines` (`:318` is the `nodes.ndim == 1`
  branch): index 0 returns the nodes, indices 1/2 return zeros.
- `directional.py:444-452` `angular_frame` docstring: *"a slab's polar `μ`
  embeds as `(μ, 0, 0)` — the SAME column-stacked embedding
  `spherical_harmonics` uses internally"*.
- `directional.py:428-442` `spherical_harmonics` feeds
  `axis_cosines(0), axis_cosines(1), axis_cosines(2)`.

(All three re-derived against the **working tree** after the mid-audit edit;
`axis_cosines` was at `:273` at HEAD.)

So the slab mirror **is σ_x** in the tree's own 3-D embedding, and `Z2` **is
σ_z**. They are different matrices; they agree only on inputs where both are
True.

### 3.6 The gap is LIVE on a shipped rule (no asymmetric fixture needed)

```
lebedev(9)           Z2.is_invariant=True   s_x=True  s_y=True  s_z=True
level_symmetric(8)   Z2.is_invariant=True   s_x=True  s_y=True  s_z=True
product(4,8)         Z2.is_invariant=True   s_x=True  s_y=True  s_z=True
product(4,3)         Z2.is_invariant=True   s_x=False s_y=True  s_z=True   <-- !!
```

**`product(4, 3)` is NOT closed under σ_x but IS closed under σ_z, and `Z2`
answers `True`.** An odd-`n_phi` product rule already discriminates the two
planes on the tree today. (This is the same odd-`n_φ` defect
`_compute_sphere_reflection_partners` is *currently being* hardened against —
`directional.py:173-181`, **uncommitted**, see drift notice — so the tree is
paying for this exact distinction one layer up, right now, without naming it.)

### 3.7 The docstring claim at `:816-818` is false

> *"Single reflection — pick z → −z as a canonical representative. (Any single
> reflection works; the choice is convention.)"*

Falsified by §3.4 and §3.6: `product(4,3)` and any asymmetric μ-set give
different answers for σ_x and σ_z. The choice is **not** convention.

### 3.8 Verdict

- **Not (a) as-shipped**: no production route reaches the wrong matrix (two
  independent guards: 1-D dispatch + stage-0 ordering).
- **Not (b) either**: "harmless because the 1-D branch never reaches
  `_realized_ops`" is true but is an *accident of dispatch*, not a designed
  invariant. Nothing in the code says "slab μ must not reach the 3-D branch";
  `_check_invariance` `:706` decides it by a support-string membership test.
- **(c)**: the tag is **overloaded across two branches with incompatible
  semantics** (plane-free vs σ_z-specific), the overload is masked by a
  degenerate fixture (symmetric μ ⇒ all planes agree), and it is already
  *observably* wrong on a shipped S² rule (`product(4,3)`, §3.6). Naming the
  plane is what makes the two branches say the same thing.

---

## 4. `_reflections(axis)` and `_vertical_mirrors(n)` (question 4)

### `_reflections(axis)` — `symmetry.py:888-899`

```python
def _reflections(axis: str) -> list[np.ndarray]:
    M = np.eye(3)
    if axis == "x":   M[0, 0] = -1.0
    elif axis == "y": M[1, 1] = -1.0
    elif axis == "z": M[2, 2] = -1.0
    else: raise ValueError(f"axis must be x/y/z, got {axis!r}")
    return [M]
```

**Yes — the plumbing already exists.** It takes `axis: str`, accepts exactly
`"x"` / `"y"` / `"z"`, raises on anything else, and returns a **one-element
list** (a generating set of one). Every one of the three axes is exercised by
the probe. `σ_y` is already constructible — only the *enum member naming it* is
missing. The name is plural (`_reflections`) but the return is always a
singleton; that is the residue of it being written for the `Dnh` composition at
`:520` / `:863`.

Callers: `symmetry.py:511` (`Z2`), `:520` (`Dnh` realization), `:819` (`Z2`
3-D arm), `:836` (`Dinfh` σ_h conjunct), `:863` (`Dnh` 3-D arm). All five pass
the literal `"z"`.

### `_vertical_mirrors(n)` — `symmetry.py:958-990`

Returns the `n` vertical mirror planes of `D_nh` (planes containing z). The
k-th **plane** is at `kπ/n` from the x-axis, so the **normal** is at
`kπ/n + π/2` (`:986`). Docstring records a fixed bug: placing the *normal* at
`kπ/n` rotated every plane by π/2 — invisible for even `n`, wrong for odd `n`,
and it made `Dnh(n).is_invariant(product(4,n))` read False at n = 3,5,7.

**This is directly relevant to the mirror-plane question**: which planes `D_nh`
contains is `n`-dependent, so `Dnh(n) ⊇ Mirror(axis)` is a genuinely
`(n, axis)`-dependent fact — measured in §5C, and it is exactly what breaks
`test_dnh_reflection_in_dnh`.

---

## 5. The subgroup lattice (question 5) — MEASURED

### 5A. Three of the five `Z2` edges are ALREADY DEAD CODE

`_contains` (`:608-622`) decides **finite × finite by computed matrix
containment** and only falls back to `_NAMED_LATTICE` when one side is
continuous. Measured, per edge:

```
(Trivial        subset Z2            )  finite=True/True   -> COMPUTED (table entry is DEAD)
(Z2             subset IcosahedralIh )  finite=True/True   -> COMPUTED (table entry is DEAD)
(Z2             subset OctahedralOh  )  finite=True/True   -> COMPUTED (table entry is DEAD)
(Z2             subset Dinfh         )  finite=True/False  -> table-read (LIVE)
(Z2             subset O3            )  finite=True/False  -> table-read (LIVE)
```

**So a parameterized mirror family needs only TWO lattice facts, not five:**
`Mirror(a) ⊆ D_∞h` and `Mirror(a) ⊆ O(3)`. Everything else against a finite
group is computed for free. (`Mirror(a) ⊆ D_∞h` is true for all three axes:
D_∞h contains σ_h *and* every σ_v.)

### 5B. But `_NAMED_LATTICE` is `set[tuple[_NamedSubgroup, _NamedSubgroup]]`

`:203` — the type is **enum-to-enum**. A parameterized tag cannot be an
element. `Cn`/`Dnh` are handled *outside* the table entirely, in the
`isinstance`-dispatch arms of `_contains` (`:628-686`). **That is the template.**

### 5C. Measured: what a live `Mirror(axis)` prototype gets for FREE

Monkeypatched a `@dataclass(frozen=True) class Mirror: axis: str` with only
`_realized_ops` + the two invariance arms extended (no lattice edits at all):

```
Oh.contains(Mirror('x'))   = True
Oh.contains(Mirror('y'))   = True
Oh.contains(Mirror('z'))   = True
Dnh(1).contains(Mirror('x')) = False      Dnh(1).contains(Mirror('z')) = True
Dnh(2).contains(Mirror('x')) = True       Dnh(2).contains(Mirror('z')) = True
Dnh(3).contains(Mirror('x')) = False      Dnh(3).contains(Mirror('z')) = True
Dnh(4).contains(Mirror('x')) = True       Dnh(4).contains(Mirror('z')) = True
Mirror('z').contains(Trivial)   = True
Mirror('x').contains(Mirror('z')) = False
Z2.contains(Mirror('z'))          = True     <-- same matrix set, computed correctly
```

Every finite × finite relation is **already right with zero lattice work**,
including the `n`-dependent `Dnh` answers §4 predicts.

Note `Z2.contains(Mirror('z')) == True` while `Mirror('z') == Z2` is `False`
(§6) — the group-theoretic identity is computed correctly but the value-object
identity is not. Two spellings of one group.

### 5D. What a new tag gets WRONG, silently

```
O3      .contains(Mirror('x')) = False     (Z2 answer: True)   <-- WRONG
Dinfh   .contains(Mirror('x')) = False     (Z2 answer: True)   <-- WRONG
SO3     .contains(Mirror('x')) = False     (Z2 answer: False)  -- right by luck
SO2     .contains(Mirror('x')) = False     (Z2 answer: False)  -- right by luck
```

`_contains` ends with a bare `return False` at `:686`. **A new tag type gets a
wrong-but-silent answer, not an error.** `O3 ⊉ Mirror` is flatly false and would
break the walk's soundness precondition (`maximal_invariance_groups` docstring
`:1390-1395`: *"the walk therefore requires a correct lattice"*).

Two of the four are "right by luck" — SO2/SO3 are proper-rotation groups so
`False` is the correct answer, reached by the wrong route.

### 5E. Does the walk special-case parameterized types? — NO, and that IS the template

`maximal_invariance_groups` (`:1373-1440`) never mentions `Cn`/`Dnh`. It:

1. calls `candidate_groups(measure)` (`:1326-1348`) — **the only place the
   parameter is enumerated.** Named entries always; `Cn(d)`/`Dnh(d)` for each
   **divisor** of `_distinct_azimuths(nodes)`. That bound (`:1303-1323`) is what
   turns two infinite families into a finite candidate set.
2. seeds `stack = _maximal(cands)` and descends via `maximal_subgroups`, which
   is *derived from* `contains` (`:1360-1370`), never tabulated.
3. keys `visited` on **`repr(group)`** (`:1430`) and `_GROUP_CACHE` on
   **`repr(tag)`** (`:586`).

The `repr` keying is a live hazard: the fallback branch `f"SubgroupOfO3({tag!r})"`
(`:346`) happens to distinguish `Mirror('x')/('y')/('z')` because the dataclass
repr carries the field — measured `repr collision = False`. But that safety is
**accidental**; any custom `__repr__` that drops the axis silently collapses the
three mirrors into one cache/visited key.

---

## 6. The `Cn`/`Dnh` template — exactly what a new parameterized member implements

Checklist derived from `grep -n "Cn\b\|Dnh\b" orpheus/numerics/symmetry.py`
(48 hits, deduped to 11 structural obligations). **`Cn`/`Dnh` are frozen
dataclasses, NOT enum members**, and they are **never imported outside
`symmetry.py`** (verified: every consumer goes through
`SubgroupOfO3.Cn(n)` / `.Dnh(n)`), so a new family is fully internal.

| # | obligation | `Cn` at | `Dnh` at | notes for a `Mirror(axis)` |
|---|---|---|---|---|
| 1 | **frozen dataclass** with `__post_init__` validation | `:132-151` | `:154-176` | `axis: str` (or an `Axis` enum) validated to x/y/z. `_reflections` already raises on bad input — but validate at construction, not at realization |
| 2 | join the `SubgroupTag` alias | `:182` | `:182` | `SubgroupTag = "_NamedSubgroup \| Cn \| Dnh \| Mirror"` (a plain string alias today) |
| 3 | `__init__` type annotation | `:309` | `:309` | same union |
| 4 | **classmethod constructor** on the façade | `:318-321` | `:323-326` | `SubgroupOfO3.Mirror(axis)`. Note the classmethod **shadows the dataclass name** inside the class body — that is the existing idiom |
| 5 | **`__repr__` arm** | `:342-343` | `:344-345` | REQUIRED — the fallback at `:346` works but is ugly (`SubgroupOfO3(Mirror(axis='x'))`) and is what `visited`/`_GROUP_CACHE` key on |
| 6 | **`.name` arm** | `:354-355` | `:356-357` | REQUIRED — `.name` is **user-facing**: `registry.py:906` selector-rejection message, `symmetry.py:1500-1503` `singular_set` ValueError, `registry.py:635` `NotImplementedError`. Un-extended it emits `Mirror(axis='x')`. Suggest `sigma_x` / `σ_x` |
| 7 | **`_realized_ops` arm** | `:517-518` | `:519-520` | `return _reflections(tag.axis)` — one line, plumbing exists |
| 8 | **`_contains` arms — only the MIXED finite/continuous ones** | `:629-634`, `:646-661`, `:669-681` | `:636-643`, `:663-666`, `:672-685` | finite × finite is computed (§5C). Needed: `outer ∈ {Dinfh, O3} ⇒ True`, `outer ∈ {SO2, SO3} ⇒ False`, and `Mirror ⊇ named` ⇒ `inner is Trivial`. **Mind the bare `return False` at `:686`** |
| 9 | **`_check_invariance_1d` arm** | `:749-752` (rotational ⇒ True) | `:769-770` (⇒ reflection test) | **the semantically load-bearing one.** Fallthrough at `:772` is `False`. `Dnh` is the precedent: an `isinstance` line routing to `_is_reflection_invariant_1d` |
| 10 | **`_check_invariance_3d` arm** | `:855-859` | `:861-864` | `_orbit_closure(nodes, weights, _reflections(tag.axis), atol) is not None` |
| 11 | **`candidate_groups` enumeration** | `:1345-1346` | `:1345-1347` | the parameter must be **finitely enumerable**. `Cn`/`Dnh` use divisors of `_distinct_azimuths`; a mirror's parameter set is just `{x, y, z}` — **finite by construction, no bound needed**. This is strictly easier than `Cn`/`Dnh` |

**NOT needed** (contrast with the named entries):
- no `ClassVar` declaration (`:300-307`) — those are for singletons
- no `_NAMED_LATTICE` edges (`:203-240`) — the set is typed enum-to-enum
- no `SubgroupOfO3.X = _from_named(...)` install (`:1277-1284`)
- no `_named_contains` change (`:243-253`)

**Equality/hash come for free**: `__eq__`/`__hash__` (`:330-336`) delegate to
`self._tag == other._tag`, and a frozen dataclass gives structural equality.
Measured: `Mirror('x') == Mirror('x')` True, `== Mirror('z')` False. **But
`Mirror('z') == Z2` is `False`** even though they are the same matrix set — if
`Z2` survives alongside a mirror family, the tree carries **two unequal
spellings of one group**, and `contains` will say they contain each other while
`==` says they differ.

---

## 7. Things that surprised me

1. **Three of the five `Z2` lattice edges are dead code.** `_contains` computes
   finite × finite from matrices (`:621`); the static table is only consulted
   when a side is continuous. The table still *reads* as five load-bearing
   facts. A parameterized family needs **two** relations, not five.

2. **The tree is paying for the σ_x/σ_z distinction RIGHT NOW, one layer up,
   without naming it — in an UNCOMMITTED edit that landed under me mid-audit.**
   `_compute_sphere_reflection_partners` (`directional.py:157-215`) is being
   rewritten to build a **per-axis** certificate precisely because odd-`n_phi`
   products are closed under σ_z but not σ_x. That function IS
   `Mirror(axis).is_invariant` + `OrbitCertificate`, hand-spelled — and its new
   docstring cites the campaign's L-013 lesson *"a predicate that internally
   builds the permutation is the missing primitive"*. **The mirror-plane type is
   that same lesson one level down, and the two changes are in flight at the same
   time.** Worth reconciling before either lands: `directional.py` now imports
   `_orbit_closure` (a private) from `symmetry.py` to spell per-axis mirrors,
   which is exactly what a public `Mirror(axis)` would supply.

3. **The gap is observable on a SHIPPED rule with no synthetic fixture.**
   `product(4, 3)`: `σ_x = False`, `σ_z = True`, `Z2.is_invariant = True`.

4. **A symmetric μ-set structurally cannot see the bug.** Every existing gate
   uses Gauss-Legendre or Lebedev — sets that satisfy *all three* mirrors — so
   the axis question is vv-Mode-7 fixture-blind. The asymmetric control (which
   already exists in the tree, `test_symmetry.py:1193-1197`) is the discriminator,
   but it is only ever run through the **1-D** path where the plane does not exist.

5. **`_check_invariance_3d` does not call `_realized_ops`.** It re-spells
   `_reflections("z")` inline at `:819` and again at `:836`. Two independent
   spellings of one realization — exactly the drift `_realized_ops`'s own
   docstring (`:488-505`) says it exists to prevent ("*single-sources the operator
   sets from the very functions `_check_invariance_3d` applies*"). It
   single-sources `contains`, not the 3-D invariance arm.

6. **`test_dnh_reflection_in_dnh` (`test_symmetry.py:205-209`) asserts a
   generalization its own prose makes false.** *"`Z_2` (a single reflection) sits
   inside every `D_nh`"* — measured False for σ_x at n = 1, 3. It passes only
   because `Z2` is σ_z, the horizontal mirror `D_nh` carries by definition.

7. **A new tag type fails SILENTLY, not loudly.** `_contains` ends in
   `return False` (`:686`) and `_check_invariance_1d` in `return False` (`:772`).
   Measured: `O3.contains(Mirror('x')) = False` — flatly wrong, no exception. The
   two rotation-group answers are right *by luck*. Given `_realized_ops` returns
   `None` for the continuous groups and the module already shipped two false
   containment claims (`Z2 ⊆ SO3`, `Dnh ⊆ O2`), a `raise` on an unhandled tag
   pair would be the ratchet the existing law-tests
   (`test_computed_containment_obeys_the_order_relation_laws`, `:728`) already
   half-provide.

8. **`docs/theory/foundations/discrete_measures.rst:461` is already stale** — it
   still lists `O(2)` among the eight named entries, renamed `Dinfh` on
   2026-08-02. Pre-existing, unrelated to this carve, and **not caught by
   `-W`** (inline literal, not a Python-domain xref).

9. **`repr` is load-bearing infrastructure here, not cosmetics.** `_GROUP_CACHE`
   keys on `repr(tag)` (`:586`) and the walk's `visited` keys on `repr(group)`
   (`:1430`). The generic fallback happens to carry the axis; a "tidier" custom
   repr that drops it would silently collapse the three mirrors into one cache
   entry. And `tests/numerics/test_measure.py:599` asserts the repr **string**
   literally.
