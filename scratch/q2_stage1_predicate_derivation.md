# Q2 — Deriving the Stage-1 quadrature-selection predicate

Frame-attack deliverable (`cross-domain-attacker`, 2026-08-02). Brief: *derive* the
mathematically honest compatibility predicate for
`orpheus/numerics/quadrature/registry.py` Stage 1, given that the current one is
unsatisfiable once `SubgroupOfO3.is_invariant` is fixed.

## 0. Provenance — what is measured, what is derived

**No shell was provisioned in this dispatch** (no Bash tool), so nothing here was
run. Every claim is one of:

- `[M-user]` — measured by the requester and stated in the brief as verified.
- `[M-plan]` — a measured block already in `.claude/plans/quadrature_machinery_campaign.md`.
- `[D]` — derived here, with the derivation shown inline so it is checkable by
  reading.

Code facts are read directly from the worktree at
`orpheus/numerics/quadrature/registry.py`, `orpheus/numerics/symmetry.py`,
`orpheus/numerics/quadrature/rules_product.py`, `rules_1d.py` (branch
`main` working tree, 2026-08-02). `symmetry.py` was **read only**, never edited.

---

## 1. The answer, up front

**The predicate's SHAPE is right and both of its ARGUMENTS are wrong.**

`G_required ⊆ G_rule` is correct — it is an exact statement about the existence of
an ordinate **permutation**, and there is no approximate version of it. What is
wrong is:

1. **The left argument.** `GEOMETRY_GROUPS` records the *continuous* symmetry group
   of the **continuous** problem. The predicate needs the **discrete residual acting
   on the angular fiber** — a group that is finite **by theorem**, not by
   approximation. Framing **(D) is correct**, and it is forced, not chosen.
2. **The right argument.** `invariance_group` is a declared **constant** where the
   object is a **function of the rule's parameters**, and — one level deeper — a
   **computable property of the constructed measure** that nothing currently checks.

The derived predicate has **three** conjuncts where the code has one, and the two
new ones are not decoration: one of them (Stage 0) is currently enforced only by
accident, and the other (Sobolev) is the campaign's own gap 3 arriving from the
symmetry side.

**Verdict on `product` for `cylinder`: ADMITTED for even `n_φ`, REJECTED for odd
`n_φ`.** The odd-`n_φ` rejection **is ERR-042**, which the derived predicate
reproduces as a consequence rather than as a separate hand-written guard. That
reproduction is the strongest single piece of evidence that the derivation is the
right one.

---

## 2. Setup and notation

- `Q` — a quadrature: nodes `{x_i}`, weights `{w_i}`, on a domain `D` (either `S²`
  or the polar-cosine interval `[−1,1]`).
- `I[f] = ∫_D f dσ` — the exact integral, with `σ` the `O(3)`-invariant surface
  measure on `S²` (or Lebesgue on `[−1,1]`).
- `E[f] = Q[f] − I[f]` — the error functional.
- `g_*Q` — the pushforward of `Q` under `g ∈ O(3)`: nodes `{g x_i}`, weights `{w_i}`.
- `Sym(Q) := { g ∈ O(3) : g_*Q = Q }` — the **maximal invariance group** of the
  rule. A property of the node/weight set; not a tag.
- `V = P_d(S²)` — the rule's declared exactness space (spherical polys of degree
  ≤ d = `span{Y_ℓ^m : ℓ ≤ d}`), `Π_V` the `L²` orthogonal projection onto it.

Ordinate convention as in the tree: `Ω = (µ_x, µ_y, µ_z)`; for the cylinder
`η = µ_x` (radial), `ξ = µ_y` (tangential), `µ = µ_z` (axial).
`σ_x, σ_y, σ_z` denote reflections in the planes `x=0, y=0, z=0` respectively
(so `σ_y : ξ → −ξ`).

---

## 3. Framing (A), the error functional — **CONFIRMED**, and the reason matters

Claim under test: the requirement is `E[f∘g] = E[f]` for all `g ∈ G_geom`.

`[D]` `I` is `O(3)`-invariant, so `I[f∘g] = I[f]` and

```
E[f∘g] − E[f]  =  Q[f∘g] − Q[f]  =  ⟨ g_*Q − Q , f ⟩ .
```

Quantified over **all** `f` in a dense subspace of `C(S²)`, this vanishes iff
`g_*Q = Q` — i.e. **exact node-set closure with matched weights**. So (A) forces
closure. **The requester's reading is CONFIRMED, not refuted.**

And then, `[D]`: for a *finite* node set containing at least one node off the
rotation axis, no continuous group can be in `Sym(Q)` — the `SO(2)` orbit of an
off-axis node is an infinite circle, which a finite set cannot contain. Hence
`Sym(Q)` is a **closed subgroup of `O(3)` with no continuous part**, i.e. finite
(it embeds faithfully into `S_N` on the nodes, because an orthogonal map fixing a
spanning set of ℝ³ is the identity, and any real angular rule spans). Therefore:

> **THEOREM A.** For any quadrature with finitely many nodes not all on one axis,
> `Sym(Q)` is a **finite** subgroup of `O(3)`. Consequently a predicate
> `G ⊆ Sym(Q)` is satisfiable **only** for finite `G`.

`GEOMETRY_GROUPS` currently supplies `SO(2)` — continuous — for three of its four
entries. **The predicate is not merely failing; it is being fed an argument of the
wrong kind.** That is the finding, and it is what points at (D).

### 3.1 Why restricting the quantifier does **not** rescue (A)

The natural repair is to quantify over a smaller test space `W` instead of all `f`.
It does not help, and the reason is instructive.

`[D]` Since `1 ∈ V` (any `d ≥ 0`), `I[(Id−Π_V)f] = 0`, and since `Q` is exact on `V`,

```
E[f] = Q[f] − I[f] = Q[(Id − Π_V) f] .            (★)
```

**The error functional IS the quadrature applied to the aliased-out part.** That is
the aliasing structure the brief asks about, in one line.

Now `Π_V` commutes with the `O(3)` action (`V` is `O(3)`-invariant and the `L²`
inner product is `O(3)`-invariant), so `E[f∘g] = Q[((Id−Π_V)f)∘g]`. Requiring
`E[f∘g] = E[f]` therefore bites **entirely on `V^⊥`**, where there is no exactness
protection at all — and `(g_*Q − Q)` already annihilates `V`, so annihilating `V^⊥`
too means annihilating everything. Closure again.

> **The exactness of the rule is exactly what makes (A) collapse to closure:** the
> part of `L²` where the rule is accurate is invisible to the symmetry test, so the
> test is decided entirely by the part where it is not.

The only way out is a `W` that is *neither* `V` (vacuous) *nor* dense (closure) —
a genuinely band-limited `W ⊂ V^⊥ ⊕ V` with `dim W < ∞`. That gives a finite list
of linear conditions on the aliasing coefficients: an *approximate equivariance*.
**It is the wrong object for this solver**, for a reason that has nothing to do
with integration accuracy — see §6.1.

---

## 4. Framing (B), the exactness space — vacuous as stated; **non-vacuous as Sobolev**

`[D]` As stated it is vacuous, exactly as the brief suspects: `P_d(S²)` is
`O(3)`-invariant for **every** `d` (rotations preserve polynomial degree), so
`Q[f∘g] = I[f∘g] = I[f] = Q[f]` holds for every rule of degree ≥ d and every
`g ∈ O(3)`. A test every candidate passes carries zero information.

**But the aliasing reading the brief proposes is right, and it has a name.**

> **THEOREM B (isotypic annihilation).** Let `G` be finite with `G ⊆ Sym(Q)` and
> let `I` be `G`-invariant. Then `E` annihilates **every non-trivial isotypic
> component** of `L²(S²)` — at every polynomial degree, arbitrarily far above `d`.
>
> *Proof.* `G ⊆ Sym(Q)` gives `Q[f∘g] = Q[f]`, and `I[f∘g] = I[f]`, so `E[f∘g] = E[f]`.
> Let `f` lie in the `χ`-isotypic component with `χ ≠ 1`; then the trivial-isotype
> projector `P₁f = (1/|G|) Σ_g f∘g = 0`. Hence
> `E[f] = (1/|G|) Σ_g E[f∘g] = E[P₁ f] = E[0] = 0`. ∎

**Corollary (Sobolev 1962).** A `G`-invariant cubature formula has degree `t` **iff**
it is exact on the `G`-**invariant** polynomials of degree ≤ `t`. (S. L. Sobolev,
*Cubature formulas on the sphere which are invariant under transformations of finite
rotation groups*, Dokl. Akad. Nauk SSSR **146**(2), 1962, 310–313; Soviet Math.
Dokl. **3**, 1307–1310.) This is the founding theorem of invariant-cubature theory
and the reason Lebedev rules are constructed by solving only the `O_h`-invariant
moment conditions.

**This is the honest non-vacuous (B), and it is a different STAGE, not a rival
predicate:**

- Stage 1 (the G stage) is about the **group action**: an exact permutation
  condition. It is not an accuracy statement.
- Stage 2 (the V stage) is then, *by Sobolev*, only a claim about
  `L²(D)^{G}` — the invariant subspace. **That is precisely the campaign's gap 3
  ("`degree_of_exactness` carries no subspace") arriving from the symmetry side,
  derived rather than bolted on.**

Two consequences that land directly on open campaign items:

- **T6 / the `+2.94`.** The folded (quotient) measure's exactness claim is `d`
  **on `L²(S²)^G` only**. `[M-plan]` `full sphere ∫ξ = +4.4e-16` vs
  `quotient ∫ξ = +5.99`. Theorem B says the full-sphere rule integrates `ξ`
  exactly **because `ξ` is `σ_y`-odd, not because `ξ` is low degree** — the fold
  removes that protection because `ξ ∉ L²(S²)^{⟨σ_y⟩}`. The `+2.94` is the domain
  check firing, confirmed from theory, not just from measurement.
- **Nothing is lost by folding.** Sobolev says the full rule and the folded rule
  carry the *same* claim on the invariants. The fold is lossless on `L²(S²)^G`
  and meaningless off it. One theorem states both halves.

---

## 5. Framing (C), the maximal representable finite subgroup — **REFUTED, twice**

**(C1) Not well-defined.** `SO(2)` has **no** largest finite subgroup: its finite
subgroups are `{C_n}_{n≥1}`, a directed family with no maximum (`C_n ⊂ C_{2n}`
forever). "The largest finite subgroup of `G_geom` the discretization can
represent" is undefined without naming the discretization.

**(C2) Circular once you name it, and vacuous as a result.** Naming it makes the
requirement `C_{n_φ(Q)} ⊆ Sym(Q)` — **true for every product rule by
construction**, for the same reason (B)-naive is vacuous. And `n_φ` is an *output*
of selection (it comes from `degree_of_exactness_for`), so defining the admission
criterion in terms of it defines the filter in terms of the thing being filtered.

**(C3) Pointed at the wrong half of the group even after de-circularizing.** The
1-D cylinder does **not** require the discrete solution to be invariant under
rotation by `2π/n_φ` about `ẑ`: that rotation **moves the spatial point**
`(r,0,0)` to a different point, and the reduced problem has exactly one spatial
point per `r`. The `SO(2)_ẑ` symmetry has already been **consumed** by the
dimensional reduction — its infinitesimal generator **is** the α-redistribution
term (campaign T9). Requiring the ordinate set to be `C_{n_φ}`-closed is requiring
a symmetry that is not a symmetry of the reduced problem. That `product` happens to
satisfy it is a coincidence of the construction, not a physics requirement.

**Salvageable kernel (keep this half):** (C)'s instinct that `G_rule` should be
*computed from the node set* rather than declared is right, and it survives into
§9. Discard "largest representable subgroup of `G_geom`"; keep "maximal invariance
group of `Q`".

---

## 6. Framing (D), the geometry table — **DERIVED**, and it is the answer

The brief's hypothesis is that the table should record *"the symmetry group the
solution is required to exhibit on the ANGULAR variable"*. That is correct. Here is
the derivation that produces it, and it is a standard construction: the **slice /
isotropy (little-group) decomposition** of an equivariant PDE under a symmetry
reduction.

### 6.1 Why the requirement is a permutation, not an integral identity

The SN discrete unknowns are `(ψ_i)_{i=1..N}`, one per ordinate. A symmetry `g` of
the continuous problem acts on phase space as `(x,Ω) ↦ (gx, R_gΩ)`. For the
**discrete** system to inherit it, the map must be expressible on the unknowns:

```
ψ_i(x)  ↦  ψ_{π_g(i)}(g⁻¹x)          for a permutation π_g of the ordinate index
```

and the scattering sum `Σ_j w_j Σ_s(Ω_i·Ω_j) ψ_j` is `g`-equivariant iff
`w_{π_g(i)} = w_i`. Those two conditions are *exactly* `g ∈ Sym(Q)`.

> **There is no approximate version.** A permutation of a finite index set either
> exists or it does not. This is why §3.1's "approximate equivariance on a
> band-limited `W`" is the wrong object: it would buy small *integration* error on a
> few modes while the discrete system still has **no `g`-action at all**, and the
> discrete solution would still be generically asymmetric. `[M-plan]` T4 measures
> that asymmetry at `1.2e-1` when the action is broken elsewhere.

### 6.2 The decomposition — which part of `G_geom` goes where

Let `G` be the isometry group of the problem **after unfolding every specular
boundary by the method of images** (so a reflective face contributes its deck
transformation to `G`; this is the boundary campaign's ruling — *G is the deck
transformation, the mirror is an ambient isometry* — reused verbatim). Let `G⁰` be
its identity component (translations, `SO(2)` rotations) and `Γ = G/G⁰` the
discrete residual. Fix a representative point `x₀` of the reduced spatial domain.
Three parts, three fates, **none discarded**:

| part of `G` | what it does | price / requirement |
|---|---|---|
| `G⁰` acting **freely** at `x₀` | reduces the SPATIAL dimension | if its fiber action is non-trivial: the **redistribution term** `α` (a connection). Trivial fiber action (translations) ⟹ no α. **This is campaign T9.** |
| `G⁰_{x₀}` (continuous **isotropy** at `x₀`) | acts on the ANGULAR FIBER only, continuously | the angular domain is the **quotient** `S²/G⁰_{x₀}`. **Stage 0.** |
| `Γ` (the DISCRETE residual — mirrors, the `r=0` fold, specular deck maps) | acts on base and/or fiber discretely | its image in `Aut(S²/G⁰_{x₀})` must be realized as an **ordinate permutation**. **Stage 1.** |

Define

> **`G_ang(geometry) := image of Γ in Aut(S²/G⁰_{x₀})`** — the **required angular
> symmetry group**.

> **THEOREM D.** `G_ang` is **always finite**. *Proof.* `Γ = G/G⁰` is discrete by
> construction; its image in `O(3)` (compact) is a discrete subgroup of a compact
> group, hence finite. ∎

Theorem A said the predicate is satisfiable only for finite `G`; Theorem D says the
correctly-derived requirement is always finite. **Together they say the predicate
is always satisfiable in principle — the unsatisfiability was entirely an artifact
of feeding it the wrong group.**

### 6.3 Why the current table records the wrong half

`GEOMETRY_GROUPS` records (roughly) `G⁰` — the part that is *consumed* by the
reduction. That part is precisely the part that **cannot** constrain a finite node
set (Theorem A), and its constraint has already been paid elsewhere in the code, as
the α term. The table records the bill that was already settled and omits the one
still outstanding.

---

## 7. The derived predicate

For geometry `g` with angular domain `D_g = S²/G⁰_{x₀}` and required group
`G_ang(g)`, and a rule `Q` with parameters `p`:

```
Stage 0  (DOMAIN)      domain(Q)  ==  D_g
Stage 1  (SYMMETRY)    G_ang(g)  ⊆  Sym(Q(p))              [ ⊆ in O(3); exact ]
Stage 2  (EXACTNESS)   deg(Q)  ≥  d   as a claim on  L²(D_g)^{G_ang(g)}   [ Sobolev ]
Stage 3  (STRUCTURAL)  F_req ⊆ F_Q                          [ unchanged ]
Stage 4  (COST)        argmin n(Q)                          [ unchanged ]
```

Three notes on shape:

- **Stage 0 is currently enforced only by accident.** `gauss_legendre_on_mu`
  returns a measure on `[−1,1]`; `product`/`lebedev`/`level_symmetric` return
  measures on `S²`. Nothing in `select_quadrature` compares `measure.support` to
  the geometry's angular domain. `[D]` Today `select_quadrature("slab", 3)` has
  `ProductQuadrature` — an `S²` rule — as a **live candidate** (it passes Stage 1
  because `SO2 ⊆ SO2`), and loses only on node count (`2` vs `8`). The domain
  conjunct is real and missing.
- **Stage 0 and Stage 1 are the two halves of ONE criterion** — the continuous and
  discrete parts of the same isotropy. This is the same one-criterion structure the
  campaign already has in T9 (`n_levels > 1 ⟺ the point isotropy group is finite`).
- **`Sym(Q)` must be the MAXIMAL group**, not a nominal tag (§9). Under-claiming
  `Sym(Q)` makes the predicate *stricter* and rejects valid rules — which is exactly
  what happens to `GaussLegendre1D` under the corrected table (§11).

---

## 8. The geometry table under the derived predicate

Derivation per entry. `x₀` = a generic point of the **reduced** spatial domain.

### slab (normal `ẑ`)

- `G⁰` = translations in `x,y` (free, trivial fiber action ⟹ **no α** ✓) **and**
  `SO(2)_ẑ` rotations, which **fix** `x₀` ⟹ continuous isotropy.
- Stage 0: `D = S²/SO(2)_ẑ = [−1,1]` in `µ_z`. ✓ This is *why* the slab rule is
  `gauss_legendre_on_mu` and not an `S²` rule — derived, not assumed.
- `Γ`: the vertical mirrors are in `G⁰`'s closure and act **trivially** on `µ_z`.
  The non-trivial element is `σ_z : µ → −µ`, contributed by the slab-face specular
  deck map (and by the midplane symmetry when present).
- **`G_ang(slab) = Z₂` (`µ → −µ`).**

### sphere (1-D radial)

- `x₀ = (r,0,0)`, `r̂ = x̂`. `G = O(3)`; `G⁰_{x₀} = SO(2)_{r̂}` (continuous).
  The freely-acting rotations reduce 3-D → 1-D with **non-trivial** fiber action
  ⟹ the sphere's `α` term ✓.
- Stage 0: `D = S²/SO(2)_{r̂} = [−1,1]` in `µ_r`. ✓
- `Γ`: the `r = 0` centre fold. The reduced domain `r ∈ [0,R]` is the quotient of a
  diameter by `x ↦ −x`, whose fiber action is `µ_r → −µ_r`. Present in **every**
  solid sphere. (Mirrors containing `r̂` act trivially on `µ_r`.)
- **`G_ang(sphere) = Z₂` (`µ_r → −µ_r`).**

### cylinder (1-D radial, infinite `z`)

- `x₀ = (r,0,0)`, axis `ẑ`. `G⁰` = `ℝ_z` translations (free, trivial fiber ⟹ no
  contribution) and `SO(2)_ẑ` (free at `x₀` for `r>0`, **non-trivial** fiber action
  ⟹ the cylinder's `α` ✓). **`G⁰_{x₀}` is trivial.**
- Stage 0: `D = S²` — the **full** sphere. ✓ This is why the cylinder needs a 2-D
  angular rule, derived from the same criterion that gave the slab a 1-D one.
- `Γ`, three generators:
  - `σ_y` (plane `y=0`) — **fixes** `x₀` ⟹ isotropy. This is campaign **T4**
    (`ψ` even in `ξ`), here derived as an isotropy element.
  - `σ_z` (plane `z=0`) — **fixes** `x₀` (its `z` is 0) ⟹ isotropy, from
    `z`-uniformity.
  - `σ_x` (plane `x=0`) — from the `r=0` centreline fold and/or an outer specular
    face. Note `C₂(ẑ) = σ_xσ_y`, so recording `σ_x` and `σ_y` yields the
    centreline map **for free** — which is campaign T8's `σ_xσ_y = C₂(ẑ)`
    identity, appearing here as group closure rather than as a coincidence.
- `⟨σ_x, σ_y, σ_z⟩ = (ℤ₂)³`, order 8 = **`D_{2h}`**.
- **`G_ang(cylinder) = D_{2h}`** (with the `r = 0` fold or a specular radial face —
  i.e. every 1-D cylinder ORPHEUS runs).
  **`= C_{2v} = ⟨σ_y, σ_z⟩`, order 4**, for a vacuum-vacuum annulus with no `r=0`.
  The distinction is the decisive discriminator — see §10.

### cartesian2d

- The problem is `z`-uniform (otherwise it is 3-D), so `σ_z` **fixes every point**
  of the 2-D plane ⟹ isotropy, but **discrete**, so `G⁰_{x₀}` is trivial.
- Stage 0: `D = S²`. ✓ (Halving to `µ_z > 0` is then the Stage-1 group being *used*,
  not a separate convention.)
- `Γ`: `σ_z` (isotropy, from `z`-uniformity) plus `σ_x, σ_y` (coordinate-aligned
  specular faces).
- **`G_ang(cartesian2d) = D_{2h}`.**
  It is **NOT `O_h`**: `O_h` additionally demands the coordinate **permutations**
  `x↔y↔z`. `x↔y` is a symmetry only of a *square* domain with identical BCs on all
  four faces; `x↔z` is never a symmetry — `z` is the uniform direction. Even
  `D_{4h}` over-claims (square only).

### The table

| geometry | Stage-0 domain | `G_ang` | in-tree spelling |
|---|---|---|---|
| `slab` | `[−1,1]` (`= S²/SO(2)_ẑ`) | `Z₂` | `SubgroupOfO3.Z2` |
| `sphere` | `[−1,1]` (`= S²/SO(2)_{r̂}`) | `Z₂` | `SubgroupOfO3.Z2` |
| `cylinder` | `S²` | `D_{2h}` | `SubgroupOfO3.Dnh(2)` |
| `cartesian2d` | `S²` | `D_{2h}` | `SubgroupOfO3.Dnh(2)` |

Two observations worth keeping:

- **`cylinder` and `cartesian2d` land on the same group for structurally different
  reasons** — for the cylinder two of the three mirrors are spatial isotropy and one
  is a boundary deck map; for cartesian2d it is one and two. `D_{2h}` = "the
  coordinate sign-flip group" is the generic requirement of any coordinate-aligned
  discretization, and it is axis-labelling-agnostic (it contains all three mirrors),
  which sidesteps `SubgroupOfO3`'s missing principal-axis parameter entirely. That
  is a second, independent reason to prefer it over `C_{2v}` (which has principal
  axis `x̂` and is therefore *not* expressible in the current axis-free lattice).
- **The docstring's stated rationale for `O_h` is structurally inverted.** It reads
  *"Conservatively tagging `O_h` (rather than just `D_{4h}`) means the selector
  accepts any `O_h`-invariant rule."* In `G_geom ⊆ G_rule`, enlarging `G_geom`
  makes the predicate **stricter**. Tagging `O_h` is the *most restrictive*
  available choice, not a conservative one — it accepts **only** `O_h`-invariant
  rules. `[D]` It is why `product` is rejected for `cartesian2d` today
  (`_contains(Dnh, Oh)` returns `False` at `symmetry.py:488-495`).

### Where the BC-dependence honestly lives

`σ_x` for the cylinder and `σ_x, σ_y` for cartesian2d come from **boundary
conditions**, which a static geometry table cannot see. Two honest options:

- **Now:** record `D_{2h}` for both — the sound over-approximation. It rejects odd
  `n_φ`, which no production configuration uses and which ERR-042 already forbids.
- **Right shape:** `required_angular_group(geometry, *, specular_normals, includes_axis)
  -> SubgroupOfO3`, with the table holding only the geometry-intrinsic isotropy
  (`slab: Trivial`, `sphere: Trivial`, `cylinder: C_{2v}`, `cartesian2d: Z₂(σ_z)`)
  and the BC generators joined in. This is strictly better and it is what the
  boundary campaign's `G`/`R` split already implies. It is not required to fix
  Stage 1 today.

---

## 9. `invariance_group`: constant field or function of parameters?

**A function — and, one level down, a computed property with a verification gate.**

`[D]` The maximal group of each registered rule:

| rule | `Sym(Q)` | depends on params? | current tag | verdict |
|---|---|---|---|---|
| `gauss_legendre_on_mu(n)` | `D_{∞h}` about `ẑ` (≈ the tree's `O2`, which includes `σ_z` at `symmetry.py:673-675`) | no | `SO2` | **under-claim** — misses `σ_z`, the only element that matters |
| `lebedev_sphere(order)` | `O_h` | no | `O_h` | correct |
| `level_symmetric_sn(N)` | `O_h` | no | `O_h` | correct (equal weights within `O_h` orbits still give `O_h`; #327 falsifies the *degree*, not the group) |
| `product_mu_phi(n_µ, n_φ)` | `D_{n_φ h}` | **YES** | `SO2` | **false**, `[M-user]` |

Three constants and one function — and **the parameter-dependent one is exactly the
one that shipped a false tag.** That is not a coincidence: a `SubgroupOfO3` *field*
cannot express `D_{n_φ h}`, so the author wrote the nearest available untruth. This
is the **third** instance in this campaign of the same defect — *the type erases a
distinction the prose states*:

1. `degree_of_exactness: int` cannot say *which space* (T2 / T12b) → #327 shipped.
2. `invariance_group: SubgroupOfO3` cannot say *which parameters* → the `SO2` lie.
3. `parameters: dict[str, type]` cannot say *which admissible set* (T12d gap 2).

**Recommendation:** mirror the existing pattern exactly —

```python
invariance_group_for: Callable[[dict[str, Any]], SubgroupOfO3]
```

with the three constants written `lambda _p: SubgroupOfO3.OctahedralOh`. Uniform
callable, same shape as `degree_of_exactness_for`, no special case.

### 9.1 The deeper half: the declared group is a CLAIM, and nothing checks it

**`is_invariant` is never called by `select_quadrature`.** Stage 1 compares two
*tags* (`registry.py:677`). The declared `invariance_group` is therefore unverified
at every point of its life — which is exactly the #327 mechanism transplanted from
the V axis to the G axis. The campaign's own thesis fires again: *an unattached tag
cannot be checked against anything, so a false claim ships and a selector consumes
it.*

`[D]` `Sym(Q)` is **computable** and cheap, because the candidate lattice is finite
and small: take the maximal element of `{Trivial, Z₂, C_n, D_{nh}, O_h, I_h, O2, …}`
that `is_invariant` accepts, with `n` restricted to divisors of `N`. This is the
`_orbit_closure` certificate of campaign lesson **L-013** wearing its widened return
type — the permutation `π_M` it already computes and discards.

So the layering is:

- `invariance_group_for(params)` — the cheap pre-instantiation **claim** the
  selector reads (it must not instantiate every candidate).
- `measure.symmetry_group()` — the **truth**, computed from nodes + weights.
- a property test asserting they agree **exactly** (not `⊆`, not `⊇`) for every
  registered rule over a parameter sweep. **That test is what would have caught the
  `SO2` tag.**

### 9.2 A consequence for Q5's RANGE/SPACING/RULE/OFFSET factorization

`[M-plan]` §4 records that the azimuthal **offset** *"provably does not change
exactness."* True — and `[D]` **the offset does change `Sym(Q)`.** With
`φ_m = δ + 2πm/n`, the vertical mirror planes of `Sym(Q)` sit at angles
`δ + kπ/n`, so

```
D_{2h} ⊆ Sym(product with offset δ)   ⟺   n even  AND  δ ≡ 0  (mod π/n).
```

> **The offset is exactness-invisible and `Sym`-visible. Stage 1 is the ONLY
> conjunct in the whole machinery that can see it.** An offset chosen "for free"
> because it does not move the degree silently destroys the geometry's mirror.

This is the direct analogue of the Γ-centered-vs-shifted Monkhorst–Pack choice in
solid-state Brillouin-zone integration, where the shift is chosen precisely by which
point-group elements survive it. It is a concrete, load-bearing input to Q5 that the
exactness axis cannot supply.

---

## 10. The decisive discriminators

Cases where the surviving candidates give **different** answers, so the choice is
settled by mathematics.

### D-1 — `product(n_µ, n_φ)` at ODD `n_φ`, for `cylinder`. **The sharpest one.**

| predicate | answer at `n_φ` odd | answer at `n_φ` even |
|---|---|---|
| current (`SO2 ⊆ SO2`) | ADMIT | ADMIT |
| (C) `C_{n_φ} ⊆ Sym(Q)` | ADMIT (tautology) | ADMIT (tautology) |
| (D) with `G_ang = C_{2v}` | ADMIT | ADMIT |
| **(D) with `G_ang = D_{2h}`** | **REJECT** | **ADMIT** |

`[D]` The derivation: `D_{2h}` needs `σ_x : φ → π − φ`, i.e.
`π − 2πm/n = 2π(n/2 − m)/n`, which is a node **iff `n` is even**. Identically,
`C₂(ẑ) = ` rotation by `π = 2π(n/2)/n` needs `n` even. `C_{2v} = ⟨σ_y, σ_z⟩` needs
only `σ_y : φ → −φ`, `−2πm/n = 2π(n−m)/n` ✓ **for every `n`**, and `σ_z : µ_z → −µ_z`
✓ (GL nodes symmetric).

> **`G_ang(cylinder) = D_{2h}` reproduces ERR-042 — "even `n_φ`" — as a
> CONSEQUENCE of the selection predicate.** The current predicate cannot express it
> at all; ERR-042 exists as a separate hand-written guard elsewhere in the tree.
> `[M-plan]` T8 records that ORPHEUS "checks the wrong group element and agrees with
> the right one only because `σ_y`-closure is automatic" — the derived predicate
> checks the right one.

This also discriminates **within** (D): `C_{2v}` vs `D_{2h}` is decided by whether
`r = 0` (or a specular radial face) is in the domain. It is a real physical
question with a definite answer per configuration, not a taste call.

### D-2 — an ASYMMETRIC 1-D rule for `slab` / `sphere`. **Live for Q3.**

`[D]` A Gauss–**Radau** rule on `[−1,1]` (one endpoint constrained; campaign Q3 adds
exactly this node-constraint axis) has non-symmetric nodes ⟹ `Sym = Trivial` on the
`µ` axis.

| predicate | Gauss-Radau for `slab` |
|---|---|
| current (`GEOMETRY_GROUPS["slab"] = SO2`) | **ADMIT** — and `_check_invariance_1d` short-circuits `SO2` on a 1-D measure to `True` (`symmetry.py:592-599`), so even calling `is_invariant` would not catch it |
| **derived (`Z₂`)** | **REJECT** — `_is_reflection_invariant_1d` fails on the asymmetric node set |

The current predicate has **zero teeth** on the one requirement the slab actually
has. Gauss-Lobatto (symmetric, includes `±1`) correctly passes both — so the test
discriminates Radau from Lobatto, which is the real distinction in Q3's
node-constraint sum type.

### D-3 — `product` for `cartesian2d`.

| predicate | answer |
|---|---|
| current (`O_h ⊆ D_{n_φ h}`) | **REJECT** (`_contains(Dnh, Oh) → False`) |
| **derived (`D_{2h} ⊆ D_{n_φ h}`)** | **ADMIT** for even `n_φ` |

A live behavioural change: 2-D Cartesian gains product quadrature, which is correct —
`D_{2h}` is what a rectangular 2-D mesh with coordinate-aligned specular faces
actually requires, and `O_h` was demanding a square-with-`x↔z` symmetry that no
2-D problem has.

### D-4 — the isotypic-exactness property (Theorem B), directly.

`[D]` Take `f = µ_x³ µ_y` (degree 4), which is `σ_y`-**odd** and `C₂(ẑ)`-**even**,
with `I[f] = 0`. Compare `product(2, 4)` (declared degree `min(3,3) = 3`, so `f` is
**above** the exactness space) against a `φ`-offset variant `φ_m = δ + 2πm/4`,
`δ = 0.3`:

- un-offset: `σ_y ∈ Sym(Q)` ⟹ Theorem B forces `Q[f] = 0` to machine precision,
  **despite `deg f > d`**.
- offset: `σ_y ∉ Sym(Q)`; the `sin 4φ` harmonic in `f` aliases (`4 ≡ 0 mod 4`) ⟹
  `Q[f] = O(sin 4δ) ≠ 0`.

This discriminates framing (B)-naive (which predicts nothing above degree `d`) from
(B)-Sobolev (which predicts **exact zero at every degree** on non-trivial isotypes).
It is also the cleanest possible demonstration that the symmetry conjunct buys
something the degree conjunct cannot.

---

## 11. Prerequisites — the derived predicate does NOT drop in

`[D]` Three in-tree defects become **live** the moment `GEOMETRY_GROUPS` carries
`Z₂` and `Dnh(2)`. Each is currently invisible because no `Dnh` and no `Z2` is ever
passed as a geometry group.

### P-1 — `_contains` says `D_{2h} ⊄ O_h`. **False, and it would reject Lebedev.**

`symmetry.py:515-518`: `outer ∈ _NamedSubgroup`, `inner ∈ Dnh` ⟹
`return outer in (O2, O3)`. So `Dnh(2).is_subgroup_of(OctahedralOh)` is `False`.
Mathematically `D_{2h} = (ℤ₂)³ ⊂ O_h` (order 8 | 48) ✓. Under the derived table this
would **reject Lebedev and level-symmetric for both `cylinder` and `cartesian2d`** —
a severe regression.

Required fix (narrow and exact): `O_h ⊇ D_{nh}` iff `n ∈ {1, 2, 4}` (no 6-fold
axis; `D_{3h}` needs a `σ_h` perpendicular to a body diagonal, which `O_h` lacks),
and `O_h ⊇ C_n` iff `n ∈ {1, 2, 3, 4}` — the case `symmetry.py:498-513` explicitly
declines to encode. Both carry an implicit "principal axis is a coordinate axis"
convention that must be documented.

### P-2 — `GaussLegendre1D`'s tag must change to `O2` in the same commit.

`[D]` `Z2 ⊄ SO2` in the lattice (`symmetry.py:199-205`, comment: *"Z2 ⊂ everything
except SO2"*), and that is mathematically right — a continuous proper-rotation group
has no reflection. So with `GEOMETRY_GROUPS["slab"] = Z2` and
`GaussLegendre1D.invariance_group = SO2`, **the slab's own canonical rule is
rejected.** The tag must become `O2` (the tree's `O2` includes `σ_z`, see
`symmetry.py:673-675`), which is also its true maximal group.

> **Table and tags change together. Changing either alone breaks selection
> completely.** This is the blast-radius statement for the carve.

### P-3 — `_vertical_mirrors` convention is rotated by `π/(2n)` relative to `Sym(product)`.

`[D]` `symmetry.py:787-799` places mirror **normals** at `kπ/n` ⟹ mirror **planes**
at `kπ/n + π/2`. At `k=0` that is the plane `x=0` (`σ_x`, `φ → π−φ`). The product
rule always owns the plane `y=0` (`φ → −φ`) and owns `x=0` only for even `n`.
For **even** `n` the two plane-sets coincide (`π/2 = (n/2)(π/n)` shifts the set onto
itself); for **odd** `n` they differ, so `is_invariant(Dnh(n_φ))` on
`product(·, odd n_φ)` returns `False` even though `Sym(product) = D_{n_φ h}`.
**Latent** (production `n_φ ∈ {8,16,32}`), **live** the moment an odd `n_φ` is
tested or `Dnh` is used as a claimed rule group. Flagged, not fixed — `symmetry.py`
is being carved concurrently.

---

## 12. Does the derived predicate ADMIT or REJECT `product` for `cylinder`?

**ADMITS, for even `n_φ`. REJECTS for odd `n_φ`.** Derivation in §10 D-1.
Production runs `n_φ ∈ {8, 16, 32}` ⟹ **admitted**.

### And it does NOT catch #326 — which is correct, not a shortfall

`[M-plan]` T4 measures `product(4,8)` homogeneous at `1.19e-1` asymmetry in `ξ`,
and states the reason: *"The product rule with equal azimuthal weights preserves the
symmetry, so the semi-discrete problem inherits it exactly."* The quadrature **is**
`σ_y`-invariant. #326 is a defect in the **solver** (the ordering / α-recursion,
T3 / T7), not in the quadrature's symmetry.

> A Stage-1 predicate that "caught" #326 would be blaming the quadrature for a sweep
> bug, and the resulting fix would be wrong. The derived predicate admits `product`
> and stays silent on #326 — that is the right behaviour and it should be stated
> explicitly so nobody reads the admission as a miss.

What the derived predicate **does** contribute to #326: it makes
`G_ang(cylinder) = D_{2h} ⊆ Sym(Q)` a checked precondition, and `D_{2h}` is exactly
the group whose element `C₂(ẑ) = σ_xσ_y` is the `r=0` centreline map that T8 shows
the tree currently implements as `σ_x` — *"the tree is right by a coincidence whose
precondition it has not taken."* The predicate is the place that precondition gets
taken.

---

## 13. First tests — each one can RED

Five tests, each stating a property the derived frame predicts that a naive
implementation **violates**.

1. **`test_declared_invariance_group_is_exactly_the_maximal_one`** — for every
   registered spec and every parameter point in a sweep, assert
   `invariance_group_for(p) == measure.symmetry_group()` (**equality**, not
   containment). REDs today on `ProductQuadrature` at every `(n_µ, n_φ)` and on
   `GaussLegendre1D` at every `n`. *This is the test that would have caught the
   `SO2` tag; it is the G-axis twin of the test that catches #327 on the V axis.*
2. **`test_odd_azimuthal_count_is_rejected_for_cylinder`** — register
   `product(4, 7)`; assert `select_quadrature("cylinder", …)` rejects it with a
   Stage-1 reason. REDs today (accepted). *Discriminates `D_{2h}` from `C_{2v}` and
   from every current candidate; it is ERR-042 arriving through the predicate.*
3. **`test_asymmetric_1d_rule_is_rejected_for_slab`** — register a Gauss-Radau rule
   on `[−1,1]`; assert rejection for `slab` and `sphere`; assert Gauss-**Lobatto**
   (symmetric) is **accepted**. REDs today (Radau accepted, and would be *chosen*
   when cheaper). *Negative test with a matched positive control, so it cannot pass
   by rejecting everything.*
4. **`test_isotypic_components_are_integrated_exactly_above_the_declared_degree`** —
   `Q = product(2,4)` (declared `d = 3`); assert
   `|Q[µ_x³µ_y]| < 1e-15` (degree 4, `σ_y`-odd) and that the `δ = 0.3`
   `φ`-offset variant gives `|Q[µ_x³µ_y]| > 1e-3`. *Discriminates Sobolev from
   degree-only reasoning, and pins §9.2's offset finding.*
5. **`test_sphere_rule_is_rejected_for_a_cylinder_domain`** — assert
   `select_quadrature("slab", d)` never returns an `S²`-supported measure, by
   registering a cheap `S²` rule (`product(1,2)`, 2 nodes) that would otherwise win
   Stage 4 on cost. REDs today. *Puts teeth on Stage 0, which is currently
   enforced only by an accident of node counts.*

---

## 14. Summary of answers to the five questions

1. **The derived predicate** — §7. `G_ang(g) ⊆ Sym(Q(p))`, where `G_ang` is the
   image on the angular fiber of the problem's **discrete symmetry residual**
   `Γ = G/G⁰` (after unfolding specular boundaries), and `Sym(Q)` is the rule's
   **maximal** invariance group. Preceded by a domain conjunct
   `domain(Q) = S²/G⁰_{x₀}` and followed by a Sobolev-restricted exactness conjunct.
   Named results: **slice / isotropy (little-group) decomposition** of an
   equivariant reduction for the left argument; **Sobolev's invariant-cubature
   theorem (1962)** for the exactness conjunct; **Theorem B (isotypic annihilation)**
   for the aliasing structure.
2. **Decisive discriminator** — §10, four of them. Sharpest: odd `n_φ` for
   `cylinder` (D-1), where the derived predicate reproduces ERR-042 and every other
   candidate admits.
3. **Constant or function** — **function**, `invariance_group_for(params)`, mirroring
   `degree_of_exactness_for`; and one level deeper, a **computed** property of the
   measure with an equality gate against the declared claim (§9, §9.1). The
   parameter-dependence is carried by exactly one of four rules — the one that
   shipped the falsehood.
4. **`GEOMETRY_GROUPS`** — §8. `slab: Z₂`, `sphere: Z₂`, `cylinder: D_{2h}`,
   `cartesian2d: D_{2h}`; plus a Stage-0 domain column (`[−1,1]`, `[−1,1]`, `S²`,
   `S²`). `GaussLegendre1D`'s tag must move `SO2 → O2` in the same commit (P-2), and
   `_contains(O_h, D_{nh})` must be fixed first (P-1).
5. **`product` for `cylinder`** — **ADMITTED** for even `n_φ` (production), rejected
   for odd. It does **not** and **should not** catch #326 (§12).

---

## 15. Cross-method pollination

**Current method class:** SN angular discretisation / quadrature selection.

- **Solid-state Brillouin-zone integration (Monkhorst–Pack / Chadi–Cohen special
  points).** Trigger: a finite point set on a compact domain, chosen by a point
  group, folded to an irreducible wedge with orbit-size weights. Same construction
  as campaign T5's orbifold measure. It brings a ready-made result ORPHEUS needs:
  the **Γ-centered vs shifted mesh** decision is made *precisely* by which
  point-group elements the shift preserves — the exact content of §9.2, and the
  known MP pathology (a mesh of lower symmetry than the crystal breaks the symmetry
  of the computed density) is the failure Stage 1 prevents.
  *First test:* §13 item 4's offset arm.
- **Quantum-chemistry DFT grids (symmetry-adapted Lebedev).** Trigger: a molecular
  point group + an angular grid, with the standard practice being "choose a grid
  whose group contains the molecule's group, then evaluate only the
  symmetry-unique points." Identical predicate, identical fold, decades of
  production evidence. *First test:* assert the folded-rule node count equals
  Burnside's `(N+F)/2`, which `[M-plan]` T5 already measures (20/8/12 and 36/8/28).
- **MoC (in-tree, live second consumer).** `moc/geometry.py:222-229` puts the
  azimuthal grid on `[0,π)` — MoC **already took the Stage-0 quotient**
  `S¹/⟨−I⟩ = RP¹` and never named it. Trigger: the same three-stage predicate with
  `G⁰_{x₀}` trivial and `−I ∈ Γ`. A `G_ang` that serves only SN would be serving one
  caller. *First test:* assert MoC's azimuthal measure passes
  `G_ang(moc_2d) ⊆ Sym(Q)` with `G_ang ∋ −I`, and that its `_reflected_azi_index`
  agrees `array_equal` with the orbit certificate's permutation (0 ULP).

---

## 16. UNEXPLORED — frames checked, not triggered

- **Wiener–Hopf / Chandrasekhar H-functions** — wrong solver family; this is a
  discretisation-selection question, not a half-space transport kernel.
- **Homology / chain complex** — the word "boundary" appears (specular deck maps),
  but `∂∘∂ ≠ 0`: two reflections compose to a non-trivial rotation. No `d² = 0`, no
  homology payoff.
- **Category theory / functorial cubature** — the compositional win here
  (quotient ∘ restriction, group-indexed families) is already fully captured
  concretely by the isotypic decomposition and the orbit certificate. No abstract
  lever produces a test the concrete frame does not.
- **Tensor networks / MPO** — the product rule is a rank-1 tensor product
  (`GL(µ) ⊗ trapezoid(φ)`) with no bond-dimension knob. Degenerate, not a network.
- **Differential geometry / connection coefficients** — fires on the α
  redistribution term (already established, campaign T3/T9) but **not** on the
  selection predicate: this question has no curvature term to redistribute, only a
  group action to represent.
- **Spherical designs / `t`-designs** — genuinely adjacent (equal-weight cubature
  from group orbits) and would be the right frame if the campaign were *constructing*
  new rules. It does not fire on *selecting* among existing ones, because a design's
  defining property is an exactness claim, i.e. Stage 2, already covered by Sobolev.
  Trigger to revisit: if #327 is fixed by **rebuilding** `level_symmetric` rather
  than renaming it (campaign Q8), the design/orbit-construction literature is the
  correct source.
- **Koksma–Hlawka / QMC discrepancy** — no trigger: the requirement derived here is
  an exact permutation condition, not an integration-error bound. Revisit at #128
  (Hardy–Krause variation space), which is a Stage-2 question.
