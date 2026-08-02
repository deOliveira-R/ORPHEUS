# Issue #326 — frame attack on the proposed quadrature-quotient machinery

**Cross-domain-attacker, structural detection.** Design analysis only; no production
code written. Branch `refactor/operator-strategy-layers`.

## Verdict in one paragraph

The framing is **not wrong**. "Half-range is a quotient, not a restriction" is
correct and is the right root. But the proposal is **under-resolved in three places
and mis-aimed in one**, and the mis-aim is the load-bearing one:

1. **Claim 2 fuses two objects that live in different categories.** The
   *measure-side* quotient is **not** a `Frame` — a frame's analysis lands in a
   **coefficient** space; a quotient's output is another **measure**. The
   *basis-side* restriction **is** a genuine `Basis`/`GalerkinFrame`
   (symmetry-adapted spherical harmonics). Split them; do not fuse them.
2. **Claim 5 mis-names the gap.** A level is **not an orbit** — it is a **fiber of
   an invariant**, and `product` and `level_symmetric` fiber over **two different
   invariants** (`μ_z` signed vs `|μ_z|`) inside **one** type. The missing
   primitive is one layer below `LevelStructure`: **the group ACTION on nodes** (the
   permutation representation), which the tree already computes **three times**
   independently and throws away twice.
3. **Claim 1's precondition is one conjunct of three**, and the field it names is
   documented in-tree as an upper bound that is literally false at finite `n_φ`.
4. **Claims 3 and 4 survive intact** and can be upgraded from observation to
   theorem (Burnside + orbit-stabilizer + character count reproduce every measured
   number, including the trapezoid `[w,2w,…,2w,w]` and the `+2.94`).

The single highest-value thing missing from all five claims: **the orbifold
singular set** — the fixed-point locus `ξ = 0` — is currently re-derived by ad-hoc
`ε`-tests at **six** unconnected sites, and every reported pathology in #326 lives
on it.

## Evidence discipline

I have **no shell in this dispatch**, so I produced **no new measurements**.

| Tag | Meaning |
|---|---|
| `[C]` | Code fact, `file:line`, read directly from the working tree this session (L-005: worktree, not Nexus) |
| `[M†]` | Measured in a **prior** dispatch — cited to `scratch/issue326_halfrange_map.md` or `scratch/issue326_alpha_theory.md` by section. Not re-run. |
| `[T]` | Derived here; checkable by hand, no code needed |

Where a `[T]` reproduces an `[M†]` number exactly, it is flagged `[T ✓ M†]` — that
agreement is the strongest evidence in this document, because the derivation and
the measurement are structurally independent.

---

# STRUCTURAL FEATURES

**Mathematical objects**

1. A compact group `O(3)` with a named finite/continuous subgroup lattice, already
   modelled as a type (`SubgroupOfO3`, `orpheus/numerics/symmetry.py:100-458`) `[C]`.
2. A discrete measure on `S²` with `support`, `weights`, `invariance_group`,
   `degree_of_exactness` (`orpheus/numerics/measure.py`) `[C]`.
3. Two *different* reduced phase spaces glued into one type: `LevelStructure`
   = `(n_levels, level_indices, level_mu)` (`rules_sphere.py:111-136`) `[C]`.
4. A `Basis`/`FrameBase` hierarchy with analysis/reconstruction faces and a
   discipline-as-type ladder `FrameBase → PetrovGalerkinFrame → GalerkinFrame`
   (`orpheus/numerics/frame.py:113-380`) `[C]`.
5. A boundary algebra `γ₋ψ = R G γ₊ψ + q` in which **`G` is explicitly named the
   deck transformation of a quotient of the domain**
   (`orpheus/geometry/boundary/__init__.py:19-46`) `[C]`.

**Symmetries present**

6. 1-D infinite cylinder spatial group: `(ℝ_z ⋊ ⟨σ_z⟩) × O(2)_ẑ`. Translations act
   **trivially** on directions; rotations act **non-trivially** `[T]`.
7. Isotropy (stabilizer) subgroup at a radial point `(r,0,0)`:
   `H = ⟨σ_y⟩ × ⟨σ_z⟩ ≅ Z₂ × Z₂ = C_{2v}` with principal axis **x̂ = the radial
   direction** `[T]`. `σ_y` is the ξ-mirror; `σ_z` is the μ_z-mirror.
8. `σ_x σ_y = C₂(ẑ)` — an exact group relation, not an approximation `[T]`.
9. `product_mu_phi`'s measure is **`C_{n_φ}`-invariant, not `SO(2)`-invariant**; the
   module docstring says so verbatim and calls the `SO2` tag "a conservative upper
   bound" (`rules_product.py:24-35`, `:147`) `[C]`.

**Iterative / differential / integral structure**

10. Per-level first-order recurrence `α_{m+1} = α_m − w_m η_m` — a cumulative
    quadrature of `dξ/dω = η` along a 1-D chain (`reduced_operator.py:778-786`) `[C]`.
11. Two-sided Dirichlet closure `α_{1/2} = α_{M+1/2} = 0` (BMC Eq. 50) `[M† alpha §1c]`.
12. Weighted-diamond angular closure with `τ ∈ [½,1]` clamp
    (`pole_angular_closure.py:588-675`) `[C]`.
13. A `Y_ℓ^m` moment integral against the same measure (`directional.py:401-458`) `[C]`.

**Boundary handling**

14. Angular closure `α = 0` at both `ξ = 0` planes — a **zero-flux condition at the
    two ends of the reduced angular interval** `[T]`.
15. Spatial `r = 0` centreline handled as `reflection_index("x")`
    (`loss_representation/__init__.py:3294, :3470, :4189, :4593`) `[M† map §3a]`,
    while the physical continuation is `C₂(ẑ)`.
16. `ERR-042` rejects odd `n_φ` on a cylinder with
    `BoundaryGeometryMapNotMeasurePreservingError` from
    `orpheus/geometry/boundary/_specular.py:154` `[M† alpha §5]`.

**Scale separation** — none relevant. No thick/thin limit in play.

**Where the elegance detector fires** — see next section.

---

# ELEGANCE DETECTOR HITS

**Smell #4 (symmetry present in the problem, absent in the method) — PRIMARY.**
The physical isotropy group of the 1-D cylinder is `C_{2v}` (order 4). `product`
quotients by **nothing** (its level is the signed-`μ_z` fiber, containing both
`ξ` signs); `level_symmetric` groups by `|μ_z|` but its level still contains **both**
signs of `μ_z` **and** both signs of `ξ` (`rules_sphere.py:211-215`) `[C]`. Neither
rule realizes the quotient. This is the issue.

**Smell #16 Shape 1 (two code paths for one operator) — THREE independent
implementations of the SAME object: "the permutation of quadrature nodes induced by
a group element."**

| Site | Implementation | What it does with the permutation |
|---|---|---|
| `orpheus/numerics/symmetry.py:904-954` `_orbit_closure` | quantised-key dict + brute-force nearest-neighbour fallback | computes the matched index `j` per `(i, M)` and **throws it away** — returns `bool` |
| `orpheus/numerics/quadrature/directional.py:170-187` `_compute_sphere_reflection_partners` | `_find_reflections` on the three axis reflections | **keeps** it, as `Quadrature.reflection_partners` |
| `orpheus/moc/geometry.py:222-229` `_reflected_azi_index` | `argmin(abs(phi − (π−phi)))` + hand-rolled wrap | keeps it, per-azimuthal-index |

`[C]` all three. `_orbit_closure` is *literally the fold, minus the output*. This is
the missing primitive, and it is missing **because the return type was chosen as
`bool`**.

**Smell #16 Shape 2 (one quantity, two incompatible representations).**
`LevelStructure` is populated by two producers keyed on **different invariants**:
`rules_product.py:139` fibers over signed `mu_x`-sorted signed-`μ_z` levels;
`rules_sphere.py:213` fibers over `|μ_z|` (`np.abs(np.abs(mu_z_arr) − level_mu[p])`)
`[C]`. Every downstream consumer reads `level_indices` and assumes one meaning.
This is why `τ_raw` is `[0,1,0,1,…]` for one and `[1,½,½,0]` for the other
`[M† alpha F3]` — **two different residual group orders (2 vs 4) wearing one type.**

**Smell #12 ("we had to add a small ε") — six sites, one locus.**
The fixed-point set of `σ_y` (the `ξ = 0` nodes) is re-detected by hand at:
`augmented_mesh.py:815` (`eps < tau_level[0] < 1−eps`),
`pole_angular_closure.py:1012` (`abs(mu[cand] − mu[m0]) > 1e-14`),
`loss_representation/__init__.py:2902` (`|η| < 1e-15` degenerate positions),
`rules_sphere.py:212` (`tol = 1e-12` level membership),
`symmetry.py:952` (`atol` weight match),
`directional.py:477` (`|μ_axis| < 1e-15` octant-zero label) `[C]`.
Six epsilons for one geometric object: **the orbifold singular locus**.

**Smell #9 (boundary handling as a special case bolted onto interior logic).**
`Quadrature.n_levels / level_indices / level_mu` (`directional.py:488-523`) are three
`if self.level_structure is None: return <degenerate>` branches `[C]`. The slab and
sphere are not "no levels" — they are **one** level whose fiber the stabilizer
collapses to a point `[T]`. The null-branch is the degenerate case spelled as absence.

**Smell #3 ("magic correction term with hand-waved derivation").** Half-fires. The
α term's derivation is now clean `[M† alpha §2]`, but *why the cylinder has one and
the slab does not* is nowhere stated. Frame D below supplies the one-line criterion.

Six smells. The threshold is two.

---

# CLAIM-BY-CLAIM ADJUDICATION

## Claim 1 — "half-range is the quotient of a G-invariant measure; `is_invariant` is that precondition"

**Confirmed on the quotient half. Refuted on the precondition half — three ways.**

**Confirmed.** The discriminator between *restriction* and *quotient* is already
spelled out in the tree, in two docstrings, and it is **total mass**:

* `DiscreteMeasure.restrict` (`measure.py:624-668`) — "Drops atoms outside `E`;
  weights of kept atoms are **preserved**" ⟹ `Σw` **falls** `[C]`.
* `DiscreteMeasure.partition_by` (`measure.py:674-697`) — "**preserves total mass**"
  `[C]`.

The measured fold **preserves `Σw`** `[M† map §2d, "weight sum preserved to 3e-15"]`.
So it is categorically not a `restrict`. **Confirmed.**

**Refutation 1 — the word is already spent, on the other operation.**
`restrict`'s own docstring says: *"Used for **half-range SN sweeps**
(`E = {μ > 0}`)"* (`measure.py:632-633`) `[C]`. The slab's `μ>0` half-range **is** a
restriction (mass drops; no group identifies the two halves; the Marshak partial
current in `acceleration.rst` is that sense `[M† map §1b]`). So **"half-range" names
two structurally different operations in one tree.** Per L-012's delete-it-and-ask
check: delete the slab half-range ⟹ you lose a partial current (a different, still
well-posed functional). Delete the cylinder fold ⟹ ψ is no longer even in ξ ⟹
**ill-posed**. Different animals. The new object must **not** be called half-range.
Name it the **fold / quotient**; reserve half-range for the restriction.

**Refutation 2 — `is_invariant` is one conjunct of three, and the sufficient one is
missing.** The fold is valid iff:

| # | Condition | Where it lives today |
|---|---|---|
| (a) the **measure** is `G`-invariant | `SubgroupOfO3.is_invariant` `[C symmetry.py:396]` — exists |
| (b) the **operator** is `G`-equivariant (`A∘g = g∘A`) | **nowhere** — no `commutes_with` on any operator `[C]` |
| (c) the **data** (source, inflow trace, BC) is `G`-invariant | **nowhere**; named as an obligation only in prose `[M† map §3a]` |

(a) alone certifies that the *quadrature* folds. (b)+(c) certify that the *solve*
folds. Without (b), `conjugate`-style descent **silently returns the symmetrized
operator** `(1/|G|)Σ_g g A g⁻¹` instead of raising — a wrong number with a valid
type. That is the landmine, and Claim 1 as written does not name it.

**Refutation 3 — the field the claim leans on is documented in-tree as not
literally true.** `rules_product.py:24-35` `[C]`, verbatim: *"The maximal
continuous-rotation symmetry the discrete rule respects is `SO(2)` only when
`n_φ → ∞`; for finite `n_φ` the rule is `C_{n_φ}`-invariant strictly. For the
algebra-of-record, we tag `SO(2)` as a conservative upper bound."* And
`is_invariant`'s own docstring (`symmetry.py:410-419`) `[C]`: for continuous groups
the check *"uses a representative orbit … this is **necessary but not sufficient**."*
`_so2_representatives()` returns **four** rotations (`symmetry.py:769-779`) `[C]`.
So `is_invariant(SO2, product_measure) == True` certifies `C₄`-closure, not `SO(2)`.

For the `ξ`-mirror `Z₂` the check is exact (a single reflection, matched exactly), so
the fold is safe — but a design that **gates on
`measure.invariance_group.contains(G)`** would be gating on a declared upper bound.
Gate on the **action**, not the tag: the quotient exists iff the permutation exists.

**And `is_invariant` is not even the checking face of the quotient — it is a boolean
SHADOW of the primitive the quotient needs.** `is_invariant(G, μ)` answers *"is `μ`
`G`-invariant?"*. The quotient needs *"give me the orbit partition"*, which requires
the **action**. `_orbit_closure` (`symmetry.py:904-954`) computes the matched index
`j` for every `(i, M)` pair and **returns `bool`** `[C]`. So Claim 1's "the machinery
has the checking face and lacks the quotient" understates the situation in the
project's favour: the machinery **computes the quotient and discards it**, because
the return type was chosen as `bool`. That is a strictly better starting point than
the claim asserts, and it changes the size of the fix (M2).

### Which `n_φ` admits which mirror — the claim's existence condition is on the WRONG group element

`[T]` For the equispaced left-endpoint grid `φ_k = 2πk/n_φ`:

| element | action on `φ` | grid closed? | fixed points per level |
|---|---|---|---|
| `σ_y` (the **fold's** mirror, `ξ → −ξ`) | `φ → −φ` (`k → n_φ − k`) | **always**, any `n_φ` | `φ = 0` always; `φ = π` iff `n_φ` even |
| `σ_x` (the **centreline's** mirror as coded, `η → −η`) | `φ → π − φ` | **iff `n_φ` even** | `φ = π/2, 3π/2` iff `4 ∣ n_φ` |
| `C₂(ẑ)` (the **physical** centreline map) | `φ → φ + π` | **iff `n_φ` even** | none |

**So the `ξ`-fold's existence is automatic and unconditional; it is the CENTRELINE
map whose existence needs even `n_φ`.** ERR-042 — *"`reflection_index('x')` does not
preserve the direction-cosine measure `w·|μ_x|`"* (`_specular.py:154` via
`reflective.py:171`) `[M† alpha §5]` — is therefore **not** the fold's existence
condition. It is the condition that the cover admits the **centreline deck
transformation**, and the group relation `C₂(ẑ) = σ_x σ_y` `[T]` explains why
checking `σ_x` gets the right answer: `σ_y`-closure is automatic, so `σ_x`-closure
and `C₂(ẑ)`-closure are equivalent. The two rows agree in the table for exactly that
reason.

Consequences, both design-relevant:

1. **Gating the fold on any invariance predicate is doubly wrong** — the fold's
   condition is automatic, and the condition that is *not* automatic belongs to a
   different group element, in a different subsystem (M7).
2. **The check ORPHEUS should be making is on `C₂(ẑ)`, not `σ_x`.** They coincide
   only because `σ_y`-closure holds; the moment a non-product azimuthal rule (or a
   folded measure) enters, `σ_x` stops being the physical map and the equivalence
   silently lapses. Same disease as Q2(iii), detected at the guard instead of the
   operator.

## Claim 2 — "restriction and lift are ONE object, and `Frame` is its shape"

**Half confirmed, half refuted. The refuted half is the one that decides the API.**

See **Q3** below for the full adjudication. Summary: the *pair* structure is right;
`FrameBase` is the right shape for the **basis** side and the **wrong base** for the
**measure** side. And the claim's closing sentence — *"on the G-invariant subspace
both composites are the identity, so it is an isomorphism"* — is **false as stated
on the full space and true only after restricting the domain**, which is exactly the
distinction the type must carry (Frame D, first test).

## Claim 3 — "the two options dissolve"

**(a) fold-the-algorithm vs fold-the-state: CONFIRMED, with a stronger reason than
given, and one cost not named.**

The user's reason is "the state lives on the quotient and the lift is applied at the
seam". The structural reason is sharper: **(a1) makes the physical invariant
`ψ(η,ξ) = ψ(η,−ξ)` an invariant maintained by a hand-written mirror-scatter
(`cache.py` row, `[M† map §3a]`), i.e. Smell #16 Shape 2 — one quantity in `|G|`
representations bridged by an index copy.** The entire point of the fix is to make
the asymmetry *unspellable*; (a1) makes it *re-established after every sweep*. So
there is no design freedom. Confirmed.

**The unnamed cost:** the lift is **not** free at the seam, because `R∘M = Π ≠ I`.
Any full-sphere datum entering the folded solve must be **decomposed**, and its
`V₋` component must **raise**, not be silently dropped. `Π` applied silently is
Claim 1 Refutation 2's landmine in its data-side form.

**(b) fold-existing-nodes vs Hébert-midpoint: CONFIRMED, and the group theory
supplies the missing structural fact.**

Both candidates are quotients — **of two different covers**. The discriminator is
whether the cover's node set **meets the singular locus**:

| Cover | Meets `ξ = 0`? | Orbit types | Quotient measure | `τ_raw[0]` |
|---|---|---|---|---|
| (A) fold existing (`φ` left-endpoint grid, even `n_φ`) | **yes** (`φ = 0, π`) | 2 fixed + free pairs | **orbifold** — `w` on fixed, `2w` on free | `0` `[M† map §2a]` |
| (B) Hébert midpoint (`φ` offset by `Δφ/2`) | **no** | all free | uniform `2w` | `0.2195` `[M† map §2a]` |

`[T]` The weight rule is **orbit-stabilizer**: `W_orbit = w · |G|/|Stab|`. Fixed
points have `Stab = G` ⟹ `W = w`; free points have `Stab = {e}` ⟹ `W = 2w`. The
trapezoid pattern `[w, 2w, …, 2w, w]` measured for (A) `[M† map §2a]` is therefore
**not a quadrature choice — it is the orbifold measure**, forced. `[T ✓ M†]`

**And the R12a trichotomy is a function of this, not of the "rule family."**
`τ_raw[0] = 0` ⟺ the level's most-inward node **is** the level edge ⟺ the cover meets
the singular locus. That upgrades the map doc's D2 finding ("the trichotomy table
presents covering-artefacts as rule properties") from an observation to a theorem
with a named cause. Choosing between (A) and (B) is choosing a **cover**, which is a
quadrature-accuracy question. Confirmed exactly as claimed.

## Claim 4 — "the SH moment break is a domain check, not a cost"

**CONFIRMED, and it is literally the trivial/sign isotypic decomposition. Three
independent numerical confirmations, all derivable by hand.**

`[T]` For `Z₂ = ⟨σ_y⟩`, `L²(S²) = V₊ ⊕ V₋` with projectors `P± = (1 ± σ_y^*)/2`.
`V₊` = ξ-even = the trivial isotypic component; `V₋` = ξ-odd = the sign isotypic
component. The folded measure is a valid cubature **on `V₊` only**; on `V₋` the
folded pairing is not "inaccurate", it evaluates a functional that is identically
zero on the space the fold lives in. Exactly as claimed.

**Confirmation 1 — orbit count is Burnside.** `#orbits = (1/|G|) Σ_g |Fix(g)| =
(N + F)/2` where `F` = number of `ξ = 0` nodes = `2 per level × n_μ`. `[T]`

| rule | `N` | `F` | Burnside `(N+F)/2` | measured `[M† Claim-1 table]` |
|---|---|---|---|---|
| `product(4,8)` | 32 | 8 | **20** | 20 (8 fixed + 12 pairs) ✓ |
| `product(4,16)` | 64 | 8 | **36** | 36 (8 fixed + 28 pairs) ✓ |

`[T ✓ M†]` Burnside counts orbits; **orbit-stabilizer weights them** — two different
theorems, and the code needs the second one.

**Confirmation 2 — the parity table is a character count.** `[T]` On degree-`ℓ` real
harmonics, `σ_y` sends `φ → −φ`: `cos mφ` even (`ℓ+1` of them), `sin mφ` odd (`ℓ`).
So `χ_ℓ(σ_y) = (ℓ+1) − ℓ = 1`, and `dim V₊^{(ℓ)} = ((2ℓ+1) + 1)/2 = ℓ+1`,
`dim V₋^{(ℓ)} = ℓ`. At `ℓ=1`: **2 even, 1 odd** — exactly the measured
`l=1 s0 EVEN | l=1 s1 EVEN | l=1 s2 ODD` `[M† map §2d]`. `[T ✓ M†]`

**Confirmation 3 — the `+2.94` is a broken Gram, and the theory predicts its
sign.** `[T]` On the folded measure the full SH basis is **no longer orthogonal**:
`⟨Y_0, Y_1^ξ⟩_fold = Σ_k W_k ξ_k` is a sum of non-negative terms with at least one
strictly positive ⟹ **strictly positive**, where the full-sphere value is `0` by
mirror cancellation. The measured `φ_1^ξ : −1.3e-16 → +2.94` `[M† map §2d]` **is that
cross-Gram leaking.** `[T ✓ M†]` This matters operationally: the fold does not merely
make one moment wrong, it **makes `frame.gram` non-diagonal**, which per L-010's
corollary means `GramStructure` becomes a lie and `project` must **raise**, not
silently divide.

**One sharpening the claim does not carry: the surviving-moment count depends on
WHICH subgroup you quotient by, so the trial-space restriction cannot be
hardcoded.** `[T]` `dim V_{A₁}^{(ℓ)}` for the full `C_{2v}`:
`χ_ℓ(e) = 2ℓ+1`, `χ_ℓ(σ_y) = χ_ℓ(σ_z) = 1`, `χ_ℓ(C₂ˣ) = (−1)^ℓ`, so
`dim V_{A₁}^{(ℓ)} = ((2ℓ+1) + 1 + 1 + (−1)^ℓ)/4`. At `ℓ = 1` that is **1**, versus
**2** for the `⟨σ_y⟩`-only fold. So a 4-folded rule and a 2-folded rule on the *same
geometry* have **different trial spaces** (the 4-fold additionally kills the axial
current `φ_1^{μ_z}`, which is physically zero for a 1-D cylinder). The restriction
must be **derived from the quotient group carried by the measure**, not written into
the basis.

## Claim 5 — "the real machinery gap is `LevelStructure`"

**Partially confirmed, and MIS-AIMED. Two corrections, one of which flips the
target.**

**Correction 1 — a level is a FIBER OF AN INVARIANT, not an orbit.** `[T]` A level is
the preimage of a value of `Ω ↦ μ_z` (for `product`) or `Ω ↦ |μ_z|` (for `LS`,
`rules_sphere.py:213`) `[C]`. For the *continuous* `SO(2)_ẑ` the fiber **is** the
orbit, so the claim is true for `product` under `SO(2)`. It is **false for `LS`**:
`|μ_z|`'s fiber is **two** `SO(2)`-orbits. And it is false for the *discrete* rule
in general: a fiber is a union of `C_{n_φ}`-orbits (one for `product`; several
`O_h`-orbits for `LS`).

**This distinction is load-bearing.** If "level = orbit", folding a level looks like
quotienting by `SO(2)`. It is not: `SO(2)_ẑ` is **not in the isotropy group** — it
moves the spatial point. The ξ-fold quotients by `⟨σ_y⟩`, which acts **within** the
fiber. Two different groups doing two different jobs (Frame D).

**Correction 2 — the primitive is one layer below `LevelStructure`.** You cannot
fold without the **group ACTION on nodes**: `g·Ω_i = Ω_{σ_g(i)}`. That permutation
is computed three times in the tree (Smell #16 table above) and discarded twice.
`LevelStructure` gaining fields is a **consequence**; the action is the **cause**.
The elegance move, in dependency order:

```
SubgroupOfO3.action(measure) -> NodePermutation | None      # THE missing primitive
  ├── is_invariant(measure)      ≡ action(g) is not None ∀g  # today: a bool that
  │                                                          #   discards the perm
  ├── Quadrature.reflection_index(axis) ≡ action(σ_axis)     # today: a 2nd impl
  ├── MOCMesh._reflected_azi_index      ≡ action(σ)          # today: a 3rd impl
  ├── orbits(action) -> partition                            # ≡ partition_by, exists
  └── measure.quotient(G) = pushforward(rep).consolidate()   # `consolidate` MISSING
```

**Claim 5 is nevertheless right that `LevelStructure` is defective — for a reason it
does not give.** Its two producers fiber over **different invariants** in **one
type** (Smell #16 Shape 2 above). That is a bigger defect than "bare index arrays
with a docstring convention", and it is the direct cause of the two different
`τ_raw` fingerprints (`[0,1,0,1,…]` vs `[1,½,½,0]`) `[M† alpha F3]`.

---

# ANSWER TO Q2 — THE SUBGROUP CHAIN, AND WHETHER IT COMPOSES

## The groups, in the order they enter

`[T]` Write the phase space as `(x, Ω) ∈ ℝ³ × S²`, acted on by an isometry `g` as
`(gx, R_g Ω)`. Three groups act, and **they enter at different stages and do
different jobs**:

| # | Group | Acts on | Reduces | Price paid on the fiber |
|---|---|---|---|---|
| 1 | `ℝ_z` (axial translation) | base only (`R_g = I`) | 3-D → 2-D | **none** — flat connection |
| 2 | `SO(2)_ẑ` (azimuthal rotation) | base **and** fiber | 2-D → 1-D | **the α-redistribution term** |
| 3 | `H = ⟨σ_y⟩ × ⟨σ_z⟩ = C_{2v}` (isotropy at `(r,0,0)`) | fiber only | fiber `S² → S²/H` | **none** — this is the free lunch |

Group 2's stabilizer at the representative point `(r,0,0)` is `⟨σ_y⟩`; group 1's
reflection partner `σ_z` supplies the second generator. Their product is
`C₂(x̂)`, so `H ≅ Z₂×Z₂` with the **radial** direction as principal axis.

## Does the chain compose? Three separate answers.

**(i) Inside the fiber: YES, a genuine tower, and the two folds COMMUTE.** `[T]`
`{e} ◁ ⟨σ_y⟩ ◁ C_{2v}`. `C_{2v}` is abelian, both subgroups are normal, so
`(S²/⟨σ_y⟩)/⟨σ_z⟩ ≅ S²/C_{2v}` — quotienting in either order gives the same
quadrant `{ξ ≥ 0, μ_z ≥ 0}`. So **`level_symmetric`'s 4-to-1 needs no new object.
It is the SAME object with `G = C_{2v}` instead of `G = ⟨σ_y⟩`.** Direct answer to
the question asked.

**(ii) Between the level partition and the fold: NO — they are not two rungs of one
ladder.** `[T]` The level partition is the fibration by the **invariant** of group 2
(`μ_z`, which is the Casimir conserved by the α-redistribution — the redistribution
is a rotation *about ẑ*, so it moves `ω` and fixes `μ_z`). The fold is a quotient by
group 3. Group 2 ⊄ group 3 and group 3 ⊄ group 2 (`σ_y ∈ O(2)_ẑ` but `SO(2)_ẑ ⊄ H`,
and `σ_z ∉ O(2)_ẑ`). They compose as a **fibered** structure, not a tower:

> reduced angular space = the fiber `S²`, foliated by the group-2 invariant `μ_z`
> into leaves, each leaf then quotiented by the group-3 stabilizer.

**Claim 5's "levels and half-ranges are the same operation at different subgroups of
a chain" is REFUTED in exactly this sense.** The right statement is: *a level is a
leaf of the foliation by the connection's conserved coordinate; a half-range is the
isotropy quotient of a leaf.* Two objects, not one.

**(iii) Between the ANGULAR fold and the SPATIAL `r`-fold: they do NOT compose
independently — they are linked by `σ_x σ_y = C₂(ẑ)`.** `[T]` This is the sharpest
finding in this section.

* The `r = 0` centreline is the deck transformation of the **spatial** quotient
  (a diameter folded to a radius). Physically, a particle crossing the axis has both
  `η → −η` **and** `ξ → −ξ`: the deck transformation on directions is `C₂(ẑ)`.
* The code uses `reflection_index("x")` = `σ_x` = `(η → −η, ξ → +ξ)`
  `[M† map §3a]`, which differs from `C₂(ẑ)` on **56 of 64 ordinates** at
  `product(4,16)`, the 8 agreeing ones being exactly the `ξ = 0` self-mirror
  pairs `[M† map §3a]`.
* `[T]` `C₂(ẑ) = σ_x ∘ σ_y` exactly. **On the quotient by `⟨σ_y⟩`, `σ_y ≡ id`, so
  `C₂(ẑ) ≡ σ_x`.**

⟹ **The tree uses the right permutation for the wrong reason, and the coincidence
holds only modulo the fold it has not taken.** The centreline map is well-defined on
directions **only after** the angular quotient; today it is defined on the cover by
an operation that is not the physical continuation. The `[M†]`-measured 56/64
disagreement is that mismatch made visible, and the map doc's observation that it
"dissolves" under the fold is now a **theorem**, not an empirical convenience.

**Consequence for the machinery:** the spatial deck transformation `G` (boundary
campaign B3) and the angular deck transformation `σ_y` are **not independent
factors** — the spatial one is only well-typed as a map `S²/⟨σ_y⟩ → S²/⟨σ_y⟩`. The
quotient object therefore has to be visible to the boundary layer, not private to
the quadrature.

## Why the redistribution exists at all — a one-line criterion

`[T]` **A geometry has an angular redistribution term iff its spatial-reduction
group acts non-trivially on the direction fiber.**

| geometry | reduction group | acts on fiber? | redistribution | fiber quotient `S²/H` |
|---|---|---|---|---|
| slab | translations `ℝ²` | trivially | **none** | `[−1,1]`, `H = O(2)_x̂` (continuous ⟹ collapses) |
| 2-D Cartesian | translation `ℝ_z` | trivially | **none** | hemisphere, `H = ⟨σ_z⟩` |
| sphere | `SO(3)` | non-trivially | `α` in `μ` | `[−1,1]`, `H = O(2)_r̂` |
| **cylinder** | `SO(2)_ẑ × ℝ_z` | non-trivially | `α` in `ω` | **2-D quadrant**, `H = C_{2v}` (finite) |

This single criterion explains three things the corpus currently states separately:
why the slab has no `α` (translation ⟹ flat connection); why the sphere has exactly
one level (`H` continuous ⟹ the fiber quotient is 1-D and there is nothing left to
label); and why the cylinder needs levels at all (`H` **finite** ⟹ the fiber
quotient is **2-D** ⟹ one coordinate is redistributed and the other is a conserved
label). **`n_levels > 1` ⟺ the isotropy group is finite.** That is the honest
definition of `LevelStructure`.

---

# ANSWER TO Q3 — IS `Frame` THE RIGHT BASE, OR ONLY THE RIGHT SHAPE?

**Only the right SHAPE, and the fit is superficial in the direction that matters.
The correct move is to SPLIT the object in two — one half genuinely is a `Frame`,
the other half genuinely is not.**

## The half that is NOT a Frame: the measure-side quotient

Four discriminators, in increasing order of decisiveness.

**D1 — codomain category.** A `FrameBase` maps a **measure space → a coefficient
space**: `analysis: measure_space → test_space`, where `test_space = test.space` is a
`Basis`'s `FunctionSpace` (`frame.py:172-180`) `[C]`. The quotient's output must be
**another `DiscreteMeasure`** — nodes = orbit representatives, weights = orbit sums,
support = `S²/G` — because every downstream consumer (the α recursion, `τ`, the level
march, the *next* fold) integrates against those weights as a **quadrature**.
`FrameBase` has no machinery to emit a measure.

> **The strongest counter-argument to this, stated honestly, and why it does not
> carry.** The fold *is* expressible as a frame: an `IndicatorBasis` over the orbit
> partition, Galerkin, whose analysis `(Mf)_k = Σ_{i∈O_k} w_i f_i` is the orbit
> integral and whose folded weights are literally `W = M𝟙` — the analysis of the
> constant function `[T]`. And the project's homogenisation frame *already* maps a
> fine field to what is effectively another mesh's nodal space, so "codomain must be
> a coefficient space" is weaker than it looks. **Expressible is not the same as
> homed.** The difference: homogenisation's coarse *volumes* come from the coarse
> `Mesh`, not from the frame — the frame emits values. Here the emitted **weights
> ARE the reduced space's quadrature rule**, consumed by the α recursion and the τ
> closure as a measure, and the reduced space's `support` (`S²/G`) is a *new space*
> the frame would also have to fabricate. So D1 survives, but on the narrower and
> more defensible ground: *a frame emits coefficients over an existing space; the
> quotient must emit a new measure over a new space.* D2–D4 are the decisive ones.

**D2 — the verb set does not transfer.** `frame.conjugate(A) = R∘A∘M` is a
change of basis, defined for **any** `A`. Descending an operator to a quotient is
defined **only for `G`-equivariant `A`**, and for non-equivariant `A` the same
formula silently returns the **group average** `(1/|G|)Σ_g gAg⁻¹`. A frame has no
precondition to check; the quotient's entire correctness is one. Putting the
quotient under `FrameBase` inherits a verb whose contract it cannot honour.

**D3 — the existing frame-native verb already exists and already refuses this
case.** `DiscreteMeasure.pushforward` (`measure.py:547-618`) `[C]` is the
measure-category morphism. Its docstring says, verbatim: *"For a non-invertible `φ`,
`φ_*μ` may have atoms collapsing onto the same target point — the caller's `g` will
then see those weights summed implicitly through the integration step.
**Mathematically valid; numerically the node array will contain duplicates with
separate weights.**"* (`:573-577`) `[C]`

So `pushforward(orbit_representative_map)` **already produces the quotient measure —
in an unreduced representation with duplicate atoms.** The missing verb is not a
frame; it is **`consolidate()`**: merge coincident atoms, summing weights. That is
one method, and it makes `quotient(G) = pushforward(rep_map).consolidate()`. The
codebase has already written the hard half and documented the missing half as a
known gap. **This is a far smaller and far more honest object than a Frame
subclass.**

**D4 — the ownership discriminator (project memory: "an operator owns its frame IFF
the frame is its eigenbasis").** Apply it: **no operator owns the quotient.** The
group is in the **commutant of every equivariant operator at once** — streaming,
collision, scattering, fission, and the BCs are all `C_{2v}`-equivariant. By Schur,
they are *all* simultaneously block-diagonal in the isotypic decomposition. A frame
owned by everything is owned by nothing: the quotient is a property of the
**symmetry**, i.e. of the *problem*, not of an operator. Under the project's own
discriminator, it fails the frame-ownership test.

## The half that IS a Frame: the basis-side isotypic restriction

`[T]` The trial-space restriction to `V₊` (or `V_{A₁}`) **is** a genuine `Basis`:
**symmetry-adapted spherical harmonics** (Altmann; Bradley & Cracknell — a standard,
literature-backed object in crystallography and computational chemistry). It has
exactly the `Basis` contract: `evaluate(nodes) → table` with `dim V_{A₁}` columns
instead of `(L+1)²`. Bound to the folded measure it is a **`GalerkinFrame`** (test
IS trial), and its Gram is the full-sphere Gram restricted to the invariant block:

`[T]` for `Y_a, Y_b ∈ V₊`, `Σ_k W_k Y_a(Ω_k)Y_b(Ω_k) = Σ_i w_i Y_a(Ω_i)Y_b(Ω_i)`
because invariant functions are constant on orbits and `W_k = Σ_{i∈orbit} w_i`.
**Exactly diagonal, exactly the parent value.** So the fold does not degrade the
frame at all — it deletes the columns that were structurally zero.

**Precise verdict:**

```
NOT a Frame:  measure.quotient(G)  ->  DiscreteMeasure       # measure category
              (= pushforward ∘ consolidate; the MISSING verb is consolidate)

IS  a Frame:  SymmetryAdaptedHarmonicBasis(G, L)  ->  Basis
              GalerkinFrame(that_basis, folded_measure)      # unchanged hierarchy
```

The two are related by **one law**, and that law is the thing worth typing:
`analysis` through the symmetry-adapted frame on the **folded** measure equals
`analysis` through the plain SH frame on the **full** measure, restricted to the
invariant slots. That is a fail-able equivalence, not an abstraction.

## The site where this must land, and why it structurally cannot today

`Quadrature.angular_frame(L)` (`directional.py:449-458`) `[C]`:

```python
s2_measure = DiscreteMeasure(
    nodes=np.column_stack([...]), weights=self.weights, support=SPACE_SPHERE,
)
return GalerkinFrame(SphericalHarmonicBasis(L=L), s2_measure)
```

Two structural facts:

1. It **rebuilds** the measure and **drops `invariance_group` entirely** (the kwarg
   is not passed) `[C]`. The frame's measure forgets the symmetry, so the frame
   *cannot* consult it even in principle.
2. It **hardcodes the full `SphericalHarmonicBasis`** `[C]`. On a folded measure
   that is a silently-wrong frame with a non-diagonal Gram (Claim 4, Confirmation 3).

The docstring's own naming tripwire (`:439-447`) `[C]` anticipates *"a second
angular basis … parametrises the basis — `angular_frame(basis=...)`, a SIGNATURE
change."* **The symmetry-adapted basis is that second basis, and it arrives via the
measure's group rather than via a caller argument** — so the parametrisation the
docstring predicts should be **derived**, not passed.

---

# FRAME CANDIDATES

## Frame A — Orbifolds / quotient by a finite group action (A.1 topology row)

**Trigger.** Feature 7 (isotropy subgroup `C_{2v}` acting on the fiber) + Feature 14
(zero-flux closure at both ends of the reduced interval) + Smell #12 (six ε-tests
detecting one geometric set). Two of the three are code facts `[C]`.

**Reformulation.** The reduced angular domain is the **orbifold** `S²/H`, not a
manifold. Concretely:

* `[T]` **Points**: orbits. `#orbits = (N + F)/2` (Burnside), `F = |{ξ = 0}| = 2n_μ`.
* `[T]` **Weights**: orbit-stabilizer, `W = w·|G|/|Stab|` — `w` on the singular
  locus, `2w` off it. This **derives** construction (A)'s trapezoid weights.
* `[T]` **Singular locus**: `Σ = {ξ = 0} = Fix(σ_y)`, a first-class *object*, not an
  ε-test. It is simultaneously: the `φ ∈ {0,π}` nodes, the self-partners of
  `reflection_index('y')`, the `|η| ≈ 0` degenerate positions, the level-edge
  endpoints, the `α = 0` closure points, and the half-weight nodes. **Six sites,
  one set** `[C]`.
* `[T]` **Boundary condition**: the α closure `α_{1/2} = α_{M+1/2} = 0` is the
  **zero-flux (reflective) condition at the orbifold boundary** — no current can
  cross a mirror plane. Identical in kind to a spatial reflective BC. And `α = −Wξ`
  vanishes on `Σ` **automatically**, because `Σ` is the mirror plane, so the closure
  is a *theorem about the quotient*, not a modelling convention.

**Elegance payoff.**

* *Structure-exposing*: the α-dome's two-sided zero closure — currently a cited
  convention (BMC Eq. 50) with an unexplained asymmetry (the cylinder path has **no**
  closure assert where the sphere does, `reduced_operator.py:778-786` vs `:688-695`
  `[C]`) — becomes a **derived** consequence of the mirror plane. It also names the
  singular locus once instead of six times.
* *Structurally-simpler*: retires six ε-tests and the `_edge_seed_stencil` tie-skip
  loop (`pole_angular_closure.py:1012-1030`, which the map doc shows becomes dead
  code `[M† map §3a]`).
* *Algorithmic-advantage*: `|G|`-fold reduction in angular unknowns
  (`1.6×` at `n_φ=8`, `1.88×` at 32 for the `⟨σ_y⟩` fold `[M† map §7]`; `4×`
  asymptotically for the full `C_{2v}`).

**First test (discriminates).** Build the fold with **uniform** doubled weights (the
naive "halve the nodes, double the weights" fold — the implementation this frame
claims is wrong) and assert `Σ W_k == Σ w_i` and `∫Y_0² == 4π`. **The naive fold
FAILS both** — it over-counts the two `ξ=0` nodes per level by `w` each, so
`ΣW = Σw + 2n_μ·w`. The orbifold weight rule PASSES. This is fail-able, targets the
specific dropped term (the stabilizer order), and cannot be passed by an
implementation that ignores fixed points.

**Structural attack on current.** The cylinder path has no closure assertion at all
(`reduced_operator.py:778-786` `[C]`), and the sphere's assertion is
telescoping-invariant `[M† alpha §4]`. Under this frame the closure is not an
assertion to add — it is a **property of the quotient's boundary** which a correct
fold cannot violate. The current formulation cannot state that, because it has no
object representing the boundary.

**Precedent.** Extends `scripts/validated_unified_geometry.md` (slab/annulus/hollow
sphere as one manifold-with-boundary) from the **base** to the **fiber**. Same frame,
new axis.

## Frame B — Representation theory / isotypic decomposition (A.2 group-theory row)

**Trigger.** Feature 7 + the measured parity table `[M† map §2d]` + Smell #4. The
question "is the ξ-even/ξ-odd split literally the trivial/sign isotypic
decomposition of `Z₂`?" is answered **yes, literally** `[T]`.

**Reformulation.** `L²(S²) = ⊕_χ V_χ` over irreps of `H`. For `H = ⟨σ_y⟩`: two
1-D irreps (trivial, sign) ⟹ `V₊ ⊕ V₋`. For `H = C_{2v}` (abelian, order 4): four
1-D irreps `A₁, A₂, B₁, B₂` ⟹ a 4-way parity split by `(ξ-parity, μ_z-parity)`. The
physical solution lives in `A₁`. Every `H`-equivariant operator is
**block-diagonal** across isotypic components (Schur). The fold is the projection
`P_{A₁}`; the folded measure is `A₁`'s representation in orbit coordinates.

Dimension counts, both derived and both checkable `[T]`:
`dim V₊^{(ℓ)} = ℓ+1` (⟨σ_y⟩); `dim V_{A₁}^{(ℓ)} = ((2ℓ+1)+1+1+(−1)^ℓ)/4` (`C_{2v}`).

**Elegance payoff.**

* *Structure-exposing*: the reported "break" is the **sign-isotypic component**, and
  the code's `-1.3e-16` for `φ_1^ξ` today is not a numerical accident — it is the
  statement `V₋ ⊥ V₊` in the full-sphere Gram. The `+2.94` is that orthogonality
  destroyed. Both numbers get a cause.
* *Expressive*: replaces "which `(ℓ,m)` slots survive?" (a per-case lookup) with a
  character formula that covers every subgroup uniformly, including
  `level_symmetric`'s 4-fold **with no new machinery**.
* *Algorithmic-advantage*: `dim V_{A₁} ≈ (L+1)²/|H|` — the moment vector shrinks by
  the group order too, not just the ordinate vector. And Schur block-diagonality is
  exactly the statement that the scattering operator never couples `V₊` to `V₋`, so
  no cross-block work is ever done.

**First test (discriminates).** Two halves, and the negative half is the real one.
*Positive*: for a `V₊` field, `folded_frame.analysis(f) == full_frame.analysis(f)[V₊
slots]` to round-off. *Negative*: on the folded measure assert
`⟨Y_0, Y_1^ξ⟩ = Σ_k W_k ξ_k > 0` while the full-sphere value is `0` — i.e. **assert
that `frame.gram` is NOT diagonal for the un-restricted basis**, and therefore that
constructing `GalerkinFrame(SphericalHarmonicBasis(L), folded_measure)` **RAISES**.
An implementation that folds the measure and leaves the basis alone PASSES nothing
here; it is exactly what produces the `+2.94`.

**Structural attack on current.** `angular_frame` builds its measure **without**
`invariance_group` (`directional.py:451-457` `[C]`), so the symmetry is erased at the
one site that would need it. The frame exposes that the truncation set `(ℓ,m)` is
being chosen by *degree alone* when the correct index set is *(degree, irrep)* — a
missing axis on the moment index, not a missing check.

## Frame C — Deck transformations / covering spaces, already built in this tree for the SPATIAL half (A.1 topology row)

**Trigger.** Feature 5 `[C]` — the boundary module states the *exact* concept, with a
membership test, in the spatial phase space; Feature 8 `[T]` — the group relation
linking the two; Smell #16 Shape 4 (a third path about to be written for machinery
that exists).

**Reformulation.** `orpheus/geometry/boundary/__init__.py:19-46` `[C]` already
carries, verbatim: *"`G : Γ₊ → Γ₋` is the **deck transformation** … the composition
operator of a measure-preserving bijection of the boundary phase space"*, plus the
membership criterion *"`G` is the deck transformation of an **actual quotient of the
physical domain**"* and the theorem *"**exactly one of** `G`, `R` **is
non-trivial**"*.

The `ξ`-mirror is a deck transformation of an actual quotient of the **angular**
phase space, and it satisfies the *same* two-part test:

* multiplicativity: `σ_y(ψφ) = (σ_yψ)(σ_yφ)` ✓ (it is a relabeling) `[T]`;
* actual-quotient: `S²/⟨σ_y⟩` is a genuine quotient of the direction sphere ✓ `[T]`.

⟹ by the campaign's own theorem, **`G` non-trivial ⟹ `R` trivial**: the angular
fold's closure carries **no constitutive kernel**. That is precisely `α = 0` — a pure
zero-flux geometric condition, no albedo, no re-emission. The theorem *predicts* the
BMC closure.

**Elegance payoff.**

* *Structure-exposing*: makes the α closure a **boundary law of the angular
  quotient** and puts it in the same algebra as the spatial reflective BC — one
  concept, one vocabulary, across two phase-space factors.
* *Structurally-simpler*: the B3 membership test (multiplicative ∧ actual-quotient)
  transfers verbatim; no second criterion needs inventing. The `G.domain_face`
  discipline (which `Γ₊` a law consumes) has an exact angular analogue (which leaf's
  `ω`-endpoint).
* *Algorithmic-advantage*: none directly. Two criteria, not three — a moderate match,
  but the *reuse* argument is strong.

**First test (discriminates).** Assert `C₂(ẑ) ∘ σ_y = σ_x` **as node permutations**,
bit-identically, on `product(4,16)` — i.e. `perm(C2z)[perm(sigma_y)] == perm(sigma_x)`
with `array_equal` (0 ULP, per L-002). Then assert that **before** the fold,
`perm(C2z) != perm(sigma_x)` on **56 of 64** ordinates, agreeing on exactly the 8
`ξ = 0` self-partners `[M† map §3a]` — and that **after** the fold they are
`array_equal`. An implementation that keeps using `σ_x` on the unfolded cover cannot
pass the first half; one that folds but keeps `σ_x` un-re-derived cannot explain the
second.

**Structural attack on current.** The centreline map `reflection_index("x")` is the
**physical** deck transformation `C₂(ẑ)` only modulo `σ_y` `[T]`. The current code is
right by a coincidence whose precondition (the fold) it has not taken. That is not a
style point: it means the `r = 0` continuation and the angular fold are **one
decision**, and splitting them across `geometry/boundary/` and
`numerics/quadrature/` guarantees they can drift.

## Frame D — Associated fiber bundle with a connection (A.1 differential-geometry row, sharpened)

**Trigger.** Feature 6 + Feature 10 (a first-order transport term generated purely by
the reduction, not by the physics) + Smell #3. This frame is already validated in
this agent's library for the Variant-α family; what is **new** here is the identity
of the fiber.

**Reformulation.** The reduced 1-D cylindrical phase space is the **associated
bundle** `P ×_H (S²/H)` where `P → [0,R]` is the principal `O(2)_ẑ`-bundle over the
radius and `H = ⟨σ_y⟩` is the isotropy. Three consequences, all `[T]`:

1. **`α` is the connection**, i.e. the derivative of the moving frame `(r̂, φ̂)` — the
   `−(1/r)∂(ξψ)/∂φ` term is parallel transport of the ordinate direction, exactly the
   validated `candidate_cylindrical_connection.md` match, now with its *existence
   criterion*: **a reduction has a connection term iff its group acts non-trivially
   on the fiber** (Q2 table).
2. **The typical fiber is `S²/H`, not `S²`.** This is the new content: the bundle's
   fiber is *already* the quotient. Modelling the fiber as the full `S²` is modelling
   the **`|H|`-fold cover** of the physical fiber.
3. **`μ_z` is the connection's conserved coordinate** (the redistribution is a
   rotation about `ẑ`), so the fiber foliates into leaves labelled by `μ_z`. **A
   level is a leaf.** `n_levels > 1` ⟺ `H` is finite ⟺ the quotient fiber has
   dimension ≥ 2.

**Elegance payoff.**

* *Structure-exposing*: unifies slab / sphere / cylinder / 2-D Cartesian under one
  criterion (Q2 table), and turns `level_structure is None` from a null-branch into
  the **`H` continuous** case. It answers "why does the cylinder have levels and the
  sphere not?" from the group, not from the code.
* *Expressive*: `LevelStructure` gets an honest signature — `(conserved coordinate,
  its value per leaf, the leaf's ordinate chain in **redistribution** order, the
  leaf's **quotient** weights, the leaf's singular set)`. Every field is derived, none
  conventional. The docstring convention that #326 turns on (`sorted by increasing η`)
  becomes **`ordered along the connection`**, which is `ω` — the variable the α
  recursion actually integrates `[M† alpha §3]` — and the ordering ambiguity
  dissolves because `η` is injective on the quotient leaf `[M† map §2a]`.
* *Algorithmic-advantage*: none new beyond Frame A's.

**First test (discriminates).** Ordering blindness. Assert that the per-level
ordinate chain is **the sort along the connection coordinate `ω`**, and that this
is `argsort`-tie-break-independent. On the **cover** the sort key `η` has exact ties
(`min|Δη| = 0` `[M† map §2a]`) and the answer moves by 2e-3–3.7e-3 on the mirror
criterion `[M† alpha F3]`; on the **quotient leaf** `min|Δη| = 0.239` at `n_φ = 8`
`[M† map §2a]` and no tie exists. So: assert `sort(key=η)` and `sort(key=ω)` produce
`array_equal` index chains. **On the cover this FAILS** (they interleave —
`φ/π : 1.00, 1.25, 0.75, 1.50, …` `[M† alpha F2]`); on the quotient it passes. An
implementation that folds the state but keeps sorting the cover cannot pass.

**Structural attack on current.** `level_indices` is documented as *"sorted by
increasing `η`"* (`rules_sphere.py:126-129` `[C]`, and
`curvilinear_one_group.rst:229-232` `[M† alpha §1b]`) — a sort by a **coordinate**
where the math requires a march along the **connection**. The two coincide only where
`η` is monotone in `ω`, i.e. **only on the quotient** `[T]`. The current formulation
has no object for "the connection's coordinate", so it cannot state the requirement
it needs; it states a proxy that is true on a half-circle and false on a circle.

## Frame E — REFUTED: the quotient as a `PetrovGalerkinFrame`

**Trigger checked**: the map doc proposes *"a folded analysis face with an unfolded
reconstruction face is exactly `PetrovGalerkinFrame`; the symmetrized basis restores
`GalerkinFrame`"* `[M† map, HEADLINE]`.

**Refuted, and the reason is structural, not stylistic.** `[T]` Petrov-Galerkin means
**test ≠ trial** — and in this project's ratified reading, test ≠ trial is a
**solution weighting** (flux-weighted homogenisation, spectrum-weighted
condensation), which is why the discipline is a type (`frame.py:41-56` `[C]`). A
symmetry fold carries **no solution weight**: `test` and `trial` are both the
symmetry-adapted harmonics, and the Gram stays diagonal at exactly its parent value
(Q3). Calling the fold Petrov-Galerkin would put a **pure geometric identification**
in the slot reserved for **solution-dependent weighting** — the same category error
the B3.0 campaign fixed when it moved the Lambertian out of the geometry slot
(`white.py:118-124` `[C]`).

The fold is **Galerkin on a smaller space**, exactly as the DSA `ℓ=0` sub-block frame
is Galerkin on a smaller space (project memory: `angular_frame(0)`, Π = P∘M
W-self-adjoint under the **plain** measure, NOT PG). Same verdict, same reason,
second sighting. The map doc's "unfolded reconstruction face" is not a Petrov-Galerkin
test space — it is the **lift `R` of the quotient**, which lives in the measure
category (Q3, D1), not on the frame's test side.

---

# WHAT IS STILL MISSING — none of the five claims names these

Ordered by value. M1 is the one I would act on first.

## M1 — The orbifold SINGULAR SET is a first-class object with six ad-hoc detectors

Every reported pathology in #326 lives on `Σ = Fix(σ_y) = {ξ = 0}`, and no code
names it. Enumerated `[C]` unless marked:

| Site | Spelling | What it is really testing |
|---|---|---|
| `augmented_mesh.py:815` | `eps < tau_level[0] < 1.0 − eps` (R12a) | is the leaf's first node **on** `Σ`? |
| `pole_angular_closure.py:1012` | `abs(mu[cand] − mu[m0]) > 1e-14` | skip the mirror partner (a **free orbit** collapsed on the cover) |
| `loss_representation/__init__.py:2902` | `|η| < 1e-15` degenerate positions | the `φ = π/2` nodes — `Σ` of the *other* mirror |
| `rules_sphere.py:212` | `tol = 1e-12` level membership | leaf membership under `|μ_z|` |
| `directional.py:477` | `|μ_axis| < 1e-15` ⟹ octant label `0` | **`Σ` for all three mirrors at once** — the tree already has a partition entry for it |
| `symmetry.py:952` | `atol` weight match | orbit-partner weight equality |

`directional.py:464-482` is the tell: `octants` **already** partitions by sign with a
zero-label class for the on-plane nodes `[C]`. That zero-label class **is** the
singular locus, already computed, already a `DiscreteMeasurePartition`, never
consumed as such. The elegance move is to name `Σ` once (as
`orbit_type(node) ∈ {free, fixed}`) and let all six sites read it.

*Why this beats the five claims:* it is the smallest change with the largest reach,
it needs no new type hierarchy, and it is the direct cause of the R12a
round-off-decides-a-data-layout defect `[M† alpha §5a]` — a **floating-point
comparison deciding whether a level allocates a first-class state block**.

## M2 — `consolidate()`: the one missing measure verb

`pushforward` exists and its docstring documents the atom-merging case as valid but
unreduced (`measure.py:573-577` `[C]`). `quotient(G) = pushforward(rep).consolidate()`.
One method. No Frame subclass, no new hierarchy. **This is the whole of "the
machinery has the checking face and lacks the quotient" — and it is one verb, not a
class.**

## M3 — `degree_of_exactness` has no space attached

`restrict` and `pushforward` both **drop** `degree_of_exactness` (`measure.py:616-617,
:652-654, :667`) `[C]`. Dropping is honest but throws away a real, provable claim:
the folded rule is **exact at the parent degree on `V₊`** and **meaningless on `V₋`**.
A consumer integrating a `V₋` function against the folded measure gets a wrong number
with a `None`-exactness contract that never warned it. The type wants
`(degree, subspace)` — and the subspace is exactly the irrep label from Frame B. Note
this is *not* speculative: the whole `+2.94` failure mode is a consumer integrating a
`V₋` function against a folded measure.

## M4 — `SubgroupOfO3` has containment but no QUOTIENT lattice operation

`contains` / `is_subgroup_of` exist (`symmetry.py:343-392`) `[C]`; there is no
"what group survives the fold". The quotient of a `C_{n_φ}`-invariant rule by
`⟨σ_y⟩` retains **no** azimuthal rotation (the fold breaks it), so the folded measure
must carry `Trivial` — and today `restrict`/`pushforward` set `invariance_group=None`
(`:616`, `:666`) `[C]`, conflating **"trivial"** with **"unknown"**. That conflation
is load-bearing at exactly one place: the *second* fold. Folding by `⟨σ_z⟩` after
`⟨σ_y⟩` requires knowing that `σ_z` survived, and `None` cannot say so. **Without
this, the tower in Q2(i) is not composable in code even though it composes in
math.**

## M5 — There is no operator-side equivariance predicate, so the descent cannot be gated

Claim 1's precondition (b). Nothing in the operator algebra answers "does `A` commute
with `g`?". Without it, the fold is gated on the *measure* while the *correctness* is
a property of the **operator and the data**. The concrete exposure: nothing forbids a
user handing a `ξ`-odd `AngularSourceSink` or a `ξ`-odd prescribed inflow to a
cylinder `[M† map §3a]`, and the folded path would **symmetrize it silently**.
Per project memory, the invariant test for a *name* is the invariant run against the
**violator** — so this owes a negative test: a deliberately `ξ`-odd source must make
the folded solve **RAISE**, not converge to a plausible wrong answer.

## M6 — "half-range" is one word for two operations (Claim 1, Refutation 1)

`restrict`'s docstring already claims the phrase for the slab `μ>0` case
(`measure.py:632-633`) `[C]`. If the fold ships under the same word, every future
reader has to disambiguate by context. Per L-012, the naming decision belongs in the
design, not after it: **fold / quotient** for the mass-preserving orbit collapse;
**half-range** stays with the mass-dropping directional restriction.

## M7 — The `ξ`-fold and the `r = 0` centreline are ONE decision, split across two packages

Q2(iii). `σ_x` is the physical continuation only modulo `σ_y` `[T]`. Today the
angular fold would land in `numerics/quadrature/` and the centreline map lives in
`geometry/boundary/` + `loss_representation/` `[C]`. Nothing links them, so a future
change to either silently breaks the relation `C₂(ẑ) = σ_x σ_y`. The fold must be
**visible to the boundary layer** (Frame C), which is an argument against making it a
private detail of `LevelStructure`.

## M8 — The two `LevelStructure` producers fiber over different invariants

Smell #16 Shape 2, restated as a deliverable: `product` fibers over signed `μ_z`
(`rules_product.py:126-140`), `level_symmetric` over `|μ_z|` (`rules_sphere.py:213`)
`[C]`. One type, two meanings, and every downstream consumer assumes one. **Fix this
before folding anything**, per the clean-before-extend rule — otherwise the fold
grows a third arm on a seam that already has two.

---

# CROSS-METHOD POLLINATION

**Current method class:** SN (curvilinear / cylindrical), angular discretisation +
quadrature machinery.

## Borrowing 1 — from MoC: the azimuthal quotient is ALREADY LIVE, and hand-rolled

**Trigger.** `orpheus/moc/geometry.py:222-229` `[C]`:

```python
def _reflected_azi_index(phi: np.ndarray, azi_index: int) -> int:
    """Both vertical and horizontal reflections map phi -> pi - phi."""
```

Two structural facts, both `[T]` on top of that `[C]`:

1. MoC's azimuthal grid lives on `[0, π)` (the wrap arithmetic at `:226-228` is
   modulo `π`, not `2π`). A track *line* is direction-agnostic, so MoC's angular
   space is **already the quotient `S¹/⟨−1⟩ = RP¹`** — MoC took a quotient and never
   named it.
2. `φ → π − φ` is the deck transformation of a **second, further** quotient (the
   mirror). It is the **third** independent implementation of the node permutation
   (Smell #16 table), with its own `argmin(abs(...))` matcher and its own wrap.

**Reformulation.** `MOCMesh`'s reflection table is `SubgroupOfO3.action(σ)` on an
azimuthal measure — the *same* primitive M1/M2 need. And MoC's `[0,π)` grid is
`azimuthal_measure.quotient(⟨−I⟩)`.

**Payoff.** *Structurally-simpler* (three implementations → one) and
*structure-exposing* (MoC's `[0,π)` convention stops being a ray-tracing folklore
convention and becomes the same object as SN's `ω ∈ [0,π]`). **This is the answer to
"is this SN-specific?" — it is not.** MoC is a second, already-shipping consumer,
and it needs the *same* two things (the action, and the quotient of a 1-D angular
measure).

**First test.** Assert `MOCMesh`'s reflection table equals
`action(σ)` on the MoC azimuthal measure, `array_equal` (0 ULP). Then the negative:
feed an azimuthal grid for which `π − φ` is **not** a node (an odd/asymmetric
grid) and assert the action is **absent** (raises / returns `None`) where today
`argmin(abs(...))` returns the **nearest** index silently. That silent
nearest-match is a live latent defect of exactly ERR-042's kind, one package over.

## Borrowing 2 — from computational chemistry / crystallography (adjacent to SN via A.2)

**Trigger.** Feature 7 + the need for a restricted trial space (Frame B). This is the
**named literature object**: symmetry-adapted spherical harmonics / projection-operator
basis construction (Altmann; Bradley & Cracknell; the same lineage `symmetry.py:59-73`
already cites for Hamermesh and Stiefel-Fässler `[C]`).

**Reformulation.** `SymmetryAdaptedHarmonicBasis(group, L)` built by the standard
projection operator `P_χ = (dim χ/|G|) Σ_g χ(g)* g`. For the abelian `H` here, all
irreps are 1-D and `P_{A₁} = (1/|H|) Σ_h h` — i.e. **the basis columns are just the
even-parity harmonics**, selected by the character formula, not by a hand-written
`(ℓ,m)` skip list.

**Payoff.** *Expressive* — one formula covers `⟨σ_y⟩`, `C_{2v}`, and any future
lattice group (`D_{6h}` for hex, which the `symmetry.py:36-42` docstring already
names as the trigger for extending the module `[C]`). *Structure-exposing* — makes
the moment index `(degree, irrep)` rather than `degree`, which is M3.

**First test.** `dim(SymmetryAdaptedHarmonicBasis(⟨σ_y⟩, L=1)) == 2` and
`dim(SymmetryAdaptedHarmonicBasis(C_2v, L=1)) == 1`, matching the character formula.
Discriminates because a naive implementation that "drops the sine harmonics" gets the
`⟨σ_y⟩` case right and the `C_{2v}` case **wrong** (it would give 2, not 1 — it
misses that `μ_z`-odd must also go).

## Borrowing 3 — from the boundary campaign (SN internal, B3)

Covered as **Frame C**. Trigger: `geometry/boundary/__init__.py:19-46` `[C]`. The
membership test (multiplicative **necessary**, actual-quotient **sufficient**,
⟹ exactly one of `G`, `R` non-trivial) transfers verbatim to the angular fiber and
**predicts** `α = 0` as a pure-geometry closure with no constitutive kernel. First
test as given in Frame C.

## Borrowing 4 — from homogenisation / condensation (SN internal, the landed frame campaign)

**Trigger.** L-010: a collapse's morphism is fixed by **what is conserved**. The fold
conserves `∫_{S²} f dμ` for `f ∈ V₊`.

**Reformulation.** The fold is the **marginalize** morphism `M` on the *weights*
(`W_k = Σ_orbit w_i` — mass-preserving, no `G⁻¹`) composed with the **average**
morphism `G⁻¹M` on the *values* (constant on orbits ⟹ the average is the value). The
two coincide **only** on `V₊` — which is another statement of Claim 4. So the fold
is not "a projection with a weight"; it is the *degenerate* case where marginalize
and average agree, and it agrees **exactly on the space the fold is defined on**.

**Payoff.** *Structure-exposing*: gives the fold the same vocabulary as the landed
`project` / `analysis` split, so the new object does not invent a third collapse
verb. **This is the argument that the fold does NOT need `normalize=` as a flag** —
it needs the invariant-subspace domain restriction, and then the flag is moot.

**First test.** Assert `fold(f).sum_against_weights() == full(f).sum_against_weights()`
for `f ∈ V₊` (mass preserved) **and** that the same assertion **FAILS** for
`f ∈ V₋`, where the full-sphere integral is `0` and the folded one is not. The
failure is the domain statement, and a design that silently returns a number for
`f ∈ V₋` is the thing being refuted.

---

# UNEXPLORED

Frames checked against the feature list that did **not** fire, with the structural
reason each failed. Recorded so the next session does not re-attack them.

* **Homology / chain complex / de Rham** — no `∂² = 0`. The fold's `M`/`R` pair is a
  section-retraction (`R∘M = Π ≠ I`, `M∘R = I`), a **projector** pair, not a
  differential. Standing refutation (L-001), re-confirmed here: `Π` is idempotent,
  not nilpotent.
* **Tensor networks / MPO** — bond dimension is `|H| = 2` or `4`, **fixed and not a
  truncation knob**. Degenerate rank; not a network (L-001). Would fire only if a
  *chain* of quotients with a tunable depth appeared, which the `C_{2v}` lattice
  forbids (it has depth 2, maximum).
* **Category theory / groupoids** — the "action groupoid `S² ⋊ H`" is the formally
  correct object, and it buys **nothing here**: `H` is finite abelian, all orbits are
  free or fully-fixed, and the concrete frames (orbifold + isotypic + deck
  transformation) already carry every lever. Per L-001, name the concrete frame.
* **Crystallographic / Bloch theory** — needs a *lattice translation* group with a
  Brillouin zone. `C_{2v}` is a point group with no translations. Fires for the hex/
  full-core lattice work the `symmetry.py:36-42` docstring anticipates, **not** for
  this.
* **Wiener-Hopf / half-space `H`-functions** — the word "half-range" is a false
  friend. Wrong solver family (Chandrasekhar line), structurally incompatible with a
  sweep (L-001). It is also the *restriction* sense of "half-range" (M6), not the
  quotient sense.
* **Asymptotic / homogenisation (A.6)** — no scale separation. The fold is exact at
  every `n_φ`; there is no small parameter.
* **Approximation theory (quadrature accuracy)** — genuinely fires, but **on
  construction (A)-vs-(B) only**, which Claim 3(b) correctly amputates from the
  symmetry question. It is answerable on quadrature grounds (which cover integrates
  `V₊` better) and must not be entangled with the fold decision.
* **Krylov / spectral (A.5)** — the fold changes `n_dof` and therefore
  `restart`/`to_flat` sizing (the ERR-053 family, flagged `[M† map §3b]`), but that is
  a *consequence* to test, not a frame that explains anything about the fold.
* **Measure theory / Radon-Nikodym (A.4)** — half-fires and is worth one line: the
  quotient's weight rule **is** the Radon-Nikodym derivative of the pushforward
  against the orbit-counting measure, and `pushforward`'s docstring explicitly says
  the Jacobian is **not** applied (`measure.py:565-571` `[C]`). But the discrete
  orbit-stabilizer rule (Frame A) is the concrete, implementable form; the RN framing
  adds vocabulary without a new lever. Recorded because it is the natural place a
  future session will reach when generalising beyond finite groups (a *continuous*
  quotient — the slab's `S²/O(2)` — **does** need the `2π` orbit-volume Jacobian,
  and that is exactly the un-named `W` normalisation constant).

---

# SELF-CORRECTION

Two errors in the first draft of this document, caught and reissued. Recorded
because a frame-attack whose own claims are not stress-tested is worth nothing.

**1 — I attributed ERR-042 to the wrong group element.** The first draft asserted
that ERR-042 (odd `n_φ` rejected on a cylinder) is *"verbatim the statement that the
`ξ`-fold's quotient does not exist"*, i.e. the fold's own existence condition
enforced at the wrong layer. **That is false.** `[T]` The equispaced grid
`{2πk/n_φ}` is closed under `φ → −φ` (the fold's mirror `σ_y`) for **every** `n_φ`
— the `ξ`-fold exists unconditionally. ERR-042 checks `σ_x` (`φ → π − φ`), which
needs even `n_φ`, and by `C₂(ẑ) = σ_x σ_y` that is the **centreline** deck
transformation's existence, not the fold's. Corrected in Claim 1, and the corrected
version is a *stronger* finding: it says the check ORPHEUS makes is on the wrong
element and agrees with the right one only because `σ_y`-closure is automatic.

**2 — I over-stated D1 in the Q3 verdict.** The first draft claimed a frame's
codomain is categorically not a measure. The homogenisation frame is a near-miss
counter-example (its codomain is effectively a coarse mesh's nodal space), and the
fold *is* expressible as an `IndicatorBasis` Galerkin frame. Reissued with the
counter-argument stated in full and D1 narrowed to its defensible form, with the
verdict re-anchored on D2 (unhonourable `conjugate` contract), D3 (`pushforward`
already exists; one verb missing), and D4 (no operator owns it — the group is in
every equivariant operator's commutant).

Trap caught in both cases: asserting a structural claim from a *pattern match* on
vocabulary ("mirror", "measure→measure") without running the group element / the
counter-example. Neither is a hedging failure; both are detection failures, which is
worse in this role.
