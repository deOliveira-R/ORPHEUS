---
name: quadrature-symmetry-quotient-frames
description: Issue #326 — the angular half-range is an ORBIFOLD quotient by the point-isotropy group; a level is a FIBER of an invariant (not an orbit); the measure-side quotient is NOT a Frame but the basis-side isotypic restriction IS; redistribution exists iff the spatial-reduction group acts on the angular fiber.
metadata:
  type: project
---

Frame attack on the proposed #326 "quadrature machinery design" (user rejected an
options-menu framing: *"Should half range be a restriction over full-range? Do we need
some abstraction class for Levels? For symmetry? Map what is missing in machinery
design."*). Deliverable `scratch/issue326_quotient_frames.md`. Grounded on
`refactor/operator-strategy-layers` worktree reads (L-005); no shell this dispatch, so
every `[M]` is cited from the two prior scratch investigations, and the derivations
below were checked against those numbers by hand.

**Why:** the fix ships as machinery, so the durable output is the STRUCTURE, not the
patch. **How to apply:** reach for this on any angular-discretisation, quadrature-
selection, or geometry-reduction question; the last three items generalise well past SN.

## The governing structure — three groups, three jobs (do NOT collapse them)

For a 1-D infinite cylinder, phase space `(x,Ω)`, isometry acting as `(gx, R_gΩ)`:

| Group | Acts on | Reduces | Price on the fiber |
|---|---|---|---|
| `ℝ_z` translation | base only (`R_g=I`) | 3-D→2-D | **none** (flat connection) |
| `SO(2)_ẑ` rotation | base **and** fiber | 2-D→1-D | **the α redistribution** |
| `H=⟨σ_y⟩×⟨σ_z⟩=C_{2v}` isotropy at `(r,0,0)` | fiber only | `S²→S²/H` | **none** — the free lunch |

**THE CRITERION (generalises):** a geometry has an angular **redistribution term iff
its spatial-reduction group acts non-trivially on the direction fiber.** Translations
⟹ flat ⟹ slab/2-D-Cartesian have none; rotations ⟹ connection ⟹ sphere/cylinder do.
And **`n_levels > 1` ⟺ the isotropy group is FINITE** (continuous `H` collapses the
quotient fiber to a point ⟹ one level; finite `H` leaves a 2-D quotient ⟹ one
coordinate redistributed, the other a conserved LABEL). One criterion explains
slab-has-no-α, sphere-has-one-level, cylinder-needs-levels.

## Adjudications that are theorems, not opinions

- **A level is a FIBER OF AN INVARIANT, not an orbit.** For continuous `SO(2)` the
  fiber = the orbit (so "level = orbit" is true for `product`), but `level_symmetric`
  fibers over `|μ_z|` = **two** `SO(2)`-orbits (`rules_sphere.py:213`) while `product`
  fibers over **signed** `μ_z` (`rules_product.py:126-140`). **One type, two
  invariants** — Smell #16 shape 2, and the direct cause of the two different `τ_raw`
  fingerprints (`[0,1,0,1,…]` vs `[1,½,½,0]` = residual group order 2 vs 4).
- **Level-partition and half-range-fold are NOT two rungs of one chain.** The level is
  the fibration by the *invariant* of the spatial-reduction group; the fold is a
  quotient by the *isotropy* group. Fibered structure, not a tower. **Inside the fiber
  the tower IS genuine and the two folds COMMUTE** (`C_{2v}` abelian ⟹
  `(S²/⟨σ_y⟩)/⟨σ_z⟩ ≅ S²/C_{2v}`) ⟹ `level_symmetric`'s 4-to-1 needs **no new object**,
  same object with a bigger `G`.
- **`σ_x σ_y = C₂(ẑ)` links the SPATIAL and ANGULAR quotients.** The `r=0` centreline's
  physical deck transformation is `C₂(ẑ)` (`η→−η` AND `ξ→−ξ`); the code uses `σ_x`.
  They differ on 56/64 ordinates at `product(4,16)` and coincide **only on the
  quotient** by `⟨σ_y⟩`. ⟹ the tree is right by a coincidence whose precondition (the
  fold) it has not taken; the two quotients are ONE decision split across two packages.
- **Which `n_φ` admits which mirror** — `φ→−φ` (the FOLD's `σ_y`) closes for **every**
  `n_φ`; `φ→π−φ` (`σ_x`) and `φ→φ+π` (`C₂ẑ`) need **even** `n_φ`. So **ERR-042 is NOT
  the fold's existence condition** — it is the CENTRELINE's, and it gets the right
  answer only because `σ_y`-closure is automatic. (I got this backwards in draft 1;
  see the deliverable's SELF-CORRECTION.)
- **Orbifold weights are forced, not chosen.** Orbit-stabilizer `W = w·|G|/|Stab|` ⟹
  `w` on the mirror-plane (fixed) nodes, `2w` off it ⟹ **the trapezoid
  `[w,2w,…,2w,w]` IS the orbifold measure**, not a quadrature taste. Burnside
  `#orbits=(N+F)/2` reproduces the measured 20/8/12 and 36/8/28 exactly. The R12a
  `τ_raw` trichotomy is then a function of **whether the cover meets the singular
  locus**, not of the "rule family".
- **The ξ-even/ξ-odd split IS the trivial/sign isotypic decomposition of `Z₂`.**
  `dim V₊^{(ℓ)} = ℓ+1`, `dim V₋^{(ℓ)} = ℓ` (character `χ_ℓ(σ_y)=1`) reproduces the
  measured `ℓ=1: 2 even, 1 odd`. For `C_{2v}`,
  `dim V_{A₁}^{(ℓ)} = ((2ℓ+1)+1+1+(−1)^ℓ)/4` ⟹ `ℓ=1` gives **1**, not 2 ⟹ **a 4-fold
  and a 2-fold rule on the SAME geometry have DIFFERENT trial spaces** ⟹ the
  restriction must be DERIVED from the group the measure carries, never hardcoded.
  The measured `+2.94` is the **cross-Gram `⟨Y_0,Y_1^ξ⟩ = ΣW_kξ_k > 0`** — the fold
  makes `frame.gram` non-diagonal, so `GramStructure` becomes a lie (L-010 corollary).

## The API verdict — SPLIT the object; do NOT fuse it into `Frame`

```
NOT a Frame: measure.quotient(G) -> DiscreteMeasure     # measure category
             = pushforward(rep_map).consolidate()       # `consolidate` is THE gap
IS   a Frame: SymmetryAdaptedHarmonicBasis(G,L) -> Basis
             GalerkinFrame(that, folded_measure)        # hierarchy unchanged
```

Discriminators, decisive first: **(D2)** `conjugate(A)=R∘A∘M` has no equivariance
precondition, so on a non-equivariant `A` it silently returns the group average
`(1/|G|)Σ gAg⁻¹` — a wrong number with a valid type. **(D3)** `DiscreteMeasure.
pushforward` already exists and its docstring already documents the atom-merging case
as *"mathematically valid; numerically the node array will contain duplicates"* — the
missing verb is **one method, `consolidate()`**, not a class hierarchy. **(D4)** no
operator owns it: the group is in the commutant of **every** equivariant operator at
once (→ lessons L-009 sharpening). **(D1, narrowed)** the fold IS expressible as an
`IndicatorBasis` Galerkin frame (`W = M𝟙`) and homogenisation is a near-miss
counter-example — D1 only survives as *"a frame emits coefficients over an existing
space; the quotient must emit a new measure over a NEW space."*

**REFUTED with reason: the fold is NOT `PetrovGalerkinFrame`.** In this project PG
means a **solution weighting**; a symmetry fold carries none, and its Gram stays
diagonal at exactly the parent value. Calling it PG repeats the category error B3.0
fixed when it moved the Lambertian out of the geometry slot. It is **Galerkin on a
smaller space** — 2nd sighting of that verdict (1st = the DSA `ℓ=0` sub-block frame).

## What is STILL missing (the highest-value half — none of it was in the proposal)

1. **The orbifold SINGULAR SET `Σ = Fix(σ_y) = {ξ=0}` has SIX ad-hoc ε-detectors and
   no name** (`augmented_mesh.py:815`, `pole_angular_closure.py:1012`,
   `loss_representation:2902`, `rules_sphere.py:212`, `directional.py:477`,
   `symmetry.py:952`). Every #326 pathology lives on it. `Quadrature.octants` ALREADY
   computes it as the zero-sign-label partition entry and nobody consumes it as such.
   Smallest change, largest reach.
2. **The group ACTION on nodes is implemented THREE times** and discarded twice:
   `symmetry.py:904 _orbit_closure` (returns `bool`!), `directional.py:170
   _compute_sphere_reflection_partners`, `moc/geometry.py:222 _reflected_azi_index`.
   THE missing primitive; `LevelStructure` gaining fields is a consequence.
3. **`degree_of_exactness` has no SUBSPACE attached** — the folded rule is exact at
   the parent degree **on `V₊`** and meaningless on `V₋`; `restrict`/`pushforward`
   just drop it, losing a provable claim. That gap IS the `+2.94` failure mode.
4. **`SubgroupOfO3` has containment but no quotient/residual-group operation**, and
   `restrict`/`pushforward` set `invariance_group=None`, conflating *trivial* with
   *unknown* — which makes the (mathematically genuine) second fold **non-composable
   in code**.
5. **No operator-side equivariance predicate** ⟹ the descent cannot be gated; a ξ-odd
   source would be silently symmetrized. Owes a NEGATIVE test (ξ-odd source ⟹ RAISE).
6. **"half-range" is one word for two operations** — `restrict`'s own docstring claims
   it for the slab `μ>0` case. Discriminator is **total mass**: restrict DROPS mass,
   quotient PRESERVES it. Reserve half-range for the restriction; the new object is
   the **fold/quotient**.

## Cross-method: NOT SN-specific — MoC is a live second consumer

`moc/geometry.py:222-229` — MoC's azimuthal grid is on `[0,π)`, i.e. MoC **already
took the quotient `S¹/⟨−I⟩ = RP¹`** and never named it, and `φ→π−φ` is a second
quotient's deck transformation implemented with its own `argmin(abs(...))` matcher +
hand-rolled wrap. Its latent defect is exactly ERR-042's kind: on a grid where `π−φ`
is not a node it returns the NEAREST index silently. A design serving only SN would
be a design serving one caller.

**Pollination ACQUIRE:** symmetry-adapted spherical harmonics (Altmann; Bradley &
Cracknell) — the projection-operator basis construction `P_χ=(dim χ/|G|)Σχ(g)*g`.
Same literature lineage `symmetry.py:59-73` already cites (Hamermesh, Stiefel-Fässler),
and it is what `D_{6h}`/hex will need.
