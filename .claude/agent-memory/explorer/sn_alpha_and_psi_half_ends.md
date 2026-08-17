---
name: sn-alpha-and-psi-half-ends
description: SN curvilinear — alpha is one-sided (far end is an FP residual, never hard-set; the bare-assert check was PROMOTED to a raising guard on 2026-08-12 and now covers BOTH arms, but it is O_h-blind); the psi-half block DOES march both decoupled ends but only the -1 leg seeds the angular recurrence, and the two ends' answers disagree by ~8% on any non-flat fixture with zero gates comparing them.
metadata:
  type: reference
---

Durable structure behind the Q6.5 "α at both ends" question. Measured 2026-08-12
on `refactor/operator-strategy-layers` @ `adb73fd5`; full map (with `file:line`)
was `scratch/q65_alpha_ends_map.md`. The SHAPE survives churn — re-derive
addresses via Nexus.

## The α recursion is ONE-SIDED, and its far end is an accident

- Owner is the **geometry factory**, not the closure: `spherical_streaming` /
  `cylindrical_streaming` in `orpheus/geometry/reduced_operator.py`.
  `alpha = np.zeros(M+1)` then `alpha[m+1] = alpha[m] - w[m]*eta[m]` — the ONLY
  structural zero is index 0. There is **no** assignment to the last index
  anywhere in `orpheus/`, no second recursion seeded at the far end, no
  renormalisation.
- `α_{M+1/2}` is therefore a floating-point residual: `[M]` ≤ `2.8e-16` on every
  shipped rule (sphere GL N=2…128: `0` to `8.2e-17`, ≈ 1 ULP of the dome peak;
  `folded_product`: `0` exactly for M ≤ 3, `1.1e-16` typical). It never drifts.
- ⛔ **The one production check WAS disabled in the canonical run** —
  `spherical_streaming` carried `assert abs(alpha[N]) < 1e-12`, a bare `assert`
  that `python -O -m pytest` strips; the cylinder arm had nothing. **Verified by
  experiment (2026-08-12)**: under `-O`, a measure whose α closes at `−0.6`
  raised nothing.
  ✅ **FIXED `bea6a367` (2026-08-12): `_assert_alpha_dome_closes(alpha, coord=,
  level=)` in `orpheus/geometry/reduced_operator.py` RAISES, on BOTH arms,
  `_ALPHA_CLOSURE_ATOL`.** `spherical_streaming` calls it right after
  `alpha_dome(mu, w)`.
  ⚠ **But it is structurally blind to a wrong-DOMAIN rule.** `α_{M+1/2} = −Σ w μ`
  is the rule's FIRST MOMENT, and every `O_h`-invariant cubature has it exactly
  zero: `[M]` 2026-08-13, `Σ w·μ_x = 0.000e+00` bit-exact on `level_symmetric(4/6/8)`,
  `product(4,8)`, `folded_product(4,8)`, and `-5.6e-17` on `lebedev(9)` ⟹ the dome
  CLOSES on all five and the guard passes. It catches a truncated/mis-ordered
  measure, never a 3-D rule handed to a 1-D arm (that is P3's job — #336).
- Four test gates assert the far end, all tolerance-based (`atol` 1e-13/1e-14),
  and all **telescoping-invariant** (L-014) ⟹ none is evidence about the
  ordinate enumeration. The only cylinder gate parameterises `folded_product`
  (4,8) and (4,6) — i.e. M ∈ {3,4}, the bit-exact end of the range.
- Contrast worth remembering: the **cell partition** producer
  (`angular_cell_edges_per_level`) HARD-SETS both cylinder endpoints
  (`edge_omega[0]=π`, `edge_omega[M]=0.0`) and the sphere's near end
  (`mu_edge[0]=-1.0`). So the partition closes exactly where α does not, and the
  two are computed in different modules from different formulas.

## The angular march is one-sided everywhere

`_psi_half_grid_single_level` (matvec) and the walk's `psi_angle` thread (solve)
both run `for m in range(M)` / `enumerate(ordinates_in_level)` from a single
index-0 seed. The ONLY descending loop is `angular_adjoint`'s reverse-mode AD
retrace. No meet-in-the-middle logic exists. (Reed & Lathrop's two-marches-meeting
prescription is cited in the tree only for the SPATIAL weight and the τ lineage,
never for α.)

## ⭐ BOTH decoupled ends ARE marched — but only one is consumed

The ψ½ state block (`RadialCharacteristicField`, R13) carries **two legs per
carried level**, `−1` (μ=−1 / ω=π) and `+1` (μ=+1 / ω=0), each solved by the SAME
`carlson_inward_sweep_from_source` engine (the `+1` leg = the same call on
reversed data — orientation is DATA, never a flag), chained by the pole
continuation ψ½⁺(0)=ψ½⁻(0).

**Only `cells(p, -1)` ever reaches the angular recurrence** — `precompute_psi_state`,
`RadialCharacteristicSeeding.apply` (`A_AB`) and its transpose all read the `−1`
leg exclusively. In the forward value path the `+1` leg reaches the bulk through
exactly ONE channel: its r=R outflow corner, reflected into the `−1` inflow
corner on a **specular** outer face. On a **vacuum** face that channel is severed
and the `+1` leg influences nothing.

Symmetrically, the recurrence's terminal face `faces[:, M, :]` (= ψ_{M+1/2}) is
produced and **never read**: `cell_contribution` consumes
`upstream_per_ordinate = faces[:, :-1, :]`, the adjoint never seeds index M, and
`downstream_per_ordinate`/`downstream()` have zero production consumers.
Consistent with `c_out[M−1] = α_{M+1/2}/τ ≈ 0`.

## ⭐ `[M]` The two ends disagree by percent, and NOTHING measures it

`faces[:, -1, :]` and `cells(p, +1)` are the same physical quantity (cell-centred
ψ at μ=+1 / ω=0, same shape). Measured, `get_mixture("A","1g")`, R=2, 10 cells,
`inner_tol=1e-12`:

| fixture | relative gap |
|---|---|
| sphere GL8, **vacuum** | **8.06e-02** |
| sphere GL8, **reflective** (flat solution) | 1.83e-12 |
| cyl `folded_product(4,8)`, vacuum | **2.2–2.6e-02** per level |

The gap does NOT vanish under refinement: sphere `8.3e-2 → 5.1e-2` over 32×
h-refinement (ratio ≈1.0–1.2, not 2) and GROWS with N (3.98e-2 at N=4 →
8.76e-2 at N=32, plateauing); cylinder INCREASES with h-refinement,
`2.1e-2 → 3.2e-2`.

⟹ **Every flat/reflective curvilinear fixture in the tree is Mode-7 blind to
this** (the flat solution makes both routes agree to round-off). The vacuum
fixture is the discriminator. No test names both objects; the α gates are
telescoping-invariant; the seed gates check `faces[0]`, never `faces[-1]`.
This is an *unreconciled redundancy*, not a wrong shipped number — because
nothing reads the terminal face.

## The σ_y-folded cylinder arc (Q5)

`folded_product(n_mu, n_phi)`: `n_levels = n_mu`, **`M = n_phi/2` per level**,
all ordinates `ξ > 0`. Nodes are strictly-interior arc MIDPOINTS
(`ω/π = [0.875, 0.625, 0.375, 0.125]` at n_phi=8) — this is candidate **(B)**
(Hébert-midpoint) of the #326 half-range map, NOT (A) fold-existing-nodes.
Neither ω=π nor ω=0 carries an ordinate; both are hard-set partition edges.
The level's two endpoint directions are the two senses of ONE ray
(`ξ=0`, `η=∓sinθ_p`) — the exact analogue of the sphere's diameter, realised as
the per-level path rescale `dr/|η_start|`. ⟹ **every folded-cylinder level is
R12a-CARRYING** (measured `radial_characteristic_levels == (0,1,2,3)`), so the
whole two-leg ψ½ machinery — and the gap above — runs `n_mu` times per solve.
The sphere is the `M = N`, single-level case of the same structure.

Contrast: the REFUSED full-circle `Quadrature.product(4,8)` closes α **bit-exactly**
(the level covers the arc twice, η pairs identically) — one more reason α[−1]≈0
is not evidence about the enumeration.
