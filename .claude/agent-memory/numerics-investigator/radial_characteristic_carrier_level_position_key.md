---
name: radial-characteristic-carrier-level-position-key
description: The ψ½ (RadialCharacteristicSpace) carrier is keyed by LEVEL POSITION (p_idx), not the walk's mu_level_idx `level`; the direct-march cells_view(p_idx) is NOT a bug — gate-source ≡ space-key invariant. Gates the RayOp/A_BB carve.
metadata:
  type: project
---

# ψ½ carrier is keyed by LEVEL POSITION — `cells_view(p_idx)` is correct

**Verdict (2026-07-07, branch `refactor/sn-walk-unification`):** the direct ψ½
solve in `orpheus/sn/loss_representation/__init__.py` (~L4078-4117) passing
`p_idx` (the `enumerate` position) to `RadialCharacteristicSpace.cells_view(buf,
level, sign)` is **NOT a bug** — not even latent. The RayOp / A_BB
(`RadialCharacteristicOperator`) carve MUST key by **`p_idx` (level POSITION)**,
NOT the loop's `level`.

**Why:** `RadialCharacteristicSpace.levels` are **level POSITION indices**
(`enumerate(raw)` in `augmented_mesh.py:826-829`, `raw =
morel_montry_tau_raw_per_level`), the SAME coordinate that (a) keys the space
slots, (b) is the sweep's `p_idx`, and (c) indexes `level_ordinates_list[p_idx]`.
They are NOT sparse μ-ordinate indices. The gate `if p_idx in seed_levels`
(L4081) and the key validator `if level not in self.levels` (`_slot_key`,
`radial_characteristic_space.py:378`) both read the SAME
`mesh.radial_characteristic_levels` tuple (`seed_levels =
frozenset(radial_characteristic_levels)` L3696; `space =
for_levels(radial_characteristic_levels)` L852). ⟹ the gate can only admit a
`p_idx ∈ space.levels`, so `cells_view` NEVER crashes AND always hits the
processed level's own slot — for ANY carrying config (contiguous, sparse, or
multi-level).

**The loop's `level` is the WRONG key.** In the sweep, `level` is the walk's
`mu_level_idx` (L4162): `None` for the sphere (`levels=[None]` L3974), the int
level index for the cylinder. The proposed "fix" (pass `level` not `p_idx`)
would CRASH the sphere (`None not in (0,)`). `p_idx` is deliberate.

**Empirical (verbatim, `nx=4 ng=2`):** the concern's premise (multi-level /
sparse-value carrier) does not occur — R12a trichotomy:

| mesh | sweep n_levels | carrying | space.levels | τ_raw[0] per level |
|---|---|---|---|---|
| sphere GL4/GL8 | 1 | `(0,)` | `(0,)` | 0.399 / 0.392 ∈ (0,1) — CARRIES |
| cyl product 2×4 / 4×8 | 2 / 4 | `()` | `None` | 0.0 (dead, #229 rank-dup ψ₀) |
| cyl LS S4 / S8 | 2 / 4 | `()` | `None` | 1.0 (dead, thread wt (1−τ)=0) |

Sphere always exactly `(0,)` (one M-M level = whole quadrature); cylinder NEVER
carries; Cartesian no curvature. So p_idx==0==key on the only carrying geometry.

**Diagnostic (34 green under `-O`, real teeth):**
`derivations/diagnostics/diag_p_idx_vs_level_radial_characteristic.py` — the
gate-source ≡ space-key invariant + R12a trichotomy across a 16-quadrature
battery; a KEEPER (general property; reddens if a future quadrature carries >1
level or if the two sources ever diverge). Existing coverage: the carrier A1
gate (`tests/sn/mesh/test_radial_characteristic_carrier.py`) already pins
sphere-carries-1 / cyl-cart-absent; the p_idx/level coordination invariant was
the gap this diagnostic fills.

Sibling rulings: [[starting-direction-metric-gauge-derivation]] (the block's
G_sd metric), [[curvilinear-inverse-seed-taxonomy]] (is `(L+C).solve` an honest
inverse), [[coupled-block-operator-numerics]] (the A_BB block in the planned
CoupledOperator).
