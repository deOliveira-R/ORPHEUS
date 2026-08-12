# Q6.5 — the α ends and the ψ½ ends: what ORPHEUS actually realises

Explorer map, 2026-08-12. Branch `refactor/operator-strategy-layers`,
HEAD `adb73fd5`. **Line numbers are at that HEAD** — the *structure* below is
durable, the addresses are re-derivable via Nexus. `orpheus/`, `tests/`, `docs/`
were CLEAN at open and at close (`git diff --stat` empty on those paths), so no
intra-dispatch drift (L-007/L-012).

Probes (uncommitted, `scratch/`): `_q65_alpha_probe.py`, `_q65_ends_probe.py`,
`_q65_final_checks.py`, `_q65_gap_refine.py`.

---

## The one-paragraph answer

ORPHEUS realises **α_{1/2} = 0 as an axiom and α_{M+1/2} = 0 as an
ACCIDENT** — the recursion is strictly one-sided, forward from the seed, and the
far end is whatever `Σ wμ` happens to round to (measured `≤ 2.8e-16`
everywhere, never hard-set, never re-imposed). The **angular march is
one-sided too** (m ascending, in the matvec and in the solve; the only
descending loop is the reverse-mode AD retrace). But ORPHEUS **does** carry a
second decoupled-end object: the ψ½ state block has TWO legs, `−1` (μ = −1 /
ω = π) and `+1` (μ = +1 / ω = 0), each solved by the SAME `carlson_inward_sweep_from_source`
radial DD engine, chained by the pole continuation ψ½⁺(0) = ψ½⁻(0). **The two
ends are never reconciled**: the recurrence marched from the `−1` end produces a
terminal face ψ_{M+1/2} that is the same physical quantity as the directly
marched `+1` leg, and on a vacuum sphere they disagree by **8 %** (measured).
Nothing reads the recurrence's terminal face, so this is an unreconciled
redundancy today, not a wrong number.

---

## Q1 — Where are the α coefficients computed?

**Owner: `orpheus/geometry/reduced_operator.py` — the two streaming factories.
Not the closure.** (The closure owns τ; the geometry factory owns α + ΔA/w. The
split is the #236 Step-C ruling, restated at `reduced_operator.py:791-796`.)

### Sphere — `spherical_streaming`, `orpheus/geometry/reduced_operator.py:779-786`

```python
alpha = np.zeros(N + 1)
for n in range(N):
    alpha[n + 1] = alpha[n] - w[n] * mu[n]

# Verify GL antisymmetry: α_{N+1/2} ≈ 0
assert abs(alpha[N]) < 1e-12, (
    f"GL antisymmetry violated: α_{{N+1/2}} = {alpha[N]:.2e}"
)
```

`α_{n+1/2} = α_{n−1/2} − w_n μ_n`, `α_{1/2} = 0` by the `np.zeros` allocation
(index 0 is never written). Cited in the body comment as Hébert §3.9.4
Eqs. 3.423–3.424, ORPHEUS factor-of-2-absorbed normalisation.

### Cylinder — `cylindrical_streaming`, `orpheus/geometry/reduced_operator.py:885-893`

```python
alpha_per_level: list[np.ndarray] = []
for level_idx in angular_measure.level_indices:
    eta = angular_measure.mu_x[level_idx]
    w   = angular_measure.weights[level_idx]
    M   = len(level_idx)
    alpha = np.zeros(M + 1)
    for m in range(M):
        alpha[m + 1] = alpha[m] - w[m] * eta[m]
    alpha_per_level.append(alpha)
```

Identical arithmetic per level, with η = `mu_x` (radial cosine) in place of μ.
**No assert on this arm.** Provenance comment (`:877-884`): Hébert §3.9.3, credited
onward to Alcouffe & O'Dell 1986 and Lathrop & Carlson JCP 1:173 (1966); "cited
by Reed & Lathrop 1970 (their ref. 7)". Neither primary is local/read.

### Consumers

Read once at closure construction, `orpheus/sn/sweep/pole_angular_closure.py:1368`
(sphere) / `:1381` (cylinder), and immediately baked into the two coefficients
that are the only things production ever sees (`:1394-1408`):

```
c_out[m] = α_{m+1/2} / τ_m                                # the ψ_m diagonal
c_in[m]  = (1−τ_m)/τ_m · α_{m+1/2} + α_{m−1/2}            # the ψ_{m−1/2} numerator
```

A **second, independent** transcription of the same recursion lives in the
derivation/analysis module `orpheus/derivations/discrete/sn/angular_differencing.py:206
`alpha_dome`` (deliberately un-single-sourced — the structural-independence
constraint spelled out at `pole_angular_closure.py:614-619`).

---

## Q2 — Is α_{M+1/2} = 0 exact, enforced, or a floating-point residual?

**A floating-point residual, on BOTH geometries.** It is:

- **not hard-set** — the array is `np.zeros(M+1)` and every index `1..M` is
  overwritten by the forward recursion. Grep for an assignment to `alpha[-1]` /
  `alpha[M]` / `alpha[N]`: none exists anywhere in `orpheus/`.
- **not closed two-sidedly** — there is no second recursion seeded at the far
  end, and no renormalisation/telescoping correction.
- **not derived from the partition** — `angular_cell_edges_per_level`
  (`pole_angular_closure.py:837-990`) DOES hard-set both ends of the *cell
  partition* (sphere `mu_edge[0] = -1.0` at `:944`; cylinder
  `edge_omega[0] = np.pi`, `edge_omega[M] = 0.0` at `:980-982`), but α is not
  computed from that partition — it runs its own accumulation in the geometry
  factory. So the partition closes exactly at the far end while α does not.

### `[M]` measured last-α (probe `_q65_alpha_probe.py`, `_q65_final_checks.py`)

Sphere, `Quadrature.gauss_legendre(N)`:

| N | α[−1] | dome peak | rel. to peak | ε·peak |
|---|---|---|---|---|
| 2 | `+0.0` (exact) | 0.5774 | 0 | — |
| 4 | `−5.5511e-17` | 0.5213 | 1.07e-16 | 1.16e-16 |
| 8 | `+1.3878e-17` | 0.5058 | 2.74e-17 | 1.12e-16 |
| 16 | `+2.0817e-17` | 0.5015 | 4.15e-17 | 1.11e-16 |
| 32 | `+1.3878e-17` | 0.5004 | 2.77e-17 | 1.11e-16 |
| 64 | `+8.1532e-17` | 0.5001 | 1.63e-16 | 1.11e-16 |
| 128 | `−1.7185e-17` | 0.5000 | 3.44e-17 | 1.11e-16 |

Cylinder, `Quadrature.folded_product(n_mu, n_phi)` — worst |α[−1]| over levels:

| rule | N | n_levels | M/level | max\|α[−1]\| | peak |
|---|---|---|---|---|---|
| (4, 4) | 8 | 4 | 2 | `0.0` exact | 1.3624 |
| (4, 6) | 12 | 4 | 3 | `0.0` exact | 1.1124 |
| (4, 8) | 16 | 4 | 4 | `1.110e-16` | 1.2587 |
| (4, 16) | 32 | 4 | 8 | `1.110e-16` | 1.2345 |
| (4, 32) | 64 | 4 | 16 | `2.776e-16` | 1.2286 |
| (8, 16) | 64 | 8 | 8 | `1.110e-16` | 0.7177 |
| (8, 32) | 128 | 8 | 16 | `5.551e-17` | 0.7142 |

⟹ the residual is at or just above 1 ULP of the dome peak; it never drifts.
Small-M cases (`M ≤ 3`, and every `M = 2`) close bit-exactly, because the level's
η values pair antisymmetrically with a constant weight and the partial sums cancel
identically.

### Gates that assert it

Four, all tolerance-based, none exact:

| gate | assertion | note |
|---|---|---|
| `tests/sn/sweep/curvilinear/test_pole_angular_closure.py:143-159` `test_alpha_endpoints_zero_gauss_legendre` | `alpha[0] == 0.0` (exact) AND `abs(alpha[N]) < 1e-13`, N ∈ {4,8,16} | hand-rebuilt recursion, GL only |
| `tests/sn/sweep/curvilinear/test_sph_sweep_regression.py:48-60` `test_alpha_boundary_conditions` | `assert_allclose(alpha_half[-1], 0.0, atol=1e-14)`, N ∈ {4,8,16,32} | reads `reduced.alpha_half` |
| `tests/sn/primitives/test_quadrature.py:239-259` `test_alpha_boundary_zero` | `assert_allclose(alpha[-1], 0.0, atol=1e-13)` per level | **the only cylinder gate**; `folded_product` (4,8) and (4,6) only — i.e. `M = 4` and `M = 3`, both near the bit-exact end of the table above |
| `tests/sn/sweep/core/test_sweep_regression.py:227-231` | `assert_allclose(alpha_half[-1], 0.0, atol=1e-14)` | sphere |

⛔ **The one production-side check is disabled in the canonical test run.**
`reduced_operator.py:784` is a bare `assert` in `orpheus/`, and
`.claude/rules/vv-testing.md` makes **`python -O -m pytest`** canonical — `-O`
strips it. Measured directly: under `python -O`, feeding `spherical_streaming` a
measure whose α closes at `−0.6` raises **nothing** and returns the operator.
(Under plain `python` it raises.) There is no cylinder counterpart at any
optimisation level.

⚠ Also note the gates are **telescoping-invariant** (L-014): α[−1] is a
permutation-invariant partial-sum endpoint, so none of the four can see an
ordinate-ordering error, and none is evidence for the recursion's *enumeration*.

---

## Q3 — The angular march direction

**Strictly one-sided, m ascending, everywhere. No two-ended / meet-in-the-middle
logic exists in the tree.**

### The matvec (apply) path — the kernel

`orpheus/sn/sweep/pole_angular_closure.py:1176-1187`, `_psi_half_grid_single_level`
(the SINGLE source of the recurrence algebra; the mesh-bound
`_psi_half_grid_for_level` at `:1445` and the public
`compute_psi_half_per_level` at `:1233` both delegate here):

```python
ng, M, nx = psi_level.shape
psi_half = np.empty((ng, M + 1, nx), dtype=psi_level.dtype)
if psi_half_seed is None:
    psi_half[:, 0, :] = 0.0
else:
    psi_half[:, 0, :] = psi_half_seed
for m in range(M):
    tau_m = tau_level[m]
    psi_half[:, m + 1, :] = (
        psi_level[:, m, :] - (1.0 - tau_m) * psi_half[:, m, :]
    ) / tau_m
return psi_half
```

One `for m in range(M)`. The only seed slot is index 0.

### The solve (sweep/walk) path — the same direction, a different body

`orpheus/sn/loss_representation/__init__.py:4206-4314`. Per level:

```python
psi_angle = np.zeros((ng, nx))                    # the M-M thread buffer
for m_local, global_n in enumerate(ordinates_in_level):
    ...
    psi_a_in_chain = psi_angle[:, chain].copy()
    ...
    psi_angle_out_chain_p = (
        geom.tau_inv[global_n] * psi_avg_chain_p
        - geom.mm_a_in_coeff[global_n] * psi_a_in_chain
    )
    psi_angle[:, chain] = psi_angle_out_chain_p
```

`enumerate(ordinates_in_level)` — ascending m again, threading one buffer forward.
(Note `psi_angle = np.zeros(...)`: on a CARRYING level the walk starts the thread
at **exactly zero** because A_BB is solved up front and substituted; comment at
`:4193-4205`.) The per-cell slow path (`DiamondDifference.update`,
`orpheus/transport/spatial/diamond.py:226-230`) runs the same
`(ψ̄ − (1−τ)ψ_in)/τ` step with the same `angular_upstream → outgoing_angular_state`
threading.

### The only descending loop, and why it is not a second march

`pole_angular_closure.py:1673-1677` in `angular_adjoint`:

```python
psi_half_bar = np.zeros((ng, M + 1, nx))
psi_half_bar[:, :M, :] += upstream_bar
for m in range(M - 1, -1, -1):
    tau_m = tau_p[m]
    phb = psi_half_bar[:, m + 1, :]
    psi_bar[:, level_idx[m], :] += phb / tau_m
    psi_half_bar[:, m, :] += -((1.0 - tau_m) / tau_m) * phb
```

This is the reverse-mode AD retrace of the ascending forward (`Lᵀ`), not a
physical march. Note `psi_half_bar[:, :M, :]` — index `M` receives **no**
cotangent, the exact mirror of Q4's finding that `faces[:, M, :]` is never read.

⛔ **Rejected candidate.** "There is a two-ended march hidden in the ψ½ block"
— the ψ½ block genuinely has two legs, but they are marches in **RADIUS** at
fixed μ = ∓1, not marches in ANGLE; the α-cascade is untouched by them
(see Q4). Grep for `meet|two-ended|both ends|bidirectional|opposite end` over
`orpheus/sn/` + `docs/theory/methods/sn/` returns nothing on this topic.

---

## Q4 — The starting-direction ψ½ seed, and the OTHER decoupled end ⭐

**Yes — the other decoupled end IS computed, by the same engine. No — it is not
used to close the angular march, and the two ends' answers are never compared.**

### (a) Both legs exist and both are directly marched

`orpheus/sn/operators/radial_characteristic.py:517-551`,
`RadialCharacteristicOperator.solve` (= `A_BB⁻¹`, the exact direct inverse):

```python
dr_path   = dr / start_cosines[level]
q_minus   = comp.interior.cells(level, -1)
q_plus    = comp.interior.cells(level, +1)
corner_in = comp.boundary.corner(level, -1)
# inward starting-direction leg: r = R inflow corner → pole
cells_minus, pole_face = carlson_inward_sweep_from_source(
    q_minus, sigma, dr_path, corner_in,
)
# outward leg: SAME engine on reversed data, entering at the pole-continued
# face (ψ½⁺(0) = ψ½⁻(0)) and exiting at r = R.
cells_plus_rev, corner_out = carlson_inward_sweep_from_source(
    q_plus[:, ::-1], sigma[:, ::-1], dr_path[::-1], pole_face,
)
```

So the μ = +1 / ω = 0 end (where α_{M+1/2} = 0 also decouples the balance) has
an exact analogue of the μ = −1 march. Orientation is carried by **data
reversal**, never a flag — the engine is one body
(`orpheus/sn/sweep/psi_half_angle_seed.py:110-185`, and its docstring `:124-135`
spells the reversed-call recipe explicitly). The residual/transpose siblings
(`psi_half_angle_seed.py:305-334`, `:365-397`) carry both legs symmetrically.
The `+1` leg is first-class STATE by ruling R13
(`orpheus/numerics/spaces/radial_characteristic_space.py:53-61`: the flat layout
is `for sign in (-1, +1): cells (ng,nx), corner (ng,)`).

### (b) But only the `−1` leg reaches the angular recurrence

Exhaustive `.cells(` consumer scan over `orpheus/`:

| site | reads |
|---|---|
| `pole_angular_closure.py:1551` `precompute_psi_state` — **the M-M seed** | `cells(p, -1)` ONLY |
| `radial_characteristic.py:850-852` `RadialCharacteristicSeeding.apply` (`A_AB`, ray → bulk) | `comp.interior` via `precompute_psi_state` ⟹ `cells(p, -1)` ONLY |
| `radial_characteristic.py:931` `…Seeding.apply_transpose` | writes `cells(p, -1)` only; its own docstring `:887`: *"The `+1` leg and corners stay zero"* |
| `loss_representation/__init__.py:3150` (walk comment) | *"its `.cells(p, -1)` read — the M-M recurrence seed"* |
| `radial_characteristic.py:1131-1143` `…Reconstruction.apply` (moment fold → q½) | writes BOTH signs (it is the SOURCE side) |
| `boundary.py:1058-1063` `_reflect_corner` | `out.corner(level, -1) = seed.corner(level, +1)` — the specular closure |

⟹ In the **forward value path** the `+1` leg influences the bulk through exactly
one channel: its **r = R outflow corner**, reflected into the `−1` leg's inflow
corner on a specular outer face. On a **vacuum** outer face that channel is
severed and the `+1` leg is solved but reaches nothing.

### (c) The recurrence's own terminal face is produced and never read

`_MMHalfGrid` (`pole_angular_closure.py:502-601`) exposes
`upstream_per_ordinate = faces[:, :-1, :]` — its docstring `:569-571` states the
trailing face is *"not consumed as anyone's upstream"*. `cell_contribution`
(`:1593`) reads `upstream_per_ordinate`. The `downstream_per_ordinate` /
`downstream(m)` accessors have **zero production consumers** (grep: only
`tests/sn/sweep/curvilinear/test_compute_psi_half_per_level.py:177-212`, accessor
smoke tests). Structurally consistent: `c_out[M−1] = α_{M+1/2}/τ_{M−1}` is
`≈ 0` (measured `+2.28e-17` on GL8), so the last ordinate's balance has no
ψ_{M+1/2} dependence to begin with.

### (d) ⭐ `[M]` The two ends DISAGREE, and no gate compares them

`faces[:, -1, :]` (= ψ_{M+1/2}, the recurrence marched M steps from the `−1`
seed) and `cells(p, +1)` (= the directly marched μ = +1 leg) are cell-centred
ψ at the same direction, same cells, same shape `(ng, nx)`. Probe
`_q65_ends_probe.py` (1-group `get_mixture("A","1g")`, R = 2, 10 cells,
`inner_tol=1e-12`):

| fixture | \|faces[0] − cells(−1)\| (seed check) | \|faces[−1] − cells(+1)\| | relative |
|---|---|---|---|
| sphere GL8, **vacuum** | `0.0` | `3.539e-03` | **8.06e-02** |
| sphere GL8, **reflective** | `0.0` | `3.633e-14` | 1.83e-12 |
| cyl `folded_product(4,8)`, vacuum, level 0/3 | `0.0` | `4.110e-04` | **2.58e-02** |
| cyl `folded_product(4,8)`, vacuum, level 1/2 | `0.0` | `3.211e-04` | **2.19e-02** |

The **reflective** fixture's discrete solution is flat (`k_eff = 1.5`,
infinite-medium), so it agrees to round-off — every flat/reflective curvilinear
fixture in the tree is **Mode-7 blind** to this. The vacuum fixture is the
discriminator.

`[M]` The gap does **not** vanish under refinement (`_q65_gap_refine.py`):

```
SPHERE GL8, vacuum, R=2 — h-refinement          CYL folded(4,8), vacuum — h-refinement
  nx=  5   8.3332e-02                             nx=  5   2.1351e-02
  nx= 10   8.0617e-02  ratio 1.03                 nx= 10   2.5833e-02  ratio 0.83
  nx= 20   7.5555e-02  ratio 1.07                 nx= 20   2.8772e-02  ratio 0.90
  nx= 40   6.8151e-02  ratio 1.11                 nx= 40   3.0776e-02  ratio 0.93
  nx= 80   5.9820e-02  ratio 1.14                 nx= 80   3.1915e-02  ratio 0.96
  nx=160   5.1322e-02  ratio 1.17

SPHERE nx=20, vacuum — N-refinement:  N=4 → 3.98e-02, N=8 → 7.56e-02,
                                      N=16 → 8.45e-02, N=32 → 8.76e-02
```

Sphere: decays far slower than first order (ratio ≈ 1.0–1.2 per halving) and
**grows** with N, plateauing near 8.8e-2. Cylinder: **increases** with
refinement toward ≈ 3.2e-2. I am **not** claiming which side is right, nor that
they are obliged to agree at any order — only that they are two independent
discrete answers for one continuous quantity, that they differ by percent, and
that **nothing in the tree measures it** (no test names both objects; the α
gates are telescoping-invariant; the seed gates check `faces[0]`, never
`faces[−1]`).

### Rejected candidates (structural reason each failed)

1. **"α_{M+1/2} is hard-set to 0."** No assignment to the last index exists;
   only `np.zeros(M+1)` (which the loop overwrites from index 1) and the
   accumulation.
2. **"The recursion is closed two-sidedly (over-determined, à la Reed & Lathrop
   'either α_{m+1/2} or α_{m−1/2} is known')."** One `for m in range(M)`, one
   seed. R&L is cited in the tree only for the SPATIAL weight and the τ
   lineage (`reduced_operator.py:181, 207`; `pole_angular_closure.py:140, 158`;
   `curvilinear_one_group.rst:1948, 2051, 2058, 2114, 2117`) — never for a
   two-ended α closure.
3. **"The `+1` ψ½ leg is the far-end angular seed."** `precompute_psi_state`
   and `A_AB.apply` read `cells(p, -1)` exclusively; `A_AB.apply_transpose`'s
   own docstring says the `+1` leg stays zero.
4. **"`faces[:, M, :]` feeds the redistribution fold."** The fold consumes
   `c_in`/`c_out`, into which ψ_{M+1/2} has already been eliminated by
   substituting the recurrence; the matvec reads `upstream_per_ordinate`
   (`faces[:, :-1, :]`) and the adjoint never seeds index `M`.
5. **"The `angular_adjoint` descending loop is a second march."** It is the
   reverse-mode transpose of the ascending forward — same coefficients, reverse
   program order.

---

## Q5 — Cylinder vs sphere (the σ_y-folded arc)

`[M]` `Quadrature.folded_product(n_mu, n_phi)`, probe `_q65_alpha_probe.py`:

- **n_levels = `n_mu`** (one per polar level, i.e. per `μ_z = cos θ_p`).
- **M = `n_phi / 2` per level** — measured: `(·,4) → M=2`, `(·,6) → M=3`,
  `(·,8) → M=4`, `(·,12) → M=6`, `(·,16) → M=8`, `(·,32) → M=16`. Every level
  of a given rule has the same M. Every ordinate has `ξ = mu_y > 0` (the fold
  keeps one half-circle).
- **Node positions are strictly INTERIOR midpoints of the arc.** At
  `n_phi = 8`: `ω/π = [0.875, 0.625, 0.375, 0.125]`; at `n_phi = 16`:
  `[0.9375, …, 0.0625]`; at `n_phi = 6`: `[0.8333, 0.5, 0.1667]`. Stored
  η-ASCENDING ⟺ ω-DESCENDING, so the march runs ω: π → 0.
  **Neither ω = π nor ω = 0 carries an ordinate** — this is the Hébert-midpoint
  half-range family (candidate (B) of the #326 map), not the fold-existing-nodes
  family (A).
- **The two angular endpoints are ω = π (η = −sin θ_p) and ω = 0
  (η = +sin θ_p)**, and they are **hard-set** in the cell-partition producer,
  `pole_angular_closure.py:979-983`:

  ```python
  edge_omega = np.empty(M + 1)
  edge_omega[0] = np.pi
  edge_omega[1:M] = 0.5 * (omega[:-1] + omega[1:])
  edge_omega[M] = 0.0
  edges.append(sin_theta * np.cos(edge_omega))
  ```

  (Same producer refuses a full-circle level at `:967-978` — "a FULL-CIRCLE
  level covers each arc twice".)

**Does the fold change the endpoint decoupling? No — it makes the cylinder
MORE like the sphere, not less.**

- **Both** endpoints decouple exactly as the sphere's diameter does:
  α[0] = 0 by construction and α[M] = 0 to `≤ 2.8e-16` (table in Q2).
- The two endpoint directions of level *p* are the two senses of **ONE ray**:
  `ξ = 0`, `μ_z = cos θ_p`, `η = ∓ sin θ_p` — a ray through the axis in the
  (r, z) plane, exactly analogous to the sphere's diameter. Stated in
  `psi_half_angle_seed.py:118-124` ("the sphere's diameter (μ = −1) or a folded
  cylinder level's ξ = 0 ray (η = −sin θ_p, Q5.6)") and realised as the per-level
  path-length rescale `dr / |η_start|` (`radial_characteristic.py:144-169`
  `_march_start_cosines`; `mu_start_per_level = −sin θ_p` at
  `reduced_operator.py:913-917`).
- ⟹ **EVERY folded-cylinder level is R12a-CARRYING** — measured
  `sn.radial_characteristic_levels == (0, 1, 2, 3)` for `folded_product(4,8)`,
  versus `(0,)` (the single level) for the sphere. So the folded cylinder runs
  `n_mu` copies of the whole two-leg ψ½ machinery, and the Q4 two-end
  disagreement above is present on **every** one of them.
- The sphere is the M = N single-level case of the same structure
  (`pole_angular_closure.py:1367`: `self.level_indices = (np.arange(N),)`).

### One asymmetry worth flagging (sphere vs cylinder, and partition vs α)

| object | sphere | cylinder (folded) |
|---|---|---|
| cell partition, near end | `mu_edge[0] = -1.0` hard-set | `edge_omega[0] = π` hard-set |
| cell partition, far end | `mu_edge[N]` accumulated (`Σw`; in-tree measured residual `0` at N∈{2,4,8,32,64}, `+2.2e-16` at N∈{12,16}, `+4.4e-16` at N=48 — `tests/sn/sweep/test_angular_cell_partition.py:272-278`) | `edge_omega[M] = 0.0` **hard-set** ⟹ `η_{M+1/2} = +sin θ` exact |
| α, near end | 0 by allocation | 0 by allocation |
| α, far end | accumulated residual | accumulated residual |
| production far-end check | `assert` at `:784`, **stripped by `-O`** | **none** |

So the cylinder's *partition* closes exactly at both ends while its *α* does
not, and the two objects are computed in different modules from different
formulas (`α_{m+1/2} = −W_p ξ(ω_{m+1/2})` is the closed form implied by the
partition — see `.claude/agent-memory/explorer/cylindrical_sn_level_order_sensitivity.md`
— but production accumulates instead). Reporting the fact only; no fix proposed.

---

## Could-not-determine / out of scope

- Whether ψ_{M+1/2} and `cells(p,+1)` *ought* to agree at any particular order
  is a numerical-analysis question I did not attempt — the two are different
  discretizations (M-M weighted-diamond-in-angle vs DD-in-radius) of the same
  continuous ψ(r, +1). I report the measured gap and the absence of any gate,
  nothing more.
- Whether the ≈ 8 % gap is present on a **heterogeneous** or multi-group
  fixture: not measured (single-region, 1-group only).
- Whether `cells(p,+1)`'s own accuracy is affected by the fact that its source
  `q_plus` is the `(+1)^ℓ` Legendre fold while `q_minus` is `(−1)^ℓ`
  (`radial_characteristic.py:1131-1143`): not investigated — at ℓ=0
  (isotropic, this probe) the two folds are identical, so the probe cannot
  separate them. A P1+ fixture would.

## Bonus: the two non-folded cylinder families (contrast, `_q65_alpha_probe.py` §C/§D)

| rule | levels | M/level | α[−1] |
|---|---|---|---|
| `Quadrature.product(4, 8)` (full circle, refused by `angular_cell_edges_per_level`) | 4 | 8 | `0.0` **exact**, all levels |
| `Quadrature.level_symmetric(4)` (does carry `level_structure`) | 2 | 16, 8 | `1.110e-16`, `0.0` |

The full-circle rule closes α bit-exactly because its level covers the arc
twice and the η values pair identically — i.e. the "cleanest" α closure in the
tree belongs to the rule the partition producer now **refuses**
(`pole_angular_closure.py:967-978`). Another reason α[−1] ≈ 0 is not evidence
of anything about the enumeration (L-014).
