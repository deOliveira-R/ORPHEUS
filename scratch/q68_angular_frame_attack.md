# Q68 — is the 1-D η-march the right structure for the cylinder's angular variable?

**Frame attack, `cross-domain-attacker`, 2026-08-12.** Branch-verified reads
(lessons L-005): every `file:line` below was read directly from the worktree at
`/Users/rodrigo/git/nuclear/ORPHEUS`, **not** from Nexus. No shell in this
dispatch — every number I state is either cited from an in-tree `[M]` record
(with its source) or **derived symbolically here and marked `[D]`**. Nothing is
measured by me.

**Verdict, one paragraph.** The (η-march + weighted-diamond closure)
decomposition is **not** the native discretisation of this operator. The native
one is named, it is uniform across sphere and cylinder, and ORPHEUS already
ships every primitive it needs: **the curvilinear angular-redistribution
operator is the first-order factor of the Jacobi–Sturm–Liouville operator
belonging to the level's own generating measure, and the level's ordinates are
that measure's Gauss nodes.** In that basis the operator is **tridiagonal**
(cost O(M), same as the march), its truncation is **exact at the collocation
nodes** (the (M+1)-th orthogonal polynomial vanishes at the M Gauss nodes),
particle balance is a **theorem at every truncation order**, its transpose is a
**matrix transpose**, and there is **no seed, no τ, no α-dome, no `MarchStart`,
no carrying-level admission gate**. Separately and independently: the march's
own frame diagnosis is that the continuous angular operator on a level admits
**no boundary condition at all** (its flux function vanishes at both arc ends),
so the seed is a discretisation artifact and the discrete system is *one
condition short and two data long*.

---

## STRUCTURAL FEATURES

Enumerated, no narrative. `[D]` = derived here; `[M-cite]` = an in-tree measured
record with its source.

**S1 — the phase space per level is a rectangle.** At fixed polar cosine
`μ_z`, the level's domain is `[0,R] × arc(ω ∈ [0,π])`. Structurally identical to
the sphere's `(r,μ)` rectangle already named in
`.claude/agent-memory/cross-domain-attacker/psi_half_seed_angular_trace_frames.md`.

**S2 — `μ_z` is EXACTLY conserved; the redistribution is exactly 1-D.** The
cylindrical streaming operator is `Ω·∇ = η ∂_r − (ξ/r) ∂_ω` with
`η = sinθ cos ω`, `ξ = sinθ sin ω`. **There is no `∂/∂μ_z` term.** The
polar cosine is invariant along every straight line in a 1-D cylinder, because
`ẑ`-translation is an exact symmetry. Corroborated structurally by the #326
attack's own table (`quadrature_symmetry_quotient_frames.md`): the
spatial-reduction group is `SO(2)_ẑ`, which acts on the fiber by rotating `ω`
only — *"one coordinate redistributed, the other a conserved LABEL"*.

**S3 — the redistribution term is a CONSERVATIVE first-order operator whose
FLUX FUNCTION vanishes at both endpoints of the arc.** From
`docs/theory/methods/sn/curvilinear_one_group.rst §balance-curvilinear` Step 1,
the discretised object is `∂_ω(ξψ)` with `ξ(ω) = sinθ sin ω`, and
`ξ(0) = ξ(π) = 0`.

**S4 — `α` IS the flux function, sampled at the cell edges.** `[D]` With equal
azimuthal weights `w_k = W/M` on an arc of total weight `W`, the recursion
`α_{m+1/2} = α_{m−1/2} − w_m η_m` (`orpheus/geometry/reduced_operator.py:712`)
integrates to

```
α(ω) = (W/π) · sinθ · sin ω = (W/π) · ξ(ω)                     [D]
```

exactly reproducing the in-tree `[M-cite]` `.claude/plans/q64_tau_partition_memo.md`
§9bis.1 premise **P-B**: *"α ∝ ξ(e_arc) to 1e-16; κ is an AMPLITUDE, not a
displacement."* ⟹ the "α-dome" is the graph of `sin ω`, and
`α_{1/2} = α_{M+1/2} = 0` is `sin 0 = sin π = 0`. It is not a "dome closing"; it
is the flux function's two zeros.

**S5 — the characteristic invariant of the coupled `(r,ω)` flow is the IMPACT
PARAMETER.** `[D]` Along a straight lab-frame line, `r sin ω = p = const`. The
`(r,ω)` characteristics are the level sets of `p`. `ω = π` is a repelling fixed
point of `dω/ds = −(sinθ/r) sin ω`; `ω = 0` is attracting.

**S6 — the arc carries an exact `Z₂`, and the closure family is
ANTI-equivariant under it.** `[D]` Under `R: ω ↦ π−ω`: `ξ → ξ` (α even),
`∂_ω → −∂_ω`, so `R A R = −A` for `A = ∂_ω(ξ·)`. The in-tree measured
`τ(π−ω) = 1−τ(ω)` and `Π_m (1−τ_m)/τ_m = 1` (brief `[M]`) are the discrete
shadow of exactly this. The **march** is not `R`-anti-equivariant: the seed slot
sits at one end only.

**S7 — the fold is already taken; the realisable basis is cosine-only.**
`folded_product` = `GL(n_μ) × trapezoid(n_φ, STAGGERED)` quotiented by `σ_y`
(`orpheus/numerics/quadrature/directional.py:681-694`). ψ is even in ξ; the
level's function space on the arc is `span{cos mω}`.

**S8 — ⭐ THE LEVEL'S ORDINATES ARE CHEBYSHEV–GAUSS NODES.** `[D]`
`STAGGERED = ½` (`orpheus/numerics/quadrature/rules_circle.py:105`) gives
`φ_m = 2π(m+½)/n_φ`. The `σ_y` fold keeps the `ξ>0` half, `M = n_φ/2` nodes at

```
ω_k = (k + ½)·π/M ,  k = 0 … M−1   ⟹   x_k ≡ η_k/sinθ = cos((2k+1)π/(2M))
```

which is **exactly** the Chebyshev–Gauss (first-kind) node set of order `M` —
the roots of `T_M`. The folded weights are equal, i.e. exactly the
Chebyshev–Gauss weights `π/M`. Degree of exactness `2M−1` against the weight
`(1−x²)^{−1/2}`. Denominator, per `plan-authoring` §2: this covers
**1 of 1** cylinder rule families ORPHEUS admits — `assert_carrying_quadrature`
(`orpheus/sn/sweep/pole_angular_closure.py:749`) refuses full-circle
`NODE_ALIGNED` products and level-symmetric rules outright.

**S9 — the SPHERE is the same construction with a different classical family.**
`[D]` GL(N) nodes are the roots of `P_N`; the flux function is `(1−μ²)`. In both
geometries the flux function is `(1−x²)·w(x)` with `w` the level's Gaussian
generating weight (Legendre `w=1 ⟹ 1−x²`; Chebyshev-1 `w=(1−x²)^{−1/2} ⟹
(1−x²)^{1/2}`) — i.e. the **Jacobi–Sturm–Liouville** factor
`(1−x)^{α+1}(1+x)^{β+1}` at `α=β=0` and `α=β=−½`. `orpheus/numerics/generating_measure.py:85`
already carries exactly this parameterisation (*"jacobi(-1/2,-1/2) is
Chebyshev-1"*).

**S10 — one free parameter per angular cell, TWO conditions on it.** `τ_m` must
be (a) the barycentric coordinate that forces Morel–Montry's `β = 0`
diffusion-limit consistency, and (b) `½ + O(w)` for second-order pointwise
angular truncation (Reed & Lathrop Eqs. 15/16). `[M-cite]`
`scratch/q68_curvilinear_angular_differencing_survey.md` HEADLINE: Lathrop 2000
Eq. (30) — *"only with `μ_m = μ̄` (`δ = 0`) is the truncation order
`O(Δμ²)`"*, where `δ = 2τ − 1`. One parameter, two conditions ⟹ over-determined
by construction.

**S11 — Smell #16 Shape 1, at maximum confidence: the angular closure is the
SPATIAL scheme registry, hand-rolled.** Character-for-character:

| | expression |
|---|---|
| `orpheus/transport/spatial/scheme.py:1197` | `return (psi_bar - (1.0 - w) * face_in) / w` |
| `orpheus/sn/sweep/pole_angular_closure.py:1184` | `(psi_level[:, m, :] - (1.0 - tau_m) * psi_half[:, m, :]) / tau_m` |
| `scheme.py:1215-1216` (VJP) | `avg_cot = face_out_bar / w` ; `in_cot = (-(1.0 - w)/w) * face_out_bar` |
| `pole_angular_closure.py:1676-1677` | `psi_bar[…] += phb / tau_m` ; `psi_half_bar[:, m, :] += -((1.0 - tau_m)/tau_m) * phb` |

Same expression, same reduction tree, both the forward AND the transpose.
**And the same module family already imports the shared primitive for the OTHER
march**: `orpheus/sn/sweep/psi_half_angle_seed.py:105` binds
`DiscretizationSchemeBase.outgoing_face_from_average_transpose` for the radial
seed chain. Half the collapse has landed; the angular half has not.

**S12 — the discarded second datum.** `RadialCharacteristicOperator.solve`
(`orpheus/sn/operators/radial_characteristic.py:526-538`) marches BOTH legs:
the inward `ω=π` leg into `cells(p,−1)` and, on reversed data through the SAME
engine, the outward `ω=0` leg into `cells(p,+1)`. Only the first is consumed:
`radial_characteristic.py:814` — *"the `cells(p,-1)` leg is read (the recurrence
seed); the `+1` leg and corners are `A_BB`-internal."* By S5, the two legs are
the **two halves of ONE characteristic** (`p = 0`, the central chord); the
pole-face continuation in the code already says so.

**S13 — boundary handling as a special case grafted onto interior logic.** The
whole `MarchStart` / `march_start_structure_per_level` / `non_carrying_levels` /
`assert_carrying_quadrature` / `edge_extrapolated_seed` / `_carrying_levels`
apparatus (`pole_angular_closure.py:622-800`, `1456-1507`) exists to supply and
police ONE datum at the arc end.

**S14 — no bond dimension.** `[M]` (brief) the recurrence never mixes cells;
`psi_half[:, m+1, :]` reads `psi_level[:, m, :]` at the same spatial index.

**S15 — cost asymmetry.** The march is a **Python `for m in range(M)` loop**
(`pole_angular_closure.py:1182`) over `(ng, nx)` slices. Its adjoint is a second
Python loop (`:1673`). `M ∈ [4, 32]` for shipped rules.

---

## ELEGANCE DETECTOR HITS

Seven fired. Two or more warrants a probe; seven is a native-frame-not-found
signature.

- **Smell #3 (magic correction term with hand-waved derivation)** — the α-dome
  is documented as "how the angular flux redistributes due to curvature" and its
  endpoint zero as a "closure contract" with a `1e-12` tolerance guard
  (`reduced_operator.py:762-820`). Per S4 it is `sin ω`, and its endpoint zero is
  `sin 0 = 0`. A guard, a tolerance constant and a raise exist for an identity.
- **Smell #4 (symmetry present in the problem, absent in the method)** — S6. The
  exact `Z₂` anti-equivariance of the closure family is broken by the one-sided
  march.
- **Smell #5 (convergence rate suboptimal by a known margin)** — S10 + the
  brief's `[M]` "plain diamond beats the principled weighting by 2.5–3.0× at
  finer azimuthal orders". Lathrop Eq. (30) names the margin as an *order*, not
  a factor.
- **Smell #7 ("we picked N because it works")** — nothing selects `ω = π` over
  `ω = 0` as the march start. Both are available (S12).
- **Smell #9 (boundary handling as a special case added to interior logic)** —
  S13, and it is the largest single block of the module.
- **Smell #13 ("this converges because we checked")** — `A(M) = 2.41…9.44` is a
  measured amplification with no spectral statement about the composite
  operator, and the existing neutral-stability gate `Π(1−τ)/τ = 1` is
  `[M-cite]` **passed by all four τ variants including "shuffled"**
  (`scratch/q68_flux_dip_discriminator.md` §1) — the gate cannot see what
  separates them.
- **Smell #16 Shape 1 (structurally distinct paths to ONE operator)** — S11, at
  maximum confidence: identical expressions, identical reduction trees, one of
  the two already consuming the shared source in a sibling file.

**Candidate NEW smell (2nd sighting → promote):** *"a discretisation
manufactures a boundary condition the continuous operator does not admit."*
First sighting was the sphere's ψ½ (recorded in
`psi_half_seed_angular_trace_frames.md` as *metric-invisible-yet-active DOF* —
"the angular advective flux `(1−μ²)` VANISHES at μ=−1"). This is the second, on
the cylinder, and it supplies the DEEPER cause that the first sighting only
saw the symptom of: **the weight is zero because the flux function is zero, and
a vanishing flux function means the endpoint is characteristic — the continuous
problem has no boundary condition there at all.** TELL: an elaborate
admission/classification apparatus (`MarchStart`, carrying/non-carrying,
`assert_*`) around a datum whose continuous counterpart carries zero measure AND
zero flux.

---

## FRAME CANDIDATES

Ranked. Rank = elegance payoff × cheapness of the discriminating test.

---

### ⭐ FRAME 1 (rank 1) — Jacobi–Sturm–Liouville ladder on the level's own Gauss nodes

*(A.3 harmonic analysis + A.2 Clebsch-Gordan/ladder + A.3 approximation theory,
firing together. This is the native frame.)*

**Trigger.** S8 + S9: the level's ordinates ARE the Gauss nodes of a classical
orthogonal family, and the redistribution's flux function IS that family's
Sturm–Liouville factor `(1−x²)w(x)`. Also S3 (conservative first-order operator)
and S7 (the realisable basis is exactly that family).

**Reformulation.** Let `x = η/sinθ_p ∈ [−1,1]` and `A = ∂_ω(ξ ·)`. Expand ψ on
the level in the family's own polynomials. `[D]`, both derivations one line of
trigonometry / one Legendre recurrence:

*Cylinder level* — basis `{T_k(x)} = {cos kω}`, `k = 0…M−1`:

```
A cos(kω) = ½(k+1)·cos((k+1)ω)  −  ½(k−1)·cos((k−1)ω)          [D]
```

*Sphere* — basis `{P_ℓ(μ)}`, `ℓ = 0…N−1`:

```
∂_μ[(1−μ²)P_ℓ] = ( ℓ(ℓ−1)·P_{ℓ−1} − (ℓ+1)(ℓ+2)·P_{ℓ+1} ) / (2ℓ+1)   [D]
```

Both are **tridiagonal with zero diagonal** — a pure raising/lowering ladder,
structurally the same object as the `ℓ ↔ ℓ±1` Clebsch–Gordan coupling that the
PN streaming recurrence already is (lessons L-009).

Four properties drop out, each a theorem rather than a check:

1. **Exact particle balance at EVERY truncation.** `[D]` The `j=0` row is
   identically zero: the down-step coefficient vanishes at `k=1`
   (`−½(k−1)|_{k=1} = 0`; `ℓ(ℓ−1)|_{ℓ=1} = 0`), and no other mode reaches `j=0`.
   So `Σ_k w_k (Aψ)_k = 0` for every ψ and every truncation order. This is
   Lathrop 2000's stated design requirement (*"telescopes the angular-derivative
   term away, giving exact particle balance — a property of the analytic
   transport equation, which is very useful to preserve"*) obtained
   structurally instead of by an α-recursion plus a runtime closure guard.
2. **The truncation is FREE.** `[D]` The top mode's up-step produces
   `T_M` (resp. `P_N`) — which **vanishes at all `M` (resp. `N`) collocation
   nodes, because they are its own roots.** Dropping it introduces exactly zero
   error at the nodes. There is no ad-hoc top-mode closure to invent.
3. **No seed, no τ, no α, no `MarchStart`.** The operator is
   `A_ord = C⁻¹ L C` with `C` the family's Vandermonde (DCT-II on the cylinder,
   the Legendre Vandermonde on the sphere — the latter is already in the tree as
   the angular-frame machinery, `harmonic_frame_ownership_funk_hecke`). It
   depends only on the quadrature, so it is precomputed **in the slot
   `_c_in_per_level` / `_c_out_per_level` occupy today**
   (`pole_angular_closure.py:1394-1409`).
4. **The transpose is `A_ord.T`.** The 20-line hand-rolled reverse recurrence at
   `pole_angular_closure.py:1661-1690` is deleted.

**Elegance payoff (4 of 4 criteria).**

- *Structure-exposing.* The redistribution stops being "a correction term with a
  recursion and a dome" and becomes **the first-order factor of the Jacobi
  operator of the level's own measure**. One sentence covers sphere and cylinder
  and predicts which classical family each uses from its quadrature alone. It
  also explains, for free, why the fold's cosine-only subspace is *invariant*
  under redistribution (the ladder is closed on it).
- *Expressive.* Sphere and cylinder become one body with a different
  `generating_measure` — the module's existing per-geometry `if coord is
  SPHERICAL / CYLINDRICAL` branches in `angular_cell_edges_per_level`
  (`:940-990`) and in `MorelMontryAngularSweep.__init__` (`:1361-1392`) collapse.
- *Structurally-simpler.* Deleted, not moved: `angular_cell_edges_per_level`
  (~150 lines incl. two refuted-convention essays), `morel_montry_tau_per_level`,
  `_assert_tau_within_unit_interval`, `alpha_dome`, `_assert_alpha_dome_closes`
  + `_ALPHA_CLOSURE_ATOL`, `MarchStart`, `march_start_structure_per_level`,
  `non_carrying_levels`, `assert_carrying_quadrature`, `edge_extrapolated_seed`,
  `_edge_seed_stencil`, `_psi_half_grid_single_level`, and the
  `angular_adjoint` reverse recurrence.
- *Algorithmic-advantage.* One precomputed `einsum` (or DCT pair) replaces a
  sequential Python loop over `M`; the transient amplification `A(M) = 9.44`
  has no counterpart; the truncation error at fixed `M` is spectral for
  ω-smooth ψ instead of `O(δΔω + Δω²)`.

**Cost — stated honestly, three items.**

1. **The sweep loses angular triangularity.** Today `c_out[m]` sits in the
   denominator and `c_in[m]·ψ_{m−1/2}` in the upstream numerator, so ordinate
   `m`'s balance needs only edges below `m` — a forward substitution. A dense
   `A_ord` makes each cell's angular problem an `M×M` solve. The matrix depends
   on `(cell, group)` through `Σ_t V` and the geometry, so it is one
   precomputable factorisation per `(cell, group)` and an `O(M²)`
   back-substitution per sweep visit. For 1-D radial meshes (`nx ≤ 320`,
   `M ≤ 32`) that is cheap. **It would not transplant to 2-D/3-D cylinders**
   without re-derivation — say so before anyone generalises it.
2. **Gibbs.** A spectral angular basis rings on a ψ that is non-smooth in ω —
   an annulus's tangent-ray kink, a directional boundary source. The march does
   not (though `A(M)=9.44` and the negative-ψ̂ observation say it is not
   monotone either). Mitigation is structural: land it as another
   `PoleAngularClosureBase` registry member (`key="jacobi_ladder"`), selectable,
   so the comparison against M-M on the same problem IS the experiment.
3. **⛔ The published failure mode applies verbatim.** Larsen & Morel 2010
   §1.5.1 `[M-cite]` via `scratch/q68_curvilinear_angular_differencing_survey.md`:
   *"the starting-direction flux is computed … with greater accuracy than the
   other directions; hence, significant accuracy is actually lost if the
   starting-direction flux plays no role in the angular derivative treatment."*
   That is exactly why linear-discontinuous-in-angle LOST to weighted diamond. A
   seedless spectral scheme discards both accurate legs. **This is the frame's
   sharpest risk and it has a named fix** — see Frame 2's synthesis.

**First test (discriminating; no solver, one script).** Extend the existing
predicate ladder in `orpheus/derivations/discrete/sn/angular_differencing.py`
(which already carries P0–P4) with **P5: exactness on `cos 2ω`**.

- Probe the production angular operator to a dense `M×M` matrix (unit-ordinate
  probes through `precompute_psi_state` + `cell_contribution`).
- Build `A_ord = C⁻¹ L C` from the ladder above.
- Apply both to `ψ_k = cos(kω)` sampled at the level's nodes, `k = 0,1,2,3`, and
  compare against the analytic `∂_ω(sinθ sin ω · cos kω)` at the same nodes.

**Prediction.** `k=0` and `k=1` (constant, and affine in η): **both exact** —
M-M by BMC-43's defining property, the ladder by construction. `k=2`: the ladder
is exact to machine precision (`array_equal`-class, ≤ few ULP); M-M carries a
residual `O(Δω²)‖ψ‖`.
**What would falsify it.** If M-M is *also* exact on `cos 2ω`, its
reconstruction is richer than affine-in-η and the frame's claimed advantage
evaporates. That is the implementation this test would PASS which I am claiming
is inferior — so the test can fail, and it costs one script and no solver.

**Structural attack on current.** The shipped nodes are the optimal collocation
points of a spectrally-accurate polynomial-in-η method, and the code interpolates
between them **piecewise-linearly** and then differences the result. The
information-theoretic content of the level's `M` ordinates is `2M−1` degrees of
exactness; the M-M closure extracts one (affine).

**Precedent.** None in `scripts/`. Nearest in-tree kin:
`harmonic_frame_ownership_funk_hecke.md` (an operator OWNS the frame that is its
eigenbasis) — this is the same rule applied to *streaming* rather than
*scattering*, and the answer is the ladder basis, not an eigenbasis, because the
operator is skew rather than self-adjoint.

---

### ⭐ FRAME 2 (rank 2) — degenerate (Fichera) advection: the operator admits NO boundary condition; the scheme is one condition short and two data long

*(A.3 functional analysis / trace theory; A.5 dynamical systems.)*

**Trigger.** S3 — the flux function `ξ = sinθ sin ω` vanishes at both ends of
the arc. S12 — two endpoint data exist and one is discarded.

**Reformulation.** The weak form on a level, for any test function φ:

```
∫₀^π ∂_ω(ξψ)·φ dω  =  [ξψφ]₀^π  −  ∫₀^π ξψ·φ' dω  =  −∫₀^π ξ ψ φ' dω     [D]
```

**The boundary term vanishes identically, for every ψ and every φ.** The Fichera
function `b·n = −ξ n_ω` is zero at both ends: neither endpoint is inflow, neither
is outflow, both are characteristic. This is a degenerate first-order operator
(Fichera / Oleĭnik–Radkevič class), **not** an initial-value problem.

⟹ Three consequences, in order of bite:

1. **`ψ_{1/2}` is not physical boundary data.** It is a datum the *marching*
   discretisation demands and the continuous operator does not supply. The
   entire S13 apparatus is a discretisation artifact wearing a physics name.
2. **What the endpoints DO supply is a COMPATIBILITY CONDITION, and there are
   two of them.** Setting `ξ = 0` in the balance leaves a plain radial ODE at
   `η = −sinθ` (`ω=π`) and at `η = +sinθ` (`ω=0`). Those are exactly the two
   legs `carlson_inward_sweep_from_source` already marches (S12). **They are two
   ends of ONE characteristic** (S5, `p = 0`), which is why the code's own
   `pole_face` continues one into the other.
3. **The discrete count.** `[D]` Unknowns per level: `M` cells + `M+1` edges
   = `2M+1`. Equations: `M` balances + `M` closures = `2M`. **Under-determined by
   exactly one** — while the continuous problem is well-posed with none. So the
   scheme is *one condition short*, and *two data long*. Today it consumes one
   datum, discards the other, and nothing checks that they agree.

**⟹ the answer to "two-point BVP vs IVP" is NEITHER.** It is a
no-boundary-condition problem with two compatibility conditions and a
one-condition-deficient discretisation. Posing it as a BVP is closer to right
than posing it as an IVP, but the honest object is a **constrained
(tau-method) Galerkin problem** — which is precisely how Frame 1's risk-3
resolves:

> **The synthesis.** Expand ψ on the level in `{T_k}` (Frame 1). Impose the two
> independently-marched leg values as the two endpoint constraints
> `Σ_k (−1)^k ĉ_k = ψ(ω=π)` and `Σ_k ĉ_k = ψ(ω=0)`, dropping the two
> highest-wavenumber balance equations. **This is Lanczos' tau method.** It
> uses BOTH accurate legs (answering Larsen & Morel's objection head-on), it is
> exactly `Z₂`-symmetric because the two constraints are each other's `R`-images
> (Frame 3), and it needs no seed *slot* — the legs enter as constraints, not as
> a state block.
> ⚠ **NAME HAZARD.** The codebase already overloads `tau` three ways
> (`pole_angular_closure.py:999` — angular closure weight, optical depth,
> critical half-thickness). A "tau method" would be a fourth. Spell it
> `endpoint_constrained` / `lanczos_constraint`, never `tau`.

**Elegance payoff (3 of 4).** *Structure-exposing*: names the seed as an
artifact and the endpoint data as compatibility conditions, which is what
dissolves the S13 apparatus. *Structurally-simpler*: `MarchStart`'s two-fact
trichotomy and `assert_carrying_quadrature` have nothing left to police.
*Algorithmic-advantage*: **available today at zero cost** — the discarded leg is
a free a-posteriori error estimator for the level's angular discretisation.
(*Expressive*: no separate gain beyond Frame 1's.)

**First test (cheapest in this document — both data already exist).**
On a converged `folded_product` cylinder solve, form per level, per cell, per
group

```
R_p = ψ_{M+1/2}^{marched}  −  radial_characteristic.cells(p, +1)
```

— the recurrence's terminal `ω=0` edge minus the independently-marched `ω=0`
radial-ODE leg. Same shape `(ng, nx)`, both live in the tree today.

**Prediction.** `‖R_p‖` converges at **order 1** in `Δω` (matching Lathrop
Eq. 30 for the shipped `δ ≠ 0`), and it tracks the anisotropic-MMS angular error
within a constant factor. Second half, the sharper one: monkeypatch `τ ≡ ½` —
the brief's `[M]` says the MMS error improves 2.5–3.0×; **`‖R_p‖` must improve
comparably.**
**What would falsify it.** `‖R_p‖` converging at order 2, or moving *against*
the MMS error under the τ≡½ mutation. Either outcome refutes the far-endpoint
compatibility defect as the controlling mechanism, and that is worth knowing
before anyone builds Frame 1.

**Structural attack on current.** An admission gate
(`assert_carrying_quadrature`) refuses quadratures on the grounds that they
cannot supply a datum the continuous problem does not have — while the second
copy of that datum, which the same operator computes, is documented as
`A_BB`-internal and thrown away.

---

### FRAME 3 (rank 3) — `Z₂` anti-equivariance: the discriminator the existing gate is missing

*(A.2 group theory / representation theory.)*

**Trigger.** S6, and the `[M-cite]` fact from
`scratch/q68_flux_dip_discriminator.md` §1 that **all four τ variants (shipped,
diamond, reversed, shuffled) pass the tree's `Π(1−τ)/τ = 1` neutral-stability
gate** — *"whatever separates them here is invisible to that gate."*

**Reformulation.** The exact per-level redistribution operator satisfies
`R A R = −A` (S6). A discretisation should inherit it. The closure family does;
the seed-injected march does not, because the injection `e₀ ⊗ seed` breaks the
reflection. The `Z₂`-anti-equivariant scheme is the projection
`½(A_h − R A_h R)` — which by lessons L-009's third outcome is **Galerkin on a
symmetry sub-block, NOT Petrov-Galerkin** (a group identification, not a
solution weighting — third sighting of that verdict after the DSA `ℓ=0`
sub-block and the #326 fold).

**Consequence with a physical name.** `R` acts on `cos(mω)` with sign `(−1)^m`,
so the `Z₂`-odd sector is `span{cos(mω): m odd}` — which contains `m=1`, i.e.
`cos ω ∝ η`. **A scheme that is not `Z₂`-anti-equivariant can leak error into
the odd sector, and the odd sector's leading member is the radial current
`J_r`.** An anisotropic MMS is exactly the fixture that sees this.

**Elegance payoff (2 of 4).** *Structure-exposing*: names an exact discrete
symmetry the scheme silently breaks, and names which physical moment the
breakage contaminates. *Algorithmic-advantage*: supplies a gate that
discriminates where the shipped gate is blind. (Not structurally-simpler; the
symmetrised march costs 2×. Not more expressive on its own.)

**First test.** Probe the production per-level angular operator to a dense
`M×M` matrix `A_h`; form the reflection permutation `R` (`m ↦ M−1−m`, exact on
the folded equispaced arc); measure `‖R A_h R + A_h‖ / ‖A_h‖`. Run it on the
same four τ variants the flux-dip harness already builds.

**Prediction.** Closure-only part (zero seed, coefficient path): residual
`~1e-16`. Full seed-injected operator: `O(1)`, concentrated in rows `0` and
`M−1`. Shuffled and reversed τ **redden this probe while passing the existing
`Π(1−τ)/τ = 1` gate**.
**What would falsify it — and the honest caveat.** ⚠ **At `M = 2` the probe
cannot separate reversed from shipped**: the condition is `τ_m + τ_{M+1−m} = 1`,
which for `M=2` is `τ₁+τ₂=1`, satisfied by both orderings (shipped
`[0.42265, 0.57735]`, reversed `[0.57735, 0.42265]` — same sum). The test is
only meaningful at `M ≥ 4`. If all four variants still pass at `M ≥ 4`, the
probe adds nothing over the existing gate and the frame contributes nothing
operational.

**Structural attack on current.** The tree's only structural gate on the angular
scheme tests a **scalar product** (`Π(1−τ)/τ = 1`) where the underlying fact is
an **operator identity** (`R A R = −A`). The scalar is the `(0,M)` entry of the
operator statement; three of the four measured τ variants exploit exactly the
gap between them.

---

### FRAME 4 (rank 4) — degrees of freedom per angular cell (approximation theory / FEM-in-angle)

*(A.3 approximation theory; A.6 asymptotic analysis for the two-scale reading.)*

**Trigger.** S10 — one parameter, two conditions.

**Reformulation.** The angular discretisation is a **P0 finite-volume scheme**
(one unknown per angular cell) with a hand-designed reconstruction. Its
reconstruction carries one parameter; its accuracy conditions number two. The
frame's rule: **raise the local polynomial degree; do not tune the
reconstruction.** The two conditions then live at different orders of the
diffusive scaling parameter and both can be met.

The published instances are already catalogued in
`scratch/q68_curvilinear_angular_differencing_survey.md` §Q1.0, all sphere:

| scheme | DOF/cell | order | seed behaviour |
|---|---|---|---|
| weighted diamond (shipped) | 1 | `O(δΔμ + Δμ²)` | consumes the seed |
| linear-discontinuous in angle (Walters–Morel 1991) | 2 | `O(Δμ²)` pointwise, **4th** on absorption/leakage | **decouples from the seed** ⟹ loses |
| quadratic-continuous (Lathrop 2000) | 2 | 4th measured | consumes it |
| **hybrid**: quadratic-continuous in the FIRST cell, LD in the rest | 2 | *"superior accuracy relative to the weighted-diamond scheme"* (Larsen–Morel) | consumes it |

⭐ **The structural rhyme.** ORPHEUS's **spatial** axis already has exactly this
registry — `DiscretizationSchemeBase` with `diamond_difference` and
`linear_discontinuous` members, an `is_affine_scannable` trait, and a shared
`affine_scan_coefficients` contract. The **angular** axis has one hard-coded
scheme. *The angular axis is the spatial axis's registry, un-built.*

**Elegance payoff (2 of 4).** *Structure-exposing*: renames "the closure has a
floor" as "the cell is one DOF short", which is a countable claim.
*Expressive*: the angular scheme becomes a registry key like the spatial one.
(Not simpler — it adds unknowns. Algorithmic-advantage is conditional: 2× the
arithmetic per angular cell buys ≥2 orders, so it wins only past a crossover
that must be measured, per Lathrop's own refusal to cost it.)

**First test — and this is the one that decides issue #235.** On the anisotropic
curvilinear MMS at `nx = 320` (the instrument
`.claude/plans/q64_tau_partition_memo.md` §4.5 records as WORKING, range to
`9e-6`), sweep `n_φ = 8/16/32/64` for (a) shipped M-M τ and (b) `τ ≡ ½`, and
**fit the ORDER, not the ratio.**

- **Hypothesis A** (Lathrop Eq. 30 — `δ ≠ 0` demotes the angular truncation to
  first order): (a) fits order ≈ 1; (b) fits order ≈ 2.
- **Hypothesis B** (#235 as filed — a genuine 2-D `(η,φ)` closure gap): **both**
  flatline with no order.

The predictions are mutually exclusive. ⚠ The brief's `[M]` "2.5–3.0× at finer
azimuthal orders" is **consistent with A but is not evidence for it** — a ratio
is not an order, and this is the measurement the campaign has not made.

**Structural attack on current.** The floor is being attributed to a missing
angular *dimension* when the counting says it is a missing angular *degree of
freedom*, and one order-fit separates the two.

---

### FRAME 5 (rank 5) — Smell #16 Shape 1: single-source the angular closure onto the spatial scheme primitives

*(Part C Smell #16, Shape 1 — the promoted kernel.)*

**Trigger.** S11 — identical expressions, identical reduction trees, forward AND
transpose, with the sibling file already importing the shared source.

**Reformulation.** Replace the two hand-rolled bodies:

```python
# pole_angular_closure.py:1184  →
psi_half[:, m + 1, :] = DiscretizationSchemeBase.outgoing_face_from_average(
    psi_level[:, m, :], psi_half[:, m, :], w=tau_level[m],
)

# pole_angular_closure.py:1676-1677  →
avg_cot, in_cot = DiscretizationSchemeBase.outgoing_face_from_average_transpose(
    phb, tau_m,
)
psi_bar[:, level_idx[m], :] += avg_cot
psi_half_bar[:, m, :] += in_cot
```

**What becomes automatic or unspellable.** (a) The angular transpose stops being
hand-maintained — it becomes the same `w`-generic VJP the spatial reverse walks
ride (#311). (b) `LinearDiscontinuous`-in-angle becomes *reachable*: Frame 4's
Lathrop scheme 4 is already a registered member of the spatial family
(`w = 1/(1+k)` per the docstring at `scheme.py:1193`) and has simply never been
instantiated on the angular axis. (c) A `step`-in-angle and a future
`quadratic_continuous`-in-angle become registry members rather than new code
paths. (d) The angular march is revealed as a **cumulative-product affine scan**,
so `CumprodScan` applies and the Python `for m in range(M)` loop (S15) can be
vectorised.

**Bounded scope, stated up front.** The two axes differ in one real way: the
spatial closure solves for `ψ̄` from a source, whereas the angular closure takes
`ψ_m` as an *input*. So the shared primitive is the **closure pair**
(`cell_average` / `outgoing_face_from_average` + their VJPs), **not** the full
`affine_scan_coefficients` triple. Do not over-reach the collapse.

**Elegance payoff (3 of 4).** *Structure-exposing*: the angular march and the
spatial march are one closure family on two axes — which is exactly the framing
Frame 4 needs to be buildable. *Structurally-simpler*: two bodies become zero.
*Algorithmic-advantage*: conditional, via `CumprodScan`.

**First test.** `np.array_equal` — **0 ULP, not `allclose`** (lessons L-002) — on
the full half-angle grid and the full cotangent for every shipped sphere and
cylinder fixture.
**Prediction.** PASSES bit-identically; the expressions are character-identical
with the same reduction tree.
**Why it can still fail, and why that is the point.** If it REDs, the two
"identical" formulas were compiled to different reduction trees (e.g. one hoists
`(1−w)/w`), which means the existing `angular_adjoint` is **not** the exact
Euclidean transpose of the forward at ULP level — and the #208 dagger claim over
the angular block is weaker than the docstring at `:1599` advertises. Either
outcome is information; a green is a free single-sourcing and a red is a
correctness finding.

**Structural attack on current.** One module family imports
`outgoing_face_from_average_transpose` for the radial seed chain
(`psi_half_angle_seed.py:105`, with a comment explaining that it is the single
source) and re-implements the same two lines for the angular chain 1500 lines
away in a sibling file.

---

### FRAME 6 (rank 6) — differential geometry: α is a Christoffel symbol; the chart where it vanishes is MoC

*(A.1 differential geometry / exterior calculus. Diagnostic + one instrument.)*

**Trigger.** S5 — the `(r,ω)` flow has the conserved quantity `p = r sin ω`;
S4 — α is the flux function of a term that appears in one chart and not another.

**Reformulation.** The redistribution is the connection coefficient of the
moving `(e_r, e_φ)` frame along a straight lab-frame trajectory. Change chart
`(r,ω) → (p,s)` with `p` the impact parameter and `s` the signed chord
coordinate (`r = √(p²+s²)`, `sin ω = p/r`). In `(p,s)` the streaming operator is
`∂/∂s` alone: **the connection is flat, `α ≡ 0` identically**, no angular
coupling, no dome, no τ, no seed. `p` is constant along trajectories, so the
`(p,s)` chart is the normal coordinates of the trajectory manifold — **and it is
MoC.**

Two theorems drop out for free:
- **The two "independent" endpoints are one characteristic.** `ω = π` and
  `ω = 0` on a level are the two ends of the `p = 0` central chord (S5). The
  code's `pole_face` continuation already encodes this; the frame names it.
- **The `(η,φ)` axis pair in #235 is the wrong pair.** By S2 there is no
  angle-angle coupling to discretise: `μ_z` is exactly conserved and the levels
  are exactly decoupled. The coupling that genuinely exists is **`(r, ω)`** —
  space-angle, along `p = r sin ω`.

**Elegance payoff (2 of 4).** *Structure-exposing*: α is chart-dependent
bookkeeping, and the whole S13 apparatus is chart-dependent with it.
*Algorithmic-advantage*: only via the instrument below — full adoption is
"replace SN with MoC", which is not on the table (SN's fixed-ordinate structure
is what makes scattering a per-cell moment sum).

**First test — this frame's real deliverable is a MISSING INSTRUMENT.** The
memo `.claude/plans/q64_tau_partition_memo.md` §0a lists **three of five
instruments as dead** and the trajectory-resolvent cross-check as
**reference-limited at `≈3e-2`** (refining the SN side moves it the wrong way —
lessons L-049). The `(p,s)` chart supplies a replacement with no reference limit:

> **A purely-absorbing cylinder with the closed-form chord solution.**
> Homogeneous, radius `R`, `Σ_s = 0`, uniform isotropic `Q`, vacuum. For ordinate
> `(θ_p, ω_m)` at radius `r`, the exact angular flux is `[D]`
> ```
> ψ(r, θ_p, ω_m) = (Q / 4πΣ_t)·( 1 − exp(−Σ_t·L/ sinθ_p) ) ,
> L = r·cos ω_m + √(R² − r² sin²ω_m)
> ```
> Elementary: no quadrature, no iteration, no scattering, no reference solver.
> Sanity checks `[D]`: `ω=π ⟹ L = R−r` ✓; `ω=0 ⟹ L = R+r`, the full central
> chord ✓.

Refine `nx` until the residual plateaus (spatially converged — this is Lathrop
2000's own construction, where he integrated the angular ODEs to 6 digits to
kill the spatial term), then read the residual as **pure angular-closure error,
per ordinate, per level.**
**Prediction.** The residual is largest at the ARC ENDS, where `τ → ¼, ¾` puts
the node furthest from its own cell midpoint (Reed & Lathrop Eqs. 15/16), and it
converges at **order 1** in `Δω`.
**What would falsify it.** A residual that is flat across the arc, or order 2.
Either refutes the R&L-16 / Lathrop-30 diagnosis on the cylinder as transplanted.

**Structural attack on current.** The campaign's own instrument audit says its
best cross-check is reference-limited, and the reference-free instrument is one
elementary formula in the chart the geometry is already written in.

---

## CROSS-METHOD POLLINATION

**Current method class:** SN — curvilinear 1-D, the angular-differencing
sub-problem specifically.

**B1 — from MoC (in-tree).** *Trigger:* S5, the shared characteristic invariant.
Two concrete borrowings. (i) The `(p,s)` chart as the **α-free oracle** — Frame
6's purely-absorbing chord solution, which is exactly the object MoC integrates.
(ii) `moc/geometry.py:222` `_reflected_azi_index` is **the same reflection
permutation `R`** that Frame 3's anti-equivariance probe needs, and MoC's
azimuthal grid is already on `[0,π)`, i.e. MoC has already taken the same
quotient (recorded in `quadrature_symmetry_quotient_frames.md`). Shared
primitive, two consumers, currently three independent implementations of the
group action.
*First test:* the same `array_equal` witness test lessons L-013 prescribes —
assert `_reflected_azi_index`, `_compute_sphere_reflection_partners` and
`_orbit_closure`'s internal matching agree at 0 ULP on the folded arc.

**B2 — from the SPATIAL scheme registry (in-tree, same package).** Frame 5. This
is the strongest pollination in the document because both endpoints already
exist and one already imports the other.

**B3 — from PN / spherical-harmonics algebra.** *Trigger:* S9 — the
redistribution is a `ℓ ↔ ℓ±1` ladder, structurally identical to the
Clebsch–Gordan coupling that makes streaming the `ℓ=1` tensor operator in the PN
recurrence (lessons L-009). The cross-method statement: **ORPHEUS already treats
the scattering kernel spectrally (Funk–Hecke, the owned Galerkin angular frame)
and the streaming redistribution by finite differences.** Frame 1 is "give
streaming the same treatment scattering already has."
*First test:* Frame 1's P5 exactness probe.

**B4 — ⭐ from Fokker–Planck angular differencing (ACQUIRE).** *Trigger:* the FP
operator `∂_μ[(1−μ²)∂_μ ·]` has **the same degenerate endpoints and the same
Sturm–Liouville factor** as the redistribution operator, one differential order
up. Morel's FP angular difference scheme is built to be exact on the Legendre
eigenbasis and to preserve the leading angular moments — i.e. it is the
published precedent for Frame 1's construction, on the sibling operator.
`scratch/q68_curvilinear_angular_differencing_survey.md` §Q1.3 explicitly flags
*"the one adjacent case (Morel's Fokker–Planck angular scheme)"* and defers it to
a §Q1.5 **that has not been written** — the survey file ends at §Q1.4. **This is
the single highest-value open literature lead for this question.**
*Search terms for the `literature-researcher` pull:* Morel 1985 "An improved
Fokker-Planck angular differencing scheme" NSE; plus
`Galerkin quadrature` + `moment-preserving angular differencing` +
`Sturm-Liouville angular discretization transport`. My own single web probe
found only *spatial* Chebyshev-collocation RTE work — the angular-side
construction is **not placed**, and I am not claiming it is novel.

**B5 — from diffusion / DSA — ⛔ REFUTED, and the refutation is worth keeping.**
*Tempting trigger:* Frame 4's two-conditions-one-parameter conflict looks like
the DSA pattern (keep a high-order operator AND a consistent low-order one, so
neither has to carry both duties). **It does not fire: DSA changes the
ITERATION, not the FIXED POINT.** The `β ≠ 0` diffusion-limit defect is a
property of the converged transport solution, so it cannot be delegated to an
accelerator. The genuine resolution is Frame 4's (add a DOF to the angular
cell, not an operator to the iteration). Consistent with lessons L-007's
corollary: the mixed/saddle-point apparatus attaches to the diffusion member,
never to the sweep.

---

## UNEXPLORED / REJECTED — each with its one-line STRUCTURAL reason

First-class output. These are the frames the next session should **not**
re-attack.

- **Miller's algorithm / Gautschi minimal-solution (march the stable
  direction).** ⛔ **REFUTED by the measured `Π(1−τ)/τ = 1`.** `[D]` With
  `G(m) = Π_{k≤m}|(1−τ_k)/τ_k|`, the forward march amplifies a seed error by
  `G(m)` and the backward march by `G(M)/(G(M)/G(m)) = G(m)` — **the same
  profile**, because `G(M) = 1` exactly. A recurrence with unit end-to-end gain
  is *reversible*: there is no dominant/minimal solution pair, so Miller's
  algorithm has no purchase and reversing the march direction buys **nothing**.
  This is a high-prior, very tempting frame (the amplification `A(M)=9.44` reads
  like a classic unstable recurrence); the exact anti-symmetry of τ is what
  kills it. See SELF-CORRECTION.
- **Two-point BVP as such.** Not refuted but **superseded**: the problem has no
  boundary conditions (Frame 2), so "BVP" is the wrong noun. The right object is
  a constrained Galerkin problem with two *compatibility* conditions.
- **Tensor networks / MPO.** ⛔ Bond dimension **1** — S14, `[M]` the recurrence
  never mixes cells. A rank-1 chain is degenerate, not a network (lessons
  L-001). Do not promote until an `N ≥ 3` truncation knob actually ships.
- **de Rham cohomology / DEC / FEEC.** ⛔ The ω-grid is a **1-D** complex, so
  `H¹ = 0` trivially and the only content DEC could supply is telescoping —
  which the α-recursion already has by construction (S4). No `∂² = 0` structure
  beyond the trivial one.
- **Symplectic geometry / structure-preserving integrators.** ⛔ The conserved
  quantity IS `p = r sin ω` (S5), and the frame in which a discretisation
  conserves it by construction IS the characteristic chart — i.e. Frame 6, i.e.
  MoC. No **independent** lever; the symplectic reading and the connection
  reading have the same single consequence.
- **Category theory / functorial lifting.** ⛔ The compositional win it gestures
  at — one angular scheme family across sphere and cylinder — is captured
  concretely by Frame 1 (the generating measure indexes the family) and Frame 5
  (the existing `RegistryMixin`). No functor or law produces a test here.
- **Information theory / MaxEnt closures.** ⛔ MaxEnt fires on a *moment*
  closure — an angular shape unknown from a few moments. Here the shape is fully
  resolved by `M` ordinates and the closure is a *reconstruction between
  resolved nodes*. There is no entropy deficit to maximise.
- **Homogenization / multiscale in ω.** ⛔ No scale separation on a level: every
  ω-cell has the same width by construction (S8, equal weights). The genuine
  scale separation in this problem is thick-vs-thin **optical** (Frame 4's
  two-order reading), not angular.
- **Hierarchical matrices / low-rank compression.** ⛔ `M ∈ [4,32]`, and the
  operator is *tridiagonal* in the right basis (Frame 1). Nothing to compress.
- **Krylov / operator-based acceleration.** ⛔ No trigger *for the closure
  question*: the angular march lives inside the matvec, not inside a solve. (A
  Krylov question exists for the outer iteration; it is a different dispatch.)
- **Wiener–Hopf factorization.** ⛔ Wrong solver family (lessons L-001) — native
  to the Chandrasekhar/H-function half-space line, structurally incompatible
  with a sweep formulation, and keeping the two families independent is itself a
  V&V requirement.
- **Bloch / crystallographic.** ⛔ Not independent: the ω-grid's periodic
  structure is exactly what makes the DCT available, i.e. it *is* Frame 1.
- **Step-characteristic / exponential differencing in angle.** ⛔ Not a
  no-trigger — a **measured dead end in the primary literature**.
  `scratch/q68_curvilinear_angular_differencing_survey.md` §Q1.3: step-in-angle
  exists (Lathrop Eqs. 7–8), is `O(Δμ)`, decouples from the seed, and Lathrop's
  hybrid step+diamond *"displayed first-order error behavior in calculating
  absorption"*. Every "characteristic/exponential" hit in the corpus is
  **spatial**. The brief's question *"does the α-weighting already secretly
  encode a step characteristic?"* answers **no**: the exact solution of the
  angular advection along the arc is a rescaling by `ξ`, and the weighted
  diamond is a *polynomial* (barycentric) closure, not an exponential one — and
  by S3 there is no exponential attenuation in ω to encode, because the operator
  is purely conservative (no zeroth-order sink in the ω direction).
- **Adjoint / control theory.** ⛔ The angular adjoint already exists
  (`angular_adjoint`), and Frame 1 would make it a matrix transpose. No separate
  lever for the accuracy question.

---

## SELF-CORRECTION

Two, both caught before emission.

1. **Nearly shipped the Miller/Gautschi "march the stable direction" claim as a
   win.** The reasoning path was: `A(M)` grows to 9.44, the forward march
   amplifies early and damps late, therefore the backward march is
   unconditionally non-amplifying and the code already computes the backward
   seed. That is **wrong**: `G(M) = 1` exactly makes the recurrence reversible,
   so both directions carry the identical amplification profile `G(m)`. Caught by
   working the suffix product out rather than reasoning from the monotone τ
   profile. Recorded in REJECTED with its structural reason, because it is a
   high-prior trap that the measured `A(M) = 9.44` actively baits.
2. **Nearly asserted that the `(η,φ)` framing of #235 is simply wrong.** It is
   wrong *as an angle-angle stencil* (S2 — `μ_z` is exactly conserved), but two
   legitimate two-variable readings survive and are recorded in Frames 4 and 6:
   the CHART question (τ is barycentric in η while the cells and nodes are
   uniform in ω — a genuine η-vs-ω tension), and the space-angle `(r,ω)`
   coupling along `p = r sin ω`. Asserting the flat negative would have thrown
   both away. Trap caught: converting a structural theorem about one axis pair
   into a dismissal of the issue that raised it.

Neither is a register violation; the reclassification held (no balance, no
hedging on trigger matches, no closing pleasantry). Both are content errors that
the "a first test must be able to fail" discipline surfaced — in case 1 by
forcing me to write down what the test would measure.

---

## ADDENDUM — a stale memory this attack refutes

`.claude/agent-memory/cross-domain-attacker/psi_half_seed_angular_trace_frames.md`
carries the ruling *"Seed is a GENUINE independent DOF ONLY on SPHERE … On
CYLINDER the level's inflow edge η=−sinθ COINCIDES with the first azimuthal
ordinate ⟹ the seed IS ψ₀ ⟹ a DEAD DOF if carried."* **That is superseded.**
Q5.6 landed the `σ_y` fold; `folded_product` levels are strictly-interior arcs,
`MarchStart.on_edge_node` is False on every one of them, and
`assert_carrying_quadrature` now *requires* every cylinder level to carry a
first-class ψ½ seed. The memory file has been updated.

---

## ONE-PARAGRAPH RECOMMENDATION

Run **Frame 2's first test** first — it is free (both arrays exist today) and it
tells you whether the far-endpoint compatibility defect is the controlling error
term. Run **Frame 4's order-fit** second — it costs four MMS solves and it
*decides issue #235*, separating "the weighted diamond is first-order in angle"
from "we need a 2-D closure", which no measurement in the campaign currently
separates. Land **Frame 5** whenever convenient — it is a bit-identity
single-sourcing with a correctness finding as its only possible downside. Build
**Frame 1** only after Frame 2 and Frame 4 have reported, because if the order-fit
says Hypothesis A then Frame 4's cheaper hybrid (a second DOF in the first
angular cell) may capture most of the win, and if it says Hypothesis B then
Frame 1 is the only frame in this document that can deliver it. Frame 6's
purely-absorbing chord oracle is worth building regardless — the campaign's
instrument audit says its best angular instrument is reference-limited, and this
one is a closed-form expression.
