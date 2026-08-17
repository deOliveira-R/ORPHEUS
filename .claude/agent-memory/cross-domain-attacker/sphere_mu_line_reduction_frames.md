---
name: sphere-mu-line-reduction-frames
description: Issue #336 — the 1-D spherical μ-line reduction of a 3-D cubature IS the continuous-isotropy quotient (pushforward is only its mechanism); GAUSS OPTIMALITY refutes the whole REDUCE branch; the admissibility guard is CDF-INTERNALITY, independent of every moment condition; the Stage-0 gate is landed and unwired.
metadata:
  type: project
---

Frame attack on GitHub #336's fork: REFUSE a non-μ-line rule at the spherical
`SNMesh`, or REDUCE it (derive the induced μ-line rule). User's framing offered
for attack: "REDUCE is a pushforward of the angular measure under Ω ↦ Ω·r̂, and
`DiscreteMeasure.quotient(G)` and it are two instances of ONE thing."

Grounded on direct worktree reads (L-005); **NO SHELL this dispatch** — every
number below is `[D]` (derived by hand from the tree's own formulae) or
`[M-cite]`. **Why:** the fork decides whether machinery gets built at all.
**How to apply:** reach for this on any dimensional-reduction-of-a-quadrature,
angular-collapse, or "should we auto-convert the user's rule" question.

## The decisive finding — GAUSS OPTIMALITY refutes REDUCE as a value path

The quotient space is **1-D with a Lebesgue reference** (Archimedes: the
pushforward of uniform-`S²` under `Ω↦Ω·r̂` is `2π·Leb` on `[-1,1]`, no
Jacobian). On such a space the maximal-degree theorem is absolute: `n` nodes ⟹
degree ≤ `2n−1`, attained **uniquely** by Gauss–Legendre. So every possible
reduction of every possible cubature to `n` distinct μ values is **≤ GL(n)**,
with equality **iff the reduced rule IS GL(n)**. ⟹ REDUCE's entire image is
(a) exactly `gauss_legendre(n)` — precisely when the parent is a `product` rule
whose μ-factor is GL, i.e. the reduce is the tensor projection — or (b) a rule
strictly dominated by a rule the user could have typed instead. There is no
third case, and it is a theorem, not a measurement.

`[D]` worked: reduce(`level_symmetric` S4) = 4 nodes `μ = ±0.3500212,
±0.8688903`, `w = (⅓, ⅔, ⅔, ⅓)` (mass 2 after ÷2π), **degree 5**, τ =
`(0.3933, 0.4750, 0.5250, 0.6067)` — fully ADMISSIBLE, all τ strictly interior,
dome closes. Against `gauss_legendre(4)`: same 4 nodes' cost, degree **7**. So
24 ordinates buy a rule 2 degrees worse than the 4-node rule already shipped.
(The degree is 5, not the parent's advertised 3, because LS's own defining
moment condition `15u²−10u+1 = 0`, `u = μ₁² = (5−√10)/15`, IS the μ-marginal's
degree-5 condition — the inheritance is `≥`, and here it is STRICT.)

⟹ **verdict: REFUSE on the value path; keep the reduction as (i) the DIAGNOSTIC
that writes the refusal message ("your rule induces a 4-node degree-5 μ-line
rule; use `gauss_legendre(4)`, degree 7") and (ii) a VERIFICATION oracle.** Bound
the claim honestly: dominated in the DEGREE metric the tree's own
`ExactnessClaim` records; a rule tuned for a non-polynomial class is outside it.

## The frame verdicts

- **Pushforward is the MECHANISM, not the SPECIFICATION.** `pushforward(φ)` takes
  an arbitrary callable, so it cannot say *which* φ; the sharper frame is
  **quotient by the CONTINUOUS ISOTROPY of the spatial-reduction group**, which
  derives φ (μ = Ω·r̂ is the complete invariant of `SO(2)_r̂`; `S²/SO(2) ≅
  [-1,1]` with the two poles as the orbifold's fixed points). This is **already
  in the tree**: `registry.AngularSymmetry(continuous_isotropy, discrete_residual)`
  with `support` DERIVED (`SO2 ⟹ SPACE_INTERVAL_M11`) and `admits_domain(measure)
  = (measure.support == self.support)`. One field should derive TWO things: the
  target space (landed) and the **morphism onto it** (missing). That kills the
  free-callable hazard by construction.
- **Archimedes' hat-box theorem is the license for "weights unchanged".** It is
  the ONLY reason the pushforward needs no Jacobian here, and it is the
  discriminator for a wrong φ: pushing forward under the polar ANGLE θ gives
  `sinθ dθ`, not `dθ`. μ is the equal-area coordinate; θ is not.
- **Disintegration is the ADJOINT half, and it names the round trip.** `R∘P = I`
  on the μ-line, `P∘R = Π₀` the azimuthal-average projector — a Galerkin frame
  (L-009), whose basis-side face is the `SO(3)↓SO(2)` `m=0` branching, i.e. the
  Legendre/zonal basis, i.e. **`angular_frame(0)`'s big brother** — already
  landed, do not re-mint.
- **REFUTED — the user's information-asymmetry argument.** "A quotient exploits a
  symmetry (no loss); a marginalization discards azimuthal dependence (loss)" is
  false as a discriminator: both are EXACT on the invariant subspace and lossy
  off it, and `DiscreteMeasure.quotient`'s own docstring already says so of the
  ξ-odd integrand. Same loss, same honesty.
- **The REAL discriminator: the fiber identification's KIND.** `quotient`'s is an
  INTEGER permutation certified by a group action (`orbit_certificate`);
  a continuous isotropy admits no finite node permutation (L-014 theorem A), so
  the reduce's fibers must be recovered from FLOAT coincidence
  (`consolidate`'s `round(x/atol)` bucketing — non-transitive at bucket
  boundaries, the ERR-042 hazard class). ⟹ the reduce must consume the integer
  partition the rule already stores (`level_indices`) instead — which
  immediately inherits the KNOWN #326 defect that `level_symmetric` fibers over
  `|μ_z|` (two `SO(2)` orbits) while `product` fibers over signed `μ_z`.
- **Unification: TRUE, and already landed.** Both are
  `pushforward(φ).consolidate()` — `quotient` documents exactly that composite.
  Do NOT mint a `Reduction` type or a `mode=` on `quotient`; add a third
  φ-builder. What legitimately DIFFERS per φ is which derived claim survives:
  `quotient` must drop `exactness` (its reference `φ_*λ` is an unnamed type), the
  reduce must KEEP it (the reference is `LEGENDRE`, a shipped type, degree ≥
  parent's).

## Admissibility: the guard the measure frame cannot see

The sphere's angular cell edges are CUMULATIVE WEIGHTS from `μ_{1/2} = −1`
(`angular_cell_edges_per_level`), so `τ_m = (μ_m − e_m)/w_m = ½ + δ_m/w_m` with
`δ_m = μ_m −` cell midpoint. Hence:

- `α`-dome closure `α_{M+1/2} = −Σw_mμ_m = 0` **is FREE** — it is degree-1
  exactness, inherited by the pushforward from any parent with degree ≥ 1, **at
  every axis r̂**. A design that demands parent σ_r-symmetry for it is over-strict.
- `τ ∈ [0,1]` ⟺ `|μ_m − mid_m| ≤ w_m/2` ⟺ **the rule's CDF staircase straddles
  the Lebesgue CDF at every node** (call it INTERNALITY / interlacing; the
  cumulative-weight partition is the monotone transport map from `dμ` to the
  rule, and internality says it moves no node out of its own cell). This is
  **implied by NO moment condition**. `[D]` witness: the Chebyshev-node
  equal-weight rule at N=8 (`μ = ±0.98079, ±0.83147, ±0.55557, ±0.19509`,
  `w = ¼`) closes the dome exactly and gives `τ₁ = −0.3259`.
- A node ON the isotropy's FIXED POINT (`μ = ±1`) is the reduce's version of
  `quotient`'s singular-set condition — same condition, and `[D]` `lebedev(26)`
  violates it: reduced μ = `{−1, −1/√2, −1/√3, 0, +1/√3, +1/√2, +1}`, so
  `μ₀ = −1` ⟹ `on_edge_node` ⟹ `τ₀ = 0` ⟹ non-carrying march.

## The gate that is BUILT and UNWIRED (L-013's sibling sub-shape)

`AngularSymmetry.admits_domain` answers #336 correctly TODAY and is TESTED
(`tests/numerics/test_registry.py:683 assert not slab.admits_domain(lebedev)`) —
its only production caller is `select_quadrature` (`registry.py:915`), the
**advisory/selection** path. The `SNMesh` **construction** path never asks. Same
for `assert_carrying_quadrature`: `march_start_structure_per_level` already has a
`SPHERICAL` arm (`mu[0] == −1.0` / `mu[0] == mu[1]`), and the guard's message is
cylinder-only prose. ⟹ REFUSE is ~one call site + a coord-parameterized message +
the internality conjunct. The proposal's scale ("build the reduction machinery")
is an order of magnitude above the gap.

## Cross-method

SAME operation, parameterised by the isotropy `H` (one table, `AngularSymmetry`
already its home; `registry.py:656`'s `NotImplementedError` is literally the
seam): `H = SO(2)_r̂` slab + 1-D sphere · `H = C_{2v}` the cylinder fold
(finite, landed) · `H = SO(3)` the P1/diffusion + `iso` + CP angular average
(quotient = a point; basis side already `angular_frame(0)`) · `H = ⟨−I⟩`
even-parity/SAAF on `RP²` and MoC's azimuthal `[0,π) = RP¹` (finite — but the
lattice carries no `C_i` member today, so this one needs a new
`_NamedSubgroup`).

⛔ **MoC's polar quadrature (Tabuchi–Yamamoto, #355) merely RHYMES** — the 2-D
geometry's reduction group is TRANSLATION `ℝ_z`, which acts trivially on the
direction fiber, so **no isotropy is spent and the fiber is not reduced**; the
polar split is a Fubini TENSOR factorization of a rule that stays on `S²`, and
its reference is `sinθ dθ` against an exponential integrand — so the "no
Jacobian" theorem fails. Its true relative is the ANALYTIC pushforward:
Bickley–Naylor `Ki_n(x) = ∫₀^{π/2} e^{−x/sinθ} sin^{n−1}θ dθ` IS this
marginalization in closed form, which makes it a verification target for a
discrete polar collapse.

⛔ Homogenization/condensation flux weighting is a **solution-weighted
Petrov–Galerkin** collapse (L-010); the reduce carries no solution weight.

## The slab is a FREE ORACLE (the best first test in the attack)

A raw 3-D rule on the **slab** arm is not ill-posed at all — it redundantly
solves the same 1-D problem per duplicate `μ_x` and sums with the correct total
weight, so the answer is CORRECT. ⟹ `solve(slab, Q).flux` vs
`solve(slab, reduce(Q)).flux` is an end-to-end equivalence gate that exists
today, and it RED-s on exactly the three errors that matter: a missing `1/2π`
(off by 2π), a wrong axis (a different converged rule), a mishandled
multiplicity. Same criterion explains why the sphere has admissibility
conditions and the slab has none: the connection term exists iff the reduction
group acts on the fiber (the #326 criterion).
