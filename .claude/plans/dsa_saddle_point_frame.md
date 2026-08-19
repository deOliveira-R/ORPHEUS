# DSA Consistency ↔ Saddle-Point / Mixed-FEM / CFD Pressure–Velocity Coupling

**A cross-domain frame memo.** Detection task, not a review. Each
hypothesis gets a structural verdict — VALIDATED / REFINED / REFUTED —
with page-cited evidence, a concrete reformulation, and a first test
that can RED. Written to be cited from #312 (LD-consistent arm), #314
(2-D DSA), #294 (diffusion follow-ups), #200 (Krylov preconditioner),
and a new k-skeleton design issue.

Author: cross-domain-attacker. Grounded on the merged DSA #2 campaign
(`docs/theory/methods/sn/acceleration.rst`, the R4 ruling), the D2
characterization (`.claude/plans/archive/dsa_d2_characterization.md`), the
literature memo (`.claude/plans/archive/dsa_literature_memo.md`), and the
local literature sidecars in `scratch/literature_ocr/`.

> **Sign-erratum guard (carried from the campaign):** Alcouffe (1977)'s
> printed discrete pairs (17)/(23) carry sign errata. This memo never
> transcribes them; the transcription reference is Larsen (1982) part I.

---

## 0. One-paragraph thesis

DSA consistency and mixed-FEM inf-sup stability are the **same
theorem** wearing two vocabularies. The transport pair (scalar flux
φ, current J) is a discrete **saddle-point / mixed system** of
**Darcy type** (not Stokes type). "Consistent DSA" is precisely
"the low-order operator is the **Schur complement of a compatible
discrete pairing**," and this is already how the ORPHEUS codebase
states it (`acceleration.rst`: the low-order operator is "the Schur
complement of a two-moment (ℓ≤1) Galerkin triple product
R₁·A_high·P₁ on the *assembled* DD operator"). Reed's divergence is a
**staggering mismatch between two individually-compatible pairings**,
NOT a primal checkerboard — transport is *born staggered* (DD is the
MAC grid of transport) and structurally cannot exhibit the collocated
checkerboard that motivates CFD stabilization. The frame validates 6
of the 8 hypotheses, refines 1 (Darcy not Stokes), and turns 1 into a
precise refutation-with-reason (no primal checkerboard). It hands the
codebase four concrete pollination imports and a k-skeleton design axis.

---

## 1. Structural features (extracted, enumerated)

From the transport within-group problem + its DSA low-order operator:

- **Mixed / first-order system.** Two coupled fields: scalar flux
  φ (0th angular moment) and current J (1st angular moment), coupled
  by (i) a *balance* law `∇·J + Σa φ = q` and (ii) a *constitutive*
  law (Fick) `D⁻¹J + ∇φ = 0`. This is a `[[A, Bᵀ],[B, C]]` block
  system — a **saddle point**.
- **Operator-block content.** The (J,J) block is a **mass matrix**
  `D⁻¹` (algebraic, local, zeroth-order — trivially invertible). The
  (φ,φ) block is the removal `Σa = Σt(1−c)` (zeroth-order, ≥0,
  **vanishing as c→1**). The off-diagonals are discrete `∇·` and `∇`.
- **Two staggerings coexist in the codebase.** (a) The DSA-derived
  low-order: φ on the **(K+1) edges / 0-skeleton**, currents at
  **cell centers**; (b) the landed diffusion module
  (`orpheus.diffusion`, #290): φ at **cell centers / d-cells**,
  compact face-difference currents at **(d−1)-faces**. These are
  *opposite* staggerings (D2 Part A).
- **Reduction / Schur structure.** The low-order operator is *derived*
  by moment-reduction of the assembled high-order operator, not chosen
  — a `R₁·A_high·P₁` Galerkin triple product with an ℓ=1 Schur
  elimination (`acceleration.rst`, sn-dsa-restriction-prolongation).
- **Fixed-point + Krylov outer structure.** SI+DSA ≡ preconditioned
  Richardson; Krylov+DSA ≡ GMRES with `M = (I+C)∘(L+C)⁻¹`
  (`acceleration.rst`, sn-dsa-both-postures; Century §1.9).
- **A vanishing corrector.** The DSA correction source is the
  scattering residual `c(φ^{l+1/2} − φ^l)`, which is **identically
  zero at the fixed point** — so a wrong accelerator degrades the
  *rate*, never the *answer* (`acceleration.rst`,
  sn-dsa-correction-vanishes).
- **Non-monotone dispersion defect.** The consistent low-order's
  off-diagonal **flips sign** at thick cells (σₜh≳√(4D/σR)·σₜ),
  deliberately abandoning M-matrix/monotone form (D2 Part B).
- **Boundary trace as half-range moment.** The Marshak BC is the
  half-range (`|Ω·n|w`-weighted) ℓ=0 moment closure — "the boundary
  analog of Fick" (`acceleration.rst`, sn-dsa-boundary-rows).
- **Skeleton-homed cochain.** The face flux is a `(d−1)`-cochain
  (`wavefront_cochain.rst`, DEC frame); the DSA scalar unknown is a
  0-cochain (edges/corners); currents live on d-cells.
- **An angular induced-mesh sibling.** Curvilinear SN carries an
  **angular** conservation law `(1−µ²)/r ∂_µ` on a half-angle-face
  induced mesh with no-flux poles (`radial_characteristic_space.py`).

## 2. Elegance-detector hits

- **Smell #16, shape 1 (two paths to one discrete operator over
  different storage conventions).** The DSA low-order (edge-homed) and
  the diffusion-module operator (cell-homed) are two discretizations
  of "the diffusion operator for φ." The R4 ruling already resolved
  this as *coexist-by-defining-law* (they are different operators, not
  a twin), but the frame explains WHY it is not a Smell-#16 collapse:
  they are **two non-isomorphic inf-sup-stable pairings**, and only
  one is the Schur complement of the *sweep's* pairing. The
  collapse-onto-one-primary-object move (L-003) here is **onto the
  reduction**: the consistent operator is not a free choice of
  staggering, it is forced by "reduce-the-discrete."
- **Smell: a repeated conditional / staggering choice that is really a
  missing type.** "Which skeleton does this DoF live on?" is currently
  answered implicitly per-operator (edge vs cell) rather than as a
  declared property of the field. This is the k-skeleton design axis
  (hypothesis g).
- **No smell fired on the corrector-vanishes property** — that is not
  a code smell, it is the load-bearing correctness-safety theorem, and
  it is the transport analog of the CFD "converged pressure-correction
  → 0" fact (see 4b, 4e).

## 3. The core correspondence dictionary

The frame is a functor between three columns. **Read the columns as
the same block operator with different operator-content in the blocks.**

| Transport / DSA | Mixed-FEM (Darcy / mixed-Poisson) | CFD incompressible (Stokes/NS) |
|---|---|---|
| scalar flux φ (ℓ=0 moment) | pressure/potential p | pressure p |
| current J (ℓ=1 moment) | flux u ∈ H(div) | velocity u |
| balance `∇·J + Σa φ = q` | mass `∇·u + c·p = f` | continuity `∇·u = 0` |
| Fick `D⁻¹J + ∇φ = 0` | Darcy `K⁻¹u + ∇p = 0` | momentum `−νΔu + ∇p = f` |
| (J,J)-block: **mass** `D⁻¹` | (u,u)-block: **mass** `K⁻¹` | (u,u)-block: **Laplacian** `−νΔ` |
| (φ,φ)-block: removal `Σa≥0`, →0 as c→1 | reaction `c≥0` | **0** (pure constraint) |
| consistent DSA low-order | Schur complement `S = BA⁻¹Bᵀ + C` | pressure Schur complement |
| DD staggering (φ edges, J cells) | RT0 / covolume compatible pair | MAC staggered grid |
| inconsistent (Reed) cell-centered | inf-sup-INCOMPATIBLE pair | collocated equal-order |
| SI + DSA | (stationary Uzawa/Schur iteration) | SIMPLE pressure-correction |
| Krylov + DSA preconditioner | block-preconditioned GMRES/MINRES | Newton–Krylov / GMRES-SIMPLE |
| M4S / partially-consistent DSA | (stabilized nearly-compatible pair) | Rhie–Chow / Brezzi–Pitkäranta |

**The load-bearing refinement** (developed in 4a): the correct import
column is the **middle one (Darcy / mixed-Poisson / H(div))**, not the
right one (Stokes). The (J,J) block is a *mass matrix*, not a
Laplacian, and the (φ,φ) block is a *nonzero (or vanishing) reaction*,
not an exact zero. That distinction changes which stable-pairing
catalog and which inf-sup theory apply.

---

## 4. Per-hypothesis verdicts

### 4a. Mixed-pair correspondence (φ,J) ≅ (p,u) — VALIDATED, then REFINED to Darcy

**Trigger:** first-order mixed system with a balance law + a
constitutive law over two co-defined fields (structural feature 1) —
the saddle-point signature `[[A,Bᵀ],[B,C]]`.

**VALIDATED core.** The P1/diffusion mixed form IS a discrete
saddle-point system, and the (φ,J) ↔ (p,u) map is exact at the level
of *roles*: φ/p is the potential the constraint acts on; J/u is the
flux the constitutive law produces; balance ≅ continuity;
Fick ≅ the constitutive (Darcy/momentum) relation. The requirement
that the (φ,J) discrete spaces form a **compatible (inf-sup / LBB)
pair** is the same requirement in both domains. This is not analogy —
`acceleration.rst` already writes the low-order operator as a Schur
complement of a two-field Galerkin reduction, which is the *definition*
of a mixed method.

**REFINED — the operator content is Darcy, not Stokes.** The
hypothesis text pairs against "incompressible (p,u)" and "momentum."
That is the wrong import column, and the difference is load-bearing,
not cosmetic:

1. **The (J,J) block is a mass matrix `D⁻¹`, not a vector Laplacian
   `−νΔ`.** Consequence: the Schur complement `S = B·A⁻¹·Bᵀ + C` is a
   genuine *sparse local* diffusion operator (because A⁻¹ is local for
   a mass matrix), so the ℓ=1 (current) block can be **Schur-eliminated
   in closed form** — which is exactly Larsen's four-step step 4
   (`acceleration.rst`, eliminate f₁ onto f₀). In Stokes, A⁻¹ (inverse
   Laplacian) is *nonlocal*, the pressure Schur complement is
   spectrally a **mass/identity** operator (not a diffusion operator),
   and no closed-form elimination exists. **The transport current
   elimination is only possible because the constitutive block is
   algebraic.**
2. **The (φ,φ) block is `Σa = Σt(1−c) ≥ 0`, not an exact zero.** In
   pure Stokes the pressure is a Lagrange multiplier with `C=0`; the
   pressure is determined *only* through the inf-sup condition and
   admits a constant null mode. In transport, `Σa > 0` (c<1)
   **regularizes** the saddle point — the φ-block is coercive, so even
   a marginal pairing is damped. This is a partial structural reason
   transport tolerates staggerings CFD would call unstable (see 4f).
   **But** as c→1 (the DSA-relevant limit), `Σa→0` and the system
   approaches the pure-constraint (Stokes-like) structure — which is
   exactly the limit where acceleration is hard and inf-sup
   compatibility becomes load-bearing.

**Elegance payoff:**
- *Structure-exposing:* names the transport mixed system as a member
  of the **H(div) / mixed-Poisson** family, which makes "RT0-P0 is the
  stable lowest-order pair" an *importable theorem* rather than a
  transport-specific discovery. It also exposes that DD's edge-φ /
  cell-J arrangement is a **covolume** pairing (0-skeleton ↔ dual
  0-skeleton), a *different* stable family than RT0 — see 4g.
- *Structurally-simpler:* one saddle-point object subsumes "balance +
  Fick + boundary Marshak" as `[[A,Bᵀ],[B,C]]` + a trace block, which
  is the #208 biproduct algebra already in the codebase (L-011: the
  coupled block system is a re-association of an existing biproduct,
  not a new object).

**First test (discriminates Darcy from Stokes import):** take the
lowest-order **Stokes-stable** pair (Taylor–Hood P2-P1: continuous
piecewise-quadratic velocity, continuous piecewise-linear pressure)
and use it to assemble a transport low-order operator; measure its DSA
spectral radius on the D2 fixture (uniform 1g slab, c∈{0.9,0.99},
σₜh∈{0.1…30}). **Prediction:** Taylor–Hood is NOT an RT0/covolume
pairing and does not reproduce the sweep's staggering, so it is an
*inconsistent* accelerator — ρ degrades/diverges at thick cells like
the cell-centered operator (D2 Part C, ρ up to 54.7), *not* the flat
0.14–0.18 of the derived edge operator. A naive "any inf-sup-stable
pair works" claim would predict it stays bounded — that is the
falsifiable divergence. (The deeper point the test pins: inf-sup
stability is necessary for a *standalone* solve but NOT sufficient for
a *consistent accelerator* — consistency additionally requires the
pairing to be the sweep's own reduction. See 4c.)

**Structural attack on current:** the R4 ruling records *that* the
diffusion module's operator diverges as an accelerator, with the
measured numbers. The frame supplies the *why in importable form*: the
diffusion module is a valid H(div)-mixed (RT0-family) discretization
chosen for standalone accuracy; it is inf-sup stable **and** wrong as
an accelerator, which proves inf-sup stability ≠ consistency. Naming
this makes the R4 coexistence a theorem, not a measured coincidence.

### 4b. Reed's divergence re-diagnosed — VALIDATED (staggering mismatch), and it is NOT a primal checkerboard

**Trigger:** two individually-well-posed discrete operators (the DD
sweep's internal diffusion limit; the cell-centered diffusion operator)
whose thick-cell limits **disagree at O(1)** (structural feature 3 +
the D2 divergence table).

**VALIDATED.** The hypothesis is correct and is stated almost verbatim
in the codebase. `acceleration.rst` (sn-dsa-consistency-is-derived,
citing Alcouffe 1977 p. 348): "the DD P1-moment equation couples the
**cell-center current** to **cell-edge scalar-flux differences** — a
staggered Fick's law — while the cell-centered diffusion operator
evaluates the current at cell edges. The two discrete diffusion limits
disagree by O(1) at finite h, so the correction source stays finite
exactly where it should vanish, and the accelerator overshoots into
instability." That is a **staggering mismatch between two operators**,
each individually a valid diffusion discretization:
- The cell-centered diffusion operator is a perfectly stable,
  monotone, M-matrix standalone solver (it solves diffusion problems
  correctly).
- The DD sweep's *implied* diffusion is the edge-φ / cell-J operator.
- Neither is unstable alone; **their composition in the SI outer loop**
  is what diverges (D2 Part C: cell-centered-as-accelerator hits
  ρ=54.7 at c=0.99, σₜh=30, while the derived edge operator is flat
  ≤0.181 over the whole grid).

**The CFD analogue is exact and is the strong form of the hypothesis:**
a **projection / SIMPLE-type pressure-correction** whose pressure-
Poisson operator is assembled on a staggering inconsistent with the
momentum discretization. If the discrete `div` and `grad` used in the
pressure equation are not the adjoint pair implied by the momentum
step, the "projection" does not actually project onto the discrete
divergence-free space the momentum operator sees — the corrected
velocity is not discretely divergence-free in the pressure solver's
metric, the correction source does not vanish at steady state, and the
outer loop stalls or diverges. This is the well-documented
"inconsistent pressure Laplacian" failure of naive fractional-step /
SIMPLE schemes. **Same mechanism, same fix** (make the pressure
operator the Schur complement of the *momentum* discretization — 4c).

**REFUTED sub-claim, promoted to its own hypothesis (4f):** Reed is
NOT "a primal checkerboard instability." A checkerboard is a *null
mode of a single collocated operator*. Reed's divergence is an
*inter-operator* thick-cell disagreement; both operators are
nonsingular. The hypothesis correctly says "NOT a primal checkerboard"
— the frame confirms it and (4f) shows transport structurally cannot
have the primal form.

**Elegance payoff:**
- *Structure-exposing:* recasts "Reed diverges" from a historical fact
  into an instance of the general **inconsistent-Schur-complement**
  failure that CFD, mixed-FEM, and DSA all share. The corrector-
  vanishes property (feature 6) is the shared safety net: in all three,
  a *converged* answer is correct regardless of the corrector, because
  the correction source → 0 at the fixed point; only the *rate* is at
  risk. This is why in all three domains you can get away with an
  approximate/mismatched corrector at the cost of iterations — until
  the mismatch pushes an eigenvalue past 1 and the fixed point diverges.
- *Algorithmic-advantage:* the diagnosis tells you the fix is not
  "stabilize φ" but "derive the low-order as the sweep's own
  reduction" — which is a different (and cheaper, closed-form)
  operation than adding a stabilization term.

**First test (discriminates staggering-mismatch from a hypothetical
primal instability):** run the cell-centered diffusion operator as a
**standalone solver** on a thick-cell (σₜh=30) diffusion problem and
verify it converges to the correct monotone answer with a well-
conditioned SPD matrix (no oscillation, no null mode); then run the
*same* operator as a DSA accelerator on the *same* mesh and watch ρ>1.
**Prediction:** the operator is stable standalone and divergent as an
accelerator — proving the pathology is in the *pairing with the sweep*,
not in the operator's own spectrum. A "primal checkerboard" hypothesis
predicts the operator would show a null mode / spurious oscillation
*standalone* — it does not. (The D2 Part C numbers already contain the
accelerator half; the discriminator adds the standalone-stability half.)

**Structural attack on current:** none needed — the campaign already
holds this diagnosis. The frame's contribution is to make it portable:
the same test discriminates a mismatched pressure-Poisson operator in
any future coupled/multiphysics solver ORPHEUS builds (#294 diffusion
follow-ups, any future CFD-adjacent coupling).

### 4c. Consistency-by-reduction = inherited compatibility — VALIDATED (this is the deepest and most elegant node)

**Trigger:** the low-order operator is *computed* as a moment-reduction
/ Schur complement of the assembled high-order operator, not chosen
(structural feature 4; `acceleration.rst` "reduce-the-discrete, never
discretize-the-reduced").

**VALIDATED, and it is the load-bearing structural identity of the
whole frame.** The claim: *the Schur complement of a well-posed
(compatible) discrete pairing inherits the compatible staggering for
free; choosing the low-order independently re-opens the compatibility
problem.* This is exactly right and the codebase proves it:

- The consistent low-order operator is `Schur[ R₁·A_high·P₁ ]` — the
  ℓ≤1 Galerkin reduction of the *assembled* DD operator followed by
  elimination of the current (ℓ=1) block (`acceleration.rst`,
  sn-dsa-restriction-prolongation; L-009 corollary: the "consistent
  low-order operator IS the Galerkin coarse operator R·A_high·P
  post-composed with a Schur complement of the retained-but-closed
  moments — Fick = odd-block Schur; Marshak = incoming-partial-current
  Schur").
- **Why inheritance holds:** the Schur complement of a saddle point
  with an inf-sup constant β has its own coercivity bounded below by β²
  (Brezzi theory). If A_high is a compatible discretization, its
  moment-reduction is automatically compatible — the staggering is not
  re-chosen, it is *projected through*. This is why the derived edge
  operator needs **no harmonic-mean interpolation** (D2 Part A: "an
  edge unknown straddles a material-homogeneous cell") while the
  independently-chosen cell operator **requires** the harmonic mean
  (its unknowns straddle material faces): the reduction inherits the
  right unknown-home; the independent choice must repair the wrong one.
- **The exactness anchor.** D2's S2 anchor is the theorem's sharpest
  witness: for S2-GL the moment system closes *exactly* (the angular
  space IS span{1,μ}), so the ℓ≤1 reduction is the **exact** Schur
  complement, and consistent DSA converges **in one iteration**
  (measured ρ ≈ 3×10⁻¹⁵). An exact Schur complement is an exact
  solver — the spectral-radius-0 limit of "inheritance." This is the
  discriminating number that a "choose the low-order independently"
  scheme can never hit (it has no reason to be the exact reduction at
  S2).

**Elegance payoff (this node hits all four criteria):**
- *Structure-exposing:* "consistency" — the single most-abused word in
  the DSA literature — is revealed to be a one-word name for "the
  low-order is a Schur complement of the high-order." The mystique
  dissolves into a standard mixed-method identity.
- *Expressive:* the same sentence covers interior (Fick = odd-block
  Schur), boundary (Marshak = incoming-partial-current Schur), and the
  P1 arm (the ℓ=1 update is the *retained* Schur block). Three
  historically-separate constructions become one operation over three
  measures.
- *Structurally-simpler:* removes the degree of freedom "which
  diffusion operator?" — there is only one, the reduction. The R4
  coexistence with `orpheus.diffusion` is then not a twin but a
  *different reduction target* (standalone accuracy vs sweep-matching).
- *Algorithmic-advantage:* the reduction is closed-form (mass-block
  elimination, 4a), so the consistent operator costs no more to build
  than the inconsistent one — the "expensive consistency" folklore is
  false for the mass-block (Darcy) case.

**First test (discriminates reduction-inheritance from independent
choice):** build the low-order operator two ways for the SAME S4-GL
fixture — (route α) the moment-reduction of the assembled DD operator;
(route β) the *analytic* P1 diffusion operator `−∇·D∇ + Σa` discretized
with the *same nodal layout* (edges) as route α, using the continuum
`D=1/(3Σt)` and lumped removal. **Prediction:** route α carries
h-baked closure coefficients (the ¼(1,2,1) mass, the σR-fold damping
`D=1/[3(Σt−Σs1)]`, γ_N≈0.2606 not ¼) and is flat-ρ; route β uses the
"nice" continuum constants and degrades — even though it lives on the
*same edges*. This isolates the claim: it is not the *skeleton home*
alone that gives consistency (β has the right home and still fails),
it is that the coefficients are the *projection of the high-order
operator*. A test that only compared cell-home vs edge-home would pass
route β falsely; comparing same-home reduced-vs-analytic is the
discriminator. (The campaign's `test_dd_instance_coefficients` and
`test_low_order_matches_reference_builder` pin route α entry-for-entry;
route β is the negative control this frame recommends adding to #314's
2-D planning as the "same-skeleton-wrong-coefficients" control.)

**Structural attack on current:** the current build proves route α ≡
Larsen (27). It does not yet ship the same-skeleton route-β negative
control that isolates *coefficient-inheritance* from *skeleton-home*.
The two are conflated in the derived-vs-landed comparison (which
differs in BOTH home and coefficients). The frame demands the
one-variable-at-a-time control.

### 4d. Rhie–Chow analogue = M4S / partial-consistency family — VALIDATED, with a three-way sharpening

**Trigger:** a stabilization that keeps a *convenient collocation* and
adds an *exactly-scaled compact-coupling correction* to recover
compatibility (the M4S "modified four-step"; the partially-consistent
DSA family; the D2 Part B sign-flip requirement).

**VALIDATED, refined into a clean three-way correspondence.** The
hypothesis pairs Rhie–Chow ↔ M4S. The frame sharpens this into a
*three-tier* ladder that both domains share exactly:

| Tier | Transport DSA | CFD pressure–velocity |
|---|---|---|
| **fully consistent (born staggered)** | derived edge operator (Schur of the sweep) | MAC staggered grid (exact div/grad adjoints) |
| **partially consistent (stabilized collocation)** | M4S / partially-consistent DSA | Rhie–Chow collocated / Brezzi–Pitkäranta / PSPG |
| **inconsistent (naive)** | Reed cell-centered | naive collocated equal-order (checkerboard) |

The middle tier is where the Rhie–Chow ≅ M4S pairing lives, and both
share the defining property: **the correction must be *exactly scaled*
by the operator being stabilized.** In Rhie–Chow, the face-velocity
correction is scaled by the momentum-equation diagonal `1/a_P` — get
the scaling wrong and you either under-stabilize (checkerboard returns)
or over-stabilize (grid-dependent artificial dissipation that does not
vanish under refinement). In DSA, the partially-consistent closure
must reproduce the sweep's dispersion *including its defects*.

**The D2 sign-flip is the sharp witness of "exactly scaled," and it is
counterintuitive in exactly the way the frame predicts.** D2 Part B:
the consistent ¼(1,2,1) removal mass **flips the off-diagonal sign** at
σₜh≳√(4D/σR)·σₜ, abandoning the M-matrix/monotone form; the lumped
M-matrix variant is the *inconsistent* class. Structural reading:
**the accelerator's job is to match the operator being accelerated,
NOT to be a good elliptic solver.** A "better-behaved" (M-matrix,
monotone) diffusion operator is a *worse* accelerator because it does
not reproduce the DD sweep's own (non-monotone) thick-cell dispersion.
The CFD twin: a symmetric "nice" pressure Laplacian is a *worse*
Rhie–Chow stabilizer for a convection-dominated momentum operator than
one that inherits the momentum operator's (non-symmetric,
non-M-matrix) coefficients — because the stabilization must cancel the
*actual* odd-even decoupling the momentum interpolation produces, not
an idealized one.

**The dose-sensitivity twin is exact.** The hypothesis names
McCoy–Larsen Table II (a partially-consistent closure diverges at a
threshold that *migrates with the inconsistency*). `acceleration.rst`
(sn-dsa-what-was-tried) records both failure modes precisely:
- a cell-*edge* inconsistent scheme "degrades gracefully toward ρ→ρ_SI
  — convergent but ineffective, it performs like source iteration"
  (Adams–Larsen eq. 3.43);
- a *partially* consistent scheme (right operator shape, wrong closure
  constant) "diverges at a threshold that merely *moves* to thicker
  cells as the inconsistency shrinks" (McCoy–Larsen Table II).

This IS the stabilization-parameter story: too little Rhie–Chow /
Brezzi–Pitkäranta stabilization → checkerboard returns; the threshold
migrates continuously with the parameter. "Inconsistency has no safe
dose" is the transport statement of "there is no stabilization
parameter that is uniformly safe across cell Reynolds/optical
thickness."

**Elegance payoff:**
- *Structure-exposing:* unifies "M4S," "partially-consistent DSA,"
  "Rhie–Chow," and "Brezzi–Pitkäranta stabilization" as one object —
  *an exactly-scaled compact correction on a convenient collocation* —
  and names the shared failure surface (dose-threshold migration).
- *Expressive:* the three-tier ladder replaces a scatter of named
  schemes with a single axis (consistency level) × two domains.

**First test (discriminates exact-scaling from a constant stabilizer):**
implement a *partially*-consistent low-order for the D2 fixture by
replacing the derived σR-fold denominator `1/[1+(3/2)ρ(Σt−Σs0)h]` (the
WD `ρ`-damping that bakes the closure weight) with a *fixed* constant
(e.g. the continuum `ρ=0`, i.e. the DD member) applied to a WD sweep
(`α≠0`). **Prediction:** the partially-consistent operator converges
for thin cells and *diverges past a threshold σₜh\* that moves to
larger σₜh as the mismatch (`α`) shrinks* — the McCoy–Larsen Table II
signature — whereas the fully-derived WD operator stays flat at all
σₜh. An "any reasonable constant works" implementation predicts a
uniform bounded ρ; the migrating threshold is the RED. (This is the
weighted-diamond partial-consistency negative control the campaign
lists as a "three discoveries" artifact — the frame confirms it is the
Rhie–Chow dose-sensitivity twin and recommends it as the reusable
negative-control *pattern* for #312/#314.)

**Structural attack on current:** the P0/P1 arms are fully consistent;
the codebase carries the WD partial-consistency case only as a
*negative control*, not as a shipped scheme. The frame says: when #312
(LD arm) forces a *partially*-consistent M4S (full LD-consistency is
"structurally unspellable without the M4S reduction," `acceleration.rst`
sn-dsa-honest-scope, R5a), the Rhie–Chow literature is the design
oracle for *how much* stabilization and *how to scale it* — import it
rather than re-deriving the dose threshold empirically.

### 4e. Krylov as iteration-level stabilization — VALIDATED, and the two failure modes are GENUINELY DISTINCT

**Trigger:** two different "eigenvalue-near-or-past-1" situations
(inconsistency divergence; consistent-but-discontinuity-degraded) both
"rescued" by a Krylov wrap (structural feature 5; Century §1.9.2–1.9.3).

**VALIDATED — the two modes are distinct mechanisms, and the Century
chapter names both precisely.** This is the hypothesis's real payload:
keep them apart. The frame confirms they are different and supplies the
Century chapter's own two-mechanism vocabulary.

**Mode (i) — inconsistency divergence (Reed-class).** As a *fixed
point*, an inconsistent/unstable accelerator has an iteration-matrix
eigenvalue with modulus >1 (high-frequency amplification) → diverges.
Recast as a Krylov *preconditioner*, that same operator is — in the
Century chapter's exact words (p. 81) — "analogous to **unstable
acceleration schemes that strongly attenuate the low-frequency Fourier
error modes but weakly amplify the high-frequency** ones"; as a
preconditioner it "moves the smallest eigenvalues away from zero but
also significantly increases and spreads out (unclusters) the largest
eigenvalues," and CG/GMRES can still net-decrease the effective
condition number → **converges**. **Bound (the frame insists on the
Century caveat):** this rescue is NOT unconditional — "if the
amplification factors for the high-frequency modes are *too large*,
the preconditioner may increase the largest eigenvalues too much …
resulting in a net decrease in convergence rate" (p. 81). So Krylov
converts *mild* inconsistency-divergence into slow convergence; *wild*
inconsistency can still defeat it. Mechanism: **spectral-radius
irrelevance** — GMRES cares about the whole spectrum, not the modulus
of the largest eigenvalue.

**Mode (ii) — consistent-DSA degradation at material discontinuities.**
Even a *consistent* DSA leaves a **few** eigenvalues of the
iteration matrix near 1 in multi-D with strong cross-section
discontinuities — the interface/boundary-layer modes that the *local
diffusion model* cannot represent (the exact transport Schur complement
is *nonlocal* at a discontinuity; the P1/diffusion truncation misses it
in a few localized modes). As a fixed point, this is fatal: Century
p. 83, "**only one eigenvalue of the iteration matrix very near unity
can ruin the performance of an acceleration scheme**." As a Krylov
method it is benign: "**but a single anomalously large eigenvalue
associated with the system being solved will generally have little
effect on a Krylov method**" — GMRES deflates the handful of interface
outliers in a few extra iterations. Hence Century p. 83: "the
deficiency of DSA in multidimensions with large discontinuities … is
**strongly diminished** when it is recast as a preconditioned Krylov
method [124]." Mechanism: **outlier tolerance** — a few bad eigenvalues
cost a few iterations, not convergence.

**Why they are genuinely distinct (the frame's structural claim):**
- Mode (i) is a **discretization-consistency** failure — the low-order
  operator is the *wrong discrete reduction* of the high-order operator
  (wrong staggering/coefficients). It is fixable *without Krylov* by
  getting the reduction right (4c); Krylov merely tolerates the wrong
  one.
- Mode (ii) is a **model-adequacy** failure — the low-order operator is
  the *right discrete reduction of the wrong continuous model*
  (diffusion is not the transport Schur complement at a sharp
  interface). It is *NOT* fixable by consistency alone; the diffusion
  model is intrinsically blind to the interface modes. Krylov is the
  *only* cheap fix (short of a nonlocal low-order model).

They differ in *what is wrong* (the discretization vs the model), in
*whether consistency fixes it* (yes vs no), and in *how Krylov helps*
(spectral-radius irrelevance vs outlier deflation). Conflating them
would predict that fixing consistency cures the multi-D discontinuity
problem — it does not; that is precisely the [124] result.

**Local corroboration of mode (ii):** the Warsa–Wareing–Morel–McGhee–
Lehoucq (2004) k-eigenvalue Krylov paper (local sidecar) and the
Century chapter's [124] pointer are the two local secondaries for the
NSE 147:218 primary (which is NOT local — flagged in §8). The Century
chapter additionally supplies the LD-on-tets fact (p. 83): the SAPD
(self-adjoint-positive-definite) property of Peierls' operator "is
**not preserved** with the linear-discontinuous spatial discretization
scheme on unstructured tetrahedral meshes [124]," and a discontinuous
diffusion discretization makes the DSA-preconditioned operator
"nonsymmetric or symmetric-indefinite," precluding CG. This is a third,
distinct degradation (definiteness loss) — see 4f and pollination for
#312.

**Elegance payoff:**
- *Structure-exposing:* replaces "Krylov helps DSA" with two named,
  disjoint mechanisms (spectral-radius-irrelevance for consistency
  failures; outlier-deflation for model-adequacy failures), each with
  its own bound.
- *Algorithmic-advantage:* tells the implementer *which* failure a
  given problem has and therefore *whether* fixing consistency will
  help. For #200 (Krylov preconditioner, now the re-enabled slot) and
  #314 (2-D DSA, where discontinuity-degradation appears), this is the
  difference between "spend effort on a more consistent 2-D operator"
  (helps mode i, wasted on mode ii) and "wrap in GMRES" (helps both,
  but is the *only* help for mode ii).

**First test (discriminates the two modes):** on a 1-D two-material
slab with a **strong cross-section discontinuity** (optically thick
absorber against a thick scatterer), measure the DSA iteration-matrix
spectrum two ways — (A) with the *consistent* derived edge operator,
(B) with an *inconsistent* cell-centered operator. **Prediction:**
- (B) shows a *band* of eigenvalues pushed past 1 (mode i, bulk
  amplification) → SI+DSA diverges, Krylov+DSA converges slowly but
  monotonically (net-conditioning still improves unless the band is
  extreme);
- (A) shows the bulk spectrum tightly clustered near 0 but **a few
  isolated eigenvalues creeping toward 1 localized at the material
  interface** (mode ii) → SI+DSA convergence *degrades* (ρ→ ~1) but
  never diverges; Krylov+DSA restores fast convergence by deflating the
  few outliers.
The discriminator is the *shape* of the near-1 part of the spectrum: a
**band** (mode i, consistency) vs **isolated interface outliers**
(mode ii, model adequacy). An implementation that treats them as one
mechanism predicts the same spectral shape for both — the differing
shapes are the RED. (This is a direct, cheap 1-D precursor to the
multi-D [124] experiment and is worth adding to #314's evidence plan.)

**Structural attack on current:** the campaign is 1-D slab, where mode
(ii) barely appears (material-discontinuity degradation is a genuinely
multi-D effect — `acceleration.rst` sn-dsa-honest-scope defers 2-D).
The codebase therefore has *no* test that exhibits mode (ii), and the
D13 iteration-count table shows the Krylov wrap "reducing counts by
only one or two" (the *consistent-scheme* regime, Century p. 110–111)
— which is exactly the regime where Krylov looks unimpressive. The
frame's warning for #314: **the Krylov posture will look nearly useless
in 1-D and become essential in 2-D-with-discontinuities**, and the
reason is mode (ii), which 1-D cannot show. Do not conclude from the
1-D "+1 or +2 iterations" that the Krylov slot (#200) is low-value.

### 4f. Primal-checkerboard hunt — REFUTED, with a precise structural reason

**Trigger checked:** a null mode of a single collocated operator (an
odd-even/checkerboard mode that produces zero discrete gradient and is
therefore unconstrained by the momentum/constitutive equation).

**REFUTED as a phenomenon of the transport schemes ORPHEUS carries.**
There is no genuine primal null-mode checkerboard in DD, in the derived
edge operator, or in the diffusion module. The structural reason is
threefold and each part is falsifiable:

1. **Transport is *born staggered* — DD is the MAC grid of transport.**
   The DD scheme homes the scalar flux on the 0-skeleton (edges/corners)
   and the current on d-cells: the *maximally-separated* skeleton pair
   (structural feature 3). The collocated checkerboard is a pathology of
   putting φ and J on the *same* skeleton with a null-space-admitting
   interpolation. DD never collocates — the staggering that CFD must
   *engineer* (MAC) or *stabilize into* (Rhie–Chow), transport gets *for
   free* from the angular-moment structure (the current is literally the
   ℓ=1 moment, homed where the sweep produces it). So the arrangement
   that would host a checkerboard is structurally absent.
2. **The removal term `Σa ≥ 0` coercivizes the φ-block** (4a, refinement
   2). A checkerboard null mode requires the (φ,φ) block to be a pure
   constraint (`C=0`); transport's is `Σa = Σt(1−c)`, strictly positive
   for c<1. The φ-block has no kernel. (Caveat, honest: as c→1, `Σa→0`
   and the coercivity vanishes — but part (1) still holds, so even the
   pure-scatterer limit is staggered-stable, which is exactly why the
   D2 table shows the *derived* operator flat even at c=0.99.)
3. **What thick cells actually produce is a *dispersion* defect, not a
   *null-mode* defect.** The DD amplification factor per cell is
   `(1−σₜh/2)/(1+σₜh/2) → −1` as σₜh→∞: the solution *oscillates in
   sign* (the classic DD negative-flux problem), but the operator is
   **nonsingular** and the solution is **uniquely determined**. The
   derived edge operator's off-diagonal **sign flip** at thick cells (D2
   Part B) is the same category: non-monotone (loses M-matrix form) but
   still invertible. **Oscillatory ≠ null.** A checkerboard is an
   *unconstrained* mode (kernel); a DD oscillation is a *determined*
   oscillation (no kernel).

**The single closest thing to a "spurious mode" in transport is a
different defect: definiteness loss for LD-on-tets.** Century p. 83
(citing [124]): the SAPD property "is **not preserved** with the
linear-discontinuous spatial discretization scheme on unstructured
tetrahedral meshes," and a discontinuous diffusion discretization makes
the DSA-preconditioned operator "nonsymmetric or symmetric-indefinite."
That is a **loss of definiteness** (the operator gains negative
eigenvalues), NOT a kernel — the operator stays nonsingular; it just
can no longer use CG (must use GMRES). So even the worst multi-D DFEM
case is an indefiniteness defect, not a checkerboard.

**What the thick-diffusion-limit asymptotics papers actually establish
(be precise — this is the hypothesis's explicit ask).** Larsen–Morel–
Miller (1987) and Larsen–Morel (1989) establish **diffusion-limit
accuracy/existence**, NOT well-posedness or null-mode absence. Verbatim
scope (LMM 1987 abstract, p. 283): "If the result of either expansion
is a legitimate diffusion description for either the cell-average or
cell-edge fluxes, then we say that the appropriate flux has the
appropriate diffusion limit; otherwise, it does not… differencing
schemes that do have a particular diffusion limit are substantially
more accurate, in the regime described by the limit, than those that do
not." So LMM answers *"does the discretized SN operator limit to a
legitimate discrete diffusion operator (and for which flux — cell-
average or cell-edge)?"* — an **accuracy** question — and it is exactly
the theoretical backbone of the *staggering* story (the cell-average
vs cell-edge distinction IS the skeleton-home distinction of 4g).
Null-mode absence follows from the staggered arrangement + coercive
removal (points 1–2 above), which LMM does not address and does not
need to. Do not over-cite LMM as a stability/null-mode result; it is a
*diffusion-limit-accuracy* result that *underpins why the DD-implied
diffusion (edge-homed) is a legitimate diffusion discretization* — the
precondition for consistent DSA, not a statement about checkerboards.

**Elegance payoff:**
- *Structure-exposing:* names the precise category of each transport
  "instability" — dispersion defect (DD oscillation, edge sign flip),
  definiteness defect (LD-on-tets), inter-operator mismatch (Reed) —
  and shows that NONE is the CFD primal null mode. The parallel to CFD
  is a parallel of *saddle-point roles*, not of *failure modes*: CFD's
  signature failure (collocated checkerboard) has no transport
  counterpart, and transport's signature failure (inter-operator
  consistency, Reed) is the one CFD calls "inconsistent projection."
- *Structurally-simpler:* removes a false worry from the #314 2-D design
  — there is no primal-stabilization (Brezzi–Pitkäranta-on-φ) term to
  add, because there is no φ null mode to kill. The 2-D consistency work
  is entirely inter-operator (get the corner/cell staggering right),
  never intra-operator (stabilize a collocated φ).

**First test (discriminates "no primal checkerboard" from a claimed
one):** assemble two low-order operators on the same 1-D thick-cell
(σₜh=30) mesh — (A) the DD-implied edge operator (born staggered);
(B) a *deliberately collocated* transport-diffusion operator with φ and
J both at cell centers and central-difference coupling (the transport
analog of naive collocated CFD). Compute the null space of each.
**Prediction:** (A) has a **trivial** null space (nonsingular) at all
σₜh — it may produce an oscillatory solution but never an unconstrained
mode; (B) exhibits the **checkerboard null mode** (a ±1 alternating φ
in the kernel). This proves the primal form exists in transport *only
if you deliberately collocate*, which the shipped schemes never do. A
"transport has a checkerboard like CFD" claim predicts a zero
eigenvalue in (A) — there is none; that absence is the RED against the
claim. (Bonus discriminator for LD-on-tets, when #312 lands: assemble
the LD-DFEM DSA operator on a tet mesh and check the *sign* of its
eigenvalues — predict indefinite, NOT singular; a "null mode" claim
predicts a zero eigenvalue, the frame predicts negative ones.)

**Structural attack on current:** the campaign's honest-scope note
correctly defers 2-D/DFEM. The frame adds a *design constraint* for
when they land: do not import CFD *pressure-stabilization* (the
Brezzi–Pitkäranta / PSPG φ-term) into 2-D DSA — it addresses a null
mode transport does not have. Import instead the CFD *consistent-
projection* machinery (get div/grad adjoint and the corner staggering
right), which addresses the mismatch transport *does* have.

### 4g. The k-skeleton design frame — VALIDATED (DEC endorses skeleton-home as the axis), with three structural warnings

**Trigger:** a discrete field whose home is a specific mesh skeleton
(0-cochain φ on edges/corners; (d−1)-cochain flux on faces; d-cochain
density in cells), and a design question about typing that home
(structural features 8–9; `wavefront_cochain.rst` DEC frame).

**VALIDATED — DEC/cochain endorses "DoF home on the mesh k-skeleton"
as the right generalization axis.** In discrete exterior calculus every
discrete field *is* a k-cochain and its home skeleton is intrinsic to
its physics: potentials are 0-cochains, gradients/circulations are
1-cochains, fluxes are (d−1)-cochains, densities are d-cochains. The
codebase already ratified one instance — the interior face flux is a
`(d−1)`-cochain (`wavefront_cochain.rst`). The DSA scalar unknown is a
**0-cochain** (edges in 1-D, Alcouffe **corners** in 2-D — the
literature memo §1.8: "the flux moments must be evaluated on the
spatial mesh cell **corners**… the vertex-centered ancestor of the
multi-D 9-point consistent stencil"). Currents live on d-cells. So the
DSA pair is the **(0-skeleton, d-skeleton)** — maximally separated.

**The frame's sharpening: there are TWO distinct DEC arrangements for
diffusion, and the codebase carries both — they are not a twin.**
- **Covolume / nodal-dual arrangement:** φ on primal 0-cells (vertices/
  corners), J on d-cells (= dual 0-cells / cell centers), coupled
  through the primal-edge / dual-edge duality. This is the **DSA-derived**
  operator. It is the transport analog of the *finite-difference /
  covolume / mimetic* family.
- **RT0 / mixed-FEM arrangement:** φ on d-cells, J on (d−1)-faces
  (H(div)). This is the **diffusion module** operator. The classic
  lowest-order mixed-Poisson pair.

Both are inf-sup-stable, both are legitimate DEC discretizations,
neither is collocated — which is why (per 4b/4f) neither has a
checkerboard and both are stable *standalone*. The DSA one is **forced**
by the sweep's reduction (the sweep produces corner/edge scalar info
and cell-center currents — LMM's "cell-edge flux diffusion limit"),
which is why it, and not the RT0 one, is the consistent accelerator.
This makes the R4 "coexist by defining law" ruling a *theorem about two
different DEC arrangements*, not a tolerated duplication.

**Three structural warnings for the planned design (typed 0-skeleton
field for f₀ + the A_edge operator now; general machinery later):**

1. **The DSA f₀ EARNS a type — it passes the persistence test the
   retired `WavefrontFlux` failed.** `wavefront_cochain.rst` S6.4(f)
   retired the typed interior-face carrier with the doctrine: *"an
   octant transient cannot be a persistent mesh-bound field; persistent
   iterate state earns a type, sweep transients stay native."* The DSA
   f₀ is the **accumulated additive correction to φ across SI
   iterations** — mesh-bound, 0-skeleton-homed, and *persistent
   iterate-state* (it is the low-order solve's output, consumed by the
   next sweep). It is on the *earns-a-type* side of the doctrine, not
   the transient side. **Cross-check:** the angular ψ½ carrier
   (`radial_characteristic_space.py`, #282 route (a)) earned a type by
   the *same* argument — "promoted from a lagged solver-internal
   estimate to typed state so the sweep is a genuine triangular solve"
   — i.e. persistent back-edge state, not a transient. Two independent
   precedents agree: **persistence, not cochain-degree, is the
   type-earning criterion.** Type f₀.

2. **DEFER the general k-cochain-field machinery — there is exactly ONE
   *persistent* typed cochain, not two.** The tempting abstraction is
   "a generic `TypedCochainField[k]` on a named skeleton." Resist it now
   (L-004; coding-elegance defer-until-≥2). The apparent second instance
   — the interior face-flux (d−1)-cochain — does **not** count, because
   S6.4(f) demoted it to a *transient* (it stays native numpy in the
   per-octant walk). So the persistent-typed-cochain count is **1** (f₀
   once built), not 2. Build the concrete `f₀`-on-0-skeleton field +
   `A_edge` now (concrete architecture + imminent consumer — the
   build-now side of L-004 / feedback-defer-only-when-vague); mint the
   general machinery only when a **second persistent** cochain lands.
   The named trigger that flips the verdict: the **2-D DSA corner field
   (#314)** — a second persistent 0-cochain on a different mesh — OR a
   nodal-DG flux state. When either lands, the (1-D edge f₀, 2-D corner
   f₀) pair is two *persistent* instances of "0-cochain field on the
   0-skeleton," and the general machinery has its ≥2 justification.

3. **The 1-D degeneracy hides the design axis — do not over-commit the
   API to it.** In 1-D the 0-skeleton (vertices) *coincides* with the
   (d−1)-faces (also vertices): d−1 = 0. So a 1-D-only `A_edge` **cannot
   distinguish** "0-cochain home" from "(d−1)-cochain home" — the whole
   skeleton-home axis is *invisible* in 1-D (D2 Part A already notes the
   derived and landed operators live on "different spaces" but in 1-D
   both spaces are vertex-sets). The skeleton-home parameter becomes
   load-bearing only at 2-D, where Alcouffe's **corners (0-skeleton)**
   are genuinely distinct from **edges ((d−1)-faces)**. **Consequence:**
   type f₀ as a 0-cochain, but do NOT bake the 1-D vertex-set identity
   into the field's API (e.g. do not let f₀ silently share the face-
   trace's layout just because they collide in 1-D). This mirrors the
   R12a degeneracy discipline in `radial_characteristic_space.py`
   (τ_raw ∈ {0,1} are degenerate cases that hide the generic structure;
   the predicate is stated on the *generic* τ_raw ∈ (0,1)). State the
   0-cochain home on the *generic* (2-D) footing so #314 inherits it,
   not on the 1-D collision.

**Elegance payoff:**
- *Structure-exposing:* "which staggering?" becomes "which skeleton-home
  pair?", a declared property of the field, and the two coexisting
  operators become two named DEC arrangements (covolume vs RT0) rather
  than a suspicious duplication.
- *Structurally-simpler:* the skeleton-home axis unifies the spatial
  face-flux cochain, the DSA 0-cochain, and (via 4h) the angular
  half-angle-face cochain under one design vocabulary — but *without*
  minting the machinery until the count justifies it.

**First test (discriminates "type f₀ now, defer the machinery" from
"build the general machinery now"):** attempt to make a single generic
`TypedCochainField[k]` serve BOTH the DSA f₀ (k=0, persistent) and the
interior face flux (k=d−1). **Prediction:** it fails the S6.4(f)
constraint — the face flux has no *persistent mesh-bound* home in the
per-octant walk (it is a rolling frontier), so forcing it into a
persistent typed field either resurrects the exact carrier that was
already retired for cause, or the "generic" field needs a persistence
flag that the face-flux instance sets to False — at which point the
abstraction has one real (persistent) client and is premature. The
RED is the resurrection of a retired type. A "build the machinery now"
plan predicts two clean clients; the retired-transient is the
counterexample. (Positive control: the SAME test at 2-D-DSA-landing —
f₀-edge-1D and f₀-corner-2D are BOTH persistent, the machinery gets two
clean clients, the abstraction is justified. The test flips from RED to
GREEN exactly when the second persistent cochain arrives — that is the
named trigger.)

**Structural attack on current:** the plan's instinct (type f₀ + A_edge
now, general machinery at the second locus) is correct and the frame
*confirms* it — but the codebase does not yet state the persistence
criterion as the gate on cochain-typing. `wavefront_cochain.rst` states
it for the face flux (S6.4(f)); `radial_characteristic_space.py` states
it for ψ½ (#282); the DSA f₀ will be the *third* independent application
of the same rule. That recurrence is itself the signal that
"persistent-iterate-state ⇒ typed cochain; sweep-transient ⇒ native"
is a promotable project-level doctrine, not three coincidences.

### 4h. The angular axis is the same cochain complex — VALIDATED, and it unifies the consistency story across axes

**Trigger:** an induced mesh on the angular axis (half-angle faces from
cumulative weights), a conservation law on it (the redistribution
`(1−µ²)/r ∂_µ`), and no-flux boundaries (poles µ=±1 where (1−µ²)=0)
(structural feature 10; `radial_characteristic_space.py`).

**VALIDATED — the angular axis is a full cochain complex with its own
diffusion-limit-consistency story, structurally identical to the
spatial one.** The mapping is exact:

| Spatial axis | Angular axis (curvilinear SN) |
|---|---|
| cells (d-cells) | µ-bins (angular cells) |
| faces (d−1 cochain) | half-angle faces µ_{m±½} (the ψ½ carriers) |
| streaming `Ω·∇` = spatial conservation law | redistribution `(1−µ²)/r ∂_µ` = angular conservation law |
| vacuum/reflective spatial BC | no-flux poles µ=±1 ((1−µ²)=0) |
| non-uniform / staggered mesh | τ_raw ≠ ½ off-centered half-angle faces |
| LMM (1987): spatial diffusion-limit accuracy | Bailey–Morel–Chang (2010): **angular** diffusion-limit accuracy |

**The angular axis has its OWN "consistency = preserve the diffusion
limit" requirement, and the codebase already implements the consistent
choice.** Bailey–Morel–Chang (2010) abstract (p. 149): the diffusion
equation is derived from transport "through a **Galerkin approximation
based upon an angular trial space that is linear in the direction
cosines**" — i.e. the ℓ≤1 angular reduction, *the same Galerkin/Schur
reduction as spatial DSA* (4c), but on the µ-axis. Morel–Montry chose
the weighted-diamond **angular** weights precisely "to ensure that the
Galerkin diffusion approximation was preserved, which **eliminated the
discrete-ordinates flux dip**." BMC then prove the equivalence the
whole frame turns on: "**preserving the Galerkin diffusion
approximation is equivalent to preserving the asymptotic diffusion
limit to first order**," and step/diamond angular schemes preserve it
only to *leading* order (so their flux dip shrinks but does not vanish),
while Morel–Montry weighted diamond preserves it to *first* order (full
consistency). **The codebase's `τ_raw` Morel–Montry weights
(`radial_characteristic_space.py`) ARE this consistent angular scheme.**

The structural identity is therefore complete:
- **the angular "flux dip" is the angular analog of the spatial
  thick-cell oscillation** (4f) — a dispersion defect that vanishes in
  the diffusion limit for a diffusion-limit-preserving scheme;
- **consistency on the angular axis = the ℓ≤1 angular Galerkin
  reduction preserving the diffusion limit to first order** = the *same
  operation* as spatial DSA consistency (4c), on a different induced
  mesh;
- **the poles are no-flux angular boundaries**: the three pole-vanishing
  quantities in `radial_characteristic_space.py` (M1 the moment/output
  weight = 0; M2 the angular through-flux coefficient (1−µ²)|_{µ=±1} =
  α_{1/2} = 0; M3 …) are the angular-BC statement — the redistribution
  flux vanishes at the poles exactly as the spatial streaming flux
  vanishes at a vacuum face with no incident current.

**This is the deep unification the frame delivers:** the "k-skeleton
DoF home + diffusion-limit consistency" machinery is **axis-agnostic**.
Spatial (LMM 1987) and angular (BMC 2010) discretizations are the *same
cochain machinery on different induced meshes*, and BOTH carry the
identical "preserve the diffusion limit / Galerkin ℓ≤1 reduction"
consistency requirement. The half-angle faces are (d−1)-cochains on the
angular 1-mesh; the redistribution is the angular codifferential; the
poles are no-flux angular BCs. (This extends the already-validated
`psi_half_seed_angular_trace_frames.md` result — the ψ½ seed is the
*angular-inflow trace* of the (r,µ) rectangle, dual to the spatial
inflow trace — from "the seed is an angular trace" to "the whole
angular axis is a cochain complex with a consistency theorem.")

**Elegance payoff:**
- *Structure-exposing:* the redistribution term stops being "a
  curvilinear correction" and becomes "the angular conservation law of
  an induced-mesh cochain complex," with the poles as its no-flux BC and
  the Morel–Montry weights as its consistency closure. (Adjacent to
  L-001/L-008's Christoffel note: the `(1−µ²)/r ∂_µ` redistribution is
  where the differential-geometry connection-coefficient frame *does*
  fire; the cochain frame here is the *conservation-law* face of the
  same term.)
- *Expressive:* one sentence — "consistency = the ℓ≤1 Galerkin reduction
  preserves the diffusion limit to first order" — now covers spatial DSA
  (4c), angular differencing (this), and any future acceleration on
  either axis.
- *Structurally-simpler:* unifies the spatial and angular discretization
  frameworks under one cochain vocabulary, which is a *latent asset* for
  the deferred curvilinear DSA (`acceleration.rst` sn-dsa-honest-scope:
  "curvilinear — no stability theory exists"): the angular
  consistency theory (BMC) is the missing half of that stability theory.

**First test (discriminates angular-consistency from angular-scheme-
indifference):** on a spherical, optically-thick, scattering-dominated
(diffusive) problem, measure the **flux-dip magnitude** at the center
for (A) the Morel–Montry weighted-diamond angular scheme (the codebase's
τ_raw ∈ (0,1)) vs (B) plain diamond angular differencing (τ = ½
unweighted). **Prediction:** (A) preserves the angular diffusion limit
to first order → **no flux dip** as the diffusion limit is approached;
(B) preserves it only to leading order → a **flux dip that shrinks but
persists** at finite thickness (BMC's confirmed conjecture). The
discriminator is the *rate* at which the dip vanishes with optical
thickness: first-order-consistent kills it, leading-order-only leaves a
residual. A "the angular differencing does not affect the diffusion
limit" claim predicts equal dips — the differing dip decay is the RED.
This is the **angular analog of the D2 spatial-consistency spectral-
radius table** (D2 Part C), and it is the concrete first experiment for
any future curvilinear-DSA work (#229/#193/#9).

**Structural attack on current:** the curvilinear angular discretization
in the codebase is *already* the consistent (Morel–Montry) choice, but
the codebase does not yet NAME it as "the angular sibling of DSA
consistency." The R12a τ_raw predicate is stated as a
*layout/rank* fact (which levels carry a ψ½ block); the frame reveals it
is *also* the angular-consistency closure (the weights that preserve the
angular diffusion limit). That dual reading is a latent asset: when
curvilinear DSA is eventually attempted, the angular half of its
stability theory is already implemented and merely needs to be
recognized as such.

---

## 5. Where the parallel honestly BREAKS

The frame is strong, but four breaks are load-bearing — naming them is
what keeps it from being over-applied.

**Break 1 — the saddle point is the DIFFUSION REDUCTION, not the
transport operator. (The most important break.)** In CFD the saddle
point `[[A,Bᵀ],[B,C]]` **is** the thing you solve (the Stokes/NS
system is the top-level problem). In transport the saddle point is the
**ℓ≤1 moment reduction** used to *accelerate* a different primary
solve. The primary transport operator is either the characteristic
sweep `L⁻¹` (triangular, NOT a saddle point) or the Peierls integral
operator `I − PL⁻¹Σs` (a **compact perturbation of the identity**, NOT
a saddle point — Century p. 82, eq. 1.219). The mixed/saddle structure
appears **only** when you truncate to P1 — which is exactly why
`orpheus.diffusion` is the resolvent-backbone **exception** (L-007):
diffusion is elliptic-self-adjoint (⇒ mixed/saddle when first-ordered),
while SN/MoC/CP are characteristic-triangular (⇒ sweeps, never saddle
points). **The saddle-point frame fires on the low-order/diffusion
member and NOWHERE else in the transport stack.** The backbone and this
frame agree from one principle: elliptic ⇒ saddle-point-when-mixed;
characteristic ⇒ triangular-sweep. Do not attack the sweep with mixed-
FEM inf-sup theory — it has no saddle to stabilize.

**Break 2 — the failure-mode parallel breaks (no primal checkerboard).**
The correspondence is of saddle-point **roles** (φ↔p, J↔u), not of
**failure modes**. CFD's signature instability — the collocated primal
checkerboard (a φ/p null mode) — has **no transport counterpart** (4f):
transport is born staggered. Transport's signature instability — Reed's
inter-operator consistency divergence — is the one CFD calls
"inconsistent projection," and it is NOT CFD's *famous* failure. So the
two domains share the saddle-point skeleton and the *consistency*
failure, but their *canonical* instabilities are different objects. A
naive transfer that expects a checkerboard in transport (and reaches
for pressure stabilization) is importing a fix for a disease transport
does not have.

**Break 3 — Darcy, not Stokes: the current elimination is closed-form.**
Because the (J,J) block is a **mass matrix** (not a Laplacian), the ℓ=1
current is eliminated **exactly and locally** (Larsen's four-step,
step 4) — the transport Schur complement IS a sparse diffusion
operator. In incompressible NS the (u,u) block is the vector Laplacian,
`A⁻¹` is nonlocal, and the pressure Schur complement is spectrally a
**mass/identity** operator, not a diffusion operator — there is no
closed-form velocity elimination. So the parallel to *incompressible
NS specifically* breaks at the elimination step; the parallel to
**Darcy / mixed-Poisson** holds exactly. Import the H(div)/mixed-Poisson
branch, never the Stokes velocity-pressure branch (4a).

**Break 4 — the angular axis has no CFD sibling.** CFD is a spatial-only
saddle point. Transport's phase space carries an **angular axis** with
its own cochain complex and its own diffusion-limit-consistency story
(4h). The DEC/mixed machinery *does* extend to the angular axis (that
is the 4h unification), but CFD's specific pressure–velocity vocabulary
does not — there is no "angular continuity equation" in fluids. So the
frame's angular half is a transport-only structure; only its *cochain*
generalization (not its *CFD* instantiation) reaches it.

**Break 5 (minor) — the constraint strength is tunable (1−c), not
fixed.** The φ-block `Σa = Σt(1−c)` interpolates from a coercive
reaction (c<1) to a pure constraint (c=1). CFD incompressibility is
always the hard-constraint end. So transport's inf-sup requirement
*tightens continuously as c→1*, whereas CFD's is uniform. The
compatibility problem is c-parametrized in transport — which is exactly
why the D2 stability tables sweep c and why the worst case is the
reflective pure-scatterer (c=1, `Σa=0`, feature 2 caveat).

## 6. Pollination candidates (CFD / mixed-FEM machinery → ORPHEUS)

Enumerated with trigger, reformulation, and a discriminating first test.
Cited targets: #312 (LD arm), #314 (2-D DSA), #294 (diffusion
follow-ups), #200 (Krylov preconditioner).

**P1 — Discrete inf-sup (LBB) test as a fast *necessary* screen for
future (φ,J) pairs [→ #314].** *Trigger:* a new low-order mixed pair
proposed on a new mesh/skeleton. *Borrowing:* compute the discrete
inf-sup constant β (smallest generalized singular value of B relative
to the A- and C-norms — the standard mixed-FEM eigenproblem) to reject
incompatible pairs cheaply, BEFORE the full spectral-radius sweep.
*Caveat (from 4a/4c — this is the load-bearing refinement):* inf-sup
stability is **necessary but NOT sufficient** for a consistent
accelerator — a Taylor–Hood-stable pair is inf-sup-stable and still
diverges as an accelerator because it is not the sweep's reduction. So
import the inf-sup test as a *fast negative screen* (β→0 ⇒ reject), not
as a green light. *First test:* compute β for (Alcouffe corner-φ,
cell-J) vs a deliberately-collocated (cell-φ, cell-J) 2-D pair;
predict β bounded-below for the covolume pair, β→0 under refinement for
the collocated one — the collapsing β is the RED that flags the
checkerboard-prone pair before any transport run.

**P2 — NDA (nonlinear diffusion acceleration) as the nonlinear sibling
[→ #294, or a new NDA issue]. (The strongest pollination candidate.)**
*Trigger:* a regime where *linear* consistency is hard or impossible —
voids (`D=1/(3Σt)→∞` as Σt→0), strong material discontinuities, or
unstructured meshes (the [124] regime). *Borrowing:* Hammer–Morel–Wang
(2019) NDA updates the low-order operator's closure (a
consistency/`\hat{D}` term) **from the high-order solution each
iteration**, so the low-order reproduces the high-order balance
**exactly at convergence, regardless of the low-order staggering**. It
is the transport analog of the CFD **SIMPLE family** (a nonlinear
pressure-correction) and of deferred-correction. It is the answer to
the [124] multi-D discontinuity degradation that does *not* require a
Krylov wrap: enforce consistency **nonlinearly** instead of
structurally. It also natively handles **voids**, where linear
diffusion DSA fails — a real reactor-physics need (gaps, coolant
channels, rod followers). *First test (discriminates NDA from linear
DSA):* build a *deliberately wrong-staggered* low-order operator and use
it (a) as a linear DSA corrector and (b) as an NDA low-order.
*Prediction:* (a) diverges (D2 Part C, the wrong staggering); (b)
**still converges and accelerates** because the nonlinear consistency
term corrects the balance regardless of staggering — and the converged
NDA low-order flux equals the high-order flux *exactly* (a discrete-
balance identity linear DSA lacks). The RED for a "NDA is just DSA with
extra steps" claim is that (b) converges where (a) diverges.

**P3 — Rhie–Chow / Brezzi–Pitkäranta stabilization *scaling* recipe for
the M4S LD arm [→ #312].** *Trigger:* #312's LD-consistent DSA is
"structurally unspellable without the M4S reduction" (partially-
consistent by necessity). *Borrowing:* the CFD collocated-stabilization
literature is the design oracle for *how much* stabilization and *how
to scale it* — the Rhie–Chow rule is "scale the compact correction by
the momentum-operator diagonal `1/a_P`," the exact-scaling requirement
(4d). *Caveat (from 4f):* import it as an **exactly-scaled compact
coupling correction** (the Rhie–Chow scaling), NOT as a φ-null-mode
stabilizer (transport has no null mode). *First test:* the McCoy–Larsen
Table II dose-threshold migration (4d's negative control) — an
incorrectly-scaled M4S closure diverges at a σₜh threshold that migrates
with the mis-scaling; the correctly-scaled one is flat. A constant
(un-scaled) stabilizer predicts a uniform bounded ρ — the migrating
threshold is the RED.

**P4 — MINRES for the symmetric-indefinite LD-on-tets operator, and
block-triangular saddle-point preconditioning [→ #200, #312].**
*Trigger:* Century p. 83 — the LD-DFEM DSA-preconditioned operator is
**symmetric-indefinite** (SAPD not preserved on tets), precluding CG.
*Borrowing:* the saddle-point Krylov literature (Elman–Silvester–Wathen)
says the natural method for a symmetric-**indefinite** system is
**MINRES** (not GMRES, not CG), and the natural preconditioner is
block-diagonal/block-triangular with the Schur complement approximated
by the diffusion operator. The codebase's `M=(I+C)(L+C)⁻¹` is already a
block-triangular preconditioner; the import is the *method choice*
(MINRES) for the indefinite case and the block-preconditioner *theory*
for choosing the Schur approximation. *First test:* on the LD-on-tets
case (when #312 lands), MINRES with block-diagonal DSA preconditioning
converges where CG stalls/breaks (CG requires SPD; the operator is not)
— the CG breakdown is the RED that forces the MINRES import.

**P5 — the covolume / mimetic-finite-difference branch (not RT0) for
the 2-D corner operator [→ #314].** *Trigger:* the DSA arrangement is
**covolume** (0-skeleton φ, dual-0 J — 4g), not RT0. *Borrowing:* the
mimetic-finite-difference / covolume literature (Hyman–Shashkov;
Lipnikov–Manzini–Shashkov) specializes in exactly the primal-dual
(0-skeleton, dual-0) arrangement and supplies the discrete
Green's-identity / div-grad-adjointness the RT0 machinery does not
directly give for the covolume pair. This is a *refinement* of "import
mixed-FEM": import the covolume branch specifically. *First test:* the
2-D corner operator must satisfy the mimetic discrete Green's identity
(discrete `⟨div J, φ⟩ = −⟨J, grad φ⟩` on the primal-dual complex, up to
the boundary trace) — a non-mimetic 2-D corner stencil breaks it, and
that broken identity is the RED (it is the 2-D generalization of the
R/P adjoint-consistency the campaign already pins in 1-D via the
W-self-adjoint `angular_frame(0)`).

**P6 — the consistent-projection (adjoint div/grad) discipline for 2-D
DSA [→ #314].** *Trigger:* 2-D low-order assembly. *Borrowing:* the
fractional-step/SIMPLE lesson (4b) — the discrete `div` and `grad` in
the low-order must be the **adjoint pair implied by the high-order**,
or the projection does not project. *First test:* verify R = P* (up to
the trace measure) on the 2-D corner/cell pair; a non-adjoint pair
breaks particle conservation (the `4π/2` normalization bug the campaign
already flags in 1-D) — the conservation defect is the RED.

## 7. Elegance assessment — the k-skeleton design question

**Decision: build the concrete typed field now; defer the general
machinery to the second persistent instance. This IS the elegant path
— building the machinery now would score *lower*, not higher.**

| What | Verdict | Four-criteria score (if built at the right time) |
|---|---|---|
| Typed 0-cochain `f₀` field + `A_edge` operator | **BUILD NOW** | structure-exposing (names the 0-skeleton home) + structurally-simpler (retires the ad-hoc edge array) — concrete, imminent consumer, passes the persistence test |
| General `TypedCochainField[k]`-on-a-skeleton machinery | **DEFER** to the 2nd persistent cochain (#314 corner-f₀ or nodal-DG) | structure-exposing + expressive (unifies spatial/angular/DSA) + structurally-simpler (2 DEC arrangements named); NO algorithmic-advantage (organizational win) → 3/4, but only *after* the count justifies it |

The elegance reasoning, stated as a decision rather than a taste:

- **`f₀` earns a type by the persistence criterion** (4g warning 1), the
  same criterion that earned ψ½ a type (#282) and *retired* the
  interior-face carrier (S6.4(f)). Three independent applications of one
  rule ⇒ the rule is real: **persistent-iterate-state ⇒ typed cochain;
  sweep-transient ⇒ native**. This is a promotable project doctrine
  (candidate for `coding-elegance` / a rule), not three coincidences.
- **The general machinery has exactly ONE persistent client today** (f₀),
  because the face-flux cochain is a *transient* and does not count
  (4g warning 2). Building the abstraction now is the premature-
  abstraction anti-pattern, and it would *fail its own first test*
  (4g: forcing the transient face-flux into a persistent typed field
  resurrects a retired type). **The elegant move is to WAIT** — and the
  waiting has a *named trigger* (the 2-D corner field #314), not an
  open-ended "someday."
- **The 1-D degeneracy is an elegance trap** (4g warning 3): 1-D
  collapses the 0-skeleton onto the (d−1)-faces, so a 1-D-tuned API
  would encode a *coincidence* as structure. State the 0-cochain home on
  the generic 2-D footing.

Net: the k-skeleton frame is a **high-value organizational unification**
(it retires the "which staggering?" degree of freedom and unifies three
cochain instances), but its algorithmic payoff is zero and its
abstraction is premature by one instance. The elegant realization is
therefore *staged*: concrete field now, machinery at #314. This is
precisely the `feedback-defer-only-when-architecture-vague` +
`defer-until-≥2` discipline, with the persistence criterion as the
sharpened counting rule.

## 8. Acquisition candidates

**Hard requirement (the frame depends on it and it is NOT local):**
- **Warsa, Wareing & Morel, "Krylov iterative methods and the degraded
  effectiveness of DSA for multidimensional SN calculations in problems
  with material discontinuities," NSE 147:218 (2004)** — the [124]
  primary behind the entire 4e mode-(ii) argument. The two local
  secondaries (the Century chapter's [124] summary + the WWM *147:26*
  k-eigenvalue paper) carry the *claims* but not the primary's data,
  proofs, or the precise conditions under which consistent DSA degrades.
  Needed before #314 relies on the mode-(ii) prediction quantitatively.
  **Do NOT acquire online (per brief); flag for the user to add to
  `scratch/literature/`.**

**Soft (acquire only if the named pollination is pursued):**
- **Boffi, Brezzi & Fortin, *Mixed Finite Element Methods and
  Applications* (2013)** — the canonical inf-sup/LBB reference; needed
  for P1 (the inf-sup screen) and to state the Darcy-vs-Stokes stable-
  pairing catalogs precisely. A math reference, not transport.
- **Elman, Silvester & Wathen, *Finite Elements and Fast Iterative
  Solvers* (2014)** — the canonical saddle-point block-preconditioner +
  MINRES reference; for P4 (#200/#312).
- **Rhie & Chow (1983), AIAA J. 21:1525** (or Ferziger–Perić,
  *Computational Methods for Fluid Dynamics*, the collocated-FV chapter)
  — the momentum-consistent interpolation scaling recipe; for P3 (#312).
- **Lipnikov, Manzini & Shashkov, "Mimetic finite difference method,"
  JCP 257 (2014)** — the covolume/MFD review; for P5 (#314's corner
  operator).

**Already local (pointers, not acquisitions):** Adams–Larsen (2002)
§IV.D.4 (the 9-point multi-D consistent stencil, for #314); Alcouffe
(1977) §II.C (the 2-D corner-centered derivation, Eqs. 41–43);
Hammer–Morel–Wang (2019) (NDA, for P2); Bailey–Morel–Chang (2010) (the
angular consistency theory, for 4h and curvilinear DSA).

## 9. UNEXPLORED — frames checked, no trigger fired

Each carries the one-line STRUCTURAL reason it did not fire (so the next
session does not re-attack a dead frame).

- **Homology / chain complex** — the DEC *cochain* structure (k-cochains,
  the codifferential δ) is load-bearing (4g/4h), but *homology* is not:
  no `∂²=0` payoff, no closed-vs-exact distinction fires. The boundary
  trace + its extension are a **dagger adjoint pair, not a differential**
  (L-001). Use the cochain-degree/skeleton-home axis; drop the
  cycles/boundaries machinery.
- **Category theory / double category** — the coupling/block-re-
  association win is already captured concretely by the #208 **biproduct**
  algebra (`Mat∘Mat≅Mat`, L-011). No abstract-nonsense lever beyond the
  biproduct is needed; a functor/2-cell restatement adds ceremony, not a
  test.
- **Tensor networks / MPO** — no bond-dimension trigger. The saddle point
  is a fixed 2×2 (3×3 with the trace block) biproduct — bond-dimension
  degenerate, not a rank-N chain with a truncation knob (L-001).
- **Multigrid** — NOT a foreign frame here: it is the NATIVE one.
  Consistent DSA IS the Galerkin coarse-grid operator `R·A_high·P` +
  Schur complement (L-009 corollary; Century p. 81 says DSA-type
  preconditioners "are analogous to the coarse-grid operators used in
  multigrid"). Already internal; nothing to import.
- **von Neumann / Fourier stability** — the native analysis tool the
  campaign already uses (the ρ<0.2247c bound). Not a cross-domain
  import; listed to mark it in-frame.
- **Optimal control / Riccati** — the saddle point is a static KKT
  (min-max) system; it can be read as the stationarity of a
  complementary-energy functional, but there is NO temporal/dynamic
  recursion, so no Riccati lever. The Schur-complement reading (4c)
  already extracts everything the KKT reading would.
- **Symplectic / Hamiltonian** — the `[[A,Bᵀ],[B,−C]]` sign pattern
  superficially resembles a symplectic block form, but there is no
  symplectic 2-form and no phase-space flow being preserved; the "saddle"
  is a stationarity, not a Hamiltonian flow. No conserved symplectic
  structure fires.
- **Wiener–Hopf / H-function** — wrong solver family (L-001): native to
  the half-space Chandrasekhar line, structurally incompatible with the
  sweep + diffusion-reduction formulation. The DSA low-order is a
  moment reduction, not a half-space factorization.
- **Domain decomposition / Schwarz** — adjacent but untriggered by *this*
  hypothesis set: it fires on "how to parallelize the sweep across
  subdomains / couple subdomain interfaces," which is orthogonal to the
  consistency/staggering question. Would become relevant for a parallel-
  sweep or multi-block 2-D DSA (#314) *implementation*, not for its
  *formulation*.
- **Differential geometry / Christoffel** — fires ONLY on the curvilinear
  angular redistribution `(1−µ²)/r ∂_µ` (L-001/L-008), which appears in
  4h — but there its *conservation-law* (cochain) face is the load-
  bearing one for the consistency question; the connection-coefficient
  face is a *separate* frame for the streaming-redistribution term, not
  for DSA consistency. Correctly scoped to 4h's redistribution note,
  not to the saddle-point frame.
