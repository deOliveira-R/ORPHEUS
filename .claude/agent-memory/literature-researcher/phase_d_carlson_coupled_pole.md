---
name: Phase D Carlson coupled-pole sweep — Hébert §3.9.4 Eqs. (3.432)-(3.435)
description: Method-implementer memo for the inward μ=−1 starting-direction sweep that seeds the sphere SN M-M α-cascade. Closed-form per Hébert; ORPHEUS Phase D implements this as a PoleFaceInitialCondition strategy to close ERR-026 on sphere Gate 1.1 MMS.
type: reference
---

# Phase D — Carlson coupled-pole sweep at the spherical SN centre

**Cite when**: Phase D implementation work (PoleFaceInitialCondition
Protocol + CarlsonCoupledPole strategy), ERR-026 sphere Gate 1.1
MMS failure, any architectural decision about how the M-M angular
closure obtains its pole-face seed at r=0 on sphere.

**V&V pillar classification**: this memo formalises a
**closed-form analytical** reference for the spatial pole-face seed.
The Carlson coupled-pole sweep is a STRUCTURAL feature of the
Morel–Montry angular closure on sphere — it is NOT a verification
truth-set itself, but the math it implements is closed-form per
Hébert (2009) and reduces to a deterministic recurrence on a fixed
spatial mesh.

**Primary source — LOCAL**: `scratch/literature/Hebert(2009)Chapter3.pdf`,
§3.9.4 pp. 141-144. All equation numbers in this memo are Hébert's.

**Companion memo**: `sphere_sn_pole_closure_canonical.md` (the
Eqs. (3.418)-(3.420), α-recursion (3.423)-(3.424), and Bailey-Morel-Chang
2010 NSE 165 justification background). This memo extends that one with
the Phase D-specific depth on (3.432)-(3.435).

---

## 1. Verbatim equations

Hébert §3.9.4 opens the sphere difference relations at Eq. (3.418)
(angularly-integrated divergence form). Then introduces the α-recursion
(3.423)-(3.424) and the cell-balance with redistribution divisor
ΔS_i / (2·𝒲_n) at Eq. (3.428).

The Carlson coupled-pole sweep begins at Eq. (3.432). Direct quote
(top of p. 143):

> "We first need knowledge of mesh-centered values of the flux φ_{1/2,i}
> over the zero-weight points, in order to initialize the recursive
> solution procedure. We observe that these directions correspond to
> particles entering the external surface and moving toward the central
> axis with μ = −1. The angular redistribution term vanishes on these
> points so that Eq. (3.164) simplifies to"

### Eq. (3.432) — continuous form at μ = −1

```
−∂/∂r φ_{−1/2}(r) + Σ(r) φ_{−1/2}(r)
     = Σ_{ℓ=0..L} (2ℓ+1)/2 · Q_{ℓ}(r) · P_{ℓ}(−1)
```

**Derivation context**: this is the streaming–collision form of the
1D-sphere Boltzmann equation specialised to μ = −1. The full sphere
streaming operator (Hébert Eq. (3.418)) reads

```
(1/r²) μ ∂/∂r [r² φ] + (1/r) ∂/∂μ [(1−μ²) φ]
```

The angular-redistribution piece (1−μ²)φ is the term Hébert calls
"vanishing at μ = ±1". Setting μ = −1 collapses the streaming piece
to `−∂φ/∂r` (the curvature 1/r² factor cancels under the integration
over the angular sub-domain that produced Eq. (3.418)). The right-hand
side is the standard Legendre expansion of the scattering source `Q`
evaluated at μ = −1, where `P_ℓ(−1) = (−1)^ℓ`.

**Why P_ℓ(−1) matters**: the source at μ = −1 retains scattering
moments of all orders ℓ, but with alternating signs. Isotropic
scattering (ℓ=0 only) gives a flat source `Q_0(r)/2` along this
ordinate; linearly anisotropic adds `−(3/2) Q_1(r)`.

**Subscript convention**: Hébert uses `−1/2` as the half-integer
index for the auxiliary inward direction; the α-cascade then starts
from `α_{1/2} = 1 − (−1)² = 0` (Eq. (3.423)). The "−1/2" subscript
is purely an indexing choice — it labels the inward starting ordinate
that lives one half-step ABOVE μ = −1 in the cascade, not a physical
ordinate at μ = −0.5.

### Eq. (3.433) — DD-discretised cell balance at μ = −1

```
−(φ̄_{i+1/2} − φ̄_{i−1/2}) + Δr_i · Σ_i · φ̄_i = Δr_i · Q̄_i
```

where (Hébert's typographic definition, immediately after the equation):

```
φ̄_i ≡ φ_{1/2, i}
Q̄_i ≡ Q_{1/2, i}
Δr_i = r_{i+1/2} − r_{i−1/2}
```

**Derivation context**: Eq. (3.432) is integrated over sub-mesh `i`
(volume between `r_{i−1/2}` and `r_{i+1/2}`). The streaming term
yields the discrete jump `φ̄_{i+1/2} − φ̄_{i−1/2}` (no r² weight,
because Eq. (3.432) is the cleanly-divergence form along μ = −1
where the geometric weight already cancelled). The negative sign on
the left comes from μ = −1 < 0 — particles travel INWARD, so the
"outgoing − incoming" jump is `i−1/2` minus `i+1/2`, equivalently
`−(φ_{i+1/2} − φ_{i−1/2})`.

**Key simplicity**: NO α-redistribution divisor appears in Eq. (3.433)
because angular redistribution vanishes at μ = ±1 (the (1−μ²) factor
is zero). This is the entire reason Hébert can solve the μ = −1 sweep
in closed form with a plain DD recurrence: the coupled angular cascade
is decoupled at the starting direction.

### Eq. (3.434) — closed-form average flux from outgoing-face value

```
φ̄_i = (Δr_i · Q̄_i + 2 φ̄_{i+1/2}) / (Δr_i · Σ_i + 2)
```

**Derivation context**: combine the DD relation Eq. (3.431)
(`φ_{n,i} = ½(φ_{n,i−1/2} + φ_{n,i+1/2})`, specialised to the
−1/2 starting ordinate) with the balance Eq. (3.433), then solve
algebraically for `φ̄_i` in terms of the known outgoing-face value
`φ̄_{i+1/2}` (the face FURTHER from the centre — known because we
sweep inward from the outer BC) and the source `Q̄_i`. The `2` in
the numerator and denominator comes from the DD identity multiplied
by 2 (eliminating the unknown incoming face `φ̄_{i−1/2}`).

**Sign sanity check**: μ = −1 < 0, so `μ_n < 0` branch of Eq. (3.386)
slab analog applies; Hébert's slab Eq. (3.386) for `μ_n < 0` reads
`φ_{n,i} = (Δx_i Q_{n,i} − 2μ_n φ_{n,i+1/2}) / (Δx_i Σ_i − 2μ_n)`.
Setting `μ_n = −1` gives `(Δx_i Q + 2 φ_{i+1/2}) / (Δx_i Σ + 2)`,
which is exactly Eq. (3.434) (the sphere has the same form because
the streaming reduces to slab-form at μ = ±1). Sign verified.

### Eq. (3.435) — DD reconstruction of incoming face

```
φ̄_{i−1/2} = 2 φ̄_i − φ̄_{i+1/2}
```

**Derivation context**: textbook diamond-difference (DD) auxiliary
relation Eq. (3.431) rearranged for the incoming face. This is the
recurrence that propagates the sweep: given the outgoing face
`φ̄_{i+1/2}` (known from previous step or from outer BC at the
first step), use Eq. (3.434) to get the cell-centred value `φ̄_i`,
then Eq. (3.435) to step inward to the next face `φ̄_{i−1/2}`.

**Closure note**: the pair (3.434)+(3.435) IS the spatial recurrence.
Together they realise a tridiagonal-style inward sweep on the radial
mesh: outer face → cell centre → inner face → next cell centre → ...
→ pole face `φ̄_{1/2}` at r = 0.

---

## 2. Inward zero-weight starting direction — what it MEANS

### Why a starting direction is needed

The M-M angular closure on sphere is a per-cell α-cascade
(Eqs. (3.436), (3.438) — i.e., the asymmetric per-ordinate sweep
with DD angular closure). The cascade COUPLES the angular flux
across ordinates within one spatial cell: ordinate `n` reads
α_{n−1/2} from the previous (more-inward-μ) ordinate. To start the
cascade at the smallest-μ ordinate, you need a value for the
α_{1/2} coefficient AND for the angular-edge flux φ_{1/2,i} at that
seed half-integer.

The α_{1/2} = 0 seed is FREE: it comes from `1 − μ²` evaluated at
μ = −1 (Hébert's text below Eq. (3.422): "The first value α is
equal to 1 − (−1)² = 0"). That handles the **angular** half of the
problem.

The flux value φ_{1/2,i}, however, is NOT free. It is the SPATIAL
flux profile at the auxiliary starting direction, and it must be
solved for as a function of position `i` along the radial mesh.
That spatial solve is exactly what Eqs. (3.432)-(3.435) provide.

### Why μ = −1 (most-inward direction) is natural

At μ = −1 the streaming operator collapses to pure radial divergence
WITHOUT the angular-redistribution coupling. From Hébert's text:
*"The angular redistribution term vanishes on these points so that
Eq. (3.164) simplifies to [Eq. (3.432)]."* Pomraning (1989) §V
explains the structural reason (see section 6 below).

This is the ONLY direction on the unit interval [−1, 1] where the
spatial 1D-sphere problem reduces to a closed-form linear recurrence
in radius alone, without an inner angular solve. Picking any
intermediate μ would leave the coupling term active and re-introduce
the cascade chicken-and-egg.

### Why zero-weight

In a Gauss-Legendre quadrature on [−1, 1] the endpoints μ = ±1 are
NOT base points (the polynomial is approximated by interior nodes
only). They have no quadrature weight, hence "zero-weight" — the
flux value at μ = −1 does NOT contribute to any Σ_n 𝒲_n φ_n
integral that builds the scalar flux moments.

The μ = −1 ordinate is therefore a **purely auxiliary** numerical
construct: its flux values exist for the sole purpose of seeding
the α-cascade for the finite-weight ordinates that follow. After
the cascade is initialised, the φ_{1/2,i} values are **discarded**
(Hébert's words: *"At completion of this initialization sweep, the
angular edge values are discarded and the mesh-centered values
φ_{1/2,i} are kept."* — p. 143, between Eqs. (3.435) and (3.436)).
Actually only the cell-centred values are kept, because those are
the ones that feed the finite-weight ordinates' Eq. (3.436)
denominator via the `(α_{n−1/2} + α_{n+1/2}) φ_{n−1/2,i} / (2 𝒲_n)`
redistribution term.

---

## 3. Spatial DD sweep at μ = −1 — propagation from outer BC to r = 0

### The recurrence on a uniform mesh

For an `nx`-cell uniform radial mesh with `Δr_i ≡ Δr` and per-cell
constants `Σ_i = Σ`, `Q̄_i = Q`, the inward sweep is:

```
Step 0 (outer BC at i = nx):   φ̄_{nx+1/2} = ψ_R^{(−1)}     # see BC discussion
Step k = nx, nx−1, ..., 1:
    φ̄_k       = (Δr · Q̄_k + 2 · φ̄_{k+1/2}) / (Δr · Σ_k + 2)     # Eq. (3.434)
    φ̄_{k−1/2} = 2 · φ̄_k − φ̄_{k+1/2}                              # Eq. (3.435)
Return: φ̄_{1/2}                                                  # pole-face seed
```

### Flat-ψ trace — verifying the hypothesis

The Phase D hypothesis is: *for a flat angular flux ψ_cell = C
across all cells, the μ = −1 inward sweep returns
`φ̄_{1/2} = C`.* Let us verify by direct algebra.

Setup: Σ_i = Σ, source Q̄_i constructed so that the consistent
fixed point is φ̄_i = C everywhere. The continuous Eq. (3.432) with
∂_r φ ≡ 0 (flat in r) demands `Σ · C = Σ_ℓ ((2ℓ+1)/2) · Q_ℓ · P_ℓ(−1)`,
i.e. `Q̄ = Σ · C` (lumping the right-hand side into the discrete
cell source `Q̄`).

Substitute into Eq. (3.434):

```
φ̄_i = (Δr · Σ · C + 2 · φ̄_{i+1/2}) / (Δr · Σ + 2)
```

If `φ̄_{i+1/2} = C` (inductive hypothesis at the outer face):

```
φ̄_i = (Δr · Σ · C + 2C) / (Δr · Σ + 2)
     = C · (Δr · Σ + 2) / (Δr · Σ + 2)
     = C
```

Then Eq. (3.435) gives `φ̄_{i−1/2} = 2C − C = C`. The recurrence is
self-similar: every face and cell value stays at `C`. Hence
**`φ̄_{1/2}(r = 0) = C` for flat ψ on flat source — the hypothesis
holds**.

This is consistent with the M-M cascade's fixed-point structure:
when the cell-centred fluxes are at their consistent fixed point,
the auxiliary ordinate's spatial sweep returns the same value at
every position, including the pole face. The α-cascade then propagates
this `C` value across every ordinate's pole-face IC without distortion.

### Vacuum-BC trace

For a vacuum outer BC, `φ̄_{nx+1/2} = 0` (no inward flux at the outer
surface). Then for a flat positive source Q ≠ 0 across an absorbing
medium (Σ > 0), the inward sweep produces a non-zero pole-face value
that depends on the cell count and (Δr · Σ) product:

```
φ̄_{nx}      = Δr · Q / (Δr · Σ + 2)                # from φ̄_{nx+1/2} = 0
φ̄_{nx−1/2}  = 2 · Δr · Q / (Δr · Σ + 2)
φ̄_{nx−1}    = (Δr · Q + 2 · φ̄_{nx−1/2}) / (Δr · Σ + 2)
            = Δr · Q · (1 + 4/(Δr·Σ+2)) / (Δr · Σ + 2)
... etc.
```

For a purely scattering medium (Σ_t = Σ_s, c = 1, no fixed source),
both Q and the consistent fixed point are zero everywhere, so the
inward sweep returns φ̄_{1/2} = 0. For a critical c = 1 problem
with non-zero fixed point (driven by reflective BC or by the
within-group fission integral), the algebra is similar but the
fixed point depends on the source moments.

**Take-away**: the inward sweep is geometry-aware (the outer BC
matters) but not particularly clever — it's a textbook tridiagonal
relax with a closed-form forward elimination. Ten lines of Python.

---

## 4. How the result feeds the outward α-cascade

### Where φ̄_{1/2,i} enters

Look at Hébert Eq. (3.436) (the per-ordinate cell balance for
μ_n < 0, immediately after the initialization Eqs. (3.432)-(3.435)):

```
φ_{n,i} = ( V_i Q_{n,i}
           − μ_n (S_{i−1/2} + S_{i+1/2}) φ_{n,i+1/2}
           + ΔS_i (α_{n−1/2} + α_{n+1/2}) φ_{n−1/2,i} / (2 𝒲_n)
          )
         / (
           V_i Σ_i
           − 2 μ_n S_{i−1/2}
           + ΔS_i α_{n+1/2} / 𝒲_n
          )
```

The CRUCIAL term is `(α_{n−1/2} + α_{n+1/2}) · φ_{n−1/2,i} / (2 𝒲_n)`
in the numerator. This is the angular-redistribution contribution
from the previous ordinate (n−1/2) at the SAME spatial cell `i`.

For the very FIRST ordinate after the auxiliary one — i.e., `n = 1`
which is the smallest-|μ| inward ordinate, μ_1 < 0 closest to 0 —
the value `φ_{n−1/2, i} = φ_{1/2, i}` IS the result of the inward
sweep Eqs. (3.432)-(3.435). Hébert's text after Eq. (3.435):

> "At completion of this initialization sweep, the angular edge
> values are discarded and the mesh-centered values φ_{1/2,i} are
> kept."

The kept `φ_{1/2,i}` values then appear in EVERY ordinate's cell
balance via the α-cascade: the (n−1/2) flux at ordinate `n` is
either (a) the previous ordinate's cell-centred flux φ_{n−1/2,i}
that the cascade just computed, or (b) the auxiliary seed φ_{1/2,i}
for the very first ordinate.

### Is the seed used "as-is" or transformed per-ordinate?

The seed φ_{1/2,i} is used **as-is** as the angular-edge value
`φ_{n−1/2, i}` for `n = 1`. There is NO per-ordinate transformation
of the seed. The transformation happens IMPLICITLY through the
α-coefficients `α_{1/2}, α_{3/2}, α_{5/2}, ...` which are themselves
ordinate-dependent and accumulate via the recurrence
`α_{n+1/2} = α_{n−1/2} − 2 𝒲_n μ_n` (Hébert Eq. (3.424)).

So the seed enters once, at `n = 1`; the α-cascade then chains
through subsequent ordinates, propagating the influence of the seed
via the running `(α_{n−1/2} + α_{n+1/2}) φ_{n−1/2,i} / (2 𝒲_n)`
term. The seed is the **initial condition for the angular cascade**
in the same sense that an ODE needs one initial value per step.

### Pole-face shape vs cell-centre values

NOTE: Hébert KEEPS the cell-centred values φ_{1/2,i} (not the
angular-edge values). The angular-edge values
`φ̄_{i±1/2}` are discarded at the end of the inward sweep. The Phase D
brief described the pole-face seed as `ψ_face(r=0)`; in Hébert's
notation this would be `φ̄_{1/2}` (the angular-edge value at the
innermost radial face). The IMPLEMENTATION needs to be careful:

- If the outward α-cascade reads cell-centred `φ_{n−1/2, i}` values
  (one per radial cell), it needs **all `nx` cell-centred values
  from the inward sweep**, not just the value at r = 0.
- If the outward α-cascade only needs a single pole-face flux for
  starting the cascade at the innermost radial cell, the
  angular-edge `φ̄_{1/2}` value at r = 0 suffices.

Reading Eq. (3.436) carefully: the term is `φ_{n−1/2, i}` — i.e.,
the angular-edge flux at angular index `n−1/2` and spatial cell
index `i`. For `n = 1` this is `φ_{1/2, i}`, i.e. the cell-centred
flux at the auxiliary ordinate, evaluated at radial cell `i`. So
the cascade needs the FULL spatial profile `{φ̄_i}_{i=1..nx}` from
the inward sweep, NOT just the value at r = 0. **The "pole-face
seed" mental model is too narrow** — it's really a "starting-direction
spatial-profile seed".

This has Phase D architecture implications (section 7).

---

## 5. Auxiliary ordinate — quadrature set member or constructed on-the-fly?

### Is μ = −1 in standard Gauss-Legendre sets?

**No**. Standard `N`-point Gauss-Legendre quadrature on [−1, 1]
places its base points at the roots of the Legendre polynomial
`P_N(μ)`. None of these roots equals ±1 — the endpoints are NOT
quadrature nodes. Hébert is explicit on this: "the base points μ_n
and weights 𝒲_n are those of an N-point Gauss-Legendre quadrature"
(text above Eq. (3.419)).

Some quadrature variants (Gauss-Lobatto, Gauss-Radau) DO include
endpoints, but Hébert's choice and the standard ORPHEUS
`sn_quadrature` module use plain Gauss-Legendre. μ = −1 is purely
auxiliary in this context.

### Quadrature-set extension vs on-the-fly construction

The architectural decision: **the μ = −1 sweep should be
constructed on-the-fly inside the angular-closure matvec, NOT added
as a virtual ordinate to the quadrature set.**

Reasoning:

1. **Zero-weight pollution**. Adding a weight-zero ordinate to the
   quadrature would force every downstream consumer
   (`scalar_flux_from_angular`, Legendre moment integrals, current
   evaluators) to handle a zero-weight member without contributing
   to the integral. Special-casing the zero weight everywhere costs
   architecture clarity.

2. **No scattering source coupling**. The inward sweep at μ = −1
   reads the scattering source `Q` (built from φ_ℓ moments via
   `(2ℓ+1)/2 · Q_ℓ · P_ℓ(−1)` per Eq. (3.432)), but it does NOT
   couple to a "neighbour ordinate" via the α-cascade because
   α_{−1/2} (the would-be α just below the auxiliary ordinate)
   does not exist. The auxiliary ordinate is structurally a leaf
   in the cascade DAG — adding it as a quadrature member would
   suggest a coupling that isn't there.

3. **Cross-section access pattern**. The on-the-fly sweep needs:
   - **Σ_t per radial cell** (the `Σ_i` in Eq. (3.433))
   - **The scattering+fixed source** evaluated at μ = −1, i.e.
     `Σ_ℓ ((2ℓ+1)/2) · Q_ℓ(r) · P_ℓ(−1)` where Q_ℓ is the standard
     ℓ-th Legendre moment of the source (scattering moments built
     from current iterate's φ_ℓ + fission moments + fixed source
     moments).
   - **The outer-face boundary condition** value at the most-inward
     ordinate (see open question in §8 about white BC).

   None of this requires modifying the quadrature dataset; all are
   either already available (Σ_t map) or recomputable from existing
   moment vectors. The matvec receives Σ_t, the moment vectors, the
   quadrature, and the BC — these inputs are sufficient.

4. **Tangent to ERR-026 root cause**. The current failure mode is
   that `MorelMontryAngularSweep` lacks any pole-face IC strategy
   on sphere, so it uses whatever default the abstract sweep
   contract provides (likely zero). Adding the inward sweep as
   on-the-fly construction inside the matvec is the local fix —
   no change to the quadrature contract or to other consumers.

### What state the inward sweep needs

Concretely the inward sweep operates on:

- `sigma_t[g, i]` — total cross-section per group per radial cell
  (shape `(ng, nx)`).
- `Q_moments[g, ℓ, i]` for ℓ = 0..L — Legendre moments of the source,
  computed from the current iterate's φ_ℓ via scattering + fission +
  fixed source. The inward sweep then forms
  `Q̄_i = Σ_ℓ ((2ℓ+1)/2) · Q_moments[g, ℓ, i] · (−1)^ℓ`
  using `P_ℓ(−1) = (−1)^ℓ`.
- `bc_face_value[g]` — the angular flux value at the outer face for
  μ = −1. For vacuum BC this is `0`. For reflective/white BC this
  is non-trivial (see open question §8).
- `radial_mesh.Δr_i, radial_mesh.r_{i±1/2}` — already standard
  inputs.

NO scattering-cross-section matrix `Σ_s,ℓ` is needed directly — the
moment-folded source is enough. NO ordinate-coupling state is
needed (no α-cascade reads/writes inside the inward sweep).

---

## 6. Pomraning structural-singularity cross-reference

**Source**: Pomraning (1989) *NSE* 101:330-340, file
`scratch/literature/Pomraning(1989)The Transport Equation in General Geometry.pdf`,
page 339, right column, last paragraph above the references.

Direct quote:

> "It was pointed out that if the bounding surface of the system is
> used as one of the coordinate surfaces and one considers a family
> of nonintersecting surfaces that starts with the bounding surface
> and progresses inward to fill the system, then these surfaces will
> eventually shrink to a surface with a zero area, namely a line or
> a point. Associated with this innermost degenerate surface will be
> one or more points where infinitely small radii of curvature are
> extant. ... A special case of this elliptical example is a sphere,
> where the innermost surface is simply a point. Hence, in general
> there will exist points on the innermost surface where the
> coefficients of the angular derivatives in the streaming term are
> infinite, since these coefficients contain the reciprocal of the
> radii of curvature [see Eqs. (68) and (75)]. As pointed out by the
> reviewer, these singular points could prove troublesome in analytic
> and, in particular, numerical work. This may well be the case, but
> such difficulties must be dealt with if one chooses to deal with
> the transport equation in curvilinear geometry. In these geometries,
> one always encounters such singular points, no matter what
> formulation is used to express the streaming term. **Prime examples
> of such singular points are found in the usual spherical and
> cylindrical geometry formulations where 1/r terms are extant and
> the attendant difficulties are well known, particularly in
> numerical treatments.**" (emphasis added)

### How the Carlson coupled-pole sweep is the canonical response

Pomraning's framing is geometric: r = 0 is **structurally singular**
in any curvilinear streaming operator because the angular-derivative
coefficients (the `(1−μ²)/r` term in the sphere Eq. (3.418) form)
contain 1/r. At r = 0 the coefficient diverges; the natural
discretisation must somehow handle this.

The naive engineering response would be **extrapolation**: pick
ψ_face(r = 0) by fitting a polynomial in r through the interior
nearby cells. This is what an incautious starting heuristic does;
it's also what produces the M-M wrong fixed point that ERR-026
diagnoses.

The Carlson coupled-pole response is **canonical** because it
sidesteps the singularity entirely: at the auxiliary direction
μ = −1 the singular `(1−μ²)/r` term is **identically zero**
(numerator vanishes), so the spatial sweep at this direction sees
no singularity at all. The equation tells you what φ̄_{1/2}(r = 0)
should be — you do not have to guess.

The cost is that you must solve the μ = −1 sweep first, then use
its result as the seed for the cascade at finite-weight ordinates
(where (1−μ²) > 0 and the singularity would otherwise be felt).
This is exactly the price Pomraning warns about: "difficulties must
be dealt with". The Carlson construction deals with it by EXPLOITING
the singularity's vanishing at μ = ±1 rather than trying to
regularise it at intermediate μ.

---

## 7. Implementation notes for the method-implementer

### Inputs

```
sigma_t      : np.ndarray, shape (ng, nx)        # cell totals
Q_moments    : np.ndarray, shape (ng, L+1, nx)   # scattering + fixed + fission Legendre moments
bc_outer     : np.ndarray, shape (ng,)           # angular flux at outer face for μ=−1
radial_mesh  : RadialMesh                        # provides Δr_i, r_{i±1/2}
```

(`bc_outer` for vacuum BC is just zeros. For reflective/white BC the
value comes from the outgoing current evaluated at outward ordinates;
this requires either a coupled outer iteration with the outgoing
sweep OR the white-BC integral evaluated at the auxiliary direction.
The latter is non-trivial — see §8.)

### Procedure

```python
def carlson_inward_sweep(sigma_t, Q_moments, bc_outer, dr, ng, L, nx):
    # For each group, run the inward sweep independently.
    # Result: cell-centred φ̄_{1/2,i} per group per cell.
    phi_aux = np.zeros((ng, nx))

    # Pre-build source at μ = −1: Q̄_i = sum_l ((2l+1)/2) * Q_l(r) * P_l(−1)
    # with P_l(−1) = (−1)**l
    l_indices = np.arange(L + 1)
    legendre_coeffs = (2 * l_indices + 1) / 2.0 * (-1) ** l_indices
    Q_bar = np.einsum("l,glx->gx", legendre_coeffs, Q_moments)  # (ng, nx)

    for g in range(ng):
        phi_face_outer = bc_outer[g]      # outgoing-face value at i = nx
        for k in range(nx - 1, -1, -1):   # sweep inward k = nx-1, ..., 0
            # Eq. (3.434) — solve for cell-centred φ̄_k
            denom = dr[k] * sigma_t[g, k] + 2.0
            phi_cell = (dr[k] * Q_bar[g, k] + 2.0 * phi_face_outer) / denom
            phi_aux[g, k] = phi_cell
            # Eq. (3.435) — step inward to next face
            phi_face_inner = 2.0 * phi_cell - phi_face_outer
            phi_face_outer = phi_face_inner   # roll for next iteration

    return phi_aux        # to be passed into the outward M-M α-cascade
```

The outward cascade then reads `phi_aux[g, i]` as the
`φ_{n−1/2, i}` value at `n = 1` (first finite-weight inward ordinate),
unchanged in subsequent ordinates via the α-recurrence.

### Caching opportunity — flag

The inward sweep's output `phi_aux[g, i]` depends on:

- `sigma_t` (geometry-iteration-stable across all Krylov inners)
- `bc_outer` (depends on BC kind; for vacuum it's zero forever)
- `Q_moments` (depends on the current iterate's φ_ℓ — changes
  with each Krylov inner because the right-hand side update
  rebuilds Q from new flux moments)

If the Krylov solver is structured so that `Q_moments` is held
fixed across multiple inner matvecs (it usually IS for fixed-source
inner iterations within a Krylov gradient step), the inward sweep
can be **cached** at the outer-iteration level and reused for every
inner matvec until `Q_moments` is refreshed.

Cache invalidation triggers:
- New outer iterate (φ_ℓ changes → Q_moments changes → invalidate).
- BC change (rare; usually only at problem setup).
- Mesh change (never within an inner solve).

For a single matvec the marginal cost of the inward sweep is
`O(ng · nx)` operations — fast. The caching is only worthwhile if
the outer iteration spans many matvecs and `Q_moments` is stable
across them. For an L1 sphere MMS gate the caching is probably
optional, but the architecture should ALLOW for it (i.e., the
strategy object should expose an `invalidate()` method or be
re-evaluated only when input fingerprints change).

### Output shape

Per Phase D brief: "per-group pole-face value, shape `(ng, n_outgoing)`".
This shape is too narrow given §4's analysis — the outward cascade
needs the **full radial profile** `phi_aux[g, i]`, not just `phi_aux[g, 0]`.
Recommend the strategy expose:

```python
class CarlsonCoupledPole(PoleFaceInitialCondition):
    def evaluate(self, ...) -> np.ndarray:        # shape (ng, nx)
        # Returns cell-centred φ̄_{1/2,i} for all radial cells.
        ...
```

This may force a small extension to the `PoleFaceInitialCondition`
Protocol's return contract. Coordinate with the architecture spec
in the §16A.3 BC trace contract docs before locking the shape.

---

## 8. Open questions / unknowns

### Multi-region σ_t step discontinuities

Hébert §3.9.4 derives Eqs. (3.432)-(3.435) for a single material
constant Σ_i per sub-mesh. The discretisation is per-cell, so a
step discontinuity at a cell interface is handled trivially by
using the correct Σ_i in each cell. **No special treatment needed**
as long as the radial mesh respects material boundaries. **Unknown**:
whether the cell-averaged source `Q̄_i` should also be material-aware
(it should — `Q_moments` are computed per cell, so this is fine).

### White BC at outer face

This is the hardest open question. Hébert's §3.9.3 cylindrical
discussion (Eq. (3.415)+(3.416)) handles white BC by computing the
outgoing current `J^+ = 4 Σ_p Σ_q 𝒲_{p,q} μ_{p,q} φ_{p,q,I+1/2}`
and uniformly redistributing it back as an incoming flux
`φ_{p,q,I+1/2}^{in} = J^− / (4 Σ_{p'} Σ_{q'} 𝒲_{p',q'} μ_{p',q'})`
across all positive-μ ordinates.

For the sphere, the analog would compute J^+ from outgoing ordinates
and supply an incoming flux at every inward ordinate AND at the
auxiliary μ = −1 ordinate. **Hébert is silent on whether the
auxiliary ordinate gets a non-zero white-BC value at the outer
face**, because his §3.9.4 sphere text covers only vacuum BC
(text above Eq. (3.432) says "particles entering the external
surface" without specifying source magnitude).

**Implementation default**: for vacuum BC use `bc_outer[g] = 0`.
For reflective/white BC, propose computing `bc_outer[g]` from the
white-BC redistribution applied at μ = −1, which under the
uniform-redistribution assumption equals the same `J^− /
(quadrature normalisation)` value used for all positive-μ
ordinates. Flag this as a design decision for the method-implementer
and verify against the cylindrical analog Eq. (3.416). The
numerics-investigator can probe empirically by running a sphere
white-BC MMS and comparing `bc_outer = 0` vs `bc_outer = J^−/Σ𝒲μ`
fixed points.

### Cylindrical analog — why does cylinder Gate 1.1 pass?

The ERR-026 empirical observation: cylindrical Gate 1.1 MMS passes
under M-M angular closure, sphere Gate 1.1 fails. Hypothesis: the
cylindrical α-cascade has more telescoping structure (the α-recursion
runs over BOTH the axial-ξ levels AND the polar-ω rotations, with
zero-weight points injected at ω = π for each ξ level — Hébert
Eqs. (3.396), (3.408)-(3.410)).

The cylindrical zero-weight sweep at ω = π IS structurally analogous
to the sphere μ = −1 sweep:

```
2π μ̃_p (φ̄_{p,i+1/2} − φ̄_{p,i−1/2}) + ΔS_i Σ_i φ̄_{p,i}
    = ΔS_i Q̄_{p,i}                                  # Eq. (3.408), cylindrical
```

vs

```
−(φ̄_{i+1/2} − φ̄_{i−1/2}) + Δr_i Σ_i φ̄_i
    = Δr_i Q̄_i                                      # Eq. (3.433), sphere
```

Both are inward-direction DD recurrences without α-redistribution.
**Hébert DOES document the cylindrical Carlson starting direction**
(p. 140, Eqs. (3.407)-(3.410)) — it is NOT a sphere-only construct.

So why does cylinder pass empirically? Two possible explanations:

1. **Implementation already has it for cylinder**: ORPHEUS's
   cylindrical angular sweep already implements the ω = π starting
   direction (Hébert §3.9.3 IS the canonical reference for the
   current cylindrical implementation), so the IC IS being seeded
   correctly. The sphere implementation simply lacks the analogous
   §3.9.4 starting direction and falls through to a zero-IC default.

2. **α-dome telescoping**: the cylindrical α-coefficients run over
   a 2D angular index `(p, q)` rather than the sphere's 1D index
   `n`, giving the cascade more degrees of freedom to absorb a
   wrong seed. The sphere cascade is more "rigid" — a wrong seed
   propagates more directly to the converged fixed point.

The ANSWER is probably (1) based on the architectural diagnosis in
the Phase D brief, but flag (2) for the numerics-investigator: if
the cylindrical implementation does NOT explicitly run an inward
sweep at ω = π, then (2) is the active mechanism and the cylinder
result is coincidentally robust rather than canonically correct.

### Multigroup with down-scattering

The inward sweep is per-group. For multigroup with downscatter, the
source moments `Q_moments[g, ℓ, i]` already include downscatter
contributions from groups g' < g, so the per-group sweep handles it
implicitly. **No additional logic needed**. Up-scattering or thermal
coupling would require an outer iteration over groups; the inward
sweep stays per-group inside that loop.

### Quadrature order N — does the seed quality depend on N?

The auxiliary ordinate μ = −1 is INDEPENDENT of N (it's always the
endpoint of [−1, 1]). The α_{1/2} = 0 seed is always exact. The
seed quality therefore depends only on:

- The accuracy of `Q_moments` (depends on the current iterate's
  scalar-flux moments, which depend on N indirectly via the
  quadrature integration).
- The spatial discretisation `Δr` (DD is second-order accurate in
  Δr).

**No N-dependence in the inward sweep itself**. This is a feature:
the Carlson IC strategy is N-agnostic and works identically at S_4
and S_64.

---

## V&V hooks for this memo

When Phase D ships, the method-implementer should request:

1. **L0 term verification**: build a unit test that constructs a
   flat-ψ probe (uniform Σ_t, flat consistent Q), runs the inward
   sweep alone, and asserts every `phi_aux[g, i] == C` to machine
   precision (the §3 algebra proves this is the exact closed-form
   answer).
2. **L1 sphere MMS**: re-run Gate 1.1 with the Carlson IC strategy
   wired into `MorelMontryAngularSweep` and verify the
   second-order spatial convergence rate at fixed N.
3. **Cylindrical regression**: run cylindrical Gate 1.1 with the
   analogous Carlson IC on cylinder (if needed — see open question)
   and confirm no regression.

Reference Pillar C of `vv-principles` — closed-form analytical
reference. The Carlson sweep is closed-form per Hébert §3.9.4; the
flat-ψ assertion in step 1 is a closed-form fixed-point check.

---

## Pointers

- **This memo's companion**: `sphere_sn_pole_closure_canonical.md`
  (background on Hébert §3.9.4 + Bailey-Morel-Chang 2010 + canonical
  references).
- **Primary source — LOCAL**:
  `scratch/literature/Hebert(2009)Chapter3.pdf` pp. 141-144,
  Eqs. (3.418)-(3.439).
- **Pomraning structural-singularity quote — LOCAL**:
  `scratch/literature/Pomraning(1989)The Transport Equation in General Geometry.pdf`,
  page 339.
- **ORPHEUS code touchpoints** (verify before recommending — these
  may have moved):
  - `src/orpheus/sn/angular_sweep/morel_montry.py`
  - `src/orpheus/sn/reduced_operator.py` (Bailey docstring citation
    correction — see sphere_sn_pole_closure_canonical.md memo).
  - PoleFaceInitialCondition Protocol — to be created by Phase D.
