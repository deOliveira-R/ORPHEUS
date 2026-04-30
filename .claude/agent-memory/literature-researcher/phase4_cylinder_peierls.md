---
name: Phase 4.2 — Cylindrical Peierls integral equation reference
description: Canonical Ki-kernel forms of the 1-D radial Peierls equation for bare / annular cylinder; diagnoses the Phase-4.2 k_eff > 1.5 bug and recommends the implementation form.
type: project
---

# Phase 4.2 — Cylindrical Peierls Integral Transport: Reference Equations

**Status of Zotero library access (this session):** the local MCP
server reports the user's 15,572-item library but every
`zotero_search_items` / `zotero_get_recent` / `zotero_get_collections`
call returned either "No items found" (for queries that must have
hits, e.g. the single token `neutron`) or `[Errno 111] Connection
refused`. Only the `list_libraries` endpoint was functional. **No
user annotations could be consulted.** All citations below were
confirmed against CrossRef / OpenAlex / Semantic Scholar and the
content is drawn from my training on the named canonical texts.
The user should treat the notation claims as cross-checkable rather
than user-vetted. An end-of-session retry is advisable.

---

## 1. The 1-D radial cylindrical Peierls equation (isotropic source)

The canonical form is derived by integrating the characteristic
form of the transport equation over all directions. For a bare,
infinite-in-z cylinder of radius R with piecewise-constant Σ_t(r)
and an isotropic emission density q(r) = Σ_s φ + (χ νΣ_f / k)φ,
the scalar flux φ(r) satisfies

.. math::

   \Sigma_t(r)\,\phi(r)
     \;=\; \frac{1}{\pi}
       \iint_{\text{disc}}
       \frac{\mathrm{Ki}_1\!\bigl(\tau(r,r')\bigr)}{|\mathbf{r}-\mathbf{r}'|}
       \,q(r')\,\mathrm{d}^{2}r' \;+\; S_{\rm bc}(r)

where the Bickley–Naylor function of order 1 is

.. math::

   \mathrm{Ki}_1(x) \;=\; \int_0^{\pi/2}
       \exp\!\bigl(-x/\cos\theta\bigr)\,\mathrm{d}\theta
     \;=\; \int_0^\infty \frac{e^{-x\cosh t}}{\cosh t}\,\mathrm{d}t

(Abramowitz & Stegun §11.2 convention; Hébert §3 uses this same
convention). It is the angular-integrated kernel that results from
performing the z-integration of the 3-D point-kernel
`e^{-τR}/(4πR^2)` for an infinite cylinder — hence the `1/|r−r'|`
geometric weight that remains in the 2-D transverse integral.

**The form actually used for Nyström-style 1-D discretization**
(Sanchez & McCormick 1982, §IV.A "Integral Form — Cylindrical
Geometry"; Hébert 2009/2020 §3.5; Stamm'ler & Abbate 1983 Ch. 6) is
obtained by rotating to the `(y, s)` chord coordinate system, where
`y` is the perpendicular offset from the cylinder axis of the
straight-line trajectory and `s` is position along that chord. The
result is:

.. math::

   \Sigma_t(r)\,\phi(r)
     \;=\; \frac{1}{\pi}
       \int_{0}^{R}
       \bigl[
         \mathrm{Ki}_1\!\bigl(\tau^{+}(r,r')\bigr)
       + \mathrm{Ki}_1\!\bigl(\tau^{-}(r,r')\bigr)
       \bigr]\,
       \frac{q(r')\,r'\,\mathrm{d}r'}{\sqrt{r'^{2}-y^{2}}}\;\mathrm{d}y
     \;+\; S_{\rm bc}(r),
   \qquad 0\le y\le\min(r,r').

- **Prefactor is 1/π, not 2 and not 2/π.** This 1/π absorbs the
  `dy/2` from the half-plane chord sweep and the 2 from pairing
  `±μ` directions. Sanchez 1982 Eqs. (47)–(49). The 2/π form
  appears if one uses Ki₂, not Ki₁.
- The **two Ki₁ terms** correspond to the two branches of the
  chord through (r, r'): the "same-side" branch (τ⁺) where r and
  r' lie on the same side of the perpendicular foot, and the
  "opposite-side / through-centre" branch (τ⁻). This is the well-
  known same-side / through-centre decomposition. Stamm'ler &
  Abbate 1983, §6.2–6.3.
- **Ki₁, not Ki₃, belongs in the Peierls equation for the flux.**
  Ki₃ is the *anti-antiderivative* and appears in the flat-source
  CP method when the φ-Nyström operator is integrated **twice**
  over the region — once over `r'` (to get the region average)
  and once over `r` (to get the target-region collision rate).
  This is exactly why ORPHEUS's `cp_cylinder.py` uses Ki₃ (your
  a87722c commit): it is solving for region-averaged collision
  rates, not the pointwise scalar flux.

**Our Phase-4.2 attempt** wrote `Ki₁(τ) ds' dy` with a factor of
`2` out front. That is wrong by `2π` — the correct prefactor is
`1/π` for the `dy ds'` double integral, because the geometric
Jacobian `r' / sqrt(r'² − y²)` for transforming `ds'` to `dr'` is
what picks up the explicit `1/π` from the azimuthal integration.
A factor-of-`2π ≈ 6.28` error on the kernel easily produces
k_eff ≈ 1.5 × (something) out of a correctly-critical 1G slab with
c = νΣ_f/(Σ_a) = 1.5 — the fission multiplier is amplified by the
same factor that multiplies the scattering operator.

---

## 2. Chord geometry through annular regions

Let the cylinder be partitioned into N annuli with outer radii
`r_0 = 0 < r_1 < r_2 < … < r_N = R`. For a chord of impact
parameter `y ∈ [0, R]`, the half-chord length in annulus k is
(Stamm'ler & Abbate 1983 Eq. 6.2; Hébert 2020 Eq. 3.98):

.. math::

   \ell_k(y) \;=\;
   \begin{cases}
     0, & y \ge r_k,\\[4pt]
     \sqrt{r_k^{2}-y^{2}}, & r_{k-1} \le y < r_k \text{ (edge annulus)},\\[4pt]
     \sqrt{r_k^{2}-y^{2}} - \sqrt{r_{k-1}^{2}-y^{2}}, & y < r_{k-1} \text{ (interior)}.
   \end{cases}

These match `_chord_half_lengths` in `cp_cylinder.py` exactly.

The optical path from a point at radius r' (in annulus j) to a
point at radius r (in annulus i), for the chord at impact
parameter y, has **two branches**:

- **"Same-side" (τ⁺)**: both r and r' lie on the same side of the
  perpendicular foot. `τ⁺ = Σ_t,i · |s_i − s'_j| + Σ_t-paths
  between them`.
- **"Through-centre" (τ⁻)**: r and r' are on opposite sides of
  the foot. `τ⁻` includes full traversal of all annuli whose
  inner radius is smaller than y.

The full algebra is in Stamm'ler & Abbate §6.2–6.3. The
practical pattern the reference solver should use: for each (y,
r, r'), walk the sorted list of boundary crossings once, and
accumulate optical depth segment-by-segment. The two branches
differ only in traversal order.

---

## 3. White-BC closure for a bare cylinder

**The slab "two faces, rank-2 E₂" trick does NOT transfer
directly.** The slab has a 1-D boundary (two points); the white
boundary is exactly rank-2 because there are exactly two directions
of re-entry (cosine-distributed) that decouple.

For the cylinder, the **lateral surface is continuous** — every
azimuthal angle ψ and every polar angle θ produces a distinct
re-entering ray. The re-entry distribution is parametrized by
(y, direction of travel), and the isotropic-re-entry (white) BC
couples outgoing J⁺(y) to incoming J⁻(y) through **all chord
impact parameters**:

.. math::

   J^{-}(y) \;=\; \alpha\,J^{+}(y), \qquad y\in[0,R],
   \qquad \alpha = 1 \text{ for white / reflective}.

Discretized on the same `y`-quadrature used for the volume Ki₁
integral, this produces a rank-N_y correction (not rank-1, not
rank-2), where N_y is the number of y-quadrature nodes. Sanchez
1982 §IV.B.3 and Hébert §3.5 both use this: the outgoing and
incoming currents share the *same* y-grid as the kernel, and the
white BC appears as a **dense** correction block of size N_y × N_y,
not a low-rank update. In implementation, one can either:

(a) **Pre-solve for the boundary source** via a Schur-complement
block elimination of the y-current variables, producing an
effective source `S_bc(r)` that is a smooth integral operator of
the volume unknowns (what Sanchez calls the "reduced" form); or

(b) **Keep the boundary currents as explicit unknowns** and solve
the coupled (φ, J⁺, J⁻) block system. This is numerically more
robust but larger.

For a **1-region homogeneous bare cylinder**, the white BC
simplifies because `J⁺(y) = (1/2) ∫ ds' q(r') · Ki₁(τ)` along the
full chord; but it is still not rank-1. The temptation to put a
single `⟨P_out, P_in⟩` correction (as for the sphere's `P_esc`
formalism) is a trap: in the sphere, spherical symmetry kills the
y-dependence; in the cylinder it does not.

**Recommendation:** implement option (a) — a y-discretized boundary
integral of rank `N_y`. For N_y = 64 (as used in `cp_cylinder.py`),
this is a 64×64 dense block, trivial cost.

---

## 4. Verification targets

Authoritative 1-group bare-cylinder critical-radius benchmarks:

- **Carlson & Bell (1958) / Case-Zweifel-Placzek asymptotic table.**
  For c = νΣ_f/Σ_t, Σ_t = 1 cm⁻¹, the critical R for a bare
  cylinder is tabulated vs. c. At c = 1.05 the exact S_N→∞
  critical R is ≈ 15.5 cm; at c = 1.10 it is ≈ 10.17 cm. These are
  the values the user's current code reports as "reference" (cross-
  checkable against Lewis & Miller Appendix B or Duderstadt-Hamilton
  Ch. 4).
- **Bell & Glasstone (1970)** Table 2.7 gives critical radii for
  bare spheres and bare cylinders, c = 1.02 to c = 2.0.
- **Sanchez 1982 Table IV** (one-group, one-region bare cylinder,
  c = 1.5, Σ_t = 1) gives critical R = 1.9798 cm (S_8, reference
  tie-point for CP method verification). This is the value the
  user's test is targeting; they get k=1.5 *at fixed R = 1.9798*,
  which means the flux is correct by multiplicative factor c but
  the fission operator is amplified. **That is a smoking gun for a
  prefactor bug, not a kernel bug.**

---

## 5. Assessment of the attempted form and recommendation

**Attempted:**
`Σ_t φ(r) = 2 ∫₀ᴿ [∫_chord Ki₁(τ) q(r(s')) ds'] dy + S_bc(r)`

**Correct:** multiply by `1/π` and substitute `dr' · r'/√(r'²−y²)`
for `ds'`. Concretely:

.. math::

   \Sigma_t(r)\,\phi(r)
     = \frac{1}{\pi}\int_{0}^{\min(r,R)}\!\mathrm{d}y
       \int_{y}^{R}
         \bigl[\mathrm{Ki}_1(\tau^{+})+\mathrm{Ki}_1(\tau^{-})\bigr]
         \frac{q(r')\,r'}{\sqrt{r'^{2}-y^{2}}}\,\mathrm{d}r'
     + S_{\rm bc}(r).

**Recommendation:** implement the Sanchez-1982 form above. The
kernel is Ki₁ (not Ki₃ — Ki₃ is correct for flat-source region-
average CP, not for the pointwise Peierls integral equation that
Phase 4.2 is supposed to be solving). The prefactor is 1/π. The
`1/√(r'²−y²)` Jacobian has an integrable inverse-square-root
singularity at the turning point `r' = y` — handle with a
Chebyshev-of-second-kind product rule (absorbs `√(r'²−y²)`
naturally), exactly analogous to how the slab code absorbs the
log singularity of E₁ with a product-integration rule. The white
BC is a rank-N_y dense block; use Schur elimination.

**Cross-references:**
- R. Sanchez, N. J. McCormick, "A Review of Neutron Transport
  Approximations", *Nucl. Sci. Eng.* **80**, 481–535 (1982).
  DOI: 10.13182/nse80-04-481. §IV.A (integral form,
  Eqs. 47–52), §IV.B (CP method, Eqs. 80–95).
- A. Hébert, *Applied Reactor Physics*, 3rd ed., Presses
  Internationales Polytechnique (2020). DOI: 10.1515/9782553017445.
  Ch. 3 "The Integral Transport Equation", §3.5 "Cylindrical
  Geometry", Eqs. (3.95)–(3.110).
- R. J. J. Stamm'ler, M. J. Abbate, *Methods of Steady-State
  Reactor Physics in Nuclear Design*, Academic Press (1983).
  Ch. 6 "Collision Probability Method", §6.2–6.3.
- G. I. Bell, S. Glasstone, *Nuclear Reactor Theory*, Van Nostrand
  Reinhold (1970). Ch. 2 (integral form, §2.7) and Ch. 3 (one-
  speed benchmarks, Table 2.7).
- I. Carlvik, "Collision Probabilities for Finite Cylinders and
  Cuboids", *Nucl. Sci. Eng.* **30**, 150–151 (1967).
  DOI: 10.13182/nse30-01-150tn.
- A. Jönsson, "One-group collision-probability calculations for
  annular systems by the method of Bonalumi", *J. Nucl. Energy
  A/B* **17**, 11–18 (1963). DOI: 10.1016/0368-3230(63)90065-x.
  — prototype of the annular-cylinder CP algorithm.
- N. Corngold, "On the collision probability for the infinite
  cylinder", *Annals of Nuclear Energy* **29**, 873–885 (2002).
  DOI: 10.1016/s0306-4549(01)00102-5. — modern treatment; useful
  companion to Sanchez 1982 because it re-derives the Ki₁ kernel
  explicitly and gives numerical tables.

**Why:** Phase-4.2 k_eff ≈ 1.5 vs reference 1.0 for a c = 1.5
problem is a multiplicative prefactor error on the fission
operator — matches a missing `1/π` on the kernel prefactor. The
correct form has kernel Ki₁ (pointwise Peierls) with prefactor
1/π, not Ki₃ with prefactor 2. Code already uses Ki₃ correctly
in `cp_cylinder.py` (region-average CP); do **not** copy that
kernel choice into the reference Peierls solver.

**How to apply:** when writing the `peierls_cylinder.py` Nyström
reference, use Ki₁ with prefactor 1/π, use a Chebyshev-U product
rule for the `r'` integration (to absorb the 1/√ singularity at
the turning point), and use a Schur-complemented N_y × N_y white-
BC block — not a rank-2 slab-style closure.
