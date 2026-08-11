# Q5.6.4 attempt 2 — adversarial QA review of the C1→C7 chain

**Reviewer:** `qa` sub-agent, 2026-08-11. Host env, `.venv/bin/python`.
**Subject:** `.claude/plans/q64_tau_partition_memo.md` §9bis, read against
`scratch/q64_attempt2_probe_outputs.md`,
`orpheus/sn/sweep/pole_angular_closure.py`, and
`orpheus/derivations/discrete/sn/angular_differencing.py`.

**Probes written for this review** (in `$CLAUDE_JOB_DIR/tmp`, ephemeral —
regenerate from the descriptions below):
`qa_c1_symbolic.py`, `qa_c6_harmonics.py`, `qa_c5_amplification.py`,
`qa_c2c3_kappa.py`, `qa_probeC_live.py`, and the two mutation plugins
`_qamut_c7_omega.py` / `_qamut_c7_faithful.py`.

> ⭐ **Headline.** **C1 SURVIVES** (independently re-derived symbolically, and
> confirmed verbatim by Hébert Eq. 3.157 / 3.393). **C2 survives.** **C3 is
> FALSE AS STATED** — its universal quantifier is an even-`M` parity artefact
> of the four orders probed; at odd `M` (`n_φ = 6, 10, 14, …`, all legal rules)
> **0 of `M+1`** edges fail. **C4 survives and is the memo's strongest
> finding.** **C5 is mis-scoped in BOTH directions** — the end-to-end
> amplification is exactly `1.0` for chord, arc AND `τ≡½`, only the retired
> absorber *dissipates* (40×), and the positivity "cost" is inverted (`τ≡½` is
> the SAFEST derived candidate). **C6 BREAKS**: its `ξ` leg is not a P1 mode
> here (`J_φ ≡ 0` by reflection symmetry; BMC Eq. 1 carries one cosine) and on
> a folded rule `quad.mu_y` samples `|ξ|`, not `ξ` — the `ξ`-WEIGHTING is
> exonerated but the BASIS is not, and under the physically-allowed family
> `{cos mω}` the ranking **INVERTS**: the LANDED arc convention becomes best,
> by ≈2× on *every* harmonic. **C7 therefore does not follow**, and §9bis.9's
> literature defence is refuted by BMC naming `τ = ½` "the diamond scheme" in
> their own cylinder section (Eq. 53) and proving the diamond only
> leading-order-accurate.
>
> ⭐⭐ **And the memo's PENDING decisive row is now FILLED (§6).** `τ ≡ 
> ½` measures **`1.0181e-01`** — it clears the `1.2e-1` gate but is **1.54×
> worse** than the retired chord+absorber (`6.5934e-02`, reproduced exactly)
> and only 1.24× better than the LANDED arc (`1.2676e-01`, also reproduced
> exactly). The four candidates' ranking is in **perfect rank correlation with
> the recurrence's error amplification** and **inverted against closure
> accuracy**, so that metric cannot adjudicate a chart at all. The unexamined
> third option — fix the starting-direction SEED and keep the arc closure — is
> the only candidate that could recover the accuracy without giving up BMC's
> exactness property, and it is one cheap experiment away.

---

## 1. C1 — the conservative cylindrical form

### 1.1 Verdict on the algebra: CORRECT, independently

I re-derived it by a route that touches neither the memo nor ORPHEUS:
build `Ω` in **lab Cartesian** components, transport along a straight ray,
and differentiate `ψ(r, ω)` by the chain rule (`qa_c1_symbolic.py`).

```
dr/ds     = sin(theta)*cos(Phi - phi_s)                    residual vs eta      : 0
domega/ds = -sin(theta)*sin(Phi - phi_s)/r                 residual vs -xi/r    : 0
```

with `ω := Φ − φ_s` (Φ = the fixed LAB azimuth of `Ω`'s in-plane
projection), giving the non-conservative operator
`Ω·∇ψ = η ∂ψ/∂r − (ξ/r) ∂ψ/∂ω`. Then, symbolically:

| candidate form | `− (non-conservative form)` |
|---|---|
| `(η/r)∂(rψ)/∂r − (1/r)∂(ξψ)/∂ω`  ← **C1** | **`0`** |
| `(η/r)∂(rψ)/∂r − (ξ/r)∂ψ/∂ω` (plausible-wrong) | `+ψ·sinθ·cosω/r` — a spurious `(η/r)ψ` |
| angular flux coefficient `η` not `ξ` (plausible-wrong) | non-zero |

The mechanism is one identity, `∂ξ/∂ω = ∂(sinθ sin ω)/∂ω = sinθ cos ω = η`,
which is exactly what makes the `(η/r)ψ` from the radial term cancel the
`(η/r)ψ` from the angular term. **C1's displayed equation is right.**

### 1.2 The face coefficient IS `ξ` at the face, exactly

By the fundamental theorem of calculus, over one angular cell
`[ω_a, ω_b]` times a polar cell of `μ_z`-measure `W`:

```
∫∫ −(1/r) ∂(ξψ)/∂ω  =  −(W/r) [ ξ(ω_b)ψ(ω_b) − ξ(ω_a)ψ(ω_a) ]     (exact)
```

so the face coefficient is `W·ξ(ω_face)` with **no `κ` and no extra
geometric factor**. C1's second sentence is confirmed.

`[M]` And `W` is identified: on `folded_product(4, n_φ)` level 0 the
per-ordinate weight is exactly `w_m = w_gl · Δω` with
`w_gl = w_m/Δω = 0.695709690` at **every** `n_φ ∈ {6,8,10,16,32,64}`
(`qa_c2c3_kappa.py`), and `0.695709690 = 2 × 0.347854845` = twice the
`gauss_legendre(4)` weight — the factor 2 being the σ_y fold. So probe B1's
"the fitted `w_gl` … is the level's polar weight" is CORRECT and the
constant is the *folded* polar weight.

### 1.3 ⛔ FINDING C1-a — one clause of C1 is a DEFINITION wearing a derivation

C1 says: *"so the half-angle faces **sit at** the geometric arc edges and
the exact coefficient there is `w_gl·ξ(ω_edge)`"*. Only the second half is
derived. The conservative form says nothing about **where** the faces are —
it says that *at whatever `ω` you choose as a face*, the exact coefficient is
`w_gl·ξ(that ω)`. Where the faces go is a free choice (it is the partition,
which is the whole thing under debate). Splicing the two halves into one
"so" is the §2 `[M]`-scope defect the memo's own companion rule
(`plan-authoring` §2) warns about, and it is load-bearing here because
§9bis.2 uses "α ∝ ξ(e_arc)" to *reinstate P-C ("the PARTITION was right")*.

The honest reinstatement is narrower and is still enough:
`α`'s recursion sums `w_m η_m = w_gl · Δω · η_m`, which is the **midpoint
rule** for `∫_cell η dω` over the cell *centred on the ordinate with width
Δω* — and on an equispaced-ω folded level that cell IS the arc cell. So
`α`'s construction does implicitly use the arc partition, via the weights.
That argument does not need C1's "faces sit at" clause.

### 1.4 ⛔ FINDING C1-b — C2's corroboration of C1 is PROCEDURAL, not structural

The brief asked whether C2's corroboration is genuine or circular. **It is
not independent**, for a sharp reason:

* C1's exact face coefficient is `ξ(ω_face)`, i.e. the **antiderivative in ω
  of `η`**.
* `α`'s recursion accumulates `−w_m η_m`, i.e. a **quadrature of the same
  integral**.

So "α ∝ ξ(e_arc) up to κ" is a statement about the *quadrature error of one
integral*, and it is true **whatever the correct face coefficient is** —
because `α`'s construction never consults C1. If C1 were wrong, production
`α` would still measure `κ·w_gl·ξ(e_arc)`. Both sides share the upstream
identity `dξ/dω = η`: this is the ERR-032 / `vv-principles` #7 shape
("two derivations agree because they share the identity, not because either
is right").

Worse, the *pre-existing* gate for C2 —
`tests/sn/sweep/curvilinear/test_alpha_closed_form.py:132`
`test_production_alpha_equals_xi_at_the_half_angle_boundaries`
(`@l1`, `@verifies("alpha-cylindrical")`) — states its own reference as
*"Hebert 2009 Eq. 3.399 (alpha IS the tangential cosine at a real half-angle
boundary)"*. That is **C1's claim**. Citing C2 as evidence for C1 is
therefore comparing a claim against its own reference.

⟹ **C1 survives**, but on §1.1's derivation, **not** on C2. Anyone rewriting
the memo must not carry "C2 corroborates C1" forward as independent evidence.

`VERDICT: C1's equation and its face-coefficient claim are CORRECT (independently re-derived symbolically; the face coefficient is exactly W·ξ(ω_face), no κ). TWO defects in how C1 is STATED: (a) "the faces sit at the arc edges" is a definitional choice spliced into a derivation with a "so", and it is the clause §9bis.2 leans on; (b) C2's corroboration is procedurally-not-structurally independent (both sides ride dξ/dω = η, and C2's own committed gate cites C1's claim as ITS reference) — so C2 must not be carried as evidence for C1.`

---

## 2. C6 — the "P1 closure defect at the arc faces" instrument

### 2.1 (c) first, because it is the cleanest: the RANKING IS STABLE under reweighting

`qa_c6_harmonics.py` §C6-B, `folded_product(4,64)` level 0, the memo's own
two legs (`η`, `ξ`), five weightings:

| convention | W0 unweighted max | W1 `ξ(e)`-wtd (memo's) | W2 `\|η(e)\|`-wtd | W3 uniform L2 | W4 sum of legs |
|---|---|---|---|---|---|
| (1) chord | 6.637e-01 | 3.374e-01 | 1.684e-01 | 4.698e-01 | 6.637e-01 |
| (1c) chord+absorber | 1.631e-02 | 5.469e-04 | 8.291e-03 | 3.400e-03 | 1.814e-02 |
| (2) arc / ω-mid **LANDED** | 1.415e-01 | 7.157e-02 | 4.458e-02 | 1.127e-01 | 1.415e-01 |
| (C) `τ ≡ ½` | **1.226e-03** | **4.045e-04** | **6.234e-04** | **8.805e-04** | **1.839e-03** |

`τ ≡ ½` wins **all five**. So the specific charge in the brief — "does the
`ξ(e_arc)` weighting bias the ranking?" — is **REFUTED**: it does not. The
`ξ` weighting is also defensible on its own terms (`α ∝ ξ(e_arc)` is the
coefficient the face value is multiplied by in the balance).

**But that is the wrong axis.** The bias is not in the weighting. It is in
the **basis**.

### 2.2 ⛔ (a) FINDING C6-a — `ξ` is NOT a mode of this problem, and on the FOLDED rule it is not even `ξ`

Two independent reasons, both measured.

**Reason 1 — the reflection symmetry kills the `ξ` coefficient.** For a 1-D
cylinder, reflection across the plane spanned by `e_r` and `e_z` maps
`ξ → −ξ`, `η → η`, `μ_z → μ_z`, and leaves geometry, materials and
(isotropic/Legendre) scattering invariant. Hence `ψ(r, η, ξ, μ_z)` is
**even in `ξ`**, so `J_φ = ∫ξψ dΩ ≡ 0`. In the P1/diffusion limit

```
ψ = (1/4π)(φ + 3 J·Ω),      J = J_r e_r      ⟹      ψ|level = A + B·η
```

with the `ξ` coefficient **identically zero**. So "a P1/diffusion-limit flux
is affine in the direction cosines" is true of P1 in general 3-D and
**false as a licence to feed `ξ` here**: in this geometry the diffusion limit
has no `ξ` component to reproduce.

**Reason 2 — on a folded level, `quad.mu_y` is `|ξ|`, not `ξ`.** `[M]`

| rule | `min ξ` | `Σ w ξ` | `Σ w \|ξ\|` |
|---|---|---|---|
| `folded_product(2,8)` | **+0.3125** | **+6.7029** | +6.7029 |
| `product(2,8)` (unfolded) | −0.8165 | **0.0000** | +6.1927 |

Every folded node has `ξ > 0`, so the array the probe feeds is the sampling
of `|ξ| = sinθ|sin ω|` — an **even** function of `ω`, not a linear function
of `Ω`. And production only ever sees folded levels: the unfolded full-circle
level is REFUSED at
`pole_angular_closure.py:890-901` (the monotone-arc guard).

`[M]` Fourier-cosine content on the arc (`qa_c6_harmonics.py` §C6-C,
`c_m` for `m = 0…6`):

| function | `c_0` | `c_1` | `c_2` | `c_3` | `c_4` | `c_5` | `c_6` |
|---|---|---|---|---|---|---|---|
| `η/sinθ = cos ω` | 0 | **+1.0000** | 0 | 0 | 0 | 0 | 0 |
| `ξ/sinθ = sin ω` | +0.6366 | 0 | −0.4244 | 0 | −0.0849 | 0 | −0.0364 |
| `ω` (what `τ≡½` is exact on) | +1.5708 | −1.2732 | 0 | −0.1415 | 0 | −0.0509 | 0 |

`cos ω` is **exactly one harmonic — the P1 mode**. `sin ω` (and `ω`) spread
over *all* harmonics with an `m⁻²` tail; as functions of `η` they are
`sqrt(sin²θ − η²)` and `arccos(η/sinθ)`, both **sqrt-singular at the level
endpoints**. No first-order reconstruction in `η` can be uniformly accurate
on them, which is the entire content of the arc convention's poor `Δξ`
column — it is generic behaviour on a singular function, not a partition
defect.

### 2.3 (a) The right local model, and the ranking under it — the RANKING INVERTS

Because `ψ` on a folded level is even in `ω`, its angular dependence is a
**cosine series**, equivalently a **Chebyshev series in `η/sinθ`**:

```
ψ|level = Σ_m c_m cos(m ω) = Σ_m c_m T_m(η / sinθ)
```

so the honest test family is `{cos mω}`, not `{η, ξ}`. `[M]`
`qa_c6_harmonics.py` §C6-A, `folded_product(4, n_φ)` level 0,
`max|ψ̂(e) − f(e)|` over the level's faces, seeded with the TRUE start-face
value:

`n_φ = 64` (`M = 32`):

| mode | chord | chord+absorber | **arc LANDED** | `τ ≡ ½` |
|---|---|---|---|---|
| `cos 0` (isotropic) | 6.7e-16 | 0 | 8.9e-16 | 0 |
| **`cos 1` = the P1 mode** | **5.4e-15** | 3.605e-03 | **2.1e-15** | 2.412e-03 |
| `cos 2` | 4.887e-03 | 1.435e-02 | **4.854e-03** | 9.677e-03 |
| `cos 3` | 1.171e-02 | 3.205e-02 | **1.132e-02** | 2.188e-02 |
| `cos 4` | 2.185e-02 | 5.636e-02 | **1.994e-02** | 3.918e-02 |
| `ξ` (memo's leg) | 6.637e-01 | 1.631e-02 | 1.415e-01 | **6.131e-04** |
| `ω` (affine-in-ω) | 1.306e+00 | 3.206e-02 | 2.777e-01 | **2.1e-15** |

**max over the physical modes `m = 0…4`:**

| convention | score | rank |
|---|---|---|
| **(2) arc / ω-mid LANDED** | **1.994e-02** | **1st** |
| (1) chord | 2.185e-02 | 2nd |
| (C) `τ ≡ ½` | 3.918e-02 | 3rd |
| (1c) chord + absorber | 5.636e-02 | 4th |

⚠ **Configuration, honestly scoped.** Read the ladder at `n_φ ≥ 16` only. At
`n_φ = 8` the level has `M = 4` cells, so `cos 4ω` sampled at the 4 ω-midpoints
**aliases** — `[M]` `τ ≡ ½` scores `2.220e-15` on it, which is a degeneracy of
the sampling, not accuracy. The ordering `arc < chord < τ≡½ < absorber` holds at
`n_φ = 16, 32, 64` (max over `m = 0…4`, in that convention order:
`4.430e-1 / 5.161e-1 / 8.284e-1 / 8.313e-1` · `8.584e-2 / 9.940e-2 / 1.648e-1 /
2.083e-1` · `1.994e-2 / 2.185e-2 / 3.918e-2 / 5.636e-2`). The only place chord
edges arc is `cos 2` at `n_φ = 32`, by 1.5 % — both are P2-in-η, so they are
expected to be close.

⟹ **The ranking is EXACTLY REVERSED.** Under the memo's `{η, ξ}` basis
`τ ≡ ½` is best and the LANDED arc worst-but-one; under the physically
realisable basis the LANDED arc is best and `τ ≡ ½` 2× worse. And the
discrimination is sharpest exactly at the P1 mode, where the arc convention
is **exact** (2e-15) and `τ ≡ ½` is not (2.4e-3).

The memo diagnosed its own circularity one level too shallow. Its §9bis.6
says *"the convention made one basis function exact and was then certified
by the instrument measuring that same basis function"* — true. But the
symmetric charge lands on C6: **`τ ≡ ½` is exact on `ω` (2.1e-15), and C6's
`ξ` leg is the basis function closest to `ω`'s span** (both spread over all
harmonics with an `m⁻²` tail). C6 does not escape the circularity; it swaps
which convention benefits from it.

### 2.4 (b) On the `ξ`-weighting, and on the "no single τ can be exact for both"

* The weighting is **fine** (§2.1) and defensible.
* The `sin p + sin q = sin(p+q)` argument is **sound but vacuous**: it proves
  no single `τ` is exact for both `η` and `ξ`, which only matters if both are
  required. §2.2 shows `ξ` is not required — it is not a P1 mode and is not
  even a linear function on the folded rule. The genuine statement is much
  weaker: *no two-point linear closure is exact beyond one mode*, and the
  design question is **which one mode**. In this geometry the answer the
  physics gives is `cos ω`, and both derived partitions (chord and arc)
  already deliver it.

`VERDICT: C6 BREAKS. Its stated justification smuggles in a mode that this geometry forbids and this quadrature cannot represent: J_φ ≡ 0 by the ξ→−ξ reflection symmetry, so the P1 limit at a level is affine in η ALONE, and on a folded level `quad.mu_y` samples |ξ| (measured Σwξ = +6.703 ≠ 0, min ξ = +0.3125) — an even, infinite-harmonic, sqrt-singular-in-η function. The ξ(e_arc) WEIGHTING is exonerated (ranking stable across 5 weightings), but the BASIS is not: under the physically-allowed family {cos mω, m=0..4} the ranking INVERTS — arc/LANDED 1.994e-02 (best), chord 2.185e-02, τ≡½ 3.918e-02, chord+absorber 5.636e-02 (worst) — and at the P1 mode itself arc is exact (2e-15) while τ≡½ is 2.4e-03. C6 is not a first-order criterion; it is a higher-order criterion that happens to favour the convention exact on ω.`

---

## 2bis. ⭐⭐ §9bis.9 LANDED MID-REVIEW, and it is where C6/C7's real defence lives

The memo grew 721 → 879 lines while this review was in flight (`8db88596`,
"the literature says tau == 1/2 for the cylinder, and our Hebert citation is
false three ways"). §9bis.9 is the strongest support C7 has, so it is reviewed
here rather than left out.

### 2bis.1 Every OCR-based FACT in §9bis.9 is CONFIRMED — I re-read the sidecar

`scratch/literature_ocr/Hebert(2009)Chapter3.md`, my own reading:

| §9bis.9 claim | my check | verdict |
|---|---|---|
| §3.9.3 = cylindrical, §3.9.4 = spherical | headers at lines 2796 / 2946 | ✅ |
| Eq. 3.431 (sphere) is the *diamond* | line 3010: *"the **diamond differencing scheme** is written"* `φ_{n,i} = ½(φ_{n−1/2,i}+φ_{n+1/2,i})` | ✅ |
| Eq. 3.439 (sphere) = `2φ_n − φ_{n−1/2}` | line 3052 verbatim | ✅ |
| Eq. 3.406 (cylinder) is the *diamond* | line 2876: same sentence, azimuthal index `q` | ✅ |
| Eq. 3.414 (cylinder) = `2φ_q − φ_{q−1/2}` | line 2919 verbatim | ✅ |
| Hébert states no τ | no τ in either section | ✅ |

⟹ **the production docstring's "Morel–Montry weighted-DD angular recurrence of
Hébert Eqs. 3.437/3.439" is indeed false three ways, and that is a REAL DEFECT
that must be fixed under every outcome.** §9bis.9 earns that finding outright.

### 2bis.2 ⛔ FINDING 9bis.9-a — but the literature CONCLUSION is refuted by BMC's own text

§9bis.9 concludes: *"`τ ≡ ½` on the cylinder is not a foreign method (P-F's
charge); it is our own predicate on the rule our own primary source
recommends."* Three measured objections, all from the sidecars.

**(i) BMC name `τ = ½` "the diamond scheme" IN THEIR CYLINDER SECTION, and the
paper exists to condemn it.** `[M]`
`scratch/literature_ocr/Bailey-Morel-Chang(2010)….md`:

* line ~479, under the **R-Z / cylindrical** Eq. (53): *"As in the
  spherical-geometry case, `τ_{m,n} = 1` gives the step scheme and
  `τ_{m,n} = ½` gives the **diamond scheme**."*
* line 95 (abstract) and line 801 (conclusions): *"any general weighted diamond
  scheme, **including the step and diamond schemes**, preserves the diffusion
  limit **to leading order**"* — while *"only the weighted diamond relationship
  of Morel and Montry preserves the asymptotic diffusion limit **to first
  order**"*, and *"If a numerical method does not preserve the diffusion limit
  through first order, we conclude that it is **not as accurate**."*

So C7 proposes, for the cylinder, exactly the scheme BMC's paper is written to
classify as first-order-deficient — in BMC's own cylindrical section. "Not a
foreign method" is true; "not the deficient one" is false.

**(ii) BMC's τ for the CYLINDER is barycentric in the RADIAL COSINE.** `[M]`
BMC Eq. (74) (printed p. 160) and Eq. (75) are written in `μ_{m±1/2,n}` — the
radial cell-edge cosines from the Eq. (52) recursion — and **ω never appears**.
Independently confirmed by the `literature-researcher` sweep
(`scratch/q64_attempt2_lit_check.md`). So the literature does not merely fail to
endorse the ω chart; it writes the cylinder τ in the cosine.

**(iii) The predicate the sources actually STATE is chart-FREE, and it selects
P2-in-η.** `[M]` BMC under Eq. (43):

> *"Equation (43) implies that Eq. (15) will **exactly relate the cell-edge and
> cell-center fluxes when the angular flux assumes the linear form defined by
> Eq. (1)**."*

and in the conclusions (line ~810): *"Preservation of the Galerkin
approximation requires that a scheme be able to accurately represent an angular
flux solution with a **linear dependence in the direction cosines**."*

`[M]` BMC **Eq. (1)** is `ψ_m = φ/4π + (3/4π) J_r μ_m` — **the radial current
alone, one cosine**. That is a chart-free exactness property, and §2.3 measures
which convention delivers it: arc/chord `2e-15`, `τ ≡ ½` `2.4e-03`.

⟹ **`lessons-L48` ("take the PREDICATE, not the recipe") applied correctly
selects the LANDED arc convention.** §9bis.9 invokes L48 and then substitutes a
*different* predicate — "barycentric in the variable the cells are equal in" —
which no source states. BMC's recipe (Eq. 52 cumulative-weight edges) is indeed
inadmissible on our rule (`[M]` P-H, τ outside `[0,1]`, worse with refinement);
their **predicate** is admissible, and the arc partition realises it.

### 2bis.3 ⛔ FINDING 9bis.9-b — the Hébert appeal is ARM-ASYMMETRIC, unjustified

`[M]` Hébert uses the SAME diamond in BOTH arms (Eq. 3.406 cylinder = Eq. 3.431
sphere, same sentence, same scheme). ORPHEUS's sphere **deliberately rejects**
it: `[M]` `morel_montry_tau_per_level` on `gauss_legendre(N)` gives
τ ∈ `[0.3992, 0.6008]` (N=4) … `[0.3897, 0.6103]` (N=64), i.e. `max|τ−½|` grows
`0.1008 → 0.1103` — never `½`.

So "Hébert prescribes `τ≡½` for the cylinder" carries **zero
arm-discriminating information**: he prescribes it for the sphere too, and the
tree overrules him there on BMC/Lathrop's authority. Adopting him on one arm and
not the other needs an independent reason.

§9bis.9's reason is that the ω-midpoint **nodes** satisfy the diffusion moment
condition exactly (`[M]` 2nd-moment rel err `≤2.1e-16`) whereas cosine-midpoint
nodes give Lathrop's `⅔(1−1/N²)`. That measurement is real and useful — **but it
is a NODE property, and probe D already measured P1 as τ-BLIND.** It establishes
that the cylinder's *quadrature* is fine; it says nothing about the *closure
weight*. Offering it as the reason to change τ is level conflation — the
`vv-principles` #9 shape, with "the nodes are already right" standing in for
"the closure is right".

### 2bis.4 Facts from the literature sweep that support the memo, recorded fairly

* **C1 is confirmed VERBATIM** by Hébert **Eq. (3.157)** (his letters swapped:
  `μ` radial, `η` azimuthal, `ξ` axial), with the non-conservative twin at
  (3.155) and the enabling identity `∂η/∂ω = μ` printed — i.e. exactly my §1.1
  `∂ξ/∂ω = η`. Corroborated by Bell & Glasstone p. 58 + Table 1.2 (the
  "infinite cylinder, axial symmetry" row) and BMC Eq. (48).
* **The face coefficient is `ξ` at the face and nothing else** — Hébert
  **Eq. (3.393)**: `W_p[η_{q+½}φ_{q+½} − η_{q−½}φ_{q−½}]`. ✅ C1(b).
* ⭐ **The sharpest reconciliation, and it dissolves §3's whole fight: TWO
  DIFFERENT cosines live at every azimuthal face.** The **azimuthal** one
  (`α/W_p`, whose value is fixed by the conservation recursion) is the
  *streaming coefficient*; the **radial** one (BMC Eq. 52) is what *τ* is
  barycentric in. C1(b) and BMC Eq. (74) are not in conflict — they name
  different numbers at the same face. This is the correct replacement for §3's
  "recursion-defined vs geometric edge" framing AND for §9bis.2's reinstatement.
* **No source prescribes any clamp** — BMC under Eq. (53): *"the weighting
  factors can take on any value between zero and one"*; Eq. (47) + Fig. 1 give
  `τ₁ = 1 − 1/√3 ≈ 0.42265`. ✅ the absorber retirement stands on the
  literature, independent of the partition choice.
* ⚠ **`τ ≥ ½ ⟺ non-amplifying` is in NO source** — an ORPHEUS observation.
  §9bis.9 already records this honestly; §5 below shows it is also mis-scoped.
* ⛔ **A published typo that would invert C1**: BMC printed p. 156 defines
  `η = sinθ cos ω` where it must be `sin ω`. Taken literally it makes the
  azimuthal-face coefficient the radial cosine — from the very paper ORPHEUS
  cites. Do not "correct" the code to match it.
* ⛔ **The primary source for ORPHEUS's cylinder cell-edge construction has
  never been read**: Alcouffe & O'Dell (Hébert ref. [36]) is unresolved after
  7 queries / 4 databases, and Morel & Montry (1984) TTSP **13**(5) 615–633,
  DOI `10.1080/00411458408211661`, is not local. Acquiring them is the user's
  call, and M&M is where a clamp/positivity argument could legitimately live.
* ⚠ **The geometric-edge closed-form τ FLIPS SIGN with march orientation**:
  `½ − ½cot(ω)tan(Δω/4)` if the index runs with ω (Hébert), `½ + ½…` if it runs
  with the radial cosine (BMC/ORPHEUS). `[M]` both agree to `8.9e-16` as a SET,
  reversed in order — **ORPHEUS's `+` is correct, and a level-symmetric fixture
  cannot see the flip** (the `lessons` B5 shape). Whichever way Q5.6.4 lands,
  that is an ungated convention.

`VERDICT: §9bis.9's OCR facts are all CONFIRMED by independent re-reading, and its production-citation finding (Hébert 3.437/3.439 cited for a both-arms module; they are unweighted; Hébert states no τ) is a genuine defect that must be fixed under EVERY outcome. Its literature CONCLUSION is REFUTED: (i) BMC name τ=½ "the diamond scheme" in their own cylinder section (Eq. 53) and their thesis is that the diamond preserves the diffusion limit only to LEADING order; (ii) BMC's cylinder τ (Eq. 74/75) is barycentric in the RADIAL COSINE and ω never appears; (iii) the predicate the sources STATE — "exactly relate cell-edge and cell-center fluxes when ψ has the linear form of Eq. (1)", Eq. (1) = φ/4π + 3J_r μ/4π — is chart-FREE and selects P2-in-η. The Hébert appeal is additionally arm-asymmetric (he prescribes the same diamond for the sphere, which the tree rejects: measured sphere τ ∈ [0.3897, 0.6103]), and the reason offered for the asymmetry is a NODE property that probe D already measured as τ-blind.`

---

## 3. Consumers that break under C7's return-type change

Full three-search audit delegated to `explorer`; the in-process swap probes
below are mine. Every direct consumer wants the **COSINE**. There are exactly
5 direct call sites (4 production, 1 test) plus the `__all__` export.

### 3.1 The three risk classes

| class | consumer (`file:line`) | what happens |
|---|---|---|
| **guarded / LOUD** | `pole_angular_closure.py:1014` `morel_montry_tau_per_level` | **RAISES.** `[M]` its own P3 guard fires: `τ_raw[0] = 4.598010630209084 ∉ [0,1]` on `folded_product(4,8)` level 0 (baseline `[0.2599, 0.4588, 0.5412, 0.7401]`). Every cylindrical solve dies at closure construction. |
| **unguarded, GATED** | `angular_differencing.py:301` `contamination_beta` (BMC β₁) | **silently wrong.** `[M]` `[-6.9e-18, 0, 0, -6.9e-18] → [-0.4436, -2.8462, -2.8462, -0.4436]`. Caught by `tests/sn/primitives/test_quadrature.py::TestL0TermVerification::test_contamination_beta_cylindrical` (`[M]` reds at `β_max = 2.85e+00`). |
| **unguarded, UNGATED** | `angular_differencing.py:339` `alpha_defect_beta` (Lathrop β₂) | **silent garbage.** `1 − e²` fed radians: `[M]` `[-0.7416,…] → [+8.8696, 4.8083, 1.8303, −0.1265, −1.0]`; `8.8696 = π² − 1`. **Zero test consumers tree-wide.** |
| **unguarded, UNGATED** | `angular_differencing.py:378` `nu_closure_residual` | **`inf`.** `[M]` `e[-1] == 0.0` *exactly* in ω (`edge_omega[M] = 0.0`), and the body's last line is `nu / float(e[-1])`. Confirmed by my own probe: `1.0/0.0 = inf`, all 4 levels. No raise. **This is the memo's own headline discriminator.** |

⟹ C7's sentence *"consumers that need the radial cosine (BMC's β functionals)
convert there, where the cosine is the functional's own variable"* is a
**Pattern-7 violation** (`coding-elegance`): it moves a convention from the
producer to N consumers, and 2 of the 3 landing sites have **no gate at all**.

### 3.2 Latent twins that would NOT move with the producer

* ⚠ **A production cosine twin, live**: `orpheus/geometry/reduced_operator.py:871`
  `cylindrical_streaming` independently produces `mu_start_per_level = −sinθ`,
  which *is* `edges[p][0]` (`sinθ·cos π`). Its live consumer is
  `pole_angular_closure.py:1414` `_edge_seed_stencil`, which mixes it with
  `quad.mu_x`. Under C7 the producer's `edges[p][0]` becomes `π` while this
  stays `−sinθ`. It does not break — it **desyncs**, i.e. C7 *re-creates* a
  two-representations-of-one-boundary twin, the exact Cardinal-Rule-2 defect
  the carve was undertaken to remove.
* `tests/sn/sweep/curvilinear/test_alpha_closed_form.py:318` hand-rolls
  `omega_half = 0.5*(omega[:-1] + omega[1:])` — the **only** angle-native
  consumer, and it lacks the `ω ∈ {π, 0}` endpoints, so it would silently
  diverge from the producer rather than track it. It is also why that file
  stayed **fully green** under my naive swap (α is τ-independent and its edges
  are its own).

### 3.3 ⛔ FINDING 3-a — the single-source partition producer has NO value gate

`[M]` **no test anywhere asserts the producer's returned values** — no
snapshot, no literal, no monotonicity check, no `[-1,1]` range check, no
`Σ Δμ = 2` check. The only two tests that call it directly
(`test_mms_ordering_blindness.py:302,312`) assert a `pytest.raises` refusal.
The invariants ARE written down — `pole_angular_closure.py:790-795`,
`docs/theory/foundations/structured_geometry.rst:295` — **in prose only**.

### 3.4 ⛔ FINDING 3-b — the FAITHFUL C7 change reds only 2 of 136 rows

Mutation-verified, in-process (`_qamut_c7_faithful.py`: cylinder partition AND
P2 node both in ω ⟹ τ ≡ ½, sphere untouched, the three `angular_differencing`
consumers given C7's promised conversion-at-the-consumer fix). Canonical
`python -O -m pytest`, 7 τ-implicated files:

```
BASELINE                136 passed in 11.89s
FAITHFUL C7             134 passed, 2 failed in 12.33s
```

Both reds are the same value pin:
`tests/sn/sweep/curvilinear/test_tau_producer_equivalence.py::test_cyl_tau_equals_the_ANALYTIC_closed_form_not_the_chord_convention[8]` and `[16]`.

**Everything else is blind**, and these are the gates a reader would assume
cover it:

| gate | why it stays green |
|---|---|
| `test_tau_arc_wellposedness::test_the_folded_tau_is_bounded_with_the_reversal_identity` | asserts `0.25 ≤ τ ≤ 0.75`; `½` is inside |
| `test_quadrature::TestL0TermVerification::test_per_ordinate_flat_flux_consistency[CYLINDRICAL]` | flat flux is τ-blind (any τ reproduces a constant) |
| `test_contamination_beta_cylindrical` | green once the consumer converts back — it is partition-loaded, not τ-loaded |
| `test_alpha_closed_form` (all rows) | α is τ-independent, and its edges are hand-rolled |
| `test_mms_ordering_blindness` ladder | τ-blind by construction |
| `test_march_start_structure` | integer facts off `ξ=0` / `η₀=η₁`, never consults τ |

⟹ **a first-order change to a production angular closure is caught by exactly
ONE gate in the SN suite**, and that gate's reference is the very closed form
being retired. Its replacement under C7 would be "τ ≡ ½", which on **every**
shipped cylinder rule is a constant — so it could not distinguish a derived ½
from a hardcoded ½ (see §4.3).

### 3.5 Docs — one equation of record goes present-tense-false

`docs/theory/foundations/structured_geometry.rst:291-302` is the doctrine home
(`:label: angular-cell-partition`) and is written with `\mu_{m+1/2}` on the LHS;
:284-286 says *"the partition is the `(M+1)` cell edges **in the radial
direction cosine**"*; :426-428 derives P3-is-a-theorem in η. Also
`curvilinear_one_group.rst:478-486`, `verification/sn.rst:1257/1424/1829`, and
the docstrings of `reduced_operator.py:113-141`, `angular_differencing.py:6-8`
and `106-119`, plus `pole_angular_closure.py:764-861` itself.

⚠ Confirmed: **`pole_angular_closure` is in NO `automodule`** (51 in the corpus;
the six under `orpheus.sn` are `solver`, `mesh.augmented_mesh`,
`loss_representation`, `operators.streaming`, `operators.boundary`,
`operators.radial_characteristic`). Neither is `geometry.reduced_operator` nor
`derivations.discrete.sn.angular_differencing`. So **no Sphinx severity can see
their Python-domain refs**, and `tools/check_docstring_xrefs.py` — `[M]` today
`949 files, 15252 roles, 11274 decidable, DEAD TARGETS: 0, EXIT=0` — is the only
gate, and it tests TARGETS, not prose truth. **Grep is the whole audit.**

`VERDICT: every direct consumer wants the COSINE. The change is not silent — morel_montry_tau_per_level RAISES (P3 guard, measured τ_raw[0] = 4.598) so all cylindrical production dies loudly. But three unguarded consumers in the derivations module must each be converted by hand, and TWO OF THREE HAVE NO GATE: alpha_defect_beta returns π²−1 garbage, and nu_closure_residual — the memo's own headline discriminator — returns `inf` because edge_omega[M] == 0.0 exactly. C7 also RE-CREATES a partition twin (reduced_operator.py:871's mu_start_per_level = −sinθ vs the producer's π, consumed live by _edge_seed_stencil:1414). And two findings the brief did not ask for: (3-a) the single-source producer has NO value gate at all — its invariants exist only in prose; (3-b) the FAITHFUL C7 change reds only 2 of 136 rows, both the same value pin whose reference is the formula being retired.`

---

## 4. Does the sphere move?

### 4.1 The claim is TRUE; the evidence offered for it is a tautology

`[M]` my faithful-C7 mutation touched only the cylinder branch and **zero
sphere-labelled rows moved** (of the 21 reds under the naive swap and 2 under
the faithful one, none is a sphere row). Consistent with the audit: no test
asserts the producer's returned values, so no sphere gate or snapshot *can*
move.

But probe E1's `max|P2_in_march − production| = 0.000e+00` at N = 4/8/16/32/64
is **not a measurement** — on the sphere the march variable *is* the radial
cosine, so `P2_in_march` and `production` are the *same expression*, and the
sphere branch of the function is not edited. The row is
signature-tautological (`vv-principles` Mode 8, third fires-but-cannot-fail
class): no input can make it fail. The honest evidence is "the sphere branch is
unchanged", which is a code fact, not a probe.

### 4.2 ⛔ FINDING 4-a — "one rule, both arms" is a re-description, not a unification

C7's own stated benefit is that it makes *"P2 carries no geometry, so ONE body
serves both arms"* true. It does not, in the sense that matters. `[M]`:

| arm | edge rule | are the ordinates the cell midpoints? | resulting τ |
|---|---|---|---|
| sphere | **cumulative WEIGHT** in μ (BMC Eq. 12) | no — `[M]` τ ∈ `[0.3992, 0.6008]` (N=4) … `[0.3897, 0.6103]` (N=64) | non-trivial |
| cylinder | **MIDPOINT** in ω | yes, by construction | **≡ ½** |

The surviving asymmetry is the **EDGE RULE**, and it is untouched by C7. Worse,
the asymmetry is load-bearing in both directions: applying the cylinder's rule
to the sphere would give τ ≡ ½ there (contradicting BMC Eq. 12 and the
committed sphere gates), and applying the sphere's rule to the cylinder is the
already-refuted P-H. **C7 moves the chart; it does not remove the fork.** After
C7 the shared body would read its input in a *different physical quantity per
arm* — a bare-`ndarray` return whose UNITS depend on a tag, which is
`coding-elegance` anti-pattern #13 plus a Pattern-7 regression, not a repair.

### 4.3 ⛔ FINDING 4-b — on every shipped cylinder rule, C7's τ is a CONSTANT

`[M]` `folded_product` is the **only** cylinder-admissible factory:
`product`, `gauss_legendre`, `level_symmetric` are all refused by
`assert_carrying_quadrature`; `lebedev` does not build. And every
`folded_product(n_mu, n_φ)` is ω-**equispaced** (`[M]` `ptp(Δω) < 1e-13` at
`n_mu ∈ {2,4}`, `n_φ ∈ {6,8,16}`).

So probe E2's *"½ is a **consequence** on equispaced arcs, never an
assertion"* is true and unfalsifiable-in-production: the only witness the memo
offers for non-triviality is a **hand-built** non-equispaced arc
(`ω = [2.9, 2.2, 1.9, 1.0, 0.35]`) that no factory can emit. Any gate asserting
"τ is derived, not hardcoded" would be green forever on production rules — the
same class as E1. If C7 lands, that gate must be built on the synthetic arc and
must SAY so, or the tree will carry a derived-ness claim with no production
witness.

`VERDICT: the sphere genuinely does not move (mutation-verified: cylinder-only swap, zero sphere rows red; and no test anywhere asserts the producer's values, so no sphere gate or snapshot CAN move). But probe E1's 0.000e+00 is signature-tautological, not a measurement — same variable, same unedited code. Two further findings: (4-a) "one rule, both arms" is a re-description — the surviving asymmetry is the EDGE RULE (cumulative-weight in μ vs midpoint in ω), it is untouched, and it is load-bearing in both directions; after C7 the shared body's return carries different UNITS per arm. (4-b) folded_product is the ONLY cylinder-admissible factory and every instance is ω-equispaced, so C7's cylinder τ is a CONSTANT on 100% of shipped configurations and its "derived, not hardcoded" defence rests on a hand-built arc no factory can emit.`

---

## 5. What C5 costs — and C5 is mis-scoped in BOTH directions

Production runs exactly C5's recurrence
(`pole_angular_closure.py:1101-1105`), so the per-step factor `|(1−τ_m)/τ_m|`
is the right error-propagation quantity. Beyond that, two corrections.

### 5.1 ⛔ FINDING 5-a — the END-TO-END factor is exactly 1.0 for chord, arc AND ½

The memo's B2 reports the **worst PARTIAL product**. `[M]` `qa_c5_amplification.py`:

| `n_φ` | | chord | chord+absorber | arc LANDED | `τ ≡ ½` |
|---|---|---|---|---|---|
| 64 | worst partial | 4.07e+01 | **1.000000** | 9.44e+00 | **1.000000** |
| 64 | **END-TO-END** | **1.000000** | **2.4549e-02** | **1.000000** | **1.000000** |
| 32 | END-TO-END | 1.000000 | 4.9127e-02 | 1.000000 | 1.000000 |
| 16 | END-TO-END | 1.000000 | 9.8491e-02 | 1.000000 | 1.000000 |
| 8 | END-TO-END | 1.000000 | 1.9891e-01 | 1.000000 | 1.000000 |

Reason: both derived partitions are level-antisymmetric,
`τ(π−ω) = 1 − τ(ω)`, so `Π_m (1−τ_m)/τ_m` **telescopes to exactly 1**. The
amplification is a *transient bulge* in the level interior (`[M]` a unit seed
error peaks at `|Δ| = 6.68` at face 8/16 for arc, `20.36` for chord, and
arrives at the last face with `|Δ| = 1.000000` for chord, arc and `τ≡½` alike).

⟹ **`τ ≡ ½` and chord+absorber do NOT tie.** The absorber is the *only*
convention that **DAMPS** — by `40×` at `n_φ = 64`, and the damping grows
`∝ M`. So the mechanism that distinguishes the memo's empirical winner is
**dissipation, not non-amplification** — and `τ ≡ ½` has none of it.

**Consequence for C7's empirical case:** §9bis.8's PENDING row asks whether
`τ≡½` reproduces the `6.593e-02` chord+absorber anchor. The mechanism above
predicts **it will not** — that number was bought with a 40× seed-error
damping that `τ≡½` does not have. C7 cannot inherit the absorber's accuracy by
inheriting its `min τ = ½`.

### 5.2 ⛔ FINDING 5-b — the positivity "cost" is INVERTED: τ≡½ is the SAFEST derived candidate

The destabilising coefficient is `(1−τ)/τ`, largest for the **smallest** τ. So
a convention with `τ < ½` entries is *strictly more* positivity-exposed than
`τ ≡ ½`. `[M]` `folded_product(4,32)` level 0, five strictly-positive angular
profiles fed through the production kernel with a positive seed:

| profile | chord | chord+absorber | arc LANDED | `τ ≡ ½` |
|---|---|---|---|---|
| flat | ok | ok | ok | ok |
| `1 + 0.9 η/sinθ` (P1) | ok `1.04e-1` | ok | ok | ok |
| beam `exp(4cos ω)` | ok `1.87e-2` | ok | ok | ok |
| **shadow `exp(−6cos ω)`** | **NEG ×7, min −230** | NEG ×6, min −23.3 | **NEG ×7, min −77.2** | **NEG ×6, min −24.2** |
| step-ish | ok | ok | ok | ok |

⟹ C5's caveat (*"`τ≡½` … has the diamond scheme's usual positivity exposure.
`τ > ½` damps. That is a real trade to state"*) states the trade **backwards**
for the candidates actually on the table: `τ ≡ ½` is the *least* exposed of the
four. The absolute exposure is real and shared; `τ → 1` (step) would be the
positivity-optimal and `τ → 0` unbounded.

### 5.3 Is there a fixture where it bites? — and the coverage gap

**No gate covers the angular recurrence's positivity.** `[M]` the two
positivity gates on the curvilinear path are both on the **scalar flux of a
converged solve**, not on `ψ̂`:

* `tests/sn/sweep/curvilinear/test_282_direct_seed_fixed_point.py:265`
  `test_ciii_coarse_sphere_fixed_source_finite_positive` — `assert np.all(flux >= 0.0)`
  on `sol.scalar_flux`, **sphere only**, 16-cell, both inner drivers.
* `tests/sn/sweep/curvilinear/test_w1_clamp_silent_on_flat.py:357`
  `test_unclamped_sphere_flux_strictly_positive` — again **sphere**, again the
  scalar flux.

Grepping the half-angle grid (`psi_half`, `_MMHalfGrid`, `compute_psi_half_per_level`,
`.upstream`) across `tests/` returns **no positivity assertion on `ψ̂`** — the
consumers are algebraic-identity, coupling-structure, and transpose tests. So:

1. There is **no cylindrical** flux-positivity gate at all.
2. There is **no `ψ̂`-positivity gate** on either arm, so the mechanism C5 names
   is entirely ungated today — for the LANDED convention as much as for `τ≡½`.
3. §5.2's `min ψ̂ = −77.2` (arc) vs `−24.2` (`τ≡½`) means a `ψ̂`-positivity gate,
   if written, would grade the LANDED convention **worse**. It is a genuine
   `τ≡½` argument — just not the one C5 makes.

`VERDICT: C5's algebra is right (production runs that exact recurrence) but the claim is MIS-SCOPED IN BOTH DIRECTIONS. (5-a) The END-TO-END product telescopes to exactly 1.0 for chord, arc AND τ≡½ (τ(π−ω) = 1−τ(ω)); the amplification is a transient interior bulge. The ONLY damping convention is chord+absorber (end-to-end 2.45e-02 at n_φ=64, ∝M) — so "non-amplification" is NOT what distinguishes the empirical winner, dissipation is, and τ≡½ has none. This predicts τ≡½ will NOT reproduce §4.6's 6.593e-02 anchor. (5-b) The positivity cost is INVERTED: τ≡½ is the LEAST exposed derived candidate (measured min ψ̂ = −24.2 vs arc −77.2 vs chord −230 on a steep-shadow profile), because (1−τ)/τ grows as τ falls. As for a biting fixture: NO gate covers the angular recurrence's positivity on either arm — the two curvilinear positivity gates are both on the SPHERE's converged SCALAR flux (test_282_direct_seed_fixed_point.py:265, test_w1_clamp_silent_on_flat.py:357), and there is no cylindrical flux-positivity gate at all.`

---
## 6. ⭐⭐ THE MEMO'S §9bis.8 PENDING ROW — FILLED IN, and it decides the matter

`qa_probeC_live.py`: a **live** solve of `cyl_2g_3reg_folded_4x8_dd_n40`
(`max_outer=200, keff_tol=1e-11, flux_tol=1e-10, max_inner=3000,
inner_tol=1e-13` — the snapshot generator's own settings) under each τ
convention, flux shape L∞-normalised per group and compared against the
trajectory-resolvent Variant-α reference, exactly as
`test_phase_e_trajectory_resolvent_flux_shape_crosscheck` does. Gate `1.2e-1`.

`[M]` 2026-08-11:

| convention | `k_eff` | per-group max \|Δφ_norm\| | overall | gate `1.2e-1` |
|---|---|---|---|---|
| (1) chord | 1.2308955887 | `[0.14409, 0.01305]` | **1.4409e-01** | FAIL |
| (1c) chord + absorber **[OLD PROD]** | 1.2302082296 | `[0.06593, 0.01335]` | **6.5934e-02** | PASS |
| (2) arc / ω-mid **[LANDED]** | 1.2310212586 | `[0.12676, 0.01341]` | **1.2676e-01** | FAIL |
| (C) `τ ≡ ½` **[C7]** | 1.2313562779 | `[0.10181, 0.01507]` | **1.0181e-01** | PASS |

**The probe is validated by the memo's own criterion.** It reproduces both
§4.6 anchors to the printed digits — `6.5934e-02` vs `6.593e-02` and
`1.2676e-01` vs `1.268e-01` — and its two `k_eff` values are exactly the
re-baseline pair recorded at
`tests/sn/verification/analytical/test_phase_c_crosscheck.py:214`
(`1.2302082296342958 → 1.2310212585879858`).

### 6.1 ⛔ C7 does NOT recover the old accuracy — §5.1's prediction confirmed

`τ ≡ ½` lands at **`1.0181e-01`**: **1.54× worse** than the old production
(`6.593e-02`), only **1.24× better** than the LANDED arc, and it clears the
`1.2e-1` gate with **15 % margin**. It buys back roughly a quarter of the
regression. The 40× seed-error dissipation the absorber supplied
(§5.1: end-to-end product `2.45e-02`) is not available to `τ ≡ ½`, whose
end-to-end product is exactly `1.000000` — and the number says so.

### 6.2 ⭐⭐ The empirical ordering is rank-ordered by the RECURRENCE, and ANTI-correlated with the CLOSURE

Put the three instruments beside the answer:

| convention | transient bulge `max Π\|(1−τ)/τ\|` | end-to-end product | physical-harmonic defect (§2.3) | **live flux-shape** |
|---|---|---|---|---|
| chord + absorber | **1.00** | **2.45e-02** (damps 40×) | 5.636e-02 (**worst**) | **6.59e-02** (best) |
| `τ ≡ ½` | **1.00** | 1.000000 | 3.918e-02 | 1.02e-01 |
| arc LANDED | 9.44 | 1.000000 | **1.994e-02** (**best**) | 1.27e-01 |
| chord | 40.7 | 1.000000 | 2.185e-02 | 1.44e-01 (worst) |

* **Rank correlation with the recurrence's amplification: PERFECT.** The order
  `{absorber, ½} < arc < chord` is exactly the bulge order `{1.00, 1.00} <
  9.44 < 40.7`, with the absorber's extra end-to-end dissipation breaking the
  tie against `τ ≡ ½`.
* **Rank correlation with closure accuracy: INVERTED at the top.** The
  convention that is *best* on every physical harmonic (arc, `1.994e-02`) is
  *third* empirically; the convention that is *worst* (absorber, `5.636e-02`)
  is *first*.

⟹ **The `1.2e-1` metric on this fixture measures the angular recurrence's
error amplification, not the closure's truncation accuracy.** That is the
quantitative confirmation of §5.1, and it reframes the whole decision:

1. Nothing here grades a *chart*. Both C6's and C7's framing are orthogonal
   to what the number responds to.
2. The trade is **structural and unavoidable on this rule**: P2-in-η on the arc
   partition FORCES `τ = ½(1 + cot ω · tan(Δω/4))`, which dips to `¼` at the
   inward end of every level, which IS the bulge. You cannot have BMC's
   exactness property and a non-amplifying recurrence simultaneously.
3. So C7 is a legitimate *engineering* choice — **"prefer a non-amplifying
   recurrence over a first-order-accurate closure"** — and it must be argued
   and ratified on exactly that ground, with `1.0181e-01` (not `6.593e-02`) as
   its number. It must NOT be argued as "the chart was wrong" (§2, refuted) or
   as "the literature prescribes ½" (§2bis, refuted — BMC name it the diamond
   and prove it leading-order-only).
4. ⭐ **And the third option is still unexamined and is the only one that could
   dominate:** if the bulge's payload is the starting-direction seed error,
   then improving the seed removes the arc convention's penalty **while keeping
   BMC's exactness**. `[M]` a unit seed error peaks at `|Δ| = 6.68` at face
   8/16 under the arc convention and arrives at the level end at exactly
   `1.000000` — so the interior faces, which is where every cell's balance
   reads `α_{m±1/2}ψ̂_{m±1/2}`, carry up to 6.7× the seed error. **The
   experiment: hold the LANDED convention and vary only
   `MorelMontryAngularSweep.psi_half_seed` (zero / edge-extrapolated /
   Carlson-inward), re-running this exact probe.** If `1.2676e-01` falls toward
   `6.6e-02`, the fix is the seed and Q5.6.4 should be closed with NO τ change.
   If it is seed-insensitive, the diagnosis is refuted and item 3 above is the
   decision.

`VERDICT: the PENDING row is filled and the probe is validated against both of the memo's anchors. τ ≡ ½ measures 1.0181e-01 — it passes the 1.2e-1 gate but recovers only ~24 % of the regression and is 1.54× worse than the retired chord+absorber. The empirical ordering is in PERFECT rank correlation with the recurrence's error amplification and INVERTED against closure accuracy, so this metric cannot adjudicate a chart. C7 is defensible only as an explicit stability-over-first-order-accuracy trade at 1.0181e-01, and the seed experiment must be run first because it is the only candidate that could recover the accuracy without giving up BMC's exactness property.`

---

## Does the chain survive?

### Link-by-link

| link | verdict | the measurement that decides it |
|---|---|---|
| **C1** conservative form + `ξ` face coefficient | ✅ **SURVIVES** | independent SymPy derivation from lab-Cartesian components (`residual 0`; two plausible-wrong variants non-zero), plus Hébert **Eq. (3.157)** verbatim and **Eq. (3.393)**'s face term `W_p[η_{q±½}φ_{q±½}]`. ⚠ two *stated-evidence* defects (§1.3, §1.4) |
| **C2** `α = κ·w_gl·ξ(e_arc)` | ✅ **SURVIVES** | `[M]` rel spread `0.0 … 4.0e-15` at `n_φ = 6…64`; `w_gl = w_m/Δω = 0.695709690 = 2×GL(4)`; and κ derived symbolically as the midpoint-vs-exact ratio of `∫_cell η dω`. Already gated (`test_alpha_closed_form.py:132`) |
| **C3** "no real solution ⟹ candidate (3) is ill-posed at every order" | ⛔ **FALSE AS STATED** | `[M]` the failing edge is the one AT `ω = π/2`, which is an edge only for **even `M`**. At odd `M` — `n_φ = 6, 10, 14, 18, 26, 34, 66`, all legal `folded_product` rules — **0 of `M+1`** edges fail (`max κ·sin ω_arc = 0.9069 / 0.9669 / 0.9832 / 0.9898 / 0.9951 / 0.9972 / 0.9992 < 1`). The four probed orders 8/16/32/64 are all even `M`. The conclusion may still be right on other grounds; **the stated evidence is a parity artefact** |
| **C4** the published criteria are τ-blind | ✅ **SURVIVES — the memo's strongest finding** | corroborated: `diffusion_limit_c` reads only `quad.mu_x`+weights (no edges, no τ); BMC β and Lathrop β consume the EDGES but never τ. ⭐ **Sharpening:** the reason is structural — in BMC's construction **τ is not a free parameter, it is DETERMINED by the partition** (Eq. 43 = barycentric on the Eq. 52 edges), so their β constrains the *partition* and τ follows. "Which τ?" is the wrong question; "which partition, and which exactness condition on it?" is the right one |
| **C5** `τ ≥ ½ ⟺ non-amplifying` | ⚠ **TRUE BUT MIS-SCOPED, both directions** | `[M]` end-to-end product is exactly `1.000000` for chord, arc AND `τ≡½` (`τ(π−ω)=1−τ(ω)` telescopes); only chord+absorber damps (`2.45e-02` at `n_φ=64`, `∝M`). And positivity is INVERTED — `τ≡½` is the least exposed derived candidate |
| **C6** the "honest τ-discriminating instrument" | ⛔ **BREAKS** | the `ξ` leg is not a P1 mode here (`J_φ ≡ 0` by the `ξ→−ξ` reflection; BMC **Eq. (1)** = `φ/4π + 3J_r μ/4π`, one cosine; BMC **Eqs. (61)–(62)** write `Ω = μ e_r + ξ e_z` with **no azimuthal component**), and on a folded rule `quad.mu_y` is `\|ξ\|` (`[M]` `Σwξ = +6.703 ≠ 0`, `min ξ = +0.3125`). Under the physically-allowed basis the ranking **INVERTS** |
| **C7** τ = barycentric in the march variable | ⛔ **DOES NOT FOLLOW** | it rests on C6; and its own literature defence (§9bis.9) is refuted by BMC naming `τ=½` "the diamond scheme" in their **cylinder** section (Eq. 53) and proving the diamond only leading-order-accurate |

### ⛔ The link that breaks, and the measurement that breaks it

**C6.** One table, `folded_product(4,64)` level 0, `max|ψ̂(e) − f(e)|`:

| mode family | chord | chord+absorber | **arc LANDED** | `τ ≡ ½` |
|---|---|---|---|---|
| **the memo's basis** `max(η, ξ)` | 6.637e-01 | 1.631e-02 | 1.415e-01 | **1.226e-03** ← memo's winner |
| **the physical basis** `max_{m=0..4} cos mω` | 2.185e-02 | 5.636e-02 | **1.994e-02** ← winner | 3.918e-02 |

And the inversion is **not** an artefact of the arc convention being exact on
the mode it was built for — it wins on **every** harmonic, by ≈2×:

| harmonic | arc LANDED | `τ ≡ ½` | ratio |
|---|---|---|---|
| `cos 1` (the P1 mode) | **2.1e-15** | 2.412e-03 | ∞ |
| `cos 2` | **4.854e-03** | 9.677e-03 | 2.0× |
| `cos 3` | **1.132e-02** | 2.188e-02 | 1.9× |
| `cos 4` | **1.994e-02** | 3.918e-02 | 2.0× |

`τ ≡ ½`'s own exactness mode is `ω` (`[M]` `2.1e-15`), and `ω` has Fourier
content on **every** harmonic with an `m⁻²` tail — it is not a physical mode
and no source names it. `cos ω` is exactly one harmonic and it is BMC's
Eq. (1).

### ⭐ What the evidence actually says, taken together

The measurements are consistent, and they do **not** say what either attempt
concluded:

1. **The arc/ω-mid partition is right** (§9bis.2, and I concur — via §1.3's
   narrower argument: `α`'s recursion is the midpoint rule *on the arc cells*).
2. **`τ` = barycentric in the RADIAL COSINE on that partition is right** — it
   is BMC's stated exactness property (Eq. 43's own sentence + Eq. (1)), it is
   chart-free, it is the *only* candidate exact on the diffusion-limit mode,
   and it is ≈2× better on every physical harmonic. **That is what the tree
   already ships.**
3. **So the empirical penalty is NOT closure truncation.** The convention that
   wins empirically (chord+absorber, `6.593e-02`) is the WORST on the physical
   basis (`5.636e-02`, last of four) and the only one that **dissipates**
   (end-to-end `2.45e-02`, i.e. 40× seed-error damping at `n_φ=64`). The
   LANDED convention is neutral end-to-end but carries a **9.44× transient
   bulge** that multiplies whatever the starting-direction seed gets wrong.
4. ⟹ **The falsifiable diagnosis: the cylinder's angular error is
   seed-and-transient dominated, not closure-truncation dominated.** The
   absorber was buying accuracy by damping a seed error, which is why retiring
   it made the MMS floor 1.8–3.4× worse (§4.1) *with no partition change at
   all* — a fact §9bis.5 already isolates and then attributes to the wrong
   mechanism.

**The experiment that decides it, and it is cheap:** hold the LANDED
convention fixed and vary ONLY the starting-direction seed
(`MorelMontryAngularSweep.psi_half_seed`: zero vs edge-extrapolated vs Carlson
inward). If the arc-vs-absorber gap collapses as the seed improves, the
diagnosis holds and the fix is in the seed, not in τ. If the gap is
seed-insensitive, the diagnosis is refuted and the transient bulge is not the
carrier — in which case the honest position is that the tree ships the
*better closure* and a *worse number*, which is `feedback_principled_over_bit_identical`
territory and a ratification decision, not a re-pose.

### ⛔ What must be fixed regardless of the outcome

1. **The production Hébert citation** (`pole_angular_closure.py` module
   docstring, lines 10–14) is false three ways — §2bis.1, confirmed
   independently. This is a Cardinal-Rule-3 defect and is outcome-independent.
2. **The single-source partition producer has no value gate** (§3.3) — its
   invariants (`[-1,1]`, `Σ Δμ = 2`, monotone, endpoints `∓sinθ`) exist only in
   prose, at `pole_angular_closure.py:790-795` and
   `structured_geometry.rst:295`.
3. **The cylinder τ has exactly ONE catcher** (§3.4) and its reference is the
   formula under debate. `alpha_defect_beta` and `nu_closure_residual` have
   **no test consumer at all**.
4. **There is no `ψ̂`-positivity gate on either arm, and no cylindrical
   flux-positivity gate at all** (§5.3).
5. **The march-orientation sign of the closed-form τ is ungated** —
   `½ ± ½cot(ω)tan(Δω/4)`; `[M]` the two agree as a SET to `8.9e-16` and differ
   only in order, so a level-symmetric fixture cannot see a flip.
6. **`Alcouffe & O'Dell` — the primary source for ORPHEUS's cylinder cell-edge
   construction — has never been read** (unresolved after 7 queries / 4
   databases), and `Morel & Montry (1984)` TTSP **13**(5) 615–633
   (DOI `10.1080/00411458408211661`) is not local. M&M is the one place a
   clamp/positivity argument could legitimately live. Acquiring them is the
   user's call.
7. ⛔ **Do not "correct" the code to BMC's printed p. 156 coordinate line**
   (`η = sinθ cos ω`, must be `sin ω`) — taken literally it inverts C1, and it
   is in the paper the sphere arm cites.

### Refuted candidates, with their structural reason (per `process-discipline`)

* **"The `ξ(e_arc)` weighting biases C6's ranking"** — rejected: the ranking is
  stable across five weightings (§2.1). The bias is in the **basis**, not the
  weight.
* **"C6 should use `ψ` affine in `ω`"** — rejected: `ω` has Fourier content on
  every harmonic with an `m⁻²` tail and is `arccos(η/sinθ)` in the cosine chart
  (infinite derivative at the level endpoints). It is not a physical mode; it is
  merely the mode `τ≡½` happens to be exact on.
* **"C6 should use `ψ` affine in `ξψ` (the conservative flux)"** — rejected: the
  closure reconstructs `ψ` at the face and the face coefficient `ξ(e)` is
  applied *afterwards*, exactly, by C1's FTC identity. Demanding exactness on
  the product would double-count a factor the scheme already carries exactly.
* **"Candidate (3) (τ on α's own recursion-defined edges)"** — rejected on a
  *different* ground than the memo's: it is well-posed at odd `M` (§C3), but the
  κ displacement is `O(Δω²)` and is a **quadrature** artefact of the midpoint
  rule, so building the partition on it would trade BMC's exactness property for
  a second-order perturbation of the same partition. The literature sweep's
  "two different cosines live at every azimuthal face" is the correct dissolution
  of the whole recursion-vs-geometric framing.
* **"The sphere might move"** — rejected: mutation-verified cylinder-only swap,
  zero sphere rows red, and no test anywhere asserts the producer's values.
* **"`test_phase_e_trajectory_resolvent_flux_shape_crosscheck` is
  structurally blind to a τ change"** — rejected: it reads a frozen snapshot
  (`_load_snapshot_scalar_flux`), *but* the snapshot tracks the convention
  because the carve re-baselines it (`39b46a31` "re-baseline the TWO cylinder
  artifacts the ω-partition moved"). ⚠ The residual hazard is real though: the
  gate can only see a τ change if someone remembers to re-baseline, and the
  anchor it is compared against was re-captured BY the carve
  (`coding-standards`' re-baselined-anchor clause).
