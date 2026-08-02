# Issue #326 — what the cylindrical α-recursion IS, and what it assumes about ordering

**Read-only investigation.** Branch `refactor/operator-strategy-layers` @ `b0a003b4`.
`git status --short` at open and close: only `.claude/` memory + skill edits are dirty;
**no production file is mid-edit**, so every `file:line` below is stable at HEAD.

Claims marked `[M]` are **measured** this session with `.venv/bin/python`.
Claims marked `[D]` are **doc/docstring assertions** — reported as CLAIMS and separately
adjudicated against the implementation (`vv-principles` L33). Literature is paraphrased
with page cites; the OCR sidecar is the search surface, the scan is the SSOT — the two
load-bearing literature pages were **spot-verified against the rendered scan**.

---

## HEADLINE — the doc and the code agree with each other and BOTH disagree with the literature

### F1. The α's have a CLOSED FORM. The user's hypothesis is CONFIRMED, and it is the literature's own *definition* of α.

Hébert (2009) §3.9.3 does not merely state the recursion — it derives it, and in doing so
**names what α is**. Hébert Eq. (3.399), verbatim (short quote, p. 72 of the PDF = printed
p. 138):

> Defining $\alpha_{p,q\pm 1/2} = \mathcal{W}_p\,\eta_{p,q\pm 1/2}$, the recurrence
> relation simplifies to  $\alpha_{p,q+1/2} = \alpha_{p,q-1/2} + \mathcal{W}_{p,q}\,\mu_{p,q}$

Here Hébert's `η` is the **tangential (azimuthal) direction cosine** and `q±1/2` are
**real half-angle boundaries in the azimuthal angle ω**. Two equations earlier he obtains it
by exactly the argument the hypothesis proposes — Eq. (3.398) is `η_{q+1/2} = η_{q−1/2} +
(𝒲_{p,q}/𝒲_p) μ_{p,q}`, i.e. the discrete quadrature of `dη/dω = μ` over the angular cell.

**Translating to ORPHEUS's spelling.** ORPHEUS `η ≡ μ_x` (radial) is Hébert's `μ`; ORPHEUS
`ξ ≡ μ_y` (tangential) is Hébert's `η` (`docs/theory/methods/sn/index.rst:308-313`, and
`rules_product.py:17-19`). So:

> **α_{m+1/2} = −W_p · ξ(φ_{m+1/2}) = −W_p · sinθ_p · sin φ_{m+1/2}**

with the integration constant fixed by `α_{1/2} = 0` at the starting boundary `φ = π`
(the most-inward direction, where `ξ = 0`). The hypothesis — "α is `−ξ` sampled at the
half-angle boundaries, up to the starting constant" — is **exactly right**, and it is
Hébert Eq. (3.399) modulo the ORPHEUS sign/normalization crosswalk already tabulated at
`docs/theory/conventions/normalization.rst:150` (`normalization-alpha-crosswalk`).

Sanity check against ORPHEUS's *own* continuous form: `docs/theory/methods/sn/index.rst`
`:label: transport-cylindrical` writes the redistribution as `−(1/r)∂(ξψ)/∂φ`. Integrating
that over one angular cell telescopes *exactly* to `−[ξψ]_{m−1/2}^{m+1/2}`, so the
coefficient multiplying `ψ_{m±1/2}` **is** `−ξ` at the boundary. There is no modelling
freedom here — it is the fundamental theorem of calculus.

**⇒ The ordering question is decidable from the math.** An ordering is correct iff the
resulting `α_{m+1/2}` sequence is `−W_p ξ` evaluated at a set of **contiguous,
monotonically-advancing** azimuthal boundaries. Any other ordering still telescopes (so
flat-flux consistency and closure survive — see F4) but the α's stop being the tangential
cosine at a real boundary, and `ψ_{m±1/2}` stops being a half-angle *face* flux.

### F2. ORPHEUS's cylindrical level spans the FULL circle `[0, 2π)`; the literature's spans a HALF range. Neither the theory page nor any docstring mentions this.

`[M]` `product_mu_phi` samples `φ` at `n_phi` left-endpoints of `[0, 2π)`
(`orpheus/numerics/quadrature/rules_product.py:115`) and puts **all `n_phi`** of them into
one level (`:129-140`). So an ORPHEUS level is **two literature levels glued** — a mirror
pair `(φ, 2π−φ)` for every distinct `η`.

Both primary sources are explicit that a level is a half-range:

* **Hébert §3.9.3, PDF p. 72 (printed p. 138)** — "The integration must be done over two
  octants … Each axial level contains `2ℱ(p)` base points in interval **0 ≤ ω_{p,q} ≤ π**".
  The remaining three-quarters of the sphere is recovered by the **factor 4** in the moment
  sum, Eq. (3.401).
* **Bailey–Morel–Chang (2010) Eq. (52), PDF p. 10 (printed p. 157)** — the level's cell-edge
  *radial* cosines run `μ_{1/2,n} = −√(1−ξ_n²)` → `μ_{M+1/2,n} = +√(1−ξ_n²)`, "where the
  weights in Eq. (52) are normalized on each level to sum to **2√(1−ξ_n²)**".
  That normalization **is** the statement that a ξ-level covers the radial-cosine range
  **exactly once**. A full circle covers it **twice**.

`[M]` The consequence, for `n_mu=2, n_phi=8`, level 0 — the as-built per-level azimuth
sequence produced by `np.argsort(mu_x[level_arr])`:

```
position :   0     1     2     3     4     5     6     7
φ/π      : 1.00  1.25  0.75  1.50  0.50  1.75  0.25  0.00
ξ = μ_y  :  0     −     +     −     +     −     +     0
```

The sort **interleaves the two half-circles**. It is not a monotone march in `φ` on either
branch, and it is not a monotone march in `φ` around the full circle. The cumulative sum is
therefore not a cumulative integral in `ω` — which is precisely the property the α's are
*defined by*.

`[M]` Level-symmetric is worse: an LS "level" contains **both signs of μ_z as well**
(`Quadrature.level_symmetric(sn_order=4)`, level 0: `mu_z ∈ {−0.408, +0.408}`, 16 ordinates
over **4 distinct η** — a 4-fold degeneracy from `(±μ_z) × (±ξ)`). LS6 level 0: 24 ordinates,
18 duplicates. LS8 level 0: 32 ordinates, 24 duplicates.

### F3. The tie-break does NOT move α. It moves **τ** — and with it, WHICH ordinate is step-differenced and which is diamond-differenced.

`[M]` Swapping the two members of a tied mirror pair leaves the α array **bit-identical**
(both partners carry the same `η`, so `α ← α − wη` is permutation-blind within a pair).
`[M]` `τ_raw` *per position* is likewise unchanged. What changes is the **assignment**:

```
n_mu=2, n_phi=8, level 0
  η   : -0.816  -0.577  -0.577  -0.000  +0.000  +0.577  +0.577  +0.816
  τ   :   0.5     1.0     0.5     1.0     0.5     1.0     0.5     1.0
  ords:    4       5       3       6       2       7       1       0     (as built)
  ords:    4       3       5       6       2       7       1       0     (pair 1↔2 swapped)
```

`[M]` **Exactly half the ordinates of every even-`n_phi` product level get `τ = 1` (STEP
differencing in angle) and half get `τ = 0.5` (DIAMOND), strictly alternating** — measured at
`(n_mu,n_phi)` = (2,4), (2,8), (4,16); and at every LS order the pattern is
`[1, 0.5, 0.5, 0.5]` repeating. The tie-break decides **which member of each mirror pair is
step-differenced and which is diamond-differenced**. The two partners *should* carry identical
flux (the 1-D cylinder flux is even in `ξ`), but the M-M recurrence chains them sequentially
and so does not preserve that symmetry — so the assignment is observable.

**Scope note (do not over-claim).** A parallel `numerics-investigator` dispatch measured, on
the same branch, that an *intra-pair swap alone* permutes the per-ordinate angular flux
**exactly** (`max|perm(base) − flip| = 9e-16`), so the scalar flux under an **isotropic**
source moves by only 1–2 ULP; it moves ~20 % under a `ξ`-dependent source, and it breaks the
`ξ → −ξ` mirror symmetry of the angular flux by 2e-3–3.7e-3. Issue #326's 0.6–7.2 % is a
*legacy-vs-`argsort`* comparison — a **broader** reordering than an intra-pair swap. So the
honest statement is: **the `τ` step/diamond split is the mechanism of the ordering's
observability, and the symmetry break is its direct measurable signature; the specific
0.6–7.2 % figure additionally involves the non-tie part of the legacy reorder.**

The `τ_raw ∈ {0, 1}` trichotomy that the codebase documents as a structural fact
(`pole_angular_closure.py:556-568`, `:657-662`) is therefore not an incidental quirk of the
product rule — it is the **direct fingerprint of the doubled covering**. `τ_raw = 0` means the
node sits on the lower edge of its own angular cell; `τ_raw = 1` means it sits on the upper
edge; both are the degenerate signature of assigning *zero-width* angular cells to the
members of a tied pair.

---

## 1. The definition of record

### 1a. In code

| Role | Location | Body |
|---|---|---|
| **cylinder α-recursion** | `orpheus/geometry/reduced_operator.py:778-786` | `α[m+1] = α[m] − w[m]·η[m]`, per level, `α[0] = 0` |
| sphere α-recursion (sibling) | `orpheus/geometry/reduced_operator.py:688-695` | same with `μ`, **plus a closure assert** |
| the `τ` weight | `orpheus/sn/sweep/pole_angular_closure.py:588-611` (raw) / `:664-675` (clamped) | η-midpoint edges, `±sinθ` endpoints, `[½,1]` clamp |
| `c_in` / `c_out` from `α`,`τ` | `orpheus/sn/sweep/pole_angular_closure.py:924-935` | `c_out = α_out/τ`, `c_in = (1−τ)/τ·α_out + α_in` |
| the half-angle recurrence | `orpheus/sn/sweep/pole_angular_closure.py:690-722` | `ψ_{m+1/2} = (ψ_m − (1−τ)ψ_{m−1/2})/τ` |
| **the ordering** | `orpheus/numerics/quadrature/rules_product.py:137-140` | `np.argsort(mu_x[level_arr])` |
| same for level-symmetric | `orpheus/numerics/quadrature/rules_sphere.py:214` | `np.argsort(mu_x[idx])` |

Verbatim, the recursion of record:

```python
# orpheus/geometry/reduced_operator.py:775-786
    # Per-level azimuthal redistribution coefficients
    # Hébert §3.9.4 (cylindrical analog): α_{m+1/2} = α_{m-1/2} − w_m · η_m
    # Ordinates are ordered by increasing η within each level.
    alpha_per_level: list[np.ndarray] = []
    for level_idx in angular_measure.level_indices:
        eta = angular_measure.mu_x[level_idx]
        w = angular_measure.weights[level_idx]
        M = len(level_idx)
        alpha = np.zeros(M + 1)
        for m in range(M):
            alpha[m + 1] = alpha[m] - w[m] * eta[m]
        alpha_per_level.append(alpha)
```

### 1b. In the theory corpus

`docs/theory/methods/sn/curvilinear_one_group.rst:207-257`, section *"The Alpha
Redistribution Coefficients"*:

* `:label: alpha-recursion` (`:214-219`) — `α_{n+1/2} = α_{n−1/2} − w_n μ_n`, `α_{1/2} = 0`.
* `:label: alpha-cylindrical` (`:234-241`) — the per-level `η` form.
* `[D]` `:229-232` — *"On level p, the ordinates are **sorted by increasing η** (radial
  direction cosine), and the recursion uses η instead of μ."* **This is the only statement of
  the ordering requirement anywhere in the corpus.**
* `[D]` `:240-241` — *"Each level's α values form an independent dome from `η = −sinθ` to
  `η = +sinθ`."*
* `[D]` `:245` — *"`α_{n+1/2} ≥ 0` for all n (non-negative dome)."*

**The page states the recursion only. It does NOT state the closed form**, does not name α
as the tangential cosine, does not derive the α's from `dξ/dφ`, and does not mention the
half-range/full-circle question at all. `[M]` A grep of `docs/theory/methods/sn/` and
`docs/theory/conventions/` for *half-range*, *two octants*, *full circle*, *mirror pair*,
*even in ξ* returns **zero** hits in the curvilinear pages (the only `half-range` hits are in
`acceleration.rst`, about Marshak currents — unrelated).

### 1c. In the literature

| Source | Statement | Where |
|---|---|---|
| **Hébert (2009) §3.9.3** | `α_{p,q±1/2} ≡ 𝒲_p η_{p,q±1/2}`; recursion Eq. (3.399); derived from flat-flux, Eq. (3.397)–(3.398); level = `0 ≤ ω ≤ π`, two octants, ×4 in Eq. (3.401); `η_{p,2ℱ(p)+1/2} = 0` (Eq. 3.396); starting direction is `ω = π`, "moving toward the central axis" (Eq. 3.407) | sidecar `Hebert(2009)Chapter3.md:2800-2895`; PDF pp. 71-73 |
| **Bailey–Morel–Chang (2010) Eq. (50)** | `α_{m+1/2,n} = α_{m−1/2,n} − μ_{m,n} w_{m,n}`, **`α_{1/2,n} = α_{M+1/2,n} = 0`** | PDF p. 9 (printed 156) — **visually verified** |
| **BMC (2010) Eq. (52)** | `μ_{1/2,n} = −√(1−ξ_n²)`, `μ_{M+1/2,n} = +√(1−ξ_n²)`, level weights normalized to sum to `2√(1−ξ_n²)` | PDF p. 10 (printed 157) — **visually verified** |
| **BMC (2010) Eqs. (53)/(54)** | the weighted-diamond-in-angle closure and the `ψ_{m+1/2} = ψ_m/τ − ((1−τ)/τ)ψ_{m−1/2}` recursion ORPHEUS runs | same page |

Two transcription notes, both **in the published paper**, not OCR artifacts (confirmed
against the rendered scan): BMC prints `η = sin θ cos ω` in the Fig. 2 caption text where it
must mean `sin θ sin ω` (`μ` is already `sin θ cos ω`), and both Eq. (50) and Eq. (52) print
the recursion with the same index on both sides (`α_{m+1/2,n} = α_{m+1/2,n} − …`) where
`m−1/2` is meant.

### 1d. ⚠ The citation at `rules_product.py:43` is the **wrong paper**, and the equation number does not say what the docstring claims.

```
# orpheus/numerics/quadrature/rules_product.py:40-44
* Bailey, T.S., Adams, M.L., Yang, B., and Zika, M.R. (2009). "A
  piecewise linear discontinuous finite element spatial discretization
  of the transport equation." *Annals of Nuclear Energy* **35**,
  1929-1936. Eq. 50 — α-recursion convention used by the cylindrical
  sweep that consumes this product rule.
```

and again in the function docstring at `:66-68`: *"Per-level indexing lists are sorted by
increasing η = μ_x to match the cylindrical-sweep convention from Bailey et al. (2009)
Eq. 50."*

Three separate defects in that one citation:

1. **Wrong paper.** `docs/theory/methods/sn/curvilinear_one_group.rst:1834-1852` already
   carries a section titled *"Citation correction"* stating that Bailey-Adams-Yang-Zika is
   *"a piecewise-linear FE diffusion paper unrelated to Sn"* and that the intended reference
   is **Bailey, Morel & Chang (2010)**. That correction lists the four modules it fixed —
   `reduced_operator`, `loss_representation`, `diamond`, `pole_angular_closure`.
   `numerics/quadrature/rules_product.py` was created later and **carries the retracted
   citation forward**. This is a live doc↔doc contradiction inside the tree.
2. **Wrong content.** BMC (2010) **Eq. (50)** is `α_{1/2,n} = α_{M+1/2,n} = 0` — the
   two-sided **closure condition**. The recursion is on the same line; the *ordering* is
   nowhere in it.
3. **The claim it is used to justify is not in the source.** Nothing in BMC Eq. (50) says
   "sorted by increasing η". What BMC *does* say (Eq. 52) is that the level's cell-edge
   cosines march monotonically from `−√(1−ξ²)` to `+√(1−ξ²)` **once** — which the full-circle
   ORPHEUS level cannot satisfy.

---

## 2. What the α's are, mathematically — the derivation

Start from ORPHEUS's own continuous form (`docs/theory/methods/sn/index.rst`,
`:label: transport-cylindrical`):

```
(η/r) ∂(rψ)/∂r  −  (1/r) ∂(ξψ)/∂φ  +  Σ_t ψ  =  Q/W
η = sinθ cos φ   (radial),   ξ = sinθ sin φ   (tangential),   μ_z = cos θ
```

Integrate the redistribution term over ordinate `m`'s angular cell `[φ_{m−1/2}, φ_{m+1/2}]`:

```
−∫ ∂(ξψ)/∂φ dφ  =  −[ ξ_{m+1/2} ψ_{m+1/2} − ξ_{m−1/2} ψ_{m−1/2} ]      (exact — FTC)
                =  α_{m+1/2} ψ_{m+1/2} − α_{m−1/2} ψ_{m−1/2}     with  α ≡ −ξ
```

so **α is `−ξ` at the half-angle boundary, by definition, exactly**. The recursion follows
by differencing the closed form and applying the midpoint/left-endpoint rule for `dξ/dφ = η`:

```
α_{m+1/2} − α_{m−1/2} = −(ξ_{m+1/2} − ξ_{m−1/2}) = −∫_{φ_{m−1/2}}^{φ_{m+1/2}} η dφ ≈ −w_m η_m
```

which is `alpha-recursion` verbatim. The starting constant: `α_{1/2} = −ξ(φ_{1/2}) = 0` at
`φ_{1/2} = π`, because `sin π = 0`. **The seed condition `α_{1/2} = 0` is not a convention —
it is the statement that the starting boundary is the plane `ξ = 0`.** Identically, `α = 0`
at `φ = 0`/`2π`, which is why the recursion closes at the *other* `ξ = 0` plane (Hébert
Eq. 3.396, BMC Eq. 50).

**The closed form is EXACT, not merely `O(Δφ)`-consistent.** For the equispaced product rule
the weight `w_m = w_gl·Δφ` is constant within a level, so the recursion's partial sum is a
Dirichlet kernel. With `ω_j = ω_0 + jΔφ` (marched in `ω`):

```
α_k = −w_gl Δφ sinθ Σ_{j<k} cos ω_j
    = −w_gl Δφ sinθ · [sin(kΔφ/2)/sin(Δφ/2)] · cos(ω_0 + (k−1)Δφ/2)
    = −w_gl · κ · [ ξ(ω_{k−1/2}) − ξ(ω_{−1/2}) ],      κ ≡ Δφ / (2 sin(Δφ/2)) → 1
```

(the product-to-sum step gives `A+B = ω_k − Δφ/2 = ω_{k−1/2}` and `A−B = −ω_{−1/2}`). So
`α` **is** `−ξ` at the *midpoint* half-angle boundaries, exactly, with a single discrete
prefactor `κ` that tends to 1 as `Δφ → 0`. Two consequences: (i) the half-angle boundaries
are `ω_j ± Δφ/2` — each node sits at the CENTRE of its angular cell (so the level's first
boundary is `ω_{−1/2} = ω_0 − Δφ/2`, *not* at a node); (ii) the identity is a pointwise
oracle, which is what makes §6's cheap gate possible. Independently derived here and pinned
by the parallel diagnostic `derivations/diagnostics/diag_326_alpha_closed_form.py` (created
2026-08-01 by `numerics-investigator`), whose Result A is this identity and whose Result B is
that the same closed form **fails by O(1)** under the production ascending-`η` ordering.

Since ORPHEUS's `w_m` carries the polar GL weight `w_gl,p` as a factor and the balance
divides by `ΔA/w_m` with the *same* `w_m`, the level-constant scale cancels and the ORPHEUS
array is `α = −w_{gl,p}·ξ` — i.e. Hébert's `𝒲_p η`, matching the crosswalk row exactly.

### 2a. The measured consequence: the interleaved α is the correct α, *sampled at the wrong places*

`[M]` For `n_mu=2, n_phi=8`, level 0, `w = π/4`, `sinθ = 0.816497`:

```
ORPHEUS (full circle, η-sorted), 9 entries:
  [0, 0.641275, 1.094725, 1.548175, 1.548175, 1.548175, 1.094725, 0.641275, 2.2e-16]

Half-range reconstruction (ξ ≤ 0 branch, φ: π → 2π, mirror-paired weights doubled,
the two ξ=0 self-mirror endpoints at single weight), 6 entries:
  [0, 0.641275, 1.548175, 1.548175, 0.641275, 1.1e-16]

ORPHEUS entries at the mirror-PAIR boundaries (indices 0,1,3,5,7,8):
  [0, 0.641275, 1.548175, 1.548175, 0.641275, 2.2e-16]      ← identical
```

So the ORPHEUS α, **read at the mirror-pair boundaries, is exactly the correct half-range
α**. The extra entries (`1.094725` at positions 2 and 6) are *intra-pair* values that
correspond to **no real angular boundary at all** — they are `α` half-way through a
zero-width angular cell. This is the cleanest available statement of the defect: the
recursion is not wrong, the *sampling grid* has spurious points interleaved into it, and
`ψ_{m±1/2}` at those points is a face flux on a face that does not exist.

---

## 3. What the recursion assumes about ordering

**The requirement, stated exactly:** the ordinates must be enumerated so that their
half-angle boundaries `φ_{1/2} < φ_{3/2} < … < φ_{M+1/2}` form a **contiguous monotone
partition of the azimuthal integration interval**, with ordinate `m` interior to cell
`[φ_{m−1/2}, φ_{m+1/2}]`. The integration variable is **`φ` (equivalently Hébert's `ω`),
not `η`.**

* **Does the theory page state the required ordering?** `[D]` Yes, but as *"sorted by
  increasing η"* (`curvilinear_one_group.rst:229-232`) — the **wrong variable**. It is a
  transcription of what the code does, not of what the math requires. It carries no
  derivation and no citation for the ordering specifically.
* **Is "increasing η" equivalent to "increasing φ"?** Only where `η = sinθ cos φ` is
  monotone in `φ` — i.e. on a **half-circle**. `[M]` The ORPHEUS level spans `[0, 2π)`
  (`rules_product.py:115`, and the measured `φ/π` sequence in F2), where the key is 2-to-1
  and the sort interleaves. **On the production quadrature the two are NOT equivalent.**
* **Does the level span `[0, 2π)` or a half range?** `[M]` **`[0, 2π)`.** Read from
  `rules_product.py:115` (`np.linspace(0, 2π, n_phi, endpoint=False)`) and `:129-136` (all
  `n_phi` indices appended to one level), and measured directly: level 0 of `(2, 8)` contains
  `φ/π ∈ {0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75}` with both signs of `ξ`.

**Marching the full circle monotonically in `φ` is a legitimate closed loop, but it does not
produce a non-negative dome.** `α = −W sinθ sin φ` is `≥ 0` on `φ ∈ [π, 2π]` and `≤ 0` on
`φ ∈ [0, π]`; it closes at zero after `2π` of travel (consistent with `Σwη = 0`). So there
are exactly **three** self-consistent designs, and ORPHEUS implements none of them:

| Design | Level content | Ordering | α | Denominator sign |
|---|---|---|---|---|
| **(A) half-range** (Hébert / Alcouffe–O'Dell / BMC) | one half-circle, weights doubled | monotone in `ω` ≡ monotone in `η` | `= −Wξ`, **≥ 0**, closes 0→0 | safe |
| **(B) full circle, φ-monotone** | full circle | monotone in `φ` (wraps) | `= −Wξ`, **sign-changing** | `c_out = α_out/τ < 0` on half the ordinates |
| **(C) full circle, η-sorted** ← **ORPHEUS today** | full circle | interleaved | ≥ 0, but sampled at fictitious boundaries; alternating `τ ∈ {0.5, 1}` | safe, but the scheme is half-step/half-diamond and tie-break-dependent |

Note that (A) is not merely "the literature's convention" — it is the only one of the three in
which `ψ_{m±1/2}` is a face flux on an actual angular cell face AND the `[½,1]` τ-clamp and
the non-negativity stability argument (`curvilinear_one_group.rst:250-251`) are both true
statements rather than accidents.

---

## 4. The closure condition — where it is asserted, and why it cannot adjudicate

**Confirmed: closure is Mode-12 blind to the ordering, exactly as the brief predicts.**
`α_{M+1/2} = α_{1/2} − Σ_m w_m η_m` is a **telescoping sum**; `Σ w η` is a sum over a set and
is invariant under **every** permutation of that set. No closure test can see the ordering.
`[M]` Directly measured on top of that: the α array itself is bit-identical under the
intra-pair swap, so even a *pointwise* α-array comparison is blind.

Where closure is asserted:

| Site | Geometry | Assertion | Blind to ordering? |
|---|---|---|---|
| `orpheus/geometry/reduced_operator.py:693-695` | **sphere only** | `assert abs(alpha[N]) < 1e-12` "GL antisymmetry violated" | yes (telescoping) |
| `orpheus/geometry/reduced_operator.py:778-786` | **cylinder** | **none — there is no closure assert on the cylinder path at all** | — |
| `orpheus/derivations/discrete/sn/balance.py:343-378` (`verify_alpha_closure`, wired as `test_alpha_closure` at `:401`) | both | `Σ w μ = 0` (GL) and `Σ w η = 0` (product, `n_phi ∈ {8,16}`) | yes — **and doubly so**: it rebuilds `η` from a fresh `np.linspace(0, 2π, n_phi)` and never touches `level_indices`, so it cannot see the sort even in principle |
| `tests/sn/sweep/curvilinear/test_pole_angular_closure.py:142-188` — `TestAlphaRecursionIdentities` (3 tests: endpoints-zero, signed-step, dome-non-negative) | **`Quadrature.gauss_legendre` only** | endpoints, step rule, non-negativity | yes — **and the fixture never instantiates the cylindrical `product` rule**, so the cylinder α-dome has no non-negativity or closure coverage at all |

### 4a. Other things that look like they gate the ordering but structurally cannot

* **The flat-flux / streaming-equilibrium identity**
  (`curvilinear_one_group.rst:259-311`, `:label: streaming-equilibrium`). The proof at
  `:280-290` uses only `α_{n+1/2} − α_{n−1/2} = −w_n μ_n` — the **step** of the recursion,
  which holds term-by-term under any permutation. `index.rst:1003-1013` already flags the
  sibling hazard in the corpus's own words: *"the identity `Σ_n w_n(α_{n+1/2} − α_{n−1/2}) =
  0` annihilates per-ordinate redistribution errors that cancel in the sum."* Flat flux is
  blind by construction — and issue #326's own measurement agrees: the 1G/2G **homogeneous**
  rows sit at `1e-14`–`1e-15`.
* **Bit-identity regression tests** — `tests/geometry/test_reduced_operator.py:155-190`
  (`test_alpha_per_level_bit_identical`) compares `reduced.alpha_per_level` against
  `sn_mesh.reduced.alpha_per_level`: two routes to the *same* array from the *same* sorted
  input. This is `lessons L13`'s "frozen-RHS" test: both sides move together, so it is immune.
* **The `-1e-13` dome non-negativity guard** — even if it were extended to the cylinder, `[M]`
  the interleaved full-circle α **is** non-negative (measured: `min α = 0`), so the guard
  would pass on the defective ordering.
* **`beta_first_order_consistent = True`** (`pole_angular_closure.py:839-846`) — the BMC
  Eq. (41) functional `β = Σ_m μ_m[α_{m+1/2}μ_{m+1/2} − α_{m−1/2}μ_{m−1/2}]` is also a
  telescoping angular sum. It is a claim about the *τ formula*, not about the enumeration.

**Net:** every existing gate on the α machinery is either telescoping-invariant, fixture-
restricted to Gauss–Legendre, or a frozen-RHS route-equivalence. Issue #326's acceptance
criterion *"the chosen ordering is adjudicated against a structurally-independent reference"*
is correct — but F1 now gives a **cheaper and stronger** adjudicator than an MMS: the closed
form. See §6.

---

## 5. The starting direction / Carlson seed

**Which ordinate is the starting direction, and how is it selected.**

* The level's starting-direction **edge** is `mu_start_per_level[p] = −sinθ_p`
  (`orpheus/geometry/reduced_operator.py:805-809`) — the `η = −sinθ` level boundary, matching
  BMC Eq. (52) `μ_{1/2,n} = −√(1−ξ_n²)` and Hébert's `ω = π` ("moving toward the central
  axis", Eq. 3.407).
* Whether the level carries **first-class ψ½ state** is the R12a predicate
  `τ_raw[0] ∈ (0,1)` strictly (`pole_angular_closure.py:553-562`,
  `SNMesh.radial_characteristic_levels`). Otherwise the seed is
  `MorelMontryAngularSweep.edge_extrapolated_seed` (`:982-1030`), a 2-point extrapolation
  through `m0 = argsort(mu)[0]` and the first ordinate with `|Δη| > 1e-14` (`:1012-1030`).
* `orpheus/transport/radial_characteristic_field.py:314` re-derives it as
  `most_inward = ords[int(np.argmin(mu[ords]))]`.

**Is the starting direction well-defined under a tied minimum? — `[M]` NO, and the tie is
present for every odd `n_phi`:**

```
n_phi   mu_start=−sinθ    most-inward NODE η   #ordinates at the min   mu_start IS a node?
   3     −0.816496581       −0.408248290              2                    False
   5     −0.816496581       −0.660559610              2                    False
   7     −0.816496581       −0.735638000              2                    False
   9     −0.816496581       −0.767255812              2                    False
  11     −0.816496581       −0.783422732              2                    False
  15     −0.816496581       −0.798654172              2                    False
   8     −0.816496581       −0.816496581              1                    True
  16     −0.816496581       −0.816496581              1                    True
```

For odd `n_phi`, `φ = π` is not sampled, so the most-inward direction is the **mirror pair**
`φ = π ± π/n_phi` — two ordinates at the identical `η = −sinθ cos(π/n_phi)`. `np.argmin`
returns the first minimum, so it resolves to `ords[0]` — but *which global ordinate is at
`ords[0]`* is decided by the same under-determined `argsort` tie-break. The brief's
`argmin == ords[0]` observation is confirmed, and it is a **second, dependent** tie-break, not
an independent one: it inherits the first one's ambiguity rather than adding a new one.
Collapsing `:314` to `ords[0]` (issue acceptance criterion 3) is correct but does **not**
remove the ambiguity — it only makes the single remaining tie-break visible.

> ### ⚠ DOWNGRADE — odd `n_phi` is UNREACHABLE for the cylinder at HEAD
>
> `[M]` `SNMesh(Mesh1D(coord=CYLINDRICAL), Quadrature.product(n_mu=2, n_phi=odd), mats)`
> **raises** before any of the above can bite:
>
> ```
> n_phi= 3: REJECTED BoundaryGeometryMapNotMeasurePreservingError
> n_phi= 4: BUILDS   radial_characteristic_levels = ()
> n_phi= 5: REJECTED BoundaryGeometryMapNotMeasurePreservingError
> n_phi= 6: BUILDS   radial_characteristic_levels = ()
> n_phi= 7: REJECTED …   n_phi= 8: BUILDS ()   n_phi= 9: REJECTED …   n_phi=16: BUILDS ()
> ```
>
> `reflection_index('x') does not preserve the direction-cosine measure w·|μ_x|` (ERR-042
> guard, `orpheus/geometry/boundary/_specular.py:154` via `reflective.py:171`) — an odd-`n_phi`
> equispaced-`φ` grid is not closed under the `x`-mirror, and the cylinder's `r = 0`
> centreline is intrinsically specular, so no outer-BC choice avoids it.
>
> **⇒ The tied-minimum starting-direction case (§5) and the round-off carrying predicate
> (§5a) are LATENT, not live.** They matter only if the guard is ever relaxed or an odd-`n_phi`
> rule reaches the cylinder by another route. `[M]` For every *buildable* (even) `n_phi`,
> `radial_characteristic_levels == ()` — no cylinder level carries ψ½ state, so the whole
> `radial_characteristic*` seed path is **sphere-only in practice**, and the cylinder always
> takes `edge_extrapolated_seed`. This corrects the emphasis of issue #326's point (5): the
> `argmin`/`ords[0]` duplication at `radial_characteristic_field.py:314` is real
> single-source-of-truth debt, but it is not a *live numerical* ambiguity today.

### 5a. `[M]` A second, independent under-determination: the R12a carrying predicate fires on round-off for odd `n_phi`

```
n_phi= 3: tau_raw[0] = 0.9999999999999993   strictly in (0,1)?  True   → level CARRIES ψ½ state
n_phi= 5: tau_raw[0] = 0.9999999999999993   strictly in (0,1)?  True   → level CARRIES ψ½ state
n_phi= 7: tau_raw[0] = 0.9999999999999987   strictly in (0,1)?  True   → level CARRIES ψ½ state
n_phi= 8: tau_raw[0] = 0.0                  strictly in (0,1)?  False  → non-carrying
n_phi=16: tau_raw[0] = 0.0                  strictly in (0,1)?  False  → non-carrying
```

The exact value should be **1.0**: `τ_raw[0] = (η_0 − η_edge_0)/(η_edge_1 − η_edge_0)` with
`η_edge_0 = −sinθ` and `η_edge_1 = ½(η_0 + η_1) = η_0` for the tied pair. It comes out
`1 − 7e-16` because `0.5*(η_0+η_1)` is not bit-identical to `η_0` when the two tied cosines
differ in the last bits. So for odd `n_phi` the level allocates a first-class
`RadialCharacteristicInteriorField` block — a **data-layout** decision, not just a coefficient
— purely on floating-point noise.

Under #325's algebraically-generated nodes the tie becomes exact, `η_edge_1 == η_0`
bit-exactly, `τ_raw[0] == 1.0`, and the level would **flip from carrying to non-carrying**.
`[M]` The LS family already sits at `τ_raw[0] = 1.0` exactly and is non-carrying, so LS is
the free oracle for the post-#325 behaviour. **Given the DOWNGRADE above this is latent, not
live** — but it is worth pinning, because it is a *predicate on a floating-point comparison
that decides a data-layout*, and that pattern will recur wherever exact ties are introduced.

---

## 6. What this means for #326 — the decision it enables

**The ordering does not need to be settled by convention, and does not need an MMS.** F1
gives a pointwise oracle: `α_{m+1/2}` must equal `−W_p ξ` at a real half-angle boundary.
Concretely, three verdicts follow directly from the evidence above:

1. **`np.lexsort((mu_y, mu_x))` (issue #326's proposed fix) makes the order *deterministic*
   but does not make it *correct*.** It picks, for each mirror pair, which partner is
   step-differenced — a choice with no mathematical content, because both members occupy the
   same zero-width angular cell. The issue already says this ("necessary… not sufficient");
   this report supplies the reason: the α's it produces are still sampled at boundaries that
   do not exist.
2. **The mathematically determined fix is the half-range level**, matching Hébert §3.9.3
   (`0 ≤ ω ≤ π`, ×4) and BMC Eq. (52) (level weights summing to `2√(1−ξ²)`). Under it the
   ordering by increasing `η` becomes *identical* to the ordering by `ω`, the tie set is
   **empty** (`η` is injective on a half-circle), `α` is `−Wξ ≥ 0` closing 0→0 at both `ξ = 0`
   planes, `τ_raw` stops being pinned to `{0, 1}`, and **the whole #326 defect class becomes
   unspellable** rather than guarded. `[M]` The evidence that this is not a re-derivation
   from scratch: the ORPHEUS α *already equals* the half-range α at the pair boundaries (§2a).
3. **A necessary regression gate exists and is cheap**: assert `α_{m+1/2} == −W_p ξ_{m+1/2}`
   pointwise (or, weaker but sufficient to bite, that the per-level `τ` is not the
   alternating `{0.5, 1}` pattern). Both are ordering-sensitive by construction — unlike every
   gate inventoried in §4. Additionally: `[M]` #326's criterion "invariant under
   `np.argsort(kind=)`" is **necessary but weak** — it would pass on ordering (C) with a
   stable sort, which is still the wrong grid.

### Open items this investigation surfaced that are NOT in #326's acceptance list

* `orpheus/numerics/quadrature/rules_product.py:40-44` and `:66-68` carry the **retracted**
  Bailey-Adams-Yang-Zika 2009 citation that `curvilinear_one_group.rst:1834-1852` explicitly
  corrects, and attribute to "Eq. 50" a claim (the η-ordering convention) that Eq. 50 does
  not make. Doc↔doc contradiction, live at HEAD.
* `docs/theory/methods/sn/curvilinear_one_group.rst:229-241` asserts the ordering
  (*"sorted by increasing η"*) and the dome shape (*"from η = −sinθ to η = +sinθ"*) with no
  derivation and no mention that the level spans the full circle. Under the current
  implementation both sentences are true only by the accident that the interleaved sort
  happens to be monotone in `η`.
* `orpheus/geometry/reduced_operator.py:778-786` has **no closure assertion** where its
  spherical sibling at `:688-695` does. Weak (telescoping ⇒ blind to ordering) but the
  asymmetry is unexplained.
* `[M]` The **level-symmetric** family (`rules_sphere.py:214`) is degenerate **4-fold**, not
  2-fold: an LS "level" contains both signs of `μ_z` *and* both signs of `ξ` (LS4 level 0 =
  16 ordinates over 4 distinct `η`; LS8 level 0 = 32 over 8). #326's acceptance item
  "`rules_sphere.py:214` audited for the same defect" is therefore an understatement — it is
  the same defect *squared*, and the `±μ_z` grouping is a separate question (is a `|μ_z|`
  level the intended `ξ`-level of BMC Eq. 51, which indexes signed levels?).
* `[M]` §5a: the R12a carrying predicate is round-off-decided for odd `n_phi` and would flip
  under #325 — **latent** (odd `n_phi` is rejected at SNMesh construction; see the §5
  downgrade box).
* `[M]` **Odd `n_phi` product rules cannot be used for the cylinder at all** (ERR-042 specular
  guard). Not a defect, but it bounds the reachable configuration space for #326's fix and
  should be stated in the issue: `n_phi` is effectively even-only, so every production
  cylinder level has an exact 2-fold `ξ`-mirror pairing with two self-mirror `ξ = 0` members
  (`φ = 0` and `φ = π`) — the exact structure the half-range decomposition (§3 design A)
  needs, with the two `ξ = 0` members carrying half weight in each half.

---

## 7. Reconciliation with the parallel #326 work (L-012)

`[M]` HEAD moved from `b0a003b4` → `1659d756` during this dispatch (the boundary campaign
B0–B3.4c; **none of those commits touch the α, quadrature, or closure paths** — every
`file:line` in this report was re-verified against the final tree). Two untracked artifacts
from a parallel `numerics-investigator` dispatch appeared mid-investigation:

* `derivations/diagnostics/diag_326_alpha_closed_form.py` — pins the closed form (Result A)
  and its O(1) failure under the production ordering (Result B). **Independently agrees with
  F1/§2 here**; I re-derived the Dirichlet identity and the `κ` prefactor before reading it.
* `derivations/diagnostics/diag_326_azimuthal_mirror_symmetry.py` — pins the `ψ(ξ) = ψ(−ξ)`
  mirror invariant, measures production breaking it at ~1e-3, and attributes it to the same
  `τ ∈ {1, ½}` split. **This is the strongest available adjudicator** and it is
  structurally independent in the `vv-principles` sense; it also correctly notes that the
  curvilinear MMS **cannot** see the defect (both ansatzes are functions of `η` and `ξ²`,
  i.e. they live inside the symmetric sector). That last point supersedes issue #326's
  proposed acceptance route ("the curvilinear MMS / L1 suite run under each candidate
  ordering") — **the MMS is Mode-7 blind here; the mirror invariant is not.**
* `.claude/agent-memory/explorer/cylindrical_sn_level_order_sensitivity.md` — a memory note
  written by that dispatch into this agent's memory directory. Its key added fact, which this
  report did not have: **the 1-D cylindrical sweep never reads `ξ = μ_y` at all** (grep-verified
  absent from `orpheus/transport/`), so mirror partners are numerically indistinguishable to
  the sweep and differ only through the per-ordinate source and inflow data.

**Where this report goes beyond the parallel work:** the literature adjudication (Hébert
§3.9.3 `0 ≤ ω ≤ π` + BMC Eq. 52's `Σw̄ = 2√(1−ξ²)` normalization ⇒ the level is a *half*
range — §1c/F2, both spot-verified against the scan); the three-design taxonomy (§3); the
`rules_product.py:43` citation defect (§1d); the gate inventory showing every existing α gate
is telescoping-invariant or GL-fixture-restricted (§4); the LS 4-fold degeneracy including the
`±μ_z` grouping (F2); and the odd-`n_phi` unreachability (§5 downgrade box).
