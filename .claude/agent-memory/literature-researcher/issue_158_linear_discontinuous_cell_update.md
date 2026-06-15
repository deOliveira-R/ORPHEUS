---
name: issue-158-linear-discontinuous-cell-update
description: Exact LD (linear-discontinuous) SN spatial cell-update math for ORPHEUS — slab moment system from Larsen-Morel 1989 JCP 83 Eqs (4.1)-(4.3) (canonical, local PDF); diffusion-limit acceptance criterion from LMM-1987 JCP 69 Table I + p.321; curvilinear LD is UNPUBLISHED → derivation path given. update/residual split mapped to CellVisit/UpstreamState/CellResult.
metadata:
  type: project
---

# Linear-Discontinuous (LD) SN cell-update — literature extraction for Issue #158

**Why:** ORPHEUS is adding LD as the first higher-order,
diffusion-limit-consistent occupant of the swappable per-cell
spatial-strategy contract (`orpheus/sn/spatial/cell_update.py`),
sibling to the shipping `DiamondDifference`. LD is chosen over Step
because LMM-1987 proves Step has NO valid diffusion limit while LD
has all four.

**How to apply:** This memo is the implementer's contract. The slab
case is fully published (transcribe Larsen-Morel 1989 Eqs. 4.1-4.3
verbatim). The curvilinear case is NOT published — derive it by
fusing the slab LD moment system with the in-house curvilinear-DD
M-M recipe (`docs/theory/discrete_ordinates.rst` Eq. `dd-solve`).
Code sets `is_linear=True`, `is_positivity_preserving=False`,
`is_affine_scannable=False` (two coupled face/slope moments → not a
single-upstream affine scan).

---

## 0. Canonical sources (all LOCAL in `scratch/literature/`)

| Tag | Citation | DOI (CrossRef-verified) | Role |
|-----|----------|--------------------------|------|
| **LM-1989** | Larsen, E. W. & Morel, J. E. (1989). *Asymptotic solutions of numerical transport problems in optically thick, diffusive regimes II.* **JCP 83(1):212-236.** | `10.1016/0021-9991(89)90229-5` | **THE slab-LD moment system** — §IV "Linear Discontinuous Methods", Eqs. (4.1a-c), (4.2a-c), (4.3a-e). Also proves LD's full diffusion limit (Eqs. 4.4-4.10). |
| **LMM-1987** | Larsen, E. W., Morel, J. E. & Miller, W. F. (1987). *Asymptotic solutions of numerical transport problems in optically thick, diffusive regimes.* **JCP 69(2):283-324.** | `10.1016/0021-9991(87)90170-7` | The diffusion-limit acceptance CRITERION (four limits, Table I p.287) + the "why LD not Step" verdict (p.321 §IX). |
| **Hébert-2009** | Hébert, A. *Applied Reactor Physics* Ch.3 §3.9. | — | In-house DD curvilinear recipe (M-M closure). **Verified: contains NO LD content (0/122 pages).** Curvilinear LD is genuinely unpublished. |
| **LewMil-1984** | Lewis, E. E. & Miller, W. F. *Computational Methods of Neutron Transport* §5.3 / §6. | — | Cited by the code docstrings for LD/positivity; NOT in the local corpus. Use LM-1989 as the primary, LewMil as the cross-reference for the positivity-fixup menu. |

⚠ **The OCR'd fulltext mangles ψ→`$` and the slope ψ̂→`*`.** The
equation IMAGES (read via `Read pages=11,12`) are authoritative; the
transcriptions below are from the images, not the text extract.

---

## 1. SLAB (Cartesian) LD — published, exact

### 1.1 In-cell representation (the unknown set)

LD represents ψ as a **linear function within each cell** via TWO
spatial moments per (group, ordinate): the cell-average flux and the
cell-average *slope*. From LM-1989 Eqs. (4.1a-c), cell `j`, ordinate `m`,
half-width `h_j/2` (h_j = cell width):

```
(4.1a)  ψ_mj     = (1/h_j)   ∫ ψ_m(x) dx                  [cell-average flux]
(4.1b)  ψ̂_mj     = (6/h_j²)  ∫ (x − x_j) ψ_m(x) dx        [cell-average SLOPE moment]
(4.1c)  ψ_m(x) ≈ ψ_mj + (2/h_j)(x − x_j) ψ̂_mj ,   |x − x_j| ≤ h_j/2   [linear reconstruction]
```

Source moments, identical structure (LM-1989 Eqs. 4.2a-c):

```
(4.2a)  Q_j   = (1/h_j)  ∫ Q(x) dx
(4.2b)  Q̂_j   = (6/h_j²) ∫ (x − x_j) Q(x) dx
(4.2c)  Q(x) ≈ Q_j + (2/h_j)(x − x_j) Q̂_j
```

So the per-cell unknown set is **{ψ̄ , ψ̂}** = {cell-average, slope},
NOT two edge values. This is the Lewis-Miller "flat + linear Legendre
moment" representation. The two edge values are *derived* from the
moments via (4.1c) evaluated at x = x_{j±1/2}.

### 1.2 The 2×2 moment system (LM-1989 Eq. 4.3a-b)

I quote the asymptotic-scaled form verbatim (σ_T/ε is leading-order
total XS, εσ_Aj is the scattering XS; the **un-scaled physical form**
for ORPHEUS is in §1.4):

```
(4.3a)  (μ_m/h_j)(ψ_{m,j+1/2} − ψ_{m,j−1/2}) + (σ_T/ε) ψ_mj
            = (σ_T/ε − ε σ_Aj)·(1/2) Σ_n ψ_nj w_n + (ε Q_j / 2)

(4.3b)  (μ_m/(θ h_j))(ψ_{m,j+1/2} + ψ_{m,j−1/2} − 2 ψ_mj) + (σ_T/ε) ψ̂_mj
            = (σ_T/ε − ε σ_Aj)·(1/2) Σ_n ψ̂_nj w_n + (ε Q̂_j / 2)
```

- **(4.3a)** = cell-balance: ∫ over cell `j` of the SN equation
  μ dψ/dx + σ_T ψ = source. The streaming difference
  (ψ_{j+1/2} − ψ_{j−1/2})/h_j is the EXACT integral of dψ/dx (the
  fundamental theorem of calculus over the cell) — this is the
  zeroth (flat) moment equation, identical in form to the DD balance.
- **(4.3b)** = first (slope) moment: multiply the SN equation by
  (x − x_j), integrate over cell `j`. The combination
  (ψ_{j+1/2} + ψ_{j−1/2} − 2ψ_mj)/h_j is the discrete first moment of
  the streaming term. **θ is a free parameter** (LM-1989 sets θ = 1/3
  for the exact-within-SN LD; θ free is used in their asymptotic
  analysis, Ref. [16]). Within the discrete-ordinates approximation,
  (4.3a)+(4.3b) with **θ = 1/3** are EXACT.

### 1.3 The discontinuous (upwind) closure — Eq. (4.3c)

**THE key LD feature.** The closure relating the cell moments to the
*downstream face value* is upwind/discontinuous (LM-1989 Eq. 4.3c):

```
(4.3c)  ψ_mj ± ψ̂_mj = ψ_{m, j±1/2} ,    μ_m ≷ 0
```

i.e. expanded (this is just (4.1c) evaluated at the OUTFLOW face):

```
  μ_m > 0  (sweep +x):  ψ_mj + ψ̂_mj = ψ_{m,j+1/2}   ← OUTFLOW = right edge
  μ_m < 0  (sweep −x):  ψ_mj − ψ̂_mj = ψ_{m,j−1/2}   ← OUTFLOW = left edge
```

**Which face value feeds downstream:** the OUTFLOW edge value
`ψ̄ + sign(μ)·ψ̂` (reconstruction (4.1c) at the downstream face) is what
the next cell consumes as its INFLOW. The INFLOW edge value (upstream
face) is taken as the upwind neighbour's outflow — i.e. the face flux
is **single-valued and = the upwind cell's reconstructed outflow**.
This is the "discontinuous" part: the cell's OWN reconstruction at the
inflow face (ψ̄ − sign(μ)·ψ̂) is generally NOT equal to the inflow value
fed in; LD does not enforce continuity, only upwinding.

Boundary closures (LM-1989 Eqs. 4.3d-e):
```
(4.3d)  ψ_{m,1/2}   = f_m ,   μ_m > 0   (left-incident)
(4.3e)  ψ_{m,J+1/2} = g_m ,   μ_m < 0   (right-incident)
```

### 1.4 Un-scaled physical form for ORPHEUS

Drop the ε-asymptotic bookkeeping. Let `Σ_t` = total XS, and let the
per-ordinate volumetric **source already include scattering + fixed**
(ORPHEUS convention: `update`/`residual` receive `source` = the full
weight-normalised per-ordinate RHS — see contract §3 and the
`cell_update.py:399-408` docstring). Multiply through by h_j (cell
volume in slab, A=1). Two coupled equations per (group, ordinate):

**Average-moment (balance):**
```
  μ_m (ψ_{out} − ψ_{in}) + Σ_t h_j ψ̄ = Q̄ h_j
```
where ψ_out, ψ_in are the OUTFLOW/INFLOW face values (sweep-resolved).

**Slope-moment:**
```
  (μ_m/θ)(ψ_{out} + ψ_{in} − 2ψ̄) + Σ_t h_j ψ̂ = Q̂ h_j ,   θ = 1/3
```

Substitute the upwind closure (4.3c): for μ_m > 0,
`ψ_out = ψ̄ + ψ̂`, `ψ_in` = given inflow; for μ_m < 0,
`ψ_out = ψ̄ − ψ̂`. Writing `|μ|` and using the sweep-resolved
"out/in" convention (ORPHEUS pre-resolves sign — strategy never sees
sign(μ)), the **outflow reconstruction is**:

```
  ψ_out = ψ̄ + ψ̂      (slope sign already folded into "downstream" by sweep direction)
```

and the inflow self-reconstruction `ψ̄ − ψ̂` is discarded (replaced by
the upstream cell's `ψ_in`). Substituting `ψ_out = ψ̄ + ψ̂` into the
two moment equations gives the **2×2 linear system in (ψ̄, ψ̂)**:

```
  [ |μ| + Σ_t h_j        |μ|       ] [ ψ̄ ]   [ Q̄ h_j + |μ| ψ_in   ]
  [                                ] [    ] = [                     ]
  [ |μ|/θ           |μ|/θ + Σ_t h_j] [ ψ̂ ]   [ Q̂ h_j + (|μ|/θ) ψ_in]
```

Derivation of the matrix entries (μ>0, ψ_out = ψ̄+ψ̂):
- balance: |μ|(ψ̄+ψ̂ − ψ_in) + Σ_t h_j ψ̄ = Q̄ h_j
  → (|μ|+Σ_t h_j)ψ̄ + |μ| ψ̂ = Q̄ h_j + |μ| ψ_in  ✓ (row 1)
- slope: (|μ|/θ)(ψ̄+ψ̂ + ψ_in − 2ψ̄) + Σ_t h_j ψ̂ = Q̂ h_j
  → (|μ|/θ)(ψ̂ − ψ̄ + ψ_in) + Σ_t h_j ψ̂ = Q̂ h_j
  → −(|μ|/θ)ψ̄ + (|μ|/θ + Σ_t h_j)ψ̂ = Q̂ h_j − (|μ|/θ)ψ_in

⚠ **Sign caution on the slope row** — the slope-moment row sign
depends on the sign convention chosen for ψ̂ and the orientation of
(x − x_j). The matrix above is the μ>0 instance; the μ<0 instance
flips the off-diagonal slope coupling. **Because ORPHEUS pre-resolves
sweep direction (the strategy sees only |μ| and "downstream"), the
implementer must verify the slope-row sign against a hand-worked
2-cell μ>0 / μ<0 pair** — the round-trip `residual(update(...))≈0`
gate (`cell_update.py:446`) is the catch-net. The DEFINITIVE entries
are best regenerated symbolically from (4.3a-c) with SymPy rather
than hand-transcribed.

### 1.5 SOLVE (`update`) vs APPLY (`residual`) split

- **`update` (solve):** form the 2×2 matrix `M` and RHS `b` above
  (per group, per ordinate; vectorise over group), solve `M·[ψ̄,ψ̂]ᵀ = b`
  (closed-form 2×2 inverse — `det = (|μ|+Σh)(|μ|/θ+Σh) ± |μ|²/θ`),
  then reconstruct `ψ_out = ψ̄ + ψ̂`. Return
  `CellResult(cell_average_flux=ψ̄, outgoing_spatial_flux=ψ_out, ...)`.
  **The slope ψ̂ is an in-cell moment, NOT a CellResult field today**
  — see Gaps §5 (CellResult has no slope slot; needs extension or
  side-channel).
- **`residual` (apply):** given probe `(ψ̄, ψ̂)`, evaluate
  `r = M·[ψ̄,ψ̂]ᵀ − b` (the 2-vector per-group residual) and
  reconstruct `ψ_out = ψ̄ + ψ̂`. At the solved `(ψ̄,ψ̂)`, `r = 0` to FP.
  The residual is **2-component** (average + slope) where DD's is
  scalar — this changes the matvec contract shape (Gaps §5).

---

## 2. CURVILINEAR (sphere + cylinder) LD — UNPUBLISHED → DERIVE

**Explicit statement: there is NO published curvilinear-LD SN
cell-update with the M-M angular closure in the local corpus or in
the open literature I searched.** LMM-1987 (p.321) analyses LD only
in SLAB geometry and defers even the slab-LD asymptotic details to
"a future article" (= LM-1989, slab only). LM-1989 itself is
slab-only. Hébert Ch.3 §3.9 is DD-only (verified 0/122 pages). An
OpenAlex search for "linear discontinuous spherical/curvilinear
diffusion limit" returned no matching reactor-physics paper. The
LMM-1987 closing promise to "extend the above analysis to curvilinear
and multidimensional geometries" (p.322) was, to my search, never
published as a curvilinear-LD *cell-update derivation*.

### 2.1 Derivation path (combine slab-LD + curvilinear-DD)

Fuse **TWO local recipes**:

1. **Slab LD moment structure** = LM-1989 Eqs. (4.1)-(4.3): the
   two-moment (ψ̄, ψ̂) representation, the slope-moment equation (4.3b)
   with θ=1/3, the upwind closure (4.3c).
2. **Curvilinear DD cell-balance with M-M angular closure** =
   in-house `docs/theory/discrete_ordinates.rst`:
   - balance Eq. `balance-general` (label `balance-general`, the
     curvilinear balance with the `(ΔA/w)[α ψ]` angular-redistribution
     term),
   - the WDD/M-M angular closure Eqs. `wdd-closure` + `dd-solve`
     (the τ-weighted `ψ_{n,i} = (1−τ)ψ_{n−1/2} + τ ψ_{n+1/2}`),
   - the assembled denominator
     `denom = 2|μ| A_down + (ΔA/w)(α_out/τ) + Σ_t V`.

**The fusion:** LD changes ONLY the SPATIAL closure (DD's
`ψ̄ = ½(ψ_in+ψ_out)` → LD's two-moment `ψ_out = ψ̄ + ψ̂` + slope
equation). The **angular term is orthogonal** — it acts on the
half-angle index (n±1/2), not the spatial face index — so the M-M
angular closure and the τ machinery **carry over UNCHANGED** into the
LD average-moment (balance) equation. Concretely:

**LD-curvilinear AVERAGE-moment (balance)** — take Eq. `balance-general`,
integrate over the cell (it already IS the cell-integrated balance),
keep the M-M angular substitution exactly as `dd-solve` does:
```
  |μ|(A_out ψ_out − A_in ψ_in)
    + (ΔA/w)(α_{n+½} ψ_{n+½} − α_{n−½} ψ_{n−½})
    + Σ_t V ψ̄  =  Q̄ V
```
Substitute the M-M angular closure (`ψ_{n+½} = (ψ̄ − (1−τ)ψ_{n−½})/τ`)
AND the LD spatial closure (`ψ_out = ψ̄ + ψ̂`). This produces the LD
analogue of `dd-solve`: the same `c_out = α_{n+½}/τ`,
`c_in = (1−τ)/τ·α_{n+½} + α_{n−½}` M-M constants, but with the
streaming-face term carrying `(ψ̄ + ψ̂)` instead of DD's `2ψ̄`.

**LD-curvilinear SLOPE-moment** — multiply the curvilinear SN
equation by the spatial moment weight (x−x_j → for sphere/cyl the
radial moment (r − r_j)) and integrate. **This is the unpublished,
load-bearing derivation step.** The curvature term
`(1−μ²)/r ∂ψ/∂μ` (sphere) / `−(1/r)∂(ξψ)/∂φ` (cyl) acquires a SLOPE
contribution under the (r−r_j) weighting that has NO published form.
Two viable routes:
   - **Route A (recommended): SymPy.** Set up the curvilinear SN
     balance + (r−r_j)-weighted moment symbolically; substitute the
     M-M angular closure and LD spatial closure; collect into a 2×2
     (ψ̄, ψ̂) system. Verify it reduces to (i) slab LD (4.3) when
     A_in=A_out=1, ΔA=0, α=0, τ=1, and (ii) curvilinear DD `dd-solve`
     when ψ̂≡0. These two reductions are the correctness oracle.
   - **Route B: drop the slope-moment curvature coupling** (set the
     (r−r_j)-weighted curvature integral to its DD value, i.e. apply
     the slope only to the spatial streaming and collision terms). This
     is an APPROXIMATION — cheaper, likely adequate for a first cut,
     but NOT diffusion-limit-proven in curvilinear geometry.

### 2.2 How τ enters

τ enters the LD AVERAGE-moment denominator **exactly as it does in
DD** — via the M-M constants `c_out = α_out/τ` (denominator) and
`c_in = (1−τ)/τ·α_out + α_in` (upstream numerator). The DD curvilinear
denominator is:
```
  denom_DD = 2|μ| A_down + (ΔA/w)(α_out/τ) + Σ_t V
```
The LD AVERAGE-moment denominator differs only in the
spatial-streaming coefficient (the LD closure ψ_out=ψ̄+ψ̂ contributes
`|μ| A_out` to the ψ̄ diagonal and `|μ| A_out` to the ψ̂ off-diagonal,
vs DD's `2|μ| A_down` lumped onto ψ̄). The `(ΔA/w)(α_out/τ)` term and
the τ-clamp policy (cylinder clamps τ∈(½,1]; sphere unclamped — see
theory page lines 937-948) are **inherited unchanged**. So:
```
  LD avg-moment ψ̄-diagonal:  |μ| A_out + (ΔA/w)(α_out/τ) + Σ_t V
  LD avg-moment ψ̂-coupling:  |μ| A_out
  LD slope-moment:           [DERIVE — curvature-weighted, Route A]
```

---

## 3. Diffusion-limit consistency condition (LMM-1987 / LM-1989)

### 3.1 The acceptance criterion (the "right LD")

LMM-1987 §III defines **FOUR asymptotic diffusion limits** =
{thick, intermediate} cell scaling × {cell-edge, cell-average} flux:

- **Two scalings** (LMM-1987 Eq. 3.11, p.294): with cell width
  `Δx_j = ε^l h_j`, `h_j = O(1)`:
  - `l = 0` → **thick** (cell ~ one scale-length),
  - `l = 1` → **intermediate** (cell ~ one mean-free-path — the
    neutronics regime),
  - `l ≥ 2` → thin (truncation-error regime; all schemes pass).
- **Two unknowns** per cell: cell-average ψ̄ and cell-edge ψ.

A scheme **passes** a limit if, as ε→0, the discretised SN equations
limit to a *legitimate* discretised diffusion equation (one that, as
Δx→0, yields the correct analytic diffusion equation) for that
flux/scaling combination (LMM-1987 p.295).

**Table I (LMM-1987 p.287) verdict:**
| Scheme | Thick edge | Thick avg | Interm. edge | Interm. avg |
|--------|:----------:|:---------:|:------------:|:-----------:|
| Diamond (DD) | yes | maybe* | yes | yes |
| Step | maybe† | maybe† | no | no |
| Lund-Wilson | no | maybe‡ | no | no |
| Castor | no | maybe‡ | no | no |
| "New" | yes | yes | yes | yes |

(* DD thick-average yes only if (σ_a h)_{j+½}=const; † Step yes only
if σ_a=σ=0 and (σh)=const; ‡ Lund-Wilson/Castor yes only if
(σh)=const. — qualifiers from LMM-1987 p.287.)

### 3.2 The "why LD not Step" verdict

LMM-1987 **p.321 §IX (verbatim sense):** schemes with THREE unknowns
per cell — *"the linear discontinuous, linear moments, and linear
characteristic methods... possess **all four diffusion limits**."*
Step (Table I) fails BOTH intermediate limits → produces inaccurate
solutions on optically-coarse diffusive meshes even with 10³ cells
(LMM-1987 p.321). **This is the load-bearing reference for choosing
LD over Step as the first higher-order occupant.**

### 3.3 The LD proof (LM-1989 §IV, the actual limit verification)

LM-1989 carries out the asymptotic expansion on the LD system (4.3):
- ansatz: expand ψ̄, ψ̂, ψ_edge in powers of ε (Eq. 2.3 form);
- **O(ε⁻¹) equations (4.4a-b):** force ψ̄ and ψ̂ to be ISOTROPIC at
  leading order — `ψ_mj^(0) = ½ φ_j^(0)`, `ψ̂_mj^(0) = ½ φ̂_j^(0)`
  (Eqs. 4.5a-b);
- **O(ε⁰) equations (4.6)-(4.10):** the solvability condition (Eq. 4.7,
  `Σ_n μ_n w_n ψ_{n,j+½}^(0) = 0` = zero leading-order current
  divergence) yields a DISCRETISED diffusion equation for the
  cell-average AND a correct edge relation — establishing all four
  limits AND **accurate diffusion boundary conditions** (LM-1989's
  central result: LD meets BOTH the diffusion-limit condition AND the
  boundary-condition accuracy that the "new" scheme fails, p.222).

**Acceptance test for the ORPHEUS LD implementation:** an MMS or
thick-diffusive benchmark (optically thick cells, c→1) must show the
LD cell-average flux converging to the diffusion solution as the cell
optical thickness stays O(1) — the DD scheme will show the
characteristic O(1)-error degradation that LD removes. This is the
L1 verification gate for "this is the *right* LD."

---

## 4. Notation crosswalk → ORPHEUS dataclasses

| LM-1989 / theory symbol | Meaning | ORPHEUS field |
|---|---|---|
| `ψ_mj` (4.1a) | cell-average flux, ordinate m, cell j | `CellResult.cell_average_flux` (ng,) |
| `ψ̂_mj` (4.1b) | cell-average **slope** moment | **NO FIELD — see Gaps §5** (extend `CellResult`, or carry in a side state) |
| `ψ_{m,j±1/2}` (4.3c) | cell-edge (face) flux | inflow=`UpstreamState.spatial_upstream`; outflow=`CellResult.outgoing_spatial_flux` |
| `μ_m` | direction cosine (slab); `η` radial cosine (cyl), `μ` (sph) | `streaming_terms.abs_mu` (sweep pre-resolves sign) |
| `h_j` | slab cell width = volume (A=1) | `streaming_terms.volume` (slab: V=Δx) |
| `Q_j` (4.2a) | cell-average source moment | `source` arg (ng,) — weight-normalised per contract |
| `Q̂_j` (4.2b) | cell **slope** source moment | **NO ARG — see Gaps §5** (slope source must be supplied) |
| `Σ_t` | total XS (un-scaled σ_T) | `total_xs` arg (ng,) |
| `θ` (=1/3) | slope-moment free parameter | LD strategy constant (`theta=1/3`) |
| `A_in, A_out` | cell face areas (curv.) | `streaming_terms.face_area_inner/outer`; downstream=`visit.face_area_downstream` |
| `ΔA = A_out − A_in` | curvature area diff | derived from `face_area_inner/outer`; `streaming_terms.delta_A_over_w` = ΔA/w |
| `α_{n±½}` | angular-redistribution coeffs | `streaming_terms.alpha_in/alpha_out` |
| `τ_n` | Morel-Montry angular weight | `streaming_terms.tau_mm` |
| `V_i` | cell volume (curv.) | `streaming_terms.volume` |
| `c_out = α_out/τ` | M-M out constant | computed inline (see `cell_balance.py:313`) |
| `c_in = (1−τ)/τ·α_out + α_in` | M-M in constant | computed inline (`cell_balance.py:314`) |
| `ψ_{n−½,i}` | upstream angular half-flux (curv.) | `UpstreamState.angular_upstream` (ng,) |
| `ψ_{n+½,i}` | downstream angular half-flux | `CellResult.outgoing_angular_state` (ng,) |

**Sweep-direction note:** ORPHEUS pre-resolves sign(μ) in the
orchestrator (`CellVisit` carries `face_area_downstream`; the
strategy sees only |μ| and an "into/out-of" view — see
`cell_update.py:239-247`). So the LD upwind closure (4.3c)'s
`μ_m ≷ 0` branch collapses to "ψ_out = ψ̄ + ψ̂" (the sign is folded
into which physical face is "downstream"). The implementer must
confirm the SLOPE sign convention survives this folding (see §1.4
caution + Gaps).

---

## 5. Gaps / risks for the implementer

1. **`CellResult`/`UpstreamState` have NO slope slot.** LD's second
   moment ψ̂ is a genuine in-cell unknown that (a) must be returned
   from `update`, (b) must be supplied as a slope SOURCE `Q̂` to both
   `update` and `residual`, and (c) makes `residual` 2-component
   (average + slope) where DD's is scalar. **This is the biggest
   structural item** — the cell-update contract was shaped around
   DD's single-moment shape. Options: extend `CellResult` with an
   optional `in_cell_moments` field; add a `source_slope` arg; or a
   richer probe type for `residual`. Flag to architecture review
   (Cardinal Rule 2) BEFORE coding — this touches the Protocol.

2. **`is_affine_scannable = False` is mandatory.** The contract
   docstring (`cell_update.py:357-362`, `diamond.py:135`) already
   anticipates this: LD "couples two face moments and is therefore
   NOT affine-scannable." LD MUST route through the DAG wavefront
   schedule, not CumprodScan/ScanMarch. Do NOT implement
   `affine_scan_coefficients` for LD.

3. **Curvilinear LD slope-moment is UNPUBLISHED.** §2.1 Route A
   (SymPy derivation with the slab-LD + curvilinear-DD double
   reduction as oracle) is the principled path. Do NOT hand-transcribe
   a curvature-weighted slope term — there is no source to transcribe
   from. Dispatch `method-implementer` with the SymPy derivation as a
   Branch-1 task; the two reduction oracles (→slab-LD, →curv-DD) are
   the L1 cross-checks.

4. **Positivity: LD is NOT positivity-preserving.** Set
   `is_positivity_preserving = False`. LD can produce negative
   cell-average AND negative edge fluxes in thin/source-steep cells
   (Lewis-Miller §5.3/§6 — the standard counter-example). Standard
   fixups: **set-to-zero** (clip negatives, simplest) and **Larsen's
   fixup** (the lumped-LD / "set negative slope to make edge zero"
   variant). **For a first cut, ship WITHOUT fixup** — match the
   `is_positivity_preserving=False` flag, document the limitation, and
   defer the fixup to a follow-up issue. The fixup makes the scheme
   NON-linear (`is_linear=False`) which would break the clean 2×2
   linear-solve `update`/`residual` symmetry; keep the first cut
   linear.

5. **θ convention.** LM-1989 uses θ=1/3 for the SN-exact LD; θ is a
   free parameter in their asymptotic analysis (Ref. [16] in LM-1989).
   Use θ=1/3. Some texts (and the "lumped LD" / LLD variant) use a
   mass-lumped slope equation that changes the off-diagonal — DO NOT
   conflate. The standard (consistent) LD is θ=1/3.

6. **Slope-row sign (§1.4 caution).** The off-diagonal slope coupling
   sign depends on the ψ̂ orientation and the (x−x_j) weight sign under
   ORPHEUS's sweep-pre-resolution. Regenerate the 2×2 entries
   symbolically from (4.3a-c); validate with the round-trip
   `residual(update(...))≈0` gate AND a hand-worked μ>0/μ<0 2-cell
   pair (the sign error is invisible to the round-trip alone if it is
   consistent between update and residual).

7. **Lewis-Miller §6 not in corpus.** The code docstrings cite
   LewMil-1984 §5.3/§6 for LD; that PDF is NOT local. LM-1989 (which
   IS local) is the superior primary source — it gives the exact
   moment system AND the diffusion-limit proof. Cross-reference
   LewMil only for the positivity-fixup menu (item 4); do not block on
   acquiring it.
