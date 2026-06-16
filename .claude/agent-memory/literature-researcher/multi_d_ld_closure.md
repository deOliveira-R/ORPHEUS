---
name: multi-d-ld-closure
description: Multi-dimensional (2-D/3-D Cartesian) LD spatial closure for SN — the per-cell discretization to extend ORPHEUS's 1-D LD to N-D. CRITICAL finding: the Cartesian-cell multi-D analog of LD is the BILINEAR/TRILINEAR DG-P1 (UBLD), NOT the simplex-P1 LD; Adams 2001 proved simplex-LD FAILS the thick diffusion limit on quadrilaterals while bilinear DFE PASSES. Cell unknowns = 2^d corner coeffs (4 in 2-D, 8 in 3-D), 4×4/8×8 dense Galerkin system. Lumping (FLBLD/LLD, Wareing 2001) suppresses (not guarantees) positivity. Sources: Maginot-Ragusa-Morel 2016 NSE (UBLD Eqs.1-12), Adams 2001 NSE (thick-diffusive verdict), Borgers-Larsen-Adams 1992 JCP (2-D LD diffusion-limit), Wareing-McGhee-Morel-Pautz 2001 NSE (3-D tetra LLD).
metadata:
  type: project
---

# Multi-dimensional LD spatial closure for SN — literature extraction (D5b)

**Why:** ORPHEUS must extend its shipping 1-D slab/Cartesian LD
(`orpheus/sn/spatial/linear_discontinuous.py`, Issue #158) to N-D
Cartesian (`cell_kernel_batch`/`residual_kernel_batch` for
`len(s_axes) > 1`). The 1-D code raises `NotImplementedError` when
`len(s_axes) != 1` with the note "a multi-D LD is bilinear (an
independent slope per axis) — deferred (#158 Increment D)". This memo
is the implementer's contract for that increment.

**How to apply:** This memo resolves the central design fork the
1-D code comment only gestures at. **The multi-D analog of LD on a
Cartesian/quadrilateral cell is NOT the simplex-P1 "average + one
slope per axis" object — it is the BILINEAR (2-D) / TRILINEAR (3-D)
discontinuous finite element (UBLD).** This distinction is
load-bearing for the diffusion limit (Item 5) and is the reason the
1-D Schur-scalar machinery does NOT generalize as a simple
per-axis tensor product. Read Item 1 before anything else.

---

## 0. Canonical sources

| Tag | Citation | DOI | Role |
|-----|----------|-----|------|
| **MRM-2016** | Maginot, P. G., Ragusa, J. C. & Morel, J. E. (2016). *Non-negative Methods for Bilinear Discontinuous Differencing of the S_N Equations on Quadrilaterals.* **NSE 185(1):17-42** (also LLNL-JRNL-697017 / ANS MC2015 conf. preprint). | `10.13182/nse16-38` | **The cleanest modern UBLD/FLBLD derivation** — §2 "Moment Equations" Eqs. (1)-(12): the 2-D Cartesian bilinear weak form, the 4 basis functions, the 4×4 system, the lumping verdict. The MC2015 conf. preprint PDF is on OSTI (`osti.gov/pages/servlets/purl/1343843`); read it directly (the contents API blob doesn't parse). |
| **Adams-2001** | Adams, M. L. (2001). *Discontinuous Finite Element Transport Solutions in Thick Diffusive Problems.* **NSE 137(3):298-333.** | `10.13182/nse00-41` | **THE thick-diffusion-limit verdict for multi-D DFE.** Establishes which DFE schemes survive in arbitrary geometry; proves **simplex-LD FAILS the thick-diffusion limit on quadrilaterals, bilinear DFE PASSES** (cited as [3] throughout MRM-2016). 180 citations. The single most important reference for the ORPHEUS `ld-thick-diffusive` tripwire in multi-D. |
| **BLA-1992** | Börgers, C., Larsen, E. W. & Adams, M. L. (1992). *The asymptotic diffusion limit of a linear discontinuous discretization of a two-dimensional linear transport equation.* **JCP 98(2):285-300.** | `10.1016/0021-9991(92)90143-m` | The original 2-D LD asymptotic diffusion-limit analysis (the multi-D successor to LMM-1987/LM-1989). Establishes the leading-order behavior of LD on 2-D rectangular cells. Predates Adams-2001's quadrilateral generalization. |
| **WMMP-2001** | Wareing, T. A., McGhee, J. M., Morel, J. E. & Pautz, S. D. (2001). *Discontinuous Finite Element S_N Methods on Three-Dimensional Unstructured Grids.* **NSE 138(3):256-268.** | `10.13182/nse138-256` | **The standard production 3-D LLD** — lumped linear discontinuous on tetrahedral meshes. The FLBLD/LLD lumping recipe + the diffusion-limit + transport-regime robustness analysis. Cited as the [7] FLBLD reference. |
| **LM-1989** | Larsen, E. W. & Morel, J. E. (1989). *Asymptotic solutions… II.* **JCP 83(1):212-236.** | `10.1016/0021-9991(89)90229-5` | The 1-D slab-LD moment system (Eqs. 4.1-4.3) — ORPHEUS's CURRENT implementation. LOCAL in `scratch/literature/`. The d=1 reduction oracle. See `[[issue-158-linear-discontinuous-cell-update]]`. |
| **LMM-1987** | Larsen, Morel & Miller (1987). *Asymptotic solutions…* **JCP 69(2):283-324.** | `10.1016/0021-9991(87)90170-7` | The diffusion-limit acceptance CRITERION (four limits). SLAB ONLY (the curvilinear/multi-D promise was never published as a cell-update). LOCAL. |
| **LewMil-1984** | Lewis, E. E. & Miller, W. F. *Computational Methods of Neutron Transport.* §5.3 (LD), §6 (DFE, multi-D). | — | Textbook DG-P1/DFE treatment + positivity-fixup menu. NOT in local corpus. |

⚠ **Local corpus has ONLY the 1-D papers** (LM-1989, LMM-1987). The
multi-D papers (Adams-2001, BLA-1992, MRM-2016, WMMP-2001) are NOT in
`scratch/literature/`. MRM-2016's MC2015 conference preprint is the
one openly readable full-text (OSTI). **Suggest the user add
Adams-2001, MRM-2016 (NSE), and WMMP-2001 to Zotero** — these are the
three canonical multi-D LD/LLD anchors and the Zotero MCP server was
down this session (could not check the library or annotations).

---

## 1. The cell unknowns in d dimensions — THE DESIGN FORK

**The naïve generalization ("ψ̄ + one slope per axis" = simplex-P1,
1+d moments) is the WRONG object on a Cartesian cell.** There are TWO
distinct multi-D "linear discontinuous" objects, and they are not
the same scheme:

### (a) Simplex-P1 LD (the literal "1 + d slopes" basis)
Basis = `{1, x, y[, z]}` — a true affine function, `1+d` moments per
cell. This is what you get by triangulating/tetrahedralizing the cell
into simplices and putting a P1 (linear) DFE on each simplex. It is
the natural object on UNSTRUCTURED triangle/tetra meshes (this is
what WMMP-2001 do in 3-D, after subdividing each cell into
simplices). **On a single Cartesian/quadrilateral cell, simplex-LD
FAILS the thick diffusion limit (Adams-2001, [3]).** Do NOT build the
2-D/3-D Cartesian ORPHEUS kernel as a tensor of 1-D slopes if the
target is a quadrilateral — it is not diffusion-limit-consistent.

### (b) Bilinear / trilinear DG-P1 (UBLD) — the RIGHT object on Cartesian cells
Basis = the **tensor product of 1-D linear Lagrange functions** —
`{1, x, y, xy}` in 2-D (bilinear, **Q1**), `{1, x, y, z, xy, xz, yz,
xyz}` in 3-D (trilinear). The natural unknowns are the **2^d corner
(vertex) values** of the cell (4 in 2-D, 8 in 3-D), one per
(group, ordinate). MRM-2016 Eqs. (8a-d) give the 2-D reference-element
basis explicitly:

```
(8a)  B_1(s,t) = (1-s)/2 · (1-t)/2     [vertex (-1,-1)]
(8b)  B_2(s,t) = (s+1)/2 · (1-t)/2     [vertex ( 1,-1)]
(8c)  B_3(s,t) = (s+1)/2 · (t+1)/2     [vertex ( 1, 1)]
(8d)  B_4(s,t) = (1-s)/2 · (t+1)/2     [vertex (-1, 1)]
```

with the trial space (MRM-2016 Eq. 9):
```
(9)   ψ̃_UBLD(s,t) = Σ_{j=1}^{4} ψ_{j,UBLD} B_j(s,t)
```

**On a Cartesian cell, bilinear DG-P1 PRESERVES the thick diffusion
limit (Adams-2001).** The `xy` cross term is exactly what simplex-P1
lacks and is the reason bilinear survives where simplex-LD fails.

### Crosswalk to ORPHEUS's 1-D moment language
The 1-D LD `{ψ̄, ψ̂}` (average + slope, the **P1 Legendre** basis on
`[-1,1]` with `θ=1/3` the L2 normalization) and the corner-value
basis are **two coordinate systems for the same linear function**.
In 1-D: `ψ̄ = ½(ψ_L + ψ_R)`, `ψ̂ = ½(ψ_R − ψ_L)` (Legendre ↔ nodal is
an invertible 2×2). The natural d-dimensional moment basis is the
**tensor-product Legendre** set:
- 2-D: `{1, P₁(x)} ⊗ {1, P₁(y)}` = `{1, x, y, xy}` → moments
  `{ψ̄, ψ̂_x, ψ̂_y, ψ̂_xy}` — **4 moments, NOT 3.** The 4th
  (cross/bilinear) moment `ψ̂_xy` is the moment ORPHEUS's "one slope
  per axis" intuition omits.
- 3-D: `{1,P₁(x)}⊗{1,P₁(y)}⊗{1,P₁(z)}` = `{1,x,y,z,xy,xz,yz,xyz}` →
  **8 moments.**

**So the corrected answer to the brief's Item 1:** the per-cell
unknown count for the diffusion-limit-consistent Cartesian scheme is
**2^d** (4 in 2-D, 8 in 3-D) per (group, ordinate) — NOT `1+d`. The
`1+d` simplex count is correct only for an unstructured-simplex mesh
(WMMP-2001), where it does preserve the limit because the simplex
geometry is different from the quadrilateral. **For ORPHEUS's
Cartesian N-D extension, build the tensor-product (bilinear/trilinear)
basis, 2^d moments.**

ORPHEUS's `CellResult` carries only `cell_average_flux` + one
`outgoing_spatial_flux`; it has no slot for `d` (or `2^d − 1`)
slope/cross moments. This is the same contract-extension gap flagged
for 1-D in `[[issue-158-linear-discontinuous-cell-update]]` §5, now
multiplied: the per-cell unknown vector is `2^d`-long and each
downstream face carries a `2^{d-1}` in-face moment set (see Item 3).

---

## 2. The per-cell linear system (the "bilinear Schur")

The d-dimensional UBLD per-cell system is the Galerkin weak form of
the SN streaming-collision operator, with the upwind flux on the
INCOMING faces. MRM-2016 §2 gives the 2-D Cartesian derivation
verbatim:

**Mono-energetic SN transport equation (MRM-2016 Eq. 1):**
```
(1)   Ω⃗·∇ψ(x,y,Ω⃗) + σ_t(x,y) ψ(x,y,Ω⃗) = S(x,y,Ω⃗)
```

**The i-th spatial moment (Galerkin, basis B_i, cell K) (Eq. 5):**
```
(5)   ∫_K B_i [Ω⃗·∇ψ + σ_t ψ] dx dy = ∫_K B_i S dx dy
```

**Integration by parts → the weak form with the SURFACE (upwind)
term (MRM-2016 Eq. 6):**
```
(6)   (Ω⃗·) ∮_{∂K} n⃗ B_i ψ dℓ  −  ∫_K ψ [Ω⃗·(∇B_i)] dx dy
                                  + σ_t ∫_K B_i ψ dx dy  =  ∫_K B_i S dx dy
```

The three volumetric/surface operators are, in matrix form per cell:
- **Mass matrix** `M_ij = ∫_K B_i B_j` (the `σ_t` collision term).
- **Gradient (stiffness) matrix** `G_ij = ∫_K B_i (Ω⃗·∇B_j)` — the
  volumetric streaming term (the `−∫ ψ Ω⃗·∇B_i` becomes `+G` after the
  IBP sign, coupling all 2^d moments).
- **Surface matrix** `F` (the `∮ n⃗ B_i ψ` face term) — on each face,
  split into OUTFLOW (`Ω⃗·n⃗ > 0`, implicit, couples the cell's own
  unknowns) and INFLOW (`Ω⃗·n⃗ < 0`, **upwind**: the incoming face
  value is the upstream neighbour's outflow trace → goes to RHS).

**The assembled per-cell system (MRM-2016 Eq. 12):**
```
(12)  A_UBLD ψ⃗_UBLD = R⃗_BLD
```
- `A_UBLD` is a **2^d × 2^d non-symmetric DENSE matrix** (4×4 in 2-D,
  8×8 in 3-D). `A = G + F_out + σ_t M` (gradient + outflow-surface +
  collision-mass).
- `ψ⃗_UBLD` = the `2^d`-vector of bilinear/trilinear coefficients.
- `R⃗_BLD` = `2^d`-vector whose i-th element = "the i-th interpolatory
  spatial moment of the volumetric source within the cell + the
  upwinded angular-flux contributions from the cell edges" (MRM-2016
  verbatim, p.5). i.e. `R = M·S⃗_moments + F_in·ψ_in_traces`.

**Mapping to ORPHEUS notation (the "moment-balance" form the brief
asks for):** writing `s_a = |μ_a|/Δ_a` (ORPHEUS's per-axis streaming,
the raw `g`), `Σ_t`, and the tensor-Legendre moment vector
`ψ⃗ = (ψ̄, ψ̂_x, ψ̂_y, ψ̂_xy)` in 2-D:

- **0th moment (i ↔ B≡1)** = the cell balance:
  `Σ_a s_a·(face-flux jump on axis a) + Σ_t·ψ̄ = S̄` — the direct
  multi-axis generalization of the 1-D balance row
  `|μ|(ψ_out − ψ_in) + Σ_t h ψ̄ = Q̄ h`.
- **1st moments (i ↔ x, y)** = the per-axis slope rows, each
  structurally the 1-D slope row `(|μ|/θ)(ψ_out+ψ_in−2ψ̄) + Σ_t h ψ̂ =
  Q̂ h` ON THAT AXIS, **but additionally coupled to the cross moment
  ψ̂_xy through the transverse face integral** (this coupling is the
  multi-D content absent in 1-D).
- **cross moment (i ↔ xy)** = the bilinear row — has NO 1-D analog;
  it is closed by the `xy` basis function's gradient and surface
  integrals. This row is what makes the system 4×4 rather than 3×3
  and is the diffusion-limit-load-bearing term.

The `θ=1/3` Legendre normalization survives per-axis (it is the
diagonal of the 1-D mass matrix in the Legendre basis); in the
tensor-product basis the mass matrix is the Kronecker product of the
1-D mass matrices, so the diagonal weights are `θ^{(number of active
axes in that moment)}` — `1` for ψ̄, `θ` for ψ̂_x/ψ̂_y, `θ²` for ψ̂_xy.
**Recommended: do NOT hand-transcribe the 4×4/8×8 entries.** Assemble
the three matrices (M, G, F) as Kronecker products of the 1-D
LD operators (which ORPHEUS already has, verified) and let the
linear-algebra do the coupling — this is the elegant, single-source
route and it makes the d=1 reduction to the existing
`linear_discontinuous.py` a Kronecker-with-empty-axes identity.
SymPy-verify against (i) the 1-D LD per-axis reduction and (ii) the
"exact on a bilinear flux ψ = a + bx + cy + dxy" oracle (the multi-D
analog of the 1-D "exact on linear-in-x" oracle).

---

## 3. The upwind face reconstruction in multi-D

In 1-D the outgoing face is the DG upwind trace `ψ_out = ψ̄ + ψ̂`
(LM-1989 Eq. 4.3c). In multi-D, **each outflow face is a (d−1)-
dimensional object that carries the cell's average PLUS the in-plane
(transverse) slope** evaluated at that face's coordinate:

- The OUTGOING face flux on the downstream face normal to axis `a` is
  the trace of the cell's bilinear/trilinear function restricted to
  that face: `ψ_out^{(a)}(transverse coords) = ψ̃(s_a = +1, other
  coords)`. In 2-D, the downstream-x face carries a 1-D LINEAR profile
  in y: `ψ_out^x(y) = (ψ̄ + ψ̂_x) + (ψ̂_y + ψ̂_xy)·P₁(y)` — i.e. a
  `2^{d-1}`-moment face object (average + transverse slope), NOT a
  single scalar.
- The INCOMING face on axis `a` is set to the upstream neighbour's
  outgoing face trace (upwind / discontinuous — no continuity
  enforced). The inflow contributes to `R⃗_BLD` via the surface matrix
  `F_in` integrated against each test function `B_i` over that face.

**Contract consequence:** ORPHEUS's interface-flux cochain
(`WavefrontFlux` / the moving-frontier face storage) must carry, per
face, a `2^{d-1}`-moment object (2-D: average + 1 transverse slope;
3-D: average + 2 transverse slopes + 1 transverse cross), not a
scalar face value. This is a genuine contract widening over the DD
face (which is a single scalar) and over the current 1-D LD face
(also scalar, because the 1-D downstream face is 0-dimensional). Flag
to architecture review (Cardinal Rule 2) — the face-cochain type is
the load-bearing extension.

---

## 4. Lumping (LLD / FLBLD) — verdict: NEEDED for production robustness, OPTIONAL for first correctness

**The difference (mass-matrix lumping):** UBLD uses the consistent
(dense) mass matrix `M_ij = ∫_K B_i B_j`. **Lumping** replaces `M`
(and, in the "fully lumped" FLBLD, also the surface and gradient
matrices) with a diagonal (row-summed) approximation:
`M^lumped_ii = Σ_j M_ij`, off-diagonals zeroed. MRM-2016 (p.4) and
WMMP-2001 give the recipe.

- **Mass-lumped BLD** (MRM-2016 ref [3]): diagonalizes the collision
  coupling. **Suppresses but does NOT guarantee** non-negative angular
  flux. Cheaper linear solve (the cell system becomes more
  diagonally dominant).
- **FLBLD = Fully Lumped BLD** (Wareing et al., MRM-2016 ref [7] =
  WMMP-2001): lumps the mass, surface, AND gradient matrices. More
  robust than mass-only lumping. **MRM-2016 verbatim: FLBLD is "the
  least accurate method, with significantly more numerical diffusion
  than the Petrov-Galerkin scheme and both fix-ups"** — robustness is
  bought with numerical diffusion.
- **SCB (Simple Corner Balance)** (Adams, MRM-2016 ref [8]): an
  algebraically-equivalent reformulation of FLBLD that gives the same
  zone-average solution; more robust on thick meshes but
  "significantly less accurate than UBLD on optically thin or
  intermediate cells."

**Positivity verdict (answers the brief's `is_positivity_preserving`
question):** **NO form of bilinear DFE — not UBLD, not mass-lumped,
not FLBLD — is STRICTLY positivity-preserving.** MRM-2016's entire
premise is that lumping only *suppresses* negatives. Strict
non-negativity requires a NON-LINEAR scheme (their BCSZ
Petrov-Galerkin, or a flux fix-up). So:
- Keep `is_positivity_preserving = False` for the linear multi-D LD
  (UBLD and FLBLD alike). Lumping does NOT change the flag to True.
- A strictly-positive multi-D LD would be non-linear
  (`is_linear = False`), breaking the clean linear cell solve — defer
  to a follow-up issue exactly as the 1-D positivity fix-up was
  deferred.

**Which to ship:** for the FIRST multi-D correctness cut, **UBLD
(unlumped)** is the right target — it is the most accurate, it is the
object the diffusion-limit proofs are stated for, and it keeps the
clean linear `A ψ = R` solve. Add FLBLD/SCB later as a robustness
option for void/grazing-incidence problems (where UBLD oscillates and
goes negative, MRM-2016 Fig. 6). Do NOT make FLBLD the default —
its numerical diffusion is large.

---

## 5. The diffusion limit — the central multi-D verdict

This is where multi-D departs sharply from the 1-D story, and it is
the reason Item 1's basis choice matters.

- **1-D (current ORPHEUS):** full LD (with slope source) preserves
  all four diffusion limits; Step has none (LMM-1987 Table I,
  LM-1989 §IV). See `[[issue-158-linear-discontinuous-cell-update]]`.
- **Multi-D, SIMPLEX-P1 LD on a QUADRILATERAL: FAILS the thick
  diffusion limit (Adams-2001, [3]).** This is the headline result.
  The literal "average + one slope per axis" object does NOT survive
  the thick-diffusive limit on Cartesian/quadrilateral cells.
- **Multi-D, BILINEAR DFE (UBLD) on a QUADRILATERAL: PRESERVES the
  thick diffusion limit (Adams-2001, [3]; BLA-1992 for the 2-D
  rectangular asymptotic analysis).** The bilinear `xy` cross term is
  what restores the limit. MRM-2016 introduction states this
  explicitly: *"While Adams showed that LD on quadrilateral meshes
  fails in the thick diffusion limit, Adams also demonstrated that
  bilinear discontinuous finite elements can maintain the thick
  diffusion limit on quadrilaterals."*
- **Multi-D, SIMPLEX-P1 LD on a tetrahedral/triangle mesh: PRESERVES
  the limit** (WMMP-2001) — but this is a different mesh topology
  (true simplices), not a Cartesian cell.

**Expected behavior for the ORPHEUS `ld-thick-diffusive` tripwire in
multi-D:** if ORPHEUS builds the N-D Cartesian kernel as the
**tensor-product bilinear/trilinear (UBLD)** object (Item 1b), the
multi-D `ld-thick-diffusive` MMS gate is expected to **PASS** (LD
recovers the diffusion solution on optically-thick Cartesian cells),
provided the slope SOURCE moments are threaded (the same Increment-C
requirement as 1-D — the flat-source cut is O(h²) but not
diffusion-limit-consistent). If instead ORPHEUS builds the naïve
simplex "1+d slope" object on Cartesian cells, the gate is expected to
**FAIL** — and that failure would be *correct physics* (Adams-2001),
not a bug. **The tripwire's expected outcome is therefore a function
of the basis choice — document this in the xfail rationale.**

---

## 6. Notation crosswalk → ORPHEUS

| Source symbol | Meaning | ORPHEUS |
|---|---|---|
| `B_j(s,t)` (MRM Eq. 8) | bilinear Lagrange basis, vertex j | tensor-Legendre moment basis (recommended) — Kronecker of 1-D LD basis |
| `ψ_{j,UBLD}` (Eq. 9) | bilinear coefficient (vertex value) | the `2^d`-moment cell unknown vector; ψ̄ = (B≡1) component |
| `A_UBLD` (Eq. 12) | 2^d×2^d cell matrix = G + F_out + σ_t M | the multi-D analog of `_LDCellTerms`; assemble as Kronecker of 1-D operators |
| `R⃗_BLD` (Eq. 12) | source moments + upwind inflow | `(2^d, ng)` moment source + face-cochain inflow traces |
| `Ω⃗·n⃗ ≷ 0` | outflow / inflow face split | ORPHEUS sweep pre-resolution (per axis: `s_a = |μ_a|/Δ_a`) |
| `M` (mass) | `∫ B_i B_j` collision term | Kronecker of 1-D LD mass; diagonal weights `1, θ, θ, θ²` (2-D) |
| `G` (gradient/stiffness) | `∫ B_i Ω·∇B_j` streaming | Kronecker; carries the `s_a` per axis |
| `F_in/F_out` (surface) | upwind face integral | the `2^{d-1}`-moment face cochain (Item 3) |
| `σ_t`, `S` | total XS, total source | `total_xs`, `source` (now `(2^d, ng)`) |
| simplex-P1 `{1,x,y,z}` | WMMP-2001 tetra basis | NOT the Cartesian object — use only on a genuine simplex mesh |

**Source-moment generalization (the brief's `(1+d, ng)` question):**
the moment source is **`(2^d, ng)`**, NOT `(1+d, ng)` — ordered as
the tensor-Legendre moments `(S̄, Ŝ_x, Ŝ_y, Ŝ_xy)` in 2-D,
`(S̄, Ŝ_x, Ŝ_y, Ŝ_z, Ŝ_xy, Ŝ_xz, Ŝ_yz, Ŝ_xyz)` in 3-D. The 1-D
`(2, ng)` source is the d=1 case (`S̄, Ŝ_x`). The current 1-D code's
`_ld_source_moments` splits `(2, ng)` → (average, slope); the multi-D
generalization splits `(2^d, ng)` → the full tensor-moment set.

---

## 7. Gaps / risks for the D5b implementer

1. **The naïve `1+d` simplex basis is a CORRECTNESS TRAP on Cartesian
   cells** (Adams-2001 — fails thick diffusion limit). Build the
   `2^d` tensor-product bilinear/trilinear (UBLD) object. The 1-D
   code's own comment ("a multi-D LD is bilinear") is correct and
   under-stated — "bilinear" means the **4-moment** Q1 object, not
   a 3-moment affine object.

2. **`CellResult`/face-cochain contract widening.** Cell unknown is
   `2^d`-long; each face carries `2^{d-1}` moments. Both exceed the
   current scalar slots. This is the biggest structural item —
   architecture review BEFORE coding (same gap as 1-D, now
   dimensionally amplified). The Schur-scalar collapse that let 1-D LD
   ride the existing scalar contract does NOT cleanly survive: the
   multi-D cross moment couples the per-axis slopes, so the cell solve
   is a genuine `2^d × 2^d` dense solve, not a sequence of scalar
   Schur reductions. (You CAN still Schur-eliminate all non-average
   moments to a scalar `ψ̄` balance, but the eliminated block is
   `(2^d−1)×(2^d−1)` and couples axes — verify whether the
   affine-scannable property survives; likely it does NOT in multi-D
   because the face object is no longer scalar, defeating the
   single-upstream affine recurrence.)

3. **Assemble via Kronecker products, do NOT hand-transcribe.** The
   M/G/F matrices are Kronecker products of the verified 1-D LD
   operators. This is the single-source, elegant route; it makes the
   d=1 reduction an identity and the d=2/3 cases fall out. SymPy-verify
   against the "exact on a bilinear/trilinear flux" oracle.

4. **Slope-source threading (Increment C dependency).** The
   diffusion-limit consistency (Item 5) needs the full `2^d` source
   moments including the scattering-source slopes — the same
   global-moment-contract dependency as 1-D. The flat-source multi-D
   UBLD is O(h²) but NOT diffusion-limit-consistent until the moment
   source iterate is threaded.

5. **Positivity stays False; lumping is a separate later option**
   (Item 4). Do not conflate FLBLD with a positivity guarantee.

6. **Zotero was DOWN this session** — could not check the user's
   library or annotations for Adams-2001 / MRM-2016 / WMMP-2001 /
   BLA-1992. No user-annotation notation oracle available; the
   notation map above is from the published PDFs (MRM-2016 read
   directly). Re-check Zotero annotations when the server is back —
   if the user highlighted the UBLD Eq.(6)/(12) cluster, that is the
   canonical form to match.
