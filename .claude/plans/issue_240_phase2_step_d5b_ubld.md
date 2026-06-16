# #240 Phase 2 Step D5b — N-dim Linear-Discontinuous on the DAG wavefront (the UBLD)

> **Durable in-repo recovery anchor** (project rule: plans live in ORPHEUS/.claude/).
> Parent: `.claude/plans/issue_240_phase2_step_d_homing.md` (§D5). Subsumes **#158 Increment D / #38**.
> Branch `feature/sn-space-angle-tier2` (off `main`@`cba6d2f`), NOT pushed/merged.
> **STATUS (2026-06-16): NOT STARTED. D5-0 + D5a DONE (the routing substrate is honest).**
> D5b is the THIRD and genuinely-new-math thread of the D5 campaign. It needs **its OWN
> design pass** (the cell/face contract widening is an architecture decision — Cardinal Rule 2)
> BEFORE implementation. The D5 design pass (test-architect spec + literature formulation +
> cross-domain-attacker frame check) is DONE — the three memos below are the authoritative inputs.

## ⭐⭐⭐ THE LOAD-BEARING FINDING (do not lose this across compaction)

The plan's original framing — "LD's multi-D closure is bilinear = average + **one slope per axis**
(`1+d` moments)" — is **WRONG for Cartesian cells**. The literature pass (Adams 2001, NSE 137;
Maginot-Ragusa-Morel 2016 = **MRM-2016**, NSE 185; Börgers-Larsen-Adams 1992, JCP 98) established a
**design fork**:

* **Simplex-P1 LD** (`{1, x, y}`, `1+d` moments) — the literal "one slope per axis." On a
  **quadrilateral/Cartesian cell it FAILS the thick-diffusion limit** (Adams 2001). It is correct
  ONLY on unstructured triangle/tetra meshes (WMMP-2001).
* **Bilinear/trilinear DG-P1 = UBLD** (`{1, x, y, xy}` in 2-D, `{1,x,y,z,xy,xz,yz,xyz}` in 3-D —
  the **tensor product of 1-D linear Lagrange functions**, `2^d` moments) — **PRESERVES the
  thick-diffusion limit on quadrilaterals**. The `xy` **cross moment** is exactly what simplex-P1
  lacks and is the reason bilinear survives.

**⇒ D5b builds the tensor-product UBLD (`2^d` moments), NOT the simplex `1+d`.** Building the naive
`1+d` object would ship a SILENT physics bug (a scheme that fails the diffusion limit but passes
every non-diffusion test). The `ld-thick-diffusive` tripwire is the catcher: it flips xfail→PASS on
UBLD, stays FAIL on simplex (correct physics — document the basis-dependence).

**Per (group, ordinate): 4 moments in 2-D, 8 in 3-D.** Moment basis = tensor-product Legendre
`{ψ̄, ψ̂_x, ψ̂_y, ψ̂_xy}` (2-D). The current 1-D LD `{ψ̄, ψ̂}` (θ=1/3 Legendre) is the d=1 reduction.

## Authoritative inputs (READ FIRST — these three memos are the contract)

1. **Formulation** — `.claude/agent-memory/literature-researcher/multi_d_ld_closure.md` (the per-cell
   UBLD system, the Kronecker-assembly recommendation, the face-contract widening, the lumping +
   diffusion-limit verdicts, the citations). **THE implementer's contract for the math.**
2. **Verification** — `.claude/agent-memory/test-architect/d5_nd_polymorphism_verification.md` §"D5b"
   (the gates D5b.0–D5b.5) + §3 (the MMS design) + §1 (the coverage matrix) + §4 (the sequence).
3. **Frame check** — `.claude/agent-memory/cross-domain-attacker/d5_trait_and_mms_frames.md` Frame 2
   (the MMS activation analysis + **the slope-driver strengthening** — see §MMS below).

## The math — the per-cell UBLD system (MRM-2016 §2, Eqs. 1–12)

Galerkin weak form of `Ω·∇ψ + Σ_t ψ = S` with upwind flux on the INCOMING faces, per cell K:

```
(12)   A_UBLD · ψ⃗ = R⃗            A = G + F_out + Σ_t·M       (2^d × 2^d, non-symmetric, DENSE)
```
- `M_ij = ∫_K B_i B_j` — mass (collision), tensor-product Legendre → diagonal weights `θ^{#active axes}`
  (`1` for ψ̄, `θ` for ψ̂_x/ψ̂_y, `θ²` for ψ̂_xy).
- `G_ij = ∫_K B_i (Ω·∇B_j)` — gradient/streaming (couples all `2^d` moments).
- `F` — surface: OUTFLOW (`Ω·n>0`, implicit, into `A`) + INFLOW (`Ω·n<0`, **upwind**: upstream
  neighbour's outflow trace → into `R`).
- `R⃗ = M·S⃗_moments + F_in·ψ_in_traces` (the source moments + the upwinded inflow-face traces).

Moment-balance reading (ORPHEUS notation, `s_a=|μ_a|/Δ_a` raw `g`):
- **0th (ψ̄)** = cell balance — the multi-axis generalization of the 1-D `|μ|(ψ_out−ψ_in)+Σ_t h ψ̄=Q̄h`.
- **1st (ψ̂_x, ψ̂_y)** = per-axis slope rows — each the 1-D slope row ON THAT AXIS,
  **PLUS coupling to ψ̂_xy through the transverse face integral** (the multi-D content).
- **cross (ψ̂_xy)** = the bilinear row — NO 1-D analog; makes the system 4×4 not 3×3; the
  diffusion-limit-load-bearing term.

## ⭐ The architecture — the elegant build (Kronecker assembly) + the contract widening

**The build (single-source, Pattern 5 "build the primitive not the product"):** assemble `M, G, F`
as **Kronecker products of the verified 1-D LD operators** (which ORPHEUS already has). The 2-D mass
matrix is `M_1D ⊗ M_1D`, etc. This (a) single-sources the 1-D math, (b) makes the **d=1 reduction a
Kronecker-with-empty-axes identity** (the existing `linear_discontinuous.py` IS the d=1 case),
(c) the diffusion-limit-load-bearing `xy` coupling falls out of the linear algebra rather than
hand-transcription. **Do NOT hand-transcribe the 4×4 / 8×8 entries.** SymPy-verify the assembly
against (i) the 1-D per-axis reduction and (ii) the "exact on a bilinear flux `ψ=a+bx+cy+dxy`" oracle
(the multi-D analog of the 1-D "exact on linear-in-x").

**⚠ THE CONTRACT WIDENING (the architecture-review item — Cardinal Rule 2):**
- The per-cell unknown is a **`2^d`-vector**, not a scalar `cell_average_flux`. `CellResult` carries
  `cell_average_flux` + ONE `outgoing_spatial_flux` today — no slot for the `2^d−1` slope/cross moments.
- Each downstream face is a **`2^{d-1}`-moment object** (average + transverse slope), NOT a scalar.
  In 2-D the downstream-x face carries a linear-in-y profile
  `ψ_out^x(y) = (ψ̄+ψ̂_x) + (ψ̂_y+ψ̂_xy)·P₁(y)`. The interface-flux cochain (`WavefrontFlux` /
  the moving-frontier face storage) must carry per face a `2^{d-1}`-moment object.
- This is a genuine widening over the DD face (scalar) AND the current 1-D LD face (scalar, because
  the 1-D downstream face is 0-dimensional). **The face-cochain type is the load-bearing extension —
  it is the main design-pass decision.** The 1-D Schur-scalar collapse likely does NOT survive to
  multi-D.

**Where it plugs in (LD rides the wavefront, NOT the scan-march):**
- Close the `linear_discontinuous.py:430-436` raise (`len(s_axes)!=1`) by supplying the multi-axis
  bilinear kernel in `cell_kernel_batch`/`residual_kernel_batch`.
- LD's 2-D path is `FullFieldWavefront`/`MovingFrontierWindow` — ALREADY pure scheme-delegators
  (`sweep_graph.py:849/890` → `scheme.cell_kernel_batch`/`residual_kernel_batch`). So closing the
  kernel raise + widening the face cochain makes LD run N-D on the wavefront with (ideally) no
  loss-rep change.
- LD stays the **structural EXCL on the scan-march** (`transverse_coupling_is_facewise=False`,
  landed in D5-0): the `affine_scan_coefficients` raise on non-neutral curvature is NOT inverted —
  only the DAG-kernel `len(s_axes)!=1` raise inverts.

**Lumping:** ship **UBLD (unlumped)** first (most accurate; the object the diffusion-limit proofs are
stated for; clean linear `Aψ=R` solve). FLBLD/SCB (Wareing 2001 / Adams) = a later robustness option
for void/grazing (never default — large numerical diffusion). **`is_positivity_preserving=False`
STAYS** — NO linear bilinear DFE (UBLD or FLBLD) is strictly positive; strict positivity needs a
NON-linear scheme (defer, like the 1-D positivity fix-up).

## Current code surface (file:line — verify at pickup; D-series moved things)

- `orpheus/sn/spatial/linear_discontinuous.py`: `_kernel_terms` (`:411`, **the raise `:430-436`**),
  `cell_kernel_batch` (`:450`), `residual_kernel_batch` (`:475`), `_ld_source_moments` (`:200`,
  requires `(2,ng)` → generalize to `(2^d, ng)`), `_LDCellTerms` (`:223`), `_schur_terms` (`:295`),
  `theta=1/3` (`:286`). `is_affine_scannable=True` (`:278`), `transverse_coupling_is_facewise`
  inherits `False` (the D5-0 EXCL).
- `orpheus/sn/spatial/scheme.py`: the `cell_kernel_batch`/`residual_kernel_batch` contract; `CellResult`;
  the reconstruction staticmethods (`source_emission`/`cell_average`/`outgoing_face_from_average`).
- `orpheus/sn/sweep_graph.py:849/890`: the scheme-delegator call sites (the walk passes the cell to
  `scheme.cell_kernel_batch`/`residual_kernel_batch`); the `WavefrontFlux` / `_MovingFrontier` face
  storage = the contract-widening target.
- `orpheus/sn/spatial/diamond.py:293-401`: DD's d-generic kernel (the explicit `zip(s_axes,psi_in)`
  left-fold) — the STRUCTURAL shape model, though LD's is a coupled `2^d` Schur solve, not a left-fold.
- The MMS home: `derivations/continuous/mms/sn.py` (`SN2DCartesianMMSCase` `:642` is the mirror;
  `_default_hetero_2d_xs_functions` `:950`); `tests/sn/verification/mms/test_mms_ld_slab.py` (1-D, the
  pattern); NEW `tests/sn/verification/mms/test_mms_ld_2d.py`.

## Verification — the D5b gates (test-architect spec §2-D5b; new capability, MMS-gated, NOT bit-id)

- **D5b.0 — INVERT the raise pin.** `test_linear_discontinuous.py::TestLDKernel::
  test_cell_kernel_batch_rejects_multi_d` (`:361-367`) asserts the `d=1` raise → REPLACE with
  `test_cell_kernel_batch_admits_multi_d` (the `(1+2)×…`/`2^d` system runs, returns
  `(psi_avg, (out_x, out_y))`, `len(faces)==2`). A false gate if left.
- **D5b.1 — round-trip (foundation).** `residual_kernel_batch` at the `psi_bar` `cell_kernel_batch`
  solves for vanishes. NON-flat per-axis `psi_in` (the bilinear slope per axis ACTIVE); feed back the
  FULL solved `2^d` state (partial feed passes spuriously). `n_groups∈{1,2}` het. `atol≈1e-12` (the
  `2^d` solve is a few division ULP deeper than DD's `1e-13` — document the bound). Catches #4
  slope-elimination Schur index drift, #2 x-slope↔y-slope swap.
- **D5b.2 — multi-D LD MMS O(h²) (L1, the structurally-independent ref).** NEW
  `test_mms_ld_2d.py`. Driver `solve_sn_fixed_source(..., Mesh2D(...), scheme=LinearDiscontinuous())`
  (`scheme=` threads end-to-end, `solver.py:1972`→`_as_sn_mesh`); the 2-D LD mesh routes to
  `FullFieldWavefront`. Ladder `nx=ny∈{20,40,80}` (or non-square) → `orders[-1]>1.95 ∧ all>1.85` +
  the VALUE band vs `phi_exact` (rate≠correctness). `@l1 @slow @verifies("ld-cartesian-2d",
  "transport-cartesian-2d")` — **`ld-cartesian-2d` is a NEW label D6 mints.** The MMS = §MMS below.
- **D5b.3 — two-paths FFW≡MFW (foundation, the headline invariant).** Same LD via two DAG schedules
  (both pure delegators → same `cell_kernel_batch`; agreement proves the walk storage doesn't perturb
  LD's cell math). Mode-9 config (the §MMS stress ansatz: het, 2G asym, μ-non-trivial, NON-SQUARE
  `nx≠ny`). `rtol≈1e-9` (SAFETY×conv_tol) + a single-sweep `assert_array_almost_equal_nulp`. ⚠ VERIFY
  each leg ran its rep (`default_for` is the expected wavefront; force MFW for the oracle leg) — a
  two-paths gate where both legs secretly ran the same rep is a silent false green.
- **D5b.4 — Krylov≡SI matvec twin (foundation).** LD forward matvec (`residual_kernel_batch`) ≡ LD SI
  sweep (`cell_kernel_batch`) on the 2-D wavefront. ⚠ **HAZARD:** the matvec reconstructs the outgoing
  face — the multi-D LD must supply a WORKING `outgoing_face_from_average` in BOTH directions per axis.
  **A sweep without a verified matvec is half a carve (L14 leg-3 standoff).** `rtol=1e-9,atol=1e-11`.
- **D5b.5 — DD≠LD routing-flip (foundation, the cross-thread guard).** Drive `solve_sn_fixed_source`
  twice (same mesh/materials/source), `scheme=DiamondDifference()` vs `LinearDiscontinuous()`; assert
  the converged fluxes DIFFER at coarse mesh (`not np.allclose(phi_dd, phi_ld, rtol=1e-3)`). **This is
  the catcher that LD is GENUINELY computing LD, not DD** — if the D5-0 misroute regressed, this sees
  DD≡LD and FAILS. (No `verifies` — a "these are different schemes" contract.)
- **Coverage-matrix delta (definition-of-done):** LD `FullFieldWavefront d=2` + `MovingFrontierWindow
  d=2` flip `NOT-YET → T`. **d=3 LD stays `NOT-YET`** (out of D5b scope — flag/track, #227-adjacent).
- The `ld-thick-diffusive` tripwire (1-D today, `test_mms_ld_slab.py:227`): a multi-D analog flips
  xfail→PASS on UBLD (stays FAIL on simplex — the basis-choice catcher).

## The multi-D LD MMS design (vv Mode-7 override MANDATORY — the hardest gate to get right)

**FORBID the existing isotropic 2-D MMS** (`SN2DCartesianMMSCase`, `phi=sin·sin`): it is angularly
flat → NULLS the per-axis slope rows → tests NOTHING about the bilinear coupling (the exact thing D5b
introduces). AI has no SymPy derivation cost → pick the stress ansatz.

**Trial (per ordinate n, group g):** `ψ_{n,g}(x,y) = [ A_g + μ_{x,n}·B_g + μ_{y,n}·C_g ] / W`
with mixed-scale spatial drivers. The test-architect base ansatz:
```
A_g = a0 + a1·sin(πx/Lx)sin(πy/Ly) + a2·cos(2πx/Lx)cos(3πy/Ly)
B_g = b0 + b1·sin(2πx/Lx)sin(πy/Ly)        # x-slope driver
C_g = c0 + c1·sin(πx/Lx)sin(2πy/Ly)        # y-slope driver
```
**⭐ cross-domain-attacker STRENGTHENING (Frame 2 — MANDATORY, the base ansatz is INSUFFICIENT):**
`B` and `C` as written are **x↔y reflections of each other** → a same-sign slope-row transcription
bug (the LIKELY bug, shared code path) **cancels in the measured flux → false green.** Fix: add
x↔y-**broken** cross-harmonics to the SLOPE drivers `B` and `C` THEMSELVES (not only `A`), with `B`
and `C` **NOT reflection-related**, e.g. `B += b2·cos(πx)cos(2πy)`, `C += c2·sin(3πx)cos(πy)`. KEEP
`a0>0` (boundary non-vanishing — the BC closure must satisfy non-trivial interior), 2G asymmetric
downscatter, NON-SQUARE domain. (A first test: weak vs strengthened ansatz against a deliberately
same-sign-flipped LD slope row — weak still converges O(h²) = false green; strengthened breaks it.)

**Activated/nulled (declare in the test, Mode 7):** ACTIVATES streaming(both axes) + collision +
both per-axis slopes + the cross moment + BC + group coupling; NULLS nothing.

**Derivation (Branch 1, L11):** a NEW `SN2DCartesianLDStressMMSCase` in `derivations/continuous/mms/
sn.py` — substitute into `μ_x ∂_x ψ + μ_y ∂_y ψ + Σ_t ψ = (1/W)(Σ_s^T φ + Q^ext)`, solve
SYMBOLICALLY (SymPy) for `Q^ext_n`. `φ_g = ∫ψ dμ = A_g` (the μ_x·B+μ_y·C terms integrate to 0 over a
symmetric quadrature); the streaming derivative carries the full μ-weighted ansatz. Pin the symbolic
source with a `@foundation` derive-test BEFORE consuming it. Structurally independent of the LD kernel.

**⚠ The Q̂≠0 slope-SOURCE posture depends on whether D5b threads the MOMENT SOURCE** (see open question
below): if D5b is flat-source (Q̂=0, like 1-D Increment B), the MMS verifies the slope-UNKNOWN sign
only (document "slope-SOURCE deferred to the moment-source increment"); if D5b threads the moment
source, the MMS MUST supply Q̂≠0 per axis (the strengthened B/C make this natural — the MMS is
slope-source-READY). **The method-implementer declares the posture; the MMS Q̂ follows it.**

## Implementation sequence (test-architect §4, adapted — D5-0/D5a already DONE)

`[Phase 0 lit-formulation ✅ DONE] ∥ [Phase 1 cross-domain MMS frame ✅ DONE]` → D5-0 ✅ → D5a ✅ →
**D5b (here)** → D6. The two proactive dispatches that FEED D5b are already done (the memos above).
At D5b pickup:

0. **DESIGN PASS FIRST (the contract widening).** This is the architecture decision the user must
   shape: the `CellResult`/face-cochain (`WavefrontFlux`) type extension (scalar → `2^d` cell /
   `2^{d-1}` face moments); whether the 1-D scalar path unifies with or stays parallel to the multi-D
   moment path; where the `2^d×2^d` per-cell solve lives. Likely a `cross-domain-attacker` or
   `Plan`-agent pass + user sign-off. ENTER PLAN MODE.
1. **The Kronecker M/G/F assembly + the per-cell `2^d` solve** (method-implementer, against the
   literature contract; SymPy-verify the d=1 reduction + the bilinear-exactness oracle).
2. **Widen the face cochain + wire `cell_kernel_batch`/`residual_kernel_batch` for `len(s_axes)>1`;
   close the raise; INVERT the D5b.0 pin.**
3. **The MMS module** (SymPy `SN2DCartesianLDStressMMSCase`, strengthened ansatz) + the D5b.1–D5b.5
   gates. elegance-enforcer + qa per sub-step.
4. **(deferred)** FLBLD/SCB robustness option (a follow-up issue); d=3 LD (a follow-up issue).
5. **D6** (archivist): mint `ld-cartesian-2d` + link the orphan labels (`loss-rep-facewise-separable`,
   `loss-rep-scanmarch-solve-affine`, `loss-rep-scanmarch-apply-residual`, `ld-cartesian-1d`,
   `ld-slab`); expand the LD theory `.. todo:` stub (`discrete_ordinates.rst:1498`) into the rich
   UBLD derivation + the Step/DD/LD↔advection narrative.

## OPEN DESIGN QUESTIONS (resolve in the design pass — these need the user)

1. **The face-cochain type (THE decision):** how do `CellResult` + `WavefrontFlux` carry the `2^d`
   cell / `2^{d-1}` face moments? A typed moment-vector? Does the 1-D scalar path reduce as the d=1
   case (one type, Kronecker-empty-axes) or stay a parallel optimization? (The Kronecker build wants
   ONE type that reduces; the perf path may want the 1-D scalar kept — the "construct general, select
   narrow" principle of `sn_sweep_strategy.md`.)
2. **D5b vs 1-D Increment C (#37, task #37, PENDING):** Inc C threads the global moment-source
   (`φ̂` in the iterate). The UBLD `2^d` solve is inherently a full-moment solve (it solves for ψ̂_xy
   etc.), so its SOURCE wants moment components (`S̄, Ŝ_x, Ŝ_y, Ŝ_xy`). **Does D5b REQUIRE Inc C first
   (the moment-source iterate), or can D5b ship flat-source (Q̂=0 scattering, the moments solved as
   unknowns) and defer the slope-SOURCE to Inc C?** This determines the MMS Q̂ posture AND the
   sequencing (does #37 block #38/D5b?). Likely: D5b can ship the UBLD CELL SOLVE flat-source (the
   diffusion limit is about the CLOSURE, gated by `ld-thick-diffusive`), with the moment-source
   scattering as a paired/follow-on increment. **Confirm with the user / a design pass.**
3. **UBLD-only vs UBLD+FLBLD:** ship unlumped UBLD first (decided); FLBLD/SCB = a tracked follow-up.
4. **d=2 only vs d=2+d=3:** D5b targets d=2 (the matrix is 4×4, the MMS is 2-D). d=3 (8×8, trilinear)
   is deferred (no 3-D quadrature production path; #227-adjacent). The Kronecker assembly is d-generic
   by construction, so d=3 is "widen the assembly + a 3-D MMS" later.

## Citations (algebra-of-record)
- **MRM-2016** Maginot, Ragusa & Morel (2016), NSE 185(1):17-42, `10.13182/nse16-38` — the UBLD/FLBLD
  derivation, Eqs. (1)-(12), basis (8a-d). OSTI preprint `osti.gov/pages/servlets/purl/1343843`.
- **Adams-2001** Adams (2001), NSE 137(3):298-333, `10.13182/nse00-41` — the thick-diffusion-limit
  verdict (simplex FAILS / bilinear PASSES on quadrilaterals). THE `ld-thick-diffusive` reference.
- **BLA-1992** Börgers, Larsen & Adams (1992), JCP 98(2):285-300, `10.1016/0021-9991(92)90143-m` —
  the original 2-D LD asymptotic diffusion-limit analysis.
- **WMMP-2001** Wareing, McGhee, Morel & Pautz (2001), NSE 138(3):256-268, `10.13182/nse138-256` —
  the 3-D LLD/FLBLD on tetra (the simplex object; the lumping recipe).
- Local (d=1 reduction oracles): Larsen-Morel 1989 (JCP 83), Larsen-Morel-Miller 1987 (JCP 69) in
  `scratch/literature/`.

## Recovery
Read this file + the three memos (top of "Authoritative inputs") + the Step-D homing plan §D5. The
D5-0/D5a commit chain: `4b465b7`(D5-0 fix)→`f7e0324`/`9ddd296`(D5-0 chore/docs)→`66dbd9a`(D5a)→
`377ccb3`/`6197351`(D5a docs/chore). Related: #38 (= this), #37 (Inc C, the moment-source question),
#242 (the deferred DD curvilinear-diagonal merge), #227 (d≥3 kernels). EXIT after D5b+D6: the
Spatial×Angular tensor-product campaign (the NEXT campaign — extract the curvilinear angular
redistribution into a distinct AngularScheme) + ff-merge `feature/sn-space-angle-tier2`→`main`.
