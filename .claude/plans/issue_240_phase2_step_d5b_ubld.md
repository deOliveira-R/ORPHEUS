# #240 Phase 2 Step D5b — N-dim Linear-Discontinuous on the DAG wavefront (the UBLD)

> **Durable in-repo recovery anchor** (project rule: plans live in ORPHEUS/.claude/).
> Parent: `.claude/plans/issue_240_phase2_step_d_homing.md` (§D5). Subsumes **#158 Increment D / #38**.
> Branch `feature/sn-space-angle-tier2` (off `main`@`cba6d2f`), NOT pushed/merged.
> **STATUS (2026-06-16): D5b-S2 DONE + committed (`495af60`, reviewed). Next = D5b-S3 (fold Inc C — the φ̂ spatial-moment iterate).**
> S1 (the unified per-cell UBLD primitive) is complete on both branches:
> **S1 commit chain** (branch `feature/sn-space-angle-tier2`, NOT pushed/merged):
> `7abba5d` (design plan, pre-session) → `cb84b7b` (Branch-1 SymPy algebra-of-record + oracles)
> / `9382ace` (chore) → `69b19c9` (Branch-2 production numpy primitive `orpheus/sn/spatial/_ubld.py`
> + the LD Rule-of-Three collapse: `_schur_terms`/`_kernel_terms`/`affine_scan_coefficients`
> single-source through `d1_closed_form`) / `effdfc1` (chore) → `c8489e4` (status) →
> `3567567` (hindsight `_xV` extraction) / `e080b10` (audit chore).
> **Verification:** bit-identical DD negative control (513 strict-gate pass, no golden `.npy` moved);
> LD principled ~1-ULP re-baseline (no golden moved); the d-generic dense `assemble_ubld`/
> `per_cell_solve` primitive is built + d=2-exact-on-bilinear verified (the S2 d≥2 path) but NOT yet
> wired (the `linear_discontinuous.py` `len(s_axes)!=1` raise STAYS). ERR-060 caught (dropped
> |μ_axis| inflow factor). elegance PASS + qa SUPPORTED both branches.
> **Hindsight audit (post-S1, independent elegance-enforcer): PASS — the code landed in the right
> place.** `d1_closed_form`/`D1ClosedForm` correctly live in `_ubld.py` (co-located with the dense
> reduction they're proven `==` to — NOT moved into `LinearDiscontinuous`); the dense primitive
> ahead of its S2 consumer is the justified reference-oracle exception (it caught ERR-060); the
> sympy↔numpy split is the justified algebra-of-record (do NOT single-source the 1-D factors —
> defeats structural independence; `balance.py` confirmed clean of a 3rd copy). ONE real duplicated
> expression FIXED (`3567567`: `schur_xV`/`scan_xV` shared `g·V`/`V·eff_denom` → `D1ClosedForm._xV`).
> ⚠ The ERR-060 `error_catalog.md` entry is written in the working tree but UNCOMMITTED
> (`.claude/skills/` is in the standing forbidden set — land via the instruction-architecture flow;
> the `catches("ERR-060")` marker IS committed in the test).
>
> ✅ **S2 (#60) DONE + committed (`495af60`; chore records in the follow-up `chore(claude)`).**
> The d≥2 UBLD LD kernel runs on the DAG wavefront: `cell_kernel_batch`/`residual_kernel_batch`
> fork on `len(s_axes)` — d=1 stays on `d1_closed_form`→CumprodScan (no perf regression), d≥2
> assembles `A=G+F_out+Σ_t·M` (Kronecker, `_ubld_system`) + the per-axis upwind inflow
> (`_ubld_inflow`→`assemble_inflow_axis`) and solves the batched `2^d` system (`per_cell_solve`);
> outgoing faces = the tensor-Legendre downstream-node trace (sum of `o_a∈{0,1}` blocks,
> `_ubld_outgoing_faces`). The contract is indexed by ONE trait
> `spatial_basis_per_axis` (DD/Step=1, LD=2; per-cell=`per_axis**d`, per-face=`per_axis**(d-1)`);
> the `is_multi_moment` (per_axis>1) gate widens the interior face cochain
> (`_MovingFrontier._win`, `FullFieldWavefront._octant_face_cochain`) + the moment-reducing emit
> (`_CellSolve`/`_CellResidual`) — DD/Step **byte-identical** (no length-1 axis appended). `Q_cells`
> scalar source is lifted to slot 0 (`AVERAGE_MOMENT`) inside `_ubld_system` (no `_SolveOperands.Q`
> carrier change this step). Single sources minted in `_ubld.py`: `AVERAGE_MOMENT` (slot-0 layout)
> + `face_moment_tail` (the two storage policies). The D5b.0 pin INVERTED
> (`test_cell_kernel_batch_admits_multi_d`).
> **Gates green:** DD bit-id strict (sweep/core+solve+cartesian_2d) **562 / 2 skip / 4 xfail** no
> golden moved; D5b.1 round-trip (full `2^d`, non-flat inflow, het); D5b.4 kernel matvec twin (faces
> both directions); D5b.3 two-paths FFW≡MFW (NON-SQUARE, forced reps); D5b.5 DD≠LD routing-flip; the
> carry-forward d=2 numpy↔symbolic A==A CELL pin + `test_d2_exact_on_bilinear` (the genuine ERR-060
> inflow catcher); 2-D LD MMS O(h²) smoke (vacuum, `@l1` NO `verifies` — Mode-7 honest); LD
> spatial+primitive+symbolic 40; Sphinx clean (matrix regen). **Reviews:** elegance PASS-WITH-NITS
> (orphaned `is_affine_scannable` docstring; `AVERAGE_MOMENT`/`face_moment_tail` single-source; gate
> invariant comment) ALL FIXED; qa SUPPORTED-WITH-CONCERNS (dropped a misattributed
> `catches("ERR-060")` off the cell-A==A pin; `_ubld_inflow` `|None` seed) ALL FIXED.
> ⚠ The method-implementer DIED on an API-overload AFTER completing the impl+tests+crosswalk+docs
> (NOT mid-logic); the main agent triaged the tree (imports+gates green), ran the reviews, fixed all
> nits, and committed. The crosswalk `.claude/plans/issue_240_d5b_s2_crosswalk.md` is the design
> record (moment ordering `[bar,ŷ,x̂,x̂y]`, the Q lift, the domain-edge minimal touch).
>
> ⏭ **S3 (#61) ENTRY POINT — fold Increment C (the φ̂ spatial-moment iterate; closes the diffusion
> limit for d=1 #37 AND d≥2).** PROACTIVE `test-architect` dispatch FIRST (cross-subsystem carve —
> the existing memo covers D5b.0–.5 but NOT the spatial-moment-iterate FP-invariance; write the
> crosswalk, L17). Sites: `_CellSolve.cell` accumulate `φ̂` (today reduces on the average moment
> only — `psi_avg[..., AVERAGE_MOMENT]`); the between-sweep flux field gains a
> `(n_spatial_moments, ng, *spatial)` slot; `_SolveOperands.Q` grows the trailing `2^d` moment axis
> (the genuine MMS/scattering moment-source carrier S2 deferred — S2 lifts a SCALAR source inside the
> kernel); `ScatteringOperator.apply` lifts over the spatial-moment axis (`Σ_s⊗I`); close the d≥2
> Krylov 2-D LD matvec raise (`_CellResidual.cell` — it RAISES today, the loud-fail S2 left). DD/Step
> untouched (n=1 → `Σ_s⊗I_1`). **Gates:** the multi-D `ld-thick-diffusive` tripwire flips xfail→PASS
> + the 1-D #37 thick-diffusive flips xfail→PASS; **vv Mode-9 FP-invariance** (the moment-source
> splitting MUST NOT move the converged FP — verify on an ANISOTROPIC/heterogeneous config + a
> DIAGONAL cubature, NOT the isotropic-reflective box). **S4 (#62)** then mints the strengthened
> Mode-7 stress-ansatz MMS + the non-vanishing-boundary `BoundaryFlux` moment trace + the
> `@verifies("ld-cartesian-2d")` flux-shape claim (S2 left the smoke unverified-by-design).
> **S2 carry-forwards STILL OPEN for D6:** replace production line-number citations in the
> `ld_ubld.py` oracle comments with symbol citations; mint/link the UBLD orphan equation family
> (`ld-ubld-n-spatial-moments`, `ld-ubld-divv-scale-free-kernel`, + the S1 set, + `ld-cartesian-2d`).
> The design pass (plan mode) resolved the two open architecture questions WITH THE USER and
> EXPANDED the scope — see "⭐⭐ THE TWO DECISIONS" below. The approved implementation plan is
> `.claude/plans/mellow-swinging-breeze.md` (the S1–S5 sub-step sequence + the unified
> `n_spatial_moments` lattice). Tasks: #58 (umbrella, in_progress) → #59 (S1, in_progress) →
> #60 (S2) → #61 (S3, the Inc C fold) → #62 (S4, MMS) → #55 (D6). D5b now subsumes BOTH
> #37 (Inc C, FOLDED IN) and #38 (Inc D, = this).
> D5-0 + D5a DONE (the routing substrate is honest). The D5 design pass (test-architect spec +
> literature formulation + cross-domain-attacker frame check) is DONE — the three memos below
> are the authoritative inputs.

## ⭐⭐ THE TWO DECISIONS (design pass 2026-06-16 — do not lose across compaction)

**Decision 1 — FOLD Increment C (#37) into D5b.** D5b ships the d-generic *global
moment-source iterate* (the scattering-slope feedback `φ̂`), closing the diffusion limit for
BOTH d=1 (#37) AND d ≥ 2 in ONE arc. (User chose "Fold Inc C into D5b" over "defer".)

**Decision 2 — UNIFY the per-cell algebra WITHOUT a 1-D perf regression** (user: "why can't we
have unification while keeping performance?"). Single-source the UBLD per-cell algebra as ONE
d-generic **primitive** (Kronecker-assembled M/G/F + moment elimination); *derive* the d=1
`affine_scan_coefficients` from it (the closed-form `2×2` Schur = the affine scan recurrence) so
d=1 keeps riding **CumprodScan** (cumprod-vectorized; zero 1-D regression; L16-safe); d ≥ 2
rides the wavefront via a **batched `2^d×2^d` solve** (`np.linalg.solve` over the anti-diagonal
cell stack — vectorized across cells, NOT a per-cell Python loop). Do NOT retire the scan path —
re-derive it from the primitive.

**The unifying architecture — ONE `n_spatial_moments` reduction lattice.** A class-level scheme
trait `n_spatial_moments: ClassVar[int]` (DD/Step = 1, LD-d = `2^d`) indexes the whole contract:
cell unknown (`ψ̄` → `2^d`), face payload (scalar → trailing `2^{d-1}` moment axis), iterate flux
(`φ̄` → `2^d` spatial moments), scattering (`Σ_s` → `Σ_s ⊗ I_{spatial-moment}`), source build
(`(ng,)` → `(2^d, ng)`). DD/Step at `n=1` are the trivial reduction — the carve MUST keep them
**bit-identical** by construction (the backward-compat invariant gating every sub-step).

**Two structural findings that de-risked the carve (the two explorer passes):**
- `WavefrontFlux`/`InteriorFaceSpace` were RETIRED at S6.4(f). The interior face cochain is now
  raw per-axis ndarray tuples, and the DAG walk (`SweepDependencyGraph.walk_full`/`walk_windowed`)
  is ALREADY moment-agnostic and ALREADY exercises the per-axis multi-face collection (DD 2-D
  drives it). ⇒ NO new face type — append a trailing moment axis to the per-axis ndarrays; the
  walk plumbing carries it untouched; the KERNEL owns the `2^d×2^d` solve.
- Today ONLY `φ̄` exists anywhere (iterate, `ScatteringOperator`, production cell kernel are all
  single-spatial-moment; the LD `(2,ng)` two-moment path is reachable only from a UNIT TEST). ⇒
  Inc C *introduces* a spatial-slope contract; it does not extend one that already flows. Minimal
  sites: scalar-flux reconstruction (`_CellSolve.cell`, accumulate `φ̂`); the iterate flux field
  (new `(n_spatial_moments, ng, *spatial)` slot); `ScatteringOperator.apply` (lift over the
  spatial-moment axis); `Q_cells` (grow a slope-source row).

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

## DESIGN QUESTIONS — RESOLVED (design pass 2026-06-16)

1. **The face-cochain type — DISSOLVED by the code's current state.** `WavefrontFlux`/
   `InteriorFaceSpace` are RETIRED; the face cochain is raw per-axis ndarrays and the walk is
   already moment-agnostic. RESOLUTION: append a trailing `2^{d-1}` moment axis to the per-axis
   face payload; the kernel owns the `2^d` solve; the 1-D scalar path is the d=1 (trailing-dim-1)
   reduction. NO new type. (See "⭐⭐ THE TWO DECISIONS" + the `n_spatial_moments` lattice above.)
2. **D5b vs Increment C — USER CHOSE TO FOLD INC C IN.** D5b ships the d-generic moment-source
   iterate (`φ̂`); the multi-D thick-diffusive tripwire PASSES within D5b, and the same contract
   closes 1-D #37. Sequenced as S2 (flat scattering, external/MMS moment source threaded) → S3
   (the scattering-slope feedback `Σ_s·φ̂`). So #37 does NOT block S1/S2; it IS S3.
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
