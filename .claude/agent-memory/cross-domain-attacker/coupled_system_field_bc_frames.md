---
name: coupled-system-field-bc-frames
description: VERDICT on the refined "augmented SN = coupled 2×2 [[A_AA,A_AB],[A_BA,A_BB]], system=field+BC" proposal (ψ½ ray as System B, transport as System A). The 2×2 is NOT a new object — it is the existing 3-block composite biproduct re-associated (Mat₂∘Mat₂≅Mat₄); no general CoupledOperator; the three cited instances are three coupling KINDS (linear-triangular / R⊣P-Galerkin / nonlinear) so defer-until-2-KINDS is UNMET.
metadata:
  type: project
---

# Coupled 2×2 augmented-SN — "system = field + BC" frame verdict

Branch `refactor/sn-walk-unification` @ `4081c0d` (dirty, C2 landed); branch-verified
reads (L-005), active Nexus graph IS the branch. Supersedes the naive `FullField ⊕
RayOp` block-DIAGONAL direct sum (correctly rejected earlier). This is the SHARPENED
proposal: a coupled 2×2 block MATRIX with off-diagonals. Second memo on this composite
after [[bordered-vs-biproduct-facefield-frames]] (which ruled biproduct-NOT-bordered).

## The proposal (user, 2026-07-07)
- System A (transport): field `BulkField`, BC = spatial trace; `A_AA = L+C−S−B`.
- System B (ray/ψ½): field codim-1 `FaceField` (ψ½ radial profile), BC = r=R corner
  (the seed); `A_BB = RayOp` (banded radial DD march, Hébert 3.434–3.435).
- Coupled `[[A_AA,A_AB],[A_BA,A_BB]]`: A_AB = ψ½ seeds bulk M-M recurrence; A_BA = bulk
  moments → ray source (fold), "currently WELDED inline into S/F, a super-smell."
- Claims: (1) the seed IS a BC (r=R inflow = Dirichlet datum of the ray ODE); (2) C2
  stands (ray FIELD = codim-1 FaceField); refinement = ray is a SYSTEM (FaceField+BC),
  a general `CoupledOperator` sits above.

## THE CRUX — the 2×2 is NOT a new object (Mat₂∘Mat₂ ≅ Mat₄)
The biproduct `⊕` is coherently associative, so the coupled 2×2-of-subsystems is the
EXISTING 3-block composite biproduct `V_bulk⊕V_trace⊕V_sd` (full_field_space.py) RE-
PARTITIONED (group V_A=bulk⊕trace, V_B=sd). The ψ½ block promoted from "third ⊕
summand" to "second subsystem" is a FREE RE-ASSOCIATION. The user CONVERGED onto the
biproduct frame prior memos ruled native — did NOT diverge to a new object. The "different
object than block-DIAGONAL ⊕" is true only vs the STRAWMAN (the uncoupled direct sum);
it is the SAME as the biproduct-with-off-diagonals (A_bs seed, A_ss=−B) already established.
The G-adjoint composes block-wise FOR FREE: G = diag(G_A, G_B) block-diag per subsystem
(G_sd=V_cell post-ERR-067), so `A† = [[A_AA†, G_A⁻¹A_BAᵀG_B],[G_B⁻¹A_ABᵀG_A, A_BB†]]` is
`_AdjointOperator`'s A_bs↔A_sbᵀ law read at the 2-subsystem granularity — NOTHING new.
DO NOT mint a `CoupledOperator` type: it is a VIEW (redundant) or a twin (Smell #16). A
coupled system is an `OperatorSum` of block-placed leaves over the block index.
Discriminator: exhibit a LINEAR coupled system expressible 2×2-of-subsystems but NOT flat
3-block — impossible; every candidate is flat-re-expressible (⇒ view) or nonlinear (⇒ not
a LinearOperator ⇒ refutes multiphysics-as-instance).

## Q1 system=field+BC — clean, recursive, native name = biproduct-carrier + boundary-row
Native object = `(V_bulk⊕V_∂, A)` with the BC = A's boundary block-ROW (A_sb,A_ss). NOT
a comma/slice category (adds no morphism the biproduct+row lacks — UNEXPLORED, L-001).
RECURSES, and System B recurses FURTHER than the proposal: the ray is a 1-D TWO-POINT BVP
— boundary = r=R Dirichlet (inward leg, = global outer BC | μ=−1, RIGOROUS) + r=0 pole-
continuity ψ⁺(0)=ψ⁻(0) (outward leg, an internal Ω→−Ω REFLECTION, psi_half_angle_seed.py:126).
"BC = the r=R corner" captures ONLY the inward Dirichlet; A_BB = [inward march]→[r=0
reflection]→[outward march], the r=0 map an internal within-B off-diagonal. Sharpened
seed-as-BC test: drop the r=0 continuity (seed outward from 0) ⇒ sphere cold-residual REDs.

## Q3+Q5 defer-vs-build — NOT justified; three instances are three KINDS
| instance | off-diagonals | metric | solve | linear |
|---|---|---|---|---|
| ψ½ (now) | seed/fold (metric-adjoint) | block-diag PSD | forward-subst (triangular) | yes |
| DSA (#2) | R/P (Galerkin R⊣P) | block-diag PSD | two-way iterative (Richardson/Krylov) | yes |
| multiphysics | ∂Σ/∂T | — | Newton/Picard fixed point | NO |
No two share a KIND ⇒ defer-until-≥2-KINDS UNMET. Minimal correct object (build now):
(1) widen `BlockRole{BULK,BOUNDARY,FULL}` (operator.py:204, the 2-block enum FREEZE) to
the mesh-enumerated block SET — `_join_block_roles` is ALREADY set-union; = the FaceField
design §5.3 block-list, now with the Mat(𝒞)-N-way justification; (2) NAME A_BA — reify
[S/F moment source]∘[pole fold] as a block leaf (`fold_moments_to_starting_direction` R14
factored, the COMPOSITION welded — Smell #16 shape-4); (3) triangularity = SN-walk
REFINEMENT mixin (#284 forward-subst), never a base property.
OVER/UNDER-reach: PSD block-diag metric — OK ψ½/DSA, EXCLUDES RQI=KKT/saddle (zero corner,
INDEFINITE, eigenvalue.py:545 — stays OUTSIDE, `SaddlePointOperator` defer, trigger
mixed-form #294). Triangularity — ψ½-only. Linearity — UNDER-reaches multiphysics (nonlinear,
NOT a Mat(𝒞) block matrix; a fixed-point combinator — DROP from the justification set).
DSA's R⊣P belongs OUTSIDE the ψ½ triangular algebra but is the DEFINING trigger for the
general coupled-ITERATIVE machinery when it lands (building it from ψ½ bakes in forward-subst
DSA violates). THREE couplings, THREE homes: ψ½=biproduct(exists), DSA=coupled-iterative(defer),
RQI=saddle-point(defer, different trigger).

## Q4 guard subsumption — YES, but via the SAME mechanism already scoped
Presence = fibration over the mesh (blocks(mesh) via R12a τ_raw∈(0,1)); "System B ∈
partition" ⟺ seed fiber non-trivial. Retires the 7 guards (`_require_starting_direction`
loss_representation/__init__.py:213 +3 callers, +4 refuse). BUT this IS the FaceField
design §5.3 mesh-enumerated block list (`PhaseSpaceCarrier`, C3 pending) — the coupled
frame supplies the NAME (fibration), NOT new power. Negative test: a composite whose block
set ≠ blocks(mesh) must be UNCONSTRUCTABLE (no runtime guard reachable).

## Pollination — the weld IS an un-named Schur elimination
A_BB cheaply invertible (triangular) ⇒ static condensation `A_AA − A_AB A_BB⁻¹ A_BA`.
"Welding A_BA into S" IS this Schur complement done implicitly — which is WHY it was
tempting, and it trades the explicit ψ½ state for an implicit condensation (the exact
"clean bulk vs fold-in" tension the FaceField design adjudicated toward EXPLICIT). Test:
welded-S bulk `== A_AA − A_AB A_BB⁻¹ A_BA` (array_equal) proves the weld = the condensation.

## Latent watch — metric-blindness on the coupled adjoint
`carlson_inward_sweep_transpose` is documented "EUCLIDEAN adjoint." The coupled seed off-
diagonal G-adjoint needs the G_sd=V_cell weighting (`G_sd⁻¹A_ABᵀG_A`). Not firing today
(gates fed zero seed) but the ERR-067 family — a nonzero-seed reciprocity gate on the
V_cell metric is the catcher.

## REFUTED / UNEXPLORED (this problem class)
Comma/slice category (no morphism beyond biproduct+row); homology (ι_*ι*=id + reflections
compose ≠0, no ∂²=0); operad/PROP (dagger-biproduct already IS it); MPO (⊕ fixed rank, no
bond knob); nonlinear-operator/fixed-point (RIGHT frame for multiphysics but NOT this linear
block machinery — flagged so it is not folded back in).
