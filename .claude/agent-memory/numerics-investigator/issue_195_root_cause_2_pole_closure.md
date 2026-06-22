---
name: issue-195-root-cause-2-pole-closure
description: Issue #195 — curvilinear-isotropic MMS plateau is a REGRESSION (Phase-D was O(h²)), root cause = the inward/outward pole-face DD seed reads the cell-CENTRE flux ψ(r=Δr/2) instead of the r=0 regularity value ψ(0), seeding an undamped odd-even sawtooth that poisons the conservative streaming domain-wide. Fix = one seed line in _compute_LpC (matvec, operator.py:410) + its sweep twin (loss_representation.py:2094).
type: project
---

# Issue #195 — VERDICT: pole-face DD seed regression (single line, twin sites)

**Date**: 2026-06-12 (round 2, term-decomposition + archaeology). **Branch**:
`main` @ `a7a67d8` (HOST, `.venv/bin/python`).
Diagnostics: `derivations/diagnostics/diag_195_probe{1,2,3,4,5}_*.py`.

## Headline (supersedes round-1's "no regressor, seed always there")

Round 1 concluded "root cause (2) pole closure, no regression." **Round 2
corrects this: there IS a regression, and the seed line is the defect.**

- **Phase-D head `7d1cdd8` (Carlson era): O(h²)** — errs `[6.393, 0.6337, 0.1149]`,
  orders `[3.335, 2.463]` (SI ≡ Krylov, bit-identical). Reproduces #195's claim.
- **Phase-C `1ea1cad` (pre-Carlson): PLATEAU** `[0.083, 0.095, 0.098]` — the SI/sweep
  was always broken before Carlson.
- **Current `main`: PLATEAU** `[0.0376, 0.0396, 0.0409]`, orders ≈ 0 — REGRESSED
  back to broken. The Phase-D Carlson consistency was LOST in the unification.

⚠ **METHODOLOGY TRAP (cost me 2 contaminated runs):** ORPHEUS installs `orpheus`
as an EDITABLE package via a custom meta-path finder
(`__editable___orpheus_0_1_0_finder`) that hardcodes the MAIN checkout and
OVERRIDES `PYTHONPATH`. Worktree archaeology with `PYTHONPATH=$WT` silently runs
MAIN's code. **Correct isolation:** `cd $WT` then strip the editable finder from
`sys.meta_path` AND `sys.path.insert(0, $WT)` inside the script; assert
`orpheus.__file__.startswith($WT)` before trusting any result.

## Part 1 — THE NAMED TERM (the deliverable)

**The inconsistent term: the discrete spatial STREAMING, corrupted by an
undamped odd-even DD sawtooth seeded by the WRONG pole-face initial condition.**

### Term table (sphere isotropic, nc=160, mid cell 80, r=2.5156; matvec units /V)

Reconstruction matches the actual matvec bit-for-bit (recon_err = 4.4e-16).

| ordinate | term | discrete | continuous | mismatch |
|---|---|---|---|---|
| μ=+0.989 | streaming | +0.0794 | −0.0031 (μA'/W) | +0.0824 |
| | redistribution | −0.4544 | 0 (isotropic) | −0.4544 |
| | collision | +0.5000 | +0.5000 | 0 |
| | source q | +0.2469 | +0.2469 | 0 |
| | **NET (L·ψ−q)/V** | **−0.122** | 0 | −0.122 |

Continuous **curvature** `μ·2A/r/W = +0.393`. The DISCRETE streaming carries only
`+0.082` of it (**21%**); the DISCRETE redistribution removes `−0.454` (**116%**).
`CURVATURE_in_stream + REDIST = +0.082 − 0.454 = −0.372 ≠ 0` (μ>0); `+0.098` (μ<0)
— direction-ASYMMETRIC. Identical at nc=160 AND 320 (mesh-independent).

### Mechanism (the smoking gun, probe-4 face-flux trace)

At mid cell the DD face fluxes have drifted from the reference:
`ψ_ref(r_in)=0.5000` but `DD ψ_face_in=0.5049`; `ψ_ref(r_out)=0.49990` but
`DD ψ_face_out=0.49502`. With the REFERENCE faces, conservative streaming
`μ(A_out·ψ_out − A_in·ψ_in)/V = 0.3902 ≈ continuous 0.3902` (PERFECT). With the
DD faces, streaming = 0.0794 (21%). The DD face drift is a **constant-amplitude
±4.91e-03 odd-even SAWTOOTH** superimposed on the reference — same amplitude at
cell 0 and cell 80 (mesh-independent).

**Origin of the sawtooth:** the outward sweep's pole-face seed is
`pole_face_seed = psi_view[:, :, 0]` (operator.py:410) = the innermost CELL-CENTRE
flux `ψ(r=Δr/2)`, but the inner face is at r=0 where the regularity value is
`ψ(0)`. The seed error is the half-cell offset `ψ(Δr/2) − ψ(0) = 4.909e-03`
(= sin(π·Δr/2/R)/W for the sin ansatz). The DD streaming recurrence has
attenuation `a = 2|μ|A_total/denom − 1 = 0.945 ≈ 1` (cells optically thin,
`Σ_t·V/(2|μ|A_down)=0.016`), so the seed error propagates UNDAMPED as the
sawtooth. The conservative streaming `A_out·ψ_out − A_in·ψ_in` amplifies the
±ε face sawtooth by `~A/V ≈ 32×` (A_in≈A_out≈80, nearly equal, both huge) into
an O(0.3) streaming error per ordinate — domain-wide.

Geometry coefficients are PERFECT: `A_face = 4πr²` exact, `V = 4π/3(r_out³−r_in³)`
exact, `(A_out−A_in)/V = 0.79502 ≈ 2/r = 0.79503` (4 digits). The MMS source is
correct (non-conservative `μA'/W`, consistent with the continuum). The
discretisation IS consistent. **Only the pole-face seed is wrong.**

### Sanity anchors reconciled

- **(a) sign-antisymmetric residual (−0.158 mid vs +0.158 outer):** the sawtooth
  is odd-even (`[-0.15,+0.15,-0.15,...]`), so the angle-integrated residual
  alternates sign cell-to-cell; the |res| envelope is ~constant, the sign flips.
- **(b) solution err PLATEAUS at 0.04 while residual is 0.24:** the operator has
  a clean fixed point (Krylov→2e-14) that is the WRONG solution; the wrong-FP
  error magnitude is set by the residual but the L2 NORM of the solution error
  is set by how the wrong-FP spreads (pole-concentrated + bulk sawtooth).
- **(c) pole ±55 per-ordinate vs continuous ±0.31:** the pole cell has A_in=0 but
  the sawtooth seed enters there first with full amplitude; `dA/w~1/r` amplifies.

### The fix (probe-5, VERIFIED on the real apply path, sphere AND cylinder)

Replace the pole-face seed with the WDD-consistent r=0 extrapolation:
`ψ(0) ≈ 1.5·ψ_cell[0] − 0.5·ψ_cell[1]`, leaving cell VALUES untouched.

Operator residual `(L+C−S−B)·ψ_ref − q`, scalar max|res|:
```
seed = ψ_cell[0]  (PRODUCTION/BUG):  0.236  0.236  0.236   (nc 40/80/160, sphere)
seed = 0  (= ψ(0) for sin ansatz):   3e-15  5e-15  1e-14   MACHINE ZERO
seed = 1.5·c0 − 0.5·c1 (general):    7.3e-4 1.8e-4 4.6e-5  O(h²) (ratio ~4)
```
Cylinder real-apply-path with the patch: `4.6e-4 → 1.2e-4 → 2.9e-5` (O(h²)).
The zero-seed gives machine zero because sin(0)=0 hits ψ(0) exactly; the
extrapolation is the production-ready general form (works for non-vacuum a0>0).

⚠ **Probe-5's hand-rolled CYLINDER reconstruction is INCOMPLETE** — its sweep loop
iterates only ±μ-masked ordinates and SKIPS the degenerate η≈0 (pure-azimuthal)
ordinates, which carry 25% of the cylinder weight. The REAL operator handles
their collision (`LC.apply` gives Σ_t·ψ on them, verified). So probe-5 asserts
on the SPHERE; the cylinder fix is verified via the real apply path (the cylinder
probe-3 residual is a clean ±0.15 odd-even sawtooth — same pole-seed mechanism).

## Part 2 — ARCHAEOLOGY VERDICT

**The Carlson coupled-pole consistency (Phase D) was lost in the curvilinear
matvec rebuild, regressing the bulk streaming sawtooth back.** Bisect (sphere
iso MMS SI ladder nx 20/40/80, GOOD = O(h²), BAD = plateau, properly isolated):

```
7d1cdd8  Phase-D head            GOOD  [6.39, 0.63, 0.11]  orders [3.34, 2.46]
2f24ae4  apply→unified matvec    GOOD  [6.39, 0.63, 0.11]   (NOT the regressor, as orchestrator said)
975edc5  retire legacy curv mv   GOOD
b667f8e  R-1 Step 0 pole hist    GOOD
c93355c  unify sweep+matvec      GOOD
c4fa030  R-1 Step 4 G0 native    GOOD  ← LAST KNOWN GOOD
... (124-commit Depth-B D-H / Wave-T window) ...
90e7d4e  Wave T T.4c curv leaf   BAD   [0.0376, 0.0396, 0.0409]
ad813fd  delete _matvec_unified  BAD
3a79ab3  walk → LossRepresent.   BAD   ← same plateau as main
main     a7a67d8                 BAD
```

**Regression window: `c4fa030`..`90e7d4e`** (the Depth-B D-H typed-field migration
+ Wave-T `M_spatial`/`M_angular_redist` leaf carve that rebuilt the curvilinear
matvec). The exact single commit was not pinned (Depth-B intermediates have API
churn that breaks `solve_sn_fixed_source`; worktree-add friction); the brief
time-boxed Part 2 and Part-1 already names the fix. **Phase-D's GOOD state had
the IDENTICAL seed line** (`psi_face_in = fi[:, outgoing_mask, i0, 0]`, operator.py
Phase-D) — so the seed line is not literally NEW; what changed is the
redistribution/closure that USED to cancel the sawtooth in the bulk. At Phase-D
the bulk was CLEAN (mid err 0.0015 at nx=80, only the pole spiked, L2=0.115);
on main the bulk carries the sawtooth too. The cleanest framing: the Carlson
M-M redistribution at Phase-D produced a half-angle flux that compensated the
streaming sawtooth; the rebuilt curvilinear `precompute_psi_state` /
`cell_contribution` on main does NOT (main's M-M ψ_ang_in = 0.5777 vs iso 0.5000
at mid cell — a 15.5% error, the compensation gone). **The surgical fix is the
seed line regardless** (it makes the operator admit ψ_ref independent of the
redistribution path, probe-5 proves O(h²) collapse).

## Part 3 — #206 CONNECT/DISCONNECT: DISCONNECTED

`tests/sn/sweep/curvilinear/test_unified_matvec_cylinder.py::test_unified_cylinder_matches_hand_reference`
is `xfail` (27 cases). **#206 is a DIFFERENT bug.** Its hand reference
(`_hand_reference_cyl_matvec`, line 183) seeds the outward pole with
`psi_face_in = psi_g_first[:, n_g, 0]` — the SAME cell-centre seed as production.
So both sides of #206's comparison share the (potential) seed choice; it cancels
in their difference. #206's ~18% divergence is the bool-mask-scatter / per-ordinate
routing mismatch (ERR-049 / SIG1, the `legacy_proxy_matvec` vs per-ordinate hand
ref), NOT the #195 pole-seed. The two are independent.

## Part 4 — FIX PROPOSAL

**What:** the pole-face DD seed for the OUTWARD (μ>0) curvilinear sweep must be the
r=0 regularity value, not the innermost cell-centre value.
**Where (twin sites — Pattern 2, fix BOTH):**
1. matvec: `orpheus/sn/operator.py:410` — `pole_face_seed = psi_view[:, :, 0].copy()`
   → `pole_face_seed = (1.5 * psi_view[:, :, 0] - 0.5 * psi_view[:, :, 1]).copy()`
2. SI sweep: `orpheus/sn/loss_representation.py:2094` —
   `psi_in = ig_values[global_n, :, 0]`
   → `psi_in = 1.5 * ig_values[global_n, :, 0] - 0.5 * ig_values[global_n, :, 1]`
   (when `ig_values is not None and nx >= 2`; cold-start `np.zeros(ng)` stays).
**Better (elegance Pattern 7 + 2):** factor the pole-face extrapolation into ONE
named helper (`pole_face_regularity_seed(psi, sn_mesh) -> (N, ng)`) consumed by
BOTH sites — the seed is a PRODUCER concern (the r=0 regularity closure), not a
per-consumer re-decision. This also lets the slab path (`face_inner`) and the
curvilinear path share the "pole/inner seed" abstraction.

**Caveat to verify before landing:** the extrapolation `1.5·c0 − 0.5·c1` is
O(h²)-accurate but uses cells 0,1; for a 1-cell mesh it degenerates. Also confirm
it does not over/under-shoot for the anisotropic cases (B(r) with B(0)=0 still
extrapolates cleanly). Consider whether the EXACT pole-face value is recoverable
from the M-M half-angle closure (the Carlson seed already computes a pole-face
half-angle flux — the seed might be derivable from `psi_state` consistently,
which is likely what Phase-D effectively did).

**Gates that prove the fix:**
- probe-1 ladder collapse: sphere/cyl iso L2 → O(h²) (orders > 1.9).
- probe-3 residual → O(h²)-consistent (max|res| shrinks ~4× per halving), NOT
  the 0.236/0.150 plateau.
- probe-5 (this round): extrap-seed operator residual O(h²) (already GREEN as a
  standalone; promote with the real source fix in place to flip it to a true
  regression gate).
- the 6 xfail markers (2 iso + 2 aniso sphere/cyl + the magnitude variants) flip
  to xpass.
- flat-flux gates stay GREEN (`test_per_ordinate_flat_flux_consistency`): the
  extrapolation on a FLAT field gives `1.5·c − 0.5·c = c` = the cell value, so
  flat-flux behaviour is UNCHANGED (the seed fix is a no-op on flat ψ — explains
  why the flat-flux gate never caught this; Mode 7 + H2).
- Cartesian bit-identity: the slab path is untouched (curvature == "cartesian"
  branch unchanged), so all slab snapshots stay bit-identical.

ERR-026 stays PARTIAL until the source fix lands (flat-flux identity closed since
Phase D; non-flat pole-seed consistency is THIS fix). Recommend a new ERR-NNN
manifestation entry (failure mode #3 missing-factor / #5 index — the seed reads
cell-0 value as a face value, a half-cell index/quantity confusion) once landed.

## Methodology lessons (round 2)

1. **Worktree archaeology under an editable-install finder is a trap.** Verify
   `orpheus.__file__.startswith($WT)` BEFORE trusting any era result; strip the
   editable meta-path finder. My first 2 Phase-C/D runs ran MAIN and falsely
   "confirmed no regression."
2. **Re-measure the issue's premise WITH proper isolation.** Round 1 (this same
   trap) wrongly concluded "no regressor." The O(h²) Phase-D claim was REAL.
3. **Term decomposition must reconstruct the ACTUAL operator (recon_err < 1e-15),
   not an idealised hand version.** My first probe-4 used a clean DD telescoping
   + isotropic ψ_ang idealisation and gave WRONG per-ordinate values; mirroring
   `_compute_LpC` exactly (imposed cell-avg + M-M `precompute_psi_state`) made it
   bit-identical and the decomposition trustworthy.
4. **The decisive MMS-plateau probe is "swap the suspect seed, watch the residual
   collapse."** Probe-5's `seed_mode ∈ {cell0, zero, extrap}` pinned the single
   line: cell0 → 0.236, zero → 1e-14, extrap → O(h²). No solve needed.
5. **Hand-rolled sweep reconstructions miss degenerate ordinates.** The cylinder's
   spurious −0.25 was my probe skipping η≈0 ordinates (25% weight); the real
   operator was fine. Always cross-check the hand sweep against `LC.apply` on a
   degenerate ordinate before believing a geometry-specific extra defect.
6. **The flat-flux gate is blind by construction:** `1.5·c − 0.5·c = c`, so the
   seed fix is a no-op on flat ψ — which is exactly why the bug hid (Mode 7 + H2).
