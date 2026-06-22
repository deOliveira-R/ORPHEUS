---
name: issue-240-d5b-s3-diffusion-limit
description: #240 D5b-S3 thick-cell LD diffusion-limit FIXED (ERR-061). Root cause = the per-ordinate LD slope moment ψ̂_n was stored in the per-ordinate SWEEP frame (downstream-positive), but the angular reduction φ̂=Σ_n w_n ψ̂_n (feeding the isotropic scattering slope source Σ_s·φ̂) assumed the GLOBAL-x frame. For backward ordinates (μ<0) the sweep frame is sign-flipped, so forward+backward slopes partially CANCEL → φ̂ ~6× under-driven → LD lost the diffusion limit (nx=4: 38.9% off DD). NOT a source-weighting bug, NOT the M-normalization, NOT the dense d=1 assembly — all individually correct. Fix = octant_moment_frame_signs involution (∏_a octant_sign_a^{o_a}) applied via _reframe at the _CellSolve/_CellResidual octant boundary (source/probe global→sweep IN, moment/residual sweep→global OUT; outgoing FACE stays sweep-frame). DD byte-identical (per_axis=1 → None). Post-fix nx=4 rel 38.9%→4.1%.
metadata:
  type: project
---

# #240 D5b-S3 — thick-cell LD diffusion limit (ERR-061, FIXED 2026-06-17)

**Branch** `feature/sn-space-angle-tier2`. **NOT committed** (main agent commits).
Host env `.venv/bin/python`; canonical `python -O -m pytest`; NEVER all `tests/sn` (#212).

## ROOT CAUSE — the LD slope moment was stored in the SWEEP frame, not the GLOBAL frame

The per-cell LD kernel produces/consumes the `2^d` moment vector in the
per-ordinate SWEEP frame (each axis oriented so the downstream face is at local
`+1`). The iterate `φ̂` + its isotropic scattering source `Σ_s·φ̂` live in the
GLOBAL-x frame, where `φ̂ = Σ_n w_n ψ̂_n` sums the slope across ordinates of BOTH
sweep directions. For an ordinate sweeping in the NEGATIVE global direction on
axis `a` (`octant_signs[a] == -1`) the sweep coordinate is the REVERSE of the
global coordinate, so the slope (P₁) moment is sign-FLIPPED between frames. The
producer (`_CellSolve.cell` emit) stored the raw sweep-frame slope; the consumer
(`integrate_angular` / `S.apply`) summed it as global-frame → backward ordinates'
opposite-signed slopes CANCELLED the forward ones → `φ̂` ~6× under-driven →
the scattering slope source `Σ_s·φ̂` too weak → LD lost the diffusion limit.

**Smoking gun:** at a cell with positive global-x gradient, forward ordinate
ψ̂_n = +0.048 but backward ψ̂_n = −0.028 (opposite signs; should both be +). The
angular sum cancels to a tiny symmetric φ̂ while the physical gradient is
anti-symmetric.

**Failure mode #1 (sign flip) + #6 (convention drift, sweep-frame producer vs
global-frame consumer).** ERR-061.

## WHAT WAS RULED OUT (all individually CORRECT — the brief's suspects)

- **Σ_s⊗I + M-normalization** (the PRIME suspect): the M-normed `(L+C)ψ` slope
  row IS the raw ÷V slope balance; subtracting raw `S.apply(ψ)` is consistent
  (`A·ψ = M·S·ψ + M·q` — the scattering IS effectively mass-weighted). NOT a
  units bug.
- **The dense d=1 UBLD assembly** (the SECOND suspect): the dense `A` slope row
  is the analytic LM-1989 2×2 slope row × θ (the Galerkin mass-weighting — CORRECT);
  the dense per-cell SOLVE matches the analytic 2×2 exactly. NOT an assembly bug.
- **The scattering slope source magnitude:** `S.apply` produces `Σ_s·φ̂/W` at
  full strength on the slope row (ratio 1.0). CORRECT.
- **The matvec round-trip / SI≡Krylov (GATE 3):** vanishes to 1e-16 — but this
  is the TRAP. It proves SI and Krylov solve the SAME operator; it CANNOT prove
  that operator is diffusion-consistent (vv §5). The frame bug is a WRONG FIXED
  POINT that self-consistency is structurally blind to.

## THE DECISIVE EVIDENCE — structural independence (vv §1/§6)

A from-scratch LD-SN slab solver (direct LM-1989 Eq 4.3 cell 2×2 + source
iteration, NO ORPHEUS kernel) reproduced ORPHEUS's WRONG value BIT-FOR-BIT
(nx=4: 1.4717 both) when it summed sweep-frame slopes, and RECOVERED the
analytic diffusion VALUE (nx=4 rel 2.3%) when it stored global-frame slopes.
This pinned the root cause to the slope-moment frame INDEPENDENT of ORPHEUS's
code — the agreement between ORPHEUS and the independent solver proved they
shared the SAME bug, and the independent fix proved the fix class. (Also
confirmed against LM-1989 §IV Eq 4.9b discrete diffusion continuity, which the
pre-fix FP violated, and Eq 4.10c φ̂⁰=½(Φ_{j+½}−Φ_{j−½}) — the slope IS the
edge-gradient.)

## THE FIX (single-sourced involution; DD byte-identical)

`octant_moment_frame_signs(octant_signs, per_axis)` in `orpheus/sn/spatial/_ubld.py`
= `∏_a (octant_sign_a)^{o_a}` over the tensor-Legendre Kronecker layout
`[o_0,…,o_{d-1}]` (average sign-invariant; per-axis slope flips once if that
axis sweeps backward; d=2 cross `x̂y` flips when an ODD number of its axes
reverse). The map is an INVOLUTION (its own inverse).

Applied via the `_reframe` helper (`orpheus/sn/sweep_graph.py`) at the octant
boundary in BOTH cell ops:
- **`_CellSolve.cell`:** source `Q_cells` global→sweep IN; emitted `psi_avg`
  sweep→global OUT.
- **`_CellResidual.cell`:** probe + source global→sweep IN; residual sweep→global OUT.
- **The d=1 matvec** (`loss_representation.py` `_apply_walk._sweep_direction`):
  per-direction frame signs on the probe + residual.
- The OUTGOING FACE (`psi_out`) stays SWEEP-frame — it propagates along the
  wavefront, never crosses into the global-frame iterate.
- `_reframe` guards on `arr.shape[-1] != frame_signs.shape[0]` → a FLAT scalar
  source (matvec zero / flat external — only the sign-invariant average moment)
  passes untouched (NOT broadcast into a spurious moment axis). DD/Step
  (`per_axis==1` → `None`) byte-identical.

**Files:** `orpheus/sn/spatial/_ubld.py` (`octant_moment_frame_signs` — returns
`None` for per_axis==1, the single-source no-op convention), `orpheus/sn/sweep_graph.py`
(`_reframe` + `_CellSolve`/`_CellResidual` `moment_frame_signs` field),
`orpheus/sn/loss_representation.py` (`_moment_frame_signs` thin delegate + 4 d≥2
construction sites + the d=1 `_apply_walk` calling the primitive directly).

**Elegance review (elegance-enforcer) = PASS-WITH-NITS, both nits ADDRESSED:**
- DO-NOW (the twin path): the d=1 `_frame_signs` closure was a smaller copy of
  the per_axis→None decision `_moment_frame_signs` owns → COLLAPSED. Both now
  call `octant_moment_frame_signs((signs), per_axis)` directly; the primitive
  owns the `None`-for-DD convention (Pattern 2 — single source).
- RECORD (the typed-predicate debt): `_reframe`'s `arr.shape[-1] != frame_signs.shape[0]`
  is the 4th value-based "does this carry the spatial-moment axis?" shape probe
  (with `_n_face_moments!=1`, the retired `len(s_axes)>1`, the pure_z
  `q.ndim>sig.ndim+1`). Filed **#246** (replace all 4 with one typed
  `SpatialMomentSpace` predicate; S4 trigger = `n_diag==2^d` collision on a 4×4
  d=2 mesh — the shape probe could silently mis-fire a sign flip once the
  non-vacuum moment trace lands).

## RESULTS (L12 paste-back)

```
# DIFFUSION-LIMIT REPRO — pre-fix → post-fix
            pre-fix                post-fix
nx=  4   DD=2.4069 LD=1.4717 rel=0.389   →   LD=2.3080 rel=0.041
nx= 16   DD=2.3608 LD=2.1752 rel=0.079   →   LD=2.3559 rel=0.002
nx= 64   DD=2.3597 LD=2.3415 rel=0.008   →   LD=2.3594 rel=0.000
# 2-D analog: rel 8.4% → 1.7% → 0.4% across n=4/8/16 (vacuum on all 4 sides;
#   the residual 8.4% at n=4 is the 2-D boundary-layer fraction, not a bug)
```
GATE 4 (DD/Step byte-id, `tests/sn/sweep/core tests/sn/solve -W error::DriftWarning`)
= **513 passed / 1 skipped / 4 xfailed** — IDENTICAL pre/post (zero drift).
GATE 3 (2-D LD MMS: O(h²) + FFW≡MFW + DD≠LD) = 3 passed.
LD foundation + primitive = 59 passed; spatial = 68; transport+numerics = 848;
operators = 420 (excl. the 7 documented pre-existing curvilinear reds).

## GATES (promoted, mutation-verified)

- `tests/sn/verification/mms/test_mms_ld_slab.py::test_ld_thick_diffusive_limit`
  (1G, flipped xfail→PASS + Mode-8 migrate to `np.testing.assert_array_less`) +
  `::test_ld_thick_diffusive_limit_2g` (2G-het, Mode-6 group-coupled slope).
- `tests/sn/spatial/test_ld_slope_frame.py::test_ld_slope_moment_global_frame_consistency`
  (frame consistency) + `::test_independent_ld_global_frame_recovers_diffusion`
  (`@foundation`, structural-independence ground — NO catches marker, it does
  not exercise `_reframe`).
- All `catches("ERR-061")` markers MUTATION-VERIFIED: re-dropping `_reframe`
  (return arr) turns 3 of them RED under `-O` (the independent one stays GREEN
  by design — it is the independent ground, not an ORPHEUS-path catcher).

## OWED-2 (the d=1 CumprodScan SI path) — STILL OPEN, my fix DICTATES its form

The d=1 SI/scan path (`CumprodScan`) still crashes on the moment source
(broadcast `(16,1,20,2)` vs `(1,1,nx)`) — PRE-EXISTING (confirmed by git-stash;
the prior closeout's 3 scan-broadcast reds). My fix is on the Krylov/matvec
path (the diffusion-limit priority). **My fix DICTATES the scan's slope-source
form:** `scan_xV` must gain an `s_hat` slope-source arg AND the scan must apply
the SAME frame sign per sweep direction — for the backward scan the slope source
is `-Σ_s·φ̂/W` (global→sweep IN) and the reconstructed ψ̂ is sign-flipped back
(sweep→global OUT), exactly as `_reframe` does for d=1. Single-source the scan's
frame handling through `octant_moment_frame_signs((direction_sign,), 2)`.

## METHODOLOGY LESSONS (durable)

1. **A per-ordinate spatial-moment quantity produced in a direction-dependent
   SWEEP frame MUST be lifted to the global frame BEFORE the angular reduction
   that sums it across ordinates** — producer and consumer must agree on the
   frame or forward/backward ordinates cancel a quantity that should reinforce.
   The angular reduction is the discriminator (H2/H3: flat flux nulls the slope;
   round-trip/conservation are telescoping-degenerate to the frame error).
2. **Matvec-self-consistency (SI≡Krylov, round-trip≈0) is necessary but NEVER
   sufficient for a moment-iterate fold.** It proves the operator is internally
   consistent, not that its FP is physically correct. Gate the converged VALUE
   against a structurally-independent reference (continuous diffusion + an
   independent from-scratch kernel), never the round-trip. (The brief named this
   trap; it held.)
3. **The structural-independence ground that CRACKED it:** a from-scratch solver
   reproducing the wrong value bit-for-bit is not a coincidence — it localizes
   the bug to the SHARED math (the 2×2 + the frame convention), and the
   independent FIX (global-frame slope) confirms the fix class without touching
   the production code. Build the independent kernel when "all components are
   individually correct but the FP is wrong."
4. **When the literature says a scheme IS consistent but your faithful
   implementation isn't, the bug is a CONVENTION between two correct pieces, not
   in either piece.** Every LM-1989 Eq (4.3) cell-level check passed; the defect
   was the sweep↔global frame at the producer-consumer seam.
