---
name: issue-341-boundary-gs-rate
description: "#341 boundary G-S vs Jacobi SI rate — the mechanism (DD's d-1 undamped channels void Varga), the refutation of ndim as the predicate, and the standing recommendation (keep the default, own the octant ORDER)"
metadata:
  type: project
---

# #341 — boundary Gauss-Seidel vs Jacobi: mechanism found, `ndim` refuted

**Status 2026-08-09: investigated, NOT changed.** Full record with every configuration:
`scratch/issue_341_gs_jacobi_mechanism.md`; 8 probes at `scratch/probe_341_*.py`.
Behavioural lessons distilled to [[lessons-L19]] — read those first; this file holds only
what a session RESUMING #341 needs.

**Why:** the issue recorded the d=2/d=3 sweep-count inversion as measured-but-unexplained
and left the default open. **How to apply:** if asked to act on #341, start from §11 of the
record, not from the issue body's face-count story (refuted).

## The three rulings

1. **Mechanism.** Zero leakage ⟹ `B` is the whole iteration and its octant action is the
   hypercube `Q_d`. The multi-D DD face transmission is `Σ = (2/D)·1wᵀ − I` — one damped
   eigenvalue `1 − 2Σ_tV/D` plus **`d−1` eigenvalues exactly `−1`** (zero-cell-average
   sawtooth, invisible to absorption). Non-negativity fails ⟹ **Varga's comparison
   theorem no longer forbids `ρ_GS > ρ_J`**, and the fold pattern decides the sign.
2. ⛔ **`ndim` is NOT the predicate** — falsified both ways by direct spectral measurement
   (3 d=2 fixtures where G-S loses, up to 2.86×; 3 d=3 where it wins, down to 0.58×).
   Σ_t, mesh, aspect ratio, quadrature order and `c` all move the sign at fixed `ndim`;
   at `c ≥ 0.99` (the real lattice regime) the inversion is gone. **Do not branch a
   default on it.**
3. ⭐ **The dominant lever is unowned: the octant sweep ORDER.** `SweepSchedule.gauss_seidel`
   inherits `Quadrature.octants`' lexicographic enumeration. The fold is exactly
   `L_a` = the constant-sign suffix run of that order; all 25 achievable d=3 patterns were
   measured and separate exactly on **`LOSES ⟺ max_a L_a > Σ_{b≠a} L_b`**, spanning
   `n_GS/n_J ∈ [0.771, 1.892]`. The shipped `(4,2,1)` is 24th of 25. ⚠ the law is exact on
   ONE d=3 fixture and does NOT transfer to d=2 — sweep before encoding it anywhere.

## The recommendation as it stands (not landed)

* **Keep `solve_sn_fixed_source`'s `"gauss_seidel"` default; add no dimension branch.**
  G-S's leverage and its pathology occupy the same zero-leakage corner; in that corner the
  measured upside (6.3× at d=2 (10,10)) exceeds the downside (1.7× shipped-order d=3).
  `solve_sn` already defaults to `"jacobi"`, so production eigenvalue exposure is nil.
* **Two present-tense-FALSE docstring claims to fix** (the real deliverable):
  `sweep_schedule.py` module docstring's `ρ_GS ≈ ρ_J²` (wrong in both directions —
  measured 0.9064 vs 0.9286 at d=2, 0.9855 vs 0.9514 at d=3), and `_select_si_splitting`'s
  **"the regular splitting"** (Varga-regular requires `M⁻¹ ≥ 0, N ≥ 0`; DD supplies `d−1`
  channels at `−1`, and that failure is exactly what licenses the inversion).
* **Own the octant order** — its own `module:sn` / `type:improvement` issue.

## Two findings worth their own tickets

* **`A = L+C−S−B` is EXACTLY singular on an all-reflective box** — a trace-only,
  zero-cell-average null space. `ker(G−I)` is annihilated by the increment, so the SI
  residual stop cannot see it and the returned state keeps a frozen null component (the
  11.26 % trace deviation of `scratch/d3_absorber_diagnosis.md` §4). Benign for cell
  averages / scalar flux; nothing warns.
* **Two promotable probes** (`scratch/issue_341_..._mechanism.md` §13):
  `probe_341_iteration_spectrum.py` → a `slow` SN rate gate that measures ρ WITHOUT the
  stopping test (immune to the #340 truncation class); `probe_341_octant_model.py` → a
  sub-second `foundation` gate on the two model theorems. The survey probes are
  print-only and must NOT be promoted.
