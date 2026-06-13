---
name: curvilinear-tau-clamp-vs-pole-floor
description: τ-clamp is NOT the dominant curvilinear-SN error floor — the sphere pole-CELL spatial closure (first-order at r→0) is; cylinder differs (clamp severe but unfixable by partition). REFINED 2026-06-13b — pole O(h) is largely an MMS comparison artifact + literature-accepted INHERENT first-order; NOT cleanly fixable by a local closure; INVISIBLE to the production volume-weighted L2 gate.
metadata:
  type: project
---

## ⭐ REFINED VERDICT (2026-06-13b — pole-cell FIXABILITY study, supersedes the §"Recommended approach Priority 1" below)

The prior §"Recommended approach" called the sphere pole O(h) a "high-leverage
OPEN first-order defect — file an issue + fix the closure." **Deeper probing
DOWNGRADES this to WONTFIX/document-only.** The pole O(h) decomposes into THREE
parts, none of which warrants a code fix:

1. **~75% MMS COMPARISON ARTIFACT (not a solver bug at all).** The production
   spherical MMS evaluates the source at the MIDPOINT `mesh.centers` and compares
   `phi_solver` against `phi_exact(midpoint)`. But the spherical DD discrete
   unknown IS the cell-VOLUME-average `(4π/V)∫r²φ dr` (Hébert 2009 Eq. 3.430 —
   the unknown is DEFINED as the shell average, not a point value; diamond
   Eq. 3.431 relates it to the FACE fluxes). Under r²dr weighting the
   volume-average and the midpoint point-value differ by O(h) at the pole cell
   (r_lo=0 maximally non-uniform weight; volume-centroid at ¾h not ½h). Using
   the SHELL-AVERAGED source (Hébert 3.430) AND comparing vs the shell-volume-
   average drops the pole err ~4× (0.0212→0.00497) — diag_18.

2. **~25% GENUINE first-order solver residual — but LITERATURE-ACCEPTED INHERENT,
   NOT a bug.** Even the fully-consistent FV-MMS (shell-avg source + shell-avg
   ref) leaves the pole at clean O(h^1.00). ROOT CAUSE (diag_20/28/29): the
   diamond spatial closure `ψ̄=½(ψ_in+ψ_out)` over-predicts the pole OUTER FACE
   flux by EXACTLY +50% (mesh-independent rel_err=0.5000; diag_20) — because at
   r_lo=0, A(0)=0 so diamond gives `ψ_out=2ψ̄`, while the true face is `A(h)`
   and `2·⟨A⟩_vol = 2·¾A(h) = 1.5·A(h)`. DEEPER: the conservative BALANCE itself
   (E1) is inconsistent at the pole — fed EXACT ψ̄ + EXACT inflow it solves for an
   outer face −46% wrong (diag_28), residual/V plateaus at ~0.08 mesh-independent
   (diag_29). This is because A_in=0 degenerates the streaming surface integral
   while V~h³. **Hébert §3.9.4 + Stacey §9.9 BOTH use exactly this plain-diamond
   + Carlson-starting-direction + symmetry scheme at the central cell with NO
   special O(h²) closure, and NEITHER flags reduced order there** (literature-
   researcher: negative Tier-2 search too). First-order at the SINGLE pole cell is
   the accepted, unflagged behavior of the standard scheme.

3. **NOT cleanly fixable by a local closure.** Tested the volume-weighted linear
   reconstruction `ψ̄=β·ψ_out+(1−β)·ψ_in`, β=¾ at pole (the value that makes ψ̄
   O(h³)-consistent vs ⟨A⟩_vol at EXACT faces, diag_21). FAITHFUL end-to-end
   validation (production sweep+cache monkeypatched, β=½-identity verified to
   3e-16, diag_25/26): β=¾ does **NOT** restore O(h²) — pole stays O(h),
   magnitude slightly WORSE (0.0050→0.0106), and full-mesh β degrades the
   interior. Closure-consistency-at-exact-faces ≠ fixed-point accuracy (the
   propagated face flux couples back through the balance). A genuine fix would
   require a non-local higher-order central-cell reconstruction the canon does
   not provide.

**INVISIBLE TO THE PRODUCTION GATE + NEGLIGIBLE DOWNSTREAM.** The production
`test_sn_spherical_mms_converges_second_order` uses a VOLUME-WEIGHTED L2 norm
(`√Σ V·diff²`); the pole O(h) at ONE cell of V~h³ contributes √V~h^1.5 →
h^2.5 to L2 — subdominant. Both midpoint AND volume-avg L2 references converge
clean O(h².00) (diag_30); only the L∞ (pole) is O(h). keff: reflective sphere
= k_inf=1.875 EXACT mesh-independent; vacuum sphere keff converges monotone to
~1.78590 at O(h^1.48) (combined pole+outer-face first-order; increments 2e-5 at
nx=160 — far below engineering tol; diag_31).

**CYLINDER (task 4 — IDENTICAL underlying defect, MASKED differently).** Cyl pole
vs MIDPOINT is O(h²) (1.94→1.98, matches old diag_12) BUT vs VOLUME-AVERAGE is
O(h) (0.99). Same diamond inconsistency, but cyl's r dr (linear) weight puts the
volume-centroid at ⅔h while diamond `½A(h)` ≈ midpoint `A(h/2)` — the midpoint
comparison the test uses is accidentally O(h²) for the cylinder. So the cyl pole
is NOT "clean O(h²)" (the old §below overstated it); it is the same O(h)
volume-average defect, masked by the midpoint comparison + linear weight. Cyl
global L2 also O(h².00).

**RECOMMENDATION: WONTFIX the closure + DOCUMENT + add an L∞/pole characterization
gate (xfail-or-band, NOT a fix gate).** No `orpheus/` change. The local closure
does not fix it; the global gate is green; downstream keff impact is negligible;
the literature accepts first-order at the single pole cell. Optional principled
hygiene (separate, low value): switch the spherical MMS SOURCE to shell-averaged
(Hébert 3.430) and the reference to the shell-volume-average — removes the ~75%
comparison artifact and makes the MMS self-consistent with the FV unknown, but
does NOT change the solver and the gate already passes. METHODOLOGY: the decisive
probes were (i) E_test = E_artifact(midpoint−volavg, O(h)) + E_true(solver−volavg)
decomposition (diag_16); (ii) the discrete-balance residual fed EXACT fields
(diag_28/29) isolating the inconsistency to the BALANCE not just the closure;
(iii) the FAITHFUL production-sweep monkeypatch with a β=½-identity regression
guard (diag_25/26) PROVING the local-closure fix fails end-to-end (the offline
closure-consistency win in diag_21 did NOT transfer). Scripts diag_16..31 in job
tmp `84fd66f8/tmp/negclosure`.


Offline-feasibility characterization (plan-mode, read-only) of negative
half-angle fluxes + τ-clamp in curvilinear (sphere/cyl) SN angular closure.
Premise studied: "the M-M τ-clamp `max(0.5,min(1,τ_raw))` is THE angular
error floor; can a positivity-preserving-AND-exact closure replace it?"

**VERDICT: the premise is correct OFFLINE but FALSIFIED end-to-end.** One
level deeper than the τ-clamp reveals the clamp is NOT the binding
constraint. Geometry-split:

**SPHERE.** Dominant solution error = the **pole-CELL spatial closure**,
not the τ-clamp. Decisive controls:
- diag_11 ISOTROPIC MMS (τ-clamp provably SILENT — flat-in-μ threads
  exactly clamped OR unclamped) shows the SAME pole floor as aniso
  (iso 3.8e-2→5.9e-3, aniso 3.3e-2→5.0e-3, nx 20→160).
- diag_14: pole cell (r→0) converges **O(h^0.91→0.97)=first-order** and is
  100% of total error; interior O(h^2.00); outer O(h^1.96) (the #168
  outer-face fix HOLDS). So the floor is a first-order r→0 pole-cell
  truncation, an OPEN undocumented defect (≠ #168 outer-face CLOSED,
  ≠ #229 angular-thread, ≠ τ-clamp).
- diag_07/diag_10: unclamping the sphere end-to-end improves err only
  ~3% (3.32e-2→3.23e-2); the residual ~6.3e-3 PLATEAUS in quadrature
  order (n_ord 8→64 unchanged) AND in clamp — neither knob touches it.
- diag_15: unclamped sphere SI CONVERGES on every stress config (thick
  abs, near-vacuum, c=0.999, S64), positive φ, no NaN; clamp costs a few
  SI iters (39→61 low-scatter) but is dispensable for STABILITY.
τ_raw range sphere GL = [0.39,0.61] (never 0); amplification (1−τ)/τ ≤ 1.56
unclamped, clamp bounds it ≤1 at τ≥0.5.

**CYLINDER.** OPPOSITE structure. Pole closure is FINE: diag_12 ISOTROPIC
cyl converges clean O(h²) (8.6e-4→1.5e-5, ratio 3.9). The aniso floor IS
angular-thread (#229: scales ~2× per n_mu doubling) and the τ-clamp is
SEVERE — but UNFIXABLE by unclamping or partition tweaks:
- diag_02: product/LS quadratures put the most-inward azimuthal ordinate
  `eta[0]` EXACTLY on `-sin_theta` ⟹ `τ_raw[0]=0` EXACTLY ⟹ unclamped
  recurrence divides by 0 (NaN). Not "near zero" — bit-exactly zero.
- diag_13: product quads have DUPLICATE azimuthal η (±φ symmetry pairs);
  the M-M thread marches in η only, so a linear-in-η field is NOT
  threadable exactly (structural mismatch: full angular variation is
  (η,φ) but the thread sees η alone). No partition (midpoint/cumweight/
  ord-interior) gives τ_raw∈[0.5,1] with bounded edges; cumweight is exact
  on LS but needs τ_raw∈[−4.5,5.5] (edges outside the level). α-dome
  telescoping is partition-INDEPENDENT (α=−cumsum(w·η)) so balance is
  always preserved; partition is a pure accuracy/positivity knob.

**NEGATIVITY (Q1).** diag_03: on EVERY realistic converged solve (smooth
MMS, homog eigenvalue keff=1.0, thick absorber) ZERO negative half-fluxes,
clamped or unclamped. EVERY clamp activation is SPURIOUS (100%: 160/320/
80/240 activations, 0 protective). diag_06: half-flux negativity that DOES
appear is in EARLY SI iterates and is INHERITED from a negative INPUT ψ
(early-iterate transient, psi_min=−0.107 at iter 1→+0.028 converged); the
clamp barely reduces it (88→78) and never eliminates it. diag_05: on
RANDOM rough positive ψ the unclamped recurrence goes negative ~99.8% (to
−12) because (1−τ)/τ>1 over-subtracts — but such rough fields do NOT arise
in converged SMOOTH iterates. So the clamp's positivity value is real only
against pathological roughness that the physics never produces; on physical
fields it is a pure accuracy tax.

**FEASIBILITY (Q5).** A positivity-preserving-AND-exact angular closure is
feasible+contained for the SPHERE (drop/relax the τ floor; SI stays
stable) — but it would NOT move the needle because the pole-cell spatial
error dominates. The HIGH-LEVERAGE fix is the **sphere pole-cell first-order
spatial closure**, not the τ-clamp. The cylinder needs a deeper angular
redesign (the η-thread can't represent (η,φ) variation; partition tweaks
fail), tracked by #229.

**METHODOLOGY LESSON (the durable one):** an OFFLINE-isolated error (the
user's proven τ-clamp thread error 1e-3→1e-15) can be a SECOND-order
contributor masked by a LARGER end-to-end error. Always run the
end-to-end swap (diag_07) AND a clamp-SILENT control (isotropic MMS,
diag_11) before accepting "X is THE floor." The isotropic-MMS control is
the structurally-independent ground that silences the suspected mechanism
by construction — it converted "τ-clamp is the floor" into "pole closure
is the floor" in one run. (vv H4 convergence-rate≠value + the
attribution discipline of L4/L5.)

Scripts: `/Users/rodrigo/.claude/jobs/84fd66f8/tmp/negclosure/diag_01..15`.
Promotion candidates: diag_11 (iso-vs-aniso pole floor control), diag_14
(per-cell rate pin) → `tests/sn/` as the OPEN pole-cell-O(h) regression
gate once a fix lands. Relates to #229 (OPEN, aniso quad retune), ERR-026
(PARTIAL), #195/#196 (CLOSED — fixed the seed/fixed-point, NOT the
pole-cell truncation), #168 (CLOSED, outer-face, distinct).
