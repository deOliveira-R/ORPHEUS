# Curvilinear aniso SN: root-cause program (#229 + #9 + new pole-cell defect)

## ⭐ STATUS 2026-06-13 — ALL 4 INVESTIGATIONS COMPLETE; HONEST REFRAME

The investigations (literature ×2, numerics ×1, test-architect ×1) revealed
**there is no solver bug to fix** — both candidate "defects" are inherent,
literature-accepted discretization properties. The deliverable is rigorous
understanding + documentation + #9 real coverage + honest #229 retune, NOT a
solver rewrite. Reconciled facts:

- **W0 ✅ DONE + committed `6c712ee`** (aniso-MMS symbolic guard 4D→3D).
- **W2 (pole-cell) = WONTFIX + document + L∞ characterization gate. ZERO code.**
  numerics verdict: the sphere pole-cell O(h) is ~75% an **MMS comparison
  artifact** (DD unknown is the r²dr shell-volume-average per Hébert 3.430; gate
  compares vs `phi_exact(midpoint)`) + ~25% genuine-but-**inherent** first-order
  residual (diamond β=½ at the A(0)=0 singular face; β=¾ closure does NOT restore
  O(h²) end-to-end — faithfully validated). Production gate GREEN (global L2
  O(h²), pole √V→h^2.5 diluted); keff impact negligible. **Cylinder shares the
  same defect, MASKED by the midpoint comparison** (not sphere-asymmetric).
  Literature (Hébert §3.9.4, Stacey §9.9): canon uses plain diamond at the
  central cell, no special O(h²) closure; O(h²)-at-origin needs a nodal/LD/SC
  scheme (Wu-Xie-Fischer 1999 NSE99-A2095) — major, out of scope.
- **#229 floor = the inherent ANGULAR half-angle-thread INTERPOLATION floor**
  (scales with quadrature, **independent of the τ-clamp**). Sphere needs **S32**
  for clean full-ladder O(h²); cylinder structurally blocked (2-D (η,φ) closure
  needed; τ_raw=0 + duplicate azimuthal η). The ORIGINAL #229 framing was right.
- **W1 (τ-clamp unclamp) ✅ DONE + COMMITTED `b2d8a6d` 2026-06-13 (qa PASS +
  elegance PASS-WITH-NITS resolved).** Doc-honesty fixes folded in: `cell_balance.py`
  stale "clamped for curvilinear" comment corrected (sphere unclamped / cyl clamped)
  + a `.. note::` caveat on the `:label: morel-montry-clamp` docstring block (the
  equation-of-record still shows the clip + cites "Eq. 5" — full restructure +
  duplicate-label collapse [also in `structured_geometry.rst:253` w/ DRIFTED
  citation] + `matrix.rst:741` DEFERRED to W5, captured in task #17). The clamp is mis-cited (Bailey says τ∈[0,1], recommends τ_raw=Eq.43
  exact-on-linear; Hébert uses τ=½) and 100% spurious on physical fields. User
  DECISION: remove it (sphere only) + S32 gate. IMPLEMENTED:
  `reduced_operator.spherical_streaming` `tau_mm[n]=max(0.5,min(1,τ_raw))`→`τ_raw`
  (inlined, `tau_raw` name dropped — no clamped/raw distinction left) + Bailey
  Eq.43 comment; CYLINDER keeps clamp (structural `τ_raw=0` ÷0 block) + comment.
  Mixed signature confirmed: cleans coarse rate (1.978→1.995) yet RAISES S16 fine
  floor (7.3e-4→1.2e-3, fortuitous-cancellation removed); NOT bit-id on iso
  (~1e-12 FP-tail). Regen `sphere_2g_3reg_dd_n40` (k 1.380766→1.381001); other 2
  sphere snaps FP-tail only (no regen). Gate files (test-architect): clamp-silent
  1a/1b/4 + 3 W1 @slow aniso gates; dropped the staged discriminator xfail.
  **⚠ NEW BUG FOUND + FIXED (same class as W0, rank-d carve `6ae3da8` fallout):**
  the 2 legacy strict-xfail tests `test_sn_{sph,cyl}_aniso_mms_converges_second_
  order` used STALE index `scalar_flux.values[:,0]` (was `(nx,ng)`, now `(ng,nx)`)
  → extracted cell-0 broadcast vs the nx-profile → garbage L2 ~14 DIVERGING. They
  were FALSE XFAILS (failed, so strict-xfail satisfied, but for a reason unrelated
  to the stated "#229 floor"). Fixed `[:,0]`→`[0,:]`; both now xfail for the REAL
  reason (sphere orders [1.995,1.999,1.407] err 1.4e-3; cyl ~1.9e-2 floor) — W3
  will rename/migrate. Broad sweep under `-O`: **225 passed / 31 xfailed / 0
  failed** (15:42; xfails = 20× #206 cyl-matvec + 4 aniso-floor + misc, all
  expected). ⚠ Mode-8 NON-issue resolved: pytest REWRITES asserts in collected
  test modules → bare asserts FIRE under `-O` (verified empirically); the `-O`
  warning only covers NON-test-module asserts.
- **W4 (#9 P1 scattering) = REAL new coverage, ready** (test-architect spec
  validated: L0 per-ordinate 3.9e-15 sphere / 6.2e-15 cyl; L1 directional
  eigenvalue P1-lowers-keff het sphere Δ=−1.40e-2 + leakage-monotone control).

REVISED actionable program: W0 ✅ → W1 ✅ (`b2d8a6d`) → W2 ✅ → W4 ✅ → W3 ✅
→ W5 docs (NEXT) → W6 verify+commit. Diagnostics:
`/Users/rodrigo/.claude/jobs/84fd66f8/tmp/negclosure/diag_01..31`.
Agent memos: `.claude/agent-memory/{numerics-investigator,test-architect,literature-researcher}/`.

## ⭐ PROGRESS 2026-06-13 (cont.) — W2/W3/W4 DONE, GREEN, qa-pending; NOT yet committed

- **W2 ✅ (#233 filed)** — pole-cell O(h) characterization gate
  `tests/sn/verification/mms/test_curvilinear_pole_cell_characterization.py` (4
  tests, test-architect). GUARANTEE: global L2 O(h²) under midpoint + Hébert-3.430
  shell-volume-average (dual structurally-indep ref). CHARACTERIZATION: pole L∞
  order **lower-bounded `>0.8`** (a future LD/nodal fix is NOT blocked — vv #5/#17)
  + interior O(h²) + cyl pole-vs-volavg O(h) masked-by-midpoint. `verifies
  ("dd-curvilinear-scalar")` on the guarantee only. **#233 filed** (WONTFIX-for-DD;
  fix = LD #6 / cell-updates #158 / nodal; cross-ref'd). **#233 stays OPEN** (tracks
  the future higher-order scheme; W5 does NOT close it).
- **W4 ✅ (#9)** — P1 path-II scattering NEW coverage. L0
  `tests/sn/verification/mms/test_curvilinear_aniso_scattering_p1.py` (per-ordinate
  `S₁−S₀` vs SH-indep hand-ref 4.7e-15 sphere / 5.6e-15 cyl + neg-control). L1 in
  `tests/sn/eigenvalue/test_keff_curvilinear.py::TestSphereP1DirectionalEigenvalue`
  (het-vac sphere P1-lowers-keff Δ=1.40e-2 + leakage-monotone R4>R25>0). ⚠ leakage
  rows need `bc=BC.vacuum` (helpers default reflective→zero-leakage→vacuous).
  `verifies("pn-scatter","flux-moments")` — first CURVILINEAR exercise of these.
  L2 deferred (subsumed). 
- **W3 ✅ (#229)** — all 5 xfails REMOVED, 6 labels GREEN (audit **exit 0**).
  Sphere `converges_second_order` → coarse `orders[:2]>1.9` + band [1e-8,5e-3]
  (err[-1]=1.40e-3); cyl → floor band [1e-3,5e-2] (no rate). Sphere phase-C test
  RETIRED (label on the W1 S32 gate; aggressive retirement). Cyl phase-C REPURPOSED
  → `test_cyl_aniso_floor_scales_with_quadrature`. ⭐ **CORRECTION: the cyl floor is
  AZIMUTHAL (`n_phi`), NOT `n_mu` (the memo mislabeled).** Measured nx=80: n_phi
  8→16→32 = 1.90e-2→7.37e-3→3.10e-3 (scales 2.58×/2.38×); n_mu flat. Physical:
  η=sinθ·cosφ, M-M threads azimuth φ per polar μ-level. prescribed_inflow sphere
  test (`test_mms_prescribed_inflow.py`): dropped xfail + rate gate, band [1e-8,5e-3]
  + kept value `assert_allclose(2e-2)`.
- ⚠ **3rd FALSE-XFAIL found** (after W1's 2 stale-index): the prescribed_inflow
  sphere "divergence" (negative L2 orders) is NOT a bug — solution CORRECT to 0.26%
  (max rel 2.6e-3 @nx160); L2 floors at ~2.4e-3 (#229) + the negative orders are the
  **#233 pole-cell O(h) creep** as the interior converges. CONFIRMED W1-innocent
  (diverges identically pre-W1 @ `6c712ee`). The xfail reason ("band fails") was
  wrong; now asserts value+band, no rate.
- Combined `-O` verify: **48 passed / 0 failed** (1:20). qa dispatched (W2+W3+W4).
- **NEXT after qa+commits: W5 archivist** (comprehensive `discrete_ordinates.rst`
  + error_catalog + close #229/#9; the DEFERRED W1 doc surfaces in task #17) → W6.

## Context

What began as "close #229 (aniso curvilinear MMS floor) + #9 (verify aniso
curvilinear scattering)" became a root-cause investigation (explorer +
literature-researcher + numerics-investigator). The investigation **falsified
the single-floor premise** and resolved the "#229 floor" into **three distinct
errors**, settled by a norm reconciliation (the test-architect measured
volume-weighted L2; the numerics-investigator measured pointwise/L∞):

| # | Error | Dominant norm | Quadrature-scaling? | Status |
|---|---|---|---|---|
| 1 | **Sphere pole-cell spatial closure** O(h) at r→0 | L∞ / pointwise central flux (suppressed in L2 by V∝r²; invisible in k_eff) | no (spatial) | **NEW open defect** |
| 2 | **Sphere angular τ-clamp floor** ~7e-4@S16 | volume-weighted L2 at fine mesh | yes | **clean fix (unclamp)** |
| 3 | **Cylinder angular floor** | both | yes | **structurally blocked** |

**The τ-clamp `max(0.5,min(1.0,τ_raw))` is an artificial, mis-cited patch
(triple-confirmed):** Bailey-Morel-Chang say τ∈[0,1] and recommend *exactly our
unclamped* `τ_raw` (their Eq. 42/43) as the unique weight exact-on-linear-in-μ;
Hébert uses pure diamond, no clamp. On every physical converged solve the clamp
is **100% spurious** (prevents no negativity); unclamped sphere SI is stable on
all stress configs (thick absorber, near-vacuum, c=0.999, S64) with positive
scalar flux. Proven: unclamped threads a linear-in-μ field to 2e-16; clamped
floors at 1e-3 (sphere) / 1e-1 (cylinder).

**Cylinder is structurally unfixable by unclamp/partition**: product/LS
quadratures put the most-inward azimuthal ordinate exactly on −sinθ
(`τ_raw=0` → division by zero) and carry duplicate azimuthal η that a 1-D
η-thread cannot represent — it needs a 2-D (η,φ) angular closure (out of scope).

**The new discovery (W2):** the sphere **central flux is only first-order**
(O(h^0.95) at the r→0 pole cell; ~60% local relative error at nx=160), proven by
the isotropic-MMS control (clamp-silent → same floor in L∞) + per-cell rate pin
(`diag_11`, `diag_14`). Invisible to the L2 MMS gates and to k_eff. Distinct from
#168 (outer face, CLOSED) and ERR-058 (seed, CLOSED). No tracking issue yet.

Decision (user): **full program this session**, compacting at natural phase
boundaries and recovering context to keep going. Fix the whole thing.

Evidence: literature memo (this session), numerics report
`/Users/rodrigo/.claude/plans/mellow-swinging-breeze-agent-a2a8e982fb0b89620.md`,
diagnostics in `/Users/rodrigo/.claude/jobs/84fd66f8/tmp/negclosure/diag_01..15`,
test-architect spec `.claude/agent-memory/test-architect/curvilinear_aniso_229_9_verification.md`.

## The one architectural subtlety (W1)

A *dynamic* negative-flux fixup (τ depends on ψ) would make the operator
**nonlinear**, breaking the **linear Krylov matvec twin** and the SI≡Krylov
twin identity (Pattern-2 discipline, Cardinal Rule 2). Since the fixup is never
needed on physical fields, the principled W1 is: **drop the clamp, use the
linear unclamped `τ_raw`** in single-source `reduced_operator` (both twins
inherit it) — **sphere only**. Cylinder keeps its clamp (the `τ_raw=0`
structural block) with an explanatory comment. No nonlinear fixup. The iso
solution is unchanged (clamp silent on flat-in-μ → iso gate bit-identical); the
aniso solution improves (floor removed).

## Plan (sequenced; ⟂ = parallelizable; ▣ = natural compaction boundary)

### W0 — pre-existing RED fix [prereq, trivial]
`tests/derivations/test_sn_mms_anisotropic_symbolic.py:259,313`:
`Q_numerical[:, 0, :, 0]` → `[:, 0, :]` (source is 3-D `(N,ng=1,nx)` since the
rank-d carve `6ae3da8`; pure test-index migration). Verify the 2 `@foundation`
symbolic guards green.

### W2-investigate — pole-cell fixability [⟂ dispatch EARLY, long pole]
numerics-investigator + literature-researcher (parallel): is the sphere
pole-cell O(h) fixable to O(h²) (the WDD pole-face IC / starting-direction
*spatial* closure at r→0), or inherent to curvilinear WDD at the coordinate
singularity? Deliverable: fixability verdict + (if fixable) closure form +
blast radius. **Branch:** fixable→W2-fix; inherent→document + L∞/pole gate.
File the issue **either way**.

### W1 — sphere τ-clamp → unclamped τ_raw [solver win]  ▣ after W1
1. **test-architect FIRST** (proactive trigger — closure change crossing
   iso/aniso + touching the ERR-058 seed): FP-invariance gate (iso-MMS
   **bit-identical**; aniso-MMS L2 → clean O(h²); SI≡Krylov holds) + identify
   which curvilinear snapshots legitimately shift (they improve).
2. `reduced_operator.spherical_streaming:644`: `max(0.5,min(1.0,τ_raw))`→`τ_raw`
   + comment (Bailey Eq. 42/43 = exact-on-linear M-M weight; clamp mis-cited/
   over-conservative). Cylinder `:760` keeps clamp + structural-τ_raw=0 comment
   + #229 xref.
3. Regenerate shifted sphere snapshots (justified: more accurate). Verify iso
   bit-identical, aniso O(h²), SI≡Krylov. qa + elegance-enforcer.

### W2-fix-or-document — per the W2-investigate verdict  ▣ after W2
Fixable: pole-face WDD in `operator.py` + sweep twin `spatial/sweep.py` (MEDIUM;
both twins; ~6 curvilinear snapshots regen) + qa + elegance. Else: L∞/pole-cell
convergence gate + thorough doc. Promote `diag_11`/`diag_14` → `tests/sn/`.

### W4 — #9 path-II P1 scattering coverage [⟂ independent]
Per the validated test-architect spec: **L0** per-ordinate P1 source through a
curvilinear sweep — `S₁.apply−S₀.apply` vs structurally-independent hand-ref,
**per-ordinate** (not weight-summed), sphere (3.9e-15) + cyl (6.2e-15), with
`_require(max|S₁−S₀|>1e-6)` negative control. **L1** aniso-vs-iso curvilinear
directional eigenvalue: P1 forward-peaked **lowers** k_eff (het sphere
Δ=−1.40e-2) + leakage-monotone control (`(P0−P1)|R=4 > (P0−P1)|R=25 > 0`);
2-group asymmetric (fuel `get_mixture("A","2g")` + mod `get_mixture("C","2g")`).
**L2 deferred** (subsumed by L0+L1; closeout note).

### W3 — #229 gate retune [after W1]
- **Sphere**: W1 → clean L2 O(h²) → assert plain O(h²) + band (no two-claim
  split). Remove xfail.
- **Cylinder**: structurally blocked → accept floor. NO spatial-rate gate; wide
  ERR-026 re-floor catcher (`1e-3<err<5e-2`, `@catches("ERR-026")`) + verified
  floor-scaling (`err(n_mu16)<err(n_mu8)/2.0`). Rename label
  `sn-mms-cylindrical-aniso-spatial-convergence`→`...-angular-floor` (honest).
  Migrate all 6 labels to GREEN tests.
- `test_mms_prescribed_inflow_sphere_activates_redistribution`: drop strict-xfail
  (sphere clean post-W1); keep converged-value `assert_allclose(2e-2)`.

### W5 — docs + closeouts [archivist]  ▣ before W6
Comprehensive Sphinx (`discrete_ordinates.rst`): the three errors + **norm
distinction** (L2 vs L∞) + τ-clamp vindication (Bailey Eq.42/43, mis-citation,
unclamp) + pole-cell defect + cylinder 2-D-closure structural limit + #9
P1-in-curvilinear (L0 trick + L1 direction). `error_catalog.md`: τ-clamp finding;
ERR-026 → pole-cell surviving manifestation; new pole-cell ERR. Close #229
(sphere fixed, cylinder accepted), #9 (W4); open pole-cell issue; update cluster
memory + Lathrop (2000) NSE 134(3) acquisition note.

### W6 — verification + commit
Branch `fix/curvilinear-aniso-pole-and-clamp` (split if cleaner). audit exit 0;
pyright ratchet (W1/W2 touch orpheus/ — verify no new errors, baseline only if
legit); full curvilinear suite + SI≡Krylov + eigenvalue green; qa +
elegance-enforcer on all solver code; Sphinx clean. `Closes #229`, `Closes #9`.

## Critical files
- `orpheus/geometry/reduced_operator.py:636-644` (sphere τ), `:740-761` (cyl τ)
- `orpheus/sn/spatial/pole_angular_closure.py:692-725` (recurrence kernel)
- `orpheus/sn/operator.py` + `orpheus/sn/spatial/sweep.py` (W2 pole-face twins)
- `tests/sn/verification/mms/test_curvilinear_aniso_convergence.py` (W3)
- `tests/sn/verification/analytical/test_mms_prescribed_inflow.py` (W3)
- `tests/derivations/test_sn_mms_anisotropic_symbolic.py` (W0)
- `orpheus/derivations/continuous/mms/sn.py` (MMS builders, read-only)
- new: `tests/sn/eigenvalue/` (W4 L1), `tests/sn/verification/mms/` (W4 L0, W2 gate)

## Verification
- `python -O -m pytest` curvilinear suite green; iso-MMS bit-identical (W1 control);
  aniso sphere O(h²); SI≡Krylov; #9 L0 per-ordinate ~1e-15 + L1 directional.
- `python -m tests._harness.audit` exit 0 (labels migrated, no orphans/phantoms).
- pyright ratchet exit 0. Sphinx builds clean.
- qa + elegance-enforcer sign-off on every solver-code change.

## Compaction protocol (user-requested)
Natural ▣ boundaries: after W1, after W2, before W5. At each, compact + recover
from this plan file + the cluster memory `project_curvilinear_sn_cluster.md` +
GitHub issue state. The plan file is the durable backbone; update its phase
checkboxes as we go (copy a durable version into `ORPHEUS/.claude/plans/`).
