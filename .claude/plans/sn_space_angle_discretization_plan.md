# SN (space ⊗ angle) discretization — development plan (cold-start)

**Status (updated 2026-06-14, session 3):**
- ⭐ **§2 PIVOTED to #206.** Built the `CellUpdate` seam capability + cache delegation (the §2 steps 1–2,
  DONE + verified, uncommitted on `refactor/sn-cellupdate-seam-slab`). Then the dependency audit revealed
  the 1-D helpers (`_run_1d_sweep`/`_sweep_1d_unified` + the operator's `_compute_LpC`) are a **SHARED 1-D
  scan core** used by BOTH `CumprodScan` AND `ScanMarch` (1-D degeneration) — so "full unification" = a
  **shared-1-D-scan-frame extraction + moving the matvec off the operator** (the 1-D half of S6.3b,
  entangled with `_compute_decomposition`) = **#206**, NOT a fold-into-CumprodScan. User chose to pivot to
  the full #206 carve. **The §2 closure-routing becomes #206 Phase A.** Carve plan =
  `.claude/plans/issue_206_carve.md`; once #206 lands, RESUME here at Tier 2 (#158/#233 spatial accuracy).
  §2-seam detail = `.claude/plans/issue_236_s2_seam_carve.md` (its DONE block + the closure-routing reference).

**Status (updated 2026-06-13, session 2):**
- ✅ **Tier 0 (#234) DONE + merged ff-only to main + pushed (`eab05ab`).** All
  curvilinear-SN `[Bailey2009]_` cites → `[BaileyMorelChang2010]_`; orphan entry +
  temp warning deleted; M-M weight cite corrected to **Eq. 43** (β cite generalized —
  no confirmable B-M-C eq number); Sphinx clean (only the 1 pre-existing `paramref`).
- ✅ **Tier 1 MAP DONE (explorer)** + **#219 home decided (user): IN-PLACE on
  `SNMesh`/`_run_1d_sweep`; DEFER #219** (its MethodSpace layer is purely aspirational —
  nothing built; ⚠ unrelated `SNMethodSpace` BC-realizer struct shares the name → future
  rename hazard). **§1 reframed (user):** drop the enum idea — the two ABCs
  (`CellUpdateBase` vs `PoleAngularClosureBase`) already disambiguate LD-spatial(#158) ≠
  LD-angular(#6). Real gaps to evaluate: **(1) the LD occupant doesn't exist on either
  axis; (2) genuine independent space⊗angle instantiation.**
- ⭐ **EMPIRICAL SEAM FINDING (session 2, probe `$CLAUDE_JOB_DIR/tmp/probe_cellupdate_seam.py`
  — PROMOTION CANDIDATE → permanent "production honors injected cell_update" gate).**
  Injected a call-counting `DiamondDifference` subclass; ran `transport_sweep` on
  slab/sphere/cylinder/2-D-Cartesian with a non-zero source. **NO production default
  honors the spatial seam:** `CumprodScan` (ALL 1-D) and `ScanMarch` (2-D) inline DD →
  ZERO `cell_update` calls. Positive control: forced `MovingFrontierWindow` → `cell_kernel_batch`
  fired 56× (seam reachable + works; spy valid). So the seam is honored ONLY by the
  explicit-select `MovingFrontierWindow`/`FullFieldWavefront` (window opt + oracle), NOT
  the production default. **⟹ "independent SPATIAL instantiation" is wired (`SNMesh(cell_update=…)`)
  but INERT in every production sweep** (vs the ANGULAR axis, which IS live in production —
  the curvilinear sweep reads+uses `sn_mesh.pole_angular_closure`). This refines #236 §2:
  the gap is the production strategies, and the seam machinery already exists in the
  window/oracle as the in-tree template. **Recommended §2 approach:** "construct general,
  select narrow, specialize on measured cost" (user's own principle / the existing
  `sn_sweep_strategy.md` SweepStrategy polymorphism) — DD keeps its fast inlined scan;
  a non-DD spatial scheme routes through the seam-honoring per-cell `cell_update` path
  (the curvilinear `_run_1d_sweep` already calls `cell_update.update` for the degenerate
  cyl ordinate at `loss_representation.py:2127` — generalize that to the non-DD fast path).
  HIGH-VALUE §2 TARGET = the 1-D/curvilinear `CumprodScan` path (unblocks #233/#158-in-curvilinear);
  2-D `ScanMarch` bypass is SECONDARY (the window already gives a seam-honoring 2-D path).

Created 2026-06-13 at the close of the curvilinear-aniso #229/#9 program. This is
the FORWARD plan; the COMPLETED predecessor program is
`.claude/plans/curvilinear_aniso_pole_clamp_program.md` (W0–W6, merged `db1b779`).

## Recovery context (where things stand)

- `main` @ `db1b779`, synced with origin, working tree clean (only the 3
  never-stage untracked files: `.claude/plans/r1_phase_a_dim_agnostic_ultraplan.md`,
  `derivations/diagnostics/diag_s69_scanmarch_vs_window_bench.py`,
  `scratch/literature/`). NEVER `git add` those three.
- Predecessor program DONE + merged: curvilinear aniso τ-clamp unclamp (W1),
  #9 P1-scattering coverage (W4), #229 retune (W3), pole-cell characterization
  + ERR-059 (W2/W5). Issues **#229/#9/#230 CLOSED**.
- This plan was produced by a paired literature + numerics evaluation of the
  user's claim that spatial and angular discretization "are the same but don't
  have to be." **Verdict: TRUE-WITH-CAVEATS** — the correct structure is a
  tensor product `(spatial-scheme ⊗ angular-scheme)`; ORPHEUS already has the
  two registries but the product is not fully realized.
- Agent memos (cold-start reading, committed in `db1b779`):
  - `.claude/agent-memory/literature-researcher/space_angle_discretization_separability.md`
  - `.claude/agent-memory/numerics-investigator/sn_space_angle_discretization_coupling.md`

## The goal

Realize the `(spatial-scheme ⊗ angular-scheme)` discretization product so spatial
and angular schemes are **independently selectable** and **validly combinable**,
and use it to lift the two curvilinear accuracy floors (the #233 pole-cell
spatial O(h) and the #229/#235 angular floor). Tracking umbrella: **#236**.

## The findings that shape the plan (load-bearing)

1. **Cartesian is already a separable product.** The Cartesian kernel
   (`DiamondDifference.cell_kernel_batch`, `orpheus/sn/spatial/diamond.py`) takes
   no τ / no angular term; empirically additive (cross-term ≈ 0, O(h²)
   bit-identical across N=4/8/16).
2. **Curvilinear couples space+angle through SHARED CODE.** The Morel–Montry
   angular weight τ enters the *spatial* cell-balance denominator
   (`orpheus/sn/spatial/cell_balance.py:262-343`:
   `denom = 2|μ|·A_down + (ΔA/w)·(α_out/τ) + Σ_t·V`). The W1 *angular* τ-clamp
   and the #233 *spatial* pole-cell genuinely share one recurrence.
   Coupling sites: `orpheus/sn/loss_representation.py` `_run_1d_sweep`
   `:2142-2150` (angle→space) and `:2167-2176` (space→angle).
3. **The coupling is an implementation binding, NOT a math necessity.** Hébert
   2009 Eq. 3.431/3.437 reuse one diamond weight τ=½ for both closures; Bailey-
   Morel-Chang 2010 + Lathrop 2000 deliberately UNBIND the angular weight. W1
   (`b2d8a6d`) already exposed the angular knob.
4. **Two registries already exist (the product is partly built):**
   - Spatial: `CellUpdateBase(RegistryMixin)` (`orpheus/sn/spatial/cell_update.py`),
     injected `SNMesh(cell_update=...)` (`geometry.py:271`), default
     `DiamondDifference()` — **ONE occupant (DD)**; Step/LD/EC are docstring
     stubs + `NotImplementedError`.
   - Angular: `PoleAngularClosureBase(RegistryMixin)`
     (`orpheus/sn/spatial/pole_angular_closure.py`), injected
     `SNMesh(pole_angular_closure=...)` (`geometry.py:460`) — Identity +
     Morel-Montry + swappable `psi_half_seed`. **Genuinely selectable.**
5. **THE SEAM BLOCKER:** the curvilinear fast path INLINES DD's closed-form
   scan + hardcodes the DD closure `0.5*(in+out)` (`loss_representation.py:2004,
   2167`); `cell_update.update` is only called for the degenerate cyl-axis
   ordinate (`:2113`). ⟹ a non-DD spatial scheme in curved geometry is BLOCKED
   until the curvilinear sweep dispatches through the `CellUpdate` Protocol.
   Cartesian already honors the seam.
6. **THE INTERFERENCE (gating):** in curvilinear the #233 spatial O(h) floor and
   the #229 azimuthal angular floor do NOT add — `E(h,N) ≈ max(E_space(h),
   E_angle(N))` (measured cross-term 0.37; sphere spatial error plateaus at S8,
   only pays off at S32). ⟹ the curvilinear accuracy payoff of a spatial fix and
   an angular fix is COUPLED — fixing one axis without the other caps the gain.
   Cartesian has no such coupling. (Mirrors the literature's "both diffusion-
   limit conditions jointly required.") Separability probes:
   `$CLAUDE_JOB_DIR/tmp/diag_sep_*` (promotion candidates → the Tier-3 gate).
7. **The diffusion-limit consistency FACTORIZES** into a spatial condition
   (Larsen-Morel-Miller 1987) + an angular condition (Bailey-Morel-Chang 2010,
   β=0 Eq. 41-42), each per-axis-provable but jointly required.
8. **LD is double-named across axes** (literature does this too): spatial-LD =
   a `CellUpdate` subclass (#158); angular-LD = a `PoleAngularClosure` subclass
   (#6). The registry MUST disambiguate: `SpatialScheme.LINEAR_DISCONTINUOUS` ≠
   `AngularScheme.LINEAR_DISCONTINUOUS` (never one `LD` enum). #233 (spatial
   pole-cell) is fixed by a SPATIAL scheme (#158 / nodal), NOT #6 (angular).

## The ranked plan (dependency + importance)

```
Tier 0  #234  bib-cite migration ........... independent, trivial (warm-up)
Tier 1  #236 §1+§2  LD-disambig + curvilinear CellUpdate seam ... FOUNDATIONAL
          │  (unblocks ALL curvilinear spatial-scheme work)
   ┌──────┴───────────────────────────────────┐
Tier 2a #158 → #233   spatial accuracy        Tier 2b #235 (+#1/#3/#6) angular
   spatial-LD / step-characteristic              cylinder 2-D (η,φ) closure
   resolves the pole-cell O(h)→O(h²)             lifts the azimuthal floor
   (Cartesian standalone; sphere gated)          (depends #236 §3 sub-axis)
   └──────┬───────────────────────────────────┘
          │  (CURVILINEAR needs BOTH — the interference gating)
Tier 3  #236 §4+§5  diffusion_limit_consistent pair-guard + separability gate
```

**Critical path:** `#234` → `#236 §1+§2 (seam)` → { `#158/#233` ∥ `#235/#1/#3/#6` } → `#236 §4+§5`.

### Tier 0 — #234 (warm-up, no deps)
Migrate the ~6 `[Bailey2009]_` body citations in `docs/theory/` to
`[BaileyMorelChang2010]_` (the corrected entry + a temporary `.. warning::`
already exist from W5). Remove the warning; rebuild Sphinx clean. ~30 min.

### Tier 1 — #236 §1 + §2 (foundational; the enabler)
- §1: split `SpatialScheme` / `AngularScheme` so LD (and any shared name) is
  axis-typed. Cheap, foundational.
- §2: make `_run_1d_sweep` dispatch the spatial closure through the `CellUpdate`
  Protocol instead of inlining DD's scan/closure (`loss_representation.py:2004,
  2113,2167`); the Pattern-2 matvec twin must follow. **This is the blocker.**
  Verify bit-identical for the DD occupant (no behavior change), then a non-DD
  spatial scheme becomes constructible for sphere/cylinder.
- ⚠ Confirm whether this product should live INSIDE #219's MethodSpace
  foundational layer BEFORE the seam refactor (build it in the right home).

### Tier 2a — #158 → #233 (spatial accuracy; depends Tier 1 for curvilinear)
Add a higher-order spatial cell update (spatial-LD or step-characteristic) as a
`CellUpdate` occupant → lifts the pole-cell O(h)→O(h²). **Reference oracle now
LOCAL:** Larsen-Morel-Miller 1987 (which spatial schemes limit to a valid
diffusion discretization). Standalone Cartesian value; curvilinear-sphere gated
by the angular floor (plan with 2b). The #233 characterization gate
(`test_curvilinear_pole_cell_characterization.py`) is the target — its pole
order is lower-bounded only, so it goes GREEN-at-higher-order when fixed.

### Tier 2b — #235 (+#1/#3/#6) (angular accuracy; depends Tier 1 §3)
Cylinder 2-D (η,φ) angular closure → lifts the #229 azimuthal floor. Supported
by #1 (Gauss azimuthal quadrature) / #3 (φ-edges); #6 (angular-LD) is an
alternative higher-order angular closure. **Reference oracle now LOCAL:**
Lathrop 2000 (benchmarks diamond / weighted-diamond / linear-continuous /
linear-discontinuous / quadratic-continuous angular schemes vs analytic
two-region sphere — the candidate menu + verification oracle, ODE-solved
WITHOUT spatial differencing, so it isolates the angular axis).

### Tier 3 — #236 §4 + §5 (verify + guard)
- §4: `diffusion_limit_consistent` as a property of the (spatial, angular) PAIR
  (encode LMM-spatial + B-M-C-angular conditions; flag bad pairings). The
  diffusion-limit pair is now LOCAL (LMM 1987 Part I + II, B-M-C 2010).
- §5: promote `$CLAUDE_JOB_DIR/tmp/diag_sep_*` to a permanent separability gate
  (Cartesian additive / curvilinear gating).

## Issue map
- NEW: **#236** (architecture umbrella — DETAILED, the cold-start spec),
  **#233** (pole-cell/ERR-059, spatial), **#234** (bib), **#235** (cyl 2-D angular).
- EXISTING (commented w/ axis-disambiguation): **#6** (angular-LD), **#158**
  (spatial cell updates).
- SUPPORT: **#1** (cyl Gauss azimuthal), **#3** (φ-edges).
- SIBLINGS: **#205** (field storage×role typing — same ×-discipline), **#219**
  (MethodSpace foundational layer — likely the home for #236).
- ORTHOGONAL (separate tracks, not this thrust): #2 (DSA), #200 (block precond),
  #8 (Case's method reference).

## References (ALL NOW LOCAL in `scratch/literature/`, added 2026-06-13)
- **Larsen, Morel & Miller (1987)** "Asymptotic solutions of numerical transport
  problems in optically thick, diffusive regimes" JCP 69(2):283-324, DOI
  10.1016/0021-9991(87)90170-7 — the SPATIAL diffusion-limit condition (Tier 2a/3).
- **Larsen & Morel (1989) Part II** JCP 83(2):212-236, DOI
  10.1016/0021-9991(89)90229-5 — companion.
- **Lathrop (2000)** "A comparison of angular difference schemes for one-
  dimensional spherical geometry S_N equations" NSE 134(3):239-264, DOI
  10.13182/NSE00-A2114 — the ANGULAR scheme menu + verification oracle (Tier 2b).
- **Bailey-Morel-Chang (2010)** NSE 165(2):149-169 (already local) — the angular
  diffusion-limit β-condition.
- **Hébert (2009)** Ch.3 §3.9 (already local) — the shared-weight coupling.

## Acceptance criteria (the product is "done")
- Two axis-disambiguated scheme registries; ≥2 spatial occupants constructible
  (DD + one of LD/SC), ≥2 angular (DD-thread + M-M, already there).
- Curvilinear sweep dispatches the spatial closure through `CellUpdate` (DD
  bit-identical); a non-DD spatial scheme runs on sphere/cylinder.
- `diffusion_limit_consistent` pair-query exists + tested on a good and a bad pair.
- Separability characterization gate green + permanent.
- The #233 pole-cell gate improves (pole order > 1.8) once a higher-order spatial
  scheme lands; the cylinder floor lifts once #235 lands.
- Sphinx documents the product + the interference constraint (extend the
  `sn-curvilinear-aniso-norm-reconciliation` section).
