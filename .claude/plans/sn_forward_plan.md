# SN forward plan — open work after D5b (the post-N-D-polymorphism backlog)

> **Durable in-repo recovery anchor** (project rule: plans live in `ORPHEUS/.claude/`).
> Consolidated **2026-06-17** from a plan triage that archived 19 completed/superseded
> plans (see the "Archived plans index" at the bottom). This is the SINGLE forward anchor
> for the SN spatial/angular discretization arc; per-effort detail lives in the archived
> plans under `.claude/plans/archive/`.
> Evergreen north-star kept active alongside this file:
> `neutron_transport_grand_report_v3.md` (the native-architecture design reference).

---

## Where we are — the banked foundation (all on `origin/main` @ `928d361`)

The SN solver's **spatial-closure polymorphism** and **N-D foundation** are DONE and
merged. Nothing below needs revisiting:

- **#206** — the matvec/sweep WDD recurrence unification + the 1-D matvec moved off the
  operator into `_OneDimScanWalk`. Merged, closed.
- **#240 Phase 2** — principled-over-bit-identical + the three-layer separation
  (agnostic `DiscretizationScheme` · SN loss-representation · `InvertibleOperator`),
  Steps A–D. The spatial closure is now a clean polymorphic family
  (`DiamondDifference` + `LinearDiscontinuous`, N-D, trait-dispatched, **zero** scheme
  `isinstance`). Merged.
- **The N-D layout + sweep-strategy campaign** (#220 / #222 / #225) — d=1/2/3 Cartesian
  admitted end-to-end (axis-native `from_axes`, no `Mesh3D`); the first-class
  `SweepStrategy` / `Compatibility` / `supports()` polymorphic selection. Merged.
- **#240 D5b** — the multi-D **UBLD** Linear-Discontinuous closure (S1→S4 + D6): the
  `n_spatial_moments` reduction lattice, the unified all-d moment matvec, the
  diffusion-limit closure (**ERR-061**, the sweep/global-frame involution), the
  strengthened `ld-cartesian-2d` stress MMS, and the full UBLD theory narrative. Merged.
- **The curvilinear-aniso program** (#229 / #195 / #168) — closed; the residual
  *structural* work is #233 / #235 (below).

The product `(spatial-scheme ⊗ angular-scheme)` is therefore **half-built**: the SPATIAL
factor and the Cartesian corner are realized; the remaining work concentrates on the
ANGULAR factor and the curvilinear corner.

---

## PRIMARY next effort — #236: realize the (spatial ⊗ angular) discretization product

> **STATUS (2026-06-17): IN PROGRESS — Phase 1 (the pairing-validity surface, the headline)
> DONE + committed on `feature/sn-spatial-angular-product` (NOT pushed; 4 commits
> `52966c3`→`5a01670`).** Detailed campaign plan (Phase 2 τ-relocation carve + Phase 3
> separability gate, with the τ producer/consumer map) =
> `.claude/plans/issue_236_spatial_angular_product.md`. The sub-task framing below is the
> ORIGINAL scoping; Phase 1 covered ST4 (1a) + ST2 (1b) + confirmed ST1 already satisfied.

**Readiness (assessed 2026-06-17, post-D5b).** The issue's "Current state" table predates
the #240 Step D / D5a / D5b campaign and is stale. At current HEAD:

- ✅ **The spatial factor is a real polymorphic family** (DD + LD both concrete, N-D) —
  the issue's *"DD only — Step/LD/EC are docstring examples + NotImplementedError"* is
  obsolete. Acceptance criterion "≥2 spatial occupants + documented path" is met for
  Cartesian.
- ✅ **The Cartesian seam is realized** (sub-task 2, Cartesian) — the 1-D scan and the
  2-D wavefront route the spatial closure through the scheme/coefficient model
  (`cartesian_scan_coefficients` / `source_emission`), not the hardcoded `0.5*(in+out)`
  the issue describes (the #239/D5a coefficient-model lift closed that).
- ✅ **The selection + pairing-validity *pattern* is proven** — D5b's
  `SweepStrategy`/`Compatibility`/`supports()` (with curvilinear-LD as a *declared
  structural exclusion*) is exactly the mechanism sub-task 4's pair-guard needs.
- ✅ **The typed-space separation is concrete** — `SpatialMomentSpace` (S3-A0) vs
  `SphericalHarmonicSpace` make the tensor framing real; the **spatial** diffusion-limit
  factor is solved (ERR-061), one half of sub-task 4's pair-property.

**What remains (the harder, angular half — D5b was entirely spatial):**

1. **Sub-task 1 — disambiguate LD-spatial vs LD-angular in the type system.** Now
   *substantive* (spatial-LD is a real symbol, not a docstring): `SpatialScheme.LINEAR_DISCONTINUOUS`
   ≠ `AngularScheme.LINEAR_DISCONTINUOUS`. Cheap opener. (Groundwork in task #21.)
2. **Sub-task 3 — extract the curvilinear angular redistribution into a selectable
   `AngularScheme`.** The Morel–Montry τ-thread currently enters the *spatial*
   cell-balance denominator; decouple the (quadrature-set, redistribution-closure) split.
   This is the architectural heart and is genuinely new.
3. **Sub-task 2, curvilinear** — route a non-DD spatial scheme through the curvilinear
   `_OneDimScanWalk` (today curvilinear-LD is a declared exclusion).
4. **Sub-task 4 — the pair-level `diffusion_limit_consistent` query** — the spatial
   condition (have: DD, LD/ERR-061) × the angular β-condition (Bailey-Morel-Chang),
   *both jointly required*. Encode as a pairing guard.
5. **Sub-task 5 — the space×angle separability characterization gate** (promote the
   numerics separability probes: Cartesian additive cross-term ≈ 0, curvilinear gating
   `E ≈ max(E_space, E_angle)`).

**The interference constraint (load-bearing for sequencing):** in curvilinear the spatial
floor (#233) and the azimuthal floor (#235) **gate** each other (`E ≈ max(E_space, E_angle)`),
so a curvilinear spatial fix and an angular fix pay off *together*. ⟹ #236's sub-task 3 and
#235 are naturally co-prioritized; Cartesian has no such coupling (harvest freely).

**Evidence memos (cold-start reading before planning the campaign):**
- `.claude/agent-memory/literature-researcher/space_angle_discretization_separability.md`
- `.claude/agent-memory/numerics-investigator/sn_space_angle_discretization_coupling.md`

**First move when starting #236:** dispatch the **explorer** for a precise
readiness/scoping pass (done-vs-remaining against the 5 sub-tasks + acceptance criteria at
current HEAD) — the issue snapshot is stale. Then plan the angular-factor campaign.

**Sibling architecture:** #205 (field storage × role typing — same `×`-product discipline);
#219 (MethodSpace ABC + builder registry — the product likely lives within that layer).

---

## Open follow-ups (each tracked by a GitHub issue)

| Issue | Title (short) | Origin / why deferred |
|---|---|---|
| **#247** | D5b-S5: consume moment-resolved external source + boundary trace | The slope-SOURCE half of the LM-1989 trap is unverified — the external Q̂ is zeroed and the scattering channel exercises-but-doesn't-constrain the slope-source sign (vv **Mode 10**). Needs a moment-resolved external-source entry + a slope-source-sign-sensitive gate + a moment-resolved boundary trace. |
| **#233** | curvilinear pole-cell spatial closure is O(h) at r→0 | Inherent to plain diamond; needs a higher-order curvilinear spatial scheme. Coupled to #235 (interference). |
| **#235** | cylinder needs a 2-D (η,φ) angular closure | The M-M 1-D η-thread cannot represent the cylinder's genuine (η,φ) variation; the azimuthal floor. The angular-closure enabler #236 sub-task 3 unblocks. |
| **#158 / #6** | remaining spatial occupants: `Step`, `ExponentialCharacteristic`; LD angular FE | `Step` is a near-stub (`scheme.py`); EC + angular-LD not started. The LD *spatial* increments (A/B/C/D) are DONE via D5b. |
| **#227** | d≥3 sweep-kernel perf | ScanMarch scan(x)∘march(y,z) generalization + MovingFrontierWindow widening (measurement-gated). |
| **#228** | Mode-8: the d=3 `_build_omega_dot_n` dead fail-loud guard | `getattr(q,"mu_z",None) is None` never fires (mu_z is always-present); discriminate all-zero-cosines + add a test. |
| **#215** | σ_r-sweep `A_wg.solve` is isotropic-limit-only | Ships 46–56% errors on anisotropic flux; the real scattering-rate win is consistent DSA (#2) / Krylov, not the σ_r fold. |
| **#246 / #245** | typed `SpatialMomentSpace` shape-probe predicate; relocate `AVERAGE_MOMENT`/`face_moment_tail` to `numerics` | D5b cleanups (S4-deadline #246 passed — fold into #236 or a cleanup pass). |

(Deeper architectural backlog also open: #2 DSA, #5 TSA, #11 GMRES, #200 block-inverse
preconditioner, #236-sibling #205/#219 — not in the immediate critical path.)

---

## Process hand-offs (the instruction-architecture flow — NOT a code commit)

These are intentionally **uncommitted on `main`** in the working tree and land via the
instruction-architecture flow (the `.claude/skills/*` are in the forbidden-to-commit set):

- **vv-principles Mode 10** ("activated-but-unconstrained") — added to
  `.claude/skills/vv-principles/SKILL.md` by the D5b-S4 qa review.
- **ERR-060 / ERR-061 / ERR-062** — appended to
  `.claude/skills/vv-principles/error_catalog.md` (the `@catches(...)` markers ARE
  committed in the tests; the catalog prose awaits harvest).
- **Phase-0 favoritism hand-off** — the #240 Phase-0 skill-text directive (bit-identity
  favoritism removal) was a HAND-OFF to the instruction-architecture flow, not landed on
  the (now-merged) feature branch.

---

## Deferred / low-priority

- **Capability-matrix generator (tooling)** — `capability_matrix_framework.md` steps 5–8
  were deferred "for later sessions" (the Nystrom/FN matrix generator). No active trigger;
  revisit if the FN method work resumes. (Plan archived.)

---

## Archived plans index (where the detail lives — `.claude/plans/archive/`)

Completed/superseded plans moved to `archive/` in this triage, with their residual:

| Archived plan | Status | Residual (→ issue) |
|---|---|---|
| `issue_240_phase2_step_d5b_ubld.md` | D5b complete + merged | #247 |
| `issue_240_d5b_s2_crosswalk.md` / `issue_240_d5b_s3_crosswalk.md` | D5b-S2/S3 design records, done | — |
| `issue_240_phase2_step_d_homing.md` / `issue_240_phase2_layer_separation.md` | Step D / layer separation done | — |
| `next_principled_polymorphism.md` | #240 Phase 0–2 done | EXIT → this plan (#236) + Phase-0 hand-off |
| `issue_206_carve.md` | #206 merged, closed | — |
| `issue_236_s2_seam_carve.md` | superseded by #206 | curvilinear seam → #236 sub-task 2 |
| `issue_158_spatial_cellupdate_carve.md` | LD increments done | Step/EC → #158 / #6 |
| `sn_space_angle_discretization_plan.md` | Tier-2 parent; LD spatial done | #236 / #233 / #235 → this plan |
| `curvilinear_aniso_pole_clamp_program.md` | program complete (#229/#195/#168 closed) | #233 / #235 |
| `cylinder_mr_variant_alpha_verification.md` | #168 closed | — |
| `sn_sweep_strategy.md` / `nd_layout_foundation.md` / `nd_foundation.md` | N-D campaign complete (#220/#222/#225) | #227 / #228 |
| `sn_development_sequence.md` / `phase2_o2a_crosswalk.md` | Wave-O → Phase 5 + #208 done | — |
| `trajectory_resolvent_hindsight_refactor.md` | one-time hindsight review (2026-05) | — |
| `capability_matrix_framework.md` | tooling, steps 1–4 done | steps 5–8 (deferred, above) |
