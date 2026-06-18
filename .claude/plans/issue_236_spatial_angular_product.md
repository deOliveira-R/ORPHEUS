# #236 — Realize the (spatial ⊗ angular) discretization product

> **Durable in-repo recovery anchor** (project rule: plans live in `ORPHEUS/.claude/`).
> Supersedes the session-scoped `~/.claude/plans/mellow-swinging-breeze.md` (the original
> plan-mode file). Phase 1 (the headline) is DONE + committed; Phase 2 (the τ-relocation
> carve, the architectural heart) + Phase 3 (separability gate) remain.
> Branch `feature/sn-spatial-angular-product` off `main` @ `3ac96b4`, **NOT pushed**.

## Context

The SN discretization is a tensor product of two independently-selectable axes — a SPATIAL
closure (`DiscretizationScheme`: Diamond Difference + Linear-Discontinuous, N-D) and an
ANGULAR factor (quadrature × `PoleAngularClosure`: Morel–Montry + Identity). Some properties
of the discretization are properties of the *pair* (chiefly the thick-diffusion limit, which
factorizes into a spatial LMM-1987 condition AND an angular BMC-2010 condition). This campaign
makes the angular factor genuinely first-class and the pairing-validity explicit/queryable/honest.

**Scope boundary (out — tracked research tail):** a *curvilinear*-LD spatial closure is
mathematically UNPUBLISHED → acceptance criterion 2 is blocked on RESEARCH (#158 curvilinear
arm / #6), NOT engineering. This campaign builds the framework + the honest exclusion.

**Recovery reading:** the explorer scoping memo
`.claude/agent-memory/explorer/issue_236_readiness_scoping.md` (done-vs-remaining at HEAD,
file:line) + the equation-pinned lit memo
`.claude/agent-memory/literature-researcher/diffusion_limit_consistency_per_scheme_booleans.md`.

---

## Phase 1 — pairing-validity surface (the headline) ✅ DONE + committed (bit-identical)

- **1a (ST4)** `52966c3` (feat) + `12c08ce` (chore): `diffusion_limit_consistent: bool` on the
  `DiscretizationScheme` Protocol/base (DD=True LMM-1987 Eq.4.24; LD=True Larsen–Morel 1989 II
  Eq.4.16 — needs the slope source, landed D5b-S3; Step=False-when-built LMM Eq.5.20) +
  `beta_first_order_consistent: bool` on `PoleAngularClosureBase` (M-M=True BMC-2010 Eq.42;
  Identity=True vacuously, Cartesian α≡0) + the joint predicate
  `pair_diffusion_limit_consistent(scheme, closure)` in NEW `orpheus/sn/spatial/pairing.py`
  (pure AND; the Cartesian collapse is encoded by `IdentityAngularClosure.beta=True`, NOT a
  branch). Corrected the stale `LinearDiscontinuous` header `.. warning::`. @foundation gate
  (positive + negative truth-table). Conformance mocks + genuine-bool teeth gained the 6th trait.
- **1b (ST2)** `c9f1ca5` (feat) + `5a01670` (chore): `supports_curvilinear: bool` (DD=True;
  LD=False — curvilinear LD unpublished) + a single-source `_curvilinear_capability(mesh)` gate
  in `loss_representation.py` composed into `CumprodScan.supports` + `ScanMarch.supports` (1-D
  arm) + a top-of-`default_for` check. Curvilinear-LD now rejected at SELECTION (frontend
  `Compatibility` / fail-fast `IncompatibleRepresentation`), not passed-then-raised-mid-sweep.
  Defense-in-depth preserved (kernel `_require_slab` + CollisionCache `affine_scan_coefficients`
  raises). `TestHonestCurvilinearSchemeSelection` + the `_fake`/conformance mocks/teeth gained
  the 7th trait.

**Both: elegance-enforcer PASS + qa SUPPORTED (gates mutation-proven). ALL GREEN** — 76
pairing+spatial, 45 dispatch+conformance, 3 curvilinear-LD guard, 460 sweep/core, 7 slab-LD MMS,
165 curvilinear sweep. **Pre-existing reds** in `tests/sn/operators` (#232 ymin 2-D BC +
curvilinear/vacuum bit-identity snapshot fragility) — NOT introduced here (the changes are
numerically inert: grep-confirmed no production reader of the new traits). **ST1 already
satisfied** (two separate registries; no single LD enum). **#158 commented** (forward-note:
register Step with `diffusion_limit_consistent=False` + LMM Eq.5.20 when built).

---

## NEXT — Phase 2 (ST3): τ is an angular property, not a streaming-geometry field

**The architectural heart.** The Bailey–Morel–Chang angular weight τ is currently computed in
the streaming-GEOMETRY factory and baked into the SPATIAL cell-balance. Relocate τ-PRODUCTION
onto the redistribution closure so the angular factor owns its defining weight; consumers read
an angular-supplied τ. The (quadrature × closure) pair IS the "AngularScheme" — realized as the
existing two `SNMesh` injection points with the closure owning its weight (**NO bundling-wrapper
class** — Pattern 5/elegance: build the primitive, not the product).

**⚠ This is a MULTI-SITE carve (not a single-site move) crossing the spatial↔angular boundary
(L17).** The proactive **`test-architect` dispatch is the FIRST action of Phase 2** — it owns the
convention crosswalk + the 4-leg verification spec BEFORE any code. The τ producer/consumer map
(verified 2026-06-17, confirm line numbers at pickup):

| Role | Site | Note |
|---|---|---|
| **PRODUCER** | `geometry/reduced_operator.py:688` | `tau_mm = (mu − mu_edge)/dmu` — this IS BMC Eq.43 (the M-M weight); sphere |
| **PRODUCER** | `geometry/reduced_operator.py:796` | `tau_mm_per_level` — cylinder |
| **PRODUCER** | `geometry/reduced_operator.py:495` | slab synthetic neutral `tau_mm = 1.0` (Identity closure) |
| consumer | `sn/spatial/cell_balance.py:306` | `tau = st.tau_mm` → spatial denom (`c_out=α_out/τ`, `c_in`) |
| consumer | `sn/spatial/diamond.py:230, :305` | DD's angular-closure branch (update + residual) |
| consumer | `sn/spatial/sweep_cache.py:298` | `tau[global_n] = st0.tau_mm` — the CollisionCache precompute |
| consumer | `sn/spatial/pole_angular_closure.py:660, :669` | `MorelMontryAngularSweep._tau_per_level = reduced.tau_mm` — the closure RECEIVES τ today |

**The move:** the M-M closure PRODUCES τ (from the α's / μ-grid it already has, BMC Eq.42/43);
Identity produces the neutral τ=1; the streaming factory stops pre-baking `tau_mm`; every
consumer reads the angular-supplied τ. **The test-architect crosswalk is authoritative for the
exact producer→consumer convention at each boundary** (do NOT over-design here).

**Verification (BIT-IDENTICAL — same τ math, new owner):** the closure-produced τ must equal
`(mu − mu_edge)/dmu` exactly. Gate the L14 4-leg standoff (sweep ≡ matvec ≡ structurally-
independent reference, under mesh refinement, on sphere + cylinder); the `DriftWarning` strict
gate on the curvilinear paths; the per-ordinate residual (L27 — NEVER weight-summed for an
angular-redistribution operator; a scalar residual is blind to a wrong angular closure).

**Sequence:** proactive `test-architect` (crosswalk + spec) → `method-implementer` (the carve)
→ `elegance-enforcer` + `qa` → resolve nits → commit (feat + chore) → archivist documents the
τ-ownership relocation + the product narrative on `theory/discrete_ordinates` (Cardinal Rule 3).

## THEN — Phase 3 (ST5): separability characterization gate

Promote the four `diag_sep_*` probes (space_angle / cyl / slab / slab_iso; construction
documented in `.claude/agent-memory/numerics-investigator/sn_space_angle_discretization_coupling.md`)
from the ephemeral job-dir into a permanent gate under `tests/sn/verification/` (analogous to the
#233 floor test). Asserts Cartesian additive cross-term ≈ 0 (separable); curvilinear
`E ≈ max(E_space, E_angle)` (the interference/gating constraint). Independent — lands any time.

---

## Verification (route-arounds + canonical invocation)

- Canonical `python -O -m pytest …` (host → `.venv/bin/python`). Bit-id/principled gate:
  `-W "error::tests.sn.regression._regression_assert.DriftWarning"`.
- Dirs per sub-step: `tests/sn/operators tests/sn/spatial tests/sn/sweep/core
  tests/sn/sweep/cartesian_2d tests/sn/solve` + the relevant MMS + (Phase 2) `tests/sn/sweep/curvilinear`
  (SLOW — ~13 min). **NEVER all `tests/sn`** (#212 `continuous_get` hang).
- Pre-existing reds to route around (NOT regressions): `tests/sn/operators` #232 ymin + the
  curvilinear/vacuum bit-identity snapshot fragility (8e15-ULP-on-near-zero signature).

## Execution discipline

- Per sub-step: implement → `elegance-enforcer` + `qa` → resolve nits → commit (one
  `feat`/`refactor` + a `chore(claude)` records commit). **Proactive `test-architect` FIRST at
  Phase 2 start.** Commit per sub-step is authorized; **pushing requires an explicit ask** (branch
  is 4 ahead of `3ac96b4`, not pushed).
- EXPLICIT paths only — NEVER `.claude/skills/*`, `CLAUDE.md`/`.claude/rules`/`.claude/hooks`,
  `docs/_build/`, the parallel-session `.claude/agent-memory/elegance-enforcer/MEMORY.md` +
  `nexus_runtime_overlay.md`, or the 3 forbidden untracked (`r1_phase_a_dim_agnostic_ultraplan.md`,
  `diag_s69_scanmarch_vs_window_bench.py`, `scratch/literature/`). `.claude/agent-memory/*`
  (explorer/qa/lit-researcher/etc., excluding the contaminated elegance one) writes ARE legitimate
  (chore). Forbid `git checkout/restore/stash` on tracked paths in every sub-agent brief (L28).
  Trailer `Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>`.

## Recovery STATE

Branch 4 ahead of `3ac96b4` (`52966c3`→`12c08ce`→`c9f1ca5`→`5a01670`), NOT pushed. Phase 1 done;
Phase 2/3 remain. Memory anchor = [[project-issue-236-spatial-angular-product]]. The FIRST action
resuming Phase 2 = dispatch `test-architect` with the τ producer/consumer map above.
