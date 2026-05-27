# Plan: Harmonic Moment-Space Typing, Operator Adjoint Correctness, and the `numerics`/`transport` Re-architecture

Status: **Phase 1 + pre-Phase-3 stabilisation COMPLETE (2026-05-26). Order
chosen: Phase 3 first, Phase 2 after.** Ready to start **P3.1 (import-linter
test)** — first step in the REVISED sequencing (P3.1 → P3.5 → P3.0 → P3.2 →
P3.3 → P3.4 → P3.6).
Branch: `refactor/moment-space-and-layering`, off `refactor/sn-operator-algebra`.
Worktree: `.claude/worktrees/moment-space-and-layering/`. This document is
self-contained by intent — every decision below is stated with its rationale so
a session with no conversational history can execute it.

---

## ⚡ Pickup checklist (read first after `/clear`)

1. **Re-load the cardinal context** via the session-start protocol (vv-principles,
   subagent-handoff-protocol, coding-elegance, `.claude/lessons.md`,
   `docs/development.rst`).
   - **CRITICAL**: Run `nexus briefing --db docs/_build/html/_nexus/graph.db`
     via Bash to confirm the Nexus graph is loaded. If the DB is missing,
     rebuild Sphinx via `.venv/bin/python -m sphinx -q docs docs/_build/html`
     BEFORE dispatching any sub-agent that depends on graph queries. The MCP
     server only spawns at session boot, so the graph must exist at startup
     — mid-session rebuild does NOT reload the MCP tool fleet (use the
     `nexus` CLI via Bash as a workaround within the session).
   - **Standalone scripts gotcha**: anything under `derivations/diagnostics/`
     run via `python <path>` will silently load the venv's
     `pip install -e .` resolution (usually points at the MAIN repo, not the
     worktree). Either prepend `PYTHONPATH=.` or use `sys.path.insert(0, str(_REPO))`
     at the top of the script. Pytest is unaffected (uses `<rootdir>`).
2. **Check this plan's outcome sections** (§§0.3 + 0.4 below) for the actual
   commit sequence + lessons learned. §0.4 covers the pre-Phase-3
   stabilisation work that fixed two pre-existing failures and logged
   ERR-052.
3. **Read the four prep memos** in `.claude/agent-memory/`:
   - `qa/phase1_moment_space_review.md` — Phase 1 QA review (Os + Cs analysis)
   - `method-implementer/phase1_self_review.md` — Phase 1 implementer self-review
   - `test-architect/phase3_verification_plan.md` — **Phase 3 verification plan**
     (P3.1 skeleton verbatim, FORBIDDEN_EDGES, WHITELIST, per-step gates, risk
     register). This is the load-bearing artefact for Phase 3 execution.
   - `numerics-investigator/vacuum_bc_eigenvalue_divergence.md` — ERR-052
     investigation memo (the hypothesis cascade that led to the fix).
4. **Order is RESOLVED**: Phase 3 first, Phase 2 after (user choice
   2026-05-26). See §"Phase 2 vs Phase 3 ordering" below for the trade-offs
   captured at decision time. Phase 2 ships as a clean tail addition into a
   reorganised `numerics/spaces/` package.
5. **Verify the baseline still green** via PROGRESSIVE CHUNKS (NOT the full
   SN suite — issue #198 documents memory growth to 50 GB that makes the
   full suite a bad gate). Target Phase 1's actual touch surface first:
   ```
   .venv/bin/python -m pytest tests/numerics/ \
       tests/sn/test_fixed_source_g1.py \
       tests/sn/test_scattering_operator.py \
       tests/sn/test_boundary_conditions.py \
       tests/sn/test_angular_average_operator.py \
       tests/sn/l1_analytical/ \
       -q --no-header
   ```
   Expected: ~750 PASS, 5 xfailed (pre-existing — sphere-4eg-krylov per
   issue #200; cylinder-2eg / cylinder-4eg variants).
6. **Proceed to P3.1** per §P3.1 below. The verification plan at
   `.claude/agent-memory/test-architect/phase3_verification_plan.md` §1.7
   has the verbatim 76-LoC test skeleton ready to copy into
   `tests/test_layer_imports.py`. Test-architect + explorer cross-verified
   zero production-side violations; the 3 derivations-only edges are
   pre-whitelisted with `RETIRE_IN_PX` retirement triggers. P3.1 lands
   GREEN and gates every subsequent P3.x commit.

---

---

## 0. Provenance and anchors

This plan resolves three coupled problems surfaced while auditing the harmonic
projection layer:

1. **ERR-039** — `HarmonicMomentReconstruction` (R) and the operator the code labels
   the "Hilbert adjoint Π*" differ by a factor `(2ℓ+1)` that is carried as a free
   literal in one class and is absent in the other, glued together by two prose
   warning boxes (`numerics/projection.py:382-447`, `sn/scattering.py:525-657`). The
   factor is the spherical-harmonic Gram matrix; it has no typed home.
2. **A `.T`/`.H` mislabel** — `HarmonicMomentProjection.apply_transpose`
   (`numerics/projection.py:342-379`) returns the naked synthesis `S₀` but its
   docstring calls it the Hilbert adjoint. `S₀` is the adjoint only under an
   unstated Euclidean metric assumption; it is not the representation transpose the
   generic `_AdjointOperator` machinery expects. The mislabel is itself a cause of
   ERR-039.
3. **Package misfiling** — `HarmonicMomentField` lives in `sn/`
   (`sn/harmonic_moment_field.py`) though it is not SN-specific;
   `spherical_harmonics.py` has no conceptual home; `numerics/` mixes pure
   mathematics with transport-flavored objects; `numerics/iteration.py`'s
   `KEigenvalue` fuses a Problem with a Solver.

Grand Report v3 anchors used throughout: §2 (core ontology), §3.3 (ABC use),
§5.3 Space hierarchy, §5.4 Basis hierarchy, §5.5 Field hierarchy, §5.7 Operator
hierarchy, §6.1 Space dunders, §6.3 `.T` vs `.H`, §10 PN, §19 harmonic projection,
§21.1 Eigenproblem hierarchy, §23 Problem classes, §24 Solver hierarchy.

Cardinal rules this plan obeys: correctness first; architecture explicit; surgical
edits to living documents; every object carries named V&V tests; fresh-context
legibility. Per the issues-as-notebook rule, assign GitHub issue numbers to each
Phase/step at implementation time — this document uses `P{n}.{m}` labels only.

### 0.1 Pre-flight reconnaissance (commissioned 2026-05-26)

The plan's P1.0 demands a consumer inventory before any edit. That memo lives at
`.claude/agent-memory/explorer/moment_space_plan_consumer_inventory.md` (483 lines,
file:line-cited throughout). Every claim in §§ 1, 1.5, and the Phase sections below
that names a file:line site is sourced from that memo. The memo is the empirical
backing for this plan; this document is the executable specification.

### 0.3 Phase 1 outcome (2026-05-26)

Phase 1 shipped as 11 commits on `refactor/moment-space-and-layering`:

```
1ab6233  P0    docs(plan): moment-space + layering plan — P0 reconnaissance landed
0eb9cf3  P1.1  feat(numerics): SphericalHarmonicBasis — canonical home for SH evaluator + Gram
f27ba47  P1.2  feat(numerics): SphericalHarmonicSpace — typed home for the SH Gram matrix
6ba9102  P1.3  refactor(numerics): MomentProjection + ReconstructionOperator — single (2ℓ+1) source
36159f5  P1.4  refactor(numerics): .T = w_n · S_0 ; .H = g_C · S_0 — ERR-039 root-cause fix
28d9275  P1.5  test(numerics): SphericalHarmonicSpace V&V test suite — ERR-039 endpoint
be92b92  P1.6  docs(numerics): collapse ERR-039 warnings + retire assert_galerkin_idempotency
f54295c  P1.7  refactor(sn): retire dead _build_rhs_{cartesian,spherical,cylindrical}
c5be4b0  P1.8a refactor(numerics): API simplification (cached_property + collapse factories + foundation marker)
1f1e85e  P1.8b docs(vv): ERR-039 round-2 update + ERR-051 + vv-principles anti-pattern #11
ea02ab5  P1.8b docs(vv): vv-principles anti-pattern #11 (followup commit; SKILL.md split from prior commit)
```

**The ERR-039 endpoint is reached:**

- The four operators ERR-039 conflated are SEPARATELY TYPED:
  `S_0` (basis.synthesize), `Π^T = w_n · S_0` (MomentProjection.apply_transpose),
  `Π* = g_C · S_0` (MomentProjection.H), `R = (2ℓ+1) · S_0`
  (HarmonicMomentReconstruction).
- The `(2ℓ+1)` literal lives in EXACTLY ONE production place:
  `SphericalHarmonicBasis.addition_theorem_factor`
  (basis/spherical_harmonic_basis.py:152). Two grep hits in other files are
  array-stride arithmetic (different semantics — verified).
- The `g_C = diag(4π/(2ℓ+1))` metric lives in EXACTLY ONE production place:
  `SphericalHarmonicSpace.inner_product_weights` (computed via
  `SphericalHarmonicBasis.metric_per_ell`).
- The broken `GalerkinProjection.assert_galerkin_idempotency` (asserted
  `Π R = I` instead of `Π R = 4π · I`) was retired alongside its sole caller
  test → **ERR-051** logged.
- `_build_rhs_{cartesian, spherical, cylindrical}` retired (dead code; the
  cartesian variant carried the inline `(2ℓ+1)` duplicate).
- `_AdjointOperator.apply` generalised via `_broadcast_for_leading_axes` helper
  to handle leading-axis metric broadcasting (operator.py:197-237).
- 13 new V&V tests in `tests/numerics/test_spherical_harmonic_space.py`, all
  `@pytest.mark.catches("ERR-039")` on the 5 math-identity tests.

**ERR catalog actions (2026-05-26):**

- **ERR-039 entry updated** with "Round 2 fix" block — Wave 0 fix made
  `apply_transpose` return bare `S_0` which was ALSO mislabeled; Phase 1
  installed the W and g_C metrics, separating `.T` from `.H`.
- **NEW ERR-051** — `assert_galerkin_idempotency` shipped wrong invariant
  (Π R = I vs Π R = 4π · I under no-prefactor SH convention); hidden by a
  negative-only test that fed a non-orthogonal Y so the wrong invariant
  produced the expected failure.
- **NEW anti-pattern #11** in `vv-principles/SKILL.md` — contract-validation
  methods MUST have positive AND negative tests; negative-only is a discipline
  failure that hides wrong-invariant claims indefinitely.

**Phase 1 lessons learned (durable insights):**

1. **Single source of truth for convention literals.** The (2ℓ+1) factor was
   in two places (R.from_Y constructor + sn/solver.py:930 inline); both got
   absorbed into the typed-space chain. Pattern: convention literals belong
   on the canonical TYPE, not at the consumer (Pattern 7 producer-side).
2. **LS quadratures are NOT exact for arbitrary SH integration.** They are
   optimised for moment integration in the SN transport equation; using them
   to test ⟨Y_l, Y_l'⟩ products at L≥2 gives ~24% diagonal error. Use Lebedev
   or product Gauss-Legendre × Chebyshev for SH-Gram tests. Recorded in user
   memory as `feedback_no_tolerance_relaxation` ("don't relax — use or
   implement an exact rule").
3. **Negative-only tests for validation methods are a discipline failure.**
   ERR-051's bug surface was the TEST, not the method per se. The test
   pattern "feed a broken case, assert raises" validates the raising behaviour
   but NOT the invariant claim. Anti-pattern #11 codifies the fix.
4. **Dead code is invisible noise that hides bugs.** Three `_build_rhs_*`
   functions were dead in production (the `TestNoLegacyMachineryInCallPath`
   sentinel already enforced their absence from the Krylov path), yet they
   carried a hidden `(2ℓ+1)` duplicate. Retiring them was free correctness +
   16% smaller solver.py.
5. **Reconnaissance memos pay 5-10× in implementation efficiency.** The
   explorer's consumer inventory (CC.1-CC.8) shaped Phase 1 in load-bearing
   ways: GroupFlux dropped (didn't exist), IsotropicSource added to P3.3,
   `numerics/iteration.py` confirmed L1-clean (P3.4 greenfield not carve).
6. **Cross-checks must be structurally independent, not just procedurally
   different.** The unit-vector cross-check in
   `test_R_equals_2l_plus_1_times_S0` is bit-identical to the einsum because
   single-term input eliminates the FP-reduction-order ambiguity — that's
   structural independence done right (L11).

### 0.4 Pre-Phase-3 stabilisation outcome (2026-05-26)

Between Phase 1's commit `8e5a263` and the start of Phase 3, a baseline
green-gate run on the inherited branch (`refactor/sn-operator-algebra@62994ad`)
surfaced two pre-existing failures that Phase 1 had inherited but not caused.
Per `feedback_fix_bugs_immediately`, the user directed fix-on-this-branch
BEFORE Phase 3 begins. Three additional commits shipped:

```
e3522fe  fix(numerics): ERR-052 — power iteration flux renormalisation
                        (subcritical eigenvalue underflow)
4b5ba90  test(sn): refresh test_z_on_gauss_legendre_raises regex for
                   current production message
6f88402  docs(memory): Phase 3 prep memos + ERR-052 investigation
```

**ERR-052 (the load-bearing fix):** `orpheus.numerics.eigenvalue.power_iteration`
(legacy, consumed by `solve_sn`) and `orpheus.numerics.iteration.KEigenvalue.solve`
(canonical going forward, P3.4 destination) both omitted the textbook per-step
flux renormalisation. For supercritical operators the iterate grew geometrically
but converged in 3-4 outers before overflow — works by accident. For subcritical
operators (the failing test's 2cm vacuum slab) the iterate decayed geometrically,
underflowed to denormalised FP (~1e-43) within ~60 outers, and the keff ratio
`production/absorption` became a 0/0 noise fixed point that `converged()` then
falsely accepted.

Fix: per-step renormalisation to **unit production rate** `∫νΣ_f·φ·dV = 1`
(NOT L2 norm — the canonical reactor-physics convention, scale-invariant in
keff, physically meaningful, rescalable to absolute flux at target power via
one multiplication by `P_target/κ`). Implemented as:
- `SNSolver.compute_production_rate(flux) -> float` (composes existing
  `compute_group_production_rate.sum()`)
- `EigenvalueSolver` Protocol gains optional `compute_production_rate`;
  `power_iteration` uses `getattr` for backward compat
- `KEigenvalue` gains a `production_estimator` constructor parameter
  (defaulting to `np.sum(F.apply(ψ))` for operator-triple generic case)

**Test-config tuning:** `tests/sn/l1_analytical/test_kinf_homogeneous.py::_TIGHT_KW`
bumped `max_inner: 300 → 1000`. The new iteration trajectory perturbs the GMRES
initial guess each outer step, so the unpreconditioned Krylov subspace re-builds
from a non-warmed state; sphere-2eg-krylov empirically needs ~600 GMRES iter per
outer (vs ~50 pre-fix) to hit `inner_tol=1e-12`. All 28 (coord × ng × inner_solver)
L1 variants pass at the original tight tolerances (`rtol=1e-10` keff, `rtol=1e-9`
spectrum) post-bump. Production `solve_sn(max_inner=200)` default is untouched.
Posted as comment on issue #200 (block-inverse preconditioner) with full
empirical evidence — the ERR-052 fix exposed that preconditioning is structurally
load-bearing for tight L1 gates, not just a performance nice-to-have.

**Second pre-existing failure (test_z_on_gauss_legendre_raises):** the test
asserted `pytest.raises(ValueError, match="mu_z")` against
`AngularAverageOperator.from_quadrature(quad_GL, axis="z", ...)`. The production
code does raise, but with a different message: the GL adapter exposes
`mu_z = zeros(N)` (not absent), so the from_quadrature path falls through past
the early "requires a quadrature with mu_z" guard and into the downstream
"no outgoing ordinates on axis='z'" check. Both messages encode the same
defect; the test now pins the actual one.

**Phase 3 prep memos** (commissioned 2026-05-26, all under
`.claude/agent-memory/`):
- `test-architect/phase3_verification_plan.md` — P3.1 import-linter test
  specification (FORBIDDEN_EDGES + WHITELIST + 76-LoC skeleton),
  layer assignment proposal for all 14 packages, per-step verification
  gates (P3.2–P3.6), AngularFlux CC.4 design-note resolution, risk register,
  concurrence on the REVISED sequencing.
- `explorer/import_walk.py` + `explorer/phase3_import_edges.json` —
  ast-walker inventory of cross-package import edges (200 edges total).
  Independently verifies the test-architect's "zero production violations"
  claim; the only 3 forbidden edges are in `derivations/` (pre-whitelisted).
- `numerics-investigator/vacuum_bc_eigenvalue_divergence.md` — full ERR-052
  hypothesis cascade. The sub-agent was sandbox-blocked from Bash/Write but
  the static analysis narrowed the suspect set enough for the main agent
  to reproduce, discriminate, and fix.

**Durable insights from the stabilisation:**

1. **Power iteration MUST renormalise.** The textbook formula isn't a
   nice-to-have. Without it, subcritical eigenproblems underflow, supercritical
   may overflow, and only "iteration count small enough" hides the issue.
   Unit-production-rate is the principled convention (`ψ` is then a power-density
   basis directly rescalable to any reactor power).
2. **Bit-identity in tight tests can be a happy accident.** The
   sphere-2eg-krylov L1 test was getting `rtol=1e-10` keff because the initial
   guess `ones(...)` happened to be the eigenmode for homogeneous reflective
   sphere (uniform flux). Once the iteration was principled, the bit-identity
   coincidence broke; the right action was to give the inner solver the budget
   it actually needs (`max_inner=1000`), NOT to relax the test tolerance.
3. **Standalone-script sys.path resolution silently picks the venv-installed
   orpheus, NOT the worktree.** Took ~30 min of head-scratching during the
   ERR-052 fix to spot that `python derivations/diagnostics/diag_*.py` was
   running against the MAIN repo's `pip install -e .`, not the worktree's
   code. Always self-prepend the repo root to `sys.path` in standalone diags.
   Pytest is fine because `<rootdir>` auto-discovery via `pyproject.toml`
   uses the worktree.
4. **Nexus MCP must be verified at session start.** This session's prep
   work would have been substantially faster if the Nexus graph had been
   loaded — the explorer wrote a 200-LoC ast walker reinventing what
   `nexus migration --from orpheus.sn` gives for free. The MCP server only
   spawns at boot; the graph DB must exist at the configured path before
   the session starts. Mid-session rebuild populates the DB but does NOT
   reload the MCP tool fleet (the `nexus` CLI via Bash works as fallback).

### 0.2 Resolved open decisions (2026-05-26)

The four "open decisions" listed at the foot of the original draft are now closed
by user direction. Folded into the relevant sections below; preserved here as a
compact resolution table:

| # | Decision | Resolution | Rationale |
|---|---|---|---|
| 1 | L2 layer name | `transport/` | Matches existing vocabulary (`sn/`, `cp/`, `moc/`, `mc/` are transport-method packages); self-describing rule: "contains no solver and no discretization." |
| 2 | `transport/` flat vs subdivided | Flat until ~15 modules | Defer subdivision; current scope keeps it under that bar. |
| 3 | `kinetics/` re-examination scope | **In this plan as P3.6** | Avoids two coupled refactors landing weeks apart. |
| 4 | `assert_orthonormal` test name | **Retire `assert_galerkin_idempotency` entirely** (memo §A3 / CC.5: asserts wrong invariant; only called against a deliberately-broken pair to make a contract test fail). The genuine `Π R = 4π·I` identity is already pinned by `tests/numerics/test_projection_operators.py:194-231`. P1.6 deletes both the method and its sole caller. |

---

## 1. The ERR-039 root cause (reference for implementers)

Four operators are in play between the discrete-ordinate angular space and the
spherical-harmonic coefficient space. All four are `S₀` post-multiplied by a
diagonal; conflating them is the bug.

- `S₀` — naked synthesis: `c ↦ Σ_ℓm Y_ℓ^m(Ω_n) c_ℓ^m`. This is what
  `HarmonicMomentProjection.apply_transpose` currently returns
  (`projection.py:342-379`).
- `Πᵀ` — the representation transpose of `Π`. Since `Π_{ℓm,n} = w_n Y_ℓ^m(Ω_n)`,
  `Πᵀ = w_n · S₀`. The current `apply_transpose` is missing the `w_n`.
- `Π.H` — the Hilbert adjoint with the quadrature-weight metric `W` on the angular
  side and the SH Gram metric `g_C` on the coefficient side:
  `Π.H = g_V⁻¹ ∘ Πᵀ ∘ g_W = g_C · S₀`.
- `R` — the addition-theorem reconstruction: `R = (2ℓ+1) · S₀ = 4π · g_C⁻¹ · S₀`.

`g_C` is the spherical-harmonic Gram matrix: under the project's no-prefactor SH
convention, `⟨Y_ℓ^m, Y_ℓ'^m'⟩_{L²(S²)} = (4π/(2ℓ+1)) δ`, so `g_C = diag(4π/(2ℓ+1))`
per ℓ and `(2ℓ+1) = 4π · g_C⁻¹`. The `(2ℓ+1)` literal in
`HarmonicMomentReconstruction` (`projection.py:465`) is `g_C⁻¹` wearing a disguise.

Fix in one sentence: give the SH basis a `Basis` object whose `mass_matrix` is the
single definition of `g_C`; derive `R`, `Π.H`, and `Πᵀ` from it; replace the prose
warnings with executable identities.

---

## 1.5 Plan gaps closed by reconnaissance (memo CC.1–CC.8)

The explorer memo's cross-cutting findings change the shape of Phases 1 and 3 in
ways the original draft did not anticipate. The summary table below; each row is
folded into the Phase section it affects. **NEVER** treat these as optional —
they are load-bearing additions discovered by the reconnaissance and confirmed
against the live codebase.

| ID | Finding | Where folded |
|----|---------|--------------|
| CC.1 | `GroupFlux` does not exist as a class (`docs/theory/index_convention.rst:739` only). | P3.3 — drop from migration list. |
| CC.2 | `IsotropicSource` + `PerOrdinateSource` (`sn/sources.py:80-426`) are L2 vocabulary with the same `mesh.{ng,nx,ny}` profile as `HarmonicMomentField`. | P3.3 — add to migration list. |
| CC.3 | `_build_rhs_cartesian` carries an inline `(2*l+1)` literal at `sn/solver.py:862, 884-898, 930` — a duplicate of the addition-theorem factor. Pre-Wave-1 algebra, parallel to `R·Λ·M`. | **NEW P1.7** — retire this duplicate to honour P1.3's "exactly one place" claim. |
| CC.4 | `AngularFlux.from_flat_with_traces` (`sn/angular_flux.py:378, 420`) imports `build_equation_map_with_traces` + siblings from `sn/operator.py` — L2 leaking to L3 via the B1'' eq-map. | P3.3 — explicit design choice required; see "AngularFlux design note" in P3.3. |
| CC.5 | `assert_galerkin_idempotency` (`projection.py:153-191`) ships broken — sole caller (`tests/numerics/test_projection_operators.py:368-381`) feeds a non-orthogonal Y so the wrong invariant FAILS as a contract test. Never run against the real (Π, R) pair. | P1.6 — retire method + delete sole caller. Sibling 4π·I test stays (`tests/numerics/test_projection_operators.py:194-231`). |
| CC.6 | `KrylovAcceleration` (`numerics/iteration.py:459-697`) is a third sibling alongside `KEigenvalue` and `SourceIteration`. Heavily used by SN production (`sn/solver.py:618, 659, 1452, 1484`). | P3.4 — explicit inclusion as a third Solver in the family. |
| CC.7 | `numerics/iteration.py` is already L1-clean (no hard `orpheus.sn` imports; cross-package flow purely duck-typed via `_is_ravellable:156-169`). P3.4 has **no import-graph carve** — it is greenfield. | P3.4 — re-scoped as greenfield, plus retirement of `numerics/eigenvalue.py` legacy `power_iteration`. |
| CC.8 | `numerics/eigenvalue.py:168` ships `power_iteration` (legacy) still imported/exported by `numerics/__init__.py:3, 76`. Duplicates `KEigenvalue`. | P3.4 — explicit retirement target. |

---

## 2. Naming decisions (Grand-Report-grounded)

| Concept | Name | Grand Report | Replaces / current code |
|---|---|---|---|
| SH coefficient space | `SphericalHarmonicSpace` (`FunctionSpace`) | §5.3 | (no current object; ad-hoc) |
| SH basis | `SphericalHarmonicBasis` (`Basis`) | §5.4 | `numerics/spherical_harmonics.py:95` |
| The Gram / metric | `MomentMassMatrix` | §10 | the free `two_l_plus_one` literal (`projection.py:432, 465`) |
| Coefficient field | `HarmonicMomentField` | §5.5 | already correct name (currently mis-located in `sn/`) |
| Projection operator | `MomentProjection` | §5.7 | `HarmonicMomentProjection` |
| Reconstruction operator | `ReconstructionOperator` subtype | §5.7 | `HarmonicMomentReconstruction` |
| Operator codomain attr | `codomain` | §3.3, §6.3 | `range` (shadows builtin) |

`SphericalHarmonicMomentSpace` (the §10 PN `MethodSpace`) is **out of scope** — it
is PN-solver work, not part of this plan.

### 2.1 Anti-recommendations (per `lessons-L13`)

**NEVER** do the following in Phase 1; **instead** do the named alternative.

1. **NEVER** create a new `SHCoefficientField` parallel to `HarmonicMomentField`
   — **instead** extend the existing `HarmonicMomentField` to consume
   `SphericalHarmonicSpace.from_L(L)` as its `domain`; the field already exists
   (`sn/harmonic_moment_field.py:70-292`).
2. **NEVER** create a per-`L` variant of `SphericalHarmonicSpace`
   (`SphericalHarmonicSpaceL0`, `SphericalHarmonicSpaceL1`, …) — **instead** make
   `L` a field on the dataclass with `shape == (L+1, 2L+1)` so spaces of equal
   `L` compare equal via the inherited `(name, shape)` identity.
3. **NEVER** modify `HarmonicMomentReconstruction.apply` (the einsum at
   `projection.py:474-477`) — **instead** rewire only its CONSTRUCTION path: the
   new `from_spherical_harmonic_space(cls, space, Y)` constructor sources the
   `two_l_plus_one` array from `space.addition_theorem_factor`. The runtime path
   stays bit-identical; the convention literal moves into the space.
4. **NEVER** make `apply_transpose` return the bare `S₀` for back-compat
   ("two paths, one for old callers, one for new") — **instead** rewire `.T` to
   return the true representation transpose (`w_n · S₀`) and update the two test
   sites (`tests/numerics/test_projection_operators.py:268, 289`) that pinned the
   old (wrong) contract. Per memo §A1.a: zero production callers, test-only
   rewire.
5. **NEVER** introduce `SphericalHarmonicBasis` as a free module-level function
   set — **instead** make it a proper `Basis` class (`§5.4`) with `evaluate`,
   `mass_matrix`, `project`, `reconstruct` methods. The functional API
   `evaluate_real_sh` becomes a thin re-export shim during transition and deletes
   in P3.2.
6. **NEVER** retire `assert_galerkin_idempotency` without also deleting its sole
   test caller (`tests/numerics/test_projection_operators.py:368-381`,
   `TestAssertGalerkinIdempotencyMethod.test_method_signals_violation`) — **per**
   `lessons-L20` (retirement = test migration); the test exists only to
   exercise a method that ships broken on the real pair. Both delete together.

---

## Phase 1 — Harmonic moment space and adjoint correctness ✅ COMPLETE (2026-05-26)

See §0.3 above for the commit sequence (1ab6233..ea02ab5), Phase 1 outcome, and
lessons learned. The step-by-step narrative below remains as the executable
specification — useful for code archeology, V&V audit, and as the template for
Phase 2's space-algebra commits.

Goal: introduce `SphericalHarmonicBasis` and `SphericalHarmonicSpace`, make the SH
metric a single typed object, correct the `.T`/`.H` split, and replace the ERR-039
warning boxes with tests. New files are placed at their final architectural location
(see Phase 3) so Phase 1 adds no debt.

### P1.0 — Consumer inventory ✅ COMPLETE

Memo at `.claude/agent-memory/explorer/moment_space_plan_consumer_inventory.md`.
Folded into §§ 1, 1.5 above and the Phase 1 / 3 steps below.

### P1.1 — `SphericalHarmonicBasis` (the `Basis`, §5.4) ✅ commit 0eb9cf3

Create `orpheus/numerics/basis/__init__.py` and
`orpheus/numerics/basis/spherical_harmonic_basis.py`. New `numerics/basis/` package
(its first member; the Basis family will grow — §5.4 lists nine basis types).

`SphericalHarmonicBasis` implements the §5.4 `Basis` interface:
- `evaluate(directions) -> Array[N, L+1, 2L+1]` — relocate the body of
  `evaluate_real_sh` (`spherical_harmonics.py:95-170`) here as a method.
- `mass_matrix(measure) -> MomentMassMatrix` — the discrete Gram
  `Σ_n w_n Y_ℓ^m Y_ℓ'^m'`. For an exact quadrature this equals `diag(4π/(2ℓ+1))`.
- `project(field, measure)` — produces / is the `MomentProjection` operator.
- `reconstruct(coefficients)` — produces / is the reconstruction operator.

The SH normalization convention and `g_C` live here and nowhere else.

`spherical_harmonics.py` becomes a thin re-export shim of
`SphericalHarmonicBasis.evaluate` (preserves the `from orpheus.numerics.spherical_harmonics import evaluate_real_sh` API for the per-quadrature delegators —
`numerics/quadrature/directional.py:75, 393` and the four sn-quadrature wrappers).
The shim deletes in P3.2 when those delegators are rewired.

### P1.2 — `SphericalHarmonicSpace` and `MomentMassMatrix` ✅ commit f27ba47

Create `orpheus/numerics/spaces/__init__.py` and
`orpheus/numerics/spaces/spherical_harmonic_space.py`. The introduction of
`numerics/spaces/` is the seed for the P3.2 reorganization; `space.py` /
`trace_space.py` migrate to that package in P3.2.

`SphericalHarmonicSpace(FunctionSpace)`, frozen dataclass, same subclassing pattern
as `TraceSpace` (`trace_space.py:251-263`):

- field `L: int = 0` (needs a default — dataclass inheritance forbids a no-default
  field after `inner_product_weights`'s default; `InflowTraceSpace` sets the same
  precedent at `trace_space.py:301-306`).
- classmethod `from_L(cls, L)`: `name="spherical_harmonic_space"`,
  `shape=(L+1, 2L+1)`, `inner_product_weights` = the `MomentMassMatrix` diagonal
  broadcast to `(L+1, 2L+1)` — row ℓ holds `4π/(2ℓ+1)` in the `2ℓ+1` valid slots,
  `0` in the `|m|>ℓ` padding.
- `__post_init__` asserts `shape == (L+1, 2L+1)`.
- `__eq__`/`__hash__` inherit `(name, shape)` from `FunctionSpace`
  (`space.py:129-135`); `shape` already encodes `L`, so equal-`L` spaces compare
  equal and the metric is correctly excluded from identity.
- derived properties: `metric_per_ell -> (L+1,)` array `4π/(2ℓ+1)`;
  `addition_theorem_factor -> (L+1,)` array `2ℓ+1`.

The space's metric **is** `SphericalHarmonicBasis.mass_matrix`; the space holds it,
the basis defines it.

### P1.3 — `MomentProjection` / `ReconstructionOperator` rewire ✅ commit 6ba9102

In `orpheus/numerics/projection.py`:

- Rename `HarmonicMomentProjection -> MomentProjection` (§5.7). Add a `codomain`
  property returning `SphericalHarmonicSpace.from_L(self.L)`.
- Reparent `HarmonicMomentReconstruction` under a `ReconstructionOperator` base
  (§5.7). Replace its free `two_l_plus_one` array (`projection.py:430-432, 450-465`)
  with a value sourced from the space: add
  `from_spherical_harmonic_space(cls, space, Y)` setting
  `two_l_plus_one = space.addition_theorem_factor`. Keep `from_Y` as a thin
  back-compat shim so the single production call site (`sn/scattering.py:633`)
  does not break.
- The `(2ℓ+1)` literal must exist in exactly one place after this step —
  `SphericalHarmonicSpace.addition_theorem_factor`. (Modulo the test fixtures at
  `tests/numerics/test_projection_operators.py:173, 293` which legitimately
  re-derive the array as the test spec; those are the test as the spec.)
- The legacy duplicate at `sn/solver.py:862, 884-898, 930` is handled in P1.7.

### P1.4 — `.T` / `.H` correction in `operator.py` (§6.3) ✅ commit 36159f5

§6.3 mandates: `.T` is the representation transpose; `.H` is the Hilbert adjoint
including domain/codomain inner products.

- Establish `A.T` as the genuine representation transpose across the operator
  algebra; `A.H` (the existing `_AdjointOperator`) composes `g_V⁻¹ ∘ A.T ∘ g_W`.
- Fix `MomentProjection`: its `.T` must return the true representation transpose
  `Πᵀ = w_n · S₀`, not the bare `S₀`. With `codomain` carrying `MomentMassMatrix`
  and the angular domain carrying quadrature weights, `M.H` then computes correctly
  via the generic machinery as `g_C · S₀`.
- The bare `S₀` (naked synthesis) remains available as an explicit basis method
  (`SphericalHarmonicBasis.synthesize`), not mislabeled as an adjoint.
- This is a behaviour change to `apply_transpose`; per memo §A1.a, **zero
  production callers** — only the two test sites
  `tests/numerics/test_projection_operators.py:268, 289` need updating.

**Technical risk — layout broadcast.** For `A.H` to be correct, metric arrays must
broadcast against the *leading* axis (quadrature weights on the ordinate axis,
`MomentMassMatrix` on the `(ℓ,m)` axes) while the projection operators stay
layout-agnostic in trailing cell/group axes. `space.py:213-225` already solves
the fixed-layout case by pre-shaping weights. Phase 1 resolves the
projection-operator case concretely: the space carries the metric shaped for axis-0
broadcast and `_AdjointOperator` broadcasts trailing. The general case is Phase 2.

### P1.5 — V&V test ✅ commit 28d9275

Create `tests/numerics/test_spherical_harmonic_space.py`. The ERR-039 prose becomes
executable:

- `SphericalHarmonicSpace.from_L(L).inner_product_weights` equals `4π/(2ℓ+1)`
  per ℓ.
- `∫_{S²} Y_ℓ^m Y_ℓ'^m'` against a fine Lebedev rule ≈ `MomentMassMatrix` — pins
  the SH convention.
- `R.apply(c) == (2ℓ+1) * S₀(c)` for random `c` — the ERR-039 identity.
- `M.apply(R.apply(c)) == 4π * c` on band-limited `c` — the genuine
  Π·R = 4π·I identity. Cross-references the existing
  `tests/numerics/test_projection_operators.py:194-231` sibling test that
  already pins this for the (Π, R) pair on a Lebedev quadrature.
- `M.H` computed generically equals `g_C ⊙ S₀`.

All marked `@pytest.mark.catches("ERR-039")` + `@pytest.mark.l1` +
`@pytest.mark.verifies(<label>)` where `<label>` matches the new
`docs/theory/spherical_harmonics.rst` equations added in P1.6.

### P1.6 — Docstring surgery + retire `assert_galerkin_idempotency` ✅ commit be92b92

Three actions in one step:

1. Replace the two ERR-039 warning boxes (`projection.py:402-422`,
   `sn/scattering.py:555-560`) with one line each: `R` is the addition-theorem
   reconstruction `= 4π·g⁻¹·S₀`, `g` = `MomentMassMatrix`; identity pinned by
   `test_spherical_harmonic_space.py`; `R` is not the Hilbert adjoint. Keep the
   ERR-039 reference.
2. Correct `MomentProjection`'s former-`apply_transpose` docstring
   (`projection.py:342-378`): `.T` is the representation transpose;
   `.H` is the metric-aware adjoint.
3. **Retire `assert_galerkin_idempotency`** per CC.5 and Resolution #4:
   - Delete the method body (`projection.py:153-191`).
   - Delete the sole caller test
     `TestAssertGalerkinIdempotencyMethod.test_method_signals_violation`
     (`tests/numerics/test_projection_operators.py:368-381`). The genuine
     `Π R = 4π·I` invariant is already pinned by the sibling test at
     `:194-231`.
   - Delete the `GalerkinProjection` ABC method declaration of
     `assert_galerkin_idempotency` (the docstring section at
     `projection.py:148-151` referencing it).

**Honest caveat.** `R` and `M.H` share a `(SphericalHarmonicSpace -> angular)`
signature; `domain`/`codomain` typing alone cannot make a mixup a type error.
Phase 1 achieves: one definition of the metric, the relationship as a test, the
warnings collapsed to a signpost. They are different classes
(`ReconstructionOperator` vs `_AdjointOperator`); that is the available type-level
separation.

### P1.7 — Retire `sn/solver.py:930` inline (2ℓ+1) duplicate (NEW per CC.3) ✅ commit f54295c

(Scope expanded principle: all three `_build_rhs_{cartesian, spherical,
cylindrical}` retired together, not just the cartesian variant carrying the
literal. All were dead code per `TestNoLegacyMachineryInCallPath`; aggregate
retirement per `feedback_aggressive_retirement`.)

### P1.8 — Phase 1 polish: post-QA-review elegance + catalog + skill (NEW; landed 2026-05-26)

QA's independent review surfaced four blockers + two catalog updates + one
skill self-improvement. Shipped as three commits:

- **P1.8a** (commit c5be4b0): `@cached_property` on
  `MomentProjection.{domain, codomain, range}` (avoid per-matvec allocation
  in Krylov inner loop); collapsed `from_spherical_harmonic_space` factories
  on both projection and reconstruction classes into one canonical factory
  per class (`from_measure(measure, L, *, Y=None)` and `from_Y(Y)` —
  `coding-elegance` Pattern 5); tagged 3 API/type tests
  `@pytest.mark.foundation` (was @l0/@l1; software invariants don't carry
  math-identity labels); deleted vacuous
  `test_from_spherical_harmonic_space_roundtrip` test.
- **P1.8b** (commit 1f1e85e): ERR-039 entry update + new ERR-051 in the
  V&V catalog.
- **P1.8b followup** (commit ea02ab5): anti-pattern #11 in
  `vv-principles/SKILL.md` (contract-validation methods need positive +
  negative tests).

The hand-rolled Pℓ source construction in `_build_rhs_cartesian` at
`orpheus/sn/solver.py:862, 884-898, 930` carries the addition-theorem factor inline:

```python
qS += (2 * l + 1) * (sig_s[mid][l].T @ SUM) / sum_w   # line 930
```

This site predates the Wave-1 `R·Λ·M` algebra. After P1.3 lands the canonical
`MomentProjection` + `ReconstructionOperator`, this site is a twin path.

Action per `feedback_unify_after_two_instances` + `lessons-L21` (sweep and matvec
are different applications of the same operator):

- Migrate `_build_rhs_cartesian` to consume `ReconstructionOperator` and
  `MomentProjection` from the canonical pipeline (the same primitives
  `ScatteringOperator.build_aniso_source` already uses, `sn/scattering.py:627-650`).
- After the migration, verify no remaining production consumer of the legacy
  inline form. If `_build_rhs_*` is fully supersedable by
  `ScatteringOperator.build_aniso_source`, retire it entirely
  (consult dependency audit before deletion per `lessons-L20`).
- Tests: the SN regression suite at `tests/sn/regression/snapshots/` must continue
  passing bit-identical to within `(iteration count × ULP)` per
  `vv-principles` §"Bit-identity vs principled-equivalence" criterion 3 (the
  reordered reduction may FP-drift; the value must agree under
  structurally-independent reference).

The plan's "exactly one place" claim depends on this step landing.

---

## Phase 2 vs Phase 3 ordering (RESOLVED 2026-05-26: Phase 3 first)

User chose Phase 3 first. The trade-off rationale captured at decision time
is preserved below for the next session's reference. **Phase 3 is the active
work.**

Both phases WILL be done. Only the order is open. The choice has trade-offs:

**Phase 2 first (smaller, ~2-3 commits):**

- P2.1 introduces `DualSpace` + `Space.dual()` that generalises the inline
  `_broadcast_for_leading_axes` helper added in P1.4 (commit 36159f5) into a
  proper Space method. The helper currently has only ONE caller
  (`_AdjointOperator.apply`); the existing project rule
  (`feedback_unify_after_two_instances`) argues for deferring the
  abstraction until a second caller arrives. Phase 2 doing it now is a
  deliberate rule-departure justified by Grand Report v3 §5.3 + §6.1
  scheduling.
- P2.2 introduces `TensorProductSpace` (`Space.__mul__`) + `DirectSumSpace`
  (`Space.__add__`). The plan's original deferral note —
  "sequence them when those consumers are scheduled" — was written when
  issues #172 (`SumOfTensorProductsOperator` §15.2 streaming-form L1 test)
  and #173 (`TensorProductOperator` adjoint distributivity test) were
  unscheduled. They are still unscheduled. **If Phase 2 is done before
  Phase 3, P2.2 ships infrastructure with no immediate test consumer** —
  the unification rule violation is intentional but creates orphan code.
- Phase 2 commits are clean additions; no consumer rewires needed.

**Phase 3 first (larger, ~7+ commits):**

- P3.1 (import-linter) lands first and lights up violations as failures —
  every subsequent P3.x commit is gated by the linter going green.
- **P3.5 (range→codomain rename) should sequence EARLY** — eliminates the
  dual-name overhead I introduced in MomentProjection (Phase 1 carries
  `range` + `codomain` because the framework still reads `range`). Doing
  P3.5 first means subsequent moves (P3.2, P3.3) work against a single
  canonical attribute name. **Recommended order: P3.1 → P3.5 → P3.0 →
  P3.2 → P3.3 → P3.4 → P3.6.**
- Phase 3 fixes layering NOW (immediate value); Phase 2 adds infrastructure
  for future consumers (deferred value).
- **CC.7 confirmed**: P3.4 is greenfield (no SN import in
  numerics/iteration.py); much smaller than originally framed.

**Recommended pick:** **Phase 3 first**, then Phase 2 as a clean codetail
addition. Reason: Phase 3 closes more debt (dual-naming, dead-code retirement
audit, kinetics restructure) and unblocks Phase 2 by giving it a clean
`numerics/spaces/` package to land into. Phase 2 then ships without scope
contortions.

**User has final say** — the plan supports either order with no rewriting.

---

## Phase 2 — Space algebra

Goal: the general `Space` compositional structure (§5.3, §6.1). ERR-039 does not
require this; the generic adjoint and tensor/streaming work do. Phase 1 added
the leading-axis broadcast helper inline in `_AdjointOperator`; Phase 2 promotes
it to a typed Space method.

### P2.1 — `DualSpace`, `Space.dual()`

`DualSpace` as a concrete `Space` (§5.3). `Space.dual()` (§6.1) returns it; for an
inner-product-bearing space it is isomorphic to the space with a covariance tag and
the inverse metric. This generalizes the axis-0 broadcast fix from P1.4 to arbitrary
spaces and makes `A.H` derivable for any operator between metric-bearing spaces.

Tests (per `vv-principles` §"structural independence"):

- `dual(dual(S)) == S` for every space the project ships.
- `(A @ B).H == B.H @ A.H` for non-self-adjoint operator pairs.
- `assert_adjoint_consistency(A, V, W, atol)` — the ⟨A x, y⟩_W = ⟨x, A.H y⟩_V
  identity over random samples.

### P2.2 — `TensorProductSpace`, `DirectSumSpace`

`Space.__mul__ -> TensorProductSpace` (independent axes; §6.1, §5.3),
`Space.__add__ -> DirectSumSpace` (coupled fields, e.g. `FluxSpace + PrecursorSpace`;
§22.4, §6.1). These unblock the §15 streaming-as-tensor-product and the
delayed-precursor block structure; sequence them when those consumers are scheduled.

Tests:

- `(V * W).shape == V.shape + W.shape` (tuple concatenation).
- `(V + W).shape[axis]` equals the sum along the coupled axis.
- `TensorProductOperator` adjoint distributivity (the existing issue
  #173 anchor; this lands its prerequisite).
- `SumOfTensorProductsOperator` §15.2 streaming-form (issue #172 anchor).

---

## Phase 3 — The `numerics` / `transport` re-architecture

Goal: organize the package by a single criterion, enforced by a test. This Phase
subsumes the earlier "numerics-internal reorg" and the `range -> codomain` rename.

### P3.0 — The criterion and the layer contract (documentation)

Record this as the framing section of the Phase 3 issue and in `docs/`.

**Criterion.** A module's home is the lowest-knowledge layer whose vocabulary
suffices to define it. Imports flow only from more-knowledge to less-knowledge.

| Layer | Knows | Packages |
|---|---|---|
| L1 mathematics | functional analysis, linear algebra, measure theory; no neutrons | `numerics/` |
| (input) | geometry; nuclear data | `geometry/`, `data/` |
| L2 transport vocabulary | the transport equation's objects; method-agnostic | `transport/` |
| L3 discretization | one method's machinery | `sn/`, `pn/`, `moc/`, `cp/`, `mc/`, `diffusion/` |
| L4 orchestration | wiring a run | thin driver / entry point |

`Problem` and `Solver` are not a layer. They are math-object families that recur at
L1 (generic: `Eigenproblem`, `PowerIteration`, `TimeStepper`), L2 (transport-named:
`CriticalityProblem`, `AlphaEigenproblem`, `FixedSourceProblem`), and L3
(method-specific: `SweepPreconditionedSolver`, DSA) — exactly as `Field` and
`Operator` do.

Sphinx page: extend or create `docs/architecture/layering.rst` with this table +
the criterion + the import-linter test that enforces it.

### P3.1 — The import-linter test (written FIRST, before any move)

Add `tests/test_layer_imports.py` asserting the import contract: no module in
`numerics/` imports `transport/` or any method package; no module in `transport/`
imports a method package; no method package imports another. This is the safety
net for every subsequent move and converts the layout into a named, enforced
invariant (§2 philosophy applied to the package graph).

**STATUS (2026-05-26):** The test-architect verification plan at
`.claude/agent-memory/test-architect/phase3_verification_plan.md` §1.7
contains the **verbatim 76-LoC test skeleton** including `FORBIDDEN_EDGES`
(all 14 packages assigned) and `WHITELIST` (3 derivations-only entries with
`RETIRE_IN_PX` retirement-trigger comments). The explorer's ast walker
(`.claude/agent-memory/explorer/import_walk.py`) independently confirmed
**zero production-side L1/L2/L3 violations**, so P3.1 lands GREEN out of
the gate — the migration to-do list is empty for production code, and the
3 whitelist entries ARE the migration to-do list for `derivations/` (each
will retire as its containing module moves to method-side tests or to an
external benchmark harness).

This rewrites the original framing: "Expect it to fail initially — its
initial failures are the migration's to-do list" was correct conceptually
but empirically the project is already structurally clean (the discipline
was enforced by earlier waves; the linter just makes it executable). The
real Phase 3 work is the file moves (P3.2/P3.3) and the conceptual carve
(P3.4 Problem/Solver split, P3.6 kinetics dissolution), NOT import-graph
surgery.

Implementation: copy the verbatim skeleton from the test-architect memo
into `tests/test_layer_imports.py`, run `pytest tests/test_layer_imports.py
-q --no-header`, expect all parametrised modules to PASS. Commit as
`test(harness): P3.1 — import-linter test for L0/L1/L2/L3 contract`. Two
tolerances per the original plan: (a) `TYPE_CHECKING` imports of L3-only
types inside L1/L2 modules are allowed (string annotations don't create
runtime edges); (b) the 3 WHITELIST entries each carry `RETIRE_IN_PX`
retirement-trigger comments so future contributors can audit them.

### P3.2 — `numerics/` internal reorganization by math-object family

Mechanical; behaviour bit-identical. Restructure `numerics/` to:

```
numerics/
  spaces/         function_space.py (= old space.py),
                  trace_space.py,
                  spherical_harmonic_space.py     (from P1.2)
  basis/          spherical_harmonic_basis.py     (from P1.1)
  measures/       discrete_measure.py (= old measure.py),
                  quadrature/
  operators/      operator.py, projection.py
  solvers/        iteration.py split (per P3.4)
  registry.py     (stays flat — cross-cutting, §4)
```

Per Cardinal Rule 2: import-linter test from P3.1 stays green at every commit.

P3.2 deletes the back-compat `spherical_harmonics.py` shim (introduced in P1.1) and
rewires the four sn-quadrature delegators (`sn/quadrature.py` —
`GaussLegendre1D`, `LebedevSphere`, `LevelSymmetricSN`, `ProductQuadrature`) plus
the `numerics/quadrature/directional.py:75, 393` site to consume the new
`SphericalHarmonicBasis.evaluate` API.

### P3.3 — Introduce the `transport/` layer; migrate fields + sources

Create `orpheus/transport/`. Migrate the following L2 vocabulary out of `sn/`:

| Source | Destination | Notes |
|---|---|---|
| `sn/harmonic_moment_field.py:HarmonicMomentField` | `transport/fields/harmonic_moment_field.py` | Clean — uses only `mesh.{ng,nx,ny}` per memo §B1; `SpatialGroupMesh` Protocol replaces `SNMesh` import. |
| `sn/scalar_flux.py:ScalarFlux` | `transport/fields/scalar_flux.py` | Cleanest of the four (memo §B3) — only `SNMesh` TYPE_CHECKING import. |
| `sn/angular_flux.py:AngularFlux` | `transport/fields/angular_flux.py` | MIXED per memo §B2 — see "AngularFlux design note" below. |
| `sn/sources.py:IsotropicSource` | `transport/sources/isotropic.py` | Per CC.2 — analogous profile to `HarmonicMomentField`. |
| `sn/sources.py:PerOrdinateSource` | `transport/sources/per_ordinate.py` | Per CC.2 — `quad.weights.sum()` read is geometry-agnostic (`from_isotropic:332`). |
| ~~`GroupFlux`~~ | — | **DROPPED** per CC.1 — class does not exist. |

Keep `transport/` flat until it exceeds ~15 modules (Resolution #2). The
above is 5 files; sub-packaging `fields/` and `sources/` IS justified by the
existing semantic grouping but introduces no further hierarchy.

**Forced dependency cleanup (per CC.4).** `HarmonicMomentField` currently carries
`SNMesh` but reads only `mesh.ng/nx/ny` (verified memo §B1). In `transport/` it
must depend on a generic `SpatialGroupMesh` Protocol:

```python
class SpatialGroupMesh(Protocol):
    ng: int
    nx: int
    ny: int
```

`SNMesh` continues to satisfy this Protocol by structural typing — no consumer
change needed at the SN call site.

**AngularFlux design note (CC.4 resolution).** `AngularFlux.from_flat_with_traces`
at `sn/angular_flux.py:378, 420` imports three helpers (`build_equation_map_with_traces`,
`solution_to_angular_flux_with_traces`, `pack_with_traces`) from `sn/operator.py`.
These are SN B1''-eq-map primitives — L3-discretization machinery, not L2 vocabulary.

**Decision (folds into P3.3):** keep the pure `AngularFlux` algebra (storage,
dunders, `integrate_angular`) in `transport/fields/angular_flux.py`; the
B1''-aware adapters become an L3 extension shipped from
`sn/angular_flux_b1pp_adapter.py` (or, equivalently, a `classmethod` injected via
the existing `SNMesh.zeros_angular_flux` factory that wraps the L2 base with the
SN-specific adapters at construction time). This is option (a) in the memo's
§B2: "thin `sn/angular_flux_b1pp_adapter.py` (L3 layer) that augments the L2
`AngularFlux` with the `from_flat_with_traces` factory at the SN package boundary."

Existing tests `tests/sn/test_angular_flux_with_boundary.py` and the
`tests/numerics/test_iteration_angular_flux.py:35` (which imports `AngularFlux`
through a cross-package path) drive the boundary placement — the L2 base must
satisfy them without the eq-map adapter.

**Optional follow-up (not blocking):** Split a generic harmonic coefficient field
(`numerics/`) from the thin transport specialization (`transport/`) — the memo's
§B1 design note suggests `l_block`, `isotropic_part`, `truncate` are operations
on a graded SH coefficient space (L1 mathematics), while `scalar_flux()` and the
mesh binding are transport-specific. Treat as a follow-up, not a blocker.

### P3.4 — Problem / Solver split (greenfield per CC.7)

Per memo §CC.7: `numerics/iteration.py` is already L1-clean. P3.4 has **no
import-graph carve**; the dependency does not exist. The cost is purely the
conceptual carve, plus the `numerics/eigenvalue.py:168` retirement.

**Greenfield introductions:**

- `transport/problems/criticality.py:CriticalityProblem` — declarative
  `A_loss ψ = (1/k) F ψ` (§21.1, §23). Transport-named.
- `transport/problems/fixed_source.py:FixedSourceProblem` — declarative
  `A_loss ψ = q` (§23).
- `transport/problems/alpha_eigen.py:AlphaEigenproblem` — sibling of
  `CriticalityProblem` (k scales production `F`; α scales the time-mass `T`;
  both `Eigenproblem`s differing only in which operator the eigenvalue
  multiplies).
- `transport/problems/initial_value.py:InitialValueProblem` — for the time-domain
  transient (§22). The kinetics work in P3.6 consumes this.

**Solver migrations (numerics/ → numerics/solvers/):**

- `numerics/iteration.py:SourceIteration` → `numerics/solvers/source_iteration.py`.
  Generic, consumes `(L, S, F, q)` operator triple.
- `numerics/iteration.py:KrylovAcceleration` → `numerics/solvers/krylov.py`. **Per
  CC.6** — the third sibling alongside SI and Power. Heavily consumed by
  `sn/solver.py:618, 659, 1452, 1484`.
- `numerics/iteration.py:KEigenvalue` → `numerics/solvers/power_iteration.py` as
  `PowerIteration` (the algorithm; the problem is the new `CriticalityProblem`).
  Renamed for clarity.
- `numerics/eigenvalue.py:power_iteration` (legacy free function) → **RETIRED** per
  CC.8. Remove from `numerics/__init__.py:3, 76`. Rewire `sn/solver.py` outer
  eigenvalue path through the new `CriticalityProblem` + `PowerIteration`
  composition.

**Construction sites that become Problem-construction sites** (per memo §C2):

- `sn/solver.py:499, 539` (SI path) — builds `(LC = StreamingOperator +
  CollisionOperator, ScatteringOperator, ZeroOperator)` triple, becomes
  `FixedSourceProblem(L=LC, S=scattering, F=ZeroOperator).solve(SourceIteration)`.
- `sn/solver.py:618, 659, 1452, 1484` (Krylov path) — analogous Problem
  construction with `KrylovAcceleration` as the solver.

Tests: `tests/numerics/test_iteration.py` + cross-package
`tests/numerics/test_iteration_angular_flux.py:35, 249` continue passing.

### P3.5 — `range -> codomain` rename

Rename the operator codomain attribute `range -> codomain` across `operator.py` and
all consumers, matching §3.3/§6.3 and removing the builtin shadow. Pure mechanical
rename; do it as the last step of Phase 3 so it rebases cleanly over the moves.

Use Nexus `rename` (`mcp__nexus__rename`) for the safe graph-aware rename; fall
back to `grep -rn '\\brange\\b' orpheus/` + sed only if Nexus is silent on
attribute references.

### P3.6 — `kinetics/` restructure (Resolution #3 — in this plan)

Under the new ontology `kinetics/` dissolves:

- **PointKinetics** is a 0-D reduced-order *model* (a peer of `diffusion/`). Keeps
  its package: `kinetics/point_kinetics/` or, more aligned with the layering rule,
  promoted to `models/point_kinetics/` if a `models/` super-package exists; lacking
  that, retain at `kinetics/point_kinetics/` and document the peer relationship.
- **Space-time transient** is an `InitialValueProblem` (`transport/problems/`) +
  a `TimeStepper` (`numerics/solvers/time_stepping.py`), i.e. Problem + Solver
  instances, not a package. Migrate the time-stepping algorithms (BDF1, BDF2,
  Crank-Nicolson, etc.) into `numerics/solvers/time_stepping.py`; the transient
  driver becomes a thin wiring at L4.

Audit `orpheus/kinetics/` to enumerate all current symbols:

- Production callers of each: drive the retirement/migration order per
  `lessons-L20`.
- Tests: rewire to the new homes per `feedback_retirement_means_test_migration`.

Sequencing within P3.6: (i) extract time-stepping primitives into
`numerics/solvers/time_stepping.py`; (ii) introduce `InitialValueProblem` in
`transport/problems/`; (iii) thin the existing kinetics drivers to consume both;
(iv) retire any duplicated abstractions; (v) the test suite stays green at each
step.

### Sequencing within Phase 3 (REVISED post-Phase-1 per QA learnings)

**Original sequencing** (from the pre-Phase-1 draft): P3.1 → P3.0 → P3.2 →
P3.3 → P3.4 → P3.5 → P3.6.

**REVISED sequencing** (recommended):

1. **P3.1** — import-linter test (FIRST, always). Expect failures; they ARE
   the migration to-do list.
2. **P3.5** — `range → codomain` rename. **NEW POSITION**: do this BEFORE the
   reorganisation moves (P3.2/P3.3). Reason: Phase 1's `MomentProjection`
   carries `range` + `codomain` because the framework still reads `range`;
   subsequent file moves and Problem/Solver splits should land against a
   single canonical attribute name, not a transitional dual-name.
3. **P3.0** — documentation of the criterion + layer table.
4. **P3.2** — `numerics/` internal reorganisation (`spaces/`, `basis/`,
   `measures/`, `operators/`, `solvers/`).
5. **P3.3** — introduce `transport/` layer; migrate fields + sources. Per
   CC.4: the `AngularFlux.from_flat_with_traces` design decision lands here
   (thin SN adapter at L3, L2 base in `transport/`).
6. **P3.4** — Problem/Solver split. **Re-scoped per CC.7**: greenfield, no
   carve. `numerics/iteration.py` is already L1-clean; only the
   `numerics/eigenvalue.py:power_iteration` legacy retires alongside.
7. **P3.6** — `kinetics/` restructure. PointKinetics becomes a peer of
   diffusion/; space-time transient = InitialValueProblem + TimeStepper.

Each step: behaviour-identical, green against the linter and the full suite
before the next. No big-bang move.

### Phase 3 — additional notes from Phase 1 learnings

**ReconstructionOperator ABC rule-departure (QA A6 — record here):** Phase 1
shipped `ReconstructionOperator` ABC with ONE concrete subclass
(`HarmonicMomentReconstruction`), borderline violating
`feedback_unify_after_two_instances`. The deferral rule's load-bearing
exception is Grand Report v3 §5.7's explicit sibling-of-`ProjectionOperator`
structure: the ABC ships eagerly because the report anticipates ≥2
instances (addition-theorem reconstruction + future PN method-space
reconstruction). This is recorded so a future implementer doesn't
inadvertently revert the ABC under the deferral rule.

**`from_Y` shim retirement (QA A2 — defer to Phase 3 close-out):**
`HarmonicMomentReconstruction.from_Y` is the canonical factory post-P1.8a
(the only one). It does NOT retire; the original QA framing of "from_Y is
the back-compat shim that retires" was inverted in P1.8a — from_Y is the
natural API and `from_spherical_harmonic_space` retired instead. No further
action.

**Smell items the user may want to revisit during Phase 3:**

- **QA A4** (`_broadcast_for_leading_axes` module-level free function with
  one caller) — defer until a second caller arises; P2.1 promotes it to a
  `Space` method if it lands then.
- **QA A8** (constructor collapse on MomentProjection / HarmonicMomentReconstruction) —
  done in P1.8a.

---

## V&V summary

| Phase | Named assertions |
|---|---|
| P1 | SH Gram = `4π/(2ℓ+1)`; quadrature Gram ≈ `MomentMassMatrix`; `R = (2ℓ+1)·S₀`; `M·R = 4π·I`; `M.H = g_C·S₀`; legacy `_build_rhs_*` inline `(2ℓ+1)` retires to single canonical source (P1.7); `assert_galerkin_idempotency` + sole caller deleted (P1.6) |
| P2 | `dual(dual(S)) == S`; `(A@B).H == B.H@A.H`; `assert_adjoint_consistency` on tensor/sum spaces; issue #172 / #173 anchors satisfied |
| P3 | import-linter contract holds at every commit; full suite bit-identical across each mechanical move (FP drift bounded by `iteration_count × ULP` per `vv-principles`); `power_iteration` legacy retired; kinetics restructure leaves test suite green |

---

## Explicitly out of scope / deferred

- `SphericalHarmonicMomentSpace` the PN `MethodSpace` (§10) — PN-solver work.
- §7 `MethodSpace` / `SpatialMesh` restructuring beyond the `SpatialGroupMesh`
  Protocol that P3.3 minimally requires.
- The generic-vs-transport `HarmonicMomentField` split — noted in P3.3
  "Optional follow-up", not in this plan.
- The B1'' equation-map machinery's relocation — P3.3 keeps it at L3 via the
  thin SN adapter; a future Phase 4 could promote it if a second method needs
  the same flat-vector pattern.

---

## Anti-recommendations summary (per `lessons-L13`, §2.1 above)

Six explicit "do NOT" rules ship with this plan; an implementer who reaches for
any of them should STOP and re-read §2.1.

1. NEVER create a parallel `SHCoefficientField` — extend `HarmonicMomentField`.
2. NEVER per-`L` variant `SphericalHarmonicSpace` — make `L` a field.
3. NEVER modify `HarmonicMomentReconstruction.apply` — only its construction.
4. NEVER make `apply_transpose` return both old and new — rewire fully; tests
   update.
5. NEVER ship `SphericalHarmonicBasis` as free functions — it's a `Basis` class.
6. NEVER retire `assert_galerkin_idempotency` without deleting its sole caller.

---

## Pointers

- Reconnaissance memo: `.claude/agent-memory/explorer/moment_space_plan_consumer_inventory.md`
- Grand Report v3 (architecture target): `.claude/plans/neutron_transport_grand_report_v3.md`
- `coding-elegance` skill (the operational guide for Cardinal Rule 2):
  `.claude/skills/coding-elegance/SKILL.md`
- `vv-principles` skill (the V&V hierarchy + bit-identity vs principled equivalence):
  `.claude/skills/vv-principles/SKILL.md`
- `subagent-handoff-protocol` skill (orchestration rules; test-architect dispatched
  proactively per its "Proactive triggers" table):
  `.claude/skills/subagent-handoff-protocol/SKILL.md`
- ERR-039 entry: `.claude/skills/vv-principles/error_catalog.md`
- Lessons referenced: L13 (briefs name existing types), L17 (convention crosswalk),
  L18 (Pattern 7 at producer), L20 (retirement requires dependency audit),
  L21 (sweep + matvec share one strategy)
