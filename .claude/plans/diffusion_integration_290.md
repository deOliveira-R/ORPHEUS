# Diffusion integration campaign — issue #290

## Context

`orpheus/diffusion/` is a 443-line MATLAB-port island (raw `(2, n_cells)` ndarrays, scipy
BiCGSTAB, own `TwoGroupXS`/`CoreGeometry`, string BCs, hardcoded ng=2 down-scatter flip
trick) predating the operator-algebra framework. Campaign #290 integrates it at SN's level
of integration or better: it is the genuine second consumer of scalar flux + scalar-flux
composites (nails down the typed-field architecture), it builds the in-algebra `A_diff`
that DSA for SN (#2, HIGH) consumes, and the diffusion pyright ignore lifts as a by-product.

Verified ground truth (explorer, at `d2a2a0c` = main): every #290 claim CONFIRMED, zero
drift. Design frames (cross-domain-attacker): the island's `_matvec` is literally
`(L+C−S)φ` hand-inlined and its fission source is the rank-1 `χ⊗νΣf/k` — the algebra
exists there in fused form. Verification plan (test-architect): full rewire map produced;
legacy `TwoGroupXS.transport` encodes Σ_tr directly, so the bit-identical Mixture bridge
keeps every analytic L1 anchor unmoved.

**Execution mode: surgical.** Main agent writes directly, user steers step-by-step with
AskUserQuestion checkpoints at phase boundaries. NO `method-implementer`. Dispatches
allowed/encouraged: `test-architect`, `qa`, `elegance-enforcer`, `archivist`, `explorer`.

## Context compaction protocol (harness `/compact` points)

Marked `⏸ COMPACTION Cn` below. Any phase boundary is a valid compaction point; the
marked ones are the recommended defaults. Before EACH compaction the main agent MUST:
1. **Commit everything** — no uncommitted state whose rationale lives only in
   conversation (agent-memory/plan edits included).
2. **Checkpoint the durable plan** (`ORPHEUS/.claude/plans/diffusion_integration_290.md`):
   append a `## CHECKPOINT Cn` block — phase statuses (done @ commit-hash / next),
   deviations from this plan + their rationale, open questions, mid-flight user rulings.
3. Then tell the user it is safe to `/compact`.
After compaction, the main agent re-anchors from the durable plan + `git log` (trust
git, not the summary — process-discipline rule), re-reads only the files the next phase
touches, and continues. Sub-agent memos live in `.claude/agent-memory/*` if re-grounding
is needed — do not re-dispatch for already-answered questions.

## User rulings (2026-07-03, this session — binding)

1. **Full composite mirror** — diffusion operators act on `FullField(bulk=ScalarFlux,
   boundary=scalar trace)`. The boundary operator is named **B** (project letter
   convention: L streaming-representation [= leakage here], C collision, S scattering,
   F fission, B boundary). Family: loss `A = L + C − S − B`, gain `F/k`.
2. **Trace carrier = partial currents (J⁺, J⁻)** — the ℓ=0 half-range moments under the
   SAME `|Ω·n|·w` metric as SN's TraceSpace. BCs are the albedo family `J⁻ = 𝒜·J⁺`.
3. **Vacuum = zero incoming current (J⁻ = 0), i.e. the Marshak realization (𝒜 = 0).**
   The current zero-flux-Dirichlet-called-"vacuum" is stale/unfaithful naming. Zero-flux
   Dirichlet is acceptable as a TEST condition but must be its own honestly-named law
   (𝒜 = −1 in the (J⁺,J⁻) basis). The analytic sine references survive mathematically
   intact — re-attributed to the zero-flux law, never renamed-away.
4. **Data seam: bit-identical encoding NOW** (`SigT = legacy transport, SigS1 = 0`;
   references unmoved), **physical P1 re-baseline AFTER** diffusion is essentially ready
   (file as follow-up issue at close-out) — prove the implementation correct first.

## Design (converged from the three pre-carve dispatches)

- **Operator family** (all typed on the scalar composite): `L` = leakage leaf `−∇·D∇`
  (FD assembly on `MaterialMesh`, harmonic-mean face D ≡ arithmetic-mean face Σ_tr —
  proven bit-identical); `C` = existing `MultiplicationOperator(σ_t)`; `S` = scalar
  composite wrapper over the SHARED K_iso kernels (`IsotropicScattering` +
  `IsotropicN2N` — never re-implement `Σ_s0ᵀφ`); `F` = the shared rank-1 `χ⊗νΣf`;
  `B` = `DiffusionBoundaryOperator`, the A_ss albedo-family block `J⁻ = 𝒜J⁺`.
  Removal Σ_r is a THEOREM (C−S in-group cancellation), not an input.
- **Resolvent**: explicit `MatrixInverseOperator(L + C − S − B)` — the taxonomy-blessed
  dense direct inverse (its docstring's worked example is literally `(−S)+(L+C)`;
  homogeneous-solver precedent `K = MatrixInverseOperator(loss) @ F`). Do NOT route
  through the structure-keyed `.inverse()` (Green splitting diverges for fine-mesh
  elliptic). Engines: `direct_eigenvalue` (dense `K = A⁻¹F`) + shared `power_iteration`
  via the 5-method `EigenvalueSolver` protocol — cross-gated at 1e-10.
- **Discretization**: cell-centered FD as-is (it IS RT0-with-mass-lumping; the mixed-form
  lift is a documented seam, trigger = off-face material interfaces needing O(h²)).
  Interior face currents stay condensed; only boundary (J⁺,J⁻) are trace DOFs.
- **ng-generic by construction** (K_iso kills the `sig_s[::-1]` flip trick) — unblocks
  #33/#34 structurally.

## Phases

Branch: `feature/diffusion-integration` off `main` @ `d2a2a0c`. First actions: copy this
plan to `ORPHEUS/.claude/plans/diffusion_integration_290.md` (durable home), create
local campaign tasks, and save the user-feedback memory (long-campaign plans carry
explicit compaction points + checkpoint protocol). Conventional commits per phase;
ff-merge at end.

### P1 — Data seam (`orpheus/data/macro_xs/mixture.py`)
- `Mixture.transport_xs` property: `Σ_tr = SigT − rowsum(SigS[1])`, absent-P1-order →
  `Σ_tr = Σ_t` EXACTLY (correct isotropic limit, not a fallback). `Mixture.
  diffusion_coefficient = 1/(3Σ_tr)`.
- Foundation tests: arithmetic pin (xs_library region with `sig_s1`), legacy bridge
  (`D == [1.52835091, 0.42462845]` bit-identical under the encoding), P0-only limit.

### P2 — Scalar trace substrate (numerics spaces + transport fields)
- Mint the scalar trace space on `MaterialMesh` (FunctionSpace sibling; do NOT generalize
  the angular `TraceSpace` — its quadrature coupling is SN's type-safety). Face layout:
  boundary faces × ng, carrying (J⁺, J⁻).
- Mint the `PartialCurrent` trace leaf (flux-role; source-sink/residual role siblings only
  as operator codomains require — mirror SN's role-triple shape, verify at contact).
- Widen `FullFieldSpace.from_blocks` (`full_field_space.py:196-215`), `FullField`
  boundary/mesh annotations (`full_field.py:126,238-254`) to admit the scalar trace.
- Fix ALL `scalar_flux.py` staleness (module docstring contradiction :32-62, class
  docstring `mesh: SNMesh` :122-125, fixed-shape claim — the anticipated second consumer
  has arrived).
- Intrinsic-property tests (foundation): `J = J⁺ − J⁻`; `φ_b = 2(J⁺ + J⁻)`; positivity
  `J± ≥ 0` for PHYSICAL laws (𝒜 ∈ [0,1]); albedo-family identities (vacuum J⁻=0 ⟺
  Marshak; reflective J⁺=J⁻; zero-flux 𝒜=−1 ⟹ φ_b=0); field algebra + mesh-identity guard.
- **L17 crosswalk FIRST**: write `.claude/plans/diffusion_crosswalk.md` before code —
  (J⁺,J⁻) ↔ (φ_b, J) dictionary `J± = φ/4 ± J/2`, outward-normal sign convention per
  face, group ordering, composite flat layout, per-block conventions.

⏸ **COMPACTION C1** — substrate landed (data seam + trace substrate committed; the
pre-carve agent memos fully absorbed into code + crosswalk + this plan).

### P3 — Geometry law + realizer (#182)
- New `ZeroFluxBoundary` Dirichlet law in `geometry/boundary/` (the mathematical
  idealization; realizer-refused where unrealizable). Audit vacuum-law docstrings for the
  faithful "zero incoming" semantics.
- `DiffusionMethodSpace` (mirror `SNMethodSpace.for_face`) +
  `DiffusionBoundaryRealizer.realize(law, method_space)` replacing the stub
  (`orpheus/diffusion/boundary_realizer.py`): vacuum→𝒜=0 (Marshak), reflective→𝒜=1,
  albedo(α)→𝒜=α, zero-flux→𝒜=−1, white→𝒜=1 (coincides with reflective at P1 —
  document). Composition laws via the walker.
- Flip `test_boundary_realizer_stub.py` negative test → positive realization tests.
  Closes #182.

### P4 — Operator family (the core carve)
- `LeakageOperator` (L): FD `−∇·D∇` on the composite; reads D from P1 seam; emits
  boundary outflow J⁺ into the trace block (mirror SN StreamingOperator's block
  structure). `DiffusionBoundaryOperator` (B): the realized A_ss albedo block.
- S/F: verify carrier-genericity of existing `ScatteringOperator`/`FissionOperator`
  at contact; if SN-shaped, mint scalar-composite siblings wrapping the SAME kernels
  (K_iso, rank-1 dyad) — the kernels are the single source of truth.
- Compose `A = L + C − S − B`; `MatrixInverseOperator(A)` resolvent over composite flat.
- Foundation gates (Mode-12 companions to the k value gates): **object-level stencil
  gate** `A.as_matrix() ≡` independently hand-posed FD matrix on 3–4 cells × 2G,
  mutation-verified RED under D-face swap / Σ_r-vs-Σ_t confusion / scatter transpose;
  per-group SPD (reflective & zero-flux); M-matrix sign pattern; column-sum conservation
  `1ᵀ(C−S) = Σ_a`; `L @ constant = 0` under reflective.

⏸ **COMPACTION C2** — operator family + realizer + foundation gates committed (P3+P4
detail no longer needed in context; APIs are in the tree).

### P5 — Solver + engines
- Modern `solve_diffusion_1d(materials: dict[int, Mixture], mesh, ...)`: mesh+laws →
  `MaterialMesh` → trace → realized B → operators → engines. 3-zone slab = `Mesh1D` with
  3-valued `mat_ids` (replaces `CoreGeometry`).
- `DiffusionSolver` implements the `EigenvalueSolver` protocol boundary with EXACT inner
  solve (the matrix inverse; no scattering inner iteration), driven by the ONE
  `power_iteration`; `direct_eigenvalue` cross-engine gate at 1e-10; adopt
  `compute_production_rate`/`IntegratedReactionRate` (#270 diffusion arm; retires the
  legacy un-normalised window + `fi /= max` conditioning note).
- 3-group discriminator test (kills the flip trick by construction). Demo
  `examples/diffusion/demo_diffusion_1d.py` rewire.

### P6 — Retirements + suite rewire (test migration, not deletion)
- Retire: `TwoGroupXS`, `CoreGeometry`, `BC_REGISTRY` strings + `_resolve_bcs`, scipy
  `A_op` + BiCGSTAB path, `print()` in `converged`. Re-run the 3-search audit per target
  at execution (explorer's table is the baseline; includes `docs/` grep — dangling
  Sphinx refs render silently).
- Derivations (the un-listed 4th consumer): rewire to `derive_2rg_continuous`
  (transcendental successor), DELETE `derive_2rg` Richardson cache + `solver_cases` +
  the "dif" lazy-loader entry; DELETE the H4 self-referencing
  `test_spatial_convergence_reflected` + duplicate `test_spatial_convergence_bare`.
- Continuous-reference suite: REWIRE 4 solver-driving tests to the new API under the
  **zero-flux law** (math unchanged); KEEP 3 reference-only tests; sharpen
  `test_properties.py` (vacuum test now asserts J⁻=0 semantics).
- NEW vacuum(Marshak) references: extend the transfer-matrix derivation with Robin faces
  if cheap; else property-level gates (J⁻=0 satisfied; k strictly between zero-flux and
  reflective) + file the analytic-Marshak-reference follow-up.
- MMS #93: manufactured heterogeneous-D, multigroup-coupled fixed-source gate (Mode-7
  override — the simple single-group sine NULLS the D-gradient and group-coupling terms);
  new `diffusion-mms` theory label. Closes #93 if landed in full.

⏸ **COMPACTION C3** — solver + retirements + suite rewire committed (the island is
gone; only P7's self-contained side-carve + P8 close-out remain).

### P7 — TransportMethod trigger (the recorded 2nd-witness forward trigger)
- Diffusion's functional realizer is the second functional realizer — the recorded
  trigger (memory: `realize_recursively → TransportMethod`) FIRES: mint the
  `TransportMethod` Protocol for BOTH witnesses (homogenization `material_mesh.py:21` +
  boundary-realizer seam), move the walker `sn/ → geometry/boundary/`, registry-route
  leaf dispatch. Execute per `.claude/plans/realize_recursively_move_spec.md` **after L24
  re-characterization** (its premises were verified on an older branch; the fresh-process
  registration gates §2b + layer-import tripwire H2 are the load-bearing hazards).
- USER may strike/defer this phase at plan review without affecting P1–P6.

⏸ **COMPACTION C4** — TransportMethod carve committed (P8 needs only the plan's
close-out list + the archivist brief; docs work is dispatch-heavy, not context-heavy).

### P8 — Docs + pyright + close-out
- Archivist dispatch: `docs/theory/diffusion_1d.rst` overhaul (operator-family
  architecture section, Development-history changelog seeded per the doc template, Key
  Facts update, **BC-law semantics rewrite** — vacuum vs zero-flux naming per ruling 3,
  fix the malformed `:class:` ref at :63), `verification.rst:87` wrong module name,
  `docs/api/diffusion_1d.rst` re-point, `operator_algebra.rst` DSA cross-refs to the
  now-real `A_diff`, stale TODOs (`test_diffusion.py:9-13`, derivations `diffusion.py`).
- Pyright: lift `"orpheus/diffusion"` from `pyproject.toml:119` ignore; add
  `diffusion: 0` to `tests/_harness/pyright_baseline.json`; CLI-verify.
- Full verification (below); ff-merge to main; delete branch.
- Close #290 (+#182, +#93 if landed) via trailers; comment on #2 (A_diff ready, seams
  named), #279 (diffusion arm landed), #270 (diffusion arm landed), #33/#34 (structurally
  unblocked); file follow-ups: physical-P1 re-baseline (ruling 4 part 2), analytic
  Marshak reference (if deferred), mixed-form/RT0 lift seam.

## DSA (#2) seams left behind (named, not built)
- `A_diff` contract: `LinearOperator` on the scalar composite, invertible via
  `MatrixInverseOperator`, consumes the moment-0 reduction of `evaluate_residual`'s
  output (`sn/solver.py:223-289`); correction → 0 at convergence.
- The SN→diffusion boundary restriction = the `M_half` half-range moment under the
  SHARED `G_s` metric (why ruling 2 matters).
- ⚠ #215 trap (`scattering.py:1305-1320`, vv Mode 9): any future accelerator
  FP-invariance gate MUST run an anisotropic config — never the isotropic reflective box.

## Verification (end-to-end)
- Canonical: `.venv/bin/python -O -m pytest` SERIAL (`-p no:xdist --timeout=300
  -p no:cacheprovider`) — tests/diffusion + tests/numerics + tests/transport +
  tests/geometry walls per phase; full `-m "not slow"` tree before merge (0 reds
  expected — the 7-red era is over; any red is a regression).
- Cross-engine: `direct_eigenvalue` ≡ `power_iteration` at 1e-10 on the analytic cases.
- Mutation-verify every new gate's teeth (in-process monkeypatch, NEVER `git checkout`
  on uncommitted files).
- `sphinx-build -W docs docs/_build/html` exit 0 (+ explicit grep for retired-symbol
  cross-refs — `-W` misses coincidentally-resolving stale refs); Nexus auto-reloads.
- `npx pyright` — total ≤ 1 (the accepted #288 transport residual) with diffusion at 0.
- `python -m tests._harness.audit` clean.

---

## CHECKPOINT C1 (2026-07-03, session 1)

**Statuses**: P1 DONE @ `836f424` (data seam; 5 gates; legacy D bit-identical
`[1.52835091, 0.42462845]`). P2 DONE @ `78d1431` (scalar trace substrate).
Branch `feature/diffusion-integration`, 2 commits ahead of main @ `d2a2a0c`.

**P2 deviations from plan (all rationale-carrying, none scope-changing):**
1. `TraceField` locus base EXTRACTED from `BoundaryField` (plan said "mint the
   leaf"; `FullField`'s isinstance check + single-substrate honesty forced the
   two-family hierarchy: `TraceField → {BoundaryField angular, ScalarTraceField
   scalar}`). `_trace_space_of` classmethod hook = the per-family trace source.
2. `PartialCurrentDisplacement` minted (NOT in plan): the #208 torsor registry
   (`Displacement._BY_REP`, keyed by the shared Rep = direct Field base) demands
   a displacement sibling the moment `state ⊖ state` runs — drivers will
   subtract composites, so the scalar role grid had to complete now. Rep key =
   `ScalarTraceField` (why the family ABC exists).
3. Four #289-style seam parses added (`fission.py` composite arm,
   `multiplication_operator.py` apply+solve arms, `sn/solver.evaluate_residual`
   trace legs): pre-existing ANGULAR-ONLY composite arms leaned statically on
   the SN-typed boundary slot; they now parse the family loudly and read the
   mesh off the parsed (SNMesh-typed) leaf. The scalar arms land in P4 through
   these exact seams.
4. `FullField.mesh` now types `MaterialMesh` (method-generic); `TimedFullField`
   overrides widened to `TraceField` (Liskov).

**Convention contract**: `.claude/plans/diffusion_crosswalk.md` (slot layout
`(2, ng, *face_spatial)`; rows J⁺=0/J⁻=1 owned by `ScalarTraceSpace`; albedo
family incl. zero-flux 𝒜=−1; face-local outward-normal signs; area metric).

**Verification state**: `npx pyright orpheus/` = **1** (exactly the accepted
#288 transport residual — baseline unchanged); `tests/{transport,numerics,data}`
1361 green; `tests/{sn,geometry}` 2212 green (serial `-O`) — TraceField
extraction behavior-neutral for SN. Intrinsic suite:
`tests/transport/fields/test_partial_current.py` (15 tests incl. torsor laws).

**Next**: P3 — `ZeroFluxBoundary` law (geometry), `DiffusionMethodSpace`,
`DiffusionBoundaryRealizer` functional (#182): vacuum→𝒜=0, reflective→𝒜=1,
albedo(α)→𝒜=α, zero-flux→𝒜=−1, white→𝒜=1 (P1-coincident, document); flip the
stub test. Then P4 operators (read the crosswalk + agent memos first:
explorer `campaign_290_diffusion_integration_map.md` §B/§C for the typing
shape + mirror pattern; attacker `diffusion_integration_frames.md` Q2/Q3).

**Open questions**: none blocking.

---

## P2.5 — Trace naming coherence (USER RULING 5, 2026-07-03 post-C1; execute FIRST after compaction, before P3)

**The vocabulary partition (binding):** **Trace = the SPACE world; Boundary = the
FIELD world** (matches the composite's `bulk × boundary` intuition), both
family-qualified; **role tokens uniform across families** — the SN trace is
angular-resolved current data under `|Ω·n|w`, so there is no principled ground
for "Flux" in one family and "Current" in the other (user). The state leaf is
`<Family>BoundaryFlux` in BOTH families (the HarmonicMomentFlux precedent: Flux
= the ROLE; the storage — ψ_Γ per-ordinate vs the (J⁺,J⁻) pair — is docstring
content). The rep-keyed `Displacement._BY_REP` pairing mechanism stays as-is;
the renames make its keys read coherently.

**Rename table (single `refactor(transport)` commit, bit-identical):**

| Current (post-P2) | New |
|---|---|
| `TraceSpace` + string `"sn_trace"` | `AngularTraceSpace` + `"angular_trace"` |
| `ScalarTraceSpace` + `"scalar_trace"` | (stays) |
| `TraceField` (locus base) | `BoundaryField` (the locus base RECLAIMS the name, now method-agnostic) |
| `BoundaryField` (angular family) | `AngularBoundaryField` |
| `ScalarTraceField` | `ScalarBoundaryField` |
| `BoundaryFlux` / `BoundarySourceSink` / `BoundaryResidual` / `BoundaryDisplacement` | `AngularBoundaryFlux` / `AngularBoundarySourceSink` / `AngularBoundaryResidual` / `AngularBoundaryDisplacement` |
| `PartialCurrent` | `ScalarBoundaryFlux` (CONFIRMED post-analysis, ruling 5b: Flux = the state ROLE whatever the rep stores — the HarmonicMomentFlux precedent; the DSA restriction then reads state→state in the operator algebra. The partial-current vocabulary LIVES ON in the accessors `outflow_view`/`inflow_view`/`net_current` + docstrings) |
| accessor `boundary_scalar_flux(face)` | `p1_boundary_scalar_flux(face)` (the φ_Γ-misreading fence: `φ_Γ = 2(J⁺+J⁻)` holds ONLY under the P1 closure — the name carries the model-dependence at the one site φ_Γ is derived; `net_current` stays unprefixed because `J = J⁺−J⁻` is EXACT for any angular distribution) |
| `PartialCurrentDisplacement` | `ScalarBoundaryDisplacement` |
| `SNMesh.trace` | `SNMesh.angular_trace` (Q2 Option 1 — pairs with the inherited `scalar_trace`; kills the same-object polysemy before DSA) |
| `MaterialMesh.scalar_trace` | (stays) |
| module files | follow the classes: `trace_space.py`→`angular_trace_space.py`, `partial_current*.py`→`scalar_boundary_flux.py`/`scalar_boundary_displacement.py`, `boundary_flux.py`→`angular_boundary_flux.py`, etc. + test files |

**Also in the same commit:** role-grid docstrings rewritten to the honest 3-axis
grid — locus {Bulk, Boundary(field)/Trace(space)} × family {Angular, Scalar,
Moment} × role {Flux, SourceSink, Residual, Displacement} ("Boundary" stops
being listed as a fourth rep family); `FullField`/`TimedFullField` isinstance
messages + the two `test_timed_full_field` match-pins flip to the reclaimed
`BoundaryField`; crosswalk file updated; the two `#290 P2` memos' class names
(MEMORY.md campaign entry) updated.

**Hazards:** (a) `"sn_trace"` string pins in tests — grep; (b) docs/ xrefs
render silently on stale targets — run the 3-search sweep (graph callers +
text-grep code/tests/docs + direct constructors) per rename target; (c)
`Displacement._BY_REP` keys derive from direct bases — verify
`sibling_of(AngularBoundaryFlux) is AngularBoundaryDisplacement` and the scalar
pair in a probe; (d) `BoundaryFlux` is referenced throughout sn/ + tests —
mechanical, nexus `rename` assists; (e) gates after: full serial walls
(transport/numerics/data + sn/geometry) + `npx pyright orpheus/` = 1 + Sphinx
`-W` if docs touched.

**Why now:** zero external consumers of the P2 scalar names exist yet, and
P3/P4 (realizer, method space, both B operators) multiply references to exactly
these names — clean-before-extend.

**Ruling 5b rationale (2026-07-03 boundary-convention analysis — feed to the P8
archivist for the theory page):** the boundary object is flux-VALUED but
current-MEASURED (transport trace theorem: traces live in
L²(Γ, |Ω·n̂| dΩ dA) — the cosine weight is not optional). SN stores the values
(raw per-ordinate ψ — the sweep's upwind data, uncosined) and keeps the measure
in the space metric `G_s`; the scalar family stores the already-contracted
value·measure products — its components genuinely ARE currents, which is also
why its space metric correctly degenerated to bare AREA. The BCs never
reconstruct a flux: the albedo family `J⁻ = 𝒜J⁺` is closed in current
variables; within the pair `J = J⁺−J⁻` is EXACT (any angular distribution)
while `φ_Γ = 2(J⁺+J⁻)` is P1-closure-dependent (hence the `p1_` accessor
prefix). The UNIFYING frame that justifies ONE family pattern: both families
store the inflow/outflow decomposition of their model's boundary state — SN
per-ordinate (BCs ordinate-diagonal via reflection pairing), scalar
per-half-range (BCs albedo-diagonal); the (J⁺,J⁻) rows are the same structural
split as SN's inflow/outflow ordinate selectors, post-integration. Basis note:
the Marshak pair (crossing-exact bookkeeping, φ ± 2J direction) was chosen over
the P1 CHARACTERISTIC basis φ ± √3J (the Mark conditions ψ(∓1/√3) = 0 — the S2
Gauss points): conservation-natural over PDE-characteristic, right for BC
bookkeeping + DSA current-consistency.

### P2.5 status (2026-07-03, session 2 — EXECUTED)

The full rename table + role-grid rewrites + message/pin flips landed as ONE
`refactor(transport)` commit (hash = this commit). Deviations from spec: NONE.
Scope notes recorded for the retirement audit: `SNMethodSpace.trace` KEPT
(per-method protocol surface, no polysemy — P3's `DiffusionMethodSpace` mirrors
it; only `SNMesh.trace` was the polysemous member); `FullFieldSpace.trace_space`
attr KEPT (generic over both families); untracked user diag scripts in
`derivations/diagnostics/` untouched (still spell old names — user scratch).
Gates: transport+numerics+data wall 1361 passed; sn+geometry wall + `sphinx -W`
+ CLI pyright = 1 (#288 residual) + `tests._harness.audit` EXIT=0 (9 MISSING
ERR entries pre-existing, 54/63) — full results in the session log.
Probes: `Displacement._BY_REP[AngularBoundaryField] is AngularBoundaryDisplacement`
+ scalar pair ✓; hierarchy `{Angular,Scalar}BoundaryField < BoundaryField(Field)`,
scalar NOT under angular ✓; torsor `ScalarBoundaryFlux ⊖ → ScalarBoundaryDisplacement`
✓; accessor fence `p1_boundary_scalar_flux` present / old name gone ✓; composite
admits both families through the reclaimed `BoundaryField` annotation ✓;
`SNMesh.angular_trace` live, `.trace` attribute gone, space name string
`"angular_trace"` ✓.
Tooling hazards hit + recovered (for future sweeps): (1) git-pathspec `*`
matches ACROSS `/` in `git ls-files` — a "top-level only" glob returned the
whole tree and re-ran the NON-IDEMPOTENT rename pipeline (double-bumped the
reclaimed `BoundaryField` → hierarchy collapse); recovered via
`git restore -- <trees>` (index held post-`git mv` pristine content), then ONE
pass over ONE definitive list. Rename pipelines with reclaim steps are
NEVER idempotent — never re-run over transformed files. (2) zsh builtin
`echo` interprets `\b` as a BACKSPACE byte when appending a perl rule —
the appended article rule silently never matched; use `print -r` /
heredocs for regex text.

## CHECKPOINT C1.5 (2026-07-03, session 2 — post-P2.5, pre-P3)

- **Statuses:** P1 @ `836f424`, P2 @ `78d1431`, C1 @ `5355057`, rulings 5/5b
  @ `c1321b0`+`1353030`, **P2.5 @ `1cd8d32`** (full rename table + role
  grids + gates — see "P2.5 status" block above), diagnostics triage @
  `6d2035c` (7 session probes promoted to tracked post-P2.5-name
  diagnostics, 6 captured ones retired; #282/#283/#276 carry the tracked
  probe paths as comments).
- **Deviations from plan:** none. Scope notes (SNMethodSpace.trace kept,
  FullFieldSpace.trace_space kept) in the P2.5 status block.
- **Open questions:** none for P3 — its spec is self-contained (§P3:
  ZeroFluxBoundary law in geometry/boundary/, DiffusionMethodSpace
  mirroring SNMethodSpace.for_face, functional DiffusionBoundaryRealizer
  replacing the stub, albedo table vacuum→𝒜=0 / reflective→𝒜=1 /
  albedo(α)→𝒜=α / zero-flux→𝒜=−1 / white→𝒜=1-documented-P1-coincident,
  flip test_boundary_realizer_stub.py negative → positive, Closes #182).
- **Re-anchor for the fresh context:** this plan + `git log --oneline -8`
  + the L17 crosswalk (`.claude/plans/diffusion_crosswalk.md`, now in
  post-P2.5 vocabulary). The live trace vocabulary is in the P2.5 status
  block and MEMORY.md's campaign entry. Branch
  `feature/diffusion-integration`, 8 ahead of main @ `d2a2a0c`.
