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

### P3 status (2026-07-03, session 3 — EXECUTED; hash = this commit)

Spec executed in full (law + method space + functional realizer + stub-test
flip; `Closes #182` fires at push). Deviations from the spec letter: NONE.
Design decisions beyond the spec letter (recorded for P4/P7):

1. **`stamp_boundary_role` extracted to `geometry/boundary/_realizer.py`**
   (exported): the BlockRole.BOUNDARY instance-stamp became a 2-producer
   shared concept at the second functional realizer — SN's module-private
   `_as_boundary` retired, all 9 SN call sites rewired (Cardinal Rule 2).
2. **"Composition laws via the walker" = a `realizer` parameter on
   `realize_recursively`** (default `None` → `SNBoundaryRealizer()`, every
   SN call site bit-unchanged). The diffusion path composes via
   `realize_recursively(tree, dms, realizer=DiffusionBoundaryRealizer())`.
   This is deliberately the MINIMAL bridge: the full generalization
   (walker moved to `geometry/boundary/`, registry-routed leaves, typed
   `MethodSpace` Protocol) stays P7 behind the L24 re-characterization
   gate. The walker's deferral comment now records the fired trigger.
3. **Realizer shape = two named maps** (`law → 𝒜` via
   `_partial_current_albedo`, `𝒜 → operator` via `_albedo_operator`) —
   ruling 2 as code: every diffusion BC IS an albedo-family scalar; the
   structure-keyed collapse (0→Zero, 1→Identity, else 𝒜·I) mirrors the SN
   albedo branch exactly.
4. **Refusals**: SN refuses `ZeroFluxBoundary` with a dedicated
   physics-message branch (negative angular inflow unrepresentable →
   redirect to VacuumInflow / diffusion realizer); diffusion refuses
   `PeriodicBoundary` (opposite-face trace-block wrap — a P4
   boundary-operator-assembly seam, not a per-face albedo) and
   `PrescribedInflow` (rank-0 AFFINE law = boundary SOURCE, not the linear
   B operator — the P5 fixed-source arm; same operator/source split SN
   keeps un-stamped).
5. **Vacuum-docstring audit result**: the geometry law docstrings were
   already faithful ("no incoming flux, irrespective of what leaks out");
   THE unfaithful naming is the island's `_diff_bc_vacuum`
   (`orpheus/diffusion/solver.py:109-114` — "Zero-flux (Dirichlet φ=0)"
   registered under `"vacuum"`) — retired + re-attributed in P6.
   `docs/theory/boundary_conditions.rst` "one functional + four stubs"
   prose is now stale → P8 archivist brief (module docstrings already
   updated in-tree).

Gates: diffusion+geometry 320 / transport+numerics+data 1361 / sn 1937 —
all green (serial `-O`); CLI pyright = 1 (the accepted #288 residual);
`sphinx -W` exit 0 (matrix.rst regenerated, 6141 collected); harness audit
EXIT=0 (9 MISSING ERR pre-existing, 54/63). NEXT = P4 (operator family; read
the crosswalk + explorer memo §B/§C + attacker memo Q2/Q3 first per C1).

### P4 status (2026-07-03, session 3 — EXECUTED; hash = this commit)

Spec executed in full: `A = L + C − S − B` composes on the scalar composite,
`MatrixInverseOperator` resolvent works over the composite flat, all
foundation gates + the Mode-12 stencil gate landed with 4 RED mutations.
Design decisions beyond the spec letter:

1. **`FlattenedOperator` minted** (`orpheus/numerics/flat_operator.py`) —
   the composite→flat serialization bridge (`flat ∘ A ∘ unflat` over the
   duck-typed `to_flat`/`from_flat` pair). The resolvent spelling is
   `MatrixInverseOperator(FlattenedOperator(A, template))`; `as_matrix`
   derives its basis dimension for FREE because `FullFieldSpace.shape` IS
   the flat direct-sum dim. Unparameterized (ndarray-carrier family
   convention — the TypeVar bounds reject bare ndarrays).
2. **S = composite arms ON the K_iso kernels** — no wrapper type minted
   (type-vs-property: `OperatorSum` already expresses `S_iso + S_n2n`).
   One shared `_scalar_composite_source` helper; each kernel's composite
   arm routes through its OWN bare `apply` (single source). Composite
   `apply_transpose` REFUSED with a #281 pointer (no consumer yet).
3. **C scalar arms via the degenerate 1-ordinate engine lift**
   (`engine.apply(values[None])[0]`) — bit-identical to the direct
   multiply, ONE broadcast engine + its frozen spectrum gate for both
   families (probe-verified). F's composite arm widened to ScalarFlux
   bulk, reusing its ScalarFlux branch verbatim.
4. **`ScalarBoundarySourceSink` minted** (consumer-driven: the L/B
   codomains demanded it now) — thin role leaf, SCALAR_FLUX_UNITS
   ("the trace is all-current"), closing the loss sum as
   `ScalarSourceSink ⊕ ScalarBoundarySourceSink` composites.
5. **`MaterialXSField.diffusion_coefficient`** per-cell channel (the P1
   seam's spatial gather, single-source) +
   **`MaterialMesh.scalar_full_field_space`** (bulk metric = V_cell
   broadcast over groups; method-agnostic — an SNMesh inherits BOTH
   composites, exactly what the DSA restriction #2 wants).
6. **L's blocks mirror SN streaming exactly**: bulk = conservative
   `(A·J)/V` divergence with condensed interior currents (series
   half-cell resistance = the RT0/harmonic form ≡ the island's
   arithmetic-mean-of-Σtr in reals); edge rows read the trace's net
   outward current; trace block = outflow-definition defect
   (`J⁺ − c_φφ_e − c_JJ⁻`, from Fick + the P1 dictionary,
   ρ = h/(2D), c_φ = 1/(ρ+2), c_J = (ρ−2)/(ρ+2)) + the inflow IDENTITY,
   so `(L−B)` inflow rows read `J⁻ − 𝒜J⁺`. Coord-general by
   construction (face_labels-driven; a curvilinear pole is not a face →
   zero flow slot). Multi-D refused at construction. The two
   discretization kernels are module-level pure functions
   (`_interior_conductance`, `_boundary_closure`) — mutation-patchable.
7. **Verified equivalences** (session probes + gates): the
   Schur-condensed zero-flux closure ≡ the island's `φ/(0.5·dz)` vacuum
   arm at 1e-12; `(L − B_reflective)` annihilates `[φ₀, φ₀/4]` (bulk
   bit-exact, closure row ULP); matrix action ≡ typed action at 1e-12.
8. **Gate file** `tests/diffusion/test_operators.py` (32 tests):
   stencil gate `[A] ≡ hand-posed` (plain-loop independent realization;
   heterogeneous 2-mat, non-uniform 4-cell, 2G, ASYMMETRIC Σs) across
   zero-flux/reflective+Marshak/albedo(0.3) configs; 4 mutations RED
   (D-face pairing swap, Σ_a-for-Σ_t, scatter transpose via the
   MaterialXSField kernel swap, closure c_J sign flip); column-sum
   theorem `1ᵀ(C−S) = Σ_a` (removal DERIVED); M-matrix sign pattern;
   per-group SPD of the bulk Schur complement **under the volume
   metric** (`diag(V)@block` — the raw FD row-scaling is V-similar on
   non-uniform h; the metric IS the composite bulk weights); resolvent
   M-materialise both ways + seed independence; substrate + arm pins
   (each composite arm ≡ its own bare kernel).

Gates: diffusion+geometry 352 / transport+numerics+data 1361 / sn 1937 —
all green (serial `-O`); CLI pyright = 1 (#288 residual; the 2
FlattenedOperator TypeVar errors fixed by unparameterizing); `sphinx -W`
exit 0; audit EXIT=0 (54/63 pre-existing). NEXT = P5 (solver + engines:
`solve_diffusion_1d` on mesh+laws → MaterialMesh → realized B → operators
→ `direct_eigenvalue` ⊗ `power_iteration` cross-gated 1e-10 +
`compute_production_rate` #270 arm + 3-group discriminator + demo rewire).

## CHECKPOINT C2 (2026-07-03, session 3 — post-P4)

- **Statuses:** P1 `836f424` / P2 `78d1431` / P2.5 `1cd8d32` / P3
  `6672e7a` / **P4 = this commit** (statuses above). Branch
  `feature/diffusion-integration`, 11 ahead of main @ `d2a2a0c`.
- **Deviations:** none scope-changing; the 8 design decisions in the P4
  status block (FlattenedOperator mint, K_iso-arms-not-wrapper, engine
  lift, ScalarBoundarySourceSink, D channel, scalar composite space, L
  block conventions, V-metric SPD).
- **Open questions:** none for P5 — its spec is §P5 + the P4 status
  block's NEXT line. The island (`solver.py`) is UNTOUCHED so far;
  P5 builds the modern path beside it, P6 retires it.
- **Re-anchor:** this plan + `git log --oneline -12` + the crosswalk.
  The P4 API surface: `orpheus.diffusion` exports `LeakageOperator`,
  `DiffusionBoundaryOperator` (+ P3's realizer/method-space);
  `orpheus.numerics.flat_operator.FlattenedOperator`;
  `mesh.scalar_full_field_space`; `mat_xs.diffusion_coefficient`;
  `ScalarBoundarySourceSink` in transport.source_sinks.

### P5 status (2026-07-03, session 3 — EXECUTED; hash = this commit)

**DONE.** The modern solver + engines, built BESIDE the untouched island
(P6 retires it). New module `orpheus/diffusion/k_eigenvalue.py`:
`DiffusionSolver` (the `EigenvalueSolver` + `ProductionRateSolver`
implementer), modern `DiffusionResult` (typed composite flux + current
profile + history), `solve_diffusion_1d(materials, mesh, ...)` driver on
the ONE `power_iteration`. Design decisions (the P5 record — full
rationale in the module docstring):

1. **Protocol carrier = the flat composite vector** — the trace (J⁺,J⁻)
   is honestly part of the eigenvector; conversion at exactly two sites
   (`unflatten` + the template-frozen FlattenedOperator).
2. **Exact inner solve** — `MatrixInverseOperator(FlattenedOperator(A,
   template))`, one eager LU, one back-substitution per outer; NO
   scattering inner iteration exists.
3. **k-update = the integrated eigenvalue relation** `k = P(ψ)/⟨1,
   (Aψ)_bulk⟩_V` — decomposes to absorption + leakage by the P4
   column-sum theorem + divergence telescoping; derived through THE loss
   operator so no term can be forgotten (the island hand-added leakage;
   SN's P/A omits it — flagged, see follow-up note below).
4. **(n,2n) loss-side** (homogeneous convention): S = IsotropicScattering
   + IsotropicN2N; F and `compute_production_rate` are νΣf ONLY.
5. **BCs resolve from the MESH's declarations** — tag→law→realize loop
   over `face_labels` mirroring `SNMesh._resolve_bcs` (None → reflective
   default); registry: vacuum/reflective/albedo/zero_flux; "white"
   DELIBERATELY absent (≡ reflective at P1). This is the 3rd spelling of
   the per-method resolution pattern — P7's TransportMethod unifies.
6. **Normalization**: returned mode carries ∫νΣf·φ dV = 1 (#270 arm;
   `fi /= max` conditioning + hardcoded e_per_fission power window
   retired from the modern path).
7. **`LeakageOperator.face_currents(psi)`** added: the current-profile
   reconstruction (interior condensed + boundary trace, axis-signed);
   `apply` now CONSUMES it (single source — flow = areas × currents,
   bit-identical, P4 stencil gates re-verified green); dead
   `_FaceClosure.area` field retired.
8. **Package exports flipped to modern** — island importable ONLY via
   `orpheus.diffusion.solver` submodule for one merge cycle (verified:
   all island consumers — tests, derivations `diffusion.py:152`, demo —
   already import the submodule).

Gates (`tests/diffusion/test_solver.py`, 22 tests): cross-engine
`power_iteration` ≡ `direct_eigenvalue` |Δk| < 1e-10 + full composite
eigenvector, 4 BC configs (foundation); **L2 infinite-medium** reflective
diffusion k ≡ homogeneous k∞ at 1e-11, 2G + **3G asymmetric (the
flip-trick discriminator — skip-scatter 1→3, split χ; ng-generic by
construction)**; **legacy bridge** modern ≡ island on CORE1D at 1e-8
(island driven through power_iteration directly — its driver's
max_iter=200 does NOT converge the PWR dominance ratio; the observed
9.5e-7 first-run delta was island non-convergence, not discretization);
trace semantics per law at the solution (vacuum J⁻=0, albedo J⁻=αJ⁺,
reflective J_net=0, zero-flux φ_Γ=0, all LU-exact 1e-12·scale) + the
asymmetric-BC face-mapping pin; balance identity P/k = absorption +
leakage at 1e-10; k monotone in 𝒜 (zero_flux < vacuum < albedo(0.6) <
reflective); production-normalization pin; ProductionRateSolver
isinstance; refusals (white with supported-list, albedo-sans-param,
multi-D — L constructs FIRST in the solver so ITS 1-D refusal fires, not
the trace's AttributeError).

**Mutation teeth probed in-process** (monkeypatch + revert, not
committed): M1 volume-weight drop in the k-update → balance RED 4.3e-1 +
cross-engine RED 2.7e-2; M2 registry swap vacuum→reflective → vacuum
trace gate RED at full scale; controls ~1e-15. (The fission /k scale is
spectrally invisible under production renormalization — genuinely
non-load-bearing, no gate possible or needed; the cross-engine gate is
the committed catcher for ALL protocol plumbing.)

Demo `examples/diffusion/demo_diffusion_1d.py` rewired to the modern API
(BC("zero_flux") — the ruling-3 honest name for the MATLAB "vacuum"):
**keff = 1.022173 = the MATLAB reference at print precision** (391
outers; the island's 200-cap driver never converged this).

Gates: diffusion+geometry 374 / transport+numerics+data 1361 / sn 1937 —
all green (serial `-O`); CLI pyright = 1 (#288 residual only); `sphinx
-W` exit 0; audit EXIT=0.

Follow-up flags for later phases: (a) P6 — `git mv k_eigenvalue.py →
solver.py` after the island deletes (family naming parity with
sn/cp/homogeneous); the bridge test class DELETES with the island. (b)
SN's `compute_keff = P/A` omits leakage — verified real-but-latent (all
tight SN k anchors are reflective; vacuum configs smoke-banded only) —
**FILED as #291** (NOT #290 scope).

NEXT = P6 (retirements + suite rewire + MMS #93; ⏸ C3 after).

---

## CHECKPOINT C3 (2026-07-03, session 4 — post-P6)

- **Statuses:** P1 `836f424` / P2 `78d1431` / P2.5 `1cd8d32` / P3
  `6672e7a` / P4 `db14643` / P5 `9470266` / **P6 = this commit**
  (carries `Closes #93` — fires at push, like #182 in P3). Branch
  `feature/diffusion-integration`; user holds pushes. Remaining: P7
  (user may strike/defer) → P8 close-out.
- **Re-anchor:** this plan + `git log --oneline -14` + the crosswalk.
  Post-P6 API surface: `orpheus.diffusion.solver` IS the modern module
  (island gone); `mixture_from_diffusion_tables` in
  `derivations/common/xs_library.py` is the ONE ruling-4 encoder
  (tests + demo consume it); `tests/diffusion/` = boundary_realizer /
  continuous_reference / mms / operators / properties / solver.

### P6 status (EXECUTED; hash = this commit)

**DONE — the island is gone; the modern path carries every anchor.**

Retirements (3-search audit re-run: graph + grep code/tests/docs +
constructors — consumers matched the explorer table exactly, plus
docstring-only mentions):

1. `orpheus/diffusion/solver.py` (island: `CoreGeometry`, `TwoGroupXS`,
   legacy result/solver/driver, `BC_REGISTRY` strings + `_resolve_bcs`,
   scipy `A_op` + BiCGSTAB, `print` in `converged`, hardcoded
   `e_per_fission`) — deleted; `git mv k_eigenvalue.py → solver.py`
   (family naming parity; `__init__`/docstrings/api page updated;
   diffusion's BC_REGISTRY exit leaves the legacy family at 3:
   cp/mc/moc).
2. Derivations legacy 2-region tier: `derive_2rg` (Richardson, ran the
   island on cache miss) + `solver_cases` deleted;
   `_richardson_cache.{py,json}` deleted WHOLE (sole consumer was
   derive_2rg; the `moc_*` JSON entries were already-orphaned data —
   moc.py never imports the cache).
3. `reference_values.py`: the ENTIRE `_load_solver_cases` machinery
   retired (deviation — plan ordered only the "dif" entry, but after
   diffusion's exit NO module defines `solver_cases`; the loop was dead
   code). Every registry entry is now analytical/semi-analytical.
4. Tests: `test_diffusion.py` deleted whole (H4 self-referencing
   reflected + duplicate bare); richardson cross-check test deleted
   (its migration-window job complete);
   `test_solver.py::TestLegacyBridge` deleted with the island (its
   1e-8 equivalence evidence recorded in the solver module docstring).

Rewires:

- 4 L1 anchors → modern API under `BC("zero_flux")` at UNCHANGED
  tolerances — all green (the math never moved; ruling 3 executed as
  re-attribution: ProblemSpec strings "vacuum"→"zero_flux", helper
  rename `_solve_2region_vacuum→zero_flux_eigenvalue`, docstrings,
  test rename `test_2region_zero_flux_boundaries_satisfied`).
- `mixture_from_diffusion_tables` minted in xs_library (ruling-4
  encoding, ONE source; ng-generic chain-downscatter; demo's inline
  bridge + the tests' copies retired onto it; demo re-verified
  keff = 1.022173 = MATLAB at print precision, 391 outers).
- `test_properties.py` sharpened: Marshak gate (BC("vacuum") ⇒ J⁻=0 at
  1e-12·scale + boundary-cell flux strictly positive — vacuum ≢
  zero-flux), positivity, symmetry (+ mirror-equal outward net
  currents pin).
- **Data correction (deviation):** derivations `_REFL_XS`
  `chi=[1,0]→[0,0]` — a placeholder on a NON-fissile region that the
  honest Mixture guard refuses; inert in every consumer
  (χ⊗νΣf ≡ 0 with production ≡ 0) ⇒ reference values bit-identical.
- Docs (P6-minimal; P8 archivist owns the overhaul): api page
  re-pointed at the 4 modern modules; verification.rst Richardson
  section rewritten as an explicit HISTORICAL note (P6 deleted the
  machinery it described in current tense — not deferring a lie);
  `richardson-diffusion` label kept as the technique's record (now an
  orphan equation in the audit — P8 may rule further);
  `docs/_generated/` regenerated (NOTE: git-ignored + MANUALLY
  generated via `python -m orpheus.derivations.generate_rst` — the
  on-disk artifacts were stale from a pre-P6 run; sphinx does NOT
  regenerate them).

**MMS #93 landed in full** (`tests/diffusion/test_mms.py` + theory
section `diffusion-mms-section` + label `diffusion-mms`): Mode-7
override executed — heterogeneous per-group D(x) (every cell its own
material ⇒ every interior face exercises the conductance
interpolation) + multigroup-coupled distinct shapes; SymPy-exact
forcing; exact-resolvent fixed-source solve; layout self-check via
`FullField.from_flat` round-trip. Orders **2.004 / 2.001 / 2.000**;
COMMITTED Mode-10 controls: flattened-D ×144 and coupling-sign-flip
×166 above the clean floor (3.4e-3 @ n=40). Production-side probe M3
(`tmp/probe_p6_mms_mutation.py`, monkeypatch `_interior_conductance`
→ one-sided): order collapses to 1.121/1.034/1.009, finest 3.4e-3 —
RED under both gate assertions; control clean.

**Marshak-reference decision: else-branch** (plan §P6) — the analytic
transfer-matrix Robin extension (rows `pL/4 + JL/2 = 0`,
`pR/4 − JR/2 = 0` + law-aware physical validation) DEFERRED; the
property gates required by the else-branch already exist (P5
trace-semantics J⁻=0 LU-exact, k-ordering zero_flux < vacuum <
albedo(0.6) < reflective, + the sharpened Marshak property gate).
File the analytic-Marshak-reference follow-up at P8 close-out
(already on P8's list).

Gates: walls diffusion+geometry+**derivations** 1921 (37 min — the
derivations suite's own weight incl. the `6d2035c` promoted probes;
exit 0) / transport+numerics+data 1361 / sn 1937, all green serial
`-O`; sphinx `-W` exit 0 (post-regenerate); CLI pyright = 1 (the
accepted #288 residual; streamed LSP noise ignored per #226); audit
exit 0 (`diffusion-mms` = 2 tests); demo = MATLAB reference.

**Open questions:** none blocking P7/P8. P7 remains user-strikeable;
its trigger evidence is now THREE spellings of tag→law→realize
resolution (homogenization, SN, diffusion `_resolve_bcs`).

NEXT = P7 (TransportMethod mint — execute per
`realize_recursively_move_spec.md` after L24 re-characterization;
USER may strike/defer) → P8 (docs overhaul + pyright lift + close-out:
close #290/#182/#93, comment #2/#279/#270/#33/#34, file follow-ups:
physical-P1 re-baseline, analytic Marshak reference, RT0 lift seam).
⏸ C4 after P7 (or straight to P8 if struck).

---

## P7 REDESIGN (2026-07-03, session 4 — USER RULING at the design discussion; SUPERSEDES the §P7 section above)

**Ruling.** `TransportMethod` is a Protocol implemented by the
**method-mesh layer** — `SNMesh` and a new `DiffusionMesh` — NOT by
stateless per-method singletons and NOT by `DiffusionSolver`. Diffusion
is missing its `augmented_mesh.py`; **create `DiffusionMesh` FIRST
(P7a)**, augmenting `MaterialMesh` with the diffusion-specific
behavior, then mint the Protocol over the two now-symmetric
method-meshes (P7b). This matches witness #1's original wording
(`material_mesh.py` module docstring: "a method-specific mesh IS a
MaterialMesh that adds behavior … conforming structurally to a future
TransportMethod Protocol … minted when a 2nd method-mesh exists") —
fix the missing piece, don't contort the abstraction around its
absence.

**The P4 contradiction this resolves.** `material_mesh.py`'s module
docstring declares the axis (MaterialMesh = mesh+materials DATA; the
method layer — quadrature, stencil, **boundary trace**, closures — is
BEHAVIOR), yet the P4-added property at `material_mesh.py:479` claims
`scalar_full_field_space` is "Method-agnostic (like scalar_trace)".
Same file, opposite claims. P4 had no diffusion method-mesh to put
them on. DiffusionMesh reclaims both — **MaterialMesh should not know
what a trace is** (SNMesh knows angular traces; DiffusionMesh knows
scalar traces). Type-vs-property check passes: 2 non-isomorphic
augmentations (quadrature-machinery vs trace+realized-BCs), real
morphisms (promotion, BC realization, solver/operator/field
consumption), and "a diffusion phase space with unresolved BCs"
becomes unrepresentable.

### P7a — the DiffusionMesh carve (surgical; do FIRST)

**STATUS: DONE @ `738e355` (2026-07-03, session 5).** Executed per the
carve list below, zero deviations of substance. Notes for P7b: (1) the
solver's `_resolve_bcs`/`_law_from_tag` moved onto `DiffusionMesh`
with the SN split kept (`_resolve_bcs` loop + `_resolve_one` +
`_law_from_tag`); the two remaining per-method loops are now
shape-identical for the P7b shared body. (2) Extra admission gate
beyond the plan: mesh-less promotion (the infinite-medium 1-cell
carrier) is refused with a points-at-homogeneous-solver message —
`self.areas` would otherwise die with a misleading AttributeError.
(3) `docs/api/geometry.rst`'s BC-resolution paragraph + diffusion
table row were island-stale AND invalidated again by the solver→mesh
move — rewritten now (method-mesh registries); the theory page's
island-era BC prose (:62-93) left for P8's planned BC-semantics
rewrite. (4) Gates: walls 379/1362/1937 serial -O (sn count ≡ P6);
demo = MATLAB 1.022173; pyright CLI = 1 (#288 only); sphinx -W 0;
audit 0; 3-search audit clean (grep + Nexus graph + built HTML).

New `orpheus/diffusion/augmented_mesh.py` (flat file — name-parity
with `sn/mesh/augmented_mesh.py`; subpackage only when diffusion/
grows), `class DiffusionMesh(MaterialMesh)`:

1. **`from_material_mesh(mm) -> DiffusionMesh`** classmethod (the
   promotion arm). NO extra params — BCs come from the axes' own `BC`
   tags (contrast `SNMesh.from_material_mesh(mm, quadrature, scheme)`;
   honest asymmetry, promotion signatures don't unify and shouldn't).
2. Construction mirrors `SNMesh.__init__` order: the **1-D refusal
   moves here** (from LeakageOperator-by-ordering — DELETE the P5
   "L constructs FIRST" hack + comment, solver.py:249-253); build
   `scalar_trace` (constructor body from the reclaimed
   `MaterialMesh.scalar_trace`, material_mesh.py:424-456; NAME KEPT —
   locus-family parity `SNMesh.angular_trace`/`DiffusionMesh.scalar_trace`);
   build **`full_field_space`** (reclaimed from
   `MaterialMesh.scalar_full_field_space` :458-494, **RENAMED** for
   `SNMesh.full_field_space` parity — the future Protocol member);
   realize **`bc: dict[str, LinearOperator]`** at construction
   (`SNMesh.bc` parity, augmented_mesh.py:452) — move
   `_resolve_bcs` + `_law_from_tag` + `BOUNDARY_OPERATOR_REGISTRY`
   OFF `DiffusionSolver` (solver.py:230-327).
3. **MaterialMesh reclamation**: DELETE `scalar_trace` +
   `scalar_full_field_space` properties (material_mesh.py:424-494);
   fix the module/property docstrings (drop the "method-agnostic"
   P4 claims; note the axis restoration). Blast radius = the grep
   captured 2026-07-03 (session 4): diffusion solver.py:254/:296,
   operators.py:318/:322/:444/:461/:465, method_space.py docstrings,
   transport/fields/_bases.py:988-992 + the 3 scalar leaf modules'
   docstring contracts, tests/diffusion/test_operators.py (+:613-617),
   tests/transport/fields/test_scalar_boundary_flux.py:50-242.
4. **Solver diet** (solver.py): `__init__(mesh: DiffusionMesh)`;
   `space = mesh.full_field_space`; `DiffusionBoundaryOperator(mesh)`
   reads `mesh.bc` (drops the realized-dict param); driver
   `solve_diffusion_1d`: `MaterialMesh(mesh1d, materials)` →
   `DiffusionMesh.from_material_mesh` → solver (PUBLIC SIGNATURE
   UNCHANGED). Solver = operators + engines only (the SNSolver split).
5. **Operators** take `DiffusionMesh` (IS-A keeps MaterialMesh-typed
   internals working); `domain`/`codomain` read `mesh.full_field_space`.
6. **Scalar field family** (transport L2): repoint the TYPE_CHECKING
   mesh annotation `MaterialMesh → DiffusionMesh` — the ANGULAR-family
   precedent (`_bases.py:685` annotates SNMesh from L2; the layer gate
   TC-exempts numerics/transport); duck read `mesh.scalar_trace`
   unchanged at runtime; update the `MaterialMesh.scalar_trace`
   docstring refs.
7. `DiffusionMethodSpace.for_face`: drop the explicit `trace=` param
   (read off the mesh) — adopted.
8. **Tests**: test_operators.py + test_scalar_boundary_flux.py
   construct DiffusionMesh; TestProtocolAndRefusals `pytest.raises`
   moves to MESH construction (albedo-sans-param, unsupported-kind,
   multi-D — sharper surface); test_solver.py direct
   `DiffusionSolver(MaterialMesh(...))` sites → DiffusionMesh; suites
   driven through `solve_diffusion_1d` unchanged. NEW intrinsic tests:
   DiffusionMesh construction laws (bc realized per-face, default
   reflective, trace/space identity, promotion round-trip,
   refusals-at-construction).
9. **Gates (bit-identical carve — zero value drift expected):** the
   full P4–P6 diffusion wall re-runs green (stencil, cross-engine,
   trace semantics LU-exact, L2 k∞, balance, MMS orders
   2.004/2.001/2.000, demo = MATLAB 1.022173); walls
   diffusion+geometry / transport+numerics+data / sn; CLI pyright = 1;
   sphinx -W 0; audit 0. Re-run the 3-search audit on the two
   reclaimed properties at execution.

### P7b — TransportMethod Protocol (after P7a; the original §P7 payload, updated)

- Protocol in `transport/method.py` (L2; conformance is STRUCTURAL —
  sn/diffusion never import it; consumers type against it via
  TYPE_CHECKING). **Instance surface only** (promotion stays
  per-method-typed classmethods). Member list minted at carve time
  from what the witnesses actually consume — candidates:
  `full_field_space`, `bc`, the law-admission registry surface,
  method-space-for-face; `name` for error text ONLY (never dispatch).
- ONE shared `resolve_boundary_conditions` generic body replaces the
  twin `SNMesh._resolve_bcs`/`DiffusionMesh` loops (the 3rd spelling
  died at P7a when the solver's copy moved to the mesh).
- Walker `realize_recursively` moves `sn/boundary/realizer.py:317` →
  `geometry/boundary/`, typed against the existing `BoundaryRealizer`
  Protocol (geometry stays method-blind).
- **`BoundaryRealizerRegistry` + string keys + the MoC/MC/CP
  NotImplementedError stub realizers DISSOLVE** (aggressive
  retirement; you hold the method-mesh → you have its realizer). The
  banked spec's dominant hazards (H1/H3 registration timing, §2b
  fresh-process gate) are import-side-effect artifacts — the dissolve
  DELETES the hazard class; §2b becomes moot.
- Banked verification structure reused (§1 composition wall rewires
  its import; §3 layer-imports tripwire = the structural gate; H4
  doc-xref grep). L24 re-characterization: re-measure the spec's
  pre-state table (stale: pyright 412→1, the 7-red era is over, the
  walker's file moved in the sn/ reorg).
- Homogenization witness #1 discharges: `material_mesh.py:21`'s
  "future TransportMethod Protocol" comment updates to point at the
  real Protocol.

⏸ **COMPACTION C4** after P7a+P7b commit (P8 as before).
