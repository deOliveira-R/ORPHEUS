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
