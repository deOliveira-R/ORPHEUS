# THE projection / reconstruction machinery — the unified Frame campaign

> **⏹ CAMPAIGN COMPLETE (2026-07-26) — P1 through P7 ALL LANDED; this file is the archaeology.**
> P1 `00f9b76` (Closes #268) · P2 `cc6d022` · P3 `79e2ac8` · P4 `74378d5`+`ed1e14d`+`83aaa8a` ·
> P4.5 (the two-param split, commits W-C…W-G) · P5 (condensation, #274 `68ceb9a`) ·
> **P6 `2cbe9da7..e301b675` (2026-07-26; #281 CLOSED)** · **P7 = the docs capstone, delivered
> INCREMENTALLY by the P3/P5/P6 archivist doc-passes** — the 2026-07-26 gap audit found the full
> charter present at the articulation standard (`docs/theory/foundations/frame.rst` self-identifies
> as the P7 capstone; every item Feynman-grade; `-E -W` clean; all cross-refs resolve). **Ruling
> kept from the audit: condensation stays a SECTION of frame.rst (one PG page, space + energy as
> sibling sections sharing the machinery), never a standalone duplicate page.**
>
> ⚠ The "Tasks #46–#52" numbering below is the ORIGINATING SESSION'S internal task list, NOT
> GitHub issues (GitHub #46–#52 are Thermal-Hydraulics/Kinetics). The real trackers were the
> phase issues #268 / #226 / #274 / #281 / #275 — all closed or standing on their own.
> Durable rulings live in `[[project-frame-projection-machinery]]` (memory) + the frame.rst /
> operator_algebra.rst theory pages; the per-phase records below are the historical log.

## Context — why this campaign exists

The campaign has been deferred three times under three names (`typed_carrier_grid_carve` →
`frame_basis_carve` → the homogenization carve), each landing **foundation** but never the **payoff**.
The payoff, in the user's words: a consumer writes **"define Frame, analyse, done"** — a couple of
lines, **nothing hand-rolled** (no membership matrices, no per-channel einsums, no inline Gram-inverse
divides in a method body). Today that is violated everywhere:

- **Scattering** hand-chains `R.apply(Λ.apply(M·ψ))` across two methods; the typed composed operator
  (`kernel = R∘Λ∘M`) exists but is **verification-only** ("semantic, not a rewiring"); a `cast(LinearOperator, Λ)`
  papers over `LegendreMomentScattering`'s `None` spaces.
- **Homogenize** gathers cross-section channels one-by-one and reassembles `Mixture`s by hand (1-D only).
- **Condensation** doesn't exist.
- The Frame **faces consume/emit bare `np.ndarray`** — every consumer unwraps `.values` and re-wraps.
- The `projection.py` Galerkin/PG discipline ABCs are **dead** (zero subclasses).

**Two load-bearing decisions (user, this session) reshape the design:**

1. **Discipline IS a type** (not a property). The hierarchy is **Petrov-Galerkin-as-base, Galerkin-as-
   specialization** (Liskov-correct: Galerkin *is-a* PG with `test is trial`, strengthening promises).

2. **The physics keystone — homogenization & condensation are Petrov-Galerkin.** The "Galerkin in
   L²(φV)" reading folds the *solution* into the *metric*; that is a forward-flux, reaction-rate-only
   convenience. The eigenvalue-consistent homogenization reactor physics actually requires preserves the
   **bilinear** functional `⟨φ*, Σφ⟩` (perturbation theory: k is stationary w.r.t. the adjoint-weighted
   residual), giving `Σ_R = ∫_R φ*Σφ / ∫_R φ*φ` with **test = φ*·1_R ≠ trial = φ·1_R ⟹ `Π* ≠ R`** —
   irreducibly Petrov-Galerkin. The forward case (`φ*=φ`) is its Galerkin degenerate (= my landed
   `homogenize`). **Consequence:** the Measure carries the **axis** (angular/spatial/energy) + the
   **fixed L² metric**, and **NEVER the discipline**; the discipline + the solution-weighting (φ/φ*)
   live in the **test side = the Frame type**. Keeping the flux OUT of the measure and IN the test side
   is exactly what makes the eigenvalue case expressible.

## The architecture — the Frame discipline-type hierarchy

```
FrameBase  (abstract)            trial basis + measure + table + reconstruction + conjugate
│                                 — the discipline-FREE mechanics (reconstruction & table do not
│                                   depend on the test side)
└─ PetrovGalerkinFrame(FrameBase)   explicit TEST side (a test Basis, or a solution-weighting on the
   │                                trial basis); analysis M = test*·W; project = coefficient
   │                                extraction; Π* ≠ R (general)
   ├─ GalerkinFrame(PetrovGalerkinFrame)   test IS trial; STRENGTHENS: Π* = R (up to tightness),
   │                                       symmetric Gram, diagonal ⟹ project overrides to the
   │                                       reciprocal fast path. Assertable specialisation.
   │                                       ← the angular SH projection (the ONLY pure-Galerkin consumer)
   └─ LeastSquaresFrame(PetrovGalerkinFrame)   test = A·trial; DENSE SPD Gram ⟹ a real solve
                                              (CAP_SOLVE coefficient space). SEAM — zero consumers; the
                                              expressivity boundary the `(basis,measure)` pair can't fake.
```

**The unified API (every consumer collapses to ≤3 lines):**
- `frame.analyze(f) -> coeffs` — M (the un-normalized moment; scattering's projection). EXISTS.
- `frame.reconstruct(c) -> values` — R (dual synthesis). EXISTS.
- `frame.conjugate(A) -> LinearOperator` — **R∘A∘M as a typed `OperatorProduct`, the production path**
  (NEW; replaces scattering's hand-chain). `frame.reconstruct_after(A) = R∘A` is the nameable
  sub-operator the windowed arm re-enters on.
- `frame.project(f) -> coeffs` — the coefficient extraction `G⁻¹ M f` (NEW; the homogenize/condense
  verb). For every real consumer the Gram is **diagonal** (disjoint indicators / SH-orthogonal) ⟹
  `apply_inverse_metric` (a reciprocal); the dense solve is the LeastSquares seam only.
- **The test side carries the solution-weighting**: `PetrovGalerkinFrame(trial, measure, test=<weighted>)`.
  Forward homogenization: `test = φ·1_R`. Eigenvalue homogenization: `test = φ*·1_R`. The MEASURE stays
  the pure geometric `volume_measure`; the flux is the TEST, never folded into the metric.

**Typed-carrier-aware** (the original "operators speak typed carriers" goal): the faces consume/emit
typed `Field`s (`analyze(AngularFlux) -> HarmonicMomentFlux`), not bare `.values`. Runtime
`FunctionSpace` domain/codomain already type the composition end-to-end — the two-param `LinearOperator[Din,Cout]`
type split is **orthogonal** (a separate #226 concern), so this carve does NOT require it.

**Where the Gram-inverse lives (decided):** keep the per-basis split (SH bakes the analytic `2ℓ+1` into
`reconstruct`, measure-free; Indicator defers `1/m_R` to the space metric). The weighted projection's
`G_w` is a *different* object (the state-weighted frame operator, owned by `project`); no single-source
conflict, and the SH `reconstruct` path stays UNTOUCHED (0-ULP canary safe).

## Consumers after the carve (the acceptance test)

```python
# angular projection (M)            — GalerkinFrame(SH, quadrature)
moments = quad.angular_frame(L).analyze(angular_flux)                 # typed in/out

# scattering  S = (1/W)·R∘Λ∘M       — GalerkinFrame.conjugate, PRODUCTION (cast + hand-chain retired)
kernel  = frame.conjugate(legendre_scattering) / W                   # the matvec IS this operator

# forward homogenization            — PetrovGalerkinFrame, test = flux-weighted indicators
coarse_xs = PetrovGalerkinFrame(coarse.indicator_basis(),
                                fine.volume_measure,
                                test=fine.flux_weighting()).project(fine_xs_field)

# eigenvalue-consistent homog.       — same, test = adjoint·forward weighting  (needs φ*, Phase 6)
coarse_xs = PetrovGalerkinFrame(coarse.indicator_basis(), fine.volume_measure,
                                test=fine.adjoint_forward_weighting()).project(fine_xs_field)

# energy condensation                — same shape, energy axis
coarse_xs = PetrovGalerkinFrame(EnergyGrid(coarse_groups).indicator_basis(),
                                energy_measure, test=spectrum_weighting()).project(fine_xs_field)
```

## Phases (each: green + CLI-pyright ≤ baseline + ff-merge; surgical, main-agent-direct)

- **P1 ✓ LANDED `00f9b76` — the discipline-type hierarchy + projection.py reconciliation.** Mint `FrameBase`,
  `PetrovGalerkinFrame`, `GalerkinFrame`; make the current `Frame(basis, measure)` callers construct a
  `GalerkinFrame` (test=trial — behaviour-identical, the realised discipline today). Wire the faces to
  the `AnalysisOperator`/`ReconstructionOperator` operator roles (the open frame_basis_carve §0 follow-up).
  **Retire `GalerkinProjection`/`PetrovGalerkinProjection`** (discipline now a Frame type, #268) + rewrite
  `galerkin_projection.rst` to the PG-base architecture. Gate: 0-ULP scattering canary; all numerics+sn green.
- **P2 ✓ LANDED `cc6d022` — `conjugate` + scattering composed-operator-as-production.** Build `frame.conjugate(A)` /
  `reconstruct_after(A)` (typed `OperatorProduct`); give `LegendreMomentScattering` real
  `domain=codomain=basis_space` (additive leaf ⟹ **retire the `cast`**); make scattering's `apply` arms
  call the conjugation (retire `build_aniso_source`/`_aniso_source_from_moment_values` hand-chains; the
  windowed arm uses `reconstruct_after`). The 0-ULP `test_scattering_kernel_crosscheck` flips from
  twin-guard to **definition**. Gate: 0-ULP canary; sn/operators+solve+sweep green.
- **P3 ✓ LANDED `79e2ac8` — `project` verb + re-frame homogenization as PG (honest test side).**
  `FrameBase.project(f)=G⁻¹Mf` + cached `frame.gram` (diagonal cross-Gram `G_R=(M∘R)(𝟙)=∫_R w dV` on the
  coefficient space; `apply_inverse_metric` Moore–Penrose zeroes empty regions). Minted
  `WeightedIndicatorBasis` (the PG test side: wraps trial indicator + nodal weight, folds it into
  `analyze` via a **leading-aligned broadcast** → source-group weighting for matrices; `evaluate`→plain
  membership; synthesis side RAISES = test-only, never half-consistent).
  `MaterialXSField.project_through(sigma_frame, emission_frame)` projects the XS field as ONE object —
  the field owns the channel→weighting taxonomy, routing Σ (reaction rate) → `sigma_frame`, χ (emission
  rate, production-weighted `p=Σ_g νΣ_f φ_g`) → `emission_frame`. **TWO conserved functionals ⟹ TWO PG
  frames** (the user-chosen design; flux-vs-emission); inline `_gather_*`/`_collapse_*` retired.
  `Solution.homogenize` re-expressed + made **n-D** (1-D guard dropped): `Mesh1D`/`Mesh2D` each gained
  polymorphic `indicator_basis()` + `with_distinct_cell_ids()` (no Mesh1D/Mesh2D reconstruction
  dispatch — the user-preferred clean primitive over isinstance). GATES: 1-D **bit-identical**
  (homogenization gate 11/11 unchanged numbers); 2-D rate preservation **2.8e-17** through a real
  `solve_sn`; **Mode-11 sentinel re-pointed** to `WeightedIndicatorBasis.analyze`+`FrameBase.project`
  (teeth CONFIRMED — the old Galerkin metric-fold gives identical numbers yet never constructs the
  weighted test ⟹ value gates blind, only the sentinel reds); new project-verb cross-Gram + Π*≠R
  discriminator + `test_weighted_indicator_basis.py` intrinsic tests; broad numerics+sn 7-and-only-7
  baseline +51 passes; pyright +0; Sphinx -W clean (`discrete_ordinates.rst` homogenization section
  rewritten Galerkin→PG by archivist; retired wrong `sn-homogenization-galerkin-equals-petrov` label,
  preserved `sn-homogenization-rate-preservation`).
- **P4 ✓ LANDED `74378d5`+`ed1e14d`+`83aaa8a` — carrier-grid completion + HarmonicFrame typed seam.**
  Commit 1 (`74378d5`, P4a/b/c): renamed `HarmonicMomentField → HarmonicMomentFlux` (module too) +
  retired the `orpheus.sn` shim (0 importers); minted **`HarmonicMomentSourceSink(MomentField)`** (bare
  MomentField, no FluxRole — `source+source` CLOSED, `SCALAR_RATE_UNITS`); **lifted** the shared
  moment-space machinery (`L`/`spatial_moments`, `_phase_space_shape`, `from_mesh_and_L`/`zeros_*`/
  `from_ndarray`, `_check_partner`) onto `MomentField` (clean-before-extending, the 2nd-leaf trigger);
  built **`HarmonicFrame(GalerkinFrame)`** in `orpheus/transport/frames/` — the typed seam (the casting
  MUST be transport-side: carriers share `Field` in numerics but castability=`mesh`+`from_mesh` lives in
  transport `BulkField`, numerics can't import carriers). Role-polymorphic `analyse`/`reconstruct`
  (singledispatchmethod): AngularFlux↔HarmonicMomentFlux, AngularSourceSink↔HarmonicMomentSourceSink.
  Generic numerics faces UNTOUCHED. Commit 2 (`ed1e14d`, option-2, USER-STEERED "explicit>implicit,
  principled>bit-identical"): `scattering.frame`→HarmonicFrame; `Λ.apply(HarmonicMomentFlux)` is now the
  explicit **role-changing edge** → returns `HarmonicMomentSourceSink` (the ndarray endomorphic arm — the
  `R∘Λ∘M` kernel `OperatorProduct` view — unchanged; both route the single
  `apply_legendre_scattering_moments` kernel); the windowed in-scatter arm = `frame.reconstruct(Λ.apply(φ))/W`
  (explicit typed Λ-then-R, replacing fused `reconstruct_after(Λ)`) — gives the new leaf +
  `HarmonicFrame.reconstruct` their FIRST production consumer. The hot `kernel = conjugate(Λ)` (full-angular
  `build_aniso_source`) STAYS the fused composed op (P2's 0-ULP canary); `_aniso_source_from_moment_values`
  (= `reconstruct_after(Λ)`) UNTOUCHED as the canary oracle. Commit 3 (`83aaa8a`, P4f docs, Archivist):
  `operator_algebra.rst` carrier-grid diagram + HarmonicFrame layering rationale + explicit-vs-fused note.
  GATES: 0-ULP crosscheck green; Phase-5a guard (windowed==AngularFlux arm) bit-identical; Λ test rewired
  to source role; **Mode-11 sentinel** (windowed arm constructs HarmonicMomentSourceSink → typed edge
  executes, no bypass); 178 scattering-unit + 782 solve/sweep/eigenvalue (incl 2G-het keff) + new-leaf
  intrinsic + HarmonicFrame tests, 7-and-only-7 baseline; Sphinx -W clean; pyright deferred (#226).
  **Two-param operator split NOT done (orthogonal, not needed) — confirmed.**
- **P4.5 — the two-param `LinearOperator[Domain, Codomain]` split (the keystone).** FULL SPEC:
  `.claude/plans/archive/p4_5_two_param_operator_split.md` (cold-pickup-ready). Generalize the operator algebra from
  single-param endomorphic `LinearOperator[V]` (`apply(x:V)->V`) to two-param `LinearOperator[Din, Cout]`
  (`apply(x:Din)->Cout`), so the non-endomorphic transport operators (S/F/L mapping **Flux → SourceSink**) are
  honestly typed. This (1) RETIRES the `@overload` "not-an-endomorphism" confession stacks (`scattering.py:1388`,
  `fission.py:466` — Python's weak spot); (2) carrier-types the kernel end-to-end (`kernel.apply(AngularFlux)->
  AngularSourceSink` directly), **dissolving the P4 option-2 fused/explicit asymmetry** (both arms uniformly
  typed); (3) unifies S/F/L under one `Operator[Flux,SourceSink]` contract (the affine-map linear-part structure;
  fission = rank-1 degenerate of scattering's `R∘Λ∘M` frame). Gating decision RESOLVED (user 2026-06-24): **raise
  the floor to Py≥3.13** (`requires-python>=3.13` + pin pyright `pythonVersion=3.13`; native PEP-696, NO
  `typing_extensions`). Sub-phases a(numerics split foundation)/b(S/F/L honest
  types + retire confessions)/c(carrier-typed kernel)/d(projection ABCs)/e(docs). pyright net DOWN is the
  DELIVERABLE (this IS the deferred #226 generic-Protocol cleanup). Runtime-zero (typing + carrier wraps);
  0-ULP canary the load-bearing gate. Proactive test-architect before b/c.
- **P5 — energy condensation (greenfield).** Mint `EnergyGrid` + its `indicator_basis(coarse_groups)`
  (same `IndicatorBasis` class, energy axis) + the energy base-measure (bin-width/lethargy) +
  `EnergyGroupSpace`; `Solution.condense(grid) -> dict[int, Mixture]` as a `PetrovGalerkinFrame`
  consumer (mesh-DECOUPLED — the asymmetry law). Gate: spectrum-weighted condensation rate-preservation
  + a within-group-varying-spectrum discriminator (the energy analog of P3's).
- **P6 — eigenvalue-consistent (adjoint-weighted) homogenization — ✅ DONE (2026-07-26, B0–B4
  on `feature/sn-adjoint-weighted-homogenization`; merge gated on the user checkpoint).**
  The ratified `adjoint:` parameter LANDED on `Solution.homogenize` AND `Solution.condense`
  (`None` ⇒ forward, bit-identical). **Derivation-first (user-ruled)**: the algebra of record
  `orpheus/derivations/common/homogenization.py` (T0 keystone · T1 vector pair weight φ\*⊙φ —
  the '(φ→φ\*)' shorthand was a TRAP, three doc instances corrected · **T1b the exact ANGULAR
  Σt pairing ρ=Σ w ψ\*ψ, user-ruled option 2** · T2 matrix per-pair sink×source · T3 fission
  mixed-fold + canonical convex χ · **T4 worth-exact ⟺ balance BROKEN** (never assert_balanced
  on adjoint-collapsed Mixtures) · T5 forward=O(ε) · **T6 the B&G Ch. 6 energy convention**
  (plain-flux carrier + flux-weighted-average adjoint carrier Ψ†; T6a exact-k at true spectra;
  T6r row-scaling freedom; T6b O(ε²)-vs-O(ε)) — all EXACT, pinned, and PROOF-WELDED to Sum-form
  display equations the doc generator renders (`generate_homogenization_collapse_rules`, the
  DERIVATION-FRAGMENT-GENERATOR's first instance, user-ruled). Five-morphism
  `MaterialXSField.project_through_bilinear`; `Mixture.condense(adjoint_spectrum=)`;
  `SNMesh.is_same_phase_space` (the two-entry pairing tier). Gates: §4.0 0-ULP pins · C1
  full-taxonomy hand rules + T4 live pins (both axes) · **C2 = the SAME-MESH XS-replacement
  CONTRAST ladder (the coarse-resolve comparative FAILED measurably — spec §4 addendum g);
  measured fwd 2.05/2.01 vs adj 6.08/9.24** · C3/C5 Mode-11 weight capture · C4 B&G hand rules
  (addendum h: the sink axis GAINS the adjoint carrier — 'marginalize frozen' was pre-B&G) ·
  Cχ + T6a-live. All mutation teeth probed RED (addendum i; the φ\*≡φ substitution is C1+C3's
  catch, not C2's — recorded). Reviews: enforcer 0 MUST-FIX (findings applied); qa SOUND
  (the tautological-guard catch → vv-principles Mode-8 extension). B&G acquired+OCR'd
  mid-campaign (637 pp sidecar); literature memo Source E carries the carrier verdict.
  Docs: frame.rst landed-state narrative + verification/sn.rst P6 rows (archivist);
  ch15 API block flipped. Gates at close: Sphinx `-E -W` 0 · audit 0 (orphans 0/315,
  ERR 69/69) · pyright floor 1 (#288).
- **P7 — docs + theory — ✅ DONE (audited 2026-07-26; delivered incrementally by the P3/P5/P6
  doc-passes, no dedicated pass needed).** The full charter is present at the articulation
  standard in `docs/theory/foundations/frame.rst` (the self-identified P7 capstone): the
  discipline hierarchy with BOTH rejected alternatives + the Liskov argument + the LS seam;
  `Π*=R` vs `≠R`; measure-is-metric-not-discipline (+ the energy twin); the composed-operator
  verbs (`conjugate`/`reconstruct_after`/`project`) with the Funk-Hecke `S=R∘Λ∘M` story and the
  eigenbasis-ownership ruling; the P6 adjoint-weighted taxonomy slice (T0 keystone → the five
  GENERATED collapse rules → T4 → the B&G energy convention); the coherent energy-condensation
  treatment + the asymmetry law. Cross-cutting rulings homed in `operator_algebra.rst` (double
  category, affine collision split, carrier census). Audit: `-E -W` clean, all cross-refs
  resolve, no load-bearing ruling is plan/issue/docstring-only. Condensation deliberately a
  frame.rst SECTION, not a standalone page (Cardinal Rule 2).

## Critical files
- **Hierarchy/API:** `orpheus/numerics/frame.py` (the type split + `conjugate`/`project`),
  `orpheus/numerics/projection.py` (retire Galerkin/PG ABCs; keep Analysis/Reconstruction roles),
  `orpheus/numerics/basis/` (test-side hooks), `orpheus/numerics/space.py` (EnergyGroupSpace, the
  CAP_SOLVE coeff space seam).
- **Consumers:** `orpheus/sn/scattering.py` (conjugation production + Λ spaces + cast retirement),
  `orpheus/sn/solution.py` (`homogenize` re-frame + `condense`), `orpheus/transport/mesh/material_xs_field.py`
  (`project_through`), the `EnergyGrid` (new), the adjoint-flux carrier (new, P6).
- **Carriers:** `orpheus/transport/fields/harmonic_moment_field.py` (rename), `source_sinks/` (new leaf).

## Verification (per phase + end-to-end)
- `.venv/bin/python -O -m pytest tests/numerics tests/sn/operators tests/sn/solve tests/sn/sweep tests/sn/eigenvalue tests/sn/verification -q -rfE --timeout=120 -p no:cacheprovider`
  — the **0-ULP `test_scattering_kernel_crosscheck`** is the load-bearing canary (any ULP shift in P1/P2 = STOP);
  the homogenization rate gate + mutation probes UNCHANGED through P3; the **7-and-only-7 pre-existing reds** baseline.
- `npx --no-install pyright --outputjson orpheus/` — net reds **down** (cast + ABC retirements); CLI is the oracle.
- Each discipline TYPE ships its intrinsic-property test (the project standard): `GalerkinFrame` asserts
  `Π*=R`; `PetrovGalerkinFrame` asserts the cross-Gram extraction + `Π*≠R`; the P6 adjoint-vs-forward
  discriminator proves the type is load-bearing (different answer).
- Proactive `test-architect` before P2/P3/P6 (operator-algebra carves); `Archivist` for P1/P7 docs;
  `elegance-enforcer` after each consumer migration.

## Scope / discipline
- Surgical, main-agent-direct, user-steered; NO `method-implementer`. Per-phase ff-merge; `main` always green.
- **Build the SEAM, not the speculation:** `LeastSquaresFrame` (dense-Gram solve) is designed-for but NOT
  built (zero consumers); the test side admits `test: Basis | weighting`, not `| Operator`, until an LS
  consumer arrives. P6's adjoint-flux is the real ≥2nd PG instance that justifies the whole type axis.
- **Durable principles banked:** flux-as-test-weighting (NOT flux-as-measure — the solution emits a
  weighting on the test side, the discretization emits the measure); discipline-is-a-Frame-type
  (PG-base, Galerkin-specialisation); measure = axis + metric, never discipline; composed-operator-as-production.
