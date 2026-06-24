# THE projection / reconstruction machinery — the unified Frame campaign

> **APPROVED 2026-06-24** (user). Branch `refactor/operator-inverse-algebra`. Surgical, main-agent-direct,
> user-steered; NO `method-implementer`. Tasks **#46–#52** (P1–P7); P1 (#46) unblocks all. Memory:
> `[[project-frame-projection-machinery]]` (the load-bearing rulings — read it FIRST; it OVERRIDES the
> discipline-as-property conclusion in `[[project-homogenization-condensation]]` / `[[project-frame-basis-carve]]`).
>
> **STATUS (2026-06-24): P1 LANDED `00f9b76` (Closes #268); P2 LANDED `cc6d022` (Refs #226);
> P3 LANDED `79e2ac8` (Refs #226); P4 LANDED `74378d5`+`ed1e14d`+`83aaa8a` (Refs #226).
> NEXT = P4.5 — the two-param `Operator[Domain,Codomain]` split (detailed spec:
> `.claude/plans/p4_5_two_param_operator_split.md`), then P5 (#50).**
> P4.5 is INSERTED before P5/P6 by user directive (2026-06-24): polish the operator algebra so the later
> phases don't build on a type that must further change. It is the keystone the post-P4 review surfaced —
> the uniform `Operator[Flux,SourceSink]` contract — and the genuine fix for the #226 generic-Protocol
> `@overload`/`None`-space confessions (pyright net DOWN becomes the deliverable, not a deferral).
> The pyright ratchet is a GUIDE, not a gate, DURING the campaign — the faces/Λ `ndarray`-vs-`Vector`
> generic-Protocol gap (the cast retirement exposed it; pre-existing) is #226, deferred to ONE
> post-campaign pyright pass (user, 2026-06-24).
> Reconcile against git before resuming.

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
  `.claude/plans/p4_5_two_param_operator_split.md` (cold-pickup-ready). Generalize the operator algebra from
  single-param endomorphic `LinearOperator[V]` (`apply(x:V)->V`) to two-param `LinearOperator[Din, Cout]`
  (`apply(x:Din)->Cout`), so the non-endomorphic transport operators (S/F/L mapping **Flux → SourceSink**) are
  honestly typed. This (1) RETIRES the `@overload` "not-an-endomorphism" confession stacks (`scattering.py:1388`,
  `fission.py:466` — Python's weak spot); (2) carrier-types the kernel end-to-end (`kernel.apply(AngularFlux)->
  AngularSourceSink` directly), **dissolving the P4 option-2 fused/explicit asymmetry** (both arms uniformly
  typed); (3) unifies S/F/L under one `Operator[Flux,SourceSink]` contract (the affine-map linear-part structure;
  fission = rank-1 degenerate of scattering's `R∘Λ∘M` frame). Gating decision FIRST: PEP-696 defaults need Py≥3.13
  or `typing_extensions` (floor is 3.11) — AskUserQuestion. Sub-phases a(numerics split foundation)/b(S/F/L honest
  types + retire confessions)/c(carrier-typed kernel)/d(projection ABCs)/e(docs). pyright net DOWN is the
  DELIVERABLE (this IS the deferred #226 generic-Protocol cleanup). Runtime-zero (typing + carrier wraps);
  0-ULP canary the load-bearing gate. Proactive test-architect before b/c.
- **P5 — energy condensation (greenfield).** Mint `EnergyGrid` + its `indicator_basis(coarse_groups)`
  (same `IndicatorBasis` class, energy axis) + the energy base-measure (bin-width/lethargy) +
  `EnergyGroupSpace`; `Solution.condense(grid) -> dict[int, Mixture]` as a `PetrovGalerkinFrame`
  consumer (mesh-DECOUPLED — the asymmetry law). Gate: spectrum-weighted condensation rate-preservation
  + a within-group-varying-spectrum discriminator (the energy analog of P3's).
- **P6 — eigenvalue-consistent (adjoint-weighted) homogenization — the headline PG consumer.**
  DEPENDS on a wired **adjoint flux φ\*** (operator `.H` makes it constructible — `A_loss.H.solve(R)` —
  but it is NOT yet a `Solution`-level capability; wire `Solution.adjoint()` / an `AdjointFlux` carrier
  first). Then `test = φ*·1_R` (or the bilinear `φ*φ`), giving the eigenvalue-preserving Σ_R. The genuine
  `Π* ≠ R` PG instance; ships the GalerkinFrame-vs-PetrovGalerkinFrame discriminating test (a case where
  the adjoint-weighted answer ≠ the forward-weighted answer ⟹ the discipline type is load-bearing).
  Gate: k-eff preservation under coarse re-solve (≥2G heterogeneous) vs the fine reference.
- **P7 — docs + theory.** The PG-base projection theory page (the discipline hierarchy, `Π*=R` vs `≠R`,
  the measure-is-metric-not-discipline ruling, the conjugation primitive), the adjoint-weighted
  homogenization derivation, the condensation page. Sphinx `-W` clean throughout.

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
