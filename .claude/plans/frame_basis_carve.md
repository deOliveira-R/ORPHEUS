# Build the `Frame` + `Basis` ABC; unify projection/reconstruction as a discrete frame

## §0. STATUS & cold-pickup — read FIRST (Phases A–D + measure-completeness COMPLETE, 2026-06-23)

**Branch** `refactor/operator-inverse-algebra` (NOT yet merged to main). **Reconcile against git
first** (`git log --oneline -8`). Surgical, main-agent-direct, NO `method-implementer`. Canonical
tests `.venv/bin/python -O -m pytest`; CLI `npx --no-install pyright --outputjson orpheus/…` is the
type oracle (the streamed `<new-diagnostics>` "could not resolve …" is #226 LSP lag — verify
with the CLI). NO `# type: ignore`.

**Landed (committed): Phases A–D ALL DONE.**
- `c65154a` **Phase A** — `Basis` ABC (`orpheus/numerics/basis/base.py`); `SphericalHarmonicBasis(Basis)`.
- `4995d6a` **Phase B** — `orpheus/numerics/frame.py` + the basis weighted contractions. Additive.
- `92e8842` **Phase C1** — scattering kernel composes `frame.analysis`/`frame.reconstruction`; 2/3
  casts retired (Λ keeps 1 = #226). Bit-identical (0-ULP canary).
- `f926b0e` **Phase C2** — in-sweep moment carrier `moment_projection`→`moment_frame: Frame`;
  `_MomentWindowedResolvent` holds `scattering_op.frame`; added **`Quadrature.angular_frame(L)`**
  (single home of the S² embedding + SH-frame binding). Consumers read order off
  `frame.table.shape`/`bulk_values.shape` (NOT `.basis.L` — sweep stays basis-agnostic + type-clean).
- `501d443` **Phase C3** — DELETE `MomentProjection`+`HarmonicMomentReconstruction`+`two_l_plus_one`;
  rename `ProjectionOperator`→`AnalysisOperator` (keep `ReconstructionOperator`/`GalerkinProjection`/
  `PetrovGalerkinProjection`); all consumers + 7 doc pages (archivist) + tests rewired;
  `test_projection_operators.py` retired (laws moved to `test_spherical_harmonic_space.py`).
- `1271490` **Phase D** — the reconstruction face earns `R.H` (adjoint-for-free, symmetric with
  analysis). New `Basis.reconstruct_transpose` (the **measure-free** representation transpose
  `R^⊤ = d_k·S_0^⊤` — NO `w_n`, asymmetric with `analyze_transpose` whose forward bakes weights in);
  `SphericalHarmonicBasis` impl `(2ℓ+1)·Σ_n Y_ℓ^m v_n`; `_FrameReconstruction.apply_transpose` +
  `CAP_APPLY_TRANSPOSE` → `R.H = (2ℓ+1)²/4π · Σ_n w_n Y_ℓ^m v_n` falls out of `_AdjointOperator`.
  3 new tests (symmetric with the analysis-face trio): bit-identical per-term-fold + unit-vector pin,
  both adjoint identities by the DEFINING inner-product law, R.H vs independent closed-form. qa SUPPORTED
  (adjoint re-derived 3 ways, all gates mutation-confirmed to bite, Mode-11 cleared); 0-ULP canary holds;
  prod pyright clean; Sphinx -W clean. (`743e954` = matrix regen.) **Did NOT build `R.H.apply_transpose`**
  (the `_AdjointOperator` stub still raises — zero consumers, per plan).

**Landed (committed): MEASURE-COMPLETENESS** (after Phase C; the user reviewed the Frame for naming +
information-completeness, this is the agreed 4-phase result — plan `.claude/plans/reactive-moseying-cake.md`,
now COMPLETE). `df700aa..d35003e`, pushed.
- `b2ed887` **SH naming** — prose `order`→`degree` (Y_ℓ^m: ℓ=degree, m=order; ~8 docstrings were a
  latent convention bug). **Param `L` kept** (highest-signal in both `ℓ≤L` and `P_L`). `angular_frame`
  docstring tripwire: it names the permanent **axis**, so a 2nd angular basis is an `angular_frame(basis=…)`
  *signature* change, NOT a rename.
- `163b236` **typed `invariance_group`** — `DiscreteMeasure.invariance_group: SubgroupOfO3 | None`
  (was `str|None`, set inconsistently); angular factories pass the typed singleton (`OctahedralOh`/`SO2`).
  TYPE_CHECKING forward-ref (symmetry imports FROM measure → runtime import would cycle).
- `fee5ff8` **the measure OWNS its space (symmetric domain/codomain)** — field `space`(str)→**`support`**
  (the continuous σ-algebra tag: `"S^2"`/`"[-1,1]"`/`"spatial_R1"`/`"cells"`). New **`measure.space` →
  `FunctionSpace`** (the induced discrete-L² DOMAIN, `(N,)`+weights metric), mirroring `basis.space`
  (codomain). `Frame.measure_space` is now just `return self.measure.space` (ad-hoc fabrication gone).
  Bit-identical (metric VALUE unchanged; only its source/name). ~37 `space=`→`support=` sites migrated.
- `de11437` **`measure.phase` derived** — `Literal["angular","spatial","energy"]`, a phase-space FACTOR
  (position×direction×energy), NOT a literal: **angular** iff `invariance_group is not None` (the
  O(3)-subgroup symmetry IS the signature — Erlangen); **spatial** iff a `"spatial…"`/`"cells"` support;
  else **raises** (energy + bare generic rules undetermined — a slab μ∈[-1,1] is geometrically
  indistinguishable from a spatial interval, so the physical identity must be supplied). New
  `tests/numerics/test_measure_phase.py` (11 foundation tests) pins the derivation + the per-category seam.
- `d35003e` — pre-existing `:paramref:` Sphinx error in `mesh.py` (surfaced when the measure work
  invalidated the doc cache) → plain ``origin`` literal. Sphinx `-W` clean.
- **DURABLE RULING (user):** the three factors (angular/spatial/energy) are GENUINELY DIFFERENT objects
  (compact O(3)-sphere vs unbounded energy half-line vs bounded mesh domain). The asymmetry IS the signal
  to type them **per-category, minted as each is first Frame-bound** — NOT a premature uniform abstraction.
  Angular is the worked instance (symmetry-derived); spatial/energy stay support-tag-recognised until bound.
- **DEFERRED** (per type-vs-property / #263 / "mint per-category as bound"): typed per-category support
  hierarchy (`AngularSupport`/`SpatialSupport`/`EnergySupport`); the energy measure + half-line support;
  `SPACE_*`→`SUPPORT_*` const rename; a `Frame.__post_init__` phase-mismatch guard (needs `Basis.phase` +
  a 2nd axis); the **Galerkin/PG discipline** field (mint at the 2nd genuinely-oblique frame — currently
  `CAP_SOLVE`-decidable — mint at Phase E+F's first oblique frame). These do NOT block Phase E+F.

**OPEN follow-up flagged to user (NOT yet decided):** the 4 `projection.py` ABCs are now
implementation-free (frame faces subclass `LinearOperatorMixin` directly). Either (a) make the frame
faces subclass the role ABCs, or (b) retire `projection.py` entirely. Also `Quadrature.angular_frame`
placement is open for review.

**As-built interface (refines the Design below — this is the source of truth):**
- **`Basis` ABC** (`basis/base.py`): `evaluate(points)→table` (the ONLY points-taking method);
  then **table-based contractions** (the Frame caches the table once, delegates — L16 perf guard):
  `synthesize(coefficients, table)` (naked `S₀`), `analyze(values, table, weights)` (M),
  `analyze_transpose(coefficients, table, weights)` (Mᵀ), `reconstruct(coefficients, table)` (R),
  `reconstruct_transpose(values, table)` (Rᵀ = `d_k·S₀ᵀ`, **measure-free** — Phase D);
  plus `mass_matrix(measure)` (diagnostic, measure-based) and the property **`space`** (the
  coefficient `FunctionSpace`; NOT `coefficient_space`/`basis_space` — those are redundant on the
  basis). Abstract-method args are positional-only so concrete bases keep domain names.
- **`SphericalHarmonicBasis`** implements all (the M/R einsums moved here verbatim, **bit-identical**;
  `reconstruct` reads `(2ℓ+1)` live from `addition_theorem_factor`; `space` = lazy
  `SphericalHarmonicSpace.from_L(L)`).
- **`Frame(basis, measure)`** (`frame.py`): cached `table`, `measure_space` (domain; weights metric),
  `basis_space` (codomain = `basis.space`), and the faces `frame.analysis` (`_FrameAnalysis`) /
  `frame.reconstruction` (`_FrameReconstruction`) — private frozen dataclasses subclassing
  `LinearOperatorMixin` (unparametrised), `capabilities` a plain unannotated class attr (the
  `block_role` override pattern). BOTH faces now carry `apply`+`apply_transpose`+
  `CAP_APPLY_TRANSPOSE` (→ `.H` free) — analysis since Phase B, reconstruction since Phase D.
- **Naming decisions (durable, in memory):** faces `analysis`/`reconstruction` (NOT
  `analysis`/`synthesis` — "synthesis"=`T*`=`S₀` is the naked basis primitive, a DIFFERENT object;
  see `[[feedback-high-signal-names]]`). The shared `S₀` is conceptual only — production keeps
  FUSED einsums (factored `w⊙S₀` drifts 3.6e-15, breaks 0-ULP).

**Verified at Phase B:** faces `np.array_equal` to `MomentProjection`/`HarmonicMomentReconstruction`;
`frame.analysis.H == M.H`; `reconstruction @ analysis` composes through `OperatorProduct` with real
spaces (no cast). 669 numerics + 129 sn eigenvalue/scattering green incl. the **0-ULP
`test_scattering_kernel_crosscheck`**. **Pre-existing unrelated failures (NOT ours, verified via
stash):** 8 in `tests/sn/operators/` (2D-mesh `mu_y` cosines; sphere streaming snapshots).

**NEXT = Phase E+F (#31)** — the **solution-as-measure** seam: `ScalarFlux.flux_volume_measure()`
(or equivalent) emits the flux·volume PG weighting as a `DiscreteMeasure` the SOLUTION produces
(keystone ruling below: the PG weightings are `Measure`s the solution emits, unifying with the
discretisation's geometric quadratures — `mesh.volume_measure` already emits one, flux·volume·Σ rates
already route through it). This is where the first genuinely-oblique (Petrov-Galerkin) measure lands —
and thus where the deferred **Galerkin/PG discipline** field is finally minted (the `support`/`phase`
per-category machinery from measure-completeness is the substrate). Then final docs polish (archivist:
a frame theory page + the reconstruction-adjoint equation labels the Phase-D tests deferred) and **EXIT
to #226** (pyright burn-down). Gate: numerics + sn/operators green; 0-ULP canary; CLI pyright ≤ baseline.

Gate discipline (every phase, held all of C): the **0-ULP `test_scattering_kernel_crosscheck`
canary** + the windowed-moments oracle (`test_2d_anisotropic_windowing`) + the structural anchors
(`test_spherical_harmonic_space`) + net-zero new CLI-pyright reds. Pre-existing unrelated reds (NOT
ours, stash-verified): 7 `tests/sn/operators/` (2D `mu_y`/sphere snapshots).

Tasks: #27/#28/#29/#30 DONE → #31 (E+F). Memory: `[[project-frame-basis-carve]]`,
`[[feedback-high-signal-names]]`.

---

## Context

This started as the "typed-carrier grid" carve and the narrower "unify M+R" refactor, but a
deep architecture discussion with the user reframed it. The findings:

1. **Partition (RESOLVED, separate work):** the `numerics(L1) | transport(L2) | method(L3)`
   three-folder split is **correct**, keyed by *coupling* (the decidable test: *what concrete
   type does a primitive NAME?*), not by abstraction altitude. `Vector`/`Generic[V]` is the
   load-bearing genericity mechanism, vindicated by the multi-physics thought experiment
   (one operator algebra must serve transport + fuel + thermal-hydraulics). Two consequences,
   both **future work, NOT this carve**: (a) the transport-state-space primitives currently
   mis-filed in numerics (`BlockRole`, `FullFieldSpace`, `TraceSpace`'s partial-current
   metric, the per-leaf unit table) move to transport; (b) `Collision`/`Fission`/`Scattering`
   are transport physics shared by all methods — only `StreamingOperator` (the method's
   Feynman-Kac view of the path integral) + its `LossRepresentation` are method-specific.

2. **The projection machinery (THIS carve):** projection ↔ reconstruction is a **discrete
   frame** (harmonic-analysis frame theory). A frame = a **`(basis, measure)` binding** that
   emits the operational `(R,Π)` pair — **analysis** `M=T` (nodes → coefficients) and
   **reconstruction** `R` (coefficients → nodes, the dual-frame synthesis); the naked
   synthesis `T*=S₀` is the shared basis primitive (`basis.synthesize`). The **basis fixes the codomain**
   (coefficient space + Gram — disciplineless trial side); the **measure fixes the domain**
   (nodal space + weights) **and the Galerkin-vs-Petrov-Galerkin discipline** (the test side).
   So `coefficient_space` is NOT a third parameter — it's derived from the basis (the
   `two_l_plus_one` smell one level up).

   **The Frame-vs-property line is CHOICE-DEPENDENCE, not iso-vs-non-iso** (GitHub #263's
   uplift criterion, surfaced #257 S9 — this CORRECTS the earlier "iso → property, non-iso →
   frame" framing). A representation earns a **frame + two first-class types iff** a dual
   basis coexists connected by a map that (a) **depends on a CHOICE** — a quadrature /
   node-set / measure, not a fixed formula — (b) is **actually modeled & applied** (carries
   truncation error, has an adjoint, lives in the operator algebra), and (c) whose
   **forbid-mixing gate guards a real bug**. **The integrating quadrature is what breaks the
   canonical isomorphism** — that break justifies the second type + the frame (the connecting
   morphism `M`/`R`). A **canonical** (fixed, choice-free) map — e.g. the LD modal
   `edge = avg ± slope/2`, where the coefficients ARE the function — is a **property/accessor
   on ONE type**, NOT a frame and NOT a second type (uplifting it = a "naming leaf").
   **Within the choice-dependent (frame) cases, iso-vs-non-iso is only a CAPABILITY:** an
   invertible square-Vandermonde map advertises `.inverse()`/`CAP_SOLVE` (the iso **nodal-DG**
   case); a lossy band-limiting map (`R∘M` = projector ≠ I, `N>(L+1)²`, SH angular) is a
   section/retraction (non-iso). **Rule:** the Frame is the SINGLE mechanism for ALL
   choice-dependent change-of-basis — both iso and non-iso, no twin path (a perf shortcut is
   allowed only as a frame-pinned oracle); a **canonical** map stays a property. This answers
   "use the Frame for both" (yes, for both iso & non-iso *choice-dependent* frames) AND why
   the spatial LD moment is correctly a property today (canonical map, clause-1 fails — and
   confirmed: the spatial Legendre basis is *operator-embedded* in the LD cell-solve, not a
   collapsible twin).

3. **Decision (user, no defer):** build the generic **`Frame`** + the **`Basis` ABC** with
   `SphericalHarmonicBasis` as the first concrete basis; unify the SH projection/reconstruction
   into a frame; give reconstruction its spaces (adjoint-for-free); delete the derived
   `two_l_plus_one` twin. The base nature is understood, so deferral is no longer warranted.

4. **Solution-as-measure (keystone principle):** the Petrov-Galerkin weightings (energy
   *spectrum*, spatial *flux·volume*) ARE measures the **solution** emits — unifying with the
   geometric quadratures the *discretization* emits, under one `Measure` type. This is
   **half-built**: `mesh.volume_measure` already emits a `DiscreteMeasure` and the
   flux·volume·Σ rates already route through it. Recorded as a first-class principle; the
   provider build is scoped below.

**Name:** `Frame` — confirmed by the cross-domain-attacker as the precise math object (a
*discrete frame*; the SH case is *4π-tight*). Rejected `Transform`/`Bijection` (imply
invertibility — false) and `GalerkinTransform` (its own "Transform misleads" argument).
Faces: **analysis/projection** and **synthesis/reconstruction** (both vocabularies in docs).

## What this supersedes

- The typed-carrier grid plan (tasks #22–26) and the M+R-only unification plan. The grid's
  friction was a *symptom* of the partition + the missing frame abstraction; it is retired.
- The earlier "rename `MomentProjection`→`HarmonicMomentProjection`" / "keep
  `HarmonicMomentReconstruction`" naming fork: **mooted** — the faces become generic
  frame-backed views (`frame.projection` / `frame.reconstruction`); the SH-specific operator
  class names retire. Harmonic-ness lives in the basis, not the operator names.

## Design

### 1. `Basis` ABC (un-defer) — `orpheus/numerics/basis/`
Promote the informal contract (`basis/__init__.py` l.15-18 deferral) to a real ABC. Abstract
surface (already satisfied by `SphericalHarmonicBasis`):
- `evaluate(directions) -> table` (the `Φ(node, mode)` tabulation; SH: `(N, L+1, 2L+1)`)
- `mass_matrix(measure) -> Gram` (reconcile the `discrete_mass_matrix` docstring-name drift)
- `metric_per_ell` / the **continuous Gram diagonal** `g_C` (coefficient-space metric)
- `addition_theorem_factor` / the **dual-frame reconstruction weight** (`= bound · g_C⁻¹`)
- `synthesize(coefficients, directions) -> nodes` (the naked `S₀`)

`SphericalHarmonicBasis(Basis)` = first concrete. **Justification for un-deferring on one
concrete** (overrides `feedback_unify_after_two_instances`): the `Frame` is the *forcing
consumer* — it binds an abstract basis and applies a non-identity morphism (project/reconstruct)
across the interface; the interface is math-rigid and grand-report-specified (§l.567-573:
`evaluate`/`mass_matrix`/`project`/`reconstruct`). State this in the ABC docstring.

### 2. `Frame(basis, measure)` — the binder — `orpheus/numerics/frame.py` (new, L1)
Frozen dataclass, two fields `basis: Basis`, `measure: DiscreteMeasure`. Everything derived:
- `table` (cached) = `basis.evaluate(measure.nodes)`  — exposes `.Y` / `.L` delegators
  (CONSTRAINT: the in-sweep walk reads `.L`/`.Y` off the analysis face — keep them as plain
  attribute reads, `operator.py:1108`, `loss_representation.py:3428/3431`).
- `measure_space` (= the operator-algebra `domain`) = `FunctionSpace` from the measure (nodes
  + `weights` as the metric).
- `basis_space` (= the operator-algebra `codomain`) = the basis's space + Gram (for SH:
  `SphericalHarmonicSpace` derived from the basis — already true; `from_L` builds the metric
  from `basis.metric_per_ell`). **`basis space`/`measure space` are the conceptual names** (the
  inputs name the spaces they produce); `.domain`/`.codomain` are the operator-algebra aliases.
- **frame operator** `S = T*T` = the metric (`= 4π·I` for the tight SH case) — **derived from
  the basis, never hard-coded** (a future Riesz/PN basis has `M∘R = I`, not `4π·I`).

### 3. The two operator faces (composable `LinearOperator` views) — `analysis` / `reconstruction`
**Naming (high-signal, corrected):** the faces are the operational `(R, Π)` pair —
`frame.analysis` (`M`=`T`) and `frame.reconstruction` (`R`, the dual-frame synthesis). NOT
`projection` (in a `Frame` "projection" is the idempotent **projector** `P = R∘M` — the
codebase already calls `R∘M` the "band-limited projector" — so reusing it for `M` collides;
`analysis`=`T` is collision-free). NOT `synthesis` for the second face — in frame theory
"synthesis" is the **naked** operator `T*=S₀=Σ Y c` (no weights, no dual factor), which is a
**basis primitive** (`basis.synthesize`), one level below the faces; the face is the
*reconstruction* `R=(2ℓ+1)·S₀`. Face↔method names line up: `analysis↔analyze`,
`reconstruction↔reconstruct`. Per the approved Phase-B fork, the **basis owns the weighted
einsums** (`analyze`/`reconstruct`) — the `Frame` is layout-agnostic and delegates.
Operator role types → `AnalysisOperator` (rename `ProjectionOperator`); **keep**
`ReconstructionOperator`. The **Galerkin / Petrov-Galerkin** discipline survives as the
frame's **canonical-dual vs oblique-dual** property (kept); FEM "projection" stays in docstrings.

**The shared `S₀` (the elegance):** three diagonal-weighted variants share ONE named naked
synthesis `S₀` = `basis.synthesize`:
`analysis.apply_transpose = wₙ⊙S₀` · `analysis.H = g_C⊙S₀` · `reconstruction.apply = (2ℓ+1)⊙S₀`.
- `frame.analysis` → `M`=`T` (`AnalysisOperator`): `apply` = `basis.analyze(x, measure)` (SH
  einsum `"n,nlm,n...->lm..."`, **verbatim for 0-ULP**); `apply_transpose(c) = measure.weights
  ⊙ basis.synthesize(c, measure.nodes)` (= `wₙ·S₀`, generic via the primitive); `domain`=measure
  space, `codomain`=basis space, `CAP_APPLY_TRANSPOSE` → `.H` free.
- `frame.reconstruction` → `R` (`ReconstructionOperator`): `apply` = `basis.reconstruct(c,
  measure.nodes)` (SH einsum `"nlm,l,lm...->n..."`), with the `(2ℓ+1)` factor **read from
  `basis.addition_theorem_factor`, NOT stored** — deletes the `two_l_plus_one` field, makes
  `g_C·(2ℓ+1)=4π` true by construction. **Gains** `domain`=basis space, `codomain`=measure
  space, `apply_transpose`, `CAP_APPLY_TRANSPOSE` → `R.H` free via the metric-aware
  `_AdjointOperator` (Phase D; do NOT build `R.H.apply_transpose` — stub raises, zero consumers).
- Faces accessed via `frame.analysis`/`frame.reconstruction`; view classes frame-backed
  (private `_FrameAnalysis`/`_FrameReconstruction`). `MomentProjection` /
  `HarmonicMomentReconstruction` **retire** into the basis methods + these faces (test
  migration, not deletion).

### 4. Solution-as-measure (principle + thin seam, NOT the full PG build)
- **Record** the principle in `measure.py` / the Frame docs: geometric measures come from the
  *discretization*, physical measures (flux·volume, spectrum) from the *solution*; both are
  `DiscreteMeasure`; the frame's analysis weighting is either.
- **Seam (cheap, latent consumers exist):** add `ScalarFlux.flux_volume_measure()` (sibling of
  `integrate_angular`) = `replace(mesh.volume_measure, weights = volumes * flux.values)`. Legal
  L2→L1 (transport field returns a numerics `Measure`). OPTIONAL in this carve — include only
  if it lands clean; otherwise document the seam and defer with the PG frames.
- **Out of scope (future, with energy condensation / spatial homogenization):** the
  `PetrovGalerkinProjection` concrete frames, the spectrum-measure provider, and the coarse
  codomains (`RegionSpace`/`EnergyGroupSpace` — currently docstring-only). The frame is built
  Galerkin (the realized discipline); the oblique-dual PG path is the documented seam.

## Migration (phased, each phase green + ff-merge)

- **A. `Basis` ABC** — formalize the interface, `SphericalHarmonicBasis(Basis)`, reconcile the
  `mass_matrix` name. Gate: basis + projection + SH-space tests green.
- **B. `Frame` + faces** — build `Frame`, the two view faces, derive spaces, derive the dual
  factor, delete `two_l_plus_one`. Behavior-preserving; faces bit-identical to today's M/R.
  Gate: `tests/numerics/test_projection_operators.py`, `test_spherical_harmonic_space.py`
  (Π R=4π·I, Π*=g_C·S₀, R=(2ℓ+1)·S₀) green.
- **C. Consumer migration** — repoint `scattering.py` kernel (`OperatorProduct` over
  `frame.analysis`/`.reconstruction`; **retire the `cast()`s** — spaces are now present, and
  Λ carries `None` spaces so the `OperatorProduct` guard short-circuits), `build_aniso_source`,
  `_aniso_source_from_moment_values`, and the in-sweep `.L`/`.Y` reads (`solver.py`). Gate:
  `tests/sn/operators/test_scattering_kernel_crosscheck.py` **0-ULP** (load-bearing canary),
  sn/operators + sn/solve + sn/sweep green.
- **D. Reconstruction adjoint-for-free** — `R.apply_transpose` + `CAP_APPLY_TRANSPOSE`; update
  `TestCapabilities`. Gate: numerics + sn/operators green; CLI pyright reds ≤ baseline (cast
  retirement should reduce them).
- **E. (optional) `ScalarFlux.flux_volume_measure()`** seam + record the principle.
- **F. Docs** — frame theory + (R,Π) vocabulary on the new module; reconcile the pervasive
  stale `HarmonicMomentProjection` name in `docs/theory/galerkin_projection.rst` etc., the
  phantom `from_spherical_harmonic_space`, the `HarmonicMomentReconstruction(Y, two_l_plus_one)`
  ctor signature. Sphinx build clean.

## Verification
- `.venv/bin/python -O -m pytest tests/numerics tests/sn/operators tests/sn/solve tests/sn/sweep tests/sn/eigenvalue tests/sn/verification -q` — Galerkin invariants, **0-ULP scattering
  kernel crosscheck** (the canary: any ULP shift = STOP), keff (≥2G heterogeneous), MMS green.
  Tolerances unchanged (einsums reorganized through the frame, not altered).
- `npx --no-install pyright --outputjson orpheus/` — net reds **down** (cast retirement); no
  new reds. Re-baseline the ratchet.
- `grep` confirms the frame is real: `Frame`, `Basis` ABC, `SphericalHarmonicBasis(Basis)`, no
  `two_l_plus_one` field, faces via `frame.analysis`/`.reconstruction`, spaces `basis_space`/
  `measure_space`.

## Scope boundaries (explicit)
- **IN:** `Basis` ABC + `SphericalHarmonicBasis`; generic `Frame`; the two faces; delete
  `two_l_plus_one`; R adjoint-for-free; kernel/consumer migration + cast retirement; tests +
  docs; record solution-as-measure + (optional) the flux-volume seam.
- **OUT (future, named):** the numerics→transport primitive moves (BlockRole etc.); the
  Collision/Fission/Scattering→transport + Streaming/LossRepresentation→method movement
  (**#261**); the Petrov-Galerkin frames (energy condensation, spatial homogenization) +
  their coarse codomains + the spectrum-measure provider; the §10 PN Riesz-basis second
  instance; the **spatial Frame extraction** — operator-extract the LD Legendre basis from
  the UBLD cell-solve + re-derive the closure on top of a `Frame(LegendreBasis,
  spatial_measure)` (entangled with the per-cell hot path). **Trigger (per #263):** a method
  introducing a **committed point-value COLLOCATION basis** coexisting with the modal moments
  via a node-dependent Vandermonde — i.e. **nodal-DG SN (Gauss-Lobatto) or Lagrange-FEM**,
  which mints `SpatialMomentField` + its collocation dual. **NOT triggered** by NEM/ANM
  "nodal" diffusion or modal-FEM (canonical/modal → property). The faces' existing
  capability-set mechanism is the seam — a future iso frame adds `CAP_SOLVE`/`.inverse()`
  additively; **no iso logic is built this carve** (SH is non-iso → section/retraction only).
- **Surgical discipline:** main-agent-direct, user-steered, NO `method-implementer`;
  `test-architect`/`explorer`/`qa`/`elegance-enforcer` available. CLI pyright is the oracle
  (not the LSP stream). NO `# type: ignore`. Per-phase ff-merge.
