# Build the `Frame` + `Basis` ABC; unify projection/reconstruction as a discrete frame

## §0. STATUS & cold-pickup — read FIRST (updated 2026-06-23, mid-carve)

**Branch** `refactor/operator-inverse-algebra`. **Reconcile against git first** (`git log --oneline -5`). Surgical, main-agent-direct, NO `method-implementer`. Canonical tests
`.venv/bin/python -O -m pytest`; CLI `npx --no-install pyright --outputjson orpheus/…` is the
type oracle (the streamed `<new-diagnostics>` "could not resolve …" is #226 LSP lag — verify
with the CLI). NO `# type: ignore`.

**Landed (committed):**
- `c65154a` **Phase A** — `Basis` ABC (`orpheus/numerics/basis/base.py`); `SphericalHarmonicBasis(Basis)`.
- `4995d6a` **Phase B** — `orpheus/numerics/frame.py` + the basis weighted contractions. Additive;
  M/R untouched.

**As-built interface (refines the Design below — this is the source of truth):**
- **`Basis` ABC** (`basis/base.py`): `evaluate(points)→table` (the ONLY points-taking method);
  then **table-based contractions** (the Frame caches the table once, delegates — L16 perf guard):
  `synthesize(coefficients, table)` (naked `S₀`), `analyze(values, table, weights)` (M),
  `analyze_transpose(coefficients, table, weights)` (Mᵀ), `reconstruct(coefficients, table)` (R);
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
  `block_role` override pattern). analysis carries `apply`+`apply_transpose`+`CAP_APPLY_TRANSPOSE`
  (→ `.H` free); reconstruction `apply`-only (Phase D adds its transpose).
- **Naming decisions (durable, in memory):** faces `analysis`/`reconstruction` (NOT
  `analysis`/`synthesis` — "synthesis"=`T*`=`S₀` is the naked basis primitive, a DIFFERENT object;
  see `[[feedback-high-signal-names]]`). The shared `S₀` is conceptual only — production keeps
  FUSED einsums (factored `w⊙S₀` drifts 3.6e-15, breaks 0-ULP).

**Verified at Phase B:** faces `np.array_equal` to `MomentProjection`/`HarmonicMomentReconstruction`;
`frame.analysis.H == M.H`; `reconstruction @ analysis` composes through `OperatorProduct` with real
spaces (no cast). 669 numerics + 129 sn eigenvalue/scattering green incl. the **0-ULP
`test_scattering_kernel_crosscheck`**. **Pre-existing unrelated failures (NOT ours, verified via
stash):** 8 in `tests/sn/operators/` (2D-mesh `mu_y` cosines; sphere streaming snapshots).

**NEXT = Phase C (consumer migration + retirement)** — see Migration §C below, sharpened:
1. **Build the Frame in `ScatteringOperator`**: `Frame(SphericalHarmonicBasis(self.scattering_order),
   DiscreteMeasure(nodes=<quadrature direction-cosines (N,3)>, weights=self.weights, space="S^2"))`.
   Its `frame.table` MUST equal `self.Y` bit-identically (same nodes + basis) — assert it.
2. **Kernel** (`scattering.py` ~629-664): `OperatorProduct(frame.reconstruction, OperatorProduct(Λ,
   frame.analysis))`; **retire the 3 `cast(LinearOperator,…)`**. NUANCE: `Λ`
   (`LegendreMomentScattering`) is still unparametrised `LinearOperatorMixin` with `None` spaces — it
   may still need ONE cast or its own minimal typing; check the CLI pyright after, keep ≤1 cast only
   if unavoidable (note #226 if so).
3. **`build_aniso_source`** (~949-951): `frame.analysis.apply(angular_flux.values)`.
   **`_aniso_source_from_moment_values`** (~537-541): `frame.reconstruction.apply(scatter.apply(...))`.
4. **In-sweep walk** reads `M.L`/`M.Y` (`operator.py:1108`, `loss_representation.py:3428/3431`,
   `solver.py` ~535-548): migrate the CONSUMER to read `frame.basis.L` / `frame.table` (adapt the
   consumer; do NOT mimic `.L`/`.Y` on the generic face).
5. **Retire** `MomentProjection` + `HarmonicMomentReconstruction` (+ the `two_l_plus_one` field) into
   the basis/frame. Run the 3-search retirement audit (graph callers + grep code/tests/docs + direct
   ctors). Test migration: behavioral→rewire to `frame`/`basis`; characterization (0-ULP crosscheck,
   `Π R=4π I`, windowed≡flat)→keep+rewire. ~25 test refs + docs (`galerkin_projection.rst` etc.).
   The ABCs `ProjectionOperator`→`AnalysisOperator` rename (keep `ReconstructionOperator`) folds in here.
6. **Gate every step on the 0-ULP `test_scattering_kernel_crosscheck` canary** + sn/operators+solve+sweep.

Tasks: #29 (Phase C, in_progress) → #30 (Phase D: reconstruction `apply_transpose`+`CAP_APPLY_TRANSPOSE`
→ `R.H`) → #31 (E+F: solution-as-measure seam + docs + EXIT to #226). Memory:
`[[project-frame-basis-carve]]`, `[[feedback-high-signal-names]]`.

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
