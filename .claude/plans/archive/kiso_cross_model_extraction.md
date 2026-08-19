# K_iso — extract the model-independent isotropic energy operator (+ SN adjoint S†)

> Campaign **#276** (adjoint SN transport), reframes task **#119** (closes **#118** A2b).
> Branch `feature/sn-adjoint-transport`. Surgical, **main-agent-direct, NO method-implementer**;
> per-gate cadence (implement → green gate + forward-safety → present → commit). Verification spec:
> `.claude/plans/archive/a2_kiso_verification.md` (test-architect, 12 gates). **Durable copy** of the
> plan-mode plan `reactive-moseying-cake.md` (which is ephemeral).

## STATUS (2026-06-29 — read first; reconcile vs git, NOT this prose)

Branch `feature/sn-adjoint-transport`, ahead of `main` @ `68ceb9a`. **P4 (meshless homogeneous closeout)
COMPLETE — P4-A…P4-F all committed (NOT pushed; user holds pushes).** NEXT = **P5 (#124)** docs/theory
closeout (archivist) + file CP/MoC/diffusion/MC K_iso-migration issues, then a **MEMORY.md distillation**
(retire the merged P5/P5.5/group-convention "Active work" entries — process-discipline merge-time rule).
Campaign commits: `a8bb027` K_iso P1 · `dcea43a` P2 (SN fwd→K_iso) · `15185e5` P3 (S†, **closes #118**) ·
`a799e92` P4-A (#267 bulk-field retype) · `8a4a094` P4-B (`MaterialMesh.from_materials`) ·
`fe9d4c6` P4-C (homogeneous refound + n2n correction + de-vacuum + docs) ·
`6611f58` P4-C′ (homogeneous A/F from the operators; `MultiplicationOperator` meshless arm) ·
`54b2781` P4-D (`direct_eigenvalue` mint + homogeneous eig/rates rewire) ·
`79e12b4` P4-E (`rayleigh_quotient_iteration` — bordered RQI, #277) ·
`29366b0` P4-F (EnergyGrid geometric diagnostics + homogeneous/MoC migration; MC twin → #278) ·
`195e590` chore (V&V matrix regen).

- **P1 DONE `a8bb027`.** `IsotropicScattering` (Σ_s0) + `IsotropicN2N` (2Σ_2n) — model-shared, sn-free.
- **P2 DONE `dcea43a`.** SN forward iso source → `ScatteringOperator.isotropic_kernel`. 0-ULP.
- **P3 DONE `15185e5` (closes #118).** `S† = (1/W)·full_scatter_kernel.apply_transpose`; reciprocity-pinned.
- **P4-A DONE `a799e92` = the #267 bulk-field mesh retype.** `BulkField.mesh: SNMesh→MaterialMesh`;
  `AngularField`/`MomentField` narrow back to `SNMesh` (covariant frozen-field override); 5 honest casts
  remain (`FullField.bulk: BulkField` consumers — the `FullField.bulk: AngularField` follow-up, #267
  comment `4829138742`, removes them). pyright 410.
- **P4-B DONE `8a4a094` (#126).** `MaterialMesh.from_materials(materials)` — meshless single-cell
  (`AxisMesh(edges=[0,1])`, `mesh=None`, `mat_map=[0]`, ng derived). Sidesteps the `from_axes` clash
  (own name + own trivial axis). Intrinsic test (meshless/1-cell/id-0/unit-volume/xs-field-resolves).
- **P4-C DONE `fe9d4c6` (#127).** Refound `solve_homogeneous_infinite` on **direct λ_max(A⁻¹F)**;
  retired the `HomogeneousSolver` power-iteration shell (kept the public entry + `power_iteration` for
  diffusion). **(n,2n) CONVENTION CORRECTION:** the bespoke added `2·colsum(Σ2)` to PRODUCTION — a
  double-count drift from the theory page + both oracles (n2n is loss-side, in A only; production =
  νΣf). Was vacuous (all fixtures Σ2=0). **De-vacuumed:** new `homo_2eg_n2n` case (asymmetric Σ2≠0,
  oracle k=1.6532258064516119); mutation-proven (double-count Δ0.43/0.45, drop-n2n Δ0.62,
  un-transposed Δ0.16 all RED; Σ2=0 stays green). Fixed `homogeneous.rst:19` doc-drift + the body's
  retired-convention labels + a pre-existing 2g worked-example error (archivist; Sphinx -E -W exit 0;
  all V&V labels kept-and-repointed, no test edits). VACUUM bit-id pre-flight passed.
- **P4-C′ DONE `6611f58`.** **A and F now ASSEMBLED FROM THE TRANSPORT OPERATORS, not np.diag/np.outer**
  (user: "use both operators"). `A = _as_dense(MultiplicationOperator.from_mesh(σt_field) − (IsoScat +
  IsoN2N), ng)`, `F = _as_dense(FissionOperator.from_solver_data(mat_xs))`; `_as_dense` =
  apply-to-basis-vectors (works on the `C − K_iso` OperatorSum). **The one real gap:
  `MultiplicationOperator.apply` was not meshless-safe** (wrapped the meshless engine in FullField +
  cast SNMesh) → gave it a bare-`ndarray` singledispatch arm (`coefficient.values * x`) mirroring
  FissionOperator; FullField arm UNCHANGED (206 callers, 16 tests bit-identical-green); ndarray-only
  scope (ScalarFlux arm = false symmetry, deferred). F needed NO change (was already meshless — just
  bypassed). k∞ bit-identical to oracle (Δ=0 all 4). `TestMeshlessBareArm` (cross-arm consistency the
  load-bearing pin) mutation-proven; A-oracle pin + Mode-11 sentinel. pyright 410; broad sweep 915.

- **P4-D DONE `54b2781` (#128).** Retired the LAST hand-rolling in `solve_homogeneous_infinite`.
  **Decision (user, AskUserQuestion 2026-06-29):** the eigenvalue is NOT routed through `KEigenvalue` —
  the 0-D infinite-medium k∞ is an EXACT dense eigenproblem and `KEigenvalue`/`power_iteration` are
  iterative (approximate-where-exact + would make the 1e-12 analytical gate dominance-ratio-fragile).
  Minted **`numerics.eigenvalue.direct_eigenvalue(A, F) → (k, φ)`** = exact `λ_max(A⁻¹F)`, the direct
  (non-iterative) sibling of `power_iteration` (free function, symmetric; fail-loud on a complex dominant
  = Perron–Frobenius). **VERIFIED against a TRANSPORT-UNRELATED pure-math gate** (`tests/numerics/test_eigenvalue.py`,
  hand-derived closed-form eigenpairs, 24 gates, mutation-proven) per the user's directive — once verified
  this way it is trusted machinery whose independence as a verification tool comes from the CONSUMER
  assembling A/F differently (fused oracle vs operator-algebra), NOT the eigensolver. Rates →
  **`IntegratedReactionRate`** (production=νΣf; absorption now via the intended
  `absorption_cross_section_field` consumer — `cell_xs.py:58` `sig_a[i]=m.absorption_xs` ⟹ bit-identical;
  0-ULP). k∞ bit-identical (test-pinned). Gates: pure-math 24/24 + homogeneous 10/10 (Mode-11 + bit-id
  pins); broad sweep 1248; pyright 410; test files clean. `_as_dense`/`_assemble_loss_matrix` kept as the
  transport-layer operator→dense adapter.

- **P4-E DONE `79e12b4` (#129, issue #277).** Minted the THIRD eigenvalue engine —
  **`numerics.eigenvalue.rayleigh_quotient_iteration(A, F, *, v0, tol, max_iter) → (k, φ)`** = bordered
  Rayleigh-Quotient Iteration (the user's "eigenvector-of-a-previous-iteration as a row" = the
  augmented-Newton system `[F−kA, −Av; vᵀ, 0]`). Polishes an estimate to the NEAREST eigenpair of A⁻¹F,
  superlinearly (locally quadratic). Verified against the now-trusted `direct_eigenvalue` oracle: eigenpair
  match + residual law + QUADRATIC-convergence teeth (error squares each step) + the NEAREST-eigenvalue
  contract (subdominant warm-start → subdominant, documented `.. warning::`). The meshed-operator
  integration (the real payoff) + adjoint φ* vehicle (RQI on the daggered operator) = **#277** (follow-on,
  NOT P4-E). Gates: 30/30 (24 pure-math + 6 RQI); numerics 725; pyright 410.
- **P4-F DONE `29366b0` (#130, issue #278).** Folded the homogeneous (+ MoC plot) energy-grid diagnostics
  onto **`EnergyGrid`** — the arithmetic group midpoint was WRONG on a log axis. Added
  `EnergyGrid.energy_widths` (ΔE) + `lethargy_widths` (Δu); `representative_energy` = GEOMETRIC
  √(E_up·E_lo) already existed. Renamed `HomogeneousResult` fields `eg_mid`→`representative_energy` (now
  geometric), `de`→`energy_widths`, `du`→`lethargy_widths`; solver reads `mix.energy_grid`. Migrated
  `plot_spectrum` + `plot_moc_spectra` (the inline twin). PRINCIPLED change (plot abscissa → geometric
  centre; flux unchanged), not bit-identical. Intrinsic geometry gate (geometric-not-arithmetic via AM-GM)
  + homogeneous eg-block wiring/None tests (the eg-block was previously UNTESTED). MC twin (separate
  `xs.eg_mid`/`xs.du` source) → **#278**. Gates: energy_grid+homogeneous 49; broad sweep 376; pyright 410;
  Sphinx -W; V&V matrix regenerated `195e590`.
- **NEXT = P5 (#124):** docs/theory closeout (archivist — the homogeneous theory page's eigenvalue story is
  now 3 engines + the geometric-centre note; the SN adjoint S† changelog) + file the CP/MoC/diffusion/MC
  K_iso-migration issues. **Then a MEMORY.md distillation** (over size limit; retire the merged
  P5/P5.5/group-convention "Active work" entries per the process-discipline merge-time rule). Downstream:
  the adjoint chain A3 (#114 solve_transpose + `.H` CAP_SOLVE) / A4 (#115 φ*) / A5 (#116 carrier) / A6
  (#117 docs).

> **NOTE:** the original P4 "migrate homogeneous onto `K_iso.as_dense()`" + its R1 `Mixture`-vs-
> `MaterialXSField` FINDING (below) are SUPERSEDED by the reframe — the meshless `MaterialMesh`
> dissolves the bare-`Mixture` problem (a homogeneous medium IS a 1-region `MaterialMesh`).

## REFRAMED P4 — meshless homogeneous on the transport operator algebra

**Vision (user):** solve homogeneous k∞ with the SAME transport operators SN uses —
`(C − S − F/k)φ = 0` with **L (streaming) DROPPED** (zero in an infinite medium). Building C/S/F on
a 0-D meshless phase space is the forcing function that reveals how to instantiate K_iso without a mesh.

**Grounding (2 explorers, 2026-06-28):**
- The operators are meshless-ready: **C** (`MultiplicationOperator`) broadcast-multiply, **K_iso**
  (`IsotropicScattering+IsotropicN2N`) bare-ndarray in/out, **F** (`FissionOperator`) pure rank-1
  einsum (`volume_measure` NOT on the matvec path). `space=None` skips the composition guard.
  The full angular `ScatteringOperator` needs a `Quadrature` — but homogeneous is isotropic ⟹ **K_iso,
  not ScatteringOperator** (no quadrature).
- **The meshless phase space = a 1-cell `MaterialMesh`** (`spatial_shape=(1,)`, `mat_map=None` → single
  region, `mesh=None`). Buildable today via the axis-primary path; exposes `volume_measure` (k∞
  extraction works). `axes=()` truly-0-D BREAKS at a `reduce()` empty-product — use the 1-cell.
- **k∞ extraction** already typed: `IntegratedReactionRate` = `R_{νΣf}/R_{Σa}` (the φ†=1 Rayleigh
  quotient, also the P6 foundation). The analytical oracle (`derivations/common/eigenvalue.py`) +
  the production solver share `A = diag(σt) − Σs0ᵀ − 2Σ2ᵀ`, `F = outer(χ,νσf)` → the operator-composed
  A/F must reproduce k∞ to **1e-12** (`tests/homogeneous/test_homogeneous.py::test_kinf_exact`).

**Phases:**
- **A — DONE `a799e92`** (#267 bulk-field mesh retype; see STATUS). The honest typing that lets
  `MaterialXSField` live on a meshless `MaterialMesh`.
- **B (#126) — NEXT:** mint `MaterialMesh.from_materials(materials)` (single-region 1-cell;
  `material_mesh.py:226` un-defer; `ng` derived — no twin). Intrinsic test: `volume_measure` present,
  `mat_map=[0]`, `material_xs_field()` resolves. Surgical, main-agent-direct.
- **C (#127):** refound `homogeneous/solver.py`: `mat_xs = MaterialMesh.from_materials({0:mix})
  .material_xs_field()`; build `C`, `K_iso`, `F`; loss `A = (C − K_iso)` (L dropped); k∞ =
  dense-A-from-operators `λ_max(A⁻¹F)` (matches oracle) — **sub-decision at the C checkpoint:** that vs
  apply-based power-iteration. GATE: k∞ **1e-12** (1eg/2eg/4eg) + production-rate-100 stay green; fix
  `homogeneous.rst:19` doc-drift (claims `numpy.eigvals` λ_max; code power-iterates + `spsolve`).

**#267 cast follow-up (separate, optional):** the 5 honest casts (`FullField.bulk: BulkField`
consumers) would vanish by narrowing **`FullField.bulk: AngularField`** — broader core-type change,
200+ consumers, own verified pass. Tracked on #267 comment `4829138742`.

## Context — why K_iso, and why NOT a new iso-frame

P6 needs an adjoint flux φ\*; #276 builds S† via the exact discrete transpose. The investigation
(explorer + cross-domain-attacker, `iso_source_frame_conjugation_unification.md`) flipped three premises:

1. **The adjoint S† is ALREADY free + committed.** `full_scatter_kernel` (`scattering.py:850`) *is*
   `frame.conjugate(Λ_{ℓ≥0}+N2N)`; `full_scatter_kernel.apply_transpose` is free via
   `OperatorProduct.apply_transpose` (pinned by `test_full_kernel_euclidean_reciprocity`). So
   `S† = (1/W)·full_scatter_kernel.apply_transpose` works today (P3 just wires it).
2. **A new `ConstantBasis`/iso-frame is the WRONG object** — Y₀⁰=1 IS the ℓ=0 spherical harmonic, so
   an "iso frame" is the ℓ=0 sub-block of the existing harmonic frame; a separate basis FORKS `R∘M`
   and is degenerate (M₀=R₀=Id) in diffusion/CP. The fission dyad is already the normal form.
3. **The real win is model-independence by RELOCATION** — the energy op `Σ_s0ᵀφ` / `2Σ_2nᵀφ` /
   `χ⟨νΣf,φ⟩` is the SAME per-cell `einsum("fg,f->g")` re-implemented **6× / 5×** (SN / CP×3 / MoC /
   homogeneous / diffusion / MC). Only the angular *wrapper* differs. Cardinal Rule 2 over 6 sites.

**Bound (do NOT over-reach):** iso-source = `frame.conjugate` (Galerkin); homogenization/condensation =
`frame.project` = G⁻¹M (Petrov-Galerkin). Different verbs — they do NOT unify.

## The design

```
IsotropicScattering(mat_xs)   # Σ_s0  : apply φ↦Σ_s0ᵀφ ;  apply_transpose χ↦Σ_s0χ ;  dense_per_material→Σ_s0ᵀ
IsotropicN2N(mat_xs)          # 2Σ_2n : apply φ↦2Σ_2nᵀφ ; apply_transpose χ↦2Σ_2nχ ; dense_per_material→2Σ_2nᵀ
```
- **SN forward (P2 — fast, energy op typed):** `φ=integrate(ψ); iso=(IsotropicScattering+IsotropicN2N).apply(φ);
  src=broadcast(iso)/W + aniso`. Replaces the inline `add_iso_source`/`add_n2n_source` in
  `_assemble_per_ordinate_source` (scattering.py ≈885); SAME `mat_xs` verbs ⇒ 0-ULP; aniso
  (`build_aniso_source`, ℓ≥1) UNCHANGED.
- **SN adjoint S† (P3 — closes #118):** `ScatteringOperator.apply_transpose =
  (1/W)·full_scatter_kernel.apply_transpose`; flip `capabilities`; retire the "no apply_transpose"
  confession. Asymmetry (fast-path forward, frame transpose) is principled + oracle-pinned.
- **Cross-model (the payoff):** operators in `transport/operators/` (sn-free); migrate consumers
  incrementally — see the P4 FINDING.

`full_scatter_kernel` stays the permanent oracle (forward + transpose).

## ⚠ P4 FINDING (test-architect, HIGH risk — resolve before P4)

`as_dense` (`dense_per_material`) has **ONE** real consumer: **homogeneous**. But
`homogeneous/solver.py` holds a bare `Mixture` — **no mesh, no `MaterialXSField`** (:87-95) — and the
P1 operators take `mat_xs` (mesh-bound). And **diffusion is NOT a clean consumer** (`solver.py:240` is
matrix-FREE + 2G-hardcoded). So P4 needs EITHER a `Mixture`/`from_materials`-constructible path on the
operators (the genuine model-independence — recommended, small) OR defer P4 + file it (don't force an
SN-mesh import into homogeneous, which would falsify "sn-free"). **Decide at P4 (likely a small
`from_materials(materials, mat_map?)` ctor; homogeneous is 1 material / 0-D).** `dense_per_material`
shape (`dict[int,(ng,ng)]`) is CONFIRMED right (not block-diagonal).

## Phases (each: green gate + forward-safety → present → commit)

- **P1 — DONE `a8bb027`** (operators + scalar transpose verbs + 16 gates).
- **P2 — route the SN forward through K_iso.** `_assemble_per_ordinate_source` → K_iso.apply; aniso
  unchanged. **LD gate re-pin:** keep `_assemble_per_ordinate_source` the typed combiner so
  `test_mms_ld_2d.py:921` still bites; sentinel-confirm it executes (Mode-11), mutation-verify RED.
  GATES (`a2_kiso_verification.md` #3,#5,#10,#12): SN forward 0-ULP bit-identity; aniso 0-ULP canary;
  keff (≥2G het)+SI≡Krylov; LD re-pin.
- **P3 — SN adjoint S† (closes #118).** GATES #4,#11: oracle equivalence (transpose, P0+P1+LD);
  free-transpose Euclidean reciprocity per-group (group-flip/omit-n2n RED); capability + S†-routing sentinel.
- **P4 — cross-model validation via homogeneous** (gated by the FINDING above). GATES #6,#7,#8:
  as_dense≡apply; homogeneous keff invariance (2eg+4eg, closed-form ground); VACUUM bit-id (Sig2=0 ⇒
  as_dense ≡ legacy SigS0ᵀ).
- **P5 — docs + file follow-on migrations.** scattering.py docstrings; file `type:improvement` issues
  for CP (`cp/solver.py`×3) / MoC (`moc/core.py`) / diffusion (`diffusion/solver.py`, n2n=0 + the
  matrix-free obstacle) / MC (`mc/solver.py:369`) with explorer file:line + layout (cell-leading vs
  ng-leading) + keff production-fold notes. Promote `diag_276_full_scatter_kernel_ld_trailing_axis.py`.
  Theory = A6 (archivist).

## Critical files

- **done (P1):** `orpheus/transport/operators/isotropic_scattering.py`,
  `orpheus/transport/mesh/material_xs_field.py` (transpose verbs), `tests/sn/operators/test_isotropic_scattering.py`.
- **P2/P3:** `orpheus/transport/operators/scattering.py` (`_assemble_per_ordinate_source` ≈885 → K_iso;
  `apply_transpose` via `full_scatter_kernel` ≈850; capabilities flip); `tests/sn/verification/mms/test_mms_ld_2d.py:921`
  (re-pin); `tests/sn/operators/test_scattering_adjoint.py` (extend).
- **P4:** `orpheus/homogeneous/solver.py:93` + `tests/homogeneous/test_homogeneous.py` (keff).
- **reuse only:** `orpheus/numerics/operator.py` (OperatorSum/Product/outer propagate apply_transpose).

## Scope / discipline

Surgical, main-agent-direct, NO method-implementer. Commit/push/merge ONLY when asked; stage explicitly
(no `git add -A` — `.claude/*` policy-uncommitted); trailer `Co-Authored-By: Claude Opus 4.8 (1M context)
<noreply@anthropic.com>`; no `# type: ignore` (`cast` OK); NEVER `git checkout/restore/stash` on
uncommitted files (`git reset` index-only safe). **Bit-identity** where the SN iso path is re-expressed;
**principled-equiv** (~1e-12) for oracle cross-check + homogeneous re-baseline. Lands #119 (reframed) +
#118 (S†) + the cross-model K_iso foundation; unblocks A3 (#114). Canonical gate:
`.venv/bin/python -O -m pytest <paths> -m "not slow" -q -rfE --timeout=300 -p no:xdist -p no:cacheprovider`.
