# A5 carrier-verification memo — the sphere φ*-shape row + the carrier testability axis

**Author:** test-architect (chartered input for #276 A5). **Status:** design-time; the
USER rules the carrier shape — this memo is the *testability* input, NOT the ruling.
**Discipline:** no production/test code written; gates specified, mutations named, teeth
proven able-to-red.

Pre-reads grounded by `file:line`: `tests/sn/solve/test_sn_adjoint_certification.py`,
`tests/sn/solve/test_sn_adjoint_entries.py`, `tests/sn/_test_helpers.py`,
`tests/sn/operators/test_loss_transpose_solve.py` (the forward-probe dense-Mᵀ pattern),
`tests/sn/mesh/test_radial_characteristic_carrier.py:264-352` (the ray `G_sd = V_cell`
hand-gauge), `orpheus/sn/solver.py:2211-2601` (the entries + posing + packaging),
`orpheus/numerics/operator.py:1204-1227` (`_AdjointOperator.apply` = `G⁻¹∘Aᵀ∘G`),
`orpheus/numerics/iteration.py:1000-1194` (KEigenvalue triple),
`.claude/plans/archive/p6_adjoint_verification_spec.md:357-420` (#281 C1/C2/C3).

---

## 0. Executive summary (findings + recommendations)

**Ask 1 (PRIMARY) — the sphere φ*-shape row.** RECOMMENDED design: a **defining-law
residual** on the FULL coupled daggered vector, referenced against a **dense
forward-probed** `(A_loss, F)` pair and a **raw-data G-metric**. It is strictly stronger
and more robust than the task's candidate (direct dense-left-eigenvector comparison),
because it is normalization-free, sign-free, degeneracy-robust, and never inverts the
(singular partial-current-trace) metric.

- The coupled operator's domain IS a flat vector space (`template.to_flat().size` +
  `CoupledField.from_flat`) — so the dense forward-probe is FEASIBLE (the task's
  fallback is not needed). Probe `resolvent.apply`, `gain.apply`, `F_posed.apply`
  (all FORWARD — structurally independent of the `.H` reverse-scan under test) via
  `from_flat`/`to_flat`; `A_loss_dense = resolvent_dense − gain_dense`.
- The SN `.H`-daggered ψ* is the **G-metric** adjoint: `A.H = G⁻¹Aᵀ G`
  (`operator.py:1204-1227`), so `v := G·ψ*` solves the true forward-Euclidean-transpose
  adjoint eigenproblem `A_lossᵀ v = (1/k) Fᵀ v`. The row asserts
  `‖A_loss_denseᵀ (G·ψ*) − (1/k) F_denseᵀ (G·ψ*)‖` is at the iteration floor — the
  coupled-sphere generalization of the shipped ∞-medium row
  `test_sn_adjoint_dominant_mode_is_the_true_left_vector` (certification `:413-437`).
- **Cover the RAY (System B) too** — the sphere leg exists to exercise the coupled
  System-A/System-B daggered structure; the full-coupled residual spans it. System-A-only
  is genuinely weaker (blind to a ray-block-isolated transpose error).
- **The Mode-10/Mode-12 tooth = the `G_sd` metric drop** (ERR-067 lineage): mutate the
  ray state-metric `V_cell → Euclidean(1)`. k stays EXACTLY equal (any invertible G gives
  a metric-similar daggered pencil — the k rows are *designed-green* on a metric error),
  while the vector residual (built with the CORRECT `V_cell` G) goes O(1) RED. This is the
  ideal tooth: k-blind, vector-visible, and it exercises the exact sphere-specific
  System-B daggered metric.
- **Fixture:** a small HETEROGENEOUS 2-region sphere (≥2G, asymmetric SigS) — the
  certification's homogeneous `_sphere()` weakens the daggered-vector materiality
  (homogeneous ⇒ spatial φ*≈φ; the adjoint content collapses toward angular-only).
- **Location/marker:** extend the certification file, new class `TestP14SphereAdjointVector`,
  `@pytest.mark.l1`; `np.testing`/`require` only (Mode 8). Runtime target ≤ ~30 s
  (keep the mesh minimal); tag `@pytest.mark.slow` only if measured beyond.

**Ask 2 (SECONDARY) — the carrier testability axis.** The decidable type-minting criterion
(`coding-standards.md`) FAILS for φ*: it is *isomorphic* to φ (same `Solution`/`ScalarFlux`
storage — one realization) and no non-identity morphism is applied to φ* alone that a
distinct type would make unrepresentable. The one latent consumer (#281 adjoint-weighted
`homogenize`) does NOT yet exist and, when it lands, its C1/C2/C3 gates consume φ* as an
**explicit test-weight**, requiring NO production-side role discriminator on `Solution`.
Therefore the testability axis FAVORS **(c) the landed unmarked `Solution`** now, with a
role **property (a)** minted only if #281 P6-B2 rules a role-dispatch API (a choice the
"repeated conditional is a missing type" caution weighs against). The **TYPE (b)** is
weakest: it demands the most net-new gates (constructor/guard) for a concept the criterion
says is not a type, at the highest migration/retirement cost. Full per-shape surface in §3.

---

## 1. System under test and the exact gap (A4, merged @ 065a0e5d)

`solve_sn_adjoint` (`solver.py:2341`) poses the daggered eigenproblem
`A_loss.H ψ* = (1/k) F.H ψ*` as `KEigenvalue(resolvent.H, gain.H, F_posed.H)` →
unchanged `power_iteration`. On a carrying (sphere) mesh, `_adjoint_posing_parts`
(`solver.py:2255-2338`) returns the COUPLED pieces: `resolvent` = the within-group joint
grid `M = [[LC, Seeding],[None, march]]`, `gain` = the coupled `N` (S+B), `F_posed`
= `[[F, ∅],[RadialCharacteristicEmission(F.kernel), Zero]]`. So `A_loss = resolvent − gain`
and `F_posed` are operators on the COUPLED space (System-A `FullField` ⊕ System-B ray).

**What the certification already certifies on the sphere** (`test_sphere_k_equality`,
`:170-188`): `k_adj == k_fwd` ONLY. Its own honest-scope note (`:174-179`) says the coupled
daggered flux SHAPE has no independent reference on this leg.

**Why k can never certify the sphere daggered VECTOR (Mode 12, as measured in the A4
sweep).** `eig(A.H) = eig(A)` by construction, so every k-level functional is
designed-green on the whole adjoint mutation class. Sharper still for the metric: `.H`
= `G⁻¹Aᵀ G` is metric-SIMILAR to `Aᵀ` for ANY invertible `G`, so `k_adj` is EXACTLY
invariant even under a WRONG G-metric — the daggered posing carries a free parameter (G)
that no eigenvalue gate can see. The ∞/slab vector rows (P1.4 spectrum, P1.5 bi-orthogonality)
do NOT lift to the sphere: the closed-form `(Aᵀ)⁻¹Fᵀ` spectrum reference is a 0-D
energy-only object with a trivial metric; the sphere adds streaming, curvilinear angular
coupling, the ψ½ ray, and a NON-trivial block G-metric (bulk `V·w_n` ⊕ trace `|Ω·n|·w_n`
⊕ ray `V_cell`). **The sphere daggered vector is currently un-referenced — that is the A5
carve-in gap.**

---

## 2. Ask 1 — the sphere φ*-shape row

### 2.1 The reference: dense forward-probe of `(A_loss, F)` + a raw-data G-metric

Three structurally-independent grounds, none of which touches the `.H`/reverse-scan/
transpose machinery under test (the L11 / anti-R1 discipline; the same posture as the G2
dense-`Mᵀ` oracle in `test_loss_transpose_solve.py`, which probes the FORWARD `apply`):

1. **Dense forward matrices via `from_flat`/`to_flat` column-probe.** For each unit vector
   `e_i` (`i` over `n = template.to_flat().size`): `x = CoupledField.from_flat(e_i, template)`,
   column = `op.apply(x).to_flat()`. Build `resolvent_dense`, `gain_dense`, `F_dense`;
   set `A_loss_dense = resolvent_dense − gain_dense`. The `from_flat`/`to_flat` protocol
   makes this layout-agnostic — it works identically for the coupled `CoupledField` and the
   non-carrying `FullField`, so NO manual seed-DOF bookkeeping is needed (a genuine
   simplification over `_probe_augmented_matrix_one_group`, which is one-group and hand-lays
   the seed block). `resolvent.apply` is the FORWARD walk (`_apply_walk`); the reverse walk
   (`loss_action_transpose`) is the thing being verified — different code path, so the probe
   is independent.
2. **Euclidean transpose via numpy `.T`** — `A_loss_dense.T`, `F_dense.T` — NOT the
   production `apply_transpose`. (Same idea as P1.5's `A.T @ spec`, `:431`.)
3. **The G-metric from raw mesh data** — `g_bulk_measure(sn)` (`V·w_n`), per-face
   `g_trace_cosine_weight` (`|Ω·n|·w_n`), and the ray `V_cell` gauge hand-built from
   `sn.volumes` exactly as `test_radial_characteristic_carrier.py:334-343` does — assembled
   into a flat `g_diag` in `CoupledField.to_flat()` order (interior ⊕ boundary ⊕ ray-interior
   ⊕ ray-boundary). This replicates `_AdjointOperator`'s `apply_metric` (`operator.py:1223`)
   from raw data, so a metric bug in production is caught rather than inherited. RECOMMEND a
   shared helper `g_coupled_diagonal(sn)` in `_test_helpers.py` (extends the existing g_inner
   cluster; anti-R1 — never read the production Gram).

**Probe-validation (anti-vacuity leg, cheap and mandatory).** The dense probe is only a
trustworthy reference if it reproduces the FORWARD physics. Assert the dense forward
dominant eigenpair matches the SN forward solve: `k_dense = ρ(np.linalg.solve(A_loss_dense,
F_dense))` equals `solve_sn(sphere).keff` to ~1e-8. (Optionally, the dense forward
eigenvector matches `solve_sn`'s forward flat flux up to scale — a bonus, fiddlier
normalization.) This pins that `A_loss_dense`/`F_dense` are assembled correctly before they
are transposed for the adjoint reference.

### 2.2 The comparison: the coupled defining-law residual (WHY it is the right granularity)

`_AdjointOperator.apply` realizes `A.H y = G_dom⁻¹ · apply_transpose(G_cod · y)`
(`operator.py:1204-1227`), i.e. `A.H = G⁻¹Aᵀ G` with `G` owned by the function space (block
`G = diag(G_A, G_sd)`, `G_A` = bulk `V·w_n` ⊕ trace `|Ω·n|·w_n`, `G_sd` = ray `V_cell`).
The daggered eigenproblem `A_loss.H ψ* = (1/k) F.H ψ*` therefore rearranges (left-multiply
by G, using `A_loss.H = (A_loss)ᵀ` conjugated by G) to

```
A_lossᵀ (G·ψ*)  =  (1/k) Fᵀ (G·ψ*).
```

So `v := G·ψ*` is the textbook forward-Euclidean-transpose adjoint eigenvector. The row:

```
psi_flat = psi_star.to_flat()                 # the SUT's converged daggered iterate
v        = g_diag * psi_flat                   # G·ψ*  (elementwise — G is diagonal)
resid    = A_loss_dense.T @ v - (F_dense.T @ v) / k_adj
rel      = ||resid|| / ||(F_dense.T @ v) / k_adj||
np.testing.assert_allclose(rel, 0.0, atol=TOL, ...)      # see §2.5
```

**Why the residual, not a direct eigenvector comparison** (the task's candidate):

- **Normalization/sign-free.** The residual is invariant to the power-iteration's arbitrary
  scale and sign of ψ* — no sign-fixing, no ℓ² renorm, no clustered-eigenvector ambiguity.
- **Never inverts G.** The trace metric `|Ω·n|·w_n` is SINGULAR at grazing ordinates (the
  "pseudo-inverse on the singular partial-current trace", `operator.py:1218`). The residual
  only MULTIPLIES by G — a direct comparison `ψ* ?= G⁻¹v` would hit that pseudo-inverse and
  be ill-posed on the trace block.
- **The G-map is load-bearing and self-checking.** If `.H` erroneously dropped the metric
  (used Euclidean transpose), ψ*_wrong would solve `A_lossᵀ ψ* = (1/k)Fᵀ ψ*` and the residual
  built with the CORRECT G would red (because `G·ψ*_wrong` no longer solves it). Conversely,
  testing the law WITHOUT the G-map would falsely-red the correct G-adjoint ψ*. So the
  raw-data G is exactly the discriminator between "right adjoint" and "metric-dropped adjoint".
- **Dominance is covered by the k rows.** The residual passes for ANY adjoint eigenvector
  (dominant or sub-dominant); a sub-dominant ψ* would shift `k_adj ≠ k_fwd` and the existing
  sphere k row catches it. Residual (vector-correct) + k row (dominant) are airtight together
  — the Mode-12 "pair every k row with a vector gate" discipline, closed in both directions.

**Granularity: full coupled ψ* (per-ordinate System-A bulk + trace + System-B ray), NOT a
reduced scalar φ*.** Reason: on a vacuum/heterogeneous sphere the daggered content lives in
the ANGULAR structure (importance is μ-asymmetric: inward-directed neutrons matter more —
the μ-reversal the sphere leg exists to test) and in the ray block. A scalar reduction
`φ* = Σ_n w_n ψ*_n` integrates over μ and would wash the μ-reversal out; it also drops
System B entirely. The residual on the full flat vector is per-ordinate/cell-resolved and
spans the ray — the strongest available scope.

### 2.3 Fixture — a small HETEROGENEOUS sphere (materiality + runtime)

The certification `_sphere()` (`:116-122`) is HOMOGENEOUS single-region 2G. That weakens
this row: on a homogeneous bare sphere the forward and adjoint SPATIAL shapes coincide
(both the fundamental mode), so the daggered-vector materiality collapses toward
angular-only, and the mutation tooth's visibility (Mode 10) weakens with it. RECOMMENDED
fixture for THIS row:

- **2-region sphere** (fuel core + reflector — e.g. `mat_ids = [0,0,0,1,1,1]`, Mixture
  A/B 2g), `coord=SPHERICAL`, `bc_right=vacuum`, `bc_left` reflective at r=0 (the pole).
  Heterogeneity makes φ* differ from φ SPATIALLY as well as angularly — the daggered vector
  is genuinely non-trivial (vv §0.6 / H2: homogeneous nulls spatial-distribution content).
- **≥2G with asymmetric SigS** (Mixture A 2g already asymmetric — the P1.3 S† mutation
  asserts `not allclose(sigs, sigs.T)`, `:230-233`). Assert the materiality precondition
  in-row: the forward and adjoint fluxes differ (`not allclose(φ, φ*)`) so the row is not a
  self-adjoint dud (mirrors P1.4's `ψ*_cf ≠ φ_cf` guard, `:285-289`).
- **Minimal size** — 6 cells, GL n_ordinates ≤ 8 — to bound `n = template.to_flat().size`
  (bulk `N·ng·nx` + trace + ray) at a few hundred, keeping the O(n) probe (≈ 3n forward
  coupled matvecs) well under the runtime target.

### 2.4 Runtime and feasibility

- Feasible: the coupled domain is a flat space (`to_flat`/`from_flat`), so column-probing
  with unit vectors is well-defined. The task's "infeasible → fallback" branch does not fire.
  **MEASURED (2026-07-25 probe on the recommended fixture):** the 2-region GL-8 6-cell sphere
  carries System B; `template` is a `CoupledField` of flat dimension **n = 140**;
  `CoupledField.from_flat(e_i, template)` → `resolvent.apply(x).to_flat()` roundtrips to a
  140-vector column; `resolvent.is_invertible` is True and `gain`/`F_posed` expose `apply`.
  So the whole reference is `3n = 420` cheap forward block-matvecs + two `140×140` dense
  solves — trivially fast; the daggered solve (the SUT) dominates the row's wall time.
- Cost: `3n` forward coupled block-matvecs (cheap, no inner solve) + the daggered solve
  (the SUT) + two dense `n×n` linear-algebra ops (`n ~ few hundred` → trivial). Estimated
  ~5–20 s for the recommended minimal mesh. Keep it a single fixture; if measured beyond
  ~30 s, add `@pytest.mark.slow`. Do NOT grow the mesh to chase precision — tighten the
  solve tolerances instead (§2.5).

### 2.5 Tolerance (justified from iteration + probe precision, never tuned-to-pass)

The measurand is the RELATIVE defining-law residual
`rel = ‖A_lossᵀ(G·ψ*) − (1/k)Fᵀ(G·ψ*)‖ / ‖(1/k)Fᵀ(G·ψ*)‖`. Its floor is set by three
independent, dimensionally-explainable sources, NOT by fitting:

- **The daggered iterate's convergence.** ψ* is converged to the power-iteration's relative
  `flux_tol`; the residual inherits that distance-to-the-eigenvector, amplified by the
  spectral-gap factor `1/(1−dominance_ratio)` (Trefethen–Bau §27) and the operator condition
  number. For a small bare sphere the amplification is O(1–10).
- **The dense probe.** `A_loss_dense`/`F_dense` are exact to ~1e-13 (machine-precision unit
  probes), negligible vs the iterate floor.
- **The G-metric.** Raw-data exact.

RECOMMENDATION: run THIS row's daggered solve at the certification's tight knobs or tighter
(`keff_tol=1e-10, flux_tol=1e-9, inner_tol=1e-11, max_inner≥500`) so the residual floor sits
near `~1e-9`; gate `atol = 1e-7` on `rel` (matching the shipped ∞-medium defining-law row
`test_sn_adjoint_dominant_mode_is_the_true_left_vector`, `:432-436`, which uses `atol=1e-7`).
The mutation (§2.6) drives `rel` to O(1) — a ≥6-order margin over the gate, so the gate is
never near its own noise. **Measure-then-gate discipline:** the implementer MUST print the
converged `rel` and confirm it is ≥100× below `1e-7`; if it is not, the fix is to tighten the
SOLVE, never to relax the gate (vv rule: a tolerance is a contract).

### 2.6 The Mode-10 / Mode-12 mutation tooth — the `G_sd` metric drop (k-blind, vector-visible)

The row's teeth MUST be proven able-to-red (`§0.5` standing discipline). The ideal tooth is
the one the task asks for: reddens THIS row while staying green on k.

**PRIMARY tooth — drop the ray state-metric `G_sd = V_cell → Euclidean(1)`.** In-process,
monkeypatch the ray space's `inner_product_weights` (the surface
`test_radial_characteristic_carrier.py:297-301` reads: `interior.inner_product_weights` /
`boundary.inner_product_weights`) to all-ones (Euclidean) BEFORE `solve_sn_adjoint` builds
the posing, so the daggered `.H` (`operator.py:1223` `apply_metric`) sees the wrong
System-B metric. Then:

- **k stays EXACTLY equal** (`k_adj == k_fwd`): `A.H = G'⁻¹Aᵀ G'` is metric-similar to `Aᵀ`
  for ANY invertible `G'`, so `eig` is untouched — a DIRECT demonstration that the k rows
  are designed-green on a metric error (assert `abs(k_mut − k_fwd) < 1e-9` INSIDE the tooth).
- **the vector residual reds O(1)**: the daggered ψ*_mut now carries the wrong ray metric,
  and the residual built with the CORRECT raw-data `V_cell` G no longer vanishes — assert
  `rel_mut > 1e-3` (huge margin over the `1e-7` gate).

This is the strongest possible tooth for this row: it is EXACTLY k-blind (proving the row
sees content the k rows cannot), and it exercises the precise sphere-specific System-B
daggered metric — the ERR-067 / diag_gsd family, where the "ghost `G_sd = 0`" was a real
production defect (1.3e-2) that put the seed in `null(G)` and made `A.H` a wrong adjoint for
any nonzero seed (`test_radial_characteristic_carrier.py:274-278`). Note the ghost value 0 is
singular — use Euclidean(1) (invertible) so the k-invariance is exact; a 0-metric would break
the similarity and could shift/NaN k, muddying the "k-blind" demonstration.

*Injection caveat (flag for the implementer):* confirm the `.H` metric reads THROUGH the
patched surface at solve time — if the space caches a Gram at construction, patch the
class-level property (so the fresh `solve_sn_adjoint` build sees it) rather than an instance.
The certification file already uses the instance-dispatch-lambda `monkeypatch.setattr(Class,
method, lambda self, ...: ...)` idiom for the singledispatch descriptors (`:213-215`); the
metric surface is a space property, so patch it at the space class.

**SECONDARY teeth (coarser, NOT k-blind — supplementary Mode-10 evidence).** `F† → F`
(`monkeypatch FissionOperator.apply_transpose = lambda self, x: self.apply(x)`, the
`:213-215` idiom) and `S† → S` (`ScatteringOperator.apply_transpose`) both reddened the
residual — but they ALSO shift k on this heterogeneous sphere (a one-factor transpose drop is
not a pencil similarity, `:336-346`), so they are not "k-blind". Include ONE of them as a
second Mode-10 tooth to show the residual reddens on a leaf-transpose drop as well, but label
it explicitly as k-visible (the sphere k row also catches it). The `G_sd` drop is the
distinguishing tooth; a leaf-drop is corroboration.

### 2.7 Mode-8, location, marker

- **Mode 8:** `np.testing.assert_allclose` / `pytest.fail` (via the file's `require`) only —
  the canonical `python -O` strips bare `assert` to a no-op. The file already banned bare
  asserts (`:42-45`); follow it.
- **Location:** extend `tests/sn/solve/test_sn_adjoint_certification.py`, a NEW class
  `TestP14SphereAdjointVector` sited beside the P1.4 spectrum section (it is the sphere
  flux-shape claim, the sibling of the ∞-medium P1.4/P1.5 vector rows). Reuse the file's
  `_solve_fwd`/`_solve_adj`/`_mix_arrays`/`require`; add the 2-region sphere builder and the
  `g_coupled_diagonal` helper (the latter in `_test_helpers.py` beside `g_bulk_measure`).
- **Marker:** `@pytest.mark.l1` (a flux-shape claim against a structurally-independent dense
  reference — an L1 equation-verification, exactly like P1.4, NOT foundation). The file
  currently partitions 10 l1 + 2 foundation; this adds to the l1 set. Optionally
  `@pytest.mark.verifies("<adjoint-eigenproblem label>")` — link to the daggered-eigenproblem
  equation label when A6/ch15 mints the adjoint theory page (`methods/sn/adjoint.rst`); until
  then leave it un-`verifies`'d with a comment (the `tests/_harness/` registry tolerates an
  l1 row with no equation link — it surfaces as a coverage note, not an orphan).

### 2.8 Honest scope + the direct-eigenvector alternative + feasibility fallback

- **Honest scope of the residual row:** it certifies that the FULL coupled daggered vector
  ψ* (System-A bulk + trace + System-B ray) is the correct discrete G-adjoint eigenvector —
  it solves the true forward-Euclidean-transpose adjoint eigenproblem under the independently
  built G-metric, up to the daggered normalization (to which it is invariant). It does NOT
  independently pin the eigenVALUE (k rides the existing sphere k row) and does NOT compare
  the daggered vector against a *different-structural-angle* semi-analytical sphere flux
  (none exists — streaming + curvilinear coupling admit no closed-form `(Aᵀ)⁻¹Fᵀ` on the
  sphere; the dense forward-probe IS the structurally-independent ground).
- **The direct-eigenvector alternative (the task's candidate), for the record:** build
  `K_dense = solve(A_loss_dense, F_dense)`, take the dense dominant LEFT eigenvector
  `u` (= dominant right eigenvector of `K_dense.T` via `scipy.linalg.eig`); the SN daggered
  ψ* relates to it by `G·ψ* ∝ A_loss_denseᵀ·u`-family (the textbook `v = G·ψ*`), so a direct
  comparison must (i) apply `G⁻¹` — hitting the singular trace pseudo-inverse — and (ii)
  sign/scale-fix a possibly-clustered eigenvector. Both are avoided by the residual. Offer the
  direct comparison ONLY as an optional strengthening leg (it adds "ψ* is the DOMINANT mode",
  which the residual+k-row pair already establishes); it is not required and carries more
  numerical fragility.
- **Feasibility fallback (not needed here, stated per the charter):** were the coupled domain
  NOT a flat probeable space, the best feasible alternative would be exactly this
  defining-law residual computed against an INDEPENDENTLY assembled dense `(A, F)` on the
  coupled space plus a bi-orthogonality row — i.e. the residual formulation IS the robust
  fallback, and it happens to also be the best primary. The coupled domain here IS flat
  (`to_flat`/`from_flat`), so the dense probe stands.

---

## 3. Ask 2 — per-shape verification surface for the φ* carrier ruling

**Decidable criterion (`coding-standards.md` "Type vs property").** Mint a TYPE iff
(a) ≥2 non-isomorphic realizations of the concept AND (b) a non-identity morphism is
actually applied to it; else it is a PROPERTY — and a property needs an actual
discriminating consumer to exist. Applied to φ*: it has ONE realization (a `Solution`
whose `scalar_flux` is a `ScalarFlux`, angular a `TimedFullField` — byte-for-byte the same
storage as the forward φ; `solver.py:2592-2600` routes BOTH through the shared
`_package_solution`), and no morphism is applied to φ* *alone* — the adjoint-weighting
`⟨φ*, Σφ⟩` is a bilinear applied to the PAIR at a call site. So condition (a) FAILS and
(b) is a call-site role, not a type morphism. **The criterion points at property/landed,
not type.**

**The one latent consumer — #281 adjoint-weighted `homogenize`** (`solution.py:368`,
currently the Galerkin degenerate φ*=φ; the frame comment `:461-471` self-describes it).
Its gates (`p6_adjoint_verification_spec.md:364-431`) ALL consume φ* as an **explicit
weight array**, not a typed/flagged object:
- **C1** — `Σ_R = ⟨φ*, Σφ⟩_R / ⟨φ*, φ⟩_R` vs a per-region loop with φ* as the test weight;
  the frame takes the scalar-flux `.values`.
- **C2** — adjoint-weighted homogenize → coarse-resolve keff-gap shrinks; mutation = swap
  φ* for φ as the test weight.
- **C3** — Mode-11 sentinel asserts the captured test-basis weight `np.array_equal` the
  ADJOINT flux (not the forward). It discriminates by VALUE, never by a role flag/type.

So the discriminator #281 actually needs is "a way to PASS the adjoint flux as the test
weight" (P6-B2's deferred parameter-vs-sibling-method API), NOT a self-describing carrier.
This is the load-bearing testability fact for the ruling.

### (a) Property — role marker / optional adjoint-partner field on `Solution`

- **Migrating gates:** `test_solution_packaging_contract` (entries `:159-197`) gains a
  role-stamp assertion (`fwd.role=="forward"`, `adj.role=="adjoint"`). All k/spectrum/
  residual rows unchanged (they read `keff`/`scalar_flux`/`angular_flux`).
- **NEW gates the shape demands:** (1) a role-stamp test — `solve_sn`→forward,
  `solve_sn_adjoint`/`_fixed_source`→adjoint; (2) IF (and only if) `homogenize` is
  redesigned to DISPATCH on the role, a Mode-11 sentinel that the role selects the
  PG-adjoint path. But #281 C1/C2/C3 as specified do NOT dispatch on a role — so this new
  gate exists only under a role-dispatch API redesign.
- **Production discriminator required by a planned gate?** NO. A role property IS a
  discriminator, but no gate on the horizon requires it; adding it now installs a field with
  no discriminating consumer — which `coding-standards` explicitly forbids ("a property needs
  an actual discriminating consumer to exist"). It also invites the `coding-elegance`
  "repeated conditional is a missing type" / boolean-flag-dispatch smell if homogenize ever
  branches `if weight.role == "adjoint"`.
- **Migration risk:** LOW-MEDIUM. Backward-compatible (optional, default forward/None), but
  the field is dead weight until a role-dispatch consumer lands, and a role flag that a
  consumer branches on is a design smell.

### (b) Type — `AdjointFlux` / `AdjointSolution`

- **Migrating gates:** ALL adjoint gates add `isinstance(adj, AdjointSolution)`;
  `test_solution_packaging_contract` reconciles its `TimedFullField`/`AngularFlux`/
  `ScalarFlux` interior assertions with the new wrapper; every certification row that reads
  `adj.keff`/`adj.scalar_flux`/`adj.angular_flux` needs those accessors on the new type
  (inherit `Solution` or re-expose).
- **NEW gates the shape demands:** constructor/guard tests (the daggered entries return it;
  it can't be built from a forward solve); an illegal-state test (a forward-only operation
  refuses an `AdjointSolution`, IF any such operation exists); homogenize type-guards the
  test weight to `AdjointSolution`. This is the ONLY shape that makes "pass φ where φ* is
  required" a TYPE error caught in production.
- **Production discriminator required?** The type IS the discriminator — but the type-minting
  criterion FAILS (φ* isomorphic to φ; no morphism on φ* alone). Per `coding-standards`,
  minting it here is "theatrics: ceremony + a conversion seam without making any illegal
  state unrepresentable" UNLESS the USER judges the forward/adjoint mix-up a real, recurring
  hazard worth an illegal-states-unrepresentable guard (`coding-elegance` Pattern 4). That is
  a CORRECTNESS judgment, not a testability one — I flag it; the USER rules it.
- **Migration risk:** HIGHEST. Touches every adjoint gate's type assertions, the return
  annotation, and the downstream homogenize signature; introduces a `Solution↔AdjointSolution`
  conversion seam; and if a later review rules "property after all", retiring the type is a
  full churn cycle (retirement + test migration).

### (c) Landed state — the unmarked second `Solution`, role at the API/call-site

- **Migrating gates:** NONE — every A4 gate (entries + certification) already asserts on the
  unmarked `Solution`, and §2's sphere row does too.
- **NEW gates the shape demands:** NONE for the carrier. #281's C1/C2/C3 work with
  `forward_solution.homogenize(coarse, adjoint=adj_solution)` (or the sibling-method form):
  the call site knows the role (it called `solve_sn_adjoint`), and C3's sentinel checks the
  weight VALUES, not a flag.
- **Production discriminator required?** NO. C1/C2/C3 discriminate by WHICH ARGUMENT the
  flux is passed as and by its VALUES — (c) satisfies the entire planned gate table.
- **Migration risk:** LOWEST now. Standing cost: "role at the call-site can drift" — a caller
  could pass a forward `Solution` as the `adjoint=` weight and no production guard catches it
  (the bare-carrier hazard, `coding-elegance` #13). In TEST, C1's value-difference gate + C3's
  Mode-11 sentinel catch a wrong weight; production stays unguarded until/unless a type is
  judged warranted.

### Testability-axis recommendation (INPUT — the USER rules the shape)

The testability axis FAVORS **(c) the landed unmarked `Solution` now**, because NO gate on
the horizon — neither A4's nor #281's C1/C2/C3 — requires a production-side role
discriminator on `Solution`: the adjoint-weighted consumer takes φ* as an explicit
test-weight argument and discriminates by VALUE. Mint the **property (a)** ONLY IF #281
P6-B2 rules a role-DISPATCH API (a choice `coding-elegance`'s "repeated conditional is a
missing type" weighs against — the cleaner realization is an explicit `adjoint=` argument
with no branch). The **type (b)** is weakest on the testability axis: the type-minting
criterion's two conditions are unmet (φ* is isomorphic to φ; no morphism on φ* alone), it
demands the most net-new gates for a concept the criterion says is not a type, and it carries
the highest migration/retirement cost — it is warranted ONLY if the USER additionally judges
"forward flux where importance is required" a real illegal state worth a type guard, which is
a correctness ruling I flag but do not make.

**Couple this ruling to the #281 P6-B2 API decision** — they are the same seam. And note the
FLIP trigger for a future type: a SECOND non-isomorphic realization of the adjoint flux
(e.g. an angular adjoint current, or a fixed-source importance map with distinct storage)
plus a morphism applied to it would satisfy the criterion and earn the type; until then it is
premature (`coding-elegance` defer-until-≥2).

---

## 4. Self-improvement check (fired before delivery)

- **New failure mode?** NO new row for `vv-principles`. The sphere row's core lever —
  "`k_adj` is EXACTLY invariant under ANY invertible G-metric (metric-similarity of
  `A.H = G⁻¹AᵀG`), so a metric error is k-blind but vector-visible" — is an INSTANCE of
  Mode 12 (the ERR-067 metric-repair family already lives in Mode-12's "SECOND closure
  mechanism"). The sharpening worth carrying is: **the G-metric is a free parameter of the
  daggered posing that no eigenvalue gate can ever see, so the daggered VECTOR row is the
  sole catcher for a metric bug, and the metric-drop is the canonical k-blind vector tooth.**
  This extends my adjoint-verification recipe family (L17 reverse-scan / L18 carrier-ψ½-metric
  / L19 swap-law / the residue spec) — I will sharpen the lessons memory with the
  coupled-adjoint-vector recipe (defining-law residual + raw-data G + `G_sd`-drop tooth),
  NOT the ephemeral A5 state.
- **Plan rejection?** N/A (design-time; no rejection to log yet).

