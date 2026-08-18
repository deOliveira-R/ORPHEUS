---
name: numerical-bug-signatures
description: 'Use when triaging a "wrong answer" or "all tests pass but..." report in a numerical solver. Provides recognition signatures for known bug classes (cylindrical-DD divergence under refinement, MOC weight cancellation, scattering matrix transpose convention, quadrature weight sums, alpha-dome positivity), each linked to the test that catches it and the ERR-NNN entry in `docs/theory/verification/error_catalog.rst`. Examples: "k_eff diverges as I refine", "homogeneous case passes but multi-region fails", "transport sweep gives wrong answer in curvilinear geometry only". Preloaded by qa, numerics-investigator, and test-architect.'
---

# Numerical Bug Signatures — recognition catalog for plausible-wrong solver bugs

This skill is a *recognition* catalog. Plausible substitution errors —
sign flips, missing factors, transposed operands, wrong recurrences,
index drift, convention drift — are the dominant failure mode of
AI-generated numerical code. They produce answers that *look right*
on the most-run tests (homogeneous, single-group, low-precision)
and only declare themselves under specific stress configurations
(heterogeneous, multi-group, mesh-refined).

The catalog below pairs each signature with its catching test and
its ERR-NNN entry in `docs/theory/verification/error_catalog.rst`. Use it to
**recognise** a reported symptom — not to apply blind fixes.

## When to use

- A solver passes its primary tests but a user report says "wrong
  answer" or "diverges under refinement" or "matches in 1G fails in
  2G."
- You see a code change that touches angular weights, scattering
  matrices, curvilinear geometry, or quadrature normalization, and
  you want to know which regression test must run before merge.
- You are designing the verification matrix for a new solver and
  want to seed it with the failure modes that have already escaped
  review at least once in this codebase.
- A V&V audit reports an equation as "implemented" but the catching
  test is implicit — this catalog tells you which test was supposed
  to catch which class of bug.

## When NOT to use

- The bug is catastrophically wrong (NaN, negative k_eff, off by
  orders of magnitude). Use direct traceback / `nexus-debugging`
  instead — these signatures are for subtle, plausible-wrong bugs.
- The solver disagrees with reference by 1–10 % and the suspect
  region has many interacting factors. Use `probe-cascade` to
  isolate the offending factor first; come back here only after
  the failing factor has been pinned to angular-weight,
  scattering-convention, or geometry-coefficient code.
- You are doing fresh exploration of unknown code. Use
  `nexus-exploring` first; this catalog assumes you already know
  which solver/sweep/operator is suspect.
- You want to add a *new* bug pattern. See "Add-a-signature
  protocol" at the end of this file.

## Recognition workflow

```
1. Match the reported symptom to "Symptom" lines below.
2. Read the "Mechanism" to confirm the candidate signature is plausible.
3. Run the "Diagnostic probe" — the cheapest test that confirms or rules it out.
4. If confirmed: the "Catching test" is the regression that should now fail;
   the "Catalog entry" gives the historical context and lesson.
5. If ruled out: re-read "Symptom" lines for the next candidate.
```

A symptom that matches **two** signatures at once (e.g. heterogeneous
keff drifts AND curvilinear flux spike at r=0) is usually a single
compound bug — fix one, the other often persists or transforms.
ERR-006 was such a case (α recursion + ΔA/w simultaneously).

---

## Reading the 2-D fingerprint — sign × magnitude scaling

Errors carry diagnostic information beyond magnitude. **Sign direction**
(under-prediction vs over-prediction), **magnitude scaling with a
parameter** (mesh refinement, σ_t, geometry parameter), and **regime
gating** (rank ≥ 2 only, heterogeneous only, σ_t·R ≥ 10 only) together
form a 2-D fingerprint that pins the bug class before debugger steps.
**Read fingerprints before opening mpmath.**

Worked examples:

- **Sign positive, magnitude `+57%` with rank ≥ 2 only, scaling as
  `(ρ_max/R)²`** → fingerprint of a missing Jacobian. The sign tells
  you "extra flux," the rank-gating tells you "rank-dependent term,"
  the scaling tells you "geometry-coupled."
- **Sign negative, growing with mesh refinement, homogeneous-exact** →
  fingerprint of a redistribution-term error (sign flip or missing
  ΔA/w). See ERR-006.
- **Sign mixed, oscillating with iteration** → fingerprint of an
  iteration-scheme bug (BC ordering, normalisation), not a discrete
  operator bug. See ERR-003, ERR-004.

When the fingerprint matches a catalogued signature below, use the
"Catching test" line. When the fingerprint is new, the
[../probe-cascade/SKILL.md](../probe-cascade/SKILL.md) skill isolates
the offending factor.

---

## Signature 1: Curvilinear sweep divergence under refinement

- **Symptom:** k_eff or fixed-source flux *diverges* (not converges)
  as the radial mesh is refined. Cell-0 flux in spherical/cylindrical
  geometry shows 35–50 % error and grows with `nx`. Homogeneous
  case is exact; heterogeneous keff drifts (e.g. 1.15 → 0.90 → 0.52
  → 0.25 across refinements).
- **Mechanism:** Curvilinear SN sweeps depend on two coupled angular
  redistribution coefficients: the α-recursion (which advects the
  half-angle ψ across the level) and the geometry factor `ΔA_i/w_m`
  (which scales redistribution by the per-ordinate area swept).
  Either bug alone is invisible to homogeneous problems because
  flat flux makes the redistribution term identically zero for every
  ordinate. Only a heterogeneous, mesh-refined run drives the
  redistribution term out of cancellation.  The WDD angular closure
  used by sweeps converges to a stable but incorrect non-flat
  fixed point — global balance is satisfied, the spatial profile is
  wrong.
- **Diagnostic probe:** Per-ordinate flat-flux residual test.
  Set Q uniform, Σ_t uniform, ψ = Q/Σ_t for every ordinate, then
  check that streaming + redistribution = 0 *per ordinate*, not
  just summed. Conservation always passes; per-ordinate exposes
  the bug.
- **Catching test:**
  - `tests/sn/test_quadrature.py::TestL0TermVerification::test_per_ordinate_flat_flux_consistency`
    (the L0-SN-003 definitive curvilinear diagnostic — flat ψ must
    satisfy streaming + redistribution = 0 per ordinate; tagged
    `@pytest.mark.catches("ERR-006", "ERR-007")`).
  - `tests/sn/test_sweep_operator_inconsistency.py::test_spherical_sweep_vs_bicgstab_flat_flux`
    (catches the WDD-closure variant — sweep diverges, BiCGSTAB
    is exact).
- **Catalog entry:** ERR-006 (α recursion + ΔA/w) and ERR-026
  (WDD curvilinear sweep wrong fixed point); ERR-007 is the
  BiCGSTAB-operator twin of ERR-006.
- **Why it hides:** Homogeneous eigenvalue is exact (flat flux ⇒ zero
  redistribution); 1-group is degenerate (k = νΣ_f/Σ_a is shape-
  independent); particle balance is exact (telescoping is by
  construction); flux non-negativity passes (no NaN until many
  iterations); single sweep is finite. ERR-006 hid behind 20
  passing tests.

---

## Signature 2: MOC angular weight cancellation in homogeneous-only tests

- **Symptom:** MOC k_eff is correct for homogeneous cases (any group
  count) but catastrophically wrong for heterogeneous: e.g. 1-group
  2-region pin cell gives k_eff = 1.344 vs CP reference 0.902.
- **Mechanism:** The MOC scalar-flux update (Boyd et al. 2014, Eq. 45)
  carries factors `4π · ω_a · ω_p · t_s · sin(θ_p)`. The `4π`
  comes from the angular flux → scalar flux integral; the
  `sin(θ_p)` comes from the 2D-segment to 3D-path projection.
  Both factors multiply `δψ` (the angular-flux *increment* across
  a track segment). For converged homogeneous boundary fluxes,
  `ψ_in = Q/Σ_t` everywhere, so `δψ ≡ 0` and the entire weight
  factor is irrelevant. Only spatially non-uniform ψ activates
  the correction term, exposing wrong weights.
- **Diagnostic probe:** Inject a non-trivial boundary flux into a
  pure-scatterer single-sweep test and compare scalar flux against
  the closed-form analytical value. Or: run a 1G 2-region pin
  cell against a CP reference.
- **Catching test:**
  `tests/moc/test_verification.py::TestL0EquilibriumFlux::test_pure_scatterer_equilibrium_single_sweep`.
- **Catalog entry:** ERR-019.
- **Why it hides:** All three homogeneous eigenvalue tests passed
  to machine precision because `δψ = 0` annihilates the weight-
  bearing term. 1G k = ν Σ_f/Σ_a, 2G/4G eigenvalue is still
  weight-independent for uniform medium. Heterogeneity is the
  only stressor.

---

## Signature 3: Scattering matrix transpose convention drift

- **Symptom:** Multi-group keff is catastrophically wrong (e.g. 2.06
  vs 1.0) but 1-group is exact. Or: 2-group keff is wrong but the
  ratio of group fluxes is the *swap* of the expected ratio.
- **Mechanism:** ORPHEUS convention is `Mixture.SigS[l][g_from, g_to]`
  — rows are source groups, columns are sink groups. The scattering
  source is `Q_g = Σ_g' SigS[g', g] · φ_g' = (SigS^T @ φ)_g`. A
  vectorised rewrite that uses `phi @ SigS^T` (instead of the
  algebraically correct `phi @ SigS`) silently double-transposes,
  which is invisible whenever `SigS = SigS^T` (1-group self-scatter,
  symmetric matrices).
- **Diagnostic probe:** L0-SN-009 — hand-calculate the scattering
  source for a 2-group asymmetric scattering matrix, compare to
  code output term-by-term.
- **Catching test:**
  - `tests/sn/test_quadrature.py::TestL0TermVerification::test_scattering_source_magnitude`
    (L0-SN-009 — hand calc against `SigS^T @ φ`, tagged
    `@pytest.mark.catches("ERR-002")`).
  - `tests/mc/test_properties.py::test_sigs_orientation_g0_to_g1`
    and `tests/mc/test_gaps.py` — explicit ERR-002-pattern guards
    in the MC verification suite.
- **Catalog entry:** ERR-002.
- **Why it hides:** Symmetric scattering matrices (any 1-group
  problem, any isotropic-scatter homogeneous problem) make
  `SigS = SigS^T` and the bug invisible. Asymmetric inputs are
  mandatory for the test to discriminate.

---

## Signature 4: Quadrature-dependent constant hardcoded

- **Symptom:** A solver works for one quadrature family (e.g.
  Lebedev) and diverges or oscillates for another (e.g. Gauss-
  Legendre). Streaming-equilibrium fixed-source test gives
  flux that does NOT equal `Q/Σ_t` under the failing quadrature.
- **Mechanism:** Different quadrature families have different weight
  sums:
  - Gauss-Legendre on `[−1, 1]`: `Σ w = 2`.
  - Lebedev / Level-Symmetric / Product (full sphere): `Σ w = 4π`.
  Any code path that bakes a literal `4π` (or `2`) into RHS
  normalization, source scaling, or angular-to-scalar conversion
  is implicitly assuming one family. The L0-SN-001 streaming
  equilibrium test (`φ → Q/Σ_t`) is the universal probe — it must
  hold for *every* quadrature.
- **Diagnostic probe:** Run streaming-equilibrium test
  (`Q` uniform, `Σ_t` uniform, vacuum BC) through the suspect
  solver path with each quadrature family in turn. Wrong
  normalization makes φ deviate from Q/Σ_t by a constant factor
  (the ratio of assumed weight sum to actual).
- **Catching test:**
  - `tests/sn/test_quadrature.py::TestWeightSums::test_gl_weights_sum_to_2`,
    `::test_lebedev_weights_sum_to_4pi`,
    `::test_level_symmetric_weights_sum_to_4pi`,
    `::test_product_weights_sum_to_4pi` (the weight-sum
    invariants — necessary preconditions for any L0-SN-001
    streaming-equilibrium check).
  - Plus a streaming-equilibrium test (`φ → Q/Σ_t`) run against
    the affected solver path under each quadrature family.
- **Catalog entry:** ERR-004 (BiCGSTAB hardcoded 4π — caught when
  GL quadrature was first tried). ERR-025 is a related case where
  the missing `1/W = 1/Σ w` normalization in the 1D cumprod
  recurrence cancelled exactly with a sign error in the numerator
  for GL — visible only at material interfaces.
- **Why it hides:** The bug is invisible whenever the quadrature
  in use happens to match the hardcoded constant. Initial
  development with one family gives full passing test suites;
  the bug surfaces when a second family is introduced.

---

## Signature 6: Log-singular kernel diagonal truncation

- **Symptom:** A semi-analytical reference solver built from a
  log-singular Fredholm kernel (Peierls / Bickley-Naylor / `E_1`
  family) stalls at 1–10 % accuracy under uniform Gauss-Legendre
  quadrature, instead of reaching the moment-floor accuracy of its
  upstream pieces. Refining the quadrature improves the answer
  *slowly* — but the convergence-rate fingerprint is **not** the
  classical `n^{-1/2}` Schneider endpoint signature.
- **Mechanism:** When the angular integral is pre-integrated
  analytically (as in F_N Path A.i), the residual kernel is
  `E_1(τ) = -γ_E - log(τ) + R(τ)` — log-singular at τ = 0. A
  μ-quadrature `Σ w_k/μ_k · exp(-τ/μ_k)` evaluates `E_1(τ)` to
  machine precision off-diagonal, but **silently saturates at finite
  `~ 2·log(n_μ)` at the diagonal** instead of diverging like
  `E_1(0+) = +∞`. The diagonal entries of the discretized kernel
  matrix are therefore **finite-but-wrong by an `n_μ`-dependent
  amount** that does not vanish with refinement. The system error
  signature is `||φ - φ_n||_∞ ~ log(n) / n` — first-order with a
  logarithmic correction, NOT `n^{-1/2}`.
- **Discriminator (the 2-D fingerprint):**
  - **`err · n / log(n)` is constant** (slow drift across 4× n range).
  - **`err · sqrt(n)` decays by ~3×** across the same range.
  - The first invariant is the textbook signature of log-singular
    kernel diagonal truncation. The second rules out Schneider
    endpoint singularity in the *solution*.
  - The discriminator is essential because both bug classes look
    similar at low-`n` (1–10 % error, slowly improving) — the
    rate fingerprint is the only cheap distinction.
- **Diagnostic probe:** Run the suspect reconstruction at
  n = 16, 32, 64, 128, 256, plot `||err||_∞ · n / log(n)` and
  `||err||_∞ · sqrt(n)`. The first ≈ constant + the second
  decaying → confirms Signature 6.
- **Catching test:**
  - `tests/derivations/test_path_ai_legacy_plain_gl_signature.py`
    (4 tests pinning the failure-mode signature: log-decomposition
    foundation, diagonal-truncation `2·log(n_μ)` scaling,
    off-diagonal machine-precision convergence, first-order-with-log
    rate). Tagged `@pytest.mark.catches("ERR-036")`.
  - `tests/derivations/test_atkinson_product_nystrom.py` — pinned
    fix verification (n_panels=64 → 3.5e-4; 256 → 1.1e-5).
- **Catalog entry:** ERR-036.
- **Why it hides:** Plain Gauss-Legendre integrates *smooth*
  functions with spectral accuracy. The diagonal saturation is a
  *one-point* defect that contributes O(`log(n)/n`) to the system
  norm — barely worse than first-order, easy to dismiss as
  "expected for a reasonable but uninformed quadrature." Off-
  diagonal machine-precision agreement at `τ > 0` further hides
  the bug by making the kernel matrix look correct elsewhere.
- **Fix family:** Atkinson product-Nyström. Decompose the kernel
  as `K(t,s) = L(t,s)·log|t-s| + smooth(t,s)`; integrate the
  log-singular subkernel via closed-form weights
  `F_k(s;t) = ∫sᵏ log|t-s| ds` for k = 0, 1, 2 (product-Simpson);
  Gauss-Legendre handles the smooth remainder. With graded mesh
  `t_j = (j/n)^q, q = 4`, recovers full O(h⁴ log h) convergence.
  See `docs/theory/references/fn_method.rst` and Atkinson 1976 / 1997.

---

## Signature 7: Quadrature endpoint pole-cancellation slow convergence

- **Symptom:** A semi-analytical evaluation involving an integrand
  with an algebraic pole (typically at `μ = 1`) algebraically
  cancelled by an opposing factor (typically `(c·tanh⁻¹μ)² ∼
  log²(1−μ)` growth) returns a 1–2 % wrong answer that **slowly
  improves with `mp.dps` but never reaches expected precision**.
  The error decay is monotone (3.3% at dps=15 → 1.6% at dps=35) —
  hallmark of a *quadrature* problem, not a *precision* problem.
- **Mechanism:** `mp.quad` (and `scipy.integrate.quad`) handle
  bounded integrands with smooth derivatives spectacularly well,
  but algebraic pole cancellation is *not* in the smooth class
  even when the integrand is mathematically bounded. The library
  cancels the pole *via local refinement around the singularity*,
  paying a polynomial cost in `dps` per digit recovered. The
  asymptote is the wrong value because the local refinement
  underestimates the contribution of the cancelled neighbourhood.
- **Discriminator (the convergence-with-dps fingerprint):**
  - Run the same solver at increasing `mp.dps` (15 → 25 → 35 → 50).
  - **Log-error decays monotonically** with dps.
  - **The decay rate is sub-linear** — adding 10 dps does NOT add
    10 digits of accuracy.
  - **The asymptote is finite, not zero** (extrapolated `lim_{dps→∞}
    err > 0`) — the smoking gun.
  - Contrast: a pure precision problem reaches machine-zero error
    at finite `dps`. A pure quadrature endpoint bug saturates at
    a non-zero asymptote.
- **Diagnostic probe:** Drop `dps ∈ {15, 25, 35, 50, 100}` into the
  same call; tabulate the log-error. If the slope of
  `log(err) vs dps` is shallow and the values plateau, you have
  Signature 7.
- **Catching test:**
  - `tests/derivations/test_case_method_z0.py::test_atalay_z0_table1_isotropic[*]`
    — 11 cases parametrised over c ∈ [0.1, 0.99]. Pre-fix: 1.5–2 %
    error at dps=15. Post-fix: 6–7 digits agreement.
    Tagged `@pytest.mark.catches("ERR-037")`.
- **Catalog entry:** ERR-037.
- **Why it hides:** "It converges with more precision" looks like
  the system is doing the right thing — it's just slow. The bug is
  invisible at low precision (looks like normal discretization
  error) and the user typically stops increasing dps once the
  *change* per step gets small (declaring "converged"), without
  comparing against an external reference. Without an external
  truth value, the wrong asymptote looks correct.
- **Fix family:** Substitution-based pole cancellation. For a pole
  at μ = 1 with `1/(1−μ²)` factor, substitute `μ = tanh(t)` mapping
  `(0, 1) → (0, ∞)` with Jacobian `sech²(t) = 1 − μ²` that cancels
  the pole *exactly* under change-of-variables. After substitution,
  the integrand is smooth at the new endpoint t → ∞ and the
  exponential decay makes the integral converge spectrally.
  Single-line code change; `mp.dps = 15` then suffices for 6–7
  digits.
- **Generalisation:** any algebraic-pole-cancellation integrand
  benefits from change-of-variables that builds the cancelling
  factor into the Jacobian. The integrand becomes "no pole" rather
  than "cancelled pole", and standard quadrature works.

---

## Signature 5: Curvilinear α-dome non-positivity

- **Symptom:** Curvilinear SN solver produces NaN after first
  iteration. Or: alpha values inspected on a level show negative
  entries.
- **Mechanism:** The α-recursion in curvilinear quadrature is a
  cumulative sum that should produce a non-negative dome on each
  level: α(0) = α(N) = 0 with positive interior. A sign error,
  wrong cosine choice (mu_y vs mu_x), or wrong cumulative direction
  produces negative α, which propagates as `1/α` or `√α` in the
  WDD closure and yields NaN.
- **Diagnostic probe:** `assert np.all(alpha >= -1e-14)` on every
  level of `SNMesh.alpha_per_level` before the first sweep.
- **Catching test:**
  `tests/sn/test_quadrature.py::TestAlphaRedistribution::test_alpha_dome_non_negative`,
  `::test_alpha_boundary_zero`, and
  `::test_spherical_alpha_dome_non_negative`.
- **Catalog entry:** Uncatalogued — pattern only. (The
  α-dome positivity check is a *necessary* invariant for ERR-006
  to be impossible; it does not have a standalone ERR entry.)
- **Why it hides:** It does not hide for long — NaN propagation is
  catastrophic. The signature exists to ensure the failure mode
  is caught at the *coefficient* level before it reaches the
  sweep, where the diagnostic message would be far less clear.

---

## Signature 8: Diverges-with-refinement masquerade (unconverged inner solve)

- **Symptom:** The observable (k_eff, flux) *diverges* (or drifts
  monotonically) as the mesh is refined `h → 0` (or as the angular
  order is raised), with the textbook fingerprint of an inconsistent
  discretization — error GROWING under refinement. But the
  fixed-source / homogeneous sanity cases are clean, and the per-
  ordinate flat-flux residual (Signature 1's probe) passes. The
  discretization looks correct yet the answer worsens with `nx`.
- **Mechanism:** The error is NOT in the discrete operator — it is an
  **unconverged inner solve** masquerading as a discretization
  inconsistency. A library convergence flag was discarded (scipy
  `gmres`/`bicgstab` return `(x, info)`; `info > 0` means "did not
  converge in `maxiter`" and is silently dropped), OR a hardcoded
  Krylov `restart` / `maxiter` cap truncates the inner iteration
  below convergence. Refining the mesh raises the inner system's
  condition number, so the *same* iteration cap delivers a
  *progressively less-converged* inner solve — the outer observable
  drifts because each refinement level is solved to a different
  (degrading) inner residual, not because the discrete operator is
  inconsistent. The "diverges with `h`" curve is an artefact of the
  cap, not the stencil.
- **Discriminator (rules it out vs Signature 1):**
  - **Tighten `inner_tol`: the fixed point does NOT move.** If the
    converged answer is insensitive to `inner_tol` across decades,
    the inner solve was already converged → it IS a discretization
    bug (Signature 1). If the answer *moves* when `inner_tol`
    tightens, the inner solve was the culprit (Signature 8).
  - **Sweep `restart` (and `maxiter`): the answer SNAPS to the exact
    value** once the cap clears the iteration count the system
    actually needs. A discretization bug ignores `restart` entirely;
    Signature 8 is gated by it.
  - **Capture and assert the discarded `info` flag.** `info > 0` on
    any refinement level is the smoking gun — the level "converged"
    to `maxiter`, not to tolerance.
- **Diagnostic probe:** At the finest failing mesh, (a) re-solve at
  `inner_tol ∈ {1e-6, 1e-9, 1e-12}` and tabulate the observable —
  drift across the column confirms Signature 8; (b) sweep
  `restart ∈ {30, 100, 300}` (and lift `maxiter`) and watch the
  observable snap to a stable value; (c) instrument the scipy call
  to `assert info == 0` per level.
- **Catching test:** speculative — no dedicated regression test yet.
  The structural fix is to **never discard the `info` flag**: route
  every scipy iterative-solver return through a checked wrapper that
  raises on `info != 0`, and pin it with a deliberately under-`maxiter`
  case that MUST raise. (When promoted, file the ERR-NNN and link the
  test path here.)
- **Catalog entry:** Uncatalogued — pattern only (numerics-investigator
  lessons L10/L9: "diverges with refinement + a discarded library
  info-flag = an unconverged inner solve, not a discretization bug";
  bound the solver cost before declaring a hang).
- **Why it hides:** It wears Signature 1's costume exactly — error
  growing under refinement is the canonical "inconsistent
  discretization" tell, so the investigator anchors on the stencil
  and never checks the inner residual. The fixed-source and per-
  ordinate probes pass (the operator IS consistent), which
  *reinforces* the wrong hypothesis ("the operator is fine, so the
  bug must be subtle"). The discarded `info` flag is invisible
  unless you instrument the solver call.

---

## Signature 9: ρ-blind stopping (increment understates the true error)

- **Symptom:** An iterative solver "reports converged" (its stopping
  test is satisfied, residual prints look small) but the answer is
  wrong by a factor that **GROWS as the scattering ratio `c → 1`**
  (or as the dominance ratio `ρ → 1` — near-critical, highly-
  scattering, optically-thick reflective). At `c = 0.5` the answer is
  fine; at `c = 0.99` it is visibly off; the error tracks `1/(1−c)`.
- **Mechanism:** The stopping test measures the **iterate increment**
  `‖Δψ‖ = ‖ψⁿ − ψⁿ⁻¹‖` and declares convergence when it drops below
  `tol`. But for a linear fixed-point iteration with spectral radius
  `ρ`, the true error relates to the increment by
  `‖ψⁿ − ψ*‖ ≈ ‖Δψ‖ / (1 − ρ)`. As `ρ → 1` (which happens precisely
  as `c → 1` for source iteration), the increment **understates** the
  true error by the factor `1/(1−ρ)`, which diverges. The solver
  stops with a small increment while still far from the fixed point —
  a *false* convergence, undetectable by the increment-based test.
- **Discriminator:**
  - **The error scales as `1/(1−c)` (or `1/(1−ρ)`)**, not as a fixed
    offset — sweep `c ∈ {0.5, 0.9, 0.99, 0.999}` and the relative
    error climbs in lock-step with `1/(1−c)`.
  - **Measure the RESIDUAL `r = Aψ − q`, not the increment.** A
    residual-based stop is `ρ`-honest: `‖r‖` does NOT carry the
    `1/(1−ρ)` amplification, so it reports the true distance to the
    fixed point. If switching the stop from `‖Δψ‖` to `‖r‖` moves the
    converged answer (and the move grows with `c`), Signature 9 is
    confirmed.
- **Diagnostic probe:** Re-run the near-critical case with the
  residual `r = Aψ − q` computed and asserted below `tol` (NOT the
  increment). Tabulate converged-`ψ` vs `c ∈ {0.9, 0.99, 0.999}`
  under both stops; the increment-stop answer drifts with `c`, the
  residual-stop answer is stable. (See the FluxDisplacement /
  AngularResidual diagnostic catalogue — the increment lives on the
  displacement, the residual is the rate-density `r = Aψ − q`.)
- **Catching test:** speculative — no dedicated regression test yet.
  The structural fix is `‖Aψ − q‖`-based stopping (or the a-posteriori
  correction `‖Δψ‖ / (1 − ρ̂)` with an estimated `ρ̂`); pin it with a
  high-`c` (`c ≥ 0.99`) fixed-source case whose exact answer is known
  (`φ = Q/Σ_a` for an infinite homogeneous absorber), asserting the
  converged value to tolerance. ORPHEUS precedent: #208's
  `FluxDisplacement` carries the `‖Δψ‖/(1−ρ)` a-posteriori true-error
  estimate for exactly this fix.
- **Catalog entry:** Uncatalogued — pattern only (numerics-investigator
  lesson L11: "measure the residual `r = Aψ − q` for a ρ-honest stop,
  not the increment `‖Δψ‖`").
- **Why it hides:** The stopping test is *self-consistent* — the
  increment really did get small, so every printed diagnostic looks
  healthy. The error is invisible at the moderate-`c` configs that
  dominate the test suite (`c ≤ 0.9` keeps `1/(1−ρ)` modest); only
  the near-critical / thick-reflective corner activates the
  amplification. Without an external truth value at high `c`, the
  false-converged answer looks plausible.

---

## Signature 10: Stale-snapshot huge-ULP red on ONE geometry (live-correct / frozen-stale)

- **Symptom:** A regression / golden-snapshot test reddens with an
  **enormous ULP distance** (thousands to millions of ULP, not the
  handful expected from an FP-reduction-tree change) on a **single
  geometry** (e.g. SPHERE), while the **sibling geometries** (slab,
  cylinder, Cartesian) using the same code path pass clean. The
  failure is localized to one `.npy` baseline, not spread across the
  suite.
- **Mechanism:** The PRODUCTION code is **live-correct**; the
  **SNAPSHOT is stale**. A legitimate, verified change to the
  production path (a new closure, a re-baselined reduction, a fixed
  bug) updated the live value, but the frozen `.npy` baseline for that
  one geometry was never regenerated — so the test compares today's
  correct output against yesterday's superseded one. The single-
  geometry locality is the tell: a real new bug in a shared code path
  would redden the siblings too; a stale baseline reddens only the one
  snapshot that was missed during re-baselining.
- **Discriminator:**
  - **Sibling geometries on the same code path PASS.** A genuine bug
    in shared logic cannot be that selective — if slab/cylinder are
    green and only sphere is red, suspect the baseline, not the code.
  - **Blob-hash the `.npy`** and `git log` the commit that last wrote
    it; check whether a later, verified production change (the
    "diverging apply-commit") post-dates the baseline. If the code
    moved after the snapshot froze, the snapshot is stale.
  - **The ULP magnitude is wrong for an FP-noise story** — a
    reduction-tree change drifts by `O(reduction-depth × ULP)` (tens
    to a few thousand ULP, bounded per `vv-principles` §bit-identity
    criterion 3). Millions of ULP means the *values genuinely differ*,
    which is a stale baseline (or a real bug — the sibling-pass test
    discriminates).
- **Diagnostic probe:** (a) Run the failing snapshot's siblings — if
  they pass, Signature 10 is likely. (b) Re-baseline the failing
  `.npy` against a **structurally-independent** reference (closed-form
  `k_∞ = νΣ_f/Σ_a`, MMS, or higher-precision recompute) — **NEVER**
  old-vs-new ULP, which only proves "differs from the stale value,"
  not "is correct." (c) If the independent reference agrees with the
  live value, the snapshot was stale → regenerate it.
- **Catching test:** speculative — this is a *triage* signature, not
  a single regression gate. The defense is process: re-baseline every
  geometry's snapshot in the same commit that changes the shared path,
  and pin the live value against a structurally-independent reference
  (not the prior snapshot) so a stale baseline cannot masquerade as a
  bug. (When a concrete instance is catalogued, file the ERR-NNN.)
- **Catalog entry:** Uncatalogued — pattern only (qa finding; cf.
  the SPHERE stale-snapshot cluster — multiple frozen baselines went
  main-baseline-red while production stayed live-correct).
- **Why it hides:** A red regression test reflexively reads as "you
  broke something" — the investigator chases the production code
  instead of the baseline. The huge ULP distance *amplifies* the
  alarm (looks catastrophic). The fix is to invert the default
  reaction for a single-geometry huge-ULP red: check the siblings and
  the baseline's provenance FIRST, and re-baseline against a
  structurally-independent reference — old-vs-new ULP can never tell
  you which side is wrong.

---

## Cross-cutting hygiene rules

These rules are invariants implied by the signatures above. The QA
agent enforces them at review time; numerics-investigator and
test-architect should use them as a fast filter when triaging
"all tests pass" claims.

### H1: 1-group eigenvalue tests are degenerate

`k = νΣ_f / Σ_a` regardless of flux shape. Every signature in this
catalog has at least one branch that hides under 1-group testing
(ERR-001, ERR-002, ERR-006, ERR-007, ERR-019, ERR-025, ERR-026 all
explicitly survived 1G suites). Any verification claim *must*
include a multi-group (≥ 2G) heterogeneous test.

### H2: Homogeneous eigenvalue is degenerate to redistribution

Flat ψ makes every redistribution / α-recursion / weight-cancel-
ling term identically zero. Tests that pass on homogeneous
material prove only that the *non-redistribution* parts of the
operator are correct. Multi-region or non-uniform Σ_t is
mandatory.

### H3: Conservation is necessary, never sufficient

Total particle balance and global conservation are *telescoping*
sums — they hold by construction even when per-ordinate or per-
group balance is wrong. The per-ordinate flat-flux residual
(Signature 1's diagnostic probe) is the canonical example: ERR-006
satisfied global conservation to machine precision while having
20–50 % per-cell flux error.

### H4: Convergence rate is not convergence value

A solver can show clean O(h) or O(h²) convergence to the *wrong*
asymptote when the reference is built from the same buggy code
(self-referencing Richardson). Convergence-rate tests must use
an external reference (analytical, MMS, or independently-verified
solver) — see the T3 dead-end pattern in
`docs/theory/methods/diffusion_1d.rst`.

### H5: Test count is not coverage

Twenty passing tests do not bound coverage of the failure modes
above. Demand a *heterogeneous, multi-group, mesh-refinement
convergence test* before accepting any "all tests pass" claim.

---

## Add-a-signature protocol

To add a new signature to this catalog:

1. **Confirm the bug class is recurrent.** A one-off bug goes in
   `docs/theory/verification/error_catalog.rst` only. A signature here must be a
   *class* of plausible-wrong errors that the test authors have
   already missed at least once.
2. **Log the underlying ERR-NNN.** If the bug instance does not
   yet have a catalog entry, add one to
   `docs/theory/verification/error_catalog.rst` first using the standard template
   (failure mode, date, solver, bug, impact, how it hid, L0 test,
   lesson).
3. **Identify the catching test by exact pytest path.** No
   "regression test exists" gestures — give
   `tests/<file>.py::<class>::<function>` so the signature can be
   re-run on demand. If no test catches it, the signature is
   speculative — tag it as such, do *not* fabricate a test path.
4. **Name the test classes blind to it.** Every signature must
   list explicitly which test classes (1G, homogeneous,
   conservation, single-iteration, etc.) are blind to the bug.
   This is what makes it a *signature* rather than a defect log.
5. **Write the diagnostic probe as the cheapest test.** The probe
   should be runnable in seconds, not minutes. If the cheapest
   probe is itself a slow integration test, factor it down.
6. **Cross-link.** Add the ERR-NNN to the signature; add the
   signature to the ERR-NNN entry's "Lesson" section as
   "→ numerical-bug-signatures Signature N."

A signature whose ERR-NNN entry is later removed should be
removed from this catalog at the same time. The catalog is
canon; signatures here are a derived view.

---

## Reference: ERR-NNN cross-index

| Signature | Catalog ERR    | Failure mode | Catching test (primary)                                                                              |
| --------- | -------------- | ------------ | ---------------------------------------------------------------------------------------------------- |
| 1         | ERR-006        | #2 + #3      | `tests/sn/test_quadrature.py::TestL0TermVerification::test_per_ordinate_flat_flux_consistency`       |
| 1         | ERR-007        | #3           | same — BiCGSTAB-operator variant covered by `@catches("ERR-006","ERR-007")`                          |
| 1         | ERR-026        | #6           | `tests/sn/test_sweep_operator_inconsistency.py::test_spherical_sweep_vs_bicgstab_flat_flux`          |
| 2         | ERR-019        | #3           | `tests/moc/test_verification.py::TestL0EquilibriumFlux::test_pure_scatterer_equilibrium_single_sweep`|
| 3         | ERR-002        | #2           | `tests/sn/test_quadrature.py::TestL0TermVerification::test_scattering_source_magnitude`              |
| 4         | ERR-004        | #4           | `tests/sn/test_quadrature.py::TestWeightSums::test_gl_weights_sum_to_2` + streaming-equilibrium      |
| 4         | ERR-025        | #3 + #4      | `tests/sn/test_cartesian.py::test_heterogeneous_absolute_keff`                                       |
| 5         | (uncatalogued) | —            | `tests/sn/test_quadrature.py::TestAlphaRedistribution::test_alpha_dome_non_negative`                 |
| 6         | ERR-036        | #3 + #4      | `tests/derivations/test_path_ai_legacy_plain_gl_signature.py` + `test_atkinson_product_nystrom.py`   |
| 7         | ERR-037        | #4           | `tests/derivations/test_case_method_z0.py::test_atalay_z0_table1_isotropic`                          |
| 8         | (uncatalogued) | —            | speculative — checked-`info` wrapper raising on `info != 0` + under-`maxiter` case that MUST raise   |
| 9         | (uncatalogued) | —            | speculative — `‖Aψ − q‖`-based stop pinned by high-`c` (`c ≥ 0.99`) fixed-source vs `φ = Q/Σ_a`      |
| 10        | (uncatalogued) | —            | speculative (triage) — sibling-pass discriminator + re-baseline vs structurally-independent reference|
