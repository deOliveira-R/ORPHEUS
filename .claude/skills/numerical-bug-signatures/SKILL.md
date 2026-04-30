---
name: numerical-bug-signatures
description: 'Use when triaging a "wrong answer" or "all tests pass but..." report in a numerical solver. Provides recognition signatures for known bug classes (cylindrical-DD divergence under refinement, MOC weight cancellation, scattering matrix transpose convention, quadrature weight sums, alpha-dome positivity), each linked to the test that catches it and the ERR-NNN entry in `tests/l0_error_catalog.md`. Examples: "k_eff diverges as I refine", "homogeneous case passes but multi-region fails", "transport sweep gives wrong answer in curvilinear geometry only". Preloaded by qa, numerics-investigator, and test-architect.'
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
its ERR-NNN entry in `tests/l0_error_catalog.md`. Use it to
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
`docs/theory/diffusion_1d.rst`.

### H5: Test count is not coverage

Twenty passing tests do not bound coverage of the failure modes
above. Demand a *heterogeneous, multi-group, mesh-refinement
convergence test* before accepting any "all tests pass" claim.

---

## Add-a-signature protocol

To add a new signature to this catalog:

1. **Confirm the bug class is recurrent.** A one-off bug goes in
   `tests/l0_error_catalog.md` only. A signature here must be a
   *class* of plausible-wrong errors that the test authors have
   already missed at least once.
2. **Log the underlying ERR-NNN.** If the bug instance does not
   yet have a catalog entry, add one to
   `tests/l0_error_catalog.md` first using the standard template
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
