# Platform-independent quadrature and k∞: clearing the macOS-update reds on `main`

## Context

On 2026-10-01 `tests/gates -m "not slow"` failed 26 tests on a clean `main` (`0a5a23fa`). The same suite passed 13 159 tests at `8f9300b2` on 2026-09-30, and no Python has changed since.

The numerics-investigator localised all 26 failures to the system linear-algebra library, Accelerate, which changed when macOS 27.0.1 was installed on 2026-09-30 at 21:50 (memo: `scratch/platform_drift/memo.md`). No production code is wrong.
- **24 failures come from the Gauss rules.** The rules are built by Golub–Welsch `np.linalg.eigh` in `GeneratingMeasure.gauss` (`orpheus/numerics/generating_measure.py:329`), and they moved by 1–2 ULP.
- **2 failures come from the homogeneous k∞.** It is computed with `np.linalg.eig` (`orpheus/numerics/eigenvalue.py:577`), which moved by 1 ULP.
- **7 of the 26 are test-design defects:** their verdict is rounding noise, and a correct input fails them in 3 to 8 of 10 random ±3-ULP jitters.

**The user's rulings (2026-10-01):**
- **Order:** diagnose first, then remedy; the reflective cleanup waits until this lands.
- **Correctly rounded Gauss rules:** each rule is computed once in high precision and rounded once, so every platform gets the same rule. The pins downstream are then re-captured once.
- **k∞ from the rank-one fission operator:** k∞ = ⟨νΣf, A⁻¹χ⟩, with no dense eigen-solve.
- **The homogeneous byte pin is re-posed:** an exact rational reference for the float inputs, with a ULP bound derived from the LU error bound, never a hand-picked number.

Mode: W2 remedy on a branch `fix/platform-independent-quadrature`. The main agent writes the code. The test-architect specifies the gates first, qa and the elegance-enforcer review, and the archivist writes the docs.

The parked reflective cleanup is untouched. Its two gate files are in `scratch/boundary_ontology/reflective_gates/pending/`, and it resumes on a fresh branch after this merges.

## The measured population (explorers, 2026-10-01)

- **Gauss map** (`scratch/platform_drift/gauss_map.md`):
  - **Families:** one generic `gauss` body serves five constructors: `LEGENDRE`, `CHEBYSHEV_T`/`U`, `HERMITE`, `jacobi(a,b)` and `laguerre(a)`. Every recurrence is closed form: rational in (k, a, b), plus one transcendental zeroth moment.
  - **Who reaches it:** production reaches only `LEGENDRE`, through `gauss_legendre_on_mu`.
  - **mpmath:** 1.3.0, present only through sympy and not declared in `pyproject.toml`. `orpheus/numerics` does not import it yet.
  - **Calls:** a spy over 2 269 tests counted 1 715 `gauss` calls covering 116 distinct (family, n) pairs. Nothing caches them, and the returned arrays are writeable.
  - **Effect of correctly rounded Legendre alone:** it fixes 6 of the reds and turns 4 green pins red, leaving 26 rule-dependent reds to re-capture or redesign. About 100 frozen artefacts in 15 groups are fed by a GL rule, produced by 7 generator scripts plus the `--capture-baseline` flag in 5 test modules.
- **Homogeneous map** (the explorer's report; plan mode blocked its file, so its content goes to `scratch/platform_drift/homogeneous_map.md` at step 0):
  - **F is rank one** in 8 of 8 shipped cases, and in every posable problem: `Mixture.chi` is a single vector, and `FissionKernel` requires 1-D factors.
  - **The byte pin** (`tests/gates/homogeneous/test_byte_stability.py`, fixture `_fixtures/cs1_prewiring.json`) compares, for 8 mixtures, k as hex, the flux as raw bytes, and `sig_prod`/`sig_abs` as hex.
  - **The rank-one route moves k** in 5 of 8 cases by 1–2 ULP. The stored record is correctly rounded for the decimal inputs, not the float ones.

## Step 0 — record the maps (no code)

- Save the homogeneous map to `scratch/platform_drift/homogeneous_map.md`.
- Create the branch.
- File two issues, each with `module:`/`level:`/`type:` labels:
  1. **The rank-one χ in the data layer.** `production_weighted_chi` averages per-isotope χ at flat flux. The exact operator is X Nᵀ, and its k is λ_max(Nᵀ A⁻¹ X), a K×K problem.
  2. **The remaining platform seams.** The level-symmetric weights use `np.linalg.solve` (`rules_sphere.py:713`), the azimuthal nodes use the system `np.cos`/`np.sin` (`roots_of_unity.py:246`), and `derivations/` holds 56 numpy `leggauss` calls. None drifted this time; all are the same class of exposure.

## Step 1 — the gates first (test-architect)

The test-architect specifies these gates before any code. Each new gate carries its first red, measured on HEAD.

- **(a) Correctly rounded rules.** For every family and a range of n, the rule equals an independent mpmath computation at a higher precision, rounded once, bit for bit:
  - Legendre, against `mpmath.gauss_quadrature`;
  - a Jacobi pair and Laguerre, against their mpmath equivalents.
  - Symmetric families come out exactly symmetric.
  - First red: today's Golub–Welsch rules are off by up to 445 ULP at n = 64.
- **(b) A platform-independence witness.** A byte fingerprint (a hash) of the GL rules for n ∈ {2, 4, 8, 16, 64}, frozen. It is a producer pin placed in front of every snapshot downstream (the `vv-principles` producer-pin rule), so the next library drift reds this one gate with a clear name rather than 26 snapshot gates.
- **(c) The k∞ reference.** The exact rational k and flux of the float inputs, computed by `Fraction` elimination in a reference module under `orpheus/derivations/` (the L0 home of exact references), and asserted against production within a derived ULP bound: the forward-error bound of an LU solve, a backward error of order G·u times the condition number κ(A), evaluated per case. The bound's derivation is written into the test's docstring.
  - Mutation: perturbing χ by 1 ULP reds it.
  - This gate replaces the byte comparison of k and flux in `test_byte_stability`. The `sig_prod`/`sig_abs` rows follow the same contract.
- **(d) Redesign the noise gates,** whatever their fate under correctly rounded rules:
  - `test_diamond::TestBitIdenticalCurvilinear`
  - the SPECULAR "blind case"
  - `test_angular_bulk_space` G1.3
  - the two MMS restriction legs
  - `psi_half[64]`
  - `test_boundary`'s GL-4 table (rtol 1e-15)
  - `test_legendre_basis` (37 ULP against a 32-ULP tolerance under correctly rounded rules)

  For each, the gate pins the invariant it was written for (association-independent or with a derived bound) rather than an accident of rounding. No tolerance is widened: where a bound changes, it is derived and its derivation is stated.
- **(e) The call-path spy** (`test_homogeneous.py:231`) and `test_R5` (exact equality on A_2g, 2 ULP off on the rank-one route) are re-posed onto the rank-one route.

## Step 2 — correctly rounded Gauss rules (commit 1)

- **`GeneratingMeasure.gauss`** computes the nodes and weights in mpmath and rounds each once to float64.
  - **The recurrence:** each family's recurrence coefficients are evaluated in mpmath from the same closed form (one definition). The float `recurrence(n)` is the rounded image of the mpmath one, so the two cannot drift. Its other consumers are checked by the census at implementation time.
  - **The near-tie check:** the rule is computed at two working precisions (about 40 and 60 digits), and the two rounded results must agree bit for bit, or the construction raises.
  - **Symmetry:** the mirror averaging becomes an assertion. A correctly rounded rule for an even weight is exactly symmetric, because rounding commutes with negation.
  - **Zeroth moment:** the renormalisation retires. It would un-round correctly rounded weights. `[M]` the weights sum to exactly 2.0 for n ∈ {4, …, 128}; gate (a) records the sum per family.
- **Cache:** `functools.lru_cache` over (family, n) stores the rounded `(nodes, weights)` tuples, and each call builds a fresh `DiscreteMeasure`, so no caller can mutate shared arrays. There are 116 distinct pairs per suite at about 0.08 s each.
- **The dependency:** `mpmath` is declared in `pyproject.toml`. L1 (`numerics/`) is mathematics only, and a high-precision arithmetic library is admissible there. The layer gate has no third-party rule to change.
- **Docstrings:** `gauss_legendre_on_mu` drops its "1–4 ULP vs leggauss" claim. The `generating-measure-golub-welsch` text states that Golub–Welsch is the algorithm and the result is correctly rounded. `jacobi(0,0)` against Legendre stays bit-identical, because both now take one correctly rounded path; gate (a) checks it.
- **Unchanged:** `gauss_chebyshev`, which only tests call, rides the same body.

## Step 3 — the homogeneous k∞ from the rank-one fission (commit 2)

- **`solve_homogeneous_infinite`** (`orpheus/homogeneous/solver.py:389-467`) spells k and the flux from the production operator's two factors, χ and νΣf. `IsotropicFission` is already a `RankOneOperator`, so the code reads as the mathematics:
  - u = A⁻¹χ, one solve through the existing `MatrixInverseOperator`;
  - k = ⟨νΣf, u⟩;
  - φ = 100·u/k, the existing gauge.

  `dominant_eigenpair` leaves this path. It stays in `numerics/` for its other callers.
- A non-rank-one F is not spellable here, because `FissionKernel` requires 1-D factors. The issue from step 0 records that the data layer is where rank one is decided.
- **Theory page** `docs/theory/foundations/infinite_medium.rst` (the archivist):
  - :1738–1748 argues for the dense inverse. It is rewritten to the rank-one derivation, with the proof that the dominant eigenvalue of A⁻¹χνΣfᵀ is ⟨νΣf, A⁻¹χ⟩, from the rank-one structure and the positive cone.
  - :787 and :1759–1760 are already false and are fixed.
  - The `.. implements::` blocks at :699 and :1366 are re-pointed.

## Step 4 — the ruled re-capture (commit 3)

- **Which artefacts:** only those whose gates are red after steps 2–3 are re-captured, each by its own generator (the 7 scripts and the `--capture-baseline` flag; list in `gauss_map.md` §5). A GL-fed pin that stays green within its stated contract is left as it is.
- **Positive control first, per generator (X1):** run the generator with the recovered pre-update GL-4 rule (`scratch/platform_drift/pd_plugin.py`, `PD_MODE=old4`) and show that it reproduces the stored artefact. A generator that cannot reproduce its own artefact is not the instrument, and that is reported, never re-captured over.
- **The `vv-principles` bit-identity conditions, each stated in the commit body:**
  1. Principled: the correctly rounded rule is a named, platform-independent quantity.
  2. Verified against an independent reference: gate (a), mpmath at a higher precision, plus the published GL-4 and GL-8 tables.
  3. The drift is explained: for each artefact, the maximum ULP and relative change, traced to the rule's change.
- **The homogeneous fixture** `cs1_prewiring.json` is not re-captured: its k and flux rows are re-posed by gate (c).

## Step 5 — docs, review, merge

- **The archivist:**
  - the homogeneous theory page (step 3);
  - the quadrature theory page on the correctly rounded construction, the two-precision check and the producer fingerprint;
  - an evidence entry: a bit pin on platform LAPACK output is pinning the platform; the instance with its dates.

  An ERR entry is not owed, since no production defect was found; the seven noise gates are test-design defects and go on the V&V anti-pattern evidence page.
- **Review:** qa and the elegance-enforcer, in parallel. qa adds a mutation that re-introduces Golub–Welsch output; gate (b) and gate (a) must red.
- **Checks:** `dead_references` reads 0; Sphinx `-W` is clean.
- **Merge:** after the full `tests/gates -m "not slow"` run, merge `--ff-only`, push, and watch CI.
- **Agent memory:** the pending memory files from today's agents are committed at close-out, uncurated.

## Verification

- **Gates:** (a)–(e) each have a first red on HEAD and are green after the change.
- **Full suite:** `.venv/bin/python -O -m pytest tests/gates -m "not slow"`, serial, run from a detached worktree driver, so no tracked file is edited while it runs.
  - Target: 0 failures in the main tree, against the baseline of 26 (the 4 `test_write_guards` failures are an artefact of the worktree path).
  - Every changed pin is accounted for in the step-4 table.
- **A second platform:** CI runs on Linux with OpenBLAS. Gate (b)'s fingerprint must be green there too, which is the cross-platform witness the rulings promise.
- **Performance:** the suite's wall time is compared with the baseline (38 min); the cache must keep the mpmath cost under about 30 s.

Sizing `[R]`: about one session; the re-capture and the noise-gate redesign are the bulk.
