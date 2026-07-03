# Verification Plan — `outer(Vector, Functional)` dyad carve + `ReactionRateFunctional`

Surgical operator-algebra carve crossing `numerics` ↔ `transport`,
refactoring the 0-ULP-canary-pinned `RankOneOperator` primitive. PLAN
MODE — this file is the test SPECIFICATION; no files under `tests/`
are written here.

## 0. Grounding findings (read the real code first)

These facts SHAPE the plan; several contradict the brief's framing and
must be honored:

1. **`RankOneOperator`'s SOLE production consumer is
   `FissionOperator.kernel`** (`operator.py:1694`, consumed at
   `fission.py:306`). Confirmed by `grep -rn RankOneOperator orpheus/`:
   every other hit is a docstring. **The brief's "scattering path that
   touches `RankOneOperator`" DOES NOT EXIST** — `ScatteringOperator.kernel`
   is `R∘Λ∘M` (an `OperatorProduct`, `scattering.py:649`), which never
   instantiates `RankOneOperator`. **Consequence for Part 1:** there is
   NO second production canary to preserve. The "dual canary" is
   actually (a) the fission crosscheck `test_fission_kernel_crosscheck.py`
   B.1/B.2 (0-ULP), and (b) the L0 primitive suite
   `test_tensor_product_operator.py::TestRankOneOperator` (the dyad's own
   unit gate). The plan treats THESE TWO as the dual canary and explicitly
   records that the scattering canary is vacuous.

2. **The current `RankOneOperator.apply` body** (`operator.py:1809-1812`):
   ```python
   inner = (self.right * x_arr).sum(axis=self.axis, keepdims=True)
   return self.left * inner
   ```
   The `inner` line IS `Functional.evaluate`. `left=χ`, `right=νΣf`,
   `axis=0`. The carve makes `right` a typed `Functional` whose
   `evaluate(x)` reproduces that exact reduction.

3. **`ProductionRateFunctional`** (`production_rate_functional.py`)
   ALREADY exists and ALREADY does `(νΣf · φ).sum(axis=0, keepdims=True)`
   (line 151) — byte-identical to `inner`. The carve GENERALIZES it to
   `ReactionRateFunctional(Σx)` (production = `Σx=νΣf`, absorption =
   `Σx=Σa`) and RETIRES `ProductionRateFunctional`. `FissionOperator.production_rate`
   (`fission.py:311`) currently returns `ProductionRateFunctional(...)`;
   it rewires to `ReactionRateFunctional(...)`.

4. **`Functional` Protocol** (`numerics/functional.py`) is the structural
   base: `evaluate(x) -> R`, contravariant input, covariant result,
   disjoint from `LinearOperator`. `ReactionRateFunctional` conforms
   structurally (no inheritance needed). `IntegralKernelOperator` Protocol
   (`integral_kernel_operator.py`) requires `FissionOperator` to keep
   exposing both `kernel` AND `production_rate`.

5. **Closed-form reference EXISTS**: `kinf_and_spectrum_homogeneous`
   (`derivations/common/eigenvalue.py:48`) returns `(k∞, φ*)` for any
   group count via `λ_max(A⁻¹F)`, `A = diag(Σt) − SigSᵀ`, `F = χ⊗νΣf`.
   This is the structurally-independent ground (an `np.linalg.eig` of the
   transfer matrix — NOT another ORPHEUS transport solver, NOT the rank-1
   path).

6. **EXACT closed-form values** (computed, not from memory):

   | mixture | k∞ | flux spectrum φ* (ℓ²-normed) | ⟨νΣf,φ*⟩ | ⟨Σa,φ*⟩ |
   |---|---|---|---|---|
   | A 2g | 1.8750000000 | `[0.70710678, 0.70710678]` | 0.15909903 | 0.08485281 |
   | A 4g | 1.4877619048 | `[0.6399506, 0.4571076, 0.2468381, 0.5662039]` | 0.19041882 | 0.12799012 |

   At φ*, `k∞ = ⟨νΣf,φ*⟩/⟨Σa,φ*⟩` EXACTLY (verified, `np.isclose` True
   both rows). **CONFIG-BLINDNESS TRAP (critical):** mixture A **2g** has
   a coincidentally FLAT spectrum (`[0.707,0.707]`), so the *flat-flux*
   ratio `Σνσf/ΣΣa = 1.875` equals k∞ there — a flat 2g test would be
   blind to flux-shape errors. Mixture A **4g** has a genuinely non-flat
   spectrum, and the flat ratio (1.4909) DIFFERS from true k∞ (1.4878).
   **Therefore the independent-pinning test MUST use the converged φ*,
   NOT a flat φ, and MUST include the 4g case** so a flux-shape-blind
   error is caught. (Mixture A 2g `SigP=[0.025,0.2]`, `χ=[1,0]`,
   `SigT=[0.5,1.0]`, `Siga=[0.02,0.1]`, asymmetric `SigS`.)

7. **`Mixture` accessors** (real names): `SigP` (= νΣf production),
   `absorption_xs` (= Σa), `SigT`, `SigS[0]` (P0, `[from,to]`), `chi`.

---

## 1.5 Claim-layer + pillar gate (vv-principles §1.5, done BEFORE the matrix)

| Part | Claim layer | Pillar | Structural ground | OK? |
|---|---|---|---|---|
| 1 (bit-identity) | equivalence/de-risk (NOT correctness) | bit-identity inheritance | the PRE-carve verified value | ✓ (0-ULP `array_equal`) |
| 2 (k∞ upgrade) | **eigenvalue** | **closed-form** (`λ_max(A⁻¹F)`) | `np.linalg.eig` of transfer matrix — independent of rank-1 path | ✓ (closed-form proves eigenvalue) |
| 3 (dyad law) | flux-shape / algebraic | closed-form (hand outer product) | hand-computed `v·(w·x)` | ✓ |
| 4 (mutations) | gate-teeth proof | n/a (mutation) | the gate's own RED | ✓ |

Part 2 pairs an **eigenvalue claim** with a **closed-form** pillar —
legal (MMS would NOT be). No row pairs eigenvalue with MMS. All chains
terminate in `np.linalg.eig`, structurally independent of the SUT.

---

## PART 1 — Dual-canary bit-identity preservation (Mode-11 sentinel)

**Goal:** prove the dyad refactor is 0-ULP-identical on every path that
touches `RankOneOperator`, AND that each canary's execution path
ACTUALLY enters the new functional-as-row-factor code (a fired sentinel,
not a green twin).

### 1.A The two real canaries (the brief's "dual canary", corrected)

- **Canary α — fission crosscheck** (`tests/sn/operators/test_fission_kernel_crosscheck.py`):
  - B.1 `test_apply_matches_hand_derived_emission` — `F.apply ≡`
    `hand_derived_fission_emission` (structurally-independent Python
    double-loop, `rtol=1e-13`). This is the CORRECTNESS anchor and is
    pillar-independent of the rank-1 reduction. MUST stay green unchanged.
  - B.2 `test_chi_times_production_rate_is_fission_apply` — `χ ·
    production_rate.evaluate(φ) == F.apply(φ)` at 0 ULP. This is the
    procedural-twin equivalence row the brief flags. **After the carve it
    becomes `χ · ReactionRateFunctional(νΣf).evaluate(φ) == F.apply(φ)`**
    — still 0-ULP by construction. Keep as an EQUIVALENCE row (not
    correctness; B.1 carries correctness). See Part 5 for the
    twin→upgrade migration.
- **Canary β — the L0 primitive suite**
  (`tests/numerics/test_tensor_product_operator.py::TestRankOneOperator`,
  9 rows incl. `test_apply_matches_hand_computed_outer_product`,
  `test_apply_multi_dim_spatial_parametrisation` at `nulp=4`,
  `test_compose_with_identity_via_tensor_product` 0-ULP). These pin the
  dyad's own apply against a hand outer product. MUST stay green; if the
  carve changes the CONSTRUCTOR signature (`RankOneOperator(left, right,
  axis)` → `outer(reconstruction, functional)`) these rows REWIRE to the
  new factory (Part 5).
- **Canary γ — scattering: VACUOUS.** Recorded explicitly:
  `ScatteringOperator` does NOT touch `RankOneOperator`. No scattering
  bit-identity gate is needed or possible for this primitive. (If the
  implementer's plan claims a scattering canary, REJECT it as a phantom —
  Mode-11: a gate whose call graph never reaches the changed line.)

### 1.B Mode-11 sentinel strategy (prove the canary EXECUTES the new code)

The carve introduces a NEW private reader: the dyad's `apply` now calls
`functional.evaluate(x)` as the row-factor (instead of the inline
`(right*x).sum(...)`). Per vv-principles Mode-11 "Sharpening (NEW private
adapter/reader)", the gold-standard proof is a **pytest in-process
wrap-sentinel**, not a file-write:

**Sentinel test** (`tests/sn/operators/test_fission_kernel_crosscheck.py`,
new row OR a new `tests/numerics/test_rank_one_functional_sentinel.py`):

```
def test_fission_apply_routes_through_functional_evaluate(monkeypatch, solver_4g):
    """Mode-11: F.apply MUST enter ReactionRateFunctional.evaluate.

    Wrap the functional's evaluate in-process; assert the counter > 0
    after F.apply. A green twin that routed around the new reader
    (inline reduction) leaves the counter at 0 and reddens this gate.
    """
    calls = []
    real_evaluate = ReactionRateFunctional.evaluate
    def counting_evaluate(self, x, /):
        calls.append(1)
        return real_evaluate(self, x)
    monkeypatch.setattr(ReactionRateFunctional, "evaluate", counting_evaluate)

    op = solver_4g.fission_op
    phi = _asymmetric_phi(...)
    op.apply(phi)            # the production matvec arm
    require(len(calls) > 0,
            "F.apply did not route through ReactionRateFunctional.evaluate "
            "— the dyad's row-factor is not the new functional (Mode-11: "
            "green twin routing around the rewired reader).")
```

- This is STRICTLY STRONGER than a file probe: it runs in-process on the
  SAME class the production path instantiates, so a routed-around path
  (the implementer left the inline `inner` line) cannot fake the wrap.
- **Decision point for the implementer (flag to user):** does
  `RankOneOperator.apply` call `functional.evaluate`, or does
  `FissionOperator.kernel` now build `outer(χ, ReactionRateFunctional(νΣf))`
  where `outer` is a free function returning a rank-1 `LinearOperator`
  whose `apply` calls `functional.evaluate`? The sentinel target
  (`ReactionRateFunctional.evaluate` vs `RankOneOperator.apply` calling a
  `.evaluate`) depends on this. The sentinel above assumes the row-factor
  IS the functional and its `evaluate` is on the apply call-path. **If the
  fused `kernel` keeps the raw `RankOneOperator(χ, νΣf)` and the functional
  is only a SEMANTIC reading (like today's `production_rate`), then the
  matvec does NOT route through `evaluate` and the sentinel would be a
  FALSE requirement** — exactly the trap the brief #1 of the prompt
  describes. The plan REQUIRES the implementer to declare which, and the
  sentinel is written to match. (If kernel stays raw rank-1, Part 1's
  Mode-11 sentinel instead targets the dyad: prove `kernel.apply` enters
  the `outer`-built rank-1 `apply`, and Part 2 carries the real coverage.)

### 1.C Bit-identity assertion (0 ULP)

For BOTH canaries, the assertion is `np.testing.assert_array_equal`
(0 ULP) for the `apply ≡ apply` / `χ·evaluate ≡ apply` rows, and
`assert_array_almost_equal_nulp(nulp=4)` ONLY where einsum-vs-`.sum`
reduction-tree drift is already accepted (the existing
`test_apply_multi_dim_spatial_parametrisation`). The carve must NOT
introduce a new reduction tree: `ReactionRateFunctional.evaluate` uses
the identical `(Σx · x).sum(axis=0, keepdims=True)` primitive as the
inline `inner`, so 0-ULP holds by construction. **If the implementer
vectorizes differently and ULP drift appears, REJECT** unless the three
vv-principles bit-identity criteria are met and the contract is narrowed
with a documented justification (drift ≤ reduction-depth × ULP).

---

## PART 2 — The structural-independence UPGRADE (headline)

**Goal:** REPLACE the retired procedural twin (B.2's `χ·evaluate ≡ apply`,
same numpy primitive both sides → vv L11 weak) with a
structurally-independent oracle: each `ReactionRateFunctional` pinned
against the closed-form `k∞` decomposition, INDEPENDENTLY.

**File target:** `tests/transport/test_reaction_rate_functional.py` (new;
the generalization of `test_production_rate_functional.py`, which retires
— Part 5). `@pytest.mark.foundation` (software-invariant + L0 value
correctness; no theory `:label:`). Could ALSO carry an `l1` row for the
eigenvalue claim — see 2.D.

### 2.A The ansatz and the closed-form expected values

Homogeneous infinite medium, reflective (∞-medium) → flux is the
dominant eigenvector φ* of `A⁻¹F`. Build φ* from
`kinf_and_spectrum_homogeneous(SigT, SigS[0], SigP, chi)`. The
fiberwise-over-space functional, summed over the (degenerate, 1-cell or
uniform-N-cell) spatial axis, reduces to the dot product `⟨Σx, φ*⟩`.

Two `ReactionRateFunctional`s, pinned INDEPENDENTLY (so a shared-factor
error the RATIO masks is caught):

```
ReactionRateFunctional(νΣf).evaluate(φ*).sum()  ==  PROD_EXPECTED   (independent)
ReactionRateFunctional(Σa ).evaluate(φ*).sum()  ==  ABSN_EXPECTED   (independent)
k∞_recovered = PROD/ABSN  ==  kinf_closed_form                       (composed)
```

| mixture | PROD_EXPECTED ⟨νΣf,φ*⟩ | ABSN_EXPECTED ⟨Σa,φ*⟩ | k∞ |
|---|---|---|---|
| A 2g | 0.15909903 (≈0.025·0.707+0.2·0.707) | 0.08485281 | 1.8750000000 |
| A 4g | 0.19041882 | 0.12799012 | 1.4877619048 |

**Why pin EACH independently (not just the ratio):** a bug that scales
BOTH `νΣf` and `Σa` by the same factor (e.g. a wrong volume measure, a
shared `mat_xs` accessor returning a mis-scaled field) leaves the RATIO
`k∞` correct but moves PROD and ABSN. The independent pins catch it; the
ratio-only pin is blind. This is the headline upgrade's teeth: B.2's old
twin only ever checked `χ·evaluate ≡ apply` (same primitive both sides),
which is invariant to ANY shared error in the functional itself.

**Construction of the functional reference (structural independence):**
`ReactionRateFunctional` is built directly with a `CrossSectionField`
wrapping `Σx`; φ* comes from `np.linalg.eig`. The expected scalars
`PROD_EXPECTED` etc. are computed from `νΣf @ φ*` IN THE TEST via the
hand dot product — NOT via the functional under test. The functional's
`.evaluate(φ*).sum()` is the SUT output; the `νΣf @ φ*` is the reference.
These share NO reduction primitive at the spatial-fiber level if the test
writes the reference as an explicit `float(sum(nusf[g]*phi[g] for g in
range(ng)))` Python loop (mirroring `hand_derived_fission_emission`'s
discipline). **MUST write the reference as a Python-loop dot, not
`νΣf @ φ*`** (a `@` would share numpy's reduction with the functional's
`.sum` — procedural-twin risk recurs).

### 2.B Anti-pattern honoring (config-blindness declarations)

- **≥2 groups (Cardinal Rule / anti-#3):** BOTH cases are ≥2G. The 1g
  mixtures are FORBIDDEN for the eigenvalue claim (k=νΣf/Σa is
  flux-shape-independent → degenerate).
- **Non-flat flux (anti-#4 / H2):** the 4g case has a genuinely non-flat
  φ* (`[0.64,0.46,0.25,0.57]`), so `⟨Σa,φ*⟩` is NOT `Σa·(flat)` — a
  flux-shape error in how the functional contracts groups is activated.
  The 2g case is ADMITTED only because its independent PROD/ABSN pins
  still constrain the per-group νΣf and Σa values (even though its
  spectrum is coincidentally flat); but the 4g case is MANDATORY and
  carries the flux-shape teeth.
- **Asymmetric inputs (Mode-2):** mixture A 2g has asymmetric `SigP=[0.025,
  0.2]`, `χ=[1,0]`, asymmetric `SigS`; 4g `SigP=[0.014,0.026,0.125,0.245]`
  is monotone-asymmetric. A νΣf↔φ or group-axis-swap in the functional
  produces a different contraction.
- **Heterogeneous companion (anti-#4):** add a THIRD row — a 2-region
  fuel(A)+absorber(C) heterogeneous SN solve, pinning
  `solver.fission_op.production_rate.evaluate(φ_converged).sum()` against
  `solver.compute_production_rate(φ_converged)` (the existing
  volume-weighted SN path) — see 2.D. This activates the spatial
  redistribution the homogeneous case nulls.

### 2.C Mode-7 / Mode-10 activation declarations (per case)

| Case | ACTIVATES | NULLS | Mode-10 risk? |
|---|---|---|---|
| A 2g ∞-medium (functional dot) | per-group νΣf & Σa magnitudes; group contraction axis | spatial redistribution (1 region); flux-shape (φ* flat coincidentally) | none — independent pins constrain each group value O(1) |
| A 4g ∞-medium (functional dot) | per-group νΣf & Σa; group contraction; **flux-shape** (φ* non-flat) | spatial redistribution | none — non-flat φ* makes a group-weight error O(1) in the dot |
| A+C heterogeneous SN (2.D) | spatial redistribution; volume measure; full F.apply path | — | the volume-measure term: see 2.D mutation |

No term in Part 2 is activated-but-unconstrained: each functional's
output is the DOMINANT quantity in its pin (a dot product), so a
sign/magnitude error in any group is O(1) in the measured scalar. (This
is unlike an MMS sub-floor forcing — here the functional value IS the
measurand.)

### 2.D The eigenvalue + heterogeneous rows (the `l1` upgrade)

Two additional rows lift the claim from "functional value" to
"eigenvalue, in the full solver":

- **Row E1 (l1, eigenvalue):** run the SN k-eigenvalue solver on a
  homogeneous reflective box of mixture A (2g AND 4g), get `solver.keff`
  and the converged `φ`. Assert:
  1. `solver.keff == kinf_closed_form` (rtol ≈ 3e-3 — the ~0.3% transport
     correction from the true SN value vs ∞-medium closed form; this is a
     CROSS-CHECK bound, not a precision target, per AGENT.md §2).
  2. `solver.fission_op.production_rate.evaluate(φ).sum()` reproduces
     `solver.compute_production_rate(φ)` to the FP-non-associativity bound
     (these MAY differ in volume weighting — see note). This pins that the
     NEW `ReactionRateFunctional`-backed `production_rate` is wired into
     the operator the solver actually uses.
- **Row E2 (foundation, heterogeneous):** 2-region fuel(A)+absorber(C)
  SN solve; assert `production_rate.evaluate(φ).sum()` agrees with an
  independent volume-weighted hand computation of `∫νΣf·φ dV`. Activates
  the spatial term the ∞-medium nulls.

**NOTE on volume weighting (a real seam to flag):** `ReactionRateFunctional.evaluate`
returns the per-cell DENSITY (no volume measure — `production_rate_functional.py:61`),
whereas `SNSolver.compute_production_rate` folds in `volume_measure`
(`solver.py:1095`). So `production_rate.evaluate(φ).sum()` ≠
`compute_production_rate(φ)` UNLESS volumes are uniform and unit. **The
plan flags this to the implementer:** Part 5's LATER fold of
`compute_keff`/`compute_production_rate` onto `ReactionRateFunctional`
MUST decide where the volume measure lives (in the functional, or in a
`∫·dV` wrapper). Row E1.2 is only a 0-ULP pin if volumes are unit;
otherwise it is an rtol agreement gate with the volume factor made
explicit. **This is a Mode-10-adjacent trap:** if the test uses a
unit-volume mesh, the volume-measure term is NULLED and a bug in volume
weighting is invisible — Row E2 (heterogeneous, non-unit volumes if the
mesh has them) is the isolating case that activates it.

---

## PART 3 — The `outer(Vector, Functional)` intrinsic-property gate

**Goal:** test the dyad's DEFINING LAW, not just its usage. (Memory
`feedback_test_intrinsic_properties.md`: every math-bearing type ships a
test of its defining laws.)

**File target:** `tests/numerics/test_outer_dyad.py` (new), or extend
`tests/numerics/test_tensor_product_operator.py::TestRankOneOperator`.
`@pytest.mark.l0`.

The laws (write each as a row):

1. **Dyad action law:** `outer(v, w).apply(x) == v * w.evaluate(x)` for
   arbitrary `v` (Vector/ndarray), `w` (Functional), `x`. Reference: the
   RHS computed via the functional's `evaluate` AND, independently, via a
   hand `v * float(w_array @ x)` Python-loop dot. Assert 0 ULP vs the
   `v*evaluate` form; `nulp=4` vs the hand-loop form (reduction-tree).
2. **Rank-1 structure:** `outer(v,w).apply(x)` is proportional to `v` for
   every `x` (the output lies in `span{v}`). Test: for two distinct
   inputs `x1, x2` with `w.evaluate(x1) != 0 != w.evaluate(x2)`, assert
   `apply(x1)/w.evaluate(x1) == apply(x2)/w.evaluate(x2) == v` (elementwise,
   `nulp`-bounded). This is the *defining* rank-1 property — the column
   space is 1-D.
3. **`RankOneOperator` ≡ `outer` of its factors** (the equivalence the
   carve preserves): `RankOneOperator(v, w_arr, axis=0).apply(x) ==
   outer(v, ReactionRateFunctional(w_field)).apply(x)` at 0 ULP. This is
   the bridge row proving the new factory reproduces the old constructor.
   (If `RankOneOperator` is RETIRED in favor of `outer`, this row instead
   pins `outer` against the hand outer product and the *old* rows migrate
   per Part 5.)
4. **Capability law:** `outer(v,w).capabilities == frozenset({CAP_APPLY})`
   — rank-1 is non-invertible by construction (no `CAP_SOLVE`). Negative
   row: `outer(v,w).solve(...)` raises `MissingCapability`.
5. **Linearity in the functional argument** (the functional is linear in
   `x`): `outer(v,w).apply(a*x1 + b*x2) == a*outer(v,w).apply(x1) +
   b*outer(v,w).apply(x2)` (`nulp`-bounded) — confirms the row-factor is a
   genuine linear functional, not an affine/nonlinear map.

**Structural-independence note:** the references for rows 1–3 are hand
outer products / hand dots (Python loops where 0-ULP is claimed), NOT the
functional under test. Row 2's rank-1 check is reference-free (an internal
consistency law) — acceptable because it tests a STRUCTURAL invariant
(column space dimension), not a value (vv-principles L11 permits internal
structural laws; the value rows 1,3 carry the independent reference).

---

## PART 4 — Mutation probes (prove the gates BITE)

Every gate credited as evidence MUST be shown ABLE to RED under the
canonical `python -O -m pytest`. Reversion is **in-process monkeypatch**
(prod files are uncommitted on `refactor/operator-inverse-algebra` —
`git checkout`/`restore`/`stash` would DESTROY uncommitted work, per
`.claude/rules/process-discipline.md`). NEVER `git checkout` a held file.

| # | Mutation (in-process) | Expected RED gate(s) | Confirms |
|---|---|---|---|
| M1 | `monkeypatch.setattr` the functional's `nu_sigma_f.values` to a scaled copy (`×1.5`) | Part 2 PROD pin (independent); Part 1 B.1 hand-loop; Part 2 E2 | Σx in the functional is constrained |
| M2 | `monkeypatch` `FissionOperator.chi` (the reconstruction `v`/χ) to a permuted spectrum | Part 1 B.1; Part 2 E1; Part 3 row 2 (rank-1 v) | the reconstruction χ is constrained |
| M3 | swap production↔absorption: build `ReactionRateFunctional(Σa)` where `νΣf` is expected (monkeypatch `mat_xs.fission_production_field` to return the absorption field) | Part 2 PROD pin (0.159 vs 0.085 → RED); Part 2 k∞ ratio | production vs absorption XS distinguished |
| M4 | mutate the dyad's contraction axis: `monkeypatch` the functional's reduction `axis` (or construct `outer` with `axis=1`) | Part 3 row 1 (dyad action); Part 1 β multi-dim; Part 2 (wrong contraction → wrong dot) | the contraction axis is constrained |
| M5 | **(Mode-11 self-check)** make `RankOneOperator.apply` (or `outer.apply`) use the inline `(right*x).sum()` and BYPASS `functional.evaluate` | Part 1.B sentinel (counter stays 0 → RED) | the canary EXECUTES the new reader |
| M6 | **(Mode-8 self-check)** confirm every Part-2/3 gate uses `np.testing.*`/`require`/`pytest.fail` — grep the new test files for `^\s*assert ` ; any bare assert under `-O` is a NO-OP | n/a (audit) | gates fire under `-O` |

**Mutation-validation protocol (the implementer/qa runs this, NOT in plan
mode):** for each Mx, apply the monkeypatch INSIDE a throwaway test
(`@pytest.mark.skip` in committed form), run `python -O -m pytest <file> -x`,
CONFIRM the named gate (and ONLY the intended class of gate) reddens, then
remove. A `catches("ERR-NNN")` is NOT claimed unless the EXACT documented
bug reddens THIS test. (No ERR-NNN exists for this carve yet — if the carve
surfaces a real bug, file one per the Log-every-caught-bug directive.)

**Mode-10 explicit check (M3 sharpening):** swapping production↔absorption
must move PROD by O(1) (0.159→0.085, a 47% change) — well above any
floor — so the term is constrained, not merely activated. The 4g case
amplifies this (different PROD/ABSN ratio). If a swap left PROD unmoved
(e.g. because the test used a mixture where νΣf≈Σa), the term would be
activated-but-unconstrained → REJECT that fixture.

---

## PART 5 — Test-migration map

For each existing test that touches the carved surface — classify
behavioral-contract (REWIRE) / API-smoke (DELETE) / procedural-twin
(RETIRE→replaced by Part 2). Inventory via `grep -rn "ProductionRateFunctional\|RankOneOperator" tests/`.

| Test file / row | Classification | Action |
|---|---|---|
| `tests/transport/test_production_rate_functional.py` (S5 correctness, `.evaluate` vs hand dot) | behavioral contract | **REWIRE** to `ReactionRateFunctional` (rename the SUT; same `evaluate` contract). Its B-rows become the production-XS instance of the generalized `test_reaction_rate_functional.py`. |
| `tests/transport/test_functional_category.py` (`ProductionRateFunctional` satisfies `Functional`, not `LinearOperator`) | behavioral contract (Protocol membership) | **REWIRE** the foil to `ReactionRateFunctional`. The category law is unchanged. |
| `tests/transport/test_integral_kernel_category.py` (`ProductionRateFunctional` foil = field→scalar, NOT a Kernel; line 145, 312) | behavioral contract | **REWIRE** foil class name → `ReactionRateFunctional`. |
| `tests/transport/_functional_helpers.py` / `_integral_kernel_helpers.py` (`require_production_rate_property`, probes for `ProductionRateFunctional`) | helper (SUT-probe) | **REWIRE** the import-probe names; keep `hand_derived_fission_emission` (structurally-independent ref — KEEP verbatim). |
| `tests/sn/operators/test_fission_kernel_crosscheck.py` B.2 (`χ·production_rate.evaluate ≡ F.apply`, 0-ULP) | **procedural twin** (vv L11 weak) | **RETIRE the correctness-pretension; KEEP as a demarcated EQUIVALENCE row** (it still de-risks that the new property reproduces the arm), with B.1 (hand loop) carrying correctness AND Part 2 (closed-form k∞) carrying the structural-independence upgrade. Update SUT name to `ReactionRateFunctional`. Do NOT delete — it is the Mode-11 "reads the property off the live operator" check; just relabel it equivalence-not-correctness (it already is, lines 176-217). |
| `tests/sn/operators/test_fission_kernel_crosscheck.py` B.1 (hand-loop correctness) | behavioral contract (correctness anchor) | **KEEP unchanged** (rtol 1e-13; structurally independent). |
| `tests/sn/operators/test_fission_operator.py::TestRankOneTensorProductKernel` (`kernel.left is emission_spectrum`; `kernel.apply ≡ apply`) | behavioral contract | **REWIRE if** `kernel` now builds `outer(χ, ReactionRateFunctional(νΣf))` — the `kernel.left is χ` identity may change to a `kernel.reconstruction is χ` / `kernel.functional.nu_sigma_f` identity. The `kernel.apply ≡ apply` single-source row stays (0-ULP). |
| `tests/numerics/test_tensor_product_operator.py::TestRankOneOperator` (9 L0 rows) | behavioral contract (primitive) | **REWIRE the constructor-signature rows** if `RankOneOperator(left,right,axis)` → `outer(reconstruction, functional)`. The apply-law rows migrate to Part 3. If `RankOneOperator` is RETIRED, delete the constructor-rejection rows (API-smoke for a gone symbol) and move the apply rows to `test_outer_dyad.py`. |
| any `RankOneOperator` API-smoke (`test_constructor_rejects_*`, `test_negative_axis_normalised`) | API-smoke for the OLD constructor | **DELETE if the old constructor is retired**; else KEEP. |

**Net:** `ProductionRateFunctional` symbol retires; its behavioral tests
REWIRE to `ReactionRateFunctional`. The B.2 procedural twin is downgraded
to a demarcated equivalence row and the STRUCTURAL-INDEPENDENCE coverage
moves to Part 2 (closed-form k∞, independent PROD/ABSN pins). No coverage
is lost: every retired/twin row's real contract is re-homed.

---

## Mode-10 / activated-but-unconstrained flags (summary)

1. **Volume-measure term (Row E1.2 / E2):** `ReactionRateFunctional.evaluate`
   carries NO volume measure; `compute_production_rate` does. On a
   unit-volume mesh the volume term is NULLED and a volume-weighting bug is
   invisible. **Isolating fix:** Row E2 uses a heterogeneous mesh with
   non-unit cell volumes (or the LATER `∫·dV` wrapper fold makes the
   measure explicit and testable). The implementer MUST declare where the
   volume measure lives before E1.2 can be a 0-ULP pin vs an rtol gate.
2. **Flat-spectrum 2g coincidence:** mixture A 2g's φ* is flat, so its
   k∞ ratio is flux-shape-blind. **Isolating fix:** the 4g case (non-flat
   φ*) is MANDATORY and carries the flux-shape teeth; 2g is kept only for
   its independent per-group PROD/ABSN pins.
3. **Mode-11 sentinel target ambiguity:** if the fused `kernel` keeps a
   raw `RankOneOperator(χ,νΣf)` and the functional is only a semantic
   reading (today's `production_rate` pattern), the matvec does NOT route
   through `evaluate` — then a sentinel requiring `evaluate` to fire is a
   FALSE gate. **Resolution:** the implementer declares whether the
   row-factor IS the functional (matvec calls `evaluate`) or the
   functional is a parallel reading. The plan ships the sentinel matched
   to that decision (1.B).

## Open decisions to put to the user (AskUserQuestion candidates)

- Does the carved fused `kernel` matvec ROUTE THROUGH `Functional.evaluate`
  (then Part 1.B sentinel + M5 are live), or is `ReactionRateFunctional`
  the semantic reading only (then Part 2 carries coverage and the sentinel
  targets the `outer` apply)?
- Where does the volume measure live after the LATER `compute_keff` fold
  (in the functional, or a `∫·dV` wrapper)? Determines whether Row E1.2 is
  0-ULP or rtol.
- Is `RankOneOperator` RETIRED in favor of a free `outer(v, functional)`,
  or KEPT and re-expressed? Determines the Part-5 delete-vs-rewire split
  for the 9 `TestRankOneOperator` rows.
