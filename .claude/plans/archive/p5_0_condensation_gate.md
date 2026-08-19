# P5.0 — Energy condensation gate design

**Mission**: design the verification gate for **energy condensation**
(GitHub #274) BEFORE the SUT (`EnergyGrid`, `GroupCondensation`,
`Mixture.condense`, `Solution.condense`) is built. Branch
`feature/sn-energy-condensation` @ `42427e8`. This is a write-tests-only
deliverable: the implementer (main agent, surgical-direct) builds
P5.1–P5.3 from this gate.

Energy condensation is the **energy-axis transpose** of the spatial
homogenization that already landed (P3). The same Petrov-Galerkin frame
machinery is reused verbatim; the only new types are `EnergyGrid` (the
energy analogue of a coarse `Mesh`) and a `GroupCondensation`
(partition-map) value object. The asymmetry law: `condense` returns
**portable** few-group XS (`dict[int, Mixture]`, mesh-DECOUPLED), vs
`homogenize` which is mesh-COUPLED (`MaterialMesh`).

---

## 0. Ground-truth reconciliation (verified against the real files)

Every API the gates reference, with `file:line`, re-confirmed this
session. Where the dispatch brief's stated ground truth needed
sharpening, it is flagged **[FLAG]**.

| Claim | Verified | Location |
| --- | --- | --- |
| `FrameBase.project(field) = G⁻¹ M f` (one line: `gram.apply_inverse_metric(analysis.apply(field))`) | ✅ | `orpheus/numerics/frame.py:275`, body line 292 |
| `PetrovGalerkinFrame(basis, measure, test_basis)`; `test_basis` REQUIRED | ✅ | `frame.py:295`, field at `:326` |
| `GalerkinFrame(basis, measure)` — test IS trial | ✅ | `frame.py:333`, `__init__` at `:353` |
| `LeastSquaresFrame` does NOT exist | ✅ | only `FrameBase`/`PetrovGalerkinFrame`/`GalerkinFrame` in `frame.py:111 __all__` |
| `IndicatorBasis(edges_per_axis=(coarse_edges,))` | ✅ | `orpheus/numerics/basis/indicator_basis.py:103`, field `:123` |
| `.evaluate(points) -> (N, n_coarse)` one-hot via `searchsorted(side="right")`; assumes **ASCENDING** edges | ✅ | `indicator_basis.py:142`, `searchsorted` at `:174` |
| `WeightedIndicatorBasis(indicator, weight)` — the TEST side; weight rides `analyze`'s trailing axes | ✅ | `orpheus/numerics/basis/weighted_indicator_basis.py:90`, `analyze` `:140`, `_weighted` `:153` |
| `Mixture` dataclass fields `SigC/SigL/SigF/SigP/SigT` (NG,), `SigS: list[csr]`, `Sig2: csr`, `chi:(NG,)`, `eg:(NG+1,)` | ✅ | `orpheus/data/macro_xs/mixture.py:53-61` |
| `Mixture.SigP` = νΣ_f (no `nubar` field); condense `SigP` | ✅ | `mixture.py:56`, `is_producing` keys on `SigP` `:73` |
| `absorption_xs` is a derived `@property` (no standalone field) | ✅ | `mixture.py:88` |
| `assert_balanced(atol=1e-9)`; residual = `|SigT − (SigC+SigL+SigF+total_scattering_xs+n2n_out_xs)|` | ✅ | `mixture.py:147`, residual `:118` |
| `eg` docstring wrongly says "Galerkin projection" | ✅ **[FLAG]** | `mixture.py:48-50` — implementer fixes to **Petrov-Galerkin** |
| `Solution.scalar_flux.values` is `(ng, *spatial)` GROUP-FIRST | ✅ | `orpheus/sn/solution.py:167`, `reaction_rate_density` `:313` |
| `Solution.homogenize` template (the spatial sibling) | ✅ | `solution.py:317` |
| `MaterialXSField.project_through(sigma_frame, emission_frame)` | ✅ | `orpheus/transport/mesh/material_xs_field.py:329` |
| `DiscreteMeasure(nodes, weights, support)`; the energy axis is a COUNTING measure (`weights=ones`) | ✅ | `orpheus/numerics/measure.py:124` |
| `FunctionSpace.apply_inverse_metric` — Moore–Penrose, zeroes empty regions | ✅ | `orpheus/numerics/space.py:232` |
| WIMS draft (untracked, present): `G69_EMAX/G172_EMAX` descending, `CONDENSE_172_TO_69` (1-based ranges), `validate()`, `boundary_mismatch_report()` | ✅ | `orpheus/data/wims_d_iaea_group_structures.py` (lines 52/73/126/138/152) |
| `EGB421.txt` present | ✅ | `orpheus/data/EGB421.txt`; **422 ascending boundaries = 421 groups** (1e-5 … 2e7 eV) |
| SUT absent: `EnergyGrid`, `GroupCondensation`, `Mixture.condense`, `Solution.condense` | ✅ | `grep` returns nothing across `orpheus/` + `tests/` |

### Discrepancies / flags vs the dispatch brief

1. **[FLAG — sharper trap]** The descending-vs-`searchsorted` trap is
   NOT "the partition breaks". It is **silent coarse-group-ORDER
   REVERSAL**. I verified concretely (realistic descending energy
   boundaries): feeding *descending* coarse edges to `IndicatorBasis`
   still produces row-sums == 1 (a "valid" partition, because
   `searchsorted` + `clip` always returns *some* in-range bin), but the
   coarse **column index is reversed** — fast fine groups land in the
   thermal coarse column. So a naive descending-edge feed yields coarse
   groups in REVERSED (thermal-first) order, breaking the canonical
   group-0-fastest convention end-to-end while every shape/partition
   check passes green. This is a Mode-6 (convention-drift) trap with a
   silent signature. **The `GroupCondensation` intrinsic-property test
   (below) MUST pin the coarse-group ORDER, not merely the partition
   completeness.** The implementer should build `GroupCondensation` from
   an **orientation-explicit** partition map (the WIMS
   `CONDENSE_172_TO_69` style: 1-based contiguous fast-first ranges),
   OR reverse edges to ascending AND reverse the resulting coarse axis
   back — and the test pins whichever choice produces fast-first coarse
   groups.

2. **[FLAG — G3 is ONE-ULP, not bit-identical]** I verified the
   scattering 2-axis collapse numerically (`sink = SigS @ T`, then
   `source = frame.project(sink)`). It reproduces the hand-summed
   in-scatter rate oracle to **2.2e-16 (one ULP)**, NOT bit-identically
   — the `@T` matmul reduction tree differs from the explicit
   group-by-group sum. This is FP-non-associativity, principled-
   equivalent per `vv-principles §bit-identity` (drift = reduction-depth
   × ULP). **The G3 rate gate MUST use `np.testing.assert_allclose` with
   a small `rtol`/`atol` (1e-12), NEVER exact `==`.** The vector
   channels (G1) ARE 0-ULP exact (single contraction), so they may use
   `==` — but I keep `assert_allclose(atol=1e-12)` uniformly for safety.

3. **[CORRECTED 2026-06-27 by user — real 421-group XS DO exist]** The
   original draft of this flag claimed "no 421-group XS library exists;
   synthesize + skipif." **That was wrong** (fresh-context error:
   `*.h5` is gitignored, so a git-only check reads "absent" — but the
   12 isotope `.h5` caches are present locally, regenerated from the
   `orpheus/data/micro_xs/*.GXS` sources by `convert_gxs_to_hdf5.py`,
   and ARE the production 421-group library). Verified this session:
   `from orpheus.data.macro_xs.recipes import pwr_like_mix, uo2_fuel` →
   real 421-group `Mixture`s, `ng=421`, `eg.shape=(422,)` DESCENDING
   (`eg[0]=2e7` eV → `eg[-1]=1e-5` eV), `pwr_like_mix().assert_balanced()
   == OK`. So **G5-legB loads REAL production XS** (`pwr_like_mix()` /
   `uo2_fuel()`), condenses to WIMS-69, and asserts shapes, balance, and
   χ fast-peaked on real data — NO synthesis, NO skipif. EGB421.txt is
   the grid oracle (it matches the real mixtures' `eg` reversed); the
   cross sections are the `.h5`/`recipes`. The G5-legB skeleton must be
   rewritten to the real-data path at P5.3 (the synthesize+skipif version
   in `tests/sn/test_condensation.py` is the to-be-replaced placeholder).

4. **[note]** The WIMS draft's `boundary_mismatch_report()` already
   reports **19** boundaries differing >0.05% (largest: 69-g1 ceiling
   10 MeV vs 172-g1 19.64 MeV at 96.4%, because 172-groups 1–5 sit above
   the 69-group ceiling). G5's derivation-validation leg reports these
   rather than asserting zero mismatch.

5. **[note]** `Quadrature.angular_frame(L)` returns a `GalerkinFrame`
   (`orpheus/numerics/quadrature/directional.py:85`) — the existing
   pure-Galerkin frame instantiation. Condensation uses
   `PetrovGalerkinFrame` (test = flux-weighted indicator), matching
   `Solution.homogenize` (`solution.py:427`).

---

## 1. The math the gate pins (per material, descending / fast-first)

Fine groups `g ∈ G`, per-material spectrum `φ_g`, partition membership
`T[g, G] ∈ {0,1}` (each fine group → exactly one coarse group,
contiguous). Verified numerically against the live frame this session.

- **Vector channels** `Σ ∈ {SigT, SigC, SigL, SigF, SigP}`:
  `Σ_G = (Σ_{g∈G} φ_g Σ_g) / (Σ_{g∈G} φ_g) = sigma_frame.project(Σ)`
  where `sigma_frame = PetrovGalerkinFrame(indicator, counting_measure,
  WeightedIndicatorBasis(indicator, φ))`. The counting measure has
  `weights = ones` (w=1), so `analyze` numerator = `Σ_{g∈G} φ_g Σ_g`
  and the diagonal Gram = `Φ_G = Σ_{g∈G} φ_g`. **0-ULP exact** vs hand
  sum.

- **Matrix channels** `SigS[ℓ], Sig2` `[g_from, g_to]`: **2-axis
  collapse** — sink `g_to` **summed** (`Σ_sink = Σ_fine @ T`, shape
  `(n_fine_from, n_coarse_to)`), then source `g_from` **flux-averaged**
  (`Σ_coarse = sigma_frame.project(Σ_sink)`). **This is the structural
  difference from `homogenize`**, which flux-weights BOTH axes
  (`sigma_frame.project(_gather_legendre(l))` at
  `material_xs_field.py:403`). Condensation sums the sink. **One-ULP
  exact** vs the hand-summed in-scatter rate.

- **chi**: **pure birth-group SUM** `χ_G = χ @ T` — NOT frame-projected.
  Verified: `χ @ T` preserves `Σχ = 1`; a flux-weighted projection would
  give `Σ ≈ 0.44 ≠ 1`, destroying the simplex. (Differs from spatial
  `homogenize`, which production-weights χ ACROSS CELLS via
  `emission_frame`. In condensation there is one material, so χ collapses
  by pure birth-group sum.)

- **Balance**: the condensed Mixture passes `assert_balanced()` iff the
  fine one does, because every removal channel collapses with the SAME
  per-coarse-group source-flux weight `Φ_G` (vectors) / `Φ_G` on the
  source axis (matrices), so the definitional identity
  `SigT = SigC + SigL + SigF + rowsum(SigS0) + rowsum(Sig2)` is
  preserved coarse-group-by-coarse-group. (The sink-sum on matrices is
  exactly what the rowsum needs: `rowsum(SigS0_coarse)[G] =
  Σ_{G'} project(SigS0 @ T)[G,G']`, and the same `Φ_G` divides SigT.)

---

## 2. Claim-layer + pillar declaration (vv-principles §1.5 gate)

Every gate row, with its claim layer, reference pillar, and structural-
independence ground. **No eigenvalue claims anywhere** — condensation is
a data-reduction operation, not a solve, so there is no `k` to verify;
the pillar question is "does the reduction preserve the rate functional",
answered by closed-form hand-summation.

| Gate | Claim layer | Pillar | Structurally-independent ground |
| --- | --- | --- | --- |
| G1 rate preservation | **flux-shape / rate-functional** (NOT eigenvalue) | closed-form | explicit per-fine-group Python hand-sum (`vv L11`) — NOT a `project` re-call |
| G2 within-group-varying discriminator | flux-shape | closed-form | hand-sum φ-weighted vs hand-computed arithmetic-average (the two must DIFFER) |
| G3 scattering 2-axis | flux-shape (matrix-rate) | closed-form | hand-summed in-scatter rate `Σ_{g∈G}Σ_{g'∈G'} φ_g Σ_s[g,g']` |
| G4 balance regression | software invariant | foundation | `assert_balanced` positive+negative pair (`vv #11`) |
| G5 WIMS + real-421 | partition correctness + integration | closed-form (partition) / integration (shapes) | the published Table 11.3 `CONDENSE_172_TO_69`; the containing-interval rule re-derived independently |
| G6 Mode-11 sentinel | execution-coverage | n/a (sentinel) | in-process monkeypatch counter on the rewired readers |
| IP `EnergyGrid` | type-invariant | foundation | the defining laws (strict-descending, positive, partition-complete) |
| IP `GroupCondensation` | type-invariant | foundation | containing-interval (one fine → exactly one coarse, contiguous, fast-first order) |

**Structural-independence terminus**: every correctness oracle (G1, G2,
G3) is a hand-coded Python `sum`/`for` over fine groups — it does NOT
call `frame.project` or a second `condense` path. A `project`-vs-`project`
check is twin-consistency (`vv L4`), FORBIDDEN as the correctness oracle.
This is the energy-axis copy of the spatial `homogenize` gate's discipline
(`tests/sn/test_homogenization.py:14-21`).

---

## 3. The gates (regime · oracle · mutation · Mode-7 declaration)

### G1 — Rate preservation (the L0/L1 correctness anchor)

- **Regime**: ≥2 coarse groups, ≥2 fine groups per coarse group
  (≥4 fine total), heterogeneous balanced fissile `Mixture`, a fine
  spectrum φ that **varies within each coarse group** (so the weight is
  load-bearing, see G2). Fast-first / descending eg.
- **Oracle** (structurally independent): for each vector channel and
  each coarse group `G`, `R_G^fine = Σ_{g∈G} φ_g Σ_g` computed by
  EXPLICIT Python summation. Coarse value ← `Mixture.condense` (i.e.
  `sigma_frame.project`); independent target ← the hand-sum. The coarse
  rate is `R_G = Σ_G · Φ_G` with `Φ_G = Σ_{g∈G} φ_g` (the group-sum).
  Assert `R_G == R_G^fine` to `atol=1e-12`.
- **vv level**: L1 (equation, against a closed-form rate functional);
  the per-channel isolation makes it L0-flavored (each `Σ` channel
  pinned independently).
- **AI failure modes defended**: #3 missing factor (a dropped φ or Φ_G
  breaks the rate), #5 index error (wrong fine→coarse membership), #2
  variable swap (wrong channel attr).
- **Mutation that RED**: drop the source flux-weight (condense with φ≡1
  / arithmetic average) → on a within-group-varying spectrum the rate is
  visibly wrong (this IS G2(b)); use the wrong `Φ_G` normalization (e.g.
  count instead of φ-sum) → rate off by a factor.

### G2 — Within-group-varying-spectrum discriminator (Mode-7 guard)

- **Regime**: a fine spectrum that VARIES within ≥1 coarse group (NOT
  flat). Concretely φ on the 4 fine groups = `[1, 4, 2, 0.5]` with
  coarse `G0={g0,g1}, G1={g2,g3}` — both coarse groups have a 4× / 4×
  within-group spread.
- **Oracle**: hand-compute BOTH (a) the φ-weighted collapse
  `Σ_G = (Σ φ_g Σ_g)/Σ φ_g` and (b) the arithmetic-average collapse
  `Σ_G^arith = mean_{g∈G} Σ_g`. Assert (a) preserves the rate (G1
  holds) AND (b) does NOT (the φ-weighted and arithmetic values are
  numerically distinct, and the rate computed from the arithmetic
  average ≠ `R_G^fine`).
- **vv level**: L1. **Mode-7 declaration**: the ansatz φ=`[1,4,2,0.5]`
  **ACTIVATES** within-group spectral variation in BOTH coarse groups;
  it **NULLS** nothing (a flat φ would null the weighting — Mode-7 trap —
  making flux-weighted ≡ arithmetic, so a wrong unweighted collapse
  accidentally agrees). The test asserts the activation: it checks the
  φ-weighted and arithmetic values genuinely differ (`discriminated`
  flag, mirroring `test_homogenization.py:283`), failing loudly if the
  fixture is too flat.
- **AI failure modes defended**: #3 missing factor (the whole point —
  proves the φ-weight is not droppable), Mode-7 (the test cannot be
  fooled by a flat config).

### G3 — Scattering 2-axis collapse mutation-probe

- **Regime**: an asymmetric (downscatter-heavy) fine `SigS[0]` 4×4, the
  same within-group-varying φ, ≥2 coarse groups. Asymmetry is mandatory
  (`vv` Signature 3 — a symmetric `SigS` makes the swap invisible).
- **Oracle** (structurally independent): the in-scatter rate
  `R_{G→G'} = Σ_{g∈G} Σ_{g'∈G'} φ_g Σ_s[g,g']` hand-summed (double Python
  loop). Assert `Φ_G · Σ_{s,coarse}[G,G'] == R_{G→G'}` to `atol=1e-12`
  (**one-ULP**, NOT exact — flag #2 above).
- **vv level**: L1. **Mode-7 declaration**: the asymmetric `SigS` +
  varying φ **ACTIVATE** the source/sink distinction AND the source
  flux-weight; flat φ or symmetric SigS would null them.
- **Mutations that MUST turn the gate RED** (each verified numerically
  this session to produce a *different* coarse matrix):
  - **(i) swap sink/source axes** — collapse `SigS.T @ T` then project
    (i.e. flux-weight the sink-group, sum the source). Differs → RED.
  - **(ii) sum BOTH axes** — `T.T @ SigS @ T` (drop the source
    flux-weight entirely). Differs → RED.
  - **(iii) project BOTH axes** — flux-weight the sink too
    (`project(project(SigS).T).T`). Differs → RED. (This is exactly the
    `homogenize` behavior, which is WRONG for condensation — so this
    mutation guards against "implementer copied homogenize verbatim".)
- **AI failure modes defended**: #2 variable swap (g_from↔g_to), #3
  missing factor (the sink-sum vs source-avg asymmetry), Mode-10/11
  (the mutations are value-RED, not merely structural).

### G4 — `assert_balanced` regression (vv #11 positive+negative)

- **Regime**: ≥2 fine groups per coarse group, ≥2 coarse groups (vv L2 —
  NEVER 1-group), heterogeneous channels (NON-zero SigL AND Sig2 so the
  full identity is exercised — every physical fixture has them zero, per
  `test_mixture_xs_balance.py:88`), balanced fine `Mixture`.
- **Oracle / assertions**:
  - **positive leg**: condense a balanced fine Mixture → the condensed
    Mixture passes `assert_balanced()` (MUST NOT raise).
  - **negative leg**: a deliberately rate-breaking mutant condense (e.g.
    condense SigT with the WRONG normalization, or hand-build a condensed
    Mixture whose SigT is off by ≫ atol) → `assert_balanced()` MUST raise.
    The broken instance is built BY HAND (not by perturbing the real
    condense), so the test pins the INVARIANT, not the raising behavior
    (`vv #11`).
- **vv level**: foundation (software invariant, no theory `:label:`).
- **AI failure modes defended**: #6 convention drift (a wrong channel
  weighting breaks balance), and the `vv #11` negative-only anti-pattern.

### G5 — WIMS Table-11.3 derivation-validation + real 421→69

- **Leg A (derivation-validation)**: derive the 172→69 partition by the
  **containing-interval rule** (each fine 172-group's representative
  energy → the coarse 69-group whose `[lower, upper)` contains it) and
  assert it reproduces the published `CONDENSE_172_TO_69`, **within the
  documented boundary tolerance** (`boundary_mismatch_report` lists 19
  non-coincident boundaries — the two structures were defined
  independently). The test REPORTS mismatches (does not silently pass):
  it asserts the partition agrees on the coincident-boundary groups and
  collects the known-19 as expected non-coincidences, failing if a NEW
  mismatch appears.
- **Leg B (real 421→69)**: load the real EGB421 grid (422 edges = 421
  groups, ascending → canonical descending), **synthesize** a fine
  421-group spectrum (fast-peaked, e.g. fission-spectrum-like) + simple
  balanced per-group XS on that grid (**[FLAG 3]**: no real 421-group XS
  data exists), condense to WIMS-69 via the analogous partition, assert
  shapes (`(69,)` vectors, `(69,69)` matrices), `assert_balanced()`, and
  χ **fast-peaked** (the bulk of χ in the low coarse-group indices, since
  group 0 = fast post-flip). `skipif` cleanly if a real 421-group XS
  loader ever lands.
- **vv level**: L1 (partition derivation) + L2 (real-structure
  integration). **≥2G throughout** (69G and 172G/421G — never 1G).
- **AI failure modes defended**: #5 index error (off-by-one in the 1-based
  partition ranges), #6 convention drift (the fast-first ordering on a
  real grid).

### G6 — Mode-11 execution sentinel

- **Regime**: any valid ≥2G condense.
- **Oracle / assertion**: in-process `monkeypatch`-count the rewired
  production readers and assert each fires `> 0`:
  - `Mixture.condense` (the new per-material verb) — count > 0;
  - `FrameBase.project` (the coefficient-extraction verb on the new
    path) — count > 0;
  - `WeightedIndicatorBasis.analyze` (the TEST-side reader — confirms the
    flux-weight is on the test side, not folded elsewhere) — count > 0.
  Sentinel via `np.testing` / `pytest.fail` assertion on the counter —
  **NEVER a bare `assert`** (Mode-8: `-O` strips it). This is the exact
  `homogenize` Mode-11 sentinel pattern (`test_homogenization.py:310`).
- **vv level**: foundation (execution-coverage).
- **AI failure modes defended**: Mode-11 (the named gate green but never
  calls the rewired reader — e.g. condense routes around `project` with a
  direct matmul), Mode-8 (the sentinel itself must fire under `-O`).

### IP — Intrinsic-property tests for the new types

Per the user's standing rule (every math-bearing type ships a test of its
DEFINING laws):

- **`EnergyGrid`** — strict-DESCENDING boundaries (the #265 invariant
  slice), all-positive energies, partition completeness (N boundaries → N−1
  groups with no gaps/overlaps). Positive (a valid grid constructs) AND
  negative (an ascending / non-monotone / non-positive grid raises). vv
  **Mode-8**: assertions via `np.testing`/`pytest.raises`, never bare
  `assert`.
- **`GroupCondensation`** — the containing-interval law: every fine group
  → exactly one coarse group; the per-coarse fine-group sets are
  contiguous (no gaps/overlaps); **the coarse-group ORDER is fast-first**
  (the [FLAG 1] orientation pin — a condensation built on raw descending
  edges would silently reverse this). Positive (a clean partition) AND
  negative (a gap / overlap / mis-ordered map raises or is detected).

---

## 4. Anti-recommendation compliance (the 5 the brief forbade)

1. ✅ Correctness oracle = hand-summed rate (G1/G2/G3), NEVER a
   `frame.project` re-call or second condense path (that is `vv L4`
   twin-consistency).
2. ✅ chi asserted as `χ @ T` (pure sum) + `Σχ` preservation, in a test
   SEPARATE from the projected channels — NOT routed through
   `frame.project`.
3. ✅ The discriminator (G2) + G3 use a within-group-VARYING spectrum
   `[1,4,2,0.5]` (never flat / single-group / homogeneous); the
   activation is asserted (`discriminated` flag).
4. ✅ Every always-on sentinel/canary (G6, G4 negative, IP) uses
   `np.testing.*` / `pytest.raises` / `pytest.fail` — no bare `assert`
   (Mode-8 under `-O`).
5. ✅ No new projection primitive / no `condense`-specific frame
   subclass: condense COMPOSES the existing `PetrovGalerkinFrame.project`;
   G6 pins that it routes through it.

---

## 5. Mode-7 activated/nulled summary (the config-blindness ledger)

| Gate | ACTIVATES | NULLS (and why it would matter) |
| --- | --- | --- |
| G1 | within-group spectral variation, every vector channel | — (flat φ would null the weighting → G2 catches that) |
| G2 | within-group spectral variation (the whole point) | the test FAILS if φ is too flat (`discriminated` guard) |
| G3 | source/sink asymmetry + source flux-weight (asymmetric SigS + varying φ) | symmetric SigS would null the swap-detection; flat φ would null the source-weight |
| G4 | the +SigL and +rowsum(Sig2) identity terms (NON-zero SigL/Sig2 fixture) | a zero-SigL/Sig2 fixture leaves those terms untested |
| G5 | fast-first ordering on a real grid; the 1-based partition ranges | a 1-group case would be degenerate (vv Cardinal) — ≥2G throughout |

**Minimum-catch confirmation**: the suite includes a within-group-varying
heterogeneous multi-group condense (G1+G2), the source/sink asymmetry
probe (G3), and a real-structure integration (G5) — activating every term
the condensation is most likely to get wrong (the source flux-weight, the
sink-sum-vs-source-avg asymmetry, and the fast-first group order).

---

## 6. Pre-implementation skeleton wiring

All three test files collect clean PRE-implementation via a flip-to-live
guard (mirroring `test_mixture_xs_balance.py:55-62`):

```python
_HAS_CONDENSE = hasattr(Mixture, "condense")        # or EnergyGrid import guard
_needs = pytest.mark.skipif(not _HAS_CONDENSE, reason="P5.0: SUT not yet implemented")
```

The **oracle math is implemented in-test** (the hand-sums), so the gates
go GREEN on contact with the SUT — no stub returns. Gates that exercise
the SUT carry `@_needs`; the pure-oracle / convention checks (the WIMS
partition derivation, the descending-edge trap demonstration) run
unconditionally (they need no SUT).

**Expected paste-back state**: xfails/skips (SUT absent), NOT errors —
collection MUST succeed.

---

## 7. vv-level → AI-failure-mode coverage matrix

| AI failure mode | Defending gate(s) |
| --- | --- |
| #1 sign flip | (n/a — condensation has no sign-bearing redistribution term; balance G4 would catch a sign error in a channel) |
| #2 variable swap (g_from↔g_to, wrong channel) | G3(i), G1 (per-channel) |
| #3 missing factor (dropped φ / Φ_G / sink-sum) | G1, G2, G3(ii) |
| #4 wrong recursion | (n/a — no recursion in condensation) |
| #5 index error (fine→coarse membership, 1-based ranges) | G1, G5, IP `GroupCondensation` |
| #6 convention drift (descending-edge order reversal, fast-first) | IP `GroupCondensation`, IP `EnergyGrid`, G5 |
| #7 MMS/flat-config blindness | G2 (within-group-varying activation guard) |
| #8 compiled-out assertion | every sentinel uses `np.testing`/`pytest.fail`/`raises` |
| #10 activated-but-unconstrained | G3 mutations are VALUE-red (verified differ), not merely structural |
| #11 gate-never-executes-rewired-path | G6 (monkeypatch counter on `condense`/`project`/`analyze`) |

---

## 8. Paste-back (L12) — verbatim pytest summary

Invocation:
`.venv/bin/python -O -m pytest tests/data/test_energy_grid.py
tests/data/test_mixture_condense.py tests/sn/test_condensation.py -q -rfE
--timeout=300 -p no:xdist -p no:cacheprovider`

```
9 passed, 39 skipped, 1 warning in 0.54s
```

(48 tests collected cleanly, 0 errors/failures. The 1 warning is the
benign `-O` "assertions not in test modules" notice — every test uses
`np.testing.*`/`pytest.fail`/`raises`, NOT bare `assert`, so nothing is
silently stripped, Mode-8-safe.)

**Split**: **9 passed** = the structurally-independent oracles + convention
traps that bite NOW with no SUT (both orientation-trap demonstrations, the
two G2 Mode-7 activation guards, the chi sum-vs-project oracle, the WIMS
partition + exact-derivation + boundary-non-coincidence legs, the EGB421
grid well-formedness). **39 skipped** = every SUT-exercising gate
(flip-to-live `@_needs`), GREEN-on-contact (oracle math implemented
in-test). The unconditional oracles' TEETH were verified in-process:

- **G3** — the design collapse (sink-summed/source-averaged) preserves the
  in-scatter rate; the homogenize-style **project-both-axes** collapse does
  NOT (rate breaks) → `test_in_scatter_rate_preserved` REDS on a
  homogenize-copy bug (the single most-likely implementer error). All three
  G3 mutations produce a numerically different coarse matrix → the `_differs`
  assertions fire.
- **Orientation trap** — descending coarse edges give coarse assignment
  `[1,1,0,0]` (thermal-first REVERSAL); fast-first would be `[0,0,1,1]`, so
  the assertion pins the reversal (non-tautological).
- **chi** — the flux-projected χ sums to 0.38 ≠ 1 for the fixture → the
  sum-vs-project `pytest.fail` guard is live.
- **G5** — the containing-interval derivation reproduces `CONDENSE_172_TO_69`
  at 172/172 EXACTLY; the 19 documented boundary non-coincidences are real
  (worst rel ≈ 0.96, the 172-above-69-ceiling) and pinned separately.

**Fixture-bug signal for the implementer**: writing the oracle helpers
caught a real trap — `make_mixture` (`xs_library.py:56`) hardcodes
`SigL = np.zeros(ng)` and silently DROPS any `sig_l` argument. A balanced
fixture with NON-zero SigL (needed for the G4 full-identity leg) MUST be
built DIRECTLY via the `Mixture(...)` constructor, NOT through
`make_mixture` (else SigT includes a SigL the Mixture zeros → imbalance).
The SUT's `Mixture.condense` must likewise carry SigL through the collapse
(it is a vector channel — condense it with the others); the G4 positive leg
+ `tests/data/test_mixture_condense.py::_balanced_fissile_4g` (hand-built,
SigL ≠ 0) pin this.

---

# ADDENDUM — fractional overlap (F1–F7), the NON-NESTED case (#274)

**Mission**: EXTEND the nested gate for FRACTIONAL-overlap re-binning. The
production case (421-group → WIMS-69/172) is **NON-NESTED**: a coarse
boundary can fall INSIDE a fine group, which must apportion a FRACTION of
its rate to each coarse group it overlaps. `GroupCondensation.table`
becomes **FRACTIONAL** `T[g,G]∈[0,1]`, partition-of-unity (rows sum to 1);
nested → one-hot (the regression-safe degenerate). All G1–G6 nested gates
are KEPT unchanged (they validate the degenerate case).

## Ground-truth reconciliation (verified against the real files this session)

| Claim | Verified | Location |
| --- | --- | --- |
| `EnergyGrid` / `GroupCondensation` (NESTED, one-hot) HAVE landed | ✅ | `orpheus/data/energy_grid.py:76` / `:153` |
| `GroupCondensation.table` is currently ONE-HOT (argmax of `coarse_of_fine`) | ✅ | `energy_grid.py:220-233` |
| span guard present (`coarse ceiling ≤ fine ceiling`) | ✅ | `energy_grid.py:184` |
| **`Mixture.condense` does NOT exist yet** | ✅ | `grep` empty in `mixture.py` + `solution.py` |
| **Fractional API absent** (`WithinGroupSpectrum`, `InverseEnergySpectrum`, `within_group` field, fractional `table`, `locally_interpolated`, upscaling guard) | ✅ | `grep` empty across `orpheus/` |
| `pwr_like_mix()` / `uo2_fuel()` are REAL 421-group balanced Mixtures | ✅ | `orpheus/data/macro_xs/recipes.py:57` / `:137`; verified `ng=421`, `eg.shape=(422,)` DESCENDING (`2e7…1e-5` eV), `assert_balanced()` OK |
| `WIMS_69` / `WIMS_172` are real `EnergyGrid` targets | ✅ | `wims_d_iaea_group_structures.py:146-147`; `WIMS_69` ceiling `1e7` floor `0`, `WIMS_172` ceiling `1.964e7` floor `1e-5` |
| 421 ceiling (`2e7`) > WIMS-69 ceiling (`1e7`) — span guard PASSES for 421→69 | ✅ | the over-ceiling 421 groups clamp into coarse 0 |

### [FLAG — the load-bearing F7 sharpening] "χ fast-peaked" is NOT `argmax==0`

The real `pwr_like_mix()` χ **peaks at 421-group index 61 (~1.15 MeV)**, NOT
group 0 — the 2e7-eV grid ceiling puts ~60 groups of small-but-nonzero
high-energy tail above the fission peak. **99.99 % of χ is above 1 keV**;
the χ-peak energy lands in **WIMS-69 coarse group ~4** (0 = fastest of 69).
So an `argmax(chi_coarse) == 0` assertion would **FALSE-RED on the real
spectrum**. F7 instead pins the physically-correct **cumulative-mass**
property: `chi[:69//2].sum() > 0.5` (the bulk of fission emission sits in
the fast half). The old synthesize-leg's `argmax==0` was only valid because
its synthetic χ was a literal step function in the fastest 10 % — it would
have been a brittle trap on real data. This is the config-blindness lesson
applied to a real fixture: the convenient synthetic spectrum nulled the
exact property (the high-energy tail) the real χ exercises.

## The gates added (each: regime · structurally-independent oracle · the RED mutation · proven teeth)

The oracle table `T` is computed **INDEPENDENTLY** of the SUT (hand-coded
energy overlap × the within-group model), NEVER read back via `cond.table`
(vv L11). The **STRADDLE FIXTURE** (the discriminating geometry): fine
`[1e6,1e4,1e2,1e0,1e-2]` (4 groups), coarse `[1e6,1e1,1e-2]` (2 groups) —
the `1e1`-eV cut falls INSIDE fine group g2 `[1e0,1e2)`, splitting it 1/E as
`T[2]=[0.5,0.5]` (lethargy ratio `ln(10)/ln(100)=0.5`). A non-straddle makes
fractional ≡ one-hot, so EVERY fractional gate uses THIS fixture.

| Gate | File · class | Regime | Structurally-independent oracle | RED mutation (the documented bug) |
| --- | --- | --- | --- | --- |
| **F1** straddle rate-preservation | `test_mixture_condense.py::TestF1StraddleRatePreservation` | straddle, within-group-varying φ=[1,3,2,0.5], 5 vector channels + scattering | `σ_G·Φ_G == Σ_g T[g,G]φ_g σ_g` with hand-built 1/E `T` | one-hot (g2 wholly to one coarse) → different rate split |
| **F2** partition of unity | `test_energy_grid.py::TestF2PartitionOfUnity` | straddle (fractional rows) + nested (one-hot) | `table.sum(axis=1)==1`; nested `array_equal` to the {0,1} one-hot | one-hot → g2 row is `{0,1}` (the "g2 row is one-hot" guard fires) |
| **F3** nested degeneracy | `test_energy_grid.py::TestF3NestedDegeneracy` | nested 8→3, coarse edges ON fine edges | fractional `table` == the contiguous fast-first one-hot | a fractional impl that ALWAYS splits (never snaps to {0,1}) |
| **F4** within-group discriminator | `test_energy_grid.py::TestF4WithinGroupModelOracle` (oracle) + `test_mixture_condense.py::TestF4ModelDiscriminatorSUT` (SUT) | straddle | 1/E fraction == `ln(hi_ov/lo_ov)/ln(hi_g/lo_g)`; flat-energy gives 90/99 ≠ 0.5 | one-hot ignores w(E) → 1/E ≡ flat-energy (model not load-bearing) |
| **F5** upscaling guard | `test_energy_grid.py::TestF5UpscalingGuard` | 64 fine → 200 "coarse" | `coarse.n_groups > fine.n_groups → ValueError`; +positive control (8 coarse OK) | no-op guard that never raises (or always raises) |
| **F6** local-interpolation report | `test_energy_grid.py::TestF6LocalInterpolationReport` | nested (empty) + straddle (non-empty) | `locally_interpolated` empty nested / ⊇ fractional columns straddle | one-hot interpolates nothing → empty report on a straddle |
| **F7** real 421→WIMS-69 succeeds | `test_condensation.py::test_real_pwr_421_to_wims69_condensation_succeeds` | REAL `pwr_like_mix()` → `WIMS_69`, NON-NESTED | no zero columns (`SigT>0` ∀ coarse), `assert_balanced`, χ fast-half-mass>0.5, `locally_interpolated` non-empty | one-hot leaves EMPTY coarse groups (3 empty under one-hot →69) → `SigT>0` fails |

**F7 replaces** the old `test_real_421_to_coarse_condensation_end_to_end`
synthesize-and-skip placeholder: NO synthesis, NO skipif — it is the L2
real-structure integration leg using the production library.

## TEETH PROVEN (the L4 discipline — every gate confirmed ABLE to RED)

Because the fractional SUT is not yet landed, the SUT-exercising gates were
validated via an **in-process reference fractional impl** (a correct
fractional `GroupCondensation`/`Mixture.condense` + a BUGGY one-hot
variant), monkeypatched in, gates run under both:

- **CORRECT mode** (fractional reference): all 19 F-class SUT tests **GREEN**
  (no false-RED on a valid impl) + F7 **GREEN** on real `pwr_like_mix()`.
- **ONEHOT-BUG mode** (the documented containing-interval bug): **11 RED** —
  F1 (all 5 vector channels + value-oracle + scattering), F2
  (`test_fractional_rows_sum_to_one`), F4 (both legs), F6
  (`test_straddle_reports_interpolated`); + **F7 RED** on exactly *"every
  WIMS-69 coarse group must be populated (SigT>0) — the empty-group bug must
  be gone"*. F3/F5 + the nested legs of F2/F6 stay **GREEN** (the bug only
  touches the straddle, not the degenerate case — correct discrimination).
- The **unconditional oracle legs** (F4 oracle-half: lethargy ratio, 1/E ≠
  flat-energy, partition-of-unity for any positive w) bite directly: a
  one-hot g2 split `[0,1] ≠ [0.5,0.5]` reds F4 leg-a; the validation harness
  was deleted after use (temp `/tmp` files, no repo artifact).

vv Mode-8 confirmed: every assertion is `np.testing.*`/`pytest.fail`/`raises`
(the run's only warning is the benign `-O` "assertions not in test modules"
notice).

## Paste-back (L12) — verbatim pytest summary

Invocation:
`.venv/bin/python -O -m pytest tests/data/test_energy_grid.py
tests/data/test_mixture_condense.py tests/sn/test_condensation.py -q -rfE
--timeout=300 -p no:xdist -p no:cacheprovider`

```
30 passed, 40 skipped, 1 warning in 0.62s
```

(70 tests collected cleanly, 0 errors/failures. **30 passed** = the
structurally-independent oracles + convention traps that bite NOW with no
fractional SUT — the 9 original nested oracles, the new F4 oracle-half
(3 legs: lethargy ratio, 1/E≠flat-energy, partition-of-unity), the EnergyGrid
intrinsic + orientation-trap legs, the WIMS partition derivation, the EGB421
grid well-formedness. **40 skipped** = every SUT-exercising gate
(flip-to-live `@_needs`/`@_needs_fractional`), GREEN-on-contact once the
implementer lands `Mixture.condense` + the fractional `within_group` API.)

## API surface the gates test against (`file:line`) + awkwardness flags

The gates probe the fractional API by **signature** (`"within_group" in
inspect.signature(GroupCondensation).parameters`) so they skip cleanly on
the current nested-only SUT and flip live the moment the field lands.

| API tested | Expected `file:line` (to be implemented) | Test reference |
| --- | --- | --- |
| `WithinGroupSpectrum` Protocol + `InverseEnergySpectrum` (`integrated_weight(lo,hi)=ln(hi/lo)`, floor-guarded) | `orpheus/data/energy_grid.py` (new) | F4 oracle (`_inv_e_weight`), F1/F2/F6/F7 pass `InverseEnergySpectrum()` |
| `GroupCondensation(fine, coarse, within_group=...)` field | `energy_grid.py:153` (add field) | all fractional gates |
| `GroupCondensation.table` FRACTIONAL `[0,1]`, rows sum 1 | `energy_grid.py:220` (generalize) | F1/F2/F3 |
| `GroupCondensation.locally_interpolated` (coarse cols with a straddle contribution) | `energy_grid.py` (new cached_property) | F6, F7 |
| upscaling guard `coarse.n_groups > fine.n_groups → ValueError` | `energy_grid.py:183` `__post_init__` | F5 |
| `Mixture.condense(condensation, spectrum) -> Mixture` | `mixture.py` (new) | F1/F4-SUT/F7 |
| `from_grids(fine, coarse, within_group=?)` | `energy_grid.py:271` (forward the kw) | gates try both spellings |

**Awkwardness flags for the implementer (I will reconcile):**

1. **`from_grids` must forward `within_group`** — the gates try
   `from_grids(fine, coarse, within_group=...)` first, then fall back to the
   direct ctor on `TypeError`. If `from_grids` does NOT forward the kw, the
   fall-back keeps the gates green, but the CANONICAL spelling at call sites
   (`from_grids`) would silently drop the model → a non-default model would
   be ignored. **Recommend `from_grids` forward `within_group` explicitly.**

2. **Over-ceiling / under-floor fine groups must CLAMP (preserve χ mass),
   not zero-row.** The 421-grid ceiling (`2e7`) exceeds WIMS-69's (`1e7`),
   so ~60 fast fine groups fall above the coarse ceiling. They must clamp
   wholly into coarse group 0 (matching the existing `coarse_of_fine` clamp
   at `energy_grid.py:217`), NOT get an all-zero `table` row — else their χ
   mass is dropped and `Σχ_coarse < 1` (F7's simplex assertion catches this;
   it surfaced in the reference impl during teeth-validation). The fractional
   `table` builder must apply the same clamp the one-hot `coarse_of_fine`
   does.

3. **`locally_interpolated` semantics** — the gates assert it lists the
   coarse columns that received a FRACTIONAL (straddle) contribution
   (F6: `fractional_cols ⊆ locally_interpolated`). If the implementer prefers
   "the FINE groups that straddled" instead of "the COARSE groups affected",
   F6 would need a one-line adjustment — flag which convention you pick and I
   will align F6.

4. **`Mixture.condense` must route the fractional `table` through BOTH
   axes of the scattering collapse** (sink summed with `T`, source
   flux-averaged with `T` — `R[G,G']=Σ_g Σ_g' T[g,G]φ_g Σ_s[g,g'] T[g',G']`).
   F1's `test_scattering_fractional_rate_preserved` pins the full
   double-fractional rate; the existing G3 (nested) pins the one-hot sink-sum.
   The G6 Mode-11 sentinel (condense routes through `FrameBase.project` +
   `WeightedIndicatorBasis.analyze`) still applies — the fractional table is
   the trial `IndicatorBasis.evaluate` generalization, the PG frame machinery
   is unchanged (per the cross-domain-attacker verdict).
