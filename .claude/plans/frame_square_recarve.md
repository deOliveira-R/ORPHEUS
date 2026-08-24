# The frame square — minted faces, metric truth, and the Riesz floor (CS4b frame re-carve)

**⚠ INSTALL ON APPROVAL (step 0):** this file was written in plan mode outside the
repo. The first execution act is installing it as
`ORPHEUS/.claude/plans/frame_square_recarve.md`, replacing the campaign plan §5
⛔ FRAME-FACTORY block's charter text with a pointer + the ruling digest, and
updating the topic memory's ▶ NEXT. The user will compact context between plan
approval and execution — **this file is the resume surface; trust it over any
summary.**

---

## Context — why

The S4-amendment ruled "an operator is not an operator without its two spaces"
and landed it at the root. Its A1 phase bound the field spaces **to the frame** —
a proxy: it bound the *factory* when the ruling was about the *products*. The
2026-08-23/24 design session (5 pressure-test rounds + the desktop master-report
reconciliation) converged the correct shape and **measured** a latent metric
defect on the way. This plan executes the converged design: the frame reverts to
a shared `(basis, measure)` factory that **mints bound typed operator faces**,
and the moment space's inner product becomes the one Parseval certifies —
the precondition of S6's G6.3/G6.8 adjoint gates.

Branch: `feature/cs1-energy-space` (BRANCH-HOLD stands; the ≥90-min pre-merge
gate runs once, later, on the user's go). Canonical runner:
`.venv/bin/python -O -m pytest -p no:randomly` SERIAL. Surgical posture: main
agent writes, user steers.

## The ruled floor (all ✅ user, 2026-08-23/24 unless noted)

1. **The frame is the shared operator factory** — identity reverts to
   `(basis, measure)`; it mints precisely-typed bound faces, taking the one
   input it cannot know (the field space) at mint time. A1's derivation logic,
   content-equality admission, and gates carry over into the mint (completion,
   not revert).
2. **Generators in the family, everything else compositional.** Minted faces
   are `AnalysisOperator`/`ReconstructionOperator` members
   (`orpheus/numerics/projection.py` — whose own docstring anticipates exactly
   this: "`AnalysisOperator[AngularFlux, HarmonicMomentFlux]`"). Every other
   analysis-shaped doer is a *word*: `[M]` census 2026-08-23 — 1 in-family
   realization pair, 6 out-of-family doers (A1 verbs; `BulkAnalysisOperator`;
   `RadialCharacteristicReconstruction`; `frame.project`; loss-representation
   in-sweep accumulation; DSA's ℓ=1 table row).
3. **Primal/dual lives in carrier generics** — no `Primal*`/`Dual*` class
   names; four faces = instances of two classes; dual faces **derived** (shared
   kernel now, `dual()` functor at CS4c), never independently implemented.
   Precisely-typed faces are the algebra citizens.
4. **F2 = β-local Riesz**: consumers see ordinary algebra;
   `A.H = A.domain.riesz_raise @ A.dual() @ A.codomain.riesz_lower` is the
   *definition*, with space-minted legs. Legs derive from the amendment's
   spaces-demand — no new abstract members. **Execution of legs/dual-typing is
   CS4c**, not this plan (recorded debts below).
5. **F1 = Proposal B** (typed `DualSpace`, primal back-reference so
   `(V*)* is V` by reference), scoped to genuinely-functional ends — the pivot
   identification at the angular/nodal level stays. Desktop ruling (relayed,
   confirmed): *"covariance is representation-bound, therefore
   construction-bound, therefore never a property of Data."* Execution at CS4c.
6. **The discrimination stack**: axes / covariance (Gram pair,
   Parseval-checkable) / role+regularity tags. Flux-vs-source is **tags** — the
   Gram's refusal to separate them is the Riesz-isometry theorem, not a gap.
7. **The triple survives** (bulk, Γ ladder, `FullFieldSpace ⊕` as posed
   codomain). `WithTrace` (per-axis, post ray-discussion) adds the *fourth*
   object (solution-bulk W) — **Phase-S territory, recorded only**.

## The measured foundation `[M]` (probes in `scratch/probe_f1_parseval*.py`, 2026-08-24)

- **Sphere (LS S4, L=1; exactness floor 1e-16):** stored moment metric
  `4π/(2ℓ+1)` is the **wrong side** for the carried covariant moments —
  Parseval ratio 118.7 stored vs **1.000…** with `(2ℓ+1)/4π`. The machinery is
  self-consistent (adjoint identity 5.5e-16): the defect is *what is stored*,
  not how it is used.
- **`M.H = R/W` to 5.6e-17** with the correct metric — the frame square closes
  with one scalar; the shipped `1/W` kernel prefactor **is** that scalar. Under
  the stored metric, `M.H` is off the physical adjoint by exactly
  `(4π/(2ℓ+1))²` per ℓ.
- **The theorem (exact, unconditional):** `φ = Mψ = (YᵀWY)c` identically, so
  the Parseval metric is the **inverse of the frame's discrete Gram** — a
  property of *(basis ⊗ measure)*. Witness that a basis constant cannot be
  right: **slab** GL has total weight **2** (not 4π), live slots `[1,1,3]` per
  degree, and a **non-diagonal** live Gram (off-diag ~1.15) — no diagonal
  candidate satisfies Parseval there.
- **Exposure:** chains are immune (interior metrics cancel — the amendment's
  composite law; the `RΛM` kernel is embedding-invariant via the addition
  theorem). **End-of-chain adjoints are exposed = exactly S6's gates.** No
  current production consumer found (estimate, not exhaustive — re-check at
  F-0 with `grep -rn "\.H" orpheus/ | grep -i moment`).
- Factor homes today (3): the space metric
  (`spherical_harmonic_space.py` `from_L`), the basis's live `(2ℓ+1)` in
  `reconstruct` (`spherical_harmonic_basis.py:329` — applied at coefficient
  level, NOT baked into the table), the kernel's `1/W`.
  `InverseMetricOperator` (`numerics/operator.py:3275`) is the raise leg with
  dishonest endomorphic typing — retype at CS4c.

---

## Step F-0 — the metric truth

**Goal.** The codomain metric a frame's faces expose is the one Parseval
certifies, so an end-exposed `.H` on a face is the physical Hilbert adjoint.

**Proposed means** (mechanics; verify shapes at execution):

1. `FrameBase` (`orpheus/numerics/frame.py`) gains a cached
   `discrete_gram` — the full K×K coefficient-side Gram
   `einsum(w, table, table)` (K = coefficient count; tiny — O(N·K²) once per
   frame) — and a diagonality verdict (off-diagonal magnitude relative to the
   diagonal scale; declare via the existing `GramStructure` vocabulary /
   a frame-level property, mirroring `FrameBase.gram`'s DENSE-refusal pattern).
2. **Diagonal arm** (sphere cubatures, indicator frames): `basis_space` (and
   through the existing override chain, `test_space`) becomes
   `replace(basis.space, inner_product_weights=<inverse diagonal>)`. The
   `is`-identity `test_space is basis_space` is preserved by the existing
   GalerkinFrame override; `(name, shape)` equality is metric-blind `[M]`, so
   no consumer's `==` moves.
3. **Non-diagonal arm** (slab): keep the basis's continuum metric, record
   `parseval unavailable — non-diagonal discrete Gram` on the frame (loud,
   readable), and scope the Parseval gate to diagonal frames. The honest
   matrix-metric home is the CS4c leg-as-operator machinery (debt recorded).
4. The **basis keeps** `metric_per_ell` = `4π/(2ℓ+1)` — it is the continuum
   Gram and `project`/`gram`'s cross-Gram vocabulary (two jobs now separated:
   basis space = continuum; frame-dressed space = discrete Parseval).
   The `gram` row-sum probe path is untouched (it probes `MR` live).
5. Re-key the stored-metric adjoint **docstrings**:
   `analyze_transpose`/`reconstruct_transpose` in
   `spherical_harmonic_basis.py` (the "`M* = g_C · S₀`" /
   "`(2ℓ+1)²/4π`" formulas change to the Parseval forms, `M.H = R/W` on
   diagonal frames), and any `test_frame.py` pins of those values
   (enumerate by running the suite — expected small).

**Done when** (checkable):
- New gate (promote the probe): for ≥2 sphere quadrature families and L ∈
  {1, 2}, band-limited ψ satisfies
  `‖Mψ‖_codomain == ‖ψ‖_W` (rtol 1e-12) **and** `M.H(y) ≈ R(y)/W`
  (the closure pin). §6c witness: the gate REDS under the pre-repair stored
  metric — run once as an in-process mutation, record in the gate docstring.
- Slab frame reports its non-diagonal verdict; the gate skips it *visibly*
  (an explicit skip naming the CS4c debt, not a silent absence).
- Full fast set green; `npx pyright orpheus/` = 0.

**Hazards.** (a) The dressing is FrameBase-general — indicator-frame tests
(`tests/numerics/test_indicator_basis.py`, `test_weighted_indicator_basis.py`)
may pin face-adjoint values; their Parseval metric is `1/V_R` (same theorem) —
re-key, don't scope-narrow, unless reds cascade beyond faces (then narrow to
the angular mint and record the general item against S2). (b) The transport
carriers' moment space is a **product**; check how `FunctionSpace.__mul__`
composes `inner_product_weights` before assuming the SH factor's dressing
reaches the product metric (one read; F-1's mint is where it must hold).
(c) 0-ULP canary: the kernel `conjugate` chain reads tables and weights, never
the codomain metric — assert bit-identity on the aniso cross-check test anyway.

## Step F-1 — the mint

**Goal.** The typed faces are bound operators in the role-ABC family, minted by
the shared frame; the frame is shareable again (frames equal ⟺ shared
projection); the space-less kernel path is legal again (A1 symptom 3 repaired);
S6's adjoint gates have their object.

**Proposed means:**

1. **Two face classes** in `orpheus/transport/frames/` (new module or inside
   `harmonic_frame.py`): `HarmonicAnalysisOperator(AnalysisOperator[D, C])`
   and `HarmonicReconstructionOperator(ReconstructionOperator[D, C])` —
   constructor `(frame, domain_space→derived codomain, carrier pair)`;
   `apply` admits exactly its carrier by content-equality against the bound
   spaces (A1's admission text survives); kernels delegate to the frame's
   numerics faces (shared table — the "derived, not independent" guarantee at
   the array level); `apply_transpose` wired so `.H` works through
   `_AdjointOperator` on the F-0 metrics.
2. **Mint verbs on `HarmonicFrame`**, one per face-with-a-consumer
   (defer-until-consumer, per the `BulkAnalysisOperator` docstring's own
   precedent): `flux_analysis_on(space)` (consumers: windowing, loss-repr,
   S6 gates) and `source_reconstruction_on(space)` (consumer: the windowed
   in-scatter arm, `scattering.py:1256`). The other two faces are documented
   as mint-ready, unminted. `_derive_moment_space` moves into the mint path
   and dresses the SH factor with the frame's F-0 metric. The A1 refusal
   ("needs the posed composite space") relocates from
   `ScatteringOperator.frame` to the mint call.
3. **`HarmonicFrame` reverts**: fields `angular_space`/`moment_space` and the
   `analyse`/`reconstruct` verbs **retire** (aggressive retirement + test
   migration); `__init__` back to `(basis, measure)`; `from_galerkin(frame)`
   keeps only the SH-narrowing job. Identity returns to the table's identity.
4. **`ScatteringOperator`** (`transport/operators/scattering.py`): the `frame`
   property drops the space-refusal (`:518-521`) and the `angular_space=`
   argument (`:522`); mints cached faces at binding; the windowed arm re-keys —
   `:1243` `angular_target = self.frame.angular_space` →
   the minted face's `.codomain`-side read, `:1256` `self.frame.reconstruct(…)`
   → the minted `source_reconstruction` face's `apply`. The frame accessor
   itself STAYS (the fused `conjugate` kernel — CS4c re-expresses it).
5. **`BulkAnalysisOperator`** (`sn/operators/windowing.py:71`): re-based as the
   spelled ⊕-word — constructor takes the **minted flux-analysis face** (not
   the frame); `codomain` = `FullFieldSpace.from_blocks(<moment bulk>, <trace>)`
   (`full_field_space.py:217` — the `None` debt dies); docstring re-described
   as `lift(M) ⊕ Id_trace`; the "must be the SCATTERING operator's own angular
   frame" single-source clause re-keys to the shared mint. `WindowedSweep`
   fusion untouched (evaluation strategy). Call sites: `sn/solver.py:898`,
   `tests/sn/operators/test_space_content_witnesses.py:200`,
   `tests/sn/solve/test_2d_anisotropic_windowing.py:295`.
6. **§6b call-site set** (complete, `[M]` grepped 2026-08-24): production =
   `scattering.py:522/:1243/:1256`, `windowing.py`, `solver.py:898`. Tests =
   `tests/transport/frames/test_harmonic_frame.py` (TestBinding re-keys:
   bound-instance → minted-face; admission refusals → face admission; derived
   content-`==` → face codomain), the two windowing test sites above,
   `tests/sn/operators/test_scattering_operator.py`,
   `tests/transport/test_reaction_rate_functional.py` (verify which verb it
   touches at execution). One step, one landing — the set is small enough that
   splitting would leave the tree broken between halves.

**Done when:**
- Two mints from one frame share the table (`face_a.frame is face_b.frame`,
  table `is`-shared) and two frames over the same `(basis, measure)` compare
  equal — the shareability gate.
- A wrong-space carrier is refused by the face (witness: any posed-space
  mismatch — constructible today).
- An S6-precursor gate consumes a minted face's `.H` and pins the `M.H = R/W`
  relation through the *transport* face (not just the numerics face).
- The aniso scattering canary is bit-identical; full fast set green; pyright 0.

## ⏹ LANDING LEDGER (2026-08-23)

| step | status | commit |
|---|---|---|
| Step 0 (bookkeeping) | ✅ LANDED | `5539168f` |
| F-0 (metric truth) | ✅ LANDED | `0317373d` |
| F-1 (the mint) | ✅ LANDED | `3dfea889` |

**Verification state at landing:** targeted sets green at the landed state —
numerics FULL + transport FULL (2959) at F-0; transport/frames (20),
sn/operators (1246, incl. the 0-ULP kernel crosscheck), sn/solve windowing
equivalence + mutation gates, space-content witnesses at F-1; pyright 0;
Sphinx `-W` clean; nexus `dead_references` 0.

**⏹ FULL-SUITE VERIFICATION COMPLETE 2026-08-24 02:02** (per-tree observable
run, `scratch/_f01_suite_driver.sh` → `scratch/_f01_full_suite.log`; the two
prior attempts died opaque — harness reap at ~90 min, then pytest's
capture+block-buffering hiding all progress — both recorded in #405):

| tree | result | wall |
|---|---|---|
| transport | ✅ 545 | 35 s |
| sn | 3296 ✅ + **2 F, both adjudicated** (below) | 2 h 22 m 44 s |
| mc | ✅ 55 (2 xfail) | 40 m 15 s |
| derivations | ⏹ DECLARED SKIP at ~145/1728 (rc=143) | 3 h 19 m spent |
| cross_method | ✅ 89 | 58 s |
| cp | ⏹ DECLARED SKIP at 66/154 (65 ✅ + one >52-min GC-bound test; rc=143) | 53 m spent |
| data | ✅ 219 | 4 s |
| diffusion | ✅ 113 | 2 s |
| geometry | ✅ 792 (4 skip) | 8 s |
| homogeneous | ✅ 50 | 2 s |
| moc | ✅ 124 | 2 m 03 s |
| numerics | ✅ 2418 (3 skip) — incl. the full Parseval suite | 5 m 26 s |
| catch-all (loose + _harness + _mutation) | ✅ 408 (5 xfail) | 34 s |

**The two sn reds, adjudicated:** (1) `test_si_gate_dispatch` — F-1's one
true blast-radius miss: duck-typed surrogates of the raised mint contract
(a stub spells the interface as a KWARG — invisible to every symbol grep;
only the suite finds it). Fixed `4aa7f951` (4/4 green). (2) `phase_e`
flux-shape crosscheck on `cyl_2g_3reg_folded_4x8_dd_n40` — **pre-existing,
`[M]` bit-identical**: re-run at pre-F-0 `5539168f` in a worktree fails with
the same 16-digit numbers (0.1267629835542986 > the 0.12 bound). Filed
**#404** (fold-era snapshot re-capture suspected; phase_d sibling passes).
Blocks the eventual pre-merge gate; not a frame-campaign defect.

**The declared skips are not campaign exposure**: derivations and cp hold
zero frame consumers (the §6b + `.H` sweeps), and both are structurally
unrunnable in a serial gate at their current shape — the measured per-tree
map and the tiering/caching/xdist decisions are **#405**. Every tree that
CAN see F-0/F-1 ran green at the final committed state.

**Execution deltas vs the proposed means (all refinements, no goal changes):**

- **F-0 negative leg is a PERMANENT gate**, not a run-once mutation: the
  pre-repair continuum metric is re-installed in-process (cached-property
  dict pre-seed) and Parseval asserted to FAIL (>10×; measured ~118 on the
  probe's seed). Stronger than the plan's §6c ask.
- **F-0 side finding:** LS4/LS8 measure DIAGONAL **and degree-exact to
  ~2e-15 at L=2** — the `test_mass_matrix_under_multiple_quadratures`
  comment "LS_8 has a 24% diagonal error at L=2, no LS order exact" was
  present-tense-FALSE (likely pre-#327); refuted in place. Parseval needs
  only DIAGONAL (any values); the closure additionally needs d·G = W —
  every shipped sphere family measures exact on both.
- ⚠ **The probes no longer reproduce their own headlines post-repair**
  (they read `frame.test_space.inner_product_weights` as "stored", which
  is now the DRESSED metric → their "stored" row prints 1.000; and the
  118.7 ratio is draw-dependent — archivist re-measured 81.4/65.2 on other
  seeds; the draw-independent statement is the per-ℓ `(4π/(2ℓ+1))²`
  factor). The live witnesses are the `test_parseval_*` gates; the probes
  stay as the DISCOVERY record only.
- **F-1 faces subclass the role ABCs UNPARAMETERIZED** + their own
  `Generic[AngularFieldT, MomentFieldT]` (the ZeroMorphism precedent): the
  carriers deliberately fail the numerics `Vector` protocol's endomorphic
  arithmetic (a source's `+` returns a union), so the per-carrier precision
  rides the face classes' own generics. `[M]` pyright reported the bound
  violation; 0 errors under this shape.
- **`angular_target` reads `self._interior_space` directly** (not the
  minted face's `.codomain`): same value by construction, and it avoids
  minting a source-reconstruction face on the L=0 path that never applies
  one.
- **A §6b near-miss caught by its gate:** `WindowedSweep.apply` read
  `self.p.frame` (an attribute read INSIDE windowing.py, beside the
  constructor sites the audit enumerated) — the windowing equivalence gates
  reddened on it immediately; re-keyed to `self.p.face.frame`. The
  Mode-11 counter-spy's `monkeypatch.setattr` on the retired verb failed
  LOUDLY (the moved-attribute hazard the project guards) and re-keyed to
  the minted face's `apply`.
- **Recorded wart (GitHub issue to file):** `HarmonicFrame` equality over
  the SAME measure object works (tuple-identity shortcut), but two
  CONTENT-equal distinct `DiscreteMeasure` instances make frame `==` RAISE
  (dataclass `__eq__` hits ndarray ambiguity) — pre-existing FrameBase
  behaviour, surfaced by the shareability gate's design.

## Recorded debts (do NOT execute here)

- **CS4c**: Riesz legs first-class (`riesz_lower/raise` space-minted;
  `InverseMetricOperator` retyped/absorbed); `DualSpace` (F1-B);
  `_AdjointOperator` retired into the leg composition (keep-and-prove until
  then); Λ's factor collection (`(2ℓ+1)/4π` → one leg); matrix-metric for the
  slab's non-diagonal Gram; `dual()` functor + the declared-predicates-become-
  theorems gates; the mint-and-forget of S's frame accessor.
- **Phase S records** (master report): per-axis `WithTrace` descriptor; D8
  trace-measure criterion (witnesses already shipped: `bea6a367` α-endpoint
  guard; the Moore–Penrose masking); RAY-CHAR's remaining half; S2
  `HarmonicAxis` refinement — *the metric is induced (frame-owned), not a
  basis constant, and not always diagonal*.
- **Opportunistic words** (record in issues, execute when touched): ray-fold =
  `ι ∘ broadcast ∘ R₂(Legendre over ray nodes)`; DSA's row =
  `restriction_{ℓ≤1} ∘ M`; `project` legs naming.

## Step 0 (on approval, before compaction) — bookkeeping

1. Install this file at `ORPHEUS/.claude/plans/frame_square_recarve.md`.
2. Campaign plan §5: the ⛔ FRAME-FACTORY block gains
   `✅ RULED 2026-08-24 — design converged + measured; execution plan =
   frame_square_recarve.md` and drops its "do not implement" banner.
3. Master report: append the session deltas — the Parseval measurements as a
   new tightness-family instance; the S2 refinement; the desktop deltas
   (per-axis WithTrace, D8, RAY-CHAR, the Riesz/dual ruling + the
   construction-bound clause); `git add` the report (currently untracked).
4. Topic memory ▶ NEXT → "execute frame_square_recarve.md F-0 → F-1".
5. Commit the plan-layer changes (docs/plans only), push.

## Verification (end-to-end)

- Per-step: the new Parseval gate (with its documented pre-repair red), the
  shareability + admission + S6-precursor gates, the aniso 0-ULP canary.
- Suite: `.venv/bin/python -O -m pytest -p no:randomly` (fast set, SERIAL);
  `npx pyright orpheus/` (scope exactly `orpheus/` — 0 errors baseline).
- The probes stay in `scratch/` as the measurement record; the gate cites them.
