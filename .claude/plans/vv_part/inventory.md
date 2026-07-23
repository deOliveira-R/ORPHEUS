# V&V documentation estate — inventory (design basis for task #10)

> Read-only inventory, 2026-07-22, tree @ branch `docs/sn-doc-architecture`
> (graph built @ `d6276013`). Line numbers current at final read; the durable
> claims are the structure, the counts drift with every matrix regeneration.
> Explorer agent; incremental write — sections appear in work order.

## 0. Executive shape (the durable claim)

The estate is **three homes + one generated artifact + two machinery layers**:

1. `docs/theory/verification.rst` (661 ln) — the ORIGINAL 2026-Q1 "Verification
   Suite" page: SymPy-derivation doctrine, per-solver derivation walkthroughs
   (each `.. include::`-ing a build-generated `../_generated/*.rst` fragment),
   tolerance rationale, unit-test property lists. Sits INSIDE the theory part
   tree; catalog marks it "temp → V&V part".
2. `docs/verification/` (index 161 + matrix 1,408 generated + reference_solutions
   543 ln) — the 2026-Q2 "V&V infrastructure" home: vocabulary discipline
   (L0–L4), three pillars, operator-form taxonomy, ContinuousReferenceSolution
   contract, kernel primitives, and the auto-generated matrix.
3. `docs/testing/` (index 20 + architecture 633 + cross_method 345 ln) — the
   test-HARNESS home: the L0..L3 ladder + foundation bucket taxonomy, marker
   conventions, audit CLI, `:vv-status:` sentinel, cross-method regression
   protocol (L4 ruling).

They are "un-linked" in the precise sense: **no `:doc:` role crosses between
home 1 and homes 2/3 in either direction**; the only bridges are (a) raw path
strings (`docs/testing/architecture.rst` in prose — invisible to the Sphinx
`-W` gate), (b) global `:ref:` labels (`vv-foundation-tests`,
`vv-status-documented`, `vv-vocabulary`) which DO resolve project-wide, and
(c) `docs/index.rst` placing `testing/` + `verification/` as toctree siblings
under one caption while `theory/verification.rst` lives in a different part.

Level-taxonomy ownership is FRAGMENTED across four definition sites (§7).

---

## 1. The three homes

### 1a. `docs/theory/verification.rst` — 661 lines, label `theory-verification`

Heading tree (`[=]` h2 under the h1 title, `[-]` h3, `[~]` h4):

| line | lvl | heading | one-line characterization |
|---|---|---|---|
| 4 | h1 | Verification Suite | |
| 12 | h2 | Overview | The self-contained-SymPy-derivation doctrine; carries the STALE absolute "No cross-verification" claim (see contradiction below) |
| 29 | h2 | Architecture | `orpheus/derivations/` → {tests/, docs/_generated/} single-source flow; library-vs-workbench (`scratch/derivations/`) split |
| 61 | h3 | Reference-Values Registry (Eager vs Lazy) | `reference_values.py` two registries (legacy `VerificationCase` + `ContinuousReferenceSolution`); eager/lazy tiers; historical note on the RETIRED `_LAZY_LOADERS` Richardson side-channel (#290 P6) |
| 128 | h3 | Cross-Section Library | Synthetic A/B/C/D region library (`xs_library.py`), {1G,2G,4G} variants, consistency relation |
| 143 | h4 | P1 scattering anisotropy | Per-region μ̄ table with nuclide-analogue rationale |
| 181 | h2 | Verification Methodology | XS + analytical eigenvalue + principled tolerance; per-method tolerance-rationale table |
| 239 | h2 | Reference Case Types | Pillar-like trichotomy: Analytical (242) / Semi-Analytical (254) / Richardson-Extrapolated (269, labels `richardson-extrapolation` + eq `richardson-extrapolation-formula`) |
| 291 | h2 | Reference Cases | one line + `.. include:: ../_generated/verification_table.rst` |
| 297 | h2 | Homogeneous Infinite Medium | k_inf matrix-eigenvalue derivation + `.. include:: ../_generated/homogeneous_derivation.rst` |
| 318 | h2 | Discrete Ordinates (SN) | flat-flux homogeneous derivation + `.. include:: ../_generated/sn_derivation.rst` |
| 348 | h2 | Slab Collision Probability | E3 kernel sketch + `.. include:: ../_generated/cp_slab_derivation.rst` |
| 367 | h2 | Cylindrical Collision Probability | Ki3/Ki4 sketch + `.. include:: ../_generated/cp_cylinder_derivation.rst` |
| 388 | h3 | The 9-Case Grid (per geometry) | label `nine-case-cp-grid`; 3×3 groups×regions design rationale; RAW PATH to testing home at :437 |
| 440 | h2 | Method of Characteristics (MOC) | characteristic ODE + `.. include:: ../_generated/moc_derivation.rst` |
| 465 | h2 | Monte Carlo | collision-probability k derivation + `.. include:: ../_generated/mc_derivation.rst` |
| 483 | h2 | Diffusion (Buckling Eigenvalue) | bare-slab buckling + `.. include:: ../_generated/diffusion_derivation.rst` |
| 498 | h3 | 2-Region Fuel + Reflector: the Richardson Reference (RETIRED) | label `diffusion-2region-richardson`; historical record of the #290 P6 retirement (H4 anti-pattern), successor `derive_2rg_continuous` |
| 550 | h2 | Unit Tests | per-solver property lists (CP row-sum/reciprocity, SN GL weights, diffusion positivity, MOC balance, MC geometry) |
| 610 | h2 | Convergence Studies | SN spatial O(h²)/angular spectral, diffusion O(h²) — prose only |
| 642 | h2 | Running the Tests | STALE run-book: claims "56 tests"/"73 tests" (actual suite: 6,652) |

Cross-references OUT (toward other homes / theory pages):

- **Real roles**: `:ref:` `nine-case-cp-grid` (:83, internal), `diffusion-2region-richardson`
  (:124, internal), `:doc:` `/theory/methods/diffusion_1d` (:510), and the
  seealso block :633–639 → `theory-homogeneous`, `theory-collision-probability`,
  `theory-discrete-ordinates`, `theory-method-of-characteristics`,
  `theory-monte-carlo`, `theory-cross-section-data` + internal `richardson-extrapolation`.
- **RAW PATH STRING (verified)**: line **437**: "used by the test harness (see
  ``docs/testing/architecture.rst``)" — the ONLY pointer from this page to the
  testing home, inert literal text, invisible to `-W`. (Corpus plan said :436;
  drifted by 1.)
- **NO reference of any kind** to `docs/verification/` (matrix, index, or
  reference_solutions). Fully un-linked toward home 2.

Six `.. include:: ../_generated/*.rst` fragments (verification_table,
homogeneous/sn/cp_slab/cp_cylinder/moc/mc/diffusion derivations) — the page is
a SHELL around build-generated content; any relocation must preserve the
`_generated` relative paths or re-point them.

### 1b. `docs/verification/` — 3 files, 2,112 lines

**`index.rst` — 161 ln, label `verification-index`.**

| line | lvl | heading | characterization |
|---|---|---|---|
| 4 | h1 | Verification — V&V infrastructure & matrix | |
| 12 | h2 | Key Facts | read-first block; carries the L4-QUALIFIED no-cross-solver statement (the correct sibling of the theory page's absolute claim); names the three artefacts + error catalog/`catches` contract |
| 47 | h2 | Evidence classes — the three pillars and their boundaries | label `verification-evidence-classes`; closed-form / MMS / semi-analytical capability table (MMS ✗ eigenvalues); two-step correctness ladder; ERR-038 paper-floor example |
| 104 | h2 | Error catalog vs. paper-floor evidence | label `verification-error-catalog-vocabulary`; code-bug vs reference-precision-floor evidence classes |
| 130 | h2 | Three-meanings taxonomy (where this verification suite consumes it) | label `verification-three-meanings-cross-link`; α/β/γ Green's-function constructions; triple-match = highest-confidence L1 |
| 154 | h2 | Pages | toctree: reference_solutions, matrix |

Xrefs OUT: all REAL roles — `:ref:` `vv-vocabulary` (:22 → reference_solutions),
`:doc:` `reference_solutions` (:25), `matrix` (:29,:53), **`/theory/references/index`**
(:32,:133 — a real `:doc:` INTO the theory part), `:ref:` `theory-trajectory-resolvent`
(:138), `theory-singular-eigenfunction` (:142). Raw out-of-docs paths:
`.claude/skills/vv-principles/error_catalog.md` (:38),
`.claude/scratch/sanchez_chandrasekhar_gap.md` (:151). No link to `docs/testing/`.

**`matrix.rst` — 1,408 ln, GENERATED** (see §2 for generator). h2 sections with
current numbers: V&V level distribution (:12), Tagging source (:26), Module ×
level grid (:42, ~390 module rows), Equation coverage (:437, label→test-count
table), Orphan equations (:718, **254** orphans), Documented-only equations
(:978, **321** labels), L0 error-catalog coverage (:1305, ERR-001.. each with
catching-test count), Unmarked tests (:1384, **39** in 9 files).

Xrefs OUT: real `:ref:` `vv-status-documented` (:981) and `vv-foundation-tests`
(:1393) — both resolve into `testing/architecture.rst` (global labels). RAW
PATH STRINGS (verified): ``docs/theory/*.rst`` glob prose (:440),
``docs/testing/architecture.rst`` at :981 and :1392 (corpus plan said :836,:1220
— drifted after the 2026-07-22 regenerations). The raw-path+`:ref:` hybrid at
:981/:1392 comes from the generator's hardcoded strings
(`tools/verification/generate_matrix.py:221-222`, :263-264).

**`reference_solutions.rst` — 543 ln, label `verification-reference-solutions`.**

| line | lvl | heading | characterization |
|---|---|---|---|
| 3 | h1 | Reference Solutions — the ORPHEUS Verification Contract | the "binding contract" page |
| 23 | h3* | Vocabulary Discipline | **label `vv-vocabulary`** — Oberkampf/Roache definitions; verification=L0/L1/L2, validation=L3, benchmark=L4-forbidden-word; NOTE: points at "the V&V-level taxonomy in the project ``CLAUDE.md``" (:38–40) — a doc→CLAUDE.md dependency |
| 68 | h3 | Operator-Form Taxonomy | label `operator-form-taxonomy`; 6 disjoint operator-form tags (homogeneous / differential-sn / differential-moc / diffusion / integral-peierls / stochastic-transport) |
| 115 | h3 | Three meanings of "Green's function" | label `verification-greens-three-meanings`; α/β/γ realizations + pillar per meaning |
| 159 | h3 | The ContinuousReferenceSolution contract | frozen-dataclass API: phi/psi callables, phi_on_mesh, ProblemSpec, Provenance, operator-form tag |
| 190 | h3 | Registry and lookup | `get` (legacy VerificationCase) vs `continuous_get` (Phase-0) parallel registries, same key namespace |
| 214 | h3 | Kernel primitives E_n and Ki_n | **7 equation labels DEFINED HERE, outside docs/theory/** (`en-definition`, `en-kernel-special-values`, `en-kernel-derivative`, `en-kernel-integral`, `kin-definition`, `kin-kernel-special-values`, `kin-kernel-derivative`) — audit-glob visibility must be checked before any page move (§3) |
| 313 | h4 | Legacy naming discrepancy in BickleyTables (historical) | retired off-by-one Ki_n naming record (Phase B.4, #94) |
| 381 | h3 | Pillar-2 reference hardening — Atkinson product Nyström | label `verification-pillar-2-hardening`; the 3-step hardening pattern; ERR-036/037/038 evidence-class triage |
| 434 | h3 | Verification campaign — audit and migration plan | label `verification-campaign-migration`; per-module tier audit table; **the T1/T2/T2.5/T3-BANNED/T4-BANNED tier taxonomy** (a second, older level-taxonomy home) |
| 519 | h3 | References | Oberkampf, Roache, Case-Zweifel, Davison, A&S, Sood, Ganapol |

(* the page uses `=`-underline for h1 and `-` for its sections — one level
shallower than it renders in the toctree.)

Xrefs OUT: real `:ref:` `reference-solvers-three-meanings` (:122 → theory
references part), `:doc:` `matrix` (:152, :410), internal `operator-form-taxonomy`
(:182). Raw: ``CLAUDE.md`` (:39). No link to `docs/testing/` or to
`theory/verification.rst`.

### 1c. `docs/testing/` — 3 files, 998 lines

**`index.rst` — 20 ln.** h1 "Testing & V&V". Defines the four-level ladder in
one paragraph (L0 term / L1 equation / L2 integration / L3 validation + marker
spelling) — the most compact ladder statement in the corpus. toctree:
architecture, cross_method. No outbound roles.

**`architecture.rst` — 633 ln.** The harness bible.

| line | lvl | heading | characterization |
|---|---|---|---|
| 1 | h1 | Test-Harness Architecture | |
| 8 | h2 | Motivation | **THE L0..L3 ladder definition table** (:17–43) + foundation-bucket preview + the 1G-degeneracy warning (ERR-006 war story) |
| 53 | h2 | Design principles | one-source-of-truth markers; no DSL; ref()-inheritance; Nexus-native traceability; central audit; enforcement mode; pyright ratchet (#226) |
| 107 | h2 | Authoring a test | raw `pytest.mark.*` = the dominant convention |
| 120 | h3 | Raw pytest.mark decorators | l0/verifies/catches worked example; docstring `:math:` role → Nexus edge |
| 169 | h3 | Optional: the verify sugar layer | `tests/_harness/verify.py` decorators — "**not currently used by any test in the tree**" (matrix confirms: source `verify` = 0) |
| 194 | h3 | Parametrize over matching cases | `vv_cases(level=, method=, geometry=)` helper |
| 219 | h3 | Inheriting through ref() | conftest stamps level+verifies from `VerificationCase`; `level_source="case"` |
| 237 | h2 | Precedence order (most specific wins) | explicit > class-name > func-name > case > unmarked; highest-level tiebreak; foundation sorts below physics |
| 263 | h2 | The audit CLI | `python -m tests._harness.audit`; sample output is a FROZEN 2026-04 snapshot (519 tests / 3 orphans / 29 documented / 22 ERR — actual: 6,652 / 254 / 321 / 69+); flags `--json --untagged --gaps --strict`; "no CI yet, run by hand before every merge" |
| 357 | h2 | Foundation tests — software invariants outside the L0..L3 ladder | **label `vv-foundation-tests`** — the foundation-bucket taxonomy, when/when-not, tiebreak, audit reporting |
| 456 | h2 | Documented-only equations (:vv-status:) | **label `vv-status-documented`** — the sentinel syntax + 3 legitimate cases + fail-closed rules (see §3) |
| 520 | h2 | Selecting tests at runtime | marker expressions |
| 540 | h2 | tests/_harness package layout | file map (STALE: omits `predicates.py`, `pyright_ratchet.py`, `pyright_baseline.json`; `meshes.py` still described as stub) |
| 565 | h2 | Nexus integration | `:math:` docstring role → `references` edge; 3 requirements for the edge |
| 591 | h2 | Contributor checklist | 7-item checklist for new tests |

Xrefs OUT: only internal `:ref:`s (`vv-foundation-tests`, `vv-status-documented`)
+ `:func:`/`:mod:`/`:class:` roles into `tests._harness` and `orpheus.derivations`.
RAW PATHS: ``docs/theory/*.rst`` glob prose at :331, :409, :461, :473, :578, :598;
`.claude/skills/vv-principles/error_catalog.md` at :49, :83, :333, :619;
`tests/conftest.py` :129. **Zero references to `docs/verification/`** (not even
to the matrix it feeds) and zero `:doc:` to any theory page.

**`cross_method.rst` — 345 ln.** The cross-method regression protocol.

| line | lvl | heading | characterization |
|---|---|---|---|
| 1 | h1 | Cross-Method Regression Protocol | |
| 8 | h2 | Motivation | 4 continuous-reference solver families; CrossMethodCase / SolverAdapter protocol |
| 57 | h2 | V&V level mapping | **the L4 ruling**: "two solvers agreeing" is L4 (informational, not registered as a marker); shipped gates tagged L1 under the structurally-independent-pair convention (L1+L1 independent ⟹ L1-strength) |
| 92 | h2 | Reference contamination — the agreement-tolerance discipline | `agreement_tolerance = max(truth tolerances)`; anti-pattern #6 |
| 117 | h2 | Pillar tags — what each truth supports | closed-form / semi-analytical / MMS / **ancillary RESERVED-and-REJECTED** (foundation gate blocks L4-only truths) |
| 141 | h2 | Truth traceability | `truth_source` primary-citation field; Sood/KLL/GS/NM examples |
| 164 | h2 | Coverage matrix (current) | frozen "as of 2026-05-03": 17 cases, 84 collected tests |
| 215 | h2 | Adding a new solver | adapter workflow |
| 255 | h2 | Adding a new case | case workflow |
| 275 | h2 | Architectural seam — relationship to wave3 meta-registry | migration seam note |
| 290 | h2 | Multi-group cross-method coverage gap (acknowledged) | 1G-degeneracy honesty + lift paths |
| 323 | h2 | References | Sood, KLL, GS-1979, NM-1980, Sanchez-1986, Siewert-Thomas + skill pointers |

Xrefs OUT: ZERO `:doc:`/`:ref:` roles. Only `:mod:`/`:class:`/`:func:` roles
(tests.cross_method, orpheus.derivations.continuous) + raw `.claude/` paths
(vv-principles skill ×4, algebra-of-record, plans/wave3, scratch assessment,
agent-memory KLL memo). Fully un-linked from all doc pages.

### 1d. The soft contradiction (verified, quoted)

- `docs/theory/verification.rst:18–20`: "**No cross-verification** (comparing
  one solver against another) **is used** — each solver's verification stands on
  its own, as if every other solver were deleted." — absolute, no L4 qualifier.
- `docs/verification/index.rst:19–21` (the sibling that supplies the
  qualifier): "No cross-solver 'benchmarking' **stands in for** verification —
  that is level **L4 (informational)**, strictly distinct from L0–L3
  (correctness)."
- Shipped practice contradicts the absolute reading twice over:
  `docs/testing/cross_method.rst` documents a LIVE cross-method protocol (17
  cases / 84 tests) whose structurally-independent-pair gates are tagged
  L1-strength (:78–87), and `docs/theory/methods/method_of_characteristics.rst:743,1181`
  records a heterogeneous MOC-vs-CP cross-verification historically catching a
  bug. The theory page's sentence is the 2026-Q1 doctrine frozen before the
  L4/L1-strength refinement existed; the verification/index + cross_method pair
  is the current ruling.

Home ages (git): theory/verification.rst born 2026-04-03 (`23b90b08`), last
touched 2026-07-22; verification/index 2026-04-13, reference_solutions
2026-04-15 (last content touch 2026-05-04 — the most frozen page);
testing/architecture 2026-04-13, cross_method 2026-05-04 (both last touched
2026-07-13).

---

## 2. The matrix generator

- **Path**: `tools/verification/generate_matrix.py` (299 ln; package
  `tools/verification/` also holds `generate_capability_matrices.py`).
  Generator logic last changed 2026-05-01 (`4b3b8b92`, a path chore; last
  functional change `e431cdc0` 2026-04-14, foundation bucket).
- **Consumes**: `python -m tests._harness.audit --json` as a subprocess
  (`generate_matrix.py:47-56`). The audit in turn (a) runs
  `pytest --collect-only` to populate `tests._harness.registry.TEST_REGISTRY`
  via the `tests/conftest.py` collection hook; (b) scans
  `docs/theory/**/*.rst` (RECURSIVE `rglob` — `audit.py:339` — so the 2026-07
  part-tree restructure is fully covered) for `:label:` lines and
  `.. vv-status: <label> documented` sentinels; (c) scans ALL of `docs/`
  (minus `_build`) for labels, for the phantom-verifies gate (`audit.py:361-384`);
  (d) regexes `ERR-\d{3}` IDs out of
  `.claude/skills/vv-principles/error_catalog.md` (69 entries, ERR-001..069).
- **Emits**: `docs/verification/matrix.rst` (default `DEFAULT_OUT`,
  `generate_matrix.py:44`; single optional positional arg overrides). The page
  is TRACKED in git and committed after builds.
- **Invocation**: `docs/conf.py` `setup(app)` connects
  `_regenerate_verification_matrix` on **`builder-inited`** — i.e. EVERY
  `sphinx-build` regenerates the matrix before source collection (conf.py
  ~:120-138, closes #79). Also runnable manually
  (`python -m tools.verification.generate_matrix`). No Makefile target, no CI
  (architecture.rst:347: "There is no CI yet, so the harness is run by hand
  before every merge").
- **Last run vs file**: matrix.rst is CLEAN in git status and was last
  committed 2026-07-22 three times (`56f0705b` explicit "regenerate after the
  G1 label backfill", then `9f305683`, `77a5cb16` — regens riding Phase-H/#304
  commits). Generator (2026-05-01) far predates the content — the page moves
  with every build, the tool is stable.
- **Orphan semantics**: orphan = a `:label:` under `docs/theory/**` that is
  (i) not named by any `@pytest.mark.verifies("label")` on any collected test
  and (ii) not excluded by a `documented` sentinel
  (`audit.py:224-235`). Rendered as a flat alphabetized bullet list
  (matrix.rst:718+). **Current: 254 orphans** of 823 total theory labels − 321
  documented = 502 testable (49.4% of testable labels are orphan). The
  catalog's "243 after Phase G" is 11 stale — Phase H's root page added
  labeled equations and #304 renamed ~39 labels.
- **Phantom gate (inverse orphan, issue #224)**: verifies-targets with NO
  matching `:label:` anywhere under `docs/` — computed in the audit JSON and
  `--strict` (`audit.py:108-126,495-501`) but **NOT rendered on the matrix
  page** (the generator never reads `payload["phantom_verifies"]`). A silent
  reporting gap for the consolidation to close.
- **Consolidation constraint (load-bearing)**: the `--theory-dir` default
  `docs/theory` is the ORPHAN-GATE BOUNDARY. Moving a page INTO the theory
  tree puts its labels under the gate (e.g. reference_solutions.rst's 7
  kernel labels — `en-definition`, `kin-definition`, … — currently sit
  outside it, covered only by the phantom gate; the audit docstring
  explicitly blesses verifies-targets on verification pages,
  `audit.py:119-123`). Moving `theory/verification.rst` OUT removes its
  labels (`richardson-extrapolation-formula`, `richardson-diffusion`) from
  the gate. Any re-homing must either keep the V&V part inside
  `docs/theory/` or teach `--theory-dir` the new root.
- Sibling generator: `_regenerate_capability_matrices` (same `builder-inited`
  hook) emits `docs/theory/_<pkg>_capability_matrix.inc.rst` per
  `orpheus.derivations.continuous` package exposing `capability_rows()`.
  **Latent bug**: its exception path still calls `app.warn` — removed in
  Sphinx 9.1 — so a capability-generator failure would raise AttributeError;
  the matrix hook's identical path was fixed (`40e9ccc`) but this one wasn't.

---

## 3. The vv-status / verifies() machinery

### 3a. `vv-status` — a comment SENTINEL, not a directive

- There is **no Sphinx directive**: no `docs/_ext/` exists, `conf.py` has no
  implementation, and the vendored nexus checkout
  (`/Users/rodrigo/git/sphinxcontrib-nexus`) has ZERO occurrences of
  `vv-status`. The syntax `.. vv-status: <label> <status>` is a plain RST
  comment Sphinx silently strips; the ONLY parser is
  `tests/_harness/audit.py:337`:
  `re.compile(r"^\.\.\s+vv-status:\s+(\S+)\s+documented\s*$")`.
- **Recognized schema**: exactly ONE status — `documented`. Other words are
  "reserved for future use and silently ignored" (architecture.rst:507-510).
- **Wild usage vs schema** (measured): 348 sentinel lines total in
  `docs/theory/**` (0 elsewhere in docs) = 324 strict `documented` + **24
  inert lines carrying OTHER statuses** the parser ignores: `tested` ×13
  (peierls.rst ×4, peierls_nystrom.rst ×9), `verified` ×6+
  (trajectory_resolvent.rst ×5, sn/verification.rst:3660), `implemented` ×1
  (foundations/frame.rst:1722), plus 4 more of the same families. All 24
  inert labels have test coverage (0 intersect the current orphan list), so
  they are harmless shadow annotations — but they PROVE authors want the
  richer {documented, tested, verified, implemented} vocabulary the parser
  doesn't implement. Direct design input for the three-layer V&V part.
- **3 STALE `documented` sentinels** name labels that no longer exist
  (dropped by the `documented &= labels` fail-closed intersection,
  `audit.py:357`): `emission-kernels-btd`, `functional-category`,
  `tensor-network-shape-table` (casualties of the G1/P10/#304 label-rename
  campaigns). Hence 324 sentinels → 321 effective documented labels.
- **Doc-vs-code mismatch**: architecture.rst:504-506 rules "the sentinel must
  appear in the same RST file as the `:label:`; cross-file sentinels are not
  supported" — but `audit.py` accumulates labels and sentinels into GLOBAL
  sets before intersecting, so a cross-file sentinel WOULD work. Code
  outranks contract; the doc overstates.
- Per-file sentinel density (top): peierls_nystrom 47, frame 34,
  operator_algebra 33, trajectory_resolvent 30, peierls 27, discretization
  26, fn_method 20, singular_eigenfunction 10, thermal_hydraulics 9,
  sn/curvilinear_one_group 9, spherical_harmonics 9, coupled_block_operator
  9, boundary_conditions 9 … (grep `^\.\. vv-status:` for the full list).
- **Nexus connection: NONE via vv-status.** Nexus links tests↔equations via
  docstring `:math:`label`` roles → graph edges (architecture.rst:565-589,
  nexus ≥ 0.6). The sentinel lives only in the audit/matrix pipeline. The two
  pipelines are parallel and independent — a consolidation-relevant twin-path.

### 3b. verifies() / catches() / the registry

- **Registry** (`tests/_harness/registry.py`, 95 ln): frozen dataclass
  `TestMetadata(nodeid, file, level, level_source, equations, catches,
  case_names, slow)`; `TEST_REGISTRY: dict[nodeid, TestMetadata]` cleared and
  refilled each collection by the `tests/conftest.py` hook
  (`pytest_collection_modifyitems`); helpers `by_level()`, `by_equation()`.
  `VVLevel = Literal["L0","L1","L2","L3","foundation"]` — **L4 deliberately
  not registered** (cross_method.rst:67-68).
- **Level resolution precedence** (conftest, documented architecture.rst:237):
  explicit marker > `TestL<N>` class name > `test_l<N>_*` func name >
  inherited from `VerificationCase` via `ref()`/case param > unmarked.
  Highest-level tiebreak; foundation loses to any physics level.
- **Current usage counts** (grep, 2026-07-22): `pytest.mark.verifies(` — 454
  occurrences across 173 test files; `pytest.mark.catches(` — 160
  occurrences naming 69 distinct ERR tags = **69/69 catalog entries have a
  catcher** (matrix flags a miss as publication-blocker; none missing).
  Matrix tagging-source row: explicit 6,534 / class-name 46 / case 33 /
  verify 0 / unmarked 39.
- **`verify` sugar layer** (`tests/_harness/verify.py`, 141 ln):
  `verify.l0(equations=[...], catches=[...], slow=...)` decorators +
  `vv_cases(level=, method=, geometry=, n_groups=)` parametrizer over
  `VerificationCase` entries (zero-match raises). Confirmed **used by ZERO
  tests** (matrix source `verify` = 0) — supported dead weight; a
  consolidation candidate for either promotion or retirement.
- Harness dir actually contains: `__init__.py, verify.py, registry.py,
  audit.py, predicates.py (122 ln), meshes.py, xs.py, pyright_ratchet.py,
  pyright_baseline.json` — architecture.rst:540's layout map omits the last
  three (stale).

---

## 4. Per-method verification content to absorb

### 4a. The SN Verification chapter — `docs/theory/methods/sn/verification.rst`

**4,680 lines** (the catalog's 2,548 predates absorption growth — premise
stale). One h1 "Verification" + one h2 "Numerical Sensitivities" (:4604).
h3-level map (extents to next h3):

- :13 MMS 1-D slab · :183 heterogeneous 2G continuous-Σ slab · :379 2-D
  Cartesian separable MMS · :513 2-D 2G heterogeneous · :605 P1 anisotropic
  MMS · :683 curvilinear isotropic (radial DD-closure probe) · :807
  curvilinear anisotropic (redistribution probe) · :1021 the curvilinear
  aniso floor reconciled W1–W5 (τ-clamp vindication, pole-cell O(h) #233 /
  ERR-059, session trail)
- **:1701–3057 "The 2-D Cartesian LD stress MMS (D5b-S4)" — EXACTLY 1,356
  lines: this IS the "1,356-line 2-D LD stress MMS page" the task brief
  sought.** It was absorbed here (catalog row 7: former
  `linear_discontinuous.rst` split — UBLD core → `cartesian_multid.rst`, the
  stress MMS → `verification.rst`, both @ `eb96f424`). Sub-blocks: Mode-7
  ansatz design, Leg A slope-moment source (tensor-Legendre projection CRUX),
  Leg B boundary transverse-face-slope (`boundary_face_layout` CRUX,
  rank-discriminated inflow lift), the coherent boundary promise +
  property-vs-type seam (#263), Mode-10 mutation tables.
- :3058 composite fixed-source API (q = bulk ⊕ boundary) · :3347 non-vacuum
  prescribed-inflow MMS (P1-in-μ ansatz, affine-BC-to-RHS, Mode-7 map, Mode-9
  splitting invariance) · :4191 heterogeneous eigenvalue via Case
  singular-eigenfunction · :4483 homogeneous infinite medium · :4524
  heterogeneous convergence · :4556 why 1-group is degenerate · :4581
  spatial/angular convergence · :4588 property tests
- :4604 h2 Numerical Sensitivities: keff sensitivity table (421-group PWR
  slab), sources of variation, matching the MATLAB result.

Catalog disposition (sn_split_catalog.md:122): "Verification" + "Numerical
Sensitivities" marked **"(temp → V&V part)"** — this chapter is the
designated bulk donor to task #10.

### 4b. Other methods — Verification is an in-page h2, never a separate page

| method | file | verification content |
|---|---|---|
| CP | `methods/collision_probability.rst` (2,215 ln) | h2 Verification :1894; h3 Eigenvalue Verification Cases :1911; h3 Extended Verification (CP-20260405-005) :1964–1990 |
| MOC | `methods/method_of_characteristics.rst` (1,483 ln) | h2 Verification :881 (incl. "Why Homogeneous … Insufficient" :968, "Heterogeneous Cross-Verification" :990 — the page that practices what theory/verification.rst:18 denies, "Convergence Properties" :1018–1134); separate h2 "MMS Verification" :1340 + Convergence evidence :1437–1459 |
| MC | `methods/monte_carlo.rst` (1,315 ln) | h2 Analytical Verification :1091; h2 Verification Suite :1188–1217 |
| Diffusion | `methods/diffusion_1d.rst` (1,457 ln) | h2 Verification :1170–1253 |
| SN curvilinear-mg | `methods/sn/curvilinear_multigroup.rst` | h2 Verification :328–390 |

Additionally, verification content is WOVEN through the SN sub-book as
per-chapter gate sections (not to be uprooted blindly):
`loss_representation.rst` :476, :798, :1306 (FullFieldWavefront oracle),
:1742, :2588 (bit-identity vs principled-equivalence);
`cartesian_multid.rst` :835, :2086, :3225 (fuller-view oracle), :3822
(mutation redefinition); `curvilinear_numerics.rst` :2716, :2753 (#196 gate);
`slab_one_group.rst` :756 (verification hooks); `solver.rst` :493. These are
chapter-local contracts; the catalog's rebalance ruling (Phase 2: "ZERO MOVED
— machinery/driver/ABC keep contract-dense prose") applies.

---

## 5. Top-level IA wiring

`docs/index.rst` toctree structure (verified :16-58):

```
caption "Theory & Derivations"      → theory/index
caption "Architecture"              → architecture/index
caption "Testing & Verification"    → testing/index, verification/index   ← the two homes as SIBLINGS
caption "Development"               → development
caption "API Reference"             → api/*
```

Root prose (:8-14) routes "the V&V infrastructure and the auto-generated
verification matrix live at :doc:`verification/index`" — no mention of
`testing/` or `theory/verification` in prose.

`docs/theory/index.rst` (112 ln): the parts table (conventions / foundations
/ methods / references) has NO verification row; the references-part row
(:54-59) says references supply "the … truth values the verification suite
consumes". The V&V-flavored identity claim is :12-20 (Hébert+Stacey grep for
"verification|benchmark|manufactured solution" = zero hits — the corpus's
differentiator). The production/reference split rationale (:61-74) ends ":
Mixing the two … corrupts the V&V vocabulary." `theory/verification.rst`
rides in the **"Cross-cutting" toctree (:91-96) alongside `glossary`** — the
ruled rump position, with the catalog row marking it "temp → V&V part"
(sn_split_catalog.md:122) and task #10 named at :947 "(three-layer V&V
part)".

Inbound `:doc:`/`:ref:` links TO the homes (the catalog's "0 inbound" for
testing/architecture is FALSE):

- `testing/architecture` ← `development.rst:258` + **`theory/references/sood_registry.rst:694`**
  ("See :doc:`/testing/architecture` for the V&V level taxonomy") + testing/index toctree.
- `verification/*` ← index.rst:14 prose, development.rst:259,
  methods/index.rst:118, collision_probability.rst:720+:2112,
  diffusion_1d.rst:1199, sn/verification.rst:9, curvilinear_numerics.rst:73,
  references/index.rst:202+:356+:358 — the references part is the heaviest
  consumer.
- `theory/verification.rst` ← only `references/index.rst:360`
  (`:ref:`theory-verification``) + the theory/index Cross-cutting toctree.

---

## 6. Nexus verification state (graph @ `d6276013`)

`verification_audit(include_tests=true, group_by="module")` returned 246 KB;
I parsed 100% of it programmatically and extracted aggregates (did not
manually read the per-row dump).

Summary counts: **verified 628 · tested 912 · implemented 57 · documented 145
· orphan_code 7,773** over total_equations 830 (Nexus's count; the audit
CLI's own theory-label census is 823 — close, different scan rules).
**tests_declared 2,574 vs tests_inferred 60,283** — the "verified" total
rides overwhelmingly on heuristic edges; declared (marker/directive/registry)
evidence is ~4% of test-equation links. `orphan_code` counts code nodes with
no equation linkage — dominated by infrastructure functions, not a per-se gap.

Gap list: 202 (= 57 implemented-but-untested + 145 documented-only).
`group_by=module` yields only 3 coarse buckets (unassigned 145 / orpheus 44 /
tests 13), so the useful clustering is by theory page — top 12 of 37 pages:

```
21 foundations/frame            10 methods/collision_probability   8 methods/sn/slab_one_group
20 foundations/discretization   10 foundations/path_integral        7 foundations/operator_tensor_network
14 foundations/operator_algebra  9 thermal_hydraulics               7 references/fn_method
11 foundations/infinite_medium   8 foundations/operator_adjoint     7 methods/monte_carlo
```

Implemented-not-tested (the actionable slice, 57): collision_probability 6,
frame 5, operator_algebra 5, sn/slab_one_group 5, discretization 4,
fn_method 4, singular_eigenfunction 4, operator_tensor_network 3, …

stale_pages: 54 — dominated by `api/*` autodoc shells (git-timestamp
staleness, low signal for this campaign).

---

## 7. Who defines the ladder today (fragmented ownership)

Definition sites of L0..L3 / foundation / L4, in authority order:

1. **`docs/testing/architecture.rst:8-51`** — the canonical in-docs TABLE
   (L0 term / L1 equation / L2 integration / L3 validation) + the
   foundation-bucket taxonomy (:357, label `vv-foundation-tests`) + tiebreak
   rules. This is what other pages cite (sood_registry.rst:694).
2. `docs/testing/index.rst:4-11` — the compact one-paragraph restatement.
3. `docs/verification/reference_solutions.rst:23-63` (label `vv-vocabulary`)
   — maps verification=L0/L1/L2, validation=L3, benchmark=L4-forbidden-word;
   **delegates to "the V&V-level taxonomy in the project ``CLAUDE.md``"
   (:37-40) — a STALE pointer**: CLAUDE.md no longer carries the taxonomy
   (it lives in the `vv-principles` skill).
4. `docs/testing/cross_method.rst:57-91` — the L4 semantics ruling
   ("informational, parallel to the ladder"; L1-strength convention for
   structurally-independent pairs); cites
   `.claude/skills/vv-principles` §"V&V level taxonomy" as its authority.
5. `docs/verification/index.rst:18-22` — Key-Facts one-liner (L4
   informational, strictly distinct from L0–L3).
6. Code: `tests/_harness/registry.py:37` `VVLevel` Literal (L0..L3 +
   foundation; no L4) — the enforceable definition.
7. Out-of-docs: `.claude/skills/vv-principles/` — the authority both
   cross_method and the CLAUDE.md-pointer actually resolve to today.

Also coexisting: the **T1/T2/T2.5/T3-banned/T4-banned reference-TIER
taxonomy** (reference_solutions.rst:434-516, campaign-era) and the
**three-pillar evidence taxonomy** (verification/index.rst:47 +
cross_method.rst:117). The three-layer V&V part inherits FOUR overlapping
classification systems (ladder, pillars, tiers, operator-forms) with no
single page relating them.

---

## 8. Hygiene punch-list surfaced by this inventory (for task #10 scoping)

1. `theory/verification.rst:18-20` — absolute "No cross-verification" needs
   the L4 qualifier (§1d).
2. `theory/verification.rst:642-661` — "Running the Tests" claims 56/73
   tests; suite is 6,652.
3. Raw-path bridges (invisible to `-W`): theory/verification.rst:437;
   matrix.rst:981+:1392 — the matrix ones are HARDCODED IN THE GENERATOR
   (`generate_matrix.py:221-222,:263-264`); fix the tool, not the page.
4. `architecture.rst:263-355` — frozen 2026-04 audit sample (519 tests / 3
   orphans / 29 documented / 22 ERR; actual 6,652 / 254 / 321 / 69);
   :540 package-layout map missing predicates/pyright_ratchet/baseline.
5. `cross_method.rst:164-213` — coverage matrix frozen "as of 2026-05-03".
6. 3 stale `documented` sentinels naming dead labels; 24 inert sentinels
   using unimplemented statuses {tested, verified, implemented} (§3a).
7. Phantom-verifies gate computed but never rendered on the matrix page (§2).
8. `reference_solutions.rst:37-40` — stale CLAUDE.md taxonomy pointer (§7).
9. `architecture.rst:504-506` same-file sentinel rule not enforced by code (§3a).
10. conf.py `_regenerate_capability_matrices` still calls removed `app.warn`
    in its failure path (§2).
11. 39 unmarked tests in 9 files (matrix.rst:1384 lists them; the largest:
    `tests/numerics/test_assembled_operator.py` ×21).
12. `verify` sugar layer: supported, documented, used by zero tests —
    promote or retire (§3b).
