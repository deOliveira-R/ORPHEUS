# Test Architect — Lessons (hot digest)

Read at the START of every dispatch. This is the **behavioral index**: one rule
per entry, imperative and standalone. War stories, measured numbers and
`file:line` detail live in **`lessons_archive.md`** (sections L1–L38) — open only
the section a pointer names. The failure-mode TAXONOMY (Modes 1–12, three
pillars, anti-patterns #1–#17) lives in the preloaded **`vv-principles`** skill:
**cite it, never restate it.** Reference inventory + XS mixtures → `AGENT.md`;
per-carve RECIPES → the topic memos indexed by `MEMORY.md` §3.

**THE SPINE.** A plan is done not when the tests pass but when, for EVERY gate:
(a) a named mutation reddens it under `python -O`; (b) the reference is
structurally INDEPENDENT of the SUT; (c) the regime ACTIVATES the term the bug
lives in. Standing directives: `AGENT.md` §0.5 / §0.6 / §1.5.

---

## 1. Gates that cannot red — one rule, many disguises

**RULE: name the mutation that reddens each gate and RUN it before crediting the
gate as evidence.** `vv` Mode 8 catalogues seven shapes (compiled-out,
tautological, signature-tautological, misattributed-xfail, self-satisfied
`raises`, skip-swallowed, decayed `catches`) with detection recipes — read them
THERE. Below: only the shapes vv lacks, plus the repair recipes.

- **An xfail's flip-edit MUST touch a statement whose VALUE the production change
  determines.** I shipped one whose prescribed flip made it a character-for-
  character duplicate of the live flip-proof beside it; production landed
  mid-session and the row did not move. Discriminator: *diff the xfail body
  against its own flip-proof* — textually equal after the edit ⟹ ceremony. Repair:
  state the claim against the production ANSWER, not a hard-coded value. → `L33`
- **A gate can stay green while its REASON becomes false.** `B` was documented
  block-DIAGONAL; a wrap made it block-STRUCTURED; all three asserting rows
  stayed green because all three sit on a now-special-case fixture. When a phase
  falsifies a structural claim, **grep the claim's WORDS in `tests/`, not only
  its symbols**, re-scope those rows in the SAME change, and give the new
  structure its own positive gate. → `L33`
- **A pin naming a "legacy"/"reference"/"adapter" counterpart must be checked
  that the counterpart still EXISTS and is not the SUT under another name.** Two
  mechanisms, same `X == X` end state: a later-added DELEGATION (survivor calls
  the other side — replacing the production function with random nodes left it
  green) `L34e`; or a RETIREMENT that deletes the counterpart and re-points the
  comparison at the successor, where the survivor is the *caller* (`SNMesh`
  stopped owning `_setup_spherical` and now calls the very factory the test
  compares it to — a fully-garbaged factory left all 47 tests in the file green,
  while 29 gates elsewhere reddened). **Two probes, in this order, before any
  mutation battery: (1) is the "other side" literally the SAME OBJECT? (`is`, 5
  seconds — `face_areas` was, via a shared `cached_property`, so that leg was
  `array_equal(x, x)`); (2) garbage the ONE shared producer in EVERY module
  binding and see if the file notices.** Then re-scope rather than delete —
  garbaging only the CONSUMER-side binding reddened exactly the 15 cases, which
  is the wiring claim they honestly carry. → `L34e`
- **A docstring that names "the surviving pins" is a CLAIM to measure, not a
  handoff to trust** — including one you wrote an hour ago. A named anchor can be
  blind for a *structural* reason: the τ producer-equivalence gate was cited as
  pinning the connection-coefficient math, but a refactor had moved τ to the
  angular closure (a function of `(μ, w)` alone), so it passes untouched under a
  fully-garbaged geometry factory. Corollary: an L0 identity that "covers" a term
  may RECOMPUTE the production array instead of reading it (`dA/w` recomputed
  from `dA` and `w`), pinning the LAW while blind to the ARRAY — check which.
- **Retiring a runtime guard that had NO negative test makes its replacement's
  teeth NET-NEW, not migrated.** Grep `pytest.raises(match=<guard msg>)` before
  crediting a mechanism-swap as behavior-identical; if nothing asserts the old
  raise, WRITE the negative test the guard never had — and `match=` the SPECIFIC
  message, since a downstream crash on the same input satisfies a bare
  `pytest.raises`. → `L4`
- **When the fix is "raise the project's TYPED error instead of the builtin", a
  `pytest.raises(<builtin>)` gate is GREEN BEFORE the fix and after** — ORPHEUS's
  `BoundaryError(ValueError)` means `except ValueError` catches the bare error
  today while `except BoundaryError` provably does not. One word, total loss of
  the gate, and it reads correct in review. Name the SUBCLASS, and require the
  gate's pre-fix state to be RED (ship it `xfail(strict=True)` so the fix deletes
  the marker). → `L36f`
- **A `pytest.raises` on a refusal is teeth-less without MESSAGE legs, KEYED to
  the argument that triggers it.** A generic message keeps `exc.value.law` true;
  a blanket "the message names both completions" pins a FALSE reason on the α=0
  row, whose defect is different. Always pair with a positive control — else an
  arm that refuses everything also passes. → `L31`
- **REPAIRING blind gates (a different design problem from writing them).** Repair
  a decayed gate by re-posing onto a REGIME-INDEPENDENT mechanism, never by
  driving the fixture back into the regime: ERR-052's catastrophe was UNREACHABLE
  in its own fixture at any depth, so the repair asserts the output convention the
  fix establishes (`∫νΣ_f φ dV = 1`), true at every outer count — check
  reachability BEFORE trying to reach it, record the answer, and never compute the
  reference with the routine that ESTABLISHES it. **⛔ A brief's own proposed
  matvec row can be tautological on the very fixture the phase is built on:** on
  the P4 MMS slab BOTH faces declare `PrescribedInflow`, which P3 collapsed onto
  the zero morphism, so `[M]` `|B(x)| = 0.0` for a random `x` and `B(0)=0` /
  `B(2x)=2B(x)` hold with both sides structurally zero — no input can red them.
  **Before writing any linearity/homogeneity/additivity row, measure `|Op(x)|`
  on a random `x` and require it `> 0` as a committed ACTIVATION leg**; if it is
  zero, the honest gate is the STRUCTURAL claim ("this IS the zero morphism") on
  that fixture plus the linearity row on a fixture where the operator is
  non-trivial (here: prescribed on one face, REFLECTIVE on the other →
  `|B(x)| = 1.320`). → `L40c`. Repair a tautological-`raises`
  file guard-by-guard with the guard named at `file:line` in each docstring; the
  acceptance measurement is the PER-GUARD red table (0/14 → 12/14), every residual
  miss categorically out of scope and said so. Pointing tests at production entry
  points is itself an audit of production. Never let a repair reduce the catch
  rate — re-run the AUDITOR's own harness, not a re-implementation; the
  per-mutation red COUNTS are the signal. → `L28`

## 2. Harness discipline — the instrument lies before the code does

**RULE: an all-blind mutation verdict is a broken harness until a positive
control proves otherwise** — `vv` anti-#17 carries the rule and both ORPHEUS
instances (a privately-loaded test-module copy; a summary parser defeated by ANSI
codes). → `L34d`, `L35l`

- **Make the harness ASSERT its own installation — printing a banner is not
  enough, because nobody reads a banner in a captured `$out`.** Two more
  false-"0 caught" verdicts in one session, both in the SAFE-LOOKING direction:
  a per-field shell loop that silently DROPPED the `-p <plugin>` flag (six
  fields, six clean zeros — exposed only by contradicting an earlier all-fields
  run that reported 29 reds on the same set), and a per-field loop over
  ~80 s-per-case tests that TIMED OUT and again reported "0 red". Repairs, both
  cheap: `raise RuntimeError` inside the plugin when it rebinds 0 symbols, and
  grep the banner COUNT into the result line so a missing instrument is visible
  per row. A mutated run that gets *slower* (garbage destroys convergence:
  3.4 s → 80 s) will blow a timeout sized on the baseline — budget mutation runs
  off the MUTATED cost, never the green one.
- **MEASURE the reachable subset before rationing the battery — a brief's cost
  estimate names a DIRECTORY, not the tests the carve can reach.** A "≈5.5 min"
  budget measured **9.40 s** for the subset that mattered (35×), turning a
  16-mutation battery from "ration it" into ≈6 min. An over-stated cost silently
  shrinks the battery, which is the same loss as a blind gate and harder to see.
  Also re-measure the pre-declared RED: if the brief names one known failure,
  confirm it is the ONLY one, so every later red is attributable with no triage.
  → `L36e`
- **`git status` + `ls -la` mid-run, not just at the start: on a shared tree
  another agent may be rewriting the very code you are measuring.** The same
  `sed` range returned different prose 20 minutes apart; a second agent was
  running its own tree-wide mutation battery concurrently. ⭐ **Worse form,
  measured: a concurrent write can DEMOTE your own measurement to a
  value-compared-with-itself.** Dispatched to plan a carve "before" it landed, I
  probed production at 22:39 and again at 22:41; the implementer shipped the
  change in between, so my "bound production vs unbound production" bit-identity
  leg was bound-vs-bound and `array_equal` was true for free. Same end state as
  the lesson-#4 rewire demotion, arriving through the CLOCK. Defences: build the
  control EXPLICITLY (`np.asarray(x.perm).copy()`, then assert
  `control.domain is None` INSIDE the test so it cannot silently become bound),
  re-`git diff` before every claim, and treat a brief's "I am about to write X"
  as "X may already be on disk" — re-`ls` and re-score. Report the correction
  with `file:line` + evidence instead of racing the writer.

- **Run the teeth harness over your OWN new module before delivering it.** It
  flagged a gate I had just written — an activation guard that was a theorem
  about parity, true for every input, surviving authoring and a green run.
  → `L34d`
- **Measure "before" with `git show HEAD:<file> > <tmp sibling>`, never
  `git checkout`/`restore`/`stash`** — this tree carries uncommitted state; the
  sibling runs under the SAME mutation, collection and fixtures. Same rule for
  the committed side of any generated-artifact diff: a concurrent Sphinx build
  regenerating a file mid-session manufactured a confident, entirely WRONG
  finding. On a shared tree `ls -la`/mtime is part of the measurement. → `L28`
- **Never baseline against a tree being written** (`32 failed` mid-edit vs
  `9 failed` settled), and run `git status --porcelain -- orpheus/` at the START
  and before every claim — five of seven plan steps were already on disk when one
  plan began. **Ship a plan's gates as a RUNNABLE dry-run module**, not only
  fenced code: transcribing exposed two errors no reading would (a perturbation
  form `B(base+e) − B(base)` that is not bit-exact; a negative whose refusal had
  moved layers). → `L31`, `L33`
- **A mutation that reds by RAISING rather than comparing has attributed
  nothing** (a shape guard fired first) — name the attributable catcher
  separately. A mutation that fails to red because it is INAPPLICABLE to the
  fixture shape is not "no teeth". → `L31`, `L25`
- **Constructing a break-exactly-ONE-invariant mutant is a design problem.**
  `np.roll(arange(N),1)` breaks measure AND sign AND involution so the earliest
  check fires; no ODD cycle can ever isolate an involution. The working
  construction was `π ∘ σ`, needing a degenerate measure class the default
  quadrature lacks — carry the fixture PER ROW. → `L31`
- **A migration's own regression proof is the DELTA, not the green** (wide suite
  `40 failed → 12 failed`, same 12 out-of-scope rows). Pair with `-rs`. → `L30`

## 3. Config blindness — the ORPHEUS fixture facts

Generic rule: `AGENT.md` §0.6, `vv` §H2 / anti-#3 / anti-#4. Below is the
project-specific inventory of builders that SILENTLY null a channel — check each
against a concrete row before trusting a green.

- **`make_mixture` nulls TWO channels**: `sig_2` defaults to all-zero (and is
  zero on every `xs_library` A/B/C/D fixture), and there is NO `sig_l` parameter
  — it hardcodes `SigL = zeros(ng)`. Any (n,2n) or (n,α) term is identically
  nulled and its test goes vacuously green; a "balanced" fixture built through it
  is IMBALANCED by exactly `sig_l`. Build `Mixture(...)` DIRECTLY. → `L1`
- **The neighbouring SN operator fixtures carry `placeholder_materials`** (SigS /
  χ / νΣf all zero) ⟹ `F` is the ZERO operator and its reciprocity row is the
  tautology `0 == 0`. A "reuse the existing fixture" brief is a hypothesis, not
  an instruction: measure it, and record WHY each neighbour was rejected. → `L26`
- **A non-fissile mixture has NO eigenvalue** — a `solve_sn` snapshot on a
  moderator mixture is `k = 0/abs → nan`, a silent dead test. Reformulate as
  fixed-source, corroborate vs `φ = (diagΣ_t − Σ_s0ᵀ)⁻¹Q`. → `L7`
- **A SYNTHETIC fixture can null a property the REAL data exercises**, making a
  synthetic-only assertion FALSE-RED on production data (`argmax(χ)==0` read off a
  synthetic step-function χ is wrong for `pwr_like_mix()`). Pin a cumulative /
  inequality property, never a brittle exact index. → `L1`
- **Production exercises the shared mechanism on a DEGENERATE slice, so the
  general term is never activated — MANUFACTURE the activating case and make the
  load-bearing mutation the generalized term, RED on the manufactured fixture and
  GREEN on production's.** That asymmetry IS the evidence. Instances: a single
  seed level makes `pos ≡ 0`, killing the level term of an offset formula in every
  pin (`2*pos → 1*pos`) `L20`; the S/F arms feed ℓ=0 ONLY, so an S/F-only gate is
  blind to `P_ℓ(±1)` for ℓ≥1 (the same mutation on iso-only input stays 0.0)
  `L22`; a single-draw probe nulls a two-face law (use a SECOND seed and assert
  the probes differ) `L32`; a slab is the degenerate two-face case for any partner
  map, so a 2-D companion is mandatory — a cross-axis partner is SHAPE-LEGAL while
  the index sets differ `L33`.
- **MMS traps (generic bias: `vv` Mode 7).** An MMS fixed-source is INHERENTLY
  anisotropic — streaming manufactures a ℓ=1 source even for an isotropic trial,
  so a fold documenting "production sources are isotropic" is FALSIFIED by the
  MMS's own source and the EIGENVALUE path is blind to the break: verify a fold's
  MOMENT REACH ≥ the MMS source's anisotropy before trusting it as an absorb-gate
  `L18`. And every prior SN MMS ansatz VANISHED at the boundary, so the
  prescribed-inflow `q.boundary ≠ 0` path was NEVER exercised — the fix needs a
  non-vanishing-at-face ansatz with `a0 > 0` load-bearing `L3`. **That fix
  LANDED** as the §4.6 family (`build_slab_{,2g_}nonvacuum_mms_case`,
  `build_sphere_nonvacuum_mms_case`, `build_2d_cartesian_ld_stress_mms_case`) —
  `a₀=0.5`, `b₀=0.3`, anisotropic `(A_g + μ_n B_g)/W`, `[M]` 39 % angular swing
  over Γ₋. **Do NOT re-derive it; re-route it.** → `L40a`
- **⛔⛔ The vv "override the simplification bias — pick high frequency / mixed
  scales" rule is SCOPED to a SPATIAL-DISCRETIZATION claim. For a
  boundary/source-CHANNEL claim it makes the gate WORSE.** `[M]` `n_wavelengths`
  1.5 → 4.5 on the §4.6 slab multiplied the bulk truncation `L2` **×16** while
  the boundary-source error was untouched ⟹ the gate's signal-to-noise for the
  boundary claim strictly degraded. **The strengthening axis must be the one the
  claim lives on** — here the trace's ANGULAR content (`b₀/a₀`), not `k`.
  → `L40b`
- **⛔⛔ A type that spells two ACTIONS apart still has ONE method that conflates
  them, and the sibling type's guard is what hides it.** `RigidMotion` spells
  `on_points` vs `on_directions` so that applying an affine map to a direction is
  unwriteable — but `permutes()` matches `on_points`, so `[M]` it returns `None`
  for EVERY deck element carrying a translation (the wrap, a seated mirror, a
  seated rotation, a glide) and the brief's contract was unrunnable on the one arm
  the step existed to build. It had never bitten because the *sibling* type's
  guard forbids the affine part. **When a plan says "call `X.permutes(...)`",
  check WHICH action the method uses against the SEMANTIC CLASS of the arguments
  the NEW consumer will pass** — and gate on the OBSERVABLE (two motions sharing a
  linear part ⟹ bit-identical output), never on the spelling, so the gate survives
  whichever way the fix lands. → `L41a`
- **⛔ The fixture a brief NAMES for the keystone can be the degenerate one, and
  the degeneracy is invisible until you compute the intermediate the keystone
  asserts on.** `[M]` `product(4,4)`'s rotation-deck local permutation is exactly
  `arange` — the shape a wrong implementation hard-codes — while `product(4,8)`
  gives `[1,2,0,…]`. Sweep the family, pick the row where the asserted quantity is
  structurally non-trivial, and KEEP the degenerate one as a labelled control that
  says it proves nothing. → `L41b`
- **⛔ An error class that manifests only as a REFUSAL cannot be credited to the
  value row — ship a SECOND, in-range in-class mutation.** π-vs-π⁻¹ on any
  non-involution deck is out of range *by a theorem*, so the obvious mutation reds
  18 rows entirely by raising. Reversing the local assignment (right set, wrong
  assignment) reds the same rows by COMPARING, 40 of them. Both, or the value
  rows' catcher status is unproven. → `L41c`
- **⛔ When a defect was closed STRUCTURALLY, the obvious mutation reds NOTHING
  and the `catches` marker looks unearned.** `[M]` the textbook ERR-073 mutation
  (bare `argmin`) reds 0 of 78, because the `Permutation` TYPE refuses a
  non-bijection at construction one frame in; deleting the TYPE's clause reds
  exactly 1. Target the type's invariant, not the consumer — and say so in the
  marker's docstring, because that is a stronger claim than "the consumer checks".
  Fixture note: duplicate a node AND its partner so `|Γ₊| = |Γ₋|` and the extent
  guard cannot fire first. → `L41d`
- **⛔ An INTROSPECTING test adapter written to survive an unknown signature
  INFLATES the battery once the signature lands** — `inspect.signature` runs at
  test-module import, i.e. AFTER the plugin installs the mutation, so every
  `**kw` wrapper made it `pytest.fail` and unrelated rows "red": `[M]` 55/60/55
  reported, 23/11/27 true. Anti-#17 in the *flattering* direction. Retiring the
  adapter once the signature exists is both the elegance fix and the harness fix.
  → `L41e`
- **⛔ Grep where a PARAMETER is READ before designing a mutation around it.** A
  `domain_face`-style argument that reaches only an f-string is not a lever —
  `[M]` overriding it changed the binding not at all (0 reds); the defect had to
  be injected at the object the binding is actually derived from. → `L41f`
- **⛔ Bound an EXCLUSIVITY claim by running the sibling module, not by
  reasoning.** "The only row in the tree that can see X" was measured two-part:
  exclusive for the overlap/gap class tree-wide, but the sibling DOES catch a
  one-sided relaxation. Put the measured table in the docstring. → `L41g`
- **⛔ A convergence+value pair can be blind in a band where the error DECREASES.**
  `[M]` scaling the declared `q` by `(1+ε)`: at `ε ≤ 3e-4` the perturbation
  partially CANCELS the `O(h²)` truncation, so `L2(80)` drops to `0.6–0.8×` the
  honest value and BOTH `orders > 1.9` AND any `rtol` value row stay green; the
  order gate first reds at `ε = 5e-4`. So "value + rate" is not a floor — name
  the band and put the real keystone at the tier the defect lives on. → `L40d`

## 4. Reference, claim layer, and the proactive refutation

- **RULE: write the (claim-layer, pillar, truth-source) triple per gate BEFORE
  drafting it**, forcing the structural-independence cross-check from a DIFFERENT
  angle (`vv` §pillars / anti-#5,#6,#7). ORPHEUS residual: **no mesh-independent
  transport eigenvalue reference exists here** — heterogeneous refs are
  diffusion-based (~0.3 % gap) or self-referencing. Diffusion-eigenvalue is a
  cross-check with an explicit tolerance, NEVER a precision target (Issue #8).
  → `L2`
- **Two-anchor template for a pure-refactor carve:** a committed snapshot
  ("didn't move" = bit-id inheritance) is necessary-NOT-sufficient — ULP distance
  cannot tell you the pre-carve value was right. Pair with a closed-form value
  anchor (`Q/Σ_t`, `k_inf`). → `L2`
- **RULE (identity-level): the highest-value output of a proactive dispatch is
  REFUTING the plan's optimistic premises with a MEASUREMENT, before the ink
  dries.** Measured false: "bit-identical" (exact-arithmetic-only, ≤1 ULP in
  IEEE); "clean O(h²) at S16" (an interpolation floor the carve does not touch);
  "improves on flat at the boundary" (the correctly-consumed slope made it
  WORSE); "this bare-`ndarray` arm is DEAD" (a RUNTIME trace refuted it — an
  argument annotated `T` is the strongest reason to suspect the `T` arm is LIVE);
  "N pyright errors clear" (never trust a count — assert the residual verbatim);
  "the same fold applies uniformly across N solvers" (read all N bodies: the
  SN/CP/diffusion `keff` DENOMINATORS are different physics). State the
  refutation IN the plan so the implementer ships the achievable carve. → `L10`
- **Measure the proposed ACCEPTANCE CRITERION as a probe before any gate is
  written.** An AC shaped "changing X must not touch Y" is usually already true
  BY SIGNATURE, so it is unfalsifiable from the first commit — and a falsifier
  check PASSES on it, giving false confidence (`vv` Mode 8, signature-
  tautological). Gate the SIGNATURE; demote the value row to a regression floor.
  → `L24`
- **Before delivering any PRE-carve plan, `ls` the target path.** A module was
  written mid-session (untracked, mtime mid-plan); re-checking re-scored 10 of 15
  API findings as DISCHARGED and surfaced 3 gates the landed design earns. The
  reconciliation outvalues the original matrix. → `L35i`, `L33`
- **⭐ Re-measure the RED COUNT at the END of a planning dispatch, not only the
  start.** An 8-red migration debt I measured mid-plan went to `555 passed`
  before I finished (four test files migrated + a brand-new module appeared).
  A migration table whose status column is stale by minutes is worse than none —
  ship it with an explicit "READ the reconciliation section FIRST" banner, keep
  the verdict/action columns as the audit of record, and make the **residual gap
  list, measured absent by grep at the end**, the deliverable. → `L39`
- **⭐⭐ A brief's "central risk, ALREADY REALISED in the tree" is a claim to
  audit, not a premise — and it comes with an ENUMERATION that is usually
  short.** A brief said *every* SN MMS ansatz vanishes on both faces and named
  "all four case families"; the module holds **12** case classes and **four
  builders are non-vacuum by design**, with an anisotropic ansatz, SymPy
  provenance and a live L1 consumer — and the module's own header states the
  brief's claim as the thing it was written to FIX. **Enumerate the population
  yourself (`dir(module)`), never the brief's list.** When the refutation lands,
  the phase usually collapses from *build a new reference* to *re-route the
  existing one* — which is also the Pattern-2-correct answer. → `L40a`
- **⭐ The keystone's ORACLE choice decides whether it catches anything — same
  assertion shape, 8 orders of sensitivity apart.** `γ₋ψ == spec.evaluate(...)`
  is self-consistency: under a magnitude mutation BOTH sides move, green for
  every ε (it still catches delivery-COUNT). `γ₋ψ == the manufactured value
  recomputed from the reference OBJECT` reds at `ε ≈ 3e-12`. For any
  "the answer satisfies the declared condition" gate, ask **which side is the
  thing under test** — if the answer is "both", it is not a gate. → `L40e`

## 5. Tolerance is a claim — choose it per law, from measurement

- **RULE: bit-exactness is EARNED PER LAW; measure before choosing the
  assertion.** On ONE type: identity 500/500 bit-exact; associativity 500/500 on
  signed permutations but 0/500 on general rotations; `g∘g⁻¹` 0/500; a seated
  reflection's PERMUTATION exactly `arange(n)[::-1]` while its coordinates land
  5.6e-17 off. A uniform choice is a false red or a thrown-away gate. → `L35h`
- **State the law in the direction that IS a float theorem, and normalise a
  residual that scales with its input.** `on_points − on_directions == t` is NOT
  exact (`fl(a+t) − a ≠ t`); `on_points == on_directions + t` IS bit-exact
  6000/6000 because it recomputes the same expression — and is the stronger
  assertion. An ABSOLUTE `atol` for a residual scaling `O(ops × ‖t‖ × eps)` reds
  on correct code for large draws and is too loose for small: normalise by
  `max(1, ‖desired‖_∞)`. → `L35j`
- **Re-derive every tolerance from structure; retire inherited `nulp` folklore.**
  Gathers and α-folds are reduction-depth 0 ⟹ `array_equal` (a tolerance there
  would admit the bug); an `n`-term positive-summand contraction vs a `tensordot`
  is `κ=1` ⟹ `|Γ₊|·ε` — and the probe being non-negative is WHY `κ=1`; say so.
  → `L32`
- **Regression-snapshot tolerance is the CLAIM, not a magic floor**: iterative →
  `SAFETY(10) × conv_tol` read OFF the run config (the SoT shared by generator
  and test); direct → `nulp(reduction_depth)`; bit-identity enforced by `-W
  error::DriftWarning` LAYERED on top ("the gate passes" ≠ "bit-identical").
  Corollary: an ITERATED end-to-end snapshot CANNOT be the bit-identity gate for a
  zero-numerical-change refactor — committed iterated snapshots already drift
  1000s–100000s ULP from cross-run FP jitter, so descend to a single-step DIRECT
  snapshot on a fixed-seed random heterogeneous ≥2G ψ with non-zero inflow,
  captured pre-carve via `--capture-baseline`. Recipe:
  `feedback_regression_tolerance_design.md`. → `L7`
- **⛔ Never assert TIGHTER than the type's own construction invariant** — split
  into a row on the type's promise and a stronger row on the constructors' actual
  quality (`vv` anti-#16). → `L35g`
- **Gating a MORE-accurate implementation against the less-accurate one it
  replaces gates it against the error it exists to remove.** "≤1 ULP vs `np.cos`"
  FAILED at 3.75 ulp — but `np.cos(2πp/q)` is not the true value; against
  100-digit mpmath the new code is 0.57 ulp and the legacy 3.72. State the
  criterion against the arbitrary-precision value, and honestly ("within 0.57 ulp
  everywhere", never "always closer"). → `L34c`

## 6. Carve archetypes — where the load-bearing gate lives, by carve shape

**Meta-rule: the keystone is decided by whether the carve INHERITS a verified
predecessor.** Wrapping / re-expressing something verified ⟹ keystone is bit-id
INHERITANCE (necessary-NOT-sufficient — always paired with an independent value
anchor). Nothing to inherit ⟹ the keystone must be structurally independent.

- **Axis-transpose / mirror of a shipped reduction:** the implementer will
  (correctly) copy the template, and the ONE rule that genuinely FLIPS is where
  the copy goes wrong. Diff the two production bodies, enumerate every flipping
  rule, and make the load-bearing mutation the UN-flipped rule. Keep
  shared-machinery checks light. → `L11`
- **Fast path folded INTO a composed form:** principled-EQUIV (`rtol≈1e-14`),
  never 0-ULP — different reduction tree. The UNCHANGED sibling kernel MUST STAY
  `array_equal`; its red is the correct signal the aniso path got re-routed. "The
  transpose falls out free" reduces to the ONE missing leaf + a Mode-11 wrap; and
  Euclidean `Aᵀ` ≠ the metric adjoint `.H`. → `L12`
- **Operator-taxonomy family (capability-string → typed operator; first ITERATIVE
  inverse; retiring the coexisting mechanism).** ADDITIVE step: keystone =
  `array_equal` inheritance + the EXISTING closed-form anchor (do not mint a new
  reference); the runtime query form must be a derived `@property`, NEVER
  `isinstance(Protocol)` — `runtime_checkable` is class-uniform, so a
  half-adjointable composite passes; grep the LITERAL string read, not only the
  CAP constant `L13`. FIRST ITERATIVE inverse: no legacy `.solve` to inherit ⟹ the
  keystone is a dense-LU anchor + the name-earning invariant, NEVER old-vs-new
  ULP; once the composite defines `solve := inverse().apply`, `inverse().apply ≡
  solve` is a TAUTOLOGY; a raise-on-non-converge gate MUST test the TRUE residual,
  not the increment (Signature 9) `L14`. TERMINAL step: the coexistence-era
  faithfulness scaffold DELETES with the mechanism — design its PERMANENT
  structural-contract successor FIRST, migrate by mechanical RULE + a completeness
  RE-GREP (never a fixed table), retire ATOMICALLY `L15`.
- **New consumption mode of a shipped algebra (`assemble()`):** the whole point is
  the ONE-SOURCE proof — a sign-flip in the SHARED coefficient source must red
  BOTH the new gates AND the existing sweep/matvec suites; if only the new one
  reds, a twin path exists → STOP, fix, log ERR-NNN. Sparse-order ≠ apply-order ⟹
  no gate is 0-ULP. Never gate a derived SCALAR (`keff`) — Mode-12. → `L16`
- **Relocation:** behavior-free BY CONSTRUCTION ⟹ argue AGAINST a new snapshot;
  the walls are existing DriftWarning suites + Sphinx `-W` + `grep -rn
  "<old.path>"` = ZERO. But a relocation moving BOTH the SUT and its reference
  leaves self-referential `array_equal` canaries GREEN even if values shifted —
  the genuine proof is a FROZEN pre-carve baseline. → `L16/L17`
- **WRAP over an already-verified engine:** count-spy (EXACT expected count) +
  bit-id reference; structural independence from the ENGINE's OWN closed form,
  never a hand re-execution of its recurrence (ERR-032). An internal transform you
  CANNOT monkeypatch is pinned by a spy on the call ARGS plus a non-vacuity check
  (a width reversal is BLIND on a uniform mesh). → `L21`
- **UN-WELD (one closure hand-rolled at N sites → one source):** the centrepiece
  is single-source ROUTING — a Mode-11 wrap-counter asserting BOTH consumers enter
  it, counted EXACT, not `> 0`. A transpose hand-coded as a bare constant needs
  its OWN single-source gate. → `L22`
- **Correction→0 accelerator (DSA/TSA):** the property PARTITIONS the failure
  surface and FP-invariance is structurally BLIND to the machinery half — it
  catches exactly ONE of eight canonical errors. The object and rate gates are
  LOAD-BEARING, not supplements; draw the value/rate partition table FIRST.
  → `L23`
- **N-DOF separation campaign:** hunt the WELD (one record spanning two stages)
  and gate it by object `is`-identity across the strategy ladder, with the
  already-green arm as CONTROL — no value gate can see it, because both splittings
  reconstruct the same `A`. Enumerate the UNSPELLABLE states separately from the
  red ones (`FunctionSpace.__eq__` refuses same-shape different-NAME domains):
  unspellable needs a phase-ordering constraint, not a tolerance. → `L24`
- **DOMAIN NARROWING:** "will this gate go tautological?" has a THIRD answer — the
  gate BREAKS (its reference expression feeds the narrowed operator the wrong
  shape), so always simulate at the REALIZER, not the call site. The teeth CHANGE
  even when the assertion survives: the original bug becomes UNSPELLABLE and only
  the write-target family remains. The new index remap appears at ≥2 sites whose
  discriminating fixtures are COMPLEMENTARY (1-D covers one, 2-D the other) — ship
  both with activation guards. → `L29`
- **MIGRATING a narrowing's inherited surface (the phase after):** a "cannot be
  posed on the narrowed operator" gate is usually a RECIPROCITY gate in disguise
  — find the mirror object the tree already builds (the opposite-face sibling)
  before reaching for `xfail`. When a narrowing retires a private CLASSIFIER,
  exactly ONE fixture usually discriminates the replacement: name it + an
  activation guard. → `L30`
- **TYPE-COLLAPSE (N types → 1 parameterised type):** the information moves from
  the TYPE to a FIELD, so every `isinstance` / class-set gate over the collapsing
  family **stays green and stops discriminating** — inventory them FIRST and
  re-pose each onto the PARAMETER *in the same commit*. Free companion:
  `type(A().f) is type(B().f)` (`is`, not `isinstance`) — the only gate that
  catches "collapsed" into a base class + two subclasses. The guard named after
  the concept is usually necessary-NOT-sufficient, and the missing clause's
  inhabitants do the DEFERRED type's job (involution admits half-turn AND
  inversion, both of which map a face to its OPPOSITE): enumerate the property's
  CONJUGACY CLASSES against the semantic claim, then build a two-way witness pair
  with DISTINCT `match=` strings. Prefer the guard that makes the blind parameter
  UNSPELLABLE over one that needs a gate — a rigid motion's TRANSLATION is
  bit-identically invisible to every angular functional (`on_directions` drops
  it; measured identical at offset −17), so `is_linear` closes it by type and
  also implies involution for free. → `L36`
- **OPTIONAL→MANDATORY BINDING (a `None` default that silently disables a
  derivation):** nothing on the forward path changes, so a bit-id keystone is cheap
  and WORTHLESS — the keystone is a metric-sensitive reciprocity gate written
  BEFORE the metric lands, so it goes RED→GREEN (written after it can only be
  green). **MEASURED how worthless: THREE wrong bindings — dropped, SWAPPED
  (`domain`↔`codomain`), collapsed-to-one-space — each produced ZERO new reds
  across 1668 tests and 1252 constructions of the bound operator.** ⭐ **Design
  the battery around the SWAP**: it survives the extent guard (`|Γ₊|==|Γ₋|` on
  every shipped fixture), survives the refusal flag (the two spaces still
  DIFFER, so an `is_involution`-style gate stays correctly False), and changes
  no arithmetic — the ONLY catcher is an `is`-identity row naming WHICH space is
  WHICH end. Decay has THREE flavours: **DIE** (the gate can no longer CONSTRUCT
  its subject), **DECAY** (green tautology — re-pose onto "the space is the RIGHT
  one", `is` not `==`, since `FunctionSpace.__eq__` is `(name,shape)`),
  **INVERT** (the gate pins the degradation as the contract). And
  `assert_array_equal` on ANY `.H` of a newly-bound operator breaks at **2 nulp**
  (`(g·x)/g` is not an IEEE identity) — grep them all first. → `L37`
- **A binding added at a LEAF may not survive to the object the producer
  RETURNS — measure at the tier the CONSUMER sees, not at the construction
  site.** A `TensorProductOperator`/`&` wrapper derived no `domain`/`codomain`
  from its factors, so a leaf bound `Γ₊→Γ₋` reached the realizer's output as
  `None`/`None` — and the campaign's next step ("route the composition through
  `@` so the check FIRES") could not fire, because one `None` short-circuits the
  composability check. Two rules: **(a)** before crediting "the check now
  fires", compose the object PRODUCTION hands out, not the leaf you just bound;
  **(b)** an already-COMMITTED sibling step probably has the identical hole —
  check it, and ship the gap as a `strict=True` xfail naming the later step, not
  as scope creep. → `L38`
- **DELIVERY-COUNT / "is this term applied, and how many times?" carve: gate the
  governing EQUATION evaluated on the CONVERGED answer, at the tier where the
  posed system's rows ARE that equation** (a trace / interface / boundary DOF).
  For the affine BC `γ₋ψ = Lγ₊ψ + q` with `L=0`, `γ₋ψ|_f` read off the solution
  gave bit-exact `2q` (double) / `q` (single) / `0` (lost) — three distinguishable
  readings, ONE assertion, no reference solver, no tolerance, and independent of
  mesh/quadrature/materials. It beat all three candidates the brief offered.
  ⛔ **`superposition in q` is a Mode-12 NON-CATCHER for a doubled source**: the
  double IS `q→2q`, still exactly affine in `q`, so `φ(a)=φ(0)+a·s` holds for any
  `s` — a "linearity in the parameter" gate can never see a wrong CONSTANT
  multiplying that parameter. ⚠ And the OBVIOUS mutation is out-of-class
  (anti-#18): re-instating the affine operator breaks LINEARITY, so solve-level
  reds come from that, not the count — the in-class mutation is `q + q` in the
  SOURCE channel, and *a gate reddening under the affine reinstatement but not
  under `q+q` is not a single-delivery catcher*. Second in-class mutation `L := I`
  (linear!) is invisible to the `B(0)==0` gate and visible only to the trace gate.
  → `L39`
- **A "not a live bug — it is fenced" claim is a HYPOTHESIS; apply the fence's own
  predicate to the object the CONSUMER holds.** A `block_role=None` stamp fenced
  nothing (`_face_laws` collects every law with no role filter, `|B(0)|=2.5`), and
  the consequence was a HARD RAISE on Krylov (`‖Aψ−q‖/‖q‖=1.718` — an affine
  operator breaks Arnoldi) on BOTH sides of the phase blamed for it. A stamp only
  fences if somebody reads it: grep the consumer for the filter. → `L39`
- **A measured number living in a COMMENT is not a gate.** The campaign's central
  measurement (`|B(0)| = 2.5`) existed in the tree ONLY as prose in a test file.
  When auditing a landed carve, grep the measurement's NUMBER and ask whether any
  `assert` consumes it. → `L39`
- **An agreement-between-two-siblings row is NOT an is-it-correct row** — a
  type-collapse carve naturally produces `prescribed.domain is vacuum.domain`,
  which stays GREEN if BOTH are bound swapped (and the extent guard passes too,
  since `|Γ₊|==|Γ₋|`). Pair every "the two collapsed things agree" row with one
  naming the RIGHT value. → `L39`, `L37`
- **When a Mode-12 fixture's SUBJECT is retired, ask whether the successor is
  also hand-constructible** — if yes the claim migrates and usually GAINS legs
  (the zero morphism added a transpose direction + two space identities the
  apply-only affine operator never had). Two production docstrings having
  independently written the same Mode-12 argument is the signal it is
  load-bearing. → `L39`
- **A rewire's demotion test: is the retired symbol a SOURCE of the expected
  value or a FORWARDER of it?** A one-line adapter (`op.apply` whose whole body
  was `spec.evaluate(shape)`) can be inlined for free — independence lived in
  *who produces each side*, not in the adapter. The trap is then "make the new
  oracle generic" by reading the SUT's own derivation. → `L39`
- **A sentinel encoding TWO states makes the discriminating gate UNWRITABLE — say
  "this gate cannot exist; here is the TYPE that makes it exist".** `domain=None`
  means BOTH "space-generic by mathematics" (an identity is the identity on every
  space) AND "nobody bothered"; no assertion can separate one runtime value from
  itself, and that is exactly WHY the degradation is silent. → `L37`
- **Run the POSITIVE CONTROL before writing the gate and read which files it
  MISSES** — the absent reds map the blind region exactly, and are far more
  persuasive in a plan than a code citation. (`apply_inverse_metric := identity`
  reddened 30 bulk gates and **zero** boundary ones, because `domain is None`
  short-circuits the call.) Enumerate the WHOLE family's Mode-12 stabiliser, not
  the one operator: 4 of 5 shipped SN boundary laws are metric-blind, so "test the
  reflective BC harder" is a provable non-catcher. → `L37`
- **Before accepting a BIT-IDENTITY acceptance line, ask which REDUCTIONS the
  change reorders.** Zero ⟹ `array_equal` is honest. Any ⟹ the line is
  arithmetically IMPOSSIBLE (a 49-term reduction under permutation is bit-identical
  25 % of the time) — re-scope to a permutation that reorders no addition (a
  packing/relabel), never loosen. And "ZERO diff in file X" is a one-shot
  commit-scoped `git diff` check + a PERMANENT AST vocabulary gate, never a pytest
  gate. → `L37`
- **When a step says "make X mandatory", find WHERE the optional default is
  DECLARED** — on a shared base it is not scoped to your tier whatever the heading
  says (a "boundary" mandate reached the homogeneous solver, 12 test-local
  subclasses in 6 files, and a helper feeding 21 call sites: ~20 tests vs ~150).
  A committed strict-xfail set is the todo list, but **read it BY ARM** — 12 of 21
  flipped, the other 8 were a different phase. → `L37`
- **An `A ≡ B` theorem holding BY SHARED BODY is designed-GREEN under a body
  bug** — the design's own justification is why the gate cannot verify. It catches
  ARGUMENT drift only; the catcher is an independent-expression anchor written
  from raw data. The `≡` DOES have teeth where production reads the law's FACTORS
  (predicates interrogating them share nothing). → `L31`
- **Gating a REVERSE solve (transpose-solve, swap law `A.H.inverse() ≡
  A.inverse().H`):** the keystone is a FORWARD-only G-reciprocity — its arithmetic
  never calls the reverse path, so it is structurally independent BY CONSTRUCTION.
  `b` MUST be bulk-only source-carried (a random FULL `b` falsely reds even
  UNMUTATED, ~1e-1: free boundary/seed DOFs lie outside the range). A predicate
  flip MUST propagate to the capability-survival CONTRACT in the SAME landing —
  grep every `is_adjointable is False` on the flipped type. → `L19`
- **Adjoint / metric gates:** a `.H` reciprocity gate is blind exactly when
  `[G, A] = 0` — compute the COMMUTATOR at design time ("non-uniform mesh" is a
  proxy, wrong both ways; `vv` Mode 12 carries the criterion). Leaves commuting
  ALGEBRAICALLY need a second, metric-agnostic mutation. And reciprocity is a
  CONSISTENCY check, not a correctness one — forward and transpose wrong the SAME
  way reciprocate at 1e-16; pair with an object-level SUPPORT gate. → `L26/L33`
- **Carrier augmentation (a new block/DOF):** PROVE the block is CONSUMED before
  crediting any gate — zero its source and the solve MUST move. A carrier DOF's
  Hilbert metric is set by its OPERATOR ROLE, not its angular-integration weight;
  conflating them gives a ghost metric and a Mode-12 false-green. Where the DOF is
  weightless, assert-UNMOVED first — a silent re-capture discards the invariance.
  → `L18`
- **Perf gate for composition-over-fusion:** the catcher is a leaf-kernel call
  COUNT, not wall clock — and "must not scale with `n_cells`" is TOO COARSE:
  tabulate arity against EVERY axis first (one path was invariant in nx, order and
  groups but linear in ny; another already per-cell). The regression it catches is
  EXACTLY value-identical — that measurement is what promotes the count to
  catcher. A perf baseline is a (number, FIXTURE) PAIR: own the sizes, fingerprint
  them, never source them from a shared `_config` correctness may retune.
  → `L24/L25`

## 7. Snapshots, generators, and exactness

- **RULE: a snapshot generator that calls production and freezes its output is
  SELF-REFERENTIAL** — it says `production == a recording of production`, detects
  change and certifies nothing, and breaks on every signature change. INVERT it:
  compute the reference from the law's EQUATION (never by transcribing the
  implementation) and freeze THAT. Precondition: the expression must be TOTAL — a
  recording pins every bit, a derived reference only what the expression
  determines. Then the FROZEN FILE is the only thing standing between a wrong
  expression and a green gate (a recomputing harness lets generator and production
  drift TOGETHER through one shared expression) — make that structural, not
  documentary: an AST gate asserting the generator imports nothing from the
  realization layer, the harness pulls only the case registry, artefacts on disk ==
  registered cases. → `L32`
- **When a completion supersedes a retired spelling, INHERIT a frozen artefact
  generated by the SIBLING law rather than regenerating** — it predates every line
  under test, so re-baseline criterion 2 holds by construction. → `L31`
- **Making a SYMMETRY EXACT is NOT a ≤1-ULP change: exactness manufactures TIES.**
  Values moved 1.06e-14; the downstream `argsort` ordering changed in 36 of 36
  configurations and the end-to-end flux by 1.008 % — twelve orders apart, under
  one justification sentence. Checklist: (1) grep consumers for
  `argsort`/`argmin`/`sort`/`unique`/`set(`/dict-keying ON the quantity made
  exact; (2) is the sort key INJECTIVE? (non-injective ⟹ the ordering was never
  determined by the physics — a latent defect to be RULED on); (3) does the
  ambiguity converge away? (flat-in-refinement ⟹ defect); (4) SPLIT the commits —
  the ordering ruling lands first, alone. And **"the level is sorted by η" is the
  wrong functional**: sortedness is invariant under permuting equal elements, so
  two orderings differing by 1.8 % are BOTH sorted. Gate the full INDEX TUPLE
  against an independently-constructed one, plus a `kind=`-invariance row
  (quicksort/stable/heapsort/mergesort must agree bit-identically) — the
  operational proof the key is injective. → `L34`
- **A nearest-neighbour partner search is not automatically a bug — measure the
  MARGIN first.** "NN matching over a noisy set must be mis-pairing" was REFUTED
  (separation 5.0e-3 vs a 1e-16 perturbation). Keep the search rather than
  replacing it with an index formula: the sibling family has no formula, so the
  replacement would mint a twin path. And **`ref[ref] == id` passes on a
  residual-0.94 garbage map** — any self-inverse pairing satisfies an involution
  law (cf. `face_opposite → identity`); the only functional outside that
  stabiliser is the RESIDUAL. → `L34b`, `L33`

## 8. Verifying a pure-math PRIMITIVE (a group / algebra type)

- **The pillars differ from a solver's: no MMS row, no semi-analytical row.**
  Every row is closed-form; the structurally independent grounds are SymPy under
  an EXPLICIT unit parameterisation (imposing `Σnᵢ²=1` by `subs` after expansion
  does NOT fire), an external impl with a DIFFERENT ALGORITHM (quaternion vs
  Rodrigues), the Lie definition `expm(θ(vuᵀ−uvᵀ))` (dimension-generic, ~4e-14 ⟹
  gate at 1e-12), published tables, and EXACT INTEGER arithmetic — the last needs
  no reference at all and is the strongest class available. → `L35a`
- **The group-action HOMOMORPHISM `π(g∘h) = π(g)∘π(h)` is the deepest cheap
  gate** — integers only, pinning composition order, the row-vs-column convention
  and `π` vs `π⁻¹` at once (wrong order violates it on 102 of 144 pairs). VACUOUS
  on an abelian fixture. A checker returning a `bool` makes the law unaskable;
  returning the PERMUTATION makes it free. → `L35b`
- **A gate that builds its own reference through the mutated path can only see
  errors that break self-consistency** — `x @ Q` for `x @ Qᵀ` reddened the
  homomorphism but NOT the positive `permutes` gate, because a transposed action
  still maps a group-invariant set onto itself. Pair every "the map is correct"
  row with a COMPOSITION law. → `L35k`
- **An involution / order law is Mode-12 BLIND to the AFFINE part** — every
  `t ∈ span(n̂)` gives an involution while the fixed plane moves by 0.37–1.48. Gate
  the FIXED SET, never the order, for anything affine. Generalisation: for any
  `(linear, affine)` decomposition, enumerate which laws factor through the linear
  part alone — designed-green on the affine half. → `L35c`
- **A `G`-preserved weighted point set has a `G`-FIXED centroid** (3-line proof) —
  48/48 seated elements preserve it, 1/48 unseated do. That row converts "where do
  we put the origin?" from a modelling choice into a computed fact. → `L35d`
- **Bijectivity and the match WINDOW are INDEPENDENT failure modes** — a set
  off-symmetry by 1e-9 certifies under a 1e-7 window with a perfectly injective π.
  The window is a first-class correctness parameter: an explicit ARGUMENT (a module
  constant makes the "window bites" gate signature-tautological), defaulting to the
  set's minimum pairwise separation. → `L35e`
- **A `-> T | None` collapses N guards into one value; isolate from the INPUT
  side.** A negative gate asserting `is None` proves only that SOME guard fired.
  Build inputs that pass all guards but one (an unequal-weight antipodal pair
  isolates the weight guard — VACUOUS on every equal-weight fixture, i.e. every
  shipped quadrature). And **when a shipped guard rejects your test input, suspect
  the FIXTURE first** — all three reds during authoring were mine (a single
  Gram–Schmidt pass is not orthonormal enough for a 1e-12 gate; use QR).
  → `L35f`, `L35j`

## 9. Pointers

- **Characterization vs guarantee:** GUARANTEE tests carry `verifies(...)` and
  assert what IS correct; CHARACTERIZATION tests carry NO `verifies(...)` and
  bound a limitation ONE-SIDED (no upper bound, so a future fix keeps them green).
  To pin a floor a fix claims to remove, measure the floor's SCALING with the
  OTHER axis — `err(S32) < err(S16)/2` is falsifiable where "the floor is gone" is
  not. An out-of-scope defect gets a POSITIVE assert-the-defect gate with a loud
  message, NOT an xfail (xfail flips silently). → `L5`, `L16`
- **Mode-10 sub-floor terms:** producer-threading at machine precision + a
  consumed-flip ≫ tol + a no-op control; where NO isolating regime exists the
  absence of a value-improvement leg is the CORRECT signature (`vv` Mode 10).
  → `L6`
- **V&V tagging idioms** (`foundation` must not carry `verifies`; module
  `pytestmark`; `slow` on params not functions) → `feedback_vv_tagging.md` (`L9`).
  **Cross-method agreement infra** (reuse the registry schema; `max(tol_a,
  tol_b)`; tag L1 not L4; truth values MUST trace to primary citations —
  transcribing from memory invented two) → `feedback_cross_method_protocol.md`
  (`L8`). **Per-carve RECIPES** → `MEMORY.md` §3.
